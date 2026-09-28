module ome_ex
  use constants_math
  use parser_input_file, &
    only: nf, e1, e2, eta, nw, iflag_write_exk, cache_ome_read, cache_ome_write
  use parser_wannier90_tb, &
    only: material_name, norb
  use parser_optics_xatu_dim, &
    only: npointstotal, vcell, &
          norb_ex, norb_ex_cut, nv_ex, nc_ex, nband_ex, e_ex, fk_ex, &
          get_ex_index_first, print_exciton_wf, &
          rkxvector, rkyvector, rkzvector
  use exciton_envelopes, &
    only: fk_ex_der, get_fk_ex_der_k
  implicit none
  logical, save :: inter_terms_ready = .false.
!   public :: inter_terms_ready

  complex(8), allocatable :: xme_ex(:,:)
  complex(8), allocatable :: vme_ex(:,:)
  ! HANDOFF 8.45: qme_ex_inter1/2, qme_ex_inter, yme_ex_inter1/2, yme_ex_inter and vme_ex_inter1/2
  ! used to exist as eight separate (3,N,N) arrays, plus six THREAD-PRIVATE copies of them inside the
  ! parallel k-loop. Nothing ever read them individually -- they were summed at the end into
  ! xme_ex_inter = (yme1+yme2)+(qme1+qme2) and vme_ex_inter = vme1+vme2 -- so the six accumulation
  ! sites now add straight into xme_ex_inter and vme_ex_inter. At N = 1875 with 32 threads that is
  ! 32.5 GiB of scratch replaced by 10.6 GiB.
  complex(8), allocatable :: xme_ex_inter(:,:,:)
  complex(8), allocatable :: vme_ex_inter(:,:,:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Builds the excitonic optical matrix elements from the Xatu envelopes and the
  !! single-particle ones: X_n and P_n (ground state to exciton, Eqs. B2a/B2b of
  !! Taghizadeh & Pedersen, PRB 97, 205432 (2018)) and, for iflag_norder = 2, the
  !! inter-exciton X_nm and P_nm as well.
  !! This k-loop is the dominant cost of an excitonic run and does not depend on which
  !! Response is being computed, so Cache_ome_ex can skip it entirely; the read hook sits
  !! before any allocation because the peak memory is the thread-private (3,N,N)
  !! accumulators inside the loop, not the result.
  !! NOTE the V arrays are the BARE momentum P, not the Heisenberg momentum Pi; anything
  !! needing Pi must build it from X (Pi_n = -i E_n X_n). On hBN they differ by 34%.
  !! @param iflag_norder  1 = ground-state-to-exciton elements only, 2 = also inter-exciton.
  !! @return void
  subroutine get_ome_ex(iflag_norder)
    use omp_lib
    implicit none
    
    integer, intent(in) :: iflag_norder
    
    integer :: kmoment


    integer :: ibz, nn, nnp, nj, nbasis, ic, iv
    integer :: u_exk
    logical :: do_write_exk
    logical :: cache_hit

    complex(8), allocatable :: vme_ex_k(:,:)
    real(8),    allocatable :: ek(:,:)
    complex(8), allocatable :: xme_ex_band(:,:,:,:)
    complex(8), allocatable :: vme_ex_band(:,:,:,:)
    complex(8), allocatable :: berry_eigen_ex_band(:,:,:,:)

    ! PATCH: per-thread accumulators/working arrays for the parallel k-loop
    ! below. Everything that used to be allocated ONCE outside the k-loop
    ! (and shared/written directly by every call) is now allocated ONCE
    ! PER THREAD, inside the parallel region, since get_ome_inter_ex_sum_k
    ! and get_ome_gs_ex_sum_k now write into these thread-private buffers
    ! rather than module-level globals directly.
    complex(8), allocatable :: vme_ex_t(:,:), xme_ex_t(:,:)
    complex(8), allocatable :: vme_ex_k_t(:,:)
    complex(8), allocatable :: xme_ex_inter_t(:,:,:), vme_ex_inter_t(:,:,:)
    integer,    allocatable :: i_ex_table(:,:)
    complex(8), allocatable :: F_cv(:,:)
    complex(8), allocatable :: FcvH(:,:)
    complex(8), allocatable :: D_c(:,:)
    complex(8), allocatable :: A_c(:,:)
    complex(8), allocatable :: B_cc(:,:)
    complex(8), allocatable :: Y_cc(:,:)
    complex(8), allocatable :: B_vv(:,:)
    complex(8), allocatable :: Y_vv(:,:)
    complex(8), allocatable :: Uc(:,:)
    complex(8), allocatable :: Wv(:,:)
    complex(8), allocatable :: mid_cc(:,:)
    complex(8), allocatable :: mid_vv(:,:)
    complex(8), allocatable :: out_v(:,:)
    complex(8), allocatable :: out_y(:,:)
    complex(8), allocatable :: out_q(:,:)

    write(*,*) '6. Entering ome_ex'
    
    inter_terms_ready = .false.   ! ADD near the top
    
    ! PATCH: do_write_exk is now driven by the input-file flag
    ! (Write_ex_kresolved) rather than being hardcoded. Kept ANDed with
    ! iflag_norder==1 since k-resolved linear output only makes sense when
    ! computing the linear-order matrix elements in the first place --
    ! toggling the flag on for iflag_norder==2 has no effect, by design.
    do_write_exk = iflag_write_exk .and. (iflag_norder == 1)
    
    allocate(ek(npointstotal, nband_ex))
    allocate(vme_ex_band(npointstotal, 3, nband_ex, nband_ex))

    ! xme_ex_band allocated and zeroed UNCONDITIONALLY, not only for
    ! iflag_norder==2: get_ome_gs_ex_sum_k reads it on every call
    ! regardless of iflag_norder, with no guard of its own -- for
    ! SECOND-ORDER OME CACHE (HANDOFF 8.44). A hit skips the whole k-loop below, which on hBN_N75
    ! with norb_ex_cut = 1875 measured 1300 s and a 43 GB peak -- the peak being the six thread-private
    ! (3,N,N) inter-exciton accumulators, ~1 GB per thread across 32 threads. read_ome_ex_second fills
    ! xme_ex, vme_ex, xme_ex_inter, vme_ex_inter and sets inter_terms_ready, i.e. everything this
    ! routine leaves behind for iflag_norder = 2; every other array it allocates is local and freed
    ! before it returns, so an early return leaves no dangling state. Not used when the k-resolved
    ! linear file is requested, since that is produced inside the loop.
    if (iflag_norder == 2 .and. cache_ome_read .and. .not. iflag_write_exk) then
      call read_ome_ex_second(cache_hit)
      if (cache_hit) return
    end if

    ! iflag_norder==1 an unallocated read here was undefined behaviour.
    allocate(xme_ex_band(npointstotal, 3, nband_ex, nband_ex))
    xme_ex_band = (0.0d0, 0.0d0)

    if (iflag_norder == 1) then
      ek = 0.0d0
      vme_ex_band = (0.0d0, 0.0d0)
      write(*,*) '   Reading optical matrix elements (sp)...'
      call read_ome_sp_linear(iflag_norder, npointstotal, nband_ex, vme_ex_band, ek)
    end if

    if (iflag_norder == 2) then
      allocate(berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex))
      ek = 0.0d0
      vme_ex_band          = (0.0d0, 0.0d0)
      berry_eigen_ex_band  = (0.0d0, 0.0d0)
      write(*,*) '   Reading optical matrix elements (sp)...'
      ! gen_der_ex_band and shift_vector_ex_band are omitted on purpose: nothing in the excitonic
      ! path reads them, so the reader discards those records instead of filling whole-mesh arrays
      ! (HANDOFF 8.45). sigma_second_sp still asks for them, and gets them.
      call read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, &
                                  berry_eigen_ex_band = berry_eigen_ex_band, &
                                  vme_ex_band = vme_ex_band, ek = ek)
      call get_ome_sp_xme_ex_band(ek, vme_ex_band, xme_ex_band)
    end if

    allocate(vme_ex(3, norb_ex_cut))
    allocate(xme_ex(3, norb_ex_cut))
    vme_ex = (0.0d0, 0.0d0)
    xme_ex = (0.0d0, 0.0d0)

    if (do_write_exk) then
      allocate(vme_ex_k(3, norb_ex_cut))
      vme_ex_k = (0.0d0, 0.0d0)
      u_exk = 77
      call write_ome_ex_linear_kresolved_init(u_exk, material_name, norb_ex_cut)
    end if

    if (iflag_norder == 2) then
      allocate(xme_ex_inter(3, norb_ex_cut, norb_ex_cut))
      allocate(vme_ex_inter(3, norb_ex_cut, norb_ex_cut))
      xme_ex_inter = (0.0d0, 0.0d0)
      vme_ex_inter = (0.0d0, 0.0d0)

      allocate(fk_ex_der(3, norb_ex, norb_ex_cut))
      call get_fk_ex_der_k()

      nbasis = nc_ex * nv_ex
      ! PATCH: i_ex_table, F_cv, FcvH, D_c, A_c, B_cc, Y_cc, B_vv, Y_vv,
      ! Uc, Wv, mid_cc, mid_vv, out_v, out_y, out_q are NO LONGER allocated
      ! here as single shared copies -- they are allocated once per thread
      ! inside the parallel region below, since get_ome_inter_ex_sum_k is
      ! now called concurrently by multiple threads on different ibz.
    end if

    if (iflag_norder == 1) write(*,*) '   Evaluating excitonic OMEs for linear conductivity...'
    if (iflag_norder == 2) write(*,*) '   Evaluating excitonic OMEs for nonlinear conductivity...'

    ! PATCH: single parallel region wrapping the whole k-point loop,
    ! instead of get_ome_gs_ex_sum_k / get_ome_gs_ex_kresolved each opening
    ! (and tearing down) their own !$omp parallel do once per ibz. That
    ! meant paying thread-team spawn/join overhead npointstotal times for
    ! a small amount of work each time (norb_ex_cut*nc_ex*nv_ex iterations),
    ! leaving the CPU under-loaded. Now the team is spawned once; each
    ! thread processes a dynamically-scheduled subset of k-points from
    ! start to finish, accumulating into thread-private totals that are
    ! combined once at the end (same pattern as sigma_second_sp's
    ! get_sigma_shift_sp). The k-resolved file write is kept in k-order via
    ! schedule(dynamic) + ordered, same pattern as print_sigma_second_ex.
    
    
    ! Add a shared/threadprivate moment counter before the parallel region
    kmoment = -1
    
    ! DEFAULT(NONE) added 2026-09-24 (audit): an unlisted variable would silently become SHARED,
    ! which is how the Wfull/Wblk_chk race reached production earlier this session. See HANDOFF 8.35.
    !$omp parallel default(none) &
    !$omp   private(ibz, ic, iv, vme_ex_t, xme_ex_t, vme_ex_k_t, &
    !$omp           xme_ex_inter_t, vme_ex_inter_t, &
    !$omp           i_ex_table, F_cv, FcvH, D_c, A_c, &
    !$omp           B_cc, Y_cc, B_vv, Y_vv, Uc, Wv, mid_cc, mid_vv, &
    !$omp           out_v, out_y, out_q)&
    ! read-only inputs:
    !$omp   shared(npointstotal, nbasis, norb_ex_cut, nv_ex, nc_ex, nband_ex, nf, iflag_norder) &
    !$omp   shared(rkxvector, rkyvector, rkzvector) &
    !$omp   shared(ek, vme_ex_band, xme_ex_band, berry_eigen_ex_band) &
    ! accumulators: every update below is inside the !$omp critical section:
    !$omp   shared(vme_ex, xme_ex, xme_ex_inter, vme_ex_inter) &
    ! k-resolved output, written only at this thread's own ibz:
    !$omp   shared(u_exk, do_write_exk)&
    !$omp   shared(kmoment)

    allocate(vme_ex_t(3, norb_ex_cut)); vme_ex_t = (0.0d0, 0.0d0)
    allocate(xme_ex_t(3, norb_ex_cut)); xme_ex_t = (0.0d0, 0.0d0)
    if (do_write_exk) then
      allocate(vme_ex_k_t(3, norb_ex_cut))
    end if

    ! i_ex_table(ic,iv) depends only on ibz (not on iflag_norder): computed once per
    ! k-point in the loop below and passed into all three consumers, which previously
    ! each rebuilt it independently for the same ibz via get_ex_index_first (code
    ! review #8, 2026-09-23).
    allocate(i_ex_table(nc_ex, nv_ex))

    if (iflag_norder == 2) then
      allocate(xme_ex_inter_t(3, norb_ex_cut, norb_ex_cut)); xme_ex_inter_t = (0.0d0,0.0d0)
      allocate(vme_ex_inter_t(3, norb_ex_cut, norb_ex_cut)); vme_ex_inter_t = (0.0d0,0.0d0)

      allocate(F_cv (norb_ex_cut, nbasis))
      allocate(FcvH (nbasis,      norb_ex_cut))
      allocate(D_c  (nbasis,      norb_ex_cut))
      allocate(A_c  (nbasis,      norb_ex_cut))
      allocate(B_cc (nc_ex,       nc_ex))
      allocate(Y_cc (nc_ex,       nc_ex))
      allocate(B_vv (nv_ex,       nv_ex))
      allocate(Y_vv (nv_ex,       nv_ex))
      allocate(Uc   (nc_ex,       norb_ex_cut))
      allocate(Wv   (nv_ex,       norb_ex_cut))
      allocate(mid_cc(nc_ex,      norb_ex_cut))
      allocate(mid_vv(nv_ex,      norb_ex_cut))
      allocate(out_v(norb_ex_cut, norb_ex_cut))
      allocate(out_y(norb_ex_cut, norb_ex_cut))
      allocate(out_q(norb_ex_cut, norb_ex_cut))
    end if
    
    
    !$omp do schedule(dynamic) ordered
    do ibz = 1, npointstotal
      write(*,*) '   OME (ex): k-point', ibz, '/', npointstotal
!       call percentage_index(ibz, npointstotal, kmoment)   ! REPLACES the per-ibz write(*,*)

      ! Built ONCE per k-point (fix, code review #8): get_ome_gs_ex_sum_k, get_ome_gs_ex_kresolved
      ! and get_ome_inter_ex_sum_k each used to rebuild this independently via get_ex_index_first for
      ! the very same ibz (2-3x redundant work per k-point); i_ex_table(ic,iv) does not depend on
      ! which of the three consumers is asking, only on ibz, so it is shared between them below.
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          call get_ex_index_first(nf, nv_ex, nc_ex, 0, ibz, i_ex_table(ic,iv), ic, iv)
        end do
      end do

      if (iflag_norder == 1 .or. iflag_norder == 2) &
        call get_ome_gs_ex_sum_k(ibz, i_ex_table, vme_ex_band, xme_ex_band, vme_ex_t, xme_ex_t)

      if (do_write_exk) then
        call get_ome_gs_ex_kresolved(ibz, i_ex_table, vme_ex_band, vme_ex_k_t)
        !$omp ordered
        call write_ome_ex_linear_kresolved_point(u_exk, rkxvector(ibz), rkyvector(ibz), &
                                                  rkzvector(ibz), norb_ex_cut, vme_ex_k_t)
        !$omp end ordered
      end if

      if (iflag_norder == 2) &
        call get_ome_inter_ex_sum_k(                                        &
               ibz, xme_ex_band, vme_ex_band, berry_eigen_ex_band,         &
               nbasis, i_ex_table, F_cv, FcvH, D_c, A_c,                   &
               B_cc, Y_cc, B_vv, Y_vv, Uc, Wv, mid_cc, mid_vv,             &
               out_v, out_y, out_q,                                        &
               xme_ex_inter_t, vme_ex_inter_t)
    end do
    !$omp end do

    !$omp critical
      vme_ex = vme_ex + vme_ex_t
      xme_ex = xme_ex + xme_ex_t
      if (iflag_norder == 2) then
        xme_ex_inter = xme_ex_inter + xme_ex_inter_t
        vme_ex_inter = vme_ex_inter + vme_ex_inter_t
      end if
    !$omp end critical

    deallocate(vme_ex_t, xme_ex_t)
    if (do_write_exk) deallocate(vme_ex_k_t)
    deallocate(i_ex_table)
    if (iflag_norder == 2) then
      deallocate(xme_ex_inter_t, vme_ex_inter_t)
      deallocate(F_cv, FcvH, D_c, A_c)
      deallocate(B_cc, Y_cc, B_vv, Y_vv, Uc, Wv, mid_cc, mid_vv)
      deallocate(out_v, out_y, out_q)
    end if

    !$omp end parallel

    if (do_write_exk) then
      call write_ome_ex_linear_kresolved_close(u_exk)
      deallocate(vme_ex_k)
      write(*,*) '   k-resolved linear excitonic OMEs written (omeexk)'
    end if

    if (iflag_norder == 2) then
      ! xme_ex_inter and vme_ex_inter were accumulated directly in the k-loop (HANDOFF 8.45);
      ! no recombination step is needed.
      deallocate(fk_ex_der)
      inter_terms_ready = .true.   ! ADD here
    end if

    write(*,*) '   Optical matrix elements (ex) have been evaluated'
    ! Both messages used to print unconditionally, so a second-order run announced a write that never
    ! happened -- write_ome_ex_linear only fires for iflag_norder == 1 (HANDOFF 8.44).
    if (iflag_norder == 1) then
      call write_ome_ex_linear(vme_ex)
      write(*,*) '   Optical matrix elements (ex, N->GS) written'
    end if
    if (iflag_norder == 2 .and. cache_ome_write) call write_ome_ex_second()

    deallocate(ek, vme_ex_band, xme_ex_band)
    if (iflag_norder == 2) deallocate(berry_eigen_ex_band)

  end subroutine get_ome_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! x_nm = -i*p_nm/(E_n-E_m) for n != m, the standard single-particle position-from-momentum
  ! relation (feeds Y_cc/Y_vv in get_ome_inter_ex_sum_k, and the ground-state xme_ex_t sum in
  ! get_ome_gs_ex_sum_k). Guard FIXED 2026-09-23 (code review, get_ome_inter_ex_sum_k
  ! investigation): previously gated only by a magnitude clip on the RESULT (|x_nm| > 20), an
  ! undocumented, unphysical threshold of exactly the kind already proven wrong for
  ! a magnitude clip on a derivative quantity (ome_sp.f90, code review finding #6) -- calibrated empirically
  ! on MoSe2's real nv_ex=2 data (34-orbital window, 4900 k-points): 12.3% of ALL k-points had
  ! the v12-v13 pair (near-degenerate, median gap 43.7 meV, min 3.6e-8 eV) x_nm component
  ! either explode (up to 1.3e4) or get force-zeroed by the magnitude clip, an unprincipled mix
  ! of false positives (a moderate, physically large-but-valid x_nm for a perfectly resolved
  ! gap, wrongly zeroed because its magnitude happened to exceed 20) and false negatives (an
  ! unstable x_nm from a barely-resolved gap, kept because its magnitude happened to land under
  ! 20). Replaced with the SAME gap-based eps_deg=1.0d-4 Ha (2.7 meV) criterion used throughout
  ! sigma_second_sp.f90 (shift_shiftvector, get_shg_intens_sp) for "is this pair
  ! near-degenerate, hence gauge-ambiguous": x_nm is large but VALID for any resolved
  ! (non-degenerate) gap, however small, and only genuinely unreliable when the two
  ! single-particle bands are degenerate to within eps_deg, where p_nm itself becomes
  ! gauge-dependent within the degenerate subspace (same root cause class as ome_sp.f90's
  ! near-degeneracy issues elsewhere). No magnitude clip is applied to the result.
  subroutine get_ome_sp_xme_ex_band(ek, vme_ex_band, xme_ex_band)
    implicit none
    real(8),    intent(in)  :: ek(npointstotal, nband_ex)
    complex(8), intent(in)  :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    integer    :: ibz, i, j, nj
    real(8)    :: de
    complex(8), parameter :: ci = (0.0d0, 1.0d0)
    real(8), parameter :: eps_deg = 1.0d-4   ! same degeneracy window as sigma_second_sp.f90

    xme_ex_band = (0.0d0, 0.0d0)
    !$omp parallel do collapse(2) default(none) schedule(static) private(ibz, nj, i, j, de) &
    !$omp   shared(npointstotal, nband_ex, ek, vme_ex_band, xme_ex_band)
    do ibz = 1, npointstotal
      do nj = 1, 3
        do i = 1, nband_ex
          do j = 1, nband_ex
            de = ek(ibz,i) - ek(ibz,j)
            if (abs(de) >= eps_deg) then
              xme_ex_band(ibz,nj,i,j) = -ci / de * vme_ex_band(ibz,nj,i,j)
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine get_ome_sp_xme_ex_band

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PATCH (both fixes applied here):
  !  1. i_ex_table(ic,iv) does not depend on nn -- it was previously
  !     recomputed via get_ex_index_first inside the "do nn" loop, i.e.
  !     norb_ex_cut times per (ic,iv) for an answer that never changes.
  !     FURTHER FIXED 2026-09-23 (code review #8): i_ex_table also does not
  !     depend on which of the three ome_ex consumers (this routine,
  !     get_ome_gs_ex_kresolved, get_ome_inter_ex_sum_k) is asking, only on
  !     ibz -- each used to rebuild it independently for the same k-point.
  !     Now built ONCE per ibz by the caller (get_ome_ex) and passed in as
  !     an intent(in) argument, shared by all three.
  !  2. The internal !$omp parallel do is REMOVED: this subroutine is now
  !     called from inside the single parallel region opened in
  !     get_ome_ex, once per ibz, by whichever thread owns that ibz -- an
  !     inner parallel region here would either be ignored (nesting
  !     disabled, the common default) or spawn a costly nested team.
  !     Results are accumulated into vme_ex_t/xme_ex_t (intent(inout)),
  !     which are thread-private buffers owned by the calling thread,
  !     rather than into the module-level vme_ex/xme_ex directly -- since
  !     multiple threads now run this concurrently for different ibz, and
  !     module-level vme_ex/xme_ex are shared across all of them.
  subroutine get_ome_gs_ex_sum_k(ibz, i_ex_table, vme_ex_band, xme_ex_band, vme_ex_t, xme_ex_t)
    implicit none
    integer,    intent(in)    :: ibz
    integer,    intent(in)    :: i_ex_table(nc_ex, nv_ex)
    complex(8), intent(in)    :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in)    :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(inout) :: vme_ex_t(3, norb_ex_cut)
    complex(8), intent(inout) :: xme_ex_t(3, norb_ex_cut)
    integer :: nn, ic, iv, nj

    do nn = 1, norb_ex_cut
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          do nj = 1, 3
            vme_ex_t(nj,nn) = vme_ex_t(nj,nn) + fk_ex(i_ex_table(ic,iv),nn)*vme_ex_band(ibz,nj,iv,nv_ex+ic)
            xme_ex_t(nj,nn) = xme_ex_t(nj,nn) + fk_ex(i_ex_table(ic,iv),nn)*xme_ex_band(ibz,nj,iv,nv_ex+ic)
          end do
        end do
      end do
    end do
  end subroutine get_ome_gs_ex_sum_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PATCH: same fixes as get_ome_gs_ex_sum_k -- i_ex_table hoisted out
  ! of the nn loop, then (code review #8, 2026-09-23) hoisted out of this
  ! subroutine entirely: it is built once per ibz by the caller and passed
  ! in as an intent(in) argument, shared with get_ome_gs_ex_sum_k and
  ! get_ome_inter_ex_sum_k rather than rebuilt separately here. The
  ! internal !$omp parallel do was also removed since this now runs inside
  ! the single outer parallel region (one thread, one ibz, at a time).
  ! vme_ex_k remains a per-call intent(out) result (already caller-private
  ! in get_ome_ex via vme_ex_k_t), so no accumulator restructuring was
  ! needed here beyond removing the nested OMP directive.
  subroutine get_ome_gs_ex_kresolved(ibz, i_ex_table, vme_ex_band, vme_ex_k)
    implicit none
    integer,    intent(in)  :: ibz
    integer,    intent(in)  :: i_ex_table(nc_ex, nv_ex)
    complex(8), intent(in)  :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: vme_ex_k(3, norb_ex_cut)
    integer :: nn, ic, iv, nj

    vme_ex_k = (0.0d0, 0.0d0)

    do nn = 1, norb_ex_cut
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          do nj = 1, 3
            vme_ex_k(nj,nn) = vme_ex_k(nj,nn) + fk_ex(i_ex_table(ic,iv),nn)*vme_ex_band(ibz,nj,iv,nv_ex+ic)
          end do
        end do
      end do
    end do
  end subroutine get_ome_gs_ex_kresolved

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get_ome_inter_ex_sum_k
!
! IMPORTANT: the (ic,icp) contraction against B_cc/Y_cc must be carried out
! separately for EACH value of iv, and the (iv,ivp) contraction against
! B_vv/Y_vv must be carried out separately for EACH value of ic, and only
! THEN summed. iv (resp. ic) is a shared "spectator" index tying together
! both matrix-element factors in the original sum
!
!   V1_{nn,nn'} = sum_iv  sum_{ic,icp} conjg(f_{(ic,iv),nn}) B_cc(ic,icp) f_{(icp,iv),nn'}
!   V2_{nn,nn'} = -sum_ic sum_{iv,ivp} conjg(f_{(ic,iv),nn}) B_vv(iv,ivp) f_{(ic,ivp),nn'}
!
! It is NOT valid to pre-sum fk_ex over iv (resp. ic) and then form a single
! bilinear contraction from the pre-summed vectors: that silently introduces
! spurious "cross" terms with iv /= iv' (resp. ic /= ic') that do not belong
! in the physical sum. Each iv (resp. ic) therefore requires its own small
! ZGEMM contraction, accumulated into the output with beta = 1.
!
! PATCH: qme_ex_inter1/2, yme_ex_inter1/2, vme_ex_inter1/2 are now
! intent(inout) THREAD-PRIVATE accumulators passed in by the caller
! (get_ome_ex, inside the parallel region) instead of module-level
! globals written directly. This subroutine is now called concurrently
! by multiple threads for different ibz; writing straight into the
! module-level qme_ex_inter1 etc. (as the previous version did) would be
! a data race the moment more than one thread is active. The working
! arrays (F_cv, D_c, B_cc, ...) were already intent(inout) dummies here,
! and are now allocated per-thread by the caller rather than once
! globally, for the same reason.
! FURTHER PATCH (code review #8, 2026-09-23): i_ex_table changed from
! intent(inout), built internally here via get_ex_index_first at the top of
! every call, to intent(in): it depends only on ibz, not on which consumer
! is asking, so get_ome_ex now builds it once per k-point and shares it
! between this routine, get_ome_gs_ex_sum_k and get_ome_gs_ex_kresolved.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_inter_ex_sum_k(                                      &
               ibz, xme_ex_band, vme_ex_band, berry_eigen_ex_band,       &
               nbasis, i_ex_table, F_cv, FcvH, D_c, A_c,                 &
               B_cc, Y_cc, B_vv, Y_vv, Uc, Wv, mid_cc, mid_vv,           &
               out_v, out_y, out_q,                                      &
               xme_ex_inter_t, vme_ex_inter_t)
    implicit none

    integer,    intent(in)    :: ibz
    integer,    intent(in)    :: nbasis
    complex(8), intent(in)    :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in)    :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in)    :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
    integer,    intent(in)    :: i_ex_table(nc_ex, nv_ex)
    complex(8), intent(inout) :: F_cv(norb_ex_cut, nbasis)
    complex(8), intent(inout) :: FcvH(nbasis, norb_ex_cut)
    complex(8), intent(inout) :: D_c(nbasis, norb_ex_cut)
    complex(8), intent(inout) :: A_c(nbasis, norb_ex_cut)
    complex(8), intent(inout) :: B_cc(nc_ex, nc_ex)
    complex(8), intent(inout) :: Y_cc(nc_ex, nc_ex)
    complex(8), intent(inout) :: B_vv(nv_ex, nv_ex)
    complex(8), intent(inout) :: Y_vv(nv_ex, nv_ex)
    complex(8), intent(inout) :: Uc(nc_ex, norb_ex_cut)
    complex(8), intent(inout) :: Wv(nv_ex, norb_ex_cut)
    complex(8), intent(inout) :: mid_cc(nc_ex, norb_ex_cut)
    complex(8), intent(inout) :: mid_vv(nv_ex, norb_ex_cut)
    complex(8), intent(inout) :: out_v(norb_ex_cut, norb_ex_cut)
    complex(8), intent(inout) :: out_y(norb_ex_cut, norb_ex_cut)
    complex(8), intent(inout) :: out_q(norb_ex_cut, norb_ex_cut)
    ! Two accumulators, not six (HANDOFF 8.45). The Y and Q pieces both belong to X_nm and the two
    ! V pieces both belong to P_nm; they were only ever summed, so they are summed here.
    complex(8), intent(inout) :: xme_ex_inter_t(3, norb_ex_cut, norb_ex_cut)
    complex(8), intent(inout) :: vme_ex_inter_t(3, norb_ex_cut, norb_ex_cut)

    integer    :: ic, iv, icp, ivp, nj, i_ex, idx, nn
    real(8)    :: berry_shift
    complex(8), parameter :: ci    = (0.0d0, 1.0d0)
    complex(8), parameter :: cone  = (1.0d0, 0.0d0)
    complex(8), parameter :: czero = (0.0d0, 0.0d0)

    !--- Step 1: F_cv(nn,idx) / FcvH(idx,nn) over the full (ic,iv) basis ---
    ! (i_ex_table is now built once by the caller, get_ome_ex, and passed in; see PATCH above)
    ! (used only for the qme_ex_inter1/2 terms, which contract over the
    !  FULL basis index and do not have the spectator-index restriction)
    do ic = 1, nc_ex
      do iv = 1, nv_ex
        idx  = (ic-1)*nv_ex + iv
        i_ex = i_ex_table(ic,iv)
        do nn = 1, norb_ex_cut
          F_cv(nn, idx) = fk_ex(i_ex, nn)
        end do
      end do
    end do
    do nn = 1, norb_ex_cut
      do idx = 1, nbasis
        FcvH(idx, nn) = conjg(F_cv(nn, idx))
      end do
    end do

    do nj = 1, 3

      !================================================================
      ! vme_ex_inter1 / yme_ex_inter1:
      !   sum_iv  Uc(iv)^H * B_cc * Uc(iv)   (and with Y_cc for yme)
      ! where Uc(iv)(ic,nn) = fk_ex(i_ex_table(ic,iv), nn)
      !================================================================
      ! B_cc holds the velocity matrix elements (feeds vme_ex_inter1);
      ! Y_cc holds the POSITION matrix elements xme_ex_band (feeds
      ! yme_ex_inter1) — these are physically distinct quantities, Y_cc
      ! is NOT simply B_cc with the diagonal removed. xme_ex_band already
      ! has an exactly-zero diagonal by construction (see
      ! get_ome_sp_xme_ex_band), so no further zeroing is strictly needed;
      ! it is kept below only for clarity/parity with the reference code.
      do ic = 1, nc_ex
        do icp = 1, nc_ex
          B_cc(ic,icp) = vme_ex_band(ibz, nj, nv_ex+ic, nv_ex+icp)
          Y_cc(ic,icp) = xme_ex_band(ibz, nj, nv_ex+ic, nv_ex+icp)
        end do
      end do
      do ic = 1, nc_ex
        Y_cc(ic,ic) = czero
      end do

      out_v = czero
      out_y = czero
      do iv = 1, nv_ex
        do ic = 1, nc_ex
          Uc(ic,:) = fk_ex(i_ex_table(ic,iv), :)
        end do

        ! mid_cc = B_cc * Uc ; out_v += Uc^H * mid_cc
        call zgemm('N', 'N', nc_ex, norb_ex_cut, nc_ex, &
                    cone,  B_cc, nc_ex, &
                           Uc,   nc_ex, &
                    czero, mid_cc, nc_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nc_ex, &
                    cone, Uc,     nc_ex, &
                          mid_cc, nc_ex, &
                    cone, out_v,  norb_ex_cut)

        ! mid_cc = Y_cc * Uc ; out_y += Uc^H * mid_cc
        call zgemm('N', 'N', nc_ex, norb_ex_cut, nc_ex, &
                    cone,  Y_cc, nc_ex, &
                           Uc,   nc_ex, &
                    czero, mid_cc, nc_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nc_ex, &
                    cone, Uc,     nc_ex, &
                          mid_cc, nc_ex, &
                    cone, out_y,  norb_ex_cut)
      end do
      vme_ex_inter_t(nj,:,:) = vme_ex_inter_t(nj,:,:) + out_v
      xme_ex_inter_t(nj,:,:) = xme_ex_inter_t(nj,:,:) + out_y

      !================================================================
      ! vme_ex_inter2 / yme_ex_inter2:
      !   -sum_ic  Wv(ic)^H * B_vv * Wv(ic)  (and with Y_vv for yme)
      ! where Wv(ic)(iv,nn) = fk_ex(i_ex_table(ic,iv), nn)
      !================================================================
      ! Same distinction as above: B_vv is the velocity matrix (feeds
      ! vme_ex_inter2), Y_vv is the POSITION matrix xme_ex_band (feeds
      ! yme_ex_inter2) — not a zero-diagonal copy of B_vv.
      !
      ! NOTE the band-index order here is TRANSPOSED relative to B_cc/Y_cc above
      ! (B_vv(iv,ivp) = vme(ivp,iv), not vme(iv,ivp)): this is NOT a copy-paste error.
      ! Matching Taghizadeh & Pedersen, PRB 97, 205432 (2018) Eq. (B3a) term-by-term to
      ! this matrix form forces exactly this transpose for the valence ("hole") block --
      ! the electron/hole asymmetry inherent to the BSE envelope formalism. Verified
      ! analytically 2026-09-23, and empirically: X_nm Hermiticity holds to 9.6e-14 on
      ! MoSe2 (nv_ex=2, the first dataset where this transpose is not a no-op; hBN's own
      ! nv_ex=1 test data can never exercise it). Do not "fix" this without re-deriving
      ! Eq. (B3a) first.
      do ivp = 1, nv_ex
        do iv = 1, nv_ex
          B_vv(iv,ivp) = vme_ex_band(ibz, nj, ivp, iv)
          Y_vv(iv,ivp) = xme_ex_band(ibz, nj, ivp, iv)
        end do
      end do
      do iv = 1, nv_ex
        Y_vv(iv,iv) = czero
      end do

      out_v = czero
      out_y = czero
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          Wv(iv,:) = fk_ex(i_ex_table(ic,iv), :)
        end do

        ! mid_vv = B_vv * Wv ; out_v += Wv^H * mid_vv
        call zgemm('N', 'N', nv_ex, norb_ex_cut, nv_ex, &
                    cone,  B_vv, nv_ex, &
                           Wv,   nv_ex, &
                    czero, mid_vv, nv_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nv_ex, &
                    cone, Wv,     nv_ex, &
                          mid_vv, nv_ex, &
                    cone, out_v,  norb_ex_cut)

        ! mid_vv = Y_vv * Wv ; out_y += Wv^H * mid_vv
        call zgemm('N', 'N', nv_ex, norb_ex_cut, nv_ex, &
                    cone,  Y_vv, nv_ex, &
                           Wv,   nv_ex, &
                    czero, mid_vv, nv_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nv_ex, &
                    cone, Wv,     nv_ex, &
                          mid_vv, nv_ex, &
                    cone, out_y,  norb_ex_cut)
      end do
      vme_ex_inter_t(nj,:,:) = vme_ex_inter_t(nj,:,:) - out_v
      xme_ex_inter_t(nj,:,:) = xme_ex_inter_t(nj,:,:) - out_y

      !================================================================
      ! qme_ex_inter1 += i * FcvH^T * D_c
      !================================================================
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx = (ic-1)*nv_ex + iv
          D_c(idx, :) = fk_ex_der(nj, i_ex_table(ic,iv), :)
        end do
      end do
      call zgemm('T', 'N', norb_ex_cut, norb_ex_cut, nbasis, &
                  ci,    FcvH, nbasis,       &
                         D_c,  nbasis,        &
                  czero, out_q, norb_ex_cut)
      xme_ex_inter_t(nj,:,:) = xme_ex_inter_t(nj,:,:) + out_q

      !================================================================
      ! qme_ex_inter2 += i * FcvH^T * A_c,  A_c = -i * fk_ex * berry_shift
      ! (the "-i" is essential — dropping it rotates this term by 90 deg)
      !================================================================
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx  = (ic-1)*nv_ex + iv
          i_ex = i_ex_table(ic,iv)
          berry_shift = dble(berry_eigen_ex_band(ibz,nj,nv_ex+ic,nv_ex+ic)) &
                      - dble(berry_eigen_ex_band(ibz,nj,iv,iv))
          A_c(idx, :) = -ci * fk_ex(i_ex, :) * berry_shift
        end do
      end do
      call zgemm('T', 'N', norb_ex_cut, norb_ex_cut, nbasis, &
                  ci,    FcvH, nbasis,       &
                         A_c,  nbasis,        &
                  czero, out_q, norb_ex_cut)
      xme_ex_inter_t(nj,:,:) = xme_ex_inter_t(nj,:,:) + out_q

    end do ! nj

  end subroutine get_ome_inter_ex_sum_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SECOND-ORDER EXCITONIC OME CACHE  (draft -- to be pasted into module ome_ex)
!
! Motivation: get_ome_ex's k-loop is the dominant cost of an excitonic second-order run and is repeated
! identically for every Response branch. Measured on hBN_N75 with norb_ex_cut = 1875: 1300 s and a 43 GB
! peak per branch, against ~5 GB for the conductivity kernel that follows. Four branches (shg,
! electrooptic, rectification, shift_shiftvector) therefore spend ~87 min recomputing one 22 min object.
!
! What is cached: everything get_ome_ex produces for iflag_norder = 2 -- xme_ex(3,N), vme_ex(3,N),
! xme_ex_inter(3,N,N), vme_ex_inter(3,N,N). Only the X arrays are needed by the second-order kernels
! (method B uses them directly; method A rebuilds Pi from them, Pi_n = -i E_n X_n,
! Pi_nm = i(E_n-E_m) X_nm), but the V arrays are cached too so that a hit is a FAITHFUL substitute for
! the routine rather than a partial one -- leaving them unallocated or zeroed would turn any future use
! into silent garbage. NOTE they are the BARE momentum P, not Pi (CLAUDE.md); caching them does not make
! them safe to use in a second-order formula.
!
! Format: unformatted stream. Text would be ~5x larger and slow to parse. Total payload is
! 2*3*N^2 + 2*3*N complex(8): 338 MB at N = 1875, 1.24 GB at N = 3600, 3.0 GB at N = 5625.
!
! Staleness: the header carries the identifying parameters AND the full exciton energy list, which is a
! strong fingerprint of the Xatu solution. It does NOT fingerprint the Wannier90 file, so a changed
! tight-binding model with an unchanged exciton spectrum would go undetected. The cache is therefore
! OPT-IN (see the call-site patch below) rather than automatic.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ome_ex_cache_name(fname)
    implicit none
    character(len=*), intent(out) :: fname
    fname = 'ome_second_ex_'//trim(material_name)//'.omeex2'
  end subroutine ome_ex_cache_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write xme_ex and xme_ex_inter. Called after the k-loop has filled them.
  subroutine write_ome_ex_second()
    implicit none
    character(len=14), parameter :: MAGIC = 'OPTICX-OMEEX2 '
    integer,           parameter :: VERSION = 1
    character(len=256) :: fname
    integer :: u, ios, m

    if (.not. allocated(xme_ex) .or. .not. allocated(xme_ex_inter) .or. &
        .not. allocated(vme_ex) .or. .not. allocated(vme_ex_inter)) then
      write(*,*) '   WARNING (write_ome_ex_second): excitonic OMEs not allocated, nothing cached.'
      return
    end if

    call ome_ex_cache_name(fname)
    open(newunit=u, file=trim(fname), form='unformatted', access='stream', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) '   WARNING (write_ome_ex_second): cannot open '//trim(fname)//', not cached.'
      return
    end if

    write(u) MAGIC, VERSION
    write(u) len_trim(material_name)
    write(u) trim(material_name)
    write(u) npointstotal, nv_ex, nc_ex, norb_ex_cut
    write(u) e_ex(1:norb_ex_cut)                      ! fingerprint of the Xatu solution
    write(u) xme_ex(1:3, 1:norb_ex_cut)
    write(u) vme_ex(1:3, 1:norb_ex_cut)
    do m = 1, norb_ex_cut                             ! one contiguous slab per m, so a smaller
      write(u) xme_ex_inter(1:3, 1:norb_ex_cut, m)    ! norb_ex_cut can be sliced out on read
    end do
    do m = 1, norb_ex_cut
      write(u) vme_ex_inter(1:3, 1:norb_ex_cut, m)
    end do
    close(u)
    write(*,'(A,A,A,F8.2,A)') '   Second-order excitonic OMEs cached to ', trim(fname), ' (', &
         (6.0d0*dble(norb_ex_cut)**2 + 6.0d0*dble(norb_ex_cut))*16.0d0/1048576.0d0, ' MB)'
  end subroutine write_ome_ex_second

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Try to fill xme_ex and xme_ex_inter from the cache. ok = .false. means "compute them normally";
! the caller must then run the k-loop as before. Never stops the run: a bad cache is a miss, not a
! failure, EXCEPT for a dimension/energy mismatch, which is reported loudly because silently reusing
! elements from a different grid or exciton set would be a wrong answer rather than a slow one.
  subroutine read_ome_ex_second(ok)
    implicit none
    logical, intent(out) :: ok
    character(len=14), parameter :: MAGIC = 'OPTICX-OMEEX2 '
    integer,           parameter :: VERSION = 1
    real(8),           parameter :: ETOL = 1.0d-10
    character(len=256) :: fname
    character(len=14)  :: magic_in
    character(len=256) :: mat_in
    integer :: u, ios, ver, nlen, np_in, nv_in, nc_in, nex_in, m
    real(8),    allocatable :: e_in(:)
    complex(8), allocatable :: buf(:,:)
    logical :: there

    ok = .false.
    call ome_ex_cache_name(fname)
    inquire(file=trim(fname), exist=there)
    if (.not. there) return

    open(newunit=u, file=trim(fname), form='unformatted', access='stream', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) return

    read(u, iostat=ios) magic_in, ver
    if (ios /= 0 .or. magic_in /= MAGIC .or. ver /= VERSION) then
      write(*,*) '   Cache '//trim(fname)//' has a different format, ignoring it.'
      close(u); return
    end if
    read(u, iostat=ios) nlen
    if (ios /= 0 .or. nlen < 1 .or. nlen > len(mat_in)) then; close(u); return; end if
    mat_in = ''
    read(u, iostat=ios) mat_in(1:nlen)
    read(u, iostat=ios) np_in, nv_in, nc_in, nex_in
    if (ios /= 0) then; close(u); return; end if

    if (trim(mat_in) /= trim(material_name) .or. np_in /= npointstotal .or. &
        nv_in /= nv_ex .or. nc_in /= nc_ex) then
      write(*,*) '   ERROR (read_ome_ex_second): '//trim(fname)//' was written for a DIFFERENT system.'
      write(*,'(A,A,A,I0,A,I0,A,I0)') '          cached: ', trim(mat_in), '  nk = ', np_in, &
           '  nv = ', nv_in, '  nc = ', nc_in
      write(*,'(A,A,A,I0,A,I0,A,I0)') '          wanted: ', trim(material_name), '  nk = ', npointstotal, &
           '  nv = ', nv_ex, '  nc = ', nc_ex
      write(*,*) '          Delete it or move it aside; refusing to mix exciton sets.'
      close(u); stop 1
    end if
    if (nex_in < norb_ex_cut) then
      write(*,'(A,I0,A,I0,A)') '   Cache holds only ', nex_in, ' excitons but ', norb_ex_cut, &
           ' are requested; recomputing.'
      close(u); return
    end if

    allocate(e_in(nex_in))
    read(u, iostat=ios) e_in
    if (ios /= 0) then; deallocate(e_in); close(u); return; end if
    if (maxval(abs(e_in(1:norb_ex_cut) - e_ex(1:norb_ex_cut))) > ETOL) then
      write(*,*) '   ERROR (read_ome_ex_second): cached exciton energies differ from the current ones'
      write(*,'(A,ES12.4)') '          max|dE| = ', &
           maxval(abs(e_in(1:norb_ex_cut) - e_ex(1:norb_ex_cut)))
      write(*,*) '          The .eigval/.states files have changed. Delete '//trim(fname)//'.'
      deallocate(e_in); close(u); stop 1
    end if
    deallocate(e_in)

    if (allocated(xme_ex))       deallocate(xme_ex)
    if (allocated(vme_ex))       deallocate(vme_ex)
    if (allocated(xme_ex_inter)) deallocate(xme_ex_inter)
    if (allocated(vme_ex_inter)) deallocate(vme_ex_inter)
    allocate(xme_ex(3, norb_ex_cut), vme_ex(3, norb_ex_cut))
    allocate(xme_ex_inter(3, norb_ex_cut, norb_ex_cut), vme_ex_inter(3, norb_ex_cut, norb_ex_cut))

    ! Both arrays were written for the CACHED count nex_in, so read full slabs and take the corner.
    allocate(buf(3, nex_in))

    read(u, iostat=ios) buf                                  ! X_n
    if (ios /= 0) go to 900
    xme_ex(:,:) = buf(:, 1:norb_ex_cut)
    read(u, iostat=ios) buf                                  ! P_n
    if (ios /= 0) go to 900
    vme_ex(:,:) = buf(:, 1:norb_ex_cut)

    ! Both slab loops MUST run over the CACHED count nex_in, not over norb_ex_cut: the slabs are laid
    ! out back to back, so stopping early would leave the file positioned inside the unread remainder
    ! of X_nm and the next array would be read from the wrong offset. (Caught by the round-trip test.)
    do m = 1, nex_in                                         ! X_nm, one slab of 3*nex_in per m
      read(u, iostat=ios) buf
      if (ios /= 0) go to 900
      if (m <= norb_ex_cut) xme_ex_inter(:, :, m) = buf(:, 1:norb_ex_cut)
    end do
    do m = 1, nex_in                                         ! P_nm
      read(u, iostat=ios) buf
      if (ios /= 0) go to 900
      if (m <= norb_ex_cut) vme_ex_inter(:, :, m) = buf(:, 1:norb_ex_cut)
    end do
    deallocate(buf)

    close(u)
    inter_terms_ready = .true.
    ok = .true.
    write(*,'(A,I0,A,A)') '   Second-order excitonic OMEs read from cache (', norb_ex_cut, &
         ' excitons) -- exciton k-loop skipped: ', trim(fname)
    return

900 continue                                                 ! truncated or corrupt: treat as a miss
    write(*,*) '   Cache '//trim(fname)//' is truncated or unreadable, recomputing.'
    deallocate(buf)
    deallocate(xme_ex, vme_ex, xme_ex_inter, vme_ex_inter)
    close(u)
    ok = .false.
    return
    write(*,'(A,I0,A,A)') '   Second-order excitonic OMEs read from cache (', norb_ex_cut, &
         ' excitons) -- k-loop skipped: ', trim(fname)
  end subroutine read_ome_ex_second

  subroutine write_ome_ex_linear(vme_ex)
    implicit none
    complex(8), intent(in) :: vme_ex(3, norb_ex_cut)
    integer :: nn, nj
    open(10, file='ome_linear_ex_'//trim(material_name)//'.omeex')
    write(10,*) 1
    do nn = 1, norb_ex_cut
      write(10,*) nn, (dble(vme_ex(nj,nn)), dimag(vme_ex(nj,nn)), nj=1,3)
    end do
    close(10)
  end subroutine write_ome_ex_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_ome_ex_linear_kresolved_init(unitno, mat_name, norb_ex_cut)
    implicit none
    integer,          intent(in) :: unitno, norb_ex_cut
    character(len=*), intent(in) :: mat_name
    open(unitno, file='ome_linear_ex_k_'//trim(mat_name)//'.omeexk', status='replace')
    write(unitno,*) 1
    write(unitno,*) norb_ex_cut
  end subroutine write_ome_ex_linear_kresolved_init

  subroutine write_ome_ex_linear_kresolved_point(unitno, kx, ky, kz, norb_ex_cut, vme_ex_k)
    implicit none
    integer,    intent(in) :: unitno, norb_ex_cut
    real(8),    intent(in) :: kx, ky, kz
    complex(8), intent(in) :: vme_ex_k(3, norb_ex_cut)
    integer :: nn
    write(unitno,*) kx, ky, kz
    do nn = 1, norb_ex_cut
      write(unitno,*) nn, dble(vme_ex_k(1,nn)), dimag(vme_ex_k(1,nn)), &
                          dble(vme_ex_k(2,nn)), dimag(vme_ex_k(2,nn)), &
                          dble(vme_ex_k(3,nn)), dimag(vme_ex_k(3,nn))
    end do
  end subroutine write_ome_ex_linear_kresolved_point

  subroutine write_ome_ex_linear_kresolved_close(unitno)
    implicit none
    integer, intent(in) :: unitno
    close(unitno)
  end subroutine write_ome_ex_linear_kresolved_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ome_sp_linear(iflag_norder, npointstotal, nband_ex, vme_ex_band, ek)
    implicit none
    integer,    intent(in)  :: iflag_norder, npointstotal, nband_ex
    real(8),    intent(out) :: ek(npointstotal, nband_ex)
    complex(8), intent(out) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    integer :: ibz, i, j, iflag_r
    real(8) :: a1, a2, a3, b1, b2, b3, b4, b5, b6
    open(10, file='ome_linear_sp_'//trim(material_name)//'.omesp')
    read(10,*) iflag_r
    do ibz = 1, npointstotal
      read(10,*) a1, a2, a3, (ek(ibz,j), j=1,nband_ex)
      do i = 1, nband_ex
        do j = 1, nband_ex
          read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
          vme_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
          vme_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
          vme_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
        end do
      end do
    end do
    close(10)
  end subroutine read_ome_sp_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, &
                                    berry_eigen_ex_band, gen_der_ex_band, &
                                    shift_vector_ex_band, vme_ex_band, ek, &
                                    vme_abs_der_ex_band, vme_abs_der_found)
   implicit none
   integer, intent(in)  :: iflag_norder, npointstotal, nband_ex
   ! optional: derivative of |v|, appended after the per-k records; absent in older files
   real(8),    intent(out), optional :: vme_abs_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
   logical,    intent(out), optional :: vme_abs_der_found
   integer :: ios_vd
   real(8), intent(out) :: ek(npointstotal, nband_ex)
   complex(8), intent(out) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
   complex(8), intent(out) :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
   ! OPTIONAL (HANDOFF 8.45): the shift vector and the sum-rule generalized derivative. Records for
   ! them sit between the per-k records the excitonic path does need, so they must still be READ --
   ! but a caller that never looks at them can omit the arrays and let the reader discard each record
   ! into a one-k-point scratch. get_ome_ex does exactly that: it used to allocate, zero and fill
   ! both over the whole mesh and never read a value. At nband_ex = 88 with 3600 k-points that is
   ! 6 GB of allocation and the matching disk traffic for nothing.
   real(8),    intent(out), optional :: shift_vector_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
   complex(8), intent(out), optional :: gen_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
   real(8)    :: sv_skip(3, 3, nband_ex, nband_ex)
   complex(8) :: gd_skip(3, 3, nband_ex, nband_ex)
   integer :: ibz, iflag_r, npts_r, nband_r
   real(8), allocatable :: rkx_r(:), rky_r(:), rkz_r(:)

   ios_vd = -1
   open(10, file='ome_nonlinear_sp_'//trim(material_name)//'.omesp', &
        form='unformatted', access='stream', status='old')

   read(10) iflag_r
   read(10) npts_r, nband_r
   ! The .omesp file is unformatted stream I/O with no self-describing record boundaries: if it
   ! was generated for a different Bandlist or k-grid than the one now being requested (e.g. a
   ! stale file reused with OME_sp = none after changing the input), the per-k reads below would
   ! silently consume the wrong number of bytes per record instead of erroring. Validate the
   ! file's own recorded size against what the caller actually needs before reading further.
   if (npts_r /= npointstotal .or. nband_r /= nband_ex) then
      write(*,*) 'ERROR (read_ome_sp_nonlinear): .omesp file was generated for a different'
      write(*,*) '       grid/Bandlist (file has',npts_r,'k-points,',nband_r,'bands; expected', &
                  npointstotal,nband_ex,'). Regenerate it with OME_sp = nonlinear.'
      stop 1
   end if
   allocate(rkx_r(npts_r), rky_r(npts_r), rkz_r(npts_r))
   read(10) rkx_r, rky_r, rkz_r

   do ibz = 1, npointstotal
      read(10) ek(ibz,:)
      read(10) vme_ex_band(ibz,:,:,:)
      read(10) berry_eigen_ex_band(ibz,:,:,:)
      if (present(shift_vector_ex_band)) then
         read(10) shift_vector_ex_band(ibz,:,:,:,:)
      else
         read(10) sv_skip
      end if
      if (present(gen_der_ex_band)) then
         read(10) gen_der_ex_band(ibz,:,:,:,:)
      else
         read(10) gd_skip
      end if
   end do

   if (present(vme_abs_der_ex_band)) then
      read(10, iostat=ios_vd) vme_abs_der_ex_band
      if (ios_vd /= 0) vme_abs_der_ex_band = 0.0d0
      if (present(vme_abs_der_found)) vme_abs_der_found = (ios_vd == 0)
   end if


   close(10)
  end subroutine read_ome_sp_nonlinear
  
  !!!!!
  
  ! DO NOT RENABLE WITHOUT UPDATING THE WRITE ROUTINE AS WELL
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, &
!                                     berry_eigen_ex_band, gen_der_ex_band, &
!                                     shift_vector_ex_band, vme_ex_band, ek)
!     implicit none
!     integer, intent(in)  :: iflag_norder, npointstotal, nband_ex
!     real(8), intent(out) :: ek(npointstotal, nband_ex)
!     real(8), intent(out) :: shift_vector_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
!     complex(8), intent(out) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
!     complex(8), intent(out) :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
!     complex(8), intent(out) :: gen_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
!     integer :: ibz, i, j, nj, iflag_r
!     real(8) :: a1, a2, a3, b1, b2, b3, b4, b5, b6
!     open(10, file='ome_nonlinear_sp_'//trim(material_name)//'.omesp')
!     read(10,*) iflag_r
!     do ibz = 1, npointstotal
!       read(10,*) a1, a2, a3, (ek(ibz,j), j=1,nband_ex)
!       do i = 1, nband_ex
!         do j = 1, nband_ex
!           read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
!           vme_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
!           vme_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
!           vme_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
!           read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
!           berry_eigen_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
!           berry_eigen_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
!           berry_eigen_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
!           do nj = 1, 3
!             read(10,*) a1, a2, a3, b1, b2, b3
!             shift_vector_ex_band(ibz,nj,1,i,j) = b1
!             shift_vector_ex_band(ibz,nj,2,i,j) = b2
!             shift_vector_ex_band(ibz,nj,3,i,j) = b3
!           end do
!           do nj = 1, 3
!             read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
!             gen_der_ex_band(ibz,nj,1,i,j) = cmplx(b1,b2,8)
!             gen_der_ex_band(ibz,nj,2,i,j) = cmplx(b3,b4,8)
!             gen_der_ex_band(ibz,nj,3,i,j) = cmplx(b5,b6,8)
!           end do
!         end do
!       end do
!     end do
!     close(10)
!   end subroutine read_ome_sp_nonlinear

end module ome_ex