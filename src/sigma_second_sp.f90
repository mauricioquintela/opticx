module sigma_second_sp
  use constants_math
  use parser_input_file, &
  only:nf,e1,e2,eta,nw,response_text,broadening_type_text
  use parser_wannier90_tb, &
  only:material_name
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex,e_ex,fk_ex, &
  get_ex_index_first,print_exciton_wf, & !routines
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use ome_ex, &
  only:read_ome_sp_nonlinear !routine
  implicit none

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_second_sp(nwp,nwq)
    implicit none 
    !in/out
    integer :: nwp,nwq
    
    !here
    integer :: iflag_norder
    integer :: ibz,j
    
    ! derivative of |v| (numerical, gauge invariant), read from the .omesp file; heap allocated
    real*8, allocatable :: vme_abs_der_ex_band(:,:,:,:,:)
    logical :: vme_abs_der_found
    ! gauge-fixed (parallel-transported) complex generalized derivative of v, (v^c_nm);k^a; needed by SHG only

    !energies and vme in k-mesh and auxiliary arrays (sp)
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    
    !energies and VME (ex)
    dimension wp(nw)
    dimension sigma_w_sp(3,3,nw)

    real*8 ek
    real*8 wp
    real*8 shift_vector_ex_band
    complex*16 vme_ex_band
    complex*16 berry_eigen_ex_band
    complex*16 gen_der_ex_band   
    complex*16 sigma_w_sp 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '10. Entering sigma_second_sp'

    !read matrix elements from file
    allocate(vme_abs_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
    call read_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,berry_eigen_ex_band, &
                                   gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek, &
                                   vme_abs_der_ex_band,vme_abs_der_found)
    write(*,*) '    Optical matrix elements (sp) have been read from file'
    if (response_text == 'shift_shiftvector' .and. .not. vme_abs_der_found) then
      write(*,*) 'ERROR (sigma_second_sp): the .omesp file has no derivative of |v| (it was written by an'
      write(*,*) '       older version). Regenerate it with OME_sp = nonlinear; shift_shiftvector needs it.'
      stop 1
    end if

    !compute shift conductivity
    if (nwp.eq.1 .and. nwq.eq.(-1)) then
      call get_sigma_shift_sp(npointstotal,nband_ex,berry_eigen_ex_band, &
                        gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek,vme_abs_der_ex_band)
      !write(*,*) 'The optical response',response_text,'has been evaluated'
    end if




    deallocate(vme_abs_der_ex_band)

  end subroutine get_sigma_second_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine get_sigma_shift_sp(npointstotal, nband_ex, berry_eigen_ex_band, &
                               gen_der_ex_band, shift_vector_ex_band, vme_ex_band, ek, vme_abs_der_ex_band)
  implicit none

  integer,    intent(in) :: nband_ex, npointstotal
  real*8,     intent(in) :: ek(npointstotal, nband_ex)
  real*8,     intent(in) :: shift_vector_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
  complex*16, intent(in) :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
  complex*16, intent(in) :: gen_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
  complex*16, intent(in) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
  real*8,     intent(in) :: vme_abs_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)

  real*8,     allocatable :: wp(:)
  real*8,     allocatable :: shift_vector_w(:,:,:)
  complex*16, allocatable :: sigma_w_sp(:,:,:,:)
  real*8 :: eta2

  integer :: ibz, i, j, nj, njp

  ! Per-thread private arrays (allocated inside PARALLEL block)
  real*8,     allocatable :: e_nband(:)
  complex*16, allocatable :: vme_nband(:,:,:)
  real*8,     allocatable :: shift_vector_nband(:,:,:,:)
  complex*16, allocatable :: gen_der_nband(:,:,:,:)
  real*8,     allocatable :: vme_abs_der_nband(:,:,:,:)
  real*8,     allocatable :: shift_vector_w_t(:,:,:)
  complex*16, allocatable :: sigma_w_sp_t(:,:,:,:)

  ! ---------------------------------------------------------------------------
  ! Allocate FIRST, then initialize
  ! ---------------------------------------------------------------------------
  allocate(wp(nw))
  allocate(sigma_w_sp(3, 3, 3, nw))
  allocate(shift_vector_w(3, 3, nw))

  call initialize_sigma_second_arrays(nw, wp, eta2, sigma_w_sp)

  shift_vector_w = 0.0d0   ! initialize doesn't touch this one

  write(*,*) '    Evaluating shift conductivity (sp)...'
  if (response_text == 'shift_sumrule') then
    write(*,*) '    WARNING: shift_sumrule is unreliable unless the band window is large (remote bands):'
    write(*,*) '             it is ~0 for two-band models and disagrees with shift_shiftvector on GeS.'
  end if

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP   SHARED(npointstotal, nband_ex, ek, nw, vme_ex_band, &
  !$OMP          shift_vector_ex_band, gen_der_ex_band, vme_abs_der_ex_band, &
  !$OMP          wp, eta2, sigma_w_sp, shift_vector_w) &
  !$OMP   PRIVATE(ibz, i, j, nj, njp, &
  !$OMP           e_nband, vme_nband, shift_vector_nband, &
  !$OMP           gen_der_nband, vme_abs_der_nband, &
  !$OMP           shift_vector_w_t, sigma_w_sp_t)

  allocate(e_nband(nband_ex))
  allocate(vme_nband(3, nband_ex, nband_ex))
  allocate(shift_vector_nband(3, 3, nband_ex, nband_ex))
  allocate(gen_der_nband(3, 3, nband_ex, nband_ex))
  allocate(vme_abs_der_nband(3, 3, nband_ex, nband_ex))
  allocate(sigma_w_sp_t(3, 3, 3, nw))
  allocate(shift_vector_w_t(3, 3, nw))

  vme_abs_der_nband    = 0.0d0
  sigma_w_sp_t     = (0.0d0, 0.0d0)
  shift_vector_w_t = 0.0d0

  !$OMP DO SCHEDULE(DYNAMIC)
  do ibz = 1, npointstotal

    do i = 1, nband_ex
      e_nband(i) = ek(ibz, i)
      do j = 1, nband_ex
        do nj = 1, 3
          vme_nband(nj, i, j) = vme_ex_band(ibz, nj, i, j)
          do njp = 1, 3
            shift_vector_nband(nj, njp, i, j) = shift_vector_ex_band(ibz, nj, njp, i, j)
            gen_der_nband(nj, njp, i, j)      = gen_der_ex_band(ibz, nj, njp, i, j)
            vme_abs_der_nband(nj, njp, i, j)      = vme_abs_der_ex_band(ibz, nj, njp, i, j)
          end do
        end do
      end do
    end do

    call get_shift_intens_sp(nband_ex, nw, e_nband, vme_nband, &
         shift_vector_nband, gen_der_nband, vme_abs_der_nband, &
         shift_vector_w_t, wp, eta2, sigma_w_sp_t)

  end do
  !$OMP END DO

  !$OMP CRITICAL
    sigma_w_sp     = sigma_w_sp     + sigma_w_sp_t
    shift_vector_w = shift_vector_w + shift_vector_w_t
  !$OMP END CRITICAL

  deallocate(e_nband, vme_nband, shift_vector_nband, &
             gen_der_nband, vme_abs_der_nband, &
             sigma_w_sp_t, shift_vector_w_t)

  !$OMP END PARALLEL

  call print_sigma_second_sp(nw, wp, sigma_w_sp, shift_vector_w)
  write(*,*) '    Shift conductivity (sp) has been printed'

  deallocate(wp, sigma_w_sp, shift_vector_w)

end subroutine get_sigma_shift_sp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift-current kernel at one k-point: the IPA expression Eq. (9) of Esteve-Paredes
  !! et al., npj Comput. Mater. 11, 13 (2025), accumulated with the -pi/2 prefactor and a
  !! normalised delta.
  !! Im[I^abc + I^acb] has two pieces. The shift-vector piece -rho_b rho_c cos(th_b-th_c)
  !! (R^{a,b} + R^{a,c}) is the familiar one; the amplitude-gradient piece
  !! sin(th_b-th_c)(rho_b d_a rho_c - rho_c d_a rho_b) vanishes for b = c but not otherwise,
  !! and omitting it left hBN's yxy at 0.324x the exact value with C3 violated. It is built
  !! from d_a|v^c| rather than the complex derivative, which gauge jumps corrupt.
  !! @param e_nband, vme_nband            Energies and velocity elements at this k-point.
  !! @param shift_vector_nband            R^{a,b}_nm.
  !! @param vme_abs_der_nband             d_a|v^c_nm|, for the amplitude-gradient term.
  !! @param wp, eta2                      Frequency grid and broadening.
  !! @param sigma_w_sp, shift_vector_w    Accumulated onto.
  !! @return void
  subroutine get_shift_intens_sp(nband_ex, nw, e_nband, vme_nband, &
     shift_vector_nband, gen_der_nband, vme_abs_der_nband, &
     shift_vector_w, wp, eta2, sigma_w_sp)
  implicit none
 
  ! ---------------------------------------------------------------------------
  ! Dummy arguments
  ! ---------------------------------------------------------------------------
  integer,    intent(in)    :: nband_ex, nw
  real*8,     intent(in)    :: e_nband(nband_ex)
  complex*16, intent(in)    :: vme_nband(3, nband_ex, nband_ex)
  real*8,     intent(in)    :: shift_vector_nband(3, 3, nband_ex, nband_ex)
  complex*16, intent(in)    :: gen_der_nband(3, 3, nband_ex, nband_ex)
  real*8,     intent(in)    :: vme_abs_der_nband(3, 3, nband_ex, nband_ex)
  real*8,     intent(in)    :: wp(nw)
  real*8,     intent(in)    :: eta2
  real*8,     intent(inout) :: shift_vector_w(3, 3, nw)     ! per-thread accumulator
  complex*16, intent(inout) :: sigma_w_sp(3, 3, 3, nw)      ! per-thread accumulator
 
  ! ---------------------------------------------------------------------------
  ! Local variables  (all stack-allocated, automatically private per call)
  ! ---------------------------------------------------------------------------
  integer    :: iw, nj, njp, njpp, nn, nnp
  real*8     :: factor1, fnn, fnnp, delta_nnp
  complex*16 :: abc(3, nband_ex, nband_ex)
  complex*16 :: shift1, shift2, shift
  complex*16 :: rb, rc
  real*8     :: dw, gb, gc
  real*8, parameter :: amp_tiny = 1.0d-10
  ! same value and meaning as clip_threshold in ome_sp.f90 (used there to zero garbage shift vectors / Berry
  ! connections): g below has the dimension of a length like R^{a,b}, so it gets the same treatment
  real*8, parameter :: g_clip = 50.0d0
  ! energy window (Hartree, 2.7 meV) below which two bands count as degenerate; the code paper's SI Note 7 finds the
  ! shift current stable for such thresholds between 1e-2 and 1e-4 eV
  real*8, parameter :: eps_deg = 1.0d-4
  logical :: degen(nband_ex)
  complex*16, parameter :: ci = (0.0d0, 1.0d0)
  ! Frequency-independent kernel for the interband term, computed ONCE per (nn,nnp) pair below and
  ! reused across all nw frequencies (fix, code review #7, 2026-09-23): shift1/shift2 (shift_sumrule)
  ! and shift/rb/rc/gb/gc (shift_shiftvector) depend only on nn, nnp and the cartesian indices, never
  ! on wp(iw) -- only delta_nnp does. The old code recomputed the full 3x3x3 kernel from scratch at
  ! every one of the nw frequencies (400 in the validation runs, up to 30000 in production -- see
  ! CLAUDE.md "Performance context"), i.e. up to 30000x more work than necessary for a quantity that
  ! does not depend on frequency at all.
  complex*16 :: shift_kernel(3, 3, 3)
 
  ! ---------------------------------------------------------------------------
  ! Module-level read-only parameters used below:
  !   nv_ex               -- valence band count
  !   npointstotal        -- total k-points (for normalisation)
  !   vcell               -- unit cell volume
  !   broadening_type_text -- 'gaussian' or 'lorentzian'
  !   response_text        -- 'shift_sumrule', 'shift_shiftvector', etc.
  ! These are read-only globals; reading them from multiple threads is safe.
  ! ---------------------------------------------------------------------------
 
  ! bands that are (near-)degenerate with another band of the window: their matrix elements are arbitrary mixtures there
  do nn = 1, nband_ex
    degen(nn) = .false.
    do nnp = 1, nband_ex
      if (nnp /= nn .and. abs(e_nband(nn) - e_nband(nnp)) < eps_deg) degen(nn) = .true.
    end do
  end do

  do nn = 1, nband_ex

    ! Fermi occupation: valence = 1, conduction = 0
    if (nn .le. nv_ex) then
      fnn = 1.0d0
    else
      fnn = 0.0d0
    end if

    do nnp = 1, nband_ex

      if (nnp .le. nv_ex) then
        fnnp = 1.0d0
      else
        fnnp = 0.0d0
      end if

      factor1 = fnn - fnnp

      ! -----------------------------------------------------------------------
      ! Frequency-independent kernel for this (nn,nnp) pair, computed once and
      ! reused for every iw below (see the shift_kernel declaration comment).
      ! -----------------------------------------------------------------------
      shift_kernel = (0.0d0, 0.0d0)
      if (fnn .ne. fnnp) then
        do nj = 1, 3
          do njp = 1, 3
            do njpp = 1, 3

              ! SIGN CONVENTION (fixed 2026-09-22): the accumulated prefactor is -pi/2, i.e.
              ! sigma = (i*pi/2V) sum f_nn' (I^{abc}+I^{acb}) delta(w - w_nn'), paper Eq. 9 of Esteve-Paredes
              ! et al., npj Comput. Mater. 11, 13 (e = 1 in a.u.). It used to be +pi/2, the opposite sign; verified
              ! against the excitonic code's IPA limit (synthetic non-interacting excitons, hBN), which follows
              ! paper Eq. 10 (see HANDOFF.md section 8.7).
              if (response_text == 'shift_sumrule') then
                ! WARNING (unreliable, kept as is): the sum-rule generalised derivative needs a large band window
                ! (remote bands). For a 2-band model its imaginary part vanishes identically (result ~ 0), and on
                ! GeS with the full 27-band window it disagrees strongly with 'shift_shiftvector' (HANDOFF.md 8.7).
                ! Sum-rule form using generalised derivative
                shift1 = -complex(0.0d0, 1.0d0) / (e_nband(nn) - e_nband(nnp)) * &
                           vme_nband(njp,  nn, nnp) * gen_der_nband(njpp, nj, nnp, nn)
                shift2 = -complex(0.0d0, 1.0d0) / (e_nband(nn) - e_nband(nnp)) * &
                           vme_nband(njpp, nn, nnp) * gen_der_nband(njp,  nj, nnp, nn)
                shift  = -complex(0.0d0, 1.0d0) * (shift1 + shift2)

                shift_kernel(nj, njp, njpp) = shift

              end if

              if (response_text == 'shift_shiftvector') then
                ! Shift-vector form (Nagaosa, 10.1103/PhysRevX.10.041041) + amplitude-gradient term (completed 2026-09-22).
                ! The shift-vector expression alone,  -(R^{a,b}_{n'n} - R^{a,c}_{nn'}) v^c_{nn'} v^b_{n'n} / w_{nn'}^2,
                ! is only part of  Im[I^{abc}+I^{acb}], I^{abc}_{nn'} = r^b_{nn'} r^{c;a}_{n'n}. With
                ! r^b_{nn'} = rho_b e^{i th_b}:
                !   Im[I^{abc}+I^{acb}] = -rho_b rho_c cos(th_b-th_c) (R^{a,b}+R^{a,c})              (shift-vector part)
                !                         + sin(th_b-th_c) (rho_b d_a rho_c - rho_c d_a rho_b)        (amplitude gradient)
                ! The second term vanishes for b = c but not for b /= c (hBN: yxy came out 0.324x the correct value).
                ! It is added here as  shift_B = -i [ r^b conj(r^c) g_c + r^c conj(r^b) g_b ],
                ! g_c = d_a ln rho_c = d_a|v^c_{nn'}| / |v^c_{nn'}| - (v^a_nn - v^a_n'n')/w_{nn'}, and only needs the
                ! derivative of the MODULUS of v, which is gauge invariant (the derivative of the complex v is
                ! corrupted by gauge jumps between separately diagonalised k-points; the shift vector above
                ! keeps the clipping it always had). The term is dropped where |v^c| ~ 0 (phase undefined) and,
                ! like the shift vector, where |g| > 50 bohr (near-degenerate bands).
                shift = -(shift_vector_nband(nj, njp,  nnp, nn) - &
                           shift_vector_nband(nj, njpp, nn,  nnp)) * &
                         vme_nband(njpp, nn, nnp) * vme_nband(njp, nnp, nn) / &
                         (e_nband(nn) - e_nband(nnp))**2
                dw = e_nband(nn) - e_nband(nnp)
                rb = -ci * vme_nband(njp,  nn, nnp) / dw
                rc = -ci * vme_nband(njpp, nn, nnp) / dw
                gb = 0.0d0
                gc = 0.0d0
                if (abs(vme_nband(njp,  nn, nnp)) > amp_tiny) gb = vme_abs_der_nband(nj, njp,  nn, nnp) / &
                      abs(vme_nband(njp,  nn, nnp)) - dble(vme_nband(nj, nn, nn) - vme_nband(nj, nnp, nnp)) / dw
                if (abs(vme_nband(njpp, nn, nnp)) > amp_tiny) gc = vme_abs_der_nband(nj, njpp, nn, nnp) / &
                      abs(vme_nband(njpp, nn, nnp)) - dble(vme_nband(nj, nn, nn) - vme_nband(nj, nnp, nnp)) / dw
                ! Essential/accidental degeneracies (e.g. the zone boundary of a non-symmorphic lattice) make the
                ! finite-difference derivatives garbage (|g| ~ 1e4-1e5). The whole pair (both b and c halves of the
                ! term) is then dropped, like the shift vector is zeroed by clip_threshold in ome_sp.f90; keeping
                ! only one half would leave an uncancelled remainder.
                if (abs(gb) > g_clip .or. abs(gc) > g_clip .or. degen(nn) .or. degen(nnp)) then
                  gb = 0.0d0
                  gc = 0.0d0
                end if
                shift = shift - ci * ( rb*conjg(rc)*gc + rc*conjg(rb)*gb )

                shift_kernel(nj, njp, njpp) = shift

              end if

              if (response_text == 'shift_gender') then
                ! Placeholder: numerical generalised derivative (Toni's paper)
                ! TODO: implement
              end if

            end do  ! njpp
          end do  ! njp
        end do  ! nj
      end if  ! fnn .ne. fnnp

      ! -----------------------------------------------------------------------
      ! Frequency loop: only delta_nnp depends on iw; the kernel above does not.
      ! -----------------------------------------------------------------------
      do iw = 1, nw

        ! Broadening lineshape
        if (trim(broadening_type_text) == 'gaussian') then
          delta_nnp = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * &
            exp(-0.5d0/(eta2**2) * (wp(iw) - e_nband(nn) + e_nband(nnp))**2)
        else if (trim(broadening_type_text) == 'lorentzian') then
          delta_nnp = 1.0d0/pi * aimag(1.0d0 / (wp(iw) - e_nband(nn) + &
            e_nband(nnp) - complex(0.0d0, eta2)))
        else
          ! Default to lorentzian
          delta_nnp = 1.0d0/pi * aimag(1.0d0 / (wp(iw) - e_nband(nn) + &
            e_nband(nnp) - complex(0.0d0, eta2)))
        end if

        do nj = 1, 3
          do njp = 1, 3

            ! Shift vector spectral function (2-index, always computed)
            shift_vector_w(nj, njp, iw) = shift_vector_w(nj, njp, iw) + &
              1.0d0/(dble(npointstotal)*vcell) * factor1 * &
              shift_vector_nband(nj, njp, nn, nnp) * delta_nnp

            if (fnn .ne. fnnp) then
              do njpp = 1, 3
                sigma_w_sp(nj, njp, njpp, iw) = sigma_w_sp(nj, njp, njpp, iw) - &
                  0.5d0*pi / (dble(npointstotal)*vcell) * factor1 * shift_kernel(nj, njp, njpp) * delta_nnp
              end do
            end if

          end do  ! njp
        end do  ! nj

      end do  ! iw

    end do  ! nnp
  end do  ! nn
 
end subroutine get_shift_intens_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  







  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_sigma_second_arrays(nw, wp, eta2, sigma_w)
    use omp_lib
    implicit none

    integer,    intent(in)  :: nw
    real(8),    intent(out) :: wp(nw)
    real(8),    intent(out)  :: eta2          ! written back as out below
    complex(8), intent(out) :: sigma_w(3,3,3,nw)

    real(8) :: wrange
    integer :: i

    ! scalar outputs that were previously implicit via host association
    real(8) :: eta2_local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Serial on purpose (audit 2026-09-24): both loops are O(nw) of trivial memory traffic, run once
    ! per calculation. Spawning a thread team costs more than the work itself. HANDOFF 8.35.
    do i = 1, nw
      wp(i)      = 0.0d0
      sigma_w(:,:,:,i) = (0.0d0, 0.0d0)
    end do

    wrange = e2 - e1

    do i = 1, nw
      wp(i) = (e1 + wrange / dble(nw) * dble(i-1)) / 27.211385d0
    end do

    eta2 = eta / 27.211385d0   ! scalar — no parallelism needed

  end subroutine initialize_sigma_second_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_second_sp(nw, wp, sigma_w_sp, shift_vector_w)
    use omp_lib
    implicit none

    integer,     intent(in) :: nw
    real(8),     intent(in) :: wp(nw)
    complex(8),  intent(in) :: sigma_w_sp(3,3,3,nw)
    real(8),     intent(in) :: shift_vector_w(3,3,nw)

    integer :: iw
    real(8) :: feps

    ! Unit-conversion factor is loop-invariant — compute once
    feps = 6.623618d-03 * 1.0d+06 * (27.211386d0**(-2)) &
         * 5.291772d-11 * 1.0d+09

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(90,  file='shift_sp_lengthgauge_'//trim(material_name)//'.dat')
    open(100, file='shift_vector.dat')

    ! --- OPTION A: simple parallel loop with ordered I/O ---
    ! Ordered writes preserve frequency ordering in the output files.
    ! serial on purpose (audit 2026-09-24): this loop is pure file I/O and every iteration was
    ! inside !$omp ordered, which serialises it completely -- the parallel wrapper only added
    ! thread spawn and synchronisation cost. HANDOFF 8.35.
    do iw = 1, nw

      write(90,*) wp(iw)*27.211385d0, &
        realpart(feps*sigma_w_sp(1,1,1,iw)), realpart(feps*sigma_w_sp(1,1,2,iw)), &
        realpart(feps*sigma_w_sp(1,1,3,iw)), realpart(feps*sigma_w_sp(1,2,1,iw)), &
        realpart(feps*sigma_w_sp(1,2,2,iw)), realpart(feps*sigma_w_sp(1,2,3,iw)), &
        realpart(feps*sigma_w_sp(1,3,1,iw)), realpart(feps*sigma_w_sp(1,3,2,iw)), &
        realpart(feps*sigma_w_sp(1,3,3,iw)), realpart(feps*sigma_w_sp(2,1,1,iw)), &
        realpart(feps*sigma_w_sp(2,1,2,iw)), realpart(feps*sigma_w_sp(2,1,3,iw)), &
        realpart(feps*sigma_w_sp(2,2,1,iw)), realpart(feps*sigma_w_sp(2,2,2,iw)), &
        realpart(feps*sigma_w_sp(2,2,3,iw)), realpart(feps*sigma_w_sp(2,3,1,iw)), &
        realpart(feps*sigma_w_sp(2,3,2,iw)), realpart(feps*sigma_w_sp(2,3,3,iw)), &
        realpart(feps*sigma_w_sp(3,1,1,iw)), realpart(feps*sigma_w_sp(3,1,2,iw)), &
        realpart(feps*sigma_w_sp(3,1,3,iw)), realpart(feps*sigma_w_sp(3,2,1,iw)), &
        realpart(feps*sigma_w_sp(3,2,2,iw)), realpart(feps*sigma_w_sp(3,2,3,iw)), &
        realpart(feps*sigma_w_sp(3,3,1,iw)), realpart(feps*sigma_w_sp(3,3,2,iw)), &
        realpart(feps*sigma_w_sp(3,3,3,iw))

      write(100,*) wp(iw)*27.211385d0, &
        shift_vector_w(1,1,iw), shift_vector_w(1,2,iw), shift_vector_w(1,3,iw), &
        shift_vector_w(2,1,iw), shift_vector_w(2,2,iw), shift_vector_w(2,3,iw), &
        shift_vector_w(3,1,iw), shift_vector_w(3,2,iw), shift_vector_w(3,3,iw)

    end do

    close(90)
    close(100)

  end subroutine print_sigma_second_sp

!!!!
  ! Writes the single-particle SHG conductivity sigma^{abc}(2*omega; omega, omega) to
  ! shg_sp_lengthgauge_<material>.dat. Columns: hbar*omega (eV) -- the FUNDAMENTAL (driving) photon energy,
  ! then for a=x,y,z; b=x,y,z; c=x,y,z
  ! (c fastest): Re, Im. Units: uA*nm/V^2, same a.u.->SI factor as the shift/excitonic-SHG conductivities.
  ! AXIS CONVENTION CHANGED 2026-09-24 (HANDOFF 8.33): column 1 used to be 2*hbar*omega. It is now
  ! hbar*omega, matching every other spectrum opticx writes (sigma_first_sp/ex, the shift conductivity),
  ! so that a two-photon resonance of an excitation at energy E appears at hbar*omega = E/2 and a
  ! one-photon resonance at hbar*omega = E. Files produced before this date carry the old 2*hbar*omega
  ! axis and will be mis-plotted by a factor of 2 if read with the new convention.
  ! No spin degeneracy factor g is included, consistent with the rest of opticx; the overall sign of e is
  ! not examined (see get_sigma_shg_sp).
  subroutine print_shg_second_sp(nw, wp, sigma_shg)
    implicit none
    integer,    intent(in) :: nw
    real(8),    intent(in) :: wp(nw)
    complex(8), intent(in) :: sigma_shg(3,3,3,nw)
    integer :: iw, ia, ib, ic
    real(8) :: feps

    feps = (6.623618d-03)*(1.0d+06)*(27.211386d0**(-2))*(5.291772d-11)*(1.0d+09)

    open(92, file='shg_sp_lengthgauge_'//trim(material_name)//'.dat')
    write(92,'(A)') '# hbar*omega(eV) [FUNDAMENTAL, not 2*hbar*omega] | sigma^{abc}(2w;w,w) (Re, Im) in uA nm/V^2, abc = xxx,xxy,xxz,xyx,...,zzz'

    do iw = 1, nw
      write(92,'(ES18.10,54ES18.10)') wp(iw)*27.211385d0, &
        ( ( ( real(feps*sigma_shg(ia,ib,ic,iw)), aimag(feps*sigma_shg(ia,ib,ic,iw)), &
              ic=1,3 ), ib=1,3 ), ia=1,3 )
    end do

    close(92)

  end subroutine print_shg_second_sp
!!!!

end module sigma_second_sp

