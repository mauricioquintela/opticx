module ome_sp
   use constants_math
   use parser_wannier90_tb, &
      only:material_name,nR,nRvec,norb,R,shop,hhop,rhop_c
   use parser_optics_xatu_dim, &
      only:npointstotal,rkxvector,rkyvector,rkzvector, &
      nband_ex,nband_index,nv_ex,nc_ex

   implicit none

   ! PATCH: dimensionality flags computed ONCE (via set_active_flags) instead
   ! of being recomputed inside every call to get_vme_kernels_ome /
   ! get_berry_eigen_fourpoint. Read-only after initialization, safe to
   ! share across OMP threads without being PRIVATE.
   logical, save :: active_x = .false., active_y = .false., active_z = .false.
   logical, save :: active_flags_set = .false.

   ! PATCH: named clipping threshold instead of a magic literal duplicated
   ! in two places (berry_eigen and shift_vector divergence clipping).
   real(8), parameter :: clip_threshold = 50.0d0

   ! ------------------------------------------------------------------------
   ! Degenerate-multiplet basis prescription of Sipe & Shkrebtii, PRB 61, 5337
   ! (2000): Eq. (A4) and the text below their Eq. (58). At a k point where two
   ! or more bands are degenerate, LAPACK returns an ARBITRARY orthonormal basis
   ! of the degenerate subspace, in which the velocity operator is generally NOT
   ! diagonal. Interband position elements are then built as
   !     r_nm = -i v_nm / (E_n - E_m)                 [get_ome_sp_xme_ex_band]
   ! i.e. a nonzero, basis-arbitrary numerator over a vanishing denominator ->
   ! unbounded noise. Sipe & Shkrebtii's prescription removes the ambiguity:
   !   "The bands at such a degeneracy point should be chosen such that the
   !    v^b_ps(k)=0 for p/=s at that point; then the r^b_ps(k) vanish at such a
   !    point for p/=s."
   ! and Eq. (A4): label the bands so that  E(t).v_mn(k) = 0  for n/=m, i.e. the
   ! basis a k.p analysis would select on moving away from the degeneracy along
   ! the field direction. apply_a4_rotation below implements exactly this.
   !
   ! CAVEAT (stated, not hidden): the prescription is DIRECTION DEPENDENT (it is
   ! E.v that is diagonalised), and v^x, v^y do not commute within a multiplet in
   ! general, so no single basis satisfies it for every Cartesian component at
   ! once. a4_dir selects the direction that is made exact. For D3h monolayer
   ! TMDs the independent SHG component is sigma_xxx (with
   ! sigma_xxx = -sigma_xyy = -sigma_yyx = -sigma_yxy), so a4_dir = 1 makes the
   ! component of interest exact, and the C3 relations then serve as a check.
   !
   ! a4_eps is the multiplet-grouping threshold. Unlike the eps_deg/clip_threshold
   ! cutoffs elsewhere -- which DISCARD data and whose results therefore depend
   ! strongly on the cutoff (HANDOFF 8.21: a 40x change over one decade) -- this
   ! one only selects which bands get re-labelled, discarding nothing, so results
   ! should be INSENSITIVE to it over a broad range. That insensitivity is the
   ! validation test for this code path.
   ! a4_mode selects how bands are grouped into multiplets:
   !   'fixed'    -- consecutive gap < a4_eps (a hard energy threshold, in Hartree). Material
   !                 dependent in the bad sense: 1e-3 Ha = 27 meV is numerical noise next to a
   !                 4 eV TMD bandwidth but a real physical scale in a 27-band GeS window, where
   !                 it grouped genuinely distinct bands and moved sigma_xxx by 86x (HANDOFF 8.26).
   !   'adaptive' -- group whenever the amplification |v_nm| / |E_n - E_m| exceeds a4_rmax, i.e.
   !                 exactly when the position element r_nm = -i v_nm/dE that this pair would
   !                 produce is larger than the code already considers trustworthy. This is
   !                 material-adaptive by construction (it scales with the material's own |v|,
   !                 hence with its hopping) and needs no energy scale to be chosen by hand.
   !                 a4_rmax is deliberately the SAME position scale as clip_threshold, which is
   !                 the value at which get_berry_eigen_fourpoint already discards berry_eigen.
   character(len=8), save :: a4_mode = 'fixed'   
   logical, save     :: a4_enabled = .true.
   integer, save     :: a4_dir     = 1          ! Cartesian direction made exact (1=x,2=y,3=z)
   real(8), save     :: a4_eps     = 1.0d-4     ! 'fixed' mode threshold, Hartree (see HANDOFF 8.26)
   real(8), save     :: a4_rmax    = clip_threshold  ! 'adaptive' mode: max trusted |r| (bohr)

   ! Per-k rotation matrix, restricted to the exciton band window, saved so that Xatu's exciton
   ! envelopes can be carried into the SAME basis (see rotate_fk_ex_to_a4_basis). Without this the
   ! single-particle quantities are rotated while fk_ex is not, which moves the excitonic linear
   ! response peak by 0.32 eV and changes the spectrum by 42-124% (HANDOFF 8.29).
   complex*16, allocatable, save :: a4_W(:,:,:)      ! (npointstotal, nband_ex, nband_ex)
   logical, save :: a4_W_ready = .false.
   logical, save :: a4_W_straddle_warned = .false.
   ! ------------------------------------------------------------------------
   
      real(8), allocatable, save :: Rx_global(:), Ry_global(:), Rz_global(:)
      logical, save :: R_cache_set = .false.


contains
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_R_cache()
            implicit none
            integer :: iRp
            if (.not. active_flags_set) call set_active_flags()
                  allocate(Rx_global(nR), Ry_global(nR), Rz_global(nR))
            do iRp = 1, nR
                  if (active_x) then
                  Rx_global(iRp) = dble(nRvec(iRp,1))*R(1,1) + dble(nRvec(iRp,2))*R(2,1) + dble(nRvec(iRp,3))*R(3,1)
                  else
                  Rx_global(iRp) = 0.0d0
                  end if
                  if (active_y) then
                  Ry_global(iRp) = dble(nRvec(iRp,1))*R(1,2) + dble(nRvec(iRp,2))*R(2,2) + dble(nRvec(iRp,3))*R(3,2)
                  else
                  Ry_global(iRp) = 0.0d0
                  end if
                  if (active_z) then
                  Rz_global(iRp) = dble(nRvec(iRp,1))*R(1,3) + dble(nRvec(iRp,2))*R(2,3) + dble(nRvec(iRp,3))*R(3,3)
                  else
                  Rz_global(iRp) = 0.0d0
                  end if
            end do
            R_cache_set = .true.
      end subroutine set_R_cache
      
   subroutine set_active_flags()
      implicit none
      active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
      active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
      active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
      active_flags_set = .true.
   end subroutine set_active_flags
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> Builds the single-particle optical matrix elements on the k-mesh and writes them
!! to ome_{linear,nonlinear}_sp_<material>.omesp for later runs to read back.
!! For iflag_norder = 2 it also produces the Berry connections, the shift vector and
!! the derivative of |v|, i.e. everything the second-order kernels consume.
!! @param iflag_norder  1 = linear (energies and v only), 2 = second order.
!! @return void
subroutine get_ome_sp(iflag_norder)
      implicit none

      integer iflag_norder
      integer ibz
      integer i,j,ii,jj,nj

      complex*16, allocatable :: skernel(:,:), hkernel(:,:)
      complex*16, allocatable :: sderkernel(:,:,:), hderkernel(:,:,:)
      complex*16, allocatable :: akernel(:,:,:)

      real(8)    :: e(norb)
      complex*16 :: hk_ev(norb,norb)
      complex*16 :: vme(norb,norb,3)

      complex*16 :: abc(norb,norb,3)
      real(8)    :: shift_vector(norb,norb,3,3)
      complex*16 :: berry_eigen1(norb,norb,3), berry_eigen2(norb,norb,3)
      complex*16 :: berry_eigen(norb,norb,3)

      integer :: kmoment

      complex*16, allocatable :: gen_der(:,:,:,:), vme_der(:,:,:,:), vme_der_pt(:,:,:,:)
      real(8),    allocatable :: vme_abs_der(:,:,:,:)
      complex*16, allocatable :: gd1(:,:,:,:), gd2(:,:,:,:), gd3(:,:,:,:)
      complex*16, allocatable :: hk_ev_neigh(:,:,:), vme_neigh(:,:,:,:)

      ! NEW: scratch for the zgemm-based basis change inside get_vme_eigen_ome.
      ! Allocated once per thread (same lifetime/placement as gd1,gd2,gd3),
      ! reused across every call to get_vme_eigen_ome made by this thread —
      ! both the "main" call below and the 7 calls inside
      ! get_berry_eigen_fourpoint.
      complex*16, allocatable :: M1(:,:), T1(:,:)
      complex*16 :: Wfull(norb,norb)   ! total basis change at this k (Eq. A4 rotation + rephasing)
      complex*16 :: Wblk_chk(nband_ex,nband_ex)

      real(8), allocatable :: vme_der_phase(:,:,:,:)

      real(8),    allocatable :: ek(:,:)
      complex*16, allocatable :: vme_ex_band(:,:,:,:)
      complex*16, allocatable :: berry_eigen_ex_band(:,:,:,:)
      complex*16, allocatable :: gen_der_ex_band(:,:,:,:,:)
      real(8),    allocatable :: shift_vector_ex_band(:,:,:,:,:)
      real(8),    allocatable :: vme_abs_der_ex_band(:,:,:,:,:)   ! (ibz, deriv a, pol c, i, j): d_a |v^c_{ij}|
      complex*16, allocatable :: vme_der_pt_ex_band(:,:,:,:,:)    ! (ibz, deriv a, pol c, i, j): (v^c_ij);k^a, gauge-fixed

      real(8) :: rkx,rky,rkz
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '5. Entering ome_sp'

      call set_active_flags()
      if (.not. R_cache_set) call set_R_cache()

      allocate(vme_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(ek(npointstotal,nband_ex))
      allocate(berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(vme_abs_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(vme_der_pt_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      if (allocated(a4_W)) deallocate(a4_W)
      allocate(a4_W(npointstotal,nband_ex,nband_ex)); a4_W=(0.0d0,0.0d0)
      vme_abs_der_ex_band=0.0d0
      vme_der_pt_ex_band=0.0d0
      gen_der_ex_band=0.0d0
      shift_vector_ex_band=0.0d0
      berry_eigen_ex_band=0.0d0
      vme_ex_band=0.0d0
      ek=0.0d0

      write(*,*) '   Calculating optical matrix elements (sp): sampling BZ...'

      kmoment = -1

      ! DEFAULT(NONE) added 2026-09-24 (audit): without it an unlisted variable silently becomes
      ! SHARED and the compiler says nothing. That is exactly how Wfull/Wblk_chk became a data race
      ! earlier in this session -- the symptom appeared three tests away, as wrong excitonic-IPA
      ! scales in check_sp_shift. With DEFAULT(NONE) any future omission is a compile error.
      !$OMP PARALLEL DEFAULT(NONE) PRIVATE(rkx,rky,rkz,ibz,i,j,ii,jj,nj), &
      !$OMP PRIVATE(hkernel,skernel,sderkernel,hderkernel,akernel), &
      !$OMP PRIVATE(hk_ev,e,vme), &
      !$OMP PRIVATE(abc,gen_der,gd1,gd2,gd3), &
      !$OMP PRIVATE(vme_der,vme_abs_der,shift_vector,berry_eigen1,berry_eigen2,berry_eigen), &
      !$OMP PRIVATE(hk_ev_neigh,vme_neigh,vme_der_phase,vme_der_pt)&
      !$OMP PRIVATE(M1,T1)&                                    ! NEW
      !$OMP PRIVATE(Wfull,Wblk_chk)&                           ! NEW: per-k basis change (Eq. A4)
      !$OMP   SHARED(kmoment), &
      ! read-only inputs:
      !$OMP   SHARED(norb,npointstotal,nband_ex,nband_index,iflag_norder), &
      !$OMP   SHARED(rkxvector,rkyvector,rkzvector), &
      ! outputs, written only at the disjoint index ibz owned by this thread:
      !$OMP   SHARED(ek,vme_ex_band,berry_eigen_ex_band,gen_der_ex_band), &
      !$OMP   SHARED(shift_vector_ex_band,vme_abs_der_ex_band,vme_der_pt_ex_band,a4_W), &
      ! double-checked flag: unguarded read is benign, the write is inside !$omp critical:
      !$OMP   SHARED(a4_W_straddle_warned)

      allocate(skernel(norb,norb), hkernel(norb,norb))
      allocate(sderkernel(norb,norb,3), hderkernel(norb,norb,3))
      allocate(akernel(norb,norb,3))
      allocate(gen_der(norb,norb,3,3), vme_der(norb,norb,3,3), vme_abs_der(norb,norb,3,3))
      allocate(vme_der_pt(norb,norb,3,3))
      allocate(gd1(norb,norb,3,3), gd2(norb,norb,3,3), gd3(norb,norb,3,3))
      allocate(hk_ev_neigh(norb,norb,7), vme_neigh(norb,norb,3,7))
      allocate(vme_der_phase(norb,norb,3,3))
      allocate(M1(norb,norb), T1(norb,norb))                    ! NEW

      !$OMP DO SCHEDULE(STATIC)
      do ibz=1,npointstotal
            ! Progress, not a transcript -- see the matching note in ome_ex.f90. This loop printed one
            ! line per k-point, 2025 of them on a 45x45 run, and it is also the loop that emits the
            ! Eq. (A4) straddle WARNING below: a real diagnostic was being hidden inside its own
            ! progress output.
            if (mod(ibz, max(1, npointstotal/10)) == 0) &
              write(*,*) '   Optical matrix elements (sp): k-point',ibz,'/',npointstotal
            rkx=rkxvector(ibz)
            rky=rkyvector(ibz)
            rkz=rkzvector(ibz)

            call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
                  hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
                  hk_ev,e,vme,M1,T1,Wfull)                       ! CHANGED: M1,T1,Wfull appended

            ! Store the exciton-window block of the basis change so fk_ex can follow (8.29).
            do j = 1, nband_ex
               do i = 1, nband_ex
                  a4_W(ibz,i,j) = Wfull(nband_index(i), nband_index(j))
               end do
            end do
            ! A multiplet straddling the window edge would make that block non-unitary, and the
            ! envelope transformation incomplete. Detect it once rather than fail silently.
            if (.not. a4_W_straddle_warned) then
               Wblk_chk = a4_W(ibz,:,:)   ! contiguous copy: avoids an array temporary at the call
               if (a4_block_is_nonunitary(nband_ex, Wblk_chk)) then
                  !$omp critical
                  if (.not. a4_W_straddle_warned) then
                     write(*,*) '   WARNING (ome_sp): an Eq. (A4) multiplet straddles the Bandlist edge;'
                     write(*,*) '            the exciton-window block of the rotation is not unitary and'
                     write(*,*) '            the fk_ex transformation is incomplete at some k points.'
                     write(*,*) '            Widen Bandlist so degenerate partners are inside the window.'
                     a4_W_straddle_warned = .true.
                  end if
                  !$omp end critical
               end if
            end if

            if (iflag_norder.eq.2) then
                  call get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
                  call get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der,vme_abs_der, &
                        shift_vector,berry_eigen1,berry_eigen2,berry_eigen, &
                        hk_ev_neigh,vme_neigh, &
                        skernel,hkernel,sderkernel,hderkernel,akernel,vme_der_phase, &
                        M1,T1,vme_der_pt)                                ! CHANGED: M1,T1,vme_der_pt appended
            end if

            do i=1,nband_ex
                  ii=nband_index(i)
                  ek(ibz,i)=e(ii)
                  do nj=1,3
                        do j=1,nband_ex
                              jj=nband_index(j)
                              vme_ex_band(ibz,nj,i,j)=vme(ii,jj,nj)
                              if (iflag_norder.eq.2) then
                                    shift_vector_ex_band(ibz,nj,1,i,j)=shift_vector(ii,jj,nj,1)
                                    shift_vector_ex_band(ibz,nj,2,i,j)=shift_vector(ii,jj,nj,2)
                                    shift_vector_ex_band(ibz,nj,3,i,j)=shift_vector(ii,jj,nj,3)
                                    gen_der_ex_band(ibz,nj,1,i,j)=gen_der(ii,jj,nj,1)
                                    gen_der_ex_band(ibz,nj,2,i,j)=gen_der(ii,jj,nj,2)
                                    gen_der_ex_band(ibz,nj,3,i,j)=gen_der(ii,jj,nj,3)
                                    ! gauge-invariant derivative of the modulus, vme_abs_der(ii,jj,deriv,pol)
                                    vme_abs_der_ex_band(ibz,nj,1,i,j)=vme_abs_der(ii,jj,nj,1)
                                    vme_abs_der_ex_band(ibz,nj,2,i,j)=vme_abs_der(ii,jj,nj,2)
                                    vme_abs_der_ex_band(ibz,nj,3,i,j)=vme_abs_der(ii,jj,nj,3)
                                    berry_eigen_ex_band(ibz,nj,i,j)=berry_eigen(ii,jj,nj)
                                    ! gauge-fixed (parallel-transported) complex derivative, (v^c_ij);k^a
                                    vme_der_pt_ex_band(ibz,nj,1,i,j)=vme_der_pt(ii,jj,nj,1)
                                    vme_der_pt_ex_band(ibz,nj,2,i,j)=vme_der_pt(ii,jj,nj,2)
                                    vme_der_pt_ex_band(ibz,nj,3,i,j)=vme_der_pt(ii,jj,nj,3)
                              end if
                        end do
                  end do
            end do
      end do
      !$OMP END DO

      deallocate(skernel,hkernel,sderkernel,hderkernel,akernel)
      deallocate(gen_der,vme_der,vme_abs_der,gd1,gd2,gd3,hk_ev_neigh,vme_neigh,vme_der_phase)
      deallocate(vme_der_pt)
      deallocate(M1,T1)                                          ! NEW
      !$OMP END PARALLEL

      a4_W_ready = .true.

      write(*,*) '   Writing optical matrix elements (sp) into file'
      if (iflag_norder.eq.1) then
         call write_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
      end if
      if (iflag_norder.eq.2) then
         call write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek, &
            gen_der_ex_band,shift_vector_ex_band,berry_eigen_ex_band,vme_abs_der_ex_band,vme_der_pt_ex_band)
      end if
      write(*,*) '   Optical matrix elements (sp) have been written in file'

      deallocate(vme_ex_band,ek,berry_eigen_ex_band,gen_der_ex_band,shift_vector_ex_band,vme_abs_der_ex_band, &
                 vme_der_pt_ex_band)
end subroutine get_ome_sp
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> Finite-difference k-derivatives at one k-point, using the six nearest neighbours
!! (+-dk along each axis) plus the centre. Returns the Berry connection split into its
!! two pieces -- berry_eigen1 = i<u_n|d u_n> and berry_eigen2 = <n|A|n>, the Wannier-centre
!! part -- the shift vector R^{a,b}_nm, and the derivatives of v.
!! The eigenvector phase from diagoz is arbitrary at each k, so a derivative of the raw
!! complex v is corrupted by gauge jumps; only the MODULUS derivative vme_abs_der is
!! safe without further work, which is why the shift-vector route uses it.
!! @param rkx,rky,rkz  Central k-point, bohr^-1.
!! @param norb         Number of Wannier orbitals.
!! @param vme_der      d/dk of the velocity matrix elements (gauge-jump prone).
!! @param vme_abs_der  d/dk of |v|, gauge invariant.
!! @param shift_vector R^{a,b}_nm, Eq. (10) of Esteve-Paredes et al., npj Comput. Mater. 11, 13 (2025).
!! @return void
subroutine get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der,vme_abs_der, &
      shift_vector,berry_eigen1,berry_eigen2,berry_eigen, &
      hk_ev_neigh,vme_neigh, &
      skernel,hkernel,sderkernel,hderkernel,akernel,vme_der_phase,&
      M1,T1,vme_der_pt)
      implicit none

      integer norb
      integer nn,nnp
      integer ialpha,ialphap
      integer nj,njp

      dimension hkernel(norb,norb),skernel(norb,norb)
      dimension sderkernel(norb,norb,3),hderkernel(norb,norb,3)
      dimension akernel(norb,norb,3)

      dimension e(norb)

      dimension berry_eigen1(norb,norb,3),berry_eigen2(norb,norb,3)
      dimension berry_eigen(norb,norb,3)
      dimension vme_der(norb,norb,3,3)
      dimension vme_abs_der(norb,norb,3,3)
      real*8 vme_abs_der
      dimension vme_der_phase(norb,norb,3,3)

      complex*16 :: hk_ev_neigh(norb,norb,7)
      complex*16 :: vme_neigh(norb,norb,3,7)
      complex*16 :: M1(norb,norb), T1(norb,norb)                 ! NEW dummy args

      ! NEW: parallel-transported (gauge-fixed) direct generalized derivative of v, (v^c_nm);k^a with
      ! index order (nn,nnp,deriv a,pol c), matching vme_der's own convention. Unlike vme_der above (a raw
      ! finite difference of vme_neigh from INDEPENDENTLY diagonalized neighbouring k-points, corrupted by
      ! gauge jumps -- see HANDOFF.md 8.8), vme_neigh is first phase-aligned per band to the central point's
      ! eigenvectors (a per-k-point, per-band U(1) parallel transport) before differencing. Needed for the
      ! single-particle SHG generalized-derivative term (get_shg_intens_sp in sigma_second_sp.f90), which
      ! needs the full COMPLEX derivative, not just |v|'s (vme_abs_der), so the shift-vector-style
      ! modulus-only workaround does not apply here.
      complex*16 :: vme_der_pt(norb,norb,3,3)
      complex*16, allocatable :: vme_neigh_pt(:,:,:,:)
      complex*16, allocatable :: overlap(:), phase_corr(:)
      integer :: ineigh

      dimension shift_vector(norb,norb,3,3)

      real*8 rkx,rky,rkz,rkx_neigh,rky_neigh,rkz_neigh
      real*8 e
      real*8 vme_der_phase
      real*8 shift_vector
      real*8 ph1,ph2,ph3,ph4,ph5,ph6

      complex*16 hkernel,akernel,skernel,sderkernel,hderkernel
      complex*16 aux1,aux2,aux3,aux4,aux5,aux6
      complex*16 vme_der
      complex*16 berry_eigen1,berry_eigen2,berry_eigen
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! active_x/y/z are module-level, computed once in set_active_flags().

      vme_der=0.0d0
      vme_abs_der=0.0d0
      shift_vector=0.0d0
      vme_der_phase=0.0d0
      hk_ev_neigh=0.0d0
      vme_neigh=0.0d0
      berry_eigen1=0.0d0
      berry_eigen2=0.0d0
      berry_eigen=0.0d0
      vme_der_pt=0.0d0
      
      ! --- central point (index 7) — MOVED to the top of the routine ---
      call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(:,:,7),e,vme_neigh(:,:,:,7),M1,T1)

      ! --- x-direction neighbours ---
      if (active_x) then
            rkx_neigh=rkx-dk; rky_neigh=rky; rkz_neigh=rkz
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,1),e,vme_neigh(:,:,:,1),M1,T1)

            rkx_neigh=rkx+dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,3),e,vme_neigh(:,:,:,3),M1,T1)
      else
            hk_ev_neigh(:,:,1)   = hk_ev_neigh(:,:,7)     ! reuse central point, NO recomputation
            vme_neigh(:,:,:,1)   = vme_neigh(:,:,:,7)
            hk_ev_neigh(:,:,3)   = hk_ev_neigh(:,:,7)
            vme_neigh(:,:,:,3)   = vme_neigh(:,:,:,7)
      end if

      ! --- y-direction neighbours ---
      if (active_y) then
            rkx_neigh=rkx; rky_neigh=rky-dk; rkz_neigh=rkz
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,2),e,vme_neigh(:,:,:,2),M1,T1)

            rky_neigh=rky+dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,4),e,vme_neigh(:,:,:,4),M1,T1)
      else
            hk_ev_neigh(:,:,2)   = hk_ev_neigh(:,:,7)     ! reuse central point, NO recomputation
            vme_neigh(:,:,:,2)   = vme_neigh(:,:,:,7)
            hk_ev_neigh(:,:,4)   = hk_ev_neigh(:,:,7)
            vme_neigh(:,:,:,4)   = vme_neigh(:,:,:,7)
      end if
      
      
      ! --- z-direction neighbours ---
         if (active_z) then
            rkx_neigh=rkx; rky_neigh=rky; rkz_neigh=rkz-dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,5),e,vme_neigh(:,:,:,5),M1,T1)

            rkz_neigh=rkz+dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,6),e,vme_neigh(:,:,:,6),M1,T1)
      else
            hk_ev_neigh(:,:,5)   = hk_ev_neigh(:,:,7)     ! reuse central point, NO recomputation
            vme_neigh(:,:,:,5)   = vme_neigh(:,:,:,7)
            hk_ev_neigh(:,:,6)   = hk_ev_neigh(:,:,7)
            vme_neigh(:,:,:,6)   = vme_neigh(:,:,:,7)
      end if

      ! --- parallel transport: phase-align each neighbour's eigenvectors, band by band, to the central
      ! point (index 7) before they are used in any derivative. overlap(n) = <psi_n(k0)|psi_n(kneigh)>;
      ! dividing by its own phase removes the arbitrary, independent-diagonalization phase jump (the
      ! source of the gauge-jump garbage in vme_der above), leaving only the smooth physical variation.
      ! Not applied to hk_ev_neigh/vme_neigh themselves (berry_eigen/shift_vector already validated as is);
      ! only vme_neigh_pt (used for vme_der_pt) gets it.
      allocate(vme_neigh_pt(norb,norb,3,7))
      allocate(overlap(norb), phase_corr(norb))
      vme_neigh_pt(:,:,:,7) = vme_neigh(:,:,:,7)
      do ineigh = 1, 6
         do nn = 1, norb
            overlap(nn) = sum(conjg(hk_ev_neigh(:,nn,7)) * hk_ev_neigh(:,nn,ineigh))
            if (abs(overlap(nn)) > 1.0d-12) then
               phase_corr(nn) = conjg(overlap(nn)) / abs(overlap(nn))
            else
               phase_corr(nn) = (1.0d0, 0.0d0)     ! near-orthogonal (degenerate/ill-conditioned): no correction possible
            end if
         end do
         do nn = 1, norb
            do nnp = 1, norb
               vme_neigh_pt(nn,nnp,:,ineigh) = conjg(phase_corr(nn)) * phase_corr(nnp) * vme_neigh(nn,nnp,:,ineigh)
            end do
         end do
      end do
      deallocate(overlap, phase_corr)

      do nn=1,norb
         do nnp=1,norb
            do ialpha=1,norb
               do ialphap=1,norb
                  !x-dir
                  if (active_x) then
                  aux1=(hk_ev_neigh(ialphap,nnp,3)-hk_ev_neigh(ialphap,nnp,1))/(2.0d0*dk)
                  berry_eigen1(nn,nnp,1)=berry_eigen1(nn,nnp,1)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(ialpha,nn,7))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(nn,nnp,1)=berry_eigen2(nn,nnp,1)+ &
                  conjg(hk_ev_neigh(ialpha,nn,7))*hk_ev_neigh(ialphap,nnp,7)*akernel(ialpha,ialphap,1)

                  !y-dir
                  if (active_y) then
                  aux1=(hk_ev_neigh(ialphap,nnp,4)-hk_ev_neigh(ialphap,nnp,2))/(2.0d0*dk)
                  berry_eigen1(nn,nnp,2)=berry_eigen1(nn,nnp,2)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(ialpha,nn,7))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(nn,nnp,2)=berry_eigen2(nn,nnp,2)+ &
                  conjg(hk_ev_neigh(ialpha,nn,7))*hk_ev_neigh(ialphap,nnp,7)*akernel(ialpha,ialphap,2)

                  !z-dir
                  if (active_z) then
                  aux1=(hk_ev_neigh(ialphap,nnp,6)-hk_ev_neigh(ialphap,nnp,5))/(2.0d0*dk)
                  berry_eigen1(nn,nnp,3)=berry_eigen1(nn,nnp,3)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(ialpha,nn,7))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(nn,nnp,3)=berry_eigen2(nn,nnp,3)+ &
                  conjg(hk_ev_neigh(ialpha,nn,7))*hk_ev_neigh(ialphap,nnp,7)*akernel(ialpha,ialphap,3)
               end do
            end do

            do nj=1,3
                  berry_eigen(nn,nnp,nj)=berry_eigen1(nn,nnp,nj)+berry_eigen2(nn,nnp,nj)
                  if (abs(berry_eigen(nn,nnp,nj)).gt.clip_threshold) then
                        berry_eigen(nn,nnp,nj)=0.0d0
                  end if

                  aux1=vme_neigh(nn,nnp,nj,1); aux3=vme_neigh(nn,nnp,nj,3)
                  aux2=vme_neigh(nn,nnp,nj,2); aux4=vme_neigh(nn,nnp,nj,4)
                  aux5=vme_neigh(nn,nnp,nj,5); aux6=vme_neigh(nn,nnp,nj,6)

                  if (active_x) then
                        vme_der(nn,nnp,1,nj)=(aux3-aux1)/(2.0d0*dk)
                        vme_abs_der(nn,nnp,1,nj)=(abs(aux3)-abs(aux1))/(2.0d0*dk)   ! gauge invariant
                  end if
                  if (active_y) then
                        vme_der(nn,nnp,2,nj)=(aux4-aux2)/(2.0d0*dk)
                        vme_abs_der(nn,nnp,2,nj)=(abs(aux4)-abs(aux2))/(2.0d0*dk)
                  end if
                  if (active_z) then
                        vme_der(nn,nnp,3,nj)=(aux6-aux5)/(2.0d0*dk)
                        vme_abs_der(nn,nnp,3,nj)=(abs(aux6)-abs(aux5))/(2.0d0*dk)
                  end if

                  ! parallel-transported (gauge-fixed) raw derivative; the -i*(xi_nn-xi_mm)*v gauge-covariant
                  ! correction is added below, once berry_eigen's diagonal is available (same point as vme_der's).
                  if (active_x) vme_der_pt(nn,nnp,1,nj) = (vme_neigh_pt(nn,nnp,nj,3)-vme_neigh_pt(nn,nnp,nj,1))/(2.0d0*dk)
                  if (active_y) vme_der_pt(nn,nnp,2,nj) = (vme_neigh_pt(nn,nnp,nj,4)-vme_neigh_pt(nn,nnp,nj,2))/(2.0d0*dk)
                  if (active_z) vme_der_pt(nn,nnp,3,nj) = (vme_neigh_pt(nn,nnp,nj,6)-vme_neigh_pt(nn,nnp,nj,5))/(2.0d0*dk)

                  call get_phase(vme_neigh(nn,nnp,nj,1),ph1)
                  call get_phase(vme_neigh(nn,nnp,nj,3),ph3)
                  call get_phase(vme_neigh(nn,nnp,nj,2),ph2)
                  call get_phase(vme_neigh(nn,nnp,nj,4),ph4)
                  call get_phase(vme_neigh(nn,nnp,nj,5),ph5)
                  call get_phase(vme_neigh(nn,nnp,nj,6),ph6)

                  if (active_x) then
                        vme_der_phase(nn,nnp,1,nj)=(ph3-ph1)/(2.0d0*dk)
                  end if
                  if (active_y) then
                        vme_der_phase(nn,nnp,2,nj)=(ph4-ph2)/(2.0d0*dk)
                  end if
                  if (active_z) then
                        vme_der_phase(nn,nnp,3,nj)=(ph6-ph5)/(2.0d0*dk)
                  end if
            end do
         end do
      end do

      !shift vector
      do nn=1,norb
         do nnp=1,norb
            do nj=1,3
               do njp=1,3
                  ! Shift vector R^{a,b}_nm = -(d_a theta_b) + (xi^a_nn - xi^a_mm), Nagaosa
                  ! (PRX 10.1103/PhysRevX.10.041041), same convention as its consumer in
                  ! sigma_second_sp.f90 (get_shift_intens_sp), which cites this equation at its
                  ! point of use; cited here too at its point of definition.
                  shift_vector(nn,nnp,nj,njp)=-vme_der_phase(nn,nnp,nj,njp) &
                        +(realpart(berry_eigen(nn,nn,nj))-realpart(berry_eigen(nnp,nnp,nj)))
                  if (abs(shift_vector(nn,nnp,nj,njp)).gt.clip_threshold) then
                        shift_vector(nn,nnp,nj,njp)=0.0d0
                  end if
                  vme_der(nn,nnp,nj,njp)=vme_der(nn,nnp,nj,njp) &
                        -complex(0.0d0,1.0d0)*vme_neigh(nn,nnp,njp,7) &
                        *(realpart(berry_eigen(nn,nn,nj))-realpart(berry_eigen(nnp,nnp,nj)))

                  ! same gauge-covariant correction (Eq. 6b), using the phase-aligned vme_neigh_pt at the
                  ! central point (identical to vme_neigh there, PT is a no-op at dk=0).
                  ! CORRECTED 2026-09-24: use berry_eigen2 (the Wannier-centre part <n|A|n>) ONLY, not the
                  ! full berry_eigen = berry_eigen1 + berry_eigen2. berry_eigen1 is i<u_n|d u_n> evaluated in
                  ! the RAW eigenvector gauge, but the derivative above is taken in the PARALLEL-TRANSPORTED
                  ! gauge, where that part of the connection is zero by construction. Adding it back
                  ! double-counted it and left the generalized derivative non-covariant: hBN's sp SHG kept a
                  ! 20% C3 (D3h) violation that did NOT converge with the k-grid (20.03/20.07/20.07% at
                  ! 30/60/120). With berry_eigen2 alone the violation collapses to ~1e-3 and does converge.
                  vme_der_pt(nn,nnp,nj,njp)=vme_der_pt(nn,nnp,nj,njp) &
                        -complex(0.0d0,1.0d0)*vme_neigh_pt(nn,nnp,njp,7) &
                        *(realpart(berry_eigen2(nn,nn,nj))-realpart(berry_eigen2(nnp,nnp,nj)))
                  ! NO magnitude clip here (removed 2026-09-23, code review finding #6): unlike
                  ! berry_eigen/shift_vector, clip_threshold=50 is the wrong kind of guard for this
                  ! quantity. Calibrated empirically on GeS's full 27-band window (2.5M band-pair
                  ! samples): vme_der_pt is intrinsically heavy-tailed even for clearly non-degenerate
                  ! pairs (gap >= 50 meV: p99=76, p99.9=2235, max=2.0e5 -- real curvature, not garbage),
                  ! so a fixed threshold of 50 wrongly zeroed 1.19% of legitimate far-from-degenerate
                  ! pairs while catching at most 0.04% of genuinely near-degenerate ones. vme_der_pt's
                  ! only consumer (get_shg_intens_sp, term 2) already drops any pair with
                  ! |E_n-E_m| < eps_deg (2.7 meV) before ever using this array -- a strictly better,
                  ! energy-gap-based criterion for "is this pair near-degenerate" than a magnitude clip.
               end do
            end do
         end do
      end do

      deallocate(vme_neigh_pt)

end subroutine get_berry_eigen_fourpoint
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_phase(aux1,ph)
      real*8 ph
      complex*16 aux1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (abs(aux1).lt.1.0d-05) then
         ph=0.0d0
      else
         ph=aimag(log(aux1))
      end if
   end subroutine get_phase
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! PATCH: gd1, gd2, gd3 are now dummy arguments again (caller-owned,
   ! allocated once per thread), NOT locals declared inside this
   ! subroutine -- see get_ome_sp for the reasoning.
   subroutine get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
      implicit none

      integer norb,norb_inter_cut
      integer nn,nnp,nnpp
      integer nj,njp

      dimension e(norb)
      dimension vme(norb,norb,3)
      dimension gen_der(norb,norb,3,3)
      dimension gd1(norb,norb,3,3)
      dimension gd2(norb,norb,3,3)
      dimension gd3(norb,norb,3,3)
      dimension abc(norb,norb,3)

      real*8 e
      complex*16 vme,gen_der,abc,gd1,gd2,gd3
      
      real*8 :: de_nnnp, inv_de_nnnp
      complex*16 :: ci_over_de
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      abc=0.0d0
      gd1=0.0d0
      gd2=0.0d0
      gd3=0.0d0
      gen_der=0.0d0
      norb_inter_cut=norb
      do nn=1,norb
         do nnp=1,norb
            do nj=1,3
               if (abs(e(nn)-e(nnp)).lt.1.0d-05) then
                  abc(nn,nnp,nj)=0.0d0
               else
                  abc(nn,nnp,nj)=-complex(0.0d0,1.0d0)*vme(nn,nnp,nj)/(e(nn)-e(nnp))
               end if
               do njp=1,3
                  de_nnnp = e(nn) - e(nnp)
                  if (abs(de_nnnp).lt.1.0d-05) then
                     gd1(nn,nnp,nj,njp)=0.0d0
                     gd2(nn,nnp,nj,njp)=0.0d0
                     gd3(nn,nnp,nj,njp)=0.0d0
                  else
                     inv_de_nnnp = 1.0d0 / de_nnnp
                     ci_over_de   = complex(0.0d0,1.0d0) * inv_de_nnnp
                     gd1(nn,nnp,nj,njp)=complex(0.0d0,1.0d0)*(inv_de_nnnp**2)* &
                        (vme(nn,nnp,nj)*(vme(nn,nn,njp)-vme(nnp,nnp,njp))+ &
                        vme(nn,nnp,njp)*(vme(nn,nn,nj)-vme(nnp,nnp,nj)))

                     gd2(nn,nnp,nj,njp)=0.0d0
                     gd3(nn,nnp,nj,njp)=0.0d0
                     do nnpp=1,norb_inter_cut
                        if (abs(e(nnpp)-e(nn)).lt.1.0d-05 .or. abs(e(nnpp)-e(nnp)).lt.1.0d-05) then
                              cycle
                        else
                              gd2(nn,nnp,nj,njp)=gd2(nn,nnp,nj,njp)+ &
                              ci_over_de*(vme(nn,nnpp,nj)*vme(nnpp,nnp,njp)/(e(nnpp)-e(nnp)))
                              gd3(nn,nnp,nj,njp)=gd3(nn,nnp,nj,njp)+ &
                              ci_over_de*(-vme(nn,nnpp,njp)*vme(nnpp,nnp,nj)/(e(nn)-e(nnpp)))
                        end if
                     end do
                  end if
                  gen_der(nn,nnp,nj,njp)=gd1(nn,nnp,nj,njp)+gd2(nn,nnp,nj,njp) &
                     +gd3(nn,nnp,nj,njp)
               end do
            end do
         end do
      end do

   end subroutine get_gen_der_sumrule
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
      implicit none
      integer :: iounit10
      integer iflag_norder
      integer npointstotal,nband_ex
      integer ibz
      integer nj,i,j

      dimension ek(npointstotal,nband_ex)
      dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)

      real*8 ek
      complex*16 vme_ex_band
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(newunit=iounit10,file='ome_linear_sp_'//trim(material_name)//'.omesp')
      write(iounit10,*) iflag_norder
      do ibz=1,npointstotal
         write(iounit10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(ek(ibz,j),j=1,nband_ex)
         do i=1,nband_ex
            do j=1,nband_ex
               write(iounit10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz), &
                  (realpart(vme_ex_band(ibz,nj,i,j)),aimag(vme_ex_band(ibz,nj,i,j)), nj=1,3)
            end do
         end do
      end do
      close(iounit10)
   end subroutine write_ome_sp_linear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek, &
      gen_der_ex_band,shift_vector_ex_band,berry_eigen_ex_band,vme_abs_der_ex_band,vme_der_pt_ex_band)
      implicit none
      integer :: iounit10
      integer iflag_norder,npointstotal,nband_ex,ibz

      dimension ek(npointstotal,nband_ex)
      dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
      dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
      dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
      dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)

      real*8 ek, shift_vector_ex_band
      complex*16 vme_ex_band, berry_eigen_ex_band, gen_der_ex_band
      real*8 vme_abs_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
      complex*16 vme_der_pt_ex_band(npointstotal,3,3,nband_ex,nband_ex)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(newunit=iounit10, file='ome_nonlinear_sp_'//trim(material_name)//'.omesp', &
           form='unformatted', access='stream', status='replace')

      write(iounit10) iflag_norder
      write(iounit10) npointstotal, nband_ex
      write(iounit10) rkxvector, rkyvector, rkzvector

      do ibz=1,npointstotal
         write(iounit10) ek(ibz,:)
         write(iounit10) vme_ex_band(ibz,:,:,:)
         write(iounit10) berry_eigen_ex_band(ibz,:,:,:)
         write(iounit10) shift_vector_ex_band(ibz,:,:,:,:)
         write(iounit10) gen_der_ex_band(ibz,:,:,:,:)
      end do
      ! appended AFTER the per-k records so that files written before these arrays existed stay readable
      write(iounit10) vme_abs_der_ex_band
      write(iounit10) vme_der_pt_ex_band

      close(iounit10)
   end subroutine write_ome_sp_nonlinear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
      hkernel,hderkernel,akernel)
      implicit none

      integer norb
      integer ialpha
      integer ialphap
      integer iRp
      integer nj

      dimension skernel(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(norb,norb,3)
      dimension hderkernel(norb,norb,3)
      dimension akernel(norb,norb,3)

      real(8) Rx,Ry,Rz
      real(8) rkx,rky,rkz

      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 phase,factor
      complex*16 hderhop_x,hderhop_y,hderhop_z
      complex*16 sderhop_x,sderhop_y,sderhop_z
      
      !real(8)    :: Rx_arr(nR), Ry_arr(nR), Rz_arr(nR)
      complex*16 :: factor_arr(nR)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! active_x/y/z now module-level (see set_active_flags),
      ! not recomputed on every call.

      hkernel=0.0d0
      hderkernel=0.0d0
      skernel=0.0d0
      sderkernel=0.0d0
      akernel=0.0d0

      
      do iRp = 1, nR
            factor_arr(iRp) = exp(complex(0.0d0, rkx*Rx_global(iRp) + rky*Ry_global(iRp) + rkz*Rz_global(iRp)))
      end do

      
      do ialpha=1,norb
         do ialphap=1,ialpha
            do iRp=1,nR
               

                  Rx = Rx_global(iRp); 
                  Ry = Ry_global(iRp); 
                  Rz = Rz_global(iRp)   ! still needed below for hderhop_x etc.
                  factor = factor_arr(iRp)
               
                  hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+factor*hhop(iRp,ialpha,ialphap)


                  skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
                  factor*shop(iRp,ialpha,ialphap)

                  hderhop_x=complex(0.0d0,Rx)*hhop(iRp,ialpha,ialphap)
                  hderhop_y=complex(0.0d0,Ry)*hhop(iRp,ialpha,ialphap)
                  hderhop_z=complex(0.0d0,Rz)*hhop(iRp,ialpha,ialphap)

                  sderhop_x=complex(0.0d0,Rx)*shop(iRp,ialpha,ialphap)
                  sderhop_y=complex(0.0d0,Ry)*shop(iRp,ialpha,ialphap)
                  sderhop_z=complex(0.0d0,Rz)*shop(iRp,ialpha,ialphap)

                  sderkernel(ialpha,ialphap,1)=sderkernel(ialpha,ialphap,1)+factor*sderhop_x
                  sderkernel(ialpha,ialphap,2)=sderkernel(ialpha,ialphap,2)+factor*sderhop_y
                  sderkernel(ialpha,ialphap,3)=sderkernel(ialpha,ialphap,3)+factor*sderhop_z

                  hderkernel(ialpha,ialphap,1)=hderkernel(ialpha,ialphap,1)+factor*hderhop_x
                  hderkernel(ialpha,ialphap,2)=hderkernel(ialpha,ialphap,2)+factor*hderhop_y
                  hderkernel(ialpha,ialphap,3)=hderkernel(ialpha,ialphap,3)+factor*hderhop_z

                  akernel(ialpha,ialphap,1)=akernel(ialpha,ialphap,1)+ &
                        factor*(rhop_c(1,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_x)
                  akernel(ialpha,ialphap,2)=akernel(ialpha,ialphap,2)+ &
                        factor*(rhop_c(2,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_y)
                  akernel(ialpha,ialphap,3)=akernel(ialpha,ialphap,3)+ &
                        factor*(rhop_c(3,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_z)
            end do

            if (ialphap < ialpha) then
               do nj=1,3
                  hkernel(ialphap,ialpha)=conjg(hkernel(ialpha,ialphap))
                  skernel(ialphap,ialpha)=conjg(skernel(ialpha,ialphap))
                  
                  sderkernel(ialphap,ialpha,nj)=conjg(sderkernel(ialpha,ialphap,nj))
                  hderkernel(ialphap,ialpha,nj)=conjg(hderkernel(ialpha,ialphap,nj))
                  akernel(ialphap,ialpha,nj)=conjg(akernel(ialpha,ialphap,nj))+ &
                     complex(0.0d0,1.0d0)*conjg(sderkernel(ialpha,ialphap,nj))
               end do
            end if

         end do
      end do
   end subroutine get_vme_kernels_ome
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    subroutine get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
!       hk_ev,e,vme)
!       implicit none
! 
!       integer :: norb
!       integer :: ialpha
!       integer :: ialphap
!       integer :: nj
!       integer :: nn,nnp
! 
!       dimension skernel(norb,norb)
!       dimension hkernel(norb,norb)
!       dimension sderkernel(3,norb,norb)
!       dimension hderkernel(3,norb,norb)
!       dimension akernel(3,norb,norb)
! 
!       dimension vjseudoa(3,norb,norb)
!       dimension vjseudob(3,norb,norb)
!       dimension e(norb)
!       dimension hk_ev(norb,norb)
!       dimension vme(3,norb,norb)
! 
!       real*8 e
!       complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
!       complex*16 hk_ev,vjseudoa,vjseudob,vme
!       complex*16 amu,amup
! 
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       e=0.0d0
!       hk_ev(:,:)=hkernel(:,:)
!       vme=0.0d0
!       vjseudoa=0.0d0
!       vjseudob=0.0d0
      
      !vjseudoa(nj,:,:) = U^H · hderkernel(nj,:,:) · U
      !vjseudob(nj,nn,nnp) = i·[ e(nn)·T1(nj,nn,nnp) − e(nnp)·conjg(T1(nj,nnp,nn)) ]
!       do nn=1,norb
!          do nnp=1,nn
!             do ialpha=1,norb
!                do ialphap=1,norb
!                   amu=hk_ev(ialpha,nn)
!                   amup=hk_ev(ialphap,nnp)
!                   do nj=1,3
!                      vjseudoa(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+ &
!                         conjg(amu)*amup*hderkernel(nj,ialpha,ialphap)
!                      vjseudob(nj,nn,nnp)=vjseudob(nj,nn,nnp)+conjg(amu)*amup* &
!                         (e(nn)*akernel(nj,ialpha,ialphap)-e(nnp)*conjg(akernel(nj,ialphap,ialpha)))* &
!                         complex(0.0d0,1.0d0)
!                   end do
!                end do
!             end do
!             do nj=1,3
!                vme(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+vjseudob(nj,nn,nnp)
!                if (nnp < nn) then
!                   vme(nj,nnp,nn)=conjg(vme(nj,nn,nnp))
!                   vjseudoa(nj,nnp,nn)=conjg(vjseudoa(nj,nn,nnp))
!                   vjseudob(nj,nnp,nn)=conjg(vjseudob(nj,nn,nnp))
!                end if
!             end do
!          end do
!       end do
!    end subroutine get_vme_eigen_ome

subroutine get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
      hk_ev,e,vme,M1,T1,Wout)
      implicit none

      complex*16, intent(out), optional :: Wout(norb,norb)
      complex*16 :: Wsave(norb,norb)
      integer :: norb, nj, nn, nnp

      dimension skernel(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(norb,norb,3)      ! <-- reordered
      dimension hderkernel(norb,norb,3)      ! <-- reordered
      dimension akernel(norb,norb,3)         ! <-- reordered
      dimension e(norb)
      dimension hk_ev(norb,norb)
      dimension vme(norb,norb,3)             ! <-- reordered

      real*8 e
      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 hk_ev,vme
      complex*16 M1(norb,norb), T1(norb,norb)
      complex*16, parameter :: cone=(1.0d0,0.0d0), czero=(0.0d0,0.0d0), ci=(0.0d0,1.0d0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      e=0.0d0
      call diagoz(norb,e,hkernel)
      hk_ev(:,:)=hkernel(:,:)
      call phase_eigvec_nk(norb,hk_ev)

      call build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)

      ! Sipe & Shkrebtii Eq. (A4): re-label near-degenerate multiplets so that the
      ! velocity along a4_dir is diagonal within each one, then rebuild vme in that
      ! basis. No-op when no multiplet is within a4_eps (e.g. hBN's 2-band model).
      ! Capture the TOTAL basis change, not just the multiplet rotation: phase_eigvec_nk runs again
      ! after the rotation and multiplies each column by a phase, so hk_ev_new = hk_ev_old * (W*D).
      ! Rather than track W and D separately, project: with orthonormal Wannier functions the
      ! eigenvectors satisfy hk_ev^H hk_ev = I, so hk_ev_old^H hk_ev_new IS the change-of-basis
      ! matrix exactly, phases included.
      if (present(Wout)) then
         Wsave = hk_ev
         if (apply_a4_rotation(norb,e,hk_ev,vme)) then
            call phase_eigvec_nk(norb,hk_ev)
            call build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)
         end if
         call zgemm('C','N',norb,norb,norb, cone, Wsave, norb, hk_ev, norb, czero, Wout, norb)
      else
         if (apply_a4_rotation(norb,e,hk_ev,vme)) then
            call phase_eigvec_nk(norb,hk_ev)
            call build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)
         end if
      end if

   end subroutine get_vme_eigen_ome
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! vme in the eigenbasis: U^H (dH/dk) U plus the Berry-connection correction
   ! i[e_n A_nm - e_m A_mn^*] (Esteve-Paredes & Palacios, SciPost Phys. Core 6,
   ! 002 (2023), Eq. 16). Factored out of get_vme_eigen_ome so it can be re-run
   ! after the Eq. (A4) multiplet rotation without duplicating the kernel.
   subroutine build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)
      implicit none
      integer,    intent(in)    :: norb
      complex*16, intent(in)    :: hderkernel(norb,norb,3), akernel(norb,norb,3)
      complex*16, intent(in)    :: hk_ev(norb,norb)
      real(8),    intent(in)    :: e(norb)
      complex*16, intent(inout) :: M1(norb,norb), T1(norb,norb)
      complex*16, intent(out)   :: vme(norb,norb,3)
      integer :: nj, nn, nnp
      complex*16, parameter :: cone=(1.0d0,0.0d0), czero=(0.0d0,0.0d0), ci=(0.0d0,1.0d0)

      vme=0.0d0
      do nj=1,3
         ! Both slices below are genuinely contiguous (norb*norb contiguous
         ! elements each) — no compiler-inserted temporary, at compile time OR
         ! at runtime under -fcheck=array-temps.
         call zgemm('N','N',norb,norb,norb, cone, hderkernel(:,:,nj), norb, hk_ev, norb, czero, M1, norb)
         call zgemm('C','N',norb,norb,norb, cone, hk_ev, norb, M1, norb, czero, vme(:,:,nj), norb)

         call zgemm('N','N',norb,norb,norb, cone, akernel(:,:,nj), norb, hk_ev, norb, czero, M1, norb)
         call zgemm('C','N',norb,norb,norb, cone, hk_ev, norb, M1, norb, czero, T1, norb)

         do nnp=1,norb
            do nn=1,norb
               vme(nn,nnp,nj) = vme(nn,nnp,nj) + ci*( e(nn)*T1(nn,nnp) - e(nnp)*conjg(T1(nnp,nn)) )
            end do
         end do
      end do
   end subroutine build_vme_blocks
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Sipe & Shkrebtii PRB 61, 5337 (2000), Eq. (A4): within each near-degenerate
   ! multiplet, re-label (rotate) the bands so that the velocity operator along
   ! a4_dir is diagonal, i.e. v_ps = 0 for p /= s inside the multiplet. Then the
   ! interband position r_ps = -i v_ps / (E_p - E_s) has a vanishing numerator
   ! there instead of an arbitrary one, which is what makes it well defined.
   !
   ! Returns .true. if any multiplet was found (so the caller knows vme must be
   ! rebuilt from the rotated eigenvectors). Rotating hk_ev and recomputing vme
   ! is used rather than transforming vme in place, because the eigenvector phase
   ! convention (phase_eigvec_nk) has to be re-imposed after the rotation and it
   ! is far less error-prone to redo the two zgemms than to track the phases
   ! through vme by hand. A per-band U(1) rephasing maps v_ps -> e^{-i th_p}
   ! e^{i th_s} v_ps, so it preserves the zeros the rotation just created.
   logical function apply_a4_rotation(norb, e, hk_ev, vme, Wout) result(rotated)
      implicit none
      integer,    intent(in)    :: norb
      real(8),    intent(in)    :: e(norb)
      complex*16, intent(inout) :: hk_ev(norb,norb)
      complex*16, intent(in)    :: vme(norb,norb,3)
      ! Optional: the accumulated basis change, identity outside multiplets, such that
      ! hk_ev_new = hk_ev_old * Wout. Needed to carry fk_ex into the same basis.
      complex*16, intent(out), optional :: Wout(norb,norb)

      integer :: nlo, nhi, nd, i, j
      complex*16, allocatable :: Vblk(:,:), col(:,:)
      real(8),    allocatable :: wblk(:)

      rotated = .false.
      if (present(Wout)) then
         Wout = (0.0d0,0.0d0)
         do i = 1, norb
            Wout(i,i) = (1.0d0,0.0d0)
         end do
      end if
      if (.not. a4_enabled) return

      nlo = 1
      do while (nlo <= norb)
         ! Extend the multiplet while consecutive bands count as degenerate by the selected
         ! criterion. 'fixed': a hard energy gap (material-dependent in the bad sense).
         ! 'adaptive': the pair would produce |r| = |v|/|dE| above a4_rmax, i.e. larger than
         ! the code already trusts anywhere else -- no hand-chosen energy scale, and it
         ! follows the material's own velocity scale automatically.
         nhi = nlo
         do while (nhi < norb)
            if (.not. pair_is_degenerate(e(nhi+1)-e(nhi), vme(nhi,nhi+1,a4_dir))) exit
            nhi = nhi + 1
         end do
         nd = nhi - nlo + 1

         if (nd > 1) then
            allocate(Vblk(nd,nd), wblk(nd), col(norb,nd))
            ! velocity block along a4_dir, explicitly Hermitised before zheev
            do j = 1, nd
               do i = 1, nd
                  Vblk(i,j) = 0.5d0*( vme(nlo+i-1, nlo+j-1, a4_dir) &
                                    + conjg(vme(nlo+j-1, nlo+i-1, a4_dir)) )
               end do
            end do
            call diagoz(nd, wblk, Vblk)      ! Vblk <- eigenvectors W of the block
            ! hk_ev(:, nlo:nhi) <- hk_ev(:, nlo:nhi) * W
            col = hk_ev(:, nlo:nhi)
            hk_ev(:, nlo:nhi) = matmul(col, Vblk)
            if (present(Wout)) Wout(nlo:nhi, nlo:nhi) = Vblk
            deallocate(Vblk, wblk, col)
            rotated = .true.
         end if

         nlo = nhi + 1
      end do
   end function apply_a4_rotation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Degeneracy test for one adjacent band pair. See a4_mode at the top of the module.
   logical function pair_is_degenerate(de, v) result(deg)
      implicit none
      real(8),    intent(in) :: de      ! E_{n+1} - E_n (Hartree)
      complex*16, intent(in) :: v       ! v_{n,n+1} along a4_dir
      if (a4_mode == 'fixed') then
         deg = (abs(de) < a4_eps)
      else
         ! |r| = |v|/|dE| > a4_rmax  <=>  |v| > a4_rmax*|dE|, written multiplicatively so a
         ! numerically zero gap cannot divide by zero.
         deg = (abs(v) > a4_rmax*abs(de))
      end if
   end function pair_is_degenerate
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Carry Xatu's exciton envelopes into the SAME basis the Eq. (A4) rotation put the
   ! single-particle states in. With |X> = sum_cv A_cv c^dagger_c c_v |0> and the single-particle
   ! states relabelled by |c~> = sum_c |c> W_{c c~}, invariance of |X> requires
   !
   !       A~ = (W_c)^H . A . W_v
   !
   ! with W_c, W_v the conduction- and valence-manifold blocks of W. Xatu and opticx are verified to
   ! start from the SAME basis (HANDOFF 8.29: with the rotation off, opticx reproduces the trusted
   ! MoSe2_ex.dat peak at 1.0705 eV with a scale of 0.959), so opticx's own W is the right matrix to
   ! apply -- Xatu's basis does not have to be discovered.
   !
   ! Must be called after get_ome_sp (which fills a4_W) and before anything consumes fk_ex.
   ! Is the exciton-window block of the rotation non-unitary? True when a multiplet straddles the
   ! Bandlist edge, i.e. a rotated band's degenerate partner lies outside the window.
   logical function a4_block_is_nonunitary(nb, Wblk) result(bad)
      implicit none
      integer,    intent(in) :: nb
      complex*16, intent(in) :: Wblk(nb,nb)
      complex*16 :: G(nb,nb)
      integer :: i, j
      G = matmul(conjg(transpose(Wblk)), Wblk)
      bad = .false.
      do j = 1, nb
         do i = 1, nb
            if (i == j) then
               if (abs(G(i,j)-(1.0d0,0.0d0)) > 1.0d-8) bad = .true.
            else
               if (abs(G(i,j)) > 1.0d-8) bad = .true.
            end if
         end do
      end do
   end function a4_block_is_nonunitary
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine rotate_fk_ex_to_a4_basis()
      use parser_optics_xatu_dim, only: fk_ex, norb_ex_cut, fk_ex_basis_ok, fk_ex_loaded
      implicit none
      integer :: ibz, ic, icp, iv, ivp, n, idx, idxp
      complex*16 :: Wc(nc_ex,nc_ex), Wv(nv_ex,nv_ex)
      complex*16, allocatable :: A(:,:), tmp(:,:)
      real(8) :: dev
      logical :: any_rot

      ! Envelopes not read yet: get_exciton_data deferred them because a second-order cache hit may
      ! make them unnecessary. Return WITHOUT setting fk_ex_basis_ok -- if the cache then misses,
      ! get_ome_ex loads them and calls this routine again, and only that second call may validate
      ! them. Marking them usable here would let an unrotated fk_ex through on a cache miss, which is
      ! exactly the O(1) error the flag exists to prevent (HANDOFF 8.46).
      if (.not. fk_ex_loaded) return

      ! Rotation switched off: the single-particle states were never rotated, so fk_ex as read from
      ! Xatu already matches them and is fit to use.
      if (.not. a4_enabled) then
         fk_ex_basis_ok = .true.
         return
      end if
      ! No rotation matrices: fk_ex CANNOT be carried over, so it stays marked unusable. This is only
      ! a warning and not a stop, because a second-order OME cache hit (HANDOFF 8.44) never reads
      ! fk_ex at all -- get_ome_ex issues the hard error at the point where it commits to the k-loop.
      if (.not. a4_W_ready) then
         write(*,*) '   WARNING (ome_sp): the Eq. (A4) rotation is active but the per-k rotation'
         write(*,*) '            matrices are unavailable (OME_sp = none reads matrix elements from'
         write(*,*) '            file and does not rebuild them). fk_ex CANNOT be carried into the'
         write(*,*) '            rotated basis, so every excitonic quantity would be computed with'
         write(*,*) '            mismatched bases. Regenerate with OME_sp = nonlinear.'
         return
      end if

      allocate(A(nc_ex,nv_ex), tmp(nc_ex,nv_ex))
      any_rot = .false.
      do ibz = 1, npointstotal
         ! blocks of W on the conduction and valence manifolds of the exciton window
         do icp = 1, nc_ex
            do ic = 1, nc_ex
               Wc(ic,icp) = a4_W(ibz, nv_ex+ic, nv_ex+icp)
            end do
         end do
         do ivp = 1, nv_ex
            do iv = 1, nv_ex
               Wv(iv,ivp) = a4_W(ibz, iv, ivp)
            end do
         end do

         ! Skip k points where this is the identity (the common case).
         dev = 0.0d0
         do ic = 1, nc_ex
            do icp = 1, nc_ex
               if (ic == icp) then
                  dev = max(dev, abs(Wc(ic,icp)-(1.0d0,0.0d0)))
               else
                  dev = max(dev, abs(Wc(ic,icp)))
               end if
            end do
         end do
         do iv = 1, nv_ex
            do ivp = 1, nv_ex
               if (iv == ivp) then
                  dev = max(dev, abs(Wv(iv,ivp)-(1.0d0,0.0d0)))
               else
                  dev = max(dev, abs(Wv(iv,ivp)))
               end if
            end do
         end do
         if (dev < 1.0d-12) cycle
         any_rot = .true.

         do n = 1, norb_ex_cut
            do ic = 1, nc_ex
               do iv = 1, nv_ex
                  idx = nc_ex*nv_ex*(ibz-1) + nv_ex*(ic-1) + iv
                  A(ic,iv) = fk_ex(idx, n)
               end do
            end do
            ! tmp = Wc^H . A ; A~ = tmp . Wv
            tmp = matmul(conjg(transpose(Wc)), A)
            A   = matmul(tmp, Wv)
            do ic = 1, nc_ex
               do iv = 1, nv_ex
                  idxp = nc_ex*nv_ex*(ibz-1) + nv_ex*(ic-1) + iv
                  fk_ex(idxp, n) = A(ic,iv)
               end do
            end do
         end do
      end do
      deallocate(A, tmp)
      ! Reached only with a4_W filled by get_ome_sp, so fk_ex is now in the rotated basis (or no k
      ! point needed rotating, which is the same thing).
      fk_ex_basis_ok = .true.
      if (any_rot) write(*,*) '   Exciton envelopes carried into the Eq. (A4) rotated basis'
   end subroutine rotate_fk_ex_to_a4_basis
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine phase_eigvec_nk(norb,hk_ev)
      implicit none

      integer norb
      integer i,j,ii

      dimension hk_ev(norb,norb)

      real*8 :: arg
      complex*16 :: aux1,hk_ev,factor
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,norb
         aux1=0.0d0
         do i=1,norb
            aux1=aux1+hk_ev(i,j)
         end do
         arg=atan2(aimag(aux1),realpart(aux1))
         factor=exp(complex(0.0d0,-arg))
         do ii=1,norb
            hk_ev(ii,j)=hk_ev(ii,j)*factor
         end do
      end do

   end subroutine phase_eigvec_nk
end module ome_sp