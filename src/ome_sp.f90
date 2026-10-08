module ome_sp
   use constants_math
   use parser_input_file, only:iflag_orthonormal, response_text, sp_covariant, nf
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
   ! Response = shift_covariant: U(n)-covariant generalised derivative over blocks of
   ! (near-)degenerate bands. Bands are grouped into blocks at the central k (consecutive gap < cov_tol);
   ! every neighbour eigenbasis is parallel-transported BLOCK-wise (unitary/polar part of the block overlap)
   ! and r^{c;a} = d_a r^c - i [xi^a_B, r^c], with xi_B the intra-block (non-Abelian) Wannier-centre
   ! connection. Built from eigensystems WITHOUT the Eq. (A4) rotation (the lineshape term in
   ! get_shift_intens_cov needs the energy eigenbasis). Results are insensitive to cov_tol over 1-30 meV
   ! (MoS2: identical for 1, 3, 10 and 30 meV).
   real(8), parameter :: cov_tol = 1.1d-4                  ! Hartree (3 meV)
   ! pairs closer than this are degenerate to numerical precision (see get_shift_intens_cov); also the
   ! grouping used by the block scramble test hook
   real(8), parameter :: cov_exact_tol = 1.0d-6           ! Hartree (27 ueV)
   complex*16, allocatable, save :: cov_rgen(:,:,:,:,:)   ! (npointstotal, a, c, i, j): r^{c;a}_ij, window
   complex*16, allocatable, save :: cov_vraw(:,:,:,:)     ! (npointstotal, c, i, j): v^c_ij, eigenbasis, window
   integer,    allocatable, save :: cov_blk(:,:)          ! (npointstotal, i): block label of window band i
   integer, save :: cov_straddle = 0                      ! k-points with a block crossing the Bandlist edge
   ! Response = shg_covariant only, window restricted:
   complex*16, allocatable, save :: cov_nb_T(:,:,:,:)     ! (npointstotal, p=1..6, i, j): neighbour eigenbasis -> centre frame
   complex*16, allocatable, save :: cov_nb_r(:,:,:,:,:)   ! (npointstotal, p=1..7, c, i, j): off-block r in p's eigenbasis
   real(8),    allocatable, save :: cov_nb_e(:,:,:)       ! (npointstotal, p=1..7, i): energies at p
   complex*16, allocatable, save :: cov_xi(:,:,:,:)       ! (npointstotal, c, i, j): Hermitian block connection, centre
   ! BLOCK SCRAMBLE TEST HOOK (OPTICX_BLOCK_SCRAMBLE set): a k-dependent unitary is applied inside every
   ! group of numerically degenerate bands (gap < cov_exact_tol) at the evaluated k-point, right after
   ! diagonalisation, on the shift_covariant path only -- i.e. the eigensolver "returned another basis".
   ! Results must not move (beyond O(splitting)). Never for production.
   logical, save :: block_scramble_on = .false.

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
   ! strongly on the cutoff (a 40x change over one decade, measured) -- this
   ! one only selects which bands get re-labelled, discarding nothing, so results
   ! should be INSENSITIVE to it over a broad range. That insensitivity is the
   ! validation test for this code path.
   ! a4_mode selects how bands are grouped into multiplets:
   !   'fixed'    -- consecutive gap < a4_eps (a hard energy threshold, in Hartree). Material
   !                 dependent in the bad sense: 1e-3 Ha = 27 meV is numerical noise next to a
   !                 4 eV TMD bandwidth but a real physical scale in a 27-band GeS window, where
   !                 it grouped genuinely distinct bands and moved sigma_xxx by 86x.
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
   real(8), save     :: a4_eps     = 1.0d-4     ! 'fixed' mode threshold, Hartree
   real(8), save     :: a4_rmax    = clip_threshold  ! 'adaptive' mode: max trusted |r| (bohr)

   ! Per-k rotation matrix, restricted to the exciton band window, saved so that Xatu's exciton
   ! envelopes can be carried into the SAME basis (see rotate_fk_ex_to_a4_basis). Without this the
   ! single-particle quantities are rotated while fk_ex is not, which moves the excitonic linear
   ! response peak by 0.32 eV and changes the spectrum by 42-124%.
   complex*16, allocatable, save :: a4_W(:,:,:)      ! (npointstotal, nband_ex, nband_ex)
   logical, save :: a4_W_ready = .false.
   logical, save :: a4_W_straddle_warned = .false.
   ! Exciton-window eigenvectors and Wannier-position connection on the k-mesh, for the covariant
   ! derivative of the exciton envelopes. Filled in the same loop and the same basis
   ! as a4_W (after the Eq. (A4) rotation and phase_eigvec_nk), so they describe exactly the states
   ! fk_ex is expressed in once rotate_fk_ex_to_a4_basis has run.
   !   win_evec(:,i,k)  = c_i(k), column nband_index(i) of hk_ev
   !   win_sevec(:,i,k) = S(k) c_i(k)  (S = identity for orthonormal Wannier functions)
   !   win_xi2(k,a,i,j) = c_i(k)^H A^a(k) c_j(k), the Wannier-position part of the Berry connection
   !                      (berry_eigen2 is its diagonal), off-diagonal window elements included
   complex*16, allocatable, save :: win_evec(:,:,:), win_sevec(:,:,:)
   complex*16, allocatable, save :: win_xi2(:,:,:,:)
   logical, save :: win_ready = .false.
   ! Tagged tail of the nonlinear .omesp (format 2, 2026-10-06). After the per-k records and the two
   ! derivative arrays the file holds OMESP_TAG, a flags word and norb, then the optional sections in this
   ! order: a4_W (bit 0), the window states win_evec/win_sevec/win_xi2 (bit 1), the shift_covariant arrays
   ! (bit 2), the covariant second-order arrays (bit 3). Storing a4_W and the window states lets a run with
   ! OME_sp = none rebuild the basis the exciton envelopes must be carried into, which it otherwise cannot.
   ! Files without the tag (older versions) are still read; they simply lack those sections.
   integer, parameter :: OMESP_TAG = 1330464562
   ! GAUGE TEST HOOK, for tools/check_gauge_covariance.py only. With the environment
   ! variable OPTICX_GAUGE_TEST set (non-empty, not '0'), every eigenvector is multiplied by a smooth
   ! but rapidly varying phase exp(i theta_n(k)) after the Eq. (A4) rotation and phase_eigvec_nk, and
   ! vme is rebuilt in that basis; the phase is also captured in Wout, so fk_ex follows it. Physical
   ! results must not move. theta changes by ~1 rad between neighbours of a 30x30 hBN mesh but is
   ! smooth on the dk = 1e-6 scale of the four-point derivatives (|grad theta| <= 30 bohr, below
   ! clip_threshold). gt_k is the k of the last get_vme_kernels_ome call on this thread, which is
   ! always the k of the get_vme_eigen_ome call that follows it.
   logical, save :: gauge_test_on = .false.
   real(8), save :: gt_amp = 1.5d0                 ! amplitude of theta; OPTICX_GAUGE_TEST=x scales it by x
   real(8), save :: gt_k(3) = 0.0d0
   !$omp threadprivate(gt_k)
   ! ------------------------------------------------------------------------
   
      real(8), allocatable, save :: Rx_global(:), Ry_global(:), Rz_global(:)
      logical, save :: R_cache_set = .false.

   ! Bloch sums as matrix products: the LOWER triangle (a >= b) of H(R), S(R) and r(R), packed once by
   ! set_R_cache into (nR, ntri) matrices (r: (nR, 3*ntri)), so that the sums over R for several k-points
   ! are one zgemm (get_vme_kernels_batch). Reading the hoppings once per stencil instead of once per
   ! k-point is what matters: the sums are memory bound (MoS2, 34 orbitals, 333 cells: 31 MB per k).
   integer, save :: ntri = 0
   integer, allocatable, save :: tri_a(:), tri_b(:)
   complex*16, allocatable, save :: hl_pack(:,:), sl_pack(:,:), rl_pack(:,:)
   ! Per-thread stencil cache: the kernels at a k-point and its +-dk neighbours, filled by
   ! fill_kernel_cache with one batched product and served by get_vme_kernels_ome on an EXACT k match
   ! (the centre is requested 4 times and each neighbour twice per k-point of a covariant run). A miss
   ! computes the k-point directly, so the cache can never return kernels of another k.
   integer, save :: kc_n = 0
   real(8), save :: kc_k(3,7)
   complex*16, allocatable, save :: kc_s(:,:,:), kc_h(:,:,:), kc_sd(:,:,:,:), kc_hd(:,:,:,:), kc_a(:,:,:,:)
   !$omp threadprivate(kc_n, kc_k, kc_s, kc_h, kc_sd, kc_hd, kc_a)


contains
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_R_cache()
            implicit none
            integer :: iRp, t, a, b, nj
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
            ! lower triangle of the hoppings, packed column by column (b outer, a = b..norb inner)
            if (allocated(hl_pack)) deallocate(hl_pack, sl_pack, rl_pack, tri_a, tri_b)
            ntri = norb*(norb+1)/2
            allocate(hl_pack(nR,ntri), sl_pack(nR,ntri), rl_pack(nR,3*ntri), tri_a(ntri), tri_b(ntri))
            t = 0
            do b = 1, norb
               do a = b, norb
                  t = t + 1
                  tri_a(t) = a; tri_b(t) = b
                  hl_pack(:,t) = hhop(:,a,b); sl_pack(:,t) = shop(:,a,b)
                  do nj = 1, 3
                     rl_pack(:,(nj-1)*ntri+t) = rhop_c(nj,:,a,b)
                  end do
               end do
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
      ! Response = shift_covariant: per-thread scratch for get_rgen_blocks
      logical :: cov_on
      complex*16, allocatable :: rgen_k(:,:,:,:), vraw_k(:,:,:)
      integer,    allocatable :: lab_k(:)
      logical :: shgcov_on
      complex*16, allocatable :: Tn_k(:,:,:), rr_k(:,:,:,:), xi_k(:,:,:)
      real(8),    allocatable :: ee_k(:,:)
      integer :: nstr
      ! Filling check: OptiX occupies bands by index, so Nfermi must put the Fermi level in a gap at EVERY
      ! k-point. Highest energy of band nf, lowest of band nf+1, smallest direct gap, over the mesh (Hartree).
      real(8) :: fill_vmax, fill_cmin, fill_gapmin
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '5. Entering ome_sp'
      block
         character(len=16) :: gt_val
         integer :: gt_len
         call get_environment_variable('OPTICX_GAUGE_TEST', gt_val, gt_len)
         gauge_test_on = (gt_len > 0 .and. trim(gt_val) /= '0')
         if (gauge_test_on) then
            gt_amp = 1.5d0
            block
               real(8) :: scale
               integer :: ios
               read(gt_val, *, iostat=ios) scale
               if (ios == 0) gt_amp = 1.5d0*scale
            end block
         end if
         if (gauge_test_on) write(*,*) '   GAUGE TEST (OPTICX_GAUGE_TEST set): eigenvector phases are being'// &
            ' scrambled on purpose; results must be unchanged. Never use for production.'
      end block

      call set_active_flags()
      if (.not. R_cache_set) call set_R_cache()

      ! Block-covariant data, stored in the .omesp for ANY second-order Response so that one file serves every
      ! branch with OME_sp = none: cov_rgen/cov_vraw/cov_blk for the shift current (shift_covariant) and the
      ! neighbour data cov_nb_*/cov_xi for the covariant method-B kernel (shg, electrooptic, rectification,
      ! general). Both come from the same get_rgen_blocks call per k-point.
      shgcov_on = (iflag_norder == 2)
      cov_on = (iflag_norder == 2)
      block
         character(len=16) :: bs_val
         integer :: bs_len
         call get_environment_variable('OPTICX_BLOCK_SCRAMBLE', bs_val, bs_len)
         block_scramble_on = (bs_len > 0 .and. trim(bs_val) /= '0')
      end block
      if (block_scramble_on .and. cov_on) write(*,*) '   BLOCK SCRAMBLE TEST (OPTICX_BLOCK_SCRAMBLE set): states inside'// &
         ' exactly degenerate groups are rotated on purpose; results must be unchanged. Never use for production.'
      if (allocated(cov_rgen)) deallocate(cov_rgen)
      if (allocated(cov_vraw)) deallocate(cov_vraw)
      if (allocated(cov_blk))  deallocate(cov_blk)
      if (allocated(cov_nb_T)) deallocate(cov_nb_T)
      if (allocated(cov_nb_r)) deallocate(cov_nb_r)
      if (allocated(cov_nb_e)) deallocate(cov_nb_e)
      if (allocated(cov_xi))   deallocate(cov_xi)
      cov_straddle = 0
      if (cov_on) then
         allocate(cov_rgen(npointstotal,3,3,nband_ex,nband_ex), cov_vraw(npointstotal,3,nband_ex,nband_ex))
         allocate(cov_blk(npointstotal,nband_ex))
         cov_rgen = (0.0d0,0.0d0); cov_vraw = (0.0d0,0.0d0); cov_blk = 0
      end if
      if (shgcov_on) then
         allocate(cov_nb_T(npointstotal,6,nband_ex,nband_ex), cov_nb_r(npointstotal,7,3,nband_ex,nband_ex))
         allocate(cov_nb_e(npointstotal,7,nband_ex), cov_xi(npointstotal,3,nband_ex,nband_ex))
         cov_nb_T = (0.0d0,0.0d0); cov_nb_r = (0.0d0,0.0d0); cov_nb_e = 0.0d0; cov_xi = (0.0d0,0.0d0)
      end if
      nstr = 0
      fill_vmax = -huge(1.0d0); fill_cmin = huge(1.0d0); fill_gapmin = huge(1.0d0)

      allocate(vme_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(ek(npointstotal,nband_ex))
      allocate(berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(vme_abs_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(vme_der_pt_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      if (allocated(a4_W)) deallocate(a4_W)
      allocate(a4_W(npointstotal,nband_ex,nband_ex)); a4_W=(0.0d0,0.0d0)
      win_ready = .false.
      if (allocated(win_evec))  deallocate(win_evec)
      if (allocated(win_sevec)) deallocate(win_sevec)
      if (allocated(win_xi2))   deallocate(win_xi2)
      if (iflag_norder.eq.2) then
         allocate(win_evec(norb,nband_ex,npointstotal), win_sevec(norb,nband_ex,npointstotal))
         allocate(win_xi2(npointstotal,3,nband_ex,nband_ex))
         win_evec=(0.0d0,0.0d0); win_sevec=(0.0d0,0.0d0); win_xi2=(0.0d0,0.0d0)
      end if
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
      ! SHARED and the compiler says nothing. That is exactly how Wfull/Wblk_chk once became a data
      ! race -- the symptom appeared three tests away, as wrong excitonic-IPA scales in check_sp_shift.
      ! With DEFAULT(NONE) any future omission is a compile error.
      !$OMP PARALLEL DEFAULT(NONE) PRIVATE(rkx,rky,rkz,ibz,i,j,ii,jj,nj), &
      !$OMP PRIVATE(hkernel,skernel,sderkernel,hderkernel,akernel), &
      !$OMP PRIVATE(hk_ev,e,vme), &
      !$OMP PRIVATE(abc,gen_der,gd1,gd2,gd3), &
      !$OMP PRIVATE(vme_der,vme_abs_der,shift_vector,berry_eigen1,berry_eigen2,berry_eigen), &
      !$OMP PRIVATE(hk_ev_neigh,vme_neigh,vme_der_phase,vme_der_pt)&
      !$OMP PRIVATE(M1,T1)&                                    ! NEW
      !$OMP PRIVATE(Wfull,Wblk_chk)&                           ! NEW: per-k basis change (Eq. A4)
      !$OMP PRIVATE(rgen_k,vraw_k,lab_k)&                      ! shift_covariant scratch
      !$OMP PRIVATE(Tn_k,rr_k,ee_k,xi_k)&                      ! shg_covariant scratch
      !$OMP   SHARED(cov_on,cov_rgen,cov_vraw,cov_blk) REDUCTION(+:nstr) &
      !$OMP   REDUCTION(max:fill_vmax) REDUCTION(min:fill_cmin,fill_gapmin) SHARED(nf) &
      !$OMP   SHARED(shgcov_on,cov_nb_T,cov_nb_r,cov_nb_e,cov_xi) &
      !$OMP   SHARED(kmoment), &
      ! read-only inputs:
      !$OMP   SHARED(norb,npointstotal,nband_ex,nband_index,iflag_norder), &
      !$OMP   SHARED(rkxvector,rkyvector,rkzvector), &
      ! outputs, written only at the disjoint index ibz owned by this thread:
      !$OMP   SHARED(ek,vme_ex_band,berry_eigen_ex_band,gen_der_ex_band), &
      !$OMP   SHARED(shift_vector_ex_band,vme_abs_der_ex_band,vme_der_pt_ex_band,a4_W), &
      !$OMP   SHARED(win_evec,win_sevec,win_xi2), &
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
      if (cov_on) allocate(rgen_k(norb,norb,3,3), vraw_k(norb,norb,3), lab_k(norb))
      if (shgcov_on) allocate(Tn_k(norb,norb,6), rr_k(norb,norb,3,7), ee_k(norb,7), xi_k(norb,norb,3))

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

            ! second order: the derivative routines below need this k-point and its +-dk neighbours
            ! several times each; build them in one pass and serve the repeats from the cache
            if (iflag_norder == 2) call fill_kernel_cache(rkx,rky,rkz)
            call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
                  hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
                  hk_ev,e,vme,M1,T1,Wfull)                       ! CHANGED: M1,T1,Wfull appended
            if (nf >= 1 .and. nf < norb) then
               fill_vmax = max(fill_vmax, e(nf)); fill_cmin = min(fill_cmin, e(nf+1))
               fill_gapmin = min(fill_gapmin, e(nf+1) - e(nf))
            end if

            ! Store the exciton-window block of the basis change so fk_ex can follow.
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
                  ! Window states for the covariant envelope derivative. Must come BEFORE
                  ! get_berry_eigen_fourpoint, which reuses skernel/akernel for the neighbour points.
                  do i = 1, nband_ex
                        ii = nband_index(i)
                        win_evec(:,i,ibz)  = hk_ev(:,ii)
                        win_sevec(:,i,ibz) = matmul(skernel, hk_ev(:,ii))
                  end do
                  do nj = 1, 3
                        do j = 1, nband_ex
                              jj = nband_index(j)
                              do i = 1, nband_ex
                                    ii = nband_index(i)
                                    win_xi2(ibz,nj,i,j) = dot_product(hk_ev(:,ii), matmul(akernel(:,:,nj), hk_ev(:,jj)))
                              end do
                        end do
                  end do
                  if (cov_on) then
                        ! own kernel calls: must come before get_berry_eigen_fourpoint like the window block above
                        if (shgcov_on) then
                              call get_rgen_blocks(rkx,rky,rkz,norb,rgen_k,vraw_k,lab_k, &
                                    skernel,sderkernel,hkernel,hderkernel,akernel,M1,T1, &
                                    Tn_out=Tn_k,rr_out=rr_k,ee_out=ee_k,xi_out=xi_k)
                              do j = 1, nband_ex
                                    jj = nband_index(j)
                                    cov_nb_e(ibz,:,j) = ee_k(jj,:)
                                    do i = 1, nband_ex
                                          ii = nband_index(i)
                                          cov_nb_T(ibz,:,i,j) = Tn_k(ii,jj,:)
                                          do nj = 1, 3
                                                cov_nb_r(ibz,:,nj,i,j) = rr_k(ii,jj,nj,:)
                                          end do
                                          cov_xi(ibz,:,i,j) = xi_k(ii,jj,:)
                                    end do
                              end do
                        else
                              call get_rgen_blocks(rkx,rky,rkz,norb,rgen_k,vraw_k,lab_k, &
                                    skernel,sderkernel,hkernel,hderkernel,akernel,M1,T1)
                        end if
                        do j = 1, nband_ex
                              jj = nband_index(j)
                              cov_blk(ibz,j) = lab_k(jj)
                              do i = 1, nband_ex
                                    ii = nband_index(i)
                                    cov_vraw(ibz,:,i,j) = vraw_k(ii,jj,:)
                                    cov_rgen(ibz,:,:,i,j) = rgen_k(ii,jj,:,:)
                              end do
                        end do
                        ! a block with members both inside and outside the window breaks the covariance
                        ! of the window sum: count such k-points (reported after the loop)
                        if (any([(count(lab_k == lab_k(nband_index(j))) /= &
                                  count(cov_blk(ibz,:) == lab_k(nband_index(j))), j = 1, nband_ex)])) nstr = nstr + 1
                        ! get_rgen_blocks overwrote the kernels with the last neighbour's; restore the centre (a copy from the stencil cache)
                        call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
                  end if
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
      if (cov_on) deallocate(rgen_k, vraw_k, lab_k)
      if (shgcov_on) deallocate(Tn_k, rr_k, ee_k, xi_k)
      !$OMP END PARALLEL
      ! A filling that is not insulating gives transitions of vanishing energy, which the second-order
      ! formulas turn into 1/omega^2 poles (GeS at Nfermi = 21: sigma^xyy 135 -> 10566), and an ill-defined
      ! occupation for the linear response. Warn; the input may be intentional (e.g. a deliberate test).
      if (nf >= 1 .and. nf < norb) then
         if (fill_vmax > fill_cmin .or. fill_gapmin < 1.0d-3/27.211385d0) then
            write(*,'(a,i0,a)') '    WARNING (ome_sp): Nfermi = ', nf, ' does not put the Fermi level in a gap on this mesh.'
            write(*,'(a,i0,a,f10.4,a,i0,a,f10.4,a)') '             band ', nf, ' reaches ', fill_vmax*27.211385d0, &
                 ' eV and band ', nf+1, ' comes down to ', fill_cmin*27.211385d0, ' eV'
            write(*,'(a,f10.4,a,f10.4,a)') '             (indirect gap ', (fill_cmin - fill_vmax)*27.211385d0, &
                 ' eV, smallest direct gap ', fill_gapmin*27.211385d0, ' eV).'
            write(*,'(a)') '             OptiX occupies bands by index, so this filling describes a metal or semimetal;'
            write(*,'(a)') '             second-order responses then contain poles at vanishing transition energy. Check Nfermi.'
         end if
      end if
      if (cov_on) then
         cov_straddle = nstr
         if (nstr > 0) then
            write(*,*) '   WARNING (shift_covariant): at',nstr,'k-points a block of degenerate bands crosses the'
            write(*,*) '            Bandlist edge; the window sum is not gauge covariant there. Widen Bandlist.'
         end if
      end if

      a4_W_ready = .true.
      if (iflag_norder.eq.2) win_ready = .true.

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
      ! gauge jumps: on GeS 2% of its entries exceed 50, max 5e4), vme_neigh is first phase-aligned per band to the central point's
      ! eigenvectors (a per-k-point, per-band U(1) parallel transport) before differencing. Needed for the
      ! single-particle SHG generalized-derivative term (get_shg_intens_sp in sigma_second_sp.f90), which
      ! needs the full COMPLEX derivative, not just |v|'s (vme_abs_der), so the shift-vector-style
      ! modulus-only workaround does not apply here.
      complex*16 :: vme_der_pt(norb,norb,3,3)
      complex*16, allocatable :: vme_neigh_pt(:,:,:,:)
      complex*16, allocatable :: overlap(:), phase_corr(:)
      integer :: ineigh
      ! Overlap matrix S(k0) at the CENTRAL point: the neighbour kernel calls below overwrite skernel, and in
      ! a non-orthonormal basis both the parallel-transport overlap and berry_eigen1 need S at k0.
      complex*16 :: s_center(norb,norb)
      complex*16, allocatable :: sc_center(:,:)
      complex*16, allocatable :: bwork(:,:), bdiff(:,:)
      complex*16, parameter :: zone=(1.0d0,0.0d0), zzero=(0.0d0,0.0d0), zi=(0.0d0,1.0d0)
      real*8, allocatable :: vme_der_phase_pt(:,:,:,:)   ! phase derivative of vme_neigh_pt (shift vector)

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
      s_center = skernel
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
      ! Not applied to hk_ev_neigh/vme_neigh themselves (berry_eigen and vme_der still use the raw gauge);
      ! only vme_neigh_pt gets it, used for vme_der_pt and, since 2026-10-05, the shift vector.
      allocate(vme_neigh_pt(norb,norb,3,7))
      allocate(overlap(norb), phase_corr(norb))
      if (.not. iflag_orthonormal) then
         allocate(sc_center(norb,norb))
         sc_center = matmul(s_center, hk_ev_neigh(:,:,7))
      end if
      vme_neigh_pt(:,:,:,7) = vme_neigh(:,:,:,7)
      do ineigh = 1, 6
         do nn = 1, norb
            ! <psi_n(k0)|psi_n(k')> = c_n(k0)^H S(k0) c_n(k'). Without S the transported phase, i.e. the Berry
            ! connection, is wrong in a non-orthonormal basis: it put the sp SHG 24% off on an on-site
            ! non-orthogonal rewrite of hBN (same physics). S = 1 for orthonormal Wannier functions.
            if (iflag_orthonormal) then
               overlap(nn) = sum(conjg(hk_ev_neigh(:,nn,7)) * hk_ev_neigh(:,nn,ineigh))
            else
               overlap(nn) = sum(conjg(sc_center(:,nn)) * hk_ev_neigh(:,nn,ineigh))
            end if
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
      if (allocated(sc_center)) deallocate(sc_center)

      ! Berry connection as matrix products (was a norb^4 scalar loop, ~1e9 iterations per k at norb = 144,
      ! essentially the whole cost of get_ome_sp for large models):
      !   berry_eigen1(:,:,a) = i C0^H S0 D_a,   D_a = central difference of the neighbour eigenvectors
      !   berry_eigen2(:,:,a) = C0^H A_a C0
      ! with C0 = hk_ev_neigh(:,:,7) and S0 = s_center. akernel is used exactly as before (the kernel left by
      ! the last neighbour call, O(dk) from the central one), so results only change at round-off.
      allocate(bwork(norb,norb), bdiff(norb,norb))
      do nj = 1, 3
         call zgemm('N','N',norb,norb,norb, zone, akernel(:,:,nj), norb, hk_ev_neigh(:,:,7), norb, zzero, bwork, norb)
         call zgemm('C','N',norb,norb,norb, zone, hk_ev_neigh(:,:,7), norb, bwork, norb, zzero, berry_eigen2(:,:,nj), norb)
      end do
      do nj = 1, 3
         if (nj == 1 .and. .not. active_x) cycle
         if (nj == 2 .and. .not. active_y) cycle
         if (nj == 3 .and. .not. active_z) cycle
         select case (nj)
            case (1); bdiff = (hk_ev_neigh(:,:,3) - hk_ev_neigh(:,:,1))/(2.0d0*dk)
            case (2); bdiff = (hk_ev_neigh(:,:,4) - hk_ev_neigh(:,:,2))/(2.0d0*dk)
            case (3); bdiff = (hk_ev_neigh(:,:,6) - hk_ev_neigh(:,:,5))/(2.0d0*dk)
         end select
         call zgemm('N','N',norb,norb,norb, zone, s_center, norb, bdiff, norb, zzero, bwork, norb)
         call zgemm('C','N',norb,norb,norb, zi, hk_ev_neigh(:,:,7), norb, bwork, norb, zzero, berry_eigen1(:,:,nj), norb)
      end do
      deallocate(bwork, bdiff)
      allocate(vme_der_phase_pt(norb,norb,3,3)); vme_der_phase_pt = 0.0d0

      do nn=1,norb
         do nnp=1,norb

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

                  ! Phase derivative as arg(v(k+dk) v*(k-dk)) / (2 dk), NOT (arg v(k+dk) - arg v(k-dk)) / (2 dk).
                  ! FIXED 2026-10-04: the difference of two principal values jumps by 2 pi
                  ! whenever v straddles the negative real axis, which the real phase convention of
                  ! phase_eigvec_nk makes common (hBN: many v^y elements lie exactly on it). The jump gave
                  ! ~3e6 bohr, the clip below zeroed the shift vector, and that k point's contribution was
                  ! lost: sigma^xyy came out 6.5% low against the exact IPA on hBN 30x30, and moved by 9.1%
                  ! under an arbitrarily small change of eigenvector phases. Same 1e-5 magnitude guard as
                  ! get_phase.
                  if (active_x) vme_der_phase(nn,nnp,1,nj) = &
                        phase_diff(vme_neigh(nn,nnp,nj,3), vme_neigh(nn,nnp,nj,1))/(2.0d0*dk)
                  if (active_y) vme_der_phase(nn,nnp,2,nj) = &
                        phase_diff(vme_neigh(nn,nnp,nj,4), vme_neigh(nn,nnp,nj,2))/(2.0d0*dk)
                  if (active_z) vme_der_phase(nn,nnp,3,nj) = &
                        phase_diff(vme_neigh(nn,nnp,nj,6), vme_neigh(nn,nnp,nj,5))/(2.0d0*dk)
                  ! Same phase derivative from the PARALLEL-TRANSPORTED neighbours; this is the one the
                  ! shift vector uses (see the shift-vector comment below).
                  if (active_x) vme_der_phase_pt(nn,nnp,1,nj) = &
                        phase_diff(vme_neigh_pt(nn,nnp,nj,3), vme_neigh_pt(nn,nnp,nj,1))/(2.0d0*dk)
                  if (active_y) vme_der_phase_pt(nn,nnp,2,nj) = &
                        phase_diff(vme_neigh_pt(nn,nnp,nj,4), vme_neigh_pt(nn,nnp,nj,2))/(2.0d0*dk)
                  if (active_z) vme_der_phase_pt(nn,nnp,3,nj) = &
                        phase_diff(vme_neigh_pt(nn,nnp,nj,6), vme_neigh_pt(nn,nnp,nj,5))/(2.0d0*dk)
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
                  ! CHANGED 2026-10-05: evaluated in the parallel-transported gauge, i.e. the
                  ! phase derivative of vme_neigh_pt plus the Wannier-centre connection berry_eigen2 ONLY
                  ! (in that gauge Re i<u_n|d u_n> = 0 by construction, as for vme_der_pt, where adding it
                  ! anyway left a 20% D3h violation on hBN that did not converge with the k-grid). It used to be the raw-gauge phase derivative plus the FULL berry_eigen,
                  ! whose diagonal carries the raw gauge's phase gradient and is zeroed by the clip at
                  ! |berry_eigen| > clip_threshold: where the clip fired on A_nn but not on the matching term
                  ! of arg(v) the gauge part no longer cancelled and the shift vector became garbage just
                  ! below the clip. On SnTe (C2v) that broke the mirror m_y: odd-y components up to 8% of
                  ! the largest, forbidden fraction 2.89% -> 0.01% after this change. hBN changes by 3e-12.
                  shift_vector(nn,nnp,nj,njp)=-vme_der_phase_pt(nn,nnp,nj,njp) &
                        +(realpart(berry_eigen2(nn,nn,nj))-realpart(berry_eigen2(nnp,nnp,nj)))
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

      deallocate(vme_neigh_pt, vme_der_phase_pt)

end subroutine get_berry_eigen_fourpoint
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> U(n)-covariant generalised derivative of the interband position for Response = shift_covariant
   !! (Esteve-Paredes et al., npj Comput. Mater. 11, 13 (2025), Eq. 9; validated against an independent NumPy
   !! implementation of the same construction to ~1% on MoS2 and against the exact Eq. 9 on hBN to 1e-7).
   !! Seven plain eigensystems (centre + six neighbours at +-dk, NO Eq. (A4) rotation). Blocks = runs of
   !! consecutive bands with gap < cov_tol at the centre. Off-block position r^c_nm = -i v^c_nm/(E_n - E_m)
   !! at every point; each neighbour's r is carried into the centre's frame by the block-diagonal unitary
   !! T_B = Y X^H, X S Y^H = svd(U_B^H S0 U'_B), i.e. U'_B -> U'_B T_B, which makes every block overlap
   !! Hermitian positive (parallel transport). Then r^{c;a} = d_a r^c - i [xi^a, r^c] with xi^a the
   !! Hermitian part of the block-diagonal Wannier-centre connection U^H A_a U at the centre (the transported
   !! part i U^H dU has no Hermitian part in that gauge, as for vme_der_pt).
   !! Covariant under any unitary inside every block, so exactly degenerate bands need no cut.
   !! @param rgen  (norb,norb,a,c): r^{c;a}_nm, zero inside blocks.
   !! @param vraw  (norb,norb,c): velocity in the centre's eigenbasis (no A4).
   !! @param lab   block label of every band at the centre.
   subroutine get_rgen_blocks(rkx,rky,rkz,norb,rgen,vraw,lab,skernel,sderkernel,hkernel,hderkernel,akernel,M1,T1, &
                              Tn_out,rr_out,ee_out,xi_out)
      implicit none
      integer, intent(in) :: norb
      real(8), intent(in) :: rkx,rky,rkz
      complex*16, intent(out) :: rgen(norb,norb,3,3), vraw(norb,norb,3)
      integer, intent(out) :: lab(norb)
      ! Optional, for Response = shg_covariant: per neighbour p (order above) the block transport
      ! eigenbasis -> centre frame (Tn_out), the off-block position in that point's OWN eigenbasis (rr_out, p = 7 is
      ! the centre) and its energies (ee_out), and the Hermitian block-diagonal Wannier-centre connection (xi_out).
      complex*16, intent(out), optional :: Tn_out(norb,norb,6), rr_out(norb,norb,3,7), xi_out(norb,norb,3)
      real(8),    intent(out), optional :: ee_out(norb,7)
      complex*16 :: skernel(norb,norb),sderkernel(norb,norb,3),hkernel(norb,norb),hderkernel(norb,norb,3)
      complex*16 :: akernel(norb,norb,3), M1(norb,norb), T1(norb,norb)
      complex*16, allocatable :: ev(:,:,:), vv(:,:,:,:), rr(:,:,:,:), a_c(:,:,:), xi(:,:,:), O(:,:), T(:,:)
      complex*16, allocatable :: rt(:,:,:,:), tmp(:,:), s_c(:,:)
      real(8),    allocatable :: ee(:,:)
      real(8) :: kp(3), dkv(3,6)
      integer :: p, n, m, a, c, nblk, s, ip
      integer :: bstart(norb), bsize(norb)
      complex*16, parameter :: ci = (0.0d0,1.0d0)
      logical :: act(3)

      allocate(ev(norb,norb,7), vv(norb,norb,3,7), rr(norb,norb,3,7), ee(norb,7))
      allocate(a_c(norb,norb,3), xi(norb,norb,3), O(norb,norb), T(norb,norb), rt(norb,norb,3,2), tmp(norb,norb))
      allocate(s_c(norb,norb))
      act = (/ active_x, active_y, active_z /)
      if (present(Tn_out)) then               ! inactive directions keep the identity
         Tn_out = (0.0d0,0.0d0)
         do n = 1, norb
            Tn_out(n,n,:) = (1.0d0,0.0d0)
         end do
      end if
      ! neighbour order as get_berry_eigen_fourpoint: 1 x-, 3 x+, 2 y-, 4 y+, 5 z-, 6 z+, 7 centre
      dkv = 0.0d0
      dkv(1,1) = -dk; dkv(1,3) = dk; dkv(2,2) = -dk; dkv(2,4) = dk; dkv(3,5) = -dk; dkv(3,6) = dk
      do p = 7, 1, -1                            ! centre first: its kernels are kept (A_c below)
         if (p < 7) then
            a = merge(1, merge(2, 3, p == 2 .or. p == 4), p == 1 .or. p == 3)
            if (.not. act(a)) cycle
            kp = (/ rkx, rky, rkz /) + dkv(:,p)
         else
            kp = (/ rkx, rky, rkz /)
         end if
         call get_vme_kernels_ome(kp(1),kp(2),kp(3),norb,skernel,sderkernel,hkernel,hderkernel,akernel)
         if (p == 7) then
            a_c = akernel; s_c = skernel        ! S(k0): overlap of the block transport (non-orthonormal basis)
         end if
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
               ev(:,:,p),ee(:,p),vv(:,:,:,p),M1,T1,skip_a4=.true.,scramble_blocks=(block_scramble_on .and. p == 7))
      end do

      ! blocks at the centre
      nblk = 1; bstart(1) = 1; lab(1) = 1
      do n = 2, norb
         if (ee(n,7) - ee(n-1,7) >= cov_tol) then
            nblk = nblk + 1; bstart(nblk) = n
         end if
         lab(n) = nblk
      end do
      do s = 1, nblk
         if (s < nblk) then
            bsize(s) = bstart(s+1) - bstart(s)
         else
            bsize(s) = norb - bstart(s) + 1
         end if
      end do

      ! off-block position at every point (centre's partition)
      rr = (0.0d0,0.0d0)
      do p = 1, 7
         if (p < 7) then
            a = merge(1, merge(2, 3, p == 2 .or. p == 4), p == 1 .or. p == 3)
            if (.not. act(a)) cycle
         end if
         do c = 1, 3
            do m = 1, norb
               do n = 1, norb
                  if (lab(n) /= lab(m)) rr(n,m,c,p) = -ci*vv(n,m,c,p)/(ee(n,p) - ee(m,p))
               end do
            end do
         end do
      end do
      vraw = vv(:,:,:,7)

      ! xi^a: Hermitian part of the block-diagonal Wannier-centre connection at the centre
      xi = (0.0d0,0.0d0)
      do a = 1, 3
         tmp = matmul(conjg(transpose(ev(:,:,7))), matmul(a_c(:,:,a), ev(:,:,7)))
         do m = 1, norb
            do n = 1, norb
               if (lab(n) == lab(m)) xi(n,m,a) = 0.5d0*(tmp(n,m) + conjg(tmp(m,n)))
            end do
         end do
      end do

      rgen = (0.0d0,0.0d0)
      do a = 1, 3
         if (.not. act(a)) then
            ! Non-periodic direction (z of a 2D model): there is no k derivative, and the generalised derivative is
            ! the connection part alone, r^{c;a} = -i [xi^a, r^c] (Quintela & Pedersen, PRB 107, 235416 (2023),
            ! Eq. (15)). Skipping it set every shift current ALONG z to zero; with it, buckled hBN's sigma^zxx and
            ! sigma^zzz equal the excitonic shift-current route to 1e-4.
            do c = 1, 3
               rgen(:,:,a,c) = -ci*(matmul(xi(:,:,a), rr(:,:,c,7)) - matmul(rr(:,:,c,7), xi(:,:,a)))
            end do
            cycle
         end if
         do s = 1, 2                                ! s = 1: +dk, s = 2: -dk
            ip = merge(merge(3, 1, s == 1), merge(merge(4, 2, s == 1), merge(6, 5, s == 1), a == 2), a == 1)
            ! <psi_n(k0)|psi_m(k')> = c_n(k0)^H S(k0) c_m(k') (S = 1 for orthonormal Wannier functions)
            if (iflag_orthonormal) then
               O = matmul(conjg(transpose(ev(:,:,7))), ev(:,:,ip))
            else
               O = matmul(conjg(transpose(ev(:,:,7))), matmul(s_c, ev(:,:,ip)))
            end if
            call block_transport(norb, nblk, bstart, bsize, O, T)
            if (present(Tn_out)) Tn_out(:,:,ip) = T
            do c = 1, 3
               tmp = matmul(rr(:,:,c,ip), T)
               rt(:,:,c,s) = matmul(conjg(transpose(T)), tmp)
            end do
         end do
         do c = 1, 3
            tmp = (rt(:,:,c,1) - rt(:,:,c,2))/(2.0d0*dk)
            rgen(:,:,a,c) = tmp - ci*(matmul(xi(:,:,a), rr(:,:,c,7)) - matmul(rr(:,:,c,7), xi(:,:,a)))
         end do
      end do
      if (present(rr_out)) rr_out = rr
      if (present(ee_out)) ee_out = ee
      if (present(xi_out)) xi_out = xi
      ! intra-block elements are not defined by this construction
      do m = 1, norb
         do n = 1, norb
            if (lab(n) == lab(m)) rgen(n,m,:,:) = (0.0d0,0.0d0)
         end do
      end do
      deallocate(ev, vv, rr, ee, a_c, xi, O, T, rt, tmp, s_c)
   end subroutine get_rgen_blocks

   !> Block-diagonal parallel transport: for every block B, T_B = Y X^H from svd(O_B) = X S Y^H, so that
   !! (U_B^H U'_B) T_B = X S X^H is Hermitian positive. 1x1 blocks reduce to the per-band phase conj(o)/|o|.
   subroutine block_transport(norb, nblk, bstart, bsize, O, T)
      implicit none
      integer, intent(in) :: norb, nblk, bstart(norb), bsize(norb)
      complex*16, intent(in) :: O(norb,norb)
      complex*16, intent(out) :: T(norb,norb)
      complex*16, allocatable :: Ob(:,:), Xs(:,:), Yh(:,:), work(:)
      real(8), allocatable :: sv(:), rwork(:)
      integer :: s, nb, i0, info, lwork
      T = (0.0d0,0.0d0)
      do s = 1, nblk
         nb = bsize(s); i0 = bstart(s)
         if (nb == 1) then
            if (abs(O(i0,i0)) > 1.0d-12) then
               T(i0,i0) = conjg(O(i0,i0))/abs(O(i0,i0))
            else
               T(i0,i0) = (1.0d0,0.0d0)
            end if
            cycle
         end if
         allocate(Ob(nb,nb), Xs(nb,nb), Yh(nb,nb), sv(nb), rwork(5*nb), work(1))
         Ob = O(i0:i0+nb-1, i0:i0+nb-1)
         call zgesvd('A','A',nb,nb,Ob,nb,sv,Xs,nb,Yh,nb,work,-1,rwork,info)
         lwork = max(1, int(real(work(1))))
         deallocate(work); allocate(work(lwork))
         Ob = O(i0:i0+nb-1, i0:i0+nb-1)
         call zgesvd('A','A',nb,nb,Ob,nb,sv,Xs,nb,Yh,nb,work,lwork,rwork,info)
         if (info /= 0) then
            write(*,*) 'ERROR (block_transport): zgesvd failed, info =', info
            stop 1
         end if
         T(i0:i0+nb-1, i0:i0+nb-1) = conjg(transpose(matmul(Xs, Yh)))      ! (X Y^H)^H = Y X^H
         deallocate(Ob, Xs, Yh, sv, rwork, work)
      end do
   end subroutine block_transport

   !> arg(a b*), branch-safe phase difference of two nearby complex numbers; 0 if either is below
   !! the 1e-5 magnitude floor get_phase uses.
   pure real(8) function phase_diff(a, b)
      complex*16, intent(in) :: a, b
      if (abs(a) < 1.0d-05 .or. abs(b) < 1.0d-05) then
         phase_diff = 0.0d0
      else
         phase_diff = atan2(aimag(a*conjg(b)), dble(a*conjg(b)))
      end if
   end function phase_diff
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
      ! Appended after the per-k records (older readers stop before it): the Eq. (A4) rotation of the
      ! exciton window, so that OME_sp = none can still carry Xatu's envelopes into this basis.
      if (allocated(a4_W) .and. a4_W_ready) then
         write(iounit10,'(a)') '#A4W'
         write(iounit10,*) a4_W
      end if
      close(iounit10)
   end subroutine write_ome_sp_linear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek, &
      gen_der_ex_band,shift_vector_ex_band,berry_eigen_ex_band,vme_abs_der_ex_band,vme_der_pt_ex_band)
      implicit none
      integer :: iounit10, flags
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
      ! Tagged tail (see OMESP_TAG): flags say which optional sections follow, in a fixed order.
      flags = 0
      if (allocated(a4_W) .and. a4_W_ready) flags = ibset(flags, 0)
      if (allocated(win_evec) .and. win_ready) flags = ibset(flags, 1)
      if (allocated(cov_rgen)) flags = ibset(flags, 2)
      if (allocated(cov_nb_T)) flags = ibset(flags, 3)
      write(iounit10) OMESP_TAG, flags, norb
      if (btest(flags, 0)) write(iounit10) a4_W
      if (btest(flags, 1)) then
         write(iounit10) win_evec
         write(iounit10) win_sevec
         write(iounit10) win_xi2
      end if
      if (btest(flags, 2)) then
         write(iounit10) cov_rgen
         write(iounit10) cov_vraw
         write(iounit10) cov_blk
      end if
      if (btest(flags, 3)) then
         write(iounit10) cov_nb_T
         write(iounit10) cov_nb_r
         write(iounit10) cov_nb_e
         write(iounit10) cov_xi
      end if

      close(iounit10)
   end subroutine write_ome_sp_nonlinear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> S(k), H(k), their k-derivatives and the position kernel A(k) at one k-point. Served from the
   !! per-thread stencil cache (fill_kernel_cache) when k matches a cached point exactly, otherwise computed.
   subroutine get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
      hkernel,hderkernel,akernel)
      implicit none
      integer, intent(in) :: norb
      real(8), intent(in) :: rkx, rky, rkz
      complex*16, intent(out) :: skernel(norb,norb), hkernel(norb,norb)
      complex*16, intent(out) :: sderkernel(norb,norb,3), hderkernel(norb,norb,3), akernel(norb,norb,3)
      complex*16, allocatable :: s1(:,:,:), h1(:,:,:), sd1(:,:,:,:), hd1(:,:,:,:), a1(:,:,:,:)
      real(8) :: kp(3,1)
      integer :: i

      gt_k = (/ rkx, rky, rkz /)       ! gauge test hook: k of the eigenproblem that follows
      do i = 1, kc_n
         if (kc_k(1,i) == rkx .and. kc_k(2,i) == rky .and. kc_k(3,i) == rkz) then
            skernel = kc_s(:,:,i); hkernel = kc_h(:,:,i)
            sderkernel = kc_sd(:,:,:,i); hderkernel = kc_hd(:,:,:,i); akernel = kc_a(:,:,:,i)
            return
         end if
      end do
      allocate(s1(norb,norb,1), h1(norb,norb,1), sd1(norb,norb,3,1), hd1(norb,norb,3,1), a1(norb,norb,3,1))
      kp(:,1) = (/ rkx, rky, rkz /)
      call get_vme_kernels_batch(1, kp, s1, h1, sd1, hd1, a1)
      skernel = s1(:,:,1); hkernel = h1(:,:,1)
      sderkernel = sd1(:,:,:,1); hderkernel = hd1(:,:,:,1); akernel = a1(:,:,:,1)
      deallocate(s1, h1, sd1, hd1, a1)
   end subroutine get_vme_kernels_ome

   !> Fill this thread's stencil cache: the k-point and its +-dk neighbours along the active directions,
   !! with the k values exactly as get_berry_eigen_fourpoint and get_rgen_blocks form them.
   subroutine fill_kernel_cache(rkx, rky, rkz)
      implicit none
      real(8), intent(in) :: rkx, rky, rkz
      if (.not. allocated(kc_s)) allocate(kc_s(norb,norb,7), kc_h(norb,norb,7), kc_sd(norb,norb,3,7), &
                                          kc_hd(norb,norb,3,7), kc_a(norb,norb,3,7))
      kc_n = 1
      kc_k(:,1) = (/ rkx, rky, rkz /)
      if (active_x) then
         kc_k(:,kc_n+1) = (/ rkx-dk, rky, rkz /); kc_k(:,kc_n+2) = (/ rkx+dk, rky, rkz /); kc_n = kc_n + 2
      end if
      if (active_y) then
         kc_k(:,kc_n+1) = (/ rkx, rky-dk, rkz /); kc_k(:,kc_n+2) = (/ rkx, rky+dk, rkz /); kc_n = kc_n + 2
      end if
      if (active_z) then
         kc_k(:,kc_n+1) = (/ rkx, rky, rkz-dk /); kc_k(:,kc_n+2) = (/ rkx, rky, rkz+dk /); kc_n = kc_n + 2
      end if
      call get_vme_kernels_batch(kc_n, kc_k(:,1:kc_n), kc_s, kc_h, kc_sd, kc_hd, kc_a)
   end subroutine fill_kernel_cache

   !> Bloch sums at nk k-points in one pass over the packed hoppings (set_R_cache):
   !! W(R,:) = f, i Rx f, i Ry f, i Rz f with f = e^{i k.R}, so [H; dH/dk] = W^T H(R), [S; dS/dk] = W^T S(R),
   !! sum_R f r(R) = f^T r(R), and A_j = sum_R f r_j(R) + i dS/dk_j. Only the lower triangle is summed; the
   !! upper one is completed from it by Hermitian conjugation, A_ba = conj(A_ab) + i conj(dS_ab/dk).
   subroutine get_vme_kernels_batch(nk, kp, skb, hkb, sdb, hdb, akb)
      implicit none
      integer, intent(in) :: nk
      real(8), intent(in) :: kp(3,nk)
      complex*16, intent(out) :: skb(norb,norb,*), hkb(norb,norb,*), sdb(norb,norb,3,*), hdb(norb,norb,3,*), akb(norb,norb,3,*)
      complex*16, allocatable :: W(:,:), F(:,:), OH(:,:), OS(:,:), OR(:,:)
      complex*16, parameter :: zone = (1.0d0,0.0d0), zzero = (0.0d0,0.0d0), ci = (0.0d0,1.0d0)
      integer :: iRp, ik, t, a, b, nj, m

      allocate(W(nR,4*nk), F(nR,nk), OH(4*nk,ntri), OS(4*nk,ntri), OR(nk,3*ntri))
      do ik = 1, nk
         m = 4*(ik-1)
         do iRp = 1, nR
            F(iRp,ik) = exp(cmplx(0.0d0, kp(1,ik)*Rx_global(iRp) + kp(2,ik)*Ry_global(iRp) + kp(3,ik)*Rz_global(iRp), 8))
            W(iRp,m+1) = F(iRp,ik)
            W(iRp,m+2) = cmplx(0.0d0, Rx_global(iRp), 8)*F(iRp,ik)
            W(iRp,m+3) = cmplx(0.0d0, Ry_global(iRp), 8)*F(iRp,ik)
            W(iRp,m+4) = cmplx(0.0d0, Rz_global(iRp), 8)*F(iRp,ik)
         end do
      end do
      call zgemm('T','N', 4*nk, ntri,   nR, zone, W, nR, hl_pack, nR, zzero, OH, 4*nk)
      call zgemm('T','N', 4*nk, ntri,   nR, zone, W, nR, sl_pack, nR, zzero, OS, 4*nk)
      call zgemm('T','N', nk,   3*ntri, nR, zone, F, nR, rl_pack, nR, zzero, OR, nk)
      do ik = 1, nk
         m = 4*(ik-1)
         do t = 1, ntri
            a = tri_a(t); b = tri_b(t)
            hkb(a,b,ik) = OH(m+1,t); skb(a,b,ik) = OS(m+1,t)
            do nj = 1, 3
               hdb(a,b,nj,ik) = OH(m+1+nj,t); sdb(a,b,nj,ik) = OS(m+1+nj,t)
               akb(a,b,nj,ik) = OR(ik,(nj-1)*ntri+t) + ci*OS(m+1+nj,t)
            end do
            if (b < a) then
               hkb(b,a,ik) = conjg(hkb(a,b,ik)); skb(b,a,ik) = conjg(skb(a,b,ik))
               do nj = 1, 3
                  hdb(b,a,nj,ik) = conjg(hdb(a,b,nj,ik)); sdb(b,a,nj,ik) = conjg(sdb(a,b,nj,ik))
                  akb(b,a,nj,ik) = conjg(akb(a,b,nj,ik)) + ci*conjg(sdb(a,b,nj,ik))
               end do
            end if
         end do
      end do
      deallocate(W, F, OH, OS, OR)
   end subroutine get_vme_kernels_batch
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
      hk_ev,e,vme,M1,T1,Wout,skip_a4,scramble_blocks)
      implicit none

      complex*16, intent(out), optional :: Wout(norb,norb)
      ! skip_a4 = .true.: plain eigensystem (no Eq. (A4) rotation), used by the shift_covariant path,
      ! which also gets the block scramble test hook. Not combined with Wout.
      logical, intent(in), optional :: skip_a4
      ! scramble_blocks = .true.: apply the block scramble test hook to THIS diagonalisation (only the
      ! central point of get_rgen_blocks: at the +-dk neighbours the pairs are genuinely split and their
      ! eigenvectors fixed, so rotating them would inject an error no eigensolver makes)
      logical, intent(in), optional :: scramble_blocks
      complex*16 :: Wsave(norb,norb)
      integer :: norb, nj, nn, nnp

      dimension skernel(norb,norb)
      dimension s_work(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(norb,norb,3)      ! <-- reordered
      dimension hderkernel(norb,norb,3)      ! <-- reordered
      dimension akernel(norb,norb,3)         ! <-- reordered
      dimension e(norb)
      dimension hk_ev(norb,norb)
      dimension vme(norb,norb,3)             ! <-- reordered

      real*8 e
      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel,s_work
      complex*16 hk_ev,vme
      complex*16 M1(norb,norb), T1(norb,norb)
      complex*16, parameter :: cone=(1.0d0,0.0d0), czero=(0.0d0,0.0d0), ci=(0.0d0,1.0d0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      e=0.0d0
      if (iflag_orthonormal) then
         call diagoz(norb,e,hkernel)
         hk_ev(:,:)=hkernel(:,:)
      else if (.not. iflag_orthonormal) then
         s_work(:,:) = skernel(:,:)
         call diagoz_gen(norb,e,hkernel,s_work)
         hk_ev(:,:)=hkernel(:,:)
      end if
      call phase_eigvec_nk(norb,hk_ev)

      call build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)

      if (present(skip_a4)) then
         if (skip_a4) then
            if (present(scramble_blocks)) then
               if (scramble_blocks) call scramble_degenerate_groups()
            end if
            if (gauge_test_on) call scramble_gauge()
            return
         end if
      end if

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
         if (gauge_test_on) call scramble_gauge()
         ! Change of basis old -> new. With orthonormal Wannier functions it is Wsave^H hk_ev; in a
         ! non-orthonormal basis the eigenvectors are S-orthonormal, so it is Wsave^H S hk_ev.
         if (iflag_orthonormal) then
            call zgemm('C','N',norb,norb,norb, cone, Wsave, norb, hk_ev, norb, czero, Wout, norb)
         else
            call zgemm('N','N',norb,norb,norb, cone, skernel, norb, hk_ev, norb, czero, M1, norb)
            call zgemm('C','N',norb,norb,norb, cone, Wsave, norb, M1, norb, czero, Wout, norb)
         end if
      else
         if (apply_a4_rotation(norb,e,hk_ev,vme)) then
            call phase_eigvec_nk(norb,hk_ev)
            call build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)
         end if
         if (gauge_test_on) call scramble_gauge()
      end if

   contains

      ! Gauge test hook (see gauge_test_on): theta_n(k) = 1.5 sin(20 khat_n . k + n), khat_n a
      ! band-dependent in-plane direction plus a small z part, so neighbouring bands get unrelated phases.
      subroutine scramble_gauge()
         integer :: n
         real(8) :: th
         do n = 1, norb
            th = gt_amp*sin(20.0d0*(gt_k(1)*cos(0.7d0*n) + gt_k(2)*sin(0.7d0*n) + 0.3d0*gt_k(3)) + dble(n))
            hk_ev(:,n) = hk_ev(:,n)*cmplx(cos(th), sin(th), 8)
         end do
         call build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)
      end subroutine scramble_gauge

      ! Block scramble hook (scramble_blocks): inside every group of numerically degenerate bands apply a
      ! k-dependent unitary built from Givens rotations (deterministic in gt_k, so thread-safe).
      subroutine scramble_degenerate_groups()
         integer :: lo, hi, i, j
         real(8) :: th, ph
         complex*16 :: c1(norb), c2(norb)
         lo = 1
         do while (lo <= norb)
            hi = lo
            do while (hi < norb)
               if (e(hi+1) - e(hi) >= cov_exact_tol) exit
               hi = hi + 1
            end do
            do i = lo, hi - 1
               do j = i + 1, hi
                  th = 1.3d0*sin(17.0d0*(gt_k(1) + 1.7d0*gt_k(2)) + dble(3*i + j))
                  ph = 2.1d0*cos(13.0d0*(gt_k(2) - 0.9d0*gt_k(1)) + dble(i + 5*j))
                  c1 = hk_ev(:,i); c2 = hk_ev(:,j)
                  hk_ev(:,i) = cos(th)*c1 - sin(th)*cmplx(cos(ph), sin(ph), 8)*c2
                  hk_ev(:,j) = sin(th)*cmplx(cos(ph), -sin(ph), 8)*c1 + cos(th)*c2
               end do
            end do
            lo = hi + 1
         end do
         call build_vme_blocks(norb,hderkernel,akernel,hk_ev,e,M1,T1,vme)
      end subroutine scramble_degenerate_groups

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
   ! start from the SAME basis (with the rotation off, opticx reproduces the trusted
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
   !> With OME_sp = none, read the exciton-window basis data that get_ome_sp would otherwise have built:
   !! the Eq. (A4) rotation a4_W (both orders) and, for norder = 2, the window states win_evec, win_sevec,
   !! win_xi2 used by the covariant X_nm. They are stored in the .omesp the excitonic stage reads (the linear
   !! one for norder = 1, the nonlinear one for norder = 2). Files written before 2026-10-06 do not contain
   !! them; then nothing is loaded, a4_W_ready / win_ready stay false, and get_ome_ex stops with a message.
   subroutine load_omesp_basis(norder)
      implicit none
      integer, intent(in) :: norder
      integer :: u, ios, iflag_r, npts_r, nb_r, tag, flags, norb_r, ibz, i
      integer(8) :: pos, nk, nb
      character(len=16) :: line
      logical :: ex
      character(len=:), allocatable :: fname

      if (norder == 1) then
         fname = 'ome_linear_sp_'//trim(material_name)//'.omesp'
         inquire(file=fname, exist=ex)
         if (.not. ex) return
         open(newunit=u, file=fname, status='old', action='read')
         read(u,*,iostat=ios) iflag_r
         do ibz = 1, npointstotal*(1 + nband_ex*nband_ex)
            if (ios /= 0) exit
            read(u,*,iostat=ios)
         end do
         line = ''
         if (ios == 0) read(u,'(a)',iostat=ios) line
         if (ios == 0 .and. trim(line) == '#A4W') then
            if (allocated(a4_W)) deallocate(a4_W)
            allocate(a4_W(npointstotal,nband_ex,nband_ex))
            read(u,*,iostat=ios) a4_W
            a4_W_ready = (ios == 0)
         end if
         close(u)
      else
         fname = 'ome_nonlinear_sp_'//trim(material_name)//'.omesp'
         inquire(file=fname, exist=ex)
         if (.not. ex) return
         open(newunit=u, file=fname, form='unformatted', access='stream', status='old', action='read')
         read(u,iostat=ios) iflag_r
         if (ios == 0) read(u,iostat=ios) npts_r, nb_r
         if (ios /= 0 .or. npts_r /= npointstotal .or. nb_r /= nband_ex) then
            close(u); return
         end if
         ! skip the fixed part: header, k-list, per-k records (ek, v, Berry, shift vector, sum-rule
         ! derivative), d|v|/dk and the parallel-transported dv/dk
         nk = npointstotal; nb = nband_ex
         pos = 1 + 4 + 8 + 3*nk*8 + nk*(nb*8 + 2*3*nb*nb*16 + 9*nb*nb*8 + 9*nb*nb*16) &
               + nk*9*nb*nb*8 + nk*9*nb*nb*16
         read(u, pos=pos, iostat=ios) tag, flags, norb_r
         if (ios == 0 .and. tag == OMESP_TAG) then
            if (btest(flags, 0)) then
               if (allocated(a4_W)) deallocate(a4_W)
               allocate(a4_W(npointstotal,nband_ex,nband_ex))
               read(u, iostat=ios) a4_W
               a4_W_ready = (ios == 0)
            end if
            if (btest(flags, 1) .and. norb_r == norb .and. ios == 0) then
               if (allocated(win_evec))  deallocate(win_evec)
               if (allocated(win_sevec)) deallocate(win_sevec)
               if (allocated(win_xi2))   deallocate(win_xi2)
               allocate(win_evec(norb,nband_ex,npointstotal), win_sevec(norb,nband_ex,npointstotal))
               allocate(win_xi2(npointstotal,3,nband_ex,nband_ex))
               read(u, iostat=ios) win_evec
               if (ios == 0) read(u, iostat=ios) win_sevec
               if (ios == 0) read(u, iostat=ios) win_xi2
               win_ready = (ios == 0)
            end if
         end if
         close(u)
      end if
      if (a4_W_ready) write(*,*) '   Eq. (A4) basis of the exciton window read from '//trim(fname)
      if (win_ready)  write(*,*) '   Exciton-window states (covariant X_nm) read from '//trim(fname)
   end subroutine load_omesp_basis

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
      ! exactly the O(1) error the flag exists to prevent.
      if (.not. fk_ex_loaded) return

      ! Rotation switched off: the single-particle states were never rotated, so fk_ex as read from
      ! Xatu already matches them and is fit to use.
      if (.not. a4_enabled) then
         fk_ex_basis_ok = .true.
         return
      end if
      ! No rotation matrices: fk_ex CANNOT be carried over, so it stays marked unusable. This is only
      ! a warning and not a stop, because a second-order OME cache hit never reads
      ! fk_ex at all -- get_ome_ex issues the hard error at the point where it commits to the k-loop.
      if (.not. a4_W_ready) then
         write(*,*) '   WARNING (ome_sp): the Eq. (A4) rotation is active but the per-k rotation'
         write(*,*) '            matrices are unavailable (OME_sp = none, and the .omesp file was'
         write(*,*) '            written by a version that did not store them). fk_ex CANNOT be carried'
         write(*,*) '            into the rotated basis. Regenerate the .omesp once with OME_sp set.'
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