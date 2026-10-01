module sigma_second_sp
  use constants_math
  use parser_input_file, &
  only:nf,e1,e2,eta,nw,response_text,broadening_type_text, &
  freq_ratio,e1b,e2b,nwb,two_freq_grid
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
    complex*16, allocatable :: vme_der_pt_ex_band(:,:,:,:,:)
    logical :: vme_der_pt_found

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
    allocate(vme_der_pt_ex_band(npointstotal,3,3,nband_ex,nband_ex))
    call read_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,berry_eigen_ex_band, &
                                   gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek, &
                                   vme_abs_der_ex_band,vme_abs_der_found, &
                                   vme_der_pt_ex_band,vme_der_pt_found)
    write(*,*) '    Optical matrix elements (sp) have been read from file'
    if (response_text == 'shift_shiftvector' .and. .not. vme_abs_der_found) then
      write(*,*) 'ERROR (sigma_second_sp): the .omesp file has no derivative of |v| (it was written by an'
      write(*,*) '       older version). Regenerate it with OME_sp = nonlinear; shift_shiftvector needs it.'
      stop 1
    end if
    if ((response_text == 'shg' .or. response_text == 'electrooptic' .or. &
         response_text == 'rectification' .or. response_text == 'general') &
        .and. .not. vme_der_pt_found) then
      write(*,*) 'ERROR (sigma_second_sp): the .omesp file has no gauge-fixed generalized derivative of v'
      write(*,*) '       (it was written by an older version). Regenerate it with OME_sp = nonlinear; shg needs it.'
      stop 1
    end if

    !compute shift conductivity
    if (nwp.eq.1 .and. nwq.eq.(-1)) then
      call get_sigma_shift_sp(npointstotal,nband_ex,berry_eigen_ex_band, &
                        gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek,vme_abs_der_ex_band)
    end if

    !compute shg susceptibility (the w_q = w_p branch; kept as its own entry point so the
    ! long-validated SHG output file and path are untouched)
    if (nwp.eq.1 .and. nwq.eq.1 .and. response_text == 'shg') then
      call get_sigma_shg_sp(npointstotal,nband_ex,vme_ex_band,ek,vme_der_pt_ex_band)
    end if

    ! general second-order sigma^{abc}(w_p+w_q; w_p, w_q), Eq. (A3a). nwq = 0 is the electro-optic
    ! (Pockels) branch sigma(w; w, 0); nwq = -1 optical rectification sigma(0; w, -w); 'general' uses
    ! Frequency_ratio or the Energy_variables_2 grid. HANDOFF 8.36.
    if (response_text == 'electrooptic') then
      freq_ratio = 0.0d0; two_freq_grid = .false.
      call get_sigma_general_sp(npointstotal,nband_ex,vme_ex_band,ek,vme_der_pt_ex_band,'electrooptic')
    else if (response_text == 'rectification') then
      ! This branch DOES reproduce the shift current (HANDOFF 8.36c). Verified on hBN over 1-16 eV:
      ! both this and shift_shiftvector peak at 8.57 eV, both have the same eta dependence (peak ratio
      ! 1.378 vs 1.373 when eta halves), and across the whole resonant region
      !     sigma_A3a / sigma_shift = -0.2479 +- 0.011  (eta = 0.05; scatter shrinks as eta -> 0),
      ! i.e. a constant factor of -1/4. The sign and the 1/4 are convention (the shift code carries the
      ! paper's -pi/2 and the (I^abc + I^acb) form; a real field also brings the usual factor of 2), NOT
      ! a physics difference. An earlier note here claimed this branch was an "O(eta) regulator residue
      ! with no resonant structure" -- that was WRONG, drawn from a 1.5-6.0 eV scan that never reached
      ! hBN's 7.25 eV gap. Retracted; see HANDOFF 8.36c.
      write(*,*) '    NOTE: the rectification branch reproduces the shift current up to a constant'
      write(*,*) '          factor of -1/4 (convention; HANDOFF 8.36c). For the production shift'
      write(*,*) '          current use Response = shift_shiftvector.'
      freq_ratio = -1.0d0; two_freq_grid = .false.
      call get_sigma_general_sp(npointstotal,nband_ex,vme_ex_band,ek,vme_der_pt_ex_band,'rectification')
    else if (response_text == 'general') then
      call get_sigma_general_sp(npointstotal,nband_ex,vme_ex_band,ek,vme_der_pt_ex_band,'general')
    end if


    deallocate(vme_abs_der_ex_band)
    deallocate(vme_der_pt_ex_band)

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
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Single-particle SHG sigma^{lam,alpha,beta}(2*omega; omega, omega). Method A, length gauge, direct
! generalized derivative (Taghizadeh, Hipolito & Pedersen, PRB 96, 195413 (2017), Eq. A3a), specialised to
! a cold insulator (occupations f_n fixed by band index, k-independent: the two Eq. A3a terms containing
! df/dk then vanish identically, leaving only the two coded below):
!
!   term 1 (l != n and l != m -- these are the only exclusions Eq. A3a states; n = m is NOT excluded, see
!           the comment in the l-loop below):
!     p^lam_nm (g^alpha_ln p^beta_ml - g^alpha_ml p^beta_ln) / (E_ml E_ln (2*hw - E_mn))
!   term 2 (n != m):
!     -p^lam_nm/(2*hw - E_mn) * (g^alpha_mn/E_mn)_{;k^beta}
!
! g^alpha_MN(w) = (f_N-f_M) p^alpha_MN/(hw-E_MN), E_MN = E_M-E_N, hw = hbar*omega + i*eta (every frequency
! carries its own +i*eta; omega_2 = omega_p+omega_q -> 2*hw, same SHG convention as the excitonic code).
! (g^alpha_mn/E_mn);k^beta is expanded by the product rule (f_nm, E_mn are k-independent/real for a T=0
! insulator, so their generalized derivative reduces to an ordinary one, Eq. 6a of the paper):
!   (g^alpha_mn/E_mn);k^beta = f_nm [ dphi/dE(E_mn) (dE_mn/dk^beta) p^alpha_mn + phi(E_mn) (p^alpha_mn);k^beta ]
!   phi(E) = 1/[(hw-E) E],  dphi/dE = (1/hw)[1/(hw-E)^2 - 1/E^2],  dE_mn/dk^beta = p^beta_mm - p^beta_nn
! NOTE the index order p^alpha_MN (not p^alpha_nm): g^alpha_mn = f_nm p^alpha_mn/(hw - E_mn), paper Eq. (10).
! This matters, and is not cosmetic: term 2 is p^lambda_nm times g^alpha_mn, so the two momentum factors carry
! OPPOSITE band-index order and their gauge phases cancel, e^{i(phi_m-phi_n)} e^{i(phi_n-phi_m)} = 1. Written
! with p^alpha_nm in both places (as this routine did until 2026-09-24) term 2 scales as e^{2i(phi_m-phi_n)},
! i.e. the SHG conductivity becomes GAUGE DEPENDENT and therefore unphysical. See the corrected note below.
! (p^alpha_mn);k^beta is the gauge-fixed (parallel-transported) generalized derivative computed in
! get_berry_eigen_fourpoint (ome_sp.f90) and read here as vme_der_pt_ex_band -- NOT gen_der_ex_band, the
! sum-rule form used by shift_sumrule: that form gives an unreliable, near-zero result on small band
! windows (verified independently this session: an independent NumPy implementation of a SECOND,
! differently-derived SHG formula, Rashkeev/Lambrecht/Segall based on Aversa & Sipe, shows the SAME
! near-zero pathology when its generalized derivative is evaluated via the sum rule instead of a direct
! k-derivative -- see HANDOFF.md). C_ee = C_ie = 1/4 in these units (e=hbar=m=1; spin factor g dropped, as
! elsewhere in opticx), combined with the code's usual 1/(Nk*V) discretisation of the BZ sum (paper Eq. A4
! and the Sigma_k -> A/(2pi)^D convention below it). Validated against an independent NumPy evaluation of
! this same equation (tools/ipa_shg_numpy.py) and cross-checked against a second, independently-derived
! formula (Rashkeev, Lambrecht & Segall, arXiv:cond-mat/9709185). CORRECTED 2026-09-24: the large residual C3
! error previously seen on hBN's 2-orbital model was NOT a 2-band basis-truncation effect -- that explanation
! is retracted. It was the p^alpha_nm/p^alpha_mn index slip described above, which both this routine and
! tools/ipa_shg_numpy.py shared (so the two "independent" implementations agreed on the same error, and the
! synthetic-data cross-check could not see it: it is a gauge-covariance defect, invisible to any test that
! does not vary the gauge). With the index corrected, hBN's 2-band SHG is exactly gauge invariant and
! satisfies D3h to 6e-6 -- hBN IS a valid validation target for this formula. Near-degenerate bands are dropped
! (eps_deg) exactly like shift_shiftvector; vme_der_pt itself is separately clipped at the source
! (clip_threshold in ome_sp.f90). No zgemm/frequency-chunking restructuring yet (see get_shg_intens_ex for
! that pattern); fine for the band windows used so far, a candidate for later optimisation on large windows.
  subroutine get_sigma_shg_sp(npointstotal,nband_ex,vme_ex_band,ek,vme_der_pt_ex_band)
    implicit none
    !in/out
    integer,    intent(in) :: nband_ex,npointstotal
    complex*16, intent(in) :: vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    real*8,     intent(in) :: ek(npointstotal,nband_ex)
    complex*16, intent(in) :: vme_der_pt_ex_band(npointstotal,3,3,nband_ex,nband_ex)

    real*8,     allocatable :: wp(:)
    complex*16, allocatable :: sigma_w_sp(:,:,:,:)
    complex*16, allocatable :: sigma_raw(:,:,:,:)
    real*8 :: eta2
    integer :: nj,njp,njpp,iw
    complex*16, allocatable :: hwp(:), hwsum(:)

    real*8,     allocatable :: e_nband(:)
    complex*16, allocatable :: vme_nband(:,:,:)
    complex*16, allocatable :: vme_der_pt_nband(:,:,:,:)
    complex*16, allocatable :: sigma_w_sp_t(:,:,:,:)
    integer :: ibz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(wp(nw))
    allocate(sigma_w_sp(3,3,3,nw))
    allocate(hwp(nw), hwsum(nw))
    call initialize_sigma_second_arrays(nw, wp, eta2, sigma_w_sp)

    ! SHG is the omega_q = omega_p branch of the general kernel (HANDOFF 8.36). Both frequencies carry
    ! +i*eta, so the outer pole is 2*omega + 2i*eta -- the paper's convention, and the one CLAUDE.md
    ! records as required for method A == method B ("build omega_2 as the SUM of the two complex
    ! frequencies, never as an independent single eta").
    do iw = 1, nw
      hwp(iw)   = cmplx(wp(iw), eta2, 8)
      hwsum(iw) = hwp(iw) + hwp(iw)
    end do

    write(*,*) '    Evaluating SHG susceptibility (sp)...'

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP   SHARED(npointstotal, nband_ex, ek, vme_ex_band, vme_der_pt_ex_band, nw, hwp, hwsum, sigma_w_sp) &
    !$OMP   PRIVATE(ibz, e_nband, vme_nband, vme_der_pt_nband, sigma_w_sp_t, nj, njp)

    allocate(e_nband(nband_ex))
    allocate(vme_nband(3,nband_ex,nband_ex))
    allocate(vme_der_pt_nband(3,3,nband_ex,nband_ex))
    allocate(sigma_w_sp_t(3,3,3,nw))
    sigma_w_sp_t = (0.0d0,0.0d0)

    !$OMP DO SCHEDULE(DYNAMIC)
    do ibz = 1, npointstotal
      e_nband(:) = ek(ibz,:)
      do nj = 1,3
        vme_nband(nj,:,:) = vme_ex_band(ibz,nj,:,:)
        do njp = 1,3
          vme_der_pt_nband(nj,njp,:,:) = vme_der_pt_ex_band(ibz,nj,njp,:,:)
        end do
      end do
      call get_shg_intens_sp(nband_ex, nw, e_nband, vme_nband, vme_der_pt_nband, hwp, hwsum, sigma_w_sp_t)
    end do
    !$OMP END DO

    !$OMP CRITICAL
      sigma_w_sp = sigma_w_sp + sigma_w_sp_t
    !$OMP END CRITICAL

    deallocate(e_nband, vme_nband, vme_der_pt_nband, sigma_w_sp_t)
    !$OMP END PARALLEL

    ! symmetrise over the two (identical-frequency) field indices -- Eq. A3a is not manifestly symmetric
    ! under alpha<->beta (text below Eq. A4); same treatment as the excitonic SHG driver (get_sigma_shg_ex)
    allocate(sigma_raw(3,3,3,nw))
    sigma_raw = sigma_w_sp
    do njpp = 1,3
      do njp = 1,3
        do nj = 1,3
          sigma_w_sp(nj,njp,njpp,:) = 0.5d0*(sigma_raw(nj,njp,njpp,:) + sigma_raw(nj,njpp,njp,:))
        end do
      end do
    end do
    deallocate(sigma_raw)

    call print_shg_second_sp(nw, wp, sigma_w_sp)
    write(*,*) '    SHG susceptibility (sp) has been printed'

    deallocate(wp, sigma_w_sp)
  end subroutine get_sigma_shg_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENERAL second-order single-particle conductivity sigma^{abc}(w_p + w_q; w_p, w_q), Taghizadeh 2017
! Eq. (A3a) for a cold intrinsic semiconductor. Added 2026-09-24 (HANDOFF 8.36).
!
! Frequency pairs POSITIONALLY with the Cartesian index: alpha <-> w_p, beta <-> w_q (Eq. 9a at first
! order; term 3's d f_nm/(hbar w_p d k^alpha); the third-order g^gamma <-> w_s). Eq. (A3a) is therefore
! NOT symmetric under alpha<->beta alone when w_p /= w_q; the physical tensor is the intrinsic-permutation
! average over the PAIRS (alpha,w_p) <-> (beta,w_q) (paper, text below Eq. A4):
!
!     sigma_sym^{l,a,b}(wp,wq) = 1/2 [ sigma^{l,a,b}(wp,wq) + sigma^{l,b,a}(wq,wp) ]
!
! Implemented as two passes of the same kernel over the BZ sharing one hwsum: pass 1 with hwp = hbar*wp,
! pass 2 with hwp = hbar*wq, then the index transpose on pass 2. For w_p = w_q (SHG) the two passes are
! identical and the average degenerates to the plain b<->c swap the old SHG driver did. Every frequency
! carries +i*eta, hence hwsum = hbar(wp+wq) + 2i*eta (paper convention; HANDOFF 8.36 records why the DC
! shift code's "eta only on w_p" rule does NOT transfer to Eq. A3a).
  subroutine get_sigma_general_sp(npointstotal, nband_ex, vme_ex_band, ek, vme_der_pt_ex_band, tag)
    implicit none
    integer,    intent(in) :: nband_ex, npointstotal
    complex*16, intent(in) :: vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    real*8,     intent(in) :: ek(npointstotal,nband_ex)
    complex*16, intent(in) :: vme_der_pt_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    character(len=*), intent(in) :: tag

    complex*16, allocatable :: hwp(:), hwq(:), hwsum(:)
    real*8,     allocatable :: wpg(:), wqg(:)
    complex*16, allocatable :: sigA(:,:,:,:), sigB(:,:,:,:), sigS(:,:,:,:)
    real*8  :: eta2, wrange, wrangeb, wa, wb
    integer :: nfreq, iw, iwb, idx, nj, njp, njpp
    logical :: same_freq

    eta2 = eta/27.211385d0
    if (two_freq_grid) then
      nfreq = nw*nwb
    else
      nfreq = nw
    end if
    allocate(hwp(nfreq), hwq(nfreq), hwsum(nfreq), wpg(nfreq), wqg(nfreq))

    wrange = e2 - e1
    if (two_freq_grid) then
      wrangeb = e2b - e1b
      idx = 0
      do iw = 1, nw
        wa = (e1 + wrange/dble(nw)*dble(iw-1))/27.211385d0
        do iwb = 1, nwb
          wb = (e1b + wrangeb/dble(nwb)*dble(iwb-1))/27.211385d0
          idx = idx + 1
          wpg(idx) = wa
          wqg(idx) = wb
        end do
      end do
    else
      do iw = 1, nw
        wpg(iw) = (e1 + wrange/dble(nw)*dble(iw-1))/27.211385d0
        wqg(iw) = freq_ratio*wpg(iw)
      end do
    end if

    ! Every frequency carries +i*eta, so hwsum = hbar(w_p + w_q) + 2i*eta. The 2i*eta is NOT optional
    ! and is not "killed by the signs" at w_q = -w_p: term 1 of Eq. (A3a) ALLOWS n = m (CLAUDE.md is
    ! emphatic -- dropping it broke C3 by 80%), and there E_mn = 0, so the outer denominator is
    ! hbar*w_sum - E_mn = hbar*w_sum. Setting w_sum = 0 exactly makes that 0/0: verified, it produces
    ! NaN at 300 of 300 frequencies (HANDOFF 8.36). The +2i*eta is what regularises it.
    do idx = 1, nfreq
      hwp(idx)   = cmplx(wpg(idx), eta2, 8)
      hwq(idx)   = cmplx(wqg(idx), eta2, 8)
      hwsum(idx) = hwp(idx) + hwq(idx)
      if (abs(hwsum(idx)) < 1.0d-30) then
        write(*,*) 'ERROR (get_sigma_general_sp): hbar(w_p+w_q) is exactly zero at a grid point.'
        write(*,*) '       Term 1 of Eq. (A3a) includes n = m, where E_mn = 0, so the outer'
        write(*,*) '       denominator would be 0/0. Use a nonzero eta.'
        stop 1
      end if
    end do

    same_freq = (.not. two_freq_grid) .and. (abs(freq_ratio - 1.0d0) < 1.0d-12)

    allocate(sigA(3,3,3,nfreq), sigB(3,3,3,nfreq), sigS(3,3,3,nfreq))

    write(*,*) '    Evaluating second-order conductivity (sp): ', trim(tag)
    write(*,'(A,I0,A)') '        ', nfreq, ' frequency pairs'

    call run_second_kernel_sp(npointstotal, nband_ex, vme_ex_band, ek, vme_der_pt_ex_band, &
                              nfreq, hwp, hwsum, sigA)
    if (same_freq) then
      sigB = sigA
    else
      call run_second_kernel_sp(npointstotal, nband_ex, vme_ex_band, ek, vme_der_pt_ex_band, &
                                nfreq, hwq, hwsum, sigB)
    end if

    do nj = 1,3
      do njp = 1,3
        do njpp = 1,3
          sigS(nj,njp,njpp,:) = 0.5d0*(sigA(nj,njp,njpp,:) + sigB(nj,njpp,njp,:))
        end do
      end do
    end do

    call print_second_general_sp(nfreq, wpg, wqg, sigS, tag)
    deallocate(hwp, hwq, hwsum, wpg, wqg, sigA, sigB, sigS)
  end subroutine get_sigma_general_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! One pass of the Eq. (A3a) kernel over the BZ for a given (hwp, hwsum) pair list.
  subroutine run_second_kernel_sp(npointstotal, nband_ex, vme_ex_band, ek, vme_der_pt_ex_band, &
                                  nfreq, hwp, hwsum, sigma)
    implicit none
    integer,    intent(in)  :: nband_ex, npointstotal, nfreq
    complex*16, intent(in)  :: vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    real*8,     intent(in)  :: ek(npointstotal,nband_ex)
    complex*16, intent(in)  :: vme_der_pt_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    complex*16, intent(in)  :: hwp(nfreq), hwsum(nfreq)
    complex*16, intent(out) :: sigma(3,3,3,nfreq)

    real*8,     allocatable :: e_nband(:)
    complex*16, allocatable :: vme_nband(:,:,:), vme_der_pt_nband(:,:,:,:), sigma_t(:,:,:,:)
    integer :: ibz, nj, njp

    sigma = (0.0d0,0.0d0)

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP   SHARED(npointstotal, nband_ex, ek, vme_ex_band, vme_der_pt_ex_band, nfreq, hwp, hwsum, sigma) &
    !$OMP   PRIVATE(ibz, e_nband, vme_nband, vme_der_pt_nband, sigma_t, nj, njp)
    allocate(e_nband(nband_ex))
    allocate(vme_nband(3,nband_ex,nband_ex))
    allocate(vme_der_pt_nband(3,3,nband_ex,nband_ex))
    allocate(sigma_t(3,3,3,nfreq))
    sigma_t = (0.0d0,0.0d0)
    !$OMP DO SCHEDULE(DYNAMIC)
    do ibz = 1, npointstotal
      e_nband(:) = ek(ibz,:)
      do nj = 1,3
        vme_nband(nj,:,:) = vme_ex_band(ibz,nj,:,:)
        do njp = 1,3
          vme_der_pt_nband(nj,njp,:,:) = vme_der_pt_ex_band(ibz,nj,njp,:,:)
        end do
      end do
      call get_shg_intens_sp(nband_ex, nfreq, e_nband, vme_nband, vme_der_pt_nband, hwp, hwsum, sigma_t)
    end do
    !$OMP END DO
    !$OMP CRITICAL
      sigma = sigma + sigma_t
    !$OMP END CRITICAL
    deallocate(e_nband, vme_nband, vme_der_pt_nband, sigma_t)
    !$OMP END PARALLEL
  end subroutine run_second_kernel_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Columns: hbar*w_p (eV), hbar*w_q (eV), then a=x,y,z; b=x,y,z; c=x,y,z (c fastest): Re, Im.
! Units uA nm/V^2. Both frequencies are printed so the file is self-describing in ratio and 2D-grid mode.
  subroutine print_second_general_sp(nfreq, wpg, wqg, sigma, tag)
    implicit none
    integer :: iounit93
    integer,    intent(in) :: nfreq
    real*8,     intent(in) :: wpg(nfreq), wqg(nfreq)
    complex*16, intent(in) :: sigma(3,3,3,nfreq)
    character(len=*), intent(in) :: tag
    integer :: iw, ia, ib, ic
    real*8  :: feps
    feps = sigma2_au_to_si
    open(newunit=iounit93, file='second_'//trim(tag)//'_lengthgauge_'//trim(material_name)//'.dat')
    write(iounit93,'(A)') '# hbar*w_p(eV) hbar*w_q(eV) | sigma^{abc}(w_p+w_q;w_p,w_q) (Re,Im) uA nm/V^2, abc=xxx,xxy,...,zzz'
    do iw = 1, nfreq
      write(iounit93,'(2ES18.10,54ES18.10)') wpg(iw)*27.211385d0, wqg(iw)*27.211385d0, &
        ( ( ( real(feps*sigma(ia,ib,ic,iw)), aimag(feps*sigma(ia,ib,ic,iw)), ic=1,3 ), ib=1,3 ), ia=1,3 )
    end do
    close(iounit93)
  end subroutine print_second_general_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! GENERALISED 2026-09-24 to arbitrary (omega_p, omega_q); see HANDOFF 8.36.
  ! Taghizadeh 2017 Eq. (A3a) is written for arbitrary frequencies. Specialised to a cold intrinsic
  ! semiconductor only two of its four terms survive, and in BOTH of them omega_q appears ONLY through the
  ! sum omega_p + omega_q:
  !    term 1: p^lam_nm (g^a_ln p^b_ml - g^a_ml p^b_ln) / (E_ml E_ln [hbar(wp+wq) - E_mn])
  !    term 2: -p^lam_nm / (hbar(wp+wq) - E_mn) * (g^a_mn/E_mn);k^b
  ! because p^b carries no frequency and the generalized derivative ;k^b carries none either -- the two
  ! terms that DO use hbar*omega_p on its own are exactly the d f_n/dk ones that vanish at T=0.
  ! The frequency pairs POSITIONALLY with the Cartesian index (alpha<->omega_p, beta<->omega_q): confirmed
  ! by Eq. (9a) at first order, by term 3's d f_nm/(hbar*omega_p d k^alpha), and by the third-order
  ! g^gamma_mn = f_nm p^gamma_mn/(hbar*omega_s - E_mn) with (wp,wq,ws)<->(alpha,beta,gamma).
  ! So this kernel needs exactly TWO complex frequencies per grid point:
  !    hwp(iw)   = hbar*omega_p + i*eta          -> enters g^alpha only
  !    hwsum(iw) = hbar*(omega_p+omega_q) + 2i*eta -> the outer pole only
  ! Every second-order process is one branch (HANDOFF 8.36):
  !    SHG (2w;w,w)             hwp = w+ie,  hwsum = 2w+2ie
  !    electro-optic (w;w,0)    hwp = w+ie,  hwsum =  w+2ie
  !    rectification (0;w,-w)   hwp = w+ie,  hwsum =  0+2ie
  !    general (w1+w2;w1,w2)    hwp = w1+ie, hwsum = w1+w2+2ie
  !> Single-particle second-order kernel at one k-point: Eq. (A3a) of Taghizadeh, Hipolito
  !! & Pedersen, PRB 96, 195413 (2017), specialised to a cold insulator so the two terms
  !! carrying df/dk drop out. Two survive: term 1 (l /= n, l /= m, but n = m IS allowed --
  !! it is the intraband/group-velocity piece) and term 2, the generalized-derivative piece.
  !! Term 2 pairs p^lambda_nm with g^alpha_MN: the two momentum factors must carry OPPOSITE
  !! band-index order or their gauge phases fail to cancel and sigma becomes gauge dependent.
  !! @param e_nband, vme_nband     Energies and velocity elements at this k-point.
  !! @param vme_der_pt_nband       Gauge-fixed generalized derivative (p^a_nm);k^b.
  !! @param hwp                    Complex hbar*omega_p.
  !! @param hwsum                  Complex hbar*(omega_p + omega_q); never zero (0/0 in term 1).
  !! @param sigma_w_sp             Accumulated onto.
  !! @return void
  subroutine get_shg_intens_sp(nband_ex, nw, e_nband, vme_nband, vme_der_pt_nband, hwp, hwsum, sigma_w_sp)
    implicit none
    integer,    intent(in)    :: nband_ex, nw
    real*8,     intent(in)    :: e_nband(nband_ex)
    complex*16, intent(in)    :: vme_nband(3,nband_ex,nband_ex)
    complex*16, intent(in)    :: vme_der_pt_nband(3,3,nband_ex,nband_ex)
    complex*16, intent(in)    :: hwp(nw)      ! hbar*omega_p   + i*eta
    complex*16, intent(in)    :: hwsum(nw)    ! hbar*(wp + wq) + 2i*eta
    complex*16, intent(inout) :: sigma_w_sp(3,3,3,nw)

    integer    :: iw, nj, njp, njpp, nn, nnp, nl
    real*8     :: fnn, fnnp, fnl, fmn2, Emn, Eml, Eln, dEdk(3)
    complex*16 :: hw, denom2, phi, dphi, gd_h
    complex*16 :: g_ln(3), g_ml(3), bracket
    real*8, parameter :: eps_deg = 1.0d-4   ! same degeneracy window as shift_shiftvector (SI Note 7)
    real*8 :: pref   ! C_ee = C_ie = 1/4 in a.u. (e=hbar=m=1, spin g dropped), plus the usual 1/(Nk*V) BZ-sum
                      ! discretisation (paper Eq. A4 and Sigma_k -> A/(2pi)^D) used throughout this module
    pref = 0.25d0/(dble(npointstotal)*vcell)

    do iw = 1, nw
      hw = hwp(iw)          ! enters g^alpha only

      ! ===================== term 2 (n != m) =====================
      do nn = 1, nband_ex
        fnn = 0.0d0; if (nn.le.nv_ex) fnn = 1.0d0
        do nnp = 1, nband_ex
          if (nnp.eq.nn) cycle
          fnnp = 0.0d0; if (nnp.le.nv_ex) fnnp = 1.0d0
          fmn2 = fnn - fnnp                     ! f_nm = f_n - f_m
          if (fmn2.eq.0.0d0) cycle
          Emn = e_nband(nnp) - e_nband(nn)
          if (abs(Emn).lt.eps_deg) cycle
          denom2 = hwsum(iw) - Emn
          phi  = 1.0d0/((hw-Emn)*Emn)
          dphi = (1.0d0/hw)*(1.0d0/(hw-Emn)**2 - 1.0d0/Emn**2)
          do njpp = 1,3
            dEdk(njpp) = dble(vme_nband(njpp,nnp,nnp) - vme_nband(njpp,nn,nn))
          end do
          do nj = 1,3
            do njp = 1,3
              do njpp = 1,3
                ! g^alpha_mn carries p^alpha_MN, i.e. (nnp,nn) -- the OPPOSITE index order to the
                ! p^lambda_nm = vme_nband(nj,nn,nnp) vertex below. Both factors must be present with
                ! opposite order for term 2 to be gauge invariant (see the header note).
                gd_h = fmn2*( dphi*dEdk(njpp)*vme_nband(njp,nnp,nn) &
                             + phi*vme_der_pt_nband(njpp,njp,nnp,nn) )
                sigma_w_sp(nj,njp,njpp,iw) = sigma_w_sp(nj,njp,njpp,iw) &
                    - pref * vme_nband(nj,nn,nnp)/denom2 * gd_h
              end do
            end do
          end do
        end do
      end do

      ! ===================== term 1 (l != n, l != m; n = m allowed) =====================
      do nn = 1, nband_ex
        do nnp = 1, nband_ex
          Emn = e_nband(nnp) - e_nband(nn)
          denom2 = hwsum(iw) - Emn
          do nl = 1, nband_ex
            if (nl.eq.nn .or. nl.eq.nnp) cycle          ! l = n or l = m: E_ln or E_ml would vanish
            Eml = e_nband(nnp) - e_nband(nl)
            Eln = e_nband(nl)  - e_nband(nn)
            if (abs(Eml).lt.eps_deg .or. abs(Eln).lt.eps_deg) cycle
            fnl  = 0.0d0; if (nl.le.nv_ex)  fnl  = 1.0d0
            fnn  = 0.0d0; if (nn.le.nv_ex)  fnn  = 1.0d0
            fnnp = 0.0d0; if (nnp.le.nv_ex) fnnp = 1.0d0
            if ((fnn-fnl).eq.0.0d0 .and. (fnl-fnnp).eq.0.0d0) cycle   ! g_ln and g_ml both zero
            do njp = 1,3
              g_ln(njp) = (fnn-fnl)*vme_nband(njp,nl,nn)/(hw-Eln)
              g_ml(njp) = (fnl-fnnp)*vme_nband(njp,nnp,nl)/(hw-Eml)
            end do
            do nj = 1,3
              do njp = 1,3
                do njpp = 1,3
                  bracket = g_ln(njp)*vme_nband(njpp,nnp,nl) - g_ml(njp)*vme_nband(njpp,nl,nn)
                  sigma_w_sp(nj,njp,njpp,iw) = sigma_w_sp(nj,njp,njpp,iw) &
                      + pref * vme_nband(nj,nn,nnp)*bracket / (Eml*Eln*denom2)
                end do
              end do
            end do
          end do
        end do
      end do

    end do
  end subroutine get_shg_intens_sp


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
    integer :: iounit100
    integer :: iounit90

    integer,     intent(in) :: nw
    real(8),     intent(in) :: wp(nw)
    complex(8),  intent(in) :: sigma_w_sp(3,3,3,nw)
    real(8),     intent(in) :: shift_vector_w(3,3,nw)

    integer :: iw
    real(8) :: feps

    ! Unit-conversion factor is loop-invariant — compute once
    feps = sigma2_au_to_si

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(newunit=iounit90,  file='shift_sp_lengthgauge_'//trim(material_name)//'.dat')
    open(newunit=iounit100, file='shift_vector.dat')

    ! --- OPTION A: simple parallel loop with ordered I/O ---
    ! Ordered writes preserve frequency ordering in the output files.
    ! serial on purpose (audit 2026-09-24): this loop is pure file I/O and every iteration was
    ! inside !$omp ordered, which serialises it completely -- the parallel wrapper only added
    ! thread spawn and synchronisation cost. HANDOFF 8.35.
    do iw = 1, nw

      write(iounit90,*) wp(iw)*27.211385d0, &
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

      write(iounit100,*) wp(iw)*27.211385d0, &
        shift_vector_w(1,1,iw), shift_vector_w(1,2,iw), shift_vector_w(1,3,iw), &
        shift_vector_w(2,1,iw), shift_vector_w(2,2,iw), shift_vector_w(2,3,iw), &
        shift_vector_w(3,1,iw), shift_vector_w(3,2,iw), shift_vector_w(3,3,iw)

    end do

    close(iounit90)
    close(iounit100)

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
    integer :: iounit92
    integer,    intent(in) :: nw
    real(8),    intent(in) :: wp(nw)
    complex(8), intent(in) :: sigma_shg(3,3,3,nw)
    integer :: iw, ia, ib, ic
    real(8) :: feps

    feps = sigma2_au_to_si

    open(newunit=iounit92, file='shg_sp_lengthgauge_'//trim(material_name)//'.dat')
    write(iounit92,'(A)') '# hbar*omega(eV) [FUNDAMENTAL, not 2*hbar*omega] | sigma^{abc}(2w;w,w) (Re, Im) in uA nm/V^2, abc = xxx,xxy,xxz,xyx,...,zzz'

    do iw = 1, nw
      write(iounit92,'(ES18.10,54ES18.10)') wp(iw)*27.211385d0, &
        ( ( ( real(feps*sigma_shg(ia,ib,ic,iw)), aimag(feps*sigma_shg(ia,ib,ic,iw)), &
              ic=1,3 ), ib=1,3 ), ia=1,3 )
    end do

    close(iounit92)

  end subroutine print_shg_second_sp
!!!!

end module sigma_second_sp

