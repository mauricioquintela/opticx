module sigma_second_ex
  use constants_math
  use parser_input_file, &
    only:e1,e2,eta,nw,response_text,broadening_type_text, &
    freq_ratio,e1b,e2b,nwb,two_freq_grid,build_freq_pairs
  use parser_wannier90_tb, &
    only:material_name
  use parser_optics_xatu_dim, &
    only:npointstotal,vcell,norb_ex_cut,nv_ex,nc_ex,nband_ex,e_ex,fk_ex
  use ome_ex, &
    only:e_ex,xme_ex,vme_ex,xme_ex_inter,vme_ex_inter,inter_terms_ready
  use sigma_second_sp, &
    only:initialize_sigma_second_arrays
  implicit none

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine get_sigma_second_ex(nwp,nwq)
    implicit none

    integer :: nwp,nwq
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    write(*,*) '11. Entering sigma_second_ex'
	write(*,*) '    Optical matrix elements (ex) will not be read from file, but passed from ome_ex'
    !compute shift conductivity
    if (nwp.eq.1 .and. nwq.eq.(-1)) then
      call get_sigma_shift_ex()
    end if

    !compute shg susceptibility (omega_p = omega_q)
    if (nwp.eq.1 .and. nwq.eq.1 .and. (response_text == 'shg' .or. response_text == 'shg_covariant')) then
      call get_sigma_shg_ex()
    end if

    ! General second order at arbitrary (omega_p, omega_q), causal broadening everywhere. get_sigma_general_ex
    ! uses Eq. (B1b) (position, method B), except where omega_1 + omega_2 = 0 is in range (rectification and 2D
    ! maps containing that line): there Eq. (B1a) (method A) with term 3 on the bare exciton current, which
    ! resolves the injection current on the mesh.
    if (response_text == 'electrooptic') then
      freq_ratio = 0.0d0; two_freq_grid = .false.
      call get_sigma_general_ex('electrooptic')
    else if (response_text == 'rectification') then
      freq_ratio = -1.0d0; two_freq_grid = .false.
      call get_sigma_general_ex('rectification')   ! whole causal sigma(0; w, -w), method A (shift current: Response = shift)
    else if (response_text == 'general') then
      call get_sigma_general_ex('general')
    end if


  end subroutine get_sigma_second_ex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Excitonic shift conductivity sigma^{abc}(0; omega, -omega).
  ! Conventions (code paper, Esteve-Paredes et al., npj Comput. Mater. 11, 13, main Eqs. 8-11 and SI Note 5):
  !  * omega_p = omega + i*eta, omega_q = -omega_p, omega_2 = 0 exactly; positive sign, sigma = +pi*e^3/(hbar V) *
  !    sum Re[S_NN'] delta(hbar*omega - E_N) (Eq. 10);
  !  * S -> Re[S] (SI Note 5: degenerate exciton pairs are not time-reversal eigenstates);
  !  * the output is symmetrised over the field indices, sigma^{abc} -> (sigma^{abc} + sigma^{acb})/2. The
  !    paper's IPA formula (Eq. 9) is written with I^{abc} + I^{acb}; the SI states that the real part is
  !    symmetric under b<->c because the current must be real, and Eq. 8 contracts sigma with
  !    eps_b eps_c, so only the b<->c-symmetric part enters for linearly polarised light. The
  !    permutation symmetrisation of Taghizadeh 2017 (App. A) reduces to this at DC.
  ! The real part of the result is printed.
  subroutine get_sigma_shift_ex()
    implicit none
    
	!here
	integer :: iw, nj, njp, njpp
    dimension :: wp(nw)
    dimension :: sigma_w_ex(3,3,3,nw)

    real*8 :: wp,eta2
    complex*16 :: sigma_w_ex
    complex(8), allocatable :: sigma_raw(:,:,:,:)   ! heap, not stack
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !initialize conductivity arrays
    
    call initialize_sigma_second_arrays(nw,wp,eta2,sigma_w_ex)
	write(*,*) '    Evaluating shift conductivity (ex)...'
    
    if (.not. inter_terms_ready) then
      write(*,*) 'ERROR (sigma_second_ex): xme_ex_inter/vme_ex_inter not '// &
                'populated — get_ome_ex must be called with iflag_norder=2 first.'
      stop 1
    end if
    call get_shift_intens_ex_matrix(wp,eta2,sigma_w_ex)
    ! symmetrise over the field indices (b,c)
    allocate(sigma_raw(3,3,3,nw))
    sigma_raw = sigma_w_ex
    do njpp = 1, 3
      do njp = 1, 3
        do nj = 1, 3
          sigma_w_ex(nj,njp,njpp,:) = 0.5d0*(sigma_raw(nj,njp,njpp,:) + sigma_raw(nj,njpp,njp,:))
        end do
      end do
    end do
    deallocate(sigma_raw)
	!print shift conductivity (ex)
	call print_sigma_second_ex(nw,wp,sigma_w_ex)
    write(*,*) '    Shift conductivity (ex) has been printed'
  end subroutine get_sigma_shift_ex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Excitonic SHG driver. Method B (Taghizadeh & Pedersen, PRB 97, 205432, Eq. B1b): position
  ! matrix elements only, frequency-chunked zgemm. Method B rather than A because vme_ex is the
  ! bare momentum P_n, not the Heisenberg momentum Pi_n = -i E_n X_n that Eq. B1a requires
  ! (see get_shg_intens_ex_matrix_methodB). Needs xme_ex_inter, i.e. the excitonic matrix
  ! elements must have been computed with OME_ex = nonlinear.
  ! The result is symmetrised over the two field indices, sigma^{abc} = sigma^{acb}: for
  ! omega_p = omega_q this is the intrinsic permutation symmetry, and the raw kernel only has it
  ! after symmetrisation (the term-3 pairing is not symmetric in b<->c).
  subroutine get_sigma_shg_ex()
    implicit none

    integer    :: nj, njp, njpp
    real(8)    :: wp(nw), eta2
    complex(8) :: sigma_shg(3,3,3,nw)
    complex(8), allocatable :: sigma_raw(:,:,:,:)   ! heap, not stack: 27*nw complex numbers

    call initialize_sigma_second_arrays(nw,wp,eta2,sigma_shg)
    write(*,*) '    Evaluating SHG susceptibility (ex)...'

    if (.not. inter_terms_ready) then
      write(*,*) 'ERROR (sigma_second_ex): xme_ex_inter/vme_ex_inter not '// &
                'populated — get_ome_ex must be called with iflag_norder=2 first.'
      stop 1
    end if

    allocate(sigma_raw(3,3,3,nw))
    call get_shg_intens_ex_matrix_methodB(wp,eta2,sigma_raw)
    do njpp = 1, 3
      do njp = 1, 3
        do nj = 1, 3
          sigma_shg(nj,njp,njpp,:) = 0.5d0*(sigma_raw(nj,njp,njpp,:) + sigma_raw(nj,njpp,njp,:))
        end do
      end do
    end do
    deallocate(sigma_raw)
    ! PHYSICAL SIGN: Eq. (B1b) as written by Taghizadeh & Pedersen already contains the electron charge,
    ! while sigma2_au_to_si applies e^3 = -1 on top; a real-time propagation of hBN (non-interacting limit) shows every
    ! causal-prescription second-order output came out with the opposite sign. Magnitudes (2018 convention) were right.
    sigma_shg = -sigma_shg
    call print_shg_second_ex(nw,wp,sigma_shg)
    write(*,*) '    SHG susceptibility (ex) has been printed'
  end subroutine get_sigma_shg_ex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  subroutine get_shift_intens_ex(wp, eta2, sigma_w_ex)
  use omp_lib
  implicit none

  real(8),    intent(in)    :: wp(nw), eta2
  complex(8), intent(inout) :: sigma_w_ex(3,3,3,nw)

  integer     :: nn, nnp, nj, njp, njpp
  integer     :: mode   ! 1 = gaussian, 2 = lorentzian
  complex(8)  :: s1, s2, s3
  complex(8), allocatable :: sigma_w_ex_t(:,:,:,:)
  complex(8), allocatable :: d1_arr(:), d2_arr(:), d3_arr(:), d4_arr(:)

  if (trim(broadening_type_text) == 'gaussian') then
    mode = 1
  else
    mode = 2   ! lorentzian, also the default fallback
  end if

  sigma_w_ex = (0.0d0, 0.0d0)

  !$omp parallel default(none) &
  !$omp   shared(mode, eta2, wp, nw, norb_ex_cut, npointstotal, vcell, sigma_w_ex) &
  !$omp   private(nn, nnp, nj, njp, njpp, s1, s2, s3, &
  !$omp           sigma_w_ex_t, d1_arr, d2_arr, d3_arr, d4_arr)

  allocate(sigma_w_ex_t(3,3,3,nw))
  sigma_w_ex_t = (0.0d0, 0.0d0)
  allocate(d1_arr(nw), d2_arr(nw), d3_arr(nw), d4_arr(nw))

  !$omp do schedule(static)
  do nn = 1, norb_ex_cut
    do nnp = 1, norb_ex_cut

      ! d1..d4 depend only on (nn,nnp) and the frequency grid wp(:) —
      ! NOT on (nj,njp,njpp). Computed once here per (nn,nnp) pair,
      ! across the whole frequency axis, instead of 27 times (once per
      ! Cartesian index triple) inside the loop below.
      call get_shift_kernel_ex_dfactors(mode, eta2, wp, nn, nnp, &
                                         d1_arr, d2_arr, d3_arr, d4_arr)

      do nj = 1, 3
        do njp = 1, 3
          do njpp = 1, 3

            call get_shift_kernel_ex_static(mode, nj, njp, njpp, nn, nnp, s1, s2, s3)

            if (mode == 1) then
              sigma_w_ex_t(nj,njp,njpp,:) = sigma_w_ex_t(nj,njp,njpp,:) &
                + ( s1*(-cmplx(0.0d0,1.0d0,8)*pi*d1_arr(:)) &
                  + s2*(-cmplx(0.0d0,1.0d0,8)*pi*d2_arr(:)) &
                  + s3*(-pi**2*d3_arr(:)*d4_arr(:)) ) &
                  / (dble(npointstotal) * vcell)
            else
              sigma_w_ex_t(nj,njp,njpp,:) = sigma_w_ex_t(nj,njp,njpp,:) &
                - ( s1*d1_arr(:) + s2*d2_arr(:) + s3*d3_arr(:) ) &
                  / (dble(npointstotal) * vcell)
            end if

          end do
        end do
      end do

    end do
  end do
  !$omp end do

  !$omp critical
    sigma_w_ex = sigma_w_ex + sigma_w_ex_t
  !$omp end critical

  deallocate(sigma_w_ex_t, d1_arr, d2_arr, d3_arr, d4_arr)
  !$omp end parallel

end subroutine get_shift_intens_ex
  
  !!!!!!!
  
  
!> Excitonic shift conductivity, zgemm-restructured form of get_shift_intens_ex.
!! Identical arithmetic (the two agree to 1e-16 in both broadening modes; guarded by
!! tests/test_shift_intens_ex_matrix.f90), but the (n,m) double sum becomes a matrix
!! product over the exciton index, which took the target run from ~2 h to ~19 min.
!! Frequencies are processed in chunks of nw_chunk to bound the workspace.
!! @param wp      Frequency grid, Hartree.
!! @param eta2    Broadening, Hartree.
!! @param sigma_w_ex  Accumulated onto, symmetrised over (b,c) by the caller.
!! @return void
subroutine get_shift_intens_ex_matrix(wp, eta2, sigma_w_ex)
  implicit none

  real(8),    intent(in)    :: wp(nw), eta2
  complex(8), intent(inout) :: sigma_w_ex(3,3,3,nw)

  integer, parameter :: nw_chunk = 3000
  integer :: nj, njp, njpp, nn, nnp
  integer :: mode
  integer :: iw0, iw1, nw_this, ichunk, nchunks
  complex(8), parameter :: ci = (0.0d0,1.0d0), czero=(0.0d0,0.0d0), cone=(1.0d0,0.0d0)

  ! ---- frequency-INDEPENDENT (lorentzian only): full norb_ex_cut length ----
  complex(8), allocatable :: Afac1(:), Afac2(:)

  ! ---- (nj,njp,njpp)-INDEPENDENT, chunk-sized, rebuilt once per chunk ----
  complex(8), allocatable :: gauss1(:,:), gauss2(:,:), gauss3(:,:), gauss4(:,:)
  complex(8), allocatable :: Bfac1(:,:), Bfac2(:,:), Cfac3(:,:), Dfac3(:,:)

  ! ---- per-(index-pair) scratch, chunk-sized, reused across every pair
  !      and every chunk ----
  complex(8), allocatable :: Mmat1(:,:), Mmat2(:,:)
  complex(8), allocatable :: Bmat1(:,:), Bmat2(:,:), Bmat3(:,:), Bmat4(:,:)
  complex(8), allocatable :: Wmat1(:,:), Wmat2(:,:), Wmat3(:,:), Wmat4(:,:)
  complex(8), allocatable :: Avec1(:), Avec2(:), Avec3(:), Avec4(:)
  complex(8), allocatable :: Amat(:,:), Amat2(:,:)
  complex(8), allocatable :: term_total(:)
  complex(8), allocatable :: PiV(:)      ! Pi_n^a = -i E_n X_n^a, rebuilt per Cartesian index

  if (trim(broadening_type_text) == 'gaussian') then
    mode = 1
  else
    mode = 2   ! lorentzian, also the default fallback
  end if

  sigma_w_ex = (0.0d0, 0.0d0)

  allocate(Mmat1(norb_ex_cut,norb_ex_cut), Mmat2(norb_ex_cut,norb_ex_cut))
  allocate(Bmat1(norb_ex_cut,nw_chunk), Bmat2(norb_ex_cut,nw_chunk))
  allocate(Bmat3(norb_ex_cut,nw_chunk), Bmat4(norb_ex_cut,nw_chunk))
  allocate(Wmat1(norb_ex_cut,nw_chunk), Wmat2(norb_ex_cut,nw_chunk))
  allocate(Wmat3(norb_ex_cut,nw_chunk), Wmat4(norb_ex_cut,nw_chunk))
  allocate(Avec1(norb_ex_cut), Avec2(norb_ex_cut), Avec3(norb_ex_cut), Avec4(norb_ex_cut))
  allocate(Amat(norb_ex_cut,nw_chunk), Amat2(norb_ex_cut,nw_chunk))
  allocate(term_total(nw_chunk), PiV(norb_ex_cut))

  if (mode == 1) then
    allocate(gauss1(norb_ex_cut,nw_chunk), gauss2(norb_ex_cut,nw_chunk))
!     allocate(gauss3(norb_ex_cut,nw_chunk), gauss4(norb_ex_cut,nw_chunk))
    allocate(gauss4(norb_ex_cut,nw_chunk))   ! gauss3 removed
  else
    allocate(Afac1(norb_ex_cut), Afac2(norb_ex_cut))
    allocate(Bfac1(norb_ex_cut,nw_chunk), Bfac2(norb_ex_cut,nw_chunk))
    !allocate(Cfac3(norb_ex_cut,nw_chunk), Dfac3(norb_ex_cut,nw_chunk))
    allocate(Dfac3(norb_ex_cut,nw_chunk))
    Afac1(:) = 1.0d0/(-e_ex(:))      ! omega_2 = 0 exactly: no eta on the omega_2 pole
    Afac2(:) = 1.0d0/( e_ex(:))
  end if

  nchunks = (nw + nw_chunk - 1) / nw_chunk

  do ichunk = 1, nchunks
    iw0     = (ichunk-1)*nw_chunk + 1
    iw1     = min(iw0 + nw_chunk - 1, nw)
    nw_this = iw1 - iw0 + 1

    ! ---- per-chunk, (nj,njp,njpp)-independent broadening factors ----
    if (mode == 1) then
      do nn = 1, norb_ex_cut
        gauss1(nn,1:nw_this) = 1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)* &
            exp(-0.5d0/(eta2**2)*(-wp(iw0:iw1)-e_ex(nn))**2)
        gauss2(nn,1:nw_this) = 1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)* &
            exp(-0.5d0/(eta2**2)*(-wp(iw0:iw1)+e_ex(nn))**2)
        !gauss3(nn,1:nw_this) = gauss2(nn,1:nw_this)
        gauss4(nn,1:nw_this) = 1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)* &
            exp(-0.5d0/(eta2**2)*( wp(iw0:iw1)-e_ex(nn))**2)
      end do
    else
      do nn = 1, norb_ex_cut
        ! omega_q = -(omega + i*eta): all poles carry z = omega + i*eta (see get_shift_kernel_ex_dfactors)
        Bfac1(nn,1:nw_this) = 1.0d0/(-wp(iw0:iw1)-e_ex(nn)-cmplx(0.0d0,eta2,8))
        Bfac2(nn,1:nw_this) = 1.0d0/(-wp(iw0:iw1)+e_ex(nn)-cmplx(0.0d0,eta2,8))
        Dfac3(nn,1:nw_this) = 1.0d0/( wp(iw0:iw1)-e_ex(nn)+cmplx(0.0d0,eta2,8))
        !Cfac3(nn,1:nw_this) = 1.0d0/( e_ex(nn)-wp(iw0:iw1)+cmplx(0.0d0,eta2,8)) !entirely identical to Bfac2
      end do
    end if

    ! =================================================================
    ! TERM 1 + TERM 2: the zgemm depends only on (njp,njpp). Loop those
    ! outermost; nj enters only in the cheap reduction below.
    ! =================================================================
    do njp = 1, 3
      do njpp = 1, 3

        if (mode == 1) then
          do nnp = 1, norb_ex_cut
            Bmat1(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) * (-ci*pi) * gauss1(nnp,1:nw_this)
            Bmat2(nnp,1:nw_this) = xme_ex(njpp,nnp)        * (-ci*pi) * gauss2(nnp,1:nw_this)
          end do
          Mmat1 = xme_ex_inter(njp,:,:)
          Mmat2 = conjg(Mmat1)

          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat1, norb_ex_cut, &
                      Bmat1, norb_ex_cut, czero, Wmat1, norb_ex_cut)
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat2, norb_ex_cut, &
                      Bmat2, norb_ex_cut, czero, Wmat2, norb_ex_cut)

          do nj = 1, 3
            ! Pi_n = -i E_n X_n  =>  -Pi_n/E_n = i X_n  and  conj(Pi_n)/E_n = i X_n^*
            Avec1(:) = ci*xme_ex(nj,:)
            Avec2(:) = ci*conjg(xme_ex(nj,:))
            term_total(1:nw_this) = matmul(Avec1, Wmat1(:,1:nw_this)) &
                                   + matmul(Avec2, Wmat2(:,1:nw_this))
            sigma_w_ex(nj,njp,njpp,iw0:iw1) = sigma_w_ex(nj,njp,njpp,iw0:iw1) &
                + term_total(1:nw_this) / (dble(npointstotal)*vcell)
          end do

        else   ! lorentzian
          do nnp = 1, norb_ex_cut
            Bmat1(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) * Bfac1(nnp,1:nw_this)
            Bmat2(nnp,1:nw_this) = xme_ex(njpp,nnp)        * Bfac1(nnp,1:nw_this)
            Bmat3(nnp,1:nw_this) = xme_ex(njpp,nnp)        * Bfac2(nnp,1:nw_this)
            Bmat4(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) * Bfac2(nnp,1:nw_this)
          end do
          Mmat1 = xme_ex_inter(njp,:,:)
          Mmat2 = conjg(Mmat1)

          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat1, norb_ex_cut, &
                      Bmat1, norb_ex_cut, czero, Wmat1, norb_ex_cut)
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat2, norb_ex_cut, &
                      Bmat2, norb_ex_cut, czero, Wmat2, norb_ex_cut)
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat2, norb_ex_cut, &
                      Bmat3, norb_ex_cut, czero, Wmat3, norb_ex_cut)
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat1, norb_ex_cut, &
                      Bmat4, norb_ex_cut, czero, Wmat4, norb_ex_cut)

          do nj = 1, 3
            PiV(:)   = -ci*e_ex(:)*xme_ex(nj,:)          ! Pi_n^a = -i E_n X_n^a
            Avec1(:) = PiV(:)        * Afac1(:)   ! term1, A
            Avec2(:) = conjg(PiV(:)) * Afac1(:)   ! term1, A*
            Avec3(:) = conjg(PiV(:)) * Afac2(:)   ! term2, A
            Avec4(:) = PiV(:)        * Afac2(:)   ! term2, A*

            ! (s - s*)/2 = i*Im(s): the S -> Re[S] projection (see get_shift_kernel_ex_static)
            term_total(1:nw_this) = &
                ( matmul(Avec1, Wmat1(:,1:nw_this)) - matmul(Avec2, Wmat2(:,1:nw_this)) &
                + matmul(Avec3, Wmat3(:,1:nw_this)) - matmul(Avec4, Wmat4(:,1:nw_this)) ) / 2.0d0

            sigma_w_ex(nj,njp,njpp,iw0:iw1) = sigma_w_ex(nj,njp,njpp,iw0:iw1) &
                - term_total(1:nw_this) / (dble(npointstotal)*vcell)
          end do
        end if

      end do
    end do

    ! =================================================================
    ! TERM 3: the zgemm depends only on (nj,njpp). Loop those outermost;
    ! njp enters only in the cheap reduction below.
    ! =================================================================
    do nj = 1, 3
      do njpp = 1, 3

        if (mode == 1) then
          do nnp = 1, norb_ex_cut
            Bmat1(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) * gauss4(nnp,1:nw_this)
          end do
          do nnp = 1, norb_ex_cut
            Mmat1(:,nnp) = ci*(e_ex(:)-e_ex(nnp))*xme_ex_inter(nj,:,nnp)   ! Pi_nm = i (E_n-E_m) X_nm
          end do
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat1, norb_ex_cut, &
                      Bmat1, norb_ex_cut, czero, Wmat1, norb_ex_cut)

          do njp = 1, 3
            do nn = 1, norb_ex_cut
              !Amat(nn,1:nw_this) = -xme_ex(njp,nn) * (-pi**2) * gauss3(nn,1:nw_this) !gauss3 is simply a copy of gauss2
              Amat(nn,1:nw_this) = -xme_ex(njp,nn) * (-pi**2) * gauss2(nn,1:nw_this)
            end do
            term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat1(:,1:nw_this), dim=1)
            sigma_w_ex(nj,njp,njpp,iw0:iw1) = sigma_w_ex(nj,njp,njpp,iw0:iw1) &
                + term_total(1:nw_this) / (dble(npointstotal)*vcell)
          end do

        else   ! lorentzian
          do nnp = 1, norb_ex_cut
            Bmat1(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) * Dfac3(nnp,1:nw_this)
            Bmat2(nnp,1:nw_this) = xme_ex(njpp,nnp)        * Dfac3(nnp,1:nw_this)
          end do
          do nnp = 1, norb_ex_cut
            Mmat1(:,nnp) = ci*(e_ex(:)-e_ex(nnp))*xme_ex_inter(nj,:,nnp)   ! Pi_nm = i (E_n-E_m) X_nm
          end do
          Mmat2 = conjg(Mmat1)

          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat1, norb_ex_cut, &
                      Bmat1, norb_ex_cut, czero, Wmat1, norb_ex_cut)
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat2, norb_ex_cut, &
                      Bmat2, norb_ex_cut, czero, Wmat2, norb_ex_cut)

          do njp = 1, 3
            do nn = 1, norb_ex_cut
              !Amat(nn,1:nw_this)  = -xme_ex(njp,nn)        * Cfac3(nn,1:nw_this)
              Amat(nn,1:nw_this)  = -xme_ex(njp,nn)        * Bfac2(nn,1:nw_this)
              !Amat2(nn,1:nw_this) = -conjg(xme_ex(njp,nn)) * Cfac3(nn,1:nw_this) !Cfac3 is Bfac2
              Amat2(nn,1:nw_this) = -conjg(xme_ex(njp,nn)) * Bfac2(nn,1:nw_this)
            end do
            term_total(1:nw_this) = &
                ( sum(Amat(:,1:nw_this)*Wmat1(:,1:nw_this), dim=1) &
                - sum(Amat2(:,1:nw_this)*Wmat2(:,1:nw_this), dim=1) ) / 2.0d0
            sigma_w_ex(nj,njp,njpp,iw0:iw1) = sigma_w_ex(nj,njp,njpp,iw0:iw1) &
                - term_total(1:nw_this) / (dble(npointstotal)*vcell)
          end do
        end if

      end do
    end do

  end do   ! ichunk

  if (mode == 1) then
    !deallocate(gauss1, gauss2, gauss3, gauss4)
    deallocate(gauss1, gauss2, gauss4)
  else
    !deallocate(Afac1, Afac2, Bfac1, Bfac2, Cfac3, Dfac3)
    deallocate(Afac1, Afac2, Bfac1, Bfac2, Dfac3)
  end if
  deallocate(Mmat1, Mmat2, Bmat1, Bmat2, Bmat3, Bmat4)
  deallocate(Wmat1, Wmat2, Wmat3, Wmat4)
  deallocate(Avec1, Avec2, Avec3, Avec4, Amat, Amat2, term_total, PiV)

end subroutine get_shift_intens_ex_matrix
  
  !!!!
  
  ! Frequency-independent part: s1, s2, s3 only.
  ! Frequency-independent matrix-element products of the excitonic shift conductivity
  ! (Esteve-Paredes et al., npj Comput. Mater. 11, 13, SI Eq. 9), written with POSITION matrix
  ! elements only (equivalent to method B of Taghizadeh & Pedersen, PRB 97, 205432):
  !   V_0N  -> Pi_n  = -i E_n X_n        (Heisenberg momentum, NOT the bare momentum in vme_ex)
  !   V_NN' -> Pi_nm =  i (E_n - E_m) X_nm
  ! (the SI uses exactly these identities to reduce Eq. 9 to the position-only Eq. 10).
  ! Under time reversal s1,s2,s3 are purely imaginary (s = -i E S with S real). In the
  ! lorentzian branch only the imaginary part of s is kept, as i*Im(s), i.e. S -> Re[S] (the
  ! SI's prescription); multiplied by the complex pole factor this gives the absorptive
  ! (delta-like) line shape when the real part is printed. The gaussian branch, which
  ! multiplies by -i*pi*delta, needs no projection.
  !> The frequency-INDEPENDENT factors of the three terms of the excitonic shift kernel,
  !! i.e. Eq. (B1a) of Taghizadeh & Pedersen, PRB 97, 205432 (2018) with Pi built from X
  !! (Pi_n = -i E_n X_n, Pi_nm = i(E_n - E_m) X_nm) -- the identities the code paper's SI
  !! uses to get its Eq. (10) from Eq. (9).
  !! In the Lorentzian branch only i*Im(s) is kept: the S -> Re[S] projection of SI Note 5,
  !! which is what makes the line shape absorptive rather than dispersive.
  !! Split out from the frequency loop because these depend only on (n,m) and the Cartesian
  !! triple, so they are evaluated once per exciton pair instead of once per frequency.
  !! @param mode              1 = gaussian, 2 = lorentzian.
  !! @param nj, njp, njpp     Cartesian indices a, b, c.
  !! @param nn, nnp           Exciton indices n, m.
  !! @param s1, s2, s3        The three static factors.
  !! @return void
  subroutine get_shift_kernel_ex_static(mode, nj, njp, njpp, nn, nnp, s1, s2, s3)
    implicit none
    integer,    intent(in)  :: mode, nj, njp, njpp, nn, nnp
    complex(8), intent(out) :: s1, s2, s3
    integer :: nj1, nj2, nj3
    complex(8), parameter :: ci = (0.0d0,1.0d0)
    complex(8) :: pi_n, pi_nm

    nj1 = nj; nj2 = njp; nj3 = njpp

    pi_n  = -ci*e_ex(nn)*xme_ex(nj1,nn)                                  ! Pi_n^a = -i E_n X_n^a
    pi_nm =  ci*(e_ex(nn)-e_ex(nnp))*xme_ex_inter(nj1,nn,nnp)            ! Pi_nm^a = i (E_n-E_m) X_nm^a

    if (mode == 1) then   ! gaussian
      s1 = -pi_n/e_ex(nn) * xme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))         ! a→V_0N, b→R_NN'
      s2 =  conjg(pi_n)/e_ex(nn) * conjg(xme_ex_inter(nj2,nn,nnp)) * xme_ex(nj3,nnp)
      s3 = -xme_ex(nj2,nn) * pi_nm * conjg(xme_ex(nj3,nnp))

    else                    ! lorentzian
      s1 = pi_n * xme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))
      s2 = conjg(pi_n) * conjg(xme_ex_inter(nj2,nn,nnp)) * xme_ex(nj3,nnp)
      s3 = -xme_ex(nj2,nn) * pi_nm * conjg(xme_ex(nj3,nnp))
      ! keep only i*Im(s)  (= (s - s*)/2): the S -> Re[S] projection; frequency-independent,
      ! so it is folded in here rather than redone inside the iw loop.
      s1 = cmplx(0.0d0, aimag(s1), 8)
      s2 = cmplx(0.0d0, aimag(s2), 8)
      s3 = cmplx(0.0d0, aimag(s3), 8)
    end if
  end subroutine get_shift_kernel_ex_static
  
  subroutine get_shift_kernel_ex_dfactors(mode, eta2, wp, nn, nnp, d1, d2, d3, d4)
  implicit none
  integer,    intent(in)  :: mode, nn, nnp
  real(8),    intent(in)  :: eta2, wp(nw)
  complex(8), intent(out) :: d1(nw), d2(nw), d3(nw), d4(nw)

  if (mode == 1) then   ! gaussian
    ! omegap = wp(:), omegaq = -wp(:), omega2 = 0.0d0  (as in the original)
    d1(:) = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(-wp(:)-e_ex(nnp))**2)
    d2(:) = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(-wp(:)+e_ex(nnp))**2)
    d3(:) = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(-wp(:)+e_ex(nn))**2)
    d4(:) = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(wp(:)-e_ex(nnp))**2)
  else                    ! lorentzian
    ! DC shift, paper convention: omega_p = omega + i*eta, omega_q = -omega_p, so omega_2 = 0
    ! EXACTLY (no eta on omega_2) and every pole carries the same complex frequency z = omega + i*eta:
    ! (omega_q -+ E) = -(z +- E). This makes term 3 vanish under time reversal (d3 symmetric under
    ! n <-> n') and gives +pi*sum Re[S]*delta_eta(omega-E_N), the sign of Eq. (10) of the code paper.
    d1(:) = 1.0d0 / ( cmplx(-e_ex(nn),0.0d0,8) * (-wp(:)-e_ex(nnp)-cmplx(0.0d0,eta2,8)) )
    d2(:) = 1.0d0 / ( cmplx( e_ex(nn),0.0d0,8) * (-wp(:)+e_ex(nnp)-cmplx(0.0d0,eta2,8)) )
    d3(:) = 1.0d0 / ( ( e_ex(nn)-wp(:)-cmplx(0.0d0,eta2,8)) * ( wp(:)-e_ex(nnp)+cmplx(0.0d0,eta2,8)) )
    d4(:) = (0.0d0, 0.0d0)   ! unused in lorentzian mode; present only so the interface is uniform
  end if

end subroutine get_shift_kernel_ex_dfactors
  
  ! Frequency-dependent part: d1..d4 and their combination with s1,s2,s3.
subroutine get_shift_kernel_ex_freq(mode, eta2, omegap, omegaq, omega2, &
                                     nn, nnp, s1, s2, s3, shift_kernel_ex)
  implicit none
  integer,    intent(in)  :: mode, nn, nnp
  real(8),    intent(in)  :: eta2, omegap, omegaq, omega2
  complex(8), intent(in)  :: s1, s2, s3
  complex(8), intent(out) :: shift_kernel_ex
  complex(8) :: d1, d2, d3, d4, aux1, aux2, aux3

  if (mode == 1) then   ! gaussian
    d1 = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(omegaq-e_ex(nnp))**2)
    d2 = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(omegaq+e_ex(nnp))**2)
    d3 = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(omegaq+e_ex(nn))**2)
    d4 = 1.0d0/eta2 * 1.0d0/sqrt(2.0d0*pi) * exp(-0.5d0/(eta2**2)*(omegap-e_ex(nnp))**2)

    ! sign chosen so that the resonant term is +pi*Re[S]*delta (paper Eq. 10); the lorentzian
    ! branch gets the same sign from omega_q = -omega_p - i*eta (see get_shift_kernel_ex_dfactors)
    aux1 = s1 * (-complex(0.0d0,1.0d0)*pi*d1)
    aux2 = s2 * (-complex(0.0d0,1.0d0)*pi*d2)
    aux3 = s3 * (-pi**2*d3*d4)
    shift_kernel_ex = +(aux1+aux2+aux3)
  else                    ! lorentzian
    ! omega_p = omegap + i*eta, omega_q = -omega_p (needs omegaq = -omegap), omega_2 = omega2 = 0 exactly
    d1 = 1.0d0 / ((omega2-e_ex(nn)) * (omegaq-e_ex(nnp)-cmplx(0.0d0,eta2,8)))
    d2 = 1.0d0 / ((omega2+e_ex(nn)) * (omegaq+e_ex(nnp)-cmplx(0.0d0,eta2,8)))
    d3 = 1.0d0 / ((omegaq+e_ex(nn)-cmplx(0.0d0,eta2,8)) * (omegap-e_ex(nnp)+cmplx(0.0d0,eta2,8)))

    aux1 = s1 * d1
    aux2 = s2 * d2
    aux3 = s3 * d3
    shift_kernel_ex = -(aux1+aux2+aux3)
  end if

end subroutine get_shift_kernel_ex_freq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! SHG 


! Frequency-independent matrix-element products. No mode branching --
! the SHG response stays fully complex, so there's no aimag()-forcing
! step the way the excitonic DC Lorentzian branch needed.
! NOTE (method A, Eq. B1a): this takes Pi_n and Pi_nm from vme_ex / vme_ex_inter, which on
! real data are the BARE momentum elements P_n, P_nm (Pi_n = P_n - i F_n, Eq. 10 of the paper),
! so on real data it does not reproduce method B. It is correct only if the caller supplies
! Pi_n = -i E_n X_n and Pi_nm = i (E_n - E_m) X_nm. Production runs use method B.
subroutine get_shg_kernel_ex_static(nj, njp, njpp, nn, nnp, s1, s2, s3)
  implicit none
  integer,    intent(in)  :: nj, njp, njpp, nn, nnp
  complex(8), intent(out) :: s1, s2, s3

  s1 =  vme_ex(nj,nn)        * xme_ex_inter(njp,nn,nnp)        * conjg(xme_ex(njpp,nnp))
  s2 =  conjg(vme_ex(nj,nn)) * conjg(xme_ex_inter(njp,nn,nnp)) * xme_ex(njpp,nnp)
  s3 = -xme_ex(njpp,nn)      * vme_ex_inter(nj,nn,nnp)         * conjg(xme_ex(njp,nnp))   ! Eq. (A9b) pairing

end subroutine get_shg_kernel_ex_static

! Frequency-dependent combination, evaluated at a single omega (SHG: ωp=ωq=omega).
subroutine get_shg_kernel_ex_freq(eta2, omega, nn, nnp, s1, s2, s3, shg_kernel)
  implicit none
  real(8),    intent(in)  :: eta2, omega
  integer,    intent(in)  :: nn, nnp
  complex(8), intent(in)  :: s1, s2, s3
  complex(8), intent(out) :: shg_kernel
  complex(8) :: om2c, omqc, ompc, d1, d2, d3

  omqc = cmplx(omega, eta2, 8)         ! omega_q + i*eta
  ompc = cmplx(omega, eta2, 8)         ! omega_p + i*eta (SHG: numerically == omega_q, kept distinct for clarity)
  om2c = ompc + omqc                   ! omega_2 = omega_p + omega_q -> 2*omega + 2i*eta (PRB 97, 205432: omega -> omega+i*eta for every frequency)

  d1 = 1.0d0 / ( (om2c - e_ex(nn)) * (omqc - e_ex(nnp)) )
  d2 = 1.0d0 / ( (om2c + e_ex(nn)) * (omqc + e_ex(nnp)) )
  d3 = 1.0d0 / ( (omqc + e_ex(nn)) * (ompc - e_ex(nnp)) )

  shg_kernel = -( s1*d1 + s2*d2 + s3*d3 )

end subroutine get_shg_kernel_ex_freq

subroutine get_shg_intens_ex(wp, eta2, sigma_shg)
  implicit none
  real(8),    intent(in)    :: wp(nw), eta2
  complex(8), intent(inout) :: sigma_shg(3,3,3,nw)

  integer     :: nj, njp, njpp, nn, nnp, iw
  complex(8)  :: s1, s2, s3, kernel

  sigma_shg = (0.0d0, 0.0d0)

  do nn = 1, norb_ex_cut
    do nnp = 1, norb_ex_cut
      do nj = 1, 3
        do njp = 1, 3
          do njpp = 1, 3
            call get_shg_kernel_ex_static(nj, njp, njpp, nn, nnp, s1, s2, s3)
            do iw = 1, nw
              call get_shg_kernel_ex_freq(eta2, wp(iw), nn, nnp, s1, s2, s3, kernel)
              ! the kernel already carries the leading sign of Eq. (B1a), so accumulate it with +
              sigma_shg(nj,njp,njpp,iw) = sigma_shg(nj,njp,njpp,iw) &
                  + kernel / (dble(npointstotal)*vcell)
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine get_shg_intens_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_shg_kernel_ex_static_methodB(nj, njp, njpp, nn, nnp, s1, s2, s3)
  implicit none
  integer,    intent(in)  :: nj, njp, njpp, nn, nnp
  complex(8), intent(out) :: s1, s2, s3

  s1 = xme_ex(nj,nn)        * xme_ex_inter(njp,nn,nnp)        * conjg(xme_ex(njpp,nnp))
  s2 = conjg(xme_ex(nj,nn)) * conjg(xme_ex_inter(njp,nn,nnp)) * xme_ex(njpp,nnp)
  s3 = xme_ex(njpp,nn)      * xme_ex_inter(nj,nn,nnp)         * conjg(xme_ex(njp,nnp))   ! Eq. (A9b) pairing

end subroutine get_shg_kernel_ex_static_methodB

subroutine get_shg_kernel_ex_freq_methodB(eta2, omega, nn, nnp, s1, s2, s3, shg_kernel)
  implicit none
  real(8),    intent(in)  :: eta2, omega
  integer,    intent(in)  :: nn, nnp
  complex(8), intent(in)  :: s1, s2, s3
  complex(8), intent(out) :: shg_kernel
  complex(8) :: omega_c, omega_2c, d1, d2, d3

  omega_c  = cmplx(omega, eta2, 8)         ! omega_p = omega_q, + i*eta
  omega_2c = omega_c + omega_c             ! omega_2 = omega_p + omega_q -> 2*omega + 2i*eta (must match method A)

  d1 = 1.0d0 / ( (omega_2c - e_ex(nn)) * (omega_c - e_ex(nnp)) )
  d2 = 1.0d0 / ( (omega_2c + e_ex(nn)) * (omega_c + e_ex(nnp)) )
  d3 = 1.0d0 / ( (omega_c  + e_ex(nn)) * (omega_c - e_ex(nnp)) )

  shg_kernel = cmplx(0.0d0,1.0d0,8)*omega_2c * ( s1*d1 + s2*d2 - s3*d3 )

end subroutine get_shg_kernel_ex_freq_methodB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Method A (Eq. B1a), zgemm version. Same caveat as get_shg_kernel_ex_static: it reads
! Pi from vme_ex/vme_ex_inter. Use get_shg_intens_ex_matrix_methodB for production.
subroutine get_shg_intens_ex_matrix(wp, eta2, sigma_shg)
  implicit none
  real(8),    intent(in)    :: wp(nw), eta2
  complex(8), intent(inout) :: sigma_shg(3,3,3,nw)

  integer, parameter :: nw_chunk = 2000
  integer :: nj, njp, njpp, nn, nnp
  integer :: iw0, iw1, nw_this, ichunk, nchunks
  complex(8), parameter :: czero=(0.0d0,0.0d0), cone=(1.0d0,0.0d0)

  complex(8), allocatable :: omega_c(:), omega2_c(:)         ! (nw_chunk) omega+i*eta and omega_2 = 2*(omega+i*eta)
  complex(8), allocatable :: Mmat(:,:)
  complex(8), allocatable :: Bmat(:,:), Wmat(:,:), Amat(:,:)
  complex(8), allocatable :: term_total(:)

  sigma_shg = (0.0d0, 0.0d0)

  allocate(omega_c(nw_chunk), omega2_c(nw_chunk))
  allocate(Mmat(norb_ex_cut,norb_ex_cut))
  allocate(Bmat(norb_ex_cut,nw_chunk), Wmat(norb_ex_cut,nw_chunk), Amat(norb_ex_cut,nw_chunk))
  allocate(term_total(nw_chunk))

  nchunks = (nw + nw_chunk - 1) / nw_chunk

  do ichunk = 1, nchunks
    iw0     = (ichunk-1)*nw_chunk + 1
    iw1     = min(iw0 + nw_chunk - 1, nw)
    nw_this = iw1 - iw0 + 1

    omega_c(1:nw_this)  = cmplx(    wp(iw0:iw1), eta2, 8)   ! omega_p and omega_q (SHG: same value)
    omega2_c(1:nw_this) = omega_c(1:nw_this) + omega_c(1:nw_this)   ! omega_2 = omega_p + omega_q -> 2*omega + 2i*eta

    ! ============ TERM 1: zgemm depends only on (njp,njpp) ============
    do njp = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) / (omega_c(1:nw_this) - e_ex(nnp))
        end do
        Mmat = xme_ex_inter(njp,:,:)
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)

        do nj = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = vme_ex(nj,nn) / (omega2_c(1:nw_this) - e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          sigma_shg(nj,njp,njpp,iw0:iw1) = sigma_shg(nj,njp,njpp,iw0:iw1) &
              - term_total(1:nw_this) / (dble(npointstotal)*vcell)
        end do
      end do
    end do

    ! ============ TERM 2: zgemm depends only on (njp,njpp) ============
    do njp = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = xme_ex(njpp,nnp) / (omega_c(1:nw_this) + e_ex(nnp))
        end do
        Mmat = conjg(xme_ex_inter(njp,:,:))
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)

        do nj = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = conjg(vme_ex(nj,nn)) / (omega2_c(1:nw_this) + e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          sigma_shg(nj,njp,njpp,iw0:iw1) = sigma_shg(nj,njp,njpp,iw0:iw1) &
              - term_total(1:nw_this) / (dble(npointstotal)*vcell)
        end do
      end do
    end do

    ! ============ TERM 3: unaffected, ω2 never appears here ============
    do nj = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) / (omega_c(1:nw_this) - e_ex(nnp))
        end do
        Mmat = vme_ex_inter(nj,:,:)
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)

        do njp = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = -xme_ex(njp,nn) / (omega_c(1:nw_this) + e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          ! Eq. (A9b) field pairing, as in get_second_intens_ex_methodA: stored at (a, b = njpp, c = njp)
          sigma_shg(nj,njpp,njp,iw0:iw1) = sigma_shg(nj,njpp,njp,iw0:iw1) &
              - term_total(1:nw_this) / (dble(npointstotal)*vcell)
        end do
      end do
    end do

  end do

  deallocate(omega_c, omega2_c, Mmat, Bmat, Wmat, Amat, term_total)

end subroutine get_shg_intens_ex_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Method B (Taghizadeh & Pedersen, PRB 97, 205432, Eq. B1b), frequency-chunked zgemm
! version. Uses ONLY the position matrix elements xme_ex, xme_ex_inter (no velocity
! elements), so it does not depend on what vme_ex holds. This is the production
! SHG routine: vme_ex/vme_ex_inter are the bare momentum P_n/P_nm, not the Heisenberg
! momentum Pi_n=-iE_nX_n needed by Eq. B1a, so method A on those arrays is not the paper's
! method A.
!   sigma^{abc} = (i*omega_2/(Nk*V)) * sum_{nm} [  X_n^a X_nm^b X_m^c* /((w2-E_n)(w-E_m))
!                                                + X_n^a* X_nm^b* X_m^c /((w2+E_n)(w+E_m))
!                                                - X_n^b X_nm^a X_m^c* /((w+E_n)(w-E_m)) ]
! with w = omega+i*eta and omega_2 = w + w = 2*omega + 2i*eta (paper convention).
! Same conventions as get_shg_kernel_ex_*_methodB, which this routine must reproduce.
! GENERALISED 2026-09-24 to arbitrary (omega_p, omega_q); see.
! Taghizadeh & Pedersen PRB 97, 205432 (2018) Eq. (B1b), method B (position elements only):
!
!   sigma^B(2) = +C_ee (i hbar w2) sum_nm [  X_n X_nm X*_m / ((hw2 - E_n)(hw_q - E_m))
!                                          + X*_n X*_nm X_m / ((hw2 + E_n)(hw_q + E_m))
!                                          - X_n X_nm X*_m / ((hw_q + E_n)(hw_p - E_m)) ]
!
! Unlike the single-particle Eq. (A3a) -- where omega_q only ever entered through the sum
! -- here omega_p appears ON ITS OWN, in the third term. So this kernel needs THREE complex frequencies:
! hwp, hwq and hw2 = hwp + hwq. Note also that the third term of Eq. (A9b) is U_n O_nm U*_m, i.e. the
! OBSERVABLE sits on the inter-exciton element there while in terms 1 and 2 it sits on the n element;
! the existing index assignment (Mmat = xme_ex_inter(nj) in term 3 but (njp) in terms 1-2) is therefore
! correct and is NOT a copy-paste slip.
!
! IMPORTANT: the prefactor is i*hbar*w2, so at w2 = 0 this expression vanishes IDENTICALLY, and on the causal
! DC line (w2 = 2i eta) it resolves the injection current only with the k-mesh. The general driver uses
! method A there. The shift current in the convention of npj Comput. Mater. 11, 13 (2025) (omega_q = -omega_p,
! omega_2 = 0, Eqs. 10-11) is get_shift_intens_ex.
subroutine get_second_intens_ex_methodB(nfreq, hwp, hwq, hw2, sigma_shg)
  implicit none
  integer,    intent(in)    :: nfreq
  complex(8), intent(in)    :: hwp(nfreq), hwq(nfreq), hw2(nfreq)
  complex(8), intent(inout) :: sigma_shg(3,3,3,nfreq)

  integer, parameter :: nw_chunk = 2000
  integer :: nj, njp, njpp, nn, nnp
  integer :: iw0, iw1, nw_this, ichunk, nchunks
  complex(8), parameter :: ci=(0.0d0,1.0d0), czero=(0.0d0,0.0d0), cone=(1.0d0,0.0d0)

  complex(8), allocatable :: omega_c(:), omega2_c(:), omegap_c(:), pref(:)
  complex(8), allocatable :: Mmat(:,:)
  complex(8), allocatable :: Bmat(:,:), Wmat(:,:), Amat(:,:)
  complex(8), allocatable :: term_total(:)

  sigma_shg = (0.0d0, 0.0d0)

  allocate(omega_c(nw_chunk), omega2_c(nw_chunk), omegap_c(nw_chunk), pref(nw_chunk))
  allocate(Mmat(norb_ex_cut,norb_ex_cut))
  allocate(Bmat(norb_ex_cut,nw_chunk), Wmat(norb_ex_cut,nw_chunk), Amat(norb_ex_cut,nw_chunk))
  allocate(term_total(nw_chunk))

  nchunks = (nfreq + nw_chunk - 1) / nw_chunk

  do ichunk = 1, nchunks
    iw0     = (ichunk-1)*nw_chunk + 1
    iw1     = min(iw0 + nw_chunk - 1, nfreq)
    nw_this = iw1 - iw0 + 1

    omega_c(1:nw_this)  = hwq(iw0:iw1)      ! omega_q: terms 1,2 (the X*_m/X_m pole) and term 3's +E_n
    omegap_c(1:nw_this) = hwp(iw0:iw1)      ! omega_p: term 3's (hw_p - E_m) only
    omega2_c(1:nw_this) = hw2(iw0:iw1)      ! omega_2 = omega_p + omega_q
    pref(1:nw_this)     = ci*omega2_c(1:nw_this) / (dble(npointstotal)*vcell)

    ! ============ TERM 1: zgemm depends only on (njp=b, njpp=c) ============
    do njp = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) / (omega_c(1:nw_this) - e_ex(nnp))
        end do
        Mmat = xme_ex_inter(njp,:,:)
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)
        do nj = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = xme_ex(nj,nn) / (omega2_c(1:nw_this) - e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          sigma_shg(nj,njp,njpp,iw0:iw1) = sigma_shg(nj,njp,njpp,iw0:iw1) &
              + pref(1:nw_this)*term_total(1:nw_this)
        end do
      end do
    end do

    ! ============ TERM 2: zgemm depends only on (njp=b, njpp=c) ============
    do njp = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = xme_ex(njpp,nnp) / (omega_c(1:nw_this) + e_ex(nnp))
        end do
        Mmat = conjg(xme_ex_inter(njp,:,:))
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)
        do nj = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = conjg(xme_ex(nj,nn)) / (omega2_c(1:nw_this) + e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          sigma_shg(nj,njp,njpp,iw0:iw1) = sigma_shg(nj,njp,njpp,iw0:iw1) &
              + pref(1:nw_this)*term_total(1:nw_this)
        end do
      end do
    end do

    ! ============ TERM 3 (enters with a minus sign): zgemm depends only on (nj=a, njpp=c) ============
    do nj = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          ! Eq. (B1b) term 3 denominator is (hbar*omega_P - E_m) -- omega_p, NOT omega_q.
          Bmat(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) / (omegap_c(1:nw_this) - e_ex(nnp))
        end do
        Mmat = xme_ex_inter(nj,:,:)
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)
        do njp = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = xme_ex(njp,nn) / (omega_c(1:nw_this) + e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          ! FIELD PAIRING (Eq. A9b): in term 3, U_n carries the w_q field and U*_m the w_p field, so the X_n
          ! factor (with hw_q + E_n) is the c index and the X*_m factor (with hw_p - E_m) the b index. Stored at
          ! (a, b = njpp, c = njp). Only the b<->c antisymmetric part depends on this: with the indices the other way
          ! round the injection current (circular light, buckled hBN) had the opposite sign to a real-time propagation.
          sigma_shg(nj,njpp,njp,iw0:iw1) = sigma_shg(nj,njpp,njp,iw0:iw1) &
              - pref(1:nw_this)*term_total(1:nw_this)
        end do
      end do
    end do

  end do

  deallocate(omega_c, omega2_c, omegap_c, pref, Mmat, Bmat, Wmat, Amat, term_total)

end subroutine get_second_intens_ex_methodB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENERAL excitonic second order via METHOD A, Taghizadeh & Pedersen PRB 97, 205432 (2018) Eq. (B1a):
!
!   sigma^A(2) = -C_ee sum_nm [  Pi_n X_nm X*_m / ((hw2 - E_n)(hw_q - E_m))
!                              + Pi*_n X*_nm X_m / ((hw2 + E_n)(hw_q + E_m))
!                              - X_n Pi_nm X*_m / ((hw_q + E_n)(hw_p - E_m)) ]
!
! WHY THIS EXISTS: unlike Eq. (B1b), method A has NO i*hbar*omega_2 prefactor, so it does
! NOT vanish at omega_2 = 0. With term 3 on the bare current (bare_term3) it is the route to the EXCITONIC
! OPTICAL RECTIFICATION on the causal DC line, injection current included.
!
! Pi is built from X INSIDE this routine, Pi_n = -i E_n X_n and Pi_nm = i(E_n - E_m) X_nm (paper Eq. 4a/4b;
! Taghizadeh & Pedersen 2018, Eq. 10: with an e-h interaction Pi_n = P_n - i F_n /= P_n). It deliberately does
! NOT read vme_ex: the ground-to-exciton BARE momentum P_n differs from Pi_n by ~34% on hBN, and terms 1-2 need
! Pi. Term 3 reads the bare exciton-exciton current vme_ex_inter only when bare_term3 is set (see below).
!
! Index roles follow Eq. (A9b) exactly, as in method B: in terms 1 and 2 the observable sits on the n
! element (Amat carries a, Mmat carries b); in term 3 it sits on the INTER-exciton element (Mmat carries a,
! Amat carries b).
!> Excitonic second-order conductivity by method A, Eq. (B1a) of Taghizadeh & Pedersen,
!! PRB 97, 205432 (2018): the Heisenberg-momentum (Pi) observable. Pi is built internally
!! from X, never read from vme_ex, which holds the bare momentum P.
!! Used where omega_1 + omega_2 = 0 is in range: Eq. (B1b) carries an overall i*hbar*omega_2 (= -2 eta there,
!! causal broadening), and its term 3 resolves the injection current only with the k-mesh.
!! @param nfreq             Number of (omega_p, omega_q) pairs.
!! @param hwp, hwq, hw2     Complex hbar*omega_p, hbar*omega_q and their sum.
!! @param sigma_shg         Result, not yet symmetrised over the field indices.
!! @param bare_term3        Optional, default .false.: term 3 (the exciton-exciton populations and coherences,
!!                          Eq. A9b term 3) reads the BARE current P_nm (vme_ex_inter) instead of
!!                          Pi_nm = i(E_n - E_m) X_nm. Pi has no diagonal and vanishes for degenerate pairs, so
!!                          with it the injection current converges only with the k-mesh (buckled hBN 75x75,
!!                          eta = 0.025 eV: 10% of its converged weight), and in the non-interacting limit it
!!                          leaves a mesh artifact in the shift part (0.90 of the exact value at 75x75,
!!                          eta = 0.1 eV). With P the non-interacting limit equals the single-particle
!!                          rectification (least squares 1.000) and the result is mesh-converged at 75x75.
!! @return void
subroutine get_second_intens_ex_methodA(nfreq, hwp, hwq, hw2, sigma_shg, bare_term3)
  implicit none
  integer,    intent(in)    :: nfreq
  complex(8), intent(in)    :: hwp(nfreq), hwq(nfreq), hw2(nfreq)
  complex(8), intent(inout) :: sigma_shg(3,3,3,nfreq)
  logical,    intent(in), optional :: bare_term3
  logical :: use_bare

  integer, parameter :: nw_chunk = 2000
  integer :: nj, njp, njpp, nn, nnp
  integer :: iw0, iw1, nw_this, ichunk, nchunks
  complex(8), parameter :: ci=(0.0d0,1.0d0), czero=(0.0d0,0.0d0), cone=(1.0d0,0.0d0)
  complex(8), allocatable :: omega_c(:), omega2_c(:), omegap_c(:)
  complex(8), allocatable :: Mmat(:,:), Bmat(:,:), Wmat(:,:), Amat(:,:), term_total(:)
  real(8) :: prefA

  sigma_shg = (0.0d0, 0.0d0)
  prefA = -1.0d0/(dble(npointstotal)*vcell)      ! -C_ee, C_ee = 1 in these units (as in method B)
  use_bare = .false.
  if (present(bare_term3)) use_bare = bare_term3

  allocate(omega_c(nw_chunk), omega2_c(nw_chunk), omegap_c(nw_chunk))
  allocate(Mmat(norb_ex_cut,norb_ex_cut))
  allocate(Bmat(norb_ex_cut,nw_chunk), Wmat(norb_ex_cut,nw_chunk), Amat(norb_ex_cut,nw_chunk))
  allocate(term_total(nw_chunk))

  nchunks = (nfreq + nw_chunk - 1) / nw_chunk
  do ichunk = 1, nchunks
    iw0     = (ichunk-1)*nw_chunk + 1
    iw1     = min(iw0 + nw_chunk - 1, nfreq)
    nw_this = iw1 - iw0 + 1
    omega_c(1:nw_this)  = hwq(iw0:iw1)
    omegap_c(1:nw_this) = hwp(iw0:iw1)
    omega2_c(1:nw_this) = hw2(iw0:iw1)

    ! ---------------- TERM 1:  Pi_n X_nm X*_m / ((hw2 - E_n)(hw_q - E_m)) ----------------
    do njp = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) / (omega_c(1:nw_this) - e_ex(nnp))
        end do
        Mmat = xme_ex_inter(njp,:,:)
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)
        do nj = 1, 3
          do nn = 1, norb_ex_cut
            ! Pi^a_n = -i E_n X^a_n
            Amat(nn,1:nw_this) = (-ci*e_ex(nn)*xme_ex(nj,nn)) / (omega2_c(1:nw_this) - e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          sigma_shg(nj,njp,njpp,iw0:iw1) = sigma_shg(nj,njp,njpp,iw0:iw1) &
              + prefA*term_total(1:nw_this)
        end do
      end do
    end do

    ! ---------------- TERM 2:  Pi*_n X*_nm X_m / ((hw2 + E_n)(hw_q + E_m)) ----------------
    do njp = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = xme_ex(njpp,nnp) / (omega_c(1:nw_this) + e_ex(nnp))
        end do
        Mmat = conjg(xme_ex_inter(njp,:,:))
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)
        do nj = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = conjg(-ci*e_ex(nn)*xme_ex(nj,nn)) / (omega2_c(1:nw_this) + e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          sigma_shg(nj,njp,njpp,iw0:iw1) = sigma_shg(nj,njp,njpp,iw0:iw1) &
              + prefA*term_total(1:nw_this)
        end do
      end do
    end do

    ! ---------------- TERM 3 (minus):  X_n Pi_nm X*_m / ((hw_q + E_n)(hw_p - E_m)) ----------------
    do nj = 1, 3
      do njpp = 1, 3
        do nnp = 1, norb_ex_cut
          Bmat(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) / (omegap_c(1:nw_this) - e_ex(nnp))
        end do
        if (use_bare) then
          ! bare exciton-exciton current P^a_nm (Eq. A10 with o = v); see bare_term3 above
          Mmat = vme_ex_inter(nj,:,:)
        else
          ! Pi^a_nm = i (E_n - E_m) X^a_nm
          do nnp = 1, norb_ex_cut
            do nn = 1, norb_ex_cut
              Mmat(nn,nnp) = ci*(e_ex(nn)-e_ex(nnp))*xme_ex_inter(nj,nn,nnp)
            end do
          end do
        end if
        call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat, norb_ex_cut, &
                    Bmat, norb_ex_cut, czero, Wmat, norb_ex_cut)
        do njp = 1, 3
          do nn = 1, norb_ex_cut
            Amat(nn,1:nw_this) = xme_ex(njp,nn) / (omega_c(1:nw_this) + e_ex(nn))
          end do
          term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat(:,1:nw_this), dim=1)
          ! FIELD PAIRING (Eq. A9b): in term 3, U_n carries the w_q field and U*_m the w_p field, so the X_n
          ! factor (with hw_q + E_n) is the c index and the X*_m factor (with hw_p - E_m) the b index. Stored at
          ! (a, b = njpp, c = njp). Only the b<->c antisymmetric part depends on this: with the indices the other way
          ! round the injection current (circular light, buckled hBN) had the opposite sign to a real-time propagation.
          sigma_shg(nj,njpp,njp,iw0:iw1) = sigma_shg(nj,njpp,njp,iw0:iw1) &
              - prefA*term_total(1:nw_this)
        end do
      end do
    end do
  end do

  deallocate(omega_c, omega2_c, omegap_c, Mmat, Bmat, Wmat, Amat, term_total)
end subroutine get_second_intens_ex_methodA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SHG entry point, kept with its original signature so the long-validated driver and
! tests/test_shg_real_data.f90 are untouched. It is now just the omega_q = omega_p branch of
! get_second_intens_ex_methodB: hw_p = hw_q = omega + i*eta and hw2 = 2*omega + 2i*eta, which is
! the paper's convention: omega_2 is the SUM of the two complex frequencies.
subroutine get_shg_intens_ex_matrix_methodB(wp, eta2, sigma_shg)
  implicit none
  real(8),    intent(in)    :: wp(nw), eta2
  complex(8), intent(inout) :: sigma_shg(3,3,3,nw)
  complex(8), allocatable :: hwp(:), hwq(:), hw2(:)
  integer :: iw
  allocate(hwp(nw), hwq(nw), hw2(nw))
  do iw = 1, nw
    hwp(iw) = cmplx(wp(iw), eta2, 8)
    hwq(iw) = hwp(iw)
    hw2(iw) = hwp(iw) + hwq(iw)
  end do
  call get_second_intens_ex_methodB(nw, hwp, hwq, hw2, sigma_shg)
  deallocate(hwp, hwq, hw2)
end subroutine get_shg_intens_ex_matrix_methodB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENERAL excitonic second-order conductivity sigma^{abc}(w_p+w_q; w_p, w_q), Taghizadeh & Pedersen
! PRB 97, 205432 (2018), Eqs. (B1a)/(B1b). Mirrors get_sigma_general_sp: the physical tensor is the
! intrinsic-permutation average over the PAIRS (alpha,w_p) <-> (beta,w_q),
!     sigma_sym^{a,b,c}(wp,wq) = 1/2 [ sigma^{a,b,c}(wp,wq) + sigma^{a,c,b}(wq,wp) ],
! implemented as two passes with (hwp,hwq) swapped, the second transposed in (b,c). For w_p = w_q the
! two passes coincide and it degenerates to the plain b<->c swap get_sigma_shg_ex already did.
!
! EVERY point, the DC line included, uses the CAUSAL prescription: each frequency carries +i*eta and
! hw_2 = hw_p + hw_q. So Response = rectification is the whole sigma(0; w, -w): the symmetric real part
! (on resonance the shift-current-like response, plus off-resonant terms) and the antisymmetric
! imaginary part, which holds the injection current. It is continuous with the anti-diagonal of a 2D map
! (identical, when the map contains w_2 = 0 points). The SHIFT CURRENT itself, in the convention of
! npj Comput. Mater. 11, 13 (w_q = -(w + i eta), w_2 = 0, Eqs. 9-11), is Response = shift, a separate route
! (get_sigma_shift_ex). The two differ through Eq. (A9b) term 3, the exciton populations and coherences,
! which that convention cancels by construction: on buckled hBN (75x75 and 90x90) the symmetric real part
! of the rectification is 0.53 of the shift current, mesh-converged; the difference comes from pairs of
! distinct bound excitons, while in the continuum and in the non-interacting limit the two agree.
!> Driver for the excitonic two-frequency response sigma^abc(w1+w2; w1, w2): builds the
!! frequency pairs, picks method A or B, symmetrises and writes the result.
!! Method A (term 3 on the bare current) is used wherever w_2 = 0 is in range -- the 1D w_q = -w_p scan
!! (rectification) and 2D maps containing such points; method B elsewhere.
!! @param tag  'electrooptic', 'rectification' or 'general'; names the output file.
!! @return void
subroutine get_sigma_general_ex(tag)
  implicit none
  character(len=*), intent(in) :: tag
  integer :: nfreq, idx, nj, njp, njpp
  real(8) :: eta2
  real(8),    allocatable :: wpg(:), wqg(:)
  complex(8), allocatable :: hwp(:), hwq(:), hw2(:)
  complex(8), allocatable :: sigA(:,:,:,:), sigB(:,:,:,:), sigS(:,:,:,:)
  logical :: same_freq, use_methodA

  if (.not. inter_terms_ready) then
    write(*,*) 'ERROR (sigma_second_ex): xme_ex_inter not populated -- get_ome_ex needs iflag_norder=2.'
    stop 1
  end if

  eta2 = eta/27.211385d0
  call build_freq_pairs(nfreq, wpg, wqg)

  allocate(hwp(nfreq), hwq(nfreq), hw2(nfreq))
  do idx = 1, nfreq
    hwp(idx) = cmplx(wpg(idx), eta2, 8)
    hwq(idx) = cmplx(wqg(idx), eta2, 8)
    hw2(idx) = hwp(idx) + hwq(idx)
  end do

  ! METHOD SELECTION. Where w_2 = 0 is in range (the DC line, w_2 = 2i eta there) method A, Eq. (B1a), with
  ! term 3 on the bare exciton-exciton current (bare_term3 in get_second_intens_ex_methodA): there term 3
  ! carries the injection current, and Pi_nm = i(E_n-E_m) X_nm, which has no diagonal, resolves it only
  ! with the k-mesh. Method B, Eq. (B1b), elsewhere: it is equivalent to method A with Pi and the
  ! long-validated production route of the resonant branches, where term 3 is a small part (0.1% of the
  ! SHG on buckled hBN).
  use_methodA = (.not. two_freq_grid) .and. (abs(freq_ratio + 1.0d0) < 1.0d-8)
  if (two_freq_grid) then
    if (minval(abs(wpg+wqg)) < 1.0d-8) use_methodA = .true.
  end if

  same_freq = (.not. two_freq_grid) .and. (abs(freq_ratio - 1.0d0) < 1.0d-12)
  allocate(sigA(3,3,3,nfreq), sigB(3,3,3,nfreq), sigS(3,3,3,nfreq))
  sigA = (0.0d0,0.0d0); sigB = (0.0d0,0.0d0)

  write(*,*) '    Evaluating second-order conductivity (ex): ', trim(tag)
  if (use_methodA) then
    write(*,'(A,I0,A,I0,A)') '        ', nfreq, ' frequency pairs, ', norb_ex_cut, &
      ' excitons  [method A, Eq. (B1a), term 3 on the bare current: omega_2 = 0 is in range]'
    write(*,'(A)') '        causal prescription (every frequency + i*eta); Re = full DC response, Im = injection channel'
    call get_second_intens_ex_methodA(nfreq, hwp, hwq, hw2, sigA, bare_term3=.true.)
    call get_second_intens_ex_methodA(nfreq, hwq, hwp, hw2, sigB, bare_term3=.true.)
  else
    write(*,'(A,I0,A,I0,A)') '        ', nfreq, ' frequency pairs, ', norb_ex_cut, &
                             ' excitons  [method B, Eq. (B1b)]'
    call get_second_intens_ex_methodB(nfreq, hwp, hwq, hw2, sigA)
    if (same_freq) then
      sigB = sigA
    else
      call get_second_intens_ex_methodB(nfreq, hwq, hwp, hw2, sigB)
    end if
  end if

  ! Both field orderings are evaluated explicitly. On the DC line the swapped pass is NOT the complex
  ! conjugate of the first (they differ by O(1) on buckled hBN), but the symmetrised sum obeys the reality
  ! of the current: its Re is b<->c symmetric and its Im antisymmetric to 1e-15 relative.
  do nj = 1,3
    do njp = 1,3
      do njpp = 1,3
        sigS(nj,njp,njpp,:) = 0.5d0*(sigA(nj,njp,njpp,:) + sigB(nj,njpp,njp,:))
      end do
    end do
  end do

  ! PHYSICAL SIGN, for the reason given in get_sigma_shg_ex (fixed by the real-time simulation).
  sigS = -sigS

  call print_second_general_ex(nfreq, wpg, wqg, sigS, tag)
  deallocate(wpg, wqg, hwp, hwq, hw2, sigA, sigB, sigS)
end subroutine get_sigma_general_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Same layout as the single-particle general writer: hbar*w_p (eV), hbar*w_q (eV), then 54 (Re, Im).
subroutine print_second_general_ex(nfreq, wpg, wqg, sigma, tag)
  implicit none
  integer :: iounit94
  integer,    intent(in) :: nfreq
  real(8),    intent(in) :: wpg(nfreq), wqg(nfreq)
  complex(8), intent(in) :: sigma(3,3,3,nfreq)
  character(len=*), intent(in) :: tag
  integer :: iw, ia, ib, ic
  real(8) :: feps
  feps = sigma2_au_to_si
  open(newunit=iounit94, file='second_ex_'//trim(tag)//'_lengthgauge_'//trim(material_name)//'.dat')
  write(iounit94,'(A)') '# hbar*w_p(eV) hbar*w_q(eV) | excitonic sigma^{abc}(w_p+w_q;w_p,w_q) (Re,Im) uA nm/V^2, abc=xxx,...,zzz'
  do iw = 1, nfreq
    write(iounit94,'(2ES18.10,54ES18.10)') wpg(iw)*27.211385d0, wqg(iw)*27.211385d0, &
      ( ( ( real(feps*sigma(ia,ib,ic,iw)), aimag(feps*sigma(ia,ib,ic,iw)), ic=1,3 ), ib=1,3 ), ia=1,3 )
  end do
  close(iounit94)
end subroutine print_second_general_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_second_ex(nw, wp, sigma_w_ex)
    use omp_lib
    implicit none
    integer :: iounit90

    integer,    intent(in) :: nw
    real(8),    intent(in) :: wp(nw)
    complex(8), intent(in) :: sigma_w_ex(3,3,3,nw)

    integer :: iw
    real(8) :: feps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !feps=-6.623618d-03/(27.21138**2)*1.0d06 !%go from au to (\mu A /V^2)*Angstrongs
    !d=2.6d0 !thickness in angstrongs for MoS2
    !d=3.28d0 !thickness in angstrongs for h-BN
    !feps=feps/(d/0.52917721067121d0) 
    feps = sigma2_au_to_si !%go from au to (\mu A /V^2)*nm
   
    open(newunit=iounit90, file='shift_ex_lengthgauge_'//trim(material_name)//'.dat')

    ! serial on purpose (audit 2026-09-24): this loop is pure file I/O and every iteration was
    ! inside !$omp ordered, which serialises it completely -- the parallel wrapper only added
    ! thread spawn and synchronisation cost..
    do iw = 1, nw
      write(iounit90,*) wp(iw)*27.211385d0, &
        realpart(feps*sigma_w_ex(1,1,1,iw)), realpart(feps*sigma_w_ex(1,1,2,iw)), &
        realpart(feps*sigma_w_ex(1,1,3,iw)), realpart(feps*sigma_w_ex(1,2,1,iw)), &
        realpart(feps*sigma_w_ex(1,2,2,iw)), realpart(feps*sigma_w_ex(1,2,3,iw)), &
        realpart(feps*sigma_w_ex(1,3,1,iw)), realpart(feps*sigma_w_ex(1,3,2,iw)), &
        realpart(feps*sigma_w_ex(1,3,3,iw)), realpart(feps*sigma_w_ex(2,1,1,iw)), &
        realpart(feps*sigma_w_ex(2,1,2,iw)), realpart(feps*sigma_w_ex(2,1,3,iw)), &
        realpart(feps*sigma_w_ex(2,2,1,iw)), realpart(feps*sigma_w_ex(2,2,2,iw)), &
        realpart(feps*sigma_w_ex(2,2,3,iw)), realpart(feps*sigma_w_ex(2,3,1,iw)), &
        realpart(feps*sigma_w_ex(2,3,2,iw)), realpart(feps*sigma_w_ex(2,3,3,iw)), &
        realpart(feps*sigma_w_ex(3,1,1,iw)), realpart(feps*sigma_w_ex(3,1,2,iw)), &
        realpart(feps*sigma_w_ex(3,1,3,iw)), realpart(feps*sigma_w_ex(3,2,1,iw)), &
        realpart(feps*sigma_w_ex(3,2,2,iw)), realpart(feps*sigma_w_ex(3,2,3,iw)), &
        realpart(feps*sigma_w_ex(3,3,1,iw)), realpart(feps*sigma_w_ex(3,3,2,iw)), &
        realpart(feps*sigma_w_ex(3,3,3,iw))
    end do

    close(iounit90)

  end subroutine print_sigma_second_ex

!!!!
  ! Writes the excitonic SHG conductivity sigma^{abc}(2*omega; omega, omega) to
  ! shg_ex_lengthgauge_<material>.dat.
  ! Columns: hbar*omega (eV) -- the FUNDAMENTAL (driving) photon energy, NOT 2*hbar*omega --
  ! then for a=x,y,z; b=x,y,z; c=x,y,z (c fastest): Re, Im.
  ! AXIS CONVENTION CHANGED 2026-09-24, same change and same reasons as
  ! print_shg_second_sp: a two-photon resonance of an exciton at energy E_N now appears at
  ! hbar*omega = E_N/2, a one-photon resonance at hbar*omega = E_N. Files written before this
  ! date carry the old 2*hbar*omega axis.
  ! Units: uA*nm/V^2, same a.u. -> SI factor as the shift conductivity (both are 2D second-order
  ! conductivities, j = sigma E E, and the kernel is in Hartree a.u.: X in Bohr, vcell in Bohr^2).
  ! No spin degeneracy factor g is included, consistent with the rest of the code
  ! (Taghizadeh & Pedersen's C_ee contains g = 2). The overall sign convention of e is not
  ! examined here.
  subroutine print_shg_second_ex(nw, wp, sigma_shg)
    implicit none
    integer :: iounit91
    integer,    intent(in) :: nw
    real(8),    intent(in) :: wp(nw)
    complex(8), intent(in) :: sigma_shg(3,3,3,nw)
    integer :: iw, ia, ib, ic
    real(8) :: feps

    feps = sigma2_au_to_si

    open(newunit=iounit91, file='shg_ex_lengthgauge_'//trim(material_name)//'.dat')
    write(iounit91,'(A)') '# hbar*omega(eV) [FUNDAMENTAL, not 2*hbar*omega] | sigma^{abc}(2w;w,w) (Re, Im) in uA nm/V^2, abc = xxx,xxy,xxz,xyx,...,zzz'

    do iw = 1, nw
      write(iounit91,'(ES18.10,54ES18.10)') wp(iw)*27.211385d0, &
        ( ( ( real(feps*sigma_shg(ia,ib,ic,iw)), aimag(feps*sigma_shg(ia,ib,ic,iw)), &
              ic=1,3 ), ib=1,3 ), ia=1,3 )
    end do

    close(iounit91)

  end subroutine print_shg_second_ex
!!!!

end module sigma_second_ex

