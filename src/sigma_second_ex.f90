module sigma_second_ex
  use constants_math
  use parser_input_file, &
    only:e1,e2,eta,nw,response_text,broadening_type_text
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
      !write(*,*) 'The optical response',response_text,'has been evaluated'
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
    
    !call the subroutine to compute the shift conductivity
    if (.not. inter_terms_ready) then
      write(*,*) 'ERROR (sigma_second_ex): xme_ex_inter/vme_ex_inter not '// &
                'populated — get_ome_ex must be called with iflag_norder=2 first.'
      stop 1
    end if
!     call get_shift_intens_ex(wp,eta2,sigma_w_ex)
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_second_ex(nw, wp, sigma_w_ex)
    use omp_lib
    implicit none

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
    feps=(6.623618d-03)*(1.0d+06)*(27.211386**(-2))*(5.291772d-11)*(1.0d+09) !%go from au to (\mu A /V^2)*nm
   
    open(90, file='shift_ex_lengthgauge_'//trim(material_name)//'.dat')

    ! serial on purpose (audit 2026-09-24): this loop is pure file I/O and every iteration was
    ! inside !$omp ordered, which serialises it completely -- the parallel wrapper only added
    ! thread spawn and synchronisation cost. HANDOFF 8.35.
    do iw = 1, nw
      write(90,*) wp(iw)*27.211385d0, &
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

    close(90)

  end subroutine print_sigma_second_ex

!!!!
  ! Writes the excitonic SHG conductivity sigma^{abc}(2*omega; omega, omega) to
  ! shg_ex_lengthgauge_<material>.dat.
  ! Columns: hbar*omega (eV) -- the FUNDAMENTAL (driving) photon energy, NOT 2*hbar*omega --
  ! then for a=x,y,z; b=x,y,z; c=x,y,z (c fastest): Re, Im.
  ! AXIS CONVENTION CHANGED 2026-09-24 (HANDOFF 8.33), same change and same reasons as
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
    integer,    intent(in) :: nw
    real(8),    intent(in) :: wp(nw)
    complex(8), intent(in) :: sigma_shg(3,3,3,nw)
    integer :: iw, ia, ib, ic
    real(8) :: feps

    feps = (6.623618d-03)*(1.0d+06)*(27.211386d0**(-2))*(5.291772d-11)*(1.0d+09)

    open(91, file='shg_ex_lengthgauge_'//trim(material_name)//'.dat')
    write(91,'(A)') '# hbar*omega(eV) [FUNDAMENTAL, not 2*hbar*omega] | sigma^{abc}(2w;w,w) (Re, Im) in uA nm/V^2, abc = xxx,xxy,xxz,xyx,...,zzz'

    do iw = 1, nw
      write(91,'(ES18.10,54ES18.10)') wp(iw)*27.211385d0, &
        ( ( ( real(feps*sigma_shg(ia,ib,ic,iw)), aimag(feps*sigma_shg(ia,ib,ic,iw)), &
              ic=1,3 ), ib=1,3 ), ia=1,3 )
    end do

    close(91)

  end subroutine print_shg_second_ex
!!!!

end module sigma_second_ex

