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
  subroutine get_sigma_shift_ex()
    implicit none
    
	!here
	integer :: iw
    dimension :: wp(nw)
    dimension :: sigma_w_ex(3,3,3,nw)

    real*8 :: wp,eta2
    complex*16 :: sigma_w_ex
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
                - ( s1*(-cmplx(0.0d0,1.0d0,8)*pi*d1_arr(:)) &
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
  allocate(term_total(nw_chunk))

  if (mode == 1) then
    allocate(gauss1(norb_ex_cut,nw_chunk), gauss2(norb_ex_cut,nw_chunk))
    allocate(gauss3(norb_ex_cut,nw_chunk), gauss4(norb_ex_cut,nw_chunk))
  else
    allocate(Afac1(norb_ex_cut), Afac2(norb_ex_cut))
    allocate(Bfac1(norb_ex_cut,nw_chunk), Bfac2(norb_ex_cut,nw_chunk))
    allocate(Cfac3(norb_ex_cut,nw_chunk), Dfac3(norb_ex_cut,nw_chunk))
    Afac1(:) = 1.0d0/(-e_ex(:)+cmplx(0.0d0,eta2,8))
    Afac2(:) = 1.0d0/( e_ex(:)+cmplx(0.0d0,eta2,8))
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
        gauss3(nn,1:nw_this) = gauss2(nn,1:nw_this)
        gauss4(nn,1:nw_this) = 1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)* &
            exp(-0.5d0/(eta2**2)*( wp(iw0:iw1)-e_ex(nn))**2)
      end do
    else
      do nn = 1, norb_ex_cut
        Bfac1(nn,1:nw_this) = 1.0d0/(-wp(iw0:iw1)-e_ex(nn)+cmplx(0.0d0,eta2,8))
        Bfac2(nn,1:nw_this) = 1.0d0/(-wp(iw0:iw1)+e_ex(nn)+cmplx(0.0d0,eta2,8))
        Dfac3(nn,1:nw_this) = 1.0d0/( wp(iw0:iw1)-e_ex(nn)+cmplx(0.0d0,eta2,8))
        Cfac3(nn,1:nw_this) = 1.0d0/( e_ex(nn)-wp(iw0:iw1)+cmplx(0.0d0,eta2,8))
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
            Avec1(:) = -vme_ex(nj,:) / e_ex(:)
            Avec2(:) = conjg(vme_ex(nj,:)) / e_ex(:)
            term_total(1:nw_this) = matmul(Avec1, Wmat1(:,1:nw_this)) &
                                   + matmul(Avec2, Wmat2(:,1:nw_this))
            sigma_w_ex(nj,njp,njpp,iw0:iw1) = sigma_w_ex(nj,njp,njpp,iw0:iw1) &
                - term_total(1:nw_this) / (dble(npointstotal)*vcell)
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
            Avec1(:) = vme_ex(nj,:)        * Afac1(:)   ! term1, A
            Avec2(:) = conjg(vme_ex(nj,:)) * Afac1(:)   ! term1, A*
            Avec3(:) = conjg(vme_ex(nj,:)) * Afac2(:)   ! term2, A
            Avec4(:) = vme_ex(nj,:)        * Afac2(:)   ! term2, A*

            term_total(1:nw_this) = &
                ( matmul(Avec1, Wmat1(:,1:nw_this)) - matmul(Avec2, Wmat2(:,1:nw_this)) &
                + matmul(Avec3, Wmat3(:,1:nw_this)) - matmul(Avec4, Wmat4(:,1:nw_this)) ) / (2.0d0*ci)

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
          Mmat1 = vme_ex_inter(nj,:,:)
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat1, norb_ex_cut, &
                      Bmat1, norb_ex_cut, czero, Wmat1, norb_ex_cut)

          do njp = 1, 3
            do nn = 1, norb_ex_cut
              Amat(nn,1:nw_this) = -xme_ex(njp,nn) * (-pi**2) * gauss3(nn,1:nw_this)
            end do
            term_total(1:nw_this) = sum(Amat(:,1:nw_this)*Wmat1(:,1:nw_this), dim=1)
            sigma_w_ex(nj,njp,njpp,iw0:iw1) = sigma_w_ex(nj,njp,njpp,iw0:iw1) &
                - term_total(1:nw_this) / (dble(npointstotal)*vcell)
          end do

        else   ! lorentzian
          do nnp = 1, norb_ex_cut
            Bmat1(nnp,1:nw_this) = conjg(xme_ex(njpp,nnp)) * Dfac3(nnp,1:nw_this)
            Bmat2(nnp,1:nw_this) = xme_ex(njpp,nnp)        * Dfac3(nnp,1:nw_this)
          end do
          Mmat1 = vme_ex_inter(nj,:,:)
          Mmat2 = conjg(Mmat1)

          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat1, norb_ex_cut, &
                      Bmat1, norb_ex_cut, czero, Wmat1, norb_ex_cut)
          call zgemm('N','N', norb_ex_cut, nw_this, norb_ex_cut, cone, Mmat2, norb_ex_cut, &
                      Bmat2, norb_ex_cut, czero, Wmat2, norb_ex_cut)

          do njp = 1, 3
            do nn = 1, norb_ex_cut
              Amat(nn,1:nw_this)  = -xme_ex(njp,nn)        * Cfac3(nn,1:nw_this)
              Amat2(nn,1:nw_this) = -conjg(xme_ex(njp,nn)) * Cfac3(nn,1:nw_this)
            end do
            term_total(1:nw_this) = &
                ( sum(Amat(:,1:nw_this)*Wmat1(:,1:nw_this), dim=1) &
                - sum(Amat2(:,1:nw_this)*Wmat2(:,1:nw_this), dim=1) ) / (2.0d0*ci)
            sigma_w_ex(nj,njp,njpp,iw0:iw1) = sigma_w_ex(nj,njp,njpp,iw0:iw1) &
                - term_total(1:nw_this) / (dble(npointstotal)*vcell)
          end do
        end if

      end do
    end do

  end do   ! ichunk

  if (mode == 1) then
    deallocate(gauss1, gauss2, gauss3, gauss4)
  else
    deallocate(Afac1, Afac2, Bfac1, Bfac2, Cfac3, Dfac3)
  end if
  deallocate(Mmat1, Mmat2, Bmat1, Bmat2, Bmat3, Bmat4)
  deallocate(Wmat1, Wmat2, Wmat3, Wmat4)
  deallocate(Avec1, Avec2, Avec3, Avec4, Amat, Amat2, term_total)

end subroutine get_shift_intens_ex_matrix
  
  !!!!
  
  ! Frequency-independent part: s1, s2, s3 only.
  subroutine get_shift_kernel_ex_static(mode, nj, njp, njpp, nn, nnp, s1, s2, s3)
    implicit none
    integer,    intent(in)  :: mode, nj, njp, njpp, nn, nnp
    complex(8), intent(out) :: s1, s2, s3
    integer :: nj1, nj2, nj3

    nj1 = nj; nj2 = njp; nj3 = njpp

    if (mode == 1) then   ! gaussian
      s1 = -vme_ex(nj1,nn)/e_ex(nn) * xme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))         ! a→V_0N, b→R_NN' — matches term1
      s2 =  conjg(vme_ex(nj1,nn))/e_ex(nn) * conjg(xme_ex_inter(nj2,nn,nnp)) * xme_ex(nj3,nnp)  ! matches term2
!       s3 = -xme_ex(nj1,nn) * vme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))                ! a→R_0N, b→V_NN'  ← WRONG
      s3 = -xme_ex(nj2,nn) * vme_ex_inter(nj1,nn,nnp) * conjg(xme_ex(nj3,nnp))
      
    else                    ! lorentzian
      s1 = vme_ex(nj1,nn) * xme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))         ! a→V_0N, b→R_NN' — matches term1
      s2 = conjg(vme_ex(nj1,nn)) * conjg(xme_ex_inter(nj2,nn,nnp)) * xme_ex(nj3,nnp)  ! matches term2
!       s3 = -xme_ex(nj1,nn) * vme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))      ! a→R_0N, b→V_NN'  ← WRONG
      s3 = -xme_ex(nj2,nn) * vme_ex_inter(nj1,nn,nnp) * conjg(xme_ex(nj3,nnp))
      ! PATCH: the lorentzian mode's real-part transform (cmplx(aimag(s),0))
      ! is also frequency-independent, so it's folded in here rather than
      ! redone inside the iw loop.
      s1 = cmplx(aimag(s1), 0.0d0, 8)
      s2 = cmplx(aimag(s2), 0.0d0, 8)
      s3 = cmplx(aimag(s3), 0.0d0, 8)
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
    d1(:) = 1.0d0 / ( (-e_ex(nn)+cmplx(0.0d0,eta2,8)) * (-wp(:)-e_ex(nnp)+cmplx(0.0d0,eta2,8)) )
    d2(:) = 1.0d0 / ( ( e_ex(nn)+cmplx(0.0d0,eta2,8)) * (-wp(:)+e_ex(nnp)+cmplx(0.0d0,eta2,8)) )
    d3(:) = 1.0d0 / ( ( e_ex(nn)-wp(:)+cmplx(0.0d0,eta2,8)) * ( wp(:)-e_ex(nnp)+cmplx(0.0d0,eta2,8)) )
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

    aux1 = s1 * (-complex(0.0d0,1.0d0)*pi*d1)
    aux2 = s2 * (-complex(0.0d0,1.0d0)*pi*d2)
    aux3 = s3 * (-pi**2*d3*d4)
    shift_kernel_ex = -(aux1+aux2+aux3)
  else                    ! lorentzian
    d1 = 1.0d0 / ((omega2-e_ex(nn) +cmplx(0.0d0,eta2,8)) * (omegaq-e_ex(nnp)+cmplx(0.0d0,eta2,8)))
    d2 = 1.0d0 / ((omega2+e_ex(nn) +cmplx(0.0d0,eta2,8)) * (omegaq+e_ex(nnp)+cmplx(0.0d0,eta2,8)))
    ! d3/aux3 intentionally not computed: same as the original, aux3 was
    ! excluded from the final sum (`shift_kernel_ex=-(aux1+aux2)`).
    ! now computed
    d3 = 1.0d0 / ((e_ex(nn)-omegap +cmplx(0.0d0,eta2,8)) * (omegap-e_ex(nnp)+cmplx(0.0d0,eta2,8)))
    
    aux1 = s1 * d1
    aux2 = s2 * d2
    aux3 = s3 * d3                                                                                    ! <-- added
    shift_kernel_ex = -(aux1+aux2+aux3)
  end if

end subroutine get_shift_kernel_ex_freq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine get_shift_kernel_ex(eta2, nj, njp, njpp, nn, nnp, &
!                                   omegap, omegaq, omega2, shift_kernel_ex)
!     implicit none
! 
!     integer,    intent(in)  :: nj, njp, njpp, nn, nnp
!     real(8),    intent(in)  :: eta2, omegap, omegaq, omega2
!     complex(8), intent(out) :: shift_kernel_ex
! 
!     integer    :: isym, ilorentzian, ihuang, imine
!     integer    :: nj1, nj2, nj3
!     complex(8) :: omegaq_c, omegap_q, omega2_c
!     complex(8) :: s1, s2, s3, eta2p
!     complex(8) :: aux1, aux2, aux3, aux4, aux5, aux6
!     complex(8) :: d1, d2, d3, d4, d5, d6, d7, d8
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!     isym = 0
!     if (isym .eq. 0) then
!       nj1 = nj
!       nj2 = njp
!       nj3 = njpp
! 
!       if (trim(broadening_type_text) == 'gaussian') then
!         imine       = 1
!         ilorentzian = 0
!       else if (trim(broadening_type_text) == 'lorentzian') then
!         imine       = 0
!         ilorentzian = 1
!       else
!         imine       = 1
!         ilorentzian = 0
!       end if
! 
!       ihuang = 0
! 
!       if (imine .eq. 1) then
! 
!         d1=1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta2**2)*(omegaq-e_ex(nnp))**2)
!         d2=1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta2**2)*(omegaq+e_ex(nnp))**2)
!         d3=1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta2**2)*(omegaq+e_ex(nn))**2)
!         d4=1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta2**2)*(omegap-e_ex(nnp))**2)
! 
!         !First version: 2023. WORKING
!         s1=-vme_ex(nj1,nn)/e_ex(nn)*xme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp))   
!         s2=conjg(vme_ex(nj1,nn))/e_ex(nn)*conjg(xme_ex_inter(nj2,nn,nnp))*xme_ex(nj3,nnp)
!         s3=-xme_ex(nj1,nn)*vme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp)) 
! 
!         aux1=s1*(-complex(0.0d0,1.0d0)*pi*d1)
!         aux2=s2*(-complex(0.0d0,1.0d0)*pi*d2)
!         aux3=s3*(-pi**2*d3*d4)
! 
!         shift_kernel_ex=-(aux1+aux2+aux3)
!         !s_kernel=-(aux1+aux2) !+aux3)
!         !s_kernel=-aux2 !+aux3)
! 
!         !full R. WORKING
!         !s1=xme_ex(nj1,nn)*xme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp))
!         !s2=conjg(xme_ex(nj1,nn))*conjg(xme_ex_inter(nj2,nn,nnp))*xme_ex(nj3,nnp)
!         !s3=-xme_ex(nj1,nn)*vme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp))
!         !s3=-xme_ex(nj1,nn)*vme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp))
!         !aux1=s1*(pi*d1)
!         !aux2=s2*(pi*d2)
!         !aux3=s3*(-pi**2*d3*d4)
!         !s_kernel=-(aux1+aux2)
! 
!         !full V
!         !s1=(complex(0.0d0,1.0d0)*vme_ex(nj1,nn)/e_ex(nn))*xme_ex_inter(nj2,nn,nnp) &
!         !*(-complex(0.0d0,1.0d0)*conjg(vme_ex(nj3,nn))/e_ex(nn))
!         !s2=(-complex(0.0d0,1.0d0)*conjg(vme_ex(nj1,nn))/e_ex(nn))*conjg(xme_ex_inter(nj2,nn,nnp)) &
!         !*(complex(0.0d0,1.0d0)*vme_ex(nj3,nnp)/e_ex(nnp))
!         !s3=-xme_ex(nj1,nn)*vme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp))
!         !aux1=s1*(pi*d1)
!         !aux2=s2*(pi*d2)
!         !aux3=s3*(-pi**2*d3*d4)
!         !shift_kernel_ex=-(aux2)
!             
!         !s1=-vme_ex(nj1,nn)/e_ex(nn)*vme_ex_inter(nj2,nn,nnp)*conjg(vme_ex(nj3,nnp))
!         !s2=conjg(vme_ex(nj1,nn))/e_ex(nn)*conjg(vme_ex_inter(nj2,nn,nnp))*vme_ex(nj3,nnp)
!         !s3=-vme_ex(nj1,nn)*vme_ex_inter(nj2,nn,nnp)*conjg(vme_ex(nj3,nnp))
! 
!         !s1=-1.0d0/omegap**2*vme_ex(nj1,nn)*vme_ex_inter(nj2,nn,nnp)*conjg(vme_ex(nj3,nnp))   
!         !s2=-1.0d0/omegap**2*conjg(vme_ex(nj1,nn))*conjg(vme_ex_inter(nj2,nn,nnp))*vme_ex(nj3,nnp)
!         !s1=-vme_ex(nj1,nn)/e_ex(nn)*xme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp))
!         !s2=conjg(vme_ex(nj1,nn))/e_ex(nn)*conjg(vme_ex_inter(nj2,nn,nnp))*xme_ex(nj3,nnp)
!         !s3=+1.0d0/omegap**2*vme_ex(nj1,nn)*vme_ex_inter(nj2,nn,nnp)*conjg(vme_ex(nj3,nnp))
!       end if
!       
!       !Taghizadeh 2018 with broadening (10.1103/PhysRevB.97.205432)
!       !tag and pedersen B1a
!       if (ilorentzian .eq. 1) then
!         d1 = 1.0d0 / ((omega2-e_ex(nn) +cmplx(0.0d0,eta2,8)) * (omegaq-e_ex(nnp)+cmplx(0.0d0,eta2,8)))
!         d2 = 1.0d0 / ((omega2+e_ex(nn) +cmplx(0.0d0,eta2,8)) * (omegaq+e_ex(nnp)+cmplx(0.0d0,eta2,8)))
!         !d3 = 1.0d0 / ((omegaq+e_ex(nn) +cmplx(0.0d0,eta2,8)) * (omegap-e_ex(nnp)+cmplx(0.0d0,eta2,8)))
! 
!         s1 = vme_ex(nj1,nn)        * xme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))
!         s2 = conjg(vme_ex(nj1,nn)) * conjg(xme_ex_inter(nj2,nn,nnp)) * xme_ex(nj3,nnp)
!         !s3 = -xme_ex(nj1,nn)       * vme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))
! 
!         s1 = cmplx(aimag(s1), 0.0d0, 8)
!         s2 = cmplx(aimag(s2), 0.0d0, 8)
!         !s3 = cmplx(aimag(s3), 0.0d0, 8)
! 
!         aux1 = s1 * d1
!         aux2 = s2 * d2
!         !aux3 = s3 * d3
!         
!         !s_kernel=-(aux1+aux2+aux3)
!         shift_kernel_ex=-(aux1+aux2) !+aux3)
!         
!         !s_kernel=-aimag((aux1+aux2))
!         !Louie .eq. tag and pedersen B1b
!         !d1=1.0d0/((omega2-e_ex(nn)+complex(0.0d0,eta2))*(omegaq-e_ex(nnp)+complex(0.0d0,eta2)))
!         !d2=1.0d0/((omega2+e_ex(nn)+complex(0.0d0,eta2))*(omegaq+e_ex(nnp)+complex(0.0d0,eta2)))
!         !d3=1.0d0/((omegaq+e_ex(nn)+complex(0.0d0,eta2))*(omegap-e_ex(nnp)+complex(0.0d0,eta2)))
!         !s1=xme_ex(nj1,nn)*xme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp))   
!         !s2=conjg(xme_ex(nj1,nn))*conjg(xme_ex_inter(nj2,nn,nnp))*xme_ex(nj3,nnp)
!         !s3=-xme_ex(nj1,nn)*xme_ex_inter(nj2,nn,nnp)*conjg(xme_ex(nj3,nnp)) 
!         !aux1=s1*d1		  
!         !aux2=s2*d2
!         !aux3=s3*d3
!         !s_kernel=aux1+aux2+aux3
!       end if
! 
!       if (ihuang .eq. 1) then
!         d1 = 1.0d0 / (omegap + e_ex(nnp) + eta2)
!         d2 = 1.0d0 / (omegap - e_ex(nnp) + eta2)
!         d3 = 1.0d0 / (omegap - e_ex(nn)  + eta2)
!         d4 = 1.0d0 / (omegap - e_ex(nnp) + eta2)
!         d5 = 1.0d0 / (-omegap + e_ex(nnp) - eta2)
!         d6 = 1.0d0 / (-omegap - e_ex(nnp) + eta2)
!         d7 = 1.0d0 / (-omegap - e_ex(nn)  + eta2)
!         d8 = 1.0d0 / (-omegap - e_ex(nnp) + eta2)
! 
!         aux1 =  cmplx(0.0d0,1.0d0,8)*xme_ex(nj1,nn)*conjg(xme_ex_inter(nj2,nn,nnp))*conjg(xme_ex(nj3,nnp))*d1
!         aux2 =  cmplx(0.0d0,1.0d0,8)*conjg(xme_ex(nj1,nn))*xme_ex_inter(nj2,nn,nnp)*xme_ex(nj3,nnp)*d2
!         aux3 = -xme_ex(nj2,nn)*vme_ex_inter(nj1,nn,nnp)*conjg(xme_ex(nj3,nnp))*d3*d4
! 
!         !b<--->c and \omega<---->-\omega
!         aux4 =  cmplx(0.0d0,1.0d0,8)*xme_ex(nj1,nn)*conjg(xme_ex_inter(nj3,nn,nnp))*conjg(xme_ex(nj2,nnp))*d5
!         aux5 =  cmplx(0.0d0,1.0d0,8)*conjg(xme_ex(nj1,nn))*xme_ex_inter(nj3,nn,nnp)*xme_ex(nj2,nnp)*d6
!         aux6 = -xme_ex(nj3,nn)*vme_ex_inter(nj3,nn,nnp)*conjg(xme_ex(nj2,nnp))*d7*d8
! 
!         !s_kernel=(aux1+aux2+aux3+aux4+aux5+aux6)	
!         shift_kernel_ex = -(aux1 + aux2 + aux3 + aux4 + aux5 + aux6)
!       end if
! 
!     end if
! 
!   end subroutine get_shift_kernel_ex

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

    !$omp parallel do schedule(static) ordered private(iw)
    do iw = 1, nw
      !$omp ordered
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
      !$omp end ordered
    end do
    !$omp end parallel do

    close(90)

  end subroutine print_sigma_second_ex



end module sigma_second_ex

