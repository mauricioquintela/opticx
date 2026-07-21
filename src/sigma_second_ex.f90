module sigma_second_ex
  use constants_math
	use parser_input_file, &
	only:e1,e2,eta,nw,response_text,broadening_type_text
  use parser_wannier90_tb, &
  only:material_name
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex, &
  e_ex,fk_ex
  use ome_ex, &
  only:e_ex,xme_ex,vme_ex, &
  xme_ex_inter,vme_ex_inter
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
    call get_shift_intens_ex(wp,eta2,sigma_w_ex)
	!print shift conductivity (ex)
	call print_sigma_second_ex(nw,wp,sigma_w_ex)
    write(*,*) '    Shift conductivity (ex) has been printed'
  end subroutine get_sigma_shift_ex


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! subroutine get_shift_intens_ex(wp, eta2, sigma_w_ex)
!     use omp_lib
!     implicit none
! 
!     real(8),    intent(in)    :: wp(nw), eta2
!     complex(8), intent(inout) :: sigma_w_ex(3,3,3,nw)
! 
!     integer     :: iw, nn, nnp, nj, njp, njpp
!     real(8)     :: omegap, omegaq, omega2
!     complex(8)  :: shift_kernel_ex1, shift_kernel_ex
!     ! sigma_w_ex_intra is per-iw — declare as a local scalar accumulator
!     ! PATCH: computation of this term is commented out below (dead output —
!     ! never returned, printed, or combined with sigma_w_ex). Left declared
!     ! and zeroed so it's a one-line change to re-enable if the intraband
!     ! ("nn==nnp") contribution turns out to be needed later.
!     !complex(8)  :: sigma_w_ex_intra(3,3,3,nw)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!     sigma_w_ex       = (0.0d0, 0.0d0)
!     !sigma_w_ex_intra = (0.0d0, 0.0d0)
! 
!     ! Each iw is fully independent: omegap/omegaq/omega2 depend only on iw,
!     ! and sigma_w_ex(nj,njp,njpp,iw) is written only by thread owning iw.
!     !$omp parallel do schedule(dynamic) &
!     !$omp   private(iw, nn, nnp, nj, njp, njpp, &
!     !$omp           omegap, omegaq, omega2, &
!     !$omp           shift_kernel_ex1, shift_kernel_ex)
!     do iw = 1, nw
!       omegap = wp(iw)
!       omegaq = -wp(iw)
!       omega2 = 0.0d0
! 
!       do nn = 1, norb_ex_cut
!         do nj = 1, 3
!           do njp = 1, 3
!             do njpp = 1, 3
! 
!               ! PATCH: intraband (nn==nnp) term commented out — it was
!               ! computed at real cost (one extra get_shift_kernel_ex call
!               ! per (nn,nj,njp,njpp,iw), i.e. an additional
!               ! O(norb_ex_cut*27*nw) kernel evaluations) but never used
!               ! downstream: not returned from this subroutine, not printed,
!               ! not combined with sigma_w_ex. Commenting it out removes
!               ! that dead work.
!               !
!               ! call get_shift_kernel_ex(eta2, nj, njp, njpp, nn, nn, &
!               !                          omegap, omegaq, omega2, shift_kernel_ex)
!               ! sigma_w_ex_intra(nj,njp,njpp,iw) = sigma_w_ex_intra(nj,njp,njpp,iw) &
!               !   + 1.0d0 / (dble(npointstotal) * vcell) * shift_kernel_ex
! 
!               do nnp = 1, norb_ex_cut
!                 call get_shift_kernel_ex(eta2, nj, njp, njpp, nn, nnp, &
!                                          omegap, omegaq, omega2, shift_kernel_ex1)
!                 sigma_w_ex(nj,njp,njpp,iw) = sigma_w_ex(nj,njp,njpp,iw) &
!                   + 1.0d0 / (dble(npointstotal) * vcell) * shift_kernel_ex1
!               end do
! 
!             end do
!           end do
!         end do
!       end do
!     end do
!     !$omp end parallel do
! 
!   end subroutine get_shift_intens_ex

  subroutine get_shift_intens_ex(wp, eta2, sigma_w_ex)
    use omp_lib
    implicit none

    real(8),    intent(in)    :: wp(nw), eta2
    complex(8), intent(inout) :: sigma_w_ex(3,3,3,nw)

    integer     :: iw, nn, nnp, nj, njp, njpp
    integer     :: mode   ! 1 = gaussian, 2 = lorentzian
    real(8)     :: omegap, omegaq, omega2
    complex(8)  :: s1, s2, s3
    complex(8)  :: shift_kernel_ex1
    complex(8), allocatable :: sigma_w_ex_t(:,:,:,:)

    ! PATCH: broadening-mode string comparison hoisted OUT of the hot path.
    ! Previously re-evaluated (trim + string compare) on every one of the
    ! nw*norb_ex_cut^2*27 calls to get_shift_kernel_ex; now done exactly once.
    if (trim(broadening_type_text) == 'lorentzian') then
      mode = 2
    else
      mode = 1   ! gaussian, also the default fallback
    end if

    sigma_w_ex = (0.0d0, 0.0d0)

    !$omp parallel default(none) &
    !$omp   shared(mode, eta2, wp, nw, norb_ex_cut, npointstotal, vcell, sigma_w_ex) &
    !$omp   private(nn, nj, njp, njpp, nnp, iw, omegap, omegaq, omega2, &
    !$omp           s1, s2, s3, shift_kernel_ex1, sigma_w_ex_t)

    allocate(sigma_w_ex_t(3,3,3,nw))
    sigma_w_ex_t = (0.0d0, 0.0d0)

    !$omp do schedule(dynamic)
    do nn = 1, norb_ex_cut
      do nj = 1, 3
        do njp = 1, 3
          do njpp = 1, 3
            do nnp = 1, norb_ex_cut

              ! PATCH: frequency-independent physics computed ONCE per
              ! (nj,njp,njpp,nn,nnp) instead of nw times.
              call get_shift_kernel_ex_static(mode, nj, njp, njpp, nn, nnp, s1, s2, s3)

              do iw = 1, nw
                omegap = wp(iw)
                omegaq = -wp(iw)
                omega2 = 0.0d0

                call get_shift_kernel_ex_freq(mode, eta2, omegap, omegaq, omega2, &
                                              nn, nnp, s1, s2, s3, shift_kernel_ex1)

                sigma_w_ex_t(nj,njp,njpp,iw) = sigma_w_ex_t(nj,njp,njpp,iw) &
                  + 1.0d0 / (dble(npointstotal) * vcell) * shift_kernel_ex1
              end do

            end do
          end do
        end do
      end do
    end do
    !$omp end do

    !$omp critical
      sigma_w_ex = sigma_w_ex + sigma_w_ex_t
    !$omp end critical

    deallocate(sigma_w_ex_t)
    !$omp end parallel

  end subroutine get_shift_intens_ex
  
  ! Frequency-independent part: s1, s2, s3 only.
  subroutine get_shift_kernel_ex_static(mode, nj, njp, njpp, nn, nnp, s1, s2, s3)
    implicit none
    integer,    intent(in)  :: mode, nj, njp, njpp, nn, nnp
    complex(8), intent(out) :: s1, s2, s3
    integer :: nj1, nj2, nj3

    nj1 = nj; nj2 = njp; nj3 = njpp

    if (mode == 1) then   ! gaussian
      s1 = -vme_ex(nj1,nn)/e_ex(nn) * xme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))
      s2 =  conjg(vme_ex(nj1,nn))/e_ex(nn) * conjg(xme_ex_inter(nj2,nn,nnp)) * xme_ex(nj3,nnp)
      s3 = -xme_ex(nj1,nn) * vme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))
    else                    ! lorentzian
      s1 = vme_ex(nj1,nn) * xme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))
      s2 = conjg(vme_ex(nj1,nn)) * conjg(xme_ex_inter(nj2,nn,nnp)) * xme_ex(nj3,nnp)
      s3 = -xme_ex(nj1,nn) * vme_ex_inter(nj2,nn,nnp) * conjg(xme_ex(nj3,nnp))
      ! PATCH: the lorentzian mode's real-part transform (cmplx(aimag(s),0))
      ! is also frequency-independent, so it's folded in here rather than
      ! redone inside the iw loop.
      s1 = cmplx(aimag(s1), 0.0d0, 8)
      s2 = cmplx(aimag(s2), 0.0d0, 8)
      s3 = cmplx(aimag(s3), 0.0d0, 8)
    end if
  end subroutine get_shift_kernel_ex_static
  
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
    aux1 = s1 * d1
    aux2 = s2 * d2
    shift_kernel_ex = -(aux1+aux2)
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

