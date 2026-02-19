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
    call read_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,berry_eigen_ex_band, &
                                   gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    write(*,*) '    Optical matrix elements (sp) have been read from file'
    
    !compute shift conductivity
    if (nwp.eq.1 .and. nwq.eq.(-1)) then
      call get_sigma_shift_sp(npointstotal,nband_ex,berry_eigen_ex_band, &
                        gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
      !write(*,*) 'The optical response',response_text,'has been evaluated'
    end if
    
    !compute shg susceptibility
    if (nwp.eq.1 .and. nwq.eq.1) then   
      !call get_sigma_shg(npointstotal,nband_ex,berry_eigen_ex_band, &
                        !gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    end if


  end subroutine get_sigma_second_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine get_sigma_shift_sp(npointstotal,nband_ex,berry_eigen_ex_band, &
                        gen_der_ex_band,shift_vector_ex_band,vme_ex_band,ek)
    implicit none
    !in/out
    integer :: nband_ex,npointstotal

    dimension ek(npointstotal,nband_ex)
    dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)

    real*8 :: ek
    complex*16 :: berry_eigen_ex_band
    complex*16 :: gen_der_ex_band
    real*8 :: shift_vector_ex_band
    complex*16 :: vme_ex_band

    !here
    integer :: ibz,i,ideg,j,nj,njp,njpp,nn,nnp
    integer :: idegg
    dimension e_nband(nband_ex)
    dimension vme_nband(3,nband_ex,nband_ex)
    dimension shift_vector_nband(3,3,nband_ex,nband_ex)
    dimension gen_der_nband(3,3,nband_ex,nband_ex)
    dimension vme_der_nband(3,3,nband_ex,nband_ex)
    dimension shift_vector_w(3,3,nw)
    dimension wp(nw)
    dimension sigma_w_sp(3,3,3,nw)

    real*8 :: e_nband
    complex*16 :: vme_nband
    real*8 :: shift_vector_nband
    complex*16 :: gen_der_nband
    complex*16 :: vme_der_nband
    real*8 :: shift_vector_w
    real*8 :: wp,eta2
    complex*16 :: sigma_w_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    !initialize conductivity arrays
    call initialize_sigma_second_arrays(nw,wp,eta2,sigma_w_sp)
    write(*,*) '    Evaluating shift conductivity (sp)...'
    !k-sampling and filliing frequency arrays
    do ibz=1,npointstotal
      !write(*,*) 'point:',ibz,npointstotal       
      !fill auxiliary arrays
      do i=1,nband_ex
        e_nband(i)=ek(ibz,i)

        do j=1,nband_ex
	        do nj=1,3
	          vme_nband(nj,i,j)=vme_ex_band(ibz,nj,i,j)
		        do njp=1,3
              shift_vector_nband(nj,njp,i,j)=shift_vector_ex_band(ibz,nj,njp,i,j)
              gen_der_nband(nj,njp,i,j)=gen_der_ex_band(ibz,nj,njp,i,j)
              !vme_der(nj,njp,i,j)=vme_der_ex_band(ibz,nj,njp,i,j)		  
            end do
          end do
        end do	 
      end do

      !deg discard: to be implemented
	    ideg=1
	    if (idegg.eq.1) then
        continue
		  end if
      call get_shift_intens_sp(nband_ex,nw,e_nband,vme_nband,shift_vector_nband,gen_der_nband,vme_der_nband, &
           shift_vector_w,wp,eta2,sigma_w_sp)
	  end do
    
    !get excitonic frequency tensor
    !write(*,*) 'broadening excitons peaks...'
    !print conductivity tensor
    !write(*,*) 'Printing sigma first...'
    call print_sigma_second_sp(nw,wp,sigma_w_sp,shift_vector_w)
    write(*,*) '    Shift conductivity (sp) has been printed'
  end subroutine get_sigma_shift_sp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_shift_intens_sp(nband_ex,nw,e_nband,vme_nband,shift_vector_nband,gen_der_nband,vme_der_nband, &
             shift_vector_w,wp,eta2,sigma_w_sp)
    implicit none
  
    !in/out
    integer :: nband_ex
    integer :: nw
  
    dimension :: e_nband(nband_ex)
    dimension :: vme_nband(3,nband_ex,nband_ex)
  
    dimension :: shift_vector_nband(3,3,nband_ex,nband_ex)
    dimension :: gen_der_nband(3,3,nband_ex,nband_ex)
    dimension :: vme_der_nband(3,3,nband_ex,nband_ex)
    dimension :: shift_vector_w(3,3,nw)
    dimension :: wp(nw)
    dimension :: sigma_w_sp(3,3,3,nw)

    real*8 :: e_nband
    complex*16 :: vme_nband
    real*8 :: shift_vector_nband
    complex*16 :: gen_der_nband
    complex*16 :: vme_der_nband
    real*8 :: shift_vector_w
    real*8 :: wp,eta2
    complex*16 :: sigma_w_sp

    !here
    integer :: iw,nj,njp,njpp,nn,nnp,i
    dimension :: abc(3,nband_ex,nband_ex)
    
    real*8 :: factor1,fnn,fnnp,delta_nnp
    complex*16 abc
    complex*16 shift1,shift2,shift
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    do iw=1,nw 
      do nn=1,nband_ex	
      !fermi disrtibution
        if (nn.le.nv_ex) then
          fnn=1.0d0
        else
          fnn=0.0d0
        end if		    
        do nnp=1,nband_ex
          !fermi distribution
          if (nnp.le.nv_ex) then
            fnnp=1.0d0
          else
            fnnp=0.0d0
          end if   
  	      factor1=fnn-fnnp
          if (trim(broadening_type_text) == 'gaussian') then
            delta_nnp = 1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta2**2)*(wp(iw)-e_nband(nn)+e_nband(nnp))**2)
          else if (trim(broadening_type_text) == 'lorentzian') then
            delta_nnp = 1.0d0/pi*aimag(1.0d0/(wp(iw)-e_nband(nn)+&
              e_nband(nnp)-complex(0.0d0,eta2)))
          else
            delta_nnp = 1.0d0/eta2*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta2**2)*(wp(iw)-e_nband(nn)+e_nband(nnp))**2)
          end if
  	      do nj=1,3
  	        do njp=1,3		
  		        shift_vector_w(nj,njp,iw)=shift_vector_w(nj,njp,iw)+ &
  		        1.0d0/(dble(npointstotal)*vcell)*factor1*shift_vector_nband(nj,njp,nn,nnp)*delta_nnp  
  		        do njpp=1,3		
                if (fnn.eq.fnnp) then
                  shift=0.0d0
                else		                   
                  if (response_text  == 'shift_sumrule') then
                    !normal sumrule usage for generalized derivative
  			            shift1=-complex(0.0d0,1.0d0)/(e_nband(nn)-e_nband(nnp))*vme_nband(njp,nn,nnp) &
  				                 *gen_der_nband(njpp,nj,nnp,nn)
  			            shift2=-complex(0.0d0,1.0d0)/(e_nband(nn)-e_nband(nnp))*vme_nband(njpp,nn,nnp) &
  				                 *gen_der_nband(njp,nj,nnp,nn)	
                    shift=-complex(0.0d0,1.0d0)*(shift1+shift2) 
                    
                    sigma_w_sp(nj,njp,njpp,iw)=sigma_w_sp(nj,njp,njpp,iw)+ &
  		              0.5d0*pi/(dble(npointstotal)*vcell)*factor1*shift*delta_nnp	
                  end if      
                  if (response_text  == 'shift_shiftvector') then
                    !NAGAOSA. WITH TRS by now. Shift vector	
  		              shift=-(shift_vector_nband(nj,njp,nnp,nn)-shift_vector_nband(nj,njpp,nn,nnp)) &
  			            *vme_nband(njpp,nn,nnp)*vme_nband(njp,nnp,nn)/(e_nband(nn)-e_nband(nnp))**2	  	
  		              
                    sigma_w_sp(nj,njp,njpp,iw)=sigma_w_sp(nj,njp,njpp,iw)+ &
  		              0.5d0*pi/(dble(npointstotal)*vcell)*factor1*shift*delta_nnp	
                  end if   
                  if (response_text  == 'shift_gender') then
                    !here I will put the numerical generalized derivative (see Toni's paper)
                  end if               	
  	            end if			  
                
  		        end do			
  	        end do
  	      end do
        end do  
      end do
    end do
  end subroutine get_shift_intens_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine get_sigma_shg_sp(npointstotal,nband_ex)
    implicit none
    !in/out
    integer :: nband_ex,npointstotal
    !here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  end subroutine get_sigma_shg_sp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine initialize_sigma_second_arrays(nw,wp,eta2,sigma_w)
    implicit none
  
    integer :: nw
    integer :: i

    dimension :: wp(nw)
    dimension :: sigma_w(3,3,3,nw)

    real(8) :: wrange,wp,eta2
    complex(8) :: sigma_w
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    wp=0.0d0
    sigma_w=0.0d0
    wrange=e2-e1
    do i=1,nw
      wp(i)=(e1+wrange/dble(nw)*dble(i-1))/27.211385d0
    end do  
    eta2=eta/27.211385d0 !change units to hartree units

  
  end subroutine initialize_sigma_second_arrays

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_second_sp(nw,wp,sigma_w_sp,shift_vector_w)
    implicit none
 
    !in/out
    integer :: iw
    integer :: nw
    dimension :: wp(nw)
    dimension :: sigma_w_sp(3,3,3,nw)
    dimension :: shift_vector_w(3,3,nw)
    
    real*8 :: wp
    complex*16 :: sigma_w_sp
    real*8 :: shift_vector_w
    
    !here
    real*8 :: feps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !write frequency dependent conductivity	  
    open(90,file='shift_sp_lengthgauge_'//trim(material_name)//'.dat')
    open(100,file='shift_vector.dat')
    do iw=1,nw
      !feps=-6.623618d-03/(27.21138**2)*1.0d06 !%go from au to (\mu A /V^2)*Angstrongs
      !d=2.6d0 !thickness in angstrongs for MoS2
      !d=3.28d0 !thickness in angstrongs for h-BN
      !feps=feps/(d/0.52917721067121d0) 
      feps=(6.623618d-03)*(1.0d+06)*(27.211386**(-2))*(5.291772d-11)*(1.0d+09) !%go from au to (\mu A /V^2)*nm	
      write(90,*) wp(iw)*27.211385d0,&
            realpart(feps*sigma_w_sp(1,1,1,iw)), &
		        realpart(feps*sigma_w_sp(1,2,2,iw)), &
		        realpart(feps*sigma_w_sp(2,1,1,iw)), &
		  	    realpart(feps*sigma_w_sp(2,2,2,iw)), & 
            realpart(feps*sigma_w_sp(1,1,3,iw)), &
		        realpart(feps*sigma_w_sp(2,2,3,iw)), &
		        realpart(feps*sigma_w_sp(3,1,1,iw)), &
		  	    realpart(feps*sigma_w_sp(3,2,2,iw)) 	    
      write(100,*) wp(iw)*27.211385d0,&
            shift_vector_w(1,1,iw),shift_vector_w(1,2,iw), &
            shift_vector_w(2,1,iw),shift_vector_w(2,2,iw)
    end do
    close(90)
    close(95)
    close(100)
  end subroutine print_sigma_second_sp
end module sigma_second_sp

