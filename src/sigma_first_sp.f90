module sigma_first_sp
  use constants_math
  use parser_input_file, &
  only:iflag_xatu,nf,e1,e2,eta,nw,broadening_type_text
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex, &
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use ome_ex, &
  only:read_ome_sp_linear !routine

  implicit none
  
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_first_sp()
    implicit none 
    integer iflag_norder
    integer :: ibz,j
    
    !energies and vme in k-mesh and auxiliary arrays (sp)
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension vme_nband(3,nband_ex,nband_ex)
    dimension e_nband(nband_ex)
    
    dimension wp(nw)
    dimension sigma_w_sp(3,3,nw)
    
    real*8 wp,eta1
    real*8 ek,e_nband
    complex*16 vme_ex_band,vme_nband 
    complex*16 sigma_w_sp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '8. Entering sigma_first_sp'
    !initialize sp arrays
    vme_ex_band=0.0d0
    ek=0.0d0   
    e_nband=0.0d0
    vme_nband=0.0d0
    !initialize ex arrays
    
    !read optical matrix elements from file
    write(*,*) '   Reading optical matrix elements...'
    call read_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
    !allocate conductivity arrays
    call fill_allocate_sigma_arrays(eta1,nw,wp,sigma_w_sp)
    
    write(*,*) '   Evaluating linear conductivity (sp)...'
 
    do ibz=1,npointstotal
      !write(*,*) 'sigma_sp point:',ibz,npointstotal       
      !fill auxiliary arrays
      e_nband(:)=ek(ibz,:)
      vme_nband(:,:,:)=vme_ex_band(ibz,:,:,:)
      !fill sigma(w) for a given k point
      call get_kubo_intens_sp(nband_ex,npointstotal,vcell,e_nband,vme_nband,nw,wp,eta1,sigma_w_sp)
    end do
    !close(50)

    !print conductivity tensor
    write(*,*) '   Printing sigma first...'
    call print_sigma_first_sp(nw,wp,sigma_w_sp)

  end subroutine get_sigma_first_sp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine fill_allocate_sigma_arrays(eta1,nw,wp,sigma_w)
  implicit none
  
  integer :: nw
  integer :: i

  dimension :: wp(nw)
  dimension :: sigma_w(3,3,nw)

  real(8) :: wrange,wp,eta1
  complex(8) :: sigma_w
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  wp=0.0d0
  sigma_w=0.0d0
  wrange=e2-e1
  do i=1,nw
    wp(i)=(e1+wrange/dble(nw)*dble(i-1))/27.211385d0
  end do  
  eta1=eta/27.211385d0 !change units to hartree units
  
  end subroutine fill_allocate_sigma_arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_first_sp(nw,wp,sigma_w_sp)
    implicit none
    integer :: iw
    integer :: nw
    dimension :: wp(nw)
    dimension :: sigma_w_sp(3,3,nw)
    
    real*8 :: wp,feps
    complex*16 :: sigma_w_sp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !write frequency dependent conductivity	  
    open(50,file='sigma_first_sp_real_'//trim(material_name)//'.dat')
    open(55,file='sigma_first_sp_imag_'//trim(material_name)//'.dat')

    do iw=1,nw
      feps=1.0d0 !use atomic units
      write(50,*) wp(iw)*27.211385d0, &
        realpart(feps*sigma_w_sp(1,1,iw)), &
        realpart(feps*sigma_w_sp(1,2,iw)), &
        realpart(feps*sigma_w_sp(1,3,iw)), &
        realpart(feps*sigma_w_sp(2,1,iw)), &
        realpart(feps*sigma_w_sp(2,2,iw)), &
        realpart(feps*sigma_w_sp(2,3,iw)), &
        realpart(feps*sigma_w_sp(3,1,iw)), &
        realpart(feps*sigma_w_sp(3,2,iw)), &
        realpart(feps*sigma_w_sp(3,3,iw))
  
      write(55,*) wp(iw)*27.211385d0,aimag(feps*sigma_w_sp(1,1,iw)), &
          aimag(feps*sigma_w_sp(1,2,iw)), &
          aimag(feps*sigma_w_sp(1,3,iw)), &
          aimag(feps*sigma_w_sp(2,1,iw)), &
          aimag(feps*sigma_w_sp(2,2,iw)), &
          aimag(feps*sigma_w_sp(2,3,iw)), &
          aimag(feps*sigma_w_sp(3,1,iw)), &
          aimag(feps*sigma_w_sp(3,2,iw)), &
          aimag(feps*sigma_w_sp(3,3,iw))	
      
    end do

    close(50)
    close(55)

  end subroutine print_sigma_first_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_kubo_intens_sp(nband_ex,npointstotal,vcell,e,vme,nw,wp,eta1,sigma_w_sp)
    implicit none
    integer :: iw,nn,nnp,nj,njp
    integer :: nw
    integer :: nband_ex,npointstotal

    dimension :: vme(3,nband_ex,nband_ex)
    dimension :: e(nband_ex)
    dimension :: sigma_w_sp(3,3,nw)

    dimension :: wp(nw)

    real*8 :: fnn,fnnp,factor1,vcell
    real*8 :: wp,eta1
    real*8 :: e
    
    complex*16 vme,delta_nnp
    complex*16 vme_prod,sigma_w_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  

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
          !DEDICE PREFACTOR WITH OCCUPATION
          if (abs(fnn-fnnp).lt.0.1d0) then 
            factor1=0.0d0         
          else
            factor1=(fnn-fnnp)/(e(nn)-e(nnp))
          end if
          ! Broadening: gaussian or lorentzian based on parser input
          if (trim(broadening_type_text) == 'gaussian') then
            delta_nnp = -1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta1**2)*(wp(iw)-e(nn)+e(nnp))**2)
          else if (trim(broadening_type_text) == 'lorentzian') then
            delta_nnp = 1.0d0/pi*aimag(1.0d0/(-wp(iw)+e(nn)-e(nnp)+&
              complex(0.0d0,eta1)))
          else
            delta_nnp = -1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta1**2)*(wp(iw)-e(nn)+e(nnp))**2)
          end if

          !save oscillator stregths
          do nj=1,3
            do njp=1,3   
              vme_prod=vme(nj,nn,nnp)*vme(njp,nnp,nn) 
              sigma_w_sp(nj,njp,iw)=sigma_w_sp(nj,njp,iw)+ &
    	        pi/(dble(npointstotal)*vcell)*factor1*vme_prod*delta_nnp 
            end do
          end do
    	  	  
        end do
      end do
    end do  

  end subroutine get_kubo_intens_sp
end module sigma_first_sp

