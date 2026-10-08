module sigma_first_ex
  use constants_math
  use parser_input_file, &
  only:iflag_xatu,nf,e1,e2,eta,nw,broadening_type_text
  use parser_wannier90_tb, &
  only:material_name,norb
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  norb_ex_cut,nv_ex,nc_ex,nband_ex, &
  e_ex,fk_ex, &
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use ome_ex, &
  only:read_ome_sp_linear !routine
  use sigma_first_sp, &
  only:fill_allocate_sigma_arrays
  
  implicit none

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_sigma_first_ex()
    implicit none 
    integer iflag_norder
    integer :: ibz,j
    
    !energies and vme in k-mesh and auxiliary arrays (sp)
    dimension ek(npointstotal,nband_ex)
    dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
    dimension vme_nband(3,nband_ex,nband_ex)
    dimension e_nband(nband_ex)
    
    !energies and VME (ex)
    dimension e_ex(norb_ex_cut)
    dimension vme_ex(3,norb_ex_cut)

    dimension wp(nw)
    dimension sigma_w_sp(3,3,nw),sigma_w_ex(3,3,nw)
    
    real*8 wp,eta1
    real*8 ek,e_nband
    real*8 e_ex
    complex*16 vme_ex_band,vme_nband 
    complex*16 vme_ex
    complex*16 sigma_w_sp,sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '9. Entering sigma_first_ex'
    !initialize ex arrays
    vme_ex=0.0d0
    call read_ome_ex_linear(vme_ex)

    !allocate conductivity arrays
    call fill_allocate_sigma_arrays(eta1,nw,wp,sigma_w_ex)
    
    write(*,*) '   Evaluating linear conductivity (ex)...'
    !get excitonic frequency tensor
    call get_kubo_intens_ex(vme_ex,nw,wp,eta1,sigma_w_ex)

    !print conductivity tensor
    write(*,*) '   Printing sigma first (ex)...'
    call print_sigma_first_ex(nw,wp,sigma_w_ex)

  end subroutine get_sigma_first_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ome_ex_linear(vme_ex)
    implicit none
    integer :: iounit10
    integer nn,nkaka
    dimension vme_ex(3,norb_ex_cut)
    complex*16 vme_ex

    real*8 :: a1,a2,a3,a4,a5,a6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(newunit=iounit10,file='ome_linear_ex_'//trim(material_name)//'.omeex') 
    read(iounit10,*)     
    do nn=1,norb_ex_cut
      read(iounit10,*) nkaka,a1,a2,a3,a4,a5,a6
      vme_ex(1,nn)=complex(a1,a2)
      vme_ex(2,nn)=complex(a3,a4)
      vme_ex(3,nn)=complex(a5,a6)
    end do
    close(iounit10)
  end subroutine read_ome_ex_linear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FLAG (physics, unchanged for now): vme_ex is the BARE momentum matrix element
  ! P_n = sum psi * v_cv (Taghizadeh & Pedersen, PRB 97, 205432, Eq. B2a), so this is
  ! sigma = pi/(Nk V) |P_n|^2/E_n delta(w-E_n), the bare-momentum ("C'") response. It is identical to
  ! Xatu's skubo_w.f90 (checked on hBN, 60x60 grid, eta = 0.08 eV: agreement 3e-7), but the paper's
  ! methods A/B use the Heisenberg momentum Pi_n = -i E_n X_n = P_n - i F_n (Eq. 10; F_n != 0 with
  ! e-h interaction), i.e. E_n |X_n|^2 instead of |P_n|^2/E_n. On hBN this makes the first peak
  ! (P/Pi)^2 = 2.34 times too large (3.147 vs 1.345 a.u.; 2.1-2.15 above 6.6 eV; the paper's Fig. 1
  ! shows the same ~2.1x C' vs A-D). The full-omega method B (Eq. 4b) matches the Pi form to 3e-4.
  ! Not switched yet: needs X_n, which is only filled when the nonlinear matrix elements are
  ! requested.
  subroutine get_kubo_intens_ex(vme_ex,nw,wp,eta1,sigma_w_ex)
    implicit none
    !dimension skubo_ex_int(3,3,norb_ex_cut)
    
    dimension wp(nw),sigma_w_ex(3,3,nw)
    dimension vme_ex(3,norb_ex_cut)
    
    integer nw
    integer iw,nn,nj,njp
    real*8 delta_n_ex
    real*8 wp,eta1
    
    complex*16 :: vme_ex
    complex*16 :: skubo_ex_int, sigma_w_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
    skubo_ex_int=0.0d0
    sigma_w_ex=0.0d0
    !$omp parallel do default(none) schedule(static) private(nn,nj,njp,iw,skubo_ex_int,delta_n_ex) &
    !$omp   shared(npointstotal, norb_ex_cut, nw, vcell, e_ex, vme_ex, wp, eta1, broadening_type_text) &
    !$omp   reduction(+:sigma_w_ex)
    do nn=1,norb_ex_cut
      do nj=1,3
        do njp=1,3
            
            !N integrand
            skubo_ex_int=pi/(dble(npointstotal)*vcell) &
            *conjg(vme_ex(nj,nn))*vme_ex(njp,nn)/e_ex(nn)   !pick the correct order of operators
            
            
          do iw=1,nw 
            !at a given frequency
            !delta function
            if (trim(broadening_type_text) == 'gaussian') then
              !delta_n_ex = pi*1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta1**2)*(wp(iw)-e_ex(nn))**2)
              delta_n_ex = 1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta1**2)*(wp(iw)-e_ex(nn))**2)
            else if (trim(broadening_type_text) == 'lorentzian') then
              delta_n_ex = 1.0d0/pi*aimag(1.0d0/(wp(iw)-e_ex(nn)-complex(0.0d0,eta1)))
                
            else
              !delta_n_ex = pi*1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta1**2)*(wp(iw)-e_ex(nn))**2)
!              delta_n_ex = 1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta1**2)*(wp(iw)-e_ex(nn))**2)
              delta_n_ex = 1.0d0/pi*aimag(1.0d0/(wp(iw)-e_ex(nn)-complex(0.0d0,eta1)))
            end if
            
            !sigma_w
            sigma_w_ex(nj,njp,iw)=sigma_w_ex(nj,njp,iw)+skubo_ex_int*delta_n_ex
          
          end do
        end do
      end do
    end do  

  end subroutine get_kubo_intens_ex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_sigma_first_ex(nw,wp,sigma_w_ex)
    implicit none
    integer :: iounit55
    integer :: iounit50
    integer :: iw
    integer :: nw
    dimension :: wp(nw)
    dimension :: sigma_w_ex(3,3,nw)
    
    real*8 :: wp,feps
    complex*16 :: sigma_w_ex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    open(newunit=iounit50,file='sigma_first_ex_real_'//trim(material_name)//'.dat')
    open(newunit=iounit55,file='sigma_first_ex_imag_'//trim(material_name)//'.dat')
    
    feps=1.0d0 !use atomic units
    ! serial on purpose (audit 2026-09-24): this loop is pure file I/O and every iteration was
    ! inside !$omp ordered, which serialises it completely -- the parallel wrapper only added
    ! thread spawn and synchronisation cost..
    do iw=1,nw
      write(iounit50,*) wp(iw)*27.211385d0, &
        realpart(feps*sigma_w_ex(1,1,iw)), &
        realpart(feps*sigma_w_ex(1,2,iw)), &
        realpart(feps*sigma_w_ex(1,3,iw)), &
        realpart(feps*sigma_w_ex(2,1,iw)), &
        realpart(feps*sigma_w_ex(2,2,iw)), &
        realpart(feps*sigma_w_ex(2,3,iw)), &
        realpart(feps*sigma_w_ex(3,1,iw)), &
        realpart(feps*sigma_w_ex(3,2,iw)), &
        realpart(feps*sigma_w_ex(3,3,iw))
  
      write(iounit55,*) wp(iw)*27.211385d0, &
          aimag(feps*sigma_w_ex(1,1,iw)), &
          aimag(feps*sigma_w_ex(1,2,iw)), &
          aimag(feps*sigma_w_ex(1,3,iw)), &
          aimag(feps*sigma_w_ex(2,1,iw)), &
          aimag(feps*sigma_w_ex(2,2,iw)), &
          aimag(feps*sigma_w_ex(2,3,iw)), &
          aimag(feps*sigma_w_ex(3,1,iw)), &
          aimag(feps*sigma_w_ex(3,2,iw)), &
          aimag(feps*sigma_w_ex(3,3,iw))	
    end do

    close(iounit50)
    close(iounit55)

  end subroutine print_sigma_first_ex

end module sigma_first_ex

