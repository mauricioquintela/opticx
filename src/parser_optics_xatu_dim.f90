module parser_optics_xatu_dim
  use constants_math
  use parser_wannier90_tb, &
    only:material_name,R,nRvec,norb !variables
  use parser_input_file, &
    only:xatu_eigval_filepath_in,xatu_states_filepath_in, & !filepaths
      ndim,npointstotal_sq, & !variables
      iflag_xatu,nf,nband_index,norb_ex_cut, & 
      cache_ome_read,iflag_ome_ex_text,iflag_write_exk, & !for the fk_ex deferral, below
      read_line_numbers_int !subroutine
  implicit none

  integer :: nv_ex,nc_ex
  integer :: npointstotal
  integer :: norb_ex
  integer :: norb_ex_band
  integer :: nband_ex
  integer :: naux
  integer :: j

  real(8) G,vcell
  real(8) rkxvector,rkyvector,rkzvector
  real(8) auxr1
  real(8) e_ex
  complex*16 fk_ex

  ! .true. once fk_ex is known to sit in the SAME single-particle basis the response routines use:
  ! either ome_sp carried it into the Eq. (A4) rotated basis (rotate_fk_ex_to_a4_basis), or that
  ! rotation is switched off and there is nothing to carry. Anything that consumes fk_ex must refuse
  ! to run while this is .false., because the answer is then wrong by O(1) and not by a little
  ! (HANDOFF 8.46). It lives here, next to fk_ex, so that ome_sp (which sets it) and ome_ex (which
  ! checks it) share it without either module having to use the other.
  logical :: fk_ex_basis_ok = .false.

  ! .true. once the exciton envelopes have actually been read from the .states file. Reading them is
  ! the second most expensive part of start-up (measured 16.8 s of a 22 s fixed cost on ReS2 with 2700
  ! excitons: 1.9 GB of ASCII), and a second-order OME cache HIT never touches them -- the cached
  ! matrix elements are what the k-loop would have produced from fk_ex. So the read is DEFERRED when a
  ! hit is possible and done lazily by load_fk_ex() if the cache turns out to miss, which costs exactly
  ! what it costs today. Anything that consumes fk_ex must check this first.
  logical :: fk_ex_loaded = .false.

  dimension G(3,3)

  allocatable rkxvector(:)
  allocatable rkyvector(:)
  allocatable rkzvector(:)
  allocatable fk_ex(:,:)
  allocatable e_ex(:)
  allocatable auxr1(:)
	  
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here we define some BZ variables, either by reading the output
! of Xatu or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_optics_xatu_dim()
  implicit none  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) '3. Entering parser_optics_xatu_dim'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Reminder:norb_ex_cut is given as input
  !get band and grid dimensions: from XATU-output or opticx-input
  if (iflag_xatu .eqv. .true.) then
    call get_exciton_dim()
    !nband_ex=nc_ex+nv_ex
  else
    !number of total k-points
    if (ndim==1) npointstotal=npointstotal_sq
    if (ndim==2) npointstotal=npointstotal_sq**2
    if (ndim==3) npointstotal=npointstotal_sq**3
    norb_ex_band=nv_ex*nc_ex
    norb_ex=norb_ex_band*npointstotal   
  end if
  !calculate nv_ex and nc_ex from the array of bands
  nband_ex=size(nband_index,dim=1)
  nv_ex=0
  do j=1,nband_ex
    if (nband_index(j).le.0) nv_ex=nv_ex+1
  end do
  nc_ex=nband_ex-nv_ex

  !change syntax for band counting
  !XATU: ...-1 0 1 2... to explicit band count
  !opticx: ...nf-1,nf,nf+1...
  nband_index(:)=nband_index(:)+nf

  ! Validate every Bandlist entry (shifted by Nfermi above) against the actual orbital count.
  ! Without this, an out-of-range entry silently indexes e(:)/vme(:,:,:) out of bounds in
  ! get_ome_sp (ome_sp.f90) -- wrong energies/matrix elements with no diagnostic, or a crash
  ! only under -fcheck=all with no message identifying the actual misconfigured entry.
  do j=1,nband_ex
    if (nband_index(j)<1 .or. nband_index(j)>norb) then
      write(*,*) 'ERROR: Bandlist entry',j,'resolves to band',nband_index(j), &
                  ', outside the valid range 1..',norb,'(check Bandlist and Nfermi).'
      stop 1
    end if
  end do

  !allocate grid and exciton arrays
  allocate (rkxvector(npointstotal))
  allocate (rkyvector(npointstotal))
  allocate (rkzvector(npointstotal))
  allocate (e_ex(norb_ex_cut))
  allocate (fk_ex(norb_ex,norb_ex_cut))
  
  !get reciprocal lattice vectors
  call get_reciprocal_vectors()  
  
  !fill exciton and other arrays: from XATU-output or opticx-input
  if (iflag_xatu .eqv. .true.) then
    call get_exciton_data() !get grid and exciton wavefunctions
    write(*,*) '   Exciton data has been read from XATU output'
  else
    call get_grid()
    !get grid and exciton variables set to zero if XATU interface is not requested
    fk_ex=0.0d0
    e_ex=0.0d0
  end if
  write(*,*) "   Grid and band parameters have been set"
  
end subroutine get_optics_xatu_dim   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! This subroutine prints a part of the exciton wavefunction or 
! the total exciton probability density
subroutine print_exciton_wf(isum,iv,ic,nn)
  implicit none
  integer :: iounit10
  integer isum,iv,ic,nn
  integer iv_s,ic_s
  integer iright,i_ex_nn
  integer ibz
  real*8 prob_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  iright=0
  open(newunit=iounit10,file='exciton_wf.dat')
  do ibz=1,npointstotal     
    if (isum.eq.1) then
      prob_k=0.0d0
      do iv_s=1,nv_ex
        do ic_s=1,nc_ex
          call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic_s,iv_s)
          prob_k=prob_k+abs(fk_ex(i_ex_nn,nn))**2
        end do
      end do
      write(iounit10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),prob_k
    else
      call get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex_nn,ic,iv)
      write(iounit10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),abs(fk_ex(i_ex_nn,nn))
    end if
  end do
  close(iounit10)
end subroutine print_exciton_wf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!this subroutine has not been updated in 2025
subroutine get_ex_index_first(nf,nv_ex,nc_ex,iright,ibz,i_ex,ic,iv)
  implicit none
  integer nf,nv_ex,nc_ex,ibz,i_ex,ic,iv
  integer iright
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !get band indeces (respect the fermi level) from A_cv index
    if (iright.eq.1) then
      iv=(nf-nv_ex)+i_ex-int((i_ex-1)/nv_ex)*nv_ex-nf
      ic=(nf+1)+int((i_ex-1)/nv_ex)-nf
    end if
    !get A_cv index from band indeces
    if (iright.ne.1) then
      i_ex=nc_ex*nv_ex*(ibz-1)+nv_ex*(ic-1)+iv
    end if
end subroutine get_ex_index_first

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_exciton_dim()
  implicit none
  integer, parameter :: max_scan = 1000
  dimension nband_index_aux1(max_scan)
  dimension nband_index_aux2(max_scan)

  integer :: nband_ex_aux
  integer :: nband_index_aux1,nband_index_aux2
  integer :: nband_ex_aux1,nband_ex_aux2
  integer :: iexit
  integer :: i,j,naux,npointstotal_sq
  integer :: hdr1, hdr2, ios
  integer :: iounit_bands, iounit_nk
  integer :: nv_ex_local, nc_ex_local, norb_ex_band_local

  real(8) aux1
  character(len=:), allocatable :: file2open
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !This part gets 'nband_ex' and 'nband_index(nband_ex)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !save number of valence bands
  file2open=trim(xatu_states_filepath_in)
  open(newunit=iounit_bands,file=file2open)
  read(iounit_bands,*) 
  do i=1,max_scan/2
    read(iounit_bands,*) aux1,aux1,aux1,nband_index_aux1(i)
    if (i.gt.1) then
      do j=1,i-1
        if (nband_index_aux1(i).eq.nband_index_aux1(j)) then
          nband_ex_aux1=i-1
          goto 128
        end if
      end do
    end if
  end do
  write(*,*) 'ERROR (get_exciton_dim): no repeated valence-band index found within ', &
            max_scan/2, ' lines.'
  stop 1
  128   continue

  !save number of conduction bands
  ! rewind, NOT a second open: re-opening a unit already connected to the SAME file leaves the file
  ! POSITION untouched (F2008 9.5.6.1, and gfortran behaves that way -- verified), so the open that
  ! used to sit here never did the rewind it looks like it does. The scan below happened to survive
  ! that because the basis list is periodic, so a stride-nv scan finds the same repeat from any
  ! offset -- checked against all four .states files. Relying on that is not worth a saved line.
  rewind(iounit_bands)
  read(iounit_bands,*) 
  do i=1,max_scan/2
    read(iounit_bands,*) aux1,aux1,aux1,naux,nband_index_aux2(i) 
    if (nband_ex_aux1.gt.1) then
      do j=1,nband_ex_aux1-1
        read(iounit_bands,*)
      end do
    end if
    
    if (i.gt.1) then
      do j=1,i-1
        if (nband_index_aux2(i).eq.nband_index_aux2(j)) then
          nband_ex_aux2=i-1
          goto 129
        end if
      end do
    end if
  end do
  write(*,*) 'ERROR (get_exciton_dim): no repeated conudction-band index found within ', &
            max_scan/2, ' lines.'
  stop 1
  129   continue
  close(iounit_bands)	 

  nband_ex_aux=nband_ex_aux1+nband_ex_aux2
  allocate(nband_index(nband_ex_aux))
  do i=1,nband_ex_aux1
    nband_index(i)=nband_index_aux1(i)-nf+1
  end do
  do i=nband_ex_aux1+1,nband_ex_aux
    nband_index(i)=nband_index_aux2(i-nband_ex_aux1)-nf+1
  end do

  nv_ex_local=0
  do i=1,nband_ex_aux
    if (nband_index(i).le.0) nv_ex_local=nv_ex_local+1
  end do
  nc_ex_local=nband_ex_aux-nv_ex_local
  norb_ex_band_local=nv_ex_local*nc_ex_local

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !get nk
  file2open=trim(xatu_eigval_filepath_in)
  open(newunit=iounit_nk,file=file2open)
  read(iounit_nk,*,iostat=ios) hdr1
  if (ios == 0) then
    read(iounit_nk,*,iostat=ios) hdr2
    if (ios == 0) then
      npointstotal_sq = hdr1
      naux = hdr2
    else
      ! New xatu Result::writeEigenvalues format: first line = exciton basis size (naux)
      naux = hdr1
      if (norb_ex_band_local > 0) then
        npointstotal = naux / norb_ex_band_local
      else
        npointstotal = 0
      end if
      ! derive npointstotal_sq based on ndim
      if (ndim == 1) then
        npointstotal_sq = npointstotal
      else if (ndim == 2) then
        npointstotal_sq = int(sqrt(dble(npointstotal)) + 0.5d0)
      else if (ndim == 3) then
        npointstotal_sq = int((dble(npointstotal))**(1.0d0/3.0d0) + 0.5d0)
      end if
    end if
  else
    rewind(iounit_nk)
    read(iounit_nk,*) npointstotal_sq
    read(iounit_nk,*) naux
  end if
  close(iounit_nk)

  !get N_BSE=nv_ex*nc_ex*nk variables
  if (npointstotal == 0) then
    if (ndim==1) npointstotal=npointstotal_sq
    if (ndim==2) npointstotal=npointstotal_sq**2
    if (ndim==3) npointstotal=npointstotal_sq**3
  end if
  if (npointstotal > 0) then
    norb_ex_band = int(naux / npointstotal)
  else
    norb_ex_band = 0
  end if
  norb_ex = norb_ex_band * npointstotal

end subroutine get_exciton_dim
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_reciprocal_vectors()
  implicit none
  real(8) cx,cy,cz
  real*8 det      
  logical :: active_x, active_y, active_z
  active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
  active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
  active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  G=0.0d0

  ! 1D
  if ( ndim == 1 ) then
        
    if (active_z) then
      G(3,3)=2.0d0*pi*(R(3,3))**(-1.0d0)
      vcell=sqrt(R(3,1)**2+R(3,2)**2+R(3,3)**2)
      
    elseif (active_y) then
      G(2,2)=2.0d0*pi*(R(2,2))**(-1.0d0)
      vcell=sqrt(R(2,1)**2+R(2,2)**2+R(2,3)**2)
      
    else
      G(1,1)=2.0d0*pi*(R(1,1))**(-1.0d0)
      vcell=sqrt(R(1,1)**2+R(1,2)**2+R(1,3)**2)
    endif
  
  ! 2D
  elseif ( ndim == 2 ) then
    if (active_y .and. active_z) then
      G(2,2)=2.0d0*pi*(-R(2,3)*R(3,2)+R(2,2)*R(3,3))**(-1.0d0) &
          *(R(3,3))
      G(2,3)=2.0d0*pi*(-R(2,2)*R(3,3)+R(2,3)*R(3,2))**(-1.0d0) &
          *(R(3,2))
      G(3,2)=2.0d0*pi*(-R(2,2)*R(3,3)+R(2,3)*R(3,2))**(-1.0d0) &
          *(R(2,3))
      G(3,3)=2.0d0*pi*(-R(2,3)*R(3,2)+R(2,2)*R(3,3))**(-1.0d0) &
          *(R(2,2))
      call crossproduct(R(3,1),R(3,2),R(3,3),R(2,1),R(2,2),R(2,3),cx,cy,cz)
      
    elseif (active_x .and. active_z) then
      G(1,1)=2.0d0*pi*(-R(1,3)*R(3,1)+R(1,1)*R(3,3))**(-1.0d0) &
          *(R(3,3))
      G(1,3)=2.0d0*pi*(-R(1,1)*R(3,3)+R(1,3)*R(3,1))**(-1.0d0) &
          *(R(3,1))
      G(3,1)=2.0d0*pi*(-R(1,1)*R(3,3)+R(1,3)*R(3,1))**(-1.0d0) &
          *(R(1,3))
      G(3,3)=2.0d0*pi*(-R(1,3)*R(3,1)+R(1,1)*R(3,3))**(-1.0d0) &
          *(R(1,1))
      call crossproduct(R(1,1),R(1,2),R(1,3),R(3,1),R(3,2),R(3,3),cx,cy,cz)
      
    else
      G(1,1)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
          *(-R(2,2))
      G(1,2)=2.0d0*pi*(-R(1,1)*R(2,2)+R(1,2)*R(2,1))**(-1.0d0) &
          *R(2,1)            
      G(2,1)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
          *(-R(1,2))
      G(2,2)=2.0d0*pi*(-R(2,1)*R(1,2)+R(2,2)*R(1,1))**(-1.0d0) &
          *R(1,1) 
      call crossproduct(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),cx,cy,cz)
    endif
    
    vcell=sqrt(cx**2+cy**2+cz**2)
    
  ! 3D
  elseif ( ndim == 3 ) then

    ! determinant of 3x3 matrix
    det = R(1,3)*R(2,2)*R(3,1) - R(1,2)*R(2,3)*R(3,1) - R(1,3)*R(2,1)*R(3,2) &
        + R(1,1)*R(2,3)*R(3,2) + R(1,2)*R(2,1)*R(3,3) - R(1,1)*R(2,2)*R(3,3)


    G(1,1)=2.0d0*pi*(det)**(-1.0d0)*( R(2,3)*R(3,2) - R(2,2)*R(3,3) )
    G(1,2)=2.0d0*pi*(det)**(-1.0d0)*(-R(2,3)*R(3,1) + R(2,1)*R(3,3) )
    G(1,3)=2.0d0*pi*(det)**(-1.0d0)*( R(2,2)*R(3,1) - R(2,1)*R(3,2) )
    G(2,1)=2.0d0*pi*(det)**(-1.0d0)*(-R(1,3)*R(3,2) + R(1,2)*R(3,3) )
    G(2,2)=2.0d0*pi*(det)**(-1.0d0)*( R(1,3)*R(3,1) - R(1,1)*R(3,3) )
    G(2,3)=2.0d0*pi*(det)**(-1.0d0)*(-R(1,2)*R(3,1) + R(1,1)*R(3,2) )
    G(3,1)=2.0d0*pi*(det)**(-1.0d0)*( R(1,3)*R(2,2) - R(1,2)*R(2,3) )
    G(3,2)=2.0d0*pi*(det)**(-1.0d0)*(-R(1,3)*R(2,1) + R(1,1)*R(2,3) )
    G(3,3)=2.0d0*pi*(det)**(-1.0d0)*( R(1,2)*R(2,1) - R(1,1)*R(2,2) )
  
    call crossproduct(R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3),cx,cy,cz)
    ! norm of triple product
    vcell=abs(R(1,1)*cx+R(1,2)*cy+R(1,3)*cz)
    
  endif
    
end subroutine get_reciprocal_vectors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
subroutine get_exciton_data()
  implicit none
  integer j,nkaka,nread
  integer ib,ibz,ibz_sum,jind  
  real(8) auxr1

  dimension auxr1(2*norb_ex)
  character(len=:), allocatable :: file2open
  integer :: header1, ios
  integer :: iounit_eex, iounit_kgrid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  
  !get energies
  file2open=trim(xatu_eigval_filepath_in)
  open(newunit=iounit_eex,file=file2open) 
  read(iounit_eex,*,iostat=ios) header1
  read(iounit_eex,*,iostat=ios)
  if (ios == 0) then
    read(iounit_eex,*,iostat=ios) nkaka
    if (ios == 0) then
      nread = 0
      do j=1,norb_ex_cut
        read(iounit_eex,*,iostat=ios) e_ex(j)
        if (ios /= 0) exit
        nread = nread + 1
      end do
      ! Asking for more excitons than Xatu wrote used to run off the end of the file and die inside
      ! load_fk_ex with a bare "Fortran runtime error: End of file" and a backtrace, naming neither the
      ! keyword nor the limit. The count is right here, so say so.
      if (nread < norb_ex_cut) then
        write(*,*) 'ERROR (get_exciton_data): Exciton_cutoff =', norb_ex_cut, 'but only', nread
        write(*,*) '       exciton energies are present in ', trim(file2open)
        write(*,*) '       Lower Exciton_cutoff to at most', nread, ', or rerun Xatu with a larger -n.'
        stop 1
      end if
    else
      rewind(iounit_eex)
      read(iounit_eex,*)
      read(iounit_eex,*) nkaka,(e_ex(j), j=1,norb_ex_cut)
    end if
  else
    rewind(iounit_eex)
    read(iounit_eex,*)
    read(iounit_eex,*) nkaka,(e_ex(j), j=1,norb_ex_cut)
  end if
  close(iounit_eex)

  file2open=trim(xatu_states_filepath_in)
    open(newunit=iounit_kgrid,file=file2open)	  	  
    read(iounit_kgrid,*) 

    !reading k-mesh
    do ibz=1,npointstotal
        read(iounit_kgrid,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz)
        do ib=1,norb_ex_band-1
          read(iounit_kgrid,*) 
      end do
    end do
  close(iounit_kgrid)

  !Please I like to work in atomic units!  	
  e_ex=e_ex/27.211385d0
  rkxvector=rkxvector*0.52917721067121d0 
  rkyvector=rkyvector*0.52917721067121d0 
  rkzvector=rkzvector*0.52917721067121d0 

  ! The envelopes themselves are read here only if something is going to use them. A second-order OME
  ! cache read may make them unnecessary; in that case load_fk_ex() is called later, by get_ome_ex, if
  ! and only if the cache misses. The predicate is deliberately coarse -- it asks whether a hit is
  ! POSSIBLE, not whether it will happen -- because a wrong guess costs only the same read, later.
  if (.not. (cache_ome_read .and. iflag_ome_ex_text == 'nonlinear' .and. .not. iflag_write_exk)) then
    call load_fk_ex()
  else
    write(*,*) '   Exciton envelopes not read yet: a second-order cache may make them unnecessary'
  end if

end subroutine get_exciton_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the exciton envelopes fk_ex from the .states file. Split out of get_exciton_data so that it can
! be skipped when a second-order OME cache hit will supply the excitonic matrix elements directly, and
! run on demand if that cache misses. Re-opens the file and skips the k-mesh block (norb_ex lines),
! which is cheap next to the norb_ex_cut wavefunctions that follow.
subroutine load_fk_ex()
  implicit none
  integer :: iounit10
  integer j,ib,ibz,ibz_sum,jind,nkaka,ios
  real(8) auxr1
  dimension auxr1(2*norb_ex)
  character(len=:), allocatable :: file2open

  if (fk_ex_loaded) return
  nkaka = 0

  file2open=trim(xatu_states_filepath_in)
  open(newunit=iounit10,file=file2open)
  read(iounit10,*)
  do ibz=1,npointstotal          ! skip the k-mesh block read by get_exciton_data
    do ib=1,norb_ex_band
      read(iounit10,*)
    end do
  end do

  ibz_sum=0
  write(*,*) '   Reading exciton wavefunctions...'
  do ibz=1,norb_ex_cut
    ibz_sum=ibz_sum+1
    if (abs(dble(ibz)/dble(norb_ex_cut))*100.0d0-100.0d0 .lt. 5.0d0) then
      call percentage_index(ibz_sum,norb_ex_cut,nkaka)
    end if
    read(iounit10,*,iostat=ios) (auxr1(j),j=1,2*norb_ex)
    if (ios /= 0) then
      write(*,*) 'ERROR (load_fk_ex): Exciton_cutoff =', norb_ex_cut, 'but ', trim(file2open)
      write(*,*) '       holds only', ibz-1, 'exciton wavefunctions (the .eigval file listed more).'
      write(*,*) '       Lower Exciton_cutoff, or rerun Xatu writing the states you asked for.'
      stop 1
    end if
    do j=1,norb_ex
      jind=2*j-1
      fk_ex(j,ibz)=complex(auxr1(jind),auxr1(jind+1))
    end do
  end do
  close(iounit10)
  fk_ex_loaded = .true.

end subroutine load_fk_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_grid()
  implicit none
 
  integer :: k,i1,i2,i3
  real(8) :: step,offset,r1,r2,r3
  logical :: active_x, active_y, active_z
  active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
  active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
  active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Same mesh as Xatu (Lattice::brillouinZoneMesh, Monkhorst-Pack, Gamma-centred): along each reciprocal axis
  ! u = c/n - 1/2, c = 0..n-1, plus 1/2 added to c for odd n (u = (c+1/2)/n - 1/2, symmetric, no zone edge). So the
  ! spacing is 1/n and one zone edge (-1/2, even n only) is present once; it used to be 1/(n-1) with both edges
  ! included (a non-periodic mesh with an O(1/n) error, e.g. 12% in the hBN shift conductivity at n = 30).
  ! The first reciprocal axis varies fastest, as in Xatu. Verified against the Xatu .states files for n = 30, 45, 60.
  step=1.0d0/dble(npointstotal_sq)
  offset=0.0d0
  if (mod(npointstotal_sq,2) == 1) offset=0.5d0

  ! 1D
  if ( ndim == 1 ) then
 
    if ( active_z ) then
      k=1
      do i1=1,npointstotal_sq
        r1=-0.5d0+(dble(i1-1)+offset)*step
        rkxvector(k)=0.0d0
        rkyvector(k)=0.0d0
        rkzvector(k)=r1*G(3,3)  
        k=k+1       
      end do
      
    elseif ( active_y ) then
      k=1
      do i1=1,npointstotal_sq
        r1=-0.5d0+(dble(i1-1)+offset)*step
        rkxvector(k)=0.0d0
        rkyvector(k)=r1*G(2,2)  
        rkzvector(k)=0.0d0
        k=k+1       
      end do
      
    else
      k=1
      do i1=1,npointstotal_sq
        r1=-0.5d0+(dble(i1-1)+offset)*step
        rkxvector(k)=r1*G(1,1) 
        rkyvector(k)=0.0d0 
        rkzvector(k)=0.0d0
        k=k+1       
      end do 
    endif
  
  ! 2D
  ! Loop order convention: the interpolation routine get_fk_ex_k_interp
  ! uses rk1 as the FAST (unit-stride) index and rk2 as the SLOW index:
  !   ibz = nblock + int((rk2+0.5)/slice)*nside
  ! so the grid must be laid out with rk1 varying in the inner loop
  ! and rk2 in the outer loop.
  elseif ( ndim == 2 ) then
            
    if ( active_y .and. active_z ) then
      ! active axes: rk2 (fast) and rk3 (slow)
      k=1
      do i2=1,npointstotal_sq        ! rk3 slow
        r2=-0.5d0+(dble(i2-1)+offset)*step
        do i1=1,npointstotal_sq      ! rk2 fast
          r1=-0.5d0+(dble(i1-1)+offset)*step
          rkxvector(k)=0.0d0
          rkyvector(k)=r1*G(2,2)+r2*G(3,2) 
          rkzvector(k)=r1*G(2,3)+r2*G(3,3) 
          k=k+1       
        end do
      end do
      
    elseif ( active_x .and. active_z ) then
      ! active axes: rk1 (fast) and rk3 (slow)
      k=1
      do i2=1,npointstotal_sq        ! rk3 slow
        r2=-0.5d0+(dble(i2-1)+offset)*step
        do i1=1,npointstotal_sq      ! rk1 fast
          r1=-0.5d0+(dble(i1-1)+offset)*step
          rkxvector(k)=r1*G(1,1)+r2*G(3,1) 
          rkyvector(k)=0.0d0
          rkzvector(k)=r1*G(1,3)+r2*G(3,3) 
          k=k+1       
        end do
      end do
      
    else
      ! active axes: rk1 (fast) and rk2 (slow)
      k=1
      do i2=1,npointstotal_sq        ! rk2 slow
        r2=-0.5d0+(dble(i2-1)+offset)*step
        do i1=1,npointstotal_sq      ! rk1 fast
          r1=-0.5d0+(dble(i1-1)+offset)*step
          rkxvector(k)=r1*G(1,1)+r2*G(2,1) 
          rkyvector(k)=r1*G(1,2)+r2*G(2,2) 
          rkzvector(k)=0.0d0
          k=k+1       
        end do
      end do 
    endif
        
  ! 3D
  ! Interpolation convention:
  !   nblock = int((rk1+0.5)/slice)
  !          + int((rk2+0.5)/slice)*nside
  !          + int((rk3+0.5)/slice)*nside^2 + 1
  ! so rk1 is fastest, rk2 middle, rk3 slowest.
  ! Loop order: i3 outermost, i2 middle, i1 innermost.
  else  
    k=1
    do i3=1,npointstotal_sq          ! rk3 slow
      r3=-0.5d0+(dble(i3-1)+offset)*step
      do i2=1,npointstotal_sq        ! rk2 middle
        r2=-0.5d0+(dble(i2-1)+offset)*step
        do i1=1,npointstotal_sq      ! rk1 fast
          r1=-0.5d0+(dble(i1-1)+offset)*step
          rkxvector(k)=r1*G(1,1)+r2*G(2,1)+r3*G(3,1)
          rkyvector(k)=r1*G(1,2)+r2*G(2,2)+r3*G(3,2)
          rkzvector(k)=r1*G(1,3)+r2*G(2,3)+r3*G(3,3)
          k=k+1      
        end do
      end do
    end do 
    
  endif
  
end subroutine get_grid
 
 
end module parser_optics_xatu_dim
 
