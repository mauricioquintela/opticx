module parser_input_file
  implicit none
  private
  public :: material_name_in
  public :: xatu_eigval_filepath_in
  public :: xatu_states_filepath_in
  public :: iflag_ome_sp_text
  public :: iflag_ome_ex_text
  public :: response_text
  ! Second-order two-frequency control (2026-09-24, HANDOFF 8.36).
  ! sigma^{abc}(w_p + w_q; w_p, w_q). Two ways to say what w_q is:
  !   Frequency_ratio r      -> w_q = r * w_p, scanning w_p over Energy_variables.
  !                             r = 1 is SHG, r = 0 electro-optic, r = -1 optical rectification.
  !   Energy_variables_2     -> an INDEPENDENT w_q grid (e1b e2b nwb); its presence selects the full
  !                             2D (w_p, w_q) map and overrides Frequency_ratio.
  public :: freq_ratio, e1b, e2b, nwb, two_freq_grid
  public :: cache_ome_read, cache_ome_write
  public :: build_freq_pairs
  public :: iflag_xatu
  public :: iflag_ome_sp
  public :: iflag_ome_ex
  public :: ndim,nf,npointstotal_sq
  public :: e1,e2,eta,nw
  public :: get_input_file
  public :: nband_index
  public :: norb_ex_cut
  public :: broadening_type_text
  public :: iflag_write_exk
  public :: read_line_numbers_int !subroutine

  character(len=1000) :: material_name_in
  character(len=100) :: filename_input
  character(len=100) :: iflag_xatu_text
  character(len=100) :: iflag_ome_sp_text
  character(len=100) :: iflag_ome_ex_text
  character(len=100) :: broadening_type_text
  character(len=100) :: iflag_write_exk_text
  character(len=1000) :: xatu_eigval_filepath_in
  character(len=1000) :: xatu_states_filepath_in
  character(len=100) :: response_text
  real(8) :: freq_ratio = 1.0d0
  ! Opt-in cache of the second-order excitonic OMEs (HANDOFF 8.44). Both default .false.: the cache
  ! header fingerprints the exciton solution but NOT the Wannier90 model, so reuse is a choice.
  ! Read and write are INDEPENDENT -- a large run may be worth reading back but too big to store
  ! (the payload is 6*N^2 complex(8): 338 MB at N = 1875, 3.0 GB at N = 5625, 9.6 GB at N = 10000).
  logical :: cache_ome_read  = .false.
  logical :: cache_ome_write = .false.
  character(len=100) :: cache_ome_ex_text = 'false'
  real(8) :: e1b = 0.0d0, e2b = 0.0d0
  integer :: nwb = 0
  logical :: two_freq_grid = .false.

  logical :: iflag_xatu
  logical :: iflag_ome_sp
  logical :: iflag_ome_ex
  logical :: iflag_write_exk

  integer :: ndim
  integer :: nf
  integer :: npointstotal_sq
  integer :: norb_ex_cut
  integer :: nband_index
  integer :: nw
  real(8) :: e1,e2,eta

  allocatable :: nband_index(:)

  contains
    function to_lower(str) result(lower_str)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str)) :: lower_str
      integer :: i, ic
      
      do i = 1, len(str)
        ic = iachar(str(i:i))
        if (ic >= iachar('A') .and. ic <= iachar('Z')) then
          lower_str(i:i) = achar(ic + 32)
        else
          lower_str(i:i) = str(i:i)
        end if
      end do
    end function to_lower
    
    subroutine get_input_file()
      implicit none
      integer :: iounit10
      integer, allocatable :: narray(:) 
      integer :: num_values, ios
      character(len=1000) :: line
      character(len=100) :: param_name
      logical :: ndim_found, material_found, xatu_found, bandlist_found
      logical :: ncells_found, nfermi_found, ome_sp_found, ome_ex_found
      logical :: response_found, energy_found, exciton_found
      logical :: write_exk_found   ! PATCH: was declared-but-unused in a comment
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '1. Entering parser_input_file'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! Initialize flags
      ndim_found = .false.
      material_found = .false.
      xatu_found = .false.
      bandlist_found = .false.
      ncells_found = .false.
      nfermi_found = .false.
      ome_sp_found = .false.
      ome_ex_found = .false.
      response_found = .false.
      energy_found = .false.
      exciton_found = .false.
      write_exk_found = .false.
      ! default broadening
      !broadening_type_text = 'gaussian'
      broadening_type_text = 'lorentzian'
      iflag_write_exk_text = 'false'
      
      call get_command_argument(1,filename_input)
      open(newunit=iounit10,file=adjustl(filename_input))
      
      ! Read file sequentially and process parameters based on their labels
      do
        read(iounit10,'(A)',iostat=ios) line
        if (ios /= 0) exit  ! End of file
        
        line = adjustl(line)
        
        ! Check if this is a comment line (parameter label)
        if (line(1:1) == '#') then
          ! Extract parameter name and read corresponding value
          param_name = adjustl(line(3:))  ! Remove "# " prefix
          
          if (index(param_name, 'Periodic dimensions') > 0) then
            read(iounit10,*) ndim
            ndim_found = .true.
            
          else if (index(param_name, 'Wannier90_filename') > 0) then
            read(iounit10,'(A)') material_name_in
            material_found = .true.
            
          else if (index(param_name, 'Xatu_interface') > 0) then
            read(iounit10,*) iflag_xatu_text
            xatu_found = .true.
            
            if (iflag_xatu_text == 'true') then
              iflag_xatu = .true.
              ! Read the eigval and states file paths that follow
              read(iounit10,'(A)') xatu_eigval_filepath_in
              read(iounit10,'(A)') xatu_states_filepath_in
            else if (iflag_xatu_text == 'false') then
              iflag_xatu = .false.
            else
              write(*,*) 'Error: Invalid value in Xatu_interface. Expected "true" or "false".'
              stop
            end if
            
          else if (index(param_name, 'Exciton_cutoff') > 0) then
            read(iounit10,*) norb_ex_cut
            exciton_found = .true.
            
          else if (index(param_name, 'Bandlist') > 0) then
            call read_line_numbers_int(iounit10, narray, num_values)
            bandlist_found = .true.
            
          else if (index(param_name, 'Ncells') > 0) then
            read(iounit10,*) npointstotal_sq
            ncells_found = .true.
            
          else if (index(param_name, 'Nfermi') > 0) then
            read(iounit10,*) nf
            nfermi_found = .true.
            
          else if (index(param_name, 'OME_sp') > 0 .or. index(param_name, 'OME_SP') > 0) then
            read(iounit10,*) iflag_ome_sp_text
            ome_sp_found = .true.
            
          else if (index(param_name, 'OME_ex') > 0 .or. index(param_name, 'OME_EX') > 0) then
            read(iounit10,*) iflag_ome_ex_text
            ome_ex_found = .true.
            
          else if (index(param_name, 'Response') > 0) then
            read(iounit10,*) response_text
            response_found = .true.
            
          else if (index(param_name, 'Write_ex_kresolved') > 0) then
            read(iounit10,*) iflag_write_exk_text
            write_exk_found = .true.
            
          else if (index(param_name, 'Energy_variables_2') > 0) then
            read(iounit10,*) e1b, e2b, nwb
            two_freq_grid = .true.
          else if (index(param_name, 'Cache_ome_ex') > 0) then
            read(iounit10,'(A)') cache_ome_ex_text
            cache_ome_ex_text = to_lower(adjustl(cache_ome_ex_text))
            select case (trim(cache_ome_ex_text))
              case ('read')
                cache_ome_read = .true.;  cache_ome_write = .false.
              case ('write')
                cache_ome_read = .false.; cache_ome_write = .true.
              case ('true', 'readwrite', 'both', '.true.')
                cache_ome_read = .true.;  cache_ome_write = .true.
              case ('false', 'off', 'none', '.false.')
                cache_ome_read = .false.; cache_ome_write = .false.
              case default
                write(*,*) 'ERROR (parser_input_file): Cache_ome_ex = "'// &
                           trim(cache_ome_ex_text)//'" is not recognised.'
                write(*,*) '       Valid: read, write, readwrite (= true, both), off (= false, none).'
                stop 1
            end select
          else if (index(param_name, 'Frequency_ratio') > 0) then
            read(iounit10,*) freq_ratio
          else if (index(param_name, 'Energy_variables') > 0) then
            read(iounit10,*) e1, e2, eta, nw
            energy_found = .true.
          else if (index(param_name, 'Broadening_type') > 0 .or. index(param_name,'Broadening')>0) then
            read(iounit10,'(A)') broadening_type_text
            broadening_type_text = adjustl(broadening_type_text)
            broadening_type_text = to_lower(broadening_type_text)
          
          end if
        end if
      end do
      
      close(iounit10)
      
      ! Handle bandlist case: allocate nband_index if bandlist was found
      if (bandlist_found) then
        allocate(nband_index(num_values))
        ! narray may be over-allocated (its size, ncount, is estimated by counting spaces in the
      ! raw line, which over-counts on a leading space or a doubled inter-number space); the
      ! successfully-parsed tokens are always packed at the front, narray(1:num_values), so copy
      ! only that slice -- copying the whole (possibly larger) array crashes under -fcheck=all
      ! ("Array bound mismatch") or silently truncates the Bandlist without it.
      nband_index(:) = narray(1:num_values)
      end if
      
      ! Set npointstotal_sq to 0 if using xatu interface
      if (iflag_xatu) then
        npointstotal_sq = 0
      end if
      
      ! Declare flags from text strings
      if (iflag_ome_sp_text == 'true') then
        iflag_ome_sp = .true.
      else
        iflag_ome_sp = .false.
      end if
      
      if (iflag_ome_ex_text == 'true') then
        iflag_ome_ex = .true.
      else
        iflag_ome_ex = .false.
      end if
      
      if (iflag_write_exk_text == 'true') then
        iflag_write_exk = .true.
      else
        iflag_write_exk = .false.
      end if

      
      write(*,*) '   Input file has been read'
    end subroutine get_input_file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This routine reads a line of numbers into an array
    subroutine read_line_numbers_int(iounit,narray,num_values)
    implicit none
    allocatable :: narray(:)

    integer, intent(in) :: iounit   ! the caller's open unit; was a bare 10 (audit 2026-09-30)
    integer :: ios,ncount,i
    integer :: narray
    integer :: num_values
    integer :: temp_num 
    integer :: istart,iposition
    character(len=1000) :: line
 
    !Read the line of text
    read(iounit,'(A)', iostat=ios) line
    if (ios /= 0) then
      print *, 'Error reading file'
    stop
    end if

    ncount=0
    do i=1,len_trim(line)
      if (line(i:i) == ' ') ncount=ncount+1
    end do
    ncount=ncount+1  ! One more than the number of spaces

    !Allocate the array based on the number of values
    allocate(narray(ncount))

    !Reset the number of values counter
    num_values=0
    istart=1
    !Now, sequentially extract numbers from the line
    do i=1,ncount
      !Find the next number in the line
      read(line(istart:),*,iostat=ios) temp_num
      if (ios == 0) then
        num_values=num_values+1
        narray(num_values)=temp_num
        !Move the starting  position to the next number
        iposition=scan(line(istart:),' ')
        if (iposition>0) then
          istart=istart+iposition
        end if
      end if
    end do

    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Builds the list of (omega_p, omega_q) pairs for a general second-order run, in Hartree.
  ! Two modes (HANDOFF 8.36):
  !   Energy_variables_2 present -> the full 2D grid, nfreq = nw*nwb, omega_p the SLOW index.
  !   otherwise                  -> omega_q = freq_ratio * omega_p over the Energy_variables grid.
  ! Shared by the single-particle and excitonic drivers so the two can never drift apart.
  subroutine build_freq_pairs(nfreq, wpg, wqg)
    implicit none
    integer,              intent(out) :: nfreq
    real(8), allocatable, intent(out) :: wpg(:), wqg(:)
    integer :: iw, iwb, idx
    real(8) :: wrange, wrangeb, wa, wb_
    if (two_freq_grid) then
      nfreq = nw*nwb
    else
      nfreq = nw
    end if
    allocate(wpg(nfreq), wqg(nfreq))
    wrange = e2 - e1
    if (two_freq_grid) then
      wrangeb = e2b - e1b
      idx = 0
      do iw = 1, nw
        wa = (e1 + wrange/dble(nw)*dble(iw-1))/27.211385d0
        do iwb = 1, nwb
          wb_ = (e1b + wrangeb/dble(nwb)*dble(iwb-1))/27.211385d0
          idx = idx + 1
          wpg(idx) = wa
          wqg(idx) = wb_
        end do
      end do
    else
      do iw = 1, nw
        wpg(iw) = (e1 + wrange/dble(nw)*dble(iw-1))/27.211385d0
        wqg(iw) = freq_ratio*wpg(iw)
      end do
    end if
  end subroutine build_freq_pairs

end module parser_input_file