module parser_input_file
  implicit none
  private
  public :: material_name_in
  public :: xatu_eigval_filepath_in
  public :: xatu_states_filepath_in
  public :: iflag_ome_sp_text
  public :: iflag_ome_ex_text
  public :: response_text
  ! Second-order two-frequency control (2026-09-24).
  ! sigma^{abc}(w_p + w_q; w_p, w_q). Two ways to say what w_q is:
  !   Frequency_ratio r      -> w_q = r * w_p, scanning w_p over Energy_variables.
  !                             r = 1 is SHG, r = 0 electro-optic, r = -1 optical rectification.
  !   Energy_variables_2     -> an INDEPENDENT w_q grid (e1b e2b nwb); its presence selects the full
  !                             2D (w_p, w_q) map and overrides Frequency_ratio.
  public :: freq_ratio, e1b, e2b, nwb, two_freq_grid
  public :: cache_ome_read, cache_ome_write
  public :: xnm_covariant, basis_repair
  public :: sp_covariant
  public :: build_freq_pairs
  public :: iflag_xatu
  public :: iflag_xatu_h5, xatu_h5_filepath_in
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
  public :: iflag_orthonormal_text
  public :: iflag_orthonormal
  public :: kpath_nv, kpath_frac, kpath_count, kpath_labels, kpath_has_labels
  public :: second_order_response

  character(len=1000) :: material_name_in
  character(len=100) :: filename_input
  character(len=100) :: iflag_xatu_text
  character(len=100) :: iflag_ome_sp_text
  character(len=100) :: iflag_ome_ex_text
  character(len=100) :: broadening_type_text
  character(len=100) :: iflag_write_exk_text
  character(len=1000) :: xatu_eigval_filepath_in
  character(len=1000) :: xatu_states_filepath_in
  ! Xatu exciton archive (xatu -H): one '<label>.h5' path in place of the .eigval/.states pair
  character(len=1000) :: xatu_h5_filepath_in = ''
  character(len=100) :: response_text
  character(len=100) :: iflag_orthonormal_text
  ! Band-structure path (Kpath): vertices in reduced coordinates along the reciprocal lattice vectors and,
  ! for each, the number of points to the next vertex (1 = this vertex alone, then a jump; the last vertex
  ! carries 1 to be included). kpath_nv = 0: no path given, bands.f90 uses its default for the lattice.
  integer :: kpath_nv = 0
  real(8), allocatable :: kpath_frac(:,:)
  integer, allocatable :: kpath_count(:)
  character(len=16), allocatable :: kpath_labels(:)
  logical :: kpath_has_labels = .false.
  real(8) :: freq_ratio = 1.0d0
  ! Opt-in cache of the second-order excitonic OMEs. Both default .false.: the cache
  ! header fingerprints the exciton solution but NOT the Wannier90 model, so reuse is a choice.
  ! Read and write are INDEPENDENT -- a large run may be worth reading back but too big to store
  ! (the payload is 6*N^2 complex(8): 338 MB at N = 1875, 3.0 GB at N = 5625, 9.6 GB at N = 10000).
  logical :: cache_ome_read  = .false.
  logical :: cache_ome_write = .false.
  character(len=100) :: cache_ome_ex_text = 'false'
  ! How the intraband part of the inter-exciton position X_nm is evaluated.
  !   covariant          (default) discrete covariant derivative of the exciton envelopes, transported
  !                      between grid neighbours with the window overlap matrices <u_n(k)|u_m(k+b)>;
  !                      invariant under any per-k unitary rotation of the window bands.
  !   finite_difference  the original form: plain central difference of the envelopes plus the
  !                      Berry-connection diagonal and r_nm = -i v_nm/(E_n-E_m) inside the window.
  !                      Correct only where the band gauge is smooth between grid neighbours.
  logical :: xnm_covariant = .true.
  ! Repair of the exciton-envelope basis at k-points where window bands are EXACTLY degenerate
  !. Xatu writes the envelopes in whatever basis its diagonaliser picked inside a
  ! degenerate block, which opticx cannot know; the repair chooses, per block, the unitary that makes
  ! the envelopes consistent with their covariantly transported neighbours (needs Xnm_derivative =
  ! covariant, second order). Keyword Exciton_basis_repair = true | false.
  logical :: basis_repair = .true.
  ! Ex_rectification: OBSOLETE. The excitonic rectification is always the whole causal sigma(0; w, -w); the
  ! shift current in the convention of npj Comput. Mater. 11, 13 is Response = shift. 'causal' is accepted
  ! with a note, 'shift' stops (it would otherwise silently return a different quantity).
  character(len=100) :: ex_rect_text = 'causal'
  character(len=100) :: basis_repair_text = 'true'
  character(len=100) :: xnm_derivative_text = 'covariant'
  ! Single-particle second-order method, keyword Sp_method = covariant | per_band. covariant (DEFAULT):
  ! block-covariant generalised derivative, no degeneracy cut -- Response = shift -> shift_covariant, and shg,
  ! electrooptic, rectification, general use the covariant method-B kernel. per_band: the older per-band routes
  ! (shift -> shift_shiftvector; Taghizadeh 2017 Eq. A3a for the others), kept for comparison.
  logical :: sp_covariant = .true.
  character(len=100) :: sp_method_text = 'covariant'
  real(8) :: e1b = 0.0d0, e2b = 0.0d0
  integer :: nwb = 0
  logical :: two_freq_grid = .false.

  logical :: iflag_xatu
  logical :: iflag_xatu_h5 = .false.
  logical :: iflag_ome_sp
  logical :: iflag_ome_ex
  logical :: iflag_orthonormal = .true.   ! default for any caller that does not parse an input file
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
  !> True for the Response values evaluated by the second-order routines (optical_response::second_order).
  ! A Xatu exciton archive is recognised by its extension, .h5 or .hdf5 (any case).
  logical function is_h5_path(path)
    character(len=*), intent(in) :: path
    character(len=len(path)) :: low
    integer :: i, n
    low = path
    do i = 1, len(low)
      if (low(i:i) >= 'A' .and. low(i:i) <= 'Z') low(i:i) = achar(iachar(low(i:i)) + 32)
    end do
    n = len_trim(low)
    is_h5_path = .false.
    if (n >= 3) is_h5_path = low(n-2:n) == '.h5'
    if (n >= 5) is_h5_path = is_h5_path .or. low(n-4:n) == '.hdf5'
  end function is_h5_path

  logical function second_order_response()
    select case (trim(response_text))
      case ('shift_sumrule', 'shift_shiftvector', 'shift_gender', 'shift_covariant', 'shg', 'shg_covariant', &
            'electrooptic', 'rectification', 'general')
        second_order_response = .true.
      case default
        second_order_response = .false.
    end select
  end function second_order_response

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
      logical :: response_found, energy_found, exciton_found, iflag_orthonormal_found
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
      norb_ex_cut = 0               ! Exciton_cutoff absent: every exciton in the Xatu files (get_exciton_dim)
      iflag_orthonormal_found = .false.
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
              ! Read the file paths that follow: either one Xatu HDF5 archive (.h5/.hdf5, written by
              ! 'xatu -H') or the .eigval and .states text files, in that order.
              read(iounit10,'(A)') xatu_eigval_filepath_in
              xatu_eigval_filepath_in = adjustl(xatu_eigval_filepath_in)
              if (is_h5_path(xatu_eigval_filepath_in)) then
                iflag_xatu_h5 = .true.
                xatu_h5_filepath_in = xatu_eigval_filepath_in
                xatu_eigval_filepath_in = ''
              else
                read(iounit10,'(A)') xatu_states_filepath_in
              end if
            else if (iflag_xatu_text == 'false') then
              iflag_xatu = .false.
            else
              write(*,*) 'Error: Invalid value in Xatu_interface. Expected "true" or "false".'
              stop 1
            end if
            
          else if (index(param_name, 'Exciton_cutoff') > 0) then
            read(iounit10,*) norb_ex_cut
            exciton_found = .true.
            if (norb_ex_cut < 0) then
              write(*,*) 'Error: Exciton_cutoff must be positive (omit it to use every exciton in the Xatu files).'
              stop 1
            end if
            
          else if (index(param_name, 'Bandlist') > 0) then
            call read_line_numbers_int(iounit10, narray, num_values)
            bandlist_found = .true.

          else if (index(param_name, 'Kpath_labels') > 0) then
            call read_kpath_labels(iounit10)

          else if (index(param_name, 'Kpath') > 0) then
            call read_kpath(iounit10)
            
          else if (index(param_name, 'Ncells') > 0) then
            read(iounit10,*) npointstotal_sq
            ncells_found = .true.
            
          else if (index(param_name, 'Nfermi') > 0) then
            read(iounit10,*) nf
            nfermi_found = .true.
            
          else if (index(param_name, 'OME_sp') > 0 .or. index(param_name, 'OME_SP') > 0) then
            read(iounit10,*) iflag_ome_sp_text
            ome_sp_found = .true.
            
          else if (index(param_name, 'Orthonormal') > 0) then
            read(iounit10,*) iflag_orthonormal_text   ! iounit10 (newunit): unit 10 is not connected
            iflag_orthonormal_found = .true. 
            
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
          else if (index(param_name, 'Ex_rectification') > 0) then
            read(iounit10,'(A)') ex_rect_text
            ex_rect_text = to_lower(adjustl(ex_rect_text))
            select case (trim(ex_rect_text))
              case ('causal')
                write(*,*) 'NOTE (parser_input_file): Ex_rectification is obsolete; the excitonic rectification'
                write(*,*) '     is always the causal sigma(0; w, -w). The keyword is ignored.'
              case ('shift')
                write(*,*) 'ERROR (parser_input_file): Ex_rectification = shift has been removed. The excitonic'
                write(*,*) '      rectification is now the whole causal sigma(0; w, -w), injection current included;'
                write(*,*) '      for the shift current (npj Comput. Mater. 11, 13 convention) use Response = shift.'
                stop 1
              case default
                write(*,*) 'ERROR (parser_input_file): Ex_rectification = "'// &
                           trim(ex_rect_text)//'" is not recognised (and the keyword is obsolete).'
                stop 1
            end select
          else if (index(param_name, 'Exciton_basis_repair') > 0) then
            read(iounit10,'(A)') basis_repair_text
            basis_repair_text = to_lower(adjustl(basis_repair_text))
            select case (trim(basis_repair_text))
              case ('true', '.true.', 'on', 'yes')
                basis_repair = .true.
              case ('false', '.false.', 'off', 'no')
                basis_repair = .false.
              case default
                write(*,*) 'ERROR (parser_input_file): Exciton_basis_repair = "'// &
                           trim(basis_repair_text)//'" is not recognised. Valid: true, false.'
                stop 1
            end select
          else if (index(param_name, 'Xnm_derivative') > 0) then
            read(iounit10,'(A)') xnm_derivative_text
            xnm_derivative_text = to_lower(adjustl(xnm_derivative_text))
            select case (trim(xnm_derivative_text))
              case ('covariant')
                xnm_covariant = .true.
              case ('finite_difference', 'plain')
                xnm_covariant = .false.
              case default
                write(*,*) 'ERROR (parser_input_file): Xnm_derivative = "'// &
                           trim(xnm_derivative_text)//'" is not recognised.'
                write(*,*) '       Valid: covariant (default), finite_difference (= plain).'
                stop 1
            end select
          else if (index(param_name, 'Sp_method') > 0) then
            read(iounit10,'(A)') sp_method_text
            sp_method_text = to_lower(adjustl(sp_method_text))
            select case (trim(sp_method_text))
              case ('covariant')
                sp_covariant = .true.
              case ('per_band')
                sp_covariant = .false.
              case default
                write(*,*) 'ERROR (parser_input_file): Sp_method = "'//trim(sp_method_text)//'" is not recognised.'
                write(*,*) '       Valid: covariant (default), per_band.'
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
      
      
      if (iflag_orthonormal_found .and. iflag_orthonormal_text == 'false') then
        iflag_orthonormal = .false.
      else if (iflag_orthonormal_found .and. iflag_orthonormal_text == 'true') then
        iflag_orthonormal = .true.
      else if (iflag_orthonormal_found) then
              write(*,*) 'ERROR: Invalid value in Orthonormal. Expected "true" or "false".'
              stop 1
      else if (.not. iflag_orthonormal_found) then
        iflag_orthonormal = .true.
      end if
      
      if (iflag_write_exk_text == 'true') then
        iflag_write_exk = .true.
      else
        iflag_write_exk = .false.
      end if

      
      ! Response = shift: the recommended single-particle shift current, chosen by Sp_method.
      ! The explicit names (shift_covariant, shift_shiftvector, ...) keep selecting exactly what they name;
      ! shg_covariant always selects the covariant SHG.
      if (trim(response_text) == 'shift') then
        if (sp_covariant) then
          response_text = 'shift_covariant'
        else
          response_text = 'shift_shiftvector'
        end if
      end if
      if (trim(response_text) == 'shg_covariant') sp_covariant = .true.

      write(*,*) '   Input file has been read'
    end subroutine get_input_file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Read the Kpath block: one vertex per line, "f1 f2 f3 n" (reduced coordinates along the reciprocal
    !! lattice vectors, n = number of points from this vertex to the next; n = 1 puts the vertex alone and
    !! jumps to the next one; the last vertex carries 1 to be included, 0 to be left out). The block ends
    !! at a blank line, the next '#' line (pushed back) or the end of the file.
    !! @param iounit  The open input file.
    !! @return void
    subroutine read_kpath(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer, parameter :: maxv = 1000
    real(8) :: f(3, maxv)
    integer :: n(maxv), ios, nv
    character(len=1000) :: line
    nv = 0
    do
      read(iounit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) exit
      if (line(1:1) == '#') then
        backspace(iounit)
        exit
      end if
      if (nv == maxv) then
        write(*,*) 'ERROR (parser_input_file): Kpath has more than', maxv, 'vertices.'
        stop 1
      end if
      nv = nv + 1
      read(line,*,iostat=ios) f(1,nv), f(2,nv), f(3,nv), n(nv)
      if (ios /= 0) then
        write(*,*) 'ERROR (parser_input_file): Kpath line "'//trim(line)//'" is not "f1 f2 f3 npoints".'
        stop 1
      end if
    end do
    if (nv < 2) then
      write(*,*) 'ERROR (parser_input_file): Kpath needs at least two vertices.'
      stop 1
    end if
    if (any(n(1:nv-1) < 1) .or. n(nv) < 0) then
      write(*,*) 'ERROR (parser_input_file): Kpath point counts must be >= 1 (1 = jump to the next vertex);'
      write(*,*) '       the last vertex takes 1 (included) or 0 (left out).'
      stop 1
    end if
    kpath_nv = nv
    allocate(kpath_frac(3,nv), kpath_count(nv))
    kpath_frac = f(:,1:nv); kpath_count = n(1:nv)
    end subroutine read_kpath
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Read the Kpath_labels line: one label per Kpath vertex, separated by spaces (e.g. "G M K G").
    !! @param iounit  The open input file.
    !! @return void
    subroutine read_kpath_labels(iounit)
    implicit none
    integer, intent(in) :: iounit
    character(len=1000) :: line
    character(len=16) :: tok(1000)
    integer :: ios, i, nt, istart
    read(iounit,'(A)',iostat=ios) line
    if (ios /= 0) return
    nt = 0; i = 1
    line = adjustl(line)
    do while (i <= len_trim(line))
      if (line(i:i) == ' ') then
        i = i + 1; cycle
      end if
      istart = i
      do while (i <= len_trim(line))
        if (line(i:i) == ' ') exit
        i = i + 1
      end do
      nt = nt + 1; tok(nt) = line(istart:i-1)
    end do
    if (nt > 0) then
      if (allocated(kpath_labels)) deallocate(kpath_labels)
      allocate(kpath_labels(nt)); kpath_labels = tok(1:nt); kpath_has_labels = .true.
    end if
    end subroutine read_kpath_labels
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
    stop 1
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
  ! Two modes:
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