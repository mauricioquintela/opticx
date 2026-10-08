!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reader of the Xatu exciton archive <label>.h5 written by 'xatu -H' (layout 'xatu-excitons', version 1).
! Reads only what opticx takes from the text .eigval/.states pair: the band window, the k mesh, the exciton
! energies and the resonant envelopes A_vck(n) (dataset 'state' of each /excitons/<NNNN> group). Values are
! returned in the archive's units (eV, Angstrom^-1); the caller converts them exactly as for the text files.
! The archive stores them at full double precision, where the text files carry 7-8 significant digits.
!
! The basis ordering is CHECKED, not assumed: row r (0-based) of /basis must be (valence band v_iv,
! conduction band c_ic, k index ik) with r = ik*nv*nc + ic*nv + iv, which is how fk_ex is indexed.
!
! Built only with 'make HDF5=1' (preprocessor flag OPTICX_HDF5); without it every entry point stops with a
! message naming the fix.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xatu_h5
#ifdef OPTICX_HDF5
  use hdf5
  use iso_c_binding
#endif
  use iso_fortran_env, only: int8, int64
  implicit none
  private

  public :: xatu_h5_info, xatu_h5_header, xatu_h5_energies, xatu_h5_kpoints, xatu_h5_states

  type :: xatu_h5_info
    integer :: nk = 0            ! k points of the BSE mesh
    integer :: ncells = 0        ! mesh points per reciprocal axis
    integer :: submesh = 1       ! Xatu submesh factor (1 = full mesh)
    integer :: nv = 0, nc = 0    ! valence and conduction bands in the exciton basis
    integer :: fermi_level = 0   ! absolute 0-based index of the highest valence band (Xatu convention)
    integer :: nexc = 0          ! exciton states in the archive
    integer :: dim = 0           ! exciton basis size, nk*nv*nc
    logical :: tda = .true.      ! Tamm-Dancoff approximation in the Xatu run
    integer, allocatable :: vbands(:), cbands(:)   ! absolute 0-based band numbers, in basis order
  end type xatu_h5_info

#ifdef OPTICX_HDF5
  character(len=*), parameter :: FORMAT_NAME = 'xatu-excitons'
  integer, parameter :: FORMAT_VERSION = 1

  ! Frees the buffer HDF5 allocates for a variable-length string (not wrapped by every Fortran build).
  interface
    integer(c_int) function h5free_memory_c(buf) bind(C, name='H5free_memory')
      import :: c_ptr, c_int
      type(c_ptr), value :: buf
    end function h5free_memory_c
  end interface
#endif

contains

#ifdef OPTICX_HDF5

!---------------------------------------------------------------------------------------------------------
! Header: format, band window, mesh and exciton count, plus the basis ordering check.
!---------------------------------------------------------------------------------------------------------
subroutine xatu_h5_header(path, info)
  character(len=*), intent(in) :: path
  type(xatu_h5_info), intent(out) :: info
  integer(hid_t) :: file
  integer(int64), allocatable :: basis(:,:), ivec(:)
  integer :: r, ik, ic, iv, np, err
  character(len=:), allocatable :: fmt

  call open_archive(path, file)

  fmt = read_attr_string(file, 'format', path)
  if (fmt /= FORMAT_NAME) then
    write(*,*) 'ERROR (xatu_h5): ', trim(path), ' is not a Xatu exciton archive (format "', fmt, &
               '", expected "', FORMAT_NAME, '").'
    stop 1
  end if
  if (read_attr_int(file, 'format_version', path) /= FORMAT_VERSION) then
    write(*,'(a,i0,a,i0,a)') ' ERROR (xatu_h5): '//trim(path)//' has format version ', &
         read_attr_int(file, 'format_version', path), '; this opticx reads version ', FORMAT_VERSION, '.'
    stop 1
  end if

  info%nk          = read_attr_int(file, 'n_kpoints', path)
  info%ncells      = read_attr_int(file, 'ncells', path)
  info%submesh     = read_attr_int(file, 'submesh_factor', path)
  info%fermi_level = read_attr_int(file, 'fermi_level', path)
  info%nexc        = read_attr_int(file, 'n_excitons', path)
  info%dim         = read_attr_int(file, 'exciton_basis_dim', path)
  info%tda         = read_attr_bool(file, 'tamm_dancoff', path)
  call read_attr_ints(file, 'valence_bands', path, ivec)
  info%vbands = int(ivec)
  call read_attr_ints(file, 'conduction_bands', path, ivec)
  info%cbands = int(ivec)
  info%nv = size(info%vbands)
  info%nc = size(info%cbands)
  np = info%nv*info%nc

  if (info%dim /= info%nk*np) then
    write(*,'(a,i0,a,i0,a,i0,a,i0)') ' ERROR (xatu_h5): inconsistent archive, exciton_basis_dim = ', info%dim, &
         ' but n_kpoints*nv*nc = ', info%nk, '*', info%nv, '*', info%nc
    stop 1
  end if

  ! The envelopes are only in the archive when Xatu ran with -c.
  if (.not. link_exists(file, 'basis') .or. .not. link_exists(file, 'excitons/'//state_name(1, info%nexc))) then
    write(*,*) 'ERROR (xatu_h5): ', trim(path), ' holds no exciton states. Rerun Xatu with -c -H.'
    stop 1
  end if

  allocate(basis(3, info%dim))
  call read_dataset_int64(file, 'basis', path, basis)
  do r = 0, info%dim - 1
    ik = r/np
    ic = mod(r, np)/info%nv
    iv = mod(r, info%nv)
    if (basis(1,r+1) /= info%vbands(iv+1) .or. basis(2,r+1) /= info%cbands(ic+1) .or. basis(3,r+1) /= ik) then
      write(*,'(a,i0,a,3i8,a,3i8)') ' ERROR (xatu_h5): exciton basis row ', r, ' is (v, c, k) =', basis(:,r+1), &
           ', expected', info%vbands(iv+1), info%cbands(ic+1), ik
      write(*,*) '       opticx indexes the envelopes as k slowest, then conduction, then valence fastest.'
      stop 1
    end if
  end do

  call h5fclose_f(file, err)
end subroutine xatu_h5_header

!---------------------------------------------------------------------------------------------------------
! The n lowest exciton energies (eV).
!---------------------------------------------------------------------------------------------------------
subroutine xatu_h5_energies(path, n, e)
  character(len=*), intent(in) :: path
  integer, intent(in) :: n
  real(8), intent(out) :: e(n)
  integer(hid_t) :: file
  real(8), allocatable, target :: eall(:)
  integer :: nexc, err

  call open_archive(path, file)
  nexc = read_attr_int(file, 'n_excitons', path)
  if (n > nexc) then
    write(*,'(a,i0,a,i0,a)') ' ERROR (xatu_h5): ', n, ' excitons requested but '//trim(path)//' holds ', nexc, '.'
    stop 1
  end if
  allocate(eall(nexc))
  call read_dataset_real(file, 'summary/energies', path, c_loc(eall), nexc)
  e = eall(1:n)
  call h5fclose_f(file, err)
end subroutine xatu_h5_energies

!---------------------------------------------------------------------------------------------------------
! The k mesh, k(1:3, ik) in Angstrom^-1 (/kpoints is (nk, 3) in C order, i.e. (3, nk) here).
!---------------------------------------------------------------------------------------------------------
subroutine xatu_h5_kpoints(path, nk, k)
  character(len=*), intent(in) :: path
  integer, intent(in) :: nk
  real(8), intent(out), target :: k(3, nk)
  integer(hid_t) :: file
  integer :: err

  call open_archive(path, file)
  call read_dataset_real(file, 'kpoints', path, c_loc(k), 3*nk)
  call h5fclose_f(file, err)
end subroutine xatu_h5_kpoints

!---------------------------------------------------------------------------------------------------------
! Resonant envelopes of the n lowest excitons: fk(:, i) = /excitons/<i>/state.
!---------------------------------------------------------------------------------------------------------
subroutine xatu_h5_states(path, dim, n, fk)
  character(len=*), intent(in) :: path
  integer, intent(in) :: dim, n
  complex(8), intent(out), target :: fk(dim, n)
  integer(hid_t) :: file, ctype, dset, space
  integer(hsize_t) :: npts
  integer :: i, nexc, err
  integer(size_t) :: dsize
  type(c_ptr) :: buf

  call open_archive(path, file)
  nexc = read_attr_int(file, 'n_excitons', path)
  if (n > nexc) then
    write(*,'(a,i0,a,i0,a)') ' ERROR (xatu_h5): ', n, ' excitons requested but '//trim(path)//' holds ', nexc, '.'
    stop 1
  end if

  ! complex(8) is (re, im) in memory: the archive's compound {r, i} of two doubles.
  dsize = 8
  call h5tcreate_f(H5T_COMPOUND_F, 2*dsize, ctype, err)
  call check(err, 'create the complex type', path)
  call h5tinsert_f(ctype, 'r', 0_size_t, H5T_NATIVE_DOUBLE, err)
  call check(err, 'create the complex type', path)
  call h5tinsert_f(ctype, 'i', dsize, H5T_NATIVE_DOUBLE, err)
  call check(err, 'create the complex type', path)

  do i = 1, n
    call h5dopen_f(file, 'excitons/'//state_name(i, nexc)//'/state', dset, err)
    call check(err, 'open excitons/'//state_name(i, nexc)//'/state', path)
    call h5dget_space_f(dset, space, err)
    call h5sget_simple_extent_npoints_f(space, npts, err)
    call h5sclose_f(space, err)
    if (npts /= dim) then
      write(*,'(a,i0,a,i0)') ' ERROR (xatu_h5): state '//state_name(i, nexc)//' has ', npts, &
           ' coefficients, expected ', dim
      stop 1
    end if
    buf = c_loc(fk(1, i))
    call h5dread_f(dset, ctype, buf, err)
    call check(err, 'read excitons/'//state_name(i, nexc)//'/state', path)
    call h5dclose_f(dset, err)
  end do

  call h5tclose_f(ctype, err)
  call h5fclose_f(file, err)
end subroutine xatu_h5_states

!---------------------------------------------------------------------------------------------------------
! Helpers
!---------------------------------------------------------------------------------------------------------
subroutine open_archive(path, file)
  character(len=*), intent(in) :: path
  integer(hid_t), intent(out) :: file
  integer :: err
  logical :: exists

  inquire(file=trim(path), exist=exists)
  if (.not. exists) then
    write(*,*) 'ERROR (xatu_h5): Xatu archive not found: ', trim(path)
    stop 1
  end if
  call h5open_f(err)
  call check(err, 'initialise the HDF5 library', path)
  call h5eset_auto_f(0, err)        ! failures are reported by check(), not by the HDF5 error stack
  call h5fopen_f(trim(path), H5F_ACC_RDONLY_F, file, err)
  call check(err, 'open the file (is it an HDF5 file?)', path)
end subroutine open_archive

subroutine check(err, what, path)
  integer, intent(in) :: err
  character(len=*), intent(in) :: what, path
  if (err < 0) then
    write(*,*) 'ERROR (xatu_h5): could not ', what, ' in ', trim(path)
    stop 1
  end if
end subroutine check

! Group name of state i: 1-based, zero padded to max(4, digits of nexc), as Xatu writes it.
function state_name(i, nexc) result(name)
  integer, intent(in) :: i, nexc
  character(len=:), allocatable :: name
  character(len=32) :: digits
  integer :: width
  write(digits, '(i0)') nexc
  width = max(4, len_trim(digits))
  write(digits, '(i0.'//itoa(width)//')') i
  name = trim(digits)
end function state_name

function itoa(i) result(s)
  integer, intent(in) :: i
  character(len=:), allocatable :: s
  character(len=16) :: buf
  write(buf, '(i0)') i
  s = trim(buf)
end function itoa

logical function link_exists(file, name)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name
  integer :: err, i
  ! h5lexists_f fails on a missing intermediate group, so test each level of the path.
  link_exists = .false.
  do i = 1, len(name)
    if (name(i:i) == '/') then
      call h5lexists_f(file, name(1:i-1), link_exists, err)
      if (err < 0 .or. .not. link_exists) then
        link_exists = .false.
        return
      end if
    end if
  end do
  call h5lexists_f(file, name, link_exists, err)
  if (err < 0) link_exists = .false.
end function link_exists

integer function read_attr_int(file, name, path)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name, path
  integer(int64), allocatable :: v(:)
  call read_attr_ints(file, name, path, v)
  read_attr_int = int(v(1))
end function read_attr_int

subroutine read_attr_ints(file, name, path, v)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name, path
  integer(int64), allocatable, target, intent(out) :: v(:)
  integer(hid_t) :: attr, space
  integer(hsize_t) :: npts
  integer :: err
  type(c_ptr) :: buf

  call open_attr(file, name, path, attr)
  call h5aget_space_f(attr, space, err)
  call h5sget_simple_extent_npoints_f(space, npts, err)
  call h5sclose_f(space, err)
  allocate(v(npts))
  buf = c_loc(v)
  call h5aread_f(attr, h5kind_to_type(int64, H5_INTEGER_KIND), buf, err)
  call check(err, 'read attribute '//name, path)
  call h5aclose_f(attr, err)
end subroutine read_attr_ints

! Booleans are the int8 enum {FALSE, TRUE}: read them with their own type into an int8.
logical function read_attr_bool(file, name, path)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name, path
  integer(hid_t) :: attr, ftype
  integer(int8), target :: v
  integer :: err
  type(c_ptr) :: buf

  call open_attr(file, name, path, attr)
  call h5aget_type_f(attr, ftype, err)
  buf = c_loc(v)
  call h5aread_f(attr, ftype, buf, err)
  call check(err, 'read attribute '//name, path)
  call h5tclose_f(ftype, err)
  call h5aclose_f(attr, err)
  read_attr_bool = (v /= 0)
end function read_attr_bool

! Variable-length UTF-8 string attribute.
function read_attr_string(file, name, path) result(s)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name, path
  character(len=:), allocatable :: s
  integer(hid_t) :: attr, ftype
  type(c_ptr), target :: cstr(1)
  character(kind=c_char), pointer :: chars(:)
  integer :: err, n
  type(c_ptr) :: buf

  call open_attr(file, name, path, attr)
  call h5aget_type_f(attr, ftype, err)
  buf = c_loc(cstr)
  call h5aread_f(attr, ftype, buf, err)
  call check(err, 'read attribute '//name, path)
  call c_f_pointer(cstr(1), chars, [4096])
  n = 0
  do while (n < 4096)
    if (chars(n+1) == c_null_char) exit
    n = n + 1
  end do
  allocate(character(len=n) :: s)
  s = transfer(chars(1:n), s)
  err = h5free_memory_c(cstr(1))
  call h5tclose_f(ftype, err)
  call h5aclose_f(attr, err)
end function read_attr_string

subroutine open_attr(file, name, path, attr)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name, path
  integer(hid_t), intent(out) :: attr
  logical :: exists
  integer :: err
  call h5aexists_f(file, name, exists, err)
  if (err < 0 .or. .not. exists) then
    write(*,*) 'ERROR (xatu_h5): attribute "', name, '" missing from ', trim(path), &
               ' (not written by a current Xatu -H?)'
    stop 1
  end if
  call h5aopen_f(file, name, attr, err)
  call check(err, 'open attribute '//name, path)
end subroutine open_attr

subroutine read_dataset_real(file, name, path, buf, n)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name, path
  type(c_ptr), intent(in) :: buf
  integer, intent(in) :: n
  integer(hid_t) :: dset, space
  integer(hsize_t) :: npts
  integer :: err
  type(c_ptr) :: p

  if (.not. link_exists(file, name)) then
    write(*,*) 'ERROR (xatu_h5): dataset ', name, ' missing from ', trim(path)
    stop 1
  end if
  call h5dopen_f(file, name, dset, err)
  call check(err, 'open dataset '//name, path)
  call h5dget_space_f(dset, space, err)
  call h5sget_simple_extent_npoints_f(space, npts, err)
  call h5sclose_f(space, err)
  if (npts /= n) then
    write(*,'(a,i0,a,i0)') ' ERROR (xatu_h5): dataset '//name//' has ', npts, ' values, expected ', n
    stop 1
  end if
  p = buf
  call h5dread_f(dset, H5T_NATIVE_DOUBLE, p, err)
  call check(err, 'read dataset '//name, path)
  call h5dclose_f(dset, err)
end subroutine read_dataset_real

subroutine read_dataset_int64(file, name, path, v)
  integer(hid_t), intent(in) :: file
  character(len=*), intent(in) :: name, path
  integer(int64), intent(out), target :: v(:,:)
  integer(hid_t) :: dset, space
  integer(hsize_t) :: npts
  integer :: err
  type(c_ptr) :: buf

  call h5dopen_f(file, name, dset, err)
  call check(err, 'open dataset '//name, path)
  call h5dget_space_f(dset, space, err)
  call h5sget_simple_extent_npoints_f(space, npts, err)
  call h5sclose_f(space, err)
  if (npts /= size(v)) then
    write(*,'(a,i0,a,i0)') ' ERROR (xatu_h5): dataset '//name//' has ', npts, ' values, expected ', size(v)
    stop 1
  end if
  buf = c_loc(v)
  call h5dread_f(dset, h5kind_to_type(int64, H5_INTEGER_KIND), buf, err)
  call check(err, 'read dataset '//name, path)
  call h5dclose_f(dset, err)
end subroutine read_dataset_int64

#else

subroutine xatu_h5_header(path, info)
  character(len=*), intent(in) :: path
  type(xatu_h5_info), intent(out) :: info
  call no_hdf5(path)
end subroutine xatu_h5_header

subroutine xatu_h5_energies(path, n, e)
  character(len=*), intent(in) :: path
  integer, intent(in) :: n
  real(8), intent(out) :: e(n)
  e = 0.0d0
  call no_hdf5(path)
end subroutine xatu_h5_energies

subroutine xatu_h5_kpoints(path, nk, k)
  character(len=*), intent(in) :: path
  integer, intent(in) :: nk
  real(8), intent(out) :: k(3, nk)
  k = 0.0d0
  call no_hdf5(path)
end subroutine xatu_h5_kpoints

subroutine xatu_h5_states(path, dim, n, fk)
  character(len=*), intent(in) :: path
  integer, intent(in) :: dim, n
  complex(8), intent(out) :: fk(dim, n)
  fk = (0.0d0, 0.0d0)
  call no_hdf5(path)
end subroutine xatu_h5_states

subroutine no_hdf5(path)
  character(len=*), intent(in) :: path
  write(*,*) 'ERROR: ', trim(path), ' is a Xatu HDF5 archive, but this opticx was built without HDF5.'
  write(*,*) '       Rebuild with "make clean && make HDF5=1", or give the .eigval and .states text files.'
  stop 1
end subroutine no_hdf5

#endif

end module xatu_h5
