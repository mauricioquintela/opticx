module ome
  use parser_input_file, &
  only:iflag_ome_sp_text,iflag_ome_ex_text,cache_ome_read,second_order_response,iflag_xatu
  use ome_sp, &
  only:get_ome_sp, rotate_fk_ex_to_a4_basis, load_omesp_basis
  use ome_ex, &
  only:get_ome_ex, get_ome_ex_from_cache
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Map an OME_sp / OME_ex keyword to the order the corresponding get_ome_* routine wants.
! Returns 0 for 'none', meaning "read it from file instead, do not compute".
  integer function ome_order(keyword)
    implicit none
    character(len=*), intent(in) :: keyword
    select case (trim(keyword))
      case ('linear');    ome_order = 1
      case ('nonlinear'); ome_order = 2
      case default;       ome_order = 0
    end select
  end function ome_order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome()
    implicit none
    integer :: norder_sp, norder_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '4. Entering ome'
    !here we compute and write into files optical matrix elements
    norder_sp = ome_order(iflag_ome_sp_text)
    norder_ex = ome_order(iflag_ome_ex_text)

    !single-particle matrix elements
    if (norder_sp > 0) then
      call get_ome_sp(norder_sp)
    else
      write(*,*) '   Optical matrix elements (sp) will be read from file'
      ! The basis the exciton envelopes must be carried into is stored in the .omesp (since 2026-10-06)
      if (norder_ex > 0) call load_omesp_basis(norder_ex)
    end if

    ! Carry Xatu's exciton envelopes into the same basis the Eq. (A4) rotation put the
    ! single-particle states in, BEFORE any excitonic quantity consumes fk_ex. Without this the
    ! two are in different bases inside every near-degenerate multiplet.
    if (norder_ex > 0) call rotate_fk_ex_to_a4_basis()

    !excitonic matrix elements
    if (norder_ex > 0) then
      call get_ome_ex(norder_ex)
    else if (iflag_xatu .and. second_order_response()) then
      ! OME_ex = none on a second-order run: the excitonic elements come from the Cache_ome_ex file, as the
      ! single-particle ones come from the .omesp (the inter-exciton elements are stored nowhere else)
      if (.not. cache_ome_read) then
        write(*,*) 'ERROR (ome): OME_ex = none on a second-order run reads the excitonic matrix elements from the'
        write(*,*) '       Cache_ome_ex file: set Cache_ome_ex = read (or readwrite), or OME_ex = nonlinear.'
        stop 1
      end if
      call get_ome_ex_from_cache()
    else
      write(*,*) '   Optical matrix elements (ex) will be read from file'
    end if

  end subroutine get_ome
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ome
