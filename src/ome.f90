module ome
  use parser_input_file, &
  only:iflag_xatu, &
  iflag_ome_sp_text,iflag_ome_ex_text
  use ome_sp, &
  only:get_ome_sp, rotate_fk_ex_to_a4_basis
  use ome_ex, &
  only:get_ome_ex
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome()
    implicit none
    integer iflag_norder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !write(*,*) iflag_ome_sp
    !write(*,*) iflag_ome_ex
    write(*,*) '4. Entering ome'
    !here we compute and write into files optical matrix elements
    !single-particle matrix elements
    if (iflag_ome_sp_text == 'linear') then
      iflag_norder=1
      call get_ome_sp(iflag_norder)     
    end if
    if (iflag_ome_sp_text == 'nonlinear') then
      iflag_norder=2
      call get_ome_sp(iflag_norder)  
    end if
    if (iflag_ome_sp_text == 'none' ) then
      write(*,*) '   Optical matrix elements (sp) will be read from file'   
    end if 
    ! Carry Xatu's exciton envelopes into the same basis the Eq. (A4) rotation put the
    ! single-particle states in, BEFORE any excitonic quantity consumes fk_ex. Without this the
    ! two are in different bases inside every near-degenerate multiplet (HANDOFF 8.29).
    if (iflag_ome_ex_text == 'linear' .or. iflag_ome_ex_text == 'nonlinear') then
      call rotate_fk_ex_to_a4_basis()
    end if

    !excitonic matrix elements
    if (iflag_ome_ex_text == 'linear') then
      iflag_norder=1
      call get_ome_ex(iflag_norder)       
    end if
    if (iflag_ome_ex_text == 'nonlinear' ) then
      iflag_norder=2
      call get_ome_ex(iflag_norder)  
    end if
    if (iflag_ome_ex_text == 'none' ) then
      write(*,*) '   Optical matrix elements (ex) will be read from file'   
    end if

  end subroutine get_ome
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
end module ome