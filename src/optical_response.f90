module optical_response
  use parser_input_file, &
  only:response_text
  use parser_input_file, &
  only:iflag_xatu
  use sigma_first_sp, &
  only:get_sigma_first_sp
  use sigma_first_ex, &
  only:get_sigma_first_ex
  use sigma_second_sp, &
  only:get_sigma_second_sp
  use sigma_second_ex, &
  only:get_sigma_second_ex
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_optical_response()
    implicit none
    integer :: nwp,nwq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    !write(*,*) iflag_ome_sp
    !write(*,*) iflag_ome_ex
    write(*,*) '7. Entering optical_response'
    !optical response selection
    if (response_text  == 'none') then
      write(*,*) '   No optical response has been evaluated'  
    end if
    if (response_text  == 'absorbance') then
      write(*,*) '   Optical response: absorbance'  
      !Evaluate sigma_first single-particle
      call get_sigma_first_sp()
      !Evaluate sigma_first exciton
      if (iflag_xatu .eqv. .true.) then
        call get_sigma_first_ex()
      end if
      
    end if
    if (response_text  == 'shift_sumrule' .or. response_text  == 'shift_shiftvector' &
        .or. response_text  == 'shift_gender') then
      write(*,*) '   Optical response: shift conductivity' 
      nwp=1
      nwq=-1
      call get_sigma_second_sp(nwp,nwq) 
      if (iflag_xatu .eqv. .true.) then
        call get_sigma_second_ex(nwp,nwq) 
      end if
    end if
    if (response_text /= 'none' .and. response_text /= 'absorbance' .and. &
        response_text /= 'shift_sumrule' .and. response_text /= 'shift_shiftvector' .and. &
        response_text /= 'shift_gender') then
      write(*,*) 'ERROR (optical_response): unknown Response = "'//trim(response_text)//'".'
      write(*,*) '       Valid: none, absorbance, shift_sumrule, shift_shiftvector, shift_gender.'
      write(*,*) '       Note it is case sensitive. Nothing would have been computed; stopping.'
      stop 1
    end if
    write(*,*) 'The optical response has been evaluated'
  end subroutine get_optical_response
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
end module optical_response