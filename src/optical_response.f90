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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) '7. Entering optical_response'
    !optical response selection. The valid keyword list lives HERE, once: an unknown Response
    !falls into `case default` and stops, so adding a branch cannot leave the validator stale
    !(it used to be a separate five-line negated conjunction repeating all nine names).
    select case (trim(response_text))

      case ('none')
        write(*,*) '   No optical response has been evaluated'

      case ('absorbance')
        write(*,*) '   Optical response: absorbance'
        call get_sigma_first_sp()
        if (iflag_xatu) call get_sigma_first_ex()

      case ('shift_sumrule', 'shift_shiftvector', 'shift_gender')
        call second_order(1, -1, 'shift conductivity')

      case ('shg')
        call second_order(1, 1, 'shg susceptibility')

      case ('electrooptic')
        call second_order(1, 0, 'electro-optic (Pockels) susceptibility, sigma(w; w, 0)')

      case ('rectification')
        ! Excitonic branch goes through METHOD A, Eq. (B1a): it has no i*hbar*omega_2 prefactor and
        ! so stays finite at omega_2 = 0, where method B vanishes identically (HANDOFF 8.38).
        call second_order(1, 1, 'optical rectification, sigma(0; w, -w) via Eq. (A3a)')

      case ('general')
        call second_order(1, 1, 'general second order, sigma(w_p+w_q; w_p, w_q)')

      case default
        write(*,*) 'ERROR (optical_response): unknown Response = "'//trim(response_text)//'".'
        write(*,*) '       Valid: none, absorbance, shift_sumrule, shift_shiftvector, shift_gender,'
        write(*,*) '              shg, electrooptic, rectification, general.'
        write(*,*) '       Note it is case sensitive. Nothing would have been computed; stopping.'
        stop 1

    end select
    write(*,*) 'The optical response has been evaluated'

  contains
    ! Every second-order branch does the same three things and differs only in the frequency pair
    ! and the label, so they are one routine rather than five copies.
    subroutine second_order(nwp, nwq, label)
      integer,          intent(in) :: nwp, nwq
      character(len=*), intent(in) :: label
      write(*,*) '   Optical response: '//label
      call get_sigma_second_sp(nwp, nwq)
      if (iflag_xatu) call get_sigma_second_ex(nwp, nwq)
    end subroutine second_order

  end subroutine get_optical_response
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module optical_response
