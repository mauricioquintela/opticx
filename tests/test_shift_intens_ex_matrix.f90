program test_shift_intens_ex_matrix
  use constants_math,         only: pi
  use parser_input_file,      only: nw, broadening_type_text
  use parser_optics_xatu_dim, only: norb_ex_cut, npointstotal, vcell, e_ex
  use ome_ex,                  only: xme_ex, vme_ex, xme_ex_inter, vme_ex_inter
  use sigma_second_ex
  implicit none

  integer :: i, j, iw, mode_idx
  real(8), allocatable :: wp(:)
  real(8) :: eta2
  complex(8), allocatable :: sigma_old(:,:,:,:), sigma_new(:,:,:,:)
  real(8) :: max_abs_diff, max_rel_diff, denom
  character(len=16) :: mode_name

  norb_ex_cut  = 9
  nw           = 14
  npointstotal = 100
  vcell        = 1.0d0
  eta2         = 0.05d0

  allocate(e_ex(norb_ex_cut))
  allocate(xme_ex(3,norb_ex_cut), vme_ex(3,norb_ex_cut))
  allocate(xme_ex_inter(3,norb_ex_cut,norb_ex_cut), vme_ex_inter(3,norb_ex_cut,norb_ex_cut))
  allocate(wp(nw))
  allocate(sigma_old(3,3,3,nw), sigma_new(3,3,3,nw))

  do i = 1, norb_ex_cut
    e_ex(i) = 1.5d0 + 0.31d0*dble(i)
    do j = 1, 3
      xme_ex(j,i) = cmplx(0.10d0*i + 0.01d0*j, 0.05d0*i - 0.02d0*j, 8)
      vme_ex(j,i) = cmplx(0.20d0*i - 0.03d0*j, 0.07d0*i + 0.04d0*j, 8)
    end do
  end do
  do i = 1, norb_ex_cut
    do j = 1, norb_ex_cut
      xme_ex_inter(:,i,j) = cmplx(0.05d0*(i+j)+0.011d0, 0.02d0*(i-j)-0.007d0, 8)
      vme_ex_inter(:,i,j) = cmplx(0.08d0*(i-j)+0.013d0, 0.03d0*(i+j)-0.009d0, 8)
    end do
  end do
  do iw = 1, nw
    wp(iw) = 0.4d0 + 0.23d0*dble(iw)
  end do

  do mode_idx = 1, 2
    if (mode_idx == 1) then
      broadening_type_text = 'gaussian';  mode_name = 'gaussian'
    else
      broadening_type_text = 'lorentzian'; mode_name = 'lorentzian'
    end if

    call get_shift_intens_ex(wp, eta2, sigma_old)
    call get_shift_intens_ex_matrix(wp, eta2, sigma_new)

    max_abs_diff = maxval(abs(sigma_old - sigma_new))
    denom        = maxval(abs(sigma_old)) + 1.0d-30
    max_rel_diff = max_abs_diff / denom

    print *, trim(mode_name), ': max abs diff =', max_abs_diff, &
             ', max rel diff =', max_rel_diff

    if (max_rel_diff > 1.0d-9) then
      print *, 'FAIL: ', trim(mode_name), ' mode disagrees beyond tolerance'
      stop 1
    end if
  end do

  print *, 'PASS: matrix-based get_shift_intens_ex_matrix agrees with the reference &
           &implementation in both broadening modes.'

end program test_shift_intens_ex_matrix