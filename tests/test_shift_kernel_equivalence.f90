program test_shift_kernel_equivalence
  use constants_math,        only: pi
  use parser_input_file,     only: nw, broadening_type_text
  use parser_optics_xatu_dim, only: norb_ex_cut, npointstotal, vcell, e_ex
  use ome_ex,                 only: xme_ex, vme_ex, xme_ex_inter, vme_ex_inter
  use sigma_second_ex
  implicit none

  integer, parameter :: nn_test(3)  = [1, 3, 5]
  integer, parameter :: nnp_test(3) = [1, 4, 5]
  integer, parameter :: iw_test(3)  = [1, 3, 6]

  integer     :: nn, nnp, nj, njp, njpp, i, j, iw, mode
  real(8)     :: omegap, omegaq, omega2, eta2
  real(8), allocatable :: wp(:)
  complex(8)  :: s1, s2, s3, ref_val, new_val
  complex(8), allocatable :: d1_arr(:), d2_arr(:), d3_arr(:), d4_arr(:)
  integer     :: nfail
  character(len=16) :: mode_name

  ! --- minimal dummy problem setup: 6 "excitons", 11 frequency points ---
  norb_ex_cut = 6
  nw          = 11
  npointstotal = 100      ! unused by the kernel routines, set for completeness
  vcell        = 1.0d0    ! unused by the kernel routines, set for completeness
  eta2         = 0.05d0

  allocate(e_ex(norb_ex_cut))
  allocate(xme_ex(3,norb_ex_cut))
  allocate(vme_ex(3,norb_ex_cut))
  allocate(xme_ex_inter(3,norb_ex_cut,norb_ex_cut))
  allocate(vme_ex_inter(3,norb_ex_cut,norb_ex_cut))
  allocate(wp(nw))

  ! Deterministic "random-looking" dummy data (no physical meaning) —
  ! just needs to be generic enough to exercise every term non-trivially.
  do i = 1, norb_ex_cut
    e_ex(i) = 1.5d0 + 0.3d0*dble(i)
    do nj = 1, 3
      xme_ex(nj,i) = cmplx(0.10d0*i + 0.01d0*nj, 0.05d0*i - 0.02d0*nj, 8)
      vme_ex(nj,i) = cmplx(0.20d0*i - 0.03d0*nj, 0.07d0*i + 0.04d0*nj, 8)
      do j = 1, norb_ex_cut
        xme_ex_inter(nj,i,j) = cmplx(0.05d0*(i+j) + 0.01d0*nj, 0.02d0*(i-j), 8)
        vme_ex_inter(nj,i,j) = cmplx(0.08d0*(i-j) - 0.02d0*nj, 0.03d0*(i+j), 8)
      end do
    end do
  end do
  do iw = 1, nw
    wp(iw) = 0.5d0 + 0.2d0*dble(iw)
  end do

  ! --- run the comparison for both broadening modes ---
  nfail = 0
  do mode = 1, 2
    if (mode == 1) then
      broadening_type_text = 'gaussian'
      mode_name = 'gaussian'
    else
      broadening_type_text = 'lorentzian'
      mode_name = 'lorentzian'
    end if

    allocate(d1_arr(nw), d2_arr(nw), d3_arr(nw), d4_arr(nw))

    do i = 1, 3
      nn  = nn_test(i)
      nnp = nnp_test(i)

      call get_shift_kernel_ex_dfactors(mode, eta2, wp, nn, nnp, &
                                         d1_arr, d2_arr, d3_arr, d4_arr)

      do j = 1, 3
        iw = iw_test(j)
        omegap = wp(iw); omegaq = -wp(iw); omega2 = 0.0d0

        do nj = 1, 3
          do njp = 1, 3
            do njpp = 1, 3

              call get_shift_kernel_ex_static(mode, nj, njp, njpp, nn, nnp, s1, s2, s3)

              ! reference: old scalar per-iw path
              call get_shift_kernel_ex_freq(mode, eta2, omegap, omegaq, omega2, &
                                             nn, nnp, s1, s2, s3, ref_val)

              ! new: array-based path, evaluated at the same iw
              if (mode == 1) then
                new_val = -( s1*(-cmplx(0.0d0,1.0d0,8)*pi*d1_arr(iw)) &
                           + s2*(-cmplx(0.0d0,1.0d0,8)*pi*d2_arr(iw)) &
                           + s3*(-pi**2*d3_arr(iw)*d4_arr(iw)) )
              else
                new_val = -( s1*d1_arr(iw) + s2*d2_arr(iw) + s3*d3_arr(iw) )
              end if

              if (abs(new_val - ref_val) > 1.0d-10 * (abs(ref_val)+1.0d-30)) then
                nfail = nfail + 1
                print '(A,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,A,2ES14.6,A,2ES14.6)', &
                  'MISMATCH ['//trim(mode_name)//'] nn,nnp,nj,njp,njpp,iw=', &
                  nn, nnp, nj, njp, njpp, iw, '  ref=', ref_val, '  new=', new_val
              end if

            end do
          end do
        end do
      end do
    end do

    deallocate(d1_arr, d2_arr, d3_arr, d4_arr)
  end do

  if (nfail == 0) then
    print *, 'PASS: all', 2*3*3*27*3, 'comparisons agree to 1e-10 relative tolerance.'
  else
    print *, 'FAIL:', nfail, 'mismatches found — see above.'
    stop 1
  end if

  deallocate(e_ex, xme_ex, vme_ex, xme_ex_inter, vme_ex_inter, wp)

end program test_shift_kernel_equivalence