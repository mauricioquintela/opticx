! Real-data check of the excitonic shift conductivity (Lorentzian branch) against the paper's
! position-only formula, evaluated independently here.
! Esteve-Paredes et al., npj Comput. Mater. 11, 13, Eq. (10)/(11) and SI Note 5:
!   sigma^{abc}(0;w,-w) = pi/(Nk V) * sum_{N,N'} Re S^{abc}_{NN'} delta_eta(w - E_N),
!   S^{abc}_{NN'} = X^a_N X^b_{NN'} X^c*_{N'},
! symmetrised over the field indices (b,c) (IPA Eq. 9 is written with I^{abc} + I^{acb}; SI: the real part is
! symmetric under b<->c). Convention: omega_p = omega + i*eta, omega_q = -omega_p, omega_2 = 0.
! What this catches (none of the scalar-vs-matrix tests can): a dispersive instead of absorptive line
! shape (zero crossing at E_N), a wrong overall sign, and a non-vanishing term 3 (with +i*eta on every pole it
! cancelled ~80% of the first peak). The reference contains only the resonant delta term; the code also has
! the (non-resonant) term 1 and finite-eta pole structure, so the agreement is not round-off: measured
! max|diff|/max|ref| = 1.7e-2 for hBN, 30x30 grid, 100 excitons. The tolerance keeps ~3x margin and is
! well below the error of the previous convention (order 1, see the control run in HANDOFF.md section 8).
! Run through `make run_test_shift_real`.
program test_shift_real_data
  use constants_math,         only: pi
  use parser_input_file,      only: get_input_file, material_name_in, nw, norb_ex_cut
  use parser_wannier90_tb,    only: wannier90_get
  use parser_optics_xatu_dim, only: get_optics_xatu_dim, e_ex, npointstotal, vcell
  use bands,                  only: get_energy_bands
  use ome,                    only: get_ome
  use ome_ex,                 only: xme_ex, xme_ex_inter, inter_terms_ready
  use sigma_second_sp,        only: initialize_sigma_second_arrays
  use sigma_second_ex,        only: get_shift_intens_ex_matrix
  implicit none
  real(8), parameter :: tol_ref = 5.0d-2
  integer :: nfail, n, a, b, c, i, j, iw, ipeak
  real(8) :: eta2, err, sc, l, sN, feps
  real(8), allocatable :: wp(:), ref(:,:,:,:), sre(:,:,:,:)
  complex(8), allocatable :: sig(:,:,:,:), sigs(:,:,:,:), dummy(:,:,:,:)

  nfail = 0
  call get_input_file()
  call wannier90_get(material_name_in)
  call get_optics_xatu_dim()
  call get_energy_bands()
  call get_ome()
  if (.not. inter_terms_ready) then
    print *, 'FAIL: excitonic inter-exciton matrix elements not available (need OME_ex = nonlinear)'
    stop 1
  end if
  n = norb_ex_cut
  allocate(wp(nw), dummy(3,3,3,nw), sig(3,3,3,nw), sigs(3,3,3,nw), ref(3,3,3,nw), sre(3,3,3,n))
  call initialize_sigma_second_arrays(nw, wp, eta2, dummy)

  call get_shift_intens_ex_matrix(wp, eta2, sig)
  do c = 1, 3
    do b = 1, 3
      do a = 1, 3
        sigs(a,b,c,:) = 0.5d0*(sig(a,b,c,:) + sig(a,c,b,:))
      end do
    end do
  end do

  ! sre(a,b,c,N) = sum_N' (1/2) Re[S^{abc}_{NN'} + S^{acb}_{NN'}]
  sre = 0.0d0
  do i = 1, n
    do j = 1, n
      do c = 1, 3
        do b = 1, 3
          do a = 1, 3
            sre(a,b,c,i) = sre(a,b,c,i) + 0.5d0*real( xme_ex(a,i)*xme_ex_inter(b,i,j)*conjg(xme_ex(c,j)) &
                                                    + xme_ex(a,i)*xme_ex_inter(c,i,j)*conjg(xme_ex(b,j)) )
          end do
        end do
      end do
    end do
  end do
  ref = 0.0d0
  do iw = 1, nw
    do i = 1, n
      l = (eta2/pi)/((wp(iw)-e_ex(i))**2 + eta2**2)                 ! Lorentzian delta, 1/Hartree
      ref(:,:,:,iw) = ref(:,:,:,iw) + pi/(dble(npointstotal)*vcell) * l * sre(:,:,:,i)
    end do
  end do
  print '(A,I0,A,I0,A,I0,A,F6.3,A)', 'real-data shift test: ', n, ' excitons, ', npointstotal, ' k-points, ', nw, &
        ' frequencies, eta = ', eta2*27.211385d0, ' eV'

  sc = maxval(abs(ref))
  ! spectra for tools/plot_test_outputs.py (uA nm/V^2)
  feps = 6.623618d-03*1.0d+06*(27.211386d0**(-2))*5.291772d-11*1.0d+09
  open(77, file='shift_real_data_spectra.dat')
  write(77,'(A)') '# kind: shift_real_data'
  write(77,'(A)') '# columns: E_eV xxx_code xxx_ref xyy_code xyy_ref yxy_code yxy_ref'
  do iw = 1, nw
    write(77,'(7ES16.8)') wp(iw)*27.211385d0, feps*real(sigs(1,1,1,iw)), feps*ref(1,1,1,iw), &
      feps*real(sigs(1,2,2,iw)), feps*ref(1,2,2,iw), feps*real(sigs(2,1,2,iw)), feps*ref(2,1,2,iw)
  end do
  close(77)

  err = maxval(abs(real(sigs) - ref))/sc
  call report('shift (Lorentzian, b<->c-symmetrised) vs paper Eq. 11: max|diff| / max|ref|', err, tol_ref)

  ! the maximum of the dominant component must be an absorptive peak at an exciton energy, with the
  ! sign of Eq. 11 (a dispersive line shape puts the extrema at E_N +- eta and a zero crossing at E_N)
  ipeak = maxloc(abs(ref(1,1,1,:)), dim=1)
  sN = minval(abs(wp(ipeak) - e_ex(1:n)))
  call report('xxx peak at an exciton energy: |w_peak - E_N| / eta', sN/eta2, 0.6d0)
  if (real(sigs(1,1,1,ipeak))*ref(1,1,1,ipeak) <= 0.0d0) then
    print '(A)', '  FAIL   xxx at its peak has the wrong sign relative to Eq. 11'
    nfail = nfail + 1
  else
    print '(A)', '  PASS   xxx at its peak has the sign of Eq. 11'
  end if

  if (nfail == 0) then
    print '(/,A)', 'ALL TESTS PASSED'
  else
    print '(/,A,I0,A)', 'FAILED: ', nfail, ' check(s)'
    stop 1
  end if

contains
  subroutine report(label, val, tol)
    character(*), intent(in) :: label
    real(8), intent(in) :: val, tol
    if (val < tol) then
      print '(A,ES10.3,A,ES8.1,A,A)', '  PASS   err=', val, '  (tol ', tol, ')  ', label
    else
      print '(A,ES10.3,A,ES8.1,A,A)', '  FAIL   err=', val, '  (tol ', tol, ')  ', label
      nfail = nfail + 1
    end if
  end subroutine report
end program test_shift_real_data
