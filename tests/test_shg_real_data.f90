! Real-data SHG checks on a Xatu exciton set (hBN). Runs the same pipeline as the main program
! (input file given as command-line argument 1), then checks, with tolerances:
!  1. X_nm is Hermitian: X_nm = X_mn^*  (a periodic k-derivative in get_fk_ex_der_k guarantees it;
!     zeroing the zone-edge points does not).
!  2. The production method-B zgemm routine reproduces the scalar method-B kernels.
!  3. Method A with the Heisenberg momentum built from position elements,
!        Pi_n = -i E_n X_n,   Pi_nm = i (E_n - E_m) X_nm,
!     agrees with method B (Taghizadeh & Pedersen, PRB 97, 205432, Eq. 8):
!       (a) diagonal components xxx, yyy, zzz: to round-off (needs X_nm Hermitian);
!       (b) the full (b,c)-symmetrised tensor: only to the discretisation error of the k-derivative.
!           The proof needs commuting position components, [X^a,X^b] = 0, which the finite-difference
!           Q_nm satisfies only approximately; measured for hBN, 100-150 excitons: 1.8e-3 (30x30 grid),
!           8.6e-4 (45x45), 5.7e-4 (60x60), i.e. ~1/N^2. The tolerance is set with margin above the
!           30x30 value and is still ~50x below the 53% error of the bare-P method A.
! Also prints (INFO, not a pass/fail) how far method A is from B when it is fed the momentum
! elements vme_ex/vme_ex_inter that the code actually computes: those are the BARE momentum
! P_n, P_nm, not Pi, so that combination is not the paper's method A.
! Run through `make run_test_shg_real` (creates a scratch directory; the pipeline writes files).
program test_shg_real_data
  use parser_input_file,      only: get_input_file, material_name_in, nw, norb_ex_cut
  use parser_wannier90_tb,    only: wannier90_get
  use parser_optics_xatu_dim, only: get_optics_xatu_dim, e_ex, npointstotal, vcell
  use bands,                  only: get_energy_bands
  use ome,                    only: get_ome
  use ome_ex,                 only: xme_ex, vme_ex, xme_ex_inter, vme_ex_inter, inter_terms_ready
  use sigma_second_sp,        only: initialize_sigma_second_arrays
  use sigma_second_ex,        only: get_shg_intens_ex_matrix, get_shg_intens_ex_matrix_methodB, &
                                    get_shg_kernel_ex_static_methodB, get_shg_kernel_ex_freq_methodB
  implicit none
  real(8), parameter :: tol_herm = 1.0d-12
  real(8), parameter :: tol_matrix_vs_scalar = 1.0d-12
  real(8), parameter :: tol_A_vs_B_diag = 1.0d-10
  real(8), parameter :: tol_A_vs_B_full = 1.0d-2
  integer :: nfail, n, a, b, c, i, j, nn, nnp, iw
  real(8) :: eta2, err, sc, feps
  real(8), allocatable :: wp(:)
  complex(8), parameter :: iu = (0.0d0,1.0d0)
  complex(8), allocatable :: sA(:,:,:,:), sB(:,:,:,:), sBs(:,:,:,:), sBm(:,:,:,:), sAs(:,:,:,:), sAP(:,:,:,:), &
                             dummy(:,:,:,:), vsave(:,:), vintsave(:,:,:), sAPs(:,:,:,:)
  complex(8) :: s1, s2, s3, kk

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
  allocate(wp(nw), dummy(3,3,3,nw), sA(3,3,3,nw), sB(3,3,3,nw), sBs(3,3,3,nw), sBm(3,3,3,nw), sAs(3,3,3,nw), &
           sAP(3,3,3,nw), vsave(3,n), vintsave(3,n,n), sAPs(3,3,3,nw))
  call initialize_sigma_second_arrays(nw, wp, eta2, dummy)
  print '(A,I0,A,I0,A,I0,A,F6.3,A)', 'real-data SHG test: ', n, ' excitons, ', npointstotal, ' k-points, ', nw, &
        ' frequencies, eta = ', eta2*27.211385d0, ' eV'

  ! ---- 1. Hermiticity of X_nm -------------------------------------------------------------
  err = 0.0d0; sc = 0.0d0
  do a = 1, 3
    do i = 1, n
      do j = 1, n
        err = max(err, abs(xme_ex_inter(a,i,j) - conjg(xme_ex_inter(a,j,i))))
        sc  = max(sc, abs(xme_ex_inter(a,i,j)))
      end do
    end do
  end do
  call report('X_nm Hermitian: max|X_nm - X_mn^*| / max|X_nm|', err/sc, tol_herm)

  ! ---- method B, production (matrix) routine ---------------------------------------------------
  call get_shg_intens_ex_matrix_methodB(wp, eta2, sBm)

  ! ---- 2. matrix B vs scalar B kernels -----------------------------------------------------------
  sB = (0.0d0, 0.0d0)
  do nn = 1, n
    do nnp = 1, n
      do a = 1, 3
        do b = 1, 3
          do c = 1, 3
            call get_shg_kernel_ex_static_methodB(a, b, c, nn, nnp, s1, s2, s3)
            do iw = 1, nw
              call get_shg_kernel_ex_freq_methodB(eta2, wp(iw), nn, nnp, s1, s2, s3, kk)
              sB(a,b,c,iw) = sB(a,b,c,iw) + kk/(dble(npointstotal)*vcell)
            end do
          end do
        end do
      end do
    end do
  end do
  sc = maxval(abs(sBm))
  call report('scalar B vs matrix B: max|diff| / max|sigma_B|', maxval(abs(sB-sBm))/sc, tol_matrix_vs_scalar)

  ! ---- 3. method A with Pi built from X vs method B --------------------------------------------
  vsave = vme_ex(:,1:n); vintsave = vme_ex_inter(:,1:n,1:n)              ! keep the code's own arrays
  ! INFO first: method A on the momentum elements the code computes (bare P)
  call get_shg_intens_ex_matrix(wp, eta2, sAP)
  call sym_bc(sAP, sAs); call sym_bc(sBm, sBs)
  sAPs = sAs                                                              ! bare-P method A, kept for the plot
  print '(A,ES10.3,A)', '  INFO   method A fed the code''s vme_ex (bare momentum P) vs B, symmetrised: ', &
        maxval(abs(sAs-sBs))/maxval(abs(sBs)), '  (not a test: P /= Pi)'
  do a = 1, 3
    do i = 1, n
      vme_ex(a,i) = -iu*e_ex(i)*xme_ex(a,i)
      do j = 1, n
        vme_ex_inter(a,i,j) = iu*(e_ex(i)-e_ex(j))*xme_ex_inter(a,i,j)
      end do
    end do
  end do
  call get_shg_intens_ex_matrix(wp, eta2, sA)
  vme_ex(:,1:n) = vsave; vme_ex_inter(:,1:n,1:n) = vintsave
  call sym_bc(sA, sAs)
  err = 0.0d0
  do a = 1, 3
    err = max(err, maxval(abs(sA(a,a,a,:) - sBm(a,a,a,:))))
  end do
  call report('method A (Pi from X) vs B, diagonal xxx/yyy/zzz, raw', err/maxval(abs(sBm)), tol_A_vs_B_diag)
  call report('method A (Pi from X) vs B, (b,c)-symmetrised, full tensor', &
              maxval(abs(sAs-sBs))/maxval(abs(sBs)), tol_A_vs_B_full)

  ! spectra for tools/plot_test_outputs.py (uA nm/V^2, hbar*omega axis -- the FUNDAMENTAL, matching
  ! opticx's own shg_*_lengthgauge files since the 2026-09-24 axis change, HANDOFF 8.33)
  feps = 6.623618d-03*1.0d+06*(27.211386d0**(-2))*5.291772d-11*1.0d+09
  open(77, file='shg_real_data_spectra.dat')
  write(77,'(A)') '# kind: shg_real_data'
  write(77,'(A)') '# columns: Ew_eV xxx_B_re xxx_B_im xxx_A_Pi_re xxx_A_Pi_im xxx_A_P_re xxx_A_P_im xyy_B_re xyy_B_im xyy_A_Pi_re xyy_A_Pi_im'
  do iw = 1, nw
    write(77,'(11ES16.8)') wp(iw)*27.211385d0, &
      feps*real(sBs(1,1,1,iw)),  feps*aimag(sBs(1,1,1,iw)),  feps*real(sAs(1,1,1,iw)),  feps*aimag(sAs(1,1,1,iw)), &
      feps*real(sAPs(1,1,1,iw)), feps*aimag(sAPs(1,1,1,iw)), &
      feps*real(sBs(1,2,2,iw)),  feps*aimag(sBs(1,2,2,iw)),  feps*real(sAs(1,2,2,iw)),  feps*aimag(sAs(1,2,2,iw))
  end do
  close(77)

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

  subroutine sym_bc(s, ss)
    complex(8), intent(in)  :: s(3,3,3,nw)
    complex(8), intent(out) :: ss(3,3,3,nw)
    integer :: ii, jj, kk2
    do ii = 1, 3
      do jj = 1, 3
        do kk2 = 1, 3
          ss(ii,jj,kk2,:) = 0.5d0*(s(ii,jj,kk2,:) + s(ii,kk2,jj,:))
        end do
      end do
    end do
  end subroutine sym_bc
end program test_shg_real_data
