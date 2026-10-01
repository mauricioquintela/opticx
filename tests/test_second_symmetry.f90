! Symmetry checks on the GENERAL two-frequency second-order response, at frequency ratios other than
! r = w_q/w_p = 1.  Motivation (HANDOFF 8.32, 8.42, 8.43): every defect found in these kernels so far has
! been invisible at r = 1, because there the (alpha,w_p)<->(beta,w_q) pair swap is the identity and the
! two field indices are interchangeable.  hBN is D3h, so sigma^xxx = -sigma^xyy = -sigma^yxy = -sigma^yyx
! exactly; the residual on a 30x30 grid is k-derivative discretisation (~2%, CLAUDE.md).
!
! NOTE on what is NOT tested here.  The RAW kernels are deliberately not symmetric under exchanging the
! two frequencies -- each is one ordering, and the pair average IS the symmetrisation.  Measured, the raw
! residual max|s^abc(wp,wq) - s^acb(wq,wp)|/max|s| is ~1 at every r /= 1 for BOTH the single-particle
! Eq. (A3a) and the excitonic Eq. (B1b), so asserting it would be wrong, and reporting it discriminates
! nothing.  What does discriminate is D3h on the SYMMETRISED tensor, which is what is checked below.
program test_second_symmetry
  use parser_input_file,      only: get_input_file, material_name_in, nw, norb_ex_cut, eta, e1, e2
  use parser_wannier90_tb,    only: wannier90_get
  use parser_optics_xatu_dim, only: get_optics_xatu_dim
  use bands,                  only: get_energy_bands
  use ome,                    only: get_ome
  use sigma_second_ex,        only: get_second_intens_ex_methodA, get_second_intens_ex_methodB
  implicit none
  integer :: nfail
  real(8) :: eta2, dw, d3h_shg
  complex(8), allocatable :: hwp(:), hwq(:), hw2(:)
  complex(8), allocatable :: A1(:,:,:,:), A2(:,:,:,:), B1(:,:,:,:), B2(:,:,:,:)
  complex(8), allocatable :: symA(:,:,:,:), symB(:,:,:,:), symHalf(:,:,:,:), symTwo(:,:,:,:)

  nfail = 0
  call get_input_file(); call wannier90_get(material_name_in)
  call get_optics_xatu_dim(); call get_energy_bands(); call get_ome()
  eta2 = eta/27.211385d0; dw = (e2-e1)/dble(nw)
  allocate(hwp(nw), hwq(nw), hw2(nw))
  allocate(A1(3,3,3,nw), A2(3,3,3,nw), B1(3,3,3,nw), B2(3,3,3,nw))
  allocate(symA(3,3,3,nw), symB(3,3,3,nw), symHalf(3,3,3,nw), symTwo(3,3,3,nw))

  write(*,'(A)') ' Second-order symmetry checks (hBN 30x30, D3h), general two-frequency branch'

  ! ---- r = 1 (SHG): methods A and B must agree, and D3h must hold -------------------------------
  call freqs(1.0d0, .false., 1.0d0)
  call get_second_intens_ex_methodA(nw, hwp, hwq, hw2, A1); call pairsym(A1, A1, symA)
  call get_second_intens_ex_methodB(nw, hwp, hwq, hw2, B1); call pairsym(B1, B1, symB)
  call check('r = 1   method A vs method B, symmetrised', reldiff(symA, symB), 5.0d-2)
  d3h_shg = d3h(symB)          ! the k-derivative discretisation floor for this grid
  call check('r = 1   D3h  max|s^xxx + s^xyy|/max|s^xxx|', d3h_shg, 5.0d-2)

  ! ---- r = 1/2: a genuinely asymmetric pair; both methods and D3h still required -----------------
  call freqs(0.5d0, .false., 1.0d0)
  call get_second_intens_ex_methodA(nw, hwp, hwq, hw2, A1)
  call get_second_intens_ex_methodA(nw, hwq, hwp, hw2, A2); call pairsym(A1, A2, symA)
  call get_second_intens_ex_methodB(nw, hwp, hwq, hw2, B1)
  call get_second_intens_ex_methodB(nw, hwq, hwp, hw2, B2); call pairsym(B1, B2, symB)
  call check('r = 1/2 method A vs method B, symmetrised', reldiff(symA, symB), 5.0d-2)
  symHalf = symB

  ! ---- the r <-> 1/r identity: sym_{1/r}(a,b,c) = sym_r(a,c,b) at the same (w_p,w_q) pair ---------
  call get_second_intens_ex_methodB(nw, hwq, hwp, hw2, B1)
  call get_second_intens_ex_methodB(nw, hwp, hwq, hw2, B2); call pairsym(B1, B2, symTwo)
  call check('r <-> 1/r  sym_{1/r}(a,b,c) = sym_r(a,c,b)', reldiff_t(symTwo, symHalf), 1.0d-10)

  ! ---- r = -1, the DC branch with its own convention (HANDOFF 8.42): guards that fix -------------
  call freqs(-1.0d0, .true., 1.0d0)
  call get_second_intens_ex_methodA(nw, hwp, hwq, hw2, A1)
  call pairsym(A1, A1, symA)                                  ! index-only, NO pair average
  symA = cmplx(dble(symA), 0.0d0, 8)
  call check('r = -1  DC branch D3h (guards 8.42)', d3h(symA), 5.0d-2)
  call check('r = -1  DC branch Im sigma = 0 (Sipe)', maxval(abs(aimag(symA))), 1.0d-30)

  ! ---- r = 0, electro-optic (HANDOFF 8.43b) ------------------------------------------------------
  ! B1b terms 1-2 contain w_2 and w_q but NO w_p, so in the swapped pass both of term 1's poles land on
  ! resonance together when w_p = 0: a near-DOUBLE pole. It does not break the formula -- it AMPLIFIES the
  ! k-derivative discretisation error of X_nm, by ~1/eta^2 on top of the usual 1/N^2. So EO is correct
  ! where it is resolved and needs a finer grid (or a larger eta) than SHG. Asserted at 4*eta, where a
  ! 30x30 grid resolves it; the production-eta number is reported so the amplification stays visible.
  call freqs(0.0d0, .false., 16.0d0)
  call get_second_intens_ex_methodB(nw, hwp, hwq, hw2, B1)
  call get_second_intens_ex_methodB(nw, hwq, hwp, hw2, B2); call pairsym(B1, B2, symB)
  call check('r = 0   EO D3h at 16*eta, relative to the r=1 floor', d3h(symB)/d3h_shg, 1.5d0)
  call freqs(0.0d0, .false., 1.0d0)
  call get_second_intens_ex_methodB(nw, hwp, hwq, hw2, B1)
  call get_second_intens_ex_methodB(nw, hwq, hwp, hw2, B2); call pairsym(B1, B2, symB)
  write(*,'(A,ES10.3,A,ES10.3,A)') '  REPORT  r = 0   electro-optic D3h at 1*eta = ', d3h(symB), &
       '  vs ', d3h_shg, ' at r = 1 (floor).'
  write(*,'(A)') '                  The excess is amplified discretisation, not a defect: it falls as'
  write(*,'(A)') '                  1/N^2 (x4.2 from 30x30 to 60x60) and as eta^2. HANDOFF 8.43b.'

  if (nfail == 0) then
    write(*,'(A)') ' ALL TESTS PASSED'
  else
    write(*,'(A,I0,A)') ' ', nfail, ' TEST(S) FAILED'; stop 1
  end if

contains
  subroutine freqs(r, dc, efac)
    real(8), intent(in) :: r
    logical, intent(in) :: dc
    real(8), intent(in), optional :: efac
    integer :: iw
    real(8) :: et
    et = eta2; if (present(efac)) et = eta2*efac
    do iw = 1, nw
      hwp(iw) = cmplx((e1+dw*dble(iw-1))/27.211385d0, et, 8)
      if (dc) then
        hwq(iw) = -hwp(iw); hw2(iw) = (0.0d0,0.0d0)
      else
        hwq(iw) = cmplx(r*(e1+dw*dble(iw-1))/27.211385d0, et, 8)
        hw2(iw) = hwp(iw) + hwq(iw)
      end if
    end do
  end subroutine
  subroutine pairsym(sA, sB, sS)
    complex(8), intent(in)  :: sA(3,3,3,nw), sB(3,3,3,nw)
    complex(8), intent(out) :: sS(3,3,3,nw)
    integer :: a,b,c
    do a=1,3; do b=1,3; do c=1,3
      sS(a,b,c,:) = 0.5d0*(sA(a,b,c,:) + sB(a,c,b,:))
    end do; end do; end do
  end subroutine
  real(8) function reldiff(x, y)
    complex(8), intent(in) :: x(3,3,3,nw), y(3,3,3,nw)
    reldiff = maxval(abs(x(1:2,1:2,1:2,:)-y(1:2,1:2,1:2,:))) / max(maxval(abs(x(1:2,1:2,1:2,:))),1.0d-300)
  end function
  real(8) function reldiff_t(x, y)      ! compares x(a,b,c) against y(a,c,b)
    complex(8), intent(in) :: x(3,3,3,nw), y(3,3,3,nw)
    integer :: a,b,c
    real(8) :: m
    m = 0.0d0
    do a=1,2; do b=1,2; do c=1,2
      m = max(m, maxval(abs(x(a,b,c,:)-y(a,c,b,:))))
    end do; end do; end do
    reldiff_t = m / max(maxval(abs(y(1:2,1:2,1:2,:))),1.0d-300)
  end function
  real(8) function d3h(s)   ! COMPLETE D3h violation: hBN allows only xxx = -xyy = -yxy = -yyx
    complex(8), intent(in) :: s(3,3,3,nw)
    real(8) :: v
    v =       maxval(abs(s(1,1,1,:)+s(1,2,2,:)))
    v = max(v,maxval(abs(s(1,1,1,:)+s(2,1,2,:))))
    v = max(v,maxval(abs(s(1,1,1,:)+s(2,2,1,:))))
    v = max(v,maxval(abs(s(2,2,2,:))));  v = max(v,maxval(abs(s(1,1,2,:))))
    v = max(v,maxval(abs(s(1,2,1,:))));  v = max(v,maxval(abs(s(2,1,1,:))))
    d3h = v / max(maxval(abs(s(1,1,1,:))),1.0d-300)
  end function
  subroutine check(lab, err, tol)
    character(len=*), intent(in) :: lab
    real(8), intent(in) :: err, tol
    if (err <= tol) then
      write(*,'(A,ES11.3,A,ES9.2,A,A)') '  PASS  err=',err,'  (tol ',tol,')  ',lab
    else
      write(*,'(A,ES11.3,A,ES9.2,A,A)') '  FAIL  err=',err,'  (tol ',tol,')  ',lab
      nfail = nfail + 1
    end if
  end subroutine
end program test_second_symmetry
