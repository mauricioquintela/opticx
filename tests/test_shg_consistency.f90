! Consolidated SHG consistency test (replaces test_shg_equivalence / _identity_check /
! _pointwise_kernel / _pointwise_kernel_eta).
!
! What is (and is not) expected to hold -- all derived from Taghizadeh & Pedersen PRB 97, 205432:
!  * sigma^A (B1a) == sigma^B (B1b) requires  omega_2 = omega_p + omega_q  as COMPLEX numbers
!    (paper: omega -> omega + i*eta for every frequency  =>  omega_2 -> 2*omega + 2i*eta).
!  * It requires physical matrix elements: commuting position components [R^a,R^b]=0, X_nm defined
!    relative to the ground state (paper Eq. B2b), and Pi_n=-iE_nX_n, Pi_nm=i(E_n-E_m)X_nm.
!    Time-reversal symmetry is NOT needed (complex Hermitian data works).
!  * Off-diagonal tensor components only agree after symmetrising the two field indices (b<->c).
!    Diagonal components (xxx,yyy,zzz) agree without it and even for arbitrary Hermitian data.
program test_shg_consistency
  use parser_input_file,      only: nw
  use parser_optics_xatu_dim, only: norb_ex_cut, npointstotal, vcell, e_ex
  use ome_ex,                 only: xme_ex, vme_ex, xme_ex_inter, vme_ex_inter
  use sigma_second_ex
  implicit none
  integer, parameter :: nex = 10, dim = nex+1
  real(8), parameter :: tol = 1.0d-12
  integer :: nfail, icplx, ie, a, b, c, i, j, k, iw
  real(8) :: etas(2), wp(14), eta2, rr(dim,3), r(3,nex,2), r2(3,nex,nex), ri(3,nex,nex)
  real(8) :: hr(dim,dim), hi(dim,dim), err, sc
  complex(8) :: V(dim,dim), Rm(dim,dim,3), iu, dotp
  complex(8), allocatable :: sA(:,:,:,:), sB(:,:,:,:), sAm(:,:,:,:), sAs(:,:,:,:), sBs(:,:,:,:), sBm(:,:,:,:)

  iu = cmplx(0d0,1d0,8); nfail = 0
  nw = 14; norb_ex_cut = nex; npointstotal = 100; vcell = 1d0
  allocate(e_ex(nex), xme_ex(3,nex), vme_ex(3,nex), xme_ex_inter(3,nex,nex), vme_ex_inter(3,nex,nex))
  allocate(sA(3,3,3,nw), sB(3,3,3,nw), sAm(3,3,3,nw), sAs(3,3,3,nw), sBs(3,3,3,nw), sBm(3,3,3,nw))
  call random_seed(put=(/ (31*i+7, i=1,64) /))
  do i=1,nex; e_ex(i) = 1.5d0 + 0.31d0*i + 0.02d0*mod(7*i,5); end do
  do iw=1,nw; wp(iw) = 0.4d0 + 0.23d0*iw; end do
  etas = (/ 0.0d0, 0.05d0 /)

  ! ---------- Part 1: physical data (commuting positions) --------------------------------
  do icplx = 0, 1
    call random_number(hr); call random_number(hi); hr = hr-0.5d0; hi = (hi-0.5d0)*dble(icplx)
    V = cmplx(hr,hi,8)
    do j=1,dim
      do k=1,j-1
        dotp = sum(conjg(V(:,k))*V(:,j)); V(:,j) = V(:,j) - dotp*V(:,k)
      end do
      V(:,j) = V(:,j)/sqrt(sum(abs(V(:,j))**2))
    end do
    call random_number(rr); rr = rr-0.5d0
    do a=1,3
      do i=1,dim; do j=1,dim
        Rm(i,j,a) = sum(V(i,:)*rr(:,a)*conjg(V(j,:)))
      end do; end do
    end do
    do a=1,3
      do i=1,nex
        xme_ex(a,i) = Rm(1,i+1,a)
        do j=1,nex
          xme_ex_inter(a,i,j) = Rm(i+1,j+1,a)
          if (i==j) xme_ex_inter(a,i,j) = xme_ex_inter(a,i,j) - Rm(1,1,a)   ! paper Eq. (B2b)
        end do
      end do
    end do
    call set_velocity_from_position()
    do ie = 1, 2
      eta2 = etas(ie)
      call get_shg_intens_ex(wp, eta2, sA)
      call get_shg_intens_ex_matrix(wp, eta2, sAm)
      call B_direct(wp, eta2, sB)
      call get_shg_intens_ex_matrix_methodB(wp, eta2, sBm)
      call sym_bc(sA, sAs); call sym_bc(sB, sBs)
      sc = maxval(abs(sA))
      call report('scalar A vs matrix A',                    maxval(abs(sA-sAm))/sc, icplx, eta2)
      call report('scalar B vs matrix B (production routine)', maxval(abs(sB-sBm))/sc, icplx, eta2)
      call report('A vs B, (b,c)-symmetrised, full tensor',  maxval(abs(sAs-sBs))/sc, icplx, eta2)
      err = 0d0
      do a=1,3; err = max(err, maxval(abs(sA(a,a,a,:)-sB(a,a,a,:)))); end do
      call report('A vs B, diagonal xxx/yyy/zzz, raw',       err/sc, icplx, eta2)
      if (icplx == 1 .and. ie == 2) then      ! complex-Hermitian data, eta = 0.05: spectrum for tools/plot_test_outputs.py
        open(77, file='shg_consistency_spectra.dat')
        write(77,'(A)') '# kind: shg_consistency'
        write(77,'(A)') '# columns: omega xxx_B_re xxx_B_im xxx_A_scalar_re xxx_A_scalar_im xxx_A_matrix_re xxx_A_matrix_im'
        do iw = 1, nw
          write(77,'(7ES16.8)') wp(iw), real(sB(1,1,1,iw)), aimag(sB(1,1,1,iw)), real(sA(1,1,1,iw)), aimag(sA(1,1,1,iw)), &
                                real(sAm(1,1,1,iw)), aimag(sAm(1,1,1,iw))
        end do
        close(77)
      end if
    end do
  end do

  ! ---------- Part 2: arbitrary Hermitian data: only the diagonal must agree ---------------
  call random_number(r); call random_number(r2); call random_number(ri)
  r = r-0.5d0; r2 = r2-0.5d0; ri = ri-0.5d0
  do i=1,nex; do a=1,3; xme_ex(a,i) = cmplx(r(a,i,1), r(a,i,2), 8); end do; end do
  do i=1,nex; do j=1,nex
    xme_ex_inter(:,i,j) = cmplx(0.5d0*(r2(:,i,j)+r2(:,j,i)), 0.5d0*(ri(:,i,j)-ri(:,j,i)), 8)
  end do; end do
  call set_velocity_from_position()
  eta2 = 0.05d0
  call get_shg_intens_ex(wp, eta2, sA); call B_direct(wp, eta2, sB)
  err = 0d0; sc = 0d0
  do a=1,3
    err = max(err, maxval(abs(sA(a,a,a,:)-sB(a,a,a,:)))); sc = max(sc, maxval(abs(sA(a,a,a,:))))
  end do
  call report('arbitrary Hermitian data: diagonal A vs B', err/sc, 2, eta2)

  ! ---------- Part 3: per-pair closed-form identity, ANY (a,b,c) ---------------------------
  ! kernelA - kernelB = -i*s1B/(z-Em) - i*s2B/(z+Em) + i*s3B*[1/(z-Em)+1/(z+En) + (z2-2z)*d3]
  block
    integer :: nn, nnp
    complex(8) :: s1A,s2A,s3A,s1B,s2B,s3B,kA,kB,z,z2,d3,pred
    real(8) :: om
    om = 1.7d0; eta2 = 0.07d0
    z = cmplx(om,eta2,8); z2 = cmplx(2*om,2*eta2,8)     ! paper convention (what the module now uses)
    err = 0d0
    do nn=1,nex; do nnp=1,nex; do a=1,3; do b=1,3; do c=1,3
      call get_shg_kernel_ex_static(a,b,c,nn,nnp,s1A,s2A,s3A)
      call get_shg_kernel_ex_static_methodB(a,b,c,nn,nnp,s1B,s2B,s3B)
      call get_shg_kernel_ex_freq(eta2,om,nn,nnp,s1A,s2A,s3A,kA)
      call get_shg_kernel_ex_freq_methodB(eta2,om,nn,nnp,s1B,s2B,s3B,kB)
      d3 = 1d0/((z+e_ex(nn))*(z-e_ex(nnp)))
      pred = -iu*s1B/(z-e_ex(nnp)) - iu*s2B/(z+e_ex(nnp)) &
           + iu*s3B*(1d0/(z-e_ex(nnp)) + 1d0/(z+e_ex(nn)) + (z2-2d0*z)*d3)
      err = max(err, abs(kA-kB-pred))
    end do; end do; end do; end do; end do
    call report('per-pair identity kA-kB = closed form (all a,b,c,n,m)', err, 2, eta2)
  end block

  if (nfail == 0) then
    print '(/,A)', 'ALL TESTS PASSED'
  else
    print '(/,A,I0,A)', 'FAILED: ', nfail, ' check(s)'; stop 1
  end if

contains
  subroutine report(label, val, icplx, e)
    character(*), intent(in) :: label
    real(8), intent(in) :: val, e
    integer, intent(in) :: icplx
    character(len=9) :: tag
    if (icplx==0) tag = 'real     '
    if (icplx==1) tag = 'complex  '
    if (icplx==2) tag = 'generic  '
    if (val < tol) then
      print '(A,1X,A9,A,F6.3,A,ES9.2,A,A)', '  PASS ', tag, ' eta=', e, '  err=', val, '   ', label
    else
      print '(A,1X,A9,A,F6.3,A,ES9.2,A,A)', '  FAIL ', tag, ' eta=', e, '  err=', val, '   ', label
      nfail = nfail + 1
    end if
  end subroutine
  subroutine set_velocity_from_position()
    integer :: a,i,j
    do a=1,3; do i=1,nex; vme_ex(a,i) = -iu*e_ex(i)*xme_ex(a,i); end do; end do
    do a=1,3; do i=1,nex; do j=1,nex
      vme_ex_inter(a,i,j) = iu*(e_ex(i)-e_ex(j))*xme_ex_inter(a,i,j)
    end do; end do; end do
  end subroutine
  subroutine sym_bc(s, ss)
    complex(8), intent(in)  :: s(3,3,3,nw)
    complex(8), intent(out) :: ss(3,3,3,nw)
    integer :: i,j,k
    do i=1,3; do j=1,3; do k=1,3; ss(i,j,k,:) = 0.5d0*(s(i,j,k,:)+s(i,k,j,:)); end do; end do; end do
  end subroutine
  subroutine B_direct(wp, eta2, sB)
    real(8), intent(in) :: wp(nw), eta2
    complex(8), intent(out) :: sB(3,3,3,nw)
    integer :: nj,njp,njpp,nn,nnp,iw
    complex(8) :: s1,s2,s3,kk
    sB = (0d0,0d0)
    do nn=1,norb_ex_cut; do nnp=1,norb_ex_cut
      do nj=1,3; do njp=1,3; do njpp=1,3
        call get_shg_kernel_ex_static_methodB(nj,njp,njpp,nn,nnp,s1,s2,s3)
        do iw=1,nw
          call get_shg_kernel_ex_freq_methodB(eta2,wp(iw),nn,nnp,s1,s2,s3,kk)
          sB(nj,njp,njpp,iw) = sB(nj,njp,njpp,iw) + kk/(dble(npointstotal)*vcell)   ! same sign as driver
        end do
      end do; end do; end do
    end do; end do
  end subroutine
end program
