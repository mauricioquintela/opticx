module constants_math
  implicit none
  real(8), parameter :: pi=3.14159265358979323846d0
  real(8), parameter :: dk=1.0d-6
  ! Second-order conductivity, atomic units -> uA nm / V^2. Was written out by hand at five
  ! sites in sigma_second_sp/ex; one of them had lost the d0 on the Hartree value and so ran
  ! 2.0e-8 relative high (audit 2026-09-30).
  ! Electron charge convention: e = -|e| (2026-10-04). Every second-order
  ! conductivity is proportional to e^3 (npj Comput. Mater. 11, 13, Eqs. 9-11; Taghizadeh et al.
  ! PRB 96, 195413, Eq. A3a) and all second-order kernels are written with e = 1, so the sign of
  ! e enters here, once, for all six second-order output files (shift, SHG, general; sp and ex).
  ! With e = -|e| opticx reproduces the published MoS2 and GeS shift spectra of the npj paper,
  ! which were computed with this convention. First-order responses carry e^2 and are unaffected.
  real(8), parameter :: e_charge_au = -1.0d0
  real(8), parameter :: sigma2_au_to_si = e_charge_au**3 * &
       (6.623618d-03)*(1.0d+06)*(27.211386d0**(-2))*(5.291772d-11)*(1.0d+09)
  
  contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine percentage_index(kacum,ktotal,kmoment)
    implicit none
    integer :: kacum,ktotal,kmoment
    integer :: npercentage,nrest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    npercentage=int(dble(kacum)/dble(ktotal)*100.0d0)
    nrest=mod(npercentage,10)
    if (nrest.eq.0) then
      if (kmoment.ne.npercentage) then
        write(*,*) '   Percentage of loop:',npercentage,' %'
      end if
      kmoment=npercentage
    end if   
  end subroutine percentage_index
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine crossproduct(ax,ay,az,bx,by,bz,cx,cy,cz)
    implicit none
    real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cx=ay*bz-az*by
    cy=az*bx-ax*bz
    cz=ax*by-ay*bx     
  end subroutine crossproduct


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NAME:         diagoz
!   INPUTS:       h matrix to diagonalize
!                 n dimension of h
!   OUTPUTS:      w; eigenvalues of h
!                 h;  gives eigenvectors by columns as output
!   DESCRIPTION:  this subroutine uses Lapack libraries to diagonalize
!                 an hermitian complex matrix.
!   
!     Juan Jose Esteve-Paredes                28.11.2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!> Diagonalises a complex Hermitian matrix in place (LAPACK zheev, 'V','U').
!! Used for H(k) throughout: the eigenvector phase it returns is arbitrary and
!! k-dependent, which is why anything differentiating in k must fix the gauge first
!! (see get_berry_eigen_fourpoint in ome_sp.f90).
!! @param n  Matrix dimension.
!! @param w  On exit, the n eigenvalues in ascending order.
!! @param h  On entry the matrix; on exit its eigenvectors, column n holding |psi_n>.
!! @return void
subroutine diagoz(n,w,h)
  implicit none
  integer, intent(in) :: n
  real(8), intent(out) :: w(n)
  complex*16, intent(inout) :: h(n,n)
  integer :: INFO, LWORK
  real(8) :: RWORK(3*n-2)
  complex*16 :: WORK_QUERY(1)
  complex*16, allocatable :: WORK(:)
  character*1 :: JOBZ, UPLO

  JOBZ='V'; UPLO='U'
  call zheev(JOBZ, UPLO, n, h, n, w, WORK_QUERY, -1, RWORK, INFO)
  LWORK = max(2*n, int(dble(WORK_QUERY(1))))
  allocate(WORK(LWORK))
  call zheev(JOBZ, UPLO, n, h, n, w, WORK, LWORK, RWORK, INFO)
  deallocate(WORK)
end subroutine diagoz
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NAME:         diagoz_gen
!   INPUTS:       h: matrix to diagonalize
!                 s: overlap matrix 
!                 n: dimension of h and s
!   OUTPUTS:      w; eigenvalues of h
!                 h;  gives eigenvectors by columns as output
!   DESCRIPTION:  this subroutine solves the generalized eigenvalue problem H*v = e*S*v using 
!                 LAPACK zhegv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine diagoz_gen(n,w,h,s)
    implicit none 
    integer n,INFO,LWORK,ITYPE
    dimension w(n)
    dimension RWORK(3*n-2)
    dimension h(n,n)
    dimension s(n,n)

    real(8) w
    real(8) RWORK
    complex*16 h,s
    complex*16 WORK(2*n)
    character*1 JOBZ,UPLO

    ITYPE=1
    JOBZ='V'
    UPLO='U'
    LWORK=2*n

    call zhegv(ITYPE,JOBZ,UPLO,n,h,n,s,n,w,WORK,LWORK,RWORK,INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: Generalized eigenvalue problem failed. zhegv failed with INFO =', INFO
      write(*,*) '       (INFO > n means the overlap matrix S(k) is not positive definite.)'
      stop 1
    end if
  end
end module constants_math

