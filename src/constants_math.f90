module constants_math
  implicit none
  real(8), parameter :: pi=3.14159265358979323846d0
  real(8), parameter :: dk=1.0d-6
  
  !real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
  !private :: ax,ay,az,bx,by,bz,cx,cy,cz
  !public :: ax,ay,az,bx,by,bz,cx,cy,cz
  contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine percentage_index(kacum,ktotal,kmoment)
    implicit none
    integer :: kacum,ktotal,kmoment
    integer :: npercentage,nrest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    npercentage=int(dble(kacum)/dble(ktotal)*100.0d0)
    nrest=mod(npercentage,10)
    !write(*,*) 'Percentage of k-points read:',npercentage,' %'
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
!   subroutine diagoz(n,w,h)
!     !finding the eigenvalues of a complex matrix using LAPACK
!     implicit real*8 (a-h,o-z)
!     !declarations, notice double precision
!     integer n,INFO,LWORK
!     dimension w(n)
!     dimension RWORK(3*n-2)
!     dimension h(n,n)
! 
!     real(8) w
!     complex*16 h
!     complex*16 WORK(2*n)
!     character*1 JOBZ,UPLO
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !find the solution using the LAPACK routine ZGEEV
!     JOBZ='V'
!     UPLO='U'
!     LWORK=2*n
!             
!     call zheev(JOBZ, UPLO, n, h, n, w, WORK, LWORK, RWORK, INFO)
!   end
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

    ! Query optimal LWORK first
    call zheev(JOBZ, UPLO, n, h, n, w, WORK_QUERY, -1, RWORK, INFO)
    LWORK = max(2*n, int(dble(WORK_QUERY(1))))
    allocate(WORK(LWORK))

    call zheev(JOBZ, UPLO, n, h, n, w, WORK, LWORK, RWORK, INFO)
    deallocate(WORK)
  end subroutine diagoz
  
end module constants_math

