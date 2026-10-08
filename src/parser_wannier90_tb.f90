module parser_wannier90_tb
   use parser_input_file, only:read_line_numbers_int, iflag_orthonormal
   implicit none
   private
   public :: norb,nR,n1,n2,n3,nRvec
   public :: R,shop,hhop,rhop_c
   public :: material_name
   public :: wannier90_get

   integer norb
   integer nR
   integer n1,n2,n3
   integer nRvec

   real(8) R
   complex*16 hhop,rhop_c,shop
   character(1000) material_name

   dimension   R(3,3)
   integer, allocatable :: Degen(:)
   allocatable nRvec(:,:)
   allocatable hhop(:,:,:)
   allocatable rhop_c(:,:,:,:)
   allocatable shop(:,:,:)

contains

!--------------------------------------------------------------------
   pure function to_lower(str) result(out)
!> Convert a string to lower‑case (portable; no compiler extension)
      character(len=*), intent(in) :: str
      character(len=len(str))      :: out
      integer :: i
      out = str
      do i = 1, len(str)
         if (out(i:i) >= 'A' .and. out(i:i) <= 'Z') &
            out(i:i) = char(iachar(out(i:i)) + 32)
         end do
      end function to_lower
      !--------------------------------------------------------------------
      subroutine wannier90_get(material_name_in)
         implicit none
         ! -----------------------------------------------------------------
         character(len=*), intent(in) :: material_name_in    ! full path read from input.txt
         ! -----------------------------------------------------------------
         integer          :: fp, iR, i, ialpha, ialphap
         integer          :: nkk1, nkk2, nRzero
         real(8)          :: a1,a2,a3,a4,a5,a6
         character(len=:), allocatable :: file2open, basename
         integer          :: p, ext_pos
         integer :: num_chunks
         integer :: rem
         ! -----------------------------------------------------------------
         write(*,*) '2. Entering parser_wannier90_tb'

      ! === 1.  Use the path exactly as supplied ========================
      file2open = trim(material_name_in)

      ! === 2.  Derive clean material name (no dir, no _tb.dat) =========
      p = max( index(file2open,'/',back=.true.),  &
         index(file2open,'\',back=.true.) )   ! works on Win/Linux
      if (p == 0) then
         basename = file2open
      else
         basename = file2open(p+1:)
      end if

      ext_pos = len_trim(basename) - len('_tb.dat') + 1
      if ( ext_pos > 0 .and. to_lower(basename(ext_pos:)) == '_tb.dat' ) then
         basename = basename(:ext_pos-1)
      end if
      material_name = adjustl(basename)

      ! === 3.  Open Wannier90 TB file ==================================
      open(unit=fp, file=file2open, action='read', status='old')
      read(fp,*)
      read(fp,*) R(1,1),R(1,2),R(1,3)
      read(fp,*) R(2,1),R(2,2),R(2,3)
      read(fp,*) R(3,1),R(3,2),R(3,3)
      read(fp,*) norb
      read(fp,*) nR

      !allocate nR, h-,r-, and s-hoppings
      allocate (nRvec(nR,3))
      allocate (hhop(nR,norb,norb))
      allocate (shop(nR,norb,norb))
      allocate (rhop_c(3,nR,norb,norb))
      allocate (Degen(nR))

      num_chunks = nR / 15
      rem = MOD(nR, 15)
      do i = 1, num_chunks
         read(fp, *) Degen((i - 1) * 15 + 1:(i - 1) * 15 + 15)
      end do
      if (rem > 0) then
         read(fp, *) Degen(num_chunks * 15 + 1:num_chunks * 15 + rem)
      end if

      !get the hopping matrices
      do iR=1,nR
         read(fp,*) nRvec(iR,:)
         do ialphap=1,norb
            do ialpha=1,norb
               read(fp,*) nkk1,nkk2,a1,a2
               hhop(iR,nkk1,nkk2)=complex(a1,a2)/Degen(iR)
            end do
         end do
         read(fp,*)
      end do
      
      !locate the (0,0,0) element of nRvec
      nRzero=0
      do iR=1,nR
         if (nRvec(iR,1)==0 .and. nRvec(iR,2)==0 .and. nRvec(iR,3)==0) then
            nRzero=iR
            exit
         end if
      end do
      
      if (nRzero==0) then
         write(*,*) 'ERROR (parser_wannier90_tb): no R=(0,0,0) cell found in the _tb.dat file.'
         write(*,*) '       A physical tight-binding Hamiltonian must include the on-site cell.'
         stop 1
      end if
      

      !get rhoppings
      do iR=1,nR
         read(fp,*) !nRvec is already strored
         do ialphap=1,norb
            do ialpha=1,norb
               read(fp,*) nkk1,nkk2,a1,a2,a3,a4,a5,a6
               rhop_c(1,iR,nkk1,nkk2)=complex(a1,a2)
               rhop_c(2,iR,nkk1,nkk2)=complex(a3,a4)
               rhop_c(3,iR,nkk1,nkk2)=complex(a5,a6)
            end do
         end do
         ! Blank separator after each position block. A plain Wannier90 _tb.dat ends right after the LAST
         ! block with no blank line (every _tb.dat in this repo does), so reading one there hits EOF; it is
         ! only present when an overlap section follows (Orthonormal = false).
         if (iR /= nR .or. .not. iflag_orthonormal) read(fp,*) !blank line
      end do

      !get shoppings 
      if (iflag_orthonormal) then
         shop = 0.0d0
         do ialpha=1,norb
            shop(nRzero,ialpha,ialpha)=1.0d0
         end do
      else 
         do iR=1,nR
            read(fp,*)
            do ialphap=1,norb
               do ialpha=1,norb
                  read(fp,*) nkk1,nkk2,a1,a2
                  shop(iR,nkk1,nkk2)=complex(a1,a2)
               end do
            end do
            if (iR /= nR) read(fp,*)   ! no blank line is required after the last overlap block
         end do
      end if
      close(fp)

      !get orthogonal overlap: this variable is a reminiscent
      !of the interface with the original crystal interface.
      !I maintain the overlap matrix though

!       !wannier functions are orthonormal
!       shop=0.0d0
!       do ialpha=1,norb
!          shop(nRzero,ialpha,ialpha)=1.0d0
!       end do


      !APPLY BIAS BY HAND
      !do iR=1,nR
      !do ialpha=1,norb
      !do ialphap=1,norb
      !hhop(iR,ialpha,ialphap)=hhop(iR,ialpha,ialphap)-0.02d0*rhop_c(3,iR,ialpha,ialphap)
      !end do
      !end do
      !end do
      !do ialpha=1,norb
      !hhop(nRzero,ialpha,ialpha)=hhop(nRzero,ialpha,ialpha)-0.1d0*rhop_c(3,nRzero,ialpha,ialpha)
      !end do

      ! Hermiticity check and repair of H, S and r, in file units, before the unit conversion
      call hermitise_tb(.not. iflag_orthonormal)

      !convert units: to Hartree and bohrs
      hhop=hhop/27.211385d0
      rhop_c=rhop_c/0.52917721067121d0
      R=R/0.52917721067121d0
      write(*,*) '   Wannier hamiltonian has been read'
   end subroutine wannier90_get


   !--------------------------------------------------------------------
   !> Hermiticity of the tight-binding blocks. A physical model has H_ij(R) = conj(H_ji(-R)), S likewise,
   !! and, translating bra and ket, <0i|r|Rj> = conj(<0j|r|-R,i>) + R S_ij(R) (the R S term vanishes for orthonormal
   !! Wannier functions). get_vme_kernels_ome reads only the LOWER triangle and fills the upper one by conjugation, so a
   !! block that breaks this relation would be replaced, silently, by an orbital-ORDER-dependent Hermitian completion that
   !! is not covariant under the crystal symmetries (MoS2: the symmetry-forbidden injection current x3.3).
   !! Here, in file units (eV, Angstrom) and before anything uses the blocks:
   !!  * a block stored as its lower triangle only (strict upper triangle identically zero, e.g. SnTe's H and S) is
   !!    completed by Hermiticity -- never averaged, which would halve its off-diagonal elements;
   !!    (a block whose strict upper triangle is identically zero cannot be told apart from such storage);
   !!  * otherwise, if the defect exceeds herm_rtol * max|X|, a WARNING is printed and X is replaced by its Hermitian part
   !!    [X(R) + X(-R)^+]/2 (positions: [r(R) + r(-R)^+ + R S(R)]/2);
   !!  * a block already Hermitian to herm_rtol is left bit-identical.
   subroutine hermitise_tb(nonorth)
      logical, intent(in) :: nonorth
      integer :: imR(nR), iR, jR, c
      complex*16, allocatable :: x4(:,:,:,:), add(:,:,:,:)

      imR = 0
      do iR = 1, nR
         do jR = 1, nR
            if (all(nRvec(jR,:) == -nRvec(iR,:))) then
               imR(iR) = jR
               exit
            end if
         end do
      end do
      if (any(imR == 0)) write(*,'(a,i0,a)') '    WARNING (parser_wannier90_tb): ', count(imR == 0), &
         ' cell(s) R without a -R partner; their Hermiticity cannot be checked or repaired.'

      allocate(x4(1,nR,norb,norb))
      x4(1,:,:,:) = hhop
      call hermitian_fix('H(R)', 'eV', x4, imR)
      hhop = x4(1,:,:,:)
      if (nonorth) then
         x4(1,:,:,:) = shop
         call hermitian_fix('S(R)', '', x4, imR)
         shop = x4(1,:,:,:)
      end if
      deallocate(x4)

      allocate(add(3,nR,norb,norb))
      do iR = 1, nR
         do c = 1, 3
            add(c,iR,:,:) = dot_product(dble(nRvec(iR,:)), R(:,c)) * shop(iR,:,:)
         end do
      end do
      if (nonorth) then
         call hermitian_fix('position matrices r(R)', 'Angstrom', rhop_c, imR, add)
      else
         call hermitian_fix('position matrices r(R)', 'Angstrom', rhop_c, imR)   ! R S(R) = 0 here
      end if
      deallocate(add)
   end subroutine hermitise_tb

   !> See hermitise_tb. X(c,R,i,j); the relation is X(c,R,i,j) = conj(X(c,-R,j,i)) + add(c,R,i,j).
   subroutine hermitian_fix(label, units, X, imR, add)
      character(len=*), intent(in)       :: label, units
      complex*16, intent(inout)          :: X(:,:,:,:)
      integer, intent(in)                :: imR(:)
      complex*16, intent(in), optional   :: add(:,:,:,:)
      real(8), parameter :: herm_rtol = 1.0d-8
      character(len=1), parameter :: comp(3) = (/'x','y','z'/)
      complex*16, allocatable :: T(:,:,:,:)
      integer :: nc, iR, i, j, w(4)
      real(8) :: xmax, dmax
      logical :: upper_zero, lower_nonzero

      nc = size(X,1)
      xmax = maxval(abs(X))
      if (xmax == 0.0d0) return
      upper_zero = .true.; lower_nonzero = .false.
      do j = 1, norb
         do i = 1, norb
            if (i < j .and. any(X(:,:,i,j) /= (0.0d0,0.0d0))) upper_zero = .false.
            if (i > j .and. any(X(:,:,i,j) /= (0.0d0,0.0d0))) lower_nonzero = .true.
         end do
      end do
      allocate(T(nc,nR,norb,norb))
      if (upper_zero .and. lower_nonzero) then
         call build_target()
         do j = 1, norb
            do i = 1, j - 1
               X(:,:,i,j) = T(:,:,i,j)
            end do
         end do
         write(*,'(4x,a)') label//' is stored as its lower triangle only: upper triangle completed by Hermiticity.'
      end if
      call build_target()
      dmax = maxval(abs(X - T))
      if (dmax > herm_rtol*xmax) then
         w = maxloc(abs(X - T))
         write(*,'(4x,a,es10.3,1x,a,a,es9.2,a,3i4,a,2i5)') 'WARNING (parser_wannier90_tb): '//label// &
            ' not Hermitian: max defect ', dmax, trim(units), ' (', dmax/xmax, ' of max|X|) at R =', nRvec(w(2),:), &
            ', orbitals', w(3), w(4)
         if (nc == 3) write(*,'(10x,a)') 'component '//comp(w(1))
         if (present(add)) then
            write(*,'(10x,a)') 'replaced by its Hermitian part [r(R) + r(-R)^+ + R S(R)]/2.'
         else
            write(*,'(10x,a)') 'replaced by its Hermitian part [X(R) + X(-R)^+]/2.'
         end if
         write(*,'(10x,a)') 'Fix the model file to silence this warning.'
         X = 0.5d0*(X + T)
      else
         write(*,'(4x,a,es9.2,1x,a)') label//' Hermitian to ', dmax, trim(units)
      end if
      deallocate(T)
   contains
      subroutine build_target()
         integer :: ii, jj, kR
         do kR = 1, nR
            if (imR(kR) == 0) then
               T(:,kR,:,:) = X(:,kR,:,:)
               cycle
            end if
            do jj = 1, norb
               do ii = 1, norb
                  T(:,kR,ii,jj) = conjg(X(:,imR(kR),jj,ii))
                  if (present(add)) T(:,kR,ii,jj) = T(:,kR,ii,jj) + add(:,kR,ii,jj)
               end do
            end do
         end do
      end subroutine build_target
   end subroutine hermitian_fix

end module parser_wannier90_tb
