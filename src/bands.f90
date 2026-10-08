module bands
  use constants_math
  use parser_wannier90_tb, &
  only:material_name,nR,nRvec,norb,R,shop,hhop
  use parser_optics_xatu_dim, &
  only:G
  use parser_input_file, only: iflag_orthonormal, nf, nband_index, &
       kpath_nv, kpath_frac, kpath_count, kpath_labels, kpath_has_labels
  implicit none

  ! the band-structure path: Cartesian k (bohr^-1), accumulated length, and the vertices for the header
  real(8), allocatable :: rkxvector_path(:), rkyvector_path(:), rkzvector_path(:), rklengthvector_path(:)
  integer :: nv_path = 0
  real(8), allocatable :: vfrac(:,:)
  integer, allocatable :: vcount(:), vrow(:)
  character(len=16), allocatable :: vlabel(:)

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Band structure along a path, written to bands_<material>.dat. The path is the Kpath block of the
  !! input (vertices in reduced coordinates along the reciprocal lattice vectors, each with the number of
  !! points to the next vertex; 1 = the vertex alone, then a jump), or, without one, a default for the
  !! lattice: Gamma-M-K-Gamma for a hexagonal 2D lattice, Gamma-X-S-Y-Gamma for any other 2D lattice
  !! (plus -Z in 3D), Gamma-X in 1D.
  !! @return void
  subroutine get_energy_bands()
    implicit none
    if (kpath_nv > 0) then
      nv_path = kpath_nv
      allocate(vfrac(3,nv_path), vcount(nv_path), vlabel(nv_path), vrow(nv_path))
      vfrac = kpath_frac; vcount = kpath_count
      call default_labels()
      if (kpath_has_labels) then
        if (size(kpath_labels) == nv_path) then
          vlabel = kpath_labels
        else
          write(*,'(a,i0,a,i0,a)') '    WARNING (bands): Kpath_labels has ', size(kpath_labels), &
               ' labels for ', nv_path, ' vertices; labels ignored.'
        end if
      end if
    else
      call default_path()
    end if
    call get_path()
    call get_eigenenergies()
  end subroutine get_energy_bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine default_labels()
    implicit none
    integer :: i
    do i = 1, nv_path
      write(vlabel(i),'(a,i0)') 'P', i
    end do
  end subroutine default_labels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Default path from the active reciprocal vectors (non-zero rows of G), about 300 points in all,
  !! distributed in proportion to the segment lengths.
  subroutine default_path()
    implicit none
    integer :: act(3), na, i, a1, a2, a3
    real(8) :: b1(3), b2(3), cosang, kc(3,3), lk(3), lb, vtx(3,6), seg, total
    character(len=16) :: lab(6)
    integer :: nvtx
    na = 0
    do i = 1, 3
      if (sum(abs(G(i,:))) > 1.0d-12) then
        na = na + 1; act(na) = i
      end if
    end do
    vtx = 0.0d0
    if (na == 1) then
      a1 = act(1)
      nvtx = 2; lab(1:2) = ['G', 'X']
      vtx(a1,2) = 0.5d0
    else if (na == 2) then
      a1 = act(1); a2 = act(2)
      b1 = G(a1,:); b2 = G(a2,:)
      cosang = dot_product(b1,b2)/(norm2(b1)*norm2(b2))
      if (abs(norm2(b1) - norm2(b2)) < 1.0d-3*norm2(b1) .and. abs(abs(cosang) - 0.5d0) < 1.0d-3) then
        ! hexagonal: K is the zone corner, the candidate of length |b|/sqrt(3)
        kc(:,1) = (b1 + b2)/3.0d0; kc(:,2) = (2.0d0*b1 + b2)/3.0d0; kc(:,3) = (b1 + 2.0d0*b2)/3.0d0
        lb = norm2(b1)/sqrt(3.0d0)
        do i = 1, 3
          lk(i) = abs(norm2(kc(:,i)) - lb)
        end do
        nvtx = 4; lab(1:4) = ['G', 'M', 'K', 'G']
        vtx(a1,2) = 0.5d0
        select case (minloc(lk, dim=1))
          case (1); vtx(a1,3) = 1.0d0/3.0d0; vtx(a2,3) = 1.0d0/3.0d0
          case (2); vtx(a1,3) = 2.0d0/3.0d0; vtx(a2,3) = 1.0d0/3.0d0
          case (3); vtx(a1,3) = 1.0d0/3.0d0; vtx(a2,3) = 2.0d0/3.0d0
        end select
      else
        nvtx = 5; lab(1:5) = ['G', 'X', 'S', 'Y', 'G']
        vtx(a1,2) = 0.5d0
        vtx(a1,3) = 0.5d0; vtx(a2,3) = 0.5d0
        vtx(a2,4) = 0.5d0
      end if
    else
      a1 = act(1); a2 = act(2); a3 = act(3)
      nvtx = 6; lab(1:6) = ['G', 'X', 'S', 'Y', 'G', 'Z']
      vtx(a1,2) = 0.5d0
      vtx(a1,3) = 0.5d0; vtx(a2,3) = 0.5d0
      vtx(a2,4) = 0.5d0
      vtx(a3,6) = 0.5d0
    end if
    nv_path = nvtx
    allocate(vfrac(3,nv_path), vcount(nv_path), vlabel(nv_path), vrow(nv_path))
    vfrac = vtx(:,1:nvtx); vlabel = lab(1:nvtx)
    total = 0.0d0
    do i = 1, nvtx - 1
      total = total + norm2(matmul(vfrac(:,i+1) - vfrac(:,i), G))
    end do
    do i = 1, nvtx - 1
      seg = norm2(matmul(vfrac(:,i+1) - vfrac(:,i), G))
      vcount(i) = max(2, nint(300.0d0*seg/total))
    end do
    vcount(nvtx) = 1
  end subroutine default_path
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Points of the path: vcount(i) >= 2 gives vcount(i) points from vertex i towards vertex i+1 (vertex
  !! i+1 itself starts the next segment); vcount(i) = 1 gives vertex i alone and a jump to vertex i+1,
  !! across which the path length does not advance; the last vertex is included if its count is >= 1.
  subroutine get_path()
    implicit none
    integer :: i, j, np, idx
    real(8) :: k(3), kprev(3), s
    logical :: jump
    np = 0
    do i = 1, nv_path - 1
      np = np + vcount(i)
    end do
    if (vcount(nv_path) >= 1) np = np + 1
    allocate(rkxvector_path(np), rkyvector_path(np), rkzvector_path(np), rklengthvector_path(np))
    idx = 0; s = 0.0d0; jump = .false.
    do i = 1, nv_path
      if (i == nv_path) then
        if (vcount(i) < 1) exit
        call add_point(vfrac(:,i))
        vrow(i) = idx
        exit
      end if
      vrow(i) = idx + 1
      do j = 0, vcount(i) - 1
        call add_point(vfrac(:,i) + dble(j)/dble(vcount(i))*(vfrac(:,i+1) - vfrac(:,i)))
      end do
      jump = (vcount(i) == 1)
    end do
  contains
    subroutine add_point(f)
      real(8), intent(in) :: f(3)
      k = matmul(f, G)                         ! k = sum_i f_i G(i,:)
      if (idx > 0 .and. .not. jump) s = s + norm2(k - kprev)
      jump = .false.
      idx = idx + 1
      rkxvector_path(idx) = k(1); rkyvector_path(idx) = k(2); rkzvector_path(idx) = k(3)
      rklengthvector_path(idx) = s
      kprev = k
    end subroutine add_point
  end subroutine get_path
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_eigenenergies()
    implicit none
    integer :: iounit10, ialpha, ialphap, iRp, ibz, j, np
    real(8) :: e(norb), rkx, rky, rkz
    real(8), allocatable :: eall(:,:), Rxyz(:,:)
    complex(8), allocatable :: f(:)
    complex(8) :: skernel(norb,norb), hkernel(norb,norb), s_work(norb,norb)
    real(8) :: vbm, cbm, svbm, scbm, dgap, sdgap
    character(len=2048) :: line

    np = size(rkxvector_path)
    allocate(eall(norb, np), f(nR), Rxyz(nR,3))
    do iRp = 1, nR
      do j = 1, 3
        Rxyz(iRp,j) = dble(nRvec(iRp,1))*R(1,j)+dble(nRvec(iRp,2))*R(2,j)+dble(nRvec(iRp,3))*R(3,j)
      end do
    end do
    do ibz = 1, np
      rkx = rkxvector_path(ibz); rky = rkyvector_path(ibz); rkz = rkzvector_path(ibz)
      ! Bloch sums as one product over R: [H S](k) = f^T [H(R) S(R)], f = e^{i k.R}; the lower triangle is
      ! kept and the upper one completed by conjugation, as in the matrix-element kernels
      do iRp = 1, nR
        f(iRp) = exp(cmplx(0.0d0, rkx*Rxyz(iRp,1)+rky*Rxyz(iRp,2)+rkz*Rxyz(iRp,3), 8))
      end do
      call zgemv('T', nR, norb*norb, (1.0d0,0.0d0), hhop, nR, f, 1, (0.0d0,0.0d0), hkernel, 1)
      call zgemv('T', nR, norb*norb, (1.0d0,0.0d0), shop, nR, f, 1, (0.0d0,0.0d0), skernel, 1)
      do ialpha = 1, norb
        do ialphap = 1, ialpha - 1
          hkernel(ialphap,ialpha) = conjg(hkernel(ialpha,ialphap))
          skernel(ialphap,ialpha) = conjg(skernel(ialpha,ialphap))
        end do
      end do
      if (iflag_orthonormal) then
        call diagoz(norb,e,hkernel)
      else
        s_work = skernel            ! zhegv overwrites its S argument
        call diagoz_gen(norb,e,hkernel,s_work)
      end if
      eall(:,ibz) = e*27.211385d0
    end do

    open(newunit=iounit10,file='bands_'//trim(material_name)//'.dat')
    write(iounit10,'(a)') '# bands along a k-path. Columns: kx ky kz (bohr^-1, Cartesian), s (path length, bohr^-1),'
    write(iounit10,'(a)') '# E_1 ... E_norb (eV, eigenvalues of the Wannier model, not shifted)'
    write(iounit10,'(a)') '# vertex  label  row  s  (reduced coordinates along the reciprocal lattice vectors)'
    do j = 1, nv_path
      if (j < nv_path .or. vcount(nv_path) >= 1) then
        write(iounit10,'(a,i4,2x,a,i7,f14.8,3f12.7)') '# vertex', j, vlabel(j), vrow(j), &
             rklengthvector_path(vrow(j)), vfrac(:,j)
      end if
    end do
    if (allocated(nband_index)) then
      write(line,'(a,i0,a)') '# Nfermi = ', nf, '; Bandlist window bands:'
      do j = 1, size(nband_index)
        write(line,'(a,1x,i0)') trim(line), nband_index(j)
      end do
      write(iounit10,'(a)') trim(line)
    end if
    if (nf >= 1 .and. nf < norb) then
      vbm = maxval(eall(nf,:)); cbm = minval(eall(nf+1,:))
      svbm = rklengthvector_path(maxloc(eall(nf,:),dim=1)); scbm = rklengthvector_path(minloc(eall(nf+1,:),dim=1))
      dgap = minval(eall(nf+1,:) - eall(nf,:)); sdgap = rklengthvector_path(minloc(eall(nf+1,:) - eall(nf,:),dim=1))
      write(iounit10,'(a,f12.6,a,f12.6,a,f12.6,a,f12.6)') '# along the path: VBM (band Nfermi) ', vbm, ' eV at s = ', &
           svbm, '; CBM ', cbm, ' eV at s = ', scbm
      write(iounit10,'(a,f12.6,a,f12.6,a,f12.6)') '# along the path: gap ', cbm - vbm, ' eV; smallest direct gap ', &
           dgap, ' eV at s = ', sdgap
    end if
    do ibz = 1, np
      write(iounit10,*) rkxvector_path(ibz), rkyvector_path(ibz), rkzvector_path(ibz), &
           rklengthvector_path(ibz), (eall(j,ibz), j=1,norb)
    end do
    close(iounit10)
    write(*,'(a,i0,a,i0,a)') '    Band structure: ', np, ' k-points along ', nv_path, ' vertices written to bands_' &
         //trim(material_name)//'.dat'
  end subroutine get_eigenenergies

end module bands
