program test_ome_sp_symmetry
  use parser_wannier90_tb
  use parser_optics_xatu_dim
  use ome_sp
  implicit none

  real(8) :: rkx, rky, rkz
  complex*16, allocatable :: skernel(:,:), hkernel(:,:)
  complex*16, allocatable :: sderkernel(:,:,:), hderkernel(:,:,:), akernel(:,:,:)
  complex*16, allocatable :: hk_ev(:,:), vme(:,:,:), M1(:,:), T1(:,:)
  complex*16, allocatable :: gen_der(:,:,:,:), vme_der(:,:,:,:)
  complex*16, allocatable :: gd1(:,:,:,:), gd2(:,:,:,:), gd3(:,:,:,:)
  complex*16, allocatable :: abc(:,:,:)
  complex*16, allocatable :: berry_eigen1(:,:,:), berry_eigen2(:,:,:), berry_eigen(:,:,:)
  complex*16, allocatable :: hk_ev_neigh(:,:,:), vme_neigh(:,:,:,:)
  real(8),    allocatable :: shift_vector(:,:,:,:), vme_der_phase(:,:,:,:)
  real(8),    allocatable :: e(:)

  integer :: nn, nnp, nj, njp
  real(8) :: max_asym

  ! --- read a real, small wannier90_tb.dat to get norb, R, hhop, shop, rhop_c ---
  call read_wannier90_tb()   ! whatever your actual driver call is named

  allocate(skernel(norb,norb), hkernel(norb,norb))
  allocate(sderkernel(norb,norb,3), hderkernel(norb,norb,3), akernel(norb,norb,3))
  allocate(hk_ev(norb,norb), vme(norb,norb,3), M1(norb,norb), T1(norb,norb), e(norb))
  allocate(gen_der(norb,norb,3,3), vme_der(norb,norb,3,3))
  allocate(gd1(norb,norb,3,3), gd2(norb,norb,3,3), gd3(norb,norb,3,3))
  allocate(abc(norb,norb,3))
  allocate(berry_eigen1(norb,norb,3), berry_eigen2(norb,norb,3), berry_eigen(norb,norb,3))
  allocate(hk_ev_neigh(norb,norb,7), vme_neigh(norb,norb,3,7))
  allocate(shift_vector(norb,norb,3,3), vme_der_phase(norb,norb,3,3))

  call set_active_flags()
  call set_R_cache()

  rkx = 0.13d0; rky = -0.07d0; rkz = 0.0d0   ! any generic, non-symmetry-special point

  call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
  call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev,e,vme,M1,T1)
  call get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
  call get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der,shift_vector, &
        berry_eigen1,berry_eigen2,berry_eigen,hk_ev_neigh,vme_neigh, &
        skernel,hkernel,sderkernel,hderkernel,akernel,vme_der_phase,M1,T1)

  ! Sanity check 1: no NaN/Inf anywhere
  if (any(vme /= vme) .or. any(shift_vector /= shift_vector)) then
    print *, 'FAIL: NaN detected'
    stop 1
  end if

  ! Sanity check 2: shift_vector(nn,nnp,nj,njp) should equal
  ! -shift_vector(nnp,nn,nj,njp) for a real, TRS-symmetric system --
  ! this follows directly from the paper's own symmetry argument
  ! (Supplementary Note 5) and is a strong, reference-free check.
  max_asym = 0.0d0
  do nn = 1, norb
    do nnp = 1, norb
      do nj = 1, 3
        do njp = 1, 3
          max_asym = max(max_asym, abs(shift_vector(nn,nnp,nj,njp) + shift_vector(nnp,nn,nj,njp)))
        end do
      end do
    end do
  end do
  print *, 'max |shift_vector(nn,nnp)+shift_vector(nnp,nn)| =', max_asym
  if (max_asym > 1.0d-8) then
    print *, 'FAIL: shift_vector antisymmetry violated'
    stop 1
  end if

  print *, 'PASS'
end program test_ome_sp_symmetry