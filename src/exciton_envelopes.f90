module exciton_envelopes
  use constants_math
  use parser_optics_xatu_dim, &
    only: G, npointstotal, rkxvector, rkyvector, rkzvector, &
          norb_ex, norb_ex_cut, norb_ex_band, fk_ex
  use parser_input_file, &
    only: ndim
  use parser_wannier90_tb, &
    only: nRvec
  implicit none

  complex(8), allocatable :: fk_ex_der(:,:,:)

  ! PATCH: dimensionality flags computed once instead of being recomputed
  ! (via NORM2 over nRvec) inside get_fk_ex_der_k, get_fk_ex_k_interp, and
  ! get_k_kc separately -- the latter two are called on every finite-
  ! difference neighbour evaluation, i.e. deep inside the hottest loop.
  logical, save :: active_x = .false., active_y = .false., active_z = .false.
  logical, save :: active_flags_set = .false.

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_active_flags()
    implicit none
    integer :: n_active
    active_x = (norm2(real(nRvec(:,1))) /= 0.0d0)
    active_y = (norm2(real(nRvec(:,2))) /= 0.0d0)
    active_z = (norm2(real(nRvec(:,3))) /= 0.0d0)
    active_flags_set = .true.

    ! PATCH: sanity-check ndim against the actual number of periodic
    ! directions instead of trusting them to agree silently.
    n_active = count([active_x, active_y, active_z])
    if (n_active /= ndim) then
      write(*,*) 'WARNING (exciton_envelopes): ndim = ', ndim, &
                 ' but ', n_active, ' lattice direction(s) appear active.', &
                 ' Finite-difference stencil selection may be inconsistent.'
    end if
  end subroutine set_active_flags

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Flat index of the nearest neighbour of point ibz along one grid direction, with
  ! periodic wrap-around. step = index stride of that direction, cnt = number of points
  ! along it, pos = 1-based position of ibz along it, sgn = +1 (next) or -1 (previous).
  pure integer function k_neighbor(ibz, step, cnt, pos, sgn) result(nb)
    integer, intent(in) :: ibz, step, cnt, pos, sgn
    nb = ibz + sgn*step
    if (sgn > 0 .and. pos == cnt) nb = nb - cnt*step
    if (sgn < 0 .and. pos == 1)   nb = nb + cnt*step
  end function k_neighbor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> k-derivative of the exciton envelope functions, d/dk psi^(n)_cvk, by central
  !! differences on the Monkhorst-Pack mesh. Neighbours are taken periodically
  !! (k_neighbor wraps the zone edge); zeroing the edge points instead makes the
  !! inter-exciton position matrix X_nm non-Hermitian by 12.6% on a 30x30 grid.
  !! Feeds the Q-part of X_nm, Eq. (B3b) of Taghizadeh & Pedersen, PRB 97, 205432 (2018).
  !! @return void (fills the module array fk_ex_der)
  subroutine get_fk_ex_der_k()
    use omp_lib
    implicit none

    real(8) :: rk1vector(npointstotal)
    real(8) :: rk2vector(npointstotal)
    real(8) :: rk3vector(npointstotal)

    complex(8), allocatable :: fk_ex_col(:)
    complex(8), allocatable :: fk_der_col(:,:)

    integer  :: ibz, nn, j, j_aux, nj
    real(8)  :: xp_bz, yp_bz, zp_bz, xc_bz, yc_bz, zc_bz
    real(8)  :: dum1,dum2,dum3,dum4,dum5,dum6

    ! --- NEW: grid bookkeeping (replaces per-point get_k_kc/interp1 calls) ---
    integer  :: nside
    real(8)  :: slice1, slice2, slice3          ! grid spacing along each active crystal direction
    real(8)  :: dRkdK(3,3)                      ! dRkdK(i,a) = d(rk_i)/d(k_a), constant, computed once
    complex(8) :: d_rk(3)                       ! central-difference df/d(rk_i), i=1,2,3
    integer  :: i_fast, i_mid, i_slow           ! grid-index position along each active direction

    ! PATCH: compute once, up front, instead of locally on every call.
    if (.not. active_flags_set) call set_active_flags()

    do ibz = 1, npointstotal
      call get_k_kc(G, rkxvector(ibz), rkyvector(ibz), rkzvector(ibz), &
                    rk1vector(ibz), rk2vector(ibz), rk3vector(ibz), &
                    xp_bz, yp_bz, zp_bz, xc_bz, yc_bz, zc_bz)
    end do

    ! ------------------------------------------------------------------
    ! Jacobian d(rk_i)/d(k_a): get_k_kc's Cartesian->crystal map is exactly
    ! linear (no additive offset -- verified from its formulas, which are
    ! homogeneous in xp,yp,zp), so evaluating it at each Cartesian unit
    ! vector gives an exact Jacobian column, with NO finite-difference
    ! epsilon and NO approximation error. This replaces, once per run
    ! instead of millions of times, exactly the coordinate-transform work
    ! that used to be redone inside get_k_kc for every (j,nn,ibz,direction)
    ! finite-difference neighbour.
    ! ------------------------------------------------------------------
    call get_k_kc(G, 1.0d0,0.0d0,0.0d0, dRkdK(1,1),dRkdK(2,1),dRkdK(3,1), &
                  dum1,dum2,dum3,dum4,dum5,dum6)
    call get_k_kc(G, 0.0d0,1.0d0,0.0d0, dRkdK(1,2),dRkdK(2,2),dRkdK(3,2), &
                  dum1,dum2,dum3,dum4,dum5,dum6)
    call get_k_kc(G, 0.0d0,0.0d0,1.0d0, dRkdK(1,3),dRkdK(2,3),dRkdK(3,3), &
                  dum1,dum2,dum3,dum4,dum5,dum6)

    ! ------------------------------------------------------------------
    ! Grid geometry. This assumes the same flat-index convention already
    ! relied upon by get_fk_ex_k_interp elsewhere in this module: crystal
    ! direction 1 is fastest-varying (index step 1), direction 2 is next
    ! (index step nside), direction 3 is slowest (index step nside*nside),
    ! restricted to whichever directions are active for this ndim. If that
    ! convention ever changes (e.g. a different k-mesh generator upstream),
    ! this routine and get_fk_ex_k_interp must be updated together.
    ! ------------------------------------------------------------------
    if (ndim == 1) then
      nside  = npointstotal
      slice1 = 0.0d0; slice2 = 0.0d0; slice3 = 0.0d0
      if (active_x) slice1 = rk1vector(2) - rk1vector(1)
      if (active_y) slice2 = rk2vector(2) - rk2vector(1)
      if (active_z) slice3 = rk3vector(2) - rk3vector(1)
    else if (ndim == 2) then
      nside = nint(sqrt(dble(npointstotal)))
      if (active_x .and. active_y) then
        slice1 = rk1vector(2)     - rk1vector(1)
        slice2 = rk2vector(nside+1) - rk2vector(1)
      else if (active_x .and. active_z) then
        slice1 = rk1vector(2)     - rk1vector(1)
        slice3 = rk3vector(nside+1) - rk3vector(1)
      else
        slice2 = rk2vector(2)     - rk2vector(1)
        slice3 = rk3vector(nside+1) - rk3vector(1)
      end if
    else  ! ndim == 3
      nside  = nint(dble(npointstotal)**(1.0d0/3.0d0))
      slice1 = rk1vector(2)              - rk1vector(1)
      slice2 = rk2vector(nside+1)        - rk2vector(1)
      slice3 = rk3vector(nside*nside+1)  - rk3vector(1)
    end if

    write(*,*) '   Evaluating exciton envelope function derivative with respect to k...'

    !$omp parallel do collapse(2) default(none) schedule(dynamic) &
    !$omp private(j, nn, ibz, nj, j_aux, &
    !$omp         d_rk, i_fast, i_mid, i_slow, &
    !$omp         fk_ex_col, fk_der_col) &
    !$omp   shared(npointstotal, norb_ex_band, norb_ex_cut, nside, ndim, drkdk) &
    !$omp   shared(fk_ex, fk_ex_der, active_x, active_y, active_z, slice1, slice2, slice3)
    do j = 1, norb_ex_band
      do nn = 1, norb_ex_cut

        allocate(fk_ex_col(npointstotal))
        allocate(fk_der_col(3, npointstotal))
        fk_ex_col  = (0.0d0, 0.0d0)
        fk_der_col = (0.0d0, 0.0d0)

        do ibz = 1, npointstotal
          j_aux = (ibz-1)*norb_ex_band + j
          fk_ex_col(ibz) = fk_ex(j_aux, nn)
        end do

        do ibz = 1, npointstotal
          d_rk = (0.0d0, 0.0d0)

          ! Central differences with PERIODIC wrap-around. The mesh covers one full
          ! reciprocal cell and, in the lattice gauge (H(k+G)=H(k)), the smooth-gauge
          ! envelope is periodic, so zone-edge points use their neighbours across the
          ! boundary. (Skipping them left the derivative operator non-antisymmetric,
          ! which made X_nm non-Hermitian.)
          if (ndim == 1) then

            if (active_z) then
              d_rk(3) = (fk_ex_col(k_neighbor(ibz,1,npointstotal,ibz, 1)) &
                       - fk_ex_col(k_neighbor(ibz,1,npointstotal,ibz,-1))) / (2.0d0*slice3)
            else if (active_y) then
              d_rk(2) = (fk_ex_col(k_neighbor(ibz,1,npointstotal,ibz, 1)) &
                       - fk_ex_col(k_neighbor(ibz,1,npointstotal,ibz,-1))) / (2.0d0*slice2)
            else
              d_rk(1) = (fk_ex_col(k_neighbor(ibz,1,npointstotal,ibz, 1)) &
                       - fk_ex_col(k_neighbor(ibz,1,npointstotal,ibz,-1))) / (2.0d0*slice1)
            end if

          else if (ndim == 2) then

            i_fast = mod(ibz-1, nside) + 1          ! step 1
            i_slow = (ibz-1)/nside + 1               ! step nside
            if (active_x .and. active_y) then       ! fast = direction 1, slow = direction 2
              d_rk(1) = (fk_ex_col(k_neighbor(ibz,1,    nside,i_fast, 1)) &
                       - fk_ex_col(k_neighbor(ibz,1,    nside,i_fast,-1))) / (2.0d0*slice1)
              d_rk(2) = (fk_ex_col(k_neighbor(ibz,nside,nside,i_slow, 1)) &
                       - fk_ex_col(k_neighbor(ibz,nside,nside,i_slow,-1))) / (2.0d0*slice2)
            else if (active_x .and. active_z) then  ! fast = direction 1, slow = direction 3
              d_rk(1) = (fk_ex_col(k_neighbor(ibz,1,    nside,i_fast, 1)) &
                       - fk_ex_col(k_neighbor(ibz,1,    nside,i_fast,-1))) / (2.0d0*slice1)
              d_rk(3) = (fk_ex_col(k_neighbor(ibz,nside,nside,i_slow, 1)) &
                       - fk_ex_col(k_neighbor(ibz,nside,nside,i_slow,-1))) / (2.0d0*slice3)
            else                                    ! fast = direction 2, slow = direction 3
              d_rk(2) = (fk_ex_col(k_neighbor(ibz,1,    nside,i_fast, 1)) &
                       - fk_ex_col(k_neighbor(ibz,1,    nside,i_fast,-1))) / (2.0d0*slice2)
              d_rk(3) = (fk_ex_col(k_neighbor(ibz,nside,nside,i_slow, 1)) &
                       - fk_ex_col(k_neighbor(ibz,nside,nside,i_slow,-1))) / (2.0d0*slice3)
            end if

          else  ! ndim == 3

            i_fast = mod(ibz-1, nside) + 1
            i_mid  = mod((ibz-1)/nside, nside) + 1
            i_slow = (ibz-1)/(nside*nside) + 1
            d_rk(1) = (fk_ex_col(k_neighbor(ibz,1,          nside,i_fast, 1)) &
                     - fk_ex_col(k_neighbor(ibz,1,          nside,i_fast,-1))) / (2.0d0*slice1)
            d_rk(2) = (fk_ex_col(k_neighbor(ibz,nside,      nside,i_mid,  1)) &
                     - fk_ex_col(k_neighbor(ibz,nside,      nside,i_mid, -1))) / (2.0d0*slice2)
            d_rk(3) = (fk_ex_col(k_neighbor(ibz,nside*nside,nside,i_slow, 1)) &
                     - fk_ex_col(k_neighbor(ibz,nside*nside,nside,i_slow,-1))) / (2.0d0*slice3)

          end if

          ! Chain rule: df/dk_a = sum_i (df/drk_i) * d(rk_i)/dk_a
          do nj = 1, 3
            fk_der_col(nj,ibz) = d_rk(1)*dRkdK(1,nj) + d_rk(2)*dRkdK(2,nj) + d_rk(3)*dRkdK(3,nj)
          end do

        end do  ! ibz

        do ibz = 1, npointstotal
          j_aux = (ibz-1)*norb_ex_band + j
          fk_ex_der(:, j_aux, nn) = fk_der_col(:, ibz)
        end do

        deallocate(fk_ex_col, fk_der_col)

      end do  ! nn
    end do  ! j
    !$omp end parallel do

  end subroutine get_fk_ex_der_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Central-difference stencil of the k-mesh, with the SAME geometry as get_fk_ex_der_k (crystal
  !! direction 1 fastest, periodic neighbours): for each point, the flat index of its +/- neighbour
  !! along each crystal direction, 1/(2*h) for that direction (0 if inactive, and then the neighbour
  !! is the point itself), and the exact crystal->Cartesian Jacobian dRkdK(i,a) = d(rk_i)/d(k_a).
  !! d f/d k_a = sum_i dRkdK(i,a) * inv2h(i) * (f(nbp(i,k)) - f(nbm(i,k))). Used by the covariant
  !! envelope derivative, which has to transport f before differencing it.
  subroutine get_kgrid_stencil(nbp, nbm, inv2h, dRkdK)
    implicit none
    integer, intent(out) :: nbp(3, npointstotal), nbm(3, npointstotal)
    real(8), intent(out) :: inv2h(3), dRkdK(3,3)
    real(8) :: rk1(npointstotal), rk2(npointstotal), rk3(npointstotal)
    real(8) :: dum1,dum2,dum3,dum4,dum5,dum6, h(3)
    integer :: ibz, nside, i_fast, i_mid, i_slow, d1, d2

    if (.not. active_flags_set) call set_active_flags()
    do ibz = 1, npointstotal
      call get_k_kc(G, rkxvector(ibz), rkyvector(ibz), rkzvector(ibz), rk1(ibz), rk2(ibz), rk3(ibz), &
                    dum1,dum2,dum3,dum4,dum5,dum6)
      nbp(:,ibz) = ibz
      nbm(:,ibz) = ibz
    end do
    call get_k_kc(G, 1.0d0,0.0d0,0.0d0, dRkdK(1,1),dRkdK(2,1),dRkdK(3,1), dum1,dum2,dum3,dum4,dum5,dum6)
    call get_k_kc(G, 0.0d0,1.0d0,0.0d0, dRkdK(1,2),dRkdK(2,2),dRkdK(3,2), dum1,dum2,dum3,dum4,dum5,dum6)
    call get_k_kc(G, 0.0d0,0.0d0,1.0d0, dRkdK(1,3),dRkdK(2,3),dRkdK(3,3), dum1,dum2,dum3,dum4,dum5,dum6)

    h = 0.0d0
    if (ndim == 1) then
      d1 = 1
      if (active_y) d1 = 2
      if (active_z) d1 = 3
      select case (d1)
        case (1); h(1) = rk1(2) - rk1(1)
        case (2); h(2) = rk2(2) - rk2(1)
        case (3); h(3) = rk3(2) - rk3(1)
      end select
      do ibz = 1, npointstotal
        nbp(d1,ibz) = k_neighbor(ibz,1,npointstotal,ibz, 1)
        nbm(d1,ibz) = k_neighbor(ibz,1,npointstotal,ibz,-1)
      end do
    else if (ndim == 2) then
      nside = nint(sqrt(dble(npointstotal)))
      if (active_x .and. active_y) then
        d1 = 1; d2 = 2; h(1) = rk1(2) - rk1(1); h(2) = rk2(nside+1) - rk2(1)
      else if (active_x .and. active_z) then
        d1 = 1; d2 = 3; h(1) = rk1(2) - rk1(1); h(3) = rk3(nside+1) - rk3(1)
      else
        d1 = 2; d2 = 3; h(2) = rk2(2) - rk2(1); h(3) = rk3(nside+1) - rk3(1)
      end if
      do ibz = 1, npointstotal
        i_fast = mod(ibz-1, nside) + 1
        i_slow = (ibz-1)/nside + 1
        nbp(d1,ibz) = k_neighbor(ibz,1,    nside,i_fast, 1)
        nbm(d1,ibz) = k_neighbor(ibz,1,    nside,i_fast,-1)
        nbp(d2,ibz) = k_neighbor(ibz,nside,nside,i_slow, 1)
        nbm(d2,ibz) = k_neighbor(ibz,nside,nside,i_slow,-1)
      end do
    else
      nside = nint(dble(npointstotal)**(1.0d0/3.0d0))
      h(1) = rk1(2) - rk1(1); h(2) = rk2(nside+1) - rk2(1); h(3) = rk3(nside*nside+1) - rk3(1)
      do ibz = 1, npointstotal
        i_fast = mod(ibz-1, nside) + 1
        i_mid  = mod((ibz-1)/nside, nside) + 1
        i_slow = (ibz-1)/(nside*nside) + 1
        nbp(1,ibz) = k_neighbor(ibz,1,          nside,i_fast, 1)
        nbm(1,ibz) = k_neighbor(ibz,1,          nside,i_fast,-1)
        nbp(2,ibz) = k_neighbor(ibz,nside,      nside,i_mid,  1)
        nbm(2,ibz) = k_neighbor(ibz,nside,      nside,i_mid, -1)
        nbp(3,ibz) = k_neighbor(ibz,nside*nside,nside,i_slow, 1)
        nbm(3,ibz) = k_neighbor(ibz,nside*nside,nside,i_slow,-1)
      end do
    end if
    inv2h = 0.0d0
    where (h /= 0.0d0) inv2h = 1.0d0/(2.0d0*h)
  end subroutine get_kgrid_stencil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Difference directions of the covariant envelope derivative of X_nm (ome_ex build_xcov): the Cartesian gradient is
  !! grad f = sum_j cfd(j,:) D_j f, with D_j f = (4/3)(T(k+s_j) - T(k-s_j)) - (1/6)(T(k+2s_j) - T(k-2s_j)) the transported
  !! difference along mesh step s_j (fourth order; the 1/2 of a central difference is in cfd).
  !!
  !! Default: the crystal axes, cfd = dRkdK * 1/(2h). On a 2D HEXAGONAL mesh (equal steps 60 or 120 degrees apart) the
  !! third star direction, s_3 = s_1 - s_2 (60 degrees) or s_1 + s_2 (120 degrees), is a mesh step too, and the gradient
  !! is taken from all three with the least-squares weights (sum_j c_j c_j^T)^-1 c_j = 2 c_j / (3 |c|^2) (c_j the Cartesian
  !! steps): the operator is then mapped onto itself by C3. The two-axis stencil is not, and for exciton envelopes that
  !! the mesh resolves poorly (excited and near-gap excitons, a few mesh points wide) it makes X_nm anisotropic: MoSe2
  !! with spin-orbit coupling, 87x87, lowest 2078 excitons, excitonic SHG C3 residual 38.6% with the axes, 1.7% with the
  !! star; a 120x120 mesh does not cure it (30.6%), because the near-gap states only get narrower. In the non-interacting
  !! limit the star is as accurate or slightly more (hBN 30x30, excitonic vs single-particle SHG 2.6e-3 -> 2.1e-3).
  !! Mesh steps are compared to 1e-4 (the lattice of a _tb.dat file has about 6 significant figures).
  subroutine get_xnm_directions(ndir, dirp, dirm, cfd, hex)
    implicit none
    integer, intent(out) :: ndir
    integer, intent(out) :: dirp(3, npointstotal), dirm(3, npointstotal)
    real(8), intent(out) :: cfd(3,3)
    logical, intent(out) :: hex
    integer, allocatable :: nbp(:,:), nbm(:,:)
    real(8) :: inv2h(3), dRkdK(3,3), B(2,2), Bi(2,2), det, c(3,3), l1, l2, lm, lp
    integer :: idim, d(3), ca(2), nca, j, k

    allocate(nbp(3, npointstotal), nbm(3, npointstotal))
    call get_kgrid_stencil(nbp, nbm, inv2h, dRkdK)
    ndir = 0; cfd = 0.0d0; hex = .false.
    do idim = 1, 3
      if (inv2h(idim) == 0.0d0) cycle
      ndir = ndir + 1; d(ndir) = idim
      dirp(ndir,:) = nbp(idim,:); dirm(ndir,:) = nbm(idim,:)
      cfd(ndir,:) = dRkdK(idim,:)*inv2h(idim)
    end do
    if (ndim /= 2 .or. ndir /= 2) return

    ! Cartesian steps of the two mesh axes: in-plane block B (crystal d(i) vs periodic Cartesian axis ca(j)), c = h B^-1 e_d
    nca = 0
    if (active_x) then; nca = nca + 1; ca(nca) = 1; end if
    if (active_y) then; nca = nca + 1; ca(nca) = 2; end if
    if (active_z) then; nca = nca + 1; ca(nca) = 3; end if
    if (nca /= 2) return
    do j = 1, 2
      do k = 1, 2
        B(j,k) = dRkdK(d(j), ca(k))
      end do
    end do
    det = B(1,1)*B(2,2) - B(1,2)*B(2,1)
    if (det == 0.0d0) return
    Bi = reshape([B(2,2), -B(2,1), -B(1,2), B(1,1)], [2,2])/det
    c = 0.0d0
    do j = 1, 2
      c(ca, j) = Bi(:, j)/(2.0d0*inv2h(d(j)))
    end do
    l1 = norm2(c(:,1)); l2 = norm2(c(:,2)); lm = norm2(c(:,1) - c(:,2)); lp = norm2(c(:,1) + c(:,2))
    if (abs(l1 - l2) > 1.0d-4*l1) return
    if (abs(lm - l1) < 1.0d-4*l1) then
      c(:,3) = c(:,1) - c(:,2)
      do k = 1, npointstotal
        dirp(3,k) = nbp(d(1), nbm(d(2), k)); dirm(3,k) = nbm(d(1), nbp(d(2), k))
      end do
    else if (abs(lp - l1) < 1.0d-4*l1) then
      c(:,3) = c(:,1) + c(:,2)
      do k = 1, npointstotal
        dirp(3,k) = nbp(d(1), nbp(d(2), k)); dirm(3,k) = nbm(d(1), nbm(d(2), k))
      end do
    else
      return
    end if
    ndir = 3; hex = .true.
    do j = 1, 3
      cfd(j,:) = 0.5d0*c(:,j)*2.0d0/(3.0d0*l1*l1)
    end do
    deallocate(nbp, nbm)
  end subroutine get_xnm_directions

  !> .true. when get_xnm_directions uses the hexagonal star (recorded in the excitonic cache).
  logical function xnm_hex_mesh() result(hex)
    implicit none
    integer, allocatable :: dirp(:,:), dirm(:,:)
    integer :: ndir
    real(8) :: cfd(3,3)
    allocate(dirp(3, npointstotal), dirm(3, npointstotal))
    call get_xnm_directions(ndir, dirp, dirm, cfd, hex)
    deallocate(dirp, dirm)
  end function xnm_hex_mesh


subroutine interp1(npointstotal, rk1v, rk2v, rk3v, fk_col, rk1, rk2, rk3, result)
  implicit none
  integer,    intent(in)  :: npointstotal
  real(8),    intent(in)  :: rk1v(npointstotal), rk2v(npointstotal), rk3v(npointstotal)
  complex(8), intent(in)  :: fk_col(npointstotal)
  real(8),    intent(in)  :: rk1, rk2, rk3
  complex(8), intent(out) :: result
  call get_fk_ex_k_interp_1d(npointstotal, rk1v, rk2v, rk3v, fk_col, rk1, rk2, rk3, result)
end subroutine interp1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rank-1 counterpart of get_fk_ex_k_interp, operating on a single
! (npointstotal)-length column instead of the full (npointstotal,norb_ex_cut)
! array. Avoids the per-call (npointstotal,1) copy that interp1 previously
! had to build to satisfy get_fk_ex_k_interp's rank-2 explicit-shape dummy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_fk_ex_k_interp_1d(npointstotal, rk1vector, rk2vector, rk3vector, &
                                   fk_ex_1d, rk1, rk2, rk3, fk_ex_k_interp)
  implicit none

  integer,    intent(in)  :: npointstotal
  real(8),    intent(in)  :: rk1vector(npointstotal), rk2vector(npointstotal), rk3vector(npointstotal)
  complex(8), intent(in)  :: fk_ex_1d(npointstotal)
  real(8),    intent(in)  :: rk1, rk2, rk3
  complex(8), intent(out) :: fk_ex_k_interp

  integer :: nside, nblock
  integer :: ibz_q11,  ibz_q21,  ibz_q31,  ibz_q12,  ibz_q22,  ibz_q32,  ibz_q13,  ibz_q23,  ibz_q33
  integer :: ibz_q111, ibz_q211, ibz_q121, ibz_q221, ibz_q112, ibz_q212, ibz_q122, ibz_q222
  real(8) :: slice
  real(8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3, x, y, z

  if (.not. active_flags_set) call set_active_flags()

  if (ndim == 1) then
    nside = nint(dble(npointstotal))
    slice = 1.0d0 / dble(nside-1)
    if (active_z) then
      nblock = int((rk3+0.5d0)/slice) + 1
      ibz_q13 = nblock;  ibz_q33 = nblock+1
      if (ibz_q13 < 1 .or. ibz_q33 > npointstotal) then
        write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
      end if
      z1 = rk3vector(ibz_q13);  z3 = rk3vector(ibz_q33);  z = rk3
      fk_ex_k_interp = fk_ex_1d(ibz_q13)*(z3-z)/(z3-z1) + fk_ex_1d(ibz_q33)*(z-z1)/(z3-z1)
    elseif (active_y) then
      nblock = int((rk2+0.5d0)/slice) + 1
      ibz_q12 = nblock;  ibz_q22 = nblock+1
      if (ibz_q12 < 1 .or. ibz_q22 > npointstotal) then
        write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
      end if
      y1 = rk2vector(ibz_q12);  y2 = rk2vector(ibz_q22);  y = rk2
      fk_ex_k_interp = fk_ex_1d(ibz_q12)*(y2-y)/(y2-y1) + fk_ex_1d(ibz_q22)*(y-y1)/(y2-y1)
    else
      nblock = int((rk1+0.5d0)/slice) + 1
      ibz_q11 = nblock;  ibz_q21 = nblock+1
      if (ibz_q11 < 1 .or. ibz_q21 > npointstotal) then
        write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
      end if
      x1 = rk1vector(ibz_q11);  x2 = rk1vector(ibz_q21);  x = rk1
      fk_ex_k_interp = fk_ex_1d(ibz_q11)*(x2-x)/(x2-x1) + fk_ex_1d(ibz_q21)*(x-x1)/(x2-x1)
    end if

  elseif (ndim == 2) then
    nside = nint(sqrt(dble(npointstotal)))
    slice = 1.0d0 / dble(nside-1)
    if (active_x .and. active_y) then
      nblock  = int((rk1+0.5d0)/slice) + 1 + int((rk2+0.5d0)/slice)*nside
      ibz_q11 = nblock;  ibz_q21 = nblock+1
      ibz_q12 = nblock+nside;  ibz_q22 = nblock+nside+1
      if (ibz_q11 < 1 .or. ibz_q22 > npointstotal) then
        write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
      end if
      x1=rk1vector(ibz_q11); x2=rk1vector(ibz_q21)
      y1=rk2vector(ibz_q11); y2=rk2vector(ibz_q12)
      x=rk1; y=rk2
      fk_ex_k_interp = 1.0d0/((x2-x1)*(y2-y1)) * &
        ( fk_ex_1d(ibz_q11)*(x2-x)*(y2-y) + fk_ex_1d(ibz_q21)*(x-x1)*(y2-y) &
        + fk_ex_1d(ibz_q12)*(x2-x)*(y-y1) + fk_ex_1d(ibz_q22)*(x-x1)*(y-y1) )
    elseif (active_x .and. active_z) then
      nblock  = int((rk1+0.5d0)/slice) + 1 + int((rk3+0.5d0)/slice)*nside
      ibz_q11 = nblock;  ibz_q31 = nblock+1
      ibz_q13 = nblock+nside;  ibz_q33 = nblock+nside+1
      if (ibz_q11 < 1 .or. ibz_q33 > npointstotal) then
        write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
      end if
      x1=rk1vector(ibz_q11); x3=rk1vector(ibz_q31)
      z1=rk3vector(ibz_q11); z3=rk3vector(ibz_q13)
      x=rk1; z=rk3
      fk_ex_k_interp = 1.0d0/((x3-x1)*(z3-z1)) * &
        ( fk_ex_1d(ibz_q11)*(x3-x)*(z3-z) + fk_ex_1d(ibz_q31)*(x-x1)*(z3-z) &
        + fk_ex_1d(ibz_q13)*(x3-x)*(z-z1) + fk_ex_1d(ibz_q33)*(x-x1)*(z-z1) )
    else
      nblock  = int((rk2+0.5d0)/slice) + 1 + int((rk3+0.5d0)/slice)*nside
      ibz_q22 = nblock;  ibz_q32 = nblock+1
      ibz_q23 = nblock+nside;  ibz_q33 = nblock+nside+1
      if (ibz_q22 < 1 .or. ibz_q33 > npointstotal) then
        write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
      end if
      y2=rk2vector(ibz_q22); y3=rk2vector(ibz_q32)
      z2=rk3vector(ibz_q22); z3=rk3vector(ibz_q23)
      y=rk2; z=rk3
      fk_ex_k_interp = 1.0d0/((y3-y2)*(z3-z2)) * &
        ( fk_ex_1d(ibz_q22)*(y3-y)*(z3-z) + fk_ex_1d(ibz_q32)*(y-y2)*(z3-z) &
        + fk_ex_1d(ibz_q23)*(y3-y)*(z-z2) + fk_ex_1d(ibz_q33)*(y-y2)*(z-z2) )
    end if

  else  ! ndim == 3
    nside = nint(dble(npointstotal)**(1.0d0/3.0d0))
    slice = 1.0d0 / dble(nside-1)
    nblock = int((rk1+0.5d0)/slice) + int((rk2+0.5d0)/slice)*nside &
           + int((rk3+0.5d0)/slice)*nside*nside + 1
    ibz_q111=nblock;           ibz_q211=nblock+1
    ibz_q121=nblock+nside;     ibz_q221=nblock+nside+1
    ibz_q112=nblock+nside*nside;       ibz_q212=nblock+nside*nside+1
    ibz_q122=nblock+nside*nside+nside; ibz_q222=nblock+nside*nside+nside+1
    if (ibz_q111 < 1 .or. ibz_q222 > npointstotal) then
      write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
    end if
    x1=rk1vector(ibz_q111); x2=rk1vector(ibz_q211)
    y1=rk2vector(ibz_q111); y2=rk2vector(ibz_q121)
    z1=rk3vector(ibz_q111); z2=rk3vector(ibz_q112)
    x=rk1; y=rk2; z=rk3
    fk_ex_k_interp = 1.0d0/((x2-x1)*(y2-y1)*(z2-z1)) * &
      ( fk_ex_1d(ibz_q111)*(x2-x)*(y2-y)*(z2-z) &
      + fk_ex_1d(ibz_q211)*(x-x1)*(y2-y)*(z2-z) &
      + fk_ex_1d(ibz_q121)*(x2-x)*(y-y1)*(z2-z) &
      + fk_ex_1d(ibz_q221)*(x-x1)*(y-y1)*(z2-z) &
      + fk_ex_1d(ibz_q112)*(x2-x)*(y2-y)*(z-z1) &
      + fk_ex_1d(ibz_q212)*(x-x1)*(y2-y)*(z-z1) &
      + fk_ex_1d(ibz_q122)*(x2-x)*(y-y1)*(z-z1) &
      + fk_ex_1d(ibz_q222)*(x-x1)*(y-y1)*(z-z1) )
  end if

end subroutine get_fk_ex_k_interp_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_fk_ex_k_interp(npointstotal, rk1vector, rk2vector, rk3vector, &
                                  norb_ex, norb_ex_cut, fk_ex, rk1, rk2, rk3, &
                                  nn, fk_ex_k_interp)
    implicit none

    integer,    intent(in)  :: npointstotal, norb_ex, norb_ex_cut, nn
    real(8),    intent(in)  :: rk1vector(npointstotal), rk2vector(npointstotal), rk3vector(npointstotal)
    complex(8), intent(in)  :: fk_ex(npointstotal, norb_ex_cut)
    real(8),    intent(in)  :: rk1, rk2, rk3
    complex(8), intent(out) :: fk_ex_k_interp

    integer :: nside, nblock
    integer :: ibz_q11,  ibz_q21,  ibz_q31,  ibz_q12,  ibz_q22,  ibz_q32,  ibz_q13,  ibz_q23,  ibz_q33
    integer :: ibz_q111, ibz_q211, ibz_q121, ibz_q221, ibz_q112, ibz_q212, ibz_q122, ibz_q222
    real(8) :: slice
    real(8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3, x, y, z

    ! PATCH: no longer recomputed here -- module-level flags.
    if (.not. active_flags_set) call set_active_flags()

    if (ndim == 1) then
      nside = nint(dble(npointstotal))
      slice = 1.0d0 / dble(nside-1)
      if (active_z) then
        nblock = int((rk3+0.5d0)/slice) + 1
        ibz_q13 = nblock;  ibz_q33 = nblock+1
        if (ibz_q13 < 1 .or. ibz_q33 > npointstotal) then
          write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
        end if
        z1 = rk3vector(ibz_q13);  z3 = rk3vector(ibz_q33);  z = rk3
        fk_ex_k_interp = fk_ex(ibz_q13,nn)*(z3-z)/(z3-z1) + fk_ex(ibz_q33,nn)*(z-z1)/(z3-z1)
      elseif (active_y) then
        nblock = int((rk2+0.5d0)/slice) + 1
        ibz_q12 = nblock;  ibz_q22 = nblock+1
        if (ibz_q12 < 1 .or. ibz_q22 > npointstotal) then
          write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
        end if
        y1 = rk2vector(ibz_q12);  y2 = rk2vector(ibz_q22);  y = rk2
        fk_ex_k_interp = fk_ex(ibz_q12,nn)*(y2-y)/(y2-y1) + fk_ex(ibz_q22,nn)*(y-y1)/(y2-y1)
      else
        nblock = int((rk1+0.5d0)/slice) + 1
        ibz_q11 = nblock;  ibz_q21 = nblock+1
        if (ibz_q11 < 1 .or. ibz_q21 > npointstotal) then
          write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
        end if
        x1 = rk1vector(ibz_q11);  x2 = rk1vector(ibz_q21);  x = rk1
        fk_ex_k_interp = fk_ex(ibz_q11,nn)*(x2-x)/(x2-x1) + fk_ex(ibz_q21,nn)*(x-x1)/(x2-x1)
      end if

    elseif (ndim == 2) then
      nside = nint(sqrt(dble(npointstotal)))
      slice = 1.0d0 / dble(nside-1)
      if (active_x .and. active_y) then
        nblock  = int((rk1+0.5d0)/slice) + 1 + int((rk2+0.5d0)/slice)*nside
        ibz_q11 = nblock;  ibz_q21 = nblock+1
        ibz_q12 = nblock+nside;  ibz_q22 = nblock+nside+1
        if (ibz_q11 < 1 .or. ibz_q22 > npointstotal) then
          write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
        end if
        x1=rk1vector(ibz_q11); x2=rk1vector(ibz_q21)
        y1=rk2vector(ibz_q11); y2=rk2vector(ibz_q12)
        x=rk1; y=rk2
        fk_ex_k_interp = 1.0d0/((x2-x1)*(y2-y1)) * &
          ( fk_ex(ibz_q11,nn)*(x2-x)*(y2-y) + fk_ex(ibz_q21,nn)*(x-x1)*(y2-y) &
          + fk_ex(ibz_q12,nn)*(x2-x)*(y-y1) + fk_ex(ibz_q22,nn)*(x-x1)*(y-y1) )
      elseif (active_x .and. active_z) then
        nblock  = int((rk1+0.5d0)/slice) + 1 + int((rk3+0.5d0)/slice)*nside
        ibz_q11 = nblock;  ibz_q31 = nblock+1
        ibz_q13 = nblock+nside;  ibz_q33 = nblock+nside+1
        if (ibz_q11 < 1 .or. ibz_q33 > npointstotal) then
          write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
        end if
        x1=rk1vector(ibz_q11); x3=rk1vector(ibz_q31)
        z1=rk3vector(ibz_q11); z3=rk3vector(ibz_q13)
        x=rk1; z=rk3
        fk_ex_k_interp = 1.0d0/((x3-x1)*(z3-z1)) * &
          ( fk_ex(ibz_q11,nn)*(x3-x)*(z3-z) + fk_ex(ibz_q31,nn)*(x-x1)*(z3-z) &
          + fk_ex(ibz_q13,nn)*(x3-x)*(z-z1) + fk_ex(ibz_q33,nn)*(x-x1)*(z-z1) )
      else
        nblock  = int((rk2+0.5d0)/slice) + 1 + int((rk3+0.5d0)/slice)*nside
        ibz_q22 = nblock;  ibz_q32 = nblock+1
        ibz_q23 = nblock+nside;  ibz_q33 = nblock+nside+1
        if (ibz_q22 < 1 .or. ibz_q33 > npointstotal) then
          write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
        end if
        y2=rk2vector(ibz_q22); y3=rk2vector(ibz_q32)
        z2=rk3vector(ibz_q22); z3=rk3vector(ibz_q23)
        y=rk2; z=rk3
        fk_ex_k_interp = 1.0d0/((y3-y2)*(z3-z2)) * &
          ( fk_ex(ibz_q22,nn)*(y3-y)*(z3-z) + fk_ex(ibz_q32,nn)*(y-y2)*(z3-z) &
          + fk_ex(ibz_q23,nn)*(y3-y)*(z-z2) + fk_ex(ibz_q33,nn)*(y-y2)*(z-z2) )
      end if

    else  ! ndim == 3
      nside = nint(dble(npointstotal)**(1.0d0/3.0d0))
      slice = 1.0d0 / dble(nside-1)
      nblock = int((rk1+0.5d0)/slice) + int((rk2+0.5d0)/slice)*nside &
             + int((rk3+0.5d0)/slice)*nside*nside + 1
      ibz_q111=nblock;           ibz_q211=nblock+1
      ibz_q121=nblock+nside;     ibz_q221=nblock+nside+1
      ibz_q112=nblock+nside*nside;       ibz_q212=nblock+nside*nside+1
      ibz_q122=nblock+nside*nside+nside; ibz_q222=nblock+nside*nside+nside+1
      if (ibz_q111 < 1 .or. ibz_q222 > npointstotal) then
        write(*,*) 'interpolation went wrong'; fk_ex_k_interp = (0.0d0,0.0d0); return
      end if
      x1=rk1vector(ibz_q111); x2=rk1vector(ibz_q211)
      y1=rk2vector(ibz_q111); y2=rk2vector(ibz_q121)
      z1=rk3vector(ibz_q111); z2=rk3vector(ibz_q112)
      x=rk1; y=rk2; z=rk3
      fk_ex_k_interp = 1.0d0/((x2-x1)*(y2-y1)*(z2-z1)) * &
        ( fk_ex(ibz_q111,nn)*(x2-x)*(y2-y)*(z2-z) &
        + fk_ex(ibz_q211,nn)*(x-x1)*(y2-y)*(z2-z) &
        + fk_ex(ibz_q121,nn)*(x2-x)*(y-y1)*(z2-z) &
        + fk_ex(ibz_q221,nn)*(x-x1)*(y-y1)*(z2-z) &
        + fk_ex(ibz_q112,nn)*(x2-x)*(y2-y)*(z-z1) &
        + fk_ex(ibz_q212,nn)*(x-x1)*(y2-y)*(z-z1) &
        + fk_ex(ibz_q122,nn)*(x2-x)*(y-y1)*(z-z1) &
        + fk_ex(ibz_q222,nn)*(x-x1)*(y-y1)*(z-z1) )
    end if

  end subroutine get_fk_ex_k_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_k_kc(G, xp, yp, zp, xc, yc, zc, xp_bz, yp_bz, zp_bz, xc_bz, yc_bz, zc_bz)
    implicit none
    real(8), intent(in)  :: G(3,3)
    real(8), intent(in)  :: xp, yp, zp
    real(8), intent(out) :: xc, yc, zc, xp_bz, yp_bz, zp_bz, xc_bz, yc_bz, zc_bz
    real(8) :: det

    ! PATCH: no longer recomputed here -- module-level flags.
    if (.not. active_flags_set) call set_active_flags()

    if (ndim == 1) then
      if (active_z) then
        xc=0.0d0; yc=0.0d0; zc=zp/G(3,3)
        xc_bz=0.0d0; yc_bz=0.0d0; zc_bz=zc-dble(int(zc/0.5d0))
        xp_bz=0.0d0; yp_bz=0.0d0; zp_bz=zc_bz*G(3,3)
      elseif (active_y) then
        xc=0.0d0; yc=yp/G(2,2); zc=0.0d0
        xc_bz=0.0d0; yc_bz=yc-dble(int(yc/0.5d0)); zc_bz=0.0d0
        xp_bz=0.0d0; yp_bz=yc_bz*G(2,2); zp_bz=0.0d0
      else
        xc=xp/G(1,1); yc=0.0d0; zc=0.0d0
        xc_bz=xc-dble(int(xc/0.5d0)); yc_bz=0.0d0; zc_bz=0.0d0
        xp_bz=xc_bz*G(1,1); yp_bz=0.0d0; zp_bz=0.0d0
      end if
    elseif (ndim == 2) then
      if (active_x .and. active_y) then
        xc=(xp*G(2,2)-G(2,1)*yp)/(G(1,1)*G(2,2)-G(2,1)*G(1,2))
        yc=(G(1,1)*yp-xp*G(1,2))/(G(1,1)*G(2,2)-G(2,1)*G(1,2))
        zc=0.0d0
        xc_bz=xc-dble(int(xc/0.5d0)); yc_bz=yc-dble(int(yc/0.5d0)); zc_bz=0.0d0
        xp_bz=xc_bz*G(1,1)+yc_bz*G(2,1); yp_bz=xc_bz*G(1,2)+yc_bz*G(2,2); zp_bz=0.0d0
      elseif (active_x .and. active_z) then
        xc=(xp*G(3,3)-G(3,1)*zp)/(G(1,1)*G(3,3)-G(3,1)*G(1,3))
        yc=0.0d0
        zc=(G(1,1)*zp-xp*G(1,3))/(G(1,1)*G(3,3)-G(3,1)*G(1,3))
        xc_bz=xc-dble(int(xc/0.5d0)); yc_bz=0.0d0; zc_bz=zc-dble(int(zc/0.5d0))
        xp_bz=xc_bz*G(1,1)+zc_bz*G(3,1); yp_bz=0.0d0; zp_bz=xc_bz*G(1,3)+zc_bz*G(3,3)
      else
        xc=0.0d0
        yc=(yp*G(3,3)-G(3,2)*zp)/(G(2,2)*G(3,3)-G(3,2)*G(2,3))
        zc=(G(2,2)*zp-yp*G(2,3))/(G(2,2)*G(3,3)-G(3,2)*G(2,3))
        xc_bz=0.0d0; yc_bz=yc-dble(int(yc/0.5d0)); zc_bz=zc-dble(int(zc/0.5d0))
        xp_bz=0.0d0; yp_bz=yc_bz*G(2,2)+zc_bz*G(3,2); zp_bz=yc_bz*G(2,3)+zc_bz*G(3,3)
      end if
    else
      det = -G(1,3)*G(2,2)*G(3,1) + G(1,2)*G(2,3)*G(3,1) + G(1,3)*G(2,1)*G(3,2) &
            -G(1,1)*G(2,3)*G(3,2) - G(1,2)*G(2,1)*G(3,3) + G(1,1)*G(2,2)*G(3,3)
!       xc=(xp*(G(2,2)*G(3,3)-G(2,3)*G(3,2))+yp*(G(1,3)*G(3,2)-G(1,2)*G(3,3))+zp*(G(1,2)*G(2,3)-G(1,3)*G(2,2)))/det
!       yc=(xp*(G(2,3)*G(3,1)-G(2,1)*G(3,3))+yp*(G(1,1)*G(3,3)-G(1,3)*G(3,1))+zp*(G(1,3)*G(2,1)-G(1,1)*G(2,3)))/det
!       zc=(xp*(G(2,1)*G(3,2)-G(2,2)*G(3,1))+yp*(G(1,2)*G(3,1)-G(1,1)*G(3,2))+zp*(G(1,1)*G(2,2)-G(1,2)*G(2,1)))/det
      xc=(xp*(G(2,2)*G(3,3)-G(2,3)*G(3,2))+yp*(G(2,3)*G(3,1)-G(2,1)*G(3,3))+zp*(G(2,1)*G(3,2)-G(2,2)*G(3,1)))/det
      yc=(xp*(G(1,3)*G(3,2)-G(1,2)*G(3,3))+yp*(G(1,1)*G(3,3)-G(1,3)*G(3,1))+zp*(G(1,2)*G(3,1)-G(1,1)*G(3,2)))/det
      zc=(xp*(G(1,2)*G(2,3)-G(1,3)*G(2,2))+yp*(G(1,3)*G(2,1)-G(1,1)*G(2,3))+zp*(G(1,1)*G(2,2)-G(1,2)*G(2,1)))/det

      xc_bz=xc-dble(int(xc/0.5d0)); yc_bz=yc-dble(int(yc/0.5d0)); zc_bz=zc-dble(int(zc/0.5d0))
      xp_bz=xc_bz*G(1,1)+yc_bz*G(2,1)+zc_bz*G(3,1)
      yp_bz=xc_bz*G(1,2)+yc_bz*G(2,2)+zc_bz*G(3,2)
      zp_bz=xc_bz*G(1,3)+yc_bz*G(2,3)+zc_bz*G(3,3)
    end if

  end subroutine get_k_kc

end module exciton_envelopes
