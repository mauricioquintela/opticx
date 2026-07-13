module ome_ex
  use constants_math
  use parser_input_file, &
    only: nf, e1, e2, eta, nw
  use parser_wannier90_tb, &
    only: material_name, norb
  use parser_optics_xatu_dim, &
    only: npointstotal, vcell, &
          norb_ex, norb_ex_cut, nv_ex, nc_ex, nband_ex, e_ex, fk_ex, &
          get_ex_index_first, print_exciton_wf, &
          rkxvector, rkyvector, rkzvector
  use exciton_envelopes, &
    only: fk_ex_der, get_fk_ex_der_k
  implicit none

  ! ex ground-state matrix elements
  complex(8), allocatable :: xme_ex(:,:)        ! (3, norb_ex_cut)
  complex(8), allocatable :: vme_ex(:,:)        ! (3, norb_ex_cut)

  ! ex inter-state matrix elements
  complex(8), allocatable :: qme_ex_inter1(:,:,:)
  complex(8), allocatable :: qme_ex_inter2(:,:,:)
  complex(8), allocatable :: qme_ex_inter(:,:,:)
  complex(8), allocatable :: yme_ex_inter1(:,:,:)
  complex(8), allocatable :: yme_ex_inter2(:,:,:)
  complex(8), allocatable :: yme_ex_inter(:,:,:)
  complex(8), allocatable :: xme_ex_inter(:,:,:)
  complex(8), allocatable :: vme_ex_inter1(:,:,:)
  complex(8), allocatable :: vme_ex_inter2(:,:,:)
  complex(8), allocatable :: vme_ex_inter(:,:,:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_ex(iflag_norder)
    implicit none
    integer, intent(in) :: iflag_norder

    integer :: ibz, nn, nnp, nj
    integer :: u_exk
    logical :: do_write_exk

    complex(8), allocatable :: vme_ex_k(:,:)

    real(8),    allocatable :: ek(:,:)
    complex(8), allocatable :: xme_ex_band(:,:,:,:)
    complex(8), allocatable :: vme_ex_band(:,:,:,:)
    complex(8), allocatable :: berry_eigen_ex_band(:,:,:,:)
    complex(8), allocatable :: gen_der_ex_band(:,:,:,:,:)
    real(8),    allocatable :: shift_vector_ex_band(:,:,:,:,:)

    write(*,*) '6. Entering ome_ex'

    do_write_exk = (iflag_norder == 1)

    ! --- allocate band-space arrays ---
    allocate(ek(npointstotal, nband_ex))
    allocate(vme_ex_band(npointstotal, 3, nband_ex, nband_ex))

    if (iflag_norder == 1) then
      ek = 0.0d0
      vme_ex_band = (0.0d0, 0.0d0)
      write(*,*) '   Reading optical matrix elements (sp)...'
      call read_ome_sp_linear(iflag_norder, npointstotal, nband_ex, vme_ex_band, ek)
    end if

    if (iflag_norder == 2) then
      allocate(berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex))
      allocate(gen_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex))
      allocate(shift_vector_ex_band(npointstotal, 3, 3, nband_ex, nband_ex))
      allocate(xme_ex_band(npointstotal, 3, nband_ex, nband_ex))
      ek = 0.0d0
      vme_ex_band = (0.0d0, 0.0d0)
      berry_eigen_ex_band = (0.0d0, 0.0d0)
      gen_der_ex_band = (0.0d0, 0.0d0)
      shift_vector_ex_band = 0.0d0
      write(*,*) '   Reading optical matrix elements (sp)...'
      call read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, &
                                  berry_eigen_ex_band, gen_der_ex_band, &
                                  shift_vector_ex_band, vme_ex_band, ek)
      xme_ex_band = (0.0d0, 0.0d0)
      call get_ome_sp_xme_ex_band(ek, vme_ex_band, xme_ex_band)
    end if

    ! --- allocate ex-ome arrays ---
    allocate(vme_ex(3, norb_ex_cut))
    allocate(xme_ex(3, norb_ex_cut))
    vme_ex = (0.0d0, 0.0d0)
    xme_ex = (0.0d0, 0.0d0)

    if (do_write_exk) then
      allocate(vme_ex_k(3, norb_ex_cut))
      vme_ex_k = (0.0d0, 0.0d0)
      u_exk = 77
      call write_ome_ex_linear_kresolved_init(u_exk, material_name, norb_ex_cut)
    end if

    if (iflag_norder == 2) then
      allocate(qme_ex_inter1(3, norb_ex_cut, norb_ex_cut))
      allocate(qme_ex_inter2(3, norb_ex_cut, norb_ex_cut))
      allocate(qme_ex_inter (3, norb_ex_cut, norb_ex_cut))
      allocate(yme_ex_inter1(3, norb_ex_cut, norb_ex_cut))
      allocate(yme_ex_inter2(3, norb_ex_cut, norb_ex_cut))
      allocate(yme_ex_inter (3, norb_ex_cut, norb_ex_cut))
      allocate(vme_ex_inter1(3, norb_ex_cut, norb_ex_cut))
      allocate(vme_ex_inter2(3, norb_ex_cut, norb_ex_cut))
      allocate(vme_ex_inter (3, norb_ex_cut, norb_ex_cut))
      qme_ex_inter1 = (0.0d0, 0.0d0)
      qme_ex_inter2 = (0.0d0, 0.0d0)
      qme_ex_inter  = (0.0d0, 0.0d0)
      yme_ex_inter1 = (0.0d0, 0.0d0)
      yme_ex_inter2 = (0.0d0, 0.0d0)
      yme_ex_inter  = (0.0d0, 0.0d0)
      vme_ex_inter1 = (0.0d0, 0.0d0)
      vme_ex_inter2 = (0.0d0, 0.0d0)
      vme_ex_inter  = (0.0d0, 0.0d0)

      allocate(fk_ex_der(3, norb_ex, norb_ex_cut))
      call get_fk_ex_der_k()
    end if

    if (iflag_norder == 1) write(*,*) '   Evaluating excitonic OMEs for linear conductivity...'
    if (iflag_norder == 2) write(*,*) '   Evaluating excitonic OMEs for nonlinear conductivity...'

    ! --- k-space integration ---
    do ibz = 1, npointstotal
      write(*,*) '   OME (ex): k-point', ibz, '/', npointstotal

      if (iflag_norder == 1 .or. iflag_norder == 2) &
        call get_ome_gs_ex_sum_k(ibz, ek, xme_ex_band, vme_ex_band)

      if (do_write_exk) then
        call get_ome_gs_ex_kresolved(ibz, ek, vme_ex_band, vme_ex_k)
        call write_ome_ex_linear_kresolved_point(u_exk, rkxvector(ibz), rkyvector(ibz), &
                                                  rkzvector(ibz), norb_ex_cut, vme_ex_k)
      end if

      if (iflag_norder == 2) &
        call get_ome_inter_ex_sum_k(ibz, ek, xme_ex_band, vme_ex_band, berry_eigen_ex_band)
    end do

    if (do_write_exk) then
      call write_ome_ex_linear_kresolved_close(u_exk)
      deallocate(vme_ex_k)
      write(*,*) '   k-resolved linear excitonic OMEs written (omeexk)'
    end if

    ! --- combine inter-state matrix elements ---
    if (iflag_norder == 2) then
      allocate(xme_ex_inter(3, norb_ex_cut, norb_ex_cut))
      ! Pure array arithmetic — compiler will vectorize/parallelise automatically
      qme_ex_inter  = qme_ex_inter1 + qme_ex_inter2
      yme_ex_inter  = yme_ex_inter1 + yme_ex_inter2
      xme_ex_inter  = yme_ex_inter  + qme_ex_inter
      vme_ex_inter  = vme_ex_inter1 + vme_ex_inter2
    end if

    write(*,*) '   Optical matrix elements (ex) have been evaluated'

    if (iflag_norder == 1) call write_ome_ex_linear(vme_ex)
    write(*,*) '   Optical matrix elements (ex, N->GS) written'
    write(*,*) '   Optical matrix elements (ex, N->N'') not printed in this version'

    ! --- cleanup band arrays ---
    deallocate(ek, vme_ex_band)
    if (iflag_norder == 2) then
      deallocate(berry_eigen_ex_band, gen_der_ex_band, shift_vector_ex_band, xme_ex_band)
      deallocate(fk_ex_der)
    end if

  end subroutine get_ome_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_sp_xme_ex_band(ek, vme_ex_band, xme_ex_band)
    implicit none
    real(8),    intent(in)  :: ek(npointstotal, nband_ex)
    complex(8), intent(in)  :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)

    integer  :: ibz, i, j, nj
    real(8)  :: de
    complex(8), parameter :: ci = (0.0d0, 1.0d0)

    xme_ex_band = (0.0d0, 0.0d0)

    !$omp parallel do collapse(2) schedule(static) &
    !$omp private(ibz, nj, i, j, de)
    do ibz = 1, npointstotal
      do nj = 1, 3
        do i = 1, nband_ex
          do j = 1, nband_ex
            de = ek(ibz,i) - ek(ibz,j)
            if (abs(de) > 1.0d-6) then
              xme_ex_band(ibz,nj,i,j) = -ci / de * vme_ex_band(ibz,nj,i,j)
              if (abs(xme_ex_band(ibz,nj,i,j)) > 20.0d0) &
                xme_ex_band(ibz,nj,i,j) = (0.0d0, 0.0d0)
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine get_ome_sp_xme_ex_band

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_gs_ex_sum_k(ibz, ek, xme_ex_band, vme_ex_band)
    implicit none
    integer,    intent(in) :: ibz
    real(8),    intent(in) :: ek(npointstotal, nband_ex)
    complex(8), intent(in) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in) :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)

    integer :: nn, ic, iv, i_ex_nn, nj

    !$omp parallel do schedule(static) private(nn, ic, iv, i_ex_nn, nj)
    do nn = 1, norb_ex_cut
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          call get_ex_index_first(nf, nv_ex, nc_ex, 0, ibz, i_ex_nn, ic, iv)
          do nj = 1, 3
            vme_ex(nj,nn) = vme_ex(nj,nn) + fk_ex(i_ex_nn,nn) * vme_ex_band(ibz,nj,iv,nv_ex+ic)
            xme_ex(nj,nn) = xme_ex(nj,nn) + fk_ex(i_ex_nn,nn) * xme_ex_band(ibz,nj,iv,nv_ex+ic)
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine get_ome_gs_ex_sum_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_gs_ex_kresolved(ibz, ek, vme_ex_band, vme_ex_k)
    implicit none
    integer,    intent(in)  :: ibz
    real(8),    intent(in)  :: ek(npointstotal, nband_ex)
    complex(8), intent(in)  :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: vme_ex_k(3, norb_ex_cut)

    integer :: nn, ic, iv, i_ex_nn, nj

    vme_ex_k = (0.0d0, 0.0d0)

    !$omp parallel do schedule(static) private(nn, ic, iv, i_ex_nn, nj)
    do nn = 1, norb_ex_cut
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          call get_ex_index_first(nf, nv_ex, nc_ex, 0, ibz, i_ex_nn, ic, iv)
          do nj = 1, 3
            vme_ex_k(nj,nn) = vme_ex_k(nj,nn) + fk_ex(i_ex_nn,nn) * vme_ex_band(ibz,nj,iv,nv_ex+ic)
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine get_ome_gs_ex_kresolved

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_inter_ex_sum_k(ibz, ek, xme_ex_band, vme_ex_band, berry_eigen_ex_band)
    !--------------------------------------------------------------------------
    ! Restructured to replace the O(norb_ex_cut^2 * nc_ex * nv_ex) explicit
    ! loop with ZGEMM calls, reducing cost from O(N^2 * B^2) to O(N^2 * B)
    ! where N = norb_ex_cut, B = nc_ex or nv_ex.
    !
    ! Strategy: for each k-point, precompute a flat index table i_ex(ic,iv)
    ! then build coefficient matrices F_c, F_v of shape (norb_ex_cut, nc*nv)
    ! and band matrices B_cc, B_vv of shape (nc, nc) or (nv, nv), and
    ! express all six accumulators as ZGEMM products.
    !
    !   vme_ex_inter1(nj,:,:) += F_c^dag · B_cc(nj) · F_c   [ZGEMM x2 per nj]
    !   vme_ex_inter2(nj,:,:) -= F_v^dag · B_vv(nj) · F_v   [ZGEMM x2 per nj]
    !   yme_ex_inter1/2: same, off-diagonal mask applied to B before GEMM
    !   qme_ex_inter1(nj,:,:) += i * F_c^dag · D_c(nj)       [ZGEMM, 1-sided]
    !   qme_ex_inter2(nj,:,:) += i * F_c^dag · A(nj)         [ZGEMM, 1-sided]
    !
    ! ZGEMM is multithreaded internally by OpenBLAS; the nj=1..3 loop and the
    ! k-point loop (external) provide additional coarse-grained parallelism.
    !--------------------------------------------------------------------------
    implicit none
    integer,    intent(in) :: ibz
    real(8),    intent(in) :: ek(npointstotal, nband_ex)
    complex(8), intent(in) :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in) :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)

    ! --- local working arrays ---
    integer :: ic, iv, icp, ivp, nj, i_ex, idx
    integer :: nbasis           ! nc_ex * nv_ex

    ! Flat index lookup: i_ex_table(ic,iv) = i_ex_nn for this ibz
    integer, allocatable :: i_ex_table(:,:)        ! (nc_ex, nv_ex)

    ! Coefficient matrices — shape (norb_ex_cut, nbasis), column = (ic,iv) pair
    ! F_left  = conjg(fk_ex) for left  factor (appears as F^dag in ZGEMM)
    ! F_right =       fk_ex  for right factor
    ! They are the same array viewed differently; we keep one and pass with 'C'/'N'.
    complex(8), allocatable :: F_cv(:,:)           ! (norb_ex_cut, nbasis)

    ! fk_ex_der reshaped: D_c(nj, nn, idx) = fk_ex_der(nj, i_ex_table(ic,iv), nnp)
    ! We build a (nbasis, norb_ex_cut) matrix for each nj
    complex(8), allocatable :: D_c(:,:)            ! (nbasis, norb_ex_cut)

    ! Berry diagonal shift matrix for qme_ex_inter2:
    ! A(nj, nn, idx) = fk_ex(i_ex, nnp) * (berry_c - berry_v)
    complex(8), allocatable :: A_c(:,:)            ! (nbasis, norb_ex_cut)

    ! Band matrices (nc_ex x nc_ex) and (nv_ex x nv_ex)
    complex(8), allocatable :: B_cc(:,:), B_vv(:,:)   ! for vme
    complex(8), allocatable :: Y_cc(:,:), Y_vv(:,:)   ! for yme (off-diagonal only)

    ! Intermediate matrix for two-step GEMM: tmp = F^dag · B  or  F^dag · D
    complex(8), allocatable :: tmp_cc(:,:)         ! (norb_ex_cut, nc_ex)
    complex(8), allocatable :: tmp_vv(:,:)         ! (norb_ex_cut, nv_ex)

    complex(8), parameter :: ci   = (0.0d0, 1.0d0)
    complex(8), parameter :: cone = (1.0d0, 0.0d0)
    complex(8), parameter :: czero= (0.0d0, 0.0d0)

    real(8) :: berry_shift

    nbasis = nc_ex * nv_ex

    allocate(i_ex_table(nc_ex, nv_ex))
    allocate(F_cv(norb_ex_cut, nbasis))
    allocate(D_c(nbasis, norb_ex_cut))
    allocate(A_c(nbasis, norb_ex_cut))
    allocate(B_cc(nc_ex, nc_ex))
    allocate(B_vv(nv_ex, nv_ex))
    allocate(Y_cc(nc_ex, nc_ex))
    allocate(Y_vv(nv_ex, nv_ex))
    allocate(tmp_cc(norb_ex_cut, nc_ex))
    allocate(tmp_vv(norb_ex_cut, nv_ex))

    !------------------------------------------------------------------
    ! Step 1: precompute flat index table for this k-point (serial,
    !         cheap — nc_ex * nv_ex calls, typically << 1000)
    !------------------------------------------------------------------
    do ic = 1, nc_ex
      do iv = 1, nv_ex
        call get_ex_index_first(nf, nv_ex, nc_ex, 0, ibz, i_ex_table(ic,iv), ic, iv)
      end do
    end do

    !------------------------------------------------------------------
    ! Step 2: build F_cv(nn, idx) = fk_ex(i_ex_table(ic,iv), nn)
    !         where idx = (ic-1)*nv_ex + iv  (column-major over iv)
    !------------------------------------------------------------------
    do ic = 1, nc_ex
      do iv = 1, nv_ex
        idx = (ic-1)*nv_ex + iv
        F_cv(:, idx) = fk_ex(i_ex_table(ic,iv), :)
      end do
    end do

    !------------------------------------------------------------------
    ! Step 3: loop over Cartesian components and accumulate via ZGEMM
    !------------------------------------------------------------------
    do nj = 1, 3

      !--- 3a. vme_ex_inter1: sum_{ic,icp,iv} conjg(F_c_nn) * V_cc(ic,icp) * F_c_nnp
      !        = F_cv^dag · B_cc · F_cv  (summed over iv implicitly via F_cv columns)
      !
      !   We contract iv first by building B_cc as sum over iv of the outer product,
      !   then do two GEMMs.  Since vme_ex_band(ibz,nj,ic,icp) is independent of iv,
      !   and F_cv sums over (ic,iv) columns, we fold iv into F by summing:
      !
      !     F_c_sum(nn, ic) = sum_iv F_cv(nn, (ic-1)*nv_ex+iv)   [collapse iv]
      !
      !   then  vme_ex_inter1(nj,:,:) += F_c_sum^dag · B_cc · F_c_sum
      !
      !   This is exact because B_cc does not depend on iv.
      !------------------------------------------------------------------

      ! Build F_c_sum by folding the iv dimension into tmp_cc
      tmp_cc = czero
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx = (ic-1)*nv_ex + iv
          tmp_cc(:, ic) = tmp_cc(:, ic) + F_cv(:, idx)
        end do
      end do
      ! tmp_cc(nn, ic) = sum_iv fk_ex(i_ex(ic,iv), nn)

      ! B_cc(ic,icp) = vme_ex_band(ibz, nj, nv_ex+ic, nv_ex+icp)
      do ic = 1, nc_ex
        do icp = 1, nc_ex
          B_cc(ic,icp) = vme_ex_band(ibz, nj, nv_ex+ic, nv_ex+icp)
        end do
      end do

      ! Y_cc = B_cc with diagonal zeroed (for yme, off-diagonal only)
      Y_cc = B_cc
      do ic = 1, nc_ex
        Y_cc(ic,ic) = czero
      end do

      ! Two-step ZGEMM:  mid = conjg(tmp_cc)^T · B_cc,  result += mid · tmp_cc^T
      ! i.e. (norb_ex_cut x nc_ex)^dag · (nc_ex x nc_ex) · (nc_ex x norb_ex_cut)
      ! Step (i): mid(nn, icp) = sum_ic conjg(tmp_cc(nn,ic)) * B_cc(ic,icp)
      !           => ZGEMM('C','N', norb_ex_cut, nc_ex, nc_ex, ...)
      ! Step (ii): result(nn,nnp) += sum_icp mid(nn,icp) * tmp_cc(nnp,icp)
      !           => ZGEMM('N','C', norb_ex_cut, norb_ex_cut, nc_ex, ...)
      !
      ! We reuse tmp_vv as scratch for the (norb_ex_cut x nc_ex) intermediate.
      ! Rename for clarity: use B_vv as the intermediate here (it is nc>=nv anyway).

      ! vme_ex_inter1
      call zgemm('C','N', norb_ex_cut, nc_ex, nc_ex, &
                  cone,  tmp_cc, norb_ex_cut, &
                         B_cc,   nc_ex, &
                  czero, tmp_vv, norb_ex_cut)
      call zgemm('N','C', norb_ex_cut, norb_ex_cut, nc_ex, &
                  cone,  tmp_vv, norb_ex_cut, &
                         tmp_cc, norb_ex_cut, &
                  cone,  vme_ex_inter1(nj,:,:), norb_ex_cut)

      ! yme_ex_inter1  (same but with Y_cc instead of B_cc)
      call zgemm('C','N', norb_ex_cut, nc_ex, nc_ex, &
                  cone,  tmp_cc, norb_ex_cut, &
                         Y_cc,   nc_ex, &
                  czero, tmp_vv, norb_ex_cut)
      call zgemm('N','C', norb_ex_cut, norb_ex_cut, nc_ex, &
                  cone,  tmp_vv, norb_ex_cut, &
                         tmp_cc, norb_ex_cut, &
                  cone,  yme_ex_inter1(nj,:,:), norb_ex_cut)

      !--- 3b. vme_ex_inter2: -(F_v_sum^dag · B_vv · F_v_sum) in valence space
      !        B_vv(iv,ivp) = vme_ex_band(ibz,nj,ivp,iv)  [note index order]

      ! Build F_v_sum: fold ic into columns
      ! F_v_sum(nn, iv) = sum_ic fk_ex(i_ex(ic,iv), nn)
      tmp_vv = czero
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx = (ic-1)*nv_ex + iv
          tmp_vv(:, iv) = tmp_vv(:, iv) + F_cv(:, idx)
        end do
      end do

      do iv = 1, nv_ex
        do ivp = 1, nv_ex
          B_vv(iv,ivp) = vme_ex_band(ibz, nj, ivp, iv)
        end do
      end do

      Y_vv = B_vv
      do iv = 1, nv_ex
        Y_vv(iv,iv) = czero
      end do

      ! Intermediate: (norb_ex_cut x nv_ex), reuse tmp_cc columns 1:nv_ex
      call zgemm('C','N', norb_ex_cut, nv_ex, nv_ex, &
                  cone,  tmp_vv, norb_ex_cut, &
                         B_vv,   nv_ex, &
                  czero, tmp_cc, norb_ex_cut)
      call zgemm('N','C', norb_ex_cut, norb_ex_cut, nv_ex, &
                  -cone, tmp_cc, norb_ex_cut, &
                         tmp_vv, norb_ex_cut, &
                  cone,  vme_ex_inter2(nj,:,:), norb_ex_cut)

      call zgemm('C','N', norb_ex_cut, nv_ex, nv_ex, &
                  cone,  tmp_vv, norb_ex_cut, &
                         Y_vv,   nv_ex, &
                  czero, tmp_cc, norb_ex_cut)
      call zgemm('N','C', norb_ex_cut, norb_ex_cut, nv_ex, &
                  -cone, tmp_cc, norb_ex_cut, &
                         tmp_vv, norb_ex_cut, &
                  cone,  yme_ex_inter2(nj,:,:), norb_ex_cut)

      !--- 3c. qme_ex_inter1: i * F_cv^dag · D_c(nj)
      !        D_c(idx, nnp) = fk_ex_der(nj, i_ex_table(ic,iv), nnp)
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx = (ic-1)*nv_ex + iv
          D_c(idx, :) = fk_ex_der(nj, i_ex_table(ic,iv), :)
        end do
      end do

      ! qme_ex_inter1(nj,nn,nnp) += i * sum_idx conjg(F_cv(nn,idx)) * D_c(idx,nnp)
      !   = i * (F_cv^dag · D_c)  shape (norb_ex_cut x norb_ex_cut)
      call zgemm('C','N', norb_ex_cut, norb_ex_cut, nbasis, &
                  ci,    F_cv, norb_ex_cut, &
                         D_c,  nbasis, &
                  cone,  qme_ex_inter1(nj,:,:), norb_ex_cut)

      !--- 3d. qme_ex_inter2: i * F_cv^dag · A_c(nj)
      !        A_c(idx, nnp) = fk_ex(i_ex, nnp) * (berry_cc - berry_vv) diagonal
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx = (ic-1)*nv_ex + iv
          i_ex = i_ex_table(ic,iv)
          berry_shift = dble( berry_eigen_ex_band(ibz,nj,nv_ex+ic,nv_ex+ic) &
                            - berry_eigen_ex_band(ibz,nj,iv,iv) )
          A_c(idx, :) = -fk_ex(i_ex, :) * berry_shift
        end do
      end do

      call zgemm('C','N', norb_ex_cut, norb_ex_cut, nbasis, &
                  ci,    F_cv, norb_ex_cut, &
                         A_c,  nbasis, &
                  cone,  qme_ex_inter2(nj,:,:), norb_ex_cut)

    end do ! nj

    deallocate(i_ex_table, F_cv, D_c, A_c, B_cc, B_vv, Y_cc, Y_vv, tmp_cc, tmp_vv)

  end subroutine get_ome_inter_ex_sum_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_ome_ex_linear(vme_ex)
    implicit none
    complex(8), intent(in) :: vme_ex(3, norb_ex_cut)
    integer :: nn, nj

    open(10, file='ome_linear_ex_'//trim(material_name)//'.omeex')
    write(10,*) 1
    do nn = 1, norb_ex_cut
      write(10,*) nn, (dble(vme_ex(nj,nn)), dimag(vme_ex(nj,nn)), nj=1,3)
    end do
    close(10)
  end subroutine write_ome_ex_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_ome_ex_linear_kresolved_init(unitno, mat_name, norb_ex_cut)
    implicit none
    integer,          intent(in) :: unitno, norb_ex_cut
    character(len=*), intent(in) :: mat_name

    open(unitno, file='ome_linear_ex_k_'//trim(mat_name)//'.omeexk', status='replace')
    write(unitno,*) 1
    write(unitno,*) norb_ex_cut
  end subroutine write_ome_ex_linear_kresolved_init

  subroutine write_ome_ex_linear_kresolved_point(unitno, kx, ky, kz, norb_ex_cut, vme_ex_k)
    implicit none
    integer,    intent(in) :: unitno, norb_ex_cut
    real(8),    intent(in) :: kx, ky, kz
    complex(8), intent(in) :: vme_ex_k(3, norb_ex_cut)
    integer :: nn

    write(unitno,*) kx, ky, kz
    do nn = 1, norb_ex_cut
      write(unitno,*) nn, dble(vme_ex_k(1,nn)), dimag(vme_ex_k(1,nn)), &
                          dble(vme_ex_k(2,nn)), dimag(vme_ex_k(2,nn)), &
                          dble(vme_ex_k(3,nn)), dimag(vme_ex_k(3,nn))
    end do
  end subroutine write_ome_ex_linear_kresolved_point

  subroutine write_ome_ex_linear_kresolved_close(unitno)
    implicit none
    integer, intent(in) :: unitno
    close(unitno)
  end subroutine write_ome_ex_linear_kresolved_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ome_sp_linear(iflag_norder, npointstotal, nband_ex, vme_ex_band, ek)
    implicit none
    integer, intent(in)  :: iflag_norder, npointstotal, nband_ex
    real(8), intent(out) :: ek(npointstotal, nband_ex)
    complex(8), intent(out) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)

    integer :: ibz, i, j, iflag_r
    real(8) :: a1, a2, a3, b1, b2, b3, b4, b5, b6

    open(10, file='ome_linear_sp_'//trim(material_name)//'.omesp')
    read(10,*) iflag_r
    do ibz = 1, npointstotal
      read(10,*) a1, a2, a3, (ek(ibz,j), j=1,nband_ex)
      do i = 1, nband_ex
        do j = 1, nband_ex
          read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
          vme_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
          vme_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
          vme_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
        end do
      end do
    end do
    close(10)
  end subroutine read_ome_sp_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, &
                                    berry_eigen_ex_band, gen_der_ex_band, &
                                    shift_vector_ex_band, vme_ex_band, ek)
    implicit none
    integer, intent(in)  :: iflag_norder, npointstotal, nband_ex
    real(8), intent(out) :: ek(npointstotal, nband_ex)
    real(8), intent(out) :: shift_vector_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: gen_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)

    integer :: ibz, i, j, nj, iflag_r
    real(8) :: a1, a2, a3, b1, b2, b3, b4, b5, b6

    open(10, file='ome_nonlinear_sp_'//trim(material_name)//'.omesp')
    read(10,*) iflag_r
    do ibz = 1, npointstotal
      read(10,*) a1, a2, a3, (ek(ibz,j), j=1,nband_ex)
      do i = 1, nband_ex
        do j = 1, nband_ex
          read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
          vme_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
          vme_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
          vme_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
          read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
          berry_eigen_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
          berry_eigen_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
          berry_eigen_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
          do nj = 1, 3
            read(10,*) a1, a2, a3, b1, b2, b3
            shift_vector_ex_band(ibz,nj,1,i,j) = b1
            shift_vector_ex_band(ibz,nj,2,i,j) = b2
            shift_vector_ex_band(ibz,nj,3,i,j) = b3
          end do
          do nj = 1, 3
            read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
            gen_der_ex_band(ibz,nj,1,i,j) = cmplx(b1,b2,8)
            gen_der_ex_band(ibz,nj,2,i,j) = cmplx(b3,b4,8)
            gen_der_ex_band(ibz,nj,3,i,j) = cmplx(b5,b6,8)
          end do
        end do
      end do
    end do
    close(10)
  end subroutine read_ome_sp_nonlinear

end module ome_ex