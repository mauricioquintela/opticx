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

  complex(8), allocatable :: xme_ex(:,:)
  complex(8), allocatable :: vme_ex(:,:)
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

    integer :: ibz, nn, nnp, nj, nbasis
    integer :: u_exk
    logical :: do_write_exk

    complex(8), allocatable :: vme_ex_k(:,:)
    real(8),    allocatable :: ek(:,:)
    complex(8), allocatable :: xme_ex_band(:,:,:,:)
    complex(8), allocatable :: vme_ex_band(:,:,:,:)
    complex(8), allocatable :: berry_eigen_ex_band(:,:,:,:)
    complex(8), allocatable :: gen_der_ex_band(:,:,:,:,:)
    real(8),    allocatable :: shift_vector_ex_band(:,:,:,:,:)

    ! Working arrays for get_ome_inter_ex_sum_k — allocated once outside k-loop
    integer,    allocatable :: i_ex_table(:,:)
    complex(8), allocatable :: F_cv(:,:)
    complex(8), allocatable :: FcvH(:,:)
    complex(8), allocatable :: D_c(:,:)
    complex(8), allocatable :: A_c(:,:)
    complex(8), allocatable :: B_cc(:,:)
    complex(8), allocatable :: Y_cc(:,:)
    complex(8), allocatable :: B_vv(:,:)
    complex(8), allocatable :: Y_vv(:,:)
    complex(8), allocatable :: Uc(:,:)
    complex(8), allocatable :: Wv(:,:)
    complex(8), allocatable :: mid_cc(:,:)
    complex(8), allocatable :: mid_vv(:,:)
    complex(8), allocatable :: out_v(:,:)
    complex(8), allocatable :: out_y(:,:)
    complex(8), allocatable :: out_q(:,:)

    write(*,*) '6. Entering ome_ex'
    do_write_exk = (iflag_norder == 1)

    allocate(ek(npointstotal, nband_ex))
    allocate(vme_ex_band(npointstotal, 3, nband_ex, nband_ex))

    ! PATCH: xme_ex_band is now allocated and zeroed UNCONDITIONALLY, not
    ! only inside the iflag_norder==2 branch. get_ome_gs_ex_sum_k reads
    ! xme_ex_band(ibz,nj,iv,nv_ex+ic) on every call regardless of
    ! iflag_norder, with no guard of its own -- for iflag_norder==1 that
    ! was previously a read of an UNALLOCATED allocatable array (undefined
    ! behaviour, most likely a crash). Allocating it here and leaving it
    ! zeroed for the linear case makes that read well-defined; the value
    ! it contributes (xme_ex accumulation) is simply never written out or
    ! otherwise used when iflag_norder==1, so this changes no output.
    allocate(xme_ex_band(npointstotal, 3, nband_ex, nband_ex))
    xme_ex_band = (0.0d0, 0.0d0)

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
      ! PATCH: xme_ex_band allocation moved above -- removed the duplicate
      ! `allocate(xme_ex_band(...))` that used to live here.
      ek = 0.0d0
      vme_ex_band          = (0.0d0, 0.0d0)
      berry_eigen_ex_band  = (0.0d0, 0.0d0)
      gen_der_ex_band      = (0.0d0, 0.0d0)
      shift_vector_ex_band = 0.0d0
      write(*,*) '   Reading optical matrix elements (sp)...'
      call read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, &
                                  berry_eigen_ex_band, gen_der_ex_band, &
                                  shift_vector_ex_band, vme_ex_band, ek)
      ! xme_ex_band already zeroed above; overwritten here with real values.
      call get_ome_sp_xme_ex_band(ek, vme_ex_band, xme_ex_band)
    end if

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

      nbasis = nc_ex * nv_ex
      allocate(i_ex_table(nc_ex,  nv_ex))
      allocate(F_cv (norb_ex_cut, nbasis))
      allocate(FcvH (nbasis,      norb_ex_cut))
      allocate(D_c  (nbasis,      norb_ex_cut))
      allocate(A_c  (nbasis,      norb_ex_cut))
      allocate(B_cc (nc_ex,       nc_ex))
      allocate(Y_cc (nc_ex,       nc_ex))
      allocate(B_vv (nv_ex,       nv_ex))
      allocate(Y_vv (nv_ex,       nv_ex))
      allocate(Uc   (nc_ex,       norb_ex_cut))
      allocate(Wv   (nv_ex,       norb_ex_cut))
      allocate(mid_cc(nc_ex,      norb_ex_cut))
      allocate(mid_vv(nv_ex,      norb_ex_cut))
      allocate(out_v(norb_ex_cut, norb_ex_cut))
      allocate(out_y(norb_ex_cut, norb_ex_cut))
      allocate(out_q(norb_ex_cut, norb_ex_cut))
    end if

    if (iflag_norder == 1) write(*,*) '   Evaluating excitonic OMEs for linear conductivity...'
    if (iflag_norder == 2) write(*,*) '   Evaluating excitonic OMEs for nonlinear conductivity...'

    do ibz = 1, npointstotal
      write(*,*) '   OME (ex): k-point', ibz, '/', npointstotal

      if (iflag_norder == 1 .or. iflag_norder == 2) &
        call get_ome_gs_ex_sum_k(ibz, ek, vme_ex_band, xme_ex_band)

      if (do_write_exk) then
        call get_ome_gs_ex_kresolved(ibz, ek, vme_ex_band, vme_ex_k)
        call write_ome_ex_linear_kresolved_point(u_exk, rkxvector(ibz), rkyvector(ibz), &
                                                  rkzvector(ibz), norb_ex_cut, vme_ex_k)
      end if

      if (iflag_norder == 2) &
        call get_ome_inter_ex_sum_k(                                        &
               ibz, ek, xme_ex_band, vme_ex_band, berry_eigen_ex_band,     &
               nbasis, i_ex_table, F_cv, FcvH, D_c, A_c,                   &
               B_cc, Y_cc, B_vv, Y_vv, Uc, Wv, mid_cc, mid_vv,             &
               out_v, out_y, out_q)
    end do

    if (do_write_exk) then
      call write_ome_ex_linear_kresolved_close(u_exk)
      deallocate(vme_ex_k)
      write(*,*) '   k-resolved linear excitonic OMEs written (omeexk)'
    end if

    if (iflag_norder == 2) then
      allocate(xme_ex_inter(3, norb_ex_cut, norb_ex_cut))
      qme_ex_inter = qme_ex_inter1 + qme_ex_inter2
      yme_ex_inter = yme_ex_inter1 + yme_ex_inter2
      xme_ex_inter = yme_ex_inter  + qme_ex_inter
      vme_ex_inter = vme_ex_inter1 + vme_ex_inter2
      deallocate(i_ex_table, F_cv, FcvH, D_c, A_c)
      deallocate(B_cc, Y_cc, B_vv, Y_vv, Uc, Wv, mid_cc, mid_vv)
      deallocate(out_v, out_y, out_q)
      deallocate(fk_ex_der)
    end if

    write(*,*) '   Optical matrix elements (ex) have been evaluated'
    if (iflag_norder == 1) call write_ome_ex_linear(vme_ex)
    write(*,*) '   Optical matrix elements (ex, N->GS) written'
    write(*,*) '   Optical matrix elements (ex, N->N'') not printed in this version'

    ! PATCH: xme_ex_band is now unconditionally allocated above, so it
    ! must be unconditionally deallocated here too (previously this was
    ! only deallocated when iflag_norder==2).
    deallocate(ek, vme_ex_band, xme_ex_band)
    if (iflag_norder == 2) &
      deallocate(berry_eigen_ex_band, gen_der_ex_band, shift_vector_ex_band)

  end subroutine get_ome_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_sp_xme_ex_band(ek, vme_ex_band, xme_ex_band)
    implicit none
    real(8),    intent(in)  :: ek(npointstotal, nband_ex)
    complex(8), intent(in)  :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(out) :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    integer    :: ibz, i, j, nj
    real(8)    :: de
    complex(8), parameter :: ci = (0.0d0, 1.0d0)

    xme_ex_band = (0.0d0, 0.0d0)
    !$omp parallel do collapse(2) schedule(static) private(ibz, nj, i, j, de)
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
  subroutine get_ome_gs_ex_sum_k(ibz, ek, vme_ex_band, xme_ex_band)
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
            vme_ex(nj,nn) = vme_ex(nj,nn) + fk_ex(i_ex_nn,nn)*vme_ex_band(ibz,nj,iv,nv_ex+ic)
            xme_ex(nj,nn) = xme_ex(nj,nn) + fk_ex(i_ex_nn,nn)*xme_ex_band(ibz,nj,iv,nv_ex+ic)
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
            vme_ex_k(nj,nn) = vme_ex_k(nj,nn) + fk_ex(i_ex_nn,nn)*vme_ex_band(ibz,nj,iv,nv_ex+ic)
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine get_ome_gs_ex_kresolved

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get_ome_inter_ex_sum_k
!
! IMPORTANT: the (ic,icp) contraction against B_cc/Y_cc must be carried out
! separately for EACH value of iv, and the (iv,ivp) contraction against
! B_vv/Y_vv must be carried out separately for EACH value of ic, and only
! THEN summed. iv (resp. ic) is a shared "spectator" index tying together
! both matrix-element factors in the original sum
!
!   V1_{nn,nn'} = sum_iv  sum_{ic,icp} conjg(f_{(ic,iv),nn}) B_cc(ic,icp) f_{(icp,iv),nn'}
!   V2_{nn,nn'} = -sum_ic sum_{iv,ivp} conjg(f_{(ic,iv),nn}) B_vv(iv,ivp) f_{(ic,ivp),nn'}
!
! It is NOT valid to pre-sum fk_ex over iv (resp. ic) and then form a single
! bilinear contraction from the pre-summed vectors: that silently introduces
! spurious "cross" terms with iv /= iv' (resp. ic /= ic') that do not belong
! in the physical sum. Each iv (resp. ic) therefore requires its own small
! ZGEMM contraction, accumulated into the output with beta = 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_ome_inter_ex_sum_k(                                      &
               ibz, ek, xme_ex_band, vme_ex_band, berry_eigen_ex_band,   &
               nbasis, i_ex_table, F_cv, FcvH, D_c, A_c,                 &
               B_cc, Y_cc, B_vv, Y_vv, Uc, Wv, mid_cc, mid_vv,           &
               out_v, out_y, out_q)
    implicit none

    integer,    intent(in)    :: ibz
    integer,    intent(in)    :: nbasis
    real(8),    intent(in)    :: ek(npointstotal, nband_ex)
    complex(8), intent(in)    :: xme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in)    :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
    complex(8), intent(in)    :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
    integer,    intent(inout) :: i_ex_table(nc_ex, nv_ex)
    complex(8), intent(inout) :: F_cv(norb_ex_cut, nbasis)
    complex(8), intent(inout) :: FcvH(nbasis, norb_ex_cut)
    complex(8), intent(inout) :: D_c(nbasis, norb_ex_cut)
    complex(8), intent(inout) :: A_c(nbasis, norb_ex_cut)
    complex(8), intent(inout) :: B_cc(nc_ex, nc_ex)
    complex(8), intent(inout) :: Y_cc(nc_ex, nc_ex)
    complex(8), intent(inout) :: B_vv(nv_ex, nv_ex)
    complex(8), intent(inout) :: Y_vv(nv_ex, nv_ex)
    complex(8), intent(inout) :: Uc(nc_ex, norb_ex_cut)
    complex(8), intent(inout) :: Wv(nv_ex, norb_ex_cut)
    complex(8), intent(inout) :: mid_cc(nc_ex, norb_ex_cut)
    complex(8), intent(inout) :: mid_vv(nv_ex, norb_ex_cut)
    complex(8), intent(inout) :: out_v(norb_ex_cut, norb_ex_cut)
    complex(8), intent(inout) :: out_y(norb_ex_cut, norb_ex_cut)
    complex(8), intent(inout) :: out_q(norb_ex_cut, norb_ex_cut)

    integer    :: ic, iv, icp, ivp, nj, i_ex, idx, nn
    real(8)    :: berry_shift
    complex(8), parameter :: ci    = (0.0d0, 1.0d0)
    complex(8), parameter :: cone  = (1.0d0, 0.0d0)
    complex(8), parameter :: czero = (0.0d0, 0.0d0)

    !--- Step 1: index table ---
    do ic = 1, nc_ex
      do iv = 1, nv_ex
        call get_ex_index_first(nf, nv_ex, nc_ex, 0, ibz, i_ex_table(ic,iv), ic, iv)
      end do
    end do

    !--- Step 2: F_cv(nn,idx) / FcvH(idx,nn) over the full (ic,iv) basis ---
    ! (used only for the qme_ex_inter1/2 terms, which contract over the
    !  FULL basis index and do not have the spectator-index restriction)
    do ic = 1, nc_ex
      do iv = 1, nv_ex
        idx  = (ic-1)*nv_ex + iv
        i_ex = i_ex_table(ic,iv)
        do nn = 1, norb_ex_cut
          F_cv(nn, idx) = fk_ex(i_ex, nn)
        end do
      end do
    end do
    do nn = 1, norb_ex_cut
      do idx = 1, nbasis
        FcvH(idx, nn) = conjg(F_cv(nn, idx))
      end do
    end do

    do nj = 1, 3

      !================================================================
      ! vme_ex_inter1 / yme_ex_inter1:
      !   sum_iv  Uc(iv)^H * B_cc * Uc(iv)   (and with Y_cc for yme)
      ! where Uc(iv)(ic,nn) = fk_ex(i_ex_table(ic,iv), nn)
      !================================================================
      ! B_cc holds the velocity matrix elements (feeds vme_ex_inter1);
      ! Y_cc holds the POSITION matrix elements xme_ex_band (feeds
      ! yme_ex_inter1) — these are physically distinct quantities, Y_cc
      ! is NOT simply B_cc with the diagonal removed. xme_ex_band already
      ! has an exactly-zero diagonal by construction (see
      ! get_ome_sp_xme_ex_band), so no further zeroing is strictly needed;
      ! it is kept below only for clarity/parity with the reference code.
      do ic = 1, nc_ex
        do icp = 1, nc_ex
          B_cc(ic,icp) = vme_ex_band(ibz, nj, nv_ex+ic, nv_ex+icp)
          Y_cc(ic,icp) = xme_ex_band(ibz, nj, nv_ex+ic, nv_ex+icp)
        end do
      end do
      do ic = 1, nc_ex
        Y_cc(ic,ic) = czero
      end do

      out_v = czero
      out_y = czero
      do iv = 1, nv_ex
        do ic = 1, nc_ex
          Uc(ic,:) = fk_ex(i_ex_table(ic,iv), :)
        end do

        ! mid_cc = B_cc * Uc ; out_v += Uc^H * mid_cc
        call zgemm('N', 'N', nc_ex, norb_ex_cut, nc_ex, &
                    cone,  B_cc, nc_ex, &
                           Uc,   nc_ex, &
                    czero, mid_cc, nc_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nc_ex, &
                    cone, Uc,     nc_ex, &
                          mid_cc, nc_ex, &
                    cone, out_v,  norb_ex_cut)

        ! mid_cc = Y_cc * Uc ; out_y += Uc^H * mid_cc
        call zgemm('N', 'N', nc_ex, norb_ex_cut, nc_ex, &
                    cone,  Y_cc, nc_ex, &
                           Uc,   nc_ex, &
                    czero, mid_cc, nc_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nc_ex, &
                    cone, Uc,     nc_ex, &
                          mid_cc, nc_ex, &
                    cone, out_y,  norb_ex_cut)
      end do
      vme_ex_inter1(nj,:,:) = vme_ex_inter1(nj,:,:) + out_v
      yme_ex_inter1(nj,:,:) = yme_ex_inter1(nj,:,:) + out_y

      !================================================================
      ! vme_ex_inter2 / yme_ex_inter2:
      !   -sum_ic  Wv(ic)^H * B_vv * Wv(ic)  (and with Y_vv for yme)
      ! where Wv(ic)(iv,nn) = fk_ex(i_ex_table(ic,iv), nn)
      !================================================================
      ! Same distinction as above: B_vv is the velocity matrix (feeds
      ! vme_ex_inter2), Y_vv is the POSITION matrix xme_ex_band (feeds
      ! yme_ex_inter2) — not a zero-diagonal copy of B_vv.
      do ivp = 1, nv_ex
        do iv = 1, nv_ex
          B_vv(iv,ivp) = vme_ex_band(ibz, nj, ivp, iv)
          Y_vv(iv,ivp) = xme_ex_band(ibz, nj, ivp, iv)
        end do
      end do
      do iv = 1, nv_ex
        Y_vv(iv,iv) = czero
      end do

      out_v = czero
      out_y = czero
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          Wv(iv,:) = fk_ex(i_ex_table(ic,iv), :)
        end do

        ! mid_vv = B_vv * Wv ; out_v += Wv^H * mid_vv
        call zgemm('N', 'N', nv_ex, norb_ex_cut, nv_ex, &
                    cone,  B_vv, nv_ex, &
                           Wv,   nv_ex, &
                    czero, mid_vv, nv_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nv_ex, &
                    cone, Wv,     nv_ex, &
                          mid_vv, nv_ex, &
                    cone, out_v,  norb_ex_cut)

        ! mid_vv = Y_vv * Wv ; out_y += Wv^H * mid_vv
        call zgemm('N', 'N', nv_ex, norb_ex_cut, nv_ex, &
                    cone,  Y_vv, nv_ex, &
                           Wv,   nv_ex, &
                    czero, mid_vv, nv_ex)
        call zgemm('C', 'N', norb_ex_cut, norb_ex_cut, nv_ex, &
                    cone, Wv,     nv_ex, &
                          mid_vv, nv_ex, &
                    cone, out_y,  norb_ex_cut)
      end do
      vme_ex_inter2(nj,:,:) = vme_ex_inter2(nj,:,:) - out_v
      yme_ex_inter2(nj,:,:) = yme_ex_inter2(nj,:,:) - out_y

      !================================================================
      ! qme_ex_inter1 += i * FcvH^T * D_c
      !================================================================
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx = (ic-1)*nv_ex + iv
          D_c(idx, :) = fk_ex_der(nj, i_ex_table(ic,iv), :)
        end do
      end do
      call zgemm('T', 'N', norb_ex_cut, norb_ex_cut, nbasis, &
                  ci,    FcvH, nbasis,       &
                         D_c,  nbasis,        &
                  czero, out_q, norb_ex_cut)
      qme_ex_inter1(nj,:,:) = qme_ex_inter1(nj,:,:) + out_q

      !================================================================
      ! qme_ex_inter2 += i * FcvH^T * A_c,  A_c = -i * fk_ex * berry_shift
      ! (the "-i" is essential — dropping it rotates this term by 90 deg)
      !================================================================
      do ic = 1, nc_ex
        do iv = 1, nv_ex
          idx  = (ic-1)*nv_ex + iv
          i_ex = i_ex_table(ic,iv)
          berry_shift = dble(berry_eigen_ex_band(ibz,nj,nv_ex+ic,nv_ex+ic)) &
                      - dble(berry_eigen_ex_band(ibz,nj,iv,iv))
          A_c(idx, :) = -ci * fk_ex(i_ex, :) * berry_shift
        end do
      end do
      call zgemm('T', 'N', norb_ex_cut, norb_ex_cut, nbasis, &
                  ci,    FcvH, nbasis,       &
                         A_c,  nbasis,        &
                  czero, out_q, norb_ex_cut)
      qme_ex_inter2(nj,:,:) = qme_ex_inter2(nj,:,:) + out_q

    end do ! nj

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
    integer,    intent(in)  :: iflag_norder, npointstotal, nband_ex
    real(8),    intent(out) :: ek(npointstotal, nband_ex)
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
   integer :: ibz, iflag_r, npts_r, nband_r
   real(8), allocatable :: rkx_r(:), rky_r(:), rkz_r(:)

   open(10, file='ome_nonlinear_sp_'//trim(material_name)//'.omesp', &
        form='unformatted', access='stream', status='old')

   read(10) iflag_r
   read(10) npts_r, nband_r
   allocate(rkx_r(npts_r), rky_r(npts_r), rkz_r(npts_r))
   read(10) rkx_r, rky_r, rkz_r

   do ibz = 1, npointstotal
      read(10) ek(ibz,:)
      read(10) vme_ex_band(ibz,:,:,:)
      read(10) berry_eigen_ex_band(ibz,:,:,:)
      read(10) shift_vector_ex_band(ibz,:,:,:,:)
      read(10) gen_der_ex_band(ibz,:,:,:,:)
   end do

   close(10)
  end subroutine read_ome_sp_nonlinear
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine read_ome_sp_nonlinear(iflag_norder, npointstotal, nband_ex, &
!                                     berry_eigen_ex_band, gen_der_ex_band, &
!                                     shift_vector_ex_band, vme_ex_band, ek)
!     implicit none
!     integer, intent(in)  :: iflag_norder, npointstotal, nband_ex
!     real(8), intent(out) :: ek(npointstotal, nband_ex)
!     real(8), intent(out) :: shift_vector_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
!     complex(8), intent(out) :: vme_ex_band(npointstotal, 3, nband_ex, nband_ex)
!     complex(8), intent(out) :: berry_eigen_ex_band(npointstotal, 3, nband_ex, nband_ex)
!     complex(8), intent(out) :: gen_der_ex_band(npointstotal, 3, 3, nband_ex, nband_ex)
!     integer :: ibz, i, j, nj, iflag_r
!     real(8) :: a1, a2, a3, b1, b2, b3, b4, b5, b6
!     open(10, file='ome_nonlinear_sp_'//trim(material_name)//'.omesp')
!     read(10,*) iflag_r
!     do ibz = 1, npointstotal
!       read(10,*) a1, a2, a3, (ek(ibz,j), j=1,nband_ex)
!       do i = 1, nband_ex
!         do j = 1, nband_ex
!           read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
!           vme_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
!           vme_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
!           vme_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
!           read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
!           berry_eigen_ex_band(ibz,1,i,j) = cmplx(b1,b2,8)
!           berry_eigen_ex_band(ibz,2,i,j) = cmplx(b3,b4,8)
!           berry_eigen_ex_band(ibz,3,i,j) = cmplx(b5,b6,8)
!           do nj = 1, 3
!             read(10,*) a1, a2, a3, b1, b2, b3
!             shift_vector_ex_band(ibz,nj,1,i,j) = b1
!             shift_vector_ex_band(ibz,nj,2,i,j) = b2
!             shift_vector_ex_band(ibz,nj,3,i,j) = b3
!           end do
!           do nj = 1, 3
!             read(10,*) a1, a2, a3, b1, b2, b3, b4, b5, b6
!             gen_der_ex_band(ibz,nj,1,i,j) = cmplx(b1,b2,8)
!             gen_der_ex_band(ibz,nj,2,i,j) = cmplx(b3,b4,8)
!             gen_der_ex_band(ibz,nj,3,i,j) = cmplx(b5,b6,8)
!           end do
!         end do
!       end do
!     end do
!     close(10)
!   end subroutine read_ome_sp_nonlinear

end module ome_ex