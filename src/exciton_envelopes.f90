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
    real(8)  :: rk1, rk2, rk3, rkxp, rkyp, rkzp
    real(8)  :: rkxp_bz, rkyp_bz, rkzp_bz, rk1_bz, rk2_bz, rk3_bz
    complex(8) :: fk_for, fk_back

    ! PATCH: compute once, up front, instead of locally on every call.
    if (.not. active_flags_set) call set_active_flags()

    do ibz = 1, npointstotal
      call get_k_kc(G, rkxvector(ibz), rkyvector(ibz), rkzvector(ibz), &
                    rk1vector(ibz), rk2vector(ibz), rk3vector(ibz), &
                    xp_bz, yp_bz, zp_bz, xc_bz, yc_bz, zc_bz)
    end do

    write(*,*) '   Evaluating exciton envelope function derivative with respect to k...'

    !$omp parallel do collapse(2) schedule(dynamic) &
    !$omp private(j, nn, ibz, nj, j_aux, &
    !$omp         rkxp, rkyp, rkzp, &
    !$omp         rk1, rk2, rk3, &
    !$omp         rkxp_bz, rkyp_bz, rkzp_bz, &
    !$omp         rk1_bz, rk2_bz, rk3_bz, &
    !$omp         fk_for, fk_back, &
    !$omp         fk_ex_col, fk_der_col)
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
          rkxp = rkxvector(ibz)
          rkyp = rkyvector(ibz)
          rkzp = rkzvector(ibz)
          call get_k_kc(G, rkxp, rkyp, rkzp, rk1, rk2, rk3, &
                        rkxp_bz, rkyp_bz, rkzp_bz, rk1_bz, rk2_bz, rk3_bz)

          if (ndim == 1) then

            if (active_z) then
              if (abs(abs(rk3)-0.5d0) < 1.0d-5) cycle
              call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)+dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)-dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(3,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            elseif (active_y) then
              if (abs(abs(rk2)-0.5d0) < 1.0d-5) cycle
              call get_k_kc(G,rkxp,rkyvector(ibz)+dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxp,rkyvector(ibz)-dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(2,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            else
              if (abs(abs(rk1)-0.5d0) < 1.0d-5) cycle
              call get_k_kc(G,rkxvector(ibz)+dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxvector(ibz)-dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(1,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            end if

          elseif (ndim == 2) then

            if (active_x .and. active_y) then
              if (abs(abs(rk1)-0.5d0) < 1.0d-5 .or. abs(abs(rk2)-0.5d0) < 1.0d-5) cycle
              call get_k_kc(G,rkxvector(ibz)+dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxvector(ibz)-dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(1,ibz) = (fk_for-fk_back)/(2.0d0*dk)
              call get_k_kc(G,rkxp,rkyvector(ibz)+dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxp,rkyvector(ibz)-dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(2,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            elseif (active_x .and. active_z) then
              if (abs(abs(rk1)-0.5d0) < 1.0d-5 .or. abs(abs(rk3)-0.5d0) < 1.0d-5) cycle
              call get_k_kc(G,rkxvector(ibz)+dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxvector(ibz)-dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(1,ibz) = (fk_for-fk_back)/(2.0d0*dk)
              call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)+dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)-dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(3,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            else
              if (abs(abs(rk2)-0.5d0) < 1.0d-5 .or. abs(abs(rk3)-0.5d0) < 1.0d-5) cycle
              call get_k_kc(G,rkxp,rkyvector(ibz)+dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxp,rkyvector(ibz)-dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(2,ibz) = (fk_for-fk_back)/(2.0d0*dk)
              call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)+dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
              call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)-dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
              call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
              fk_der_col(3,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            end if

          else  ! ndim == 3

            if (abs(abs(rk1)-0.5d0) < 1.0d-5 .or. &
                abs(abs(rk2)-0.5d0) < 1.0d-5 .or. &
                abs(abs(rk3)-0.5d0) < 1.0d-5) cycle
            call get_k_kc(G,rkxvector(ibz)+dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
            call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
            call get_k_kc(G,rkxvector(ibz)-dk,rkyp,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
            call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
            fk_der_col(1,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            call get_k_kc(G,rkxp,rkyvector(ibz)+dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
            call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
            call get_k_kc(G,rkxp,rkyvector(ibz)-dk,rkzp,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
            call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
            fk_der_col(2,ibz) = (fk_for-fk_back)/(2.0d0*dk)
            call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)+dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
            call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_for)
            call get_k_kc(G,rkxp,rkyp,rkzvector(ibz)-dk,rk1,rk2,rk3,rkxp_bz,rkyp_bz,rkzp_bz,rk1_bz,rk2_bz,rk3_bz)
            call interp1(npointstotal,rk1vector,rk2vector,rk3vector,fk_ex_col,rk1,rk2,rk3,fk_back)
            fk_der_col(3,ibz) = (fk_for-fk_back)/(2.0d0*dk)

          end if
        end do ! ibz

        do ibz = 1, npointstotal
          j_aux = (ibz-1)*norb_ex_band + j
          do nj = 1, 3
            fk_ex_der(nj, j_aux, nn) = fk_der_col(nj, ibz)
          end do
        end do

        deallocate(fk_ex_col, fk_der_col)

      end do ! nn
    end do ! j
    !$omp end parallel do

  end subroutine get_fk_ex_der_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PATCH: fk_col was passed directly as a rank-1 actual argument to a
  ! rank-2 explicit-shape dummy in get_fk_ex_k_interp. Since
  ! get_fk_ex_k_interp is a module procedure (implicit explicit interface),
  ! Fortran requires matching ranks -- the original would likely fail to
  ! compile ("Rank mismatch in argument") or misbehave under lax argument
  ! checking. Fixed by explicitly reshaping into a genuine rank-2 (n,1)
  ! temporary before the call.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine interp1(npointstotal, rk1v, rk2v, rk3v, fk_col, rk1, rk2, rk3, result)
    implicit none
    integer,    intent(in)  :: npointstotal
    real(8),    intent(in)  :: rk1v(npointstotal), rk2v(npointstotal), rk3v(npointstotal)
    complex(8), intent(in)  :: fk_col(npointstotal)
    real(8),    intent(in)  :: rk1, rk2, rk3
    complex(8), intent(out) :: result
    complex(8) :: fk_col2(npointstotal,1)
    fk_col2(:,1) = fk_col
    call get_fk_ex_k_interp(npointstotal, rk1v, rk2v, rk3v, &
                             1, 1, fk_col2, rk1, rk2, rk3, 1, result)
  end subroutine interp1

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
        xc=(xp*G(3,3)-G(3,1)*yp)/(G(1,1)*G(3,3)-G(3,1)*G(1,3))
        yc=0.0d0
        zc=(G(1,1)*yp-xp*G(1,3))/(G(1,1)*G(3,3)-G(3,1)*G(1,3))
        xc_bz=xc-dble(int(xc/0.5d0)); yc_bz=0.0d0; zc_bz=zc-dble(int(zc/0.5d0))
        xp_bz=xc_bz*G(1,1)+yc_bz*G(3,1); yp_bz=0.0d0; zp_bz=xc_bz*G(1,3)+yc_bz*G(3,3)
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
      xc=(xp*(G(2,2)*G(3,3)-G(2,3)*G(3,2))+yp*(G(1,3)*G(3,2)-G(1,2)*G(3,3))+zp*(G(1,2)*G(2,3)-G(1,3)*G(2,2)))/det
      yc=(xp*(G(2,3)*G(3,1)-G(2,1)*G(3,3))+yp*(G(1,1)*G(3,3)-G(1,3)*G(3,1))+zp*(G(1,3)*G(2,1)-G(1,1)*G(2,3)))/det
      zc=(xp*(G(2,1)*G(3,2)-G(2,2)*G(3,1))+yp*(G(1,2)*G(3,1)-G(1,1)*G(3,2))+zp*(G(1,1)*G(2,2)-G(1,2)*G(2,1)))/det
      xc_bz=xc-dble(int(xc/0.5d0)); yc_bz=yc-dble(int(yc/0.5d0)); zc_bz=zc-dble(int(zc/0.5d0))
      xp_bz=xc_bz*G(1,1)+yc_bz*G(2,1)+zc_bz*G(3,1)
      yp_bz=xc_bz*G(1,2)+yc_bz*G(2,2)+zc_bz*G(3,2)
      zp_bz=xc_bz*G(1,3)+yc_bz*G(2,3)+zc_bz*G(3,3)
    end if

  end subroutine get_k_kc

end module exciton_envelopes