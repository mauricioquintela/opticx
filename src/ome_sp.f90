module ome_sp
   use constants_math
   use parser_wannier90_tb, &
      only:material_name,nR,nRvec,norb,R,shop,hhop,rhop_c
   use parser_optics_xatu_dim, &
      only:npointstotal,rkxvector,rkyvector,rkzvector, &
      nband_ex,nband_index,nv_ex,nc_ex

   implicit none

   ! PATCH: dimensionality flags computed ONCE (via set_active_flags) instead
   ! of being recomputed inside every call to get_vme_kernels_ome /
   ! get_berry_eigen_fourpoint. Read-only after initialization, safe to
   ! share across OMP threads without being PRIVATE.
   logical, save :: active_x = .false., active_y = .false., active_z = .false.
   logical, save :: active_flags_set = .false.

   ! PATCH: named clipping threshold instead of a magic literal duplicated
   ! in two places (berry_eigen and shift_vector divergence clipping).
   real(8), parameter :: clip_threshold = 50.0d0
   
      real(8), allocatable, save :: Rx_global(:), Ry_global(:), Rz_global(:)
      logical, save :: R_cache_set = .false.


contains
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_R_cache()
            implicit none
            integer :: iRp
            if (.not. active_flags_set) call set_active_flags()
                  allocate(Rx_global(nR), Ry_global(nR), Rz_global(nR))
            do iRp = 1, nR
                  if (active_x) then
                  Rx_global(iRp) = dble(nRvec(iRp,1))*R(1,1) + dble(nRvec(iRp,2))*R(2,1) + dble(nRvec(iRp,3))*R(3,1)
                  else
                  Rx_global(iRp) = 0.0d0
                  end if
                  if (active_y) then
                  Ry_global(iRp) = dble(nRvec(iRp,1))*R(1,2) + dble(nRvec(iRp,2))*R(2,2) + dble(nRvec(iRp,3))*R(3,2)
                  else
                  Ry_global(iRp) = 0.0d0
                  end if
                  if (active_z) then
                  Rz_global(iRp) = dble(nRvec(iRp,1))*R(1,3) + dble(nRvec(iRp,2))*R(2,3) + dble(nRvec(iRp,3))*R(3,3)
                  else
                  Rz_global(iRp) = 0.0d0
                  end if
            end do
            R_cache_set = .true.
      end subroutine set_R_cache
      
   subroutine set_active_flags()
      implicit none
      active_x = (NORM2(real(nRvec(:,1))) /= 0.0d0)
      active_y = (NORM2(real(nRvec(:,2))) /= 0.0d0)
      active_z = (NORM2(real(nRvec(:,3))) /= 0.0d0)
      active_flags_set = .true.
   end subroutine set_active_flags
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_ome_sp(iflag_norder)
      implicit none

      integer iflag_norder
      integer ibz
      integer i,j,ii,jj,nj

      complex*16, allocatable :: skernel(:,:), hkernel(:,:)
      complex*16, allocatable :: sderkernel(:,:,:), hderkernel(:,:,:)
      complex*16, allocatable :: akernel(:,:,:)

      real(8)    :: e(norb)
      complex*16 :: hk_ev(norb,norb)
      complex*16 :: vme(norb,norb,3)

      complex*16 :: abc(norb,norb,3)
      real(8)    :: shift_vector(norb,norb,3,3)
      complex*16 :: berry_eigen1(norb,norb,3), berry_eigen2(norb,norb,3)
      complex*16 :: berry_eigen(norb,norb,3)

      integer :: kmoment

      complex*16, allocatable :: gen_der(:,:,:,:), vme_der(:,:,:,:)
      complex*16, allocatable :: gd1(:,:,:,:), gd2(:,:,:,:), gd3(:,:,:,:)
      complex*16, allocatable :: hk_ev_neigh(:,:,:), vme_neigh(:,:,:,:)

      ! NEW: scratch for the zgemm-based basis change inside get_vme_eigen_ome.
      ! Allocated once per thread (same lifetime/placement as gd1,gd2,gd3),
      ! reused across every call to get_vme_eigen_ome made by this thread —
      ! both the "main" call below and the 7 calls inside
      ! get_berry_eigen_fourpoint.
      complex*16, allocatable :: M1(:,:), T1(:,:)

      real(8), allocatable :: vme_der_phase(:,:,:,:)

      real(8),    allocatable :: ek(:,:)
      complex*16, allocatable :: vme_ex_band(:,:,:,:)
      complex*16, allocatable :: berry_eigen_ex_band(:,:,:,:)
      complex*16, allocatable :: gen_der_ex_band(:,:,:,:,:)
      real(8),    allocatable :: shift_vector_ex_band(:,:,:,:,:)

      real(8) :: rkx,rky,rkz
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '5. Entering ome_sp'

      call set_active_flags()
      if (.not. R_cache_set) call set_R_cache()

      allocate(vme_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(ek(npointstotal,nband_ex))
      allocate(berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex))
      allocate(gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      allocate(shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex))
      gen_der_ex_band=0.0d0
      shift_vector_ex_band=0.0d0
      berry_eigen_ex_band=0.0d0
      vme_ex_band=0.0d0
      ek=0.0d0

      write(*,*) '   Calculating optical matrix elements (sp): sampling BZ...'

      kmoment = -1

      !$OMP PARALLEL PRIVATE(rkx,rky,rkz,ibz,i,j,ii,jj,nj), &
      !$OMP PRIVATE(hkernel,skernel,sderkernel,hderkernel,akernel), &
      !$OMP PRIVATE(hk_ev,e,vme), &
      !$OMP PRIVATE(abc,gen_der,gd1,gd2,gd3), &
      !$OMP PRIVATE(vme_der,shift_vector,berry_eigen1,berry_eigen2,berry_eigen), &
      !$OMP PRIVATE(hk_ev_neigh,vme_neigh,vme_der_phase)&
      !$OMP PRIVATE(M1,T1)&                                    ! NEW
      !$OMP   SHARED(kmoment)

      allocate(skernel(norb,norb), hkernel(norb,norb))
      allocate(sderkernel(norb,norb,3), hderkernel(norb,norb,3))
      allocate(akernel(norb,norb,3))
      allocate(gen_der(norb,norb,3,3), vme_der(norb,norb,3,3))
      allocate(gd1(norb,norb,3,3), gd2(norb,norb,3,3), gd3(norb,norb,3,3))
      allocate(hk_ev_neigh(norb,norb,7), vme_neigh(norb,norb,3,7))
      allocate(vme_der_phase(norb,norb,3,3))
      allocate(M1(norb,norb), T1(norb,norb))                    ! NEW

      !$OMP DO SCHEDULE(STATIC)
      do ibz=1,npointstotal
            write(*,*) '   Optical matrix elements (sp): k-point',ibz,'/',npointstotal
            rkx=rkxvector(ibz)
            rky=rkyvector(ibz)
            rkz=rkzvector(ibz)

            call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
                  hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
                  hk_ev,e,vme,M1,T1)                             ! CHANGED: M1,T1 appended

            if (iflag_norder.eq.2) then
                  call get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
                  call get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der, &
                        shift_vector,berry_eigen1,berry_eigen2,berry_eigen, &
                        hk_ev_neigh,vme_neigh, &
                        skernel,hkernel,sderkernel,hderkernel,akernel,vme_der_phase, &
                        M1,T1)                                          ! CHANGED: M1,T1 appended
            end if

            do i=1,nband_ex
                  ii=nband_index(i)
                  ek(ibz,i)=e(ii)
                  do nj=1,3
                        do j=1,nband_ex
                              jj=nband_index(j)
                              vme_ex_band(ibz,nj,i,j)=vme(ii,jj,nj)
                              if (iflag_norder.eq.2) then
                                    shift_vector_ex_band(ibz,nj,1,i,j)=shift_vector(ii,jj,nj,1)
                                    shift_vector_ex_band(ibz,nj,2,i,j)=shift_vector(ii,jj,nj,2)
                                    shift_vector_ex_band(ibz,nj,3,i,j)=shift_vector(ii,jj,nj,3)
                                    gen_der_ex_band(ibz,nj,1,i,j)=gen_der(ii,jj,nj,1)
                                    gen_der_ex_band(ibz,nj,2,i,j)=gen_der(ii,jj,nj,2)
                                    gen_der_ex_band(ibz,nj,3,i,j)=gen_der(ii,jj,nj,3)
                                    berry_eigen_ex_band(ibz,nj,i,j)=berry_eigen(ii,jj,nj)
                              end if
                        end do
                  end do
            end do
      end do
      !$OMP END DO

      deallocate(skernel,hkernel,sderkernel,hderkernel,akernel)
      deallocate(gen_der,vme_der,gd1,gd2,gd3,hk_ev_neigh,vme_neigh,vme_der_phase)
      deallocate(M1,T1)                                          ! NEW
      !$OMP END PARALLEL

      write(*,*) '   Writing optical matrix elements (sp) into file'
      if (iflag_norder.eq.1) then
         call write_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
      end if
      if (iflag_norder.eq.2) then
         call write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek, &
            gen_der_ex_band,shift_vector_ex_band,berry_eigen_ex_band)
      end if
      write(*,*) '   Optical matrix elements (sp) have been written in file'

      deallocate(vme_ex_band,ek,berry_eigen_ex_band,gen_der_ex_band,shift_vector_ex_band)
end subroutine get_ome_sp
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der, &
      shift_vector,berry_eigen1,berry_eigen2,berry_eigen, &
      hk_ev_neigh,vme_neigh, &
      skernel,hkernel,sderkernel,hderkernel,akernel,vme_der_phase,&
      M1,T1)               
      implicit none

      integer norb
      integer nn,nnp
      integer ialpha,ialphap
      integer nj,njp

      dimension hkernel(norb,norb),skernel(norb,norb)
      dimension sderkernel(norb,norb,3),hderkernel(norb,norb,3)
      dimension akernel(norb,norb,3)

      dimension e(norb)

      dimension berry_eigen1(norb,norb,3),berry_eigen2(norb,norb,3)
      dimension berry_eigen(norb,norb,3)
      dimension vme_der(norb,norb,3,3)
      dimension vme_der_phase(norb,norb,3,3)

      complex*16 :: hk_ev_neigh(norb,norb,7)
      complex*16 :: vme_neigh(norb,norb,3,7)
      complex*16 :: M1(norb,norb), T1(norb,norb)                 ! NEW dummy args

      dimension shift_vector(norb,norb,3,3)

      real*8 rkx,rky,rkz,rkx_neigh,rky_neigh,rkz_neigh
      real*8 e
      real*8 vme_der_phase
      real*8 shift_vector
      real*8 ph1,ph2,ph3,ph4,ph5,ph6

      complex*16 hkernel,akernel,skernel,sderkernel,hderkernel
      complex*16 aux1,aux2,aux3,aux4,aux5,aux6
      complex*16 vme_der
      complex*16 berry_eigen1,berry_eigen2,berry_eigen
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! active_x/y/z are module-level, computed once in set_active_flags().

      vme_der=0.0d0
      shift_vector=0.0d0
      vme_der_phase=0.0d0
      hk_ev_neigh=0.0d0
      vme_neigh=0.0d0
      berry_eigen1=0.0d0
      berry_eigen2=0.0d0
      berry_eigen=0.0d0
      
      ! --- central point (index 7) — MOVED to the top of the routine ---
      call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(:,:,7),e,vme_neigh(:,:,:,7),M1,T1)

      ! --- x-direction neighbours ---
      if (active_x) then
            rkx_neigh=rkx-dk; rky_neigh=rky; rkz_neigh=rkz
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,1),e,vme_neigh(:,:,:,1),M1,T1)

            rkx_neigh=rkx+dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,3),e,vme_neigh(:,:,:,3),M1,T1)
      else
            hk_ev_neigh(:,:,1)   = hk_ev_neigh(:,:,7)     ! reuse central point, NO recomputation
            vme_neigh(:,:,:,1)   = vme_neigh(:,:,:,7)
            hk_ev_neigh(:,:,3)   = hk_ev_neigh(:,:,7)
            vme_neigh(:,:,:,3)   = vme_neigh(:,:,:,7)
      end if

      ! --- y-direction neighbours ---
      if (active_y) then
            rkx_neigh=rkx; rky_neigh=rky-dk; rkz_neigh=rkz
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,2),e,vme_neigh(:,:,:,2),M1,T1)

            rky_neigh=rky+dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,4),e,vme_neigh(:,:,:,4),M1,T1)
      else
            hk_ev_neigh(:,:,2)   = hk_ev_neigh(:,:,7)     ! reuse central point, NO recomputation
            vme_neigh(:,:,:,2)   = vme_neigh(:,:,:,7)
            hk_ev_neigh(:,:,4)   = hk_ev_neigh(:,:,7)
            vme_neigh(:,:,:,4)   = vme_neigh(:,:,:,7)
      end if
      
      
      ! --- z-direction neighbours ---
         if (active_z) then
            rkx_neigh=rkx; rky_neigh=rky; rkz_neigh=rkz-dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
            call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,5),e,vme_neigh(:,:,:,5),M1,T1)

            rkz_neigh=rkz+dk
            call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel,sderkernel,hkernel,hderkernel,akernel)
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel,hk_ev_neigh(:,:,6),e,vme_neigh(:,:,:,6),M1,T1)
      else
            hk_ev_neigh(:,:,5)   = hk_ev_neigh(:,:,7)     ! reuse central point, NO recomputation
            vme_neigh(:,:,:,5)   = vme_neigh(:,:,:,7)
            hk_ev_neigh(:,:,6)   = hk_ev_neigh(:,:,7)
            vme_neigh(:,:,:,6)   = vme_neigh(:,:,:,7)
      end if

      do nn=1,norb
         do nnp=1,norb
            do ialpha=1,norb
               do ialphap=1,norb
                  !x-dir
                  if (active_x) then
                  aux1=(hk_ev_neigh(ialphap,nnp,3)-hk_ev_neigh(ialphap,nnp,1))/(2.0d0*dk)
                  berry_eigen1(nn,nnp,1)=berry_eigen1(nn,nnp,1)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(ialpha,nn,7))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(nn,nnp,1)=berry_eigen2(nn,nnp,1)+ &
                  conjg(hk_ev_neigh(ialpha,nn,7))*hk_ev_neigh(ialphap,nnp,7)*akernel(ialpha,ialphap,1)

                  !y-dir
                  if (active_y) then
                  aux1=(hk_ev_neigh(ialphap,nnp,4)-hk_ev_neigh(ialphap,nnp,2))/(2.0d0*dk)
                  berry_eigen1(nn,nnp,2)=berry_eigen1(nn,nnp,2)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(ialpha,nn,7))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(nn,nnp,2)=berry_eigen2(nn,nnp,2)+ &
                  conjg(hk_ev_neigh(ialpha,nn,7))*hk_ev_neigh(ialphap,nnp,7)*akernel(ialpha,ialphap,2)

                  !z-dir
                  if (active_z) then
                  aux1=(hk_ev_neigh(ialphap,nnp,6)-hk_ev_neigh(ialphap,nnp,5))/(2.0d0*dk)
                  berry_eigen1(nn,nnp,3)=berry_eigen1(nn,nnp,3)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(ialpha,nn,7))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(nn,nnp,3)=berry_eigen2(nn,nnp,3)+ &
                  conjg(hk_ev_neigh(ialpha,nn,7))*hk_ev_neigh(ialphap,nnp,7)*akernel(ialpha,ialphap,3)
               end do
            end do

            do nj=1,3
                  berry_eigen(nn,nnp,nj)=berry_eigen1(nn,nnp,nj)+berry_eigen2(nn,nnp,nj)
                  if (abs(berry_eigen(nn,nnp,nj)).gt.clip_threshold) then
                        berry_eigen(nn,nnp,nj)=0.0d0
                  end if

                  aux1=vme_neigh(nn,nnp,nj,1); aux3=vme_neigh(nn,nnp,nj,3)
                  aux2=vme_neigh(nn,nnp,nj,2); aux4=vme_neigh(nn,nnp,nj,4)
                  aux5=vme_neigh(nn,nnp,nj,5); aux6=vme_neigh(nn,nnp,nj,6)

                  if (active_x) then
                        vme_der(nn,nnp,1,nj)=(aux3-aux1)/(2.0d0*dk)
                  end if
                  if (active_y) then
                        vme_der(nn,nnp,2,nj)=(aux4-aux2)/(2.0d0*dk)
                  end if
                  if (active_z) then
                        vme_der(nn,nnp,3,nj)=(aux6-aux5)/(2.0d0*dk)
                  end if

                  call get_phase(vme_neigh(nn,nnp,nj,1),ph1)
                  call get_phase(vme_neigh(nn,nnp,nj,3),ph3)
                  call get_phase(vme_neigh(nn,nnp,nj,2),ph2)
                  call get_phase(vme_neigh(nn,nnp,nj,4),ph4)
                  call get_phase(vme_neigh(nn,nnp,nj,5),ph5)
                  call get_phase(vme_neigh(nn,nnp,nj,6),ph6)

                  if (active_x) then
                        vme_der_phase(nn,nnp,1,nj)=(ph3-ph1)/(2.0d0*dk)
                  end if
                  if (active_y) then
                        vme_der_phase(nn,nnp,2,nj)=(ph4-ph2)/(2.0d0*dk)
                  end if
                  if (active_z) then
                        vme_der_phase(nn,nnp,3,nj)=(ph6-ph5)/(2.0d0*dk)
                  end if
            end do
         end do
      end do

      !shift vector
      do nn=1,norb
         do nnp=1,norb
            do nj=1,3
               do njp=1,3
                  shift_vector(nn,nnp,nj,njp)=-vme_der_phase(nn,nnp,nj,njp) &
                        +(realpart(berry_eigen(nn,nn,nj))-realpart(berry_eigen(nnp,nnp,nj)))
                  if (abs(shift_vector(nn,nnp,nj,njp)).gt.clip_threshold) then
                        shift_vector(nn,nnp,nj,njp)=0.0d0
                  end if
                  vme_der(nn,nnp,nj,njp)=vme_der(nn,nnp,nj,njp) &
                        -complex(0.0d0,1.0d0)*vme_neigh(nn,nnp,njp,7) &
                        *(realpart(berry_eigen(nn,nn,nj))-realpart(berry_eigen(nnp,nnp,nj)))
               end do
            end do
         end do
      end do

end subroutine get_berry_eigen_fourpoint
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_phase(aux1,ph)
      real*8 ph
      complex*16 aux1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (abs(aux1).lt.1.0d-05) then
         ph=0.0d0
      else
         ph=aimag(log(aux1))
      end if
   end subroutine get_phase
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! PATCH: gd1, gd2, gd3 are now dummy arguments again (caller-owned,
   ! allocated once per thread), NOT locals declared inside this
   ! subroutine -- see get_ome_sp for the reasoning.
   subroutine get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
      implicit none

      integer norb,norb_inter_cut
      integer nn,nnp,nnpp
      integer nj,njp

      dimension e(norb)
      dimension vme(norb,norb,3)
      dimension gen_der(norb,norb,3,3)
      dimension gd1(norb,norb,3,3)
      dimension gd2(norb,norb,3,3)
      dimension gd3(norb,norb,3,3)
      dimension abc(norb,norb,3)

      real*8 e
      complex*16 vme,gen_der,abc,gd1,gd2,gd3
      
      real*8 :: de_nnnp, inv_de_nnnp
      complex*16 :: ci_over_de
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      abc=0.0d0
      gd1=0.0d0
      gd2=0.0d0
      gd3=0.0d0
      gen_der=0.0d0
      norb_inter_cut=norb
      do nn=1,norb
         do nnp=1,norb
            do nj=1,3
               if (abs(e(nn)-e(nnp)).lt.1.0d-05) then
                  abc(nn,nnp,nj)=0.0d0
               else
                  abc(nn,nnp,nj)=-complex(0.0d0,1.0d0)*vme(nn,nnp,nj)/(e(nn)-e(nnp))
               end if
               do njp=1,3
                  de_nnnp = e(nn) - e(nnp)
                  if (abs(de_nnnp).lt.1.0d-05) then
                     gd1(nn,nnp,nj,njp)=0.0d0
                     gd2(nn,nnp,nj,njp)=0.0d0
                     gd3(nn,nnp,nj,njp)=0.0d0
                  else
                     inv_de_nnnp = 1.0d0 / de_nnnp
                     ci_over_de   = complex(0.0d0,1.0d0) * inv_de_nnnp
                     gd1(nn,nnp,nj,njp)=complex(0.0d0,1.0d0)*(inv_de_nnnp**2)* &
                        (vme(nn,nnp,nj)*(vme(nn,nn,njp)-vme(nnp,nnp,njp))+ &
                        vme(nn,nnp,njp)*(vme(nn,nn,nj)-vme(nnp,nnp,nj)))

                     gd2(nn,nnp,nj,njp)=0.0d0
                     gd3(nn,nnp,nj,njp)=0.0d0
                     do nnpp=1,norb_inter_cut
                        if (abs(e(nnpp)-e(nn)).lt.1.0d-05 .or. abs(e(nnpp)-e(nnp)).lt.1.0d-05) then
                              cycle
                        else
                              gd2(nn,nnp,nj,njp)=gd2(nn,nnp,nj,njp)+ &
                              ci_over_de*(vme(nn,nnpp,nj)*vme(nnpp,nnp,njp)/(e(nnpp)-e(nnp)))
                              gd3(nn,nnp,nj,njp)=gd3(nn,nnp,nj,njp)+ &
                              ci_over_de*(-vme(nn,nnpp,njp)*vme(nnpp,nnp,nj)/(e(nn)-e(nnpp)))
                        end if
                     end do
                  end if
                  gen_der(nn,nnp,nj,njp)=gd1(nn,nnp,nj,njp)+gd2(nn,nnp,nj,njp) &
                     +gd3(nn,nnp,nj,njp)
               end do
            end do
         end do
      end do

   end subroutine get_gen_der_sumrule
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
      implicit none
      integer iflag_norder
      integer npointstotal,nband_ex
      integer ibz
      integer nj,i,j

      dimension ek(npointstotal,nband_ex)
      dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)

      real*8 ek
      complex*16 vme_ex_band
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(10,file='ome_linear_sp_'//trim(material_name)//'.omesp')
      write(10,*) iflag_norder
      do ibz=1,npointstotal
         write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(ek(ibz,j),j=1,nband_ex)
         do i=1,nband_ex
            do j=1,nband_ex
               write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz), &
                  (realpart(vme_ex_band(ibz,nj,i,j)),aimag(vme_ex_band(ibz,nj,i,j)), nj=1,3)
            end do
         end do
      end do
      close(10)
   end subroutine write_ome_sp_linear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek, &
      gen_der_ex_band,shift_vector_ex_band,berry_eigen_ex_band)
      implicit none
      integer iflag_norder,npointstotal,nband_ex,ibz

      dimension ek(npointstotal,nband_ex)
      dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
      dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
      dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
      dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)

      real*8 ek, shift_vector_ex_band
      complex*16 vme_ex_band, berry_eigen_ex_band, gen_der_ex_band
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(10, file='ome_nonlinear_sp_'//trim(material_name)//'.omesp', &
           form='unformatted', access='stream', status='replace')

      write(10) iflag_norder
      write(10) npointstotal, nband_ex
      write(10) rkxvector, rkyvector, rkzvector

      do ibz=1,npointstotal
         write(10) ek(ibz,:)
         write(10) vme_ex_band(ibz,:,:,:)
         write(10) berry_eigen_ex_band(ibz,:,:,:)
         write(10) shift_vector_ex_band(ibz,:,:,:,:)
         write(10) gen_der_ex_band(ibz,:,:,:,:)
      end do

      close(10)
   end subroutine write_ome_sp_nonlinear
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
      hkernel,hderkernel,akernel)
      implicit none

      integer norb
      integer ialpha
      integer ialphap
      integer iRp
      integer nj

      dimension skernel(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(norb,norb,3)
      dimension hderkernel(norb,norb,3)
      dimension akernel(norb,norb,3)

      real(8) Rx,Ry,Rz
      real(8) rkx,rky,rkz

      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 phase,factor
      complex*16 hderhop_x,hderhop_y,hderhop_z
      complex*16 sderhop_x,sderhop_y,sderhop_z
      
      !real(8)    :: Rx_arr(nR), Ry_arr(nR), Rz_arr(nR)
      complex*16 :: factor_arr(nR)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! active_x/y/z now module-level (see set_active_flags),
      ! not recomputed on every call.

      hkernel=0.0d0
      hderkernel=0.0d0
      skernel=0.0d0
      sderkernel=0.0d0
      akernel=0.0d0

      
      do iRp = 1, nR
            factor_arr(iRp) = exp(complex(0.0d0, rkx*Rx_global(iRp) + rky*Ry_global(iRp) + rkz*Rz_global(iRp)))
      end do

      
      do ialpha=1,norb
         do ialphap=1,ialpha
            do iRp=1,nR
               

                  Rx = Rx_global(iRp); 
                  Ry = Ry_global(iRp); 
                  Rz = Rz_global(iRp)   ! still needed below for hderhop_x etc.
                  factor = factor_arr(iRp)
               
                  hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+factor*hhop(iRp,ialpha,ialphap)


                  skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
                  factor*shop(iRp,ialpha,ialphap)

                  hderhop_x=complex(0.0d0,Rx)*hhop(iRp,ialpha,ialphap)
                  hderhop_y=complex(0.0d0,Ry)*hhop(iRp,ialpha,ialphap)
                  hderhop_z=complex(0.0d0,Rz)*hhop(iRp,ialpha,ialphap)

                  sderhop_x=complex(0.0d0,Rx)*shop(iRp,ialpha,ialphap)
                  sderhop_y=complex(0.0d0,Ry)*shop(iRp,ialpha,ialphap)
                  sderhop_z=complex(0.0d0,Rz)*shop(iRp,ialpha,ialphap)

                  sderkernel(ialpha,ialphap,1)=sderkernel(ialpha,ialphap,1)+factor*sderhop_x
                  sderkernel(ialpha,ialphap,2)=sderkernel(ialpha,ialphap,2)+factor*sderhop_y
                  sderkernel(ialpha,ialphap,3)=sderkernel(ialpha,ialphap,3)+factor*sderhop_z

                  hderkernel(ialpha,ialphap,1)=hderkernel(ialpha,ialphap,1)+factor*hderhop_x
                  hderkernel(ialpha,ialphap,2)=hderkernel(ialpha,ialphap,2)+factor*hderhop_y
                  hderkernel(ialpha,ialphap,3)=hderkernel(ialpha,ialphap,3)+factor*hderhop_z

                  akernel(ialpha,ialphap,1)=akernel(ialpha,ialphap,1)+ &
                        factor*(rhop_c(1,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_x)
                  akernel(ialpha,ialphap,2)=akernel(ialpha,ialphap,2)+ &
                        factor*(rhop_c(2,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_y)
                  akernel(ialpha,ialphap,3)=akernel(ialpha,ialphap,3)+ &
                        factor*(rhop_c(3,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_z)
            end do

            if (ialphap < ialpha) then
               do nj=1,3
                  hkernel(ialphap,ialpha)=conjg(hkernel(ialpha,ialphap))
                  skernel(ialphap,ialpha)=conjg(skernel(ialpha,ialphap))
                  
                  sderkernel(ialphap,ialpha,nj)=conjg(sderkernel(ialpha,ialphap,nj))
                  hderkernel(ialphap,ialpha,nj)=conjg(hderkernel(ialpha,ialphap,nj))
                  akernel(ialphap,ialpha,nj)=conjg(akernel(ialpha,ialphap,nj))+ &
                     complex(0.0d0,1.0d0)*conjg(sderkernel(ialpha,ialphap,nj))
               end do
            end if

         end do
      end do
   end subroutine get_vme_kernels_ome
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    subroutine get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
!       hk_ev,e,vme)
!       implicit none
! 
!       integer :: norb
!       integer :: ialpha
!       integer :: ialphap
!       integer :: nj
!       integer :: nn,nnp
! 
!       dimension skernel(norb,norb)
!       dimension hkernel(norb,norb)
!       dimension sderkernel(3,norb,norb)
!       dimension hderkernel(3,norb,norb)
!       dimension akernel(3,norb,norb)
! 
!       dimension vjseudoa(3,norb,norb)
!       dimension vjseudob(3,norb,norb)
!       dimension e(norb)
!       dimension hk_ev(norb,norb)
!       dimension vme(3,norb,norb)
! 
!       real*8 e
!       complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
!       complex*16 hk_ev,vjseudoa,vjseudob,vme
!       complex*16 amu,amup
! 
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       e=0.0d0
!       call diagoz(norb,e,hkernel)
!       hk_ev(:,:)=hkernel(:,:)
!       call phase_eigvec_nk(norb,hk_ev)
!       vme=0.0d0
!       vjseudoa=0.0d0
!       vjseudob=0.0d0
      
      !vjseudoa(nj,:,:) = U^H · hderkernel(nj,:,:) · U
      !vjseudob(nj,nn,nnp) = i·[ e(nn)·T1(nj,nn,nnp) − e(nnp)·conjg(T1(nj,nnp,nn)) ]
!       do nn=1,norb
!          do nnp=1,nn
!             do ialpha=1,norb
!                do ialphap=1,norb
!                   amu=hk_ev(ialpha,nn)
!                   amup=hk_ev(ialphap,nnp)
!                   do nj=1,3
!                      vjseudoa(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+ &
!                         conjg(amu)*amup*hderkernel(nj,ialpha,ialphap)
!                      vjseudob(nj,nn,nnp)=vjseudob(nj,nn,nnp)+conjg(amu)*amup* &
!                         (e(nn)*akernel(nj,ialpha,ialphap)-e(nnp)*conjg(akernel(nj,ialphap,ialpha)))* &
!                         complex(0.0d0,1.0d0)
!                   end do
!                end do
!             end do
!             do nj=1,3
!                vme(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+vjseudob(nj,nn,nnp)
!                if (nnp < nn) then
!                   vme(nj,nnp,nn)=conjg(vme(nj,nn,nnp))
!                   vjseudoa(nj,nnp,nn)=conjg(vjseudoa(nj,nn,nnp))
!                   vjseudob(nj,nnp,nn)=conjg(vjseudob(nj,nn,nnp))
!                end if
!             end do
!          end do
!       end do
!    end subroutine get_vme_eigen_ome

subroutine get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
      hk_ev,e,vme,M1,T1)
      implicit none

      integer :: norb, nj, nn, nnp

      dimension skernel(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(norb,norb,3)      ! <-- reordered
      dimension hderkernel(norb,norb,3)      ! <-- reordered
      dimension akernel(norb,norb,3)         ! <-- reordered
      dimension e(norb)
      dimension hk_ev(norb,norb)
      dimension vme(norb,norb,3)             ! <-- reordered

      real*8 e
      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 hk_ev,vme
      complex*16 M1(norb,norb), T1(norb,norb)
      complex*16, parameter :: cone=(1.0d0,0.0d0), czero=(0.0d0,0.0d0), ci=(0.0d0,1.0d0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      e=0.0d0
      call diagoz(norb,e,hkernel)
      hk_ev(:,:)=hkernel(:,:)
      call phase_eigvec_nk(norb,hk_ev)
      vme=0.0d0

      do nj=1,3
         ! Both slices below are now genuinely contiguous (norb*norb
         ! contiguous elements each) — no compiler-inserted temporary,
         ! at compile time OR at runtime under -fcheck=array-temps.
         call zgemm('N','N',norb,norb,norb, cone, hderkernel(:,:,nj), norb, hk_ev, norb, czero, M1, norb)
         call zgemm('C','N',norb,norb,norb, cone, hk_ev, norb, M1, norb, czero, vme(:,:,nj), norb)

         call zgemm('N','N',norb,norb,norb, cone, akernel(:,:,nj), norb, hk_ev, norb, czero, M1, norb)
         call zgemm('C','N',norb,norb,norb, cone, hk_ev, norb, M1, norb, czero, T1, norb)

         do nnp=1,norb
            do nn=1,norb
               vme(nn,nnp,nj) = vme(nn,nnp,nj) + ci*( e(nn)*T1(nn,nnp) - e(nnp)*conjg(T1(nnp,nn)) )
            end do
         end do
      end do

   end subroutine get_vme_eigen_ome
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine phase_eigvec_nk(norb,hk_ev)
      implicit none

      integer norb
      integer i,j,ii

      dimension hk_ev(norb,norb)

      real*8 :: arg
      complex*16 :: aux1,hk_ev,factor
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,norb
         aux1=0.0d0
         do i=1,norb
            aux1=aux1+hk_ev(i,j)
         end do
         arg=atan2(aimag(aux1),realpart(aux1))
         factor=exp(complex(0.0d0,-arg))
         do ii=1,norb
            hk_ev(ii,j)=hk_ev(ii,j)*factor
         end do
      end do

   end subroutine phase_eigvec_nk
end module ome_sp