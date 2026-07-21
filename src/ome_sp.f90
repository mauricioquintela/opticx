module ome_sp
   use constants_math
   use parser_wannier90_tb, &
      only:material_name,nR,nRvec,norb,R,shop,hhop,rhop_c
   use parser_optics_xatu_dim, &
      only:npointstotal,rkxvector,rkyvector,rkzvector, &
      nband_ex,nband_index,nv_ex,nc_ex

   implicit none

   allocatable vme_ex_band(:,:,:,:)
   allocatable ek(:,:)
   allocatable gen_der_ex_band(:,:,:,:,:)
   allocatable shift_vector_ex_band(:,:,:,:,:)
   allocatable berry_eigen_ex_band(:,:,:,:)

   real*8 ek
   real*8 shift_vector_ex_band
   complex*16 vme_ex_band
   complex*16 gen_der_ex_band
   complex*16 berry_eigen_ex_band

   ! PATCH: dimensionality flags computed ONCE (via set_active_flags) instead
   ! of being recomputed inside every call to get_vme_kernels_ome /
   ! get_berry_eigen_fourpoint (was previously duplicated in two places and
   ! recomputed npointstotal*8 times). Read-only after initialization, so
   ! safe to share across OMP threads without being PRIVATE.
   logical, save :: active_x = .false., active_y = .false., active_z = .false.
   logical, save :: active_flags_set = .false.

   ! PATCH: named clipping threshold instead of a magic literal duplicated
   ! in two places (berry_eigen and shift_vector divergence clipping).
   real(8), parameter :: clip_threshold = 50.0d0

contains
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

      ! PATCH: these were `allocatable` + pre-loop `allocate`, then marked
      ! OMP PRIVATE. Per the OpenMP/Fortran spec that gives each thread an
      ! UNALLOCATED private copy at the start of the parallel region (the
      ! pre-loop allocate does not carry over) -> undefined behavior on any
      ! read/write inside the loop. Switching to automatic (non-allocatable)
      ! arrays sized by norb fixes this: each thread gets a correctly-sized
      ! private instance automatically, with no allocate/deallocate needed.
      complex*16 :: skernel(norb,norb), hkernel(norb,norb)
      complex*16 :: sderkernel(3,norb,norb), hderkernel(3,norb,norb)
      complex*16 :: akernel(3,norb,norb)

      real(8)    :: e(norb)
      complex*16 :: hk_ev(norb,norb)
      complex*16 :: vme(3,norb,norb)

      complex*16 :: abc(3,norb,norb)
      complex*16 :: gd1(3,3,norb,norb), gd2(3,3,norb,norb), gd3(3,3,norb,norb)
      complex*16 :: gen_der(3,3,norb,norb)
      complex*16 :: vme_der(3,3,norb,norb)
      real(8)    :: shift_vector(3,3,norb,norb)
      complex*16 :: berry_eigen1(3,norb,norb), berry_eigen2(3,norb,norb)
      complex*16 :: berry_eigen(3,norb,norb)

      real(8) :: rkx,rky,rkz
      ! PATCH: removed unused `ecomplex`, `naux1`, `ibz_sum` (declared/
      ! allocated in the original but never read or incremented).
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) '5. Entering ome_sp'

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

      ! PATCH: compute lattice-dimensionality flags exactly once, up front,
      ! instead of inside the hot k-point loop.
      call set_active_flags()

      write(*,*) '   Calculating optical matrix elements (sp): sampling BZ...'

      !$OMP PARALLEL DO PRIVATE(rkx,rky,rkz), &
      !$OMP PRIVATE(hkernel,skernel,sderkernel,hderkernel,akernel), &
      !$OMP PRIVATE(hk_ev,e,vme), &
      !$OMP PRIVATE(abc,gen_der,gd1,gd2,gd3), &
      !$OMP PRIVATE(vme_der,shift_vector,berry_eigen1,berry_eigen2,berry_eigen)
      do ibz=1,npointstotal
         write(*,*) '   Optical matrix elements (sp): k-point',ibz,'/',npointstotal
         rkx=rkxvector(ibz)
         rky=rkyvector(ibz)
         rkz=rkzvector(ibz)

         call get_vme_kernels_ome(rkx,rky,rkz,norb,skernel,sderkernel, &
            hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev,e,vme)

         if (iflag_norder.eq.2) then
            call get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
            call get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der, &
               shift_vector,berry_eigen1,berry_eigen2,berry_eigen)
         end if

         do i=1,nband_ex
            ii=nband_index(i)
            ek(ibz,i)=e(ii)
            do nj=1,3
               do j=1,nband_ex
                  jj=nband_index(j)
                  vme_ex_band(ibz,nj,i,j)=vme(nj,ii,jj)

                  if (iflag_norder.eq.2) then
                     shift_vector_ex_band(ibz,nj,1,i,j)=shift_vector(nj,1,ii,jj)
                     shift_vector_ex_band(ibz,nj,2,i,j)=shift_vector(nj,2,ii,jj)
                     shift_vector_ex_band(ibz,nj,3,i,j)=shift_vector(nj,3,ii,jj)

                     gen_der_ex_band(ibz,nj,1,i,j)=gen_der(nj,1,ii,jj)
                     gen_der_ex_band(ibz,nj,2,i,j)=gen_der(nj,2,ii,jj)
                     gen_der_ex_band(ibz,nj,3,i,j)=gen_der(nj,3,ii,jj)

                     berry_eigen_ex_band(ibz,nj,i,j)=berry_eigen(nj,ii,jj)
                  end if
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      write(*,*) '   Writing optical matrix elements (sp) into file'
      if (iflag_norder.eq.1) then
         call write_ome_sp_linear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek)
      end if
      if (iflag_norder.eq.2) then
         call write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek, &
            gen_der_ex_band,shift_vector_ex_band,berry_eigen_ex_band)
      end if
      write(*,*) '   Optical matrix elements (sp) have been written in file'
   end subroutine get_ome_sp
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_berry_eigen_fourpoint(rkx,rky,rkz,norb,vme_der, &
      shift_vector,berry_eigen1,berry_eigen2,berry_eigen)
      implicit none

      ! PATCH: removed local `integer nR,norb` shadowing the module-level
      ! `nR` from parser_wannier90_tb. Only `norb` is actually a dummy
      ! argument here; `nR` was an unused, silently-shadowing landmine.
      integer norb
      integer nn,nnp
      integer ialpha,ialphap
      integer nj,njp

      dimension hkernel(norb,norb),skernel(norb,norb)
      dimension sderkernel(3,norb,norb),hderkernel(3,norb,norb)

      ! PATCH: removed unused `hk_alpha`.
      dimension hk_ev(norb,norb),e(norb)
      dimension akernel(3,norb,norb)

      dimension vjseudoa(3,norb,norb),vjseudob(3,norb,norb),vme(3,norb,norb)
      dimension berry_eigen1(3,norb,norb),berry_eigen2(3,norb,norb)
      dimension berry_eigen(3,norb,norb)
      dimension vme_der(3,3,norb,norb)
      dimension vme_der_phase(3,3,norb,norb)
      ! 7 neighbour slots: 1/3 = x-,x+ ; 2/4 = y-,y+ ; 5/6 = z-,z+ ; 7 = central
      dimension hk_ev_neigh(7,norb,norb)
      dimension vme_neigh(7,3,norb,norb)

      dimension shift_vector(3,3,norb,norb)

      real*8 rkx,rky,rkz,rkx_neigh,rky_neigh,rkz_neigh
      real*8 e
      real*8 vme_der_phase
      real*8 shift_vector
      real*8 ph1,ph2,ph3,ph4,ph5,ph6

      complex*16 hkernel,akernel,skernel,sderkernel,hderkernel
      complex*16 hk_ev,vjseudoa,vjseudob,vme
      complex*16 aux1,aux2,aux3,aux4,aux5,aux6
      complex*16 vme_der
      complex*16 berry_eigen1,berry_eigen2,berry_eigen
      complex*16 hk_ev_neigh,vme_neigh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! PATCH: active_x/y/z are now module-level, computed once in
      ! set_active_flags() -- no longer recomputed here on every one of the
      ! ~8 calls per k-point.

      vme=0.0d0
      hk_ev_neigh=0.0d0
      vme_neigh=0.0d0
      berry_eigen1=0.0d0
      berry_eigen2=0.0d0
      berry_eigen=0.0d0
      vme_der=0.0d0
      shift_vector=0.0d0
      vme_der_phase=0.0d0

      ! --- x-direction neighbours ---
      if (active_x) then
         rkx_neigh=rkx-dk; rky_neigh=rky; rkz_neigh=rkz
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(1,:,:),e,vme_neigh(1,:,:,:))

         rkx_neigh=rkx+dk
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(3,:,:),e,vme_neigh(3,:,:,:))
      else
         ! PATCH: previously called get_vme_kernels_ome + get_vme_eigen_ome
         ! TWICE at the same unshifted point to fill both neighbour slots,
         ! even though both results are gated to zero by active_x anyway.
         ! Now: compute once and copy into both slots (or just leave them
         ! at the pre-zeroed default -- they're never used when active_x is
         ! false, since the difference below is skipped). Left as an
         ! explicit single call + copy for clarity / in case of future use.
         rkx_neigh=rkx; rky_neigh=rky; rkz_neigh=rkz
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(1,:,:),e,vme_neigh(1,:,:,:))
         hk_ev_neigh(3,:,:) = hk_ev_neigh(1,:,:)
         vme_neigh(3,:,:,:) = vme_neigh(1,:,:,:)
      end if

      ! --- y-direction neighbours ---
      if (active_y) then
         rkx_neigh=rkx; rky_neigh=rky-dk; rkz_neigh=rkz
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(2,:,:),e,vme_neigh(2,:,:,:))

         rky_neigh=rky+dk
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(4,:,:),e,vme_neigh(4,:,:,:))
      else
         rkx_neigh=rkx; rky_neigh=rky; rkz_neigh=rkz
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(2,:,:),e,vme_neigh(2,:,:,:))
         hk_ev_neigh(4,:,:) = hk_ev_neigh(2,:,:)
         vme_neigh(4,:,:,:) = vme_neigh(2,:,:,:)
      end if

      ! --- z-direction neighbours ---
      if (active_z) then
         rkx_neigh=rkx; rky_neigh=rky; rkz_neigh=rkz-dk
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(5,:,:),e,vme_neigh(5,:,:,:))

         rkz_neigh=rkz+dk
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(6,:,:),e,vme_neigh(6,:,:,:))
      else
         rkx_neigh=rkx; rky_neigh=rky; rkz_neigh=rkz
         call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
         call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(5,:,:),e,vme_neigh(5,:,:,:))
         hk_ev_neigh(6,:,:) = hk_ev_neigh(5,:,:)
         vme_neigh(6,:,:,:) = vme_neigh(5,:,:,:)
      end if

      ! --- central point (index 7) ---
      rkx_neigh=rkx; rky_neigh=rky; rkz_neigh=rkz
      call get_vme_kernels_ome(rkx_neigh,rky_neigh,rkz_neigh,norb,skernel, &
            sderkernel,hkernel,hderkernel,akernel)
      call get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
            hk_ev_neigh(7,:,:),e,vme_neigh(7,:,:,:))

      do nn=1,norb
         do nnp=1,norb
            do ialpha=1,norb
               do ialphap=1,norb
                  !x-dir
                  if (active_x) then
                     aux1=(hk_ev_neigh(3,ialphap,nnp)-hk_ev_neigh(1,ialphap,nnp))/(2.0d0*dk)
                     berry_eigen1(1,nn,nnp)=berry_eigen1(1,nn,nnp)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(7,ialpha,nn))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(1,nn,nnp)=berry_eigen2(1,nn,nnp)+ &
                     conjg(hk_ev_neigh(7,ialpha,nn))*hk_ev_neigh(7,ialphap,nnp)*akernel(1,ialpha,ialphap)

                  !y-dir
                  if (active_y) then
                     aux1=(hk_ev_neigh(4,ialphap,nnp)-hk_ev_neigh(2,ialphap,nnp))/(2.0d0*dk)
                     berry_eigen1(2,nn,nnp)=berry_eigen1(2,nn,nnp)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(7,ialpha,nn))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(2,nn,nnp)=berry_eigen2(2,nn,nnp)+ &
                     conjg(hk_ev_neigh(7,ialpha,nn))*hk_ev_neigh(7,ialphap,nnp)*akernel(2,ialpha,ialphap)

                  !z-dir
                  if (active_z) then
                     aux1=(hk_ev_neigh(6,ialphap,nnp)-hk_ev_neigh(5,ialphap,nnp))/(2.0d0*dk)
                     berry_eigen1(3,nn,nnp)=berry_eigen1(3,nn,nnp)+ &
                        complex(0.0d0,1.0d0)*conjg(hk_ev_neigh(7,ialpha,nn))*skernel(ialpha,ialphap)*aux1
                  end if
                  berry_eigen2(3,nn,nnp)=berry_eigen2(3,nn,nnp)+ &
                     conjg(hk_ev_neigh(7,ialpha,nn))*hk_ev_neigh(7,ialphap,nnp)*akernel(3,ialpha,ialphap)
               end do
            end do

            do nj=1,3
               berry_eigen(nj,nn,nnp)=berry_eigen1(nj,nn,nnp)+berry_eigen2(nj,nn,nnp)
               if (abs(berry_eigen(nj,nn,nnp)).gt.clip_threshold) then
                  berry_eigen(nj,nn,nnp)=0.0d0
               end if

               aux1=vme_neigh(1,nj,nn,nnp); aux3=vme_neigh(3,nj,nn,nnp)
               aux2=vme_neigh(2,nj,nn,nnp); aux4=vme_neigh(4,nj,nn,nnp)
               aux5=vme_neigh(5,nj,nn,nnp); aux6=vme_neigh(6,nj,nn,nnp)

               if (active_x) then
                  vme_der(1,nj,nn,nnp)=(aux3-aux1)/(2.0d0*dk)
               end if
               if (active_y) then
                  vme_der(2,nj,nn,nnp)=(aux4-aux2)/(2.0d0*dk)
               end if
               if (active_z) then
                  vme_der(3,nj,nn,nnp)=(aux6-aux5)/(2.0d0*dk)
               end if

               call get_phase(vme_neigh(1,nj,nn,nnp),ph1)
               call get_phase(vme_neigh(3,nj,nn,nnp),ph3)
               call get_phase(vme_neigh(2,nj,nn,nnp),ph2)
               call get_phase(vme_neigh(4,nj,nn,nnp),ph4)
               call get_phase(vme_neigh(5,nj,nn,nnp),ph5)
               call get_phase(vme_neigh(6,nj,nn,nnp),ph6)

               if (active_x) then
                  vme_der_phase(1,nj,nn,nnp)=(ph3-ph1)/(2.0d0*dk)
               end if
               if (active_y) then
                  vme_der_phase(2,nj,nn,nnp)=(ph4-ph2)/(2.0d0*dk)
               end if
               if (active_z) then
                  vme_der_phase(3,nj,nn,nnp)=(ph6-ph5)/(2.0d0*dk)
               end if
               ! Note: vme_der / vme_der_phase were initialised to 0.0d0
               ! above, so inactive directions are explicitly zero rather
               ! than left as uninitialised stack memory (see the old 2D
               ! version, where the z-slot of vme_der_phase was read
               ! without ever being assigned).
            end do
         end do
      end do

      !shift vector
      do nn=1,norb
         do nnp=1,norb
            do nj=1,3
               do njp=1,3
                  shift_vector(nj,njp,nn,nnp)=-vme_der_phase(nj,njp,nn,nnp) &
                     +(realpart(berry_eigen(nj,nn,nn))-realpart(berry_eigen(nj,nnp,nnp)))
                  if (abs(shift_vector(nj,njp,nn,nnp)).gt.clip_threshold) then
                     shift_vector(nj,njp,nn,nnp)=0.0d0
                  end if
                  vme_der(nj,njp,nn,nnp)=vme_der(nj,njp,nn,nnp) &
                     -complex(0.0d0,1.0d0)*vme_neigh(7,njp,nn,nnp) &
                     *(realpart(berry_eigen(nj,nn,nn))-realpart(berry_eigen(nj,nnp,nnp)))
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
   subroutine get_gen_der_sumrule(norb,vme,e,abc,gen_der,gd1,gd2,gd3)
      ! PATCH: was `implicit real*8 (a-h,o-z)`, i.e. any mistyped/new
      ! variable name silently became an implicitly-typed real or integer
      ! instead of a compile-time error. Switched to implicit none like
      ! every other routine in the module.
      implicit none

      integer norb,norb_inter_cut
      integer nn,nnp,nnpp
      integer nj,njp

      dimension e(norb)
      dimension vme(3,norb,norb)
      dimension gen_der(3,3,norb,norb)
      dimension gd1(3,3,norb,norb)
      dimension gd2(3,3,norb,norb)
      dimension gd3(3,3,norb,norb)
      dimension abc(3,norb,norb)

      real*8 e
      ! PATCH: removed unused aux1,aux2,aux3.
      complex*16 vme,gen_der,abc,gd1,gd2,gd3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
                  abc(nj,nn,nnp)=0.0d0
               else
                  abc(nj,nn,nnp)=-complex(0.0d0,1.0d0)*vme(nj,nn,nnp)/(e(nn)-e(nnp))
               end if
               do njp=1,3
                  if (abs(e(nn)-e(nnp)).lt.1.0d-05) then
                     gd1(nj,njp,nn,nnp)=0.0d0
                     gd2(nj,njp,nn,nnp)=0.0d0
                     gd3(nj,njp,nn,nnp)=0.0d0
                  else
                     gd1(nj,njp,nn,nnp)=complex(0.0d0,1.0d0)/((e(nn)-e(nnp))**2)* &
                        (vme(nj,nn,nnp)*(vme(njp,nn,nn)-vme(njp,nnp,nnp))+ &
                        vme(njp,nn,nnp)*(vme(nj,nn,nn)-vme(nj,nnp,nnp)))

                     gd2(nj,njp,nn,nnp)=0.0d0
                     gd3(nj,njp,nn,nnp)=0.0d0
                     do nnpp=1,norb_inter_cut
                        if (abs(e(nnpp)-e(nn)).lt.1.0d-05 .or. &
                           abs(e(nnpp)-e(nnp)).lt.1.0d-05) then
                           gd2(nj,njp,nn,nnp)=gd2(nj,njp,nn,nnp)+0.0d0
                           gd3(nj,njp,nn,nnp)=gd3(nj,njp,nn,nnp)+0.0d0
                        else
                           gd2(nj,njp,nn,nnp)=gd2(nj,njp,nn,nnp)+ &
                              complex(0.0d0,1.0d0)/(e(nn)-e(nnp))* &
                              (vme(nj,nn,nnpp)*vme(njp,nnpp,nnp)/(e(nnpp)-e(nnp)))

                           gd3(nj,njp,nn,nnp)=gd3(nj,njp,nn,nnp)+ &
                              complex(0.0d0,1.0d0)/(e(nn)-e(nnp))* &
                              (-vme(njp,nn,nnpp)*vme(nj,nnpp,nnp)/(e(nn)-e(nnpp)))
                        end if
                     end do
                  end if
                  gen_der(nj,njp,nn,nnp)=gd1(nj,njp,nn,nnp)+gd2(nj,njp,nn,nnp) &
                     +gd3(nj,njp,nn,nnp)
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
   ! PATCH: unformatted stream I/O, bulk array writes/reads instead of
   ! 8*nband_ex^2 scalar-by-scalar formatted writes/reads per k-point.
   ! See prior discussion for the matching read_ome_sp_nonlinear.
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
   
!       subroutine write_ome_sp_nonlinear(iflag_norder,npointstotal,nband_ex,vme_ex_band,ek,gen_der_ex_band, &
!       shift_vector_ex_band,berry_eigen_ex_band)
!       implicit none
!       integer iflag_norder
!       integer npointstotal,nband_ex
!       integer ibz
!       integer nj,i,j
!  
!       dimension ek(npointstotal,nband_ex)
!       dimension vme_ex_band(npointstotal,3,nband_ex,nband_ex)
!       dimension berry_eigen_ex_band(npointstotal,3,nband_ex,nband_ex)
!       dimension gen_der_ex_band(npointstotal,3,3,nband_ex,nband_ex)
!       dimension shift_vector_ex_band(npointstotal,3,3,nband_ex,nband_ex)
!  
!       real*8 ek
!       real*8 shift_vector_ex_band
!       complex*16 vme_ex_band
!       complex*16 berry_eigen_ex_band
!       complex*16 gen_der_ex_band
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       open(10,file='ome_nonlinear_sp_'//trim(material_name)//'.omesp')
!       write(10,*) iflag_norder
!       do ibz=1,npointstotal
!          write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(ek(ibz,j),j=1,nband_ex)
!          do i=1,nband_ex
!             do j=1,nband_ex
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz), &
!                   (realpart(vme_ex_band(ibz,nj,i,j)),aimag(vme_ex_band(ibz,nj,i,j)), nj=1,3)
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(berry_eigen_ex_band(ibz,nj,i,j)), &
!                   aimag(berry_eigen_ex_band(ibz,nj,i,j)), nj=1,3)
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(shift_vector_ex_band(ibz,1,nj,i,j), nj=1,3)
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(shift_vector_ex_band(ibz,2,nj,i,j), nj=1,3)
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(shift_vector_ex_band(ibz,3,nj,i,j), nj=1,3)
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(gen_der_ex_band(ibz,1,nj,i,j)), &
!                   aimag(gen_der_ex_band(ibz,1,nj,i,j)), nj=1,3)
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(gen_der_ex_band(ibz,2,nj,i,j)), &
!                   aimag(gen_der_ex_band(ibz,2,nj,i,j)), nj=1,3)
!                write(10,*) rkxvector(ibz),rkyvector(ibz),rkzvector(ibz),(realpart(gen_der_ex_band(ibz,3,nj,i,j)), &
!                   aimag(gen_der_ex_band(ibz,3,nj,i,j)), nj=1,3)
!             end do
!          end do
!       end do
!       close(10)
!    end subroutine write_ome_sp_nonlinear
   
   
   
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
      dimension sderkernel(3,norb,norb)
      dimension hderkernel(3,norb,norb)
      dimension akernel(3,norb,norb)

      real(8) Rx,Ry,Rz
      real(8) rkx,rky,rkz

      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 phase,factor
      ! PATCH: hderhop/sderhop were (3,nR,norb,norb) automatic (stack)
      ! arrays, but only the current-iRp slice was ever used before being
      ! overwritten -- a large unnecessary stack allocation, risky under
      ! OpenMP where each thread needs its own copy. Replaced with plain
      ! scalars, computed and consumed within the same iRp iteration.
      complex*16 hderhop_x,hderhop_y,hderhop_z
      complex*16 sderhop_x,sderhop_y,sderhop_z
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! PATCH: active_x/y/z now module-level (see set_active_flags),
      ! no longer recomputed on every call.

      hkernel=0.0d0
      hderkernel=0.0d0
      skernel=0.0d0
      sderkernel=0.0d0
      akernel=0.0d0

      do ialpha=1,norb
         do ialphap=1,ialpha
            do iRp=1,nR
               if (active_x) then
                  Rx = dble(nRvec(iRp,1))*R(1,1) + dble(nRvec(iRp,2))*R(2,1) + dble(nRvec(iRp,3))*R(3,1)
               else
                  Rx = 0.0d0
               end if
               if (active_y) then
                  Ry = dble(nRvec(iRp,1))*R(1,2) + dble(nRvec(iRp,2))*R(2,2) + dble(nRvec(iRp,3))*R(3,2)
               else
                  Ry = 0.0d0
               end if
               if (active_z) then
                  Rz = dble(nRvec(iRp,1))*R(1,3) + dble(nRvec(iRp,2))*R(2,3) + dble(nRvec(iRp,3))*R(3,3)
               else
                  Rz = 0.0d0
               end if

               phase=complex(0.0d0,rkx*Rx+rky*Ry+rkz*Rz)
               factor=exp(phase)

               hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+ &
                  factor*hhop(iRp,ialpha,ialphap)
               skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
                  factor*shop(iRp,ialpha,ialphap)

               hderhop_x=complex(0.0d0,Rx)*hhop(iRp,ialpha,ialphap)
               hderhop_y=complex(0.0d0,Ry)*hhop(iRp,ialpha,ialphap)
               hderhop_z=complex(0.0d0,Rz)*hhop(iRp,ialpha,ialphap)

               sderhop_x=complex(0.0d0,Rx)*shop(iRp,ialpha,ialphap)
               sderhop_y=complex(0.0d0,Ry)*shop(iRp,ialpha,ialphap)
               sderhop_z=complex(0.0d0,Rz)*shop(iRp,ialpha,ialphap)

               sderkernel(1,ialpha,ialphap)=sderkernel(1,ialpha,ialphap)+factor*sderhop_x
               sderkernel(2,ialpha,ialphap)=sderkernel(2,ialpha,ialphap)+factor*sderhop_y
               sderkernel(3,ialpha,ialphap)=sderkernel(3,ialpha,ialphap)+factor*sderhop_z

               hderkernel(1,ialpha,ialphap)=hderkernel(1,ialpha,ialphap)+factor*hderhop_x
               hderkernel(2,ialpha,ialphap)=hderkernel(2,ialpha,ialphap)+factor*hderhop_y
               hderkernel(3,ialpha,ialphap)=hderkernel(3,ialpha,ialphap)+factor*hderhop_z

               akernel(1,ialpha,ialphap)=akernel(1,ialpha,ialphap)+ &
                  factor*(rhop_c(1,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_x)
               akernel(2,ialpha,ialphap)=akernel(2,ialpha,ialphap)+ &
                  factor*(rhop_c(2,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_y)
               akernel(3,ialpha,ialphap)=akernel(3,ialpha,ialphap)+ &
                  factor*(rhop_c(3,iRp,ialpha,ialphap)+complex(0.0d0,1.0d0)*sderhop_z)
            end do

            ! PATCH: the original mirrored (ialpha,ialphap) -> (ialphap,ialpha)
            ! for ialphap=1,ialpha, which INCLUDES the diagonal ialphap==ialpha.
            ! On the diagonal that overwrites e.g. sderkernel(nj,ialpha,ialpha)
            ! with its own conjugate before akernel's diagonal term reads it,
            ! making the diagonal formula inconsistent with the off-diagonal
            ! one and flipping any floating-point imaginary residue. Only
            ! mirror strict off-diagonal pairs; the diagonal element was
            ! already fully computed by the accumulation loop above.
            if (ialphap < ialpha) then
               do nj=1,3
                  hkernel(ialphap,ialpha)=conjg(hkernel(ialpha,ialphap))
                  skernel(ialphap,ialpha)=conjg(skernel(ialpha,ialphap))
                  sderkernel(nj,ialphap,ialpha)=conjg(sderkernel(nj,ialpha,ialphap))
                  hderkernel(nj,ialphap,ialpha)=conjg(hderkernel(nj,ialpha,ialphap))
                  akernel(nj,ialphap,ialpha)=conjg(akernel(nj,ialpha,ialphap))+ &
                     complex(0.0d0,1.0d0)*conjg(sderkernel(nj,ialpha,ialphap))
               end do
            end if

         end do
      end do
   end subroutine get_vme_kernels_ome
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! NOTE (unchanged, flagged for review): skernel/sderkernel are received
   ! as dummy arguments but never referenced below -- diagonalization is a
   ! standard Hermitian eigenproblem (`diagoz(norb,e,hkernel)`), not the
   ! generalized form (H c = E S c) that a non-orthogonal Wannier basis
   ! would require. If shop/skernel is not effectively the identity, this
   ! is a physics bug, not a style issue -- please confirm the basis is
   ! orthogonal before relying on this routine.
   subroutine get_vme_eigen_ome(norb,skernel,sderkernel,hkernel,hderkernel,akernel, &
      hk_ev,e,vme)
      implicit none

      integer :: norb
      integer :: ialpha
      integer :: ialphap
      integer :: nj
      integer :: nn,nnp

      dimension skernel(norb,norb)
      dimension hkernel(norb,norb)
      dimension sderkernel(3,norb,norb)
      dimension hderkernel(3,norb,norb)
      dimension akernel(3,norb,norb)

      dimension vjseudoa(3,norb,norb)
      dimension vjseudob(3,norb,norb)
      dimension e(norb)
      dimension hk_ev(norb,norb)
      dimension vme(3,norb,norb)

      real*8 e
      complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
      complex*16 hk_ev,vjseudoa,vjseudob,vme
      complex*16 amu,amup

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      e=0.0d0
      call diagoz(norb,e,hkernel)
      hk_ev(:,:)=hkernel(:,:)
      call phase_eigvec_nk(norb,hk_ev)
      vme=0.0d0
      vjseudoa=0.0d0
      vjseudob=0.0d0

      do nn=1,norb
         do nnp=1,nn
            do ialpha=1,norb
               do ialphap=1,norb
                  amu=hk_ev(ialpha,nn)
                  amup=hk_ev(ialphap,nnp)
                  do nj=1,3
                     vjseudoa(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+ &
                        conjg(amu)*amup*hderkernel(nj,ialpha,ialphap)
                     vjseudob(nj,nn,nnp)=vjseudob(nj,nn,nnp)+conjg(amu)*amup* &
                        (e(nn)*akernel(nj,ialpha,ialphap)-e(nnp)*conjg(akernel(nj,ialphap,ialpha)))* &
                        complex(0.0d0,1.0d0)
                  end do
               end do
            end do
            ! PATCH: same diagonal double-conjugation issue as in
            ! get_vme_kernels_ome -- `do nnp=1,nn` includes nnp==nn, so the
            ! unconditional mirror below was applying conjg() to the
            ! diagonal element too. Diagonal terms are computed directly
            ! by the accumulation above and need no mirroring.
            do nj=1,3
               vme(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+vjseudob(nj,nn,nnp)
               if (nnp < nn) then
                  vme(nj,nnp,nn)=conjg(vme(nj,nn,nnp))
                  vjseudoa(nj,nnp,nn)=conjg(vjseudoa(nj,nn,nnp))
                  vjseudob(nj,nnp,nn)=conjg(vjseudob(nj,nn,nnp))
               end if
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