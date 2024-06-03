!#define OLD_SCHEME
#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_fxy_disp
#ifdef MPI
use mpi
#endif
use mod_grid
! === Dispersive ===============================================================
use mod_params, only : max_step, min_depth
#ifdef CONV_CHECK
use mod_params, only : conv_val
#endif
! ==============================================================================
#ifdef MPI
use mod_mpi
#endif
use mod_nest
use mod_params, only : VEL
! === Limiter with max Froude number. ==========================================
use mod_params, only : froude_lim
! ==============================================================================
#ifdef BANKFILE
use mod_params, only : broken_rate
#endif
#ifdef NORMALMODE
use mod_normalmode, only : calc_nm_P
use mod_params, only : m_rho
#endif
implicit none

contains

#ifndef MPI
#ifndef CONV_CHECK
#ifndef NORMALMODE
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,nlon,nlat,gflag,fg,cg)
#else
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,nlon,nlat,gflag,fg,cg,istep)
#endif
#else
#ifndef NORMALMODE
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,nlon,nlat,gflag,fg,cg,conv_step)
#else
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,nlon,nlat,gflag,fg,cg,conv_step,istep)
#endif
#endif
#else
#ifndef CONV_CHECK
#ifndef NORMALMODE
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,joff,nlon,nlat,gflag,bflag,fg,cg)
#else
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,joff,nlon,nlat,gflag,bflag,fg,cg,istep)
#endif
#else
#ifndef NORMALMODE
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,joff,nlon,nlat,gflag,bflag,fg,cg,conv_step)
#else
   subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,joff,nlon,nlat,gflag,bflag,fg,cg,conv_step,istep)
#endif
#endif
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: joff
#endif
#ifdef CONV_CHECK
      integer(kind=4), intent(out) :: conv_step
#endif
#ifdef NORMALMODE
      integer(kind=4), intent(in) :: istep
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: nm_P, nm_P0, nm_P1
      integer(kind=4), pointer, dimension(:,:) :: nm_ind
#endif

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz, p
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: rote, g, gts, theta, sint, cfc
      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy
! === Dispersive ===============================================================
      type(data_grids), target, intent(inout) :: fg, cg
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
      integer(kind=4) :: ist_, jst_
#endif
      integer(kind=4), intent(in) :: gflag
      integer(kind=4) :: ist, ind, jst, jnd
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx_old, fy_old
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: yu, yv ! RHS
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: cu, cv ! Coefficient
      integer(kind=4) :: conv                                 ! Convergence loop
#ifdef CONV_CHECK
      real(kind=REAL_BYTE) :: diff, diffx, diffy, diffxy
#ifdef MPI
      integer(kind=4) :: ierr
#endif
      real(kind=REAL_BYTE) :: tmpx, tmpy
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: tdx, tdy
#endif
      real(kind=REAL_BYTE) :: c0, c1, c2, c3, c4
#ifdef __FUJITSU
      real(kind=REAL_BYTE), allocatable, dimension(:) :: ta0, ta1
#endif
#ifdef LESS_CC
      integer(kind=4), parameter :: check = 10
#endif
! ==============================================================================

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz
! === Dispersive ===============================================================
      yu => wfld%yu
      yv => wfld%yv
      cu => wfld%cu
      cv => wfld%cv
      cu = 0.0d0
      cv = 0.0d0
      fx_old => wfld%fx_old
      fy_old => wfld%fy_old
! ==============================================================================
#ifdef CONV_CHECK
      tdx => wfld%tdx
      tdy => wfld%tdy
#endif
#ifdef NORMALMODE
      nm_ind => wfld%nm_ind
      nm_P   => wfld%nm_P
      nm_P0  => wfld%nm_P0
      nm_P1  => wfld%nm_P1
      call calc_nm_P(nlon, nlat, nm_ind, nm_P, nm_P0, nm_P1, istep)
#ifdef MPI
      call exchange_edges_P(fg)
#endif
#endif

      dz => dfld%dz
      dx => dfld%dx
      dy => dfld%dy

! === Dispersive ===============================================================
!$omp parallel
#ifndef MPI
!$omp do private(i)
      do j = 0, nlat
         do i = 0, nlon
            fx_old(i,j) = fx(i,j)
         end do
      end do
!$omp do private(i)
      do j = 0, nlat
         do i = 0, nlon
            fy_old(i,j) = fy(i,j)
         end do
      end do
#else
!$omp single
      ist_ = 1
      jst_ = 1
      if(iand(bflag, WEST_BOUND)  /= 0) ist_ = 0
      if(iand(bflag, NORTH_BOUND) /= 0) jst_ = 0
!$omp end single
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            fx_old(i,j) = fx(i,j)
         end do
      end do
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            fy_old(i,j) = fy(i,j)
         end do
      end do
#endif
!$omp end parallel
! ==============================================================================
      rote = 6371.0d0 * 1000.0d0
      g = 9.8d0
      gts = g*dt/(dth*rote)

!$omp parallel do private(theta, sint, cfc, i)
      do j = 1, nlat-1
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         cfc = gts/sint
         do i = 1, nlon-1
            if(dz(i,j) > zap) then
               fx(i,j) = fx(i,j) - dx(i,j)*cfc*(hz(i+1,j)  -hz(i,j))
               fy(i,j) = fy(i,j) - dy(i,j)*gts*(hz(i,  j+1)-hz(i,j))
#ifdef NORMALMODE
               fx(i,j) = fx(i,j) - dx(i,j)*cfc/g/m_rho*(nm_P(i+1,j)  -nm_P(i,j))
               fy(i,j) = fy(i,j) - dy(i,j)*gts/g/m_rho*(nm_P(i,  j+1)-nm_P(i,j))
#endif
            else
               fx(i,j) = zap
               fy(i,j) = zap
            end if
         end do
      end do
! === Dispersive ===============================================================
      jst = 1
      jnd = nlat - 1
      ist = 1
      ind = nlon - 1

      if(gflag == 1) then
         jst = 2
         jnd = nlat - 1
         ist = 2
         ind = nlon - 1
      end if

      jnd = jnd - jst + 1
      ind = ind - ist + 1
#ifdef MPI
      if(gflag == 1) then
         jst = 2
         jnd = nlat - 1
         ist = 2
         ind = nlon - 1
         if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 2
         if(iand(bflag, EAST_BOUND)  /= 0) ind = nlon - 2
      else
         jst = 1
         jnd = nlat - 1
         ist = 1
         ind = nlon - 1
         if(iand(bflag, NORTH_BOUND) == 0) jst = 2
         if(iand(bflag, WEST_BOUND)  == 0) ist = 2
      end if
#endif

      rote = 6371.0d0 * 1000.0d0
      g = 9.8d0
      gts = g*dt/(dth*rote)

#ifndef __FUJITSU
!$omp parallel private(conv)
#else
!$omp parallel private(conv, ta0, ta1)
      allocate(ta0(ist:ind))
      allocate(ta1(ist:ind))
#endif
!$omp do private(theta, sint, cfc, i)
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         cfc = gts/sint
         do i = ist, ind
            if(dz(i,j) > min_depth) then
               cu(i,j) = ((1.0d0/(rote*dth*sint))**2/3.0d0)*(0.5d0*(dz(i+1,j)+dz(i,j)))**2
               fx(i,j) = fx(i,j) - cu(i,j)*( &
               ! U_n
               &   fx_old(i+1,j) - 2.0d0*fx_old(i,j) + fx_old(i-1,j) &
               ! V_n
               & + (fy_old(i+1,j  )-fy_old(i,  j  ))*sin(theta+dth/2.0d0) &
               & - (fy_old(i+1,j-1)-fy_old(i,  j-1))*sin(theta-dth/2.0d0) &
               & )
               cv(i,j) = ((1.0d0/(rote*dth))**2/3.0d0)*(0.5d0*(dz(i,j+1)+dz(i,j)))**2
               fy(i,j) = fy(i,j)  - cv(i,j)*( &
               ! V_n
               &   fy_old(i  ,j+1)*sin(theta+dth)/sin(theta+dth/2.0d0) &
               & - fy_old(i  ,j  )*sint/sin(theta-dth/2.0d0)           &
               & - fy_old(i  ,j  )*sint/sin(theta+dth/2.0d0)           &
               & + fy_old(i  ,j-1)*sin(theta-dth)/sin(theta-dth/2.0d0) &
               ! U_n
               & + (fx_old(i  ,j+1)-fx_old(i-1,j+1))/sin(theta+dth/2.0d0) &
               & - (fx_old(i,  j  )-fx_old(i-1,j  ))/sin(theta-dth/2.0d0) &
               & )
            end if
         end do
      end do

#ifndef MPI
!$omp do private(i)
      do j = 0, nlat
         do i = 0, nlon
            yu(i,j) = fx(i,j)
            fx(i,j) = fx_old(i,j)
            yv(i,j) = fy(i,j)
            fy(i,j) = fy_old(i,j)
         end do
      end do
      if(gflag /= 1) then
!$omp single
         call interp2fine(VEL,cg,fg)
!$omp end single
!$omp do private(i)
         do j = 0, nlat
            do i = 0, nlon
               fx_old(i,j) = fx(i,j)
               fy_old(i,j) = fy(i,j)
            end do
         end do
      end if
#else
!$omp single
      call exchange_edges_disp_fx(fg,1)
      call exchange_edges_disp_fy(fg,1)
!$omp end single
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            yu(i,j) = fx(i,j)
            fx(i,j) = fx_old(i,j)
            yv(i,j) = fy(i,j)
            fy(i,j) = fy_old(i,j)
         end do
      end do
      if(gflag /= 1) then
!$omp single
         call interp2fine(VEL,cg,fg)
         call exchange_edges_disp_fx(fg,1)
         call exchange_edges_disp_fy(fg,1)
!$omp end single
!$omp do private(i)
         do j = jst_, nlat
            do i = ist_, nlon
               fx_old(i,j) = fx(i,j)
               fy_old(i,j) = fy(i,j)
            end do
         end do
      end if
#endif
#ifdef CONV_CHECK
!$omp do private(i, tmpx, tmpy)
      do j = 1, nlat
         do i = 1, nlon
            tmpx = 0.5d0*(dz(i+1,j)+dz(i,j)+hz(i+1,j)+hz(i,j))
            tmpy = 0.5d0*(dz(i,j+1)+dz(i,j)+hz(i,j+1)+hz(i,j))
            if(abs(tmpx) > 0.0d0) then
               tmpx = 1.0d0/tmpx
            else
               tmpx = 0.0d0
            end if
            if(abs(tmpy) > 0.0d0) then
               tmpy = 1.0d0/tmpy
            else
               tmpy = 0.0d0
            end if
            tdx(i,j) = tmpx
            tdy(i,j) = tmpy
         end do
      end do
#endif
      do conv = 1, max_step
!$omp do private(theta, i, c0, c1)
         do j = jst, jnd
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            c0 = sin(theta+0.5d0*dth)
            c1 = sin(theta-0.5d0*dth)
            do i = ist, ind
#ifndef __FUJITSU
               fx(i,j) = (yu(i,j) + cu(i,j)*( &
               ! V_n+1
               & + (fy_old(i+1,j  )-fy_old(i,  j  ))*c0 &
               & - (fy_old(i+1,j-1)-fy_old(i,  j-1))*c1 &
               ! U_n+1
               & + fx_old(i+1,j)+fx_old(i-1,j)))/(1.0d0+2.0d0*cu(i,j))
#else
               ta0(i) = yu(i,j) + cu(i,j)*( &
               & + (fy_old(i+1,j  )-fy_old(i,  j  ))*c0 &
               & - (fy_old(i+1,j-1)-fy_old(i,  j-1))*c1 &
               & + fx_old(i+1,j)+fx_old(i-1,j))
            end do
            do i = ist, ind
               ta1(i) = 1.0d0/(1.0d0+2.0d0*cu(i,j))
            end do
            do i = ist, ind
               fx(i,j) = ta0(i)*ta1(i)
#endif
            end do
         end do
#ifdef MPI
!$omp single
         call exchange_edges_disp_fx(fg,mod(conv,2))
!$omp end single
#endif
!$omp do private(theta, i, c0, c1, c2, c3, c4)
         do j = jst, jnd
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            c0 = 1.0d0/sin(theta+0.5d0*dth)
            c1 = 1.0d0/sin(theta-0.5d0*dth)
            c2 = sin(theta+dth)/sin(theta+0.5d0*dth)
            c3 = sin(theta-dth)/sin(theta-0.5d0*dth)
            c4 = sin(theta)*(1.0d0/sin(theta-0.5d0*dth)+1.0d0/sin(theta+0.5d0*dth))
            do i = ist, ind
#ifndef __FUJITSU
               fy(i,j) = (yv(i,j) + cv(i,j)*( &
               ! U_n+1
               & + (fx(i  ,j+1)-fx(i-1,j+1))*c0 &
               & - (fx(i,  j  )-fx(i-1,j  ))*c1 &
               ! V_n+1
               & + fy_old(i,j+1)*c2+fy_old(i,j-1)*c3))/(1.0d0+cv(i,j)*c4)
#else
               ta0(i) = yv(i,j) + cv(i,j)*( &
               & + (fx(i  ,j+1)-fx(i-1,j+1))*c0 &
               & - (fx(i,  j  )-fx(i-1,j  ))*c1 &
               & + fy_old(i,j+1)*c2+fy_old(i,j-1)*c3)
            end do
            do i = ist, ind
               ta1(i) = 1.0d0/(1.0d0+cv(i,j)*c4)
            end do
            do i = ist, ind
               fy(i,j) = ta0(i)*ta1(i)
#endif
            end do
         end do
#ifdef MPI
!$omp single
         call exchange_edges_disp_fy(fg,mod(conv,2))
!$omp end single
#endif
#ifdef CONV_CHECK
#ifdef LESS_CC
         if(mod(conv,check) == 0) then
#endif
!$omp single
         diffxy = 0.0d0
!$omp end single
!$omp do private(i, diff, diffx, diffy) reduction(max:diffxy)
         do j = jst, jnd
            do i = ist, ind
               diffx = abs((fx(i,j)-fx_old(i,j))*tdx(i,j))
               diffy = abs((fy(i,j)-fy_old(i,j))*tdy(i,j))
               diff = max(diffx,diffy)
               diffxy = max(diffxy,diff)
            end do
         end do
#ifdef MPI
!$omp single
#ifndef MULTI
         call MPI_Allreduce(MPI_IN_PLACE, diffxy, 1, REAL_MPI, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
         call MPI_Allreduce(MPI_IN_PLACE, diffxy, 1, REAL_MPI, MPI_MAX, MPI_MEMBER_WORLD, ierr)
#endif
!$omp end single
#endif
         if(diffxy < conv_val) then
            exit
         endif
#ifdef LESS_CC
         end if
#endif
#endif
!$omp single
         p => fx_old
         fx_old => fx
         fx => p
         p => fy_old
         fy_old => fy
         fy => p
!$omp end single
      end do
#ifdef CONV_CHECK
!$omp single
      conv_step = conv
!$omp end single
#endif
#ifndef MPI
!$omp do private(i)
      do j = 0, nlat
         do i = 0, nlon
            fx_old(i,j) = fx(i,j)
            fy_old(i,j) = fy(i,j)
         end do
      end do
#else
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            fx_old(i,j) = fx(i,j)
            fy_old(i,j) = fy(i,j)
         end do
      end do
#endif
#ifdef __FUJITSU
      deallocate(ta0)
      deallocate(ta1)
#endif
!$omp end parallel
! ==============================================================================

      return
   end subroutine fxy_rwg_disp

#ifndef MPI
#ifndef CONV_CHECK
#ifndef NORMALMODE
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,nlon,nlat, &
                             gflag,smallh,fg,cg)
#else
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,nlon,nlat, &
                             gflag,smallh,fg,cg,istep)
#endif
#else
#ifndef NORMALMODE
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,nlon,nlat, &
                             gflag,smallh,fg,cg,conv_step)
#else
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,nlon,nlat, &
                             gflag,smallh,fg,cg,conv_step,istep)
#endif
#endif
#else
#ifndef CONV_CHECK
#ifndef NORMALMODE
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,joff,nlon,nlat, &
                             gflag,smallh,bflag,fg,cg)
#else
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,joff,nlon,nlat, &
                             gflag,smallh,bflag,fg,cg,istep)
#endif
#else
#ifndef NORMALMODE
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,joff,nlon,nlat, &
                             gflag,smallh,bflag,fg,cg,conv_step)
#else
   subroutine fxynl_rwg_disp(wfld,dfld,ffld,ifz,cfs,cfl,cflag,dt,th0,dth,joff,nlon,nlat, &
                             gflag,smallh,bflag,fg,cg,conv_step,istep)
#endif
#endif
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(inout) :: dfld
      real(kind=REAL_BYTE), target, dimension(nlon,nlat), intent(in) :: ffld
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(in) :: ifz
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(in) :: ifz
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(in) :: ifz
#endif
! ==============================================================================
      real(kind=REAL_BYTE), intent(in) :: cfs, cfl
      integer(kind=4), intent(in) :: cflag
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat, gflag
      real(kind=REAL_BYTE), intent(in) :: smallh
#ifdef MPI
      integer(kind=4), intent(in) :: joff, bflag
#endif
#ifdef CONV_CHECK
      integer(kind=4), intent(out) :: conv_step
#endif

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, ddx, ddy, dz, p
      real(kind=REAL_BYTE) :: dtds, gdtds, bcf, theta, sint, invst, gdtdss, crls
      real(kind=REAL_BYTE) :: dvdx, dvdy, fric, dh, advc
#ifdef OLD_SCHEME
      real(kind=REAL_BYTE) :: cf
#else
      real(kind=REAL_BYTE) :: cf, bcf2, cf2, fric2
#endif

      real(kind=REAL_BYTE), parameter :: rote = 6371000.0d0
      real(kind=REAL_BYTE), parameter :: g = 9.8d0
      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
      real(kind=REAL_BYTE), parameter :: half = 0.5d0
      real(kind=REAL_BYTE), parameter :: quart = 0.25d0
      real(kind=REAL_BYTE), parameter :: cfreq = 7.27220522d-05 ! angular freq of earth's rotation (rad/s)

      integer(kind=4) :: i, j
      integer(kind=4) :: jst, jnd, ist, ind

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx_old, fy_old
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz_old
      real(kind=REAL_BYTE) :: fxbar, fybar
      real(kind=REAL_BYTE) :: ddx_tmp, ddy_tmp
! === Dispersive ===============================================================
      type(data_grids), target, intent(inout) :: fg, cg
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: yu, yv ! RHS
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: cu, cv ! Coefficient
      integer(kind=4) :: conv                                 ! Convergence loop
#ifdef CONV_CHECK
      real(kind=REAL_BYTE) :: diff, diffx, diffy, diffxy
#ifdef MPI
      integer(kind=4) :: ierr
#endif
      real(kind=REAL_BYTE) :: tmpx, tmpy
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: tdx, tdy
#endif
#ifndef UPWIND3
      real(kind=REAL_BYTE) :: tmp0, tmp1
#else
      real(kind=REAL_BYTE) :: tmp0, tmp1, tmp2, tmp3
#endif
      real(kind=REAL_BYTE), parameter :: small_dxy = 0.01d0
#ifdef MPI
      integer(kind=4) :: ist_, jst_
#endif
      real(kind=REAL_BYTE) :: c0, c1, c2, c3, c4
#ifdef __FUJITSU
      real(kind=REAL_BYTE), allocatable, dimension(:) :: ta0, ta1
#endif
#ifdef LESS_CC
      integer(kind=4), parameter :: check = 10
#endif
! ==============================================================================
! === Limiter with max Froude number. ==========================================
      real(kind=REAL_BYTE) :: d, lim
! ==============================================================================
#ifdef BANKFILE
      integer(kind=4), pointer, dimension(:,:) :: ir, brokenx, brokeny
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: btx, bty, dx, dy

      real(kind=REAL_BYTE) :: zhigh, zlow, discharge, dhigh
      real(kind=REAL_BYTE), parameter :: GX = 1.0d-5, GY = 1.0d-10
#endif
#ifdef NORMALMODE
      integer(kind=4), intent(in) :: istep
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: nm_P, nm_P0, nm_P1
      integer(kind=4), pointer, dimension(:,:) :: nm_ind
#endif

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz
      fx_old => wfld%fx_old
      fy_old => wfld%fy_old
      hz_old => wfld%hz_old
! === Dispersive ===============================================================
      yu => wfld%yu
      yv => wfld%yv
      cu => wfld%cu
      cv => wfld%cv
      cu = 0.0d0
      cv = 0.0d0
! ==============================================================================

#ifndef BANKFILE
      ddx => dfld%dx
      ddy => dfld%dy
#else
      if(allocated(wfld%ir)) then
         ir  => wfld%ir

         ddx => dfld%dxbx
         ddy => dfld%dyby

         btx => wfld%btx
         dx => dfld%dx
         brokenx => wfld%brokenx
         bty => wfld%bty
         dy => dfld%dy
         brokeny => wfld%brokeny
      end if
#endif
      dz  => dfld%dz
#ifdef CONV_CHECK
      tdx => wfld%tdx
      tdy => wfld%tdy
#endif
#ifdef NORMALMODE
      nm_ind => wfld%nm_ind
      nm_P   => wfld%nm_P
      nm_P0  => wfld%nm_P0
      nm_P1  => wfld%nm_P1
      call calc_nm_P(nlon, nlat, nm_ind, nm_P, nm_P0, nm_P1, istep)
#ifdef MPI
      call exchange_edges_P(fg)
#endif
#endif

      dtds = dt/(dth*rote)
      gdtds = g*dtds

      jst = 1
      jnd = nlat - 1
      ist = 1
      ind = nlon - 1

      ! if coarse grid use non-linear adjacent to boundary for absorbing stability
      if(gflag == 1) then
#ifndef UPWIND3
         jst = 3
         jnd = nlat - 2
         ist = 3
         ind = nlon - 2
#else
         jst = 4
         jnd = nlat - 3
         ist = 4
         ind = nlon - 3
#endif
      end if

#ifndef __FUJITSU
!$omp parallel private(conv)
#else
!$omp parallel private(conv, ta0, ta1)
#endif
      ! do non-linear calc within interior of grid
!$omp single
#ifdef MPI
      if(gflag == 1) then
         jst = 2
         jnd = nlat - 1
         ist = 2
         ind = nlon - 1
#ifndef UPWIND3
         if(iand(bflag, NORTH_BOUND) /= 0) jst = 3
         if(iand(bflag, WEST_BOUND)  /= 0) ist = 3
         if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 2
         if(iand(bflag, EAST_BOUND)  /= 0) ind = nlon - 2
#else
         if(iand(bflag, NORTH_BOUND) /= 0) jst = 4
         if(iand(bflag, WEST_BOUND)  /= 0) ist = 4
         if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 3
         if(iand(bflag, EAST_BOUND)  /= 0) ind = nlon - 3
#endif
      else
         jst = 1
         jnd = nlat - 1
         ist = 1
         ind = nlon - 1
         if(iand(bflag, NORTH_BOUND) == 0) jst = 2
         if(iand(bflag, WEST_BOUND)  == 0) ist = 2
      end if
#endif
!$omp end single
#ifndef MPI
!$omp do private(i)
      do j = -1, nlat+1
         do i = -2, nlon+1
            fx_old(i,j) = fx(i,j)
         end do
      end do
!$omp do private(i)
      do j = -2, nlat+1
         do i = -1, nlon+1
            fy_old(i,j) = fy(i,j)
         end do
      end do
!$omp do private(i)
      do j = -1, nlat+2
         do i = -1, nlon+2
            hz_old(i,j) = hz(i,j)
         end do
      end do
#else
!$omp single
      ist_ = -1
      jst_ = 0
      if(iand(bflag, WEST_BOUND)  /= 0) ist_ = -2
      if(iand(bflag, NORTH_BOUND) /= 0) jst_ = -1
!$omp end single
!$omp do private(i)
      do j = jst_, nlat+1
         do i = ist_, nlon+1
            fx_old(i,j) = fx(i,j)
         end do
      end do
!$omp single
      ist_ = 0
      jst_ = -1
      if(iand(bflag, WEST_BOUND)  /= 0) ist_ = -1
      if(iand(bflag, NORTH_BOUND) /= 0) jst_ = -2
!$omp end single
!$omp do private(i)
      do j = jst_, nlat+1
         do i = ist_, nlon+1
            fy_old(i,j) = fy(i,j)
         end do
      end do
!$omp single
      ist_ = 0
      jst_ = 0
      if(iand(bflag, WEST_BOUND)  /= 0) ist_ = -1
      if(iand(bflag, NORTH_BOUND) /= 0) jst_ = -1
!$omp end single
!$omp do private(i)
      do j = jst_, nlat+2
         do i = ist_, nlon+2
            hz_old(i,j) = hz(i,j)
         end do
      end do
#endif

! ==============================================================================
#ifndef UPWIND3
!$omp do &
!$omp private(theta, sint, invst, gdtdss, crls, i, fybar, &
!$omp         dvdx, dvdy, advc, tmp0, tmp1)
#else
!$omp do &
!$omp private(theta, sint, invst, gdtdss, crls, i, fybar, &
!$omp         dvdx, dvdy, advc, tmp0, tmp1, tmp2, tmp3)
#endif
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         invst = 1.0d0/sint
         gdtdss = gdtds*invst
         crls = 2.0d0*cflag*cfreq*cos(theta)*dt

         do i = ist, ind
            fybar = quart*(fy_old(i,j) + fy_old(i+1,j) + fy_old(i,j-1) + fy_old(i+1,j-1))

            if(ifz(i,j)+ifz(i+1,j) > 0) then ! both wet
#ifndef UPWIND3
               if(fx_old(i,j) >= zap) then
                  tmp0 = half*(dz(i+1,j)+dz(i,j)+hz_old(i+1,j)+hz_old(i,j))
                  tmp1 = half*(dz(i-1,j)+dz(i,j)+hz_old(i-1,j)+hz_old(i,j))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fx_old(i,j) /= zap .and. fx_old(i-1,j) /= zap) then
                     dvdx = fx_old(i,j)**2/tmp0 - fx_old(i-1,j)**2/tmp1
                  else
                     dvdx = zap
                  end if
               else
                  tmp0 = half*(dz(i+2,j)+dz(i+1,j)+hz_old(i+2,j)+hz_old(i+1,j))
                  tmp1 = half*(dz(i+1,j)+dz(i,  j)+hz_old(i+1,j)+hz_old(i,  j))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fx_old(i+1,j) /= zap .and. fx_old(i,j) /= zap) then
                     dvdx = fx_old(i+1,j)**2/tmp0 - fx_old(i,j)**2/tmp1
                  else
                     dvdx = zap
                  end if
               end if
               if(fybar >= zap) then
                  tmp0 = half*(dz(i+1,j  )+dz(i,j  )+hz_old(i+1,j  )+hz_old(i,j  ))
                  tmp1 = half*(dz(i+1,j-1)+dz(i,j-1)+hz_old(i+1,j-1)+hz_old(i,j-1))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fx_old(i,j  ) /= zap .and. fx_old(i,  j-1) /= zap .and. &
                     fy_old(i,j  ) /= zap .and. fy_old(i+1,j  ) /= zap .and. &
                     fy_old(i,j-1) /= zap .and. fy_old(i+1,j-1) /= zap .and. &
                     fy_old(i,j-2) /= zap .and. fy_old(i+1,j-2) /= zap) then
                     dvdy = fx_old(i,j  )*(fy_old(i,j  )+fy_old(i+1,j  )+fy_old(i,j-1)+fy_old(i+1,j-1))/tmp0 &
                          - fx_old(i,j-1)*(fy_old(i,j-1)+fy_old(i+1,j-1)+fy_old(i,j-2)+fy_old(i+1,j-2))/tmp1
                  else
                     dvdy = zap
                  end if
               else
                  tmp0 = half*(dz(i+1,j+1)+dz(i,j+1)+hz_old(i+1,j+1)+hz_old(i,j+1))
                  tmp1 = half*(dz(i+1,j  )+dz(i,j  )+hz_old(i+1,j  )+hz_old(i,j  ))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fx_old(i,j+1) /= zap .and. fx_old(i,  j  ) /= zap .and. &
                     fy_old(i,j+1) /= zap .and. fy_old(i+1,j+1) /= zap .and. &
                     fy_old(i,j  ) /= zap .and. fy_old(i+1,j  ) /= zap .and. &
                     fy_old(i,j-1) /= zap .and. fy_old(i+1,j-1) /= zap) then
                     dvdy = fx_old(i,j+1)*(fy_old(i,j+1)+fy_old(i+1,j+1)+fy_old(i,j  )+fy_old(i+1,j  ))/tmp0 &
                          - fx_old(i,j  )*(fy_old(i,j  )+fy_old(i+1,j  )+fy_old(i,j-1)+fy_old(i+1,j-1))/tmp1
                  else
                     dvdy = zap
                  end if
               end if
               advc = dtds*(dvdx*invst + quart*dvdy)
#else
               if(fx_old(i,j) >= zap) then
                  tmp0 = 0.5d0*(dz(i-2,j)+dz(i-1,j)+hz_old(i-2,j)+hz_old(i-1,j))
                  tmp1 = 0.5d0*(dz(i-1,j)+dz(i,  j)+hz_old(i-1,j)+hz_old(i,  j))
                  tmp2 = 0.5d0*(dz(i,  j)+dz(i+1,j)+hz_old(i,  j)+hz_old(i+1,j))
                  tmp3 = 0.5d0*(dz(i+1,j)+dz(i+2,j)+hz_old(i+1,j)+hz_old(i+2,j))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fx_old(i-2,j) /= zap .and. fx_old(i-1,j) /= zap .and. &
                     fx_old(i,  j) /= zap .and. fx_old(i+1,j) /= zap) then
                     dvdx =         fx_old(i-2,j)**2/tmp0 &
                            - 6.0d0*fx_old(i-1,j)**2/tmp1 &
                            + 3.0d0*fx_old(i,  j)**2/tmp2 &
                            + 2.0d0*fx_old(i+1,j)**2/tmp3
                  else
                     dvdx = zap
                  end if
               else
                  tmp0 = 0.5d0*(dz(i+2,j)+dz(i+3,j)+hz_old(i+2,j)+hz_old(i+3,j))
                  tmp1 = 0.5d0*(dz(i+1,j)+dz(i+2,j)+hz_old(i+1,j)+hz_old(i+2,j))
                  tmp2 = 0.5d0*(dz(i,  j)+dz(i+1,j)+hz_old(i,  j)+hz_old(i+1,j))
                  tmp3 = 0.5d0*(dz(i-1,j)+dz(i,  j)+hz_old(i-1,j)+hz_old(i,  j))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fx_old(i+2,j) /= zap .and. fx_old(i+1,j) /= zap .and. &
                     fx_old(i,  j) /= zap .and. fx_old(i-1,j) /= zap) then
                     dvdx = -       fx_old(i+2,j)**2/tmp0 &
                            + 6.0d0*fx_old(i+1,j)**2/tmp1 &
                            - 3.0d0*fx_old(i,  j)**2/tmp2 &
                            - 2.0d0*fx_old(i-1,j)**2/tmp3
                  else
                     dvdx = zap
                  end if
               end if
               if(fybar >= zap) then
                  tmp0 = 0.5d0*(dz(i,j-2)+dz(i+1,j-2)+hz_old(i,j-2)+hz_old(i+1,j-2))
                  tmp1 = 0.5d0*(dz(i,j-1)+dz(i+1,j-1)+hz_old(i,j-1)+hz_old(i+1,j-1))
                  tmp2 = 0.5d0*(dz(i,j  )+dz(i+1,j  )+hz_old(i,j  )+hz_old(i+1,j  ))
                  tmp3 = 0.5d0*(dz(i,j+1)+dz(i+1,j+1)+hz_old(i,j+1)+hz_old(i+1,j+1))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fx_old(i,j-2) /= zap .and. fx_old(i,  j-1) /= zap .and. &
                     fx_old(i,j  ) /= zap .and. fx_old(i,  j+1) /= zap .and. &
                     fy_old(i,j+1) /= zap .and. fy_old(i+1,j+1) /= zap .and. &
                     fy_old(i,j  ) /= zap .and. fy_old(i+1,j  ) /= zap .and. &
                     fy_old(i,j-1) /= zap .and. fy_old(i+1,j-1) /= zap .and. &
                     fy_old(i,j-2) /= zap .and. fy_old(i+1,j-2) /= zap .and. &
                     fy_old(i,j-3) /= zap .and. fy_old(i+1,j-3) /= zap) then
                     dvdy =         fx_old(i,j-2)*(fy_old(i,j-3)+fy_old(i+1,j-3)+fy_old(i,j-2)+fy_old(i+1,j-2))/tmp0 &
                            - 6.0d0*fx_old(i,j-1)*(fy_old(i,j-2)+fy_old(i+1,j-2)+fy_old(i,j-1)+fy_old(i+1,j-1))/tmp1 &
                            + 3.0d0*fx_old(i,j  )*(fy_old(i,j-1)+fy_old(i+1,j-1)+fy_old(i,j  )+fy_old(i+1,j  ))/tmp2 &
                            + 2.0d0*fx_old(i,j+1)*(fy_old(i,j  )+fy_old(i+1,j  )+fy_old(i,j+1)+fy_old(i+1,j+1))/tmp3
                  else
                     dvdy = zap
                  end if
               else
                  tmp0 = 0.5d0*(dz(i,j+2)+dz(i+1,j+2)+hz_old(i,j+2)+hz_old(i+1,j+2))
                  tmp1 = 0.5d0*(dz(i,j+1)+dz(i+1,j+1)+hz_old(i,j+1)+hz_old(i+1,j+1))
                  tmp2 = 0.5d0*(dz(i,j  )+dz(i+1,j  )+hz_old(i,j  )+hz_old(i+1,j  ))
                  tmp3 = 0.5d0*(dz(i,j-1)+dz(i+1,j-1)+hz_old(i,j-1)+hz_old(i+1,j-1))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fx_old(i,j+2) /= zap .and. fx_old(i,  j+1) /= zap .and. &
                     fx_old(i,j  ) /= zap .and. fx_old(i,  j-1) /= zap .and. &
                     fy_old(i,j+2) /= zap .and. fy_old(i+1,j+2) /= zap .and. &
                     fy_old(i,j+1) /= zap .and. fy_old(i+1,j+1) /= zap .and. &
                     fy_old(i,j  ) /= zap .and. fy_old(i+1,j  ) /= zap .and. &
                     fy_old(i,j-1) /= zap .and. fy_old(i+1,j-1) /= zap .and. &
                     fy_old(i,j-2) /= zap .and. fy_old(i+1,j-2) /= zap) then
                     dvdy = -       fx_old(i,j+2)*(fy_old(i,j+1)+fy_old(i+1,j+1)+fy_old(i,j+2)+fy_old(i+1,j+2))/tmp0 &
                            + 6.0d0*fx_old(i,j+1)*(fy_old(i,j  )+fy_old(i+1,j  )+fy_old(i,j+1)+fy_old(i+1,j+1))/tmp1 &
                            - 3.0d0*fx_old(i,j  )*(fy_old(i,j-1)+fy_old(i+1,j-1)+fy_old(i,j  )+fy_old(i+1,j  ))/tmp2 &
                            - 2.0d0*fx_old(i,j-1)*(fy_old(i,j-2)+fy_old(i+1,j-2)+fy_old(i,j-1)+fy_old(i+1,j-1))/tmp3
                  else
                     dvdy = zap
                  end if
               end if
               advc = dtds/6.0d0*(dvdx*invst + quart*dvdy)
#endif
               fx(i,j) = fx_old(i,j)                                                           &
                       - half*gdtdss*(dz(i+1,j)+dz(i,j)+hz_old(i+1,j)+hz_old(i,j))*(hz_old(i+1,j)-hz_old(i,j)) &
#ifndef NORMALMODE
                       - advc - crls*fybar
#else
                       - advc - crls*fybar - half*(dz(i+1,j)+dz(i,j))*dtds*invst/m_rho*(nm_P(i+1,j)-nm_P(i,j))
#endif
            end if
         end do
      end do
! ==============================================================================
!$omp do private(i)
      do j = jst, jnd
         do i = ist, ind
            if(ifz(i,j)+ifz(i+1,j) < 0) then ! both dry
               fx(i,j) = zap
            end if
         end do
      end do
! ==============================================================================
!$omp do private(theta, sint, invst, gdtdss, i, dh)
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         invst = 1.0d0/sint
         gdtdss = gdtds*invst

         do i = ist, ind
            if((ifz(i,j) > 0) .and. (ifz(i+1,j) < 0)) then ! wet west, dry east
               dh = dz(i+1,j) + hz_old(i,j)
               if(dh <= zap) then
                  fx(i,j) = zap
               else if((-dz(i+1,j)-hz_old(i,j)) > 10.0d0) then
                  fx(i,j) = - gdtdss*dh*10.0d0
               else if((-dz(i+1,j)-hz_old(i,j)) < -10.0d0) then
                  fx(i,j) = - gdtdss*dh*(-10.0d0)
               else
                  fx(i,j) = - gdtdss*dh*(-dz(i+1,j)-hz_old(i,j))
               end if
               if(dh > zap) then
                  if((-dz(i+1,j)-hz_old(i,j)) < 0.0d0) then
                     hz(i+1,j) = -dz(i+1,j)
                  end if
               end if
            end if
         end do
      end do
! ==============================================================================
!$omp do private(theta, sint, invst, gdtdss, i, dh)
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         invst = 1.0d0/sint
         gdtdss = gdtds*invst

         do i = ist, ind
            if((ifz(i,j) < 0) .and. (ifz(i+1,j) > 0)) then ! dry west, wet east
               dh = dz(i,j) + hz_old(i+1,j)
               if(dh <= zap) then
                  fx(i,j) = zap
               else if((hz_old(i+1,j)+dz(i,j)) > 10.0d0) then
                  fx(i,j) = - gdtdss*dh*10.0d0
               else if((hz_old(i+1,j)+dz(i,j)) < -10.0d0) then
                  fx(i,j) = - gdtdss*dh*(-10.0d0)
               else
                  fx(i,j) = - gdtdss*dh*(hz_old(i+1,j)+dz(i,j))
               end if
               if(dh > zap) then
                  if((hz_old(i+1,j)+dz(i,j)) > 0.0d0) then
                     hz(i,j) = -dz(i,j)
                  end if
               end if
            end if
         end do
      end do
! ==============================================================================
#ifdef BANKFILE
      if(allocated(wfld%ir)) then
!$omp do private(i, zhigh, zlow, discharge, dhigh)
         do j = jst, jnd
            do i = ist, ind
               if((ir(i,j) == 1) .or. (ir(i,j) == 3)) then
                  dhigh = min(dz(i,j), dz(i+1,j))
! === CRITICAL! Almost all line-data will be ignored!!! ========================
!                 if(ddx(i,j) > -dhigh) then
                  if(ddx(i,j) < dhigh) then
! ==============================================================================
                     ! === CAL. OF DISCHANGE OF OVERFLOW ===
                     if(hz(i,j) + ddx(i,j) > hz(i+1,j) + ddx(i,j)) then
                        zhigh = hz(i,j) + ddx(i,j)
                        zlow = hz(i+1,j) + ddx(i,j)
                     else
                        zhigh = hz(i+1,j) + ddx(i,j)
                        zlow = hz(i,j) + ddx(i,j)
                     end if
                     if(zhigh < GX) then
                        fx(i,j) = zap
                     else
                        if(zhigh*0.66667d0 < zlow) then
                           discharge = 4.029d0*zlow*sqrt(zhigh - zlow) ! 4.029 = u'*(9.8*2)*0.5, u' = 2.6u
                        else
                           discharge = 1.55d0*zhigh**1.5d0 ! 1.55 = u*(9.8*2)*0.5, u = 0.35
                        end if
                        if(hz(i+1,j) + ddx(i,j) > hz(i,j) + ddx(i,j)) discharge = -discharge
                        fx(i,j) = discharge
                        if(brokenx(i,j) == 0) then
                           btx(i,j) = broken_rate*(btx(i,j) + dhigh) - dhigh
                           brokenx(i,j) = 1
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
#endif
#ifndef UPWIND3
!$omp do &
!$omp private(theta, sint, invst, crls, i, fxbar, &
!$omp         dvdx, dvdy, advc, tmp0, tmp1)
#else
!$omp do &
!$omp private(theta, sint, invst, crls, i, fxbar, &
!$omp         dvdx, dvdy, advc, tmp0, tmp1, tmp2, tmp3)
#endif
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         invst = 1.0d0/sint
         crls = 2.0d0*cflag*cfreq*cos(theta)*dt

         do i = ist, ind
            fxbar = quart*(fx_old(i,j) + fx_old(i-1,j) + fx_old(i-1,j+1) + fx_old(i,j+1))

            if(ifz(i,j)+ifz(i,j+1) > 0) then ! both wet
#ifndef UPWIND3
               if(fxbar >= zap) then
                  tmp0 = half*(dz(i,  j+1)+dz(i,  j)+hz_old(i,  j+1)+hz_old(i,  j))
                  tmp1 = half*(dz(i-1,j+1)+dz(i-1,j)+hz_old(i-1,j+1)+hz_old(i-1,j))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fy_old(i,  j) /= zap .and. fy_old(i-1,j  ) /= zap .and. &
                     fx_old(i,  j) /= zap .and. fx_old(i,  j+1) /= zap .and. &
                     fx_old(i-1,j) /= zap .and. fx_old(i-1,j+1) /= zap .and. &
                     fx_old(i-2,j) /= zap .and. fx_old(i-2,j+1) /= zap) then
                     dvdx = (fx_old(i,  j)+fx_old(i,  j+1)+fx_old(i-1,j)+fx_old(i-1,j+1))*fy_old(i,  j)/tmp0 &
                          - (fx_old(i-1,j)+fx_old(i-1,j+1)+fx_old(i-2,j)+fx_old(i-2,j+1))*fy_old(i-1,j)/tmp1
                  else
                     dvdx = zap
                  end if
               else
                  tmp0 = half*(dz(i+1,j+1)+dz(i+1,j)+hz_old(i+1,j+1)+hz_old(i+1,j))
                  tmp1 = half*(dz(i,  j+1)+dz(i,  j)+hz_old(i,  j+1)+hz_old(i,  j))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fy_old(i+1,j) /= zap .and. fy_old(i,  j  ) /= zap .and. &
                     fx_old(i+1,j) /= zap .and. fx_old(i+1,j+1) /= zap .and. &
                     fx_old(i,  j) /= zap .and. fx_old(i,  j+1) /= zap .and. &
                     fx_old(i-1,j) /= zap .and. fx_old(i-1,j+1) /= zap) then
                     dvdx = (fx_old(i+1,j)+fx_old(i+1,j+1)+fx_old(i,  j)+fx_old(i,  j+1))*fy_old(i+1,j)/tmp0 &
                          - (fx_old(i,  j)+fx_old(i,  j+1)+fx_old(i-1,j)+fx_old(i-1,j+1))*fy_old(i,  j)/tmp1
                  else
                     dvdx = zap
                  end if
               end if
               if(fy_old(i,j) >= zap) then
                  tmp0 = half*(dz(i,j+1)+dz(i,j)+hz_old(i,j+1)+hz_old(i,j))
                  tmp1 = half*(dz(i,j-1)+dz(i,j)+hz_old(i,j-1)+hz_old(i,j))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fy_old(i,j) /= zap .and. fy_old(i,j-1) /= zap) then
                     dvdy = fy_old(i,j)**2/tmp0 - fy_old(i,j-1)**2/tmp1
                  else
                     dvdy = zap
                  end if
               else
                  tmp0 = half*(dz(i,j+2)+dz(i,j+1)+hz_old(i,j+2)+hz_old(i,j+1))
                  tmp1 = half*(dz(i,j+1)+dz(i,j  )+hz_old(i,j+1)+hz_old(i,j  ))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     fy_old(i,j+1) /= zap .and. fy_old(i,j) /= zap) then
                     dvdy = fy_old(i,j+1)**2/tmp0 - fy_old(i,j)**2/tmp1
                  else
                     dvdy = zap
                  end if
               end if
               advc = dtds*(quart*dvdx*invst + dvdy)
#else
               if(fxbar >= zap) then
                  tmp0 = 0.5d0*(dz(i-2,j)+dz(i-2,j+1)+hz_old(i-2,j)+hz_old(i-2,j+1))
                  tmp1 = 0.5d0*(dz(i-1,j)+dz(i-1,j+1)+hz_old(i-1,j)+hz_old(i-1,j+1))
                  tmp2 = 0.5d0*(dz(i,  j)+dz(i,  j+1)+hz_old(i,  j)+hz_old(i,  j+1))
                  tmp3 = 0.5d0*(dz(i+1,j)+dz(i+1,j+1)+hz_old(i+1,j)+hz_old(i+1,j+1))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fy_old(i-2,j) /= zap .and. fy_old(i-1,j  ) /= zap .and. &
                     fy_old(i,  j) /= zap .and. fy_old(i+1,j  ) /= zap .and. &
                     fx_old(i+1,j) /= zap .and. fx_old(i+1,j+1) /= zap .and. &
                     fx_old(i,  j) /= zap .and. fx_old(i,  j+1) /= zap .and. &
                     fx_old(i-1,j) /= zap .and. fx_old(i-1,j+1) /= zap .and. &
                     fx_old(i-2,j) /= zap .and. fx_old(i-2,j+1) /= zap .and. &
                     fx_old(i-3,j) /= zap .and. fx_old(i-3,j+1) /= zap) then
                     dvdx =         fy_old(i-2,j)*(fx_old(i-3,j)+fx_old(i-2,j)+fx_old(i-3,j+1)+fx_old(i-2,j+1))/tmp0 &
                            - 6.0d0*fy_old(i-1,j)*(fx_old(i-2,j)+fx_old(i-1,j)+fx_old(i-2,j+1)+fx_old(i-1,j+1))/tmp1 &
                            + 3.0d0*fy_old(i,  j)*(fx_old(i-1,j)+fx_old(i,  j)+fx_old(i-1,j+1)+fx_old(i,  j+1))/tmp2 &
                            + 2.0d0*fy_old(i+1,j)*(fx_old(i,  j)+fx_old(i+1,j)+fx_old(i,  j+1)+fx_old(i+1,j+1))/tmp3
                  else
                     dvdx = zap
                  end if
               else
                  tmp0 = 0.5d0*(dz(i+2,j)+dz(i+2,j+1)+hz_old(i+2,j)+hz_old(i+2,j+1))
                  tmp1 = 0.5d0*(dz(i+1,j)+dz(i+1,j+1)+hz_old(i+1,j)+hz_old(i+1,j+1))
                  tmp2 = 0.5d0*(dz(i,  j)+dz(i,  j+1)+hz_old(i,  j)+hz_old(i,  j+1))
                  tmp3 = 0.5d0*(dz(i-1,j)+dz(i-1,j+1)+hz_old(i-1,j)+hz_old(i-1,j+1))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fy_old(i+2,j) /= zap .and. fy_old(i+1,j  ) /= zap .and. &
                     fy_old(i,  j) /= zap .and. fy_old(i-1,j  ) /= zap .and. &
                     fx_old(i+2,j) /= zap .and. fx_old(i+2,j+1) /= zap .and. &
                     fx_old(i+1,j) /= zap .and. fx_old(i+1,j+1) /= zap .and. &
                     fx_old(i,  j) /= zap .and. fx_old(i,  j+1) /= zap .and. &
                     fx_old(i-1,j) /= zap .and. fx_old(i-1,j+1) /= zap .and. &
                     fx_old(i-2,j) /= zap .and. fx_old(i-2,j+1) /= zap) then
                     dvdx = -       fy_old(i+2,j)*(fx_old(i+1,j)+fx_old(i+2,j)+fx_old(i+1,j+1)+fx_old(i+2,j+1))/tmp0 &
                            + 6.0d0*fy_old(i+1,j)*(fx_old(i,  j)+fx_old(i+1,j)+fx_old(i,  j+1)+fx_old(i+1,j+1))/tmp1 &
                            - 3.0d0*fy_old(i,  j)*(fx_old(i-1,j)+fx_old(i,  j)+fx_old(i-1,j+1)+fx_old(i,  j+1))/tmp2 &
                            - 2.0d0*fy_old(i-1,j)*(fx_old(i-2,j)+fx_old(i-1,j)+fx_old(i-2,j+1)+fx_old(i-1,j+1))/tmp3
                  else
                     dvdx = zap
                  end if
               end if
               if(fy_old(i,j) >= zap) then
                  tmp0 = 0.5d0*(dz(i,j-2)+dz(i,j-1)+hz_old(i,j-2)+hz_old(i,j-1))
                  tmp1 = 0.5d0*(dz(i,j-1)+dz(i,j  )+hz_old(i,j-1)+hz_old(i,j  ))
                  tmp2 = 0.5d0*(dz(i,j  )+dz(i,j+1)+hz_old(i,j  )+hz_old(i,j+1))
                  tmp3 = 0.5d0*(dz(i,j+1)+dz(i,j+2)+hz_old(i,j+1)+hz_old(i,j+2))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fy_old(i,j-2) /= zap .and. fy_old(i,j-1) /= zap .and. &
                     fy_old(i,j  ) /= zap .and. fy_old(i,j+1) /= zap) then
                     dvdy =         fy_old(i,j-2)**2/tmp0 &
                            - 6.0d0*fy_old(i,j-1)**2/tmp1 &
                            + 3.0d0*fy_old(i,j  )**2/tmp2 &
                            + 2.0d0*fy_old(i,j+1)**2/tmp3
                  else
                     dvdy = zap
                  end if
               else
                  tmp0 = 0.5d0*(dz(i,j+2)+dz(i,j+3)+hz_old(i,j+2)+hz_old(i,j+3))
                  tmp1 = 0.5d0*(dz(i,j+1)+dz(i,j+2)+hz_old(i,j+1)+hz_old(i,j+2))
                  tmp2 = 0.5d0*(dz(i,j  )+dz(i,j+1)+hz_old(i,j  )+hz_old(i,j+1))
                  tmp3 = 0.5d0*(dz(i,j-1)+dz(i,j  )+hz_old(i,j-1)+hz_old(i,j  ))
                  if(tmp0 > small_dxy .and. tmp1 > small_dxy .and. &
                     tmp2 > small_dxy .and. tmp3 > small_dxy .and. &
                     fy_old(i,j+2) /= zap .and. fy_old(i,j+1) /= zap .and. &
                     fy_old(i,j  ) /= zap .and. fy_old(i,j-1) /= zap) then
                     dvdy = -       fy_old(i,j+2)**2/tmp0 &
                            + 6.0d0*fy_old(i,j+1)**2/tmp1 &
                            - 3.0d0*fy_old(i,j  )**2/tmp2 &
                            - 2.0d0*fy_old(i,j-1)**2/tmp3
                  else
                     dvdy = zap
                  end if
               end if
               advc = dtds/6.0d0*(quart*dvdx*invst + dvdy)
#endif
               fy(i,j) = fy_old(i,j)                                                          &
                       - half*gdtds*(dz(i,j+1)+dz(i,j)+hz_old(i,j+1)+hz_old(i,j))*(hz_old(i,j+1)-hz_old(i,j)) &
#ifndef NORMALMODE
                       - advc + crls*fxbar
#else
                       - advc + crls*fxbar - half*(dz(i,j+1)+dz(i,j))*dtds/m_rho*(nm_P(i,j+1)-nm_P(i,j))
#endif
            end if
         end do
      end do
! ==============================================================================
!$omp do private(i)
      do j = jst, jnd
         do i = ist, ind
            if(ifz(i,j)+ifz(i,j+1) < 0) then ! both dry
               fy(i,j) = zap
            end if
         end do
      end do
! ==============================================================================
!$omp do private(i, dh)
      do j = jst, jnd
         do i = ist, ind
            if((ifz(i,j) > 0) .and. (ifz(i,j+1) < 0)) then ! wet north, dry south
               dh = dz(i,j+1) + hz_old(i,j)
               if(dh <= zap) then
                  fy(i,j) = zap
               else if((-dz(i,j+1)-hz_old(i,j)) > 10.0d0) then
                  fy(i,j) = - gdtds*dh*10.0d0
               else if((-dz(i,j+1)-hz_old(i,j)) < -10.0d0) then
                  fy(i,j) = - gdtds*dh*(-10.0d0)
               else
                  fy(i,j) = - gdtds*dh*(-dz(i,j+1)-hz_old(i,j))
               end if
               if(dh > zap) then
                  if((-dz(i,j+1)-hz_old(i,j)) < 0.0d0) then
                     hz(i,j+1) = -dz(i,j+1)
                  end if
               end if
            end if
         end do
      end do
! ==============================================================================
!$omp do private(i, dh)
      do j = jst, jnd
         do i = ist, ind
            if((ifz(i,j) < 0) .and. (ifz(i,j+1) > 0)) then ! dry north, wet south
               dh = dz(i,j) + hz_old(i,j+1)
               if(dh <= zap) then
                  fy(i,j) = zap
               else if((hz_old(i,j+1)+dz(i,j)) > 10.0d0) then
                  fy(i,j) = - gdtds*dh*10.0d0
               else if((hz_old(i,j+1)+dz(i,j)) < -10.0d0) then
                  fy(i,j) = - gdtds*dh*(-10.0d0)
               else
                  fy(i,j) = - gdtds*dh*(hz_old(i,j+1)+dz(i,j))
               end if
               if(dh > zap) then
                  if((hz_old(i,j+1)+dz(i,j)) > 0.0d0) then
                     hz(i,j) = -dz(i,j)
                  end if
               end if
            end if
         end do
      end do
! ==============================================================================
#ifdef BANKFILE
      if(allocated(wfld%ir)) then
!$omp do private(i, zhigh, zlow, discharge, dhigh)
         do j = jst, jnd
            do i = ist, ind
               if((ir(i,j+1) == 2) .or. (ir(i,j+1) == 3)) then
                  dhigh = min(dz(i,j), dz(i,j+1))
! === CRITICAL! Almost all line-data will be ignored!!! ========================
!                 if(ddy(i,j) > -dhigh) then
                  if(ddy(i,j) < dhigh) then
! ==============================================================================
                     ! === CAL. OF DISCHANGE OF OVERFLOW ===
                     if(hz(i,j) + ddy(i,j) > hz(i,j+1) + ddy(i,j)) then
                        zhigh = hz(i,j) + ddy(i,j)
                        zlow = hz(i,j+1) + ddy(i,j)
                     else
                        zhigh = hz(i,j+1) + ddy(i,j)
                        zlow = hz(i,j) + ddy(i,j)
                     end if
                     if(zhigh < GX) then
                        fy(i,j) = zap
                     else
                        if(zhigh*0.66667d0 < zlow) then
                           discharge = 4.029d0*zlow*sqrt(zhigh - zlow) ! 4.029 = u'*(9.8*2)*0.5, u' = 2.6u
                        else
                           discharge = 1.55d0*zhigh**1.5d0 ! 1.55 = u*(9.8*2)*0.5, u = 0.35
                        end if
                        if(hz(i,j+1) + ddy(i,j) > hz(i,j) + ddy(i,j)) discharge = -discharge
                        fy(i,j) = discharge
                        if(brokeny(i,j) == 0) then
                           bty(i,j) = broken_rate*(bty(i,j) + dhigh) - dhigh
                           brokeny(i,j) = 1
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
#endif
!$omp do private(theta, sint, i, fybar, cf, bcf, fric, ddx_tmp, d, lim)
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         do i = ist, ind
            fybar = quart*(fy_old(i,j) + fy_old(i+1,j) + fy_old(i,j-1) + fy_old(i+1,j-1))
            if(ifz(i,j)+ifz(i+1,j) > 0) then ! both wet
               ddx_tmp = half*(hz_old(i+1,j) + hz_old(i,j) + dz(i+1,j) + dz(i,j))
               if(ddx_tmp > smallh) then
                  if(ffld(i,j) == 0.0d0) then ! case of NO_FRICTION_FILE_GIVEN
                     if(half*(dz(i+1,j)+dz(i,j)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else ! case of FRICTION_FILE_GIVEN
                     cf = ffld(i,j)
                  end if
                  if(cf > 0.0d0) then ! non-dimensional friction coefficent
                     bcf = cf
                  else  ! Manning's roughness coefficent
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  ! explicit
!                 fric = dt*bcf*fx_old(i,j)*sqrt(fx_old(i,j)*fx_old(i,j) + fybar*fybar)/ddx_tmp**2
#ifdef OLD_SCHEME
                  ! semi-implicit added on 14/04/2019, Baba
                  fric = dt*bcf*half*sqrt(fx_old(i,j)*fx_old(i,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(i,j) = (fx(i,j)-fric*fx_old(i,j))/(1.0d0+fric)
#else
                  ! combined semi-implicit added on 14/04/2019, Baba
!                 fric = dt*bcf*half*sqrt(fx_old(i,j)*fx_old(i,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  ! simple semi-implicit added on 16/05/2023, Minami
                  fric = dt*bcf*sqrt(fx_old(i,j)*fx_old(i,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(i,j) = fx(i,j)/(1.0d0+fric)
#endif
               else
!                  fric = zap
                  fx(i,j) = zap
               end if
!               fx(i,j) = fx(i,j) - fric
! === Limiter with max Froude number. ==========================================
! === To prevent sqrt of negative numbers. =====================================
!              d = half*(dz(i+1,j)+dz(i,j)+hz_old(i+1,j)+hz_old(i,j))
               d = max(0.0d0, half*(dz(i+1,j)+dz(i,j)+hz_old(i+1,j)+hz_old(i,j)))
! ==============================================================================
               lim = froude_lim*d*sqrt(9.8d0*d)
               if(fx(i,j) > lim) then
                  fx(i,j) = lim
               else if(fx(i,j) < -lim) then
                  fx(i,j) = -lim
               end if
! ==============================================================================
#ifdef BANKFILE
            end if
         end do
      end do

      if(allocated(wfld%ir)) then
!$omp do private(theta, sint, i)
         do j = jst, jnd
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            sint = sin(theta)
            do i = ist, ind
               if(ifz(i,j)+ifz(i+1,j) > 0) then ! both wet
                  if(dz(i,j) > min_depth) then
                     if((ir(i,j) /= 1) .and. (ir(i,j) /= 3)) then
                        cu(i,j) = ((1.0d0/(rote*dth*sint))**2/3.0d0)*(0.5d0*(dz(i+1,j)+dz(i,j)))**2
                        fx(i,j) = fx(i,j) - cu(i,j)*( &
                           ! U_n
                           &   fx_old(i+1,j) - 2.0d0*fx_old(i,j) + fx_old(i-1,j) &
                           ! V_n
                           & + (fy_old(i+1,j  )-fy_old(i,  j  ))*sin(theta+dth/2.0d0) &
                           & - (fy_old(i+1,j-1)-fy_old(i,  j-1))*sin(theta-dth/2.0d0) &
                           & )
                     end if
                  end if
               end if
            end do
         end do
      else
!$omp do private(theta, sint, i)
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         do i = ist, ind
            if(ifz(i,j)+ifz(i+1,j) > 0) then ! both wet
#endif
! === Dispersive ===============================================================
               if(dz(i,j) > min_depth) then
                  cu(i,j) = ((1.0d0/(rote*dth*sint))**2/3.0d0)*(0.5d0*(dz(i+1,j)+dz(i,j)))**2
                  fx(i,j) = fx(i,j) - cu(i,j)*( &
                          ! U_n
                          &   fx_old(i+1,j) - 2.0d0*fx_old(i,j) + fx_old(i-1,j) &
                          ! V_n
                          & + (fy_old(i+1,j  )-fy_old(i,  j  ))*sin(theta+dth/2.0d0) &
                          & - (fy_old(i+1,j-1)-fy_old(i,  j-1))*sin(theta-dth/2.0d0) &
                          & )
               end if
! ==============================================================================
            end if
         end do
      end do
! ==============================================================================
#ifdef BANKFILE
      end if
#endif
!$omp do private(theta, sint, i, fxbar, cf, bcf, fric, ddy_tmp, d, lim)
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         do i = ist, ind
            fxbar = quart*(fx_old(i,j) + fx_old(i-1,j) + fx_old(i-1,j+1) + fx_old(i,j+1))
            if(ifz(i,j)+ifz(i,j+1) > 0) then ! both wet
               ddy_tmp = half*(hz_old(i,j+1) + hz_old(i,j) + dz(i,j+1) + dz(i,j))
               if(ddy_tmp > smallh) then
                  if(ffld(i,j) == 0.0d0) then ! case of NO_FRICTION_FILE_GIVEN
                     if(half*(dz(i,j+1)+dz(i,j)) > 0.0d0) then ! sea
                        cf = cfs
                     else ! land
                        cf = cfl
                     end if
                  else ! case of FRICTION_FILE_GIVEN
                     cf = ffld(i,j)
                  end if
                  if(cf > 0.0d0) then ! non-dimensional coefficent
                     bcf = cf
                  else !  Manning's roughness coefficent
                     bcf = cf*cf*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  ! explicit
!                 fric = dt*bcf*fy_old(i,j)*sqrt(fy_old(i,j)*fy_old(i,j) + fxbar*fxbar)/ddy_tmp**2
#ifdef OLD_SCHEME
                  ! semi-implicit added on 14/04/2019, Baba
                  fric = dt*bcf*half*sqrt(fy_old(i,j)*fy_old(i,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(i,j) = (fy(i,j)-fric*fy_old(i,j))/(1.0d0+fric)
#else
                  ! combinedsemi-implicit added on 14/04/2019, Baba
!                 fric = dt*bcf*half*sqrt(fy_old(i,j)*fy_old(i,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  ! simple semi-implicit added on 14/04/2019, Baba
                  fric = dt*bcf*sqrt(fy_old(i,j)*fy_old(i,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(i,j) = fy(i,j)/(1.0d0+fric)
#endif
               else
                  fy(i,j) = zap
!                  fric = zap
               end if
!               fy(i,j) = fy(i,j) - fric
! === Limiter with max Froude number. ==========================================
! === To prevent sqrt of negative numbers. =====================================
!              d = half*(dz(i,j+1)+dz(i,j)+hz_old(i,j+1)+hz_old(i,j))
               d = max(0.0d0, half*(dz(i,j+1)+dz(i,j)+hz_old(i,j+1)+hz_old(i,j)))
! ==============================================================================
               lim = froude_lim*d*sqrt(9.8d0*d)
               if(fy(i,j) > lim) then
                  fy(i,j) = lim
               else if(fy(i,j) < -lim) then
                  fy(i,j) = -lim
               end if
! ==============================================================================
#ifdef BANKFILE
            end if
         end do
      end do

      if(allocated(wfld%ir)) then
!$omp do private(theta, sint, i)
         do j = jst, jnd
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            sint = sin(theta)
            do i = ist, ind
               if(ifz(i,j)+ifz(i,j+1) > 0) then ! both wet
                  if(dz(i,j) > min_depth) then
                     if((ir(i,j+1) /= 2) .and. (ir(i,j+1) /= 3)) then
                        cv(i,j) = ((1.0d0/(rote*dth))**2/3.0d0)*(0.5d0*(dz(i,j+1)+dz(i,j)))**2
                        fy(i,j) = fy(i,j) - cv(i,j)*( &
                           ! V_n
                           &   fy_old(i  ,j+1)*sin(theta+dth)/sin(theta+dth/2.0d0) &
                           & - fy_old(i  ,j  )*sint/sin(theta-dth/2.0d0)           &
                           & - fy_old(i  ,j  )*sint/sin(theta+dth/2.0d0)           &
                           & + fy_old(i  ,j-1)*sin(theta-dth)/sin(theta-dth/2.0d0) &
                           ! U_n
                           & + (fx_old(i  ,j+1)-fx_old(i-1,j+1))/sin(theta+dth/2.0d0) &
                           & - (fx_old(i,  j  )-fx_old(i-1,j  ))/sin(theta-dth/2.0d0) &
                           & )
                     end if
                  end if
               end if
            end do
         end do
      else
!$omp do private(theta, sint, i)
      do j = jst, jnd
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint = sin(theta)
         do i = ist, ind
            if(ifz(i,j)+ifz(i,j+1) > 0) then ! both wet
#endif
! === Dispersive ===============================================================
               if(dz(i,j) > min_depth) then
                  cv(i,j) = ((1.0d0/(rote*dth))**2/3.0d0)*(0.5d0*(dz(i,j+1)+dz(i,j)))**2
                  fy(i,j) = fy(i,j) - cv(i,j)*( &
                          ! V_n
                          &   fy_old(i  ,j+1)*sin(theta+dth)/sin(theta+dth/2.0d0) &
                          & - fy_old(i  ,j  )*sint/sin(theta-dth/2.0d0)           &
                          & - fy_old(i  ,j  )*sint/sin(theta+dth/2.0d0)           &
                          & + fy_old(i  ,j-1)*sin(theta-dth)/sin(theta-dth/2.0d0) &
                          ! U_n
                          & + (fx_old(i  ,j+1)-fx_old(i-1,j+1))/sin(theta+dth/2.0d0) &
                          & - (fx_old(i,  j  )-fx_old(i-1,j  ))/sin(theta-dth/2.0d0) &
                          & )
               end if
! ==============================================================================
            end if
         end do
      end do
! ==============================================================================
#ifdef BANKFILE
      end if
#endif
#ifdef MPI
      if(iand(bflag, NORTH_BOUND) == 0) then
         j = 1
!$omp do private(dh)
         do i = ist, ind
            if(ifz(i,j) > 0 .and. ifz(i,j+1) < 0) then ! wet north, dry south
               dh = dz(i,j+1) + hz_old(i,j)
               if(dh > zap) then
                  if((-dz(i,j+1)-hz_old(i,j)) < 0.0d0) then
                     hz(i,j+1) = -dz(i,j+1)
                  end if
               end if
            end if
         end do
      end if

      if(iand(bflag, WEST_BOUND)  == 0) then
         i = 1
!$omp do private(dh)
         do j = jst, jnd
            if(ifz(i,j) > 0 .and. ifz(i+1,j) < 0) then ! wet west, dry east
               dh = dz(i+1,j) + hz_old(i,j)
               if(dh > zap) then
                  if((-dz(i+1,j)-hz_old(i,j)) < 0.0d0) then
                     hz(i+1,j) = -dz(i+1,j)
                  end if
               end if
            end if
         end do
      end if
#endif
! === Dispersive ===============================================================
#ifndef MPI
!$omp do private(i)
      do j = 0, nlat
         do i = 0, nlon
            yu(i,j) = fx(i,j)
            fx(i,j) = fx_old(i,j)
            yv(i,j) = fy(i,j)
            fy(i,j) = fy_old(i,j)
         end do
      end do
      if(gflag /= 1) then
!$omp single
         call interp2fine(VEL,cg,fg)
!$omp end single
!$omp do private(i)
         do j = 0, nlat
            do i = 0, nlon
               fx_old(i,j) = fx(i,j)
               fy_old(i,j) = fy(i,j)
            end do
         end do
      end if
#else
!$omp single
      ist_ = 1
      jst_ = 1
      if(iand(bflag, WEST_BOUND)  /= 0) ist_ = 0
      if(iand(bflag, NORTH_BOUND) /= 0) jst_ = 0
      call exchange_edges_disp_fx(fg,1)
      call exchange_edges_disp_fy(fg,1)
!$omp end single
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            yu(i,j) = fx(i,j)
            fx(i,j) = fx_old(i,j)
            yv(i,j) = fy(i,j)
            fy(i,j) = fy_old(i,j)
         end do
      end do
      if(gflag /= 1) then
!$omp single
         call interp2fine(VEL,cg,fg)
         call exchange_edges_disp_fx(fg,1)
         call exchange_edges_disp_fy(fg,1)
!$omp end single
!$omp do private(i)
         do j = jst_, nlat
            do i = ist_, nlon
               fx_old(i,j) = fx(i,j)
               fy_old(i,j) = fy(i,j)
            end do
         end do
      end if
#endif
#ifdef CONV_CHECK
!$omp do private(i, tmpx, tmpy)
      do j = 1, nlat
         do i = 1, nlon
            tmpx = 0.5d0*(dz(i+1,j)+dz(i,j)+hz_old(i+1,j)+hz_old(i,j))
            tmpy = 0.5d0*(dz(i,j+1)+dz(i,j)+hz_old(i,j+1)+hz_old(i,j))
            if(abs(tmpx) > 0.0d0) then
               tmpx = 1.0d0/tmpx
            else
               tmpx = 0.0d0
            end if
            if(abs(tmpy) > 0.0d0) then
               tmpy = 1.0d0/tmpy
            else
               tmpy = 0.0d0
            end if
            tdx(i,j) = tmpx
            tdy(i,j) = tmpy
         end do
      end do
#endif
#ifdef __FUJITSU
      allocate(ta0(ist:ind))
      allocate(ta1(ist:ind))
#endif
      do conv = 1, max_step
!$omp do private(theta, i, c0, c1)
         do j = jst, jnd
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            c0 = sin(theta+0.5d0*dth)
            c1 = sin(theta-0.5d0*dth)
            do i = ist, ind
#ifndef __FUJITSU
               fx(i,j) = (yu(i,j) + cu(i,j)*( &
                       ! V_n+1
                       & + (fy_old(i+1,j  )-fy_old(i,  j  ))*c0 &
                       & - (fy_old(i+1,j-1)-fy_old(i,  j-1))*c1 &
                       ! U_n+1
                       & + fx_old(i+1,j)+fx_old(i-1,j)))/(1.0d0+2.0d0*cu(i,j))
#else
               ta0(i) = yu(i,j) + cu(i,j)*( &
               & + (fy_old(i+1,j  )-fy_old(i,  j  ))*c0 &
               & - (fy_old(i+1,j-1)-fy_old(i,  j-1))*c1 &
               & + fx_old(i+1,j)+fx_old(i-1,j))
            end do
            do i = ist, ind
               ta1(i) = 1.0d0/(1.0d0+2.0d0*cu(i,j))
            end do
            do i = ist, ind
               fx(i,j) = ta0(i)*ta1(i)
#endif
            end do
         end do
#ifdef MPI
!$omp single
         call exchange_edges_disp_fx(fg,mod(conv,2))
!$omp end single
#endif
!$omp do private(theta, i, c0, c1, c2, c3, c4)
         do j = jst, jnd
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            c0 = 1.0d0/sin(theta+0.5d0*dth)
            c1 = 1.0d0/sin(theta-0.5d0*dth)
            c2 = sin(theta+dth)/sin(theta+0.5d0*dth)
            c3 = sin(theta-dth)/sin(theta-0.5d0*dth)
            c4 = sin(theta)*(1.0d0/sin(theta-0.5d0*dth)+1.0d0/sin(theta+0.5d0*dth))
            do i = ist, ind
#ifndef __FUJITSU
               fy(i,j) = (yv(i,j) + cv(i,j)*( &
                       ! U_n+1
                       & + (fx(i  ,j+1)-fx(i-1,j+1))*c0 &
                       & - (fx(i,  j  )-fx(i-1,j  ))*c1 &
                       ! V_n+1
                       & + fy_old(i,j+1)*c2+fy_old(i,j-1)*c3))/(1.0d0+cv(i,j)*c4)
#else
               ta0(i) = yv(i,j) + cv(i,j)*( &
               & + (fx(i  ,j+1)-fx(i-1,j+1))*c0 &
               & - (fx(i,  j  )-fx(i-1,j  ))*c1 &
               & + fy_old(i,j+1)*c2+fy_old(i,j-1)*c3)
            end do
            do i = ist, ind
               ta1(i) = 1.0d0/(1.0d0+cv(i,j)*c4)
            end do
            do i = ist, ind
               fy(i,j) = ta0(i)*ta1(i)
#endif
            end do
         end do
#ifdef MPI
!$omp single
         call exchange_edges_disp_fy(fg,mod(conv,2))
!$omp end single
#endif
#ifdef CONV_CHECK
#ifdef LESS_CC
         if(mod(conv,check) == 0) then
#endif
!$omp single
         diffxy = 0.0d0
!$omp end single
!$omp do private(i, diff, diffx, diffy) reduction(max:diffxy)
         do j = jst, jnd
            do i = ist, ind
               diffx = abs((fx(i,j)-fx_old(i,j))*tdx(i,j))
               diffy = abs((fy(i,j)-fy_old(i,j))*tdy(i,j))
               diff = max(diffx,diffy)
               diffxy = max(diffxy,diff)
            end do
         end do
#ifdef MPI
!$omp single
#ifndef MULTI
         call MPI_Allreduce(MPI_IN_PLACE, diffxy, 1, REAL_MPI, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
         call MPI_Allreduce(MPI_IN_PLACE, diffxy, 1, REAL_MPI, MPI_MAX, MPI_MEMBER_WORLD, ierr)
#endif
!$omp end single
#endif
         if(diffxy < conv_val) then
            exit
         endif
#ifdef LESS_CC
         end if
#endif
#endif
!$omp single
         p => fx_old
         fx_old => fx
         fx => p
         p => fy_old
         fy_old => fy
         fy => p
!$omp end single
      end do
#ifdef CONV_CHECK
!$omp single
      conv_step = conv
!$omp end single
#endif
#ifndef MPI
!$omp do private(i)
      do j = 0, nlat
         do i = 0, nlon
            fx_old(i,j) = fx(i,j)
            fy_old(i,j) = fy(i,j)
         end do
      end do
#else
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            fx_old(i,j) = fx(i,j)
            fy_old(i,j) = fy(i,j)
         end do
      end do
#endif
!$omp single
      fx => wfld%fx
      fy => wfld%fy
!$omp end single
#ifdef __FUJITSU
      deallocate(ta0)
      deallocate(ta1)
#endif
! ==============================================================================

      if(gflag == 1) then
#ifdef OLD_SCHEME
         ! do linear calc along edges (i=0,nlon-2;j=0,nlat-2) of grid
#else
         ! do non-linear calc without advection term along edges (i=0,nlon-2;j=0,nlat-2) of grid
#endif
         ! note that fx,fy updates are moved in 1 grid along maximal edges

         ! north
#ifndef MPI
!$omp single
         sint = sin(th0)
#ifndef OLD_SCHEME
         crls = 2.0d0*cflag*cfreq*cos(th0)*dt
#endif
#else
         if(iand(bflag, NORTH_BOUND) /= 0) then
!$omp single
            ist = 1
            ind = nlon - 1
            if(iand(bflag, EAST_BOUND) /= 0) ind = nlon - 1
            sint = sin(th0+joff*dth)
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(th0+joff*dth)*dt
#endif
#endif
         gdtdss = gdtds/sint
!$omp end single
#ifdef OLD_SCHEME
!$omp do
#else
!$omp do  private(fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
#ifndef MPI
         do i = 1, nlon-1
#else
         do i = ist, ind
#endif
#ifdef OLD_SCHEME
            if(dz(i,1) > zap) then
               fx(i,1) = fx(i,1) - (0.5d0*(dz(i+1,1)+dz(i,1)))*gdtdss*(hz_old(i+1,1)-hz_old(i,1))
               fy(i,1) = fy(i,1) - (0.5d0*(dz(i,  2)+dz(i,1)))*gdtds *(hz_old(i,  2)-hz_old(i,1))
#else
            ddx_tmp = half*(hz_old(i+1,1) + hz_old(i,1) + dz(i+1,1) + dz(i,1))
            if(dz(i,1) > zap .AND. dz(i+1,1) > zap .AND. ddx_tmp > zap) then
               fybar = half*(fy_old(i,1) + fy_old(i+1,1))
               fx(i,1) = fx_old(i,1) - ddx_tmp*gdtdss*(hz_old(i+1,1)-hz_old(i,1))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(i,1) == 0.0d0) then 
                     if(half*(dz(i+1,1)+dz(i,1)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(i,1)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(i,1)*fx_old(i,1)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(i,1) = fx(i,1)/(1.0d0+fric)
               end if
#endif
            else
               fx(i,1) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(i,2) + hz_old(i,1) + dz(i,2) + dz(i,1))
            if(dz(i,1) > zap .AND. dz(i,2) > zap .AND. ddy_tmp > zap) then
               if (i > 1) then
                  fxbar = quart*(fx_old(i,1) + fx_old(i-1,1) + fx_old(i-1,2) + fx_old(i,2))
               else
                  fxbar = half *(fx_old(i,1) + fx_old(i,2))
               end if
               fy(i,1) = fy_old(i,1) - ddy_tmp*gdtds*(hz_old(i,2)-hz_old(i,1))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(i,1) == 0.0d0) then 
                     if(half*(dz(i,2)+dz(i,1)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(i,1)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(i,1)*fy_old(i,1)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(i,1) = fy(i,1)/(1.0d0+fric2)
               end if
            else
#endif
               fy(i,1) = zap
            end if
         end do
#ifdef MPI
         end if
#endif
#ifndef MPI
!$omp single
         sint = sin(th0 + dth)
#ifndef OLD_SCHEME
         crls = 2.0d0*cflag*cfreq*cos(th0 + dth)*dt
#endif
#else
         if(iand(bflag, NORTH_BOUND) /= 0) then
!$omp single
            ist = 1
            ind = nlon - 1
            if(iand(bflag, EAST_BOUND) /= 0) ind = nlon - 1
            sint = sin(th0+(joff+1)*dth)
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(th0+(joff+1)*dth)*dt
#endif
#endif
         gdtdss = gdtds/sint
!$omp end single
#ifdef OLD_SCHEME
!$omp do
#else
!$omp do  private(fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
#ifndef MPI
         do i = 1, nlon-1
#else
         do i = ist, ind
#endif
#ifdef OLD_SCHEME
            if(dz(i,2) > zap) then
               fx(i,2) = fx(i,2) - (0.5d0*(dz(i+1,2)+dz(i,2)))*gdtdss*(hz_old(i+1,2)-hz_old(i,2))
               fy(i,2) = fy(i,2) - (0.5d0*(dz(i,  3)+dz(i,2)))*gdtds *(hz_old(i,  3)-hz_old(i,2))
#else
            ddx_tmp = half*(hz_old(i+1,2) + hz_old(i,2) + dz(i+1,2) + dz(i,2))
            if(dz(i,2) > zap .AND. dz(i+1,2) > zap .AND. ddx_tmp > zap) then
               fybar = quart*(fy_old(i,2) + fy_old(i+1,2) + fy_old(i,1) + fy_old(i+1,1))
               fx(i,2) = fx_old(i,2) - ddx_tmp*gdtdss*(hz_old(i+1,2)-hz_old(i,2))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(i,2) == 0.0d0) then 
                     if(half*(dz(i+1,2)+dz(i,2)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(i,2)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(i,2)*fx_old(i,2)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(i,2) = fx(i,2)/(1.0d0+fric)
               end if
#endif
            else
               fx(i,2) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(i,3) + hz_old(i,2) + dz(i,3) + dz(i,2))
            if(dz(i,2) > zap .AND. dz(i,3) > zap .AND. ddy_tmp > zap) then
               if (i > 1) then
                  fxbar = quart*(fx_old(i,2) + fx_old(i-1,2) + fx_old(i-1,3) + fx_old(i,3))
               else
                  fxbar =  half*(fx_old(i,2) + fx_old(i,3))
               end if
               fy(i,2) = fy_old(i,2) - ddy_tmp*gdtds*(hz_old(i,3)-hz_old(i,2))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(i,2) == 0.0d0) then 
                     if(half*(dz(i,3)+dz(i,2)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(i,2)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(i,2)*fy_old(i,2)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(i,2) = fy(i,2)/(1.0d0+fric2)
               end if
            else
#endif
               fy(i,2) = zap
            end if
         end do
#ifdef MPI
         end if
#endif
#ifdef UPWIND3
#ifndef MPI
!$omp single
         sint = sin(th0 + 2*dth)
#ifndef OLD_SCHEME
         crls = 2.0d0*cflag*cfreq*cos(th0 + 2*dth)*dt
#endif
#else
         if(iand(bflag, NORTH_BOUND) /= 0) then
!$omp single
            ist = 1
            ind = nlon - 1
            if(iand(bflag, EAST_BOUND) /= 0) ind = nlon - 1
            sint = sin(th0+(joff+2)*dth)
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(th0+(joff+2)*dth)*dt
#endif
#endif
         gdtdss = gdtds/sint
!$omp end single
#ifdef OLD_SCHEME
!$omp do
#else
!$omp do  private(fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
#ifndef MPI
         do i = 1, nlon-1
#else
         do i = ist, ind
#endif
#ifdef OLD_SCHEME
            if(dz(i,3) > zap) then
               fx(i,3) = fx(i,3) - (0.5d0*(dz(i+1,3)+dz(i,3)))*gdtdss*(hz_old(i+1,3)-hz_old(i,3))
               fy(i,3) = fy(i,3) - (0.5d0*(dz(i,  4)+dz(i,3)))*gdtds *(hz_old(i,  4)-hz_old(i,3))
#else
            ddx_tmp = half*(hz_old(i+1,3) + hz_old(i,3) + dz(i+1,3) + dz(i,3))
            if(dz(i,3) > zap .AND. dz(i+1,3) > zap .AND. ddx_tmp > zap) then
               fybar = quart*(fy_old(i,3) + fy_old(i+1,3) + fy_old(i,2) + fy_old(i+1,2))
               fx(i,3) = fx_old(i,3) - ddx_tmp*gdtdss*(hz_old(i+1,3)-hz_old(i,3))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(i,3) == 0.0d0) then 
                     if(half*(dz(i+1,3)+dz(i,3)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(i,3)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(i,2)*fx_old(i,3)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(i,3) = fx(i,3)/(1.0d0+fric)
               end if
            else
            fx(i,3) = zap
            end if
            ddy_tmp = half*(hz_old(i,4) + hz_old(i,3) + dz(i,4) + dz(i,3))
            if(dz(i,3) > zap .AND. dz(i,4) > zap .AND. ddy_tmp > zap) then
               if (i > 1) then
                  fxbar = quart*(fx_old(i,3) + fx_old(i-1,3) + fx_old(i-1,4) + fx_old(i,4))
               else
                  fxbar = half*(fx_old(i,3) + fx_old(i,4))
               end if
               fy(i,3) = fy_old(i,3) - ddy_tmp*gdtds*(hz_old(i,4)-hz_old(i,3))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(i,3) == 0.0d0) then 
                     if(half*(dz(i,4)+dz(i,3)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(i,3)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(i,3)*fy_old(i,3)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(i,3) = fy(i,3)/(1.0d0+fric2)
               end if
#endif
            else
#ifdef OLD_SCHEME
               fx(i,3) = zap
#endif
               fy(i,3) = zap
            end if
         end do
#ifdef MPI
         end if
#endif
#endif

         ! south
#ifndef MPI
!$omp single
         sint = sin(th0 + (nlat-2)*dth)
#ifndef OLD_SCHEME
         crls = 2.0d0*cflag*cfreq*cos(th0 + (nlat-2)*dth)*dt
#endif
#else
         if(iand(bflag, SOUTH_BOUND) /= 0) then
!$omp single
            ist = 1
            ind = nlon - 1
            if(iand(bflag, WEST_BOUND) /= 0) ist = 1
            if(iand(bflag, EAST_BOUND) /= 0) ind = nlon - 1
            sint = sin(th0 + (nlat-2+joff)*dth)
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(th0 + (nlat-2+joff)*dth)*dt
#endif
#endif
         gdtdss = gdtds/sint
!$omp end single
#ifdef OLD_SCHEME
!$omp do
#else
!$omp do  private(fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
#ifndef MPI
         do i = 1, nlon-1
#else
         do i = ist, ind
#endif
#ifdef OLD_SCHEME
            if(dz(i,nlat-1) > zap) then
               fx(i,nlat-1) = fx(i,nlat-1) - (0.5d0*(dz(i+1,nlat-1)+dz(i,nlat-1)))*gdtdss*(hz_old(i+1,nlat-1)-hz_old(i,nlat-1))
               fy(i,nlat-1) = fy(i,nlat-1) - (0.5d0*(dz(i,  nlat)  +dz(i,nlat-1)))*gdtds *(hz_old(i,  nlat)  -hz_old(i,nlat-1))
#else
            ddx_tmp = half*(hz_old(i+1,nlat-1) + hz_old(i,nlat-1) + dz(i+1,nlat-1) + dz(i,nlat-1))
            if(dz(i,nlat-1) > zap .AND. dz(i+1,nlat-1) > zap .AND. ddx_tmp > zap) then
               fybar = quart*(fy_old(i,nlat-1) + fy_old(i+1,nlat-1) + fy_old(i,nlat-2) + fy_old(i+1,nlat-2))
               fx(i,nlat-1) = fx_old(i,nlat-1) - ddx_tmp*gdtdss*(hz_old(i+1,nlat-1)-hz_old(i,nlat-1))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(i,nlat-1) == 0.0d0) then 
                     if(half*(dz(i+1,nlat-1)+dz(i,nlat-1)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(i,nlat-1)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(i,nlat-1)*fx_old(i,nlat-1)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(i,nlat-1) = fx(i,nlat-1)/(1.0d0+fric)
               end if
#endif
            else
               fx(i,nlat-1) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(i,nlat) + hz_old(i,nlat-1) + dz(i,nlat) + dz(i,nlat-1))
            if(dz(i,nlat-1) > zap .AND. dz(i,nlat) > zap .AND. ddy_tmp > zap) then
               if (i > 1) then
                  fxbar = quart*(fx_old(i,nlat-1) + fx_old(i-1,nlat-1) + fx_old(i-1,nlat) + fx_old(i,nlat))
               else
                  fxbar = half *(fx_old(i,nlat-1) + fx_old(i,nlat))
               end if
               fy(i,nlat-1) = fy_old(i,nlat-1) - ddy_tmp*gdtds *(hz_old(i,nlat)-hz_old(i,nlat-1))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(i,nlat-1) == 0.0d0) then 
                     if(half*(dz(i,nlat)+dz(i,nlat-1)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(i,nlat-1)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(i,nlat-1)*fy_old(i,nlat-1)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(i,nlat-1) = fy(i,nlat-1)/(1.0d0+fric2)
               end if
            else
#endif
               fy(i,nlat-1) = zap
            end if
         end do
#ifdef MPI
         end if
#endif
#ifdef UPWIND3
#ifndef MPI
!$omp single
         sint = sin(th0 + (nlat-3)*dth)
#ifndef OLD_SCHEME
         crls = 2.0d0*cflag*cfreq*cos(th0 + (nlat-3)*dth)*dt
#endif
#else
         if(iand(bflag, SOUTH_BOUND) /= 0) then
!$omp single
            ist = 1
            ind = nlon - 1
            if(iand(bflag, WEST_BOUND) /= 0) ist = 1
            if(iand(bflag, EAST_BOUND) /= 0) ind = nlon - 1
            sint = sin(th0 + (nlat-3+joff)*dth)
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(th0 + (nlat-3+joff)*dth)*dt
#endif
#endif
         gdtdss = gdtds/sint
!$omp end single
#ifdef OLD_SCHEME
!$omp do
#else
!$omp do  private(fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
#ifndef MPI
         do i = 1, nlon-1
#else
         do i = ist, ind
#endif
#ifdef OLD_SCHEME
            if(dz(i,nlat-2) > zap) then
               fx(i,nlat-2) = fx(i,nlat-2) - (0.5d0*(dz(i+1,nlat-2)+dz(i,nlat-2)))*gdtdss*(hz_old(i+1,nlat-2)-hz_old(i,nlat-2))
               fy(i,nlat-2) = fy(i,nlat-2) - (0.5d0*(dz(i,  nlat-1)+dz(i,nlat-2)))*gdtds *(hz_old(i,  nlat-1)-hz_old(i,nlat-2))
#else
            ddx_tmp = half*(hz_old(i+1,nlat-2) + hz_old(i,nlat-2) + dz(i+1,nlat-2) + dz(i,nlat-2))
            if(dz(i,nlat-2) > zap .AND. dz(i+1,nlat-2) > zap .AND. ddx_tmp > zap) then
               fybar = quart*(fy_old(i,nlat-2) + fy_old(i+1,nlat-2) + fy_old(i,nlat-3) + fy_old(i+1,nlat-3))
               fx(i,nlat-2) = fx_old(i,nlat-2) - ddx_tmp*gdtdss*(hz_old(i+1,nlat-2)-hz_old(i,nlat-2))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(i,nlat-2) == 0.0d0) then 
                     if(half*(dz(i+1,nlat-2)+dz(i,nlat-2)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(i,nlat-2)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(i,nlat-2)*fx_old(i,nlat-2)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(i,nlat-2) = fx(i,nlat-2)/(1.0d0+fric)
               end if
#endif
            else
               fx(i,nlat-2) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(i,nlat-1) + hz_old(i,nlat-2) + dz(i,nlat-1) + dz(i,nlat-2))
            if(dz(i,nlat-2) > zap .AND. dz(i,nlat-1) > zap .AND. ddy_tmp > zap) then
               if (i > 1) then
                  fxbar = quart*(fx_old(i,nlat-2) + fx_old(i-1,nlat-2) + fx_old(i-1,nlat-1) + fx_old(i,nlat-1))
               else
                  fxbar = half *(fx_old(i,nlat-2) + fx_old(i,nlat-1))
               end if
               fy(i,nlat-2) = fy_old(i,nlat-2) - ddy_tmp*gdtds *(hz_old(i,nlat-1)-hz_old(i,nlat-2))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(i,nlat-2) == 0.0d0) then 
                     if(half*(dz(i,nlat)+dz(i,nlat-2)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(i,nlat-2)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(i,nlat-2)*fy_old(i,nlat-2)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(i,nlat-2) = fy(i,nlat-2)/(1.0d0+fric2)
               end if
            else
#endif
               fy(i,nlat-2) = zap
            end if
         end do
#ifdef MPI
         end if
#endif
#endif

#ifndef MPI
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
#ifndef UPWIND3
         do j = 2, nlat-2
#else
         do j = 4, nlat-3
#endif
            theta = th0 + (j-1)*dth
#else
         if(iand(bflag, EAST_BOUND) /= 0) then
!$omp single
            jst = 1
            jnd = nlat - 1
#ifndef UPWIND3
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 2
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 2
#else
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 4
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 3
#endif
!$omp end single
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
            do j = jst, jnd
               theta = th0 + (j+joff-1)*dth
#endif
            sint = sin(theta)
            gdtdss = gdtds/sint
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(theta)*dt
#endif

            ! east
#ifdef OLD_SCHEME
            if(dz(nlon-1,j) > zap) then
               fx(nlon-1,j) = fx(nlon-1,j) - (0.5d0*(dz(nlon,  j)  +dz(nlon-1,j)))*gdtdss*(hz_old(nlon,  j)  -hz_old(nlon-1,j))
               fy(nlon-1,j) = fy(nlon-1,j) - (0.5d0*(dz(nlon-1,j+1)+dz(nlon-1,j)))*gdtds *(hz_old(nlon-1,j+1)-hz_old(nlon-1,j))
#else
            ddx_tmp = half*(hz_old(nlon,j) + hz_old(nlon-1,j) + dz(nlon,j) + dz(nlon-1,j))
            if(dz(nlon-1,j) > zap .AND. dz(nlon,j) > zap .AND. ddx_tmp > zap) then
               if (j > 1) then
                  fybar = quart*(fy_old(nlon-1,j) + fy_old(nlon,j) + fy_old(nlon-1,j-1) + fy_old(nlon,j-1))
               else
                  fybar = half*(fy_old(nlon-1,j) + fy_old(nlon,j))
               end if
               fx(nlon-1,j) = fx_old(nlon-1,j) - ddx_tmp*gdtdss*(hz_old(nlon,j)-hz_old(nlon-1,j))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(nlon-1,j) == 0.0d0) then 
                     if(half*(dz(nlon,j)+dz(nlon-1,j)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(nlon-1,j)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(nlon-1,j)*fx_old(nlon-1,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(nlon-1,j) = fx(nlon-1,j)/(1.0d0+fric)
               end if
#endif
            else
               fx(nlon-1,j) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(nlon-1,j+1) + hz_old(nlon-1,j) + dz(nlon-1,j+1) + dz(nlon-1,j))
            if(dz(nlon-1,j) > zap .AND. dz(nlon-1,j+1) > zap .AND. ddy_tmp > zap) then
               fxbar = quart*(fx_old(nlon-1,j) + fx_old(nlon-2,j) + fx_old(nlon-2,j+1) + fx_old(nlon-1,j+1))
               fy(nlon-1,j) = fy_old(nlon-1,j) - ddy_tmp*gdtds *(hz_old(nlon-1,j+1)-hz_old(nlon-1,j))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(nlon-1,j) == 0.0d0) then 
                     if(half*(dz(nlon-1,j+1)+dz(nlon-1,j)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(nlon-1,j)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(nlon-1,j)*fy_old(nlon-1,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(nlon-1,j) = fy(nlon-1,j)/(1.0d0+fric2)
               end if
            else
#endif
               fy(nlon-1,j) = zap
            end if
#ifdef MPI
            end do
         end if

         if(iand(bflag, WEST_BOUND) /= 0) then
!$omp single
            jst = 1
            jnd = nlat - 1
#ifndef UPWIND3
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 2
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 2
#else
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 4
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 3
#endif
!$omp end single
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
            do j = jst, jnd
               theta = th0 + (j+joff-1)*dth
               sint = sin(theta)
               gdtdss = gdtds/sint
#ifndef OLD_SCHEME
               crls = 2.0d0*cflag*cfreq*cos(theta)*dt
#endif
#endif

            ! west
#ifdef OLD_SCHEME
            if(dz(1,j) > zap) then
               fx(1,j) = fx(1,j) - (0.5d0*(dz(2,j)  +dz(1,j)))*gdtdss*(hz_old(2,j)  -hz_old(1,j))
               fy(1,j) = fy(1,j) - (0.5d0*(dz(1,j+1)+dz(1,j)))*gdtds *(hz_old(1,j+1)-hz_old(1,j))
#else
            ddx_tmp = half*(hz_old(2,j) + hz_old(1,j) + dz(2,j) + dz(1,j))
            if(dz(1,j) > zap .AND. dz(2,j) > zap .AND. ddx_tmp > zap) then
               if (j > 1) then
                  fybar = quart*(fy_old(1,j) + fy_old(2,j) + fy_old(1,j-1) + fy_old(2,j-1))
               else
                  fybar = half*(fy_old(1,j) + fy_old(2,j))
               end if
               fx(1,j) = fx_old(1,j) - ddx_tmp*gdtdss*(hz_old(2,j)-hz_old(1,j))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(1,j) == 0.0d0) then 
                     if(half*(dz(2,j)+dz(1,j)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(1,j)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(1,j)*fx_old(1,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(1,j) = fx(1,j)/(1.0d0+fric)
               end if
#endif
            else
               fx(1,j) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(1,j+1) + hz_old(1,j) + dz(1,j+1) + dz(1,j))
            if(dz(1,j) > zap .AND. dz(1,j+1) > zap .AND. ddy_tmp > zap) then
               fxbar = half*(fx_old(1,j) + fx_old(1,j+1))
               fy(1,j) = fy_old(1,j) - ddy_tmp*gdtds *(hz_old(1,j+1)-hz_old(1,j))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(1,j) == 0.0d0) then 
                     if(half*(dz(1,j+1)+dz(1,j)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(1,j)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(1,j)*fy_old(1,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(1,j) = fy(1,j)/(1.0d0+fric2)
               end if
            else
#endif
               fy(1,j) = zap
            end if
         end do
#ifdef MPI
         end if
#endif
#ifndef MPI
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
#ifndef UPWIND3
         do j = 2, nlat-2
#else
         do j = 4, nlat-3
#endif
            theta = th0 + (j-1)*dth
            sint = sin(theta)
            gdtdss = gdtds/sint
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(theta)*dt
#endif
            ! west
#ifdef OLD_SCHEME
            if(dz(2,j) > zap) then
               fx(2,j) = fx(2,j) - (0.5d0*(dz(3,j)  +dz(2,j)))*gdtdss*(hz_old(3,j)  -hz_old(2,j))
               fy(2,j) = fy(2,j) - (0.5d0*(dz(2,j+1)+dz(2,j)))*gdtds *(hz_old(2,j+1)-hz_old(2,j))
#else
            ddx_tmp = half*(hz_old(3,j) + hz_old(2,j) + dz(3,j) + dz(2,j))
            if(dz(2,j) > zap .AND. dz(3,j) > zap .AND. ddx_tmp > zap) then
               fybar = quart*(fy_old(2,j) + fy_old(3,j) + fy_old(2,j-1) + fy_old(3,j-1))
               fx(2,j) = fx_old(2,j) - ddx_tmp*gdtdss*(hz_old(3,j)-hz_old(2,j))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(2,j) == 0.0d0) then 
                     if(half*(dz(3,j)+dz(2,j)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(2,j)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(2,j)*fx_old(2,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(2,j) = fx(2,j)/(1.0d0+fric)
               end if
#endif
            else
               fx(2,j) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(2,j+1) + hz_old(2,j) + dz(2,j+1) + dz(2,j))
            if(dz(2,j) > zap .AND. dz(2,j+1) > zap .AND. ddy_tmp > zap) then
               fxbar = quart*(fx_old(2,j) + fx_old(1,j) + fx_old(1,j+1) + fx_old(2,j+1))
               fy(2,j) = fy_old(2,j) - ddy_tmp*gdtds *(hz_old(2,j+1)-hz_old(2,j))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(2,j) == 0.0d0) then 
                     if(half*(dz(2,j+1)+dz(2,j)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(2,j)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(2,j)*fy_old(2,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(2,j) = fy(2,j)/(1.0d0+fric2)
               end if
            else
#endif
               fy(2,j) = zap
            end if
         end do

#else
         if(iand(bflag, WEST_BOUND) /= 0) then
!$omp single
            jst = 1
            jnd = nlat - 1
#ifndef UPWIND3
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 2
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 2
#else
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 4
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 3
#endif
!$omp end single
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
            do j = jst, jnd
               theta = th0 + (j+joff-1)*dth
               sint = sin(theta)
               gdtdss = gdtds/sint
#ifndef OLD_SCHEME
               crls = 2.0d0*cflag*cfreq*cos(theta)*dt
#endif
               ! west
#ifdef OLD_SCHEME
               if(dz(2,j) > zap) then
                  fx(2,j) = fx(2,j) - (0.5d0*(dz(3,j)  +dz(2,j)))*gdtdss*(hz_old(3,j)  -hz_old(2,j))
                  fy(2,j) = fy(2,j) - (0.5d0*(dz(2,j+1)+dz(2,j)))*gdtds *(hz_old(2,j+1)-hz_old(2,j))
#else
            ddx_tmp = half*(hz_old(3,j) + hz_old(2,j) + dz(3,j) + dz(2,j))
            if(dz(2,j) > zap .AND. dz(3,j) > zap .AND. ddx_tmp > zap) then
               if (j > 1) then
                  fybar = quart*(fy_old(2,j) + fy_old(3,j) + fy_old(2,j-1) + fy_old(3,j-1))
#endif
               else
#ifdef OLD_SCHEME
                  fx(2,j) = zap
                  fy(2,j) = zap
#else
                  fybar = half*(fy_old(2,j) + fy_old(3,j))
               end if
               fx(2,j) = fx_old(2,j) - ddx_tmp*gdtdss*(hz_old(3,j)-hz_old(2,j))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(2,j) == 0.0d0) then 
                     if(half*(dz(3,j)+dz(2,j)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(2,j)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(2,j)*fx_old(2,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(2,j) = fx(2,j)/(1.0d0+fric)
               end if
            else
               fx(2,j) = zap
            end if
            ddy_tmp = half*(hz_old(2,j+1) + hz_old(2,j) + dz(2,j+1) + dz(2,j))
            if(dz(2,j) > zap .AND. dz(2,j+1) > zap .AND. ddy_tmp > zap) then
               fxbar = quart*(fx_old(2,j) + fx_old(1,j) + fx_old(1,j+1) + fx_old(2,j+1))
               fy(2,j) = fy_old(2,j) - ddy_tmp*gdtds *(hz_old(2,j+1)-hz_old(2,j))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(2,j) == 0.0d0) then 
                     if(half*(dz(2,j+1)+dz(2,j)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(2,j)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(2,j)*fy_old(2,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(2,j) = fy(2,j)/(1.0d0+fric2)
#endif
               end if
#ifndef OLD_SCHEME
            else
               fy(2,j) = zap
            end if
#endif
            end do
         end if
#endif
#ifdef UPWIND3
#ifndef MPI
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
         do j = 4, nlat-3
            theta = th0 + (j-1)*dth
#else
         if(iand(bflag, EAST_BOUND) /= 0) then
!$omp single
            jst = 1
            jnd = nlat - 1
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 4
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 3
!$omp end single
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
            do j = jst, jnd
               theta = th0 + (j+joff-1)*dth
#endif
            sint = sin(theta)
            gdtdss = gdtds/sint
#ifndef OLD_SCHEME
            crls = 2.0d0*cflag*cfreq*cos(theta)*dt
#endif

            ! east
#ifdef OLD_SCHEME
            if(dz(nlon-2,j) > zap) then
               fx(nlon-2,j) = fx(nlon-2,j) - (0.5d0*(dz(nlon-1,j)  +dz(nlon-2,j)))*gdtdss*(hz_old(nlon-1,j)  -hz_old(nlon-2,j))
               fy(nlon-2,j) = fy(nlon-2,j) - (0.5d0*(dz(nlon-2,j+1)+dz(nlon-2,j)))*gdtds *(hz_old(nlon-2,j+1)-hz_old(nlon-2,j))
#else
            ddx_tmp = half*(hz_old(nlon-1,j) + hz_old(nlon-2,j) + dz(nlon-1,j) + dz(nlon-2,j))
            if(dz(nlon-1,j) > zap .AND. dz(nlon-2,j) > zap .AND. ddx_tmp > zap) then
               if (j > 1) then
                  fybar = quart*(fy_old(nlon-2,j) + fy_old(nlon-1,j) + fy_old(nlon-2,j-1) + fy_old(nlon-1,j-1))
               else
                  fybar = half*(fy_old(nlon-2,j) + fy_old(nlon-1,j))
               end if
               fx(nlon-2,j) = fx_old(nlon-2,j) - ddx_tmp*gdtdss*(hz_old(nlon-2,j)-hz_old(nlon-1,j))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(nlon-2,j) == 0.0d0) then 
                     if(half*(dz(nlon,j)+dz(nlon-2,j)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(nlon-2,j)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(nlon-2,j)*fx_old(nlon-2,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(nlon-2,j) = fx(nlon-2,j)/(1.0d0+fric)
               end if
#endif
            else
               fx(nlon-2,j) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(nlon-2,j+1) + hz_old(nlon-2,j) + dz(nlon-2,j+1) + dz(nlon-2,j))
            if(dz(nlon-2,j) > zap .AND. dz(nlon-2,j+1) > zap .AND. ddy_tmp > zap) then
               fxbar = quart*(fx_old(nlon-2,j) + fx_old(nlon-3,j) + fx_old(nlon-3,j+1) + fx_old(nlon-2,j+1))
               fy(nlon-2,j) = fy_old(nlon-2,j) - ddy_tmp*gdtds *(hz_old(nlon-2,j+1)-hz_old(nlon-2,j))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(nlon-2,j) == 0.0d0) then 
                     if(half*(dz(nlon-2,j+1)+dz(nlon-2,j)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(nlon-2,j)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(nlon-2,j)*fy_old(nlon-2,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(nlon-2,j) = fy(nlon-2,j)/(1.0d0+fric2)
               end if
            else
#endif
               fy(nlon-2,j) = zap
            end if
#ifdef MPI
            end do
         end if

         if(iand(bflag, WEST_BOUND) /= 0) then
!$omp single
            jst = 1
            jnd = nlat - 1
            if(iand(bflag, NORTH_BOUND) /= 0) jst = 4
            if(iand(bflag, SOUTH_BOUND) /= 0) jnd = nlat - 3
!$omp end single
#ifdef OLD_SCHEME
!$omp do private(theta, sint, gdtdss)
#else
!$omp do private(theta, sint, gdtdss, crls, fxbar, fybar, cf, cf2, bcf, bcf2, fric, fric2, ddx_tmp, ddy_tmp)
#endif
            do j = jst, jnd
               theta = th0 + (j+joff-1)*dth
               sint = sin(theta)
               gdtdss = gdtds/sint
#ifndef OLD_SCHEME
               crls = 2.0d0*cflag*cfreq*cos(theta)*dt
#endif
#endif

            ! west
#ifdef OLD_SCHEME
            if(dz(3,j) > zap) then
               fx(3,j) = fx(3,j) - (0.5d0*(dz(4,j)  +dz(3,j)))*gdtdss*(hz_old(4,j)  -hz_old(3,j))
               fy(3,j) = fy(3,j) - (0.5d0*(dz(3,j+1)+dz(3,j)))*gdtds *(hz_old(3,j+1)-hz_old(3,j))
#else
            ddx_tmp = half*(hz_old(4,j) + hz_old(3,j) + dz(4,j) + dz(3,j))
            if(dz(3,j) > zap .AND. dz(4,j) > zap .AND. ddx_tmp > zap) then
               if (j > 1) then
                  fybar = quart*(fy_old(3,j) + fy_old(4,j) + fy_old(3,j-1) + fy_old(4,j-1))
               else
                  fybar = half*(fy_old(3,j) + fy_old(4,j))
               end if
               fx(3,j) = fx_old(3,j) - ddx_tmp*gdtdss*(hz_old(4,j)-hz_old(3,j))-crls*fybar
               if(ddx_tmp < 100.0d0) then
                  if(ffld(3,j) == 0.0d0) then 
                     if(half*(dz(4,j)+dz(3,j)) > 0.0d0) then
                        cf = cfs
                     else
                        cf = cfl
                     end if
                  else 
                     cf = ffld(3,j)
                  end if
                  if(cf > 0.0d0) then 
                     bcf = cf
                  else  
                     bcf = cf*cf*9.8d0*ddx_tmp**(-1.0d0/3.0d0)
                  end if
                  fric = dt*bcf*sqrt(fx_old(3,j)*fx_old(3,j)+fybar*fybar)/ddx_tmp/ddx_tmp
                  fx(3,j) = fx(3,j)/(1.0d0+fric)
               end if
#endif
            else
               fx(3,j) = zap
#ifndef OLD_SCHEME
            end if
            ddy_tmp = half*(hz_old(3,j+1) + hz_old(3,j) + dz(3,j+1) + dz(3,j))
            if(dz(3,j) > zap .AND. dz(3,j+1) > zap .AND. ddy_tmp > zap) then
               fxbar = quart*(fx_old(3,j) + fx_old(2,j) + fx_old(2,j+1) + fx_old(3,j+1))
               fy(3,j) = fy_old(3,j) - ddy_tmp*gdtds *(hz_old(3,j+1)-hz_old(3,j))+crls*fxbar
               if(ddy_tmp < 100.0d0) then
                  if(ffld(3,j) == 0.0d0) then 
                     if(half*(dz(3,j+1)+dz(3,j)) > 0.0d0) then 
                        cf2 = cfs
                     else 
                        cf2 = cfl
                     end if
                  else 
                     cf2 = ffld(3,j)
                  end if
                  if(cf2 > 0.0d0) then
                     bcf2 = cf2
                  else 
                     bcf2 = cf2*cf2*9.8d0*ddy_tmp**(-1.0d0/3.0d0)
                  end if
                  fric2 = dt*bcf2*sqrt(fy_old(3,j)*fy_old(3,j)+fxbar*fxbar)/ddy_tmp/ddy_tmp
                  fy(3,j) = fy(3,j)/(1.0d0+fric2)
               end if
            else
#endif
               fy(3,j) = zap
            end if
         end do
#ifdef MPI
         end if
#endif
#endif
      end if
!$omp end parallel

      return
   end subroutine fxynl_rwg_disp

end module mod_fxy_disp
