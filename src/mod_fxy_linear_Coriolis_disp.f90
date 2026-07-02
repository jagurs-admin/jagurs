#ifdef TIMER
#   define TIMER_START(a) call start_timer(a)
#   define TIMER_STOP(a)  call stop_timer(a)
#else
#   define TIMER_START(a)
#   define TIMER_STOP(a)
#endif
#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_fxy_Coriolis_disp
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
use mod_params, only : timenest
#ifdef MPI
use mod_params, only : VEL
#endif
#ifdef TIMER
use mod_timer
#endif
implicit none

contains

! === Coriolis force is supported on linear calc. ==============================
!  subroutine fxy_rwg_disp(wfld,dfld,dt,th0,dth,joff,nlon,nlat,gflag,bflag,fg,cg,conv_step)
   subroutine fxy_rwg_Coriolis_disp(wfld,dfld,ifz,cflag,dt,th0,dth,joff,nlon,nlat,gflag,bflag,fg,cg,conv_step)
! ==============================================================================
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(in) :: ifz
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(in) :: ifz
#endif
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
      integer(kind=4), intent(in) :: joff
      integer(kind=4), intent(out) :: conv_step
! === Coriolis force is supported on linear calc. ==============================
      integer(kind=4), intent(in) :: cflag
! ==============================================================================

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz, p
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: rote, g, gts, theta, sint, cfc
      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
! === Coriolis force is supported on linear calc. ==============================
      real(kind=REAL_BYTE), parameter :: quart = 0.25d0
      real(kind=REAL_BYTE), parameter :: cfreq = 7.27220522d-05 ! angular freq of earth's rotation (rad/s)
      real(kind=REAL_BYTE) :: crls, fxbar, fybar
! ==============================================================================
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy
! === Dispersive ===============================================================
      type(data_grids), target, intent(inout) :: fg, cg
      integer(kind=4), intent(in) :: bflag
#ifdef MPI
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
! ==============================================================================
      real(kind=REAL_BYTE) :: c0, c1, c2, c3, c4
#ifdef __FUJITSU
      real(kind=REAL_BYTE), allocatable, dimension(:) :: ta0, ta1
#endif
#ifdef LESS_CC
      integer(kind=4), parameter :: check = 10
#endif
      integer(kind=4) :: numneststeps, neststepvel, neststep
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx_new, fy_new
      real(kind=REAL_BYTE) :: rate

      if(timenest /= 1) then
         fx => wfld%fx
         fy => wfld%fy
         hz => wfld%hz
      end if

      numneststeps = fg%my%numneststeps
      neststepvel  = fg%my%neststepvel

      if(neststepvel == 0) then
         if(timenest == 1) then
            fx_old => wfld%fx_a
            fx     => wfld%fx_b
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = lbound(fx,2), ubound(fx,2)
               do i = lbound(fx,1), ubound(fx,1)
                  fx_old(i,j) = fx(i,j)
               end do
            end do
            fy_old => wfld%fy_a
            fy     => wfld%fy_b
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = lbound(fy,2), ubound(fy,2)
               do i = lbound(fy,1), ubound(fy,1)
                  fy_old(i,j) = fy(i,j)
               end do
            end do
            hz     => wfld%hz_b
         end if
! === Dispersive ===============================================================
         yu => wfld%yu
         yv => wfld%yv
         cu => wfld%cu
         cv => wfld%cv
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(cu,2), ubound(cu,2)
            do i = lbound(cv,1), ubound(cv,1)
               cu(i,j) = 0.0d0
               cv(i,j) = 0.0d0
            end do
         end do
         fx_old => wfld%fx_old
         fy_old => wfld%fy_old
! ==============================================================================
#ifdef CONV_CHECK
         tdx => wfld%tdx
         tdy => wfld%tdy
#endif

         dz => dfld%dz
         dx => dfld%dx
         dy => dfld%dy

! === Dispersive ===============================================================
#ifndef MPI
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 0, nlat
            do i = -1, nlon
               fx_old(i,j) = fx(i,j)
            end do
         end do
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = -1, nlat
            do i = 0, nlon
               fy_old(i,j) = fy(i,j)
            end do
         end do
#else
#ifndef USE_GPU
!$omp single
#endif
         ist_ = 0
         jst_ = 1
         if(iand(bflag, WEST_BOUND)  /= 0) ist_ = -1
         if(iand(bflag, NORTH_BOUND) /= 0) jst_ = 0
#ifndef USE_GPU
!$omp end single
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = jst_, nlat
            do i = ist_, nlon
               fx_old(i,j) = fx(i,j)
            end do
         end do
#ifndef USE_GPU
!$omp single
#endif
         ist_ = 1
         jst_ = 0
         if(iand(bflag, WEST_BOUND)  /= 0) ist_ = 0
         if(iand(bflag, NORTH_BOUND) /= 0) jst_ = -1
#ifndef USE_GPU
!$omp end single
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = jst_, nlat
            do i = ist_, nlon
               fy_old(i,j) = fy(i,j)
            end do
         end do
#endif
! ==============================================================================
         rote = 6371.0d0 * 1000.0d0
         g = 9.8d0
         gts = g*dt/(dth*rote)

! === Coriolis force is supported on linear calc. ==============================
         jst = 1
         jnd = nlat - 1
         ist = 1
         ind = nlon - 1
         if(gflag == 1) then
            jst = 2
            jnd = nlat - 2
            ist = 2
            ind = nlon - 2
         end if
#ifdef MPI
         if(iand(bflag, NORTH_BOUND) == 0) jst = 2
         if(iand(bflag, WEST_BOUND)  == 0) ist = 2
         if(iand(bflag, SOUTH_BOUND) == 0) jnd = nlat - 1
         if(iand(bflag, EAST_BOUND)  == 0) ind = nlon - 1
#endif
! ==============================================================================
! === Coriolis force is supported on linear calc. ==============================
!!$omp parallel do private(theta, sint, cfc, i)
!     do j = 1, nlat-1
#ifndef USE_GPU
!$omp parallel do private(theta, sint, cfc, i, crls, fybar, fxbar)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint,cfc,i,crls,fybar,fxbar)
#endif
         do j = jst, jnd
#ifdef USE_GPU
            do i = ist, ind
#endif
! ==============================================================================
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            sint = sin(theta)
            cfc = gts/sint
! === Coriolis force is supported on linear calc. ==============================
!        do i = 1, nlon-1
            crls = 2.0d0*cflag*cfreq*cos(theta)*dt
#ifndef USE_GPU
            do i = ist, ind
#endif
! ==============================================================================
               if(dz(i,j) > zap) then
! === Coriolis force is supported on linear calc. ==============================
!              fx(i,j) = fx(i,j) - dx(i,j)*cfc*(hz(i+1,j)  -hz(i,j))
!              fy(i,j) = fy(i,j) - dy(i,j)*gts*(hz(i,  j+1)-hz(i,j))
               if(ifz(i,j)+ifz(i+1,j) > 0) then
                  fybar = quart*(fy_old(i,j) + fy_old(i+1,j) + fy_old(i,j-1) + fy_old(i+1,j-1))
                  fx(i,j) = fx(i,j) - dx(i,j)*cfc*(hz(i+1,j)  -hz(i,j)) - crls*fybar
               else
                  fx(i,j) = zap
               endif
               if(ifz(i,j)+ifz(i,j+1) > 0) then
                  fxbar = quart*(fx_old(i,j) + fx_old(i-1,j) + fx_old(i,j+1) + fx_old(i-1,j+1))
                  fy(i,j) = fy(i,j) - dy(i,j)*gts*(hz(i,  j+1)-hz(i,j)) + crls*fxbar
               else
                  fy(i,j) = zap
               endif
! ==============================================================================
               else
                  fx(i,j) = zap
                  fy(i,j) = zap
               end if
            end do
         end do
! === Coriolis force is supported on linear calc. ==============================
         if(gflag == 1) then
#ifndef MPI
               theta = th0
#else
            if(iand(bflag, NORTH_BOUND) /= 0) then
               theta = th0 + joff*dth
#endif
               sint = sin(theta)
               cfc = gts/sint
#ifndef MPI
#ifndef USE_GPU
!$omp parallel do
#else
!$omp target teams distribute parallel do
#endif
               do i = 1, nlon-1
#else
               ist = 1
               ind = nlon - 1
               if(iand(bflag, WEST_BOUND)  == 0) ist = 2
#ifndef USE_GPU
!$omp parallel do
#else
!$omp target teams distribute parallel do
#endif
               do i = ist, ind
#endif
                  if(dz(i,1) > zap) then
                     fx(i,1) = fx(i,1) - dx(i,1)*cfc*(hz(i+1,1)-hz(i,1))
                     fy(i,1) = fy(i,1) - dy(i,1)*gts*(hz(i,  2)-hz(i,1))
                  else
                     fx(i,1) = zap
                     fy(i,1) = zap
                  end if
               end do
#ifndef MPI
               theta = th0 + (nlat-2)*dth
#else
            end if
            if(iand(bflag, SOUTH_BOUND) /= 0) then
               theta = th0 + (nlat+joff-2)*dth
#endif
               sint = sin(theta)
               cfc = gts/sint
#ifndef MPI
#ifndef USE_GPU
!$omp parallel do
#else
!$omp target teams distribute parallel do
#endif
               do i = 1, nlon-1
#else
               ist = 1
               ind = nlon - 1
               if(iand(bflag, WEST_BOUND)  == 0) ist = 2
#ifndef USE_GPU
!$omp parallel do
#else
!$omp target teams distribute parallel do
#endif
               do i = ist, ind
#endif
                  if(dz(i,nlat-1) > zap) then
                     fx(i,nlat-1) = fx(i,nlat-1) - dx(i,nlat-1)*cfc*(hz(i+1,nlat-1)-hz(i,nlat-1))
                     fy(i,nlat-1) = fy(i,nlat-1) - dy(i,nlat-1)*gts*(hz(i,  nlat)  -hz(i,nlat-1))
                  else
                     fx(i,nlat-1) = zap
                     fy(i,nlat-1) = zap
                  end if
               end do
#ifndef MPI
#ifndef USE_GPU
!$omp parallel do private(theta, sint, cfc)
#else
!$omp target teams distribute parallel do private(theta,sint,cfc)
#endif
               do j = 2, nlat-2
                  theta = th0 + (j-1)*dth
#else
            end if
            if(iand(bflag, WEST_BOUND) /= 0) then
               jst = 2
               jnd = nlat - 2
               if(iand(bflag, SOUTH_BOUND) == 0) jnd = nlat - 1
#ifndef USE_GPU
!$omp parallel do private(theta, sint, cfc)
#else
!$omp target teams distribute parallel do private(theta,sint,cfc)
#endif
               do j = jst, jnd
                  theta = th0 + (j+joff-1)*dth
#endif
                  sint = sin(theta)
                  cfc = gts/sint
                  if(dz(1,j) > zap) then
                     fx(1,j) = fx(1,j) - dx(1,j)*cfc*(hz(2,j)  -hz(1,j))
                     fy(1,j) = fy(1,j) - dy(1,j)*gts*(hz(1,j+1)-hz(1,j))
                  else
                     fx(1,j) = zap
                     fy(1,j) = zap
                  end if
               end do
#ifndef MPI
#ifndef USE_GPU
!$omp parallel do private(theta, sint, cfc)
#else
!$omp target teams distribute parallel do private(theta,sint,cfc)
#endif
               do j = 2, nlat-2
                  theta = th0 + (j-1)*dth
#else
            end if
            if(iand(bflag, EAST_BOUND) /= 0) then
               jst = 2
               jnd = nlat - 2
               if(iand(bflag, SOUTH_BOUND) == 0) jnd = nlat - 1
#ifndef USE_GPU
!$omp parallel do private(theta, sint, cfc)
#else
!$omp target teams distribute parallel do private(theta,sint,cfc)
#endif
               do j = jst, jnd
                  theta = th0 + (j+joff-1)*dth
#endif
                  sint = sin(theta)
                  cfc = gts/sint
                  if(dz(nlon-1,j) > zap) then
                     fx(nlon-1,j) = fx(nlon-1,j) - dx(nlon-1,j)*cfc*(hz(nlon,  j)  -hz(nlon-1,j))
                     fy(nlon-1,j) = fy(nlon-1,j) - dy(nlon-1,j)*gts*(hz(nlon-1,j+1)-hz(nlon-1,j))
                  else
                     fx(nlon-1,j) = zap
                     fy(nlon-1,j) = zap
                  end if
               end do
#ifdef MPI
            end if
#endif
         end if
! ==============================================================================
! === Coriolis force is supported on linear calc. ==============================
#ifdef MPI
         call exchange_edges(VEL,fg)

         ist_ = 1
         jst_ = 1
         if(iand(bflag, WEST_BOUND)  /= 0) ist_ = 0
         if(iand(bflag, NORTH_BOUND) /= 0) jst_ = 0
#endif
! ==============================================================================
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
#ifndef USE_GPU
!$omp parallel private(conv)
#endif
#else
!$omp parallel private(conv, ta0, ta1)
         allocate(ta0(ist:ind))
         allocate(ta1(ist:ind))
#endif
#ifndef USE_GPU
!$omp do private(theta, sint, cfc, i)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint,cfc,i)
#endif
         do j = jst, jnd
#ifdef USE_GPU
            do i = ist, ind
#endif
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            sint = sin(theta)
            cfc = gts/sint
#ifndef USE_GPU
            do i = ist, ind
#endif
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
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 0, nlat
            do i = 0, nlon
               yu(i,j) = fx(i,j)
               fx(i,j) = fx_old(i,j)
               yv(i,j) = fy(i,j)
               fy(i,j) = fy_old(i,j)
            end do
         end do
         if(gflag /= 1) then
#ifndef USE_GPU
!$omp single
#endif
            neststep = cg%my%neststepvel / fg%my%numneststeps
            call interp2fine(VEL,cg,fg,neststep)
#ifndef USE_GPU
!$omp end single
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = 0, nlat
               do i = 0, nlon
                  fx_old(i,j) = fx(i,j)
                  fy_old(i,j) = fy(i,j)
               end do
            end do
         end if
#else
#ifndef USE_GPU
!$omp single
#endif
         call exchange_edges_disp_fx(fg,1)
         call exchange_edges_disp_fy(fg,1)
#ifndef USE_GPU
!$omp end single
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = jst_, nlat
            do i = ist_, nlon
               yu(i,j) = fx(i,j)
               fx(i,j) = fx_old(i,j)
               yv(i,j) = fy(i,j)
               fy(i,j) = fy_old(i,j)
            end do
         end do
         if(gflag /= 1) then
#ifndef USE_GPU
!$omp single
#endif
            neststep = cg%my%neststepvel / fg%my%numneststeps
            call interp2fine(VEL,cg,fg,neststep)
            call exchange_edges_disp_fx(fg,1)
            call exchange_edges_disp_fy(fg,1)
#ifndef USE_GPU
!$omp end single
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = jst_, nlat
               do i = ist_, nlon
                  fx_old(i,j) = fx(i,j)
                  fy_old(i,j) = fy(i,j)
               end do
            end do
         end if
! ==============================================================================
#endif
#ifdef CONV_CHECK
#ifndef USE_GPU
!$omp do private(i, tmpx, tmpy)
#else
!$omp target teams distribute parallel do collapse(2) private(i,tmpx,tmpy)
#endif
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
#ifndef USE_GPU
!$omp master
#endif
TIMER_START('___ld:reg01')
#ifndef USE_GPU
!$omp end master
#endif
#ifndef USE_GPU
!$omp do private(theta, i, c0, c1)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,i,c0,c1)
#endif
            do j = jst, jnd
#ifdef USE_GPU
               do i = ist, ind
#endif
#ifndef MPI
               theta = th0 + (j-1)*dth
#else
               theta = th0 + (j+joff-1)*dth
#endif
               c0 = sin(theta+0.5d0*dth)
               c1 = sin(theta-0.5d0*dth)
#ifndef USE_GPU
               do i = ist, ind
#endif
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
#ifndef USE_GPU
!$omp master
#endif
TIMER_STOP('___ld:reg01')
#ifndef USE_GPU
!$omp end master
#endif
#ifdef MPI
#ifndef USE_GPU
!$omp single
#endif
            call exchange_edges_disp_fx(fg,mod(conv,2))
#ifndef USE_GPU
!$omp end single
#endif
#endif
#ifndef USE_GPU
!$omp master
#endif
TIMER_START('___ld:reg02')
#ifndef USE_GPU
!$omp end master
#endif
#ifndef USE_GPU
!$omp do private(theta, i, c0, c1, c2, c3, c4)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,i,c0,c1,c2,c3,c4)
#endif
            do j = jst, jnd
#ifdef USE_GPU
               do i = ist, ind
#endif
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
#ifndef USE_GPU
               do i = ist, ind
#endif
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
#ifndef USE_GPU
!$omp master
#endif
TIMER_STOP('___ld:reg02')
#ifndef USE_GPU
!$omp end master
#endif
#ifdef MPI
#ifndef USE_GPU
!$omp single
#endif
            call exchange_edges_disp_fy(fg,mod(conv,2))
#ifndef USE_GPU
!$omp end single
#endif
#endif
#ifdef CONV_CHECK
#ifdef LESS_CC
            if(mod(conv,check) == 0) then
#endif
#ifndef USE_GPU
!$omp master
#endif
TIMER_START('___ld:reg03')
#ifndef USE_GPU
!$omp end master
#endif
#ifndef USE_GPU
!$omp single
#endif
            diffxy = 0.0d0
#ifndef USE_GPU
!$omp end single
!$omp do private(i, diff, diffx, diffy) reduction(max:diffxy)
#else
!$omp target teams distribute parallel do collapse(2) private(i,diff,diffx,diffy) reduction(max:diffxy)
#endif
            do j = jst, jnd
               do i = ist, ind
                  diffx = abs((fx(i,j)-fx_old(i,j))*tdx(i,j))
                  diffy = abs((fy(i,j)-fy_old(i,j))*tdy(i,j))
                  diff = max(diffx,diffy)
                  diffxy = max(diffxy,diff)
               end do
            end do
#ifndef USE_GPU
!$omp master
#endif
TIMER_STOP('___ld:reg03')
#ifndef USE_GPU
!$omp end master
#endif
#ifdef MPI
#ifndef USE_GPU
!$omp single
#endif
            call MPI_Allreduce(MPI_IN_PLACE, diffxy, 1, REAL_MPI, MPI_MAX, __MPICOMM__, ierr)
#ifndef USE_GPU
!$omp end single
#endif
#endif
            if(diffxy < conv_val) then
               exit
            endif
#ifdef LESS_CC
            end if
#endif
#endif
#ifndef USE_GPU
!$omp single
#endif
            p => fx_old
            fx_old => fx
            fx => p
            p => fy_old
            fy_old => fy
            fy => p
#ifndef USE_GPU
!$omp end single
#endif
         end do
#ifdef CONV_CHECK
#ifndef USE_GPU
!$omp single
#endif
         conv_step = conv
#ifndef USE_GPU
!$omp end single
#endif
#endif
#ifndef MPI
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 0, nlat
            do i = 0, nlon
               fx_old(i,j) = fx(i,j)
               fy_old(i,j) = fy(i,j)
            end do
         end do
#else
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
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
#ifndef USE_GPU
!$omp end parallel
#endif
! ==============================================================================
         if(timenest == 1) then
            fx     => wfld%fx
            fx_new => wfld%fx_b
            fy     => wfld%fy
            fy_new => wfld%fy_b
#ifdef MPI
            call exchange_edges_b(VEL,fg)
#endif
            rate = dble(neststepvel + 1)/dble(numneststeps)
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = lbound(fx,2), ubound(fx,2)
               do i = lbound(fx,1), ubound(fx,1)
                  fx(i,j) = fx_old(i,j)*(1.0d0 - rate) + fx_new(i,j)*rate
               end do
            end do
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = lbound(fy,2), ubound(fy,2)
               do i = lbound(fy,1), ubound(fy,1)
                  fy(i,j) = fy_old(i,j)*(1.0d0 - rate) + fy_new(i,j)*rate
               end do
            end do
         end if
      else if(neststepvel == numneststeps - 1) then
         fx     => wfld%fx
         fx_new => wfld%fx_b
         fy     => wfld%fy
         fy_new => wfld%fy_b
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(fx,2), ubound(fx,2)
            do i = lbound(fx,1), ubound(fx,1)
               fx(i,j) = fx_new(i,j)
            end do
         end do
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(fy,2), ubound(fy,2)
            do i = lbound(fy,1), ubound(fy,1)
               fy(i,j) = fy_new(i,j)
            end do
         end do
      else
         fx     => wfld%fx
         fx_old => wfld%fx_a
         fx_new => wfld%fx_b
         fy     => wfld%fy
         fy_old => wfld%fy_a
         fy_new => wfld%fy_b
         rate = dble(neststepvel + 1)/dble(numneststeps)
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(fx,2), ubound(fx,2)
            do i = lbound(fx,1), ubound(fx,1)
               fx(i,j) = fx_old(i,j)*(1.0d0 - rate) + fx_new(i,j)*rate
            end do
         end do
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(fy,2), ubound(fy,2)
            do i = lbound(fy,1), ubound(fy,1)
               fy(i,j) = fy_old(i,j)*(1.0d0 - rate) + fy_new(i,j)*rate
            end do
         end do
      end if

      return
   end subroutine fxy_rwg_Coriolis_disp

end module mod_fxy_Coriolis_disp
