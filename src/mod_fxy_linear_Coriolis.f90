#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_fxy_Coriolis
use mod_grid
implicit none

contains

#ifndef MPI
! === Coriolis force is supported on linear calc. ==============================
!  subroutine fxy_rwg(wfld,dfld,dt,th0,dth,nlon,nlat)
   subroutine fxy_rwg_Coriolis(wfld,dfld,cflag,dt,th0,dth,nlon,nlat,gflag)
! ==============================================================================
#else
! === Coriolis force is supported on linear calc. ==============================
!  subroutine fxy_rwg(wfld,dfld,dt,th0,dth,joff,nlon,nlat)
   subroutine fxy_rwg_Coriolis(wfld,dfld,cflag,dt,th0,dth,joff,nlon,nlat,gflag,bflag)
! ==============================================================================
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: joff
#endif
! === Coriolis force is supported on linear calc. ==============================
      integer(kind=4), intent(in) :: cflag
! ==============================================================================

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: rote, g, gts, theta, sint, cfc
      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
! === Coriolis force is supported on linear calc. ==============================
      real(kind=REAL_BYTE), parameter :: quart = 0.25d0
      real(kind=REAL_BYTE), parameter :: cfreq = 7.27220522d-05 ! angular freq of earth's rotation (rad/s)
      real(kind=REAL_BYTE) :: crls, fxbar, fybar
! ==============================================================================
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy
! === Coriolis force is supported on linear calc. ==============================
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
      integer(kind=4) :: ist_, jst_
#endif
      integer(kind=4), intent(in) :: gflag
      integer(kind=4) :: ist, ind, jst, jnd
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx_old, fy_old
! ==============================================================================

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz
! === Coriolis force is supported on linear calc. ==============================
      fx_old => wfld%fx_old
      fy_old => wfld%fy_old
! ==============================================================================

      dz => dfld%dz
      dx => dfld%dx
      dy => dfld%dy

! === Coriolis force is supported on linear calc. ==============================
#ifndef MPI
!$omp do private(i)
      do j = 0, nlat
         do i = -1, nlon
            fx_old(i,j) = fx(i,j)
         end do
      end do
!$omp do private(i)
      do j = -1, nlat
         do i = 0, nlon
            fy_old(i,j) = fy(i,j)
         end do
      end do
#else
!$omp single
      ist_ = 0
      jst_ = 1
      if(iand(bflag, WEST_BOUND)  /= 0) ist_ = -1
      if(iand(bflag, NORTH_BOUND) /= 0) jst_ = 0
!$omp end single
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            fx_old(i,j) = fx(i,j)
         end do
      end do
!$omp single
      ist_ = 1
      jst_ = 0
      if(iand(bflag, WEST_BOUND)  /= 0) ist_ = 0
      if(iand(bflag, NORTH_BOUND) /= 0) jst_ = -1
!$omp end single
!$omp do private(i)
      do j = jst_, nlat
         do i = ist_, nlon
            fy_old(i,j) = fy(i,j)
         end do
      end do
#endif
! ==============================================================================
      rote = 6371.0d0 * 1000.0d0
      g = 9.8d0
      gts = g * dt /(dth * rote)

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
!$omp parallel do private(theta, sint, cfc, i, crls, fybar, fxbar)
      do j = jst, jnd
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
         do i = ist, ind
! ==============================================================================
            if(dz(i,j) > zap) then
! === Coriolis force is supported on linear calc. ==============================
!              fx(i,j) = fx(i,j) - dx(i,j)*cfc*(hz(i+1,j)  -hz(i,j))
!              fy(i,j) = fy(i,j) - dy(i,j)*gts*(hz(i,  j+1)-hz(i,j))
               fybar = quart*(fy_old(i,j) + fy_old(i+1,j) + fy_old(i,j-1) + fy_old(i+1,j-1))
               fxbar = quart*(fx_old(i,j) + fx_old(i-1,j) + fx_old(i,j+1) + fx_old(i-1,j+1))
               fx(i,j) = fx(i,j) - dx(i,j)*cfc*(hz(i+1,j)  -hz(i,j)) - crls*fybar
               fy(i,j) = fy(i,j) - dy(i,j)*gts*(hz(i,  j+1)-hz(i,j)) + crls*fxbar
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
!$omp parallel do
            do i = 1, nlon-1
#else
            ist = 1
            ind = nlon - 1
            if(iand(bflag, WEST_BOUND)  == 0) ist = 2
!$omp parallel do
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
!$omp parallel do
            do i = 1, nlon-1
#else
            ist = 1
            ind = nlon - 1
            if(iand(bflag, WEST_BOUND)  == 0) ist = 2
!$omp parallel do
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
!$omp parallel do private(theta, sint, cfc)
            do j = 2, nlat-2
               theta = th0 + (j-1)*dth
#else
         end if
         if(iand(bflag, WEST_BOUND) /= 0) then
            jst = 2
            jnd = nlat - 2
            if(iand(bflag, SOUTH_BOUND) == 0) jnd = nlat - 1
!$omp parallel do private(theta, sint, cfc)
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
!$omp parallel do private(theta, sint, cfc)
            do j = 2, nlat-2
               theta = th0 + (j-1)*dth
#else
         end if
         if(iand(bflag, EAST_BOUND) /= 0) then
            jst = 2
            jnd = nlat - 2
            if(iand(bflag, SOUTH_BOUND) == 0) jnd = nlat - 1
!$omp parallel do private(theta, sint, cfc)
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

      return
   end subroutine fxy_rwg_Coriolis

end module mod_fxy_Coriolis
