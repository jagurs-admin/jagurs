#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_hxy
use mod_grid
#ifndef CARTESIAN
! === Density ==================================================================
use mod_params, only : with_density
! ==============================================================================
#endif
implicit none

contains

#ifndef MPI
   subroutine hxy_rwg(wfld,dt,th0,dth,nlon,nlat)
#else
   subroutine hxy_rwg(wfld,dt,th0,dth,joff,nlon,nlat)
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: joff
#endif

#ifdef CARTESIAN
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz
#else
! === Density ==================================================================
!     real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, m_rhoC
! ==============================================================================
#endif
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: rote, dtds, theta, sint1, sint2, cfac1

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz
#ifndef CARTESIAN
! === Density ==================================================================
      if(with_density == 1) m_rhoC => wfld%m_rhoC
! ==============================================================================
#endif

      rote = 6371.0d0 * 1000.0d0
      dtds = dt/(dth*rote)

#ifndef CARTESIAN
! === Density ==================================================================
      if(with_density == 1) then
!$omp parallel do private(theta, sint1, sint2, cfac1, i)
         do j = 2, nlat-1
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            sint1 = sin(theta)
            sint2 = sin(theta + dth)
            cfac1 = dtds/sint1
            do i = 2, nlon-1
               hz(i,j) = hz(i,j) &
                  - m_rhoC(i,j)*cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
            end do
         end do
      else
! ==============================================================================
#endif
!$omp parallel do private(theta, sint1, sint2, cfac1, i)
      do j = 2, nlat-1
#ifndef MPI
         theta = th0 + (j-1)*dth
#else
         theta = th0 + (j+joff-1)*dth
#endif
         sint1 = sin(theta)
         sint2 = sin(theta + dth)
         cfac1 = dtds/sint1
         do i = 2, nlon-1
            hz(i,j) = hz(i,j) - cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
         end do
      end do
#ifndef CARTESIAN
! === Density ==================================================================
      end if
! ==============================================================================
#endif

      return
   end subroutine hxy_rwg

#ifndef MPI
   subroutine hxynl_rwg(wfld,dfld,ifz,dt,th0,dth,nlon,nlat,smallh)
#else
   subroutine hxynl_rwg(wfld,dfld,ifz,dt,th0,dth,joff,nlon,nlat,smallh,bflag)
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(inout) :: dfld
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(inout) :: ifz
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(inout) :: ifz
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(inout) :: ifz
#endif
! ==============================================================================
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: smallh
#ifdef MPI
      integer(kind=4), intent(in) :: joff, bflag
#endif

#ifdef CARTESIAN
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz
#else
! === Density ==================================================================
!     real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz, m_rhoC
! ==============================================================================
#endif
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: rote, dtds, theta, sint1, sint2, cfac1

      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
      real(kind=REAL_BYTE), parameter :: half = 0.5d0
#ifdef MPI
      integer(kind=4) :: ist, jst, ind, jnd
#endif

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz

      dz  => dfld%dz
#ifndef CARTESIAN
! === Density ==================================================================
      if(with_density == 1) m_rhoC => wfld%m_rhoC
! ==============================================================================
#endif

      rote = 6371.0d0 * 1000.0d0
      dtds = dt/(dth*rote)

#ifndef CARTESIAN
! === Density ==================================================================
      if(with_density == 1) then
!$omp parallel
#ifndef MPI
!$omp do private(theta, sint1, sint2, cfac1, i)
         do j = 2, nlat-1
            theta = th0 + (j-1)*dth
#else
!$omp single
         ist = 2
         ind = nlon
         jst = 2
         jnd = nlat

         if(iand(bflag, EAST_BOUND)  /= 0) ind = ind - 1
         if(iand(bflag, SOUTH_BOUND) /= 0) jnd = jnd - 1
!$omp end single

!$omp do private(theta, sint1, sint2, cfac1, i)
         do j = jst, jnd
            theta = th0 + (j+joff-1)*dth
#endif
            sint1 = sin(theta)
            sint2 = sin(theta + dth)
            cfac1 = dtds/sint1
#ifndef MPI
            do i = 2, nlon-1
#else
            do i = ist, ind
#endif
               hz(i,j) = hz(i,j) &
                  - m_rhoC(i,j)*cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
            end do
         end do

      ! check for wet-or-dry
!$omp do private(i)
         do j = 1, nlat
            do i = 1, nlon
               if(dz(i,j) + hz(i,j) > smallh) then
                  ifz(i,j) = 1
               else
! === Flood Change =============================================================
                  if(ifz(i,j) == 1) hz(i,j) = -dz(i,j)
! ==============================================================================
                  ifz(i,j) = -1
               end if
            end do
         end do
!$omp end parallel
      else
! ==============================================================================
#endif
!$omp parallel
#ifndef MPI
!$omp do private(theta, sint1, sint2, cfac1, i)
      do j = 2, nlat-1
         theta = th0 + (j-1)*dth
#else
!$omp single
      ist = 2
      ind = nlon
      jst = 2
      jnd = nlat

      if(iand(bflag, EAST_BOUND)  /= 0) ind = ind - 1
      if(iand(bflag, SOUTH_BOUND) /= 0) jnd = jnd - 1
!$omp end single

!$omp do private(theta, sint1, sint2, cfac1, i)
      do j = jst, jnd
         theta = th0 + (j+joff-1)*dth
#endif
         sint1 = sin(theta)
         sint2 = sin(theta + dth)
         cfac1 = dtds/sint1
#ifndef MPI
         do i = 2, nlon-1
#else
         do i = ist, ind
#endif
            hz(i,j) = hz(i,j) - cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
         end do
      end do

      ! check for wet-or-dry
!$omp do private(i)
      do j = 1, nlat
         do i = 1, nlon
            if(dz(i,j) + hz(i,j) > smallh) then
               ifz(i,j) = 1
            else
! === Flood Change =============================================================
               if(ifz(i,j) == 1) hz(i,j) = -dz(i,j)
! ==============================================================================
               ifz(i,j) = -1
            end if
         end do
      end do
!$omp end parallel
#ifndef CARTESIAN
! === Density ==================================================================
      end if
! ==============================================================================
#endif

      return
   end subroutine hxynl_rwg

end module mod_hxy
