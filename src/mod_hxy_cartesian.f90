#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_hxy_cartesian
use mod_grid
implicit none

contains

   subroutine hxy_rwg(wfld,dt,dxdy,nlon,nlat)
      type(wave_arrays), target, intent(inout) :: wfld
      real(kind=REAL_BYTE), intent(in) :: dt, dxdy
      integer(kind=4), intent(in) :: nlon, nlat

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: dtds

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz

      dtds = dt/dxdy

!$omp parallel do private(i)
      do j = 2, nlat-1
         do i = 2, nlon-1
            hz(i,j) = hz(i,j) - dtds*(fx(i,j) - fx(i-1,j) + fy(i,j) - fy(i,j-1))
         end do
      end do

      return
   end subroutine hxy_rwg

#ifndef MPI
   subroutine hxynl_rwg(wfld,dfld,ifz,dt,dxdy,nlon,nlat,smallh)
#else
   subroutine hxynl_rwg(wfld,dfld,ifz,dt,dxdy,nlon,nlat,smallh,bflag)
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
      real(kind=REAL_BYTE), intent(in) :: dt, dxdy
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: smallh
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
#endif

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: dtds

      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
      real(kind=REAL_BYTE), parameter :: half = 0.5d0
#ifdef MPI
      integer(kind=4) :: ist, jst, ind, jnd
#endif

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz

      dz  => dfld%dz

      dtds = dt/dxdy

!$omp parallel
#ifndef MPI
!$omp do private(i)
      do j = 2, nlat-1
         do i = 2, nlon-1
#else
!$omp single
      ist = 2
      ind = nlon
      jst = 2
      jnd = nlat

      if(iand(bflag, EAST_BOUND)  /= 0) ind = ind - 1
      if(iand(bflag, SOUTH_BOUND) /= 0) jnd = jnd - 1
!$omp end single

!$omp do private(i)
      do j = jst, jnd
         do i = ist, ind
#endif
            hz(i,j) = hz(i,j) - dtds*(fx(i,j) - fx(i-1,j) + fy(i,j) - fy(i,j-1))
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

      return
   end subroutine hxynl_rwg

end module mod_hxy_cartesian
