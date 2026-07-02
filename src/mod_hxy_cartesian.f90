#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_hxy_cartesian
use mod_grid
use mod_params, only : timenest
implicit none

contains

   subroutine hxy_rwg(wfld,dt,dxdy,dth,joff,nlon,nlat,fg)
      type(data_grids), target, intent(inout) :: fg
      type(wave_arrays), target, intent(inout) :: wfld
      real(kind=REAL_BYTE), intent(in) :: dt, dxdy
      integer(kind=4), intent(in) :: nlon, nlat
! === Dummy arguments ==========================================================
      real(kind=REAL_BYTE), intent(in) :: dth
      integer(kind=4), intent(in) :: joff
! ==============================================================================

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: dtds
      integer(kind=4) :: numneststeps, neststephgt

      if(timenest /= 1) then
         fx => wfld%fx
         fy => wfld%fy
         hz => wfld%hz
      end if

      numneststeps = fg%my%numneststeps
      neststephgt  = fg%my%neststephgt

      if(neststephgt == 0) then
         if(timenest == 1) then
            fx     => wfld%fx_b
            fy     => wfld%fy_b
            hz     => wfld%hz_b
         end if

         dtds = dt/dxdy

#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 2, nlat-1
            do i = 2, nlon-1
               hz(i,j) = hz(i,j) - dtds*(fx(i,j) - fx(i-1,j) + fy(i,j) - fy(i,j-1))
            end do
         end do
      end if

      return
   end subroutine hxy_rwg

   subroutine hxynl_rwg(wfld,dfld,ifz,dt,dxdy,dth,joff,nlon,nlat,smallh,bflag,fg)
      type(data_grids), target, intent(inout) :: fg
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
      integer(kind=4), intent(in) :: bflag
! === Dummy arguments ==========================================================
      real(kind=REAL_BYTE), intent(in) :: dth
      integer(kind=4), intent(in) :: joff
! ==============================================================================

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, dz
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: dtds

      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
      real(kind=REAL_BYTE), parameter :: half = 0.5d0
#ifdef MPI
      integer(kind=4) :: ist, jst, ind, jnd
#endif
      integer(kind=4) :: numneststeps, neststephgt

      if(timenest /= 1) then
         fx => wfld%fx
         fy => wfld%fy
         hz => wfld%hz
      end if

      numneststeps = fg%my%numneststeps
      neststephgt  = fg%my%neststephgt

      if(neststephgt == 0) then
         if(timenest == 1) then
            fx     => wfld%fx_b
            fy     => wfld%fy_b
            hz     => wfld%hz_b
         end if

         dz  => dfld%dz

         dtds = dt/dxdy

#ifndef USE_GPU
!$omp parallel
#endif
#ifndef MPI
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 2, nlat-1
            do i = 2, nlon-1
#else
#ifndef USE_GPU
!$omp single
#endif
         ist = 2
         ind = nlon
         jst = 2
         jnd = nlat

         if(iand(bflag, EAST_BOUND)  /= 0) ind = ind - 1
         if(iand(bflag, SOUTH_BOUND) /= 0) jnd = jnd - 1
#ifndef USE_GPU
!$omp end single

!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = jst, jnd
            do i = ist, ind
#endif
               hz(i,j) = hz(i,j) - dtds*(fx(i,j) - fx(i-1,j) + fy(i,j) - fy(i,j-1))
            end do
         end do

         ! check for wet-or-dry
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
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
#ifndef USE_GPU
!$omp end parallel
#endif
      end if

      return
   end subroutine hxynl_rwg

end module mod_hxy_cartesian
