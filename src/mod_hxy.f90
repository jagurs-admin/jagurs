#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_hxy
use mod_grid
use mod_params, only : timenest
#ifndef CARTESIAN
! === Density ==================================================================
use mod_params, only : with_density
! ==============================================================================
#endif
implicit none

contains

   subroutine hxy_rwg(wfld,dt,th0,dth,joff,nlon,nlat,fg)
      type(data_grids), target, intent(inout) :: fg
      type(wave_arrays), target, intent(inout) :: wfld
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
      integer(kind=4), intent(in) :: joff

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
#ifndef USE_GPU
!$omp parallel do private(theta, sint1, sint2, cfac1, i)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint1,sint2,cfac1,i)
#endif
            do j = 2, nlat-1
#ifdef USE_GPU
               do i = 2, nlon-1
#endif
#ifndef MPI
               theta = th0 + (j-1)*dth
#else
               theta = th0 + (j+joff-1)*dth
#endif
               sint1 = sin(theta)
               sint2 = sin(theta + dth)
               cfac1 = dtds/sint1
#ifndef USE_GPU
               do i = 2, nlon-1
#endif
                  hz(i,j) = hz(i,j) &
                     - m_rhoC(i,j)*cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
               end do
            end do
         else
! ==============================================================================
#endif
#ifndef USE_GPU
!$omp parallel do private(theta, sint1, sint2, cfac1, i)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint1,sint2,cfac1,i)
#endif
         do j = 2, nlat-1
#ifdef USE_GPU
            do i = 2, nlon-1
#endif
#ifndef MPI
            theta = th0 + (j-1)*dth
#else
            theta = th0 + (j+joff-1)*dth
#endif
            sint1 = sin(theta)
            sint2 = sin(theta + dth)
            cfac1 = dtds/sint1
#ifndef USE_GPU
            do i = 2, nlon-1
#endif
               hz(i,j) = hz(i,j) - cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
            end do
         end do
#ifndef CARTESIAN
! === Density ==================================================================
         end if
! ==============================================================================
#endif
      end if

      return
   end subroutine hxy_rwg

   subroutine hxynl_rwg(wfld,dfld,ifz,dt,th0,dth,joff,nlon,nlat,smallh,bflag,fg)
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
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: smallh
      integer(kind=4), intent(in) :: joff, bflag

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
#ifndef USE_GPU
!$omp parallel
#endif
#ifndef MPI
#ifndef USE_GPU
!$omp do private(theta, sint1, sint2, cfac1, i)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint1,sint2,cfac1,i)
#endif
            do j = 2, nlat-1
#ifdef USE_GPU
               do i = 2, nlon-1
#endif
               theta = th0 + (j-1)*dth
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
#endif

#ifndef USE_GPU
!$omp do private(theta, sint1, sint2, cfac1, i)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint1,sint2,cfac1,i)
#endif
            do j = jst, jnd
#ifdef USE_GPU
               do i = ist, ind
#endif
               theta = th0 + (j+joff-1)*dth
#endif
               sint1 = sin(theta)
               sint2 = sin(theta + dth)
               cfac1 = dtds/sint1
#ifndef USE_GPU
#ifndef MPI
               do i = 2, nlon-1
#else
               do i = ist, ind
#endif
#endif
                  hz(i,j) = hz(i,j) &
                     - m_rhoC(i,j)*cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
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
         else
! ==============================================================================
#endif
#ifndef USE_GPU
!$omp parallel
#endif
#ifndef MPI
#ifndef USE_GPU
!$omp do private(theta, sint1, sint2, cfac1, i)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint1,sint2,cfac1,i)
#endif
         do j = 2, nlat-1
#ifdef USE_GPU
            do i = 2, nlon-1
#endif
            theta = th0 + (j-1)*dth
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
#endif

#ifndef USE_GPU
!$omp do private(theta, sint1, sint2, cfac1, i)
#else
!$omp target teams distribute parallel do collapse(2) private(theta,sint1,sint2,cfac1,i)
#endif
         do j = jst, jnd
#ifdef USE_GPU
            do i = ist, ind
#endif
            theta = th0 + (j+joff-1)*dth
#endif
            sint1 = sin(theta)
            sint2 = sin(theta + dth)
            cfac1 = dtds/sint1
#ifndef USE_GPU
#ifndef MPI
            do i = 2, nlon-1
#else
            do i = ist, ind
#endif
#endif
               hz(i,j) = hz(i,j) - cfac1*(fx(i,j) - fx(i-1,j) + fy(i,j)*sint2 - fy(i,j-1)*sint1)
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
#ifndef CARTESIAN
! === Density ==================================================================
         end if
! ==============================================================================
#endif
      end if

      return
   end subroutine hxynl_rwg

end module mod_hxy
