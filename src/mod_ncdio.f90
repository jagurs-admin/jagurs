#include "real.h"
module mod_ncdio
#ifdef NCDIO
use mod_grid
! === Arrival time =============================================================
!use mod_params, only : velgrd_flag, start_date, missing_value, VEL, HGT, speedgrd_flag
use mod_params, only : velgrd_flag, start_date, missing_value, VEL, HGT, &
                       speedgrd_flag, check_arrival_time
! ==============================================================================
! === Multiple rupture =========================================================
use mod_params, only : multrupt
! ==============================================================================
#if defined(MPI) && defined(ONEFILE)
use mod_onefile, only : onefile_gather_array_ncdio
#endif
implicit none
include 'netcdf.inc'
integer(kind=4), private :: stat

contains

#ifndef MPI
   subroutine open_file(base, dgrid)
#else
   subroutine open_file(base, dgrid, myrank)
#endif
      character(len=256), intent(in) :: base
      type(data_grids), target, intent(inout) :: dgrid
#ifdef MPI
      integer(kind=4), intent(in) :: myrank
#endif

      integer(kind=4) :: nx, ny
#ifndef SKIP_MAX_VEL
      integer(kind=4), pointer :: ncid, idid, mhid, mvid, hzid, vxid, vyid, speedid
#else
      integer(kind=4), pointer :: ncid, idid, mhid, hzid, vxid, vyid, speedid
#endif
#ifdef HZMINOUT
      integer(kind=4), pointer :: minhid
#endif
! === Arrival time =============================================================
      integer(kind=4), pointer :: atid
! ==============================================================================
      integer(kind=4), pointer, dimension(:) :: start, count

      integer(kind=4) :: xid, yid, tid
! === Multiple rupture =========================================================
      integer(kind=4) :: nruptid, nrupt
! ==============================================================================
      integer(kind=4), dimension(3) :: vdims
! === Multiple rupture =========================================================
      integer(kind=4), dimension(3) :: vdims_id
! ==============================================================================
      character(len=512) :: fname
#ifdef MPI
#ifndef ONEFILE
      integer(kind=4), dimension(0:5) :: mpi_id
      integer(kind=4), dimension(2) :: tmp
#endif
#endif
#ifndef MPI
      real(kind=REAL_BYTE), pointer :: mlat0, mlon0
#else
      real(kind=REAL_BYTE), pointer :: glat0, glon0
      integer(kind=4), pointer :: ix, iy, totalNx, totalNy
#endif
      real(kind=REAL_BYTE), pointer :: dxdy
      integer(kind=4), pointer :: isize, ist, ien, jsize, jst, jen
      real(kind=8) :: lat_north, lat_south, lon_east, lon_west, dx, dy
      real(kind=8), allocatable, dimension(:) :: lon, lat
      integer(kind=4) :: lonid, latid
      integer(kind=4), pointer :: timeid, stepid
      integer(kind=4) :: i, j
      character(len=128) :: att

#ifndef MPI
      write(fname,'(a,a)') trim(base), '.nc'

      nx = dgrid%my%nx
      ny = dgrid%my%ny

      dgrid%my%ncdio%isize = nx
      dgrid%my%ncdio%ist = 1
      dgrid%my%ncdio%ien = nx

      dgrid%my%ncdio%jsize = ny
      dgrid%my%ncdio%jst = 1
      dgrid%my%ncdio%jen = ny
#else
#ifndef ONEFILE
      write(fname,'(a,a,i6.6,a)') trim(base), '.', myrank, '.nc'

      nx = dgrid%my%ixend - dgrid%my%ix + 1
      ny = dgrid%my%iyend - dgrid%my%iy + 1

      dgrid%my%ncdio%isize = nx
      dgrid%my%ncdio%ist = 1 + (dgrid%my%ix - dgrid%my%kx)
      dgrid%my%ncdio%ien = dgrid%my%ncdio%ist + nx - 1

      dgrid%my%ncdio%jsize = ny
      dgrid%my%ncdio%jst = 1 + (dgrid%my%iy - dgrid%my%ky)
      dgrid%my%ncdio%jen = dgrid%my%ncdio%jst + ny - 1
#else
      nx = dgrid%my%totalNx
      ny = dgrid%my%totalNy

      dgrid%my%ncdio%isize = dgrid%my%nx
      dgrid%my%ncdio%ist = 1
      dgrid%my%ncdio%ien = dgrid%my%nx

      dgrid%my%ncdio%jsize = dgrid%my%ny
      dgrid%my%ncdio%jst = 1
      dgrid%my%ncdio%jen = dgrid%my%ny

      if(myrank == 0) then
      write(fname,'(a,a)') trim(base), '.nc'
#endif
#endif

      ncid => dgrid%my%ncdio%ncid

      idid => dgrid%my%ncdio%idid
      mhid => dgrid%my%ncdio%mhid
#ifdef HZMINOUT
      minhid => dgrid%my%ncdio%minhid
#endif
#ifndef SKIP_MAX_VEL
      mvid => dgrid%my%ncdio%mvid
#endif
! === Arrival time =============================================================
      atid => dgrid%my%ncdio%atid
! ==============================================================================

      hzid => dgrid%my%ncdio%hzid
      vxid => dgrid%my%ncdio%vxid
      vyid => dgrid%my%ncdio%vyid
      speedid => dgrid%my%ncdio%speedid

      start => dgrid%my%ncdio%start
      count => dgrid%my%ncdio%count

! === For EXTREMELY numerous initial_displacement. =============================
!     stat = nf_create(trim(fname), ior(NF_CLOBBER, NF_64BIT_OFFSET), ncid)
      stat = nf_create(trim(fname), ior(NF_CLOBBER, NF_NETCDF4), ncid)
! ==============================================================================
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      stat = nf_def_dim(ncid, 'lon', nx, xid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_def_dim(ncid, 'lat', ny, yid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
! === Multiple rupture =========================================================
      nrupt = 1
      if(multrupt == 1) nrupt = dgrid%nrupt
      stat = nf_def_dim(ncid, 'nrupt', nrupt, nruptid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
! ==============================================================================
      stat = nf_def_dim(ncid, 'time', NF_UNLIMITED, tid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

#ifdef MPI
#ifndef ONEFILE
      stat = nf_def_dim(ncid, 'mpi_params', 2, mpi_id(0))
      vdims(1) = mpi_id(0)
      stat = nf_def_var(ncid, 'total_xy',   NF_INT, 1, vdims, mpi_id(1))
      stat = nf_def_var(ncid, 'nprocs',     NF_INT, 1, vdims, mpi_id(2))
      stat = nf_def_var(ncid, 'myrank',     NF_INT, 1, vdims, mpi_id(3))
      stat = nf_def_var(ncid, 'nprocs_xy',  NF_INT, 1, vdims, mpi_id(4))
      stat = nf_def_var(ncid, 'myrank_xy',  NF_INT, 1, vdims, mpi_id(5))
#endif
#endif

      vdims(1) = xid
      vdims(2) = yid
      vdims(3) = tid
! === Multiple rupture =========================================================
      vdims_id(1) = xid
      vdims_id(2) = yid
      vdims_id(3) = nruptid
! ==============================================================================

      timeid => dgrid%my%ncdio%timeid
      stepid => dgrid%my%ncdio%stepid

      stat = nf_def_var(ncid, 'lon', NF_DOUBLE, 1, vdims(1), lonid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Longitude'
      stat = nf_put_att_text(ncid, lonid, 'long_name', len(trim(att)), att)
      att = 'Degrees'
      stat = nf_put_att_text(ncid, lonid, 'units', len(trim(att)), att)

      stat = nf_def_var(ncid, 'lat', NF_DOUBLE, 1, vdims(2), latid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Latitude'
      stat = nf_put_att_text(ncid, latid, 'long_name', len(trim(att)), att)
      att = 'Degrees'
      stat = nf_put_att_text(ncid, latid, 'units', len(trim(att)), att)

      stat = nf_def_var(ncid, 'time', NF_DOUBLE, 1, vdims(3), timeid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Time'
      stat = nf_put_att_text(ncid, timeid, 'long_name', len(trim(att)), att)
! === For ArcGIS ===============================================================
!     att = 'Seconds'
! === Support start date in time attribute. ====================================
!     att = 'seconds since 2000-01-01 00:00:00'
      att = 'seconds since ' // trim(start_date)
! ==============================================================================
! ==============================================================================
      stat = nf_put_att_text(ncid, timeid, 'units', len(trim(att)), att)

      stat = nf_def_var(ncid, 'step', NF_INT, 1, vdims(3), stepid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Model step'
      stat = nf_put_att_text(ncid, stepid, 'long_name', len(trim(att)), att)
      att = 'Steps'
      stat = nf_put_att_text(ncid, stepid, 'units', len(trim(att)), att)

! === Multiple rupture =========================================================
!     stat = nf_def_var(ncid, 'initial_displacement', NF_REAL, 2, vdims, idid)
      stat = nf_def_var(ncid, 'initial_displacement', NF_REAL, 3, vdims_id, idid)
! ==============================================================================
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Initial displacement'
      stat = nf_put_att_text(ncid, idid, 'long_name', len(trim(att)), att)
      att = 'Meters'
      stat = nf_put_att_text(ncid, idid, 'units', len(trim(att)), att)

      stat = nf_def_var(ncid, 'max_height', NF_REAL, 2, vdims, mhid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Maximum wave height'
      stat = nf_put_att_text(ncid, mhid, 'long_name', len(trim(att)), att)
      att = 'Meters'
      stat = nf_put_att_text(ncid, mhid, 'units', len(trim(att)), att)
! === For negative max. height =================================================
!     stat = nf_put_att_real(ncid, mhid, 'missing_value', NF_REAL, 1, &
      stat = nf_put_att_real(ncid, mhid, '_FillValue', NF_REAL, 1, &
                             real(missing_value))
! ==============================================================================
#ifdef HZMINOUT
      stat = nf_def_var(ncid, 'min_height', NF_REAL, 2, vdims, minhid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Minimum wave height'
      stat = nf_put_att_text(ncid, minhid, 'long_name', len(trim(att)), att)
      att = 'Meters'
      stat = nf_put_att_text(ncid, minhid, 'units', len(trim(att)), att)
! === For negative min. height =================================================
!     stat = nf_put_att_real(ncid, minhid, 'missing_value', NF_REAL, 1, &
      stat = nf_put_att_real(ncid, minhid, '_FillValue', NF_REAL, 1, &
                             real(missing_value))
! ==============================================================================
#endif

#ifndef SKIP_MAX_VEL
      stat = nf_def_var(ncid, 'max_velocity', NF_REAL, 2, vdims, mvid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Maximum flow speed'
      stat = nf_put_att_text(ncid, mvid, 'long_name', len(trim(att)), att)
      att = 'Meters/Second'
      stat = nf_put_att_text(ncid, mvid, 'units', len(trim(att)), att)
#endif
! === Arrival time =============================================================
      if(check_arrival_time == 1) then
         stat = nf_def_var(ncid, 'arrival_time', NF_REAL, 2, vdims, atid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         att = 'Arrival time'
         stat = nf_put_att_text(ncid, atid, 'long_name', len(trim(att)), att)
         att = 'seconds since ' // trim(start_date)
         stat = nf_put_att_text(ncid, atid, 'units', len(trim(att)), att)
         stat = nf_put_att_real(ncid, atid, '_FillValue', NF_REAL, 1, &
                                real(missing_value))
      end if
! ==============================================================================

      stat = nf_def_var(ncid, 'wave_height', NF_REAL, 3, vdims, hzid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Wave height'
      stat = nf_put_att_text(ncid, hzid, 'long_name', len(trim(att)), att)
      att = 'Meters'
      stat = nf_put_att_text(ncid, hzid, 'units', len(trim(att)), att)
! === Wave height should be missing value on dry cell. =========================
      stat = nf_put_att_real(ncid, hzid, '_FillValue', NF_REAL, 1, &
                             real(missing_value))
! ==============================================================================

      if(velgrd_flag == 1) then
         stat = nf_def_var(ncid, 'velocity_x', NF_REAL, 3, vdims, vxid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         att = 'Flow velocity (X-component)'
         stat = nf_put_att_text(ncid, vxid, 'long_name', len(trim(att)), att)
         att = 'Meters/Second'
         stat = nf_put_att_text(ncid, vxid, 'units', len(trim(att)), att)

         stat = nf_def_var(ncid, 'velocity_y', NF_REAL, 3, vdims, vyid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         att = 'Flow velocity (Y-component)'
         stat = nf_put_att_text(ncid, vyid, 'long_name', len(trim(att)), att)
         att = 'Meters/Second'
         stat = nf_put_att_text(ncid, vyid, 'units', len(trim(att)), att)
      end if
      if(speedgrd_flag == 1) then
         stat = nf_def_var(ncid, 'speed', NF_REAL, 3, vdims, speedid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         att = 'Flow speed'
         stat = nf_put_att_text(ncid, speedid, 'long_name', len(trim(att)), att)
         att = 'Meters/Second'
         stat = nf_put_att_text(ncid, speedid, 'units', len(trim(att)), att)
      end if

      stat = nf_enddef(ncid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      start(1) = 1
      start(2) = 1
      start(3) = 0
      count(1) = nx
      count(2) = ny
      count(3) = 1

#ifdef MPI
#ifndef ONEFILE
      tmp(1) = dgrid%my%totalNx
      tmp(2) = dgrid%my%totalNy
      stat = nf_put_var_int(ncid, mpi_id(1), tmp)

      tmp(1) = dgrid%my%px * dgrid%my%py
      tmp(2) = 0
      stat = nf_put_var_int(ncid, mpi_id(2), tmp)

      tmp(1) = myrank
      tmp(2) = 0
      stat = nf_put_var_int(ncid, mpi_id(3), tmp)

      tmp(1) = dgrid%my%px
      tmp(2) = dgrid%my%py
      stat = nf_put_var_int(ncid, mpi_id(4), tmp)

      tmp(1) = dgrid%my%rx
      tmp(2) = dgrid%my%ry
      stat = nf_put_var_int(ncid, mpi_id(5), tmp)
#endif
#endif

#ifdef MPI
      totalNx => dgrid%my%totalNx
      totalNy => dgrid%my%totalNy
      ix => dgrid%my%ix
      iy => dgrid%my%iy
#endif
      isize => dgrid%my%ncdio%isize
      ist   => dgrid%my%ncdio%ist
      ien   => dgrid%my%ncdio%ien
      jsize => dgrid%my%ncdio%jsize
      jst   => dgrid%my%ncdio%jst
      jen   => dgrid%my%ncdio%jen

#ifndef MPI
      mlon0 => dgrid%my%mlon0
      mlat0 => dgrid%my%mlat0
      dxdy  => dgrid%my%dh

#ifndef CARTESIAN
      lat_north = 90.0d0 - mlat0/60.0d0
      lon_west  = mlon0/60.0d0
      lat_south = 90.0d0 - (mlat0/60.0d0 + (ny-1.0d0)*dxdy)
      lon_east  = mlon0/60.0d0 + (nx-1.0d0)*dxdy
      if(lon_west > 180.0d0) lon_west = lon_west - 360.0d0
      if(lon_east > 180.0d0) lon_east = lon_east - 360.0d0
#else
      lon_west  = mlon0
      lat_south = mlat0
      lon_east  = lon_west  + (nx-1.0d0)*dxdy
      lat_north = lat_south + (ny-1.0d0)*dxdy
#endif
      dx = dxdy
      dy = dx
#else
      glon0 => dgrid%my%glon0
      glat0 => dgrid%my%glat0
      dxdy  => dgrid%my%dh

#ifndef CARTESIAN
      lat_north = 90.0d0 - glat0/60.0d0
      lon_west  = glon0/60.0d0
      lat_south = 90.0d0 - (glat0/60.0d0 + (totalNy-1.0d0)*dxdy)
      lon_east  = glon0/60.0d0 + (totalNx-1.0d0)*dxdy
      if(lon_west > 180.0d0) lon_west = lon_west - 360.0d0
      if(lon_east > 180.0d0) lon_east = lon_east - 360.0d0
#else
      lon_west  = glon0
      lat_south = glat0
      lon_east  = lon_west  + (totalNx-1.0d0)*dxdy
      lat_north = lat_south + (totalNy-1.0d0)*dxdy
#endif
      dx = dxdy
      dy = dx
#endif

#if !defined(MPI) || !defined(ONEFILE)
      allocate(lon(isize))
      allocate(lat(jsize))
#else
      allocate(lon(totalNx))
      allocate(lat(totalNy))
#endif

#ifndef MPI
      do i = ist, ien
         lon(i-ist+1) = lon_west + dx*(i-1)
      end do

      do j = jst, jen
         lat(j-jst+1) = lat_north - dy*(jen-j)
      end do
#else
#ifndef ONEFILE
      do i = ix, ix+isize-1
         lon(i-ix+1) = lon_west + dx*(i-1)
      end do

      do j = iy, iy+jsize-1
         lat(j-iy+1) = lat_north - dy*(2*iy+jsize-j-2)
      end do
#else
      do i = 1, totalNx
         lon(i) = lon_west + dx*(i-1)
      end do

      do j = 1, totalNy
         lat(j) = lat_north - dy*(totalNy-j)
      end do
#endif
#endif

      stat = nf_put_var_double(ncid, lonid, lon)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_put_var_double(ncid, latid, lat)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      deallocate(lon)
      deallocate(lat)
#if defined(MPI) && defined(ONEFILE)
      end if
#endif

      return
   end subroutine open_file

! === Conversion from flux to velocity should be done right after calc. ========
!  subroutine write_snapshot(dgrid, t, istep)
#ifndef MPI
   subroutine write_snapshot(dgrid, t, istep, mode, linear_flag)
#else
#ifndef ONEFILE
   subroutine write_snapshot(dgrid, t, istep, mode, bflag, linear_flag)
#else
   subroutine write_snapshot(dgrid, t, istep, mode, bflag, myrank, linear_flag)
#endif
#endif
! ==============================================================================
      type(data_grids), target, intent(inout) :: dgrid
      real(kind=REAL_BYTE), intent(in) :: t
      integer(kind=4), intent(in) :: istep
! === Conversion from flux to velocity should be done right after calc. ========
      integer(kind=4), intent(in) :: mode
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
#ifdef ONEFILE
      integer(kind=4), intent(in) :: myrank
#endif
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: linear_flag

      integer(kind=4), pointer, dimension(:) :: start, count
      integer(kind=4), pointer :: ncid, timeid, stepid, hzid, vxid, vyid, speedid
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, fx, fy
      integer(kind=4), pointer, dimension(:,:) :: wod
      integer(kind=4), pointer :: isize, ist, ien, jsize, jst, jen
      real(kind=4), allocatable, dimension(:,:) :: tmp
      real(kind=8) :: time
#if !defined(MPI) || !defined(ONEFILE)
      integer(kind=4) :: i, j, i_, j_
#else
      integer(kind=4) :: i, j
#endif
! === Flux -> Flow speed =======================================================
! === Conversion from flux to velocity should be done right after calc. ========
!     real(kind=REAL_BYTE) :: td
      real(kind=REAL_BYTE) :: tdxm, tdxp, tdym, tdyp
      integer(kind=4) :: im, jm, ip, jp
! ==============================================================================
      real(kind=REAL_BYTE), parameter :: td_min = 0.01d0 ! 1 cm
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dz
! ==============================================================================
      real(kind=REAL_BYTE) :: tx, ty, speed
#if defined(MPI) && defined(ONEFILE)
      real(kind=4), allocatable, dimension(:,:) :: tmp_all
#endif

      start => dgrid%my%ncdio%start
      count => dgrid%my%ncdio%count

      ncid => dgrid%my%ncdio%ncid
      timeid => dgrid%my%ncdio%timeid
      stepid => dgrid%my%ncdio%stepid
      hzid => dgrid%my%ncdio%hzid
      vxid => dgrid%my%ncdio%vxid
      vyid => dgrid%my%ncdio%vyid
      speedid => dgrid%my%ncdio%speedid

      hz => dgrid%wave_field%hz
      fx => dgrid%wave_field%fx
      fy => dgrid%wave_field%fy
      wod => dgrid%wod_flags
! === Flux -> Flow speed =======================================================
      dz => dgrid%depth_field%dz
! ==============================================================================

      isize => dgrid%my%ncdio%isize
      ist   => dgrid%my%ncdio%ist
      ien   => dgrid%my%ncdio%ien
      jsize => dgrid%my%ncdio%jsize
      jst   => dgrid%my%ncdio%jst
      jen   => dgrid%my%ncdio%jen

      if(((velgrd_flag == 0) .and. (mode == HGT)) .or. &
         ((velgrd_flag == 1) .and. (mode == VEL)) .or. &
         ((speedgrd_flag == 1) .and. (mode == VEL))) then
      start(3) = start(3) + 1

      time = t
#if defined(MPI) && defined(ONEFILE)
      if(myrank == 0) then
#endif
      stat = nf_put_vara_double(ncid, timeid, start(3), count(3), time)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_put_vara_int(ncid, stepid, start(3), count(3), istep)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
#if defined(MPI) && defined(ONEFILE)
      end if
#endif
! === Conversion from flux to velocity should be done right after calc. ========
      end if
! ==============================================================================

      allocate(tmp(isize,jsize))

! === Conversion from flux to velocity should be done right after calc. ========
      if(mode == HGT) then
! ==============================================================================
      do j = 1, jsize
         do i = 1, isize
#if !defined(MPI) || !defined(ONEFILE)
            i_ = i + ist - 1
            j_ = j + jst - 1
            if(wod(i_, j_) == 1) then
               tmp(i, jsize-j+1) = hz(i_, j_)
            else
! === Wave height should be missing value on dry cell. =========================
!              tmp(i, jsize-j+1) = 0.0d0
               tmp(i, jsize-j+1) = missing_value
! ==============================================================================
            end if
#else
            if(wod(i,j) == 1) then
               tmp(i,j) = hz(i,j)
            else
               tmp(i,j) = missing_value
            end if
#endif
         end do
      end do
#if !defined(MPI) || !defined(ONEFILE)
      stat = nf_put_vara_real(ncid, hzid, start, count, tmp)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
#else
      if(myrank == 0) then
         allocate(tmp_all(dgrid%my%totalNx,dgrid%my%totalNy))
      else
         allocate(tmp_all(1,1))
      end if
      call onefile_gather_array_ncdio(tmp,tmp_all,dgrid)
      if(myrank == 0) then
         stat = nf_put_vara_real(ncid, hzid, start, count, tmp_all)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
      deallocate(tmp_all)
#endif
! === Conversion from flux to velocity should be done right after calc. ========
      end if
! ==============================================================================

! === Conversion from flux to velocity should be done right after calc. ========
      if(mode == VEL) then
! ==============================================================================
      if(velgrd_flag == 1) then
! === Conversion from flux to velocity should be done right after calc. ========
         if(linear_flag == 0) then
            do j = 1, jsize
               do i = 1, isize
#if !defined(MPI) || !defined(ONEFILE)
                  i_ = i + ist - 1
                  j_ = j + jst - 1
#ifndef MPI
                  im = max(1,  i_-1)
                  ip = min(ien,i_+1)
#else
                  im = i_ - 1
                  ip = i_ + 1
                  if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                  if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
#endif
                  tdxm = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(im,j_) + hz(im,j_))
                  tdxp = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(ip,j_) + hz(ip,j_))
                  if(wod(im,j_) == 1 .and. tdxm > td_min .and. &
                     wod(ip,j_) == 1 .and. tdxp > td_min .and. &
                     wod(i_,j_) == 1) then
                     tmp(i, jsize-j+1) = 0.5d0*(fx(i_,j_)/tdxp + fx(im,j_)/tdxm)
                  else
                     tmp(i, jsize-j+1) = 0.0d0
                  end if
#else
                  im = i - 1
                  ip = i + 1
                  if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                  if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
                  tdxm = 0.5d0*(dz(i,j) + hz(i,j) + dz(im,j) + hz(im,j))
                  tdxp = 0.5d0*(dz(i,j) + hz(i,j) + dz(ip,j) + hz(ip,j))
                  if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                     wod(ip,j) == 1 .and. tdxp > td_min .and. &
                     wod(i,j) == 1) then
                     tmp(i,j) = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
                  else
                     tmp(i,j) = 0.0d0
                  end if
#endif
               end do
            end do
         else
            do j = 1, jsize
               do i = 1, isize
#if !defined(MPI) || !defined(ONEFILE)
                  i_ = i + ist - 1
                  j_ = j + jst - 1
#ifndef MPI
                  im = max(1,  i_-1)
                  ip = min(ien,i_+1)
#else
                  im = i_ - 1
                  ip = i_ + 1
                  if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                  if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
#endif
                  tdxm = 0.5d0*(dz(i_,j_) + dz(im,j_))
                  tdxp = 0.5d0*(dz(i_,j_) + dz(ip,j_))
                  if(wod(im,j_) == 1 .and. tdxm > td_min .and. &
                     wod(ip,j_) == 1 .and. tdxp > td_min .and. &
                     wod(i_,j_) == 1) then
                     tmp(i, jsize-j+1) = 0.5d0*(fx(i_,j_)/tdxp + fx(im,j_)/tdxm)
                  else
                     tmp(i, jsize-j+1) = 0.0d0
                  end if
#else
                  im = i - 1
                  ip = i + 1
                  if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                  if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
                  tdxm = 0.5d0*(dz(i,j) + + dz(im,j))
                  tdxp = 0.5d0*(dz(i,j) + + dz(ip,j))
                  if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                     wod(ip,j) == 1 .and. tdxp > td_min .and. &
                     wod(i,j) == 1) then
                     tmp(i,j) = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
                  else
                     tmp(i,j) = 0.0d0
                  end if
#endif
               end do
            end do
         end if
! ==============================================================================
#if !defined(MPI) || !defined(ONEFILE)
         stat = nf_put_vara_real(ncid, vxid, start, count, tmp)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
#else
         if(myrank == 0) then
            allocate(tmp_all(dgrid%my%totalNx,dgrid%my%totalNy))
         else
            allocate(tmp_all(1,1))
         end if
         call onefile_gather_array_ncdio(tmp,tmp_all,dgrid)
         if(myrank == 0) then
            stat = nf_put_vara_real(ncid, vxid, start, count, tmp_all)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         end if
         deallocate(tmp_all)
#endif

! === Conversion from flux to velocity should be done right after calc. ========
         if(linear_flag == 0) then
            do j = 1, jsize
               do i = 1, isize
#if !defined(MPI) || !defined(ONEFILE)
                  i_ = i + ist - 1
                  j_ = j + jst - 1
#ifndef MPI
                  jm = max(1,  j_-1)
                  jp = min(jen,j_+1)
#else
                  jm = j_ - 1
                  jp = j_ + 1
                  if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                  if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
#endif
                  tdym = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(i_,jm) + hz(i_,jm))
                  tdyp = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(i_,jp) + hz(i_,jp))
                  if(wod(i_,jm) == 1 .and. tdym > td_min .and. &
                     wod(i_,jp) == 1 .and. tdyp > td_min .and. &
                     wod(i_,j_) == 1) then
                     ! Invert latitude!!!
                     tmp(i, jsize-j+1) = -0.5d0*(fy(i_,j_)/tdyp + fy(i_,jm)/tdym)
                  else
                     tmp(i, jsize-j+1) = 0.0d0
                  end if
#else
                  jm = j - 1
                  jp = j + 1
                  if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                  if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
                  tdym = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jm) + hz(i,jm))
                  tdyp = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jp) + hz(i,jp))
                  if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                     wod(i,jp) == 1 .and. tdyp > td_min .and. &
                     wod(i,j) == 1) then
                     ! Invert latitude!!!
                     tmp(i,j) = -0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
                  else
                     tmp(i,j) = 0.0d0
                  end if
#endif
               end do
            end do
         else
            do j = 1, jsize
               do i = 1, isize
#if !defined(MPI) || !defined(ONEFILE)
                  i_ = i + ist - 1
                  j_ = j + jst - 1
#ifndef MPI
                  jm = max(1,  j_-1)
                  jp = min(jen,j_+1)
#else
                  jm = j_ - 1
                  jp = j_ + 1
                  if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                  if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
#endif
                  tdym = 0.5d0*(dz(i_,j_) + dz(i_,jm))
                  tdyp = 0.5d0*(dz(i_,j_) + dz(i_,jp))
                  if(wod(i_,jm) == 1 .and. tdym > td_min .and. &
                     wod(i_,jp) == 1 .and. tdyp > td_min .and. &
                     wod(i_,j_) == 1) then
                     ! Invert latitude!!!
                     tmp(i, jsize-j+1) = -0.5d0*(fy(i_,j_)/tdyp + fy(i_,jm)/tdym)
                  else
                     tmp(i, jsize-j+1) = 0.0d0
                  end if
#else
                  jm = j - 1
                  jp = j + 1
                  if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                  if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
                  tdym = 0.5d0*(dz(i,j) + dz(i,jm))
                  tdyp = 0.5d0*(dz(i,j) + dz(i,jp))
                  if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                     wod(i,jp) == 1 .and. tdyp > td_min .and. &
                     wod(i,j) == 1) then
                     ! Invert latitude!!!
                     tmp(i,j) = -0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
                  else
                     tmp(i,j) = 0.0d0
                  end if
#endif
               end do
            end do
         end if
! ==============================================================================
#if !defined(MPI) || !defined(ONEFILE)
         stat = nf_put_vara_real(ncid, vyid, start, count, tmp)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
#else
         if(myrank == 0) then
            allocate(tmp_all(dgrid%my%totalNx,dgrid%my%totalNy))
         else
            allocate(tmp_all(1,1))
         end if
         call onefile_gather_array_ncdio(tmp,tmp_all,dgrid)
         if(myrank == 0) then
            stat = nf_put_vara_real(ncid, vyid, start, count, tmp_all)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         end if
         deallocate(tmp_all)
#endif
      end if
         if(speedgrd_flag == 1) then
            if(linear_flag == 0) then
               do j = 1, jsize
                  do i = 1, isize
                     tx = 0.0d0

#if !defined(MPI) || !defined(ONEFILE)
                     i_ = i + ist - 1
                     j_ = j + jst - 1
#ifndef MPI
                     im = max(1,  i_-1)
                     ip = min(ien,i_+1)
#else
                     im = i_ - 1
                     ip = i_ + 1
                     if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                     if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
#endif
                     tdxm = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(im,j_) + hz(im,j_))
                     tdxp = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(ip,j_) + hz(ip,j_))
                     if(wod(im,j_) == 1 .and. tdxm > td_min .and. &
                        wod(ip,j_) == 1 .and. tdxp > td_min .and. &
                        wod(i_,j_) == 1) then
                        tx = 0.5d0*(fx(i_,j_)/tdxp + fx(im,j_)/tdxm)
                     else
                        tx = 0.0d0
                     end if

#ifndef MPI
                     jm = max(1,  j_-1)
                     jp = min(jen,j_+1)
#else
                     jm = j_ - 1
                     jp = j_ + 1
                     if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                     if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
#endif
                     tdym = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(i_,jm) + hz(i_,jm))
                     tdyp = 0.5d0*(dz(i_,j_) + hz(i_,j_) + dz(i_,jp) + hz(i_,jp))
                     if(wod(i_,jm) == 1 .and. tdym > td_min .and. &
                        wod(i_,jp) == 1 .and. tdyp > td_min .and. &
                        wod(i_,j_) == 1) then
                        ! Invert latitude!!!
                        ty = -0.5d0*(fy(i_,j_)/tdyp + fy(i_,jm)/tdym)
                     else
                        ty = 0.0d0
                     end if

                     speed = sqrt(tx**2 + ty**2)
                     tmp(i, jsize-j+1) = speed
#else
                     im = i - 1
                     ip = i + 1
                     if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                     if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
                     tdxm = 0.5d0*(dz(i,j) + hz(i,j) + dz(im,j) + hz(im,j))
                     tdxp = 0.5d0*(dz(i,j) + hz(i,j) + dz(ip,j) + hz(ip,j))
                     if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                        wod(ip,j) == 1 .and. tdxp > td_min .and. &
                        wod(i,j) == 1) then
                        tx = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
                     else
                        tx = 0.0d0
                     end if

                     jm = j - 1
                     jp = j + 1
                     if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                     if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
                     tdym = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jm) + hz(i,jm))
                     tdyp = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jp) + hz(i,jp))
                     if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                        wod(i,jp) == 1 .and. tdyp > td_min .and. &
                        wod(i,j) == 1) then
                        ! Invert latitude!!!
                        ty = -0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
                     else
                        ty = 0.0d0
                     end if

                     speed = sqrt(tx**2 + ty**2)
                     tmp(i,j) = speed
#endif
                  end do
               end do
            else
               do j = 1, jsize
                  do i = 1, isize
                     tx = 0.0d0

#if !defined(MPI) || !defined(ONEFILE)
                     i_ = i + ist - 1
                     j_ = j + jst - 1
#ifndef MPI
                     im = max(1,  i_-1)
                     ip = min(ien,i_+1)
#else
                     im = i_ - 1
                     ip = i_ + 1
                     if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                     if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
#endif
                     tdxm = 0.5d0*(dz(i_,j_) + dz(im,j_))
                     tdxp = 0.5d0*(dz(i_,j_) + dz(ip,j_))
                     if(wod(im,j_) == 1 .and. tdxm > td_min .and. &
                        wod(ip,j_) == 1 .and. tdxp > td_min .and. &
                        wod(i_,j_) == 1) then
                        tx = 0.5d0*(fx(i_,j_)/tdxp + fx(im,j_)/tdxm)
                     else
                        tx = 0.0d0
                     end if

#ifndef MPI
                     jm = max(1,  j_-1)
                     jp = min(jen,j_+1)
#else
                     jm = j_ - 1
                     jp = j_ + 1
                     if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                     if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
#endif
                     tdym = 0.5d0*(dz(i_,j_) + dz(i_,jm))
                     tdyp = 0.5d0*(dz(i_,j_) + dz(i_,jp))
                     if(wod(i_,jm) == 1 .and. tdym > td_min .and. &
                        wod(i_,jp) == 1 .and. tdyp > td_min .and. &
                        wod(i_,j_) == 1) then
                        ! Invert latitude!!!
                        ty = -0.5d0*(fy(i_,j_)/tdyp + fy(i_,jm)/tdym)
                     else
                        ty = 0.0d0
                     end if

                     speed = sqrt(tx**2 + ty**2)
                     tmp(i, jsize-j+1) = speed
#else
                     im = i - 1
                     ip = i + 1
                     if(iand(bflag, WEST_BOUND) /= 0) im = max(1,  im)
                     if(iand(bflag, EAST_BOUND) /= 0) ip = min(ien,ip)
                     tdxm = 0.5d0*(dz(i,j) + dz(im,j))
                     tdxp = 0.5d0*(dz(i,j) + dz(ip,j))
                     if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                        wod(ip,j) == 1 .and. tdxp > td_min .and. &
                        wod(i,j) == 1) then
                        tx = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
                     else
                        tx = 0.0d0
                     end if

                     jm = j - 1
                     jp = j + 1
                     if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,  jm)
                     if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(jen,jp)
                     tdym = 0.5d0*(dz(i,j) + dz(i,jm))
                     tdyp = 0.5d0*(dz(i,j) + dz(i,jp))
                     if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                        wod(i,jp) == 1 .and. tdyp > td_min .and. &
                        wod(i,j) == 1) then
                        ! Invert latitude!!!
                        ty = -0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
                     else
                        ty = 0.0d0
                     end if

                     speed = sqrt(tx**2 + ty**2)
                     tmp(i,j) = speed
#endif
                  end do
               end do
            end if

#if !defined(MPI) || !defined(ONEFILE)
            stat = nf_put_vara_real(ncid, speedid, start, count, tmp)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
#else
            if(myrank == 0) then
               allocate(tmp_all(dgrid%my%totalNx,dgrid%my%totalNy))
            else
               allocate(tmp_all(1,1))
            end if
            call onefile_gather_array_ncdio(tmp,tmp_all,dgrid)
            if(myrank == 0) then
               stat = nf_put_vara_real(ncid, speedid, start, count, tmp_all)
               if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
            end if
            deallocate(tmp_all)
#endif
         end if
! === Conversion from flux to velocity should be done right after calc. ========
      end if
! ==============================================================================

      deallocate(tmp)

      return
   end subroutine write_snapshot

! === Multiple rupture =========================================================
!  subroutine write_initial_displacement(dgrid)
#if !defined(MPI) || !defined(ONEFILE)
   subroutine write_initial_displacement(dgrid, irupt)
#else
   subroutine write_initial_displacement(dgrid, myrank, irupt)
#endif
! ==============================================================================
      type(data_grids), target, intent(inout) :: dgrid
! === Multiple rupture =========================================================
      integer(kind=4), intent(in) :: irupt
! ==============================================================================
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4), intent(in) :: myrank
#endif
      integer(kind=4), pointer :: idid
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: zz

      idid => dgrid%my%ncdio%idid
      zz => dgrid%zz

! === Multiple rupture =========================================================
!     call write_array(dgrid, zz, idid)
#if !defined(MPI) || !defined(ONEFILE)
      call write_array(dgrid, zz, idid, irupt)
#else
      call write_array(dgrid, zz, idid, myrank, irupt)
#endif
! ==============================================================================

      return
   end subroutine write_initial_displacement

#if !defined(MPI) || !defined(ONEFILE)
   subroutine write_max_height(dgrid)
#else
   subroutine write_max_height(dgrid, myrank)
#endif
      type(data_grids), target, intent(inout) :: dgrid
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4), intent(in) :: myrank
#endif
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hzmax
      integer(kind=4), pointer :: mhid

      mhid => dgrid%my%ncdio%mhid
      hzmax => dgrid%hzmax

#if !defined(MPI) || !defined(ONEFILE)
      call write_array(dgrid, hzmax, mhid)
#else
      call write_array(dgrid, hzmax, mhid, myrank)
#endif

      return
   end subroutine write_max_height
! === Arrival time =============================================================
#if !defined(MPI) || !defined(ONEFILE)
   subroutine write_arrival_time(dgrid)
#else
   subroutine write_arrival_time(dgrid, myrank)
#endif
      type(data_grids), target, intent(inout) :: dgrid
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4), intent(in) :: myrank
#endif
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: arrival_time
      integer(kind=4), pointer :: atid

      atid => dgrid%my%ncdio%atid
      arrival_time => dgrid%wave_field%arrival_time

#if !defined(MPI) || !defined(ONEFILE)
      call write_array(dgrid, arrival_time, atid)
#else
      call write_array(dgrid, arrival_time, atid, myrank)
#endif

      return
   end subroutine write_arrival_time
! ==============================================================================
#ifdef HZMINOUT
#if !defined(MPI) || !defined(ONEFILE)
   subroutine write_min_height(dgrid)
#else
   subroutine write_min_height(dgrid, myrank)
#endif
      type(data_grids), target, intent(inout) :: dgrid
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4), intent(in) :: myrank
#endif
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hzmin
      integer(kind=4), pointer :: minhid

      minhid => dgrid%my%ncdio%minhid
      hzmin => dgrid%hzmin

#if !defined(MPI) || !defined(ONEFILE)
      call write_array(dgrid, hzmin, minhid)
#else
      call write_array(dgrid, hzmin, minhid, myrank)
#endif

      return
   end subroutine write_min_height
#endif

#ifndef SKIP_MAX_VEL
#if !defined(MPI) || !defined(ONEFILE)
   subroutine write_max_velocity(dgrid)
#else
   subroutine write_max_velocity(dgrid, myrank)
#endif
      type(data_grids), target, intent(inout) :: dgrid
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4), intent(in) :: myrank
#endif
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: vmax
      integer(kind=4), pointer :: mvid

      vmax => dgrid%vmax
      mvid => dgrid%my%ncdio%mvid

#if !defined(MPI) || !defined(ONEFILE)
      call write_array(dgrid, vmax, mvid)
#else
      call write_array(dgrid, vmax, mvid, myrank)
#endif

      return
   end subroutine write_max_velocity
#endif

   ! NOTE: write_array cannot be used for fx, fy.
   !       Because they are NOT 1-origin.
! === Multiple rupture =========================================================
!  subroutine write_array(dgrid, array, id)
#if !defined(MPI) || !defined(ONEFILE)
   subroutine write_array(dgrid, array, id, irupt)
#else
   subroutine write_array(dgrid, array, id, myrank, irupt)
#endif
! ==============================================================================
      type(data_grids), target, intent(inout) :: dgrid
      real(kind=REAL_BYTE), dimension(:,:), intent(in) :: array
      integer(kind=4), intent(in) :: id
! === Multiple rupture =========================================================
      integer(kind=4), intent(in), optional :: irupt
! ==============================================================================
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4), intent(in) :: myrank
#endif

      integer(kind=4), pointer :: ncid, isize, ist, ien, jsize, jst, jen
      integer(kind=4), pointer, dimension(:) :: start, count
      real(kind=4), allocatable, dimension(:,:) :: tmp
      integer(kind=4) :: i, j
! === Multiple rupture =========================================================
      integer(kind=4), dimension(3) :: start_id
! ==============================================================================
#if defined(MPI) && defined(ONEFILE)
      real(kind=4), allocatable, dimension(:,:) :: tmp_all
#endif

      ncid  => dgrid%my%ncdio%ncid
      isize => dgrid%my%ncdio%isize
      ist   => dgrid%my%ncdio%ist
      ien   => dgrid%my%ncdio%ien
      jsize => dgrid%my%ncdio%jsize
      jst   => dgrid%my%ncdio%jst
      jen   => dgrid%my%ncdio%jen

      start => dgrid%my%ncdio%start
! === Multiple rupture =========================================================
      if(present(irupt)) then
         start_id = start
         start_id(3) = irupt
      end if
! ==============================================================================
      count => dgrid%my%ncdio%count

      allocate(tmp(isize,jsize))

      do j = 1, jsize
         do i = 1, isize
#if !defined(MPI) || !defined(ONEFILE)
            tmp(i,jsize-j+1) = array(i+ist-1,j+jst-1)
#else
            tmp(i,j) = array(i,j)
#endif
         end do
      end do
! === Multiple rupture =========================================================
      if(present(irupt)) then
#if !defined(MPI) || !defined(ONEFILE)
         stat = nf_put_vara_real(ncid, id, start_id, count, tmp)
#else
         if(myrank == 0) then
            allocate(tmp_all(dgrid%my%totalNx,dgrid%my%totalNy))
         else
            allocate(tmp_all(1,1))
         end if
         call onefile_gather_array_ncdio(tmp,tmp_all,dgrid)
         if(myrank == 0) then
            stat = nf_put_vara_real(ncid, id, start_id, count, tmp_all)
         end if
         deallocate(tmp_all)
#endif
      else
! ==============================================================================
#if !defined(MPI) || !defined(ONEFILE)
      stat = nf_put_vara_real(ncid, id, start, count, tmp)
#else
      if(myrank == 0) then
         allocate(tmp_all(dgrid%my%totalNx,dgrid%my%totalNy))
      else
         allocate(tmp_all(1,1))
      end if
      call onefile_gather_array_ncdio(tmp,tmp_all,dgrid)
       if(myrank == 0) then
         stat = nf_put_vara_real(ncid, id, start, count, tmp_all)
      end if
      deallocate(tmp_all)
#endif
! === Multiple rupture =========================================================
      end if
! ==============================================================================
#if defined(MPI) && defined(ONEFILE)
      if(myrank == 0) then
#endif
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
#if defined(MPI) && defined(ONEFILE)
      end if
#endif

      deallocate(tmp)

      return
   end subroutine write_array

#if !defined(MPI) || !defined(ONEFILE)
   subroutine close_file(dgrid)
#else
   subroutine close_file(dgrid, myrank)
#endif
      type(data_grids), target, intent(inout) :: dgrid
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4), intent(in) :: myrank
#endif
      integer(kind=4), pointer :: ncid

      ncid => dgrid%my%ncdio%ncid

#if defined(MPI) && defined(ONEFILE)
      if(myrank == 0) then
#endif
      stat = nf_close(ncid)
      if(stat /= NF_NOERR) write(0, '(a)') nf_strerror(stat)
#if defined(MPI) && defined(ONEFILE)
      end if
#endif

      return
   end subroutine close_file
#endif

end module mod_ncdio

