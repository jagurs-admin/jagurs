#include "real.h"
module mod_mygmt_gridio
! === For negative max. height =================================================
use mod_params, only : missing_value, fill_value, fill_value_is_nan
! ==============================================================================
use mod_params, only : id_cf17, grd_nan2zero
implicit none
include 'netcdf.inc'

contains

! === For negative max. height =================================================
!  subroutine mygmt_grdio_d(z,x0,x1,y0,y1,dx,dy,zmin,zmax,nx,ny,filename)
   subroutine mygmt_grdio_d(z,x0,x1,y0,y1,dx,dy,zmin,zmax,nx,ny,filename, &
                            formatid,flag_missing_value)
! ==============================================================================
      real(kind=REAL_BYTE), dimension(nx,ny), intent(in) :: z
      real(kind=8), intent(in) :: x0, x1, y0, y1, dx, dy
      real(kind=REAL_BYTE), intent(in) :: zmin, zmax
      integer(kind=4), intent(in) :: nx, ny
      character(len=512), intent(in) :: filename
! === For negative max. height =================================================
      logical, optional, intent(in) :: flag_missing_value
! ==============================================================================

      integer(kind=4) :: ncid, side_dim, xysize_dim ! netCDF id dimensions ids
      integer(kind=4) :: xysize_len, side_len ! 2-dimension lengths
      integer(kind=4) :: x_range_id, y_range_id, z_range_id ! variable ids
      integer(kind=4) :: spacing_id, z_id, dimension_id
      integer(kind=4) :: x_range_dims(2), y_range_dims(2), z_range_dims(2)
      integer(kind=4) :: spacing_dims(2), dimension_dims(2), z_dims(2)
      real(kind=8) :: z_scale_factor(1), z_add_offset(1)
      integer(kind=4) :: z_node_offset(1)
      integer(kind=4) :: stat

      !** set these ***
      real(kind=8), allocatable, dimension(:) :: x_range, y_range, z_range, spacing
      integer(kind=4), allocatable, dimension(:) :: dimension
! === Write buffer must be 4 byte. =============================================
      real(kind=4), allocatable, dimension(:,:) :: z_tmp
! ==============================================================================
      logical :: missing_value_is_available
      integer(kind=4), intent(in) :: formatid
      real(kind=8), allocatable, dimension(:) :: tmp
#if !defined(__NEC__) && !defined(__GFORTRAN__) && !defined(__amdflang__)
      real(kind=4), parameter :: NaN = transfer(Z'FFFFFFFF', 0.e0)
#else
      real(kind=4), parameter :: NaN = Z'FFFFFFFF'
#endif
      integer(kind=4) :: x_dim, y_dim, i
      integer(kind=4) :: j

      xysize_len = nx*ny
      side_len = 2 ! 2 sides

      allocate(x_range(side_len))
      allocate(y_range(side_len))
      allocate(z_range(side_len))
      allocate(spacing(side_len))
      allocate(dimension(side_len))

      if(formatid == nf_format_classic) then
      ! enter define mode
      stat = nf_create(filename, NF_CLOBBER, ncid)

      ! define 2 dimensions return dimension id's
      stat = nf_def_dim(ncid, 'side', side_len, side_dim)
      stat = nf_def_dim(ncid, 'xysize', xysize_len, xysize_dim)

      ! define variables
      x_range_dims(1) = side_dim
      y_range_dims(1) = side_dim
      z_range_dims(1) = side_dim
      spacing_dims(1) = side_dim
      dimension_dims(1) = side_dim
      z_dims(1) = xysize_dim

      stat = nf_def_var(ncid, 'x_range', NF_DOUBLE, 1, x_range_dims, x_range_id )
      stat = nf_def_var(ncid, 'y_range', NF_DOUBLE, 1, y_range_dims, y_range_id )
      stat = nf_def_var(ncid, 'z_range', NF_DOUBLE, 1, z_range_dims, z_range_id )
      stat = nf_def_var(ncid, 'spacing', NF_DOUBLE, 1, spacing_dims, spacing_id )
      stat = nf_def_var(ncid, 'dimension', NF_INT, 1, dimension_dims, dimension_id )
      stat = nf_def_var(ncid, 'z', NF_REAL, 1, z_dims, z_id)

      ! assign attributes
      stat = nf_put_att_text(ncid, x_range_id, 'units', 19, 'Longitude (Degrees)')
      stat = nf_put_att_text(ncid, y_range_id, 'units', 18, 'Latitude (Degrees)')
      stat = nf_put_att_text(ncid, z_range_id, 'units', 8, '(meters)')
      z_scale_factor(1) = 1
      stat = nf_put_att_double(ncid, z_id, 'scale_factor', NF_DOUBLE, 1, z_scale_factor)
      z_add_offset(1) = 0
      stat = nf_put_att_double(ncid, z_id, 'add_offset', NF_DOUBLE, 1, z_add_offset)
! === For negative max. height =================================================
      missing_value_is_available = .false.
      if(present(flag_missing_value)) then
         if(flag_missing_value) then
            missing_value_is_available = .true.
         end if
      end if

      if(missing_value_is_available) then
         stat = nf_put_att_real(ncid, z_id, '_FillValue', NF_REAL, 1, &
                                real(fill_value))
      end if
! ==============================================================================
      ! Burbidge - changed this to normal    z_node_offset[0] = 1;
      z_node_offset(1) = 0 ! 1=pixel 0=node registration
      stat = nf_put_att_int(ncid, z_id, 'node_offset', NF_INT, 1, z_node_offset)
      stat = nf_put_att_text(ncid, NF_GLOBAL, 'title', 29, 'tsunami linear FD calculation')
      stat = nf_put_att_text(ncid, NF_GLOBAL, 'source', 17, 'ctsunami_rwg-v1.4')

      ! leave define mode
      stat = nf_enddef(ncid)

      x_range(1) = dble(x0)
      x_range(2) = dble(x1)
      y_range(1) = dble(y0)
      y_range(2) = dble(y1)
      z_range(1) = dble(zmin)
      z_range(2) = dble(zmax)
      spacing(1) = dble(dx)
      spacing(2) = dble(dy)
      dimension(1) = nx
      dimension(2) = ny

      ! store x y and z range spacing
      stat = nf_put_var_double(ncid, x_range_id, x_range)
      stat = nf_put_var_double(ncid, y_range_id, y_range)
      stat = nf_put_var_double(ncid, z_range_id, z_range)
      stat = nf_put_var_double(ncid, spacing_id, spacing)
      ! store dimension
      stat = nf_put_var_int(ncid, dimension_id, dimension)
      else
         ! enter define mode
         stat = nf_create(filename, NF_NETCDF4, ncid)

         ! define 2 dimensions return dimension id's
         if(formatid /= id_cf17) then
            stat = nf_def_dim(ncid, 'x', nx, x_dim)
            stat = nf_def_dim(ncid, 'y', ny, y_dim)
         else
            stat = nf_def_dim(ncid, 'lon', nx, x_dim)
            stat = nf_def_dim(ncid, 'lat', ny, y_dim)
         end if

         ! define variables
         z_dims(1) = x_dim
         z_dims(2) = y_dim

         if(formatid /= id_cf17) then
            stat = nf_def_var(ncid, 'x', NF_DOUBLE, 1, (/x_dim/), x_range_id )
            stat = nf_def_var(ncid, 'y', NF_DOUBLE, 1, (/y_dim/), y_range_id )
         else
            stat = nf_def_var(ncid, 'lon', NF_DOUBLE, 1, (/x_dim/), x_range_id )
            stat = nf_def_var(ncid, 'lat', NF_DOUBLE, 1, (/y_dim/), y_range_id )
         end if
         stat = nf_def_var(ncid, 'z', NF_REAL, 2, z_dims, z_id )

         ! assign attributes
         if(formatid /= id_cf17) then
            stat = nf_put_att_text(ncid, x_range_id, 'long_name', 1, 'x')
         else
            stat = nf_put_att_text(ncid, x_range_id, 'long_name', 1, 'longitude')
         end if
         x_range(1) = dble(x0)
         x_range(2) = dble(x1)
         stat = nf_put_att_double(ncid, x_range_id, 'actual_range', NF_DOUBLE, 2, x_range)

         if(formatid /= id_cf17) then
            stat = nf_put_att_text(ncid, y_range_id, 'long_name', 1, 'y')
         else
            stat = nf_put_att_text(ncid, y_range_id, 'long_name', 1, 'latitude')
         end if
         y_range(1) = dble(y0)
         y_range(2) = dble(y1)
         stat = nf_put_att_double(ncid, y_range_id, 'actual_range', NF_DOUBLE, 2, y_range)

         stat = nf_put_att_text(ncid, z_id, 'long_name', 1, 'z')
!        stat = nf_put_att_real(ncid, z_id, '_FillValue', NF_REAL, 1, NaN)
! === For negative max. height =================================================
         missing_value_is_available = .false.
         if(present(flag_missing_value)) then
            if(flag_missing_value) then
               missing_value_is_available = .true.
            end if
         end if

         if(missing_value_is_available) then
            stat = nf_put_att_real(ncid, z_id, '_FillValue', NF_REAL, 1, &
                                   real(fill_value))
         else
            stat = nf_put_att_real(ncid, z_id, '_FillValue', NF_REAL, 1, NaN)
         end if
! ==============================================================================
         z_range(1) = dble(zmin)
         z_range(2) = dble(zmax)
         stat = nf_put_att_double(ncid, z_id, 'actual_range', NF_DOUBLE, 2, z_range)

         ! leave define mode
         stat = nf_enddef(ncid)

         ! store x y and z range spacing
         allocate(tmp(nx))
#ifdef USE_GPU
!$omp target teams distribute parallel do
#endif
         do i = 1, nx
            tmp(i) = dble(x0) + dble(dx)*(i - 1)
         end do
         stat = nf_put_var_double(ncid, x_range_id, tmp)
         deallocate(tmp)

         allocate(tmp(ny))
#ifdef USE_GPU
!$omp target teams distribute parallel do
#endif
         do i = 1, ny
            tmp(i) = dble(y0) + dble(dy)*(i - 1)
         end do
         stat = nf_put_var_double(ncid, y_range_id, tmp)
         deallocate(tmp)
      end if

      ! store z
! === Write buffer must be 4 byte. =============================================
      allocate(z_tmp(nx,ny))
      if(fill_value_is_nan == 0) then
         if(formatid == nf_format_classic) then
#ifndef USE_GPU
            z_tmp = z
#else
!$omp target teams distribute parallel do collapse(2) private(i)
            do j = lbound(z_tmp,2), ubound(z_tmp,2)
               do i = lbound(z_tmp,1), ubound(z_tmp,1)
                  z_tmp(i,j) = z(i,j)
               end do
            end do
#endif
         else
#ifdef USE_GPU
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = 1, ny
               do i = 1, nx
                  z_tmp(i,ny-j+1) = z(i,j)
               end do
            end do
         end if
      else
         if(formatid == nf_format_classic) then
#ifdef USE_GPU
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = lbound(z_tmp,2), ubound(z_tmp,2)
               do i = lbound(z_tmp,1), ubound(z_tmp,1)
                  if(z(i,j) == missing_value) then
                     z_tmp(i,j) = real(fill_value)
                  else
                     z_tmp(i,j) = z(i,j)
                  end if
               end do
            end do
         else
#ifdef USE_GPU
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = 1, ny
               do i = 1, nx
                  if(z(i,j) == missing_value) then
                     z_tmp(i,ny-j+1) = real(fill_value)
                  else
                     z_tmp(i,ny-j+1) = z(i,j)
                  end if
               end do
            end do
         end if
      end if
      stat = nf_put_var_real(ncid, z_id, z_tmp)
      deallocate(z_tmp)
! ==============================================================================
      stat = nf_close(ncid)

      deallocate(x_range)
      deallocate(y_range)
      deallocate(z_range)
      deallocate(spacing)
      deallocate(dimension)

      return
   end subroutine mygmt_grdio_d

   subroutine read_gmt_grd_hdr(infilename, nx, ny, dx, dy, west, east, south, north, zmin, zmax, nxorg, nyorg, formatid, read_header)
      use mod_params, only : RUDEF, IUDEF
#ifndef __GFORTRAN__
      character(len=256), intent(in) :: infilename
#else
      character(len=*), intent(in) :: infilename
#endif
      integer(kind=4), intent(out) :: nx, ny
      real(kind=REAL_BYTE), intent(out) :: dx, dy, west, east, south, north, zmin, zmax

      integer(kind=4), intent(out) :: nxorg, nyorg
#ifndef PIXELIN
      integer(kind=4) :: ncid, err
      integer(kind=4) :: side_id, xysize_id, side_len, xysize_len, side_dim, xysize_dim
      real(kind=8), allocatable, dimension(:) :: x_range, y_range, z_range, spacing
      integer(kind=4), allocatable, dimension(:) :: dimension
      integer(kind=4) :: x_range_id, y_range_id, z_range_id, spacing_id, dimension_id
#else
      character(len=256) :: descfile
      namelist /desc/ west, east, south, north, dx, dy, zmin, zmax, nx, ny
      real(kind=REAL_BYTE) :: westorg, eastorg, southorg, northorg
#endif
      integer(kind=4), intent(out) :: formatid
      integer(kind=4) :: node_offset
      logical, intent(in), optional :: read_header
      logical :: check_offset

      !*** open the file and inquire about the dimensions ***
#ifndef PIXELIN
      err = nf_open(infilename, NF_NOWRITE, ncid)
      if(err /= NF_NOERR) then
         write(0,'(a,i0,a,a)') 'netcdf err=', err, ' nonexistent file=', trim(infilename)
         stop
      else if(err == NF_ENOTNC) then
         write(0,'(a,i0,a,a)') 'netcdf err=', err, ' non-netCDF file=', trim(infilename)
         stop
      end if

      err = nf_inq_format(ncid, formatid)

      if(formatid == nf_format_classic) then
      err = nf_inq_dimid(ncid, 'side', side_id)
      err = nf_inq_dimid(ncid, 'xysize', xysize_id)

      !** get the dimensions ***
      err = nf_inq_dimlen(ncid, side_id, side_len)
      err = nf_inq_dimlen(ncid, xysize_id, xysize_len)
      side_dim = int(side_len)
      xysize_dim = int(xysize_len)
      else
         side_dim = 2
      end if

      !*** allocate space for variables ***
      allocate(x_range(side_dim))
      allocate(y_range(side_dim))
      allocate(z_range(side_dim))
      allocate(spacing(side_dim))
      allocate(dimension(side_dim))

      if(formatid == nf_format_classic) then
      !*** get the variable id's - these variable names are common in GMT ***
      err = nf_inq_varid(ncid, 'x_range', x_range_id)
      err = nf_inq_varid(ncid, 'y_range', y_range_id)
      err = nf_inq_varid(ncid, 'z_range', z_range_id)
      err = nf_inq_varid(ncid, 'spacing', spacing_id)
      err = nf_inq_varid(ncid, 'dimension', dimension_id)

      !*** get the variables ***
      err = nf_get_var_double(ncid, x_range_id, x_range)
      err = nf_get_var_double(ncid, y_range_id, y_range)
      err = nf_get_var_double(ncid, z_range_id, z_range)
      err = nf_get_var_double(ncid, spacing_id, spacing)
      err = nf_get_var_int(ncid, dimension_id, dimension)
      else
         check_offset = .true.
         if(present(read_header)) then
            if(read_header) check_offset = .false.
         end if

         err = nf_get_att_int(ncid, nf_global, 'node_offset', node_offset)
         if((node_offset == 1) .and. (err == 0) .and. check_offset) then
            write(0,'(a)') 'Error! Grid file with pixel node registration is NOT available!!!'
            stop
         end if

         err = nf_inq_varid(ncid, 'x', x_range_id)
         if(err == nf_enotvar) then
            err = nf_inq_varid(ncid, 'lon', x_range_id)
            if(err == nf_enotvar) then
               write(0,'(a,a,a)') 'Error! The format of GMT file "', trim(infilename), '" is invalid!'
               stop
            end if
            formatid = id_cf17
         end if

         if(formatid /= id_cf17) then
            err = nf_get_att_double(ncid, x_range_id, 'actual_range', x_range)

            err = nf_inq_varid(ncid, 'y', y_range_id)
            err = nf_get_att_double(ncid, y_range_id, 'actual_range', y_range)

            err = nf_inq_varid(ncid, 'z', z_range_id)
            err = nf_get_att_double(ncid, z_range_id, 'actual_range', z_range)

            err = nf_inq_dimid(ncid, 'x', x_range_id)
            err = nf_inq_dimlen(ncid, x_range_id, dimension(1))

            err = nf_inq_dimid(ncid, 'y', y_range_id)
            err = nf_inq_dimlen(ncid, y_range_id, dimension(2))
         else
            err = nf_get_att_double(ncid, x_range_id, 'actual_range', x_range)

            err = nf_inq_varid(ncid, 'lat', y_range_id)
            err = nf_get_att_double(ncid, y_range_id, 'actual_range', y_range)

            err = nf_inq_varid(ncid, 'z', z_range_id)
            err = nf_get_att_double(ncid, z_range_id, 'actual_range', z_range)

            err = nf_inq_dimid(ncid, 'lon', x_range_id)
            err = nf_inq_dimlen(ncid, x_range_id, dimension(1))

            err = nf_inq_dimid(ncid, 'lat', y_range_id)
            err = nf_inq_dimlen(ncid, y_range_id, dimension(2))
         end if

         spacing(1) = (x_range(2) - x_range(1))/(dimension(1) - 1)
         spacing(2) = (y_range(2) - y_range(1))/(dimension(2) - 1)
         !spacing(1) = (x_range(2) - x_range(1))/dimension(1)
         !spacing(2) = (y_range(2) - y_range(1))/dimension(2)
      end if
      err = nf_close(ncid)

      !*** move values into return variables ***
#ifndef REAL_DBLE
      west  = x_range(1)
      east  = x_range(2)
      south = y_range(1)
      north = y_range(2)
      zmin  = z_range(1)
      zmax  = z_range(2)
      dx = spacing(1)
      dy = spacing(2)
      nx = dimension(1)
      ny = dimension(2)
#else
! === To avoid numerical error on MPI decomposition. ===========================
      west  = real(x_range(1))
      east  = real(x_range(2))
      south = real(y_range(1))
      north = real(y_range(2))
      zmin  = real(z_range(1))
      zmax  = real(z_range(2))
      dx = real(spacing(1))
      dy = real(spacing(2))
      nx = real(dimension(1))
      ny = real(dimension(2))
! ==============================================================================
#endif

      !*** free space ***/
      deallocate(x_range)
      deallocate(y_range)
      deallocate(z_range)
      deallocate(spacing)
      deallocate(dimension)
#else
      descfile = trim(infilename) // '.desc'
      open(1,file=trim(descfile),action='read',status='old',form='formatted',err=100)
      west = RUDEF; east = RUDEF; south = RUDEF; north = RUDEF
      dx   = RUDEF; dy   = RUDEF; nx    = IUDEF; ny    = IUDEF
      read(1,desc)
      if((west == RUDEF) .or. (east == RUDEF) .or. (south == RUDEF) .or. (north == RUDEF) .or. &
         (dx   == RUDEF) .or. (dy   == RUDEF) .or. (nx    == IUDEF) .or. (ny    == IUDEF)) goto 101
      nxorg = nx
      nyorg = ny
      nx = (nxorg-1)/3; nx = 3*nx+1
      ny = (nyorg-1)/3; ny = 3*ny+1
      westorg = west; eastorg = east; northorg = north; southorg = south
      west  = west  + 1.5d0*dx
      east  = east  - 1.5d0*dx
      north = north - 1.5d0*dy
      south = south + 1.5d0*dy
      close(1)
      write(6,'(a,a,a)') '[pixcel format] 1 grid on each edge is ignored on "', trim(infilename), """!"
      write(6,'(a,i6,a,i6)') '[pixcel format] Grid size is dealt as ', nx, ' x ', ny
      write(6,'(a,i6,a,i6,a)') '[pixcel format]    instead of it''s original gird size ', nxorg, ' x ', nyorg, '.'
      write(6,'(a,f0.3,a,f0.3,a)') '[pixcel format] West edge is changed from ', westorg, ' to ', west, '.'
      write(6,'(a,f0.3,a,f0.3,a)') '[pixcel format] East edge is changed from ', eastorg, ' to ', east, '.'
      write(6,'(a,f0.3,a,f0.3,a)') '[pixcel format] South edge is changed from ', southorg, ' to ', south, '.'
      write(6,'(a,f0.3,a,f0.3,a)') '[pixcel format] North edge is changed from ', northorg, ' to ', north, '.'
#ifdef PIXCELOUT
      write(6,'(a)') '[pixcel format] Note that the value on truncated area is zero!'
#endif
#endif

      return
#ifdef PIXELIN
 100  continue
      write(0,'(a,a,a)') 'Error! Descriptor file "', trim(descfile), '" does not exist!'
      stop
 101  continue
      write(0,'(a,a,a)') 'Error! Invalid description on descriptor file "', trim(descfile), '"!'
      stop
#endif
   end subroutine read_gmt_grd_hdr

   subroutine read_gmt_grd(infilename,z,nx,ny,formatid)
#ifndef __GFORTRAN__
      character(len=256), intent(in) :: infilename
#else
      character(len=*), intent(in) :: infilename
#endif
      real(kind=REAL_BYTE), dimension(nx,ny), intent(inout) :: z
!     real(kind=REAL_BYTE), intent(inout) :: z
      integer(kind=4), intent(in) :: nx, ny

      integer(kind=4) :: err, ncid, z_id
! === Read buffer must be 4 byte. ==============================================
      real(kind=4), allocatable, dimension(:,:) :: z_tmp
! ==============================================================================
      integer(kind=4), intent(in) :: formatid
      integer(kind=4) :: i, j, ierr

      !*** open file again ***
      err = nf_open(infilename, NF_NOWRITE, ncid)
      if(err /= NF_NOERR) then
         write(0,'(a,i0,a,a)') 'netcdf err=', err, ' nonexistent file=', trim(infilename)
         stop
      else if(err == NF_ENOTNC) then
         write(0,'(a,i0,a,a)') 'netcdf err=', err, ' non-netCDF file=', trim(infilename)
         stop
      end if

      !*** inquire about z-data and read into float array
      !    assumes space is already allocated ***
      err = nf_inq_varid(ncid, 'z', z_id)
! === Read buffer must be 4 byte. ==============================================
      allocate(z_tmp(nx,ny))
      err = nf_get_var_real(ncid, z_id, z_tmp)
      ierr = 0
      if(formatid == nf_format_classic) then
         ! NaN is modified into zero.
#ifndef USE_GPU
#ifndef __NEC__
!$omp parallel do private(i)
#else
!$omp parallel do private(i) reduction(+:ierr)
#endif
#else
!$omp target teams distribute parallel do collapse(2) private(i) reduction(+:ierr)
#endif
         do j = lbound(z,2), ubound(z,2)
            do i = lbound(z,1), ubound(z,1)
               if(z_tmp(i,j) /= z_tmp(i,j)) then
                  z_tmp(i,j) = 0.0d0
#if !defined(__NEC__) && !defined(USE_GPU)
!$omp critical
                  ierr = 1
!$omp end critical
#else
                  ierr = ierr + 1
#endif
               end if
               z(i,j) = z_tmp(i,j)
            end do
         end do
      else
         ! NaN is modified into zero.
#ifndef USE_GPU
#ifndef __NEC__
!$omp parallel do private(i)
#else
!$omp parallel do private(i) reduction(+:ierr)
#endif
#else
!$omp target teams distribute parallel do collapse(2) private(i) reduction(+:ierr)
#endif
         do j = 1, ny
            do i = 1, nx
               if(z_tmp(i,j) /= z_tmp(i,j)) then
                  z_tmp(i,j) = 0.0d0
#if !defined(__NEC__) && !defined(USE_GPU)
!$omp critical
                  ierr = 1
!$omp end critical
#else
                  ierr = ierr + 1
#endif
               end if
               z(i,ny-j+1) = z_tmp(i,j)
            end do
         end do
      end if
      deallocate(z_tmp)

      if((ierr /= 0) .and. (grd_nan2zero == 0)) then
         write(0,'(a)')     '======================================================='
         write(0,'(a)')     '======================================================='
         write(0,'(a)')     '=== ERROR!!!'
         write(0,'(a,a,a)') '=== Input file "', trim(infilename), '" includes NaN!'
         write(0,'(a)')     '=== Check the file or specify "grd_nan2zero = 1"'
         write(0,'(a)')     '=== to replace NaNs into zeros.'
         write(0,'(a)')     '======================================================='
         write(0,'(a)')     '======================================================='
         stop
      end if
! ==============================================================================

      !*** close file and return pointer to array ***
      err = nf_close(ncid)

      return
   end subroutine read_gmt_grd

end module mod_mygmt_gridio
