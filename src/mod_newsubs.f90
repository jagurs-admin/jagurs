#include "real.h"
module mod_newsubs
use mod_grid
use mod_rwg
use mod_mygmt_gridio, only : read_gmt_grd_hdr, read_gmt_grd, mygmt_grdio_d
#ifdef MPI
use mod_mpi_fatal
! === Support multiple ruptures. ===============================================
use mod_params, only : multrupt
! ==============================================================================
#endif
! === For ensemble =============================================================
#ifdef MULTI
use mod_multi, only : members_dir
#endif
! ==============================================================================
implicit none

contains

#ifdef MPI
   subroutine add_suffix_to_grid_filenames(grids, ngrid, suffix)
      type(data_grids), target, dimension(ngrid), intent(inout) :: grids
      integer(kind=4), intent(in) :: ngrid
      character(len=128), intent(in) :: suffix

      integer(kind=4) :: slen = 0
      integer(kind=4) :: flen = 0
      integer(kind=4) :: maxlen = 128
      character(len=256), pointer :: filename
      integer(kind=4) :: ig

      slen = len_trim(suffix)

      do ig = 1, ngrid
         ! bath_file
         filename => grids(ig)%my%bath_file
         flen = len_trim(filename)
         if(flen + slen + 1 > maxlen) then
             write(0,'(a,a,a,a,a,i0,a)') 'add_suffix_to_grid_filenames(newsubs.o): too long filename with suffix ( file ', &
                trim(filename), ' suffix ', trim(suffix), ' max ', maxlen, ' characters )'
             call fatal_error(220)
         end if
         filename = trim(filename) // trim(suffix)

         ! disp_file
! === Support multiple ruptures. ===============================================
         if(multrupt /= 1) then
! ==============================================================================
         filename => grids(ig)%my%disp_file
         if(trim(filename) /= 'NO_DISPLACEMENT_FILE_GIVEN') then
            flen = len_trim(filename)
            if(flen + slen + 1 > maxlen) then
                write(0,'(a,a,a,a,a,i0,a)') 'add_suffix_to_grid_filenames(newsubs.o): too long filename with suffix ( file ', &
                   trim(filename), ' suffix ', trim(suffix), ' max ', maxlen, ' characters )'
                call fatal_error(221)
            end if
            filename = trim(filename) // trim(suffix)
         end if
! === Support multiple ruptures. ===============================================
         end if
! ==============================================================================

         ! wod_file
         filename => grids(ig)%wod_file
         if(trim(filename) /= 'NO_WETORDRY_FILE_GIVEN') then
            flen = len_trim(filename)
            if(flen + slen + 1 > maxlen) then
                write(0,'(a,a,a,a,a,i0,a)') 'add_suffix_to_grid_filenames(newsubs.o): too long filename with suffix ( file ', &
                   trim(filename), ' suffix ', trim(suffix), ' max ', maxlen, ' characters )'
                call fatal_error(222)
            end if
            filename = trim(filename) // trim(suffix)
         end if

         ! bcf_file
         filename => grids(ig)%bcf_file
         if(trim(filename) /= 'NO_FRICTION_FILE_GIVEN') then
            flen = len_trim(filename)
            if(flen + slen + 1 > maxlen) then
                write(0,'(a,a,a,a,a,i0,a)') 'add_suffix_to_grid_filenames(newsubs.o): too long filename with suffix ( file ', &
                   trim(filename), ' suffix ', trim(suffix), ' max ', maxlen, ' characters )'
                call fatal_error(223)
            end if
            filename = trim(filename) // trim(suffix)
         end if
      end do

      return
   end subroutine add_suffix_to_grid_filenames
#endif
   subroutine read_bathymetry_gmt_grdhdr(fname,nx,ny,dxdy,mlon0,mlat0)
      character(len=256), intent(in) :: fname
      integer(kind=4), intent(out) :: nx, ny
      real(kind=REAL_BYTE), intent(out) :: dxdy, mlon0, mlat0
      integer(kind=4) :: nlon, nlat
      real(kind=REAL_BYTE) :: dx, dy, west, east, south, north, zmin, zmax

      call read_gmt_grd_hdr(fname,nlon,nlat,dx,dy,west,east,south,north,zmin,zmax)

      write(6,'(/,8x,a,a)') 'read_bathymetry_gmt_grdhdr(newsub.o): file name=', trim(fname)
      write(6,'(8x,a,i0,a,i0)') 'nx=', nlon, ' ny=', nlat
      write(6,'(8x,a,f0.6,a,f0.6)') 'dx=', dx, ' dy=', dy
      write(6,'(8x,a,f0.3,a,f0.3)') 'zmin=', zmin, ' zmax=', zmax
      write(6,'(8x,a,f0.3,a,f0.3,a,f0.3,a,f0.3)') 'west=', west, ' east=', east, ' south=', south, ' north=', north

      nx = nlon
      ny = nlat

      dxdy = dx ! assumes dx=dy

      if(abs(dx - dy) >= 0.00001d0) then
         write(0,'(a)') 'read_bathymetry_gmt_grdhdr(newsubs.o): dxdy not same!'
#ifndef MPI
         stop
#else
         call fatal_error(204)
#endif
      end if

#ifndef CARTESIAN
      !*** colongitude ***
      if(west < 0.0d0) then
         mlon0 = (360.0d0+west)*60.0d0
      else
         mlon0 = west*60.0d0
      end if

      !*** colatitude ***
      mlat0 = (90.0d0-north)*60.0d0
#else
      mlon0 = west
      mlat0 = south
#endif

      write(6,'(8x,a,f0.3,a,f0.3,a,f0.6)') 'mlon0=', mlon0, ' mlat0=', mlat0, ' dxdy=', dxdy

      return
   end subroutine read_bathymetry_gmt_grdhdr

   subroutine read_friction_gmt_grd(fname,ffld,nx,ny)
      character(len=256), intent(in) :: fname
      real(kind=REAL_BYTE), target, dimension(nx,ny), intent(inout) :: ffld
      integer(kind=4), intent(in) :: nx, ny

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: bcf

      bcf => ffld

      if(trim(fname) == 'NO_FRICTION_FILE_GIVEN') then
         write(6,'(8x,a)') trim(fname)
         bcf = 0.0d0
      else
         write(6,'(8x,a,a)') 'FRICTION_FILE_GIVEN:', trim(fname)
         call read_gmt_grd(fname,bcf,nx,ny)
      end if

      return
   end subroutine read_friction_gmt_grd

   subroutine read_bathymetry_gmt_grd(fname,dfld,nx,ny,linear)
      character(len=256), intent(in) :: fname
      type(depth_arrays), target, intent(inout) :: dfld
      integer(kind=4), intent(in) :: nx, ny, linear
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy, dz
      real(kind=REAL_BYTE) :: zap = 0.0d0, zmin, zmax
#ifndef MPI
      integer(kind=4) ::  imin, jmin, imax, jmax
#endif
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dz_tmp

      dx => dfld%dx
      dy => dfld%dy
      dz => dfld%dz

      allocate(dz_tmp(nx,ny))
      call read_gmt_grd(fname,dz_tmp,nx,ny)
      do j = 1, ny
         do i = 1, nx
            dz(i,j) = dz_tmp(i,j)
         end do
      end do
#ifndef MPI
      call minmax_rwg(nx,ny,dz_tmp,zmax,zmin,imin,jmin,imax,jmax)
#else
      call minmax_rwg(nx,ny,dz_tmp,zmax,zmin)
#endif
      deallocate(dz_tmp)
      write(6,'(8x,a)') 'read_bathymetry_gmt_grd(newsubs.o):'
#ifndef MPI
      write(6,'(8x,a,f0.3,a,i0,a,i0,a,f0.3,a,i0,a,i0)') &
         'min=', zmin, ' imin=', imin, ' jmin=', jmin, &
         ' max=', zmax, ' imax=', imax, ' jmax=', jmax
#else
      write(6,'(8x,a,f0.3,a,f0.3)') 'min=', zmin, ' max=', zmax
#endif

      if(linear == 1) then
         do j = 1, ny
            do i = 1, nx-1
               if(dz(i+1,j) < zap .or. dz(i,j) < zap) then
                  dx(i,j) = zap
               else
                  dx(i,j) = 0.5d0*(dz(i+1,j)+dz(i,j))
               end if
            end do
            if(dz(nx,j) < zap) then
               dx(nx,j) = zap
            else 
               dx(nx,j) = dz(nx,j)
            end if
         end do
         do i = 1, nx
            do j = 1, ny-1
               if(dz(i,j+1) < zap .or. dz(i,j) < zap) then
                  dy(i,j) = zap
               else
                  dy(i,j) = 0.5d0*(dz(i,j+1)+dz(i,j))
               end if
            end do
            if(dz(i,ny) < zap) then
               dy(i,ny) = zap
            else
               dy(i,ny) = dz(i,ny)
            end if
         end do
      else
         do j = 1, ny
            do i = 1, nx-1
               dx(i,j) = 0.5d0*(dz(i+1,j)+dz(i,j))
            end do
            dx(nx,j) = dz(nx,j)
         end do
         do i = 1, nx
            do j = 1, ny-1
               dy(i,j) = 0.5d0*(dz(i,j+1)+dz(i,j))
            end do
            dy(i,ny) = dz(i,ny)
         end do
      end if

      return
   end subroutine read_bathymetry_gmt_grd

#ifndef DIROUT
! === For negative max. height =================================================
!  subroutine maxgrd_write_gmt(hzmax,nlon,nlat,mlon0,mlat0,dxdy,fname)
   subroutine maxgrd_write_gmt(hzmax,nlon,nlat,mlon0,mlat0,dxdy,fname, &
                               flag_missing_value)
! ==============================================================================
#else
! === For negative max. height =================================================
!  subroutine maxgrd_write_gmt(hzmax,nlon,nlat,mlon0,mlat0,dxdy,dirname,fname)
   subroutine maxgrd_write_gmt(hzmax,nlon,nlat,mlon0,mlat0,dxdy,dirname,fname, &
                               flag_missing_value)
! ==============================================================================
#endif
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(inout) :: hzmax
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: mlon0, mlat0, dxdy
      character(len=512), intent(in) :: fname
! === For negative max. height =================================================
      logical, intent(in) :: flag_missing_value
! ==============================================================================

      real(kind=REAL_BYTE) :: zmax, zmin
#ifndef MPI
      integer(kind=4) :: imax, jmax, imin, jmin
#endif
      real(kind=8) :: lat_north, lat_south, lon_east, lon_west, dx, dy
#ifdef DIROUT
      character(len=128), intent(in) :: dirname
      character(len=128) :: command = ''
      character(len=512) :: fname_dir = ''
      write(command,'(a,a)') 'mkdir -p ', trim(dirname)
      call system(trim(command))
#endif

      !*** calculate GMT grdinfo -R{west}/{east}/{south}/{west}
      !    corners -I{dx}/{dy} -N{nlon}/{nlat} ***/
#ifndef MPI
! === For negative max. height =================================================
!     call minmax_rwg(nlon,nlat,hzmax,zmax,zmin,imin,jmin,imax,jmax)
      call minmax_rwg(nlon,nlat,hzmax,zmax,zmin,imin,jmin,imax,jmax, &
                      flag_missing_value)
! ==============================================================================
#else
! === For negative max. height =================================================
!     call minmax_rwg(nlon,nlat,hzmax,zmax,zmin)
      call minmax_rwg(nlon,nlat,hzmax,zmax,zmin,flag_missing_value)
! ==============================================================================
#endif
#ifndef CARTESIAN
      lat_north = 90.0d0 - mlat0/60.0d0
      lon_west  = mlon0/60.0d0
      ! Burbidge - Changed the following to be in normal node registration
      ! NB 1.0 is used instead of 1 to improve output precision to GMT tolerable
      lat_south = 90.0d0 - (mlat0/60.0d0 + (nlat-1.0d0)*dxdy)
      lon_east  = mlon0/60.0d0 + (nlon-1.0d0)*dxdy
! === DEBUG by tkato 2015/02/02 ================================================
!     if(lon_west > 180.0d0) lon_west = lon_west - 360.0d0
! ==============================================================================
#else
      lon_west  = mlon0
      lat_south = mlat0
      lon_east  = lon_west  + (nlon-1.0d0)*dxdy
      lat_north = lat_south + (nlat-1.0d0)*dxdy
#endif
      dx = dxdy
      dy = dx
! === To add max velocity output. by tkato 2012/10/02 ==========================
      write(6,'(a)') 'Max wave hight'
! ==============================================================================
#ifndef CARTESIAN
      write(6,'(a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i0,a,i0,a,a)') &
#else
      write(6,'(a,f9.3,a,f9.3,a,f0.1,a,f0.1,a,f0.1,a,f0.1,a,i0,a,i0,a,a)') &
#endif
         'min=', zmin, ' max =', zmax, ' (', lon_west, '/', lon_east, '/', &
         lat_south, '/', lat_north, ') nx=', nlon, ' ny=', nlat, ' ', trim(fname)

#ifdef DIROUT
      fname_dir = trim(dirname) // '/' // trim(fname)
#endif
      !*** write a GMT grid file using libnetcdf.a routines ***
      call mygmt_grdio_d(hzmax, lon_west, lon_east, lat_south, lat_north, &
#ifndef DIROUT
! === For negative max. height =================================================
!                        dx, dy, zmin, zmax, nlon, nlat, fname)
                         dx, dy, zmin, zmax, nlon, nlat, fname, flag_missing_value)
! ==============================================================================
#else
! === For negative max. height =================================================
!                        dx, dy, zmin, zmax, nlon, nlat, fname_dir)
                         dx, dy, zmin, zmax, nlon, nlat, fname_dir, flag_missing_value)
! ==============================================================================
#endif

      return
   end subroutine maxgrd_write_gmt
! === To add max velocity output. by tkato 2012/10/02 ==========================
#ifndef DIROUT
   subroutine maxgrd_v_write_gmt(vmax,nlon,nlat,mlon0,mlat0,dxdy,fname)
#else
   subroutine maxgrd_v_write_gmt(vmax,nlon,nlat,mlon0,mlat0,dxdy,dirname,fname)
#endif
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(inout) :: vmax
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: mlon0, mlat0, dxdy
      character(len=512), intent(in) :: fname

      real(kind=REAL_BYTE) :: zmax, zmin
#ifndef MPI
      integer(kind=4) :: imax, jmax, imin, jmin
#endif
      real(kind=8) :: lat_north, lat_south, lon_east, lon_west, dx, dy
#ifdef DIROUT
      character(len=128), intent(in) :: dirname
      character(len=128) :: command = ''
      character(len=512) :: fname_dir = ''
      write(command,'(a,a)') 'mkdir -p ', trim(dirname)
      call system(trim(command))
#endif

      !*** calculate GMT grdinfo -R{west}/{east}/{south}/{west}
      !    corners -I{dx}/{dy} -N{nlon}/{nlat} ***/
#ifndef MPI
      call minmax_rwg(nlon,nlat,vmax,zmax,zmin,imin,jmin,imax,jmax)
#else
      call minmax_rwg(nlon,nlat,vmax,zmax,zmin)
#endif
#ifndef CARTESIAN
      lat_north = 90.0d0 - mlat0/60.0d0
      lon_west  = mlon0/60.0d0
      ! Burbidge - Changed the following to be in normal node registration
      ! NB 1.0 is used instead of 1 to improve output precision to GMT tolerable
      lat_south = 90.0d0 - (mlat0/60.0d0 + (nlat-1.0d0)*dxdy)
      lon_east  = mlon0/60.0d0 + (nlon-1.0d0)*dxdy
! === DEBUG by tkato 2015/02/02 ================================================
!     if(lon_west > 180.0d0) lon_west = lon_west - 360.0d0
! ==============================================================================
#else
      lon_west  = mlon0
      lat_south = mlat0
      lon_east  = lon_west  + (nlon-1.0d0)*dxdy
      lat_north = lat_south + (nlat-1.0d0)*dxdy
#endif
      dx = dxdy
      dy = dx
      write(6,'(a)') 'Max velocity'
#ifndef CARTESIAN
      write(6,'(a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i0,a,i0,a,a)') &
#else
      write(6,'(a,f9.3,a,f9.3,a,f0.1,a,f0.1,a,f0.1,a,f0.1,a,i0,a,i0,a,a)') &
#endif
         'min=', zmin, ' max =', zmax, ' (', lon_west, '/', lon_east, '/', &
         lat_south, '/', lat_north, ') nx=', nlon, ' ny=', nlat, ' ', trim(fname)

#ifdef DIROUT
      fname_dir = trim(dirname) // '/' // trim(fname)
#endif
      !*** write a GMT grid file using libnetcdf.a routines ***
      call mygmt_grdio_d(vmax, lon_west, lon_east, lat_south, lat_north, &
#ifndef DIROUT
                         dx, dy, zmin, zmax, nlon, nlat, fname)
#else
                         dx, dy, zmin, zmax, nlon, nlat, fname_dir)
#endif

      return
   end subroutine maxgrd_v_write_gmt
! ==============================================================================

#ifndef MPI
   subroutine dump_gmt_nl(wfld,dfld,tfld,nlon,nlat,wod, &
! === Speed output. ============================================================
!                         mlat0,mlon0,dxdy,t,istep,base,velgrd_flag,mode)
                          mlat0,mlon0,dxdy,t,istep,base,mode)
! ==============================================================================
#else
   subroutine dump_gmt_nl(wfld,dfld,tfld,nlon,nlat,wod, &
! === Speed output. ============================================================
!                         mlat0,mlon0,dxdy,t,istep,myrank,base,velgrd_flag,mode,bflag)
                          mlat0,mlon0,dxdy,t,istep,myrank,base,mode,bflag)
! ==============================================================================
#endif
! === Speed output. ============================================================
!     use mod_params, only : VEL, HGT
      use mod_params, only : VEL, HGT, velgrd_flag, speedgrd_flag
! ==============================================================================
      type(wave_arrays), target, intent(in) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
      real(kind=REAL_BYTE), target, dimension(nlon,nlat), intent(in) :: tfld
      integer(kind=4), intent(in) :: nlon, nlat
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(in) :: wod
#endif
! ==============================================================================
      real(kind=REAL_BYTE), intent(in) :: mlat0, mlon0, dxdy, t
      integer(kind=4), intent(in) :: istep
      character(len=256), intent(in) :: base
#ifdef MPI
      integer(kind=4), intent(in) :: myrank
#endif
! === Conversion from flux to velocity should be done right after calc. ========
      integer(kind=4), intent(in) :: mode
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
#endif
! ==============================================================================

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz, tp, fx, fy
      character(len=512) :: fname

      real(kind=REAL_BYTE) :: zmax, zmin
      real(kind=8) :: lat_north, lat_south, lon_east, lon_west, dx, dy
#ifndef MPI
      integer(kind=4) :: imax, jmax, imin, jmin
#endif
      integer(kind=4) :: i, j

      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
! === Flux -> Flow speed =======================================================
! === Conversion from flux to velocity should be done right after calc. ========
!     real(kind=REAL_BYTE) :: td
      real(kind=REAL_BYTE) :: tdxm, tdxp, tdym, tdyp
      integer(kind=4) :: im, jm, ip, jp
! ==============================================================================
      real(kind=REAL_BYTE), parameter :: td_min = 0.01d0 ! 1 cm
! ==============================================================================
! === Speed output. ============================================================
      real(kind=REAL_BYTE) :: tx, ty, speed
! ==============================================================================
#ifdef DIROUT
      character(len=128) :: dirname = '', command = ''
      write(dirname,'(i8.8,a)') istep, '.grd'
! === For ensemble =============================================================
#ifdef MULTI
      dirname = trim(members_dir) // trim(dirname)
#endif
! ==============================================================================
      write(command,'(a,a)') 'mkdir -p ', trim(dirname)
      call system(trim(command))
#endif
! === Conversion from flux to velocity should be done right after calc. ========
      hz => wfld%hz
      dz => dfld%dz
      tp => tfld

#ifndef CARTESIAN
      lat_north = 90.0d0 - mlat0/60.0d0
      lon_west  = mlon0/60.0d0
      lat_south = 90.0d0 - (mlat0/60.0d0 + (nlat-1.0d0)*dxdy)
      lon_east  = mlon0/60.0d0 + (nlon-1.0d0)*dxdy
! === DEBUG by tkato 2015/02/02 ================================================
!     if(lon_west > 180.0d0) lon_west = lon_west - 360.0d0
! ==============================================================================
! === DEBUG by tkato. 2015/01/23 ===============================================
!     if(lon_east > 180.0d0) lon_east = lon_east - 360.0d0
! ==============================================================================
#else
      lon_west  = mlon0
      lat_south = mlat0
      lon_east  = lon_west  + (nlon-1.0d0)*dxdy
      lat_north = lat_south + (nlat-1.0d0)*dxdy
#endif
      dx = dxdy
      dy = dx
! ==============================================================================

! === Conversion from flux to velocity should be done right after calc. ========
      if(mode == HGT) then
! ==============================================================================
#ifndef MPI
! === To support over 1M steps. by tkato 2012/11/27 ============================
!     write(fname,'(a,a,i6.6,a)') trim(base), '.', istep, '.grd'
      write(fname,'(a,a,i8.8,a)') trim(base), '.', istep, '.grd'
! ==============================================================================
#else
! === To support over 1M steps. by tkato 2012/11/27 ============================
!     write(fname,'(a,a,i6.6,a,i6.6)') trim(base), '.', istep, '.grd.', myrank
      write(fname,'(a,a,i8.8,a,i6.6)') trim(base), '.', istep, '.grd.', myrank
! ==============================================================================
#endif
! === For ensemble =============================================================
#ifndef DIROUT
#ifdef MULTI
      fname = trim(members_dir) // trim(fname)
#endif
#endif
! ==============================================================================

! === Conversion from flux to velocity should be done right after calc. ========
!     hz => wfld%hz
!     dz => dfld%dz
!     tp => tfld
! ==============================================================================

      ! check for wet-or-dry
      do j = 1, nlat
         do i = 1, nlon
            if(wod(i,j) == 1) then
               tp(i,j) = hz(i,j)
            else
               tp(i,j) = zap
            end if
         end do
      end do

      !*** GMT -R{west}/{east}/{south}/{west} corners -I{dx}/{dy} -N{nlon}/{nlat} ***
#ifndef MPI
      call minmax_rwg(nlon,nlat,tp,zmax,zmin,imin,jmin,imax,jmax)
#else
      call minmax_rwg(nlon,nlat,tp,zmax,zmin)
#endif
! === Conversion from flux to velocity should be done right after calc. ========
!#ifndef CARTESIAN
!     lat_north = 90.0d0 - mlat0/60.0d0
!     lon_west  = mlon0/60.0d0
!     lat_south = 90.0d0 - (mlat0/60.0d0 + (nlat-1.0d0)*dxdy)
!     lon_east  = mlon0/60.0d0 + (nlon-1.0d0)*dxdy
!     if(lon_west > 180.0d0) lon_west = lon_west - 360.0d0
!     if(lon_east > 180.0d0) lon_east = lon_east - 360.0d0
!#else
!     lon_west  = mlon0
!     lat_south = mlat0
!     lon_east  = lon_west  + (nlon-1.0d0)*dxdy
!     lat_north = lat_south + (nlat-1.0d0)*dxdy
!#endif
!     dx = dxdy
!     dy = dx
! ==============================================================================

! === To support over 1M steps. by tkato 2012/11/27 ============================
!     write(6,'(a,i6.6,a,f9.2,a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i4,a,i4,a,a)') &
#ifndef CARTESIAN
      write(6,'(a,i8.8,a,f9.2,a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i4,a,i4,a,a)') &
#else
      write(6,'(a,i8.8,a,f9.2,a,f9.3,a,f9.3,a,f0.1,a,f0.1,a,f0.1,a,f0.1,a,i4,a,i4,a,a)') &
#endif
! ==============================================================================
         'it=', istep, ' t=', t, ' min=', zmin, ' max=', zmax, ' (', &
         lon_west, '/', lon_east, '/', lat_south, '/', lat_north, ') nx=', &
         nlon, ' ny=', nlat, ' ', trim(fname)

#ifdef DIROUT
      fname = trim(dirname) // '/' // trim(fname)
#endif
      call mygmt_grdio_d(tp,lon_west,lon_east,lat_south,lat_north, &
                         dx,dy,zmin,zmax,nlon,nlat,fname)

! === Conversion from flux to velocity should be done right after calc. ========
      else if(mode == VEL) then
! ==============================================================================
      ! Burbidge - added velocity output GMT files
      if(velgrd_flag == 1) then
         fx => wfld%fx
         fy => wfld%fy

         ! Burbidge: Now output vx
! === Conversion from flux to velocity should be done right after calc. ========
         do j = 1, nlat
            do i = 1, nlon
#ifndef MPI
               im = max(1,   i-1)
               ip = min(nlon,i+1)
#else
               im = i - 1
               ip = i + 1
               if(iand(bflag, WEST_BOUND) /= 0) im = max(1,im)
               if(iand(bflag, EAST_BOUND) /= 0) ip = min(nlon,ip)
#endif
               tdxm = 0.5d0*(dz(i,j) + hz(i,j) + dz(im,j) + hz(im,j))
               tdxp = 0.5d0*(dz(i,j) + hz(i,j) + dz(ip,j) + hz(ip,j))
               if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                  wod(ip,j) == 1 .and. tdxp > td_min .and. &
                  wod(i,j) == 1) then
                  tp(i,j) = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
               else
                  tp(i,j) = zap
               end if
            end do
         end do
! ==============================================================================

#ifndef MPI
! === To support over 1M steps. by tkato 2012/11/27 ============================
!        write(fname,'(a,a,i6.6,a)') trim(base), '.', istep, '-vx.grd'
         write(fname,'(a,a,i8.8,a)') trim(base), '.', istep, '-vx.grd'
! ==============================================================================
         call minmax_rwg(nlon,nlat,tp,zmax,zmin,imin,jmin,imax,jmax)
#else
! === To support over 1M steps. by tkato 2012/11/27 ============================
!        write(fname,'(a,a,i6.6,a,i6.6)') trim(base), '.', istep, '-vx.grd.', myrank
         write(fname,'(a,a,i8.8,a,i6.6)') trim(base), '.', istep, '-vx.grd.', myrank
! ==============================================================================
         call minmax_rwg(nlon,nlat,tp,zmax,zmin)
#endif
! === For ensemble =============================================================
#ifndef DIROUT
#ifdef MULTI
         fname = trim(members_dir) // trim(fname)
#endif
#endif
! ==============================================================================

! === To support over 1M steps. by tkato 2012/11/27 ============================
!        write(6,'(a,i6.6,a,f9.2,a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i4,a,i4,a,a)') &
#ifndef CARTESIAN
         write(6,'(a,i8.8,a,f9.2,a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i4,a,i4,a,a)') &
#else
         write(6,'(a,i8.8,a,f9.2,a,f9.3,a,f9.3,a,f0.1,a,f0.1,a,f0.1,a,f0.1,a,i4,a,i4,a,a)') &
#endif
! ==============================================================================
            'fx - it=', istep, ' t=', t, ' min=', zmin, ' max=', zmax, ' (', &
            lon_west, '/', lon_east, '/', lat_south, '/', lat_north, ') nx=', &
            nlon, ' ny=', nlat, ' ', trim(fname)

#ifdef DIROUT
         fname = trim(dirname) // '/' // trim(fname)
#endif
         call mygmt_grdio_d(tp,lon_west,lon_east,lat_south,lat_north, &
                            dx,dy,zmin,zmax,nlon,nlat,fname)

         ! Burbidge - Now do vy
! === Conversion from flux to velocity should be done right after calc. ========
         do j = 1, nlat
            do i = 1, nlon
#ifndef MPI
               jm = max(1,   j-1)
               jp = min(nlat,j+1)
#else
               jm = j - 1
               jp = j + 1
               if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,jm)
               if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(nlat,jp)
#endif
               tdym = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jm) + hz(i,jm))
               tdyp = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jp) + hz(i,jp))
               if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                  wod(i,jp) == 1 .and. tdyp > td_min .and. &
                  wod(i,j) == 1) then
                  tp(i,j) = 0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
               else
                  tp(i,j) = zap
               end if
            end do
         end do
! ==============================================================================

#ifndef MPI
! === To support over 1M steps. by tkato 2012/11/27 ============================
!        write(fname,'(a,a,i6.6,a)') trim(base), '.', istep, '-vy.grd'
         write(fname,'(a,a,i8.8,a)') trim(base), '.', istep, '-vy.grd'
! ==============================================================================

         call minmax_rwg(nlon,nlat,tp,zmax,zmin,imin,jmin,imax,jmax)
#else
! === To support over 1M steps. by tkato 2012/11/27 ============================
!        write(fname,'(a,a,i6.6,a,i6.6)') trim(base), '.', istep, '-vy.grd.', myrank
         write(fname,'(a,a,i8.8,a,i6.6)') trim(base), '.', istep, '-vy.grd.', myrank
! ==============================================================================

         call minmax_rwg(nlon,nlat,tp,zmax,zmin)
#endif
! === For ensemble =============================================================
#ifndef DIROUT
#ifdef MULTI
         fname = trim(members_dir) // trim(fname)
#endif
#endif
! ==============================================================================

! === To support over 1M steps. by tkato 2012/11/27 ============================
!        write(6,'(a,i6.6,a,f9.2,a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i4,a,i4,a,a)') &
#ifndef CARTESIAN
         write(6,'(a,i8.8,a,f9.2,a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i4,a,i4,a,a)') &
#else
         write(6,'(a,i8.8,a,f9.2,a,f9.3,a,f9.3,a,f0.1,a,f0.1,a,f0.1,a,f0.1,a,i4,a,i4,a,a)') &
#endif
! ==============================================================================
            'fy - it=', istep, ' t=', t, ' min=', zmin, ' max=', zmax, ' (', &
            lon_west, '/', lon_east, '/', lat_south, '/', lat_north, ') nx=', &
            nlon, ' ny=', nlat, ' ', trim(fname)

#ifdef DIROUT
         fname = trim(dirname) // '/' // trim(fname)
#endif
         call mygmt_grdio_d(tp,lon_west,lon_east,lat_south,lat_north, &
                            dx,dy,zmin,zmax,nlon,nlat,fname)
      end if
! === Speed output. ============================================================
         if(speedgrd_flag == 1) then
            fx => wfld%fx
            fy => wfld%fy

!$omp parallel do private(i, im, ip, jm, jp, tdxm, tdxp, tdym, tdyp, tx, ty, speed)
            do j = 1, nlat
               do i = 1, nlon
                  tx = 0.0d0
#ifndef MPI
                  im = max(1,   i-1)
                  ip = min(nlon,i+1)
#else
                  im = i - 1
                  ip = i + 1
                  if(iand(bflag, WEST_BOUND) /= 0) im = max(1,im)
                  if(iand(bflag, EAST_BOUND) /= 0) ip = min(nlon,ip)
#endif
                  tdxm = 0.5d0*(dz(i,j) + hz(i,j) + dz(im,j) + hz(im,j))
                  tdxp = 0.5d0*(dz(i,j) + hz(i,j) + dz(ip,j) + hz(ip,j))
                  if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                     wod(ip,j) == 1 .and. tdxp > td_min .and. &
                     wod(i,j) == 1) then
                     tx = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
                  end if

                  ty = 0.0d0
#ifndef MPI
                  jm = max(1,   j-1)
                  jp = min(nlat,j+1)
#else
                  jm = j - 1
                  jp = j + 1
                  if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,jm)
                  if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(nlat,jp)
#endif
                  tdym = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jm) + hz(i,jm))
                  tdyp = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jp) + hz(i,jp))
                  if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                     wod(i,jp) == 1 .and. tdyp > td_min .and. &
                     wod(i,j) == 1) then
                     ty = 0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
                  end if

                  speed = sqrt(tx**2 + ty**2)
                  tp(i,j) = speed
               end do
            end do

#ifndef MPI
            write(fname,'(a,a,i8.8,a)') trim(base), '.', istep, '-speed.grd'
            call minmax_rwg(nlon,nlat,tp,zmax,zmin,imin,jmin,imax,jmax)
#else
            write(fname,'(a,a,i8.8,a,i6.6)') trim(base), '.', istep, '-speed.grd.', myrank
            call minmax_rwg(nlon,nlat,tp,zmax,zmin)
#endif
! === For ensemble =============================================================
#ifndef DIROUT
#ifdef MULTI
            fname = trim(members_dir) // trim(fname)
#endif
#endif
! ==============================================================================

            write(6,'(a,i8.8,a,f9.2,a,f9.3,a,f9.3,a,f5.1,a,f5.1,a,f5.1,a,f5.1,a,i4,a,i4,a,a)') &
               'speed - it=', istep, ' t=', t, ' min=', zmin, ' max=', zmax, ' (', &
               lon_west, '/', lon_east, '/', lat_south, '/', lat_north, ') nx=', &
               nlon, ' ny=', nlat, ' ', trim(fname)

#ifdef DIROUT
            fname = trim(dirname) // '/' // trim(fname)
#endif
            call mygmt_grdio_d(tp,lon_west,lon_east,lat_south,lat_north, &
                               dx,dy,zmin,zmax,nlon,nlat,fname)

         end if
! ==============================================================================
! === Conversion from flux to velocity should be done right after calc. ========
      end if
! ==============================================================================

      return
   end subroutine dump_gmt_nl

   !thomas -  initialise wfld and rupture field, don't read in rupture field
   !          for use in multiple rupture version where reading in of rupture field
   !          must occur inside of the loop
   subroutine initl_wfld(wfld,zz,nlon,nlat)
      type(wave_arrays), target, intent(inout) :: wfld
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(out) :: zz
      integer(kind=4), intent(in) :: nlon, nlat

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz

      fx = 0.0d0
      fy = 0.0d0
      hz = 0.0d0
      zz = 0.0d0

      return
   end subroutine initl_wfld

   ! thomas -  version of initl_gene2 that doesn't initialise wave and rupture fields
   !           for use in multiple rupture version where reading in of rupture field
   !           must occur inside the loop, but initialising of wave field must not
   subroutine read_rupture(zz,nlon,nlat,mlat0,mlon0,dxdy,fname,pname)
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(inout) :: zz
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: mlat0, mlon0, dxdy
      character(len=256), intent(in) :: fname, pname

      integer(kind=4) :: i, j
      integer(kind=4) :: is, ie, js, je, ii, jj ! new
#ifndef MPI
      integer(kind=4) :: imax, jmax, imin, jmin
#endif
      real(kind=REAL_BYTE) :: zmax, zmin, west, east, south, north, slat0, slon0, dx, dy
      integer(kind=4) :: snx, sny
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: zz_tmp
      real(kind=REAL_BYTE) :: east_chk, west_chk, north_chk, south_chk

      ! RWG
      if(trim(fname) == 'NO_DISPLACEMENT_FILE_GIVEN') return

      call read_gmt_grd_hdr(fname, snx, sny, dx, dy, west, east, south, north, zmin, zmax)
      ! Burbidge : Use a tolerance here instead
#ifdef MPI
      !** assumes nlon=snx, nlat=sny **
      if(nlon /= snx .or. nlat /= sny) then
         write(0,'(a)') 'read_rupture(newsub.o): nlon or nlat does not same!'
         write(0,'(a,i0,a,i0,a,i0,a,i0,a)') '(nlon,nlat)=(', nlon, ',', nlat, ') (snx,sny)=(', snx, ',', sny, ')'
         call fatal_error(210)
      end if
#endif
      !** assumes dx=dy=dxdy **
      if(abs(dxdy - dx) > 0.0001d0 .or. abs(dxdy - dy) > 0.0001d0 .or. abs(dx - dy) > 0.0001d0) then
         write(0,'(a)') 'read_rupture(newsubs.o): dxdy not same!'
         write(0,'(a,f15.10,a,f15.10,a,f15.10)') 'dx=', dx, ' dy=', dy, ' dxdy=', dxdy
#ifndef MPI
         stop
#else
         call  fatal_error(209)
#endif
      end if

      write(6,'(8x,a,a)') 'read_rupture(newsubs.o): file name=', trim(fname)
      write(6,'(8x,a,i0,a,i0,a,e14.6,a,e14.6)') 'snx=', snx, ' sny=', sny, ' dx=', dx, ' dy=', dy
      write(6,'(8x,a,f0.3,a,f0.3,a,f0.3,a,f0.3,a,f0.6,a,f0.6)') &
         'west=', west, ' east=', east, ' south=', south, ' north=', north, &
         ' zmin=', zmin, ' zmax=', zmax

#ifndef CARTESIAN
      !*** colongitude ***
      if(west < 0.0d0) then
         slon0 = (360.0d0+west)*60.0d0
      else
         slon0 = west*60.0d0
      end if

      ! Burbidge - Fix up location if bathymetry straddles the date line
      ! First try adjusting the hemisphere. If that still does not work then exit with error
      !*** colatitude ***
      slat0 = (90.0d0-north)*60.0d0
#else
      slon0 = west
      slat0 = south
#endif

      write(6,'(8x,a,f0.0,a,f0.0)') 'slon0=', slon0, ' slat0=', slat0

      !** the (i,j) dimensions of the source grid wrt the depth/wavefield grids **
      !** is ie -> snx dlon **
      !** js je -> sny dlat **
      ! Burbidge - adjusted for normal node registration
#ifndef CARTESIAN
      is = int(abs(mlon0-slon0)/(dxdy*60.0d0)) + 1
      js = int(abs(mlat0-slat0)/(dxdy*60.0d0)) + 1
#else
      is = int(abs(mlon0-slon0)/dxdy) + 1
      js = int(abs(mlat0-slat0)/dxdy) + 1
#endif
! === DEBUG by tkato 2012/11/13 ================================================
      ie = is + snx - 1
      je = js + sny - 1
! ==============================================================================

#ifndef CARTESIAN
      write(6,'(a,f0.3)') 'mlon0/60.=', mlon0/60.0d0
      east_chk = (ie-1)*dxdy+mlon0/60.0d0
      if(east_chk > 180.0d0) east_chk = east_chk - 360.0d0

      west_chk = (is-1)*dxdy+mlon0/60.0d0
      if(west_chk > 180.0d0) west_chk = west_chk - 360.0d0

      south_chk = 90.0d0-((je-1)*dxdy+mlat0/60.0d0)
      north_chk = 90.0d0-((js-1)*dxdy+mlat0/60.0d0)

      write(6,'(a,f0.3,a,f0.3,a,f0.3,a,f0.3)') 'Source Check bounds: west=', west_chk, &
         ' east=', east_chk, ' south=', south_chk, ' north=', north_chk
#else
      write(6,'(a,f0.3)') 'mlon0=', mlon0
      east_chk = (ie-1)*dxdy+mlon0
      west_chk = (is-1)*dxdy+mlon0

      south_chk = (js-1)*dxdy+mlat0
      north_chk = (je-1)*dxdy+mlat0

      write(6,'(a,f0.3,a,f0.3,a,f0.3,a,f0.3)') 'Source Check bounds: west=', west_chk, &
         ' east=', east_chk, ' south=', south_chk, ' north=', north_chk
#endif
      ! Burbidge - Changed this to stdout so redirect works
      write(6,'(a,a,i0,a,i0,a,i0,a,i0)') &
         trim(pname), ': read_rupture(): is=', is, ' ie=', ie, ' js=', js, ' je=', je

      allocate(zz_tmp(snx,sny))
      call read_gmt_grd(fname,zz_tmp,snx,sny)

      !** map zz_tmp into zz **
      j = 1
! === DEBUG by tkato 2012/11/13 ================================================
!     do jj = js, je-1
      do jj = js, je
! ==============================================================================
         i = 1
! === DEBUG by tkato 2012/11/13 ================================================
!        do ii = is, ie-1
         do ii = is, ie
! ==============================================================================
            zz(ii,jj) = zz_tmp(i,j)
            i = i + 1
         end do
         j = j + 1
      end do
      deallocate(zz_tmp) ! free temporary space

#ifndef MPI
      call minmax_rwg(nlon, nlat, zz, zmax, zmin, imin, jmin, imax, jmax )
      write(6,'(a,a,a,a,f0.6,a,i0,a,i0,a,f0.6,a,i0,a,i0)') &
         trim(pname), ': read_rupture: ', trim(fname), &
         ' min=', zmin, ' imin=', imin, ' jmin=', jmin, &
         ' max=', zmax, ' imax=', imax, ' jmax=', jmax
#else
      call minmax_rwg(nlon, nlat, zz, zmax, zmin)
      write(6,'(a,a,a,a,f0.6,a,f0.6)') &
         trim(pname), ': read_rupture: ', trim(fname), ' min=', zmin, ' max=', zmax
#endif

      return
   end subroutine read_rupture

end module mod_newsubs
