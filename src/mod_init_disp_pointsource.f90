#include "real.h"
module mod_init_disp_pointsource
use mod_grid
use mod_rwg
implicit none

#ifndef CARTESIAN
real(kind=8), parameter, private :: R       = 6371.d0
real(kind=8), parameter, private :: PI      = 3.14159265d0
real(kind=8), parameter, private :: DEG2RAD = PI / 180.d0
#endif

real(kind=8), private :: h0 = 1.0d-3, lon_o = 0.0d0, lat_o = 0.0d0

contains

   subroutine specify_pointsource_params(file)
      character(len=256), intent(in) :: file
      namelist /pointsource/ h0, lon_o, lat_o

      open(1,file=trim(file),action='read',status='old',form='formatted')
      read(1,pointsource)
      close(1)

      write(6,'(a)') '================================================================'
      write(6,'(a)') '=== Check point source initialization! [Begin] ================='
      write(6,'(a)') '================================================================'
      write(6,'(a,a)') '- Filename: ', trim(file)
#ifndef CARTESIAN
      write(6,'(a,f15.6)') '- Wave height[m] (h0):          ', h0
      write(6,'(a,f15.6)') '- Center lon.[Degrees] (lon_o): ', lon_o
      write(6,'(a,f15.6)') '- Center lat.[Degrees] (lat_o): ', lat_o
#else
      write(6,'(a,f15.6)') '- Wave height[m] (h0):    ', h0
      write(6,'(a,f15.6)') '- Center lon.[m] (lon_o): ', lon_o
      write(6,'(a,f15.6)') '- Center lat.[m] (lat_o): ', lat_o
#endif
      write(6,'(a)') '================================================================'
      write(6,'(a)') '=== Check point source initialization! [End] ==================='
      write(6,'(a)') '================================================================'

      return
   end subroutine specify_pointsource_params

   subroutine make_pointsource_rupture(dgrid)
      type(data_grids), target, intent(inout) :: dgrid

      integer(kind=4) :: nx, ny
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: zz
#ifndef MPI
      real(kind=REAL_BYTE), pointer :: mlat0, mlon0
#else
      integer(kind=4), pointer :: totalNx, totalNy, kx, ky
      real(kind=REAL_BYTE), pointer :: glat0, glon0
#endif
      real(kind=REAL_BYTE), pointer :: dxdy
      real(kind=8) :: lat_north, lat_south, lon_east, lon_west, dx, dy
      integer(kind=4) :: i, j
#ifdef MPI
      integer(kind=4) :: ig, jg
#endif
      real(kind=8) :: lon, lat
#ifndef MPI
      integer(kind=4) :: imax, jmax, imin, jmin
#endif
      real(kind=REAL_BYTE) :: zmax, zmin

      nx = dgrid%my%nx
      ny = dgrid%my%ny
      zz => dgrid%zz
#ifdef MPI
      totalNx => dgrid%my%totalNx
      totalNy => dgrid%my%totalNy
      kx => dgrid%my%kx
      ky => dgrid%my%ky
#endif

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

#ifndef MPI
      i = (lon_o - lon_west)/dx
      j = (lat_north - lat_o)/dy
#else
      ig = (lon_o - lon_west)/dx
      jg = (lat_north - lat_o)/dy
      i = ig - kx + 1
      j = jg - ky + 1
#endif
      if((i >= 1) .and. (i <= nx) .and. (j >= 1) .and. (j <= ny)) zz(i,j) = h0

#ifndef MPI
      call minmax_rwg(nx, ny, zz, zmax, zmin, imin, jmin, imax, jmax )
      write(6,'(a,f0.6,a,i0,a,i0,a,f0.6,a,i0,a,i0)') &
         'make_pointsource_rupture: min=', zmin, ' imin=', imin, ' jmin=', jmin, &
         ' max=', zmax, ' imax=', imax, ' jmax=', jmax
#else
      call minmax_rwg(nx, ny, zz, zmax, zmin)
      write(6,'(a,f0.6,a,f0.6)') 'make_pointsource_rupture: min=', zmin, ' max=', zmax
#endif

      return
   end subroutine make_pointsource_rupture

end module mod_init_disp_pointsource
