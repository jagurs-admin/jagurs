#include "real.h"
module mod_tgs
use mod_grid
implicit none

! thomas - modification of tgs_open_rwg to include better checking of tide gauge positions
!          with respect to grid extents

contains

#ifndef __SX__
   subroutine tgs_open_ct(dg,ngrid,dt,nstep,nsta,program_name,tg_station_file_name,mytgs)
      type(data_grids), dimension(:), intent(inout) :: dg
      integer(kind=4), intent(in) :: ngrid
      real(kind=REAL_BYTE), intent(in) :: dt
      integer(kind=4), intent(in) :: nstep
      integer(kind=4), intent(inout) :: nsta
      character(len=256), intent(in) :: program_name
      character(len=128), intent(in) :: tg_station_file_name
      type(tgsrwg), allocatable, dimension(:), intent(out) :: mytgs

      integer(kind=4) :: ista
      integer(kind=4), parameter :: fp = 1
#ifndef CARTESIAN
      real(kind=REAL_BYTE) ::  colat, colon
#endif
      integer(kind=4) :: ig
      character(len=128) :: str
      real(kind=REAL_BYTE) :: mindh
      integer(kind=4) :: ilon, ilat
#ifdef MPI
      integer(kind=4) :: kg, ksta
#endif

      !*** read in the latitude and longitude ****
      open(fp,file=trim(tg_station_file_name),action='read',status='old',form='formatted',err=100)
      ista = 0

      read(fp,'(a)',err=101) str
      read(str,'(i20)',err=101) nsta
      if(nsta <= 0) goto 101

      allocate(mytgs(nsta))

#ifndef MPI
      do ista = 1, nsta
#else
      ista = 1
      do ksta = 1, nsta
#endif
         read(fp,'(a)') str
         read(str,*) mytgs(ista)%geolat, mytgs(ista)%geolon, mytgs(ista)%number

#ifndef CARTESIAN
         colat = 90.0d0 - mytgs(ista)%geolat
         if(mytgs(ista)%geolon < 0.0d0) then
            colon = 360.0d0 + mytgs(ista)%geolon
         else
            colon = mytgs(ista)%geolon
         end if

         mytgs(ista)%mcolat = colat*60.0d0
         mytgs(ista)%mcolon = colon*60.0d0
#else
         mytgs(ista)%mcolat = mytgs(ista)%geolat
         mytgs(ista)%mcolon = mytgs(ista)%geolon
#endif
         mytgs(ista)%ig = 0
         mytgs(ista)%ilon = 0
         mytgs(ista)%ilat = 0

         mindh = 2.0d0*dg(1)%my%dh
#ifdef MPI
         kg = -1
#endif
         do ig = 1, ngrid
            call tgs_find_grid_coords(mytgs(ista)%mcolon,mytgs(ista)%mcolat, &
#ifndef MPI
#ifndef CARTESIAN
                                      dg(ig)%my%mlon0,dg(ig)%my%mlat0,REAL_FUNC(dg(ig)%my%dh*60.0d0), &
#else
                                      dg(ig)%my%mlon0,dg(ig)%my%mlat0,REAL_FUNC(dg(ig)%my%dh), &
#endif
#else
#ifndef CARTESIAN
                                      dg(ig)%my%glon0,dg(ig)%my%glat0,REAL_FUNC(dg(ig)%my%dh*60.0d0), &
#else
                                      dg(ig)%my%glon0,dg(ig)%my%glat0,REAL_FUNC(dg(ig)%my%dh), &
#endif
#endif
                                      ilon,ilat)
! === DEBUG by tkato 2015/03/04 ================================================
#ifdef CARTESIAN
#ifndef MPI
            ilat = dg(ig)%my%ny - ilat + 1
#else
            ilat = dg(ig)%my%totalNy - ilat + 1
#endif
#endif
! ==============================================================================
#ifndef MPI
            if(ilon < 1 .or. ilon > dg(ig)%my%nx .or. ilat < 1 .or. ilat > dg(ig)%my%ny) then
               !  tide gauge is not in extent of dg[ig]
               if(ig == 1) then
                  ! we assume that ig=1 is the ancestor of all grids and if the tide gauge is
                  ! outside of that then exit
                  write(0,'(a,i0,a)') 'Tide Gauge Point ', ista, ' appears to be outside the box.'
                  ! DB added some more information here
                  write(0,'(a,f18.6,a,f18.6)') 'Station geolon=', mytgs(ista)%geolon, ' geolat=', mytgs(ista)%geolat
                  write(0,'(a,i0,a,i0,a,i0,a,i0)') 'ilon=', ilon, ' ilat=', ilat, ' nx=', dg(ig)%my%nx, ' ny=', dg(ig)%my%ny
                  stop
               else
                  ! in this case ig is not the big grid and tide gauge is allowed to be
                  !  outside, so  continue on to the next iteration of the loop
                  cycle
               end if
            end if

            if(dg(ig)%my%dh < mindh) then
               mindh = dg(ig)%my%dh
               mytgs(ista)%ig = ig

               mytgs(ista)%ilon = ilon
               mytgs(ista)%ilat = ilat
               mytgs(ista)%dt = dt
               mytgs(ista)%nt = nstep

               mytgs(ista)%z = dg(ig)%depth_field%dz(mytgs(ista)%ilon,mytgs(ista)%ilat)
            end if
         end do
#else
            if(ilon >= 1 .and. ilon <= dg(ig)%my%totalNx .and. ilat >= 1 .and. ilat <= dg(ig)%my%totalNy) then
               if(dg(ig)%my%dh < mindh) then
                  mindh = dg(ig)%my%dh
                  kg = ig
               end if
            end if
         end do
         if(kg < 0) then
            write(0,'(a,i0,a)') 'Tide Gauge Point ', ksta, ' appears to be outside the box.'
            ! DB added some more information here
            write(0,'(a,f18.6,a,f18.6)') 'Station geolon=', mytgs(ista)%geolon, ' geolat=', mytgs(ista)%geolat
            write(0,'(a,i0,a,i0,a,i0,a,i0)') 'ilon=', ilon, ' ilat=', ilat, &
               ' nx=', dg(1)%my%totalNx, ' ny=', dg(1)%my%totalNy
            call fatal_error(303)
         end if

         ! set if tgs exists in local area
         call tgs_find_grid_coords(mytgs(ista)%mcolon,mytgs(ista)%mcolat, &
#ifndef CARTESIAN
                                   dg(kg)%my%glon0,dg(kg)%my%glat0,REAL_FUNC(dg(kg)%my%dh*60.0d0), &
#else
                                   dg(kg)%my%glon0,dg(kg)%my%glat0,REAL_FUNC(dg(kg)%my%dh), &
#endif
                                   ilon,ilat)
! === DEBUG by tkato 2015/03/04 ================================================
#ifdef CARTESIAN
         ilat = dg(kg)%my%totalNy - ilat + 1
#endif
! ==============================================================================
         if(ilon >= dg(kg)%my%ix .and. ilon <= dg(kg)%my%ixend .and. &
            ilat >= dg(kg)%my%iy .and. ilat <= dg(kg)%my%iyend) then
            mytgs(ista)%ig = kg

            mytgs(ista)%ilon = ilon - dg(kg)%my%kx + 1
            mytgs(ista)%ilat = ilat - dg(kg)%my%ky + 1
            mytgs(ista)%dt = dt
            mytgs(ista)%nt = nstep

            mytgs(ista)%z = dg(kg)%depth_field%dz(mytgs(ista)%ilon,mytgs(ista)%ilat)

            ista = ista + 1
         end if
#endif
      end do
      close(fp)
#ifdef MPI
      nsta = ista - 1
#endif

      write(6,'(a,a,i0,a,a)') trim(program_name), ': tgs_open_rwg.c: nsta=', nsta, &
         ' read from file ', trim(tg_station_file_name)

      return

100   write(0,'(a,a)') 'missing tide gauge station ', trim(tg_station_file_name)
      nsta = 0
      return

101   write(0,'(a,a)') 'missing or malformed tide gauge station file ', trim(tg_station_file_name)
      nsta = 0
      return
   end subroutine tgs_open_ct
#endif

   ! thomas - function to compute row and column position of point with respect to a grid
   !          point can be outside of grid
   subroutine tgs_find_grid_coords(pt_m_colon,  & ! colongitude of point in minutes
                                   pt_m_colat,  & ! colatitude of point in minutes
                                   grd_m_colon, & ! colongitude of NW corner of grid in minutes
                                   grd_m_colat, & ! colat of NW corner of grid in minutes
                                   dh_m,        & ! grid resolution in minutes
                                   ilon_ptr,    & ! grid position of point east-west
                                   ilat_ptr)      ! grid position of point north-south, could be -ve
      real(kind=REAL_BYTE), intent(inout) :: pt_m_colon, pt_m_colat, grd_m_colon, grd_m_colat
      real(kind=REAL_BYTE), intent(in) :: dh_m
      integer(kind=4), intent(out) :: ilon_ptr, ilat_ptr

#ifndef CARTESIAN
      ! just in case, let's standardise the colongitudes
      do while(pt_m_colon < 0.0d0)
         pt_m_colon = pt_m_colon + 360.0d0*60.0d0
      end do
      do while(pt_m_colon > 60.0d0*360.0d0)
         pt_m_colon = pt_m_colon - 60.0d0*360.0d0
      end do
      do while(grd_m_colon < 0.0d0)
         grd_m_colon = grd_m_colon + 360.0d0*60.0d0
      end do
      do while(grd_m_colon > 60.0d0*360.0d0)
         grd_m_colon = grd_m_colon - 60.0d0*360.0d0
      end do

      if(pt_m_colon < grd_m_colon - 0.5d0*dh_m) then
         ! perhaps the grid is big enough to wrap around the point
         ! so always make *ilon_ptr >= 0
         pt_m_colon = pt_m_colon + 360.0d0*60.0d0
      end if
#endif

      ilon_ptr = int((pt_m_colon - grd_m_colon)/dh_m + 0.5d0) + 1
      ilat_ptr = int((pt_m_colat - grd_m_colat)/dh_m + 0.5d0) + 1

      return
   end subroutine tgs_find_grid_coords

end module mod_tgs
