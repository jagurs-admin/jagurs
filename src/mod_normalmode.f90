#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_normalmode
#ifdef MPI
use mod_mpi_fatal
#endif
implicit none
real(kind=REAL_BYTE) :: srclon, srclat, tmax, dt
real(kind=REAL_BYTE), allocatable, dimension(:,:) :: nm_atm
integer(kind=4) :: dist_s, dist_e, sind, minmode, maxmode
real(kind=REAL_BYTE) :: Samp, SDt, OH
real(kind=REAL_BYTE), parameter :: rhow = 1.0d3 ! water density (kg/m3)
integer(kind=4) :: nstep, interval
logical :: DCflag, SYNflag
character(len=256) :: atmfile

contains

   subroutine normalmode_read_namelist(dt_model)
      real(kind=REAL_BYTE), intent(in) :: dt_model
      integer(kind=4) :: i, j
      namelist /normalmode/ srclon, srclat, dist_s, dist_e, tmax, dt, sind, minmode, maxmode, &
         Samp, SDt, DCflag, SYNflag, OH, atmfile

      open(1,file='normalmode.namelist',action='read',status='old',form='formatted')
      read(1,normalmode)
      close(1)
      write(6,*) '[Normal mode] srclon, srclat, dist_s, dist_e, tmax, dt: ', srclon, srclat, dist_s, dist_e, tmax, dt

      nstep = nint(tmax/dt) + 1
      allocate(nm_atm(dist_s:dist_e,nstep))
      open(1,file='normalmode.dat',action='read',status='old',form='unformatted')
      read(1) nm_atm
      close(1)

      if(dt_model > dt) then
         write(0,*) '[Normal mode] ERROR! "dt" on JAGURS is grater than normal-mode interval!'
         write(0,*) '[Normal mode] JAGRUS dt[s]:          ', dt_model
         write(0,*) '[Normal mode] Normalmode file dt[s]: ', dt
#ifndef MPI
         stop
#else
         call fatal_error(401)
#endif
      end if

      interval = (dt + dt_model/10.0d0)/dt_model
      write(6,*) '[Normal mode] JAGRUS dt[s]:          ', dt_model
      write(6,*) '[Normal mode] Normalmode file dt[s]: ', dt
      write(6,*) '[Normal mode] Interval [steps]:      ', interval

      return
   end subroutine normalmode_read_namelist

   subroutine deg2dist(lat1, lon1, lat2, lon2, dist)
      real(kind=REAL_BYTE), intent(in) :: lat1, lon1, lat2, lon2
      real(kind=REAL_BYTE), intent(out) :: dist
      real(kind=REAL_BYTE), parameter :: R = 6371.0d3
      real(kind=REAL_BYTE), parameter :: pi = 3.14159265d0
      real(kind=REAL_BYTE) :: phi1, phi2, lamda1, lamda2

      phi1 = lat1*pi/180.0d0
      phi2 = lat2*pi/180.0d0
      lamda1 = lon1*pi/180.0d0
      lamda2 = lon2*pi/180.0d0

      dist = R*acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lamda1 - lamda2))

      return
   end subroutine deg2dist

   subroutine make_nm_ind(nlon, nlat, mlat0, mlon0, dxdy, nm_ind)
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: mlat0, mlon0, dxdy
      integer(kind=4), intent(out) :: nm_ind(nlon, nlat)
      real(kind=REAL_BYTE) :: north, west, lat1, lon1, dist
      integer(kind=4) :: i, j, minind, maxind

      north = 90.0d0 - mlat0/60.0d0
      west  = mlon0/60.0d0
      if(west > 180.0d0) west = west - 360.0d0

      minind = dist_e
      maxind = dist_s
!$omp parallel do private(i,j,lat1,lon1,dist) reduction(min:minind) reduction(max:maxind)
      do j = 1, nlat
         lat1 = north - (j - 1)*dxdy
         do i = 1, nlon
            lon1 = west + (i - 1)*dxdy
            call deg2dist(lat1, lon1, srclat, srclon, dist)
            nm_ind(i,j) = nint(dist/1000.0d0)
            if(minind > nm_ind(i,j)) minind = nm_ind(i,j)
            if(maxind < nm_ind(i,j)) maxind = nm_ind(i,j)
         end do
      end do

      if((minind < dist_s) .or. (maxind > dist_e)) then
         write(0,*) '[Normal mode] ERROR! The range accessed is out of normal-mode range!'
         write(0,*) '[Normal mode] The range accessed: ', minind, '-', maxind
         write(0,*) '[Normal mode] Normalmode dist_s-dist_e: ', dist_s, '-', dist_e
#ifndef MPI
         stop
#else
         call fatal_error(402)
#endif
      end if

#if 0
      write(0,*) 'nlon, nlat: ', nlon, nlat
      write(0,*) 'north, west: ', north, west
      write(0,*) 'dxdy: ', dxdy
      do j = 1, nlat
         lat1 = north - (j - 1)*dxdy
         do i = 1, nlon
            lon1 = west + (i - 1)*dxdy
            write(0,'(3i8,2f8.3)') j, i, nm_ind(i,j), lon1, lat1
         end do
      end do
#endif

      return
   end subroutine make_nm_ind

   subroutine calc_nm_P(nlon, nlat, ind, P, P0, P1, istep)
      integer(kind=4), intent(in) :: nlon, nlat, istep
      integer(kind=4), intent(in) :: ind(nlon, nlat)
      real(kind=REAL_BYTE), intent(out) :: P(nlon, nlat)
      real(kind=REAL_BYTE), intent(inout) :: P0(nlon, nlat), P1(nlon, nlat)
      real(kind=REAL_BYTE) :: Ptmp(nlon, nlat)
      integer(kind=4) :: i, j, l, m
      real(kind=REAL_BYTE) :: c0, c1

      l = (istep - 1)/interval + 1
      m = mod(istep - 1, interval)

      if(((l+1>nstep).and.(m>0)).or.(l>nstep)) then
         P = 0.0d0
         return
      end if

      if(m == 0) then
         if(istep == 1) then
!$omp parallel do private(i,j)
            do j = 1, nlat
               do i = 1, nlon
                  P(i,j) = nm_atm(ind(i,j),l)
               end do
            end do
            call interp_P(nlon, nlat, P, P0)
         else
            P0 = P1
         end if

         if(l + 1 <= nstep) then
!$omp parallel do private(i,j)
            do j = 1, nlat
               do i = 1, nlon
                  P(i,j) = nm_atm(ind(i,j),l+1)
               end do
            end do
            call interp_P(nlon, nlat, P, P1)
         else
            write(6,'(a,i0,a)') '[Normal mode] Final read at step ', istep, '!'
            P1 = 0.0d0
         end if
      end if

      c0 = dble(interval - m)/dble(interval)
      c1 = dble(m)/dble(interval)
!$omp parallel do private(i,j)
      do j = 1, nlat
         do i = 1, nlon
            P(i,j) = c0*P0(i,j) + c1*P1(i,j)
         end do
      end do

      return
   end subroutine calc_nm_P

   subroutine interp_P(nlon, nlat, Pin, Pout)
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: Pin(nlon, nlat)
      real(kind=REAL_BYTE), intent(out) :: Pout(nlon, nlat)
      integer(kind=4) :: i, j

!$omp parallel do private(i,j)
      do j = 1,nlat
         do i = 1,nlon
            if((i == 1) .and. (j == 1)) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                  &
                         & 1.0d0/ 8.0d0*(Pin(i+1,j) + Pin(i,j+1)) + &
                         & 1.0d0/16.0d0*Pin(i+1,j+1)
            else if((i == nlon) .and. (j == nlat)) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                  &
                         & 1.0d0/ 8.0d0*(Pin(i-1,j) + Pin(i,j-1)) + &
                         & 1.0d0/16.0d0*Pin(i-1,j-1)
            else if((i == 1) .and. (j == nlat)) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                  &
                         & 1.0d0/ 8.0d0*(Pin(i+1,j) + Pin(i,j-1)) + &
                         & 1.0d0/16.0d0*Pin(i+1,j-1)
            else if((i == nlon) .and. (j == 1)) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                  &
                         & 1.0d0/ 8.0d0*(Pin(i-1,j) + Pin(i,j+1)) + &
                         & 1.0d0/16.0d0*Pin(i-1,j+1)
            else if(i == 1) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                               &
                         & 1.0d0/ 8.0d0*(Pin(i+1,j) + Pin(i,j+1) + Pin(i,j-1)) + &
                         & 1.0d0/16.0d0*(Pin(i+1,j+1) + Pin(i+1,j-1))
            else if(i == nlon) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                               &
                         & 1.0d0/ 8.0d0*(Pin(i-1,j) + Pin(i,j+1) + Pin(i,j-1)) + &
                         & 1.0d0/16.0d0*(Pin(i-1,j+1) + Pin(i-1,j-1))
            else if(j == 1) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                               &
                         & 1.0d0/ 8.0d0*(Pin(i+1,j) + Pin(i-1,j) + Pin(i,j+1)) + &
                         & 1.0d0/16.0d0*(Pin(i+1,j+1) + Pin(i-1,j+1))
            else if(j == nlat) then
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                               &
                         & 1.0d0/ 8.0d0*(Pin(i+1,j) + Pin(i-1,j) + Pin(i,j-1)) + &
                         & 1.0d0/16.0d0*(Pin(i+1,j-1) + Pin(i-1,j-1) )
            else
               Pout(i,j) = 1.0d0/ 4.0d0*Pin(i,j) +                                             &
                         & 1.0d0/ 8.0d0*(Pin(i+1,j) + Pin(i-1,j) + Pin(i,j+1) + Pin(i,j-1) ) + &
                         & 1.0d0/16.0d0*(Pin(i+1,j+1) + Pin(i+1,j-1) + Pin(i-1,j+1) + Pin(i-1,j-1))
            end if
         end do
      end do

      return
   end subroutine interp_P
end module mod_normalmode
