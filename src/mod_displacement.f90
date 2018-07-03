!#define QUAD
#include "real.h"
module mod_displacement
#ifdef _OPENMP
use omp_lib
#endif
use mod_grid
use mod_params, only : hzdisp_effect, fault_param_file, apply_kj_filter
use mod_mygmt_gridio, only : read_gmt_grd_hdr
#ifdef MPI
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
#endif
#ifdef MPI
use mod_mpi, only : exchange_edges_zz
#endif
implicit none

! Fault parameters.
integer(kind=4), private :: num_faults
real(kind=8), private, allocatable, dimension(:) :: xref, yref
real(kind=8), private, allocatable, dimension(:) :: length, width, depth, dip
real(kind=8), private, allocatable, dimension(:) :: us, ud, ut, strike
#ifndef CARTESIAN
real(kind=8), private, allocatable, dimension(:) :: lon_0, lat_0
real(kind=8), private, allocatable, dimension(:) :: h0, lat1
#else
real(kind=8), private, allocatable, dimension(:) :: h0
#endif

! Parameters.
#ifndef CARTESIAN
real(kind=8), private, parameter :: rote    = 6371.0d3
#endif
real(kind=8), private, parameter :: M_PI    = 3.14159265d0
real(kind=8), private, parameter :: DEG2RAD = 3.14159265d0/180.d0
real(kind=8), private, parameter :: ALP     = 0.5d0

#ifdef MPI
integer(kind=4), private :: nprocs, myrank, npx, npy, rankx, ranky
integer(kind=4), private :: MPI_X_WORLD, MPI_Y_WORLD
#endif

contains

#ifdef  MPI
   subroutine displacement_mpi_initialize(nprocs_in, myrank_in, npx_in, npy_in, rankx_in, ranky_in)
      use mpi
      integer(kind=4), intent(in)  :: nprocs_in, myrank_in, npx_in, npy_in, rankx_in, ranky_in
      integer(kind=4) :: color, key, ierr
      nprocs = nprocs_in
      myrank = myrank_in
      npx = npx_in
      npy = npy_in
      rankx = rankx_in
      ranky = ranky_in

      color = ranky
      key = rankx
#ifndef MULTI
      call MPI_comm_split(MPI_COMM_WORLD, color, key, MPI_X_WORLD, ierr)
#else
      call MPI_comm_split(MPI_MEMBER_WORLD, color, key, MPI_X_WORLD, ierr)
#endif

      color = rankx
      key = ranky
#ifndef MULTI
      call MPI_comm_split(MPI_COMM_WORLD, color, key, MPI_Y_WORLD, ierr)
#else
      call MPI_comm_split(MPI_MEMBER_WORLD, color, key, MPI_Y_WORLD, ierr)
#endif

      return
   end subroutine displacement_mpi_initialize
#endif

#ifndef MULTI
   subroutine displacement_initialize()
#else
   subroutine displacement_initialize(input_dirname)
      character(len=128), intent(in) :: input_dirname
      character(len=256) :: file_multi
#endif
#ifndef CARTESIAN
      real(kind=8), allocatable, dimension(:) :: lat, lon, rake, slip_amp
#else
      real(kind=8), allocatable, dimension(:) :: rake, slip_amp
#endif
      character(len=512) :: buf
      real(kind=8), dimension(9) :: tmp
      integer(kind=4) :: n, num_lines

      write(6,'(a,a,a)') '[displacement] Open fault parameter file ', trim(fault_param_file), '!'
#ifndef MULTI
      open(1,file=trim(fault_param_file),action='read',status='old',form='formatted')
#else
      file_multi = trim(input_dirname) // trim(fault_param_file)
      open(1,file=trim(file_multi),action='read',status='old',form='formatted')
#endif

      num_lines = 0
      num_faults = 0
      do while(.true.)
         num_lines = num_lines + 1
         read(1,'(a)',end=100) buf
         buf = adjustl(buf)
         if(buf(1:1) == '!') cycle
         read(buf,*,end=101,err=101) tmp
         num_faults = num_faults + 1
      end do
100   write(6,'(a,i6)') '[displacement] Number of faults: ', num_faults

      ! local
#ifndef CARTESIAN
      allocate(lat     (num_faults))
      allocate(lon     (num_faults))
#endif
      allocate(rake    (num_faults))
      allocate(slip_amp(num_faults))

      ! global
      allocate(depth   (num_faults))
      allocate(length  (num_faults))
      allocate(width   (num_faults))
      allocate(dip     (num_faults))
      allocate(strike  (num_faults))
      allocate(xref    (num_faults))
      allocate(yref    (num_faults))
      allocate(us      (num_faults))
      allocate(ud      (num_faults))
      allocate(ut      (num_faults))
#ifndef CARTESIAN
      allocate(lon_0   (num_faults))
      allocate(lat_0   (num_faults))
#endif
      if(apply_kj_filter == 1) then
         allocate(h0  (num_faults))
#ifndef CARTESIAN
         allocate(lat1(num_faults))
#endif
      end if

      rewind(1)
      n = 0
      do while(.true.)
         read(1,'(a)',end=200) buf
         buf = adjustl(buf)
         if(buf(1:1) == '!') cycle
         n = n + 1
#ifndef CARTESIAN
         read(buf,*) lat(n), lon(n), depth(n), length(n), width(n), &
                     dip(n), strike(n), rake(n), slip_amp(n)
#else
         read(buf,*) yref(n), xref(n), depth(n), length(n), width(n), &
                     dip(n), strike(n), rake(n), slip_amp(n)
#endif
      end do
200   continue
      close(1)

      write(6,'(a)') '+--------------------------------------------&
                     &----------------------------------------------+'
#ifndef CARTESIAN
      write(6,'(a)') '|lat       lon       depth     length    width&
                     &     dip       strike    rake      slip_amp  |'
#else
      write(6,'(a)') '|yref      xref      depth     length    width&
                     &     dip       strike    rake      slip_amp  |'
#endif
      write(6,'(a)') '+--------------------------------------------&
                     &----------------------------------------------+'
      do n = 1, num_faults
#ifndef CARTESIAN
         write(6,'(a,9f10.3,a)') '|', lat(n), lon(n), depth(n), length(n), width(n), &
#else
         write(6,'(a,2f10.1,7f10.3,a)') '|', yref(n), xref(n), depth(n), length(n), width(n), &
#endif
            dip(n), strike(n), rake(n), slip_amp(n), '|'
      end do
      write(6,'(a)') '+--------------------------------------------&
                     &----------------------------------------------+'

      do n = 1, num_faults
#ifndef CARTESIAN
         lon_0(n) = lon(n)
         lat_0(n) = lat(n)
         call mapproject_initialize(lon_0(n))
         call ll2xy(1, 1, lon(n), lat(n), xref(n), yref(n))
         call mapproject_finalize()
#endif
         us(n) = slip_amp(n)*dcos(rake(n)*DEG2RAD)*1.0d2
         ud(n) = slip_amp(n)*dsin(rake(n)*DEG2RAD)*1.0d2
         ut(n) = 0.d0*1.0d2
      end do

      ! local
#ifndef CARTESIAN
      deallocate(lat)
      deallocate(lon)
#endif
      deallocate(rake)
      deallocate(slip_amp)

      return

101   write(0,'(a,i0,a)') 'ERROR: Invalid format on line ', num_lines, ' in ''fault''!'
      stop
   end subroutine displacement_initialize

   subroutine displacement_finalize()
      ! global
      deallocate(depth)
      deallocate(length)
      deallocate(width)
      deallocate(dip)
      deallocate(strike)
      deallocate(xref)
      deallocate(yref)
      deallocate(us)
      deallocate(ud)
      deallocate(ut)
#ifndef CARTESIAN
      deallocate(lon_0)
      deallocate(lat_0)
#endif
      if(apply_kj_filter == 1) then
         deallocate(h0)
#ifndef CARTESIAN
         deallocate(lat1)
#endif
      end if
      return 
   end subroutine displacement_finalize

   subroutine displacement_calc_displacement(dg, ig)
! === Specify lower limit of depth to adopt horizontal displacement effect. ====
      use mod_params, only : min_depth_hde
! ==============================================================================
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), intent(in) :: ig

      integer(kind=4), pointer :: nlon, nlat
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dz, zz

      real(kind=8), allocatable, dimension(:,:) :: uz, uh1, uh2
#ifndef CARTESIAN
      real(kind=8), allocatable, dimension(:,:,:) :: lonlat, xy, dxy
      real(kind=8), allocatable, dimension(:,:,:) :: lltmp
      integer(kind=4), allocatable, dimension(:,:) :: mask
#else
      real(kind=8), allocatable, dimension(:,:,:) :: xy, dxy
#endif

      real(kind=REAL_BYTE) :: x_inc, y_inc, x_min, x_max, y_min, y_max, z_min, z_max
      real(kind=REAL_BYTE) :: umin, umax, ulim
#ifndef CARTESIAN
      real(kind=8) :: lat, lon, cd, sd, ct, st, dep
#else
      real(kind=8) :: cd, sd, ct, st, dep
#endif
#ifndef __SX__
      real(kind=8) :: x, y, u1, u2, u3, dumm
#else
      real(kind=8) :: x, y, u1, u2, u3
#endif
      real(kind=8) :: dist
      integer(kind=4) :: i, j, n
#ifdef MPI
      integer(kind=4) :: ierr
#endif
! === Check distance============================================================
      real(kind=8) :: p, q, xi, et, r, d, rd
      integer(kind=4) :: num_distance_zero, l_, w_
! ==============================================================================
#if defined(MPI) && defined(ONEFILE)
      integer(kind=4) :: dummynx, dummyny
#endif

      write(6,'(a)') '[displacement] Initial displacement with fault calculation.'

      nlon => dg%my%nx
      nlat => dg%my%ny
      dz   => dg%depth_field%dz
      zz   => dg%zz

      allocate(uz    (nlon,nlat))
      allocate(uh1   (nlon,nlat))
      allocate(uh2   (nlon,nlat))
#ifndef CARTESIAN
      allocate(lonlat(nlon,nlat,2))
#endif
      allocate(xy    (nlon,nlat,2))
      allocate(dxy   (nlon,nlat,2)); dxy = 0.0d0
#ifndef CARTESIAN
      allocate(lltmp (nlon,nlat,2))
      allocate(mask  (nlon,nlat))
#endif

#if !defined(MPI) || !defined(ONEFILE)
      call read_gmt_grd_hdr(dg%my%bath_file,nlon,nlat,x_inc,y_inc,x_min,x_max,y_min,y_max,z_min,z_max)
#ifdef MPI
#ifndef MULTI
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, x_min, 1, REAL_MPI, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, y_min, 1, REAL_MPI, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, x_max, 1, REAL_MPI, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, y_max, 1, REAL_MPI, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, x_min, 1, REAL_MPI, MPI_MIN, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, y_min, 1, REAL_MPI, MPI_MIN, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, x_max, 1, REAL_MPI, MPI_MAX, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, y_max, 1, REAL_MPI, MPI_MAX, MPI_MEMBER_WORLD, ierr)
#endif
#endif
#else
      if(myrank == 0) then
         call read_gmt_grd_hdr(dg%my%bath_file,dummynx,dummyny,x_inc,y_inc,x_min,x_max,y_min,y_max,z_min,z_max)
      end if
#ifndef MULTI
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(x_min, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(y_min, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(x_max, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(y_max, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
#else
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(x_min, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(y_min, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(x_max, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(y_max, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
#endif
#endif
      x_inc = real(x_inc)
      y_inc = real(y_inc)
      x_min = real(x_min)
      x_max = real(x_max)
      y_min = real(y_min)
      y_max = real(y_max)
      write(6,'(a,4f15.6)') '[displacement] x_min, x_inc, y_min, y_inc: ', x_min, x_inc, y_min, y_inc

!$omp parallel
#ifndef CARTESIAN
!$omp do private(lat, i, lon)
      do j = 1, nlat
#ifndef MPI
         lat = y_min + y_inc*(j - 1)
#else
         lat = y_min + y_inc*(dg%my%totalNy - dg%my%kyend + j - 1)
#endif
         do i = 1, nlon
#ifndef MPI
            lon = x_min + x_inc*(i - 1)
#else
            lon = x_min + x_inc*(dg%my%kx + i - 2)
#endif
            lonlat(i,nlat-j+1,1) = real(lon)
            lonlat(i,nlat-j+1,2) = real(lat)
         end do
      end do

!$omp single
#ifndef MPI
      write(6,'(a,2f15.6)') '[displacement] lon-range: ', x_min, x_min + x_inc*(nlon - 1)
      write(6,'(a,2f15.6)') '[displacement] lat-range: ', y_min, y_min + y_inc*(nlat - 1)
#else
      write(6,'(a,2f15.6)') '[displacement] lon-range: ', x_min + x_inc*(dg%my%kx - 1), x_min + x_inc*(dg%my%kx + nlon - 2)
      write(6,'(a,2f15.6)') '[displacement] lat-range: ', y_min + y_inc*(dg%my%totalNy - dg%my%kyend), &
                                                          y_min + y_inc*(dg%my%totalNy - dg%my%kyend + nlat - 1)
#endif
!$omp end single
#else
!$omp do private(y, i, x)
      do j = 1, nlat
#ifndef MPI
         y = y_min + y_inc*(j - 1)
#else
         y = y_min + y_inc*(dg%my%totalNy - dg%my%kyend + j - 1)
#endif
         do i = 1, nlon
#ifndef MPI
            x = x_min + x_inc*(i - 1)
#else
            x = x_min + x_inc*(dg%my%kx + i - 2)
#endif
            xy(i,nlat-j+1,1) = real(x)
            xy(i,nlat-j+1,2) = real(y)
         end do
      end do

!$omp single
#ifndef MPI
      write(6,'(a,2f15.6)') '[displacement] x-range: ', x_min, x_min + x_inc*(nlon - 1)
      write(6,'(a,2f15.6)') '[displacement] y-range: ', y_min, y_min + y_inc*(nlat - 1)
#else
      write(6,'(a,2f15.6)') '[displacement] x-range: ', x_min + x_inc*(dg%my%kx - 1), x_min + x_inc*(dg%my%kx + nlon - 2)
      write(6,'(a,2f15.6)') '[displacement] y-range: ', y_min + y_inc*(dg%my%totalNy - dg%my%kyend), &
                                                        y_min + y_inc*(dg%my%totalNy - dg%my%kyend + nlat - 1)
#endif
!$omp end single
#endif

!$omp do private(i)
      do j = 1, nlat
         do i = 1, nlon
            zz(i,j) = 0.0d0
         end do
      end do
!$omp end parallel

! === Check distance============================================================
      do n = 1, num_faults
#ifndef CARTESIAN
!$omp parallel do private(i, r)
         do j = 1, nlat
            do i = 1, nlon
               r = dsqrt((lonlat(i,j,1) - lon_0(n))**2 + (lonlat(i,j,2) - lat_0(n))**2)
               if(r > 30.0d0) then
                  mask(i,j) = 0
                  lltmp(i,j,1) = lon_0(n)
                  lltmp(i,j,2) = lat_0(n)
               else
                  mask(i,j) = 1
                  lltmp(i,j,1) = lonlat(i,j,1)
                  lltmp(i,j,2) = lonlat(i,j,2)
               end if
            end do
         end do

         ! TODO: Investigation of thread safety of PROJ.4 if it's high-cost.
         call mapproject_initialize(lon_0(n))
         call ll2xy(nlon, nlat, lltmp(1,1,1), lltmp(1,1,2), xy(1,1,1), xy(1,1,2))
         call mapproject_finalize()
#endif

         cd = dcos(dip(n)*DEG2RAD)
         sd = dsin(dip(n)*DEG2RAD)

         if(dabs(cd) < 1.d-3) then
            cd = 0.d0
            if(sd .gt. 0.d0) then
               sd = 1.D0
            end if
            if(sd .le. 0.d0) then
               sd = -1.D0
            end if
         end if
         dep = depth(n) + width(n)*sd
         ct  = dcos(strike(n)*DEG2RAD)
         st  = dsin(strike(n)*DEG2RAD)

         num_distance_zero = 0
!$omp parallel do private(i, x, y, p, q, xi, et, r, d, rd, l_, w_)
         do j = 1, nlat
            do i = 1, nlon
#ifndef CARTESIAN
               if(mask(i,j) == 1) then
#endif
                  x =  (xy(i,j,1)-xref(n))/1000.0d0*st + (xy(i,j,2)-yref(n))/1000.0d0*ct
                  y = -(xy(i,j,1)-xref(n))/1000.0d0*ct + (xy(i,j,2)-yref(n))/1000.0d0*st + width(n)*cd
! ------------------------------------------------------------------------------
                  p = y*cd + dep*sd
                  q = y*sd - dep*cd
                  do l_ = 0, 1
                     xi = x - l_*length(n)
                     do w_ = 0, 1
                        et = p - w_*width(n)
                        r = dsqrt(xi*xi+et*et+q*q)
                        d = et*sd - q*cd
                        rd = r + d
                        if(rd < tiny(rd)) then
!$omp critical
                           num_distance_zero = num_distance_zero + 1
!$omp end critical
                        end if
                     end do
                  end do
! ------------------------------------------------------------------------------
#ifndef CARTESIAN
               end if
#endif
            end do
         end do
#ifdef MPI
#ifndef MULTI
         call MPI_Allreduce(MPI_IN_PLACE, num_distance_zero, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
         call MPI_Allreduce(MPI_IN_PLACE, num_distance_zero, 1, MPI_INTEGER, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
#endif
         if(num_distance_zero > 0) then
            write(6,'(a)')      '[displacement] ========================= WARNING!!! BEGIN ========================='
            write(6,'(a,i5,a)') '[displacement] There exist the grid points in which distance to the source number ', n, ' is zero.'
            write(6,'(a)')      '[displacement] So, it is displaced +1.0d-6 [m] to x direction.'
            write(6,'(a)')      '[displacement] ========================= WARNING!!! END ==========================='
            xref(n) = xref(n) + 1.0d-6
         end if
      end do
! ==============================================================================
      do n = 1, num_faults
#ifndef CARTESIAN
!$omp parallel do private(i, r)
         do j = 1, nlat
            do i = 1, nlon
               r = dsqrt((lonlat(i,j,1) - lon_0(n))**2 + (lonlat(i,j,2) - lat_0(n))**2)
               if(r > 30.0d0) then
                  mask(i,j) = 0
                  lltmp(i,j,1) = lon_0(n)
                  lltmp(i,j,2) = lat_0(n)
               else
                  mask(i,j) = 1
                  lltmp(i,j,1) = lonlat(i,j,1)
                  lltmp(i,j,2) = lonlat(i,j,2)
               end if
            end do
         end do

         ! TODO: Investigation of thread safety of PROJ.4 if it's high-cost.
         call mapproject_initialize(lon_0(n))
         call ll2xy(nlon, nlat, lltmp(1,1,1), lltmp(1,1,2), xy(1,1,1), xy(1,1,2))
         call mapproject_finalize()
#endif
!$omp parallel
!$omp do private(i)
         do j = 1, nlat
            do i = 1, nlon
               uz(i,j) = 0.d0
            end do
         end do

!$omp single
         cd = dcos(dip(n)*DEG2RAD)
         sd = dsin(dip(n)*DEG2RAD)

         if(dabs(cd) < 1.d-3) then
            cd = 0.d0
            if(sd .gt. 0.d0) then
               sd = 1.D0
            end if
            if(sd .le. 0.d0) then
               sd = -1.D0
            end if
         end if
         dep = depth(n) + width(n)*sd
         ct  = dcos(strike(n)*DEG2RAD)
         st  = dsin(strike(n)*DEG2RAD)
!$omp end single

#ifndef __SX__
!$omp do private(i, x, y, u1, u2, u3, dumm)
#else
!$omp do private(i, x, y, u1, u2, u3)
#endif
         do j = 1, nlat
            do i = 1, nlon
#ifndef CARTESIAN
               if(mask(i,j) == 1) then
#endif
                  x =  (xy(i,j,1)-xref(n))/1000.0d0*st + (xy(i,j,2)-yref(n))/1000.0d0*ct
                  y = -(xy(i,j,1)-xref(n))/1000.0d0*ct + (xy(i,j,2)-yref(n))/1000.0d0*st + width(n)*cd
#ifndef __SX__
                  call srectf(ALP, x, y, dep, length(n), width(n), sd, cd, us(n), ud(n), ut(n), &
                     u1, u2, u3, dumm, dumm, dumm, dumm, dumm, dumm)
#else
                  call srectf(ALP, x, y, dep, length(n), width(n), sd, cd, us(n), ud(n), ut(n), &
                     u1, u2, u3)
#endif
                  uh1(i,j) =  u1*sin(strike(n)*DEG2RAD) - u2*cos(strike(n)*DEG2RAD)
                  uh2(i,j) = -u1*cos(strike(n)*DEG2RAD) - u2*sin(strike(n)*DEG2RAD)
                  uz(i,j)  =  u3
#ifndef CARTESIAN
               else
                  uh1(i,j) = 0.0d0
                  uh2(i,j) = 0.0d0
                  uz(i,j)  = 0.0d0
               end if
#endif
            end do
         end do

!$omp single
         ulim = 0.d0
!$omp end single
!$omp do private(i)
         do j = 1, nlat
            do i = 1, nlon
               if(dabs(uz(i,j)) < ulim) then
                  uz(i,j) = 0.d0
               end if
            end do
         end do

         ! cm to m
!$omp do private(i)
         do j = 1, nlat
            do i = 1, nlon
               uh1(i,j) =  uh1(i,j)*0.01d0
               uh2(i,j) = -uh2(i,j)*0.01d0
               uz(i,j)  =  uz(i,j) *0.01d0
            end do
         end do

         if(hzdisp_effect == 1) then
!$omp single
            write(6,'(a)') '[displacement] Horizontal displacement effect is adopted!'
!$omp end single
!$omp do private(i, dist)
            do j = 1, nlat
               do i = 2, nlon-1
#ifndef CARTESIAN
                  if((mask(i+1,j) == 1) .and. (mask(i,j) == 1) .and. (mask(i-1,j) == 1)) then
#endif
                     dist = dsqrt((xy(i+1,j,1)-xy(i-1,j,1))**2 &
                                + (xy(i+1,j,2)-xy(i-1,j,2))**2)
                     dxy(i,j,1) = (dz(i+1,j) - dz(i-1,j))/dist
#ifndef CARTESIAN
                  else
                     dxy(i,j,1) = 0.0d0
                  end if
#endif
               end do
            end do

!$omp do
            do j = 1, nlat
               dxy(1,   j,1) = dxy(2,     j,1)
               dxy(nlon,j,1) = dxy(nlon-1,j,1)
            end do

!$omp do private(i, dist)
            do j = 2, nlat-1
               do i = 1, nlon
#ifndef CARTESIAN
                  if((mask(i,j+1) == 1) .and. (mask(i,j) == 1) .and. (mask(i,j-1) == 1)) then
#endif
                     dist = dsqrt((xy(i,j+1,1)-xy(i,j-1,1))**2 &
                                + (xy(i,j+1,2)-xy(i,j-1,2))**2)
                     dxy(i,j,2) = (dz(i,j+1) - dz(i,j-1))/dist
#ifndef CARTESIAN
                  else
                     dxy(i,j,2) = 0.0d0
                  end if
#endif
               end do
            end do

!$omp do
            do i = 1, nlon
               dxy(i,1,   2) = dxy(i,2,     2)
               dxy(i,nlat,2) = dxy(i,nlat-1,2)
            end do

!$omp do private(i)
            do j = 1, nlat
               do i = 1, nlon
! === Specify lower limit of depth to adopt horizontal displacement effect. ====
                  if(dz(i,j) > min_depth_hde) then
! ==============================================================================
                  uz(i,j) = uz(i,j) + uh1(i,j)*dxy(i,j,1) + uh2(i,j)*dxy(i,j,2)
! === Specify lower limit of depth to adopt horizontal displacement effect. ====
                  end if
! ==============================================================================
               end do
            end do
         end if
!$omp end parallel

         if(apply_kj_filter == 1) then
#ifndef CARTESIAN
            if(ig == 1) call displacement_calc_h0_lat1(dg, uz, nlon, nlat, h0(n), lat1(n))
            call displacement_apply_kj_filter(dg, uz, nlon, nlat, h0(n), lat1(n))
#else
            if(ig == 1) call displacement_calc_h0(dg, uz, nlon, nlat, h0(n))
            call displacement_apply_kj_filter(dg, uz, nlon, nlat, h0(n))
#endif
         end if

!$omp parallel
!$omp do private(i)
         do j = 1, nlat
            do i = 1, nlon
               zz(i,j) = zz(i,j) + uz(i,j)
            end do
         end do
!$omp end parallel
      end do

#ifdef MPI
      call exchange_edges_zz(dg)
#endif

!$omp parallel
!$omp single
      umax = -huge(umax)
!$omp end single
!$omp do private(i) reduction(max:umax)
      do j = 1, nlat
         do i = 1, nlon
            umax = max(umax,zz(i,j))
         end do
      end do

!$omp single
      umin = huge(umin)
!$omp end single
!$omp do private(i) reduction(min:umin)
      do j = 1, nlat
         do i = 1, nlon
            umin = min(umin,zz(i,j))
         end do
      end do
!$omp end parallel
#ifdef MPI
#ifndef MULTI
      call MPI_Allreduce(MPI_IN_PLACE, umax, 1, REAL_MPI, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, umin, 1, REAL_MPI, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(MPI_IN_PLACE, umax, 1, REAL_MPI, MPI_MAX, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, umin, 1, REAL_MPI, MPI_MIN, MPI_MEMBER_WORLD, ierr)
#endif
#endif
      write(6,'(a,2f15.6)') '[displacement] umax, umin: ', umax, umin

      deallocate(uz)
      deallocate(uh1)
      deallocate(uh2)
#ifndef CARTESIAN
      deallocate(lonlat)
#endif
      deallocate(xy)
      deallocate(dxy)
#ifndef CARTESIAN
      deallocate(lltmp)
      deallocate(mask)
#endif

      return
   end subroutine displacement_calc_displacement

#ifndef CARTESIAN
   subroutine displacement_apply_kj_filter(dg, zz, nlon, nlat, h0, lat1)
#else
   subroutine displacement_apply_kj_filter(dg, zz, nlon, nlat, h0)
#endif
#ifndef __SX__
      include 'fftw3.f'
#endif
      type(data_grids), target, intent(inout) :: dg 
      real(kind=8), dimension(nlon,nlat), intent(inout) :: zz
      integer(kind=4), intent(in) :: nlon, nlat
#ifndef CARTESIAN
      real(kind=8), intent(in) :: h0, lat1
#else
      real(kind=8), intent(in) :: h0
#endif
      
! ==============================================================================
! === Initialize FFT Begin =====================================================
! ==============================================================================
      ! Local vals.
      integer(kind=4) :: m_nxg, m_nyg, N_X, N_Y, N
#ifdef MPI
      integer(kind=4) :: ix, iy
#endif
#ifndef __SX__
      integer(kind=8), allocatable, dimension(:) :: &
         xplan_forward, yplan_forward, xplan_backward, yplan_backward ! FFTW3 plan
#else
      integer(kind=4), allocatable, dimension(:,:) :: ifax_x, ifax_y
      real(kind=8), allocatable, dimension(:,:) :: trigs_x, trigs_y
      real(kind=8), allocatable, dimension(:,:) :: work_x
      complex(kind=8), allocatable, dimension(:,:) :: work_y
#endif
      complex(kind=8), pointer, dimension(:,:) :: xfftbuf, yfftbuf
      real(kind=8), pointer, dimension(:,:) :: realbuf

      integer(kind=4) :: i, j
#ifdef MPI
      integer(kind=4) :: nx0, ny0, nx1, ny1, nx2, ny2, ibias, jbias
      real(kind=8), allocatable, dimension(:) :: sendbuf1, recvbuf1
      complex(kind=8), allocatable, dimension(:) :: sendbuf2, recvbuf2
      integer(kind=4), allocatable, dimension(:) :: sendcounts1, sdispls1, recvcounts1, rdispls1
      integer(kind=4), allocatable, dimension(:) :: sendcounts2, sdispls2, recvcounts2, rdispls2

      integer(kind=4) :: ixend, iyend, kx, kxend, ky, kyend
      integer(kind=4) :: p, j_, ind, ierr
#endif
! === OpenMP ===================================================================
      integer(kind=4), allocatable, dimension(:) :: its2, ite2, itn2 ! N_X, nx2
      integer(kind=4), allocatable, dimension(:) :: jts1, jte1, jtn1 ! N_Y, ny1
      integer(kind=4) :: it, nthreads
! ==============================================================================
#ifdef __SX__
#ifndef MPI
      integer(kind=4) :: ierr
#endif
#endif
! ==============================================================================
! === Initialize FFT End =======================================================
! ==============================================================================
      real(kind=REAL_BYTE) :: x_inc, y_inc, x_min, x_max, y_min, y_max, z_min, z_max
      real(kind=8) :: m_dx, m_dy, coef_norm
      real(kind=8) :: dk_x, dk_y, k_x, k_y, k
      integer(kind=4) :: nlon_dummy, nlat_dummy
#ifdef MPI
      integer(kind=4) :: num, xst, xlen, yst, ylen, ist, i_
#endif
#ifndef CARTESIAN
      real(kind=8) :: x, lon1, lat2, lon2
#endif

      write(6,'(a)') '[displacement] Kajiura filter is applied!'

! ==============================================================================
! === Initialize FFT Begin =====================================================
! ==============================================================================
#ifndef MPI
      m_nxg = dg%my%nx
      m_nyg = dg%my%ny
#else
      m_nxg = dg%my%totalNx
      m_nyg = dg%my%totalNy
#endif

      N_X = m_nxg
      N_Y = m_nyg
      if(m_nxg > m_nyg) then
         N_X = m_nxg
         N_Y = N_X
      else
         N_Y = m_nyg
         N_X = N_Y
      end if
      N_X = 2**ceiling(dlog(dble(N_X))/dlog(2.0d0))
      N_Y = 2**ceiling(dlog(dble(N_Y))/dlog(2.0d0))
      N = N_X*N_Y

      write(6,'(a,2i6)') '[displacement] m_nxg, m_nyg: ', m_nxg, m_nyg
      write(6,'(a,2i6,i12)') '[displacement] N_X, N_Y, N: ', N_X, N_Y, N
#ifdef MPI
      ix = dg%my%ix
      iy = dg%my%iy
      ixend = dg%my%ixend
      iyend = dg%my%iyend

      kx = dg%my%kx
      ky = dg%my%ky
      kxend = dg%my%kxend
      kyend = dg%my%kyend

      ! Space grid without edges.
      nx0 = ixend - ix + 1 ! xsize.
      ny0 = iyend - iy + 1 ! ysize.
      ibias = ix - kx ! x bias to skip edges.
      jbias = iy - ky ! y bias to skip edges.

      write(6,'(a,2i6)') '[displacement] nprocs, myrank: ', nprocs, myrank
      write(6,'(a,2i6)') '[displacement] ix, ixend: ', ix, ixend
      write(6,'(a,2i6)') '[displacement] iy, iyend: ', iy, iyend
      write(6,'(a,2i6)') '[displacement] nx0, ny0: ', nx0, ny0

      ! Space grid for x-FFT.
      nx1 = dg%my%totalNx
      ny1 = ny0/npx
      if(rankx < mod(ny0,npx)) ny1 = ny1 + 1

      write(6,'(a,2i6)') '[displacement] nx1, ny1: ', nx1, ny1

      ! MPI_Alltoallv1: (nx0,ny0) -> (nx1,ny1)
      allocate(sendbuf1(nx0*ny0))
      allocate(sendcounts1(0:npx-1))
      allocate(sdispls1(0:npx-1))

      allocate(recvbuf1(nx1*ny1))
      allocate(recvcounts1(0:npx-1))
      allocate(rdispls1(0:npx-1))

      do p = 0, npx-1
         j = ny0/npx
         if(p < mod(ny0,npx)) j = j + 1
         sendcounts1(p) = nx0*j
      end do
      call MPI_Alltoall(sendcounts1, 1, MPI_INTEGER, &
                        recvcounts1, 1, MPI_INTEGER, MPI_X_WORLD, ierr)
      sdispls1(0) = 0
      rdispls1(0) = 0
      do p = 1, npx-1
         sdispls1(p) = sdispls1(p-1) + sendcounts1(p-1)
         rdispls1(p) = rdispls1(p-1) + recvcounts1(p-1)
      end do

      write(6,'(a)') '[displacement] p, sendcounts1(p), sdispls1(p), recvcounts1(p), rdispls1(p)'
      do p = 0, npx-1
         write(6,'(5i8)') p, sendcounts1(p), sdispls1(p), recvcounts1(p), rdispls1(p)
      end do

      ! Grid for y-FFT.
      nx2 = (N_X/2+1)/nprocs
      if(myrank < mod(N_X/2+1,nprocs)) nx2 = nx2 + 1
      ny2 = N_Y

      write(6,'(a,2i6)') '[displacement] nx2, ny2: ', nx2, ny2

      ! MPI_Alltoallv1: (N_X,ny1) -> (N_Y,nx2)
      allocate(sendbuf2((N_X/2+1)*ny1))
      allocate(sendcounts2(0:nprocs-1))
      allocate(sdispls2(0:nprocs-1))

      allocate(recvbuf2(N_Y*nx2))
      allocate(recvcounts2(0:nprocs-1))
      allocate(rdispls2(0:nprocs-1))

      do p = 0, nprocs-1
         i = (N_X/2+1)/nprocs
         if(p < mod(N_X/2+1,nprocs)) i = i + 1
         sendcounts2(p) = i*ny1
      end do
#ifndef MULTI
      call MPI_Alltoall(sendcounts2, 1, MPI_INTEGER, &
                        recvcounts2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoall(sendcounts2, 1, MPI_INTEGER, &
                        recvcounts2, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      sdispls2(0) = 0
      rdispls2(0) = 0
      do p = 1, nprocs-1
         sdispls2(p) = sdispls2(p-1) + sendcounts2(p-1)
         rdispls2(p) = rdispls2(p-1) + recvcounts2(p-1)
      end do

      write(6,'(a)') '[displacement] p, sendcounts2(p), sdispls2(p), recvcounts2(p), rdispls2(p)'
      do p = 0, nprocs-1
         write(6,'(5i8)') p, sendcounts2(p), sdispls2(p), recvcounts2(p), rdispls2(p)
      end do
#endif

! === OpenMP ===================================================================
#ifdef _OPENMP
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif

#ifndef __SX__
      allocate(xplan_forward(0:nthreads-1))
      allocate(yplan_forward(0:nthreads-1))
      allocate(xplan_backward(0:nthreads-1))
      allocate(yplan_backward(0:nthreads-1))
#else
      allocate(ifax_x(20,0:nthreads-1))
      allocate(ifax_y(20,0:nthreads-1))
      allocate(trigs_x(N_X*2,0:nthreads-1))
      allocate(trigs_y(N_Y*2,0:nthreads-1))
#ifndef MPI
      allocate(work_x(N_X+2,N_Y))
      allocate(work_y(N_Y,N_X/2+1))
#else
      allocate(work_x(N_X+2,ny1))
      allocate(work_y(N_Y,nx2))
#endif
#endif
      allocate(its2(0:nthreads-1))
      allocate(ite2(0:nthreads-1))
      allocate(itn2(0:nthreads-1))
      allocate(jts1(0:nthreads-1))
      allocate(jte1(0:nthreads-1))
      allocate(jtn1(0:nthreads-1))

#ifndef MPI
      do it = 0, nthreads-1
         itn2(it) = (N_X/2+1)/nthreads
         if(it < mod(N_X/2+1,nthreads)) itn2(it) = itn2(it) + 1
         jtn1(it) = N_Y/nthreads
         if(it < mod(N_Y,nthreads)) jtn1(it) = jtn1(it) + 1
      end do
#else
      do it = 0, nthreads-1
         itn2(it) = nx2/nthreads
         if(it < mod(nx2,nthreads)) itn2(it) = itn2(it) + 1
         jtn1(it) = ny1/nthreads
         if(it < mod(ny1,nthreads)) jtn1(it) = jtn1(it) + 1
      end do
#endif
      its2(0) = 1
      jts1(0) = 1
      do it = 1, nthreads-1
         its2(it) = its2(it-1) + itn2(it-1)
         jts1(it) = jts1(it-1) + jtn1(it-1)
      end do
      do it = 0, nthreads-1
         ite2(it) = its2(it) + itn2(it) - 1
         jte1(it) = jts1(it) + jtn1(it) - 1
      end do
! ==============================================================================

#ifndef __SX__
#ifndef MPI
      allocate(realbuf(N_X,N_Y))
      allocate(xfftbuf(N_X/2+1,N_Y))
      allocate(yfftbuf(N_Y,N_X/2+1))
#else
      allocate(realbuf(N_X,ny1))
      allocate(xfftbuf(N_X/2+1,ny1))
      allocate(yfftbuf(N_Y,nx2))
#endif
#else
#ifndef MPI
      allocate(realbuf(N_X+2,N_Y))
      allocate(xfftbuf(N_X/2+1,N_Y))
      allocate(yfftbuf(N_Y+1,N_X/2+1))
#else
      allocate(realbuf(N_X+16,ny1))
      allocate(xfftbuf(N_X/2+1,ny1))
      allocate(yfftbuf(N_Y+16,nx2))
#endif
#endif

      ! First touch!
!$omp parallel do private(i, j)
      do it = 0, nthreads-1
         do j = jts1(it), jte1(it)
            do i = 1, N_X
               realbuf(i,j) = 0.0d0
            end do
         end do
         do j = jts1(it), jte1(it)
            do i = 1, N_X/2+1
               xfftbuf(i,j) = dcmplx(0.0d0,0.0d0)
            end do
         end do
         do i = its2(it), ite2(it)
            do j = 1, N_Y
               yfftbuf(j,i) = dcmplx(0.0d0,0.0d0)
            end do
         end do
      end do

      do it = 0, nthreads-1
#ifndef __SX__
         call dfftw_plan_many_dft_r2c(xplan_forward(it), 1, N_X, jtn1(it), &
                                      realbuf(1,jts1(it)), 0, 1, N_X,      &
                                      xfftbuf(1,jts1(it)), 0, 1, N_X/2+1,  &
                                      FFTW_MEASURE)

         call dfftw_plan_many_dft(yplan_forward(it), 1, N_Y, itn2(it), &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,      &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,      &
                                  FFTW_FORWARD, FFTW_MEASURE)

         call dfftw_plan_many_dft_c2r(xplan_backward(it), 1, N_X, jtn1(it), &
                                      xfftbuf(1,jts1(it)), 0, 1, N_X/2+1,   &
                                      realbuf(1,jts1(it)), 0, 1, N_X,       &
                                      FFTW_MEASURE)

         call dfftw_plan_many_dft(yplan_backward(it), 1, N_Y, itn2(it), &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,       &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,       &
                                  FFTW_BACKWARD, FFTW_MEASURE)
#else
         call DFRMFB(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+16,0,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
         call ZFCMFB(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+16,0,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
#endif
      end do
! ==============================================================================
! === Initialize FFT End =======================================================
! ==============================================================================

! ==============================================================================
! === Calc. m_dx and m_dy Begin ================================================
! ==============================================================================
#ifdef MPI
      if(myrank == 0) then
#endif
      call read_gmt_grd_hdr(dg%my%bath_file,nlon_dummy,nlat_dummy,x_inc,y_inc,x_min,x_max,y_min,y_max,z_min,z_max)
#ifdef MPI
      end if
#endif
#ifdef MPI
#ifndef MULTI
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
#else
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
#endif
#endif
#ifndef CARTESIAN
      m_dy  = y_inc*DEG2RAD*rote ! Average Y space step

      lon1 = x_min
#ifdef MPI
#ifndef MULTI
      call MPI_Allreduce(MPI_IN_PLACE, lon1, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(MPI_IN_PLACE, lon1, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_MEMBER_WORLD, ierr)
#endif
#endif
      lat2 = lat1
      lon2 = lon1 + x_inc
      x = dacos(dsin(lat1*DEG2RAD)*dsin(lat2*DEG2RAD) &
              + dcos(lat1*DEG2RAD)*dcos(lat2*DEG2RAD)*dcos((lon2-lon1)*DEG2RAD))
      m_dx = rote*x
      write(6,'(a)')        '[displacement] Parameters to calc. sphere distance'
      write(6,'(a,4f15.6)') '[displacement] lat1, lon1, lat2, lon2: ', lat1, lon1, lat2, lon2
#else
      m_dx  = x_inc
      m_dy  = y_inc
#endif
      write(6,'(a,2f15.6)') '[displacement] m_dx, m_dy: ', m_dx, m_dy
! ==============================================================================
! === Calc. m_dx and m_dy End ==================================================
! ==============================================================================

! ==============================================================================
! === Apply filter Begin =======================================================
! ==============================================================================
      coef_norm = 1.0d0/(dble(N_X)*dble(N_Y))

!$omp parallel
#ifndef MPI
!$omp do private(j, i)
      do it = 0, nthreads-1
         do j = jts1(it), jte1(it)
            do i = 1, N_X
               realbuf(i,j) = 0.0d0
            end do
         end do
      end do

!$omp do private(i)
      do j = 1, nlat
         do i = 1, nlon
            realbuf(i,j) = zz(i,j)
         end do
      end do

      ! Forward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute_dft_r2c(xplan_forward(it),realbuf(1,jts1(it)),xfftbuf(1,jts1(it)))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+2,1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif

      ! Transposiion: (N_X,N_Y) -> (N_Y,N_X)
!$omp do private(i, j)
      do it = 0, nthreads-1
         do i = its2(it), ite2(it)
            do j = 1, N_Y
#ifndef __SX__
               yfftbuf(j,i) = xfftbuf(i,j)
#else
               yfftbuf(j,i) = dcmplx(realbuf(2*i-1,j),realbuf(2*i,j))
#endif
            end do
         end do
      end do

      ! Forward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_forward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+1,1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

!$omp single
      dk_x = 2.0d0*M_PI/(m_dx*dble(N_X))
      dk_y = 2.0d0*M_PI/(m_dy*dble(N_Y))
!$omp end single
!$omp do private(i, k_x, j, k_y, k)
      do it = 0, nthreads-1
         do i = its2(it), ite2(it)
            k_x = (i - 1)*dk_x
            if(i > N_X/2 + 1) k_x = (N_X - i + 1)*dk_x
            do j = 1, N_Y
               k_y = (j - 1)*dk_y
               if(j > N_Y/2 + 1) k_y = (N_Y - j + 1)*dk_y
               k = dsqrt(k_x**2 + k_y**2)
               yfftbuf(j,i) = yfftbuf(j,i)/dcosh(k*h0)*coef_norm
            end do
         end do
      end do

      ! Backward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_backward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+1,-1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

      ! Transposiion: (N_Y,N_X) -> (N_X,N_Y)
!$omp do private(j, i)
      do it = 0, nthreads-1
         do j = jts1(it), jte1(it)
            do i = 1, N_X/2+1
#ifndef __SX__
               xfftbuf(i,j) = yfftbuf(j,i)
#else
               realbuf(2*i-1,j) = dble (yfftbuf(j,i))
               realbuf(2*i,  j) = dimag(yfftbuf(j,i))
#endif
            end do
         end do
      end do

      ! Backward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute_dft_c2r(xplan_backward(it),xfftbuf(1,jts1(it)),realbuf(1,jts1(it)))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+2,-1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif

!$omp do private(i)
      do j = 1, nlat
         do i = 1, nlon
            zz(i,j) = realbuf(i,j)
         end do
      end do
#else
      ! Transposiion: (nx0,ny0) -> (N_X,ny1)
!$omp do private(i, ind)
      do j = 1, ny0
         do i = 1, nx0
            ind = (j-1)*nx0 + i
            sendbuf1(ind) = 0.0d0
         end do
      end do

!$omp do private(j_, i, i_, ind)
      do j = 1, ny0
         j_ = jbias + j
         do i = 1, nx0
            i_ = ibias + i
            ind = (j-1)*nx0 + i
            sendbuf1(ind) = zz(i_,j_)
         end do
      end do

!$omp single
      call MPI_Alltoallv(sendbuf1, sendcounts1, sdispls1, MPI_DOUBLE_PRECISION, &
                         recvbuf1, recvcounts1, rdispls1, MPI_DOUBLE_PRECISION, MPI_X_WORLD, ierr)
!$omp end single

!$omp do private(i)
      do j = 1, ny1
#ifndef __SX__
         do i = 1, N_X
#else
         do i = 1, N_X+16
#endif
            realbuf(i,j) = 0.0d0
         end do
      end do

!$omp do private(xst, xlen, j, i, ind, i_)
      do p = 0, npx-1
         xst = rdispls1(p)/ny1
         xlen = recvcounts1(p)/ny1
         do j = 1, ny1
            do i = 1, xlen
               ind = rdispls1(p) + (j-1)*xlen + i
               i_ =  xst + i
               realbuf(i_,j) = recvbuf1(ind)
            end do
         end do
      end do

      ! Forward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute_dft_r2c(xplan_forward(it),realbuf(1,jts1(it)),xfftbuf(1,jts1(it)))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+16,1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif

      ! Transposiion: (N_X,ny1) -> (N_Y,nx2)
#ifndef __SX__
!$omp do private(j, ind)
      do i = 1, N_X/2+1
         do j = 1, ny1
            ind = ny1*(i-1) + j
            sendbuf2(ind) = xfftbuf(i,j)
         end do
      end do
#else
!$omp do private(i, ind)
      do j = 1, ny1
         do i = 1, N_X/2+1
            ind = ny1*(i-1) + j
            sendbuf2(ind) = dcmplx(realbuf(2*i-1,j),realbuf(2*i,j))
         end do
      end do
#endif

!$omp single
#ifndef MULTI
      call MPI_Alltoallv(sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, &
                         recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoallv(sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, &
                         recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, MPI_MEMBER_WORLD, ierr)
#endif
!$omp end single

!$omp do private(j)
      do i = 1, nx2
#ifndef __SX__
         do j = 1, N_Y
#else
         do j = 1, N_Y+16
#endif
            yfftbuf(j,i) = dcmplx(0.0d0,0.0d0)
         end do
      end do

!$omp do private(yst, ylen, i, j, ind, j_)
      do p = 0, nprocs-1
         yst = rdispls2(p)/nx2
         ylen = recvcounts2(p)/nx2
         do i = 1, nx2
            do j = 1, ylen
               ind = rdispls2(p) + (i-1)*ylen + j
               j_ =  yst + j
               yfftbuf(j_,i) = recvbuf2(ind)
            end do
         end do
      end do

      ! Forward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_forward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+16,1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

      ! Calc. F(hz)*F(g)
!$omp single
      num = (N_X/2+1)/nprocs
      ist = myrank*num
      ist = ist + min(mod(N_X/2+1,nprocs),myrank) + 1
      dk_x = 2.0d0*M_PI/(m_dx*dble(N_X))
      dk_y = 2.0d0*M_PI/(m_dy*dble(N_Y))
!$omp end single
!$omp do private(i, i_, k_x, j, k_y, k)
      do it = 0, nthreads-1
         do i = its2(it), ite2(it)
            i_ = ist + i - 1
            k_x = (i_ - 1)*dk_x
            if(i_ > N_X/2 + 1) k_x = (N_X - i_ + 1)*dk_x
            do j = 1, N_Y
               k_y = (j - 1)*dk_y
               if(j > N_Y/2 + 1) k_y = (N_Y - j + 1)*dk_y
               k = dsqrt(k_x**2 + k_y**2)
               yfftbuf(j,i) = yfftbuf(j,i)/dcosh(k*h0)*coef_norm
            end do
         end do
      end do

      ! Backward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_backward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+16,-1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

      ! Transposiion: (N_Y,nx2) -> (N_X,ny1)
!$omp do private(yst, ylen, i, j, ind, j_)
      do p = 0, nprocs-1
         yst = rdispls2(p)/nx2
         ylen = recvcounts2(p)/nx2
         do i = 1, nx2
            do j = 1, ylen
               ind = rdispls2(p) + (i-1)*ylen + j
               j_ =  yst + j
               recvbuf2(ind) = yfftbuf(j_,i)
            end do
         end do
      end do

!$omp single
#ifndef MULTI
      call MPI_Alltoallv(recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, &
                         sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoallv(recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, &
                         sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, MPI_MEMBER_WORLD, ierr)
#endif
!$omp end single

#ifndef __SX__
!$omp do private(i, ind)
      do j = 1, ny1
         do i = 1, N_X/2+1
            ind = ny1*(i-1) + j
            xfftbuf(i,j) = sendbuf2(ind)
         end do
      end do
#else
!$omp do private(i, ind)
      do j = 1, ny1
         do i = 1, N_X/2+1
            ind = ny1*(i-1) + j
            realbuf(2*i-1,j) = dble (sendbuf2(ind))
            realbuf(2*i,  j) = dimag(sendbuf2(ind))
         end do
      end do
#endif

      ! Backward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute_dft_c2r(xplan_backward(it),xfftbuf(1,jts1(it)),realbuf(1,jts1(it)))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+16,-1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif

      ! Transposiion: (N_X,ny1) -> (nx0,ny0)
!$omp do private(xst, xlen, j, i, ind, i_)
      do p = 0, npx-1
         xst = rdispls1(p)/ny1
         xlen = recvcounts1(p)/ny1
         do j = 1, ny1
            do i = 1, xlen
               ind = rdispls1(p) + (j-1)*xlen + i
               i_ =  xst + i
               recvbuf1(ind) = realbuf(i_,j)
            end do
         end do
      end do

!$omp single
      call MPI_Alltoallv(recvbuf1, recvcounts1, rdispls1, MPI_DOUBLE_PRECISION, &
                         sendbuf1, sendcounts1, sdispls1, MPI_DOUBLE_PRECISION, MPI_X_WORLD, ierr)
!$omp end single

!$omp do private(j_, i, i_, ind)
      do j = 1, ny0
         j_ = jbias + j
         do i = 1, nx0
            i_ = ibias + i
            ind = (j-1)*nx0 + i
            zz(i_,j_) = sendbuf1(ind)
         end do
      end do
#endif
!$omp end parallel
! ==============================================================================
! === Apply filter End =========================================================
! ==============================================================================

! ==============================================================================
! === Finalize FFT Begin =======================================================
! ==============================================================================
#ifndef __SX__
      do it = 0, nthreads-1
         call dfftw_destroy_plan(xplan_forward(it))
         call dfftw_destroy_plan(yplan_forward(it))
         call dfftw_destroy_plan(xplan_backward(it))
         call dfftw_destroy_plan(yplan_backward(it))
      end do

      deallocate(xplan_forward)
      deallocate(yplan_forward)
      deallocate(xplan_backward)
      deallocate(yplan_backward)
#else
      deallocate(ifax_x)
      deallocate(ifax_y)
      deallocate(trigs_x)
      deallocate(trigs_y)
      deallocate(work_x)
      deallocate(work_y)
#endif

      deallocate(its2)
      deallocate(ite2)
      deallocate(itn2)
      deallocate(jts1)
      deallocate(jte1)
      deallocate(jtn1)

      deallocate(realbuf)
      deallocate(xfftbuf)
      deallocate(yfftbuf)

#ifdef MPI
      deallocate(sendbuf1)
      deallocate(sendcounts1)
      deallocate(sdispls1)

      deallocate(recvbuf1)
      deallocate(recvcounts1)
      deallocate(rdispls1)

      deallocate(sendbuf2)
      deallocate(sendcounts2)
      deallocate(sdispls2)

      deallocate(recvbuf2)
      deallocate(recvcounts2)
      deallocate(rdispls2)
#endif
! ==============================================================================
! === Finalize FFT End =========================================================
! ==============================================================================

      return
   end subroutine displacement_apply_kj_filter

#ifndef CARTESIAN
   subroutine displacement_calc_h0_lat1(dg, zz, nlon, nlat, h0, lat1)
#else
   subroutine displacement_calc_h0(dg, zz, nlon, nlat, h0)
#endif
      type(data_grids), target, intent(inout) :: dg 
      real(kind=8), dimension(nlon,nlat), intent(inout) :: zz
      integer(kind=4), intent(in) :: nlon, nlat
#ifndef CARTESIAN
      real(kind=8), intent(inout) :: h0, lat1
#else
      real(kind=8), intent(inout) :: h0
#endif

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dz
#ifdef QUAD
      real(kind=16) :: h0q
#ifdef MPI
      real(kind=16), allocatable, dimension(:) :: qbuf
#endif
#endif
! === get-large-area.f90 2015/06/19 ============================================
      real(kind=8), parameter :: disp_thres_ratio = 0.1d0
      real(kind=8) :: disp_max, disp_min
      integer(kind=4) :: imin, imax, jmin, jmax
! ==============================================================================
#ifdef MPI
      integer(kind=4) :: ix, iy, ixend, iyend, kx, kxend, ky, kyend, nx0, ny0, ibias, jbias, ierr
#endif
      real(kind=REAL_BYTE) :: x_inc, y_inc, x_min, x_max, y_min, y_max, z_min, z_max
      integer(kind=4) :: nlon_dummy, nlat_dummy, num, i, j

#ifdef MPI
      ix = dg%my%ix
      iy = dg%my%iy
      ixend = dg%my%ixend
      iyend = dg%my%iyend

      kx = dg%my%kx
      ky = dg%my%ky
      kxend = dg%my%kxend
      kyend = dg%my%kyend

      ! Space grid without edges.
      nx0 = ixend - ix + 1 ! xsize.
      ny0 = iyend - iy + 1 ! ysize.
      ibias = ix - kx ! x bias to skip edges.
      jbias = iy - ky ! y bias to skip edges.

      write(6,'(a,2i6)') '[displacement] nprocs, myrank: ', nprocs, myrank
      write(6,'(a,2i6)') '[displacement] ix, ixend: ', ix, ixend
      write(6,'(a,2i6)') '[displacement] iy, iyend: ', iy, iyend
      write(6,'(a,2i6)') '[displacement] nx0, ny0: ', nx0, ny0
#endif

      dz   => dg%depth_field%dz
#ifdef QUAD
#ifdef MPI
      allocate(qbuf(0:nprocs-1))
#endif
#endif

#ifdef MPI
      if(myrank == 0) then
#endif
      call read_gmt_grd_hdr(dg%my%bath_file,nlon_dummy,nlat_dummy,x_inc,y_inc,x_min,x_max,y_min,y_max,z_min,z_max)
#ifdef MPI
      end if
#endif
#ifdef MPI
#ifndef MULTI
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
#else
      call MPI_Bcast(x_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(y_inc, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
#endif
#endif

! === get-large-area.f90 2015/06/19 ============================================
      disp_max = -huge(disp_max)
      disp_min =  huge(disp_min)
!$omp parallel
!$omp do private(i) reduction(max:disp_max) reduction(min:disp_min)
      do j = 1, nlat
         do i = 1, nlon
            disp_max = max(disp_max,zz(i,j))
            disp_min = min(disp_min,zz(i,j))
         end do
      end do
!$omp single
#ifdef MPI
#ifndef MULTI
      call MPI_Allreduce(MPI_IN_PLACE, disp_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, disp_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(MPI_IN_PLACE, disp_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, disp_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_MEMBER_WORLD, ierr)
#endif
#endif

      imin =  huge(imin)
      imax = -huge(imax)
      jmin =  huge(jmin)
      jmax = -huge(jmax)
!$omp end single
      ! Case-1: all the disp is POSITIVE
      if(disp_max >= 0.0d0 .and. disp_min >= 0.0d0) then
!$omp do private(i) reduction(min:imin, jmin) reduction(max:imax, jmax)
         do j = 1, nlat
            do i = 1, nlon
               if(zz(i,j) >= disp_max*disp_thres_ratio) then
                  if(i <= imin) imin = i
                  if(i >= imax) imax = i
                  if(j <= jmin) jmin = j
                  if(j >= jmax) jmax = j
               end if
            end do
         end do
      ! Case-2: all the disp is NEGATIVE
      else if(disp_max <= 0.0d0 .and. disp_min <= 0.0d0) then
!$omp do private(i) reduction(min:imin, jmin) reduction(max:imax, jmax)
         do j = 1, nlat
            do i = 1, nlon
               if(zz(i,j) <= disp_min*disp_thres_ratio) then
                  if(i <= imin) imin = i
                  if(i >= imax) imax = i
                  if(j <= jmin) jmin = j
                  if(j >= jmax) jmax = j
               end if
            end do
         end do
      ! Case-3: disp includes both POSITIVE and NEGATIVE
      else
!$omp do private(i) reduction(min:imin, jmin) reduction(max:imax, jmax)
         do j = 1, nlat
            do i = 1, nlon
               if((zz(i,j) >= disp_max*disp_thres_ratio) .or. &
                  (zz(i,j) <= disp_min*disp_thres_ratio)) then
                  if(i <= imin) imin = i
                  if(i >= imax) imax = i
                  if(j <= jmin) jmin = j
                  if(j >= jmax) jmax = j
               end if
            end do
         end do
      end if
!$omp single
#ifndef MPI
#ifndef CARTESIAN
      lat1 = (dble(jmax) + dble(jmin))/2.0d0
      lat1 = y_min + y_inc*(lat1 - 1.0d0)
#endif
#else
      if(imin <  huge(imin)) imin = kx + imin - 1
      if(imax > -huge(imax)) imax = kx + imax - 1
      if(jmin <  huge(jmin)) jmin = ky + jmin - 1
      if(jmax > -huge(imax)) jmax = ky + jmax - 1
#ifndef MULTI
      call MPI_Allreduce(MPI_IN_PLACE, imin, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, imax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, jmin, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, jmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(MPI_IN_PLACE, imin, 1, MPI_INTEGER, MPI_MIN, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, imax, 1, MPI_INTEGER, MPI_MAX, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, jmin, 1, MPI_INTEGER, MPI_MIN, MPI_MEMBER_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, jmax, 1, MPI_INTEGER, MPI_MAX, MPI_MEMBER_WORLD, ierr)
#endif
#ifndef CARTESIAN
      lat1 = (dble(jmax) + dble(jmin))/2.0d0
      lat1 = y_min + y_inc*(lat1 - 1.0d0)
#ifndef MULTI
      call MPI_Allreduce(MPI_IN_PLACE, lat1, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(MPI_IN_PLACE, lat1, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_MEMBER_WORLD, ierr)
#endif
#endif
      imin = imin - kx + 1; imin = max(imin, ibias + 1)
      imax = imax - kx + 1; imax = min(imax, ibias + nx0)
      jmin = jmin - ky + 1; jmin = max(jmin, jbias + 1)
      jmax = jmax - ky + 1; jmax = min(jmax, jbias + ny0)
#endif
! ==============================================================================
#ifndef QUAD
      h0 = 0.0d0
#else
      h0q = 0.0q0
#endif
      num = 0
!$omp end single
#ifndef QUAD
!$omp do private(i) reduction(+:h0,num)
#else
!$omp do private(i) reduction(+:h0q,num)
#endif
! === get-large-area.f90 2015/06/19 ============================================
!     do j = 1, nlat
!        do i = 1, nlon
      do j = jmin, jmax
         do i = imin, imax
! ==============================================================================
            if(dz(i,j) > 0.0d0) then
#ifndef QUAD
               h0 = h0 + dz(i,j)
#else
               h0q = h0q + dz(i,j)
#endif
               num = num + 1
            end if
         end do
      end do
!$omp single
#ifdef MPI
#ifndef MULTI
#ifndef QUAD
      call MPI_Allreduce(MPI_IN_PLACE, h0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(h0q, 1, MPI_REAL16, qbuf, 1, MPI_REAL16, MPI_COMM_WORLD, ierr)
#endif
      call MPI_Allreduce(MPI_IN_PLACE, num, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
#ifndef QUAD
      call MPI_Allreduce(MPI_IN_PLACE, h0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#else
      call MPI_Allgather(h0q, 1, MPI_REAL16, qbuf, 1, MPI_REAL16, MPI_MEMBER_WORLD, ierr)
#endif
      call MPI_Allreduce(MPI_IN_PLACE, num, 1, MPI_INTEGER, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
#ifdef QUAD
      h0q = 0.0q0
      do i = 0, nprocs-1
         h0q = h0q + qbuf(i)
      end do
      deallocate(qbuf)
#endif
#endif
#ifdef QUAD
      h0 = h0q
#endif
      if(num == 0) then
         h0 = 0.0d0
      else
         h0 = h0/dble(num)
      end if

#ifndef CARTESIAN
      write(6,'(a,2f15.6)')  '[displacement] h0, lat1: ', h0, lat1
#else
      write(6,'(a,f15.6)')  '[displacement] h0: ', h0
#endif
!$omp end single
!$omp end parallel

      return
#ifndef CARTESIAN
   end subroutine displacement_calc_h0_lat1
#else
   end subroutine displacement_calc_h0
#endif

end module mod_displacement
