#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_dump1d
#ifdef MPI
use mpi
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
#endif
use mod_grid
implicit none

integer(kind=4) :: target_domain = 1
integer(kind=4) :: dump_interval = 100
integer(kind=4) :: space_interval = 2
integer(kind=4) :: target_j = -1 ! -1: none, 0: center
integer(kind=4) :: target_i = -1 ! -1: none, 0: center

integer(kind=4) :: nx, ny
#ifdef MPI
integer(kind=4) :: ix, ixend, iy, iyend, kx, ky
integer(kind=4) :: nprocs, myrank

real(kind=REAL_BYTE), allocatable, dimension(:) :: tmphzx, tmpdzx
real(kind=REAL_BYTE), allocatable, dimension(:) :: tmphzy, tmpdzy
#endif
real(kind=REAL_BYTE) :: dh, lon0, lat0

namelist /params/ target_domain, dump_interval, space_interval, target_j, target_i

contains

   subroutine dump1d_initialize(dgrid)
      type(data_grids), intent(inout), target, dimension(:) :: dgrid
#ifdef MPI
      integer(kind=4) :: ierr
#endif

      open(1,file='dump1d',action='read',status='old',form='formatted',err=100)
      read(1,params)
      close(1)
      go to 200
 100  continue
      write(6,'(a)') 'No dump1d file exists. Default parameters are utilized.'
 200  continue

#ifndef MPI
      nx = dgrid(target_domain)%my%nx
      ny = dgrid(target_domain)%my%ny
#else
      nx    = dgrid(target_domain)%my%totalNx
      ny    = dgrid(target_domain)%my%totalNy
      ix    = dgrid(target_domain)%my%ix
      ixend = dgrid(target_domain)%my%ixend
      iy    = dgrid(target_domain)%my%iy
      iyend = dgrid(target_domain)%my%iyend
      kx    = dgrid(target_domain)%my%kx
      ky    = dgrid(target_domain)%my%ky

#ifndef MULTI
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
#else
      call MPI_Comm_size(MPI_MEMBER_WORLD, nprocs, ierr)
      call MPI_Comm_rank(MPI_MEMBER_WORLD, myrank, ierr)
#endif
#endif
      dh = dgrid(target_domain)%my%dh

#ifndef MPI
      lon0 = dgrid(target_domain)%my%mlon0
      lat0 = dgrid(target_domain)%my%mlat0
#else
      lon0 = dgrid(target_domain)%my%glon0
      lat0 = dgrid(target_domain)%my%glat0
#endif
#ifndef CARTESIAN
      lon0  = lon0/60.0d0
      if(lon0 > 180.0d0) lon0 = lon0 - 360.0d0

      lat0 = 90.0d0 - (lat0/60.0d0 + (ny-1.0d0)*dh)
#endif

      if(target_j == 0) target_j = ny/2
      if(target_i == 0) target_i = nx/2

#ifdef MPI
      if(myrank == 0) then
#endif
      write(6,'(a)') '========== DUMP1D parameters ==============================='
      write(6,'(a,i6)') 'Target domain:   ', target_domain
      write(6,'(a,i6)') 'Interval[steps]: ', dump_interval
      write(6,'(a,i6)') 'Interval[grids]: ', space_interval
      write(6,'(a,i6)') 'j-index:         ', target_j
      write(6,'(a,i6)') 'i-index:         ', target_i
      write(6,'(a,e15.6)') 'dx/dy:           ', dh
      write(6,'(a,e15.6)') 'lon0:            ', lon0
      write(6,'(a,e15.6)') 'lat0:            ', lat0
      write(6,'(a)') '============================================================'
#ifdef MPI
      end if
#endif

#ifdef MPI
      if(target_j /= -1) then
         allocate(tmphzx(nx))
         allocate(tmpdzx(nx))
      end if

      if(target_i /= -1) then
         allocate(tmphzy(ny))
         allocate(tmpdzy(ny))
      end if
#endif

      return
   end subroutine dump1d_initialize

   subroutine dump1d_finalize()
#ifdef MPI
      if(target_j /= -1) then
         deallocate(tmphzx)
         deallocate(tmpdzx)
      end if

      if(target_i /= -1) then
         deallocate(tmphzy)
         deallocate(tmpdzy)
      end if
#endif
      return
   end subroutine dump1d_finalize

   subroutine dump1d(ig, dg, istep)
      integer(kind=4), intent(in) :: ig, istep
      type(data_grids), target, intent(inout) :: dg

      character(len=512) :: str
      integer(kind=4) :: i, j
#ifdef MPI
      integer(kind=4) :: ierr
#endif

      if((ig /= target_domain) .or. (mod(istep,dump_interval) /= 0)) return

      if(target_j /= -1) then
!        write(str,'(a,a,i8.8)') trim(members_dir), 'hzx.', istep
         write(str,'(a,i8.8)') 'hzx.', istep
#ifndef MPI
         open(998,file=trim(str),action='write',status='replace',form='formatted')
         do i = 1, nx, space_interval
            write(998,'(3f15.6)') lon0 + dh*(i-1), dg%wave_field%hz(i,target_j), -dg%depth_field%dz(i,target_j)
         end do
         close(998)
#else
         tmphzx = 0.0d0
         tmpdzx = 0.0d0
         if((target_j >= iy) .and. (target_j <= iyend)) then
            do i = max(1,ix), min(nx,ixend)
               tmphzx(i) = dg%wave_field%hz(i-kx+1,target_j-ky+1)
               tmpdzx(i) = dg%depth_field%dz(i-kx+1,target_j-ky+1)
            end do
         end if
#ifndef MULTI
         call MPI_Allreduce(MPI_IN_PLACE,tmphzx,nx,REAL_MPI,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_Allreduce(MPI_IN_PLACE,tmpdzx,nx,REAL_MPI,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
         call MPI_Allreduce(MPI_IN_PLACE,tmphzx,nx,REAL_MPI,MPI_SUM,MPI_MEMBER_WORLD,ierr)
         call MPI_Allreduce(MPI_IN_PLACE,tmpdzx,nx,REAL_MPI,MPI_SUM,MPI_MEMBER_WORLD,ierr)
#endif
         if(myrank == 0) then
            open(998,file=trim(str),action='write',status='replace',form='formatted')
            do i = 1, nx, space_interval
               write(998,'(3f15.6)') lon0 + dh*(i-1), tmphzx(i), -tmpdzx(i)
            end do
            close(998)
         end if
#endif
      end if

      if(target_i /= -1) then
!        write(str,'(a,a,i8.8)') trim(members_dir), 'hzy.', istep
         write(str,'(a,i8.8)') 'hzy.', istep
#ifndef MPI
         open(999,file=trim(str),action='write',status='replace',form='formatted')
         do j = 1, ny, space_interval
            write(999,'(3f15.6)') lat0 + dh*(j-1), dg%wave_field%hz(target_i,ny-j+1), -dg%depth_field%dz(target_i,ny-j+1)
         end do
         close(999)
#else
         tmphzy = 0.0d0
         tmpdzy = 0.0d0
         if((target_i >= ix) .and. (target_i <= ixend)) then
            do j = max(1,iy), min(ny,iyend)
               tmphzy(j) = dg%wave_field%hz(target_i-kx+1,j-ky+1)
               tmpdzy(j) = dg%depth_field%dz(target_i-kx+1,j-ky+1)
            end do
         end if
#ifndef MULTI
         call MPI_Allreduce(MPI_IN_PLACE,tmphzy,ny,REAL_MPI,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_Allreduce(MPI_IN_PLACE,tmpdzy,ny,REAL_MPI,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
         call MPI_Allreduce(MPI_IN_PLACE,tmphzy,ny,REAL_MPI,MPI_SUM,MPI_MEMBER_WORLD,ierr)
         call MPI_Allreduce(MPI_IN_PLACE,tmpdzy,ny,REAL_MPI,MPI_SUM,MPI_MEMBER_WORLD,ierr)
#endif
         if(myrank == 0) then
            open(999,file=trim(str),action='write',status='replace',form='formatted')
            do j = 1, ny, space_interval
               write(999,'(3f15.6)') lat0 + dh*(j-1), tmphzy(ny-j+1), -tmpdzy(ny-j+1)
            end do
            close(999)
         end if
#endif
      end if

      return
   end subroutine dump1d

end module mod_dump1d
