#include "real.h"
module mod_onefile
use mpi
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
use mod_grid
use mod_params
implicit none
integer(kind=4), private :: nprocs, myrank

contains

#ifdef ONEFILE
   subroutine onefile_initialize(nprocs_in, myrank_in)
      integer(kind=4), intent(in) :: nprocs_in, myrank_in
      nprocs = nprocs_in
      myrank = myrank_in
      return
   end subroutine onefile_initialize

   subroutine onefile_setparams(dg)
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), pointer :: kx, ky, kxend, kyend
      integer(kind=4), pointer :: ix, iy, ixend, iyend
      integer(kind=4), pointer :: srcount_x, srcount_y, srcount
      integer(kind=4), pointer, dimension(:) :: kx_all, ky_all, kxend_all, kyend_all
      integer(kind=4), pointer, dimension(:) :: ix_all, iy_all, ixend_all, iyend_all
      integer(kind=4) :: ierr

      kx        => dg%my%kx
      ky        => dg%my%ky
      kxend     => dg%my%kxend
      kyend     => dg%my%kyend
      ix        => dg%my%ix
      iy        => dg%my%iy
      ixend     => dg%my%ixend
      iyend     => dg%my%iyend
      srcount_x => dg%my%srcount_x
      srcount_y => dg%my%srcount_y
      srcount   => dg%my%srcount

      allocate(dg%my%kx_all(0:nprocs-1))
      allocate(dg%my%ky_all(0:nprocs-1))
      allocate(dg%my%kxend_all(0:nprocs-1))
      allocate(dg%my%kyend_all(0:nprocs-1))

      kx_all    => dg%my%kx_all
      ky_all    => dg%my%ky_all
      kxend_all => dg%my%kxend_all
      kyend_all => dg%my%kyend_all

#ifndef MULTI
      call MPI_Allgather(kx,    1, MPI_INTEGER, kx_all,    1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_Allgather(ky,    1, MPI_INTEGER, ky_all,    1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_Allgather(kxend, 1, MPI_INTEGER, kxend_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_Allgather(kyend, 1, MPI_INTEGER, kyend_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(kx,    1, MPI_INTEGER, kx_all,    1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
      call MPI_Allgather(ky,    1, MPI_INTEGER, ky_all,    1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
      call MPI_Allgather(kxend, 1, MPI_INTEGER, kxend_all, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
      call MPI_Allgather(kyend, 1, MPI_INTEGER, kyend_all, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      srcount_x = maxval(kxend_all - kx_all + 1)
      srcount_y = maxval(kyend_all - ky_all + 1)
      srcount   = srcount_x*srcount_y

      if(myrank == 0) then
         allocate(dg%my%buf_g(srcount*nprocs))
#ifdef NCDIO
         allocate(dg%my%buf_g_ncdio(srcount*nprocs))
#endif
      else
         allocate(dg%my%buf_g(1))
#ifdef NCDIO
         allocate(dg%my%buf_g_ncdio(1))
#endif
      end if

      allocate(dg%my%buf_l(srcount_x,srcount_y))
#ifdef NCDIO
      allocate(dg%my%buf_l_ncdio(srcount_x,srcount_y))
#endif

      allocate(dg%my%ix_all(0:nprocs-1))
      allocate(dg%my%iy_all(0:nprocs-1))
      allocate(dg%my%ixend_all(0:nprocs-1))
      allocate(dg%my%iyend_all(0:nprocs-1))

      ix_all    => dg%my%ix_all
      iy_all    => dg%my%iy_all
      ixend_all => dg%my%ixend_all
      iyend_all => dg%my%iyend_all

#ifndef MULTI
      call MPI_Allgather(ix,    1, MPI_INTEGER, ix_all,    1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_Allgather(iy,    1, MPI_INTEGER, iy_all,    1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_Allgather(ixend, 1, MPI_INTEGER, ixend_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_Allgather(iyend, 1, MPI_INTEGER, iyend_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(ix,    1, MPI_INTEGER, ix_all,    1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
      call MPI_Allgather(iy,    1, MPI_INTEGER, iy_all,    1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
      call MPI_Allgather(ixend, 1, MPI_INTEGER, ixend_all, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
      call MPI_Allgather(iyend, 1, MPI_INTEGER, iyend_all, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      return
   end subroutine onefile_setparams

   subroutine onefile_scatter_array(ain, aout, dg)
      real(kind=REAL_BYTE), intent(in),  dimension(:,:) :: ain
      real(kind=REAL_BYTE), intent(out), dimension(:,:) :: aout
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), pointer :: nx, ny, srcount_x, srcount
      integer(kind=4), pointer, dimension(:) :: kx_all, ky_all, kxend_all, kyend_all
      real(kind=REAL_BYTE), pointer, dimension(:) :: buf_g
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: buf_l
      integer(kind=4) :: i, j, p, ind, ierr

      nx        => dg%my%nx
      ny        => dg%my%ny
      srcount_x => dg%my%srcount_x
      srcount   => dg%my%srcount

      kx_all    => dg%my%kx_all
      ky_all    => dg%my%ky_all
      kxend_all => dg%my%kxend_all
      kyend_all => dg%my%kyend_all

      buf_g     => dg%my%buf_g
      buf_l     => dg%my%buf_l

      if(myrank == 0) then
         do p = 0, nprocs - 1
            do j = ky_all(p), kyend_all(p)
               do i = kx_all(p), kxend_all(p)
                 ind = srcount*p + (i - kx_all(p) + 1) + (j - ky_all(p))*srcount_x
                 buf_g(ind) = ain(i,j)
               end do
            end do
         end do
      end if

#ifndef MULTI
      call MPI_Scatter(buf_g, srcount, REAL_MPI, buf_l, srcount, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
#else
      call MPI_Scatter(buf_g, srcount, REAL_MPI, buf_l, srcount, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
#endif

      do j = 1, ny
         do i = 1, nx
            aout(i,j) = buf_l(i,j)
         end do
      end do

      return
   end subroutine onefile_scatter_array

   subroutine onefile_gather_array(ain, aout, dg)
      real(kind=REAL_BYTE), intent(in),  dimension(:,:) :: ain
      real(kind=REAL_BYTE), intent(out), dimension(:,:) :: aout
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), pointer :: nx, ny, srcount_x, srcount
      integer(kind=4), pointer, dimension(:) :: ix_all, iy_all, ixend_all, iyend_all
      integer(kind=4), pointer, dimension(:) :: kx_all, ky_all
      real(kind=REAL_BYTE), pointer, dimension(:) :: buf_g
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: buf_l
      integer(kind=4) :: i, j, p, ind, ierr

      nx        => dg%my%nx
      ny        => dg%my%ny
      srcount_x => dg%my%srcount_x
      srcount   => dg%my%srcount

      ix_all    => dg%my%ix_all
      iy_all    => dg%my%iy_all
      ixend_all => dg%my%ixend_all
      iyend_all => dg%my%iyend_all

      kx_all    => dg%my%kx_all
      ky_all    => dg%my%ky_all

      buf_g     => dg%my%buf_g
      buf_l     => dg%my%buf_l

      do j = 1, ny
         do i = 1, nx
            buf_l(i,j) = ain(i,j)
         end do
      end do

#ifndef MULTI
      call MPI_Gather(buf_l, srcount, REAL_MPI, buf_g, srcount, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
#else
      call MPI_Gather(buf_l, srcount, REAL_MPI, buf_g, srcount, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
#endif

      if(myrank == 0) then
         do p = 0, nprocs - 1
            do j = iy_all(p), iyend_all(p)
               do i = ix_all(p), ixend_all(p)
                 ind = srcount*p + (i - kx_all(p) + 1) + (j - ky_all(p))*srcount_x
                 aout(i,j) = buf_g(ind)
               end do
            end do
         end do
      end if

      return
   end subroutine onefile_gather_array
#ifdef NCDIO
   subroutine onefile_gather_array_ncdio(ain, aout, dg)
      real(kind=4), intent(in),  dimension(:,:) :: ain
      real(kind=4), intent(out), dimension(:,:) :: aout
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), pointer :: nx, ny, srcount_x, srcount, totalNy
      integer(kind=4), pointer, dimension(:) :: ix_all, iy_all, ixend_all, iyend_all
      integer(kind=4), pointer, dimension(:) :: kx_all, ky_all
      real(kind=4), pointer, dimension(:) :: buf_g
      real(kind=4), pointer, dimension(:,:) :: buf_l
      integer(kind=4) :: i, j, p, ind, ierr

      nx        => dg%my%nx
      ny        => dg%my%ny
      totalNy   => dg%my%totalNy
      srcount_x => dg%my%srcount_x
      srcount   => dg%my%srcount

      ix_all    => dg%my%ix_all
      iy_all    => dg%my%iy_all
      ixend_all => dg%my%ixend_all
      iyend_all => dg%my%iyend_all

      kx_all    => dg%my%kx_all
      ky_all    => dg%my%ky_all

      buf_g     => dg%my%buf_g_ncdio
      buf_l     => dg%my%buf_l_ncdio

      do j = 1, ny
         do i = 1, nx
            buf_l(i,j) = ain(i,j)
         end do
      end do

#ifndef MULTI
      call MPI_Gather(buf_l, srcount, MPI_REAL4, buf_g, srcount, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
#else
      call MPI_Gather(buf_l, srcount, MPI_REAL4, buf_g, srcount, MPI_REAL4, 0, MPI_MEMBER_WORLD, ierr)
#endif

      if(myrank == 0) then
         do p = 0, nprocs - 1
            do j = iy_all(p), iyend_all(p)
               do i = ix_all(p), ixend_all(p)
                 ind = srcount*p + (i - kx_all(p) + 1) + (j - ky_all(p))*srcount_x
                 aout(i,totalNy-j+1) = buf_g(ind)
               end do
            end do
         end do
      end if

      return
   end subroutine onefile_gather_array_ncdio
#endif
#else
#ifdef __SX__
   subroutine onefile_dummy()
   end subroutine onefile_dummy
#endif
#endif

end module mod_onefile
