#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_interpolation
#ifdef MPI
use mpi
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
#endif
use mod_grid
implicit none
#ifdef MPI
integer(kind=4), private :: nprocs, myrank, npx, npy, rankx, ranky
integer(kind=4), private :: MPI_X_WORLD, MPI_Y_WORLD
#endif

private :: spline_inerpolation, tri, linear_inerpolation

contains

#ifdef  MPI
   subroutine interpolation_mpi_initialize(nprocs_in, myrank_in, npx_in, npy_in, rankx_in, ranky_in)
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
   end subroutine interpolation_mpi_initialize
#endif

   subroutine interp2fine_init_disp(cg,fg)
      use mod_params, only : use_linear
      type(data_grids), target, intent(inout) :: cg, fg
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: buf0, buf1, buf2
      integer(kind=4) :: bigNX, bigNY, zeroIX, zeroIY, nx, ny
      integer(kind=4) :: i, j
#ifdef MPI
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: buf3, tmp2d
      real(kind=REAL_BYTE), allocatable, dimension(:) :: tmp1d
      integer(kind=4), allocatable, dimension(:) :: sendcounts0, sdispls0 ! MPI_Scatterv
      integer(kind=4), allocatable, dimension(:) :: sendcounts1, sdispls1 ! MPI_Alltoallv (1)
      integer(kind=4), allocatable, dimension(:) :: recvcounts1, rdispls1 ! MPI_Alltoallv (1)
      integer(kind=4), allocatable, dimension(:) :: sendcounts2, sdispls2 ! MPI_Alltoallv (2)
      integer(kind=4), allocatable, dimension(:) :: recvcounts2, rdispls2 ! MPI_Alltoallv (2)
      integer(kind=4) :: totalNx, totalNy, ix, iy, ixend, iyend, kx, ky
      integer(kind=4) :: nx0, ny1
      integer(kind=4) :: ic, jc, p, ierr, ist, ien, jst, jen, ind, nbx, nby, ib, jb
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: commbuf
! ==============================================================================
#endif

      bigNX = fg%my%bigNX
      bigNY = fg%my%bigNY
      zeroIX = fg%my%zeroIX
      zeroIY = fg%my%zeroIY
      nx = fg%my%nx
      ny = fg%my%ny
#ifdef MPI
      totalNx = fg%my%totalNx
      totalNy = fg%my%totalNy

      ix = fg%my%ix
      iy = fg%my%iy
      ixend = fg%my%ixend
      iyend = fg%my%iyend
      kx = fg%my%kx
      ky = fg%my%ky

      nbx = ixend - ix + 1
      nby = iyend - iy + 1

      ! Scatter tmp2d(bigNY,bigNX) -> buf0(bigNY,nx0)
      nx0 = bigNX/nprocs
      if(myrank < mod(bigNX,nprocs)) nx0 = nx0 + 1

      allocate(sendcounts0(0:nprocs-1))
      allocate(sdispls0(0:nprocs-1))
#ifndef MULTI
      call MPI_Allgather(nx0,1,MPI_INTEGER,sendcounts0,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
#else
      call MPI_Allgather(nx0,1,MPI_INTEGER,sendcounts0,1,MPI_INTEGER,MPI_MEMBER_WORLD,ierr)
#endif
      do p = 0, nprocs-1
         sendcounts0(p) = sendcounts0(p)*bigNY
      end do
      sdispls0(0) = 0
      do p = 1, nprocs-1
         sdispls0(p) = sdispls0(p-1) + sendcounts0(p-1)
      end do

      ! Transposition tmp2d(nx0,totalNy) -> tmp1d(bigNX*ny1)
      ny1 = nby/npx
      if(rankx < mod(nby,npx)) ny1 = ny1 + 1

      allocate(sendcounts1(0:nprocs-1))
      allocate(sdispls1(0:nprocs-1))
#ifndef MULTI
      call MPI_Allgather(ny1,1,MPI_INTEGER,sendcounts1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
#else
      call MPI_Allgather(ny1,1,MPI_INTEGER,sendcounts1,1,MPI_INTEGER,MPI_MEMBER_WORLD,ierr)
#endif
      do p = 0, nprocs-1
         sendcounts1(p) = sendcounts1(p)*nx0
      end do
      sdispls1(0) = 0
      do p = 1, nprocs-1
         sdispls1(p) = sdispls1(p-1) + sendcounts1(p-1)
      end do

      allocate(recvcounts1(0:nprocs-1))
      allocate(rdispls1(0:nprocs-1))
      call MPI_Alltoall(sendcounts1,1,MPI_INTEGER, &
#ifndef MULTI
                        recvcounts1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
#else
                        recvcounts1,1,MPI_INTEGER,MPI_MEMBER_WORLD,ierr)
#endif
      rdispls1(0) = 0
      do p = 1, nprocs-1
         rdispls1(p) = rdispls1(p-1) + recvcounts1(p-1)
      end do

      ! Transposition tmp2d(ny1,totalNx) -> tmp1d(nby*nbx)
      allocate(sendcounts2(0:npx-1))
      allocate(sdispls2(0:npx-1))
      call MPI_Allgather(nbx,1,MPI_INTEGER,sendcounts2,1,MPI_INTEGER,MPI_X_WORLD,ierr)
      do p = 0, npx-1
         sendcounts2(p) = sendcounts2(p)*ny1
      end do
      sdispls2(0) = 0
      do p = 1, npx-1
         sdispls2(p) = sdispls2(p-1) + sendcounts2(p-1)
      end do

      allocate(recvcounts2(0:nprocs-1))
      allocate(rdispls2(0:nprocs-1))
      call MPI_Alltoall(sendcounts2,1,MPI_INTEGER, &
                        recvcounts2,1,MPI_INTEGER,MPI_X_WORLD,ierr)
      rdispls2(0) = 0
      do p = 1, npx-1
         rdispls2(p) = rdispls2(p-1) + recvcounts2(p-1)
      end do
#endif

!$omp parallel
#ifndef MPI
!$omp single
      allocate(buf0(bigNY,bigNX))
      allocate(buf1(ny,   bigNX))
      allocate(buf2(bigNX,ny))
!$omp end single

!$omp do private(j)
      do i = 1, bigNX
         do j = 1, bigNY
            buf0(j,i) = cg%zz(i+zeroIX-1,j+zeroIY-1)
         end do
      end do

      ! Y interpolation
      if(use_linear == 1) then
!$omp do
         do i = 1, bigNX
            call linear_inerpolation(fg%my%nr, bigNY, buf0(1,i), ny, buf1(1,i))
         end do
      else
!$omp do
         do i = 1, bigNX
            call spline_inerpolation(fg%my%nr, bigNY, buf0(1,i), ny, buf1(1,i))
         end do
      end if
!$omp do private(i)
      do j = 1, ny
         do i = 1, bigNX
            buf2(i,j) = buf1(j,i)
         end do
      end do
      ! X interpolation
      if(use_linear == 1) then
!$omp do
         do j = 1, ny
            call linear_inerpolation(fg%my%nr, bigNX, buf2(1,j), nx, fg%zz(1,j))
         end do
      else
!$omp do
         do j = 1, ny
            call spline_inerpolation(fg%my%nr, bigNX, buf2(1,j), nx, fg%zz(1,j))
         end do
      end if

!$omp single
      deallocate(buf0)
      deallocate(buf1)
      deallocate(buf2)
!$omp end single
#else
      ! Collect all coarse elements to rank 0.
!$omp single
      allocate(tmp2d(bigNY,bigNX))
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
      allocate(commbuf(bigNY,bigNX))
! ==============================================================================
!$omp end single
!$omp do private(j)
      do i = 1, bigNX
         do j = 1, bigNY
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
!           tmp2d(j,i) = 0.0d0
            commbuf(j,i) = 0.0d0
! ==============================================================================
         end do
      end do
!$omp do private(ic, j, jc)
      do i = 1, bigNX
         ic = i + zeroIX - 1
         if((ic >= cg%my%ix) .and. (ic <= cg%my%ixend)) then
            do j = 1, bigNY
               jc = j + zeroIY - 1
               if((jc >= cg%my%iy) .and. (jc <= cg%my%iyend)) then
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
!                 tmp2d(j,i) = cg%zz(ic-cg%my%kx+1,jc-cg%my%ky+1)
                  commbuf(j,i) = cg%zz(ic-cg%my%kx+1,jc-cg%my%ky+1)
! ==============================================================================
               end if
            end do
         end if
      end do
!$omp single
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
!      if(myrank == 0) then
!         call MPI_Reduce(MPI_IN_PLACE, tmp2d, bigNY*bigNX, REAL_MPI, &
!#ifndef MULTI
!                         MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!#else
!                         MPI_SUM, 0, MPI_MEMBER_WORLD, ierr)
!#endif
!      else
!         call MPI_Reduce(tmp2d, tmp2d, bigNY*bigNX, REAL_MPI, &
!#ifndef MULTI
!                         MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!#else
!                         MPI_SUM, 0, MPI_MEMBER_WORLD, ierr)
!#endif
!      end if
      call MPI_Reduce(commbuf, tmp2d, bigNY*bigNX, REAL_MPI, &
#ifndef MULTI
                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
#else
                      MPI_SUM, 0, MPI_MEMBER_WORLD, ierr)
#endif
      deallocate(commbuf)
! ==============================================================================

      ! Scatter tmp2d(bigNY,bigNX) -> buf0(bigNY,nx0)
      allocate(buf0(bigNY,nx0))
      call MPI_Scatterv(tmp2d, sendcounts0, sdispls0, REAL_MPI, &
                        buf0,   bigNY*nx0,             REAL_MPI, &
#ifndef MULTI
                        0, MPI_COMM_WORLD, ierr)
#else
                        0, MPI_MEMBER_WORLD, ierr)
#endif
      deallocate(tmp2d)

      ! Y interpolation buf0(bigNY,nx0) -> buf1(totalNy,nx0)
      allocate(buf1(totalNy,nx0))
!$omp end single
      if(use_linear == 1) then
!$omp do
         do i = 1, nx0
            call linear_inerpolation(fg%my%nr, bigNY, buf0(1,i), totalNy, buf1(1,i))
         end do
      else
!$omp do
         do i = 1, nx0
            call spline_inerpolation(fg%my%nr, bigNY, buf0(1,i), totalNy, buf1(1,i))
         end do
      end if
!$omp single
      deallocate(buf0)

      ! Local transposition buf1(totalNy,nx0) -> tmp2d(nx0,totalNy)
      allocate(tmp2d(nx0,totalNy))
!$omp end single
!$omp do private(i)
      do j = 1, totalNy
         do i = 1, nx0
            tmp2d(i,j) = buf1(j,i)
         end do
      end do
!$omp single
      deallocate(buf1)
      ! Transposition tmp2d(nx0,totalNy) -> tmp1d(bigNX*ny1)
      allocate(tmp1d(bigNX*ny1))
      call MPI_Alltoallv(tmp2d, sendcounts1, sdispls1, REAL_MPI, &
                         tmp1d, recvcounts1, rdispls1, REAL_MPI, &
#ifndef MULTI
                         MPI_COMM_WORLD, ierr)
#else
                         MPI_MEMBER_WORLD, ierr)
#endif
      deallocate(tmp2d)
      ! Rearrange tmp1d(bigNX*ny1) -> buf2(bigNX,ny1)
      allocate(buf2(bigNX,ny1))
!$omp end single
      do p = 0, nprocs-1
         ist = rdispls1(p)/ny1 + 1
         ien = ist + recvcounts1(p)/ny1 - 1
!$omp do private(i, ind)
         do j = 1, ny1
            do i = ist, ien
               ind = rdispls1(p) + (ien - ist + 1)*(j - 1) + i - ist + 1
               buf2(i,j) = tmp1d(ind)
            end do
         end do
      end do
!$omp single
      deallocate(tmp1d)

      ! X interpolation buf2(i,j) -> buf3(totalNx,ny1)
      allocate(buf3(totalNx,ny1))
!$omp end single
      if(use_linear == 1) then
!$omp do
         do j = 1, ny1
            call linear_inerpolation(fg%my%nr, bigNX, buf2(1,j), totalNx, buf3(1,j))
         end do
      else
!$omp do
         do j = 1, ny1
            call spline_inerpolation(fg%my%nr, bigNX, buf2(1,j), totalNx, buf3(1,j))
         end do
      end if
!$omp single
      deallocate(buf2)

      ! Local transposition buf3(totalNx,ny1) -> tmp2d(ny1,totalNx)
      allocate(tmp2d(ny1,totalNx))
!$omp end single
!$omp do private(j)
      do i = 1, totalNx
         do j = 1, ny1
            tmp2d(j,i) = buf3(i,j)
         end do
      end do
!$omp single
      deallocate(buf3)
      ! Transposition tmp2d(ny1,totalNx) -> tmp1d(nby*nbx)
      allocate(tmp1d(nby*nbx))
      call MPI_Alltoallv(tmp2d, sendcounts2, sdispls2, REAL_MPI, &
                         tmp1d, recvcounts2, rdispls2, REAL_MPI, &
                         MPI_X_WORLD, ierr)
      deallocate(tmp2d)
!$omp end single

      ! Finally, rearrange tmp1d(nby*nbx) -> fg%zz(nx,ny)
      ib = ix - kx
      jb = iy - ky
      do p = 0, npx-1
         jst = rdispls2(p)/nbx + 1
         jen = jst + recvcounts2(p)/nbx - 1
!$omp do private(j, ind)
         do i = 1, nbx
            do j = jst, jen
               ind = rdispls2(p) + (jen - jst + 1)*(i - 1) + j - jst + 1
               fg%zz(i+ib,j+jb) = tmp1d(ind)
            end do
         end do
      end do
!$omp single
      deallocate(tmp1d)

      deallocate(sendcounts0)
      deallocate(sdispls0)
      deallocate(sendcounts1)
      deallocate(sdispls1)
      deallocate(recvcounts1)
      deallocate(rdispls1)
      deallocate(sendcounts2)
      deallocate(sdispls2)
      deallocate(recvcounts2)
      deallocate(rdispls2)
!$omp end single
#endif
!$omp end parallel
      return
   end subroutine interp2fine_init_disp 

   subroutine spline_inerpolation(r, n_in, y, n_out, yo)
      ! Arguments
      integer(kind=4), intent(in) :: r, n_in, n_out
      real(kind=REAL_BYTE), intent(in), dimension(n_in) :: y
      real(kind=REAL_BYTE), intent(out), dimension(n_out) :: yo
      ! Work
      real(kind=REAL_BYTE), dimension(n_in-2)   :: d, b
      real(kind=REAL_BYTE), dimension(n_in)     :: u
      real(kind=REAL_BYTE), dimension(n_in-1,4) :: c
      real(kind=REAL_BYTE) :: x0
      integer(kind=4) :: n, i, j, i_

      n = n_in - 2

      ! Solve ax = b
      ! a is tridiagonal matrix (dl, d, du)
      do i = 1, n
         d(i) = 2.0d0*2.0d0
         b(i) = 6.0d0*(y(i+2) - 2.0d0*y(i+1) + y(i))
      end do

      call tri(n, d, b)

      ! Make u (= d^2y/dx^2)
      u(1) = 0.0d0
      do i = 1, n
         u(i+1) = b(i)
      end do
      u(n_in) = 0.0d0

      ! Calc. coefficients
      ! S = c(i,1)*x^3 + c(i,2)*x^2 + c(i,3)*x + c(i,4)
      do i = 1, n_in-1
         c(i,1) = (u(i+1) - u(i))/6.0d0
         c(i,2) = u(i)/2.0d0
         c(i,3) = (y(i+1) - y(i)) - (2.0d0*u(i) + u(i+1))/6.0d0
         c(i,4) = y(i)
      end do

#ifdef DEBUG
      write(6,'(a)') '--------------------------------------------------------------------'
      write(6,'(a)') 'i       a              b              c              d              '
      write(6,'(a)') '--------------------------------------------------------------------'
      do i = 1, n_in-1
         write(6,'(i8,4f15.6)') i, c(i,1), c(i,2), c(i,3), c(i,4)
      end do
      write(6,'(a)') '--------------------------------------------------------------------'
#endif

      ! Interpolation
      do j = 1, n_in-1
         do i = 0, r-1
            i_ = (j-1)*r + i + 1
            x0 = dble(i)/dble(r)
            yo(i_) = c(j,1)*x0**3 + c(j,2)*x0**2 + c(j,3)*x0 + c(j,4)
         end do
      end do
      yo(n_out) = y(n_in)

      return
   end subroutine spline_inerpolation

   subroutine tri(n, d, b)
      implicit none
      integer(kind=4), intent(in) :: n
      real(kind=REAL_BYTE), intent(inout), dimension(n) :: d, b
      integer(kind=4) :: i
      do i = 1, n-1
         d(i+1) = d(i+1) - 1.0d0/d(i)
      end do
      do i = 1, n-1
         b(i+1) = b(i+1) - b(i)/d(i)
      end do
      b(n) = b(n)/d(n)
      do i = n-1, 1, -1
         b(i) = (b(i) - b(i+1))/d(i)
      end do
      return
   end subroutine tri

   subroutine linear_inerpolation(r, n_in, y, n_out, yo)
      ! Arguments
      integer(kind=4), intent(in) :: r, n_in, n_out
      real(kind=REAL_BYTE), intent(in), dimension(n_in) :: y
      real(kind=REAL_BYTE), intent(out), dimension(n_out) :: yo
      real(kind=REAL_BYTE) :: x0, x1
      integer(kind=4) :: i, j, i_
      ! Interpolation
      do j = 1, n_in-1
         do i = 0, r-1
            i_ = (j-1)*r + i + 1
            x0 = dble(r-i)/dble(r)
            x1 = dble(i)/dble(r)
            yo(i_) = x0*y(j) + x1*y(j+1)
         end do
      end do
      yo(n_out) = y(n_in)
      return
   end subroutine linear_inerpolation

! === Elastic loading with interpolation =======================================
!  subroutine interp2fine_init_disp(cg,fg)
#ifndef CARTESIAN
   subroutine interp2fine_elastic_loading(cg,fg)
! ==============================================================================
      use mod_params, only : use_linear
      type(data_grids), target, intent(inout) :: cg, fg
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: buf0, buf1, buf2
      integer(kind=4) :: bigNX, bigNY, zeroIX, zeroIY, nx, ny
      integer(kind=4) :: i, j
#ifdef MPI
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: buf3, tmp2d
      real(kind=REAL_BYTE), allocatable, dimension(:) :: tmp1d
      integer(kind=4), allocatable, dimension(:) :: sendcounts0, sdispls0 ! MPI_Scatterv
      integer(kind=4), allocatable, dimension(:) :: sendcounts1, sdispls1 ! MPI_Alltoallv (1)
      integer(kind=4), allocatable, dimension(:) :: recvcounts1, rdispls1 ! MPI_Alltoallv (1)
      integer(kind=4), allocatable, dimension(:) :: sendcounts2, sdispls2 ! MPI_Alltoallv (2)
      integer(kind=4), allocatable, dimension(:) :: recvcounts2, rdispls2 ! MPI_Alltoallv (2)
      integer(kind=4) :: totalNx, totalNy, ix, iy, ixend, iyend, kx, ky
      integer(kind=4) :: nx0, ny1
      integer(kind=4) :: ic, jc, p, ierr, ist, ien, jst, jen, ind, nbx, nby, ib, jb
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: commbuf
! ==============================================================================
#endif

      bigNX = fg%my%bigNX
      bigNY = fg%my%bigNY
      zeroIX = fg%my%zeroIX
      zeroIY = fg%my%zeroIY
      nx = fg%my%nx
      ny = fg%my%ny
#ifdef MPI
      totalNx = fg%my%totalNx
      totalNy = fg%my%totalNy

      ix = fg%my%ix
      iy = fg%my%iy
      ixend = fg%my%ixend
      iyend = fg%my%iyend
      kx = fg%my%kx
      ky = fg%my%ky

      nbx = ixend - ix + 1
      nby = iyend - iy + 1

      ! Scatter tmp2d(bigNY,bigNX) -> buf0(bigNY,nx0)
      nx0 = bigNX/nprocs
      if(myrank < mod(bigNX,nprocs)) nx0 = nx0 + 1

      allocate(sendcounts0(0:nprocs-1))
      allocate(sdispls0(0:nprocs-1))
#ifndef MULTI
      call MPI_Allgather(nx0,1,MPI_INTEGER,sendcounts0,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
#else
      call MPI_Allgather(nx0,1,MPI_INTEGER,sendcounts0,1,MPI_INTEGER,MPI_MEMBER_WORLD,ierr)
#endif
      do p = 0, nprocs-1
         sendcounts0(p) = sendcounts0(p)*bigNY
      end do
      sdispls0(0) = 0
      do p = 1, nprocs-1
         sdispls0(p) = sdispls0(p-1) + sendcounts0(p-1)
      end do

      ! Transposition tmp2d(nx0,totalNy) -> tmp1d(bigNX*ny1)
      ny1 = nby/npx
      if(rankx < mod(nby,npx)) ny1 = ny1 + 1

      allocate(sendcounts1(0:nprocs-1))
      allocate(sdispls1(0:nprocs-1))
#ifndef MULTI
      call MPI_Allgather(ny1,1,MPI_INTEGER,sendcounts1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
#else
      call MPI_Allgather(ny1,1,MPI_INTEGER,sendcounts1,1,MPI_INTEGER,MPI_MEMBER_WORLD,ierr)
#endif
      do p = 0, nprocs-1
         sendcounts1(p) = sendcounts1(p)*nx0
      end do
      sdispls1(0) = 0
      do p = 1, nprocs-1
         sdispls1(p) = sdispls1(p-1) + sendcounts1(p-1)
      end do

      allocate(recvcounts1(0:nprocs-1))
      allocate(rdispls1(0:nprocs-1))
      call MPI_Alltoall(sendcounts1,1,MPI_INTEGER, &
#ifndef MULTI
                        recvcounts1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
#else
                        recvcounts1,1,MPI_INTEGER,MPI_MEMBER_WORLD,ierr)
#endif
      rdispls1(0) = 0
      do p = 1, nprocs-1
         rdispls1(p) = rdispls1(p-1) + recvcounts1(p-1)
      end do

      ! Transposition tmp2d(ny1,totalNx) -> tmp1d(nby*nbx)
      allocate(sendcounts2(0:npx-1))
      allocate(sdispls2(0:npx-1))
      call MPI_Allgather(nbx,1,MPI_INTEGER,sendcounts2,1,MPI_INTEGER,MPI_X_WORLD,ierr)
      do p = 0, npx-1
         sendcounts2(p) = sendcounts2(p)*ny1
      end do
      sdispls2(0) = 0
      do p = 1, npx-1
         sdispls2(p) = sdispls2(p-1) + sendcounts2(p-1)
      end do

      allocate(recvcounts2(0:nprocs-1))
      allocate(rdispls2(0:nprocs-1))
      call MPI_Alltoall(sendcounts2,1,MPI_INTEGER, &
                        recvcounts2,1,MPI_INTEGER,MPI_X_WORLD,ierr)
      rdispls2(0) = 0
      do p = 1, npx-1
         rdispls2(p) = rdispls2(p-1) + recvcounts2(p-1)
      end do
#endif

!$omp parallel
#ifndef MPI
!$omp single
      allocate(buf0(bigNY,bigNX))
      allocate(buf1(ny,   bigNX))
      allocate(buf2(bigNX,ny))
!$omp end single

!$omp do private(j)
      do i = 1, bigNX
         do j = 1, bigNY
! === Elastic loading with interpolation =======================================
!           buf0(j,i) = cg%zz(i+zeroIX-1,j+zeroIY-1)
            buf0(j,i) = cg%loading%delta(i+zeroIX-1,j+zeroIY-1)
! ==============================================================================
         end do
      end do

      ! Y interpolation
      if(use_linear == 1) then
!$omp do
         do i = 1, bigNX
            call linear_inerpolation(fg%my%nr, bigNY, buf0(1,i), ny, buf1(1,i))
         end do
      else
!$omp do
         do i = 1, bigNX
            call spline_inerpolation(fg%my%nr, bigNY, buf0(1,i), ny, buf1(1,i))
         end do
      end if
!$omp do private(i)
      do j = 1, ny
         do i = 1, bigNX
            buf2(i,j) = buf1(j,i)
         end do
      end do
      ! X interpolation
      if(use_linear == 1) then
!$omp do
         do j = 1, ny
! === Elastic loading with interpolation =======================================
!           call linear_inerpolation(fg%my%nr, bigNX, buf2(1,j), nx, fg%zz(1,j))
            call linear_inerpolation(fg%my%nr, bigNX, buf2(1,j), nx, fg%loading%delta(1,j))
! ==============================================================================
         end do
      else
!$omp do
         do j = 1, ny
! === Elastic loading with interpolation =======================================
!           call spline_inerpolation(fg%my%nr, bigNX, buf2(1,j), nx, fg%zz(1,j))
            call spline_inerpolation(fg%my%nr, bigNX, buf2(1,j), nx, fg%loading%delta(1,j))
! ==============================================================================
         end do
      end if

! === Elastic loading with interpolation =======================================
!$omp do private(i)
      do j = 1, ny
         do i = 1, nx
            fg%wave_field%hz(i,j) = fg%wave_field%hz(i,j) + fg%loading%delta(i,j)
         end do
      end do
! ==============================================================================
!$omp single
      deallocate(buf0)
      deallocate(buf1)
      deallocate(buf2)
!$omp end single
#else
      ! Collect all coarse elements to rank 0.
!$omp single
      allocate(tmp2d(bigNY,bigNX))
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
      allocate(commbuf(bigNY,bigNX))
! ==============================================================================
!$omp end single
!$omp do private(j)
      do i = 1, bigNX
         do j = 1, bigNY
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
!           tmp2d(j,i) = 0.0d0
            commbuf(j,i) = 0.0d0
! ==============================================================================
         end do
      end do
!$omp do private(ic, j, jc)
      do i = 1, bigNX
         ic = i + zeroIX - 1
         if((ic >= cg%my%ix) .and. (ic <= cg%my%ixend)) then
            do j = 1, bigNY
               jc = j + zeroIY - 1
               if((jc >= cg%my%iy) .and. (jc <= cg%my%iyend)) then
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
!                 tmp2d(j,i) = cg%zz(ic-cg%my%kx+1,jc-cg%my%ky+1)
! === Elastic loading with interpolation =======================================
!                 commbuf(j,i) = cg%zz(ic-cg%my%kx+1,jc-cg%my%ky+1)
                  commbuf(j,i) = cg%loading%delta(ic-cg%my%kx+1,jc-cg%my%ky+1)
! ==============================================================================
! ==============================================================================
               end if
            end do
         end if
      end do
!$omp single
! === MPI_IN_PLACE cannot be utilized with MPI_Reduce on SX. ===================
!      if(myrank == 0) then
!         call MPI_Reduce(MPI_IN_PLACE, tmp2d, bigNY*bigNX, REAL_MPI, &
!#ifndef MULTI
!                         MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!#else
!                         MPI_SUM, 0, MPI_MEMBER_WORLD, ierr)
!#endif
!      else
!         call MPI_Reduce(tmp2d, tmp2d, bigNY*bigNX, REAL_MPI, &
!#ifndef MULTI
!                         MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!#else
!                         MPI_SUM, 0, MPI_MEMBER_WORLD, ierr)
!#endif
!      end if
      call MPI_Reduce(commbuf, tmp2d, bigNY*bigNX, REAL_MPI, &
#ifndef MULTI
                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
#else
                      MPI_SUM, 0, MPI_MEMBER_WORLD, ierr)
#endif
      deallocate(commbuf)
! ==============================================================================

      ! Scatter tmp2d(bigNY,bigNX) -> buf0(bigNY,nx0)
      allocate(buf0(bigNY,nx0))
      call MPI_Scatterv(tmp2d, sendcounts0, sdispls0, REAL_MPI, &
                        buf0,   bigNY*nx0,             REAL_MPI, &
#ifndef MULTI
                        0, MPI_COMM_WORLD, ierr)
#else
                        0, MPI_MEMBER_WORLD, ierr)
#endif
      deallocate(tmp2d)

      ! Y interpolation buf0(bigNY,nx0) -> buf1(totalNy,nx0)
      allocate(buf1(totalNy,nx0))
!$omp end single
      if(use_linear == 1) then
!$omp do
         do i = 1, nx0
            call linear_inerpolation(fg%my%nr, bigNY, buf0(1,i), totalNy, buf1(1,i))
         end do
      else
!$omp do
         do i = 1, nx0
            call spline_inerpolation(fg%my%nr, bigNY, buf0(1,i), totalNy, buf1(1,i))
         end do
      end if
!$omp single
      deallocate(buf0)

      ! Local transposition buf1(totalNy,nx0) -> tmp2d(nx0,totalNy)
      allocate(tmp2d(nx0,totalNy))
!$omp end single
!$omp do private(i)
      do j = 1, totalNy
         do i = 1, nx0
            tmp2d(i,j) = buf1(j,i)
         end do
      end do
!$omp single
      deallocate(buf1)
      ! Transposition tmp2d(nx0,totalNy) -> tmp1d(bigNX*ny1)
      allocate(tmp1d(bigNX*ny1))
      call MPI_Alltoallv(tmp2d, sendcounts1, sdispls1, REAL_MPI, &
                         tmp1d, recvcounts1, rdispls1, REAL_MPI, &
#ifndef MULTI
                         MPI_COMM_WORLD, ierr)
#else
                         MPI_MEMBER_WORLD, ierr)
#endif
      deallocate(tmp2d)
      ! Rearrange tmp1d(bigNX*ny1) -> buf2(bigNX,ny1)
      allocate(buf2(bigNX,ny1))
!$omp end single
      do p = 0, nprocs-1
         ist = rdispls1(p)/ny1 + 1
         ien = ist + recvcounts1(p)/ny1 - 1
!$omp do private(i, ind)
         do j = 1, ny1
            do i = ist, ien
               ind = rdispls1(p) + (ien - ist + 1)*(j - 1) + i - ist + 1
               buf2(i,j) = tmp1d(ind)
            end do
         end do
      end do
!$omp single
      deallocate(tmp1d)

      ! X interpolation buf2(i,j) -> buf3(totalNx,ny1)
      allocate(buf3(totalNx,ny1))
!$omp end single
      if(use_linear == 1) then
!$omp do
         do j = 1, ny1
            call linear_inerpolation(fg%my%nr, bigNX, buf2(1,j), totalNx, buf3(1,j))
         end do
      else
!$omp do
         do j = 1, ny1
            call spline_inerpolation(fg%my%nr, bigNX, buf2(1,j), totalNx, buf3(1,j))
         end do
      end if
!$omp single
      deallocate(buf2)

      ! Local transposition buf3(totalNx,ny1) -> tmp2d(ny1,totalNx)
      allocate(tmp2d(ny1,totalNx))
!$omp end single
!$omp do private(j)
      do i = 1, totalNx
         do j = 1, ny1
            tmp2d(j,i) = buf3(i,j)
         end do
      end do
!$omp single
      deallocate(buf3)
      ! Transposition tmp2d(ny1,totalNx) -> tmp1d(nby*nbx)
      allocate(tmp1d(nby*nbx))
      call MPI_Alltoallv(tmp2d, sendcounts2, sdispls2, REAL_MPI, &
                         tmp1d, recvcounts2, rdispls2, REAL_MPI, &
                         MPI_X_WORLD, ierr)
      deallocate(tmp2d)
!$omp end single

      ! Finally, rearrange tmp1d(nby*nbx) -> fg%zz(nx,ny)
      ib = ix - kx
      jb = iy - ky
      do p = 0, npx-1
         jst = rdispls2(p)/nbx + 1
         jen = jst + recvcounts2(p)/nbx - 1
!$omp do private(j, ind)
         do i = 1, nbx
            do j = jst, jen
               ind = rdispls2(p) + (jen - jst + 1)*(i - 1) + j - jst + 1
! === Elastic loading with interpolation =======================================
!              fg%zz(i+ib,j+jb) = tmp1d(ind)
               fg%loading%delta(i+ib,j+jb) = tmp1d(ind)
! ==============================================================================
            end do
         end do
      end do
! === Elastic loading with interpolation =======================================
!$omp do private(i)
      do j = 1, fg%my%ny
         do i = 1, fg%my%nx
            fg%wave_field%hz(i,j) = fg%wave_field%hz(i,j) + fg%loading%delta(i,j)
         end do
      end do
! ==============================================================================
!$omp single
      deallocate(tmp1d)

      deallocate(sendcounts0)
      deallocate(sdispls0)
      deallocate(sendcounts1)
      deallocate(sdispls1)
      deallocate(recvcounts1)
      deallocate(rdispls1)
      deallocate(sendcounts2)
      deallocate(sdispls2)
      deallocate(recvcounts2)
      deallocate(rdispls2)
!$omp end single
#endif
!$omp end parallel
      return
! === Elastic loading with interpolation =======================================
!  end subroutine interp2fine_init_disp
   end subroutine interp2fine_elastic_loading
#endif
! ==============================================================================

end module mod_interpolation
