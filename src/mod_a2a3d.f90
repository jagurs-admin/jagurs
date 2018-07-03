#define DEBUG
#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_a2a3d
use mpi
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
implicit none

integer(kind=4), private :: nprocs, myrank, max_desc, num_desc
integer(kind=4), dimension(3), private :: nprocs3d, myrank3d, MPI_XYZ_GROUP ! 1:X, 2:Y, 3:Z

type A2A3D_descriptor
   integer(kind=4) :: tot_sendx, tot_sendy, tot_sendz
   integer(kind=4), allocatable, dimension(:) :: scntx, sdispx, rcntx, rdispx
   integer(kind=4), allocatable, dimension(:) :: scnty, sdispy, rcnty, rdispy
   integer(kind=4), allocatable, dimension(:) :: scntz, sdispz, rcntz, rdispz
   real(kind=REAL_BYTE), allocatable, dimension(:) :: sbufx, rbufx
   real(kind=REAL_BYTE), allocatable, dimension(:) :: sbufy, rbufy
   real(kind=REAL_BYTE), allocatable, dimension(:) :: sbufz
   integer(kind=4), allocatable, dimension(:) :: map_sbuf2x
   integer(kind=4), allocatable, dimension(:) :: map_x2y
   integer(kind=4), allocatable, dimension(:) :: map_y2z
end type A2A3D_descriptor

type(A2A3D_descriptor), allocatable, dimension(:) :: desc

contains

   subroutine A2A3D_init(nprocs_given, myrank_given, nproc3d_given, max_desc_given)
      integer(kind=4), intent(in) :: nprocs_given, myrank_given, max_desc_given
      integer(kind=4), dimension(3), intent(in) :: nproc3d_given
      integer(kind=4), dimension(3) :: color, nprcs3d_tmp
      integer(kind=4) :: key, i, ierr

      nprocs = nprocs_given
      myrank = myrank_given
      nprocs3d = nproc3d_given
      max_desc = max_desc_given

      if(product(nprocs3d) /= nprocs) then
         if(myrank == 0) then
            write(0,'(a)') '###########################################################'
            write(0,'(a)') '# A2A3D ERROR!!!                                          #'
            write(0,'(a)') '# Product of X/Y/Z procs. is NOT equal to num. of procs.! #'
            write(0,'(a,i6,a)') '# Num. of procs.: ', nprocs, &
                                '                                  #'
            write(0,'(a,3i6,a)') '# Num. of X/Y/Z procs.: ', nprocs3d, &
                                 '                #'
            write(0,'(a)') '###########################################################'
         end if
#ifndef MULTI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         call MPI_Abort(MPI_COMM_WORLD,999,ierr)
#else
         call MPI_Barrier(MPI_MEMBER_WORLD,ierr)
         call MPI_Abort(MPI_MEMBER_WORLD,999,ierr)
#endif
         stop
      end if

      ! Incriment X -> Y -> Z
      color(1) = myrank/nprocs3d(1)
      color(2) = myrank/(nprocs3d(1)*nprocs3d(2))
      color(2) = color(2)*nprocs3d(1) + mod(myrank,nprocs3d(1))
      color(3) = mod(myrank,nprocs3d(1)*nprocs3d(2))
      key = myrank

      do i = 1, 3
#ifndef MULTI
         call MPI_comm_split(MPI_COMM_WORLD, color(i), key, MPI_XYZ_GROUP(i), ierr)
#else
         call MPI_comm_split(MPI_MEMBER_WORLD, color(i), key, MPI_XYZ_GROUP(i), ierr)
#endif
         call MPI_comm_size(MPI_XYZ_GROUP(i), nprcs3d_tmp(i), ierr)
         call MPI_comm_rank(MPI_XYZ_GROUP(i), myrank3d(i), ierr)
      end do

      allocate(desc(max_desc))
      num_desc = 0

      return
   end subroutine A2A3D_init

   subroutine A2A3D_plan(sendcnts, handler)
      integer(kind=4), dimension(0:nprocs-1), intent(in) :: sendcnts
      integer(kind=4), intent(out) :: handler

      integer(kind=4), allocatable, dimension(:) :: sendcnts_yzx
      integer(kind=4), allocatable, dimension(:) :: sendcnts_zy
      integer(kind=4), allocatable, dimension(:) :: map_x2y_count
      integer(kind=4), allocatable, dimension(:) :: map_y2z_count

      integer(kind=4) :: i, j, k, l, n, p, send_count, start_ind, ind, ierr

      num_desc = num_desc + 1
      handler = num_desc
      n = num_desc
! ========================================================================================
! === X-direction ========================================================================
! ========================================================================================
! ----------------------------------------------------------------------------------------
! --- Sort send buffer -------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
      desc(n)%tot_sendx = sum(sendcnts)

      allocate(desc(n)%map_sbuf2x(desc(n)%tot_sendx))
      allocate(desc(n)%sbufx(desc(n)%tot_sendx))

      ind = 0
      do i = 0, nprocs3d(1)-1
         do k = 0, nprocs3d(3)-1
            do j = 0, nprocs3d(2)-1
               p = k*nprocs3d(1)*nprocs3d(2)+ j*nprocs3d(1) + i
               start_ind = 0
               do l = 0, p-1
                  start_ind = start_ind + sendcnts(l)
               end do
               do l = 1, sendcnts(p)
                  ind = ind + 1
                  desc(n)%map_sbuf2x(ind) = start_ind + l
               end do
            end do
         end do
      end do
! ----------------------------------------------------------------------------------------
! --- Alltoallv parameters for X-direction -----------------------------------------------
! ----------------------------------------------------------------------------------------
      allocate(desc(n)%scntx(0:nprocs3d(1)-1))
      allocate(desc(n)%rcntx(0:nprocs3d(1)-1))
      allocate(desc(n)%sdispx(0:nprocs3d(1)-1))
      allocate(desc(n)%rdispx(0:nprocs3d(1)-1))

      do i = 0, nprocs3d(1)-1
         send_count = 0
         do k = 0, nprocs3d(3)-1
            do j = 0, nprocs3d(2)-1
               p = k*nprocs3d(1)*nprocs3d(2)+ j*nprocs3d(1) + i
               send_count = send_count + sendcnts(p)
            end do
         end do
         desc(n)%scntx(i) = send_count
      end do

      call MPI_Alltoall(desc(n)%scntx, 1, MPI_INTEGER, &
                        desc(n)%rcntx, 1, MPI_INTEGER, MPI_XYZ_GROUP(1), ierr)

      desc(n)%sdispx(0) = 0
      desc(n)%rdispx(0) = 0
      do i = 1, nprocs3d(1)-1
         desc(n)%sdispx(i) = desc(n)%sdispx(i-1) + desc(n)%scntx(i-1)
         desc(n)%rdispx(i) = desc(n)%rdispx(i-1) + desc(n)%rcntx(i-1)
      end do

      desc(n)%tot_sendy = sum(desc(n)%rcntx)

      allocate(desc(n)%rbufx(desc(n)%tot_sendy))
      allocate(desc(n)%sbufy(desc(n)%tot_sendy))
! ----------------------------------------------------------------------------------------
! --- Sort recv buffer for Y-direction ---------------------------------------------------
! ----------------------------------------------------------------------------------------
      allocate(desc(n)%map_x2y(desc(n)%tot_sendy))

      allocate(sendcnts_yzx(0:nprocs-1))
      ind = 0
      do i = 0, nprocs3d(1)-1
         do k = 0, nprocs3d(3)-1
            do j = 0, nprocs3d(2)-1
               p = k*nprocs3d(1)*nprocs3d(2)+ j*nprocs3d(1) + i
               sendcnts_yzx(ind) = sendcnts(p)
               ind = ind + 1
            end do
         end do
      end do

      allocate(map_x2y_count(0:nprocs-1))
      call MPI_Alltoall(sendcnts_yzx,  nprocs3d(2)*nprocs3d(3), MPI_INTEGER, &
                        map_x2y_count, nprocs3d(2)*nprocs3d(3), MPI_INTEGER, MPI_XYZ_GROUP(1), ierr)

      ind = 0
      do j = 0, nprocs3d(2)-1
         do k = 0, nprocs3d(3)-1
            do i = 0, nprocs3d(1)-1
               p = i*nprocs3d(2)*nprocs3d(3)+ k*nprocs3d(2) + j
               start_ind = 0
               do l = 0, p-1
                  start_ind = start_ind + map_x2y_count(l)
               end do
               do l = 1, map_x2y_count(p)
                  ind = ind + 1
                  desc(n)%map_x2y(ind) = start_ind + l
               end do
            end do
         end do
      end do

      deallocate(sendcnts_yzx)
! ========================================================================================
! === Y-direction ========================================================================
! ========================================================================================
! ----------------------------------------------------------------------------------------
! --- Alltoallv parameters for Y-direction -----------------------------------------------
! ----------------------------------------------------------------------------------------
      allocate(desc(n)%scnty(0:nprocs3d(2)-1))
      allocate(desc(n)%rcnty(0:nprocs3d(2)-1))
      allocate(desc(n)%sdispy(0:nprocs3d(2)-1))
      allocate(desc(n)%rdispy(0:nprocs3d(2)-1))

      allocate(sendcnts_zy(0:nprocs3d(2)*nprocs3d(3)-1))

      sendcnts_zy = 0
      do k = 0, nprocs3d(3)-1
         do j = 0, nprocs3d(2)-1
            do i = 0, nprocs3d(1)-1
               p = i*nprocs3d(2)*nprocs3d(3)+ k*nprocs3d(2) + j
               sendcnts_zy(j*nprocs3d(3)+k) = sendcnts_zy(j*nprocs3d(3)+k) + map_x2y_count(p)
            end do
         end do
      end do

      do j = 0, nprocs3d(2)-1
         send_count = 0
         do k = 0, nprocs3d(3)-1
            p = j*nprocs3d(3) + k
            send_count = send_count + sendcnts_zy(p)
         end do
         desc(n)%scnty(j) = send_count
      end do

      call MPI_Alltoall(desc(n)%scnty, 1, MPI_INTEGER, &
                        desc(n)%rcnty, 1, MPI_INTEGER, MPI_XYZ_GROUP(2), ierr)

      desc(n)%sdispy(0) = 0
      desc(n)%rdispy(0) = 0
      do j = 1, nprocs3d(2)-1
         desc(n)%sdispy(j) = desc(n)%sdispy(j-1) + desc(n)%scnty(j-1)
         desc(n)%rdispy(j) = desc(n)%rdispy(j-1) + desc(n)%rcnty(j-1)
      end do

      desc(n)%tot_sendz = sum(desc(n)%rcnty)

      allocate(desc(n)%rbufy(desc(n)%tot_sendz))
      allocate(desc(n)%sbufz(desc(n)%tot_sendz))
 
      deallocate(map_x2y_count)
! ----------------------------------------------------------------------------------------
! --- Sort recv buffer for Z-direction ---------------------------------------------------
! ----------------------------------------------------------------------------------------
      allocate(desc(n)%map_y2z(desc(n)%tot_sendz))
      allocate(map_y2z_count(0:nprocs3d(2)*nprocs3d(3)-1))

      call MPI_Alltoall(sendcnts_zy,   nprocs3d(3), MPI_INTEGER, &
                        map_y2z_count, nprocs3d(3), MPI_INTEGER, MPI_XYZ_GROUP(2), ierr)

      ind = 0
      do k = 0, nprocs3d(3)-1
         do j = 0, nprocs3d(2)-1
            p = j*nprocs3d(3) + k
            start_ind = 0
            do l = 0, p-1
               start_ind = start_ind + map_y2z_count(l)
            end do
            do l = 1, map_y2z_count(p)
               ind = ind + 1
               desc(n)%map_y2z(ind) = start_ind + l
            end do
         end do
      end do

      deallocate(sendcnts_zy)
! ========================================================================================
! === Z-direction ========================================================================
! ========================================================================================
! ----------------------------------------------------------------------------------------
! --- Alltoallv parameters for Z-direction -----------------------------------------------
! ----------------------------------------------------------------------------------------
      allocate(desc(n)%scntz(0:nprocs3d(3)-1))
      allocate(desc(n)%rcntz(0:nprocs3d(3)-1))
      allocate(desc(n)%sdispz(0:nprocs3d(3)-1))
      allocate(desc(n)%rdispz(0:nprocs3d(3)-1))

      do k = 0, nprocs3d(3)-1
         send_count = 0
         do j = 0, nprocs3d(2)-1
            p = j*nprocs3d(3) + k
            send_count = send_count + map_y2z_count(p)
         end do
         desc(n)%scntz(k) = send_count
      end do

      call MPI_Alltoall(desc(n)%scntz, 1, MPI_INTEGER, &
                        desc(n)%rcntz, 1, MPI_INTEGER, MPI_XYZ_GROUP(3), ierr)

      desc(n)%sdispz(0) = 0
      desc(n)%rdispz(0) = 0
      do k = 1, nprocs3d(3)-1
         desc(n)%sdispz(k) = desc(n)%sdispz(k-1) + desc(n)%scntz(k-1)
         desc(n)%rdispz(k) = desc(n)%rdispz(k-1) + desc(n)%rcntz(k-1)
      end do

      deallocate(map_y2z_count)
! ========================================================================================
! === Done!!! ============================================================================
! ========================================================================================
      return
   end subroutine A2A3D_plan

   subroutine A2A3D_finalize()
      integer(kind=4) :: n

      do n = 1, num_desc
! ========================================================================================
! === X-direction ========================================================================
! ========================================================================================
         deallocate(desc(n)%map_sbuf2x)
         deallocate(desc(n)%sbufx)
         deallocate(desc(n)%scntx)
         deallocate(desc(n)%rcntx)
         deallocate(desc(n)%sdispx)
         deallocate(desc(n)%rdispx)
         deallocate(desc(n)%rbufx)
         deallocate(desc(n)%sbufy)
         deallocate(desc(n)%map_x2y)
! ========================================================================================
! === Y-direction ========================================================================
! ========================================================================================
         deallocate(desc(n)%scnty)
         deallocate(desc(n)%rcnty)
         deallocate(desc(n)%sdispy)
         deallocate(desc(n)%rdispy)
         deallocate(desc(n)%rbufy)
         deallocate(desc(n)%sbufz)
         deallocate(desc(n)%map_y2z)
! ========================================================================================
! === Z-direction ========================================================================
! ========================================================================================
         deallocate(desc(n)%scntz)
         deallocate(desc(n)%rcntz)
         deallocate(desc(n)%sdispz)
         deallocate(desc(n)%rdispz)
! ========================================================================================
! === Done!!! ============================================================================
! ========================================================================================
      end do

      deallocate(desc)
      num_desc = 0

      return
   end subroutine A2A3D_finalize

   subroutine A2A3D_execute(ssize, sbuf, rsize, rbuf, n)
      integer(kind=4), intent(in) :: ssize, rsize, n
      real(kind=REAL_BYTE), intent(in),  dimension(ssize) :: sbuf
      real(kind=REAL_BYTE), intent(out), dimension(rsize) :: rbuf
      integer(kind=4) :: i, ierr

      do i = 1, desc(n)%tot_sendx
         desc(n)%sbufx(i) = sbuf(desc(n)%map_sbuf2x(i))
      end do

      call MPI_Alltoallv(desc(n)%sbufx, desc(n)%scntx, desc(n)%sdispx, REAL_MPI, &
                         desc(n)%rbufx, desc(n)%rcntx, desc(n)%rdispx, REAL_MPI, &
                         MPI_XYZ_GROUP(1), ierr)

      do i = 1, desc(n)%tot_sendy
         desc(n)%sbufy(i) = desc(n)%rbufx(desc(n)%map_x2y(i))
      end do

      call MPI_Alltoallv(desc(n)%sbufy, desc(n)%scnty, desc(n)%sdispy, REAL_MPI, &
                         desc(n)%rbufy, desc(n)%rcnty, desc(n)%rdispy, REAL_MPI, &
                         MPI_XYZ_GROUP(2), ierr)

      do i = 1, desc(n)%tot_sendz
         desc(n)%sbufz(i) = desc(n)%rbufy(desc(n)%map_y2z(i))
      end do

      call MPI_Alltoallv(desc(n)%sbufz, desc(n)%scntz, desc(n)%sdispz, REAL_MPI, &
                         rbuf,          desc(n)%rcntz, desc(n)%rdispz, REAL_MPI, &
                         MPI_XYZ_GROUP(3), ierr)

      return
   end subroutine A2A3D_execute

   subroutine decomp_default(n, a2a3d_decomp)
      implicit none
      integer(kind=4), intent(in) :: n
      integer(kind=4), dimension(3), intent(out) :: a2a3d_decomp
      integer(kind=4) :: sqrtn, numdiv, tmpn, i
      integer(kind=4), allocatable, dimension(:) :: divlist
      integer(kind=4), dimension(3) :: tmpa

      sqrtn = sqrt(real(n))
      allocate(divlist(sqrtn+1))

      ! Make the list of divisors 'divlist'.
      tmpn = n
      numdiv=0
      do
         sqrtn = sqrt(real(tmpn))
         do i = 2, sqrtn
            if(mod(tmpn,i) == 0) then
               tmpn = tmpn/i
               numdiv = numdiv + 1
               divlist(numdiv) = i
               exit
            end if
         end do
         if(i == sqrtn+1) then
            numdiv = numdiv + 1
            divlist(numdiv) = tmpn
            exit
         end if
      end do

      ! Make 3 factors 'tmpa' as far as close.
      tmpa = 1
      do i = numdiv, 1, -1
         if(tmpa(1) <= tmpa(2)) then
            if(tmpa(1) <= tmpa(3)) then
               tmpa(1) = tmpa(1)*divlist(i)
            else
               tmpa(3) = tmpa(3)*divlist(i)
            end if
         else
            if(tmpa(2) <= tmpa(3)) then
               tmpa(2) = tmpa(2)*divlist(i)
            else
               tmpa(3) = tmpa(3)*divlist(i)
            end if
         end if
      end do

      ! Sort 'tmpa' and make 'a2a3d_decomp'.
      a2a3d_decomp(1) = maxval(tmpa)
      a2a3d_decomp(3) = minval(tmpa)
      a2a3d_decomp(2) = n/a2a3d_decomp(1)/a2a3d_decomp(3)

      ! Check 'a2a3d_decomp'.
      if(product(a2a3d_decomp) /= n) then
         write(0,'(a)') 'ERROR! Invalid decomposition!'
         write(0,'(a,i6)') 'n = ', n
         write(0,'(a,3i6)') 'a2a3d_decomp = ', a2a3d_decomp
      end if

      deallocate(divlist)

      return
   end subroutine decomp_default

end module mod_a2a3d
