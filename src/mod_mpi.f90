#include "real.h"
module mod_mpi
use mpi
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
use mod_grid
use mod_params
implicit none

integer(kind=4) :: north_rank, south_rank, east_rank, west_rank

contains

   subroutine allocate_edges(grid)
      type(data_grids), target, intent(inout) :: grid
      integer(kind=4) :: nx, ny
      integer(kind=4) :: px, py
      integer(kind=4) :: rx, ry
      integer(kind=4) :: has_boundary

      nx = grid%my%nx
      ny = grid%my%ny
      px = grid%my%px
      py = grid%my%py
      rx = grid%my%rx
      ry = grid%my%ry
      has_boundary = grid%my%has_boundary

      north_rank = (ry-1)*px+rx
      south_rank = (ry+1)*px+rx
      east_rank = ry*px+rx+1
      west_rank = ry*px+rx-1

      if(iand(has_boundary, NORTH_BOUND) /= 0) north_rank = MPI_PROC_NULL
      if(iand(has_boundary, SOUTH_BOUND) /= 0) south_rank = MPI_PROC_NULL
      if(iand(has_boundary, EAST_BOUND)  /= 0) east_rank  = MPI_PROC_NULL
      if(iand(has_boundary, WEST_BOUND)  /= 0) west_rank  = MPI_PROC_NULL

      ! for exchange_edges
      allocate(grid%edges%fne(nx,4))
      allocate(grid%edges%fse(nx,5))
      allocate(grid%edges%fnb(nx,5))
      allocate(grid%edges%fsb(nx,4))

      allocate(grid%edges%hne(nx,3))
      allocate(grid%edges%hse(nx,2))
      allocate(grid%edges%hnb(nx,2))
      allocate(grid%edges%hsb(nx,3))

      allocate(grid%edges%fee(-1:ny+1,5))
      allocate(grid%edges%fwe(-1:ny+1,4))
      allocate(grid%edges%feb(-1:ny+1,4))
      allocate(grid%edges%fwb(-1:ny+1,5))

      allocate(grid%edges%hee(0:ny+2,2))
      allocate(grid%edges%hwe(0:ny+2,3))
      allocate(grid%edges%heb(0:ny+2,3))
      allocate(grid%edges%hwb(0:ny+2,2))

      ! for exchange_edges_dz
      allocate(grid%edges%dne(nx,3))
      allocate(grid%edges%dse(nx,2))
      allocate(grid%edges%dnb(nx,2))
      allocate(grid%edges%dsb(nx,3))

      allocate(grid%edges%dee(0:ny+2,2))
      allocate(grid%edges%dwe(0:ny+2,3))
      allocate(grid%edges%deb(0:ny+2,3))
      allocate(grid%edges%dwb(0:ny+2,2))

      ! exchange_edges_disp_f[xy]
      allocate(grid%edges%fne_d(nx,1))
      allocate(grid%edges%fse_d(nx,1))
      allocate(grid%edges%fnb_d(nx,1))
      allocate(grid%edges%fsb_d(nx,1))

      allocate(grid%edges%fee_d(ny,1))
      allocate(grid%edges%fwe_d(ny,1))
      allocate(grid%edges%feb_d(ny,1))
      allocate(grid%edges%fwb_d(ny,1))

      ! exchange_edges_wod
      allocate(grid%edges%wne_d(nx,1))
      allocate(grid%edges%wse_d(nx,1))
      allocate(grid%edges%wnb_d(nx,1))
      allocate(grid%edges%wsb_d(nx,1))

      allocate(grid%edges%wee_d(ny,1))
      allocate(grid%edges%wwe_d(ny,1))
      allocate(grid%edges%web_d(ny,1))
      allocate(grid%edges%wwb_d(ny,1))

      return
   end subroutine allocate_edges

   subroutine deallocate_edges(grid)
      type(data_grids), target, intent(inout) :: grid

      ! for exchange_edges
      deallocate(grid%edges%fne)
      deallocate(grid%edges%fse)
      deallocate(grid%edges%fnb)
      deallocate(grid%edges%fsb)

      deallocate(grid%edges%hne)
      deallocate(grid%edges%hse)
      deallocate(grid%edges%hnb)
      deallocate(grid%edges%hsb)

      deallocate(grid%edges%fee)
      deallocate(grid%edges%fwe)
      deallocate(grid%edges%feb)
      deallocate(grid%edges%fwb)

      deallocate(grid%edges%hee)
      deallocate(grid%edges%hwe)
      deallocate(grid%edges%heb)
      deallocate(grid%edges%hwb)

      ! for exchange_edges_dz
      deallocate(grid%edges%dne)
      deallocate(grid%edges%dse)
      deallocate(grid%edges%dnb)
      deallocate(grid%edges%dsb)

      deallocate(grid%edges%dee)
      deallocate(grid%edges%dwe)
      deallocate(grid%edges%deb)
      deallocate(grid%edges%dwb)

      ! exchange_edges_disp_f[xy]
      deallocate(grid%edges%fne_d)
      deallocate(grid%edges%fse_d)
      deallocate(grid%edges%fnb_d)
      deallocate(grid%edges%fsb_d)

      deallocate(grid%edges%fee_d)
      deallocate(grid%edges%fwe_d)
      deallocate(grid%edges%feb_d)
      deallocate(grid%edges%fwb_d)

      ! exchange_edges_wod
      deallocate(grid%edges%wne_d)
      deallocate(grid%edges%wse_d)
      deallocate(grid%edges%wnb_d)
      deallocate(grid%edges%wsb_d)

      deallocate(grid%edges%wee_d)
      deallocate(grid%edges%wwe_d)
      deallocate(grid%edges%web_d)
      deallocate(grid%edges%wwb_d)

      return
   end subroutine deallocate_edges

   subroutine exchange_edges(mode, grid)
      integer(kind=4), intent(in) :: mode
      type(data_grids), target, intent(inout) :: grid

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_buf
      integer(kind=4) :: nx
      integer(kind=4) :: ny
      integer(kind=4) :: has_boundary
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fy
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz

      integer(kind=4), dimension(2) :: sreqs, rreqs
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: sstat, rstat
      integer(kind=4) :: ierr = 0

      nx = grid%my%nx
      ny = grid%my%ny
      has_boundary = grid%my%has_boundary

      fx => grid%wave_field%fx
      fy => grid%wave_field%fy
      hz => grid%wave_field%hz

      !**************************************
      !*                                    *
      !* North-South Communication          *
      !*                                    *
      !**************************************
      if(mode == VEL) then
         north_edge => grid%edges%fne
         south_edge => grid%edges%fse
         north_buf  => grid%edges%fnb
         south_buf  => grid%edges%fsb

         if(iand(has_boundary, NORTH_BOUND) == 0) north_edge(1:nx,1) = fx(1:nx,2)
         if(iand(has_boundary, NORTH_BOUND) == 0) north_edge(1:nx,2) = fy(1:nx,2)
         if(iand(has_boundary, NORTH_BOUND) == 0) north_edge(1:nx,3) = fx(1:nx,3)
         if(iand(has_boundary, NORTH_BOUND) == 0) north_edge(1:nx,4) = fy(1:nx,3)

         if(iand(has_boundary, SOUTH_BOUND) == 0) south_edge(1:nx,1) = fx(1:nx,ny-1)
         if(iand(has_boundary, SOUTH_BOUND) == 0) south_edge(1:nx,4) = fx(1:nx,ny-2)
         if(iand(has_boundary, SOUTH_BOUND) == 0) then
            south_edge(1:nx,2) = fy(1:nx,ny-1)
            south_edge(1:nx,3) = fy(1:nx,ny-2)
            south_edge(1:nx,5) = fy(1:nx,ny-3)
         end if

#ifndef MULTI
         call MPI_Irecv(north_buf,  5*nx, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, rreqs(1), ierr)
         call MPI_Irecv(south_buf,  4*nx, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, rreqs(2), ierr)
         call MPI_Isend(north_edge, 4*nx, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, sreqs(1), ierr)
         call MPI_Isend(south_edge, 5*nx, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, sreqs(2), ierr)
#else
         call MPI_Irecv(north_buf,  5*nx, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, rreqs(1), ierr)
         call MPI_Irecv(south_buf,  4*nx, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, rreqs(2), ierr)
         call MPI_Isend(north_edge, 4*nx, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, sreqs(1), ierr)
         call MPI_Isend(south_edge, 5*nx, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
         call MPI_Waitall(2, sreqs, sstat, ierr)
         call MPI_Waitall(2, rreqs, rstat, ierr)

         if(iand(has_boundary, NORTH_BOUND) == 0) fx(1:nx,1) = north_buf(1:nx,1)
         if(iand(has_boundary, NORTH_BOUND) == 0) fx(1:nx,0) = north_buf(1:nx,4)
         if(iand(has_boundary, NORTH_BOUND) == 0) then
            fy(1:nx,1) = north_buf(1:nx,2)
            fy(1:nx,0) = north_buf(1:nx,3)
            fy(1:nx,-1) = north_buf(1:nx,5)
         end if
         if(iand(has_boundary, SOUTH_BOUND) == 0) fx(1:nx,ny) = south_buf(1:nx,1)
         if(iand(has_boundary, SOUTH_BOUND) == 0) fy(1:nx,ny) = south_buf(1:nx,2)
         if(iand(has_boundary, SOUTH_BOUND) == 0) fx(1:nx,ny+1) = south_buf(1:nx,3)
         if(iand(has_boundary, SOUTH_BOUND) == 0) fy(1:nx,ny+1) = south_buf(1:nx,4)
      else if(mode == HGT) then
         north_edge => grid%edges%hne
         south_edge => grid%edges%hse
         north_buf  => grid%edges%hnb
         south_buf  => grid%edges%hsb

         if(iand(has_boundary, NORTH_BOUND) == 0) then
            north_edge(1:nx,1) = hz(1:nx,2)
            north_edge(1:nx,2) = hz(1:nx,3)
            north_edge(1:nx,3) = hz(1:nx,4)
         end if
         if(iand(has_boundary, SOUTH_BOUND) == 0) then
            south_edge(1:nx,1) = hz(1:nx,ny-1)
            south_edge(1:nx,2) = hz(1:nx,ny-2)
         end if

#ifndef MULTI
         call MPI_Irecv(north_buf,  2*nx, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, rreqs(1), ierr)
         call MPI_Irecv(south_buf,  3*nx, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, rreqs(2), ierr)
         call MPI_Isend(north_edge, 3*nx, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, sreqs(1), ierr)
         call MPI_Isend(south_edge, 2*nx, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, sreqs(2), ierr)
#else
         call MPI_Irecv(north_buf,  2*nx, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, rreqs(1), ierr)
         call MPI_Irecv(south_buf,  3*nx, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, rreqs(2), ierr)
         call MPI_Isend(north_edge, 3*nx, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, sreqs(1), ierr)
         call MPI_Isend(south_edge, 2*nx, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
         call MPI_Waitall(2, sreqs, sstat, ierr)
         call MPI_Waitall(2, rreqs, rstat, ierr)

         if(iand(has_boundary, NORTH_BOUND) == 0) then
            hz(1:nx,1) = north_buf(1:nx,1)
            hz(1:nx,0) = north_buf(1:nx,2)
         end if
         if(iand(has_boundary, SOUTH_BOUND) == 0) then
            hz(1:nx,ny  ) = south_buf(1:nx,1)
            hz(1:nx,ny+1) = south_buf(1:nx,2)
            hz(1:nx,ny+2) = south_buf(1:nx,3)
         end if
      end if

      !**************************************
      !*                                    *
      !* East-West Communication            *
      !*                                    *
      !**************************************
      if(mode == VEL) then
         east_edge => grid%edges%fee
         west_edge => grid%edges%fwe
         east_buf  => grid%edges%feb
         west_buf  => grid%edges%fwb

         if(iand(has_boundary, EAST_BOUND) == 0) then
             east_edge(0:ny+1,1) = fx(nx-1,0:ny+1)
             east_edge(0:ny+1,2) = fx(nx-2,0:ny+1)
             east_edge(0:ny+1,4) = fx(nx-3,0:ny+1)
         end if
         if(iand(has_boundary, EAST_BOUND) == 0) east_edge(-1:ny+1,3) = fy(nx-1,-1:ny+1)
         if(iand(has_boundary, EAST_BOUND) == 0) east_edge(-1:ny+1,5) = fy(nx-2,-1:ny+1)
         if(iand(has_boundary, WEST_BOUND) == 0) west_edge( 0:ny+1,1) = fx(2, 0:ny+1)
         if(iand(has_boundary, WEST_BOUND) == 0) west_edge(-1:ny+1,2) = fy(2,-1:ny+1)
         if(iand(has_boundary, WEST_BOUND) == 0) west_edge( 0:ny+1,3) = fx(3, 0:ny+1)
         if(iand(has_boundary, WEST_BOUND) == 0) west_edge(-1:ny+1,4) = fy(3,-1:ny+1)

#ifndef MULTI
         call MPI_Irecv(east_buf,  4*(ny+3), REAL_MPI, east_rank, 2, MPI_COMM_WORLD, rreqs(1), ierr)
         call MPI_Irecv(west_buf,  5*(ny+3), REAL_MPI, west_rank, 3, MPI_COMM_WORLD, rreqs(2), ierr)
         call MPI_Isend(east_edge, 5*(ny+3), REAL_MPI, east_rank, 3, MPI_COMM_WORLD, sreqs(1), ierr)
         call MPI_Isend(west_edge, 4*(ny+3), REAL_MPI, west_rank, 2, MPI_COMM_WORLD, sreqs(2), ierr)
#else
         call MPI_Irecv(east_buf,  4*(ny+3), REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, rreqs(1), ierr)
         call MPI_Irecv(west_buf,  5*(ny+3), REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, rreqs(2), ierr)
         call MPI_Isend(east_edge, 5*(ny+3), REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, sreqs(1), ierr)
         call MPI_Isend(west_edge, 4*(ny+3), REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
         call MPI_Waitall(2, sreqs, sstat, ierr)
         call MPI_Waitall(2, rreqs, rstat, ierr)

         if(iand(has_boundary, EAST_BOUND) == 0) fx(nx,   0:ny+1) = east_buf( 0:ny+1,1)
         if(iand(has_boundary, EAST_BOUND) == 0) fy(nx,  -1:ny+1) = east_buf(-1:ny+1,2)
         if(iand(has_boundary, EAST_BOUND) == 0) fx(nx+1, 0:ny+1) = east_buf( 0:ny+1,3)
         if(iand(has_boundary, EAST_BOUND) == 0) fy(nx+1,-1:ny+1) = east_buf(-1:ny+1,4)
         if(iand(has_boundary, WEST_BOUND) == 0) then
             fx( 1,0:ny+1) = west_buf(0:ny+1,1)
             fx( 0,0:ny+1) = west_buf(0:ny+1,2)
             fx(-1,0:ny+1) = west_buf(0:ny+1,4)
         end if
         if(iand(has_boundary, WEST_BOUND) == 0) fy(1,-1:ny+1) = west_buf(-1:ny+1,3)
         if(iand(has_boundary, WEST_BOUND) == 0) fy(0,-1:ny+1) = west_buf(-1:ny+1,5)
      else if(mode == HGT) then
         east_edge => grid%edges%hee
         west_edge => grid%edges%hwe
         east_buf  => grid%edges%heb
         west_buf  => grid%edges%hwb

         if(iand(has_boundary, EAST_BOUND) == 0) then
            east_edge(0:ny+2,1) = hz(nx-1,0:ny+2)
            east_edge(0:ny+2,2) = hz(nx-2,0:ny+2)
         end if
         if(iand(has_boundary, WEST_BOUND) == 0) then
            west_edge(0:ny+2,1) = hz(2,0:ny+2)
            west_edge(0:ny+2,2) = hz(3,0:ny+2)
            west_edge(0:ny+2,3) = hz(4,0:ny+2)
         end if

#ifndef MULTI
         call MPI_Irecv(east_buf,  3*(ny+3), REAL_MPI, east_rank, 2, MPI_COMM_WORLD, rreqs(1), ierr)
         call MPI_Irecv(west_buf,  2*(ny+3), REAL_MPI, west_rank, 3, MPI_COMM_WORLD, rreqs(2), ierr)
         call MPI_Isend(east_edge, 2*(ny+3), REAL_MPI, east_rank, 3, MPI_COMM_WORLD, sreqs(1), ierr)
         call MPI_Isend(west_edge, 3*(ny+3), REAL_MPI, west_rank, 2, MPI_COMM_WORLD, sreqs(2), ierr)
#else
         call MPI_Irecv(east_buf,  3*(ny+3), REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, rreqs(1), ierr)
         call MPI_Irecv(west_buf,  2*(ny+3), REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, rreqs(2), ierr)
         call MPI_Isend(east_edge, 2*(ny+3), REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, sreqs(1), ierr)
         call MPI_Isend(west_edge, 3*(ny+3), REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
         call MPI_Waitall(2, sreqs, sstat, ierr)
         call MPI_Waitall(2, rreqs, rstat, ierr)

         if(iand(has_boundary, EAST_BOUND) == 0) then
            hz(nx,  0:ny+2) = east_buf(0:ny+2,1)
            hz(nx+1,0:ny+2) = east_buf(0:ny+2,2)
            hz(nx+2,0:ny+2) = east_buf(0:ny+2,3)
         end if
         if(iand(has_boundary, WEST_BOUND) == 0) then
            hz(1,0:ny+2) = west_buf(0:ny+2,1)
            hz(0,0:ny+2) = west_buf(0:ny+2,2)
         end if
      end if

      return
   end subroutine exchange_edges

   subroutine exchange_edges_dz(grid)
      type(data_grids), target, intent(inout) :: grid

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_buf
      integer(kind=4) :: nx
      integer(kind=4) :: ny
      integer(kind=4) :: has_boundary
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz

      integer(kind=4), dimension(2) :: sreqs, rreqs
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: sstat, rstat
      integer(kind=4) :: ierr = 0

      nx = grid%my%nx
      ny = grid%my%ny
      has_boundary = grid%my%has_boundary

      hz => grid%depth_field%dz

      !**************************************
      !*                                    *
      !* North-South Communication          *
      !*                                    *
      !**************************************
      north_edge => grid%edges%dne
      south_edge => grid%edges%dse
      north_buf  => grid%edges%dnb
      south_buf  => grid%edges%dsb

      if(iand(has_boundary, NORTH_BOUND) == 0) then
         north_edge(1:nx,1) = hz(1:nx,2)
         north_edge(1:nx,2) = hz(1:nx,3)
         north_edge(1:nx,3) = hz(1:nx,4)
      end if
      if(iand(has_boundary, SOUTH_BOUND) == 0) then
         south_edge(1:nx,1) = hz(1:nx,ny-1)
         south_edge(1:nx,2) = hz(1:nx,ny-2)
      end if

#ifndef MULTI
      call MPI_Irecv(north_buf,  2*nx, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  3*nx, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, 3*nx, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, 2*nx, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(north_buf,  2*nx, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  3*nx, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, 3*nx, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, 2*nx, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, NORTH_BOUND) == 0) then
         hz(1:nx,1) = north_buf(1:nx,1)
         hz(1:nx,0) = north_buf(1:nx,2)
      end if
      if(iand(has_boundary, SOUTH_BOUND) == 0) then
         hz(1:nx,ny  ) = south_buf(1:nx,1)
         hz(1:nx,ny+1) = south_buf(1:nx,2)
         hz(1:nx,ny+2) = south_buf(1:nx,3)
      end if

      !**************************************
      !*                                    *
      !* East-West Communication            *
      !*                                    *
      !**************************************
      east_edge => grid%edges%dee
      west_edge => grid%edges%dwe
      east_buf  => grid%edges%deb
      west_buf  => grid%edges%dwb

      if(iand(has_boundary, EAST_BOUND) == 0) then
         east_edge(0:ny+2,1) = hz(nx-1,0:ny+2)
         east_edge(0:ny+2,2) = hz(nx-2,0:ny+2)
      end if
      if(iand(has_boundary, WEST_BOUND) == 0) then
         west_edge(0:ny+2,1) = hz(2,0:ny+2)
         west_edge(0:ny+2,2) = hz(3,0:ny+2)
         west_edge(0:ny+2,3) = hz(4,0:ny+2)
      end if

#ifndef MULTI
      call MPI_Irecv(east_buf,  3*(ny+3), REAL_MPI, east_rank, 2, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  2*(ny+3), REAL_MPI, west_rank, 3, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, 2*(ny+3), REAL_MPI, east_rank, 3, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, 3*(ny+3), REAL_MPI, west_rank, 2, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(east_buf,  3*(ny+3), REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  2*(ny+3), REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, 2*(ny+3), REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, 3*(ny+3), REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, EAST_BOUND) == 0) then
         hz(nx,  0:ny+2) = east_buf(0:ny+2,1)
         hz(nx+1,0:ny+2) = east_buf(0:ny+2,2)
         hz(nx+2,0:ny+2) = east_buf(0:ny+2,3)
      end if
      if(iand(has_boundary, WEST_BOUND) == 0) then
         hz(1,0:ny+2) = west_buf(0:ny+2,1)
         hz(0,0:ny+2) = west_buf(0:ny+2,2)
      end if

      return
   end subroutine exchange_edges_dz
! === Displacement =============================================================
   subroutine exchange_edges_zz(grid)
      type(data_grids), target, intent(inout) :: grid

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_buf
      integer(kind=4) :: nx
      integer(kind=4) :: ny
      integer(kind=4) :: has_boundary
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: zz

      integer(kind=4), dimension(2) :: sreqs, rreqs
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: sstat, rstat
      integer(kind=4) :: ierr = 0

      nx = grid%my%nx
      ny = grid%my%ny
      has_boundary = grid%my%has_boundary

      zz => grid%zz

      !**************************************
      !*                                    *
      !* North-South Communication          *
      !*                                    *
      !**************************************
      north_edge => grid%edges%fne_d
      south_edge => grid%edges%fse_d
      north_buf  => grid%edges%fnb_d
      south_buf  => grid%edges%fsb_d

      if(iand(has_boundary, NORTH_BOUND) == 0) then
         north_edge(1:nx,1) = zz(1:nx,2)
      end if
      if(iand(has_boundary, SOUTH_BOUND) == 0) then
         south_edge(1:nx,1) = zz(1:nx,ny-1)
      end if

#ifndef MULTI
      call MPI_Irecv(north_buf,  nx, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(north_buf,  nx, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, NORTH_BOUND) == 0) then
         zz(1:nx,1) = north_buf(1:nx,1)
      end if
      if(iand(has_boundary, SOUTH_BOUND) == 0) then
         zz(1:nx,ny  ) = south_buf(1:nx,1)
      end if

      !**************************************
      !*                                    *
      !* East-West Communication            *
      !*                                    *
      !**************************************
      east_edge => grid%edges%fee_d
      west_edge => grid%edges%fwe_d
      east_buf  => grid%edges%feb_d
      west_buf  => grid%edges%fwb_d

      if(iand(has_boundary, EAST_BOUND) == 0) then
         east_edge(1:ny,1) = zz(nx-1,1:ny)
      end if
      if(iand(has_boundary, WEST_BOUND) == 0) then
         west_edge(1:ny,1) = zz(2,1:ny)
      end if

#ifndef MULTI
      call MPI_Irecv(east_buf,  ny, REAL_MPI, east_rank, 2, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, REAL_MPI, west_rank, 3, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, REAL_MPI, east_rank, 3, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, REAL_MPI, west_rank, 2, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(east_buf,  ny, REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, EAST_BOUND) == 0) then
         zz(nx,1:ny) = east_buf(1:ny,1)
      end if
      if(iand(has_boundary, WEST_BOUND) == 0) then
         zz(1,1:ny) = west_buf(1:ny,1)
      end if

      return
   end subroutine exchange_edges_zz
! ==============================================================================
! === Dispersive ===============================================================
   subroutine exchange_edges_disp_fx(grid, parity)
      type(data_grids), target, intent(inout) :: grid
      integer(kind=4), intent(in) :: parity

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_buf
      integer(kind=4) :: nx
      integer(kind=4) :: ny
      integer(kind=4) :: has_boundary
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx

      integer(kind=4), dimension(2) :: sreqs, rreqs
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: sstat, rstat
      integer(kind=4) :: ierr = 0

      nx = grid%my%nx
      ny = grid%my%ny
      has_boundary = grid%my%has_boundary

      if(parity == 1) then
         fx => grid%wave_field%fx
      else
         fx => grid%wave_field%fx_old
      end if

      !**************************************
      !*                                    *
      !* North-South Communication          *
      !*                                    *
      !**************************************
      north_edge => grid%edges%fne_d
      south_edge => grid%edges%fse_d
      north_buf  => grid%edges%fnb_d
      south_buf  => grid%edges%fsb_d

      if(iand(has_boundary, NORTH_BOUND) == 0) north_edge(1:nx,1) = fx(1:nx,2)
      if(iand(has_boundary, SOUTH_BOUND) == 0) south_edge(1:nx,1) = fx(1:nx,ny-1)

#ifndef MULTI
      call MPI_Irecv(north_buf,  nx, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(north_buf,  nx, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, NORTH_BOUND) == 0) fx(1:nx,1) = north_buf(1:nx,1)
      if(iand(has_boundary, SOUTH_BOUND) == 0) fx(1:nx,ny) = south_buf(1:nx,1)

      !**************************************
      !*                                    *
      !* East-West Communication            *
      !*                                    *
      !**************************************
      east_edge => grid%edges%fee_d
      west_edge => grid%edges%fwe_d
      east_buf  => grid%edges%feb_d
      west_buf  => grid%edges%fwb_d

      if(iand(has_boundary, EAST_BOUND) == 0) east_edge(1:ny,1) = fx(nx-1,1:ny)
      if(iand(has_boundary, WEST_BOUND) == 0) west_edge(1:ny,1) = fx(2,1:ny)

#ifndef MULTI
      call MPI_Irecv(east_buf,  ny, REAL_MPI, east_rank, 2, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, REAL_MPI, west_rank, 3, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, REAL_MPI, east_rank, 3, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, REAL_MPI, west_rank, 2, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(east_buf,  ny, REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, EAST_BOUND) == 0) fx(nx,1:ny) = east_buf(1:ny,1)
      if(iand(has_boundary, WEST_BOUND) == 0) fx(1,1:ny) = west_buf(1:ny,1)

      return
   end subroutine exchange_edges_disp_fx

   subroutine exchange_edges_disp_fy(grid, parity)
      type(data_grids), target, intent(inout) :: grid
      integer(kind=4), intent(in) :: parity

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_edge
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: north_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: south_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: east_buf
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: west_buf
      integer(kind=4) :: nx
      integer(kind=4) :: ny
      integer(kind=4) :: has_boundary
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fy

      integer(kind=4), dimension(2) :: sreqs, rreqs
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: sstat, rstat
      integer(kind=4) :: ierr = 0

      nx = grid%my%nx
      ny = grid%my%ny
      has_boundary = grid%my%has_boundary

      if(parity == 1) then
         fy => grid%wave_field%fy
      else
         fy => grid%wave_field%fy_old
      end if

      !**************************************
      !*                                    *
      !* North-South Communication          *
      !*                                    *
      !**************************************
      north_edge => grid%edges%fne_d
      south_edge => grid%edges%fse_d
      north_buf  => grid%edges%fnb_d
      south_buf  => grid%edges%fsb_d

      if(iand(has_boundary, NORTH_BOUND) == 0) north_edge(1:nx,1) = fy(1:nx,2)
      if(iand(has_boundary, SOUTH_BOUND) == 0) south_edge(1:nx,1) = fy(1:nx,ny-1)

#ifndef MULTI
      call MPI_Irecv(north_buf,  nx, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(north_buf,  nx, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, NORTH_BOUND) == 0) fy(1:nx,1) = north_buf(1:nx,1)
      if(iand(has_boundary, SOUTH_BOUND) == 0) fy(1:nx,ny) = south_buf(1:nx,1)

      !**************************************
      !*                                    *
      !* East-West Communication            *
      !*                                    *
      !**************************************
      east_edge => grid%edges%fee_d
      west_edge => grid%edges%fwe_d
      east_buf  => grid%edges%feb_d
      west_buf  => grid%edges%fwb_d

      if(iand(has_boundary, EAST_BOUND) == 0) east_edge(1:ny,1) = fy(nx-1,1:ny)
      if(iand(has_boundary, WEST_BOUND) == 0) west_edge(1:ny,1) = fy(2,1:ny)

#ifndef MULTI
      call MPI_Irecv(east_buf,  ny, REAL_MPI, east_rank, 2, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, REAL_MPI, west_rank, 3, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, REAL_MPI, east_rank, 3, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, REAL_MPI, west_rank, 2, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(east_buf,  ny, REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, EAST_BOUND) == 0) fy(nx,1:ny) = east_buf(1:ny,1)
      if(iand(has_boundary, WEST_BOUND) == 0) fy(1,1:ny) = west_buf(1:ny,1)

      return
   end subroutine exchange_edges_disp_fy
! ==============================================================================
! === Conversion from flux to velocity should be done right after calc. ========
   subroutine exchange_edges_wod(grid)
      type(data_grids), target, intent(inout) :: grid

      integer(kind=4), pointer, dimension(:,:) :: north_edge
      integer(kind=4), pointer, dimension(:,:) :: south_edge
      integer(kind=4), pointer, dimension(:,:) :: east_edge
      integer(kind=4), pointer, dimension(:,:) :: west_edge
      integer(kind=4), pointer, dimension(:,:) :: north_buf
      integer(kind=4), pointer, dimension(:,:) :: south_buf
      integer(kind=4), pointer, dimension(:,:) :: east_buf
      integer(kind=4), pointer, dimension(:,:) :: west_buf
      integer(kind=4) :: nx
      integer(kind=4) :: ny
      integer(kind=4) :: has_boundary
      integer(kind=4), pointer, dimension(:,:) :: wod

      integer(kind=4), dimension(2) :: sreqs, rreqs
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: sstat, rstat
      integer(kind=4) :: ierr = 0

      nx = grid%my%nx
      ny = grid%my%ny
      has_boundary = grid%my%has_boundary

      wod => grid%wod_flags

      !**************************************
      !*                                    *
      !* North-South Communication          *
      !*                                    *
      !**************************************
      north_edge => grid%edges%wne_d
      south_edge => grid%edges%wse_d
      north_buf  => grid%edges%wnb_d
      south_buf  => grid%edges%wsb_d

      if(iand(has_boundary, NORTH_BOUND) == 0) then
         north_edge(1:nx,1) = wod(1:nx,3)
      end if
      if(iand(has_boundary, SOUTH_BOUND) == 0) then
         south_edge(1:nx,1) = wod(1:nx,ny-2)
      end if

#ifndef MULTI
      call MPI_Irecv(north_buf,  nx, MPI_INTEGER, north_rank, 0, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, MPI_INTEGER, south_rank, 1, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, MPI_INTEGER, north_rank, 1, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, MPI_INTEGER, south_rank, 0, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(north_buf,  nx, MPI_INTEGER, north_rank, 0, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(south_buf,  nx, MPI_INTEGER, south_rank, 1, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(north_edge, nx, MPI_INTEGER, north_rank, 1, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(south_edge, nx, MPI_INTEGER, south_rank, 0, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, NORTH_BOUND) == 0) then
         wod(1:nx,0) = north_buf(1:nx,1)
      end if
      if(iand(has_boundary, SOUTH_BOUND) == 0) then
         wod(1:nx,ny+1) = south_buf(1:nx,1)
      end if

      !**************************************
      !*                                    *
      !* East-West Communication            *
      !*                                    *
      !**************************************
      east_edge => grid%edges%wee_d
      west_edge => grid%edges%wwe_d
      east_buf  => grid%edges%web_d
      west_buf  => grid%edges%wwb_d

      if(iand(has_boundary, EAST_BOUND) == 0) then
         east_edge(1:ny,1) = wod(nx-2,1:ny)
      end if
      if(iand(has_boundary, WEST_BOUND) == 0) then
         west_edge(1:ny,1) = wod(3,1:ny)
      end if

#ifndef MULTI
      call MPI_Irecv(east_buf,  ny, MPI_INTEGER, east_rank, 2, MPI_COMM_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, MPI_INTEGER, west_rank, 3, MPI_COMM_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, MPI_INTEGER, east_rank, 3, MPI_COMM_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, MPI_INTEGER, west_rank, 2, MPI_COMM_WORLD, sreqs(2), ierr)
#else
      call MPI_Irecv(east_buf,  ny, MPI_INTEGER, east_rank, 2, MPI_MEMBER_WORLD, rreqs(1), ierr)
      call MPI_Irecv(west_buf,  ny, MPI_INTEGER, west_rank, 3, MPI_MEMBER_WORLD, rreqs(2), ierr)
      call MPI_Isend(east_edge, ny, MPI_INTEGER, east_rank, 3, MPI_MEMBER_WORLD, sreqs(1), ierr)
      call MPI_Isend(west_edge, ny, MPI_INTEGER, west_rank, 2, MPI_MEMBER_WORLD, sreqs(2), ierr)
#endif
      call MPI_Waitall(2, sreqs, sstat, ierr)
      call MPI_Waitall(2, rreqs, rstat, ierr)

      if(iand(has_boundary, EAST_BOUND) == 0) then
         wod(nx+1,1:ny) = east_buf(1:ny,1)
      end if
      if(iand(has_boundary, WEST_BOUND) == 0) then
         wod(0,1:ny) = west_buf(1:ny,1)
      end if

      return
   end subroutine exchange_edges_wod
! ==============================================================================

end module mod_mpi
