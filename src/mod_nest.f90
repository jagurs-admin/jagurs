#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_nest
#ifdef MPI
use mpi
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
#endif
use mod_grid
! === 1-way nest ===============================================================
!use mod_params, only : VEL, HGT
use mod_params, only : VEL, HGT, nest_1way
! ==============================================================================
implicit none

contains

   subroutine mapgrids(mode,cg,fg,c2p_all)
      integer(kind=4), intent(in) :: mode
      type(data_grids), target, intent(inout) :: cg, fg
      integer(kind=4), intent(in) :: c2p_all
! === 1-way nest ===============================================================
!     call copy2coarse(mode,cg,fg,c2p_all)
      if(nest_1way == 0) call copy2coarse(mode,cg,fg,c2p_all)
! ==============================================================================
      call interp2fine(mode,cg,fg)
      return
   end subroutine mapgrids

#ifndef USE_ALLTOALLV
! === Config of copy from child to perent. by tkato 2012/11/15 =================
!  subroutine initl_gridmap(dg)
   subroutine initl_gridmap(dg, c2p_all)
! ==============================================================================
#else
! === USE_MPI_ALLTOALLV ========================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
!  subroutine initl_gridmap(dg, nprocs)
   subroutine initl_gridmap(dg, nprocs, c2p_all)
#ifdef SINGLE_A2A
#ifdef A2A3D
      use mod_a2a3d
#endif
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: nprocs
! === USE_MPI_ALLTOALLV ========================================================
#endif
      type(data_grids), target, intent(inout) :: dg
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      integer(kind=4), intent(in) :: c2p_all
! ==============================================================================

      type(grid_info), pointer :: parent, my
      type(interp_info), pointer :: fxo, fyo, fxi, fyi, hzi
      integer(kind=4) :: bigIX, bigIY, ix, iy, k
      real(kind=REAL_BYTE) :: lfac
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
      integer(kind=4) :: fxi_np0, fxi_np1
      integer(kind=4) :: fyi_np0, fyi_np1
#ifdef USE_ALLTOALLV
      integer(kind=4) :: ixst, iyst
      integer(kind=4) :: ixen, iyen
#ifdef SINGLE_A2A
! === SINGLE Alltoallv =========================================================
      integer(kind=4) :: i, p, sc, rc
      integer(kind=4) :: scx0, scy0, rcx0, rcy0
      integer(kind=4) :: scx1, scy1, rcx1, rcy1
! ==============================================================================
#endif
#endif
! ==============================================================================
! === copy2coarse for hz =======================================================
      type(interp_info), pointer :: hzo
! ==============================================================================

      parent => dg%parent
      my => dg%my

! === Config of copy from child to perent. by tkato 2012/11/15 =================
      if(c2p_all == 1) then
         write(6,'(a)') '*** NOTE: c2p_all is defined!'
! ==============================================================================
! === DEBUG for wave hight gap on nest boundary. 2012/10/31 ====================
!     dg%fxo%np = 2*(my%bigNY - 2)
      dg%fxo%np = (my%bigNX - 1)*(my%bigNY - 2)
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      else
         write(6,'(a)') '*** NOTE: c2p_all is NOT defined!'
! ==============================================================================
! === Reduce child -> parent copy by tkato 2012/11/13 ==========================
      dg%fxo%np = 4*(my%bigNY - 2) + 2*(my%bigNX - 5)
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      end if
! ==============================================================================
      allocate(dg%fxo%fndx(dg%fxo%np,2))
      allocate(dg%fxo%cndx0(dg%fxo%np,2))

! === Config of copy from child to perent. by tkato 2012/11/15 =================
      if(c2p_all == 1) then
! ==============================================================================
! === DEBUG for wave hight gap on nest boundary. 2012/10/31 ====================
!     dg%fyo%np = 2*(my%bigNX - 2)
      dg%fyo%np = (my%bigNY - 1)*(my%bigNX - 2)
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      else
! ==============================================================================
! === Reduce child -> parent copy by tkato 2012/11/13 ==========================
      dg%fyo%np = 4*(my%bigNX - 2) + 2*(my%bigNY - 5)
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      end if
! ==============================================================================
      allocate(dg%fyo%fndx(dg%fyo%np,2))
      allocate(dg%fyo%cndx0(dg%fyo%np,2))

      fxo => dg%fxo
      fyo => dg%fyo

#ifndef MPI
      ix = (my%nx-2) - my%nr/2 + 1
#else
      ix = (my%totalNx-2) - my%nr/2 + 1
#endif
      bigIX = my%zeroIX + my%bigNX - 2
      bigIY = my%zeroIY + 1
      k = 1
#ifndef MPI
      do iy = my%nr+1, my%ny-1, my%nr
#else
      do iy = my%nr+1, my%totalNy-1, my%nr
#endif
         fxo%fndx(k,1) = ix
         fxo%fndx(k,2) = iy
         fxo%cndx0(k,1) = bigIX
         fxo%cndx0(k,2) = bigIY
         bigIY = bigIY+1
         k = k+1
      end do

      ix = int(my%nr/2) + 1
      bigIX = my%zeroIX
      bigIY = my%zeroIY + my%bigNY - 2
#ifndef MPI
      do iy = my%ny - my%nr, my%nr, -my%nr
#else
      do iy = my%totalNy - my%nr, my%nr, -my%nr
#endif
         fxo%fndx(k,1) = ix
         fxo%fndx(k,2) = iy
         fxo%cndx0(k,1) = bigIX
         fxo%cndx0(k,2) = bigIY
         bigIY = bigIY-1
         k = k+1
      end do
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      if(c2p_all == 1) then
! ==============================================================================
! === DEBUG for wave hight gap on nest boundary. 2012/10/31 ====================
      bigIX = my%zeroIX + 1
#ifndef MPI
      do ix = int(my%nr/2) + 1 + my%nr, (my%nx-2) - my%nr/2 + 1 - my%nr, my%nr
         bigIY = my%zeroIY + 1
         do iy = my%nr+1, my%ny-1, my%nr
#else
      do ix = int(my%nr/2) + 1 + my%nr, (my%totalNx-2) - my%nr/2 + 1 - my%nr, my%nr
         bigIY = my%zeroIY + 1
         do iy = my%nr+1, my%totalNy-1, my%nr
#endif
            fxo%fndx(k,1) = ix
            fxo%fndx(k,2) = iy
            fxo%cndx0(k,1) = bigIX
            fxo%cndx0(k,2) = bigIY
            bigIY = bigIY+1
            k = k+1
         end do
         bigIX = bigIX+1
      end do
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      else
! ==============================================================================
! === Reduce child -> parent copy by tkato 2012/11/13 ==========================
#ifndef MPI
      ix = (my%nx-2) - my%nr/2 + 1 - my%nr
#else
      ix = (my%totalNx-2) - my%nr/2 + 1 - my%nr
#endif
      bigIX = my%zeroIX + my%bigNX - 3
      bigIY = my%zeroIY + 1
#ifndef MPI
      do iy = my%nr+1, my%ny-1, my%nr
#else
      do iy = my%nr+1, my%totalNy-1, my%nr
#endif
         fxo%fndx(k,1) = ix
         fxo%fndx(k,2) = iy
         fxo%cndx0(k,1) = bigIX
         fxo%cndx0(k,2) = bigIY
         bigIY = bigIY+1
         k = k+1
      end do

      ix = int(my%nr/2) + 1 + my%nr
      bigIX = my%zeroIX + 1
      bigIY = my%zeroIY + my%bigNY - 2
#ifndef MPI
      do iy = my%ny - my%nr, my%nr, -my%nr
#else
      do iy = my%totalNy - my%nr, my%nr, -my%nr
#endif
         fxo%fndx(k,1) = ix
         fxo%fndx(k,2) = iy
         fxo%cndx0(k,1) = bigIX
         fxo%cndx0(k,2) = bigIY
         bigIY = bigIY-1
         k = k+1
      end do

      iy = my%nr+1
      bigIX = my%zeroIX + 2
      bigIY = my%zeroIY + 1
#ifndef MPI
      do ix = int(my%nr/2) + 1 + 2*my%nr, (my%nx-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#else
      do ix = int(my%nr/2) + 1 + 2*my%nr, (my%totalNx-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#endif
         fxo%fndx(k,1) = ix
         fxo%fndx(k,2) = iy
         fxo%cndx0(k,1) = bigIX
         fxo%cndx0(k,2) = bigIY
         bigIX = bigIX+1
         k = k+1
      end do

#ifndef MPI
      iy = my%ny - my%nr
#else
      iy = my%totalNy - my%nr
#endif
      bigIX = my%zeroIX + 2
      bigIY = my%zeroIY + my%bigNY - 2
#ifndef MPI
      do ix = int(my%nr/2) + 1 + 2*my%nr, (my%nx-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#else
      do ix = int(my%nr/2) + 1 + 2*my%nr, (my%totalNx-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#endif
         fxo%fndx(k,1) = ix
         fxo%fndx(k,2) = iy
         fxo%cndx0(k,1) = bigIX
         fxo%cndx0(k,2) = bigIY
         bigIX = bigIX+1
         k = k+1
      end do
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      end if
! ==============================================================================
#ifdef MPI
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
      call make_a2a_params_c2c(dg, dg%fxo, nprocs)
! === USE_MPI_ALLTOALLV ========================================================
#endif
#endif

      iy = int(my%nr/2) + 1
      bigIX = my%zeroIX + 1
      bigIY = my%zeroIY
      k = 1
#ifndef MPI
      do ix = my%nr+1, my%nx-1, my%nr
#else
      do ix = my%nr+1, my%totalNx-1, my%nr
#endif
         fyo%fndx(k,1) = ix
         fyo%fndx(k,2) = iy
         fyo%cndx0(k,1) = bigIX
         fyo%cndx0(k,2) = bigIY
         bigIX = bigIX+1
         k = k+1
      end do

#ifndef MPI
      iy = (my%ny-2) - my%nr/2 + 1
#else
      iy = (my%totalNy-2) - my%nr/2 + 1
#endif
      bigIX = my%zeroIX + my%bigNX - 2
      bigIY = my%zeroIY + my%bigNY - 2
#ifndef MPI
      do ix = my%nx - my%nr, my%nr, -my%nr
#else
      do ix = my%totalNx - my%nr, my%nr, -my%nr
#endif
         fyo%fndx(k,1) = ix
         fyo%fndx(k,2) = iy
         fyo%cndx0(k,1) = bigIX
         fyo%cndx0(k,2) = bigIY
         bigIX = bigIX-1
         k = k+1
      end do
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      if(c2p_all == 1) then
! ==============================================================================
! === DEBUG for wave hight gap on nest boundary. 2012/10/31 ====================
      bigIY = my%zeroIY + 1
#ifndef MPI
      do iy = int(my%nr/2) + 1 + my%nr, (my%ny-2) - my%nr/2 + 1 - my%nr, my%nr
         bigIX = my%zeroIX + 1
         do ix = my%nr+1, my%nx-1, my%nr
#else
      do iy = int(my%nr/2) + 1 + my%nr, (my%totalNy-2) - my%nr/2 + 1 - my%nr, my%nr
         bigIX = my%zeroIX + 1
         do ix = my%nr+1, my%totalNx-1, my%nr
#endif
            fyo%fndx(k,1) = ix
            fyo%fndx(k,2) = iy
            fyo%cndx0(k,1) = bigIX
            fyo%cndx0(k,2) = bigIY
            bigIX = bigIX+1
            k = k+1
         end do
         bigIY = bigIY+1
      end do
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      else
! ==============================================================================
! === Reduce child -> parent copy by tkato 2012/11/13 ==========================
      iy = int(my%nr/2) + 1 + my%nr
      bigIX = my%zeroIX + 1
      bigIY = my%zeroIY + 1
#ifndef MPI
      do ix = my%nr+1, my%nx-1, my%nr
#else
      do ix = my%nr+1, my%totalNx-1, my%nr
#endif
         fyo%fndx(k,1) = ix
         fyo%fndx(k,2) = iy
         fyo%cndx0(k,1) = bigIX
         fyo%cndx0(k,2) = bigIY
         bigIX = bigIX+1
         k = k+1
      end do

#ifndef MPI
      iy = (my%ny-2) - my%nr/2 + 1 - my%nr
#else
      iy = (my%totalNy-2) - my%nr/2 + 1 - my%nr
#endif
      bigIX = my%zeroIX + my%bigNX - 2
      bigIY = my%zeroIY + my%bigNY - 3
#ifndef MPI
      do ix = my%nx - my%nr, my%nr, -my%nr
#else
      do ix = my%totalNx - my%nr, my%nr, -my%nr
#endif
         fyo%fndx(k,1) = ix
         fyo%fndx(k,2) = iy
         fyo%cndx0(k,1) = bigIX
         fyo%cndx0(k,2) = bigIY
         bigIX = bigIX-1
         k = k+1
      end do

      ix =  my%nr + 1
      bigIX = bigIX + 1
      bigIY = my%zeroIY + 2
#ifndef MPI
      do iy = int(my%nr/2) + 1 + 2*my%nr, (my%ny-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#else
      do iy = int(my%nr/2) + 1 + 2*my%nr, (my%totalNy-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#endif
         fyo%fndx(k,1) = ix
         fyo%fndx(k,2) = iy
         fyo%cndx0(k,1) = bigIX
         fyo%cndx0(k,2) = bigIY
         bigIY = bigIY+1
         k = k+1
      end do

#ifndef MPI
      ix =  my%nx - my%nr
#else
      ix =  my%totalNx - my%nr
#endif
      bigIX = my%zeroIX + my%bigNX - 2
      bigIY = my%zeroIY + 2
#ifndef MPI
      do iy = int(my%nr/2) + 1 + 2*my%nr, (my%ny-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#else
      do iy = int(my%nr/2) + 1 + 2*my%nr, (my%totalNy-2) - my%nr/2 + 1 - 2*my%nr, my%nr
#endif
         fyo%fndx(k,1) = ix
         fyo%fndx(k,2) = iy
         fyo%cndx0(k,1) = bigIX
         fyo%cndx0(k,2) = bigIY
         bigIY = bigIY+1
         k = k+1
      end do
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
      end if
! ==============================================================================
#ifdef MPI
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
      call make_a2a_params_c2c(dg, dg%fyo, nprocs)
! === USE_MPI_ALLTOALLV ========================================================
#ifdef SINGLE_A2A
! === SINGLE Alltoallv =========================================================
      allocate(dg%fxo%smap0(dg%fxo%snp0))
      allocate(dg%fyo%smap0(dg%fyo%snp0))
      allocate(dg%fxo%rmap0(dg%fxo%rnp0))
      allocate(dg%fyo%rmap0(dg%fyo%rnp0))
      sc = 1
      rc = 1
      scx0 = 1
      scy0 = 1
      rcx0 = 1
      rcy0 = 1
      do p = 0, nprocs-1
         do i = 1, dg%fxo%sendcnts0(p)
            dg%fxo%smap0(scx0) = sc
            scx0 = scx0 + 1
            sc = sc + 1
         end do
         do i = 1, dg%fyo%sendcnts0(p)
            dg%fyo%smap0(scy0) = sc
            scy0 = scy0 + 1
            sc = sc + 1
         end do
         do i = 1, dg%fxo%recvcnts0(p)
            dg%fxo%rmap0(rcx0) = rc
            rcx0 = rcx0 + 1
            rc = rc + 1
         end do
         do i = 1, dg%fyo%recvcnts0(p)
            dg%fyo%rmap0(rcy0) = rc
            rcy0 = rcy0 + 1
            rc = rc + 1
         end do
         dg%fyo%sendcnts0(p) = dg%fxo%sendcnts0(p) + dg%fyo%sendcnts0(p)
         dg%fyo%recvcnts0(p) = dg%fxo%recvcnts0(p) + dg%fyo%recvcnts0(p)
      end do
      dg%fyo%sdispls0(0) = 0
      dg%fyo%rdispls0(0) = 0
      do p = 1, nprocs-1
         dg%fyo%sdispls0(p) = dg%fyo%sdispls0(p-1) + dg%fyo%sendcnts0(p-1)
         dg%fyo%rdispls0(p) = dg%fyo%rdispls0(p-1) + dg%fyo%recvcnts0(p-1)
      end do
#ifdef A2A3D
      call A2A3D_plan(dg%fyo%sendcnts0, dg%fyo%handler0)
#endif
! ==============================================================================
#endif
#endif
#endif
! === copy2coarse for hz =======================================================
      if(c2p_all == 1) then
         dg%hzo%np = (my%bigNX - 2)*(my%bigNY - 2)

         allocate(dg%hzo%fndx(dg%hzo%np,2))
         allocate(dg%hzo%cndx0(dg%hzo%np,2))
         hzo => dg%hzo

         k = 1

         bigIX = my%zeroIX + 1
#ifndef MPI
         do ix = 1 + my%nr, my%nx - 1, my%nr
            bigIY = my%zeroIY + 1
            do iy = 1 + my%nr, my%ny - 1, my%nr
#else
         do ix = 1 + my%nr, my%totalNx - 1, my%nr
            bigIY = my%zeroIY + 1
            do iy = 1 + my%nr, my%totalNy - 1, my%nr
#endif
               hzo%fndx(k,1) = ix
               hzo%fndx(k,2) = iy
               hzo%cndx0(k,1) = bigIX
               hzo%cndx0(k,2) = bigIY
               bigIY = bigIY+1
               k = k+1
            end do
            bigIX = bigIX+1
         end do
#ifdef MPI
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
         call make_a2a_params_c2c(dg, dg%hzo, nprocs)
! === USE_MPI_ALLTOALLV ========================================================
#ifdef SINGLE_A2A
#ifdef A2A3D
      call A2A3D_plan(dg%hzo%sendcnts0, dg%hzo%handler0)
#endif
#endif
#endif
#endif
      end if
! ==============================================================================

      lfac = 1.0/REAL_FUNC(my%nr)

#ifndef MPI
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!     dg%fxi%np = 2*(my%nx - 1)
      fxi_np0 = (my%ny - 1)/my%nr + 1 ! for i = 0
      fxi_np1 = (my%ny - 1)/my%nr + 1 ! for i = nx
! === Upwind3 ==================================================================
!     dg%fxi%np = 2*(my%nx - 1) + fxi_np0 + fxi_np1 + (my%nx + 1) + (my%ny + 1)
      dg%fxi%np = 2*(my%nx - 1) + fxi_np0 + fxi_np1 + (my%nx + 1) + (my%ny + 1) &
                + 2*((my%nx - 1)/my%nr + 2) + fxi_np0 + my%ny
! ==============================================================================
! ==============================================================================
#else
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!     dg%fxi%np = 2*(my%totalNx - 1)
      fxi_np0 = (my%totalNy - 1)/my%nr + 1 ! for i = 0
      fxi_np1 = (my%totalNy - 1)/my%nr + 1 ! for i = nx
! === Upwind3 ==================================================================
!     dg%fxi%np = 2*(my%totalNx - 1) + fxi_np0 + fxi_np1 + (my%totalNx + 1) + (my%totalNy + 1)
      dg%fxi%np = 2*(my%totalNx - 1) + fxi_np0 + fxi_np1 + (my%totalNx + 1) + (my%totalNy + 1) &
                + 2*((my%totalNx - 1)/my%nr + 2) + fxi_np0 + my%totalNy
! ==============================================================================
! ==============================================================================
#endif
      allocate(dg%fxi%fndx(dg%fxi%np,2))
      allocate(dg%fxi%cndx0(dg%fxi%np,2))
      allocate(dg%fxi%cndx1(dg%fxi%np,2))
      allocate(dg%fxi%wt0(dg%fxi%np))
      allocate(dg%fxi%wt1(dg%fxi%np))

#ifndef MPI
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!     dg%fyi%np = 2*(my%ny - 1)
      fyi_np0 = (my%nx - 1)/my%nr + 1 ! for j = 0
      fyi_np1 = (my%nx - 1)/my%nr + 1 ! for j = ny
! === Upwind3 ==================================================================
!     dg%fyi%np = 2*(my%ny - 1) + fyi_np0 + fyi_np1 + (my%ny + 1) + (my%nx + 1)
      dg%fyi%np = 2*(my%ny - 1) + fyi_np0 + fyi_np1 + (my%ny + 1) + (my%nx + 1) &
                + 2*((my%ny - 1)/my%nr + 2) + fyi_np0 + my%nx
! ==============================================================================
! ==============================================================================
#else
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!     dg%fyi%np = 2*(my%totalNy - 1)
      fyi_np0 = (my%totalNx - 1)/my%nr + 1 ! for j = 0
      fyi_np1 = (my%totalNx - 1)/my%nr + 1 ! for j = ny
! === Upwind3 ==================================================================
!     dg%fyi%np = 2*(my%totalNy - 1) + fyi_np0 + fyi_np1 + (my%totalNy + 1) + (my%totalNx + 1)
      dg%fyi%np = 2*(my%totalNy - 1) + fyi_np0 + fyi_np1 + (my%totalNy + 1) + (my%totalNx + 1) &
                + 2*((my%totalNy - 1)/my%nr + 2) + fyi_np0 + my%totalNx
! ==============================================================================
! ==============================================================================
#endif
      allocate(dg%fyi%fndx(dg%fyi%np,2))
      allocate(dg%fyi%cndx0(dg%fyi%np,2))
      allocate(dg%fyi%cndx1(dg%fyi%np,2))
      allocate(dg%fyi%wt0(dg%fyi%np))
      allocate(dg%fyi%wt1(dg%fyi%np))

#ifndef MPI
      dg%hzi%np = 2*((my%nx - 1) + (my%ny - 1)) + 2*my%nx + 2*(my%ny + 2)
#else
      dg%hzi%np = 2*((my%totalNx - 1) + (my%totalNy - 1)) + 2*my%totalNx + 2*(my%totalNy + 2)
#endif
      allocate(dg%hzi%fndx(dg%hzi%np,2))
      allocate(dg%hzi%cndx0(dg%hzi%np,2))
      allocate(dg%hzi%cndx1(dg%hzi%np,2))
      allocate(dg%hzi%wt0(dg%hzi%np))
      allocate(dg%hzi%wt1(dg%hzi%np))

      fxi => dg%fxi
      fyi => dg%fyi
      hzi => dg%hzi

      iy = 1
      k = 1
#ifndef MPI
      do ix = 1, my%nx-1
#else
      do ix = 1, my%totalNx-1
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%wt0(k) = 1.0d0 - mod(ix+1,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do

#ifndef MPI
      iy = my%ny
      do ix = my%nx-1, 1, -1
#else
      iy = my%totalNy
      do ix = my%totalNx-1, 1, -1
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%wt0(k) = 1.0d0 - mod(ix+1,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do

      iy = 0
#ifndef MPI
      do ix = 0, my%nx
#else
      do ix = 0, my%totalNx
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx0(k,2) = my%zeroIY - 1 + (iy-1)/my%nr
         fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr
         fxi%cndx1(k,2) = my%zeroIY - 1 + (iy-1)/my%nr
         fxi%wt0(k) = 1.0d0 - mod(ix+1,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do

      ix = 0
      do iy = 1, fxi_np0
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = (iy-1)*my%nr + 1
         fxi%cndx0(k,1) = my%zeroIX - 1
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)
         fxi%cndx1(k,1) = my%zeroIX
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)
         fxi%wt0(k) = 2.0d0*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do

#ifndef MPI
      ix = my%nx
#else
      ix = my%totalNx
#endif
      do iy = 1, fxi_np1
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = (iy-1)*my%nr + 1
#ifndef MPI
         fxi%cndx0(k,1) = my%zeroIX + (my%nx-1)/my%nr - 1
#else
         fxi%cndx0(k,1) = my%zeroIX + (my%totalNx-1)/my%nr - 1
#endif
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)
#ifndef MPI
         fxi%cndx1(k,1) = my%zeroIX + (my%nx-1)/my%nr
#else
         fxi%cndx1(k,1) = my%zeroIX + (my%totalNx-1)/my%nr
#endif
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)
         fxi%wt0(k) = 1.0d0 - 2.0d0*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do

      ix = -1

      iy = 0
      fxi%fndx(k,1) = ix
      fxi%fndx(k,2) = iy
      fxi%cndx0(k,1) = my%zeroIX - 1
      fxi%cndx0(k,2) = my%zeroIY - 1
      fxi%cndx1(k,1) = my%zeroIX - 1
      fxi%cndx1(k,2) = my%zeroIY
      fxi%wt0(k) = lfac
      fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
      k = k+1

#ifndef MPI
      do iy = 1, my%ny
#else
      do iy = 1, my%totalNy
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX - 1
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%cndx1(k,1) = my%zeroIX - 1
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         fxi%wt0(k) = 1.0d0 - mod(iy+2,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do

! === Upwind3 ==================================================================
! === fx(1:nx,-1) ==============================================================
      iy = -1

      ix = 1
      fxi%fndx(k,1) = ix
      fxi%fndx(k,2) = iy
      fxi%cndx0(k,1) = my%zeroIX - 1
      fxi%cndx0(k,2) = my%zeroIY
      fxi%cndx1(k,1) = my%zeroIX - 1
      fxi%cndx1(k,2) = my%zeroIY - 1
      fxi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
      k = k+1
#ifndef MPI
      do ix = 2, my%nx, my%nr
#else
      do ix = 2, my%totalNx, my%nr
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx0(k,2) = my%zeroIY
         fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx1(k,2) = my%zeroIY - 1
         fxi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do
#ifndef MPI
      ix = my%nx
#else
      ix = my%totalNx
#endif
      fxi%fndx(k,1) = ix
      fxi%fndx(k,2) = iy
      fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr
      fxi%cndx0(k,2) = my%zeroIY
      fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr
      fxi%cndx1(k,2) = my%zeroIY - 1
      fxi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
      k = k+1
! === fx(1:nx,ny+1) ============================================================
#ifndef MPI
      iy = my%ny+1
#else
      iy = my%totalNy+1
#endif

      ix = 1
      fxi%fndx(k,1) = ix
      fxi%fndx(k,2) = iy
      fxi%cndx0(k,1) = my%zeroIX - 1
      fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr + 1
      fxi%cndx1(k,1) = my%zeroIX - 1
      fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
      fxi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
      k = k+1
#ifndef MPI
      do ix = 2, my%nx, my%nr
#else
      do ix = 2, my%totalNx, my%nr
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do
#ifndef MPI
      ix = my%nx
#else
      ix = my%totalNx
#endif
      fxi%fndx(k,1) = ix
      fxi%fndx(k,2) = iy
      fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr
      fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr + 1
      fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr
      fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
      fxi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
      k = k+1
! === fx(-2:1:ny) ==============================================================
      ix = -2
#ifndef MPI
      do iy = 1, my%ny, my%nr
#else
      do iy = 1, my%totalNy, my%nr
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX - 2
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%cndx1(k,1) = my%zeroIX - 1
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do
! === fx(nx+1:1:ny) ============================================================
#ifndef MPI
      ix = my%nx+1
      do iy = 1, my%ny
#else
      ix = my%totalNx+1
      do iy = 1, my%totalNy
#endif
         fxi%fndx(k,1) = ix
         fxi%fndx(k,2) = iy
         fxi%cndx0(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         fxi%cndx1(k,1) = my%zeroIX + (ix+1)/my%nr - 1
         fxi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         fxi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         fxi%wt1(k) = 1.0d0 - fxi%wt0(k)
         k = k+1
      end do
! ==============================================================================
#ifdef MPI
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
      ixst = my%kx
      iyst = my%ky
! === Upwind3 ==================================================================
!     if(iand(my%has_boundary, WEST_BOUND)  /= 0) ixst = -1
!     if(iand(my%has_boundary, NORTH_BOUND) /= 0) iyst = 0
      if(iand(my%has_boundary, WEST_BOUND) /= 0) ixst = -2
      if(iand(my%has_boundary, NORTH_BOUND) /= 0) iyst = -1
! ==============================================================================
      ixen = my%kxend
      iyen = my%kyend
! === Upwind3 ==================================================================
      if(iand(my%has_boundary, EAST_BOUND) /= 0) ixen = my%kxend + 1
      if(iand(my%has_boundary, SOUTH_BOUND) /= 0) iyen = my%kyend + 1
! ==============================================================================
      call make_a2a_params_i2f(dg, dg%fxi, nprocs, ixst, iyst, ixen, iyen)
! === USE_MPI_ALLTOALLV ========================================================
#endif
#endif

#ifndef MPI
      ix = my%nx
#else
      ix = my%totalNx
#endif
      k = 1
#ifndef MPI
      do iy = 1, my%ny-1
#else
      do iy = 1, my%totalNy-1
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr
         fyi%wt0(k) = 1.0d0 - mod(iy+1,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do

      ix = 1
#ifndef MPI
      do iy = my%ny-1, 1, -1
#else
      do iy = my%totalNy-1, 1, -1
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr
         fyi%wt0(k) = 1.0d0 - mod(iy+1,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do

      ix = 0
#ifndef MPI
      do iy = my%ny, 0, -1
#else
      do iy = my%totalNy, 0, -1
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX - 1 + (ix-1)/my%nr
         fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%cndx1(k,1) = my%zeroIX - 1 + (ix-1)/my%nr
         fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr
         fyi%wt0(k) = 1.0d0 - mod(iy+1,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do

      iy = 0
      do ix = 1, fyi_np0
         fyi%fndx(k,1) = (ix-1)*my%nr + 1
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)
         fyi%cndx0(k,2) = my%zeroIY - 1
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)
         fyi%cndx1(k,2) = my%zeroIY
         fyi%wt0(k) = 2.0d0*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do

#ifndef MPI
      iy = my%ny
#else
      iy = my%totalNy
#endif
      do ix = 1, fyi_np1
         fyi%fndx(k,1) = (ix-1)*my%nr + 1
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)
#ifndef MPI
         fyi%cndx0(k,2) = my%zeroIY + (my%ny-1)/my%nr - 1
#else
         fyi%cndx0(k,2) = my%zeroIY + (my%totalNy-1)/my%nr - 1
#endif
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)
#ifndef MPI
         fyi%cndx1(k,2) = my%zeroIY + (my%ny-1)/my%nr
#else
         fyi%cndx1(k,2) = my%zeroIY + (my%totalNy-1)/my%nr
#endif
         fyi%wt0(k) = 1.0d0 - 2.0d0*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do

      iy = -1

      ix = 0
      fyi%fndx(k,1) = ix
      fyi%fndx(k,2) = iy
      fyi%cndx0(k,1) = my%zeroIX - 1
      fyi%cndx0(k,2) = my%zeroIY - 1
      fyi%cndx1(k,1) = my%zeroIX
      fyi%cndx1(k,2) = my%zeroIY - 1
      fyi%wt0(k) = lfac
      fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
      k = k+1

#ifndef MPI
      do ix = 1, my%nx
#else
      do ix = 1, my%totalNx
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx0(k,2) = my%zeroIY - 1
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         fyi%cndx1(k,2) = my%zeroIY - 1
         fyi%wt0(k) = 1.0d0 - mod(ix+2,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do

! === Upwind3 ==================================================================
! === fy(-1,1:ny) ==============================================================
      ix = -1

      iy = 1
      fyi%fndx(k,1) = ix
      fyi%fndx(k,2) = iy
      fyi%cndx0(k,1) = my%zeroIX
      fyi%cndx0(k,2) = my%zeroIY - 1
      fyi%cndx1(k,1) = my%zeroIX - 1
      fyi%cndx1(k,2) = my%zeroIY - 1
      fyi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
      k = k+1
#ifndef MPI
      do iy = 2, my%ny, my%nr
#else
      do iy = 2, my%totalNy, my%nr
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX
         fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%cndx1(k,1) = my%zeroIX - 1
         fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do
#ifndef MPI
      iy = my%ny
#else
      iy = my%totalNy
#endif
      fyi%fndx(k,1) = ix
      fyi%fndx(k,2) = iy
      fyi%cndx0(k,1) = my%zeroIX
      fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr
      fyi%cndx1(k,1) = my%zeroIX - 1
      fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr
      fyi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
      k = k+1
! === fy(nx+1,1:ny) ============================================================
#ifndef MPI
      ix = my%nx+1
#else
      ix = my%totalNx+1
#endif

      iy = 1
      fyi%fndx(k,1) = ix
      fyi%fndx(k,2) = iy
      fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr + 1
      fyi%cndx0(k,2) = my%zeroIY - 1
      fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
      fyi%cndx1(k,2) = my%zeroIY - 1
      fyi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
      k = k+1
#ifndef MPI
      do iy = 2, my%ny, my%nr
#else
      do iy = 2, my%totalNy, my%nr
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do
#ifndef MPI
      iy = my%ny
#else
      iy = my%totalNy
#endif
      fyi%fndx(k,1) = ix
      fyi%fndx(k,2) = iy
      fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr + 1
      fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr
      fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
      fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr
      fyi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
      fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
      k = k+1
! === fy(1:nx,-2) ==============================================================
      iy = -2
#ifndef MPI
      do ix = 1, my%nx, my%nr
#else
      do ix = 1, my%totalNx, my%nr
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx0(k,2) = my%zeroIY - 2
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx1(k,2) = my%zeroIY - 1
         fyi%wt0(k) = 1.0d0 - mod(2,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do
! === fy(1:nx,ny+1) ============================================================
#ifndef MPI
      iy = my%ny+1
      do ix = 1, my%nx
#else
      iy = my%totalNy+1
      do ix = 1, my%totalNx
#endif
         fyi%fndx(k,1) = ix
         fyi%fndx(k,2) = iy
         fyi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         fyi%cndx0(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         fyi%cndx1(k,2) = my%zeroIY + (iy+1)/my%nr - 1
         fyi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         fyi%wt1(k) = 1.0d0 - fyi%wt0(k)
         k = k+1
      end do
! ==============================================================================
#ifdef MPI
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
      ixst = my%kx
      iyst = my%ky
! === Upwind3 ==================================================================
!     if(iand(my%has_boundary, WEST_BOUND)  /= 0) ixst = 0
!     if(iand(my%has_boundary, NORTH_BOUND) /= 0) iyst = -1
      if(iand(my%has_boundary, WEST_BOUND) /= 0) ixst = -1
      if(iand(my%has_boundary, NORTH_BOUND) /= 0) iyst = -2
! ==============================================================================
      ixen = my%kxend
      iyen = my%kyend
! === Upwind3 ==================================================================
      if(iand(my%has_boundary, EAST_BOUND) /= 0) ixen = my%kxend + 1
      if(iand(my%has_boundary, SOUTH_BOUND) /= 0) iyen = my%kyend + 1
! ==============================================================================
      call make_a2a_params_i2f(dg, dg%fyi, nprocs, ixst, iyst, ixen, iyen)
! === USE_MPI_ALLTOALLV ========================================================
#ifdef SINGLE_A2A
! === SINGLE Alltoallv =========================================================
      allocate(dg%fxi%smap0(dg%fxi%snp0))
      allocate(dg%fxi%smap1(dg%fxi%snp1))
      allocate(dg%fyi%smap0(dg%fyi%snp0))
      allocate(dg%fyi%smap1(dg%fyi%snp1))
      allocate(dg%fxi%rmap0(dg%fxi%rnp0))
      allocate(dg%fxi%rmap1(dg%fxi%rnp1))
      allocate(dg%fyi%rmap0(dg%fyi%rnp0))
      allocate(dg%fyi%rmap1(dg%fyi%rnp1))
      sc = 1
      rc = 1
      scx0 = 1
      scx1 = 1
      scy0 = 1
      scy1 = 1
      rcx0 = 1
      rcx1 = 1
      rcy0 = 1
      rcy1 = 1
      do p = 0, nprocs-1
         do i = 1, dg%fxi%sendcnts0(p)
            dg%fxi%smap0(scx0) = sc
            scx0 = scx0 + 1
            sc = sc + 1
         end do
         do i = 1, dg%fxi%sendcnts1(p)
            dg%fxi%smap1(scx1) = sc
            scx1 = scx1 + 1
            sc = sc + 1
         end do
         do i = 1, dg%fyi%sendcnts0(p)
            dg%fyi%smap0(scy0) = sc
            scy0 = scy0 + 1
            sc = sc + 1
         end do
         do i = 1, dg%fyi%sendcnts1(p)
            dg%fyi%smap1(scy1) = sc
            scy1 = scy1 + 1
            sc = sc + 1
         end do
         do i = 1, dg%fxi%recvcnts0(p)
            dg%fxi%rmap0(rcx0) = rc
            rcx0 = rcx0 + 1
            rc = rc + 1
         end do
         do i = 1, dg%fxi%recvcnts1(p)
            dg%fxi%rmap1(rcx1) = rc
            rcx1 = rcx1 + 1
            rc = rc + 1
         end do
         do i = 1, dg%fyi%recvcnts0(p)
            dg%fyi%rmap0(rcy0) = rc
            rcy0 = rcy0 + 1
            rc = rc + 1
         end do
         do i = 1, dg%fyi%recvcnts1(p)
            dg%fyi%rmap1(rcy1) = rc
            rcy1 = rcy1 + 1
            rc = rc + 1
         end do
         dg%fyi%sendcnts1(p) = dg%fxi%sendcnts0(p) + dg%fyi%sendcnts0(p) &
                             + dg%fxi%sendcnts1(p) + dg%fyi%sendcnts1(p)
         dg%fyi%recvcnts1(p) = dg%fxi%recvcnts0(p) + dg%fyi%recvcnts0(p) &
                             + dg%fxi%recvcnts1(p) + dg%fyi%recvcnts1(p)
      end do
      dg%fyi%sdispls1(0) = 0
      dg%fyi%rdispls1(0) = 0
      do p = 1, nprocs-1
         dg%fyi%sdispls1(p) = dg%fyi%sdispls1(p-1) + dg%fyi%sendcnts1(p-1)
         dg%fyi%rdispls1(p) = dg%fyi%rdispls1(p-1) + dg%fyi%recvcnts1(p-1)
      end do
#ifdef A2A3D
      call A2A3D_plan(dg%fyi%sendcnts1, dg%fyi%handler1)
#endif
! ==============================================================================
#endif
#endif
#endif

      iy = 1
      k = 1
#ifndef MPI
      do ix = 1, my%nx-1
#else
      do ix = 1, my%totalNx-1
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         hzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%wt0(k) = 1.0d0 - mod(ix-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do

#ifndef MPI
      ix = my%nx
      do iy = 1, my%ny-1
#else
      ix = my%totalNx
      do iy = 1, my%totalNy-1
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         hzi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do

#ifndef MPI
      iy = my%ny
      do ix = my%nx, 2, -1
#else
      iy = my%totalNy
      do ix = my%totalNx, 2, -1
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         hzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%wt0(k) = 1.0d0 - mod(ix-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do

      ix = 1
#ifndef MPI
      do iy = my%ny, 2, -1
#else
      do iy = my%totalNy, 2, -1
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         hzi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do

      iy = 0
#ifndef MPI
      do ix = 1, my%nx
#else
      do ix = 1, my%totalNx
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx0(k,2) = my%zeroIY - 1
         hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         hzi%cndx1(k,2) = my%zeroIY - 1
         hzi%wt0(k) = 1.0d0 - mod(ix-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do

#ifndef MPI
      iy = my%ny + 1
      do ix = my%nx, 1, -1
#else
      iy = my%totalNy + 1
      do ix = my%totalNx, 1, -1
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         hzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         hzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         hzi%wt0(k) = 1.0d0 - mod(ix-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do

      ix = 0
#ifndef MPI
      do iy = my%ny + 1, 1, -1
#else
      do iy = my%totalNy + 1, 1, -1
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX - 1
         hzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%cndx1(k,1) = my%zeroIX - 1
         hzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         hzi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do
      iy = 0
      hzi%fndx(k,1) = ix
      hzi%fndx(k,2) = iy
      hzi%cndx0(k,1) = my%zeroIX - 1
      hzi%cndx0(k,2) = my%zeroIY - 1
      hzi%cndx1(k,1) = my%zeroIX - 1
      hzi%cndx1(k,2) = my%zeroIY
      hzi%wt0(k) = lfac
      hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
      k = k+1

#ifndef MPI
      ix = my%nx + 1
#else
      ix = my%totalNx + 1
#endif
      iy = 0
      hzi%fndx(k,1) = ix
      hzi%fndx(k,2) = iy
      hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr + 1
      hzi%cndx0(k,2) = my%zeroIY
      hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
      hzi%cndx1(k,2) = my%zeroIY - 1
      hzi%wt0(k) = 1.0d0 - lfac
      hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
      k = k+1
#ifndef MPI
      do iy = 1, my%ny + 1
#else
      do iy = 1, my%totalNy + 1
#endif
         hzi%fndx(k,1) = ix
         hzi%fndx(k,2) = iy
         hzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         hzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         hzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         hzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         hzi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         hzi%wt1(k) = 1.0d0 - hzi%wt0(k)
         k = k+1
      end do

#ifdef MPI
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
      ixst = my%kx
      iyst = my%ky
! === Upwind3 ==================================================================
!     if(iand(my%has_boundary, WEST_BOUND)  /= 0) ixst = 0
!     if(iand(my%has_boundary, NORTH_BOUND) /= 0) iyst = 0
      if(iand(my%has_boundary, WEST_BOUND) /= 0) ixst = -1
      if(iand(my%has_boundary, NORTH_BOUND) /= 0) iyst = -1
! ==============================================================================
      ixen = my%kxend + 1
      iyen = my%kyend + 1
! === Upwind3 ==================================================================
      if(iand(my%has_boundary, EAST_BOUND) /= 0) ixen = my%kxend + 2
      if(iand(my%has_boundary, SOUTH_BOUND) /= 0) iyen = my%kyend + 2
! ==============================================================================
      call make_a2a_params_i2f(dg, dg%hzi, nprocs, ixst, iyst, ixen, iyen)
! === USE_MPI_ALLTOALLV ========================================================
#ifdef SINGLE_A2A
! === SINGLE Alltoallv =========================================================
      allocate(dg%hzi%smap0(dg%hzi%snp0))
      allocate(dg%hzi%smap1(dg%hzi%snp1))
      allocate(dg%hzi%rmap0(dg%hzi%rnp0))
      allocate(dg%hzi%rmap1(dg%hzi%rnp1))
      sc = 1
      rc = 1
      scx0 = 1
      scx1 = 1
      rcx0 = 1
      rcx1 = 1
      do p = 0, nprocs-1
         do i = 1, dg%hzi%sendcnts0(p)
            dg%hzi%smap0(scx0) = sc
            scx0 = scx0 + 1
            sc = sc + 1
         end do
         do i = 1, dg%hzi%sendcnts1(p)
            dg%hzi%smap1(scx1) = sc
            scx1 = scx1 + 1
            sc = sc + 1
         end do
         do i = 1, dg%hzi%recvcnts0(p)
            dg%hzi%rmap0(rcx0) = rc
            rcx0 = rcx0 + 1
            rc = rc + 1
         end do
         do i = 1, dg%hzi%recvcnts1(p)
            dg%hzi%rmap1(rcx1) = rc
            rcx1 = rcx1 + 1
            rc = rc + 1
         end do
         dg%hzi%sendcnts1(p) = dg%hzi%sendcnts0(p) + dg%hzi%sendcnts1(p)
         dg%hzi%recvcnts1(p) = dg%hzi%recvcnts0(p) + dg%hzi%recvcnts1(p)
      end do
      dg%hzi%sdispls1(0) = 0
      dg%hzi%rdispls1(0) = 0
      do p = 1, nprocs-1
         dg%hzi%sdispls1(p) = dg%hzi%sdispls1(p-1) + dg%hzi%sendcnts1(p-1)
         dg%hzi%rdispls1(p) = dg%hzi%rdispls1(p-1) + dg%hzi%recvcnts1(p-1)
      end do
#ifdef A2A3D
      call A2A3D_plan(dg%hzi%sendcnts1, dg%hzi%handler1)
#endif
! ==============================================================================
#endif
#endif
#endif
! === Don't repeat allocate/deallocate! ========================================
#ifdef MPI
#ifndef USE_ALLTOALLV
      ! for copy2coarse
      allocate(dg%fxo%fb(dg%fxo%np+dg%fyo%np))
      allocate(dg%fxo%cb(dg%fxo%np+dg%fyo%np))
      allocate(dg%hzo%fb(dg%hzo%np))
      allocate(dg%hzo%cb(dg%hzo%np))

      ! for interp2fine
      allocate(dg%fxi%fb(dg%fxi%np+dg%fyi%np))
      allocate(dg%fxi%cb(dg%fxi%np+dg%fyi%np))
      allocate(dg%hzi%fb(dg%hzi%np))
      allocate(dg%hzi%cb(dg%hzi%np))
#else
#ifndef SINGLE_A2A
      ! for copy2coarse
      allocate(dg%fxo%fb(dg%fxo%snp0))
      allocate(dg%fxo%cb(dg%fxo%rnp0))
      allocate(dg%fyo%fb(dg%fyo%snp0))
      allocate(dg%fyo%cb(dg%fyo%rnp0))
      allocate(dg%hzo%fb(dg%hzo%snp0))
      allocate(dg%hzo%cb(dg%hzo%rnp0))

      ! for interp2fine
      allocate(dg%fxi%fb0(dg%fxi%rnp0))
      allocate(dg%fxi%fb1(dg%fxi%rnp1))
      allocate(dg%fxi%cb0(dg%fxi%snp0))
      allocate(dg%fxi%cb1(dg%fxi%snp1))
      allocate(dg%fyi%fb0(dg%fyi%rnp0))
      allocate(dg%fyi%fb1(dg%fyi%rnp1))
      allocate(dg%fyi%cb0(dg%fyi%snp0))
      allocate(dg%fyi%cb1(dg%fyi%snp1))
      allocate(dg%hzi%fb0(dg%hzi%rnp0))
      allocate(dg%hzi%fb1(dg%hzi%rnp1))
      allocate(dg%hzi%cb0(dg%hzi%snp0))
      allocate(dg%hzi%cb1(dg%hzi%snp1))
#else
      ! for copy2coarse
      allocate(dg%fyo%fb(dg%fxo%snp0+dg%fyo%snp0))
      allocate(dg%fyo%cb(dg%fxo%rnp0+dg%fyo%rnp0))
      allocate(dg%hzo%fb(dg%hzo%snp0))
      allocate(dg%hzo%cb(dg%hzo%rnp0))

      ! for interp2fine
      allocate(dg%fyi%fb1(dg%fxi%rnp0+dg%fxi%rnp1+dg%fyi%rnp0+dg%fyi%rnp1))
      allocate(dg%fyi%cb1(dg%fxi%snp0+dg%fxi%snp1+dg%fyi%snp0+dg%fyi%snp1))
      allocate(dg%hzi%fb1(dg%hzi%rnp0+dg%hzi%rnp1))
      allocate(dg%hzi%cb1(dg%hzi%snp0+dg%hzi%snp1))
#endif
#endif
#endif
! ==============================================================================

      return
   end subroutine initl_gridmap

   subroutine copy2coarse(mode,cg,fg,c2p_all)
#ifdef USE_ALLTOALLV
#ifdef A2A3D
      use mod_a2a3d
#endif
#endif
      integer(kind=4), intent(in) :: mode
      type(data_grids), target, intent(inout) :: cg, fg
      integer(kind=4), intent(in) :: c2p_all

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fxc, fyc
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fxf, fyf
      type(interp_info), pointer :: fxo
      type(interp_info), pointer :: fyo
      integer(kind=4) :: k
#ifdef MPI
#ifndef USE_ALLTOALLV
      real(kind=REAL_BYTE), pointer, dimension(:) :: xfbuf, xcbuf
      integer(kind=4) :: i, j, ix, iy, ierr
#else
! === USE_MPI_ALLTOALLV ========================================================
#ifndef SINGLE_A2A
      real(kind=REAL_BYTE), pointer, dimension(:) :: xfbuf, xcbuf, yfbuf, ycbuf
#else
      real(kind=REAL_BYTE), pointer, dimension(:) :: yfbuf, ycbuf
#endif
#ifndef A2A3D
      integer(kind=4) :: i, j, ierr
#else
      integer(kind=4) :: i, j
#endif
! === USE_MPI_ALLTOALLV ========================================================
#endif
#endif
! === copy2coarse for hz =======================================================
      type(interp_info), pointer :: hzo
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hzf, hzc
#ifdef MPI
      real(kind=REAL_BYTE), pointer, dimension(:) :: zfbuf, zcbuf
#else
      integer(kind=4) :: i, j, ic, jc
#endif
      integer(kind=4), pointer, dimension(:,:) :: wodf, wodc
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dzc
      real(kind=REAL_BYTE) :: tmp
      integer(kind=4) :: i_, j_
! ==============================================================================

      fxc => cg%wave_field%fx
      fyc => cg%wave_field%fy

      fxf => fg%wave_field%fx
      fyf => fg%wave_field%fy

      fxo => fg%fxo
      fyo => fg%fyo

      if(mode == VEL) then
#ifndef MPI
!$omp parallel
         if(c2p_all == 0) then
!$omp do
            do k = 1, fxo%np
               fxc(fxo%cndx0(k,1),fxo%cndx0(k,2)) = fxf(fxo%fndx(k,1),fxo%fndx(k,2))
            end do
!$omp do
            do k = 1, fyo%np
               fyc(fyo%cndx0(k,1),fyo%cndx0(k,2)) = fyf(fyo%fndx(k,1),fyo%fndx(k,2))
            end do
         else
!$omp do private(i, j)
            do k = 1, fxo%np
               i = fxo%fndx(k,1)
               j = fxo%fndx(k,2)
               fxc(fxo%cndx0(k,1),fxo%cndx0(k,2)) = &
               &  (fxf(i-1,j-1) + fxf(i,j-1) + fxf(i+1,j-1) &
               & + fxf(i-1,j  ) + fxf(i,j  ) + fxf(i+1,j  ) &
               & + fxf(i-1,j+1) + fxf(i,j+1) + fxf(i+1,j+1))/9.0d0
            end do
!$omp do private(i, j)
            do k = 1, fyo%np
               i = fyo%fndx(k,1)
               j = fyo%fndx(k,2)
               fyc(fyo%cndx0(k,1),fyo%cndx0(k,2)) = &
               &  (fyf(i-1,j-1) + fyf(i,j-1) + fyf(i+1,j-1) &
               & + fyf(i-1,j  ) + fyf(i,j  ) + fyf(i+1,j  ) &
               & + fyf(i-1,j+1) + fyf(i,j+1) + fyf(i+1,j+1))/9.0d0
            end do
         end if
!$omp end parallel
#else
#ifndef USE_ALLTOALLV
!$omp parallel
!$omp single
         xfbuf => fxo%fb
         xcbuf => fxo%cb
!$omp end single
!$omp do
         do k = 1, fxo%np+fyo%np
            xfbuf(k) = 0.0d0
            xcbuf(k) = 0.0d0
         end do
         !*==============*
         !*  fine2buf    * must not read from edges.
         !*==============*
         if(c2p_all == 0) then
!$omp do private(ix, iy, i, j)
            do k = 1, fxo%np
               ix = fxo%fndx(k,1)
               iy = fxo%fndx(k,2)
               if(fg%my%ix <= ix .and. ix <= fg%my%ixend .and. fg%my%iy <= iy .and. iy <= fg%my%iyend) then
                  i = ix - fg%my%kx + 1 ! with egde
                  j = iy - fg%my%ky + 1 ! with edge
                  xfbuf(k) = fxf(i,j)
               end if
            end do
!$omp do private(ix, iy, i, j)
            do k = 1, fyo%np
               ix = fyo%fndx(k,1)
               iy = fyo%fndx(k,2)
               if(fg%my%ix <= ix .and. ix <= fg%my%ixend .and. fg%my%iy <= iy .and. iy <= fg%my%iyend) then
                  i = ix - fg%my%kx + 1 ! with egde 
                  j = iy - fg%my%ky + 1 ! with edge 
                  xfbuf(fxo%np+k) = fyf(i,j)
               end if
            end do
         else
!$omp do private(ix, iy, i, j)
            do k = 1, fxo%np
               ix = fxo%fndx(k,1)
               iy = fxo%fndx(k,2)
               if(fg%my%ix <= ix .and. ix <= fg%my%ixend .and. fg%my%iy <= iy .and. iy <= fg%my%iyend) then
                  i = ix - fg%my%kx + 1 ! with egde
                  j = iy - fg%my%ky + 1 ! with edge
                  xfbuf(k) = (fxf(i-1,j-1) + fxf(i,j-1) + fxf(i+1,j-1) &
                          & + fxf(i-1,j  ) + fxf(i,j  ) + fxf(i+1,j  ) &
                          & + fxf(i-1,j+1) + fxf(i,j+1) + fxf(i+1,j+1))/9.0d0
               end if
            end do
!$omp do private(ix, iy, i, j)
            do k = 1, fyo%np
               ix = fyo%fndx(k,1)
               iy = fyo%fndx(k,2)
               if(fg%my%ix <= ix .and. ix <= fg%my%ixend .and. fg%my%iy <= iy .and. iy <= fg%my%iyend) then
                  i = ix - fg%my%kx + 1 ! with egde 
                  j = iy - fg%my%ky + 1 ! with edge 
                  xfbuf(fxo%np+k) = (fyf(i-1,j-1) + fyf(i,j-1) + fyf(i+1,j-1) &
                          & + fyf(i-1,j  ) + fyf(i,j  ) + fyf(i+1,j  ) &
                          & + fyf(i-1,j+1) + fyf(i,j+1) + fyf(i+1,j+1))/9.0d0
               end if
            end do
         end if
         !*==============*
         !*  allreduce   *
         !*==============*
!$omp single
#ifndef MULTI
         call MPI_Allreduce(xfbuf, xcbuf, fxo%np+fyo%np, REAL_MPI, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
         call MPI_Allreduce(xfbuf, xcbuf, fxo%np+fyo%np, REAL_MPI, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
         if(ierr /= 0) then
            select case (ierr)
               case(MPI_ERR_BUFFER)
                  write(0,'(a)') 'MPI Error : Invalid buffer pointer'
               case(MPI_ERR_COUNT)
                  write(0,'(a)') 'MPI Error : Invalid count argument'
               case(MPI_ERR_TYPE)
                  write(0,'(a)') 'MPI Error : Invalid datatype argument'
               case(MPI_ERR_OP)
                  write(0,'(a)') 'MPI Error : Invalid operation'
               case(MPI_ERR_COMM)
                  write(0,'(a)') 'MPI Error : Invalid communicator'
               case default
                  write(0,'(a)') 'MPI Error : Unknown error'
            end select
            call fatal_error(ierr)
         end if
!$omp end single
         !*==============*
         !*  buf2coarce  * must write to edges.
         !*==============*
!$omp do private(ix, iy, i, j)
         do k = 1, fxo%np
            ix = fxo%cndx0(k,1)
            iy = fxo%cndx0(k,2)
            if(cg%my%kx <= ix .and. ix <= cg%my%kxend .and. cg%my%ky <= iy .and. iy <= cg%my%kyend) then
               i = ix - cg%my%kx + 1 ! with egde 
               j = iy - cg%my%ky + 1 ! with edge 
               fxc(i,j) = xcbuf(k)
            end if
         end do
!$omp do private(ix, iy, i, j)
         do k = 1, fyo%np
            ix = fyo%cndx0(k,1)
            iy = fyo%cndx0(k,2)
            if(cg%my%kx <= ix .and. ix <= cg%my%kxend .and. cg%my%ky <= iy .and. iy <= cg%my%kyend) then
               i = ix - cg%my%kx + 1 ! with egde 
               j = iy - cg%my%ky + 1 ! with edge 
               fyc(i,j) = xcbuf(fxo%np+k)
            end if
         end do
!$omp end parallel
#else
! === USE_MPI_ALLTOALLV ========================================================
#ifndef SINGLE_A2A
         xfbuf => fxo%fb
         xcbuf => fxo%cb
         yfbuf => fyo%fb
         ycbuf => fyo%cb
#else
         yfbuf => fyo%fb
         ycbuf => fyo%cb
#endif
         !*==============*
         !*  fine2buf    * must not read from edges.
         !*==============*
         if(c2p_all == 0) then
            do k = 1, fxo%snp0
               i = fxo%fndx0_l(k,1)
               j = fxo%fndx0_l(k,2)
#ifndef SINGLE_A2A
               xfbuf(k) = fxf(i,j)
#else
               yfbuf(fxo%smap0(k)) = fxf(i,j)
#endif
            end do
            do k = 1, fyo%snp0
               i = fyo%fndx0_l(k,1)
               j = fyo%fndx0_l(k,2)
#ifndef SINGLE_A2A
               yfbuf(k) = fyf(i,j)
#else
               yfbuf(fyo%smap0(k)) = fyf(i,j)
#endif
            end do
         else
!$omp parallel do private(i, j)
            do k = 1, fxo%snp0
               i = fxo%fndx0_l(k,1)
               j = fxo%fndx0_l(k,2)
#ifndef SINGLE_A2A
               xfbuf(k) = (fxf(i-1,j-1) + fxf(i,j-1) + fxf(i+1,j-1) &
#else
               yfbuf(fxo%smap0(k)) = (fxf(i-1,j-1) + fxf(i,j-1) + fxf(i+1,j-1) &
#endif
                       & + fxf(i-1,j  ) + fxf(i,j  ) + fxf(i+1,j  ) &
                       & + fxf(i-1,j+1) + fxf(i,j+1) + fxf(i+1,j+1))/9.0d0
            end do
!$omp parallel do private(i, j)
            do k = 1, fyo%snp0
               i = fyo%fndx0_l(k,1)
               j = fyo%fndx0_l(k,2)
#ifndef SINGLE_A2A
               yfbuf(k) = (fyf(i-1,j-1) + fyf(i,j-1) + fyf(i+1,j-1) &
#else
               yfbuf(fyo%smap0(k)) = (fyf(i-1,j-1) + fyf(i,j-1) + fyf(i+1,j-1) &
#endif
                       & + fyf(i-1,j  ) + fyf(i,j  ) + fyf(i+1,j  ) &
                       & + fyf(i-1,j+1) + fyf(i,j+1) + fyf(i+1,j+1))/9.0d0
            end do
         end if
         !*==============*
         !*  alltoallv   *
         !*==============*
#ifndef A2A3D
#ifndef MULTI
#ifndef SINGLE_A2A
         call MPI_Alltoallv(xfbuf, fxo%sendcnts0, fxo%sdispls0, REAL_MPI, &
                            xcbuf, fxo%recvcnts0, fxo%rdispls0, REAL_MPI, MPI_COMM_WORLD, ierr)
#endif
         call MPI_Alltoallv(yfbuf, fyo%sendcnts0, fyo%sdispls0, REAL_MPI, &
                            ycbuf, fyo%recvcnts0, fyo%rdispls0, REAL_MPI, MPI_COMM_WORLD, ierr)
#else
#ifndef SINGLE_A2A
         call MPI_Alltoallv(xfbuf, fxo%sendcnts0, fxo%sdispls0, REAL_MPI, &
                            xcbuf, fxo%recvcnts0, fxo%rdispls0, REAL_MPI, MPI_MEMBER_WORLD, ierr)
#endif
         call MPI_Alltoallv(yfbuf, fyo%sendcnts0, fyo%sdispls0, REAL_MPI, &
                            ycbuf, fyo%recvcnts0, fyo%rdispls0, REAL_MPI, MPI_MEMBER_WORLD, ierr)
#endif
#else
#ifndef SINGLE_A2A
         call A2A3D_execute(fxo%snp0, xfbuf, fxo%rnp0, xcbuf, fxo%handler0)
         call A2A3D_execute(fyo%snp0, yfbuf, fyo%rnp0, ycbuf, fyo%handler0)
#else
         call A2A3D_execute(fxo%snp0+fyo%snp0, yfbuf, fxo%rnp0+fyo%rnp0, ycbuf, fyo%handler0)
#endif
#endif
         !*==============*
         !*  buf2coarce  * must write to edges.
         !*==============*
         do k = 1, fxo%rnp0
            i = fxo%cndx0_l(k,1)
            j = fxo%cndx0_l(k,2)
#ifndef SINGLE_A2A
            fxc(i,j) = xcbuf(k)
#else
            fxc(i,j) = ycbuf(fxo%rmap0(k))
#endif
         end do
         do k = 1, fyo%rnp0
            i = fyo%cndx0_l(k,1)
            j = fyo%cndx0_l(k,2)
#ifndef SINGLE_A2A
            fyc(i,j) = ycbuf(k)
#else
            fyc(i,j) = ycbuf(fyo%rmap0(k))
#endif
         end do
! === USE_MPI_ALLTOALLV ========================================================
#endif
#endif
      end if
! === copy2coarse for hz =======================================================
      if(c2p_all == 1) then
         hzc => cg%wave_field%hz

         hzf => fg%wave_field%hz

         hzo => fg%hzo

         wodf => fg%wod_flags
         wodc => cg%wod_flags
         dzc  => cg%depth_field%dz

         if(mode == HGT) then
#ifndef MPI
!$omp parallel
!$omp do private(i, j, ic, jc, tmp, i_, j_)
            do k = 1, hzo%np
               i = hzo%fndx(k,1)
               j = hzo%fndx(k,2)
               ic = hzo%cndx0(k,1)
               jc = hzo%cndx0(k,2)
               if(wodc(ic,jc) == 1) then ! Only coarse domain is wet
                  tmp = 0.0d0
                  do j_ = -1, 1
                     do i_ = -1, 1
                        if(wodf(i+i_, j+j_) == 1) then ! Only fine domain is wet
                           tmp = tmp + hzf(i+i_, j+j_)
                        end if
                     end do
                  end do
                  hzc(ic,jc) = tmp/9.0d0
                  if(dzc(ic,jc) + hzc(ic,jc) > 0.01d0) then
                     wodc(ic,jc) = 1
                  else
                     wodc(ic,jc) = -1
                  end if
               end if
            end do
!$omp end parallel
#else
#ifndef USE_ALLTOALLV
!$omp parallel
!$omp single
            zfbuf => hzo%fb
            zcbuf => hzo%cb
!$omp end single
!$omp do
            do k = 1, hzo%np
               zfbuf(k) = 0.0d0
               zcbuf(k) = 0.0d0
            end do
            !*==============*
            !*  fine2buf    * must not read from edges.
            !*==============*
!$omp do private(ix, iy, i, j, tmp, i_, j_)
            do k = 1, hzo%np
               ix = hzo%fndx(k,1)
               iy = hzo%fndx(k,2)
               if(fg%my%ix <= ix .and. ix <= fg%my%ixend .and. fg%my%iy <= iy .and. iy <= fg%my%iyend) then
                  i = ix - fg%my%kx + 1 ! with egde
                  j = iy - fg%my%ky + 1 ! with edge
                  tmp = 0.0d0
                  do j_ = -1, 1
                     do i_ = -1, 1
                        if(wodf(i+i_, j+j_) == 1) then ! Only fine domain is wet
                           tmp = tmp + hzf(i+i_, j+j_)
                        end if
                     end do
                  end do
                  zfbuf(k) = tmp/9.0d0
               end if
            end do
            !*==============*
            !*  allreduce   *
            !*==============*
!$omp single
#ifndef MULTI
            call MPI_Allreduce(zfbuf, zcbuf, hzo%np, REAL_MPI, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            call MPI_Allreduce(zfbuf, zcbuf, hzo%np, REAL_MPI, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
            if(ierr /= 0) then
               select case (ierr)
                  case(MPI_ERR_BUFFER)
                     write(0,'(a)') 'MPI Error : Invalid buffer pointer'
                  case(MPI_ERR_COUNT)
                     write(0,'(a)') 'MPI Error : Invalid count argument'
                  case(MPI_ERR_TYPE)
                     write(0,'(a)') 'MPI Error : Invalid datatype argument'
                  case(MPI_ERR_OP)
                     write(0,'(a)') 'MPI Error : Invalid operation'
                  case(MPI_ERR_COMM)
                     write(0,'(a)') 'MPI Error : Invalid communicator'
                  case default
                     write(0,'(a)') 'MPI Error : Unknown error'
               end select
               call fatal_error(ierr)
            end if
!$omp end single
            !*==============*
            !*  buf2coarce  * must write to edges.
            !*==============*
!$omp do private(ix, iy, i, j)
            do k = 1, hzo%np
               ix = hzo%cndx0(k,1)
               iy = hzo%cndx0(k,2)
               if(cg%my%kx <= ix .and. ix <= cg%my%kxend .and. cg%my%ky <= iy .and. iy <= cg%my%kyend) then
                  i = ix - cg%my%kx + 1 ! with egde 
                  j = iy - cg%my%ky + 1 ! with edge 
                  if(wodc(i,j) == 1) then ! Only coarse domain is wet
                     hzc(i,j) = zcbuf(k)
                     if(dzc(i,j) + hzc(i,j) > 0.01d0) then
                        wodc(i,j) = 1
                     else
                        wodc(i,j) = -1
                     end if
                  end if
               end if
            end do
!$omp end parallel
#else
! === USE_MPI_ALLTOALLV ========================================================
            zfbuf => hzo%fb
            zcbuf => hzo%cb
            !*==============*
            !*  fine2buf    * must not read from edges.
            !*==============*
!$omp parallel do private(i, j, tmp, i_, j_)
            do k = 1, hzo%snp0
               i = hzo%fndx0_l(k,1)
               j = hzo%fndx0_l(k,2)
               tmp = 0.0d0
               do j_ = -1, 1
                  do i_ = -1, 1
                     if(wodf(i+i_, j+j_) == 1) then ! Only fine domain is wet
                        tmp = tmp + hzf(i+i_, j+j_)
                     end if
                  end do
               end do
               zfbuf(k) = tmp/9.0d0
            end do
            !*==============*
            !*  alltoallv   *
            !*==============*
#ifndef A2A3D
#ifndef MULTI
            call MPI_Alltoallv(zfbuf, hzo%sendcnts0, hzo%sdispls0, REAL_MPI, &
                               zcbuf, hzo%recvcnts0, hzo%rdispls0, REAL_MPI, MPI_COMM_WORLD, ierr)
#else
            call MPI_Alltoallv(zfbuf, hzo%sendcnts0, hzo%sdispls0, REAL_MPI, &
                               zcbuf, hzo%recvcnts0, hzo%rdispls0, REAL_MPI, MPI_MEMBER_WORLD, ierr)
#endif
#else
            call A2A3D_execute(hzo%snp0, zfbuf, hzo%rnp0, zcbuf, hzo%handler0)
#endif
            !*==============*
            !*  buf2coarce  * must write to edges.
            !*==============*
            do k = 1, hzo%rnp0
               i = hzo%cndx0_l(k,1)
               j = hzo%cndx0_l(k,2)
               if(wodc(i,j) == 1) then ! Only coarse domain is wet
                  hzc(i,j) = zcbuf(k)
                  if(dzc(i,j) + hzc(i,j) > 0.01d0) then
                     wodc(i,j) = 1
                  else
                     wodc(i,j) = -1
                  end if
               end if
            end do
! === USE_MPI_ALLTOALLV ========================================================
#endif
#endif
         end if
      end if
! ==============================================================================

      return
   end subroutine copy2coarse

   subroutine interp2fine(mode,cg,fg)
#ifdef USE_ALLTOALLV
#ifdef A2A3D
      use mod_a2a3d
#endif
#endif
      integer(kind=4), intent(in) :: mode
      type(data_grids), target, intent(inout) :: cg, fg

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fxc, fyc, hzc, dzc
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fxf, fyf, hzf, dzf
      integer(kind=4) :: k
      type(interp_info), pointer :: fxi
      type(interp_info), pointer :: fyi
      type(interp_info), pointer :: hzi
#ifdef MPI
#ifndef USE_ALLTOALLV
      real(kind=REAL_BYTE), pointer, dimension(:) :: xfbuf, xcbuf, zfbuf, zcbuf
#else
! === USE_MPI_ALLTOALLV ========================================================
#ifndef SINGLE_A2A
      real(kind=REAL_BYTE), pointer, dimension(:) :: xfbuf0, xcbuf0, yfbuf0, ycbuf0, zfbuf0, zcbuf0
      real(kind=REAL_BYTE), pointer, dimension(:) :: xfbuf1, xcbuf1, yfbuf1, ycbuf1, zfbuf1, zcbuf1
#else
      real(kind=REAL_BYTE), pointer, dimension(:) :: yfbuf1, ycbuf1, zfbuf1, zcbuf1
#endif
! === USE_MPI_ALLTOALLV ========================================================
#endif
      integer(kind=4) :: i, j, ix, iy, ierr
      integer(kind=4), dimension(2) :: ireq11, ireq12, ireq13, ireq14
      integer(kind=4), dimension(2) :: ireq21, ireq22, ireq23, ireq24
#endif
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
#ifndef MPI
      integer(kind=4) :: i, j, imod, ind0, ind1
#else
#ifndef USE_ALLTOALLV
      integer(kind=4) :: imod, ind0, ind1, ixst, iyst, ist, ien, jst, jen, &
                         shift_st, shift_en
! === Upwind3 ==================================================================
      integer(kind=4) :: ixen, iyen
! ==============================================================================
#else
      integer(kind=4) :: imod, ind0, ind1, ist, ien, jst, jen, &
                         shift_st, shift_en
#endif
! === Upwind3 ==================================================================
!     real(kind=REAL_BYTE) :: send_to_w, send_to_e, recv_from_w, recv_from_e
!     real(kind=REAL_BYTE) :: send_to_n, send_to_s, recv_from_n, recv_from_s
!     integer(kind=4) :: nx, ny, px, py, rx, ry, north_rank, south_rank, east_rank, west_rank
      real(kind=REAL_BYTE), dimension(3) :: send_to_w, send_to_e, recv_from_w, recv_from_e
      real(kind=REAL_BYTE), dimension(3) :: send_to_n, send_to_s, recv_from_n, recv_from_s
      integer(kind=4) :: px, py, rx, ry, north_rank, south_rank, east_rank, west_rank
      integer(kind=4) :: ist2, ien2, jst2, jen2, shift_st2, shift_en2
! ==============================================================================
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: stat11, stat12, stat13, stat14
      integer(kind=4), dimension(MPI_STATUS_SIZE,2) :: stat21, stat22, stat23, stat24
#endif
      real(kind=REAL_BYTE) :: lfac, fac0, fac1
! ==============================================================================
      real(kind=REAL_BYTE) :: t0, t1
! === Upwind3 ==================================================================
      integer(kind=4) :: nx, ny
! ==============================================================================
   
      fxc => cg%wave_field%fx
      fyc => cg%wave_field%fy
      hzc => cg%wave_field%hz
      dzc => cg%depth_field%dz

      fxf => fg%wave_field%fx
      fyf => fg%wave_field%fy
      hzf => fg%wave_field%hz
      dzf => fg%depth_field%dz

      fxi => fg%fxi
      fyi => fg%fyi
      hzi => fg%hzi

      if(mode == VEL) then
#ifndef MPI
!$omp parallel
!$omp do
         do k = 1, fxi%np
            fxf(fxi%fndx(k,1),fxi%fndx(k,2)) = &
               fxc(fxi%cndx0(k,1),fxi%cndx0(k,2))*fxi%wt0(k) + &
               fxc(fxi%cndx1(k,1),fxi%cndx1(k,2))*fxi%wt1(k)
         end do
!$omp do
         do k = 1, fyi%np
            fyf(fyi%fndx(k,1),fyi%fndx(k,2)) = &
               fyc(fyi%cndx0(k,1),fyi%cndx0(k,2))*fyi%wt0(k) + &
               fyc(fyi%cndx1(k,1),fyi%cndx1(k,2))*fyi%wt1(k)
         end do
!$omp end parallel
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
         lfac = 1.0d0/REAL_FUNC(fg%my%nr)

         do j = 1, fg%my%ny-1
            imod = mod(j-1,fg%my%nr)
            ind0 = j - imod
            ind1 = ind0 + fg%my%nr
            fac0 = 1.0d0 - imod*lfac
            fac1 = 1.0d0 - fac0
            fxf(0,       j) = fxf(0,       ind0)*fac0 + fxf(0,       ind1)*fac1
            fxf(fg%my%nx,j) = fxf(fg%my%nx,ind0)*fac0 + fxf(fg%my%nx,ind1)*fac1
         end do

         do i = 1, fg%my%nx-1
            imod = mod(i-1,fg%my%nr)
            ind0 = i - imod
            ind1 = ind0 + fg%my%nr
            fac0 = 1.0d0 - imod*lfac
            fac1 = 1.0d0 - fac0
            fyf(i,0)        = fyf(ind0,0       )*fac0 + fyf(ind1,0       )*fac1
            fyf(i,fg%my%ny) = fyf(ind0,fg%my%ny)*fac0 + fyf(ind1,fg%my%ny)*fac1
         end do
! ==============================================================================
         t0 = 1.0d0/REAL_FUNC(fg%my%nr)
         t1 = 1.0d0 - t0
         do i = 0, fg%my%nx
            fxf(i,0) = fxf(i,0)*t0 + fxf(i,1)*t1
         end do
         do j = 0, fg%my%ny
            fyf(0,j) = fyf(0,j)*t0 + fyf(1,j)*t1
         end do
! === Upwind3 ==================================================================
         nx = fg%my%nx
         ny = fg%my%ny

         fac0 = lfac
         fac1 = 1.0d0 - fac0
         fxf(1,      -1) = fxf(1,     -1)*fac0 + fxf(2,   -1)*fac1
         fxf(1,    ny+1) = fxf(1,   ny+1)*fac0 + fxf(2, ny+1)*fac1
         fxf(nx-1,   -1) = fxf(nx-2,  -1)*fac1 + fxf(nx,  -1)*fac0
         fxf(nx-1, ny+1) = fxf(nx-2,ny+1)*fac1 + fxf(nx,ny+1)*fac0
         fxf(nx,     -1) = fxf(nx-2,  -1)*fac0 + fxf(nx,  -1)*fac1
         fxf(nx,   ny+1) = fxf(nx-2,ny+1)*fac0 + fxf(nx,ny+1)*fac1
         do i = 2, nx - fg%my%nr
            imod = mod(i-2,fg%my%nr)
            ind0 = i - imod
            ind1 = ind0 + fg%my%nr
            fac0 = 1.0d0 - imod*lfac
            fac1 = 1.0d0 - fac0
            fxf(i,  -1) = fxf(ind0,  -1)*fac0 + fxf(ind1,  -1)*fac1
            fxf(i,ny+1) = fxf(ind0,ny+1)*fac0 + fxf(ind1,ny+1)*fac1
         end do

         do j = 1, ny - 1
            imod = mod(j-1,fg%my%nr)
            ind0 = j - imod
            ind1 = ind0 + fg%my%nr
            fac0 = 1.0d0 - imod*lfac
            fac1 = 1.0d0 - fac0
            fxf(-2,j) = fxf(-2,ind0)*fac0 + fxf(-2,ind1)*fac1
         end do

         fac0 = lfac
         fac1 = 1.0d0 - fac0
         fyf(-1,      1) = fyf(-1,     1)*fac0 + fyf(-1,   2)*fac1
         fyf(nx+1,    1) = fyf(nx+1,   1)*fac0 + fyf(nx+1, 2)*fac1
         fyf(-1,   ny-1) = fyf(-1,  ny-2)*fac1 + fyf(-1,  ny)*fac0
         fyf(nx+1, ny-1) = fyf(nx+1,ny-2)*fac1 + fyf(nx+1,ny)*fac0
         fyf(-1,     ny) = fyf(-1,  ny-2)*fac0 + fyf(-1,  ny)*fac1
         fyf(nx+1,   ny) = fyf(nx+1,ny-2)*fac0 + fyf(nx+1,ny)*fac1
         do j = 2, ny - fg%my%nr
            imod = mod(j-2,fg%my%nr)
            ind0 = j - imod
            ind1 = ind0 + fg%my%nr
            fac0 = 1.0d0 - imod*lfac
            fac1 = 1.0d0 - fac0
            fyf(-1,  j) = fyf(-1,  ind0)*fac0 + fyf(-1,  ind1)*fac1
            fyf(nx+1,j) = fyf(nx+1,ind0)*fac0 + fyf(nx+1,ind1)*fac1
         end do

         do i = 1, nx - 1
            imod = mod(i-1,fg%my%nr)
            ind0 = i - imod
            ind1 = ind0 + fg%my%nr
            fac0 = 1.0d0 - imod*lfac
            fac1 = 1.0d0 - fac0
            fyf(i,-2) = fyf(ind0,-2)*fac0 + fyf(ind1,-2)*fac1
         end do
! ==============================================================================
#else
#ifndef USE_ALLTOALLV
!$omp parallel
!$omp single
         xfbuf => fxi%fb
         xcbuf => fxi%cb
!$omp end single
!$omp do
         do k = 1, fxi%np+fyi%np
            xfbuf(k) = 0.0d0
            xcbuf(k) = 0.0d0
         end do
         !*==============*
         !*  coarse2buf  * must not read from edges
         !*==============*
!$omp do private(ix, iy, i, j)
         do k = 1, fxi%np
            ix = fxi%cndx0(k,1)
            iy = fxi%cndx0(k,2)
            if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
               i = ix - cg%my%kx + 1 ! with edges 
               j = iy - cg%my%ky + 1 ! with edges
               xcbuf(k) = xcbuf(k) + fxc(i,j)*fxi%wt0(k)
            end if
         end do
!$omp do private(ix, iy, i, j)
         do k = 1, fxi%np
            ix = fxi%cndx1(k,1)
            iy = fxi%cndx1(k,2)
            if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
               i = ix - cg%my%kx + 1 ! with edges 
               j = iy - cg%my%ky + 1 ! with edges 
               xcbuf(k) = xcbuf(k) + fxc(i,j)*fxi%wt1(k)
            end if
         end do
!$omp do private(ix, iy, i, j)
         do k = 1, fyi%np
            ix = fyi%cndx0(k,1)
            iy = fyi%cndx0(k,2)
            if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
               i = ix - cg%my%kx + 1 ! with edges 
               j = iy - cg%my%ky + 1 ! with edges 
               xcbuf(fxi%np+k) = xcbuf(fxi%np+k) + fyc(i,j)*fyi%wt0(k)
            end if
         end do
!$omp do private(ix, iy, i, j)
         do k = 1, fyi%np
            ix = fyi%cndx1(k,1)
            iy = fyi%cndx1(k,2)
            if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
               i = ix - cg%my%kx + 1 ! with edges 
               j = iy - cg%my%ky + 1 ! with edges 
               xcbuf(fxi%np+k) = xcbuf(fxi%np+k) + fyc(i,j)*fyi%wt1(k)
            end if
         end do
         !*==============*
         !*  allreduce   *
         !*==============*
!$omp single
#ifndef MULTI
         call MPI_Allreduce(xcbuf, xfbuf, fxi%np+fyi%np, REAL_MPI, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
         call MPI_Allreduce(xcbuf, xfbuf, fxi%np+fyi%np, REAL_MPI, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
         if(ierr /= 0) then
            select case (ierr)
               case(MPI_ERR_BUFFER)
                  write(0,'(a)') 'MPI Error : Invalid buffer pointer'
               case(MPI_ERR_COUNT)
                  write(0,'(a)') 'MPI Error : Invalid count argument'
               case(MPI_ERR_TYPE)
                  write(0,'(a)') 'MPI Error : Invalid datatype argument'
               case(MPI_ERR_OP)
                  write(0,'(a)') 'MPI Error : Invalid operation'
               case(MPI_ERR_COMM)
                  write(0,'(a)') 'MPI Error : Invalid communicator'
               case default
                  write(0,'(a)') 'MPI Error : Unknown error'
            end select
            call fatal_error(ierr)
         end if
!$omp end single
         !*==============*
         !*  buf2fine    * write to edges
         !*==============*
!$omp single
         ixst = fg%my%kx
         iyst = fg%my%ky
! === Upwind3 ==================================================================
!        if(iand(fg%my%has_boundary, WEST_BOUND)  /= 0) ixst = -1
!        if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) iyst = 0
         if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) ixst = -2
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) iyst = -1
         ixen = fg%my%kxend
         iyen = fg%my%kyend
         if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) ixen = fg%my%kxend+1
         if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) iyen = fg%my%kyend+1
! ==============================================================================
!$omp end single
!$omp do private(ix, iy, i, j)
         do k = 1, fxi%np
            ix = fxi%fndx(k,1)
            iy = fxi%fndx(k,2)
! === Upwind3 ==================================================================
!           if(ixst <= ix .and. ix <= fg%my%kxend .and. iyst <= iy .and. iy <= fg%my%kyend) then
            if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
               i = ix - fg%my%kx + 1 ! with edges
               j = iy - fg%my%ky + 1 ! with edges
               fxf(i,j) = xfbuf(k)
            end if
         end do
!$omp single
         ixst = fg%my%kx
         iyst = fg%my%ky
! === Upwind3 ==================================================================
!        if(iand(fg%my%has_boundary, WEST_BOUND)  /= 0) ixst = 0
!        if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) iyst = -1
         if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) ixst = -1
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) iyst = -2
         ixen = fg%my%kxend
         iyen = fg%my%kyend
         if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) ixen = fg%my%kxend+1
         if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) iyen = fg%my%kyend+1
! ==============================================================================
!$omp end single
!$omp do private(ix, iy, i, j)
         do k = 1, fyi%np
            ix = fyi%fndx(k,1)
            iy = fyi%fndx(k,2)
! === Upwind3 ==================================================================
!           if(ixst <= ix .and. ix <= fg%my%kxend .and. iyst <= iy .and. iy <= fg%my%kyend) then
            if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
               i = ix - fg%my%kx + 1 ! with edges
               j = iy - fg%my%ky + 1 ! with edges
               fyf(i,j) = xfbuf(fxi%np+k)
            end if
         end do
!$omp end parallel
#else
! === USE_MPI_ALLTOALLV ========================================================
#ifndef SINGLE_A2A
         xfbuf0 => fxi%fb0
         xfbuf1 => fxi%fb1
         xcbuf0 => fxi%cb0
         xcbuf1 => fxi%cb1
         yfbuf0 => fyi%fb0
         yfbuf1 => fyi%fb1
         ycbuf0 => fyi%cb0
         ycbuf1 => fyi%cb1
#else
         yfbuf1 => fyi%fb1
         ycbuf1 => fyi%cb1
#endif
         !*==============*
         !*  coarse2buf  * must not read from edges
         !*==============*
         do k = 1, fxi%snp0
            i = fxi%cndx0_l(k,1)
            j = fxi%cndx0_l(k,2)
#ifndef SINGLE_A2A
            xcbuf0(k) = fxc(i,j)*fxi%wt0_l(k)
#else
            ycbuf1(fxi%smap0(k)) = fxc(i,j)*fxi%wt0_l(k)
#endif
         end do
         do k = 1, fxi%snp1
            i = fxi%cndx1_l(k,1)
            j = fxi%cndx1_l(k,2)
#ifndef SINGLE_A2A
            xcbuf1(k) = fxc(i,j)*fxi%wt1_l(k)
#else
            ycbuf1(fxi%smap1(k)) = fxc(i,j)*fxi%wt1_l(k)
#endif
         end do
         do k = 1, fyi%snp0
            i = fyi%cndx0_l(k,1)
            j = fyi%cndx0_l(k,2)
#ifndef SINGLE_A2A
            ycbuf0(k) = fyc(i,j)*fyi%wt0_l(k)
#else
            ycbuf1(fyi%smap0(k)) = fyc(i,j)*fyi%wt0_l(k)
#endif
         end do
         do k = 1, fyi%snp1
            i = fyi%cndx1_l(k,1)
            j = fyi%cndx1_l(k,2)
#ifndef SINGLE_A2A
            ycbuf1(k) = fyc(i,j)*fyi%wt1_l(k)
#else
            ycbuf1(fyi%smap1(k)) = fyc(i,j)*fyi%wt1_l(k)
#endif
         end do
         !*==============*
         !*  alltoallv   *
         !*==============*
#ifndef A2A3D
#ifndef MULTI
#ifndef SINGLE_A2A
         call MPI_Alltoallv(xcbuf0, fxi%sendcnts0, fxi%sdispls0, REAL_MPI, &
                            xfbuf0, fxi%recvcnts0, fxi%rdispls0, REAL_MPI, MPI_COMM_WORLD, ierr)
         call MPI_Alltoallv(xcbuf1, fxi%sendcnts1, fxi%sdispls1, REAL_MPI, &
                            xfbuf1, fxi%recvcnts1, fxi%rdispls1, REAL_MPI, MPI_COMM_WORLD, ierr)
         call MPI_Alltoallv(ycbuf0, fyi%sendcnts0, fyi%sdispls0, REAL_MPI, &
                            yfbuf0, fyi%recvcnts0, fyi%rdispls0, REAL_MPI, MPI_COMM_WORLD, ierr)
#endif
         call MPI_Alltoallv(ycbuf1, fyi%sendcnts1, fyi%sdispls1, REAL_MPI, &
                            yfbuf1, fyi%recvcnts1, fyi%rdispls1, REAL_MPI, MPI_COMM_WORLD, ierr)
#else
#ifndef SINGLE_A2A
         call MPI_Alltoallv(xcbuf0, fxi%sendcnts0, fxi%sdispls0, REAL_MPI, &
                            xfbuf0, fxi%recvcnts0, fxi%rdispls0, REAL_MPI, MPI_MEMBER_WORLD, ierr)
         call MPI_Alltoallv(xcbuf1, fxi%sendcnts1, fxi%sdispls1, REAL_MPI, &
                            xfbuf1, fxi%recvcnts1, fxi%rdispls1, REAL_MPI, MPI_MEMBER_WORLD, ierr)
         call MPI_Alltoallv(ycbuf0, fyi%sendcnts0, fyi%sdispls0, REAL_MPI, &
                            yfbuf0, fyi%recvcnts0, fyi%rdispls0, REAL_MPI, MPI_MEMBER_WORLD, ierr)
#endif
         call MPI_Alltoallv(ycbuf1, fyi%sendcnts1, fyi%sdispls1, REAL_MPI, &
                            yfbuf1, fyi%recvcnts1, fyi%rdispls1, REAL_MPI, MPI_MEMBER_WORLD, ierr)
#endif
#else
#ifndef SINGLE_A2A
         call A2A3D_execute(fxi%snp0, xcbuf0, fxi%rnp0, xfbuf0, fxi%handler0)
         call A2A3D_execute(fxi%snp1, xcbuf1, fxi%rnp1, xfbuf1, fxi%handler1)
         call A2A3D_execute(fyi%snp0, ycbuf0, fyi%rnp0, yfbuf0, fyi%handler0)
         call A2A3D_execute(fyi%snp1, ycbuf1, fyi%rnp1, yfbuf1, fyi%handler1)
#else
         call A2A3D_execute(fxi%snp0+fxi%snp1+fyi%snp0+fyi%snp1, ycbuf1, &
                            fxi%rnp0+fxi%rnp1+fyi%rnp0+fyi%rnp1, yfbuf1, fyi%handler1)
#endif
#endif
         !*==============*
         !*  buf2fine    * write to edges
         !*==============*
         do k = 1, fxi%rnp0
            i = fxi%fndx0_l(k,1)
            j = fxi%fndx0_l(k,2)
#ifndef SINGLE_A2A
            fxf(i,j) = xfbuf0(k)
#else
            fxf(i,j) = yfbuf1(fxi%rmap0(k))
#endif
         end do
         do k = 1, fxi%rnp1
            i = fxi%fndx1_l(k,1)
            j = fxi%fndx1_l(k,2)
#ifndef SINGLE_A2A
            fxf(i,j) = fxf(i,j) + xfbuf1(k)
#else
            fxf(i,j) = fxf(i,j) + yfbuf1(fxi%rmap1(k))
#endif
         end do
         do k = 1, fyi%rnp0
            i = fyi%fndx0_l(k,1)
            j = fyi%fndx0_l(k,2)
#ifndef SINGLE_A2A
            fyf(i,j) = yfbuf0(k)
#else
            fyf(i,j) = yfbuf1(fyi%rmap0(k))
#endif
         end do
         do k = 1, fyi%rnp1
            i = fyi%fndx1_l(k,1)
            j = fyi%fndx1_l(k,2)
#ifndef SINGLE_A2A
            fyf(i,j) = fyf(i,j) + yfbuf1(k)
#else
            fyf(i,j) = fyf(i,j) + yfbuf1(fyi%rmap1(k))
#endif
         end do
! === USE_MPI_ALLTOALLV ========================================================
#endif
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
         lfac = 1.0d0/REAL_FUNC(fg%my%nr)

         nx = fg%my%nx
         ny = fg%my%ny
         px = fg%my%px
         py = fg%my%py
         rx = fg%my%rx
         ry = fg%my%ry
         north_rank = (ry-1)*px+rx
         south_rank = (ry+1)*px+rx
         east_rank = ry*px+rx+1
         west_rank = ry*px+rx-1
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) north_rank = MPI_PROC_NULL
         if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) south_rank = MPI_PROC_NULL
         if(iand(fg%my%has_boundary, EAST_BOUND)  /= 0) east_rank  = MPI_PROC_NULL
         if(iand(fg%my%has_boundary, WEST_BOUND)  /= 0) west_rank  = MPI_PROC_NULL

         ! Num. of shift indices to first element has a value.
         shift_st = mod(fg%my%nr-mod(fg%my%iy-1,fg%my%nr),fg%my%nr)
! === Upwind3 ==================================================================
         shift_st2 = mod(fg%my%nr-mod(fg%my%iy-1,fg%my%nr)+1,fg%my%nr)
! ==============================================================================
         ! Num. of shift indices to last element has a value.
         shift_en = mod(fg%my%iyend-1,fg%my%nr)
! === Upwind3 ==================================================================
!        shift_en2 = mod(fg%my%iyend,fg%my%nr)
         shift_en2 = mod(fg%my%iyend-2,fg%my%nr)
! ==============================================================================

         if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) then
            jst = 1        + shift_st
            jen = fg%my%ny - shift_en
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst = jst + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen = jen - 1
! === Upwind3 ==================================================================
            jst2 = 1        + shift_st2
            jen2 = fg%my%ny - shift_en2
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst2 = jst2 + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen2 = jen2 - 1
! ==============================================================================
#ifndef MULTI
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_n, 1, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, ireq11(1), ierr)
!           call MPI_Irecv(recv_from_s, 1, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, ireq21(1), ierr)
!           send_to_n = fxf(0,jst)
!           call MPI_Isend(send_to_n, 1, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, ireq11(2), ierr)
!           send_to_s = fxf(0,jen)
!           call MPI_Isend(send_to_s, 1, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, ireq21(2), ierr)
            call MPI_Irecv(recv_from_n, 3, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, ireq11(1), ierr)
            call MPI_Irecv(recv_from_s, 3, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, ireq21(1), ierr)
            send_to_n(1) = fxf(0, jst)
            send_to_n(2) = fxf(-2,jst)
            send_to_n(3) = fyf(-1,jst2)
            call MPI_Isend(send_to_n, 3, REAL_MPI, north_rank, 0, MPI_COMM_WORLD, ireq11(2), ierr)
            send_to_s(1) = fxf(0, jen)
            send_to_s(2) = fxf(-2,jen)
            send_to_s(3) = fyf(-1,jen2)
            call MPI_Isend(send_to_s, 3, REAL_MPI, south_rank, 0, MPI_COMM_WORLD, ireq21(2), ierr)
! ==============================================================================
#else
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_n, 1, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, ireq11(1), ierr)
!           call MPI_Irecv(recv_from_s, 1, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, ireq21(1), ierr)
!           send_to_n = fxf(0,jst)
!           call MPI_Isend(send_to_n, 1, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, ireq11(2), ierr)
!           send_to_s = fxf(0,jen)
!           call MPI_Isend(send_to_s, 1, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, ireq21(2), ierr)
            call MPI_Irecv(recv_from_n, 3, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, ireq11(1), ierr)
            call MPI_Irecv(recv_from_s, 3, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, ireq21(1), ierr)
            send_to_n(1) = fxf(0, jst)
            send_to_n(2) = fxf(-2,jst)
            send_to_n(3) = fyf(-1,jst2)
            call MPI_Isend(send_to_n, 3, REAL_MPI, north_rank, 0, MPI_MEMBER_WORLD, ireq11(2), ierr)
            send_to_s(1) = fxf(0, jen)
            send_to_s(2) = fxf(-2,jen)
            send_to_s(3) = fyf(-1,jen2)
            call MPI_Isend(send_to_s, 3, REAL_MPI, south_rank, 0, MPI_MEMBER_WORLD, ireq21(2), ierr)
! ==============================================================================
#endif
         end if

         if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) then
            jst = 1        + shift_st
            jen = fg%my%ny - shift_en
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst = jst + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen = jen - 1
! === Upwind3 ==================================================================
            jst2 = 1        + shift_st2
            jen2 = fg%my%ny - shift_en2
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst2 = jst2 + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen2 = jen2 - 1
! ==============================================================================
#ifndef MULTI
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_n, 1, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, ireq12(1), ierr)
!           call MPI_Irecv(recv_from_s, 1, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, ireq22(1), ierr)
!           send_to_n = fxf(fg%my%nx,jst)
!           call MPI_Isend(send_to_n, 1, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, ireq12(2), ierr)
!           send_to_s = fxf(fg%my%nx,jen)
!           call MPI_Isend(send_to_s, 1, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, ireq22(2), ierr)
            call MPI_Irecv(recv_from_n, 2, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, ireq12(1), ierr)
            call MPI_Irecv(recv_from_s, 2, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, ireq22(1), ierr)
            send_to_n(1) = fxf(fg%my%nx,  jst)
            send_to_n(2) = fyf(fg%my%nx+1,jst2)
            call MPI_Isend(send_to_n, 2, REAL_MPI, north_rank, 1, MPI_COMM_WORLD, ireq12(2), ierr)
            send_to_s(1) = fxf(fg%my%nx,  jen)
            send_to_s(2) = fyf(fg%my%nx+1,jen2)
            call MPI_Isend(send_to_s, 2, REAL_MPI, south_rank, 1, MPI_COMM_WORLD, ireq22(2), ierr)
! ==============================================================================
#else
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_n, 1, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, ireq12(1), ierr)
!           call MPI_Irecv(recv_from_s, 1, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, ireq22(1), ierr)
!           send_to_n = fxf(fg%my%nx,jst)
!           call MPI_Isend(send_to_n, 1, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, ireq12(2), ierr)
!           send_to_s = fxf(fg%my%nx,jen)
!           call MPI_Isend(send_to_s, 1, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, ireq22(2), ierr)
            call MPI_Irecv(recv_from_n, 2, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, ireq12(1), ierr)
            call MPI_Irecv(recv_from_s, 2, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, ireq22(1), ierr)
            send_to_n(1) = fxf(fg%my%nx,  jst)
            send_to_n(2) = fyf(fg%my%nx+1,jst2)
            call MPI_Isend(send_to_n, 2, REAL_MPI, north_rank, 1, MPI_MEMBER_WORLD, ireq12(2), ierr)
            send_to_s(1) = fxf(fg%my%nx,  jen)
            send_to_s(2) = fyf(fg%my%nx+1,jen2)
            call MPI_Isend(send_to_s, 2, REAL_MPI, south_rank, 1, MPI_MEMBER_WORLD, ireq22(2), ierr)
! ==============================================================================
#endif
         end if
         
         ! Num. of shift indices to first element has a value.
         shift_st = mod(fg%my%nr-mod(fg%my%ix-1,fg%my%nr),fg%my%nr)
! === Upwind3 ==================================================================
         shift_st2 = mod(fg%my%nr-mod(fg%my%ix-1,fg%my%nr)+1,fg%my%nr)
! ==============================================================================
         ! Num. of shift indices to last element has a value.
         shift_en = mod(fg%my%ixend-1,fg%my%nr)
! === Upwind3 ==================================================================
!        shift_en2 = mod(fg%my%ixend,fg%my%nr)
         shift_en2 = mod(fg%my%ixend-2,fg%my%nr)
! ==============================================================================

         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) then
            ist = 1        + shift_st
            ien = fg%my%nx - shift_en
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist = ist + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien = ien - 1
! === Upwind3 ==================================================================
            ist2 = 1        + shift_st2
            ien2 = fg%my%nx - shift_en2
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist2 = ist2 + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien2 = ien2 - 1
! ==============================================================================
#ifndef MULTI
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_w, 1, REAL_MPI, west_rank, 2, MPI_COMM_WORLD, ireq13(1), ierr)
!           call MPI_Irecv(recv_from_e, 1, REAL_MPI, east_rank, 2, MPI_COMM_WORLD, ireq23(1), ierr)
!           send_to_w = fyf(ist,0)
!           call MPI_Isend(send_to_w, 1, REAL_MPI, west_rank, 2, MPI_COMM_WORLD, ireq13(2), ierr)
!           send_to_e = fyf(ien,0)
!           call MPI_Isend(send_to_e, 1, REAL_MPI, east_rank, 2, MPI_COMM_WORLD, ireq23(2), ierr)
            call MPI_Irecv(recv_from_w, 3, REAL_MPI, west_rank, 2, MPI_COMM_WORLD, ireq13(1), ierr)
            call MPI_Irecv(recv_from_e, 3, REAL_MPI, east_rank, 2, MPI_COMM_WORLD, ireq23(1), ierr)
            send_to_w(1) = fyf(ist, 0)
            send_to_w(2) = fxf(ist2,-1)
            send_to_w(3) = fyf(ist,-2)
            call MPI_Isend(send_to_w, 3, REAL_MPI, west_rank, 2, MPI_COMM_WORLD, ireq13(2), ierr)
            send_to_e(1) = fyf(ien, 0)
            send_to_e(2) = fxf(ien2,-1)
            send_to_e(3) = fyf(ien,-2)
            call MPI_Isend(send_to_e, 3, REAL_MPI, east_rank, 2, MPI_COMM_WORLD, ireq23(2), ierr)
! ==============================================================================
#else
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_w, 1, REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, ireq13(1), ierr)
!           call MPI_Irecv(recv_from_e, 1, REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, ireq23(1), ierr)
!           send_to_w = fyf(ist,0)
!           call MPI_Isend(send_to_w, 1, REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, ireq13(2), ierr)
!           send_to_e = fyf(ien,0)
!           call MPI_Isend(send_to_e, 1, REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, ireq23(2), ierr)
            call MPI_Irecv(recv_from_w, 3, REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, ireq13(1), ierr)
            call MPI_Irecv(recv_from_e, 3, REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, ireq23(1), ierr)
            send_to_w(1) = fyf(ist, 0)
            send_to_w(2) = fxf(ist2,-1)
            send_to_w(3) = fyf(ist,-2)
            call MPI_Isend(send_to_w, 3, REAL_MPI, west_rank, 2, MPI_MEMBER_WORLD, ireq13(2), ierr)
            send_to_e(1) = fyf(ien, 0)
            send_to_e(2) = fxf(ien2,-1)
            send_to_e(3) = fyf(ien,-2)
            call MPI_Isend(send_to_e, 3, REAL_MPI, east_rank, 2, MPI_MEMBER_WORLD, ireq23(2), ierr)
! ==============================================================================
#endif
         end if

         if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) then
            ist = 1        + shift_st
            ien = fg%my%nx - shift_en
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist = ist + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien = ien - 1
! === Upwind3 ==================================================================
            ist2 = 1        + shift_st2
            ien2 = fg%my%nx - shift_en2
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist2 = ist2 + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien2 = ien2 - 1
! ==============================================================================
#ifndef MULTI
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_w, 1, REAL_MPI, west_rank, 3, MPI_COMM_WORLD, ireq14(1), ierr)
!           call MPI_Irecv(recv_from_e, 1, REAL_MPI, east_rank, 3, MPI_COMM_WORLD, ireq24(1), ierr)
!           send_to_w = fyf(ist,fg%my%ny)
!           call MPI_Isend(send_to_w, 1, REAL_MPI, west_rank, 3, MPI_COMM_WORLD, ireq14(2), ierr)
!           send_to_e = fyf(ien,fg%my%ny)
!           call MPI_Isend(send_to_e, 1, REAL_MPI, east_rank, 3, MPI_COMM_WORLD, ireq24(2), ierr)
            call MPI_Irecv(recv_from_w, 2, REAL_MPI, west_rank, 3, MPI_COMM_WORLD, ireq14(1), ierr)
            call MPI_Irecv(recv_from_e, 2, REAL_MPI, east_rank, 3, MPI_COMM_WORLD, ireq24(1), ierr)
            send_to_w(1) = fyf(ist,  fg%my%ny)
            send_to_w(2) = fxf(ist2,fg%my%ny+1)
            call MPI_Isend(send_to_w, 2, REAL_MPI, west_rank, 3, MPI_COMM_WORLD, ireq14(2), ierr)
            send_to_e(1) = fyf(ien,  fg%my%ny)
            send_to_e(2) = fxf(ien2,fg%my%ny+1)
            call MPI_Isend(send_to_e, 2, REAL_MPI, east_rank, 3, MPI_COMM_WORLD, ireq24(2), ierr)
! ==============================================================================
#else
! === Upwind3 ==================================================================
!           call MPI_Irecv(recv_from_w, 1, REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, ireq14(1), ierr)
!           call MPI_Irecv(recv_from_e, 1, REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, ireq24(1), ierr)
!           send_to_w = fyf(ist,fg%my%ny)
!           call MPI_Isend(send_to_w, 1, REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, ireq14(2), ierr)
!           send_to_e = fyf(ien,fg%my%ny)
!           call MPI_Isend(send_to_e, 1, REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, ireq24(2), ierr)
            call MPI_Irecv(recv_from_w, 2, REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, ireq14(1), ierr)
            call MPI_Irecv(recv_from_e, 2, REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, ireq24(1), ierr)
            send_to_w(1) = fyf(ist,  fg%my%ny)
            send_to_w(2) = fxf(ist2,fg%my%ny+1)
            call MPI_Isend(send_to_w, 2, REAL_MPI, west_rank, 3, MPI_MEMBER_WORLD, ireq14(2), ierr)
            send_to_e(1) = fyf(ien,  fg%my%ny)
            send_to_e(2) = fxf(ien2,fg%my%ny+1)
            call MPI_Isend(send_to_e, 2, REAL_MPI, east_rank, 3, MPI_MEMBER_WORLD, ireq24(2), ierr)
! ==============================================================================
#endif
         end if

         ! Num. of shift indices to first element has a value.
         shift_st = mod(fg%my%nr-mod(fg%my%iy-1,fg%my%nr),fg%my%nr)
! === Upwind3 ==================================================================
         shift_st2 = mod(fg%my%nr-mod(fg%my%iy-1,fg%my%nr)+1,fg%my%nr)
! ==============================================================================
         ! Num. of shift indices to last element has a value.
         shift_en = mod(fg%my%iyend-1,fg%my%nr)
! === Upwind3 ==================================================================
!        shift_en2 = mod(fg%my%iyend,fg%my%nr)
         shift_en2 = mod(fg%my%iyend-2,fg%my%nr)
! ==============================================================================

         if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) then
            jst = 1        + shift_st
            jen = fg%my%ny - shift_en
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst = jst + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen = jen - 1
! === Upwind3 ==================================================================
            jst2 = 1        + shift_st2
            jen2 = fg%my%ny - shift_en2
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst2 = jst2 + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen2 = jen2 - 1
! ==============================================================================

            ! Calc. others.
            do j = jst, jen-1
               iy = j + fg%my%ky - 1
               imod = mod(iy-1,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fxf(0,       j) = fxf(0,       ind0)*fac0 + fxf(0,       ind1)*fac1
! === Upwind3 ==================================================================
               fxf(-2,      j) = fxf(-2,      ind0)*fac0 + fxf(-2,      ind1)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            do j = jst2, jen2-1
               iy = j + fg%my%ky - 1
               imod = mod(iy-2,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fyf(-1,      j) = fyf(-1,      ind0)*fac0 + fyf(-1,      ind1)*fac1
            end do
! ==============================================================================

            ! Calc. north boundary elemnts.
            call MPI_Waitall(2, ireq11, stat11, ierr)
            do j = jst-shift_st, jst-1
               iy = j + fg%my%ky - 1
               imod = mod(iy-1,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fxf(0,       j) = recv_from_n*fac0 + fxf(0,       ind1)*fac1
               fxf(0,       j) = recv_from_n(1)*fac0 + fxf(0,       ind1)*fac1
               fxf(-2,      j) = recv_from_n(2)*fac0 + fxf(-2,      ind1)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) then
               do j = jst2-shift_st2, jst2-1
                  iy = j + fg%my%ky - 1
                  imod = mod(iy-2,fg%my%nr)
                  ind0 = j - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fyf(-1,      j) = recv_from_n(3)*fac0 + fyf(-1,      ind1)*fac1
               end do
            end if
! ==============================================================================
            ! Calc. south boundary elemnts.
            call MPI_Waitall(2, ireq21, stat21, ierr)
            do j = jen+1, jen+shift_en
               iy = j + fg%my%ky - 1
               imod = mod(iy-1,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fxf(0,       j) = fxf(0,       ind0)*fac0 + recv_from_s*fac1
               fxf(0,       j) = fxf(0,       ind0)*fac0 + recv_from_s(1)*fac1
               fxf(-2,      j) = fxf(-2,      ind0)*fac0 + recv_from_s(2)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) then
               do j = jen2+1, jen2+shift_en2
                  iy = j + fg%my%ky - 1
                  imod = mod(iy-2,fg%my%nr)
                  ind0 = j - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fyf(-1,      j) = fyf(-1,      ind0)*fac0 + recv_from_s(3)*fac1
               end do
            end if
! ==============================================================================
         end if

         if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) then
            jst = 1        + shift_st
            jen = fg%my%ny - shift_en
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst = jst + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen = jen - 1
! === Upwind3 ==================================================================
            jst2 = 1        + shift_st2
            jen2 = fg%my%ny - shift_en2
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) jst2 = jst2 + 1
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) jen2 = jen2 - 1
! ==============================================================================

            ! Calc. others.
            do j = jst, jen-1
               iy = j + fg%my%ky - 1
               imod = mod(iy-1,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fxf(fg%my%nx,j) = fxf(fg%my%nx,ind0)*fac0 + fxf(fg%my%nx,ind1)*fac1
            end do
! === Upwind3 ==================================================================
            do j = jst2, jen2-1
               iy = j + fg%my%ky - 1
               imod = mod(iy-2,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fyf(fg%my%nx+1,j) = fyf(fg%my%nx+1,ind0)*fac0 + fyf(fg%my%nx+1,ind1)*fac1
            end do
! ==============================================================================

            ! Calc. north boundary elemnts.
            call MPI_Waitall(2, ireq12, stat12, ierr)
            do j = jst-shift_st, jst-1
               iy = j + fg%my%ky - 1
               imod = mod(iy-1,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fxf(fg%my%nx,j) = recv_from_n*fac0 + fxf(fg%my%nx,ind1)*fac1
               fxf(fg%my%nx,  j) = recv_from_n(1)*fac0 + fxf(fg%my%nx,  ind1)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, NORTH_BOUND) == 0) then
               do j = jst2-shift_st2, jst2-1
                  iy = j + fg%my%ky - 1
                  imod = mod(iy-2,fg%my%nr)
                  ind0 = j - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fyf(fg%my%nx+1,j) = recv_from_n(2)*fac0 + fyf(fg%my%nx+1,ind1)*fac1
               end do
            end if
! ==============================================================================
            ! Calc. south boundary elemnts.
            call MPI_Waitall(2, ireq22, stat22, ierr)
            do j = jen+1, jen+shift_en
               iy = j + fg%my%ky - 1
               imod = mod(iy-1,fg%my%nr)
               ind0 = j - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fxf(fg%my%nx,j) = fxf(fg%my%nx,ind0)*fac0 + recv_from_s*fac1
               fxf(fg%my%nx,  j) = fxf(fg%my%nx,  ind0)*fac0 + recv_from_s(1)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, SOUTH_BOUND) == 0) then
               do j = jen2+1, jen2+shift_en2
                  iy = j + fg%my%ky - 1
                  imod = mod(iy-2,fg%my%nr)
                  ind0 = j - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fyf(fg%my%nx+1,j) = fyf(fg%my%nx+1,ind0)*fac0 + recv_from_s(2)*fac1
               end do
            end if
! ==============================================================================
         end if
         
         ! Num. of shift indices to first element has a value.
         shift_st = mod(fg%my%nr-mod(fg%my%ix-1,fg%my%nr),fg%my%nr)
! === Upwind3 ==================================================================
         shift_st2 = mod(fg%my%nr-mod(fg%my%ix-1,fg%my%nr)+1,fg%my%nr)
! ==============================================================================
         ! Num. of shift indices to last element has a value.
         shift_en = mod(fg%my%ixend-1,fg%my%nr)
! === Upwind3 ==================================================================
!        shift_en2 = mod(fg%my%ixend,fg%my%nr)
         shift_en2 = mod(fg%my%ixend-2,fg%my%nr)
! ==============================================================================

         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) then
            ist = 1        + shift_st
            ien = fg%my%nx - shift_en
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist = ist + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien = ien - 1
! === Upwind3 ==================================================================
            ist2 = 1        + shift_st2
            ien2 = fg%my%nx - shift_en2
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist2 = ist2 + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien2 = ien2 - 1
! ==============================================================================

            ! Calc. others.
            do i = ist, ien-1
               ix = i + fg%my%kx - 1
               imod = mod(ix-1,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fyf(i,0)        = fyf(ind0,0       )*fac0 + fyf(ind1,0       )*fac1
! === Upwind3 ==================================================================
               fyf(i,-2)       = fyf(ind0,-2      )*fac0 + fyf(ind1,-2      )*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            do i = ist2, ien2-1
               ix = i + fg%my%kx - 1
               imod = mod(ix-2,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fxf(i,-1)       = fxf(ind0,-1      )*fac0 + fxf(ind1,-1      )*fac1
            end do
! ==============================================================================

            ! Calc. west boundary elemnts.
            call MPI_Waitall(2, ireq13, stat13, ierr)
            do i = ist-shift_st, ist-1
               ix = i + fg%my%kx - 1
               imod = mod(ix-1,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fyf(i,0)        = recv_from_w*fac0 + fyf(ind1,0       )*fac1
               fyf(i,0)        = recv_from_w(1)*fac0 + fyf(ind1,0       )*fac1
               fyf(i,-2)       = recv_from_w(3)*fac0 + fyf(ind1,-2      )*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) then
               do i = ist2-shift_st2, ist2-1
                  ix = i + fg%my%kx - 1
                  imod = mod(ix-2,fg%my%nr)
                  ind0 = i - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fxf(i,-1)       = recv_from_w(2)*fac0 + fxf(ind1,-1      )*fac1
               end do
            end if
! ==============================================================================
            ! Calc. east boundary elemnts.
            call MPI_Waitall(2, ireq23, stat23, ierr)
            do i = ien+1, ien+shift_en
               ix = i + fg%my%kx - 1
               imod = mod(ix-1,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fyf(i,0)        = fyf(ind0,0       )*fac0 + recv_from_e*fac1
               fyf(i,0)        = fyf(ind0,0       )*fac0 + recv_from_e(1)*fac1
               fyf(i,-2)       = fyf(ind0,-2      )*fac0 + recv_from_e(3)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) then
               do i = ien2+1, ien2+shift_en2
                  ix = i + fg%my%kx - 1
                  imod = mod(ix-2,fg%my%nr)
                  ind0 = i - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fxf(i,-1)       = fxf(ind0,-1      )*fac0 + recv_from_e(2)*fac1
               end do
            end if
! ==============================================================================
         end if

         if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) then
            ist = 1        + shift_st
            ien = fg%my%nx - shift_en
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist = ist + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien = ien - 1
! === Upwind3 ==================================================================
            ist2 = 1        + shift_st2
            ien2 = fg%my%nx - shift_en2
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) ist2 = ist2 + 1
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) ien2 = ien2 - 1
! ==============================================================================

            ! Calc. others.
            do i = ist, ien-1
               ix = i + fg%my%kx - 1
               imod = mod(ix-1,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fyf(i,fg%my%ny) = fyf(ind0,fg%my%ny)*fac0 + fyf(ind1,fg%my%ny)*fac1
            end do
! === Upwind3 ==================================================================
            do i = ist2, ien2-1
               ix = i + fg%my%kx - 1
               imod = mod(ix-2,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
               fxf(i,fg%my%ny+1) = fxf(ind0,fg%my%ny+1)*fac0 + fxf(ind1,fg%my%ny+1)*fac1
            end do
! ==============================================================================

            ! Calc. west boundary elemnts.
            call MPI_Waitall(2, ireq14, stat14, ierr)
            do i = ist-shift_st, ist-1
               ix = i + fg%my%kx - 1
               imod = mod(ix-1,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fyf(i,fg%my%ny) = recv_from_w*fac0 + fyf(ind1,fg%my%ny)*fac1
               fyf(i  ,fg%my%ny) = recv_from_w(1)*fac0 + fyf(  ind1,fg%my%ny)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, WEST_BOUND) == 0) then
               do i = ist2-shift_st2, ist2-1
                  ix = i + fg%my%kx - 1
                  imod = mod(ix-2,fg%my%nr)
                  ind0 = i - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fxf(i,fg%my%ny+1) = recv_from_w(2)*fac0 + fxf(ind1,fg%my%ny+1)*fac1
               end do
            end if
! ==============================================================================
            ! Calc. east boundary elemnts.
            call MPI_Waitall(2, ireq24, stat24, ierr)
            do i = ien+1, ien+shift_en
               ix = i + fg%my%kx - 1
               imod = mod(ix-1,fg%my%nr)
               ind0 = i - imod
               ind1 = ind0 + fg%my%nr
               fac0 = 1.0d0 - imod*lfac
               fac1 = 1.0d0 - fac0
! === Upwind3 ==================================================================
!              fyf(i,fg%my%ny) = fyf(ind0,fg%my%ny)*fac0 + recv_from_e*fac1
               fyf(i,  fg%my%ny) = fyf(ind0,  fg%my%ny)*fac0 + recv_from_e(1)*fac1
! ==============================================================================
            end do
! === Upwind3 ==================================================================
            if(iand(fg%my%has_boundary, EAST_BOUND) == 0) then
               do i = ien2+1, ien2+shift_en2
                  ix = i + fg%my%kx - 1
                  imod = mod(ix-2,fg%my%nr)
                  ind0 = i - imod
                  ind1 = ind0 + fg%my%nr
                  fac0 = 1.0d0 - imod*lfac
                  fac1 = 1.0d0 - fac0
                  fxf(i,fg%my%ny+1) = fxf(ind0,fg%my%ny+1)*fac0 + recv_from_e(2)*fac1
               end do
            end if
! ==============================================================================
         end if
! ==============================================================================
         t0 = 1.0d0/REAL_FUNC(fg%my%nr)
         t1 = 1.0d0 - t0

         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) then
            ist = 1
            if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) ist = 0
            do i = ist, fg%my%nx
               fxf(i,0) = fxf(i,0)*t0 + fxf(i,1)*t1
            end do
         end if

         if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) then
            jst = 1
            if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) jst = 0
            do j = jst, fg%my%ny
               fyf(0,j) = fyf(0,j)*t0 + fyf(1,j)*t1
            end do
         end if
! === Upwind3 ==================================================================
         nx = fg%my%nx
         ny = fg%my%ny

         fac0 = lfac
         fac1 = 1.0d0 - fac0

         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) then
            if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) then
               fxf(1,  -1) = fxf(1, -1)*fac0 + fxf(2, -1)*fac1
               fyf(-1,  1) = fyf(-1, 1)*fac0 + fyf(-1, 2)*fac1
            end if
            if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) then
               fxf(nx-1, -1) = fxf(nx-2, -1)*fac1 + fxf(nx,  -1)*fac0
               fxf(nx,   -1) = fxf(nx-2, -1)*fac0 + fxf(nx,  -1)*fac1
               fyf(nx+1,  1) = fyf(nx+1,  1)*fac0 + fyf(nx+1, 2)*fac1
            end if
         end if

         if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) then
            if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) then
               fxf(1,    ny+1) = fxf(1,   ny+1)*fac0 + fxf(2, ny+1)*fac1
               fyf(-1,   ny-1) = fyf(-1,  ny-2)*fac1 + fyf(-1,  ny)*fac0
               fyf(-1,     ny) = fyf(-1,  ny-2)*fac0 + fyf(-1,  ny)*fac1
            end if
            if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) then
               fxf(nx-1, ny+1) = fxf(nx-2,ny+1)*fac1 + fxf(nx,ny+1)*fac0
               fxf(nx,   ny+1) = fxf(nx-2,ny+1)*fac0 + fxf(nx,ny+1)*fac1
               fyf(nx+1, ny-1) = fyf(nx+1,ny-2)*fac1 + fyf(nx+1,ny)*fac0
               fyf(nx+1,   ny) = fyf(nx+1,ny-2)*fac0 + fyf(nx+1,ny)*fac1
            end if
         end if
! ==============================================================================
#endif
      else if(mode == HGT) then
#ifndef MPI
!$omp parallel
!$omp do
         do k = 1, hzi%np
            hzf(hzi%fndx(k,1),hzi%fndx(k,2)) = &
               hzc(hzi%cndx0(k,1),hzi%cndx0(k,2))*hzi%wt0(k) + &
               hzc(hzi%cndx1(k,1),hzi%cndx1(k,2))*hzi%wt1(k)
         end do
!$omp end parallel
         t0 = 1.0d0/REAL_FUNC(fg%my%nr)
         t1 = 1.0d0 - t0
         do i = 1, fg%my%nx
! === Upwind3 ==================================================================
            hzf(i,-1        ) = hzf(i,0         )*t1 + hzf(i,1       )*t0
            hzf(i,fg%my%ny+2) = hzf(i,fg%my%ny+1)*t1 + hzf(i,fg%my%ny)*t0
! ==============================================================================
            hzf(i,0         ) = hzf(i,0         )*t0 + hzf(i,1       )*t1
            hzf(i,fg%my%ny+1) = hzf(i,fg%my%ny+1)*t0 + hzf(i,fg%my%ny)*t1
         end do
         do j = 0, fg%my%ny+1
! === Upwind3 ==================================================================
            hzf(-1,        j) = hzf(0,         j)*t1 + hzf(1,       j)*t0
            hzf(fg%my%nx+2,j) = hzf(fg%my%nx+1,j)*t1 + hzf(fg%my%nx,j)*t0
! ==============================================================================
            hzf(0,         j) = hzf(0,         j)*t0 + hzf(1,       j)*t1
            hzf(fg%my%nx+1,j) = hzf(fg%my%nx+1,j)*t0 + hzf(fg%my%nx,j)*t1
         end do
#else
#ifndef USE_ALLTOALLV
!$omp parallel
!$omp single
         zfbuf => hzi%fb
         zcbuf => hzi%cb
!$omp end single
!$omp do
         do k = 1, hzi%np
            zfbuf(k) = 0.0d0
            zcbuf(k) = 0.0d0
         end do
         !*==============*
         !*  coarse2buf  * must not read from edges
         !*==============*
!$omp do private(ix, iy, i, j)
         do k = 1, hzi%np
            ix = hzi%cndx0(k,1)
            iy = hzi%cndx0(k,2)
            if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
               i = ix - cg%my%kx + 1
               j = iy - cg%my%ky + 1
               zcbuf(k) = zcbuf(k) + hzc(i,j)*hzi%wt0(k)
            end if
         end do
!$omp do private(ix, iy, i, j)
         do k = 1, hzi%np
            ix = hzi%cndx1(k,1)
            iy = hzi%cndx1(k,2)
            if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
               i = ix - cg%my%kx + 1 ! with edges
               j = iy - cg%my%ky + 1 ! with edges
               zcbuf(k) = zcbuf(k) + hzc(i,j)*hzi%wt1(k)
            end if
         end do
         !*==============*
         !*  allreduce   *
         !*==============*
!$omp single
#ifndef MULTI
         call MPI_Allreduce(zcbuf, zfbuf, hzi%np, REAL_MPI, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
         call MPI_Allreduce(zcbuf, zfbuf, hzi%np, REAL_MPI, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
         if(ierr /= 0) then
            select case (ierr)
               case(MPI_ERR_BUFFER)
                  write(0,'(a)') 'MPI Error : Invalid buffer pointer'
               case(MPI_ERR_COUNT)
                  write(0,'(a)') 'MPI Error : Invalid count argument'
               case(MPI_ERR_TYPE)
                  write(0,'(a)') 'MPI Error : Invalid datatype argument'
               case(MPI_ERR_OP)
                  write(0,'(a)') 'MPI Error : Invalid operation'
               case(MPI_ERR_COMM)
                  write(0,'(a)') 'MPI Error : Invalid communicator'
               case default
                  write(0,'(a)') 'MPI Error : Unknown error'
            end select
            call fatal_error(ierr)
         end if

         ist = fg%my%kx
! === Upwind3 ==================================================================
!        if(iand(fg%my%has_boundary, WEST_BOUND)  /= 0) ist = 0
         if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) ist = -1
         ien = fg%my%kxend+1
         if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) ien = fg%my%kxend+2
! ==============================================================================
         jst = fg%my%ky
! === Upwind3 ==================================================================
!        if(iand(fg%my%has_boundary, NORTH_BOUND)  /= 0) jst = 0
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) jst = -1
         jen = fg%my%kyend+1
         if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) jen = fg%my%kyend+2
! ==============================================================================
!$omp end single
         !*==============*
         !*  buf2fine    * write to edges
         !*==============*
!$omp do private(ix, iy, i, j)
         do k = 1, hzi%np
            ix = hzi%fndx(k,1)
            iy = hzi%fndx(k,2)
! === Upwind3 ==================================================================
!           if(ist <= ix .and. ix <= fg%my%kxend+1 .and. jst <= iy .and. iy <= fg%my%kyend+1) then
            if(ist <= ix .and. ix <= ien .and. jst <= iy .and. iy <= jen) then
! ==============================================================================
               i = ix - fg%my%kx + 1 ! with edges
               j = iy - fg%my%ky + 1 ! with edges
               hzf(i,j) = zfbuf(k)
            end if
         end do
!$omp end parallel
#else
! === USE_MPI_ALLTOALLV ========================================================
#ifndef SINGLE_A2A
         zfbuf0 => hzi%fb0
         zfbuf1 => hzi%fb1
         zcbuf0 => hzi%cb0
         zcbuf1 => hzi%cb1
#else
         zfbuf1 => hzi%fb1
         zcbuf1 => hzi%cb1
#endif
         !*==============*
         !*  coarse2buf  * must not read from edges
         !*==============*
         do k = 1, hzi%snp0
            i = hzi%cndx0_l(k,1)
            j = hzi%cndx0_l(k,2)
#ifndef SINGLE_A2A
            zcbuf0(k) = hzc(i,j)*hzi%wt0_l(k)
#else
            zcbuf1(hzi%smap0(k)) = hzc(i,j)*hzi%wt0_l(k)
#endif
         end do
         do k = 1, hzi%snp1
            i = hzi%cndx1_l(k,1)
            j = hzi%cndx1_l(k,2)
#ifndef SINGLE_A2A
            zcbuf1(k) = hzc(i,j)*hzi%wt1_l(k)
#else
            zcbuf1(hzi%smap1(k)) = hzc(i,j)*hzi%wt1_l(k)
#endif
         end do
         !*==============*
         !*  alltoallv   *
         !*==============*
#ifndef A2A3D
#ifndef MULTI
#ifndef SINGLE_A2A
         call MPI_Alltoallv(zcbuf0, hzi%sendcnts0, hzi%sdispls0, REAL_MPI, &
                            zfbuf0, hzi%recvcnts0, hzi%rdispls0, REAL_MPI, MPI_COMM_WORLD, ierr)
#endif
         call MPI_Alltoallv(zcbuf1, hzi%sendcnts1, hzi%sdispls1, REAL_MPI, &
                            zfbuf1, hzi%recvcnts1, hzi%rdispls1, REAL_MPI, MPI_COMM_WORLD, ierr)
#else
#ifndef SINGLE_A2A
         call MPI_Alltoallv(zcbuf0, hzi%sendcnts0, hzi%sdispls0, REAL_MPI, &
                            zfbuf0, hzi%recvcnts0, hzi%rdispls0, REAL_MPI, MPI_MEMBER_WORLD, ierr)
#endif
         call MPI_Alltoallv(zcbuf1, hzi%sendcnts1, hzi%sdispls1, REAL_MPI, &
                            zfbuf1, hzi%recvcnts1, hzi%rdispls1, REAL_MPI, MPI_MEMBER_WORLD, ierr)
#endif
#else
#ifndef SINGLE_A2A
         call A2A3D_execute(hzi%snp0, zcbuf0, hzi%rnp0, zfbuf0, hzi%handler0)
         call A2A3D_execute(hzi%snp1, zcbuf1, hzi%rnp1, zfbuf1, hzi%handler1)
#else
         call A2A3D_execute(hzi%snp0+hzi%snp1, zcbuf1, hzi%rnp0+hzi%rnp1, zfbuf1, hzi%handler1)
#endif
#endif
         !*==============*
         !*  buf2fine    * write to edges
         !*==============*
         do k = 1, hzi%rnp0
            i = hzi%fndx0_l(k,1)
            j = hzi%fndx0_l(k,2)
#ifndef SINGLE_A2A
            hzf(i,j) = zfbuf0(k)
#else
            hzf(i,j) = zfbuf1(hzi%rmap0(k))
#endif
         end do
         do k = 1, hzi%rnp1
            i = hzi%fndx1_l(k,1)
            j = hzi%fndx1_l(k,2)
#ifndef SINGLE_A2A
            hzf(i,j) = hzf(i,j) + zfbuf1(k)
#else
            hzf(i,j) = hzf(i,j) + zfbuf1(hzi%rmap1(k))
#endif
         end do
! === USE_MPI_ALLTOALLV ========================================================
#endif

      t0 = 1.0d0/REAL_FUNC(fg%my%nr)
      t1 = 1.0d0 - t0

      if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) then
!$omp single
         ist = 1
!$omp end single
!$omp do
         do i = ist, fg%my%nx
! === Upwind3 ==================================================================
            hzf(i,-1) = hzf(i,0)*t1 + hzf(i,1)*t0
! ==============================================================================
            hzf(i,0) = hzf(i,0)*t0 + hzf(i,1)*t1
         end do
      end if

      if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) then
!$omp single
         ist = 1
!$omp end single
!$omp do
         do i = ist, fg%my%nx
! === Upwind3 ==================================================================
            hzf(i,fg%my%ny+2) = hzf(i,fg%my%ny+1)*t1 + hzf(i,fg%my%ny)*t0
! ==============================================================================
            hzf(i,fg%my%ny+1) = hzf(i,fg%my%ny+1)*t0 + hzf(i,fg%my%ny)*t1
         end do
      end if

      if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) then
!$omp single
         jst = 1
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) jst = 0
!$omp end single
!$omp do
         do j = jst, fg%my%ny+1
! === Upwind3 ==================================================================
            hzf(-1,j) = hzf(0,j)*t1 + hzf(1,j)*t0
! ==============================================================================
            hzf(0,j) = hzf(0,j)*t0 + hzf(1,j)*t1
         end do
      end if

      if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) then
!$omp single
         jst = 1
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) jst = 0
!$omp end single
!$omp do
         do j = jst, fg%my%ny+1
! === Upwind3 ==================================================================
            hzf(fg%my%nx+2,j) = hzf(fg%my%nx+1,j)*t1 + hzf(fg%my%nx,j)*t0
! ==============================================================================
            hzf(fg%my%nx+1,j) = hzf(fg%my%nx+1,j)*t0 + hzf(fg%my%nx,j)*t1
         end do
      end if
#endif
      end if

      return
   end subroutine interp2fine

#ifdef MPI
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
   subroutine make_a2a_params_c2c(dg, fxo, nprocs)
#ifndef SINGLE_A2A
#ifdef A2A3D
      use mod_a2a3d
#endif
#endif
      type(data_grids), target, intent(inout) :: dg
      type(interp_info), target, intent(inout) :: fxo
      integer(kind=4), intent(in) :: nprocs

      ! Pointers to type objecs.
      type(grid_info), pointer :: parent, my
      integer(kind=4), pointer :: snp0, rnp0
      integer(kind=4), pointer, dimension(:,:) :: fndx0_l, cndx0_l
      integer(kind=4), pointer, dimension(:) :: sendcnts0, sdispls0, recvcnts0, rdispls0

      ! Work vals.
      integer(kind=4) :: k, ix, iy, p, cnt, ind, ierr
      integer(kind=4) :: scnt, rcnt, total_scnt, total_rcnt

      ! Work arrays.
      integer(kind=4), allocatable, dimension(:) :: scnts, rcnts, displs
      integer(kind=4), allocatable, dimension(:) :: table, table_all

      parent => dg%parent
      my => dg%my

      snp0 => fxo%snp0
      rnp0 => fxo%rnp0

      allocate(scnts(0:nprocs-1))
      allocate(rcnts(0:nprocs-1))
      allocate(displs(0:nprocs-1))

      ! Count send elements.
      scnt = 0
      do k = 1, fxo%np
         ix = fxo%fndx(k,1)
         iy = fxo%fndx(k,2)
         if(my%ix <= ix .and. ix <= my%ixend .and. my%iy <= iy .and. iy <= my%iyend) then
            scnt = scnt + 1
         end if
      end do
#ifndef MULTI
      call MPI_Allgather(scnt, 1, MPI_INTEGER, scnts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(scnt, 1, MPI_INTEGER, scnts, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      ! Count recv elements.
      rcnt = 0
      do k = 1, fxo%np
         ix = fxo%cndx0(k,1)
         iy = fxo%cndx0(k,2)
         if(parent%kx <= ix .and. ix <= parent%kxend .and. parent%ky <= iy .and. iy <= parent%kyend) then
            rcnt = rcnt + 1
         end if
      end do
#ifndef MULTI
      call MPI_Allgather(rcnt, 1, MPI_INTEGER, rcnts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(rcnt, 1, MPI_INTEGER, rcnts, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      total_scnt = 0
      total_rcnt = 0
      do p = 0, nprocs-1
         total_scnt = total_scnt + scnts(p)
         total_rcnt = total_rcnt + rcnts(p)
      end do

! --- Make alltoallv params for send. ------------------------------------------
      allocate(table(rcnt))
      allocate(table_all(total_rcnt))

      displs(0) = 0
      do p = 1, nprocs-1
         displs(p) = displs(p-1) + rcnts(p-1)
      end do

      cnt = 0
      do k = 1, fxo%np
         ix = fxo%cndx0(k,1)
         iy = fxo%cndx0(k,2)
         if(parent%kx <= ix .and. ix <= parent%kxend .and. parent%ky <= iy .and. iy <= parent%kyend) then
            cnt = cnt + 1
            table(cnt) = k
         end if
      end do
#ifndef MULTI
      call MPI_Allgatherv(table, rcnt, MPI_INTEGER, table_all, rcnts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgatherv(table, rcnt, MPI_INTEGER, table_all, rcnts, displs, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      
      snp0 = 0
      do p = 0, nprocs-1
         do ind = displs(p)+1, displs(p)+rcnts(p)
            k = table_all(ind)
            ix = fxo%fndx(k,1)
            iy = fxo%fndx(k,2)
            if(my%ix <= ix .and. ix <= my%ixend .and. my%iy <= iy .and. iy <= my%iyend) then
               snp0 = snp0 + 1
            end if
         end do
      end do

      allocate(fxo%fndx0_l(snp0,2))
      allocate(fxo%sendcnts0(0:nprocs-1))
      allocate(fxo%sdispls0(0:nprocs-1))
      fndx0_l   => fxo%fndx0_l
      sendcnts0 => fxo%sendcnts0
      sdispls0  => fxo%sdispls0

      snp0 = 0
      do p = 0, nprocs-1
         cnt = 0
         do ind = displs(p)+1, displs(p)+rcnts(p)
            k = table_all(ind)
            ix = fxo%fndx(k,1)
            iy = fxo%fndx(k,2)
            if(my%ix <= ix .and. ix <= my%ixend .and. my%iy <= iy .and. iy <= my%iyend) then
               cnt = cnt + 1
               snp0 = snp0 + 1
               fndx0_l(snp0,1) = ix - my%kx + 1 ! with egde
               fndx0_l(snp0,2) = iy - my%ky + 1 ! with edge
            end if
         end do
         sendcnts0(p) = cnt
      end do

      sdispls0(0) = 0
      do p = 1, nprocs-1
         sdispls0(p) = sdispls0(p-1) + sendcnts0(p-1)
      end do

      deallocate(table)
      deallocate(table_all)

! --- Make alltoallv params for recv. ------------------------------------------
      allocate(table(scnt))
      allocate(table_all(total_scnt))

      displs(0) = 0
      do p = 1, nprocs-1
         displs(p) = displs(p-1) + scnts(p-1)
      end do

      cnt = 0
      do k = 1, fxo%np
         ix = fxo%fndx(k,1)
         iy = fxo%fndx(k,2)
         if(my%ix <= ix .and. ix <= my%ixend .and. my%iy <= iy .and. iy <= my%iyend) then
            cnt = cnt + 1
            table(cnt) = k
         end if
      end do
#ifndef MULTI
      call MPI_Allgatherv(table, scnt, MPI_INTEGER, table_all, scnts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgatherv(table, scnt, MPI_INTEGER, table_all, scnts, displs, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      
      rnp0 = 0
      do p = 0, nprocs-1
         do ind = displs(p)+1, displs(p)+scnts(p)
            k = table_all(ind)
            ix = fxo%cndx0(k,1)
            iy = fxo%cndx0(k,2)
            if(parent%kx <= ix .and. ix <= parent%kxend .and. parent%ky <= iy .and. iy <= parent%kyend) then
               rnp0 = rnp0 + 1
            end if
         end do
      end do

      allocate(fxo%cndx0_l(rnp0,2))
      allocate(fxo%recvcnts0(0:nprocs-1))
      allocate(fxo%rdispls0(0:nprocs-1))
      cndx0_l   => fxo%cndx0_l
      recvcnts0 => fxo%recvcnts0
      rdispls0  => fxo%rdispls0

      rnp0 = 0
      do p = 0, nprocs-1
         cnt = 0
         do ind = displs(p)+1, displs(p)+scnts(p)
            k = table_all(ind)
            ix = fxo%cndx0(k,1)
            iy = fxo%cndx0(k,2)
            if(parent%kx <= ix .and. ix <= parent%kxend .and. parent%ky <= iy .and. iy <= parent%kyend) then
               cnt = cnt + 1
               rnp0 = rnp0 + 1
               cndx0_l(rnp0,1) = ix - parent%kx + 1 ! with egde
               cndx0_l(rnp0,2) = iy - parent%ky + 1 ! with edge
            end if
         end do
         recvcnts0(p) = cnt
      end do

      rdispls0(0) = 0
      do p = 1, nprocs-1
         rdispls0(p) = rdispls0(p-1) + recvcnts0(p-1)
      end do

      deallocate(table)
      deallocate(table_all)

! ------------------------------------------------------------------------------

      deallocate(scnts)
      deallocate(rcnts)
      deallocate(displs)
#ifndef SINGLE_A2A
#ifdef A2A3D
      call A2A3D_plan(fxo%sendcnts0, fxo%handler0)
#endif
#endif

      return
   end subroutine make_a2a_params_c2c

   subroutine make_a2a_params_i2f(dg, fxi, nprocs, ixst, iyst, ixen, iyen)
#ifndef SINGLE_A2A
#ifdef A2A3D
      use mod_a2a3d
#endif
#endif
      type(data_grids), target, intent(inout) :: dg
      type(interp_info), target, intent(inout) :: fxi
      integer(kind=4), intent(in) :: nprocs
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
      integer(kind=4), intent(in) :: ixst, iyst
! ==============================================================================
      integer(kind=4), intent(in) :: ixen, iyen

      ! Pointers to type objecs.
      type(grid_info), pointer :: parent, my
      integer(kind=4), pointer :: snp0, rnp0
      integer(kind=4), pointer :: snp1, rnp1
      integer(kind=4), pointer, dimension(:,:) :: fndx0_l, cndx0_l
      integer(kind=4), pointer, dimension(:,:) :: fndx1_l, cndx1_l
      integer(kind=4), pointer, dimension(:) :: sendcnts0, sdispls0, recvcnts0, rdispls0
      integer(kind=4), pointer, dimension(:) :: sendcnts1, sdispls1, recvcnts1, rdispls1
      real(kind=REAL_BYTE), pointer, dimension(:) :: wt0_l, wt1_l

      ! Work vals.
      integer(kind=4) :: k, ix, iy, p, cnt, ind, ierr
      integer(kind=4) :: scnt, rcnt, total_scnt, total_rcnt

      ! Work arrays.
      integer(kind=4), allocatable, dimension(:) :: scnts, rcnts, displs
      integer(kind=4), allocatable, dimension(:) :: table, table_all

      parent => dg%parent
      my => dg%my

! === Section 1: ===============================================================
! === snp0, rnp0, fndx0_l, cndx0_l, sendcnts0, sdispls0, recvcnts0, rdispls0 ===
! ==============================================================================
      snp0 => fxi%snp0
      rnp0 => fxi%rnp0

      allocate(scnts(0:nprocs-1))
      allocate(rcnts(0:nprocs-1))
      allocate(displs(0:nprocs-1))

      ! Count send elements.
      scnt = 0
      do k = 1, fxi%np
         ix = fxi%cndx0(k,1)
         iy = fxi%cndx0(k,2)
         if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
            scnt = scnt + 1
         end if
      end do
#ifndef MULTI
      call MPI_Allgather(scnt, 1, MPI_INTEGER, scnts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(scnt, 1, MPI_INTEGER, scnts, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      ! Count recv elements.
      rcnt = 0
      do k = 1, fxi%np
         ix = fxi%fndx(k,1)
         iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!        if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
         if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
            rcnt = rcnt + 1
         end if
      end do
#ifndef MULTI
      call MPI_Allgather(rcnt, 1, MPI_INTEGER, rcnts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(rcnt, 1, MPI_INTEGER, rcnts, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      total_scnt = 0
      total_rcnt = 0
      do p = 0, nprocs-1
         total_scnt = total_scnt + scnts(p)
         total_rcnt = total_rcnt + rcnts(p)
      end do

! --- Make alltoallv params for send. ------------------------------------------
      allocate(table(rcnt))
      allocate(table_all(total_rcnt))

      displs(0) = 0
      do p = 1, nprocs-1
         displs(p) = displs(p-1) + rcnts(p-1)
      end do

      cnt = 0
      do k = 1, fxi%np
         ix = fxi%fndx(k,1)
         iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!        if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
         if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
            cnt = cnt + 1
            table(cnt) = k
         end if
      end do
#ifndef MULTI
      call MPI_Allgatherv(table, rcnt, MPI_INTEGER, table_all, rcnts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgatherv(table, rcnt, MPI_INTEGER, table_all, rcnts, displs, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      
      snp0 = 0
      do p = 0, nprocs-1
         do ind = displs(p)+1, displs(p)+rcnts(p)
            k = table_all(ind)
            ix = fxi%cndx0(k,1)
            iy = fxi%cndx0(k,2)
            if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
               snp0 = snp0 + 1
            end if
         end do
      end do

      allocate(fxi%cndx0_l(snp0,2))
      allocate(fxi%wt0_l(snp0))
      allocate(fxi%sendcnts0(0:nprocs-1))
      allocate(fxi%sdispls0(0:nprocs-1))
      cndx0_l   => fxi%cndx0_l
      wt0_l     => fxi%wt0_l
      sendcnts0 => fxi%sendcnts0
      sdispls0  => fxi%sdispls0

      snp0 = 0
      do p = 0, nprocs-1
         cnt = 0
         do ind = displs(p)+1, displs(p)+rcnts(p)
            k = table_all(ind)
            ix = fxi%cndx0(k,1)
            iy = fxi%cndx0(k,2)
            if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
               cnt = cnt + 1
               snp0 = snp0 + 1
               cndx0_l(snp0,1) = ix - parent%kx + 1 ! with egde
               cndx0_l(snp0,2) = iy - parent%ky + 1 ! with edge
               wt0_l(snp0) = fxi%wt0(k)
            end if
         end do
         sendcnts0(p) = cnt
      end do

      sdispls0(0) = 0
      do p = 1, nprocs-1
         sdispls0(p) = sdispls0(p-1) + sendcnts0(p-1)
      end do

      deallocate(table)
      deallocate(table_all)

! --- Make alltoallv params for recv. ------------------------------------------
      allocate(table(scnt))
      allocate(table_all(total_scnt))

      displs(0) = 0
      do p = 1, nprocs-1
         displs(p) = displs(p-1) + scnts(p-1)
      end do

      cnt = 0
      do k = 1, fxi%np
         ix = fxi%cndx0(k,1)
         iy = fxi%cndx0(k,2)
         if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
            cnt = cnt + 1
            table(cnt) = k
         end if
      end do
#ifndef MULTI
      call MPI_Allgatherv(table, scnt, MPI_INTEGER, table_all, scnts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgatherv(table, scnt, MPI_INTEGER, table_all, scnts, displs, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      
      rnp0 = 0
      do p = 0, nprocs-1
         do ind = displs(p)+1, displs(p)+scnts(p)
            k = table_all(ind)
            ix = fxi%fndx(k,1)
            iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!           if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
            if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
               rnp0 = rnp0 + 1
            end if
         end do
      end do

      allocate(fxi%fndx0_l(rnp0,2))
      allocate(fxi%recvcnts0(0:nprocs-1))
      allocate(fxi%rdispls0(0:nprocs-1))
      fndx0_l   => fxi%fndx0_l
      recvcnts0 => fxi%recvcnts0
      rdispls0  => fxi%rdispls0

      rnp0 = 0
      do p = 0, nprocs-1
         cnt = 0
         do ind = displs(p)+1, displs(p)+scnts(p)
            k = table_all(ind)
            ix = fxi%fndx(k,1)
            iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!           if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
            if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
               cnt = cnt + 1
               rnp0 = rnp0 + 1
               fndx0_l(rnp0,1) = ix - my%kx + 1 ! with egde
               fndx0_l(rnp0,2) = iy - my%ky + 1 ! with edge
            end if
         end do
         recvcnts0(p) = cnt
      end do

      rdispls0(0) = 0
      do p = 1, nprocs-1
         rdispls0(p) = rdispls0(p-1) + recvcnts0(p-1)
      end do

      deallocate(table)
      deallocate(table_all)

! ------------------------------------------------------------------------------

      deallocate(scnts)
      deallocate(rcnts)
      deallocate(displs)

! === Section 2: ===============================================================
! === snp1, rnp1, fndx1_l, cndx1_l, sendcnts1, sdispls1, recvcnts1, rdispls1 ===
! ==============================================================================
      snp1 => fxi%snp1
      rnp1 => fxi%rnp1

      allocate(scnts(0:nprocs-1))
      allocate(rcnts(0:nprocs-1))
      allocate(displs(0:nprocs-1))

      ! Count send elements.
      scnt = 0
      do k = 1, fxi%np
         ix = fxi%cndx1(k,1)
         iy = fxi%cndx1(k,2)
         if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
            scnt = scnt + 1
         end if
      end do
#ifndef MULTI
      call MPI_Allgather(scnt, 1, MPI_INTEGER, scnts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(scnt, 1, MPI_INTEGER, scnts, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      ! Count recv elements.
      rcnt = 0
      do k = 1, fxi%np
         ix = fxi%fndx(k,1)
         iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!        if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
         if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
            rcnt = rcnt + 1
         end if
      end do
#ifndef MULTI
      call MPI_Allgather(rcnt, 1, MPI_INTEGER, rcnts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgather(rcnt, 1, MPI_INTEGER, rcnts, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif

      total_scnt = 0
      total_rcnt = 0
      do p = 0, nprocs-1
         total_scnt = total_scnt + scnts(p)
         total_rcnt = total_rcnt + rcnts(p)
      end do

! --- Make alltoallv params for send. ------------------------------------------
      allocate(table(rcnt))
      allocate(table_all(total_rcnt))

      displs(0) = 0
      do p = 1, nprocs-1
         displs(p) = displs(p-1) + rcnts(p-1)
      end do

      cnt = 0
      do k = 1, fxi%np
         ix = fxi%fndx(k,1)
         iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!        if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
         if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
            cnt = cnt + 1
            table(cnt) = k
         end if
      end do
#ifndef MULTI
      call MPI_Allgatherv(table, rcnt, MPI_INTEGER, table_all, rcnts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgatherv(table, rcnt, MPI_INTEGER, table_all, rcnts, displs, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      
      snp1 = 0
      do p = 0, nprocs-1
         do ind = displs(p)+1, displs(p)+rcnts(p)
            k = table_all(ind)
            ix = fxi%cndx1(k,1)
            iy = fxi%cndx1(k,2)
            if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
               snp1 = snp1 + 1
            end if
         end do
      end do

      allocate(fxi%cndx1_l(snp1,2))
      allocate(fxi%wt1_l(snp1))
      allocate(fxi%sendcnts1(0:nprocs-1))
      allocate(fxi%sdispls1(0:nprocs-1))
      cndx1_l   => fxi%cndx1_l
      wt1_l     => fxi%wt1_l
      sendcnts1 => fxi%sendcnts1
      sdispls1  => fxi%sdispls1

      snp1 = 0
      do p = 0, nprocs-1
         cnt = 0
         do ind = displs(p)+1, displs(p)+rcnts(p)
            k = table_all(ind)
            ix = fxi%cndx1(k,1)
            iy = fxi%cndx1(k,2)
            if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
               cnt = cnt + 1
               snp1 = snp1 + 1
               cndx1_l(snp1,1) = ix - parent%kx + 1 ! with egde
               cndx1_l(snp1,2) = iy - parent%ky + 1 ! with edge
               wt1_l(snp1) = fxi%wt1(k)
            end if
         end do
         sendcnts1(p) = cnt
      end do

      sdispls1(0) = 0
      do p = 1, nprocs-1
         sdispls1(p) = sdispls1(p-1) + sendcnts1(p-1)
      end do

      deallocate(table)
      deallocate(table_all)

! --- Make alltoallv params for recv. ------------------------------------------
      allocate(table(scnt))
      allocate(table_all(total_scnt))

      displs(0) = 0
      do p = 1, nprocs-1
         displs(p) = displs(p-1) + scnts(p-1)
      end do

      cnt = 0
      do k = 1, fxi%np
         ix = fxi%cndx1(k,1)
         iy = fxi%cndx1(k,2)
         if(parent%ix <= ix .and. ix <= parent%ixend .and. parent%iy <= iy .and. iy <= parent%iyend) then
            cnt = cnt + 1
            table(cnt) = k
         end if
      end do
#ifndef MULTI
      call MPI_Allgatherv(table, scnt, MPI_INTEGER, table_all, scnts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allgatherv(table, scnt, MPI_INTEGER, table_all, scnts, displs, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      
      rnp1 = 0
      do p = 0, nprocs-1
         do ind = displs(p)+1, displs(p)+scnts(p)
            k = table_all(ind)
            ix = fxi%fndx(k,1)
            iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!           if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
            if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
               rnp1 = rnp1 + 1
            end if
         end do
      end do

      allocate(fxi%fndx1_l(rnp1,2))
      allocate(fxi%recvcnts1(0:nprocs-1))
      allocate(fxi%rdispls1(0:nprocs-1))
      fndx1_l   => fxi%fndx1_l
      recvcnts1 => fxi%recvcnts1
      rdispls1  => fxi%rdispls1

      rnp1 = 0
      do p = 0, nprocs-1
         cnt = 0
         do ind = displs(p)+1, displs(p)+scnts(p)
            k = table_all(ind)
            ix = fxi%fndx(k,1)
            iy = fxi%fndx(k,2)
! === DEBUG for wave hight gap on nest boundary. 2012/10/30 ====================
!           if(my%kx <= ix .and. ix <= my%kxend .and. my%ky <= iy .and. iy <= my%kyend) then
            if(ixst <= ix .and. ix <= ixen .and. iyst <= iy .and. iy <= iyen) then
! ==============================================================================
               cnt = cnt + 1
               rnp1 = rnp1 + 1
               fndx1_l(rnp1,1) = ix - my%kx + 1 ! with egde
               fndx1_l(rnp1,2) = iy - my%ky + 1 ! with edge
            end if
         end do
         recvcnts1(p) = cnt
      end do

      rdispls1(0) = 0
      do p = 1, nprocs-1
         rdispls1(p) = rdispls1(p-1) + recvcnts1(p-1)
      end do

      deallocate(table)
      deallocate(table_all)

! ------------------------------------------------------------------------------

      deallocate(scnts)
      deallocate(rcnts)
      deallocate(displs)
#ifndef SINGLE_A2A
#ifdef A2A3D
      call A2A3D_plan(fxi%sendcnts0, fxi%handler0)
      call A2A3D_plan(fxi%sendcnts1, fxi%handler1)
#endif
#endif

      return
   end subroutine make_a2a_params_i2f
! === USE_MPI_ALLTOALLV ========================================================
#endif
#endif

   subroutine initl_gridmap_dz(dg)
      type(data_grids), target, intent(inout) :: dg

      type(grid_info), pointer :: parent, my
      type(interp_info), pointer :: dzi
      integer(kind=4) :: ix, iy, k
      real(kind=REAL_BYTE) :: lfac

      parent => dg%parent
      my => dg%my

      lfac = 1.0/REAL_FUNC(my%nr)

#ifndef MPI
! === SOLVE WITH FLOW VOLUME! ==================================================
!     dg%dzi%np = (my%nx - 1) + (my%ny - 1)
      dg%dzi%np = 2*my%nx + 2*(my%ny + 2)
! === SOLVE WITH FLOW VOLUME! ==================================================
#else
! === SOLVE WITH FLOW VOLUME! ==================================================
!     dg%dzi%np = (my%totalNx - 1) + (my%totalNy - 1)
      dg%dzi%np = 2*my%totalNx + 2*(my%totalNy + 2)
! === SOLVE WITH FLOW VOLUME! ==================================================
#endif
      allocate(dg%dzi%fndx(dg%dzi%np,2))
      allocate(dg%dzi%cndx0(dg%dzi%np,2))
      allocate(dg%dzi%cndx1(dg%dzi%np,2))
      allocate(dg%dzi%wt0(dg%dzi%np))
      allocate(dg%dzi%wt1(dg%dzi%np))

      dzi => dg%dzi

! === SOLVE WITH FLOW VOLUME! ==================================================
#if 0
! === SOLVE WITH FLOW VOLUME! ==================================================
      iy = 0
      k = 1
#ifndef MPI
      do ix = 1, my%nx-1
#else
      do ix = 1, my%totalNx-1
#endif
         dzi%fndx(k,1) = ix
         dzi%fndx(k,2) = iy
         dzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         dzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr - 1
         dzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         dzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr - 1
         dzi%wt0(k) = 1.0d0 - mod(ix-1,my%nr)*lfac
         dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
         k = k+1
      end do

      ix = 0
#ifndef MPI
      do iy = my%ny, 2, -1
#else
      do iy = my%totalNy, 2, -1
#endif
         dzi%fndx(k,1) = ix
         dzi%fndx(k,2) = iy
         dzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr - 1
         dzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         dzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr - 1
         dzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         dzi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
         k = k+1
      end do
! === SOLVE WITH FLOW VOLUME! ==================================================
#else
      k = 1

      iy = 0
#ifndef MPI
      do ix = 1, my%nx
#else
      do ix = 1, my%totalNx
#endif
         dzi%fndx(k,1) = ix
         dzi%fndx(k,2) = iy
         dzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         dzi%cndx0(k,2) = my%zeroIY - 1
         dzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         dzi%cndx1(k,2) = my%zeroIY - 1
         dzi%wt0(k) = 1.0d0 - mod(ix-1,my%nr)*lfac
         dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
         k = k+1
      end do

#ifndef MPI
      iy = my%ny + 1
      do ix = my%nx, 1, -1
#else
      iy = my%totalNy + 1
      do ix = my%totalNx, 1, -1
#endif
         dzi%fndx(k,1) = ix
         dzi%fndx(k,2) = iy
         dzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr
         dzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         dzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         dzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         dzi%wt0(k) = 1.0d0 - mod(ix-1,my%nr)*lfac
         dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
         k = k+1
      end do

      ix = 0
#ifndef MPI
      do iy = my%ny + 1, 1, -1
#else
      do iy = my%totalNy + 1, 1, -1
#endif
         dzi%fndx(k,1) = ix
         dzi%fndx(k,2) = iy
         dzi%cndx0(k,1) = my%zeroIX - 1
         dzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         dzi%cndx1(k,1) = my%zeroIX - 1
         dzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         dzi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
         k = k+1
      end do
      iy = 0
      dzi%fndx(k,1) = ix
      dzi%fndx(k,2) = iy
      dzi%cndx0(k,1) = my%zeroIX - 1
      dzi%cndx0(k,2) = my%zeroIY - 1
      dzi%cndx1(k,1) = my%zeroIX - 1
      dzi%cndx1(k,2) = my%zeroIY
      dzi%wt0(k) = lfac
      dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
      k = k+1

#ifndef MPI
      ix = my%nx + 1
#else
      ix = my%totalNx + 1
#endif
      iy = 0
      dzi%fndx(k,1) = ix
      dzi%fndx(k,2) = iy
      dzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr + 1
      dzi%cndx0(k,2) = my%zeroIY
      dzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
      dzi%cndx1(k,2) = my%zeroIY - 1
      dzi%wt0(k) = 1.0d0 - lfac
      dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
      k = k+1
#ifndef MPI
      do iy = 1, my%ny + 1
#else
      do iy = 1, my%totalNy + 1
#endif
         dzi%fndx(k,1) = ix
         dzi%fndx(k,2) = iy
         dzi%cndx0(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         dzi%cndx0(k,2) = my%zeroIY + (iy-1)/my%nr
         dzi%cndx1(k,1) = my%zeroIX + (ix-1)/my%nr + 1
         dzi%cndx1(k,2) = my%zeroIY + (iy-1)/my%nr + 1
         dzi%wt0(k) = 1.0d0 - mod(iy-1,my%nr)*lfac
         dzi%wt1(k) = 1.0d0 - dzi%wt0(k)
         k = k+1
      end do
#endif
! === SOLVE WITH FLOW VOLUME! ==================================================
#ifdef MPI
!#ifndef USE_ALLTOALLV
      allocate(dzi%fb(dzi%np))
      allocate(dzi%cb(dzi%np))
!#endif
#endif

      return
   end subroutine initl_gridmap_dz

   subroutine interp2fine_dz(cg,fg)
      type(data_grids), target, intent(inout) :: cg, fg

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dzc, dzf
      integer(kind=4) :: k
      type(interp_info), pointer :: dzi

      integer(kind=4) :: i, j
#ifdef MPI
      real(kind=REAL_BYTE), pointer, dimension(:) :: zfbuf, zcbuf
      integer(kind=4) :: ix, iy, ierr
      integer(kind=4) :: ist, jst
! === Upwind3 ==================================================================
      integer(kind=4) :: ien, jen
! ==============================================================================
#endif
      real(kind=REAL_BYTE) :: t0, t1

      dzc => cg%depth_field%dz
      dzf => fg%depth_field%dz
      dzi => fg%dzi

#ifndef MPI
!$omp parallel
!$omp do
      do k = 1, dzi%np
         dzf(dzi%fndx(k,1),dzi%fndx(k,2)) = &
            dzc(dzi%cndx0(k,1),dzi%cndx0(k,2))*dzi%wt0(k) + &
            dzc(dzi%cndx1(k,1),dzi%cndx1(k,2))*dzi%wt1(k)
      end do

!$omp single
      t0 = 1.0d0/REAL_FUNC(fg%my%nr)
      t1 = 1.0d0 - t0
!$omp end single
!$omp do
      do i = 1, fg%my%nx
! === Upwind3 ==================================================================
         dzf(i,-1        ) = dzf(i,0         )*t1 + dzf(i,1       )*t0
         dzf(i,fg%my%ny+2) = dzf(i,fg%my%ny+1)*t1 + dzf(i,fg%my%ny)*t0
! ==============================================================================
         dzf(i,0         ) = dzf(i,0         )*t0 + dzf(i,1       )*t1
         dzf(i,fg%my%ny+1) = dzf(i,fg%my%ny+1)*t0 + dzf(i,fg%my%ny)*t1
      end do
!$omp do
      do j = 0, fg%my%ny+1
! === Upwind3 ==================================================================
         dzf(-1,        j) = dzf(0,         j)*t1 + dzf(1,       j)*t0
         dzf(fg%my%nx+2,j) = dzf(fg%my%nx+1,j)*t1 + dzf(fg%my%nx,j)*t0
! ==============================================================================
         dzf(0,         j) = dzf(0,         j)*t0 + dzf(1,       j)*t1
         dzf(fg%my%nx+1,j) = dzf(fg%my%nx+1,j)*t0 + dzf(fg%my%nx,j)*t1
      end do
!$omp end parallel
#else
!$omp parallel
!$omp single
      zfbuf => dzi%fb
      zcbuf => dzi%cb
!$omp end single
!$omp do
      do k = 1, dzi%np
         zfbuf(k) = 0.0d0
         zcbuf(k) = 0.0d0
      end do
      !*==============*
      !*  coarse2buf  * must not read from edges
      !*==============*
!$omp do private(ix, iy, i, j)
      do k = 1, dzi%np
         ix = dzi%cndx0(k,1)
         iy = dzi%cndx0(k,2)
         if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
            i = ix - cg%my%kx + 1
            j = iy - cg%my%ky + 1
            zcbuf(k) = zcbuf(k) + dzc(i,j)*dzi%wt0(k)
         end if
      end do
!$omp do private(ix, iy, i, j)
      do k = 1, dzi%np
         ix = dzi%cndx1(k,1)
         iy = dzi%cndx1(k,2)
         if(cg%my%ix <= ix .and. ix <= cg%my%ixend .and. cg%my%iy <= iy .and. iy <= cg%my%iyend) then
            i = ix - cg%my%kx + 1 ! with edges
            j = iy - cg%my%ky + 1 ! with edges
            zcbuf(k) = zcbuf(k) + dzc(i,j)*dzi%wt1(k)
         end if
      end do
      !*==============*
      !*  allreduce   *
      !*==============*
!$omp single
#ifndef MULTI
      call MPI_Allreduce(zcbuf, zfbuf, dzi%np, REAL_MPI, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(zcbuf, zfbuf, dzi%np, REAL_MPI, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
      if(ierr /= 0) then
         select case (ierr)
            case(MPI_ERR_BUFFER)
               write(0,'(a)') 'MPI Error : Invalid buffer pointer'
            case(MPI_ERR_COUNT)
               write(0,'(a)') 'MPI Error : Invalid count argument'
            case(MPI_ERR_TYPE)
               write(0,'(a)') 'MPI Error : Invalid datatype argument'
            case(MPI_ERR_OP)
               write(0,'(a)') 'MPI Error : Invalid operation'
            case(MPI_ERR_COMM)
               write(0,'(a)') 'MPI Error : Invalid communicator'
            case default
               write(0,'(a)') 'MPI Error : Unknown error'
         end select
         call fatal_error(ierr)
      end if
      !*==============*
      !*  buf2fine    * write to edges
      !*==============*
      ist = fg%my%kx
! === Upwind3 ==================================================================
!     if(iand(fg%my%has_boundary, WEST_BOUND)  /= 0) ist = 0
      if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) ist = -1
      ien = fg%my%kxend+1
      if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) ien = fg%my%kxend+2
! ==============================================================================
      jst = fg%my%ky
! === Upwind3 ==================================================================
!     if(iand(fg%my%has_boundary, NORTH_BOUND)  /= 0) jst = 0
      if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) jst = -1
      jen = fg%my%kyend+1
      if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) jen = fg%my%kyend+2
! ==============================================================================
!$omp end single
!$omp do private(ix, iy, i, j)
      do k = 1, dzi%np
         ix = dzi%fndx(k,1)
         iy = dzi%fndx(k,2)
! === Upwind3 ==================================================================
!        if(ist <= ix .and. ix <= fg%my%kxend+1 .and. jst <= iy .and. iy <= fg%my%kyend+1) then
         if(ist <= ix .and. ix <= ien .and. jst <= iy .and. iy <= jen) then
! ==============================================================================
            i = ix - fg%my%kx + 1 ! with edges
            j = iy - fg%my%ky + 1 ! with edges
            dzf(i,j) = zfbuf(k)
         end if
      end do
!$omp single
      t0 = 1.0d0/REAL_FUNC(fg%my%nr)
      t1 = 1.0d0 - t0
!$omp end single

      if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) then
!$omp single
         ist = 1
!$omp end single
!$omp do
         do i = ist, fg%my%nx
! === Upwind3 ==================================================================
            dzf(i,-1) = dzf(i,0)*t1 + dzf(i,1)*t0
! ==============================================================================
            dzf(i,0) = dzf(i,0)*t0 + dzf(i,1)*t1
         end do
      end if

      if(iand(fg%my%has_boundary, SOUTH_BOUND) /= 0) then
!$omp single
         ist = 1
!$omp end single
!$omp do
         do i = ist, fg%my%nx
! === Upwind3 ==================================================================
            dzf(i,fg%my%ny+2) = dzf(i,fg%my%ny+1)*t1 + dzf(i,fg%my%ny)*t0
! ==============================================================================
            dzf(i,fg%my%ny+1) = dzf(i,fg%my%ny+1)*t0 + dzf(i,fg%my%ny)*t1
         end do
      end if

      if(iand(fg%my%has_boundary, WEST_BOUND) /= 0) then
!$omp single
         jst = 1
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) jst = 0
!$omp end single
!$omp do
         do j = jst, fg%my%ny+1
! === Upwind3 ==================================================================
            dzf(-1,j) = dzf(0,j)*t1 + dzf(1,j)*t0
! ==============================================================================
            dzf(0,j) = dzf(0,j)*t0 + dzf(1,j)*t1
         end do
      end if

      if(iand(fg%my%has_boundary, EAST_BOUND) /= 0) then
!$omp single
         jst = 1
         if(iand(fg%my%has_boundary, NORTH_BOUND) /= 0) jst = 0
!$omp end single
!$omp do
         do j = jst, fg%my%ny+1
! === Upwind3 ==================================================================
            dzf(fg%my%nx+2,j) = dzf(fg%my%nx+1,j)*t1 + dzf(fg%my%nx,j)*t0
! ==============================================================================
            dzf(fg%my%nx+1,j) = dzf(fg%my%nx+1,j)*t0 + dzf(fg%my%nx,j)*t1
         end do
      end if
!$omp end parallel
#endif

      return
   end subroutine interp2fine_dz

end module mod_nest
