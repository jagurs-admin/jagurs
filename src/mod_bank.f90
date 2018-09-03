!#define OLDFORMAT
#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_bank
#ifdef MPI
use mpi
#endif
use mod_grid
implicit none

contains

#ifdef BANKFILE
   subroutine read_bank_file(dg)
      type(data_grids), target, intent(inout) :: dg

      integer(kind=4), pointer :: nx, ny
#if !defined(MPI) || !defined(ONEFILE)
#ifdef CARTESIAN
      real(kind=REAL_BYTE), pointer :: mlon0, mlat0, dh
#else
      real(kind=REAL_BYTE), pointer :: mlon0, mlat0
      real(kind=REAL_BYTE) :: dh
#endif
#else
#ifdef CARTESIAN
      real(kind=REAL_BYTE), pointer :: dh
#else
      real(kind=REAL_BYTE) :: dh
#endif
      integer(kind=4), pointer :: kx, kyend, totalNy
      real(kind=REAL_BYTE) :: mlon0, mlat0
#endif
      integer(kind=4), pointer, dimension(:,:) :: ir
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: btx, bty

#ifndef CARTESIAN
      real(kind=REAL_BYTE), allocatable, dimension(:) :: xin, yin
      real(kind=REAL_BYTE) :: xtmpin, ytmpin
#endif
#ifdef OLDFORMAT
      integer(kind=4), allocatable, dimension(:) :: x, y, val

#ifdef CARTESIAN
      integer(kind=4) :: n, i, j, num_lines, xtmp, ytmp, valtmp
#else
      integer(kind=4) :: n, i, j, num_lines, valtmp
#endif
#else
      integer(kind=4), allocatable, dimension(:) :: x, y, irread
      real(kind=REAL_BYTE), allocatable, dimension(:) :: height

#ifdef CARTESIAN
      integer(kind=4) :: n, i, j, num_lines, xtmp, ytmp, irtmp
#else
      integer(kind=4) :: n, i, j, num_lines, irtmp
#endif
      real(kind=REAL_BYTE) :: heighttmp
#endif
      real(kind=REAL_BYTE) :: bh

      if(dg%bank_file /= 'NO_BANK_FILE_GIVEN') then

         write(6,'(a,a,a)') '[bank] Open bank file ', trim(dg%bank_file), '!'
         open(1,file=trim(dg%bank_file),action='read',status='old',form='formatted')

         num_lines = 0
         do while(.true.)
#ifdef CARTESIAN
#ifdef OLDFORMAT
            read(1,*,end=100) xtmp, ytmp, valtmp
#else
            read(1,*,end=100) ytmp, xtmp, irtmp, heighttmp
#endif
#else
#ifdef OLDFORMAT
            read(1,*,end=100) xtmpin, ytmpin, valtmp
#else
            read(1,*,end=100) ytmpin, xtmpin, irtmp, heighttmp
#endif
#endif
            num_lines = num_lines + 1
         end do
 100     write(6,'(a,i6)') '[bank] Number of banks: ', num_lines

#ifndef CARTESIAN
         allocate(xin(num_lines))
         allocate(yin(num_lines))
#endif
         allocate(x(num_lines))
         allocate(y(num_lines))
#ifdef OLDFORMAT
         allocate(val(num_lines))
#else
         allocate(irread(num_lines))
         allocate(height(num_lines))
#endif

         rewind(1)
         do n = 1, num_lines
#ifdef CARTESIAN
#ifdef OLDFORMAT
            read(1,*) x(n), y(n), val(n)
#else
            read(1,*) y(n), x(n), irread(n), height(n)
#endif
#else
#ifdef OLDFORMAT
            read(1,*) xin(n), yin(n), val(n)
#else
            read(1,*) yin(n), xin(n), irread(n), height(n)
#endif
#endif
         end do

         nx    => dg%my%nx
         ny    => dg%my%ny
#ifdef CARTESIAN
#if !defined(MPI) || !defined(ONEFILE)
         mlon0 => dg%my%mlon0
         mlat0 => dg%my%mlat0
         dh    => dg%my%dh
#else
         kx      => dg%my%kx
         kyend   => dg%my%kyend
         totalNy => dg%my%totalNy
         dh      => dg%my%dh

         mlon0 = dg%my%mlon0 + (kx - 1)*dh
         mlat0 = dg%my%mlat0 + (totalNy - 1)*dh - (kyend - 1)*dh
#endif
#else
#if !defined(MPI) || !defined(ONEFILE)
         mlon0 => dg%my%mlon0
         mlat0 => dg%my%mlat0
         dh    =  dg%my%dh*60.0d0
#else
         mlon0 = dg%my%mlon0
         mlat0 = dg%my%mlat0
         dh    = dg%my%dh*60.0d0
#endif
#endif
#ifndef CARTESIAN
         if(mlon0 < 0.0d0) then
            do n = 1, num_lines
               xin(n) = (360.0d0+xin(n))*60.0d0
               yin(n) = (90.0d0-yin(n))*60.0d0
            end do
         else
            do n = 1, num_lines
               xin(n) = xin(n)*60.0d0
               yin(n) = (90.0d0-yin(n))*60.0d0
            end do
         end if
#endif

         ir    => dg%wave_field%ir
         btx   => dg%wave_field%btx
         bty   => dg%wave_field%bty

         do n = 1, num_lines
!           x(n) = int((x(n) - mlon0 + 0.5d0)/dh) + 1
!           y(n) = int((y(n) - mlat0 + 0.5d0)/dh) + 1
#ifdef CARTESIAN
            x(n) = floor((x(n) - mlon0 + 0.5d0)/dh) + 1
            y(n) = floor((y(n) - mlat0 + 0.5d0)/dh) + 1
            y(n) = ny - y(n) + 1
#else
            x(n) = floor((xin(n) - mlon0)/dh) + 1
            y(n) = floor((yin(n) - mlat0)/dh) + 1
#endif

#ifdef MPI
            if((x(n) >= 0) .and. (x(n) <= nx+1) .and. (y(n) >= 0) .and. (y(n) <= ny+1)) then
#endif
#ifdef OLDFORMAT
            if(val(n)/300000 /= 0) then
               i = x(n)
               j = y(n)
               ir(i,j) = 3
               bh = dble(val(n) - 300000)/100.0d0
               btx(i,j  ) = bh
               bty(i,j-1) = bh
            else if(val(n)/200000 /= 0) then
               i = x(n)
               j = y(n)
               ir(i,j) = 2
               bh = dble(val(n) - 200000)/100.0d0
               bty(i,j-1) = bh
            else if(val(n)/100000 /= 0) then
               i = x(n)
               j = y(n)
               ir(i,j) = 1
               bh = dble(val(n) - 100000)/100.0d0
               btx(i,j  ) = bh
            end if
#else
            if(irread(n) == 3) then
               i = x(n)
               j = y(n)
               ir(i,j) = 3
               bh = height(n)
               btx(i,j  ) = bh
               bty(i,j-1) = bh
            else if(irread(n) == 2) then
               i = x(n)
               j = y(n)
               ir(i,j) = 2
               bh = height(n)
               bty(i,j-1) = bh
            else if(irread(n) == 1) then
               i = x(n)
               j = y(n)
               ir(i,j) = 1
               bh = height(n)
               btx(i,j  ) = bh
            end if
#endif
#ifdef MPI
            end if
#endif
         end do

#if 0
         do j = -1, ny+1
            do i = -1, nx+1
               if(ir(i,j) == 0) then
                  if(abs(btx(i,j)) > tiny(bh)) then
                     write(0,*) 'Somthing is wrong on btx!!!', btx(i,j)
                  end if
                  if(abs(bty(i,j-1)) > tiny(bh)) then
                     write(0,*) 'Somthing is wrong on bty!!!', bty(i,j-1)
                  end if
               else
                  write(400000,'(3i6,2f15.3)') i, j, ir(i,j), btx(i,j), bty(i,j-1)
               end if
            end do
         end do
#endif

#ifndef CARTESIAN
         deallocate(xin)
         deallocate(yin)
#endif
         deallocate(x)
         deallocate(y)
#ifdef OLDFORMAT
         deallocate(val)
#else
         deallocate(irread)
         deallocate(height)
#endif

      end if

      return
   end subroutine read_bank_file

   subroutine update_bathymetory(dg)
#ifdef MPI
      use mod_mpi, only : exchange_edges_dxbxdyby
#endif
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), pointer :: nx, ny
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy, dz
      integer(kind=4), pointer, dimension(:,:) :: ir
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: btx, bty, dxbx, dyby
      integer(kind=4) :: i, j

      if(dg%bank_file /= 'NO_BANK_FILE_GIVEN') then
         nx   => dg%my%nx
         ny   => dg%my%ny

         dx   => dg%depth_field%dx
         dy   => dg%depth_field%dy
         dz   => dg%depth_field%dz

         ir   => dg%wave_field%ir
         btx  => dg%wave_field%btx
         bty  => dg%wave_field%bty
         dxbx => dg%depth_field%dxbx
         dyby => dg%depth_field%dyby

         ! update bathymetory with vertical wall (the same as HMN)
!$omp parallel do private(i)
         do j = 1, ny
            do i = 1, nx
               if((ir(i,j) == 1) .or. (ir(i,j) == 3)) then
                  dxbx(i,j) = min(-btx(i,j), dz(i,j), dz(i+1,j))
               else
                  dxbx(i,j) = dx(i,j)
               end if
               if((ir(i,j+1) == 2) .or. (ir(i,j+1) == 3)) then
                  dyby(i,j) = min(-bty(i,j), dz(i,j), dz(i,j+1))
               else
                  dyby(i,j) = dy(i,j)
               end if
            end do
         end do
#ifdef MPI
         call exchange_edges_dxbxdyby(dg)
#endif
      end if

      return
   end subroutine update_bathymetory

   subroutine save_dxdy_old(dg)
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), pointer :: nx, ny
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx_old, dx, dy_old, dy
      integer(kind=4) :: i, j

      if(dg%bank_file /= 'NO_BANK_FILE_GIVEN') then
         nx     => dg%my%nx
         ny     => dg%my%ny

         dx_old => dg%depth_field%dx_old
         dx     => dg%depth_field%dx
         dy_old => dg%depth_field%dy_old
         dy     => dg%depth_field%dy

!$omp parallel do private(i)
         do j = 1, ny
            do i = 1, nx
               dx_old(i,j) = dx(i,j)
               dy_old(i,j) = dy(i,j)
            end do
         end do
      end if

      return
   end subroutine save_dxdy_old

   subroutine update_btxbty(dg)
#ifdef MPI
      use mod_mpi, only : exchange_edges_btxbty
#endif
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), pointer :: nx, ny
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx_old, dx, dy_old, dy
      integer(kind=4), pointer, dimension(:,:) :: ir
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: btx, bty
      integer(kind=4) :: i, j

      if(dg%bank_file /= 'NO_BANK_FILE_GIVEN') then
         nx     => dg%my%nx
         ny     => dg%my%ny

         dx_old => dg%depth_field%dx_old
         dx     => dg%depth_field%dx
         dy_old => dg%depth_field%dy_old
         dy     => dg%depth_field%dy

         ir     => dg%wave_field%ir
         btx    => dg%wave_field%btx
         bty    => dg%wave_field%bty

!$omp parallel do private(i)
         do j = 1, ny
            do i = 1, nx
               if((ir(i,j) == 1) .or. (ir(i,j) == 3)) then
                  btx(i,j) = btx(i,j) - (dx(i,j) - dx_old(i,j))
               end if
               if((ir(i,j+1) == 2) .or. (ir(i,j+1) == 3)) then
                  bty(i,j) = bty(i,j) - (dy(i,j) - dy_old(i,j))
               end if
            end do
         end do
#ifdef MPI
         call exchange_edges_btxbty(dg)
#endif
      end if

      return
   end subroutine update_btxbty
#else
#ifdef __SX__
   subroutine bank_dummy()
   end subroutine bank_dummy
#endif
#endif

end module mod_bank
