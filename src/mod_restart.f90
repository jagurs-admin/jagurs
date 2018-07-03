#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_restart
use mod_grid
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
use mod_params, only : with_elastic_loading
! ==============================================================================
#endif
! === DEBUG for restart with sinwave/file-input 2016/02/18 =====================
use mod_init_disp_sinwave, only : ccount
! ==============================================================================
! === Arrival time =============================================================
use mod_params, only : check_arrival_time
! ==============================================================================
implicit none

contains

   subroutine read_restart_file(ngrid, dgrid, fname)
      integer(kind=4), intent(in) :: ngrid
      type(data_grids), target, dimension(ngrid), intent(inout) :: dgrid
      character(len=256), intent(in) :: fname

      type(data_grids),pointer :: dg
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: rarray
      real(kind=REAL_BYTE), pointer, dimension(:) :: barray
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
      real(kind=8), pointer, dimension(:,:) :: darray
! ==============================================================================
#endif
      integer(kind=4), pointer, dimension(:,:) :: iarray
      integer(kind=4) :: ig

      open(1,file=trim(fname),action='read',status='old',form='unformatted')

! === DEBUG for restart with sinwave/file-input 2016/02/18 =====================
      read(1) ccount
! ==============================================================================
      do ig = 1, ngrid
         dg => dgrid(ig)

         rarray => dg%wave_field%fx
         read(1) rarray
         rarray => dg%wave_field%fy
         read(1) rarray
         rarray => dg%wave_field%hz
         read(1) rarray

         rarray => dg%depth_field%dx
         read(1) rarray
         rarray => dg%depth_field%dy
         read(1) rarray
         rarray => dg%depth_field%dz
         read(1) rarray

         rarray => dg%zz
         read(1) rarray
         rarray => dg%hzmax
         read(1) rarray
#ifdef HZMINOUT
         rarray => dg%hzmin
         read(1) rarray
#endif
#ifndef SKIP_MAX_VEL
         rarray => dg%vmax
         read(1) rarray
#endif
         iarray => dg%wod_flags
         read(1) iarray

         barray => dg%ubnd%north
         read(1) barray
         barray => dg%ubnd%east
         read(1) barray
         barray => dg%ubnd%south
         read(1) barray
         barray => dg%ubnd%west
         read(1) barray

         barray => dg%hbnd%north
         read(1) barray
         barray => dg%hbnd%east
         read(1) barray
         barray => dg%hbnd%south
         read(1) barray
         barray => dg%hbnd%west
         read(1) barray
#ifdef BANKFILE
         rarray => dg%wave_field%btx
         read(1) rarray
         rarray => dg%wave_field%bty
         read(1) rarray
         rarray => dg%depth_field%dxbx
         read(1) rarray
         rarray => dg%depth_field%dyby
         read(1) rarray
         iarray => dg%wave_field%brokenx
         read(1) iarray
         iarray => dg%wave_field%brokeny
         read(1) iarray
#endif

#ifndef CARTESIAN
! === Elastic Loading ==========================================================
         if(with_elastic_loading == 1) then
            darray => dg%loading%defZmap1
            read(1) darray
         end if
! ==============================================================================
#endif
! === Arrival time =============================================================
         if(check_arrival_time == 1) then
            iarray => dg%wave_field%arrivedat
            read(1) iarray
         end if
! ==============================================================================
      end do

      close(1)

      return
   end subroutine read_restart_file

   subroutine write_restart_file(ngrid, dgrid, fname)
      integer(kind=4), intent(in) :: ngrid
      type(data_grids), target, dimension(ngrid), intent(inout) :: dgrid
      character(len=256), intent(in) :: fname

      type(data_grids),pointer :: dg
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: rarray
      real(kind=REAL_BYTE), pointer, dimension(:) :: barray
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
      real(kind=8), pointer, dimension(:,:) :: darray
! ==============================================================================
#endif
      integer(kind=4), pointer, dimension(:,:) :: iarray
      integer(kind=4) :: ig

      open(1,file=trim(fname),action='write',status='replace',form='unformatted')

! === DEBUG for restart with sinwave/file-input 2016/02/18 =====================
      write(1) ccount
! ==============================================================================
      do ig = 1, ngrid
         dg => dgrid(ig)

         rarray => dg%wave_field%fx
         write(1) rarray
         rarray => dg%wave_field%fy
         write(1) rarray
         rarray => dg%wave_field%hz
         write(1) rarray

         rarray => dg%depth_field%dx
         write(1) rarray
         rarray => dg%depth_field%dy
         write(1) rarray
         rarray => dg%depth_field%dz
         write(1) rarray

         rarray => dg%zz
         write(1) rarray
         rarray => dg%hzmax
         write(1) rarray
#ifdef HZMINOUT
         rarray => dg%hzmin
         write(1) rarray
#endif
#ifndef SKIP_MAX_VEL
         rarray => dg%vmax
         write(1) rarray
#endif
         iarray => dg%wod_flags
         write(1) iarray

         barray => dg%ubnd%north
         write(1) barray
         barray => dg%ubnd%east
         write(1) barray
         barray => dg%ubnd%south
         write(1) barray
         barray => dg%ubnd%west
         write(1) barray

         barray => dg%hbnd%north
         write(1) barray
         barray => dg%hbnd%east
         write(1) barray
         barray => dg%hbnd%south
         write(1) barray
         barray => dg%hbnd%west
         write(1) barray
#ifdef BANKFILE
         rarray => dg%wave_field%btx
         write(1) rarray
         rarray => dg%wave_field%bty
         write(1) rarray
         rarray => dg%depth_field%dxbx
         write(1) rarray
         rarray => dg%depth_field%dyby
         write(1) rarray
         iarray => dg%wave_field%brokenx
         write(1) iarray
         iarray => dg%wave_field%brokeny
         write(1) iarray
#endif

#ifndef CARTESIAN
! === Elastic Loading ==========================================================
         if(with_elastic_loading == 1) then
            darray => dg%loading%defZmap1
            write(1) darray
         end if
! ==============================================================================
#endif
! === Arrival time =============================================================
         if(check_arrival_time == 1) then
            iarray => dg%wave_field%arrivedat
            write(1) iarray
         end if
! ==============================================================================
      end do

      close(1)

      return
   end subroutine write_restart_file
end module mod_restart
