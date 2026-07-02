#ifdef TIMER_DETAIL
#   define TIMER_START(a) call start_timer(a)
#   define TIMER_STOP(a)  call stop_timer(a)
#else
#   define TIMER_START(a)
#   define TIMER_STOP(a)
#endif
#ifdef DBLE_MATH
#include "dble_math.h"
#endif
#include "real.h"
module mod_rwg
use mod_grid
! === For negative max. height =================================================
!use mod_params, only : VEL, HGT
! === Arrival time =============================================================
!use mod_params, only : VEL, HGT, missing_value
use mod_params, only : VEL, HGT, missing_value, check_arrival_time, check_arrival_height
! ==============================================================================
use mod_params, only : check_tt_time, check_ttt_height
! ==============================================================================
use mod_mygmt_gridio, only : read_gmt_grd
#ifdef MPI
use mod_mpi
#endif
#ifdef TIMER_DETAIL
use mod_timer
#endif
use mod_params, only : with_disp, with_abc, nxa, nya, apara
#ifndef CARTESIAN
use mod_hxy
use mod_fxy
use mod_fxy_disp
#else
use mod_hxy_cartesian
use mod_fxy_cartesian
use mod_fxy_disp_cartesian
#endif
use mod_nest
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
use mod_loading, only : loading_run
! === Elastic loading with interpolation =======================================
!use mod_params, only : with_elastic_loading
use mod_params, only : with_elastic_loading, elastic_loading_interpolation
! ==============================================================================
! ==============================================================================
! === Coriolis force is supported on linear calc. ==============================
use mod_fxy_Coriolis
use mod_fxy_Coriolis_disp
! ==============================================================================
#endif
#if defined(MPI) && defined(ONEFILE)
use mod_onefile, only : onefile_scatter_array
#endif
! === Elastic loading with interpolation =======================================
#ifndef CARTESIAN
use mod_interpolation, only : interp2fine_elastic_loading
#endif
! ==============================================================================
implicit none

contains

   subroutine maxgrd_init_rwg(hzmax,nlon,nlat)
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(out) :: hzmax
      integer(kind=4), intent(in) :: nlon, nlat
      integer :: i, j
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
! === For negative max. height =================================================
!           hzmax(i,j) = 0.0d0
            hzmax(i,j) = missing_value
! ==============================================================================
         end do
      end do
      return
   end subroutine maxgrd_init_rwg
#ifdef HZMINOUT
   subroutine mingrd_init_rwg(hzmin,nlon,nlat)
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(out) :: hzmin
      integer(kind=4), intent(in) :: nlon, nlat
      integer :: i, j
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
! === For negative min. height =================================================
!           hzmin(i,j) = 0.0d0
            hzmin(i,j) = missing_value
! ==============================================================================
         end do
      end do
      return
   end subroutine mingrd_init_rwg
#endif

   subroutine maxgrd_check_nl(hzmax,wfld,wod,nlon,nlat)
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(inout) :: hzmax
      type(wave_arrays), target, intent(in) :: wfld
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(in) :: wod
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: nlon, nlat

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz
      integer(kind=4) :: i, j

      hz => wfld%hz

      !** if wet check for hzmax **
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(wod(i,j) == 1) then
               if(hz(i,j) > hzmax(i,j)) hzmax(i,j) = hz(i,j)
            end if
         end do
      end do

      return
   end subroutine maxgrd_check_nl
#ifdef HZMINOUT
   subroutine mingrd_check_nl(hzmin,wfld,wod,nlon,nlat)
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(inout) :: hzmin
      type(wave_arrays), target, intent(in) :: wfld
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(in) :: wod
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: nlon, nlat

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz
      integer(kind=4) :: i, j

      hz => wfld%hz

      !** if wet check for hzmin **
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(wod(i,j) == 1) then
               if((hzmin(i,j) < missing_value*0.9d0) .or. (hz(i,j) < hzmin(i,j))) then
                  hzmin(i,j) = hz(i,j)
               end if
            end if
         end do
      end do

      return
   end subroutine mingrd_check_nl
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
   subroutine maxgrd_v_init_rwg(vmax,nlon,nlat)
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(out) :: vmax
      integer(kind=4), intent(in) :: nlon, nlat
      integer :: i, j
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
            vmax(i,j) = 0.0d0
         end do
      end do
      return
   end subroutine maxgrd_v_init_rwg

! === Conversion from flux to velocity should be done right after calc. ========
!  subroutine maxgrd_v_check_nl(vmax,wfld,dfld,wod,nlon,nlat)
#ifndef MPI
   subroutine maxgrd_v_check_nl(vmax,wfld,dfld,wod,nlon,nlat,linear_flag)
#else
   subroutine maxgrd_v_check_nl(vmax,wfld,dfld,wod,nlon,nlat,bflag,linear_flag)
#endif
! ==============================================================================
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(inout) :: vmax
      type(wave_arrays), target, intent(in) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(in) :: wod
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: nlon, nlat
! === Conversion from flux to velocity should be done right after calc. ========
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: linear_flag

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: vel
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz
! === Conversion from flux to velocity should be done right after calc. ========
!     real(kind=REAL_BYTE) :: td, tx, ty
      real(kind=REAL_BYTE) :: tx, ty
! ==============================================================================
      real(kind=REAL_BYTE), parameter :: td_min = 0.01d0 ! 1 cm
! === Conversion from flux to velocity should be done right after calc. ========
      integer(kind=4) :: im, jm, ip, jp
      real(kind=REAL_BYTE) :: tdxm, tdxp, tdym, tdyp
! ==============================================================================

      fx => wfld%fx
      fy => wfld%fy
      hz => wfld%hz
      dz => dfld%dz

      !** if wet check for vmax **
! === Conversion from flux to velocity should be done right after calc. ========
      if(linear_flag == 0) then
#ifndef USE_GPU
!$omp parallel do private(i, im, ip, jm, jp, tdxm, tdxp, tdym, tdyp, tx, ty, vel)
#else
!$omp target teams distribute parallel do collapse(2) private(i,im,ip,jm,jp,tdxm,tdxp,tdym,tdyp,tx,ty,vel)
#endif
         do j = 1, nlat
            do i = 1, nlon
               tx = 0.0d0
#ifndef MPI
               im = max(1,   i-1)
               ip = min(nlon,i+1)
#else
               im = i - 1
               ip = i + 1
               if(iand(bflag, WEST_BOUND) /= 0) im = max(1,im)
               if(iand(bflag, EAST_BOUND) /= 0) ip = min(nlon,ip)
#endif
               tdxm = 0.5d0*(dz(i,j) + hz(i,j) + dz(im,j) + hz(im,j))
               tdxp = 0.5d0*(dz(i,j) + hz(i,j) + dz(ip,j) + hz(ip,j))
               if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                  wod(ip,j) == 1 .and. tdxp > td_min .and. &
                  wod(i,j) == 1) then
                  tx = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
               end if

               ty = 0.0d0
#ifndef MPI
               jm = max(1,   j-1)
               jp = min(nlat,j+1)
#else
               jm = j - 1
               jp = j + 1
               if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,jm)
               if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(nlat,jp)
#endif
               tdym = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jm) + hz(i,jm))
               tdyp = 0.5d0*(dz(i,j) + hz(i,j) + dz(i,jp) + hz(i,jp))
               if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                  wod(i,jp) == 1 .and. tdyp > td_min .and. &
                  wod(i,j) == 1) then
                  ty = 0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
               end if

               vel = sqrt(tx**2 + ty**2)
               if(vel > vmax(i,j)) vmax(i,j) = vel
            end do
         end do
      else
#ifndef USE_GPU
!$omp parallel do private(i, im, ip, jm, jp, tdxm, tdxp, tdym, tdyp, tx, ty, vel)
#else
!$omp target teams distribute parallel do collapse(2) private(i,im,ip,jm,jp,tdxm,tdxp,tdym,tdyp,tx,ty,vel)
#endif
         do j = 1, nlat
            do i = 1, nlon
               tx = 0.0d0
#ifndef MPI
               im = max(1,   i-1)
               ip = min(nlon,i+1)
#else
               im = i - 1
               ip = i + 1
               if(iand(bflag, WEST_BOUND) /= 0) im = max(1,im)
               if(iand(bflag, EAST_BOUND) /= 0) ip = min(nlon,ip)
#endif
               tdxm = 0.5d0*(dz(i,j) + dz(im,j))
               tdxp = 0.5d0*(dz(i,j) + dz(ip,j))
               if(wod(im,j) == 1 .and. tdxm > td_min .and. &
                  wod(ip,j) == 1 .and. tdxp > td_min .and. &
                  wod(i,j) == 1) then
                  tx = 0.5d0*(fx(i,j)/tdxp + fx(im,j)/tdxm)
               end if

               ty = 0.0d0
#ifndef MPI
               jm = max(1,   j-1)
               jp = min(nlat,j+1)
#else
               jm = j - 1
               jp = j + 1
               if(iand(bflag, NORTH_BOUND) /= 0) jm = max(1,jm)
               if(iand(bflag, SOUTH_BOUND) /= 0) jp = min(nlat,jp)
#endif
               tdym = 0.5d0*(dz(i,j) + dz(i,jm))
               tdyp = 0.5d0*(dz(i,j) + dz(i,jp))
               if(wod(i,jm) == 1 .and. tdym > td_min .and. &
                  wod(i,jp) == 1 .and. tdyp > td_min .and. &
                  wod(i,j) == 1) then
                  ty = 0.5d0*(fy(i,j)/tdyp + fy(i,jm)/tdym)
               end if

               vel = sqrt(tx**2 + ty**2)
               if(vel > vmax(i,j)) vmax(i,j) = vel
            end do
         end do
      end if
! ==============================================================================

      return
   end subroutine maxgrd_v_check_nl
! ==============================================================================

#ifndef MPI
   subroutine outsea_rwg(wfld,dfld,hbnd,ubnd,nlon,nlat)
#else
   subroutine outsea_rwg(wfld,dfld,hbnd,ubnd,nlon,nlat,bflag)
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(inout) :: dfld
      type(boundary_arrays), target, intent(inout) :: hbnd, ubnd
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
#endif

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz
      real(kind=REAL_BYTE), pointer, dimension(:) :: hb, ub
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE), parameter :: zap = 0.0d0

      if(timenest /= 1) then
         hz => wfld%hz
      else
         hz => wfld%hz_b
      end if
      dz => dfld%dz

      hb => hbnd%north
      ub => ubnd%north
#ifndef USE_GPU
!$omp parallel
#endif
#ifdef MPI
      if(iand(bflag, NORTH_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
      do i = 1, nlon-1
         if(dz(i,1) <= zap) then
            hz(i,1) = zap
         else
            hz(i,1) = hb(i) + ub(i)*(hz(i,1)-hz(i,2))
         end if
         hb(i) = hz(i,2)
      end do
#ifdef MPI
      end if
#endif

#ifndef USE_GPU
!$omp single
#endif
      hb => hbnd%east
      ub => ubnd%east
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, EAST_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
      do j = 1, nlat-1
         if(dz(nlon,j) <= zap) then
            hz(nlon,j) = zap
         else
            hz(nlon,j) = hb(j) + ub(j)*(hz(nlon,j)-hz(nlon-1,j))
         end if
         hb(j) = hz(nlon-1,j)
      end do
#ifdef MPI
      end if
#endif

#ifndef USE_GPU
!$omp single
#endif
      hb => hbnd%south
      ub => ubnd%south
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, SOUTH_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
      do i = 1, nlon-1
         if(dz(nlon-i+1,nlat) <= zap) then
            hz(nlon-i+1,nlat) = zap
         else
            hz(nlon-i+1,nlat) = hb(i) + ub(i)*(hz(nlon-i+1,nlat)-hz(nlon-i+1,nlat-1))
         end if
         hb(i) = hz(nlon-i+1,nlat-1)
      end do
#ifdef MPI
      end if
#endif

#ifndef USE_GPU
!$omp single
#endif
      hb => hbnd%west
      ub => ubnd%west
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, WEST_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
      do j = 1, nlat-1
         if(dz(1,nlat-j+1) <= zap) then
            hz(1,nlat-j+1) = zap
         else
            hz(1,nlat-j+1) = hb(j) + ub(j)*(hz(1,nlat-j+1)-hz(2,nlat-j+1))
         end if
         hb(j) = hz(2,nlat-j+1)
      end do
#ifdef MPI
      end if
#endif
#ifndef USE_GPU
!$omp end parallel
#endif

      return
   end subroutine outsea_rwg

   subroutine tstep_grid(mode,ig,cg,fg,cf,cfl,crls,dt,smallh_xy,smallh_wod,c2p_all,conv_step,istep)
      integer(kind=4), intent(in) :: mode, ig
      type(data_grids), target, intent(inout) :: cg, fg
      real(kind=REAL_BYTE), intent(in) :: cf
! === Modification to fit pointer version! =====================================
      real(kind=REAL_BYTE), intent(in) :: cfl
! ==============================================================================
      integer(kind=4), intent(in) :: crls
      real(kind=REAL_BYTE), intent(in) :: dt, smallh_xy, smallh_wod
      integer(kind=4), intent(in) :: c2p_all
      integer(kind=4), intent(out) :: conv_step
      integer(kind=4), intent(in) :: istep

      type(wave_arrays), pointer :: wfld
      type(depth_arrays), pointer :: dfld
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: ffld
      real(kind=REAL_BYTE) :: th0, dth
      integer(kind=4) :: niz, njz, lflag
      integer(kind=4), pointer, dimension(:,:) :: wod
      type(boundary_arrays), pointer :: hbnd, ubnd
      integer(kind=4) :: rank, joff, bflag
      integer(kind=4) :: calchgt, numneststeps, neststephgt
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz_old, hz_new, hz
      real(kind=REAL_BYTE) :: rate
      integer(kind=4) :: i, j

      ! set pointers for function calls
      wfld   => fg%wave_field
      dfld   => fg%depth_field
      lflag  =  fg%my%linear_flag
      niz    =  fg%my%nx
      njz    =  fg%my%ny
#ifndef CARTESIAN
      th0    =  fg%my%th0
#else
      th0    =  fg%my%dh
#endif
      dth    =  fg%my%dth
      wod    => fg%wod_flags
      ffld   => fg%bcf_field
      hbnd   => fg%hbnd
      ubnd   => fg%ubnd
#ifdef MPI
      rank  = fg%my%ry*fg%my%px+fg%my%rx
      bflag = fg%my%has_boundary
      if(fg%my%iy == 1) then
         joff = fg%my%iy - 1
      else
         joff = fg%my%iy - 2 ! sub. edge.
      end if
#else
      rank  = 0
      bflag = 0
      joff  = 0
#endif
      if(timenest == 1) conv_step = 0

      if(mode == VEL) then
         if(lflag == 1) then
            TIMER_START('- fxy_rwg')
            if((with_disp == 0) .or. (with_disp == 2 .and. ig == 1)) then
#ifndef CARTESIAN
! === Coriolis force is supported on linear calc. ==============================
               if(crls == 0) then
! ==============================================================================
#endif
               call fxy_rwg(wfld,dfld,dt,th0,dth,joff,niz,njz,istep,fg)
#ifndef CARTESIAN
! === Coriolis force is supported on linear calc. ==============================
               else
                  call fxy_rwg_Coriolis(wfld,dfld,wod,crls,dt,th0,dth,joff,niz,njz,ig,bflag,fg)
               end if
! ==============================================================================
#endif
            else
#ifndef CARTESIAN
! === Coriolis force is supported on linear calc. ==============================
               if(crls == 0) then
! ==============================================================================
#endif
               call fxy_rwg_disp(wfld,dfld,dt,th0,dth,joff,niz,njz,ig,bflag,fg,cg,conv_step,istep)
#ifndef CARTESIAN
! === Coriolis force is supported on linear calc. ==============================
               else
                  call fxy_rwg_Coriolis_disp(wfld,dfld,wod,crls,dt,th0,dth,joff,niz,njz,ig,bflag,fg,cg,conv_step)
               end if
! ==============================================================================
#endif
            end if
            TIMER_STOP('- fxy_rwg')
         else
            TIMER_START('- fxynl_rwg')
            if((with_disp == 0) .or. (with_disp == 2 .and. ig == 1)) then
               call fxynl_rwg(wfld,dfld,ffld,wod,cf,cfl,crls,dt,th0,dth,joff,niz,njz,ig,smallh_xy,bflag,istep,fg)
            else
               call fxynl_rwg_disp(wfld,dfld,ffld,wod,cf,cfl,crls,dt,th0,dth,joff,niz,njz, &
                                   ig,smallh_xy,bflag,fg,cg,conv_step,istep)
            end if
#ifdef MPI
! === Flood Change =============================================================
            if(timenest /= 1) then
               call exchange_edges(HGT,fg)
            else
               call exchange_edges_b(HGT,fg)
            end if
! ==============================================================================
#endif
            TIMER_STOP('- fxynl_rwg')
         end if
         if(with_abc == 1) then
            if(timenest /= 1) then
               if(ig == 1) call apply_abc(mode,wfld,niz,njz)
            else
               if(ig == 1 .and. fg%my%calcvel == 1) call apply_abc(mode,wfld,niz,njz)
            end if
         end if
#ifdef MPI
! === For 9-point averaging on copy2coarse, exchange_edges is necessary! =======
         if(c2p_all == 1) then
            TIMER_START('- exch_edges_vel')
            call exchange_edges(mode,fg)
            if(timenest == 1) call exchange_edges_b(mode,fg)
            TIMER_STOP('- exch_edges_vel')
         end if
! ==============================================================================
         TIMER_START('- mapgrids_vel')
         if(ig /= 1) call mapgrids(mode,cg,fg,c2p_all)
         TIMER_STOP('- mapgrids_vel')
! === DEBUG by tkato 2013/10/10 ================================================
         if(ig /= 1) call exchange_edges(mode,cg)
         if((timenest == 1) .and. (ig /= 1)) then
            call exchange_edges_b(mode,cg)
            call exchange_edges_b(mode,fg)
         end if
! ==============================================================================
         TIMER_START('- exch_edges_vel')
         call exchange_edges(mode,fg)
         if(timenest == 1) call exchange_edges_b(mode,fg)
         TIMER_STOP('- exch_edges_vel')
! === Modification to fit pointer version! =====================================
         TIMER_START('- exch_edges_wod')
         call exchange_edges_wod(fg)
         TIMER_STOP('- exch_edges_wod')
! ==============================================================================
#endif
      else if(mode == HGT) then
         if(lflag == 1) then
            TIMER_START('- hxy_rwg')
            call hxy_rwg(wfld,dt,th0,dth,joff,niz,njz,fg)
            TIMER_STOP('- hxy_rwg')
         else
            TIMER_START('- hxynl_rwg')
            call hxynl_rwg(wfld,dfld,wod,dt,th0,dth,joff,niz,njz,smallh_xy,bflag,fg)
            TIMER_STOP('- hxynl_rwg')
         end if
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
         if(with_elastic_loading == 1) then
! === Elastic loading with interpolation =======================================
         if(elastic_loading_interpolation == 0) then
! ==============================================================================
            if((timenest == 0) .or. (fg%my%calchgt == 1)) then
               TIMER_START('- loading_run')
               call loading_run(fg)
               TIMER_STOP('- loading_run')
            end if
! === Elastic loading with interpolation =======================================
         else
            if(ig == 1) then
               if((timenest == 0) .or. (fg%my%calchgt == 1)) then
                  TIMER_START('- loading_run')
                  call loading_run(fg)
                  TIMER_STOP('- loading_run')
               end if
            else
               TIMER_START('- interp2fine_elastic_loading')
               call interp2fine_elastic_loading(cg,fg)
               TIMER_STOP('- interp2fine_elastic_loading')
            end if
         end if
! ==============================================================================
         end if
! ==============================================================================
#endif
#ifdef MPI
! === For 9-point averaging on copy2coarse, exchange_edges is necessary! =======
         if(c2p_all == 1) then
            TIMER_START('- exch_edges_hgt')
            call exchange_edges(mode,fg)
            if(timenest == 1) call exchange_edges_b(mode,fg)
            TIMER_STOP('- exch_edges_hgt')
! === DEBUG by tkato 2013/10/10 ================================================
            if(lflag == 0) call recheck_wod(wfld,dfld,wod,niz,njz,smallh_xy)
! ==============================================================================
         end if
! ==============================================================================
         TIMER_START('- mapgrids_hgt')
         if(ig /= 1) call mapgrids(mode,cg,fg,c2p_all)
         TIMER_STOP('- mapgrids_hgt')
! === DEBUG by tkato 2013/10/10 ================================================
         if(ig /= 1) call exchange_edges(mode,cg)
         if((timenest == 1) .and. (ig /= 1)) then
            call exchange_edges_b(mode,cg)
            call exchange_edges_b(mode,fg)
         end if
! ==============================================================================
! === recheck_wod should be called after outsea_rwg. by tkato 2012/09/11 =======
!        if(lflag == 0) call recheck_wod(wfld,dfld,wod,niz,njz,smallh_wod)
!        if(ig == 1) call outsea_rwg(wfld,dfld,hbnd,ubnd,niz,njz,bflag)
!        call exchange_edges(mode,fg)
!        call exchange_edges_wod(fg)
         if(timenest /= 1) then
            TIMER_START('- outsea_rwg')
            if(ig == 1) call outsea_rwg(wfld,dfld,hbnd,ubnd,niz,njz,bflag)
            TIMER_STOP('- outsea_rwg')
            TIMER_START('- exch_edges_hgt')
            if(with_abc == 1) then
               if(ig == 1) call apply_abc(mode,wfld,niz,njz)
            end if
            call exchange_edges(mode,fg)
            TIMER_STOP('- exch_edges_hgt')
         else
            calchgt      = fg%my%calchgt
            numneststeps = fg%my%numneststeps
            neststephgt  = fg%my%neststephgt

            TIMER_START('- outsea_rwg')
            if(ig == 1 .and. calchgt == 1) call outsea_rwg(wfld,dfld,hbnd,ubnd,niz,njz,bflag)
            TIMER_STOP('- outsea_rwg')
            TIMER_START('- exch_edges_hgt')
            if(with_abc == 1) then
               if(ig == 1 .and. calchgt == 1) call apply_abc(mode,wfld,niz,njz)
            end if
            call exchange_edges_b(mode,fg)
            TIMER_STOP('- exch_edges_hgt')

            if(neststephgt == numneststeps - 1) then
               hz     => wfld%hz
               hz_old => wfld%hz_a
               hz_new => wfld%hz_b
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
               do j = lbound(hz,2), ubound(hz,2)
                  do i = lbound(hz,1), ubound(hz,1)
                     hz(i,j) = hz_new(i,j)
                     hz_old(i,j) = hz_new(i,j)
                  end do
               end do
            else
               hz     => wfld%hz
               hz_old => wfld%hz_a
               hz_new => wfld%hz_b
               call exchange_edges_b(mode,fg)
               rate = dble(neststephgt + 1)/dble(numneststeps)
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
               do j = lbound(hz,2), ubound(hz,2)
                  do i = lbound(hz,1), ubound(hz,1)
                     hz(i,j) = hz_old(i,j)*(1.0d0 - rate) + hz_new(i,j)*rate
                  end do
               end do
            end if
         end if
         TIMER_START('- recheck_wod')
! === Flood Change =============================================================
!        if(lflag == 0) call recheck_wod(wfld,dfld,wod,niz,njz,smallh_wod)
         if(lflag == 0) then
            call recheck_wod(wfld,dfld,wod,niz,njz,smallh_wod)
            call exchange_edges(mode,fg)
         end if
! ==============================================================================
         TIMER_STOP('- recheck_wod')
! ==============================================================================
#endif
      end if

#ifndef MPI
      if(timenest == 1 .and. mode == HGT .and. lflag == 0) call recheck_wod(wfld,dfld,wod,niz,njz,smallh_xy)
      TIMER_START('- mapgrids')
      if(ig > 1) call mapgrids(mode,cg,fg,c2p_all)
      TIMER_STOP('- mapgrids')

! === recheck_wod should be called after outsea_rwg. by tkato 2012/09/11 =======
      if(timenest /= 1) then
         TIMER_START('- outsea_rwg')
         if(mode == HGT .and. ig == 1) call outsea_rwg(wfld,dfld,hbnd,ubnd,niz,njz)
         TIMER_STOP('- outsea_rwg')
! ==============================================================================
         if(with_abc == 1) then
            if(mode == HGT .and. ig == 1) call apply_abc(mode,wfld,niz,njz)
         end if
      else
         calchgt      = fg%my%calchgt
         numneststeps = fg%my%numneststeps
         neststephgt  = fg%my%neststephgt

         TIMER_START('- outsea_rwg')
         if(mode == HGT .and. ig == 1 .and. calchgt == 1) call outsea_rwg(wfld,dfld,hbnd,ubnd,niz,njz)
         TIMER_STOP('- outsea_rwg')

         if(with_abc == 1) then
            if(mode == HGT .and. ig == 1 .and. calchgt == 1) call apply_abc(mode,wfld,niz,njz)
         end if

         if(mode == HGT) then
            if(neststephgt == numneststeps - 1) then
               hz     => wfld%hz
               hz_old => wfld%hz_a
               hz_new => wfld%hz_b
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
               do j = lbound(hz,2), ubound(hz,2)
                  do i = lbound(hz,1), ubound(hz,1)
                     hz(i,j) = hz_new(i,j)
                     hz_old(i,j) = hz_new(i,j)
                  end do
               end do
            else
               hz     => wfld%hz
               hz_old => wfld%hz_a
               hz_new => wfld%hz_b
               rate = dble(neststephgt + 1)/dble(numneststeps)
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
               do j = lbound(hz,2), ubound(hz,2)
                  do i = lbound(hz,1), ubound(hz,1)
                     hz(i,j) = hz_old(i,j)*(1.0d0 - rate) + hz_new(i,j)*rate
                  end do
               end do
            end if
         end if
      end if

      TIMER_START('- recheck_wod')
      if(mode == HGT .and. lflag == 0) call recheck_wod(wfld,dfld,wod,niz,njz,smallh_wod)
      TIMER_STOP('- recheck_wod')
#endif

      return
   end subroutine tstep_grid

   subroutine recheck_wod(wfld,dfld,ifz,nlon,nlat,smallh)
      type(wave_arrays), target, intent(in) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(inout) :: ifz
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(inout) :: ifz
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(inout) :: ifz
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: smallh

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz
      integer(kind=4) :: i, j

      hz => wfld%hz
      dz => dfld%dz

#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(dz(i,j) + hz(i,j) > smallh) then
               ifz(i,j) = 1
            else
! === Flood Change =============================================================
               if(ifz(i,j) == 1) hz(i,j) = -dz(i,j)
! ==============================================================================
               ifz(i,j) = -1
            end if
         end do
      end do

      return
   end subroutine recheck_wod

   subroutine drise_rwg(dfld,zz,dt,tau,nlon,nlat,linear,defbathy_flag)
      type(depth_arrays), target, intent(inout) :: dfld
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(in) :: zz
      real(kind=REAL_BYTE), intent(in) :: dt, tau
      integer(kind=4), intent(in) :: nlon, nlat, linear, defbathy_flag

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy, dz
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: dev

      real(kind=REAL_BYTE), parameter :: zap = 0.0d0
      real(kind=REAL_BYTE), parameter :: half = 0.5d0

      dx => dfld%dx
      dy => dfld%dy
      dz => dfld%dz

      if(tau > 0.0d0) then
         dev = dt/tau
      else
         dev = 1.0d0
      end if

#ifndef USE_GPU
!$omp parallel
#endif
      ! Burbidge Only deform bathymetry if defbathy_flag = 1
      if(defbathy_flag == 1) then
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon
               dz(i,j) = dz(i,j) - dev*zz(i,j)
            end do
         end do
      end if

      if(linear == 1) then
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon-1
               if(dz(i+1,j) < zap .or. dz(i,j) < zap) then
                  dx(i,j) = zap
               else
                  dx(i,j) = 0.5d0*(dz(i+1,j) + dz(i,j))
               end if
            end do
#ifdef USE_GPU
         end do
!$omp target teams distribute parallel do
         do j = 1, nlat
#endif
            if(dz(nlon,j) < zap) then
               dx(nlon,j) = zap
            else
               dx(nlon,j) = dz(nlon,j)
            end if
         end do

#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat-1
            do i = 1, nlon
               if(dz(i,j+1) < zap .or. dz(i,j) < zap) then
                  dy(i,j) = zap
               else
                  dy(i,j) = 0.5d0*(dz(i,j+1) + dz(i,j))
               end if
            end do
         end do

#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do i = 1, nlon
            if(dz(i,nlat) < zap) then
               dy(i,nlat) = zap
            else
               dy(i,nlat) = dz(i,nlat)
            end if
         end do
      else
#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon-1
               dx(i,j) = half*(dz(i+1,j) + dz(i,j))
            end do
#ifdef USE_GPU
         end do
!$omp target teams distribute parallel do
         do j = 1, nlat
#endif
            dx(nlon,j) = dz(nlon,j)
         end do

#ifndef USE_GPU
!$omp do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat-1
            do i = 1, nlon
               dy(i,j) = half*(dz(i,j+1) + dz(i,j))
            end do
         end do

#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do i = 1, nlon
            dy(i,nlat) = dz(i,nlat)
         end do
      end if
#ifndef USE_GPU
!$omp end parallel
#endif

      return
   end subroutine drise_rwg

! === When "def_bathy=0", hz is changed on dry cell and it can become wet. =====
!  subroutine hrise_rwg(wfld,zz,dt,tau,nlon,nlat)
   subroutine hrise_rwg(wfld,zz,dt,tau,nlon,nlat,wod,defbathy_flag)
! ==============================================================================
      type(wave_arrays), target, intent(inout) :: wfld
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(in) :: zz
      real(kind=REAL_BYTE), intent(in) :: dt, tau
      integer(kind=4), intent(in) :: nlon, nlat
! === When "def_bathy=0", hz is changed on dry cell and it can become wet. =====
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(in) :: wod
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(in) :: wod
#endif
      integer(kind=4), intent(in) :: defbathy_flag
! ==============================================================================

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz
      integer(kind=4) :: i, j
      real(kind=REAL_BYTE) :: dev

      if(timenest /= 1) then
         hz => wfld%hz
      else
         hz => wfld%hz_b
      end if

      if(tau > 0.0d0) then
         dev = dt/tau
      else
         dev = 1.0d0
      end if

! === When "def_bathy=0", hz is changed on dry cell and it can become wet. =====
      if(defbathy_flag == 1) then
! ==============================================================================
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
            hz(i,j) = hz(i,j) + dev*zz(i,j)
         end do
      end do
! === When "def_bathy=0", hz is changed on dry cell and it can become wet. =====
      else
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon
               if(wod(i,j) == 1) hz(i,j) = hz(i,j) + dev*zz(i,j)
            end do
         end do
      end if
! ==============================================================================

      return
   end subroutine hrise_rwg

#ifndef CARTESIAN
#ifndef MPI
   subroutine boundary_rwg(dfld,hbnd,ubnd,dt,th0,dth,nlon,nlat)
#else
   subroutine boundary_rwg(dfld,hbnd,ubnd,dt,th0,dth,joff,nlon,nlat,bflag)
#endif
      type(depth_arrays), target, intent(in) :: dfld
      type(boundary_arrays), target, intent(inout) :: hbnd, ubnd
      real(kind=REAL_BYTE), intent(in) :: dt, th0, dth
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: joff, bflag
#endif

      real(kind=REAL_BYTE) :: sint, theta
      real(kind=REAL_BYTE) :: g, rote
      real(kind=REAL_BYTE), pointer, dimension(:) :: ub
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy
      integer(kind=4) :: ib, jb
      real(kind=REAL_BYTE) :: dtds, cfac
      real(kind=REAL_BYTE) :: zap = 0.0d0

      hbnd%north = 0.0d0
      hbnd%east  = 0.0d0
      hbnd%south = 0.0d0
      hbnd%west  = 0.0d0

      dx => dfld%dx
      dy => dfld%dy

      rote=6371.0d0 * 1000.0d0 !*** meters ***
      g = 9.8d0  !*** m/s^2 ***

      dtds = dt*sqrt(g)/(dth*rote)

      ! north
      ub => ubnd%north
#ifndef USE_GPU
!$omp parallel
#endif
#ifdef MPI
      if(iand(bflag, NORTH_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(cfac)
#else
!$omp target teams distribute parallel do private(cfac)
#endif
      do ib = 1, nlon
         if(dy(ib,1) > zap) then
            cfac = dtds*sqrt(dy(ib,1))
         else
            cfac = zap
         end if
         ub(ib) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do ib = 1, nlon
            ub(ib) = 1.0d0
         end do
      end if
#endif

      ! east
#ifndef USE_GPU
!$omp single
#endif
      ub => ubnd%east
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, EAST_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(theta, sint, cfac)
#else
!$omp target teams distribute parallel do private(theta,sint,cfac)
#endif
      do jb = 1, nlat
#ifndef MPI
         theta = th0 + (jb-1)*dth
#else
         theta = th0 + (jb+joff-1)*dth
#endif
         sint = sin(theta)
         if(dx(nlon-1,jb) > zap) then
            cfac = dtds*sqrt(dx(nlon-1,jb))/sint
         else
            cfac = zap
         end if
         ub(jb) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do jb = 1, nlat
            ub(jb) = 1.0d0
         end do
      end if
#endif

      ! south; be careful, here we go backwards across bottom of grid
#ifndef USE_GPU
!$omp single
#endif
      ub => ubnd%south
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, SOUTH_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(cfac)
#else
!$omp target teams distribute parallel do private(cfac)
#endif
      do ib = 1, nlon
         if(dy(nlon-ib+1,nlat-1) > zap) then
            cfac = dtds*sqrt(dy(nlon-ib+1,nlat-1))
         else
            cfac = zap
         end if
         ub(ib) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do ib = 1, nlon
            ub(ib) = 1.0d0
         end do
      end if
#endif

      ! west; be careful, here we go backwards up the left of grid
#ifndef USE_GPU
!$omp single
#endif
      ub => ubnd%west
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, WEST_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(theta, sint, cfac)
#else
!$omp target teams distribute parallel do private(theta,sint,cfac)
#endif
      do jb = 1, nlat
#ifndef MPI
         theta = th0 + (nlat-jb)*dth
#else
         theta = th0 + (nlat+joff-jb)*dth
#endif
         sint = sin(theta)
         if(dx(1,nlat-jb+1) > zap) then
            cfac = dtds*sqrt(dx(1,nlat-jb+1))/sint
         else
            cfac = zap
         end if
         ub(jb) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do jb = 1, nlat
            ub(jb) = 1.0d0
         end do
      end if
#endif
#ifndef USE_GPU
!$omp end parallel
#endif

      return
   end subroutine boundary_rwg
#else
#ifndef MPI
   subroutine boundary_rwg(dfld,hbnd,ubnd,dt,dxdy,nlon,nlat)
#else
   subroutine boundary_rwg(dfld,hbnd,ubnd,dt,dxdy,nlon,nlat,bflag)
#endif
      type(depth_arrays), target, intent(in) :: dfld
      type(boundary_arrays), target, intent(inout) :: hbnd, ubnd
      real(kind=REAL_BYTE), intent(in) :: dt, dxdy
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: bflag
#endif

      real(kind=REAL_BYTE) :: g
      real(kind=REAL_BYTE), pointer, dimension(:) :: ub
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy
      integer(kind=4) :: ib, jb
      real(kind=REAL_BYTE) :: dtds, cfac
      real(kind=REAL_BYTE) :: zap = 0.0d0

      hbnd%north = 0.0d0
      hbnd%east  = 0.0d0
      hbnd%south = 0.0d0
      hbnd%west  = 0.0d0

      dx => dfld%dx
      dy => dfld%dy

      g = 9.8d0  !*** m/s^2 ***

      dtds = dt*sqrt(g)/dxdy

      ! north
      ub => ubnd%north
#ifndef USE_GPU
!$omp parallel
#endif
#ifdef MPI
      if(iand(bflag, NORTH_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(cfac)
#else
!$omp target teams distribute parallel do private(cfac)
#endif
      do ib = 1, nlon
         if(dy(ib,1) > zap) then
            cfac = dtds*sqrt(dy(ib,1))
         else
            cfac = zap
         end if
         ub(ib) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do ib = 1, nlon
            ub(ib) = 1.0d0
         end do
      end if
#endif

      ! east
#ifndef USE_GPU
!$omp single
#endif
      ub => ubnd%east
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, EAST_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(cfac)
#else
!$omp target teams distribute parallel do private(cfac)
#endif
      do jb = 1, nlat
         if(dx(nlon-1,jb) > zap) then
            cfac = dtds*sqrt(dx(nlon-1,jb))
         else
            cfac = zap
         end if
         ub(jb) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do jb = 1, nlat
            ub(jb) = 1.0d0
         end do
      end if
#endif

      ! south; be careful, here we go backwards across bottom of grid
#ifndef USE_GPU
!$omp single
#endif
      ub => ubnd%south
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, SOUTH_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(cfac)
#else
!$omp target teams distribute parallel do private(cfac)
#endif
      do ib = 1, nlon
         if(dy(nlon-ib+1,nlat-1) > zap) then
            cfac = dtds*sqrt(dy(nlon-ib+1,nlat-1))
         else
            cfac = zap
         end if
         ub(ib) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do ib = 1, nlon
            ub(ib) = 1.0d0
         end do
      end if
#endif

      ! west; be careful, here we go backwards up the left of grid
#ifndef USE_GPU
!$omp single
#endif
      ub => ubnd%west
#ifndef USE_GPU
!$omp end single
#endif
#ifdef MPI
      if(iand(bflag, WEST_BOUND) /= 0) then
#endif
#ifndef USE_GPU
!$omp do private(cfac)
#else
!$omp target teams distribute parallel do private(cfac)
#endif
      do jb = 1, nlat
         if(dx(1,nlat-jb+1) > zap) then
            cfac = dtds*sqrt(dx(1,nlat-jb+1))
         else
            cfac = zap
         end if
         ub(jb) = (1.0d0 - cfac)/(1.0d0 + cfac)
      end do
#ifdef MPI
      else
#ifndef USE_GPU
!$omp do
#else
!$omp target teams distribute parallel do
#endif
         do jb = 1, nlat
            ub(jb) = 1.0d0
         end do
      end if
#endif
#ifndef USE_GPU
!$omp end parallel
#endif

      return
   end subroutine boundary_rwg
#endif

#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELIN
   subroutine wet_or_dry(wfld,dfld,ifz,nlon,nlat,fname,wodfld,formatid)
#else
   subroutine wet_or_dry(wfld,dfld,ifz,nlon,nlat,fname,wodfld,nxorg,nyorg,formatid)
#endif
#else
#ifndef PIXELIN
   subroutine wet_or_dry(wfld,dfld,ifz,nlon,nlat,fname,wodfld,dg,myrank,formatid)
#else
   subroutine wet_or_dry(wfld,dfld,ifz,nlon,nlat,fname,wodfld,dg,myrank,nxorg,nyorg,formatid)
#endif
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
! === Conversion from flux to velocity should be done right after calc. ========
!     integer(kind=4), dimension(nlon,nlat), intent(out) :: ifz
#ifndef MPI
      integer(kind=4), dimension(nlon,nlat), intent(out) :: ifz
#else
      integer(kind=4), dimension(0:nlon+1,0:nlat+1), intent(out) :: ifz
#endif
! ==============================================================================
      integer(kind=4), intent(in) :: nlon, nlat
      character(len=256), intent(in) :: fname
      real(kind=REAL_BYTE), target, dimension(nlon,nlat), intent(inout) :: wodfld
#ifdef PIXELIN
      integer(kind=4), intent(in) :: nxorg, nyorg
#endif

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz, wod
! === Arrival time =============================================================
      integer(kind=4), pointer, dimension(:,:) :: arrivedat
! ==============================================================================
      integer(kind=4), pointer, dimension(:,:) :: tttdat
      real(kind=REAL_BYTE) :: zap = 0.0d0
      integer(kind=4) :: i, j
#if defined(MPI) && defined(ONEFILE)
      type(data_grids), target, intent(inout) :: dg
      integer(kind=4), intent(in) :: myrank
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: wod_all
#endif
#ifdef PIXELIN
      real(kind=REAL_BYTE), allocatable, dimension(:,:) :: wodorg
#endif
      integer(kind=4), intent(in) :: formatid

      hz => wfld%hz
      dz => dfld%dz
      wod => wodfld
! === Arrival time =============================================================
      if(check_arrival_time == 1) arrivedat => wfld%arrivedat
! ==============================================================================
      if(check_tt_time == 1) tttdat => wfld%tttdat

#ifndef __NEC__
#ifndef USE_GPU
!$omp parallel
#endif
#endif
      if(trim(fname(1:11)) == 'NO_WETORDRY') then
#ifndef __NEC__
#ifndef USE_GPU
!$omp single
#endif
#endif
#if defined(MPI) && defined(ONEFILE)
         if(myrank == 0) then
#endif
         write(6,'(8x,a)') trim(fname)
#if defined(MPI) && defined(ONEFILE)
         end if
#endif
#ifndef USE_GPU
#ifndef __NEC__
!$omp end single
!$omp do private(i)
#else
!$omp parallel do private(i)
#endif
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon
               if(dz(i,j) > zap) then ! wet
                  ifz(i,j) = 1
                  hz(i,j) = zap
               else ! dry
                  ifz(i,j) = -1
                  hz(i,j) = zap
               end if
            end do
         end do
      else
#ifndef __NEC__
#ifndef USE_GPU
!$omp single
#endif
#endif
#if !defined(MPI) || !defined(ONEFILE)
         write(6,'(8x,a,a)') 'WETORDRY_FILE_GIVEN:', trim(fname)
#ifndef PIXELIN
         call read_gmt_grd(fname, wod, nlon, nlat,formatid)
#else
         allocate(wodorg(0:nxorg-1,0:nyorg-1))
         open(1,file=trim(fname),action='read',status='old',form='formatted')
         read(1,'(10f8.2)') wodorg
         close(1)
         wod(1:nlon,1:nlat) = wodorg(1:nlon,1:nlat)
         deallocate(wodorg)
#endif
#else
         if(myrank == 0) then
            allocate(wod_all(dg%my%totalNx,dg%my%totalNy))
            write(6,'(8x,a,a)') 'WETORDRY_FILE_GIVEN:', trim(fname)
#ifndef PIXELIN
            call read_gmt_grd(fname, wod_all, dg%my%totalNx,dg%my%totalNy,formatid)
#else
            allocate(wodorg(0:nxorg-1,0:nyorg-1))
            open(1,file=trim(fname),action='read',status='old',form='formatted')
            read(1,'(10f8.2)') wodorg
            close(1)
            wod_all(1:dg%my%totalNx,1:dg%my%totalNy) = wodorg(1:dg%my%totalNx,1:dg%my%totalNy)
         deallocate(wodorg)
#endif
         else
            allocate(wod_all(1,1))
         end if
         call onefile_scatter_array(wod_all,wod,dg)
         deallocate(wod_all)
#endif
#ifndef USE_GPU
#ifndef __NEC__
!$omp end single
!$omp do private(i)
#else
!$omp parallel do private(i)
#endif
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon
               if(dz(i,j) > zap .and. wod(i,j) > zap) then ! wet
                  ifz(i,j) = 1
                  hz(i,j) = zap
               else ! dry
                  ifz(i,j) = -1
                  ! hz should be initally set to be elevation for anticipated inundation area
! === Flood Change =============================================================
!                 if(wod(i,j) > -30.0d0 .and. abs(dz(i,j)) > 0.01d0) then
!                    hz(i,j) = -dz(i,j)-0.01d0
                  if(wod(i,j) < 0.0d0 .and. dz(i,j) >= 0.0d0) then
                     hz(i,j) = -dz(i,j)-1.0d-6
! ==============================================================================
                  else
                     hz(i,j) = zap
                  end if
               end if
            end do
         end do
      end if
! === Arrival time =============================================================
      if(check_arrival_time == 1) then
#ifndef USE_GPU
#ifndef __NEC__
!$omp do private(i)
#else
!$omp parallel do private(i)
#endif
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon
               if(ifz(i,j) == 1) then ! wet first
                  arrivedat(i,j) = -1
               else ! dry first
                  arrivedat(i,j) = -2
               end if
            end do
         end do
      end if
! ==============================================================================
      if(check_tt_time == 1) then
#ifndef USE_GPU
#ifndef __NEC__
!$omp do private(i)
#else
!$omp parallel do private(i)
#endif
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon
               if(ifz(i,j) == 1) then ! wet first
                  tttdat(i,j) = -1
               else ! dry first
                  tttdat(i,j) = -2
               end if
            end do
         end do
      end if
#ifndef __NEC__
#ifndef USE_GPU
!$omp end parallel
#endif
#endif

      return
   end subroutine wet_or_dry

#ifndef MPI
! === For negative max. height =================================================
!  subroutine minmax_rwg(nx,ny,a,zmax,zmin,imin,jmin,imax,jmax)
   subroutine minmax_rwg(nx,ny,a,zmax,zmin,imin,jmin,imax,jmax, &
                         flag_missing_value)
! ==============================================================================
#else
! === For negative max. height =================================================
!  subroutine minmax_rwg(nx,ny,a,zmax,zmin)
   subroutine minmax_rwg(nx,ny,a,zmax,zmin,flag_missing_value)
! ==============================================================================
      use mpi
#endif
#ifdef _OPENMP
      use omp_lib
#endif
      integer(kind=4), intent(in) :: nx, ny
      real(kind=REAL_BYTE), dimension(nx,ny), intent(in) :: a
      real(kind=REAL_BYTE), intent(out) :: zmax, zmin
#ifndef MPI
      integer(kind=4), intent(out) :: imin, jmin, imax, jmax
#endif
! === For negative max. height =================================================
      logical, optional, intent(in) :: flag_missing_value
! ==============================================================================

      real(kind=REAL_BYTE) :: vmin, vmax
      integer(kind=4) :: i, j
#ifdef MPI
      integer(kind=4) :: ierr
#endif
#ifdef _OPENMP
#ifndef __NEC__
      integer(kind=4), parameter :: mt = 256
#else
      integer(kind=4), parameter :: mt = 8
#endif
#ifndef MPI
      integer(kind=4), dimension(0:mt-1) :: imin_, jmin_, imax_, jmax_
#endif
      real(kind=REAL_BYTE), dimension(0:mt-1) :: min_, max_
      integer(kind=4) :: t
#endif
      logical :: missing_value_is_available
   
#ifndef MPI
      imin = 1
      jmin = 1
      imax = 1
      jmax = 1
#endif

#ifndef _OPENMP
! === For negative max. height =================================================
      missing_value_is_available = .false.
      if(present(flag_missing_value)) then
         if(flag_missing_value) then
            missing_value_is_available = .true.
         end if
      end if

      if(missing_value_is_available) then
! ==============================================================================
         vmin = -missing_value
         vmax =  missing_value

#ifndef USE_GPU
         do j = 1, ny
            do i = 1, nx
               if(a(i,j) /= missing_value) then
                  if(a(i,j) < vmin) then
                     vmin = a(i,j)
#ifndef MPI
                     imin = i
                     jmin = j
#endif
                  end if
                  if(a(i,j) > vmax) then
                     vmax = a(i,j)
#ifndef MPI
                     imax = i
                     jmax = j
#endif
                  end if
               end if
            end do
         end do
#else
!$omp target teams distribute parallel do collapse(2) private(i) &
!$omp reduciton(min:vmin,max:vmax)
         do j = 1, ny
            do i = 1, nx
               if(a(i,j) /= missing_value) then
                  vmin = min(vmin, a(i,j))
                  vmax = max(vmax, a(i,j))
               end if
            end do
         end do
#endif
! === For negative max. height =================================================
      else
! ==============================================================================
      vmin = a(1,1)
      vmax = a(1,1)

#ifndef USE_GPU
      do j = 1, ny
         do i = 1, nx
            if(a(i,j) < vmin) then
               vmin = a(i,j)
#ifndef MPI
               imin = i
               jmin = j
#endif
            end if
            if(a(i,j) > vmax) then
               vmax = a(i,j)
#ifndef MPI
               imax = i
               jmax = j
#endif
            end if
         end do
      end do
#else
!$omp target teams distribute parallel do collapse(2) private(i) &
!$omp reduciton(min:vmin,max:vmax)
         do j = 1, ny
            do i = 1, nx
               vmin = min(vmin, a(i,j))
               vmax = max(vmax, a(i,j))
            end do
         end do
#endif
! === For negative max. height =================================================
      end if
! ==============================================================================
#ifndef MPI
!$omp target teams distribute parallel do collapse(2) private(i)
      do j = 1, ny
         do i = 1, nx
            if(a(i,j) == vmin) then
!$omp atomic write
               imin = i
!$omp atomic write
               jmin = j
            end if
            if(a(i,j) == vmax) then
!$omp atomic write
               imax = i
!$omp atomic write
               jmax = j
            end if
         end do
      end do
#endif
#else
      missing_value_is_available = .false.
      if(present(flag_missing_value)) then
         if(flag_missing_value) then
            missing_value_is_available = .true.
         end if
      end if

      if(missing_value_is_available) then
         vmin = -missing_value
         vmax =  missing_value
         min_ = vmin
         max_ = vmax

!$omp parallel private(t)
         t = omp_get_thread_num()
!$omp do private(i)
         do j = 1, ny
            do i = 1, nx
               if(a(i,j) /= missing_value) then
                  if(a(i,j) < min_(t)) then
                     min_(t) = a(i,j)
#ifndef MPI
                     imin_(t) = i
                     jmin_(t) = j
#endif
                  end if
                  if(a(i,j) > max_(t)) then
                     max_(t) = a(i,j)
#ifndef MPI
                     imax_(t) = i
                     jmax_(t) = j
#endif
                  end if
               end if
            end do
         end do
!$omp end parallel
         do t = 0, mt-1
            if(min_(t) < vmin) then
               vmin = min_(t)
#ifndef MPI
               imin = imin_(t)
               jmin = jmin_(t)
#endif
            end if
            if(max_(t) > vmax) then
               vmax = max_(t)
#ifndef MPI
               imax = imax_(t)
               jmax = jmax_(t)
#endif
            end if
         end do
      else
         vmin = a(1,1)
         vmax = a(1,1)
         min_ = vmin
         max_ = vmax

!$omp parallel private(t)
         t = omp_get_thread_num()
!$omp do private(i)
         do j = 1, ny
            do i = 1, nx
               if(a(i,j) < min_(t)) then
                  min_(t) = a(i,j)
#ifndef MPI
                  imin_(t) = i
                  jmin_(t) = j
#endif
               end if
               if(a(i,j) > max_(t)) then
                  max_(t) = a(i,j)
#ifndef MPI
                  imax_(t) = i
                  jmax_(t) = j
#endif
               end if
            end do
         end do
!$omp end parallel
         do t = 0, mt-1
            if(min_(t) < vmin) then
               vmin = min_(t)
#ifndef MPI
               imin = imin_(t)
               jmin = jmin_(t)
#endif
            end if
            if(max_(t) > vmax) then
               vmax = max_(t)
#ifndef MPI
               imax = imax_(t)
               jmax = jmax_(t)
#endif
            end if
         end do
      end if
#endif

#ifndef MPI
      zmin = vmin
      zmax = vmax
#else
      call MPI_Allreduce(vmin, zmin, 1, REAL_MPI, MPI_MIN, __MPICOMM__, ierr)
      if(ierr /= 0) write(0,'(a)') 'MPI Error : MPI_Allreduce in minmax_rwg'
      call MPI_Allreduce(vmax, zmax, 1, REAL_MPI, MPI_MAX, __MPICOMM__, ierr)
      if(ierr /= 0) write(0,'(a)') 'MPI Error : MPI_Allreduce in minmax_rwg'
#endif

      return
   end subroutine minmax_rwg

   subroutine smooth_bath(cg,fg)
#ifdef MPI
      use mpi
#endif
      type(data_grids), target, intent(inout) :: cg, fg
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: dx, dy, dz, cdx, cdy, cdz
#ifndef MPI
      integer(kind=4) :: bigIX, bigIY, ix, iy
! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!     integer(kind=4) :: k1, k2, k3
! ==============================================================================
#else
      integer(kind=4) :: ix, iy
#endif
      real(kind=REAL_BYTE) :: lfac, x0, x1, y0, y1
      type(grid_info), pointer :: parent, my
      integer(kind=4) :: i0, j0, i1, j1, i2, j2, i3, j3
#ifdef MPI
      real(kind=REAL_BYTE), allocatable, dimension(:,:,:) :: cbuf0, cbuf1, cbuf2, cbuf3
      real(kind=REAL_BYTE), allocatable, dimension(:,:,:) :: fbuf0, fbuf1, fbuf2, fbuf3
      integer(kind=4) :: i, j, ldx, ldy, ixcbgn, ixcend, iycbgn, iycend, ixb, iyb, ierr
#endif

      parent => fg%parent
      my => fg%my

      cdx => cg%depth_field%dx
      cdy => cg%depth_field%dy
      cdz => cg%depth_field%dz

      dx => fg%depth_field%dx
      dy => fg%depth_field%dy
      dz => fg%depth_field%dz

#ifndef MPI
      bigIY = my%zeroIY
      do iy = 1, my%nr+1, my%nr
         bigIX = my%zeroIX
         do ix = 1, my%nx, my%nr
            dx(ix,iy) = cdx(bigIX,bigIY)
            dy(ix,iy) = cdy(bigIX,bigIY)
            dz(ix,iy) = cdz(bigIX,bigIY)
            bigIX = bigIX+1
         end do
         bigIY = bigIY+1
      end do

      bigIY = my%zeroIY + 2
      do iy = 2*my%nr+1, my%ny-2*my%nr, my%nr
         bigIX = my%zeroIX
         do ix = 1, my%nr+1, my%nr
            dx(ix,iy) = cdx(bigIX,bigIY)
            dy(ix,iy) = cdy(bigIX,bigIY)
            dz(ix,iy) = cdz(bigIX,bigIY)
            bigIX = bigIX+1
         end do
         bigIX = my%zeroIX + my%bigNX - 2
         do ix = my%nx-my%nr, my%nx, my%nr
            dx(ix,iy) = cdx(bigIX,bigIY)
            dy(ix,iy) = cdy(bigIX,bigIY)
            dz(ix,iy) = cdz(bigIX,bigIY)
            bigIX = bigIX+1
         end do
         bigIY = bigIY+1
      end do

      bigIY = my%zeroIY + my%bigNY - 2
      do iy = my%ny-my%nr, my%ny, my%nr
         bigIX = my%zeroIX
         do ix = 1, my%nx, my%nr
            dx(ix,iy) = cdx(bigIX,bigIY)
            dy(ix,iy) = cdy(bigIX,bigIY)
            dz(ix,iy) = cdz(bigIX,bigIY)
            bigIX = bigIX+1
         end do
         bigIY = bigIY+1
      end do

      lfac = 1.0d0/REAL_FUNC(my%nr)

      do iy = 1, 2*my%nr+1
         do ix = 1, my%nx
            i0 = ((ix-1)/my%nr)*my%nr + 1
            j0 = ((iy-1)/my%nr)*my%nr + 1
! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           i1 = mod(i0 + my%nr - 1, my%nx) + 1
!           j1 = j0 + (i0 + my%nr - 1)/my%nx
            i1 = i0 + my%nr
            j1 = j0
! ==============================================================================
            i2 = i0
            j2 = j0 + my%nr
            i3 = i1
            j3 = j1 + my%nr

! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           k1 = i1 + (j1-1)*my%nx
!           k2 = i2 + (j2-1)*my%nx
!           k3 = i3 + (j3-1)*my%nx
!           if(k1 < 1 .or. k1 > my%nx*my%ny) then
!              i1 = i0
!              j1 = j0
!           end if
!           if(k2 < 1 .or. k2 > my%nx*my%ny) then
!              i2 = i0
!              j2 = j0
!           end if
!           if(k3 < 1 .or. k3 > my%nx*my%ny) then
!              i3 = i0
!              j3 = j0
!           end if
            if(ix + my%nr > my%nx) then
               i1 = i1 - my%nr
               i3 = i3 - my%nr
            end if
            if(iy + my%nr > my%ny) then
               j2 = j2 - my%nr
               j3 = j3 - my%nr
            end if
! ==============================================================================

            x1 = mod(ix-1,my%nr)*lfac
            x0 = 1.0d0 - x1
            y1 = mod(iy-1,my%nr)*lfac
            y0 = 1.0d0 - y1

            dx(ix,iy) = x0*y0*dx(i0,j0) + x1*y0*dx(i1,j1) + x0*y1*dx(i2,j2) + x1*y1*dx(i3,j3)
            dy(ix,iy) = x0*y0*dy(i0,j0) + x1*y0*dy(i1,j1) + x0*y1*dy(i2,j2) + x1*y1*dy(i3,j3)
            dz(ix,iy) = x0*y0*dz(i0,j0) + x1*y0*dz(i1,j1) + x0*y1*dz(i2,j2) + x1*y1*dz(i3,j3)
         end do
      end do

      do iy = 2*my%nr+1, my%ny-2*my%nr
         do ix = 1, 2*my%nr+1
            i0 = ((ix-1)/my%nr)*my%nr + 1
            j0 = ((iy-1)/my%nr)*my%nr + 1
! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           i1 = mod(i0 + my%nr - 1, my%nx) + 1
!           j1 = j0 + (i0 + my%nr - 1)/my%nx
            i1 = i0 + my%nr
            j1 = j0
! ==============================================================================
            i2 = i0
            j2 = j0 + my%nr
            i3 = i1
            j3 = j1 + my%nr

! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           k1 = i1 + (j1-1)*my%nx
!           k2 = i2 + (j2-1)*my%nx
!           k3 = i3 + (j3-1)*my%nx
!           if(k1 < 1 .or. k1 > my%nx*my%ny) then
!              i1 = i0
!              j1 = j0
!           end if
!           if(k2 < 1 .or. k2 > my%nx*my%ny) then
!              i2 = i0
!              j2 = j0
!           end if
!           if(k3 < 1 .or. k3 > my%nx*my%ny) then
!              i3 = i0
!              j3 = j0
!           end if
            if(ix + my%nr > my%nx) then
               i1 = i1 - my%nr
               i3 = i3 - my%nr
            end if
            if(iy + my%nr > my%ny) then
               j2 = j2 - my%nr
               j3 = j3 - my%nr
            end if
! ==============================================================================

            x1 = mod(ix-1,my%nr)*lfac
            x0 = 1.0d0 - x1
            y1 = mod(iy-1,my%nr)*lfac
            y0 = 1.0d0 - y1
   
            dx(ix,iy) = x0*y0*dx(i0,j0) + x1*y0*dx(i1,j1) + x0*y1*dx(i2,j2) + x1*y1*dx(i3,j3)
            dy(ix,iy) = x0*y0*dy(i0,j0) + x1*y0*dy(i1,j1) + x0*y1*dy(i2,j2) + x1*y1*dy(i3,j3)
            dz(ix,iy) = x0*y0*dz(i0,j0) + x1*y0*dz(i1,j1) + x0*y1*dz(i2,j2) + x1*y1*dz(i3,j3)
         end do
         do ix = my%nx-2*my%nr, my%nx
            i0 = ((ix-1)/my%nr)*my%nr + 1
            j0 = ((iy-1)/my%nr)*my%nr + 1
! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           i1 = mod(i0 + my%nr - 1, my%nx) + 1
!           j1 = j0 + (i0 + my%nr - 1)/my%nx
            i1 = i0 + my%nr
            j1 = j0
! ==============================================================================
            i2 = i0
            j2 = j0 + my%nr
            i3 = i1
            j3 = j1 + my%nr

! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           k1 = i1 + (j1-1)*my%nx
!           k2 = i2 + (j2-1)*my%nx
!           k3 = i3 + (j3-1)*my%nx
!           if(k1 < 1 .or. k1 > my%nx*my%ny) then
!              i1 = i0
!              j1 = j0
!           end if
!           if(k2 < 1 .or. k2 > my%nx*my%ny) then
!              i2 = i0
!              j2 = j0
!           end if
!           if(k3 < 1 .or. k3 > my%nx*my%ny) then
!              i3 = i0
!              j3 = j0
!           end if
            if(ix + my%nr > my%nx) then
               i1 = i1 - my%nr
               i3 = i3 - my%nr
            end if
            if(iy + my%nr > my%ny) then
               j2 = j2 - my%nr
               j3 = j3 - my%nr
            end if
! ==============================================================================

            x1 = mod(ix-1,my%nr)*lfac
            x0 = 1.0d0 - x1
            y1 = mod(iy-1,my%nr)*lfac
            y0 = 1.0d0 - y1

            dx(ix,iy) = x0*y0*dx(i0,j0) + x1*y0*dx(i1,j1) + x0*y1*dx(i2,j2) + x1*y1*dx(i3,j3)
            dy(ix,iy) = x0*y0*dy(i0,j0) + x1*y0*dy(i1,j1) + x0*y1*dy(i2,j2) + x1*y1*dy(i3,j3)
            dz(ix,iy) = x0*y0*dz(i0,j0) + x1*y0*dz(i1,j1) + x0*y1*dz(i2,j2) + x1*y1*dz(i3,j3)
         end do
      end do

      do iy = my%ny-2*my%nr, my%ny
         do ix = 1, my%nx
            i0 = ((ix-1)/my%nr)*my%nr + 1
            j0 = ((iy-1)/my%nr)*my%nr + 1
! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           i1 = mod(i0 + my%nr - 1, my%nx) + 1
!           j1 = j0 + (i0 + my%nr - 1)/my%nx
            i1 = i0 + my%nr
            j1 = j0
! ==============================================================================
            i2 = i0
            j2 = j0 + my%nr
            i3 = i1
            j3 = j1 + my%nr

! === DEBUG for MPI result check!!! by tkato 2012/08/13 ========================
!           k1 = i1 + (j1-1)*my%nx
!           k2 = i2 + (j2-1)*my%nx
!           k3 = i3 + (j3-1)*my%nx
!           if(k1 < 1 .or. k1 > my%nx*my%ny) then
!              i1 = i0
!              j1 = j0
!           end if
!           if(k2 < 1 .or. k2 > my%nx*my%ny) then
!              i2 = i0
!              j2 = j0
!           end if
!           if(k3 < 1 .or. k3 > my%nx*my%ny) then
!              i3 = i0
!              j3 = j0
!           end if
            if(ix + my%nr > my%nx) then
               i1 = i1 - my%nr
               i3 = i3 - my%nr
            end if
            if(iy + my%nr > my%ny) then
               j2 = j2 - my%nr
               j3 = j3 - my%nr
            end if
! ==============================================================================

            x1 = mod(ix-1,my%nr)*lfac
            x0 = 1.0d0 - x1
            y1 = mod(iy-1,my%nr)*lfac
            y0 = 1.0d0 - y1
   
            dx(ix,iy) = x0*y0*dx(i0,j0) + x1*y0*dx(i1,j1) + x0*y1*dx(i2,j2) + x1*y1*dx(i3,j3)
            dy(ix,iy) = x0*y0*dy(i0,j0) + x1*y0*dy(i1,j1) + x0*y1*dy(i2,j2) + x1*y1*dy(i3,j3)
            dz(ix,iy) = x0*y0*dz(i0,j0) + x1*y0*dz(i1,j1) + x0*y1*dz(i2,j2) + x1*y1*dz(i3,j3)
         end do
      end do
#else
      lfac = 1.0d0/REAL_FUNC(my%nr)

      ldx = my%totalNx/my%nr + 1
      ldy = my%totalNy/my%nr + 1

      ixcbgn = max(my%zeroIX+my%ix   /my%nr-1, parent%ix   )
      ixcend = min(my%zeroIX+my%ixend/my%nr+1, parent%ixend)
      iycbgn = max(my%zeroIY+my%iy   /my%nr-1, parent%iy   )
      iycend = min(my%zeroIY+my%iyend/my%nr+1, parent%iyend)

      !*==============*
      !*  allocate    *
      !*==============*
      allocate(cbuf0(ldx,4,3)); cbuf0 = 0.0d0
      allocate(cbuf1(ldy,4,3)); cbuf1 = 0.0d0
      allocate(cbuf2(ldy,4,3)); cbuf2 = 0.0d0
      allocate(cbuf3(ldx,4,3)); cbuf3 = 0.0d0

      allocate(fbuf0(ldx,4,3)); fbuf0 = 0.0d0
      allocate(fbuf1(ldy,4,3)); fbuf1 = 0.0d0
      allocate(fbuf2(ldy,4,3)); fbuf2 = 0.0d0
      allocate(fbuf3(ldx,4,3)); fbuf3 = 0.0d0

      !*=================================================*
      !*                 Update dx                       *
      !*=================================================*
      !*==============*
      !*  copy2buf    *
      !*==============*
      ! (1)
      do iyb = 1, 2 ! buffer
         iy = iyb + my%zeroIY - 1 ! global
         if(parent%iy <= iy .and. iy <= parent%iyend) then
            j = iy - parent%ky + 1 ! local
            do ixb = 1, my%bigNX ! buffer
                ix = ixb + my%zeroIX - 1 ! global
                if(parent%ix <= ix .and. ix <= parent%ixend) then
                   i = ix - parent%kx + 1 ! local
                   cbuf0(ixb,iyb,1) = cdx(i,j)
                   cbuf0(ixb,iyb,2) = cdy(i,j)
                   cbuf0(ixb,iyb,3) = cdz(i,j)
                end if
            end do
         end if
      end do
      ! (2)
      do iyb = 3, my%bigNY-1 ! buffer
         iy = iyb + my%zeroIY - 1 ! global
         if(parent%iy <= iy .and. iy <= parent%iyend) then
            j = iy - parent%ky + 1 ! local
            do ixb = 1, 2 ! buffer
               ix = ixb + my%zeroIX - 1 ! global
               if(parent%ix <= ix .and. ix <= parent%ixend) then
                  i = ix - parent%kx + 1 ! local
                  cbuf1(iyb,ixb,1) = cdx(i,j)
                  cbuf1(iyb,ixb,2) = cdy(i,j)
                  cbuf1(iyb,ixb,3) = cdz(i,j)
               end if
            end do
         end if
      end do
      ! (3)
      do iyb = 3, my%bigNY-1 ! buffer
         iy = iyb + my%zeroIY - 1 ! global
         if(parent%iy <= iy .and. iy <= parent%iyend) then
            j = iy - parent%ky + 1 ! local
            do ixb = 1, 2 ! buffer
               ix = ixb + (my%zeroIX + my%bigNX - 1) - 2 ! global
               if(parent%ix <= ix .and. ix <= parent%ixend) then
                  i = ix - parent%kx + 1 ! local
                  cbuf2(iyb,ixb,1) = cdx(i,j)
                  cbuf2(iyb,ixb,2) = cdy(i,j)
                  cbuf2(iyb,ixb,3) = cdz(i,j)
               end if
            end do
         end if
      end do
      ! (4)
      do iyb = 1, 2 ! buffer
         iy = iyb + (my%zeroIY + my%bigNY - 1) -2 ! global
         if(parent%iy <= iy .and. iy <= parent%iyend) then
            j = iy - parent%ky + 1 ! local
            do ixb = 1, my%bigNX ! buffer
               ix = ixb + my%zeroIX - 1 ! global
               if(parent%ix <= ix .and. ix <= parent%ixend) then
                  i = ix - parent%kx + 1 ! local
                  cbuf3(ixb,iyb,1) = cdx(i,j)
                  cbuf3(ixb,iyb,2) = cdy(i,j)
                  cbuf3(ixb,iyb,3) = cdz(i,j)
               end if
            end do
         end if
      end do
      !*==============*
      !*  allreduce   *
      !*==============*
      call MPI_Allreduce(cbuf0, fbuf0, ldx*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      call MPI_Allreduce(cbuf1, fbuf1, ldy*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      call MPI_Allreduce(cbuf2, fbuf2, ldy*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      call MPI_Allreduce(cbuf3, fbuf3, ldx*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      !*==============*
      !*  copy2dx     *
      !*==============*
      ! (1)
      do iy = 1, my%nr+1, my%nr ! global
         do ix = 1, my%totalNx, my%nr ! global
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ixb = ix/my%nr + 1 ! buffer
               iyb = iy/my%nr + 1 ! buffer
               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = fbuf0(ixb,iyb,1)
               dy(i,j) = fbuf0(ixb,iyb,2)
               dz(i,j) = fbuf0(ixb,iyb,3)
            end if
         end do
      end do
      ! (2)
      do iy = 2*my%nr+1, my%totalNy-2*my%nr, my%nr ! global
         do ix = 1, my%nr+1, my%nr ! global
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ixb = ix/my%nr + 1 ! buffer
               iyb = iy/my%nr + 1 ! buffer
               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = fbuf1(iyb,ixb,1)
               dy(i,j) = fbuf1(iyb,ixb,2)
               dz(i,j) = fbuf1(iyb,ixb,3)
            end if
         end do
      end do
      ! (3)
      do iy = 2*my%nr+1, my%totalNy-2*my%nr, my%nr ! global
         do ix = my%totalNx - my%nr, my%totalNx, my%nr ! global
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ixb = (ix - (my%totalNx - my%nr))/my%nr + 1 ! buffer
               iyb = iy/my%nr + 1 ! buffer
               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = fbuf2(iyb,ixb,1)
               dy(i,j) = fbuf2(iyb,ixb,2)
               dz(i,j) = fbuf2(iyb,ixb,3)
            end if
         end do
      end do
      ! (4)
      do iy = my%totalNy - my%nr, my%totalNy, my%nr ! global
         do ix = 1, my%totalNx, my%nr ! global
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ixb = ix/my%nr + 1 ! buffer
               iyb = (iy - (my%totalNy - my%nr))/my%nr + 1 ! buffer
               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = fbuf3(ixb,iyb,1)
               dy(i,j) = fbuf3(ixb,iyb,2)
               dz(i,j) = fbuf3(ixb,iyb,3)
            end if
         end do
      end do
      !*==============*
      !*  clearbuf    *
      !*==============*
      cbuf0 = 0.0d0
      cbuf1 = 0.0d0
      cbuf2 = 0.0d0
      cbuf3 = 0.0d0

      fbuf0 = 0.0d0
      fbuf1 = 0.0d0
      fbuf2 = 0.0d0
      fbuf3 = 0.0d0
      !*==============*
      !*  fine2buf    *
      !*==============*
      ! (1)
      do iy = 1, 4*my%nr, my%nr
         if(my%iy <= iy .and. iy <= my%iyend) then
            j = iy - my%ky + 1 ! local
            do ix = 1, my%totalNx, my%nr
               if(my%ix <= ix .and. ix <= my%ixend) then
                  i = ix - my%kx + 1 ! local
                  ixb = ix/my%nr + 1 ! buffer
                  iyb = iy/my%nr + 1 ! buffer
                  cbuf0(ixb,iyb,1) = dx(i,j)
                  cbuf0(ixb,iyb,2) = dy(i,j)
                  cbuf0(ixb,iyb,3) = dz(i,j)
               end if
            end do
         end if
      end do
      ! (2)
      do iy = 1, my%totalNy-my%nr, my%nr
         if(my%iy <= iy .and. iy <= my%iyend) then
            j = iy - my%ky + 1 ! local
            do ix = 1, 4*my%nr, my%nr
               if(my%ix <= ix .and. ix <= my%ixend) then
                  i = ix - my%kx + 1 ! local
                  ixb = ix/my%nr + 1 ! buffer
                  iyb = iy/my%nr + 1 ! buffer
                  cbuf1(iyb,ixb,1) = dx(i,j)
                  cbuf1(iyb,ixb,2) = dy(i,j)
                  cbuf1(iyb,ixb,3) = dz(i,j)
                end if
             end do
         end if
      end do
      ! (3)
      do iy = 1, my%totalNy-my%nr, my%nr
         if(my%iy <= iy .and. iy <= my%iyend) then
            j = iy - my%ky + 1 ! local
            do ix = my%totalNx-2*my%nr, my%totalNx, my%nr
               if(my%ix <= ix .and. ix <= my%ixend) then
                  i = ix - my%kx + 1 ! local
                  ixb = (ix - (my%totalNx - 2*my%nr) + 1)/my%nr + 1 ! buffer
                  iyb = iy/my%nr + 1 ! buffer
                  cbuf2(iyb,ixb,1) = dx(i,j)
                  cbuf2(iyb,ixb,2) = dy(i,j)
                  cbuf2(iyb,ixb,3) = dz(i,j)
               end if
            end do
         end if
      end do
      ! (4)
      do iy = my%totalNy-2*my%nr, my%totalNy, my%nr
         if(my%iy <= iy .and. iy <= my%iyend) then
            j = iy - my%ky + 1 ! local
            do ix = 1, my%totalNx, my%nr
               if(my%ix <= ix .and. ix <= my%ixend) then
                  i = ix - my%kx + 1 ! local
                  ixb = ix/my%nr + 1 ! buffer
                  iyb = (iy - (my%totalNy - 2*my%nr) + 1)/my%nr + 1 ! buffer
                  cbuf3(ixb,iyb,1) = dx(i,j)
                  cbuf3(ixb,iyb,2) = dy(i,j)
                  cbuf3(ixb,iyb,3) = dz(i,j)
               end if
            end do
         end if
      end do
      !*==============*
      !*  allreduce   *
      !*==============*
      call MPI_Allreduce(cbuf0, fbuf0, ldx*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      call MPI_Allreduce(cbuf1, fbuf1, ldy*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      call MPI_Allreduce(cbuf2, fbuf2, ldy*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      call MPI_Allreduce(cbuf3, fbuf3, ldx*4*3, REAL_MPI, MPI_SUM, __MPICOMM__, ierr)
      !*==============*
      !*  buf2fine    *
      !*==============*
      ! (1)
      do iy = 1, 2*my%nr+1 ! global
         do ix = 1, my%totalNx ! global
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ! buffer coordinates
               i0 = (ix-1)/my%nr + 1
               j0 = (iy-1)/my%nr + 1
               i1 = i0 + 1
               j1 = j0
               i2 = i0
               j2 = j0 + 1
               i3 = i1
               j3 = j1 + 1
    
               if(ix + my%nr > my%totalNx) then
                  i1 = i1 - 1
                  i3 = i3 - 1
               end if
               if(iy + my%nr > my%totalNy) then
                  j2 = j2 - 1
                  j3 = j3 - 1
               end if

               x1 = mod(ix-1,my%nr)*lfac
               x0 = 1.0d0 - x1
               y1 = mod(iy-1,my%nr)*lfac
               y0 = 1.0d0 - y1

               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = x0*y0*fbuf0(i0,j0,1) + x1*y0*fbuf0(i1,j1,1) + x0*y1*fbuf0(i2,j2,1) + x1*y1*fbuf0(i3,j3,1)
               dy(i,j) = x0*y0*fbuf0(i0,j0,2) + x1*y0*fbuf0(i1,j1,2) + x0*y1*fbuf0(i2,j2,2) + x1*y1*fbuf0(i3,j3,2)
               dz(i,j) = x0*y0*fbuf0(i0,j0,3) + x1*y0*fbuf0(i1,j1,3) + x0*y1*fbuf0(i2,j2,3) + x1*y1*fbuf0(i3,j3,3)
            end if
         end do
      end do
      ! (2)
      do iy = 2*my%nr+1, my%totalNy-2*my%nr ! global
         do ix = 1, 2*my%nr+1 ! global
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ! buffer coordinates
               i0 = (ix-1)/my%nr + 1
               j0 = (iy-1)/my%nr + 1
               i1 = i0 + 1
               j1 = j0
               i2 = i0
               j2 = j0 + 1
               i3 = i1
               j3 = j1 + 1
    
               if(ix + my%nr > my%totalNx) then
                  i1 = i1 - 1
                  i3 = i3 - 1
               end if
               if(iy + my%nr > my%totalNy) then
                  j2 = j2 - 1
                  j3 = j3 - 1
               end if

               x1 = mod(ix-1,my%nr)*lfac
               x0 = 1.0d0 - x1
               y1 = mod(iy-1,my%nr)*lfac
               y0 = 1.0d0 - y1

               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = x0*y0*fbuf1(j0,i0,1) + x1*y0*fbuf1(j1,i1,1) + x0*y1*fbuf1(j2,i2,1) + x1*y1*fbuf1(j3,i3,1)
               dy(i,j) = x0*y0*fbuf1(j0,i0,2) + x1*y0*fbuf1(j1,i1,2) + x0*y1*fbuf1(j2,i2,2) + x1*y1*fbuf1(j3,i3,2)
               dz(i,j) = x0*y0*fbuf1(j0,i0,3) + x1*y0*fbuf1(j1,i1,3) + x0*y1*fbuf1(j2,i2,3) + x1*y1*fbuf1(j3,i3,3)
            end if
         end do
      end do
      ! (3)
      do iy = 2*my%nr+1, my%totalNy-2*my%nr ! global
         do ix = my%totalNx-2*my%nr, my%totalNx
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ! buffer coordinates
               i0 = (ix - (my%totalNx - 2*my%nr))/my%nr + 1
               j0 = (iy-1)/my%nr + 1
               i1 = i0 + 1
               j1 = j0
               i2 = i0
               j2 = j0 + 1
               i3 = i1
               j3 = j1 + 1
    
               if(ix + my%nr > my%totalNx) then
                  i1 = i1 - 1
                  i3 = i3 - 1
               end if
               if(iy + my%nr > my%totalNy) then
                  j2 = j2 - 1
                  j3 = j3 - 1
               end if

               x1 = mod(ix-1,my%nr)*lfac
               x0 = 1.0d0 - x1
               y1 = mod(iy-1,my%nr)*lfac
               y0 = 1.0d0 - y1

               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = x0*y0*fbuf2(j0,i0,1) + x1*y0*fbuf2(j1,i1,1) + x0*y1*fbuf2(j2,i2,1) + x1*y1*fbuf2(j3,i3,1)
               dy(i,j) = x0*y0*fbuf2(j0,i0,2) + x1*y0*fbuf2(j1,i1,2) + x0*y1*fbuf2(j2,i2,2) + x1*y1*fbuf2(j3,i3,2)
               dz(i,j) = x0*y0*fbuf2(j0,i0,3) + x1*y0*fbuf2(j1,i1,3) + x0*y1*fbuf2(j2,i2,3) + x1*y1*fbuf2(j3,i3,3)
            end if
         end do
      end do
      ! (4)
      do iy = my%totalNy-2*my%nr, my%totalNy
         do ix = 1, my%totalNx ! global
            if(my%kx <= ix .and. ix <= my%kxend .and. &
               my%ky <= iy .and. iy <= my%kyend) then
               ! buffer coordinates
               i0 = (ix-1)/my%nr + 1
               j0 = (iy - (my%totalNy - 2*my%nr))/my%nr + 1
               i1 = i0 + 1
               j1 = j0
               i2 = i0
               j2 = j0 + 1
               i3 = i1
               j3 = j1 + 1
    
               if(ix + my%nr > my%totalNx) then
                  i1 = i1 - 1
                  i3 = i3 - 1
               end if
               if(iy + my%nr > my%totalNy) then
                  j2 = j2 - 1
                  j3 = j3 - 1
               end if

               x1 = mod(ix-1,my%nr)*lfac
               x0 = 1.0d0 - x1
               y1 = mod(iy-1,my%nr)*lfac
               y0 = 1.0d0 - y1

               i = ix - my%kx + 1 ! local with edge
               j = iy - my%ky + 1 ! local with edge
               dx(i,j) = x0*y0*fbuf3(i0,j0,1) + x1*y0*fbuf3(i1,j1,1) + x0*y1*fbuf3(i2,j2,1) + x1*y1*fbuf3(i3,j3,1)
               dy(i,j) = x0*y0*fbuf3(i0,j0,2) + x1*y0*fbuf3(i1,j1,2) + x0*y1*fbuf3(i2,j2,2) + x1*y1*fbuf3(i3,j3,2)
               dz(i,j) = x0*y0*fbuf3(i0,j0,3) + x1*y0*fbuf3(i1,j1,3) + x0*y1*fbuf3(i2,j2,3) + x1*y1*fbuf3(i3,j3,3)
            end if
         end do
      end do
      !*==============*
      !*  deallocate  *
      !*==============*
      deallocate(cbuf0)
      deallocate(cbuf1)
      deallocate(cbuf2)
      deallocate(cbuf3)
      deallocate(fbuf0)
      deallocate(fbuf1)
      deallocate(fbuf2)
      deallocate(fbuf3)
#endif

      return
   end subroutine smooth_bath

! === Absorbing boundary condition =============================================
#ifndef MPI
   subroutine make_abc(abc,nlon,nlat)
#else
   subroutine make_abc(abc,nlon,nlat,totalNx,totalNy,kx,ky)
#endif
      real(kind=REAL_BYTE), dimension(nlon,nlat), intent(out) :: abc
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: totalNx, totalNy, kx, ky
      integer(kind=4) :: ig, jg
#endif
      real(kind=REAL_BYTE) :: tmp
      integer(kind=4) :: i, j, ii, jj

#ifndef USE_GPU
!$omp parallel
#endif
#ifndef MPI
#ifndef USE_GPU
!$omp do private(jj, ii, i, tmp)
#else
!$omp target teams distribute parallel do collapse(2) private(jj,ii,i,tmp)
#endif
      do j = 1, nlat
#ifndef USE_GPU
         jj = min(j, nlat+1-j)
#endif
         do i = 1, nlon
#ifdef USE_GPU
           jj = min(j, nlat+1-j)
#endif
            ii = min(i, nlon+1-i)
            if(ii <= nxa .and. jj > nya) then
               tmp = exp(-((apara*(nxa-ii))**2))
            else if(ii > nxa .and. jj <= nya) then
               tmp = exp(-((apara*(nya-jj))**2))
            else if(ii <= nxa .and. jj <= nya) then
               tmp = exp(-((apara*(nxa-ii))**2+(apara*(nya-jj))**2))
            else
               tmp = 1.0d0
            end if
            abc(i,j) = tmp
         end do
      end do
#else
#ifndef USE_GPU
!$omp do private(jj, ii, i, ig, jg, tmp)
#else
!$omp target teams distribute parallel do collapse(2) private(jj,ii,i,ig,jg,tmp)
#endif
      do j = 1, nlat
#ifndef USE_GPU
         jg = ky+j-1
         jj = min(jg, totalNy+1-jg)
#endif
         do i = 1, nlon
#ifdef USE_GPU
           jg = ky+j-1
           jj = min(jg, totalNy+1-jg)
#endif
            ig = kx+i-1
            ii = min(ig, totalNx+1-ig)
            if(ii <= nxa .and. jj > nya) then
               tmp = exp(-((apara*(nxa-ii))**2))
            else if(ii > nxa .and. jj <= nya) then
               tmp = exp(-((apara*(nya-jj))**2))
            else if(ii <= nxa .and. jj <= nya) then
               tmp = exp(-((apara*(nxa-ii))**2+(apara*(nya-jj))**2))
            else
               tmp = 1.0d0
            end if
            abc(i,j) = tmp
         end do
      end do
#endif
#ifndef USE_GPU
!$omp end parallel
#endif

      return
   end subroutine make_abc

   subroutine apply_abc(mode,wfld,nlon,nlat)
      type(wave_arrays), target, intent(inout) :: wfld
      integer(kind=4), intent(in) :: mode, nlon, nlat

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: fx, fy, hz, abc
      integer(kind=4) :: i, j

      abc => wfld%abc

      if(mode == VEL) then
         if(timenest /= 1) then
            fx => wfld%fx
            fy => wfld%fy
         else
            fx => wfld%fx_b
            fy => wfld%fy_b
         end if
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = 1, nlat
            do i = 1, nlon
               fx(i,j) = fx(i,j) * abc(i,j)
               fy(i,j) = fy(i,j) * abc(i,j)
            end do
         end do
      else
         if(timenest /= 1) then
            hz => wfld%hz
         else
            hz => wfld%hz_b
         end if
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
          do j = 1, nlat
             do i = 1, nlon
                hz(i,j) = hz(i,j) * abc(i,j)
             end do
          end do
      end if

      return
   end subroutine apply_abc
! ==============================================================================
! === Arrival time =============================================================
   subroutine check_arrival(wfld,dfld,nlon,nlat,istep)
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
      integer(kind=4), intent(in) :: nlon, nlat, istep

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz
      real(kind=REAL_BYTE) :: td
      integer(kind=4), pointer, dimension(:,:) :: arrivedat
      integer(kind=4) :: i, j

      hz        => wfld%hz
      dz        => dfld%dz
      arrivedat => wfld%arrivedat

#ifndef USE_GPU
!$omp parallel do private(i, td)
#else
!$omp target teams distribute parallel do collapse(2) private(i,td)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(arrivedat(i,j) == -1) then ! wet first
               if(hz(i,j) > check_arrival_height) arrivedat(i,j) = istep
            else if(arrivedat(i,j) == -2) then ! dry first
               td = dz(i,j) + hz(i,j)
               if(td > check_arrival_height) arrivedat(i,j) = istep
            end if
         end do
      end do

      return
   end subroutine check_arrival

   subroutine calc_arrival_time(wfld,nlon,nlat,dt)
      type(wave_arrays), target, intent(inout) :: wfld
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: dt

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: arrival_time
      integer(kind=4), pointer, dimension(:,:) :: arrivedat
      integer(kind=4) :: i, j

      arrival_time => wfld%arrival_time
      arrivedat    => wfld%arrivedat

#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(arrivedat(i,j) < 0) then
               arrival_time(i,j) = missing_value
            else
               arrival_time(i,j) = dt*arrivedat(i,j)
            end if
         end do
      end do

      return
   end subroutine calc_arrival_time
! ==============================================================================
   subroutine check_ttt(wfld,dfld,nlon,nlat,istep)
      type(wave_arrays), target, intent(inout) :: wfld
      type(depth_arrays), target, intent(in) :: dfld
      integer(kind=4), intent(in) :: nlon, nlat, istep

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz, dz
      real(kind=REAL_BYTE) :: td
      integer(kind=4), pointer, dimension(:,:) :: tttdat
      integer(kind=4) :: i, j

      hz     => wfld%hz
      dz     => dfld%dz
      tttdat => wfld%tttdat

#ifndef USE_GPU
!$omp parallel do private(i, td)
#else
!$omp target teams distribute parallel do collapse(2) private(i,td)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(tttdat(i,j) == -1) then ! wet first
               if(abs(hz(i,j)) > check_ttt_height) tttdat(i,j) = istep
            else if(tttdat(i,j) == -2) then ! dry first
               td = dz(i,j) + hz(i,j)
               if(td > check_ttt_height) tttdat(i,j) = istep
            end if
         end do
      end do

      return
   end subroutine check_ttt

   subroutine calc_tt_time(wfld,nlon,nlat,dt)
      type(wave_arrays), target, intent(inout) :: wfld
      integer(kind=4), intent(in) :: nlon, nlat
      real(kind=REAL_BYTE), intent(in) :: dt

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: tt_time
      integer(kind=4), pointer, dimension(:,:) :: tttdat
      integer(kind=4) :: i, j

      tt_time => wfld%tt_time
      tttdat  => wfld%tttdat

#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(tttdat(i,j) < 0) then
               tt_time(i,j) = missing_value
            else
               tt_time(i,j) = dt*tttdat(i,j)
            end if
         end do
      end do

      return
   end subroutine calc_tt_time

end module mod_rwg
