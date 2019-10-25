#include "real.h"
module mod_grid
#ifdef MPI
use mod_mpi_fatal
#endif
#ifdef MULTI
! === Split Dir ================================================================
use mod_multi, only : input_dirname
! ==============================================================================
#endif
implicit none
#ifdef MPI
integer(kind=4), parameter :: WEST_BOUND = 1, EAST_BOUND = 2, SOUTH_BOUND = 4, NORTH_BOUND = 8 ! kinds of boundary flags
#endif

#ifdef NCDIO
type ncdio_info
   integer(kind=4), dimension(3) :: start, count
   integer(kind=4) :: ncid, idid, mhid, mvid, hzid, vxid, vyid, speedid
! === Arrival time =============================================================
   integer(kind=4) :: atid
! ==============================================================================
#ifdef HZMINOUT
   integer(kind=4) :: minhid
#endif
   integer(kind=4) :: isize, ist, ien, jsize, jst, jen
   integer(kind=4) :: timeid, stepid
end type ncdio_info
#endif
type grid_info
   integer(kind=4) :: id = 1
   real(kind=REAL_BYTE) :: mlon0
   real(kind=REAL_BYTE) :: mlat0
   real(kind=REAL_BYTE) :: dh
   real(kind=REAL_BYTE) :: th0
   real(kind=REAL_BYTE) :: dth
   integer(kind=4) :: nx
   integer(kind=4) :: ny
#ifdef PIXELIN
   integer(kind=4) :: nxorg
   integer(kind=4) :: nyorg
#endif
   integer(kind=4) :: linear_flag
   character(len=256) :: base_name
   character(len=256) :: bath_file
   character(len=256) :: disp_file
   integer(kind=4) :: nr
#ifndef MPI
   integer(kind=4) :: bigNX
   integer(kind=4) :: bigNY
   integer(kind=4) :: zeroIX
   integer(kind=4) :: zeroIY
#else
   integer(kind=4) :: bigNX        ! change to global region value
   integer(kind=4) :: bigNY        ! change to global region value
   integer(kind=4) :: zeroIX       ! change to global region value
   integer(kind=4) :: zeroIY       ! change to global region value
   integer(kind=4) :: ix           ! x-position in global grid
   integer(kind=4) :: iy           ! y-position in global grid
   integer(kind=4) :: bigIX        ! x-position in global big grid
   integer(kind=4) :: bigIY        ! y-position in global big grid
   real(kind=REAL_BYTE) :: glon0           ! mlon0 of global region
   real(kind=REAL_BYTE) :: glat0           ! mlat0 of global region
   integer(kind=4) :: px           ! number of MPI proceses for x-direction
   integer(kind=4) :: py           ! number of MPI proceses for y-direction
   integer(kind=4) :: rx           ! MPI rank number for x-direction
   integer(kind=4) :: ry           ! MPI rank number for y-direction
   integer(kind=4) :: totalNx      ! nx of global region
   integer(kind=4) :: totalNy      ! ny of global region
   integer(kind=4) :: nbx          ! block size x
   integer(kind=4) :: nby          ! block size y
   integer(kind=4) :: ixend        ! end of local region for x-direction
   integer(kind=4) :: iyend        ! end of local region for y-direction
   integer(kind=4) :: has_boundary ! boundary flags
   integer(kind=4) :: kx           ! ix including edge
   integer(kind=4) :: ky           ! iy including edge
   integer(kind=4) :: kxend        ! ixend including edge
   integer(kind=4) :: kyend        ! iyend including edge
#endif
#ifdef NCDIO
   type(ncdio_info) :: ncdio
#endif
#ifdef CONV_CHECK
   integer(kind=4) :: nconvout
#endif
#if defined(MPI) && defined(ONEFILE)
   integer(kind=4), allocatable, dimension(:) :: kx_all, ky_all, kxend_all, kyend_all
   integer(kind=4), allocatable, dimension(:) :: ix_all, iy_all, ixend_all, iyend_all
   integer(kind=4) :: srcount_x, srcount_y, srcount
   real(kind=REAL_BYTE), allocatable, dimension(:) :: buf_g
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: buf_l
#ifdef NCDIO
   real(kind=4), allocatable, dimension(:) :: buf_g_ncdio
   real(kind=4), allocatable, dimension(:,:) :: buf_l_ncdio
#endif
#endif
end type grid_info

type interp_info
   integer(kind=4) :: np
   integer(kind=4), allocatable, dimension(:,:) :: fndx
   integer(kind=4), allocatable, dimension(:,:) :: cndx0
   integer(kind=4), allocatable, dimension(:,:) :: cndx1
   real(kind=REAL_BYTE), allocatable, dimension(:) :: wt0
   real(kind=REAL_BYTE), allocatable, dimension(:) :: wt1
#ifdef MPI
   integer(kind=4), pointer :: sendto
   integer(kind=4), pointer :: recvfm
   integer(kind=4) :: nsend
   integer(kind=4) :: nrecv
#ifdef USE_ALLTOALLV
! === USE_MPI_ALLTOALLV ========================================================
   integer(kind=4) :: snp0, rnp0, snp1, rnp1
   integer(kind=4), allocatable, dimension(:,:) :: fndx0_l, cndx0_l
   integer(kind=4), allocatable, dimension(:,:) :: fndx1_l, cndx1_l
   integer(kind=4), allocatable, dimension(:) :: sendcnts0, sdispls0, recvcnts0, rdispls0
   integer(kind=4), allocatable, dimension(:) :: sendcnts1, sdispls1, recvcnts1, rdispls1
   real(kind=REAL_BYTE), allocatable, dimension(:) :: wt0_l
   real(kind=REAL_BYTE), allocatable, dimension(:) :: wt1_l
! === USE_MPI_ALLTOALLV ========================================================
#ifdef SINGLE_A2A
! === SINGLE Alltoallv =========================================================
   integer(kind=4), allocatable, dimension(:) :: smap0, smap1
   integer(kind=4), allocatable, dimension(:) :: rmap0, rmap1
! ==============================================================================
#endif
#ifdef A2A3D
   integer(kind=4) :: handler0, handler1
#endif
#endif
! === Don't repeat allocate/deallocate! ========================================
   real(kind=REAL_BYTE), allocatable, dimension(:) :: fb, cb
#ifdef USE_ALLTOALLV
   real(kind=REAL_BYTE), allocatable, dimension(:) :: fb0, cb0
   real(kind=REAL_BYTE), allocatable, dimension(:) :: fb1, cb1
#endif
! ==============================================================================
#endif
end type interp_info

type wave_arrays
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fx
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fy
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hz
! === Flood Change =============================================================
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hz_old
! ==============================================================================
! === DEBUG for wave height gap on nest boundary. 2012/10/30 ===================
!  real(kind=REAL_BYTE), allocatable, dimension(:) :: fx_old
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fx_old
! ==============================================================================
! === DEBUG for wave height gap on nest boundary. 2012/10/30 ===================
!  real(kind=REAL_BYTE), allocatable, dimension(:) :: fy_old
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fy_old
! ==============================================================================
! === Dispersive ===============================================================
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: yu, yv ! RHS
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: cu, cv ! Coefficient
! ==============================================================================
! === Absorbing boundary condition =============================================
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: abc
! ==============================================================================
! === Optimization =============================================================
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: tdx, tdy
! ==============================================================================
#ifndef CARTESIAN
! === Density ==================================================================
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: m_rhoC
! ==============================================================================
#endif
! === Arrival time =============================================================
   integer(kind=4), allocatable, dimension(:,:) :: arrivedat
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: arrival_time
! ==============================================================================
#ifdef BANKFILE
   integer(kind=4), allocatable, dimension(:,:) :: ir
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: btx, bty
   integer(kind=4), allocatable, dimension(:,:) :: brokenx, brokeny
#endif
end type wave_arrays
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
type loading
   integer(kind=4) :: N_X, N_Y ! Size of FFT-box
#ifndef __SX__
   integer(kind=8), allocatable, dimension(:) :: &
      xplan_forward, yplan_forward, xplan_backward, yplan_backward ! FFTW3 plan
#else
   integer(kind=4), allocatable, dimension(:,:) :: ifax_x, ifax_y
   real(kind=8), allocatable, dimension(:,:) :: trigs_x, trigs_y
   real(kind=8), allocatable, dimension(:,:) :: work_x
   complex(kind=8), allocatable, dimension(:,:) :: work_y
#endif
   complex(kind=8), allocatable, dimension(:,:) :: green_out_Z
   complex(kind=8), allocatable, dimension(:,:) :: xfftbuf, yfftbuf
   real(kind=8), allocatable, dimension(:,:) :: realbuf
   real(kind=8), allocatable, dimension(:,:) :: defZmap, defZmap1
! === Elastic loading with interpolation =======================================
   real(kind=8), allocatable, dimension(:,:) :: delta
! ==============================================================================
#ifdef MPI
   integer(kind=4) :: nx0, ny0, nx1, ny1, nx2, ny2, ibias, jbias
   real(kind=8), allocatable, dimension(:) :: sendbuf1, recvbuf1
   complex(kind=8), allocatable, dimension(:) :: sendbuf2, recvbuf2
   integer(kind=4), allocatable, dimension(:) :: sendcounts1, sdispls1, recvcounts1, rdispls1
   integer(kind=4), allocatable, dimension(:) :: sendcounts2, sdispls2, recvcounts2, rdispls2
#endif
! === OpenMP ===================================================================
   integer(kind=4), allocatable, dimension(:) :: its2, ite2, itn2 ! N_X, nx2
   integer(kind=4), allocatable, dimension(:) :: jts1, jte1, jtn1 ! N_Y, ny1
! ==============================================================================
end type loading
! ==============================================================================
#endif

type depth_arrays
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dx
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dy
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dz
#ifdef BANKFILE
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dx_old
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dy_old

   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dxbx
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dyby
#endif
end type depth_arrays

type boundary_arrays
   real(kind=REAL_BYTE), allocatable, dimension(:) :: north
   real(kind=REAL_BYTE), allocatable, dimension(:) :: east
   real(kind=REAL_BYTE), allocatable, dimension(:) :: south
   real(kind=REAL_BYTE), allocatable, dimension(:) :: west
end type boundary_arrays

! === Don't repeat allocate/deallocate! ========================================
#ifdef MPI
type edges
   ! for exchange_edges
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fne, fnb ! Flux North Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fse, fsb ! Flux South Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hne, hnb ! Height North Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hse, hsb ! Height South Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fee, feb ! Flux East Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fwe, fwb ! Flux West Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hee, heb ! Height East Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hwe, hwb ! Height West Edge/Buf

   ! for exchange_edges_dz
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dne, dnb ! Depth North Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dse, dsb ! Depth South Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dee, deb ! Depth East Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: dwe, dwb ! Depth West Edge/Buf

   ! exchange_edges_disp_f[xy]
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fne_d, fnb_d ! Flux North Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fse_d, fsb_d ! Flux South Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fee_d, feb_d ! Flux East Edge/Buf
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: fwe_d, fwb_d ! Flux West Edge/Buf

   ! exchange_edges_wod
   integer(kind=4), allocatable, dimension(:,:) :: wne_d, wnb_d ! Wod North Edge/Buf
   integer(kind=4), allocatable, dimension(:,:) :: wse_d, wsb_d ! Wod South Edge/Buf
   integer(kind=4), allocatable, dimension(:,:) :: wee_d, web_d ! Wod East Edge/Buf
   integer(kind=4), allocatable, dimension(:,:) :: wwe_d, wwb_d ! Wod West Edge/Buf
end type edges
#endif
! ==============================================================================

type data_grids
   type(grid_info) :: parent
   type(grid_info) :: my
   type(wave_arrays) :: wave_field
   type(depth_arrays) :: depth_field
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: ts_field
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: zz
#ifndef REAL_DBLE
   real(kind=8), allocatable, dimension(:,:) :: zz_dp
#endif
   type(boundary_arrays) :: ubnd
   type(boundary_arrays) :: hbnd
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hzmax
#ifdef HZMINOUT
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: hzmin
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: vmax
! ==============================================================================
   integer(kind=4), allocatable, dimension(:,:) :: wod_flags
   character(len=256) :: wod_file
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: wod_field
   character(len=256) :: bcf_file
#ifdef BANKFILE
   character(len=256) :: bank_file
#endif
   real(kind=REAL_BYTE), allocatable, dimension(:,:) :: bcf_field
   type(interp_info) :: fxo
   type(interp_info) :: fyo
   type(interp_info) :: fxi
   type(interp_info) :: fyi
   type(interp_info) :: hzi
   type(interp_info) :: dzi
! === copy2coarse for hz =======================================================
   type(interp_info) :: hzo
! ==============================================================================
! === Support multiple ruptures. ===============================================
   integer(kind=4) :: nrupt, irupt, jrupt ! nrupt = number of rupture time steps
   character(len=256), allocatable, dimension(:) :: ruptgrd
! ==============================================================================
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
   type(loading) :: loading
! ==============================================================================
#endif
! === Don't repeat allocate/deallocate! ========================================
#ifdef MPI
   type(edges) :: edges
#endif
! ==============================================================================
end type data_grids

!*** new RWG DEFINITION FOR TIDE GAUGE OUTPUT ***
type tgsrwg
   real(kind=REAL_BYTE) :: geolat = 0.0d0
   real(kind=REAL_BYTE) :: geolon = 0.0d0 ! location and water depth in geographic coordinates -180/180/-90/90
   real(kind=REAL_BYTE) :: mcolat = 0.0d0
   real(kind=REAL_BYTE) :: mcolon = 0.0d0
   integer(kind=4) :: ig = 0
   integer(kind=4) :: ilon = 0 ! grid point location
   integer(kind=4) :: ilat = 0
   real(kind=REAL_BYTE) :: z = 0.0d0 ! water depth at this location
   real(kind=REAL_BYTE) :: center_lat = 0.0d0, center_lon = 0.0d0 ! center of this array
   real(kind=REAL_BYTE) :: offset = 0.0d0, az = 0.0d0, baz = 0.0d0 ! for arrays this is the distance in km from center of array
   real(kind=REAL_BYTE) :: dt = 0.0d0 ! sampling rate for this site
   integer(kind=4) :: nt = 0 ! number of points
   integer(kind=4) :: number ! tgs number in input file
   character(len=16) :: id = '' ! identifier
end type tgsrwg

!*** DEFINITION FOR TIDE GAUGE OUTPUT ***
type tgs
   integer(kind=4) :: ista
   real(kind=REAL_BYTE) :: geolat, geolon ! location and water depth in geographic coordinates -180/180/-90/90
   integer(kind=4) :: ilat, ilon ! grid point location
   real(kind=REAL_BYTE) :: z ! water depth at this location
   real(kind=REAL_BYTE) :: center_lat, center_lon ! center of this array
   real(kind=REAL_BYTE) :: offset, az, baz ! for arrays this is the distance in km from center of array
   real(kind=REAL_BYTE) :: dt ! sampling rate for this site
   integer(kind=4) :: nt ! number of points
   character(len=16) :: id ! identifier
end type tgs

#ifndef __SX__
contains

   subroutine read_grid_info(gfile,dg,ng)
      character(len=128), intent(in) :: gfile
      type(data_grids), allocatable, dimension(:), intent(out) :: dg
      integer(kind=4), intent(out) :: ng
      character(len=1024) :: buf
      integer(kind=4) :: ig
#ifdef MPI
      integer(kind=4) :: ierr = 0
#endif
#ifdef MULTI
! === Split Dir ================================================================
      character(len=128) :: tmp_filename
! ==============================================================================
#endif

! === Split Dir ================================================================
!     open(1,file=trim(gfile),action='read',status='old',form='formatted',err=100)
#ifndef MULTI
      open(1,file=trim(gfile),action='read',status='old',form='formatted',err=100)
#else
      tmp_filename = trim(input_dirname) // trim(gfile)
      open(1,file=trim(tmp_filename),action='read',status='old',form='formatted',err=100)
#endif
! ==============================================================================

      ng = 0
      do while(.true.)
         read(1,'(a)',end=101) buf
         ng = ng + 1
      end do

101   allocate(dg(ng))
      rewind(1)

      do ig = 1, ng
         read(1,'(a)') buf
         read(buf,*,end=104) dg(ig)%my%base_name, dg(ig)%parent%base_name, dg(ig)%my%linear_flag, dg(ig)%my%bath_file, &
                             dg(ig)%my%disp_file
         read(buf,*,end=105) dg(ig)%my%base_name, dg(ig)%parent%base_name, dg(ig)%my%linear_flag, dg(ig)%my%bath_file, &
                             dg(ig)%my%disp_file, dg(ig)%wod_file
         read(buf,*,end=106) dg(ig)%my%base_name, dg(ig)%parent%base_name, dg(ig)%my%linear_flag, dg(ig)%my%bath_file, &
                             dg(ig)%my%disp_file, dg(ig)%wod_file, dg(ig)%bcf_file
#ifdef BANKFILE
         read(buf,*,end=107) dg(ig)%my%base_name, dg(ig)%parent%base_name, dg(ig)%my%linear_flag, dg(ig)%my%bath_file, &
                             dg(ig)%my%disp_file, dg(ig)%wod_file, dg(ig)%bcf_file, dg(ig)%bank_file
#endif
         cycle
104      dg(ig)%my%disp_file = 'NO_DISPLACEMENT_FILE_GIVEN'
105      dg(ig)%wod_file     = 'NO_WETORDRY_FILE_GIVEN'
106      dg(ig)%bcf_file     = 'NO_FRICTION_FILE_GIVEN'
#ifdef BANKFILE
107      dg(ig)%bank_file    = 'NO_BANK_FILE_GIVEN'
#endif
      end do
#ifdef MULTI
! === Split Dir ================================================================
      do ig = 1, ng
         dg(ig)%my%bath_file = trim(input_dirname) // trim(dg(ig)%my%bath_file)
         if(dg(ig)%my%disp_file /= 'NO_DISPLACEMENT_FILE_GIVEN') then
            dg(ig)%my%disp_file = trim(input_dirname) // trim(dg(ig)%my%disp_file)
         end if
         if(dg(ig)%wod_file /= 'NO_WETORDRY_FILE_GIVEN') then
            dg(ig)%wod_file = trim(input_dirname) // trim(dg(ig)%wod_file)
         end if
         if(dg(ig)%bcf_file /= 'NO_FRICTION_FILE_GIVEN') then
            dg(ig)%bcf_file = trim(input_dirname) // trim(dg(ig)%bcf_file)
         end if
#ifdef BANKFILE
         if(dg(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
            dg(ig)%bank_file = trim(input_dirname) // trim(dg(ig)%bank_file)
         end if
#endif
      end do
! ==============================================================================
#endif

      return

100   write(0,'(a,a)') 'read_grid_info(): cant open file ', trim(gfile)
#ifndef MPI
      stop
#else
      ierr = 1
      call fatal_error(201)
#endif
   end subroutine read_grid_info
#endif

end module mod_grid
