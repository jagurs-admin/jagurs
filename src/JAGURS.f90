#ifdef TIMER
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
program JAGURS
#ifdef MPI
   use mpi
#else
#ifdef MULTI
   use mpi
#endif
#endif
#ifdef MULTI
   use mod_multi
#endif
   use mod_grid
   use mod_params
   use mod_tgs
   use mod_rwg
   use mod_nest
   use mod_newsubs
#ifdef TIMER
   use mod_timer
#endif
#ifdef MPI
   use mod_mpi_fatal
#endif
! === Support restart ==========================================================
   use mod_restart
! ==============================================================================
#ifdef NCDIO
   use mod_ncdio
#endif
#ifdef MPI
#ifdef USE_ALLTOALLV
#ifdef A2A3D
   use mod_a2a3d
#ifdef __FUJITSU
   use extfjmpi
#endif
#endif
#endif
#endif
! === Support truncation =======================================================
   use mod_truncation
! ==============================================================================
! === For MRI ==================================================================
   use mod_init_disp_gaussian
! ==============================================================================
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
   use mod_loading
! ==============================================================================
! === Density ==================================================================
   use mod_density
! ==============================================================================
#endif
! === SINWAVE ==================================================================
   use mod_init_disp_sinwave, only : specify_sinwave_params, hrise_rwg_sin
! ==============================================================================
! === Don't repeat allocate/deallocate! ========================================
#ifdef MPI
   use mod_mpi, only : allocate_edges, deallocate_edges
#endif
! ==============================================================================
! === Displacement =============================================================
   use mod_displacement, only : displacement_initialize,        &
                                displacement_calc_displacement, &
                                displacement_apply_kj_filter,   &
#ifndef CARTESIAN
                                displacement_calc_h0_lat1,      &
#else
                                displacement_calc_h0,           &
#endif
#ifndef MPI
                                displacement_finalize
#else
                                displacement_finalize,          &
                                displacement_mpi_initialize
#endif
! ==============================================================================
! === Initial displacement of child domains is given by interpolation. =========
   use mod_interpolation, only : interp2fine_init_disp
#ifdef MPI
   use mod_interpolation, only : interpolation_mpi_initialize
#endif
! ==============================================================================
#if defined(MPI) && defined(ONEFILE)
   use mod_onefile, only : onefile_initialize, onefile_setparams
#endif
#ifdef BANKFILE
   use mod_bank, only : read_bank_file, update_bathymetory, &
                        save_dxdy_old, update_btxbty
#endif
#ifdef DUMP1D
   use mod_dump1d, only : dump1d_initialize, dump1d, dump1d_finalize
#endif
   implicit none

   type(data_grids), allocatable, target, dimension(:) :: dgrid
   integer(kind=4) :: ig, ngrid, pid

   type(wave_arrays), pointer :: wave_field
   type(depth_arrays), pointer :: depth_field
   type(boundary_arrays), pointer :: ubnd, hbnd
   real(kind=REAL_BYTE), pointer, dimension(:,:) :: wod_field, bcf_field, ts_field, zz, hzmax
#ifdef HZMINOUT
   real(kind=REAL_BYTE), pointer, dimension(:,:) :: hzmin
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
#ifndef SKIP_MAX_VEL
   real(kind=REAL_BYTE), pointer, dimension(:,:) :: vmax
#endif
! ==============================================================================
   integer(kind=4), pointer, dimension(:,:) :: wod_flags
   character(len=256), pointer :: base, file_name_bathymetry, displacement_file_name, wod_file_name, bcf_file_name

   integer(kind=4) :: niz, njz
   real(kind=REAL_BYTE) :: dxdy, mlon0, mlat0, th0, dth
   integer(kind=4) :: linear_flag ! do linear calculation for any nonzero number else 0 gives nonlinear

   integer(kind=4) :: istep, nstep
   integer(kind=4) :: nsta ! total number of tide gauge stations
   type(tgsrwg), allocatable, dimension(:) :: mytgs ! tide gauge data structure see "mod_grid.f90"
   real(kind=REAL_BYTE) :: radmin
   real(kind=REAL_BYTE) :: t
#ifndef NCDIO
   character(len=512) :: str
#else
#ifdef CONV_CHECK
   character(len=512) :: str
#endif
#endif

! === Support multiple ruptures. ===============================================
   integer(kind=4), pointer :: nrupt, irupt, jrupt
   character(len=256), pointer, dimension(:) :: ruptgrd
! ==============================================================================

   integer(kind=4) :: i, j
   integer(kind=4) :: error
#ifdef __SX__
   character(len=1024) :: buf
   integer(kind=4) :: ista
#ifndef CARTESIAN
   real(kind=REAL_BYTE) :: colat, colon
#endif
   real(kind=REAL_BYTE) :: mindh
   integer(kind=4) :: ilon, ilat
#ifdef MPI
   integer(kind=4) :: kg, ksta
#endif
#endif
#ifdef MPI
   ! MPI
   integer(kind=4) :: nprocs = 0
   integer(kind=4) :: myrank = 0
   integer(kind=4) :: ierr = 0
!  integer(kind=4) :: npx = 0 ! in mod_params
!  integer(kind=4) :: npy = 0 ! in mod_params
   integer(kind=4) :: rankx = 0
   integer(kind=4) :: ranky = 0
   integer(kind=4) :: nbx = 0
   integer(kind=4) :: nby = 0
   integer(kind=4) :: ix = 0
   integer(kind=4) :: iy = 0
   integer(kind=4) :: totalNx = 0
   integer(kind=4) :: totalNy = 0
   integer(kind=4) :: has_boundary = 0
   integer(kind=4) :: ixend = 0
   integer(kind=4) :: iyend = 0
   real(kind=REAL_BYTE) :: glon0 =0.0d0
   real(kind=REAL_BYTE) :: glat0 =0.0d0
   integer(kind=4) :: kx = 0
   integer(kind=4) :: ky = 0
   integer(kind=4) :: kxend = 0
   integer(kind=4) :: kyend = 0
   character(len=128) :: suffix = ''
#endif
#ifndef SINGLE_TGS
   character(len=128), allocatable, dimension(:) :: tgsfiles
   integer(kind=4), allocatable, dimension(:) :: tgsfp
#else
   character(len=128) :: tgsfiles
   integer(kind=4) :: tgsfp
   character(len=8) :: tgsnum
#endif
   integer(kind=4) :: tgs_nx
   integer(kind=4) :: tgs_ny
   type(wave_arrays), pointer :: tgs_wfld
!  integer(kind=4) :: tgs_step
!  real(kind=REAL_BYTE) :: tgs_time
   real(kind=REAL_BYTE) :: tgs_lon
   real(kind=REAL_BYTE) :: tgs_lat
   real(kind=REAL_BYTE) :: tgs_z
   integer(kind=4) :: tgs_i
   integer(kind=4) :: tgs_j
   real(kind=REAL_BYTE) :: tgs_fx
   real(kind=REAL_BYTE) :: tgs_fy
   real(kind=REAL_BYTE) :: tgs_hz
!  integer(kind=4) :: tgstxtout = 1 ! in mod_params
!  character(len=128) :: tgstxtoutfile = 'tgs' ! in mod_params
   logical :: tgs_opened
#ifdef MPI
! === Separate stdout into for each process. ===================================
   character(len=128) :: stdout = 'stdout'
! ==============================================================================
#endif
! === Support restart ==========================================================
   integer(kind=4) :: istart
   character(len=256) :: restart_file_name
! ==============================================================================
#ifdef DIROUT
   character(len=128) :: dirname
#endif
#ifdef CONV_CHECK
   integer(kind=4) :: conv_step
#endif
! === Support truncation =======================================================
   logical :: trunc_flag = .false.
! ==============================================================================
#ifdef MPI
#ifdef USE_ALLTOALLV
#ifdef A2A3D
   integer(kind=4) :: a2a3d_num_handler
   integer(kind=4), dimension(3) :: decomp
#ifndef __FUJITSU
   logical :: a2a3d_decomp_exists
   integer(kind=4) :: a2a3d_x, a2a3d_y, a2a3d_z
   namelist /a2a3d_decomp/ a2a3d_x, a2a3d_y, a2a3d_z
#endif
#endif
#endif
#endif
! === For MRI ==================================================================
   character(len=256) :: gaussian_file_name = 'gaussian'
! ==============================================================================
#ifdef __SX__
#ifdef MULTI
! === Split Dir ================================================================
   character(len=128) :: tmp_filename
! ==============================================================================
#endif
#endif
! === SINWAVE ==================================================================
   character(len=256) :: sinwave_file_name = 'sinwave'
! ==============================================================================
#ifdef MPI
   integer(kind=4) :: ist, ien, jst, jen
#endif
#ifdef CARTESIAN
   real(kind=REAL_BYTE) :: plat0, clat0
#endif
! === Displacement =============================================================
#ifndef CARTESIAN
   real(kind=8) :: h0, lat1
#else
   real(kind=8) :: h0
#endif
! ==============================================================================
#ifdef PIXELIN
   integer(kind=4) :: nxorg, nyorg
#endif
#ifdef NFSUPPORT
   integer(kind=4) :: formatid
#endif
   TIMER_START('All')

#ifdef MULTI
   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, g_nprocs, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, g_myrank, ierr)

   call get_arguments()

#ifndef MPI
   member_id = g_myrank
#else
   member_id = g_nprocs/num_members
   member_id = g_myrank/member_id
   call MPI_Comm_Split(MPI_COMM_WORLD, member_id, g_myrank, MPI_MEMBER_WORLD, ierr)
#endif

! === Split Dir ================================================================
   write(suffix,'(i6.6)') member_ids(member_id+1)
   input_dirname = trim(input_dirname) // trim(suffix) // '/'
! ==============================================================================
   write(suffix,'(a,i6.6)') '.', member_ids(member_id+1)
   members_dir = 'member' // trim(suffix) // '/'
   gaussian_file_name = trim(gaussian_file_name) // trim(suffix)

   if(restart == 0) then
      write(command,'(a,a)') 'mkdir -p ', trim(members_dir)
      call system(trim(command))
   end if

#ifndef MPI
   stdout = trim(members_dir) // trim(stdout)
   open(6,file=stdout,action='write',status='replace',form='formatted')
#endif
#endif

#ifdef MPI
#ifndef MULTI
   call MPI_Init(ierr)
#endif
#endif
   TIMER_START('getpar')
   call getpar()
   TIMER_STOP('getpar')
! === Support truncation =======================================================
   if(max_time_i /= 0) call start_trunc()
! ==============================================================================
#ifdef MPI
#ifndef MULTI
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
#else
   call MPI_Comm_size(MPI_MEMBER_WORLD, nprocs, ierr)
   call MPI_Comm_rank(MPI_MEMBER_WORLD, myrank, ierr)
#endif
#endif
#if defined(MPI) && defined(ONEFILE)
   TIMER_START('onefile_initialize')
   call onefile_initialize(nprocs, myrank)
   TIMER_STOP('onefile_initialize')
#endif

   nstep = nint(tend/dt) + 1

   !** Burbidge Added another sanity check **
   if(nstep <= 0) then
      write(0,'(a,i0)') 'Number of steps is invalid. nstep=', nstep
#ifndef MPI
      stop
#else
      call fatal_error(107)
#endif
   end if

   radmin = acos(-1.0d0)/(180.0d0)

   TIMER_START('read_grid_info')
#ifndef __SX__
   call read_grid_info(gridfile,dgrid,ngrid)
#else
! === Split Dir ================================================================
!  open(1,file=trim(gridfile),action='read',status='old',form='formatted',err=300)
#ifndef MULTI
   open(1,file=trim(gridfile),action='read',status='old',form='formatted',err=300)
#else
   tmp_filename = trim(input_dirname) // trim(gridfile)
   open(1,file=trim(tmp_filename),action='read',status='old',form='formatted',err=300)
#endif
! ==============================================================================

   ngrid = 0
   do while(.true.)
      read(1,'(a)',end=301) buf
      ngrid = ngrid + 1
   end do

301 allocate(dgrid(ngrid))
   rewind(1)

   do ig = 1, ngrid
      read(1,'(a)') buf
      read(buf,*,end=304) dgrid(ig)%my%base_name, dgrid(ig)%parent%base_name, dgrid(ig)%my%linear_flag, dgrid(ig)%my%bath_file, &
                          dgrid(ig)%my%disp_file
      read(buf,*,end=305) dgrid(ig)%my%base_name, dgrid(ig)%parent%base_name, dgrid(ig)%my%linear_flag, dgrid(ig)%my%bath_file, &
                          dgrid(ig)%my%disp_file, dgrid(ig)%wod_file
      read(buf,*,end=306) dgrid(ig)%my%base_name, dgrid(ig)%parent%base_name, dgrid(ig)%my%linear_flag, dgrid(ig)%my%bath_file, &
                          dgrid(ig)%my%disp_file, dgrid(ig)%wod_file, dgrid(ig)%bcf_file
#ifdef BANKFILE
      read(buf,*,end=307) dgrid(ig)%my%base_name, dgrid(ig)%parent%base_name, dgrid(ig)%my%linear_flag, dgrid(ig)%my%bath_file, &
                          dgrid(ig)%my%disp_file, dgrid(ig)%wod_file, dgrid(ig)%bcf_file, dgrid(ig)%bank_file
#endif
      cycle
304   dgrid(ig)%my%disp_file = 'NO_DISPLACEMENT_FILE_GIVEN'
305   dgrid(ig)%wod_file     = 'NO_WETORDRY_FILE_GIVEN'
306   dgrid(ig)%bcf_file     = 'NO_FRICTION_FILE_GIVEN'
#ifdef BANKFILE
307   dgrid(ig)%bank_file    = 'NO_BANK_FILE_GIVEN'
#endif
   end do
#ifdef MULTI
! === Split Dir ================================================================
   do ig = 1, ngrid
      dgrid(ig)%my%bath_file = trim(input_dirname) // trim(dgrid(ig)%my%bath_file)
      if(dgrid(ig)%my%disp_file /= 'NO_DISPLACEMENT_FILE_GIVEN') then
         dgrid(ig)%my%disp_file = trim(input_dirname) // trim(dgrid(ig)%my%disp_file)
      end if
      if(dgrid(ig)%wod_file /= 'NO_WETORDRY_FILE_GIVEN') then
         dgrid(ig)%wod_file = trim(input_dirname) // trim(dgrid(ig)%wod_file)
      end if
      if(dgrid(ig)%bcf_file /= 'NO_FRICTION_FILE_GIVEN') then
         dgrid(ig)%bcf_file = trim(input_dirname) // trim(dgrid(ig)%bcf_file)
      end if
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         dgrid(ig)%bank_file = trim(input_dirname) // trim(dgrid(ig)%bank_file)
      end if
#endif
   end do
! ==============================================================================
#endif
#endif
   TIMER_STOP('read_grid_info')
#ifdef MPI
#ifdef USE_ALLTOALLV
#ifdef A2A3D
#ifndef SINGLE_A2A
   a2a3d_num_handler = (ngrid - 1)*9
#else
   a2a3d_num_handler = (ngrid - 1)*4
#endif
#ifdef __FUJITSU
! === DEBUG for A2A3D & MULTI on K-computer on 2016/02/08 ======================
#ifndef MULTI
! ==============================================================================
   ierr = FJMPI_Topology_get_shape(decomp(1), decomp(2), decomp(3))
   do i = 1, 3
     if(decomp(i) == 0) decomp(i) = 1
   end do
! === DEBUG for A2A3D & MULTI on K-computer on 2016/02/08 ======================
#else
   inquire(file='a2a3d_decomp', exist=a2a3d_decomp_exists)
   if(a2a3d_decomp_exists) then
      if(myrank == 0) write(0,'(a)') 'A2A3D decomposition file a2a3d_decomp exists!!!'
      open(1,file='a2a3d_decomp',action='read',status='old',form='formatted')
      read(1,a2a3d_decomp)
      close(1)
      decomp(1) = a2a3d_x
      decomp(2) = a2a3d_y
      decomp(3) = a2a3d_z
   else
      if(myrank == 0) write(0,'(a)') 'A2A3D decomposition file a2a3d_decomp does NOT exist!!!'
      if(myrank == 0) write(0,'(a)') 'Default decomposition is used!!!'
      call decomp_default(nprocs, decomp)
   end if
#endif
! ==============================================================================
#else
   inquire(file='a2a3d_decomp', exist=a2a3d_decomp_exists)
   if(a2a3d_decomp_exists) then
      if(myrank == 0) write(0,'(a)') 'A2A3D decomposition file a2a3d_decomp exists!!!'
      open(1,file='a2a3d_decomp',action='read',status='old',form='formatted')
      read(1,a2a3d_decomp)
      close(1)
      decomp(1) = a2a3d_x
      decomp(2) = a2a3d_y
      decomp(3) = a2a3d_z
   else
      if(myrank == 0) write(0,'(a)') 'A2A3D decomposition file a2a3d_decomp does NOT exist!!!'
      if(myrank == 0) write(0,'(a)') 'Default decomposition is used!!!'
      call decomp_default(nprocs, decomp)
   end if
#endif
   if(myrank == 0) write(0,'(a,3i6)') 'NOTE: A2A3D decomposition: ', decomp
   TIMER_START('A2A3D_init')
   call A2A3D_init(nprocs, myrank, decomp, a2a3d_num_handler)
   TIMER_STOP('A2A3D_init')
#endif
#endif
   write(suffix,'(a,i6.6)') '.', myrank
! === Separate stdout into for each process. ===================================
   stdout = trim(stdout) // trim(suffix)
#ifdef MULTI
   stdout = trim(members_dir) // trim(stdout)
#endif
   open(6,file=stdout,action='write',status='replace',form='formatted')
! ==============================================================================
#ifndef ONEFILE
   call add_suffix_to_grid_filenames(dgrid, ngrid, suffix)
#endif

   if(nprocs /= npx*npy) then
      write(0,'(a,i0,a,i0,a,i0,a,i0,a)') 'FATAL Error : Number of MPI processes (nprocs=', nprocs, &
         ') is not equal nporcx * nproxy (', npx, ' * ', npy, ' = ', npx*npy, ') !'
      call fatal_error(111)
   end if
   rankx = mod(myrank, npx)
   ranky = myrank/npx
#endif
! === For MRI ==================================================================
   if(init_disp_gaussian == 1) then
      TIMER_START('specify_gaussian_params')
      call specify_gaussian_params(gaussian_file_name)
      TIMER_STOP('specify_gaussian_params')
   end if
! ==============================================================================
! === SINWAVE ==================================================================
#ifndef MULTI
   if(init_disp_sinwave == 1) then
      TIMER_START('specify_sinwave_params')
      call specify_sinwave_params(sinwave_file_name)
      TIMER_STOP('specify_sinwave_params')
   end if
#else
   if(init_disp_sinwave == 1) then
      TIMER_START('specify_sinwave_params')
      call specify_sinwave_params(input_dirname, sinwave_file_name)
      TIMER_STOP('specify_sinwave_params')
   end if
#endif
! ==============================================================================
! === Displacement =============================================================
#ifdef MPI
   if(apply_kj_filter == 1) then
      TIMER_START('displacement_mpi_initialize')
      call displacement_mpi_initialize(nprocs, myrank, npx, npy, rankx, ranky)
      TIMER_STOP('displacement_mpi_initialize')
   end if
#endif
! ==============================================================================
! === Initial displacement of child domains is given by interpolation. =========
#ifdef MPI
! === Elastic loading with interpolation =======================================
!  if(init_disp_interpolation == 1) then
   if((init_disp_interpolation == 1) .or. (elastic_loading_interpolation == 1)) then
! ==============================================================================
      TIMER_START('interpolation_mpi_initialize')
      call interpolation_mpi_initialize(nprocs, myrank, npx, npy, rankx, ranky)
      TIMER_STOP('interpolation_mpi_initialize')
   end if
#endif
! ==============================================================================
#ifndef _WIN32
! === Output config! ===========================================================
   TIMER_START('putpar')
#ifndef MPI
   call putpar(gridfile, ngrid, dgrid)
#else
   call putpar(nprocs, gridfile, ngrid, dgrid)
#endif
   TIMER_STOP('putpar')
! ==============================================================================
#endif

   ! thomas - code to accommodate multiple ruptures
   !          in the multiple rupture case dgrid[ig].my.disp_file is interpreted as
   !          the name of a text file (a rupture list file) that lists the names
   !          of the multiple rupture files for grid ig
   !          tau is interpreted as the time step between ruptures
   if(multrupt == 1) then
#ifdef MPI
      ierr = 0
#endif
      write(6,'(a)') 'Multiple ruptures'
! === Support multiple ruptures. ===============================================
      do ig = 1, ngrid
         nrupt   => dgrid(ig)%nrupt
         irupt   => dgrid(ig)%irupt
! === Displacement =============================================================
         if((ig /= 1) .and. (init_disp_interpolation == 1)) then
            if(trim(dgrid(ig)%my%disp_file) == 'NO_DISPLACEMENT_FILE_GIVEN') then
               dgrid(ig)%nrupt = dgrid(1)%nrupt
               dgrid(ig)%irupt = dgrid(1)%irupt
               allocate(dgrid(ig)%ruptgrd(1))
               cycle
            end if
         end if
! ==============================================================================
         ! there must be a rupture list file for each domain, else error
         open(1,file=trim(dgrid(ig)%my%disp_file),action='read',status='old',form='formatted',err=100)
         ! count the rupture file names for each domain, this tells us the number of ruptures
         ! temporarily allocate some space to ruptgrd, just for the count
         allocate(dgrid(ig)%ruptgrd(1))
         nrupt = 0
         do while(.true.)
            read(1,'(a)',end=101) dgrid(ig)%ruptgrd(1)
            nrupt = nrupt + 1
         end do
101      continue
         close(1)
         write(6,'(a,i0,a,f0.6)') 'Number of ruptures=', nrupt, ' tau=', tau
         ! reallocate ruptgrd since we now know nrupt
         deallocate(dgrid(ig)%ruptgrd)
         allocate(dgrid(ig)%ruptgrd(nrupt))
         ruptgrd => dgrid(ig)%ruptgrd
         ! Fill in the rupture file names
         if(trim(dgrid(ig)%my%disp_file) == 'NO_DISPLACEMENT_FILE_GIVEN') then
            do irupt = 1, nrupt
               ruptgrd(irupt) = 'NO_DISPLACEMENT_FILE_GIVEN'
               write(6,'(a,i0,a,a)') 'grid: ', ig, ' rupture file: ', trim(ruptgrd(irupt))
            end do
         else
            open(1,file=trim(dgrid(ig)%my%disp_file),action='read',status='old',form='formatted',err=100)
            do irupt = 1, nrupt
               read(1,'(a)',err=100) ruptgrd(irupt)
#ifdef MULTI
               ruptgrd(irupt) = trim(input_dirname) // trim(ruptgrd(irupt))
#endif
#ifdef MPI
#ifndef ONEFILE
! === Displacement =============================================================
! === DEBUG for MPI & multiple rupture & interpolation! ========================
!              if(init_disp_interpolation == 0) then
               if(init_disp_fault == 0) then
! ==============================================================================
! ==============================================================================
               ruptgrd(irupt) = trim(ruptgrd(irupt)) // trim(suffix)
! === Displacement =============================================================
               end if
! ==============================================================================
#endif
#endif
               write(6,'(a,i0,a,a)') 'grid: ', ig, ' rupture file: ', trim(ruptgrd(irupt))
            end do
            close(1)
         end if
      end do
   else
      do ig = 1, ngrid
         nrupt   => dgrid(ig)%nrupt
         irupt   => dgrid(ig)%irupt
         nrupt = 1
         allocate(dgrid(ig)%ruptgrd(1))
         ruptgrd => dgrid(ig)%ruptgrd
         ruptgrd(1) = dgrid(ig)%my%disp_file
      end do
   end if
! ==============================================================================
#ifdef MPI
110 if(ierr > 0) call fatal_error(108)
#endif

#if defined(MPI) && defined(ONEFILE)
   if(myrank == 0) then
#endif
   write(6,'(a,a)') trim(program_name), ': Reading grid data...'
#if defined(MPI) && defined(ONEFILE)
   end if
#endif

   do ig = 1, ngrid
#if defined(MPI) && defined(ONEFILE)
      if(myrank == 0) then
#endif
      write(6,'(/,2x,a,i0)') '*** grid ', ig

      !*** read in bathymetry GMT netcdf grid file ***
      !*** first read the header values ***
      TIMER_START('read_bathymetry_gmt_grdhdr')
#ifndef PIXELIN
#ifndef NFSUPPORT
      call read_bathymetry_gmt_grdhdr(dgrid(ig)%my%bath_file,niz,njz,dxdy,mlon0,mlat0)
#else
      call read_bathymetry_gmt_grdhdr(dgrid(ig)%my%bath_file,niz,njz,dxdy,mlon0,mlat0,formatid)
#endif
#else
#ifndef NFSUPPORT
      call read_bathymetry_gmt_grdhdr(dgrid(ig)%my%bath_file,niz,njz,dxdy,mlon0,mlat0,nxorg,nyorg)
#else
      call read_bathymetry_gmt_grdhdr(dgrid(ig)%my%bath_file,niz,njz,dxdy,mlon0,mlat0,nxorg,nyorg,formatid)
#endif
#endif
      TIMER_STOP('read_bathymetry_gmt_grdhdr')
#if defined(MPI) && defined(ONEFILE)
      end if

#ifndef MULTI
      call MPI_Bcast(niz,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(njz,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(dxdy,  1, REAL_MPI,    0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(mlon0, 1, REAL_MPI,    0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(mlat0, 1, REAL_MPI,    0, MPI_COMM_WORLD, ierr)
#ifdef PIXELIN
      call MPI_Bcast(nxorg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(nyorg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
#ifdef NFSUPPORT
      call MPI_Bcast(formatid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
#else
      call MPI_Bcast(niz,   1, MPI_INTEGER, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(njz,   1, MPI_INTEGER, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(dxdy,  1, REAL_MPI,    0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(mlon0, 1, REAL_MPI,    0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(mlat0, 1, REAL_MPI,    0, MPI_MEMBER_WORLD, ierr)
#ifdef PIXELIN
      call MPI_Bcast(nxorg, 1, MPI_INTEGER, 0, MPI_MEMBER_WORLD, ierr)
      call MPI_Bcast(nyorg, 1, MPI_INTEGER, 0, MPI_MEMBER_WORLD, ierr)
#endif
#ifdef NFSUPPORT
      call MPI_Bcast(formatid, 1, MPI_INTEGER, 0, MPI_MEMBER_WORLD, ierr)
#endif
#endif

      totalNx = niz
      nbx = totalNx / npx + min(1, mod(totalNx, npx))
      niz = nbx
      if(rankx == npx - 1) niz = totalNx - nbx*(npx - 1)
      if(rankx <  npx - 1) niz = niz + 1
      if(rankx >  0)       niz = niz + 1
      if(niz < 3) then
         write(6,'(a)')    '############################################################'
         write(6,'(a)')    '# ERROR!!! Num. of procs. on x-direction may be too large. #'
         write(6,'(a)')    '#          Local size must be greater than 3.              #'
         write(6,'(a)')    '#          Please try smaller procx.                       #'
         write(6,'(a)')    '#                                                          #'
         write(6,'(a,i6,a)') '# Domain ID:   ', ig,      '                                      #'
         write(6,'(a,i6,a)') '# Global size: ', totalNx, '                                      #'
         write(6,'(a,i6,a)') '# procx:       ', npx,     '                                      #'
         write(6,'(a,i6,a)') '# Local size:  ', niz,     ' (must be > 3)                        #'
         write(6,'(a)') '############################################################'
         call MPI_Abort(MPI_COMM_WORLD, 999, ierr)
      end if

      totalNy = njz
      nby = totalNy / npy + min(1, mod(totalNy, npy))
      njz = nby
      if(ranky == npy - 1) njz = totalNy - nby*(npy - 1)
      if(ranky <  npy - 1) njz = njz + 1
      if(ranky >  0)       njz = njz + 1
      if(njz < 3) then
         write(6,'(a)')    '############################################################'
         write(6,'(a)')    '# ERROR!!! Num. of procs. on y-direction may be too large. #'
         write(6,'(a)')    '#          Local size must be greater than 3.              #'
         write(6,'(a)')    '#          Please try smaller procy.                       #'
         write(6,'(a)')    '#                                                          #'
         write(6,'(a,i6,a)') '# Domain ID:   ', ig,      '                                      #'
         write(6,'(a,i6,a)') '# Global size: ', totalNy, '                                      #'
         write(6,'(a,i6,a)') '# procy:       ', npy,     '                                      #'
         write(6,'(a,i6,a)') '# Local size:  ', njz,     ' (must be > 3)                        #'
         write(6,'(a)') '############################################################'
         call MPI_Abort(MPI_COMM_WORLD, 999, ierr)
      end if
#endif

      dgrid(ig)%my%nx = niz
      dgrid(ig)%my%ny = njz
      dgrid(ig)%my%dh = dxdy
      dgrid(ig)%my%mlon0 = mlon0
      dgrid(ig)%my%mlat0 = mlat0
#ifdef PIXELIN
      dgrid(ig)%my%nxorg = nxorg
      dgrid(ig)%my%nyorg = nyorg
#endif
#ifdef NFSUPPORT
      dgrid(ig)%my%formatid = formatid
#endif
 
#ifdef MPI
      ! MPI
      !       north 
      !   +---+---+---+
      !   | 0 | 1 | 2 | e
      ! w +---+---+---+ a
      ! e | 3 | 4 | 5 | s
      ! s +---+---+---+ t
      ! t | 6 | 7 | 8 |  
      !   +---+---+---+
      !       south

      dgrid(ig)%my%px = npx
      dgrid(ig)%my%py = npy
      dgrid(ig)%my%rx = rankx
      dgrid(ig)%my%ry = ranky

#ifndef MULTI
      call MPI_Allreduce(niz, totalNx, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(niz, totalNx, 1, MPI_INTEGER, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
      totalNx = totalNx/npy - (npx-1)*2
      dgrid(ig)%my%totalNx = totalNx

#ifndef MULTI
      call MPI_Allreduce(njz, totalNy, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(njz, totalNy, 1, MPI_INTEGER, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
      totalNy = totalNy/npx - (npy-1)*2
      dgrid(ig)%my%totalNy = totalNy

#ifndef MULTI
      call MPI_Allreduce(mlon0, glon0, 1, REAL_MPI, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(mlon0, glon0, 1, REAL_MPI, MPI_MIN, MPI_MEMBER_WORLD, ierr)
#endif
      dgrid(ig)%my%glon0 = glon0

#ifndef MULTI
      call MPI_Allreduce(mlat0, glat0, 1, REAL_MPI, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(mlat0, glat0, 1, REAL_MPI, MPI_MIN, MPI_MEMBER_WORLD, ierr)
#endif
      dgrid(ig)%my%glat0 = glat0

      nbx    = totalNx / npx ! x-size of local region without edges
      if(mod(totalNx, npx) /= 0) nbx = nbx + 1
      nby    = totalNy / npy ! y-size of local region without edges
      if(mod(totalNy, npy) /= 0) nby = nby + 1
      ix     = rankx*nbx + 1 ! x start global index without edges
      iy     = ranky*nby + 1 ! y start global index without edges
      ixend  = min(ix + nbx - 1, totalNx) ! x end global index without edges
      iyend  = min(iy + nby - 1, totalNy) ! y end global index without edges
      kx     = ix - 1 ! with edges
      ky     = iy - 1 ! with edges
      kxend  = ixend + 1 ! with edges
      kyend  = iyend + 1 ! with edges

      has_boundary = 0
      if(ix == 1) then
         has_boundary = ior(has_boundary, WEST_BOUND)
         kx = ix 
      end if
      if(ixend == totalNx) then
         has_boundary = ior(has_boundary, EAST_BOUND)
         kxend = ixend
      end if
      if(iy == 1) then
         has_boundary = ior(has_boundary, NORTH_BOUND)
         ky = iy
      end if
      if(iyend == totalNy) then
         has_boundary = ior(has_boundary, SOUTH_BOUND)
         kyend = iyend
      end if

      dgrid(ig)%my%nbx          = nbx
      dgrid(ig)%my%nby          = nby
      dgrid(ig)%my%ix           = ix
      dgrid(ig)%my%iy           = iy
      dgrid(ig)%my%ixend        = ixend
      dgrid(ig)%my%iyend        = iyend
      dgrid(ig)%my%kx           = kx
      dgrid(ig)%my%ky           = ky
      dgrid(ig)%my%kxend        = kxend
      dgrid(ig)%my%kyend        = kyend
      dgrid(ig)%my%has_boundary = has_boundary
#ifdef ONEFILE
      TIMER_START('onefile_setparams')
      call onefile_setparams(dgrid(ig))
      TIMER_STOP('onefile_setparams')
#endif
! === To avoid numerical error =================================================
#ifndef MULTI
      call MPI_Bcast(dgrid(ig)%my%dh, 1, REAL_MPI, 0, MPI_COMM_WORLD, ierr)
#else
      call MPI_Bcast(dgrid(ig)%my%dh, 1, REAL_MPI, 0, MPI_MEMBER_WORLD, ierr)
#endif
#ifndef CARTESIAN
      dgrid(ig)%my%mlon0 = dgrid(ig)%my%glon0 + (kx - 1)*dgrid(ig)%my%dh*60.0d0
      dgrid(ig)%my%mlat0 = dgrid(ig)%my%glat0 + (ky - 1)*dgrid(ig)%my%dh*60.0d0
#endif
! ==============================================================================
#endif
      !*** allocate and initalize the memory space ***
#ifndef MPI
      allocate(dgrid(ig)%wave_field%fx    (-2:niz+1,-1:njz+1))
      allocate(dgrid(ig)%wave_field%fy    (-1:niz+1,-2:njz+1))
      allocate(dgrid(ig)%wave_field%fx_old(-2:niz+1,-1:njz+1))
      allocate(dgrid(ig)%wave_field%fy_old(-1:niz+1,-2:njz+1))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         allocate(dgrid(ig)%wave_field%btx(-2:niz+1,-1:njz+1))
         allocate(dgrid(ig)%wave_field%brokenx(-2:niz+1,-1:njz+1))
         allocate(dgrid(ig)%wave_field%bty(-1:niz+1,-2:njz+1))
         allocate(dgrid(ig)%wave_field%brokeny(-1:niz+1,-2:njz+1))
      end if
#endif
      if((with_disp == 1) .or. (with_disp == 2 .and. ig /= 1)) then
         allocate(dgrid(ig)%wave_field%yu(0:niz,0:njz))
         allocate(dgrid(ig)%wave_field%yv(0:niz,0:njz))
      end if
#else
      ! fx, fx_old
      ist = -1
      ien = niz + 1
      if(iand(has_boundary, WEST_BOUND) /= 0) ist = -2
      jst = 0
      jen = njz + 1
      if(iand(has_boundary, NORTH_BOUND) /= 0) jst = -1
      allocate(dgrid(ig)%wave_field%fx    (ist:ien,jst:jen))
      allocate(dgrid(ig)%wave_field%fx_old(ist:ien,jst:jen))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         allocate(dgrid(ig)%wave_field%btx(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%brokenx(ist:ien,jst:jen))
      end if
#endif

      ! fy, fy_old
      ist = 0
      ien = niz + 1
      if(iand(has_boundary, WEST_BOUND) /= 0) ist = -1
      jst = -1
      jen = njz + 1
      if(iand(has_boundary, NORTH_BOUND) /= 0) jst = -2
      allocate(dgrid(ig)%wave_field%fy    (ist:ien,jst:jen))
      allocate(dgrid(ig)%wave_field%fy_old(ist:ien,jst:jen))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         allocate(dgrid(ig)%wave_field%bty(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%brokeny(ist:ien,jst:jen))
      end if
#endif
      if((with_disp == 1) .or. (with_disp == 2 .and. ig /= 1)) then
         if(iand(has_boundary, WEST_BOUND)  /= 0) then
            if(iand(has_boundary, NORTH_BOUND)  /= 0) then
               allocate(dgrid(ig)%wave_field%yu(0:niz,0:njz))
               allocate(dgrid(ig)%wave_field%yv(0:niz,0:njz))
            else
               allocate(dgrid(ig)%wave_field%yu(0:niz,njz))
               allocate(dgrid(ig)%wave_field%yv(0:niz,njz))
            end if
         else
            if(iand(has_boundary, NORTH_BOUND)  /= 0) then
               allocate(dgrid(ig)%wave_field%yu(niz,0:njz))
               allocate(dgrid(ig)%wave_field%yv(niz,0:njz))
            else
               allocate(dgrid(ig)%wave_field%yu(niz,njz))
               allocate(dgrid(ig)%wave_field%yv(niz,njz))
            end if
         end if
      end if
#endif
      if((with_disp == 1) .or. (with_disp == 2 .and. ig /= 1)) then
         allocate(dgrid(ig)%wave_field%cu(niz,njz))
         allocate(dgrid(ig)%wave_field%cv(niz,njz))
      end if
#ifndef MPI
      allocate(dgrid(ig)%wave_field%hz    (-1:niz+2,-1:njz+2))
      allocate(dgrid(ig)%wave_field%hz_old(-1:niz+2,-1:njz+2))
      allocate(dgrid(ig)%depth_field%dz   (-1:niz+2,-1:njz+2))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         allocate(dgrid(ig)%wave_field%ir   (-1:niz+2,-1:njz+2))
         allocate(dgrid(ig)%depth_field%dxbx(-1:niz+2,-1:njz+2))
         allocate(dgrid(ig)%depth_field%dyby(-1:niz+2,-1:njz+2))
      end if
#endif
#ifndef CARTESIAN
! === Density ==================================================================
      if(with_density == 1) then
         allocate(dgrid(ig)%wave_field%m_rhoC(-1:niz+2,-1:njz+2))
      end if
! ==============================================================================
#endif
#else
      ist = 0
      ien = niz + 2
      if(iand(has_boundary, WEST_BOUND) /= 0) ist = -1
      jst = 0
      jen = njz + 2
      if(iand(has_boundary, NORTH_BOUND) /= 0) jst = -1
      allocate(dgrid(ig)%wave_field%hz    (ist:ien,jst:jen))
#ifndef NONESTDEBUG
      allocate(dgrid(ig)%wave_field%noi2f (ist:ien,jst:jen))
#endif
      allocate(dgrid(ig)%wave_field%hz_old(ist:ien,jst:jen))
      allocate(dgrid(ig)%depth_field%dz   (ist:ien,jst:jen))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         allocate(dgrid(ig)%wave_field%ir   (ist:ien,jst:jen))
         allocate(dgrid(ig)%depth_field%dxbx(ist:ien,jst:jen))
         allocate(dgrid(ig)%depth_field%dyby(ist:ien,jst:jen))
      end if
#endif
#ifndef CARTESIAN
! === Density ==================================================================
      if(with_density == 1) then
         allocate(dgrid(ig)%wave_field%m_rhoC(ist:ien,jst:jen))
      end if
! ==============================================================================
#endif
#endif
      if(with_abc == 1) then
         if(ig == 1) then
            allocate(dgrid(ig)%wave_field%abc(niz,njz))
            TIMER_START('make_abc')
#ifndef MPI
            call make_abc(dgrid(ig)%wave_field%abc,niz,njz)
#else
            call make_abc(dgrid(ig)%wave_field%abc,niz,njz,          &
                          dgrid(ig)%my%totalNx,dgrid(ig)%my%totalNy, &
                          dgrid(ig)%my%kx,dgrid(ig)%my%ky)
#endif
            TIMER_STOP('make_abc')
         end if
      end if
      allocate(dgrid(ig)%depth_field%dx(niz,njz))
      allocate(dgrid(ig)%depth_field%dy(niz,njz))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         allocate(dgrid(ig)%depth_field%dx_old(niz,njz))
         allocate(dgrid(ig)%depth_field%dy_old(niz,njz))
      end if
#endif
      allocate(dgrid(ig)%wod_field(niz,njz))
      allocate(dgrid(ig)%bcf_field(niz,njz))
      allocate(dgrid(ig)%ts_field(niz,njz))
      allocate(dgrid(ig)%zz(niz,njz))
#ifndef REAL_DBLE
      allocate(dgrid(ig)%zz_dp(niz,njz))
#endif
      allocate(dgrid(ig)%hzmax(niz,njz))
#ifdef HZMINOUT
      allocate(dgrid(ig)%hzmin(niz,njz))
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
#ifndef SKIP_MAX_VEL
      allocate(dgrid(ig)%vmax(niz,njz))
#endif
! ==============================================================================
! === Conversion from flux to velocity should be done right after calc. ========
!     allocate(dgrid(ig)%wod_flags(niz,njz))
#ifndef MPI
      allocate(dgrid(ig)%wod_flags(niz,njz))
#else
      allocate(dgrid(ig)%wod_flags(0:niz+1,0:njz+1))
#endif
! ==============================================================================
! === Arrival time =============================================================
      if(check_arrival_time == 1) then
         allocate(dgrid(ig)%wave_field%arrivedat(niz,njz))
         allocate(dgrid(ig)%wave_field%arrival_time(niz,njz))
      end if
! ==============================================================================
      allocate(dgrid(ig)%hbnd%north(niz))
      allocate(dgrid(ig)%hbnd%east(njz))
      allocate(dgrid(ig)%hbnd%south(niz))
      allocate(dgrid(ig)%hbnd%west(njz))
      allocate(dgrid(ig)%ubnd%north(niz))
      allocate(dgrid(ig)%ubnd%east(njz))
      allocate(dgrid(ig)%ubnd%south(niz))
      allocate(dgrid(ig)%ubnd%west(njz))
! === Optimization =============================================================
      if((with_disp == 1) .or. (with_disp == 2 .and. ig /= 1)) then
         allocate(dgrid(ig)%wave_field%tdx(niz,njz))
         allocate(dgrid(ig)%wave_field%tdy(niz,njz))
      end if
! ==============================================================================

      !*** initialize depth and friction field ***
      dgrid(ig)%depth_field%dx = 0.0d0 ! depth field
      dgrid(ig)%depth_field%dy = 0.0d0 ! depth field
      dgrid(ig)%depth_field%dz = 0.0d0 ! depth field
      dgrid(ig)%bcf_field = 0.0d0 ! friction field
#ifdef BANKFILE
      if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
         dgrid(ig)%depth_field%dx_old = 0.0d0
         dgrid(ig)%depth_field%dy_old = 0.0d0

         dgrid(ig)%wave_field%btx = 0.0d0
         dgrid(ig)%wave_field%bty = 0.0d0
         dgrid(ig)%wave_field%ir = 0
         dgrid(ig)%depth_field%dxbx = 0.0d0
         dgrid(ig)%depth_field%dyby = 0.0d0

         if(broken_rate < 0.0d0) then
            dgrid(ig)%wave_field%brokenx = 1
            dgrid(ig)%wave_field%brokeny = 1
         else
            dgrid(ig)%wave_field%brokenx = 0
            dgrid(ig)%wave_field%brokeny = 0
         end if
      end if
#endif

      !*** set simulation parameters ***
#ifndef CARTESIAN
#ifndef MPI
      dgrid(ig)%my%th0 = dgrid(ig)%my%mlat0/60.0d0*radmin
#else
      dgrid(ig)%my%th0 = REAL_FUNC(dgrid(ig)%my%glat0/60.0d0)*radmin ! global domain value
#endif
#else
      dgrid(ig)%my%th0 = dgrid(ig)%my%mlat0
#endif
      dgrid(ig)%my%dth = dgrid(ig)%my%dh*radmin  !*** dxdy in degrees ***

      !*** set pointers for function calls ***
      base                   => dgrid(ig)%my%base_name
      file_name_bathymetry   => dgrid(ig)%my%bath_file
      displacement_file_name => dgrid(ig)%my%disp_file
      niz                    =  dgrid(ig)%my%nx
      njz                    =  dgrid(ig)%my%ny
#ifdef PIXELIN
      nxorg                  =  dgrid(ig)%my%nxorg
      nyorg                  =  dgrid(ig)%my%nyorg
#endif
      dxdy                   =  dgrid(ig)%my%dh
      mlon0                  =  dgrid(ig)%my%mlon0
      mlat0                  =  dgrid(ig)%my%mlat0
      linear_flag            =  dgrid(ig)%my%linear_flag
      th0                    =  dgrid(ig)%my%th0
      dth                    =  dgrid(ig)%my%dth
      wave_field             => dgrid(ig)%wave_field
      depth_field            => dgrid(ig)%depth_field
      ts_field               => dgrid(ig)%ts_field
      zz                     => dgrid(ig)%zz
      ubnd                   => dgrid(ig)%ubnd
      hbnd                   => dgrid(ig)%hbnd
      hzmax                  => dgrid(ig)%hzmax
#ifdef HZMINOUT
      hzmin                  => dgrid(ig)%hzmin
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
#ifndef SKIP_MAX_VEL
      vmax                   => dgrid(ig)%vmax
#endif
! ==============================================================================
      wod_flags              => dgrid(ig)%wod_flags
      wod_file_name          => dgrid(ig)%wod_file
      wod_field              => dgrid(ig)%wod_field
      bcf_file_name          => dgrid(ig)%bcf_file
      bcf_field              => dgrid(ig)%bcf_field
#ifdef MPI
      has_boundary          =  dgrid(ig)%my%has_boundary
! === Don't repeat allocate/deallocate! ========================================
      TIMER_START('allocate_edges')
      call allocate_edges(dgrid(ig))
      TIMER_STOP('allocate_edges')
! ==============================================================================
#endif

      !*** read in the bathymetry ***
      TIMER_START('read_bathymetry_gmt_grd')
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELIN
      call read_bathymetry_gmt_grd(file_name_bathymetry, depth_field, niz, njz, linear_flag)
#else
      call read_bathymetry_gmt_grd(file_name_bathymetry, depth_field, niz, njz, linear_flag, nxorg, nyorg)
#endif
#else
#ifndef PIXELIN
      call read_bathymetry_gmt_grd(file_name_bathymetry, depth_field, niz, njz, linear_flag, dgrid(ig), myrank)
#else
      call read_bathymetry_gmt_grd(file_name_bathymetry, depth_field, niz, njz, linear_flag, dgrid(ig), myrank, nxorg, nyorg)
#endif
#endif
      TIMER_STOP('read_bathymetry_gmt_grd')
#ifdef BANKFILE
      TIMER_START('read_bank_file')
      call read_bank_file(dgrid(ig))
      TIMER_STOP('read_bank_file')
#endif

      !*** read in the bottom friction coefficient ***
      TIMER_START('read_friction_gmt_grd')
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELIN
      call read_friction_gmt_grd(bcf_file_name, bcf_field, niz, njz)
#else
      call read_friction_gmt_grd(bcf_file_name, bcf_field, niz, njz, nxorg, nyorg)
#endif
#else
#ifndef PIXELIN
      call read_friction_gmt_grd(bcf_file_name, bcf_field, niz, njz, dgrid(ig), myrank)
#else
      call read_friction_gmt_grd(bcf_file_name, bcf_field, niz, njz, dgrid(ig), myrank, nxorg, nyorg)
#endif
#endif
      TIMER_STOP('read_friction_gmt_grd')

      !*** read seismic displacement ***
      ! thomas - to accommodate multiple ruptures, reading of the seismic displacements
      !          is left until the main loop now
      !          however here we intialise the wave field and the displacement field
      TIMER_START('initl_wfld')
      call initl_wfld(wave_field,zz,niz,njz)
      TIMER_STOP('initl_wfld')

#ifdef NCDIO
      TIMER_START('open_file')
#ifndef MULTI
#ifndef MPI
      call open_file(base, dgrid(ig))
#else
      call open_file(base, dgrid(ig), myrank)
#endif
#else
      temp_filename =  trim(members_dir) // trim(base)
#ifndef MPI
      call open_file(temp_filename, dgrid(ig))
#else
      call open_file(temp_filename, dgrid(ig), myrank)
#endif
#endif
      TIMER_STOP('open_file')
#endif

      !*** initialize wet or dry flags ***
      TIMER_START('wet_or_dry')
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELIN
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field)
#else
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field,nxorg,nyorg)
#endif
#else
#ifndef PIXELIN
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field,dgrid(ig),myrank)
#else
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field,dgrid(ig),myrank,nxorg,nyorg)
#endif
#endif
      TIMER_STOP('wet_or_dry')

      if(ig == 1) then ! set zero slope radiation condition for coarse grid only
         TIMER_START('boundary_rwg')
#ifndef MPI
#ifndef CARTESIAN
         call boundary_rwg(depth_field,hbnd,ubnd,dt,th0,dth,niz,njz)
#else
         call boundary_rwg(depth_field,hbnd,ubnd,dt,dxdy,niz,njz)
#endif
#else
#ifndef CARTESIAN
         call boundary_rwg(depth_field,hbnd,ubnd,dt,th0,dth,max(iy-2,0),niz,njz,has_boundary)
#else
         call boundary_rwg(depth_field,hbnd,ubnd,dt,dxdy,niz,njz,has_boundary)
#endif
#endif
         TIMER_STOP('boundary_rwg')
      end if

      TIMER_START('maxgrd_init_rwg')
      call maxgrd_init_rwg(hzmax,niz,njz)
      TIMER_STOP('maxgrd_init_rwg')
#ifdef HZMINOUT
      TIMER_START('mingrd_init_rwg')
      call mingrd_init_rwg(hzmin,niz,njz)
      TIMER_STOP('mingrd_init_rwg')
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
#ifndef SKIP_MAX_VEL
      TIMER_START('maxgrd_v_init_rwg')
      call maxgrd_v_init_rwg(vmax,niz,njz)
      TIMER_STOP('maxgrd_v_init_rwg')
#endif
! ==============================================================================
   end do

   ! find parents, initialize fields and set-up interpolation/copying mapping
   do ig = 2, ngrid
#ifdef MPI
      ierr = 0
#endif
      pid = 1
      do while(trim(dgrid(ig)%parent%base_name) /= trim(dgrid(pid)%my%base_name) .and. pid <= ngrid)
         pid = pid+1
      end do

      if(pid == ngrid+1) then
         write(0,'(a,i0,a,a,a)') '****** grid= ', ig, ' cannot find pid ', &
            trim(dgrid(ig)%parent%base_name), ', exiting...'
#ifndef MPI
         stop
#else
         ierr = ierr + 1
         call fatal_error(109)
#endif
      else
         dgrid(ig)%parent%id = pid
#ifndef MPI
      end if
#endif

      write(6,'(/,2x,a,i0,a,i0,a,a)') '*** grid ', ig, ', parent= ', dgrid(ig)%parent%id, &
         ' parent_base= ', trim(dgrid(ig)%parent%base_name)

      dgrid(ig)%parent%mlon0 = dgrid(pid)%my%mlon0
      dgrid(ig)%parent%mlat0 = dgrid(pid)%my%mlat0
      dgrid(ig)%parent%dh = dgrid(pid)%my%dh
      dgrid(ig)%parent%nx = dgrid(pid)%my%nx
      dgrid(ig)%parent%ny = dgrid(pid)%my%ny
      dgrid(ig)%my%nr = int(dgrid(ig)%parent%dh/dgrid(ig)%my%dh + 0.5d0)
#ifdef MPI
      ! MPI
      dgrid(ig)%parent%glon0 = dgrid(pid)%my%glon0
      dgrid(ig)%parent%glat0 = dgrid(pid)%my%glat0
      dgrid(ig)%parent%px = dgrid(pid)%my%px ! mpi proccess x
      dgrid(ig)%parent%py = dgrid(pid)%my%py ! mpi proccess y
      dgrid(ig)%parent%rx = dgrid(pid)%my%rx ! mpi rank x
      dgrid(ig)%parent%ry = dgrid(pid)%my%ry ! mpi rank y
      dgrid(ig)%parent%totalNx = dgrid(pid)%my%totalNx
      dgrid(ig)%parent%totalNy = dgrid(pid)%my%totalNy
      dgrid(ig)%parent%has_boundary = dgrid(pid)%my%has_boundary
      dgrid(ig)%parent%nbx    = dgrid(pid)%my%nbx
      dgrid(ig)%parent%nby    = dgrid(pid)%my%nby
      dgrid(ig)%parent%ix     = dgrid(pid)%my%ix
      dgrid(ig)%parent%iy     = dgrid(pid)%my%iy
      dgrid(ig)%parent%ixend  = dgrid(pid)%my%ixend
      dgrid(ig)%parent%iyend  = dgrid(pid)%my%iyend
! === DEBUG by tkato 2012/08/13 ================================================
      dgrid(ig)%parent%kx     = dgrid(pid)%my%kx
      dgrid(ig)%parent%ky     = dgrid(pid)%my%ky
      dgrid(ig)%parent%kxend  = dgrid(pid)%my%kxend
      dgrid(ig)%parent%kyend  = dgrid(pid)%my%kyend
! ==============================================================================
#endif

#ifndef CARTESIAN
      write(6,'(7x,a,i0,a,i0,a,f0.3,a,f0.3,a,i0)') 'niz=', dgrid(ig)%my%nx, ' njz=', dgrid(ig)%my%ny, &
         ' tend=', tend, '(sec) dt=', dt,'(sec) nstep=', nstep
      write(6,'(7x,a,f0.6,a,f0.6,a,e23.15e3)') 'th0=', dgrid(ig)%my%th0, '(rad) dxdy=', dgrid(ig)%my%dh, &
         '(deg) dth=',  dgrid(ig)%my%dth
      write(6,'(7x,a,f0.3,a,f0.3,a,f0.3,a,f0.3,a)') 'mlon0=', dgrid(ig)%my%mlon0, &
         ' (', dgrid(ig)%my%mlon0/60.0d0, ') mlat0=', dgrid(ig)%my%mlat0, ' (', dgrid(ig)%my%mlat0/60.0d0, ')'

#ifndef MPI
      mlon0 = (int((dgrid(ig)%my%mlon0-dgrid(ig)%parent%mlon0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0))*dgrid(ig)%parent%dh
      mlon0 = mlon0 + dgrid(ig)%parent%mlon0/60.0d0
      mlat0 = (int((dgrid(ig)%my%mlat0-dgrid(ig)%parent%mlat0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0))*dgrid(ig)%parent%dh
      mlat0 = mlat0 + dgrid(ig)%parent%mlat0/60.0d0
#else
      mlon0 = (int((dgrid(ig)%my%glon0-dgrid(ig)%parent%glon0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0))*dgrid(ig)%parent%dh
      mlon0 = mlon0 + dgrid(ig)%parent%glon0/60.0d0
      mlat0 = (int((dgrid(ig)%my%glat0-dgrid(ig)%parent%glat0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0))*dgrid(ig)%parent%dh
      mlat0 = mlat0 + dgrid(ig)%parent%glat0/60.0d0
#endif
#else
      write(6,'(7x,a,i0,a,i0,a,f0.3,a,f0.3,a,i0)') 'niz=', dgrid(ig)%my%nx, ' njz=', dgrid(ig)%my%ny, &
         ' tend=', tend, '(sec) dt=', dt,'(sec) nstep=', nstep
      write(6,'(7x,a,f0.6,a,f0.6,a,e23.15e3)') 'th0=', dgrid(ig)%my%th0, '(m) dxdy=', dgrid(ig)%my%dh, &
         '(m) dth=',  dgrid(ig)%my%dth
      write(6,'(7x,a,f0.3,a,f0.3)') 'mlon0=', dgrid(ig)%my%mlon0, &
         ' mlat0=', dgrid(ig)%my%mlat0

#ifndef MPI
      mlon0 = (int((dgrid(ig)%my%mlon0-dgrid(ig)%parent%mlon0)/dgrid(ig)%parent%dh + 0.5d0))*dgrid(ig)%parent%dh
      mlon0 = mlon0 + dgrid(ig)%parent%mlon0
      mlat0 = (int((dgrid(ig)%my%mlat0-dgrid(ig)%parent%mlat0)/dgrid(ig)%parent%dh + 0.5d0))*dgrid(ig)%parent%dh
      mlat0 = mlat0 + dgrid(ig)%parent%mlat0
#else
      mlon0 = (int((dgrid(ig)%my%glon0-dgrid(ig)%parent%glon0)/dgrid(ig)%parent%dh + 0.5d0))*dgrid(ig)%parent%dh
      mlon0 = mlon0 + dgrid(ig)%parent%glon0
      mlat0 = (int((dgrid(ig)%my%glat0-dgrid(ig)%parent%glat0)/dgrid(ig)%parent%dh + 0.5d0))*dgrid(ig)%parent%dh
      mlat0 = mlat0 + dgrid(ig)%parent%glat0
#endif
#endif

#ifndef CARTESIAN
      write(6,'(/,2x,a)') '*** from global origin (deg):'
#ifndef MPI
      write(6,'(7x,a,f13.8,a,f13.8)') 'mlon0= ', mlon0, ' input value= ', dgrid(ig)%my%mlon0/60.0d0
      write(6,'(7x,a,f13.8,a,f13.8)') 'mlat0= ', mlat0, ' input value= ', dgrid(ig)%my%mlat0/60.0d0
#else
      write(6,'(7x,a,f13.8,a,f13.8)') 'mlon0= ', mlon0, ' input value= ', dgrid(ig)%my%glon0/60.0d0
      write(6,'(7x,a,f13.8,a,f13.8)') 'mlat0= ', mlat0, ' input value= ', dgrid(ig)%my%glat0/60.0d0
#endif
#else
      write(6,'(/,2x,a)') '*** from global origin (meter):'
#ifndef MPI
      write(6,'(7x,a,f0.8,a,f0.8)') 'mlon0= ', mlon0, ' input value= ', dgrid(ig)%my%mlon0
      write(6,'(7x,a,f0.8,a,f0.8)') 'mlat0= ', mlat0, ' input value= ', dgrid(ig)%my%mlat0
#else
      write(6,'(7x,a,f0.8,a,f0.8)') 'mlon0= ', mlon0, ' input value= ', dgrid(ig)%my%glon0
      write(6,'(7x,a,f0.8,a,f0.8)') 'mlat0= ', mlat0, ' input value= ', dgrid(ig)%my%glat0
#endif
#endif

#ifndef CARTESIAN
#ifndef MPI
      if((mlon0 <= dgrid(ig)%my%mlon0/60.0d0 - 0.0001d0) .or. (mlon0 >= dgrid(ig)%my%mlon0/60.0d0 + 0.0001d0) .or. &
         (mlat0 <= dgrid(ig)%my%mlat0/60.0d0 - 0.0001d0) .or. (mlat0 >= dgrid(ig)%my%mlat0/60.0d0 + 0.0001d0)) then
#else
      if((mlon0 <= dgrid(ig)%my%glon0/60.0d0 - 0.0001d0) .or. (mlon0 >= dgrid(ig)%my%glon0/60.0d0 + 0.0001d0) .or. &
         (mlat0 <= dgrid(ig)%my%glat0/60.0d0 - 0.0001d0) .or. (mlat0 >= dgrid(ig)%my%glat0/60.0d0 + 0.0001d0)) then
#endif
         write(0,'(a,i0,a)') '****** grid= ', ig, ' not properly registered to coarse grid, exiting...'
#ifndef MPI
         stop
#else
         write(0,'(a)') '  *** from global origin (deg):\n'
         write(0,'(7x,a,f13.8,a,f13.8)') 'mlon0= ', mlon0, ' input value= ', dgrid(ig)%my%glon0/60.0d0
         write(0,'(7x,a,f13.8,a,f13.8)') 'mlat0= ', mlat0, ' input value= ', dgrid(ig)%my%glat0/60.0d0
         ierr = ierr + 1
         call fatal_error(110)
#endif
#else
#ifndef MPI
      if((mlon0 <= dgrid(ig)%my%mlon0 - 0.0001d0) .or. (mlon0 >= dgrid(ig)%my%mlon0 + 0.0001d0) .or. &
         (mlat0 <= dgrid(ig)%my%mlat0 - 0.0001d0) .or. (mlat0 >= dgrid(ig)%my%mlat0 + 0.0001d0)) then
#else
      if((mlon0 <= dgrid(ig)%my%glon0 - 0.0001d0) .or. (mlon0 >= dgrid(ig)%my%glon0 + 0.0001d0) .or. &
         (mlat0 <= dgrid(ig)%my%glat0 - 0.0001d0) .or. (mlat0 >= dgrid(ig)%my%glat0 + 0.0001d0)) then
#endif
         write(0,'(a,i0,a)') '****** grid= ', ig, ' not properly registered to coarse grid, exiting...'
#ifndef MPI
         stop
#else
         write(0,'(a)') '  *** from global origin (meter):\n'
         write(0,'(7x,a,f13.8,a,f13.8)') 'mlon0= ', mlon0, ' input value= ', dgrid(ig)%my%glon0
         write(0,'(7x,a,f13.8,a,f13.8)') 'mlat0= ', mlat0, ' input value= ', dgrid(ig)%my%glat0
         ierr = ierr + 1
         call fatal_error(110)
#endif
#endif
      end if

      write(6,'(/,13x,a,f0.6)') 'H= ', dgrid(ig)%parent%dh
      write(6,'(13x,a,f0.6)') 'h= ', dgrid(ig)%my%dh
      write(6,'(13x,a,f0.0)') 'H/h= ', dgrid(ig)%parent%dh/dgrid(ig)%my%dh
      write(6,'(7x,a,i0,/)') 'grid reduction factor= ', dgrid(ig)%my%nr

#ifndef MPI
      if((mod(dgrid(ig)%my%nx-1,dgrid(ig)%my%nr) /= 0) .or. (mod(dgrid(ig)%my%ny-1,dgrid(ig)%my%nr) /= 0)) then
#else
      if((mod(dgrid(ig)%my%totalNx-1,dgrid(ig)%my%nr) /= 0) .or. (mod(dgrid(ig)%my%totalNy-1,dgrid(ig)%my%nr) /= 0)) then
#endif
         write(0,'(a,i0,a)') '****** grid= ', ig, ' not properly registered to coarse grid, exiting...'
#ifndef MPI
         write(0,'(a,i0,a,i0,a)') '****** (nx - 1 = ', dgrid(ig)%my%nx-1, ' ) / (nr = ', dgrid(ig)%my%nr, ' ) must be zero'
         write(0,'(a,i0,a,i0,a)') '****** (ny - 1 = ', dgrid(ig)%my%ny-1, ' ) / (nr = ', dgrid(ig)%my%nr, ' ) must be zero'
         stop
#else
         write(0,'(a,i0,a,i0,a)') '****** (nx - 1 = ', dgrid(ig)%my%totalNx-1, ' ) / (nr = ', dgrid(ig)%my%nr, ' ) must be zero'
         write(0,'(a,i0,a,i0,a)') '****** (ny - 1 = ', dgrid(ig)%my%totalNy-1, ' ) / (nr = ', dgrid(ig)%my%nr, ' ) must be zero'
         ierr = ierr + 1
         call fatal_error(111)
#endif
      end if

#ifndef CARTESIAN
#ifndef MPI
      dgrid(ig)%my%bigNX = ((dgrid(ig)%my%nx-1)/dgrid(ig)%my%nr + 1)
      dgrid(ig)%my%bigNY = ((dgrid(ig)%my%ny-1)/dgrid(ig)%my%nr + 1)

      dgrid(ig)%my%zeroIX = int((dgrid(ig)%my%mlon0 - dgrid(ig)%parent%mlon0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0) + 1
      dgrid(ig)%my%zeroIY = int((dgrid(ig)%my%mlat0 - dgrid(ig)%parent%mlat0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0) + 1

      write(6,'(13x,a,i0,a,i0,a,i0)') 'nx= ', dgrid(ig)%my%nx, ' zeroIX= ', dgrid(ig)%my%zeroIX, ' NX= ', dgrid(ig)%my%bigNX
      write(6,'(13x,a,i0,a,i0,a,i0,/)') 'ny= ', dgrid(ig)%my%ny, ' zeroIY= ', dgrid(ig)%my%zeroIY, ' NY= ', dgrid(ig)%my%bigNY
#else
      dgrid(ig)%my%bigNX = ((dgrid(ig)%my%totalNx-1)/dgrid(ig)%my%nr + 1)
      dgrid(ig)%my%bigNY = ((dgrid(ig)%my%totalNy-1)/dgrid(ig)%my%nr + 1)
      dgrid(ig)%my%bigIX = ((dgrid(ig)%my%ix       )/dgrid(ig)%my%nr    )
      dgrid(ig)%my%bigIY = ((dgrid(ig)%my%iy       )/dgrid(ig)%my%nr    )

      dgrid(ig)%my%zeroIX = int((dgrid(ig)%my%glon0 - dgrid(ig)%parent%glon0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0) + 1
      dgrid(ig)%my%zeroIY = int((dgrid(ig)%my%glat0 - dgrid(ig)%parent%glat0)/(60.0d0*dgrid(ig)%parent%dh) + 0.5d0) + 1

      write(6,'(13x,a,i0,a,i0,a,i0)')   'totalNx= ', dgrid(ig)%my%totalNx, &
         ' zeroIX= ', dgrid(ig)%my%zeroIX, ' NX= ', dgrid(ig)%my%bigNX
      write(6,'(13x,a,i0,a,i0,a,i0,/)') 'totalNy= ', dgrid(ig)%my%totalNy, &
         ' zeroIY= ', dgrid(ig)%my%zeroIY, ' NY= ', dgrid(ig)%my%bigNY
#endif
#else
#ifndef MPI
      dgrid(ig)%my%bigNX = ((dgrid(ig)%my%nx-1)/dgrid(ig)%my%nr + 1)
      dgrid(ig)%my%bigNY = ((dgrid(ig)%my%ny-1)/dgrid(ig)%my%nr + 1)

      dgrid(ig)%my%zeroIX = int((dgrid(ig)%my%mlon0 - dgrid(ig)%parent%mlon0)/dgrid(ig)%parent%dh + 0.5d0) + 1
      plat0 = dgrid(ig)%parent%mlat0 + (dgrid(ig)%parent%ny-1)*dgrid(ig)%parent%dh
      clat0 = dgrid(ig)%my%mlat0     + (dgrid(ig)%my%ny-1)*    dgrid(ig)%my%dh
      dgrid(ig)%my%zeroIY = int((plat0 - clat0)/dgrid(ig)%parent%dh + 0.5d0) + 1

      write(6,'(13x,a,i0,a,i0,a,i0)') 'nx= ', dgrid(ig)%my%nx, ' zeroIX= ', dgrid(ig)%my%zeroIX, ' NX= ', dgrid(ig)%my%bigNX
      write(6,'(13x,a,i0,a,i0,a,i0,/)') 'ny= ', dgrid(ig)%my%ny, ' zeroIY= ', dgrid(ig)%my%zeroIY, ' NY= ', dgrid(ig)%my%bigNY
#else
      dgrid(ig)%my%bigNX = ((dgrid(ig)%my%totalNx-1)/dgrid(ig)%my%nr + 1)
      dgrid(ig)%my%bigNY = ((dgrid(ig)%my%totalNy-1)/dgrid(ig)%my%nr + 1)
      dgrid(ig)%my%bigIX = ((dgrid(ig)%my%ix       )/dgrid(ig)%my%nr    )
      dgrid(ig)%my%bigIY = ((dgrid(ig)%my%iy       )/dgrid(ig)%my%nr    )

      dgrid(ig)%my%zeroIX = int((dgrid(ig)%my%glon0 - dgrid(ig)%parent%glon0)/dgrid(ig)%parent%dh + 0.5d0) + 1
      plat0 = dgrid(ig)%parent%glat0 + (dgrid(ig)%parent%totalNy-1)*dgrid(ig)%parent%dh
      clat0 = dgrid(ig)%my%glat0     + (dgrid(ig)%my%totalNy-1)*    dgrid(ig)%my%dh
      dgrid(ig)%my%zeroIY = int((plat0 - clat0)/dgrid(ig)%parent%dh + 0.5d0) + 1

      write(6,'(13x,a,i0,a,i0,a,i0)')   'totalNx= ', dgrid(ig)%my%totalNx, &
         ' zeroIX= ', dgrid(ig)%my%zeroIX, ' NX= ', dgrid(ig)%my%bigNX
      write(6,'(13x,a,i0,a,i0,a,i0,/)') 'totalNy= ', dgrid(ig)%my%totalNy, &
         ' zeroIY= ', dgrid(ig)%my%zeroIY, ' NY= ', dgrid(ig)%my%bigNY
#endif
#endif

      TIMER_START('initl_gridmap')
#ifndef USE_ALLTOALLV
! === Config of copy from child to perent. by tkato 2012/11/15 =================
!     call initl_gridmap(dgrid(ig))
      call initl_gridmap(dgrid(ig), c2p_all)
! ==============================================================================
#else
! === USE_MPI_ALLTOALLV ========================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
!     call initl_gridmap(dgrid(ig), nprocs)
      call initl_gridmap(dgrid(ig), nprocs, c2p_all)
! ==============================================================================
! === USE_MPI_ALLTOALLV ========================================================
#endif
      TIMER_STOP('initl_gridmap')
      TIMER_START('initl_gridmap_dz')
      call initl_gridmap_dz(dgrid(ig))
      TIMER_STOP('initl_gridmap_dz')
      if(ig /= 1) then
         pid = dgrid(ig)%parent%id
         TIMER_START('interp2fine_dz')
         call interp2fine_dz(dgrid(pid),dgrid(ig))
         TIMER_STOP('interp2fine_dz')
      end if
#ifdef MPI
      TIMER_START('exchange_edges_dz')
      call exchange_edges_dz(dgrid(ig))
      TIMER_STOP('exchange_edges_dz')
#endif

      if(smooth_edges == 1) then
         TIMER_START('smooth_bath')
         call smooth_bath(dgrid(pid),dgrid(ig))
         TIMER_STOP('smooth_bath')
         if(ig /= 1) then
            pid = dgrid(ig)%parent%id
            TIMER_START('interp2fine_dz')
            call interp2fine_dz(dgrid(pid),dgrid(ig))
            TIMER_STOP('interp2fine_dz')
         end if
#ifdef MPI
         TIMER_START('exchange_edges_dz')
         call exchange_edges_dz(dgrid(ig))
         TIMER_STOP('exchange_edges_dz')
#endif
         ! recheck wet or dry since depth may have changed
         TIMER_START('wet_or_dry')
         call wet_or_dry(dgrid(ig)%wave_field,dgrid(ig)%depth_field,dgrid(ig)%wod_flags, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELIN
                         dgrid(ig)%my%nx,dgrid(ig)%my%ny,dgrid(ig)%wod_file,dgrid(ig)%wod_field)
#else
                         dgrid(ig)%my%nx,dgrid(ig)%my%ny,dgrid(ig)%wod_file,dgrid(ig)%wod_field, &
                         dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
                         dgrid(ig)%my%nx,dgrid(ig)%my%ny,dgrid(ig)%wod_file,dgrid(ig)%wod_field, &
#ifndef PIXELIN
                         dgrid(ig),myrank)
#else
                         dgrid(ig),myrank,dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
         TIMER_STOP('wet_or_dry')
      end if

#ifdef MPI
      end if
      if(ierr > 0) call fatal_error(109)
#endif
   end do

   !********************************************
   !*** load in the locations of tide gauges ***
   !********************************************

   ! thomas - now tide gauge data only dumped every itgrn steps, so should send dt*itgrn instead of dt
   !          and nstep/itgrn instead of nstep
   TIMER_START('tgs_open_ct')
! === Split Dir ================================================================
#ifdef MULTI
   tg_station_file_name = trim(input_dirname) // trim(tg_station_file_name)
#endif
! ==============================================================================
#ifndef __SX__
   call tgs_open_ct(dgrid,ngrid,dt*REAL_FUNC(itgrn),nstep/itgrn,nsta,program_name,tg_station_file_name,mytgs)
#else
   !*** read in the latitude and longitude ****
   open(1,file=trim(tg_station_file_name),action='read',status='old',form='formatted',err=500)
   ista = 0

   read(1,'(a)',err=501) str
   read(str,'(i)',err=501) nsta
   if(nsta <= 0) goto 501

   allocate(mytgs(nsta))

#ifndef MPI
   do ista = 1, nsta
#else
   ista = 1
   do ksta = 1, nsta
#endif
      read(1,'(a)') str
      read(str,*) mytgs(ista)%geolat, mytgs(ista)%geolon, mytgs(ista)%number

#ifndef CARTESIAN
      colat = 90.0d0 - mytgs(ista)%geolat
      if(mytgs(ista)%geolon < 0.0d0) then
         colon = 360.0d0 + mytgs(ista)%geolon
      else
         colon = mytgs(ista)%geolon
      end if

      mytgs(ista)%mcolat = colat*60.0d0
      mytgs(ista)%mcolon = colon*60.0d0
#else
      mytgs(ista)%mcolat = mytgs(ista)%geolat
      mytgs(ista)%mcolon = mytgs(ista)%geolon
#endif
      mytgs(ista)%ig = 0
      mytgs(ista)%ilon = 0
      mytgs(ista)%ilat = 0

      mindh = 2.0d0*dgrid(1)%my%dh
#ifdef MPI
      kg = -1
#endif
      do ig = 1, ngrid
         call tgs_find_grid_coords(mytgs(ista)%mcolon,mytgs(ista)%mcolat, &
#ifndef MPI
#ifndef CARTESIAN
                                   dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,REAL_FUNC(dgrid(ig)%my%dh*60.0d0), &
#else
                                   dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,REAL_FUNC(dgrid(ig)%my%dh), &
#endif
#else
#ifndef CARTESIAN
                                   dgrid(ig)%my%glon0,dgrid(ig)%my%glat0,REAL_FUNC(dgrid(ig)%my%dh*60.0d0), &
#else
                                   dgrid(ig)%my%glon0,dgrid(ig)%my%glat0,REAL_FUNC(dgrid(ig)%my%dh), &
#endif
#endif
                                   ilon,ilat)
! === DEBUG by tkato 2015/03/04 ================================================
#ifdef CARTESIAN
#ifndef MPI
         ilat = dgrid(ig)%my%ny - ilat + 1
#else
         ilat = dgrid(ig)%my%totalNy - ilat + 1
#endif
#endif
! ==============================================================================
#ifndef MPI
         if(ilon < 1 .or. ilon > dgrid(ig)%my%nx .or. ilat < 1 .or. ilat > dgrid(ig)%my%ny) then
            !  tide gauge is not in extent of dg[ig]
            if(ig == 1) then
                ! we assume that ig=1 is the ancestor of all grids and if the tide gauge is
                ! outside of that then exit
                write(0,'(a,i0,a)') 'Tide Gauge Point ', ista, ' appears to be outside the box.'
                ! DB added some more information here
                write(0,'(a,f,a,f)') 'Station geolon=', mytgs(ista)%geolon, ' geolat=', mytgs(ista)%geolat
                write(0,'(a,i0,a,i0,a,i0,a,i0)') 'ilon=', ilon, ' ilat=', ilat, ' nx=', dgrid(ig)%my%nx, ' ny=', dgrid(ig)%my%ny
                stop
            else
               ! in this case ig is not the big grid and tide gauge is allowed to be
               !  outside, so  continue on to the next iteration of the loop
               cycle
            end if
         end if

         if(dgrid(ig)%my%dh < mindh) then
            mindh = dgrid(ig)%my%dh
            mytgs(ista)%ig = ig

            mytgs(ista)%ilon = ilon
            mytgs(ista)%ilat = ilat
            mytgs(ista)%dt = dt
            mytgs(ista)%nt = nstep

            mytgs(ista)%z = dgrid(ig)%depth_field%dz(mytgs(ista)%ilon,mytgs(ista)%ilat)
         end if
      end do
#else
         if(ilon >= 1 .and. ilon <= dgrid(ig)%my%totalNx .and. ilat >= 1 .and. ilat <= dgrid(ig)%my%totalNy) then
            if(dgrid(ig)%my%dh < mindh) then
               mindh = dgrid(ig)%my%dh
               kg = ig
            end if
         end if
      end do
      if(kg < 0) then
         write(0,'(a,i0,a)') 'Tide Gauge Point ', ksta, ' appears to be outside the box.'
         ! DB added some more information here
         write(0,'(a,f,a,f)') 'Station geolon=', mytgs(ista)%geolon, ' geolat=', mytgs(ista)%geolat
         write(0,'(a,i0,a,i0,a,i0,a,i0)') 'ilon=', ilon, ' ilat=', ilat, &
            ' nx=', dgrid(1)%my%totalNx, ' ny=', dgrid(1)%my%totalNy
         call fatal_error(303)
      end if

      ! set if tgs exists in local area
      call tgs_find_grid_coords(mytgs(ista)%mcolon,mytgs(ista)%mcolat, &
#ifndef CARTESIAN
                                dgrid(kg)%my%glon0,dgrid(kg)%my%glat0,REAL_FUNC(dgrid(kg)%my%dh*60.0d0), &
#else
                                dgrid(kg)%my%glon0,dgrid(kg)%my%glat0,REAL_FUNC(dgrid(kg)%my%dh), &
#endif
                                ilon,ilat)
! === DEBUG by tkato 2015/03/04 ================================================
#ifdef CARTESIAN
      ilat = dgrid(kg)%my%totalNy - ilat + 1
#endif
! ==============================================================================
      if(ilon >= dgrid(kg)%my%ix .and. ilon <= dgrid(kg)%my%ixend .and. &
         ilat >= dgrid(kg)%my%iy .and. ilat <= dgrid(kg)%my%iyend) then
         mytgs(ista)%ig = kg

         mytgs(ista)%ilon = ilon - dgrid(kg)%my%kx + 1
         mytgs(ista)%ilat = ilat - dgrid(kg)%my%ky + 1
         mytgs(ista)%dt = dt
         mytgs(ista)%nt = nstep

         mytgs(ista)%z = dgrid(kg)%depth_field%dz(mytgs(ista)%ilon,mytgs(ista)%ilat)

         ista = ista + 1
      end if
#endif
   end do
   close(1)
#ifdef MPI
   nsta = ista - 1
#endif

   write(6,'(a,a,i0,a,a)') trim(program_name), ': tgs_open_rwg.c: nsta=', nsta, &
      ' read from file ', trim(tg_station_file_name)

503 continue
#endif
   TIMER_STOP('tgs_open_ct')

   !*** thomas - all this conditional on nsta > 0 ***
   if(nsta > 0) then
      ! Burbidge: Changed this to stdout so that redirect output works OK
      write(6,'(a,a,i0,a)') trim(program_name), ': nsta=', nsta, ' read sucessfully'
      ! tgs text output
      if(tgstxtout > 0) then
#ifndef SINGLE_TGS
         allocate(tgsfiles(nsta))
         allocate(tgsfp(nsta))
         do j = 1, nsta
            tgsfiles(j) = ''
            tgsfp(j) = 200 + j
         end do
         do j = 1, nsta
#ifndef MULTI
            write(tgsfiles(j),'(a,i6.6)') trim(tgstxtoutfile), mytgs(j)%number
#else
            write(tgsfiles(j),'(a,i6.6)') &
               trim(members_dir) // trim(tgstxtoutfile), mytgs(j)%number
#endif
            open(tgsfp(j),file=trim(tgsfiles(j)),action='write',status='replace',form='formatted',err=200)
            cycle
200         write(0,'(a,a,a)') 'tgs text out : file open error "', trim(tgsfiles(j)), '"'
         end do
#else
         tgsfiles = ''
         tgsfp = 200
#ifndef MPI
#ifndef MULTI
         write(tgsfiles,'(a)') trim(tgstxtoutfile)
#else
         write(tgsfiles,'(a)') trim(members_dir) // trim(tgstxtoutfile)
#endif
#else
#ifndef MULTI
         write(tgsfiles,'(a)') trim(tgstxtoutfile) // trim(suffix)
#else
         write(tgsfiles,'(a)') trim(members_dir) // trim(tgstxtoutfile) // trim(suffix)
#endif
#endif
         open(tgsfp,file=trim(tgsfiles),action='write',status='replace',form='formatted')
#endif
      end if
   else
      write(6,'(a,a)') trim(program_name), ': warning, no tide gauges read.'
   end if

   ! thomas - counter to keep track of progressive rupture in main loop
! === Support multiple ruptures. ===============================================
   do ig = 1, ngrid
      jrupt => dgrid(ig)%jrupt
      jrupt = -1
   end do
! ==============================================================================

#ifdef CONV_CHECK
#ifdef MPI
   if(myrank == 0) then
#endif
      if(with_disp == 1 .or. with_disp == 2) then
         do ig = 1, ngrid
            if(with_disp == 1 .or. ig /= 1) then
               dgrid(ig)%my%nconvout = 1000 + ig
               str = 'conv_step.' // trim(dgrid(ig)%my%base_name)
#ifdef MULTI
               str = trim(members_dir) // trim(str)
#endif
               open(dgrid(ig)%my%nconvout,file=str,action='write',status='replace',form='formatted')
            end if
         end do
      end if
#ifdef MPI
   end if
#endif
#endif
#ifndef CARTESIAN
! === Density ==================================================================
   if(with_density == 1) then
      do ig = 1, ngrid
         niz =  dgrid(ig)%my%nx
         njz =  dgrid(ig)%my%ny
         TIMER_START('density_initialize')
         call density_initialize(dgrid(ig),niz,njz)
         TIMER_STOP('density_initialize')
      end do
   end if
! ==============================================================================
! === Elastic Loading ==========================================================
   if(with_elastic_loading == 1) then
#ifdef MPI
      TIMER_START('loading_mpi_initialize')
      call loading_mpi_initialize(nprocs, myrank, npx, npy, rankx, ranky)
      TIMER_STOP('loading_mpi_initialize')
#endif
      TIMER_START('loading_getval')
#ifndef MULTI
      call loading_getval()
#else
      call loading_getval(input_dirname)
#endif
      TIMER_STOP('loading_getval')
      do ig = 1, ngrid
         TIMER_START('loading_initialize')
         call loading_initialize(dgrid(ig))
         TIMER_STOP('loading_initialize')
      end do
   end if
! ==============================================================================
#endif
#ifdef BANKFILE
   do ig = 1, ngrid
      TIMER_START('update_bathymetory')
      call update_bathymetory(dgrid(ig))
      TIMER_STOP('update_bathymetory')
   end do
#endif
#ifdef DUMP1D
   TIMER_START('dump1d_initialize')
   call dump1d_initialize(dgrid)
   TIMER_STOP('dump1d_initialize')
#endif
   !*** main loop ***
! === Support restart ==========================================================
!  do istep = 1, nstep
   if(restart == 0) then
      istart = 1
   else
      istart = restart + 1
#ifndef MPI
      write(restart_file_name, '(a,i8.8)') 'restart.', restart
#else
      write(restart_file_name, '(a,i8.8,a,i6.6)') 'restart.', restart, '.', myrank
#endif
      write(6,'(a,i0)') '[RESTART] Restart from step ', restart
      write(6,'(a,a)') '[RESTART] Restart file read: ', trim(restart_file_name)
      write(6,'(a,i0)') '[RESTART] Next step is ', istart
#ifdef MULTI
      restart_file_name = trim(members_dir) // trim(restart_file_name)
#endif
      TIMER_START('read_restart_file')
      call read_restart_file(ngrid, dgrid, restart_file_name)
      TIMER_STOP('read_restart_file')
   end if
   do istep = istart, nstep
! ==============================================================================
      t = REAL_FUNC(istep) * dt
! === Support truncation =======================================================
      if(max_time_i /= 0) then
#ifndef MPI
         trunc_flag = check_trunc(max_time_i)
#else
         if(myrank == 0) trunc_flag = check_trunc(max_time_i)
#ifndef MULTI
         call MPI_Bcast(trunc_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#else
         call MPI_Bcast(trunc_flag, 1, MPI_LOGICAL, 0, MPI_MEMBER_WORLD, ierr)
#endif
#endif
         if(trunc_flag) then
            write(6,'(a)')          '======================================================='
            write(6,'(a)')          '======================================================='
            write(6,'(a)')          '=== NOTE!!!'
            write(6,'(a,i8,a)')     '=== Execution is truncated at step ', istep-1, '!'
            write(6,'(a,i8,a)')     '=== Because elapsed time reached specified max_time.'
            write(6,'(a,a,a,i0,a)') '=== Max time: ', trim(max_time), ' (', max_time_i, ' seconds)'
            write(6,'(a)')          '======================================================='
            write(6,'(a)')          '======================================================='
            exit
         end if
      end if
! ==============================================================================

      ! thomas - modified to accommodate multiple rupture steps
! === SINWAVE ==================================================================
      if(init_disp_sinwave == 0) then
! ==============================================================================
! === Support multiple ruptures. ===============================================
      do ig = 1, ngrid
         nrupt   => dgrid(ig)%nrupt
         irupt   => dgrid(ig)%irupt
         jrupt   => dgrid(ig)%jrupt
         ruptgrd => dgrid(ig)%ruptgrd
         if(t <= nrupt*tau) then
            irupt = -1
            do while((irupt+1)*tau < t)
               irupt = irupt + 1
            end do
            irupt = irupt + 1

            if(irupt > jrupt) then
! === For MRI ==================================================================
               if(init_disp_gaussian == 0) then
! ==============================================================================
! === Initial displacement of child domains is given by interpolation. =========
               if((init_disp_interpolation /= 1) .or. (ig == 1))  then
! ==============================================================================
! === Displacement =============================================================
               if(init_disp_fault == 0) then
! ==============================================================================
               write(6,'(a,i0,a,f0.6,a,f0.6,a,f0.6,a,i0,a,i0,a,i0,a,a)') &
                  'New rupture step: istep=', istep, ' dt=', dt, ' t=', t, ' tau=', tau, ' ig=', ig, &
                  ' irupt=', irupt, ' nrupt=', nrupt, ' file=', trim(ruptgrd(irupt))
               TIMER_START('read_rupture')
               call read_rupture(dgrid(ig)%zz,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELIN
                                 dgrid(ig)%my%mlat0,dgrid(ig)%my%mlon0,dgrid(ig)%my%dh, &
                                 ruptgrd(irupt),program_name)
#else
                                 dgrid(ig)%my%dh, &
                                 ruptgrd(irupt),program_name,dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELIN
                                 dgrid(ig)%my%mlat0,dgrid(ig)%my%mlon0,dgrid(ig)%my%dh, &
                                 ruptgrd(irupt),program_name,dgrid(ig),myrank)
#else
                                 dgrid(ig)%my%dh, &
                                 ruptgrd(irupt),program_name,dgrid(ig),myrank,dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
               TIMER_STOP('read_rupture')
! === Displacement =============================================================
               else
                  if(ig == 1) then
                     if(multrupt == 1) fault_param_file = ruptgrd(irupt)
                     TIMER_START('displacement_initialize')
#ifndef MULTI
                     call displacement_initialize()
#else
                     call displacement_initialize(input_dirname)
#endif
                     TIMER_STOP('displacement_initialize')
                  end if
                  TIMER_START('displacement_calc_displacement')
                  call displacement_calc_displacement(dgrid(ig), ig)
                  TIMER_STOP('displacement_calc_displacement')
               end if
! ==============================================================================
! === Displacement =============================================================
! === DEBUG: Kajiura filter is applied 2 times when "init_disp_fault == 0". ====
!              if(apply_kj_filter == 1) then
               if((apply_kj_filter == 1) .and. (init_disp_fault == 0)) then
! ==============================================================================
#ifndef REAL_DBLE
                  dgrid(ig)%zz_dp = dgrid(ig)%zz
#endif
#ifndef CARTESIAN
#ifndef REAL_DBLE
                  if(ig == 1) then
                     TIMER_START('displacement_calc_h0_lat1')
                     call displacement_calc_h0_lat1(dgrid(ig), dgrid(ig)%zz_dp, &
#else
                  if(ig == 1) then
                     TIMER_START('displacement_calc_h0_lat1')
                     call displacement_calc_h0_lat1(dgrid(ig), dgrid(ig)%zz, &
#endif
                     dgrid(ig)%my%nx, dgrid(ig)%my%ny, h0, lat1)
                     TIMER_STOP('displacement_calc_h0_lat1')
                  end if
                  TIMER_START('displacement_apply_kj_filter')
#ifndef REAL_DBLE
                  call displacement_apply_kj_filter(dgrid(ig), dgrid(ig)%zz_dp, &
#else
                  call displacement_apply_kj_filter(dgrid(ig), dgrid(ig)%zz, &
#endif
                     dgrid(ig)%my%nx, dgrid(ig)%my%ny, h0, lat1)
                  TIMER_STOP('displacement_apply_kj_filter')
#else
#ifndef REAL_DBLE
                  if(ig == 1) then
                     TIMER_START('displacement_calc_h0')
                     call displacement_calc_h0(dgrid(ig), dgrid(ig)%zz_dp, &
#else
                  if(ig == 1) then
                     TIMER_START('displacement_calc_h0')
                     call displacement_calc_h0(dgrid(ig), dgrid(ig)%zz, &
#endif
                     dgrid(ig)%my%nx, dgrid(ig)%my%ny, h0)
                     TIMER_STOP('displacement_calc_h0')
                  end if
                  TIMER_START('displacement_apply_kj_filter')
#ifndef REAL_DBLE
                  call displacement_apply_kj_filter(dgrid(ig), dgrid(ig)%zz_dp, &
#else
                  call displacement_apply_kj_filter(dgrid(ig), dgrid(ig)%zz, &
#endif
                     dgrid(ig)%my%nx, dgrid(ig)%my%ny, h0)
                  TIMER_STOP('displacement_apply_kj_filter')
#endif
#ifdef MPI
                  call exchange_edges_zz(dgrid(ig))
#endif
               end if
! ==============================================================================
! === Initial displacement of child domains is given by interpolation. =========
               else
                  pid = dgrid(ig)%parent%id
                  TIMER_START('interp2fine_init_disp')
                  call interp2fine_init_disp(dgrid(pid),dgrid(ig))
                  TIMER_STOP('interp2fine_init_disp')
#ifdef MPI
                  call exchange_edges_zz(dgrid(ig))
#endif
               end if
! ==============================================================================
! === Displacement =============================================================
               if((init_disp_fault == 1) .and. (ig == ngrid)) call displacement_finalize()
! ==============================================================================
! === For MRI ==================================================================
               else
                  write(6,'(a,i0,a,f0.6,a,f0.6,a,f0.6,a,i0,a,i0,a,i0,a,a)') &
                     'New rupture step: istep=', istep, ' dt=', dt, ' t=', t, ' tau=', tau, ' ig=', ig, &
                     ' irupt=', irupt, ' nrupt=', nrupt, ' file=', trim(ruptgrd(irupt))
                  TIMER_START('make_gaussian_rupture')
                  call make_gaussian_rupture(dgrid(ig))
                  TIMER_STOP('make_gaussian_rupture')
               end if
! ==============================================================================
! === Multiple rupture =========================================================
#ifndef NCDIO
#ifndef MPI
               write(str,'(a,a,i3.3,a)') &
#ifndef PIXELOUT
                  trim(dgrid(ig)%my%base_name),  '.initl_disp', irupt, '.grd'
#else
                  trim(dgrid(ig)%my%base_name),  '.initl_disp', irupt, '.dat'
#endif
#else
#ifndef ONEFILE
               write(str,'(a,a,i3.3,a,a)') &
                  trim(dgrid(ig)%my%base_name),  '.initl_disp', irupt, '.grd', trim(suffix)
#else
               write(str,'(a,a,i3.3,a)') &
#ifndef PIXELOUT
                  trim(dgrid(ig)%my%base_name),  '.initl_disp', irupt, '.grd'
#else
                  trim(dgrid(ig)%my%base_name),  '.initl_disp', irupt, '.dat'
#endif
#endif
#endif
#ifndef DIROUT
#ifdef MULTI
               str =  trim(members_dir) // trim(str)
#endif
               TIMER_START('maxgrd_write_gmt_init')
               call maxgrd_write_gmt(dgrid(ig)%zz,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.false.)
#else
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.false., &
                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.false.,dgrid(ig),myrank)
#else
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.false.,dgrid(ig),myrank, &
                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
               TIMER_STOP('maxgrd_write_gmt_init')
#else
#ifndef PIXELOUT
               dirname = 'initl_disp.grd'
#else
               dirname = 'initl_disp.dat'
#endif
#ifdef MULTI
               dirname =  trim(members_dir) // trim(dirname)
#endif
               TIMER_START('maxgrd_write_gmt_init')
               call maxgrd_write_gmt(dgrid(ig)%zz,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.false.)
#else
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.false., &
                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.false.,dgrid(ig),myrank)
#else
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.false.,dgrid(ig),myrank, &
                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
               TIMER_STOP('maxgrd_write_gmt_init')
#endif
#else
               TIMER_START('write_initial_displacement')
#if !defined(MPI) || !defined(ONEFILE)
               call write_initial_displacement(dgrid(ig), irupt)
#else
               call write_initial_displacement(dgrid(ig), myrank, irupt)
#endif
               TIMER_STOP('write_initial_displacement')
#endif
! ==============================================================================
               jrupt = irupt
            end if
         end if
      end do
! ==============================================================================
! === SINWAVE ==================================================================
      end if
! ==============================================================================

      do ig = 1, ngrid
! === Support multiple ruptures. ===============================================
         nrupt   => dgrid(ig)%nrupt
         irupt   => dgrid(ig)%irupt
         jrupt   => dgrid(ig)%jrupt
         ruptgrd => dgrid(ig)%ruptgrd
! ==============================================================================
         if(t <= nrupt*tau) then
! === SINWAVE ==================================================================
! === Displacement =============================================================
!           if(init_disp_sinwave == 0) then
! === DEBUG for MPI & multiple rupture & interpolation! ========================
!           if((init_disp_sinwave == 0) .and. (init_disp_fault == 0)) then
            if((init_disp_sinwave == 0) .and. (init_disp_fault == 0) .and. &
               ((init_disp_interpolation /= 1) .or. (ig == 1))) then
! ==============================================================================
! ==============================================================================
! ==============================================================================
            write(6,'(a,i0,a,f0.6,a,f0.6,a,f0.6,a,i0,a,i0,a,a)') &
               'Computing rupture effects: time_step=', istep, ' dt=', dt, ' t=', t, &
! === Support multiple ruptures. ===============================================
               ' tau=', tau, ' ig=', ig, ' irupt=', irupt, ' file=', trim(ruptgrd(irupt))
! ==============================================================================
! === SINWAVE ==================================================================
            else
            write(6,'(a,i0,a,f0.6,a,f0.6,a,f0.6,a,i0,a,i0)') &
               'Computing rupture effects: time_step=', istep, ' dt=', dt, ' t=', t, &
               ' tau=', tau, ' ig=', ig, ' irupt=', irupt
            end if
! ==============================================================================
            TIMER_START('hrise_rwg')
! === SINWAVE ==================================================================
            if(init_disp_sinwave == 0) then
! ==============================================================================
! === When "def_bathy=0", hz is changed on dry cell and it can become wet. =====
!           call hrise_rwg(dgrid(ig)%wave_field,dgrid(ig)%zz,dt,tau,dgrid(ig)%my%nx,dgrid(ig)%my%ny)
            call hrise_rwg(dgrid(ig)%wave_field,dgrid(ig)%zz,dt,tau,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
                           dgrid(ig)%wod_flags,defbathy_flag)
! ==============================================================================
! === SINWAVE ==================================================================
            else
               if(ig == 1) then
#ifndef MPI
                  call hrise_rwg_sin(dgrid(ig)%wave_field,t,dgrid(ig)%my%nx,dgrid(ig)%my%ny)
#else
                  call hrise_rwg_sin(dgrid(ig)%wave_field,t,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
                                     dgrid(ig)%my%kx,dgrid(ig)%my%ky)
#endif
              end if
            end if
! ==============================================================================
            TIMER_STOP('hrise_rwg')
#ifdef BANKFILE
            if(defbathy_flag == 1) then
               TIMER_START('save_dxdy_old')
               call save_dxdy_old(dgrid(ig))
               TIMER_STOP('save_dxdy_old')
            end if
#endif
#ifdef MPI
            call exchange_edges(HGT,dgrid(ig))
#endif
            TIMER_START('drise_rwg')
            call drise_rwg(dgrid(ig)%depth_field,dgrid(ig)%zz,dt,tau,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
                           dgrid(ig)%my%linear_flag,defbathy_flag)
            TIMER_STOP('drise_rwg')
            if(ig /= 1) then
               pid = dgrid(ig)%parent%id
               TIMER_START('interp2fine_dz')
               call interp2fine_dz(dgrid(pid),dgrid(ig))
               TIMER_STOP('interp2fine_dz')
            end if
#ifdef MPI
            call exchange_edges_dz(dgrid(ig))
#endif
#ifdef BANKFILE
            if(defbathy_flag == 1) then
               TIMER_START('update_btxbty')
               call update_btxbty(dgrid(ig))
               TIMER_STOP('update_btxbty')
               TIMER_START('update_bathymetory')
               call update_bathymetory(dgrid(ig))
               TIMER_STOP('update_bathymetory')
            end if
#endif
         end if
! === Write snapshot on step 0 =================================================
         if(plotgrd(1) < 0 .or. plotgrd_num(ig)) then
            if(istep == 1 .and. istep >= itmap_start .and. istep <= itmap_end) then
               ! set pointers for function calls
               wave_field             => dgrid(ig)%wave_field
               depth_field            => dgrid(ig)%depth_field
               ts_field               => dgrid(ig)%ts_field
               niz                    =  dgrid(ig)%my%nx
               njz                    =  dgrid(ig)%my%ny
#ifdef PIXELIN
               nxorg                  =  dgrid(ig)%my%nxorg
               nyorg                  =  dgrid(ig)%my%nyorg
#endif
               wod_flags              => dgrid(ig)%wod_flags
               mlat0                  =  dgrid(ig)%my%mlat0
               mlon0                  =  dgrid(ig)%my%mlon0
               dxdy                   =  dgrid(ig)%my%dh
               base                   => dgrid(ig)%my%base_name
#ifdef MPI
               has_boundary          =  dgrid(ig)%my%has_boundary
#endif
               linear_flag           =  dgrid(ig)%my%linear_flag
#ifndef NCDIO
               TIMER_START('dump_gmt_nl_vel')
#ifndef MPI
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,VEL,linear_flag)
#else
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,VEL,nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,VEL,has_boundary,linear_flag)
#else
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,VEL,has_boundary,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,VEL,has_boundary,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#endif
#endif
               TIMER_STOP('dump_gmt_nl_vel')
               TIMER_START('dump_gmt_nl_hgt')
#ifndef MPI
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,HGT,linear_flag)
#else
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,HGT,nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,HGT,has_boundary,linear_flag)
#else
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,HGT,has_boundary,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,HGT,has_boundary,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#endif
#endif
               TIMER_STOP('dump_gmt_nl_hgt')
#else
               TIMER_START('write_snapshot')
#ifndef MPI
               call write_snapshot(dgrid(ig),REAL_FUNC(0),0,VEL,linear_flag)
               call write_snapshot(dgrid(ig),REAL_FUNC(0),0,HGT,linear_flag)
#else
#ifndef ONEFILE
               call write_snapshot(dgrid(ig),REAL_FUNC(0),0,VEL,has_boundary,linear_flag)
               call write_snapshot(dgrid(ig),REAL_FUNC(0),0,HGT,has_boundary,linear_flag)
#else
               call write_snapshot(dgrid(ig),REAL_FUNC(0),0,VEL,has_boundary,myrank,linear_flag)
               call write_snapshot(dgrid(ig),REAL_FUNC(0),0,HGT,has_boundary,myrank,linear_flag)
#endif
#endif
               TIMER_STOP('write_snapshot')
#endif
            end if
         end if
! ==============================================================================
! === Arrival time =============================================================
         if(check_arrival_time == 1 .and. istep == 1) then
            TIMER_START('check_arrival')
            call check_arrival(wave_field,depth_field,niz,njz,0)
            TIMER_STOP('check_arrival')
         end if
! ==============================================================================
         pid = dgrid(ig)%parent%id
         TIMER_START('tstep_grid_vel')
#ifndef CONV_CHECK
#ifndef CARTESIAN
         call tstep_grid(VEL,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dt,smallh_xy,smallh_wod,c2p_all)
#else
         call tstep_grid(VEL,ig,dgrid(pid),dgrid(ig),cf,cfl,dt,smallh_xy,smallh_wod,c2p_all)
#endif
#else
#ifndef CARTESIAN
         call tstep_grid(VEL,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dt,smallh_xy,smallh_wod,c2p_all,conv_step)
#else
         call tstep_grid(VEL,ig,dgrid(pid),dgrid(ig),cf,cfl,dt,smallh_xy,smallh_wod,c2p_all,conv_step)
#endif
#endif
         TIMER_STOP('tstep_grid_vel')
#ifdef BANKFILE
         if(dgrid(ig)%bank_file /= 'NO_BANK_FILE_GIVEN') then
#ifdef MPI
            TIMER_START('exchange_edges_btxbty')
            call exchange_edges_btxbty(dgrid(ig))
            TIMER_STOP('exchange_edges_btxbty')
#endif
            TIMER_START('update_bathymetory')
            call update_bathymetory(dgrid(ig))
            TIMER_STOP('update_bathymetory')
         end if
#endif
#ifdef CONV_CHECK
#ifdef MPI
         if(myrank == 0) then
#endif
            if((with_disp == 1) .or. (with_disp == 2 .and. ig /= 1)) then
               if(conv_step > max_step) then
                  write(dgrid(ig)%my%nconvout,'(a,i0,a,f0.6,a,i0,a)') &
                     'step=', istep, ' t=', t, ' conv_step=', conv_step, ' [NOT CONVERGED!!!]'
               else
                  write(dgrid(ig)%my%nconvout,'(a,i0,a,f0.6,a,i0)') &
                     'step=', istep, ' t=', t, ' conv_step=', conv_step
               end if
            end if
#ifdef MPI
         end if
#endif
#endif
! === Conversion from flux to velocity should be done right after calc. ========
         base                   => dgrid(ig)%my%base_name
         file_name_bathymetry   => dgrid(ig)%my%bath_file
         displacement_file_name => dgrid(ig)%my%disp_file
         niz                    =  dgrid(ig)%my%nx
         njz                    =  dgrid(ig)%my%ny
#ifdef PIXELIN
         nxorg                  =  dgrid(ig)%my%nxorg
         nyorg                  =  dgrid(ig)%my%nyorg
#endif
         dxdy                   =  dgrid(ig)%my%dh
         mlon0                  =  dgrid(ig)%my%mlon0
         mlat0                  =  dgrid(ig)%my%mlat0
         linear_flag            =  dgrid(ig)%my%linear_flag
         th0                    =  dgrid(ig)%my%th0
         dth                    =  dgrid(ig)%my%dth
         wave_field             => dgrid(ig)%wave_field
         depth_field            => dgrid(ig)%depth_field
         ts_field               => dgrid(ig)%ts_field
         zz                     => dgrid(ig)%zz
         ubnd                   => dgrid(ig)%ubnd
         hbnd                   => dgrid(ig)%hbnd
         hzmax                  => dgrid(ig)%hzmax
#ifdef HZMINOUT
         hzmin                  => dgrid(ig)%hzmin
#endif
#ifndef SKIP_MAX_VEL
         vmax                   => dgrid(ig)%vmax
#endif
         wod_flags              => dgrid(ig)%wod_flags
#ifdef MPI
         has_boundary          =  dgrid(ig)%my%has_boundary
#endif
         linear_flag           =  dgrid(ig)%my%linear_flag

#ifndef SKIP_MAX_VEL
         TIMER_START('maxgrd_v_check_nl')
#ifndef MPI
         call maxgrd_v_check_nl(vmax,wave_field,depth_field,wod_flags,niz,njz,linear_flag)
#else
         call maxgrd_v_check_nl(vmax,wave_field,depth_field,wod_flags,niz,njz,has_boundary,linear_flag)
#endif
         TIMER_STOP('maxgrd_v_check_nl')
#endif

         if(plotgrd(1) < 0 .or. plotgrd_num(ig)) then
            if(mod(istep,itmap) == 0 .and. istep >= itmap_start .and. istep <= itmap_end) then
#ifndef NCDIO
               TIMER_START('dump_gmt_nl_vel')
#ifndef MPI
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,t,istep,base,VEL,linear_flag)
#else
                                mlat0,mlon0,dxdy,t,istep,base,VEL,nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,t,istep,myrank,base,VEL,has_boundary,linear_flag)
#else
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,t,istep,myrank,base,VEL,has_boundary,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,t,istep,myrank,base,VEL,has_boundary,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#endif
#endif
               TIMER_STOP('dump_gmt_nl_vel')
#else
               TIMER_START('write_snapshot_vel')
#ifndef MPI
               call write_snapshot(dgrid(ig),t,istep,VEL,linear_flag)
#else
#ifndef ONEFILE
               call write_snapshot(dgrid(ig),t,istep,VEL,has_boundary,linear_flag)
#else
               call write_snapshot(dgrid(ig),t,istep,VEL,has_boundary,myrank,linear_flag)
#endif
#endif
               TIMER_STOP('write_snapshot_vel')
#endif
            end if
         end if
! ==============================================================================
      end do

      !*** update the wave height ***
      do ig = 1, ngrid
         pid = dgrid(ig)%parent%id
         TIMER_START('tstep_grid_hgt')
#ifndef CONV_CHECK
#ifndef CARTESIAN
         call tstep_grid(HGT,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dt,smallh_xy,smallh_wod,c2p_all)
#else
         call tstep_grid(HGT,ig,dgrid(pid),dgrid(ig),cf,cfl,dt,smallh_xy,smallh_wod,c2p_all)
#endif
#else
#ifndef CARTESIAN
         call tstep_grid(HGT,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dt,smallh_xy,smallh_wod,c2p_all,conv_step)
#else
         call tstep_grid(HGT,ig,dgrid(pid),dgrid(ig),cf,cfl,dt,smallh_xy,smallh_wod,c2p_all,conv_step)
#endif
#endif
         TIMER_STOP('tstep_grid_hgt')

! === recheck_wod should be called after outsea_rwg. by tkato 2012/09/11 =======
!#ifndef MPI
!        if(ig == 1) then
!           TIMER_START('outsea_rwg')
!           call outsea_rwg(dgrid(ig)%wave_field,dgrid(ig)%depth_field,dgrid(ig)%hbnd, &
!                           dgrid(ig)%ubnd,dgrid(ig)%my%nx,dgrid(ig)%my%ny)
!           TIMER_STOP('outsea_rwg')
!        end if
!#endif
! ==============================================================================
      end do

      !*** write outputs (plotgrd<0 => all of them, otherwise just plotgrd number ***
      do ig = 1, ngrid
         ! set pointers for function calls
         base                   => dgrid(ig)%my%base_name
         file_name_bathymetry   => dgrid(ig)%my%bath_file
         displacement_file_name => dgrid(ig)%my%disp_file
         niz                    =  dgrid(ig)%my%nx
         njz                    =  dgrid(ig)%my%ny
#ifdef PIXELIN
         nxorg                  =  dgrid(ig)%my%nxorg
         nyorg                  =  dgrid(ig)%my%nyorg
#endif
         dxdy                   =  dgrid(ig)%my%dh
         mlon0                  =  dgrid(ig)%my%mlon0
         mlat0                  =  dgrid(ig)%my%mlat0
         linear_flag            =  dgrid(ig)%my%linear_flag
         th0                    =  dgrid(ig)%my%th0
         dth                    =  dgrid(ig)%my%dth
         wave_field             => dgrid(ig)%wave_field
         depth_field            => dgrid(ig)%depth_field
         ts_field               => dgrid(ig)%ts_field
         zz                     => dgrid(ig)%zz
         ubnd                   => dgrid(ig)%ubnd
         hbnd                   => dgrid(ig)%hbnd
         hzmax                  => dgrid(ig)%hzmax
#ifdef HZMINOUT
         hzmin                  => dgrid(ig)%hzmin
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
#ifndef SKIP_MAX_VEL
         vmax                   => dgrid(ig)%vmax
#endif
! ==============================================================================
         wod_flags              => dgrid(ig)%wod_flags
! === Conversion from flux to velocity should be done right after calc. ========
#ifdef MPI
         has_boundary          =  dgrid(ig)%my%has_boundary
#endif
         linear_flag           =  dgrid(ig)%my%linear_flag
! ==============================================================================
#ifdef DUMP1D
         TIMER_START('dump1d')
         call dump1d(ig, dgrid(ig), istep)
         TIMER_STOP('dump1d')
#endif

! === Arrival time =============================================================
         if(check_arrival_time == 1) then
            call check_arrival(wave_field,depth_field,niz,njz,istep)
         end if
! ==============================================================================
         !*** check for maximum wave heights ***
         TIMER_START('maxgrd_check_nl')
         call maxgrd_check_nl(hzmax,wave_field,wod_flags,niz,njz)
         TIMER_STOP('maxgrd_check_nl')
#ifdef HZMINOUT
         TIMER_START('mingrd_check_nl')
         call mingrd_check_nl(hzmin,wave_field,wod_flags,niz,njz)
         TIMER_STOP('mingrd_check_nl')
#endif

         ! Burbidge: Stop things if the maximum wave height gets silly
         TIMER_START('error_check')
         error = 0
#ifndef __NEC__
!$omp parallel do private(i)
#else
!$omp parallel do private(i) reduction(+:error)
#endif
         do j = 1, njz
            do i = 1, niz
               if(hzmax(i,j) > 1.0d6) then
#ifndef __NEC__
!$omp critical
                  error = 1
!$omp end critical
#else
                  error = error + 1
#endif
               end if
            end do
         end do
! === Finalized if error occurred. =============================================
#ifdef MPI
#ifndef MULTI
         call MPI_Allreduce(MPI_IN_PLACE, error, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
         call MPI_Allreduce(MPI_IN_PLACE, error, 1, MPI_INTEGER, MPI_SUM, MPI_MEMBER_WORLD, ierr)
#endif
#endif
         if(error /= 0) then
            TIMER_STOP('error_check')
            write(6,'(a)')      '======================================================='
            write(6,'(a)')      '======================================================='
            write(6,'(a)')      '=== ERROR!!!'
            write(6,'(a)')      '=== Maximum wave height exceeds 1 million'
            write(6,'(a,i3,a)') '=== in domain ID ', ig, '!'
            write(6,'(a)')      '=== Probable instability!'
            write(6,'(a,i8,a)') '=== Execution is stopped at step ', istep, ' and finalized!'
            write(6,'(a)')      '======================================================='
            write(6,'(a)')      '======================================================='
            exit
         end if
! ==============================================================================
         TIMER_STOP('error_check')

         if(plotgrd(1) < 0 .or. plotgrd_num(ig)) then
            if(mod(istep,itmap) == 0 .and. istep >= itmap_start .and. istep <= itmap_end) then
#ifndef NCDIO
! === Conversion from flux to velocity should be done right after calc. ========
!              TIMER_START('dump_gmt_nl')
               TIMER_START('dump_gmt_nl_hgt')
! ==============================================================================
#ifndef MPI
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,t,istep,base,HGT,linear_flag)
#else
                                mlat0,mlon0,dxdy,t,istep,base,HGT,nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,t,istep,myrank,base,HGT,has_boundary,linear_flag)
#else
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,t,istep,myrank,base,HGT,has_boundary,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,t,istep,myrank,base,HGT,has_boundary,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#endif
#endif
! === Conversion from flux to velocity should be done right after calc. ========
!              TIMER_STOP('dump_gmt_nl')
               TIMER_STOP('dump_gmt_nl_hgt')
! ==============================================================================
#else
! === Conversion from flux to velocity should be done right after calc. ========
!              call write_snapshot(dgrid(ig),t,istep)
               TIMER_START('write_snapshot_hgt')
#ifndef MPI
               call write_snapshot(dgrid(ig),t,istep,HGT,linear_flag)
#else
#ifndef ONEFILE
               call write_snapshot(dgrid(ig),t,istep,HGT,has_boundary,linear_flag)
#else
               call write_snapshot(dgrid(ig),t,istep,HGT,has_boundary,myrank,linear_flag)
#endif
#endif
               TIMER_STOP('write_snapshot_hgt')
! ==============================================================================
#endif
            end if
         end if
      end do

! === Finalized if error occurred. =============================================
      if(error /= 0) exit
! ==============================================================================

      ! tgs text output
! === To decimate outputs by tkato. 2016/10/24 =================================
      if(mod(istep,itgrn) == 0) then
! ==============================================================================
#ifndef SINGLE_TGS
      do j = 1, nsta
         inquire(tgsfp(j), opened=tgs_opened)
         if(tgs_opened) then
            tgs_nx   =  dgrid(mytgs(j)%ig)%my%nx
            tgs_ny   =  dgrid(mytgs(j)%ig)%my%ny
            tgs_wfld => dgrid(mytgs(j)%ig)%wave_field
            tgs_lon  =  mytgs(j)%geolon
            tgs_lat  =  mytgs(j)%geolat
            tgs_z    =  mytgs(j)%z
            tgs_i    =  mytgs(j)%ilon
            tgs_j    =  mytgs(j)%ilat
            tgs_fx   =  tgs_wfld%fx(tgs_i, tgs_j)
            tgs_fy   =  tgs_wfld%fy(tgs_i, tgs_j)
            tgs_hz   =  tgs_wfld%hz(tgs_i, tgs_j)
            if(istep/itgrn == 1) then
               write(tgsfp(j), '(a,i0,a,f0.3,a,f0.3,a,f0.3,a,i0)') '> TGS No.=', mytgs(j)%number, ' lat=', tgs_lat, &
                  ' lon=', tgs_lon, ' depth=', tgs_z, ' grid_id=', mytgs(j)%ig
! === DEBUG for tgs text output by tkato. 2012/10/11 ===========================
!              write(tgsfp(j), '(a,i0,a,e,a,e,a,e,a,e)') 'step=', istep, ' t=', t, ' hz=', tgs_hz, ' fx=', tgs_fx, ' fy=', tgs_fy
! ==============================================================================
            end if
! === DEBUG for tgs text output by tkato. 2012/10/11 ===========================
#ifndef __FUJITSU
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
#ifndef REAL_DBLE
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#endif
#endif
! ==============================================================================
         end if
      end do
#else
      do j = 1, nsta
         inquire(tgsfp, opened=tgs_opened)
         if(tgs_opened) then
            tgs_nx   =  dgrid(mytgs(j)%ig)%my%nx
            tgs_ny   =  dgrid(mytgs(j)%ig)%my%ny
            tgs_wfld => dgrid(mytgs(j)%ig)%wave_field
            tgs_lon  =  mytgs(j)%geolon
            tgs_lat  =  mytgs(j)%geolat
            tgs_z    =  mytgs(j)%z
            tgs_i    =  mytgs(j)%ilon
            tgs_j    =  mytgs(j)%ilat
            tgs_fx   =  tgs_wfld%fx(tgs_i, tgs_j)
            tgs_fy   =  tgs_wfld%fy(tgs_i, tgs_j)
            tgs_hz   =  tgs_wfld%hz(tgs_i, tgs_j)
            write(tgsnum,'(a,i6.6,a)') '[', mytgs(j)%number, ']'
            if(istep/itgrn == 1) then
               write(tgsfp, '(a,i0,a,f0.3,a,f0.3,a,f0.3,a,i0)') trim(tgsnum) // '> TGS No.=', mytgs(j)%number, &
                  ' lat=', tgs_lat, ' lon=', tgs_lon, ' depth=', tgs_z, ' grid_id=', mytgs(j)%ig
            end if
#ifndef __FUJITSU
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
#ifndef REAL_DBLE
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#endif
#endif
         end if
      end do
#endif
! === To decimate outputs by tkato. 2016/10/24 =================================
      end if
! ==============================================================================

      !*** Burbidge - print a message to stdout every 100 time steps ***
      if(mod(istep,100) == 0) then
         write(6,'(a,i0,a,f0.3,a)') 'Timestep ', istep, ' completed. Time =', t, 's'
      end if
! === Support restart ==========================================================
      if(restart_interval /= 0) then
         if(mod(istep,restart_interval) == 0) then
#ifndef MPI
            write(restart_file_name, '(a,i8.8)') 'restart.', istep
#else
            write(restart_file_name, '(a,i8.8,a,i6.6)') 'restart.', istep, '.', myrank
#endif
            write(6,'(a,a)') '[RESTART] Restart file written: ', trim(restart_file_name)
#ifdef MULTI
            restart_file_name = trim(members_dir) // trim(restart_file_name)
#endif
            TIMER_START('write_restart_file')
            call write_restart_file(ngrid, dgrid, restart_file_name)
            TIMER_STOP('write_restart_file')
         end if
      end if
! ==============================================================================
   end do ! end of main loop
! === Restart output for truncation ============================================
      if(restart_interval /= 0 .and. trunc_flag) then
         istep = istep - 1
#ifndef MPI
         write(restart_file_name, '(a,i8.8)') 'restart.', istep
#else
         write(restart_file_name, '(a,i8.8,a,i6.6)') 'restart.', istep, '.', myrank
#endif
         write(6,'(a,a)') '[RESTART] Restart file written: ', trim(restart_file_name)
#ifdef MULTI
         restart_file_name = trim(members_dir) // trim(restart_file_name)
#endif
         TIMER_START('write_restart_file')
         call write_restart_file(ngrid, dgrid, restart_file_name)
         TIMER_STOP('write_restart_file')
      end if
! ==============================================================================

   !*** dump max wave height grid ***
   do ig = 1, ngrid
! === Multi-grids can be specified! ============================================
!     if(plotgrd_num < 0 .or. plotgrd_num == ig) then
      if(plotgrd(1) < 0 .or. plotgrd_num(ig)) then
! ==============================================================================
#ifndef NCDIO
#ifndef MPI
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(max_grid_file_name)
#else
#ifndef ONEFILE
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(max_grid_file_name) // trim(suffix)
#else
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(max_grid_file_name)
#endif
#endif
         TIMER_START('maxgrd_write_gmt')
#ifndef DIROUT
#ifdef MULTI
         str = trim(members_dir) // trim(str)
#endif
         call maxgrd_write_gmt(dgrid(ig)%hzmax,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
! === For negative max. height =================================================
!                              dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str)
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true., &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig),myrank)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig),myrank, &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
!===============================================================================
#else
         dirname = trim(max_grid_file_name)
#ifdef MULTI
         dirname = trim(members_dir) // trim(dirname)
#endif
         call maxgrd_write_gmt(dgrid(ig)%hzmax,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
! === For negative max. height =================================================
!                              dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str)
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true., &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig),myrank)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig),myrank, &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
!===============================================================================
#endif
         TIMER_STOP('maxgrd_write_gmt')
#else
         TIMER_START('write_max_height')
#if !defined(MPI) || !defined(ONEFILE)
         call write_max_height(dgrid(ig))
#else
         call write_max_height(dgrid(ig), myrank)
#endif
         TIMER_STOP('write_max_height')
#endif
#ifdef HZMINOUT
#ifndef NCDIO
#ifndef MPI
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(min_grid_file_name)
#else
#ifndef ONEFILE
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(min_grid_file_name) // trim(suffix)
#else
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(min_grid_file_name)
#endif
#endif
         TIMER_START('mingrd_write_gmt')
#ifndef DIROUT
#ifdef MULTI
         str = trim(members_dir) // trim(str)
#endif
         call mingrd_write_gmt(dgrid(ig)%hzmin,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
! === For negative min. height =================================================
!                              dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str)
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true., &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig),myrank)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig),myrank, &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
!===============================================================================
#else
         dirname = trim(min_grid_file_name)
#ifdef MULTI
         dirname = trim(members_dir) // trim(dirname)
#endif
         call mingrd_write_gmt(dgrid(ig)%hzmin,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
! === For negative min. height =================================================
!                              dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str)
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true., &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig),myrank)
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig),myrank, &
                               dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
!===============================================================================
#endif
         TIMER_STOP('mingrd_write_gmt')
#else
         TIMER_START('write_min_height')
#if !defined(MPI) || !defined(ONEFILE)
         call write_min_height(dgrid(ig))
#else
         call write_min_height(dgrid(ig), myrank)
#endif
         TIMER_STOP('write_min_height')
#endif
#endif
! === Arrival time =============================================================
         if(check_arrival_time == 1) then
            TIMER_START('calc_arrival_time')
            call calc_arrival_time(dgrid(ig)%wave_field,dgrid(ig)%my%nx,dgrid(ig)%my%ny,dt)
            TIMER_STOP('calc_arrival_time')
#ifndef NCDIO
#ifndef MPI
#ifndef PIXELOUT
            str = trim(dgrid(ig)%my%base_name) // '.arrival_time.grd'
#else
            str = trim(dgrid(ig)%my%base_name) // '.arrival_time.dat'
#endif
#else
#ifndef ONEFILE
            str = trim(dgrid(ig)%my%base_name) // '.arrival_time.grd' // trim(suffix)
#else
#ifndef PIXELOUT
            str = trim(dgrid(ig)%my%base_name) // '.arrival_time.grd'
#else
            str = trim(dgrid(ig)%my%base_name) // '.arrival_time.dat'
#endif
#endif
#endif
            TIMER_START('maxgrd_write_gmt_arrival')
#ifndef DIROUT
#ifdef MULTI
            str = trim(members_dir) // trim(str)
#endif
            call maxgrd_write_gmt(dgrid(ig)%wave_field%arrival_time,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.)
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true., &
                                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig),myrank)
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig),myrank, &
                                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
#else
#ifndef PIXELOUT
            dirname = 'arrival_time.grd'
#else
            dirname = 'arrival_time.dat'
#endif
#ifdef MULTI
            dirname = trim(members_dir) // trim(dirname)
#endif
            call maxgrd_write_gmt(dgrid(ig)%wave_field%arrival_time,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.)
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true., &
                                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig),myrank)
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig),myrank, &
                                  dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
#endif
            TIMER_STOP('maxgrd_write_gmt_arrival')
#else
            TIMER_START('write_arrival_time')
#if !defined(MPI) || !defined(ONEFILE)
            call write_arrival_time(dgrid(ig))
#else
            call write_arrival_time(dgrid(ig), myrank)
#endif
            TIMER_STOP('write_arrival_time')
#endif
         end if
! ==============================================================================
      end if
   end do
! === To add max velocity output. by tkato 2012/10/02 ==========================
#ifndef SKIP_MAX_VEL
   !*** dump max velocity grid ***
   do ig = 1, ngrid
! === Multi-grids can be specified! ============================================
!     if(plotgrd_num < 0 .or. plotgrd_num == ig) then
      if(plotgrd(1) < 0 .or. plotgrd_num(ig)) then
! ==============================================================================
#ifndef NCDIO
#ifndef MPI
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(vmax_grid_file_name)
#else
#ifndef ONEFILE
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(vmax_grid_file_name) // trim(suffix)
#else
         str = trim(dgrid(ig)%my%base_name) // '.' // trim(vmax_grid_file_name)
#endif
#endif
         TIMER_START('maxgrd_v_write_gmt')
#ifndef DIROUT
#ifdef MULTI
         str = trim(members_dir) // trim(str)
#endif
         call maxgrd_v_write_gmt(dgrid(ig)%vmax,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str)
#else
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str, &
                                 dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,dgrid(ig),myrank)
#else
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,dgrid(ig),myrank, &
                                 dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
#else
         dirname = trim(vmax_grid_file_name)
#ifdef MULTI
         dirname = trim(members_dir) // trim(dirname)
#endif
         call maxgrd_v_write_gmt(dgrid(ig)%vmax,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str)
#else
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str, &
                                 dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
#ifndef PIXELOUT
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,dgrid(ig),myrank)
#else
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,dgrid(ig),myrank, &
                                 dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#endif
#endif
         TIMER_STOP('maxgrd_v_write_gmt')
#else
         TIMER_START('write_max_velocity')
#if !defined(MPI) || !defined(ONEFILE)
         call write_max_velocity(dgrid(ig))
#else
         call write_max_velocity(dgrid(ig), myrank)
#endif
         TIMER_STOP('write_max_velocity')
#endif
      end if
   end do
#endif
! ==============================================================================
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
   if(with_elastic_loading == 1) then
      do ig = 1, ngrid
         TIMER_START('loading_finalize')
         call loading_finalize(dgrid(ig))
         TIMER_STOP('loading_finalize')
      end do
   end if
! ==============================================================================
#endif
#ifdef NCDIO
   do ig = 1, ngrid
! === Output only specified domains. ===========================================
! === Multi-grids can be specified! ============================================
!     if(plotgrd_num < 0 .or. plotgrd_num == ig) then
      if(plotgrd(1) < 0 .or. plotgrd_num(ig)) then
! ==============================================================================
! ==============================================================================
#if !defined(MPI) || !defined(ONEFILE)
      call close_file(dgrid(ig))
#else
      call close_file(dgrid(ig), myrank)
#endif
! === Output only specified domains. ===========================================
      end if
! ==============================================================================
   end do
#endif
#ifdef CONV_CHECK
#ifdef MPI
   if(myrank == 0) then
#endif
      if(with_disp == 1 .or. with_disp == 2) then
         do ig = 1, ngrid
            if(with_disp == 1 .or. ig /= 1) then
               close(dgrid(ig)%my%nconvout)
            end if
         end do
      end if
#ifdef MPI
   end if
#endif
#endif

   !*** close all open tide gauge files ***
   if(nsta > 0) then
#ifndef SINGLE_TGS
      deallocate(mytgs)
      ! tgs text output
      do j = 1, nsta
         inquire(tgsfp(j), opened=tgs_opened)
         if(tgs_opened) then
            close(tgsfp(j))
         end if
      end do
      deallocate(tgsfiles)
      deallocate(tgsfp)
#else
      close(tgsfp)
#endif
   end if
! === Support multiple ruptures. ===============================================
   do ig = 1, ngrid
      deallocate(dgrid(ig)%ruptgrd)
   end do
! ==============================================================================
! === Don't repeat allocate/deallocate! ========================================
#ifdef MPI
   do ig = 1, ngrid
      call deallocate_edges(dgrid(ig))
   end do
#endif
! ==============================================================================
#ifdef DUMP1D
   call dump1d_finalize()
#endif

#ifdef MPI
#ifndef MULTI
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
#else
   call MPI_Barrier(MPI_MEMBER_WORLD, ierr)
#endif
#endif
   TIMER_STOP('All')
#ifdef TIMER
   call print_timer()
#endif
#ifdef MPI
! === Separate stdout into for each process. ===================================
   close(6)
! ==============================================================================
#ifdef USE_ALLTOALLV
#ifdef A2A3D
   call A2A3D_finalize()
#endif
#endif
   call MPI_Finalize(ierr)
#else
#ifdef MULTI
   call MPI_Finalize(ierr)
#endif
#endif
   stop

#ifdef __SX__
300 write(0,'(a,a)') 'read_grid_info(): cant open file ', trim(gridfile)
#ifndef MPI
   stop
#else
   ierr = 1
   call fatal_error(201)
#endif
500 write(0,'(a,a)') 'missing tide gauge station ', trim(tg_station_file_name)
   nsta = 0
   goto 503
501 write(0,'(a,a)') 'missing or malformed tide gauge station file ', trim(tg_station_file_name)
   nsta = 0
   goto 503
#endif
100 write(0,'(a,a)') 'Error opening rupture list file ', trim(dgrid(1)%my%disp_file)
#ifndef MPI
   stop
#else
   ierr = ierr + 1
   goto 110
#endif
end program JAGURS
