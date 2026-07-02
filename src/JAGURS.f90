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
   use mod_init_disp_pointsource
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
   use mod_interpolation, only : interp2fine_init_disp, interp2fine_init_hz, &
                                 interp2fine_init_fx, interp2fine_init_fy
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
#ifdef NORMALMODE
#ifndef NM_DEF
   use mod_normalmode, only : normalmode_read_namelist, make_nm_ind, calc_nm_P
#else
   use mod_normalmode, only : normalmode_read_namelist, normalmode_set_params, calc_nm_P
#endif
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
   character(len=256), pointer :: file_name_init_hz, file_name_init_fx, file_name_init_fy

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
   integer(kind=4) :: myrank = 0
#ifdef MPI
   ! MPI
   integer(kind=4) :: nprocs = 0
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
#ifdef NORMALMODE
   real(kind=REAL_BYTE) :: tgs_P
#endif
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
   integer(kind=4) :: conv_step
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
   character(len=256) :: pointsource_file_name = 'pointsource'
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
   integer(kind=4) :: nxorg, nyorg
   integer(kind=4) :: formatid
   real(kind=REAL_BYTE), pointer, dimension(:,:) :: tmp
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
   pointsource_file_name = trim(pointsource_file_name) // trim(suffix)

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
   call MPI_Comm_size(__MPICOMM__, nprocs, ierr)
   call MPI_Comm_rank(__MPICOMM__, myrank, ierr)
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
   call read_grid_info(gridfile,dgrid,ngrid)
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
   if(init_disp_pointsource == 1) then
      TIMER_START('specify_pointsource_params')
      call specify_pointsource_params(pointsource_file_name)
      TIMER_STOP('specify_pointsource_params')
   end if
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
   if((init_disp_interpolation == 1) .or. (elastic_loading_interpolation == 1) .or. (init_val_interpolation == 1)) then
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
            if(trim(dgrid(ig)%my%disp_file(1:15)) == 'NO_DISPLACEMENT') then
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
         if(trim(dgrid(ig)%my%disp_file(1:15)) == 'NO_DISPLACEMENT') then
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
      call read_bathymetry_gmt_grdhdr(dgrid(ig)%my%bath_file,niz,njz,dxdy,mlon0,mlat0,nxorg,nyorg,formatid)
      TIMER_STOP('read_bathymetry_gmt_grdhdr')
#if defined(MPI) && defined(ONEFILE)
      end if

      call MPI_Bcast(niz,   1, MPI_INTEGER, 0, __MPICOMM__, ierr)
      call MPI_Bcast(njz,   1, MPI_INTEGER, 0, __MPICOMM__, ierr)
      call MPI_Bcast(dxdy,  1, REAL_MPI,    0, __MPICOMM__, ierr)
      call MPI_Bcast(mlon0, 1, REAL_MPI,    0, __MPICOMM__, ierr)
      call MPI_Bcast(mlat0, 1, REAL_MPI,    0, __MPICOMM__, ierr)
#ifdef PIXELIN
      call MPI_Bcast(nxorg, 1, MPI_INTEGER, 0, __MPICOMM__, ierr)
      call MPI_Bcast(nyorg, 1, MPI_INTEGER, 0, __MPICOMM__, ierr)
#endif
      call MPI_Bcast(formatid, 1, MPI_INTEGER, 0, __MPICOMM__, ierr)

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
      dgrid(ig)%my%formatid = formatid
 
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

      call MPI_Allreduce(niz, totalNx, 1, MPI_INTEGER, MPI_SUM, __MPICOMM__, ierr)
      totalNx = totalNx/npy - (npx-1)*2
      dgrid(ig)%my%totalNx = totalNx

      call MPI_Allreduce(njz, totalNy, 1, MPI_INTEGER, MPI_SUM, __MPICOMM__, ierr)
      totalNy = totalNy/npx - (npy-1)*2
      dgrid(ig)%my%totalNy = totalNy

      call MPI_Allreduce(mlon0, glon0, 1, REAL_MPI, MPI_MIN, __MPICOMM__, ierr)
      dgrid(ig)%my%glon0 = glon0

      call MPI_Allreduce(mlat0, glat0, 1, REAL_MPI, MPI_MIN, __MPICOMM__, ierr)
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
      call MPI_Bcast(dgrid(ig)%my%dh, 1, REAL_MPI, 0, __MPICOMM__, ierr)
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
      if(timenest == 1) then
         allocate(dgrid(ig)%wave_field%fx_a(-2:niz+1,-1:njz+1))
         allocate(dgrid(ig)%wave_field%fx_b(-2:niz+1,-1:njz+1))
         allocate(dgrid(ig)%wave_field%fy_a(-1:niz+1,-2:njz+1))
         allocate(dgrid(ig)%wave_field%fy_b(-1:niz+1,-2:njz+1))
         allocate(dgrid(ig)%wave_field%fx_i2f0(-2:niz+1,-1:njz+1))
         allocate(dgrid(ig)%wave_field%fx_i2f1(-2:niz+1,-1:njz+1))
         allocate(dgrid(ig)%wave_field%fy_i2f0(-1:niz+1,-2:njz+1))
         allocate(dgrid(ig)%wave_field%fy_i2f1(-1:niz+1,-2:njz+1))
      end if
#ifdef BANKFILE
      if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
      if(timenest == 1) then
         allocate(dgrid(ig)%wave_field%fx_a(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%fx_b(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%fx_i2f0(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%fx_i2f1(ist:ien,jst:jen))
      end if
#ifdef BANKFILE
      if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
      if(timenest == 1) then
         allocate(dgrid(ig)%wave_field%fy_a(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%fy_b(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%fy_i2f0(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%fy_i2f1(ist:ien,jst:jen))
      end if
#ifdef BANKFILE
      if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
      if(timenest == 1) then
         allocate(dgrid(ig)%wave_field%hz_a(-1:niz+2,-1:njz+2))
         allocate(dgrid(ig)%wave_field%hz_b(-1:niz+2,-1:njz+2))
         allocate(dgrid(ig)%wave_field%hz_i2f0(-1:niz+2,-1:njz+2))
         allocate(dgrid(ig)%wave_field%hz_i2f1(-1:niz+2,-1:njz+2))
      end if
      allocate(dgrid(ig)%depth_field%dz   (-1:niz+2,-1:njz+2))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
#ifdef NORMALMODE
      allocate(dgrid(ig)%wave_field%nm_ind(niz,njz))
      allocate(dgrid(ig)%wave_field%nm_P(niz,njz))
      allocate(dgrid(ig)%wave_field%nm_P0(niz,njz))
      allocate(dgrid(ig)%wave_field%nm_P1(niz,njz))
#endif
#else
      ist = 0
      ien = niz + 2
      if(iand(has_boundary, WEST_BOUND) /= 0) ist = -1
      jst = 0
      jen = njz + 2
      if(iand(has_boundary, NORTH_BOUND) /= 0) jst = -1
      allocate(dgrid(ig)%wave_field%hz    (ist:ien,jst:jen))
      if(timenest == 1) then
         allocate(dgrid(ig)%wave_field%hz_a(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%hz_b(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%hz_i2f0(ist:ien,jst:jen))
         allocate(dgrid(ig)%wave_field%hz_i2f1(ist:ien,jst:jen))
      end if
#ifndef NONESTDEBUG
      allocate(dgrid(ig)%wave_field%noi2f (ist:ien,jst:jen))
#endif
      allocate(dgrid(ig)%wave_field%hz_old(ist:ien,jst:jen))
      allocate(dgrid(ig)%depth_field%dz   (ist:ien,jst:jen))
#ifdef BANKFILE
      if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
#ifdef NORMALMODE
      allocate(dgrid(ig)%wave_field%nm_ind(niz,njz))
      allocate(dgrid(ig)%wave_field%nm_P(niz,njz))
      allocate(dgrid(ig)%wave_field%nm_P0(niz,njz))
      allocate(dgrid(ig)%wave_field%nm_P1(niz,njz))
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
      if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
      if(check_tt_time == 1) then
         allocate(dgrid(ig)%wave_field%tttdat(niz,njz))
         allocate(dgrid(ig)%wave_field%tt_time(niz,njz))
      end if
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
#ifndef USE_GPU
      dgrid(ig)%depth_field%dx = 0.0d0 ! depth field
      dgrid(ig)%depth_field%dy = 0.0d0 ! depth field
      dgrid(ig)%depth_field%dz = 0.0d0 ! depth field
      dgrid(ig)%bcf_field = 0.0d0 ! friction field
#else
      tmp => dgrid(ig)%depth_field%dx
!$omp target teams distribute parallel do collapse(2) private(i)
      do j = lbound(tmp,2), ubound(tmp,2)
         do i = lbound(tmp,1), ubound(tmp,1)
            tmp(i,j) = 0.0d0
         end do
      end do

      tmp => dgrid(ig)%depth_field%dy
!$omp target teams distribute parallel do collapse(2) private(i)
      do j = lbound(tmp,2), ubound(tmp,2)
         do i = lbound(tmp,1), ubound(tmp,1)
            tmp(i,j) = 0.0d0
         end do
      end do

      tmp => dgrid(ig)%depth_field%dz
!$omp target teams distribute parallel do collapse(2) private(i)
      do j = lbound(tmp,2), ubound(tmp,2)
         do i = lbound(tmp,1), ubound(tmp,1)
            tmp(i,j) = 0.0d0
         end do
      end do

      tmp => dgrid(ig)%bcf_field
!$omp target teams distribute parallel do collapse(2) private(i)
      do j = lbound(tmp,2), ubound(tmp,2)
         do i = lbound(tmp,1), ubound(tmp,1)
            tmp(i,j) = 0.0d0
         end do
      end do
#endif
#ifdef BANKFILE
      if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
      if(timenest == 1) then
         tmp => dgrid(ig)%wave_field%fx
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(tmp,2), ubound(tmp,2)
            do i = lbound(tmp,1), ubound(tmp,1)
               dgrid(ig)%wave_field%fx_a(i,j) = 0.0d0
               dgrid(ig)%wave_field%fx_b(i,j) = 0.0d0
               dgrid(ig)%wave_field%fx_i2f0(i,j) = 0.0d0
               dgrid(ig)%wave_field%fx_i2f1(i,j) = 0.0d0
            end do
         end do

         tmp => dgrid(ig)%wave_field%fy
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(tmp,2), ubound(tmp,2)
            do i = lbound(tmp,1), ubound(tmp,1)
               dgrid(ig)%wave_field%fy_a(i,j) = 0.0d0
               dgrid(ig)%wave_field%fy_b(i,j) = 0.0d0
               dgrid(ig)%wave_field%fy_i2f0(i,j) = 0.0d0
               dgrid(ig)%wave_field%fy_i2f1(i,j) = 0.0d0
            end do
         end do

         tmp => dgrid(ig)%wave_field%hz
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(tmp,2), ubound(tmp,2)
            do i = lbound(tmp,1), ubound(tmp,1)
               dgrid(ig)%wave_field%hz_a(i,j) = 0.0d0
               dgrid(ig)%wave_field%hz_b(i,j) = 0.0d0
               dgrid(ig)%wave_field%hz_i2f0(i,j) = 0.0d0
               dgrid(ig)%wave_field%hz_i2f1(i,j) = 0.0d0
            end do
         end do
      end if

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
      call read_bathymetry_gmt_grd(file_name_bathymetry, depth_field, niz, njz, linear_flag, dgrid(ig), myrank, nxorg, nyorg, dgrid(ig)%my%formatid)
      TIMER_STOP('read_bathymetry_gmt_grd')
#ifdef BANKFILE
      TIMER_START('read_bank_file')
      call read_bank_file(dgrid(ig))
      TIMER_STOP('read_bank_file')
#endif

      !*** read in the bottom friction coefficient ***
      TIMER_START('read_friction_gmt_grd')
      call read_friction_gmt_grd(bcf_file_name, bcf_field, niz, njz, dgrid(ig), myrank, nxorg, nyorg, dgrid(ig)%my%formatid)
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
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field,dgrid(ig)%my%formatid)
#else
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field,nxorg,nyorg,dgrid(ig)%my%formatid)
#endif
#else
#ifndef PIXELIN
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field,dgrid(ig),myrank,dgrid(ig)%my%formatid)
#else
      call wet_or_dry(wave_field,depth_field,wod_flags,niz,njz,wod_file_name,wod_field,dgrid(ig),myrank,nxorg,nyorg,dgrid(ig)%my%formatid)
#endif
#endif
      TIMER_STOP('wet_or_dry')

      if(timenest /= 1) then
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
      if(timenest == 1) then
         dgrid(ig)%my%nl = min(dgrid(pid)%my%nl + 1, max_nest_level)
         dgrid(ig)%my%maxnl = dgrid(ig)%my%nl
      end if

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
#ifndef NOCHECK
#ifndef MPI
         stop
#else
         write(0,'(a)') '  *** from global origin (deg):\n'
         write(0,'(7x,a,f13.8,a,f13.8)') 'mlon0= ', mlon0, ' input value= ', dgrid(ig)%my%glon0/60.0d0
         write(0,'(7x,a,f13.8,a,f13.8)') 'mlat0= ', mlat0, ' input value= ', dgrid(ig)%my%glat0/60.0d0
         ierr = ierr + 1
         call fatal_error(110)
#endif
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
#ifndef NOCHECK
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
#ifndef NOCHECK
         stop
#endif
#else
         write(0,'(a,i0,a,i0,a)') '****** (nx - 1 = ', dgrid(ig)%my%totalNx-1, ' ) / (nr = ', dgrid(ig)%my%nr, ' ) must be zero'
         write(0,'(a,i0,a,i0,a)') '****** (ny - 1 = ', dgrid(ig)%my%totalNy-1, ' ) / (nr = ', dgrid(ig)%my%nr, ' ) must be zero'
#ifndef NOCHECK
         ierr = ierr + 1
         call fatal_error(111)
#endif
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
                         dgrid(ig)%my%nx,dgrid(ig)%my%ny,dgrid(ig)%wod_file,dgrid(ig)%wod_field,dgrid(ig)%my%formatid)
#else
                         dgrid(ig)%my%nx,dgrid(ig)%my%ny,dgrid(ig)%wod_file,dgrid(ig)%wod_field, &
                         dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg)
#endif
#else
                         dgrid(ig)%my%nx,dgrid(ig)%my%ny,dgrid(ig)%wod_file,dgrid(ig)%wod_field, &
#ifndef PIXELIN
                         dgrid(ig),myrank,dgrid(ig)%my%formatid)
#else
                         dgrid(ig),myrank,dgrid(ig)%my%nxorg,dgrid(ig)%my%nyorg,dgrid(ig)%my%formatid)
#endif
#endif
         TIMER_STOP('wet_or_dry')
      end if

#ifdef MPI
      end if
      if(ierr > 0) call fatal_error(109)
#endif
   end do

   if(timenest == 1) then
      do ig = 1, ngrid-1
         dgrid(ig)%my%maxnl = dgrid(ngrid)%my%maxnl
      end do
      do ig = 1, ngrid
         dgrid(ig)%my%dt = dt*3.0**(dgrid(ig)%my%maxnl - dgrid(ig)%my%nl)
      end do

      TIMER_START('boundary_rwg')
      ig = 1
      depth_field  => dgrid(ig)%depth_field
      hbnd         => dgrid(ig)%hbnd
      ubnd         => dgrid(ig)%ubnd
      th0          =  dgrid(ig)%my%th0
      dxdy         =  dgrid(ig)%my%dh
      dth          =  dgrid(ig)%my%dth
      niz          =  dgrid(ig)%my%nx
      njz          =  dgrid(ig)%my%ny
#ifdef MPI
      iy           =  dgrid(ig)%my%iy
      has_boundary =  dgrid(ig)%my%has_boundary
#endif
#ifndef MPI
#ifndef CARTESIAN
      call boundary_rwg(depth_field,hbnd,ubnd,dgrid(ig)%my%dt,th0,dth,niz,njz)
#else
      call boundary_rwg(depth_field,hbnd,ubnd,dgrid(ig)%my%dt,dxdy,niz,njz)
#endif
#else
#ifndef CARTESIAN
      call boundary_rwg(depth_field,hbnd,ubnd,dgrid(ig)%my%dt,th0,dth,max(iy-2,0),niz,njz,has_boundary)
#else
      call boundary_rwg(depth_field,hbnd,ubnd,dgrid(ig)%my%dt,dxdy,niz,njz,has_boundary)
#endif
#endif
      TIMER_STOP('boundary_rwg')
   end if

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
   call tgs_open_ct(dgrid,ngrid,dt*REAL_FUNC(itgrn),nstep/itgrn,nsta,program_name,tg_station_file_name,mytgs)
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
#ifdef NORMALMODE
   call normalmode_read_namelist(dt)
   do ig = 1, ngrid
#ifndef NM_DEF
      call make_nm_ind(dgrid(ig)%my%nx,dgrid(ig)%my%ny,dgrid(ig)%my%mlat0,dgrid(ig)%my%mlon0,dgrid(ig)%my%dh, &
                       dgrid(ig)%wave_field%nm_ind)
#else
      call normalmode_set_params(dgrid(ig)%my%mlat0,dgrid(ig)%my%mlon0,dgrid(ig)%my%dh)
#endif
   end do
#endif
   !*** main loop ***
   if(timenest == 1) then
      do ig = 1, ngrid
         tmp => dgrid(ig)%wave_field%fx
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(tmp,2), ubound(tmp,2)
            do i = lbound(tmp,1), ubound(tmp,1)
               dgrid(ig)%wave_field%fx_a(i,j) = dgrid(ig)%wave_field%fx(i,j)
               dgrid(ig)%wave_field%fx_b(i,j) = dgrid(ig)%wave_field%fx(i,j)
            end do
         end do

         tmp => dgrid(ig)%wave_field%fy
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(tmp,2), ubound(tmp,2)
            do i = lbound(tmp,1), ubound(tmp,1)
               dgrid(ig)%wave_field%fy_a(i,j) = dgrid(ig)%wave_field%fy(i,j)
               dgrid(ig)%wave_field%fy_b(i,j) = dgrid(ig)%wave_field%fy(i,j)
            end do
         end do

         tmp => dgrid(ig)%wave_field%hz
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
         do j = lbound(tmp,2), ubound(tmp,2)
            do i = lbound(tmp,1), ubound(tmp,1)
               dgrid(ig)%wave_field%hz_a(i,j) = dgrid(ig)%wave_field%hz(i,j)
               dgrid(ig)%wave_field%hz_b(i,j) = dgrid(ig)%wave_field%hz(i,j)
            end do
         end do
      end do
   end if
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

   do ig = 1, ngrid
      file_name_init_hz => dgrid(ig)%init_hz_file
      file_name_init_fx => dgrid(ig)%init_fx_file
      file_name_init_fy => dgrid(ig)%init_fy_file
      wave_field        => dgrid(ig)%wave_field
      depth_field       => dgrid(ig)%depth_field
      niz               =  dgrid(ig)%my%nx
      njz               =  dgrid(ig)%my%ny
      linear_flag       =  dgrid(ig)%my%linear_flag
      wod_flags         => dgrid(ig)%wod_flags
#ifdef PIXELIN
      nxorg             =  dgrid(ig)%my%nxorg
      nyorg             =  dgrid(ig)%my%nyorg
#endif

      if((init_val_interpolation == 0) .or. (ig == 1)) then
         if(file_name_init_hz(1:10) /= 'NO_INIT_HZ') then
            call read_init_val_gmt_grd(file_name_init_hz, wave_field, niz, njz, linear_flag, dgrid(ig), myrank, nxorg, nyorg, dgrid(ig)%my%formatid, HGT)
#ifdef MPI
            call exchange_edges(HGT,dgrid(ig))
#endif
            if(linear_flag == 0) then
               call recheck_wod(wave_field, depth_field, wod_flags, niz, njz, smallh_xy)
#ifdef MPI
               call exchange_edges(HGT,dgrid(ig))
               call exchange_edges_wod(dgrid(ig))
#endif
            end if
         end if

         if(file_name_init_fx(1:10) /= 'NO_INIT_FX') then
            call read_init_val_gmt_grd(file_name_init_fx, wave_field, niz, njz, linear_flag, dgrid(ig), myrank, nxorg, nyorg, dgrid(ig)%my%formatid, IFX)

            if(linear_flag == 0) then
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
               do j = 1, njz
                  do i = 1, niz
                     if(wod_flags(i,j) /= 1) wave_field%fx(i,j) = 0.0d0
                  end do
               end do
            end if

#ifdef MPI
            call exchange_edges(VEL,dgrid(ig))
#endif
         end if

         if(file_name_init_fy(1:10) /= 'NO_INIT_FY') then
            call read_init_val_gmt_grd(file_name_init_fy, wave_field, niz, njz, linear_flag, dgrid(ig), myrank, nxorg, nyorg, dgrid(ig)%my%formatid, IFY)

            if(linear_flag == 0) then
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
               do j = 1, njz
                  do i = 1, niz
                     if(wod_flags(i,j) /= 1) wave_field%fy(i,j) = 0.0d0
                  end do
               end do
            end if

#ifdef MPI
            call exchange_edges(VEL,dgrid(ig))
#endif
         end if
      else
         pid = dgrid(ig)%parent%id

         call interp2fine_init_hz(dgrid(pid), dgrid(ig))

#ifdef MPI
         call exchange_edges(HGT,dgrid(ig))
#endif
         if(linear_flag == 0) then
            call recheck_wod(wave_field, depth_field, wod_flags, niz, njz, smallh_xy)
#ifdef MPI
            call exchange_edges(HGT,dgrid(ig))
            call exchange_edges_wod(dgrid(ig))
#endif
         end if

         call interp2fine_init_fx(dgrid(pid), dgrid(ig))
         call interp2fine_init_fy(dgrid(pid), dgrid(ig))
!
         if(linear_flag == 0) then
#ifndef USE_GPU
!$omp parallel do private(i)
#else
!$omp target teams distribute parallel do collapse(2) private(i)
#endif
            do j = 1, njz
               do i = 1, niz
                  if(wod_flags(i,j) /= 1) then
                     wave_field%fx(i,j) = 0.0d0
                     wave_field%fy(i,j) = 0.0d0
                  end if
               end do
            end do
         end if

#ifdef MPI
         call exchange_edges(VEL,dgrid(ig))
#endif
      end if
   end do

   do istep = istart, nstep
      do ig = 1, ngrid
         dgrid(ig)%my%numneststeps = 3**(dgrid(ig)%my%maxnl - dgrid(ig)%my%nl)
         dgrid(ig)%my%neststephgt  = mod(istep-1,                               dgrid(ig)%my%numneststeps)
         dgrid(ig)%my%neststepvel  = mod(istep-1 + dgrid(ig)%my%numneststeps/2, dgrid(ig)%my%numneststeps)
         dgrid(ig)%my%calchgt = 0
         if(dgrid(ig)%my%neststephgt == 0) dgrid(ig)%my%calchgt = 1
         dgrid(ig)%my%calcvel = 0
         if(dgrid(ig)%my%neststepvel == 0) dgrid(ig)%my%calcvel = 1
      end do
! ==============================================================================
      t = REAL_FUNC(istep) * dt
! === Support truncation =======================================================
      if(max_time_i /= 0) then
#ifndef MPI
         trunc_flag = check_trunc(max_time_i)
#else
         if(myrank == 0) trunc_flag = check_trunc(max_time_i)
         call MPI_Bcast(trunc_flag, 1, MPI_LOGICAL, 0, __MPICOMM__, ierr)
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
               if((init_disp_gaussian == 0) .and. (init_disp_pointsource == 0)) then
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
                  if(init_disp_gaussian == 1) then
                     TIMER_START('make_gaussian_rupture')
                     call make_gaussian_rupture(dgrid(ig))
                     TIMER_STOP('make_gaussian_rupture')
                  end if
                  if(init_disp_pointsource == 1) then
                     TIMER_START('make_pointsource_rupture')
                     call make_pointsource_rupture(dgrid(ig))
                     TIMER_STOP('make_pointsource_rupture')
                  end if
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
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.false.,dgrid(ig))
#else
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.false.,dgrid(ig), &
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
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.false.,dgrid(ig))
#else
                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.false.,dgrid(ig), &
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
            if(timenest /= 1) then
               call hrise_rwg(dgrid(ig)%wave_field,dgrid(ig)%zz,dt,tau,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
                              dgrid(ig)%wod_flags,defbathy_flag)
            else
               if(dgrid(ig)%my%calchgt == 1) then
                  call hrise_rwg(dgrid(ig)%wave_field,dgrid(ig)%zz,dgrid(ig)%my%dt,tau,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
                                 dgrid(ig)%wod_flags,defbathy_flag)
               end if
            end if
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
            if(timenest /= 1) then
               call exchange_edges(HGT,dgrid(ig))
            else
               call exchange_edges_b(HGT,dgrid(ig))
            end if
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
#ifdef NORMALMODE
                if(dumpp == 1) then
#ifndef NM_DEF
                   call calc_nm_P(niz, njz, wave_field%nm_ind, wave_field%nm_P, wave_field%nm_P0, wave_field%nm_P1, 1)
#else
                   call calc_nm_P(niz, njz, wave_field%nm_P, .true.)
#endif
                end if
#endif
#ifndef NCDIO
               TIMER_START('dump_gmt_nl_vel')
#ifndef MPI
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef PIXELOUT
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,VEL,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,VEL,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,VEL,has_boundary,dgrid(ig),linear_flag)
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
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,HGT,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,base,HGT,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,REAL_FUNC(0),0,myrank,base,HGT,has_boundary,dgrid(ig),linear_flag)
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
         if(check_tt_time == 1 .and. istep == 1) then
            TIMER_START('check_ttt')
            call check_ttt(wave_field,depth_field,niz,njz,0)
            TIMER_STOP('check_ttt')
         end if
         pid = dgrid(ig)%parent%id
         TIMER_START('tstep_grid_vel')
         if(timenest /= 1) then
            call tstep_grid(VEL,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dt,smallh_xy,smallh_wod,c2p_all,conv_step,istep)
         else
            call tstep_grid(VEL,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dgrid(ig)%my%dt,smallh_xy,smallh_wod,c2p_all,conv_step,istep)
         end if
         TIMER_STOP('tstep_grid_vel')
#ifdef BANKFILE
         if(dgrid(ig)%bank_file(1:7) /= 'NO_BANK') then
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
                                mlat0,mlon0,dxdy,t,istep,base,VEL,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,t,istep,base,VEL,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,t,istep,myrank,base,VEL,has_boundary,dgrid(ig),linear_flag)
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
         if(timenest /= 1) then
            call tstep_grid(HGT,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dt,smallh_xy,smallh_wod,c2p_all,conv_step,istep)
         else
            call tstep_grid(HGT,ig,dgrid(pid),dgrid(ig),cf,cfl,coriolis,dgrid(ig)%my%dt,smallh_xy,smallh_wod,c2p_all,conv_step,istep)
         end if
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
         if(check_tt_time == 1) then
            call check_ttt(wave_field,depth_field,niz,njz,istep)
         end if
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
#ifndef USE_GPU
#ifndef __NEC__
!$omp parallel do private(i)
#else
!$omp parallel do private(i) reduction(+:error)
#endif
#else
!$omp target teams distribute parallel do collapse(2) private(i) reduction(+:error)
#endif
         do j = 1, njz
            do i = 1, niz
               if(hzmax(i,j) > 1.0d6) then
#if !defined(__NEC__) && !defined(USE_GPU)
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
         call MPI_Allreduce(MPI_IN_PLACE, error, 1, MPI_INTEGER, MPI_SUM, __MPICOMM__, ierr)
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
                                mlat0,mlon0,dxdy,t,istep,base,HGT,dgrid(ig),linear_flag)
#else
                                mlat0,mlon0,dxdy,t,istep,base,HGT,dgrid(ig),nxorg,nyorg,linear_flag)
#endif
#else
               call dump_gmt_nl(wave_field,depth_field,ts_field,niz,njz,wod_flags, &
#ifndef ONEFILE
                                mlat0,mlon0,dxdy,t,istep,myrank,base,HGT,has_boundary,dgrid(ig),linear_flag)
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
#ifdef NORMALMODE
            tgs_P    =  tgs_wfld%nm_P(tgs_i, tgs_j)
#endif
            if(istep/itgrn == 1) then
               write(tgsfp(j), '(a,i0,a,f0.3,a,f0.3,a,f0.3,a,i0)') '> TGS No.=', mytgs(j)%number, ' lat=', tgs_lat, &
                  ' lon=', tgs_lon, ' depth=', tgs_z, ' grid_id=', mytgs(j)%ig
! === DEBUG for tgs text output by tkato. 2012/10/11 ===========================
!              write(tgsfp(j), '(a,i0,a,e,a,e,a,e,a,e)') 'step=', istep, ' t=', t, ' hz=', tgs_hz, ' fx=', tgs_fx, ' fy=', tgs_fy
! ==============================================================================
            end if
! === DEBUG for tgs text output by tkato. 2012/10/11 ===========================
#ifndef __FUJITSU
#ifndef NORMALMODE
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy, ' P= ', tgs_P
#endif
#else
#ifndef REAL_DBLE
#ifndef NORMALMODE
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy, ' P= ', tgs_P
#endif
#else
#ifndef NORMALMODE
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp(j), '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy, ' P= ', tgs_P
#endif
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
#ifdef NORMALMODE
            tgs_P    =  tgs_wfld%nm_P(tgs_i, tgs_j)
#endif
            write(tgsnum,'(a,i6.6,a)') '[', mytgs(j)%number, ']'
            if(istep/itgrn == 1) then
               write(tgsfp, '(a,i0,a,f0.3,a,f0.3,a,f0.3,a,i0)') trim(tgsnum) // '> TGS No.=', mytgs(j)%number, &
                  ' lat=', tgs_lat, ' lon=', tgs_lon, ' depth=', tgs_z, ' grid_id=', mytgs(j)%ig
            end if
#ifndef __FUJITSU
#ifndef NORMALMODE
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy, ' P= ', tgs_P
#endif
#else
#ifndef REAL_DBLE
#ifndef NORMALMODE
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy, ' P= ', tgs_P
#endif
#else
#ifndef NORMALMODE
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy
#else
            write(tgsfp, '(a,i0,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3,a,e23.15e3)') &
               trim(tgsnum) // 'step=', istep, ' t=', t, ' hz= ', tgs_hz, ' fx= ', tgs_fx, ' fy= ', tgs_fy, ' P= ', tgs_P
#endif
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
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig))
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig), &
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
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig))
#else
                               dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig), &
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
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig))
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig), &
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
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig))
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig), &
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
         if(check_tt_time == 1) then
            TIMER_START('calc_tt_time')
            call calc_tt_time(dgrid(ig)%wave_field,dgrid(ig)%my%nx,dgrid(ig)%my%ny,dt)
            TIMER_STOP('calc_tt_time')
#ifndef NCDIO
#ifndef MPI
#ifndef PIXELOUT
            str = trim(dgrid(ig)%my%base_name) // '.' // trim(ttt_file) // '.grd'
#else
            str = trim(dgrid(ig)%my%base_name) // '.' // trim(ttt_file) // '.dat'
#endif
#else
#ifndef ONEFILE
            str = trim(dgrid(ig)%my%base_name) // '.' // trim(ttt_file) // '.dat' // trim(suffix)
#else
#ifndef PIXELOUT
            str = trim(dgrid(ig)%my%base_name) // '.' // trim(ttt_file) // '.grd'
#else
            str = trim(dgrid(ig)%my%base_name) // '.' // trim(ttt_file) // '.dat'
#endif
#endif
#endif
            TIMER_START('maxgrd_write_gmt_ttt')
#ifndef DIROUT
#ifdef MULTI
            str = trim(members_dir) // trim(str)
#endif
            call maxgrd_write_gmt(dgrid(ig)%wave_field%tt_time,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig))
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,.true.,dgrid(ig), &
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
            dirname = trim(ttt_file) // '.grd'
#else
            dirname = trim(ttt_file) // '.dat'
#endif
#ifdef MULTI
            dirname = trim(members_dir) // trim(dirname)
#endif
            call maxgrd_write_gmt(dgrid(ig)%wave_field%tt_time,dgrid(ig)%my%nx,dgrid(ig)%my%ny, &
#if !defined(MPI) || !defined(ONEFILE)
#ifndef PIXELOUT
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig))
#else
                                  dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,.true.,dgrid(ig), &
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
            TIMER_STOP('maxgrd_write_gmt_ttt')
#else
            TIMER_START('write_tt_time')
#if !defined(MPI) || !defined(ONEFILE)
            call write_tt_time(dgrid(ig))
#else
            call write_tt_time(dgrid(ig), myrank)
#endif
            TIMER_STOP('write_tt_time')
#endif
         end if
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
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,dgrid(ig))
#else
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,str,dgrid(ig), &
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
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,dgrid(ig))
#else
                                 dgrid(ig)%my%mlon0,dgrid(ig)%my%mlat0,dgrid(ig)%my%dh,dirname,str,dgrid(ig), &
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
   call MPI_Barrier(__MPICOMM__, ierr)
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

100 write(0,'(a,a)') 'Error opening rupture list file ', trim(dgrid(1)%my%disp_file)
#ifndef MPI
   stop
#else
   ierr = ierr + 1
   goto 110
#endif
end program JAGURS
