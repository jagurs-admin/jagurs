#include "real.h"
module mod_params
#ifdef MPI
use mod_mpi_fatal
#endif
! === Support truncation =======================================================
use mod_truncation
! ==============================================================================
#ifdef MULTI
! === Split Dir ================================================================
!  use mod_multi, only : member_id, input_files
   use mod_multi, only : member_id, input_files, input_dirname
! ==============================================================================
#endif
implicit none

integer(kind=4), parameter :: INI = 0
integer(kind=4), parameter :: VEL = 1
integer(kind=4), parameter :: HGT = 2
#ifdef MPI
integer(kind=4), parameter :: DXY = 3
integer(kind=4), parameter :: DPT = 4
#endif

!*** check undefined values ***
#if !defined(PIXELIN) && defined(NONESTDEBUG)
integer(kind=4),   private, parameter :: IUDEF = 999999
real(kind=REAL_BYTE),      private, parameter :: RUDEF = 999999.9d0
#else
integer(kind=4),            parameter :: IUDEF = 999999
real(kind=REAL_BYTE),               parameter :: RUDEF = 999999.9d0
#endif
character(len=32), private, parameter :: SUDEF = ''

character(len=256) :: program_name

character(len=128) :: gridfile = SUDEF
real(kind=REAL_BYTE) :: dt = RUDEF, tau = RUDEF, tend = RUDEF
integer(kind=4) :: itmap = IUDEF
real(kind=REAL_BYTE) :: cf = RUDEF ! when cf> 0: it is a non-dimensional frictional coefficient
                           ! when cf< 0: -cf (>0) is Manning's roughness coefficient which
                           !                produces depth-dependent friction
real(kind=REAL_BYTE) :: cfl = 0.0d0 ! frictional coefficient for land
integer(kind=4) :: coriolis = IUDEF ! coriolis=1: coriolis force included
                                    ! coriolis=0: coriolis force neglected
! === Output file name should be optional. =====================================
!character(len=128) :: maxgrdfn = SUDEF, max_grid_file_name
#ifndef PIXELOUT
character(len=128) :: maxgrdfn = 'zmax.grd', max_grid_file_name
#else
character(len=128) :: maxgrdfn = 'zmax.dat', max_grid_file_name
#endif
! ==============================================================================
#ifdef HZMINOUT
! === Output file name should be optional. =====================================
!character(len=128) :: mingrdfn = SUDEF, min_grid_file_name
#ifndef PIXELOUT
character(len=128) :: mingrdfn = 'zmin.grd', min_grid_file_name
#else
character(len=128) :: mingrdfn = 'zmin.dat', min_grid_file_name
#endif
! ==============================================================================
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
! === Output file name should be optional. =====================================
!character(len=128) :: vmaxgrdfn = SUDEF, vmax_grid_file_name
#ifndef PIXELOUT
character(len=128) :: vmaxgrdfn = 'vmax.grd', vmax_grid_file_name
#else
character(len=128) :: vmaxgrdfn = 'vmax.dat', vmax_grid_file_name
#endif
! ==============================================================================
! ==============================================================================
character(len=128) :: tgstafn = SUDEF, tg_station_file_name
character(len=256) :: tgsoutfile = SUDEF
integer(kind=4) :: pointers = 1 ! if( pointers <= 0 ) don't pointers
integer(kind=4) :: smooth_edges = 0
!*** Burbidge ***
integer(kind=4) :: velmux = 1, velmux_flag ! If 1 then out velocity mux files, else do not output
integer(kind=4) :: velgrd = 1, velgrd_flag ! If 1 then out velocity grd files, else do not output
integer(kind=4) :: speedgrd = 0, speedgrd_flag ! If 1 then out speed grd files, else do not output
! === Multi-grids can be specified! ============================================
!integer(kind=4) :: plotgrd = -1, plotgrd_num ! If < 0 then output all grd files, else grdnumber equal to plotgrd_num
integer(kind=4), parameter, private :: ngrid_max = 16
integer(kind=4), dimension(ngrid_max) :: plotgrd = -1
logical, dimension(ngrid_max) :: plotgrd_num = .false.
integer(kind=4), private :: i
! ==============================================================================
! === Write snapshot on step 0 =================================================
!integer(kind=4) :: itmap_start = 1 ! Output grd files after itmap_start steps
integer(kind=4) :: itmap_start = 0 ! Output grd files after itmap_start steps
! ==============================================================================
integer(kind=4) :: itmap_end = 99999999 ! steps at Output grd files end
integer(kind=4) :: def_bathy = 1, defbathy_flag ! If 1 then deform bathy, otherwise skip that step
! === Flood Change =============================================================
!real(kind=REAL_BYTE) :: smallh_xy = 0.01d0  ! Shifted the smallh parameters to here, to help with initial instability
!real(kind=REAL_BYTE) :: smallh_wod = 0.01d0
real(kind=REAL_BYTE) :: smallh_xy = 1.0d-6  ! Shifted the smallh parameters to here, to help with initial instability
real(kind=REAL_BYTE) :: smallh_wod = 1.0d-6
! ==============================================================================
!*** thomas ***
integer(kind=4) :: muxfmt = IUDEF ! 0 = single grid format, 1 = old variable grid format, 2 = small variable grid format
real(kind=REAL_BYTE) :: htmin = 0.0d0 ! when abs(waveheight) exceeds this then recording begins
real(kind=REAL_BYTE) :: tgrt = -1.0d0 ! recording duration in seconds
integer(kind=4) :: itgrn = 1 ! tide gauge data will be dumped every itgrn steps
integer(kind=4) :: multrupt = 0 ! 0 = single rupture, 1 = multiple ruptures over time
integer(kind=4) :: smallmux ! computed to hold 0 or 1
#ifdef MPI
integer(kind=4) :: npx
integer(kind=4) :: npy
integer(kind=4), private :: procx = IUDEF, procy = IUDEF ! temp values
#endif
integer(kind=4) :: tgstxtout = 1
character(len=128) :: tgstxtoutfile = 'tgs'
! === Config of copy from child to perent. by tkato 2012/11/15 =================
integer(kind=4) :: c2p_all = 0
! ==============================================================================
! === Support restart ==========================================================
integer(kind=4) :: restart = 0, restart_interval = 0
! ==============================================================================
! === Dispersive ===============================================================
integer(kind=4) :: with_disp = 0             ! 0: without disp. 1: with disp.
! === Modified on 2013/10/01 ===================================================
!integer(kind=4) :: max_step = 1000           ! Max. convergence loop
integer(kind=4) :: max_step = 9999           ! Max. convergence loop
! ==============================================================================
real(kind=REAL_BYTE) :: min_depth = 5.0d0    ! Minimum depth dispersive is calculated
#ifdef CONV_CHECK
! === Modified on 2013/10/01 ===================================================
!real(kind=REAL_BYTE) :: conv_val = 1.0d-6    ! Conv val.
real(kind=REAL_BYTE) :: conv_val = 1.0d-8    ! Conv val.
! ==============================================================================
#endif
! ==============================================================================
! === Absorbing boundary condition =============================================
integer(kind=4) :: with_abc = 0
integer(kind=4) :: nxa = 20, nya = 20
real(kind=REAL_BYTE) :: apara = 0.055d0
! ==============================================================================
! === Support start date in time attribute. ====================================
character(len=64) :: start_date = '2000-01-01 00:00:00'
! ==============================================================================
! === Support truncation =======================================================
character(len=64) :: max_time = ''
integer(kind=4) :: max_time_i = 0
! ==============================================================================
! === For negative max. height =================================================
real(kind=REAL_BYTE), parameter :: missing_value = -1.0d10
! ==============================================================================
! === For MRI ==================================================================
integer(kind=4) :: init_disp_gaussian = 0
! ==============================================================================
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
integer(kind=4) :: with_elastic_loading = 0
real(kind=8) :: m_radius = 2000.0d0 ! Radius in which the loading effect will be estimated.
character(len=256) :: m_pyfile = 'PREM_Ggz.nc'
!character(len=256) :: m_pyfile = 'PREM_Gz.nc'
! ==============================================================================
! === Density ==================================================================
integer(kind=4) :: with_density = 0
real(kind=8) :: m_rho = 1025.5d0
real(kind=8) :: m_K = 2.2d9
! ==============================================================================
#endif
! === SINWAVE ==================================================================
integer(kind=4) :: init_disp_sinwave = 0
! ==============================================================================
! === Limiter with max Froude number. ==========================================
real(kind=REAL_BYTE) :: froude_lim = 2.0d0
! ==============================================================================
! === Displacement =============================================================
integer(kind=4) :: init_disp_fault = 0
integer(kind=4) :: hzdisp_effect = 0
integer(kind=4) :: apply_kj_filter = 0 ! Apply Kajiura filter to zz.
character(len=256) :: fault_param_file = 'fault'
! ==============================================================================
! === 1-way nest ===============================================================
integer(kind=4) :: nest_1way = 0
! ==============================================================================
! === Initial displacement of child domains is given by interpolation. =========
integer(kind=4) :: init_disp_interpolation = 0
integer(kind=4) :: use_linear = 0
! ==============================================================================
! === Specify lower limit of depth to adopt horizontal displacement effect. ====
real(kind=REAL_BYTE) :: min_depth_hde = 50.0d0
! ==============================================================================
! === Arrival time =============================================================
integer(kind=4) :: check_arrival_time = 0
real(kind=REAL_BYTE) :: check_arrival_height = 0.01d0 ! [m]
! ==============================================================================
! === Elastic loading with interpolation =======================================
integer(kind=4) :: elastic_loading_interpolation = 1
! ==============================================================================
#ifdef BANKFILE
real(kind=REAL_BYTE) :: broken_rate = -1.0d0
#endif

contains

   subroutine getpar()
#ifndef MULTI
      integer(kind=4) :: num_arg, len_arg
#endif
      character(len=128) :: arg
      namelist /params/ gridfile, dt, tau, tend, itmap, cf, cfl, coriolis, &
! === To add max velocity output. by tkato 2012/10/02 ==========================
!        maxgrdfn, tgstafn, tgsoutfile, pointers, smooth_edges, &
         maxgrdfn, vmaxgrdfn, tgstafn, tgsoutfile, pointers, smooth_edges, &
! ==============================================================================
#ifdef HZMINOUT
         mingrdfn, &
#endif
      !** Burbidge - GA's input flags **
         velmux, velgrd, plotgrd, itmap_start, itmap_end, def_bathy, smallh_xy, smallh_wod, &
      ! thomas - additional parameters for version that can write smaller mux files
         muxfmt, & ! 0=single grid format, 1=old variable grid format, 2=small variable grid format
         htmin,  & ! if muxfmt=2 then recording at  tide gauge begins when abs wave ht exceeds this
         tgrt,   & ! if muxfmt=2 then duration in seconds for recording at tide gauge
         itgrn,  & ! output tide gauge data every itgrn time steps
      ! thomas - additional parameters for multiple rupture version
#ifdef MPI
         procx, procy, tgstxtout, tgstxtoutfile, &
#else
         tgstxtout, tgstxtoutfile, &
#endif
! === Support restart ==========================================================
         restart, restart_interval, &
! ==============================================================================
! === Dispersive ===============================================================
         with_disp, max_step, min_depth, &
#ifdef CONV_CHECK
         conv_val, &
#endif
! ==============================================================================
! === Absorbing boundary condition =============================================
         with_abc, nxa, nya, apara, &
! ==============================================================================
! === Config of copy from child to perent. by tkato 2012/11/15 =================
!        multrupt ! 1 = multiple ruptures, 0 = single rupture
         multrupt, & ! 1 = multiple ruptures, 0 = single rupture
         c2p_all, &
! ==============================================================================
! === Support truncation =======================================================
         max_time, &
! ==============================================================================
! === Support start date in time attribute. ====================================
         start_date, &
! ==============================================================================
#ifndef CARTESIAN
! === For MRI ==================================================================
         init_disp_gaussian, &
! === SINWAVE ==================================================================
         init_disp_sinwave, &
! ==============================================================================
! === Elastic Loading ==========================================================
!        speedgrd
         speedgrd, &
! === Density ==================================================================
!        with_elastic_loading, m_radius, m_pyfile
         with_elastic_loading, m_radius, m_pyfile, &
! === Limiter with max Froude number. ==========================================
!        with_density, m_rho, m_K
         with_density, m_rho, m_K, &
! === Displacement =============================================================
!        froude_lim
! === 1-way nest ===============================================================
!        froude_lim, init_disp_fault, hzdisp_effect, apply_filter, fault_param_file
         froude_lim, init_disp_fault, hzdisp_effect, apply_kj_filter, fault_param_file, &
! === Initial displacement of child domains is given by interpolation. =========
!        nest_1way
! === Specify lower limit of depth to adopt horizontal displacement effect. ====
!        nest_1way, init_disp_interpolation, use_linear
         nest_1way, init_disp_interpolation, use_linear, &
#ifndef BANKFILE
! === Arrival time =============================================================
!        min_depth_hde
! === Elastic loading with interpolation =======================================
!        min_depth_hde, check_arrival_time, check_arrival_height
         min_depth_hde, check_arrival_time, check_arrival_height, &
         elastic_loading_interpolation
! ==============================================================================
! ==============================================================================
#else
! === Arrival time =============================================================
!        min_depth_hde, broken_rate
! === Elastic loading with interpolation =======================================
!        min_depth_hde, broken_rate, check_arrival_time, check_arrival_height
         min_depth_hde, broken_rate, check_arrival_time, check_arrival_height, &
         elastic_loading_interpolation
! ==============================================================================
! ==============================================================================
#endif
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
#else
! === For MRI ==================================================================
         init_disp_gaussian, &
! ==============================================================================
! === SINWAVE ==================================================================
         init_disp_sinwave, &
! === Limiter with max Froude number. ==========================================
!        speedgrd
         speedgrd, &
! === Displacement =============================================================
!        froude_lim
! === 1-way nest ===============================================================
!        froude_lim, init_disp_fault, hzdisp_effect, apply_filter, fault_param_file
         froude_lim, init_disp_fault, hzdisp_effect, apply_kj_filter, fault_param_file, &
! === Initial displacement of child domains is given by interpolation. =========
!        nest_1way
! === Specify lower limit of depth to adopt horizontal displacement effect. ====
!        nest_1way, init_disp_interpolation, use_linear
         nest_1way, init_disp_interpolation, use_linear, &
#ifndef BANKFILE
! === Arrival time =============================================================
!        min_depth_hde
! === Elastic loading with interpolation =======================================
!        min_depth_hde, check_arrival_time, check_arrival_height
         min_depth_hde, check_arrival_time, check_arrival_height, &
         elastic_loading_interpolation
! ==============================================================================
! ==============================================================================
#else
! === Arrival time =============================================================
!        min_depth_hde, broken_rate
! === Elastic loading with interpolation =======================================
!        min_depth_hde, broken_rate, check_arrival_time, check_arrival_height
         min_depth_hde, broken_rate, check_arrival_time, check_arrival_height, &
         elastic_loading_interpolation
! ==============================================================================
! ==============================================================================
#endif
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
! ==============================================================================
#endif

      call get_command_argument(0, program_name)

#ifndef MULTI
      num_arg = command_argument_count()
      if(num_arg /= 1) then
         call print_usage()
         stop
      end if

      call get_command_argument(1, arg)
      len_arg = len_trim(arg)
      if(arg(1:4) .eq. 'par=') then
         arg = trim(arg(5:len_arg))
      else
         call print_usage()
         stop
      end if
#else
! === Split Dir ================================================================
!     arg = input_files(member_id+1)
      arg = trim(input_dirname) // trim(input_files(member_id+1))
! ==============================================================================
#endif

      open(1,file=trim(arg),action='read',status='old',form='formatted')
      read(1,params)
      close(1)

      max_grid_file_name   = maxgrdfn
#ifdef HZMINOUT
      min_grid_file_name   = mingrdfn
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
      vmax_grid_file_name  = vmaxgrdfn
! ==============================================================================
      tg_station_file_name = tgstafn
      velmux_flag          = velmux
      velgrd_flag          = velgrd
      speedgrd_flag        = speedgrd
! === Multi-grids can be specified! ============================================
!     plotgrd_num          = plotgrd
      do i = 1, ngrid_max
         if(plotgrd(i) == -1) exit
         if((plotgrd(i) >= 1) .and. (plotgrd(i) <= ngrid_max)) then
            plotgrd_num(plotgrd(i)) = .true.
         end if
      end do
! ==============================================================================
      defbathy_flag        = def_bathy
#ifdef MPI
      npx = procx
      npy = procy
#endif
! === Support truncation =======================================================
      max_time_i = truncation_time(max_time)
! ==============================================================================

      call check_mandatory_params()

      call check_params()
! === Dispersive ===============================================================
      if(with_disp == 1 .or. with_disp == 2) then
         c2p_all =  1 ! Dispersive needs c2p_all
! === ABC is optional! =========================================================
!        with_abc = 1 !  Dispersive needs with_abc
! ==============================================================================
      end if
! ==============================================================================

      return
   end subroutine getpar

   subroutine check_params()
      if(cfl == 0.0d0) cfl = cf ! if cfl is not defined, use cf instead

      !*** Burbidge - Error check the input ***
      if(dt <= 0.0d0) then
         write(0,'(a)') 'Error in tsun.par file:'
         write(0,'(a)') 'Timestep is not positive!'
         write(0,'(a,f18.6)') 'dt=', dt
#ifndef MPI
         stop
#else
         call fatal_error(103)
#endif
      end if

      if(dt > tau) then
! === tau < dt is okay now =====================================================
!         write(0,'(a)') 'Error in tsun.par file:'
!         write(0,'(a)') 'Timestep is greater than uplift time!'
!         write(0,'(a,f18.6,a,f18.6)') 'dt=', dt, ' tau=', tau
!#ifndef MPI
!         stop
!#else
!         call fatal_error(104)
!#endif
         write(6,'(a)') '[NOTE] Specified uplift time "tau" is less than timestep "dt".'
         write(6,'(a)') '       So, "tau" is reset to "dt"!'
         tau = dt
! ==============================================================================
      end if

      ! thomas - check that tau is (close to) an integer multiple of dt,
      !          otherwise the full rise will not be achieved
      if(abs(tau/dt - anint(tau/dt)) > 0.001d0) then
         write(0,'(a)') 'Error in tsun.par file:'
         write(0,'(a,f18.6,a,f18.6)') 'Uplift time=', tau, ' is not a multiple of time step=', dt
#ifndef MPI
         stop
#else
         call fatal_error(105)
#endif
       end if

      if(dt > tend) then
         write(0,'(a)') 'Error in tsun.par file:'
         write(0,'(a)') 'Timestep is greater than final time!'
         write(0,'(a,f18.6,a,f18.6)') 'dt=', dt, ' tend=', tend
#ifndef MPI
         stop
#else
         call fatal_error(106)
#endif
      end if

      return
   end subroutine check_params

   subroutine check_mandatory_params()
      logical :: error = .false.
      if(gridfile == SUDEF)   then; write(0,'(a)') 'ERROR!: gridfile must be specified!';   error = .true.; end if
      if(dt == RUDEF)         then; write(0,'(a)') 'ERROR!: dt must be specified!';         error = .true.; end if
      if(tau == RUDEF)        then; write(0,'(a)') 'ERROR!: tau must be specified!';        error = .true.; end if
      if(tend == RUDEF)       then; write(0,'(a)') 'ERROR!: tend must be specified!';       error = .true.; end if
      if(itmap == IUDEF)      then; write(0,'(a)') 'ERROR!: itmap must be specified!';      error = .true.; end if
      if(cf == RUDEF)         then; write(0,'(a)') 'ERROR!: cf must be specified!';         error = .true.; end if
      if(coriolis == IUDEF)   then; write(0,'(a)') 'ERROR!: coriolis must be specified!';   error = .true.; end if
! === Output file name should be optional. =====================================
!     if(maxgrdfn == SUDEF)   then; write(0,'(a)') 'ERROR!: maxgrdfn must be specified!';   error = .true.; end if
! ==============================================================================
#ifdef HZMINOUT
! === Output file name should be optional. =====================================
!     if(mingrdfn == SUDEF)   then; write(0,'(a)') 'ERROR!: mingrdfn must be specified!';   error = .true.; end if
! ==============================================================================
#endif
! === To add max velocity output. by tkato 2012/10/02 ==========================
! === Output file name should be optional. =====================================
!     if(vmaxgrdfn == SUDEF)  then; write(0,'(a)') 'ERROR!: vmaxgrdfn must be specified!';  error = .true.; end if
! ==============================================================================
! ==============================================================================
      if(tgstafn == SUDEF)    then; write(0,'(a)') 'ERROR!: tgstafn must be specified!';    error = .true.; end if
#ifdef MPI
      if(procx == IUDEF)      then; write(0,'(a)') 'ERROR!: procx must be specified!';      error = .true.; end if
      if(procy == IUDEF)      then; write(0,'(a)') 'ERROR!: procy must be specified!';      error = .true.; end if
#endif
      if(error) stop
      return
   end subroutine check_mandatory_params

   subroutine print_usage()
      write(0,'(a)') '[Usage]'
      write(0,'(a,a)') trim(program_name), ' par=[parameter file]'
      return
   end subroutine print_usage

#ifndef MPI
   subroutine putpar(gfile, ngrid, dg)
#else
   subroutine putpar(nprocs, gfile, ngrid, dg)
#endif
#ifdef _OPENMP
      use omp_lib
#endif
      use mod_grid
#ifdef MPI
      integer(kind=4), intent(in) :: nprocs
#endif
      character(len=128), intent(in) :: gfile
      integer(kind=4), intent(in) :: ngrid
      type(data_grids), dimension(ngrid), intent(in) :: dg
      integer(kind=4) :: i, disp, wod, bcf
#ifdef BANKFILE
      integer(kind=4) :: bkf
#endif
      character(len=8) :: myname, pname

      write(6,'(a)') '============================================================'
      write(6,'(a)') '=== Check configurations! [Begin] =========================='
      write(6,'(a)') '============================================================'
! ========================================================================================
      write(6,'(a)') '[Make-time configurations]'
#ifndef CARTESIAN
      write(6,'(a)') '- Coordination: Polar (Default)'
#else
      write(6,'(a)') '- Coordination: Cartesian (-DCARTESIAN)'
#endif
#ifndef UPWIND3
      write(6,'(a)') '- Difference scheme for advection: 1st order upwind (Default)'
#else
      write(6,'(a)') '- Difference scheme for advection: 3rd order upwind (-DUPWIND3)'
#endif
#ifdef REAL_DBLE
      write(6,'(a)') '- Precision: Double precision except file I/O (-DREAL_DBLE)'
#else
#ifdef DBLE_MATH
      write(6,'(a)') '- Precision: Double precision only math functions (-DDBLE_MATH)'
#else
      write(6,'(a)') '- Precision: Single precision'
#endif
#endif
#ifndef MPI
      write(6,'(a)') '- MPI parallelization: OFF'
#else
      write(6,'(a)') '- MPI parallelization: ON (-DMPI)'
#ifdef ONEFILE
      write(6,'(a)') '   - Undecomposed input/output files: ON (-DONEFILE)'
#else
      write(6,'(a)') '   - Undecomposed input/output files: OFF (Default)'
#endif
#ifdef USE_ALLTOALLV
      write(6,'(a)') '   - Inter-domain communication: MPI_Alltoallv (-DUSE_ALLTOALLV)'
#ifdef A2A3D
      write(6,'(a)') '      - A2A3D: ON (-DA2A3D)'
#endif
#ifdef SINGLE_A2A
      write(6,'(a)') '      - SINGLE_A2A: ON (-DSINGLE_A2A)'
#endif
#else
      write(6,'(a)') '   - Inter-domain communication: MPI_Allreduce'
#endif
#endif
#ifdef _OPENMP
      write(6,'(a)') '- OpenMP parallelization: ON'
#else
      write(6,'(a)') '- OpenMP parallelization: OFF'
#endif
#ifndef TIMER
      write(6,'(a)') '- Built-in timer: OFF'
#else
#ifndef TIMER_DETAIL
      write(6,'(a)') '- Built-in timer: ON (-DTIMER)'
#else
      write(6,'(a)') '- Built-in timer: DETAIL (-DTIMER -DTIMER_DETAIL)'
#endif
#endif
#ifndef CONV_CHECK
      write(6,'(a)') '- Convergence check on dispersive: OFF'
#else
      write(6,'(a)') '- Convergence check on dispersive: ON (-DCONV_CHECK)'
#endif
#ifdef DIROUT
      write(6,'(a)') '- Snapshot output: GMT output to directories (-DDIROUT)'
#else
#ifdef NCDIO
      write(6,'(a)') '- Snapshot output: NetCDF (-DNCDIO)'
! === Support start date in time attribute. ====================================
      write(6,'(a,a)') '   - Start date: ', trim(start_date)
! ==============================================================================
#else
      write(6,'(a)') '- Snapshot output: Orignal GMT output'
#endif
#endif
#ifdef PIXELIN
      write(6,'(a)') '- Input files with pixel format: ON (-DPIXELIN)'
#endif
#ifdef PIXELOUT
      write(6,'(a)') '- Output files with pixel format: ON (-DPIXELOUT)'
#endif
! ========================================================================================
      write(6,'(/,a)') '[Run-time configurations]'
! ----------------------------------------------------------------------------------------
      write(6,'(a)') '(Model parameters)'
      write(6,'(a,e15.6)') '- Time step[s] (dt): ', dt
! === tend is [s]. =============================================================
!     write(6,'(a,e15.6)') '- Number of steps[steps] (tend): ', tend
      write(6,'(a,e15.6)') '- End time[s] (tend): ', tend
! ==============================================================================
      write(6,'(a,i9)') '- Snapshot interval[steps] (itmap): ', itmap
      write(6,'(a,i9)') '   - From [step] (itmap_start): ', itmap_start
      write(6,'(a,i9)') '   - To [step] (itmap_end):     ', itmap_end
      write(6,'(a,e15.6)') '- Time step between ruptures[s] (tau): ', tau
      write(6,'(a)') '- Non-dimensional coeff. (positive val.) or Manning''s roughtness coeff. (negative val)'
      write(6,'(a,e15.6)') '   - For sea (cf):   ', cf
      write(6,'(a,e15.6)') '   - For land (cfl): ', cfl
! === Limiter with max Froude number. ==========================================
      write(6,'(a,e15.6)') '- Max. Froude number (froude_lim): ', froude_lim
! ==============================================================================
      write(6,'(a,i3)') '- Coriolis force (coriolis=1:ON/0:OFF): ', coriolis
      write(6,'(a,i3)') '- Deformation bathymetry (def_bathy=1:ON/0:OFF): ', def_bathy
      write(6,'(a,i3)') '- Smooth edges (smooth_edges=1:ON/0:OFF): ', smooth_edges
      write(6,'(a,i3)') '- All grids are copiedy to coarse (c2p_all=1:ON/0:OFF): ', c2p_all
! === 1-way nest ===============================================================
      write(6,'(a,i3)') '- Only p2c interpolation is performed. (nest_1way=1:ON/0:OFF): ', nest_1way
! ==============================================================================
      write(6,'(a,i3)') '- Multiple ruptures (multrupt=1:ON/0:OFF): ', multrupt
! === Initial displacement of child domains is given by interpolation. =========
      write(6,'(a,i3)') '- Initial disp. with interpolation (init_disp_interpolation=1:ON/0:OFF): ', &
         init_disp_interpolation
      if(use_linear == 1) then
         write(6,'(a)')  '   - Linear interpolation is utilized.'
      else
         write(6,'(a)')  '   - 3rd order spline interpolation is utilized.'
      end if
! ==============================================================================
! === for MRI ==================================================================
      write(6,'(a,i3)') '- Initial disp. with Gaussian (init_disp_gaussian=1:ON/0:OFF): ', &
         init_disp_gaussian
! ==============================================================================
! === SINWAVE ==================================================================
      write(6,'(a,i3)') '- Initial disp. with sin wave (init_disp_sinwave=1:ON/0:OFF): ', &
         init_disp_sinwave
! ==============================================================================
! === Displacement =============================================================
      write(6,'(a,i3)') '- Initial disp. with fault calc. (init_disp_fault=1:ON/0:OFF): ', &
         init_disp_fault
      if(init_disp_fault == 1) then
         write(6,'(a,a)')  '   - Fault parameter file (fault_param_file): ', &
            trim(fault_param_file)
         write(6,'(a,i3)') '   - Adopt horizontal disp. effect (hzdisp_effect=1:ON/0:OFF): ', &
            hzdisp_effect
! === Specify lower limit of depth to adopt horizontal displacement effect. ====
         write(6,'(a,e15.6)') '      - Lower limit of depth to apply it[m] (min_depth_hde): ', min_depth_hde
! ==============================================================================
      end if
      write(6,'(a,i3)') '- Apply Kajiura filter. (apply_kj_filter=1:ON/0:OFF): ', &
         apply_kj_filter
! ==============================================================================
#ifndef CARTESIAN
! === Elastic Loading ==========================================================
      write(6,'(a,i3)') '- Elastic loading (with_elastic_loading=1:ON/0:OFF): ', with_elastic_loading
      if(with_elastic_loading == 1) then
         write(6,'(a,e15.6)') '   - Radius in which the loading effect will be estimated (m_radius): ', m_radius
         write(6,'(a,a)') '   - NetCDF file to specify Green function (m_pyfile): ', trim(m_pyfile)
! === Elastic loading with interpolation =======================================
         write(6,'(a,i3)') '   - Elastic loading with interpolation (elastic_loading_interpolation=1:ON/0:OFF): ', &
            elastic_loading_interpolation
! ==============================================================================
      end if
! ==============================================================================
! === Density ==================================================================
      write(6,'(a,i3)') '- Density (with_density=1:ON/0:OFF): ', with_density
      if(with_density == 1) then
         write(6,'(a,e15.6)') '   - The value of rho. (m_rho): ', m_rho
         write(6,'(a,e15.6)') '   - The value of K. (m_K): ', m_K
      end if
! ==============================================================================
#endif
! === Arrival time =============================================================
      write(6,'(a,i3)') '- Check arrival time (check_arrival_time=1:ON/0:OFF): ', check_arrival_time
      if(check_arrival_time == 1) then
         write(6,'(a,e15.6)') '   - Threshold height[m] (check_arrival_height): ', check_arrival_height
      end if
! ==============================================================================
! ----------------------------------------------------------------------------------------
      write(6,'(/,a)') '(Domains)'
      write(6,'(a,a)') '- Grid file: ', trim(gfile)
      write(6,'(a,i6)') '- Number of domains: ', ngrid
#ifndef BANKFILE
      write(6,'(/,a)') '   ID  Name     Parent  Linear  disp_file wod_file bcf_file'
#else
      write(6,'(/,a)') '   ID  Name     Parent  Linear  disp_file wod_file bcf_file bank_file'
#endif
      do i = 1, ngrid
         myname = ''
         myname = trim(dg(i)%my%base_name)
         pname = ''
         pname = trim(dg(i)%parent%base_name)
         disp = 1
         if(trim(dg(i)%my%disp_file) == 'NO_DISPLACEMENT_FILE_GIVEN') disp = 0
         wod = 1
         if(trim(dg(i)%wod_file) == 'NO_WETORDRY_FILE_GIVEN') wod = 0
         bcf = 1
         if(trim(dg(i)%bcf_file) == 'NO_FRICTION_FILE_GIVEN') bcf = 0
#ifndef BANKFILE
         write(6,'(a,i3,a,a,a,a,i6,i11,2i9)') '   ', i, ' ', myname, ' ', pname, dg(i)%my%linear_flag, disp, wod, bcf
#else
         bkf = 1
         if(trim(dg(i)%bank_file) == 'NO_BANK_FILE_GIVEN') bkf = 0
         write(6,'(a,i3,a,a,a,a,i6,i11,2i9,i10)') '   ', i, ' ', myname, ' ', pname, dg(i)%my%linear_flag, disp, wod, bcf, bkf
#endif
      end do
! ----------------------------------------------------------------------------------------
      write(6,'(/,a)') '(Dispersive)'
      write(6,'(a,i3)') '- Dispersive (with_disp=0:OFF/1:ON/2:Except root dom.): ', with_disp
      if(with_disp == 1 .or. with_disp == 2) then
         write(6,'(a,i6)') '   - Max. steps (max_step): ', max_step
#ifdef CONV_CHECK
         write(6,'(a,e15.6)') '   - Truncation error[m/s] (conv_val): ', conv_val
#endif
         write(6,'(a,e15.6)') '   - Min. depth despersive is available[m] (min_depth): ', min_depth
      end if
! ----------------------------------------------------------------------------------------
#ifdef BANKFILE
      write(6,'(/,a)') '(Linedata)'
      if(broken_rate < 0.0d0) then
         write(6,'(a)') '- Wall break: OFF (broken_rate=[negative number])'
      else
         write(6,'(a)') '- Wall break: ON'
         write(6,'(a,e15.6)') '   - Broken rate (broken_rate=[real number]): ', broken_rate
      end if
#endif
! ----------------------------------------------------------------------------------------
      write(6,'(/,a)') '(Absorbing boundary condition)'
      write(6,'(a,i3)') '- Absorbing boundary condition (with_abc=1:ON/0:OFF): ', with_abc
      if(with_abc == 1) then
         write(6,'(a,i6)') '   - Num. of target grids on East/West boundary (nxa):   ', nxa
         write(6,'(a,i6)') '   - Num. of target grids on North/South boundary (nya): ', nya
         write(6,'(a,e15.6)') '   - Absorbing parameter (apara): ', apara
      end if
! ----------------------------------------------------------------------------------------
      write(6,'(/,a)') '(Restart)'
      write(6,'(a,i9)') '- Restart step [step] (restart=0:OFF): ', restart
      write(6,'(a,i9)') '- Restart file interval [steps] (restartinterval=0:OFF): ', restart_interval
! ----------------------------------------------------------------------------------------
! === Support truncation =======================================================
      if(max_time_i /= 0) then
         write(6,'(/,a)') '(Truncation)'
         write(6,'(a,a,a,i0,a)') '- Max time: ', trim(max_time), ' (', max_time_i, ' seconds)'
      end if
! ==============================================================================
#ifdef MPI
      write(6,'(/,a)') '(MPI parallelization)'
      write(6,'(a,i7)') '- Number of processes: ', nprocs
      write(6,'(a,i7)') '   - East-West direction (procx):   ', npx
      write(6,'(a,i7)') '   - North-South direction (procy): ', npy
#endif
! ----------------------------------------------------------------------------------------
#ifdef _OPENMP
      write(6,'(/,a)') '(OpenMP parallelization)'
      write(6,'(a,i7)') '- Number of threads: ', omp_get_max_threads()
#endif
! ========================================================================================
      write(6,'(a)') '============================================================'
      write(6,'(a)') '=== Check configurations! [End] ============================'
      write(6,'(a)') '============================================================'
      return
   end subroutine putpar

end module mod_params
