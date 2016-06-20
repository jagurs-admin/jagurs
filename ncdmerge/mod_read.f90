module mod_read
implicit none
include 'netcdf.inc'
integer(kind=4), private :: stat, varid
! === Speed output. ============================================================
!integer(kind=4), private :: ncid, idid, mhid, mvid, hzid, vxid, vyid
integer(kind=4), private :: ncid, idid, mhid, mvid, hzid, vxid, vyid, speedid
! ==============================================================================
integer(kind=4), private :: lonid, latid, timeid, stepid
integer(kind=4), dimension(3), private :: start, count
real(kind=4) :: missing_value
logical :: flag_missing_value

integer(kind=4) :: timesteps           ! Number of time steps
integer(kind=4) :: totalNx, totalNy    ! Total domain size
integer(kind=4) :: nprocs              ! Number of processes, my rank
integer(kind=4) :: nprocs_x, nprocs_y  ! Number of processes for each direction
logical :: vel_exists = .true.         ! Velocity data exists or not
! === Speed output. ============================================================
logical :: speed_exists = .true.       ! Speed data exists or not
! ==============================================================================
! === Multiple rupture =========================================================
logical :: multrupt = .true.
integer(kind=4) :: nrupt = 1
! ==============================================================================
character(len=64) :: time_att

contains

   subroutine read_params(fname)
      character(len=256), intent(in) :: fname
      integer(kind=4), dimension(2) :: rbuf

      stat = nf_open(fname, NF_NOWRITE, ncid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      ! Read number of time steps.
      stat = nf_inq_dimid(ncid, 'time', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_dimlen(ncid, varid, timesteps)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      ! Read attribute of time step.
      stat = nf_inq_varid(ncid, 'time', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_att_text(ncid, varid, 'units', time_att)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      ! Read total domain size.
      stat = nf_inq_varid(ncid, 'total_xy', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_var_int(ncid, varid, rbuf)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      totalNx = rbuf(1)
      totalNy = rbuf(2)

      ! Read number of processes.
      stat = nf_inq_varid(ncid, 'nprocs', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_var_int(ncid, varid, rbuf)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      nprocs = rbuf(1)

      ! Read Number of processes for each direction
      stat = nf_inq_varid(ncid, 'nprocs_xy', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_var_int(ncid, varid, rbuf)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      nprocs_x = rbuf(1)
      nprocs_y = rbuf(2)

      ! Check whether velocity data exists or not.
      stat = nf_inq_varid(ncid, 'velocity_x', varid)
      if(stat == NF_ENOTVAR) then
         vel_exists = .false.
      end if
      if(stat /= NF_NOERR .and. stat /= NF_ENOTVAR) write(0,'(a)') nf_strerror(stat)
! === Speed output. ============================================================
      ! Check whether speed data exists or not.
      stat = nf_inq_varid(ncid, 'speed', varid)
      if(stat == NF_ENOTVAR) then
         speed_exists = .false.
      end if
      if(stat /= NF_NOERR .and. stat /= NF_ENOTVAR) write(0,'(a)') nf_strerror(stat)
! ==============================================================================

      stat = nf_inq_varid(ncid, 'max_height', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_att_real(ncid, varid, '_FillValue', missing_value)
      if(stat == NF_ENOTATT) then
         flag_missing_value = .false.
      else
         flag_missing_value = .true.
      end if
! === Multiple rupture =========================================================
      stat = nf_inq_dimid(ncid, 'nrupt', varid)
      if(stat /= NF_NOERR) then
         multrupt = .false.
      end if
      if(multrupt) then
         stat = nf_inq_dimlen(ncid, varid, nrupt)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! ==============================================================================
      stat = nf_close(ncid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      write(6,'(a)')     '===== Parameters of global domain ====='
      write(6,'(a,2i6)') 'totalNx, totalNy: ', totalNx, totalNy
      write(6,'(a,i6)')  'timesteps:        ', timesteps
      write(6,'(a,i6)')  'nprocs:           ', nprocs
      write(6,'(a,2i6)') 'nprocs_xy:        ', nprocs_x, nprocs_y
      write(6,'(a,l)')   'vel_exists:       ', vel_exists
! === Speed output. ============================================================
      write(6,'(a,l)')   'speed_exists:     ', speed_exists
! ==============================================================================
! === Multiple rupture =========================================================
      write(6,'(a,l)')   'multrupt:         ', multrupt
! ==============================================================================
      write(6,'(a)')     '======================================='

      return
   end subroutine read_params

   subroutine read_header(fname, lon, lat, rx, ry)
      character(len=256), intent(in) :: fname
      integer(kind=4), intent(out) :: lon, lat, rx, ry
      integer(kind=4), dimension(2) :: rbuf

      stat = nf_open(fname, NF_NOWRITE, ncid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      ! Read local lon size.
      stat = nf_inq_dimid(ncid, 'lon', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_dimlen(ncid, varid, lon)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      ! Read local lat size.
      stat = nf_inq_dimid(ncid, 'lat', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_dimlen(ncid, varid, lat)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
! === Multiple rupture =========================================================
      if(multrupt) then
         stat = nf_inq_dimid(ncid, 'nrupt', varid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         stat = nf_inq_dimlen(ncid, varid, nrupt)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! ==============================================================================

      ! Read myrank in each direction.
      stat = nf_inq_varid(ncid, 'myrank_xy', varid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_var_int(ncid, varid, rbuf)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      rx = rbuf(1)
      ry = rbuf(2)

      ! Get id for each variables.
      stat = nf_inq_varid(ncid, 'lon', lonid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_varid(ncid, 'lat', latid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_varid(ncid, 'time', timeid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_varid(ncid, 'step', stepid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      stat = nf_inq_varid(ncid, 'initial_displacement', idid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_varid(ncid, 'max_height', mhid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_inq_varid(ncid, 'max_velocity', mvid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      stat = nf_inq_varid(ncid, 'wave_height', hzid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      if(vel_exists) then
         stat = nf_inq_varid(ncid, 'velocity_x', vxid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         stat = nf_inq_varid(ncid, 'velocity_y', vyid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! === Speed output. ============================================================
      if(speed_exists) then
         stat = nf_inq_varid(ncid, 'speed', speedid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! ==============================================================================

      ! Set parameters used on nf_get_vara_real.
      start(1) = 1
      start(2) = 1
      start(3) = 0
      count(1) = lon
      count(2) = lat
      count(3) = 1

      return
   end subroutine read_header

! === Speed output. ============================================================
!  subroutine read_snapshot(hz, vx, vy, t)
!     real(kind=4), dimension(:,:), intent(out) :: hz, vx, vy
   subroutine read_snapshot(hz, vx, vy, speed, t)
      real(kind=4), dimension(:,:), intent(out) :: hz, vx, vy, speed
! ==============================================================================
      integer(kind=4), intent(in) :: t

      start(3) = t

      stat = nf_get_vara_real(ncid, hzid, start, count, hz)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      if(vel_exists) then
         stat = nf_get_vara_real(ncid, vxid, start, count, vx)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

         stat = nf_get_vara_real(ncid, vyid, start, count, vy)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! === Speed output. ============================================================
      if(speed_exists) then
         stat = nf_get_vara_real(ncid, speedid, start, count, speed)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! ==============================================================================

      return
   end subroutine read_snapshot

   subroutine read_initial_displacement(zz)
! === Multiple rupture =========================================================
!     real(kind=4), dimension(:,:), intent(out) :: zz
      real(kind=4), dimension(:,:,:), intent(out) :: zz
      integer(kind=4), dimension(3) :: start_id, count_id
      start_id = start
      start_id(3) = 1
      count_id = count
      count_id(3) = nrupt
! ==============================================================================

! === Multiple rupture =========================================================
!     stat = nf_get_vara_real(ncid, idid, start, count, zz)
      stat = nf_get_vara_real(ncid, idid, start_id, count_id, zz)
! ==============================================================================
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      return
   end subroutine read_initial_displacement

   subroutine read_max_hight(hzmax)
      real(kind=4), dimension(:,:), intent(out) :: hzmax

      stat = nf_get_vara_real(ncid, mhid, start, count, hzmax)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      return
   end subroutine read_max_hight

   subroutine read_max_velocity(vmax)
      real(kind=4), dimension(:,:), intent(out) :: vmax

      stat = nf_get_vara_real(ncid, mvid, start, count, vmax)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      return
   end subroutine read_max_velocity

   subroutine read_lonlat(lon, lat, time, step)
      real(kind=8), dimension(:), intent(out) :: lon, lat, time
      integer(kind=4), dimension(:), intent(out) :: step

      stat = nf_get_var_double(ncid, lonid, lon)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_var_double(ncid, latid, lat)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_vara_double(ncid, timeid, 1, timesteps, time)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_get_vara_int(ncid, stepid, 1, timesteps, step)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      return
   end subroutine read_lonlat

   subroutine close_read_file()
      stat = nf_close(ncid)
      if(stat /= NF_NOERR) write(0, '(a)') nf_strerror(stat)
      return
   end subroutine close_read_file

end module mod_read

