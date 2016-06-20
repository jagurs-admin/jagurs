module mod_write
use mod_read, only : totalNx, totalNy, vel_exists, timesteps, time_att, &
! === Speed output. ============================================================
!                    missing_value, flag_missing_value
                     missing_value, flag_missing_value, speed_exists
! ==============================================================================
! === Multiple rupture =========================================================
use mod_read, only : multrupt, nrupt
! ==============================================================================
implicit none
include 'netcdf.inc'
integer(kind=4) :: netcdf4 = 0 ! NetCDF-4 Fortmat for > 4GiB records.

integer(kind=4), private :: stat
! === Speed output. ============================================================
!integer(kind=4), private :: ncid, idid, mhid, mvid, hzid, vxid, vyid
integer(kind=4), private :: ncid, idid, mhid, mvid, hzid, vxid, vyid, speedid
! ==============================================================================
integer(kind=4), private :: lonid, latid, timeid, stepid
! === Multiple rupture =========================================================
integer(kind=4), private :: nruptid
! ==============================================================================
integer(kind=4), dimension(3), private :: start, count

contains

! === Partial merge support ====================================================
!  subroutine open_write_file(fname, only_step)
   subroutine open_write_file(fname, only_step, osize_x, osize_y)
! ==============================================================================
      character(len=256), intent(in) :: fname
      integer(kind=4), intent(in) :: only_step
! === Partial merge support ====================================================
      integer(kind=4), intent(in) :: osize_x, osize_y
! ==============================================================================

      integer(kind=4) :: xid, yid, tid
      integer(kind=4), dimension(3) :: vdims
! === Multiple rupture =========================================================
      integer(kind=4), dimension(3) :: vdims_id
! ==============================================================================
      character(len=128) :: att

      if(netcdf4 == 1) then
         stat = nf_create(trim(fname), ior(NF_CLOBBER, NF_NETCDF4), ncid)
      else
         stat = nf_create(trim(fname), ior(NF_CLOBBER, NF_64BIT_OFFSET), ncid)
      end if
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

! === Partial merge support ====================================================
!     stat = nf_def_dim(ncid, 'lon', totalNx, xid)
      stat = nf_def_dim(ncid, 'lon', osize_x, xid)
! ==============================================================================
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
! === Partial merge support ====================================================
!     stat = nf_def_dim(ncid, 'lat', totalNy, yid)
      stat = nf_def_dim(ncid, 'lat', osize_y, yid)
! ==============================================================================
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
! === Multiple rupture =========================================================
      if(multrupt) then
         stat = nf_def_dim(ncid, 'nrupt', nrupt, nruptid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! ==============================================================================
      stat = nf_def_dim(ncid, 'time', NF_UNLIMITED, tid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      vdims(1) = xid
      vdims(2) = yid
      vdims(3) = tid
! === Multiple rupture =========================================================
      vdims_id(1) = xid
      vdims_id(2) = yid
      vdims_id(3) = nruptid
! ==============================================================================

      stat = nf_def_var(ncid, 'lon', NF_DOUBLE, 1, vdims(1), lonid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Longitude'
      stat = nf_put_att_text(ncid, lonid, 'long_name', len(trim(att)), att)
      att = 'Degrees'
      stat = nf_put_att_text(ncid, lonid, 'units', len(trim(att)), att)

      stat = nf_def_var(ncid, 'lat', NF_DOUBLE, 1, vdims(2), latid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      att = 'Latitude'
      stat = nf_put_att_text(ncid, latid, 'long_name', len(trim(att)), att)
      att = 'Degrees'
      stat = nf_put_att_text(ncid, latid, 'units', len(trim(att)), att)

      if(only_step >= 0) then
         stat = nf_def_var(ncid, 'time', NF_DOUBLE, 1, vdims(3), timeid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         att = 'Time'
         stat = nf_put_att_text(ncid, timeid, 'long_name', len(trim(att)), att)
         att = trim(time_att)
         stat = nf_put_att_text(ncid, timeid, 'units', len(trim(att)), att)

         stat = nf_def_var(ncid, 'step', NF_INT, 1, vdims(3), stepid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         att = 'Model step'
         stat = nf_put_att_text(ncid, stepid, 'long_name', len(trim(att)), att)
         att = 'Steps'
         stat = nf_put_att_text(ncid, stepid, 'units', len(trim(att)), att)
      end if

      if(only_step <= 0) then
         if(only_step == 0 .or. only_step == -1) then
! === Multiple rupture =========================================================
            if(multrupt) then
               stat = nf_def_var(ncid, 'initial_displacement', NF_REAL, 3, vdims_id, idid)
            else
! ==============================================================================
            stat = nf_def_var(ncid, 'initial_displacement', NF_REAL, 2, vdims, idid)
! === Multiple rupture =========================================================
            end if
! ==============================================================================
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
            att = 'Initial displacement'
            stat = nf_put_att_text(ncid, idid, 'long_name', len(trim(att)), att)
            att = 'Meters'
            stat = nf_put_att_text(ncid, idid, 'units', len(trim(att)), att)
         end if
         if(only_step == 0 .or. only_step == -2) then
            stat = nf_def_var(ncid, 'max_height', NF_REAL, 2, vdims, mhid)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
            att = 'Maximum wave height'
            stat = nf_put_att_text(ncid, mhid, 'long_name', len(trim(att)), att)
            att = 'Meters'
            stat = nf_put_att_text(ncid, mhid, 'units', len(trim(att)), att)
            if(flag_missing_value) then
                stat = nf_put_att_real(ncid, mhid, '_FillValue', NF_REAL, 1, &
                                       missing_value)
            end if
         end if
         if(only_step == 0 .or. only_step == -3) then
            stat = nf_def_var(ncid, 'max_velocity', NF_REAL, 2, vdims, mvid)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
            att = 'Maximum flow speed'
            stat = nf_put_att_text(ncid, mvid, 'long_name', len(trim(att)), att)
            att = 'Meters/Second'
            stat = nf_put_att_text(ncid, mvid, 'units', len(trim(att)), att)
         end if
      end if

      if(only_step >= 0) then
         stat = nf_def_var(ncid, 'wave_height', NF_REAL, 3, vdims, hzid)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         att = 'Wave height'
         stat = nf_put_att_text(ncid, hzid, 'long_name', len(trim(att)), att)
         att = 'Meters'
         stat = nf_put_att_text(ncid, hzid, 'units', len(trim(att)), att)

         if(vel_exists) then
            stat = nf_def_var(ncid, 'velocity_x', NF_REAL, 3, vdims, vxid)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
            att = 'Flow velocity (X-component)'
            stat = nf_put_att_text(ncid, vxid, 'long_name', len(trim(att)), att)
            att = 'Meters/Second'
            stat = nf_put_att_text(ncid, vxid, 'units', len(trim(att)), att)

            stat = nf_def_var(ncid, 'velocity_y', NF_REAL, 3, vdims, vyid)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
            att = 'Flow velocity (Y-component)'
            stat = nf_put_att_text(ncid, vyid, 'long_name', len(trim(att)), att)
            att = 'Meters/Second'
            stat = nf_put_att_text(ncid, vyid, 'units', len(trim(att)), att)
         end if
! === Speed output. ============================================================
         if(speed_exists) then
            stat = nf_def_var(ncid, 'speed', NF_REAL, 3, vdims, speedid)
            if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
            att = 'Flow speed'
            stat = nf_put_att_text(ncid, speedid, 'long_name', len(trim(att)), att)
            att = 'Meters/Second'
            stat = nf_put_att_text(ncid, speedid, 'units', len(trim(att)), att)
         end if
! ==============================================================================
      end if

      stat = nf_enddef(ncid)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      start(1) = 1
      start(2) = 1
      start(3) = 0
! === Partial merge support ====================================================
!     count(1) = totalNx
!     count(2) = totalNy
      count(1) = osize_x
      count(2) = osize_y
! ==============================================================================
      count(3) = 1

      return
   end subroutine open_write_file

! === Speed output. ============================================================
!  subroutine write_snapshot(hz, vx, vy)
!     real(kind=4), dimension(:,:), intent(in) :: hz, vx, vy
   subroutine write_snapshot(hz, vx, vy, speed)
      real(kind=4), dimension(:,:), intent(in) :: hz, vx, vy, speed
! ==============================================================================

      start(3) = start(3) + 1

      stat = nf_put_vara_real(ncid, hzid, start, count, hz)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      if(vel_exists) then
         stat = nf_put_vara_real(ncid, vxid, start, count, vx)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

         stat = nf_put_vara_real(ncid, vyid, start, count, vy)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! === Speed output. ============================================================
      if(speed_exists) then
         stat = nf_put_vara_real(ncid, speedid, start, count, speed)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if
! ==============================================================================

      return
   end subroutine write_snapshot

   subroutine write_initial_displacement(zz)
! === Multiple rupture =========================================================
!     real(kind=4), dimension(:,:), intent(in) :: zz
      real(kind=4), dimension(:,:,:), intent(in) :: zz
      integer(kind=4), dimension(3) :: start_id, count_id
      start_id = start
      start_id(3) = 1
      count_id = count
      count_id(3) = nrupt
! ==============================================================================

! === Multiple rupture =========================================================
!     stat = nf_put_vara_real(ncid, idid, start, count, zz)
      stat = nf_put_vara_real(ncid, idid, start_id, count_id, zz)
! ==============================================================================
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      return
   end subroutine write_initial_displacement

   subroutine write_max_hight(hzmax)
      real(kind=4), dimension(:,:), intent(in) :: hzmax

      stat = nf_put_vara_real(ncid, mhid, start, count, hzmax)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      return
   end subroutine write_max_hight

   subroutine write_max_velocity(vmax)
      real(kind=4), dimension(:,:), intent(in) :: vmax

      stat = nf_put_vara_real(ncid, mvid, start, count, vmax)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)

      return
   end subroutine write_max_velocity

   subroutine write_lonlat(lon, lat, time, step, only_step)
      real(kind=8), dimension(:), intent(in) :: lon, lat, time
      integer(kind=4), dimension(:), intent(in) :: step
      integer(kind=4), intent(in) :: only_step

      stat = nf_put_var_double(ncid, lonid, lon)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      stat = nf_put_var_double(ncid, latid, lat)
      if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      if(only_step >= 0) then
         stat = nf_put_vara_double(ncid, timeid, 1, timesteps, time)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
         stat = nf_put_vara_int(ncid, stepid, 1, timesteps, step)
         if(stat /= NF_NOERR) write(0,'(a)') nf_strerror(stat)
      end if

      return
   end subroutine write_lonlat

   subroutine close_write_file()
      stat = nf_close(ncid)
      if(stat /= NF_NOERR) write(0, '(a)') nf_strerror(stat)
      return
   end subroutine close_write_file

end module mod_write

