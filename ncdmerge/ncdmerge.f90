program ncdmerge
use mod_read
use mod_write
implicit none
character(len=256) :: basename, fname
real(kind=4), allocatable, dimension(:,:) :: rhz, rvx, rvy ! Read buffers
real(kind=4), allocatable, dimension(:,:) :: whz, wvx, wvy ! Write buffers
! === Multiple rupture =========================================================
real(kind=4), allocatable, dimension(:,:,:) :: rid ! Read buffers
real(kind=4), allocatable, dimension(:,:,:) :: wid ! Write buffers
! ==============================================================================
! === Speed output. ============================================================
real(kind=4), allocatable, dimension(:,:) :: rspeed ! Read buffers
real(kind=4), allocatable, dimension(:,:) :: wspeed ! Write buffers
! ==============================================================================
real(kind=8), allocatable, dimension(:) :: rlon, rlat ! Read buffers
real(kind=8), allocatable, dimension(:) :: wlon, wlat ! Write buffers
real(kind=8), allocatable, dimension(:) :: time       ! Read/Write buffer
integer(kind=4), allocatable, dimension(:) :: step    ! Read/Write buffer
integer(kind=4) :: p, lon, lat, ist, jst, rx, ry, nbx, nby, t
integer(kind=4) :: i, j
! === Multiple rupture =========================================================
integer(kind=4) :: irupt
! ==============================================================================
integer(kind=4) :: only_step = 0, index_step
! === Partial merge support ====================================================
!namelist /params/ basename, netcdf4, only_step
integer(kind=4) :: root_proc = -1, xnum = -1, ynum = -1
integer(kind=4) :: root_x, root_y, merged_procs, k
integer(kind=4), allocatable, dimension(:) :: proc_list
integer(kind=4) :: osize_x, osize_y, ostart_x, ostart_y, oend_x, oend_y
namelist /params/ basename, netcdf4, only_step, root_proc, xnum, ynum
! ==============================================================================

read(5,params)
write(6,'(a,a)') 'Target file basename: ', trim(basename)
if(only_step /= 0) then
   if(only_step == -1) then
      write(6,'(a)') 'NOTE: Only initial_displacement is merged!'
   else if(only_step == -2) then
      write(6,'(a)') 'NOTE: Only max_height is merged!'
   else if(only_step == -3) then
      write(6,'(a)') 'NOTE: Only max_velocity is merged!'
   else
      write(6,'(a,i9,a)') 'NOTE: Only time step ', only_step, ' is merged!'
   end if
end if
if(netcdf4 == 1) then
   write(6,'(a)') 'NOTE: NetCDF format is NetCDF4!'
else
   write(6,'(a)') 'NOTE: NetCDF format is classic/64-bit offset!'
   write(6,'(a)') '      If target domains are very large (> 4GiB/record) and error occurs,'
   write(6,'(a)') '      try NetCDF4 Format! (Specify "netcdf4=1" in namelist.)'
end if

! First, read parameters used to merge from [basename].000000.nc.
write(fname,'(a,a,i6.6,a)') trim(basename), '.', 0, '.nc'
call read_params(fname)
! === Partial merge support ====================================================
if(root_proc /= -1) then
   if(root_proc < 0 .or. root_proc > nprocs-1) then
     write(0,'(a,i9,a)') 'ERROR: root_proc should be 0 ~ ', nprocs-1, '!'
     stop
   end if
   root_x = mod(root_proc, nprocs_x)
   root_y = root_proc/nprocs_x
   if((xnum == -1) .or. (root_x+xnum > nprocs_x)) xnum = nprocs_x - root_x
   if((ynum == -1) .or. (root_y+ynum > nprocs_y)) ynum = nprocs_y - root_y
   merged_procs = xnum*ynum
   allocate(proc_list(merged_procs))
   k = 0
   do j = root_y, root_y+ynum-1
      do i = root_x, root_x+xnum-1
          k = k + 1
          proc_list(k) = i + j*nprocs_x
      end do
   end do
   write(6,'(a)') '====== Partial merge specifications ==='
   write(0,'(a,i8)')        'Root process rank:  ', root_proc
   write(0,'(a,i6,a,i6,a)') '   (x-rank,y-rank): (', root_x, ',', root_y, ')'
   write(0,'(a,i6,a,i6)')   'Target size:         ', xnum, 'x', ynum
   write(0,'(a)')           'Target list:'
   write(0,'(4i8)') proc_list
   write(6,'(a)') '======================================='
else
   merged_procs = nprocs
   allocate(proc_list(merged_procs))
   do i = 0, nprocs-1
      proc_list(i+1) = i
   end do
end if
! ==============================================================================

! Calc. local domain size.
nbx = totalNx / nprocs_x
if(mod(totalNx, nprocs_x) /= 0) nbx = nbx + 1
nby = totalNy / nprocs_y
if(mod(totalNy, nprocs_y) /= 0) nby = nby + 1

! === Partial merge support ====================================================
if(root_proc /= -1) then
   osize_x = 0
   do i = root_x, root_x+xnum-1
      if(i == nprocs_x-1) then
         osize_x = osize_x + (totalNx-nbx*(nprocs_x-1))
      else
         osize_x = osize_x + nbx
      end if
   end do
   osize_y = 0
   do j = root_y, root_y+ynum-1
      if(j == nprocs_y-1) then
         osize_y = osize_y + (totalNy-nby*(nprocs_y-1))
      else
         osize_y = osize_y + nby
      end if
   end do
   ostart_x = nbx*root_x + 1
!  ostart_y = nby*root_y + 1
   ostart_y = max(totalNy-((root_y+ynum-1)+1)*nby+1, 1)
   oend_x = ostart_x + osize_x - 1
   oend_y = ostart_y + osize_y - 1
else
   osize_x = totalNx
   osize_y = totalNy
   ostart_x = 1
   ostart_y = 1
   oend_x = totalNx
   oend_y = totalNy
end if
! ==============================================================================

! Open write files.
if(only_step == 0) then
   write(fname,'(a,a)') trim(basename), '.nc'
else
   if(only_step == -1) then
      write(fname,'(a,a)') trim(basename), '.initial_displacement.nc'
   else if(only_step == -2) then
      write(fname,'(a,a)') trim(basename), '.max_height.nc'
   else if(only_step == -3) then
      write(fname,'(a,a)') trim(basename), '.max_velocity.nc'
   else
      write(fname,'(a,a,i8.8,a)') trim(basename), '.', only_step, '.nc'
   end if
end if
write(6,'(a,a)') 'Target file name: ', trim(fname)
! === Partial merge support ====================================================
!call open_write_file(fname, only_step)
call open_write_file(fname, only_step, osize_x, osize_y)
! ==============================================================================

! Allocate write buffers.
! === Partial merge support ====================================================
!allocate(whz(totalNx,totalNy))
!allocate(wvx(totalNx,totalNy))
!allocate(wvy(totalNx,totalNy))
!allocate(wlon(totalNx))
!allocate(wlat(totalNy))
allocate(whz(ostart_x:oend_x,ostart_y:oend_y))
! === Multiple rupture =========================================================
allocate(wid(ostart_x:oend_x,ostart_y:oend_y,nrupt))
! ==============================================================================
allocate(wvx(ostart_x:oend_x,ostart_y:oend_y))
allocate(wvy(ostart_x:oend_x,ostart_y:oend_y))
! === Speed output. ============================================================
allocate(wspeed(ostart_x:oend_x,ostart_y:oend_y))
! ==============================================================================
allocate(wlon(ostart_x:oend_x))
allocate(wlat(ostart_y:oend_y))
! ==============================================================================
allocate(time(timesteps))
allocate(step(timesteps))

! Read/Write initial_displacement, max_height and wave_height.
! === Partial merge support ====================================================
!do p = 0, nprocs-1
do k = 1, merged_procs
   p = proc_list(k)
! ==============================================================================
   ! Open local file and read parameters.
   write(fname,'(a,a,i6.6,a)') trim(basename), '.', p, '.nc'
   call read_header(fname, lon, lat, rx, ry)

   ! Allocate read buffers.
   if(only_step <= 0) then
! === Multiple rupture =========================================================
!     allocate(rhz(lon,lat))
      allocate(rid(lon,lat,nrupt))
! ==============================================================================
      allocate(rvx(lon,lat))
      allocate(rvy(lon,lat))
   end if
   allocate(rlon(lon))
   allocate(rlat(lat))

   ! Read local variables.
   if(only_step <= 0) then
! === Multiple rupture =========================================================
!     call read_initial_displacement(rhz)
      call read_initial_displacement(rid)
! ==============================================================================
      call read_max_hight(rvx)
      call read_max_velocity(rvy)
   end if
   call read_lonlat(rlon,rlat,time,step)

   ! Local vals. are copied into global write buffers.
   ist = rx*nbx + 1
   jst = max(totalNy-(ry+1)*nby+1, 1)
   if(only_step <= 0) then
      do j = 1, lat
         do i = 1, lon
! === Multiple rupture =========================================================
!           whz(i+ist-1,j+jst-1) = rhz(i,j)
! ==============================================================================
            wvx(i+ist-1,j+jst-1) = rvx(i,j)
            wvy(i+ist-1,j+jst-1) = rvy(i,j)
         end do
      end do
! === Multiple rupture =========================================================
      do irupt = 1, nrupt
         do j = 1, lat
            do i = 1, lon
               wid(i+ist-1,j+jst-1,irupt) = rid(i,j,irupt)
            end do
         end do
      end do
! ==============================================================================
   end if
   do i = 1, lon
      wlon(i+ist-1) = rlon(i)
   end do
   do j = 1, lat
      wlat(j+jst-1) = rlat(j)
   end do

   ! Deallocate read buffers.
   if(only_step <= 0) then
! === Multiple rupture =========================================================
!     deallocate(rhz)
      deallocate(rid)
! ==============================================================================
      deallocate(rvx)
      deallocate(rvy)
   end if
   deallocate(rlon)
   deallocate(rlat)

   ! Close local file.
   call close_read_file()
end do
! Write global variables.
if(only_step <= 0) then
! === Multiple rupture =========================================================
!  if(only_step == 0 .or. only_step == -1) call write_initial_displacement(whz)
   if(only_step == 0 .or. only_step == -1) call write_initial_displacement(wid)
! ==============================================================================
   if(only_step == 0 .or. only_step == -2) call write_max_hight(wvx)
   if(only_step == 0 .or. only_step == -3) call write_max_velocity(wvy)
end if
if(only_step <= -1) then
   timesteps = 0
else
   if(only_step /= 0) then
      do t = 1, timesteps
         if(step(t) == only_step) exit
      end do
      if(t > timesteps) then
         write(0,'(a,i9,a)') 'ERROR: Specified time step ', only_step, ' does NOT exist!'
         stop
      end if
      timesteps = 1
      index_step = t
      time(1) = time(index_step)
      step(1) = step(index_step)
   end if
end if
call write_lonlat(wlon,wlat,time,step,only_step)
deallocate(wlon)
deallocate(wlat)
deallocate(time)
deallocate(step)

if(only_step <= 0) then
   if(only_step == -1) then
      write(6,'(a)') 'initial_displacement has been merged!!!'
   else if(only_step == -2) then
      write(6,'(a)') 'max_height has been merged!!!'
   else if(only_step == -3) then
      write(6,'(a)') 'max_velocity has been merged!!!'
   else
      write(6,'(a)') 'initial_displacement, max_height and max_velocity have been merged!!!'
   end if
end if

! Read/Write snapshots in each time step.
do t = 1, timesteps
! === Partial merge support ====================================================
!  do p = 0, nprocs-1
   do k = 1, merged_procs
      p = proc_list(k)
! ==============================================================================
      ! Open local file and read parameters.
      write(fname,'(a,a,i6.6,a)') trim(basename), '.', p, '.nc'
      call read_header(fname, lon, lat, rx, ry)

      ! Allocate read buffers.
      allocate(rhz(lon,lat))
      allocate(rvx(lon,lat))
      allocate(rvy(lon,lat))
! === Speed output. ============================================================
      allocate(rspeed(lon,lat))
! ==============================================================================

      ! Read local variables.
      if(only_step == 0) then
! === Speed output. ============================================================
!        call read_snapshot(rhz, rvx, rvy, t)
         call read_snapshot(rhz, rvx, rvy, rspeed, t)
! ==============================================================================
      else
! === Speed output. ============================================================
!        call read_snapshot(rhz, rvx, rvy, index_step)
         call read_snapshot(rhz, rvx, rvy, rspeed, index_step)
! ==============================================================================
      end if

      ! Local vals. are copied into global write buffers.
      ist = rx*nbx + 1
      jst = max(totalNy-(ry+1)*nby+1, 1)
      do j = 1, lat
         do i = 1, lon
            whz(i+ist-1,j+jst-1) = rhz(i,j)
            wvx(i+ist-1,j+jst-1) = rvx(i,j)
            wvy(i+ist-1,j+jst-1) = rvy(i,j)
! === Speed output. ============================================================
            wspeed(i+ist-1,j+jst-1) = rspeed(i,j)
! ==============================================================================
         end do
      end do

      ! Deallocate read buffers.
      deallocate(rhz)
      deallocate(rvx)
      deallocate(rvy)
! === Speed output. ============================================================
      deallocate(rspeed)
! ==============================================================================

      ! Close local file.
      call close_read_file()
   end do
   ! Write global variables.
! === Speed output. ============================================================
!  call write_snapshot(whz, wvx, wvy)
   call write_snapshot(whz, wvx, wvy, wspeed)
! ==============================================================================
   write(6,'(a,i8,a,i8)') 'Snapshot has been merged!!! ', t, '/', timesteps
end do

! Finally, close write file.
call close_write_file()

! Deallocate write buffers.
deallocate(whz)
! === Multiple rupture =========================================================
deallocate(wid)
! ==============================================================================
deallocate(wvx)
deallocate(wvy)
! === Speed output. ============================================================
deallocate(wspeed)
! ==============================================================================
! === Partial merge support ====================================================
deallocate(proc_list)
! ==============================================================================

write(6,'(a)') 'DONE!!!'

stop
end program ncdmerge
