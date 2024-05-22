program printp
   implicit none
   real(kind=8) :: srclon, srclat, tmax, dt
   real(kind=8), allocatable, dimension(:,:) :: nm_atm
   integer(kind=4) :: dist_s, dist_e
   real(kind=8) :: distlon, distlat, t, dist, Samp, SDt, OH
   integer(kind=4) :: nstep, i, j, it, nt, nm_ind, sind, minmode, maxmode
   logical :: DCflag, SYNflag
   character(len=256) :: atmfile
   namelist /normalmode/ srclon, srclat, dist_s, dist_e, tmax, dt,sind, minmode, maxmode, Samp, SDt, DCflag, SYNflag, OH, atmfile
   namelist /params/ distlon, distlat, dist_s, dist_e, tmax, dt

   open(1,file='normalmode.namelist',action='read',status='old',form='formatted')
   read(1,normalmode)
   close(1)

   write(6,'(a)') '- Normalmode parameters in "normalmode.namelist"'
   write(6,'(a,2f15.3)') 'srclon, srclat: ', srclon, srclat
   write(6,'(a,2i12)') 'dist_s, dist_e:   ', dist_s, dist_e
   write(6,'(a,2f15.3)') 'tmax, dt:       ', tmax, dt

   nstep = nint(tmax/dt) + 1

   allocate(nm_atm(dist_s:dist_e,nstep))
   write(6,'(/,a)') '- Reading "normalmode.dat"'
   open(1,file='normalmode.dat',action='read',status='old',form='unformatted')
   read(1) nm_atm
   close(1)

   open(1,file='printp.namelist',action='read',status='old',form='formatted')
   read(1,params)
   close(1)

   write(6,'(/a)') '- Observation point is specified in "printp.namelist"'
   write(6,'(a,2f15.3)') 'distlon, distlat: ', distlon, distlat

   call deg2dist(distlat, distlon, srclat, srclon, dist)
   nm_ind = nint(dist/1000.0d0)
   write(6,'(a,i8)') 'nm_ind: ', nm_ind

   write(6,'(/a)')  '------------------------------'
   write(6,'(a)') 'Time[s]         P'
   write(6,'(a)')  '------------------------------'
   nt = nint(tmax/dt)
   do it = 0, nt
      t = it*dt
      write(6,'(f15.3,e15.3)') t, nm_atm(nm_ind,it+1)
   end do
   write(6,'(a)')  '------------------------------'

   stop
end program printp

subroutine deg2dist(lat1, lon1, lat2, lon2, dist)
   real(kind=8), intent(in) :: lat1, lon1, lat2, lon2
   real(kind=8), intent(out) :: dist
   real(kind=8), parameter :: R = 6371.0d3
   real(kind=8), parameter :: pi = 3.14159265d0
   real(kind=8) :: phi1, phi2, lamda1, lamda2

   phi1 = lat1*pi/180.0d0
   phi2 = lat2*pi/180.0d0
   lamda1 = lon1*pi/180.0d0
   lamda2 = lon2*pi/180.0d0

   dist = R*acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lamda1 - lamda2))

   return
end subroutine deg2dist
