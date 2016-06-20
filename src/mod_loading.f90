#ifndef MPI
#define FASTER
#endif

#include "real.h"
module mod_loading
#ifdef _OPENMP
use omp_lib
#endif
use mod_grid
use mod_params, only : m_radius, m_pyfile
! === For ensemble =============================================================
#ifdef MPI
#ifdef MULTI
use mod_multi, only : MPI_MEMBER_WORLD
#endif
#endif
! ==============================================================================
implicit none
integer(kind=4), private :: Nrec
real(kind=8), private, allocatable, dimension(:) :: Dist, Val

! Parameters.
real(kind=8), private, parameter :: rote = 6371.0d3
real(kind=8), private, parameter :: M_PI = 3.14159265d0

#ifdef MPI
integer(kind=4), private :: nprocs, myrank, npx, npy, rankx, ranky
integer(kind=4), private :: MPI_X_WORLD, MPI_Y_WORLD
#endif

private :: LinearInterpolate

contains

#ifdef  MPI
   subroutine loading_mpi_initialize(nprocs_in, myrank_in, npx_in, npy_in, rankx_in, ranky_in)
      use mpi
      integer(kind=4), intent(in)  :: nprocs_in, myrank_in, npx_in, npy_in, rankx_in, ranky_in
      integer(kind=4) :: color, key, ierr
      nprocs = nprocs_in
      myrank = myrank_in
      npx = npx_in
      npy = npy_in
      rankx = rankx_in
      ranky = ranky_in

      color = ranky
      key = rankx
#ifndef MULTI
      call MPI_comm_split(MPI_COMM_WORLD, color, key, MPI_X_WORLD, ierr)
#else
      call MPI_comm_split(MPI_MEMBER_WORLD, color, key, MPI_X_WORLD, ierr)
#endif

      color = rankx
      key = ranky
#ifndef MULTI
      call MPI_comm_split(MPI_COMM_WORLD, color, key, MPI_Y_WORLD, ierr)
#else
      call MPI_comm_split(MPI_MEMBER_WORLD, color, key, MPI_Y_WORLD, ierr)
#endif

      return
   end subroutine loading_mpi_initialize
#endif

   function LinearInterpolate(y1, y2, mu)
      real(kind=8), intent(in)  :: y1, y2, mu
      real(kind=8) :: LinearInterpolate
      LinearInterpolate =  y1*(1.0d0 - mu) + y2*mu
      return
   end function LinearInterpolate

#ifndef MULTI
   subroutine loading_getval()
#else
   subroutine loading_getval(input_dirname)
#endif
      include 'netcdf.inc'
#ifdef MULTI
      character(len=128), intent(in) :: input_dirname
      character(len=256) :: file_multi
#endif
      integer(kind=4) :: retval, ncid, NrecID, varid, i

      ! Read File!
      ! Open the file. NC_NOWRITE tells netCDF we want read-only access to the file.
#ifndef MULTI
      retval = nf_open(m_pyfile, NF_NOWRITE, ncid)
#else
      file_multi = trim(input_dirname) // trim(m_pyfile)
      retval = nf_open(file_multi, NF_NOWRITE, ncid)
#endif
      if(retval /= NF_NOERR) then
         write(0,'(a,i,a,a)') 'netcdf err=', retval, ' in file=', trim(m_pyfile)
         stop
      end if

      retval = nf_inq_dimid(ncid, 'Nrec', NrecID)
      if(retval /= NF_NOERR) then
         write(0,'(a,i)') 'netcdf err=', retval
         stop
      end if

      retval = nf_inq_dimlen(ncid, NrecID, Nrec)
      if(retval /= NF_NOERR) then
         write(0,'(a,i)') 'netcdf err=', retval
         stop
      end if

      allocate(Dist(Nrec))
      allocate(Val(Nrec))

      ! Get the varid of the data variable, based on its name.
      retval = nf_inq_varid(ncid, 'dist', varid)
      if(retval /= NF_NOERR) then
         write(0,'(a,i)') 'netcdf err=', retval
         stop
      end if

      ! Read the data.
      retval = nf_get_var_double(ncid, varid, Dist)
      if(retval /= NF_NOERR) then
         write(0,'(a,i)') 'netcdf err=', retval
         stop
      end if

      ! Get the varid of the data variable, based on its name.
      retval = nf_inq_varid(ncid, 'Gz', varid)
      if(retval /= NF_NOERR) then
         write(0,'(a,i)') 'netcdf err=', retval
         stop
      end if

      ! Read the data.
      retval = nf_get_var_double(ncid, varid, Val)
      if(retval /= NF_NOERR) then
         write(0,'(a,i)') 'netcdf err=', retval
         stop
      end if

      retval = nf_close(ncid)

      write(6,'(a)') '[loading] ReadNC file'
      write(6,'(a,i6)') '[loading] Nrec: ', Nrec
      write(6,'(a)') '[loading] i      Dist           Val'
      do i = 1, Nrec
         write(6,'(a,i6,2e15.6)') '[loading] ', i, Dist(i), Val(i)
      end do

      return
   end subroutine loading_getval

   subroutine loading_initialize(dg)
#ifndef __SX__
      include 'fftw3.f'
#endif
      type(data_grids), target, intent(inout) :: dg 
      
      ! Local work arrays.
      real(kind=8), allocatable, dimension(:,:) :: greenZZ      ! Global work
      real(kind=8), allocatable, dimension(:,:) :: greenZdouble ! Local work
      complex(kind=8), allocatable, dimension(:,:) :: green_in_Z

      ! Local vals.
      integer(kind=4) :: m_nxg, m_nyg, N, ix, iy, num, iradius
#ifndef __SX__
      integer(kind=8), allocatable, dimension(:) :: xgreen_planZ, ygreen_planZ ! FFTW3 plan
#else
      integer(kind=4), allocatable, dimension(:,:) :: ifax
      real(kind=8), allocatable, dimension(:,:) :: trigs
      complex(kind=8), allocatable, dimension(:,:) :: work
#endif
      real(kind=8) :: m_dx_d, m_dy_d, m_dx, m_dy, &
         xval, yval, xvalm, yvalm, res, sum, Delta, Delta_t, r, &
         xval_t, yval_t, val_, maxi

      ! Pointers to type loading in each domain.
      integer(kind=4), pointer :: N_X, N_Y
#ifndef __SX__
      integer(kind=8), pointer, dimension(:) :: &
         xplan_forward, yplan_forward, xplan_backward, yplan_backward ! FFTW3 plan
#else
      integer(kind=4), pointer, dimension(:,:) :: ifax_x, ifax_y
      real(kind=8), pointer, dimension(:,:) :: trigs_x, trigs_y
#ifndef REAL_FFT
      complex(kind=8), pointer, dimension(:,:) :: work_x, work_y
#else
      real(kind=8), pointer, dimension(:,:) :: work_x
      complex(kind=8), pointer, dimension(:,:) :: work_y
#endif
#endif
      complex(kind=8), pointer, dimension(:,:) :: green_out_Z
      complex(kind=8), pointer, dimension(:,:) :: xfftbuf, yfftbuf
#ifdef REAL_FFT
      real(kind=8), pointer, dimension(:,:) :: realbuf
#endif
      real(kind=8), pointer, dimension(:,:) :: defZmap1

      integer(kind=4) :: i, j
#ifdef MPI
      integer(kind=4), pointer :: nx0, ny0, nx1, ny1, nx2, ny2, ibias, jbias
      real(kind=8), pointer, dimension(:) :: sendbuf1, recvbuf1
      complex(kind=8), pointer, dimension(:) :: sendbuf2, recvbuf2
      integer(kind=4), pointer, dimension(:) :: sendcounts1, sdispls1, recvcounts1, rdispls1
      integer(kind=4), pointer, dimension(:) :: sendcounts2, sdispls2, recvcounts2, rdispls2

      complex(kind=8), allocatable, dimension(:) :: sendbufg, recvbufg
      integer(kind=4), allocatable, dimension(:) :: sendcountsg, sdisplsg, recvcountsg, rdisplsg

      integer(kind=4) :: ixend, iyend, kx, kxend, ky, kyend, nyg
      integer(kind=4) :: p, j_, yst, yen, ylen, ind, ierr
#endif
! === OpenMP ===================================================================
      integer(kind=4), allocatable, dimension(:) :: jtsg, jteg, jtng ! N_Y, nyg
      integer(kind=4), pointer, dimension(:) :: its2, ite2, itn2     ! N_X, nx2
      integer(kind=4), pointer, dimension(:) :: jts1, jte1, jtn1     ! N_Y, ny1
      integer(kind=4) :: it, nthreads
! ==============================================================================
#ifdef __SX__
#ifndef MPI
      integer(kind=4) :: ierr
#endif
#endif
#ifndef MPI
      m_nxg = dg%my%nx
      m_nyg = dg%my%ny
#else
      m_nxg = dg%my%totalNx
      m_nyg = dg%my%totalNy
#endif
      m_dx_d = dg%my%dh     ! X space step in degree
      m_dy_d = dg%my%dh     ! Y space step in degree
      m_dx = dg%my%dth*rote ! Average X space step
      m_dy = dg%my%dth*rote ! Average Y space step

      N_X => dg%loading%N_X
      N_Y => dg%loading%N_Y

      N_X = m_nxg
      N_Y = m_nyg
      if(m_nxg > m_nyg) then
         N_X = m_nxg
         N_Y = N_X
      else
         N_Y = m_nyg
         N_X = N_Y
      end if
      N_X = 2**ceiling(dlog(dble(N_X))/dlog(2.0d0))
      N_Y = 2**ceiling(dlog(dble(N_Y))/dlog(2.0d0))
      N = N_X*N_Y

      write(6,'(a,2i6)') '[loading] m_nxg, m_nyg: ', m_nxg, m_nyg
      write(6,'(a,3i8)') '[loading] N_X, N_Y, N: ', N_X, N_Y, N
#ifdef MPI
      ix = dg%my%ix
      iy = dg%my%iy
      ixend = dg%my%ixend
      iyend = dg%my%iyend

      kx = dg%my%kx
      ky = dg%my%ky
      kxend = dg%my%kxend
      kyend = dg%my%kyend

      nx0 => dg%loading%nx0
      ny0 => dg%loading%ny0
      nx1 => dg%loading%nx1
      ny1 => dg%loading%ny1
      nx2 => dg%loading%nx2
      ny2 => dg%loading%ny2
      ibias => dg%loading%ibias
      jbias => dg%loading%jbias

      ! Space grid without edges.
      nx0 = ixend - ix + 1 ! xsize.
      ny0 = iyend - iy + 1 ! ysize.
      ibias = ix - kx ! x bias to skip edges.
      jbias = iy - ky ! y bias to skip edges.

      write(6,'(a,2i6)') '[loading] nprocs, myrank: ', nprocs, myrank
      write(6,'(a,2i6)') '[loading] ix, ixend: ', ix, ixend
      write(6,'(a,2i6)') '[loading] iy, iyend: ', iy, iyend
      write(6,'(a,2i6)') '[loading] nx0, ny0: ', nx0, ny0

      ! Space grid for x-FFT.
      nx1 = dg%my%totalNx
      ny1 = ny0/npx
      if(rankx < mod(ny0,npx)) ny1 = ny1 + 1

      write(6,'(a,2i6)') '[loading] nx1, ny1: ', nx1, ny1

      ! MPI_Alltoallv1: (nx0,ny0) -> (nx1,ny1)
      allocate(dg%loading%sendbuf1(nx0*ny0))
      allocate(dg%loading%sendcounts1(0:npx-1))
      allocate(dg%loading%sdispls1(0:npx-1))
      sendbuf1    => dg%loading%sendbuf1
      sendcounts1 => dg%loading%sendcounts1
      sdispls1    => dg%loading%sdispls1

      allocate(dg%loading%recvbuf1(nx1*ny1))
      allocate(dg%loading%recvcounts1(0:npx-1))
      allocate(dg%loading%rdispls1(0:npx-1))
      recvbuf1    => dg%loading%recvbuf1
      recvcounts1 => dg%loading%recvcounts1
      rdispls1    => dg%loading%rdispls1

      do p = 0, npx-1
         j = ny0/npx
         if(p < mod(ny0,npx)) j = j + 1
         sendcounts1(p) = nx0*j
      end do
      call MPI_Alltoall(sendcounts1, 1, MPI_INTEGER, &
                        recvcounts1, 1, MPI_INTEGER, MPI_X_WORLD, ierr)
      sdispls1(0) = 0
      rdispls1(0) = 0
      do p = 1, npx-1
         sdispls1(p) = sdispls1(p-1) + sendcounts1(p-1)
         rdispls1(p) = rdispls1(p-1) + recvcounts1(p-1)
      end do

      write(6,'(a)') '[loading] p, sendcounts1(p), sdispls1(p), recvcounts1(p), rdispls1(p)'
      do p = 0, npx-1
         write(6,'(5i8)') p, sendcounts1(p), sdispls1(p), recvcounts1(p), rdispls1(p)
      end do

      ! Grid for y-FFT.
#ifndef REAL_FFT
      nx2 = N_X/nprocs
      if(myrank < mod(N_X,nprocs)) nx2 = nx2 + 1
#else
      nx2 = (N_X/2+1)/nprocs
      if(myrank < mod(N_X/2+1,nprocs)) nx2 = nx2 + 1
#endif
      ny2 = N_Y

      write(6,'(a,2i6)') '[loading] nx2, ny2: ', nx2, ny2

      ! MPI_Alltoallv1: (N_X,ny1) -> (N_Y,nx2)
#ifndef REAL_FFT
      allocate(dg%loading%sendbuf2(N_X*ny1))
#else
      allocate(dg%loading%sendbuf2((N_X/2+1)*ny1))
#endif
      allocate(dg%loading%sendcounts2(0:nprocs-1))
      allocate(dg%loading%sdispls2(0:nprocs-1))
      sendbuf2    => dg%loading%sendbuf2
      sendcounts2 => dg%loading%sendcounts2
      sdispls2    => dg%loading%sdispls2

      allocate(dg%loading%recvbuf2(N_Y*nx2))
      allocate(dg%loading%recvcounts2(0:nprocs-1))
      allocate(dg%loading%rdispls2(0:nprocs-1))
      recvbuf2    => dg%loading%recvbuf2
      recvcounts2 => dg%loading%recvcounts2
      rdispls2    => dg%loading%rdispls2

      do p = 0, nprocs-1
#ifndef REAL_FFT
         i = N_X/nprocs
         if(p < mod(N_X,nprocs)) i = i + 1
#else
         i = (N_X/2+1)/nprocs
         if(p < mod(N_X/2+1,nprocs)) i = i + 1
#endif
         sendcounts2(p) = i*ny1
      end do
#ifndef MULTI
      call MPI_Alltoall(sendcounts2, 1, MPI_INTEGER, &
                        recvcounts2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoall(sendcounts2, 1, MPI_INTEGER, &
                        recvcounts2, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      sdispls2(0) = 0
      rdispls2(0) = 0
      do p = 1, nprocs-1
         sdispls2(p) = sdispls2(p-1) + sendcounts2(p-1)
         rdispls2(p) = rdispls2(p-1) + recvcounts2(p-1)
      end do

      write(6,'(a)') '[loading] p, sendcounts2(p), sdispls2(p), recvcounts2(p), rdispls2(p)'
      do p = 0, nprocs-1
         write(6,'(5i8)') p, sendcounts2(p), sdispls2(p), recvcounts2(p), rdispls2(p)
      end do

! === Used only for calc. green_out_Z_tmp. =====================================
      ! Space grid for x-FFT.
      nyg = N_Y/nprocs
      if(myrank < mod(N_Y,nprocs)) nyg = nyg + 1

      ! MPI_Alltoallv1: (N_X,nyg) -> (N_Y,nx2)
#ifndef REAL_FFT
      allocate(sendbufg(N_X*nyg))
#else
      allocate(sendbufg((N_X/2+1)*nyg))
#endif
      allocate(sendcountsg(0:nprocs-1))
      allocate(sdisplsg(0:nprocs-1))

      allocate(recvbufg(N_Y*nx2))
      allocate(recvcountsg(0:nprocs-1))
      allocate(rdisplsg(0:nprocs-1))

      do p = 0, nprocs-1
#ifndef REAL_FFT
         i = N_X/nprocs
         if(p < mod(N_X,nprocs)) i = i + 1
#else
         i = (N_X/2+1)/nprocs
         if(p < mod(N_X/2+1,nprocs)) i = i + 1
#endif
         sendcountsg(p) = i*nyg
      end do
#ifndef MULTI
      call MPI_Alltoall(sendcountsg, 1, MPI_INTEGER, &
                        recvcountsg, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoall(sendcountsg, 1, MPI_INTEGER, &
                        recvcountsg, 1, MPI_INTEGER, MPI_MEMBER_WORLD, ierr)
#endif
      sdisplsg(0) = 0
      rdisplsg(0) = 0
      do p = 1, nprocs-1
         sdisplsg(p) = sdisplsg(p-1) + sendcountsg(p-1)
         rdisplsg(p) = rdisplsg(p-1) + recvcountsg(p-1)
      end do

      write(6,'(a)') '[loading] p, sendcountsg(p), sdisplsg(p), recvcountsg(p), rdisplsg(p)'
      do p = 0, nprocs-1
         write(6,'(5i8)') p, sendcountsg(p), sdisplsg(p), recvcountsg(p), rdisplsg(p)
      end do
! ==============================================================================
#endif
! === OpenMP ===================================================================
#ifdef _OPENMP
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif
      allocate(jtsg(0:nthreads-1))
      allocate(jteg(0:nthreads-1))
      allocate(jtng(0:nthreads-1))

#ifndef __SX__
      allocate(dg%loading%xplan_forward(0:nthreads-1))
      allocate(dg%loading%yplan_forward(0:nthreads-1))
      allocate(dg%loading%xplan_backward(0:nthreads-1))
      allocate(dg%loading%yplan_backward(0:nthreads-1))
#else
      allocate(dg%loading%ifax_x(20,0:nthreads-1))
      allocate(dg%loading%ifax_y(20,0:nthreads-1))
      allocate(dg%loading%trigs_x(N_X*2,0:nthreads-1))
      allocate(dg%loading%trigs_y(N_Y*2,0:nthreads-1))
#ifndef REAL_FFT
#ifndef MPI
      allocate(dg%loading%work_x(N_X,N_Y))
      allocate(dg%loading%work_y(N_Y,N_X))
#else
      allocate(dg%loading%work_x(N_X,ny1))
      allocate(dg%loading%work_y(N_Y,nx2))
#endif
#else
#ifndef MPI
      allocate(dg%loading%work_x(N_X+2,N_Y))
      allocate(dg%loading%work_y(N_Y,N_X/2+1))
#else
      allocate(dg%loading%work_x(N_X+2,ny1))
      allocate(dg%loading%work_y(N_Y,nx2))
#endif
#endif
#endif
      allocate(dg%loading%its2(0:nthreads-1))
      allocate(dg%loading%ite2(0:nthreads-1))
      allocate(dg%loading%itn2(0:nthreads-1))
      allocate(dg%loading%jts1(0:nthreads-1))
      allocate(dg%loading%jte1(0:nthreads-1))
      allocate(dg%loading%jtn1(0:nthreads-1))

#ifndef __SX__
      xplan_forward  => dg%loading%xplan_forward
      yplan_forward  => dg%loading%yplan_forward
      xplan_backward => dg%loading%xplan_backward
      yplan_backward => dg%loading%yplan_backward
#else
      ifax_x => dg%loading%ifax_x
      ifax_y => dg%loading%ifax_y
      trigs_x => dg%loading%trigs_x
      trigs_y => dg%loading%trigs_y
      work_x => dg%loading%work_x
      work_y => dg%loading%work_y
#endif
      its2 => dg%loading%its2
      ite2 => dg%loading%ite2
      itn2 => dg%loading%itn2
      jts1 => dg%loading%jts1
      jte1 => dg%loading%jte1
      jtn1 => dg%loading%jtn1
#ifndef MPI
      do it = 0, nthreads-1
         jtng(it) = N_Y/nthreads
         if(it < mod(N_Y,nthreads)) jtng(it) = jtng(it) + 1
#ifndef REAL_FFT
         itn2(it) = N_X/nthreads
         if(it < mod(N_X,nthreads)) itn2(it) = itn2(it) + 1
#else
         itn2(it) = (N_X/2+1)/nthreads
         if(it < mod(N_X/2+1,nthreads)) itn2(it) = itn2(it) + 1
#endif
         jtn1(it) = N_Y/nthreads
         if(it < mod(N_Y,nthreads)) jtn1(it) = jtn1(it) + 1
      end do
#else
      do it = 0, nthreads-1
         jtng(it) = nyg/nthreads
         if(it < mod(nyg,nthreads)) jtng(it) = jtng(it) + 1
         itn2(it) = nx2/nthreads
         if(it < mod(nx2,nthreads)) itn2(it) = itn2(it) + 1
         jtn1(it) = ny1/nthreads
         if(it < mod(ny1,nthreads)) jtn1(it) = jtn1(it) + 1
      end do
#endif
      jtsg(0) = 1
      its2(0) = 1
      jts1(0) = 1
      do it = 1, nthreads-1
         jtsg(it) = jtsg(it-1) + jtng(it-1)
         its2(it) = its2(it-1) + itn2(it-1)
         jts1(it) = jts1(it-1) + jtn1(it-1)
      end do
      do it = 0, nthreads-1
         jteg(it) = jtsg(it) + jtng(it) - 1
         ite2(it) = its2(it) + itn2(it) - 1
         jte1(it) = jts1(it) + jtn1(it) - 1
      end do
! ==============================================================================

#ifndef MPI
      allocate(greenZZ(N_X,N_Y))
#ifndef __SX__
      allocate(greenZdouble(N_X,N_Y))
#else
#ifndef REAL_FFT
      allocate(greenZdouble(N_X,N_Y))
#else
      allocate(greenZdouble(N_X+2,N_Y))
#endif
#endif
#else
      allocate(greenZZ(N_X,nyg))
#ifndef __SX__
      allocate(greenZdouble(N_X,nyg))
#else
#ifndef REAL_FFT
      allocate(greenZdouble(N_X,nyg))
#else
      allocate(greenZdouble(N_X+2,nyg))
#endif
#endif

      yst = N_Y/nprocs
      yst = yst*myrank + min(myrank, mod(N_Y,nprocs)) + 1
      yen = yst + nyg - 1
#endif

      maxi = 0.0d0
!$omp parallel
#ifndef FASTER
#ifndef MPI
!$omp do private(ix, xval, yval, xvalm, yvalm, res, Delta, r, sum, num, xval_t, yval_t, &
!$omp            Delta_t, val_, iradius) reduction(max:maxi)
      do iy = 1, N_Y
#else
!$omp do private(ix, xval, yval, xvalm, yvalm, res, Delta, r, sum, num, xval_t, yval_t, &
!$omp            Delta_t, val_, iradius) reduction(max:maxi)
      do iy = yst, yen
#endif
         do ix = 1, N_X
            if((ix <= N_X/2+1) .and. (iy <= N_Y/2+1)) then
               xval  = (      ix - 1)*m_dx_d
               yval  = (      iy - 1)*m_dy_d
               xvalm = (      ix - 1)*m_dx
               yvalm = (      iy - 1)*m_dy
            else if((ix <= N_X/2+1) .and. (iy > N_Y/2+1)) then
               xval  = (      ix - 1)*m_dx_d
               yval  = (N_Y - iy + 1)*m_dy_d
               xvalm = (      ix - 1)*m_dx
               yvalm = (N_Y - iy + 1)*m_dy
            else if((ix > N_X/2+1) .and. (iy > N_Y/2+1)) then
               xval  = (N_X - ix + 1)*m_dx_d
               yval  = (N_Y - iy + 1)*m_dy_d
               xvalm = (N_X - ix + 1)*m_dx
               yvalm = (N_Y - iy + 1)*m_dy
            else if((ix > N_X/2+1) .and. (iy <= N_Y/2+1)) then
               xval  = (N_X - ix + 1)*m_dx_d
               yval  = (      iy - 1)*m_dy_d
               xvalm = (N_X - ix + 1)*m_dx
               yvalm = (      iy - 1)*m_dy
            end if
#else
!$omp do private(ix, xval, yval, xvalm, yvalm, res, Delta, r, sum, num, xval_t, yval_t, &
!$omp            Delta_t, val_, iradius) reduction(max:maxi)
      do iy = 1, N_Y/2+1
         do ix = 1, N_X/2+1
            xval  = (ix - 1)*m_dx_d
            yval  = (iy - 1)*m_dy_d
            xvalm = (ix - 1)*m_dx
            yvalm = (iy - 1)*m_dy
#endif
            res = 0.0d0
            Delta = dacos(dcos(yval*M_PI/180.0d0)*dcos(dabs(xval*M_PI/180.0d0)))
            r = Delta*rote
            if(r < m_radius*1000.0d0) then
               sum = 0.0d0
               num = 0
               xval_t = xval
               do while(xval_t < xval + m_dx_d)
                  yval_t = yval
                  do while(yval_t < yval + m_dy_d)
                     Delta_t = dabs(dacos(dcos(yval_t*M_PI/180.0d0)*dcos(dabs(xval_t*M_PI/180.0d0))))
                     if(Delta_t < 1.0d-5) Delta_t = 1.0d-5
                     Delta_t = Delta_t*180.0d0/M_PI
                     val_ = 0.0d0
                     if(Delta_t < Dist(1)) then
                        val_ = Val(1)
                     else
                        do iradius = 1, Nrec
                           if(Dist(iradius) < Delta_t .and. Dist(iradius + 1) > Delta_t) exit
                        end do
                        val_ = LinearInterpolate(Val(iradius), Val(iradius + 1), &
                           (Delta_t - Dist(iradius))/(Dist(iradius + 1) - Dist(iradius)))
                     end if
                     sum = sum + val_
                     num = num + 1
                     yval_t = yval_t + m_dy_d/10.0d0
                  end do
                  xval_t = xval_t + m_dx_d/10.0d0
               end do
               res = sum/dble(num)*1020.0d0*m_dx*m_dy
            else
               res = 0.0d0
            end if

            if(res > 1.0d0) then
               write(6,'(a,e15.6,a,e15.6)') '[loading] UNBUG ', Delta, ' ', res
            end if

#ifndef MPI
            greenZZ(ix, iy) = res
#else
            greenZZ(ix, iy-yst+1) = res
#endif

            if(maxi < dabs(res)) then
               maxi = abs(res)
            end if
         end do
      end do
!$omp single
#ifdef MPI
#ifndef MULTI
      call MPI_Allreduce(MPI_IN_PLACE, maxi, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Allreduce(MPI_IN_PLACE, maxi, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_MEMBER_WORLD, ierr)
#endif
#endif

      write(6,'(a,e15.6,a,e15.6,a,e15.6,a)') '[loading] MAX GREEN ', maxi, ' (', m_dx, ', ', m_dy, ')'
!$omp end single
#ifdef FASTER
!$omp do private(ix)
      do iy = 1, N_Y/2+1
         do ix = 1, N_X/2-1
            greenZZ(N_X/2+ix+1,         iy) = greenZZ(N_X/2-ix+1,         iy)
         end do
      end do

!$omp do private(ix)
      do iy = 1, N_Y/2-1
         do ix = 1, N_X/2+1
            greenZZ(        ix, N_Y/2+iy+1) = greenZZ(        ix, N_Y/2-iy+1)
         end do
      end do

!$omp do private(ix)
      do iy = 1, N_Y/2-1
         do ix = 1, N_X/2-1
            greenZZ(N_X/2+ix+1, N_Y/2+iy+1) = greenZZ(N_X/2-ix+1, N_Y/2-iy+1)
         end do
      end do
#endif

!$omp do private(iy, ix)
      do it = 0, nthreads-1
         do iy = jtsg(it), jteg(it)
            do ix = 1, N_X
               greenZdouble(ix,iy) = greenZZ(ix,iy)
            end do
         end do
      end do

!$omp single
#ifndef REAL_FFT
#ifndef MPI
      allocate(green_in_Z(N_X,N_Y))
      allocate(dg%loading%green_out_Z(N_Y,N_X))
#else
      allocate(green_in_Z(N_X,nyg))
      allocate(dg%loading%green_out_Z(N_Y,nx2))
#endif
#else
#ifndef MPI
      allocate(green_in_Z(N_X/2+1,N_Y))
      allocate(dg%loading%green_out_Z(N_Y,N_X/2+1))
#else
      allocate(green_in_Z(N_X/2+1,nyg))
      allocate(dg%loading%green_out_Z(N_Y,nx2))
#endif
#endif

      green_out_Z => dg%loading%green_out_Z
!$omp end single

#ifndef REAL_FFT
!$omp do private(iy, ix)
      do it = 0, nthreads-1
         do iy = jtsg(it), jteg(it)
            do ix = 1, N_X
               green_in_Z(ix, iy) = dcmplx(greenZdouble(ix, iy),0.0d0)
            end do
         end do
      end do
#endif

#ifndef __SX__
!$omp single
      allocate(xgreen_planZ(0:nthreads-1))
      allocate(ygreen_planZ(0:nthreads-1))

      ! dfftw_plan* is NOT thread-safe!
      do it = 0, nthreads-1
#ifndef REAL_FFT
         call dfftw_plan_many_dft(xgreen_planZ(it), 1, N_X, jtng(it), &
                                  green_in_Z(1,jtsg(it)), 0, 1, N_X,  &
                                  green_in_Z(1,jtsg(it)), 0, 1, N_X,  &
                                  FFTW_FORWARD, FFTW_ESTIMATE)
#else
         call dfftw_plan_many_dft_r2c(xgreen_planZ(it), 1, N_X, jtng(it),      &
                                      greenZdouble(1,jtsg(it)), 0, 1, N_X,     &
                                      green_in_Z(1,jtsg(it)),   0, 1, N_X/2+1, &
                                      FFTW_ESTIMATE)
#endif

         call dfftw_plan_many_dft(ygreen_planZ(it), 1, N_Y, itn2(it), &
                                  green_out_Z(1,its2(it)), 0, 1, N_Y, &
                                  green_out_Z(1,its2(it)), 0, 1, N_Y, &
                                  FFTW_FORWARD, FFTW_ESTIMATE)
      end do
!$omp end single
#else
!$omp single
      allocate(ifax(20,0:nthreads-1))
#ifndef REAL_FFT
      allocate(trigs(2*N_X,0:nthreads-1))
#else
      allocate(trigs(N_X,0:nthreads-1))
#endif
#ifndef REAL_FFT
#ifndef MPI
      allocate(work(N_X,N_Y))
#else
      allocate(work(N_X,nyg))
#endif
#else
#ifndef MPI
      allocate(work(N_X/2+1,N_Y))
#else
      allocate(work(N_X/2+1,nyg))
#endif
#endif
!$omp end single
#endif

      ! Forward FFT on X-direction
#ifndef __SX__
#ifndef REAL_FFT
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(xgreen_planZ(it))
      end do
#else
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute_dft_r2c(xgreen_planZ(it),greenZdouble(1,jtsg(it)),green_in_Z(1,jtsg(it)))
      end do
#endif
#else
#ifndef REAL_FFT
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMFB(N_X,jtng(it),green_in_Z(1,jtsg(it)),1,N_X,1,ifax(1,it),trigs(1,it),work(1,jtsg(it)),ierr)
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMFB(N_X,jtng(it),greenZdouble(1,jtsg(it)),1,N_X+2,1,ifax(1,it),trigs(1,it),work(1,jtsg(it)),ierr)
      end do
#endif
!$omp single
      deallocate(ifax)
      deallocate(trigs)
      deallocate(work)
!$omp end single
#endif

#ifndef MPI
      ! Transposiion: (N_X,N_Y) -> (N_Y,N_X)
!$omp do private(i, j)
      do it = 0, nthreads-1
         do i = its2(it), ite2(it)
            do j = 1, N_Y
#ifndef __SX__
               green_out_Z(j,i) = green_in_Z(i,j)
#else
#ifndef REAL_FFT
               green_out_Z(j,i) = green_in_Z(i,j)
#else
               green_out_Z(j,i) = dcmplx(greenZdouble(2*i-1,j),greenZdouble(2*i,j))
#endif
#endif
            end do
         end do
      end do
#else
      ! Transposiion: (N_X,nyg) -> (N_Y,nx2)
!$omp single
      ind = 0
#ifndef REAL_FFT
      do i = 1, N_X
#else
      do i = 1, N_X/2+1
#endif
         do j = 1, nyg
            ind = ind + 1
#ifndef __SX__
            sendbufg(ind) = green_in_Z(i,j)
#else
#ifndef REAL_FFT
            sendbufg(ind) = green_in_Z(i,j)
#else
            sendbufg(ind) = dcmplx(greenZdouble(2*i-1,j),greenZdouble(2*i,j))
#endif
#endif
         end do
      end do
#ifndef MULTI
      call MPI_Alltoallv(sendbufg, sendcountsg, sdisplsg, MPI_DOUBLE_COMPLEX, &
                         recvbufg, recvcountsg, rdisplsg, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoallv(sendbufg, sendcountsg, sdisplsg, MPI_DOUBLE_COMPLEX, &
                         recvbufg, recvcountsg, rdisplsg, MPI_DOUBLE_COMPLEX, MPI_MEMBER_WORLD, ierr)
#endif
      yst = 0
      do p = 0, nprocs-1
         ylen = recvcountsg(p)/nx2
         do i = 1, nx2
            do j = 1, ylen
               ind = rdisplsg(p) + (i-1)*ylen + j
               j_ =  yst + j
               green_out_Z(j_,i) = recvbufg(ind)
            end do
         end do
         yst = yst + ylen
      end do
!$omp end single
#endif

      ! Forward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(ygreen_planZ(it))
      end do
#else
!$omp single
      allocate(ifax(20,0:nthreads-1))
      allocate(trigs(2*N_Y,0:nthreads-1))
#ifndef MPI
      allocate(work(N_Y,N_X))
#else
      allocate(work(N_Y,nx2))
#endif
!$omp end single
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMFB(N_Y,itn2(it),green_out_Z(1,its2(it)),1,N_Y,1,ifax(1,it),trigs(1,it),work(1,its2(it)),ierr)
      end do
#endif

#ifndef __SX__
!$omp single
      do it = 0, nthreads-1
         call dfftw_destroy_plan(xgreen_planZ(it))
         call dfftw_destroy_plan(ygreen_planZ(it))
      end do
#else
!$omp single
      deallocate(ifax)
      deallocate(trigs)
      deallocate(work)
#endif

      deallocate(jtsg)
      deallocate(jteg)
      deallocate(jtng)
#ifndef __SX__
      deallocate(xgreen_planZ)
      deallocate(ygreen_planZ)
#endif

#ifdef MPI
      deallocate(sendbufg)
      deallocate(sendcountsg)
      deallocate(sdisplsg)

      deallocate(recvbufg)
      deallocate(recvcountsg)
      deallocate(rdisplsg)
#endif

#ifndef __SX__
#ifndef REAL_FFT
#ifndef MPI
      allocate(dg%loading%xfftbuf(N_X,N_Y))
      allocate(dg%loading%yfftbuf(N_Y,N_X))
#else
      allocate(dg%loading%xfftbuf(N_X,ny1))
      allocate(dg%loading%yfftbuf(N_Y,nx2))
#endif
#else
#ifndef MPI
      allocate(dg%loading%realbuf(N_X,N_Y))
      allocate(dg%loading%xfftbuf(N_X/2+1,N_Y))
      allocate(dg%loading%yfftbuf(N_Y,N_X/2+1))
#else
      allocate(dg%loading%realbuf(N_X,ny1))
      allocate(dg%loading%xfftbuf(N_X/2+1,ny1))
      allocate(dg%loading%yfftbuf(N_Y,nx2))
#endif
#endif
#else
#ifndef REAL_FFT
#ifndef MPI
      allocate(dg%loading%xfftbuf(N_X+1,N_Y))
      allocate(dg%loading%yfftbuf(N_Y+1,N_X))
#else
      allocate(dg%loading%xfftbuf(N_X+1,ny1))
      allocate(dg%loading%yfftbuf(N_Y+1,nx2))
#endif
#else
#ifndef MPI
      allocate(dg%loading%realbuf(N_X+2,N_Y))
      allocate(dg%loading%xfftbuf(N_X/2+1,N_Y))
      allocate(dg%loading%yfftbuf(N_Y+1,N_X/2+1))
#else
      allocate(dg%loading%realbuf(N_X+2,ny1))
      allocate(dg%loading%xfftbuf(N_X/2+1,ny1))
      allocate(dg%loading%yfftbuf(N_Y+1,nx2))
#endif
#endif
#endif
#ifdef REAL_FFT
      realbuf => dg%loading%realbuf
#endif
      xfftbuf => dg%loading%xfftbuf
      yfftbuf => dg%loading%yfftbuf
!$omp end single

      ! First touch!
#ifndef REAL_FFT
!$omp do private(i, j)
      do it = 0, nthreads-1
         do j = jts1(it), jte1(it)
            do i = 1, N_X
               xfftbuf(i,j) = dcmplx(0.0d0,0.0d0)
            end do
         end do
         do i = its2(it), ite2(it)
            do j = 1, N_Y
               yfftbuf(j,i) = dcmplx(0.0d0,0.0d0)
            end do
         end do
      end do
#else
!$omp do private(i, j)
      do it = 0, nthreads-1
         do j = jts1(it), jte1(it)
            do i = 1, N_X
               realbuf(i,j) = 0.0d0
            end do
         end do
         do j = jts1(it), jte1(it)
            do i = 1, N_X/2+1
               xfftbuf(i,j) = dcmplx(0.0d0,0.0d0)
            end do
         end do
         do i = its2(it), ite2(it)
            do j = 1, N_Y
               yfftbuf(j,i) = dcmplx(0.0d0,0.0d0)
            end do
         end do
      end do
#endif

!$omp single
      ! dfftw_plan* is NOT thread-safe!
      do it = 0, nthreads-1
#ifndef __SX__
#ifndef REAL_FFT
         call dfftw_plan_many_dft(xplan_forward(it), 1, N_X, jtn1(it), &
                                  xfftbuf(1,jts1(it)), 0, 1, N_X,      &
                                  xfftbuf(1,jts1(it)), 0, 1, N_X,      &
                                  FFTW_FORWARD, FFTW_MEASURE)
#else
         call dfftw_plan_many_dft_r2c(xplan_forward(it), 1, N_X, jtn1(it), &
                                      realbuf(1,jts1(it)), 0, 1, N_X,      &
                                      xfftbuf(1,jts1(it)), 0, 1, N_X/2+1,  &
                                      FFTW_MEASURE)
#endif

         call dfftw_plan_many_dft(yplan_forward(it), 1, N_Y, itn2(it), &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,      &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,      &
                                  FFTW_FORWARD, FFTW_MEASURE)

#ifndef REAL_FFT
         call dfftw_plan_many_dft(xplan_backward(it), 1, N_X, jtn1(it), &
                                  xfftbuf(1,jts1(it)), 0, 1, N_X,       &
                                  xfftbuf(1,jts1(it)), 0, 1, N_X,       &
                                  FFTW_BACKWARD, FFTW_MEASURE)
#else
         call dfftw_plan_many_dft_c2r(xplan_backward(it), 1, N_X, jtn1(it), &
                                      xfftbuf(1,jts1(it)), 0, 1, N_X/2+1,   &
                                      realbuf(1,jts1(it)), 0, 1, N_X,       &
                                      FFTW_MEASURE)
#endif

         call dfftw_plan_many_dft(yplan_backward(it), 1, N_Y, itn2(it), &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,       &
                                  yfftbuf(1,its2(it)), 0, 1, N_Y,       &
                                  FFTW_BACKWARD, FFTW_MEASURE)
#else
#ifndef REAL_FFT
         call ZFCMFB(N_X,jtn1(it),xfftbuf(1,jts1(it)),1,N_X+1,0,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
#else
         call DFRMFB(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+2,0,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
#endif
         call ZFCMFB(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+1,0,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
#endif
      end do

      deallocate(greenZZ)
      deallocate(greenZdouble)
      deallocate(green_in_Z)
#ifndef MPI
      allocate(dg%loading%defZmap(m_nxg,m_nyg))
      allocate(dg%loading%defZmap1(m_nxg,m_nyg))
#else
      allocate(dg%loading%defZmap(nx0,ny0))
      allocate(dg%loading%defZmap1(nx0,ny0))
#endif
      defZmap1 => dg%loading%defZmap1
!$omp end single

#ifndef MPI
!$omp do private(i)
      do j = 1, m_nyg
         do i = 1, m_nxg
#else
!$omp do private(i)
      do j = 1, ny0
         do i = 1, nx0
#endif
            defZmap1(i,j) = 0.0d0
         end do
      end do
!$omp end parallel

      return
   end subroutine loading_initialize

   subroutine loading_finalize(dg)
#ifndef __SX__
      include 'fftw3.f'
#endif
      type(data_grids), target, intent(inout) :: dg 
! === OpenMP ===================================================================
      integer(kind=4) :: it, nthreads
! ==============================================================================
#ifdef _OPENMP
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif
      if(allocated(Dist)) deallocate(Dist)
      if(allocated(Val)) deallocate(Val)
      deallocate(dg%loading%green_out_Z)

#ifndef __SX__
      do it = 0, nthreads-1
         call dfftw_destroy_plan(dg%loading%xplan_forward(it))
         call dfftw_destroy_plan(dg%loading%yplan_forward(it))
         call dfftw_destroy_plan(dg%loading%xplan_backward(it))
         call dfftw_destroy_plan(dg%loading%yplan_backward(it))
      end do

      deallocate(dg%loading%xplan_forward)
      deallocate(dg%loading%yplan_forward)
      deallocate(dg%loading%xplan_backward)
      deallocate(dg%loading%yplan_backward)
#else
      deallocate(dg%loading%ifax_x)
      deallocate(dg%loading%ifax_y)
      deallocate(dg%loading%trigs_x)
      deallocate(dg%loading%trigs_y)
      deallocate(dg%loading%work_x)
      deallocate(dg%loading%work_y)
#endif

      deallocate(dg%loading%its2)
      deallocate(dg%loading%ite2)
      deallocate(dg%loading%itn2)
      deallocate(dg%loading%jts1)
      deallocate(dg%loading%jte1)
      deallocate(dg%loading%jtn1)

#ifdef REAL_FFT
      deallocate(dg%loading%realbuf)
#endif
      deallocate(dg%loading%xfftbuf)
      deallocate(dg%loading%yfftbuf)

#ifdef MPI
      deallocate(dg%loading%sendbuf1)
      deallocate(dg%loading%sendcounts1)
      deallocate(dg%loading%sdispls1)

      deallocate(dg%loading%recvbuf1)
      deallocate(dg%loading%recvcounts1)
      deallocate(dg%loading%rdispls1)

      deallocate(dg%loading%sendbuf2)
      deallocate(dg%loading%sendcounts2)
      deallocate(dg%loading%sdispls2)

      deallocate(dg%loading%recvbuf2)
      deallocate(dg%loading%recvcounts2)
      deallocate(dg%loading%rdispls2)
#endif

      deallocate(dg%loading%defZmap)
      deallocate(dg%loading%defZmap1)
      return
   end subroutine loading_finalize

   subroutine loading_run(dg)
#ifndef __SX__
      include 'fftw3.f'
#endif
      type(data_grids), target, intent(inout) :: dg 

      ! Local vals.
      real(kind=8) :: coef_norm
      integer(kind=4) :: i, j

      ! Pointers to vals. in each domain.
#ifndef MPI
      integer(kind=4), pointer :: nlon, nlat
#endif
      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz
      integer(kind=4), pointer, dimension(:,:) :: ifz

      ! Pointers to type loading in each domain.
      integer(kind=4), pointer :: N_X, N_Y
#ifndef __SX__
      integer(kind=8), pointer, dimension(:) :: &
         xplan_forward, yplan_forward, xplan_backward, yplan_backward ! FFTW3 plan
#else
      integer(kind=4), pointer, dimension(:,:) :: ifax_x, ifax_y
      real(kind=8), pointer, dimension(:,:) :: trigs_x, trigs_y
#ifndef REAL_FFT
      complex(kind=8), pointer, dimension(:,:) :: work_x, work_y
#else
      real(kind=8), pointer, dimension(:,:) :: work_x
      complex(kind=8), pointer, dimension(:,:) :: work_y
#endif
#endif
      complex(kind=8), pointer, dimension(:,:) :: green_out_Z
      complex(kind=8), pointer, dimension(:,:) :: xfftbuf, yfftbuf
#ifdef REAL_FFT
      real(kind=8), pointer, dimension(:,:) :: realbuf
#endif
      real(kind=8), pointer, dimension(:,:) :: defZmap, defZmap1
#ifdef MPI
      integer(kind=4), pointer :: nx0, ny0, nx1, ny1, nx2, ny2, ibias, jbias
      real(kind=8), pointer, dimension(:) :: sendbuf1, recvbuf1
      complex(kind=8), pointer, dimension(:) :: sendbuf2, recvbuf2
      integer(kind=4), pointer, dimension(:) :: sendcounts1, sdispls1, recvcounts1, rdispls1
      integer(kind=4), pointer, dimension(:) :: sendcounts2, sdispls2, recvcounts2, rdispls2

      integer(kind=4) :: i_, j_, p, ind, xst, xlen, yst, ylen, ierr
#endif
#ifdef __SX__
#ifndef MPI
      integer(kind=4) :: ierr
#endif
#endif
! === OpenMP ===================================================================
      integer(kind=4), pointer, dimension(:) :: its2, ite2, itn2 ! N_X, nx2
      integer(kind=4), pointer, dimension(:) :: jts1, jte1, jtn1 ! N_Y, ny1
      integer(kind=4) :: it, nthreads
! ==============================================================================

#ifndef MPI
      nlon => dg%my%nx
      nlat => dg%my%ny

      N_X => dg%loading%N_X
      N_Y => dg%loading%N_Y
#else
      nx0 => dg%loading%nx0
      ny0 => dg%loading%ny0
      nx1 => dg%loading%nx1
      ny1 => dg%loading%ny1
      nx2 => dg%loading%nx2
      ny2 => dg%loading%ny2
      ibias => dg%loading%ibias
      jbias => dg%loading%jbias

      N_X => dg%loading%N_X
      N_Y => dg%loading%N_Y

      sendbuf1    => dg%loading%sendbuf1
      sendcounts1 => dg%loading%sendcounts1
      sdispls1    => dg%loading%sdispls1

      recvbuf1    => dg%loading%recvbuf1
      recvcounts1 => dg%loading%recvcounts1
      rdispls1    => dg%loading%rdispls1

      sendbuf2    => dg%loading%sendbuf2
      sendcounts2 => dg%loading%sendcounts2
      sdispls2    => dg%loading%sdispls2

      recvbuf2    => dg%loading%recvbuf2
      recvcounts2 => dg%loading%recvcounts2
      rdispls2    => dg%loading%rdispls2
#endif

      coef_norm = 1.0d0/(dble(N_X)*dble(N_Y))

      hz => dg%wave_field%hz
      ifz => dg%wod_flags

#ifndef __SX__
      xplan_forward  => dg%loading%xplan_forward
      yplan_forward  => dg%loading%yplan_forward
      xplan_backward => dg%loading%xplan_backward
      yplan_backward => dg%loading%yplan_backward
#else
      ifax_x => dg%loading%ifax_x
      ifax_y => dg%loading%ifax_y
      trigs_x => dg%loading%trigs_x
      trigs_y => dg%loading%trigs_y
      work_x => dg%loading%work_x
      work_y => dg%loading%work_y
#endif
! === OpenMP ===================================================================
      its2 => dg%loading%its2
      ite2 => dg%loading%ite2
      itn2 => dg%loading%itn2
      jts1 => dg%loading%jts1
      jte1 => dg%loading%jte1
      jtn1 => dg%loading%jtn1
! ==============================================================================
      green_out_Z => dg%loading%green_out_Z
      xfftbuf => dg%loading%xfftbuf
      yfftbuf => dg%loading%yfftbuf
#ifdef REAL_FFT
      realbuf => dg%loading%realbuf
#endif
      defZmap => dg%loading%defZmap
      defZmap1 => dg%loading%defZmap1

#ifdef _OPENMP
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif
!$omp parallel
#ifndef MPI
!$omp do private(j, i)
      do it = 0, nthreads-1
         do j = jts1(it), jte1(it)
            do i = 1, N_X
#ifndef REAL_FFT
               xfftbuf(i,j) = dcmplx(0.0d0,0.0d0)
#else
               realbuf(i,j) = 0.0d0
#endif
            end do
         end do
      end do

!$omp do private(i)
      do j = 1, nlat
         do i = 1, nlon
            if(ifz(i,j) == 1) then
#ifndef REAL_FFT
               xfftbuf(i,j) = dcmplx(hz(i,j),0.0d0)
#else
               realbuf(i,j) = hz(i,j)
#endif
            end if
         end do
      end do

      ! Forward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
#ifndef REAL_FFT
         call dfftw_execute(xplan_forward(it))
#else
         call dfftw_execute_dft_r2c(xplan_forward(it),realbuf(1,jts1(it)),xfftbuf(1,jts1(it)))
#endif
      end do
#else
#ifndef REAL_FFT
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_X,jtn1(it),xfftbuf(1,jts1(it)),1,N_X+1,1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+2,1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif
#endif

      ! Transposiion: (N_X,N_Y) -> (N_Y,N_X)
!$omp do private(i, j)
      do it = 0, nthreads-1
         do i = its2(it), ite2(it)
            do j = 1, N_Y
#ifndef __SX__
               yfftbuf(j,i) = xfftbuf(i,j)
#else
#ifndef REAL_FFT
               yfftbuf(j,i) = xfftbuf(i,j)
#else
               yfftbuf(j,i) = dcmplx(realbuf(2*i-1,j),realbuf(2*i,j))
#endif
#endif
            end do
         end do
      end do

      ! Forward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_forward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+1,1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

!$omp do private(i, j)
      do it = 0, nthreads-1
         do i = its2(it), ite2(it)
            do j = 1, N_Y
               yfftbuf(j,i) = yfftbuf(j,i)*green_out_Z(j,i)*coef_norm
            end do
         end do
      end do

      ! Backward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_backward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+1,-1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

      ! Transposiion: (N_Y,N_X) -> (N_X,N_Y)
!$omp do private(j, i)
      do it = 0, nthreads-1
         do j = jts1(it), jte1(it)
#ifndef REAL_FFT
            do i = 1, N_X
#else
            do i = 1, N_X/2+1
#endif
#ifndef __SX__
               xfftbuf(i,j) = yfftbuf(j,i)
#else
#ifndef REAL_FFT
               xfftbuf(i,j) = yfftbuf(j,i)
#else
               realbuf(2*i-1,j) = dble (yfftbuf(j,i))
               realbuf(2*i,  j) = dimag(yfftbuf(j,i))
#endif
#endif
            end do
         end do
      end do

      ! Backward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
#ifndef REAL_FFT
         call dfftw_execute(xplan_backward(it))
#else
         call dfftw_execute_dft_c2r(xplan_backward(it),xfftbuf(1,jts1(it)),realbuf(1,jts1(it)))
#endif
      end do
#else
#ifndef REAL_FFT
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_X,jtn1(it),xfftbuf(1,jts1(it)),1,N_X+1,-1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+2,-1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif
#endif

!$omp do private(i)
      do j = 1, nlat
         do i = 1, nlon
#ifndef REAL_FFT
            defZmap(i,j) = dble(xfftbuf(i,j))
#else
            defZmap(i,j) = realbuf(i,j)
#endif
         end do
      end do

!$omp do private(i)
      do j = 1, nlat
         do i = 1, nlon
            if(ifz(i,j) == 1) then
               hz(i,j) = hz(i,j) - (defZmap(i,j)-defZmap1(i,j))
               defZmap1(i,j) = defZmap(i,j)
            end if
         end do
      end do
#else
      ! Transposiion: (nx0,ny0) -> (N_X,ny1)
!$omp do private(i, ind)
      do j = 1, ny0
         do i = 1, nx0
            ind = (j-1)*nx0 + i
            sendbuf1(ind) = 0.0d0
         end do
      end do

!$omp do private(j_, i, i_, ind)
      do j = 1, ny0
         j_ = jbias + j
         do i = 1, nx0
            i_ = ibias + i
            ind = (j-1)*nx0 + i
            if(ifz(i_,j_) == 1) then
               sendbuf1(ind) = hz(i_,j_)
            end if
         end do
      end do

!$omp single
      call MPI_Alltoallv(sendbuf1, sendcounts1, sdispls1, MPI_DOUBLE_PRECISION, &
                         recvbuf1, recvcounts1, rdispls1, MPI_DOUBLE_PRECISION, MPI_X_WORLD, ierr)
!$omp end single

!$omp do private(i)
      do j = 1, ny1
#ifndef REAL_FFT
#ifndef __SX__
         do i = 1, N_X
#else
         do i = 1, N_X+1
#endif
            xfftbuf(i,j) = dcmplx(0.0d0,0.0d0)
#else
#ifndef __SX__
         do i = 1, N_X
#else
         do i = 1, N_X+2
#endif
            realbuf(i,j) = 0.0d0
#endif
         end do
      end do

!$omp do private(xst, xlen, j, i, ind, i_)
      do p = 0, npx-1
         xst = rdispls1(p)/ny1
         xlen = recvcounts1(p)/ny1
         do j = 1, ny1
            do i = 1, xlen
               ind = rdispls1(p) + (j-1)*xlen + i
               i_ =  xst + i
#ifndef REAL_FFT
               xfftbuf(i_,j) = dcmplx(recvbuf1(ind),0.0d0)
#else
               realbuf(i_,j) = recvbuf1(ind)
#endif
            end do
         end do
      end do

      ! Forward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
#ifndef REAL_FFT
         call dfftw_execute(xplan_forward(it))
#else
         call dfftw_execute_dft_r2c(xplan_forward(it),realbuf(1,jts1(it)),xfftbuf(1,jts1(it)))
#endif
      end do
#else
#ifndef REAL_FFT
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_X,jtn1(it),xfftbuf(1,jts1(it)),1,N_X+1,1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+2,1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif
#endif

      ! Transposiion: (N_X,ny1) -> (N_Y,nx2)
#ifndef __SX__
#ifndef REAL_FFT
!$omp do private(j, ind)
      do i = 1, N_X
#else
!$omp do private(j, ind)
      do i = 1, N_X/2+1
#endif
         do j = 1, ny1
            ind = ny1*(i-1) + j
            sendbuf2(ind) = xfftbuf(i,j)
         end do
      end do
#else
!$omp do private(i, ind)
      do j = 1, ny1
#ifndef REAL_FFT
         do i = 1, N_X
#else
         do i = 1, N_X/2+1
#endif
            ind = ny1*(i-1) + j
#ifndef REAL_FFT
            sendbuf2(ind) = xfftbuf(i,j)
#else
            sendbuf2(ind) = dcmplx(realbuf(2*i-1,j),realbuf(2*i,j))
#endif
         end do
      end do
#endif

!$omp single
#ifndef MULTI
      call MPI_Alltoallv(sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, &
                         recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoallv(sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, &
                         recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, MPI_MEMBER_WORLD, ierr)
#endif
!$omp end single

!$omp do private(j)
      do i = 1, nx2
#ifndef __SX__
         do j = 1, N_Y
#else
         do j = 1, N_Y+1
#endif
            yfftbuf(j,i) = dcmplx(0.0d0,0.0d0)
         end do
      end do

!$omp do private(yst, ylen, i, j, ind, j_)
      do p = 0, nprocs-1
         yst = rdispls2(p)/nx2
         ylen = recvcounts2(p)/nx2
         do i = 1, nx2
            do j = 1, ylen
               ind = rdispls2(p) + (i-1)*ylen + j
               j_ =  yst + j
               yfftbuf(j_,i) = recvbuf2(ind)
            end do
         end do
      end do

      ! Forward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_forward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+1,1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

      ! Calc. F(hz)*F(g)
!$omp do private(i, j)
      do it = 0, nthreads-1
         do i = its2(it), ite2(it)
            do j = 1, N_Y
               yfftbuf(j,i) = yfftbuf(j,i)*green_out_Z(j,i)*coef_norm
            end do
         end do
      end do

      ! Backward FFT on Y-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
         call dfftw_execute(yplan_backward(it))
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_Y,itn2(it),yfftbuf(1,its2(it)),1,N_Y+1,-1,ifax_y(1,it),trigs_y(1,it),work_y(1,its2(it)),ierr)
      end do
#endif

      ! Transposiion: (N_Y,nx2) -> (N_X,ny1)
!$omp do private(yst, ylen, i, j, ind, j_)
      do p = 0, nprocs-1
         yst = rdispls2(p)/nx2
         ylen = recvcounts2(p)/nx2
         do i = 1, nx2
            do j = 1, ylen
               ind = rdispls2(p) + (i-1)*ylen + j
               j_ =  yst + j
               recvbuf2(ind) = yfftbuf(j_,i)
            end do
         end do
      end do

!$omp single
#ifndef MULTI
      call MPI_Alltoallv(recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, &
                         sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
#else
      call MPI_Alltoallv(recvbuf2, recvcounts2, rdispls2, MPI_DOUBLE_COMPLEX, &
                         sendbuf2, sendcounts2, sdispls2, MPI_DOUBLE_COMPLEX, MPI_MEMBER_WORLD, ierr)
#endif
!$omp end single

#ifndef __SX__
!$omp do private(i, ind)
      do j = 1, ny1
#ifndef REAL_FFT
         do i = 1, N_X
#else
         do i = 1, N_X/2+1
#endif
            ind = ny1*(i-1) + j
            xfftbuf(i,j) = sendbuf2(ind)
         end do
      end do
#else
!$omp do private(i, ind)
      do j = 1, ny1
#ifndef REAL_FFT
         do i = 1, N_X
#else
         do i = 1, N_X/2+1
#endif
            ind = ny1*(i-1) + j
#ifndef REAL_FFT
            xfftbuf(i,j) = sendbuf2(ind)
#else
            realbuf(2*i-1,j) = dble (sendbuf2(ind))
            realbuf(2*i,  j) = dimag(sendbuf2(ind))
#endif
         end do
      end do
#endif

      ! Backward FFT on X-direction
#ifndef __SX__
!$omp do
      do it = 0, nthreads-1
#ifndef REAL_FFT
         call dfftw_execute(xplan_backward(it))
#else
         call dfftw_execute_dft_c2r(xplan_backward(it),xfftbuf(1,jts1(it)),realbuf(1,jts1(it)))
#endif
      end do
#else
#ifndef REAL_FFT
!$omp do private(ierr)
      do it = 0, nthreads-1
         call ZFCMBF(N_X,jtn1(it),xfftbuf(1,jts1(it)),1,N_X+1,-1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#else
!$omp do private(ierr)
      do it = 0, nthreads-1
         call DFRMBF(N_X,jtn1(it),realbuf(1,jts1(it)),1,N_X+2,-1,ifax_x(1,it),trigs_x(1,it),work_x(1,jts1(it)),ierr)
      end do
#endif
#endif

      ! Transposiion: (N_X,ny1) -> (nx0,ny0)
!$omp do private(xst, xlen, j, i, ind, i_)
      do p = 0, npx-1
         xst = rdispls1(p)/ny1
         xlen = recvcounts1(p)/ny1
         do j = 1, ny1
            do i = 1, xlen
               ind = rdispls1(p) + (j-1)*xlen + i
               i_ =  xst + i
#ifndef REAL_FFT
               recvbuf1(ind) = dble(xfftbuf(i_,j))
#else
               recvbuf1(ind) = realbuf(i_,j)
#endif
            end do
         end do
      end do

!$omp single
      call MPI_Alltoallv(recvbuf1, recvcounts1, rdispls1, MPI_DOUBLE_PRECISION, &
                         sendbuf1, sendcounts1, sdispls1, MPI_DOUBLE_PRECISION, MPI_X_WORLD, ierr)
!$omp end single

!$omp do private(i, ind)
      do j = 1, ny0
         do i = 1, nx0
            ind = (j-1)*nx0 + i
            defZmap(i,j) = sendbuf1(ind)
         end do
      end do

      ! Update hz and defZmap1
!$omp do private(j_, i, i_)
      do j = 1, ny0
         j_ = jbias + j
         do i = 1, nx0
            i_ = ibias + i
            if(ifz(i_,j_) == 1) then
               hz(i_,j_) = hz(i_,j_) - (defZmap(i,j)-defZmap1(i,j))
               defZmap1(i,j) = defZmap(i,j)
            end if
         end do
      end do
#endif
!$omp end parallel

      return
   end subroutine loading_run

end module mod_loading
