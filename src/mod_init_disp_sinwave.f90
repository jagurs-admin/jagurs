#include "real.h"
module mod_init_disp_sinwave
use mod_grid
use mod_rwg
implicit none

real(kind=8), private :: T = 20.0d0 ! [s]
real(kind=8), private :: a = 0.03d0 ! [m]
integer(kind=4), private :: x_index = 0
integer(kind=4), private :: y_index = 0
! === FILE INPUT ===============================================================
character(len=64), private :: input_file = ''
real(kind=REAL_BYTE), private, allocatable, dimension(:) :: tlist, hlist
! === DEBUG for restart with sinwave/file-input 2016/02/18 =====================
!integer(kind=4), private :: lcount = 0, ccount = 0
integer(kind=4), private :: lcount = 0
integer(kind=4) :: ccount = 0
! ==============================================================================
! ==============================================================================

real(kind=8), private :: pi

contains

#ifndef MULTI
   subroutine specify_sinwave_params(file)
#else
   subroutine specify_sinwave_params(input_dirname, file)
#endif
      character(len=256), intent(in) :: file
#ifdef MULTI
      character(len=128), intent(in) :: input_dirname
      character(len=256) :: file_multi
#endif
! === FILE INPUT ===============================================================
!     namelist /sinwave/ T, a, x_index, y_index
      real(kind=REAL_BYTE) :: t_in, h_in
      namelist /sinwave/ T, a, x_index, y_index, input_file
! ==============================================================================

#ifndef MULTI
      open(1,file=trim(file),action='read',status='old',form='formatted')
#else
      file_multi = trim(input_dirname) // trim(file)
      open(1,file=trim(file_multi),action='read',status='old',form='formatted')
#endif
      read(1,sinwave)
      close(1)

! === FILE INPUT ===============================================================
      if(input_file == '') then
! ==============================================================================
         write(6,'(a)') '============================================================'
         write(6,'(a)') '=== Check Sin Wave initialization! [Begin] ================='
         write(6,'(a)') '============================================================'
         write(6,'(a,a)') '- Filename: ', trim(file)
         write(6,'(a,f15.6)') '- Cycle[s] (T):       ', T
         write(6,'(a,f15.6)') '- Wave height[m] (a): ', a
         write(6,'(a,i8)')    '- X-index (x_index):  ', x_index
         write(6,'(a,i8)')    '- Y-index (y_index):  ', y_index
         write(6,'(a)') '============================================================'
         write(6,'(a)') '=== Check Sin Wave initialization! [End] ==================='
         write(6,'(a)') '============================================================'

         pi = 2.0d0*dasin(1.0d0)
! === FILE INPUT ===============================================================
      else
         write(6,'(a)') '============================================================'
         write(6,'(a)') '=== Check Sin Wave initialization! [Begin] ================='
         write(6,'(a)') '============================================================'
         write(6,'(a,a)') '- Filename: ', trim(file)
         write(6,'(a,a)')     '- Input File (input_file): ', trim(input_file)
         write(6,'(a,i8)')    '- X-index (x_index):       ', x_index
         write(6,'(a,i8)')    '- Y-index (y_index):       ', y_index
         write(6,'(a)') '============================================================'
         write(6,'(a)') '=== Check Sin Wave initialization! [End] ==================='
         write(6,'(a)') '============================================================'

#ifndef MULTI
         open(1,file=trim(input_file),action='read',status='old',form='formatted')
#else
         file_multi = trim(input_dirname) // trim(input_file)
         open(1,file=trim(file_multi),action='read',status='old',form='formatted')
#endif

         lcount=0
         do
            read(1, *, end=100, err=200) t_in, h_in
            lcount = lcount+1
200         continue
         end do
100      continue

         allocate(hlist(lcount)); hlist = 0.0d0
         allocate(tlist(lcount)); tlist = 0.0d0

         rewind(1)
         lcount=0
         do
            read(1, *, end=101, err=201) t_in, h_in
            lcount = lcount+1
            hlist(lcount) = h_in
            tlist(lcount) = t_in
201         continue
         end do
101      continue
      end if
! ==============================================================================

      return
   end subroutine specify_sinwave_params

#ifndef MPI
   subroutine hrise_rwg_sin(wfld,time,nlon,nlat)
#else
   subroutine hrise_rwg_sin(wfld,time,nlon,nlat,kx,ky)
#endif
      type(wave_arrays), target, intent(inout) :: wfld
      real(kind=REAL_BYTE), intent(in) :: time
      integer(kind=4), intent(in) :: nlon, nlat
#ifdef MPI
      integer(kind=4), intent(in) :: kx, ky
#endif

      real(kind=REAL_BYTE), pointer, dimension(:,:) :: hz
      integer(kind=4) :: i, j
#ifdef MPI
      integer(kind=4) :: ig, jg
#endif

      hz => wfld%hz

! === FILE INPUT ===============================================================
      if(lcount == 0) then
! ==============================================================================
#ifndef MPI
         if(x_index /= 0) then
            i = x_index
!$omp parallel do
            do j = 1, nlat
               hz(i,j) = a*dsin(2.0d0*pi/T*time)
            end do
         end if

         if(y_index /= 0) then
            j = y_index
!$omp parallel do
            do i = 1, nlon
               hz(i,j) = a*dsin(2.0d0*pi/T*time)
            end do
         end if
#else
         if(x_index /= 0) then
            do i = 1, nlon
               ig = kx+i-1
               if(ig == x_index) then
!$omp parallel do
                  do j = 1, nlat
                     hz(i,j) = a*dsin(2.0d0*pi/T*time)
                  end do
               end if
            end do
         end if

         if(y_index /= 0) then
            do j = 1, nlat
               jg = ky+j-1
               if(jg == y_index) then
!$omp parallel do
                  do i = 1, nlon
                     hz(i,j) = a*dsin(2.0d0*pi/T*time)
                  end do
               end if
            end do
         end if
#endif
! === FILE INPUT ===============================================================
      else
         ccount = ccount + 1
         if(ccount > lcount) return
!        write(6,'(i8,3e15.6)') ccount, time, tlist(ccount), hlist(ccount)
#ifndef MPI
         if(x_index /= 0) then
            i = x_index
!$omp parallel do
            do j = 1, nlat
               hz(i,j) = hlist(ccount)
            end do
         end if

         if(y_index /= 0) then
            j = y_index
!$omp parallel do
            do i = 1, nlon
               hz(i,j) = hlist(ccount)
            end do
         end if
#else
         if(x_index /= 0) then
            do i = 1, nlon
               ig = kx+i-1
               if(ig == x_index) then
!$omp parallel do
                  do j = 1, nlat
                    hz(i,j) = hlist(ccount)
                  end do
               end if
            end do
         end if

         if(y_index /= 0) then
            do j = 1, nlat
               jg = ky+j-1
               if(jg == y_index) then
!$omp parallel do
                  do i = 1, nlon
                     hz(i,j) = hlist(ccount)
                  end do
               end if
            end do
         end if
#endif
      end if
! ==============================================================================

      return
   end subroutine hrise_rwg_sin

end module mod_init_disp_sinwave
