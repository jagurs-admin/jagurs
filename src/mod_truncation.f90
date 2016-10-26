module mod_truncation
implicit none

real(kind=8), private :: start_time, current_time
private :: wallclock

contains

   subroutine start_trunc()
      call wallclock(start_time)
      return
   end subroutine start_trunc

   logical function check_trunc(max_time)
      integer(kind=4), intent(in) :: max_time
      check_trunc = .false.
      call wallclock(current_time)
      if(current_time - start_time > max_time) check_trunc = .true.
      return
   end function check_trunc

   integer(kind=4) function truncation_time(max_time)
      character(len=64), intent(in) :: max_time
      integer(kind=4) :: hours, minutes, seconds
      integer(kind=4) :: l, c1, c2

      hours   = 0
      minutes = 0
      seconds = 0

      l = len(trim(max_time))

      if(l /= 0) then
         c1 = index(max_time(1:l),':')
         if(c1 == 0) then
            read(max_time(1:l),'(i20)') seconds
         else
            c2 = index(max_time(c1+1:l),':')
            if(c2 == 0) then
               read(max_time(1:c1-1),'(i20)') minutes
               read(max_time(c1+1:l),'(i20)') seconds
            else
               c2 = c1 + c2
               read(max_time(1:c1-1),'(i20)') hours
               read(max_time(c1+1:c2-1),'(i20)') minutes
               read(max_time(c2+1:l),'(i20)') seconds
            end if
         end if
      end if

      truncation_time = 3600*hours + 60*minutes + seconds

      return
   end function truncation_time

   subroutine wallclock(t)
      real(kind=8), intent(out) :: t
#ifndef __FUJITSU
#ifndef NO_SYSTEM_CLOCK
      integer(kind=8) :: c, c_rate
#else
! === This should be wall-clock time. Not user time. ===========================
!     real :: etime, tarray(2)
      real(kind=8) :: count, count_rate, count_max
! ==============================================================================
#endif

#ifndef NO_SYSTEM_CLOCK
      call system_clock(c, c_rate)

      t = dble(c)/dble(c_rate)
#else
! === This should be wall-clock time. Not user time. ===========================
!     t = etime(tarray)
      call system_clockx(count, count_rate, count_max)

      t = count/count_rate
! ==============================================================================
#endif
#else
      call gettod(t)
      t = t*1.0d-6
#endif

      return
   end subroutine wallclock

end module mod_truncation
