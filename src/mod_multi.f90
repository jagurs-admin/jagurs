module mod_multi
#ifdef MULTI
use mpi
implicit none

! === Split Dir ================================================================
!integer(kind=4), parameter :: max_paramset = 10
integer(kind=4), parameter :: max_paramset = 1
! ==============================================================================

integer(kind=4) :: num_members, g_nprocs, g_myrank, member_id
character(len=128) :: members_dir = '', command = ''
#ifdef NCDIO
character(len=256) :: temp_filename = ''
#endif
#ifndef MPI
character(len=128) :: suffix = '', stdout = 'stdout'
integer(kind=4) :: ierr
#else
integer(kind=4) :: MPI_MEMBER_WORLD
#endif

character(len=128), allocatable, dimension(:) :: input_files
integer(kind=4), allocatable, dimension(:)  :: member_ids
! === Split Dir ================================================================
character(len=128) :: input_dirname = 'input.'
! ==============================================================================

private :: usage

contains

   subroutine get_arguments()
      integer(kind=4) :: num_arg, len_arg, num_paramset
      character(len=128) :: program_name, arg
      logical :: invalid

      integer(kind=4), dimension(max_paramset) :: istart, iend
      character(len=128), dimension(max_paramset) :: filenames

      integer(kind=4) :: i, j, k, ist, ien, ierr

      num_arg = COMMAND_ARGUMENT_COUNT()
      call  GET_COMMAND_ARGUMENT(0, program_name)

      if((num_arg < 2) .or. (num_arg > 2*max_paramset+1) .or. (mod(num_arg, 2) /= 0)) then
         goto 200
      end if

      i = 0
      num_paramset = 0
      do while(i < num_arg)
         i = i + 1
         num_paramset = num_paramset + 1
         call  GET_COMMAND_ARGUMENT(i, arg)
         arg = ADJUSTL(arg)
         len_arg = LEN_TRIM(arg)
         invalid = .true.
         do j = 1, len_arg
            if(arg(j:j) .eq. '-') then
               invalid = .false.
               exit
             end if
         end do
         if(invalid) goto 200
#ifndef __GFORTRAN__
         read(arg(1:j-1), '(i)', err=200) istart(num_paramset)
         read(arg(j+1:len_arg), '(i)', err=200) iend(num_paramset)
#else
         read(arg(1:j-1), *, err=200) istart(num_paramset)
         read(arg(j+1:len_arg), *, err=200) iend(num_paramset)
#endif

         i = i + 1
         call  GET_COMMAND_ARGUMENT(i, filenames(num_paramset))
      end do

      if(g_myrank == 0) then
         write(0,'(a,i3)') 'Number of paramset: ', num_paramset
         write(0,'(/,a)') 'No StartID EndID    InputFileName'
         do i = 1, num_paramset
            write(0,'(i3,2i8,a,a)') i, istart(i), iend(i), ' ', trim(filenames(i))
         end do
      end if

      do i = 1, num_paramset
         if(istart(i) > iend(i)) then
            if(g_myrank == 0) then
               write(0,'(a,i3)')'Error! StartID is grather than EndID in No. ', i
            end if
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Abort(MPI_COMM_WORLD, 201, ierr)
         end if
      end do

      ist = minval(istart(1:num_paramset))
      ien = maxval(iend(1:num_paramset))
      num_members = ien - ist + 1
      if(g_myrank == 0) then
         write(0,'(/,a,i8)') 'Number of members: ',  num_members
         write(0,'(a,i8)') 'From: ', ist
         write(0,'(a,i8)') 'To:   ', ien
      end if
   
      allocate(input_files(num_members))
      allocate(member_ids(num_members))
      input_files = ''

      do i = 1, num_paramset
         do j = istart(i), iend(i)
            k = j - ist + 1
            input_files(k) = trim(filenames(i))
         end do
      end do

      if(g_myrank == 0) then
         write(0,'(/,a)') 'MemberID InputFileName'
      end if
      do i = 1, num_members
         member_ids(i) = i + ist - 1
         if(g_myrank == 0) then
            write(0,'(i8,a,a)') member_ids(i), ' ', trim(input_files(i))
         end if
      end do
   
      do i = 1, num_members
         if(input_files(i) == '') then
            if(g_myrank == 0) then
               write(0,*) 'Input file for member ID ',  i+ist-1, ' is NOT specified!!!'
            end if
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Abort(MPI_COMM_WORLD, 202, ierr)
         end if
      end do

      return

200   if(g_myrank == 0) call usage(program_name)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Abort(MPI_COMM_WORLD, 200, ierr)
 
   end subroutine get_arguments

   subroutine usage(program_name)
      character(len=128), intent(in) :: program_name

! === Split Dir ================================================================
!     write(0,'(a,a,a,/)') 'usage: ', trim(program_name), ' paramset1 [paramset2 ...]'
!     write(0,'(a)')       '   where ''paramsetN'' is the set of arguments ''sid-eid'' and ''filename''.'
!     write(0,'(a,i3,/)')  '   N must be less or equal ', max_paramset
      write(0,'(a,a,a,/)') 'usage: ', trim(program_name), ' paramset'
      write(0,'(a)')       '   where ''paramset'' is the set of arguments ''sid-eid'' and ''filename''.'
! ==============================================================================
      write(0,'(a)')       '   sid: Start ID of target membars. (integer, >= 0)'
      write(0,'(a)')       '   eid: End ID of target membars. (integer, >= sid)'
      write(0,'(a,/)')     '   filename: Input filename for members whose ID is ''sid''-''eid''. (strings)'
      write(0,'(a,/)')     '   For example,'
! === Split Dir ================================================================
!     write(0,'(a,a,a)')   '      ', trim(program_name), ' 1-100 tsun0.par 101-200 tsun1.par'
      write(0,'(a,a,a)')   '      ', trim(program_name), ' 1-100 tsun.par'
! ==============================================================================

      return
   end subroutine usage
#endif

end module mod_multi

