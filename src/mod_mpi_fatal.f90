module mod_mpi_fatal
use mpi
implicit none
contains

   subroutine fatal_error(err)
      integer(kind=4), intent(in) :: err
      integer(kind=4) :: ierr
      write(0,'(a,i0)') 'error code = ', err
      call MPI_Abort(MPI_COMM_WORLD, err, ierr)
      call MPI_Finalize(ierr)
      stop
   end subroutine fatal_error

   subroutine stop_process()
      integer(kind=4) :: ierr
      call MPI_Finalize(ierr)
      stop
   end subroutine stop_process

end module mod_mpi_fatal
