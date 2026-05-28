!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
! Timing accounting module for separating compute, IO, wait, and idle times.
!
! Each process measures its own local wall-clock time using MPI_Wtime()
! without any synchronization barriers. At the end of the time loop,
! each process writes its timing breakdown to a shared log file.
!
! Categories:
!   - compute : time in force computation, time stepping, kernels
!   - io      : time in write/send (compute nodes) or recv/write (IO nodes)
!   - wait_io : time compute nodes block in wait_all_send (HDF5_IO_NODES > 0)
!   - idle_io : time IO nodes block in MPI_Probe waiting for messages
!
!=====================================================================

module timing_accounting

  use mpi, only: MPI_Wtime

  implicit none
  private :: MPI_Wtime

  ! accumulators (seconds)
  double precision :: time_compute = 0.0d0
  double precision :: time_io      = 0.0d0
  double precision :: time_wait_io = 0.0d0
  double precision :: time_idle_io = 0.0d0

  ! scratch start-time variables
  double precision :: t_compute_start = 0.0d0
  double precision :: t_io_start      = 0.0d0
  double precision :: t_wait_start    = 0.0d0
  double precision :: t_idle_start    = 0.0d0

contains

  subroutine timing_reset()
    time_compute = 0.0d0
    time_io      = 0.0d0
    time_wait_io = 0.0d0
    time_idle_io = 0.0d0
  end subroutine timing_reset

  ! ---- compute timing ----
  subroutine timing_compute_start()
    t_compute_start = MPI_Wtime()
  end subroutine timing_compute_start

  subroutine timing_compute_stop()
    time_compute = time_compute + (MPI_Wtime() - t_compute_start)
  end subroutine timing_compute_stop

  ! ---- IO timing ----
  subroutine timing_io_start()
    t_io_start = MPI_Wtime()
  end subroutine timing_io_start

  subroutine timing_io_stop()
    time_io = time_io + (MPI_Wtime() - t_io_start)
  end subroutine timing_io_stop

  ! ---- wait-for-IO timing (compute nodes, HDF5_IO_NODES > 0) ----
  subroutine timing_wait_start()
    t_wait_start = MPI_Wtime()
  end subroutine timing_wait_start

  subroutine timing_wait_stop()
    time_wait_io = time_wait_io + (MPI_Wtime() - t_wait_start)
  end subroutine timing_wait_stop

  ! ---- idle timing (IO nodes) ----
  subroutine timing_idle_start()
    t_idle_start = MPI_Wtime()
  end subroutine timing_idle_start

  subroutine timing_idle_stop()
    time_idle_io = time_idle_io + (MPI_Wtime() - t_idle_start)
  end subroutine timing_idle_stop

  ! ---- reporting (no MPI barriers) ----
  subroutine timing_report(rank, group, is_io_node, total_elapsed)
    implicit none
    integer, intent(in) :: rank, group
    logical, intent(in) :: is_io_node
    double precision, intent(in) :: total_elapsed

    integer :: unit_number, ierr
    character(len=80) :: filename
    character(len=12) :: group_str
    character(len=12) :: rank_str
    character(len=8)  :: date_str
    character(len=10) :: time_str
    character(len=16) :: role_str
    double precision :: time_other, current_time

    ! compute unaccounted time
    time_other = total_elapsed - time_compute - time_io - time_wait_io - time_idle_io

    ! role label
    if (is_io_node) then
      role_str = 'io_server'
    else
      role_str = 'compute'
    endif

    ! one file per rank+role to avoid filename collision between io_server and
    ! compute sub-communicator when both have myrank=0 (no barriers needed)
    write(group_str, '(I0)') group
    write(rank_str,  '(I0)') rank
    filename = 'timing_acct_group' // trim(adjustl(group_str)) // &
               '_' // trim(adjustl(role_str)) // &
               '_rank' // trim(adjustl(rank_str)) // '.txt'

    ! get current wall-clock and date for the log entry
    current_time = MPI_Wtime()
    call date_and_time(date=date_str, time=time_str)

    ! each rank writes exclusively to its own file (no contention, no barriers)
    open(newunit=unit_number, file=filename, status='unknown', action='write', &
         position='append', iostat=ierr)
    if (ierr /= 0) then
      print *, 'timing_accounting: Error opening ', trim(filename), ' rank=', rank, ' ierr=', ierr
      return
    endif

    write(unit_number, '(A,I0,A,I0,A,A,A,F14.6,A,F14.6,A,F14.6,A,F14.6,A,F14.6,A,F14.6,A,F24.12,A,A,A,A)') &
      'mygroup: ', group, &
      ', myrank: ', rank, &
      ', role: ', trim(role_str), &
      ', compute (s): ', time_compute, &
      ', io (s): ', time_io, &
      ', wait_io (s): ', time_wait_io, &
      ', idle_io (s): ', time_idle_io, &
      ', other (s): ', time_other, &
      ', total (s): ', total_elapsed, &
      ', mpi_wtime (s): ', current_time, &
      ', date: ', trim(date_str), &
      ', time: ', trim(time_str)

    close(unit_number)

  end subroutine timing_report

end module timing_accounting
