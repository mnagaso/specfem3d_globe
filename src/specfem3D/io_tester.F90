module io_bandwidth
    use mpi, only: MPI_Wtime

    implicit none
    double precision :: start_time, end_time, time_delta
    integer(kind=8) :: bytes_written  ! Use 64-bit integer to avoid overflow for large I/O

  contains

    subroutine initialize_bytes_written()
      bytes_written = 0
      time_delta = 0.0
    end subroutine initialize_bytes_written

    subroutine start_timer()
        start_time = MPI_Wtime()
    end subroutine start_timer

    subroutine stop_timer()
        end_time = MPI_Wtime()
        time_delta = time_delta + (end_time - start_time)
    end subroutine stop_timer

    subroutine set_bytes_written(bytes)
      implicit none
      integer(kind=8), intent(in) :: bytes
      bytes_written = bytes
    end subroutine set_bytes_written

    subroutine set_bytes_written_from_array(element_size, num_elements)
        implicit none
        integer, intent(in):: element_size, num_elements

        ! element_size is the size of each element in bits (as storage_size returns in bits)

        ! Calculate the total bytes written and add to the existing value
        ! Use 64-bit arithmetic to avoid overflow
        bytes_written = bytes_written + int(element_size, kind=8) * int(num_elements, kind=8) / 8_8 ! Convert bits to bytes

      end subroutine set_bytes_written_from_array


    subroutine calculate_bandwidth_this_proc()
      implicit none
      double precision :: elapsed_time, bandwidth

      elapsed_time = time_delta
      if (elapsed_time > 0.0) then
        bandwidth = bytes_written / (elapsed_time * 1.0e6)  ! Bandwidth in MB/s
        print *, 'I/O Bandwidth: ', bandwidth, ' MB/s'
      else
        print *, 'Elapsed time is zero, cannot calculate bandwidth'
      endif
    end subroutine calculate_bandwidth_this_proc

    ! calculate total bandwidth across all processes
    subroutine calculate_bandwidth_all_procs()
      use specfem_par
      use shared_parameters
      implicit none
      integer :: i, unit_number, ierr
      double precision :: elapsed_time, max_elapsed_time, &
    	      bandwidth, total_bandwidth, current_time
      double precision :: bytes_written_dp, total_bytes_dp  ! Use double precision for large byte counts
      character(len=20) :: filename
      character(len=10) :: mygroup_str
      character(len=8) :: date_str
      character(len=10) :: time_str
      integer :: datetime_values(8)
      logical :: file_exists

      ! Convert mygroup to a character string
      write(mygroup_str, '(I0)') mygroup

      ! Create the filename
      filename = 'io_band_' // trim(adjustl(mygroup_str)) // '.txt'
      unit_number = 1000010010 + mygroup

      ! Check if the file exists
      inquire(file=filename, exist=file_exists)

      ! Open the file for writing
      if (myrank == 0) then
        ! create the file if it does not exist
        if (file_exists) then
          open(unit=unit_number, file=filename, status='old', action='write', position='append', iostat=ierr)
        else
          open(unit=unit_number, file=filename, status='new', action='write', iostat=ierr)
        end if
        if (ierr /= 0) then
            print*, 'Error opening file: ', filename
            stop
        end if
        write(unit_number, *) '-------------------------------------------------------------------------------'
        close(unit_number)
      end if

      call synchronize_all()

      elapsed_time = time_delta
      if (elapsed_time > 0.0) then
          bandwidth = dble(bytes_written) / (elapsed_time * 1.0d6)  ! Bandwidth in MB/s

          ! calculate total bytes written across all processes using double precision
          bytes_written_dp = dble(bytes_written)
          call sum_all_dp(bytes_written_dp, total_bytes_dp)
          ! get the maximum elapsed time across all processes
          call max_all_dp(elapsed_time, max_elapsed_time)

            ! calculate total bandwidth
            total_bandwidth = total_bytes_dp / (max_elapsed_time * 1.0d6) ! Bandwidth in MB/s

            ! Each process writes to the file sequentially
            do i = 0, NPROCTOT_VAL-1
              if (myrank == i) then
                open(unit=unit_number, file=filename, status='old', action='write', position='append', iostat=ierr)
                if (ierr /= 0) then
                  print*, 'Error opening file: ', filename
                  stop
                end if

                ! record MPI wall-clock time and system date/time for this rank's bandwidth entry
                current_time = MPI_Wtime()
                call date_and_time(date=date_str, time=time_str, values=datetime_values)

                write(unit_number, '(A, I0, A, I0, A, I0, A, F12.6, A, F12.6, A, F24.12, A, A, A, A)') &
                  'mygroup: ', mygroup, ', myrank: ', myrank, ', bytes_written: ', bytes_written, &
                  ', elapsed_time (s): ', elapsed_time, ', bandwidth: ', bandwidth, ' MB/s, mpi_wtime (s): ', current_time, &
                  ', date: ', trim(date_str), ', time: ', trim(time_str)
                close(unit_number)
              end if
              call synchronize_all()
            end do

            ! Only the root process writes the total bandwidth
            if (myrank == 0) then
              open(unit=unit_number, file=filename, status='old', action='write', position='append', iostat=ierr)
              if (ierr /= 0) then
                print*, 'Error opening file: ', filename
                stop
              end if

              ! record MPI wall-clock time and system date/time for this I/O bandwidth measurement
              current_time = MPI_Wtime()
              call date_and_time(date=date_str, time=time_str, values=datetime_values)

              write(unit_number, '(A, I0, A, F20.0, A, F12.6, A, F12.6, A, F24.12, A, A, A, A, A, I0)') &
	          'mygroup: ', mygroup, ', total_bytes_written: ', total_bytes_dp, ', max_elapsed_time (s): ', &
	          max_elapsed_time, ', total_bandwidth: ', total_bandwidth, ' MB/s, mpi_wtime (s): ', current_time, &
            ', date: ', trim(date_str), ', time: ', trim(time_str), ', ms: ', datetime_values(8)
              close(unit_number)
            end if
      end if
    end subroutine calculate_bandwidth_all_procs

  end module io_bandwidth