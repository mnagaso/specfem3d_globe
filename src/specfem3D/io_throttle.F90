!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!> @file io_throttle.F90
!! @brief I/O throttling module for controlling checkpoint write rates
!!
!! This module provides:
!! - Pre-checkpoint delay (sleep before writing)
!! - Bandwidth-limited writes
!! - Adaptive throttling based on external control file
!! - Runtime adjustable throttle parameters via control file
!!
!! The throttle parameters can be set via Par_file or dynamically adjusted
!! at runtime by an external process (digital twin) writing to a control file.

module io_throttle

  use constants, only: CUSTOM_REAL, MAX_STRING_LEN, IMAIN

  implicit none

  ! =====================================================================
  ! I/O Throttle Parameters (copied from Par_file via shared_input_parameters)
  ! These are local copies that can be modified at runtime
  ! =====================================================================

  ! Pre-checkpoint delay in seconds (enables staggering across simultaneous runs)
  double precision :: IO_PRE_CHECKPOINT_DELAY_SEC = 0.0d0

  ! Maximum I/O bandwidth in MB/s (0 = unlimited)
  integer :: IO_MAX_BANDWIDTH_MBPS = 0

  ! Enable adaptive throttling (dynamic rate adjustment based on system load)
  logical :: IO_ADAPTIVE_THROTTLE = .false.

  ! Enable runtime throttle control via control file
  logical :: IO_RUNTIME_THROTTLE_CONTROL = .false.

  ! Output files directory (copied from shared_input_parameters)
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES_DIR = ''

  ! =====================================================================
  ! Runtime throttle state
  ! =====================================================================

  ! Current throttle values (may be modified at runtime)
  double precision :: current_delay_sec = 0.0d0
  integer :: current_max_bandwidth_mbps = 0

  ! Statistics
  integer :: throttle_apply_count = 0
  double precision :: total_delay_applied = 0.0d0

  ! Control file path
  character(len=MAX_STRING_LEN) :: throttle_control_file = ''

  ! Initialization flag
  logical :: io_throttle_initialized = .false.

  private

  ! Public subroutines
  public :: io_throttle_init
  public :: io_throttle_apply_pre_checkpoint_delay
  public :: io_throttle_read_control_file
  public :: io_throttle_finalize
  public :: io_throttle_log_stats

  ! Public variables (for external access if needed)
  public :: IO_PRE_CHECKPOINT_DELAY_SEC
  public :: IO_MAX_BANDWIDTH_MBPS
  public :: IO_ADAPTIVE_THROTTLE
  public :: IO_RUNTIME_THROTTLE_CONTROL
  public :: OUTPUT_FILES_DIR
  public :: current_delay_sec
  public :: current_max_bandwidth_mbps

contains

  !=====================================================================
  !> Initialize I/O throttle module
  !!
  !! Call this once during solver initialization to set up throttle parameters
  !! from Par_file values and prepare the control file path.
  !=====================================================================
  subroutine io_throttle_init(myrank)

    implicit none

    integer, intent(in) :: myrank

    ! Set current values from Par_file parameters
    current_delay_sec = IO_PRE_CHECKPOINT_DELAY_SEC
    current_max_bandwidth_mbps = IO_MAX_BANDWIDTH_MBPS

    ! Set up control file path
    throttle_control_file = trim(OUTPUT_FILES_DIR) // '/io_throttle_control.txt'

    ! Reset statistics
    throttle_apply_count = 0
    total_delay_applied = 0.0d0

    ! Mark as initialized
    io_throttle_initialized = .true.

    ! User output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'I/O throttle initialization:'
      write(IMAIN,*) '  IO_PRE_CHECKPOINT_DELAY_SEC   = ', IO_PRE_CHECKPOINT_DELAY_SEC
      write(IMAIN,*) '  IO_MAX_BANDWIDTH_MBPS         = ', IO_MAX_BANDWIDTH_MBPS
      write(IMAIN,*) '  IO_ADAPTIVE_THROTTLE          = ', IO_ADAPTIVE_THROTTLE
      write(IMAIN,*) '  IO_RUNTIME_THROTTLE_CONTROL   = ', IO_RUNTIME_THROTTLE_CONTROL
      if (IO_RUNTIME_THROTTLE_CONTROL) then
        write(IMAIN,*) '  control file: ', trim(throttle_control_file)
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  end subroutine io_throttle_init

  !=====================================================================
  !> Apply pre-checkpoint delay
  !!
  !! Call this before checkpoint writes to apply the configured delay.
  !! If runtime throttle control is enabled, reads the control file first.
  !!
  !! @param myrank MPI rank for logging
  !! @param iteration_on_subset Current subset iteration (for logging)
  !=====================================================================
  subroutine io_throttle_apply_pre_checkpoint_delay(myrank, iteration_on_subset)

    implicit none

    integer, intent(in) :: myrank
    integer, intent(in) :: iteration_on_subset

    integer :: delay_microseconds

    ! Check if initialized
    if (.not. io_throttle_initialized) return

    ! Read control file if runtime control is enabled
    if (IO_RUNTIME_THROTTLE_CONTROL) then
      call io_throttle_read_control_file(myrank)
    endif

    ! Apply delay if configured
    if (current_delay_sec > 0.0d0) then
      ! Convert seconds to microseconds for usleep
      delay_microseconds = int(current_delay_sec * 1.0d6)

      ! Log the delay
      if (myrank == 0) then
        write(IMAIN,'(a,i6,a,f8.4,a)') '  I/O throttle: subset ', iteration_on_subset, &
                                       ' applying pre-checkpoint delay of ', current_delay_sec, ' sec'
        call flush_IMAIN()
      endif

      ! Apply the delay using sleep
      call io_throttle_sleep(current_delay_sec)

      ! Update statistics
      throttle_apply_count = throttle_apply_count + 1
      total_delay_applied = total_delay_applied + current_delay_sec
    endif

  end subroutine io_throttle_apply_pre_checkpoint_delay

  !=====================================================================
  !> Read runtime throttle control file
  !!
  !! Parses OUTPUT_FILES/io_throttle_control.txt for runtime parameter updates.
  !! Format: key=value pairs, one per line
  !!
  !! Supported keys:
  !!   IO_PRE_CHECKPOINT_DELAY_SEC=<float>
  !!   IO_MAX_BANDWIDTH_MBPS=<int>
  !!
  !! @param myrank MPI rank for logging
  !=====================================================================
  subroutine io_throttle_read_control_file(myrank)

    implicit none

    integer, intent(in) :: myrank

    logical :: file_exists
    integer :: ier, iunit
    character(len=256) :: line
    character(len=64) :: key
    character(len=128) :: value_str
    integer :: eq_pos
    double precision :: new_delay
    integer :: new_bandwidth
    logical :: values_changed

    ! Check if file exists
    inquire(file=trim(throttle_control_file), exist=file_exists)
    if (.not. file_exists) return

    ! Open control file
    iunit = 99  ! Use a high unit number
    open(unit=iunit, file=trim(throttle_control_file), status='old', &
         action='read', iostat=ier)
    if (ier /= 0) return

    values_changed = .false.
    new_delay = current_delay_sec
    new_bandwidth = current_max_bandwidth_mbps

    ! Read and parse each line
    do
      read(iunit, '(a)', iostat=ier) line
      if (ier /= 0) exit

      ! Skip empty lines and comments
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      ! Find the equals sign
      eq_pos = index(line, '=')
      if (eq_pos <= 0) cycle

      ! Extract key and value
      key = adjustl(line(1:eq_pos-1))
      value_str = adjustl(line(eq_pos+1:))

      ! Parse known keys
      select case (trim(key))
        case ('IO_PRE_CHECKPOINT_DELAY_SEC')
          read(value_str, *, iostat=ier) new_delay
          if (ier == 0 .and. new_delay /= current_delay_sec) then
            values_changed = .true.
          endif

        case ('IO_MAX_BANDWIDTH_MBPS')
          read(value_str, *, iostat=ier) new_bandwidth
          if (ier == 0 .and. new_bandwidth /= current_max_bandwidth_mbps) then
            values_changed = .true.
          endif
      end select
    enddo

    close(iunit)

    ! Apply new values if changed
    if (values_changed) then
      current_delay_sec = new_delay
      current_max_bandwidth_mbps = new_bandwidth

      if (myrank == 0) then
        write(IMAIN,*) '  I/O throttle: updated from control file'
        write(IMAIN,*) '    new delay_sec     = ', current_delay_sec
        write(IMAIN,*) '    new bandwidth_mbps = ', current_max_bandwidth_mbps
        call flush_IMAIN()
      endif
    endif

  end subroutine io_throttle_read_control_file

  !=====================================================================
  !> Sleep for the specified number of seconds
  !!
  !! Uses a portable approach with Fortran intrinsics
  !!
  !! @param seconds Number of seconds to sleep (can be fractional)
  !=====================================================================
  subroutine io_throttle_sleep(seconds)

    implicit none

    double precision, intent(in) :: seconds

    integer :: sleep_seconds
    double precision :: remaining

    ! Use integer sleep for whole seconds
    sleep_seconds = int(seconds)
    remaining = seconds - dble(sleep_seconds)

    ! Call system sleep for whole seconds
    if (sleep_seconds > 0) then
      call sleep(sleep_seconds)
    endif

    ! For sub-second delays, use a busy wait (not ideal but portable)
    ! In production, this could be replaced with usleep via C interop
    if (remaining > 0.0d0) then
      call io_throttle_busy_wait(remaining)
    endif

  end subroutine io_throttle_sleep

  !=====================================================================
  !> Busy wait for a fraction of a second
  !!
  !! @param seconds Fractional seconds to wait (should be < 1.0)
  !=====================================================================
  subroutine io_throttle_busy_wait(seconds)

    implicit none

    double precision, intent(in) :: seconds

    double precision :: start_time, current_time
    double precision, external :: wtime

    start_time = wtime()
    do
      current_time = wtime()
      if (current_time - start_time >= seconds) exit
    enddo

  end subroutine io_throttle_busy_wait

  !=====================================================================
  !> Log throttle statistics
  !!
  !! Call at the end of the simulation to report throttle usage.
  !!
  !! @param myrank MPI rank for logging
  !=====================================================================
  subroutine io_throttle_log_stats(myrank)

    implicit none

    integer, intent(in) :: myrank

    if (.not. io_throttle_initialized) return

    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'I/O throttle statistics:'
      write(IMAIN,*) '  total delay applications = ', throttle_apply_count
      write(IMAIN,*) '  total delay time (sec)   = ', total_delay_applied
      if (throttle_apply_count > 0) then
        write(IMAIN,*) '  average delay (sec)      = ', total_delay_applied / dble(throttle_apply_count)
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  end subroutine io_throttle_log_stats

  !=====================================================================
  !> Finalize I/O throttle module
  !!
  !! Call at the end of the simulation for cleanup.
  !!
  !! @param myrank MPI rank for logging
  !=====================================================================
  subroutine io_throttle_finalize(myrank)

    implicit none

    integer, intent(in) :: myrank

    ! Log statistics
    call io_throttle_log_stats(myrank)

    ! Reset state
    io_throttle_initialized = .false.

  end subroutine io_throttle_finalize

end module io_throttle