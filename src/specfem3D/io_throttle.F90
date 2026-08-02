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
!! This module provides four throttle mechanisms for simultaneous runs:
!!
!! 1. Pre-checkpoint delay (IO_PRE_CHECKPOINT_DELAY_SEC)
!!    - Staggered per simultaneous run: group i sleeps i * delta seconds
!!    - Spreads checkpoint writes across time to avoid I/O burst
!!
!! 2. Bandwidth cap (IO_MAX_BANDWIDTH_MBPS)
!!    - Enforces maximum write bandwidth per checkpoint
!!    - Adds post-write delay if writes complete too fast
!!
!! 3. Adaptive throttling (IO_ADAPTIVE_THROTTLE)
!!    - EWMA-based feedback loop on observed write times
!!    - Increases delay when contention detected, decreases when clear
!!
!! 4. Runtime control (IO_RUNTIME_THROTTLE_CONTROL)
!!    - External digital twin writes control file between checkpoints
!!    - All parameters adjustable at runtime
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

  ! Pre-checkpoint delay in seconds per group step
  ! For group i, actual delay = i * IO_PRE_CHECKPOINT_DELAY_SEC
  ! This staggers writes so simultaneous runs do not hit the FS at once
  double precision :: IO_PRE_CHECKPOINT_DELAY_SEC = 0.0d0

  ! Maximum I/O bandwidth in MB/s (0 = unlimited)
  ! When set, a post-write delay is added to enforce this bandwidth cap
  integer :: IO_MAX_BANDWIDTH_MBPS = 0

  ! Enable adaptive throttling (dynamic rate adjustment based on observed contention)
  logical :: IO_ADAPTIVE_THROTTLE = .false.

  ! Enable runtime throttle control via control file
  logical :: IO_RUNTIME_THROTTLE_CONTROL = .false.

  ! Output files directory (copied from shared_input_parameters)
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES_DIR = ''

  ! =====================================================================
  ! Runtime throttle state
  ! =====================================================================

  ! Group index for this simultaneous run (0-based; -1 if single run)
  integer :: my_run_group = -1

  ! Current throttle values (may be modified at runtime)
  double precision :: current_delay_sec = 0.0d0
  integer :: current_max_bandwidth_mbps = 0

  ! Statistics for pre-checkpoint delay
  integer :: throttle_apply_count = 0
  double precision :: total_delay_applied = 0.0d0

  ! Statistics for post-write bandwidth enforcement
  integer :: bw_enforce_count = 0
  double precision :: total_bw_delay_applied = 0.0d0

  ! =====================================================================
  ! Adaptive throttle state (EWMA-based feedback loop)
  ! =====================================================================

  ! EWMA smoothing factor (0 < alpha <= 1; higher = more responsive)
  double precision, parameter :: ADAPTIVE_ALPHA = 0.3d0

  ! Factor to scale delay adjustments
  double precision, parameter :: ADAPTIVE_INCREASE_FACTOR = 1.5d0
  double precision, parameter :: ADAPTIVE_DECREASE_FACTOR = 0.8d0

  ! Minimum adaptive delay step (seconds)
  double precision, parameter :: ADAPTIVE_MIN_STEP = 0.01d0

  ! Maximum adaptive delay (seconds) — safety cap
  double precision, parameter :: ADAPTIVE_MAX_DELAY = 10.0d0

  ! EWMA of observed write times (seconds)
  double precision :: ewma_write_time = 0.0d0

  ! Baseline write time — established from first checkpoint (no contention reference)
  double precision :: baseline_write_time = 0.0d0

  ! Number of write observations for adaptive
  integer :: adaptive_obs_count = 0

  ! Number of baseline samples to collect before starting adaptation
  integer, parameter :: ADAPTIVE_BASELINE_SAMPLES = 2

  ! Adaptive delay component (added on top of staggered delay)
  double precision :: adaptive_delay_sec = 0.0d0

  ! Control file path
  character(len=MAX_STRING_LEN) :: throttle_control_file = ''

  ! Initialization flag
  logical :: io_throttle_initialized = .false.

  private

  ! Public subroutines
  public :: io_throttle_init
  public :: io_throttle_apply_pre_checkpoint_delay
  public :: io_throttle_apply_post_write_delay
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
  public :: my_run_group

contains

  !=====================================================================
  !> Initialize I/O throttle module
  !!
  !! Call this once during solver initialization to set up throttle parameters
  !! from Par_file values and prepare the control file path.
  !!
  !! @param myrank  MPI rank for logging
  !! @param mygroup Simultaneous run group index (0-based; -1 if single run)
  !=====================================================================
  subroutine io_throttle_init(myrank, mygroup)

    implicit none

    integer, intent(in) :: myrank
    integer, intent(in) :: mygroup

    ! Store group index
    my_run_group = mygroup

    ! Set current values from Par_file parameters
    current_delay_sec = IO_PRE_CHECKPOINT_DELAY_SEC
    current_max_bandwidth_mbps = IO_MAX_BANDWIDTH_MBPS

    ! Set up control file path
    throttle_control_file = trim(OUTPUT_FILES_DIR) // '/io_throttle_control.txt'

    ! Reset statistics
    throttle_apply_count = 0
    total_delay_applied = 0.0d0
    bw_enforce_count = 0
    total_bw_delay_applied = 0.0d0

    ! Reset adaptive state
    ewma_write_time = 0.0d0
    baseline_write_time = 0.0d0
    adaptive_obs_count = 0
    adaptive_delay_sec = 0.0d0

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
      write(IMAIN,*) '  my_run_group (0-based)        = ', my_run_group
      if (my_run_group >= 0) then
        write(IMAIN,*) '  staggered delay for this group = ', &
                       dble(my_run_group) * IO_PRE_CHECKPOINT_DELAY_SEC, ' sec'
      endif
      if (IO_RUNTIME_THROTTLE_CONTROL) then
        write(IMAIN,*) '  control file: ', trim(throttle_control_file)
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  end subroutine io_throttle_init

  !=====================================================================
  !> Apply pre-checkpoint delay (staggered by group index)
  !!
  !! For simultaneous runs, group i sleeps for i * IO_PRE_CHECKPOINT_DELAY_SEC
  !! seconds before writing. Group 0 writes first, group 1 next, etc.
  !! This spreads the I/O load across time.
  !!
  !! If adaptive throttle is enabled, an additional adaptive delay is added.
  !! If runtime control is enabled, the control file is read first.
  !!
  !! @param myrank MPI rank for logging
  !! @param iteration_on_subset Current subset iteration (for logging)
  !=====================================================================
  subroutine io_throttle_apply_pre_checkpoint_delay(myrank, iteration_on_subset)

    implicit none

    integer, intent(in) :: myrank
    integer, intent(in) :: iteration_on_subset

    double precision :: staggered_delay
    double precision :: total_pre_delay

    ! Check if initialized
    if (.not. io_throttle_initialized) return

    ! Read control file if runtime control is enabled
    if (IO_RUNTIME_THROTTLE_CONTROL) then
      call io_throttle_read_control_file(myrank)
    endif

    ! Compute staggered delay based on group index
    ! Group 0: delay = 0, Group 1: delay = 1*delta, ..., Group N: delay = N*delta
    if (my_run_group > 0) then
      staggered_delay = dble(my_run_group) * current_delay_sec
    else if (my_run_group == 0) then
      ! Group 0 writes first — no stagger delay
      staggered_delay = 0.0d0
    else
      ! Single run (mygroup == -1): apply delay as-is
      staggered_delay = current_delay_sec
    endif

    ! Add adaptive delay component if enabled
    total_pre_delay = staggered_delay + adaptive_delay_sec

    ! Apply delay if > 0
    if (total_pre_delay > 0.0d0) then
      ! Log the delay
      if (myrank == 0) then
        write(IMAIN,'(a,i6,a,i3,a,f8.4,a,f8.4,a)') &
          '  I/O throttle: subset ', iteration_on_subset, &
          ' group ', my_run_group, &
          ' stagger=', staggered_delay, &
          ' adaptive=', adaptive_delay_sec, ' sec'
        call flush_IMAIN()
      endif

      ! Apply the delay using sleep
      call io_throttle_sleep(total_pre_delay)

      ! Update statistics
      throttle_apply_count = throttle_apply_count + 1
      total_delay_applied = total_delay_applied + total_pre_delay
    endif

  end subroutine io_throttle_apply_pre_checkpoint_delay

  !=====================================================================
  !> Apply post-write delay for bandwidth enforcement and adaptive feedback
  !!
  !! Called AFTER the checkpoint write completes. Performs two functions:
  !!
  !! 1. Bandwidth cap: if IO_MAX_BANDWIDTH_MBPS > 0, computes the target
  !!    write duration and sleeps for the remaining time.
  !!
  !! 2. Adaptive throttle: if IO_ADAPTIVE_THROTTLE is true, updates the
  !!    EWMA of write times and adjusts the adaptive delay for the next
  !!    checkpoint.
  !!
  !! @param myrank        MPI rank for logging
  !! @param bytes_written Total bytes written in this checkpoint (this rank)
  !! @param write_elapsed Actual wall-clock time for the write (seconds)
  !=====================================================================
  subroutine io_throttle_apply_post_write_delay(myrank, bytes_written, write_elapsed)

    implicit none

    integer, intent(in) :: myrank
    double precision, intent(in) :: bytes_written
    double precision, intent(in) :: write_elapsed

    double precision :: target_write_time, bw_delay
    double precision :: contention_ratio

    if (.not. io_throttle_initialized) return

    ! -----------------------------------------------------------------
    ! 1. Bandwidth cap enforcement
    ! -----------------------------------------------------------------
    if (current_max_bandwidth_mbps > 0 .and. bytes_written > 0.0d0) then
      ! Target time = bytes / bandwidth
      ! bytes_written is in bytes, bandwidth is in MB/s = 1e6 bytes/s
      target_write_time = bytes_written / (dble(current_max_bandwidth_mbps) * 1.0d6)

      bw_delay = target_write_time - write_elapsed
      if (bw_delay > 0.0d0) then
        if (myrank == 0) then
          write(IMAIN,'(a,f8.4,a,f8.4,a,f8.4,a)') &
            '  I/O throttle: bw cap: target=', target_write_time, &
            's actual=', write_elapsed, 's post-delay=', bw_delay, 's'
          call flush_IMAIN()
        endif

        call io_throttle_sleep(bw_delay)

        bw_enforce_count = bw_enforce_count + 1
        total_bw_delay_applied = total_bw_delay_applied + bw_delay
      endif
    endif

    ! -----------------------------------------------------------------
    ! 2. Adaptive throttle feedback
    ! -----------------------------------------------------------------
    if (IO_ADAPTIVE_THROTTLE) then
      call io_throttle_adaptive_update(myrank, write_elapsed)
    endif

  end subroutine io_throttle_apply_post_write_delay

  !=====================================================================
  !> Update adaptive throttle based on observed write time
  !!
  !! Uses EWMA to smooth write times. Compares against baseline to detect
  !! contention. Adjusts adaptive_delay_sec accordingly.
  !!
  !! Algorithm:
  !!   - First ADAPTIVE_BASELINE_SAMPLES writes establish baseline
  !!   - EWMA = alpha * new_sample + (1-alpha) * EWMA
  !!   - If EWMA > 1.5 * baseline: contention detected → increase delay
  !!   - If EWMA < 1.2 * baseline: contention easing → decrease delay
  !!   - Delay is clamped to [0, ADAPTIVE_MAX_DELAY]
  !!
  !! @param myrank       MPI rank for logging
  !! @param write_elapsed Actual write time for this checkpoint (seconds)
  !=====================================================================
  subroutine io_throttle_adaptive_update(myrank, write_elapsed)

    implicit none

    integer, intent(in) :: myrank
    double precision, intent(in) :: write_elapsed

    double precision :: old_delay
    double precision :: contention_ratio

    adaptive_obs_count = adaptive_obs_count + 1

    ! Update EWMA
    if (adaptive_obs_count == 1) then
      ewma_write_time = write_elapsed
    else
      ewma_write_time = ADAPTIVE_ALPHA * write_elapsed + &
                        (1.0d0 - ADAPTIVE_ALPHA) * ewma_write_time
    endif

    ! Collect baseline samples
    if (adaptive_obs_count <= ADAPTIVE_BASELINE_SAMPLES) then
      ! Accumulate baseline (average of first N samples)
      baseline_write_time = baseline_write_time + write_elapsed
      if (adaptive_obs_count == ADAPTIVE_BASELINE_SAMPLES) then
        baseline_write_time = baseline_write_time / dble(ADAPTIVE_BASELINE_SAMPLES)
        if (myrank == 0) then
          write(IMAIN,'(a,f8.4,a)') &
            '  I/O throttle adaptive: baseline established = ', baseline_write_time, ' sec'
          call flush_IMAIN()
        endif
      endif
      return
    endif

    ! Baseline established — adjust delay based on contention
    if (baseline_write_time <= 0.0d0) return  ! safety

    old_delay = adaptive_delay_sec

    ! Compute contention ratio = current EWMA / baseline
    contention_ratio = ewma_write_time / baseline_write_time

    if (contention_ratio > 1.5d0) then
      ! High contention: increase delay
      adaptive_delay_sec = adaptive_delay_sec + ADAPTIVE_MIN_STEP * ADAPTIVE_INCREASE_FACTOR
      if (adaptive_delay_sec > ADAPTIVE_MAX_DELAY) then
        adaptive_delay_sec = ADAPTIVE_MAX_DELAY
      endif
    else if (contention_ratio < 1.2d0) then
      ! Low contention: decrease delay
      adaptive_delay_sec = adaptive_delay_sec * ADAPTIVE_DECREASE_FACTOR
      if (adaptive_delay_sec < ADAPTIVE_MIN_STEP) then
        adaptive_delay_sec = 0.0d0
      endif
    endif
    ! Between 1.2 and 1.5: hold steady (hysteresis zone)

    if (myrank == 0 .and. adaptive_delay_sec /= old_delay) then
      write(IMAIN,'(a,f6.2,a,f8.4,a,f8.4,a)') &
        '  I/O throttle adaptive: ratio=', contention_ratio, &
        ' ewma=', ewma_write_time, ' delay=', adaptive_delay_sec, ' sec'
      call flush_IMAIN()
    endif

  end subroutine io_throttle_adaptive_update

  !=====================================================================
  !> Read runtime throttle control file
  !!
  !! Parses OUTPUT_FILES/io_throttle_control.txt for runtime parameter updates.
  !! Format: key=value pairs, one per line
  !!
  !! Supported keys:
  !!   IO_PRE_CHECKPOINT_DELAY_SEC=<float>   — base stagger step
  !!   IO_MAX_BANDWIDTH_MBPS=<int>           — bandwidth cap
  !!   IO_ADAPTIVE_THROTTLE=<true|false>     — enable/disable adaptive mode
  !!   ADAPTIVE_DELAY_SEC=<float>            — directly set adaptive delay
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
    logical :: new_adaptive
    double precision :: new_adaptive_delay
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
    new_adaptive = IO_ADAPTIVE_THROTTLE
    new_adaptive_delay = adaptive_delay_sec

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

        case ('IO_ADAPTIVE_THROTTLE')
          ! Parse boolean: .true. / true / 1 → true; else false
          value_str = adjustl(value_str)
          if (value_str(1:5) == '.true' .or. value_str(1:4) == 'true' .or. &
              value_str(1:1) == '1') then
            new_adaptive = .true.
          else
            new_adaptive = .false.
          endif
          if (new_adaptive .neqv. IO_ADAPTIVE_THROTTLE) then
            values_changed = .true.
          endif

        case ('ADAPTIVE_DELAY_SEC')
          read(value_str, *, iostat=ier) new_adaptive_delay
          if (ier == 0 .and. new_adaptive_delay /= adaptive_delay_sec) then
            values_changed = .true.
          endif
      end select
    enddo

    close(iunit)

    ! Apply new values if changed
    if (values_changed) then
      current_delay_sec = new_delay
      current_max_bandwidth_mbps = new_bandwidth
      IO_ADAPTIVE_THROTTLE = new_adaptive
      adaptive_delay_sec = new_adaptive_delay

      if (myrank == 0) then
        write(IMAIN,*) '  I/O throttle: updated from control file'
        write(IMAIN,*) '    delay_sec (base stagger) = ', current_delay_sec
        write(IMAIN,*) '    bandwidth_mbps           = ', current_max_bandwidth_mbps
        write(IMAIN,*) '    adaptive_throttle        = ', IO_ADAPTIVE_THROTTLE
        write(IMAIN,*) '    adaptive_delay_sec       = ', adaptive_delay_sec
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
      write(IMAIN,*) '  run group index              = ', my_run_group
      write(IMAIN,*) '  pre-checkpoint delays        = ', throttle_apply_count
      write(IMAIN,*) '  total pre-delay time (sec)   = ', total_delay_applied
      if (throttle_apply_count > 0) then
        write(IMAIN,*) '  average pre-delay (sec)      = ', &
                       total_delay_applied / dble(throttle_apply_count)
      endif
      write(IMAIN,*) '  bandwidth enforcement count  = ', bw_enforce_count
      write(IMAIN,*) '  total bw-delay time (sec)    = ', total_bw_delay_applied
      if (IO_ADAPTIVE_THROTTLE) then
        write(IMAIN,*) '  adaptive observations        = ', adaptive_obs_count
        write(IMAIN,*) '  final EWMA write time (sec)  = ', ewma_write_time
        write(IMAIN,*) '  baseline write time (sec)    = ', baseline_write_time
        write(IMAIN,*) '  final adaptive delay (sec)   = ', adaptive_delay_sec
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