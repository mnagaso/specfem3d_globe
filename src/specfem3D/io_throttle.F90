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
!! 5. Coordinated barriers (IO_COORDINATED_CHECKPOINT_BARRIERS)
!!    - MPI_COMM_WORLD barrier before/after each checkpoint across simultaneous runs
!!    - Resets wall-clock drift so adjacent groups stay ~n seconds apart every CKPT
!!
!! 6. Arrival-order reservation (IO_STAGGER_MODE=arrival_order)
!!    - Barrier-free shared-memory scheduler on /dev/shm
!!    - Reserves write-start slots in actual CKPT arrival order
!!
!! The throttle parameters can be set via Par_file or dynamically adjusted
!! at runtime by an external process (digital twin) writing to a control file.

module io_throttle

  use constants, only: CUSTOM_REAL, MAX_STRING_LEN, IMAIN
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS
  use, intrinsic :: iso_c_binding

  implicit none

  interface
    subroutine synchronize_all()
    end subroutine synchronize_all
    subroutine synchronize_all_world()
    end subroutine synchronize_all_world
    subroutine world_require_single_hostname()
    end subroutine world_require_single_hostname
    subroutine world_bcast_ch(buffer, countval)
      character(len=*) :: buffer
      integer :: countval
    end subroutine world_bcast_ch
    subroutine world_rank_and_size(rank, nprocs)
      integer :: rank, nprocs
    end subroutine world_rank_and_size
    subroutine bcast_all_singledp(buffer)
      double precision :: buffer
    end subroutine bcast_all_singledp
    subroutine bcast_all_singlei(buffer)
      integer :: buffer
    end subroutine bcast_all_singlei
    subroutine bcast_all_singlei8(buffer)
      integer(kind=8) :: buffer
    end subroutine bcast_all_singlei8
    double precision function wtime()
    end function wtime
  end interface

  interface
    function aos_init_c(run_key, state_dir, num_groups, world_rank, create_new) &
         bind(C, name="aos_init")
      import :: c_char, c_int
      integer(c_int) :: aos_init_c
      character(kind=c_char), intent(in) :: run_key(*)
      character(kind=c_char), intent(in) :: state_dir(*)
      integer(c_int), value :: num_groups
      integer(c_int), value :: world_rank
      integer(c_int), value :: create_new
    end function aos_init_c

    function aos_reserve_c(checkpoint_id, group_id, spacing_sec, arrival_ns, &
                           scheduled_start_ns, arrival_sequence, assigned_delay_sec) &
         bind(C, name="aos_reserve")
      import :: c_int, c_double, c_int64_t
      integer(c_int) :: aos_reserve_c
      integer(c_int), value :: checkpoint_id
      integer(c_int), value :: group_id
      real(c_double), value :: spacing_sec
      integer(c_int64_t) :: arrival_ns
      integer(c_int64_t) :: scheduled_start_ns
      integer(c_int) :: arrival_sequence
      real(c_double) :: assigned_delay_sec
    end function aos_reserve_c

    function aos_sleep_until_c(scheduled_start_ns) bind(C, name="aos_sleep_until")
      import :: c_int, c_int64_t
      integer(c_int) :: aos_sleep_until_c
      integer(c_int64_t), value :: scheduled_start_ns
    end function aos_sleep_until_c

    function aos_finalize_c(run_key, state_dir, world_rank) bind(C, name="aos_finalize")
      import :: c_char, c_int
      integer(c_int) :: aos_finalize_c
      character(kind=c_char), intent(in) :: run_key(*)
      character(kind=c_char), intent(in) :: state_dir(*)
      integer(c_int), value :: world_rank
    end function aos_finalize_c

    function aos_now_ns_c() bind(C, name="aos_now_ns")
      import :: c_int64_t
      integer(c_int64_t) :: aos_now_ns_c
    end function aos_now_ns_c

    function aos_last_error_c() bind(C, name="aos_last_error")
      import :: c_ptr
      type(c_ptr) :: aos_last_error_c
    end function aos_last_error_c
  end interface

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

  ! Coordinated checkpoint barriers across simultaneous runs
  logical :: IO_COORDINATED_CHECKPOINT_BARRIERS = .false.

  ! Spacing between adjacent run groups (0 = use IO_PRE_CHECKPOINT_DELAY_SEC)
  double precision :: IO_CHECKPOINT_SPACING_SEC = 0.0d0

  ! Runtime-selectable stagger mode
  character(len=MAX_STRING_LEN) :: IO_STAGGER_MODE = 'fixed_group'
  character(len=MAX_STRING_LEN) :: IO_ARRIVAL_ORDER_FALLBACK = 'abort'
  character(len=MAX_STRING_LEN) :: IO_ARRIVAL_ORDER_STATE_DIR = '/dev/shm'

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

  ! Effective stagger mode string used at runtime
  character(len=32) :: active_stagger_mode = 'fixed_group'

  ! Arrival-order scheduler state
  logical :: aos_enabled = .false.
  logical :: aos_initialized = .false.
  character(len=MAX_STRING_LEN) :: aos_run_key = ''
  character(len=MAX_STRING_LEN) :: aos_state_dir = '/dev/shm'
  character(len=MAX_STRING_LEN) :: aos_fallback = 'abort'
  character(len=MAX_STRING_LEN) :: aos_timing_csv = ''

  ! Last reservation / write timing (per group, for CSV)
  integer :: last_checkpoint_id = -1
  integer :: last_arrival_sequence = -1
  integer(kind=8) :: last_arrival_ns = 0
  integer(kind=8) :: last_scheduled_start_ns = 0
  integer(kind=8) :: last_write_start_ns = 0
  integer(kind=8) :: last_write_end_ns = 0
  double precision :: last_assigned_delay_sec = 0.0d0
  double precision :: last_target_spacing_sec = 0.0d0
  double precision :: last_local_sync_sec = 0.0d0
  double precision :: last_bytes_written = 0.0d0
  logical :: aos_csv_header_written = .false.

  ! Statistics for pre-checkpoint delay
  integer :: throttle_apply_count = 0
  double precision :: total_delay_applied = 0.0d0

  ! Statistics for post-write bandwidth enforcement
  integer :: bw_enforce_count = 0
  double precision :: total_bw_delay_applied = 0.0d0

  ! Statistics for coordinated barriers
  integer :: coord_barrier_pre_count = 0
  integer :: coord_barrier_post_count = 0

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
  public :: io_throttle_coordinated_post_checkpoint
  public :: io_throttle_record_write_timing
  public :: io_throttle_maybe_fsync_after_checkpoint
  public :: io_throttle_read_control_file
  public :: io_throttle_finalize
  public :: io_throttle_log_stats

  ! Public variables (for external access if needed)
  public :: IO_PRE_CHECKPOINT_DELAY_SEC
  public :: IO_MAX_BANDWIDTH_MBPS
  public :: IO_ADAPTIVE_THROTTLE
  public :: IO_RUNTIME_THROTTLE_CONTROL
  public :: IO_COORDINATED_CHECKPOINT_BARRIERS
  public :: IO_CHECKPOINT_SPACING_SEC
  public :: IO_STAGGER_MODE
  public :: IO_ARRIVAL_ORDER_FALLBACK
  public :: IO_ARRIVAL_ORDER_STATE_DIR
  public :: OUTPUT_FILES_DIR
  public :: current_delay_sec
  public :: current_max_bandwidth_mbps
  public :: my_run_group
  public :: active_stagger_mode

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

    integer :: world_rank, ier_env, aos_rc, world_nprocs
    character(len=64) :: pid_str, time_str
    character(kind=c_char), dimension(MAX_STRING_LEN+1) :: c_run_key, c_state_dir
    integer :: i

    ! Store group index
    my_run_group = mygroup

    ! Set current values from Par_file parameters
    current_delay_sec = IO_PRE_CHECKPOINT_DELAY_SEC
    current_max_bandwidth_mbps = IO_MAX_BANDWIDTH_MBPS

    ! Resolve stagger mode (legacy flags remain supported)
    call io_throttle_resolve_mode()

    ! Set up control file path
    throttle_control_file = trim(OUTPUT_FILES_DIR) // '/io_throttle_control.txt'
    aos_timing_csv = trim(OUTPUT_FILES_DIR) // '/arrival_order_timing.csv'
    aos_csv_header_written = .false.

    ! Reset statistics
    throttle_apply_count = 0
    total_delay_applied = 0.0d0
    bw_enforce_count = 0
    total_bw_delay_applied = 0.0d0
    coord_barrier_pre_count = 0
    coord_barrier_post_count = 0

    ! Reset adaptive state
    ewma_write_time = 0.0d0
    baseline_write_time = 0.0d0
    adaptive_obs_count = 0
    adaptive_delay_sec = 0.0d0

    ! Arrival-order setup
    aos_enabled = (trim(active_stagger_mode) == 'arrival_order')
    aos_initialized = .false.
    aos_fallback = trim(IO_ARRIVAL_ORDER_FALLBACK)
    if (len_trim(aos_fallback) == 0) aos_fallback = 'abort'
    aos_state_dir = trim(IO_ARRIVAL_ORDER_STATE_DIR)
    if (len_trim(aos_state_dir) == 0) aos_state_dir = '/dev/shm'

    call world_rank_and_size(world_rank, world_nprocs)

    ! Shared run key for all groups (required for posix shm name)
    call get_environment_variable('SPECFEM_AOS_RUN_KEY', aos_run_key, status=ier_env)
    if (world_rank == 0) then
      if (ier_env /= 0 .or. len_trim(aos_run_key) == 0) then
        write(pid_str, '(i0)') getpid_compat()
        write(time_str, '(i0)') int(wtime() * 1.0d6)
        aos_run_key = 'aos_' // trim(pid_str) // '_' // trim(time_str)
      endif
    endif
    call world_bcast_ch(aos_run_key, len(aos_run_key))

    if (aos_enabled) then
      ! Single-node only
      call world_require_single_hostname()

      if (index(aos_state_dir, '/dev/shm') /= 1) then
        if (myrank == 0) then
          write(IMAIN,*) 'ERROR: IO_ARRIVAL_ORDER_STATE_DIR must be under /dev/shm'
          write(IMAIN,*) '  got: ', trim(aos_state_dir)
          call flush_IMAIN()
        endif
        stop 'Arrival-order state dir must be /dev/shm'
      endif

      call io_throttle_to_c_string(aos_run_key, c_run_key)
      call io_throttle_to_c_string(aos_state_dir, c_state_dir)

      aos_rc = aos_init_c(c_run_key, c_state_dir, &
                          int(max(NUMBER_OF_SIMULTANEOUS_RUNS, 1), c_int), &
                          int(world_rank, c_int), &
                          merge(1_c_int, 0_c_int, world_rank == 0))

      if (aos_rc /= 0) then
        call io_throttle_handle_aos_failure(myrank, 'aos_init', aos_rc)
      else
        aos_initialized = .true.
      endif
    endif

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
      write(IMAIN,*) '  IO_COORDINATED_CHECKPOINT_BARRIERS = ', IO_COORDINATED_CHECKPOINT_BARRIERS
      write(IMAIN,*) '  IO_CHECKPOINT_SPACING_SEC     = ', IO_CHECKPOINT_SPACING_SEC
      write(IMAIN,*) '  active_stagger_mode           = ', trim(active_stagger_mode)
      write(IMAIN,*) '  my_run_group (0-based)        = ', my_run_group
      if (trim(active_stagger_mode) == 'arrival_order') then
        write(IMAIN,*) '  Scheduler backend: posix_shm'
        write(IMAIN,*) '  Scheduler state: ', trim(aos_state_dir), '/specfem_aos_', trim(aos_run_key)
        write(IMAIN,*) '  Fallback: ', trim(aos_fallback)
        write(IMAIN,*) '  timing csv: ', trim(aos_timing_csv)
      else if (my_run_group >= 0) then
        if (trim(active_stagger_mode) == 'coordinated') then
          write(IMAIN,*) '  coordinated stagger for group  = ', &
                         dble(my_run_group) * io_throttle_get_spacing_sec(), ' sec'
        else if (trim(active_stagger_mode) == 'fixed_group') then
          write(IMAIN,*) '  staggered delay for this group = ', &
                         dble(my_run_group) * IO_PRE_CHECKPOINT_DELAY_SEC, ' sec'
        endif
      endif
      if (IO_COORDINATED_CHECKPOINT_BARRIERS .and. IO_ADAPTIVE_THROTTLE) then
        write(IMAIN,*) '  WARNING: adaptive throttle with coordinated barriers may break spacing'
      endif
      if (IO_RUNTIME_THROTTLE_CONTROL) then
        write(IMAIN,*) '  control file: ', trim(throttle_control_file)
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! silence unused
    i = 0

  end subroutine io_throttle_init

  !=====================================================================
  subroutine io_throttle_resolve_mode()

    implicit none

    character(len=MAX_STRING_LEN) :: mode

    mode = adjustl(IO_STAGGER_MODE)
    if (len_trim(mode) == 0) then
      if (IO_COORDINATED_CHECKPOINT_BARRIERS) then
        mode = 'coordinated'
      else
        mode = 'fixed_group'
      endif
    endif

    select case (trim(mode))
      case ('none', 'no_delay', 'off')
        active_stagger_mode = 'none'
        IO_COORDINATED_CHECKPOINT_BARRIERS = .false.
      case ('fixed_group', 'fixed')
        active_stagger_mode = 'fixed_group'
        IO_COORDINATED_CHECKPOINT_BARRIERS = .false.
      case ('arrival_order', 'arrival')
        active_stagger_mode = 'arrival_order'
        IO_COORDINATED_CHECKPOINT_BARRIERS = .false.
      case ('coordinated', 'coord')
        active_stagger_mode = 'coordinated'
        IO_COORDINATED_CHECKPOINT_BARRIERS = .true.
      case default
        write(*,*) 'Unknown IO_STAGGER_MODE: ', trim(mode)
        stop 'Invalid IO_STAGGER_MODE'
    end select

  end subroutine io_throttle_resolve_mode

  !=====================================================================
  integer function getpid_compat()
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    interface
      function c_getpid() bind(C, name="getpid")
        import :: c_int
        integer(c_int) :: c_getpid
      end function c_getpid
    end interface
    getpid_compat = int(c_getpid())
  end function getpid_compat

  !=====================================================================
  subroutine io_throttle_to_c_string(fstr, cstr)
    implicit none
    character(len=*), intent(in) :: fstr
    character(kind=c_char), intent(out) :: cstr(*)
    integer :: i, n
    n = len_trim(fstr)
    do i = 1, n
      cstr(i) = fstr(i:i)
    enddo
    cstr(n+1) = c_null_char
  end subroutine io_throttle_to_c_string

  !=====================================================================
  subroutine io_throttle_handle_aos_failure(myrank, where, rc)
    implicit none
    integer, intent(in) :: myrank, rc
    character(len=*), intent(in) :: where
    type(c_ptr) :: err_ptr
    character(kind=c_char), pointer :: err_chars(:)
    character(len=512) :: err_msg
    integer :: i

    if (myrank == 0) then
      err_msg = ''
      err_ptr = aos_last_error_c()
      if (c_associated(err_ptr)) then
        call c_f_pointer(err_ptr, err_chars, [512])
        do i = 1, 512
          if (err_chars(i) == c_null_char) exit
          err_msg(i:i) = err_chars(i)
        enddo
      endif
      write(IMAIN,*) 'ERROR: arrival-order scheduler failure in ', trim(where), ' rc=', rc
      if (len_trim(err_msg) > 0) then
        write(IMAIN,*) '  aos_last_error: ', trim(err_msg)
      endif
      call flush_IMAIN()
    endif

    select case (trim(aos_fallback))
      case ('abort')
        stop 'Arrival-order scheduler failure (fallback=abort)'
      case ('fixed_group')
        if (myrank == 0) write(IMAIN,*) '  Falling back to fixed_group (explicit)'
        active_stagger_mode = 'fixed_group'
        aos_enabled = .false.
        aos_initialized = .false.
      case ('none', 'no_delay')
        if (myrank == 0) write(IMAIN,*) '  Falling back to no delay (explicit)'
        active_stagger_mode = 'none'
        aos_enabled = .false.
        aos_initialized = .false.
      case default
        stop 'Arrival-order scheduler failure: unknown fallback (silent fallback forbidden)'
    end select
  end subroutine io_throttle_handle_aos_failure

  !=====================================================================
  !> Apply pre-checkpoint delay (staggered by group index or arrival-order)
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

    last_checkpoint_id = iteration_on_subset
    last_target_spacing_sec = io_throttle_get_spacing_sec()

    if (trim(active_stagger_mode) == 'none') then
      last_assigned_delay_sec = 0.0d0
      last_arrival_sequence = -1
      last_local_sync_sec = 0.0d0
      return
    endif

    if (trim(active_stagger_mode) == 'coordinated' .and. &
        NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. my_run_group >= 0) then
      call io_throttle_coordinated_pre_checkpoint(myrank, iteration_on_subset)
      return
    endif

    if (trim(active_stagger_mode) == 'arrival_order' .and. aos_enabled) then
      call io_throttle_arrival_order_pre_checkpoint(myrank, iteration_on_subset)
      return
    endif

    ! fixed_group baseline
    if (my_run_group > 0) then
      staggered_delay = dble(my_run_group) * current_delay_sec
    else if (my_run_group == 0) then
      staggered_delay = 0.0d0
    else
      staggered_delay = current_delay_sec
    endif

    total_pre_delay = staggered_delay + adaptive_delay_sec
    last_assigned_delay_sec = total_pre_delay
    last_arrival_sequence = my_run_group
    last_local_sync_sec = 0.0d0

    if (total_pre_delay > 0.0d0) then
      if (myrank == 0) then
        write(IMAIN,'(a,i6,a,i3,a,f8.4,a,f8.4,a)') &
          '  I/O throttle: subset ', iteration_on_subset, &
          ' group ', my_run_group, &
          ' stagger=', staggered_delay, &
          ' adaptive=', adaptive_delay_sec, ' sec'
        call flush_IMAIN()
      endif

      call io_throttle_sleep(total_pre_delay)

      throttle_apply_count = throttle_apply_count + 1
      total_delay_applied = total_delay_applied + total_pre_delay
    endif

  end subroutine io_throttle_apply_pre_checkpoint_delay

  !=====================================================================
  !> Arrival-order reservation + absolute sleep (no MPI_COMM_WORLD barrier)
  !=====================================================================
  subroutine io_throttle_arrival_order_pre_checkpoint(myrank, iteration_on_subset)

    implicit none

    integer, intent(in) :: myrank
    integer, intent(in) :: iteration_on_subset

    double precision :: t0, t1, spacing_sec
    integer :: aos_rc, seq
    integer(kind=8) :: arrival_ns, scheduled_ns
    double precision :: delay_sec
    integer(c_int64_t) :: c_arrival, c_scheduled
    integer(c_int) :: c_seq
    real(c_double) :: c_delay

    if (.not. aos_initialized) then
      call io_throttle_handle_aos_failure(myrank, 'pre_checkpoint_not_init', -7)
      if (trim(active_stagger_mode) /= 'arrival_order') then
        call io_throttle_apply_pre_checkpoint_delay(myrank, iteration_on_subset)
      endif
      return
    endif

    spacing_sec = io_throttle_get_spacing_sec()
    last_target_spacing_sec = spacing_sec

    ! Ensure all ranks in this simulation group are ready before reserving
    t0 = wtime()
    call synchronize_all()
    t1 = wtime()
    last_local_sync_sec = t1 - t0

    arrival_ns = 0
    scheduled_ns = 0
    seq = 0
    delay_sec = 0.0d0

    if (myrank == 0) then
      c_arrival = 0_c_int64_t
      c_scheduled = 0_c_int64_t
      c_seq = 0_c_int
      c_delay = 0.0_c_double
      aos_rc = aos_reserve_c(int(iteration_on_subset, c_int), &
                             int(max(my_run_group, 0), c_int), &
                             real(spacing_sec, c_double), &
                             c_arrival, c_scheduled, c_seq, c_delay)
      if (aos_rc /= 0) then
        call io_throttle_handle_aos_failure(myrank, 'aos_reserve', aos_rc)
        if (trim(active_stagger_mode) /= 'arrival_order') then
          ! Recurse into non-arrival path after explicit fallback
          call io_throttle_apply_pre_checkpoint_delay(myrank, iteration_on_subset)
          return
        endif
      endif
      arrival_ns = int(c_arrival, kind=8)
      scheduled_ns = int(c_scheduled, kind=8)
      seq = int(c_seq)
      delay_sec = dble(c_delay)
    endif

    ! Broadcast reservation to all ranks in this simulation group
    call bcast_all_singlei8(arrival_ns)
    call bcast_all_singlei8(scheduled_ns)
    call bcast_all_singlei(seq)
    call bcast_all_singledp(delay_sec)

    last_arrival_ns = arrival_ns
    last_scheduled_start_ns = scheduled_ns
    last_arrival_sequence = seq
    last_assigned_delay_sec = delay_sec

    if (myrank == 0) then
      write(IMAIN,'(a,i6,a,i3,a,i3,a,f8.4,a,f8.4,a)') &
        '  I/O throttle arrival: subset ', iteration_on_subset, &
        ' group ', my_run_group, &
        ' seq=', seq, &
        ' spacing=', spacing_sec, &
        ' delay=', delay_sec, ' sec'
      call flush_IMAIN()
    endif

    aos_rc = aos_sleep_until_c(int(scheduled_ns, c_int64_t))
    if (aos_rc /= 0) then
      call io_throttle_handle_aos_failure(myrank, 'aos_sleep_until', aos_rc)
    endif

    if (delay_sec > 0.0d0) then
      throttle_apply_count = throttle_apply_count + 1
      total_delay_applied = total_delay_applied + delay_sec
    endif

  end subroutine io_throttle_arrival_order_pre_checkpoint

  !=====================================================================
  !> Optional durable flush after checkpoint write (causal experiment).
  !!
  !! Enabled when env SPECFEM_IO_FSYNC_AFTER_CHECKPOINT is 1/true/.true.
  !! Uses libc sync() so dirty pages for the filesystem are pushed before
  !! compute resumes. Cost lands in timing_acct CKPT io (outside write_elapsed).
  !! Only myrank==0 issues sync to avoid N-way stampedes; other ranks wait
  !! at the next MPI sync points of the solver.
  !=====================================================================
  subroutine io_throttle_maybe_fsync_after_checkpoint(myrank)

    implicit none

    integer, intent(in) :: myrank

    interface
      subroutine c_sync() bind(C, name="sync")
      end subroutine c_sync
    end interface

    character(len=64) :: env
    integer :: status
    double precision :: t0, t1
    double precision, external :: wtime
    logical, save :: logged_once = .false.
    logical :: enabled

    env = ''
    call get_environment_variable('SPECFEM_IO_FSYNC_AFTER_CHECKPOINT', env, status=status)
    enabled = .false.
    if (status == 0) then
      env = adjustl(env)
      if (trim(env) == '1' .or. trim(env) == 'true' .or. trim(env) == '.true.' &
          .or. trim(env) == 'TRUE' .or. trim(env) == 'yes') then
        enabled = .true.
      endif
    endif
    if (.not. enabled) return

    ! Only group-local rank0 issues sync() (4 calls with nsim=4, not 96).
    if (myrank /= 0) return

    t0 = wtime()
    call c_sync()
    t1 = wtime()

    if (.not. logged_once) then
      write(*,*) 'IO_FSYNC_AFTER_CHECKPOINT: sync() enabled (first call took ', &
                 t1 - t0, ' s)'
      logged_once = .true.
    endif

  end subroutine io_throttle_maybe_fsync_after_checkpoint

  !=====================================================================
  !> Record actual write start/end into arrival_order_timing.csv (rank 0)
  !=====================================================================
  subroutine io_throttle_record_write_timing(myrank, bytes_written, write_elapsed)

    implicit none

    integer, intent(in) :: myrank
    double precision, intent(in) :: bytes_written
    double precision, intent(in) :: write_elapsed

    integer :: iunit, ier
    integer(kind=8) :: now_ns
    character(len=64) :: host
    logical :: exist

    if (.not. io_throttle_initialized) return
    if (myrank /= 0) return

    last_bytes_written = bytes_written
    now_ns = int(aos_now_ns_c(), kind=8)

    ! Approximate write window using CLOCK_MONOTONIC around the measured elapsed
    last_write_end_ns = now_ns
    if (write_elapsed > 0.0d0) then
      last_write_start_ns = now_ns - int(write_elapsed * 1.0d9, kind=8)
    else
      last_write_start_ns = now_ns
    endif

    ! Always record for arrival_order; also useful for other modes when enabled
    if (trim(active_stagger_mode) /= 'arrival_order' .and. &
        trim(active_stagger_mode) /= 'fixed_group' .and. &
        trim(active_stagger_mode) /= 'none' .and. &
        trim(active_stagger_mode) /= 'coordinated') return

    inquire(file=trim(aos_timing_csv), exist=exist)
    iunit = 97
    if (.not. exist .or. .not. aos_csv_header_written) then
      open(unit=iunit, file=trim(aos_timing_csv), status='unknown', &
           action='write', position='rewind', iostat=ier)
      if (ier /= 0) return
      write(iunit,'(a)') &
        'run_id,checkpoint_id,group_id,arrival_sequence,'// &
        'arrival_time_ns,scheduled_start_time_ns,'// &
        'actual_write_start_time_ns,actual_write_end_time_ns,'// &
        'assigned_delay_sec,target_spacing_sec,scheduler_mode,'// &
        'local_sync_time_sec,write_elapsed_sec,bytes_written,hostname'
      aos_csv_header_written = .true.
    else
      open(unit=iunit, file=trim(aos_timing_csv), status='old', &
           action='write', position='append', iostat=ier)
      if (ier /= 0) return
    endif

    host = ''
    call get_environment_variable('HOSTNAME', host)

    write(iunit, &
      '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",es16.8,",",es16.8,",",a,",",es16.8,",",es16.8,",",es16.8,",",a)') &
      trim(aos_run_key), last_checkpoint_id, max(my_run_group, 0), last_arrival_sequence, &
      last_arrival_ns, last_scheduled_start_ns, &
      last_write_start_ns, last_write_end_ns, &
      last_assigned_delay_sec, last_target_spacing_sec, trim(active_stagger_mode), &
      last_local_sync_sec, write_elapsed, bytes_written, trim(host)

    close(iunit)

  end subroutine io_throttle_record_write_timing

  !=====================================================================
  !> Coordinated pre-checkpoint: align all groups, then stagger writes
  !=====================================================================
  subroutine io_throttle_coordinated_pre_checkpoint(myrank, iteration_on_subset)

    implicit none

    integer, intent(in) :: myrank
    integer, intent(in) :: iteration_on_subset

    double precision :: spacing_sec, staggered_delay, total_pre_delay
    double precision :: t_before, t_after

    spacing_sec = io_throttle_get_spacing_sec()

    if (my_run_group > 0) then
      staggered_delay = dble(my_run_group) * spacing_sec
    else
      staggered_delay = 0.0d0
    endif
    total_pre_delay = staggered_delay + adaptive_delay_sec

    t_before = wtime()

    ! Phase 1: align all ranks within group, then across all simultaneous runs
    call synchronize_all()
    call synchronize_all_world()

    t_after = wtime()

    if (myrank == 0) then
      write(IMAIN,'(a,i6,a,i3,a,f8.4,a,f8.4,a,f8.4,a)') &
        '  I/O throttle coord: subset ', iteration_on_subset, &
        ' group ', my_run_group, &
        ' spacing=', spacing_sec, &
        ' stagger=', staggered_delay, &
        ' barrier_wait=', t_after - t_before, ' sec'
      call flush_IMAIN()
    endif

    coord_barrier_pre_count = coord_barrier_pre_count + 1

    if (total_pre_delay > 0.0d0) then
      call io_throttle_sleep(total_pre_delay)
      throttle_apply_count = throttle_apply_count + 1
      total_delay_applied = total_delay_applied + total_pre_delay
    endif

  end subroutine io_throttle_coordinated_pre_checkpoint

  !=====================================================================
  !> Coordinated post-checkpoint: wait until all groups finish CKPT I/O
  !=====================================================================
  subroutine io_throttle_coordinated_post_checkpoint(myrank, iteration_on_subset)

    implicit none

    integer, intent(in) :: myrank
    integer, intent(in) :: iteration_on_subset

    double precision :: t_before, t_after

    if (.not. io_throttle_initialized) return
    if (.not. IO_COORDINATED_CHECKPOINT_BARRIERS) return
    if (NUMBER_OF_SIMULTANEOUS_RUNS <= 1 .or. my_run_group < 0) return

    t_before = wtime()

    call synchronize_all()
    call synchronize_all_world()

    t_after = wtime()
    coord_barrier_post_count = coord_barrier_post_count + 1

    if (myrank == 0) then
      write(IMAIN,'(a,i6,a,i3,a,f8.4,a)') &
        '  I/O throttle coord post: subset ', iteration_on_subset, &
        ' group ', my_run_group, &
        ' barrier_wait=', t_after - t_before, ' sec'
      call flush_IMAIN()
    endif

  end subroutine io_throttle_coordinated_post_checkpoint

  !=====================================================================
  !> Effective checkpoint spacing in seconds
  !=====================================================================
  function io_throttle_get_spacing_sec() result(spacing_sec)

    implicit none

    double precision :: spacing_sec

    if (IO_CHECKPOINT_SPACING_SEC > 0.0d0) then
      spacing_sec = IO_CHECKPOINT_SPACING_SEC
    else
      spacing_sec = current_delay_sec
    endif

  end function io_throttle_get_spacing_sec

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
    double precision :: new_spacing
    character(len=MAX_STRING_LEN) :: new_mode
    logical :: mode_changed
    double precision :: bare_value
    integer :: bare_ier

    ! Check if file exists
    inquire(file=trim(throttle_control_file), exist=file_exists)
    if (.not. file_exists) return

    ! Open control file
    iunit = 99  ! Use a high unit number
    open(unit=iunit, file=trim(throttle_control_file), status='old', &
         action='read', iostat=ier)
    if (ier /= 0) return

    values_changed = .false.
    mode_changed = .false.
    new_delay = current_delay_sec
    new_bandwidth = current_max_bandwidth_mbps
    new_adaptive = IO_ADAPTIVE_THROTTLE
    new_adaptive_delay = adaptive_delay_sec
    new_spacing = IO_CHECKPOINT_SPACING_SEC
    new_mode = active_stagger_mode

    ! Read and parse each line
    do
      read(iunit, '(a)', iostat=ier) line
      if (ier /= 0) exit

      ! Skip empty lines and comments
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      ! Legacy bare number => fixed_group target spacing / delay
      eq_pos = index(line, '=')
      if (eq_pos <= 0) then
        read(line, *, iostat=bare_ier) bare_value
        if (bare_ier == 0) then
          new_delay = bare_value
          new_spacing = bare_value
          new_mode = 'fixed_group'
          values_changed = .true.
          mode_changed = .true.
        endif
        cycle
      endif

      ! Extract key and value
      key = adjustl(line(1:eq_pos-1))
      value_str = adjustl(line(eq_pos+1:))

      ! Parse known keys
      select case (trim(key))
        case ('IO_PRE_CHECKPOINT_DELAY_SEC', 'target_spacing_sec', 'TARGET_SPACING_SEC')
          read(value_str, *, iostat=ier) new_delay
          if (ier == 0) then
            new_spacing = new_delay
            if (new_delay /= current_delay_sec) values_changed = .true.
          endif

        case ('IO_MAX_BANDWIDTH_MBPS')
          read(value_str, *, iostat=ier) new_bandwidth
          if (ier == 0 .and. new_bandwidth /= current_max_bandwidth_mbps) then
            values_changed = .true.
          endif

        case ('IO_ADAPTIVE_THROTTLE')
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

        case ('IO_CHECKPOINT_SPACING_SEC')
          read(value_str, *, iostat=ier) new_spacing
          if (ier == 0 .and. new_spacing /= IO_CHECKPOINT_SPACING_SEC) then
            values_changed = .true.
          endif

        case ('mode', 'IO_STAGGER_MODE')
          new_mode = adjustl(value_str)
          if (trim(new_mode) /= trim(active_stagger_mode)) then
            mode_changed = .true.
            values_changed = .true.
          endif

        case ('version', 'effective_checkpoint')
          ! accepted and ignored (metadata for Digital Twin writers)

      end select
    enddo

    close(iunit)

    ! Apply new values if changed
    if (values_changed) then
      current_delay_sec = new_delay
      current_max_bandwidth_mbps = new_bandwidth
      IO_ADAPTIVE_THROTTLE = new_adaptive
      adaptive_delay_sec = new_adaptive_delay
      IO_CHECKPOINT_SPACING_SEC = new_spacing

      if (mode_changed) then
        IO_STAGGER_MODE = new_mode
        call io_throttle_resolve_mode()
        ! Do not silently enable arrival_order without prior aos_init
        if (trim(active_stagger_mode) == 'arrival_order' .and. .not. aos_initialized) then
          if (myrank == 0) then
            write(IMAIN,*) 'ERROR: control file requested arrival_order but scheduler was not initialized'
          endif
          call io_throttle_handle_aos_failure(myrank, 'runtime_mode_arrival_order', -7)
        endif
        aos_enabled = (trim(active_stagger_mode) == 'arrival_order' .and. aos_initialized)
      endif

      if (myrank == 0) then
        write(IMAIN,*) '  I/O throttle: updated from control file'
        write(IMAIN,*) '    delay_sec / target_spacing = ', current_delay_sec
        write(IMAIN,*) '    checkpoint_spacing_sec   = ', IO_CHECKPOINT_SPACING_SEC
        write(IMAIN,*) '    stagger_mode             = ', trim(active_stagger_mode)
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
      if (IO_COORDINATED_CHECKPOINT_BARRIERS) then
        write(IMAIN,*) '  coordinated pre-barriers   = ', coord_barrier_pre_count
        write(IMAIN,*) '  coordinated post-barriers  = ', coord_barrier_post_count
      endif
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
    integer :: world_rank, world_nprocs, aos_rc
    character(kind=c_char), dimension(MAX_STRING_LEN+1) :: c_run_key, c_state_dir

    ! Log statistics
    call io_throttle_log_stats(myrank)

    if (aos_initialized) then
      call world_rank_and_size(world_rank, world_nprocs)
      call io_throttle_to_c_string(aos_run_key, c_run_key)
      call io_throttle_to_c_string(aos_state_dir, c_state_dir)
      aos_rc = aos_finalize_c(c_run_key, c_state_dir, int(world_rank, c_int))
      aos_initialized = .false.
      aos_enabled = .false.
    endif

    ! Reset state
    io_throttle_initialized = .false.

  end subroutine io_throttle_finalize

end module io_throttle