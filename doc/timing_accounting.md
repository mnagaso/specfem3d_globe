# Per-Process Timing Accounting

## Overview

`src/specfem3D/timing_accounting.F90` provides a barrier-free, per-process wall-clock
timing breakdown for the solver time loop. Each MPI rank independently accumulates time
into four categories and writes periodic reports — one line per undo-attenuation subset
(or per `NTSTEP_BETWEEN_OUTPUT_INFO` steps for non-undoatt runs).

The design constraint is **zero synchronization barriers**: all timing is done with
`MPI_Wtime()` calls local to each process, so the measurement itself has no impact on
parallel efficiency or on the phenomena being measured.

### Recording intervals

The timing is recorded **periodically** rather than as a single accumulated total:

**`iterate_time_undoatt.F90` (undo-attenuation):**
```
subset 1 start: checkpoint IO (save_forward_arrays_undoatt)
  compute it=1..200 + movie IO     → record, reset   (includes checkpoint IO)
  compute it=201..400 + movie IO   → record, reset
  ...
  compute it=601..700 + movie IO   → record, reset
subset 2 start: checkpoint IO
  compute it=701..900 + movie IO   → record, reset   (includes checkpoint IO)
  ...
```
Recording fires every `NTSTEP_BETWEEN_OUTPUT_INFO` steps and at `it_end`.
Checkpoint IO time is included in the first interval of each subset.
Any remaining steps at subset boundaries are also recorded.

**`iterate_time.F90` (no undo-attenuation):**
```
every NTSTEP_BETWEEN_OUTPUT_INFO steps → record, reset
final step                            → record, reset
```

**IO server (`hdf5_io_server.F90`):**
```
after each undo snapshot write → record, reset
final (residual after last snapshot) → record
```

Each record line shows `it_begin` and `it_end` to identify the iteration range covered.

---

## Output files

Each rank writes to its own file to avoid any concurrent-write conflicts:

```
timing_acct_group<mygroup>_<role>_rank<myrank>.txt
```

| Token | Meaning |
|---|---|
| `<mygroup>` | Inter-communicator group index (`mygroup` from `constants` module). `-1` when `HDF5_IO_NODES = 0` (no group split). |
| `<role>` | `compute` or `io_server` |
| `<myrank>` | Rank within the role's sub-communicator (`myrank` from `constants` module). |

The files are written to the simulation working directory (same directory as
`OUTPUT_FILES/`).

> **Note on sub-communicator ranks**: when `HDF5_IO_NODES > 0`, both the
> IO server and compute rank 0 have `myrank = 0` within their respective
> sub-communicators. Including the role in the filename prevents collision.

---

## Output columns

Each file contains **one line per interval** (multiple lines per simulation):

```
it_begin: B, it_end: E, mygroup: G, myrank: R, role: ROLE, compute (s): C,
io (s): I, wait_io (s): W, idle_io (s): D, other (s): O, total (s): T,
mpi_wtime (s): E, date: DATE, time: TIME
```

### `it_begin`, `it_end`

The iteration range covered by this record.

- **Compute ranks**: timestep numbers (e.g. `it_begin: 1, it_end: 100` for
  the first subset of 100 steps).
- **IO server ranks**: undo snapshot index (e.g. `it_begin: 1, it_end: 1`
  for the first snapshot). The final residual record uses `it_begin: 0, it_end: 0`.

### `compute (s)`

Wall-clock time spent inside force computation, time stepping, and kernel evaluation.

**Source — `iterate_time.F90` (lines ~119–143, ~175–242, ~257–259, ~650–651):**
```fortran
call timing_compute_start()
do istage = 1, NSTAGE_TIME_SCHEME
  call update_displ_Newmark()       ! or LDDRK variant
  call compute_forces_acoustic()
  call compute_forces_viscoelastic()
enddo
call timing_compute_stop()
```

**Source — `iterate_time_undoatt.F90` (lines ~289–313, ~398–422, ~484–508, ~606–630, ~650–653):**
```fortran
call timing_compute_start()
do istage = 1, NSTAGE_TIME_SCHEME
  call update_displ_Newmark()
  call compute_forces_acoustic()
  call compute_forces_viscoelastic()
enddo
call timing_compute_stop()
```

Also wraps `compute_kernels()` in adjoint/kernel simulations (`SIMULATION_TYPE == 3`).

Always **zero for IO server ranks** — the IO server performs no physics.

---

### `io (s)`

**For compute ranks**: wall-clock time inside IO operations — seismogram writes,
movie output, and undo-attenuation checkpoint saves/reads.

**Source — `iterate_time_undoatt.F90` (lines ~316–327, ~350–358):**
```fortran
! undo-attenuation checkpoint save (SIMULATION_TYPE 1/2):
call timing_io_start()
call save_forward_arrays_undoatt()   ! packs wavefield arrays + MPI_Isend to IO server
call timing_io_stop()

! undo-attenuation checkpoint read (SIMULATION_TYPE 3):
call timing_io_start()
call read_forward_arrays_undoatt()
call timing_io_stop()
```

**Source — `iterate_time.F90` / `iterate_time_undoatt.F90` (seismogram + movie):**
```fortran
call timing_io_start()
call write_seismograms()
if (MOVIE_SURFACE .or. MOVIE_VOLUME) call write_movie_output()
call timing_io_stop()
```

> **Why compute's `io (s)` is larger than the IO server's `io (s)`**:
> The compute rank times the *entire* `save_forward_arrays_undoatt()` call, which
> includes packing large wavefield arrays (crust_mantle, outer_core, inner_core
> displ/veloc/accel) and posting the non-blocking `MPI_Isend`.  The IO server times
> only the `recv_and_write_ford_undo()` / `write_buffered_undo_snapshot()` calls, which
> are just writing pre-packed bytes to disk and are much cheaper.

**For IO server ranks**: wall-clock time inside receive + write operations.

**Source — `hdf5_io_server.F90` (lines ~1649, ~1667, ~1684, ~1701, ~1722):**
```fortran
! undo-attenuation: receive per-field arrays and write to HDF5 buffer
call timing_io_start()
call recv_and_write_ford_undo(tag, tag_src, status, ...)
call timing_io_stop()

! flush buffered snapshot to disk + close HDF5 file
call timing_io_start()
call write_buffered_undo_snapshot(undo_use_collective)
call h5_close_file_p()
call timing_io_stop()

! surface / volume movie frames
call timing_io_start()
call write_surface_frame(...) / write_volume_frame(...)
call timing_io_stop()
```

---

### `wait_io (s)`

**Compute ranks only.** Time blocked in `wait_all_send()` after posting the non-blocking
`MPI_Isend` to the IO server. Non-zero only when `HDF5_IO_NODES > 0`.

**Source — `save_forward_arrays_hdf5.F90` (lines ~661–666), called from within
`save_forward_arrays_undoatt()` when `HDF5_ENABLED = .true.`:**
```fortran
! split timing: IO timer paused while waiting for MPI send to complete
call timing_io_stop()
call timing_wait_start()
call wait_all_send()      ! MPI_Waitall on the isend requests
call timing_wait_stop()
call timing_io_start()    ! resume for caller's stop
```

> In practice this is often ~0 because the IO server is parked in a blocking
> `MPI_Probe` and accepts the send immediately, so `MPI_Waitall` returns
> in nanoseconds.

Always **zero for IO server ranks**.

---

### `idle_io (s)`

**IO server ranks only.** Time blocked in `MPI_Probe` waiting for the next message
from compute ranks.

**Source — `hdf5_io_server.F90` (lines ~1607–1609):**
```fortran
call timing_idle_start()
call idle_mpi_io(status)   ! wraps blocking MPI_Probe(MPI_ANY_SOURCE, ...)
call timing_idle_stop()
```

On a fast compute node with light IO load, the IO server spends the vast majority of
its time here (see example results below).

Always **zero for compute ranks**.

---

### `other (s)`

Residual time not captured by any of the four categories:

```
other = total - compute - io - wait_io - idle_io
```

Accounts for simulation initialization, stability checks, VTK updates, MPI collective
operations outside the timed regions, and the seismogram flush at the end of the loop.

---

### `total (s)`

Wall-clock elapsed time for this interval, measured from the last `timing_reset()` to
the `timing_report()` call. The sum of all intervals' `total` values equals the overall
simulation wall-clock time.

**Source — `timing_accounting.F90`:**
```fortran
total_elapsed = MPI_Wtime() - t_interval_start   ! t_interval_start set by timing_reset()
```

---

### `mpi_wtime (s)`, `date`, `time`

`MPI_Wtime()` at the moment `timing_report` is called (epoch-relative wall clock),
plus the system date/time from `date_and_time()`. Useful for correlating per-rank
log entries and correlating with system-level profiling traces.

---

## Example results

Test configuration: `regional_Greece_small_LDDRK` example on Fugaku (A64FX),
`NCHUNKS=1`, `NPROC_XI=2`, `NPROC_ETA=2`, `NSTEP=700`, `UNDO_ATTENUATION=.true.`,
`SIMULATION_TYPE=1`, `MOVIE_VOLUME=.true.`, `MOVIE_COARSE=.true.`,
`NTSTEP_BETWEEN_FRAMES=50`, `NTSTEP_BETWEEN_OUTPUT_INFO=200`.
`NT_DUMP_ATTENUATION_VAL=1802` (computed by mesher), so the entire run fits in a
single undo-attenuation subset. Periodic recording produces 4 records per rank
(at steps 200, 400, 600, 700).

### Case 1: HDF5_IO_NODES=0 (no IO server, binary checkpoint)

```
timing_acct_group-1_compute_rank0.txt:
it_begin: 1, it_end: 200, ..., compute (s):  144.44, io (s):  2.91, ..., total (s):  147.39
it_begin: 201, it_end: 400, ..., compute (s):   50.08, io (s):  2.50, ..., total (s):   52.60
it_begin: 401, it_end: 600, ..., compute (s):   45.61, io (s):  2.43, ..., total (s):   48.06
it_begin: 601, it_end: 700, ..., compute (s):   22.51, io (s):  1.37, ..., total (s):   23.90
```

- **4 records per rank**: at steps 200, 400, 600, and 700 (every `NTSTEP_BETWEEN_OUTPUT_INFO` + final step).
- **First interval (1–200) has 3× more compute** than later intervals: includes
  undo-attenuation checkpoint save (`save_forward_arrays_undoatt()`) overhead at subset start,
  plus cache warm-up effects.
- **IO proportional to movie frames**: 4 frames in 200-step intervals (~2.5s IO),
  2 frames in the 100-step interval (~1.4s IO). Each frame writes `MOVIE_VOLUME_TYPE=2`
  (div+curl) fields for all GLL points.
- `wait_io` and `idle_io` are 0 (no IO server).

### Case 2: HDF5_IO_NODES=1 (1 IO server, HDF5 checkpoint)

```
timing_acct_group-1_compute_rank0.txt:
it_begin: 1, it_end: 200, ..., compute (s):  144.33, io (s):  2.53, ..., total (s):  146.90
it_begin: 201, it_end: 400, ..., compute (s):   49.95, io (s):  2.10, ..., total (s):   52.08
it_begin: 401, it_end: 600, ..., compute (s):   45.49, io (s):  2.08, ..., total (s):   47.59
it_begin: 601, it_end: 700, ..., compute (s):   22.39, io (s):  1.14, ..., total (s):   23.55

timing_acct_group-1_io_server_rank0.txt:
it_begin: 0, it_end: 0, ..., io (s):  0.13, idle_io (s):  269.87, ..., total (s):  270.19
```

- Compute ranks: same pattern as case 1, but slightly less IO time (~2.1s vs ~2.5s
  for 200-step intervals) thanks to HDF5 batching.
- IO server: 0s compute, 0.13s IO, 269.9s idle — the server spent nearly all its time
  blocked in `MPI_Probe` waiting for messages.
- IO server `it_begin: 0, it_end: 0` — the final residual record (no per-snapshot
  records because `NT_DUMP_ATTENUATION_VAL > NSTEP`, so no undo checkpoint IO occurred
  through the IO server).

---

## Files changed

| File | Change |
|---|---|
| `src/specfem3D/timing_accounting.F90` | Module: accumulators, start/stop helpers, periodic `timing_report(rank, group, is_io_node, it_start, it_end)` with `t_interval_start` tracking |
| `src/specfem3D/rules.mk` | Added `$O/timing_accounting.solverstatic.o` to OBJECTS; dependency rules for 4 consumers |
| `src/specfem3D/iterate_time.F90` | Periodic recording every `NTSTEP_BETWEEN_OUTPUT_INFO` steps and at `it_end` |
| `src/specfem3D/iterate_time_undoatt.F90` | Periodic recording every `NTSTEP_BETWEEN_OUTPUT_INFO` steps within each undo-attenuation subset, plus conditional recording at subset boundaries |
| `src/specfem3D/save_forward_arrays_hdf5.F90` | `wait_io` bracket around `wait_all_send()` inside `save_forward_arrays_undoatt_hdf5()` |
| `src/specfem3D/hdf5_io_server.F90` | Per-undo-snapshot recording after each `write_buffered_undo_snapshot()` + final residual record |
