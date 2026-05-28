# Per-Process Timing Accounting

## Overview

`src/specfem3D/timing_accounting.F90` provides a barrier-free, per-process wall-clock
timing breakdown for the solver time loop. Each MPI rank independently accumulates time
into four categories and writes a one-line report at the end of the simulation.

The design constraint is **zero synchronization barriers**: all timing is done with
`MPI_Wtime()` calls local to each process, so the measurement itself has no impact on
parallel efficiency or on the phenomena being measured.

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

Each file contains a single line:

```
mygroup: G, myrank: R, role: ROLE, compute (s): C, io (s): I, wait_io (s): W,
idle_io (s): D, other (s): O, total (s): T, mpi_wtime (s): E, date: DATE, time: TIME
```

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

`wtime() - time_start` measured from just before the main time loop to just after it
exits. Matches the value printed by `print_elapsed_time()`.

**Source — `iterate_time_undoatt.F90` (line ~694):**
```fortran
call timing_report(myrank, mygroup, .false., wtime() - time_start)
call print_elapsed_time()
```

---

### `mpi_wtime (s)`, `date`, `time`

`MPI_Wtime()` at the moment `timing_report` is called (epoch-relative wall clock),
plus the system date/time from `date_and_time()`. Useful for correlating per-rank
log entries and correlating with system-level profiling traces.

---

## Example results

**System**: Fugaku (Fujitsu A64FX), 1 node, 1 chunk  
**Example**: `EXAMPLES/regional_Greece_small_LDDRK`  
**Settings**: `NCHUNKS=1`, `NPROC_XI=2`, `NPROC_ETA=2`, `NSTEP=700`,
`UNDO_ATTENUATION=.true.`, `SIMULATION_TYPE=1`

### Case 1 — `HDF5_IO_NODES = 0` (standard binary IO, 4 compute ranks)

```
timing_acct_group-1_compute_rank0.txt:
  mygroup: -1, myrank: 0, role: compute,
    compute (s):   262.113682,  io (s):   9.274764,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.114510,  total (s): 271.502956

timing_acct_group-1_compute_rank1.txt:
  mygroup: -1, myrank: 1, role: compute,
    compute (s):   261.976416,  io (s):   9.444544,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.049303,  total (s): 271.470263

timing_acct_group-1_compute_rank2.txt:
  mygroup: -1, myrank: 2, role: compute,
    compute (s):   261.887815,  io (s):   9.537152,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.049648,  total (s): 271.474614

timing_acct_group-1_compute_rank3.txt:
  mygroup: -1, myrank: 3, role: compute,
    compute (s):   261.597632,  io (s):   9.886255,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.051348,  total (s): 271.535234
```

| category | fraction of total |
|---|---|
| compute | ~96.5 % |
| io (binary undo-att checkpoint writes) | ~3.4 % |
| wait_io | 0 % (no IO server) |
| idle_io | 0 % (no IO server) |
| other | ~0.04 % |

---

### Case 2 — `HDF5_IO_NODES = 1`, `HDF5_ENABLED = .true.` (4 compute + 1 IO server)

```
timing_acct_group-1_compute_rank0.txt:
  mygroup: -1, myrank: 0, role: compute,
    compute (s):   261.640241,  io (s):   7.850409,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.096213,  total (s): 269.586863

timing_acct_group-1_compute_rank1.txt:
  mygroup: -1, myrank: 1, role: compute,
    compute (s):   261.511526,  io (s):   8.034591,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.049422,  total (s): 269.595539

timing_acct_group-1_compute_rank2.txt:
  mygroup: -1, myrank: 2, role: compute,
    compute (s):   261.361392,  io (s):   8.184355,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.049793,  total (s): 269.595539

timing_acct_group-1_compute_rank3.txt:
  mygroup: -1, myrank: 3, role: compute,
    compute (s):   261.056227,  io (s):   8.487966,
    wait_io (s):     0.000000,  idle_io (s):   0.000000,
    other (s):       0.051347,  total (s): 269.595540

timing_acct_group-1_io_server_rank0.txt:
  mygroup: -1, myrank: 0, role: io_server,
    compute (s):     0.000000,  io (s):   0.085370,
    wait_io (s):     0.000000,  idle_io (s): 269.340074,
    other (s):       0.191127,  total (s): 269.616572
```

| rank | role | compute | io | wait_io | idle_io |
|---|---|---|---|---|---|
| 0–3 | compute | ~97 % | ~3 % (data pack + MPI_Isend) | ~0 % | 0 % |
| 0 (IO server) | io_server | 0 % | ~0.03 % (HDF5 writes) | 0 % | ~99.9 % |

**Key observations from Case 2:**

1. **IO server spends 99.9 % of its time idle** in `idle_mpi_io()` / `MPI_Probe`
   (269.3 s idle vs. 0.085 s actual HDF5 writes). The write work is negligible compared
   to the wait time between checkpoints.

2. **`wait_io ≈ 0` for compute ranks** — the IO server is parked in a blocking
   `MPI_Probe`, so it accepts the `MPI_Isend` essentially instantaneously and
   `wait_all_send()` returns in nanoseconds.

3. **Compute `io` decreased vs. Case 1** (~8.1 s vs. ~9.5 s, −15 %) — the actual disk
   write latency is hidden on the IO server while compute ranks move on. The remaining
   compute io time is the data-packing cost.

4. **Total elapsed time decreased** (~269.6 s vs. ~271.5 s) despite adding an extra
   MPI rank, confirming the IO server removes a blocking disk-write from the critical
   path.

---

## Files changed

| File | Change |
|---|---|
| `src/specfem3D/timing_accounting.F90` | New module: accumulators, start/stop helpers, `timing_report()` |
| `src/specfem3D/rules.mk` | Added `$O/timing_accounting.solverstatic.o` to OBJECTS; dependency rules for 4 consumers |
| `src/specfem3D/iterate_time.F90` | `use timing_accounting`; `timing_reset()`; compute/io brackets; `timing_report()` call |
| `src/specfem3D/iterate_time_undoatt.F90` | Same as above for undo-attenuation time loop |
| `src/specfem3D/save_forward_arrays_hdf5.F90` | `wait_io` bracket around `wait_all_send()` inside `save_forward_arrays_undoatt_hdf5()` |
| `src/specfem3D/hdf5_io_server.F90` | `timing_reset()`; `idle_io` bracket around `idle_mpi_io()`; `io` brackets around all recv/write calls; `timing_report()` call |
