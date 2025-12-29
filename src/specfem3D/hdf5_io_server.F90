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

module io_server_hdf5

  use constants, only: CUSTOM_REAL
  use shared_parameters, only: HDF5_IO_NODES, &
    IO_storage_task, IO_compute_task, H5_COL

  implicit none

#ifdef USE_HDF5

  ! public routines
  public :: initialize_io_server
  public :: finalize_io_server
  public :: do_io_start_idle
  public :: pass_info_to_io
  public :: get_info_from_comp
  public :: wait_all_send

  ! functions to prepare the offset arrays for the HDF5 output
  public :: initialize_hdf5_solver
  public :: finalize_hdf5_solver
  public :: allocate_offset_arrays
  public :: deallocate_offset_arrays
  public :: create_hdf5_files_and_datasets

  ! public parameters
  public :: dest_ionod

  ! MPI message tags
  ! undo attenuation
  public :: io_tag_ford_undo_d_cm, io_tag_ford_undo_v_cm, io_tag_ford_undo_a_cm
  public :: io_tag_ford_undo_d_oc, io_tag_ford_undo_v_oc, io_tag_ford_undo_a_oc
  public :: io_tag_ford_undo_d_ic, io_tag_ford_undo_v_ic, io_tag_ford_undo_a_ic
  public :: io_tag_ford_undo_eps_xx_cm, io_tag_ford_undo_eps_yy_cm, io_tag_ford_undo_eps_xy_cm
  public :: io_tag_ford_undo_eps_xz_cm, io_tag_ford_undo_eps_yz_cm
  public :: io_tag_ford_undo_eps_xx_ic, io_tag_ford_undo_eps_yy_ic, io_tag_ford_undo_eps_xy_ic
  public :: io_tag_ford_undo_eps_xz_ic, io_tag_ford_undo_eps_yz_ic
  public :: io_tag_ford_undo_A_rot, io_tag_ford_undo_B_rot
  public :: io_tag_ford_undo_R_xx_cm, io_tag_ford_undo_R_yy_cm, io_tag_ford_undo_R_xy_cm
  public :: io_tag_ford_undo_R_xz_cm, io_tag_ford_undo_R_yz_cm
  public :: io_tag_ford_undo_R_xx_ic, io_tag_ford_undo_R_yy_ic, io_tag_ford_undo_R_xy_ic
  public :: io_tag_ford_undo_R_xz_ic, io_tag_ford_undo_R_yz_ic
  public :: io_tag_ford_undo_neq, io_tag_ford_undo_neq1, io_tag_ford_undo_pgrav1
  public :: io_tag_ford_undo_nmsg

  ! surface movie metadata
  public :: io_tag_surf_offset, io_tag_surf_npoints

  ! surface movie data (ux, uy, uz)
  public :: io_tag_surf_ux, io_tag_surf_uy, io_tag_surf_uz

  ! volume movie metadata (movie-point offsets and total count)
  public :: io_tag_vol_offset, io_tag_vol_npoints

  ! volume movie data (strain and vector components at movie points)
  public :: io_tag_vol_strain_NN, io_tag_vol_strain_EE, io_tag_vol_strain_ZZ
  public :: io_tag_vol_strain_NE, io_tag_vol_strain_NZ, io_tag_vol_strain_EZ
  public :: io_tag_vol_vec_N,     io_tag_vol_vec_E,     io_tag_vol_vec_Z

  ! MPI requests
  public :: n_req_ford_undo
  public :: req_dump_ford_undo
  public :: n_msg_ford_undo

  ! surface movie nonblocking send requests
  public :: n_req_surf
  public :: req_dump_surf

  ! volume movie nonblocking send requests
  public :: n_req_vol
  public :: req_dump_vol

  public :: nproc_io

  ! verbosity
  public :: VERBOSE



  private

  ! MPI tags for io server implementation
  ! for forward undo att arrays
  integer :: io_tag_ford_undo_d_cm      = 100001
  integer :: io_tag_ford_undo_v_cm      = 100002
  integer :: io_tag_ford_undo_a_cm      = 100003
  integer :: io_tag_ford_undo_d_oc      = 100004
  integer :: io_tag_ford_undo_v_oc      = 100005
  integer :: io_tag_ford_undo_a_oc      = 100006
  integer :: io_tag_ford_undo_d_ic      = 100007
  integer :: io_tag_ford_undo_v_ic      = 100008
  integer :: io_tag_ford_undo_a_ic      = 100009
  integer :: io_tag_ford_undo_eps_xx_cm = 100010
  integer :: io_tag_ford_undo_eps_yy_cm = 100011
  integer :: io_tag_ford_undo_eps_xy_cm = 100012
  integer :: io_tag_ford_undo_eps_xz_cm = 100013
  integer :: io_tag_ford_undo_eps_yz_cm = 100014
  integer :: io_tag_ford_undo_eps_xx_ic = 100015
  integer :: io_tag_ford_undo_eps_yy_ic = 100016
  integer :: io_tag_ford_undo_eps_xy_ic = 100017
  integer :: io_tag_ford_undo_eps_xz_ic = 100018
  integer :: io_tag_ford_undo_eps_yz_ic = 100019
  integer :: io_tag_ford_undo_A_rot     = 100020
  integer :: io_tag_ford_undo_B_rot     = 100021
  integer :: io_tag_ford_undo_R_xx_cm   = 100022
  integer :: io_tag_ford_undo_R_yy_cm   = 100023
  integer :: io_tag_ford_undo_R_xy_cm   = 100024
  integer :: io_tag_ford_undo_R_xz_cm   = 100025
  integer :: io_tag_ford_undo_R_yz_cm   = 100026
  integer :: io_tag_ford_undo_R_xx_ic   = 100027
  integer :: io_tag_ford_undo_R_yy_ic   = 100028
  integer :: io_tag_ford_undo_R_xy_ic   = 100029
  integer :: io_tag_ford_undo_R_xz_ic   = 100030
  integer :: io_tag_ford_undo_R_yz_ic   = 100031
  integer :: io_tag_ford_undo_neq       = 100032
  integer :: io_tag_ford_undo_neq1      = 100033
  integer :: io_tag_ford_undo_pgrav1    = 100034
  integer :: io_tag_nsubset_iterations  = 100035
  integer :: io_tag_ford_undo_nmsg      = 100036

  ! surface movie metadata tags
  integer :: io_tag_surf_offset         = 110001
  integer :: io_tag_surf_npoints        = 110002

  ! surface movie data tags
  integer :: io_tag_surf_ux             = 110010
  integer :: io_tag_surf_uy             = 110011
  integer :: io_tag_surf_uz             = 110012

  ! surface movie time-stepping tags
  integer :: io_tag_surf_it_begin       = 110020
  integer :: io_tag_surf_it_end         = 110021

  ! volume movie metadata tags
  integer :: io_tag_vol_offset          = 120001
  integer :: io_tag_vol_npoints         = 120002

  ! volume movie data tags (movie-point-based)
  integer :: io_tag_vol_strain_NN       = 120010
  integer :: io_tag_vol_strain_EE       = 120011
  integer :: io_tag_vol_strain_ZZ       = 120012
  integer :: io_tag_vol_strain_NE       = 120013
  integer :: io_tag_vol_strain_NZ       = 120014
  integer :: io_tag_vol_strain_EZ       = 120015
  integer :: io_tag_vol_vec_N           = 120020
  integer :: io_tag_vol_vec_E           = 120021
  integer :: io_tag_vol_vec_Z           = 120022

  ! mpi_req dump (used in wait_all_send)
  integer :: n_req_ford_undo = 0
  integer, dimension(34) :: req_dump_ford_undo
  integer :: n_msg_ford_undo = 0

  ! surface movie nonblocking send requests
  integer :: n_req_surf = 0
  integer, dimension(3) :: req_dump_surf

  ! volume movie nonblocking send requests
  integer :: n_req_vol = 0
  integer, dimension(6) :: req_dump_vol

  ! responsible id of io node
  integer :: dest_ionod = 0
  ! number of computer nodes sending info to each io node
  integer :: nproc_io
  ! id for io node
  integer :: my_io_id

  ! mapping from IO nodes to the compute ranks they serve
  ! size: (number of IO nodes) x (max number of compute ranks per IO node)
  integer, allocatable :: io_compute_ranks(:,:)
  integer, allocatable :: io_nproc_all(:)

  ! ordered list of undo-attenuation MPI tags per snapshot
  integer, allocatable :: undo_tag_list(:)
  integer :: n_undo_tags

  ! verbose output (for debugging)
  logical, parameter :: VERBOSE = .true.


! USE_HDF5
#endif

contains

  subroutine initialize_io_server()

  ! initialization of IO server splits compute and io nodes
#ifdef USE_HDF5

  use constants, only: IMAIN,myrank,MAX_STRING_LEN
  use shared_parameters, only: NPROCTOT !,NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer :: sizeval,irank
  integer :: mpi_comm,split_comm,inter_comm
  integer :: key,io_start,comp_start
  ! test node name
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: node_names
  character(len=MAX_STRING_LEN) :: this_node_name
  integer :: node_len

  ! split comm into computation nodes and io node
  ! here we use the last HDF5_IO_NODES ranks (intra comm) as the io node
  ! for using one additional node, xspecfem3D need to be run with + HDF5_IO_NODES node
  ! thus for running mpirun, it should be like e.g.
  ! mpirun -n $((NPROC+HDF5_IO_NODES)) ./bin/xspecfem3D

  ! safety check
  if (HDF5_IO_NODES < 0) stop 'Invalid HDF5_IO_NODES, must be zero or positive'

  ! initializes
  call world_get_comm(mpi_comm)         ! mpi_comm == my_local_mpi_comm_world

  ! default MPI group, no i/o server
  call world_set_comm_inter(mpi_comm)   ! my_local_mpi_comm_inter = mpi_comm

  ! checks if anything to do
  if (HDF5_IO_NODES == 0) return

  ! select the task of this proc
  ! get the local mpi_size and rank
  call world_size(sizeval)
  call world_rank(myrank)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'HDF5 I/O server run:'
    write(IMAIN,*) '  number of dedicated HDF5 IO nodes = ',HDF5_IO_NODES
    write(IMAIN,*) '  total number of MPI processes     = ',sizeval
    write(IMAIN,*)
    write(IMAIN,*) '  separating subgroups for compute tasks with ',sizeval - HDF5_IO_NODES,'MPI processes'
    write(IMAIN,*) '                       for io tasks      with ',HDF5_IO_NODES,'MPI processes'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks if solver run with a number of MPI processes equal to (NPROC + HDF5_IO_NODES)
  if (sizeval /= NPROCTOT + HDF5_IO_NODES) then
    if (myrank == 0) then
      print *,"Error: HDF5 IO server needs ",HDF5_IO_NODES," additional process together with NPROCTOT",NPROCTOT,"processes"
      print *
      print *,"However, this simulation runs with ",sizeval,"MPI processes,"
      print *,"where instead it should run with ",NPROCTOT + HDF5_IO_NODES,"MPI processes"
      print *
      print *,"Please check your call to 'mpirun -np .. xspecfem3D' and rerun with the correct number of processes"
      print *
    endif
    stop 'Invalid number of MPI processes for HDF5_IO_NODES > 0'
  endif

  ! to select io node and compute nodes on the same cluster node.
  allocate(node_names(0:sizeval-1))

  call world_get_processor_name(this_node_name,node_len)

  ! share the node name between procs
  call gather_all_all_single_ch(this_node_name,node_names,sizeval,MAX_STRING_LEN)

  call synchronize_all()

  ! select the task of this proc
  call select_io_node(node_names,myrank,sizeval,key,io_start)

  deallocate(node_names)

  ! split communicator into compute_comm and io_comm
  call world_comm_split(mpi_comm, key, myrank, split_comm)

  ! create inter communicator and set as my_local_mpi_comm_inter
  if (IO_storage_task) then
    comp_start = 0
    call world_create_intercomm(split_comm, 0, mpi_comm, comp_start, 11111, inter_comm)
  else
    !io_start = sizeval - HDF5_IO_NODES !+dest_ionod
    call world_create_intercomm(split_comm, 0, mpi_comm, io_start, 11111, inter_comm)
  endif

  ! sets new comm_world subgroup
  call world_set_comm(split_comm)            ! such that: my_local_mpi_comm_world == split_comm

  ! use inter_comm as my_local_mpi_comm_world for all send/recv
  call world_set_comm_inter(inter_comm)      ! such that: my_local_mpi_comm_inter == inter_comm

  ! exclude io nodes from the other compute nodes
  !if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) NPROCTOT = NPROCTOT - HDF5_IO_NODES

  ! re-define myrank (within new comm_world subgroup)
  call world_rank(myrank)
  call world_size(sizeval)

  ! debug print inter_comm and if IO_compute_task
  if (VERBOSE) then
    do irank = 0, sizeval-1
      if (myrank == irank) then
        if (IO_compute_task) then
          print *, "io_server: rank ", myrank, " compute task, my_local_mpi_comm_inter = ", &
            inter_comm, " my_local_mpi_comm_world = ", split_comm
        else
          print *, "io_server: rank ", myrank, " io task, my_local_mpi_comm_inter = ", &
            inter_comm, " my_local_mpi_comm_world = ", split_comm
        endif
      endif
      call flush_stdout()
      call synchronize_all()
    enddo
  endif

  call synchronize_inter()

  ! note: IO compute tasks cover the lower part of MPI processes, i.e., with ranks in
  !       the initial 0 to (sizeval-HDF5_IO_NODES)-1 range.
  !       the io nodes with IO_storage_task set to .true. are added at the end, with ranks in
  !       the upper range (sizeval-HDF5_IO_NODES) to (sizeval-1).
  !
  !       for user output, we only have the initial rank 0 as the main process opening the output_solver.txt file.
  !       after separating the tasks and creating new subgroup, we have a rank 0 process in the compute subgroup
  !       as well as a rank 0 process in the storage subgroup.
  !
  !       from here on, we only want the process with rank 0 in the compute group write output to IMAIN.
  !
  ! user output
  if (IO_compute_task) then
    if (myrank == 0) then
      write(IMAIN,*) '  new number of processes in compute group = ',sizeval
      write(IMAIN,*) '  new addtional processes in io group      = ',HDF5_IO_NODES
      write(IMAIN,*) '  MPI i/o server setup done'
      call flush_IMAIN()
    endif
  endif

  ! allocate the offset arrays
  !call initialize_hdf5_solver()

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine initialize_io_server() called without HDF5 Support."
  print *, "       HDF5_IO_NODES > 0 requires HDF5 support and HDF5_ENABLED must be set to .true."
  print *
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'initialize_io_server() called without HDF5 compilation support'

#endif

  end subroutine initialize_io_server

  !
  !-------------------------------------------------------------------------------------------------
  !

  subroutine finalize_io_server()

#ifdef USE_HDF5

  implicit none

  ! checks if anything to do
  if (HDF5_IO_NODES == 0) return

  ! finish MPI subgroup
  if (HDF5_IO_NODES > 0) then
    ! deallocate offset arrays
    !call deallocate_offset_arrays()
    ! wait for all to finish
    call synchronize_inter()
    ! free subgroup
    call world_comm_free_inter()
  endif
#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine finalize_io_server() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'finalize_io_server() called without HDF5 compilation support'

#endif

  end subroutine finalize_io_server

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine select_io_node(node_names, myrank, sizeval, key, io_start)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer, intent(in)                         :: myrank, sizeval
  integer, intent(out)                        :: key, io_start
  character(len=MAX_STRING_LEN), dimension(0:sizeval-1), intent(in) :: node_names
  ! local parameters
  integer, dimension(sizeval) :: n_procs_on_node ! number of procs on each cluster node
  integer, dimension(:), allocatable :: n_ionode_on_cluster ! number of ionode on the cluster nodes
  integer :: i,j,c,n_cluster_node,my_cluster_id,n_rest_io,n_ionode,n_comp_node
  integer :: comp_rank_counter, dest_io_id, idx
  real(kind=CUSTOM_REAL) :: io_ratio ! dum
  character(len=MAX_STRING_LEN), dimension(sizeval) :: dump_node_names ! names of cluster nodes

  ! initialize
  my_cluster_id = -1
  n_cluster_node = 0
  n_procs_on_node(:) = 0
  dump_node_names(:) = "nan"

  ! BUG: nprocs goes wrong when 2 cluster nodes and 64 compute 8 io procs
  ! (only 62 compute nodes are assigned)

  ! cluster_node_nums = [n_1,...,n_i,...n_cn,-1,-1,...] ! n_i is the number of procs on
  ! each cluster node
  do i = 1, sizeval
    ! search the node name already found and registered in the dump_node_names array
    c = 0
    do j = 1, sizeval
      if (node_names(i-1) == dump_node_names(j)) c = j
    enddo

    ! if node_name[i-1] is not registered yet
    if (c == 0) then
      ! count the number of cluster node
      n_cluster_node = n_cluster_node + 1

      dump_node_names(n_cluster_node) = node_names(i-1) ! register the name
      n_procs_on_node(n_cluster_node) = 1 ! count up the number of procs on this cluster node

      c = n_cluster_node
    else ! if node_name[i] has already be found, count up the number of procs
      n_procs_on_node(c) = n_procs_on_node(c) + 1
    endif

    ! the id of cluster which this process belongs to
    if (i-1 == myrank) then
      my_cluster_id = c
    endif
  enddo

  ! warning when a cluster node has only one single proc.
  do i = 1, n_cluster_node
    if (n_procs_on_node(i) == 1) then
      print *
      print *, "***************************************************************************"
      print *, "IO server Warning:"
      print *, "  node name: " // trim(dump_node_names(i)) // " has only one procs."
      print *, "  this may lead a io performance issue by inter-clusternode communication."
      print *, "***************************************************************************"
      print *
    endif
  enddo

  ! select HDF5_IO_NODES of io nodes
  allocate(n_ionode_on_cluster(n_cluster_node))
  !! decide the number of io node on each cluster node
  ! at least one io node on each cluster node
  n_ionode_on_cluster(:) = 1

  ! check if the total number of io node > HDF5_IO_NODES
  if (sum(n_ionode_on_cluster) > HDF5_IO_NODES) then
    print *, "Error: HDF5_IO_NODES in Par_file is too small,"
    print *, "       at least one io node for each cluster node is necessary"
    stop 'Invalid HDF5_IO_NODES value too small'
  endif

  ! share the rest of ionodes based on the ratio of compute nodes
  n_rest_io = HDF5_IO_NODES - sum(n_ionode_on_cluster)
  if (n_rest_io /= 0) then
    ! add io nodes one by one to the cluster node where
    ! the io_node/total_node ratio is rowest
    do i = 1, n_rest_io
      io_ratio = 9999.0 !initialize at each i
      do j = 1, n_cluster_node
        if (io_ratio > real(n_ionode_on_cluster(j))/real(n_procs_on_node(j))) then
          ! dump if largest
          io_ratio = real(n_ionode_on_cluster(j))/real(n_procs_on_node(j))
          ! destination of additional io node
          c = j
        endif
      enddo
      ! add one additional io node
      n_ionode_on_cluster(c) = n_ionode_on_cluster(c) + 1
    enddo
  endif

  !! choose the io node from the last rank of each cluster node
  n_ionode = 0
  comp_rank_counter = -1

  ! allocate and initialize IO-to-compute mapping arrays
  if (HDF5_IO_NODES > 0) then
    if (.not. allocated(io_nproc_all)) then
      allocate(io_nproc_all(HDF5_IO_NODES))
    endif
    if (.not. allocated(io_compute_ranks)) then
      allocate(io_compute_ranks(HDF5_IO_NODES,sizeval))
    endif
    io_nproc_all(:) = 0
  endif
  do i = 1, n_cluster_node
    c = 0
    ! number of compute node on this cluster node
    n_comp_node = n_procs_on_node(i) - n_ionode_on_cluster(i)
    do j = 1, sizeval
      ! if this rank is on i cluster node
      if (dump_node_names(i) == node_names(j-1)) then
        c = c + 1
        if (n_comp_node < c) then
          ! j is io node
          if (j-1 == myrank) then
            ! check  if j is myrank
            IO_storage_task = .true. ! set io node flag
            IO_compute_task = .false.
            key          = 0
            my_io_id     = n_ionode ! id of io_node
            nproc_io     = n_comp_node/n_ionode_on_cluster(i) ! number of compute nodes which use this io node
            dest_ionod   = -1
            if (c-n_comp_node-1 < mod(n_comp_node,n_ionode_on_cluster(i))) nproc_io = nproc_io + 1
          endif

          ! rank of io_start
          if (n_ionode == 0) io_start = j-1

          n_ionode = n_ionode+1

        else
          ! j is compute node
          comp_rank_counter = comp_rank_counter + 1
          dest_io_id = mod(c-1,n_ionode_on_cluster(i)) + n_ionode

          if (HDF5_IO_NODES > 0) then
            idx = io_nproc_all(dest_io_id+1) + 1
            io_nproc_all(dest_io_id+1) = idx
            io_compute_ranks(dest_io_id+1,idx) = comp_rank_counter
          endif

          if (j-1 == myrank) then
            IO_storage_task = .false.
            IO_compute_task = .true.
            key          = 1
            ! set the destination of MPI communication
            dest_ionod   = dest_io_id
          endif
        endif
      endif
    enddo
  enddo

  ! debug
  if (VERBOSE) then
    if (myrank == 0) then
      print *, "io_server: n_procs_on_node", n_procs_on_node(:)
      print *, "io_server: n_ionode_on_cluster", n_ionode_on_cluster
      print *
    endif
    call flush_stdout()
    call synchronize_all()

    if (IO_storage_task) then
      print *, "io_server: io_task rank ", myrank, " nprocio ",nproc_io
      print *, "io_server: io_task rank ", myrank, " node ", trim(node_names(myrank)), " my_io_id", my_io_id
      print *
    endif
    call flush_stdout
    call synchronize_all()

    do i = 0, n_ionode-1
      if (.not. IO_storage_task) then
        if (dest_ionod == i) then
          print *, "io_server: rank ", myrank, " node ", trim(node_names(myrank)), " dest ", dest_ionod
        endif
      endif
      call flush_stdout()
      call synchronize_all()
    enddo
    call flush_stdout()
    call synchronize_all()

  endif

  deallocate(n_ionode_on_cluster)

  end subroutine select_io_node

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine do_io_start_idle()

#ifdef USE_HDF5

  use specfem_par
  use specfem_par_movie_hdf5
  use manager_hdf5
  use constants, only: myrank, my_status_size, my_status_source, my_status_tag

  use io_bandwidth

  implicit none

  integer :: status(my_status_size)
  integer :: tag, tag_src

  ! vars forward undo att arrays
  ! undo attenuation
  integer :: n_recv_msg_ford_undo, &      ! number of messages received for undo attenuation of one iteration
             max_ford_undo_out, &         ! number of iterations when IO happens for undo attenuation
             ford_undo_out_count          ! count the completed iterations for undo attenuation

  ! surface movie
  integer :: n_recv_msg_surf, &           ! number of messages received for one surface movie frame
             max_surf_frames, &           ! maximum number of surface movie frames
             surf_frame_count, &          ! count the completed surface movie frames
             n_msg_surf, &                ! number of messages per surface movie frame
             max_surf_points, &           ! maximum number of local surface points on any compute rank
             it_first_surf                ! first time step index producing a surface frame

  ! volume movie (movie-point-based types: strain/vector)
  integer :: n_recv_msg_vol, &            ! number of messages received for one volume movie frame
             max_vol_frames, &            ! maximum number of volume movie frames
             vol_frame_count, &           ! count the completed volume movie frames
             n_msg_vol, &                 ! number of messages per volume movie frame
             max_vol_points, &            ! maximum number of local volume movie points on any compute rank
             it_first_vol                 ! first time step index producing a volume frame

  ! track which frame index already has its HDF5 group/datasets created
  integer :: surf_group_frame_prepared
  integer :: vol_group_frame_prepared

  ! array for dumping the data array
  real(kind=CUSTOM_REAL), dimension(:),         allocatable :: dump_ford_undo_1d_glob
  real(kind=CUSTOM_REAL), dimension(:,:),       allocatable :: dump_ford_undo_2d_glob
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),   allocatable :: dump_ford_undo_4d
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: dump_ford_undo_5d

  ! arrays for surface movie values (per rank)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_surf_ux
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_surf_uy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_surf_uz

  ! arrays for volume movie values at movie points (per rank)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol1
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol3
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol4
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol5
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol6

  ! maximum nglob and nspec in offset arrays
  integer :: max_nglob
  integer :: max_nspec

  integer :: i_out

  ! initialize all counters
  ! undo attenuation
  n_recv_msg_ford_undo = 0 ! number of messages received for undo attenuation of one iteration
  max_ford_undo_out    = 0 ! number of iterations when IO happens for undo attenuation
  ford_undo_out_count  = 0 ! count the completed iterations for undo attenuation

  ! surface movie
  n_recv_msg_surf   = 0
  max_surf_frames   = 0
  surf_frame_count  = 0
  n_msg_surf        = 0
  max_surf_points   = 0
  it_first_surf     = 0
  surf_group_frame_prepared = -1

  ! volume movie
  n_recv_msg_vol    = 0
  max_vol_frames    = 0
  vol_frame_count   = 0
  n_msg_vol         = 0
  max_vol_points    = 0
  it_first_vol      = 0
  vol_group_frame_prepared = -1


  !
  ! initialize
  !

  ! undo attenuation
  if (UNDO_ATTENUATION .and. SAVE_FORWARD) then

    ! initialize HDF5 MPI context on IO storage tasks for undo snapshots
    if (HDF5_ENABLED .and. IO_storage_task) then
      call world_get_comm(comm)
      call world_get_info_null(info)
      call h5_initialize()
      call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)
    endif

    ! create all the HDF5 files and datasets for undo attenuation (one per snapshot)
    do i_out = 1, NSUBSET_ITERATIONS
      call create_hdf5_files_and_datasets(i_out)
    enddo

    n_msg_ford_undo   = n_msg_ford_undo*nproc_io ! multiply by the number of compute nodes for each io node
    max_ford_undo_out = NSUBSET_ITERATIONS ! multiply by the number of snapshots

    ! find the maximum number of nglob in offset_nglob_cm, _oc, _ic
    max_nglob = maxval([maxval(offset_nglob_cm), maxval(offset_nglob_oc), maxval(offset_nglob_ic)])

    ! find the maximum number of nspec in offset_nspec_cm_soa, _ic_soa, _oc_rot, _cm_att, _ic_att
    max_nspec = maxval([maxval(offset_nspec_cm_soa), maxval(offset_nspec_ic_soa), &
                        merge(maxval(offset_nspec_oc_rot), -huge(1), ROTATION_VAL), &
                        merge(maxval(offset_nspec_cm_att), -huge(1), ATTENUATION_VAL), &
                        merge(maxval(offset_nspec_ic_att), -huge(1), ATTENUATION_VAL)])

    ! allocate the dump arrays
    allocate(dump_ford_undo_1d_glob(max_nglob), &
             dump_ford_undo_2d_glob(NDIM, max_nglob), &
             dump_ford_undo_4d(NGLLX, NGLLY, NGLLZ, max_nspec), &
             dump_ford_undo_5d(NGLLX, NGLLY, NGLLZ, N_SLS, max_nspec))

  endif ! UNDO_ATTENUATION .and. SAVE_FORWARD

  ! surface movie initialization
  if (MOVIE_SURFACE) then

    ! number of messages per surface frame (ux, uy, uz from each compute rank handled by this IO node)
    n_msg_surf = 3 * nproc_io

    ! determine number of surface movie frames from time-stepping parameters
    ! frames are written whenever mod(it,NTSTEP_BETWEEN_FRAMES) == 0 within [it_begin,it_end]
    if (NTSTEP_BETWEEN_FRAMES > 0) then
      it_first_surf = ((it_begin + NTSTEP_BETWEEN_FRAMES - 1)/NTSTEP_BETWEEN_FRAMES) * NTSTEP_BETWEEN_FRAMES
      if (it_first_surf <= it_end) then
        max_surf_frames = (it_end - it_first_surf) / NTSTEP_BETWEEN_FRAMES + 1
        print *, 'io_server: surface movie frames from it=', it_first_surf, ' to it=', it_end, &
                 ' every ', NTSTEP_BETWEEN_FRAMES, ' steps: total ', max_surf_frames, ' frames'
      else
        max_surf_frames = 0
      endif
    else
      max_surf_frames = 0
    endif

    ! allocate buffers for the largest local surface chunk
    if (allocated(offset_poin)) then
      max_surf_points = maxval(offset_poin)
      if (max_surf_points > 0) then
        allocate(dump_surf_ux(max_surf_points))
        allocate(dump_surf_uy(max_surf_points))
        allocate(dump_surf_uz(max_surf_points))
      endif
    endif

  endif

  ! volume movie initialization (only for movie-point-based types: strains/vector)
  if (MOVIE_VOLUME) then

    ! MOVIE_VOLUME_TYPE 1-3: strains / time-integrated / potency (6 components)
    ! MOVIE_VOLUME_TYPE 5-6: displacement / velocity vectors (3 components)
    select case (MOVIE_VOLUME_TYPE)
    case (1,2,3)
      n_msg_vol = 6 * nproc_io
    case (5,6)
      n_msg_vol = 3 * nproc_io
    case default
      n_msg_vol = 0
    end select

    ! determine number of volume movie frames from movie parameters
    ! frames are written whenever mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0
    ! and it is within [MOVIE_START,MOVIE_STOP]
    if (NTSTEP_BETWEEN_FRAMES > 0) then
      ! first time step (within this run) that satisfies the movie sampling rule
      it_first_vol = MOVIE_START + ((max(it_begin, MOVIE_START) - MOVIE_START + NTSTEP_BETWEEN_FRAMES - 1) / NTSTEP_BETWEEN_FRAMES) * NTSTEP_BETWEEN_FRAMES

      ! limit movie duration to the actual simulated range
      if (min(it_end, MOVIE_STOP) >= it_first_vol) then
        max_vol_frames = (min(it_end, MOVIE_STOP) - it_first_vol) / NTSTEP_BETWEEN_FRAMES + 1
      else
        max_vol_frames = 0
      endif
    else
      max_vol_frames = 0
    endif

    ! allocate buffers for the largest local volume movie chunk
    if (allocated(offset_poin_vol)) then
      max_vol_points = maxval(offset_poin_vol)
      if (max_vol_points > 0) then
        allocate(dump_vol1(max_vol_points))
        allocate(dump_vol2(max_vol_points))
        allocate(dump_vol3(max_vol_points))
        allocate(dump_vol4(max_vol_points))
        allocate(dump_vol5(max_vol_points))
        allocate(dump_vol6(max_vol_points))
      endif
    endif

  endif


  ! initialize timer
  call initialize_bytes_written()

  if (VERBOSE .and. IO_storage_task) then
    print *, 'io_server: entering do_io_start_idle'
    print *, '  undo: max_ford_undo_out =', max_ford_undo_out, ' n_msg_ford_undo =', n_msg_ford_undo
    print *, '  surf: it_first_surf =', it_first_surf, ' max_surf_frames =', max_surf_frames, ' n_msg_surf =', n_msg_surf
    print *, '  vol : it_first_vol  =', it_first_vol,  ' max_vol_frames  =', max_vol_frames,  ' n_msg_vol  =', n_msg_vol
  endif

  !
  ! idling loop
  !
  do while ( (UNDO_ATTENUATION .and. SAVE_FORWARD .and. ford_undo_out_count < max_ford_undo_out) .or. &
             (MOVIE_SURFACE .and. HDF5_ENABLED .and. surf_frame_count < max_surf_frames) .or. &
             (MOVIE_VOLUME .and. HDF5_ENABLED .and. vol_frame_count  < max_vol_frames) )

    ! for surface movies, create the HDF5 group and datasets once per frame
    if (MOVIE_SURFACE) then
      if (max_surf_frames > 0 .and. n_msg_surf > 0) then
        if (surf_frame_count < max_surf_frames .and. n_recv_msg_surf == 0 .and. &
            surf_group_frame_prepared /= surf_frame_count) then
          if (myrank == 0) then
            call create_surface_frame_group(surf_frame_count, it_first_surf)
          endif
          call synchronize_all()
          if (VERBOSE .and. IO_storage_task) then
            print *, 'io_server: prepared surface frame group ', surf_frame_count, '/', max_surf_frames
          endif
          surf_group_frame_prepared = surf_frame_count
        endif
      endif
    endif

    ! for volume movies, create the HDF5 group and datasets once per frame
    if (MOVIE_VOLUME) then
      if (max_vol_frames > 0 .and. n_msg_vol > 0) then
        if (vol_frame_count < max_vol_frames .and. n_recv_msg_vol == 0 .and. &
            vol_group_frame_prepared /= vol_frame_count) then
          if (myrank == 0) then
            call create_volume_frame_group(vol_frame_count, it_first_vol)
          endif
          call synchronize_all()
          if (VERBOSE .and. IO_storage_task) then
            print *, 'io_server: prepared volume frame group ', vol_frame_count, '/', max_vol_frames
          endif
          vol_group_frame_prepared = vol_frame_count
        endif
      endif
    endif

    ! check the iteration counter


    ! waiting for a MPI message
    call idle_mpi_io(status)

    tag = status(my_status_tag)
    tag_src = status(my_status_source)

    ! debug output on the received message
    if (VERBOSE .and. IO_storage_task) then
      print *, 'io_server: rank', myrank, 'received message tag', tag, 'from rank', tag_src
      print *, '  ford_undo:', n_recv_msg_ford_undo, '/', n_msg_ford_undo
      print *, '  surf     :', n_recv_msg_surf, '/', n_msg_surf, ' frame ', surf_frame_count, '/', max_surf_frames
      print *, '  vol      :', n_recv_msg_vol,  '/', n_msg_vol,  ' frame ', vol_frame_count,  '/', max_vol_frames
    endif

    ! undo attenuation
    if (UNDO_ATTENUATION .and. SAVE_FORWARD .and. tag >= io_tag_ford_undo_d_cm .and. tag <= io_tag_ford_undo_pgrav1) then
      ! receive the data
      call recv_and_write_ford_undo(tag, tag_src, status, &
                                    dump_ford_undo_1d_glob, &
                                    dump_ford_undo_2d_glob, &
                                    dump_ford_undo_4d, &
                                    dump_ford_undo_5d, &
                                    ford_undo_out_count) ! use for the filename

      ! count 1 message received
      n_recv_msg_ford_undo = n_recv_msg_ford_undo + 1

    endif ! UNDO_ATTENUATION .and. SAVE_FORWARD

    ! surface movie
    if (MOVIE_SURFACE .and. &
        (tag == io_tag_surf_ux .or. tag == io_tag_surf_uy .or. tag == io_tag_surf_uz)) then

      call recv_and_write_surface_movie(tag, tag_src, status, &
                                        dump_surf_ux, dump_surf_uy, dump_surf_uz, &
                                        surf_frame_count, it_first_surf)

      ! count 1 message received
      n_recv_msg_surf = n_recv_msg_surf + 1

    endif

    ! volume movie (movie-point-based)
    if (MOVIE_VOLUME .and. &
        (tag == io_tag_vol_strain_NN .or. tag == io_tag_vol_strain_EE .or. tag == io_tag_vol_strain_ZZ .or. &
         tag == io_tag_vol_strain_NE .or. tag == io_tag_vol_strain_NZ .or. tag == io_tag_vol_strain_EZ .or. &
         tag == io_tag_vol_vec_N     .or. tag == io_tag_vol_vec_E     .or. tag == io_tag_vol_vec_Z)) then

      call recv_and_write_volume_movie(tag, tag_src, status, &
                                       dump_vol1, dump_vol2, dump_vol3, &
                                       dump_vol4, dump_vol5, dump_vol6, &
                                       vol_frame_count, it_first_vol)

      ! count 1 message received
      n_recv_msg_vol = n_recv_msg_vol + 1

    endif

    ! check receive counters
    if (UNDO_ATTENUATION .and. SAVE_FORWARD) then
      if (n_recv_msg_ford_undo >= n_msg_ford_undo) then
        ! increment the iteration counter
        ford_undo_out_count = ford_undo_out_count + 1
        ! reset the receive counter
        n_recv_msg_ford_undo = 0

        if (VERBOSE .and. IO_storage_task) then
          print *, 'io_server: completed undo iteration ', ford_undo_out_count, '/', max_ford_undo_out
        endif

        ! write out the times
        call calculate_bandwidth_all_procs()

        ! re-initialize timer
        call initialize_bytes_written()

      endif
    endif

    if (MOVIE_SURFACE) then
      if (max_surf_frames > 0 .and. n_msg_surf > 0) then
        if (n_recv_msg_surf >= n_msg_surf) then
          surf_frame_count = surf_frame_count + 1
          n_recv_msg_surf = 0

          if (VERBOSE .and. IO_storage_task) then
            print *, 'io_server: completed surface frame ', surf_frame_count, '/', max_surf_frames
          endif
        endif
      endif
    endif

    if (MOVIE_VOLUME) then
      if (max_vol_frames > 0 .and. n_msg_vol > 0) then
        if (n_recv_msg_vol >= n_msg_vol) then
          vol_frame_count = vol_frame_count + 1
          n_recv_msg_vol = 0

          if (VERBOSE .and. IO_storage_task) then
            print *, 'io_server: completed volume frame ', vol_frame_count, '/', max_vol_frames
          endif
        endif
      endif
    endif

  enddo
  !
  ! end of idling loop
  !

  call calculate_bandwidth_all_procs()

  if (VERBOSE .and. IO_storage_task) then
    print *, 'io_server: leaving do_io_start_idle'
    print *, '  final surf frames: ', surf_frame_count, '/', max_surf_frames
    print *, '  final vol  frames: ', vol_frame_count,  '/', max_vol_frames
  endif

  ! deallocate temporary arrays
  ! undo attenuation
  if (UNDO_ATTENUATION .and. SAVE_FORWARD) then
    deallocate(dump_ford_undo_1d_glob, &
               dump_ford_undo_2d_glob, &
               dump_ford_undo_4d, &
               dump_ford_undo_5d)
  endif

  ! surface movie buffers
  if (MOVIE_SURFACE) then
    if (allocated(dump_surf_ux)) deallocate(dump_surf_ux)
    if (allocated(dump_surf_uy)) deallocate(dump_surf_uy)
    if (allocated(dump_surf_uz)) deallocate(dump_surf_uz)
  endif

  ! volume movie buffers
  if (MOVIE_VOLUME) then
    if (allocated(dump_vol1)) deallocate(dump_vol1)
    if (allocated(dump_vol2)) deallocate(dump_vol2)
    if (allocated(dump_vol3)) deallocate(dump_vol3)
    if (allocated(dump_vol4)) deallocate(dump_vol4)
    if (allocated(dump_vol5)) deallocate(dump_vol5)
    if (allocated(dump_vol6)) deallocate(dump_vol6)
  endif


#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine do_io_start_idle() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'pass_info_to_io() called without HDF5 compilation support'

#endif

  end subroutine do_io_start_idle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_info_from_comp()

#ifdef USE_HDF5

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer, dimension(1) :: tmp_arr

  ! pass necessary information to io node

  ! forward undo att
  if (UNDO_ATTENUATION .and. SAVE_FORWARD) then
    ! receive NSUBSET_ITERATIONS
    call recv_i_inter(tmp_arr, 1, 0, io_tag_nsubset_iterations)
    NSUBSET_ITERATIONS = tmp_arr(1)

    ! receive offset arrays
    call recv_i_inter(offset_nglob_cm, NPROCTOT_VAL, 0, io_tag_ford_undo_d_cm) !!!!!!!!!!!!!!
    call recv_i_inter(offset_nglob_oc, NPROCTOT_VAL, 0, io_tag_ford_undo_d_oc)
    call recv_i_inter(offset_nglob_ic, NPROCTOT_VAL, 0, io_tag_ford_undo_d_ic)
    call recv_i_inter(offset_nspec_cm_soa, NPROCTOT_VAL, 0, io_tag_ford_undo_eps_xx_cm)
    call recv_i_inter(offset_nspec_ic_soa, NPROCTOT_VAL, 0, io_tag_ford_undo_eps_xx_ic)
    if (ROTATION_VAL) then
      call recv_i_inter(offset_nspec_oc_rot, NPROCTOT_VAL, 0, io_tag_ford_undo_A_rot)
    endif
    if (ATTENUATION_VAL) then
      call recv_i_inter(offset_nspec_cm_att, NPROCTOT_VAL, 0, io_tag_ford_undo_R_xx_cm)
      call recv_i_inter(offset_nspec_ic_att, NPROCTOT_VAL, 0, io_tag_ford_undo_R_xx_ic)
    endif
    if (FULL_GRAVITY_VAL) then
      call recv_i_inter(offset_pgrav1, NPROCTOT_VAL, 0, io_tag_ford_undo_pgrav1)
    endif

    ! receive the number of messages
    call recv_i_inter(tmp_arr, 1, 0, io_tag_ford_undo_nmsg)
    n_msg_ford_undo = tmp_arr(1)

    ! debug: log undo configuration on IO ranks
    if (VERBOSE) then
      print *, 'io_server: get_info_from_comp rank', myrank, 'NSUBSET_ITERATIONS =', NSUBSET_ITERATIONS
      print *, 'io_server: get_info_from_comp rank', myrank, 'sum(offset_nglob_cm) =', sum(offset_nglob_cm)
      print *, 'io_server: get_info_from_comp rank', myrank, 'n_msg_ford_undo =', n_msg_ford_undo
      call flush_stdout()
    endif

    ! build ordered tag list for undo-attenuation messages
    call build_undo_tag_list()

  else
    n_msg_ford_undo = 0
    NSUBSET_ITERATIONS = 0
  endif ! UNDO_ATTENUATION and SAVE_FORWARD

  ! surface movie metadata
  if (MOVIE_SURFACE) then

    ! allocate offset array on IO tasks
    if (.not. allocated(offset_poin)) then
      allocate(offset_poin(0:NPROCTOT_VAL-1))
    endif

    ! receive surface movie offsets and total number of points
    call recv_i_inter(offset_poin, NPROCTOT_VAL, 0, io_tag_surf_offset)
    call recv_i_inter(tmp_arr, 1, 0, io_tag_surf_npoints)
    npoints_surf_mov_all_proc = tmp_arr(1)

    ! receive time-stepping range for surface movie frames
    call recv_i_inter(tmp_arr, 1, 0, io_tag_surf_it_begin)
    it_begin = tmp_arr(1)
    call recv_i_inter(tmp_arr, 1, 0, io_tag_surf_it_end)
    it_end = tmp_arr(1)

  endif

  ! volume movie metadata (movie-point-based types)
  if (MOVIE_VOLUME) then

    ! allocate offset array on IO tasks
    if (.not. allocated(offset_poin_vol)) then
      allocate(offset_poin_vol(0:NPROCTOT_VAL-1))
    endif

    ! receive volume movie offsets and total number of points
    call recv_i_inter(offset_poin_vol, NPROCTOT_VAL, 0, io_tag_vol_offset)
    call recv_i_inter(tmp_arr, 1, 0, io_tag_vol_npoints)
    npoints_vol_mov_all_proc = tmp_arr(1)

  endif

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine pass_info_to_io() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'pass_info_to_io() called without HDF5 compilation support'

#endif

  end subroutine get_info_from_comp
!
!-------------------------------------------------------------------------------------------------
!

  subroutine pass_info_to_io()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer :: i_ionod
    integer :: tmp_int

    ! pass necessary information to io node

    ! forward undo att
    if (UNDO_ATTENUATION .and. SAVE_FORWARD) then
      if (myrank == 0) then

        ! count up the number of messages
        n_msg_ford_undo = 19
        if (ROTATION_VAL) n_msg_ford_undo = n_msg_ford_undo + 2
        if (ATTENUATION_VAL) n_msg_ford_undo = n_msg_ford_undo + 10
        if (FULL_GRAVITY_VAL) n_msg_ford_undo = n_msg_ford_undo + 3

        do i_ionod = 0, HDF5_IO_NODES-1

          ! send NSUBSET_ITERATIONS
          call send_i_inter((/NSUBSET_ITERATIONS/), 1, i_ionod, io_tag_nsubset_iterations)

          ! send offset arrays
          call send_i_inter(offset_nglob_cm, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_cm) !!!!!!!!!!!!!!
          call send_i_inter(offset_nglob_oc, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_oc)
          call send_i_inter(offset_nglob_ic, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_ic)
          call send_i_inter(offset_nspec_cm_soa, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_eps_xx_cm)
          call send_i_inter(offset_nspec_ic_soa, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_eps_xx_ic)
          if (ROTATION_VAL) then
            call send_i_inter(offset_nspec_oc_rot, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_A_rot)
          endif
          if (ATTENUATION_VAL) then
            call send_i_inter(offset_nspec_cm_att, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_R_xx_cm)
            call send_i_inter(offset_nspec_ic_att, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_R_xx_ic)
          endif
          if (FULL_GRAVITY_VAL) then
            call send_i_inter(offset_pgrav1, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_pgrav1)
          endif

          ! send the number of messages
          tmp_int = n_msg_ford_undo
          call send_i_inter((/tmp_int/), 1, i_ionod, io_tag_ford_undo_nmsg)

        enddo ! i_ionod
      endif ! myrank == 0

    endif ! UNDO_ATTENUATION .and. SAVE_FORWARD

    ! surface movie metadata
    if (MOVIE_SURFACE) then
      if (myrank == 0) then

        do i_ionod = 0, HDF5_IO_NODES-1

          ! send per-rank surface offsets and total number of points
          call send_i_inter(offset_poin, NPROCTOT_VAL, i_ionod, io_tag_surf_offset)
          tmp_int = npoints_surf_mov_all_proc
          call send_i_inter((/tmp_int/), 1, i_ionod, io_tag_surf_npoints)

          ! send time-stepping range for surface movie frames
          tmp_int = it_begin
          call send_i_inter((/tmp_int/), 1, i_ionod, io_tag_surf_it_begin)
          tmp_int = it_end
          call send_i_inter((/tmp_int/), 1, i_ionod, io_tag_surf_it_end)

        enddo
      endif
    endif

    ! volume movie metadata (movie-point-based types)
    if (MOVIE_VOLUME) then
      if (myrank == 0) then

        do i_ionod = 0, HDF5_IO_NODES-1

          ! send per-rank volume movie offsets and total number of movie points
          call send_i_inter(offset_poin_vol, NPROCTOT_VAL, i_ionod, io_tag_vol_offset)
          tmp_int = npoints_vol_mov_all_proc
          call send_i_inter((/tmp_int/), 1, i_ionod, io_tag_vol_npoints)

        enddo
      endif
    endif

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine pass_info_to_io() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'pass_info_to_io() called without HDF5 compilation support'

#endif

  end subroutine pass_info_to_io

!-------------------------------------------------------------------------------------------------
!

  subroutine build_undo_tag_list()

#ifdef USE_HDF5

  use specfem_par

  implicit none

  integer, parameter :: N_BASE_TAGS = 19
  integer, dimension(N_BASE_TAGS) :: base_tags
  integer, dimension(2)  :: rot_tags
  integer, dimension(10) :: att_tags
  integer, dimension(3)  :: grav_tags
  integer :: idx

  ! define base tags in the exact send order used in
  ! save_forward_arrays_undoatt_hdf5()
  base_tags = (/ &
    io_tag_ford_undo_d_cm, &
    io_tag_ford_undo_v_cm, &
    io_tag_ford_undo_a_cm, &
    io_tag_ford_undo_d_oc, &
    io_tag_ford_undo_v_oc, &
    io_tag_ford_undo_a_oc, &
    io_tag_ford_undo_d_ic, &
    io_tag_ford_undo_v_ic, &
    io_tag_ford_undo_a_ic, &
    io_tag_ford_undo_eps_xx_cm, &
    io_tag_ford_undo_eps_yy_cm, &
    io_tag_ford_undo_eps_xy_cm, &
    io_tag_ford_undo_eps_xz_cm, &
    io_tag_ford_undo_eps_yz_cm, &
    io_tag_ford_undo_eps_xx_ic, &
    io_tag_ford_undo_eps_yy_ic, &
    io_tag_ford_undo_eps_xy_ic, &
    io_tag_ford_undo_eps_xz_ic, &
    io_tag_ford_undo_eps_yz_ic /)

  rot_tags = (/ io_tag_ford_undo_A_rot, io_tag_ford_undo_B_rot /)

  att_tags = (/ &
    io_tag_ford_undo_R_xx_cm, &
    io_tag_ford_undo_R_yy_cm, &
    io_tag_ford_undo_R_xy_cm, &
    io_tag_ford_undo_R_xz_cm, &
    io_tag_ford_undo_R_yz_cm, &
    io_tag_ford_undo_R_xx_ic, &
    io_tag_ford_undo_R_yy_ic, &
    io_tag_ford_undo_R_xy_ic, &
    io_tag_ford_undo_R_xz_ic, &
    io_tag_ford_undo_R_yz_ic /)

  grav_tags = (/ io_tag_ford_undo_neq, io_tag_ford_undo_neq1, io_tag_ford_undo_pgrav1 /)

  ! total number of tags equals n_msg_ford_undo provided by rank 0
  n_undo_tags = n_msg_ford_undo

  if (allocated(undo_tag_list)) deallocate(undo_tag_list)
  if (n_undo_tags > 0) then
    allocate(undo_tag_list(n_undo_tags))
  else
    return
  endif

  ! fill in order: base, optional rotation, optional attenuation, optional gravity
  idx = 0

  undo_tag_list(1:N_BASE_TAGS) = base_tags
  idx = N_BASE_TAGS

  if (ROTATION_VAL) then
    undo_tag_list(idx+1:idx+2) = rot_tags
    idx = idx + 2
  endif

  if (ATTENUATION_VAL) then
    undo_tag_list(idx+1:idx+10) = att_tags
    idx = idx + 10
  endif

  if (FULL_GRAVITY_VAL) then
    undo_tag_list(idx+1:idx+3) = grav_tags
    idx = idx + 3
  endif

  ! safety: idx should match n_undo_tags
  if (idx /= n_undo_tags) then
    print *, 'build_undo_tag_list: mismatch between constructed tag count and n_msg_ford_undo', idx, n_undo_tags
  endif

#else

  ! no-op when built without HDF5 support

#endif

  end subroutine build_undo_tag_list

!
!-------------------------------------------------------------------------------------------------
!
! MPI communications
!
!-------------------------------------------------------------------------------------------------

  subroutine wait_all_send()

#ifdef USE_HDF5

  implicit none

  integer :: ireq

  ! forward undo arrays
  if (n_req_ford_undo /= 0) then
    ! wait till all mpi_isends are finished
    do ireq = 1,n_req_ford_undo
      call wait_req(req_dump_ford_undo(ireq))
    enddo
  endif
  n_req_ford_undo = 0

  ! surface movie
  if (n_req_surf /= 0) then
    do ireq = 1,n_req_surf
      call wait_req(req_dump_surf(ireq))
    enddo
  endif
  n_req_surf = 0

  ! volume movie
  if (n_req_vol /= 0) then
    do ireq = 1,n_req_vol
      call wait_req(req_dump_vol(ireq))
    enddo
  endif
  n_req_vol = 0

  call synchronize_all()

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine wait_all_send() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'wait_all_send() called without HDF5 compilation support'

#endif

  end subroutine wait_all_send

!
!-------------------------------------------------------------------------------------------------
!
  subroutine allocate_offset_arrays()

#ifdef USE_HDF5
    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer :: ier

    ! nglobs
    allocate(offset_nglob_cm(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nglob_cm')
    allocate(offset_nglob_oc(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nglob_oc')
    allocate(offset_nglob_ic(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nglob_ic')

    ! nelems
    allocate(offset_nspec_cm(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_cm')
    allocate(offset_nspec_oc(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_oc')
    allocate(offset_nspec_ic(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_ic')
    allocate(offset_nspec_cm_soa(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_cm_soa')
    allocate(offset_nspec_ic_soa(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_ic_soa')
    if (ROTATION_VAL) then
      allocate(offset_nspec_oc_rot(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_oc_rot')
    endif
    if (ATTENUATION_VAL) then
      allocate(offset_nspec_cm_att(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_cm_att')
      allocate(offset_nspec_ic_att(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_nspec_ic_att')
    endif

    if (FULL_GRAVITY_VAL) then
      allocate(offset_pgrav1(0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_pgrav1')
    endif

    npoints_vol_mov_all_proc_cm     = 0
    npoints_vol_mov_all_proc_oc     = 0
    npoints_vol_mov_all_proc_ic     = 0

    nspec_vol_mov_all_proc_cm     = 0
    nspec_vol_mov_all_proc_oc     = 0
    nspec_vol_mov_all_proc_ic     = 0
    nspec_vol_mov_all_proc_cm_soa = 0
    nspec_vol_mov_all_proc_ic_soa = 0
    if (ROTATION_VAL) then
      nspec_vol_mov_all_proc_oc_rot = 0
    endif
    if (ATTENUATION_VAL) then
      nspec_vol_mov_all_proc_cm_att = 0
      nspec_vol_mov_all_proc_ic_att = 0
    endif

#else

  print *
  print *, "Error: HDF5 I/O server routine allocate_offset_arrays() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'allocate_offset_arrays() called without HDF5 compilation support'

#endif

  end subroutine allocate_offset_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine deallocate_offset_arrays()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    deallocate(offset_nglob_cm)
    deallocate(offset_nglob_oc)
    deallocate(offset_nglob_ic)

    deallocate(offset_nspec_cm)
    deallocate(offset_nspec_oc)
    deallocate(offset_nspec_ic)
    deallocate(offset_nspec_cm_soa)
    deallocate(offset_nspec_ic_soa)
    if (ROTATION_VAL) then
      deallocate(offset_nspec_oc_rot)
    endif
    if (ATTENUATION_VAL) then
      deallocate(offset_nspec_cm_att)
      deallocate(offset_nspec_ic_att)
    endif
    if (FULL_GRAVITY_VAL) then
      deallocate(offset_pgrav1)
    endif

#else

  print *
  print *, "Error: HDF5 I/O server routine deallocate_offset_arrays() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'deallocate_offset_arrays() called without HDF5 compilation support'

#endif

  end subroutine deallocate_offset_arrays

!
!-------------------------------------------------------------------------------------------------
!



  subroutine initialize_hdf5_solver()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5
    use specfem_par_full_gravity, only : neq1

    implicit none

    ! allocate offset arrays
    call allocate_offset_arrays()

    call gather_all_all_singlei(NGLOB_CRUST_MANTLE, offset_nglob_cm, NPROCTOT_VAL)
    call gather_all_all_singlei(NGLOB_OUTER_CORE,   offset_nglob_oc, NPROCTOT_VAL)
    call gather_all_all_singlei(NGLOB_INNER_CORE,   offset_nglob_ic, NPROCTOT_VAL)

    call gather_all_all_singlei(NSPEC_CRUST_MANTLE, offset_nspec_cm, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_OUTER_CORE,   offset_nspec_oc, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_INNER_CORE,   offset_nspec_ic, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_CRUST_MANTLE_STR_OR_ATT, offset_nspec_cm_soa, NPROCTOT_VAL)
    call gather_all_all_singlei(NSPEC_INNER_CORE_STR_OR_ATT,   offset_nspec_ic_soa, NPROCTOT_VAL)
    if (ROTATION_VAL) then
      call gather_all_all_singlei(NSPEC_OUTER_CORE_ROTATION, offset_nspec_oc_rot, NPROCTOT_VAL)
    endif
    if (ATTENUATION_VAL) then
      call gather_all_all_singlei(NSPEC_CRUST_MANTLE_ATTENUATION, offset_nspec_cm_att, NPROCTOT_VAL)
      call gather_all_all_singlei(NSPEC_INNER_CORE_ATTENUATION,   offset_nspec_ic_att, NPROCTOT_VAL)
    endif
    if (FULL_GRAVITY_VAL) then
      call gather_all_all_singlei(neq1, offset_pgrav1, NPROCTOT_VAL)
    endif

    npoints_vol_mov_all_proc_cm = sum(offset_nglob_cm)
    npoints_vol_mov_all_proc_oc = sum(offset_nglob_oc)
    npoints_vol_mov_all_proc_ic = sum(offset_nglob_ic)

    nspec_vol_mov_all_proc_cm = sum(offset_nspec_cm)
    nspec_vol_mov_all_proc_oc = sum(offset_nspec_oc)
    nspec_vol_mov_all_proc_ic = sum(offset_nspec_ic)
    nspec_vol_mov_all_proc_cm_soa = sum(offset_nspec_cm_soa)
    nspec_vol_mov_all_proc_ic_soa = sum(offset_nspec_ic_soa)
    if (ROTATION_VAL) then
      nspec_vol_mov_all_proc_oc_rot = sum(offset_nspec_oc_rot)
    endif
    if (ATTENUATION_VAL) then
      nspec_vol_mov_all_proc_cm_att = sum(offset_nspec_cm_att)
      nspec_vol_mov_all_proc_ic_att = sum(offset_nspec_ic_att)
    endif

#else

  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine initialize_hdf5_solver() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'initialize_hdf5_solver() called without HDF5 compilation support'

#endif

  end subroutine initialize_hdf5_solver

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_hdf5_solver()

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    ! deallocate offset arrays
    call deallocate_offset_arrays()

#else

  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine finalize_hdf5_solver() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'finalize_hdf5_solver() called without HDF5 compilation support'

#endif

  end subroutine finalize_hdf5_solver

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_hdf5_files_and_datasets(i_snapshot)


#ifdef USE_HDF5
    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
#endif

    implicit none

    integer, intent(in) :: i_snapshot

#ifdef USE_HDF5
    ! forward undo arrays

    ! In multi-IO mode, undo checkpoints are written into per-IO
    ! shard files by recv_and_write_ford_undo(). To avoid creating
    ! an unused single shared save_frame_at*.h5 here, skip global
    ! pre-creation when HDF5_IO_NODES > 1.
    if (HDF5_IO_NODES > 1) then
      return
    endif

    if (UNDO_ATTENUATION .and. SAVE_FORWARD) then

      ! construct output file name for this snapshot
      write(file_name, '(a,i6.6,a)') 'save_frame_at',i_snapshot,'.h5'
      file_name = trim(LOCAL_PATH)//'/'//trim(file_name)

      if (VERBOSE) then
        print *, 'io_server: create_hdf5_files_and_datasets rank', myrank, &
                 ' snapshot', i_snapshot, ' file =', trim(file_name)
        print *, 'io_server: create_hdf5_files_and_datasets rank', myrank, &
                 ' sum(offset_nglob_cm) =', sum(offset_nglob_cm)
        call flush_stdout()
      endif

      ! all IO ranks collectively create the file with the parallel driver
      ! (MPI communicator and info are already configured via h5_set_mpi_info)
      call h5_create_file_p_collect(file_name)

      ! create datasets (collective metadata operations across IO ranks)
      call h5_create_dataset_gen('displ_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
      call h5_create_dataset_gen('veloc_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
      call h5_create_dataset_gen('accel_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
      call h5_create_dataset_gen('displ_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('veloc_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('accel_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('displ_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
      call h5_create_dataset_gen('veloc_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
      call h5_create_dataset_gen('accel_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xx_crust_mantle', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_yy_crust_mantle', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xy_crust_mantle', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xz_crust_mantle', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_yz_crust_mantle', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xx_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_yy_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xy_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xz_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_yz_inner_core', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)

      if (ROTATION_VAL) then
        call h5_create_dataset_gen('A_array_rotation', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_oc_rot)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('B_array_rotation', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_oc_rot)/), 4, CUSTOM_REAL)
      endif

      if (ATTENUATION_VAL) then
        call h5_create_dataset_gen('R_xx_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_yy_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_xy_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_xz_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_yz_crust_mantle', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)

        call h5_create_dataset_gen('R_xx_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_yy_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_xy_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_xz_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        call h5_create_dataset_gen('R_yz_inner_core', (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
      endif ! ATTENUATION_VAL

      if (FULL_GRAVITY_VAL) then
        call h5_create_dataset_gen('neq', (/NPROCTOT_VAL/), 1, 1)
        call h5_create_dataset_gen('neq1', (/NPROCTOT_VAL/), 1, 1)
        call h5_create_dataset_gen('pgrav1', (/sum(offset_pgrav1)/), 1, CUSTOM_REAL)
      endif ! FULL_GRAVITY_VAL

      ! close file collectively
      call h5_close_file_p()

    endif ! UNDO_ATTENUATION .and. SAVE_FORWARD

    ! global barrier among IO tasks before proceeding to next snapshot
    call synchronize_all()

#else

  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine create_files_and_datasets_ford_undo() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'create_files_and_datasets_ford_undo() called without HDF5 compilation support'

#endif

  end subroutine create_hdf5_files_and_datasets


!-------------------------------------------------------------------------------
!
! HDF5 io routines (only available with HDF5 compilation support)
!
!-------------------------------------------------------------------------------

#if defined(USE_HDF5)
  ! only available with HDF5 compilation support

  subroutine idle_mpi_io(status)
  ! wait for an arrival of any MPI message

    use constants, only: my_status_size

    implicit none

    integer, intent(inout) :: status(my_status_size)

    call world_probe_any_inter(status)

  end subroutine idle_mpi_io

#endif
  !
  !-------------------------------------------------------------------------------------------------
  !

  subroutine create_surface_frame_group(i_frame, it_first_surf)

#ifdef USE_HDF5

  use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
  use constants, only: myrank

    implicit none

    integer, intent(in) :: i_frame
    integer, intent(in) :: it_first_surf

    integer :: it_val
    integer :: i_io

    ! safety: this routine must only be called on rank 0 of the world
    if (myrank > 0) then
      print *, 'Error: create_surface_frame_group called from rank ', myrank
      stop 'create_surface_frame_group must be called only on rank 0'
    endif

    ! compute actual time-step index for this frame
    it_val = it_first_surf + i_frame * NTSTEP_BETWEEN_FRAMES

    ! construct group name for this time step
    group_name = "it_"//trim(i2c(it_val))

    if (HDF5_IO_NODES > 1) then
      ! multi-IO-server mode: do not touch shard files here. Each IO
      ! server will lazily create/open its own shard file and frame
      ! group when it first receives data for this frame. This avoids
      ! any cross-rank contention or file locking issues.
      return
    else
      ! single-file mode: original behaviour creating the frame in
      ! OUTPUT_FILES/movie_surface.h5
      file_name = trim(OUTPUT_FILES)//"/movie_surface.h5"
      call h5_create_or_open_file(file_name)
      call h5_create_group(group_name)
      call h5_open_group(group_name)
      call h5_create_dataset_gen_in_group("ux", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uy", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uz", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      call h5_close_group()
      call h5_close_file()
    endif

#endif

  end subroutine create_surface_frame_group

  !-------------------------------------------------------------------------------------------------
  !

  subroutine create_volume_frame_group(i_frame, it_first_vol)

#ifdef USE_HDF5

  use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
  use constants, only: myrank

    implicit none

    integer, intent(in) :: i_frame
    integer, intent(in) :: it_first_vol

    integer :: it_val
    character(len=2) :: movie_prefix2
    logical :: dset_exists
    character(len=MAX_STRING_LEN) :: dset_full_name

    ! safety: this routine must only be called on rank 0
    if (myrank > 0) then
      print *, 'Error: create_volume_frame_group called from rank ', myrank
      stop 'create_volume_frame_group must be called only on rank 0'
    endif

    ! in multi-IO mode, we let each IO rank lazily create its own
    ! shard file and frame group when it first receives data for
    ! this frame. Avoids cross-rank contention on a single file.
    if (HDF5_IO_NODES > 1) then
      return
    endif

    ! compute actual time-step index for this frame
    it_val = it_first_vol + i_frame * NTSTEP_BETWEEN_FRAMES

    ! construct file and group names
    file_name = trim(OUTPUT_FILES)//"/movie_volume.h5"
    group_name = "it_"//trim(i2c(it_val))

    ! open file in serial mode and create group/datasets
    call h5_open_file(file_name)
    call h5_open_or_create_group(group_name)

    select case (MOVIE_VOLUME_TYPE)
    case (1,2,3)
      ! strains / time-integrated / potency at movie points
      if (MOVIE_VOLUME_TYPE == 1) then
        movie_prefix2 = 'E '
      else if (MOVIE_VOLUME_TYPE == 2) then
        movie_prefix2 = 'S '
      else
        movie_prefix2 = 'P '
      endif

      ! create datasets only if they do not already exist
      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'NN'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'NN', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'EE'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'EE', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'ZZ'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'ZZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'NE'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'NE', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'NZ'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'NZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'EZ'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'EZ', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

    case (5,6)
      ! displacement / velocity vectors at movie points
      if (MOVIE_VOLUME_TYPE == 5) then
        movie_prefix2 = 'DI'
      else
        movie_prefix2 = 'VE'
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'N'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'N', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'E'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'E', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/'//trim(movie_prefix2)//'Z'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group(trim(movie_prefix2)//'Z', (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
      endif
    end select

    call h5_close_group()
    call h5_close_file()

#endif

  end subroutine create_volume_frame_group

  !-------------------------------------------------------------------------------------------------
  !

#if defined(USE_HDF5)

  subroutine recv_and_write_ford_undo(tag, tag_src, status, &
                                     dump_ford_undo_1d_glob, &
                                     dump_ford_undo_2d_glob, &
                                     dump_ford_undo_4d, &
                                     dump_ford_undo_5d, &
                                     i_snapshot)

    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
    use constants, only: CUSTOM_REAL, my_status_size, my_status_source, my_status_tag
    use io_bandwidth

    implicit none

    integer, intent(in) :: tag, tag_src
    integer, intent(in) :: status(my_status_size)
    real(kind=CUSTOM_REAL), dimension(:),         intent(inout) :: dump_ford_undo_1d_glob
    real(kind=CUSTOM_REAL), dimension(:,:),       intent(inout) :: dump_ford_undo_2d_glob
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),   intent(inout) :: dump_ford_undo_4d
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(inout) :: dump_ford_undo_5d
    integer, intent(in) :: i_snapshot

    integer :: msg_size, ista, iend, req_dummy, data_len
    integer :: neq, neq1

    integer :: data_size ! size of one data element in bytes
    logical :: use_collective
    logical :: dset_exists

    ! default datasize is for CUSTOM_REAL
    data_size = CUSTOM_REAL ! 8 for double precision, 4 for single precision

    ! use collective HDF5 only when there is a single IO node
    ! multi-IO with asynchronous message arrival breaks collective ordering
    use_collective = H5_COL .and. (HDF5_IO_NODES <= 1)

    ! get message size
    call world_get_size_msg(status, msg_size)

    ! file name to write
    !  - single-IO mode:  save_frame_atXXXXXX.h5
    !  - multi-IO mode:   save_frame_atXXXXXX.io<my_io_id>.h5
    write(file_name, '(a,i6.6,a)') 'save_frame_at',i_snapshot+1,'.h5'
    if (HDF5_IO_NODES > 1) then
      file_name = trim(LOCAL_PATH)//'/save_frame_at'//trim(i2c(i_snapshot+1))//'.io'//trim(i2c(my_io_id))//'.h5'
    else
      file_name = trim(LOCAL_PATH)//'/'//trim(file_name)
    endif

    if (VERBOSE) then
      print *, 'io_server: recv_and_write_ford_undo rank', myrank, &
               ' snapshot', i_snapshot+1, ' tag', tag, ' from src', tag_src, &
               ' open file =', trim(file_name)
      call flush_stdout()
    endif

    ! open file
    !! get MPI parameters
    !call world_get_comm(comm)
    !call world_get_info_null(info)

    !! initialize HDF5
    !call h5_initialize() ! called in initialize_mesher()
    !! set MPI
    !call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

    if (HDF5_IO_NODES > 1) then
      ! multi-IO mode: each IO rank independently opens or creates
      ! its own shard file using serial HDF5
      call h5_create_or_open_file(file_name)
    else
      if (use_collective) then
        ! open file using parallel HDF5 collectively
        call h5_open_file_p_collect(file_name)
      else
        ! open file using parallel HDF5 without collectives
        call h5_open_file_p(file_name)
      endif
    endif

    ! receive the data
    if (tag == io_tag_ford_undo_d_cm) then
      ! displ_cust_mantle
      ista = sum(offset_nglob_cm(0:tag_src-1))
      iend = sum(offset_nglob_cm(0:tag_src))
      data_len = offset_nglob_cm(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('displ_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('displ_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('displ_crust_mantle', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_v_cm) then
      ! veloc_crust_mantle
      ista = sum(offset_nglob_cm(0:tag_src-1))
      iend = sum(offset_nglob_cm(0:tag_src))
      data_len = offset_nglob_cm(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (VERBOSE) then
        print *, 'io_server: writing veloc_crust_mantle on rank', myrank, &
                 ' snapshot', i_snapshot+1, ' ista =', ista, ' data_len =', data_len, &
                 ' use_collective =', use_collective
        call flush_stdout()
      endif
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('veloc_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('veloc_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('veloc_crust_mantle', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_a_cm) then
      ! accel_crust_mantle
      ista = sum(offset_nglob_cm(0:tag_src-1))
      iend = sum(offset_nglob_cm(0:tag_src))
      data_len = offset_nglob_cm(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('accel_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('accel_crust_mantle', (/NDIM, sum(offset_nglob_cm)/), 2, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('accel_crust_mantle', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_d_oc) then
      ! displ_outer_core
      ista = sum(offset_nglob_oc(0:tag_src-1))
      iend = sum(offset_nglob_oc(0:tag_src))
      data_len = offset_nglob_oc(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('displ_outer_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('displ_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('displ_outer_core', dump_ford_undo_1d_glob(1:data_len), (/ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_v_oc) then
      ! veloc_outer_core
      ista = sum(offset_nglob_oc(0:tag_src-1))
      iend = sum(offset_nglob_oc(0:tag_src))
      data_len = offset_nglob_oc(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('veloc_outer_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('veloc_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('veloc_outer_core', dump_ford_undo_1d_glob(1:data_len), (/ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_a_oc) then
      ! accel_outer_core
      ista = sum(offset_nglob_oc(0:tag_src-1))
      iend = sum(offset_nglob_oc(0:tag_src))
      data_len = offset_nglob_oc(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('accel_outer_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('accel_outer_core', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('accel_outer_core', dump_ford_undo_1d_glob(1:data_len), (/ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_d_ic) then
      ! displ_inner_core
      ista = sum(offset_nglob_ic(0:tag_src-1))
      iend = sum(offset_nglob_ic(0:tag_src))
      data_len = offset_nglob_ic(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('displ_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('displ_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('displ_inner_core', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_v_ic) then
      ! veloc_inner_core
      ista = sum(offset_nglob_ic(0:tag_src-1))
      iend = sum(offset_nglob_ic(0:tag_src))
      data_len = offset_nglob_ic(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('veloc_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('veloc_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('veloc_inner_core', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_a_ic) then
      ! accel_inner_core
      ista = sum(offset_nglob_ic(0:tag_src-1))
      iend = sum(offset_nglob_ic(0:tag_src))
      data_len = offset_nglob_ic(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_2d_glob(:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('accel_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('accel_inner_core', (/NDIM, sum(offset_nglob_ic)/), 2, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('accel_inner_core', dump_ford_undo_2d_glob(:,1:data_len), (/0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xx_cm) then
      ! epsilondev_xx_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_xx_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_xx_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                     (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yy_cm) then
      ! epsilondev_yy_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_yy_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_yy_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                     (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xy_cm) then
      ! epsilondev_xy_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_xy_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_xy_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                     (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xz_cm) then
      ! epsilondev_xz_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_xz_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_xz_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                     (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yz_cm) then
      ! epsilondev_yz_crust_mantle
      ista = sum(offset_nspec_cm_soa(0:tag_src-1))
      iend = sum(offset_nspec_cm_soa(0:tag_src))
      data_len = offset_nspec_cm_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_yz_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_yz_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', dump_ford_undo_4d(:,:,:,1:data_len), &
                                     (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xx_ic) then
      ! epsilondev_xx_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_xx_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_xx_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xx_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                   (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yy_ic) then
      ! epsilondev_yy_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_yy_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_yy_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yy_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                   (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xy_ic) then
      ! epsilondev_xy_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_xy_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_xy_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xy_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                   (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_xz_ic) then
      ! epsilondev_xz_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_xz_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_xz_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_xz_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                   (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_eps_yz_ic) then
      ! epsilondev_yz_inner_core
      ista = sum(offset_nspec_ic_soa(0:tag_src-1))
      iend = sum(offset_nspec_ic_soa(0:tag_src))
      data_len = offset_nspec_ic_soa(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('epsilondev_yz_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('epsilondev_yz_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('epsilondev_yz_inner_core', dump_ford_undo_4d(:,:,:,1:data_len), &
                                   (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_A_rot) then
      ! A_array_rotation
      ista = sum(offset_nspec_oc_rot(0:tag_src-1))
      iend = sum(offset_nspec_oc_rot(0:tag_src))
      data_len = offset_nspec_oc_rot(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('A_array_rotation', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('A_array_rotation', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_oc_rot)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('A_array_rotation', dump_ford_undo_4d(:,:,:,1:data_len), &
                                (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_B_rot) then
      ! B_array_rotation
      ista = sum(offset_nspec_oc_rot(0:tag_src-1))
      iend = sum(offset_nspec_oc_rot(0:tag_src))
      data_len = offset_nspec_oc_rot(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_4d(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('B_array_rotation', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('B_array_rotation', &
                                     (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_oc_rot)/), 4, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('B_array_rotation', dump_ford_undo_4d(:,:,:,1:data_len), &
                                (/0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xx_cm) then
      ! R_xx_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_xx_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_xx_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xx_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yy_cm) then
      ! R_yy_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_yy_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_yy_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yy_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xy_cm) then
      ! R_xy_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_xy_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_xy_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xy_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xz_cm) then
      ! R_xz_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_xz_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_xz_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xz_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yz_cm) then
      ! R_yz_crust_mantle
      ista = sum(offset_nspec_cm_att(0:tag_src-1))
      iend = sum(offset_nspec_cm_att(0:tag_src))
      data_len = offset_nspec_cm_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_yz_crust_mantle', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_yz_crust_mantle', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_cm_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yz_crust_mantle', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xx_ic) then
      ! R_xx_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_xx_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_xx_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xx_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                  (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yy_ic) then
      ! R_yy_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_yy_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_yy_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yy_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                  (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xy_ic) then
      ! R_xy_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_xy_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_xy_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xy_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                  (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_xz_ic) then
      ! R_xz_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_xz_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_xz_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_xz_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                  (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_R_yz_ic) then
      ! R_yz_inner_core
      ista = sum(offset_nspec_ic_att(0:tag_src-1))
      iend = sum(offset_nspec_ic_att(0:tag_src))
      data_len = offset_nspec_ic_att(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_5d(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('R_yz_inner_core', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('R_yz_inner_core', &
                                     (/NGLLX, NGLLY, NGLLZ, N_SLS, sum(offset_nspec_ic_att)/), 5, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('R_yz_inner_core', dump_ford_undo_5d(:,:,:,:,1:data_len), &
                                  (/0, 0, 0, 0, ista/), use_collective)
      call stop_timer()

    else if (tag == io_tag_ford_undo_neq) then
      ! neq
      ! receive
      call irecv_i_inter((/neq/), 1, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('neq', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('neq', (/NPROCTOT_VAL/), 1, 1)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('neq', (/neq/), (/tag_src/), use_collective)
      call stop_timer()

      ! data size if for integer
      data_size = 4

    else if (tag == io_tag_ford_undo_neq1) then
      ! neq1
      ! receive
      call irecv_i_inter((/neq1/), 1, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('neq1', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('neq1', (/NPROCTOT_VAL/), 1, 1)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('neq1', (/neq1/), (/tag_src/), use_collective)
      call stop_timer()

      ! data size if for integer
      data_size = 4

    else if (tag == io_tag_ford_undo_pgrav1) then
      ! pgrav1
      ista = sum(offset_pgrav1(0:tag_src-1))
      iend = sum(offset_pgrav1(0:tag_src))
      data_len = offset_pgrav1(tag_src)
      ! receive
      call irecvv_cr_inter(dump_ford_undo_1d_glob(1:data_len), msg_size, tag_src, tag, req_dummy)
      ! write
      if (HDF5_IO_NODES > 1) then
        call h5_check_dataset_exists('pgrav1', dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen('pgrav1', (/sum(offset_pgrav1)/), 1, CUSTOM_REAL)
        endif
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab('pgrav1', dump_ford_undo_1d_glob(1:data_len), (/ista/), use_collective)
      call stop_timer()
    else
      ! unknown tag
      print *, 'Error: unknown tag in recv_and_write_ford_undo'
      stop 'Error: unknown tag in recv_and_write_ford_undo'

    end if


    ! count the bytes written
    call set_bytes_written_from_array(data_size*8, msg_size) ! converting to bits from bytes

    ! close file
    if (HDF5_IO_NODES > 1) then
      call h5_close_file()
    else
      call h5_close_file_p()
    endif

  end subroutine recv_and_write_ford_undo

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_and_write_surface_movie(tag, tag_src, status, &
                                          dump_surf_ux, dump_surf_uy, dump_surf_uz, &
                                          i_frame, it_first_surf)

    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
    use constants, only: CUSTOM_REAL, my_status_size
    use io_bandwidth

    implicit none

    integer, intent(in) :: tag, tag_src
    integer, intent(in) :: status(my_status_size)
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_surf_ux
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_surf_uy
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_surf_uz
    integer, intent(in) :: i_frame
    integer, intent(in) :: it_first_surf

    integer :: msg_size, ista, data_len, req_dummy
    integer :: it_val
    integer :: data_size
    logical :: dset_exists
    character(len=MAX_STRING_LEN) :: dset_full_name
    logical :: use_collective

    ! get message size
    call world_get_size_msg(status, msg_size)

    ! size of one data element in bytes
    data_size = CUSTOM_REAL

    ! use collective HDF5 only when there is a single IO node writing
    use_collective = H5_COL .and. (HDF5_IO_NODES <= 1)

    ! compute actual time-step index for this frame
    it_val = it_first_surf + i_frame * NTSTEP_BETWEEN_FRAMES

    ! construct file and group names
    if (HDF5_IO_NODES > 1) then
      ! multi-IO-server mode: each IO server writes to its own shard file
      file_name = trim(OUTPUT_FILES)//"/movie_surface.io"//trim(i2c(my_io_id))//".h5"
    else
      ! single IO server or no IO sharding: original movie_surface.h5 file
      file_name = trim(OUTPUT_FILES)//"/movie_surface.h5"
    endif
    group_name = "it_"//trim(i2c(it_val))

    if (VERBOSE) then
      print *, 'io_server: recv_and_write_surface_movie rank', myrank, &
               ' frame', i_frame+1, ' it =', it_val, ' tag', tag, ' src', tag_src, &
               ' file =', trim(file_name), ' group =', trim(group_name)
      call flush_stdout()
    endif

    ! open file
    if (HDF5_IO_NODES > 1) then
      ! each IO rank independently opens or creates its own shard file
      call h5_create_or_open_file(file_name)

      ! ensure group and datasets exist locally in this shard file
      call h5_open_or_create_group(group_name)

      dset_full_name = trim(group_name)//'/ux'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group('ux', (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/uy'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group('uy', (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      dset_full_name = trim(group_name)//'/uz'
      call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
      if (.not. dset_exists) then
        call h5_create_dataset_gen_in_group('uz', (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      endif

      ! stay in this group for writing; no cross-rank collectives here
    else
      if (use_collective) then
        call h5_open_file_p_collect(file_name)
      else
        call h5_open_file_p(file_name)
      endif

      ! group and datasets were already created collectively
      call h5_open_group(group_name)
    endif

    if (VERBOSE) then
      print *, 'io_server: opened surface group on rank', myrank, ' group =', trim(group_name)
      call flush_stdout()
    endif

    ! receive and write according to tag
    ista = sum(offset_poin(0:tag_src-1))
    data_len = offset_poin(tag_src)

    if (tag == io_tag_surf_ux) then
      call irecvv_cr_inter(dump_surf_ux(1:data_len), msg_size, tag_src, tag, req_dummy)
      if (VERBOSE) then
        print *, 'io_server: writing surf ux on rank', myrank, ' frame', i_frame+1, &
                 ' ista =', ista, ' len =', data_len, ' use_collective =', use_collective
        call flush_stdout()
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab_in_group("ux", dump_surf_ux(1:data_len), (/ista/), use_collective)
      call stop_timer()
    else if (tag == io_tag_surf_uy) then
      call irecvv_cr_inter(dump_surf_uy(1:data_len), msg_size, tag_src, tag, req_dummy)
      if (VERBOSE) then
        print *, 'io_server: writing surf uy on rank', myrank, ' frame', i_frame+1, &
                 ' ista =', ista, ' len =', data_len, ' use_collective =', use_collective
        call flush_stdout()
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab_in_group("uy", dump_surf_uy(1:data_len), (/ista/), use_collective)
      call stop_timer()
    else if (tag == io_tag_surf_uz) then
      call irecvv_cr_inter(dump_surf_uz(1:data_len), msg_size, tag_src, tag, req_dummy)
      if (VERBOSE) then
        print *, 'io_server: writing surf uz on rank', myrank, ' frame', i_frame+1, &
                 ' ista =', ista, ' len =', data_len, ' use_collective =', use_collective
        call flush_stdout()
      endif
      call start_timer()
      call h5_write_dataset_collect_hyperslab_in_group("uz", dump_surf_uz(1:data_len), (/ista/), use_collective)
      call stop_timer()
    else
      print *, 'Error: unknown surface movie tag in recv_and_write_surface_movie'
      stop 'Error: unknown surface movie tag in recv_and_write_surface_movie'
    end if

    ! count bytes written for this surface movie message
    call set_bytes_written_from_array(data_size*8, msg_size)

    ! close group and file
    if (HDF5_IO_NODES > 1) then
      call h5_close_group()
      call h5_close_file()
    else
      call h5_close_group()
      call h5_close_file_p()
    endif

  end subroutine recv_and_write_surface_movie

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_and_write_volume_movie(tag, tag_src, status, &
                                         dump_vol1, dump_vol2, dump_vol3, &
                                         dump_vol4, dump_vol5, dump_vol6, &
                                         i_frame, it_first_vol)

    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
    use constants, only: CUSTOM_REAL, my_status_size
    use io_bandwidth

    implicit none

    integer, intent(in) :: tag, tag_src
    integer, intent(in) :: status(my_status_size)
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol1
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol2
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol3
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol4
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol5
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol6
    integer, intent(in) :: i_frame
    integer, intent(in) :: it_first_vol

    integer :: msg_size, ista, data_len, req_dummy
    integer :: it_val
    character(len=2) :: movie_prefix2
    character(len=MAX_STRING_LEN) :: dset_name
    character(len=MAX_STRING_LEN) :: dset_full_name
    integer :: data_size
    logical :: use_collective
    logical :: dset_exists

    ! get message size
    call world_get_size_msg(status, msg_size)

    ! size of one data element in bytes
    data_size = CUSTOM_REAL

    ! use collective HDF5 only when there is a single IO node
    use_collective = H5_COL .and. (HDF5_IO_NODES <= 1)

    ! compute actual time-step index for this frame
    it_val = it_first_vol + i_frame * NTSTEP_BETWEEN_FRAMES

    ! construct file and group names
    if (HDF5_IO_NODES > 1) then
      ! multi-IO-server mode: each IO server writes to its own shard
      file_name = trim(OUTPUT_FILES)//"/movie_volume.io"//trim(i2c(my_io_id))//".h5"
    else
      ! single IO server or no IO sharding: original single file
      file_name = trim(OUTPUT_FILES)//"/movie_volume.h5"
    endif
    group_name = "it_"//trim(i2c(it_val))

    if (VERBOSE) then
      print *, 'io_server: recv_and_write_volume_movie rank', myrank, &
               ' frame', i_frame+1, ' it =', it_val, ' tag', tag, ' src', tag_src, &
               ' file =', trim(file_name), ' group =', trim(group_name)
      call flush_stdout()
    endif

    ! open file
    if (HDF5_IO_NODES > 1) then
      ! each IO rank independently opens or creates its own shard file
      call h5_create_or_open_file(file_name)
      ! and lazily creates the frame group inside it
      call h5_open_or_create_group(group_name)
    else
      if (use_collective) then
        call h5_open_file_p_collect(file_name)
      else
        call h5_open_file_p(file_name)
      endif

      ! open group collectively (group and datasets have been created beforehand)
      call h5_open_group(group_name)
    endif

    if (VERBOSE) then
      print *, 'io_server: opened volume group on rank', myrank, ' group =', trim(group_name)
      call flush_stdout()
    endif

    ista = sum(offset_poin_vol(0:tag_src-1))
    data_len = offset_poin_vol(tag_src)

    select case (MOVIE_VOLUME_TYPE)
    case (1,2,3)
      if (MOVIE_VOLUME_TYPE == 1) then
        movie_prefix2 = 'E '
      else if (MOVIE_VOLUME_TYPE == 2) then
        movie_prefix2 = 'S '
      else
        movie_prefix2 = 'P '
      endif

      if (tag == io_tag_vol_strain_NN) then
        call irecvv_cr_inter(dump_vol1(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'NN'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        if (VERBOSE) then
          print *, 'io_server: writing vol', trim(dset_name), ' on rank', myrank, ' frame', i_frame+1, &
                   ' ista =', ista, ' len =', data_len, ' use_collective =', use_collective
          call flush_stdout()
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol1(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else if (tag == io_tag_vol_strain_EE) then
        call irecvv_cr_inter(dump_vol2(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'EE'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol2(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else if (tag == io_tag_vol_strain_ZZ) then
        call irecvv_cr_inter(dump_vol3(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'ZZ'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol3(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else if (tag == io_tag_vol_strain_NE) then
        call irecvv_cr_inter(dump_vol4(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'NE'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol4(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else if (tag == io_tag_vol_strain_NZ) then
        call irecvv_cr_inter(dump_vol5(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'NZ'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol5(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else if (tag == io_tag_vol_strain_EZ) then
        call irecvv_cr_inter(dump_vol6(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'EZ'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol6(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else
        print *, 'Error: unknown volume strain tag in recv_and_write_volume_movie'
        stop 'Error: unknown volume strain tag in recv_and_write_volume_movie'
      end if

    case (5,6)
      if (MOVIE_VOLUME_TYPE == 5) then
        movie_prefix2 = 'DI'
      else
        movie_prefix2 = 'VE'
      endif

      if (tag == io_tag_vol_vec_N) then
        call irecvv_cr_inter(dump_vol1(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'N'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol1(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else if (tag == io_tag_vol_vec_N) then
        call irecvv_cr_inter(dump_vol1(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'N'
        if (VERBOSE) then
          print *, 'io_server: writing vol', trim(dset_name), ' on rank', myrank, ' frame', i_frame+1, &
                   ' ista =', ista, ' len =', data_len, ' use_collective =', use_collective
          call flush_stdout()
        endif
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol1(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else if (tag == io_tag_vol_vec_Z) then
        call irecvv_cr_inter(dump_vol3(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'Z'
        if (HDF5_IO_NODES > 1) then
          dset_full_name = trim(group_name)//'/'//trim(dset_name)
          call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
          if (.not. dset_exists) then
            call h5_create_dataset_gen_in_group(trim(dset_name), (/npoints_vol_mov_all_proc/), 1, CUSTOM_REAL)
          endif
        endif
        call start_timer()
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol3(1:data_len), (/ista/), use_collective)
        call stop_timer()
      else
        print *, 'Error: unknown volume vector tag in recv_and_write_volume_movie'
        stop 'Error: unknown volume vector tag in recv_and_write_volume_movie'
      end if

    case default
      ! other volume movie types are currently not handled by IO server
      print *, 'Error: MOVIE_VOLUME_TYPE not supported in recv_and_write_volume_movie'
      stop 'Error: MOVIE_VOLUME_TYPE not supported in recv_and_write_volume_movie'
    end select

    ! count bytes written for this volume movie message
    call set_bytes_written_from_array(data_size*8, msg_size)

    ! close group and file
    call h5_close_group()
    if (HDF5_IO_NODES > 1) then
      call h5_close_file()
    else
      call h5_close_file_p()
    endif

  end subroutine recv_and_write_volume_movie

!
!-------------------------------------------------------------------------------------------------
!

#endif

end module io_server_hdf5
