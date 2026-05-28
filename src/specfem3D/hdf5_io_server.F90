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

#ifdef USE_HDF5

  implicit none

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

  ! volume movie data (norm values for MOVIE_VOLUME_TYPE 7,8,9: displnorm, velnorm, accelnorm)
  public :: io_tag_vol_norm_cm, io_tag_vol_norm_oc, io_tag_vol_norm_ic

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

  !---------------------------------------------------------------------------
  ! MPI TAGS FOR IO SERVER IMPLEMENTATION
  !---------------------------------------------------------------------------
  ! Tags are organized by functional category:
  ! - 100001-100036: Undo attenuation (forward snapshots)
  ! - 110001-110021: Surface movie
  ! - 120001-120022: Volume movie
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! UNDO ATTENUATION TAGS (100001-100036)
  !---------------------------------------------------------------------------
  ! Crust-mantle displacement/velocity/acceleration
  integer, parameter :: io_tag_ford_undo_d_cm      = 100001
  integer, parameter :: io_tag_ford_undo_v_cm      = 100002
  integer, parameter :: io_tag_ford_undo_a_cm      = 100003

  ! Outer core displacement/velocity/acceleration
  integer, parameter :: io_tag_ford_undo_d_oc      = 100004
  integer, parameter :: io_tag_ford_undo_v_oc      = 100005
  integer, parameter :: io_tag_ford_undo_a_oc      = 100006

  ! Inner core displacement/velocity/acceleration
  integer, parameter :: io_tag_ford_undo_d_ic      = 100007
  integer, parameter :: io_tag_ford_undo_v_ic      = 100008
  integer, parameter :: io_tag_ford_undo_a_ic      = 100009

  ! Crust-mantle strain components
  integer, parameter :: io_tag_ford_undo_eps_xx_cm = 100010
  integer, parameter :: io_tag_ford_undo_eps_yy_cm = 100011
  integer, parameter :: io_tag_ford_undo_eps_xy_cm = 100012
  integer, parameter :: io_tag_ford_undo_eps_xz_cm = 100013
  integer, parameter :: io_tag_ford_undo_eps_yz_cm = 100014

  ! Inner core strain components
  integer, parameter :: io_tag_ford_undo_eps_xx_ic = 100015
  integer, parameter :: io_tag_ford_undo_eps_yy_ic = 100016
  integer, parameter :: io_tag_ford_undo_eps_xy_ic = 100017
  integer, parameter :: io_tag_ford_undo_eps_xz_ic = 100018
  integer, parameter :: io_tag_ford_undo_eps_yz_ic = 100019

  ! Rotation arrays
  integer, parameter :: io_tag_ford_undo_A_rot     = 100020
  integer, parameter :: io_tag_ford_undo_B_rot     = 100021

  ! Crust-mantle attenuation R components
  integer, parameter :: io_tag_ford_undo_R_xx_cm   = 100022
  integer, parameter :: io_tag_ford_undo_R_yy_cm   = 100023
  integer, parameter :: io_tag_ford_undo_R_xy_cm   = 100024
  integer, parameter :: io_tag_ford_undo_R_xz_cm   = 100025
  integer, parameter :: io_tag_ford_undo_R_yz_cm   = 100026

  ! Inner core attenuation R components
  integer, parameter :: io_tag_ford_undo_R_xx_ic   = 100027
  integer, parameter :: io_tag_ford_undo_R_yy_ic   = 100028
  integer, parameter :: io_tag_ford_undo_R_xy_ic   = 100029
  integer, parameter :: io_tag_ford_undo_R_xz_ic   = 100030
  integer, parameter :: io_tag_ford_undo_R_yz_ic   = 100031

  ! Gravity arrays and metadata
  integer, parameter :: io_tag_ford_undo_neq       = 100032
  integer, parameter :: io_tag_ford_undo_neq1      = 100033
  integer, parameter :: io_tag_ford_undo_pgrav1    = 100034
  integer, parameter :: io_tag_nsubset_iterations  = 100035
  integer, parameter :: io_tag_ford_undo_nmsg      = 100036

  !---------------------------------------------------------------------------
  ! SURFACE MOVIE TAGS (110001-110021)
  !---------------------------------------------------------------------------
  ! Metadata tags
  integer :: io_tag_surf_offset         = 110001
  integer :: io_tag_surf_npoints        = 110002

  ! Data tags (displacement components)
  integer :: io_tag_surf_ux             = 110010
  integer :: io_tag_surf_uy             = 110011
  integer :: io_tag_surf_uz             = 110012

  ! Time-stepping tags
  integer, parameter :: io_tag_surf_it_begin       = 110020
  integer, parameter :: io_tag_surf_it_end         = 110021

  !---------------------------------------------------------------------------
  ! VOLUME MOVIE TAGS (120001-120022)
  !---------------------------------------------------------------------------
  ! Metadata tags
  integer, parameter :: io_tag_vol_offset          = 120001
  integer, parameter :: io_tag_vol_npoints         = 120002

  ! Strain component tags (for MOVIE_VOLUME_TYPE 1-3)
  integer, parameter :: io_tag_vol_strain_NN       = 120010
  integer, parameter :: io_tag_vol_strain_EE       = 120011
  integer, parameter :: io_tag_vol_strain_ZZ       = 120012
  integer, parameter :: io_tag_vol_strain_NE       = 120013
  integer, parameter :: io_tag_vol_strain_NZ       = 120014
  integer, parameter :: io_tag_vol_strain_EZ       = 120015

  ! Vector component tags (for MOVIE_VOLUME_TYPE 5-6: displacement/velocity)
  integer, parameter :: io_tag_vol_vec_N           = 120020
  integer, parameter :: io_tag_vol_vec_E           = 120021
  integer, parameter :: io_tag_vol_vec_Z           = 120022

  ! Norm component tags (for MOVIE_VOLUME_TYPE 7-9: displnorm, velnorm, accelnorm)
  ! These send norm data for each region (crust-mantle, outer core, inner core)
  integer, parameter :: io_tag_vol_norm_cm         = 120030
  integer, parameter :: io_tag_vol_norm_oc         = 120031
  integer, parameter :: io_tag_vol_norm_ic         = 120032

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

  ! undo attenuation frame buffers (one snapshot on IO node)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_displ_cm(:,:), undo_frame_veloc_cm(:,:), undo_frame_accel_cm(:,:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_displ_oc(:),   undo_frame_veloc_oc(:),   undo_frame_accel_oc(:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_displ_ic(:,:), undo_frame_veloc_ic(:,:), undo_frame_accel_ic(:,:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_eps_xx_cm(:,:,:,:), undo_frame_eps_yy_cm(:,:,:,:), &
                                                 undo_frame_eps_xy_cm(:,:,:,:), undo_frame_eps_xz_cm(:,:,:,:), &
                                                 undo_frame_eps_yz_cm(:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_eps_xx_ic(:,:,:,:), undo_frame_eps_yy_ic(:,:,:,:), &
                                                 undo_frame_eps_xy_ic(:,:,:,:), undo_frame_eps_xz_ic(:,:,:,:), &
                                                 undo_frame_eps_yz_ic(:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_A_rot(:,:,:,:), undo_frame_B_rot(:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_R_xx_cm(:,:,:,:,:), undo_frame_R_yy_cm(:,:,:,:,:), &
                                                 undo_frame_R_xy_cm(:,:,:,:,:), undo_frame_R_xz_cm(:,:,:,:,:), &
                                                 undo_frame_R_yz_cm(:,:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_R_xx_ic(:,:,:,:,:), undo_frame_R_yy_ic(:,:,:,:,:), &
                                                 undo_frame_R_xy_ic(:,:,:,:,:), undo_frame_R_xz_ic(:,:,:,:,:), &
                                                 undo_frame_R_yz_ic(:,:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable, target :: undo_frame_pgrav1(:)
  integer,          allocatable, target :: undo_frame_neq(:), undo_frame_neq1(:)

  ! encapsulated undo state: flags plus pointer views of all buffers
  type undo_state_type
    logical :: datasets_initialized = .false.
    logical :: file_open            = .false.

    ! pointer aliases to frame buffers (storage remains in the module arrays above)
    real(kind=CUSTOM_REAL), pointer :: displ_cm(:,:)    => null()
    real(kind=CUSTOM_REAL), pointer :: veloc_cm(:,:)    => null()
    real(kind=CUSTOM_REAL), pointer :: accel_cm(:,:)    => null()
    real(kind=CUSTOM_REAL), pointer :: displ_oc(:)      => null()
    real(kind=CUSTOM_REAL), pointer :: veloc_oc(:)      => null()
    real(kind=CUSTOM_REAL), pointer :: accel_oc(:)      => null()
    real(kind=CUSTOM_REAL), pointer :: displ_ic(:,:)    => null()
    real(kind=CUSTOM_REAL), pointer :: veloc_ic(:,:)    => null()
    real(kind=CUSTOM_REAL), pointer :: accel_ic(:,:)    => null()

    real(kind=CUSTOM_REAL), pointer :: eps_xx_cm(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_yy_cm(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_xy_cm(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_xz_cm(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_yz_cm(:,:,:,:) => null()

    real(kind=CUSTOM_REAL), pointer :: eps_xx_ic(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_yy_ic(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_xy_ic(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_xz_ic(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: eps_yz_ic(:,:,:,:) => null()

    real(kind=CUSTOM_REAL), pointer :: A_rot(:,:,:,:)     => null()
    real(kind=CUSTOM_REAL), pointer :: B_rot(:,:,:,:)     => null()

    real(kind=CUSTOM_REAL), pointer :: R_xx_cm(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_yy_cm(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_xy_cm(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_xz_cm(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_yz_cm(:,:,:,:,:) => null()

    real(kind=CUSTOM_REAL), pointer :: R_xx_ic(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_yy_ic(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_xy_ic(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_xz_ic(:,:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: R_yz_ic(:,:,:,:,:) => null()

    real(kind=CUSTOM_REAL), pointer :: pgrav1(:)          => null()
    integer,          pointer :: neq(:)                   => null()
    integer,          pointer :: neq1(:)                  => null()
  end type undo_state_type

  type(undo_state_type), save :: undo_state

  integer, parameter :: UNDO_BUFFER_NONE   = 0
  integer, parameter :: UNDO_BUFFER_REAL1D = 1
  integer, parameter :: UNDO_BUFFER_REAL2D = 2
  integer, parameter :: UNDO_BUFFER_REAL4D = 3
  integer, parameter :: UNDO_BUFFER_REAL5D = 4
  integer, parameter :: UNDO_BUFFER_INT1D  = 5

  integer, parameter :: UNDO_OFFSET_NONE          = 0
  integer, parameter :: UNDO_OFFSET_NGLOB_CM      = 1
  integer, parameter :: UNDO_OFFSET_NGLOB_OC      = 2
  integer, parameter :: UNDO_OFFSET_NGLOB_IC      = 3
  integer, parameter :: UNDO_OFFSET_NSPEC_CM_SOA  = 4
  integer, parameter :: UNDO_OFFSET_NSPEC_IC_SOA  = 5
  integer, parameter :: UNDO_OFFSET_NSPEC_ROT     = 6
  integer, parameter :: UNDO_OFFSET_NSPEC_CM_ATT  = 7
  integer, parameter :: UNDO_OFFSET_NSPEC_IC_ATT  = 8
  integer, parameter :: UNDO_OFFSET_PGRAV1        = 9

  type :: undo_message_descriptor
    integer :: tag = -1
    integer :: buffer_kind = UNDO_BUFFER_NONE
    integer :: offset_kind = UNDO_OFFSET_NONE
    integer :: h5_kind = CUSTOM_REAL
    character(len=64) :: dataset_name = ''
    real(kind=CUSTOM_REAL), pointer :: dest_real1(:) => null()
    real(kind=CUSTOM_REAL), pointer :: dest_real2(:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: dest_real4(:,:,:,:) => null()
    real(kind=CUSTOM_REAL), pointer :: dest_real5(:,:,:,:,:) => null()
    integer, pointer :: dest_int(:) => null()
  end type undo_message_descriptor

  type(undo_message_descriptor), allocatable :: undo_descriptors(:)

  ! verbose output (for debugging)
  logical, parameter :: VERBOSE = .true.


! USE_HDF5
#endif

contains

  !---------------------------------------------------------------------------
  ! HDF5 FILE OPERATION HELPERS
  !---------------------------------------------------------------------------
  ! These helper subroutines consolidate common HDF5 file open/close patterns
  ! that vary based on whether we're in multi-IO-server mode or single mode.
  !---------------------------------------------------------------------------

  subroutine open_hdf5_file_for_io(file_name, use_collective, is_shard_mode)
  ! Open an HDF5 file with appropriate method based on IO mode.
  ! In shard mode (HDF5_IO_NODES > 1), each IO server creates/opens its own file.
  ! In single mode, use parallel/collective or independent access.

#ifdef USE_HDF5

    use manager_hdf5

    implicit none

    character(len=*), intent(in) :: file_name
    logical, intent(in) :: use_collective
    logical, intent(in) :: is_shard_mode

    if (is_shard_mode) then
      ! Multi-IO-server mode: each IO rank creates/opens its own shard file
      call h5_create_or_open_file(file_name)
    else
      ! Single file mode: use parallel access
      if (use_collective) then
        call h5_open_file_p_collect(file_name)
      else
        call h5_open_file_p(file_name)
      endif
    endif

#endif

  end subroutine open_hdf5_file_for_io


  subroutine close_hdf5_file_for_io(is_shard_mode)
  ! Close an HDF5 file with appropriate method based on IO mode.

#ifdef USE_HDF5

    use manager_hdf5

    implicit none

    logical, intent(in) :: is_shard_mode

    if (is_shard_mode) then
      call h5_close_file()
    else
      call h5_close_file_p()
    endif

#endif

  end subroutine close_hdf5_file_for_io


  subroutine ensure_dataset_exists(group_name, dset_name, total_points)
  ! Ensure a dataset exists in the current group, creating it if necessary.
  ! Useful for lazy dataset creation in shard mode.

#ifdef USE_HDF5

    use manager_hdf5
    use constants, only: CUSTOM_REAL, MAX_STRING_LEN

    implicit none

    character(len=*), intent(in) :: group_name
    character(len=*), intent(in) :: dset_name
    integer, intent(in) :: total_points

    character(len=MAX_STRING_LEN) :: dset_full_name
    logical :: dset_exists

    dset_full_name = trim(group_name)//'/'//trim(dset_name)
    call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
    if (.not. dset_exists) then
      call h5_create_dataset_gen_in_group(trim(dset_name), (/total_points/), 1, CUSTOM_REAL)
    endif

#endif

  end subroutine ensure_dataset_exists

  !---------------------------------------------------------------------------
  ! END HDF5 FILE OPERATION HELPERS
  !---------------------------------------------------------------------------

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

  integer function io_local_compute_count()

  use specfem_par, only: NPROCTOT_VAL

  implicit none

  io_local_compute_count = 0

  if (HDF5_IO_NODES <= 1) then
    io_local_compute_count = NPROCTOT_VAL
    return
  endif

  if (.not. allocated(io_nproc_all)) return
  if (my_io_id < 0) return
  if (my_io_id + 1 > size(io_nproc_all)) return

  io_local_compute_count = io_nproc_all(my_io_id+1)

  end function io_local_compute_count

  subroutine build_io_assignment_mask(assigned_to_this_io)

  use specfem_par, only: NPROCTOT_VAL

  implicit none

  logical, dimension(0:NPROCTOT_VAL-1), intent(out) :: assigned_to_this_io
  integer :: i, comp_rank, local_count

  assigned_to_this_io = .false.

  if (HDF5_IO_NODES <= 1) then
    assigned_to_this_io = .true.
    return
  endif

  if (.not. allocated(io_compute_ranks)) return

  local_count = io_local_compute_count()
  if (local_count <= 0) return

  do i = 1, local_count
    comp_rank = io_compute_ranks(my_io_id+1,i)
    if (comp_rank >= 0 .and. comp_rank < NPROCTOT_VAL) then
      assigned_to_this_io(comp_rank) = .true.
    endif
  enddo

  end subroutine build_io_assignment_mask

  subroutine zero_unassigned_offsets(offsets, assigned_to_this_io)

  implicit none

  integer, dimension(0:), intent(inout) :: offsets
  logical, dimension(0:), intent(in) :: assigned_to_this_io
  integer :: rank_id

  do rank_id = 0, ubound(offsets,1)
    if (.not. assigned_to_this_io(rank_id)) offsets(rank_id) = 0
  enddo

  end subroutine zero_unassigned_offsets

  subroutine apply_local_io_partition_to_undo_offsets()

  use specfem_par, only: ATTENUATION_VAL, FULL_GRAVITY_VAL, NPROCTOT_VAL, &
                         ROTATION_VAL, myrank
  use specfem_par_movie_hdf5

  implicit none

  logical, dimension(0:NPROCTOT_VAL-1) :: assigned_to_this_io

  if (HDF5_IO_NODES <= 1) return
  if (.not. IO_storage_task) return

  call build_io_assignment_mask(assigned_to_this_io)

  call zero_unassigned_offsets(offset_nglob_cm, assigned_to_this_io)
  call zero_unassigned_offsets(offset_nglob_oc, assigned_to_this_io)
  call zero_unassigned_offsets(offset_nglob_ic, assigned_to_this_io)
  call zero_unassigned_offsets(offset_nspec_cm_soa, assigned_to_this_io)
  call zero_unassigned_offsets(offset_nspec_ic_soa, assigned_to_this_io)

  if (ROTATION_VAL) then
    call zero_unassigned_offsets(offset_nspec_oc_rot, assigned_to_this_io)
  endif
  if (ATTENUATION_VAL) then
    call zero_unassigned_offsets(offset_nspec_cm_att, assigned_to_this_io)
    call zero_unassigned_offsets(offset_nspec_ic_att, assigned_to_this_io)
  endif
  if (FULL_GRAVITY_VAL) then
    call zero_unassigned_offsets(offset_pgrav1, assigned_to_this_io)
  endif

  npoints_vol_mov_all_proc_cm = sum(offset_nglob_cm)
  npoints_vol_mov_all_proc_oc = sum(offset_nglob_oc)
  npoints_vol_mov_all_proc_ic = sum(offset_nglob_ic)
  nspec_vol_mov_all_proc_cm_soa = sum(offset_nspec_cm_soa)
  nspec_vol_mov_all_proc_ic_soa = sum(offset_nspec_ic_soa)

  if (ROTATION_VAL) then
    nspec_vol_mov_all_proc_oc_rot = sum(offset_nspec_oc_rot)
  endif
  if (ATTENUATION_VAL) then
    nspec_vol_mov_all_proc_cm_att = sum(offset_nspec_cm_att)
    nspec_vol_mov_all_proc_ic_att = sum(offset_nspec_ic_att)
  endif

  if (VERBOSE) then
    print *, 'io_server: rank ', myrank, ' io shard ', my_io_id, &
             ' local undo/nglob totals = ', npoints_vol_mov_all_proc_cm, &
             npoints_vol_mov_all_proc_oc, npoints_vol_mov_all_proc_ic
    call flush_stdout()
  endif

  end subroutine apply_local_io_partition_to_undo_offsets

  subroutine apply_local_io_partition_to_surface_offsets()

  use specfem_par, only: NPROCTOT_VAL, myrank
  use specfem_par_movie_hdf5

  implicit none

  logical, dimension(0:NPROCTOT_VAL-1) :: assigned_to_this_io

  if (HDF5_IO_NODES <= 1) return
  if (.not. IO_storage_task) return
  if (.not. allocated(offset_poin)) return

  call build_io_assignment_mask(assigned_to_this_io)
  call zero_unassigned_offsets(offset_poin, assigned_to_this_io)

  npoints_surf_mov_all_proc = sum(offset_poin)

  if (VERBOSE) then
    print *, 'io_server: rank ', myrank, ' io shard ', my_io_id, &
             ' local surface points = ', npoints_surf_mov_all_proc
    call flush_stdout()
  endif

  end subroutine apply_local_io_partition_to_surface_offsets

  subroutine apply_local_io_partition_to_volume_offsets()

  use specfem_par, only: NPROCTOT_VAL, myrank
  use specfem_par_movie_hdf5

  implicit none

  logical, dimension(0:NPROCTOT_VAL-1) :: assigned_to_this_io

  if (HDF5_IO_NODES <= 1) return
  if (.not. IO_storage_task) return
  if (.not. allocated(offset_poin_vol)) return

  call build_io_assignment_mask(assigned_to_this_io)
  call zero_unassigned_offsets(offset_poin_vol, assigned_to_this_io)

  npoints_vol_mov_all_proc = sum(offset_poin_vol)

  if (VERBOSE) then
    print *, 'io_server: rank ', myrank, ' io shard ', my_io_id, &
             ' local volume points = ', npoints_vol_mov_all_proc
    call flush_stdout()
  endif

  end subroutine apply_local_io_partition_to_volume_offsets

  subroutine pack_local_int_values(global_values, local_values)

  use specfem_par, only: NPROCTOT_VAL

  implicit none

  integer, dimension(0:), intent(in) :: global_values
  integer, dimension(:), intent(out) :: local_values
  integer :: i, comp_rank, local_count

  local_values = 0

  if (HDF5_IO_NODES <= 1) then
    local_count = min(size(local_values), size(global_values))
    if (local_count > 0) local_values(1:local_count) = global_values(0:local_count-1)
    return
  endif

  if (.not. allocated(io_compute_ranks)) return

  local_count = min(size(local_values), io_local_compute_count())
  if (local_count <= 0) return

  do i = 1, local_count
    comp_rank = io_compute_ranks(my_io_id+1,i)
    if (comp_rank >= 0 .and. comp_rank < NPROCTOT_VAL) then
      local_values(i) = global_values(comp_rank)
    endif
  enddo

  end subroutine pack_local_int_values

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine do_io_start_idle()

#ifdef USE_HDF5

  use specfem_par
  use specfem_par_movie_hdf5
  use manager_hdf5
  use constants, only: myrank, my_status_size, my_status_source, my_status_tag, IMAIN

  use io_bandwidth
  use timing_accounting

  implicit none

  integer :: status(my_status_size)
  integer :: tag, tag_src

  ! timing for IO node total elapsed
  double precision :: io_time_start

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

  ! track undo snapshot file state
  logical :: undo_use_collective

  ! array for dumping the data array
  real(kind=CUSTOM_REAL), dimension(:),         allocatable :: dump_ford_undo_1d_glob
  real(kind=CUSTOM_REAL), dimension(:,:),       allocatable :: dump_ford_undo_2d_glob
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),   allocatable :: dump_ford_undo_4d
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: dump_ford_undo_5d

  ! arrays for surface movie values (per rank)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_surf_ux
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_surf_uy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_surf_uz

  ! full-frame buffers for surface movie (per IO node)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_frame_ux
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_frame_uy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_frame_uz

  ! arrays for volume movie values at movie points (per rank)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol1
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol3
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol4
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol5
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol6

  ! full-frame buffers for volume movie (per IO node)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame1
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame3
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame4
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame5
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame6

  ! per-rank receive buffers and full-frame buffers for norm movie types (7,8,9)
  ! These store NGLOB data per region (cm, oc, ic) rather than movie points
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol_norm_cm
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol_norm_oc
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dump_vol_norm_ic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame_norm_cm
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame_norm_oc
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol_frame_norm_ic

  ! maximum nglob and nspec in offset arrays
  integer :: max_nglob
  integer :: max_nspec

  integer :: i_out
  integer :: total_nglob_cm, total_nglob_oc, total_nglob_ic
  integer :: total_nspec_cm_soa, total_nspec_ic_soa
  integer :: total_nspec_oc_rot, total_nspec_cm_att, total_nspec_ic_att
  integer :: total_pgrav1

  ! initialize all counters
  !---------------------------------------------------------------------------
  ! PHASE 1: COUNTER INITIALIZATION
  !---------------------------------------------------------------------------
  call timing_reset()
  io_time_start = MPI_Wtime()

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

  ! undo file state
  undo_state%file_open      = .false.
  undo_use_collective = H5_COL .and. (HDF5_IO_NODES <= 1)

  if (SAVE_FORWARD .and. (MOVIE_SURFACE .or. MOVIE_VOLUME) .and. IO_storage_task) then
    if (my_io_id == 0) then
      write(IMAIN,*) 'Warning: SAVE_FORWARD together with MOVIE output requires buffering undo arrays on IO nodes.'
      write(IMAIN,*) '         Ensure sufficient memory is available for the HDF5 I/O server.'
      call flush_IMAIN()
    endif
  endif


  !---------------------------------------------------------------------------
  ! PHASE 2: BUFFER AND HDF5 INITIALIZATION
  !---------------------------------------------------------------------------

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

    total_nglob_cm = sum(offset_nglob_cm)
    total_nglob_oc = sum(offset_nglob_oc)
    total_nglob_ic = sum(offset_nglob_ic)
    total_nspec_cm_soa = sum(offset_nspec_cm_soa)
    total_nspec_ic_soa = sum(offset_nspec_ic_soa)
    if (ROTATION_VAL) then
      total_nspec_oc_rot = sum(offset_nspec_oc_rot)
    else
      total_nspec_oc_rot = 0
    endif
    if (ATTENUATION_VAL) then
      total_nspec_cm_att = sum(offset_nspec_cm_att)
      total_nspec_ic_att = sum(offset_nspec_ic_att)
    else
      total_nspec_cm_att = 0
      total_nspec_ic_att = 0
    endif
    if (FULL_GRAVITY_VAL) then
      total_pgrav1 = sum(offset_pgrav1)
    else
      total_pgrav1 = 0
    endif

    if (total_nglob_cm > 0) then
      if (.not. allocated(undo_frame_displ_cm)) then
        allocate(undo_frame_displ_cm(NDIM,total_nglob_cm))
        allocate(undo_frame_veloc_cm(NDIM,total_nglob_cm))
        allocate(undo_frame_accel_cm(NDIM,total_nglob_cm))
      endif
      if (.not. associated(undo_state%displ_cm)) then
        undo_state%displ_cm => undo_frame_displ_cm
        undo_state%veloc_cm => undo_frame_veloc_cm
        undo_state%accel_cm => undo_frame_accel_cm
      endif
    endif
    if (total_nglob_oc > 0) then
      if (.not. allocated(undo_frame_displ_oc)) then
        allocate(undo_frame_displ_oc(total_nglob_oc))
        allocate(undo_frame_veloc_oc(total_nglob_oc))
        allocate(undo_frame_accel_oc(total_nglob_oc))
      endif
      if (.not. associated(undo_state%displ_oc)) then
        undo_state%displ_oc => undo_frame_displ_oc
        undo_state%veloc_oc => undo_frame_veloc_oc
        undo_state%accel_oc => undo_frame_accel_oc
      endif
    endif
    if (total_nglob_ic > 0) then
      if (.not. allocated(undo_frame_displ_ic)) then
        allocate(undo_frame_displ_ic(NDIM,total_nglob_ic))
        allocate(undo_frame_veloc_ic(NDIM,total_nglob_ic))
        allocate(undo_frame_accel_ic(NDIM,total_nglob_ic))
      endif
      if (.not. associated(undo_state%displ_ic)) then
        undo_state%displ_ic => undo_frame_displ_ic
        undo_state%veloc_ic => undo_frame_veloc_ic
        undo_state%accel_ic => undo_frame_accel_ic
      endif
    endif
    if (total_nspec_cm_soa > 0) then
      if (.not. allocated(undo_frame_eps_xx_cm)) then
        allocate(undo_frame_eps_xx_cm(NGLLX,NGLLY,NGLLZ,total_nspec_cm_soa))
        allocate(undo_frame_eps_yy_cm(NGLLX,NGLLY,NGLLZ,total_nspec_cm_soa))
        allocate(undo_frame_eps_xy_cm(NGLLX,NGLLY,NGLLZ,total_nspec_cm_soa))
        allocate(undo_frame_eps_xz_cm(NGLLX,NGLLY,NGLLZ,total_nspec_cm_soa))
        allocate(undo_frame_eps_yz_cm(NGLLX,NGLLY,NGLLZ,total_nspec_cm_soa))
      endif
      if (.not. associated(undo_state%eps_xx_cm)) then
        undo_state%eps_xx_cm => undo_frame_eps_xx_cm
        undo_state%eps_yy_cm => undo_frame_eps_yy_cm
        undo_state%eps_xy_cm => undo_frame_eps_xy_cm
        undo_state%eps_xz_cm => undo_frame_eps_xz_cm
        undo_state%eps_yz_cm => undo_frame_eps_yz_cm
      endif
    endif
    if (total_nspec_ic_soa > 0) then
      if (.not. allocated(undo_frame_eps_xx_ic)) then
        allocate(undo_frame_eps_xx_ic(NGLLX,NGLLY,NGLLZ,total_nspec_ic_soa))
        allocate(undo_frame_eps_yy_ic(NGLLX,NGLLY,NGLLZ,total_nspec_ic_soa))
        allocate(undo_frame_eps_xy_ic(NGLLX,NGLLY,NGLLZ,total_nspec_ic_soa))
        allocate(undo_frame_eps_xz_ic(NGLLX,NGLLY,NGLLZ,total_nspec_ic_soa))
        allocate(undo_frame_eps_yz_ic(NGLLX,NGLLY,NGLLZ,total_nspec_ic_soa))
      endif
      if (.not. associated(undo_state%eps_xx_ic)) then
        undo_state%eps_xx_ic => undo_frame_eps_xx_ic
        undo_state%eps_yy_ic => undo_frame_eps_yy_ic
        undo_state%eps_xy_ic => undo_frame_eps_xy_ic
        undo_state%eps_xz_ic => undo_frame_eps_xz_ic
        undo_state%eps_yz_ic => undo_frame_eps_yz_ic
      endif
    endif
    if (ROTATION_VAL .and. total_nspec_oc_rot > 0) then
      if (.not. allocated(undo_frame_A_rot)) then
        allocate(undo_frame_A_rot(NGLLX,NGLLY,NGLLZ,total_nspec_oc_rot))
        allocate(undo_frame_B_rot(NGLLX,NGLLY,NGLLZ,total_nspec_oc_rot))
      endif
      if (.not. associated(undo_state%A_rot)) then
        undo_state%A_rot => undo_frame_A_rot
        undo_state%B_rot => undo_frame_B_rot
      endif
    endif
    if (ATTENUATION_VAL .and. total_nspec_cm_att > 0) then
      if (.not. allocated(undo_frame_R_xx_cm)) then
        allocate(undo_frame_R_xx_cm(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_cm_att))
        allocate(undo_frame_R_yy_cm(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_cm_att))
        allocate(undo_frame_R_xy_cm(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_cm_att))
        allocate(undo_frame_R_xz_cm(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_cm_att))
        allocate(undo_frame_R_yz_cm(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_cm_att))
      endif
      if (.not. associated(undo_state%R_xx_cm)) then
        undo_state%R_xx_cm => undo_frame_R_xx_cm
        undo_state%R_yy_cm => undo_frame_R_yy_cm
        undo_state%R_xy_cm => undo_frame_R_xy_cm
        undo_state%R_xz_cm => undo_frame_R_xz_cm
        undo_state%R_yz_cm => undo_frame_R_yz_cm
      endif
    endif
    if (ATTENUATION_VAL .and. total_nspec_ic_att > 0) then
      if (.not. allocated(undo_frame_R_xx_ic)) then
        allocate(undo_frame_R_xx_ic(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_ic_att))
        allocate(undo_frame_R_yy_ic(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_ic_att))
        allocate(undo_frame_R_xy_ic(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_ic_att))
        allocate(undo_frame_R_xz_ic(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_ic_att))
        allocate(undo_frame_R_yz_ic(NGLLX,NGLLY,NGLLZ,N_SLS,total_nspec_ic_att))
      endif
      if (.not. associated(undo_state%R_xx_ic)) then
        undo_state%R_xx_ic => undo_frame_R_xx_ic
        undo_state%R_yy_ic => undo_frame_R_yy_ic
        undo_state%R_xy_ic => undo_frame_R_xy_ic
        undo_state%R_xz_ic => undo_frame_R_xz_ic
        undo_state%R_yz_ic => undo_frame_R_yz_ic
      endif
    endif
    if (FULL_GRAVITY_VAL) then
      if (.not. allocated(undo_frame_pgrav1) .and. total_pgrav1 > 0) then
        allocate(undo_frame_pgrav1(total_pgrav1))
      endif
      if (total_pgrav1 > 0 .and. .not. associated(undo_state%pgrav1)) then
        undo_state%pgrav1 => undo_frame_pgrav1
      endif
      if (.not. allocated(undo_frame_neq)) then
        allocate(undo_frame_neq(0:NPROCTOT_VAL-1))
      endif
      if (.not. associated(undo_state%neq)) then
        undo_state%neq => undo_frame_neq
      endif
      if (.not. allocated(undo_frame_neq1)) then
        allocate(undo_frame_neq1(0:NPROCTOT_VAL-1))
      endif
      if (.not. associated(undo_state%neq1)) then
        undo_state%neq1 => undo_frame_neq1
      endif
    endif

    call bind_undo_descriptor_targets()

  endif ! UNDO_ATTENUATION .and. SAVE_FORWARD

  ! initialize HDF5 MPI context on IO storage tasks for movie output
  if ((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. HDF5_ENABLED .and. IO_storage_task) then
    call world_get_comm(comm)
    call world_get_info_null(info)
    call h5_initialize()
    call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)
  endif

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
      ! allocate full-frame buffers for all surface movie points handled by this IO node
      if (npoints_surf_mov_all_proc > 0) then
        allocate(surf_frame_ux(npoints_surf_mov_all_proc))
        allocate(surf_frame_uy(npoints_surf_mov_all_proc))
        allocate(surf_frame_uz(npoints_surf_mov_all_proc))
      endif
    endif

  endif

  ! volume movie initialization (only for movie-point-based types: strains/vector)
  if (MOVIE_VOLUME) then

    ! MOVIE_VOLUME_TYPE 1-3: strains / time-integrated / potency (6 components)
    ! MOVIE_VOLUME_TYPE 5-6: displacement / velocity vectors (3 components)
    ! MOVIE_VOLUME_TYPE 7-9: norm movies (1 component per enabled region)
    select case (MOVIE_VOLUME_TYPE)
    case (1,2,3)
      n_msg_vol = 6 * nproc_io
    case (5,6)
      n_msg_vol = 3 * nproc_io
    case (7,8,9)
      ! norm types send one message per enabled region
      n_msg_vol = 0
      if (OUTPUT_CRUST_MANTLE) n_msg_vol = n_msg_vol + nproc_io
      if (OUTPUT_OUTER_CORE) n_msg_vol = n_msg_vol + nproc_io
      if (OUTPUT_INNER_CORE) n_msg_vol = n_msg_vol + nproc_io
    case default
      n_msg_vol = 0
    end select

    ! determine number of volume movie frames from movie parameters
    ! frames are written whenever mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0
    ! and it is within [MOVIE_START,MOVIE_STOP]
    if (NTSTEP_BETWEEN_FRAMES > 0) then
      ! first time step (within this run) that satisfies the movie sampling rule
      it_first_vol = MOVIE_START + &
        ((max(it_begin, MOVIE_START) - MOVIE_START + NTSTEP_BETWEEN_FRAMES - 1) / NTSTEP_BETWEEN_FRAMES) &
        * NTSTEP_BETWEEN_FRAMES

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
      ! allocate full-frame buffers for all volume movie points handled by this IO node
      if (npoints_vol_mov_all_proc > 0) then
        allocate(vol_frame1(npoints_vol_mov_all_proc))
        allocate(vol_frame2(npoints_vol_mov_all_proc))
        allocate(vol_frame3(npoints_vol_mov_all_proc))
        allocate(vol_frame4(npoints_vol_mov_all_proc))
        allocate(vol_frame5(npoints_vol_mov_all_proc))
        allocate(vol_frame6(npoints_vol_mov_all_proc))
      endif
    endif

    ! allocate buffers for norm movie types (7,8,9) which use NGLOB data per region
    if (MOVIE_VOLUME_TYPE >= 7 .and. MOVIE_VOLUME_TYPE <= 9) then
      ! allocate per-rank receive buffers using max_nglob (set during undo attenuation init)
      if (max_nglob > 0) then
        allocate(dump_vol_norm_cm(max_nglob))
        allocate(dump_vol_norm_oc(max_nglob))
        allocate(dump_vol_norm_ic(max_nglob))
      endif
      ! allocate full-frame buffers for all NGLOB points per region
      if (OUTPUT_CRUST_MANTLE .and. npoints_vol_mov_all_proc_cm > 0) then
        allocate(vol_frame_norm_cm(npoints_vol_mov_all_proc_cm))
      endif
      if (OUTPUT_OUTER_CORE .and. npoints_vol_mov_all_proc_oc > 0) then
        allocate(vol_frame_norm_oc(npoints_vol_mov_all_proc_oc))
      endif
      if (OUTPUT_INNER_CORE .and. npoints_vol_mov_all_proc_ic > 0) then
        allocate(vol_frame_norm_ic(npoints_vol_mov_all_proc_ic))
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

  !---------------------------------------------------------------------------
  ! PHASE 3: MAIN IDLING LOOP - MESSAGE PROCESSING
  !---------------------------------------------------------------------------
  ! This loop waits for MPI messages from compute nodes and processes them:
  ! - Undo attenuation snapshots (tags 100001-100034)
  ! - Surface movie frames (tags 110010-110012)
  ! - Volume movie frames (tags 120010-120022)
  !---------------------------------------------------------------------------
  do while ( (UNDO_ATTENUATION .and. SAVE_FORWARD .and. ford_undo_out_count < max_ford_undo_out) .or. &
             (MOVIE_SURFACE .and. HDF5_ENABLED .and. surf_frame_count < max_surf_frames) .or. &
             (MOVIE_VOLUME .and. HDF5_ENABLED .and. vol_frame_count  < max_vol_frames) )

    ! for surface movies without IO server, create the HDF5 group and datasets once per frame
    if (MOVIE_SURFACE .and. HDF5_IO_NODES == 0) then
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
    call timing_idle_start()
    call idle_mpi_io(status)
    call timing_idle_stop()

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

      ! lazily open undo snapshot file on first undo message for this snapshot
      if (.not. undo_state%file_open) then
        ! current snapshot index is ford_undo_out_count (0-based)
        if (HDF5_IO_NODES > 1) then
          file_name = trim(LOCAL_PATH)//'/save_frame_at'//trim(i2c(ford_undo_out_count+1))//'.io'//trim(i2c(my_io_id))//'.h5'
          call h5_create_or_open_file(file_name)
        else
          write(file_name, '(a,i6.6,a)') 'save_frame_at',ford_undo_out_count+1,'.h5'
          file_name = trim(LOCAL_PATH)//'/'//trim(file_name)
          if (undo_use_collective) then
            call h5_open_file_p_collect(file_name)
          else
            call h5_open_file_p(file_name)
          endif
        endif
        if (HDF5_IO_NODES > 1) then
          undo_state%datasets_initialized = .false.
        else
          undo_state%datasets_initialized = .true.
        endif
        undo_state%file_open = .true.
      endif

      ! receive the data
      call timing_io_start()
      call recv_and_write_ford_undo(tag, tag_src, status, &
                                    dump_ford_undo_1d_glob, &
                                    dump_ford_undo_2d_glob, &
                                    dump_ford_undo_4d, &
                                    dump_ford_undo_5d, &
                                    ford_undo_out_count) ! use for the filename
      call timing_io_stop()

      ! count 1 message received
      n_recv_msg_ford_undo = n_recv_msg_ford_undo + 1

    endif ! UNDO_ATTENUATION .and. SAVE_FORWARD

    ! surface movie
    if (MOVIE_SURFACE .and. &
        (tag == io_tag_surf_ux .or. tag == io_tag_surf_uy .or. tag == io_tag_surf_uz)) then

      call timing_io_start()
      call recv_surface_movie(tag, tag_src, status, &
                              dump_surf_ux, dump_surf_uy, dump_surf_uz, &
                              surf_frame_ux, surf_frame_uy, surf_frame_uz)
      call timing_io_stop()

      ! count 1 message received
      n_recv_msg_surf = n_recv_msg_surf + 1

    endif

    ! volume movie (movie-point-based)
    if (MOVIE_VOLUME .and. &
        (tag == io_tag_vol_strain_NN .or. tag == io_tag_vol_strain_EE .or. tag == io_tag_vol_strain_ZZ .or. &
         tag == io_tag_vol_strain_NE .or. tag == io_tag_vol_strain_NZ .or. tag == io_tag_vol_strain_EZ .or. &
         tag == io_tag_vol_vec_N     .or. tag == io_tag_vol_vec_E     .or. tag == io_tag_vol_vec_Z)) then

      call timing_io_start()
      call recv_volume_movie(tag, tag_src, status, &
                 dump_vol1, dump_vol2, dump_vol3, &
                 dump_vol4, dump_vol5, dump_vol6, &
                 vol_frame1, vol_frame2, vol_frame3, &
                 vol_frame4, vol_frame5, vol_frame6)
      call timing_io_stop()

      ! count 1 message received
      n_recv_msg_vol = n_recv_msg_vol + 1

    endif

    ! volume movie (norm types: 7,8,9 - NGLOB data per region)
    if (MOVIE_VOLUME .and. &
        (tag == io_tag_vol_norm_cm .or. tag == io_tag_vol_norm_oc .or. tag == io_tag_vol_norm_ic)) then

      call timing_io_start()
      call recv_volume_norm_movie(tag, tag_src, status, &
                                  dump_vol_norm_cm, dump_vol_norm_oc, dump_vol_norm_ic, &
                                  vol_frame_norm_cm, vol_frame_norm_oc, vol_frame_norm_ic)
      call timing_io_stop()

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

        ! close undo snapshot file for this iteration
        if (undo_state%file_open) then
          call timing_io_start()
          call write_buffered_undo_snapshot(undo_use_collective)
          if (HDF5_IO_NODES > 1) then
            call h5_close_file()
          else
            call h5_close_file_p()
          endif
          call timing_io_stop()
          undo_state%file_open = .false.
          undo_state%datasets_initialized = .false.
        endif

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
          ! all messages for this frame have been received; write buffered frame
          call timing_io_start()
          call write_surface_frame(surf_frame_count, it_first_surf, &
                                   surf_frame_ux, surf_frame_uy, surf_frame_uz)
          call timing_io_stop()
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
          ! all messages for this frame have been received; write buffered frame
          call timing_io_start()
          select case (MOVIE_VOLUME_TYPE)
          case (1,2,3,5,6)
            call write_volume_frame(vol_frame_count, it_first_vol, &
                                    vol_frame1, vol_frame2, vol_frame3, &
                                    vol_frame4, vol_frame5, vol_frame6)
          case (7,8,9)
            call write_volume_norm_frame(vol_frame_count, it_first_vol, &
                                         vol_frame_norm_cm, vol_frame_norm_oc, vol_frame_norm_ic)
          end select
          call timing_io_stop()
          vol_frame_count = vol_frame_count + 1
          n_recv_msg_vol = 0

          ! per-frame IO bandwidth tracking for volume movie
          call calculate_bandwidth_all_procs()
          call initialize_bytes_written()

          if (VERBOSE .and. IO_storage_task) then
            print *, 'io_server: completed volume frame ', vol_frame_count, '/', max_vol_frames
          endif
        endif
      endif
    endif

  enddo
  !---------------------------------------------------------------------------
  ! END OF MAIN IDLING LOOP
  !---------------------------------------------------------------------------

  ! timing accounting report for IO node (no barriers)
  call timing_report(myrank, mygroup, .true., MPI_Wtime() - io_time_start)

  call calculate_bandwidth_all_procs()

  if (VERBOSE .and. IO_storage_task) then
    print *, 'io_server: leaving do_io_start_idle'
    print *, '  final surf frames: ', surf_frame_count, '/', max_surf_frames
    print *, '  final vol  frames: ', vol_frame_count,  '/', max_vol_frames
  endif

  !---------------------------------------------------------------------------
  ! PHASE 4: CLEANUP - DEALLOCATE TEMPORARY BUFFERS
  !---------------------------------------------------------------------------
  ! undo attenuation
  if (UNDO_ATTENUATION .and. SAVE_FORWARD) then
    deallocate(dump_ford_undo_1d_glob, &
               dump_ford_undo_2d_glob, &
               dump_ford_undo_4d, &
               dump_ford_undo_5d)
    if (allocated(undo_descriptors)) then
      deallocate(undo_descriptors)
    endif
    if (allocated(undo_tag_list)) then
      deallocate(undo_tag_list)
    endif
    n_undo_tags = 0
  endif

  ! surface movie buffers
  if (MOVIE_SURFACE) then
    if (allocated(dump_surf_ux)) deallocate(dump_surf_ux)
    if (allocated(dump_surf_uy)) deallocate(dump_surf_uy)
    if (allocated(dump_surf_uz)) deallocate(dump_surf_uz)
    if (allocated(surf_frame_ux)) deallocate(surf_frame_ux)
    if (allocated(surf_frame_uy)) deallocate(surf_frame_uy)
    if (allocated(surf_frame_uz)) deallocate(surf_frame_uz)
  endif

  ! volume movie buffers
  if (MOVIE_VOLUME) then
    if (allocated(dump_vol1)) deallocate(dump_vol1)
    if (allocated(dump_vol2)) deallocate(dump_vol2)
    if (allocated(dump_vol3)) deallocate(dump_vol3)
    if (allocated(dump_vol4)) deallocate(dump_vol4)
    if (allocated(dump_vol5)) deallocate(dump_vol5)
    if (allocated(dump_vol6)) deallocate(dump_vol6)
    if (allocated(vol_frame1)) deallocate(vol_frame1)
    if (allocated(vol_frame2)) deallocate(vol_frame2)
    if (allocated(vol_frame3)) deallocate(vol_frame3)
    if (allocated(vol_frame4)) deallocate(vol_frame4)
    if (allocated(vol_frame5)) deallocate(vol_frame5)
    if (allocated(vol_frame6)) deallocate(vol_frame6)
    ! norm movie buffers (types 7,8,9)
    if (allocated(dump_vol_norm_cm)) deallocate(dump_vol_norm_cm)
    if (allocated(dump_vol_norm_oc)) deallocate(dump_vol_norm_oc)
    if (allocated(dump_vol_norm_ic)) deallocate(dump_vol_norm_ic)
    if (allocated(vol_frame_norm_cm)) deallocate(vol_frame_norm_cm)
    if (allocated(vol_frame_norm_oc)) deallocate(vol_frame_norm_oc)
    if (allocated(vol_frame_norm_ic)) deallocate(vol_frame_norm_ic)
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
  !---------------------------------------------------------------------------
  ! Receive metadata from compute rank 0 on IO server tasks.
  ! This is the RECEIVER side of the metadata exchange.
  ! MIRROR FUNCTION of pass_info_to_io() - must be kept in sync.
  !
  ! Data exchanged:
  ! - UNDO_ATTENUATION: offset arrays, NSUBSET_ITERATIONS, n_msg_ford_undo
  ! - MOVIE_SURFACE: offset_poin, npoints_surf_mov_all_proc, time range
  ! - MOVIE_VOLUME: offset_poin_vol, npoints_vol_mov_all_proc
  !---------------------------------------------------------------------------

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

    ! In shard mode, each IO rank must reduce undo offsets to its assigned
    ! compute ranks before buffer sizing and dataset creation.
    call apply_local_io_partition_to_undo_offsets()

  else
    n_msg_ford_undo = 0
    NSUBSET_ITERATIONS = 0
  endif ! UNDO_ATTENUATION and SAVE_FORWARD

  if (MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE >= 7 .and. MOVIE_VOLUME_TYPE <= 9) then
    if (.not. (UNDO_ATTENUATION .and. SAVE_FORWARD)) then
      call recv_i_inter(offset_nglob_cm, NPROCTOT_VAL, 0, io_tag_ford_undo_d_cm)
      call recv_i_inter(offset_nglob_oc, NPROCTOT_VAL, 0, io_tag_ford_undo_d_oc)
      call recv_i_inter(offset_nglob_ic, NPROCTOT_VAL, 0, io_tag_ford_undo_d_ic)
      call apply_local_io_partition_to_undo_offsets()
    endif
  endif

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

    call apply_local_io_partition_to_surface_offsets()

  endif

  ! receive time-stepping range for movie frames (used by both surface and volume movies)
  if (MOVIE_SURFACE .or. MOVIE_VOLUME) then
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

    call apply_local_io_partition_to_volume_offsets()

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
  !---------------------------------------------------------------------------
  ! Send metadata from compute rank 0 to all IO server tasks.
  ! This is the SENDER side of the metadata exchange.
  ! MIRROR FUNCTION of get_info_from_comp() - must be kept in sync.
  !
  ! Data exchanged:
  ! - UNDO_ATTENUATION: offset arrays, NSUBSET_ITERATIONS, n_msg_ford_undo
  ! - MOVIE_SURFACE: offset_poin, npoints_surf_mov_all_proc, time range
  ! - MOVIE_VOLUME: offset_poin_vol, npoints_vol_mov_all_proc
  !---------------------------------------------------------------------------

#ifdef USE_HDF5

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer :: i_ionod
    integer :: tmp_int

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

    if (MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE >= 7 .and. MOVIE_VOLUME_TYPE <= 9) then
      if (myrank == 0 .and. .not. (UNDO_ATTENUATION .and. SAVE_FORWARD)) then
        do i_ionod = 0, HDF5_IO_NODES-1
          call send_i_inter(offset_nglob_cm, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_cm)
          call send_i_inter(offset_nglob_oc, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_oc)
          call send_i_inter(offset_nglob_ic, NPROCTOT_VAL, i_ionod, io_tag_ford_undo_d_ic)
        enddo
      endif
    endif

    ! surface movie metadata
    if (MOVIE_SURFACE) then
      if (myrank == 0) then

        do i_ionod = 0, HDF5_IO_NODES-1

          ! send per-rank surface offsets and total number of points
          call send_i_inter(offset_poin, NPROCTOT_VAL, i_ionod, io_tag_surf_offset)
          tmp_int = npoints_surf_mov_all_proc
          call send_i_inter((/tmp_int/), 1, i_ionod, io_tag_surf_npoints)

        enddo
      endif
    endif

    ! send time-stepping range for movie frames (used by both surface and volume movies)
    if (MOVIE_SURFACE .or. MOVIE_VOLUME) then
      if (myrank == 0) then

        do i_ionod = 0, HDF5_IO_NODES-1

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

  integer :: idx

  if (allocated(undo_tag_list)) deallocate(undo_tag_list)
  if (allocated(undo_descriptors)) deallocate(undo_descriptors)

  n_undo_tags = n_msg_ford_undo

  if (n_undo_tags <= 0) return

  allocate(undo_tag_list(n_undo_tags))
  allocate(undo_descriptors(n_undo_tags))

  undo_tag_list = 0
  undo_descriptors%tag = -1

  idx = 0

  call register_undo_descriptor(idx, io_tag_ford_undo_d_cm, 'displ_crust_mantle', &
                                UNDO_BUFFER_REAL2D, UNDO_OFFSET_NGLOB_CM)
  call register_undo_descriptor(idx, io_tag_ford_undo_v_cm, 'veloc_crust_mantle', &
                                UNDO_BUFFER_REAL2D, UNDO_OFFSET_NGLOB_CM)
  call register_undo_descriptor(idx, io_tag_ford_undo_a_cm, 'accel_crust_mantle', &
                                UNDO_BUFFER_REAL2D, UNDO_OFFSET_NGLOB_CM)
  call register_undo_descriptor(idx, io_tag_ford_undo_d_oc, 'displ_outer_core', &
                                UNDO_BUFFER_REAL1D, UNDO_OFFSET_NGLOB_OC)
  call register_undo_descriptor(idx, io_tag_ford_undo_v_oc, 'veloc_outer_core', &
                                UNDO_BUFFER_REAL1D, UNDO_OFFSET_NGLOB_OC)
  call register_undo_descriptor(idx, io_tag_ford_undo_a_oc, 'accel_outer_core', &
                                UNDO_BUFFER_REAL1D, UNDO_OFFSET_NGLOB_OC)
  call register_undo_descriptor(idx, io_tag_ford_undo_d_ic, 'displ_inner_core', &
                                UNDO_BUFFER_REAL2D, UNDO_OFFSET_NGLOB_IC)
  call register_undo_descriptor(idx, io_tag_ford_undo_v_ic, 'veloc_inner_core', &
                                UNDO_BUFFER_REAL2D, UNDO_OFFSET_NGLOB_IC)
  call register_undo_descriptor(idx, io_tag_ford_undo_a_ic, 'accel_inner_core', &
                                UNDO_BUFFER_REAL2D, UNDO_OFFSET_NGLOB_IC)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_xx_cm, 'epsilondev_xx_crust_mantle', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_CM_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_yy_cm, 'epsilondev_yy_crust_mantle', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_CM_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_xy_cm, 'epsilondev_xy_crust_mantle', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_CM_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_xz_cm, 'epsilondev_xz_crust_mantle', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_CM_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_yz_cm, 'epsilondev_yz_crust_mantle', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_CM_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_xx_ic, 'epsilondev_xx_inner_core', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_IC_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_yy_ic, 'epsilondev_yy_inner_core', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_IC_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_xy_ic, 'epsilondev_xy_inner_core', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_IC_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_xz_ic, 'epsilondev_xz_inner_core', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_IC_SOA)
  call register_undo_descriptor(idx, io_tag_ford_undo_eps_yz_ic, 'epsilondev_yz_inner_core', &
                                UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_IC_SOA)

  if (ROTATION_VAL) then
    call register_undo_descriptor(idx, io_tag_ford_undo_A_rot, 'A_array_rotation', &
                                  UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_ROT)
    call register_undo_descriptor(idx, io_tag_ford_undo_B_rot, 'B_array_rotation', &
                                  UNDO_BUFFER_REAL4D, UNDO_OFFSET_NSPEC_ROT)
  endif

  if (ATTENUATION_VAL) then
    call register_undo_descriptor(idx, io_tag_ford_undo_R_xx_cm, 'R_xx_crust_mantle', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_CM_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_yy_cm, 'R_yy_crust_mantle', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_CM_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_xy_cm, 'R_xy_crust_mantle', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_CM_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_xz_cm, 'R_xz_crust_mantle', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_CM_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_yz_cm, 'R_yz_crust_mantle', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_CM_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_xx_ic, 'R_xx_inner_core', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_IC_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_yy_ic, 'R_yy_inner_core', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_IC_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_xy_ic, 'R_xy_inner_core', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_IC_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_xz_ic, 'R_xz_inner_core', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_IC_ATT)
    call register_undo_descriptor(idx, io_tag_ford_undo_R_yz_ic, 'R_yz_inner_core', &
                                  UNDO_BUFFER_REAL5D, UNDO_OFFSET_NSPEC_IC_ATT)
  endif

  if (FULL_GRAVITY_VAL) then
    call register_undo_descriptor(idx, io_tag_ford_undo_neq, 'neq', &
                                  UNDO_BUFFER_INT1D, UNDO_OFFSET_NONE, data_kind=1)
    call register_undo_descriptor(idx, io_tag_ford_undo_neq1, 'neq1', &
                                  UNDO_BUFFER_INT1D, UNDO_OFFSET_NONE, data_kind=1)
    call register_undo_descriptor(idx, io_tag_ford_undo_pgrav1, 'pgrav1', &
                                  UNDO_BUFFER_REAL1D, UNDO_OFFSET_PGRAV1)
  endif

  if (idx /= n_undo_tags) then
    print *, 'build_undo_tag_list: mismatch between constructed tag count and n_msg_ford_undo', idx, n_undo_tags
  endif

#else

  ! no-op when built without HDF5 support

#endif

  end subroutine build_undo_tag_list

  subroutine register_undo_descriptor(idx, tag, dataset_name, buffer_kind, offset_kind, data_kind)

#ifdef USE_HDF5

  implicit none

  integer, intent(inout) :: idx
  integer, intent(in) :: tag, buffer_kind, offset_kind
  integer, intent(in), optional :: data_kind
  character(len=*), intent(in) :: dataset_name

  idx = idx + 1
  if (.not. allocated(undo_tag_list)) return
  if (idx > size(undo_tag_list)) return

  undo_tag_list(idx) = tag
  undo_descriptors(idx)%tag = tag
  undo_descriptors(idx)%buffer_kind = buffer_kind
  undo_descriptors(idx)%offset_kind = offset_kind
  undo_descriptors(idx)%dataset_name = trim(dataset_name)
  if (present(data_kind)) then
    undo_descriptors(idx)%h5_kind = data_kind
  else
    undo_descriptors(idx)%h5_kind = CUSTOM_REAL
  endif

#else

  ! no-op

#endif

  end subroutine register_undo_descriptor

  subroutine bind_undo_descriptor_targets()
  ! Bind undo_state pointers to descriptor destination pointers.
  ! Uses a helper subroutine to avoid 30+ case statements.

#ifdef USE_HDF5

  implicit none

  integer :: i

  if (.not. allocated(undo_descriptors)) return

  do i = 1, size(undo_descriptors)
    call bind_single_undo_descriptor(undo_descriptors(i))
  enddo

#else

  ! no-op

#endif

  end subroutine bind_undo_descriptor_targets


  subroutine bind_single_undo_descriptor(desc)
  ! Helper subroutine to bind a single descriptor's destination pointer
  ! based on its MPI tag. This centralizes the tag-to-pointer mapping.

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(inout) :: desc

  ! Map tag to appropriate undo_state pointer
  ! Crust-mantle displacement/velocity/acceleration (2D real arrays)
  if (desc%tag == io_tag_ford_undo_d_cm) then
    if (associated(undo_state%displ_cm)) desc%dest_real2 => undo_state%displ_cm
  else if (desc%tag == io_tag_ford_undo_v_cm) then
    if (associated(undo_state%veloc_cm)) desc%dest_real2 => undo_state%veloc_cm
  else if (desc%tag == io_tag_ford_undo_a_cm) then
    if (associated(undo_state%accel_cm)) desc%dest_real2 => undo_state%accel_cm

  ! Outer core displacement/velocity/acceleration (1D real arrays)
  else if (desc%tag == io_tag_ford_undo_d_oc) then
    if (associated(undo_state%displ_oc)) desc%dest_real1 => undo_state%displ_oc
  else if (desc%tag == io_tag_ford_undo_v_oc) then
    if (associated(undo_state%veloc_oc)) desc%dest_real1 => undo_state%veloc_oc
  else if (desc%tag == io_tag_ford_undo_a_oc) then
    if (associated(undo_state%accel_oc)) desc%dest_real1 => undo_state%accel_oc

  ! Inner core displacement/velocity/acceleration (2D real arrays)
  else if (desc%tag == io_tag_ford_undo_d_ic) then
    if (associated(undo_state%displ_ic)) desc%dest_real2 => undo_state%displ_ic
  else if (desc%tag == io_tag_ford_undo_v_ic) then
    if (associated(undo_state%veloc_ic)) desc%dest_real2 => undo_state%veloc_ic
  else if (desc%tag == io_tag_ford_undo_a_ic) then
    if (associated(undo_state%accel_ic)) desc%dest_real2 => undo_state%accel_ic

  ! Crust-mantle strain components (4D real arrays)
  else if (desc%tag == io_tag_ford_undo_eps_xx_cm) then
    if (associated(undo_state%eps_xx_cm)) desc%dest_real4 => undo_state%eps_xx_cm
  else if (desc%tag == io_tag_ford_undo_eps_yy_cm) then
    if (associated(undo_state%eps_yy_cm)) desc%dest_real4 => undo_state%eps_yy_cm
  else if (desc%tag == io_tag_ford_undo_eps_xy_cm) then
    if (associated(undo_state%eps_xy_cm)) desc%dest_real4 => undo_state%eps_xy_cm
  else if (desc%tag == io_tag_ford_undo_eps_xz_cm) then
    if (associated(undo_state%eps_xz_cm)) desc%dest_real4 => undo_state%eps_xz_cm
  else if (desc%tag == io_tag_ford_undo_eps_yz_cm) then
    if (associated(undo_state%eps_yz_cm)) desc%dest_real4 => undo_state%eps_yz_cm

  ! Inner core strain components (4D real arrays)
  else if (desc%tag == io_tag_ford_undo_eps_xx_ic) then
    if (associated(undo_state%eps_xx_ic)) desc%dest_real4 => undo_state%eps_xx_ic
  else if (desc%tag == io_tag_ford_undo_eps_yy_ic) then
    if (associated(undo_state%eps_yy_ic)) desc%dest_real4 => undo_state%eps_yy_ic
  else if (desc%tag == io_tag_ford_undo_eps_xy_ic) then
    if (associated(undo_state%eps_xy_ic)) desc%dest_real4 => undo_state%eps_xy_ic
  else if (desc%tag == io_tag_ford_undo_eps_xz_ic) then
    if (associated(undo_state%eps_xz_ic)) desc%dest_real4 => undo_state%eps_xz_ic
  else if (desc%tag == io_tag_ford_undo_eps_yz_ic) then
    if (associated(undo_state%eps_yz_ic)) desc%dest_real4 => undo_state%eps_yz_ic

  ! Rotation arrays (4D real arrays)
  else if (desc%tag == io_tag_ford_undo_A_rot) then
    if (associated(undo_state%A_rot)) desc%dest_real4 => undo_state%A_rot
  else if (desc%tag == io_tag_ford_undo_B_rot) then
    if (associated(undo_state%B_rot)) desc%dest_real4 => undo_state%B_rot

  ! Crust-mantle attenuation R components (5D real arrays)
  else if (desc%tag == io_tag_ford_undo_R_xx_cm) then
    if (associated(undo_state%R_xx_cm)) desc%dest_real5 => undo_state%R_xx_cm
  else if (desc%tag == io_tag_ford_undo_R_yy_cm) then
    if (associated(undo_state%R_yy_cm)) desc%dest_real5 => undo_state%R_yy_cm
  else if (desc%tag == io_tag_ford_undo_R_xy_cm) then
    if (associated(undo_state%R_xy_cm)) desc%dest_real5 => undo_state%R_xy_cm
  else if (desc%tag == io_tag_ford_undo_R_xz_cm) then
    if (associated(undo_state%R_xz_cm)) desc%dest_real5 => undo_state%R_xz_cm
  else if (desc%tag == io_tag_ford_undo_R_yz_cm) then
    if (associated(undo_state%R_yz_cm)) desc%dest_real5 => undo_state%R_yz_cm

  ! Inner core attenuation R components (5D real arrays)
  else if (desc%tag == io_tag_ford_undo_R_xx_ic) then
    if (associated(undo_state%R_xx_ic)) desc%dest_real5 => undo_state%R_xx_ic
  else if (desc%tag == io_tag_ford_undo_R_yy_ic) then
    if (associated(undo_state%R_yy_ic)) desc%dest_real5 => undo_state%R_yy_ic
  else if (desc%tag == io_tag_ford_undo_R_xy_ic) then
    if (associated(undo_state%R_xy_ic)) desc%dest_real5 => undo_state%R_xy_ic
  else if (desc%tag == io_tag_ford_undo_R_xz_ic) then
    if (associated(undo_state%R_xz_ic)) desc%dest_real5 => undo_state%R_xz_ic
  else if (desc%tag == io_tag_ford_undo_R_yz_ic) then
    if (associated(undo_state%R_yz_ic)) desc%dest_real5 => undo_state%R_yz_ic

  ! Gravity arrays (integer and real arrays)
  else if (desc%tag == io_tag_ford_undo_neq) then
    if (associated(undo_state%neq)) desc%dest_int => undo_state%neq
  else if (desc%tag == io_tag_ford_undo_neq1) then
    if (associated(undo_state%neq1)) desc%dest_int => undo_state%neq1
  else if (desc%tag == io_tag_ford_undo_pgrav1) then
    if (associated(undo_state%pgrav1)) desc%dest_real1 => undo_state%pgrav1
  endif
  ! Unknown tags are silently ignored (no action needed)

#else

  ! no-op

#endif

  end subroutine bind_single_undo_descriptor

  integer function find_undo_descriptor_index(tag)

#ifdef USE_HDF5

  implicit none

  integer, intent(in) :: tag
  integer :: i

  find_undo_descriptor_index = 0
  if (.not. allocated(undo_descriptors)) return

  do i = 1, size(undo_descriptors)
    if (undo_descriptors(i)%tag == tag) then
      find_undo_descriptor_index = i
      return
    endif
  enddo

#else

  find_undo_descriptor_index = 0

#endif

  end function find_undo_descriptor_index

  subroutine compute_hyperslab_bounds(offset_kind, tag_src, ista, data_len)

#ifdef USE_HDF5

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer, intent(in) :: offset_kind, tag_src
  integer, intent(out) :: ista, data_len

  ista = 0
  data_len = 0

  select case (offset_kind)
  case (UNDO_OFFSET_NGLOB_CM)
    call compute_prefix_sum(offset_nglob_cm, tag_src, ista, data_len)
  case (UNDO_OFFSET_NGLOB_OC)
    call compute_prefix_sum(offset_nglob_oc, tag_src, ista, data_len)
  case (UNDO_OFFSET_NGLOB_IC)
    call compute_prefix_sum(offset_nglob_ic, tag_src, ista, data_len)
  case (UNDO_OFFSET_NSPEC_CM_SOA)
    call compute_prefix_sum(offset_nspec_cm_soa, tag_src, ista, data_len)
  case (UNDO_OFFSET_NSPEC_IC_SOA)
    call compute_prefix_sum(offset_nspec_ic_soa, tag_src, ista, data_len)
  case (UNDO_OFFSET_NSPEC_ROT)
    call compute_prefix_sum(offset_nspec_oc_rot, tag_src, ista, data_len)
  case (UNDO_OFFSET_NSPEC_CM_ATT)
    call compute_prefix_sum(offset_nspec_cm_att, tag_src, ista, data_len)
  case (UNDO_OFFSET_NSPEC_IC_ATT)
    call compute_prefix_sum(offset_nspec_ic_att, tag_src, ista, data_len)
  case (UNDO_OFFSET_PGRAV1)
    call compute_prefix_sum(offset_pgrav1, tag_src, ista, data_len)
  case default
    ista = 0
    data_len = 0
  end select

#else

  ista = 0
  data_len = 0

#endif

  end subroutine compute_hyperslab_bounds

  subroutine compute_prefix_sum(array, tag_src, ista, data_len)

#ifdef USE_HDF5

  implicit none

  integer, dimension(:), intent(in) :: array
  integer, intent(in) :: tag_src
  integer, intent(out) :: ista, data_len

  if (tag_src > 0) then
    ista = sum(array(0:tag_src-1))
  else
    ista = 0
  endif
  data_len = array(tag_src)

#else

  ista = 0
  data_len = 0

#endif

  end subroutine compute_prefix_sum

  subroutine handle_real1d_descriptor(desc, tag_src, tag, msg_size, scratch, element_bytes)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  integer, intent(in) :: tag_src, tag, msg_size
  real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: scratch
  integer, intent(out) :: element_bytes
  integer :: ista, data_len, iend, req_dummy

  call compute_hyperslab_bounds(desc%offset_kind, tag_src, ista, data_len)
  iend = ista + data_len

  call irecvv_cr_inter(scratch(1:data_len), msg_size, tag_src, tag, req_dummy)
  if (data_len > 0 .and. associated(desc%dest_real1)) then
    desc%dest_real1(ista+1:iend) = scratch(1:data_len)
  endif

  element_bytes = CUSTOM_REAL

#else

  element_bytes = 0

#endif

  end subroutine handle_real1d_descriptor

  subroutine handle_real2d_descriptor(desc, tag_src, tag, msg_size, scratch, element_bytes)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  integer, intent(in) :: tag_src, tag, msg_size
  real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: scratch
  integer, intent(out) :: element_bytes
  integer :: ista, data_len, iend, req_dummy

  call compute_hyperslab_bounds(desc%offset_kind, tag_src, ista, data_len)
  iend = ista + data_len

  call irecvv_cr_inter(scratch(:,1:data_len), msg_size, tag_src, tag, req_dummy)
  if (data_len > 0 .and. associated(desc%dest_real2)) then
    desc%dest_real2(:,ista+1:iend) = scratch(:,1:data_len)
  endif

  element_bytes = CUSTOM_REAL

#else

  element_bytes = 0

#endif

  end subroutine handle_real2d_descriptor

  subroutine handle_real4d_descriptor(desc, tag_src, tag, msg_size, scratch, element_bytes)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  integer, intent(in) :: tag_src, tag, msg_size
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(inout) :: scratch
  integer, intent(out) :: element_bytes
  integer :: ista, data_len, iend, req_dummy

  call compute_hyperslab_bounds(desc%offset_kind, tag_src, ista, data_len)
  iend = ista + data_len

  call irecvv_cr_inter(scratch(:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
  if (data_len > 0 .and. associated(desc%dest_real4)) then
    desc%dest_real4(:,:,:,ista+1:iend) = scratch(:,:,:,1:data_len)
  endif

  element_bytes = CUSTOM_REAL

#else

  element_bytes = 0

#endif

  end subroutine handle_real4d_descriptor

  subroutine handle_real5d_descriptor(desc, tag_src, tag, msg_size, scratch, element_bytes)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  integer, intent(in) :: tag_src, tag, msg_size
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(inout) :: scratch
  integer, intent(out) :: element_bytes
  integer :: ista, data_len, iend, req_dummy

  call compute_hyperslab_bounds(desc%offset_kind, tag_src, ista, data_len)
  iend = ista + data_len

  call irecvv_cr_inter(scratch(:,:,:,:,1:data_len), msg_size, tag_src, tag, req_dummy)
  if (data_len > 0 .and. associated(desc%dest_real5)) then
    desc%dest_real5(:,:,:,:,ista+1:iend) = scratch(:,:,:,:,1:data_len)
  endif

  element_bytes = CUSTOM_REAL

#else

  element_bytes = 0

#endif

  end subroutine handle_real5d_descriptor

  subroutine handle_int_descriptor(desc, tag_src, tag, element_bytes)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  integer, intent(in) :: tag_src, tag
  integer, intent(out) :: element_bytes
  integer :: req_dummy
  integer :: val

  call irecv_i_inter((/val/), 1, tag_src, tag, req_dummy)
  if (associated(desc%dest_int)) then
    desc%dest_int(tag_src) = val
  endif

  element_bytes = 4

#else

  element_bytes = 0

#endif

  end subroutine handle_int_descriptor

  logical function descriptor_has_data(desc)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc

  select case (desc%buffer_kind)
  case (UNDO_BUFFER_REAL1D)
    descriptor_has_data = associated(desc%dest_real1)
  case (UNDO_BUFFER_REAL2D)
    descriptor_has_data = associated(desc%dest_real2)
  case (UNDO_BUFFER_REAL4D)
    descriptor_has_data = associated(desc%dest_real4)
  case (UNDO_BUFFER_REAL5D)
    descriptor_has_data = associated(desc%dest_real5)
  case (UNDO_BUFFER_INT1D)
    descriptor_has_data = associated(desc%dest_int)
  case default
    descriptor_has_data = .false.
  end select

#else

  descriptor_has_data = .false.

#endif

  end function descriptor_has_data

  subroutine descriptor_shape(desc, dims, rank)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  integer, dimension(5), intent(out) :: dims
  integer, intent(out) :: rank

  dims = 0
  rank = 0

  select case (desc%buffer_kind)
  case (UNDO_BUFFER_REAL1D)
    if (associated(desc%dest_real1)) then
      rank = 1
      dims(1) = size(desc%dest_real1)
    endif
  case (UNDO_BUFFER_REAL2D)
    if (associated(desc%dest_real2)) then
      rank = 2
      dims(1) = size(desc%dest_real2,1)
      dims(2) = size(desc%dest_real2,2)
    endif
  case (UNDO_BUFFER_REAL4D)
    if (associated(desc%dest_real4)) then
      rank = 4
      dims(1) = size(desc%dest_real4,1)
      dims(2) = size(desc%dest_real4,2)
      dims(3) = size(desc%dest_real4,3)
      dims(4) = size(desc%dest_real4,4)
    endif
  case (UNDO_BUFFER_REAL5D)
    if (associated(desc%dest_real5)) then
      rank = 5
      dims(1) = size(desc%dest_real5,1)
      dims(2) = size(desc%dest_real5,2)
      dims(3) = size(desc%dest_real5,3)
      dims(4) = size(desc%dest_real5,4)
      dims(5) = size(desc%dest_real5,5)
    endif
  case (UNDO_BUFFER_INT1D)
    if (associated(desc%dest_int)) then
      rank = 1
      if (HDF5_IO_NODES > 1 .and. IO_storage_task) then
        dims(1) = io_local_compute_count()
      else
        dims(1) = size(desc%dest_int)
      endif
      if (dims(1) <= 0) rank = 0
    endif
  case default
    rank = 0
  end select

#else

  dims = 0
  rank = 0

#endif

  end subroutine descriptor_shape

  integer function descriptor_element_size(desc)

#ifdef USE_HDF5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc

  if (desc%buffer_kind == UNDO_BUFFER_INT1D) then
    descriptor_element_size = 4
  else
    descriptor_element_size = CUSTOM_REAL
  endif

#else

  descriptor_element_size = 0

#endif

  end function descriptor_element_size

  subroutine ensure_descriptor_dataset(desc)

#ifdef USE_HDF5

  use manager_hdf5

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  logical :: dset_exists
  integer :: dims(5)
  integer :: rank

  if (.not. descriptor_has_data(desc)) return

  call descriptor_shape(desc, dims, rank)
  if (rank == 0) return

  call h5_check_dataset_exists(trim(desc%dataset_name), dset_exists)
  if (.not. dset_exists) then
    call h5_create_dataset_gen(trim(desc%dataset_name), dims(1:rank), rank, desc%h5_kind)
  endif

#else

  ! no-op

#endif

  end subroutine ensure_descriptor_dataset

  subroutine write_descriptor_dataset(desc, use_collective)

#ifdef USE_HDF5

  use manager_hdf5
  use io_bandwidth

  implicit none

  type(undo_message_descriptor), intent(in) :: desc
  logical, intent(in) :: use_collective
  integer :: dims(5)
  integer :: rank
  integer :: start_idx(5)
  integer :: num_elements
  integer :: elem_size
  integer, allocatable :: packed_int(:)

  if (.not. descriptor_has_data(desc)) return

  call descriptor_shape(desc, dims, rank)
  if (rank == 0) return

  start_idx(1:rank) = 0

  call start_timer()
  select case (desc%buffer_kind)
  case (UNDO_BUFFER_REAL1D)
    call h5_write_dataset_collect_hyperslab(trim(desc%dataset_name), desc%dest_real1, start_idx(1:rank), use_collective)
  case (UNDO_BUFFER_REAL2D)
    call h5_write_dataset_collect_hyperslab(trim(desc%dataset_name), desc%dest_real2, start_idx(1:rank), use_collective)
  case (UNDO_BUFFER_REAL4D)
    call h5_write_dataset_collect_hyperslab(trim(desc%dataset_name), desc%dest_real4, start_idx(1:rank), use_collective)
  case (UNDO_BUFFER_REAL5D)
    call h5_write_dataset_collect_hyperslab(trim(desc%dataset_name), desc%dest_real5, start_idx(1:rank), use_collective)
  case (UNDO_BUFFER_INT1D)
    if (HDF5_IO_NODES > 1 .and. IO_storage_task) then
      allocate(packed_int(dims(1)))
      call pack_local_int_values(desc%dest_int, packed_int)
      call h5_write_dataset_collect_hyperslab(trim(desc%dataset_name), packed_int, start_idx(1:rank), use_collective)
      deallocate(packed_int)
    else
      call h5_write_dataset_collect_hyperslab(trim(desc%dataset_name), desc%dest_int, start_idx(1:rank), use_collective)
    endif
  case default
    call stop_timer()
    return
  end select
  call stop_timer()

  elem_size = descriptor_element_size(desc)
  if (rank > 0 .and. elem_size > 0) then
    num_elements = product(dims(1:rank))
    call set_bytes_written_from_array(elem_size*8, num_elements)
  endif

#else

  ! no-op

#endif

  end subroutine write_descriptor_dataset

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

  subroutine recv_surface_movie(tag, tag_src, status, &
                                dump_surf_ux, dump_surf_uy, dump_surf_uz, &
                                surf_frame_ux, surf_frame_uy, surf_frame_uz)

    use specfem_par_movie_hdf5
    use constants, only: CUSTOM_REAL, my_status_size

    implicit none

    integer, intent(in) :: tag, tag_src
    integer, intent(in) :: status(my_status_size)
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_surf_ux
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_surf_uy
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_surf_uz
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: surf_frame_ux
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: surf_frame_uy
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: surf_frame_uz

    integer :: msg_size, ista, data_len, req_dummy

    ! get message size
    call world_get_size_msg(status, msg_size)

    ! compute start index and length within full-frame buffers
    ista = sum(offset_poin(0:tag_src-1))
    data_len = offset_poin(tag_src)

    ! If this rank has no surface movie points, we must still consume the message
    ! (which was already probed) to remove it from the MPI queue
    if (data_len <= 0) then
      call consume_empty_message_inter(tag_src, tag)
      return
    endif

    if (tag == io_tag_surf_ux) then
      call irecvv_cr_inter(dump_surf_ux(1:data_len), msg_size, tag_src, tag, req_dummy)
      surf_frame_ux(ista+1:ista+data_len) = dump_surf_ux(1:data_len)
    else if (tag == io_tag_surf_uy) then
      call irecvv_cr_inter(dump_surf_uy(1:data_len), msg_size, tag_src, tag, req_dummy)
      surf_frame_uy(ista+1:ista+data_len) = dump_surf_uy(1:data_len)
    else if (tag == io_tag_surf_uz) then
      call irecvv_cr_inter(dump_surf_uz(1:data_len), msg_size, tag_src, tag, req_dummy)
      surf_frame_uz(ista+1:ista+data_len) = dump_surf_uz(1:data_len)
    else
      print *, 'Error: unknown surface movie tag in recv_surface_movie'
      stop 'Error: unknown surface movie tag in recv_surface_movie'
    end if

  end subroutine recv_surface_movie

  !-------------------------------------------------------------------------------------------------

  subroutine recv_volume_movie(tag, tag_src, status, &
                               dump_vol1, dump_vol2, dump_vol3, &
                               dump_vol4, dump_vol5, dump_vol6, &
                               vol_frame1, vol_frame2, vol_frame3, &
                               vol_frame4, vol_frame5, vol_frame6)

    use specfem_par
    use specfem_par_movie_hdf5
    use constants, only: CUSTOM_REAL, my_status_size

    implicit none

    integer, intent(in) :: tag, tag_src
    integer, intent(in) :: status(my_status_size)
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol1
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol2
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol3
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol4
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol5
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol6
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame1
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame2
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame3
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame4
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame5
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame6

    integer :: msg_size, ista, data_len, req_dummy

    ! get message size
    call world_get_size_msg(status, msg_size)

    ista = sum(offset_poin_vol(0:tag_src-1))
    data_len = offset_poin_vol(tag_src)

    ! If this rank has no volume movie points, we must still consume the message
    ! (which was already probed) to remove it from the MPI queue
    if (data_len <= 0) then
      call consume_empty_message_inter(tag_src, tag)
      return
    endif

    select case (MOVIE_VOLUME_TYPE)
    case (1,2,3)
      if (tag == io_tag_vol_strain_NN) then
        call irecvv_cr_inter(dump_vol1(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame1(ista+1:ista+data_len) = dump_vol1(1:data_len)
      else if (tag == io_tag_vol_strain_EE) then
        call irecvv_cr_inter(dump_vol2(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame2(ista+1:ista+data_len) = dump_vol2(1:data_len)
      else if (tag == io_tag_vol_strain_ZZ) then
        call irecvv_cr_inter(dump_vol3(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame3(ista+1:ista+data_len) = dump_vol3(1:data_len)
      else if (tag == io_tag_vol_strain_NE) then
        call irecvv_cr_inter(dump_vol4(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame4(ista+1:ista+data_len) = dump_vol4(1:data_len)
      else if (tag == io_tag_vol_strain_NZ) then
        call irecvv_cr_inter(dump_vol5(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame5(ista+1:ista+data_len) = dump_vol5(1:data_len)
      else if (tag == io_tag_vol_strain_EZ) then
        call irecvv_cr_inter(dump_vol6(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame6(ista+1:ista+data_len) = dump_vol6(1:data_len)
      else
        print *, 'Error: unknown volume strain tag in recv_volume_movie'
        stop 'Error: unknown volume strain tag in recv_volume_movie'
      end if

    case (5,6)
      if (tag == io_tag_vol_vec_N) then
        call irecvv_cr_inter(dump_vol1(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame1(ista+1:ista+data_len) = dump_vol1(1:data_len)
      else if (tag == io_tag_vol_vec_E) then
        call irecvv_cr_inter(dump_vol2(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame2(ista+1:ista+data_len) = dump_vol2(1:data_len)
      else if (tag == io_tag_vol_vec_Z) then
        call irecvv_cr_inter(dump_vol3(1:data_len), msg_size, tag_src, tag, req_dummy)
        vol_frame3(ista+1:ista+data_len) = dump_vol3(1:data_len)
      else
        print *, 'Error: unknown volume vector tag in recv_volume_movie'
        stop 'Error: unknown volume vector tag in recv_volume_movie'
      end if

    case default
      print *, 'Error: MOVIE_VOLUME_TYPE not supported in recv_volume_movie'
      stop 'Error: MOVIE_VOLUME_TYPE not supported in recv_volume_movie'
    end select

  end subroutine recv_volume_movie

  !-------------------------------------------------------------------------------------------------

  subroutine recv_volume_norm_movie(tag, tag_src, status, &
                                    dump_vol_norm_cm, dump_vol_norm_oc, dump_vol_norm_ic, &
                                    vol_frame_norm_cm, vol_frame_norm_oc, vol_frame_norm_ic)
  ! Receive handler for norm movie types (7,8,9) which use NGLOB data per region
  ! These have different offset arrays (offset_nglob_cm/oc/ic) from the movie-point types

    use specfem_par
    use specfem_par_movie_hdf5
    use constants, only: CUSTOM_REAL, my_status_size

    implicit none

    integer, intent(in) :: tag, tag_src
    integer, intent(in) :: status(my_status_size)
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol_norm_cm
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol_norm_oc
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: dump_vol_norm_ic
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame_norm_cm
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame_norm_oc
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: vol_frame_norm_ic

    integer :: msg_size, ista, data_len, req_dummy

    ! get message size
    call world_get_size_msg(status, msg_size)

    if (tag == io_tag_vol_norm_cm) then
      ista = sum(offset_nglob_cm(0:tag_src-1))
      data_len = offset_nglob_cm(tag_src)
      if (data_len <= 0) then
        call consume_empty_message_inter(tag_src, tag)
        return
      endif
      call irecvv_cr_inter(dump_vol_norm_cm(1:data_len), msg_size, tag_src, tag, req_dummy)
      vol_frame_norm_cm(ista+1:ista+data_len) = dump_vol_norm_cm(1:data_len)

    else if (tag == io_tag_vol_norm_oc) then
      ista = sum(offset_nglob_oc(0:tag_src-1))
      data_len = offset_nglob_oc(tag_src)
      if (data_len <= 0) then
        call consume_empty_message_inter(tag_src, tag)
        return
      endif
      call irecvv_cr_inter(dump_vol_norm_oc(1:data_len), msg_size, tag_src, tag, req_dummy)
      vol_frame_norm_oc(ista+1:ista+data_len) = dump_vol_norm_oc(1:data_len)

    else if (tag == io_tag_vol_norm_ic) then
      ista = sum(offset_nglob_ic(0:tag_src-1))
      data_len = offset_nglob_ic(tag_src)
      if (data_len <= 0) then
        call consume_empty_message_inter(tag_src, tag)
        return
      endif
      call irecvv_cr_inter(dump_vol_norm_ic(1:data_len), msg_size, tag_src, tag, req_dummy)
      vol_frame_norm_ic(ista+1:ista+data_len) = dump_vol_norm_ic(1:data_len)

    else
      print *, 'Error: unknown volume norm tag in recv_volume_norm_movie'
      stop 'Error: unknown volume norm tag in recv_volume_norm_movie'
    end if

  end subroutine recv_volume_norm_movie

  !-------------------------------------------------------------------------------------------------

  subroutine write_surface_frame(i_frame, it_first_surf, &
                                 surf_frame_ux, surf_frame_uy, surf_frame_uz)

    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
    use constants, only: CUSTOM_REAL, myrank
    use io_bandwidth

    implicit none

    integer, intent(in) :: i_frame
    integer, intent(in) :: it_first_surf
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: surf_frame_ux
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: surf_frame_uy
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: surf_frame_uz

    integer :: it_val
    integer :: data_size
    logical :: use_collective

    if (.not. MOVIE_SURFACE) return
    if (npoints_surf_mov_all_proc <= 0) return

    data_size = CUSTOM_REAL
    ! collective I/O is only meaningful in non-IO-server runs (HDF5_IO_NODES == 0)
    use_collective = H5_COL .and. (HDF5_IO_NODES == 0)

    ! compute actual time-step index for this frame
    it_val = it_first_surf + i_frame * NTSTEP_BETWEEN_FRAMES

    ! construct file and group names
    if (HDF5_IO_NODES > 1) then
      file_name = trim(OUTPUT_FILES)//"/movie_surface.io"//trim(i2c(my_io_id))//".h5"
    else
      file_name = trim(OUTPUT_FILES)//"/movie_surface.h5"
    endif
    group_name = "it_"//trim(i2c(it_val))
    if (HDF5_IO_NODES > 0) then
      ! IO-server mode: only IO ranks touch HDF5 files, use
      ! simple (non-MPI) create/open helper to avoid MPI FAPL.
      call h5_create_or_open_file(file_name)
      call h5_open_or_create_group(group_name)

      ! Ensure datasets exist using helper
      call ensure_dataset_exists(group_name, 'ux', npoints_surf_mov_all_proc)
      call ensure_dataset_exists(group_name, 'uy', npoints_surf_mov_all_proc)
      call ensure_dataset_exists(group_name, 'uz', npoints_surf_mov_all_proc)
    else
      ! non-IO-server runs: all compute ranks participate in I/O
      call open_hdf5_file_for_io(file_name, use_collective, .false.)
      call h5_open_or_create_group(group_name)

      ! Ensure datasets exist using helper
      call ensure_dataset_exists(group_name, 'ux', npoints_surf_mov_all_proc)
      call ensure_dataset_exists(group_name, 'uy', npoints_surf_mov_all_proc)
      call ensure_dataset_exists(group_name, 'uz', npoints_surf_mov_all_proc)
    endif

    ! write full-frame datasets as single hyperslabs
    call start_timer()
    call h5_write_dataset_collect_hyperslab_in_group('ux', surf_frame_ux, (/0/), use_collective)
    call h5_write_dataset_collect_hyperslab_in_group('uy', surf_frame_uy, (/0/), use_collective)
    call h5_write_dataset_collect_hyperslab_in_group('uz', surf_frame_uz, (/0/), use_collective)
    call stop_timer()

    call set_bytes_written_from_array(data_size*8, 3*npoints_surf_mov_all_proc)

    call h5_close_group()
    call close_hdf5_file_for_io(HDF5_IO_NODES > 0)

  end subroutine write_surface_frame

  !-------------------------------------------------------------------------------------------------

  subroutine write_volume_frame(i_frame, it_first_vol, &
                                vol_frame1, vol_frame2, vol_frame3, &
                                vol_frame4, vol_frame5, vol_frame6)

    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
    use constants, only: CUSTOM_REAL, myrank
    use io_bandwidth

    implicit none

    integer, intent(in) :: i_frame
    integer, intent(in) :: it_first_vol
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame1
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame2
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame3
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame4
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame5
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame6

    integer :: it_val
    character(len=2) :: movie_prefix2
    character(len=MAX_STRING_LEN) :: dset_name
    character(len=MAX_STRING_LEN) :: dset_full_name
    integer :: data_size
    logical :: use_collective
    logical :: dset_exists

    if (.not. MOVIE_VOLUME) return
    if (npoints_vol_mov_all_proc <= 0) return

    data_size = CUSTOM_REAL
    ! collective I/O is only meaningful in non-IO-server runs (HDF5_IO_NODES == 0)
    use_collective = H5_COL .and. (HDF5_IO_NODES == 0)

    it_val = it_first_vol + i_frame * NTSTEP_BETWEEN_FRAMES

    if (HDF5_IO_NODES > 1) then
      file_name = trim(OUTPUT_FILES)//"/movie_volume.io"//trim(i2c(my_io_id))//".h5"
    else
      file_name = trim(OUTPUT_FILES)//"/movie_volume.h5"
    endif
    group_name = "it_"//trim(i2c(it_val))
    if (HDF5_IO_NODES > 0) then
      ! IO-server mode: only IO ranks touch HDF5 files, use
      ! simple (non-MPI) create/open helper to avoid MPI FAPL.
      call h5_create_or_open_file(file_name)
      call h5_open_or_create_group(group_name)
    else
      ! non-IO-server runs: all compute ranks participate in I/O
      if (use_collective) then
        call h5_open_file_p_collect(file_name)
      else
        call h5_open_file_p(file_name)
      endif
      call h5_open_group(group_name)
    endif

    select case (MOVIE_VOLUME_TYPE)
    case (1,2,3)
      if (MOVIE_VOLUME_TYPE == 1) then
        movie_prefix2 = 'E '
      else if (MOVIE_VOLUME_TYPE == 2) then
        movie_prefix2 = 'S '
      else
        movie_prefix2 = 'P '
      endif

      ! ensure datasets exist in multi-IO mode
      if (HDF5_IO_NODES > 1) then
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
      endif

      call start_timer()
      dset_name = trim(movie_prefix2)//'NN'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame1, (/0/), use_collective)
      dset_name = trim(movie_prefix2)//'EE'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame2, (/0/), use_collective)
      dset_name = trim(movie_prefix2)//'ZZ'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame3, (/0/), use_collective)
      dset_name = trim(movie_prefix2)//'NE'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame4, (/0/), use_collective)
      dset_name = trim(movie_prefix2)//'NZ'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame5, (/0/), use_collective)
      dset_name = trim(movie_prefix2)//'EZ'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame6, (/0/), use_collective)
      call stop_timer()

      call set_bytes_written_from_array(data_size*8, 6*npoints_vol_mov_all_proc)

    case (5,6)
      if (MOVIE_VOLUME_TYPE == 5) then
        movie_prefix2 = 'DI'
      else
        movie_prefix2 = 'VE'
      endif

      if (HDF5_IO_NODES > 1) then
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
      endif

      call start_timer()
      dset_name = trim(movie_prefix2)//'N'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame1, (/0/), use_collective)
      dset_name = trim(movie_prefix2)//'E'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame2, (/0/), use_collective)
      dset_name = trim(movie_prefix2)//'Z'
      call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), vol_frame3, (/0/), use_collective)
      call stop_timer()

      call set_bytes_written_from_array(data_size*8, 3*npoints_vol_mov_all_proc)

    case default
      print *, 'Error: MOVIE_VOLUME_TYPE not supported in write_volume_frame'
      stop 'Error: MOVIE_VOLUME_TYPE not supported in write_volume_frame'
    end select

    call h5_close_group()
    if (HDF5_IO_NODES > 0) then
      call h5_close_file()
    else
      call h5_close_file_p()
    endif

  end subroutine write_volume_frame

  !-------------------------------------------------------------------------------------------------

  subroutine write_volume_norm_frame(i_frame, it_first_vol, &
                                     vol_frame_norm_cm, vol_frame_norm_oc, vol_frame_norm_ic)
  ! Write handler for norm movie types (7,8,9) which output NGLOB data per region
  ! These write separate datasets:
  !   Type 7 (displnorm): reg1_displ, reg2_displ, reg3_displ
  !   Type 8 (velnorm):   reg1_veloc, reg2_veloc, reg3_veloc
  !   Type 9 (accelnorm): reg1_accel, reg2_accel, reg3_accel

    use specfem_par
    use specfem_par_movie_hdf5
    use manager_hdf5
    use constants, only: CUSTOM_REAL, myrank
    use io_bandwidth

    implicit none

    integer, intent(in) :: i_frame
    integer, intent(in) :: it_first_vol
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame_norm_cm
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame_norm_oc
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vol_frame_norm_ic

    integer :: it_val
    character(len=MAX_STRING_LEN) :: dset_full_name
    character(len=16) :: dset_suffix  ! _displ, _veloc, or _accel
    integer :: data_size, total_bytes
    logical :: use_collective
    logical :: dset_exists

    if (.not. MOVIE_VOLUME) return
    if (MOVIE_VOLUME_TYPE < 7 .or. MOVIE_VOLUME_TYPE > 9) return

    ! determine dataset suffix based on movie type
    select case (MOVIE_VOLUME_TYPE)
    case (7)
      dset_suffix = '_displ'
    case (8)
      dset_suffix = '_veloc'
    case (9)
      dset_suffix = '_accel'
    end select

    data_size = CUSTOM_REAL
    total_bytes = 0
    ! collective I/O is only meaningful in non-IO-server runs (HDF5_IO_NODES == 0)
    use_collective = H5_COL .and. (HDF5_IO_NODES == 0)

    it_val = it_first_vol + i_frame * NTSTEP_BETWEEN_FRAMES

    if (HDF5_IO_NODES > 1) then
      file_name = trim(OUTPUT_FILES)//"/movie_volume.io"//trim(i2c(my_io_id))//".h5"
    else
      file_name = trim(OUTPUT_FILES)//"/movie_volume.h5"
    endif
    group_name = "it_"//trim(i2c(it_val))

    if (HDF5_IO_NODES > 0) then
      ! IO-server mode: only IO ranks touch HDF5 files
      call h5_create_or_open_file(file_name)
      call h5_open_or_create_group(group_name)
    else
      ! non-IO-server runs: all compute ranks participate in I/O
      if (use_collective) then
        call h5_open_file_p_collect(file_name)
      else
        call h5_open_file_p(file_name)
      endif
      call h5_open_group(group_name)
    endif

    call start_timer()

    ! write crust-mantle region
    if (OUTPUT_CRUST_MANTLE .and. npoints_vol_mov_all_proc_cm > 0) then
      if (HDF5_IO_NODES > 1) then
        dset_full_name = trim(group_name)//'/reg1'//trim(dset_suffix)
        call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen_in_group('reg1'//trim(dset_suffix), (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
        endif
      endif
      call h5_write_dataset_collect_hyperslab_in_group('reg1'//trim(dset_suffix), vol_frame_norm_cm, (/0/), use_collective)
      total_bytes = total_bytes + npoints_vol_mov_all_proc_cm
    endif

    ! write outer core region
    if (OUTPUT_OUTER_CORE .and. npoints_vol_mov_all_proc_oc > 0) then
      if (HDF5_IO_NODES > 1) then
        dset_full_name = trim(group_name)//'/reg2'//trim(dset_suffix)
        call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen_in_group('reg2'//trim(dset_suffix), (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
        endif
      endif
      call h5_write_dataset_collect_hyperslab_in_group('reg2'//trim(dset_suffix), vol_frame_norm_oc, (/0/), use_collective)
      total_bytes = total_bytes + npoints_vol_mov_all_proc_oc
    endif

    ! write inner core region
    if (OUTPUT_INNER_CORE .and. npoints_vol_mov_all_proc_ic > 0) then
      if (HDF5_IO_NODES > 1) then
        dset_full_name = trim(group_name)//'/reg3'//trim(dset_suffix)
        call h5_check_dataset_exists(trim(dset_full_name), dset_exists)
        if (.not. dset_exists) then
          call h5_create_dataset_gen_in_group('reg3'//trim(dset_suffix), (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
        endif
      endif
      call h5_write_dataset_collect_hyperslab_in_group('reg3'//trim(dset_suffix), vol_frame_norm_ic, (/0/), use_collective)
      total_bytes = total_bytes + npoints_vol_mov_all_proc_ic
    endif

    call stop_timer()
    call set_bytes_written_from_array(data_size*8, total_bytes)

    call h5_close_group()
    if (HDF5_IO_NODES > 0) then
      call h5_close_file()
    else
      call h5_close_file_p()
    endif

  end subroutine write_volume_norm_frame

  !-------------------------------------------------------------------------------------------------

  subroutine recv_and_write_ford_undo(tag, tag_src, status, &
                                     dump_ford_undo_1d_glob, &
                                     dump_ford_undo_2d_glob, &
                                     dump_ford_undo_4d, &
                                     dump_ford_undo_5d, &
                                     i_snapshot)

    use specfem_par
    use specfem_par_movie_hdf5
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

    integer :: msg_size
    integer :: element_bytes
    integer :: desc_idx
    type(undo_message_descriptor) :: descriptor

    ! get message size
    call world_get_size_msg(status, msg_size)

    desc_idx = find_undo_descriptor_index(tag)
    if (desc_idx == 0) then
      print *, 'Error: unknown tag in recv_and_write_ford_undo', tag
      stop 'Error: unknown tag in recv_and_write_ford_undo'
    endif

    descriptor = undo_descriptors(desc_idx)
    element_bytes = 0

    select case (descriptor%buffer_kind)
    case (UNDO_BUFFER_REAL1D)
      call handle_real1d_descriptor(descriptor, tag_src, tag, msg_size, &
                                    dump_ford_undo_1d_glob, element_bytes)
    case (UNDO_BUFFER_REAL2D)
      call handle_real2d_descriptor(descriptor, tag_src, tag, msg_size, &
                                    dump_ford_undo_2d_glob, element_bytes)
    case (UNDO_BUFFER_REAL4D)
      call handle_real4d_descriptor(descriptor, tag_src, tag, msg_size, &
                                    dump_ford_undo_4d, element_bytes)
    case (UNDO_BUFFER_REAL5D)
      call handle_real5d_descriptor(descriptor, tag_src, tag, msg_size, &
                                    dump_ford_undo_5d, element_bytes)
    case (UNDO_BUFFER_INT1D)
      call handle_int_descriptor(descriptor, tag_src, tag, element_bytes)
    case default
      print *, 'Error: unsupported descriptor buffer kind in recv_and_write_ford_undo', descriptor%buffer_kind
      stop 'Error: unsupported descriptor buffer kind in recv_and_write_ford_undo'
    end select

    if (element_bytes == 0) element_bytes = CUSTOM_REAL

    ! Undo-snapshot bytes are accounted when the buffered descriptors are flushed
    ! to HDF5 in write_descriptor_dataset(). Counting them here would double count
    ! every buffered message in the IO-server path.

  end subroutine recv_and_write_ford_undo

  subroutine write_buffered_undo_snapshot(use_collective)

#ifdef USE_HDF5

    use specfem_par
    use manager_hdf5
    use constants, only: CUSTOM_REAL
    use io_bandwidth

    implicit none

    logical, intent(in) :: use_collective

    integer :: i

    if (.not. IO_storage_task) return
    if (.not. undo_state%file_open) return

    if (HDF5_IO_NODES > 1) then
      call ensure_undo_datasets_in_shard()
    endif

    if (.not. allocated(undo_descriptors)) return

    do i = 1, size(undo_descriptors)
      call write_descriptor_dataset(undo_descriptors(i), use_collective)
    enddo

#endif

  end subroutine write_buffered_undo_snapshot

  subroutine ensure_undo_datasets_in_shard()

#ifdef USE_HDF5

    use specfem_par
    use manager_hdf5
    use constants, only: CUSTOM_REAL

    implicit none

    integer :: i

    if (undo_state%datasets_initialized) return
    if (HDF5_IO_NODES <= 1) return
    if (.not. undo_state%file_open) return
    if (.not. allocated(undo_descriptors)) then
      undo_state%datasets_initialized = .true.
      return
    endif

    do i = 1, size(undo_descriptors)
      call ensure_descriptor_dataset(undo_descriptors(i))
    enddo

    undo_state%datasets_initialized = .true.

#endif

  end subroutine ensure_undo_datasets_in_shard

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

    ! open file and ensure the frame group/datasets exist
    if (HDF5_IO_NODES > 1) then
      ! each IO rank independently opens or creates its own shard file
      call h5_create_or_open_file(file_name)

      ! ensure group and datasets exist locally in this shard file
      call h5_open_or_create_group(group_name)
    else
      if (use_collective) then
        call h5_open_file_p_collect(file_name)
      else
        call h5_open_file_p(file_name)
      endif

      ! single IO node: lazily create or open the frame group as needed
      call h5_open_or_create_group(group_name)
    endif

    ! ensure datasets ux/uy/uz exist in this frame group
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

    ! debug: check that volume movie hyperslab fits into dataset extent
    if (VERBOSE) then
      print *, 'DEBUG recv_and_write_volume_movie: rank', myrank, ' src', tag_src, &
               ' ista =', ista, ' len =', data_len, &
               ' npoints_vol_mov_all_proc =', npoints_vol_mov_all_proc
      call flush_stdout()
    endif

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
      else if (tag == io_tag_vol_vec_E) then
        call irecvv_cr_inter(dump_vol2(1:data_len), msg_size, tag_src, tag, req_dummy)
        dset_name = trim(movie_prefix2)//'E'
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
        call h5_write_dataset_collect_hyperslab_in_group(trim(dset_name), dump_vol2(1:data_len), (/ista/), use_collective)
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
