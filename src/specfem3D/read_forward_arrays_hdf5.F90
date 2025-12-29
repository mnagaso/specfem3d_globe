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

  subroutine read_intermediate_forward_arrays_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity
  use specfem_par_movie_hdf5

#ifdef USE_HDF5
  use manager_hdf5
  use io_server_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! full gravity
  integer :: neq_read,neq1_read

  ! MPI variables
  !integer :: info, comm

  ! TODO HDF5: put offset array creation in a initialization process
  ! offset array
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_cm
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_oc
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_ic
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_soa
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic_soa
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_oc_rot
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_att
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic_att
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_pgrav1

  !character(len=MAX_STRING_LEN) :: file_name

  ! gather the offset arrays
  !call gather_all_all_singlei(size(displ_crust_mantle), offset_nglob_cm, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(displ_inner_core),   offset_nglob_oc, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(displ_outer_core),   offset_nglob_ic, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(epsilondev_xx_crust_mantle,4), offset_nspec_cm_soa, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(epsilondev_xx_inner_core,4),   offset_nspec_ic_soa, NPROCTOT_VAL)
  !if (ROTATION_VAL) then
  !  call gather_all_all_singlei(size(A_array_rotation,4), offset_nspec_oc_rot, NPROCTOT_VAL)
  !endif
  !if (ATTENUATION_VAL) then
  !  call gather_all_all_singlei(size(R_xx_crust_mantle,5), offset_nspec_cm_att, NPROCTOT_VAL)
  !  call gather_all_all_singlei(size(R_xx_inner_core,5),   offset_nspec_ic_att, NPROCTOT_VAL)
  !endif
  !if (FULL_GRAVITY_VAL) then
  !  call gather_all_all_singlei(size(pgrav1), offset_pgrav1, NPROCTOT_VAL)
  !endif

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/dump_all_arrays.h5'

  ! get MPI parameters
  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initialize HDF5
  call h5_initialize() ! called in initialize_mesher()
  ! set MPI
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! open the file
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! read the arrays
  call h5_read_dataset_collect_hyperslab('displ_crust_mantle', displ_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('displ_outer_core',   displ_outer_core,   (/0,sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('displ_inner_core',   displ_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('veloc_crust_mantle', veloc_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('veloc_outer_core',   veloc_outer_core,   (/0,sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('veloc_inner_core',   veloc_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('accel_crust_mantle', accel_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('accel_outer_core',   accel_outer_core,   (/0,sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('accel_inner_core',   accel_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', epsilondev_xx_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', epsilondev_yy_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', epsilondev_xy_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', epsilondev_xz_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', epsilondev_yz_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xx_inner_core',   epsilondev_xx_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yy_inner_core',   epsilondev_yy_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xy_inner_core',   epsilondev_xy_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xz_inner_core',   epsilondev_xz_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yz_inner_core',   epsilondev_yz_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)

  if (ROTATION_VAL) then
    call h5_read_dataset_collect_hyperslab('A_array_rotation', A_array_rotation, &
                                           (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('A_array_rotation', B_array_rotation, &
                                           (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
  endif

  if (ATTENUATION_VAL) then
    call h5_read_dataset_collect_hyperslab('R_xx_crust_mantle', R_xx_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yy_crust_mantle', R_yy_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xy_crust_mantle', R_xy_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xz_crust_mantle', R_xz_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yz_crust_mantle', R_yz_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xx_inner_core',   R_xx_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yy_inner_core',   R_yy_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xy_inner_core',   R_xy_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xz_inner_core',   R_xz_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yz_inner_core',   R_yz_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
  endif

  if (FULL_GRAVITY_VAL) then
    call h5_read_dataset_scalar_collect_hyperslab('neq', neq_read, (/myrank/), H5_COL)
    call h5_read_dataset_scalar_collect_hyperslab('neq1', neq1_read, (/myrank/), H5_COL)
    call h5_read_dataset_collect_hyperslab('pgrav1', pgrav1, (/0,sum(offset_pgrav1(0:myrank-1))/), H5_COL)
  endif

  ! close the file
  call h5_close_file_p()

#else
  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')
#endif

  end subroutine read_intermediate_forward_arrays_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity
  use specfem_par_movie_hdf5

#ifdef USE_HDF5
  use manager_hdf5
  use io_server_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! full gravity
  integer :: b_neq_read, b_neq1_read

  ! MPI variables
  !integer :: info, comm


  ! TODO HDF5: put offset array creation in a initialization process
  ! offset array
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_cm
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_oc
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_ic
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_soa
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic_soa
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_oc_rot
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_att
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic_att
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_pgrav1

  !character(len=MAX_STRING_LEN) :: file_name

  ! gather the offset arrays
  !call gather_all_all_singlei(size(displ_crust_mantle,2), offset_nglob_cm, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(displ_inner_core,2),   offset_nglob_ic, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(displ_outer_core,1),   offset_nglob_oc, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(epsilondev_xx_crust_mantle,4), offset_nspec_cm_soa, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(epsilondev_xx_inner_core,4),   offset_nspec_ic_soa, NPROCTOT_VAL)
  !if (ROTATION_VAL) then
  !  call gather_all_all_singlei(size(A_array_rotation,4), offset_nspec_oc_rot, NPROCTOT_VAL)
  !endif
  !if (ATTENUATION_VAL) then
  !  call gather_all_all_singlei(size(R_xx_crust_mantle,5), offset_nspec_cm_att, NPROCTOT_VAL)
  !  call gather_all_all_singlei(size(R_xx_inner_core,5),   offset_nspec_ic_att, NPROCTOT_VAL)
  !endif
  !if (FULL_GRAVITY_VAL) then
  !  call gather_all_all_singlei(size(pgrav1), offset_pgrav1, NPROCTOT_VAL)
  !endif

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/dump_all_arrays.h5'

  ! get MPI parameters
  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initialize HDF5
  call h5_initialize() ! called in initialize_mesher()
  ! set MPI
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! open the file
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! read the arrays
  call h5_read_dataset_collect_hyperslab('displ_crust_mantle', b_displ_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('displ_outer_core',   b_displ_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('displ_inner_core',   b_displ_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('veloc_crust_mantle', b_veloc_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('veloc_outer_core',   b_veloc_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('veloc_inner_core',   b_veloc_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('accel_crust_mantle', b_accel_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('accel_outer_core',   b_accel_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('accel_inner_core',   b_accel_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', b_epsilondev_xx_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', b_epsilondev_yy_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', b_epsilondev_xy_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', b_epsilondev_xz_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', b_epsilondev_yz_crust_mantle, &
                                         (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xx_inner_core',   b_epsilondev_xx_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yy_inner_core',   b_epsilondev_yy_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xy_inner_core',   b_epsilondev_xy_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_xz_inner_core',   b_epsilondev_xz_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
  call h5_read_dataset_collect_hyperslab('epsilondev_yz_inner_core',   b_epsilondev_yz_inner_core, &
                                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)

  if (ROTATION_VAL) then
    call h5_read_dataset_collect_hyperslab('A_array_rotation', b_A_array_rotation, &
                                           (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('A_array_rotation', b_B_array_rotation, &
                                           (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
  endif

  if (ATTENUATION_VAL) then
    call h5_read_dataset_collect_hyperslab('R_xx_crust_mantle', b_R_xx_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yy_crust_mantle', b_R_yy_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xy_crust_mantle', b_R_xy_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xz_crust_mantle', b_R_xz_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yz_crust_mantle', b_R_yz_crust_mantle, &
                                           (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xx_inner_core',   b_R_xx_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yy_inner_core',   b_R_yy_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xy_inner_core',   b_R_xy_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_xz_inner_core',   b_R_xz_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('R_yz_inner_core',   b_R_yz_inner_core, &
                                           (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
  endif

  if (FULL_GRAVITY_VAL) then
    call h5_read_dataset_scalar_collect_hyperslab('neq', b_neq_read, (/myrank/), H5_COL)
    call h5_read_dataset_scalar_collect_hyperslab('neq1', b_neq1_read, (/myrank/), H5_COL)

    ! check if array sizes match
    if (b_neq_read /= neq) then
      print *,'Error reading forward array for startrun: rank ',myrank,'has read neq =',b_neq_read,' - shoud be ',neq
      call exit_MPI(myrank,'Invalid forward array neq for startrun')
    endif
    if (b_neq1_read /= neq1) then
      print *,'Error reading forward array for startrun: rank ',myrank,'has read neq1 =',b_neq1_read,' - shoud be ',neq1
      call exit_MPI(myrank,'Invalid forward array neq1 for startrun')
    endif

    call h5_read_dataset_collect_hyperslab('pgrav1', b_pgrav1, (/0,sum(offset_pgrav1(0:myrank-1))/), H5_COL)
  endif

  ! close the file
  call h5_close_file_p()

#else
  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')
#endif

  end subroutine read_forward_arrays_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_undoatt_hdf5(iteration_on_subset_tmp)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity
  use specfem_par_movie_hdf5

#ifdef USE_HDF5
  use manager_hdf5
  use io_server_hdf5
#endif

  implicit none

  ! input parameter
  integer, intent(in) :: iteration_on_subset_tmp

#ifdef USE_HDF5

  ! full gravity
  integer :: b_neq_read, b_neq1_read

  logical :: file_exists

  ! MPI variables
  !integer :: info, comm


  ! TODO HDF5: put offset array creation in a initialization process
  ! offset array
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_cm
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_oc
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nglob_ic
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_soa
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic_soa
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_oc_rot
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_att
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_ic_att
  !integer, dimension(0:NPROCTOT_VAL-1) :: offset_pgrav1

  !character(len=MAX_STRING_LEN) :: file_name

  ! gather the offset arrays
  !call gather_all_all_singlei(size(displ_crust_mantle,2), offset_nglob_cm, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(displ_inner_core,2),   offset_nglob_ic, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(displ_outer_core,1),   offset_nglob_oc, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(epsilondev_xx_crust_mantle,4), offset_nspec_cm_soa, NPROCTOT_VAL)
  !call gather_all_all_singlei(size(epsilondev_xx_inner_core,4),   offset_nspec_ic_soa, NPROCTOT_VAL)
  !if (ROTATION_VAL) then
  !  call gather_all_all_singlei(size(A_array_rotation,4), offset_nspec_oc_rot, NPROCTOT_VAL)
  !endif
  !if (ATTENUATION_VAL) then
  !  call gather_all_all_singlei(size(R_xx_crust_mantle,5), offset_nspec_cm_att, NPROCTOT_VAL)
  !  call gather_all_all_singlei(size(R_xx_inner_core,5),   offset_nspec_ic_att, NPROCTOT_VAL)
  !endif
  !if (FULL_GRAVITY_VAL) then
  !  call gather_all_all_singlei(size(pgrav1), offset_pgrav1, NPROCTOT_VAL)
  !endif

  write(file_name, '(a,i6.6,a)') 'save_frame_at',iteration_on_subset_tmp,'.h5'
  file_name = trim(LOCAL_PATH)//'/'//trim(file_name)

  ! for multi-IO runs, check if a single shared checkpoint file
  ! exists; if not, we assume per-IO shard files and read from
  ! the shard corresponding to this rank's IO node.
  if (HDF5_IO_NODES > 1) then
    inquire(file=trim(file_name), exist=file_exists)
  else
    file_exists = .true.
  endif

  if (file_exists .or. HDF5_IO_NODES <= 1) then

    !---- Original single-file parallel HDF5 path ----

    ! get MPI parameters
    call world_get_comm(comm)
    call world_get_info_null(info)

    ! initialize HDF5
    call h5_initialize() ! called in initialize_mesher()
    ! set MPI
    call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

    ! open the file
    if (H5_COL) then
      ! open file
      call h5_open_file_p_collect(file_name)
    else
      ! open file
      call h5_open_file_p(file_name)
    endif

    ! read the arrays
    call h5_read_dataset_collect_hyperslab('displ_crust_mantle', b_displ_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('displ_outer_core',   b_displ_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('displ_inner_core',   b_displ_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('veloc_crust_mantle', b_veloc_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('veloc_outer_core',   b_veloc_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('veloc_inner_core',   b_veloc_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('accel_crust_mantle', b_accel_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('accel_outer_core',   b_accel_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('accel_inner_core',   b_accel_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', b_epsilondev_xx_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', b_epsilondev_yy_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', b_epsilondev_xy_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', b_epsilondev_xz_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', b_epsilondev_yz_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_xx_inner_core',   b_epsilondev_xx_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_yy_inner_core',   b_epsilondev_yy_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_xy_inner_core',   b_epsilondev_xy_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_xz_inner_core',   b_epsilondev_xz_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_read_dataset_collect_hyperslab('epsilondev_yz_inner_core',   b_epsilondev_yz_inner_core, &
                         (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)

    if (ROTATION_VAL) then
      call h5_read_dataset_collect_hyperslab('A_array_rotation', b_A_array_rotation, &
                                             (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('A_array_rotation', b_B_array_rotation, &
                                             (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
    endif

    if (ATTENUATION_VAL) then
      call h5_read_dataset_collect_hyperslab('R_xx_crust_mantle', b_R_xx_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_yy_crust_mantle', b_R_yy_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_xy_crust_mantle', b_R_xy_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_xz_crust_mantle', b_R_xz_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_yz_crust_mantle', b_R_yz_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_xx_inner_core',   b_R_xx_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_yy_inner_core',   b_R_yy_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_xy_inner_core',   b_R_xy_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_xz_inner_core',   b_R_xz_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_read_dataset_collect_hyperslab('R_yz_inner_core',   b_R_yz_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    endif

    if (FULL_GRAVITY_VAL) then
      call h5_read_dataset_scalar_collect_hyperslab('neq', b_neq_read, (/myrank/), H5_COL)
      call h5_read_dataset_scalar_collect_hyperslab('neq1', b_neq1_read, (/myrank/), H5_COL)

      ! check if array sizes match
      if (b_neq_read /= neq) then
        print *,'Error reading forward array for startrun: rank ',myrank,'has read neq =',b_neq_read,' - shoud be ',neq
        call exit_MPI(myrank,'Invalid forward array neq for startrun')
      endif
      if (b_neq1_read /= neq1) then
        print *,'Error reading forward array for startrun: rank ',myrank,'has read neq1 =',b_neq1_read,' - shoud be ',neq1
        call exit_MPI(myrank,'Invalid forward array neq1 for startrun')
      endif

      call h5_read_dataset_collect_hyperslab('pgrav1', b_pgrav1, (/0,sum(offset_pgrav1(0:myrank-1))/), H5_COL)
    endif

    ! close the file
    call h5_close_file_p()

  else

    !---- Sharded per-IO checkpoint path (multi-IO, no single file) ----

    ! Each rank reads its own slice from the shard corresponding to
    ! its IO node. The mapping from compute ranks to IO ids is the
    ! same as used when writing the shards.

    write(file_name, '(a,i6.6,a,a,a)') 'save_frame_at',iteration_on_subset_tmp,'.io',trim(i2c(dest_ionod)),'.h5'
    file_name = trim(LOCAL_PATH)//'/'//trim(file_name)

    ! open shard file with serial HDF5
    call h5_open_file(file_name)

    ! read the arrays (non-collective)
    call h5_read_dataset_collect_hyperslab('displ_crust_mantle', b_displ_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('displ_outer_core',   b_displ_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('displ_inner_core',   b_displ_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('veloc_crust_mantle', b_veloc_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('veloc_outer_core',   b_veloc_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('veloc_inner_core',   b_veloc_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('accel_crust_mantle', b_accel_crust_mantle, (/0,sum(offset_nglob_cm(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('accel_outer_core',   b_accel_outer_core,   (/sum(offset_nglob_oc(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('accel_inner_core',   b_accel_inner_core,   (/0,sum(offset_nglob_ic(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', b_epsilondev_xx_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', b_epsilondev_yy_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', b_epsilondev_xy_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', b_epsilondev_xz_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', b_epsilondev_yz_crust_mantle, &
                                           (/0,0,0,sum(offset_nspec_cm_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_xx_inner_core',   b_epsilondev_xx_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_yy_inner_core',   b_epsilondev_yy_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_xy_inner_core',   b_epsilondev_xy_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_xz_inner_core',   b_epsilondev_xz_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), .false.)
    call h5_read_dataset_collect_hyperslab('epsilondev_yz_inner_core',   b_epsilondev_yz_inner_core, &
                                           (/0,0,0,sum(offset_nspec_ic_soa(0:myrank-1))/), .false.)

    if (ROTATION_VAL) then
      call h5_read_dataset_collect_hyperslab('A_array_rotation', b_A_array_rotation, &
                                             (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('A_array_rotation', b_B_array_rotation, &
                                             (/0,0,0,sum(offset_nspec_oc_rot(0:myrank-1))/), .false.)
    endif

    if (ATTENUATION_VAL) then
      call h5_read_dataset_collect_hyperslab('R_xx_crust_mantle', b_R_xx_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_yy_crust_mantle', b_R_yy_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_xy_crust_mantle', b_R_xy_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_xz_crust_mantle', b_R_xz_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_yz_crust_mantle', b_R_yz_crust_mantle, &
                                             (/0,0,0,0,sum(offset_nspec_cm_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_xx_inner_core',   b_R_xx_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_yy_inner_core',   b_R_yy_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_xy_inner_core',   b_R_xy_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_xz_inner_core',   b_R_xz_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), .false.)
      call h5_read_dataset_collect_hyperslab('R_yz_inner_core',   b_R_yz_inner_core, &
                                             (/0,0,0,0,sum(offset_nspec_ic_att(0:myrank-1))/), .false.)
    endif

    if (FULL_GRAVITY_VAL) then
      call h5_read_dataset_scalar_collect_hyperslab('neq', b_neq_read, (/myrank/), .false.)
      call h5_read_dataset_scalar_collect_hyperslab('neq1', b_neq1_read, (/myrank/), .false.)

      ! check if array sizes match
      if (b_neq_read /= neq) then
        print *,'Error reading forward array for startrun: rank ',myrank,'has read neq =',b_neq_read,' - shoud be ',neq
        call exit_MPI(myrank,'Invalid forward array neq for startrun')
      endif
      if (b_neq1_read /= neq1) then
        print *,'Error reading forward array for startrun: rank ',myrank,'has read neq1 =',b_neq1_read,' - shoud be ',neq1
        call exit_MPI(myrank,'Invalid forward array neq1 for startrun')
      endif

      call h5_read_dataset_collect_hyperslab('pgrav1', b_pgrav1, (/0,sum(offset_pgrav1(0:myrank-1))/), .false.)
    endif

    ! close shard file
    call h5_close_file()

  endif

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = iteration_on_subset_tmp

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')
#endif

  end subroutine read_forward_arrays_undoatt_hdf5
