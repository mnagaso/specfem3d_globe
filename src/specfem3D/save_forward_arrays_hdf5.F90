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


  subroutine save_intermediate_forward_arrays_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

    implicit none

#ifdef USE_HDF5

  if (HDF5_IO_NODES > 0) then
    ! TODO: HDF5 IOSERVER
    print *,'Error: HDF5 IOSERVER not implemented yet for save_intermediate_forward_arrays_hdf5'
    call exit_mpi(myrank,'Error: HDF5 IOSERVER not implemented yet')

  else

    file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/dump_all_arrays.h5'

    ! get MPI parameters
    call world_get_comm(comm)
    call world_get_info_null(info)

    ! initialize HDF5
    call h5_initialize() ! called in initialize_mesher()
    ! set MPI
    call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

    ! create file and datasets by myrank==0
    if (myrank == 0) then
      call h5_create_file(file_name)

      ! create datasets
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
        call h5_create_dataset_gen('neq', (/NPROCTOT_VAL/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('neq1', (/NPROCTOT_VAL/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('pgrav1', (/sum(offset_pgrav1)/), 1, CUSTOM_REAL)
      endif ! FULL_GRAVITY_VAL

      ! close file
      call h5_close_file()
    endif ! myrank == 0

    call synchronize_all()

    ! write data from all ranks
    if (H5_COL) then
      ! open file
      call h5_open_file_p_collect(file_name)
    else
      ! open file
      call h5_open_file_p(file_name)
    endif

    ! write datasets
    call h5_write_dataset_collect_hyperslab('displ_crust_mantle',displ_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_crust_mantle',veloc_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_crust_mantle',accel_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('displ_outer_core', displ_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_outer_core', veloc_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_outer_core', accel_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('displ_inner_core', displ_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_inner_core', veloc_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_inner_core', accel_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', epsilondev_xx_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', epsilondev_yy_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', epsilondev_xy_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', epsilondev_xz_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', epsilondev_yz_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xx_inner_core', epsilondev_xx_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yy_inner_core', epsilondev_yy_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xy_inner_core', epsilondev_xy_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xz_inner_core', epsilondev_xz_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yz_inner_core', epsilondev_yz_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)

    if (ROTATION_VAL) then
      call h5_write_dataset_collect_hyperslab('A_array_rotation', A_array_rotation, &
                                              (/0, 0, 0, sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('B_array_rotation', B_array_rotation, &
                                              (/0, 0, 0, sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
    endif

    if (ATTENUATION_VAL) then
      call h5_write_dataset_collect_hyperslab('R_xx_crust_mantle', R_xx_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yy_crust_mantle', R_yy_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xy_crust_mantle', R_xy_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xz_crust_mantle', R_xz_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yz_crust_mantle', R_yz_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xx_inner_core', R_xx_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yy_inner_core', R_yy_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xy_inner_core', R_xy_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xz_inner_core', R_xz_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yz_inner_core', R_yz_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    endif ! ATTENUATION_VAL

    if (FULL_GRAVITY_VAL) then
      call h5_write_dataset_collect_hyperslab('neq', (/neq/), (/myrank/), H5_COL)
      call h5_write_dataset_collect_hyperslab('neq1', (/neq1/), (/myrank/), H5_COL)
      call h5_write_dataset_collect_hyperslab('pgrav1', pgrav1, (/sum(offset_pgrav1(0:myrank-1))/), H5_COL)
    endif ! FULL_GRAVITY_VAL

    ! close file
    call h5_close_file_p()

  endif ! HDF5_IO_NODES == 0

#else

    print *,'Error: HDF5 not enabled in this version of the code'
    print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
    call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine save_intermediate_forward_arrays_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_forward_arrays_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

#ifdef USE_HDF5

  ! Forward arrays use a simple parallel HDF5 file that is
  ! independent of the IO server. Allow this path even when
  ! HDF5_IO_NODES > 0 instead of aborting.

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/save_forward_arrays.h5'

    ! get MPI parameters
    call world_get_comm(comm)
    call world_get_info_null(info)

    ! initialize HDF5
    call h5_initialize() ! called in initialize_mesher()
    ! set MPI
    call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

    ! create file and datasets by myrank==0
    if (myrank == 0) then
      call h5_create_file(file_name)

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
      call h5_create_dataset_gen('epsilondev_xx_inner_core', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_yy_inner_core', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xy_inner_core', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_xz_inner_core', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('epsilondev_yz_inner_core', &
                                 (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic_soa)/), 4, CUSTOM_REAL)

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
        call h5_create_dataset_gen('neq', (/NPROCTOT_VAL/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('neq1', (/NPROCTOT_VAL/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('pgrav1', (/sum(offset_pgrav1)/), 1, CUSTOM_REAL)
      endif

      ! close file
      call h5_close_file()

    endif ! myrank == 0

    call synchronize_all()

    ! write data from all ranks
    if (H5_COL) then
      ! open file
      call h5_open_file_p_collect(file_name)
    else
      ! open file
      call h5_open_file_p(file_name)
    endif

    ! write datasets
    call h5_write_dataset_collect_hyperslab('displ_crust_mantle',displ_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_crust_mantle',veloc_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_crust_mantle',accel_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('displ_outer_core', displ_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_outer_core', veloc_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_outer_core', accel_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('displ_inner_core', displ_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_inner_core', veloc_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_inner_core', accel_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', epsilondev_xx_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', epsilondev_yy_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', epsilondev_xy_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', epsilondev_xz_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', epsilondev_yz_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xx_inner_core', epsilondev_xx_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yy_inner_core', epsilondev_yy_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xy_inner_core', epsilondev_xy_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xz_inner_core', epsilondev_xz_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yz_inner_core', epsilondev_yz_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    if (ROTATION_VAL) then
      call h5_write_dataset_collect_hyperslab('A_array_rotation', A_array_rotation, &
                                              (/0, 0, 0, sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('B_array_rotation', B_array_rotation, &
                                              (/0, 0, 0, sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
    endif
    if (ATTENUATION_VAL) then
      call h5_write_dataset_collect_hyperslab('R_xx_crust_mantle', R_xx_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yy_crust_mantle', R_yy_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xy_crust_mantle', R_xy_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xz_crust_mantle', R_xz_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yz_crust_mantle', R_yz_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xx_inner_core', R_xx_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yy_inner_core', R_yy_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xy_inner_core', R_xy_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xz_inner_core', R_xz_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yz_inner_core', R_yz_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    endif ! ATTENUATION_VAL

    if (FULL_GRAVITY_VAL) then
      call h5_write_dataset_collect_hyperslab('neq', (/neq/), (/myrank/), H5_COL)
      call h5_write_dataset_collect_hyperslab('neq1', (/neq1/), (/myrank/), H5_COL)
      call h5_write_dataset_collect_hyperslab('pgrav1', pgrav1, (/sum(offset_pgrav1(0:myrank-1))/), H5_COL)
    endif ! FULL_GRAVITY_VAL

    ! close file
    call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine save_forward_arrays_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_forward_arrays_undoatt_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity

  use io_bandwidth

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
  use io_server_hdf5
#endif

    implicit none

#ifdef USE_HDF5

  integer :: req_count

  ! hdf5 i/o server request index
  req_count = 1

  if (HDF5_IO_NODES > 0) then
    ! debug: log start of undo send on compute side
    print *, 'compute undo_send_start: subset', iteration_on_subset, 'rank', myrank, 'dest_ionod', dest_ionod

    ! send data to the IO server
    call isend_cr_inter(displ_crust_mantle,NDIM*offset_nglob_cm(myrank),dest_ionod, &
                        io_tag_ford_undo_d_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(veloc_crust_mantle,NDIM*offset_nglob_cm(myrank),dest_ionod, &
                        io_tag_ford_undo_v_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(accel_crust_mantle,NDIM*offset_nglob_cm(myrank),dest_ionod, &
                        io_tag_ford_undo_a_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(displ_outer_core,offset_nglob_oc(myrank),dest_ionod, &
                        io_tag_ford_undo_d_oc,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(veloc_outer_core,offset_nglob_oc(myrank),dest_ionod, &
                        io_tag_ford_undo_v_oc,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(accel_outer_core,offset_nglob_oc(myrank),dest_ionod, &
                        io_tag_ford_undo_a_oc,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(displ_inner_core,NDIM*offset_nglob_ic(myrank),dest_ionod, &
                        io_tag_ford_undo_d_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(veloc_inner_core,NDIM*offset_nglob_ic(myrank),dest_ionod, &
                        io_tag_ford_undo_v_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(accel_inner_core,NDIM*offset_nglob_ic(myrank),dest_ionod, &
                        io_tag_ford_undo_a_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_xx_crust_mantle,NGLLX*NGLLY*NGLLZ*offset_nspec_cm_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_xx_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_yy_crust_mantle,NGLLX*NGLLY*NGLLZ*offset_nspec_cm_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_yy_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_xy_crust_mantle,NGLLX*NGLLY*NGLLZ*offset_nspec_cm_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_xy_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_xz_crust_mantle,NGLLX*NGLLY*NGLLZ*offset_nspec_cm_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_xz_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_yz_crust_mantle,NGLLX*NGLLY*NGLLZ*offset_nspec_cm_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_yz_cm,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_xx_inner_core,NGLLX*NGLLY*NGLLZ*offset_nspec_ic_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_xx_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_yy_inner_core,NGLLX*NGLLY*NGLLZ*offset_nspec_ic_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_yy_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_xy_inner_core,NGLLX*NGLLY*NGLLZ*offset_nspec_ic_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_xy_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_xz_inner_core,NGLLX*NGLLY*NGLLZ*offset_nspec_ic_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_xz_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1
    call isend_cr_inter(epsilondev_yz_inner_core,NGLLX*NGLLY*NGLLZ*offset_nspec_ic_soa(myrank),dest_ionod, &
                        io_tag_ford_undo_eps_yz_ic,req_dump_ford_undo(req_count))
    req_count = req_count + 1

    if (ROTATION_VAL) then
      call isend_cr_inter(A_array_rotation,NGLLX*NGLLY*NGLLZ*offset_nspec_oc_rot(myrank),dest_ionod, &
                          io_tag_ford_undo_A_rot,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(B_array_rotation,NGLLX*NGLLY*NGLLZ*offset_nspec_oc_rot(myrank),dest_ionod, &
                          io_tag_ford_undo_B_rot,req_dump_ford_undo(req_count))
      req_count = req_count + 1
    endif

    if (ATTENUATION_VAL) then
      call isend_cr_inter(R_xx_crust_mantle,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_cm_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_xx_cm,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_yy_crust_mantle,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_cm_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_yy_cm,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_xy_crust_mantle,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_cm_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_xy_cm,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_xz_crust_mantle,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_cm_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_xz_cm,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_yz_crust_mantle,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_cm_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_yz_cm,req_dump_ford_undo(req_count))
      req_count = req_count + 1

      call isend_cr_inter(R_xx_inner_core,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_ic_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_xx_ic,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_yy_inner_core,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_ic_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_yy_ic,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_xy_inner_core,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_ic_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_xy_ic,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_xz_inner_core,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_ic_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_xz_ic,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(R_yz_inner_core,NGLLX*NGLLY*NGLLZ*N_SLS*offset_nspec_ic_att(myrank),dest_ionod, &
                          io_tag_ford_undo_R_yz_ic,req_dump_ford_undo(req_count))
      req_count = req_count + 1
    endif ! ATTENUATION_VAL

    if (FULL_GRAVITY_VAL) then
      call isend_i_inter((/neq/),NPROCTOT_VAL,dest_ionod,io_tag_ford_undo_neq,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_i_inter((/neq1/),NPROCTOT_VAL,dest_ionod,io_tag_ford_undo_neq1,req_dump_ford_undo(req_count))
      req_count = req_count + 1
      call isend_cr_inter(pgrav1,offset_pgrav1(myrank),dest_ionod,io_tag_ford_undo_pgrav1,req_dump_ford_undo(req_count))
      req_count = req_count + 1
    endif ! FULL_GRAVITY_VAL

    ! hdf5 i/o server
    ! store the number of mpi_isend reqs
    n_req_ford_undo = req_count - 1

    ! wait for current subset's sends to complete before returning to time loop
    ! (prevents MPI resource conflict between inter-communicator and compute-only communicator)
    call wait_all_send()

    ! debug: log completion of undo send on compute side
    print *, 'compute undo_send_done: subset', iteration_on_subset, 'rank', myrank, 'dest_ionod', dest_ionod, &
         'n_req_ford_undo', n_req_ford_undo

  else
    write(file_name, '(a,i6.6,a)') 'save_frame_at',iteration_on_subset,'.h5'
    file_name = trim(LOCAL_PATH)//'/'//trim(file_name)

    ! get MPI parameters
    call world_get_comm(comm)
    call world_get_info_null(info)

    ! initialize HDF5
    call h5_initialize() ! called in initialize_mesher()
    ! set MPI
    call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

    ! create file and datasets by myrank==0
    if (myrank == 0) then
      call h5_create_file(file_name)

      ! create datasets
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

      ! close file
      call h5_close_file()
    endif ! myrank == 0

    call synchronize_all()

    ! initialize timer
    call initialize_bytes_written()
    call start_timer()

    ! write data from all ranks
    if (H5_COL) then
      ! open file
      call h5_open_file_p_collect(file_name)
    else
      ! open file
      call h5_open_file_p(file_name)
    endif

    ! write datasets
    call h5_write_dataset_collect_hyperslab('displ_crust_mantle',displ_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_crust_mantle',veloc_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_crust_mantle',accel_crust_mantle,(/0,sum(offset_nglob_cm(0:myrank-1))/),H5_COL)
    call h5_write_dataset_collect_hyperslab('displ_outer_core', displ_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_outer_core', veloc_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_outer_core', accel_outer_core, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('displ_inner_core', displ_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('veloc_inner_core', veloc_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('accel_inner_core', accel_inner_core, (/0, sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xx_crust_mantle', epsilondev_xx_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yy_crust_mantle', epsilondev_yy_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xy_crust_mantle', epsilondev_xy_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xz_crust_mantle', epsilondev_xz_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yz_crust_mantle', epsilondev_yz_crust_mantle, &
                                            (/0, 0, 0, sum(offset_nspec_cm_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xx_inner_core', epsilondev_xx_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yy_inner_core', epsilondev_yy_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xy_inner_core', epsilondev_xy_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_xz_inner_core', epsilondev_xz_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('epsilondev_yz_inner_core', epsilondev_yz_inner_core, &
                                            (/0, 0, 0, sum(offset_nspec_ic_soa(0:myrank-1))/), H5_COL)

    if (ROTATION_VAL) then
      call h5_write_dataset_collect_hyperslab('A_array_rotation', A_array_rotation, &
                                              (/0, 0, 0, sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('B_array_rotation', B_array_rotation, &
                                              (/0, 0, 0, sum(offset_nspec_oc_rot(0:myrank-1))/), H5_COL)
    endif

    if (ATTENUATION_VAL) then
      call h5_write_dataset_collect_hyperslab('R_xx_crust_mantle', R_xx_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yy_crust_mantle', R_yy_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xy_crust_mantle', R_xy_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xz_crust_mantle', R_xz_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yz_crust_mantle', R_yz_crust_mantle, &
                                              (/0, 0, 0, 0, sum(offset_nspec_cm_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xx_inner_core', R_xx_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yy_inner_core', R_yy_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xy_inner_core', R_xy_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_xz_inner_core', R_xz_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('R_yz_inner_core', R_yz_inner_core, &
                                              (/0, 0, 0, 0, sum(offset_nspec_ic_att(0:myrank-1))/), H5_COL)
    endif ! ATTENUATION_VAL

    if (FULL_GRAVITY_VAL) then
      call h5_write_dataset_collect_hyperslab('neq', (/neq/), (/myrank/), H5_COL)
      call h5_write_dataset_collect_hyperslab('neq1', (/neq1/), (/myrank/), H5_COL)
      call h5_write_dataset_collect_hyperslab('pgrav1', pgrav1, (/sum(offset_pgrav1(0:myrank-1))/), H5_COL)
    endif ! FULL_GRAVITY_VAL

    ! close file
    call h5_close_file_p()

    ! end timer
    call stop_timer()
    ! calculate the total bytes written to disk
    call set_bytes_written_from_array(storage_size(displ_crust_mantle), size(displ_crust_mantle))
    call set_bytes_written_from_array(storage_size(veloc_crust_mantle), size(veloc_crust_mantle))
    call set_bytes_written_from_array(storage_size(accel_crust_mantle), size(accel_crust_mantle))
    call set_bytes_written_from_array(storage_size(displ_outer_core), size(displ_outer_core))
    call set_bytes_written_from_array(storage_size(veloc_outer_core), size(veloc_outer_core))
    call set_bytes_written_from_array(storage_size(accel_outer_core), size(accel_outer_core))
    call set_bytes_written_from_array(storage_size(displ_inner_core), size(displ_inner_core))
    call set_bytes_written_from_array(storage_size(veloc_inner_core), size(veloc_inner_core))
    call set_bytes_written_from_array(storage_size(accel_inner_core), size(accel_inner_core))
    call set_bytes_written_from_array(storage_size(epsilondev_xx_crust_mantle), size(epsilondev_xx_crust_mantle))
    call set_bytes_written_from_array(storage_size(epsilondev_yy_crust_mantle), size(epsilondev_yy_crust_mantle))
    call set_bytes_written_from_array(storage_size(epsilondev_xy_crust_mantle), size(epsilondev_xy_crust_mantle))
    call set_bytes_written_from_array(storage_size(epsilondev_xz_crust_mantle), size(epsilondev_xz_crust_mantle))
    call set_bytes_written_from_array(storage_size(epsilondev_yz_crust_mantle), size(epsilondev_yz_crust_mantle))
    call set_bytes_written_from_array(storage_size(epsilondev_xx_inner_core), size(epsilondev_xx_inner_core))
    call set_bytes_written_from_array(storage_size(epsilondev_yy_inner_core), size(epsilondev_yy_inner_core))
    call set_bytes_written_from_array(storage_size(epsilondev_xy_inner_core), size(epsilondev_xy_inner_core))
    call set_bytes_written_from_array(storage_size(epsilondev_xz_inner_core), size(epsilondev_xz_inner_core))
    call set_bytes_written_from_array(storage_size(epsilondev_yz_inner_core), size(epsilondev_yz_inner_core))

    if (ROTATION_VAL) then
      call set_bytes_written_from_array(storage_size(A_array_rotation), size(A_array_rotation))
      call set_bytes_written_from_array(storage_size(B_array_rotation), size(B_array_rotation))
    endif

    if (ATTENUATION_VAL) then
      call set_bytes_written_from_array(storage_size(R_xx_crust_mantle), size(R_xx_crust_mantle))
      call set_bytes_written_from_array(storage_size(R_yy_crust_mantle), size(R_yy_crust_mantle))
      call set_bytes_written_from_array(storage_size(R_xy_crust_mantle), size(R_xy_crust_mantle))
      call set_bytes_written_from_array(storage_size(R_xz_crust_mantle), size(R_xz_crust_mantle))
      call set_bytes_written_from_array(storage_size(R_yz_crust_mantle), size(R_yz_crust_mantle))
      call set_bytes_written_from_array(storage_size(R_xx_inner_core), size(R_xx_inner_core))
      call set_bytes_written_from_array(storage_size(R_yy_inner_core), size(R_yy_inner_core))
      call set_bytes_written_from_array(storage_size(R_xy_inner_core), size(R_xy_inner_core))
      call set_bytes_written_from_array(storage_size(R_xz_inner_core), size(R_xz_inner_core))
      call set_bytes_written_from_array(storage_size(R_yz_inner_core), size(R_yz_inner_core))
    endif

    if (FULL_GRAVITY_VAL) then
      call set_bytes_written_from_array(storage_size(neq), 1)
      call set_bytes_written_from_array(storage_size(neq1), 1)
      call set_bytes_written_from_array(storage_size(pgrav1), size(pgrav1))
    endif

    ! calculate bandwidth
    call calculate_bandwidth_all_procs()

  endif ! HDF5_IO_NODES == 0

#else

    print *,'Error: HDF5 not enabled in this version of the code'
    print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
    call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine save_forward_arrays_undoatt_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_forward_model_at_shifted_frequency_hdf5(factor_scale_relaxed_crust_mantle,factor_scale_relaxed_inner_core)

  use constants

  use specfem_par_crustmantle
  use specfem_par_innercore

#ifdef USE_HDF5
  use shared_parameters, only: R_PLANET,RHOAV,LOCAL_PATH,TRANSVERSE_ISOTROPY,H5_COL,HDF5_IO_NODES
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

  real(kind=CUSTOM_REAL),dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL) :: factor_scale_relaxed_crust_mantle
  real(kind=CUSTOM_REAL),dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL) :: factor_scale_relaxed_inner_core

#ifdef USE_HDF5

  ! local parameters
  integer :: ier
  real(kind=CUSTOM_REAL) :: scaleval1,scale_factor_r
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: muv_shifted,muh_shifted
  integer :: i,j,k,ispec

  ! debug
  logical, parameter :: OUTPUT_RELAXED_MODEL = .false.

  if (HDF5_IO_NODES > 0) then

    print *, 'IO_SERVER not supported yet for save_forward_model_at_shifted_frequency_hdf5'
    call exit_mpi(myrank,'IO_SERVER not supported yet for save_forward_model_at_shifted_frequency_hdf5')

  else

    ! file name
    file_name = trim(LOCAL_PATH)//'/'//'model_shifted.h5'

    ! create file and datasets by myrank==0
    if (myrank == 0) then
      call h5_create_file(file_name)

      if (NSPEC_CRUST_MANTLE > 0) then
        ! safety check
        if (ANISOTROPIC_3D_MANTLE_VAL) &
          call exit_mpi(myrank,'ANISOTROPIC_3D_MANTLE not supported yet for shifted model file output')

        ! group name
        write(group_name, "('reg',i1)") IREGION_CRUST_MANTLE
        call h5_create_group(group_name)
        call h5_open_group(group_name)

        ! create datasets
        call h5_create_dataset_gen_in_group('muv_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen_in_group('muh_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

        if (TRANSVERSE_ISOTROPY) then
          call h5_create_dataset_gen_in_group('vpv_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('vph_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('vsv_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('vsh_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        else ! isotropic
          call h5_create_dataset_gen_in_group('vp_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('vs_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        endif

        ! close group
        call h5_close_group()
      endif ! NSPEC_CRUST_MANTLE > 0

      if (NSPEC_INNER_CORE > 0) then
        if (ANISOTROPIC_INNER_CORE_VAL) then
          call exit_mpi(myrank,'ANISOTROPIC_INNER_CORE not supported yet for shifted model file output')
        else
          ! only isotropic inner core supported

          ! group name
          write(group_name, "('reg',i1)") IREGION_INNER_CORE
          call h5_create_group(group_name)
          call h5_open_group(group_name)

          ! create datasets
          call h5_create_dataset_gen_in_group('vp_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('vs_shifted', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_ic)/), 4, CUSTOM_REAL)
        endif

        ! close group
        call h5_close_group()
      endif

      ! output relaxed model values
      if (OUTPUT_RELAXED_MODEL) then
        ! checks
        if (.not. TRANSVERSE_ISOTROPY) stop 'Outputting relaxed model requires TRANSVERSE_ISOTROPY'

        if (NSPEC_CRUST_MANTLE > 0) then
          ! group name
          write(group_name, "('reg',i1)") IREGION_CRUST_MANTLE
          call h5_create_group(group_name)
          ! open group
          call h5_open_group(group_name)

          ! create datasets
          call h5_create_dataset_gen_in_group('muv_relaxed', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('muh_relaxed', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('kappav_relaxed', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen_in_group('kappah_relaxed', (/NGLLX, NGLLY, NGLLZ, sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        endif

      endif ! OUTPUT_RELAXED_MODEL

      ! close file
      call h5_close_file()
    endif ! myrank == 0

    call synchronize_all()

    !
    ! write data from all ranks
    !
    if (H5_COL) then
      ! open file
      call h5_open_file_p_collect(file_name)
    else
      ! open file
      call h5_open_file_p(file_name)
    endif

    ! scaling factors to re-dimensionalize units
    scaleval1 = real( sqrt(PI*GRAV*RHOAV)*(R_PLANET/1000.0d0), kind=CUSTOM_REAL)  ! velocities

    if (NSPEC_CRUST_MANTLE > 0) then
      ! open group
      write(group_name, "('reg',i1)") IREGION_CRUST_MANTLE
      call h5_open_group(group_name)

      ! uses temporary array
      allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
               muv_shifted(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
               muh_shifted(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), stat=ier)
      if (ier /= 0) stop 'Error allocating temp_store array'
      temp_store(:,:,:,:) = 0._CUSTOM_REAL

      ! safety check
      if (ANISOTROPIC_3D_MANTLE_VAL) &
        call exit_mpi(myrank,'ANISOTROPIC_3D_MANTLE not supported yet for shifted model file output')

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  shifted model files in directory: ',trim(LOCAL_PATH)
      endif

      ! user output
      if (myrank == 0) write(IMAIN,*) '  crust/mantle:'

      ! moduli (muv,muh) are at relaxed values (only Qmu implemented),
      ! scales back to have values at center frequency
      muv_shifted(:,:,:,:) = muvstore_crust_mantle(:,:,:,:)
      muh_shifted(:,:,:,:) = muhstore_crust_mantle(:,:,:,:)
      do ispec = 1,NSPEC_CRUST_MANTLE
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
                scale_factor_r = factor_scale_relaxed_crust_mantle(i,j,k,ispec)
              else
                scale_factor_r = factor_scale_relaxed_crust_mantle(1,1,1,ispec)
              endif
              ! scaling back from relaxed to values at shifted frequency
              ! (see in prepare_attenuation.f90 for how muv,muh are scaled to become relaxed moduli)
              ! muv
              muv_shifted(i,j,k,ispec) = muv_shifted(i,j,k,ispec) / scale_factor_r
              ! muh
              if (ispec_is_tiso_crust_mantle(ispec)) then
                muh_shifted(i,j,k,ispec) = muh_shifted(i,j,k,ispec) / scale_factor_r
              endif
            enddo
          enddo
        enddo
      enddo

      if (TRANSVERSE_ISOTROPY) then
        ! vpv (at relaxed values)
        temp_store(:,:,:,:) = sqrt((kappavstore_crust_mantle(:,:,:,:) &
                              + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                              * scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vpv_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        ! vph
        temp_store(:,:,:,:) = sqrt((kappahstore_crust_mantle(:,:,:,:) &
                              + FOUR_THIRDS * muh_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                              * scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vph_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        ! vsv
        temp_store(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vsv_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        ! vsh
        temp_store(:,:,:,:) = sqrt( muh_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vsh_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      else ! isotropic
        ! vp
        temp_store(:,:,:,:) = sqrt((kappavstore_crust_mantle(:,:,:,:) &
                              + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                              * scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vp_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        ! vs
        temp_store(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vs_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

      endif

      deallocate(temp_store,muv_shifted,muh_shifted)

      ! close group
      call h5_close_group()
    endif ! NSPEC_CRUST_MANTLE > 0

    if (NSPEC_INNER_CORE > 0) then
      ! open group
      write(group_name, "('reg',i1)") IREGION_INNER_CORE
      call h5_open_group(group_name)

      ! uses temporary array
      allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
               muv_shifted(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), stat=ier)
      if (ier /= 0) stop 'Error allocating temp_store array'
      temp_store(:,:,:,:) = 0._CUSTOM_REAL

      ! user output
      if (myrank == 0) write(IMAIN,*) '  inner core:'

      ! moduli (muv,muh) are at relaxed values, scale back to have shifted values at center frequency
      muv_shifted(:,:,:,:) = muvstore_inner_core(:,:,:,:)
      do ispec = 1,NSPEC_INNER_CORE
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
                scale_factor_r = factor_scale_relaxed_inner_core(i,j,k,ispec)
              else
                scale_factor_r = factor_scale_relaxed_inner_core(1,1,1,ispec)
              endif

              ! inverts to scale relaxed back to shifted factor
              ! scaling back from relaxed to values at shifted frequency
              ! (see in prepare_attenuation.f90 for how muv,muh are scaled to become relaxed moduli)
              ! muv
              muv_shifted(i,j,k,ispec) = muv_shifted(i,j,k,ispec) / scale_factor_r
            enddo
          enddo
        enddo
      enddo

      if (ANISOTROPIC_INNER_CORE_VAL) then
        call exit_mpi(myrank,'ANISOTROPIC_INNER_CORE not supported yet for shifted model file output')
      else
        ! isotropic model
        ! vp
        temp_store(:,:,:,:) = sqrt((kappavstore_inner_core(:,:,:,:) &
                              + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_inner_core(:,:,:,:)) &
                              * scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vp_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_ic(0:myrank-1))/), H5_COL)
        ! vs
        temp_store(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_inner_core(:,:,:,:) )*scaleval1
        call h5_write_dataset_collect_hyperslab_in_group('vs_shifted', temp_store, &
                                                         (/0, 0, 0, sum(offset_nspec_ic(0:myrank-1))/), H5_COL)
      endif

      deallocate(temp_store,muv_shifted)

      ! close group
      call h5_close_group()

    endif ! NSPEC_INNER_CORE > 0

    if (OUTPUT_RELAXED_MODEL) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  outputting relaxed model:'
        call flush_IMAIN()
      endif
      ! checks
      if (.not. TRANSVERSE_ISOTROPY) stop 'Outputting relaxed model requires TRANSVERSE_ISOTROPY'

      ! scaling factor to re-dimensionalize units
      ! the scale of GPa--[g/cm^3][(km/s)^2]
      scaleval1 = real( ((sqrt(PI*GRAV*RHOAV)*R_PLANET/1000.d0)**2)*(RHOAV/1000.d0), kind=CUSTOM_REAL) ! moduli GPa

      if (NSPEC_CRUST_MANTLE > 0) then
        ! open group
        write(group_name, "('reg',i1)") IREGION_CRUST_MANTLE
        call h5_open_group(group_name)

        ! muv_relaxed
        call h5_write_dataset_collect_hyperslab_in_group('muv_relaxed', muvstore_crust_mantle*scaleval1, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        ! muh_relaxed
        call h5_write_dataset_collect_hyperslab_in_group('muh_relaxed', muhstore_crust_mantle*scaleval1, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        ! kappav_relaxed
        call h5_write_dataset_collect_hyperslab_in_group('kappav_relaxed', kappavstore_crust_mantle*scaleval1, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        ! kappah_relaxed
        call h5_write_dataset_collect_hyperslab_in_group('kappah_relaxed', kappahstore_crust_mantle*scaleval1, &
                                                         (/0, 0, 0, sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

        ! close group
        call h5_close_group()

      endif ! NSPEC_CRUST_MANTLE > 0

    endif ! OUTPUT_RELAXED_MODEL

    ! close file
    call h5_close_file_p()

  endif ! HDF5_IO_NODES == 0

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(factor_scale_relaxed_crust_mantle,kind=4)
  idummy = size(factor_scale_relaxed_inner_core,kind=4)

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine save_forward_model_at_shifted_frequency_hdf5
