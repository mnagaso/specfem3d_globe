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

  ! initialize the kernels.h5 file
  subroutine initialize_kernels_hdf5()

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    logical :: file_exists
    integer :: stat

    file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

    ! erase the file if it exists
    if (myrank == 0) then
      inquire(file=file_name, exist=file_exists)
      if (file_exists) then
        ! remove the file
        open(unit=10101010, file=file_name, iostat=stat, status='old')
        if (stat == 0) then
          close(unit=10101010, status='delete')
        else
          print *,'Error: could not delete file ',file_name
          call exit_mpi(myrank,'Error: could not delete file')
        end if
      endif
    endif

    call synchronize_all()

  end subroutine initialize_kernels_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_strength_noise_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

#ifdef USE_HDF5

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    if (.not. HDF5_KERNEL_VIS) then
      call h5_create_dataset_gen('sigma_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
    else
      call h5_create_dataset_gen('sigma_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
    endif
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! write data
  if (.not. HDF5_KERNEL_VIS) then
    call h5_write_dataset_collect_hyperslab('sigma_kernel', sigma_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
  else
    call write_array3dspec_as_1d_hdf5('sigma_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      sigma_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
  endif

  call synchronize_all()

  ! close hdf5
  call h5_close_file_p()

  !

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_strength_noise_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_cm_ani_hdf5(alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
                                   betav_kl_crust_mantle,betah_kl_crust_mantle, &
                                   eta_kl_crust_mantle, &
                                   bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
                                   bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle, &
                                   Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle)


  use specfem_par
  use specfem_par_crustmantle

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

  ! input Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
      betav_kl_crust_mantle,betah_kl_crust_mantle, &
      eta_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
      bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle

#ifdef USE_HDF5

  ! check if anything to do
  if (.not. ANISOTROPIC_KL) return

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! create or open file
    call h5_create_or_open_file(file_name)
    ! create dataset

    if (.not. HDF5_KERNEL_VIS) then
      if (SAVE_TRANSVERSE_KL_ONLY) then
        call h5_create_dataset_gen('alphav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('alphah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('eta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

        call h5_create_dataset_gen('bulk_c_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

        call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

      else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
        call h5_create_dataset_gen('alphav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('alphah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

        call h5_create_dataset_gen('bulk_c_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betav_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betah_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

        call h5_create_dataset_gen('eta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

        call h5_create_dataset_gen('Gc_prime_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('Gs_prime_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

        ! check isotropic kernel
        if (.false.) then
          call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('bulk_beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        endif
        ! check anisotropic kernels
        if (.false.) then
          call h5_create_dataset_gen('A_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('C_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('L_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('N_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('F_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Gc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Gs_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Jc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Kc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Mc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Bc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Hc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Ec_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
          call h5_create_dataset_gen('Dc_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        endif

      else
        ! fully anisotropic kernels
        call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        call h5_create_dataset_gen('cijkl_kernel', (/21,NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 5, CUSTOM_REAL)

      endif

    else ! HDF5_KERNEL_VIS
      if (SAVE_TRANSVERSE_KL_ONLY) then
        call h5_create_dataset_gen('alphav_kernel',     (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('alphah_kernel',     (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('betav_kernel',      (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('betah_kernel',      (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('eta_kernel',        (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('rho_kernel',        (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

        call h5_create_dataset_gen('bulk_c_kernel',     (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betav_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betah_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

        call h5_create_dataset_gen('alpha_kernel',      (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('beta_kernel',       (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_beta_kernel',  (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

      else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
        call h5_create_dataset_gen('alphav_kernel',     (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('alphah_kernel',     (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('betav_kernel',      (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('betah_kernel',      (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

        call h5_create_dataset_gen('bulk_c_kernel',     (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betav_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('bulk_betah_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

        call h5_create_dataset_gen('eta_kernel',        (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('rho_kernel',        (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

        call h5_create_dataset_gen('Gc_prime_kernel',   (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
        call h5_create_dataset_gen('Gs_prime_kernel',   (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

      else
        ! fully anisotropic kernels
        !call h5_create_dataset_gen('rho_kernel',      (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
        !call h5_create_dataset_gen('cijkl_kernel', (/21,NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 5, CUSTOM_REAL)
      endif


    endif
    ! close file
    call h5_close_file()

  endif ! myrank == 0

  ! synchronize all
  call synchronize_all()

  ! write data from all ranks
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! write data
  if (.not. HDF5_KERNEL_VIS) then

    if (SAVE_TRANSVERSE_KL_ONLY) then
      call h5_write_dataset_collect_hyperslab('alphav_kernel', alphav_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('alphah_kernel', alphah_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('betav_kernel', betav_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('betah_kernel', betah_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('eta_kernel', eta_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

      call h5_write_dataset_collect_hyperslab('bulk_c_kernel', bulk_c_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('bulk_betav_kernel', bulk_betav_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('bulk_betah_kernel', bulk_betah_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

      call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('bulk_beta_kernel', bulk_beta_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

    else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
      call h5_write_dataset_collect_hyperslab('alphav_kernel', alphav_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('alphah_kernel', alphah_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('betav_kernel', betav_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('betah_kernel', betah_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

      call h5_write_dataset_collect_hyperslab('bulk_c_kernel', bulk_c_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('bulk_betav_kernel', bulk_betav_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('bulk_betah_kernel', bulk_betah_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

      call h5_write_dataset_collect_hyperslab('eta_kernel', eta_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

      call h5_write_dataset_collect_hyperslab('Gc_prime_kernel', Gc_prime_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('Gs_prime_kernel', Gs_prime_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

      ! check isotropic kernel
      if (.false.) then
        call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_crust_mantle, &
                                                (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_crust_mantle, &
                                                (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
        call h5_write_dataset_collect_hyperslab('bulk_beta_kernel', bulk_beta_kl_crust_mantle, &
                                                (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      endif

      ! check anisotropic kernels
      !if (.false.) then
      !  call h5_write_dataset_collect_hyperslab('A_kernel', A_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('C_kernel', C_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('L_kernel', L_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('N_kernel', N_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('F_kernel', F_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Gc_kernel', Gc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Gs_kernel', Gs_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Jc_kernel', Jc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Kc_kernel', Kc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Mc_kernel', Mc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Bc_kernel', Bc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Hc_kernel', Hc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Ec_kernel', Ec_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !  call h5_write_dataset_collect_hyperslab('Dc_kernel', Dc_kl_crust_mantle, (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      !endif

    else

      ! fully anisotropic kernels
      call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                              (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
      call h5_write_dataset_collect_hyperslab('cijkl_kernel', cijkl_kl_crust_mantle, &
                                              (/0,0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

    endif

  else ! HDF5_KERNEL_VIS
    if (SAVE_TRANSVERSE_KL_ONLY) then
      call write_array3dspec_as_1d_hdf5('alphav_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      alphav_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('alphah_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      alphah_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('betav_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      betav_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('betah_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      betah_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('eta_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      eta_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('rho_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      rho_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('bulk_c_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      bulk_c_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('bulk_betav_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      bulk_betav_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('bulk_betah_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      bulk_betah_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('alpha_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      alpha_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('beta_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      beta_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('bulk_beta_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      bulk_beta_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)

    else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
      call write_array3dspec_as_1d_hdf5('alphav_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      alphav_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('alphah_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      alphah_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('betav_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      betav_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('betah_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      betah_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('bulk_c_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      bulk_c_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('bulk_betav_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      bulk_betav_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('bulk_betah_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      bulk_betah_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('eta_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      eta_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('rho_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      rho_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('Gc_prime_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      Gc_prime_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)
      call write_array3dspec_as_1d_hdf5('Gs_prime_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                      Gs_prime_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                      ibool_crust_mantle, .false.)

    else
      ! fully anisotropic kernels ! not implemented
    endif

  endif

  call synchronize_all()

  ! close hdf5
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(alphav_kl_crust_mantle,kind=4)
  idummy = size(alphah_kl_crust_mantle,kind=4)
  idummy = size(betav_kl_crust_mantle,kind=4)
  idummy = size(betah_kl_crust_mantle,kind=4)
  idummy = size(eta_kl_crust_mantle,kind=4)
  idummy = size(bulk_c_kl_crust_mantle,kind=4)
  idummy = size(bulk_beta_kl_crust_mantle,kind=4)
  idummy = size(bulk_betav_kl_crust_mantle,kind=4)
  idummy = size(bulk_betah_kl_crust_mantle,kind=4)
  idummy = size(Gc_prime_kl_crust_mantle,kind=4)
  idummy = size(Gs_prime_kl_crust_mantle,kind=4)

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_cm_ani_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_cm_iso_hdf5(mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
                                    bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle)

  use specfem_par
  use specfem_par_crustmantle

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

  ! Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle

#ifdef USE_HDF5

  ! checks if anything to do
  if (ANISOTROPIC_KL) return

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)

    ! create dataset
    if (.not. HDF5_KERNEL_VIS) then

      call h5_create_dataset_gen('rhonotprime_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('kappa_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('mu_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

      call h5_create_dataset_gen('bulk_c_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('bulk_beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)

    else ! HDF5_KERNEL_VIS

      call h5_create_dataset_gen('rhonotprime_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('kappa_kernel',       (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('mu_kernel',          (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

      call h5_create_dataset_gen('rho_kernel',         (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('alpha_kernel',       (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('beta_kernel',        (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

      call h5_create_dataset_gen('bulk_c_kernel',      (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('bulk_beta_kernel',   (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)

    endif

    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! write data
  if (.not. HDF5_KERNEL_VIS) then
    call h5_write_dataset_collect_hyperslab('rhonotprime_kernel', rhonotprime_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('kappa_kernel', kappa_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('mu_kernel', mu_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

    call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)

    call h5_write_dataset_collect_hyperslab('bulk_c_kernel', bulk_c_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('bulk_beta_kernel', bulk_beta_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
  else
!print*, myrank, 'DEBUG 0002a'
!print*, myrank, 'DEBUG shape of rhonotprime_kl_crust_mantle', shape(rhonotprime_kl_crust_mantle)
!print*, myrank, 'DEBUG offset_nspec_cm', offset_nspec_cm
!print*, myrank, 'DEBUG offset_nglob_cm', offset_nglob_cm
!print*, myrank, 'DEBUG ibool_crust_mantle', shape(ibool_crust_mantle)
! access memory ibool_crust_mantle
!print*, myrank, 'DEBUG ibool_crust_mantle(1)', ibool_crust_mantle(1,1,1,1)
    call write_array3dspec_as_1d_hdf5('rhonotprime_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    rhonotprime_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.) !!!!!!!!!!!
    call write_array3dspec_as_1d_hdf5('kappa_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    kappa_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
    call write_array3dspec_as_1d_hdf5('mu_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    mu_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)

    call write_array3dspec_as_1d_hdf5('rho_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    rho_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
    call write_array3dspec_as_1d_hdf5('alpha_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    alpha_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
    call write_array3dspec_as_1d_hdf5('beta_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    beta_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)

    call write_array3dspec_as_1d_hdf5('bulk_c_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    bulk_c_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
    call write_array3dspec_as_1d_hdf5('bulk_beta_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    bulk_beta_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)

  endif

  call synchronize_all()

  ! close hdf5
  call h5_close_file_p()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = size(mu_kl_crust_mantle,kind=4)
  idummy = size(kappa_kl_crust_mantle,kind=4)
  idummy = size(rhonotprime_kl_crust_mantle,kind=4)
  idummy = size(bulk_c_kl_crust_mantle,kind=4)
  idummy = size(bulk_beta_kl_crust_mantle,kind=4)

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_cm_iso_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_oc_hdf5()

  use specfem_par
  use specfem_par_outercore

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

#ifdef USE_HDF5

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    if (.not. HDF5_KERNEL_VIS) then
      call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_oc)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_oc)/), 4, CUSTOM_REAL)
    else
      call h5_create_dataset_gen('rho_kernel', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('alpha_kernel', (/sum(offset_nglob_oc)/), 1, CUSTOM_REAL)
    endif

    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! write data
  if (.not. HDF5_KERNEL_VIS) then
    call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_outer_core, &
                                            (/0,0,0,sum(offset_nspec_oc(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_outer_core, &
                                            (/0,0,0,sum(offset_nspec_oc(0:myrank-1))/), H5_COL)
  else
    call write_array3dspec_as_1d_hdf5('rho_kernel', offset_nspec_oc(myrank-1), offset_nglob_oc(myrank-1), &
                                    rho_kl_outer_core, sum(offset_nglob_oc(0:myrank-1)), &
                                    ibool_outer_core, .false.)
    call write_array3dspec_as_1d_hdf5('alpha_kernel', offset_nspec_oc(myrank-1), offset_nglob_oc(myrank-1), &
                                    alpha_kl_outer_core, sum(offset_nglob_oc(0:myrank-1)), &
                                    ibool_outer_core, .false.)
  endif

  call synchronize_all()

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_oc_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_ic_hdf5()

  use specfem_par
  use specfem_par_innercore

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

#ifdef USE_HDF5


  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)

    ! create dataset
    if (.not. HDF5_KERNEL_VIS) then
      call h5_create_dataset_gen('rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_ic)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('alpha_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_ic)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('beta_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_ic)/), 4, CUSTOM_REAL)
    else
      call h5_create_dataset_gen('rho_kernel', (/sum(offset_nglob_ic)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('alpha_kernel', (/sum(offset_nglob_ic)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('beta_kernel', (/sum(offset_nglob_ic)/), 1, CUSTOM_REAL)
    endif
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! write data
  if (.not. HDF5_KERNEL_VIS) then
    call h5_write_dataset_collect_hyperslab('rho_kernel', rho_kl_inner_core, &
                                            (/0,0,0,sum(offset_nspec_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('alpha_kernel', alpha_kl_inner_core, &
                                            (/0,0,0,sum(offset_nspec_ic(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('beta_kernel', beta_kl_inner_core, &
                                            (/0,0,0,sum(offset_nspec_ic(0:myrank-1))/), H5_COL)
  else
    call write_array3dspec_as_1d_hdf5('rho_kernel', offset_nspec_ic(myrank-1), offset_nglob_ic(myrank-1), &
                                    rho_kl_inner_core, sum(offset_nglob_ic(0:myrank-1)), &
                                    ibool_inner_core, .false.)
    call write_array3dspec_as_1d_hdf5('alpha_kernel', offset_nspec_ic(myrank-1), offset_nglob_ic(myrank-1), &
                                    alpha_kl_inner_core, sum(offset_nglob_ic(0:myrank-1)), &
                                    ibool_inner_core, .false.)
    call write_array3dspec_as_1d_hdf5('beta_kernel', offset_nspec_ic(myrank-1), offset_nglob_ic(myrank-1), &
                                    beta_kl_inner_core, sum(offset_nglob_ic(0:myrank-1)), &
                                    ibool_inner_core, .false.)
  endif

  call synchronize_all()

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_ic_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_boundary_kl_hdf5()

  ! TODO: HDF5 kernel vis for boundary kernels not implemented yet

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

#ifdef USE_HDF5
  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! offset array
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_moho
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_400
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_670
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_cmb
  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec2d_icb

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  integer :: info, comm

  if (.not. SAVE_KERNELS_BOUNDARY) return

  ! gather the number of elements in each region
  call gather_all_all_singlei(NSPEC2D_MOHO, offset_nspec2d_moho, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_400, offset_nspec2d_400, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_670, offset_nspec2d_670, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_CMB, offset_nspec2d_cmb, NPROCTOT_VAL)
  call gather_all_all_singlei(NSPEC2D_ICB, offset_nspec2d_icb, NPROCTOT_VAL)

  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)
    ! create dataset
    if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
      call h5_create_dataset_gen('moho_kernel', (/0,0,sum(offset_nspec2d_moho(0:myrank-1))/), 3, CUSTOM_REAL)
    endif
    call h5_create_dataset_gen('d400_kernel', (/0,0,sum(offset_nspec2d_400(0:myrank-1))/), 3, CUSTOM_REAL)
    call h5_create_dataset_gen('d670_kernel', (/0,0,sum(offset_nspec2d_670(0:myrank-1))/), 3, CUSTOM_REAL)
    call h5_create_dataset_gen('CMB_kernel',  (/0,0,sum(offset_nspec2d_cmb(0:myrank-1))/), 3, CUSTOM_REAL)
    call h5_create_dataset_gen('ICB_kernel',  (/0,0,sum(offset_nspec2d_icb(0:myrank-1))/), 3, CUSTOM_REAL)

    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! write data
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    call h5_write_dataset_collect_hyperslab('moho_kernel', moho_kl, (/0,0,sum(offset_nspec2d_moho(0:myrank-1))/), H5_COL)
  endif
  call h5_write_dataset_collect_hyperslab('d400_kernel', d400_kl, (/0,0,sum(offset_nspec2d_400(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('d670_kernel', d670_kl, (/0,0,sum(offset_nspec2d_670(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('CMB_kernel', CMB_kl, (/0,0,sum(offset_nspec2d_CMB(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab('ICB_kernel', ICB_kl, (/0,0,sum(offset_nspec2d_ICB(0:myrank-1))/), H5_COL)

  call synchronize_all()

  ! close hdf5
  call h5_close_file_p()

#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_boundary_kl_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_kernels_Hessian_hdf5()

  use specfem_par
  use specfem_par_crustmantle

#ifdef USE_HDF5
  use manager_hdf5
  use specfem_par_movie_hdf5
#endif

  implicit none

#ifdef USE_HDF5


  ! initialize hdf5
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize() ! called in initialize_mesher()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'

  if (myrank == 0) then
    ! check if file exists
    call h5_create_or_open_file(file_name)

    ! create dataset
    if (.not. HDF5_KERNEL_VIS) then
      call h5_create_dataset_gen('hess_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('hess_rho_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('hess_kappa_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
      call h5_create_dataset_gen('hess_mu_kernel', (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec_cm)/), 4, CUSTOM_REAL)
    else ! HDF5_KERNEL_VIS
      call h5_create_dataset_gen('hess_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('hess_rho_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('hess_kappa_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen('hess_mu_kernel', (/sum(offset_nglob_cm)/), 1, CUSTOM_REAL)
    endif
    ! close file
    call h5_close_file()
  endif

  call synchronize_all()

  ! open hdf5
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  ! write data
  if (.not. HDF5_KERNEL_VIS) then
    call h5_write_dataset_collect_hyperslab('hess_kernel', hess_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('hess_rho_kernel', hess_rho_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('hess_kappa_kernel', hess_kappa_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab('hess_mu_kernel', hess_mu_kl_crust_mantle, &
                                            (/0,0,0,sum(offset_nspec_cm(0:myrank-1))/), H5_COL)
  else
    call write_array3dspec_as_1d_hdf5('hess_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    hess_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
    call write_array3dspec_as_1d_hdf5('hess_rho_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    hess_rho_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
    call write_array3dspec_as_1d_hdf5('hess_kappa_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    hess_kappa_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
    call write_array3dspec_as_1d_hdf5('hess_mu_kernel', offset_nspec_cm(myrank-1), offset_nglob_cm(myrank-1), &
                                    hess_mu_kl_crust_mantle, sum(offset_nglob_cm(0:myrank-1)), &
                                    ibool_crust_mantle, .false.)
  endif

  call synchronize_all()

  ! close hdf5
  call h5_close_file_p()


#else

  print *,'Error: HDF5 not enabled in this version of the code'
  print *, 'Please recompile with the HDF5 option enabled with the configure flag --with-hdf5'
  call exit_mpi(myrank,'Error: HDF5 not enabled in this version of the code')

#endif

  end subroutine write_kernels_Hessian_hdf5

!
!-------------------------------------------------------------------------------------------------
!
  subroutine write_mesh_info_kernel_h5()

  use specfem_par

#ifdef USE_HDF5
  use specfem_par_crustmantle, only: ibool_crust_mantle,rstore_crust_mantle
  use specfem_par_outercore, only: ibool_outer_core,rstore_outer_core
  use specfem_par_innercore, only: ibool_inner_core,rstore_inner_core
  use specfem_par_movie_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! local parameters
  integer :: ispec,i,j,k
  integer :: iglob
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE)  :: store_val3D_x_cm, store_val3D_y_cm, store_val3D_z_cm
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE)    :: store_val3D_x_oc, store_val3D_y_oc, store_val3D_z_oc
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE)    :: store_val3D_x_ic, store_val3D_y_ic, store_val3D_z_ic
  ! dummy num_ibool_3dmovie for cm oc ic
  integer, dimension(NGLOB_CRUST_MANTLE) :: num_ibool_3dmovie_cm
  integer, dimension(NGLOB_OUTER_CORE)   :: num_ibool_3dmovie_oc
  integer, dimension(NGLOB_INNER_CORE)   :: num_ibool_3dmovie_ic
  ! dummy mask_ibool_3dmovie for cm oc ic
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool_3dmovie_cm
  logical, dimension(NGLOB_OUTER_CORE)   :: mask_ibool_3dmovie_oc
  logical, dimension(NGLOB_INNER_CORE)   :: mask_ibool_3dmovie_ic

  integer, dimension(:,:), allocatable :: elm_conn_cm, elm_conn_oc, elm_conn_ic

  integer, dimension(0:NPROCTOT_VAL-1) :: offset_nspec_cm_conn, offset_nspec_oc_conn, offset_nspec_ic_conn

  integer :: nelems_3dmovie_cm, nelems_3dmovie_oc, nelems_3dmovie_ic


  ! safety check
  if (NDIM /= 3) stop 'movie volume output requires NDIM = 3'

  ! outer core and inner core is always fine mesh
  nelems_3dmovie_cm = NSPEC_CRUST_MANTLE * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
  nelems_3dmovie_oc = NSPEC_OUTER_CORE   * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
  nelems_3dmovie_ic = NSPEC_INNER_CORE   * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)

  ! allocate elm_conn
  allocate(elm_conn_cm(9,nelems_3dmovie_cm))
  allocate(elm_conn_oc(9,nelems_3dmovie_oc))
  allocate(elm_conn_ic(9,nelems_3dmovie_ic))

  ! offset arrays for element connectivity
  call gather_all_all_singlei(nelems_3dmovie_cm, offset_nspec_cm_conn, NPROCTOT_VAL)
  call gather_all_all_singlei(nelems_3dmovie_oc, offset_nspec_oc_conn, NPROCTOT_VAL)
  call gather_all_all_singlei(nelems_3dmovie_ic, offset_nspec_ic_conn, NPROCTOT_VAL)

  nspec_vol_mov_all_proc_cm_conn = sum(offset_nspec_cm_conn)
  nspec_vol_mov_all_proc_oc_conn = sum(offset_nspec_oc_conn)
  nspec_vol_mov_all_proc_ic_conn = sum(offset_nspec_ic_conn)

  !
  ! create the xyz arrays for crust and mantle (not strain or vector output)
  !

  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,1
      do j = 1,NGLLY,1
        do i = 1,NGLLX,1
          iglob           = ibool_crust_mantle(i,j,k,ispec)
          rval            = rstore_crust_mantle(1,iglob)
          thetaval        = rstore_crust_mantle(2,iglob)
          phival          = rstore_crust_mantle(3,iglob)

          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

          store_val3D_x_cm(iglob) = xval
          store_val3D_y_cm(iglob) = yval
          store_val3D_z_cm(iglob) = zval

          ! dummy num_ibool_3dmovie for cm
          num_ibool_3dmovie_cm(iglob) = iglob
          ! all the mask_ibool_3dmovie_cm are true
          mask_ibool_3dmovie_cm(iglob) = .true.
        enddo
      enddo
    enddo
  enddo

  !
  ! create the xyz arrays for outer core
  !

  do ispec = 1, NSPEC_OUTER_CORE
    do k = 1,NGLLZ,1
      do j = 1,NGLLY,1
        do i = 1,NGLLX,1
          iglob           = ibool_outer_core(i,j,k,ispec)
          rval            = rstore_outer_core(1,iglob)
          thetaval        = rstore_outer_core(2,iglob)
          phival          = rstore_outer_core(3,iglob)

          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

          store_val3D_x_oc(iglob) = xval
          store_val3D_y_oc(iglob) = yval
          store_val3D_z_oc(iglob) = zval

          ! dummy num_ibool_3dmovie for oc
          num_ibool_3dmovie_oc(iglob) = iglob
          ! all the mask_ibool_3dmovie_oc are true
          mask_ibool_3dmovie_oc(iglob) = .true.
        enddo
      enddo
    enddo
  enddo

  !
  ! create the xyz arrays for inner core
  !
  do ispec = 1, NSPEC_INNER_CORE
    do k = 1,NGLLZ,1
      do j = 1,NGLLY,1
        do i = 1,NGLLX,1
          iglob           = ibool_inner_core(i,j,k,ispec)
          rval            = rstore_inner_core(1,iglob)
          thetaval        = rstore_inner_core(2,iglob)
          phival          = rstore_inner_core(3,iglob)

          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

          store_val3D_x_ic(iglob) = xval
          store_val3D_y_ic(iglob) = yval
          store_val3D_z_ic(iglob) = zval

          ! dummy num_ibool_3dmovie for ic
          num_ibool_3dmovie_ic(iglob) = iglob
          ! all the mask_ibool_3dmovie_ic are true
          mask_ibool_3dmovie_ic(iglob) = .true.
        enddo
      enddo
    enddo
  enddo

  ! create elm_conn
  call get_conn_for_movie(elm_conn_cm, sum(offset_nglob_cm(0:myrank-1)), 1, nelems_3dmovie_cm, &
                        NGLOB_CRUST_MANTLE, NSPEC_CRUST_MANTLE, num_ibool_3dmovie_cm, mask_ibool_3dmovie_cm, ibool_crust_mantle)
  ! for outer core
  call get_conn_for_movie(elm_conn_oc, sum(offset_nglob_oc(0:myrank-1)), 1, nelems_3dmovie_oc, &
                        NGLOB_OUTER_CORE, NSPEC_OUTER_CORE, num_ibool_3dmovie_oc, mask_ibool_3dmovie_oc, ibool_outer_core)
  ! for inner core
  call get_conn_for_movie(elm_conn_ic, sum(offset_nglob_ic(0:myrank-1)), 1, nelems_3dmovie_ic, &
                        NGLOB_INNER_CORE, NSPEC_INNER_CORE, num_ibool_3dmovie_ic, mask_ibool_3dmovie_ic, ibool_inner_core)

  ! initialize h5 file for volume movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! file name and group name
  file_name = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))//'/kernels.h5'
  group_name = 'mesh'

  ! create the file, group and dataset
  if (myrank == 0) then
    call h5_create_file(file_name)
    call h5_open_or_create_group(group_name)

    ! create the dataset
    ! for crust and mantle
    call h5_create_dataset_gen_in_group('elm_conn_cm', (/9, nspec_vol_mov_all_proc_cm_conn/), 2, 1)
    call h5_create_dataset_gen_in_group('x_cm', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('y_cm', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('z_cm', (/npoints_vol_mov_all_proc_cm/), 1, CUSTOM_REAL)
    ! for outer core
    call h5_create_dataset_gen_in_group('elm_conn_oc', (/9, nspec_vol_mov_all_proc_oc_conn/), 2, 1)
    call h5_create_dataset_gen_in_group('x_oc', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('y_oc', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('z_oc', (/npoints_vol_mov_all_proc_oc/), 1, CUSTOM_REAL)
    ! for inner core
    call h5_create_dataset_gen_in_group('elm_conn_ic', (/9, nspec_vol_mov_all_proc_ic_conn/), 2, 1)
    call h5_create_dataset_gen_in_group('x_ic', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('y_ic', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group('z_ic', (/npoints_vol_mov_all_proc_ic/), 1, CUSTOM_REAL)

    ! close the group
    call h5_close_group()
    ! close the file
    call h5_close_file()

  endif

  ! synchronize
  call synchronize_all()

  ! write the data
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif

  call h5_open_group(group_name)

  ! write the data
  ! for crust and mantle
  call h5_write_dataset_collect_hyperslab_in_group('elm_conn_cm', elm_conn_cm, &
                                                   (/0, sum(offset_nspec_cm_conn(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('x_cm', store_val3D_x_cm, (/sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('y_cm', store_val3D_y_cm, (/sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('z_cm', store_val3D_z_cm, (/sum(offset_nglob_cm(0:myrank-1))/), H5_COL)
  ! for outer core
  call h5_write_dataset_collect_hyperslab_in_group('elm_conn_oc', elm_conn_oc, &
                                                   (/0, sum(offset_nspec_oc_conn(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('x_oc', store_val3D_x_oc, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('y_oc', store_val3D_y_oc, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('z_oc', store_val3D_z_oc, (/sum(offset_nglob_oc(0:myrank-1))/), H5_COL)
  ! for inner core
  call h5_write_dataset_collect_hyperslab_in_group('elm_conn_ic', elm_conn_ic, &
                                                   (/0, sum(offset_nspec_ic_conn(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('x_ic', store_val3D_x_ic, (/sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('y_ic', store_val3D_y_ic, (/sum(offset_nglob_ic(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group('z_ic', store_val3D_z_ic, (/sum(offset_nglob_ic(0:myrank-1))/), H5_COL)

  call synchronize_all()

  ! close the group
  call h5_close_group()
  ! close the file
  call h5_close_file_p()

  ! deallocate
  deallocate(elm_conn_cm)
  deallocate(elm_conn_oc)
  deallocate(elm_conn_ic)

  ! write xdmf for all timesteps
  call write_xdmf_kernel_hdf5()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = nelems_3dmovie_in

  idummy = size(num_ibool_3dmovie,kind=4)
  idummy = size(nu_3dmovie,kind=4)
  idummy = size(mask_3dmovie,kind=4)
  idummy = size(mask_ibool_3dmovie,kind=4)
  idummy = size(muvstore_crust_mantle_3dmovie,kind=4)

  print *, 'Error: HDF5 is not enabled in this version of the code.'
  print *, 'Please recompile with the HDF5 option enabled with --with-hdf5'
  stop

#endif


  end subroutine write_mesh_info_kernel_h5

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine write_xdmf_kernel_hdf5_header(nelems, nglobs, fname_h5_data_xdmf, target_unit, region_flag)

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    integer, intent(in) :: nelems, nglobs
    character(len=*), intent(in) :: fname_h5_data_xdmf
    integer, intent(in) :: target_unit
    integer, intent(in) :: region_flag ! 1: crust mantle, 2: outer core, 3: inner core

    character(len=64) :: nelm_str, nglo_str, elemconn_str, x_str, y_str, z_str

   if (region_flag == 1) then
      elemconn_str = 'elm_conn_cm'
      x_str = 'x_cm'
      y_str = 'y_cm'
      z_str = 'z_cm'
    else if (region_flag == 2) then
      elemconn_str = 'elm_conn_oc'
      x_str = 'x_oc'
      y_str = 'y_oc'
      z_str = 'z_oc'
    else if (region_flag == 3) then
      elemconn_str = 'elm_conn_ic'
      x_str = 'x_ic'
      y_str = 'y_ic'
      z_str = 'z_ic'
    else
      print *,'Error: invalid region_flag in write_xdmf_hdf5_header'
      call exit_mpi(myrank,'Error invalid region_flag in write_xdmf_hdf5_header')
    endif

    ! convert integer to string
    nelm_str = i2c(nelems)
    nglo_str = i2c(nglobs)

    ! definition of topology and geometry
    ! refer only control nodes (8 or 27) as a coarse output
    ! data array need to be extracted from full data array on GLL points
    write(target_unit,'(a)') '<?xml version="1.0" ?>'
    write(target_unit,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(target_unit,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
    write(target_unit,*) '<Domain name="mesh">'

    ! loop for writing information of mesh partitions
    write(target_unit,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm_str)//'">'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                         //trim(nelm_str)//' 9">'
    write(target_unit,*) '       '//trim(fname_h5_data_xdmf)//':/mesh/'//trim(elemconn_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '</Topology>'
    write(target_unit,*) '<Geometry GeometryType="X_Y_Z">'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
    write(target_unit,*) '       '//trim(fname_h5_data_xdmf)//':/mesh/'//trim(x_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
    write(target_unit,*) '       '//trim(fname_h5_data_xdmf)//':/mesh/'//trim(y_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo_str)//'">'
    write(target_unit,*) '       '//trim(fname_h5_data_xdmf)//':/mesh/'//trim(z_str)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '</Geometry>'

  end subroutine write_xdmf_kernel_hdf5_header

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_xdmf_kernel_hdf5_footer(target_unit)

    implicit none

    integer, intent(in) :: target_unit

    ! file finish
    write(target_unit,*) '</Domain>'
    write(target_unit,*) '</Xdmf>'

  end subroutine write_xdmf_kernel_hdf5_footer

  subroutine write_xdmf_kernel_hdf5_one_data(fname_h5, attr_name, dset_name, len_data, target_unit, value_on_node)

    use specfem_par
    use specfem_par_movie_hdf5

    implicit none

    character(len=*), intent(in) :: fname_h5
    character(len=*), intent(in) :: attr_name
    character(len=*), intent(in) :: dset_name
    integer, intent(in) :: len_data
    integer, intent(in) :: target_unit
    logical, intent(in) :: value_on_node ! values are defined on nodes or elements

    character(len=20) :: center_str

    if (value_on_node) then
      center_str = 'Node' ! len data should be the number of nodes
    else
      center_str = 'Cell' ! len data should be the number of elements
    endif

    write(target_unit,*) '<Grid Name="kernel" GridType="Uniform">'
    write(target_unit,*) '<Topology Reference="/Xdmf/Domain/Topology" />'
    write(target_unit,*) '<Geometry Reference="/Xdmf/Domain/Geometry" />'

    ! write a header for single data
    write(target_unit,*) '<Attribute Name="'//trim(attr_name)//'" AttributeType="Scalar" Center="'//trim(center_str)//'">'
    write(target_unit,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(len_data))//'">'
    write(target_unit,*) '      '//trim(fname_h5)//':/'//trim(dset_name)
    write(target_unit,*) '</DataItem>'
    write(target_unit,*) '</Attribute>'

    write(target_unit,*) '</Grid>'
  end subroutine write_xdmf_kernel_hdf5_one_data


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_xdmf_kernel_hdf5

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  ! local parameters
  integer                       :: ierr
  character(len=MAX_STRING_LEN) :: fname_xdmf_kl, fname_xdmf_kl_oc, fname_xdmf_kl_ic
  character(len=MAX_STRING_LEN) :: fname_h5_data_kl_xdmf

  ! checks if anything do, only main process writes out xdmf file
  if (myrank /= 0) return

  !
  ! write out the crust and mantle xdmf file
  !
  fname_xdmf_kl = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH)) // "/kernel_cm.xmf"
  fname_h5_data_kl_xdmf = "./kernels.h5"  ! relative to movie_volume_cm.xmf file

  ! open xdmf file
  open(unit=xdmf_kl, file=trim(fname_xdmf_kl), status='replace', action='write', iostat=ierr, recl=256)
  if (ierr /= 0) then
    print *, 'Error: could not open kernel XDMF file ', trim(fname_xdmf_kl)
    return
  endif

  call write_xdmf_kernel_hdf5_header(nspec_vol_mov_all_proc_cm_conn, npoints_vol_mov_all_proc_cm, &
                                    fname_h5_data_kl_xdmf, xdmf_kl, 1)

  ! write headers for each dataset
  if (ANISOTROPIC_KL) then
    ! anisotropic kernels
    if (SAVE_TRANSVERSE_KL_ONLY ) then
      ! alphav_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alphav_kernel', 'alphav_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! alphah_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alphah_kernel', 'alphah_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! betav_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'betav_kernel', 'betav_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! betah_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'betah_kernel', 'betah_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! eta_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'eta_kernel', 'eta_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! rho_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'rho_kernel', 'rho_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! bulk_c_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_c_kernel', 'bulk_c_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! bulk_betav_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_betav_kernel', 'bulk_betav_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node

      ! bulk_betah_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_betah_kernel', 'bulk_betah_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! alpha_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alpha_kernel', 'alpha_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! beta_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'beta_kernel', 'beta_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! bulk_beta_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_beta_kernel', 'bulk_beta_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node

    else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
      ! alphav_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alphav_kernel', 'alphav_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! alphah_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alphah_kernel', 'alphah_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! betav_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'betav_kernel', 'betav_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! betah_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'betah_kernel', 'betah_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! eta_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'eta_kernel', 'eta_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! rho_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'rho_kernel', 'rho_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! bulk_c_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_c_kernel', 'bulk_c_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! bulk_betav_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_betav_kernel', 'bulk_betav_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
      ! bulk_betah_kernel
      call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_betah_kernel', 'bulk_betah_kernel', &
                                        npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    else
      ! full anisotropic kernels
      ! not implemented yet
      call exit_MPI(myrank,'Error: full Anisotropic kernels not implemented for HDF5 yet')
      ! rho_kernel
      ! cijkl_kernel

    endif
  else
    ! isotropic kernels
    ! rhonotprime_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'rhonotprime_kernel', 'rhonotprime_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! kappa_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'kappa_kernel', 'kappa_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! mu_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'mu_kernel', 'mu_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! rho_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'rho_kernel', 'rho_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! alpha_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alpha_kernel', 'alpha_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! beta_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'beta_kernel', 'beta_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! bulk_c_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_c_kernel', 'bulk_c_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! bulk_beta_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'bulk_beta_kernel', 'bulk_beta_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
  endif

  if (SAVE_REGULAR_KL) then
    ! not implemented yet
    call exit_MPI(myrank,'Error: SAVE_REGULAR_KL not implemented for HDF5 yet')
  endif

  if (FULL_GRAVITY_VAL) then
    ! not implemented yet
    call exit_MPI(myrank,'Error: FULL_GRAVITY_VAL not implemented for HDF5 yet')
  endif

  if (NOISE_TOMOGRAPHY == 3) then
    ! sigma_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'sigma_kernel', 'sigma_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
  endif

  if (APPROXIMATE_HESS_KL) then
    ! hess_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'hess_kernel', 'hess_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! hess_rho_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'hess_rho_kernel', 'hess_rho_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! hess_kappa_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'hess_kappa_kernel', 'hess_kappa_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
    ! hess_mu_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'hess_mu_kernel', 'hess_mu_kernel', &
                                      npoints_vol_mov_all_proc_cm, xdmf_kl, .true.) ! value on node
  endif

  call write_xdmf_kernel_hdf5_footer(xdmf_kl)

  ! close xdmf file
  close(xdmf_kl, iostat=ierr)

  !
  ! write out the outer core xdmf file
  !
  if (SAVE_KERNELS_OC) then

    fname_xdmf_kl_oc = LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH)) // "/kernel_oc.xmf"
    fname_h5_data_kl_xdmf = "./kernels.h5"  ! relative to movie_volume_oc.xmf file

    ! open xdmf file
    open(unit=xdmf_kl, file=trim(fname_xdmf_kl_oc), status='replace', action='write', iostat=ierr, recl=256)
    if (ierr /= 0) then
      print *, 'Error: could not open kernel XDMF file ', trim(fname_xdmf_kl_oc)
      return
    endif

    call write_xdmf_kernel_hdf5_header(nspec_vol_mov_all_proc_oc_conn, npoints_vol_mov_all_proc_oc, &
                                       fname_h5_data_kl_xdmf, xdmf_kl, 2)

    ! write headers for each dataset
    ! rho_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'rho_kernel', 'rho_kernel', &
                                      npoints_vol_mov_all_proc_oc, xdmf_kl, .true.) ! value on node
    ! alpha_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alpha_kernel', 'alpha_kernel', &
                                      npoints_vol_mov_all_proc_oc, xdmf_kl, .true.) ! value on node

    call write_xdmf_kernel_hdf5_footer(xdmf_kl)

    ! close xdmf file
    close(xdmf_kl, iostat=ierr)

  endif


  !
  ! write out the inner core xdmf file
  !
  if (SAVE_KERNELS_IC) then
    fname_xdmf_kl_ic = trim(OUTPUT_FILES) // "/kernel_ic.xmf"
    fname_h5_data_kl_xdmf = "./kernels.h5"  ! relative to movie_volume_ic.xmf file

    ! open xdmf file
    open(unit=xdmf_kl, file=trim(fname_xdmf_kl_ic), status='replace', action='write', iostat=ierr, recl=256)
    if (ierr /= 0) then
      print *, 'Error: could not open kernel XDMF file ', trim(fname_xdmf_kl_ic)
      return
    endif

    call write_xdmf_kernel_hdf5_header(nspec_vol_mov_all_proc_ic_conn, npoints_vol_mov_all_proc_ic, &
                                    fname_h5_data_kl_xdmf, xdmf_kl, 3)

    ! write headers for each dataset
    ! rho_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'rho_kernel', 'rho_kernel', &
                                      npoints_vol_mov_all_proc_ic, xdmf_kl, .true.) ! value on node
    ! alpha_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'alpha_kernel', 'alpha_kernel', &
                                      npoints_vol_mov_all_proc_ic, xdmf_kl, .true.) ! value on node
    ! beta_kernel
    call write_xdmf_kernel_hdf5_one_data(fname_h5_data_kl_xdmf, 'beta_kernel', 'beta_kernel', &
                                      npoints_vol_mov_all_proc_ic, xdmf_kl, .true.) ! value on node

    call write_xdmf_kernel_hdf5_footer(xdmf_kl)

    ! close xdmf file
    close(xdmf_kl, iostat=ierr)

  endif


  ! other kernels not implemented yet
  if (SAVE_KERNELS_BOUNDARY) then
    ! not implemented yet
    call exit_MPI(myrank,'Error: SAVE_KERNELS_BOUNDARY not implemented for HDF5 yet')
  endif

  end subroutine write_xdmf_kernel_hdf5

#endif
