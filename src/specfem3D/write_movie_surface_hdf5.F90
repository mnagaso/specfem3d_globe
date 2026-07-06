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

subroutine movie_surface_init_hdf5()
#ifdef USE_HDF5
  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer :: ier

  allocate(offset_poin(0:NPROCTOT_VAL-1),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating offset_poin array')

  npoints_surf_mov_all_proc = 0

#else

      write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
      write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
      stop
#endif
end subroutine movie_surface_init_hdf5


subroutine movie_surface_finalize_hdf5()
#ifdef USE_HDF5
  use specfem_par_movie_hdf5

  implicit none

  deallocate(offset_poin)

#else

      write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
      write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
      stop
#endif
end subroutine movie_surface_finalize_hdf5

subroutine write_movie_surface_mesh_hdf5()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

#ifdef USE_HDF5
  use specfem_par_movie_hdf5

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_x,store_val_y,store_val_z
  integer :: ipoin,ispec2D,ispec,i,j,k,ier,iglob1,iglob2,iglob3,iglob4,npoin
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval

  call movie_surface_init_hdf5()

  ! gather npoints on each process
  call gather_all_all_singlei(nmovie_points,offset_poin,NPROCTOT_VAL)
  ! total number of points on all processes
  npoints_surf_mov_all_proc = sum(offset_poin)

  ! allocates movie surface arrays
  allocate(store_val_x(nmovie_points), &
           store_val_y(nmovie_points), &
           store_val_z(nmovie_points),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface location arrays')

  ! initialize h5 file for surface movie
  call world_get_comm(comm)
  call world_get_info_null(info)
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

  ! create file and dataset
  file_name = trim(OUTPUT_FILES)//"/movie_surface.h5"

  if (myrank == 0) then
    call h5_create_file(file_name)
    ! create group surf_coord
    call h5_create_group("surf_coord")
    call h5_open_group("surf_coord")
    ! create datasets x, y, z
    call h5_create_dataset_gen_in_group("x", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group("y", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
    call h5_create_dataset_gen_in_group("z", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)

    ! close group surf_coord
    call h5_close_group()

    ! close file
    call h5_close_file()
  endif

  ! gets coordinates of surface mesh
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)
    ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
    ! the cubed sphere, the mapping changes for different chunks,
    ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
    ! for future consideration, like in create_movie_GMT_global.f90 ...
    k = NGLLZ
    ! loop on all the points inside the element
    if (.not. MOVIE_COARSE) then
      do j = 1, NGLLY-1, 1
        do i = 1, NGLLX-1, 1
          ! stores values
          iglob1 = ibool_crust_mantle(i,j,k,ispec)
          iglob2 = ibool_crust_mantle(i+1,j,k,ispec)
          iglob3 = ibool_crust_mantle(i+1,j+1,k,ispec)
          iglob4 = ibool_crust_mantle(i,j+1,k,ispec)
          ! iglob1
          ipoin    = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob1) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob1) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob1) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
          ! iglob2
          ipoin = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob2) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob2) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob2) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
          ! iglob3
          ipoin = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob3) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob3) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob3) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
          ! iglob4
          ipoin = ipoin + 1
          rval     = rstore_crust_mantle(1,iglob4) ! radius r (normalized)
          thetaval = rstore_crust_mantle(2,iglob4) ! colatitude theta (in radian)
          phival   = rstore_crust_mantle(3,iglob4) ! longitude phi (in radian)
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val_x(ipoin) = xval
          store_val_y(ipoin) = yval
          store_val_z(ipoin) = zval
        enddo
      enddo
    else ! MOVIE_COARSE
      iglob1 = ibool_crust_mantle(1,1,k,ispec)
      iglob2 = ibool_crust_mantle(NGLLX,1,k,ispec)
      iglob3 = ibool_crust_mantle(NGLLX,NGLLY,k,ispec)
      iglob4 = ibool_crust_mantle(1,NGLLY,k,ispec)

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob1) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob1) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob1) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob2) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob2) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob2) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob3) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob3) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob3) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

      ipoin = ipoin + 1
      rval  = rstore_crust_mantle(1,iglob4) ! radius r (normalized)
      thetaval = rstore_crust_mantle(2,iglob4) ! colatitude theta (in radian)
      phival   = rstore_crust_mantle(3,iglob4) ! longitude phi (in radian)
      call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
      store_val_x(ipoin) = xval
      store_val_y(ipoin) = yval
      store_val_z(ipoin) = zval

    endif
  enddo
  npoin = ipoin
  if (npoin /= nmovie_points ) call exit_mpi(myrank,'Error number of movie points not equal to nmovie_points')

  call synchronize_all()

  ! write data to h5 file
  if (H5_COL) then
    ! open file
    call h5_open_file_p_collect(file_name)
  else
    ! open file
    call h5_open_file_p(file_name)
  endif
  call h5_open_group("surf_coord")

  ! write x, y, z
  call h5_write_dataset_collect_hyperslab_in_group("x", store_val_x, (/sum(offset_poin(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group("y", store_val_y, (/sum(offset_poin(0:myrank-1))/), H5_COL)
  call h5_write_dataset_collect_hyperslab_in_group("z", store_val_z, (/sum(offset_poin(0:myrank-1))/), H5_COL)

  ! close group and file
  call h5_close_group()
  call h5_close_file_p()

  deallocate(store_val_x,store_val_y,store_val_z)

  if (myrank == 0) then
    ! write XDMF description
    if (HDF5_IO_NODES > 1) then
      ! multi-IO mode: one XDMF file per IO shard referencing
      ! geometry in movie_surface.h5 and data in movie_surface.io<k>.h5
      call write_xdmf_surface_shards(npoints_surf_mov_all_proc)
    else
      ! single-file mode: complete XDMF referencing movie_surface.h5
      call write_xdmf_surface_complete(npoints_surf_mov_all_proc)
    endif
  endif

#else

    write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
    write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif

end subroutine write_movie_surface_mesh_hdf5


subroutine write_movie_surface_hdf5()

#ifdef USE_HDF5
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie
  use specfem_par_movie_hdf5
  use io_server_hdf5
  use timing_accounting

  implicit none

  ! local parameters
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob1,iglob2,iglob3,iglob4
  integer :: req_count

  ! by default: save velocity here to avoid static offset on displacement for movies

  ! gets coordinates of surface mesh and surface displacement
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
    ! the cubed sphere, the mapping changes for different chunks,
    ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
    ! for future consideration, like in create_movie_GMT_global.f90 ...
    k = NGLLZ

    ! loop on all the points inside the element
    if (.not. MOVIE_COARSE) then
      do j = 1, NGLLY-1, 1
        do i = 1, NGLLX-1, 1
          ! stores values
          iglob1 = ibool_crust_mantle(i,j,k,ispec)
          iglob2 = ibool_crust_mantle(i+1,j,k,ispec)
          iglob3 = ibool_crust_mantle(i+1,j+1,k,ispec)
          iglob4 = ibool_crust_mantle(i,j+1,k,ispec)

          if (MOVIE_VOLUME_TYPE == 5) then
            ! stores displacement
            ! iglob1
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob2
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob3
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob4
            ipoin = ipoin + 1
            store_val_ux(ipoin) = displ_crust_mantle(1,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = displ_crust_mantle(2,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = displ_crust_mantle(3,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
          else
            ! stores velocity
            ! iglob1
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob2
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob3
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
            ! iglob4
            ipoin = ipoin + 1
            store_val_ux(ipoin) = veloc_crust_mantle(1,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
            store_val_uy(ipoin) = veloc_crust_mantle(2,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
            store_val_uz(ipoin) = veloc_crust_mantle(3,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
          endif
        enddo
      enddo
    else ! MOVIE_COARSE
      iglob1 = ibool_crust_mantle(1,1,k,ispec)
      iglob2 = ibool_crust_mantle(NGLLX-1,1,k,ispec)
      iglob3 = ibool_crust_mantle(NGLLX-1,NGLLY-1,k,ispec)
      iglob4 = ibool_crust_mantle(1,NGLLY-1,k,ispec)

      if (MOVIE_VOLUME_TYPE == 5) then
        ! stores displacement
        ! iglob1
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob1) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob2
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob2) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob3
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob3) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob4
        ipoin = ipoin + 1
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob4) * real(scale_displ,kind=CUSTOM_REAL) ! longitude phi (in radian)
      else
        ! stores velocity
        ! iglob1
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob1) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob2
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob2) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob3
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob3) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
        ! iglob4
        ipoin = ipoin + 1
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! radius r (normalized)
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! colatitude theta (in radian)
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob4) * real(scale_veloc,kind=CUSTOM_REAL) ! longitude phi (in radian)
      endif
    endif
  enddo
  ! use IO server when dedicated HDF5 IO nodes are enabled
  if (HDF5_IO_NODES > 0) then

    ! hdf5 i/o server request index for surface movie
    req_count = 1

    ! debug: log surface movie send on compute side
    print *, 'compute surf_send: it', it, 'rank', myrank, 'dest_ionod', dest_ionod, 'npoints', ipoin

    ! send surface movie data to IO server
    call isend_cr_inter(store_val_ux,ipoin,dest_ionod, &
                        io_tag_surf_ux,req_dump_surf(req_count))
    req_count = req_count + 1
    call isend_cr_inter(store_val_uy,ipoin,dest_ionod, &
                        io_tag_surf_uy,req_dump_surf(req_count))
    req_count = req_count + 1
    call isend_cr_inter(store_val_uz,ipoin,dest_ionod, &
                        io_tag_surf_uz,req_dump_surf(req_count))
    req_count = req_count + 1

    ! store number of MPI_ISEND requests for surface movie
    n_req_surf = req_count - 1

    ! wait for current frame's sends to complete before returning to time loop
    ! (prevents MPI resource conflict between inter-communicator and compute-only communicator)
    call timing_io_stop()
    call timing_wait_start()
    call wait_all_send()
    call timing_wait_stop()
    call timing_io_start()

    ! in multi-IO mode, XDMF is generated once per shard by
    ! write_xdmf_surface_shards() and does not need per-step appends

  else

    ! initialize h5 file for surface movie
    call world_get_comm(comm)
    call world_get_info_null(info)
    call h5_initialize()
    call h5_set_mpi_info(comm, info, myrank, NPROCTOT_VAL)

    ! create file and dataset
    file_name = trim(OUTPUT_FILES)//"/movie_surface.h5"
    group_name = "it_"//trim(i2c(it))

    ! create dataset
    if (myrank == 0) then
      call h5_open_file(file_name)
      call h5_create_group(group_name)
      call h5_open_group(group_name)

      ! create datasets ux, uy, uz
      call h5_create_dataset_gen_in_group("ux", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uy", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)
      call h5_create_dataset_gen_in_group("uz", (/npoints_surf_mov_all_proc/), 1, CUSTOM_REAL)

      ! close group
      call h5_close_group()
      ! close file
      call h5_close_file()
    endif

    call synchronize_all()

    ! write data to h5 file
    if (H5_COL) then
      ! open file
      call h5_open_file_p_collect(file_name)
    else
      ! open file
      call h5_open_file_p(file_name)
    endif

    call h5_open_group(group_name)

    ! write ux, uy, uz
    call h5_write_dataset_collect_hyperslab_in_group("ux", store_val_ux, (/sum(offset_poin(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group("uy", store_val_uy, (/sum(offset_poin(0:myrank-1))/), H5_COL)
    call h5_write_dataset_collect_hyperslab_in_group("uz", store_val_uz, (/sum(offset_poin(0:myrank-1))/), H5_COL)

    ! close group and file
    call h5_close_group()
    call h5_close_file_p()

    ! Note: XDMF file is written once at initialization (complete file approach)
    ! so no per-timestep write_xdmf_surface_body call is needed

  endif

#else

    write(*,*) 'Error: HDF5 is not enabled in this version of Specfem3D_Globe.'
    write(*,*) 'Please recompile with the HDF5 option enabled with --with-hdf5'
    stop

#endif


end subroutine write_movie_surface_hdf5



!
! xdmf output routines
!
! The XDMF routines generate XML metadata files that describe the HDF5 data layout
! for visualization software (ParaView, VisIt, etc.).
!
! For surface movies, there are two modes:
!   - Single-file mode (HDF5_IO_NODES <= 1): writes movie_surface.xmf referencing movie_surface.h5
!   - Multi-IO mode (HDF5_IO_NODES > 1): writes movie_surface_io<k>.xmf for each IO shard k
!
! Both modes now use complete file generation (not incremental appending) for robustness.
!
#ifdef USE_HDF5

  subroutine write_xdmf_surface_complete(num_nodes)
  !
  ! Writes a complete XDMF file for surface movie in single-file mode.
  ! This replaces the old incremental approach (header + body) for robustness.
  ! The file describes all time frames upfront based on simulation parameters.
  !
  use specfem_par
  use specfem_par_movie_hdf5

  implicit none
  integer, intent(in) :: num_nodes

  ! local parameters
  integer :: num_elm
  integer :: i_frame, max_surf_frames
  integer :: it_first_surf, it_val
  integer :: ierr
  character(len=MAX_STRING_LEN) :: fname_xdmf_surf
  character(len=MAX_STRING_LEN) :: fname_h5_data_surf_xdmf
  character(len=32) :: it_str

  ! only main process writes out xdmf file
  if (myrank /= 0) return

  ! this routine is used only in single-file mode (no IO sharding)
  if (HDF5_IO_NODES > 1) return

  ! writeout xdmf file for surface movie
  fname_xdmf_surf = trim(OUTPUT_FILES) // "/movie_surface.xmf"
  fname_h5_data_surf_xdmf = "./movie_surface.h5"   ! relative to movie_surface.xmf file

  num_elm = num_nodes / 4

  ! determine number of surface movie frames from time-stepping parameters
  ! frames are written whenever mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0
  ! and it is within [MOVIE_START,MOVIE_STOP]
  max_surf_frames = 0
  it_first_surf = 0

  if (NTSTEP_BETWEEN_FRAMES > 0) then
    ! first time step (within this run) that satisfies the movie sampling rule
    it_first_surf = MOVIE_START + ((max(it_begin, MOVIE_START) - MOVIE_START + NTSTEP_BETWEEN_FRAMES - 1) &
                                   / NTSTEP_BETWEEN_FRAMES) * NTSTEP_BETWEEN_FRAMES

    ! limit movie duration to the actual simulated range
    if (min(it_end, MOVIE_STOP) >= it_first_surf) then
      max_surf_frames = (min(it_end, MOVIE_STOP) - it_first_surf) / NTSTEP_BETWEEN_FRAMES + 1
    else
      max_surf_frames = 0
    endif
  endif

  ! open file with error handling
  open(unit=xdmf_surf, file=trim(fname_xdmf_surf), status='replace', &
       action='write', iostat=ierr, recl=512)
  if (ierr /= 0) then
    print *, 'Error: cannot open XDMF file for surface movie: ', trim(fname_xdmf_surf)
    call exit_mpi(myrank, 'Error opening surface movie XDMF file')
  endif

  ! Write XDMF header
  write(xdmf_surf,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_surf,'(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_surf,'(a)') '<Xdmf Version="3.0">'
  write(xdmf_surf,'(a)') '<Domain Name="mesh">'

  ! Define shared topology and geometry (referenced by all time grids)
  write(xdmf_surf,'(a)') '<Topology Name="topo" TopologyType="Quadrilateral" NumberOfElements="'//trim(i2c(num_elm))//'"/>'
  write(xdmf_surf,'(a)') '<Geometry Name="geom" GeometryType="X_Y_Z">'
  write(xdmf_surf,'(a)') '  <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,'(a)') '    '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/x'
  write(xdmf_surf,'(a)') '  </DataItem>'
  write(xdmf_surf,'(a)') '  <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,'(a)') '    '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/y'
  write(xdmf_surf,'(a)') '  </DataItem>'
  write(xdmf_surf,'(a)') '  <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,'(a)') '    '//trim(fname_h5_data_surf_xdmf)//':/surf_coord/z'
  write(xdmf_surf,'(a)') '  </DataItem>'
  write(xdmf_surf,'(a)') '</Geometry>'

  ! Temporal collection grid
  write(xdmf_surf,'(a)') '<Grid Name="SurfaceMovie" GridType="Collection" CollectionType="Temporal">'

  ! Write all time frames
  do i_frame = 0, max_surf_frames - 1
    it_val = it_first_surf + i_frame * NTSTEP_BETWEEN_FRAMES
    it_str = i2c(it_val)

    write(xdmf_surf,'(a)') '  <Grid Name="it_'//trim(it_str)//'" GridType="Uniform">'
    write(xdmf_surf,'(a)') '    <Time Value="'//trim(r2c(sngl((it_val-1)*DT-t0)))//'" />'
    write(xdmf_surf,'(a)') '    <Topology Reference="/Xdmf/Domain/Topology" />'
    write(xdmf_surf,'(a)') '    <Geometry Reference="/Xdmf/Domain/Geometry" />'

    ! ux attribute
    write(xdmf_surf,'(a)') '    <Attribute Name="ux" AttributeType="Scalar" Center="Node">'
    write(xdmf_surf,'(a)') '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
    write(xdmf_surf,'(a)') '        '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/ux'
    write(xdmf_surf,'(a)') '      </DataItem>'
    write(xdmf_surf,'(a)') '    </Attribute>'

    ! uy attribute
    write(xdmf_surf,'(a)') '    <Attribute Name="uy" AttributeType="Scalar" Center="Node">'
    write(xdmf_surf,'(a)') '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
    write(xdmf_surf,'(a)') '        '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/uy'
    write(xdmf_surf,'(a)') '      </DataItem>'
    write(xdmf_surf,'(a)') '    </Attribute>'

    ! uz attribute
    write(xdmf_surf,'(a)') '    <Attribute Name="uz" AttributeType="Scalar" Center="Node">'
    write(xdmf_surf,'(a)') '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
    write(xdmf_surf,'(a)') '        '//trim(fname_h5_data_surf_xdmf)//':/it_'//trim(it_str)//'/uz'
    write(xdmf_surf,'(a)') '      </DataItem>'
    write(xdmf_surf,'(a)') '    </Attribute>'

    write(xdmf_surf,'(a)') '  </Grid>'
  enddo

  ! Close temporal collection and document
  write(xdmf_surf,'(a)') '</Grid>'
  write(xdmf_surf,'(a)') '</Domain>'
  write(xdmf_surf,'(a)') '</Xdmf>'

  close(xdmf_surf, iostat=ierr)
  if (ierr /= 0) then
    print *, 'Warning: error closing XDMF file for surface movie'
  endif

  end subroutine write_xdmf_surface_complete

#endif


#ifdef USE_HDF5

  subroutine write_xdmf_surface_shards(num_nodes)
  !
  ! Writes complete XDMF files for surface movie in multi-IO mode.
  ! Creates one XDMF file per IO shard, each referencing:
  !   - Shared geometry from movie_surface.h5
  !   - Data from movie_surface.io<k>.h5 for shard k
  !

  use specfem_par
  use specfem_par_movie_hdf5

  implicit none

  integer, intent(in) :: num_nodes

  ! local parameters
  integer :: io_id
  integer :: i_frame, max_surf_frames
  integer :: it_first_surf, it_val
  integer :: ierr
  character(len=MAX_STRING_LEN) :: fname_xdmf_surf
  character(len=MAX_STRING_LEN) :: fname_h5_geom
  character(len=MAX_STRING_LEN) :: fname_h5_data
  character(len=32) :: it_str

  ! only main rank writes XDMF
  if (myrank /= 0) return

  ! only meaningful in multi-IO mode
  if (HDF5_IO_NODES <= 1) return

  ! determine number of surface movie frames from time-stepping parameters
  ! frames are written whenever mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0
  ! and it is within [MOVIE_START,MOVIE_STOP]
  max_surf_frames = 0
  it_first_surf = 0

  if (NTSTEP_BETWEEN_FRAMES > 0) then
    ! first time step (within this run) that satisfies the movie sampling rule
    it_first_surf = MOVIE_START + ((max(it_begin, MOVIE_START) - MOVIE_START + NTSTEP_BETWEEN_FRAMES - 1) &
                                   / NTSTEP_BETWEEN_FRAMES) * NTSTEP_BETWEEN_FRAMES

    ! limit movie duration to the actual simulated range
    if (min(it_end, MOVIE_STOP) >= it_first_surf) then
      max_surf_frames = (min(it_end, MOVIE_STOP) - it_first_surf) / NTSTEP_BETWEEN_FRAMES + 1
    else
      max_surf_frames = 0
    endif
  endif

  fname_h5_geom = "./movie_surface.h5"

  do io_id = 0, HDF5_IO_NODES-1

    ! XDMF file name for this shard
    fname_xdmf_surf = trim(OUTPUT_FILES)//"/movie_surface_io"//trim(i2c(io_id))//".xmf"
    fname_h5_data   = "./movie_surface.io"//trim(i2c(io_id))//".h5"

    ! open file with error handling
    open(unit=xdmf_surf, file=trim(fname_xdmf_surf), status='replace', &
         action='write', iostat=ierr, recl=512)
    if (ierr /= 0) then
      print *, 'Error: cannot open XDMF file for surface movie shard: ', trim(fname_xdmf_surf)
      call exit_mpi(myrank, 'Error opening surface movie XDMF shard file')
    endif

    write(xdmf_surf,'(a)') '<?xml version="1.0" ?>'
    write(xdmf_surf,'(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(xdmf_surf,'(a)') '<Xdmf Version="3.0">'
    write(xdmf_surf,'(a)') '<Domain>'
    write(xdmf_surf,'(a)') '<Grid Name="SurfaceMovie" GridType="Collection" CollectionType="Temporal">'

    do i_frame = 0, max_surf_frames-1

      it_val = it_first_surf + i_frame * NTSTEP_BETWEEN_FRAMES
      it_str = i2c(it_val)

      write(xdmf_surf,'(a)') '  <Grid Name="it_'//trim(it_str)//'" GridType="Uniform">'
      write(xdmf_surf,'(a)') '    <Time Value="'//trim(r2c(sngl((it_val-1)*DT-t0)))//'" />'

      write(xdmf_surf,'(a)') '    <Topology TopologyType="Polyvertex" NumberOfElements="'//trim(i2c(num_nodes))//'"/>'

      ! Note: GeometryType X_Y_Z expects separate x, y, z arrays (not interleaved XYZ)
      write(xdmf_surf,'(a)') '    <Geometry GeometryType="X_Y_Z">'
      write(xdmf_surf,'(a)') '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                             //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
      write(xdmf_surf,'(a)') '        '//trim(fname_h5_geom)//':/surf_coord/x'
      write(xdmf_surf,'(a)') '      </DataItem>'
      write(xdmf_surf,'(a)') '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                             //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
      write(xdmf_surf,'(a)') '        '//trim(fname_h5_geom)//':/surf_coord/y'
      write(xdmf_surf,'(a)') '      </DataItem>'
      write(xdmf_surf,'(a)') '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                             //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
      write(xdmf_surf,'(a)') '        '//trim(fname_h5_geom)//':/surf_coord/z'
      write(xdmf_surf,'(a)') '      </DataItem>'
      write(xdmf_surf,'(a)') '    </Geometry>'

      write(xdmf_surf,'(a)') '    <Attribute Name="velocity" AttributeType="Vector" Center="Node">'
      write(xdmf_surf,'(a)') '      <DataItem ItemType="Function" Function="JOIN($0,$1,$2)" Dimensions="' &
                             //trim(i2c(num_nodes))//' 3">'
      write(xdmf_surf,'(a)') '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                             //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
      write(xdmf_surf,'(a)') '          '//trim(fname_h5_data)//':/it_'//trim(it_str)//'/ux'
      write(xdmf_surf,'(a)') '        </DataItem>'
      write(xdmf_surf,'(a)') '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                             //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
      write(xdmf_surf,'(a)') '          '//trim(fname_h5_data)//':/it_'//trim(it_str)//'/uy'
      write(xdmf_surf,'(a)') '        </DataItem>'
      write(xdmf_surf,'(a)') '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                             //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
      write(xdmf_surf,'(a)') '          '//trim(fname_h5_data)//':/it_'//trim(it_str)//'/uz'
      write(xdmf_surf,'(a)') '        </DataItem>'
      write(xdmf_surf,'(a)') '      </DataItem>'
      write(xdmf_surf,'(a)') '    </Attribute>'

      write(xdmf_surf,'(a)') '  </Grid>'

    enddo

    write(xdmf_surf,'(a)') '</Grid>'
    write(xdmf_surf,'(a)') '</Domain>'
    write(xdmf_surf,'(a)') '</Xdmf>'

    close(xdmf_surf, iostat=ierr)
    if (ierr /= 0) then
      print *, 'Warning: error closing XDMF file for surface movie shard: ', io_id
    endif

  enddo

  end subroutine write_xdmf_surface_shards

#endif

! Note: The old write_xdmf_surface_body subroutine has been removed.
! It used a fragile incremental file-writing approach with magic line numbers.
! The new write_xdmf_surface_complete subroutine writes the complete XDMF file
! at initialization, which is more robust and maintainable.
