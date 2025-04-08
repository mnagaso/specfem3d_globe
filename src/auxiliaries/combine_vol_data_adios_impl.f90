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

module combine_vol_data_adios_mod

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! member variables
  public
  logical is_forward ! forward array
  logical is_undo !

contains

!=============================================================================
!> Print help message.
subroutine print_usage_adios()

  implicit none
  print *, ' Usage: '
  print *, '   xcombine_data slice_list varname var_file mesh_file output_dir high/low-resolution region'
  print *
  print *, ' with'
  print *, '   slice_list   - text file containing slice numbers to combine (or use name "all" for all slices)'
  print *, '   varname      - possible varnames are: '
  print *, '                    rho, vp, vs, kappastore, mustore, alpha_kl, beta_kl, etc.'
  print *, '   var_file     - datafile that holds array, as real(kind=CUSTOM_REAL):: varname(NGLLX,NGLLY,NGLLZ,NSPEC),'
  print *, '                  (e.g. OUTPUT_FILES/kernels.bp)'
  print *, '   mesh_file    - are used to link variable to the topology (e.g. DATABASES_MPI/solver_data.bp)'
  print *, '   output_dir   - indicates where var_name.vtk will be written'
  print *, '   high/low res - give 0 for low resolution and 1 for high resolution'
  print *, '   region       - (optional) region number, only use 1 == crust/mantle, 2 == outer core, 3 == inner core'
  print *, '   iter         - (optional) iteration number'
  print *

  stop ' Reenter command line options'

end subroutine print_usage_adios

!=============================================================================
!> Interpret command line arguments

subroutine read_args_adios(arg, var_name, value_file_name, mesh_file_name, slice_list_name, &
                           outdir, ires, iregion, iiter)

  use constants, only: IIN,MAX_STRING_LEN

  implicit none
  ! Arguments
  character(len=*), intent(in) :: arg(:)
  integer, intent(out) :: ires, iregion, iiter
  character(len=*), intent(out) :: var_name, value_file_name, mesh_file_name, &
                                   outdir, slice_list_name

  ! initializes
  iregion = 0
  iiter = 0

  ! gets arguments
  if ((command_argument_count() == 6) .or. (command_argument_count() == 7) &
                                      .or. (command_argument_count() == 8)) then
    slice_list_name = arg(1)
    var_name = arg(2)
    value_file_name = arg(3)
    mesh_file_name = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
  else
    call print_usage_adios()
  endif

  if ((command_argument_count() >= 7)) then
    read(arg(7),*) iregion
  endif

  if ((command_argument_count() == 8)) then
    read(arg(8),*) iiter
  endif

  !debug
  !print *,'debug: read adios: arguments: ',trim(slice_list_name),"|",trim(var_name),"|", &
  !        trim(value_file_name),"|",trim(mesh_file_name),"|",trim(outdir),"|",ires,"|",iregion

end subroutine read_args_adios


!=============================================================================
!> Open ADIOS value and mesh files, read mode

subroutine init_adios(value_file_name, mesh_file_name, i_iter)
  use constants, only: ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE

  implicit none
  ! Parameters
  character(len=*), intent(in) :: value_file_name, mesh_file_name
  integer, intent(in) :: i_iter

  ! debug
  logical, parameter :: DEBUG = .false.
  integer :: nglob,nspec

  if (i_iter >= 0) then
    is_forward = .true.
    ! if file name contains "undoatt" then it is an undo array
    if (index(value_file_name, "undoatt") > 0) then
      is_undo = .true.
    else
      is_undo = .false.
    endif
    ! print if it is a forward array and undo array
    print * , 'is_forward = ',is_forward
    print * , 'is_undo = ',is_undo
  endif

  ! initializes adios
  call initialize_adios()

  ! initializes read method and opens mesh file (using default handle)
  call init_adios_group(myadios_group,"MeshReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,mesh_file_name)

  ! opens second adios file for reading data values
  if (.not. is_undo) then
    call init_adios_group(myadios_val_group,"ValReader")
  else
    if (.not. ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) then
      ! not supported
      stop 'Error: adios save all snapshots in one file not supported'
    else
      ! undoatt file
      call init_adios_group_undo_att(myadios_val_group,"SPECFEM3D_GLOBE_FORWARD_ARRAYS_UNDOATT")
    endif
  endif
  call open_file_adios_read(myadios_val_file,myadios_val_group,value_file_name)

  ! debug output list variables and attributs in mesh file
  if (DEBUG) then
    call show_adios_file_variables(myadios_file,myadios_group,mesh_file_name)
    call show_adios_file_variables(myadios_val_file,myadios_val_group,value_file_name)
    ! checks reading nglob/nspec values
    print *,'mesh file: ',trim(mesh_file_name)
    call read_adios_scalar(myadios_file,myadios_group,0,"reg1/nglob",nglob)
    call read_adios_scalar(myadios_file,myadios_group,0,"reg1/nspec",nspec)
    print *,'  nglob/nspec = ',nglob,"/",nspec
  endif

end subroutine init_adios


!=============================================================================
!> Open ADIOS value and mesh files, read mode

subroutine clean_adios()

  implicit none

  ! closes file with data values
  call close_file_adios_read(myadios_val_file)

  ! closes default file and finalizes read method
  call close_file_adios_read_and_finalize_method(myadios_file)

  ! finalizes
  call finalize_adios()

end subroutine clean_adios


!=============================================================================

subroutine read_scalars_adios_mesh(iproc, ir, nglob, nspec)

  implicit none
  ! Parameters
  integer, intent(in) :: iproc, ir
  integer, intent(out) :: nglob, nspec
  ! Variables
  !integer(kind=8) :: sel
  character(len=80) :: reg_name
  ! debug
  logical, parameter :: DEBUG = .false.
  !character(len=1024) :: err_message
  !integer :: ier

  ! region name
  write(reg_name,"('reg',i1,'/')") ir

  ! debug output
  if (DEBUG) then
    call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nglob",nglob)
    ! note: adios_get_scalar here retrieves the same nglob for everyone (from writer rank 0)
    !call adios_get_scalar(myadios_file, trim(reg_name)//"nglob", nglob, ier)
    !if (ier /= 0) then
    !  call adios_errmsg(err_message)
    !  print *,'Error adios: could not read parameter: ',trim(reg_name)//"nglob"
    !  print *,trim(err_message)
    !  stop 'Error adios helper read scalar'
    !endif

    call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nspec",nspec)
    !call adios_get_scalar(myadios_file, trim(reg_name)//"nspec", nspec, ier)
    !if (ier /= 0) then
    !  call adios_errmsg(err_message)
    !  print *,'Error adios: could not read parameter: ',trim(reg_name)//"nspec"
    !  print *,trim(err_message)
    !  stop 'Error adios helper read scalar'
    !endif

    print *,'read nglob/nspec: ',iproc,ir,nglob,nspec
  endif

  ! reads nglob & nspec
  ! the following adios calls allow to have different nglob/nspec values for different processes (iproc)
  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nglob",nglob)
  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nspec",nspec)

end subroutine read_scalars_adios_mesh


!=============================================================================

subroutine read_coordinates_adios_mesh(iproc, ir, nglob, nspec, &
                                       xstore, ystore, zstore, ibool)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none
  ! Parameters
  integer, intent(in) :: iproc, ir, nglob, nspec
  real(kind=CUSTOM_REAL),dimension(:), intent(inout) :: xstore, ystore, zstore
  integer, dimension(:,:,:,:), intent(inout) :: ibool
  ! Variables
  character(len=80) :: reg_name
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel_coord, sel_ibool
  integer(kind=8) :: offset_coord, offset_ibool

  write(reg_name,"('reg',i1,'/')") ir

  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "ibool/offset",offset_ibool)
  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "x_global/offset",offset_coord)

  start(1) = offset_ibool
  count(1) = NGLLX * NGLLY * NGLLZ * nspec
  call set_selection_boundingbox(sel_ibool , start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel_ibool, start, count, trim(reg_name) // "ibool/array", ibool)

  start(1) = offset_coord
  count(1) = nglob
  call set_selection_boundingbox(sel_coord, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel_coord, start, count, trim(reg_name) // "x_global/array", xstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel_coord, start, count, trim(reg_name) // "y_global/array", ystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel_coord, start, count, trim(reg_name) // "z_global/array", zstore)

  call read_adios_perform(myadios_file)

  call delete_adios_selection(sel_ibool)
  call delete_adios_selection(sel_coord)

end subroutine read_coordinates_adios_mesh


!=============================================================================
!> reads in data from ADIOS value file

subroutine read_values_adios(var_name, iproc, ir, i_iter, ires, nglob, nspec, ibool, data)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,IREGION_CRUST_MANTLE,IREGION_INNER_CORE,IREGION_OUTER_CORE

  implicit none
  ! Parameters
  character(len=*), intent(in) :: var_name
  integer, intent(in) :: iproc, ir, i_iter, ires, nspec, nglob
  integer, dimension(:,:,:,:), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(inout) :: data
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: data_tmp ! used for reading forward arrays
  ! Variables
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel
  integer(kind=8) :: offset
  character(len=128) :: data_name
  character(len=8) :: reg_name
  logical :: is_kernel
  integer :: iexist
  integer :: di, dj, dk, i, j, k, ispec, iglob

  ! note: we can  visualize
  !         wavespeed arrays (in DATABASES_MPI/model_gll.bp)
  !                          (e.g. rho,vp,vs,vph,vpv,vsh,..)
  !       or
  !         sensitivity kernels (in OUTPUT_FILES/kernels.bp)
  !                             (e.g. rho_kl,alpha_kl,beta_kl,alphah_kl,alphav_kl,betah_kl,..)
  !       or
  !         forward arrays (in DATABASES_MPI/save_forward_arrays.pb) *only with adios
  !         (e.g. rho_forward,alpha_forward,beta_forward,alphah_forward,alphav_forward,betah_forward,..)
  !       or
  !         undo arrays (in DATABASES_MPI/save_forward_arrays_undoatt.bp) *only with adios
  !       unfortunately, they have different naming conventions for different Earth regions:
  !        - wavespeed arrays: reg1/vp/, reg2/vp/, reg3/vp/, ..
  !        - kernels: alpha_kl_crust_mantle/, alpha_kl_outer_core/, alpha_kl_inner_core/,..
  ! here we try to estimate the type by the ending of the variable name given by the user,
  ! i.e. if the ending is "_kl" we assume it is a kernel name
  !
  is_kernel = .false.

  ! check if the variable name is a kernel name
  ! example: alpha_kl checks ending '_kl'
  if (len_trim(var_name) > 3) then
    if (var_name(len_trim(var_name)-2:len_trim(var_name)) == '_kl') then
      is_kernel = .true.
    endif
  endif
  ! example: alpha_kl_crust_mantle checks ending '_kl_crust_mantle'
  if (len_trim(var_name) > 16) then
    if (var_name(len_trim(var_name)-15:len_trim(var_name)) == '_kl_crust_mantle') then
      is_kernel = .true.
    endif
  endif

  ! verify if i_iter is not negative for is_forward
  if (is_forward .and. i_iter < 0) then
    print *,'Error: i_iter is negative for forward array'
    stop 'Error: i_iter is negative for forward array'
  endif

  ! determines full data array name
  if (is_kernel) then
    ! for kernel name: alpha_kl, .. , **_kl only:
    ! adds region name to get kernel_name
    ! for example: var_name = "alpha_kl" -> alpha_kl_crust_mantle
    !
    ! note: this must match the naming convention used
    !       in file save_kernels_adios.F90
    if (var_name(len_trim(var_name)-2:len_trim(var_name)) == '_kl') then
      select case (ir)
      case (IREGION_CRUST_MANTLE)
        data_name = trim(var_name) // "_crust_mantle"
      case (IREGION_OUTER_CORE)
        data_name = trim(var_name) // "_outer_core"
      case (IREGION_INNER_CORE)
        data_name = trim(var_name) // "_inner_core"
      case default
        stop 'Error wrong region code in read_values_adios() routine'
      end select
    else
      ! example: alpha_kl_crust_mantle
      data_name = trim(var_name)
    endif
    print *,'  kernel data name: ',trim(data_name)
  else if (is_forward) then
    ! expeting var_name = "R_**", "accel", "veloc" or "displ"
    ! stop if the name is not correct
    if (index(var_name, "R_") == 0 .and. index(var_name, "accel") == 0 .and. &
        index(var_name, "veloc") == 0 .and. index(var_name, "displ") == 0) then
      print *,'Error: undo array name is not correct'
      stop 'Error: undo array name is not correct'
    endif

    select case (ir)
    case (IREGION_CRUST_MANTLE)
      data_name = trim(var_name) // "_crust_mantle"
    case (IREGION_OUTER_CORE)
      data_name = trim(var_name) // "_outer_core"
    case (IREGION_INNER_CORE)
      data_name = trim(var_name) // "_inner_core"
    case default
      stop 'Error wrong region code in read_values_adios() routine'
    end select

    if (ires == 1) then ! high resolution
      ! set increments
      di = 1
      dj = 1
      dk = 1
    else if (ires == 2) then ! mid. resolution
      di = int((NGLLX-1)/2.0)
      dj = int((NGLLY-1)/2.0)
      dk = int((NGLLZ-1)/2.0)
    else ! high resolution
      di = NGLLX - 1
      dj = NGLLY - 1
      dk = NGLLZ - 1
    endif

  else
    ! for wavespeed name: rho,vp,..
    ! adds region name: var_name = "rho" -> reg1/rho
    write(reg_name,"('reg',i1,'/')") ir
    data_name = trim(reg_name) // trim(var_name)
    print *,'  data name: ',trim(data_name)
  endif

  ! gets data values
  if (is_forward) then
    ! allocate data_tmp
    allocate(data_tmp(NDIM*nglob))

    ! assumes GLL type array size (NGLLX,NGLLY,NGLLZ,nspec)
    call read_adios_array_gll_check_forward(myadios_val_file,myadios_val_group,iproc,nglob,&
                                            trim(data_name),data_tmp,iexist,int(i_iter, kind=8))

    ! recompose the array from glob based to GLL based
    ! data_tmp -> data
    do ispec = 1, nspec
      do k = 1,NGLLZ, dk
        do j = 1,NGLLY, dj
          do i = 1,NGLLX, di
            iglob = ibool(i,j,k,ispec)
            data(i,j,k,ispec) = data_tmp(iglob)
          end do
        end do
      end do
    end do
    deallocate(data_tmp)
  else
    ! reads in data offset
    call read_adios_scalar(myadios_val_file,myadios_val_group,iproc,trim(data_name) // "/offset",offset)

    ! reads in data array
    start(1) = offset
    count(1) = NGLLX * NGLLY * NGLLZ * nspec
    call set_selection_boundingbox(sel , start, count)

    call read_adios_schedule_array(myadios_val_file, myadios_val_group, sel, start, count, trim(data_name) // "/array", data)
    call read_adios_perform(myadios_val_file)

    call delete_adios_selection(sel)
  endif

end subroutine read_values_adios

end module combine_vol_data_adios_mod
