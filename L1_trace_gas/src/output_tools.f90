module output_tools
  use netcdf
  use tell_module
  use tio_module
  use OMSAO_precision_module
  use ctrlvars, only: yn_diagnostic_run, yn_refseccor, yn_scat_weights
  use OMSAO_omidata_module, only : input_vars_type, result_vars_type

  implicit none
  private

  public create_output_file, close_output_file, write_radfit_output

  type (tiof_object_type), private, target :: primary_output_file

contains

  subroutine write_coordinate_vars (obj, num_steps, num_xtrack, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    integer, intent(in) :: num_steps, num_xtrack
    integer, intent(inout) :: errstat

    integer, dimension(num_xtrack) :: dim_xtrack
    integer, dimension(num_steps) :: dim_step
    integer :: i

    if (errstat < 0) return

    ! FIXME: eventually, want: dim_step=[mirror_step_start, ..., mirror_step_end]

    dim_step = [(i, i=0,num_steps-1)]
    call tiof_put1d_i4 (obj, tempo_dim_step, 0, num_steps, dim_step, errstat)

    dim_xtrack = [(i, i=0,num_xtrack-1)]
    call tiof_put1d_i4 (obj, tempo_dim_xtrack, 0, num_xtrack, dim_xtrack, errstat)

  end subroutine write_coordinate_vars

  subroutine create_output_file (filename, num_steps, num_xtrack, num_swlevels, &
                                 errstat)
    !USE OMSAO_parameters_module, ONLY: NWAVEL_MAX, nUTCdim
    !USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
    !USE OMSAO_omidata_module,   ONLY: nclenfit, n_comm_wvl
    !USE OMSAO_variables_module, ONLY: n_fitvar_rad

    implicit none
    character (len=*), intent(in) :: filename
    integer (kind=i4), intent(in) :: num_steps, num_xtrack, num_swlevels
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file
    type (tiof_dimlist_type) :: dimlist
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord, att_latbnd, att_lonbnd
    integer, dimension(2) :: dimids_xtrack_step

    if (errstat < 0) return

    ! create a file
    call tiof_create (obj, filename, ior(nf90_clobber,nf90_netcdf4), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: creating file "//trim(filename), &
                       errstat)
      return
    endif

    ! Define a dimension list.
    call tiof_dimlist_append (dimlist, tempo_dim_step, num_steps, errstat)
    call tiof_dimlist_append (dimlist, tempo_dim_xtrack, num_xtrack, errstat)
    call tiof_dimlist_append (dimlist, tempo_dim_swt_levels, num_swlevels, errstat)
    call tiof_def_dims (obj, dimlist, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: defining dimensions in "//trim(filename), &
                       errstat)
      return
    endif

    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              [tempo_dim_xtrack, tempo_dim_step], &
                              dimids_xtrack_step, &
                              errstat)

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    ! netcdf coordinate variables:
    call tiof_varlist_append (varlist, errstat, tempo_dim_xtrack, nf90_int, &
                             dimids=[dimids_xtrack_step(1)])
    call tiof_varlist_append (varlist, errstat, tempo_dim_step, nf90_int, &
                             dimids=[dimids_xtrack_step(2)])

    ! data field variables with their attribute lists:
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(tempo_var_longitude) &
                              //' '//trim(tempo_var_latitude))
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_column_amount, &
                              nf90_float, &   ! data type in output file
                              dimids = dimids_xtrack_step,  &
                              comment = "column amount", &
                              units = "molec/cm2", &
                              valid_range = [-1.e30_r8, 1.e30_r8], &
                              attlist=att_coord &
                             )

    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = tempo_var_latitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_latitude, &
                              nf90_float, &   ! data type in output file
                              dimids = dimids_xtrack_step,  &
                              comment = "latitude at pixel center", &
                              units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              attlist=att_latbnd &
                             )

    call tiof_attlist_append (att_lonbnd, errstat, "bounds", &
                              att_text = tempo_var_latitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_longitude, &
                              nf90_float, &   ! data type in output file
                              dimids = dimids_xtrack_step,  &
                              comment = "longitude at pixel center", &
                              units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              attlist=att_lonbnd &
                             )

    call tiof_def_vars (obj, varlist, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: defining variables in "//trim(filename), &
                       errstat)
      return
    endif

    call write_coordinate_vars (obj, num_steps, num_xtrack, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: writing coordinate variables to "//trim(filename), &
                       errstat)
      return
    endif

  end subroutine create_output_file

  subroutine write_radfit_output (iline, nblock, nxtrack, &
                                  input_vars, result_vars, errstat)
    implicit none

    integer, intent(in) :: iline, nblock, nxtrack
    type (input_vars_type), intent(in) :: input_vars
    type (result_vars_type), intent(in) :: result_vars
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    if (errstat < 0) return

    call tiof_put2d_r8 (obj, tempo_var_column_amount, iline, nblock, &
                        result_vars % column_amount (1:nxtrack, 0:nblock-1), errstat)

    call tiof_put2d_r4 (obj, tempo_var_longitude, iline, nblock, &
                        input_vars % longitude (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tempo_var_latitude, iline, nblock, &
                        input_vars % latitude (1:nxtrack, 0:nblock-1), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, "write_radfit_output: failed", errstat)
      return
    endif

  end subroutine write_radfit_output

  subroutine close_output_file (errstat)
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    call tiof_close (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_error, "close_output_file failed", errstat)
    endif

  end subroutine close_output_file

end module output_tools
