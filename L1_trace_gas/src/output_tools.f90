module output_tools
  use netcdf
  use tell_module
  use tio_module
  use OMSAO_precision_module

  implicit none
  private

  public create_output_file, close_output_file, &
    write_radfit_output, write_fitting_statistics

  type (tiof_object_type), private, target :: primary_output_file

contains

  subroutine write_coordinate_vars (obj, num_steps, num_xtrack, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    integer, intent(in) :: num_steps, num_xtrack
    integer, intent(inout) :: errstat

    integer, dimension(num_xtrack) :: xtrack_indices
    integer, dimension(num_steps) :: step_indices
    integer :: i

    if (errstat < 0) return

    ! FIXME: eventually, this will be something like
    ! step_indices=[mirror_step_beg, ..., mirror_step_end]
    ! where mirror_step_beg/end are granule-specific

    step_indices = [(i, i=0,num_steps-1)]
    call tiof_put1d_i4 (obj, tempo_dim_step, 0, num_steps, step_indices, errstat)

    xtrack_indices = [(i, i=0,num_xtrack-1)]
    call tiof_put1d_i4 (obj, tempo_dim_xtrack, 0, num_xtrack, xtrack_indices, errstat)

  end subroutine write_coordinate_vars

  subroutine create_output_file (filename, num_steps, num_xtrack, num_swlevels, &
                                 errstat)
    !USE OMSAO_parameters_module, ONLY: NWAVEL_MAX, nUTCdim
    !USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
    !USE OMSAO_omidata_module,   ONLY: nclenfit, n_comm_wvl
    !USE OMSAO_variables_module, ONLY: n_fitvar_rad, n_rad_wvl
    use ctrlvars, only: yn_diagnostic_run !, yn_refseccor, yn_scat_weights
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

    ! data field variables with optional attribute lists:
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(tempo_var_longitude) &
                              //' '//trim(tempo_var_latitude))
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_column_amount, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "column amount", &
                              units = "molec/cm2", &
                              valid_range = [-1.e30_r8, 1.e30_r8], &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_column_uncert, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "column amount uncertainty", &
                              units = "molec/cm2", &
                              valid_range = [-1.e30_r8, 1.e30_r8], &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_fit_rms, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "fit rms", &
                              valid_range = [0.0_r8, 1.e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_fit_convergence_flag, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "fit convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8])

    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_time, &
                              nf90_double, &
                              dimids = [dimids_xtrack_step(2)],  &
                              comment = "exposure start time", &
                              units = "s", &
                              valid_range = [0.0_r8, 1.e30_r8])

    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = tempo_var_latitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_latitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "latitude at pixel center", &
                              units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              attlist=att_latbnd)

    call tiof_attlist_append (att_lonbnd, errstat, "bounds", &
                              att_text = tempo_var_longitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_longitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "longitude at pixel center", &
                              units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              attlist=att_lonbnd)

    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_sz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_sa_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_vz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_va_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8])

    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_main_dqf, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "main data quality flag", &
                              valid_range = [-1.0_r8, 2.0_r8])

    if (yn_diagnostic_run) then
      call tiof_varlist_append (varlist, errstat, &
                                tempo_var_fit_iteration_count, &
                                nf90_short, &
                                dimids = dimids_xtrack_step,  &
                                comment = "radiance fit iteration count", &
                                valid_range = [0.0_r8, 32767.0_r8])      
    endif

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
    use OMSAO_omidata_module, only : input_vars_type, result_vars_type
    use ctrlvars, only: yn_diagnostic_run !, yn_refseccor, yn_scat_weights
    implicit none

    integer, intent(in) :: iline, nblock, nxtrack
    type (input_vars_type), intent(in) :: input_vars
    type (result_vars_type), intent(in) :: result_vars
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    if (errstat < 0) return

    ! result_vars
    call tiof_put2d_r8 (obj, tempo_var_column_amount, iline, nblock, &
                        result_vars % column_amount (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r8 (obj, tempo_var_column_uncert, iline, nblock, &
                        result_vars % column_uncert (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r8 (obj, tempo_var_fit_rms, iline, nblock, &
                        result_vars % fit_rms (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_i2 (obj, tempo_var_fit_convergence_flag, iline, nblock, &
                        result_vars % fit_convergence_flag (1:nxtrack, 0:nblock-1), errstat)
    
    if (yn_diagnostic_run) then
      call tiof_put2d_i2 (obj, tempo_var_fit_iteration_count, iline, nblock, &
                          result_vars % fit_iteration_count (1:nxtrack, 0:nblock-1), errstat)
    endif

    ! input_vars
    call tiof_put1d_r8 (obj, tempo_var_time, iline, nblock, &
                        input_vars % time (0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tempo_var_longitude, iline, nblock, &
                        input_vars % longitude (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tempo_var_latitude, iline, nblock, &
                        input_vars % latitude (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tempo_var_sz_angle, iline, nblock, &
                        input_vars % solar_zenith (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tempo_var_sa_angle, iline, nblock, &
                        input_vars % solar_azimuth (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tempo_var_vz_angle, iline, nblock, &
                        input_vars % viewing_zenith (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tempo_var_va_angle, iline, nblock, &
                        input_vars % viewing_azimuth (1:nxtrack, 0:nblock-1), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, "write_radfit_output: failed", errstat)
      return
    endif

  end subroutine write_radfit_output

  subroutine write_fitting_statistics (stats, errstat)
    use omi_pge_fitting_aux, only : fitting_statistics_type
    implicit none

    type (fitting_statistics_type), intent(in) :: stats
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file
    type (tiof_attlist_type) :: attlist

    call tiof_attlist_append (attlist, errstat, "num_crosstrack_pixels", &
                             att_i4=[stats % num_crosstrack_pixels])
    call tiof_attlist_append (attlist, errstat, "num_scan_lines", &
                             att_i4=[stats % num_scan_lines])
    call tiof_attlist_append (attlist, errstat, "num_good_input", &
                             att_i4=[stats % num_good_input])
    call tiof_attlist_append (attlist, errstat, "num_good_output", &
                             att_i4=[stats % num_good_output])
    call tiof_attlist_append (attlist, errstat, "num_suspect_output", &
                             att_i4=[stats % num_suspect_output])
    call tiof_attlist_append (attlist, errstat, "num_bad_output", &
                             att_i4=[stats % num_bad_output])
    call tiof_attlist_append (attlist, errstat, "num_converged", &
                             att_i4=[stats % num_converged])
    call tiof_attlist_append (attlist, errstat, "num_failed_convergence", &
                             att_i4=[stats % num_failed_convergence])
    call tiof_attlist_append (attlist, errstat, "num_exceeded_iterations", &
                             att_i4=[stats % num_exceeded_iterations])
    call tiof_attlist_append (attlist, errstat, "num_out_of_bounds", &
                             att_i4=[stats % num_out_of_bounds])
    call tiof_attlist_append (attlist, errstat, "percent_good_output", &
                             att_r4=[stats % percent_good_output])
    call tiof_attlist_append (attlist, errstat, "percent_bad_output", &
                             att_r4=[stats % percent_bad_output])
    call tiof_attlist_append (attlist, errstat, "percent_suspect_output", &
                             att_r4=[stats % percent_suspect_output])

    call tiof_def_atts (obj, attlist, nf90_global, errstat)

    call tiof_put2d_i2 (obj, tempo_var_main_dqf, 0, stats % num_scan_lines, &
                        stats % quality_flag (1:stats % num_crosstrack_pixels, &
                                              0:stats % num_scan_lines-1), &
                        errstat)

    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "write_fitting_statistics: writing fitting statistics", &
                       errstat)
      return
    endif

  end subroutine write_fitting_statistics

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
