module output_tools
  use netcdf
  use tell_module
  use tio_module
  use OMSAO_precision_module
  use ctrlvars, only: yn_diagnostic_run !, yn_refseccor, yn_scat_weights

  implicit none
  private

  public create_output_file, close_output_file, &
    write_radfit_output, write_fitting_statistics, &
    write_albedo, write_gas_profile, write_scattering_weights, &
    write_amf_correction

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

  subroutine append_amf_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_object_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_xtrack_step_levels
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat < 0) return

    ! lookup dimids for relevant array shapes
    call tiof_dimlist_lookup (dimlist, &
                              [tempo_dim_xtrack, tempo_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tempo_dim_xtrack, tempo_dim_step,tempo_dim_swt_level], &
                              dimids_xtrack_step_levels, &
                              errstat)

    ! append amf variables
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_scattering_weights, &
                              nf90_double, &
                              dimids = dimids_xtrack_step_levels,  &
                              comment = "scattering weights", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_climatology_levels, &
                              nf90_double, &
                              dimids = dimids_xtrack_step_levels,  &
                              comment = "climatology levels", &
                              units = "hPa", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_gas_profile, &
                              nf90_double, &
                              dimids = dimids_xtrack_step_levels,  &
                              comment = "gas profile", &
                              units = "ppb", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_albedo, &
                              nf90_double, &
                              dimids = dimids_xtrack_step,  &
                              comment = "albedo", &
                              valid_range = [-1e30_r8, 1e30_r8])

    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_molecule_specific, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "molecule-specific air mass factor (AMF)", &
                              valid_range = [0.0_r8, 1e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_diagnostic_flag, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "diagnostic flag for molecule-specific air mass factor (AMF)", &
                              valid_range = [-2.0_r8, 13127.0_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_geometric, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "geometric air mass factor (AMF)", &
                              valid_range = [0.0_r8, 1e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_cloud_fraction, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "adjusted cloud fraction for AMF computation", &
                              valid_range = [0.0_r8, 1e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_amf_cloud_pressure, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "adjusted cloud pressure for AMF computation", &
                              valid_range = [0.0_r8, 1e30_r8])

    call tiof_def_vars (obj, varlist, errstat)

  end subroutine append_amf_vars

  subroutine append_diagnostic_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_object_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_var_xtrack_step, dimids_commwvl_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat < 0) return

    ! lookup dimids for relevant array shapes
    call tiof_dimlist_lookup (dimlist, &
                              [tempo_dim_xtrack, tempo_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tempo_dim_fitvar, tempo_dim_xtrack, tempo_dim_step], &
                              dimids_var_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tempo_dim_commwvl, tempo_dim_xtrack, tempo_dim_step], &
                              dimids_commwvl_xtrack_step, &
                              errstat)

    ! append diagnostic variables
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_fit_iteration_count, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "radiance fit iteration count", &
                              valid_range = [0.0_r8, 32767.0_r8])

    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_diag_params, &
                              nf90_float, &
                              dimids = dimids_var_xtrack_step, &
                              comment = "fit parameter", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_diag_errors, &
                              nf90_float, &
                              dimids = dimids_var_xtrack_step, &
                              comment = "fit parameter uncertainty", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_diag_correl, &
                              nf90_float, &
                              dimids = dimids_var_xtrack_step, &
                              comment = "fit parameter correlation", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)

    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_diag_measured_spectrum, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "measured spectrum", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_diag_measured_wavelengths, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "measured wavelength", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_diag_model_spectrum, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "model spectrum", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tempo_var_diag_fit_weights, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "spectrum fit weights", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)

    call tiof_def_vars (obj, varlist, errstat)

  end subroutine append_diagnostic_vars

  subroutine create_output_file (filename, num_steps, num_xtrack, num_swlevels, &
                                 errstat)
    !USE OMSAO_parameters_module, ONLY: NWAVEL_MAX, nUTCdim
    !USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
    USE OMSAO_omidata_module,   ONLY: n_comm_wvl ! nclenfit
    USE OMSAO_variables_module, ONLY: n_fitvar_rad
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
    call tiof_dimlist_append (dimlist, tempo_dim_swt_level, num_swlevels, errstat)
    if (yn_diagnostic_run) then
      call tiof_dimlist_append (dimlist, tempo_dim_fitvar, n_fitvar_rad, errstat)
      call tiof_dimlist_append (dimlist, tempo_dim_commwvl, n_comm_wvl, errstat)
    endif
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
                              tempo_var_fit_rms_residual, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "fit rms residual", &
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

    call append_amf_vars (obj, dimlist, errstat)

    if (yn_diagnostic_run) then
      call append_diagnostic_vars (obj, dimlist, errstat)
      if (errstat < 0) then
        call tell_error (tell_io_write_error, &
                         "create_output_file: defining diagnostic variables in "//trim(filename), &
                         errstat)
        return
      endif
    endif

  end subroutine create_output_file

  subroutine write_radfit_output (iline, nblock, nxtrack, n_fitvar_rad, n_rad_wvl, &
                                  input_vars, result_vars, radfit_diagnostics, &
                                  errstat)
    use OMSAO_omidata_module, only : input_vars_type, result_vars_type, &
      radfit_diagnostics_type
    implicit none

    integer, intent(in) :: iline, nblock, nxtrack, n_fitvar_rad, n_rad_wvl
    type (input_vars_type), intent(in) :: input_vars
    type (result_vars_type), intent(in) :: result_vars
    type (radfit_diagnostics_type), intent(in) :: radfit_diagnostics
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    if (errstat < 0) return

    ! result_vars
    call tiof_put2d_r8 (obj, tempo_var_column_amount, iline, nblock, &
                        result_vars % column_amount (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r8 (obj, tempo_var_column_uncert, iline, nblock, &
                        result_vars % column_uncert (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r8 (obj, tempo_var_fit_rms_residual, iline, nblock, &
                        result_vars % fit_rms_residual (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_i2 (obj, tempo_var_fit_convergence_flag, iline, nblock, &
                        result_vars % fit_convergence_flag (1:nxtrack, 0:nblock-1), errstat)

    if (yn_diagnostic_run) then
      call tiof_put2d_i2 (obj, tempo_var_fit_iteration_count, iline, nblock, &
                          result_vars % fit_iteration_count (1:nxtrack, 0:nblock-1), errstat)

      call tiof_put3d_r8 (obj, tempo_var_diag_params, iline, nblock, &
                          radfit_diagnostics % params(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj, tempo_var_diag_errors, iline, nblock, &
                          radfit_diagnostics % errors(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj, tempo_var_diag_correl, iline, nblock, &
                          radfit_diagnostics % correl(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)

      call tiof_put3d_r8 (obj, tempo_var_diag_model_spectrum, iline, nblock, &
                          radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 1, 0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj, tempo_var_diag_measured_spectrum, iline, nblock, &
                          radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 2, 0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj, tempo_var_diag_measured_wavelengths, iline, nblock, &
                          radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 3, 0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj, tempo_var_diag_fit_weights, iline, nblock, &
                          radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 4, 0:nblock-1), &
                          errstat)
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

  subroutine write_albedo (albedo, nxtrack, ntimes, errstat)
    implicit none

    real (kind=r8), dimension (1:nxtrack, 0:ntimes-1), intent(in) :: albedo
    integer, intent(in) :: nxtrack, ntimes
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    if (errstat < 0) return
    call tiof_put2d_r8 (obj, tempo_var_amf_albedo, 0, ntimes, &
                        albedo (1:nxtrack, 0:ntimes-1), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, "in write_albedo", errstat)
      return
    endif

  end subroutine write_albedo

  subroutine write_gas_profile (gas_profile, climatology_levels, &
                                nxtrack, ntimes, nlevels, errstat)
    implicit none

    real (kind=r8), dimension (1:nxtrack, 0:ntimes-1, 1:nlevels), intent(in) :: gas_profile
    real (kind=r8), dimension (1:nxtrack, 0:ntimes-1, 1:nlevels), intent(in) :: climatology_levels
    integer, intent(in) :: nxtrack, ntimes, nlevels
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    if (errstat < 0) return
    call tiof_put3d_r8 (obj, tempo_var_amf_gas_profile, 0, nlevels, &
                        gas_profile(1:nxtrack, 0:ntimes-1, 1:nlevels), errstat)
    call tiof_put3d_r8 (obj, tempo_var_amf_climatology_levels, 0, nlevels, &
                        climatology_levels(1:nxtrack, 0:ntimes-1, 1:nlevels), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, "in write_gas_profile", errstat)
      return
    endif

  end subroutine write_gas_profile

  subroutine write_scattering_weights (scattw, nxtrack, ntimes, nlevels, errstat)
    implicit none

    real (kind=r8), dimension (1:nxtrack, 0:ntimes-1, 1:nlevels), intent(in) :: scattw
    integer, intent(in) :: nxtrack, ntimes, nlevels
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    if (errstat < 0) return
    call tiof_put3d_r8 (obj, tempo_var_amf_scattering_weights, 0, nlevels, &
                        scattw (1:nxtrack, 0:ntimes-1, 1:nlevels), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, "in write_scattering_weights", errstat)
      return
    endif

  end subroutine write_scattering_weights

  subroutine write_amf_correction (pge_idx, nxtrack, ntimes, amf_corr, &
                                   amf_corr_column, amf_corr_column_uncertainty, &
                                   errstat)
    use OMSAO_omidata_module, only : amf_correction_type
    use OMSAO_indices_module, only: pge_hcho_idx, pge_gly_idx
    implicit none

    integer, intent(in) :: pge_idx, nxtrack, ntimes
    type (amf_correction_type), intent(in) :: amf_corr
    real (kind=r8), dimension(1:nxtrack,0:ntimes-1), intent(in) :: amf_corr_column
    real (kind=r8), dimension(1:nxtrack,0:ntimes-1), intent(in) :: amf_corr_column_uncertainty
    integer, intent(inout) :: errstat

    type (tiof_object_type), pointer :: obj => primary_output_file

    if (errstat < 0) return

    call tiof_put2d_i2 (obj, tempo_var_amf_diagnostic_flag, 0, ntimes, &
                        amf_corr % diagnostic_flag (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, tempo_var_amf_geometric, 0, ntimes, &
                        amf_corr % amf_geometric (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, tempo_var_amf_molecule_specific, 0, ntimes, &
                        amf_corr % amf_molecule_specific (1:nxtrack, 0:ntimes-1), errstat)

    if (pge_idx == pge_hcho_idx .or. pge_idx == pge_gly_idx) then
      call tiof_put2d_r8 (obj, tempo_var_amf_cloud_fraction, 0, ntimes, &
                          amf_corr % cloud_fraction (1:nxtrack, 0:ntimes-1), errstat)
      call tiof_put2d_r8 (obj, tempo_var_amf_cloud_pressure, 0, ntimes, &
                          amf_corr % cloud_pressure (1:nxtrack, 0:ntimes-1), errstat)
    endif

    call tiof_put2d_r8 (obj, tempo_var_column_amount, 0, ntimes, &
                        amf_corr_column (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, tempo_var_column_uncert, 0, ntimes, &
                        amf_corr_column_uncertainty (1:nxtrack, 0:ntimes-1), errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, "in write_amf_correction", errstat)
      return
    endif

  end subroutine write_amf_correction

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
