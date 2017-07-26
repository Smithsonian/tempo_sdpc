!> NetCDF output functions
!! @file
!! @note Ideally, \a use statements in this module will provide
!!       access only to subroutines, type definitions, and compile-time
!!       constants. All variables should be passed via the subroutine
!!       parameter argument lists.

module output_tools
  use netcdf
  use tell_module
  use tio_module
  use tg_names_module
  use OMSAO_precision_module
  use ctrlvars, only: yn_diagnostic_run, yn_refseccor, yn_scat_weights

  implicit none
  private

  public create_output_file, close_output_file, write_wavcal_output, &
    write_radfit_output, write_fitting_statistics, write_common_mode, &
    write_albedo, write_gas_profile, write_scattering_weights, &
    write_amf_correction, write_refspec_database, &
    write_reference_sector_corrected_column, &
    write_solar_wavecal_diagnostics, &
    write_radiance_wavecal_diagnostics, copy_pixel_corners, &
    copy_metadata, label_output_file, &
    read_geofields, read_column_results, read_cloud_params

  type (tiof_file_type), private, save, target :: primary_output_file

  ! using fill values from the original code simplifies diffing output files
  real (kind=8), private, parameter :: &
    fill_short = -30000, &
    fill_float = -1.0e30, &
    fill_double = -1.0e30_r8

contains

  subroutine write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(in) :: num_steps, num_xtrack
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(num_xtrack) :: xtrack_indices
    integer, dimension(num_steps) :: step_indices
    integer :: i, dimids(2)

    if (errstat /= 0) return

    call tiof_dimlist_lookup (dimlist, [tg_dim_xtrack, tg_dim_step], dimids, errstat)

    ! netcdf coordinate variables:
    call tiof_varlist_append (varlist, errstat, tg_dim_xtrack, nf90_int, &
                             dimids=[dimids(1)])
    call tiof_varlist_append (varlist, errstat, tg_dim_step, nf90_int, &
                             dimids=[dimids(2)])
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)

    ! mirror step indices may be overwritten later by copying from the input
    ! radiance file.
    step_indices = [(i, i=0,num_steps-1)]
    call tiof_put1d_i4 (obj, tg_dim_step, [0], [num_steps], step_indices, errstat)

    xtrack_indices = [(i, i=0,num_xtrack-1)]
    call tiof_put1d_i4 (obj, tg_dim_xtrack, [0], [num_xtrack], xtrack_indices, errstat)

  end subroutine write_coordinate_vars

  subroutine append_common_mode_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(2) :: dimids_commwvl_xtrack, dimids_xtrack_pair

    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_commwvl, tg_dim_xtrack], &
                              dimids_commwvl_xtrack, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_xtrack, tg_dim_pair], &
                              dimids_xtrack_pair, &
                              errstat)

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_common_mode_spectrum, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack,  &
                              comment = "common mode spectrum", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_common_mode_wavelengths, &
                              nf90_float, &
                              dimids = dimids_commwvl_xtrack,  &
                              comment = "common mode wavelengths", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_common_mode_count, &
                              nf90_int, &
                              dimids = [dimids_commwvl_xtrack(2)],  &
                              comment = "common mode spectrum averaging count", &
                              valid_range = [0.0_r8, 2147483647.0_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_common_mode_ccd_pixel_range, &
                              nf90_short, &
                              dimids = dimids_xtrack_pair,  &
                              comment = "first and last ccd pixel number fitted", &
                              valid_range = [0.0_r8, 1024.0_r8])

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)

  end subroutine append_common_mode_vars

  subroutine append_amf_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_xtrack_step_levels, dimsizes_xtrack_step_levels
    integer, dimension(3) :: chunksizes
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat /= 0) return

    ! lookup dimids for relevant array shapes
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_xtrack, tg_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_xtrack, tg_dim_step,tg_dim_swt_level], &
                              dimids_xtrack_step_levels, &
                              errstat, dimsizes = dimsizes_xtrack_step_levels)

    ! append amf variables
    chunksizes(1) = min(dimsizes_xtrack_step_levels(1), 128)  ! xtrack dimension
    chunksizes(2) = min(dimsizes_xtrack_step_levels(2), 128)  ! step dimension
    chunksizes(3) = dimsizes_xtrack_step_levels(3)            ! level dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_scattering_weights, &
                              nf90_double, &
                              dimids = dimids_xtrack_step_levels,  &
                              comment = "scattering weights", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_gas_profile, &
                              nf90_double, &
                              dimids = dimids_xtrack_step_levels,  &
                              comment = "gas profile", &
                              units = "ppb", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_albedo, &
                              nf90_double, &
                              dimids = dimids_xtrack_step,  &
                              comment = "albedo", &
                              valid_range = [-1e30_r8, 1e30_r8])

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_molecule_specific, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "molecule-specific air mass factor (AMF)", &
                              valid_range = [0.0_r8, 1e30_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_diagnostic_flag, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "diagnostic flag for molecule-specific air mass factor (AMF)", &
                              valid_range = [-2.0_r8, 13127.0_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_geometric, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "geometric air mass factor (AMF)", &
                              valid_range = [0.0_r8, 1e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_cloud_fraction, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "adjusted cloud fraction for AMF computation", &
                              valid_range = [0.0_r8, 1e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_cloud_pressure, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "adjusted cloud pressure for AMF computation", &
                              valid_range = [0.0_r8, 1e30_r8])

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)

  end subroutine append_amf_vars

  subroutine append_diagnostic_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(2) :: dimids_xtrack_step, dimids_refwavl_xtrack
    integer, dimension(3) :: dimids_var_xtrack_step, dimsizes_var_xtrack_step
    integer, dimension(3) :: dimids_commwvl_xtrack_step, dimsizes_commwvl_xtrack_step
    integer, dimension(3) :: dimids_refwavl_xtrack_refspec, dimsizes_refwavl_xtrack_refspec
    integer, dimension(3) :: chunksizes
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat /= 0) return

    ! lookup dimids for relevant array shapes
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_xtrack, tg_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_fitvar, tg_dim_xtrack, tg_dim_step], &
                              dimids_var_xtrack_step, &
                              errstat, dimsizes = dimsizes_var_xtrack_step)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_commwvl, tg_dim_xtrack, tg_dim_step], &
                              dimids_commwvl_xtrack_step, &
                              errstat, dimsizes = dimsizes_commwvl_xtrack_step)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_refwavl, tg_dim_xtrack, tg_dim_refspec], &
                              dimids_refwavl_xtrack_refspec, &
                              errstat, dimsizes = dimsizes_refwavl_xtrack_refspec)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_refwavl, tg_dim_xtrack], &
                              dimids_refwavl_xtrack, &
                              errstat)

    ! append diagnostic variables
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_iteration_count, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "radiance fit iteration count", &
                              valid_range = [0.0_r8, 32767.0_r8], &
                              fillvalue = fill_short)

    chunksizes(1) = dimsizes_var_xtrack_step(1)  ! var dimension
    chunksizes(2) = dimsizes_var_xtrack_step(2)  ! xtrack dimension
    chunksizes(3) = 1                            ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_params, &
                              nf90_double, &
                              dimids = dimids_var_xtrack_step, &
                              comment = "fit parameter", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_errors, &
                              nf90_double, &
                              dimids = dimids_var_xtrack_step, &
                              comment = "fit parameter uncertainty", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_correl, &
                              nf90_double, &
                              dimids = dimids_var_xtrack_step, &
                              comment = "fit parameter correlation", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_commwvl_xtrack_step(1)            ! wavelength dimension
    chunksizes(2) = min(dimsizes_commwvl_xtrack_step(2),1024)  ! xtrack dimension
    chunksizes(3) = 1                                          ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_measured_spectrum, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "measured spectrum", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_measured_wavelengths, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "measured wavelength", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_model_spectrum, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "model spectrum", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_weights, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "spectrum fit weights", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_residuals, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              comment = "spectrum fit residuals", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_param_names, &
                              nf90_string, &
                              dimids = [dimids_var_xtrack_step(1)])

    chunksizes(1) = dimsizes_refwavl_xtrack_refspec(1)              ! wavelength dimension
    chunksizes(2) = min(dimsizes_refwavl_xtrack_refspec(2), 1024)   ! xtrack dimension
    chunksizes(3) = 1                                               ! refspec dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_refspec, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack_refspec, &
                              comment = "reference spectra used in the fitting process", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_refspec_wavelength, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack, &
                              comment = "reference spectra wavelengths", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_refspec_norm, &
                              nf90_double, &
                              dimids = [dimids_refwavl_xtrack_refspec(3) ], &
                              comment = "reference spectra normalization factors", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)

  end subroutine append_diagnostic_vars

  subroutine append_column_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist, varlist_geo, varlist_qa
    type (tiof_attlist_type) :: att_coord, att_latbnd, att_lonbnd
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_corner_xtrack_step

    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_xtrack, tg_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_corner, tg_dim_xtrack, tg_dim_step], &
                              dimids_corner_xtrack_step, &
                              errstat)

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    ! data field variables with optional attribute lists:
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(tg_var_longitude) &
                              //' '//trim(tg_var_latitude))
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_column_amount, &
                              nf90_double, &
                              dimids = dimids_xtrack_step,  &
                              comment = "column amount", &
                              units = "molec/cm2", &
                              valid_range = [-1.e30_r8, 1.e30_r8], &
                              fillvalue = fill_double, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_column_uncert, &
                              nf90_double, &
                              dimids = dimids_xtrack_step,  &
                              comment = "column amount uncertainty", &
                              units = "molec/cm2", &
                              valid_range = [-1.e30_r8, 1.e30_r8], &
                              fillvalue = fill_double, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_main_dqf, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "main data quality flag", &
                              valid_range = [-1.0_r8, 2.0_r8])
    if (yn_refseccor) then
      call tiof_varlist_append (varlist, errstat, &
                                tg_var_refseccor_vertical_column, &
                                nf90_double, &
                                dimids = dimids_xtrack_step, &
                                comment = "reference sector corrected vertical_column", &
                                units = "molec/cm2", &
                                valid_range = [-1e30_r8, 1e30_r8], &
                                fillvalue = fill_double)
    endif
    call tiof_push_group (obj, tg_grp_product, errstat)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist)

    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_radfit_rms_residual, &
                              nf90_double, &
                              dimids = dimids_xtrack_step,  &
                              comment = "fit rms residual", &
                              valid_range = [0.0_r8, 1.e30_r8], &
                              fillvalue = fill_double)
    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_radfit_convergence_flag, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "fit convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8], &
                              fillvalue = fill_short)

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_def_vars (obj, varlist_qa, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist_qa)

    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_time, &
                              nf90_double, &
                              dimids = [dimids_xtrack_step(2)],  &
                              comment = "exposure start time", &
                              units = "s", &
                              valid_range = [0.0_r8, 1.e30_r8], &
                              fillvalue = fill_double)

    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = tg_var_latitude_bounds)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_latitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "latitude at pixel center", &
                              units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_latbnd)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_latitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              comment = "latitude at pixel corners (SW,SE,NE,NW)", &
                              units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              fillvalue = fill_float)

    call tiof_attlist_append (att_lonbnd, errstat, "bounds", &
                              att_text = tg_var_longitude_bounds)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_longitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "longitude at pixel center", &
                              units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_lonbnd)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_longitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              comment = "longitude at pixel corners (SW,SE,NE,NW)", &
                              units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float)

    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_sz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_sa_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_vz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_va_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_relative_azimuth, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "relative azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_terrain_height, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "terrain height", &
                              units = "m", &
                              valid_range = [-1000.0_r8, 10000.0_r8], &
                              fillvalue = fill_short)

    call tiof_push_group (obj, tg_grp_geolocation, errstat)
    call tiof_def_vars (obj, varlist_geo, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist_geo)

    call tiof_attlist_free (att_coord)
    call tiof_attlist_free (att_latbnd)
    call tiof_attlist_free (att_lonbnd)

  end subroutine append_column_vars

  subroutine append_wavcal_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist, varlist_diag, varlist_qa
    integer, dimension(1) :: dimid_xtrack
    integer, dimension(2) :: dimids_refwavl_xtrack

    call tiof_dimlist_lookup (dimlist, [tg_dim_xtrack], dimid_xtrack, errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [tg_dim_refwavl, tg_dim_xtrack], &
                              dimids_refwavl_xtrack, &
                              errstat)

    if (yn_diagnostic_run) then
      call tiof_varlist_append (varlist_diag, errstat, &
                                tg_var_solcal_wavelengths, &
                                nf90_double, &
                                dimids = dimids_refwavl_xtrack,  &
                                comment = "solar wavecal wavelengths", &
                                valid_range = [100.0_r8, 1000.0_r8], &
                                fillvalue = fill_double)
      call tiof_varlist_append (varlist_diag, errstat, &
                                tg_var_solcal_residuals, &
                                nf90_double, &
                                dimids = dimids_refwavl_xtrack,  &
                                comment = "solar wavecal residuals", &
                                valid_range = [-1e30_r8, 1e30_r8], &
                                fillvalue = fill_double)
      call tiof_varlist_append (varlist_diag, errstat, &
                                tg_var_radcal_wavelengths, &
                                nf90_double, &
                                dimids = dimids_refwavl_xtrack,  &
                                comment = "radiance wavecal wavelengths", &
                                valid_range = [100.0_r8, 1000.0_r8], &
                                fillvalue = fill_double)
      call tiof_varlist_append (varlist_diag, errstat, &
                                tg_var_radcal_residuals, &
                                nf90_double, &
                                dimids = dimids_refwavl_xtrack,  &
                                comment = "radiance wavecal residuals", &
                                valid_range = [-1e30_r8, 1e30_r8], &
                                fillvalue = fill_double)

      call tiof_push_group (obj, tg_grp_diagnostic, errstat)
      call tiof_def_vars (obj, varlist_diag, errstat)
      call tiof_pop_group (obj, errstat)
      call tiof_varlist_free (varlist_diag)
      if (errstat /= 0) return
    endif

    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_solcal_convergence_flag, &
                              nf90_short, &
                              dimids = dimid_xtrack,  &
                              comment = "solar wavelength calibration convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8], &
                              fillvalue = fill_short)
    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_radcal_convergence_flag, &
                              nf90_short, &
                              dimids = dimid_xtrack,  &
                              comment = "radiance wavelength calibration convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8], &
                              fillvalue = fill_short)
    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_radref_convergence_flag, &
                              nf90_short, &
                              dimids = dimid_xtrack,  &
                              comment = "radiance reference fit convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8], &
                              fillvalue = fill_short)
    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_radref_fit_rms, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              comment = "radiance reference fit RMS", &
                              valid_range = [0.0_r8, 1.e30_r8])

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_def_vars (obj, varlist_qa, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist_qa)

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radref_column_amount, &
                              nf90_double, &
                              dimids = dimid_xtrack, &
                              comment = "radiance reference fit column amount", &
                              units = "molec/cm2", &
                              valid_range = [-1.e30_r8, 1.e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radref_column_uncert, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              comment = "radiance reference fit column uncert", &
                              units = "molec/cm2", &
                              valid_range = [0.0_r8, 1.e30_r8])
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radref_column_xtrfit, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              comment = "radiance reference fit column XTR fit", &
                              units = "molec/cm2", &
                              valid_range = [0.0_r8, 1.e30_r8])

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist)

  end subroutine append_wavcal_vars

  subroutine append_fit_stats (obj, errstat)
    implicit none

    type (tiof_file_type), intent(in) :: obj
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist

    if (errstat /= 0) return

    call tiof_varlist_append (varlist, errstat, "num_crosstrack_pixels", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_scan_lines", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_good_input", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_good_output", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_suspect_output", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_bad_output", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_converged", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_failed_convergence", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_exceeded_iterations", nf90_int)
    call tiof_varlist_append (varlist, errstat, "num_out_of_bounds", nf90_int)
    call tiof_varlist_append (varlist, errstat, "percent_good_output", nf90_float)
    call tiof_varlist_append (varlist, errstat, "percent_bad_output", nf90_float)
    call tiof_varlist_append (varlist, errstat, "percent_suspect_output", nf90_float)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)
  end subroutine append_fit_stats

  !> Create netCDF format Level 2 product file
  !! @param[in] filename   netCDF output file name
  !! @param[in] num_steps  Number of scan steps
  !! @param[in] num_xtrack  Number of cross-track pixels
  !! @param[in] num_swlevels  Number of height levels used in AMF climatologies
  !! @param[in] n_comm_wvl  Number of common-mode wavelengths
  !! @param[in] nwavel_max  Maximum number of reference spectrum wavelengths
  !! @param[in] max_rs_idx  Maximum reference spectrum index
  !! @param[in] n_fitvar_rad  Number of radiance spectrum fit variables
  !! @param[inout]  errstat  Error status variable
  subroutine create_output_file (filename, num_steps, num_xtrack, num_swlevels, &
                                 n_comm_wvl, nwavel_max, max_rs_idx, n_fitvar_rad, &
                                 errstat)
    implicit none
    character (len=*), intent(in) :: filename
    integer (kind=i4), intent(in) :: num_steps, num_xtrack, num_swlevels, &
      n_comm_wvl, nwavel_max, max_rs_idx, n_fitvar_rad
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_dimlist_type) :: dimlist

    if (errstat /= 0) return

    obj => primary_output_file

    ! Create a file.
    call tiof_create (obj, filename, nf90_clobber, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: creating file "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_put_git_commit_hash (obj, errstat)

    ! Create default groups.
    call tiof_def_group (obj, tg_grp_product, errstat)
    call tiof_def_group (obj, tg_grp_geolocation, errstat)
    call tiof_def_group (obj, tg_grp_support_data, errstat)
    call tiof_def_group (obj, tg_grp_qa_stats, errstat)
    call tiof_def_group (obj, tg_grp_metadata, errstat)
    if (yn_diagnostic_run) then
      call tiof_def_group (obj, tg_grp_diagnostic, errstat)
    endif
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file:  defining groups in "//trim(filename), &
                       errstat)
      return
    endif

    ! Define a dimension list.
    call tiof_dimlist_append (dimlist, tg_dim_step, num_steps, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_xtrack, num_xtrack, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_corner, 4, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_swt_level, num_swlevels, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_pair, 2, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_commwvl, n_comm_wvl, errstat)
    if (yn_diagnostic_run) then
      ! For simplicity, use one dimlist for all groups
      call tiof_dimlist_append (dimlist, tg_dim_fitvar, n_fitvar_rad, errstat)
      call tiof_dimlist_append (dimlist, tg_dim_refwavl, nwavel_max, errstat)
      call tiof_dimlist_append (dimlist, tg_dim_refspec, max_rs_idx, errstat)
    endif
    call tiof_def_dims (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: defining dimensions in "//trim(filename), &
                       errstat)
      return
    endif

    call write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: writing coordinate variables to "//trim(filename), &
                       errstat)
      return
    endif

    ! Define variables roughly in order of importance to aid users
    ! viewing the file with an application like ncdump
    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call append_fit_stats (obj, errstat)
    call tiof_pop_group (obj, errstat)

    call append_column_vars (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: defining variables in "//trim(filename), &
                       errstat)
      return
    endif

    if (yn_scat_weights) then
      call tiof_push_group (obj, tg_grp_support_data, errstat)
      call append_amf_vars (obj, dimlist, errstat)
      call tiof_pop_group (obj, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_write_error, &
                         "create_output_file: defining amf variables in "//trim(filename), &
                         errstat)
        return
      endif
    endif

    call append_wavcal_vars (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: defining variables in "//trim(filename), &
                       errstat)
      return
    endif

    if (yn_diagnostic_run) then
      call tiof_push_group (obj, tg_grp_diagnostic, errstat)
      call append_common_mode_vars (obj, dimlist, errstat)
      call append_diagnostic_vars (obj, dimlist, errstat)
      call tiof_pop_group (obj, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_write_error, &
                         "create_output_file: defining diagnostic variables in "//trim(filename), &
                         errstat)
        return
      endif
    endif

    call tiof_dimlist_free (dimlist)
  end subroutine create_output_file

  !> Write a block of radiance fit results to Level 2 product file
  !! @param[in] iline   Starting scan line index of block
  !! @param[in] nblock  Number of scan lines in block
  !! @param[in] nxtrack  Number of cross-track pixels
  !! @param[in] n_fitvar_rad  Number of fit variables in radiance spectrum model
  !! @param[in] n_rad_wvl  Number of wavelengths in fitted spectral range
  !! @param[in] input_vars  Structure containing variables read from input files
  !!                        and copied to the output file unmodified.
  !!                        For now, it's mostly geolocation variables.
  !! @param[in] result_vars  Structure containing radiance fit results.
  !! @param[in] radfit_diagnostics Structure containing radiance fit diagnostics.
  !! @param[inout] errstat  Error status variable
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

    type (tiof_file_type), pointer :: obj
    real (kind=r8), dimension(1:n_rad_wvl, 1:nxtrack, 0:nblock-1) :: residuals
    real (kind=r8), dimension(:,:,:), pointer :: waves, meas, model, weights
    real (kind=r4), dimension(1:nxtrack,0:nblock-1) :: relative_azimuth
    integer :: i,j
    real (kind=r4) :: razi, sazi, vazi

    if (errstat /= 0) return

    obj => primary_output_file

    ! result_vars
    call tiof_push_group (obj, tg_grp_product, errstat)
    call tiof_put2d_r8 (obj, tg_var_column_amount, [iline,0], [nblock, nxtrack], &
                        result_vars % column_amount (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r8 (obj, tg_var_column_uncert, [iline,0], [nblock, nxtrack], &
                        result_vars % column_uncert (1:nxtrack, 0:nblock-1), errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_put2d_r8 (obj, tg_var_radfit_rms_residual, [iline,0], [nblock, nxtrack], &
                        result_vars % fit_rms_residual (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_i2 (obj, tg_var_radfit_convergence_flag, [iline,0], [nblock, nxtrack], &
                        result_vars % fit_convergence_flag (1:nxtrack, 0:nblock-1), errstat)
    call tiof_pop_group (obj, errstat)

    if (yn_diagnostic_run) then
      call tiof_push_group (obj, tg_grp_diagnostic, errstat)
      call tiof_put2d_i2 (obj, tg_var_radfit_iteration_count, [iline,0], [nblock, nxtrack], &
                          result_vars % fit_iteration_count (1:nxtrack, 0:nblock-1), errstat)

      call tiof_put3d_r8 (obj, tg_var_radfit_params, [iline,0,0], [nblock,nxtrack,n_fitvar_rad], &
                          radfit_diagnostics % params(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj, tg_var_radfit_errors, [iline,0,0], [nblock,nxtrack,n_fitvar_rad], &
                          radfit_diagnostics % errors(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj, tg_var_radfit_correl, [iline,0,0], [nblock,nxtrack,n_fitvar_rad], &
                          radfit_diagnostics % correl(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)

      model   => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 1, 0:nblock-1)
      meas    => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 2, 0:nblock-1)
      waves   => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 3, 0:nblock-1)
      weights => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 4, 0:nblock-1)

      call tiof_put3d_r8 (obj, tg_var_radfit_model_spectrum, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], model, errstat)
      call tiof_put3d_r8 (obj, tg_var_radfit_measured_spectrum, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], meas, errstat)
      call tiof_put3d_r8 (obj, tg_var_radfit_measured_wavelengths, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl],  waves, errstat)
      call tiof_put3d_r8 (obj, tg_var_radfit_weights, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], weights, errstat)

      residuals(:,:,:) = meas - model
      call tiof_put3d_r8 (obj, tg_var_radfit_residuals, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], residuals, errstat)
      call tiof_pop_group (obj, errstat)
    endif

    ! Compute relative azimuth angle
    relative_azimuth(:,:) = fill_float
    do j=0,nblock-1
      do i=1,nxtrack
        vazi = input_vars % viewing_azimuth(i,j)
        sazi = input_vars % solar_azimuth(i,j)
        if (abs(vazi) <= 360.0 .and. abs(sazi) <= 360.0) then
          razi = vazi - sazi
          if (razi > 180.0) then
            razi = razi - 360.0
          else if (razi < -180.0) then
            razi = razi + 360.0
          endif
          relative_azimuth(i,j) = razi
        endif
      enddo
    enddo

    ! input_vars
    call tiof_push_group (obj, tg_grp_geolocation, errstat)
    call tiof_put1d_r8 (obj, tg_var_time, [iline], [nblock], &
                        input_vars % time (0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tg_var_longitude, [iline,0], [nblock,nxtrack], &
                        input_vars % longitude (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tg_var_latitude, [iline,0], [nblock,nxtrack], &
                        input_vars % latitude (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tg_var_sz_angle, [iline,0], [nblock,nxtrack], &
                        input_vars % solar_zenith (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tg_var_sa_angle, [iline,0], [nblock,nxtrack], &
                        input_vars % solar_azimuth (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tg_var_vz_angle, [iline,0], [nblock,nxtrack], &
                        input_vars % viewing_zenith (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tg_var_va_angle, [iline,0], [nblock,nxtrack], &
                        input_vars % viewing_azimuth (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r4 (obj, tg_var_relative_azimuth, [iline,0], [nblock,nxtrack], &
                        relative_azimuth (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_i2 (obj, tg_var_terrain_height, [iline,0], [nblock,nxtrack], &
                        input_vars % terrain_height (1:nxtrack, 0:nblock-1), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_radfit_output: failed", errstat)
      return
    endif
  end subroutine write_radfit_output

  !> Write wavelength calibration to Level 2 product file
  !! @param[in] result_vars Structure containing radiance fit results.
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[inout] errstat  Error status variable
  subroutine write_wavcal_output (result_vars, nxtrack, errstat)
    use OMSAO_omidata_module, only : result_vars_type
    implicit none
    type(result_vars_type), intent(in) :: result_vars
    integer, intent(in) :: nxtrack
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_put1d_i2 (obj, tg_var_solcal_convergence_flag, [0], [nxtrack], &
                        result_vars % solcal_convergence_flag (1:nxtrack), errstat)
    call tiof_put1d_i2 (obj, tg_var_radcal_convergence_flag, [0], [nxtrack], &
                        result_vars % radcal_convergence_flag (1:nxtrack), errstat)
    call tiof_put1d_i2 (obj, tg_var_radref_convergence_flag, [0], [nxtrack], &
                        result_vars % radref_convergence_flag (1:nxtrack), errstat)
    call tiof_put1d_r8 (obj, tg_var_radref_fit_rms, [0], [nxtrack], &
                        result_vars % radref_fit_rms (1:nxtrack), errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put1d_r8 (obj, tg_var_radref_column_amount, [0], [nxtrack], &
                        result_vars % radref_column_amount (1:nxtrack), errstat)
    call tiof_put1d_r8 (obj, tg_var_radref_column_uncert, [0], [nxtrack], &
                        result_vars % radref_column_uncert (1:nxtrack), errstat)
    call tiof_put1d_r8 (obj, tg_var_radref_column_xtrfit, [0], [nxtrack], &
                        result_vars % radref_column_xtrfit (1:nxtrack), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_wavcal_output: failed", errstat)
      return
    endif
  end subroutine write_wavcal_output

  !> Write radiance fit QA statistics to Level 2 product file
  !! @param[in] stats  Structure containing radiance fit QA statistics
  !! @param[in] param_names  Array of fit parameter names
  !! @param[in] num_params Number of fit parameters
  !! @param[inout] errstat  Error status variable
  subroutine write_fitting_statistics (stats, param_names, num_params, errstat)
    use omi_pge_fitting_aux, only : fitting_statistics_type
    implicit none

    type (fitting_statistics_type), intent(in) :: stats
    character (len=*), dimension (:), intent(in) :: param_names
    integer, intent(in) :: num_params
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    integer :: ncp, nsl

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_put_i4 (obj, "num_crosstrack_pixels", stats % num_crosstrack_pixels, errstat)
    call tiof_put_i4 (obj, "num_scan_lines", stats % num_scan_lines, errstat)
    call tiof_put_i4 (obj, "num_good_input", stats % num_good_input, errstat)
    call tiof_put_i4 (obj, "num_good_output", stats % num_good_output, errstat)
    call tiof_put_i4 (obj, "num_suspect_output", stats % num_suspect_output, errstat)
    call tiof_put_i4 (obj, "num_bad_output", stats % num_bad_output, errstat)
    call tiof_put_i4 (obj, "num_converged", stats % num_converged, errstat)
    call tiof_put_i4 (obj, "num_failed_convergence", stats % num_failed_convergence, errstat)
    call tiof_put_i4 (obj, "num_exceeded_iterations", stats % num_exceeded_iterations, errstat)
    call tiof_put_i4 (obj, "num_out_of_bounds", stats % num_out_of_bounds, errstat)
    call tiof_put_r4 (obj, "percent_good_output", stats % percent_good_output, errstat)
    call tiof_put_r4 (obj, "percent_bad_output", stats % percent_bad_output, errstat)
    call tiof_put_r4 (obj, "percent_suspect_output", stats % percent_suspect_output, errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "write_fitting_statistics: writing fitting statistics", &
                       errstat)
      return
    endif

    ncp = stats % num_crosstrack_pixels
    nsl = stats % num_scan_lines

    call tiof_push_group (obj, tg_grp_product, errstat)
    call tiof_put2d_i2 (obj, tg_var_main_dqf, [0,0], [nsl,ncp], &
                        stats % quality_flag (1:ncp, 0:nsl-1), &
                        errstat)
    call tiof_pop_group (obj, errstat)

    if (yn_diagnostic_run) then
      call tiof_push_group (obj, tg_grp_diagnostic, errstat)
      call tiof_put1d_string (obj, tg_var_radfit_param_names, 0, num_params, &
                              param_names(1:num_params), errstat)
      call tiof_pop_group (obj, errstat)
    endif

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "write_fitting_statistics: writing fitting statistics", &
                       errstat)
      return
    endif
  end subroutine write_fitting_statistics

  !> Write AMF albedo to Level 2 product file
  !! @param[in] albedo  Surface albedo used in AMF calculation
  !! @param[in] nxtrack  Number of cross-track pixels
  !! @param[in] ntimes  Number of scans
  !! @param[inout] errstat  Error status variable
  subroutine write_albedo (albedo, nxtrack, ntimes, errstat)
    implicit none

    integer, intent(in) :: nxtrack, ntimes
    real (kind=r8), dimension (1:nxtrack, 0:ntimes-1), intent(in) :: albedo
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put2d_r8 (obj, tg_var_amf_albedo, [0,0], [ntimes,nxtrack], &
                        albedo (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_albedo", errstat)
      return
    endif
  end subroutine write_albedo

  !> Write AMF vertical gas profile climatology to Level 2 product file
  !! @param[in] gas_profile  Column density vertical profile climatology for each pixel
  !!                     [molecules/cm^2].
  !! @param[in] climatology_levels  Altitude [km] coordinate for each gas profile point.
  !! @param[in] nxtrack  Number of cross-track pixels
  !! @param[in] ntimes  Number of scans
  !! @param[in] nlevels  Number of altitudes in vertical profile climatology
  !! @param[inout] errstat  Error status variable
  subroutine write_gas_profile (gas_profile, &
                                nxtrack, ntimes, nlevels, errstat)
    implicit none

    integer, intent(in) :: nxtrack, ntimes, nlevels
    real (kind=r8), dimension (1:nxtrack, 0:ntimes-1, 1:nlevels), intent(in) :: gas_profile
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put3d_r8 (obj, tg_var_amf_gas_profile, [0,0,0], [nlevels,ntimes,nxtrack], &
                        gas_profile(1:nxtrack, 0:ntimes-1, 1:nlevels), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_gas_profile", errstat)
      return
    endif
  end subroutine write_gas_profile

  !> Write AMF scattering weights to Level 2 product file
  !! @param[in] scattw Scattering weights
  !! @param nxtrack  Number of cross-track pixels
  !! @param ntimes  Number of scans
  !! @param nlevels  Number of altitudes in vertical profile climatology
  !! @param errstat  Error status variable
  subroutine write_scattering_weights (scattw, nxtrack, ntimes, nlevels, errstat)
    implicit none

    integer, intent(in) :: nxtrack, ntimes, nlevels
    real (kind=r8), dimension (1:nxtrack, 0:ntimes-1, 1:nlevels), intent(in) :: scattw
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put3d_r8 (obj, tg_var_amf_scattering_weights, [0,0,0], [nlevels,ntimes,nxtrack], &
                        scattw (1:nxtrack, 0:ntimes-1, 1:nlevels), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_scattering_weights", errstat)
      return
    endif
  end subroutine write_scattering_weights

  !> Write AMF correction to Level 2 product file
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[in] ntimes Number of scans
  !! @param[in] amf_corr Air mass factor (AMF) correction
  !! @param[in] amf_corr_column  AMF-corrected column density
  !! @param[in] amf_corr_column_uncertainty  Uncertainty in AMF-corrected column density
  !! @param[in] yn_write_cloud_variables   If \a .true., write cloud variables
  !!                                      to Level 2 product file
  !! @param[inout] errstat  Error status variable
  subroutine write_amf_correction (nxtrack, ntimes, amf_corr, &
                                   amf_corr_column, amf_corr_column_uncertainty, &
                                   yn_write_cloud_variables, errstat)
    use OMSAO_omidata_module, only : amf_correction_type
    implicit none

    integer, intent(in) :: nxtrack, ntimes
    type (amf_correction_type), intent(in) :: amf_corr
    real (kind=r8), dimension(1:nxtrack,0:ntimes-1), intent(in) :: amf_corr_column
    real (kind=r8), dimension(1:nxtrack,0:ntimes-1), intent(in) :: amf_corr_column_uncertainty
    logical, intent(in) :: yn_write_cloud_variables
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put2d_i2 (obj, tg_var_amf_diagnostic_flag, [0,0], [ntimes,nxtrack], &
                        amf_corr % diagnostic_flag (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, tg_var_amf_geometric, [0,0], [ntimes,nxtrack], &
                        amf_corr % amf_geometric (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, tg_var_amf_molecule_specific, [0,0], [ntimes,nxtrack], &
                        amf_corr % amf_molecule_specific (1:nxtrack, 0:ntimes-1), errstat)

    if (yn_write_cloud_variables) then
      call tiof_put2d_r8 (obj, tg_var_amf_cloud_fraction, [0,0], [ntimes,nxtrack], &
                          amf_corr % cloud_fraction (1:nxtrack, 0:ntimes-1), errstat)
      call tiof_put2d_r8 (obj, tg_var_amf_cloud_pressure, [0,0], [ntimes,nxtrack], &
                          amf_corr % cloud_pressure (1:nxtrack, 0:ntimes-1), errstat)
    endif
    call tiof_pop_group (obj, errstat)

    ! Note that we're over-writing the column amount variable in the file
    ! (to which we previously wrote the slant column values).
    call tiof_push_group (obj, tg_grp_product, errstat)
    call tiof_put2d_r8 (obj, tg_var_column_amount, [0,0], [ntimes,nxtrack], &
                        amf_corr_column (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, tg_var_column_uncert, [0,0], [ntimes,nxtrack], &
                        amf_corr_column_uncertainty (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_amf_correction", errstat)
      return
    endif

  end subroutine write_amf_correction

  !> Write common-mode spectral to Level 2 product file
  !! @param[in] nxtrack  Number of cross-track pixels
  !! @param[in] n_comm_wvl  Number of common-mode wavelengths
  !! @param[in] common_mode  Structure containing common-mode spectra.
  !! @param[inout]  errstat  Error status variable
  !! @details
  !! Common-mode spectra are essentially average fit residuals.
  subroutine write_common_mode (nxtrack, n_comm_wvl, common_mode, errstat)
    use OMSAO_variables_module, only : common_mode_spectrum_type
    implicit none

    integer (kind=i4), intent(in) :: nxtrack, n_comm_wvl
    type (common_mode_spectrum_type), intent(in) :: common_mode
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_diagnostic, errstat)
    call tiof_put2d_r8 (obj, tg_var_common_mode_spectrum, [0,0], [nxtrack,n_comm_wvl], &
                        common_mode % refspecdata (1:n_comm_wvl,1:nxtrack), errstat)
    call tiof_put2d_r8 (obj, tg_var_common_mode_wavelengths, [0,0], [nxtrack,n_comm_wvl], &
                        common_mode % refspecwavs (1:n_comm_wvl,1:nxtrack), errstat)
    call tiof_put2d_i2 (obj, tg_var_common_mode_ccd_pixel_range, [0,0], [2,nxtrack], &
                        common_mode % ccdpixel (1:nxtrack,1:2), errstat)
    call tiof_put1d_i4 (obj, tg_var_common_mode_count, [0], [nxtrack], &
                        common_mode % refspeccount (1:nxtrack), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_common_mode", errstat)
      return
    endif

  end subroutine write_common_mode

  !> Write reference spectra to Level 2 product file
  !! @param[in] db  Array of reference spectrum values
  !! @param[in] db_wvl  Array of reference spectrum wavelength grids
  !! @param[in] refspec  Structure containing reference spectra on the
  !!                     original wavelength grid, with normalization
  !!                     factors during the fit.
  !! @param[in] nrefspec Number of reference spectra
  !! @param[in] npts  Number of wavelength points
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[inout]  errstat  Error status variable
  subroutine write_refspec_database (db, db_wvl, refspec, &
                                     nrefspec, npts, nxtrack, errstat)
    use OMSAO_variables_module, only : reference_spectrum_type
    implicit none
    real (kind=r8), intent(in), dimension (:,:,:) :: db
    real (kind=r8), intent(in), dimension (:,:) :: db_wvl
    type(reference_spectrum_type), intent(in), dimension(:) :: refspec
    integer (kind=i4), intent(in) :: nrefspec, npts, nxtrack
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    integer :: i

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_diagnostic, errstat)

    ! Loop avoids creation of temporary array that may exceed process address space,
    ! causing a segv.  This can happen, when nxtrack is large, e.g. 2048.
    do i=1,nrefspec
      call tiof_put3d_r8 (obj, tg_var_refspec, [i-1,0,0], [1, nxtrack, npts], &
                          db(1:npts, 1:nxtrack, i:i), errstat)
    enddo

    call tiof_put2d_r8 (obj, tg_var_refspec_wavelength, [0,0], [nxtrack, npts], &
                        db_wvl(1:npts, 1:nxtrack), errstat)
    call tiof_put1d_r8 (obj, tg_var_refspec_norm, [0], [nrefspec], &
                        refspec (1:nrefspec) % normfactor, errstat)

    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_refspec_database", errstat)
      return
    endif

  end subroutine write_refspec_database

  !> Write reference sector corrected column to Level 2 product file
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[in] ntimes Number of scans
  !! @param[in] column Reference sector corrected column density
  !! @param[inout] errstat  Error status variable
  subroutine write_reference_sector_corrected_column (nxtrack, ntimes, column, errstat)
    implicit none
    integer (kind=i4), intent(in) :: nxtrack, ntimes
    real (kind=r8), dimension(1:nxtrack,0:ntimes-1), intent(in) :: column
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_product, errstat)
    call tiof_put2d_r8 (obj, tg_var_refseccor_vertical_column, [0,0], [ntimes, nxtrack], &
                        column(1:nxtrack, 0:ntimes-1), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_reference_sector_corrected_column", &
                       errstat)
      return
    endif

  end subroutine write_reference_sector_corrected_column

  !> Write solar wavelength calibration diagnostics to Level 2 product file
  !! @param[in] nwaves  Number of wavelength points
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[in] waves  Wavelength grids
  !! @param[in] resid  Residuals from wavelength calibration best-fit.
  !! @param[inout] errstat  Error status variable
  subroutine write_solar_wavecal_diagnostics (nwaves, nxtrack, waves, resid, &
                                              errstat)
    implicit none
    integer, intent(in) :: nwaves, nxtrack
    real (kind=r8), dimension(:,:), intent(in) :: waves, resid
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_diagnostic, errstat)
    call tiof_put2d_r8 (obj, tg_var_solcal_wavelengths, [0,0], [nxtrack, nwaves], &
                        waves (1:nwaves,1:nxtrack), errstat)
    call tiof_put2d_r8 (obj, tg_var_solcal_residuals, [0,0], [nxtrack, nwaves], &
                        resid (1:nwaves,1:nxtrack), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_solar_wavecal_diagnostics", errstat)
      return
    endif

  end subroutine write_solar_wavecal_diagnostics

  !> Write radiance wavelength calibration diagnostics to Level 2 product file
  !! @param[in] nwaves  Number of wavelength points
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[in] waves  Wavelength grids
  !! @param[in] resid  Residuals from wavelength calibration best-fit.
  !! @param[inout] errstat  Error status variable
  subroutine write_radiance_wavecal_diagnostics (nwaves, nxtrack, waves, resid, &
                                                 errstat)
    implicit none
    integer, intent(in) :: nwaves, nxtrack
    real (kind=r8), dimension(:,:), intent(in) :: waves, resid
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_diagnostic, errstat)
    call tiof_put2d_r8 (obj, tg_var_radcal_wavelengths, [0,0], [nxtrack, nwaves], &
                        waves (1:nwaves,1:nxtrack), errstat)
    call tiof_put2d_r8 (obj, tg_var_radcal_residuals, [0,0], [nxtrack, nwaves], &
                        resid (1:nwaves,1:nxtrack), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_radiance_wavecal_diagnostics", errstat)
      return
    endif

  end subroutine write_radiance_wavecal_diagnostics

  !> Close Level 2 product file
  !! @param[inout] errstat  Error status variable
  subroutine close_output_file (errstat)
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    obj => primary_output_file

    call tiof_close (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, "close_output_file failed", errstat)
    endif

  end subroutine close_output_file

  subroutine copy_pixel_corners (l1bfile, rad_group, ntimes, nxtrack, &
                                 corners_copied, errstat)
    implicit none
    character (len=*), intent(in) :: l1bfile, rad_group
    integer (kind=4), intent(in) :: ntimes, nxtrack
    logical, intent(out) :: corners_copied
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_file_type) :: l1b
    integer, dimension(ntimes) :: step_indices
    real (kind=4), dimension(4,1:nxtrack,1:ntimes) :: tmp

    if (errstat /= 0) return

    corners_copied = .false.

    call tiof_open (l1bfile, l1b, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "copy_pixel_corners: opening file "//trim(l1bfile), &
                       errstat)
      return
    endif

    ! Back-compatibility:
    ! If there's a problem reading mirror step indices,
    ! return errstat=0, corners_copied=.false.
    call tell_push_queue
    call tiof_get1d_i4 (l1b, tg_dim_step, [0], [ntimes], step_indices, errstat)
    if (errstat /= 0) then
      call tell_pop_queue (1)
      errstat = 0
      return
    endif
    call tell_pop_queue (0)

    ! use pixel corners associated with the radiances being fitted
    call tiof_push_group (l1b, rad_group, errstat)

    ! Back-compatibility:
    ! If there's a problem reading pixel corners,
    ! return errstat=0, corners_copied=.false.
    call tell_push_queue
    call tiof_get3d_r4 (l1b, tg_var_latitude_bounds, [0,0,0], [ntimes, nxtrack, 4], &
                        tmp(1:4,1:nxtrack,1:ntimes), errstat)
    if (errstat /= 0) then
      call tell_pop_queue (1)
      errstat = 0
      return
    endif
    call tell_pop_queue (0)

    obj => primary_output_file

    ! copy mirror step indices from input radiance file
    call tiof_push_group (obj, "/", errstat)
    call tiof_put1d_i4 (obj, tg_dim_step, [0], [ntimes], step_indices, errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_geolocation, errstat)
    call tiof_put3d_r4 (obj, tg_var_latitude_bounds, [0,0,0], [ntimes, nxtrack, 4], &
                        tmp(1:4,1:nxtrack,1:ntimes), errstat)
    call tiof_get3d_r4 (l1b, tg_var_longitude_bounds, [0,0,0], [ntimes, nxtrack, 4], &
                        tmp(1:4,1:nxtrack,1:ntimes), errstat)
    call tiof_put3d_r4 (obj, tg_var_longitude_bounds, [0,0,0], [ntimes, nxtrack, 4], &
                        tmp(1:4,1:nxtrack,1:ntimes), errstat)

    call tiof_close (l1b, errstat)
    if (errstat /= 0) then
      call tiof_pop_group (obj, errstat)
      call tell_error (tell_io_read_error, "copy_pixel_corners: reading file "//trim(l1bfile), &
                       errstat)
      return
    endif

    call tiof_pop_group (obj, errstat)

  end subroutine copy_pixel_corners

  subroutine copy_metadata (l1bfile, errstat)
    implicit none
    character (len=*), intent(in) :: l1bfile
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj
    type (tiof_file_type) :: l1b

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_open (l1bfile, l1b, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "copy_metadata: opening file "//trim(l1bfile), &
                       errstat)
      return
    endif

    call tiof_copy_granule_ident (l1b, obj, errstat)
    call tiof_close (l1b, errstat)

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "copy_metadata: copying from "//trim(l1bfile), &
                       errstat)
    endif
  end subroutine copy_metadata

  subroutine label_output_file (label, processing_version, errstat)
    implicit none
    character (len=*), intent(in) :: label
    integer, intent(in) :: processing_version
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file
    call tiof_label_product (obj, label, processing_version, errstat)
  end subroutine label_output_file

  !> Read geolocation fields
  !! @param[in] ntimes  Number of scans
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[inout] lat   Latitude [deg]
  !! @param[inout] lon   Longitude [deg]
  !! @param[inout] sza   Solar zenith angle [deg]
  !! @param[inout] vza   Viewing zenith angle [deg]
  !! @param[inout] thgt  Terrain height
  !! @param[inout] time  Time [s]
  !! @param[inout] errstat  Error status variable
  subroutine read_geofields (ntimes, nxtrack, lat, lon, sza, vza, thgt, time, errstat)
    use OMSAO_precision_module, only : i2, i4, r4
    use tio_module
    implicit none

    integer (kind=i4), intent(in) :: ntimes, nxtrack
    ! these arrays are (1:nxtrack,0:ntimes-1)
    real (kind=r4), dimension(:,:), intent(inout) :: lat, lon, sza, vza, thgt
    real (kind=r8), dimension(:), intent(inout) :: time
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    integer (kind=i2), dimension(nxtrack,ntimes) :: i2_thgt

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_geolocation, errstat)
    call tiof_get2d_r4 (obj, tg_var_latitude, [0,0], [ntimes, nxtrack], lat(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_longitude, [0,0], [ntimes, nxtrack], lon(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_sz_angle, [0,0], [ntimes, nxtrack], sza(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_vz_angle, [0,0], [ntimes, nxtrack], vza(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_i2 (obj, tg_var_terrain_height, [0,0], [ntimes, nxtrack], i2_thgt(1:nxtrack,1:ntimes), errstat)
    call tiof_get1d_r8 (obj, tg_var_time, [0], [ntimes], time(1:ntimes), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "in read_geofields", errstat)
      return
    endif

    thgt(1:nxtrack,1:ntimes) = real(i2_thgt(1:nxtrack,1:ntimes),kind=r4)

  end subroutine read_geofields

  !> Read fit results from Level 2 product file
  !! @param[in] ntimes  Number of scans
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[inout] col  Column density
  !! @param[inout] col_unc  Uncertainty in column density
  !! @param[inout] rms  Radiance fit RMS residual
  !! @param[inout] amf  Air mass factor (AMF)
  !! @param[inout] convergence_flag  Convergence flag
  !! @param[inout] errstat  Error status variable
  subroutine read_column_results (ntimes, nxtrack, col, col_unc, rms, amf, &
                                  convergence_flag, errstat)
    use OMSAO_precision_module, only : i2, r8
    use tio_module
    implicit none

    integer (kind=i4), intent(in) :: ntimes, nxtrack
    ! these arrays are (1:nxtrack,0:ntimes-1)
    real (kind=r8), dimension(:,:), intent(inout) :: col, col_unc, rms, amf
    integer (kind=i2), dimension(:,:) :: convergence_flag
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_product, errstat)
    call tiof_get2d_r8 (obj, tg_var_column_amount, [0,0], [ntimes, nxtrack], col(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r8 (obj, tg_var_column_uncert, [0,0], [ntimes, nxtrack], col_unc(1:nxtrack,1:ntimes), errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_get2d_r8 (obj, tg_var_radfit_rms_residual, [0,0], [ntimes, nxtrack], rms(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_i2 (obj, tg_var_radfit_convergence_flag, [0,0], [ntimes, nxtrack], &
                        convergence_flag(1:nxtrack,1:ntimes), errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_get2d_r8 (obj, tg_var_amf_molecule_specific, [0,0], [ntimes, nxtrack], amf(1:nxtrack,1:ntimes), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "in read_column_results", errstat)
      return
    endif

  end subroutine read_column_results

  !> Read cloud parameters from Level 2 cloud product file.
  !! @param[in] cloud_file  Name of Level 2 cloud product file
  !! @param[in] ntimes Number of scans
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[out] cloud_fraction  Cloud fraction
  !! @param[out] cloud_top_pressure  Cloud "top" pressure where the definition
  !!                               of "top" depends on the cloud product.
  !! @param[inout] errstat  Error status variable
  subroutine read_cloud_params (cloud_file, ntimes, nxtrack, cloud_fraction, &
                                cloud_top_pressure, errstat)
    use OMSAO_parameters_module, only : r8_missval
    implicit none
    character (len=*), intent(in) :: cloud_file
    integer (kind=i4), intent(in) :: ntimes, nxtrack
    real (kind=r8), dimension (1:nxtrack,0:ntimes-1), intent(out) :: cloud_fraction, cloud_top_pressure
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: cld
    character (len=128) :: logmsg

    if (errstat /= 0) return

    call tiof_open (cloud_file, cld, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "read_cloud_params: opening file "//trim(cloud_file), &
                       errstat)
      return
    endif

    write(logmsg,'(a,i6,i6,a)')'Reading cloud file, ntimes,nxtrack=', &
      ntimes, nxtrack, ' file = '//trim(cloud_file)
    call tell_log (1, logmsg)

    call tiof_push_group (cld, "/product", errstat)
    call tiof_get2d_r8 (cld, "cloud_pressure_for_O3", [0,0], [ntimes,nxtrack], &
                        cloud_top_pressure, errstat, replace_fill=r8_missval)
    call tiof_get2d_r8 (cld, "cloud_fraction_for_O3", [0,0], [ntimes,nxtrack], &
                        cloud_fraction, errstat, replace_fill=r8_missval)
    call tiof_close (cld, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_cloud_params: reading file "//trim(cloud_file), &
                       errstat)
      return
    endif

    ! Force the cloud parameters into physical bounds, avoiding missing data
    where (cloud_fraction > r8_missval .and. cloud_fraction < 0.0_r8)
      cloud_fraction = 0.0_r8
    elsewhere (cloud_fraction > 1.0_r8)
      cloud_fraction = 1.0_r8
    endwhere
    where (cloud_top_pressure > r8_missval .and. cloud_top_pressure < 0.0_r8)
      cloud_top_pressure = 0.0_r8
    endwhere

  end subroutine read_cloud_params

end module output_tools
