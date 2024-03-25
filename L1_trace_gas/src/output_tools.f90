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
  use ctrlvars, only: yn_diagnostic_run, yn_refseccor, yn_scat_weights, &
       yn_stratrop, yn_gems
  use sao_pge_utils, only: calc_relaz_angle

  implicit none
  private

  public create_output_file, close_output_file, write_wavcal_output, &
    write_radfit_output, write_fitting_statistics, write_common_mode, &
    write_albedo, write_gas_profile, write_scattering_weights, &
    write_amf_correction, write_temperature_profile, write_refspec_database, &
    write_reference_sector_corrected_column, &
    write_solar_wavecal_diagnostics, &
    write_radiance_wavecal_diagnostics, copy_pixel_corners, &
    copy_metadata, copy_gpqf_attributes, label_output_file, &
    read_geofields, read_column_results, read_cloud_params

  type (tiof_file_type), private, save, target :: primary_output_file
  type (tiof_file_type), private, save, target :: diagnostic_output_file

  ! using fill values from the original code simplifies diffing output files
  real (kind=8), private, parameter :: &
    fill_short = -9999, &
    fill_float = -1.0e30, &
    fill_double = -1.0e30_r8

  integer, private, parameter :: molname_len = 33
  type, private :: molname_type
    character (len=molname_len) :: name, name_sans_spaces
    integer :: pge_idx
  end type
  type (molname_type) :: target_molecule

  ! These file variable names are different for NO2:
  character (len=tg_max_name_len) :: var_amf, var_amf_error
  character (len=tg_max_name_len) :: var_vertical_column, var_vertical_column_error

contains

  subroutine set_molecule_name (pge_idx, errstat)
    implicit none
    integer, intent(in) :: pge_idx
    integer, intent(inout) :: errstat

    ! index range matches pge_idx definitions from OMSAO_indices_module.f90
    integer, parameter :: beg_idx=10, end_idx=22
    character (len=molname_len), dimension(beg_idx:end_idx), &
      parameter :: mol_names = (/ &
      'chlorine dioxide                ' &
    , 'bromine monoxide                ' &
    , 'formaldehyde                    ' &
    , 'ozone                           ' &
    , 'nitrogen dioxide                ' &
    , 'sulfur dioxide                  ' &
    , 'glyoxal                         ' &
    , 'iodine monoxide                 ' &
    , 'water vapor                     ' &
    , 'nitrous acid                    ' &
    , 'collision induced oxygen complex' &
    , 'liquid water                    ' &
    , 'nitrogen dioxide                ' &
      /)
    character (len=128) :: msg
    character (len=molname_len) :: name_sans_spaces
    integer :: i, n

    if (errstat /= 0) return

    if (pge_idx < beg_idx .or. end_idx < pge_idx) then
      write (msg, *)'molecule_name: unsupported value pge_idx=',pge_idx
      call tell_error (tell_runtime_error, msg, errstat)
      return
    endif

    name_sans_spaces = mol_names (pge_idx)
    n = len_trim(name_sans_spaces)
    do i = 1, n
      if (name_sans_spaces(i:i) == ' ') then
        name_sans_spaces(i:i) = '_'
      endif
    enddo

    target_molecule % pge_idx = pge_idx
    target_molecule % name = mol_names (pge_idx)
    target_molecule % name_sans_spaces = name_sans_spaces

  end subroutine

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
                             dimids=[dimids(1)], &
                             long_name = "pixel index along slit")
    call tiof_varlist_append (varlist, errstat, tg_dim_step, nf90_int, &
                             dimids=[dimids(2)], &
                             long_name = "scan mirror position index")
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
                              long_name = "common mode spectrum", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_common_mode_wavelengths, &
                              nf90_float, &
                              dimids = dimids_commwvl_xtrack,  &
                              long_name = "common mode wavelength", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              fillvalue = fill_float)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_common_mode_count, &
                              nf90_int, &
                              dimids = [dimids_commwvl_xtrack(2)],  &
                              long_name = "common mode count", &
                              comment = "number of values averaged at each wavelength", &
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
    type (tiof_attlist_type) :: att_coord, att_amf_diag
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_levels_xtrack_step, dimsizes_levels_xtrack_step
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
                              [tg_dim_swt_level, tg_dim_xtrack, tg_dim_step], &
                              dimids_levels_xtrack_step, &
                              errstat, dimsizes = dimsizes_levels_xtrack_step)

    ! coordinates for 2D variables
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(tg_var_time) &
                              //' '//trim(tg_var_longitude) &
                              //' '//trim(tg_var_latitude))
    call tiof_attlist_append (att_amf_diag, errstat, "coordinates", &
                              att_text = trim(tg_var_time) &
                              //' '//trim(tg_var_longitude) &
                              //' '//trim(tg_var_latitude))
    call tiof_attlist_append (att_amf_diag, errstat, "flag_meanings", &
                              att_text = "geometric_AMF glint snow_correction "// &
                              "no_cloud_pressure adjusted_surface_pressure "// &
                              "adjusted_cloud_pressure no_albedo no_cloud_fraction "// &
                              "no_gas_profile no_scattering_weights AMF_disabled")
    call tiof_attlist_append (att_amf_diag, errstat, "flag_masks", &
                              att_i4 = [1, 2, 4, 8, 16, 32, 2048, 4096, 8192, 16384, 32768])
    ! append amf variables
    chunksizes(1) = dimsizes_levels_xtrack_step(1)            ! level dimension
    chunksizes(2) = min(dimsizes_levels_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = min(dimsizes_levels_xtrack_step(3), 128)  ! step dimension

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_scattering_weights, &
                              nf90_float, &
                              dimids = dimids_levels_xtrack_step,  &
                              long_name = "scattering weights", &
                              comment = "vertical profile of scattering weights", &
                              valid_min = 0.0_r8, &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes, &
                              attlist = att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_gas_profile, &
                              nf90_float, &
                              dimids = dimids_levels_xtrack_step,  &
                              long_name = "vertical profile of "//trim(target_molecule % name)//" partial column", &
                              units = "molecules/cm^2", &
                              valid_min = 0.0_r8, &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes, &
                              attlist = att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_albedo, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "surface albedo", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              attlist = att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_temperature_profile, &
                              nf90_float, &
                              dimids = dimids_levels_xtrack_step,  &
                              long_name = "air temperature", &
                              units = "K", &
                              valid_min = 0.0_r8, &
                              valid_max = 400.0_r8, &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes, &
                              attlist = att_coord)
                        
    call tiof_varlist_append (varlist, errstat, &
                              var_amf, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = trim(target_molecule % name)//" air mass factor", &
                              comment = "total "//trim(target_molecule % name)//" air mass factor (AMF) "// &
                              "calculated from surface to top of atmosphere", &
                              valid_min = 0.0_r8, &
                              fillvalue = fill_float, &
                              attlist = att_coord)
    if (.false.) then
    call tiof_varlist_append (varlist, errstat, &
                              var_amf_error, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = trim(target_molecule % name)//" air mass factor uncertainty", &
                              valid_min = 0.0_r8, &
                              fillvalue = fill_float, &
                              attlist = att_coord)
    endif
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_diagnostic_flag, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              long_name = trim(target_molecule % name)//" air mass factor diagnostic flag ", &
                              fillvalue = fill_short, &
                              attlist = att_amf_diag)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_eff_cld_frac, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "effective cloud fraction", &
                              comment = "effective cloud fraction from cloud retrieval", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = -1.0_r8, &
                              attlist = att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_cloud_fraction, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "cloud fraction", &
                              comment = "cloud radiance fraction for AMF computation", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = -1.0_r8, &
                              attlist = att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_amf_cloud_pressure, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "cloud pressure", &
                              comment = "cloud pressure for AMF computation", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1200.0_r8], &
                              fillvalue = fill_float, &
                              attlist = att_coord)
    IF (yn_stratrop) THEN
       call tiof_varlist_append (varlist, errstat, &
                                 tg_var_amf_troposphere, &
                                 nf90_float, &
                                 dimids = dimids_xtrack_step,  &
                                 long_name = trim(target_molecule % name)//" tropospheric air mass factor", &
                                 valid_min = 0.0_r8, &
                                 fillvalue = fill_float, &
                                 attlist = att_coord)
       call tiof_varlist_append (varlist, errstat, &
                                 tg_var_amf_stratosphere, &
                                 nf90_float, &
                                 dimids = dimids_xtrack_step,  &
                                 long_name = trim(target_molecule % name)//" stratospheric air mass factor", &
                                 valid_min = 0.0_r8, &
                                 fillvalue = fill_float, &
                                 attlist = att_coord)
    END IF

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)

    call tiof_attlist_free (att_coord)
    call tiof_attlist_free (att_amf_diag)

  end subroutine append_amf_vars

  subroutine append_diagnostic_vars (obj, dimlist, errstat)
    implicit none

    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_convergence_flag
    integer, dimension(1) :: dimid_xtrack
    integer, dimension(2) :: dimids_xtrack_step, dimids_refwavl_xtrack
    integer, dimension(3) :: dimids_var_xtrack_step, dimsizes_var_xtrack_step
    integer, dimension(3) :: dimids_commwvl_xtrack_step, dimsizes_commwvl_xtrack_step
    integer, dimension(3) :: dimids_refwavl_xtrack_refspec, dimsizes_refwavl_xtrack_refspec
    integer, dimension(3) :: chunksizes
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat /= 0) return

    ! lookup dimids for relevant array shapes
    call tiof_dimlist_lookup (dimlist, [tg_dim_xtrack], dimid_xtrack, errstat)
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

    call tiof_attlist_append (att_convergence_flag, errstat, "coordinates", &
                              att_text = trim(tg_var_time) &
                              //' '//trim(tg_var_longitude) &
                              //' '//trim(tg_var_latitude))
    call tiof_attlist_append (att_convergence_flag, errstat, "flag_meanings", &
                              att_text = "failed maxiter_exceeded suspect good")
    call tiof_attlist_append (att_convergence_flag, errstat, "flag_values", &
                              att_i4 = [-2,-1,0,1])

    ! append diagnostic variables
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_wavelengths, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              units = "nm", &
                              long_name = "calibrated solar spectrum wavelengths", &
                              valid_range = [100.0_r8, 1000.0_r8], &
                              fillvalue = fill_double)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_residuals, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              long_name = "solar spectrum residuals", &
                              comment = "fit residuals from solar spectrum wavelength calibration", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_convergence_flag, &
                              nf90_short, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar wavelength calibration convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8], &
                              fillvalue = fill_short, &
                              attlist=att_convergence_flag)

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radcal_wavelengths, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              units = "nm", &
                              long_name = "calibrated radiance wavelengths", &
                              comment = "fit residuals from radiance wavelength calibration", &
                              valid_range = [100.0_r8, 1000.0_r8], &
                              fillvalue = fill_double)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radcal_residuals, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              long_name = "radiance spectrum residuals", &
                              comment = "fit residuals from radiance spectrum wavelength calibration", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radcal_convergence_flag, &
                              nf90_short, &
                              dimids = dimid_xtrack,  &
                              long_name = "radiance wavelength calibration convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8], &
                              fillvalue = fill_short, &
                              attlist = att_convergence_flag)

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_iteration_count, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "radiance fit iteration count", &
                              valid_range = [0.0_r8, 32767.0_r8], &
                              fillvalue = fill_short)

    chunksizes(1) = dimsizes_var_xtrack_step(1)           ! var dimension
    chunksizes(2) = min(dimsizes_var_xtrack_step(2),128)  ! xtrack dimension
    chunksizes(3) = 1                                     ! step dimension

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_params, &
                              nf90_double, &
                              dimids = dimids_var_xtrack_step, &
                              long_name = "radiance fit parameter", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_errors, &
                              nf90_double, &
                              dimids = dimids_var_xtrack_step, &
                              long_name = "radiance fit parameter uncertainty", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_correl, &
                              nf90_double, &
                              dimids = dimids_var_xtrack_step, &
                              long_name = "radiance fit parameter correlation", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_commwvl_xtrack_step(1)            ! wavelength dimension
    chunksizes(2) = min(dimsizes_commwvl_xtrack_step(2),128)   ! xtrack dimension
    chunksizes(3) = 1                                          ! step dimension

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_measured_spectrum, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              long_name = "measured radiance spectrum", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_measured_wavelengths, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              long_name = "measured radiance wavelength", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_model_spectrum, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              long_name = "model radiance spectrum", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_weights, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              long_name = "radiance spectrum fit weights", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_radfit_residuals, &
                              nf90_double, &
                              dimids = dimids_commwvl_xtrack_step, &
                              long_name = "radiance spectrum fit residuals", &
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
    chunksizes(2) = min(dimsizes_refwavl_xtrack_refspec(2), 128)    ! xtrack dimension
    chunksizes(3) = 1                                               ! refspec dimension

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_refspec, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack_refspec, &
                              long_name = "reference spectrum", &
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
                              comment = "reference spectrum wavelengths", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              tg_var_refspec_norm, &
                              nf90_double, &
                              dimids = [dimids_refwavl_xtrack_refspec(3) ], &
                              comment = "reference spectrum normalization factors", &
                              valid_range=[-1e30_r8, 1e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_convergence_flag)

    call append_common_mode_vars (obj, dimlist, errstat)

  end subroutine append_diagnostic_vars

  subroutine append_column_vars (obj, dimlist, amf_wvl, errstat)
    use OMSAO_indices_module, only : pge_no2_idx, pge_hcho_idx, pge_o2o2_idx
    implicit none

    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    real (kind=r8), intent(in) :: amf_wvl
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist_geo, varlist_qa
    type (tiof_varlist_type), target :: varlist, varlist_supp
    type (tiof_varlist_type), target :: varlist_tmp
    type (tiof_attlist_type) :: att_coord, att_latbnd, att_lonbnd
    type (tiof_attlist_type) :: att_main_dqf, att_convergence_flag
    type (tiof_attlist_type) :: att_time
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_corner_xtrack_step

    character (len=32) :: slant_column_units
    character (len=32) :: epoch_buf

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

    ! optional attribute lists:
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(tg_var_time) &
                              //' '//trim(tg_var_longitude) &
                              //' '//trim(tg_var_latitude))

    call tiof_attlist_append (att_main_dqf, errstat, "coordinates", &
                              att_text = trim(tg_var_time) &
                              //' '//trim(tg_var_longitude) &
                              //' '//trim(tg_var_latitude))
    call tiof_attlist_append (att_main_dqf, errstat, "flag_meanings", &
                              att_text = "normal suspicious bad")
    call tiof_attlist_append (att_main_dqf, errstat, "flag_values", &
                              att_i4 = [0,1,2])

    call tiof_attlist_append (att_convergence_flag, errstat, "coordinates", &
                              att_text = trim(tg_var_time) &
                              //' '//trim(tg_var_longitude) &
                              //' '//trim(tg_var_latitude))
    call tiof_attlist_append (att_convergence_flag, errstat, "flag_meanings", &
                              att_text = "failed maxiter_exceeded suspect good")
    call tiof_attlist_append (att_convergence_flag, errstat, "flag_values", &
                              att_i4 = [-2,-1,0,1])

    if (target_molecule % pge_idx == pge_no2_idx) then
      ! For NO2, separate contributions from stratosphere/troposphere will be derived
      ! in post-processing. For this reason, selected NO2 file variable names have
      ! the word "total". These columns go to support data group.
      varlist_tmp = varlist_supp
      var_amf       = trim(tg_var_amf)//"_total"
      !var_amf_error = trim(tg_var_amf)//"_total_uncertainty"
      var_vertical_column       = trim(tg_var_vertical_column)//"_total"
      var_vertical_column_error = trim(tg_var_vertical_column)//"_total_uncertainty"
    else
      ! For all other molecules, the vertical column goes to the "product" group.
      varlist_tmp = varlist
      var_amf       =      tg_var_amf
      !var_amf_error = trim(tg_var_amf)//"_uncertainty"
      var_vertical_column       =      tg_var_vertical_column
      var_vertical_column_error = trim(tg_var_vertical_column)//"_uncertainty"
    endif

    ! data field variables with optional attribute lists:
 
    if (amf_wvl > 0.0) then
      call tiof_varlist_append (varlist_tmp, errstat, &
                                var_vertical_column, &
                                nf90_double, &
                                dimids = dimids_xtrack_step,  &
                                long_name = trim(target_molecule % name)//" vertical column", &
                                units = "molecules/cm^2", &
                                comment = trim(target_molecule % name)// &
                                " vertical column determined from fitted slant column"// &
                                " and total AMF calculated from surface to top of atmosphere", &
                                fillvalue = fill_double, &
                                attlist=att_coord)

       call tiof_varlist_append (varlist_tmp, errstat, &
                                 var_vertical_column_error, &
                                 nf90_double, &
                                 dimids = dimids_xtrack_step,  &
                                 long_name = trim(target_molecule % name)//" vertical column uncertainty", &
                                 units = "molecules/cm^2", &
                                 fillvalue = fill_double, &
                                 attlist=att_coord)
    endif
    if (target_molecule % pge_idx == pge_no2_idx) then
        varlist_supp = varlist_tmp
    else
        varlist = varlist_tmp
    endif

    call tiof_varlist_append (varlist, errstat, &
                              tg_var_main_dqf, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              long_name = "main data quality flag", &
                              valid_range = [0.0_r8, 2.0_r8], &
                              fillvalue = fill_short, &
                              attlist=att_main_dqf)
    call tiof_push_group (obj, tg_grp_product, errstat)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist)

    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_radfit_rms_residual, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "radiance fit RMS residual", &
                              valid_range = [0.0_r8, 0.01_r8], &
                              fillvalue = fill_double, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_qa, errstat, &
                              tg_var_radfit_convergence_flag, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "radiance fit convergence flag", &
                              valid_range = [-10.0_r8, 12344.0_r8], &
                              fillvalue = fill_short, &
                              attlist=att_convergence_flag)

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_def_vars (obj, varlist_qa, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist_qa)
    call tiof_attlist_free (att_convergence_flag)

    epoch_buf(:)=''
    call tiof_mktimestamp_str (0.0_r8, epoch_buf, errstat)

    call tiof_attlist_append (att_time, errstat, "calendar", &
                              att_text = "gregorian")
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_time, &
                              nf90_double, &
                              dimids = [dimids_xtrack_step(2)],  &
                              standard_name = "time", &
                              long_name = "radiance exposure start time", &
                              units = "seconds since "//trim(epoch_buf), &
                              fillvalue = fill_double, &
                              attlist=att_time)

    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = tg_var_latitude_bounds)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_latitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              standard_name = "latitude", &
                              long_name = "pixel center latitude", &
                              comment = "latitude at pixel center", &
                              units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_latbnd)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_latitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              long_name = "pixel corner latitude", &
                              comment = "latitude at pixel corners (SW,SE,NE,NW)", &
                              !units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              fillvalue = fill_float)

    call tiof_attlist_append (att_lonbnd, errstat, "bounds", &
                              att_text = tg_var_longitude_bounds)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_longitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              standard_name = "longitude", &
                              long_name = "pixel center longitude", &
                              comment = "longitude at pixel center", &
                              units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_lonbnd)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_longitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              long_name = "pixel corner longitude", &
                              comment = "longitude at pixel corners (SW,SE,NE,NW)", &
                              !units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float)

    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_sz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "solar zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_sa_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "solar azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_vz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "viewing zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_va_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "viewing azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_geo, errstat, &
                              tg_var_relative_azimuth, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "relative azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_push_group (obj, tg_grp_geolocation, errstat)
    call tiof_def_vars (obj, varlist_geo, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist_geo)

    if (target_molecule % pge_idx == pge_o2o2_idx) then
      slant_column_units = "molecules^2/cm^5"
    else
      slant_column_units = "molecules/cm^2"
    endif

    call tiof_varlist_append (varlist_supp, errstat, &
                              tg_var_fitted_slant_column, &
                              nf90_double, &
                              dimids = dimids_xtrack_step,  &
                              long_name = trim(target_molecule % name)//" fitted slant column", &
                              units = trim(slant_column_units), &
                              fillvalue = fill_double, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_supp, errstat, &
                              tg_var_fitted_slant_column_error, &
                              nf90_double, &
                              dimids = dimids_xtrack_step,  &
                              long_name = trim(target_molecule % name)//" fitted slant column uncertainty", &
                              units = trim(slant_column_units), &
                              fillvalue = fill_double, &
                              attlist=att_coord)
    if (.not.yn_gems) then
      call tiof_varlist_append (varlist_supp, errstat, &
                                tg_var_snowice_fraction, &
                                nf90_float, &
                                dimids = dimids_xtrack_step,  &
                                long_name = "fraction of pixel area covered by snow and/or ice", &
                                valid_range = [0.0_r8, 1.0_r8], &
                                fillvalue = fill_float, &
                                attlist=att_coord)
    endif
    call tiof_varlist_append (varlist_supp, errstat, &
                              tg_var_terrain_height, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              long_name = "terrain height", &
                              units = "m", &
                              valid_range = [-1000.0_r8, 10000.0_r8], &
                              fillvalue = fill_short, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_supp, errstat, &
                              tg_var_gpqf, &
                              nf90_int, &
                              dimids = dimids_xtrack_step, &
                              long_name = "ground pixel quality flag", &
                              fillvalue = fill_short, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist_supp, errstat, &
                              tg_var_surface_pressure, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "surface pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1200.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    IF (yn_stratrop) THEN
      call tiof_varlist_append (varlist_supp, errstat, &
                                tg_var_tropopause_pressure, &
                                nf90_float, &
                                dimids = dimids_xtrack_step,  &
                                long_name = "tropopause pressure", &
                                units = "hPa", &
                                valid_range = [0.0_r8, 1200.0_r8], &
                                fillvalue = fill_float, &
                                attlist=att_coord)
    END IF
    if (yn_refseccor .and. (target_molecule % pge_idx == pge_hcho_idx)) then
      call tiof_varlist_append (varlist_supp, errstat, &
                                tg_var_refsec_corr, &
                                nf90_double, &
                                dimids = dimids_xtrack_step, &
                                long_name = "reference sector correction", &
                                comment = "reference sector correction based on differences"//&
                                " between model and retrieval over reference sector", &
                                units = "molecules/cm^2", &
                                fillvalue = fill_double, &
                                attlist=att_coord)
    endif
    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_def_vars (obj, varlist_supp, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist_supp)

    call tiof_attlist_free (att_coord)
    call tiof_attlist_free (att_latbnd)
    call tiof_attlist_free (att_lonbnd)
    call tiof_attlist_free (att_main_dqf)
    call tiof_attlist_free (att_time)

  end subroutine append_column_vars

  subroutine create_diagnostic_file (product_filename, num_steps, num_xtrack, n_comm_wvl, &
                                     nwavel_max, max_rs_idx, n_fitvar_rad, &
                                     errstat)
    implicit none
    character (len=*), intent(in) :: product_filename
    integer (kind=i4), intent(in) :: num_steps, num_xtrack, n_comm_wvl, nwavel_max
    integer (kind=i4), intent(in) :: max_rs_idx, n_fitvar_rad
    integer, intent(inout) :: errstat

    integer :: dot_pos
    character (len=*), parameter :: diag_label = '_diag'
    character (len=len(product_filename)+len(diag_label)) :: filename
    type (tiof_file_type), pointer :: obj
    type (tiof_dimlist_type) :: dimlist

    if (errstat /= 0) return

    obj => diagnostic_output_file

    ! Generate the filename
    dot_pos = scan (trim(product_filename), '.', back=.true.)
    filename = product_filename (1:dot_pos-1)//diag_label//'.nc'

    ! Create the file
    call tiof_create (obj, filename, nf90_clobber, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_diagnostic_file: creating file "//trim(filename), &
                       errstat)
      return
    endif

    ! Define a dimension list.
    call tiof_dimlist_append (dimlist, tg_dim_step, num_steps, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_xtrack, num_xtrack, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_commwvl, n_comm_wvl, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_fitvar, n_fitvar_rad, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_refwavl, nwavel_max, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_refspec, max_rs_idx, errstat)
    call tiof_dimlist_append (dimlist, tg_dim_pair, 2, errstat)
    call tiof_def_dims (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_diagnostic_file: defining dimensions in "//trim(filename), &
                       errstat)
      return
    endif

    call write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_diagnostic_file: writing coordinate variables to "//trim(filename), &
                       errstat)
      return
    endif

    call append_diagnostic_vars (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_diagnostic_file: defining diagnostic variables in "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_dimlist_free (dimlist)

  end subroutine

  !> Create netCDF format Level 2 product file
  !! @param[in] filename   netCDF output file name
  !! @param[in] pge_idx    Index of target molecule [integer]
  !! @param[in] amf_wvl    AMF wavelength (>0 when AMF vars are present)
  !! @param[in] num_steps  Number of scan steps
  !! @param[in] num_xtrack  Number of cross-track pixels
  !! @param[in] num_swlevels  Number of height levels used in AMF climatologies
  !! @param[in] n_comm_wvl  Number of common-mode wavelengths
  !! @param[in] nwavel_max  Maximum number of reference spectrum wavelengths
  !! @param[in] max_rs_idx  Maximum reference spectrum index
  !! @param[in] n_fitvar_rad  Number of radiance spectrum fit variables
  !! @param[inout]  errstat  Error status variable
  subroutine create_output_file (filename, pge_idx, amf_wvl, &
                                 num_steps, num_xtrack, num_swlevels, &
                                 n_comm_wvl, nwavel_max, max_rs_idx, &
                                 n_fitvar_rad, errstat)
    implicit none
    character (len=*), intent(in) :: filename
    integer (kind=i4), intent(in) :: pge_idx, num_steps, num_xtrack, num_swlevels, &
      n_comm_wvl, nwavel_max, max_rs_idx, n_fitvar_rad
    real (kind=r8), intent(in) :: amf_wvl
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_dimlist_type) :: dimlist

    if (errstat /= 0) return

    obj => primary_output_file

    call set_molecule_name (pge_idx, errstat)
    if (errstat /= 0) return

    ! Create a file.
    call tiof_create (obj, filename, nf90_clobber, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file: creating file "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_history_append_cmdline (obj)
    call tiof_put_git_commit_hash (obj, errstat)

    ! Create default groups.
    call tiof_def_group (obj, tg_grp_product, errstat)
    call tiof_def_group (obj, tg_grp_geolocation, errstat)
    call tiof_def_group (obj, tg_grp_support_data, errstat)
    call tiof_def_group (obj, tg_grp_qa_stats, errstat)
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
    if (yn_scat_weights) then
      call tiof_dimlist_append (dimlist, tg_dim_swt_level, num_swlevels, errstat)
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

    call append_column_vars (obj, dimlist, amf_wvl, errstat)
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

    call tiof_dimlist_free (dimlist)

    if (yn_diagnostic_run) then
      call create_diagnostic_file (filename, num_steps, num_xtrack, n_comm_wvl, &
                                   nwavel_max, max_rs_idx, n_fitvar_rad, errstat)
      if (errstat /= 0) return
    endif

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

    type (tiof_file_type), pointer :: obj, obj_diag
    real (kind=r8), dimension(1:n_rad_wvl, 1:nxtrack, 0:nblock-1) :: residuals
    real (kind=r8), dimension(:,:,:), pointer :: waves, meas, model, weights
    real (kind=r4), dimension(1:nxtrack,0:nblock-1) :: relative_azimuth
    integer :: i,j

    if (errstat /= 0) return

    obj => primary_output_file
    obj_diag => null()

    ! result_vars
    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put2d_r8 (obj, tg_var_fitted_slant_column, [iline,0], [nblock, nxtrack], &
                        result_vars % column_amount (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_r8 (obj, tg_var_fitted_slant_column_error, [iline,0], [nblock, nxtrack], &
                        result_vars % column_uncert (1:nxtrack, 0:nblock-1), errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_put2d_r8 (obj, tg_var_radfit_rms_residual, [iline,0], [nblock, nxtrack], &
                        result_vars % fit_rms_residual (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_i2 (obj, tg_var_radfit_convergence_flag, [iline,0], [nblock, nxtrack], &
                        result_vars % fit_convergence_flag (1:nxtrack, 0:nblock-1), errstat)
    call tiof_pop_group (obj, errstat)

    if (yn_diagnostic_run) then
      obj_diag => diagnostic_output_file
      call tiof_put2d_i2 (obj_diag, tg_var_radfit_iteration_count, [iline,0], [nblock, nxtrack], &
                          result_vars % fit_iteration_count (1:nxtrack, 0:nblock-1), errstat)

      call tiof_put3d_r8 (obj_diag, tg_var_radfit_params, [iline,0,0], [nblock,nxtrack,n_fitvar_rad], &
                          radfit_diagnostics % params(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj_diag, tg_var_radfit_errors, [iline,0,0], [nblock,nxtrack,n_fitvar_rad], &
                          radfit_diagnostics % errors(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)
      call tiof_put3d_r8 (obj_diag, tg_var_radfit_correl, [iline,0,0], [nblock,nxtrack,n_fitvar_rad], &
                          radfit_diagnostics % correl(1:n_fitvar_rad,1:nxtrack,0:nblock-1), &
                          errstat)

      model   => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 1, 0:nblock-1)
      meas    => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 2, 0:nblock-1)
      waves   => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 3, 0:nblock-1)
      weights => radfit_diagnostics % fitspc(1:n_rad_wvl, 1:nxtrack, 4, 0:nblock-1)

      call tiof_put3d_r8 (obj_diag, tg_var_radfit_model_spectrum, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], model, errstat)
      call tiof_put3d_r8 (obj_diag, tg_var_radfit_measured_spectrum, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], meas, errstat)
      call tiof_put3d_r8 (obj_diag, tg_var_radfit_measured_wavelengths, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl],  waves, errstat)
      call tiof_put3d_r8 (obj_diag, tg_var_radfit_weights, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], weights, errstat)

      residuals(:,:,:) = meas - model
      call tiof_put3d_r8 (obj_diag, tg_var_radfit_residuals, &
                          [iline,0,0], [nblock,nxtrack,n_rad_wvl], residuals, errstat)
    endif

    ! Compute relative azimuth angle
    relative_azimuth(:,:) = fill_float
    do j=0,nblock-1
      do i=1,nxtrack
        relative_azimuth(i,j) = calc_relaz_angle( &
             input_vars%solar_azimuth(i,j), input_vars%viewing_azimuth(i,j))
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
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    if (.not.yn_gems) then
      call tiof_put2d_r4 (obj, tempo_var_snowice_fraction, [iline,0], [nblock,nxtrack], &
                          input_vars % snow_ice_fraction (1:nxtrack, 0:nblock-1), errstat)
    endif
    call tiof_put2d_i2 (obj, tg_var_terrain_height, [iline,0], [nblock,nxtrack], &
                        input_vars % terrain_height (1:nxtrack, 0:nblock-1), errstat)
    call tiof_put2d_i4 (obj, tg_var_gpqf, [iline,0], [nblock,nxtrack], &
                        input_vars % ground_pixel_quality_flag (1:nxtrack, 0:nblock-1), errstat)
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

    type (tiof_file_type), pointer :: obj, obj_diag

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_wavcal_output: failed", errstat)
      return
    endif

    if (yn_diagnostic_run) then
      obj_diag => diagnostic_output_file
      call tiof_put1d_i2 (obj_diag, tg_var_solcal_convergence_flag, [0], [nxtrack], &
                          result_vars % solcal_convergence_flag (1:nxtrack), errstat)
      call tiof_put1d_i2 (obj_diag, tg_var_radcal_convergence_flag, [0], [nxtrack], &
                          result_vars % radcal_convergence_flag (1:nxtrack), errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_write_error, &
                         "write_wavcal_output: diagnostic output failed", errstat)
        return
      endif
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

    type (tiof_file_type), pointer :: obj, obj_diag
    type (tiof_attlist_type) :: attlist
    integer :: ncp, nsl

    obj => primary_output_file
    obj_diag => null()

    call tiof_attlist_append (attlist, errstat, "num_crosstrack_pixels", att_i4 = [stats % num_crosstrack_pixels])
    call tiof_attlist_append (attlist, errstat, "num_scan_lines", att_i4 = [stats % num_scan_lines])
    call tiof_attlist_append (attlist, errstat, "num_good_input", att_i4 = [stats % num_good_input])
    call tiof_attlist_append (attlist, errstat, "num_good_output", att_i4 = [stats % num_good_output])
    call tiof_attlist_append (attlist, errstat, "num_suspect_output", att_i4 = [stats % num_suspect_output])
    call tiof_attlist_append (attlist, errstat, "num_bad_output", att_i4 = [stats % num_bad_output])
    call tiof_attlist_append (attlist, errstat, "num_converged", att_i4 = [stats % num_converged])
    call tiof_attlist_append (attlist, errstat, "num_failed_convergence", att_i4 = [stats % num_failed_convergence])
    call tiof_attlist_append (attlist, errstat, "num_exceeded_iterations", att_i4 = [stats % num_exceeded_iterations])
    call tiof_attlist_append (attlist, errstat, "num_out_of_bounds", att_i4 = [stats % num_out_of_bounds])
    call tiof_attlist_append (attlist, errstat, "percent_good_output", att_r4 = [stats % percent_good_output])
    call tiof_attlist_append (attlist, errstat, "percent_bad_output", att_r4 = [stats % percent_bad_output])
    call tiof_attlist_append (attlist, errstat, "percent_suspect_output", att_r4 = [stats % percent_suspect_output])

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_def_atts (obj, attlist, nf90_global, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_attlist_free (attlist)

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
      obj_diag => diagnostic_output_file
      call tiof_put1d_string (obj_diag, tg_var_radfit_param_names, 0, num_params, &
                              param_names(1:num_params), errstat)
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
  !! @param[in] nxtrack  Number of cross-track pixels
  !! @param[in] ntimes  Number of scans
  !! @param[in] nlevels  Number of altitudes in vertical profile climatology
  !! @param[in] apriori_source  String describing the source of the gas profile
  !! @param[inout] errstat  Error status variable
  subroutine write_gas_profile (gas_profile, &
                                nxtrack, ntimes, nlevels, &
                                apriori_source, errstat)
    implicit none

    integer, intent(in) :: nxtrack, ntimes, nlevels
    real (kind=r8), dimension (1:nlevels, 1:nxtrack, 0:ntimes-1), intent(in) :: gas_profile
    character (len=*), intent(in) :: apriori_source
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_attlist_type) :: attlist

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_attlist_append (attlist, errstat, "apriori_source", att_text = apriori_source)
    call tiof_def_atts (obj, attlist, nf90_global, errstat)
    call tiof_attlist_free (attlist)

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put3d_r4 (obj, tg_var_amf_gas_profile, [0,0,0], &
         [ntimes,nxtrack,nlevels], &
         real(gas_profile(1:nlevels,1:nxtrack, 0:ntimes-1), kind=4), errstat)
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
    real (kind=r8), dimension (1:nlevels, 1:nxtrack, 0:ntimes-1), intent(in) :: scattw
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put3d_r4 (obj, tg_var_amf_scattering_weights, [0,0,0], &
         [ntimes,nxtrack,nlevels], &
         real(scattw (1:nlevels,1:nxtrack, 0:ntimes-1), kind=4), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_scattering_weights", errstat)
      return
    endif
  end subroutine write_scattering_weights

  !> Write AMF vertical temperature profile climatology to Level 2 product file
  !! @param[in] temperature_profile  Temperature vertical profile climatology for each pixel [K]
  !! @param[in] nxtrack  Number of cross-track pixels
  !! @param[in] ntimes  Number of scans
  !! @param[in] nlevels  Number of altitudes in vertical profile climatology
  !! @param[inout] errstat  Error status variable
  subroutine write_temperature_profile (temperature_profile, &
                                        nxtrack, ntimes, nlevels, &
                                        errstat)
    implicit none
    integer, intent(in) :: nxtrack, ntimes, nlevels
    real (kind=r8), dimension (1:nlevels, 1:nxtrack, 0:ntimes-1), intent(in) :: temperature_profile
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj
    if (errstat /= 0) return
    obj => primary_output_file
    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put3d_r4 (obj, tg_var_amf_temperature_profile, [0,0,0], &
         [ntimes,nxtrack,nlevels], &
         real(temperature_profile(1:nlevels,1:nxtrack, 0:ntimes-1), kind=4), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "in write_temperature_profile", errstat)
      return
    endif
  end subroutine write_temperature_profile

  !> Write AMF correction to Level 2 product file
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[in] ntimes Number of scans
  !! @param[in] amf_corr Air mass factor (AMF) correction
  !! @param[in] amf_corr_column  AMF-corrected column density
  !! @param[in] amf_corr_column_uncertainty  Uncertainty in AMF-corrected column density
  !! @param[in] yn_write_cloud_variables   If \a .true., write cloud variables
  !!                                      to Level 2 product file
  !! @param[in] crfrc cloud radiance fraction 2D array
  !! @param[inout] errstat  Error status variable
  subroutine write_amf_correction (nxtrack, ntimes, amf_corr, &
                                   amf_corr_column, amf_corr_column_uncertainty, &
                                   yn_write_cloud_variables, crfrc, errstat)
    use OMSAO_omidata_module, only : amf_correction_type
    use OMSAO_indices_module, only : pge_no2_idx
    implicit none

    integer, intent(in) :: nxtrack, ntimes
    type (amf_correction_type), intent(in) :: amf_corr
    real (kind=r8), dimension(1:nxtrack,0:ntimes-1), intent(in) :: amf_corr_column
    real (kind=r8), dimension(1:nxtrack,0:ntimes-1), intent(in) :: amf_corr_column_uncertainty, crfrc
    logical, intent(in) :: yn_write_cloud_variables
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_attlist_type) :: attlist
    integer :: status, varid_surface_pressure

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put2d_i2 (obj, tg_var_amf_diagnostic_flag, [0,0], [ntimes,nxtrack], &
                        amf_corr % diagnostic_flag (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, var_amf, [0,0], [ntimes,nxtrack], &
                        amf_corr % amf_molecule_specific (1:nxtrack, 0:ntimes-1), errstat)

    if (yn_stratrop) then
       call tiof_put2d_r8 (obj, tg_var_amf_stratosphere, [0,0], [ntimes,nxtrack], &
                           amf_corr % amf_molecule_stratospheric (1:nxtrack, 0:ntimes-1), errstat)
       call tiof_put2d_r8 (obj, tg_var_amf_troposphere, [0,0], [ntimes,nxtrack], &
                           amf_corr % amf_molecule_tropospheric (1:nxtrack, 0:ntimes-1), errstat)
    end if

    if (yn_write_cloud_variables) then
      call tiof_put2d_r8 (obj, tg_var_amf_cloud_fraction, [0,0], &
           [ntimes,nxtrack], crfrc (1:nxtrack, 0:ntimes-1), errstat)
      call tiof_put2d_r8 (obj, tg_var_eff_cld_frac, [0,0], [ntimes,nxtrack], &
           amf_corr % cloud_fraction (1:nxtrack, 0:ntimes-1), errstat)
      !call tiof_put2d_r8 (obj, tg_var_amf_cloud_pressure, [0,0], [ntimes,nxtrack], &
      !                    amf_corr % cloud_pressure (1:nxtrack, 0:ntimes-1), errstat)
      ! FIXME: netcdf error will occur if cloud_pressure is not representable as a float
      call tiof_put2d_r4 (obj, tg_var_amf_cloud_pressure, [0,0], [ntimes,nxtrack], &
                          real(amf_corr % cloud_pressure (1:nxtrack, 0:ntimes-1),kind=r4), errstat)
    endif
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_put2d_r4 (obj, tg_var_surface_pressure, [0,0], [ntimes,nxtrack], &
                        amf_corr % surface_pressure (1:nxtrack, 0:ntimes-1), errstat)
    ! Pressure profile is parameterized by p(z) = eta_a(z) + eta_b(z) * psurf
    status = nf90_inq_varid (obj % groupid, tg_var_surface_pressure, varid_surface_pressure)
    call tiof_attlist_append (attlist, errstat, "Eta_A", att_r4 = amf_corr % eta_a)
    call tiof_attlist_append (attlist, errstat, "Eta_B", att_r4 = amf_corr % eta_b)
    call tiof_def_atts (obj, attlist, varid_surface_pressure, errstat)
    call tiof_attlist_free (attlist)

    if (yn_stratrop) then
       call tiof_put2d_r4 (obj, tg_var_tropopause_pressure, [0,0], [ntimes,nxtrack], &
                           amf_corr % tropopause_pressure (1:nxtrack, 0:ntimes-1), errstat)
    end if

    call tiof_pop_group (obj, errstat)

    if (target_molecule % pge_idx == pge_no2_idx) then
        call tiof_push_group (obj, tg_grp_support_data, errstat)
    else
        call tiof_push_group (obj, tg_grp_product, errstat)
    endif
    call tiof_put2d_r8 (obj, var_vertical_column, [0,0], [ntimes,nxtrack], &
                        amf_corr_column (1:nxtrack, 0:ntimes-1), errstat)
    call tiof_put2d_r8 (obj, var_vertical_column_error, [0,0], [ntimes,nxtrack], &
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

    obj => diagnostic_output_file

    call tiof_put2d_r8 (obj, tg_var_common_mode_spectrum, [0,0], [nxtrack,n_comm_wvl], &
                        common_mode % refspecdata (1:n_comm_wvl,1:nxtrack), errstat)
    call tiof_put2d_r8 (obj, tg_var_common_mode_wavelengths, [0,0], [nxtrack,n_comm_wvl], &
                        common_mode % refspecwavs (1:n_comm_wvl,1:nxtrack), errstat)
    call tiof_put2d_i2 (obj, tg_var_common_mode_ccd_pixel_range, [0,0], [2,nxtrack], &
                        common_mode % ccdpixel (1:nxtrack,1:2), errstat)
    call tiof_put1d_i4 (obj, tg_var_common_mode_count, [0], [nxtrack], &
                        common_mode % refspeccount (1:nxtrack), errstat)

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

    obj => diagnostic_output_file

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
    call tiof_put2d_r8 (obj, tg_var_refsec_corr, [0,0], [ntimes, nxtrack], &
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

    obj => diagnostic_output_file

    call tiof_put2d_r8 (obj, tg_var_solcal_wavelengths, [0,0], [nxtrack, nwaves], &
                        waves (1:nwaves,1:nxtrack), errstat)
    call tiof_put2d_r8 (obj, tg_var_solcal_residuals, [0,0], [nxtrack, nwaves], &
                        resid (1:nwaves,1:nxtrack), errstat)

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

    obj => diagnostic_output_file

    call tiof_put2d_r8 (obj, tg_var_radcal_wavelengths, [0,0], [nxtrack, nwaves], &
                        waves (1:nwaves,1:nxtrack), errstat)
    call tiof_put2d_r8 (obj, tg_var_radcal_residuals, [0,0], [nxtrack, nwaves], &
                        resid (1:nwaves,1:nxtrack), errstat)

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
                        tmp(1:4,1:nxtrack,1:ntimes), errstat, &
                        replace_fill=real(fill_float, kind=4))
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
                        tmp(1:4,1:nxtrack,1:ntimes), errstat, &
                        replace_fill=real(fill_float, kind=4))
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

    corners_copied = .true.

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

    call tiof_write_epoch_timestamp (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "copy_metadata: copying from "//trim(l1bfile), &
                       errstat)
    endif
  end subroutine copy_metadata

  subroutine copy_gpqf_attributes (l1bfile, rad_group, errstat)
    implicit none
    character (len=*), intent(in) :: l1bfile
    character (len=*), intent(in) :: rad_group
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj
    type (tiof_file_type) :: l1b

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_open (l1bfile, l1b, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "copy_gpqf_attributes: opening file "//trim(l1bfile), &
                       errstat)
      return
    endif

    call tiof_push_group (l1b, rad_group, errstat)

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_copy_attr (l1b, tg_var_gpqf, obj, tg_var_gpqf, &
                         (/"long_name    ", &
                          "flag_meanings", &
                          "flag_values  "/), errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_close (l1b, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
                       "copy_gpqf_attributes: copying attributes", errstat)
      return
    endif

  end subroutine copy_gpqf_attributes

  subroutine label_output_file (label, processing_version, errstat)
    implicit none
    character (len=*), intent(in) :: label
    integer, intent(in) :: processing_version
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file
    call tiof_label_product (obj, label, 2, processing_version, errstat)
  end subroutine label_output_file

  !> Read geolocation fields
  !! @param[in] l1bfile level 1 radiance filename
  !! @param[in] ntimes  Number of scans
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[inout] lat   Latitude [deg]
  !! @param[inout] lon   Longitude [deg]
  !! @param[inout] sza   Solar zenith angle [deg]
  !! @param[inout] vza   Viewing zenith angle [deg]
  !! @param[inout] saa   Solar azimuth angle [deg]
  !! @param[inout] vaa   Viewing azimuth angle [deg]
  !! @param[inout] thgt  Terrain height
  !! @param[inout] time  Time [s]
  !! @param[inout] errstat  Error status variable
  subroutine read_geofields (l1bfile, ntimes, nxtrack, lat, lon, sza, vza, saa, vaa, thgt, time, errstat)
    use OMSAO_precision_module, only : i2, i4, r4
    use OMSAO_omidata_module, only: omi_radiance_swathname
    use OMSAO_parameters_module, only : max_latitude, max_longitude, r4_missval, r8_missval
    use tio_module
    use netcdf, only: nf90_nowrite
    implicit none

    character (len=*), intent(in) :: l1bfile
    integer (kind=i4), intent(in) :: ntimes, nxtrack
    ! these arrays are (1:nxtrack,0:ntimes-1)
    real (kind=r4), dimension(:,:), intent(inout) :: lat, lon, sza, vza, &
         saa, vaa, thgt
    real (kind=r8), dimension(:), intent(inout) :: time
    integer, intent(inout) :: errstat

    type (tiof_file_type) :: obj
    integer (kind=i2), dimension(nxtrack,ntimes) :: i2_thgt

    if (errstat /= 0) return

    call tiof_open (l1bfile, obj, nf90_nowrite, errstat)
    call tiof_get1d_r8 (obj, tg_var_time, [0], [ntimes], time(1:ntimes), errstat, &
                        replace_fill=r8_missval)
    call tiof_push_group (obj, omi_radiance_swathname, errstat)
    call tiof_get2d_r4 (obj, tg_var_latitude, [0,0], [ntimes, nxtrack], lat(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_longitude, [0,0], [ntimes, nxtrack], lon(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_sz_angle, [0,0], [ntimes, nxtrack], sza(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_vz_angle, [0,0], [ntimes, nxtrack], vza(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_sa_angle, [0,0], [ntimes, nxtrack], saa(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r4 (obj, tg_var_va_angle, [0,0], [ntimes, nxtrack], vaa(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_i2 (obj, tg_var_terrain_height, [0,0], [ntimes, nxtrack], i2_thgt(1:nxtrack,1:ntimes), errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_close (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "in read_geofields", errstat)
      return
    endif

    where (abs(lat(1:nxtrack,1:ntimes)) > max_latitude)
      lat(1:nxtrack,1:ntimes) = r4_missval
    endwhere
    where (abs(lon(1:nxtrack,1:ntimes)) > max_longitude)
      lon(1:nxtrack,1:ntimes) = r4_missval
    endwhere

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

    call tiof_push_group (obj, tg_grp_support_data, errstat)
    call tiof_get2d_r8 (obj, tg_var_fitted_slant_column, [0,0], [ntimes, nxtrack], col(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_r8 (obj, tg_var_fitted_slant_column_error, [0,0], [ntimes, nxtrack], col_unc(1:nxtrack,1:ntimes), errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, tg_grp_qa_stats, errstat)
    call tiof_get2d_r8 (obj, tg_var_radfit_rms_residual, [0,0], [ntimes, nxtrack], rms(1:nxtrack,1:ntimes), errstat)
    call tiof_get2d_i2 (obj, tg_var_radfit_convergence_flag, [0,0], [ntimes, nxtrack], &
                        convergence_flag(1:nxtrack,1:ntimes), errstat)
    call tiof_pop_group (obj, errstat)

    if (yn_scat_weights) then
      call tiof_push_group (obj, tg_grp_support_data, errstat)
      call tiof_get2d_r8 (obj, var_amf, [0,0], [ntimes, nxtrack], amf(1:nxtrack,1:ntimes), errstat)
      call tiof_pop_group (obj, errstat)
    endif

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
    character (len=256) :: logmsg

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
    call tiof_get2d_r8 (cld, "cloud_pressure", [0,0], [ntimes,nxtrack], &
                        cloud_top_pressure, errstat, replace_fill=r8_missval)
    call tiof_get2d_r8 (cld, "cloud_fraction", [0,0], [ntimes,nxtrack], &
                        cloud_fraction, errstat, replace_fill=r8_missval)
    call tiof_close (cld, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_cloud_params: reading file "//trim(cloud_file), &
                       errstat)
      return
    endif

    ! Clouds are already cropped to physical bounds in cloud code, but fill
    ! value is different for cloud fraction between codes, so need to change
    ! here.
    where (cloud_fraction < 0.0_r8)
      cloud_fraction = -1.0_r8
    endwhere
    where (cloud_top_pressure > r8_missval .and. cloud_top_pressure < 0.0_r8)
      cloud_top_pressure = 0.0_r8
    endwhere

  end subroutine read_cloud_params

end module output_tools
