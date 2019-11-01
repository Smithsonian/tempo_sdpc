!> Level 2 data product output functions
!! @file
module l2_tio_class
  use netcdf
  use tell_module
  use tio_module
  use o3t_names_module
  implicit none
  private

  public l2_tio_create, l2_tio_close, &
    l2_tio_write_etc, l2_tio_write_geo, l2_tio_write_fields, &
    l2_tio_write_skipped_fields, l2_tio_write_mqf, &
    l2_tio_copy_metadata, l2_tio_label_output_file

  type (tiof_file_type), private, save, target :: primary_output_file

  ! using fill values from the original code simplifies diffing output files
  real (kind=8), private :: fill_float, fill_double
  real (kind=8), private, parameter :: fill_short = -32767

contains

  subroutine append_product_vars (obj, dimlist, have_omi_data, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    logical, intent(in) :: have_omi_data
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord
    integer, dimension(2) :: dimids_xtrack_step

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3t_dim_xtrack, o3t_dim_step], &
                              dimids_xtrack_step, &
                              errstat)

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    ! data field variables with optional attribute lists:
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3t_var_longitude) &
                              //' '//trim(o3t_var_latitude))
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_column_amount_o3, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "best total ozone solution", &
                              units = "DU", &
                              valid_range = [50.0_8, 700.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    if (have_omi_data) then
      call tiof_varlist_append (varlist, errstat, &
                              o3t_var_mqf, &
                              nf90_ubyte, &
                              dimids = [dimids_xtrack_step(2)], &
                              comment = "measurement quality flag", &
                              valid_range = [0.0_8, 254.0_8])
    endif

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_radiative_cloudfrac, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "radiative cloud fraction = fc*Ic331/Im331", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_cloudfrac_param, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "mixed LER model (cloud fraction) parameter", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_quality_flag, &
                              nf90_ushort, &
                              dimids = dimids_xtrack_step, &
                              comment = "quality flags", &
                              valid_range = [0.0_8, 65534.0_8])

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_o3_below_cloud, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "ozone below fractional cloud", &
                              units = "DU", &
                              valid_range = [0.0_8, 100.0_8], &
                              fillvalue = fill_float)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_so2_index, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "SO2 index", &
                              valid_range = [-300.0_8, 300.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_uv_aerosol_index, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "UV aerosol index", &
                              valid_range = [-30.0_8, 30.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_def_vars (obj, varlist, errstat)

    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_coord)

  end subroutine

  subroutine append_support_vars (obj, dimlist, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_layer_xtrack_step, dimids_wavel_xtrack_step
    integer, dimension(3) :: dimsizes_layer_xtrack_step, dimsizes_wavel_xtrack_step
    integer, dimension(3) :: chunksizes
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3t_dim_xtrack, o3t_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3t_dim_layer, o3t_dim_xtrack, o3t_dim_step], &
                              dimids_layer_xtrack_step, &
                              errstat, dimsizes = dimsizes_layer_xtrack_step)
    call tiof_dimlist_lookup (dimlist, &
                              [o3t_dim_wavelength, o3t_dim_xtrack, o3t_dim_step], &
                              dimids_wavel_xtrack_step, &
                              errstat, dimsizes = dimsizes_wavel_xtrack_step)

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    ! data field variables with optional attribute lists:
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3t_var_longitude) &
                              //' '//trim(o3t_var_latitude))

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_lut_wavelength, &
                              nf90_float, &
                              dimids = [dimids_wavel_xtrack_step(1)], &
                              comment = "lookup table wavelength", &
                              units = "nm", &
                              valid_range = [300.0_8, 400.0_8], &
                              fillvalue = fill_float)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_cloud_pressure, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "effective cloud pressure", &
                              units = "hPA", &
                              valid_range = [0.0_8, 1013.25_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_terrain_pressure, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "terrain pressure", &
                              units = "hPA", &
                              valid_range = [0.0_8, 1013.25_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_algorithm_flags, &
                              nf90_ubyte, &
                              dimids = dimids_xtrack_step, &
                              comment = "algorithm flags", &
                              valid_range = [0.0_8, 13.0_8])

    chunksizes(1) = dimsizes_layer_xtrack_step(1)            ! layer dimension
    chunksizes(2) = min(dimsizes_layer_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_apriori_layer_o3, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "a priori ozone profile", &
                              valid_range = [0.0_8, 125.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_radbpix_flag_accepted, &
                              nf90_ushort, &
                              dimids = dimids_xtrack_step, &
                              comment = "radiance bad pixel flag accepted", &
                              valid_range = [0.0_8, 65534.0_8])

    chunksizes(1) = dimsizes_layer_xtrack_step(1)            ! layer dimension
    chunksizes(2) = min(dimsizes_layer_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_layer_efficiency, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "algorithmic layer efficiency", &
                              valid_range = [0.0_8, 10.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_wavel_xtrack_step(1)            ! wavel dimension
    chunksizes(2) = min(dimsizes_wavel_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_dndr, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step, &
                              comment = "reflectivity sensitivity ratio", &
                              valid_range = [-200.0_8, 0.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_wavel_xtrack_step(1)            ! wavel dimension
    chunksizes(2) = min(dimsizes_wavel_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_nvalue, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step, &
                              comment = "measured N-value", &
                              valid_range = [0.0_8, 600.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_reflectivity_331, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "effective surface reflectivity at 331 nm", &
                              units = "percent", &
                              valid_range = [-15.0_8, 115.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_reflectivity_360, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "effective surface reflectivity at 360 nm", &
                              units = "percent", &
                              valid_range = [-15.0_8, 115.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    chunksizes(1) = dimsizes_wavel_xtrack_step(1)            ! wavel dimension
    chunksizes(2) = min(dimsizes_wavel_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_residual, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step, &
                              comment = "N-value residual", &
                              valid_range = [-32.0_8, 32.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_wavel_xtrack_step(1)            ! wavel dimension
    chunksizes(2) = min(dimsizes_wavel_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_residual_step1, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step, &
                              comment = "step 1 N-value residual", &
                              valid_range = [-32.0_8, 32.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_wavel_xtrack_step(1)            ! wavel dimension
    chunksizes(2) = min(dimsizes_wavel_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_residual_step2, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step, &
                              comment = "step 2 N-value residual", &
                              valid_range = [-32.0_8, 32.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_wavel_xtrack_step(1)            ! wavel dimension
    chunksizes(2) = min(dimsizes_wavel_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_sensitivity, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step, &
                              comment = "ozone sensitivity ratio, dN/dOmega", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    chunksizes(1) = dimsizes_wavel_xtrack_step(1)            ! wavel dimension
    chunksizes(2) = min(dimsizes_wavel_xtrack_step(2), 128)  ! xtrack dimension
    chunksizes(3) = 1                                        ! step dimension
    ! FIXME - choose more optimal chunk sizes
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_temp_sensitivity_ratio, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step, &
                              comment = "ozone weighted temperature sensitivity ratio, dN/dT", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_step1_o3, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "step 1 ozone solution", &
                              units = "DU", &
                              valid_range = [50.0_8, 700.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_step2_o3, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "step 2 ozone solution", &
                              units = "DU", &
                              valid_range = [50.0_8, 700.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_cal_adjustment, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step(1:2), &
                              comment = "calibration adjustment", &
                              valid_range = [-10.0_8, 10.0_8], &
                              fillvalue = fill_float)

    call tiof_def_vars (obj, varlist, errstat)

    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_coord)
  end subroutine

  subroutine append_geolocation_vars (obj, dimlist, have_omi_data, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    logical, intent(in) :: have_omi_data
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_latbnd, att_lonbnd, att_coord
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_corner_xtrack_step

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3t_dim_xtrack, o3t_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3t_dim_corner, o3t_dim_xtrack, o3t_dim_step], &
                              dimids_corner_xtrack_step, &
                              errstat)

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    ! data field variables with optional attribute lists:
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_time, &
                              nf90_double, &
                              dimids = [dimids_xtrack_step(2)],  &
                              comment = "exposure start time", &
                              units = "s", &
                              valid_range = [-5.0e9_8, 1.e10_8], &
                              fillvalue = fill_double)

    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = o3t_var_latitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_latitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "latitude at pixel center", &
                              units = "degrees_north", &
                              valid_range = [-90.0_8, 90.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_latbnd)

    call tiof_attlist_append (att_lonbnd, errstat, "bounds", &
                              att_text = o3t_var_longitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_longitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "longitude at pixel center", &
                              units = "degrees_east", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_lonbnd)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_latitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              comment = "latitude at pixel corners (sw,se,ne,nw)", &
                              units = "degrees_north", &
                              valid_range = [-90.0_8, 90.0_8], &
                              fillvalue = fill_float)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_longitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              comment = "longitude at pixel corners (sw,se,ne,nw)", &
                              units = "degrees_east", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float)

    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3t_var_longitude) &
                              //' '//trim(o3t_var_latitude))
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_sz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_sa_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_vz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_8, 70.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_va_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_relaz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "relative azimuth angle (sun + 180 - view)", &
                              units = "degrees", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_terrain_height, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "terrain height", &
                              units = "m", &
                              valid_range = [-200.0_8, 10000.0_8], &
                              fillvalue = fill_short, &
                              attlist=att_coord)

    call tiof_varlist_append (varlist, errstat, &
                              o3t_var_geoflg, &
                              nf90_ushort, &
                              dimids = dimids_xtrack_step, &
                              comment = "ground pixel quality flag", &
                              valid_range = [0.0_8, 65534.0_8])

    if (have_omi_data) then
      call tiof_varlist_append (varlist, errstat, &
                              o3t_var_xtrack_qf, &
                              nf90_ubyte, &
                              dimids = dimids_xtrack_step, &
                              comment = "cross-track quality flag", &
                              valid_range = [0.0_8, 254.0_8])
    endif

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_latbnd)
    call tiof_attlist_free (att_lonbnd)
    call tiof_attlist_free (att_coord)

  end subroutine

  subroutine write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, num_layers, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(in) :: num_steps, num_xtrack, num_layers
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(num_xtrack) :: xtrack_indices
    integer, dimension(num_steps) :: step_indices
    integer, dimension(num_layers) :: layer_indices
    integer :: i, dimids(3)

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, [o3t_dim_xtrack, o3t_dim_step, o3t_dim_layer], dimids, errstat)

    ! netcdf coordinate variables:
    call tiof_varlist_append (varlist, errstat, o3t_dim_xtrack, nf90_int, &
                             dimids=[dimids(1)])
    call tiof_varlist_append (varlist, errstat, o3t_dim_step, nf90_int, &
                             dimids=[dimids(2)])
    call tiof_varlist_append (varlist, errstat, o3t_dim_layer, nf90_int, &
                             dimids=[dimids(3)])
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)

    ! mirror step indices may be overwritten later by copying from the input
    ! radiance file
    step_indices = [(i, i=0,num_steps-1)]
    call tiof_put1d_i4 (obj, o3t_dim_step, [0], [num_steps], step_indices, errstat)

    xtrack_indices = [(i, i=0,num_xtrack-1)]
    call tiof_put1d_i4 (obj, o3t_dim_xtrack, [0], [num_xtrack], xtrack_indices, errstat)

    layer_indices = [(i, i=0,num_layers-1)]
    call tiof_put1d_i4 (obj, o3t_dim_layer, [0], [num_layers], layer_indices, errstat)

  end subroutine write_coordinate_vars

  !> Create netCDF format Level 2 product file
  !! @param[in] filename   output netCDF file name
  !! @param[in] num_steps  Number of scan steps
  !! @param[in] num_xtrack  Number of cross-strack pixels
  !! @param[in] num_layers  Number of vertical profile layers in forward model
  !! @param[in] num_wavel  Number of wavelengths
  !! @param[in] have_omi_data true = OMI, false = TEMPO
  !! @param[inout] errstat   Error status variable
  subroutine l2_tio_create (filename, num_steps, num_xtrack, &
                            num_layers, num_wavel, have_omi_data, errstat)
    implicit none
    character (len=*), intent(in) :: filename
    integer (kind=4), intent(in) :: num_steps, num_xtrack, num_layers, num_wavel
    logical, intent(in) :: have_omi_data
    integer, intent(inout) :: errstat

    integer (kind=4), external :: r8fill
    type (tiof_file_type), pointer :: obj
    type (tiof_dimlist_type) :: dimlist

    if (errstat < 0) return

    obj => primary_output_file

    if ((r8fill(fill_float) /= 0) .or. (r8fill(fill_double) /= 0)) then
      call tell_error (tell_runtime_error, &
                       "l2_tio_create: defining fill values", &
                       errstat)
      return
    endif

    ! Create a file.
    call tiof_create (obj, filename, nf90_clobber, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: creating file "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_write_epoch_timestamp (obj, errstat)

    ! Create default groups.
    call tiof_def_group (obj, o3t_grp_product, errstat)
    call tiof_def_group (obj, o3t_grp_geolocation, errstat)
    call tiof_def_group (obj, o3t_grp_support_data, errstat)
    call tiof_def_group (obj, o3t_grp_qa_stats, errstat)
    call tiof_def_group (obj, o3t_grp_metadata, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create:  defining groups in "//trim(filename), &
                       errstat)
      return
    endif

    ! Define a dimension list.
    call tiof_dimlist_append (dimlist, o3t_dim_step, num_steps, errstat)
    call tiof_dimlist_append (dimlist, o3t_dim_xtrack, num_xtrack, errstat)
    call tiof_dimlist_append (dimlist, o3t_dim_corner, 4, errstat)
    call tiof_dimlist_append (dimlist, o3t_dim_layer, num_layers, errstat)
    call tiof_dimlist_append (dimlist, o3t_dim_wavelength, num_wavel, errstat)
    call tiof_def_dims (obj, dimlist, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: defining dimensions in "//trim(filename), &
                       errstat)
      return
    endif

    call write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, num_layers, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: writing coordinate variables to "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_push_group (obj, o3t_grp_product, errstat)
    call append_product_vars (obj, dimlist, have_omi_data, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: creating product group in "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_push_group (obj, o3t_grp_geolocation, errstat)
    call append_geolocation_vars (obj, dimlist, have_omi_data, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: creating geolocation group in "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_push_group (obj, o3t_grp_support_data, errstat)
    call append_support_vars (obj, dimlist, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: creating support_data group in "//trim(filename), &
                       errstat)
      return
    endif

    !call tiof_push_group (obj, o3t_grp_qa_stats, errstat)
    !call append_qa_vars (obj, dimlist, errstat)
    !call tiof_pop_group (obj, errstat)
    !
    !call tiof_push_group (obj, o3t_grp_metadata, errstat)
    !call append_metadata_vars (obj, dimlist, errstat)
    !call tiof_pop_group (obj, errstat)

    call tiof_dimlist_free (dimlist)

    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: creating file "//trim(filename), &
                       errstat)
      return
    endif

  end subroutine

  !> Close Level 2 product file
  !! @param[inout] errstat Error status variable
  subroutine l2_tio_close (errstat)
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    obj => primary_output_file

    call tiof_close (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_error, "l2_tio_close failed", errstat)
    endif

  end subroutine

  subroutine l2_tio_copy_metadata (l1b, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: l1b
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_copy_granule_ident (l1b, obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "copying L1 metadata to output file", &
                       errstat)
    endif

  end subroutine l2_tio_copy_metadata

  subroutine l2_tio_label_output_file (label, processing_version, errstat)
    implicit none
    character (len=*), intent(in) :: label
    integer, intent(in) :: processing_version
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file
    call tiof_label_product (obj, label, processing_version, errstat)
  end subroutine l2_tio_label_output_file

  !> Write calibration adjustment variables to Level 2 product file
  !! @param[in] nwavel  Number of wavelengths
  !! @param[in] wl_com  Wavelengths
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[in] swpcr  N-value correction for given year and day
  !! @param[inout] errstat  Error status variable
  subroutine l2_tio_write_etc (nwavel, wl_com, nxtrack, swpcr, errstat)
    implicit none
    integer, intent(in) :: nwavel, nxtrack
    real (kind=4), dimension(:), intent(in) :: wl_com
    real (kind=4), dimension(:,:), allocatable, intent(in) :: swpcr
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    ! support_data group
    call tiof_push_group (obj, o3t_grp_support_data, errstat)
    call tiof_put1d_r4 (obj, o3t_var_lut_wavelength, [0], [nwavel], &
                        wl_com(1:nwavel), errstat)
    if (allocated (swpcr)) then
      call tiof_put2d_r4 (obj, o3t_var_cal_adjustment, [0,0], [nxtrack, nwavel], &
                          swpcr(1:nwavel,1:nxtrack), errstat)
    endif
    call tiof_pop_group (obj, errstat)

  end subroutine

  !> Write scan line geolocation variables to Level 2 product file
  !! @param[in] iline  Scan line index
  !! @param[in] nxtrack Number of cross-track pixels
  !! @param[in] have_omi_data logical, true = OMI, false = TEMPO
  !! @param[inout] errstat  Error status variable
  !! @details
  !! The geolocation variables are passed via the \a O3T_radgeo_class module.
  subroutine l2_tio_write_geo (iline, step_index, nxtrack, have_omi_data, &
       errstat)
    use O3T_radgeo_class
    implicit none
    integer, intent(in) :: iline, step_index, nxtrack
    logical, intent(in) :: have_omi_data
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    ! main group
    call tiof_push_group (obj, "/", errstat)
    call tiof_put1d_i4 (obj, o3t_dim_step, [iline], [1], [step_index], errstat)
    call tiof_pop_group (obj, errstat)

    ! geolocation group
    call tiof_push_group (obj, o3t_grp_geolocation, errstat)
    call tiof_put1d_r8 (obj, o3t_var_time, [iline], [1], [time], errstat)
    call tiof_put1d_ui2 (obj, o3t_var_geoflg, [iline,0], [1, nxtrack], &
                        geoflg(1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_latitude, [iline,0], [1, nxtrack], &
                        latitude(1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_longitude, [iline,0], [1, nxtrack], &
                        longitude(1:nxtrack), errstat)
    call tiof_put2d_r4 (obj, o3t_var_latitude_bounds, [iline,0,0], [1, nxtrack,4], &
                        lat_bounds(1:4,1:nxtrack), errstat)
    call tiof_put2d_r4 (obj, o3t_var_longitude_bounds, [iline,0,0], [1, nxtrack,4], &
                        lon_bounds(1:4,1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_sz_angle, [iline,0], [1, nxtrack], &
                        szenith(1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_sa_angle, [iline,0], [1, nxtrack], &
                        sazimuth(1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_vz_angle, [iline,0], [1, nxtrack], &
                        vzenith(1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_va_angle, [iline,0], [1, nxtrack], &
                        vazimuth(1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_relaz_angle, [iline,0], [1, nxtrack], &
                        phiArray(1:nxtrack), errstat)
    call tiof_put1d_i2 (obj, o3t_var_terrain_height, [iline,0], [1, nxtrack], &
                        height(1:nxtrack), errstat)
    if (have_omi_data) then
      call tiof_put1d_ui1 (obj, o3t_var_xtrack_qf, [iline,0], [1, nxtrack], &
                         anomflg(1:nxtrack), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_geo: writing to geolocation group", &
                       errstat)
      return
    endif

    ! support data group
    call tiof_push_group (obj, o3t_grp_support_data, errstat)
    call tiof_put1d_r4 (obj, o3t_var_cloud_pressure, [iline,0], [1, nxtrack], &
                        pcArray(1:nxtrack), errstat)
    call tiof_put1d_r4 (obj, o3t_var_terrain_pressure, [iline,0], [1, nxtrack], &
                        ptArray(1:nxtrack), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_geo: writing to support_data group", &
                       errstat)
      return
    endif

  end subroutine

  !> Write intermediate variables to Level 2 product file
  !! @param[in]  iline   scan line index
  !! @param[in]  ix      cross-track pixel index
  !! @param[in]  nwavel  Number of wavelengths
  !! @param[in]  nlayers  Number of vertical profile layers
  !! @param[in]  algflg   Algorithm convergence flag
  !! @param[in]  qaflags  Quality assurance flags
  !! @param[in]  radbadpixflgs  Radiance bad pixel flags
  !! @param[in]  stp1oz   Step 1 ozone column
  !! @param[in]  stp2oz   Step 2 ozone column
  !! @param[in]  stp3oz   Step 3 ozone column
  !! @param[in]  oz_cld   Ozone column below fractional cloud
  !! @param[in]  aerind   Aerosol index
  !! @param[in]  so2ind   Sulfur dioxide index
  !! @param[in]  pixsurf  Pixel surface land cover metadata
  !! @param[in]  eff      Layer efficiency factor
  !! @param[in]  aprfoz   A-priori ozone profile
  !! @param[inout] errstat  Error status flag
  subroutine l2_tio_write_fields (iline, ix, nwavel, nlayers, &
                                  algflg, qaflags, radbadpixflgs, &
                                  stp1oz, stp2oz, stp3oz, oz_cld, aerind, so2ind, &
                                  pixsurf, eff, aprfoz, errstat)
    use O3T_pixel_class
    use O3T_L2output_class, only : dndr, xnvalm, res_stp1, res_stp2, res_stp3, &
      dndomega_t, dNdT
    implicit none
    integer, intent(in) :: iline, ix, nwavel, nlayers
    integer (kind=1), intent(in) :: algflg
    integer (kind=2), intent(in) :: qaflags, radbadpixflgs
    real (kind=4), intent(in) :: stp1oz, stp2oz, stp3oz, oz_cld, aerind, so2ind
    type (O3T_pixcover_type), intent(inout) :: pixsurf
    real (kind=4), dimension(:), intent(in) :: eff, aprfoz
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    ! support data group
    call tiof_push_group (obj, o3t_grp_support_data, errstat)
    call tiof_put1d_ui1 (obj, o3t_var_algorithm_flags, [iline,ix-1], [1, 1], &
                         [algflg], errstat)
    call tiof_put1d_ui2 (obj, o3t_var_radbpix_flag_accepted, [iline,ix-1], [1, 1], &
                         [radbadpixflgs], errstat)
    call tiof_put1d_r4 (obj, o3t_var_step1_o3, [iline,ix-1], [1, 1], &
                        [stp1oz], errstat)
    call tiof_put1d_r4 (obj, o3t_var_step2_o3, [iline,ix-1], [1, 1], &
                        [stp2oz], errstat)
    call tiof_put1d_r4 (obj, o3t_var_reflectivity_331, [iline,ix-1], [1, 1], &
                        [100.0*pixsurf%ref], errstat)
    call tiof_put1d_r4 (obj, o3t_var_reflectivity_360, [iline,ix-1], [1, 1], &
                        [100.0*pixsurf%ref360], errstat)

    call tiof_put1d_r4 (obj, o3t_var_apriori_layer_o3, [iline,ix-1,0], [1,1,nlayers], &
                        aprfoz(1:nlayers), errstat)
    call tiof_put1d_r4 (obj, o3t_var_layer_efficiency, [iline,ix-1,0], [1,1,nlayers], &
                        eff(1:nlayers), errstat)

    call tiof_put1d_r4 (obj, o3t_var_nvalue, [iline,ix-1,0], [1,1,nwavel], &
                        100.0*xnvalm(1:nwavel), errstat)
    call tiof_put1d_r4 (obj, o3t_var_dndr, [iline,ix-1,0], [1,1,nwavel], &
                        dndr(1:nwavel), errstat)
    call tiof_put1d_r4 (obj, o3t_var_residual_step1, [iline,ix-1,0], [1,1,nwavel], &
                        res_stp1(1:nwavel), errstat)
    call tiof_put1d_r4 (obj, o3t_var_residual_step2, [iline,ix-1,0], [1,1,nwavel], &
                        res_stp2(1:nwavel), errstat)
    call tiof_put1d_r4 (obj, o3t_var_residual, [iline,ix-1,0], [1,1,nwavel], &
                        res_stp3(1:nwavel), errstat)
    call tiof_put1d_r4 (obj, o3t_var_sensitivity, [iline,ix-1,0], [1,1,nwavel], &
                        dndomega_t(1:nwavel), errstat)
    call tiof_put1d_r4 (obj, o3t_var_temp_sensitivity_ratio, [iline,ix-1,0], [1,1,nwavel], &
                        dNdT(1:nwavel), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_fields: writing to support_data group", &
                       errstat)
      return
    endif

    ! product group
    call tiof_push_group (obj, o3t_grp_product, errstat)
    call tiof_put1d_r4 (obj, o3t_var_column_amount_o3, [iline,ix-1], [1, 1], &
                        [stp3oz], errstat)
    call tiof_put1d_ui2 (obj, o3t_var_quality_flag, [iline,ix-1], [1, 1], &
                         [qaflags], errstat)
    call tiof_put1d_r4 (obj, o3t_var_cloudfrac_param, [iline,ix-1], [1, 1], &
                        [pixsurf%clfrac], errstat)
    call tiof_put1d_r4 (obj, o3t_var_radiative_cloudfrac, [iline,ix-1], [1, 1], &
                        [pixsurf%rcf1], errstat)
    call tiof_put1d_r4 (obj, o3t_var_o3_below_cloud, [iline,ix-1], [1, 1], &
                        [oz_cld], errstat)
    call tiof_put1d_r4 (obj, o3t_var_so2_index, [iline,ix-1], [1, 1], &
                        [so2ind], errstat)
    call tiof_put1d_r4 (obj, o3t_var_uv_aerosol_index, [iline,ix-1], [1, 1], &
                        [aerind], errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_fields: writing to product group", &
                       errstat)
      return
    endif

  end subroutine

  !> Write parameters for skipped pixels to L2 product file
  !! @param[in]  iline   scan line index
  !! @param[in]  ix      cross-track pixel index
  !! @param[in]  nwavel  Number of wavelengths
  !! @param[in]  algflg   Algorithm convergence flag
  !! @param[in]  qaflags  Quality assurance flags
  !! @param[in]  radbadpixflgs  Radiance bad pixel flags
  !! @param[inout] errstat  Error status flag
  subroutine l2_tio_write_skipped_fields (iline, ix, nwavel, &
                                          algflg, qaflags, radbadpixflgs, errstat)
    use O3T_L2output_class, only : xnvalm
    implicit none
    integer, intent(in) :: iline, ix, nwavel
    integer (kind=1), intent(in) :: algflg
    integer (kind=2), intent(in) :: qaflags, radbadpixflgs
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    ! support data group
    call tiof_push_group (obj, o3t_grp_support_data, errstat)
    call tiof_put1d_ui1 (obj, o3t_var_algorithm_flags, [iline,ix-1], [1, 1], &
                         [algflg], errstat)
    call tiof_put1d_ui2 (obj, o3t_var_radbpix_flag_accepted, [iline,ix-1], [1, 1], &
                         [radbadpixflgs], errstat)
    call tiof_put1d_r4 (obj, o3t_var_nvalue, [iline,ix-1,0], [1,1,nwavel], &
                        100.0*xnvalm(1:nwavel), errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_skipped_fields: writing to support_data group", &
                       errstat)
      return
    endif

    ! product group
    call tiof_push_group (obj, o3t_grp_product, errstat)
    call tiof_put1d_ui2 (obj, o3t_var_quality_flag, [iline,ix-1], [1, 1], &
                         [qaflags], errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_skipped_fields: writing to product group", &
                       errstat)
      return
    endif

  end subroutine

  !> Write measurement quality flag to L2 product file
  !! @param[in]  iline   scan line index
  !! @param[in]  mqf     measurement quality flag.
  !! @param[inout]  errstat  Error status variable.
  subroutine l2_tio_write_mqf (iline, mqf, errstat)
    implicit none
    integer, intent(in) :: iline
    integer (kind=1), intent(in) :: mqf
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    ! product group
    call tiof_push_group (obj, o3t_grp_product, errstat)
    call tiof_put1d_ui1 (obj, o3t_var_mqf, [iline], [1], [mqf], errstat)
    call tiof_pop_group (obj, errstat)

  end subroutine

end module
