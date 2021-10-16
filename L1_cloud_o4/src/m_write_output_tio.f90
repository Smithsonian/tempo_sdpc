!> Write output to a TEMPO-format netCDF4 file
module m_write_output_tio
  use tio_module
  use tell_module
  use netcdf, only: nf90_clobber, nf90_double, nf90_float, nf90_short, &
       nf90_uint, nf90_int
  use m_read_input_tio, only: open_tio, close_tio

  private write_coordinate_vars, write_geo_struct, write_geo_data, &
       copy_pixel_corners
  private write_product_struct, write_product_data
  private write_support_struct, write_support_data

  public create_output_file_tio

  type (tiof_file_type), private, target :: primary_output_file

  !fill values
  real (kind=8), private, parameter :: fill_bit = -128, &
       fill_short = -32768, fill_int = fill_short, fill_float=-9999.0, &
       fill_double = fill_float

contains


  !>Top-level subroutine to create and populate an L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] outfile    name of L2 netCDF file to be written
  !> @param[in] l1_file    L1 radiance file from which to copy corners
  !> @param[in] swathname  radiance swath name
  !> @param[in] nstep      along-track dimension size
  !> @param[in] nxtrack    cross-track dimension size
  !> @param     errstat    error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   April 2021
  !-----------------------------------------------------------------------
  subroutine create_output_file_tio (outfile, l1_file, swathname, &
       nstep, nxtrack, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: outfile, l1_file, swathname
    integer (kind=4), intent(in) :: nstep, nxtrack
    !output variables
    integer (kind=4), intent(inout) :: errstat
    !local variables
    type (tiof_file_type), pointer :: tio_l2obj
    type (tiof_dimlist_type) :: dimlist

    if (errstat /= 0) return

    tio_l2obj => primary_output_file

    ! Create the file
    call tiof_create (tio_l2obj, outfile, nf90_clobber, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file_tio: creating file "//trim(outfile), &
           errstat)
      return
    endif

    call tiof_put_git_commit_hash (tio_l2obj, errstat)
    !call tiof_write_epoch_timestamp (tio_l2obj, errstat)

    ! Create default groups.
    call tiof_def_group (tio_l2obj, "product", errstat)
    call tiof_def_group (tio_l2obj, "geolocation", errstat)
    call tiof_def_group (tio_l2obj, "support_data", errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file_tio: failed to define groups", &
           errstat)
      return
    endif

    ! define the dimension list
    call tiof_dimlist_append (dimlist, "mirror_step", nstep, errstat)
    call tiof_dimlist_append (dimlist, "xtrack", nxtrack, errstat)
    call tiof_dimlist_append (dimlist, "corner", 4, errstat)
    call tiof_def_dims (tio_l2obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file_tio: failed to define dims", errstat)
      return
    endif

    ! coordinate variables
    call write_coordinate_vars (tio_l2obj, dimlist, nstep, nxtrack, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file_tio: failed to write cordinate vars", &
           errstat)
      return
    endif

    ! geolocation variable definitions
    call write_geo_struct (tio_l2obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write structures", &
           errstat)
      return
    endif

    ! geolocation data
    call write_geo_data (tio_l2obj, nstep, nxtrack, errstat)
    call copy_pixel_corners (l1_file, swathname, tio_l2obj, nstep, nxtrack, &
         errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write data", &
           errstat)
      return
    endif

    !hqw addition -------------------------------------------------
    ! product variable definitions
    call write_product_struct (tio_l2obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write structures", &
           errstat)
      return
    endif

    ! product data
    call write_product_data (tio_l2obj, nstep, nxtrack, errstat)
    call copy_pixel_corners (l1_file, swathname, tio_l2obj, nstep, nxtrack, &
         errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write data", &
           errstat)
      return
    endif

    ! support variable definitions
    call write_support_struct(tio_l2obj, dimlist, errstat)
    if (errstat /= 0) then 
       call tell_error (tell_io_write_error, &
            "creat_output_file: failed to write structures", &
            errstat)
       return
    endif

    ! support data
    call write_support_data(tio_l2obj, nstep, nxtrack, errstat)
    if (errstat /= 0) then
       call tell_error (tell_io_write_error, &
            "create_output_file: filed to write data", &
            errstat)
       return
    endif

    ! global metadata
    call write_tio_glbattr(tio_l2obj,errstat)

    !------------------------------------------------------------
    call tiof_dimlist_free (dimlist)

    call close_tio (tio_l2obj, errstat)

  end subroutine create_output_file_tio

  !>Write coordinate variables into L2 netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] tio_l2obj  file object to be written into
  !> @param[in] dimlist    list of dimension parameters
  !> @param[in] nstep      along-track dimension size
  !> @param[in] nxtrack    cross-track dimension size
  !> @param     errstat    error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan   April 2021
  !-----------------------------------------------------------------------
  subroutine write_coordinate_vars (tio_l2obj, dimlist, nstep, nxtrack, &
       errstat)

    implicit none

    !input variables
    type (tiof_file_type), intent(in) :: tio_l2obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer (kind=4), intent(in) :: nstep, nxtrack
    !output variables
    integer (kind=4), intent(inout) :: errstat
    !local variables
    type (tiof_varlist_type) :: varlist
    integer, dimension(nxtrack) :: xtrack_indices
    integer (kind=4), dimension(2) :: dimids_xtrack_step
    integer :: i


    if (errstat /= 0) return

    call tiof_dimlist_lookup (dimlist, &
                              ["xtrack     ", "mirror_step"], &
                              dimids_xtrack_step, &
                              errstat)

    ! netcdf coordinate variables
    call tiof_varlist_append (varlist, errstat, "xtrack", nf90_int, &
                             dimids=[dimids_xtrack_step(1)])
    call tiof_varlist_append (varlist, errstat, "mirror_step", nf90_int, &
                             dimids=[dimids_xtrack_step(2)])
    call tiof_def_vars (tio_l2obj, varlist, errstat)
    call tiof_varlist_free (varlist)

    !write xtrack indices, mirror_step will be copied from rad file
    xtrack_indices = [(i, i=0,nxtrack-1)]
    call tiof_put1d_i4 (tio_l2obj, "xtrack", [0], [nxtrack], &
         xtrack_indices, errstat)


  end subroutine write_coordinate_vars


  !> Create the structure for the geolocation data in L2 netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] tio_l2obj   file object to be written into
  !> @param[in] dimlist     list of dimension parameters
  !> @param     errstat     error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan   April 2021
  !-----------------------------------------------------------------------
  subroutine write_geo_struct(tio_l2obj, dimlist, errstat)

    implicit none

    !input variables
    type (tiof_file_type), intent(inout) :: tio_l2obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    !output variables
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_geo, att_latbnd, att_lonbnd
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_corner_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    !define r8 kind for use in setting parameter valid ranges
    integer, parameter :: r8 = kind(1.0d0)
    character (len=32) :: epoch_buf

    if (errstat /= 0) return

    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              ["xtrack     ", "mirror_step"], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              ["corner     ", "xtrack     ", "mirror_step"], &
                              dimids_corner_xtrack_step, &
                              errstat)


    epoch_buf(:)=''
    call tiof_mktimestamp_str (0.0_r8, epoch_buf, errstat)


    ! Geolocation Fields with optional attribute lists
    call tiof_attlist_append (att_geo, errstat, "coordinates", &
                              att_text = "longitude latitude")
    call tiof_varlist_append (varlist, errstat, &
                              "time", &
                              nf90_double, &
                              dimids = [dimids_xtrack_step(2)],  &
                              comment = "radiance exposure start time", &
                              units = "seconds since "//trim(epoch_buf), &
                              valid_range = [0.0_r8, 1.0e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              fillvalue = fill_double)
    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = "latitude_bounds")
    call tiof_varlist_append (varlist, errstat, &
                              "latitude", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "latitude at pixel center", &
                              units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_latbnd)
    call tiof_attlist_append (att_lonbnd, errstat, "bounds", &
                              att_text = "longitude_bounds")
    call tiof_varlist_append (varlist, errstat, &
                              "longitude", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "longitude at pixel center", &
                              units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_lonbnd)
    call tiof_varlist_append (varlist, errstat, &
                              "latitude_bounds", &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              long_name = "pixel corner latitude", &
                              comment = "latitude at pixel corners (SW,SE,NE,NW)", &
                              units = "degrees_north", &
                              valid_range = [-90.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              "longitude_bounds", &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              long_name = "pixel corner longitude", &
                              comment = "longitude at pixel corners (SW,SE,NE,NW)", &
                              units = "degrees_east", &
                              valid_range = [-180.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle)
    call tiof_varlist_append (varlist, errstat, &
                              "solar_zenith_angle", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_geo)
    call tiof_varlist_append (varlist, errstat, &
                              "solar_azimuth_angle", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_geo)
    call tiof_varlist_append (varlist, errstat, &
                              "viewing_zenith_angle", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 90.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_geo)
    call tiof_varlist_append (varlist, errstat, &
                              "viewing_azimuth_angle", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_geo)

    call tiof_varlist_append (varlist, errstat, &
                              "relative_azimuth_angle", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "relative azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_geo)

    call tiof_varlist_append (varlist, errstat, &
                              "terrain_height", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "terrain height", &
                              units = "m", &
                              valid_range = [-100.0_r8, 10000.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_geo)
! note that using a nf90_short causes the flags to be written out wrong
! should be nf90_ushort, but the fill value is unsupported...
!    call tiof_varlist_append (varlist, errstat, &
!                              "ground_pixel_quality_flag", &
!                              nf90_uint, &
!                              dimids = dimids_xtrack_step,  &
!                              comment = "ground pixel quality flag", &
!                              valid_range = [0.0_r8, 65535.0_r8], &
!                              fillvalue = fill_int, &
!                              attlist=att_geo)

    call tiof_push_group (tio_l2obj, "geolocation", errstat)
    call tiof_def_vars (tio_l2obj, varlist, errstat)
    call tiof_pop_group (tio_l2obj, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_geo)
    call tiof_attlist_free (att_latbnd)
    call tiof_attlist_free (att_lonbnd)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_geo_struct: failed", &
           errstat)
      return
    endif

  end subroutine write_geo_struct


  !> Write geolocation data into L2 netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in]  tio_l2obj  file object to be written into
  !> @param[in]  nsteps     along-track dimension size
  !> @param[in]  nxtrack    cross-track dimension size
  !> @param      errstat    error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   April 2021
  !-----------------------------------------------------------------------
  subroutine write_geo_data(tio_l2obj, nstep, nxtrack, errstat)
    use m_vars, only: rad_Latitude, rad_Longitude, rad_Time, &
         rad_SolarZenithAngle, rad_SolarAzimuthAngle, &
         rad_ViewingZenithAngle, rad_ViewingAzimuthAngle, &
         out_RelativeAzimuthAngle, out_TerrainHeight!, &
!        rad_GroundPixelQualityFlags

    implicit none

    !input variables
    integer, intent(in) :: nxtrack, nstep
    !output variables
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: tio_l2obj


    if (errstat /= 0) return

    tio_l2obj => primary_output_file

    call tiof_push_group (tio_l2obj, "geolocation", errstat)

    call tiof_put2d_r4 (tio_l2obj, "latitude", [0,0], &
         [nstep, nxtrack], rad_Latitude, errstat)

    call tiof_put2d_r4 (tio_l2obj, "longitude", [0,0], &
         [nstep, nxtrack], rad_Longitude, errstat)

    call tiof_put1d_r8 (tio_l2obj, "time", [0], [nstep], &
                        rad_Time, errstat)

    call tiof_put2d_r4 (tio_l2obj, "solar_zenith_angle", [0,0], &
         [nstep, nxtrack], rad_SolarZenithAngle, errstat)

    call tiof_put2d_r4 (tio_l2obj, "viewing_zenith_angle", [0,0], &
         [nstep, nxtrack], rad_ViewingZenithAngle, errstat)

    call tiof_put2d_r4 (tio_l2obj, "solar_azimuth_angle", [0,0], &
         [nstep, nxtrack], rad_SolarAzimuthAngle, errstat)

    call tiof_put2d_r4 (tio_l2obj, "viewing_azimuth_angle", [0,0], &
         [nstep, nxtrack], rad_ViewingAzimuthAngle, errstat)
 
    call tiof_put2d_r4 (tio_l2obj, "relative_azimuth_angle", [0,0], &
         [nstep, nxtrack], out_RelativeAzimuthAngle, errstat)

    call tiof_put2d_r4 (tio_l2obj, "terrain_height", [0,0], &
         [nstep, nxtrack], out_TerrainHeight, errstat)

!    call tiof_put2d_ui4 (tio_l2obj, "ground_pixel_quality_flag", [0,0], &
!         [nstep, nxtrack], rad_GroundPixelQualityFlags, errstat)

    call tiof_pop_group (tio_l2obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_geo_data: failed", errstat)
      return
    endif

  end subroutine write_geo_data


  !> Write geolocation data into L2 netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in]  l1_file    Radiance input filename
  !> @param[in]  swathname  radiance input swath name
  !> @param[in]  tio_l2obj  output file object
  !> @param[in]  nstep      along-track dimension size
  !> @param[in]  nxtrack    cross-track dimension size
  !> @param      errstat    error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   April 2021
  !-----------------------------------------------------------------------
  subroutine copy_pixel_corners (l1_file, swathname, tio_l2obj, &
       nstep, nxtrack, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l1_file, swathname
    integer (kind=4), intent(in) :: nstep, nxtrack
    !output variables
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_file_type) :: tio_l1obj
    integer, dimension(nstep) :: step_indices
    real (kind=4), dimension(4,1:nxtrack,1:nstep) :: tmp_lat, tmp_lon

    type (tiof_file_type), pointer :: tio_l2obj


    if (errstat /= 0) return

    ! open L1 rad file, read corners and dimensio indices
    call open_tio (l1_file, tio_l1obj, errstat)
    if (errstat /= 0) return

    call tiof_get1d_i4 (tio_l1obj, "mirror_step", [0], [nstep], &
         step_indices, errstat)
    call tiof_push_group (tio_l1obj, trim(swathname), errstat)
    call tiof_get3d_r4 (tio_l1obj, "latitude_bounds", [0,0,0], &
         [nstep, nxtrack, 4], tmp_lat(1:4,1:nxtrack,1:nstep), errstat)
    call tiof_get3d_r4 (tio_l1obj, "longitude_bounds", [0,0,0], &
         [nstep, nxtrack, 4], tmp_lon(1:4,1:nxtrack,1:nstep), errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "copy_pixel_corners: failed read from "//trim(l1_file), &
           errstat)
      return
    endif

    call close_tio (tio_l1obj, errstat)

    !write to output file
    tio_l2obj => primary_output_file

    ! copy mirror step indices from input radiance file
    call tiof_push_group (tio_l2obj, "/", errstat)
    call tiof_put1d_i4 (tio_l2obj, "mirror_step", [0], [nstep], step_indices,&
         errstat)
    call tiof_pop_group (tio_l2obj, errstat)
    call tiof_push_group (tio_l2obj, "geolocation", errstat)
    call tiof_put3d_r4 (tio_l2obj, "latitude_bounds", [0,0,0], &
         [nstep, nxtrack, 4], tmp_lat(1:4,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (tio_l2obj, "longitude_bounds", [0,0,0], &
         [nstep, nxtrack, 4], tmp_lon(1:4,1:nxtrack,1:nstep), errstat)

    if (errstat /= 0) then
      call tiof_pop_group (tio_l2obj, errstat)
      call tell_error (tell_io_read_error, "copy_pixel_corners: write failed",&
                       errstat)
      return
    endif

    call tiof_pop_group (tio_l2obj, errstat)

  end subroutine copy_pixel_corners


!-------------------------------
! hqw addition below 
!-------------------------------

   subroutine write_product_struct(tio_l2obj, dimlist, errstat)

      implicit none

      !input variables
    type (tiof_file_type), intent(inout) :: tio_l2obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    !output variables
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_product
    integer, dimension(2) :: dimids_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    !define r8 kind for use in setting parameter valid ranges
    integer, parameter :: r8 = kind(1.0d0)
    character (len=32) :: epoch_buf

    if (errstat /= 0) return
    
    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              ["xtrack     ", "mirror_step"], &
                              dimids_xtrack_step, &
                              errstat)

    ! Product Fields with optional attributes
    call tiof_attlist_append (att_product, errstat, "coordinates", &
                              att_text = "longitude latitude")

    call tiof_varlist_append (varlist, errstat, &
                              "CloudPressure", &
                              nf90_int, &
                              dimids = dimids_xtrack_step,  &
                              comment = "cloud pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.2e3_r8], &
                              fillvalue = fill_int, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

     call tiof_varlist_append (varlist, errstat, &
                              "EffectiveCloudFraction", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "effective cloud fraction", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

     call tiof_varlist_append (varlist, errstat, &
                              "CloudRadianceFraction466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "cloud radiance fraction at 466nm", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

     call tiof_varlist_append (varlist, errstat, &
                              "CloudRadianceFraction440", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "cloud radiance fraction at 440nm", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

    call tiof_varlist_append (varlist, errstat, &
                              "ProcessingQualityFlags", &
                              nf90_int, &
                              dimids = dimids_xtrack_step,  &
                              comment = "processing quality flags", &
                              valid_range = [0.0_r8, 65535.0_r8], &
                              fillvalue = fill_int, &
                              attlist=att_product)


    call tiof_push_group (tio_l2obj, "product", errstat)
    call tiof_def_vars (tio_l2obj, varlist, errstat)
    call tiof_pop_group (tio_l2obj, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_product)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_product_struct: failed", &
           errstat)
      return
    endif

  end subroutine write_product_struct

  !>Write product data into L2 netCDF file
  !--------------------------------------------------------------------------
  subroutine write_product_data(tio_l2obj, nstep, nxtrack, errstat)
     use m_vars, only: out_EffectiveCloudFraction, out_CloudPressure, &
            out_CloudRadianceFraction466, out_CloudRadianceFraction440, &
            out_ProcessingQualityFlags

    implicit none

    !input variables
    integer, intent(in) :: nxtrack, nstep
    !output variables
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: tio_l2obj

    if (errstat /= 0) return

    tio_l2obj => primary_output_file

    call tiof_push_group (tio_l2obj, "product", errstat)

    call tiof_put2d_i2 (tio_l2obj, "CloudPressure", [0,0], &
         [nstep, nxtrack], out_CloudPressure, errstat)

    call tiof_put2d_r4 (tio_l2obj, "EffectiveCloudFraction", [0,0], &
         [nstep, nxtrack], out_EffectiveCloudFraction, errstat)

    call tiof_put2d_r4 (tio_l2obj, "CloudRadianceFraction466", [0,0], &
         [nstep, nxtrack], out_CloudRadianceFraction466, errstat)

    call tiof_put2d_r4 (tio_l2obj, "CloudRadianceFraction440", [0,0], &
         [nstep, nxtrack], out_CloudRadianceFraction440, errstat)

    call tiof_put2d_i2 (tio_l2obj, "ProcessingQualityFlags", [0,0], &
         [nstep, nxtrack], out_ProcessingQualityFlags, errstat)

    call tiof_pop_group (tio_l2obj, errstat)

    if (errstat /= 0) then
       call tell_error (tell_io_write_error,"write_product_data:failed",errstat)
       return
    endif

   end subroutine write_product_data


  !>Write support structure into L2 netCDF file
  !----------------------------------------------------------------------------
  subroutine write_support_struct(tio_l2obj, dimlist, errstat)
    use m_vars, only: name_option_SurfaceReflectivity
    implicit none

    !input variables
    type (tiof_file_type), intent(inout) :: tio_l2obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    !output variables
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_support
    integer, dimension(2) :: dimids_xtrack_step
    character(len=255) :: name466, name440
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    !define r8 kind for use in setting parameter valid ranges
    integer, parameter :: r8 = kind(1.0d0)
    character (len=32) :: epoch_buf

    if (errstat /= 0) return

    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              ["xtrack     ", "mirror_step"], &
                              dimids_xtrack_step, &
                              errstat)

    ! Product Fields with optional attributes
    call tiof_attlist_append (att_support, errstat, "coordinates", &
                              att_text = "longitude latitude")


    call tiof_varlist_append (varlist, errstat, &
                              "ReflectanceFactor466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "466nm Reflectance=(Pi*rad466)/(irr466*cos(SZA))", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

     name466 = 'GLER466'
     name440 = 'GLER440'
     if (name_option_SurfaceReflectivity .eq. 'Kleipool') then
           name466 = 'Kleipool466'
           name440 = 'Kleipool440'
     endif
           
     call tiof_varlist_append (varlist, errstat, &
                              name466, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "466nm surface reflectivity used in calculation", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

     call tiof_varlist_append (varlist, errstat, &
                              name440, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "440nm surface reflectivity used in calculation", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SurfaceLER466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "LER at 466nm calculated at TerrainPressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SurfaceLER440", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "LER at 440nm calculated at TerrainPressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SceneLER466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "LER at 466nm calculated at ScenePressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SceneLER440", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "LER at 440nm calculated at ScenePressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "ScenePressure", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "ScenePressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.5e3_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "TerrainPressure", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "Pixel Surface Pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.5e3_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SlantColumnAmountO2O2", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "O2-O2 slant column at EffectiveTemperature", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SlantColumnReferenceO2O2", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "O2-O2 slant column at 273K", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "O2O2CloudTemperature", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "effective T for SlantColumnAmountO2O2", &
                              units = "K", &
                              valid_range = [160.0_r8, 310.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SlantColumnSceneO2O2", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "O2-O2 slant column at O2O2SceneTemperature", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SlantColumnTerrainO2O2", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "O2-O2 slant column at O2O2TerrainTemperature", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "O2O2SceneTemperature", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "effective T for SlantColumnSceneO2O2", &
                              units = "K", &
                              valid_range = [160.0_r8, 310.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "O2O2TerrainTemperature", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "effective T for SlantColumnTerrainO2O2", &
                              units = "K", &
                              valid_range = [160.0_r8, 310.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SCD_MainDataQualityFlags", &
                              nf90_int, &
                              dimids = dimids_xtrack_step,  &
                              comment = "main data quality flags for SlantColumnReferenceO2O2", &
                              valid_range = [0.0_r8, 65535.0_r8], &
                              fillvalue = fill_int, &
                              attlist=att_support)


    call tiof_push_group (tio_l2obj, "support_data", errstat)
    call tiof_def_vars (tio_l2obj, varlist, errstat)
    call tiof_pop_group (tio_l2obj, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_support)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_support_struct: failed", &
           errstat)
      return
    endif

  end subroutine write_support_struct

  !>Write support data into L2 netCDF file
  !-------------------------------------------------------------------------------
  subroutine write_support_data(tio_l2obj, nstep, nxtrack, errstat)

     use m_vars, only: out_SlantColumnAmountO2O2, nasa_SlantColumnAmountO2O2,&
               out_ReflectanceFactor, out_O2O2CloudTemperature, out_TerrainPressure,&
               out_SurfaceReflectivity440, out_SurfaceReflectivity466,&
               out_SurfaceLER440, out_SurfaceLER466,&
               out_SceneLer440, out_SceneLER466, out_ScenePressure,&
               out_SlantColumnSceneO2O2, out_O2O2SceneTemperature,&
               out_SlantColumnTerrainO2O2, out_O2O2TerrainTemperature,&
               name_option_SurfaceReflectivity

     use m_vars, only: scd_mdqfl

     implicit none

    !input variables
    integer, intent(in) :: nxtrack, nstep
    !output variables
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: tio_l2obj

    !local variables
    character(len=255) :: name466, name440

    if (errstat /= 0) return

    tio_l2obj => primary_output_file

    name466 = 'GLER466'
    name440 = 'GLER440'
    if (name_option_SurfaceReflectivity .eq. 'Kleipool') then
        name466 = 'Kleipool466'
        name440 = 'Kleipool440'
    endif

    call tiof_push_group (tio_l2obj, "support_data", errstat)

    call tiof_put2d_r4 (tio_l2obj, "ReflectanceFactor466", [0,0], &
         [nstep, nxtrack], out_ReflectanceFactor, errstat)

    call tiof_put2d_r4 (tio_l2obj, name466, [0,0], &
         [nstep, nxtrack], out_SurfaceReflectivity466, errstat)

    call tiof_put2d_r4 (tio_l2obj, name440, [0,0], &
         [nstep, nxtrack], out_SurfaceReflectivity440, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SurfaceLER466", [0,0], &
         [nstep, nxtrack], out_SurfaceLER466, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SurfaceLER440", [0,0], &
         [nstep, nxtrack], out_SurfaceLER440, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SceneLER466", [0,0], &
         [nstep, nxtrack], out_SceneLER466, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SceneLER440", [0,0], &
         [nstep, nxtrack], out_SceneLER440, errstat)

    call tiof_put2d_r4 (tio_l2obj, "ScenePressure", [0,0], &
         [nstep, nxtrack], out_ScenePressure, errstat)

    call tiof_put2d_r4 (tio_l2obj, "TerrainPressure", [0,0], &
         [nstep, nxtrack], out_TerrainPressure, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SlantColumnAmountO2O2", [0,0], &
         [nstep, nxtrack], out_SlantColumnAmountO2O2, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SlantColumnReferenceO2O2", [0,0], &
         [nstep, nxtrack], nasa_SlantColumnAmountO2O2, errstat)

    call tiof_put2d_r4 (tio_l2obj, "O2O2CloudTemperature", [0,0], &
         [nstep, nxtrack], out_O2O2CloudTemperature, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SlantColumnSceneO2O2", [0,0], &
         [nstep, nxtrack], out_SlantColumnSceneO2O2, errstat)

    call tiof_put2d_r4 (tio_l2obj, "O2O2SceneTemperature", [0,0], &
         [nstep, nxtrack], out_O2O2SceneTemperature, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SlantColumnTerrainO2O2", [0,0], &
         [nstep, nxtrack], out_SlantColumnTerrainO2O2, errstat)

    call tiof_put2d_r4 (tio_l2obj, "O2O2TerrainTemperature", [0,0], &
         [nstep, nxtrack], out_O2O2TerrainTemperature, errstat)

    call tiof_put2d_i2 (tio_l2obj,"SCD_MainDataQualityFlags", [0,0], &
         [nstep, nxtrack], scd_mdqfl, errstat)

    call tiof_pop_group(tio_l2obj, errstat)

    if (errstat /= 0) then 
       call tell_error (tell_io_write_error,"write_support_data: failed",errstat)
       return
    endif

  end subroutine write_support_data

!-----------------
   subroutine write_tio_glbattr(tio_l2obj,errstat)

    use m_vars, only: gmetadata
    use netcdf, only: nf90_global, nf90_put_att

    implicit none

    include 'GetConfig.inc'

    type (tiof_file_type), pointer, intent(in) :: tio_l2obj
    integer (kind=4), intent(inout) :: errstat

    character(len=CFG_VAL_LEN) :: buf, attname
    integer (kind=4) :: tmpint
    real (kind=4) :: tmpreal

    if (errstat /= 0) return

    attname = 'StartDate' 
    buf = gmetadata%startdate
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'StartTime' 
    buf = gmetadata%starttime
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'EndDate' 
    buf = gmetadata%enddate
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'EndTime' 
    buf = gmetadata%endtime
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'scan_num' 
    tmpint = gmetadata%scan_num
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,tmpint)

    attname = 'granule_num' 
    tmpint = gmetadata%granule_num
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,tmpint)

    attname = 'Collection' 
    buf = gmetadata%omi_collection
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'Instrument' 
    buf = gmetadata%InstrumentName
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'platform' 
    buf = gmetadata%platformShortName
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'ProcessingCenter' 
    buf = gmetadata%ProcessingCenter
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'PricessingLevel'
    buf = gmetadata%ProcessingLevel
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'LeadScientist' 
    buf = gmetadata%LeadScientist
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'author_name' 
    buf = gmetadata%author_name
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

    attname = 'author_affiliation' 
    buf = gmetadata%author_affiliation
    errstat=nf90_put_att(tio_l2obj%fileid, nf90_global,attname,buf)

   end subroutine write_tio_glbattr

end module m_write_output_tio
