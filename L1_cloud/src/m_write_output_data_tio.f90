!>Subroutines to write out L2 Cloud netCDF file
module m_write_output_data_tio
  use netcdf
  use tio_module
  use tell_module
  use cld_names_module
  use m_vars, only: fill_value, fill_value_int

  implicit none
  private

  public create_output_file, close_output_file, write_coordinate_vars, &
       write_geo_struct, write_geo_data, write_cloud_struct, &
       copy_hdr_metadata, copy_pixel_corners, label_output_file!, write_cloud_data

  type (tiof_file_type), private, target :: primary_output_file

  !define fill values in correct format
  real (kind=8), private, parameter :: fill_float=fill_value, &
       fill_double=fill_value, fill_short=fill_value_int, &
       fill_int=fill_value_int

contains

  !>Write coordinate variables into L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] obj file object to be written into
  !> @param[in] dimlist list of dimension parameters
  !> @param[in] num_steps size of dimension in scan direction
  !> @param[in] num_xtrack size of dimension across scan direction
  !> @param errstat error tracking code, non-zero indiactes problem
  !> @param[in] num_wavel size of spectral dimension
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
  subroutine write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, &
       errstat, num_wavel)
    use m_vars, only: write_resid
    implicit none
    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(in) :: num_steps, num_xtrack
    integer, intent(in), optional :: num_wavel
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer, dimension(num_xtrack) :: xtrack_indices
    integer, dimension(num_steps) :: step_indices
    integer, dimension(:), allocatable :: wavel_indices
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_wavel_xtrack_step
    integer :: i

    if (errstat /= 0) return

    if (present(num_wavel)) allocate(wavel_indices(num_wavel), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "write_coordinate_vars: unable to allocate num_wavel", &
           errstat)
      return
    endif

    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              [cld_dim_xtrack, cld_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    if (write_resid) then
      call tiof_dimlist_lookup (dimlist, &
                             [cld_dim_channel, cld_dim_xtrack, cld_dim_step], &
                             dimids_wavel_xtrack_step, &
                             errstat)
    endif

    ! Make a list of variables with their dimension ids and attributes:

    ! netcdf coordinate variables:
    call tiof_varlist_append (varlist, errstat, cld_dim_xtrack, nf90_int, &
                             dimids=[dimids_xtrack_step(1)])
    call tiof_varlist_append (varlist, errstat, cld_dim_step, nf90_int, &
                             dimids=[dimids_xtrack_step(2)])
    if (write_resid) then
      call tiof_varlist_append (varlist, errstat, cld_dim_channel, nf90_int, &
                             dimids=[dimids_wavel_xtrack_step(1)])
    endif
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)

    ! FIXME: eventually, this will be something like
    ! step_indices=[mirror_step_beg, ..., mirror_step_end]
    ! where mirror_step_beg/end are granule-specific

    step_indices = [(i, i=0,num_steps-1)]
    call tiof_put1d_i4 (obj, cld_dim_step, [0], [num_steps], &
         step_indices, errstat)

    xtrack_indices = [(i, i=0,num_xtrack-1)]
    call tiof_put1d_i4 (obj, cld_dim_xtrack, [0], [num_xtrack], &
         xtrack_indices, errstat)

    if (present(num_wavel)) then
      wavel_indices = [(i, i=0,num_wavel-1)]
      call tiof_put1d_i4 (obj, cld_dim_channel, [0], [num_wavel], &
         wavel_indices, errstat)
    endif

  end subroutine write_coordinate_vars

  subroutine copy_pixel_corners (l1bfile, num_steps, num_xtrack, errstat)
    use m_vars, only: nc_swathname
    implicit none
    character (len=*), intent(in) :: l1bfile
    integer (kind=4), intent(in) :: num_steps, num_xtrack
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_file_type) :: l1b
    integer, dimension(num_steps) :: step_indices
    real (kind=4), dimension(4,1:num_xtrack,1:num_steps) :: tmp

    if (errstat /= 0) return

    call tiof_open (l1bfile, l1b, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       "copy_pixel_corners: opening file "//trim(l1bfile), &
                       errstat)
      return
    endif

    call tiof_get1d_i4 (l1b, cld_dim_step, [0], [num_steps], step_indices, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       "copy_pixel_corners: reading step indices from file "//trim(l1bfile), &
                       errstat)
      return
    endif

    ! use pixel corners associated with the radiances being fitted
    call tiof_push_group (l1b, trim(nc_swathname), errstat)

    call tiof_get3d_r4 (l1b, cld_var_latitude_bounds, [0,0,0], [num_steps, num_xtrack, 4], &
                        tmp(1:4,1:num_xtrack,1:num_steps), errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
                       "copy_pixel_corners: reading latitude bounds from file "//trim(l1bfile), &
                       errstat)
      return
    endif

    obj => primary_output_file

    ! copy mirror step indices from input radiance file
    call tiof_push_group (obj, "/", errstat)
    call tiof_put1d_i4 (obj, cld_dim_step, [0], [num_steps], step_indices, errstat)
    call tiof_pop_group (obj, errstat)

    call tiof_push_group (obj, cld_grp_geolocation, errstat)
    call tiof_put3d_r4 (obj, cld_var_latitude_bounds, [0,0,0], [num_steps, num_xtrack, 4], &
                        tmp(1:4,1:num_xtrack,1:num_steps), errstat)
    call tiof_get3d_r4 (l1b, cld_var_longitude_bounds, [0,0,0], [num_steps, num_xtrack, 4], &
                        tmp(1:4,1:num_xtrack,1:num_steps), errstat)
    call tiof_put3d_r4 (obj, cld_var_longitude_bounds, [0,0,0], [num_steps, num_xtrack, 4], &
                        tmp(1:4,1:num_xtrack,1:num_steps), errstat)

    call tiof_close (l1b, errstat)
    if (errstat /= 0) then
      call tiof_pop_group (obj, errstat)
      call tell_error (tell_io_read_error, "copy_pixel_corners: reading file "//trim(l1bfile), &
                       errstat)
      return
    endif

    call tiof_pop_group (obj, errstat)

  end subroutine copy_pixel_corners

  !>Top-level subroutine to create and populate an L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] outfile_nc name of L2 netCDF file to be written
  !> @param[in] num_steps size of dimension in scan direction
  !> @param[in] num_xtrack size of dimension across scan direction
  !> @param errstat error tracking code, non-zero indiactes problem
  !> @param[in] num_wavel size of spectral dimension
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
  subroutine create_output_file (outfile_nc, num_steps, num_xtrack, &
       errstat, num_wavel)
    use m_vars, only: write_resid
    implicit none
    character (len=*), intent(in) :: outfile_nc
    integer (kind=4), intent(in) :: num_xtrack, num_steps
    integer (kind=4), intent(in), optional :: num_wavel
    integer (kind=4), intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_dimlist_type) :: dimlist

    if (errstat /= 0) return

    obj => primary_output_file

    ! create the file
    call tiof_create (obj, outfile_nc, nf90_clobber, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: creating file "//trim(outfile_nc), &
           errstat)
      return
    endif

    call tiof_put_git_commit_hash (obj, errstat)
    call tiof_write_epoch_timestamp (obj, errstat)

    ! Create default groups.
    call tiof_def_group (obj, cld_grp_product, errstat)
    call tiof_def_group (obj, cld_grp_geolocation, errstat)
    call tiof_def_group (obj, cld_grp_support_data, errstat)
    call tiof_def_group (obj, cld_grp_qa_stats, errstat)
    if (write_resid) then
      call tiof_def_group (obj, cld_grp_diagnostic, errstat)
    endif
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "create_output_file:  defining groups in "//trim(outfile_nc), &
                       errstat)
      return
    endif

    ! define the dimension list
    call tiof_dimlist_append (dimlist, cld_dim_step, num_steps, errstat)
    call tiof_dimlist_append (dimlist, cld_dim_xtrack, num_xtrack, errstat)
    call tiof_dimlist_append (dimlist, cld_dim_corner, 4, errstat)
    if (write_resid) then
      call tiof_dimlist_append (dimlist, cld_dim_channel, num_wavel, errstat)
    endif
    call tiof_def_dims (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: defining dimensions in "//trim(outfile_nc), &
           errstat)
      return
    endif

    ! coordinate variables
    if (write_resid) then
      call write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, errstat, &
           num_wavel)
    else
      call write_coordinate_vars (obj, dimlist, num_steps, num_xtrack, errstat)
    endif
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: writing coordinate variables to "//trim(outfile_nc), &
           errstat)
      return
    endif

    ! geolocation & cloud variable definitions
    call write_geo_struct (obj, dimlist, errstat)
    call write_cloud_struct (obj, dimlist, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: writing variables to "//trim(outfile_nc), &
           errstat)
      return
    endif

    ! geolocation & cloud data
    call write_geo_data (obj, num_steps, num_xtrack, errstat)
    if (write_resid) then
      call write_cloud_data (obj, num_steps, num_xtrack, errstat, num_wavel)
    else
      call write_cloud_data (obj, num_steps, num_xtrack, errstat)
    endif
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: writing data to "//trim(outfile_nc), &
           errstat)
      return
    endif

    !metadata
    call write_metadata (errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: writing metadata to "//trim(outfile_nc), &
           errstat)
      return
    endif

    !Free dimension list
    call tiof_dimlist_free (dimlist)

    !!Close the netCDF file
    !call close_output_file (errstat)
    !if (errstat /= 0) then
    !  call tell_error (tell_io_error, &
    !       "create_output_file: unable to close file "//trim(outfile_nc), &
    !       errstat)
    !  return
    !endif

  end subroutine create_output_file

  !>Close L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param errstat error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
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

  !> Create the structure for the geolocation data in L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] obj file object to be written into
  !> @param[in] dimlist list of dimension parameters
  !> @param errstat error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
  subroutine write_geo_struct(obj, dimlist, errstat)
    use m_vars, only: write_resid
    implicit none

    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_geo, att_latbnd, att_lonbnd
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_wavel_xtrack_step
    integer, dimension(3) :: dimids_corner_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    !define r8 kind for use in setting parameter valid ranges
    integer, parameter :: r8 = kind(1.0d0)
    character (len=32) :: epoch_buf

    if (errstat /= 0) return

    ! Define dimid arrays associated with common data field shapes.
    call tiof_dimlist_lookup (dimlist, &
                              [cld_dim_xtrack, cld_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [cld_dim_corner, cld_dim_xtrack, cld_dim_step], &
                              dimids_corner_xtrack_step, &
                              errstat)
    if (write_resid) then
      call tiof_dimlist_lookup (dimlist, &
                             [cld_dim_channel, cld_dim_xtrack, cld_dim_step], &
                             dimids_wavel_xtrack_step, &
                             errstat)
    endif

    ! Make a list of variables with their dimension ids and attributes:
    epoch_buf(:)=''
    call tiof_mktimestamp_str (0.0_r8, epoch_buf, errstat)

    ! Geolocation Fields with optional attribute lists
    call tiof_attlist_append (att_geo, errstat, "coordinates", &
                              att_text = trim(cld_var_longitude) &
                              //' '//trim(cld_var_latitude))
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_time, &
                              nf90_double, &
                              dimids = [dimids_xtrack_step(2)],  &
                              comment = "radiance exposure start time", &
                              units = "seconds since "//trim(epoch_buf), &
                              valid_range = [0.0_r8, 1.0e30_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              fillvalue = fill_double)
    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = cld_var_latitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_latitude, &
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
                              att_text = cld_var_longitude_bounds)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_longitude, &
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
                              cld_var_latitude_bounds, &
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
                              cld_var_longitude_bounds, &
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
                              cld_var_sz_angle, &
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
                              cld_var_vz_angle, &
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
                              cld_var_ra_angle, &
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
                              cld_var_terr_height, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "terrain height", &
                              units = "m", &
                              valid_range = [-1000.0_r8, 10000.0_r8], &
                              fillvalue = fill_short, &
                              attlist=att_geo)
    ! note that using a nf90_short causes the flags to be written out wrong
    ! should be nf90_ushort, but the fill value is unsupported...
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_gpqf, &
                              nf90_uint, &
                              dimids = dimids_xtrack_step,  &
                              comment = "ground pixel quality flag", &
                              valid_range = [0.0_r8, 65535.0_r8], &
                              fillvalue = fill_int, &
                              attlist=att_geo)
    call tiof_push_group (obj, cld_grp_geolocation, errstat)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_pop_group (obj, errstat)
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

  !> Write geolocation data into L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] obj file object to be written into
  !> @param[in] num_steps size of dimension in scan direction
  !> @param[in] num_xtrack size of dimension across scan direction
  !> @param errstat error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
  subroutine write_geo_data(obj, num_steps, num_xtrack, errstat)
    use m_vars, only: lat, lon, time, sza, sat_zen, azimuth, terr_height, &
         geoflg

    implicit none

    integer, intent(in) :: num_xtrack, num_steps
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_push_group (obj, cld_grp_geolocation, errstat)

    call tiof_put2d_r4 (obj, cld_var_latitude, [0,0], &
         [num_steps, num_xtrack], lat(1:num_xtrack,1:num_steps), errstat)

!    where(lon(:,:) > 180.d0) lon(:,:) = lon(:,:) -360.d0
    call tiof_put2d_r4 (obj, cld_var_longitude, [0,0], &
         [num_steps, num_xtrack], lon(1:num_xtrack,1:num_steps), errstat)

    call tiof_put1d_r8 (obj, cld_var_time, [0], [num_steps], &
                        time (1:num_steps), errstat)

    call tiof_put2d_r4 (obj, cld_var_sz_angle, [0,0], &
         [num_steps, num_xtrack], sza(0:num_xtrack-1,1:num_steps), errstat)

    call tiof_put2d_r4 (obj, cld_var_vz_angle, [0,0], &
         [num_steps, num_xtrack], sat_zen(0:num_xtrack-1,1:num_steps), &
         errstat)

    call tiof_put2d_r4 (obj, cld_var_ra_angle, [0,0], &
         [num_steps, num_xtrack], azimuth(0:num_xtrack-1,1:num_steps), &
         errstat)

    call tiof_put2d_i2 (obj, cld_var_terr_height, [0,0], &
         [num_steps, num_xtrack], terr_height(1:num_xtrack,1:num_steps), &
         errstat)

    call tiof_put2d_ui4 (obj, cld_var_gpqf, [0,0], &
         [num_steps, num_xtrack], geoflg(1:num_xtrack,1:num_steps), &
         errstat)

    call tiof_pop_group (obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_geo_data: failed", errstat)
      return
    endif

  end subroutine write_geo_data

  !> Create cloud data structure in L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] obj file object to be written into
  !> @param[in] dimlist list of dimension parameters
  !> @param errstat error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
  subroutine write_cloud_struct(obj, dimlist, errstat)
    use m_vars, only: squeeze, write_fill, write_resid, cal_reflec, &
         do_mler, do_cloud_mask, cloud_mask, have_omi_data

    implicit none

    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_cld
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_wavel_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    !define r8 kind for use in setting parameter valid ranges
    integer, parameter :: r8 = kind(1.0d0)

    if (errstat /= 0) return

    ! Define dimid arrays associated with common data field shapes
    if (write_resid) then
      call tiof_dimlist_lookup (dimlist, &
                         [cld_dim_channel, cld_dim_xtrack, cld_dim_step], &
                         dimids_wavel_xtrack_step, &
                         errstat)
    endif
    call tiof_dimlist_lookup (dimlist, &
                              [cld_dim_xtrack, cld_dim_step], &
                              dimids_xtrack_step, &
                              errstat)

    ! Cloud processing fields with optional attribute lists
    call tiof_attlist_append (att_cld, errstat, "coordinates", &
                              att_text = trim(cld_var_longitude) &
                              //' '//trim(cld_var_latitude))

    !Product group
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_cld_press, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "cloud pressure for O3", &
                              units = "hPa", &
                              valid_range = [-1000.0_r8, 10000.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_cld_frac, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "cloud fraction for O3", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_proc_qf, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "processing quality flags", &
                              valid_range = [0.0_r8, 32767.0_r8], &
                              fillvalue = fill_short, &
                              attlist=att_cld)
    if (do_cloud_mask .and. allocated(cloud_mask)) then
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_cloud_mask, &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              comment = "cloud mask", &
                              valid_range = [0.0_r8, 3.0_r8], &
                              fillvalue = fill_short, &
                              attlist=att_cld)
    endif

    call tiof_push_group (obj, cld_grp_product, errstat)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist)

    !Support data group
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_chlorophyll, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "chlorophyll concentration", &
                              units = "mg/m3", &
                              valid_range = [0.0_r8, 1000.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_terr_press, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "terrain pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 2000.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_rad_cld_frac, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "radiative cloud fraction", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_wav_shift, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "wavelength shift", &
                              units = "nm", &
                              valid_range = [-100.0_r8, 100.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_surface_reflec, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "surface reflectivity climatology", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_reflec, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "reflectivity", &
                              valid_range = [-10.0_r8, 10.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    if (have_omi_data) then
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_mqf, &
                              nf90_short, &
                              dimids = [dimids_xtrack_step(2)],  &
                              comment = "measurement quality flags", &
                              valid_range = [0.0_r8, 65536.0_r8], &
                              fillvalue = fill_short, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    endif
    if (squeeze) then
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_wav_squeeze, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "wavelength squeeze", &
                              valid_range = [0.0_r8, 10.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    endif
    if (write_fill) then
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_filling_in, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "effective filling in", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    endif
    if (.not. do_mler) then
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_cld_reflec, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "cloud reflectivity", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    endif
    if (cal_reflec) then
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_didr, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "radiance (fractional) refl. sens.", &
                              valid_range = [0.0_r8, 100.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    endif

    call tiof_push_group (obj, cld_grp_support_data, errstat)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist)

    !QA_statistics group
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_resid_bias, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "residual bias", &
                              valid_range = [-1000.0_r8, 1000.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_resid_stddev, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "residual standard deviation", &
                              valid_range = [-1000.0_r8, 1000.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)
    call tiof_varlist_append (varlist, errstat, &
                              cld_var_convergence_factor, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "convergence factor", &
                              valid_range = [0.0_r8, 1e38_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)

    call tiof_push_group (obj, cld_grp_qa_stats, errstat)
    call tiof_def_vars (obj, varlist, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_varlist_free (varlist)

    !Diagnostic group
    if (write_resid) then
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_wav_resid, &
                              nf90_float, &
                              dimids = [dimids_wavel_xtrack_step(1)],  &
                              comment = "residual wavelengths", &
                              units = "nm", &
                              valid_range = [0.0_r8, 1000.0_r8], &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              fillvalue = fill_float)
      call tiof_varlist_append (varlist, errstat, &
                              cld_var_resid, &
                              nf90_float, &
                              dimids = dimids_wavel_xtrack_step,  &
                              comment = "radiance residual", &
                              units = "percent", &
                              valid_range = [-100.0_r8, 100.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_cld)

      call tiof_push_group (obj, cld_grp_diagnostic, errstat)
      call tiof_def_vars (obj, varlist, errstat)
      call tiof_pop_group (obj, errstat)
      call tiof_varlist_free (varlist)
    endif

    call tiof_attlist_free (att_cld)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_cloud_struct: failed", &
           errstat)
      return
    endif
  end subroutine write_cloud_struct

  !>Write out cloud data into L2 Cloud netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] obj file object to be written into
  !> @param[in] num_steps size of dimension in scan direction
  !> @param[in] num_xtrack size of dimension across scan direction
  !> @param errstat error tracking code, non-zero indiactes problem
  !> @param[in] num_wavel size of spectral dimension
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
  subroutine write_cloud_data(obj, num_steps, num_xtrack, errstat, &
       num_wavel)
    use m_vars, only: squeeze, write_fill, write_resid, cal_reflec, &
         do_mler, cloud_mask, rad_cld_frac, ps, shifts2, squeezes, &
         ref_clr, fill, wave_resid, resid, dIdR, reflect_cld, &
         refl, meas_qual_flg, biases2, stds2, chi_sqr2, chlorophyll, &
         cld_pres2, eff_cld_frac2, qc2, do_cloud_mask, read_he4, &
         have_omi_data

    implicit none

    integer, intent(in) :: num_xtrack, num_steps
    integer, intent(in), optional :: num_wavel
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    if (errstat /= 0) return

    obj => primary_output_file

    !Product group
    call tiof_push_group (obj, cld_grp_product, errstat)
    if (do_cloud_mask .and. allocated(cloud_mask)) then
      call tiof_put2d_i2 (obj, cld_var_cloud_mask, [0,0], &
         [num_steps, num_xtrack], cloud_mask(1:num_xtrack,1:num_steps), &
         errstat)
    endif

    call tiof_put2d_r4 (obj, cld_var_cld_press, [0,0], &
         [num_steps,num_xtrack], cld_pres2(0:num_xtrack-1,1:num_steps), &
         errstat)

    call tiof_put2d_r4 (obj, cld_var_cld_frac, [0,0], &
         [num_steps,num_xtrack], eff_cld_frac2(0:num_xtrack-1,1:num_steps), &
         errstat)

    call tiof_put2d_i2 (obj, cld_var_proc_qf, [0,0], &
         [num_steps,num_xtrack], qc2(0:num_xtrack-1,1:num_steps), errstat)
    call tiof_pop_group (obj, errstat)

    !Support_data group
    call tiof_push_group (obj, cld_grp_support_data, errstat)

    call tiof_put2d_r4 (obj, cld_var_chlorophyll, [0,0], &
         [num_steps,num_xtrack], chlorophyll(0:num_xtrack-1,1:num_steps), &
         errstat)

    if (.not. read_he4) ps=ps*1013.25
    call tiof_put2d_r4 (obj, cld_var_terr_press, [0,0], &
         [num_steps,num_xtrack], ps(0:num_xtrack-1,1:num_steps), errstat)

    call tiof_put2d_r4 (obj, cld_var_rad_cld_frac, [0,0], &
         [num_steps,num_xtrack], rad_cld_frac(0:num_xtrack-1,1:num_steps), &
         errstat)

    call tiof_put2d_r4 (obj, cld_var_surface_reflec, [0,0], &
         [num_steps,num_xtrack], ref_clr(0:num_xtrack-1,1:num_steps), errstat)

    call tiof_put2d_r4 (obj, cld_var_reflec, [0,0], &
         [num_steps,num_xtrack], refl(0:num_xtrack-1,1:num_steps), errstat)

    if (have_omi_data) then
      call tiof_put1d_i2 (obj, cld_var_mqf, [0], &
           [-1], meas_qual_flg(1:num_steps), errstat)
    endif

    call tiof_put2d_r4 (obj, cld_var_wav_shift, [0,0], &
         [num_steps,num_xtrack], shifts2(0:num_xtrack-1,1:num_steps), errstat)

    if (squeeze) then
      call tiof_put2d_r4 (obj, cld_var_wav_squeeze, [0,0], &
         [num_steps,num_xtrack], squeezes(0:num_xtrack-1,1:num_steps), errstat)
    endif

    if (write_fill) then
      call tiof_put2d_r4 (obj, cld_var_filling_in, [0,0], &
         [num_steps,num_xtrack], fill(0:num_xtrack-1,1:num_steps), errstat)
    endif

    if (.not. do_mler) then
      call tiof_put2d_r4 (obj, cld_var_cld_reflec, [0,0], &
         [num_steps,num_xtrack], reflect_cld(0:num_xtrack-1,1:num_steps), &
         errstat)
    endif

    if (cal_reflec) then
      call tiof_put2d_r4 (obj, cld_var_didr, [0,0], &
         [num_steps,num_xtrack], dIdR(0:num_xtrack-1,1:num_steps), errstat)
    endif

    call tiof_pop_group (obj, errstat)

    !QA_statistics group
    call tiof_push_group (obj, cld_grp_qa_stats, errstat)
    call tiof_put2d_r4 (obj, cld_var_resid_bias, [0,0], &
         [num_steps,num_xtrack], biases2(0:num_xtrack-1,1:num_steps), errstat)

    call tiof_put2d_r4 (obj, cld_var_resid_stddev, [0,0], &
         [num_steps,num_xtrack], stds2(0:num_xtrack-1,1:num_steps), errstat)

    call tiof_put2d_r4 (obj, cld_var_convergence_factor, [0,0], &
         [num_steps,num_xtrack], chi_sqr2(0:num_xtrack-1,1:num_steps), errstat)

    call tiof_pop_group (obj, errstat)

    !Diagnostic group
    if (write_resid) then
      call tiof_push_group (obj, cld_grp_diagnostic, errstat)
      call tiof_put1d_r4 (obj, cld_var_wav_resid, [0], &
         [num_wavel], wave_resid(1:num_wavel), errstat)
      call tiof_put3d_r4 (obj, cld_var_resid, [0,0,0], &
         [num_steps,num_xtrack,num_wavel], &
         resid(1:num_wavel,1:num_xtrack,1:num_steps), errstat)
      call tiof_pop_group (obj, errstat)
    endif

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_cloud_data: failed", errstat)
      return
    endif

  end subroutine write_cloud_data

  !>Add some basic metadata to L2 Cloud netCDF file.
  !>Proof of concept for now
  !-----------------------------------------------------------------------
  !
  !> @param errstat error tracking code, non-zero indiactes problem
  !
  !> @author E. O'Sullivan   March 2015
  !-----------------------------------------------------------------------
  subroutine write_metadata (errstat)
    use m_vars, only: qc, cloud_pres, eff_cld_frac, cld_frac_min, &
         n_good_input, n_good_output, n_missing, n_input, &
         highqual, badqual

    implicit none
    integer, intent(inout) :: errstat
    integer :: Qamissingdata , Qaboundsdata, QAPercentCloudCover, &
         PerGoodQualData, ind
    character(len=*), parameter :: expl="Flag set to Passed if "// &
         "QAPercentHighQualityData >= 80%, "// &
         "Flag set to Suspect if percent high quality data >= 20%, "//&
         "or L1B AutomaticQualityFlag not set to Passed, "//         &
         "otherwise Flag set to Failed"
    character(len=10) :: value

    type (tiof_file_type), pointer :: obj
    type (tiof_attlist_type) :: attlist

    obj => primary_output_file

    if (errstat /= 0) return

    ! calculate quality stats
    ind = count(btest(qc(:,:),2) .or. btest(qc(:,:),3))
    QAboundsdata = nint( real(ind) / real(size(cloud_pres))*100.0)

    ind = count(eff_cld_frac > cld_frac_min)
    QAPercentCloudCover= nint( real(ind) / real(size(eff_cld_frac))*100.0)

    PerGoodQualData = nint( real(n_good_output)*100.0 / real(n_good_input))
    if(PerGoodQualData >= highqual ) then
      value = "Passed"
    else if(PerGoodQualData >= badqual) then
      value = "Suspect"
    else
      value = "Failed"
    endif
    !if( trim(value) == "Passed" .and. trim(L1B_AutQualFl) /= "Passed") &
    !value = "Suspect"

    QAmissingdata = nint( real(n_missing)*100.0 /2.0/ real(n_input))
    ! n_missing counts twice in m_cloud_pres_ret

    call tiof_attlist_append (attlist, errstat, "QA_percent_missing_data", &
         att_i4=[QAmissingdata])
    call tiof_attlist_append (attlist, errstat, "QA_percent_cloud_cover", &
         att_i4=[QAPercentCloudCover])
    call tiof_attlist_append (attlist, errstat, "QA_percent_out_of_bounds", &
         att_i4=[QAboundsdata])
    call tiof_attlist_append (attlist, errstat, "automatic_quality_flag", &
         att_text=value)
    call tiof_attlist_append (attlist, errstat, &
         "automatic_quality_flag_explanation", att_text=expl)

    call tiof_push_group (obj, cld_grp_qa_stats, errstat)
    call tiof_def_atts (obj, attlist, nf90_global, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_metadata: failed", &
           errstat)
      return
    endif

  end subroutine write_metadata

  subroutine copy_hdr_metadata (l1bfile, errstat)
    implicit none
    character (len=*), intent(in) :: l1bfile
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj
    type (tiof_file_type) :: l1b

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_open (l1bfile, l1b, nf90_nowrite, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_open_error, &
                       "copy_hdr_metadata: opening file "//trim(l1bfile), &
                       errstat)
      return
    endif

    call tiof_copy_granule_ident (l1b, obj, errstat)
    call tiof_close (l1b, errstat)

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, &
                       "copy_hdr_metadata: copying from "//trim(l1bfile), &
                       errstat)
    endif

  end subroutine copy_hdr_metadata

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

end module m_write_output_data_tio
