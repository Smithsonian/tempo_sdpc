!> Write output to a TEMPO-format netCDF4 file
module m_write_output_tio
  use tio_module
  use tell_module
  use netcdf, only: nf90_clobber, nf90_write, &
    nf90_double, nf90_float, nf90_short, nf90_uint, nf90_int
  use m_read_input_tio, only: open_tio, close_tio
  use m_vars, only: iFillValue, fFillValue

  private write_coordinate_vars, write_geo_struct, write_geo_data, &
       copy_pixel_corners
  private write_product_struct, write_product_data
  private write_support_struct, write_support_data

  public update_output_file_tio
  public create_output_file_tio

  type (tiof_file_type), private, target :: primary_output_file

  !fill values
  real (kind=8), private, parameter :: fill_bit = -128, &
       fill_short = iFillValue, fill_int = iFillValue, fill_float = fFillValue, &
       fill_double = fill_float, fill_float_nines = -9999.0

contains

  !>Top-level subroutine to update an L2 O2O2 netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param[in] outfile    name of L2 netCDF file to be updated
  !> @param[in] nstep      along-track dimension size
  !> @param[in] nxtrack    cross-track dimension size
  !> @param     errstat    error tracking code, non-zero indicates problem
  !
  !> @author John Houck  Oct 2021
  !-----------------------------------------------------------------------
  subroutine update_output_file_tio (outfile, nstep, nxtrack, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: outfile
    integer (kind=4), intent(in) :: nstep, nxtrack
    !output variables
    integer (kind=4), intent(inout) :: errstat
    !local variables
    type (tiof_file_type), pointer :: tio_l2obj
    integer :: dimid_xtrack, dimid_step

    if (errstat /= 0) return

    tio_l2obj => primary_output_file

    ! Open the file
    call tiof_open (trim(adjustl(outfile)), tio_l2obj, nf90_write, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "update_output_file_tio: opening file "//trim(outfile), &
           errstat)
      return
    endif

    call tiof_history_append_cmdline (tio_l2obj)

    call tiof_inq_dimid (tio_l2obj, 'xtrack', dimid_xtrack, errstat)
    call tiof_inq_dimid (tio_l2obj, 'mirror_step', dimid_step, errstat)

    !hqw addition -------------------------------------------------
    ! product variable definitions
    call write_product_struct (tio_l2obj, dimid_xtrack, dimid_step, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write structures", &
           errstat)
      return
    endif

    ! product data
    call write_product_data (tio_l2obj, nstep, nxtrack, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write data", &
           errstat)
      return
    endif

    ! support variable definitions
    call write_support_struct(tio_l2obj, dimid_xtrack, dimid_step, errstat)
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

    call close_tio (tio_l2obj, errstat)

  end subroutine update_output_file_tio

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
    integer, dimension(2) :: dimids_xtrack_step

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
           "create_output_file: failed to write geolocation structures", &
           errstat)
      return
    endif

    ! geolocation data
    call write_geo_data (tio_l2obj, nstep, nxtrack, errstat)
    call copy_pixel_corners (l1_file, swathname, tio_l2obj, nstep, nxtrack, &
         errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write geolocation data", &
           errstat)
      return
    endif

    !hqw addition -------------------------------------------------
    ! product variable definitions
    call tiof_dimlist_lookup (dimlist, &
                              ["xtrack     ", "mirror_step"], &
                              dimids_xtrack_step, &
                              errstat)
    call write_product_struct (tio_l2obj, dimids_xtrack_step(1), &
                               dimids_xtrack_step(2), errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write product structures", &
           errstat)
      return
    endif

    ! product data
    call write_product_data (tio_l2obj, nstep, nxtrack, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed to write product data", &
           errstat)
      return
    endif

    ! support variable definitions
    call write_support_struct(tio_l2obj, dimids_xtrack_step(1), &
                              dimids_xtrack_step(2), errstat)
    if (errstat /= 0) then
       call tell_error (tell_io_write_error, &
            "creat_output_file: failed to write support structures", &
            errstat)
       return
    endif

    ! support data
    call write_support_data(tio_l2obj, nstep, nxtrack, errstat)
    if (errstat /= 0) then
       call tell_error (tell_io_write_error, &
            "create_output_file: filed to write support data", &
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
    integer, parameter :: deflate_level = 1
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
                              long_name = "radiance exposure start time", &
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
                              long_name = "latitude at pixel center", &
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
                              long_name= "longitude at pixel center", &
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
                              long_name = "solar zenith angle at pixel center", &
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
                              long_name = "solar azimuth angle at pixel center", &
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
                              long_name = "viewing zenith angle at pixel center", &
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
                              long_name = "viewing azimuth angle at pixel center", &
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
                              long_name = "relative azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_r8, 180.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
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
         rad_SolarZenithAngle, &
         rad_SolarAzimuthAngle, &
         rad_ViewingZenithAngle, &
         rad_ViewingAzimuthAngle, &
         rad_RelativeAzimuthAngle 

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
         [nstep, nxtrack], rad_RelativeAzimuthAngle, errstat)

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

   subroutine write_product_struct(tio_l2obj, dimid_xtrack, dimid_step, errstat)

      implicit none

      !input variables
    type (tiof_file_type), intent(inout) :: tio_l2obj
    integer, intent(in) :: dimid_xtrack, dimid_step
    !output variables
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_product
    integer, dimension(2) :: dimids_xtrack_step
    integer, parameter :: deflate_level = 1
    logical, parameter :: shuffle = .true.

    !define r8 kind for use in setting parameter valid ranges
    integer, parameter :: r8 = kind(1.0d0)
    !character (len=32) :: epoch_buf

    character(len=1000) :: comment1

    if (errstat /= 0) return

    dimids_xtrack_step(1) = dimid_xtrack
    dimids_xtrack_step(2) = dimid_step

    ! Product Fields with optional attributes
    call tiof_attlist_append (att_product, errstat, "coordinates", &
                              att_text = "longitude latitude")

    call tiof_varlist_append (varlist, errstat, &
                              "cloud_pressure", &
                              nf90_int, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "cloud pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.2e3_r8], &
                              fillvalue = fill_int, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

     call tiof_varlist_append (varlist, errstat, &
                              "cloud_fraction", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "effective cloud fraction at 466nm", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

     call tiof_varlist_append (varlist, errstat, &
                              "CloudRadianceFraction466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "cloud radiance fraction at 466nm", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

     call tiof_varlist_append (varlist, errstat, &
                              "CloudRadianceFraction440", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "cloud radiance fraction at 440nm", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_product)

!refer to m_cal_ocp for confirmation
!this string is too long to fit as a comment
    comment1="0: (ERROR) lat/lon/SZA/VZA/RAA error; "// &
            "1: "// &
            "2: (WARNING) pcloud replaced by pscene as 0.<ecf<min_ecf; " // &
            "3: (ERROR) input psfc/rsfc error "// &
            "4: (WARNING) pcloud replaced by pscene as snow_ice_fraction>min_snowice; " // &
            "5: (WARNING) ocp SCD iteration max_iter reached; "// &
            "6: (ERROR) SCD<0. or SCD_MainDataQualityFlag=2(bad); "// &
            "7: (WARNING) 440nm rad or irr error; "// &
            "8: (ERROR) 466nm rad or irr error; " // &
            "9: (ERROR) calculated ecf beyond normal range; "// &
            "10:(WARNING) SceneAlbedoAtTerrain.eq.'yes' skipped or SCD correction problem; " // &
            "11:(WARNING) SceneAlbedoAtTerrain.eq.'no' skpped or SCD correction problem; " // &
            "12:(ERROR) ecf calculation skipped during processing; "// &
            "13:(ERROR) ocp calculation skipped during processing; "// &
            "14:(ERROR) calculated ocp beyond normal range; "// &
            "15:(WARNING) pscene calculation skipped during processing;"

    call tiof_varlist_append (varlist, errstat, &
                              "processing_quality_flag", &
                              nf90_short, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "bitwise processing quality flag", &
                              comment = " ", &
                              valid_range = [0.0_r8, 32767.0_r8], &
                              fillvalue = fill_short, &
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

    call tiof_put2d_i2 (tio_l2obj, "cloud_pressure", [0,0], &
         [nstep, nxtrack], out_CloudPressure, errstat)

    call tiof_put2d_r4 (tio_l2obj, "cloud_fraction", [0,0], &
         [nstep, nxtrack], out_EffectiveCloudFraction, errstat)

    call tiof_put2d_r4 (tio_l2obj, "CloudRadianceFraction466", [0,0], &
         [nstep, nxtrack], out_CloudRadianceFraction466, errstat)

    call tiof_put2d_r4 (tio_l2obj, "CloudRadianceFraction440", [0,0], &
         [nstep, nxtrack], out_CloudRadianceFraction440, errstat)

    write(*,*)'FIXME: processing_quality_flag still in development'
    !!!!!!!!out_ProcessingQualityFlags(:,:) = 0
    call tiof_put2d_i2 (tio_l2obj, "processing_quality_flag", [0,0], &
         [nstep, nxtrack], out_ProcessingQualityFlags, errstat)

    call tiof_pop_group (tio_l2obj, errstat)

    if (errstat /= 0) then
       call tell_error (tell_io_write_error,"write_product_data:failed",errstat)
       return
    endif

   end subroutine write_product_data

  !>Write support structure into L2 netCDF file
  !----------------------------------------------------------------------------
  subroutine write_support_struct(tio_l2obj, dimid_xtrack, dimid_step, errstat)
    use m_vars, only: name_option_SurfaceReflectivity,run_mode
    implicit none

    !input variables
    type (tiof_file_type), intent(inout) :: tio_l2obj
    integer, intent(in) :: dimid_xtrack, dimid_step
    !output variables
    integer, intent(inout) :: errstat
    !local variables
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_support
    integer, dimension(2) :: dimids_xtrack_step
    character(len=255) :: name466, name440
    integer, parameter :: deflate_level = 1
    logical, parameter :: shuffle = .true.

    !define r8 kind for use in setting parameter valid ranges
    integer, parameter :: r8 = kind(1.0d0)
    !character (len=32) :: epoch_buf

    if (errstat /= 0) return

    ! Define dimid arrays associated with common data field shapes.
    dimids_xtrack_step(1) = dimid_xtrack
    dimids_xtrack_step(2) = dimid_step

    ! Product Fields with optional attributes
    call tiof_attlist_append (att_support, errstat, "coordinates", &
                              att_text = "longitude latitude")


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
                              long_name = "466nm surface reflectivity used in calculation", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

     call tiof_varlist_append (varlist, errstat, &
                              name440, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "440nm surface reflectivity used in calculation", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SurfaceLER466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "LER at 466nm calculated at surface pressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SurfaceLER440", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "LER at 440nm calculated at surface pressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SceneLER466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "LER at 466nm calculated at ScenePressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SceneLER440", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "LER at 440nm calculated at ScenePressure", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "ScenePressure", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "scene pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.5e3_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "TerrainPressure", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
          long_name = "terrain pressure calculated for SurfaceLER466", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.5e3_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    if (run_mode .EQ. 'production') then 

    call tiof_varlist_append (varlist, errstat, &
                              "SurfacePressure", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "Pixel Surface Pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.5e3_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    else 
    call tiof_varlist_append (varlist, errstat, &
                              "nonclipped_cloud_fraction", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                    long_name = "nonclipped 466nm effective cloud fraction", &
                              units = "no unit", &
                              valid_range = [-10.0_r8, 10.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SlantColumnAmountO2O2", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "O2-O2 SCD used for cloud pressure", &
                              comment = "O2-O2 slant column at EffectiveTemperature", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SlantColumnSceneO2O2", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "O2-O2 SCD used for scene pressure", &
                              comment = "O2-O2 slant column at O2O2SceneTemperature", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "SlantColumnTerrainO2O2", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "O2-O2 SCD used for surface pressure", &
                              comment = "O2-O2 slant column at O2O2TerrainTemperature", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

!    call tiof_varlist_append (varlist, errstat, &
!                              "ReflectanceFactor466", &
!                              nf90_float, &
!                              dimids = dimids_xtrack_step,  &
!                              long_name = "466nm Reflectance=(Pi*rad466)/(irr466*cos(SZA))", &
!                              units = "no unit", &
!                              valid_range = [0.0_r8, 1.0_r8], &
!                              fillvalue = fill_float_nines, &
!                              deflate_level = deflate_level, &
!                              shuffle = shuffle, &
!                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "rad_of_irr466", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "rad/irr at 466nm for ecf", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "cal_rad_clr", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "cal_rad_clr at 466nm for ecf", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "cal_rad_cld", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "cal_rad_cld at 466nm for ecf", &
                              units = "no unit", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "fitted_slant_column", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                          long_name = "fitted O2-O2 SCD at reference temperature" , &
                              comment = "fitted O2-O2 slant column at 273K", &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)
   
    call tiof_varlist_append (varlist, errstat, &
                              "fitted_slant_column_uncertainty", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                          long_name = "fitted reference O2-O2 SCD uncertainty" , &
                              units = "1.e43 molec^2 cm^-5", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)
    
    call tiof_varlist_append (varlist, errstat, &
                              "fit_rms_residual", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                          long_name = "fitted reference O2-O2 SCD uncertainty" , &
                              units = "unitless", &
                              valid_range = [0.0_r8, 20.0_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)
    
    call tiof_varlist_append (varlist, errstat, &
                              "surface_pressure", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "Pixel Surface Pressure", &
                              units = "hPa", &
                              valid_range = [0.0_r8, 1.5e3_r8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "terrain_height", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              long_name = "terrain height", &
                              units = "m", &
                              valid_range = [-100.0_r8, 10000.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "snow_ice_fraction", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                  long_name = "Fraction of pixel area covered by snow and/or ice", &
                              units = "unitless", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_float, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "O2O2CloudTemperature", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "effective T for SlantColumnAmountO2O2", &
                              units = "K", &
                              valid_range = [160.0_r8, 310.0_r8], &
                              fillvalue = fill_float_nines, &
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
                              fillvalue = fill_float_nines, &
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
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)

    call tiof_varlist_append (varlist, errstat, &
                              "RelativeAzimuthAngle", &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "RAA used in calculation", &
                              units = "degree", &
                              valid_range = [0._r8, 180.0_r8], &
                              fillvalue = fill_float_nines, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_support)
    endif

    call tiof_varlist_append (varlist, errstat, &
                              "SCD_MainDataQualityFlags", &
                              nf90_int, &
                              dimids = dimids_xtrack_step,  &
                   long_name = "main data quality flags for fitted_slant_column", &
                              comment = "0=normal, 1=suspicious, 2=bad", &
                              valid_range = [0.0_r8, 2.0_r8], &
                              fillvalue = -30000.0_r8, &
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
               out_SurfaceLER440, out_SurfaceLER466, out_TerrainHeight,&
               out_SceneLer440, out_SceneLER466, out_ScenePressure,&
               out_SlantColumnSceneO2O2, out_O2O2SceneTemperature,&
               out_SlantColumnTerrainO2O2, out_O2O2TerrainTemperature,&
               rad_SnowIceFraction, nasa_scduncertainty, nasa_scdrms,&
               name_option_SurfaceReflectivity, out_RelativeAzimuthAngle,&
               out_EffectiveCloudFractionNotClipped, &
               out_CloudPressureNotClipped,l2_TerrainPressure, &
               rad_of_irr466, cal_rad_clr, cal_rad_cld

     use m_vars, only: scd_mdqfl, run_mode

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

    if (run_mode .EQ. 'production') then

    call tiof_put2d_r4 (tio_l2obj, "SurfacePressure", [0,0], &
         [nstep, nxtrack], l2_TerrainPressure, errstat)

    else
    call tiof_put2d_r4 (tio_l2obj, "nonclipped_cloud_fraction", [0,0], &
         [nstep, nxtrack], out_EffectiveCloudFractionNotClipped, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SlantColumnAmountO2O2", [0,0], &
         [nstep, nxtrack], out_SlantColumnAmountO2O2, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SlantColumnSceneO2O2", [0,0], &
         [nstep, nxtrack], out_SlantColumnSceneO2O2, errstat)

    call tiof_put2d_r4 (tio_l2obj, "SlantColumnTerrainO2O2", [0,0], &
         [nstep, nxtrack], out_SlantColumnTerrainO2O2, errstat)

!    call tiof_put2d_r4 (tio_l2obj, "ReflectanceFactor466", [0,0], &
!         [nstep, nxtrack], out_ReflectanceFactor, errstat)

    call tiof_put2d_r4 (tio_l2obj, "rad_of_irr466", [0,0], &
         [nstep, nxtrack], rad_of_irr466, errstat)

    call tiof_put2d_r4 (tio_l2obj, "cal_rad_clr", [0,0], &
         [nstep, nxtrack], cal_rad_clr, errstat)

    call tiof_put2d_r4 (tio_l2obj, "cal_rad_cld", [0,0], &
         [nstep, nxtrack], cal_rad_cld, errstat)

    call tiof_put2d_r4 (tio_l2obj, "surface_pressure", [0,0], &
         [nstep, nxtrack], l2_TerrainPressure, errstat)

    call tiof_put2d_r4 (tio_l2obj, "terrain_height", [0,0], &
         [nstep, nxtrack], out_TerrainHeight, errstat)

    call tiof_put2d_r4 (tio_l2obj, "snow_ice_fraction", [0,0], &
         [nstep, nxtrack], rad_SnowIceFraction, errstat)

    call tiof_put2d_r4 (tio_l2obj, "fitted_slant_column", [0,0], &
         [nstep, nxtrack], nasa_SlantColumnAmountO2O2, errstat)

    call tiof_put2d_r4 (tio_l2obj, "fitted_slant_column_uncertainty", [0,0], &
         [nstep, nxtrack], nasa_scduncertainty, errstat)

    call tiof_put2d_r4 (tio_l2obj, "fit_rms_residual", [0,0], &
         [nstep, nxtrack], nasa_scdrms, errstat)

    call tiof_put2d_r4 (tio_l2obj, "O2O2CloudTemperature", [0,0], &
         [nstep, nxtrack], out_O2O2CloudTemperature, errstat)

    call tiof_put2d_r4 (tio_l2obj, "O2O2SceneTemperature", [0,0], &
         [nstep, nxtrack], out_O2O2SceneTemperature, errstat)

    call tiof_put2d_r4 (tio_l2obj, "O2O2TerrainTemperature", [0,0], &
         [nstep, nxtrack], out_O2O2TerrainTemperature, errstat)

    call tiof_put2d_r4 (tio_l2obj, "RelativeAzimuthAngle", [0,0], &
         [nstep, nxtrack], out_RelativeAzimuthAngle, errstat)
    endif

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
