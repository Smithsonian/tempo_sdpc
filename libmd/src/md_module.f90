!>Write archive metadata as attributes in level 2 netCDF files
module md_module

  use netcdf, only: nf90_write, nf90_global
  use tell_module
  use tio_module

  implicit none

  private
  public md_open, md_close, md_write_geo_bounds, md_write_inputs, &
       md_write_prodid, md_write_fixed

  type (tiof_file_type), private, target :: primary_output_file

contains

  !>Open level 2 netCDF file for metadata writing
  !-----------------------------------------------------------------------
  !
  !> @param[in]    l2file  netCDF output filename
  !> @param        errstat error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan  March 2018
  !
  !-----------------------------------------------------------------------
  subroutine md_open (l2file, errstat)

    implicit none

    character (len=*), intent(in) :: l2file
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: l2obj

    if (errstat /= 0) return

    l2obj => primary_output_file

    call tiof_open(l2file, l2obj, nf90_write, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_open_error, &
                       "md_open: opening file "//trim(l2file), &
                       errstat)
      return
    endif

  end subroutine md_open

  !>Close level 2 netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param errstat error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan   March 2018
  !-----------------------------------------------------------------------
  subroutine md_close (errstat)
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: l2obj

    l2obj => primary_output_file

    call tiof_close (l2obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, "md_close failed", errstat)
    endif

  end subroutine md_close

  !>Write geographic bounding polygon and centroid values
  !-----------------------------------------------------------------------
  !
  !> @param[in]  lon           Bounding polygon longitude coordinates
  !> @param[in]  lat           Bounding polygon latitude coordinates
  !> @param[in]  centroid_lon  Bounding polygon centroid longitude
  !> @param[in]  centroid_lat  Bounding polygon centroid latitude
  !> @param[in/out]  errstat error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan  March 2018
  !> @author J. Houck,  March 2019
  !-----------------------------------------------------------------------
  subroutine md_write_geo_bounds (lon, lat, centroid_lon, centroid_lat, errstat)

    implicit none

    !input variables
    real (kind=4), dimension(:), intent(in) :: lon, lat
    real (kind=4), intent(in) :: centroid_lat, centroid_lon
    integer (kind=4), intent (inout) :: errstat

    !local variables
    integer (kind=4), allocatable, dimension(:) :: seq
    type (tiof_file_type), pointer :: l2obj
    type (tiof_attlist_type) :: attlist
    integer status, i, num

    if (errstat /= 0) return

    l2obj => primary_output_file

    num = size(lon)
    allocate (seq(num), stat=status)
    if (status /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif

    seq(1:num) = (/(i, i=1,num)/)

    !write to l2obj
    call tiof_attlist_append (attlist, errstat, "centroid_mean_latitude", &
         att_r4=[centroid_lat])
    call tiof_attlist_append (attlist, errstat, "centroid_mean_longitude", &
         att_r4=[centroid_lon])
    call tiof_attlist_append (attlist, errstat, "polygon_latitudes", &
         att_r4=[lat(1:num)])
    call tiof_attlist_append (attlist, errstat, "polygon_longitudes", &
         att_r4=[lon(1:num)])
    call tiof_attlist_append (attlist, errstat, "polygon_sequence", &
         att_i4=[seq(1:num)])

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "md_write_geo_bounds: failed", &
           errstat)
      return
    endif
    call tiof_push_group (l2obj, "metadata", errstat)
    call tiof_def_atts (l2obj, attlist, nf90_global, errstat)
    call tiof_pop_group (l2obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "md_write_geo_bounds: failed", &
           errstat)
      return
    endif

  end subroutine md_write_geo_bounds

  !> Write a list of input files into metadata
  !-----------------------------------------------------------------------
  !
  !> @param[in]    ninp    number of input files
  !> @param[in]    inputs  filenames as a 1D array of strings, length ninp
  !> @param        errstat error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan  March 2018
  !-----------------------------------------------------------------------
  subroutine md_write_inputs (ninp, inputs, errstat)

    implicit none

    integer (kind=4), intent(in) :: ninp
    character (len=*), dimension(ninp), intent(in) :: inputs
    integer (kind=4), intent(inout) :: errstat

    character (len=128*ninp) :: input_list
    integer :: n
    type (tiof_attlist_type) :: attlist
    type (tiof_file_type), pointer :: l2obj

    if (errstat /= 0) return

    l2obj => primary_output_file

    input_list=''
    input_list = trim(adjustl(inputs(1)))
    do n=2,ninp
      input_list = trim(input_list)//', '//trim(adjustl(inputs(n)))
    enddo

    call tiof_attlist_append (attlist, errstat, "input_files", &
         att_text=trim(input_list))

    call tiof_push_group (l2obj, "metadata", errstat)
    call tiof_def_atts (l2obj, attlist, nf90_global, errstat)
    call tiof_pop_group (l2obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "md_write_inputs: failed", &
           errstat)
      return
    endif

  end subroutine md_write_inputs

  !> Write product id, version, production date/time to netCDF attributes
  !-----------------------------------------------------------------------
  !
  !> @param[in]    filename product file name (used as unique identifier)
  !> @param[in]    version product version number
  !> @param        errstat  error tracking code, non-zero indicates error
  !
  !> @author E. O'Sullivan  March 2018
  !-----------------------------------------------------------------------
  subroutine md_write_prodid (filename, version, errstat)

    implicit none

    character (len=*), intent(in) :: filename
    integer (kind=4), intent(in) :: version
    integer (kind=4), intent(inout) :: errstat

    type (tiof_attlist_type) :: attlist
    type (tiof_file_type), pointer :: l2obj
    character (len=8) :: date
    character (len=10) :: time
    character (len=5) :: zone
    character (len=32) :: prod_datetime
    integer :: j

    if (errstat /= 0) return

    l2obj => primary_output_file

    !FIXME Eventually we'll need a way to put current date & time in
    ! TEMPO time standard, but for now just use system time
    call date_and_time(DATE=date, TIME=time, ZONE=zone)
    prod_datetime=date(1:4)//'-'//date(5:6)//'-'//date(7:8) &
         //'T'//time(1:2)//':'//time(3:4)//':'//time(5:10)//' UTC'//zone

    !trim any leading path from the product filename
    j = index( filename, '/', BACK = .true. ) + 1

    call tiof_attlist_append (attlist, errstat, "local_granule_id", &
         att_text=trim(adjustl(filename(j:))))
    call tiof_attlist_append (attlist, errstat, "version_id", &
         att_i4=[version])
    call tiof_attlist_append (attlist, errstat, "production_date_time", &
         att_text=prod_datetime)

    call tiof_push_group (l2obj, "metadata", errstat)
    call tiof_def_atts (l2obj, attlist, nf90_global, errstat)
    call tiof_pop_group (l2obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "md_write_prodid: failed", &
           errstat)
      return
    endif

  end subroutine md_write_prodid

  !>Write fixed product-specific metadata from namelist into netCDF
  !----------------------------------------------------------------------
  !
  !> @param[in]    nlfile    filename of namelist data file to read
  !> @param        errstat   error tracking code, non-zero indicates error
  !
  !> @author E. O'Sullivan  March 2018
  !----------------------------------------------------------------------
  subroutine md_write_fixed (nlfile, errstat)

    implicit none

    ! input variables
    character (len=*), intent(in) :: nlfile
    integer (kind=4), intent(inout) :: errstat

    !local variables read from namelist
    character (len=32) :: collection_shortname, platform
    character (len=5) :: collection_version
    character (len=2048) :: abstract
    character (len=128) :: access_description
    character (len=32) :: access_value
    character (len=1024) :: keywords

    character (len=5), parameter :: inst='TEMPO'

    type (tiof_attlist_type) :: attlist
    type (tiof_file_type), pointer :: l2obj

    namelist /metadata/  collection_shortname, collection_version, platform, &
         access_description, access_value, abstract, keywords

    if (errstat /= 0) return

    l2obj => primary_output_file

    ! read from namelist
    open (888, file=nlfile, status='OLD', iostat=errstat)
    if (errstat == 0) read (888, nml=metadata, iostat=errstat)
    if (errstat == 0) close (888, iostat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "md_write_fixed: failed to read from namelist", errstat)
      return
    endif

    call tiof_attlist_append (attlist, errstat, "collection_shortname", &
         att_text=trim(adjustl(collection_shortname)))
    call tiof_attlist_append (attlist, errstat, "collection_version", &
         att_text=trim(adjustl(collection_version)))
    call tiof_attlist_append (attlist, errstat, "platform", &
         att_text=trim(adjustl(platform)))
    call tiof_attlist_append (attlist, errstat, "instrument", &
         att_text=inst)
    call tiof_attlist_append (attlist, errstat, "access_description", &
         att_text=trim(adjustl(access_description)))
    call tiof_attlist_append (attlist, errstat, "access_value", &
         att_text=trim(adjustl(access_value)))
    call tiof_attlist_append (attlist, errstat, "abstract", &
         att_text=trim(adjustl(abstract)))
    call tiof_attlist_append (attlist, errstat, "keywords", &
         att_text=trim(adjustl(keywords)))

    call tiof_push_group (l2obj, "metadata", errstat)
    call tiof_def_atts (l2obj, attlist, nf90_global, errstat)
    call tiof_pop_group (l2obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "md_write_fixed: failed", &
           errstat)
      return
    endif

  end subroutine md_write_fixed

end module md_module
