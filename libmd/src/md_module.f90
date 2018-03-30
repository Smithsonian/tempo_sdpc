!>Write archive metadata as attributes in level 2 netCDF files
module md_module

  use netcdf, only: nf90_write, nf90_global
  use tell_module
  use tio_module

  implicit none

  private
  public open_md, close_md, write_geo_bounds_md, write_inputs_md, &
       write_prodid_md, write_fixed_md

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
  subroutine open_md (l2file, errstat)

    implicit none

    character (len=*), intent(in) :: l2file
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: l2obj

    if (errstat /= 0) return

    l2obj => primary_output_file

    call tiof_open(l2file, l2obj, nf90_write, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_open_error, &
                       "open_md: opening file "//trim(l2file), &
                       errstat)
      return
    endif

  end subroutine open_md


  !>Close level 2 netCDF file
  !-----------------------------------------------------------------------
  !
  !> @param errstat error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan   March 2018
  !-----------------------------------------------------------------------
  subroutine close_md (errstat)
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: l2obj

    l2obj => primary_output_file

    call tiof_close (l2obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, "close_md failed", errstat)
    endif

  end subroutine close_md


  !>Write geographic bounding polygon and centroid values
  !-----------------------------------------------------------------------
  !
  !> @param[in]    nxtrack cross-track dimension size (e.g., 2000)
  !> @param[in]    nstep   mirror step dimension size (e.g., 128)
  !> @param[in]    lat     latitude array (nxtrack, nstep)
  !> @param[in]    lon     longitude array
  !> @param        errstat error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan  March 2018
  !
  !-----------------------------------------------------------------------
  subroutine write_geo_bounds_md (nxtrack, nstep, lat, lon, errstat)

    implicit none

    !input variables
    integer (kind=4), intent(in) :: nxtrack, nstep
    real (kind=4), dimension(nxtrack, nstep), intent(in) :: lat, lon
    integer (kind=4), intent (inout) :: errstat

    !local variables
    integer, parameter :: nintermed = 3
    integer :: npts, n, xtidx, atidx
    integer (kind=4), dimension(0:4*nintermed+3) :: xtstep, atstep
    real (kind=4), dimension(4*nintermed+4) :: polygon_lats, polygon_lons
    integer (kind=4), dimension(4*nintermed+4) :: polygon_seq
    real (kind=4) :: center_lat, center_lon
    type (tiof_attlist_type) :: attlist
    type (tiof_file_type), pointer :: l2obj

    if (errstat /= 0) return

    l2obj => primary_output_file

    !number of points in polygon
    npts = 4+(4*nintermed)

    !polygon points
    do n=0,npts/2
      xtstep(n)=0+(n-(nintermed+1))
      atstep(n)=0+n
      if(atstep(n) > (nintermed+1)) atstep(n) = nintermed+1
      if(xtstep(n) < 0) xtstep(n) = 0
    enddo
    do n=(npts/2)+1,npts-1
      xtstep(n)=atstep(npts-n)
      atstep(n)=xtstep(npts-n)
    enddo
    do n=1,npts
      xtidx=1+int((nxtrack-1)*xtstep(n-1)/(nintermed+1))
      atidx=1+int((nstep-1)*atstep(n-1)/(nintermed+1))
      polygon_lats(n)=lat(xtidx, atidx)
      polygon_lons(n)=lon(xtidx, atidx)
      polygon_seq(n)=n
    enddo

    ! Centroid longitude and latitude
    center_lon=lon(nxtrack/2,nstep/2)
    center_lat=lat(nxtrack/2,nstep/2)

    !write to l2obj
    call tiof_attlist_append (attlist, errstat, "centroid_mean_latitude", &
         att_r4=[center_lat])
    call tiof_attlist_append (attlist, errstat, "centroid_mean_longitude", &
         att_r4=[center_lon])
    call tiof_attlist_append (attlist, errstat, "polygon_latitudes", &
         att_r4=[polygon_lats])
    call tiof_attlist_append (attlist, errstat, "polygon_longitudes", &
         att_r4=[polygon_lons])
    call tiof_attlist_append (attlist, errstat, "polygon_sequence", &
         att_i4=[polygon_seq])

    call tiof_push_group (l2obj, "metadata", errstat)
    call tiof_def_atts (l2obj, attlist, nf90_global, errstat)
    call tiof_pop_group (l2obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_geo_bounds_md: failed", &
           errstat)
      return
    endif

  end subroutine write_geo_bounds_md


  !> Write a list of input files into metadata
  !-----------------------------------------------------------------------
  !
  !> @param[in]    ninp    number of input files
  !> @param[in]    inputs  filenames as a 1D array of strings, length ninp
  !> @param        errstat error tracking code, non-zero indicates problem
  !
  !> @author E. O'Sullivan  March 2018
  !-----------------------------------------------------------------------
  subroutine write_inputs_md (ninp, inputs, errstat)

    implicit none

    integer (kind=4), intent(in) :: ninp
    character (len=128), dimension(ninp), intent(in) :: inputs
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
         att_text=input_list)

    call tiof_push_group (l2obj, "metadata", errstat)
    call tiof_def_atts (l2obj, attlist, nf90_global, errstat)
    call tiof_pop_group (l2obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_inputs_md: failed", &
           errstat)
      return
    endif

  end subroutine write_inputs_md


  !> Write product id, version, production date/time to netCDF attributes
  !-----------------------------------------------------------------------
  !
  !> @param[in]    filename product file name (used as unique identifier)
  !> @param[in]    version_str  product version as a string
  !> @param        errstat  error tracking code, non-zero indicates error
  !
  !> @author E. O'Sullivan  March 2018
  !-----------------------------------------------------------------------
  subroutine write_prodid_md (filename, version_str, errstat)

    implicit none

    character (len=*), intent(in) :: filename, version_str
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
    ! FIXME - not clear how local_version_id should be defined, may well
    ! be a checksum or commit hash. For now just using a placeholder string
    call tiof_attlist_append (attlist, errstat, "local_version_id", &
         att_text=trim(adjustl(version_str)))
    call tiof_attlist_append (attlist, errstat, "production_date_time", &
         att_text=prod_datetime)

    call tiof_push_group (l2obj, "metadata", errstat)
    call tiof_def_atts (l2obj, attlist, nf90_global, errstat)
    call tiof_pop_group (l2obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "write_prodid_md: failed", &
           errstat)
      return
    endif

  end subroutine write_prodid_md


  !>Write fixed product-specific metadata from namelist into netCDF
  !----------------------------------------------------------------------
  !
  !> @param[in]    nlfile    filename of namelist data file to read
  !> @param        errstat   error tracking code, non-zero indicates error
  !
  !> @author E. O'Sullivan  March 2018
  !----------------------------------------------------------------------
  subroutine write_fixed_md (nlfile, errstat)

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
           "write_fixed_md: failed to read from namelist", errstat)
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
      call tell_error (tell_io_write_error, "write_fixed_md: failed", &
           errstat)
      return
    endif

  end subroutine write_fixed_md



end module md_module
