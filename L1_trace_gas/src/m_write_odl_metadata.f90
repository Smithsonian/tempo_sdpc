!>Function to write ODL-format metadata to an ASCII file
module m_write_odl_metadata

  use netcdf
  use tell_module
  use tio_module
  use tg_names_module
  use ISO_C_BINDING, only: C_NULL_CHAR, C_DOUBLE, C_INT, C_CHAR
  use OMSAO_indices_module, only: mcf_lun
  use OMSAO_omidata_module, only: omi_radiance_swathname
  use OMSAO_variables_module, only: mdlist_filename
  use m_pgs_include
  use md_module

  implicit none

  private
  public write_odl_metadata

contains

  !>Write ODL-format metadata to an ASCII .met output file
  !-----------------------------------------------------------------------
  !
  !> @param[in] l1bfile    L1 radiance netCDF file name
  !> @param[in] outfilnm   L2 product netCDF file name
  !> @param[in] nXtrack    number of cross-track pixels
  !> @param[in] nLines     number of mirror steps
  !> @param[in] lun_input  list of LUNs for input files
  !> @param[in] ninp       number of input files to be listed (n_lun_input)
  !> @param     errstat    error tracking code, non-zero indicates problem
  !
  ! @author E. O'Sullivan  March 2018
  !
  ! Note: as yet, this is only a proof-of-concept test to ensure we can
  !       produce ASCII format ODL metadata when processing netCDF only.
  !-----------------------------------------------------------------------
  function write_odl_metadata(l1bfile, outfilnm, nXtrack, nLines, &
       lun_input, ninp) result(errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l1bfile, outfilnm
    integer (kind=4), intent(in) :: nXtrack, nLines
    integer (kind=4), intent(in) :: ninp
    integer (kind=4), intent(in), dimension(ninp) :: lun_input
    integer :: errstat

    !local variables
    integer, parameter :: INVENTORY=2
    integer, parameter :: ninvname=5


    real (kind=4), dimension(nXtrack, nLines) :: lon, lat

    integer :: pgs_MET_setAttr_s, pgs_MET_setAttr_i, pgs_MET_setmultiAttr_s, &
         pgs_MET_setmultiAttr_d, pgs_MET_setmultiAttr_i
    integer :: pgs_met_init,pgs_met_write, pgs_pc_getreference, &
         pgs_met_sfstart, pgs_met_sfend, pgs_met_remove

    integer :: i, status, returnstatus, version, sdid, Fil_Lun, j

    character(LEN=PGSd_MET_GROUP_NAME_L), dimension(PGSd_MET_NUM_OF_GROUPS) :: GROUPS
    character(LEN=100), dimension(50) :: Objvalue
    character(LEN=100), dimension(ninp) :: InputPnt, supflnm
    character(LEN=200) :: buf, localgranuleid
    character(LEN=PGSd_MET_GROUP_NAME_L),dimension(ninvname), parameter :: &
         INVOBJ = (/                                     &
         "INPUTPOINTER                     ", &
         "RANGEENDINGDATE                  ", &
         "RANGEENDINGTIME                  ", &
         "RANGEBEGINNINGDATE               ", &
         "RANGEBEGINNINGTIME               "/)
    character (kind=C_CHAR) :: NULL = C_NULL_CHAR



    ! Additional attributes
    integer, parameter :: nadd = 2
    character (len=32), dimension(nadd) :: AddAttrNam, AddAttrVal

    ! bounding polgon / footprint parameters
    integer, parameter :: max_npts=100
    integer :: npts
    integer (kind=C_INT), dimension(max_npts) :: polygon_seq
    real (kind=C_DOUBLE), dimension(max_npts) :: polygon_lats, polygon_lons
    real (kind=4), dimension(max_npts) :: p_lats, p_lons
    logical, dimension(nXtrack, nLines) :: valid
    real (kind=4) :: center_lat, center_lon, clon, clat

    integer :: ncerr
    character (len=32) :: cov_start_string, cov_end_string

    type (tiof_file_type) :: tio_l1obj

    status = OMI_S_SUCCESS
    errstat = 0


    !
    ! TBD - here we need a section reading metadata from input radiance file
    !


    ! TBD - read metadata values from PCF file if you want to operate on them

    ! Get the start and end date and time from the input file
    ! TBD - probably want to move this to m_read_metadata_tio eventually
    call tiof_open(l1bfile, tio_l1obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "write_odl_metadata: failed to open L1 radiance file", errstat)
      return
    endif

    ncerr = nf90_get_att (tio_l1obj % fileid, nf90_global, &
         "time_coverage_start", cov_start_string)
    if (ncerr /= nf90_noerr) then
      call tell_error (tell_io_read_error, &
           "write_odl_metadata: failed to read time_coverage_start", &
           errstat)
      return
    endif

    ncerr = nf90_get_att (tio_l1obj % fileid, nf90_global, &
         "time_coverage_end", cov_end_string)
    if (ncerr /= nf90_noerr) then
      call tell_error (tell_io_read_error, &
           "write_odl_metadata: failed to read time_coverage_end", &
           errstat)
      return
    endif

    ! We're going to need to read lon and lat in YET AGAIN because the
    ! retrieval doesn't retain those arrays
    call tiof_push_group (tio_l1obj, omi_radiance_swathname, errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_latitude, [0,0], &
         [nLines, nXtrack], lat(1:nXtrack,1:nLines), errstat)
    call tiof_get2d_r4 (tio_l1obj, tg_var_longitude, [0,0], &
         [nLines, nXtrack], lon(1:nXtrack,1:nLines), errstat)
    call tiof_pop_group (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "write_odl_metadata: failed to read longitude, latitude arrays", &
           errstat)
      return
    endif


    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "write_odl_metadata: failed to close L1 radiance file", errstat)
      return
    endif

    read(cov_end_string,'(a10,1x,a8)') Objvalue(2), Objvalue(3)
    read(cov_start_string,'(a10,1x,a8)') Objvalue(4), Objvalue(5)

    ! Bounding polgon points
    ! For now, do the simplest thing, a bounding box
    call bounding_box_md(nXtrack, nLines, lat, lon, p_lats, &
         p_lons, npts, errstat)
    if (npts > max_npts) then
      call tell_error (tell_io_write_error, &
           "write_geo_bounds_md: npts in polygon exceeds max allowed", &
           errstat)
      return
    endif
    polygon_lats=p_lats
    polygon_lons=p_lons
    do i=1,npts
      polygon_seq(i) = i
    enddo

    ! Mean longitude and latitude
    where (lat.ge.-90.0d0 .and. lat.le.90.0d0 .and. &
         lon.ge.-180.0d0 .and. lon.le.180.0d0)
      valid=.true.
    elsewhere
      valid=.false.
    end where
    center_lon=sum(lon,mask=valid)/count(valid)
    center_lat=sum(lat,mask=valid)/count(valid)

    ! Centroid values classed as additional attributes
    AddAttrNam(1) = 'CENTROID_MEAN_LONGITUDE'
    AddAttrNam(2) = 'CENTROID_MEAN_LATITUDE'
    write(AddAttrVal(1),'(f7.1)') center_lon
    write(AddAttrVal(2),'(f7.1)') center_lat

    ! FIXME - at present the code only includes a very limited set of input
    ! files (RAD, IRRAD, RADREF, PREFITS). It should really include all the
    ! reference datasets used

    ! Input files
    do i=1,ninp
      Fil_Lun=lun_input(i)
      version = 1

      returnstatus = PGS_PC_GetReference( Fil_Lun, version, buf )
      if( returnstatus /= 0 ) then
        call tell_error(tell_io_read_error, &
             "write_odl_metadata: failed to read input file names", errstat)
        return
      else
        j = index( buf, '/', BACK = .true. ) + 1
        supflnm(i) = trim( buf( j:) )
      endif
    enddo

    InputPnt=supflnm

    ! Set up metadata file in memory
    returnstatus = pgs_met_init(mcf_lun, GROUPS)
    if (returnstatus /= 0 ) then
      call tell_error(tell_io_error, &
           "write_odl_metadata: failed to initialise metadata file", errstat)
      return
    endif

    j = index( outfilnm, '/', BACK = .true. ) + 1
    localgranuleid = trim(outfilnm(j:))
    returnstatus = pgs_MET_setattr_s(GROUPS(INVENTORY), "LOCALGRANULEID", &
         trim(localgranuleid))
    if (returnstatus /= 0 ) then
      call tell_error(tell_io_error, &
           "write_odl_metadata: failed to set localgranuleid", errstat)
      return
    endif


    do i=1,ninvname
      returnstatus = pgs_met_setattr_s(GROUPS(INVENTORY),trim(INVOBJ(i)), &
           Objvalue(i))
      if (returnstatus /= 0 ) then
        call tell_error(tell_io_error, &
             "write_odl_metadata: failed to set time/date ranges", errstat)
        return
      endif
    enddo

    returnstatus = pgs_MET_setmultiAttr_s(GROUPS(INVENTORY),"InputPointer", &
         ninp,InputPnt)

    if (returnstatus /= 0) then
      call tell_error(tell_io_error, &
           "write_odl_metadata: failed to set input pointer", errstat)
      return
    endif

    returnstatus = pgs_MET_setmultiAttr_i(GROUPS(INVENTORY), &
         "GRINGPOINTSEQUENCENO.1", npts, polygon_seq(1:npts))
    returnstatus = pgs_MET_setmultiAttr_d(GROUPS(INVENTORY), &
         "GRINGPOINTLATITUDE.1", npts, polygon_lats(1:npts))
    returnstatus = pgs_MET_setmultiAttr_d(GROUPS(INVENTORY), &
         "GRINGPOINTLONGITUDE.1", npts, polygon_lons(1:npts))

    if(returnstatus /= 0)then
      call tell_error(tell_io_error, &
           "write_odl_metadata: failed to set bounding polygon", errstat)
      return
    endif

    do i=1,nadd
      write(buf,*) i
      buf=adjustl(buf)
      returnstatus = pgs_met_setAttr_i(GROUPS(INVENTORY), &
           "ADDITIONALATTRIBUTENAME."//trim(buf),AddAttrNam(i))
      returnstatus = pgs_met_setAttr_i(GROUPS(INVENTORY), &
           "PARAMETERVALUE."//trim(buf),AddAttrVal(i))
      if (returnstatus /= 0) then
        call tell_error(tell_io_error, &
             "write_odl_metadata: failed to set additional attr "//trim(buf), &
             errstat)
        return
      endif
    enddo

    ! Write archive metadata attributes to netCDF file
    ! do this first since pgs_met functions apparently leave nc file open!
    call open_md(outfilnm, errstat)
    call write_geo_bounds_md(nXtrack, nLines, lat, lon, errstat)
    call write_inputs_md(ninp, InputPnt, errstat)
    call write_fixed_md(mdlist_filename,errstat)!md_namelist, errstat)
    call write_prodid_md(outfilnm,"(1)",errstat)
    call close_md(errstat)

    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "write_odl_metadata: failed writing netCDF attributes", &
           errstat)
      return
    endif

    ! write ODL-format text file
    version =1

    returnstatus = pgs_met_sfstart( trim(outfilnm), HDF5_ACC_RDWR,sdid)
    if(returnstatus /= 0) then
      call tell_error(tell_io_error, &
           "write_odl_metadata: failed to open .met file", errstat)
      return
    endif

    ! Write metadata file
    returnstatus = pgs_met_write(groups(INVENTORY), NULL, sdid)
    ! PGSMET_E_SD_SETATTR indicates HDFEOS fie not set, as we intend
    if(returnstatus /= PGSMET_E_SD_SETATTR .AND. returnstatus /= 0) then
      call tell_error (tell_io_write_error, &
           "write_odl_metadata: failed to write ArchiveMetadata", errstat)
      return
    endif

    returnstatus = pgs_met_sfend(sdid)

    if(returnstatus /=0) then
      call tell_error(tell_io_error, &
           "write_odl_metadata: failed to close .met file", errstat)
      return
    endif

    returnstatus = pgs_met_remove()




    call tell_log (1,"write_odl_metadata: success")

    return
  end function write_odl_metadata

end module m_write_odl_metadata
