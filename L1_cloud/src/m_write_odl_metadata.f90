!>Function to write ODL-format metadata to an ASCII file
module m_write_odl_metadata

  use m_LUN_set
  use m_pgs_include
  use m_vars, only: filename_in_nc, lon, lat, nXtrack, nLines
  use tell_module
  use netcdf, only: nf90_nowrite, nf90_global, nf90_noerr, nf90_get_att
  use tio_module
  use ISO_C_BINDING, only: C_NULL_CHAR, C_DOUBLE, C_INT, C_CHAR

  implicit none
  private

  public write_odl_metadata

contains

  !>Write ODL-format metadata to an ASCII .met output file
  !-----------------------------------------------------------------------
  !
  !> @param[in] outfilnm   Filename of output L2 netCDF product file
  !> @param     errstat    error tracking code, non-zero indicates problem
  !
  ! @author E. O'Sullivan  November 2017
  !
  ! Note: as yet, this is only a proof-of-concept test to ensure we can
  !       produce ASCII format ODL metadata when processing netCDF only.
  ! Also: note that there's a write_metadata subroutine in
  !       m_write_output_data_tio that writes to the netCDF file
  !-----------------------------------------------------------------------
  function write_odl_metadata(outfilnm) result(errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: outfilnm
    integer :: errstat

    !local variables
    integer, parameter :: INVENTORY=2
    integer, parameter :: ninp = 9
    integer, parameter :: ninvname=5

    integer :: pgs_MET_setAttr_s, pgs_MET_setAttr_i, pgs_MET_setmultiAttr_s, &
         pgs_MET_setmultiAttr_d, pgs_MET_setmultiAttr_i
    integer :: pgs_met_init,pgs_met_write, pgs_pc_getreference, &
         pgs_met_sfstart, pgs_met_sfend, pgs_met_remove

    integer :: i, status, returnstatus, version, sdid, Fil_Lun, j, &
         resid_id

    character(LEN=PGSd_MET_GROUP_NAME_L), dimension(PGSd_MET_NUM_OF_GROUPS) :: GROUPS
    character(LEN=100), dimension(50) :: Objvalue
    character(LEN=100), dimension(ninp) :: InputPnt,supflnm
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
    integer, parameter :: nintermed = 3 ! change dimensions and MCF if changing from 3
    integer :: npts, n, xtidx, atidx
    integer (kind=4), dimension(0:15) :: xtstep, atstep
    integer (kind=C_INT), dimension(16) :: polygon_seq
    real (kind=C_DOUBLE), dimension(16) :: polygon_lats, polygon_lons
    real (kind=4) :: center_lat, center_lon

    integer :: ncerr
    character (len=32) :: cov_start_string, cov_end_string

    type(tiof_file_type) :: tio_l1obj

    status = OMI_S_SUCCESS
    errstat = 0

    ! Bounding polgon points
    ! One position for each corner, nintermed points along each side
    ! Note that MCF specifies total number of values in polygon, so has to be
    ! updated if you change nintermed, as well as polygon variable dimensions
    npts = 4+(4*nintermed) ! number of points in polygon

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
      xtidx=1+int((nXtrack-1)*xtstep(n-1)/(nintermed+1))
      atidx=1+int((nLines-1)*atstep(n-1)/(nintermed+1))
      polygon_lats(n)=lat(xtidx, atidx)
      polygon_lons(n)=lon(xtidx, atidx)
      polygon_seq(n)=n
    enddo
    ! Mean longitude and latitude
    center_lon=lon(nXtrack/2,nLines/2)
    center_lat=lat(nXtrack/2,nLines/2)

    ! Centroid values classed as additional attributes
    AddAttrNam(1) = 'CENTROID_MEAN_LONGITUDE'
    AddAttrNam(2) = 'CENTROID_MEAN_LATITUDE'
    write(AddAttrVal(1),'(f14.9)') center_lon
    write(AddAttrVal(2),'(f14.9)') center_lat

    !
    ! TBD - here we need a section reading metadata from input radiance file
    !


    ! TBD - read metadata values from PCF file if you want to operate on them

    ! Get the start and end date and time from the input file
    ! TBD - probably want to move this to m_read_metadata_tio eventually
    call tiof_open(filename_in_nc, tio_l1obj, nf90_nowrite, errstat)
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

    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "write_odl_metadata: failed to close L1 radiance file", errstat)
      return
    endif

    read(cov_end_string,'(a10,1x,a8)') Objvalue(2), Objvalue(3)
    read(cov_start_string,'(a10,1x,a8)') Objvalue(4), Objvalue(5)

    ! Input files
    resid_id = resid_id_late
    do i=1,ninp
      if(i==1)Fil_Lun=L1B_LUN
      if(i==2)Fil_Lun=IRR1B_file
      if(i==3)Fil_Lun=terr_prs_id
      if(i==4)Fil_Lun=chl_id
      if(i==5)Fil_Lun=oc_ram_id
      if(i==6)Fil_Lun=ring_id
      if(i==7)Fil_Lun=thresh_id
      if(i==8)Fil_Lun=resid_id
      if(i==9)Fil_Lun=refl_id
      if(i==10) then
        supflnm(i)=''
        go to 50
      endif
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
50    continue
    enddo

    InputPnt=supflnm

    ! Set up metadata file in memory
    returnstatus = pgs_met_init(MCF_LUN, GROUPS)
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
         "GRINGPOINTSEQUENCENO.1", npts, polygon_seq)
    returnstatus = pgs_MET_setmultiAttr_d(GROUPS(INVENTORY), &
         "GRINGPOINTLATITUDE.1", npts, polygon_lats)
    returnstatus = pgs_MET_setmultiAttr_d(GROUPS(INVENTORY), &
         "GRINGPOINTLONGITUDE.1", npts, polygon_lons)

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
