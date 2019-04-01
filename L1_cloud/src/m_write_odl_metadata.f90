!>Function to write ODL-format metadata to an ASCII file
module m_write_odl_metadata

  use m_LUN_set
  use m_pgs_include
  use m_vars, only: filename_in_nc, nc_swathname
  use tell_module
  use netcdf, only: nf90_nowrite, nf90_global, nf90_noerr, nf90_get_att
  use tio_module
  use ISO_C_BINDING, only: C_NULL_CHAR, C_DOUBLE, C_INT, C_CHAR

  implicit none
  private

  type, public :: boundary_polygon_type
    real (kind=4), dimension(:), allocatable :: lons, lats
    integer (kind=4), dimension(:), allocatable :: seq
    real (kind=4) :: centroid_lat, centroid_lon
  end type

  public write_odl_metadata, make_bounding_polygon

contains

  subroutine make_bounding_polygon (bdry, errstat)
    implicit none

    ! input variables
    type (boundary_polygon_type), intent(inout) :: bdry
    integer, intent(inout) :: errstat

    ! local variables
    type (tiof_file_type) :: tio_l1obj
    integer :: npts, alloc_status, i

    if (errstat /= 0) return

    call tiof_open (filename_in_nc, tio_l1obj, nf90_nowrite, errstat)

    ! Bounding polygon and centroid
    call tiof_push_group (tio_l1obj, trim(nc_swathname), errstat)
    call tiof_make_lev1_bounding_polygon (tio_l1obj, bdry % lons, bdry % lats, &
                                          bdry % centroid_lon, bdry % centroid_lat, errstat)
    call tiof_pop_group (tio_l1obj, errstat)
    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
                       "write_odl_metadata: failed generating bounding polygon", &
                       errstat)
      return
    endif

    npts = size(bdry % lons)
    allocate (bdry % seq(npts), stat=alloc_status)
    if (alloc_status /= 0) then
      call tell_error (tell_malloc_error, "malloc failed", errstat)
      return
    endif
    bdry % seq(1:npts) = (/(i, i=1,npts)/)

  end subroutine

  !>Write ODL-format metadata to an ASCII .met output file
  !-----------------------------------------------------------------------
  !
  !> @param[in] outfilnm   Filename of output L2 netCDF product file
  !> @param     errstat    error tracking code, non-zero indicates problem
  !
  ! @author E. O'Sullivan  November 2017
  ! @author J. Houck, March 2019
  !
  ! Note: as yet, this is only a proof-of-concept test to ensure we can
  !       produce ASCII format ODL metadata when processing netCDF only.
  ! Also: note that there's a write_metadata subroutine in
  !       m_write_output_data_tio that writes to the netCDF file
  !-----------------------------------------------------------------------
  function write_odl_metadata (outfilnm, bdry) result(errstat)

    implicit none

    !input variables
    type(boundary_polygon_type), intent(in) :: bdry
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

    integer :: ncerr, npts
    character (len=32) :: cov_start_string, cov_end_string

    type(tiof_file_type) :: tio_l1obj

    status = OMI_S_SUCCESS
    errstat = 0

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

    ! Centroid values classed as additional attributes
    AddAttrNam(1) = 'CENTROID_MEAN_LONGITUDE'
    AddAttrNam(2) = 'CENTROID_MEAN_LATITUDE'
    write(AddAttrVal(1),'(f10.5)') bdry % centroid_lon
    write(AddAttrVal(2),'(f10.5)') bdry % centroid_lat

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

    npts = size(bdry % lons)

    returnstatus = pgs_MET_setmultiAttr_i(GROUPS(INVENTORY), &
         "GRINGPOINTSEQUENCENO.1", npts, bdry % seq(1:npts))
    returnstatus = pgs_MET_setmultiAttr_d(GROUPS(INVENTORY), &
         "GRINGPOINTLATITUDE.1", npts, real(bdry % lats(1:npts), kind=8) )
    returnstatus = pgs_MET_setmultiAttr_d(GROUPS(INVENTORY), &
         "GRINGPOINTLONGITUDE.1", npts, real(bdry % lons(1:npts), kind=8) )

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
