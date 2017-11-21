!>Function to write ODL-format metadata to an ASCII file
module m_write_metadata

  use m_LUN_set
  use m_pgs_include
  use m_vars, only: filename_in_nc
  use tell_module
  use netcdf
  use tio_module
  use ISO_C_BINDING, only: C_NULL_CHAR

  implicit none
  private

  public write_metadata

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
  !-----------------------------------------------------------------------
  function write_metadata(outfilnm) result(errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: outfilnm
    integer :: errstat

    !local variables
    integer, parameter :: INVENTORY=2
    integer, parameter :: ninp = 9
    integer, parameter :: ninvname=5

    integer :: pgs_MET_setAttr_s, pgs_MET_setmultiAttr_s
    integer :: pgs_met_init,pgs_met_write, pgs_pc_getreference, &
         pgs_met_sfstart, pgs_met_sfend, pgs_met_remove

    integer :: i, status, returnstatus, version, sdid, Fil_Lun, j, &
         resid_id

    character(LEN=PGSd_MET_GROUP_NAME_L) :: GROUPS(PGSd_MET_NUM_OF_GROUPS)
    character(LEN=100), dimension(50) :: Objvalue
    character(LEN=100), dimension(ninp) :: InputPnt,supflnm
    character(LEN=200) :: buf
    character(LEN=PGSd_MET_GROUP_NAME_L),dimension(ninvname), parameter :: &
         INVOBJ = (/                                     &
         "INPUTPOINTER                     ", &
         "RANGEENDINGDATE                  ", &
         "RANGEENDINGTIME                  ", &
         "RANGEBEGINNINGDATE               ", &
         "RANGEBEGINNINGTIME               "/)
    character (len=4) :: NULL = C_NULL_CHAR

    integer :: ncerr
    character (len=256) :: cov_start_string, cov_end_string
    type(tiof_file_type) :: tio_l1obj

    status = OMI_S_SUCCESS

    !
    ! TBD - here we need a section reading metadata from input radiance file
    !


    ! TBD - read metadata values from PCF file if you want to operate on them

    ! Get the start and end date and time from the input file
    ! TBD - probably want to move this to m_read_metadata_tio eventually
    call tiof_open(filename_in_nc, tio_l1obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "write_metadata: failed to open L1 radiance file", errstat)
      return
    endif

    ncerr = nf90_get_att (tio_l1obj % fileid, nf90_global, &
         "time_coverage_start", cov_start_string)
    if (ncerr /= nf90_noerr) then
      call tell_error (tell_io_read_error, &
           "write_metadata: failed to read time_coverage_start", &
           errstat)
      return
    endif

    ncerr = nf90_get_att (tio_l1obj % fileid, nf90_global, &
         "time_coverage_end", cov_end_string)
    if (ncerr /= nf90_noerr) then
      call tell_error (tell_io_read_error, &
           "write_metadata: failed to read time_coverage_end", &
           errstat)
      return
    endif

    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "write_metadata: failed to close L1 radiance file", errstat)
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
             "write_metadata: failed to read input file names", errstat)
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
    if(returnstatus /= 0 ) then
      call tell_error(tell_io_error, &
           "write_metadata: failed to initialise metadata file", errstat)
      return
    endif


    do i=1,ninvname
      returnstatus = pgs_met_setattr_s(GROUPS(INVENTORY),trim(INVOBJ(i)), &
           Objvalue(i))
      if(returnstatus /= 0 ) then
        call tell_error(tell_io_error, &
             "write_metadata: failed to set input filenames", errstat)
        return
      endif
    enddo

    returnstatus = pgs_MET_setmultiAttr_s(GROUPS(INVENTORY),"InputPointer", &
         ninp,InputPnt)

    if(returnstatus /=0)then
      call tell_error(tell_io_error, &
           "write_metadata: failed to set input pointer", errstat)
      return
    endif

    version =1

    returnstatus = pgs_met_sfstart( trim(outfilnm), HDF5_ACC_RDWR,sdid)

    if(returnstatus /=0) then
      call tell_error(tell_io_error, &
           "write_metadata: failed to open .met file", errstat)
      return
    endif

    ! Write metadata file
    returnstatus = pgs_met_write(groups(INVENTORY), NULL, sdid)
    ! PGSMET_E_SD_SETATTR indicates HDFEOS fie not set, as we intend
    if(returnstatus /= PGSMET_E_SD_SETATTR .AND. returnstatus /= 0) then
      call tell_error (tell_io_write_error, &
           "write_metadata: failed to write ArchiveMetadata", errstat)
      return
    endif

    returnstatus = pgs_met_sfend(sdid)

    if(returnstatus /=0) then
      call tell_error(tell_io_error, &
           "write_metadata: failed to close .met file", errstat)
      return
    endif

    returnstatus = pgs_met_remove()



    call tell_log (1,"write_metadata: success")

    return
  end function write_metadata
end module m_write_metadata
