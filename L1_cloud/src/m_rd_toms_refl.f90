module m_rd_toms_refl

  private
  public rd_toms_refl

contains

  subroutine rd_toms_refl (errstat)   

    use m_vars, ONLY: done_read_refl, iprt, lat, lon, toms_refl, ref_nmon, &
         iLine, nXtrack, ref_nlat, ref_nlon, ref_clr, ref_lats, ref_lons, &
         month, ler_sz, ler_th, ler_ph, ler354, startlat, startlon, deltlat, &
         deltlon
    use m_LUN_set
    use m_pgs_include
    use tell_module
    implicit NONE          

    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !!BOP
    !
    ! !ROUTINE: rd_toms_refl
    !
    ! !DESCRIPTION: get TOMS reflectivity climatology
    !
    ! !CALLING SEQUENCE: call rd_toms_refl (lat, lon, terr_pres)
    !
    ! !INPUT PARAMETERS: 
    !real (KIND=8), dimension(:), intent(in)  :: lat, lon
    !                     lat     : latitude
    !                     lon     : longitude
    !
    ! !OUTPUT PARAMETERS:  
    !real (KIND=8), dimension(:), intent(out) :: terr_pres
    !                     terr_pres      : terrain pressure
    !
    ! !SEE ALSO: 
    !
    ! !REVISION HISTORY:
    !
    !  01Aug07   Joiner      Original code
    !  26Mar15   O'Sullivan  Update for TEMPO
    !!EOP
    !-------------------------------------------------------------------------
    integer, intent(inout) :: errstat

    !local variables
    !================
    integer :: lun=2 
    integer :: pgs_io_gen_openf, pgs_io_gen_closef, OMI_SMF_setmsg
    integer :: status,ierr, version=1
    integer :: ipts, i, j
    real (KIND=8) :: lont, latt
    character(len=100) :: txt


    if (errstat /= 0) return

    !=======================
    !read terrain data set
    !=======================
    if (.not. done_read_refl) then
      status = pgs_io_gen_openf ( refl_id, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
      if(status.ne.0) then
        ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN, &
             'error opening reflectivity file', &
             'rd_toms_refl, module m_rd_toms_refl',2)
        errstat = -1
        return
      endif
      if (iprt > 0) print *, &
           'rd_toms_refl: opening reflectivity file, status :',status
      read (lun,*,err=200)  txt
      if (iprt > 0) print *, txt
      read (lun,*,err=200)  txt
      if (iprt > 0) print *, txt
      read (lun,*,err=200)  ref_nlon, ref_nlat, ref_nmon
      if (iprt > 0) print *, ref_nlon, ref_nlat, ref_nmon

      !Allocate memory, read in data
      allocate(ref_lats(ref_nlat), ref_lons(ref_nlon), &
           toms_refl(ref_nmon,ref_nlon,ref_nlat), stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "rd_toms_refl: failed to allocate memory", &
             errstat)
        errstat = -1
        return
      endif
      read (lun,*,err=200)  ref_lons
      read (lun,*,err=200)  ref_lats

      if (iprt > 1) print *, ref_lons  
      if (iprt > 1) print *, ref_lats
      read (lun,*,err=200)  toms_refl
      status = pgs_io_gen_closef (lun)
      if (iprt > 0) print *, &
           'rd_toms_refl: closing reflectivity file, status :',status

      status = pgs_io_gen_openf ( ler354_id, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
      if(status.ne.0) then
        ierr=OMI_SMF_setmsg(OMI_E_FILE_OPEN,'error opening ler354_cox_munk file', &
             'rd_toms_refl, module m_rd_toms_refl',2)
        errstat = -1
        return
      endif

      read (lun,*,err=201) ler_sz
      read (lun,*,err=201) ler_th
      read (lun,*,err=201) ler_ph
      read (lun,*,err=201) ler354
      status = pgs_io_gen_closef (lun)
      if (iprt > 0) print *, &
           'rd_toms_refl: closing ler354_cox_munk file, status :',status
      done_read_refl=.true.
      deltlat=ref_nlat/180.
      deltlon=ref_nlon/360.
      startlat=ref_lats(1)
      startlon=ref_lons(1)
    endif ! not done reading

    !if (iprt > 0) print *,'deltlat, lon ',deltlat, deltlon

    do ipts=1, nXtrack
      latt=lat(ipts,iLine)
      i=anint((startlat-latt)/deltlat, kind=4)+1   
      if (i <= 0) i=1
      if (i >= ref_nlat+1) i=ref_nlat
      lont=lon(ipts,iLine)   
      if(lont > 180.) then 
        lont=lont-360 
      endif
      j=anint((lont-startlon)/deltlon, kind=4)+1   
      if(j <= 0) j=1   
      if(j >= ref_nlon+1) j=1   
      ! print *, 'rd_toms_refl : ',ipts, lat(ipts,iLine), lon(ipts,iLine), &
      !   i, j
      ! print *, toms_refl(month,j,i)
      ref_clr(ipts-1,iLine)=toms_refl(month,j,i)/100.
    enddo   ! ipts

    return
200 print *, 'rd_toms_refl: error reading reflectivity file'
    toms_refl=0.
201 print *, 'rd_toms_refl: error reading Cox-Munk LER file'
    ler354=0.

  end subroutine rd_toms_refl

end module m_rd_toms_refl
