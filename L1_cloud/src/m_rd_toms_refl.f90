!>Read TOMS surface reflectivity climatology reference file
!
!-------------------------------------------------------------------------
!
! !ROUTINE: rd_toms_refl
!
! !DESCRIPTION: get TOMS reflectivity climatology
!
! !CALLING SEQUENCE: call rd_toms_refl (lat, lon, terr_pres)
!
! !INPUT PARAMETERS: 
!> @param  lat[in]      latitude
!> @param  lon[in]      longitude
!
! !OUTPUT PARAMETERS:  
!> @param  terr_pres[out]     terrain pressure
!
! !SEE ALSO: 
!
! !REVISION HISTORY:
!
!> @author  01Aug07   Joiner      Original code
!> @author  26Mar15   O'Sullivan  Update for TEMPO
!-------------------------------------------------------------------------
module m_rd_toms_refl

  private
  public rd_toms_refl

contains

  subroutine rd_toms_refl (errstat)   

    use m_vars, ONLY: done_read_refl, lat, lon, toms_refl, ref_nmon, &
         iLine, nXtrack, ref_nlat, ref_nlon, ref_clr, ref_lats, ref_lons, &
         month, ler_sz, ler_th, ler_ph, ler354, startlat, startlon, deltlat, &
         deltlon
    use m_LUN_set
    use m_pgs_include
    use tell_module
    implicit NONE          

    integer, intent(inout) :: errstat

    !local variables
    !================
    integer :: lun=2 
    integer :: pgs_io_gen_openf, pgs_io_gen_closef
    integer :: status, version=1
    integer :: ipts, i, j
    real (KIND=8) :: lont, latt
    character(len=100) :: txt, logmsg


    if (errstat /= 0) return

    !=======================
    !read terrain data set
    !=======================
    if (.not. done_read_refl) then
      status = pgs_io_gen_openf ( refl_id, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
      if(status.ne.0) then
        call tell_error(tell_io_open_error, &
             "rd_toms_refl: error opening reflectivity file", errstat)
        return
      endif
      write(logmsg,"(A49, I4)") &
           'rd_toms_refl: opening reflectivity file, status :',status
      call tell_log(1,logmsg)
      read (lun,*,err=200)  txt
      call tell_log(1,txt)
      read (lun,*,err=200)  txt
      call tell_log(1,txt)
      read (lun,*,err=200)  ref_nlon, ref_nlat, ref_nmon
      write(logmsg,"(3I6)") ref_nlon, ref_nlat, ref_nmon
      call tell_log(1,logmsg)

      !Allocate memory, read in data
      allocate(ref_lats(ref_nlat), ref_lons(ref_nlon), &
           toms_refl(ref_nmon,ref_nlon,ref_nlat), stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "rd_toms_refl: failed to allocate memory", &
             errstat)
        return
      endif
      read (lun,*,err=200)  ref_lons
      read (lun,*,err=200)  ref_lats

!      if (iprt > 1) print *, ref_lons  
!      if (iprt > 1) print *, ref_lats
      read (lun,*,err=200)  toms_refl
      status = pgs_io_gen_closef (lun)
      write(logmsg,"(A49, I4)") &
           'rd_toms_refl: closing reflectivity file, status :',status
      call tell_log(1,logmsg)

      status = pgs_io_gen_openf ( ler354_id, PGSd_IO_Gen_RSeqFrm, &
           0,lun, version)
      if(status.ne.0) then
        call tell_error(tell_io_open_error, &
             "rd_toms_refl: error opening ler354_cox_munk file", errstat)
        return
      endif

      read (lun,*,err=201) ler_sz
      read (lun,*,err=201) ler_th
      read (lun,*,err=201) ler_ph
      read (lun,*,err=201) ler354
      status = pgs_io_gen_closef (lun)
      write(logmsg,"(A52, I4)") &
           'rd_toms_refl: closing ler354_cox_munk file, status :',status
      call tell_log(1,logmsg)
      done_read_refl=.true.
      deltlat=ref_nlat/180.
      deltlon=ref_nlon/360.
      startlat=ref_lats(1)
      startlon=ref_lons(1)
    endif ! not done reading

    !if (iprt > 0) print *,'deltlat, lon ',deltlat, deltlon

    do ipts=1, nXtrack
      ! Check lat/lon within bounds, ref_clr=0 if not
      if (lat(ipts,iLine) > 90.0 .or. lat(ipts,iLine) < -90.0 .or. &
           lon(ipts,iLine) > 180.0 .or. lon(ipts,iLine) < -180.0) then
        ref_clr(ipts-1,iLine)=0.0d0
      else
        latt=lat(ipts,iLine)
        i=int(anint((startlat-latt)/deltlat),kind(i))+1
        if (i <= 0) i=1
        if (i >= ref_nlat+1) i=ref_nlat
        lont=lon(ipts,iLine)   
        if(lont > 180.) then 
          lont=lont-360 
        endif
        j=int(anint((lont-startlon)/deltlon),kind(j))+1
        if(j <= 0) j=1   
        if(j >= ref_nlon+1) j=1   
        ref_clr(ipts-1,iLine)=toms_refl(month,j,i)/100.
      endif
    enddo   ! ipts

    return

    !Note - OMCLDRR aborts at these errors but could continue
200 call tell_error(tell_io_read_error, &
         "rd_toms_refl: error reading reflectivity file", errstat)
    toms_refl=0.d0

201 call tell_error(tell_io_read_error, &
    'rd_toms_refl: error reading Cox-Munk LER file', errstat)
    ler354=0.d0

  end subroutine rd_toms_refl

end module m_rd_toms_refl
