MODULE m_anomflg
  USE PGS_PC_class
  USE OMI_SMF_class
  USE HE4_class
  USE OMI_LUN_set
  IMPLICIT NONE
  CHARACTER (LEN=256) :: anomflg_fn
  INTEGER(KIND = 4) :: ntimes,nxtrack,nlats
  PUBLIC :: get_anomflg_3
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: anomflg_3
CONTAINS

  FUNCTION GET_ANOMFLG_3(year, doy, anomflg_3) RESULT(status)
    INTEGER(KIND=4) :: status, ierr
    INTEGER(KIND=4) :: version
    CHARACTER(LEN=PGS_SMF_MAX_MSG_SIZE) :: msg,stryear,strdoy
    INTEGER(KIND=4), PARAMETER :: zero = 0
    CHARACTER(LEN=128) :: filename
    INTEGER(KIND=4) :: i
    INTEGER(KIND=4), INTENT(IN) :: year, doy
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: anomflg_3
    CHARACTER(LEN=20) :: FUNCTIONNAME = "get_anomflg_3"
    INTEGER(KIND = 4) ::  fid, swid

    INTEGER(KIND = 4) :: start1(1), stride1(1), edge1(1)
    INTEGER(KIND = 4) :: start3(3), stride3(3), edge3(3)
    INTEGER(KIND = 4), DIMENSION(:,:,:), ALLOCATABLE :: flgs
    INTEGER(KIND = 4), DIMENSION(:), ALLOCATABLE :: flg_year, flg_doy

    version = 1
    status = PGS_PC_getReference(ANOMFLG3_LUN,version, anomflg_fn)

    IF( status /= PGS_S_SUCCESS ) THEN
      WRITE( msg,'(A)' ) " Get anomaly flag file name failed "
      ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, "GET_ANOMFLG_3", zero )
      CALL EXIT(1)
    ENDIF

    fid = swopen(anomflg_fn, DFACC_READ)
    IF(fid == -1 ) THEN
      WRITE( msg,'(A)' ) " Open anomaly flag file failed "
      ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, "GET_ANOMFLG_3", zero )
      CALL EXIT(1)
    ENDIF

    swid = swattach(fid, "OMTO3 Row Anomaly Data" )

    ntimes = swdiminfo(swid, "Time")
    nxtrack = swdiminfo(swid, "Xtrack")
    nlats = swdiminfo(swid, "Lats")

    start1(1)=0
    stride1(1)=1
    edge1(1)=ntimes
    start3(1)=0
    start3(2)=0
    start3(3)=0
    stride3(1)=1
    stride3(2)=1
    stride3(3)=1
    edge3(1)=ntimes
    edge3(2)=nxtrack
    edge3(3)=nlats

    ALLOCATE(flgs(ntimes,nxtrack,nlats),anomflg_3(nxtrack,nlats), &
         flg_year(ntimes),flg_doy(ntimes))

    status=swrdfld(swid,"Year",start1,stride1,edge1,flg_year)
    status=swrdfld(swid,"DayOfYear",start1,stride1,edge1,flg_doy)
    status=swrdfld(swid,"AnomalyFlags",start3,stride3,edge3,flgs)

    status = swdetach(swid)
    status = swclose(fid)

    WRITE(stryear,'(I4)') year
    WRITE(strdoy,'(I4)') doy

    status=0
    DO i=1,ntimes
      IF (flg_year(i) == year .AND. flg_doy(i) == doy) THEN
        anomflg_3(:,:)=flgs(i,:,:)
        WRITE(msg,'(A)' ) "Anomaly flags retrieved for "//TRIM(stryear)//"/"//TRIM(strdoy)
        ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, "GET_ANOMFLG_3", zero )
        status=1
      ENDIF
    ENDDO

    IF (status /= 1) THEN
      WRITE( msg,'(A)' ) "No anomaly flags found for "//TRIM(stryear)//"/"//TRIM(strdoy)//" in"//anomflg_fn
      ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, "GET_ANOMFLG_3", zero )
      CALL EXIT(1)
    ENDIF

    !write(34,'(60I1)') anomflg_3

    DEALLOCATE(flgs,flg_year,flg_doy)
  END FUNCTION GET_ANOMFLG_3


END MODULE m_anomflg
