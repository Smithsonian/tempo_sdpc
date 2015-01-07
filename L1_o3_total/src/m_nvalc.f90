MODULE m_nvalc
  !!-----------------------------------------------------------------------
  !!This module read the orbital n-value correction ( or adjustment ) in he4
  !!and output the corrections(or adjustment) for this orbit 
  !!
  !! programmer : Shifang Luo  SSAI                       10/25/2004
  !!
  !!-----------------------------------------------------------------------
  USE PGS_PC_class
  USE OMI_SMF_class
  USE HE4_class

  IMPLICIT NONE

  INTEGER(KIND = 4), PARAMETER :: nvCORR_LUN = 420001
  CHARACTER (LEN=256) :: nvcorr_flnm
  INTEGER(KIND = 4) :: nWvc, nXtc, nOrbc
  INTEGER(KIND = 4) ::  fid, swid


  PUBLIC :: OmiNVCinfo
  PUBLIC :: OmiNvalueCorr
!  PUBLIC :: writeNVXcorr

CONTAINS

  FUNCTION OmiNVCinfo()  RESULT(status)
    !!------------------------------------------------------------------
    !! Get the dimension information from orbital N-Value correction
    !! file in he4
    !!-----------------------------------------------------------------

    CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
    INTEGER (KIND = 4), PARAMETER :: zero = 0
    INTEGER (KIND = 4) :: status, ierr, version

    !
    ! start
    !  
    version = 1
    status = PGS_PC_getReference(nvCORR_LUN,version, nvcorr_flnm)

    IF( status /= PGS_S_SUCCESS ) THEN
      WRITE( msg,'(A)' ) " Get N-value correction file name failed "
      ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, "", zero )
      CALL EXIT(1)
    ENDIF

    fid = swopen(nvcorr_flnm, DFACC_READ)
    IF(fid == -1 ) THEN
      WRITE( msg,'(A)' ) " Open N-value correction file name failed "
      ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, "", zero )
      CALL EXIT(1)
    ENDIF

    swid = swattach(fid, "OrbitNvalueCorrection" )
    nWvc = swdiminfo(swid, "nWavelc")
    nXtc = swdiminfo(swid, "nXtrack")
    nOrbc = swdiminfo(swid, "nOrb")

    status = swdetach(swid)
    status = swclose(fid)

  END FUNCTION OmiNVCinfo

  FUNCTION OmiNvalueCorr(orbit, crwl,swpcr) RESULT(status) 
    !
    !***********************************************************************
    ! OmiNvalueCorr was from omi_6ch_nvcorv1
    !
    !     omi_6ch_nvcorv1         june 2004  s. taylor SSAI
    !
    !     purpose
    !        to provide OMI Nvalue corrections at six EP  wavelengths
    !        for all 60 scan positions by
    !        interpolating from an external table for a given orbit     
    !
    !     method
    !        for the first call, read the external table of corrections.
    !        Interpolate the table in time to provide swpcr(:,:), the
    !        Nvalue corrections for each of the x wavelengths
    !        crwl(x) and scan positions for orbit.
    !
    !     variables
    !       name        type i/o  description
    !       ----        ---- ---  -----------
    !
    !      argument
    !       orbit        i*4  i   orbit number of requested corrections
    !       crwl(:)      r*4  o   wavelengths
    !       swpcr(:,:)   r*4  o   correction for given year and doy
    !
    !     calling routine       
    ! 
    !  Original Programmer: Steve Taylor.
    !  October, 2004:   Shifang Luo --   SSAI  
    !          get N-Value adjustments from he4 file
    !
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    USE UTIL_tools_class
    !
    CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
    INTEGER (KIND = 4), PARAMETER :: zero = 0
    INTEGER(KIND=4), INTENT(IN) :: orbit
    REAL(KIND = 4), DIMENSION(nWvc,nXtc),INTENT(OUT) :: swpcr
    REAL(KIND = 4), DIMENSION(nWvc), INTENT(OUT) :: crwl
    REAL(KIND = 4), DIMENSION(nWvc,nXtc,nOrbc) :: scr 
    REAL(KIND = 4), DIMENSION(nOrbc) :: norb
    REAL(KIND = 4) :: frac
    INTEGER(KIND = 4) :: start1(1), stride1(1), edge1(1)
    INTEGER(KIND = 4) :: start2(1), stride2(1), edge2(1)
    INTEGER(KIND = 4) :: start3(3), stride3(3), edge3(3)
    INTEGER(KIND=4) :: i,ji, iw, isw, is, j
    INTEGER(KIND = 4) :: status, ierr, version

    data is/0/
    !
    version = 1
    status = PGS_PC_getReference(nvCORR_LUN,version, nvcorr_flnm)

    IF( status /= PGS_S_SUCCESS ) THEN   
      WRITE( msg,'(A)' ) " Get N-value correction file name failed "
      ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, "", zero )
      CALL EXIT(1)
    ENDIF

    fid = swopen(nvcorr_flnm, DFACC_READ) 
    IF(fid == -1 ) THEN
      WRITE( msg,'(A)' ) " Open N-value correction file name failed "
      ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, "", zero )
      CALL EXIT(1)
    ENDIF

    swid = swattach(fid, "OrbitNvalueCorrection" )

    start1(1) = 0
    stride1(1) = 1
    edge1(1) = nWvc

    start2(1) = 0
    stride2(1) = 1
    edge2(1) = nOrbc

    start3(1) = 0
    start3(2) = 0
    start3(3) = 0
    stride3(1) = 1
    stride3(2) = 1
    stride3(3) = 1
    edge3(1) = nWvc
    edge3(2) = nXtc
    edge3(3) = nOrbc

    status = swrdfld(swid, "WavelengthOfCorrection",start1,stride1,edge1,&
         crwl)
    status = swrdfld(swid,"OrbitNumber", start2,stride2,edge2,&
         norb)
    status = swrdfld(swid, "NvaluesCorrection", start3,stride3,edge3,&
         scr)

    is = hunt(norb, REAL(orbit), is)  
    !!
    !! interpolate correction by orbit
    IF( is > 1) THEN
      frac = (REAL(orbit)-norb(is))/(norb(is+1)-norb(is))
      do iw = 1,nWvc
        do isw=1,nXtc
          swpcr(iw,isw) = frac*(scr(iw,isw,is+1)-scr(iw,isw,is)) + &
               scr(iw,isw,is) 
        enddo
      enddo
    ELSE
      swpcr = scr(:,:,1) 
    ENDIF
    !!
    IF( is == 0 ) THEN
      swpcr =scr(:,:, nOrbc) 
      RETURN 
    ENDIF
    !!
    status = OZT_S_SUCCESS

  END FUNCTION OmiNvalueCorr

!  FUNCTION writeNVXcorr(L2_filename,crwl, swpcr) RESULT(status)
!    !!-----------------------------------------------------------------------
!    !! Write out the N-value correction for the currect orbit to
!    !! Global Attribute in L2 output
!    !!-----------------------------------------------------------------------
!    USE HE5_class
!    CHARACTER( LEN = *), INTENT(IN) :: L2_filename
!    REAL(KIND = 4), DIMENSION(nWvc,nXtc), INTENT(IN) :: swpcr
!    REAL(KIND = 4), DIMENSION(nWvc), INTENT(IN) :: crwl
!    INTEGER (KIND=4) :: nc, numtype, SW_fileid
!    INTEGER (KIND=4) :: ierr, status
!
!    SW_fileid = he5_swopen( L2_filename, HE5F_ACC_RDWR )
!    nc =SIZE(crwl) 
!    numtype = HE5T_NATIVE_FLOAT
!    ierr = he5_ehwrglatt( SW_fileid, "WavelengthOfAdjustment", &
!         numtype, nc, crwl )
!    IF( ierr == -1 ) THEN
!      status  = - 1
!    ENDIF
!
!    nc = SIZE(swpcr)
!    numtype = HE5T_NATIVE_FLOAT
!    ierr = he5_ehwrglatt( SW_fileid, "NVXAdjustment", &
!         numtype, nc, swpcr )
!    IF( ierr == -1 ) THEN
!      status  = - 1
!    ELSE
!      status = 0
!    ENDIF
!
!  END FUNCTION writeNVXcorr


END MODULE m_nvalc

