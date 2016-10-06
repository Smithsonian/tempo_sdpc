!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_QA_class
! 
! This module contains a set of functions that set the QualityFlags for
! each output pixel in the OMTO3 product, set the MeasurementQualityFlag
! for each scan line, check to make sure the instrument configuration ID
! are compatible between input solar irradiance and Earth view radiance
! products, gather zonal minimun and maximun values, and gather statics.
! for each individaul granules.
!
!!Revision History:
! Initial version 03/26/2002  Kai Yang/UMBC
!
!!Team-unique Header:
! This software was developed by the OMI Science Team Support
! Group for the National Aeronautics and Space Administration, Goddard
! Space Flight Center, under NASA Task 916-003-1
!
!!References and Credits
! Written by 
! Kai Yang 
! University of Maryland Baltimore County
! email: Kai.Yang-1@nasa.gov
! 
!!Design Notes
!
!!END
!!****************************************************************************
MODULE O3T_QA_class
  USE L1B_class
  USE PGS_MET_class
  IMPLICIT NONE
  PUBLIC  :: O3T_setMqaL2
  PUBLIC  :: O3T_setQAflgsI
  PUBLIC  :: O3T_setQAflgsP
  PUBLIC  :: O3T_instIDsCheck
  PUBLIC  :: O3T_nMissing
  PUBLIC  :: O3T_zonalMinMax
  PUBLIC  :: O3T_QAcounter

CONTAINS
  FUNCTION O3T_setMqaL2( mqaL1B, radLMissing, instIDmismatch, &
       bit7Q, L2param )  RESULT( mqaL2 )

    implicit none

    INTEGER (KIND=2), INTENT(IN) :: mqaL1B
    TYPE (L2PARAM_T), INTENT(INOUT) :: L2param
    LOGICAL, INTENT(IN) :: radLMissing, instIDmismatch, bit7Q
    INTEGER (KIND=1) :: mqaL2
    INTEGER (KIND = 1) :: mesL_error, mesL_warning
    INTEGER (KIND = 1) :: rebinF, saaF, maneuverF
    !INTEGER (KIND = 1) :: mqaL_error

    mesL_error = L1B_extractMqa( mqaL1B, mesL_warning, &
         rebinF, saaF, maneuverF  )
    mqaL2 = 0
    IF( radLMissing      ) THEN
      mqaL2 = IBSET( mqaL2, 0  )
      L2param%NumberOfMeasurementMissing    = &
           L2param%NumberOfMeasurementMissing + 1
    ENDIF

    IF( mesL_error > 0   ) THEN
      mqaL2 = IBSET( mqaL2, 1  )
      L2param%NumberOfMeasurementError      = &
           L2param%NumberOfMeasurementError + 1
    ENDIF

    IF( mesL_warning > 0 ) THEN
      mqaL2 = IBSET( mqaL2, 2  )
      L2param%NumberOfMeasurementWarning    = &
           L2param%NumberOfMeasurementWarning + 1
    ENDIF
    IF( rebinF > 0       ) THEN
      mqaL2 = IBSET( mqaL2, 3  )
      L2param%NumberOfMeasurementRebinned   = &
           L2param%NumberOfMeasurementRebinned + 1
    ENDIF
    IF( saaF > 0         ) THEN
      mqaL2 = IBSET( mqaL2, 4  )
      L2param%NumberOfMeasurementSAA        = &
           L2param%NumberOfMeasurementSAA  + 1
    ENDIF
    IF( maneuverF > 0    ) THEN
      mqaL2 = IBSET( mqaL2, 5  )
      L2param%NumberOfMeasurementManeuver   = &
           L2param%NumberOfMeasurementManeuver + 1
    ENDIF
    IF( instIDmismatch   ) THEN
      mqaL2 = IBSET( mqaL2, 6  )
      L2param%NumberOfInstrumentSettingsError= &
           L2param%NumberOfInstrumentSettingsError + 1
    ENDIF
    IF( bit7Q            ) THEN
      mqaL2 = IBSET( mqaL2, 7  )
    ENDIF

    RETURN
  END FUNCTION O3T_setMqaL2




  FUNCTION O3T_instIDsCheck( ShortName_rad, mqa_rad, instID_rad, &
       instID_irr ) RESULT( instIDmismatch )

    implicit none

    CHARACTER(LEN=*), INTENT(IN) :: ShortName_rad
    INTEGER (KIND=2), INTENT(IN) :: mqa_rad
    INTEGER (KIND=1), INTENT(IN) :: instID_rad, instID_irr
    INTEGER(KIND=1) :: rebinF
    INTEGER(KIND=4) :: l
    LOGICAL :: instIDmismatch

    rebinF = int(IBITS( mqa_rad, 7, 1 ),kind(rebinF))
    l = LEN_TRIM( ShortName_rad )

    IF( ShortName_rad(l:l) == 'G' ) THEN
      IF( instID_irr /= 8 ) THEN
        instIDmismatch = .TRUE.
        RETURN
      ENDIF

      IF( rebinF == 0 ) THEN
        IF( instID_rad == 0 .OR. instID_rad == 1 .OR. &
             instID_rad == 2 .OR. instID_rad == 7 ) THEN
          instIDmismatch = .FALSE.
          RETURN
        ELSE
          instIDmismatch = .TRUE.
          RETURN
        ENDIF
      ELSE IF( rebinF == 1 ) THEN
        IF( (42 <= instID_rad .AND. instID_rad <= 44) .OR. &
             instID_rad == 49 .OR. &
             ( 56 <= instID_rad .AND. instID_rad <= 58 ) ) THEN
          instIDmismatch = .FALSE.
          RETURN
        ELSE
          instIDmismatch = .TRUE.
          RETURN
        ENDIF
      ELSE
        instIDmismatch = .TRUE.
        RETURN
      ENDIF
    ELSE IF( ShortName_rad(l:l) == 'Z' ) THEN
      IF( rebinF /= 0 ) THEN
        instIDmismatch = .TRUE.
        RETURN
      ENDIF

      IF( instID_irr == 50 ) THEN
        IF( (42 <= instID_rad .AND. instID_rad <= 44) .OR. &
             instID_rad == 49 ) THEN
          instIDmismatch = .FALSE.
          RETURN
        ELSE
          instIDmismatch = .TRUE.
          RETURN
        ENDIF
      ELSE IF( instID_irr == 62 ) THEN
        IF( 56 <= instID_rad .AND. instID_rad <= 58 ) THEN
          instIDmismatch = .FALSE.
          RETURN
        ELSE
          instIDmismatch = .TRUE.
          RETURN
        ENDIF
      ELSE
        instIDmismatch = .TRUE.
        RETURN
      ENDIF
    ELSE
      instIDmismatch = .TRUE.
      RETURN
    ENDIF

    instIDmismatch = .TRUE.
    RETURN
  END FUNCTION O3T_instIDsCheck

  FUNCTION O3T_nMissing( iwlArray, earthRadQFLine ) RESULT( Nmissing )

    implicit none

    INTEGER (KIND=4), DIMENSION(:), INTENT(IN) :: iwlArray
    INTEGER (KIND=4) :: nwl, nXtrack
    INTEGER (KIND=2), DIMENSION(:,:), INTENT(IN) :: earthRadQFLine
    INTEGER (KIND=4) :: iwl, iX, pixMissing, nMissing

    nwl = SIZE( iwlArray )
    nXtrack = SIZE( earthRadQFLine, 2 )
    nMissing = 0
    DO iX = 1, nXtrack
      pixMissing = 0
      DO iwl = 1, nwl
        pixMissing = pixMissing &
             + IBITS( earthRadQFLine(iwlArray(iwl), iX), 0, 1 )
      ENDDO
      IF( pixMissing > 0 ) nMissing = nMissing + 1
    ENDDO
    RETURN
  END FUNCTION O3T_nMissing




  FUNCTION O3T_setQAflgsI( QAflags, radBadPixflgs, anomflg_3, &
       algflg, L2param, descendQ, PclimQ, sza, gflg, iwlArray, &
       earthRadQF, solarIrrQF) &
       RESULT( skipit )
!  FUNCTION O3T_setQAflgsI( QAflags, radBadPixflgs, anomflg, anomflg_3, &
!       algflg, L2param, &
!       descendQ, PclimQ, sza, gflg, iwlArray, &
!       earthRadQF, solarIrrQF, earthRadPN, &
!       solarIrrPN, earthWlPN,  solarWlPN )  &
!       RESULT( skipit )

    implicit none

    INTEGER (KIND=2), INTENT(INOUT) :: QAflags, radBadPixflgs
    INTEGER (KIND=1), INTENT(INOUT) :: algflg
    !INTEGER (KIND=1), INTENT(IN) :: anomflg ! L1B flag
    !INTEGER (KIND=1) :: mask7 = 7 ! mask for bits 0-2
    INTEGER (KIND=4), INTENT(IN) :: anomflg_3 ! O3 derived flag
    TYPE (L2PARAM_T), INTENT(INOUT) :: L2param
    LOGICAL, INTENT(IN) :: descendQ, PclimQ
    REAL (KIND=4), INTENT(IN)    :: sza
    INTEGER (KIND=2), INTENT(IN) :: gflg
    INTEGER (KIND=4), DIMENSION(:), INTENT(IN) :: iwlArray
    INTEGER (KIND=4) :: nwl
    INTEGER (KIND=2), DIMENSION(:), INTENT(IN) :: earthRadQF, solarIrrQF
    !REAL (KIND=4), DIMENSION(:), INTENT(IN)    :: earthRadPN, solarIrrPN
    !REAL (KIND=4), DIMENSION(:), INTENT(IN), OPTIONAL :: earthWlPN,  &
    !     solarWlPN
    LOGICAL :: skipit, badpixQ

    INTEGER (KIND=4) :: radMissing, radBadPix, radError, radWarning
    INTEGER (KIND=4) :: irrMissing, irrError, irrWarning

    INTEGER (KIND=4) :: iwl

    nwl = SIZE( iwlArray )
    radMissing = 0
    radBadPix  = 0
    radError   = 0
    radWarning = 0
    DO iwl = 1, nwl
      radMissing = radMissing &
           + IBITS( earthRadQF(iwlArray(iwl)), 0, 1 )
      radBadPix   = radBadPix   &
           + IBITS( earthRadQF(iwlArray(iwl)), 1, 1 )
      radError   = radError   &
           + IBITS( earthRadQF(iwlArray(iwl)), 2, 1 )
      radWarning = radWarning &
           + IBITS( earthRadQF(iwlArray(iwl)), 3, 8 )
    ENDDO

    irrMissing = 0
    irrError   = 0
    irrWarning = 0
    DO iwl = 1, nwl
      irrMissing = irrMissing &
           + IBITS( solarIrrQF(iwlArray(iwl)), 0, 1 )
      irrError   = irrError   &
           + IBITS( solarIrrQF(iwlArray(iwl)), 1, 2 )
      irrWarning = irrWarning &
           + IBITS( solarIrrQF(iwlArray(iwl)), 3, 8 )
    ENDDO

    IF( descendQ )       QAflags = IBSET( QAflags, 3  )
    IF(anomflg_3 > 0)    QAflags = IBSET( QAflags, 6 )
    !   to use KMNI L1B XtrackQualityFlags instead, use:
    !   IF (IAND(anomflg, mask7) > 0)) QAflags = IBSET( QAflags, 6 )
    IF( PclimQ )         QAflags = IBSET( QAflags, 7  )
    IF( IBITS( gflg, 6, 1 ) == 1 ) &
         QAflags = IBSET( QAflags, 8  )
    IF( sza  >= 88.    ) QAflags = IBSET( QAflags, 9  )
    IF( radMissing > 0 ) QAflags = IBSET( QAflags, 10 )
    IF( radError   > 0 ) QAflags = IBSET( QAflags, 11 )
    IF( radWarning > 0 ) QAflags = IBSET( QAflags, 12 )
    IF( irrMissing > 0 ) QAflags = IBSET( QAflags, 13 )
    IF( irrError   > 0 ) QAflags = IBSET( QAflags, 14 )
    IF( irrWarning > 0 ) QAflags = IBSET( QAflags, 15 )

    IF( radMissing > 0 ) L2param%NumberOfRadianceMissing &
         = L2param%NumberOfRadianceMissing + 1
    IF( radError   > 0 ) L2param%NumberOfRadianceError &
         = L2param%NumberOfRadianceError + 1
    IF( radWarning > 0 ) L2param%NumberOfRadianceWarning &
         = L2param%NumberOfRadianceWarning + 1
    IF( irrMissing > 0 ) L2param%NumberOfIrradianceMissing &
         = L2param%NumberOfIrradianceMissing + 1
    IF( irrError   > 0 ) L2param%NumberOfIrradianceError &
         = L2param%NumberOfIrradianceError + 1
    IF( irrWarning > 0 ) L2param%NumberOfIrradianceWarning &
         = L2param%NumberOfIrradianceWarning + 1

    IF( sza >= 88.0 .OR. radMissing > 0 .OR. radError > 0 &
         .OR. irrMissing > 0 .OR. irrError > 0 ) THEN
      skipit = .TRUE.
      L2param%NumberOfSkippedSamples = L2param%NumberOfSkippedSamples + 1
      algflg = 0
    ELSE
      skipit = .FALSE.
    ENDIF

    IF( radMissing > 0 .OR. irrMissing > 0 ) THEN
      L2param%NumberOfMissingInputSamples = &
           L2param%NumberOfMissingInputSamples + 1
    ENDIF

    IF( radError > 0 .OR. irrError > 0 ) THEN
      L2param%NumberOfBadInputSamples = &
           L2param%NumberOfBadInputSamples + 1
    ENDIF

    IF( radWarning > 0 .OR. irrWarning > 0 ) THEN
      L2param%NumberOfInputWarningSamples = &
           L2param%NumberOfInputWarningSamples + 1
    ENDIF

    !! Set flag when radiance bad pixel flag from L1b is encountered
    !! This is done for tracking. The algorithm proceeds nonetheless.
    IF( radBadPix > 0 ) THEN
      DO iwl = 1, nwl
        IF( IBITS( earthRadQF(iwlArray(iwl)), 1, 1 ) > 0 ) THEN
          radBadPixflgs = IBSET( radBadPixflgs, iwlArray(iwl)-1 )
        ENDIF
      ENDDO
    ENDIF

    IF( IBITS( QAflags, 8, 4 ) + IBITS( QAflags, 13, 2 ) > 0 ) THEN
      badPixQ = .TRUE.
    ELSE
      badPixQ = .FALSE.
    ENDIF

    IF( (.NOT. descendQ) .AND. (.NOT. badPixQ ) ) THEN
      L2param%NumberOfGoodInputSamples = &
           L2param%NumberOfGoodInputSamples + 1
      IF( sza > 84. ) L2param%NumberOfLargeSZAInputSamples  &
           = L2param%NumberOfLargeSZAInputSamples + 1
    ENDIF

    RETURN
  END FUNCTION O3T_setQAflgsI




  FUNCTION O3T_setQAflgsP( irflo, irfhi, iozon, imixr, glint, maxitr, &
       algflg, resn, flg3lm, flg4lm, soilim, &
       sza, pathl, so2ind ) RESULT( errflgs )
!       sza, pathl, so2ind, ozbst ) RESULT( errflgs )

    implicit none

    INTEGER( KIND=4 ), INTENT(IN) :: irflo, irfhi, iozon, imixr
    INTEGER( KIND=1 ), INTENT(IN) :: algflg
    REAL (KIND=4), INTENT(IN)     :: sza, pathl
    REAL (KIND=4), INTENT(IN)     :: flg3lm(3), flg4lm(3), soilim
    REAL (KIND=4), DIMENSION(:), INTENT(IN)  :: resn
    REAL (KIND=4), INTENT(INOUT)  :: so2ind !, ozbst
    LOGICAL, INTENT(IN)           :: glint, maxitr
    INTEGER( KIND=1 )             :: errflgs
    INTEGER( KIND=4 )             :: i, modalg
    INTEGER( KIND=4 )             :: nwl

    ! -- set error flag to 0
    errflgs = 0

    ! -- set glint flag
    IF( glint ) errflgs = 1

    ! -- set sza flag
    IF( sza > 84. ) errflgs = 2
    !
    ! -- flag 3 for high 360 nm residue
    !
    modalg = MOD(algflg,10)
    IF( modalg == 1) THEN
      IF( abs(resn(irfhi)) >  flg3lm(1) ) errflgs = 3
    ELSE IF( modalg == 2 ) THEN
      IF( abs(resn(irfhi)) >  flg3lm(2) ) errflgs = 3
    ELSE IF( modalg == 3 ) THEN
      IF( abs(resn(irfhi)) >  flg3lm(3) ) errflgs = 3
    ENDIF
    !
    ! -- flag 4 for ozone inconsistency
    !
    IF( modalg == 1 ) THEN
      IF( abs(resn(imixr)) > flg4lm(1) ) errflgs = 4
    ELSE IF( modalg == 2 ) THEN
      IF( abs(resn(irflo)) > flg4lm(2) ) errflgs = 4
    ELSE IF(modalg == 3) THEN
      IF( abs(resn(iozon)) > flg4lm(3) ) errflgs = 4
    ENDIF
    !
    ! -- check for so2 contamination (this code uses soi from
    !    subroutine soi, which is v7 algorithm relative to 360)
    !
    IF( so2ind >= soilim .AND. pathl < 3.5 ) errflgs = 5

    IF( maxitr ) errflgs = 6
    !
    ! -- set "fatal" flag if residue has hit the limit
    !

    nwl = SIZE( resn )
    IF( nwl > 20 ) THEN  !! a large nubmer of residuals
      IF( ABS(resn(imixr)) >= 32.0 .OR.  ABS(resn(iozon)) >= 32.0 .OR. &
           ABS(resn(irflo)) >= 32.0 .OR.  ABS(resn(irfhi)) >= 32.0  ) THEN
        errflgs = 7
      ENDIF
    ELSE  !! a small number of residuals
      DO i = imixr, SIZE( resn )
        IF( abs(resn(i)) >= 32.0 ) THEN
          errflgs = 7
        ENDIF
      ENDDO
    ENDIF
  END FUNCTION O3T_setQAflgsP

  SUBROUTINE O3T_QAcounter( algflg, QAflags, QualityFlagsCounters )

    implicit none

    INTEGER(KIND=1),INTENT(IN) :: algflg
    INTEGER(KIND=2),INTENT(IN) :: QAflags
    INTEGER(KIND=4),DIMENSION(:,:), INTENT(INOUT) :: QualityFlagsCounters
    INTEGER(KIND=4) :: xalg, xqa
    INTEGER(KIND=4) :: di

    xalg = MOD( algflg, 10 )
    xqa  = MOD( IBITS(QAflags, 0, 4), 8 ) + 1

    IF( xalg > 0 .AND. xalg < 4 .AND. xqa > 0 .AND. xqa <= 8 ) THEN
      QualityFlagsCounters(xalg,xqa) = QualityFlagsCounters(xalg,xqa) + 1
    ENDIF

    IF( xalg > 0 .AND. xalg < 4 ) THEN
      DO di = 9, 16
        QualityFlagsCounters(xalg,di) = QualityFlagsCounters(xalg,di) &
             + IBITS( QAflags, di-1, 1 )
      ENDDO
    ENDIF
  END SUBROUTINE O3T_QAcounter

  SUBROUTINE O3T_zonalMinMax( QAflags, lat, ozone, L2param )

    implicit none

    INTEGER (KIND=2), INTENT(IN) :: QAflags
    REAL (KIND=4), INTENT(IN) :: lat, ozone
    TYPE (L2PARAM_T), INTENT(INOUT) :: L2param
    INTEGER (KIND=4) :: iz
    INTEGER( KIND=1 ) :: errflgs

    errflgs = int(IBITS( QAflags, 0, 4 ),kind(errflgs))
    IF( errflgs == 0 .OR. errflgs == 1 .OR. &
         errflgs == 8 .OR. errflgs == 9 ) THEN
      DO iz = 1, nZones
        IF(        L2param%ZonalLatRange(1,iz) <= lat .AND. &
             lat <= L2param%ZonalLatRange(2,iz) ) THEN
          IF( ozone < L2param%ZonalOzoneMin(iz) ) THEN
            L2param%ZonalOzoneMin(iz) = ozone
          ELSE IF ( ozone > L2param%ZonalOzoneMax(iz) ) THEN
            L2param%ZonalOzoneMax(iz) = ozone
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE O3T_zonalMinMax
END MODULE O3T_QA_class
