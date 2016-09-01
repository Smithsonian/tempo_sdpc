!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE OMI_metaData_class
! 
!  contains functions to read and write the ECS metadata 
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
! 
!!Return
! None 
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
! GEST/UMBC
! email: Kai.Yang.1@gsfc.nasa.gov
! 
!!Design Notes
!
!!END
!!****************************************************************************
MODULE OMI_metaData_class
    USE PGS_TD_class
    USE PGS_MET_class
    USE OMI_SMF_class
    USE OMSAO_errstat_module
    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, four = 4
    INTEGER (KIND=4), PARAMETER :: NATTR_MAX = 20
    INTEGER (KIND=4), PARAMETER :: Npsa_MAX = 30
    INTEGER (KIND=4), PARAMETER :: Nelm_MAX    = 20
    INTEGER (KIND=4), PARAMETER :: ORBITNUMBER_LUN = 200200 

    TYPE, PUBLIC :: ECSMETA_ITEM_T
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: Name
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: DataType
      INTEGER ( KIND = 4 ) :: GROUP
      INTEGER ( KIND = 4 ) :: NUM_VAL
      CHARACTER(LEN=20*PGSd_MET_MAX_STRING_SET_L) :: Values
    END TYPE ECSMETA_ITEM_T

    TYPE, PUBLIC :: OMIECSMETA_T
      !! Core MetaData String Field
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: RangeBeginningDate,   &
                                                  RangeBeginningTime,   &
                                                  RangeEndingDate,      &
                                                  RangeEndingTime,      &
                                                  EquatorCrossingTime,  &
                                                  EquatorCrossingDate,  &
                                                  ShortName,            &
                                                  AutomaticQualityFlag, &
                                                  ScienceQualityFlag

      !! value get from ECS Core to be compared with those in L2 mcf
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: AssociatedPlatformSN,   &
                                                  AssociatedInstrumentSN, &
                                                  AssociatedSensorSN

      ! Integers and double to be read from Core Meta
      INTEGER(KIND=4) :: orbitNumber, QAPercentMissingData
      REAL(KIND=8) :: EqCrossLon

      INTEGER (KIND=4) :: Npsa   !! Number of Product Specific Attributes

      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(Npsa_MAX) :: &
                                                psaNames
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), &
                     DIMENSION(Npsa_MAX,Nelm_MAX) ::  psaValues

    END TYPE OMIECSMETA_T

    INTEGER, PARAMETER :: INVENTORY = 2
    INTEGER, PARAMETER :: ARCHIVE   = 3

    !! string type of PSA to be copied from Core Meta Data
    CHARACTER( LEN = 1000 ), PRIVATE ::                       &
     PSAtoBeCopiedFromCore= "SouthAtlanticAnomalyCrossing" // &
                            "SolarEclipse"                 // &
                            "PathNr"                       // &
                            "StartBlockNr"                 // &
                            "EndBlockNr"                   

     !! The list below contains some PSAs that are not in the OMSO2
     !! ESDT Descriptor. The correct list is above.
     !PSAtoBeCopiedFromCore= "NrMeasurements"               // &
     !                       "NrZoom"                       // &
     !                       "NrSpatialZoom"                // &
     !                       "NrSpectralZoom"               // &
     !                       "ExpeditedData"                // &
     !                       "SouthAtlanticAnomalyCrossing" // &
     !                       "SpacecraftManeuverFlag"       // &
     !                       "SolarEclipse"                 // &
     !                       "MasterClockPeriods"           // &
     !                       "ExposureTimes"                // &
     !                       "PathNr"                       // &
     !                       "StartBlockNr"                 // &
     !                       "EndBlockNr"                   // &
     !                       "InstrumentConfigurationIDs"

    CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: &
        AutomaticQualityFlagExplanation = &
        "Flag set to Passed if QAPercentHighQualityData >= 90%, "//&
        "Flag set to Suspect if QAPercentHightQualityData >= 60% or "//&
        "or L1B does not have its AutomaticQualityFlag set to Passed, "//&
        "Flag set to Failed if QAPercentHighQualityData < 60%. "

    CHARACTER( LEN = 80 ):: inventoryMetadataName
    CHARACTER( LEN = 80 ):: archivedMetadataName
    LOGICAL (KIND=4) :: OMIMETADATA_READ = .FALSE.

    PUBLIC :: OMI_getCoreMetaData
    PUBLIC :: OMI_setCoreArchMetaData

    CONTAINS

      FUNCTION OMI_getCoreMetaData( OMI_fileLUN, versionIN, &
                                            OMIcoreMeta,  &
                                            year, month, day, jday ) &
                                            RESULT( status ) 

        use utilities, only: day_of_year

        INTEGER (KIND=4), INTENT(IN) :: OMI_fileLUN, versionIN
        TYPE (OMIECSMETA_T), INTENT(OUT) :: OMIcoreMeta
        INTEGER (KIND=4), INTENT(OUT) :: year, month, day, jday
        INTEGER (KIND=4) :: status!, iName
        INTEGER (KIND=4) :: ierr, version
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: dummyName!, dummyValue
        INTEGER (KIND=4) :: counter_r
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
        CHARACTER( LEN = 1 ) :: foo
        !INTEGER (KIND=4), EXTERNAL :: day_of_year
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), &
                     DIMENSION(Nelm_MAX) ::  tmp_psaValues


        OMIMETADATA_READ = .FALSE.

        !! CoreMetadata in OMI can sometimes written in several
        !! variations, such as "coremetadata", or "CoreMetadata.0"
        !! Here it figure out exactly how it was written in the OMI,
        !! then use the name to read the fields in the CoreMetadata.

        inventoryMetadataName = "CoreMetadata.0"

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEBEGINNINGDATE",             &
                                      OMIcoreMeta%RangeBeginningDate )
        IF( status == PGS_S_SUCCESS ) THEN
           READ( OMIcoreMeta%RangeBeginningDate , '(I4,A1,I2,A1,I2)' ) &
                 year, foo, month, foo, day
           jday = day_of_year( year, month, day )
           WRITE( msg, '(A, I5,I3,I3,I4)' ) "(y,m,d) = ", &
                       year, month, day, jday 
           ierr = OMI_SMF_setmsg( OMI_W_GENERAL, msg, &
                                  "OMI_getCoreMetaData",four )
        ELSE
           OMIcoreMeta%RangeBeginningDate = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                 "can't get RangeBeginningDate" &
                              // " from OMI file core Meta data.", &
                                 "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEBEGINNINGTIME",             &
                                      OMIcoreMeta%RangeBeginningTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%RangeBeginningTime = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                              "can't get RangeBeginningTime" &
                           // " from OMI file core Meta data.", &
                              "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEENDINGDATE",                &
                                      OMIcoreMeta%RangeEndingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%RangeEndingDate = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                            "can't get RangeEndingDate" &
                         // " from OMI file core Meta data.", &
                            "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEENDINGTIME",                &
                                      OMIcoreMeta%RangeEndingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%RangeEndingTime = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                 "can't get RangeEndingTime" &
                              // " from OMI file core Meta data.", &
                                 "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version,            &
                                      inventoryMetadataName,            & 
                                      "EQUATORCROSSINGTIME.1",          &
                                      OMIcoreMeta%EquatorCrossingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%EquatorCrossingTime = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                              "can't get EquatorCrossingTime" &
                           // " from OMI file core Meta data.", &
                              "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version,            &
                                      inventoryMetadataName,            & 
                                      "EQUATORCROSSINGDATE.1",          &
                                      OMIcoreMeta%EquatorCrossingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%EquatorCrossingDate = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                  "can't get EquatorCrossingDate" &
                               // " from OMI file core Meta data.", &
                                  "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version,            &
                                      inventoryMetadataName,            & 
                                      "SHORTNAME",                      &
                                      OMIcoreMeta%ShortName )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%ShortName = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, "can't get ShortName" &
                                  // " from OMI file core Meta data.", &
                                  "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version , &
                                      inventoryMetadataName,& 
                                     "ASSOCIATEDPLATFORMSHORTNAME.1", & 
                                      OMIcoreMeta%AssociatedPlatformSN  )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%AssociatedPlatformSN = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                             "can't get ASSOCIATEDPLATFORMSHORTNAME" &
                          // " from OMI file core Meta data.", &
                             "OMI_getCoreMetaData", zero )
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version , &
                                      inventoryMetadataName,& 
                                     "ASSOCIATEDINSTRUMENTSHORTNAME.1", &
                                      OMIcoreMeta%AssociatedInstrumentSN  )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%AssociatedInstrumentSN = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                             "can't get ASSOCIATEDINSTRUMENTSHORTNAME" &
                                  // " from OMI file core Meta data.", &
                                  "OMI_getCoreMetaData", zero )
        ENDIF

        version = versionIN
        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version , &
                                      inventoryMetadataName,& 
                                     "ASSOCIATEDSENSORSHORTNAME.1", &
                                      OMIcoreMeta%AssociatedSensorSN  )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%AssociatedSensorSN = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                  "can't get ASSOCIATEDSENSORSHORTNAME" &
                                  // " from OMI file core Meta data.", &
                                  "OMI_getCoreMetaData", zero )
        ENDIF

        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version , &
                                      inventoryMetadataName,& 
                                     "AUTOMATICQUALITYFLAG.1", &
                                      OMIcoreMeta%AutomaticQualityFlag )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%AutomaticQualityFlag = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                  "can't get AUTOMATICQUALITYFLAG" &
                                  // " from OMI file core Meta data.", &
                                  "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_GetPCAttr_s( OMI_fileLUN, version , &
                                      inventoryMetadataName,& 
                                     "SCIENCEQUALITYFLAG.1", &
                                      OMIcoreMeta%ScienceQualityFlag )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%AutomaticQualityFlag = ""
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                  "can't get SCIENCEQUALITYFLAG" &
                                  // " from OMI file core Meta data.", &
                                  "OMI_getCoreMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! Read Integer Type of OMI data fileds.
        ! 1. Get OrbitNumber
        version = versionIN
        status = PGS_MET_GetPCAttr_i( OMI_fileLUN, version , &
                                      inventoryMetadataName,& 
                                     "ORBITNUMBER.1", &  
                                      OMIcoreMeta%orbitNumber )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%orbitNumber = -1
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                 "get OrbitNumber.1 failed.", &
                                 "OMI_getCoreMetaData", zero )  
           status = OMI_E_FAILURE
           return
        ENDIF

        ! 2. Get QAPercentMissingData
        version = versionIN
        status = PGS_MET_GetPCAttr_i( OMI_fileLUN, version , &
                                      inventoryMetadataName,& 
                                     "QAPERCENTMISSINGDATA.1", &
                                      OMIcoreMeta%QAPercentMissingData )
        IF( status /= PGS_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                 "get QAPercentMissingData.1 failed.", &
                                 "OMI_getCoreMetaData", zero )  
        ENDIF

        !! Read Double Type of OMI data fileds (there is just one of them).
        version = versionIN
        status = PGS_MET_GetPCAttr_d( OMI_fileLUN, version, &
                                      inventoryMetadataName,& 
                                     "EQUATORCROSSINGLONGITUDE.1", &
                                      OMIcoreMeta%EqCrossLon )
        IF( status /= PGS_S_SUCCESS ) THEN
           OMIcoreMeta%EqCrossLon = -999.0
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                     "get EQUATORCROSSINGLONGITUDE.1 failed.", &
                     "OMI_getCoreMetaData", zero )  
        ELSE          ! write to the Log files the field values of OMI file
           WRITE(msg, '(A,F9.5)' ) "EQUATORCROSSINGLONGITUDE.1 = ", &
                                   OMIcoreMeta%EqCrossLon
           ierr = OMI_SMF_setmsg( OMI_W_GENERAL, msg, &
                                 "OMI_getCoreMetaData",four )
        ENDIF

        counter_r = 1  !! counter for reading PSA from ECS
        version = 1
        WRITE( dummyName, '(A,I1)') "ADDITIONALATTRIBUTENAME.", counter_r
        status = PGS_MET_getPCAttr_s( OMI_fileLUN, version,               &
                                      inventoryMetadataName,               &
                                      dummyName,                           &
                                      OMIcoreMeta%psaNames(counter_r) )
        DO WHILE( status == PGS_S_SUCCESS )
          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "PARAMETERVALUE.", counter_r
          ELSE IF( counter_r < 100 ) THEN
             WRITE( dummyName, '(A,I2)') "PARAMETERVALUE.", counter_r
          ELSE
             ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                    "Too many PSA in input file", &
                                    "OMI_getCoreArchivedMetaData", zero )
             status = OMI_E_FAILURE
             RETURN
          ENDIF
          
          !! empty string Array before it is read in from OMI metadata
          OMIcoreMeta%psaValues(counter_r,:) = ""
          tmp_psaValues(:) = ""
          version = 1

          ! FIXME - masking array temporary
          status = PGS_MET_getPCAttr_s( OMI_fileLUN, version,  &
                                        inventoryMetadataName,  &
                                        dummyName,              &
!                                        OMIcoreMeta%psaValues(counter_r,:) )
                                        tmp_psaValues )
          OMIcoreMeta%psaValues(counter_r,:) = tmp_psaValues

          IF( status /= PGS_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OMI_E_INPUT, "GEt PSA value failed", &
                                   "OMI_getCoreArchivedMetaData", zero )
             status = OMI_E_FAILURE
             RETURN
          ENDIF

          OMIcoreMeta%Npsa   = counter_r
          counter_r = counter_r + 1;

          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "ADDITIONALATTRIBUTENAME.", counter_r
          ELSE IF( counter_r < 100 ) THEN
             WRITE( dummyName, '(A,I2)') "ADDITIONALATTRIBUTENAME.", counter_r
          ELSE
             ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                   "Too many PSA in input file",&
                                   "OMI_getCoreMetaData", zero )
             status = OMI_E_FAILURE
             RETURN
          ENDIF
          status = PGS_MET_getPCAttr_s( OMI_fileLUN, version,               &
                                        inventoryMetadataName,               &
                                        dummyName,                           &
                                        OMIcoreMeta%psaNames(counter_r) )
        END DO

        OMIMETADATA_READ = .TRUE.
        status = OMI_S_SUCCESS
      END FUNCTION OMI_getCoreMetaData


      FUNCTION OMI_setCoreArchMetaData( L2_outFile_LUN , OMIcoreMeta, &
                                        LUNinputPointer, mcf_LUN,  &
                                        L2specificItems, ShortName ) RESULT( status ) 
        USE OMSAO_he5_module
        USE PGS_PC_class
        INTEGER (KIND=4), INTENT(IN) :: L2_outFile_LUN
        INTEGER (KIND=4), INTENT(IN) :: mcf_LUN
        TYPE (OMIECSMETA_T), INTENT(IN) :: OMIcoreMeta
        TYPE (ECSMETA_ITEM_T), DIMENSION(:), INTENT(IN) :: L2specificItems 
        INTEGER (KIND=4), DIMENSION(:), INTENT(IN) :: LUNinputPointer
        CHARACTER(LEN=*), INTENT(OUT) :: ShortName
        INTEGER (KIND=4) :: status, I4foo
        INTEGER (KIND=4) :: nL2items
        CHARACTER ( LEN = 1 ) :: delim = ','
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), &
                DIMENSION( SIZE(LUNinputPointer) ) :: inputPointer
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), &
                DIMENSION( 100 )  :: StringValues

        !! value from L2 mcf to be compared with those in CoreMeta
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) ::  &
                         MCF_AssociatedPlatformSN,   &
                         MCF_AssociatedInstrumentSN, &
                         MCF_AssociatedSensorSN

        INTEGER (KIND=4) :: orbitNumber
        INTEGER (KIND=4) :: VersionID

        CHARACTER( LEN = PGSd_PC_FILE_PATH_MAX ) :: InstrumentPlatformName
        INTEGER, EXTERNAL :: OMI_localGranuleID
        CHARACTER( LEN = PGSd_PC_FILE_PATH_MAX ) :: localGranuleID, &
                                                    L2_filename 
        INTEGER (KIND=4) :: in_s, ii, num_val, version, ierr
        INTEGER (KIND=4) :: counter_w, nStr!, nVal
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: dummyName!, dummyValue
        INTEGER (KIND=4) :: HE5id
        REAL    (KIND=4) :: R4foo
        REAL    (KIND=8) :: R8foo
        INTEGER (KIND=4), DIMENSION(100) :: I4fooArray
        REAL    (KIND=4), DIMENSION(100) :: R4fooArray
        REAL    (KIND=8), DIMENSION(100) :: R8fooArray
        CHARACTER( LEN = 28 ) :: RangeBeginningDateTime!, ProductionDateTime
        CHARACTER( LEN=PGSd_MET_GROUP_NAME_L ) :: &
             GROUPS(PGSd_MET_NUM_OF_GROUPS)
        CHARACTER( LEN=PGS_SMF_MAX_MSG_SIZE ) :: msg
        character(LEN=PGSd_MET_MAX_STRING_SET_L), dimension(Nelm_MAX) :: &
             tmp_psaValues
        integer :: idiot

        inputPointer(:) = ""
        DO ii = 1, SIZE(LUNinputPointer) 
          version = 1
          status = PGS_PC_GetReference( LUNinputPointer(ii), version, msg )
          IF( status /= PGS_S_SUCCESS ) THEN
             WRITE( msg,'(A,I8)' ) "get filename failed at LUN = ", &
                                  LUNinputPointer(ii)
             ierr = OMI_SMF_setmsg( status,msg,"OMI_setCoreArchMetaData", zero )
             status = OMI_E_FAILURE
             RETURN
          ELSE
             in_s = INDEX( msg, '/', BACK = .TRUE. ) + 1
             inputPointer(ii) = TRIM( msg( in_s:) )
          ENDIF
        ENDDO 

        !! read the Orbit Number from the PCF file, later it will
        !! be compared with the Orbit Number retrieved from CoreMetaData
        !! to make sure they are the same, otherwise this code exit 
        !! with an error. Doing so will ensure that intended input 
        !! file is correctly staged.

        status = PGS_PC_GetConfigData( ORBITNUMBER_LUN, msg )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I8)' ) "get orbit number failed at LUN = ", &
                                ORBITNUMBER_LUN  
           ierr = OMI_SMF_setmsg( status, msg, "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ELSE
           READ( msg, '(I8)'), orbitNumber
           IF( orbitNumber /= OMIcoreMeta%orbitNumber ) THEN
              WRITE( msg,'(A,I8,A,I8,A)' ) "Input file OrbitNumber = ", &
               OMIcoreMeta%orbitNumber, ", differs in PCF = ", orbitNumber, &
               ", use orbitNumber in input L2 file."
              ierr = OMI_SMF_setmsg( OMI_W_GENERAL, msg, &
                                    "OMI_setCoreArchMetaData", zero )
              orbitNumber = OMIcoreMeta%orbitNumber
           ENDIF
        ENDIF

        !! combined RangeBeginningDate and RangeBeginningTime to get
        !! RangeBeginningDateTime, which will be used in seting up
        !! the LocalGranuleID.
        WRITE( RangeBeginningDateTime, '(A)' ) &
               TRIM( OMIcoreMeta%RangeBeginningDate ) // "T" // &
               TRIM( OMIcoreMeta%RangeBeginningTime ) // "Z" 

        GROUPS(:) = ""
        status = PGS_MET_init( mcf_LUN, GROUPS )
        IF( status /= PGS_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, "PGS_MET_init failed.", &
                                  "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! Read SHORTNAME which is set in the MCF file.
        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "SHORTNAME", ShortName )
        IF( status /= PGS_S_SUCCESS ) THEN
            ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                      "PGS_MET_getSetAttr_s failed for SHORTNAME", &
                                  "OMI_setCoreArchMetaData", zero )
            status = OMI_E_FAILURE
            RETURN
        ENDIF

        !! Read VERSIONID, a integer with no more than 3 digit, again 
        !! it is set in the MCF file
        status = PGS_MET_getSetAttr_i( GROUPS(INVENTORY), &
                                      "VERSIONID", VersionID )
        IF( status /= PGS_S_SUCCESS ) THEN
            ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                      "PGS_MET_getSetAttr_s failed for VERSIONID", &
                                  "OMI_setCoreArchMetaData", zero )
            status = OMI_E_FAILURE
            RETURN
        ENDIF
                         
        !! Check the make sure the Associated Sensor stuff in MCF are the same
        !! in the input OMOcoreMeta. If different, generate warning message 
        !! in Log file, but use the values in MCF to construct Local granule ID.
        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "ASSOCIATEDPLATFORMSHORTNAME.1", & 
                                       MCF_AssociatedPlatformSN )
        IF( TRIM(OMIcoreMeta%AssociatedPlatformSN) /= &
            TRIM(MCF_AssociatedPlatformSN) ) THEN
            WRITE( msg,'(A)' ) "AssociatedPlatformSN not matched,Core:"//&
                  TRIM(OMIcoreMeta%AssociatedPlatformSN)// " MCF:" // &
                  TRIM(MCF_AssociatedPlatformSN)
            ierr = OMI_SMF_setmsg( OMI_W_GENERAL, msg, &
                                  "OMI_setCoreArchMetaData", four )
        ENDIF

        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "ASSOCIATEDINSTRUMENTSHORTNAME.1", & 
                                       MCF_AssociatedInstrumentSN )
        IF( TRIM(OMIcoreMeta%AssociatedInstrumentSN) /= &
            TRIM(MCF_AssociatedInstrumentSN) ) THEN
           WRITE( msg,'(A)') &
                  "AssociatedInstrumentShortname not matched,Core:"//&
                  TRIM(OMIcoreMeta%AssociatedInstrumentSN)//" MCF:"// &
                  TRIM(MCF_AssociatedInstrumentSN)
           ierr = OMI_SMF_setmsg( OMI_W_GENERAL, msg, &
                                  "OMI_setCoreArchMetaData", four )
        ENDIF
        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "ASSOCIATEDSENSORSHORTNAME.1", &
                                       MCF_AssociatedSensorSN )
        IF( TRIM(OMIcoreMeta%AssociatedSensorSN) /=  &
            TRIM(MCF_AssociatedSensorSN) ) THEN
           WRITE( msg,'(A)' ) &
                  "AssociatedSensorShortname not matched,Core:"//&
                  TRIM(OMIcoreMeta%AssociatedSensorSN )//" MCF:"// &
                  TRIM(MCF_AssociatedSensorSN )
           ierr = OMI_SMF_setmsg( OMI_W_GENERAL, msg, &
                                  "OMI_setCoreArchMetaData", four )
        ENDIF

        !! Construct the LocalGranuleID for L2 output
        InstrumentPlatformName =  TRIM(MCF_AssociatedInstrumentSN ) // &
                                  "-" // TRIM(MCF_AssociatedPlatformSN ) 
        status = OMI_localGranuleID( RangeBeginningDateTime, &
                                     InstrumentPlatformName, &
                                    "L2", ShortName,  orbitNumber, &
                                     VersionID, "he5", localGranuleID );
        IF( status /= OMI_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "set local granule id failed.", &
                                  "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN 
        ENDIF

        !! set the local granule ID metadata
        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), "LOCALGRANULEID", &
                                    TRIM( localGranuleID )  )
        IF( status /= PGS_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                  "PGS_MET_setAttr_s failed for LOCALGRANULEID",&
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEBEGINNINGDATE", &
                                   OMIcoreMeta%RangeBeginningDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "RANGEBEGINNINGDATE"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEBEGINNINGTIME", &
                                   OMIcoreMeta%RangeBeginningTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "RANGEBEGINNINGTIME"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEENDINGDATE", & 
                                   OMIcoreMeta%RangeEndingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                               "RANGEENDINGDATE"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEENDINGTIME", &
                                   OMIcoreMeta%RangeEndingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "RANGEENDINGTIME"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "EQUATORCROSSINGDATE.1", &
                                   OMIcoreMeta%EquatorCrossingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                               "EQUATORCROSSINGDATE.1"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "EQUATORCROSSINGTIME.1", &
                                   OMIcoreMeta%EquatorCrossingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "EQUATORCROSSINGTIME.1"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_i( GROUPS(INVENTORY), &
                                    "ORBITNUMBER.1", orbitNumber )
        IF( status /= PGS_S_SUCCESS ) THEN
         WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for ORBITNUMBER.1"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! Write the Input Pointer Metadata Field
        status = PGS_MET_SetMultiAttr_s( GROUPS(INVENTORY),&
                                         "InputPointer", &
                                        SIZE(LUNinputPointer), inputPointer )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "InputPointer"
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! Now set the 1 (double type) core metadata field 
        !! that is supposed to be set by the PGE. 
        !! Write the EQUATORCROSSINGLONGITUDE Metadata Field

        status = PGS_MET_setAttr_d( GROUPS(INVENTORY), &
                                   "EQUATORCROSSINGLONGITUDE.1", &
                                   OMIcoreMeta%EqCrossLon )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_d failed for "// &
                              "EQUATORCROSSINGLONGITUDE.1"
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, &
                                 "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! copy some PSA from OMIcoreMeta.

        counter_w = 1  !! counter for writing of L2 PSA
        DO ii = 1, OMIcoreMeta%Npsa
          IF( INDEX( TRIM(PSAtoBeCopiedFromCore), &
              TRIM(OMIcoreMeta%psaNames(ii)) ) == 0 ) CYCLE
          IF( counter_w < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "ADDITIONALATTRIBUTENAME.",counter_w
          ELSE IF( counter_w < 100 ) THEN
             WRITE( dummyName, '(A,I2)') "ADDITIONALATTRIBUTENAME.",counter_w
          ELSE
             ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                   "Too many PSA in input file",&
                                  "L2_setCoreArchMetaData", zero )
             status = OMI_E_FAILURE
             RETURN
          ENDIF
          status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                      TRIM( dummyName ), &
                                      TRIM(OMIcoreMeta%psaNames(ii)) )

          !! figure out the number of strings in the metadata, and
          !! write the exact same number of stings in the L2 metadata.
          nStr = 0
          !original counter would loop to infinity when high values of
          !psaValues were empty.
          idiot = 0
!          DO WHILE( nStr+1 <= Nelm_MAX )
          DO WHILE( idiot+1 <= Nelm_MAX )
            if (LEN_TRIM(OMIcoreMeta%psaValues(ii,nStr+1))  > 0 ) then
              nStr = nStr + 1
            endif
            idiot = idiot +1
          ENDDO

          IF( counter_w < 10 ) THEN
            WRITE( dummyName, '(A,I1)') "PARAMETERVALUE.", counter_w
          ELSE IF( counter_w < 100 ) THEN
            WRITE( dummyName, '(A,I2)') "PARAMETERVALUE.", counter_w
          ENDIF

          ! FIXME - masking array temporary
          tmp_psaValues=OMIcoreMeta%psaValues(ii,:)
          status = PGS_MET_SetMultiAttr_s( GROUPS(INVENTORY), &
                                           TRIM( dummyName ), &
                                           nStr, &
!                                           OMIcoreMeta%psaValues(ii,:) )
                                           tmp_psaValues(:))
          counter_w = counter_w + 1
        END DO

        nL2items = SIZE( L2specificItems ) 
        DO ii = 1, nL2items
          in_s = L2specificItems(ii)%GROUP          
          num_val = L2specificItems(ii)%NUM_VAL 
          IF( num_val < 1 ) THEN 
             CYCLE
          ELSE IF( num_val == 1 ) THEN 
             IF( INDEX( TRIM(L2specificItems(ii)%DataType), "STRING")>0 ) THEN 
                status = PGS_MET_setAttr_s( GROUPS(in_s), &
                                            TRIM(L2specificItems(ii)%Name), &
                                            TRIM(L2specificItems(ii)%Values) )
             ELSE IF( INDEX( TRIM(L2specificItems(ii)%DataType), &
                      "INTEGER")>0 ) THEN 
                READ( L2specificItems(ii)%Values, *) I4foo 
                status = PGS_MET_setAttr_i( GROUPS(in_s), &
                                            TRIM(L2specificItems(ii)%Name), &
                                            I4foo )
             ELSE IF( INDEX( TRIM(L2specificItems(ii)%DataType), &
                      "DOUBLE")>0 ) THEN 
                READ( L2specificItems(ii)%Values, *) R8foo 
                status = PGS_MET_setAttr_d( GROUPS(in_s), &
                                            TRIM(L2specificItems(ii)%Name), &
                                            R8foo )
             ELSE IF( INDEX( TRIM(L2specificItems(ii)%DataType), &
                      "REAL")>0 ) THEN 
                READ( L2specificItems(ii)%Values, *) R4foo 
                status = PGS_MET_setAttr_r( GROUPS(in_s), &
                                            TRIM(L2specificItems(ii)%Name), &
                                            R4foo )
             ENDIF   
          ELSE IF( num_val < 100 ) THEN 
             IF( INDEX( TRIM(L2specificItems(ii)%DataType), "INTEGER")>0 ) THEN 
                READ( L2specificItems(ii)%Values, *) I4fooArray(1:num_val) 
                status = PGS_MET_SetMultiAttr_i( GROUPS(in_s), &
                                             TRIM(L2specificItems(ii)%Name), &
                                             num_val, I4fooArray(1:num_val) )
             ELSE IF( INDEX( TRIM(L2specificItems(ii)%DataType), &
                      "DOUBLE")>0 ) THEN 
                READ( L2specificItems(ii)%Values, *) R8fooArray(1:num_val) 
                status = PGS_MET_SetMultiAttr_d( GROUPS(in_s), &
                                            TRIM(L2specificItems(ii)%Name), &
                                             num_val, R8fooArray(1:num_val) )
             ELSE IF( INDEX( TRIM(L2specificItems(ii)%DataType), &
                      "REAL")>0 ) THEN 
                READ( L2specificItems(ii)%Values, *) R4fooArray(1:num_val) 
                status = PGS_MET_SetMultiAttr_d( GROUPS(in_s), &
                                            TRIM(L2specificItems(ii)%Name), &
                                             num_val, R4fooArray(1:num_val) )
             ELSE IF( INDEX( TRIM(L2specificItems(ii)%DataType), &
                      "STRING")>0 ) THEN 
                StringValues(:) = ""
                nStr = EH_parsestrF( TRIM(L2specificItems(ii)%Values), &
                                     delim, StringValues ) 
                status = PGS_MET_SetMultiAttr_s( GROUPS(in_s), &
                                            TRIM(L2specificItems(ii)%Name), &
                                            num_val, StringValues )
             ENDIF
          ENDIF
        ENDDO

        !! Now all the Metadata are set, write them to the output file.
        !! retrieve the L2 output file name
        version = 1
        status = PGS_PC_getReference( L2_outFile_LUN, version, L2_filename )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I8)' ) "get filename from PCF file at LUN =", &
                                 L2_outFile_LUN
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, &
                                  "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! Initiate it for writing metadata 
        status = PGS_MET_sfstart( TRIM(L2_filename), HDF5_ACC_RDWR, HE5id )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_sfstart failed " // &
                               TRIM(L2_filename)
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, &
                                  "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! Write the CoreMetadata
        status = PGS_MET_write( GROUPS(INVENTORY),'CoreMetadata', HE5id )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_write CoreMetadata failed. " 
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, &
                                  "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        !! Write the Archived metadata
        status = PGS_MET_write( GROUPS(ARCHIVE),'ArchivedMetadata', HE5id )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_write ArchivedMetadata failed. " 
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, msg, &
                                  "OMI_setCoreArchMetaData", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF
        
        !! Close the L2 file 
        status = PGS_MET_SFend(HE5id)
        !! Clean up the memory associated writing metadata
        status = PGS_MET_Remove()

        status = OMI_S_SUCCESS
      END FUNCTION OMI_setCoreArchMetaData
END MODULE OMI_metaData_class
