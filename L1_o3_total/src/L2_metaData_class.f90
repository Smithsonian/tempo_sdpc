!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE L2_metaData_class
! 
! contains the function that set the L2 ECS core and archived metadata.
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
! University of Maryland Baltimore County 
! email: Kai.Yang-1@nasa.gov
! 
!!Design Notes
!
!!END
!!****************************************************************************
MODULE L2_metaData_class
    USE PGS_TD_class
    USE PGS_MET_class
    USE OMI_SMF_class
    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, four = 4
    INTEGER, PARAMETER :: INVENTORY = 2
    INTEGER, PARAMETER :: ARCHIVE   = 3

    !! string type of PSA to be copied from L1B
    CHARACTER( LEN = 1000 ), PRIVATE ::                       &
     PSAtoBeCopiedFromL1B = "NrMeasurements"               // &
                            "NrZoom"                       // &
                            "NrSpatialZoom"                // &
                            "NrSpectralZoom"               // &
                            "ExpeditedData"                // &
                            "SouthAtlanticAnomalyCrossing" // &
                            "SpacecraftManeuverFlag"       // &
                            "SolarEclipse"                 // &
                            "MasterClockPeriods"           // &
                            "ExposureTimes"                // &
                            "PathNr"                       // &
                            "StartBlockNr"                 // &
                            "EndBlockNr"                   // &
                            "InstrumentConfigurationIDs"

    CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: &
        AutomaticQualityFlagExplanation = &
        "Flag set to Passed if QAPercentHighQualityData >= 90%, "//&
        "Flag set to Suspect if QAPercentHightQualityData >= 60% or "//&
        "or input L1B do not have its AutomaticQualityFlag set to Passed, "//&
        "Flag set to Failed if QAPercentHighQualityData < 60%. "

    PUBLIC :: L2_setCoreArchMetaData

    CONTAINS
      FUNCTION L2_setCoreArchMetaData( L2_outFile_LUN , L1BRadCoreArch, &
                                       L1BIrrCoreArch,  LUNinputPointer,&
                                       mcf_LUN, L2_Parameters )         &
                                       RESULT( status ) 
        USE OMI_LUN_set
        USE PGS_PC_class
        INTEGER (KIND=4), INTENT(IN) :: L2_outFile_LUN
        INTEGER (KIND=4), INTENT(IN) :: mcf_LUN
        TYPE (L1BECSMETA_T), INTENT(IN) :: L1BRadCoreArch
        TYPE (L1BECSMETA_T), INTENT(IN) :: L1BIrrCoreArch
        INTEGER (KIND=4), DIMENSION(:), INTENT(IN) :: LUNinputPointer
        TYPE (L2PARAM_T), DIMENSION(:), INTENT(INOUT) :: L2_Parameters
        INTEGER :: npL2, nqfc 
        INTEGER (KIND=4) :: status, I4foo

        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), &
                DIMENSION( SIZE(LUNinputPointer) ) :: inputPointer

        !! value from L2 mcf to be compared with those in L1B
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
        INTEGER (KIND=4) :: in_s, ii, version, ierr, jj
        INTEGER (KIND=4) :: counter_w, nStr, nVal
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: dummyName, dummyValue
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: ShortName
        INTEGER (KIND=4) :: HE5id
        CHARACTER( LEN = 28 ) :: RangeBeginningDateTime!, ProductionDateTime
        CHARACTER( LEN=PGSd_MET_GROUP_NAME_L ) :: GROUPS(PGSd_MET_NUM_OF_GROUPS)
        CHARACTER( LEN=PGS_SMF_MAX_MSG_SIZE ) :: msg
        CHARACTER( LEN = 28 ) :: L1BRadDateTime, L1BIrrDateTime
        REAL (KIND=8) :: L1BRadTAI93, L1BIrrTAI93
        ! To avoid array temporaries, temporary arrays...
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), &
                DIMENSION(20) :: temp_L1BpsaValues
        integer (kind=4), dimension(16) :: temp_QFCounters
        real (kind=8), dimension(5) :: temp_ZLatRange
        ! Round off Zonal Ozone min & max to avoid differences
        ! in output under different compilers
        integer (kind=4), dimension(5) :: roundZOmin, roundZOmax
        real (kind=8), dimension(5) :: tempZOmin, tempZOmax

        inputPointer(:) = ""
        DO ii = 1, SIZE(LUNinputPointer) 
          version = 1
          status = PGS_PC_GetReference( LUNinputPointer(ii), version, msg )
          IF( status /= PGS_S_SUCCESS ) THEN
             WRITE( msg,'(A,I9)' ) "get filename failed at LUN = ", &
                                  LUNinputPointer(ii)
             ierr = OMI_SMF_setmsg( status,msg, "L2_setCoreArchMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ELSE
             in_s = INDEX( msg, '/', BACK = .TRUE. ) + 1
             inputPointer(ii) = TRIM( msg( in_s:) )
          ENDIF
        ENDDO 

        !! read the Orbit Number from the PCF file, later it will
        !! be compared with the Orbit Number retrieved from the 
        !! the L1B file to make sure they are the same, otherwise
        !! this code exit with an error. Doing so will ensure that
        !! intended L1B file is correctly staged.

        status = PGS_PC_GetConfigData( ORBITNUMBER_LUN, msg )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I9)' ) "get orbit number failed at LUN = ", &
                                ORBITNUMBER_LUN  
           ierr = OMI_SMF_setmsg( status, msg, "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ELSE
           READ( msg, '(I9)'), orbitNumber
           IF( orbitNumber /= L1BRadCoreArch%orbitNumber ) THEN
              WRITE( msg,'(A,I8,A,I8,A)' ) "OrbitNumber in input file = ", &
               L1BRadCoreArch%orbitNumber, ", differs in PCF = ", orbitNumber, &
               ", use orbitNumber in input L2 file."
              ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                    "L2_setCoreArchMetaData", zero )
              orbitNumber = L1BRadCoreArch%orbitNumber
           ENDIF
        ENDIF

        !! combined RangeBeginningDate and RangeBeginningTime to get
        !! RangeBeginningDateTime, which will be used in seting up
        !! the LocalGranuleID.
        WRITE( RangeBeginningDateTime, '(A)' ) &
               TRIM( L1BRadCoreArch%RangeBeginningDate ) // "T" // &
               TRIM( L1BRadCoreArch%RangeBeginningTime ) // "Z" 
        L1BRadDateTime = RangeBeginningDateTime
        WRITE( L1BIrrDateTime, '(A)' ) &
               TRIM( L1BIrrCoreArch%RangeBeginningDate ) // "T" // &
               TRIM( L1BIrrCoreArch%RangeBeginningTime ) // "Z" 

        GROUPS(:) = ""
        status = PGS_MET_init( mcf_LUN, GROUPS )
        IF( status /= PGS_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "PGS_MET_init failed.", &
                                  "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! Read SHORTNAME which is set in the MCF file.
        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "SHORTNAME", ShortName )
        IF( status /= PGS_S_SUCCESS ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "PGS_MET_getSetAttr_s failed for SHORTNAME", &
                                  "L2_setCoreArchMetaData", zero )
            status = OZT_E_FAILURE
            RETURN
        ENDIF

        !! Read VERSIONID, a integer with no more than 3 digit, again 
        !! it is set in the MCF file
        status = PGS_MET_getSetAttr_i( GROUPS(INVENTORY), &
                                      "VERSIONID", VersionID )
        IF( status /= PGS_S_SUCCESS ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "PGS_MET_getSetAttr_s failed for VERSIONID", &
                                  "L2_setCoreArchMetaData", zero )
            status = OZT_E_FAILURE
            RETURN
        ENDIF
                         
        !! Check the make sure the Associated Sensor stuff in MCF are the same
        !! in L1B, generate warning message in Log file, but use the values
        !! in MCF for constructing Local granule ID.
        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "ASSOCIATEDPLATFORMSHORTNAME.1", & 
                                       MCF_AssociatedPlatformSN )
        IF( TRIM(L1BRadCoreArch%AssociatedPlatformSN) /= &
            TRIM(MCF_AssociatedPlatformSN) ) THEN
            WRITE( msg,'(A)' ) "AssociatedPlatformSN not matched,L1B:"//&
                  TRIM(L1BRadCoreArch%AssociatedPlatformSN)// " MCF:" // &
                  TRIM(MCF_AssociatedPlatformSN)
            ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                  "L2_setCoreArchMetaData", four )
        ENDIF

        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "ASSOCIATEDINSTRUMENTSHORTNAME.1", & 
                                       MCF_AssociatedInstrumentSN )
        IF( TRIM(L1BRadCoreArch%AssociatedInstrumentSN) /= &
            TRIM(MCF_AssociatedInstrumentSN) ) THEN
           WRITE( msg,'(A)' )"AssociatedInstrumentShortname not matched,L1B:"//&
                  TRIM(L1BRadCoreArch%AssociatedInstrumentSN)//" MCF:"// &
                  TRIM(MCF_AssociatedInstrumentSN)
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                  "L2_setCoreArchMetaData", four )
        ENDIF
        status = PGS_MET_getSetAttr_s( GROUPS(INVENTORY), &
                                      "ASSOCIATEDSENSORSHORTNAME.1", &
                                       MCF_AssociatedSensorSN )
        IF( TRIM(L1BRadCoreArch%AssociatedSensorSN) /=  &
            TRIM(MCF_AssociatedSensorSN) ) THEN
           WRITE( msg,'(A)' )"AssociatedSensorShortname not matched,L1B:"//&
                  TRIM(L1BRadCoreArch%AssociatedSensorSN )//" MCF:"// &
                  TRIM(MCF_AssociatedSensorSN )
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                  "L2_setCoreArchMetaData", four )
        ENDIF

        !! Construct the LocalGranuleID for L2 output
        InstrumentPlatformName =  TRIM(MCF_AssociatedInstrumentSN ) // &
                                  "-" // TRIM(MCF_AssociatedPlatformSN ) 
        status = OMI_localGranuleID( RangeBeginningDateTime, &
                                     InstrumentPlatformName, &
                                    "L2", ShortName,  orbitNumber, &
                                     VersionID, "he5", localGranuleID );
        IF( status /= OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "set local granule id failed.", &
                                  "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN 
        ENDIF

        !! set the local granule ID metadata
        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), "LOCALGRANULEID", &
                                    TRIM( localGranuleID )  )
        IF( status /= PGS_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "PGS_MET_setAttr_s failed for LOCALGRANULEID",&
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEBEGINNINGDATE", &
                                   L1BRadCoreArch%RangeBeginningDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "RANGEBEGINNINGDATE"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEBEGINNINGTIME", &
                                   L1BRadCoreArch%RangeBeginningTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "RANGEBEGINNINGTIME"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEENDINGDATE", & 
                                   L1BRadCoreArch%RangeEndingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                               "RANGEENDINGDATE"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "RANGEENDINGTIME", &
                                   L1BRadCoreArch%RangeEndingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "RANGEENDINGTIME"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "EQUATORCROSSINGDATE.1", &
                                   L1BRadCoreArch%EquatorCrossingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                               "EQUATORCROSSINGDATE.1"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                   "EQUATORCROSSINGTIME.1", &
                                   L1BRadCoreArch%EquatorCrossingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "EQUATORCROSSINGTIME.1"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_setAttr_i( GROUPS(INVENTORY), &
                                    "ORBITNUMBER.1", orbitNumber )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for ORBITNUMBER.1"
           ierr = OMI_SMF_setmsg( status, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! Write the Input Pointer Metadata Field
        status = PGS_MET_SetMultiAttr_s( GROUPS(INVENTORY),"InputPointer", &
                                        SIZE(LUNinputPointer), inputPointer )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_s failed for "// &
                              "InputPointer"
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        npL2 = SIZE( L2_Parameters )
        DO ii = 1, npL2
          WRITE( dummyName, '(A,I1)') "PARAMETERNAME.", ii
          status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                      TRIM( dummyName ), &
                                      L2_Parameters(ii)%Name )
          IF( status /= PGS_S_SUCCESS ) THEN
             WRITE( msg,'(A)' ) "PGS_MET_setAttr_s PARAMETERNAME failed: "// &
                                TRIM( L2_Parameters(ii)%Name )
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                   "L2_setCoreArchMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          IF( TRIM( ShortName ) == "OMTO3"    .OR. &
              TRIM( ShortName ) == "EPTO3" .OR. &
              TRIM( ShortName ) == "N7TO3" .OR. &
              TRIM( ShortName ) == "GMTO3"            ) THEN
            L2_Parameters(ii)%QAPercentHighQualityData = &
               NINT( (  L2_Parameters(ii)%NumberOfGoodOutputSamples            &
                      + L2_Parameters(ii)%NumberOfGlintCorrectedSamples)*100.0/&
                     (  L2_Parameters(ii)%NumberOfGoodInputSamples             &
                      - L2_Parameters(ii)%NumberOfLargeSZAInputSamples ) )
            L2_Parameters(ii)%QAPercentMissingData =                           &
               NINT( L2_Parameters(ii)%NumberOfMissingInputSamples*100.0/      &
                     L2_Parameters(ii)%NumberOfInputSamples )
          ENDIF

          IF( L2_Parameters(ii)%QAPercentHighQualityData >= 90  ) THEN
             dummyValue = "Passed"
          ELSE IF( L2_Parameters(ii)%QAPercentHighQualityData >= 60  ) THEN
             dummyValue = "Suspect"
          ELSE 
             dummyValue = "Failed"
          ENDIF

          IF( TRIM(dummyValue) == "Passed" .AND. &
             INDEX( L1BRadCoreArch%AutomaticQualityFlag, "Passed" ) == 0 ) THEN
             dummyValue = "Suspect"
          ENDIF

          WRITE( dummyName, '(A,I1)') "AUTOMATICQUALITYFLAG.", ii
          status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                      TRIM( dummyName ), dummyValue )
          IF( status /= PGS_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "PGS_MET_setAttr_s failed "//&
                                    TRIM(dummyName)//":"// TRIM(dummyValue), &
                                   "L2_setCoreArchMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          WRITE( dummyName, '(A,I1)') "AUTOMATICQUALITYFLAGEXPLANATION.", ii
          status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                      TRIM( dummyName ), &
                                      TRIM( AutomaticQualityFlagExplanation ) )
          IF( status /= PGS_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "PGS_MET_setAttr_s failed "//&
                                    TRIM( dummyName ), &
                                   "L2_setCoreArchMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          WRITE( dummyName, '(A,I1)') "QAPERCENTMISSINGDATA.", ii
          status = PGS_MET_setAttr_i( GROUPS(INVENTORY), TRIM( dummyName), &
                                      L2_Parameters(ii)%QAPercentMissingData )

          IF( status /= PGS_S_SUCCESS ) THEN
             WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 TRIM( dummyName )
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                   "L2_setCoreArchMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF
        ENDDO

        !! Now set the 1 (double type) core metadata fields (in both
        !! L1B  and L2 ) that is supposed to be set by the PGE. 
        !! Write the EQUATORCROSSINGLONGITUDE Metadata Field

        status = PGS_MET_setAttr_d( GROUPS(INVENTORY), &
                                   "EQUATORCROSSINGLONGITUDE.1", &
                                   L1BRadCoreArch%EqCrossLon )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_d failed for "// &
                              "EQUATORCROSSINGLONGITUDE.1"
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! copy some PSA from L1B.
        
        counter_w = 1  !! counter for writing of L2 PSA
        DO ii = 1, L1BRadCoreArch%nL1Bpsa
          IF( INDEX( TRIM(PSAtoBeCopiedFromL1B), & 
              TRIM(L1BRadCoreArch%L1BpsaNames(ii)) ) == 0 ) CYCLE
          IF( counter_w < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "ADDITIONALATTRIBUTENAME.",counter_w
          ELSE IF( counter_w < 100 ) THEN
             WRITE( dummyName, '(A,I2)') "ADDITIONALATTRIBUTENAME.",counter_w
          ELSE
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                   "Too many PSA in input file",&
                                  "L2_setCoreArchMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                      TRIM( dummyName ), &
                                      TRIM(L1BRadCoreArch%L1BpsaNames(ii)) )

          !! figure out the number of strings in the metadata, and
          !! write the exact same number of stings in the L2 metadata.
          nStr = 0
          !Allows counter to go above array limit, so replace
!          DO WHILE( LEN_TRIM(L1BRadCoreArch%L1BpsaValues(ii,nStr+1))  > 0 )
!            nStr = nStr + 1
!          ENDDO
          do jj=1,20
            if (LEN_TRIM(L1BRadCoreArch%L1BpsaValues(ii,jj)) > 0) then
              nStr = nStr + 1
            endif
          enddo

          IF( counter_w < 10 ) THEN
            WRITE( dummyName, '(A,I1)') "PARAMETERVALUE.", counter_w
          ELSE IF( counter_w < 100 ) THEN
            WRITE( dummyName, '(A,I2)') "PARAMETERVALUE.", counter_w
          ENDIF

          temp_L1BpsaValues(:) = L1BRadCoreArch%L1BpsaValues(ii,:) 
          status = PGS_MET_SetMultiAttr_s( GROUPS(INVENTORY), &
                                           TRIM( dummyName ), &
                                           nStr, &
!                                           L1BRadCoreArch%L1BpsaValues(ii,:) )
                                           temp_L1BpsaValues(:))
          counter_w = counter_w + 1 
        END DO 

        status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                   "QAPercentHighQualityData", &
                                    L2_Parameters(1)%QAPercentHighQualityData )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                              "QAPercentHighQualityData"
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                 "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF
        IF( TRIM( ShortName ) == "EPTO3" .OR. &
            TRIM( ShortName ) == "N7TO3" ) THEN
           nVal = 6
           DO ii = 1, L1BRadCoreArch%nL1BGring
             IF( ii < 10 ) THEN
                WRITE( dummyName, '(A,I1)') "GRINGPOINTLONGITUDE.", ii
             ELSE IF( ii < 100 ) THEN
                WRITE( dummyName, '(A,I2)') "GRINGPOINTLONGITUDE.", ii
             ELSE
                ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                   "Too many GRings in output file",&
                                   "L2_setCoreArchMetaData", zero )
                status = OZT_E_FAILURE
                RETURN
             ENDIF
              
             status = PGS_MET_SetMultiAttr_d( GROUPS(INVENTORY), &
                                              TRIM( dummyName ), &
                                              nVal, &
                                      L1BRadCoreArch%GRingPointLongitude(ii,:) )
             IF( status /= PGS_S_SUCCESS ) THEN
                ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Set GRingLon failed", &
                                       "L2_setCoreArchMetaData", zero )
                status = OZT_E_FAILURE
                RETURN
             ENDIF
 
             IF( ii < 10 ) THEN
                WRITE( dummyName, '(A,I1)') "GRINGPOINTLATITUDE.", ii
             ELSE IF( ii < 100 ) THEN
                WRITE( dummyName, '(A,I2)') "GRINGPOINTLATITUDE.", ii
             ENDIF
 
             status = PGS_MET_SetMultiAttr_d( GROUPS(INVENTORY), &
                                              TRIM( dummyName ), &
                                              nVal, &
                                       L1BRadCoreArch%GRingPointLatitude(ii,:) )
             IF( status /= PGS_S_SUCCESS ) THEN
                ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Set GRingLat failed", &
                                       "L2_setCoreArchMetaData", zero )
                status = OZT_E_FAILURE
                RETURN
             ENDIF
 
             IF( ii < 10 ) THEN
                WRITE( dummyName, '(A,I1)') "GRINGPOINTSEQUENCENO.", ii
             ELSE IF( ii < 100 ) THEN
                WRITE( dummyName, '(A,I2)') "GRINGPOINTSEQUENCENO.", ii
             ENDIF
 
             status = PGS_MET_SetMultiAttr_d( GROUPS(INVENTORY), &
                                              TRIM( dummyName ), &
                                              nVal, &
                                     L1BRadCoreArch%GRingPointSequenceNo(ii,:) )
             IF( status /= PGS_S_SUCCESS ) THEN
                ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Set GRingSeq failed", &
                                       "L2_setCoreArchMetaData", zero )
                status = OZT_E_FAILURE
                RETURN
             ENDIF
  
             IF( ii < 10 ) THEN
                WRITE( dummyName, '(A,I1)') "EXCLUSIONGRINGFLAG.", ii
             ELSE IF( ii < 100 ) THEN
                WRITE( dummyName, '(A,I2)') "EXCLUSIONGRINGFLAG.", ii
             ENDIF
             status = PGS_MET_setAttr_s( GROUPS(INVENTORY), &
                                         TRIM( dummyName ), &
                                         L1BRadCoreArch%ExclusionGRingFlag(ii) )
             IF( status /= PGS_S_SUCCESS ) THEN
                ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Set GRingExl failed "//&
                                       TRIM( dummyName ), &
                                      "L2_setCoreArchMetaData", zero )
                status = OZT_E_FAILURE
                RETURN
             ENDIF
           ENDDO
        ENDIF

        !! set some achive metadata
        IF( TRIM( ShortName ) == "OMTO3"    .OR. &
            TRIM( ShortName ) == "EPTO3" .OR. &
            TRIM( ShortName ) == "N7TO3" .OR. &
            TRIM( ShortName ) == "GMTO3"            ) THEN
           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), "NumberOfInputSamples",&
                                       L2_Parameters(1)%NumberOfInputSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfInputSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfGoodInputSamples", &
                                      L2_Parameters(1)%NumberOfGoodInputSamples)
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfGoodInputSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfLargeSZAInputSamples", &
                                L2_Parameters(1)%NumberOfLargeSZAInputSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfLargeSZAInputSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfMissingInputSamples", &
                                 L2_Parameters(1)%NumberOfMissingInputSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfMissingInputSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfBadInputSamples", &
                                      L2_Parameters(1)%NumberOfBadInputSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfBadInputSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfInputWarningSamples", &
                                   L2_Parameters(1)%NumberOfInputWarningSamples)
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfInputWarningSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfGoodOutputSamples", &
                                     L2_Parameters(1)%NumberOfGoodOutputSamples)
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfGoodOutputSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfGlintCorrectedSamples", &
                                L2_Parameters(1)%NumberOfGlintCorrectedSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfGlintCorrectedSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfSkippedSamples", &
                                      L2_Parameters(1)%NumberOfSkippedSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfSkippedSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfStep1InvalidSamples", &
                                  L2_Parameters(1)%NumberOfStep1InvalidSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfStep1InvalidSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "NumberOfStep2InvalidSamples", &
                                  L2_Parameters(1)%NumberOfStep2InvalidSamples )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A)' ) "PGS_MET_setAttr_i failed for "// &
                                 "NumberOfStep2InvalidSamples"
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "L2_setCoreArchMetaData", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           nqfc = SIZE( L2_Parameters(1)%QualityFlagsCounters(1,:) )
           temp_QFCounters(:) = L2_Parameters(1)%QualityFlagsCounters(1,:) 
           status = PGS_MET_SetMultiAttr_i( GROUPS(ARCHIVE), &
                      "AlgorithmPath1QualityFlagCounts", &
!                      nqfc, L2_Parameters(1)%QualityFlagsCounters(1,:) )
                      nqfc, temp_QFCounters(:) )

           temp_QFCounters(:)=L2_Parameters(1)%QualityFlagsCounters(2,:)
           status = PGS_MET_SetMultiAttr_i( GROUPS(ARCHIVE), &
                      "AlgorithmPath2QualityFlagCounts", &
!                      nqfc, L2_Parameters(1)%QualityFlagsCounters(2,:) )
                      nqfc, temp_QFCounters(:) )

           temp_QFCounters(:) = L2_Parameters(1)%QualityFlagsCounters(3,:)
           status = PGS_MET_SetMultiAttr_i( GROUPS(ARCHIVE), &
                      "AlgorithmPath3QualityFlagCounts", &
!                      nqfc, L2_Parameters(1)%QualityFlagsCounters(3,:) )
                      nqfc, temp_QFCounters(:) )

           temp_ZLatRange(:) = L2_Parameters(1)%ZonalLatRange(1,:)
           status = PGS_MET_SetMultiAttr_d( GROUPS(ARCHIVE), &
                      "ZonalLatitudeMinimum", &
!                      nZones, L2_Parameters(1)%ZonalLatRange(1,:) )
                      nZones, temp_ZLatRange(:) )

           temp_ZLatRange(:) = L2_Parameters(1)%ZonalLatRange(2,:)
           status = PGS_MET_SetMultiAttr_d( GROUPS(ARCHIVE), &
                      "ZonalLatitudeMaximum", &
!                      nZones, L2_Parameters(1)%ZonalLatRange(2,:) )
                      nZones, temp_ZLatRange(:) )

           ! rounding
           roundZOmin=ANINT(L2_Parameters(1)%ZonalOzoneMin*100.0)
           tempZOmin=real(roundZOmin, kind=8)/100.0d0
           status = PGS_MET_SetMultiAttr_d( GROUPS(ARCHIVE), &
                      "ZonalOzoneMinimum", &
!                      nZones, L2_Parameters(1)%ZonalOzoneMin )
                      nZones, tempZOmin)

           roundZOmax=ANINT(L2_Parameters(1)%ZonalOzoneMax*100.0)
           tempZOmax=real(roundZOmax, kind=8)/100.0d0
           status = PGS_MET_SetMultiAttr_d( GROUPS(ARCHIVE), &
                      "ZonalOzoneMaximum", &
!                      nZones, L2_Parameters(1)%ZonalOzoneMax )
                      nZones, tempZOmax)

           IF( INDEX( TRIM(ShortName), "TOMS" ) == 0 ) THEN
              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfIrradianceMissing/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctIrradianceMissing", I4foo )
                                    
              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfIrradianceError/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctIrradianceError", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfIrradianceWarning/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                      "QAPctIrradianceWarning", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfRadianceMissing/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctRadianceMissing", I4foo )
                                    
              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfRadianceError/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctRadianceError", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfRadianceWarning/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctRadianceWarning", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfMeasurementMissing/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctMeasurementMissing", I4foo )
                                    
              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfMeasurementError/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctMeasurementError", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfMeasurementWarning/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctMeasurementWarning", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfMeasurementRebinned/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctMeasurementRebinned", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfMeasurementSAA/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctMeasurementSAA", I4foo )

              I4foo = NINT( 100.0*L2_Parameters(1)%NumberOfMeasurementManeuver/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctMeasurementManeuver", I4foo )

              IF( L2_Parameters(1)%NumberOfInstrumentSettingsError > 0 ) THEN
                 status = PGS_MET_setAttr_s( GROUPS(ARCHIVE), &
                           "INSTRUMENTCONFIGURATIONIDSNONCOMPATIBLE", &
                           "TRUE" )
              ENDIF
              
              I4foo = NINT(100.0*&
                           L2_Parameters(1)%NumberOfInstrumentSettingsError/&
                                  L2_Parameters(1)%NumberOfMeasurement )
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "QAPctInstrumentSettingsError", I4foo )

              I4foo = 0
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "SOLARPRODUCTMISSING", &
                                  L2_Parameters(1)%SolarProductMissing )
              status = PGS_TD_UTCtoTAI( L1BRadDateTime, L1BRadTAI93 )
              status = PGS_TD_UTCtoTAI( L1BIrrDateTime, L1BIrrTAI93 )
              IF( ABS( L1BRadTAI93 - L1BIrrTAI93 ) > 48*60*60 ) THEN 
                 L2_Parameters(1)%SolarProductOutOfDate = 1
                 WRITE( msg,'(A,F8.2,A)' ) "Solar Irradiance out of date by ", &
                                  (L1BRadTAI93 - L1BIrrTAI93)/(24*3600.0), &
                                   " day(s)." 
                 ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                       "L2_setCoreArchMetaData", zero )
              ELSE
                 L2_Parameters(1)%SolarProductOutOfDate = 0
              ENDIF

              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "SOLARPRODUCTOUTOFDATE", &
                                      L2_Parameters(1)%SolarProductOutOfDate )
              I4foo = 0
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "BACKUPSOLARPRODUCTUSED", &
                                  L2_Parameters(1)%BackupSolarProductUsed )

              IF( L2_Parameters(1)%NumberOfIrradianceWarning > 0 ) THEN
                 I4foo = 1
              ELSE
                 I4foo = 0
              ENDIF
              status = PGS_MET_setAttr_i( GROUPS(ARCHIVE), &
                                         "SOLARIRRADIANCEWARNING", I4foo )

              status = PGS_MET_setAttr_s( GROUPS(ARCHIVE), &
                                         "RADIANCESCIENCEQUALITYFLAG", &
                                 TRIM( L1BRadCoreArch%ScienceQualityFlag ) )
              IF( status /= PGS_S_SUCCESS ) THEN
                 ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                        "PGS_MET_setAttr_s failed "//&
                                        "RADIANCESCIENCEQUALITYFLAG", &
                                        "L2_setCoreArchMetaData", zero )
                 status = OZT_E_FAILURE
                 RETURN
              ENDIF

              status = PGS_MET_setAttr_s( GROUPS(ARCHIVE), &
                                         "IRRADIANCESCEINCEQUALITYFLAG", &
                                 TRIM( L1BIrrCoreArch%ScienceQualityFlag ) )
              IF( status /= PGS_S_SUCCESS ) THEN
                 ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                        "PGS_MET_setAttr_s failed "//&
                                        "IrradianceSceinceQualityFlag", &
                                        "L2_setCoreArchMetaData", zero )
                 status = OZT_E_FAILURE
                 RETURN
              ENDIF
           ENDIF  !!end IF( INDEX( TRIM(ShortName), "TOMS" ) == 0 ) THEN
        ENDIF !! end IF( TRIM( ShortName ) == "OMTO3"    .OR. &
              !!         TRIM( ShortName ) == "EPTO3" .OR. &
              !!         TRIM( ShortName ) == "N7TO3" .OR. &
              !!         TRIM( ShortName ) == "GMTO3"            ) THEN

        !! Now all the Metadata are set, write them to the output file.
        !! retrieve the L2 output file name
        version = 1
        status = PGS_PC_getReference( L2_outFile_LUN, version, L2_filename )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I9)' ) "get filename from PCF file at LUN =", &
                                 L2_outFile_LUN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                  "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! Initiate it for writing metadata 
        status = PGS_MET_sfstart( TRIM(L2_filename), HDF5_ACC_RDWR, HE5id )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_sfstart failed " // &
                               TRIM(L2_filename)
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                  "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! Write the CoreMetadata
        status = PGS_MET_write( GROUPS(INVENTORY),'CoreMetadata', HE5id )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_write CoreMetadata failed. " 
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                  "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! Write the Archived metadata
        status = PGS_MET_write( GROUPS(ARCHIVE),'ArchivedMetadata', HE5id )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "PGS_MET_write ArchivedMetadata failed. " 
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                  "L2_setCoreArchMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF
        
        !! Close the L2 file 
        status = PGS_MET_SFend(HE5id)
        !! Clean up the memory associated writing metadata
        status = PGS_MET_Remove()

        status = OZT_S_SUCCESS
      END FUNCTION L2_setCoreArchMetaData

END MODULE L2_metaData_class
