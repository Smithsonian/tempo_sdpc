!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE L1B_metaData_class
! 
!  contains functions to read the ECS core and archived metadata in the 
!  L1B file.
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
MODULE L1B_metaData_class
    USE OMI_SMF_class
    USE PGS_MET_class
    USE UTIL_tools_class
    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, four = 4

    CHARACTER( LEN = 80 ):: inventoryMetadataName
    CHARACTER( LEN = 80 ):: archivedMetadataName
    LOGICAL (KIND=4) :: L1BMETADATA_READ = .FALSE.

    PUBLIC :: L1B_getCoreArchivedMetaData
    PUBLIC :: L1B_getGlobalAttrname

    CONTAINS
      FUNCTION L1B_getCoreArchivedMetaData( l1b_file_lun, &
                                            L1BcoreArch,  &
                                            year, month, day, jday, &
                                            version_in ) &
                                            RESULT( status ) 
        INTEGER (KIND=4), INTENT(IN) :: l1b_file_lun
        TYPE (L1BECSMETA_T), INTENT(OUT) :: L1BcoreArch
        INTEGER (KIND=4), INTENT(OUT) :: year, month, day, jday
        INTEGER (KIND=4), INTENT(IN), OPTIONAL :: version_in
        INTEGER (KIND=4) :: status ! iName, 
        INTEGER (KIND=4) :: ierr, version
        CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: dummyName!, dummyValue
        INTEGER (KIND=4) :: counter_r
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
        CHARACTER( LEN = 1 ) :: foo
        INTEGER (KIND=4), EXTERNAL :: day_of_year

        L1BMETADATA_READ = .FALSE.

        !! CoreMetadata in L1B can sometimes written in several
        !! variations, such as "coremetadata", or "CoreMetadata.0"
        !! Here it figure out exactly how it was written in the L1B,
        !! then use the name to read the fields in the CoreMetadata.

        status = L1B_getGlobalAttrname( l1b_file_lun, &
                                        inventoryMetadataName, "core" ) 
        IF( status /= OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, "can't find coremetadata "//&
                                  " or its variation in input L1B file.",&
                                  "L1B_getCoreArchivedMetaData", zero )  
           inventoryMetadataName = "CoreMetadata.0"
        ENDIF

        version = 1
        IF( PRESENT( version_in) ) version = version_in
           
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEBEGINNINGDATE",             &
                                      L1BcoreArch%RangeBeginningDate )
        IF( status == PGS_S_SUCCESS ) THEN
           READ( L1BcoreArch%RangeBeginningDate , '(I4,A1,I2,A1,I2)' ) &
                 year, foo, month, foo, day
           jday = day_of_year( year, month, day )
           WRITE( msg, '(A, I5,I3,I3,I4)' ) "(y,m,d) = ", &
                       year, month, day, jday 
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                  "L1B_getCoreArchivedMetaData",four )
        ELSE
           L1BcoreArch%RangeBeginningDate = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "can't get RangeBeginningDate" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEBEGINNINGTIME",             &
                                      L1BcoreArch%RangeBeginningTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%RangeBeginningTime = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "can't get RangeBeginningTime" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEENDINGDATE",                &
                                      L1BcoreArch%RangeEndingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%RangeEndingDate = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "can't get RangeEndingDate" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version,            &
                                      inventoryMetadataName,            & 
                                      "RANGEENDINGTIME",                &
                                      L1BcoreArch%RangeEndingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%RangeEndingTime = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "can't get RangeEndingTime" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version,            &
                                      inventoryMetadataName,            & 
                                      "EQUATORCROSSINGTIME.1",          &
                                      L1BcoreArch%EquatorCrossingTime )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%EquatorCrossingTime = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "can't get EquatorCrossingTime" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version,            &
                                      inventoryMetadataName,            & 
                                      "EQUATORCROSSINGDATE.1",          &
                                      L1BcoreArch%EquatorCrossingDate )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%EquatorCrossingDate = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "can't get EquatorCrossingDate" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version,            &
                                      inventoryMetadataName,            & 
                                      "SHORTNAME",                      &
                                      L1BcoreArch%ShortName )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%ShortName = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "can't get ShortName" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                      inventoryMetadataName,& 
                                     "ASSOCIATEDPLATFORMSHORTNAME.1", & 
                                      L1BcoreArch%AssociatedPlatformSN  )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%AssociatedPlatformSN = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "can't get ASSOCIATEDPLATFORMSHORTNAME" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                      inventoryMetadataName,& 
                                     "ASSOCIATEDINSTRUMENTSHORTNAME.1", &
                                      L1BcoreArch%AssociatedInstrumentSN  )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%AssociatedInstrumentSN = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "can't get ASSOCIATEDINSTRUMENTSHORTNAME" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
        ENDIF

        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                      inventoryMetadataName,& 
                                     "ASSOCIATEDSENSORSHORTNAME.1", &
                                      L1BcoreArch%AssociatedSensorSN  )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%AssociatedSensorSN = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "can't get ASSOCIATEDSENSORSHORTNAME" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
        ENDIF

        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                      inventoryMetadataName,& 
                                     "AUTOMATICQUALITYFLAG.1", &
                                      L1BcoreArch%AutomaticQualityFlag )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%AutomaticQualityFlag = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "can't get AUTOMATICQUALITYFLAG" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                      inventoryMetadataName,& 
                                     "SCIENCEQUALITYFLAG.1", &
                                      L1BcoreArch%ScienceQualityFlag )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%AutomaticQualityFlag = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "can't get SCIENCEQUALITYFLAG" &
                                  // " from L1B file core Meta data.", &
                                  "L1B_getCoreArchivedMetaData", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! Read Integer Type of L1B data fileds.
        ! 1. Get OrbitNumber
        version = 1
        status = PGS_MET_GetPCAttr_i( l1b_file_lun, version , &
                                      inventoryMetadataName,& 
                                     "ORBITNUMBER.1", &  
                                      L1BcoreArch%orbitNumber )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%orbitNumber = -1
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "get OrbitNumber.1 failed.", &
                                 "L1B_getCoreArchivedMetaData", zero )  
           status = OZT_E_FAILURE
           return
        ENDIF

        ! 2. Get QAPercentMissingData
        version = 1
        status = PGS_MET_GetPCAttr_i( l1b_file_lun, version , &
                                      inventoryMetadataName,& 
                                     "QAPERCENTMISSINGDATA.1", &
                                      L1BcoreArch%QAPercentMissingData )
        IF( status /= PGS_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "get QAPercentMissingData.1 failed.", &
                                 "L1B_getCoreArchivedMetaData", zero )  
        ENDIF

        !! Read Double Type of L1B data fileds (there is just one of them).
        version = 1
        status = PGS_MET_GetPCAttr_d( l1b_file_lun, version, &
                                      inventoryMetadataName,& 
                                     "EQUATORCROSSINGLONGITUDE.1", &
                                      L1BcoreArch%EqCrossLon )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%EqCrossLon = -999.0
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                     "get EQUATORCROSSINGLONGITUDE.1 failed.", &
                     "L1B_getCoreArchivedMetaData", zero )  
        ELSE          ! write to the Log files the field values of L1B file
           WRITE(msg, '(A,F9.5)' ) "EQUATORCROSSINGLONGITUDE.1 = ", &
                                   L1BcoreArch%EqCrossLon
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                 "L1B_getCoreArchivedMetaData",four )
        ENDIF

        counter_r = 1  !! counter for reading of L1B PSA
        version = 1
        WRITE( dummyName, '(A,I1)') "ADDITIONALATTRIBUTENAME.", counter_r
        status = PGS_MET_getPCAttr_s( l1b_file_lun, version,               &
                                      inventoryMetadataName,               &
                                      dummyName,                           &
                                      L1BcoreArch%L1BpsaNames(counter_r) )
        DO WHILE( status == PGS_S_SUCCESS )
          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "PARAMETERVALUE.", counter_r
          ELSE IF( counter_r < 100 ) THEN
             WRITE( dummyName, '(A,I2)') "PARAMETERVALUE.", counter_r
          ELSE
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "Too many PSA in input file", &
                                    "L1B_getCoreArchivedMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          !! empty string Array before it is read in from L1B metadata
          L1BcoreArch%L1BpsaValues(counter_r,:) = "";
          version = 1
          status = PGS_MET_getPCAttr_s( l1b_file_lun, version,  &
                                        inventoryMetadataName,  &
                                        dummyName,              &
                                        L1BcoreArch%L1BpsaValues(counter_r,:) )
          IF( status /= PGS_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "GEt PSA value failed", &
                                   "L1B_getCoreArchivedMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          L1BcoreArch%nL1Bpsa   = counter_r
          counter_r = counter_r + 1;

          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "ADDITIONALATTRIBUTENAME.", counter_r
          ELSE IF( counter_r < 100 ) THEN
             WRITE( dummyName, '(A,I2)') "ADDITIONALATTRIBUTENAME.", counter_r
          ELSE
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Too many PSA in input file",&
                                    "O3T_getL1BmetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF
          status = PGS_MET_getPCAttr_s( l1b_file_lun, version,               &
                                        inventoryMetadataName,               &
                                        dummyName,                           &
                                        L1BcoreArch%L1BpsaNames(counter_r) )
        END DO


        L1BcoreArch%nL1BGring = 0
        counter_r = 1  !! counter for reading of L1B Gring
        version = 1
        WRITE( dummyName, '(A,I1)') "GRINGPOINTLONGITUDE.", counter_r
        status = PGS_MET_getPCAttr_d( l1b_file_lun, version,               &
                                      inventoryMetadataName,               &
                                      dummyName,                           &
                                L1BcoreArch%GRingPointLongitude(counter_r,:) )
        DO WHILE( status == PGS_S_SUCCESS )
          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "GRINGPOINTLATITUDE.", counter_r
          ELSE IF( counter_r < NRing_MAX ) THEN
             WRITE( dummyName, '(A,I2)') "GRINGPOINTLATITUDE.", counter_r
          ELSE
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "Too many GRings in input file", &
                                    "L1B_getCoreArchivedMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          version = 1
          status = PGS_MET_getPCAttr_d( l1b_file_lun, version,  &
                                        inventoryMetadataName,  &
                                        dummyName,              &
                                 L1BcoreArch%GRingPointLatitude(counter_r,:) )
          IF( status /= PGS_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Get GRingLat value failed", &
                                   "L1B_getCoreArchivedMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "GRINGPOINTSEQUENCENO.", counter_r
          ELSE IF( counter_r < NRing_MAX ) THEN
             WRITE( dummyName, '(A,I2)') "GRINGPOINTSEQUENCENO.", counter_r
          ENDIF

          version = 1
          status = PGS_MET_getPCAttr_i( l1b_file_lun, version,  &
                                        inventoryMetadataName,  &
                                        dummyName,              &
                                 L1BcoreArch%GringPointSequenceNo(counter_r,:) )
          IF( status /= PGS_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Get GRingSeq value failed", &
                                   "L1B_getCoreArchivedMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF
          
          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "EXCLUSIONGRINGFLAG.", counter_r
          ELSE IF( counter_r < NRing_MAX ) THEN
             WRITE( dummyName, '(A,I2)') "EXCLUSIONGRINGFLAG.", counter_r
          ENDIF

          version = 1
          status = PGS_MET_getPCAttr_s( l1b_file_lun, version,  &
                                        inventoryMetadataName,  &
                                        dummyName,              &
                                 L1BcoreArch%ExclusionGRingFlag(counter_r) )
          IF( status /= PGS_S_SUCCESS ) THEN
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Get GRingExc value failed", &
                                   "L1B_getCoreArchivedMetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          L1BcoreArch%nL1BGring   = counter_r
          counter_r = counter_r + 1;

          IF( counter_r < 10 ) THEN
             WRITE( dummyName, '(A,I1)') "GRINGPOINTLONGITUDE.", counter_r
          ELSE IF( counter_r < NRing_MAX ) THEN
             WRITE( dummyName, '(A,I2)') "GRINGPOINTLONGITUDE.", counter_r
          ELSE
             ierr = OMI_SMF_setmsg( OZT_E_INPUT,"Too many Gring in input file",&
                                    "O3T_getL1BmetaData", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF
          status = PGS_MET_getPCAttr_d( l1b_file_lun, version,               &
                                        inventoryMetadataName,               &
                                        dummyName,                           &
                                L1BcoreArch%GRingPointLongitude(counter_r,:) )
        END DO


        status = L1B_getGlobalAttrname( l1b_file_lun, &
                                        archivedMetadataName, "arch" ) 
        IF( status /= OZT_S_SUCCESS ) THEN
           WRITE( msg, '(A,I9, A)' ) "can't find archivemetadata"// &
                    " or its variation in input L1B file at lun = ", &
                    l1b_file_lun, ". Skip reading it."
           ierr=OMI_SMF_setmsg( OZT_W_GENERAL, msg,  &
                               "L1B_getCoreArchivedMetaData", zero )  
           L1BMETADATA_READ = .TRUE.
           status = OZT_S_SUCCESS
           RETURN
        ENDIF
        version = 1
        status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                      archivedMetadataName,   &
                                      "LongName",  &
                                      L1BcoreArch%LongName )
        IF( status /= PGS_S_SUCCESS ) THEN
           L1BcoreArch%LongName = ""
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "can't get LongName" &
                                  // " from L1B file Archived Meta data.", &
                                  "L1B_getCoreArchivedMetaData", four )
        ENDIF

        IF( INDEX( L1BcoreArch%ShortName, "OML1B" ) >  0 ) THEN
          !! Read String Type of L1B data fileds.
          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "AlgorithmBypassList",  &
                                         L1BcoreArch%AlgorithmBypassList )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%AlgorithmBypassList = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get AlgorithmBypassList" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF

          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "ProcessingMode",  &
                                        L1BcoreArch%ProcessingMode )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%ProcessingMode = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get ProcessingMode" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF

          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "OrbitData",  &
                                        L1BcoreArch%OrbitData )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%OrbitData = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get OrbitData" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF
          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "ProcessingCenter",  &
                                        L1BcoreArch%ProcessingCenter )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%ProcessingCenter = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get ProcessingCenter" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF
          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "ESDTDescriptorRevision",  &
                                        L1BcoreArch%ESDTDescriptorRevision )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%ESDTDescriptorRevision = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get ESDTDescriptorRevision" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF
          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "OPFSmearSwitchValue",  &
                                        L1BcoreArch%OPFSmearSwitchValue )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%OPFSmearSwitchValue = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get OPFSmearSwitchValue" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF
          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "OPFMeasurementStrayFlag",  &
                                        L1BcoreArch%OPFMeasurementStrayFlag )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%OPFMeasurementStrayFlag = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get OPFMeasurementStrayFlag" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF
          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "OPFVersion",  &
                                        L1BcoreArch%OPFVersion )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%OPFVersion = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get OPFVersion" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF
          version = 1
          status = PGS_MET_GetPCAttr_s( l1b_file_lun, version , &
                                        archivedMetadataName,   &
                                        "OPFValid",  &
                                        L1BcoreArch%OPFValid )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%OPFValid = ""
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "can't get OPFValid" &
                                    // " from L1B file Archived Meta data.", &
                                    "L1B_getCoreArchivedMetaData", four )
          ENDIF

          !! Read Real Type of Archived data fileds (there are two of them).
          version = 1
          status = PGS_MET_GetPCAttr_d( l1b_file_lun, version, &
                                        archivedMetadataName,&
                                       "SpacecraftMinAltitude", &
                                        L1BcoreArch%SpacecraftMinAltitude )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%SpacecraftMinAltitude = -1.0
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                       "get SpacecraftMinAltitude failed.", &
                       "L1B_getCoreArchivedMetaData", four )
          ELSE          ! write to the Log files the field values of L1B file
             WRITE(msg, '(A,E12.5)' ) "SpacecraftMinAltitude = ", &
                                      L1BcoreArch%SpacecraftMinAltitude
             ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                   "L1B_getCoreArchivedMetaData",four )
          ENDIF

          version = 1
          status = PGS_MET_GetPCAttr_d( l1b_file_lun, version, &
                                        archivedMetadataName,&
                                       "SpacecraftMaxAltitude", &
                                        L1BcoreArch%SpacecraftMaxAltitude )
          IF( status /= PGS_S_SUCCESS ) THEN
             L1BcoreArch%SpacecraftMaxAltitude = -1.0;
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                       "get SpacecraftMaxAltitude failed.", &
                       "L1B_getCoreArchivedMetaData", four )
          ELSE          ! write to the Log files the field values of L1B file
             WRITE(msg, '(A,E12.5)' ) "SpacecraftMaxAltitude = ", &
                                      L1BcoreArch%SpacecraftMaxAltitude
             ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, &
                                   "L1B_getCoreArchivedMetaData",four )
          ENDIF
        ENDIF

        L1BMETADATA_READ = .TRUE.
        status = OZT_S_SUCCESS
      END FUNCTION L1B_getCoreArchivedMetaData

      FUNCTION L1B_getGlobalAttrname( l1b_file_lun, attr_name, sub_string ) &
                                      RESULT( status ) 
        USE PGS_PC_class
        USE HDF4_class
        INTEGER (KIND=4), INTENT(IN) :: l1b_file_lun
        CHARACTER( LEN = * ), INTENT(IN), OPTIONAL :: sub_string
        CHARACTER( LEN = * ), INTENT(OUT)  :: attr_name
        INTEGER (KIND=4) :: ierr, version, status
        INTEGER (KIND=4) :: sd_id
        INTEGER (KIND=4) :: attr_index, data_type, count
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: sub_low
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: att_low
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
        CHARACTER( LEN = PGSd_PC_PATH_LENGTH_MAX ) :: L1BPathAndName

        version = 1
        status = PGS_PC_GetReference( l1b_file_lun, version, L1BPathAndName )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I9)' ) "get L1B filename failed at LUN = ", &
                                l1b_file_lun 
           ierr = OMI_SMF_setmsg( status, msg, "L1B_getGlobalAttrname",zero )
           status = OZT_E_FAILURE
           RETURN
        ELSE
           sd_id = sfstart( L1BPathAndName, DFACC_READ )  
           IF( sd_id == -1 ) THEN
              ierr = OMI_SMF_setmsg( status, "open " // &
                                     TRIM(L1BPathAndName) // "failed." , &
                                    "L1B_getGlobalAttrname", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
          
           attr_index = 0
           ierr = sfgainfo( sd_id, attr_index, attr_name, data_type, count )
           IF( PRESENT( sub_string ) ) THEN
              DO WHILE( ierr /= - 1) 
                att_low = toLowerString(attr_name)
                sub_low = toLowerString(sub_string)
                IF( INDEX( TRIM(att_low), TRIM(sub_low) ) > 0 ) THEN
                  ierr = sfend( sd_id )
                  status = OZT_S_SUCCESS
                  RETURN
                ENDIF
                attr_index = attr_index + 1
                ierr = sfgainfo( sd_id, attr_index, attr_name, data_type,count )
              ENDDO

              ierr = sfend( sd_id )
              attr_name = ""
              status = OZT_E_FAILURE
           ELSE 
              ierr = sfend( sd_id )
              IF( ierr /= -1 ) THEN
                 status = OZT_S_SUCCESS
                 RETURN
              ELSE
                 attr_name = ""
                 status = OZT_E_FAILURE
                 RETURN
              ENDIF
           ENDIF
        ENDIF

      END FUNCTION L1B_getGlobalAttrname

END MODULE L1B_metaData_class
