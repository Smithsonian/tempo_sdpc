!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE L2_attr_class
! 
!  contsins function the write the global attributes for the L2 output.
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
MODULE L2_attr_class
    USE L1B_metaData_class, ONLY :  L1B_getGlobalAttrname, &
                                    inventoryMetadataName, &
                                    archivedMetadataName
    USE OMI_SMF_class ! include PGE specific messages and OMI_SMF_setmsg
    USE PGS_MET_class
    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four =4, five = 5

    PUBLIC  :: OMI_writeGlobalAttribute
    PUBLIC  :: OMI_writeSwathAttribute
!    PRIVATE :: PGEVersion2PhaseScience

    CONTAINS

       FUNCTION OMI_writeGlobalAttribute( L2_filename, &
                                          GranuleYear, &
                                          GranuleMonth, &
                                          GranuleDay, &
                                          inputLUNs, &
                                          PSANumSampOutOfRange ) & 
                            RESULT( status )
         USE HE5_class
         USE OMI_LUN_set
         USE PGS_PC_class  ! define PGSd_PC_FILE_PATH_MAX, and PGS_PC functions
                           ! include PGS_SMF.f define PGS_SMF_MAX_MSG_SIZE
         USE PGS_TD_class
         USE UTIL_tools_class
         USE ISO_C_BINDING, ONLY: C_LONG
         CHARACTER( LEN = *), INTENT(IN) :: L2_filename
         INTEGER (KIND=4), DIMENSION(:), INTENT(IN) :: inputLUNs
         INTEGER (KIND=4), INTENT(IN) :: GranuleDay, GranuleMonth,  &
                                         GranuleYear  !!PGE
         INTEGER (KIND=4), INTENT(IN), OPTIONAL :: PSANumSampOutOfRange !! PGE
         INTEGER (KIND=4), PARAMETER :: npcfattr = 7
         !INTEGER (KIND=4) :: nc, numtype 
         INTEGER (KIND=C_LONG) :: nc, numtype 
         CHARACTER( LEN = 28 ) :: GranuleDAY0Z
         CHARACTER( LEN = 1  ) :: char
         REAL (KIND = 8 )      :: TAI93At0zOfGranule !!PGE
         CHARACTER( LEN = PGSd_MET_MAX_STRING_SET_L) :: ShortName, &
                                                        InputPGEVersion,&
                                                        InputVersions, &
                                                        OrbitData
         CHARACTER( LEN = PGSd_PC_VALUE_LENGTH_MAX ), DIMENSION(npcfattr) ::  &
             globalAttributeName_PCF =                                 &
             (/"PGEVERSION              ", "ProcessingCenter        ", &
               "InstrumentName          ", "ProcessingHost          ", &
               "ProcessLevel            ", "AuthorAffiliation       ", &
               "AuthorName              " /) 
         INTEGER (KIND=4), DIMENSION(npcfattr) :: lun

         CHARACTER( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: StringValue, strTemp
         INTEGER (KIND=4) :: version 
         INTEGER (KIND=4) :: status, ierr, di
         INTEGER (KIND=4) :: SW_fileid, SW_id

         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         SW_fileid = he5_swopen( L2_filename, HE5F_ACC_RDWR )
         IF( SW_fileid == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_swopen:"// TRIM(L2_filename) //&
                               " failed."
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeGlobalAttribute", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF
         
         nc = 1
         numtype = HE5T_NATIVE_INT
         ierr = he5_ehwrglatt( SW_fileid, "GranuleDay", &
                               numtype, nc, GranuleDay )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleDay" //&
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeGlobalAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         nc = 1
         numtype = HE5T_NATIVE_INT
         ierr = he5_ehwrglatt( SW_fileid, "GranuleMonth", &
                               numtype, nc, GranuleMonth )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleMonth" //&
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeGlobalAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         nc = 1
         numtype = HE5T_NATIVE_INT
         ierr = he5_ehwrglatt( SW_fileid, "GranuleYear", &
                               numtype, nc, GranuleYear )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleYear" //&
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeGlobalAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         IF( PRESENT( PSANumSampOutOfRange ) ) THEN
            nc = 1
            numtype = HE5T_NATIVE_INT
            ierr = he5_ehwrglatt( SW_fileid, "PSANumSampOutOfRange", &
                                  numtype, nc, PSANumSampOutOfRange )
            IF( ierr == -1 ) THEN
               WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "PSANumSampOutOfRange" //&
                                  " failed in file " // TRIM(L2_filename )
               ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                     "OMI_writeGlobalAttribute", zero )
               status = OZT_E_FAILURE
               RETURN 
            ENDIF
         ENDIF

         WRITE( UNIT = GranuleDAY0Z, FMT = '(I4.4,A1,I2.2,A1,I2.2,A)' ) &
             GranuleYear, '-', GranuleMonth, '-', GranuleDay, 'T00:00:00.000Z'
           
         status = PGS_TD_UTCtoTAI( GranuleDAY0Z, TAI93At0zOfGranule )
         IF( status /= PGS_S_SUCCESS ) THEN
            IF( status /= PGSTD_E_NO_LEAP_SECS ) THEN
               WRITE( msg,'(A)' ) "he5_ehwrglatt Time error:"// GranuleDAY0Z
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                     "OMI_writeGlobalAttribute", zero )
               status = OZT_E_FAILURE
               RETURN 
            ENDIF
         ENDIF
         
         nc = 1
         numtype = HE5T_NATIVE_DOUBLE
         ierr = he5_ehwrglatt( SW_fileid, "TAI93At0zOfGranule", &
                               numtype, nc, TAI93At0zOfGranule )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "TAI93At0zOfGranule" //&
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeGlobalAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF


         InputVersions = ""
         DO di = 1, SIZE( inputLUNs )
           status = L1B_getGlobalAttrname( inputLUNs(di), &
                                           inventoryMetadataName, "core" )
           version = 1
           status = PGS_MET_getPCAttr_s( inputLUNs(di), version , &
                                         TRIM(inventoryMetadataName), &
                                        "SHORTNAME", ShortName )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A,I9)' ) "get ShortName failed at LUN:", inputLUNs(di)
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "OMI_writeGlobalAttribute", zero )
               status = OZT_E_FAILURE
               return
           ENDIF
         
           InputVersions = TRIM(InputVersions) // " " // TRIM(ShortName) // ":"

           version = 1
           status = PGS_MET_getPCAttr_s( inputLUNs(di), version , &
                                         TRIM(inventoryMetadataName), &
                                        "PGEVERSION", InputPGEVersion )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A,I9)' ) "get InputPGEVersion failed at LUN:", &
                                   inputLUNs(di) 
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "OMI_writeGlobalAttribute", zero )
               status = OZT_E_FAILURE
               return
           ENDIF

           InputVersions = TRIM(InputVersions) // TRIM(InputPGEVersion)
         ENDDO

         !! Write InputVersions without the space in the beginning of string
         nc = LEN( TRIM( InputVersions ))
         IF( nc > 1 ) nc = nc - 1
         numtype = HE5T_NATIVE_CHAR
         ierr = he5_ehwrglatt( SW_fileid, "InputVersions", &
                               numtype, nc, InputVersions(2:) )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt InputVersions:"// &
                                InputVersions //&
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeGlobalAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF
      
         !! Read the global attribute 
         lun(1:npcfattr) = (/ PGEVERSION_LUN, PROCESSINGCENTER_LUN, &
                              INSTRUMENTNAME_LUN, PROCESSINGHOST_LUN, &
                              PROCESSLEVEL_LUN, AUTHORAFFILIAT_LUN, &
                              AUTHORNAME_LUN /)
         numtype = HE5T_NATIVE_CHAR
         DO di = 1, npcfattr
           status = PGS_PC_GetConfigData( lun(di), StringValue )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A,I9)' ) "get from PCF failed at LUN = ", lun(di)
              ierr = OMI_SMF_setmsg( status, msg, "OMI_writeGlobalAttribute", &
                                     zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           !! get rid of the double quote in the string values retrieved 
           !! from the PCF
           nc = deQuote( StringValue )
           IF( nc > 0 ) THEN
              ierr = he5_ehwrglatt( SW_fileid,                          &
                                    TRIM( globalAttributeName_PCF(di) ), &
                                    numtype, nc, StringValue )
              IF( ierr == -1 ) THEN
                  WRITE( msg,'(A)' ) "he5_ehwrglatt: "//  &
                         TRIM( globalAttributeName_PCF(di)) // &
                         " failed in file " // TRIM(L2_filename )
                  ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                        "OMI_writeGlobalAttribute", zero )
                  status = OZT_E_FAILURE
                  RETURN 
               ENDIF
           ENDIF
         ENDDO

         status = L1B_getGlobalAttrname( L1B_UV_FILE_LUN, &
                                         archivedMetadataName, "arch" )
         version = 1
         status = PGS_MET_getPCAttr_s( L1B_UV_FILE_LUN, version , &
                                       TRIM( archivedMetadataName ), &
                                      "ORBITDATA", OrbitData )
         IF( status /= PGS_S_SUCCESS ) THEN
            WRITE( msg,'(A,I9)' ) "get OrbitData failed at LUN:", L1B_UV_FILE_LUN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                    "OMI_writeGlobalAttribute", zero )
            OrbitData = "Unknown"
         ENDIF

         nc = LEN_TRIM(OrbitData) 
         numtype = HE5T_NATIVE_CHAR
         ierr = he5_ehwrglatt( SW_fileid, "OrbitData", numtype, nc, OrbitData )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt: OrbitData "//  &
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                   "OMI_writeGlobalAttribute", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         status = he5_swclose( SW_fileid )
         status = OZT_S_SUCCESS
         RETURN 
       END FUNCTION OMI_writeGlobalAttribute

       FUNCTION OMI_writeSwathAttribute( L2_filename, &
                                         L2_swathname, &
                                         NumTimes, &
                                         NumTimesSmallPixel, &
                                         EarthSunDistance, &
                                         VerticalCoordinate )  RESULT( status )
         USE HE5_class
         USE ISO_C_BINDING, ONLY: C_LONG
         CHARACTER( LEN = *), INTENT(IN) :: L2_filename, L2_swathname
         INTEGER (KIND=4), INTENT(IN)   :: NumTimes, NumTimesSmallPixel
         REAL (KIND=4), INTENT(IN)      :: EarthSunDistance
         CHARACTER( LEN=* ), INTENT(IN), OPTIONAL :: VerticalCoordinate
         INTEGER (KIND=4) :: ierr, status
         !INTEGER (KIND=4) :: count
         INTEGER (KIND=C_LONG) :: count
         INTEGER (KIND=4) :: SW_fileid, SW_id
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
 
         SW_fileid = he5_swopen( L2_filename, HE5F_ACC_RDWR )
         IF( SW_fileid == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_swopen:"// TRIM(L2_filename) //&
                               " failed."
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeSwathAttribute", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         !! attach the swath
         SW_id = he5_swattach( SW_fileid, L2_swathname )
         IF( SW_id == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_swattach:"// TRIM(L2_swathname) //&
                               " failed in file " // TRIM(L2_filename )
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                  "OMI_writeSwathAttribute", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF
         
         count  = 1
         status = he5_swwrattr( SW_id, "NumTimes", HE5T_NATIVE_INT, &
                                count, NumTimes )
         IF( status == -1 ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, &
                                   "Write Swath Attribute NumTimes failed.", &
                                   "OMI_writeSwathAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         count  = 1
         status = he5_swwrattr( SW_id, "NumTimesSmallPixel", HE5T_NATIVE_INT, &
                                count, NumTimesSmallPixel )
         IF( status == -1 ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, &
                          "Write Swath Attribute NumTimesSmallPixel failed.", &
                          "OMI_writeSwathAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         status = he5_swwrattr( SW_id, "EarthSunDistance", HE5T_NATIVE_FLOAT, &
                                count, EarthSunDistance )
         IF( status == -1 ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, &
                            "Write Swath Attribute EarthSunDistance failed.", &
                            "OMI_writeSwathAttribute", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         IF( PRESENT( VerticalCoordinate ) ) THEN
            count = LEN( TRIM(VerticalCoordinate) )
            status = he5_swwrattr( SW_id, "VerticalCoordinate", &
                                   HE5T_NATIVE_CHAR, &
                                   count, VerticalCoordinate )
            IF( status == -1 ) THEN
               ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, &
                           "Write Swath Attribute VerticalCoordinate failed.", &
                           "OMI_writeSwathAttribute", zero )
               status = OZT_E_FAILURE
               RETURN 
            ENDIF
         ENDIF

         status = he5_swdetach( SW_id )
         status = he5_swclose( SW_fileid )
         status = OZT_S_SUCCESS
       END FUNCTION OMI_writeSwathAttribute

!       FUNCTION PGEVersion2PhaseScience( InputPGEVersion ) RESULT( PS )
!         CHARACTER( LEN=PGSd_MET_MAX_STRING_SET_L ) :: PS
!         CHARACTER( LEN = * ), INTENT(IN) :: InputPGEVersion
!         INTEGER (KIND=4) :: ii, di, ierr 
!         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
!
!         !! search for the two '.' in the PGE version string like 11.22.33.44, 
!         !! and extract it the part before the second '.', in the above example
!         !! it would be 11.22. If there is only one '.' in the stirng, return
!         !! the whole string. If there is no '.', return error.
!         ii = INDEX( InputPGEVersion, '.' )
!         PS = InputPGEVersion( ii+1: )
!         di = INDEX( PS, '.' )
!         IF( di >= 1 ) THEN
!            ii = ii -1 + di
!         ELSE
!            ii = LEN( InputPGEVersion )
!         ENDIF
!
!         PS = InputPGEVersion(1:ii) 
!         IF( ii <= 1 ) THEN
!            WRITE( msg,'(A)' ) "InputPGEVersion:"//InputPGEVersion //&
!                               "does not have the right format xx.xx.xx"
!            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
!                                  "PGEVersion2PhaseScience", zero )
!            RETURN
!         ENDIF
!       END FUNCTION PGEVersion2PhaseScience

END MODULE L2_attr_class
