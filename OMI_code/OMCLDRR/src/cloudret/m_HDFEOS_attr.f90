MODULE m_HDFEOS_attr

    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four =4, five = 5
    INTEGER (KIND=4), PARAMETER:: HE5T_NATIVE_INT = 0, HE5T_NATIVE_FLOAT = 10, &
                                  HE5T_NATIVE_DOUBLE = 11, HE5T_NATIVE_CHAR = 56

    PUBLIC :: CLDRR_writeGlobalAttr
    PUBLIC :: CLDRR_writeSwathAttr
    PRIVATE :: PGEVersion2PhaseScience
    PUBLIC :: deQuote

    CONTAINS

       FUNCTION CLDRR_WriteGlobalAttr( SW_fileid, GranuleYear, &
                                          GranuleMonth, &
                                          GranuleDay ) &
                            RESULT( status )
         USE m_LUN_set
!        INCLUDE 'PGS_TD.f'
         INCLUDE 'PGS_TD_3.f'
         INCLUDE 'PGS_MET.f'
 INCLUDE 'PGS_PC.f'
 INCLUDE 'PGS_PC_9.f'
 INCLUDE 'PGS_SMF.f'
 INCLUDE 'PGS_MET_13.f'
 INCLUDE 'PGS_OMI_1900.f'
 INCLUDE 'PGS_OMCLDRR_52251.f'

    INTEGER(KIND=4), EXTERNAL :: OMI_SMF_setmsg, he5_ehwrglatt, PGS_MET_getPCAttr_s
    INTEGER(KIND=4), EXTERNAL :: PGS_PC_GetConfigData   

         INTEGER (KIND=4), INTENT(IN) :: GranuleDay, GranuleMonth,  &
                                         GranuleYear  !!PGE
! character (len=*), intent(in) :: outfile
  integer (kind=4), intent(in) :: SW_fileid

         INTEGER (KIND=4), PARAMETER :: npcfattr = 7
         INTEGER (KIND=4) :: nc, numtype 
         CHARACTER( LEN = 28 ) :: GranuleDAY0Z
         CHARACTER( LEN = 1  ) :: char
         REAL (KIND = 8 )      :: TAI93At0zOfGranule !!PGE
         CHARACTER( LEN = 200) :: ShortName, &
                                                        InputPGEVersion,&
                                                        InputVersions 
         CHARACTER( LEN = 200 ), DIMENSION(npcfattr) ::  &
             globalAttributeName_PCF =                                 &
             (/"PGEVERSION              ", "ProcessingCenter        ", &
               "InstrumentName          ", "ProcessingHost          ", &
               "ProcessLevel            ", "AuthorAffiliation       ", &
               "AuthorName              " /) 
         INTEGER (KIND=4), DIMENSION(npcfattr) :: lun

         CHARACTER( LEN = 200 ) :: StringValue, strTemp, OrbitData
         INTEGER (KIND=4) :: version 
         INTEGER (KIND=4) :: status, ierr, di
         CHARACTER( LEN = 200 ) :: msg
         INTEGER (KIND=4), EXTERNAL :: PGS_TD_UTCtoTAI

         nc = 1
         numtype = HE5T_NATIVE_INT
         ierr = he5_ehwrglatt( SW_fileid, "GranuleDay", &
                               numtype, nc, GranuleDay )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleDay" // " failed "
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                  "CLDRR_writeGlobalAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF

         nc = 1
         numtype = HE5T_NATIVE_INT
         ierr = he5_ehwrglatt( SW_fileid, "GranuleMonth", &
                               numtype, nc, GranuleMonth )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleMonth" // " failed "
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                  "CLDRR_writeGlobalAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF

         nc = 1
         numtype = HE5T_NATIVE_INT
         ierr = he5_ehwrglatt( SW_fileid, "GranuleYear", &
                               numtype, nc, GranuleYear )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleYear" // " failed "
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                  "CLDRR_writeGlobalAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF

         WRITE( UNIT = GranuleDAY0Z, FMT = '(I4.4,A1,I2.2,A1,I2.2,A)' ) &
             GranuleYear, '-', GranuleMonth, '-', GranuleDay, 'T00:00:00.000Z'
         status = PGS_TD_UTCtoTAI( GranuleDAY0Z, TAI93At0zOfGranule )
         IF( status /= PGS_S_SUCCESS ) THEN
            IF( status /= PGSTD_E_NO_LEAP_SECS ) THEN
               WRITE( msg,'(A)' ) "he5_ehwrglatt Time error:"// GranuleDAY0Z
               ierr = OMI_SMF_setmsg(omcldrr_f_hdfeos, msg, &
                                     "CLDRR_writeGlobalAttribute", zero )
               call exit(1)
               RETURN 
            ENDIF
         ENDIF
         
         nc = 1
         numtype = HE5T_NATIVE_DOUBLE
         ierr = he5_ehwrglatt( SW_fileid, "TAI93At0zOfGranule", &
                               numtype, nc, TAI93At0zOfGranule )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "TAI93At0zOfGranule" // " failed "
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                  "CLDRR_writeGlobalAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF

         InputVersions = ""
         lun(1:2) = (/ l1b_lun, irr1b_file /)
         DO di = 1, 2
           version = 1
           status = PGS_MET_getPCAttr_s( lun(di), version , "CoreMetadata.0", &
                                        "SHORTNAME", ShortName )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A,I0)' ) "get ShortName failed at LUN:", lun(di) 
              ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                    "CLDRR_writeGlobalAttribute", zero )
               call exit(1)
               return
           ENDIF
           InputVersions = TRIM(InputVersions) // " " // TRIM(ShortName) // ":"

           version = 1
           status = PGS_MET_getPCAttr_s( lun(di), version , "CoreMetadata.0", &
                                         "PGEVERSION", InputPGEVersion )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A,I0)' ) "get InputPGEVersion failed at LUN:", lun(di) 
              ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                    "CLDRR_writeGlobalAttribute", zero )
               call exit(1)
               return
           ENDIF

           InputVersions = TRIM(InputVersions) //  &
                           TRIM(PGEVersion2PhaseScience(InputPGEVersion)) 
         ENDDO

         !! Write InputVersions without the space in the beginning of string
         nc = LEN( TRIM( InputVersions ))
         IF( nc > 1 ) nc = nc - 1
         numtype = HE5T_NATIVE_CHAR
         ierr = he5_ehwrglatt( SW_fileid, "InputVersions", &
                               numtype, nc, InputVersions(2:) )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt InputVersions:"// &
                                InputVersions //" failed "
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                  "CLDRR_writeGlobalAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF
      
         !! Read the global attribute 
         lun(1:npcfattr) = (/ PGEVERSION_LUN, PROCESSINGCENTER_LUN, &
                              INSTRUMENTNAME_LUN, PROCESSINGHOST_LUN, &
                              PROCESSLEVEL_LUN, AUTHORAFFILIATION_LUN, &
                              AUTHORNAME_LUN /)
         numtype = HE5T_NATIVE_CHAR
         DO di = 1, npcfattr
           status = PGS_PC_GetConfigData( lun(di), StringValue )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(A,I0)' ) "get from PCF failed at LUN = ", lun(di)
              ierr = OMI_SMF_setmsg( status, msg, "CLDRR_writeGlobalAttribute", &
                                     zero )
              call exit(1)
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
                         TRIM( globalAttributeName_PCF(di)) // " failed "
                  ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                        "CLDRR_writeGlobalAttribute", zero )
                  call exit(1)
                  RETURN 
               ENDIF
           ENDIF
         ENDDO

  version=1
  status = pgs_met_getPCAttr_s(L1B_LUN, version , "ArchiveMetadata.0", &
                               "ORBITDATA",OrbitData)
  IF(status /= PGS_S_SUCCESS ) THEN
    ierr = OMI_SMF_setmsg( OMCLDRR_W_MET, "get OrbitData from L1B failed", &
                              "MetadataModule", 0 )
  ENDIF
         numtype = HE5T_NATIVE_CHAR
         nc = LEN(TRIM(OrbitData)) 
         ierr = he5_ehwrglatt( SW_fileid, "OrbitData", &
                               numtype, nc, TRIM(OrbitData) )
         IF( ierr == -1 ) THEN
            WRITE( msg,'(A)' ) "he5_ehwrglatt OrbitData failed "
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                  "CLDRR_writeGlobalAttribute", zero )
            call exit(1)
            RETURN
         ENDIF

         RETURN 
       END FUNCTION CLDRR_WriteGlobalAttr

       FUNCTION CLDRR_WriteSwathAttr( SW_id, fn, swn)  RESULT( status )
 INCLUDE 'PGS_SMF.f'
 INCLUDE 'PGS_MET_13.f'
 INCLUDE 'PGS_OMI_1900.f'
 INCLUDE 'PGS_OMCLDRR_52251.f'

    INTEGER, EXTERNAL :: OMI_SMF_setmsg, he5_swwrattr, he5_ehwrglatt, PGS_MET_getPCAttr_s    
    INTEGER, EXTERNAL :: swopen, swattach, swrdattr, swdetach, swclose 
  integer, parameter :: DFACC_READ = 1   

         CHARACTER( LEN=* ), INTENT(IN) :: fn, swn 
  integer (kind=4), intent(in) :: SW_id 
         INTEGER (KIND=4) :: ierr, status, swfid, swid
         INTEGER (KIND=4) :: count, NumTimes, NumTimesSmallPix
  REAL (KIND=4) :: EarthSundistance
  character(len=200) :: VerticalCoordinate = 'Total Column'

!Get NumTimes & NumTimesSmallPixel from L1B

        !! open the  swath file
        swfid = swopen( fn, DFACC_READ )
        IF( swfid < zero ) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg( OMI_E_FILE_OPEN, fn, "WriteSwathAttr", zero )
           RETURN
        ENDIF

        !! attach to the swath
        swid = swattach( swfid, swn )
        IF( swid < zero ) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg( OMI_E_SWATH_ATTACH, swn, "WriteSwathAttr", zero )
           RETURN
        ENDIF

        !! retrieve the dimension info from the swath file.
        status = swrdattr( swid, "NumTimes", NumTimes )
        If( status < zero ) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, "get NumTimes size failed", &
                                 "WriteSwathAttr", zero )
           RETURN
        ENDIF
         
        status = swrdattr( swid, "NumTimesSmallPixel", NumTimesSmallPix )
        If( status < zero ) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
                                 "get NumTimesSmallPixel size failed", &
                                 "WriteSwathAttr", zero )
           RETURN
        ENDIF

        status = swrdattr( swid, "EarthSunDistance", EarthSundistance )
        If( status < zero ) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg( OMI_E_HDFEOS, &
                                 "get EarthSundistance failed", &
                                 "WriteSwathAttr", zero )
           RETURN
        ENDIF

        !! detach and close the swath file.  No need for 
        !! error checking, for an error is unlikely to occur here.
        ierr = swdetach( swid )
        ierr = swclose( swfid )

         count  = 1
         status = he5_swwrattr( SW_id, "NumTimes", HE5T_NATIVE_INT, &
                                count, NumTimes )
         IF( status == -1 ) THEN
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, &
                                   "Write Swath Attribute NumTimes failed.", &
                                   "CLDRR_writeSwathAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF

         count  = 1
         status = he5_swwrattr( SW_id, "NumTimesSmallPixel", HE5T_NATIVE_INT, &
                                count, NumTimesSmallPix )
         IF( status == -1 ) THEN
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, &
                                   "Write Swath Attribute NumTimesSmallPix failed.", &
                                   "CLDRR_writeSwathAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF

         status = he5_swwrattr( SW_id, "EarthSunDistance", HE5T_NATIVE_FLOAT, &
                                count, EarthSundistance )
         IF( status == -1 ) THEN
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, &
                            "Write Swath Attribute EarthSundistance failed.", &
                            "CLDRR_writeSwathAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF

         count = LEN( TRIM(VerticalCoordinate) )
         status = he5_swwrattr( SW_id, "VerticalCoordinate", HE5T_NATIVE_CHAR, &
                                count, VerticalCoordinate )
         IF( status == -1 ) THEN
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, &
                           "Write Swath Attribute VerticalCoordinate failed.", &
                           "CLDRR_writeSwathAttribute", zero )
            call exit(1)
            RETURN 
         ENDIF
       END FUNCTION CLDRR_WriteSwathAttr

   
       FUNCTION PGEVersion2PhaseScience( InputPGEVersion ) RESULT( PS )
         INCLUDE 'PGS_MET.f'
 INCLUDE 'PGS_OMCLDRR_52251.f'
         CHARACTER( LEN=PGSd_MET_MAX_STRING_SET_L ) :: PS
         CHARACTER( LEN = * ), INTENT(IN) :: InputPGEVersion
         INTEGER (KIND=4) :: ii, di, ierr, OMI_SMF_setmsg 
         CHARACTER( LEN = 200 ) :: msg

         !! search for the two '.' in the PGE version string like 11.22.33.44, 
         !! and extract it the part before the second '.', in the above example
         !! it would be 11.22. If there is only one '.' in the stirng, return
         !! the whole string. If there is no '.', return error.
         ii = INDEX( InputPGEVersion, '.' )
         PS = InputPGEVersion( ii+1: )
         di = INDEX( PS, '.' )
         IF( di >= 1 ) THEN
            ii = ii -1 + di
         ELSE
            ii = LEN( InputPGEVersion )
         ENDIF

         PS = InputPGEVersion(1:ii) 
         IF( ii <= 1 ) THEN
            WRITE( msg,'(A)' ) "InputPGEVersion:"//InputPGEVersion //&
                               "does not have the right format xx.xx.xx"
            ierr = OMI_SMF_setmsg( omcldrr_f_hdfeos, msg, &
                                  "PGEVersion2PhaseScience", zero )
            RETURN
         ENDIF
       END FUNCTION PGEVersion2PhaseScience

       FUNCTION deQuote( StringValue ) RESULT( nc )
         CHARACTER(LEN=*) :: StringValue
         CHARACTER(LEN=LEN(StringValue)) :: strTemp
         INTEGER :: nc, di

         di = INDEX( StringValue, '"' )
         IF( di > 0 ) THEN
            strTemp = StringValue(di+1:)
            di = INDEX( strTemp, '"', BACK = .TRUE. )
            IF( di > 0 ) THEN
               StringValue = strTemp(1:di-1)
            ELSE
               StringValue = strTemp(1:)
            ENDIF
         ENDIF
         nc = LEN( TRIM(StringValue) )
         RETURN
       END FUNCTION deQuote

END MODULE m_HDFEOS_attr
