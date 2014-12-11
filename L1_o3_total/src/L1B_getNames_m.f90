MODULE L1B_getNames_m
    USE HE4_class
    USE PGS_PC_class  ! define PGSd_PC_FILE_PATH_MAX, and PGS_PC functions 
                      ! include PGS_SMF.f define PGS_SMF_MAX_MSG_SIZE
    USE OMI_SMF_class ! include PGE specific messages and OMI_SMF_setmsg

    IMPLICIT NONE

    INTEGER(KIND=4), PARAMETER, PRIVATE :: zero = 0
    PUBLIC  :: L1B_getNames

    CONTAINS
       FUNCTION L1B_getNames( L1B_file_LUN, numfiles, L1B_filenames, &
                              L1B_swathlist ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN) :: L1B_file_LUN
         CHARACTER( LEN = * ), DIMENSION(:), INTENT(OUT) :: L1B_filenames
         CHARACTER( LEN = * ),               INTENT(OUT) :: L1B_swathlist
         INTEGER (KIND=4), INTENT(OUT) :: numfiles
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         CHARACTER( LEN = LEN(L1B_swathlist) ) :: swathlist
         INTEGER (KIND=4) :: di, version, nswath, strbufsize
         INTEGER (KIND=4) :: status, ierr
         INTEGER (KIND=4) :: SW_fileid 

         !! Get number of L1B output files from PCF
         !! and make sure there is only one 

         status = PGS_PC_getnumberoffiles( L1B_file_LUN,  numfiles )
         IF( status /= PGS_S_SUCCESS ) THEN
            WRITE( msg,'(A,I9)' ) "can't get numfiles from PCF file at LUN =", &
                                  L1B_file_LUN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT,msg,"L1B_getNames",zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF
         !! Get the L1B file name from the PCF.
         DO di = 1, numfiles
           version = di
           status = PGS_PC_getreference( L1B_file_LUN, version, &
                                         L1B_filenames(di) )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(2(A,I9))' ) "get filename from PCF file at LUN =", &
                                      L1B_file_LUN, " Version =", version
              ierr = OMI_SMF_setmsg( OZT_E_INPUT,msg,"L1B_getNames",zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           SW_fileid = swopen( TRIM(L1B_filenames(di)), DFACC_READ )
           IF( SW_fileid == -1 ) THEN
              WRITE( msg,'(A)' ) "swopen:"// TRIM(L1B_filenames(di)) // &
                                 " failed."
              ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L1B_getNames",zero)
              status = OZT_E_FAILURE
              RETURN 
           ENDIF

           nswath = swinqswath( TRIM(L1B_filenames(di)), swathlist, &
                                strbufsize )
           IF( di == 1 ) THEN
              L1B_swathlist = TRIM(swathlist)
           ELSE IF( TRIM(L1B_swathlist) /= TRIM(swathlist) ) THEN
              WRITE( msg,'(A)' ) "file contains different swathlist:"// &
                                 TRIM( L1B_swathlist ) // "," // &
                                 TRIM(swathlist) 
              ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L1B_getNames",zero)
              status = OZT_E_FAILURE
              RETURN 
           ENDIF 

           status = swclose( SW_fileid )
         ENDDO

         status = OZT_S_SUCCESS
         RETURN 
       END FUNCTION L1B_getNames

END MODULE L1B_getNames_m
