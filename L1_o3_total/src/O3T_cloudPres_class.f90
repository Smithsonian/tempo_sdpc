!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_cloudPres_class
! 
! This module contains the allocatable arrays used to store in Cloud
! inforamtion retrived from the input L2 cloud file. 
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
MODULE O3T_cloudPres_class
    USE PGS_PC_class ! define PGSd_PC_FILE_PATH_MAX, and pgs_pc functions 
                     ! include PGS_SMF.f define PGS_SMF_MAX_MSG_SIZE
    USE OMI_L2reader_class
    IMPLICIT NONE
    INTEGER(KIND = 4) :: np_cld
    INTEGER (KIND=2) :: I2foo
    REAL (KIND=4) :: R4foo
    INTEGER (KIND=1), DIMENSION(:), ALLOCATABLE :: lineMem
    INTEGER (KIND=4), DIMENSION(0:40) :: als
    INTEGER (KIND=4) :: nXtrack_cld, nTimes_cld
    INTEGER (KIND=2), DIMENSION(:), ALLOCATABLE :: ProcessingQualityFlags, &
                                                   CloudPressureI2
    REAL (KIND=4), DIMENSION(:), ALLOCATABLE :: CloudPressureR4, CloudFraction

    TYPE (L2_generic_type) :: cld_blk

    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX), DIMENSION(20) :: CLD_filename
    CHARACTER (LEN=PGS_SMF_MAX_MSG_SIZE) :: cld_swathlist, cld_parameterlist
    CHARACTER (LEN=PGS_SMF_MAX_MSG_SIZE) :: CLDshortName

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four =4, five = 5
    REAL (KIND=4), PARAMETER :: Pc_min = 250.0, Pc_max = 1013.25

    PUBLIC  :: O3T_initCLD
    PUBLIC  :: O3T_freeCLD
    PUBLIC  :: O3T_getOMICldPress

    CONTAINS
       FUNCTION O3T_initCLD( CloudPressureSource ) RESULT( status )
         USE OMI_LUN_set      ! define Logical Unit for Input and Output
         USE OMI_SMF_class    ! include PGE specific messages and OMI_SMF_setmsg
         USE PGS_MET_class
         CHARACTER( LEN = *  ), INTENT(OUT) :: CloudPressureSource

         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND=4) :: numfiles
         INTEGER (KIND=4) :: version
         INTEGER (KIND=4) :: status, ierr
 
         !! clean up everything before proceed
         CALL O3T_freeCLD

         !! test Cloud Pressure source option
         status = PGS_PC_GetConfigData( CLOUDPRESSURESOURCE_L, msg )
         IF( status /= PGS_S_SUCCESS ) THEN    ! default ="Climatology"
            CloudPressureSource = '"Climatology"'
         ELSE
            CloudPressureSource = TRIM(msg)
            IF( TRIM(msg) /= '"OMCLDRR"' .AND. &
                TRIM(msg) /= '"OMCLDO2"' .AND. &
                TRIM(msg) /= '"Climatology"' ) THEN
               CloudPressureSource = '"Climatology"'
            ENDIF
         ENDIF

         ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, "Cloud pressure source: " // &
                                TRIM(CloudPressureSource), "O3T_initCLD", zero )

         !! dealing with L2 cloud
         IF( TRIM(CloudPressureSource) /= '"Climatology"' ) THEN
            status =  L2_getFileNames( OMCLDRR_L2_LUN, numfiles, CLD_filename, &
                                       cld_swathlist )
            IF( status /= OZT_S_SUCCESS ) THEN
               WRITE( msg,'(A,I9)' ) "can't get numfiles from PCF file at LUN =", &
                                     OMCLDRR_L2_LUN
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_initCLD", zero )
               status = OZT_E_FAILURE
               RETURN
            ELSE IF( numfiles /= one ) THEN
               WRITE( msg,'(A,I9,A,I9)' ) "numfiles =", numfiles, "not equal to 1"//&
                                        " at LUN =", OMCLDRR_L2_LUN
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_initCLD", zero )
               status = OZT_E_FAILURE
               RETURN
            ENDIF
            
            version = 1
            status = PGS_MET_GetPCAttr_s( OMCLDRR_L2_LUN, version,"CoreMetadata",&
                                         "SHORTNAME", CLDshortName )
            IF( status /= PGS_S_SUCCESS ) THEN
               status = PGS_MET_GetPCAttr_s( OMCLDRR_L2_LUN, version, &
                                            "CoreMetadata.0", &
                                            "SHORTNAME", CLDshortName )
            ENDIF
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, "Cloud File Type: " // &
                                   TRIM(CLDshortName), "O3T_initCLD", zero )
            IF( INDEX( TRIM(CloudPressureSource),TRIM(CLDshortName)) == 0 ) THEN 
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, "cloud file type mismatch:"// &
                        TRIM(CLDshortName)// " and "//TRIM(CloudPressureSource), &
                        "O3T_initCLD", zero )
               status = OZT_E_FAILURE
               RETURN
            ENDIF
   
            IF( INDEX( CLDshortName, "OMCLDRR") > 0) THEN    ! OMCLDRR
                cld_parameterlist = "ProcessingQualityFlagsforO3,"//&
                                    "CloudPressureforO3,CloudFractionforO3"
            ELSE IF( INDEX( CLDshortName, "OMCLDO2") > 0) THEN    ! OMCLDO2
                cld_parameterlist = "ProcessingQualityFlags,CloudPressure,"// &
                                    "CloudFraction"
            ELSE IF( INDEX( CLDshortName, "OMSFO3") > 0) THEN !SpectralFittingResults
                cld_parameterlist = "QualityFlags,CloudPressure,"// &
                                    "CloudFraction"
            ELSE
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Unknown cloud file type:"// &
                                      TRIM(CLDshortName), "O3T_initCLD", zero )
               status = OZT_E_FAILURE
               RETURN
            ENDIF
   
            status = L2_newBlock( cld_blk, TRIM(CLD_filename(1)), &
                                  TRIM(cld_swathlist), TRIM(cld_parameterlist) )
            np_cld = cld_blk%nFields
            als(0:np_cld) = cld_blk%accuLineSize(0:np_cld)
   
            nXtrack_cld = L2_getSWdim( cld_blk, 'nXtrack' )
            nTimes_cld  = L2_getSWdim( cld_blk, 'nTimes'  )
   
            IF( INDEX( CLDshortName, "OMCLDO2") > 0) THEN    ! OMCLDO2
               status = he5_swrdattr( cld_blk%swathID, "NumTimes", nTimes_cld )
            ENDIF
   
            WRITE(msg,*) "nXtrack=", nXtrack_cld, "nTimes =", nTimes_cld
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, TRIM(CLDshortName)//TRIM(msg), &
                                   "O3T_initCLD", zero )
   
            ! allocate memory for arrays
           
            IF( INDEX( CLDshortName, "OMCLDO2") > 0) THEN    ! OMCLDO2
               ALLOCATE( lineMem( als(np_cld ) ), &
                         ProcessingQualityFlags(nXtrack_cld),&
                         CloudFraction(nXtrack_cld), CloudPressureI2(nXtrack_cld),&
                         STAT=ierr )
            ELSE
               ALLOCATE( lineMem( als(np_cld ) ), &
                         ProcessingQualityFlags(nXtrack_cld),&
                         CloudFraction(nXtrack_cld), CloudPressureR4(nXtrack_cld),&
                         STAT=ierr )
            ENDIF
            IF( ierr /= zero ) THEN
               ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                     "radiance and geo allocation failure", &
                                     "O3T_initCLD", zero )
               status = OZT_E_FAILURE
               RETURN
            ENDIF
         ENDIF
       END FUNCTION O3T_initCLD

       SUBROUTINE O3T_getOMICldPress( iLine, pcArray, nXtrack )
         INTEGER(KIND=4), INTENT(IN) :: iLine, nXtrack
         REAL (KIND = 4), DIMENSION(:), INTENT(OUT) :: pcArray
         INTEGER(KIND=4) :: ierr, status
         
         IF( nXtrack /= nXtrack_cld ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                   "nXtrack not equal to nXtrack_cld", &
                                   "O3T_getOMICldPress", zero )
            CALL EXIT(1) 
         ENDIF

         status = L2_getLine( cld_blk, iLine, lineMem )
         ProcessingQualityFlags(:) = TRANSFER( lineMem(als(0)+1:als(1)), &
                                               I2foo, nXtrack_cld )
         IF( INDEX( CLDshortName, "OMCLDRR" ) > 0 ) THEN
            CloudPressureR4(:) = TRANSFER( lineMem(als(1)+1:als(2)), &
                                           R4foo, nXtrack_cld )
            pcArray(1:nXtrack) = CloudPressureR4(1:nXtrack_cld)
         ELSE IF( INDEX( CLDshortName, "OMCLDO2" ) > 0 ) THEN
            CloudPressureI2(:) = TRANSFER( lineMem(als(1)+1:als(2)), &
                                           I2foo, nXtrack_cld )
            pcArray(1:nXtrack) = CloudPressureI2(1:nXtrack_cld)
         ENDIF
         CloudFraction(:) = TRANSFER( lineMem(als(2)+1:als(3)), &
                                      R4foo, nXtrack_cld )
       END SUBROUTINE O3T_getOMICldPress

       SUBROUTINE O3T_freeCLD
         INTEGER (KIND=4) :: status
         IF( ALLOCATED( lineMem       ) ) DEALLOCATE( lineMem       )
         IF( ALLOCATED( CloudFraction ) ) DEALLOCATE( CloudFraction )
         IF( ALLOCATED( CloudPressureI2 ) ) DEALLOCATE( CloudPressureI2 )
         IF( ALLOCATED( CloudPressureR4 ) ) DEALLOCATE( CloudPressureR4 )
         IF( ALLOCATED( ProcessingQualityFlags ) ) &
                                          DEALLOCATE( ProcessingQualityFlags )

         ! close data block structure
         CALL L2_disposeBlock( cld_blk )
       END SUBROUTINE O3T_freeCLD

END MODULE O3T_cloudPres_class
