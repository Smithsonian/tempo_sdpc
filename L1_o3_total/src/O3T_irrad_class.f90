!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE O3T_irrad_class
! 
! read in solar irradiance from either normal or backup solar irradiance files.
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
MODULE O3T_irrad_class
    USE L1B_class    ! L1B_class contains module to read OMI L1B 
                     ! geolocation, radiance, and irradiance parameters
    USE HE4_class
    USE PGS_PC_class ! define PGSd_PC_FILE_PATH_MAX, and pgs_pc functions 
                     ! include PGS_SMF.f define PGS_SMF_MAX_MSG_SIZE
    USE OMI_SMF_class    ! include PGE specific messages and OMI_SMF_setmsg
    IMPLICIT NONE
    REAL (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: irradiance, &
                                                    irrPrecision, &
                                                    irrWavelength
    INTEGER (KIND = 2), DIMENSION(:,:), ALLOCATABLE :: irrQAflags
    INTEGER (KIND = 2) :: irrMeasurementFlags
    INTEGER (KIND = 1) :: instID_irr
    INTEGER (KIND = 1) :: irr_error = 1
    INTEGER (KIND = 1) :: irr_warning = 1, irr_any = 1
    INTEGER (KIND = 4) :: nTimes_irr, nXtrack_irr, nWavel_irr
    REAL (KIND = 4) :: EarthSunDistanceIRR
    CHARACTER( LEN = 6 ) :: IRR_FILE_TYPE  !! Normal or Backup
    INTEGER( KIND = 4 ) :: USED_L1BIRR_LUN
    LOGICAL :: NORMAL_L1BIRR_MISSING = .FALSE.
    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four = 4, five = 5
    PUBLIC  :: O3T_getIRR
    PUBLIC  :: O3T_freeIRR
    PUBLIC  :: O3T_AdjustIRREarthSun
!    PUBLIC  :: O3T_irrRepair

    CONTAINS

!!
!  In general, two solar files is stage, one is the daily solar meassurement,
!  and the other is the back-up solar irradiance file. This function first 
!  decide which of these two IRR files will be used in the processing
!  based on the quality flags in the file. Then if the optional input wl_com
!  is presnet, the solar irradiance values is interpolated to the wl_com, 
!  otherwise it is directly read from the IRR file. 
!!
       FUNCTION O3T_getIRR( IRR_filename, IRR_swathname, wl_com ) &
                            RESULT( status )
         USE OMI_LUN_set      ! define Logical Unit for Input and Output
         CHARACTER( LEN = * ), INTENT(IN) :: IRR_filename, IRR_swathname
                              !! irradiance file and swath names
         REAL (KIND = 4), DIMENSION(:), INTENT(IN), OPTIONAL :: wl_com
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND=4) :: iLine
         INTEGER (KIND=4) :: status, ierr ! , status1
         INTEGER (KIND=4) :: nWavelCoef_irr
         TYPE (L1B_radirr_type) :: irr_blk

         EarthSunDistanceIRR = L1Bga_EarthSunDistance( IRR_filename, &
                                                       IRR_swathname )
         IF( EarthSunDistanceIRR <= 0.0 ) THEN
            WRITE( msg,'(A,g16.7,A)' ) "IRR EarthSunDistance = ", &
                     EarthSunDistanceIRR, "Reset it to 1"
            EarthSunDistanceIRR = 1.0
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_getIRR", zero )
         ENDIF

         status = L1Bri_open( irr_blk, IRR_filename, IRR_swathname )
         IF( status /= OZT_S_SUCCESS ) THEN
            WRITE( msg,'(A)' ) "L1Bri_open "// TRIM(IRR_swathname) //&
                               " in file " // TRIM(IRR_filename) // " failed."
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_getIRR", zero )
            status = OZT_E_FAILURE
            RETURN
         ELSE
            WRITE( msg,'(A)' ) "L1Bri_open "// TRIM(IRR_swathname) //&
                              " in file "//TRIM(IRR_filename)//" successfully."
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, "O3T_getIRR", four )
         ENDIF

         !! obtain sizes of dimensions defined in swath
         status = L1Bri_getSWdims( irr_blk, nTimes_irr, &
                                   nXtrack_irr, nWavel_irr, &
                                   nWavelCoef_irr )
         IF( status /= OZT_S_SUCCESS ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
                                 "L1Bri_getSWdims failed.", "O3T_getIRR", zero )
            status = OZT_E_FAILURE
            RETURN
         ELSE
            WRITE( msg,'(A,4I4)' ) "(nTime, nXtrack, nWavel, nWavelCoef)=", &
                       nTimes_irr, nXtrack_irr, nWavel_irr, nWavelCoef_irr
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, "O3T_getIRR", four )
         ENDIF

         !! set nWavel_irr to the input Wl array size, 
         !! if input wavelength is used.
         IF( PRESENT( wl_com ) ) THEN
            nWavel_irr = SIZE( wl_com )
            WRITE( msg,'(A,I3)' ) "Common array of fixed wavelength used." // &
                                 "nWave_irr =", nWavel_irr
            ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, "O3T_initRAD", four )
         ENDIF

         ! allocate memory for arrays
         CALL O3T_freeIRR 
         ALLOCATE( irradiance(nWavel_irr,nXtrack_irr), &
                   irrPrecision(nWavel_irr,nXtrack_irr), &
                   irrQAflags(nWavel_irr,nXtrack_irr), &
                   irrWavelength(nWavel_irr,nXtrack_irr), STAT=ierr )
         IF( ierr /= zero ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                  "irradiance allocation failure", &
                                  "O3T_getIRR", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         iLine = 0
         DO iLine = 0, nTimes_irr-1
           IF( PRESENT( wl_com ) ) THEN
              status = L1Bri_getLineWL( irr_blk, iLine, &
                                    wl_com, nWavel_irr, &
                                    RadIrr_k = irradiance, &
                                    RadIrrPrecision_k =irrPrecision, &
                                    PixelQualityFlags_k = irrQAflags )
              IF( status /= OZT_S_SUCCESS ) THEN
                 ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
                            "L1Bri_getLineWL failed.", "O3T_getIRR", zero )
                 status = OZT_E_FAILURE
                 RETURN
              ENDIF
           ELSE
              status = L1Bri_getLine( irr_blk, iLine, &
                                    RadIrr_k = irradiance, &
                                    RadIrrPrecision_k =irrPrecision, &
                                    PixelQualityFlags_k = irrQAflags, &
                                    Wavelength_k = irrWavelength ) 
              IF( status /= OZT_S_SUCCESS ) THEN
                 ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
                               "L1Bri_getLine failed.", "O3T_getIRR", zero )
                 status = OZT_E_FAILURE
                 RETURN
              ENDIF
           ENDIF

           status = L1Bri_getMeasurementQF( irr_blk, iLine, &
                                            irrMeasurementFlags )
           status = L1Bri_getInstConfigId( irr_blk, iLine, instID_irr )
           !! No error, so no need for read another scan
           IF( irrMeasurementFlags == 0 ) EXIT
         ENDDO

         !! close data irradiance block structure
         status = L1Bri_close( irr_blk )
         IF( status /= OZT_S_SUCCESS ) THEN
            WRITE( msg,'(A)' ) "L1Bri_close failed." 
            ierr = OMI_SMF_setmsg( status, msg, "O3T_getIRR", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF
         
         !! Replace solar data with unusable QA flags with those that
         !! have usable QA flags in the adjacent scan position.
    !!!!!CALL O3T_irrRepair

       END FUNCTION O3T_getIRR

       SUBROUTINE O3T_freeIRR
         IF( ALLOCATED( irradiance    ) ) DEALLOCATE( irradiance    )
         IF( ALLOCATED( irrPrecision  ) ) DEALLOCATE( irrPrecision  )
         IF( ALLOCATED( irrQAflags    ) ) DEALLOCATE( irrQAflags    )
         IF( ALLOCATED( irrWavelength ) ) DEALLOCATE( irrWavelength )
       END SUBROUTINE O3T_freeIRR

       SUBROUTINE O3T_AdjustIRREarthSun( EarthSundistanceRAD )
         REAL (KIND = 4), INTENT(IN) :: EarthSunDistanceRAD
         REAL (KIND = 4) :: Ratio
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         IF( EarthSunDistanceRAD <= 0.0 ) THEN
            WRITE( msg,'(A,g16.7,A)' ) "Radiance EarthSunDistance = ", &
                     EarthSunDistanceRAD, ", No adjustment are made to IRR"
            RETURN
         ELSE
            Ratio = EarthSunDistanceIRR/EarthSunDistanceRAD
            IF( 0.5 < Ratio .AND. Ratio < 1.5 ) THEN
               irradiance = irradiance*Ratio**2 
            ELSE
               WRITE( msg,'(A,F10.4,A)' ) "EarthSunDistance Ratio= ", &
                     Ratio, ", out of range, no adjustment are made to IRR"
               RETURN
            ENDIF
         ENDIF
       END SUBROUTINE O3T_AdjustIRREarthSun

!       SUBROUTINE O3T_irrRepair
!         INTEGER (KIND = 4 ) :: iX, iwl
!
!         DO iX = 1, SIZE( irrQAflags, 2 )
!           DO iwl = 1, SIZE( irrQAflags, 1 )
!             IF( IBITS( irrQAflags(iwl,iX), 0, 3 ) > 0 ) THEN
!                IF( iX == 1 ) THEN
!                   irradiance(iwl,iX) = irradiance(iwl,iX+1)
!                   irrQAflags(iwl,iX) = irrQAflags(iwl,iX+1)
!                ELSE IF( iX == SIZE( irrQAflags, 2 ) ) THEN
!                   irradiance(iwl,iX) = irradiance(iwl,iX-1)
!                   irrQAflags(iwl,iX) = irrQAflags(iwl,iX-1)
!                ELSE
!                   IF( IBITS( irrQAflags(iwl,iX-1), 0, 3 ) == 0 .AND. &
!                       IBITS( irrQAflags(iwl,iX+1), 0, 3 ) == 0 ) THEN
!                      irrQAflags(iwl,iX) = IOR( irrQAflags(iwl,iX-1), &
!                                                irrQAflags(iwl,iX+1) )
!
!                      irradiance(iwl,iX) = 0.5*( irradiance(iwl,iX-1) + &
!                                                 irradiance(iwl,iX+1) )
!                   ELSE IF( IBITS( irrQAflags(iwl,iX-1), 0, 3 ) == 0 ) THEN
!                      irrQAflags(iwl,iX) = irrQAflags(iwl,iX-1)
!                      irradiance(iwl,iX) = irradiance(iwl,iX-1)
!                   ELSE IF( IBITS( irrQAflags(iwl,iX+1), 0, 3 ) == 0 ) THEN
!                      irrQAflags(iwl,iX) = irrQAflags(iwl,iX+1)
!                      irradiance(iwl,iX) = irradiance(iwl,iX+1)
!                   ENDIF
!                ENDIF
!             ENDIF
!           ENDDO
!         ENDDO
!       END SUBROUTINE O3T_irrRepair
       
END MODULE O3T_irrad_class
