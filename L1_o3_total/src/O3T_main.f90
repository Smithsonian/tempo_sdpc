!!****************************************************************************
!!F90
!
!!Description:
!
!  PROGRAM O3T_mainNVAdj
!
!  This is the main program that implements the three-step process of
!  successive improvement of the estimation of total column amount ozone.
!  Step 1: The alogrithm uses an pair of wavelengths to derive an estimation
!  of the column amount ozone. This is done by iterative deriving reflectivity
!  from 331 nm and ozone from the single B wavelength 318nm, until the N-values
!  at this pair are the same for the measurements and for the calculation using
!  the reflectivity and ozone values.
!  Step 2: Ozone and temperature climatologies are applied at all levels to
!  account for seasonal ad latitudinal variations in profile shape.
!  Step 3: The step 2 ozone estimation is modified to correct for wavelength
!  dependent effects (torpospheric aerosol and sun glint) and local upper level
!  profile shape effects.
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
PROGRAM O3T_mainNVAdj
    USE O3T_irrad_class  ! contains irradiance parameters
    USE O3T_radgeo_class ! contains radiance and geolocation parameters
    USE O3T_nval_class, ONLY: wl_com, nwl_com, nvRRS, &  ! , nvELS
                              iwl_ozone, iwl_refl_l, iwl_refl_h, iwl_mix, &
                              O3T_nval_setup, O3T_nval_dispose, &
                              O3T_getLambdaSet
    USE O3T_dndx_class, ONLY: O3T_dndx_setup, O3T_dndx_dispose ! nwl_sub,
    USE O3T_cloudPres_class
    USE O3T_apriori_class
    USE O3T_lpolycoef_class
    USE O3T_lpolyinterp_class
    USE O3T_pixel_class
    USE O3T_class
    USE O3T_QA_class
    USE O3T_const
    USE O3T_so2_class
    USE OMI_SMF_class    ! include PGE specific messages and OMI_SMF_setmsg
    USE L1B_getNames_m
    USE OMI_LUN_set
    USE O3T_omto3_fs
    USE OMI_copyHE4toHE5_class
    USE O3T_L2output_class
    USE OMI_L2writer_class
    USE L1B_metaData_class
    USE L2_metaData_class
    USE L2_attr_class
    USE m_nvalc
    USE m_anomflg

    IMPLICIT NONE
    INTEGER (KIND=4), PARAMETER :: zero = 0, one = 1!, two = 2, four = 4
    INTEGER (KIND=4), PARAMETER :: nLinesPerWrite = 100
    INTEGER (KIND=4), PARAMETER :: nGeoF = 15, nDatF = 25, nwlA = 4
    INTEGER (KIND=4), DIMENSION(5) :: dims
    CHARACTER (LEN=256) :: dimList
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: OMTO3_fn
    CHARACTER (LEN=256) :: OMTO3_swathname
    INTEGER (KIND=4) :: OMTO3_fileid, OMTO3_swid
    CHARACTER (LEN=2) :: satname
    TYPE (DFHE5_T), DIMENSION(nGeoF) :: gf_omto3
    TYPE (DFHE5_T), DIMENSION(nDatF) :: df_omto3
    TYPE (DFHE5_T), DIMENSION( 1   ) :: df_omto3wl, df_caladj
    TYPE (DFHE5_T), DIMENSION( 3   ) :: df_copy
    TYPE (L2_generic_type) :: wlblk, geoblk, datablk, calblk
    REAL (KIND=4) :: wl_cutoff = 360.0  !dndx cutoff wl in nm
    REAL (KIND=4) :: sza_p, vza_p, phi, pt, pc, lat, lon, hgt
    TYPE (O3T_lpoly_cden_type) :: LUT_cden_blk
    TYPE (O3T_lpoly_coef_type) :: coefs
    TYPE (O3T_pixgeo_type)     :: pixGEO
    TYPE (O3T_pixcover_type)   :: pixSURF
    TYPE (L1BECSMETA_T)        :: L1BUV2coreArch
    TYPE (L1BECSMETA_T)        :: L1BIRRcoreArch
    TYPE (L2PARAM_T), DIMENSION(1) :: L2_parameters
    INTEGER (KIND=4), DIMENSION(11):: LUNinputPointer
    INTEGER (KIND=4) :: mcfLUN
    REAL (KIND=4) :: doz_limit = 5.0, guesoz, stp1oz, stp2oz, stp3oz, dr
    REAL (KIND=4), DIMENSION(NLYR) :: stp2prf, eff, aprfoz
    REAL (KIND=4)  :: aerind, so2ind, soilim, pathl, oz_cld !, cloudcov = 0.0
    INTEGER (KIND=4), DIMENSION(4) :: iso2w
    REAL (KIND=4), DIMENSION(5) :: o3abs, so2abs

    REAL (KIND=4)  :: latPre, latNow
    LOGICAL :: absrfl, skipit, maxitr, descendQ, nXevenQ
    LOGICAL :: stp1oz_valid, stp2oz_valid
    INTEGER (KIND=4) :: iwl_oz, iwl_refl, iplow
    INTEGER (KIND=4), DIMENSION(nwlA) :: iwlArray
    INTEGER (KIND=1) :: algflg = 0, mqaL2 = 0, cld_errflg
    INTEGER (KIND=2) :: QAflags = 0, radBadPixflgs = 0, errflgs = 0
    INTEGER (KIND=4) :: year, month, day, jday
    INTEGER (KIND=2) :: nise_flag
    LOGICAL :: radLMissing = .FALSE., instIDmismatch = .FALSE., bit7Q = .FALSE.
    INTEGER (KIND=4) :: iwl, iX, iLine, iLine_s, iLine_b, iLine_e, iT, ii, nLw, iLat
    INTEGER (KIND=4) :: pixID
    INTEGER (KIND=4) :: status, ierr
    !INTEGER (KIND=4) :: version
    !INTEGER (KIND=4) :: orbitNumber = 999

    INTEGER (KIND=4) :: numfiles
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX), DIMENSION(100) :: L1B_filenames
    CHARACTER( LEN=2048 ) :: L1B_swathlist
    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: IRR_filename, UV_filename
    CHARACTER (LEN=2048 ) :: IRR_swathname, UV_swathname
    CHARACTER( LEN = 20  ) :: FUNCTIONNAME = "O3T_mainNVAdj"
    CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg, CloudPressureSource
    INTEGER, EXTERNAL :: OMI_pixGetCldPres
    INTEGER, EXTERNAL :: OMI_pixGetSnowIce
    INTEGER, EXTERNAL :: OMI_pixGetTerPres

    CHARACTER( LEN = PGSd_PC_FILE_PATH_MAX ) :: SnowIceSourceOP

    REAL(KIND = 4), DIMENSION(:), ALLOCATABLE :: wl_cor
                                              !! wavelengths from corr.table
    REAL(KIND = 4), DIMENSION(:,:), ALLOCATABLE :: swpcr  !! corrections
    INTEGER(KIND = 4), DIMENSION(:), ALLOCATABLE :: iwlSub
    INTEGER(KIND = 4) :: il, ih
    REAL(KIND = 4) :: frac
    REAL(KIND=4), PARAMETER :: EPSILON10=1.0E-4
    !CHARACTER(LEN=128) :: filename

    logical :: use_tio = .true.
    type (l1b_tio_type) :: rad_file_obj
    integer :: errstat
    errstat = 0

    !! test snowice source option
    status = PGS_PC_GetConfigData( SNOWICESOURCE_LUN, msg )
    IF( status /= PGS_S_SUCCESS ) THEN    ! default ="Climatology"
       SnowIceSourceOP = '"Climatology"'
    ELSE
       READ( msg, '(A)'), SnowIceSourceOP
       IF( TRIM(SnowIceSourceOP) /= '"NISE"' .AND. &
           TRIM(SnowIceSourceOP) /= '"Climatology"' ) THEN
          SnowIceSourceOP = '"Climatology"'
       ENDIF
    ENDIF

    status = O3T_nval_setup( OMTO3_NVAL_LUN, nvRRS  )
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "Read NVAL LUT failed, PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF
    status = O3T_getLambdaSet()

    !! read in DNDX table
    status = O3T_dndx_setup( wl_cutoff )
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "Read DNDX LUT failed, PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    !! calculate common part of the coefficient for interpolation
    status = O3T_lpoly_cden( LUT_cden_blk )
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "Setup denominator for coef, " // &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    status = L1B_getNames( L1B_IRR_FILE_LUN, numfiles, L1B_filenames, &
                           L1B_swathlist )
    IF( status /= OZT_S_SUCCESS ) THEN
       status = L1B_getNames( BACKUP_L1BIRR_LUN, numfiles, L1B_filenames, &
                           L1B_swathlist )
    ENDIF

    IF( numfiles /= 1 ) THEN
       WRITE( msg,'(A)' ) "Number of IRR file not equal to 1, " // &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    !! Select Irradiance file
    status = L1B_selectIRR( "UV-2", USED_L1BIRR_LUN, &
                            IRR_filename, IRR_swathname, &
                            IRR_FILE_TYPE, NORMAL_L1BIRR_MISSING,  &
                            irr_error, irr_warning, irr_any )
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "Get IRR file or Swath Name failed, "// &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF
    IF( irr_error > 0 ) THEN
       WRITE( msg,'(A)' ) "One or more MeasurementQualityflags was set to "// &
             "Eorror in irradiance file:"// TRIM(IRR_filename) // &
             ". Alternative file should be staged. PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    !! read the L1B irradiance file
    if (use_tio) then
      call o3t_tio_getirr (IRR_filename, IRR_swathname, errstat)
      if (errstat < 0) stop 1
    else
      status = O3T_getIRR(IRR_filename, IRR_swathname)
      IF( status /= OZT_S_SUCCESS ) THEN
        WRITE( msg,'(A)' ) "O3T_getIRR failed, PGE aborting, exit code = 1."
        ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
        CALL EXIT(1)
      END IF
    endif

    !! open the L1B radiance file
    status = L1B_getNames( L1B_UV_FILE_LUN, numfiles, L1B_filenames, &
                           L1B_swathlist )
    IF( numfiles /= 1 ) THEN
       WRITE( msg,'(A)' ) "Number of L1BRUG file not equal to 1, " // &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ELSE
       UV_filename = L1B_filenames(1)
    ENDIF

    IF( INDEX( TRIM(L1B_swathlist), "Earth UV-2 Swath" ) == 0 ) THEN
       WRITE( msg,'(A)' ) "file:" // TRIM(UV_filename)//"," // &
                          TRIM(L1B_swathlist)//", does not have UV-2 Swath"//&
                          ", PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ELSE
       UV_swathname = "Earth UV-2 Swath"
    ENDIF

    if (use_tio) then
      call o3t_tio_initrad (rad_file_obj, UV_filename, UV_swathname, errstat)
      if (errstat < 0) stop 1
    else
      status = O3T_initRAD( UV_filename, UV_swathname )
      IF( status /= OZT_S_SUCCESS ) THEN
        WRITE( msg,'(A)' ) "O3T_initRAD failed, PGE aborting, exit code = 1."
        ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
        CALL EXIT(1)
      END IF
    endif

    !! Earth Sun distance Adjustment of IRR values to the date of radiance
    !! measurement because IRR and RUG measurment may be quite diferent in time.
    CALL O3T_AdjustIRREarthSun( EarthSundistance )

    !! setup storage for L2 output parameters, number of wavlengths in
    !! output is determined by the fixed output wavelength grid, while nXtrack,
    !! and nTimes are determined by input L1B.
    status = O3T_initL2out( wl_com )

    !! get the irradiance QA flags and precision at the fixed output
    !! wavelength grid.
    status = L1Bri_interpWL( irrWavelength, irrQAflags, wl_com, &
                             irrQAflags_com, irrPrecision, irrPrecision_com )

    !! read metadata from L1B irradiance file
    status = L1B_getCoreArchivedMetaData( USED_L1BIRR_LUN, &
                                          L1BIRRcoreArch, &
                                          year, month, day, jday )
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "L1B_getCoreArchivedMetaData IRR failed, " // &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    !! read metadata from L1B radiance file
    status = L1B_getCoreArchivedMetaData( L1B_UV_FILE_LUN, &
                                          L1BUV2coreArch, &
                                          year, month, day, jday )
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "L1B_getCoreArchivedMetaData UV2 failed, " // &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    !! Process Global Mode data only, skip any zoom mode data
    !! Note: variables names used here are meant to be used
    !! somewhre else with specific meaning. But it is ok burrow these
    !! variables to be used here.
    ii = LEN_TRIM( L1BUV2coreArch%ShortName )
    IF( L1BUV2coreArch%ShortName(ii:ii) == 'Z' ) THEN

       WRITE( msg,'(A,I9,A)' ) "PGE skip zoom mode orbit", &
                              L1BUV2coreArch%orbitNumber,&
                              "PGE finishes normally, exit code = 0"
       ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, FUNCTIONNAME, zero )

       CALL EXIT(0)
    ENDIF

    ! DO iLine = 1, NL1Bpsa_MAX
    !   ii = LEN_TRIM( L1BUV2coreArch%L1BpsaNames(iLine) )
    !   IF( ii <= 0 ) EXIT
    !   IF( INDEX( L1BUV2coreArch%L1BpsaNames(iLine), &
    !              "InstrumentConfigurationIDs" ) > 0 ) THEN
    !      DO iX = 1, Nelm_MAX
    !        ii = LEN_TRIM( L1BUV2coreArch%L1BpsaValues(iLine,iX) )
    !        IF( ii <= 0 ) EXIT
    !        READ( L1BUV2coreArch%L1BpsaValues(iLine,iX), '(I)'), pixID
    !        IF( pixID == 0 .OR. pixID == 1 .OR. &
    !            pixID == 2 .OR. pixID == 7 ) THEN
    !          skipit = .FALSE.
    !        ELSE
    !          skipit = .TRUE.
    !          WRITE( msg,'(A,I,A)' ) "PGE skip orbit = ",  &
    !           L1BUV2coreArch%orbitNumber, "which has non global mode data," // &
    !                 " PGE finishes normally, exit code = 0"
    !          ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, FUNCTIONNAME, zero )
    !          CALL EXIT(0)
    !        ENDIF
    !      ENDDO
    !   ENDIF
    ! ENDDO

    status = O3T_initCLD( CloudPressureSource )

    ! Read anomaly flags
    status= GET_ANOMFLG_3(year,jday,anomflg_3)

    status = OmiNVCinfo()
    ALLOCATE(wl_cor(nWvc))
    ALLOCATE(iwlSub(nWvc))
    ALLOCATE(swpcr(nWvc,nXtc))
    status = OmiNvalueCorr( L1BUV2coreArch%orbitNumber, wl_cor, swpcr )

    !! read in climatology data
    status = O3T_apriori_rd()
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "Read climatology data failed, " // &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    !! setup coefficient for so2 index calculations
    satname = O3T_getSatName( L1BUV2coreArch%ShortName, OMTO3_NVAL_LUN )
    status = O3_so2_setCoef( satname, o3abs, so2abs, iso2w, wl_com, soilim )
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "Set coefficient for SOI failed, " // &
                          "PGE aborting, exit code = 1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    !! Create the L2 output file
    status = L2_createFile( OMTO3_L2_LUN, OMTO3_fn )

    !! Create or setup the swath in the L2 output file
    dimList = "nWavel,nXtrack,nTimes,nLayers,nTimesSmallPixel"
    dims(1:5) = (/ nwl_com,nXtrack_rad,nTimes_rad,NLYR,nTimesSmallPixel_rad /)
    status = L2_setupSwath( OMTO3_fn, dimList, dims, OMTO3_SWATH_LUN, &
                            OMTO3_swathname, OMTO3_fileid, OMTO3_swid )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "OMI_setupL2swath failed" // &
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    !! Write retrieval wavelengths to L2 output
    df_omto3wl(:) = (/ df_Wavelength /)            !nWavel
    status = L2_defSWdatafields( OMTO3_swid, df_omto3wl )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_defSWdatafields failed," // &
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    status = L2_newBlockW( wlblk, OMTO3_fn, OMTO3_swathname, &
                           OMTO3_fileid, OMTO3_swid, df_omto3wl, nwl_com )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_newBlockW failed,"//&
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    wlblk%data(:) = TRANSFER( wl_com(1:nwl_com),  I1, &
                              nwl_com*wlblk%elmSize(1) )

    !! Write soft-calibration adjustments to L2 output
    df_caladj(:) = (/ df_CalibrationAdjustment /) !nWavel,nXtrack
    status = L2_defSWdatafields( OMTO3_swid, df_caladj )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_defSWdatafields failed," // &
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    status = L2_newBlockW( calblk, OMTO3_fn, OMTO3_swathname, &
                           OMTO3_fileid, OMTO3_swid, df_caladj, nxtrack )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_newBlockW failed,"//&
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    calblk%data(:) = TRANSFER( swpcr(:,:),  I1, &
         nwl_com*nxtrack*calblk%elmSize(1) )

    DO iwl = 1, nWvc
      il = iwl
      il = hunt( wl_com(1:nwl_com), wl_cor(iwl), il )
      IF( il == 0 .OR. il == nwl_com ) THEN
         !! wl_cor outside wl_com range
         ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "wl_cor outside wl_com range,"//&
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
         CALL EXIT(1)
      ELSE
         ih = il+1
         frac = (wl_cor(iwl) - wl_com(il))/ (wl_com(ih) - wl_com(il) )
         IF( ABS( frac ) < EPSILON10 ) THEN
            ih = il
         ELSE IF( ABS( frac -1.0 ) < EPSILON10 ) THEN
            il = ih
         ENDIF

         IF( ih /= il ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
                               "wl_cor not mathced any in wl_com,"//&
                               "PGE aborting, exit code = 1.", &
                               FUNCTIONNAME, zero )
            CALL EXIT(1)
         ENDIF
      ENDIF
      iwlSub(iwl) = ih
    ENDDO

    iLine_b = 0
    iLine_s = iLine_b

    status = L2_writeBlock( wlblk, iLine_s, nwl_com )
    CALL L2_disposeBlockW( wlblk   )
    status = L2_writeBlock( calblk, iLine_s, 720 )
    CALL L2_disposeBlockW( calblk  )

    IF( TRIM(L1BUV2coreArch%ShortName) == "OML1BRUG" ) THEN
       df_copy(:) = (/ df_NumberSmallPixelColumns,  &
                       df_SmallPixelColumn       ,  &
                       df_InstrumentConfigurationId /)
       status = L2_defSWdatafields( OMTO3_swid, df_copy )
       IF( status /= OZT_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_defSWdatafields failed,"//&
                                 "PGE aborting, exit code = 1.", &
                                 FUNCTIONNAME, zero )
          CALL EXIT(1)
       END IF
    ENDIF

    !! now starting defining the geo fields and data fields with one
    !! "nTime" as one of its dimensions.
    gf_omto3(:) = (/ gf_GroundPixelQualityFlags, &
                     gf_Latitude               , &
                     gf_Longitude              , &
                     gf_SolarZenithAngle       , &
                     gf_SolarAzimuthAngle      , &
                     gf_ViewingZenithAngle     , &
                     gf_ViewingAzimuthAngle    , &
                     gf_RelativeAzimuthAngle   , &
                     gf_TerrainHeight          , &
                     gf_Time                   , &
                     gf_SecondsInDay           , &
                     gf_SpacecraftLatitude     , &
                     gf_SpacecraftLongitude    , &
                     gf_SpacecraftAltitude     , &
                     gf_XTrackQualityFlags     /)

    df_omto3(:) = (/ df_CloudPressure,           &
                     df_TerrainPressure,         &
                     df_AlgorithmFlags,          &
                     df_QualityFlags,            &
                     df_RadianceBadPixelFlagAccepted,&
                     df_fc,                      &
                     df_CloudFraction,           &
                     df_ColumnAmountO3,          &
                     df_O3BelowCloud,            &
                     df_Reflectivity331,         &
                     df_Reflectivity360,         &
                     df_SO2index,                &
                     df_StepOneO3,               &
                     df_StepTwoO3,               &
                     df_UVAerosolIndex,          &
                     df_dNdR,                    &
                     df_NValue,                  &
                     df_Residual,                &
                     df_ResidualStep1,           &
                     df_ResidualStep2,           &
                     df_Sensitivity,             &
                     df_dNdT,                    &
                     df_APrioriLayerO3,          &
                     df_LayerEfficiency,         &
                     df_MeasurementQualityFlags /)

    !! define the geofield in the L2 output file
    status = L2_defSWgeofields( OMTO3_swid, gf_omto3 )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_defSWgeofields failed," // &
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    !! define the datafield in the L2 output file
    status = L2_defSWdatafields( OMTO3_swid, df_omto3 )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_defSWdatafields failed," // &
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    !! Allocate the memory for storing the L2 output information for a block.
    !! The memory allocated before for the wl_com is deallocated and
    !! reallocated for the new df_omto3 data fields.
    status = L2_newBlockW( geoblk, OMTO3_fn, OMTO3_swathname, &
                           OMTO3_fileid, OMTO3_swid, gf_omto3, &
                           nLinesPerWrite )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_newBlockW for Geo failed,"//&
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    status = L2_newBlockW( datablk, OMTO3_fn, OMTO3_swathname, &
                           OMTO3_fileid, OMTO3_swid, df_omto3, &
                           nLinesPerWrite )
    IF( status /= OZT_S_SUCCESS ) THEN
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, "L2_newBlockW for Data failed,"//&
                              "PGE aborting, exit code = 1.", &
                              FUNCTIONNAME, zero )
       CALL EXIT(1)
    END IF

    !! ------- now begin processing -------

    !! determine descendQ flag in the beginning of the granule
    IF( MOD( nXtrack_rad, 2 ) == 0 ) THEN
       nXevenQ = .TRUE.
    ELSE
       nXevenQ = .FALSE.
    ENDIF

    iLine = 0
    if (use_tio) then
      call l1b_tio_getgeo (rad_file_obj, radgeo_blk, iLine, latitude, longitude, &
                           szenith, sazimuth, vzenith, vazimuth, geoflg, &
                           errstat)
      if (errstat < 0) stop 1
    else
      status = L1Bga_getLine( geo_blk, iLine, time, &
                             latitude, longitude, &
                             szenith, sazimuth, vzenith, vazimuth, &
                             height, geoflg )
    endif

    IF( nXevenQ ) THEN
       latPre = (latitude(nXtrack_rad/2-1) + latitude(nXtrack_rad/2))*0.5
    ELSE
       latPre = latitude((nXtrack_rad-1)/2)
    ENDIF

    iLine = 1
    if (use_tio) then
      call l1b_tio_getgeo (rad_file_obj, radgeo_blk, iLine, latitude, longitude, &
                           szenith, sazimuth, vzenith, vazimuth, geoflg, &
                           errstat)
      if (errstat < 0) stop 1
    else
      status = L1Bga_getLine( geo_blk, iLine, time, &
                             latitude, longitude, &
                             szenith, sazimuth, vzenith, vazimuth, &
                             height, geoflg )
    endif

    IF( nXevenQ ) THEN
       latNow = (latitude(nXtrack_rad/2-1) + latitude(nXtrack_rad/2))*0.5
    ELSE
       latNow = latitude((nXtrack_rad-1)/2)
    ENDIF
    descendQ = O3T_descendQ( latPre, latNow )

    L2_parameters(1)%NumberOfInputSamples          = nTimes_rad*nXtrack_rad
    L2_parameters(1)%NumberOfGoodInputSamples      = 0
    L2_parameters(1)%NumberOfLargeSZAInputSamples  = 0
    L2_parameters(1)%NumberOfMissingInputSamples   = 0
    L2_parameters(1)%NumberOfBadInputSamples       = 0
    L2_parameters(1)%NumberOfInputWarningSamples   = 0
    L2_parameters(1)%NumberOfGoodOutputSamples     = 0
    L2_parameters(1)%NumberOfGlintCorrectedSamples = 0
    L2_parameters(1)%NumberOfSkippedSamples        = 0
    L2_parameters(1)%NumberOfStep1InvalidSamples   = 0
    L2_parameters(1)%NumberOfStep2InvalidSamples   = 0
    L2_parameters(1)%NumberOfIrradianceMissing     = 0
    L2_parameters(1)%NumberOfIrradianceError       = 0
    L2_parameters(1)%NumberOfIrradianceWarning     = 0
    L2_parameters(1)%NumberOfRadianceMissing       = 0
    L2_parameters(1)%NumberOfRadianceError         = 0
    L2_parameters(1)%NumberOfRadianceWarning       = 0
    L2_parameters(1)%NumberOfMeasurement           = nTimes_rad
    L2_parameters(1)%NumberOfMeasurementMissing    = 0
    L2_parameters(1)%NumberOfMeasurementError      = 0
    L2_parameters(1)%NumberOfMeasurementWarning    = 0
    L2_parameters(1)%NumberOfMeasurementRebinned   = 0
    L2_parameters(1)%NumberOfMeasurementSAA        = 0
    L2_parameters(1)%NumberOfMeasurementManeuver   = 0
    L2_parameters(1)%NumberOfInstrumentSettingsError= 0
    L2_parameters(1)%ZonalOzoneMin(1:5)            = 1000.0
    L2_parameters(1)%ZonalOzoneMax(1:5)            =-1000.0
    L2_parameters(1)%ZonalLatRange(1,1:5) = (/ -90, -60, -30, 30, 60 /)
    L2_parameters(1)%ZonalLatRange(2,1:5) = (/ -60, -30,  30, 60, 90 /)

    !! processing one line at a time.
    iLine_s = iLine_b
    iLine_e = nTimes_rad-1
    DO iLine = iLine_b, iLine_e

      if (use_tio) then
        call l1b_tio_getgeo (rad_file_obj, radgeo_blk, iLine, latitude, longitude, &
                             szenith, sazimuth, vzenith, vazimuth, geoflg, &
                             errstat, anomflg)
        if (errstat < 0) stop 1
      else
        status = L1Bga_getLine( geo_blk, iLine, time, &
                               latitude, longitude, &
                               szenith, sazimuth, vzenith, vazimuth, &
                               height, geoflg, anomflg, sInD, scLat, scLon, scHgt )
        IF( status /=  OZT_S_SUCCESS ) THEN
          WRITE( msg,'(A)' ) "L1Bga_getLine failed, " // &
            "PGE aborting, exit code = 1."
          ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
          CALL EXIT(1)
        END IF
      endif

      phiArray(1:nXtrack_rad ) =  &
         (/ ( adjustDEG(180.0+sazimuth(iX)-vazimuth(iX)), iX =1,nXtrack_rad ) /)

      DO iX = 1, nXtrack_rad
        status = OMI_pixGetTerPres( latitude(iX), longitude(iX), ptArray(iX) )
        status = OMI_pixGetSnowIce( latitude(iX), longitude(iX), &
                                          year, month, day, snowIceArray(iX) )
      ENDDO

      IF( TRIM(CloudPressureSource) == '"Climatology"' ) THEN
         DO iX = 1, nXtrack_rad
           status = OMI_pixGetCldPres( latitude(iX), longitude(iX), &
                                               year, month, day, pcArray(iX) )
           !! Make sure cloud pressure is inside the range 200.0 to 1013.25 hPa.
           PclimQ(iX) = .TRUE.
           IF( pcArray(iX) < Pc_min      ) THEN
              pcArray(iX) = Pc_min
           ELSE IF( pcArray(iX) > Pc_max ) THEN
              pcArray(iX) = Pc_max
           ENDIF
           !! Do not allow cloud pressure lower than the terrain pressure.
           IF( ptArray(iX) < pcArray(iX)) pcArray(iX) = ptArray(iX)
         ENDDO
      ELSE
         CALL O3T_getOMICldPress( iLine, pcArray, nXtrack_rad )
         PclimQ(1:nXtrack_rad) = .FALSE.
         DO iX = 1, nXtrack_rad
           IF( TRIM(CloudPressureSource) == '"OMCLDO2"' ) THEN
              cld_errflg = IBITS( ProcessingQualityFlags(iX),5,7 )
           ELSE IF( TRIM(CloudPressureSource) == '"OMCLDRR"' ) THEN
              cld_errflg = IBITS( ProcessingQualityFlags(iX),0, 3 ) &
                         + IBITS( ProcessingQualityFlags(iX),6, 2 ) &
                         + IBITS( ProcessingQualityFlags(iX),14,2 )
           ENDIF
           !! if cloud pressure not good, either by its quality flag
           !! or by the cloud fraction, use climatology cloud pressure instead.

           !!  IF( cld_errflg > 0 .OR. CloudFraction(iX) < 0.15 ) THEN
           IF( cld_errflg > 0 ) THEN
              status = OMI_pixGetCldPres( latitude(iX), longitude(iX), &
                                          year, month, day, pcArray(iX) )
              PclimQ(iX) = .TRUE.
           ENDIF
           !! Make sure cloud pressure is inside the range 250.0 to 1013.25 hPa.
           IF( pcArray(iX) < Pc_min      ) THEN
              pcArray(iX) = Pc_min
           ELSE IF( pcArray(iX) > Pc_max ) THEN
              pcArray(iX) = Pc_max
           ENDIF
           !! Do not allow cloud pressure lower than the terrain pressure.
           IF( ptArray(iX) < pcArray(iX)) pcArray(iX) = ptArray(iX)
         ENDDO
      ENDIF

      iT = MOD( iLine, nLinesPerWrite ) + 1
      !! All geolocation information read from the input L1B is copied
      !! to the output file, no modification is done to geo info before
      !! written out the L2 output.
      CALL O3T_L2setGeoLine( iT, geoblk, datablk )

      if (use_tio) then
        call l1b_tio_getrad (rad_file_obj, radgeo_blk, iLine, errstat, &
                             radiance, radWavelength, radQAflags)
        if (errstat < 0) stop 1
      else
        status = L1Bri_getLine( rad_blk, iLine, &
                               RadIrr_k = radiance, &
                               RadIrrPrecision_k = radPrecision, &
                               PixelQualityFlags_k = radQAflags, &
                               Wavelength_k = radWavelength )
        IF( status /=  OZT_S_SUCCESS ) THEN
          WRITE( msg,'(A)' ) "L1Bri_getLineWL failed, " // &
            "PGE aborting, exit code = 1."
          ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
          CALL EXIT(1)
        END IF
      endif

      !! get radiance QA flags and precision at the fixed output
      !! wavelength grid.
      status = L1Bri_interpWL( radWavelength, radQAflags, wl_com, &
                               radQAflags_com, radPrecision, radPrecision_com )

      pixSURF%iwl_glint = iwl_refl_h
      pixSURF%nres_glnt_ub = glntlm        !! defined in O3T_const file

      IF( iLine > 0 ) THEN
         IF( nXevenQ ) THEN
            latNow = (latitude(nXtrack_rad/2-1) + latitude(nXtrack_rad/2))*0.5
         ELSE
            latNow = latitude((nXtrack_rad-1)/2)
         ENDIF
         descendQ = O3T_descendQ( latPre, latNow )
         latPre = latNow
      ENDIF

      iwlArray = (/ iwl_mix, iwl_ozone, iwl_refl_l, iwl_refl_h /)

      DO iX = 1, nXtrack_rad
        iwl_oz             = iwl_ozone
        skipit             = .FALSE.
        pixSURF%glint_flag = .FALSE.
        sza_p              = szenith(iX)
        vza_p              = vzenith(iX)
        phi                = phiArray(iX)
        lat                = latitude(iX)
        lon                = longitude(iX)
        hgt                = height(iX)
        pt                 = ptArray(iX)/Pc_max
        pc                 = pcArray(iX)/Pc_max
        ilat               = (int(lat)+90)/5+1

        IF (ilat > 36) ilat=36

        IF( TRIM(SnowIceSourceOP) == '"NISE"' ) THEN
           nise_flag  = IBITS( geoflg(iX), 8, 7 )
           IF( nise_flag >= 80 .AND. nise_flag <= 103  )THEN
              pixSURF%isnow = 10
           ELSE
              pixSURF%isnow = 0
           ENDIF
        ELSE IF( TRIM(SnowIceSourceOP) == '"Climatology"' ) THEN
           IF( snowIceArray(iX) >= 95.0 ) THEN
              pixSURF%isnow = 10
           ELSE
              pixSURF%isnow = 0
           ENDIF
        ENDIF

        ! check instrument configuration id
        status = L1Bri_getInstConfigId(  rad_blk, iLine, instId_rad )

        QAflags = 0
        radBadPixflgs = 0
        skipit = ( instid_rad /= 0 .AND. &
                   instid_rad /= 1 .AND. &
                   instid_rad /= 2 .AND. &
                   instid_rad /= 7 ) .OR. &
                   O3T_setQAflgsI( QAflags, radBadPixflgs, anomflg(iX), anomflg_3(iX,ilat), algflg, &
                                   L2_parameters(1), descendQ, PclimQ(iX), &
                                   sza_p, geoflg(iX), iwlArray, &
                                   radQAflags_com(:,IX), &
                                   irrQAflags_com(:,IX), &
                                   radPrecision_com(:,IX), &
                                   irrPrecision_com(:,IX) )

        !! calculates xnvalm no matter what QAflags is set for the pixel.
        status = O3T_nvalm( irrWavelength(:,iX), irradiance(:,iX),  &
                            radWavelength(:,iX),   radiance(:,iX),  &
                            wl_com(:), xnvalm(:) )
        IF( skipit .OR. status /= OZT_S_SUCCESS  ) THEN
           IF( pixSURF%isnow == 10 ) algflg = algflg + 10
           errflgs = 7
           QAflags = QAflags + errflgs
           CALL O3T_L2fillDataPix( iT, iX, nwl_com, NLYR, algflg, QAflags, &
                                   radBadPixflgs, datablk )
           CALL O3T_QAcounter( algflg, QAflags, &
                               L2_parameters(1)%QualityFlagsCounters )
           CYCLE
           IF( status /= OZT_S_SUCCESS ) THEN
              WRITE( msg,'(A,2(A,I4),A)' ) "Calculate xnalm failed at ",  &
                     "iLine =", iLine, " iX=", iX, "Skip to next pixel."
              ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
           ENDIF
        ENDIF
!!!
        xnvalm(iwlSub(1:nWvc)) = xnvalm(iwlSub(1:nWvc)) + swpcr(1:nWvc, IX)/100.
!!!

        pixID  = iX + nXtrack_rad*iLine
        pixGEO = O3T_pixgeo_set( sza_p, vza_p, phi, pt, pc, lat, lon, pixID )

        status = O3T_lpoly_coef( coefs, LUT_cden_blk, pixGEO )
        IF( status /= OZT_S_SUCCESS ) THEN
           WRITE( msg,'(A)' ) "Calculate poly coefficient failed, " // &
                              "Skip to next pixel"
           ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
           IF( pixSURF%isnow == 10 ) algflg = algflg + 10
           errflgs = 7
           QAflags = QAflags + errflgs
           CALL O3T_L2fillDataPix( iT, iX, nwl_com, NLYR, algflg, QAflags, &
                                   radBadPixflgs, datablk )
           CALL O3T_QAcounter( algflg, QAflags, &
                               L2_parameters(1)%QualityFlagsCounters )
           CYCLE
        ENDIF

        pixSURF%landsea_mask = IBITS( geoflg(iX), 0, 4 )
        IF( pixSURF%landsea_mask /= 1 .AND. pixSURF%landsea_mask /= 3 ) THEN
           pixSURF%surface_category = 0
        ELSE
           pixSURF%surface_category = 1
        ENDIF

        IF( iX == 1 ) THEN
           IF( ABS(latitude(iX) ) <=  45.0 ) THEN
              guesoz = 260
           ELSE IF(ABS(latitude(iX)) <= 75.0 ) THEN
              guesoz = 340
           ELSE
              guesoz = 360
           ENDIF
        ELSE
           guesoz = stp3oz
        ENDIF
        status = O3T_step1( iwl_oz, iwl_refl_l, iwl_refl_h, xnvalm, pixGEO, &
                           pixSURF, coefs, iwl_refl, iplow, guesoz, stp1oz, &
                           res_stp1, dndomega_t, dndr, algflg, absrfl,maxitr, &
                           stp1oz_valid, skipit, doz_limit )
        IF( skipit .OR. status /= OZT_S_SUCCESS ) THEN
           WRITE( msg,'(A,I3,I5, 2(A9,F9.3), 3(A9,L1) )' ) &
                     "O3T_step1 skipped at (iX, iL)=", iX, iLine, &
                     ",stp1oz=", stp1oz, ", refl=", pixSURF%ref, ",maxitr=", maxitr, &
                     ",stp1oz=", stp1oz_valid, ",skipit=", skipit
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, FUNCTIONNAME, one )
           IF( maxitr ) THEN
              errflgs = 6
           ELSE
              errflgs = 7
           ENDIF
           QAflags = QAflags + errflgs
           IF( .NOT. stp1oz_valid ) THEN
              L2_parameters(1)%NumberOfStep1InvalidSamples = &
              L2_parameters(1)%NumberOfStep1InvalidSamples + 1
           ENDIF
           CALL O3T_L2fillDataPix( iT, iX, nwl_com, NLYR, algflg, QAflags, &
                                   radBadPixflgs, datablk )
           CALL O3T_QAcounter( algflg, QAflags, &
                               L2_parameters(1)%QualityFlagsCounters )
           CYCLE
        ENDIF

        status = O3T_step2( iwl_oz, iwl_refl, iplow, jday, pixGEO, pixSURF, &
                           coefs, stp1oz, res_stp1, dndomega_t, dndr,absrfl, &
                           stp2oz, res_stp2, dr, stp2prf, eff, aprfoz, &
                           stp2oz_valid, skipit, dNdT )

        IF( skipit .OR. status /= OZT_S_SUCCESS ) THEN
           WRITE( msg,'(A,I3,I5,2(A9,G10.3),A8,L1)' ) &
                     "O3T_step2 skipped at (iX, iL)=(", iX, iLine, &
                     "),stp1oz=", stp1oz, ",stp2oz=", stp2oz,",skipit=", skipit
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, FUNCTIONNAME, one )
           IF( .NOT. stp2oz_valid ) THEN
              L2_parameters(1)%NumberOfStep2InvalidSamples = &
              L2_parameters(1)%NumberOfStep2InvalidSamples + 1
           ENDIF
           errflgs = 7
           QAflags = QAflags + errflgs
           CALL O3T_L2fillDataPix( iT, iX, nwl_com, NLYR, algflg, QAflags, &
                                   radBadPixflgs, datablk )
           CALL O3T_QAcounter( algflg, QAflags, &
                               L2_parameters(1)%QualityFlagsCounters )
           CYCLE
        ENDIF

        !  set aerosol index
        aerind = fill_float32
        IF( iwl_refl /= iwl_refl_h ) aerind = -res_stp2(iwl_refl_h)

        !! compute step 3 ozone
        status = O3T_step3( iwl_refl_h, iwl_refl_l, iwl_mix, pixGEO, stp2oz, &
                            f313, f360, res_stp2, dndomega_t, &
                            stp3oz, res_stp3, pathl, algflg )
        IF( status /= OZT_S_SUCCESS ) THEN
           errflgs = 7
           QAflags = QAflags + errflgs
           CALL O3T_L2fillDataPix( iT, iX, nwl_com, NLYR, algflg, QAflags, &
                                   radBadPixflgs, datablk )
           CALL O3T_QAcounter( algflg, QAflags, &
                               L2_parameters(1)%QualityFlagsCounters )
           CYCLE
        ENDIF

        ! calculate reflectance at wl_refl_h
        pixSURF%ref360 = pixSURF%ref &
                       + (res_stp3(iwl_refl_h)-res_stp3(iwl_refl_l)) &
                       / dndr(iwl_refl_h)

        ! calculate SO2 index
        IF( pathl < 3.5 ) THEN
           so2ind = O3_so2_index( wl_com, dndomega_t, res_stp3, &
                                  iso2w, o3abs, so2abs )
        ELSE
           so2ind = fill_float32
        ENDIF

        IF( pixSURF%isnow == 10 ) algflg = algflg + 10

        ! calculate Ozone below fractional cloud
        oz_cld = O3T_blwcld( stp2prf, pixGEO, pixSURF )

        ! set processing part of the QA flags: errflgs
        errflgs = O3T_setQAflgsP( iwl_refl_l, iwl_refl_h, iwl_oz, iwl_mix, &
                                  pixSURF%glint_flag, maxitr, algflg, &
                                  res_stp3, flg3lm, flg4lm, soilim, &
                                  sza_p, pathl, so2ind, stp3oz )
        ! Total QA flags is the summation of the input QA flags and the
        ! processing error flags
        QAflags = QAflags + errflgs
        CALL O3T_QAcounter( algflg, QAflags, &
                            L2_parameters(1)%QualityFlagsCounters )
        IF( IBITS( QAflags, 0, 4 ) == 0 ) THEN
           L2_parameters(1)%NumberOfGoodOutputSamples = &
           L2_parameters(1)%NumberOfGoodOutputSamples + 1
        ELSE IF( IBITS( QAflags, 0, 4 ) == 1 ) THEN
           L2_parameters(1)%NumberOfGlintCorrectedSamples = &
           L2_parameters(1)%NumberOfGlintCorrectedSamples + 1
        ENDIF

        !! store data in data fields output block
        CALL O3T_L2setDataPix( iT, iX, nwl_com, NLYR, algflg, QAflags, radBadPixflgs, &
                               stp1oz, stp2oz, stp3oz, oz_cld, aerind, so2ind,&
                               pixSURF, eff, aprfoz, datablk )
        CALL O3T_zonalMinMax( QAflags, lat, stp3oz, L2_parameters(1) )

      ENDDO ! iX end loop for Cross Track index

      status = L1Bri_getMeasurementQF( rad_blk, iLine, mqa_rad    )

      IF( nXtrack_rad == &
          O3T_nMissing( iwlArray, radQAflags(:,:) ) ) THEN
         radLmissing = .TRUE.
      ELSE
         radLmissing = .FALSE.
      ENDIF

      instIDmismatch = O3T_instIDsCheck( TRIM(L1BUV2coreArch%ShortName), &
                                         mqa_rad, instID_rad, instID_irr )
      mqaL2          = O3T_setMqaL2( mqa_rad, radLMissing, instIDmismatch, &
                                     bit7Q, L2_parameters(1) )
      CALL O3T_L2setMqaLine( iT, mqaL2, datablk )

      IF( iT == nLinesPerWrite ) THEN
         nLw = nLinesPerWrite
         IF( nLw > (iLine - iLine_s + 1 ) ) THEN
            nLw = iLine - iLine_s + 1
         ENDIF
         status = L2_writeBlock( geoblk,  iLine_s, nLw )
         status = L2_writeBlock( datablk, iLine_s, nLw )
         iLine_s = iLine_s + nLw
         iT = 0 ! When nTimes = n*nLinesPerWrite (n is an integer), iT =
                ! nLinesPerWrite at the end of the loop (not equal to zero).
                ! This would induce the the last write outside of the iLine
                ! loop. Setting iT = 0 would avoid this extra write.
      ENDIF
    ENDDO  ! iLine end loop for Along Track index

    !! write the last blocks of data (if any)
    IF( iT /= 0 ) THEN
       nLw = iT
       status = L2_writeBlock( geoblk,  iLine_s, nLw )
       status = L2_writeBlock( datablk, iLine_s, nLw )
    ENDIF

    !! writing finished, close the L2 file
    CALL L2_disposeBlockW( geoblk  )
    CALL L2_disposeBlockW( datablk )
    ierr = HE5_SWdetach( OMTO3_swid )
    ierr = HE5_SWclose( OMTO3_fileid )

    !! free the allocated memory
    CALL O3T_freeIRR
    CALL O3T_freeRAD
    CALL O3T_freeCLD
    CALL O3T_freeL2out

    !! free coefficient denominator data block
    status = O3T_lpolycden_dispose( LUT_cden_blk )
    !! free the dndx tables.
    status = O3T_dndx_dispose()
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "O3T_dndx_dispose failed, PGE aborting, exit code=1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    !! clean up space allocated for NVAL LUT
    status = O3T_nval_dispose()
    IF( status /= OZT_S_SUCCESS ) THEN
       WRITE( msg,'(A)' ) "O3T_nval_dispose failed, PGE aborting, exit code=1."
       ierr = OMI_SMF_setmsg( OZT_E_FAILURE, msg, FUNCTIONNAME, zero )
       CALL EXIT(1)
    ENDIF

    ! free memory used for anomaly flag
    DEALLOCATE(anomflg_3)

    L2_parameters(1)%Name = "Total Column Ozone"

    IF( NORMAL_L1BIRR_MISSING ) THEN
      L2_parameters(1)%SolarProductMissing = 1
    ELSE
      L2_parameters(1)%SolarProductMissing = 0
    ENDIF

    IF( IRR_FILE_TYPE == "Backup" ) THEN
       L2_parameters(1)%BackupSolarProductUsed = 1
    ELSE
       L2_parameters(1)%BackupSolarProductUsed = 0
    ENDIF

    IF( TRIM(L1BUV2coreArch%ShortName) == "OML1BRUZ" ) THEN
       mcfLUN = OMTO3Z_MCF_LUN
    ELSE
       mcfLUN = OMTO3_MCF_LUN
    ENDIF

    !! Read the input file names from the PCF file and
    !! use these files as the input pointer for the L2 output
    IF( TRIM(CloudPressureSource) == '"Climatology"' ) THEN
       LUNinputPointer(1:10) = (/ L1B_UV_FILE_LUN, USED_L1BIRR_LUN,  &
                                 O3_CLIM_LUN,     TM_CLIM_LUN,      &
                                 TERRAINPRES_LUN, CLOUDPRES_LUN,    &
                                 OMTO3_NVAL_LUN,  OMTO3_DNDX_LUN,   &
                                 nvCORR_LUN, ANOMFLG3_LUN /)
       status = L2_setCoreArchMetaData( OMTO3_L2_LUN , L1BUV2coreArch, &
                                 L1BIRRcoreArch, LUNinputPointer(1:10), &
                                 mcfLUN, L2_parameters )
    ELSE
       LUNinputPointer(1:11)= (/ L1B_UV_FILE_LUN, USED_L1BIRR_LUN,  &
                                O3_CLIM_LUN,     TM_CLIM_LUN,       &
                                TERRAINPRES_LUN, CLOUDPRES_LUN,     &
                                OMCLDRR_L2_LUN, OMTO3_NVAL_LUN,     &
                                OMTO3_DNDX_LUN,    nvCORR_LUN, ANOMFLG3_LUN /)
       status = L2_setCoreArchMetaData( OMTO3_L2_LUN , L1BUV2coreArch, &
                                L1BIRRcoreArch, LUNinputPointer(1:11), &
                                mcfLUN, L2_parameters )
    ENDIF

    !! Write global attribute
    status = OMI_writeGlobalAttribute( OMTO3_fn, year, month, day, &
                                       LUNinputPointer(1:2) )

    !! Write Swath attribute
    status = OMI_writeSwathAttribute( OMTO3_fn, OMTO3_swathname, &
                                      nTimes_rad, nTimesSmallPixel_rad, &
                                      EarthSundistance, "Total Column" )

    !! copy some L1B data fileds to L2 output, should be for OMI only.
    IF( TRIM(L1BUV2coreArch%ShortName) == "OML1BRUG" ) THEN
      msg = "NumberSmallPixelColumns,SmallPixelColumn," // &
            "InstrumentConfigurationId"
      status = OMI_readHE4fields( L1B_UV_FILE_LUN, "UV-2", msg )
      status = OMI_cpwtHE5fields( OMTO3_L2_LUN, "Column Amount O3" )
    ENDIF

    !! copy some HE4 swath attributes and write to HE5 global attributes
    !! mainly for TOMS processing.
    IF( INDEX( L1BUV2coreArch%ShortName, "TOMS" ) > 0 ) THEN
       status = OMI_readHE4swAttr( L1B_UV_FILE_LUN )
       status = OMI_cpwtHE5glAttr( OMTO3_L2_LUN    )
    ENDIF

    !! write the final message
    ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, &
                           "PGE finishes normally, exit code = 0",&
                           FUNCTIONNAME, zero )

    CALL EXIT(0)
END PROGRAM O3T_mainNVAdj
