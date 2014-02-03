MODULE omi_pge_fitting_aux

  use errormodule
  private
  public find_swathrange_by_latitude, read_latitude, &
    find_swathline_by_latitude, check_wavelength_overlap, convert_tai_to_utc, &
    find_swathline_range, &
    compute_fitting_statistics, compute_fitting_statistics_nohe5, &
    omi_set_xtrpix_range, &
    omi_set_fitting_parameters, set_input_pointer_and_versions

CONTAINS
  SUBROUTINE omi_set_fitting_parameters ( pge_idx, errstat )

    USE OMSAO_precision_module
    USE OMSAO_variables_module,  ONLY: l1b_channel
    USE OMSAO_he5_module,        ONLY: swath_base_name, pge_swath_name
    USE OMSAO_errstat_module,    ONLY: pge_errstat_ok
    USE OMSAO_omidata_module,    ONLY: omi_radiance_swathname, omi_irradiance_swathname
    USE OMSAO_indices_module,    ONLY: sao_molecule_names

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN)  :: pge_idx

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=26), PARAMETER :: modulename = 'omi_set_fitting_parameters'

    ! --------------------------
    ! Initialize OUTPUT variable
    ! --------------------------
    errstat = pge_errstat_ok

    ! ---------------------------------------------------------------------
    ! Name of solar, earthshine, and L2 swaths (normally obtained from PCF)
    ! ---------------------------------------------------------------------
    SELECT CASE ( l1b_channel )
    CASE ( 'UV1' )
      omi_radiance_swathname   = 'Earth UV-1 Swath'
      omi_irradiance_swathname = 'Sun Volume UV-1 Swath'
    CASE ( 'UV2' )
      omi_radiance_swathname   = 'Earth UV-2 Swath'
      omi_irradiance_swathname = 'Sun Volume UV-2 Swath'
    CASE ( 'VIS' )
      omi_radiance_swathname   = 'Earth VIS Swath '
      omi_irradiance_swathname = 'Sun Volume VIS Swath'
    END SELECT

    IF ( pge_swath_name(1:1) == "?" ) &
      pge_swath_name = &
      TRIM(ADJUSTL(swath_base_name(pge_idx)))//" "//TRIM(ADJUSTL(sao_molecule_names(pge_idx)))

    RETURN
  END SUBROUTINE omi_set_fitting_parameters

  SUBROUTINE compute_fitting_statistics ( &
      pge_idx, ntimes, nxtrack, xtrange, saocol, saodco, saorms, saofcf, saomqf, errstat )

    USE OMSAO_precision_module,  ONLY: i2, i4, r4, r8
    USE OMSAO_parameters_module, ONLY: &
      i2_missval, r8_missval, main_qa_good, main_qa_suspect, main_qa_bad
    use optimizer_interface_module, only: &
      opt_convergence_failed, opt_convergence_maxiter_exceeded, opt_convergence_suspect, &
      opt_convergence_good
    USE metadata_tools,  ONLY:  QAPercentMissingData, QAPercentOutofBoundsData
    USE OMSAO_he5_module,       ONLY:  &
      NrOfInputSamples, NrofGoodOutputSamples, NrofSuspectOutputSamples,        &
      NrofBadOutputSamples, NrofConvergedSamples, NrofFailedConvergenceSamples, &
      NrofExceededIterationsSamples, NrofOutofBoundsSamples, NrofMissingSamples, &
      NrofGoodInputSamples, NrofSuspectOutputSamples, NrofBadOutputSamples,      &
      NrofConvergedSamples, NrofFailedConvergenceSamples, &
      PercentGoodOutputSamples, PercentSuspectOutputSamples, &
      PercentBadOutputSamples, &
      AbsolutePercentMissingSamples
    USE OMSAO_errstat_module,   ONLY: vb_lev_screen, pge_errstat_ok
    USE OMSAO_variables_module, ONLY: verb_thresh_lev, max_good_col
    USE he5_output_tools, ONLY: he5_write_fitting_statistics
    USE metadata_tools, ONLY: set_automatic_quality_flag

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ntimes, nxtrack
    INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),     INTENT (IN) :: xtrange
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saocol, saodco, saorms
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saofcf

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saomqf

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Local variables
    ! ----------------
    INTEGER (KIND=i4) :: locerrstat, ix, it, spix, epix
    REAL    (KIND=r4) :: PercentOutofBoundsSamples
    REAL    (KIND=r8) :: fitcol_avg, rms_avg, dfitcol_avg, nfitcol
    REAL    (KIND=r8) :: col2sig, col3sig

    locerrstat = pge_errstat_ok

    ! ---------------------------------------------------------
    ! The total number of input samples is simply the number of
    ! pixels in the granule
    ! ---------------------------------------------------------
    NrofInputSamples = nxtrack*ntimes

    ! ------------------------------------------------------------------
    ! Compute all other fitting statistics variables over two nice loops
    ! ------------------------------------------------------------------
    saomqf                        = i2_missval
    NrofGoodInputSamples          = 0_i4
    NrofGoodOutputSamples         = 0_i4
    NrofSuspectOutputSamples      = 0_i4
    NrofBadOutputSamples          = 0_i4
    NrofOutOfBoundsSamples        = 0_i4
    NrofConvergedSamples          = 0_i4
    NrofFailedConvergenceSamples  = 0_i4
    NrofExceededIterationsSamples = 0_i4
    NrofMissingSamples            = 0_i4

    nfitcol    = 0.0_r8
    fitcol_avg = 0.0_r8 ; rms_avg = 0.0_r8 ; dfitcol_avg = 0.0_r8
    DO it = 0, ntimes-1

      spix = xtrange(it,1) ; epix = xtrange(it,2)
      DO ix = spix, epix

        col2sig = saocol(ix,it)+2.0_r8*saodco(ix,it)
        col3sig = saocol(ix,it)+3.0_r8*saodco(ix,it)

        ! ------------------------------------------------------
        ! The Good: Columns are postive within two sigma fitting
        !           uncertainty and the fitting has converged.
        !           For this "sweet spot" we compute the average
        !           fitting statistics.
        ! ------------------------------------------------------
        IF ( (saofcf(ix,it)   == opt_convergence_good) .AND. &
            (saocol(ix,it)      >  r8_missval                       ) .AND. &
            (ABS(saocol(ix,it)) <= max_good_col                     ) .AND. &
            (col2sig            >= 0.0_r8                           ) ) THEN

          saomqf(ix,it) = main_qa_good

          NrofGoodInputSamples  = NrofGoodInputSamples   + 1
          NrofConvergedSamples  = NrofConvergedSamples   + 1
          NrofGoodOutputSamples = NrofGoodOutputSamples + 1

          fitcol_avg  = fitcol_avg  + saocol(ix,it)
          dfitcol_avg = dfitcol_avg + saodco(ix,it)
          rms_avg     = rms_avg     + saorms(ix,it)
          nfitcol     = nfitcol     + 1.0_r8

          CYCLE
        END IF

        ! ----------------------------------------------------------
        ! The Bad: Fitting hasn't converged or columns are negative
        !          within three sigma fitting uncertainty. Note that
        !          pixels can count towards both the number of out-
        !          of bounds and the failed convergence samples.
        ! ----------------------------------------------------------
        IF ( (saofcf(ix,it)      > i2_missval .AND. saofcf(ix,it) < 0_i2) .OR. &
            (saocol(ix,it)      > r8_missval .AND. col3sig < 0.0_r8    ) ) THEN

          saomqf(ix,it) = main_qa_bad

          NrofGoodInputSamples = NrofGoodInputSamples + 1
          NrofBadOutputSamples = NrofBadOutputSamples + 1

          IF ( saocol(ix,it) > r8_missval .AND. col3sig < 0.0_r8 ) &
            NrofOutofBoundsSamples        = NrofOutofBoundsSamples        + 1
          IF ( saofcf(ix,it) == opt_convergence_failed .or. saofcf(ix,it) == opt_convergence_maxiter_exceeded ) &
            NrofFailedConvergenceSamples  = NrofFailedConvergenceSamples  + 1
          IF ( saofcf(ix,it) == opt_convergence_maxiter_exceeded)                      &
            NrofExceededIterationsSamples = NrofExceededIterationsSamples + 1

          CYCLE
        END IF

        ! ----------------------------------------------------------
        ! The Ugly: Whatever is left (outside plain missing columns)
        ! ----------------------------------------------------------
        IF ( saocol(ix,it) > r8_missval ) THEN

          IF ( (saofcf(ix,it) == opt_convergence_suspect) .OR. &
              (col2sig <  0.0_r8  .AND. col3sig >= 0.0_r8                      ) .OR. &
              (ABS(saocol(ix,it)) > max_good_col                               ) ) THEN

            saomqf(ix,it) = main_qa_suspect

            NrofGoodInputSamples     = NrofGoodInputSamples     + 1
            NrofSuspectOutputSamples = NrofSuspectOutputSamples + 1

            CYCLE
          END IF

        ELSE

          ! ----------------------------------------------------------
          ! The Missing: Not processed because of either missing input
          !              or restrictions on lat, lon, sza, etc.
          ! ----------------------------------------------------------
          NrofMissingSamples = NrofMissingSamples + 1

        END IF

      END DO
    END DO

    ! --------------------------------------------
    ! Now we can compute averages and percentages,
    ! and write out the final statistics
    ! --------------------------------------------

    IF ( nfitcol >= 1.0_r8 ) THEN
      fitcol_avg  = fitcol_avg  / nfitcol
      rms_avg     = rms_avg     / nfitcol
      dfitcol_avg = dfitcol_avg / nfitcol
    END IF

    PercentGoodOutputSamples      = 100_r4    * &
      REAL(NrofGoodOutputSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples,  KIND=r4) )

    PercentBadOutputSamples       = 100_r4        * &
      REAL(NrofBadOutputSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )

    PercentSuspectOutputSamples   =  100.0_r4         * &
      REAL(NrofSuspectOutputSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )

    PercentOutofBoundsSamples     =  100.0_r4         * &
      REAL(NrofOutofBoundsSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )

    AbsolutePercentMissingSamples = 100_r4 * &
      REAL(NrofMissingSamples, KIND=r4) / &
      MAX ( 1.0_4, REAL(NrofInputSamples, KIND=r4) )

    QAPercentMissingData     = NINT ( AbsolutePercentMissingSamples, KIND=i4 )
    QAPercentOutofBoundsData = NINT ( PercentOutofBoundsSamples,     KIND=i4 )

    ! ------------------------------------------------------------------------
    ! With the above information we can easily determine the Automatic QA Flag
    ! ------------------------------------------------------------------------
    CALL set_automatic_quality_flag ( PercentGoodOutputSamples )

    IF ( verb_thresh_lev >= vb_lev_screen ) THEN
      WRITE (*, '(A, 3(1PE15.5))')          'Col-DCol-RMS: ', fitcol_avg, dfitcol_avg, rms_avg
      WRITE (*, '(A, I7,A,I7,A,F7.1,A)')  'Statistics:   ', &
        MAX(NrofGoodOutputSamples,0), ' of ', MAX(NrofGoodInputSamples,0), ' converged - ', &
        MAX(PercentGoodOutputSamples, 0.0), '%'
      WRITE (*, '(A, F7.1)') 'Nfitcol =', nfitcol
    END IF

    CALL he5_write_fitting_statistics ( &
      pge_idx, max_good_col, nxtrack, ntimes, saomqf(1:nxtrack,0:ntimes-1), &
      fitcol_avg, dfitcol_avg, rms_avg, locerrstat                            )
    errstat = MAX ( locerrstat, errstat )

    RETURN

  END SUBROUTINE compute_fitting_statistics

  SUBROUTINE set_input_pointer_and_versions ( pge_idx )

    USE OMSAO_precision_module,    ONLY: i4
    USE OMSAO_prefitcol_module,    ONLY : yn_o3_prefit, yn_bro_prefit,&
      yn_lqh2o_prefit
    USE OMSAO_indices_module,      ONLY: &
      pge_oclo_idx, pge_bro_idx, pge_hcho_idx, pge_o3_idx,    &
      pge_gly_idx, l1b_radiance_lun, l1b_radianceref_lun, l1b_irradiance_lun, &
      o3_prefit_lun, bro_prefit_lun, lqh2o_prefit_lun,                        &
      voc_amf_luns, voc_omicld_idx, pge_h2o_idx
    USE OMSAO_he5_module,          ONLY: n_lun_inp, lun_input, input_versions
    USE OMSAO_variables_module,    ONLY: l1b_rad_filename, &
      l1b_radref_filename
    use ctrlvars, only: yn_radiance_reference, yn_solar_comp

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx

    ! ---------------
    ! Local variables
    ! ---------------
    LOGICAL :: do_radref

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=31), PARAMETER :: modulename = 'set_input_pointer_and_versions'

    ! ----------------------------------------
    ! Initialize variables returned via MODULE
    ! ----------------------------------------
    lun_input      = -1
    n_lun_inp      =  0

    ! --------------------------------------------------------------------
    ! (Almost) Common InputVersions for all PGEs: OMBRUG/OMBRVG and OMBIRR
    ! --------------------------------------------------------------------
    ! The total number of PGE input files depends on
    !
    ! (a) Whether a solar composite spectrum is being used
    ! (b) Whether a radiance reference from a granule other than the one
    !     being processed is being used.
    ! --------------------------------------------------------------------
    ! (0) Processed granule
    ! ---------------------
    n_lun_inp = 1
    lun_input(n_lun_inp) = l1b_radiance_lun
    ! --------------------------------
    ! (a) Earthshine reference granule
    ! --------------------------------
    IF ( yn_radiance_reference .AND. &
        ( TRIM(ADJUSTL(l1b_rad_filename)) /= TRIM(ADJUSTL(l1b_radref_filename))) ) THEN
      n_lun_inp = n_lun_inp + 1
      lun_input(n_lun_inp) = l1b_radianceref_lun
      do_radref = .TRUE.
    ELSE
      do_radref = .FALSE.
    END IF
    ! --------------------
    ! (b) Solar Irradiance
    ! --------------------
    IF ( .NOT. yn_solar_comp ) THEN
      n_lun_inp = n_lun_inp + 1
      lun_input(n_lun_inp) = l1b_irradiance_lun
    END IF

    ! -------------------------------------------------------------------
    ! Depending on the PGE we are running, we need to add some input LUNs
    ! -------------------------------------------------------------------
    SELECT CASE ( pge_idx )
    CASE (pge_oclo_idx)
      ! -----------------------
      ! Nothing to be done here
      ! -----------------------
    CASE (pge_bro_idx)
      ! -----------------------
      ! Nothing to be done here
      ! -----------------------
    CASE (pge_hcho_idx)
      ! -----------------
      ! Add the Cloud LUN
      ! -----------------
      n_lun_inp            = n_lun_inp + 1
      lun_input(n_lun_inp) = voc_amf_luns(voc_omicld_idx)
      ! ----------------------------------------------------------
      ! Add possibly pre-fitted OMSAO3 and OMBRO
      ! ----------------------------------------------------------
      IF ( yn_o3_prefit(1) ) THEN
        n_lun_inp            = n_lun_inp + 1
        lun_input(n_lun_inp) = o3_prefit_lun
      END IF
      IF ( yn_bro_prefit(1) ) THEN
        n_lun_inp            = n_lun_inp + 1
        lun_input(n_lun_inp) = bro_prefit_lun
      END IF
    CASE (pge_gly_idx)
      ! -----------------
      ! Add the Cloud LUN
      ! -----------------
      n_lun_inp            = n_lun_inp + 1
      lun_input(n_lun_inp) = voc_amf_luns(voc_omicld_idx)

      ! --------------------
      ! Pre-fitted lqH2O CCM
      ! --------------------
      IF ( yn_lqh2o_prefit(1) ) THEN
        n_lun_inp            = n_lun_inp + 1
        lun_input(n_lun_inp) = lqh2o_prefit_lun
      END IF
    CASE (pge_h2o_idx)
      ! -----------------
      ! Add the Cloud LUN
      ! -----------------
      n_lun_inp            = n_lun_inp + 1
      lun_input(n_lun_inp) = voc_amf_luns(voc_omicld_idx)

    CASE (pge_o3_idx)
      ! -----------------------
      ! Nothing to be done here
      ! -----------------------
    END SELECT

    ! ------------------------------------------------------------
    ! Composing the InputVersion string is more difficult, because
    ! we have to compose the pieces of information from various
    ! MetaData strings.
    ! ------------------------------------------------------------
    CALL get_input_versions ( pge_idx, do_radref, input_versions )
    input_versions = TRIM(ADJUSTL(input_versions))

    RETURN
  END SUBROUTINE set_input_pointer_and_versions

  SUBROUTINE get_input_versions (pge_idx, do_radref, input_versions )

    USE OMSAO_precision_module,  ONLY: i4
    USE OMSAO_indices_module,    ONLY: pge_hcho_idx, pge_gly_idx, pge_h2o_idx
    USE OMSAO_prefitcol_module,  ONLY: yn_o3_prefit, yn_bro_prefit, yn_lqh2o_prefit
    USE metadata_tools, ONLY: PGSd_MET_NAME_L, n_mdata_omhcho, &
      n_mdata_str, n_mdata_voc, mdata_string_fields, n_mdata_omchocho, &
      mdata_string_values, mdata_voc_fields, mdata_omhcho_fields, &
      mdata_omhcho_values, mdata_voc_values, &
      mdata_omchocho_fields, mdata_omchocho_values
    use ctrlvars, only: yn_solar_comp

    IMPLICIT NONE

    ! --------------
    ! Input variable
    ! --------------
    INTEGER,           INTENT (IN)  :: pge_idx
    LOGICAL,           INTENT (IN)  :: do_radref

    ! ---------------
    ! Output variable
    ! ---------------
    CHARACTER (LEN=*), INTENT (OUT) :: input_versions

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)               :: i
    CHARACTER (LEN=PGSd_MET_NAME_L) :: &
      rad_name, rad_version, &  ! Radiance granule
      rrf_name, rrf_version, &  ! Radiance Reference granule
      irr_name, irr_version, &  ! Irradiance granule
      cld_name, cld_version, &  ! Cloud ESDT
      bro_name, bro_version, &  ! BrO Prefit ESDT
      ooo_name, ooo_version, &  ! Ozone Prefit ESDT
      lqh2o_name, lqh2o_version   ! lqH2O Prefit ESDT
    ! --------------------------
    ! Initialize output variable
    ! --------------------------
    input_versions = ''

    ! -------------------------------------------------------
    ! Collect information on Radiance and Irradiance versions
    ! -------------------------------------------------------
    DO i = 1, n_mdata_str
      IF ( TRIM(ADJUSTL(mdata_string_fields(1,i))) == 'PGEVERSION' ) THEN
        IF ( TRIM(ADJUSTL(mdata_string_fields(3,i))) == 'l1r' ) &
          rad_version = TRIM(ADJUSTL(mdata_string_values(i)))
        IF ( TRIM(ADJUSTL(mdata_string_fields(3,i))) == 'l1R' ) &
          rrf_version = TRIM(ADJUSTL(mdata_string_values(i)))
        IF ( TRIM(ADJUSTL(mdata_string_fields(3,i))) == 'l1i' ) &
          irr_version = TRIM(ADJUSTL(mdata_string_values(i)))
      END IF
      IF ( TRIM(ADJUSTL(mdata_string_fields(1,i))) == 'ShortName' ) THEN
        IF ( TRIM(ADJUSTL(mdata_string_fields(3,i))) == 'l1r' ) &
          rad_name = TRIM(ADJUSTL(mdata_string_values(i)))
        IF ( TRIM(ADJUSTL(mdata_string_fields(3,i))) == 'l1R' ) &
          rrf_name = TRIM(ADJUSTL(mdata_string_values(i)))
        IF ( TRIM(ADJUSTL(mdata_string_fields(3,i))) == 'l1i' ) &
          irr_name = TRIM(ADJUSTL(mdata_string_values(i)))
      END IF
    END DO

    ! ------------------------------------------------------------------
    ! As per request, the full version string rather than only the first
    ! decimal digit will be appended to the InputVersions string. This
    ! makes obsolete determining the position of the first "." in the
    ! version strings. See also in block below for L2 ESTDs.
    ! ------------------------------------------------------------------
    ! Start with the radiance granule currently being processed
    ! ---------------------------------------------------------
    input_versions = TRIM(ADJUSTL(rad_name))//':'//TRIM(ADJUSTL(rad_version))
    ! ----------------------------------
    ! Add the radiance reference granule
    ! ----------------------------------
    IF ( do_radref ) &
      input_versions = TRIM(ADJUSTL(input_versions)) // ' ' // &
      TRIM(ADJUSTL(rrf_name))//':'//TRIM(ADJUSTL(rrf_version))

    ! ---------------------
    ! Add the solar granule
    ! ---------------------
    IF ( .NOT. yn_solar_comp ) &
      input_versions = TRIM(ADJUSTL(input_versions))  // ' ' // &
      TRIM(ADJUSTL(irr_name))//':'//TRIM(ADJUSTL(irr_version))

    ! ---------------------------------------------------
    ! For OMHCHO and OMCHOCHO we have to add some extras.
    ! ---------------------------------------------------
    SELECT CASE ( pge_idx )
    CASE ( pge_hcho_idx )
      ! ---------------------------------------------------
      ! In all cases we should have auxilliary cloud inputs
      ! ---------------------------------------------------
      DO i = 1, n_mdata_voc
        IF ( TRIM(ADJUSTL(mdata_voc_fields(3,i))) == 'CLD' ) THEN
          IF ( TRIM(ADJUSTL(mdata_voc_fields(1,i))) == 'PGEVERSION' ) &
            cld_version = TRIM(ADJUSTL(mdata_voc_values(i)))
          IF ( TRIM(ADJUSTL(mdata_voc_fields(1,i))) == 'ShortName'  ) &
            cld_name    = TRIM(ADJUSTL(mdata_voc_values(i)))
        END IF
      END DO

      ! ----------------------------
      ! Full Version string required
      ! ----------------------------
      input_versions = TRIM(ADJUSTL(input_versions)) // ' ' // &
        TRIM(ADJUSTL(cld_name))//':'//TRIM(ADJUSTL(cld_version))

      ! ----------------------------------------------------------
      ! Add possibly pre-fitted OMSAO3 and OMBRO, and OMCLDO2 ESDT
      ! ----------------------------------------------------------
      IF ( yn_o3_prefit(1) ) THEN
        DO i = 1, n_mdata_omhcho
          IF ( TRIM(ADJUSTL(mdata_omhcho_fields(3,i))) ==  'OOO' ) THEN
            IF ( TRIM(ADJUSTL(mdata_omhcho_fields(1,i))) == 'PGEVERSION' ) &
              ooo_version = TRIM(ADJUSTL(mdata_omhcho_values(i)))
            IF ( TRIM(ADJUSTL(mdata_omhcho_fields(1,i))) == 'ShortName'  ) &
              ooo_name    = TRIM(ADJUSTL(mdata_omhcho_values(i)))
          END IF
        END DO
        ! ----------------------------
        ! Full Version string required
        ! ----------------------------
        input_versions = TRIM(ADJUSTL(input_versions))            // ' ' // &
          TRIM(ADJUSTL(ooo_name))//':'//TRIM(ADJUSTL(ooo_version))
      END IF
      IF ( yn_bro_prefit(1) ) THEN
        DO i = 1, n_mdata_omhcho
          IF ( TRIM(ADJUSTL(mdata_omhcho_fields(3,i))) == 'BRO' ) THEN
            IF ( TRIM(ADJUSTL(mdata_omhcho_fields(1,i))) == 'PGEVERSION' ) &
              bro_version = TRIM(ADJUSTL(mdata_omhcho_values(i)))
            IF ( TRIM(ADJUSTL(mdata_omhcho_fields(1,i))) == 'ShortName'  ) &
              bro_name    = TRIM(ADJUSTL(mdata_omhcho_values(i)))
          END IF
        END DO
        ! ----------------------------
        ! Full Version string required
        ! ----------------------------
        input_versions = TRIM(ADJUSTL(input_versions))            // ' ' // &
          TRIM(ADJUSTL(bro_name))//':'//TRIM(ADJUSTL(bro_version))
      END IF

    CASE ( pge_gly_idx)
      ! ---------------------------------------------------
      ! In all cases we should have auxilliary cloud inputs
      ! ---------------------------------------------------
      DO i = 1, n_mdata_voc
        IF ( TRIM(ADJUSTL(mdata_voc_fields(3,i))) == 'CLD' ) THEN
          IF ( TRIM(ADJUSTL(mdata_voc_fields(1,i))) == 'PGEVERSION' ) &
            cld_version = TRIM(ADJUSTL(mdata_voc_values(i)))
          IF ( TRIM(ADJUSTL(mdata_voc_fields(1,i))) == 'ShortName'  ) &
            cld_name    = TRIM(ADJUSTL(mdata_voc_values(i)))
        END IF
      END DO

      ! ----------------------------
      ! Full Version string required
      ! ----------------------------
      input_versions = TRIM(ADJUSTL(input_versions)) // ' ' // &
        TRIM(ADJUSTL(cld_name))//':'//TRIM(ADJUSTL(cld_version))

      ! lqH2O prefit
      IF ( yn_lqh2o_prefit(1) ) THEN
        DO i = 1, n_mdata_omchocho
          IF ( TRIM(ADJUSTL(mdata_omchocho_fields(3,i))) == 'LQH2O' ) THEN
            IF ( TRIM(ADJUSTL(mdata_omchocho_fields(1,i))) == 'PGEVERSION' ) &
              lqh2o_version = TRIM(ADJUSTL(mdata_omchocho_values(i)))
            IF ( TRIM(ADJUSTL(mdata_omchocho_fields(1,i))) == 'ShortName'  ) &
              lqh2o_name    = TRIM(ADJUSTL(mdata_omchocho_values(i)))
          END IF
        END DO
        ! ----------------------------
        ! Full Version string required
        ! ----------------------------
        input_versions = TRIM(ADJUSTL(input_versions))            // ' ' // &
          TRIM(ADJUSTL(lqh2o_name))//':'//TRIM(ADJUSTL(lqh2o_version))
      END IF

    CASE ( pge_h2o_idx)
      ! ---------------------------------------------------
      ! In all cases we should have auxilliary cloud inputs
      ! ---------------------------------------------------
      DO i = 1, n_mdata_voc
        IF ( TRIM(ADJUSTL(mdata_voc_fields(3,i))) == 'CLD' ) THEN
          IF ( TRIM(ADJUSTL(mdata_voc_fields(1,i))) == 'PGEVERSION' ) &
            cld_version = TRIM(ADJUSTL(mdata_voc_values(i)))
          IF ( TRIM(ADJUSTL(mdata_voc_fields(1,i))) == 'ShortName'  ) &
            cld_name    = TRIM(ADJUSTL(mdata_voc_values(i)))
        END IF
      END DO

    END SELECT

    RETURN
  END SUBROUTINE get_input_versions

!UNUSED!   SUBROUTINE omi_radiance_wvl_smoothing ( nxt, nwl, omi_radiance_wavl )
!UNUSED! 
!UNUSED!     USE OMSAO_precision_module, ONLY: i4, r4, r8
!UNUSED!     USE SLATEC_davint, ONLY: dpolft
!UNUSED! 
!UNUSED!     ! ---------------
!UNUSED!     ! Input Variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4), INTENT (IN) :: nxt, nwl
!UNUSED! 
!UNUSED!     ! ------------------
!UNUSED!     ! Modified Variables
!UNUSED!     ! ------------------
!UNUSED!     REAL (KIND=r4), DIMENSION(nwl,nxt), INTENT (INOUT) :: omi_radiance_wavl
!UNUSED! 
!UNUSED!     ! ------------------------------
!UNUSED!     ! Local Variables and Parameters
!UNUSED!     ! ------------------------------
!UNUSED!     INTEGER (KIND=i4), PARAMETER       :: max_deg = 2
!UNUSED!     INTEGER (KIND=i4)                  :: i, ierr, ndeg
!UNUSED!     REAL    (KIND=r8)                  :: eps
!UNUSED!     REAL    (KIND=r8), DIMENSION (nxt) :: x, y, w, yf
!UNUSED!     REAL    (KIND=r8), DIMENSION (3*(nxt+max_deg+1)) :: a
!UNUSED! 
!UNUSED!     x(1:nxt) = (/ (REAL(i,KIND=KIND(r8)), i = 1, nxt) /)
!UNUSED!     w(1:nxt) = -1.0_r8
!UNUSED!     eps      = 0.0_r8
!UNUSED! 
!UNUSED!     DO i = 1, nwl
!UNUSED!       y(1:nxt) = omi_radiance_wavl(i,1:nxt)
!UNUSED!       CALL dpolft (nxt, x, y, w, max_deg, ndeg, eps, yf, ierr, a)
!UNUSED!       omi_radiance_wavl(i,1:nxt) = REAL(yf(1:nxt),KIND=r4)
!UNUSED!     END DO
!UNUSED! 
!UNUSED!     RETURN
!UNUSED!   END SUBROUTINE omi_radiance_wvl_smoothing

!UNUSED!   SUBROUTINE compact_fitting_spectrum ( n_fit_wvl, curr_fit_spec )
!UNUSED! 
!UNUSED!     USE OMSAO_precision_module,  ONLY: i4, r8
!UNUSED!     USE OMSAO_indices_module,    ONLY: wvl_idx, sig_idx, ccd_idx
!UNUSED!     USE OMSAO_parameters_module, ONLY: normweight
!UNUSED! 
!UNUSED!     IMPLICIT NONE
!UNUSED! 
!UNUSED!     ! ------------------
!UNUSED!     ! Modified variables
!UNUSED!     ! ------------------
!UNUSED!     INTEGER (KIND=i4),                                        INTENT (INOUT) :: n_fit_wvl
!UNUSED!     REAL    (KIND=r8), DIMENSION (wvl_idx:ccd_idx,n_fit_wvl), INTENT (INOUT) :: curr_fit_spec
!UNUSED! 
!UNUSED!     ! ---------------
!UNUSED!     ! Local variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4)                                        :: i, n_loc_wvl
!UNUSED!     REAL    (KIND=r8), DIMENSION (wvl_idx:ccd_idx,n_fit_wvl) :: loc_fit_spec
!UNUSED! 
!UNUSED!     n_loc_wvl = 0_i4  ;  loc_fit_spec = 0.0_r8
!UNUSED!     DO i = 1, n_fit_wvl
!UNUSED!       IF ( curr_fit_spec(sig_idx,i) == normweight ) THEN
!UNUSED!         n_loc_wvl = n_loc_wvl + 1
!UNUSED!         loc_fit_spec(wvl_idx:ccd_idx,n_loc_wvl) = curr_fit_spec(wvl_idx:ccd_idx,i)
!UNUSED!       END IF
!UNUSED!     END DO
!UNUSED! 
!UNUSED!     ! ---------------------------
!UNUSED!     ! Assign the output variables
!UNUSED!     ! ---------------------------
!UNUSED!     n_fit_wvl     = n_loc_wvl
!UNUSED!     curr_fit_spec = loc_fit_spec
!UNUSED! 
!UNUSED!     RETURN
!UNUSED!   END SUBROUTINE compact_fitting_spectrum

  SUBROUTINE check_wavelength_overlap ( &
      n_fitvar_rad, n_sol_wvl, irradiance_wvl, n_rad_wvl, radiance_wvl, &
      do_cycle_this_pix )

    USE OMSAO_precision_module,  ONLY: i4, r8

    IMPLICIT NONE

    ! Input variables
    INTEGER (KIND=i4),                        INTENT (IN) :: n_sol_wvl, n_rad_wvl, n_fitvar_rad
    REAL    (KIND=r8), DIMENSION (n_sol_wvl), INTENT (IN) :: irradiance_wvl
    REAL    (KIND=r8), DIMENSION (n_rad_wvl), INTENT (IN) :: radiance_wvl

    ! Output variable
    LOGICAL, INTENT (OUT) :: do_cycle_this_pix

    ! Local variables
    INTEGER (KIND=i4) :: j, n_overlap1, n_overlap2

    do_cycle_this_pix = .FALSE.

    IF ( radiance_wvl(1)         >= irradiance_wvl(n_sol_wvl) .OR. &
        radiance_wvl(n_rad_wvl) <= irradiance_wvl(1)                 ) THEN
      do_cycle_this_pix = .TRUE.
      RETURN
    END IF

    n_overlap1 = 0
    DO j = 1, n_rad_wvl
      IF ( radiance_wvl(j) >= irradiance_wvl(1)         .AND. &
          radiance_wvl(j) <= irradiance_wvl(n_sol_wvl)         ) n_overlap1 = n_overlap1 + 1
    END DO
    IF ( n_overlap1 < n_fitvar_rad ) THEN
      do_cycle_this_pix = .TRUE.
      RETURN
    END IF

    n_overlap2 = 0
    DO j = 1, n_sol_wvl
      IF ( irradiance_wvl(j) >= radiance_wvl(1)         .AND. &
          irradiance_wvl(j) <= radiance_wvl(n_rad_wvl)         ) n_overlap2 = n_overlap2 + 1
    END DO
    IF ( n_overlap2 < n_fitvar_rad ) THEN
      do_cycle_this_pix = .TRUE.
      RETURN
    END IF

    RETURN
  END SUBROUTINE check_wavelength_overlap

  SUBROUTINE omi_set_xtrpix_range (                 &
      nTimes, nXtrack, pixnum_limits, omi_binfac,  &
      omi_xtrpix_range, first_wc_pix, last_wc_pix, &
      errstat )

    USE OMSAO_precision_module,  ONLY: i1, i4
    USE OMSAO_omidata_module,    ONLY: szoom_mode, global_mode, gzoom_spix, gzoom_epix
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! =======================================================================
    !
    ! Purpose of this routine:
    !
    !    Set first and last cross track pixel position to process
    !
    ! The problem:
    !
    !    The usual approach to fit all cross-track pixel, or a subset
    !    of them that has been specified in the fitting control file,
    !    fails for granules that contain spatial zoom data, since those
    !    have the center 30 cross-track positions 16-45 stored with the
    !    rest missing.
    !
    ! The approach:
    !
    !    For each swath line we define first and last cross track pixel
    !    based on whether it is global or spatial zoom. Data are stored
    !    as usual (1:60, respectively 16:45) but the "first" and "last"
    !    cross-track pixel to process is adjusted in cases of spatial
    !    zoom mode to avoid needless checking of boundaries.
    !
    !    "Mixed Mode" cases are conceivable but pathological. They are
    !    not considered here.
    !
    !    For the radiance fits, an array with first/last pixel and shift is
    !    set up for each swath line.
    !
    ! =======================================================================

    ! ===================================================================
    !
    ! Explanation of subroutine arguments
    !
    !  nTimes .............. number of swath lines
    !  nXtrack ............. number of cross track ("XT") positions
    !  pixnum_limits ....... constraints on XT pixels
    !  omy_binfac .......... Binning factor for each scan line, either
    !                           szoom_mode or global_mode.
    !
    !  first_wc_pix ........ first XT position for wavelengh calibration
    !  last_wc_pix ......... last XT position for wavelengh calibration
    !  omi_xtrpix_range .... (nTimes,2) array of first and last XT for
    !                        radiance fitting
    !
    ! ====================================================================

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nTimes, nXtrack
    INTEGER (KIND=i4), DIMENSION(2),          INTENT (IN) :: pixnum_limits
    INTEGER (KIND=i1), DIMENSION(0:nTimes-1), INTENT (IN) :: omi_binfac

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                          INTENT (OUT) :: first_wc_pix, last_wc_pix
    INTEGER (KIND=i4), DIMENSION(0:nTimes-1,2), INTENT (OUT) :: omi_xtrpix_range

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)          :: first_pix, last_pix, i, locerrstat

    ! -----------------------
    ! Name of this subroutine
    ! -----------------------
    !CHARACTER (LEN=20), PARAMETER :: modulename = 'omi_set_xtrpix_range'

    ! ---------------------------
    ! Initialize return variables
    ! ---------------------------
    omi_xtrpix_range(0:nTimes-1,1:2) = -1
    first_wc_pix                     = -1
    last_wc_pix                      = -1

    locerrstat = pge_errstat_ok

    ! ------------------------------------------
    ! Find the range of XT pixels to process
    ! ------------------------------------------
    first_pix = MAX (                                  1, pixnum_limits(1) )
    last_pix  = MAX ( MIN ( nXtrack,  pixnum_limits(2) ), first_pix        )

    ! --------------------------------------------------------------------
    ! We go through the pixels one-by-one. Lots of redundancy, and no
    ! guarantee that we are actually doing the right thing. Proceed with
    ! fingers crossed.
    ! --------------------------------------------------------------------
    DO i = 0, nTimes -1
      SELECT CASE ( omi_binfac(i) )
      CASE ( global_mode )
        first_wc_pix = first_pix
        last_wc_pix  = last_pix
      CASE ( szoom_mode )
        first_wc_pix = MAX ( first_pix, gzoom_spix )
        last_wc_pix  = MIN ( last_pix,  gzoom_epix )
      CASE DEFAULT
        ! ---------------------------------------------------
        ! Let's hope we never reach here. This spells TROUBLE.
        ! Setting the end pixel to less than the start pixel
        ! is a feeble attempt to avoid the worst of loops.
        ! ---------------------------------------------------
        first_wc_pix = -1 ; last_wc_pix = -2
      END SELECT
      omi_xtrpix_range(i,1) = first_wc_pix
      omi_xtrpix_range(i,2) = last_wc_pix
    END DO

    errstat = MAX ( errstat, locerrstat )
    RETURN
  END SUBROUTINE omi_set_xtrpix_range

  SUBROUTINE convert_tai_to_utc ( nUTCdim, time_tai, time_utc )

    USE OMSAO_precision_module, ONLY: i2, i4, r8

    IMPLICIT NONE

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER, EXTERNAL :: PGS_TD_TAItoUTC

    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nUTCdim
    REAL    (KIND=r8), INTENT (IN) :: time_tai

    ! ---------------
    ! Output Variable
    ! ---------------
    INTEGER (KIND=i2), DIMENSION (nUTCdim), INTENT (OUT) :: time_utc

    ! --------------
    ! Local Variable
    ! --------------
    INTEGER   (KIND=i4 ) :: locerrstat
    CHARACTER (LEN=27)   :: utc_string

    ! -----------------------
    ! Convert TAI to UTC time
    ! -----------------------
    utc_string = ""
    locerrstat = PGS_TD_TAItoUTC ( time_tai, utc_string )

    ! ---------------------------------------------------------------
    ! Now we convert the UTC string to INTEGER values. We had rather
    ! write this to file as a CHAR but at least HE5 v1.6.4 has issues
    ! with that. Hence this kludge.
    ! ---------------------------------------------------------------
    ! The format of the UTC string is: YYYY-MM-DDThh:mm:ss.ddddddZ
    ! where
    !         YYYY = year             (4 characters)
    !         MM   = month            (2)
    !         DD   = day              (2)
    !         T    = "T"              (separator)
    !         hh   = hour             (2)
    !         mm   = minutes          (2)
    !         ss   = seconds          (2)
    !         d    = decimal fraction (6)
    !         Z    = "Z" (terminator)
    ! ------------------------------------------------------------------
    ! The conversion below is sort-of Low-Brow, but there doesn't seem
    ! to be a simple routine that returns "interger value of a numerical
    ! character". Hence the recourse to old F77 practices.
    ! ------------------------------------------------------------------
    READ ( utc_string( 1: 4), * ) time_utc(1)   ! YYYY
    READ ( utc_string( 6: 7), * ) time_utc(2)   ! MM
    READ ( utc_string( 9:10), * ) time_utc(3)   ! DD
    READ ( utc_string(12:13), * ) time_utc(4)   ! hh
    READ ( utc_string(15:16), * ) time_utc(5)   ! mm
    READ ( utc_string(18:19), * ) time_utc(6)   ! ss

    RETURN
  END SUBROUTINE convert_tai_to_utc

  SUBROUTINE find_swathline_range ( &
      l1bfile, l1bswath, nt, nx, l1blats, latrange, in_range, errstat )

    USE OMSAO_precision_module
    USE OMSAO_variables_module,  ONLY: pixnum_lim
    USE OMSAO_errstat_module
    USE omi_read_l1b_data, ONLY: omi_read_binning_factor

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=20), PARAMETER :: modulename = 'find_swathline_range'

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*),                            INTENT (IN) :: l1bfile, l1bswath
    INTEGER   (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL      (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: l1blats
    REAL      (KIND=r4), DIMENSION (2),           INTENT (IN) :: latrange

    ! -----------------------------
    ! Output and Modified variables
    ! -----------------------------
    LOGICAL, DIMENSION (0:nt-1), INTENT (INOUT) :: in_range
    INTEGER (KIND=i4),           INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                          :: estat, fpix, lpix, midnum, locerrstat
    REAL    (KIND=r4)                          :: midlat
    INTEGER (KIND=i4), DIMENSION (0:nt-1, 1:2) :: xtrange
    INTEGER (KIND=i1), DIMENSION (0:nt-1)      :: binfac
    LOGICAL,           DIMENSION (0:nt-1)      :: ynzoom

    locerrstat = pge_errstat_ok
    estat      = pge_errstat_ok

    ! ----------------------------------------------------------------
    ! Read preparatory arrays for determining the range of swath lines
    ! that fall within the desired latitude interval.
    ! ----------------------------------------------------------------
    CALL omi_read_binning_factor ( &
      l1bfile, l1bswath, nt, binfac(0:nt-1), ynzoom(0:nt-1), estat )

    CALL omi_set_xtrpix_range ( &
      nt, nx, pixnum_lim(3:4), binfac(0:nt-1), &
      xtrange(0:nt-1,1:2), fpix, lpix, estat    )

    ! ----------------------------------------------------------------------
    ! Determine the range of swath line numbers that go into the radiance
    ! reference spectrum. This is based either on a finite latitude interval
    ! or on a single latitude.
    ! ----------------------------------------------------------------------
    IF ( latrange(1) /= latrange(2) ) THEN
      midlat = SUM(latrange(1:2)) / 2.0_r4
    ELSE
      midlat = latrange(1)
    END IF

    CALL find_swathrange_by_latitude (                      &
      nt, nx, latrange(1), latrange(2),                  &
      l1blats(1:nx,0:nt-1), xtrange(0:nt-1,1:2), midlat, &
      midnum, in_range(0:nt-1)                        )

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE find_swathline_range

  SUBROUTINE find_swathline_by_latitude ( &
      nxrr, sline, eline, latr4, lat, xtrange, lnum, was_found )

    USE OMSAO_precision_module, ONLY: i4, r4
    USE OMSAO_parameters_module, ONLY: r4_missval
    USE OMSAO_errstat_module
    !USE L1B_Reader_class

    IMPLICIT NONE

    ! --------------------------------------------------------------------------
    ! This subroutine returns a swath line number from an OMI swath
    ! based on a given latitude. Subroutine arguments:
    !
    ! nxrr ............. Number of cross-track entries in OMI swath
    ! sline ............ Lowest swath line number
    ! eline ............ Highest swath line number
    ! latr4 ............ Latitude array for whole swath
    ! lat .............. Latitude to locate
    ! xtrange .......... Number of valid cross-track positions in swath
    !                    (possibly smaller than nx)
    ! was_found ......... TRUE if line number has been found, FALSE otherwise
    ! --------------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                 INTENT (IN) :: nxrr, sline, eline
    REAL    (KIND=r4),                                 INTENT (IN) :: lat
    INTEGER (KIND=i4), DIMENSION (sline:eline,2),      INTENT (IN) :: xtrange
    REAL    (KIND=r4), DIMENSION (1:nxrr,sline:eline), INTENT (IN) :: latr4

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: lnum
    LOGICAL,           INTENT (OUT) :: was_found

    ! ---------------
    ! Local variables
    ! ---------------
    REAL    (KIND=r4)                   :: diff, mindiff
    REAL    (KIND=r4), DIMENSION (nxrr) :: cntr4, latdiff, loclat
    INTEGER (KIND=i4)                   :: &
      j1, j2, icnt, iline, fpix, lpix, locerr

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    lnum = -1
    was_found = .FALSE.

    ! --------------------------------------------------------------------------
    ! First, start a bisection of the [0, NLINES-1] interval to find the closest
    ! match in latitude to the mipoint of the latitude regime to average.
    ! --------------------------------------------------------------------------
    j1 = sline ; j2 = eline-1  ;  icnt = 0
    FindLine: DO WHILE ( .NOT. was_found )
      icnt  = icnt + 1
      iline = (j1 + j2) / 2

      ! -----------------------------------------------------------------------
      ! Get first and last pixel.
      ! -----------------------------------------------------------------------
      fpix = xtrange(iline,1)
      lpix = xtrange(iline,2)

      IF ( iline < sline .OR. iline > eline ) THEN
        locerr = pge_errstat_error
        EXIT FindLine
      END IF

      loclat(fpix:lpix) = latr4(fpix:lpix,iline)
      cntr4 (fpix:lpix) = 1.0_r4

      WHERE ( ABS(loclat(fpix:lpix)) > 90.0_r4 )
        loclat(fpix:lpix) = r4_missval
        cntr4 (fpix:lpix) = 0.0_r4
      END WHERE

      IF ( MAXVAL(loclat(fpix:lpix)) > r4_missval ) THEN
        latdiff(fpix:lpix) = ( loclat(fpix:lpix) - lat ) * cntr4(fpix:lpix)
        diff = SUM(latdiff(fpix:lpix))/SUM(cntr4(fpix:lpix))
        IF ( diff < 0.0_r4 ) THEN
          j1 = iline
        ELSE
          j2 = iline
        END IF
        IF ( ABS(j1 - j2) <= 2 ) THEN
          lnum = (j1 + j2) / 2
          was_found     = .TRUE.
          EXIT FindLine
        END IF
      ELSE
        ! -------------------------------------------------------------
        ! Reaching here spells trouble. The only thing we can think of
        ! to do here is to alternately increase and decrease the top
        ! boundary in the hope that we hit upon something. ICNT is
        ! increased by one for each iteration, so "-1**ICNT" has
        ! alternating sign between two iterations. The line below first
        ! subtracts 1, then adds 2, then subtracts 3, asf.
        ! -------------------------------------------------------------
        j2 = j2 + icnt * (-1)**icnt
      END IF
    END DO FindLine

    ! ----------------------------------------------------------------
    ! Now we fine-tune the retrieved scan line number by checking +/-2
    ! scan lines on either side.
    ! ----------------------------------------------------------------
    IF ( was_found ) THEN

      ! -----------------------------------------------------------------------------
      ! MINDIFF will contain the smallest difference found; set to large value first.
      ! -----------------------------------------------------------------------------
      mindiff = REAL((lpix-fpix+1),KIND=r4)*ABS(r4_missval) ; j1 = lnum
      DO iline = lnum-2, lnum+2

        IF ( iline < 0 .OR. iline > eline ) CYCLE

        cntr4(fpix:lpix) = 1.0_r4
        fpix = xtrange(iline,1)
        lpix = xtrange(iline,2)
        loclat(fpix:lpix) = latr4(fpix:lpix,iline)

        WHERE ( ABS(loclat(fpix:lpix)) > 90.0_r4 )
          loclat(fpix:lpix) = r4_missval
          cntr4 (fpix:lpix) = 0.0_r4
        END WHERE

        IF ( MAXVAL(loclat(fpix:lpix)) > r4_missval ) THEN
          latdiff(fpix:lpix) = ( loclat(fpix:lpix) - lat ) * cntr4(fpix:lpix)
          diff = SUM(latdiff(fpix:lpix)*latdiff(fpix:lpix))/SUM(cntr4(fpix:lpix))
          IF ( diff < mindiff ) THEN
            j1      = iline
            mindiff = diff
          END IF
        END IF

      END DO
      lnum = j1
    END IF

    RETURN
  END SUBROUTINE find_swathline_by_latitude

  SUBROUTINE find_swathrange_by_latitude ( &
      nt, nx, latlow, latupp, latr4, xtrange, latmid, latnum, in_range )

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: r4_missval
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------------------------------------------------------------------------
    ! This subroutine returns a swath line number from an OMI swath
    ! based on a given latitude. Subroutine arguments:
    !
    ! nt .................. Number of swath lines
    ! nx .................. Number of cross-track entries in swath
    ! latlow .............. Lowest  latitude to include
    ! latupp .............. Highest latitude to include
    ! latr4 ............... Latitude array for whole swath
    ! xtrange ............. Number of valid cross-track positions in swath
    !                       (possibly smaller than nx)
    ! latmid .............. "Midpoint" latitude closest to average of latlow and latupp
    ! latnum .............. Swath line number of latmid
    ! in_range ......... TRUE if swath line latitude falls between latlow and latupp
    ! ---------------------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4),                          INTENT (IN) :: latlow, latupp, latmid
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: latr4

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),           INTENT (OUT)   :: latnum
    LOGICAL, DIMENSION (0:nt-1), INTENT (INOUT) :: in_range

    ! ---------------
    ! Local variables
    ! ---------------
    REAL    (KIND=r4), PARAMETER      :: dlat = 10.0_r4
    INTEGER (KIND=i4)                 :: iline
    REAL    (KIND=r4), DIMENSION (nx) :: cntr4, latdiff, loclat
    REAL    (KIND=r4)                 :: diff, mindiff
    INTEGER (KIND=i4)                 :: fpix, lpix
    LOGICAL                           :: is_single_lat

    ! -------------------------------------
    ! Initialize output and local variables
    ! -------------------------------------
    latnum = -1 ; mindiff = ABS(r4_missval)

    ! ------------------------------------------------------------
    ! Check whether we are working with a finite latitude interval
    ! or with a single latitude
    ! ------------------------------------------------------------
    is_single_lat = .TRUE.
    IF ( latlow /= latupp ) is_single_lat = .FALSE.

    ! --------------------------------------------------------------------------
    ! Owing to the discontiguous nature of NRT L1b storage, we can't assume that
    ! we are working with monotonously increasing latitues. The only possibility
    ! is to go throught the latitude array from start to finish and flag any
    ! swath line that falls within the desired range. So, out with the speedy
    ! bisection and in with the brute-force plowing through.
    ! --------------------------------------------------------------------------
    GetRange: DO iline = 0, nt-1

      ! -----------------------------------------------------------------------
      ! Get first and last pixel.
      ! -----------------------------------------------------------------------
      fpix = xtrange(iline,1)
      lpix = xtrange(iline,2)

      ! --------------------------------------
      ! Store the pixel slice in a local array
      ! --------------------------------------
      loclat(fpix:lpix) = latr4(fpix:lpix,iline)
      cntr4 (fpix:lpix) = 1.0_r4

      ! ---------------------------------------------------
      ! Replace and out-of-bounds entries by missing values
      ! ---------------------------------------------------
      WHERE ( ABS(loclat(fpix:lpix)) > 90.0_r4 )
        loclat(fpix:lpix) = r4_missval
        cntr4 (fpix:lpix) = 0.0_r4
      END WHERE

      ! -----------------------------------------------------------------
      ! The trick is to handle both an extended latitude interval and a
      ! single latitude value to locate. Since we have 60 cross-track
      ! postions in OMI with a somewhat slanted swath, we can't simply
      ! check for equality with a single latitude value.
      !
      ! Solution: If we are looking for a single latitude only, reject
      !           any swath lines that fall outside a 10deg window of
      !           the target. For the rest, we compute the RMS difference
      !           of the swath latitudes to the target latitude.
      ! -----------------------------------------------------------------

      IF ( is_single_lat ) THEN
        ! -----------------------------------------------------------------------
        ! Check 1: Single latitude. Skip if nothing is within the
        !          [LATMID-DLAT, LAMID+DLAT] interval
        ! -----------------------------------------------------------------------
        IF ( .NOT. (                                        &
          ANY ( loclat(fpix:lpix) >= latmid-dlat ) .AND. &
          ANY ( loclat(fpix:lpix) <= latmid+dlat )      )  ) CYCLE

        ! ----------------------------------------------------------
        ! Update search for the swath line that is closest to LATMID
        ! ----------------------------------------------------------
        latdiff(fpix:lpix) = ( loclat(fpix:lpix) - latmid ) * cntr4(fpix:lpix)
        diff = SQRT( SUM(latdiff(fpix:lpix)*latdiff(fpix:lpix)) ) / SUM(cntr4(fpix:lpix))

        IF ( diff < mindiff ) THEN
          mindiff = diff
          latnum  = iline
        END IF
      ELSE
        ! -----------------------------------------------------------------------
        ! Check 2: Finite latitude interval. Skip if nothing is within the
        !          [LATLOW, LATUPP] interval
        ! -----------------------------------------------------------------------
        IF ( .NOT. (                                    &
          ANY ( loclat(fpix:lpix) >= latlow ) .AND.  &
          ANY ( loclat(fpix:lpix) <= latupp )       )  ) CYCLE

        ! -------------------------
        ! Set in_range to .TRUE.
        ! -------------------------
        in_range(iline) = .TRUE.
      END IF

    END DO GetRange

    ! ------------------------------------------
    ! If working with a single latitude, set the
    ! corresponding line logical to .TRUE.
    ! ------------------------------------------
    IF ( is_single_lat .AND. latnum >= 0 .AND. latnum <= nt-1 ) &
      in_range(latnum) = .TRUE.

    RETURN
  END SUBROUTINE find_swathrange_by_latitude

  ! read_latitude reads all cross track latitudes for ntimes steps from tstart
  SUBROUTINE read_latitude (l1bfile, l1bswath, tstart, ntimes, latr4, errstat )

    USE OMSAO_precision_module, ONLY: i4, r4
    use l1bread, only: l1bread_open_swath, l1bread_close, L1B_Object_Type, &
      l1bread_get2d_r4

    implicit none
    CHARACTER (LEN=*),     INTENT (IN) :: l1bfile, l1bswath
    INTEGER (KIND=i4),     INTENT (IN) :: tstart, ntimes
    integer, intent(inout) :: errstat
    REAL    (KIND=r4), DIMENSION(:,:), INTENT (out) :: latr4

    type (L1B_Object_Type) :: l1bobj

    if (errstat < 0) return

    call l1bread_open_swath (l1bfile, l1bswath, l1bobj, errstat)
    if (errstat < 0) return

    if (size(latr4, 1) /= l1bobj%num_xtrack) then
      call err_message_error ("read_latitude: nxtrack dimension is not correct", errstat)
      call l1bread_close (l1bobj)
      return
    endif

    call l1bread_get2d_r4 (l1bobj, "latitude", tstart, ntimes, latr4, errstat)

    call l1bread_close (l1bobj)

    return

  end subroutine read_latitude

  SUBROUTINE compute_fitting_statistics_nohe5 ( &
      pge_idx, ntimes, nxtrack, xtrange, saocol, saodco, saorms, &
      saofcf, saomqf, errstat )

    USE OMSAO_precision_module,  ONLY: i2, i4, r4, r8
    USE OMSAO_parameters_module, ONLY: &
      i2_missval, r8_missval, main_qa_good, main_qa_suspect, main_qa_bad
    use optimizer_interface_module, only: &
      opt_convergence_failed, opt_convergence_maxiter_exceeded, opt_convergence_suspect, &
      opt_convergence_good
    USE metadata_tools, ONLY:  QAPercentMissingData, QAPercentOutofBoundsData
    USE OMSAO_he5_module,       ONLY:  &
      NrOfInputSamples, NrofGoodOutputSamples, NrofSuspectOutputSamples,        &
      NrofBadOutputSamples, NrofConvergedSamples, NrofFailedConvergenceSamples, &
      NrofExceededIterationsSamples, NrofOutofBoundsSamples, NrofMissingSamples, &
      NrofGoodInputSamples, NrofSuspectOutputSamples, NrofBadOutputSamples,      &
      NrofConvergedSamples, NrofFailedConvergenceSamples, &
      PercentGoodOutputSamples, PercentSuspectOutputSamples, &
      PercentBadOutputSamples, &
      AbsolutePercentMissingSamples
    USE OMSAO_errstat_module,   ONLY: vb_lev_screen, pge_errstat_ok
    USE OMSAO_variables_module, ONLY: verb_thresh_lev, max_good_col
    USE metadata_tools, ONLY: set_automatic_quality_flag

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ntimes, nxtrack
    INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),     INTENT (IN) :: xtrange
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saocol, saodco, saorms
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saofcf

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saomqf

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Local variables
    ! ----------------
    INTEGER (KIND=i4) :: locerrstat, ix, it, spix, epix
    REAL    (KIND=r4) :: PercentOutofBoundsSamples
    REAL    (KIND=r8) :: fitcol_avg, rms_avg, dfitcol_avg, nfitcol
    REAL    (KIND=r8) :: col2sig, col3sig

    locerrstat = pge_errstat_ok

    ! ---------------------------------------------------------
    ! The total number of input samples is simply the number of
    ! pixels in the granule
    ! ---------------------------------------------------------
    NrofInputSamples = nxtrack*ntimes

    ! ------------------------------------------------------------------
    ! Compute all other fitting statistics variables over two nice loops
    ! ------------------------------------------------------------------
    saomqf                        = i2_missval
    NrofGoodInputSamples          = 0_i4
    NrofGoodOutputSamples         = 0_i4
    NrofSuspectOutputSamples      = 0_i4
    NrofBadOutputSamples          = 0_i4
    NrofOutOfBoundsSamples        = 0_i4
    NrofConvergedSamples          = 0_i4
    NrofFailedConvergenceSamples  = 0_i4
    NrofExceededIterationsSamples = 0_i4
    NrofMissingSamples            = 0_i4

    nfitcol    = 0.0_r8
    fitcol_avg = 0.0_r8 ; rms_avg = 0.0_r8 ; dfitcol_avg = 0.0_r8
    DO it = 0, ntimes-1

      spix = xtrange(it,1) ; epix = xtrange(it,2)
      DO ix = spix, epix

        col2sig = saocol(ix,it)+2.0_r8*saodco(ix,it)
        col3sig = saocol(ix,it)+3.0_r8*saodco(ix,it)

        ! ------------------------------------------------------
        ! The Good: Columns are postive within two sigma fitting
        !           uncertainty and the fitting has converged.
        !           For this "sweet spot" we compute the average
        !           fitting statistics.
        ! ------------------------------------------------------
        IF ( (saofcf(ix,it) == opt_convergence_good) .AND. &
            (saocol(ix,it)      >  r8_missval                       ) .AND. &
            (ABS(saocol(ix,it)) <= max_good_col                     ) .AND. &
            (col2sig            >= 0.0_r8                           ) ) THEN

          saomqf(ix,it) = main_qa_good

          NrofGoodInputSamples  = NrofGoodInputSamples   + 1
          NrofConvergedSamples  = NrofConvergedSamples   + 1
          NrofGoodOutputSamples = NrofGoodOutputSamples + 1

          fitcol_avg  = fitcol_avg  + saocol(ix,it)
          dfitcol_avg = dfitcol_avg + saodco(ix,it)
          rms_avg     = rms_avg     + saorms(ix,it)
          nfitcol     = nfitcol     + 1.0_r8
          !write(*,'(1P3E15.5,2I5)') saocol(ix,it),saorms(ix,it), saodco(ix,it), ix, it

          CYCLE
        END IF

        ! ----------------------------------------------------------
        ! The Bad: Fitting hasn't converged or columns are negative
        !          within three sigma fitting uncertainty. Note that
        !          pixels can count towards both the number of out-
        !          of bounds and the failed convergence samples.
        ! ----------------------------------------------------------
        IF ( (saofcf(ix,it)      > i2_missval .AND. saofcf(ix,it) < 0_i2) .OR. &
            (saocol(ix,it)      > r8_missval .AND. col3sig < 0.0_r8    ) ) THEN

          saomqf(ix,it) = main_qa_bad

          NrofGoodInputSamples = NrofGoodInputSamples + 1
          NrofBadOutputSamples = NrofBadOutputSamples + 1

          IF ( saocol(ix,it) > r8_missval .AND. col3sig < 0.0_r8 ) &
            NrofOutofBoundsSamples        = NrofOutofBoundsSamples        + 1
          IF ( saofcf(ix,it) == opt_convergence_failed .or. saofcf(ix,it) == opt_convergence_maxiter_exceeded) &
            NrofFailedConvergenceSamples  = NrofFailedConvergenceSamples  + 1
          IF ( saofcf(ix,it) == opt_convergence_maxiter_exceeded)  &
            NrofExceededIterationsSamples = NrofExceededIterationsSamples + 1

          CYCLE
        END IF

        ! ----------------------------------------------------------
        ! The Ugly: Whatever is left (outside plain missing columns)
        ! ----------------------------------------------------------
        IF ( saocol(ix,it) > r8_missval ) THEN

          IF ( (saofcf(ix,it) == opt_convergence_suspect) .OR. &
              (col2sig <  0.0_r8  .AND. col3sig >= 0.0_r8                      ) .OR. &
              (ABS(saocol(ix,it)) > max_good_col                               ) ) THEN

            saomqf(ix,it) = main_qa_suspect

            NrofGoodInputSamples     = NrofGoodInputSamples     + 1
            NrofSuspectOutputSamples = NrofSuspectOutputSamples + 1

            CYCLE
          END IF

        ELSE

          ! ----------------------------------------------------------
          ! The Missing: Not processed because of either missing input
          !              or restrictions on lat, lon, sza, etc.
          ! ----------------------------------------------------------
          NrofMissingSamples = NrofMissingSamples + 1

        END IF

      END DO
    END DO

    ! --------------------------------------------
    ! Now we can compute averages and percentages,
    ! and write out the final statistics
    ! --------------------------------------------

    IF ( nfitcol >= 1.0_r8 ) THEN
      fitcol_avg  = fitcol_avg  / nfitcol
      rms_avg     = rms_avg     / nfitcol
      dfitcol_avg = dfitcol_avg / nfitcol
    END IF

    PercentGoodOutputSamples      = 100_r4    * &
      REAL(NrofGoodOutputSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples,  KIND=r4) )

    PercentBadOutputSamples       = 100_r4        * &
      REAL(NrofBadOutputSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )

    PercentSuspectOutputSamples   =  100.0_r4         * &
      REAL(NrofSuspectOutputSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )

    PercentOutofBoundsSamples     =  100.0_r4         * &
      REAL(NrofOutofBoundsSamples, KIND=r4) / &
      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )

    AbsolutePercentMissingSamples = 100_r4 * &
      REAL(NrofMissingSamples, KIND=r4) / &
      MAX ( 1.0_4, REAL(NrofInputSamples, KIND=r4) )

    QAPercentMissingData     = NINT ( AbsolutePercentMissingSamples, KIND=i4 )
    QAPercentOutofBoundsData = NINT ( PercentOutofBoundsSamples,     KIND=i4 )

    ! ------------------------------------------------------------------------
    ! With the above information we can easily determine the Automatic QA Flag
    ! ------------------------------------------------------------------------
    CALL set_automatic_quality_flag ( PercentGoodOutputSamples )

    IF ( verb_thresh_lev >= vb_lev_screen ) THEN
      WRITE (*, '(A, 3(1PE15.5))')          'Col-DCol-RMS: ', fitcol_avg, dfitcol_avg, rms_avg
      WRITE (*, '(A, I7,A,I7,A,F7.1,A)')  'Statistics:   ', &
        MAX(NrofGoodOutputSamples,0), ' of ', MAX(NrofGoodInputSamples,0), ' converged - ', &
        MAX(PercentGoodOutputSamples, 0.0), '%'
      WRITE (*, '(A, F7.1)') 'Nfitcol =', nfitcol
    END IF

    errstat = locerrstat
    RETURN

  END SUBROUTINE compute_fitting_statistics_nohe5

END MODULE
