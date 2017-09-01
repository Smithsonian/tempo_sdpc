!> subroutine removed from omi_pge_fitting_aux but kept here for refence
module omi_pge_fitting_aux_disabled

  use omi_pge_fitting_aux

  implicit none

contains

  SUBROUTINE DISABLED_get_input_versions (pge_idx, do_radref, input_versions )

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
  END SUBROUTINE DISABLED_get_input_versions

end module omi_pge_fitting_aux_disabled
