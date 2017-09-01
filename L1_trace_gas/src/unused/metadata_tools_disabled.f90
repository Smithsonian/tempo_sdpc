!> subroutine from metadata_tools module, disabled but kept here for reference
module metadata_tools_diabled

  use metadata_tools

  implicit none

contains

  SUBROUTINE DISABLED_init_metadata (pge_idx, errstat )

    USE OMSAO_precision_module,   ONLY: i4
    USE OMSAO_indices_module,     ONLY: &
      l1b_irradiance_lun, l1b_radiance_lun, l1b_radianceref_lun, &
      mcf_lun, md_inventory_idx, md_archive_idx,                 &
      pge_hcho_idx, pge_gly_idx, pge_h2o_idx,                    &
      voc_amf_luns, voc_omicld_idx, o3_prefit_lun, bro_prefit_lun

    USE OMSAO_errstat_module
    USE OMSAO_parameters_module, ONLY: str_missval, int16_missval, r8_missval
    USE OMSAO_variables_module,  ONLY: l1br_opf_version
    USE OMSAO_prefitcol_module,  ONLY: yn_o3_prefit, yn_bro_prefit
    USE OMSAO_he5_module,        ONLY: &
      granule_day, granule_month, granule_year, TAI93At0zOfGranule, l1b_orbitdata
    USE OMSAO_casestring_module, ONLY: upper_case
    USE sao_pge_utils, ONLY: day_of_year
    IMPLICIT NONE

    !INCLUDE 'PGS_TD_3.f'

    INTEGER (KIND=i4), INTENT(IN) ::  pge_idx
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=14), PARAMETER :: modulename = "init_metadata"

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    CHARACTER (LEN=3)                         :: mdata_typ, mdata_loc
    CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L) :: tmp_string, metadata_type
    INTEGER   (KIND=i4)                       :: &
      imd, version, estat, md_stat, locerrstat, mdata_index, jday

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER (KIND=i4), EXTERNAL :: &
      PGS_TD_UTCtoTAI,      PGS_MET_init, &
      PGS_MET_getSetAttr_s, PGS_MET_getSetAttr_i, PGS_MET_getSetAttr_d, &
      PGS_MET_getPCAttr_s,  PGS_MET_getPCAttr_i,  PGS_MET_getPCAttr_d

    ! ------------------------------------------------------------
    ! Variables for the extraction of the L1B radiance OPF version
    ! ------------------------------------------------------------
    INTEGER   (KIND=i4), PARAMETER                     :: n_ipp = 20
    INTEGER   (KIND=i4)                                :: idx
    CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_ipp) :: l1br_inputp
    CHARACTER (LEN=4)                                  :: opf_string

    ! PGS_MET_getPCAttr_i will fail if LEN is set to anything else
    !CHARACTER (LEN=PGSd_MET_NAME_L),           EXTERNAL :: upper_case

    locerrstat = pge_errstat_ok

    ! initialize this to an invalid value
    mdata_index = HUGE(1_i4)

    ! ----------------------------------------------------------------
    ! Initialize MetaData file (Example taken from PGS Toolkit Primer)
    ! ----------------------------------------------------------------
    mcf_groups = ""
    md_stat = PGS_MET_init ( mcf_lun, mcf_groups )
    CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_fatal, OMSAO_F_METINIT, &
                      modulename, vb_lev_default, errstat )
    IF ( errstat >= pge_errstat_fatal ) RETURN

    ! ---------------------------------
    ! Initialize MetaData STRING fields
    ! ---------------------------------
    DO imd = 1, n_mdata_str

      IF ( TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "mcf" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "pcf" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "l1i" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "l1r" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "l1R"        ) THEN
        version = 1
        mdata_string_values(imd) = str_missval

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_string_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_string_fields(3,imd)))

        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
          mdata_index = md_inventory_idx
          metadata_type = cmd_str//'.0'
        CASE ( "psa" )
          mdata_index   = md_inventory_idx
          metadata_type = cmd_str//'.0'
          ! Do we have any of these at all?
        CASE ( "arc" )
          mdata_index   = md_archive_idx
          metadata_type = amd_str//'.0'
        CASE DEFAULT
          ! No idea what to do here
        END SELECT

        SELECT CASE ( mdata_loc )
        CASE ( "mcf" )
          md_stat = PGS_MET_getSetAttr_s( mcf_groups(mdata_index), &
                                         TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
                                         mdata_string_values(imd) )
          IF ( TRIM(ADJUSTL(mdata_string_fields(1,imd))) == "ShortName" ) &
            mcf_shortname = TRIM(ADJUSTL(mdata_string_values(imd)))
        CASE ( "l1i" )
          md_stat = PGS_MET_getPCAttr_s( &
            l1b_irradiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
            TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
            mdata_string_values(imd) )
        CASE ( "l1r" )
          SELECT CASE ( TRIM(ADJUSTL(mdata_string_fields(2,imd))) )
          CASE ( "arc")
            md_stat = PGS_MET_getPCAttr_s( &
              l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
              TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), mdata_string_values(imd) )
          CASE ( "inv")
            md_stat = PGS_MET_getPCAttr_s( &
              l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
              TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
              mdata_string_values(imd) )
          END SELECT
          IF ( TRIM(ADJUSTL(mdata_string_fields(1,imd))) == "OrbitData" ) &
            l1b_orbitdata = TRIM(ADJUSTL(mdata_string_values(imd)))
          ! -----------------------------------------------------------
          ! Obtain YEAR, MONTH, and DAY, and convert to Day-of-the-Year
          ! -----------------------------------------------------------
          IF ( TRIM(ADJUSTL(mdata_string_fields(1,imd))) == rbd_str ) THEN
            IF ( TRIM(ADJUSTL(mdata_string_values(imd))) /= str_missval ) THEN
              READ ( mdata_string_values(imd), '(I4,1X,I2,1X,I2)') &
                granule_year, granule_month, granule_day
              jday = day_of_year ( granule_year, granule_month, granule_day )

              ! ------------------------------------------------------------
              ! Since we are here, compute Granule Time in seconds; this is
              ! one of the global file attributes that we need to write out.
              ! ------------------------------------------------------------
              tmp_string = ""
              tmp_string = TRIM(ADJUSTL(mdata_string_values(imd))) //"T00:00:00.000Z"
              estat = PGS_TD_UTCtoTAI ( TRIM(ADJUSTL(tmp_string)), TAI93At0zOfGranule )
              IF ( estat /= PGS_S_SUCCESS .AND. estat /= PGSTD_E_NO_LEAP_SECS ) &
                CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_TAI93, &
                                  modulename, vb_lev_default, errstat )
              ! ------------------------------------------------------------
            ELSE
              granule_year = 0 ; granule_month = 0 ; granule_day = 0 ; jday = 0
            END IF
          END IF

        CASE ( "l1R" ) ! Radiance reference granule
          SELECT CASE ( TRIM(ADJUSTL(mdata_string_fields(2,imd))) )
          CASE ( "inv")
            md_stat = PGS_MET_getPCAttr_s( &
              l1b_radianceref_lun, version, TRIM(ADJUSTL(metadata_type)), &
              TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
              mdata_string_values(imd) )
          CASE DEFAULT
            ! Nothing to do here
          END SELECT

        CASE DEFAULT
          ! Nothing to do here
        END SELECT

        ! ------------------------------
        ! Error check for initialization
        ! ------------------------------
        CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                          modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,imd))), &
                          vb_lev_default, errstat )

      END IF

    END DO

    ! ----------------------------------------------------
    ! Initialize MetaData STRING fields special for OMHCHO
    ! ----------------------------------------------------
    SELECT CASE ( pge_idx )
    CASE (  pge_hcho_idx )
      DO imd = 1, n_mdata_omhcho
        version = 1
        mdata_omhcho_values(imd) = str_missval

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_omhcho_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_omhcho_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "BRO" )
          IF ( yn_bro_prefit(1) ) THEN
            md_stat = PGS_MET_getPCAttr_s( &
              bro_prefit_lun, version, cmd_str, &
              TRIM(ADJUSTL(upper_case(mdata_omhcho_fields(1,imd)))), &
              mdata_omhcho_values(imd) )
            CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                              modulename//f_sep//TRIM(ADJUSTL(mdata_omhcho_fields(1,imd))), &
                              vb_lev_default, errstat )
          END IF
        CASE ( "OOO" )
          IF ( yn_o3_prefit(1) ) THEN
            md_stat = PGS_MET_getPCAttr_s( &
              o3_prefit_lun, version, cmd_str, &
              TRIM(ADJUSTL(upper_case(mdata_omhcho_fields(1,imd)))), &
              mdata_omhcho_values(imd) )
            CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                              modulename//f_sep//TRIM(ADJUSTL(mdata_omhcho_fields(1,imd))), &
                              vb_lev_default, errstat )
          END IF
        CASE DEFAULT
          ! No idea what to do here
        END SELECT
      END DO

      DO imd = 1, n_mdata_voc
        version = 1
        mdata_voc_values(imd) = str_missval

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_voc_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_voc_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "CLD" )
          md_stat = PGS_MET_getPCAttr_s( &
            voc_amf_luns(voc_omicld_idx), version, TRIM(ADJUSTL(cmd_str))//'.0', &
            TRIM(ADJUSTL(upper_case(mdata_voc_fields(1,imd)))), &
            mdata_voc_values(imd) )
          CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                            modulename//f_sep//TRIM(ADJUSTL(mdata_voc_fields(1,imd))), &
                            vb_lev_default, errstat )
        CASE DEFAULT
          ! No idea what to do here
        END SELECT
      END DO

      ! ------------------------------------------------------
      ! Initialize MetaData STRING fields special for OMCHOCHO
      ! ------------------------------------------------------
    CASE ( pge_gly_idx)
      DO imd = 1, n_mdata_omchocho
        version = 1
        mdata_voc_values(imd) = str_missval

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_voc_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_voc_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "CLD" )
          md_stat = PGS_MET_getPCAttr_s( &
            voc_amf_luns(voc_omicld_idx), version, TRIM(ADJUSTL(cmd_str))//'.0', &
            TRIM(ADJUSTL(upper_case(mdata_voc_fields(1,imd)))), &
            mdata_voc_values(imd) )
          CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                            modulename//f_sep//TRIM(ADJUSTL(mdata_voc_fields(1,imd))), &
                            vb_lev_default, errstat )
        CASE DEFAULT
          ! No idea what to do here
        END SELECT
      END DO
      ! ----------------------------------------------------------------------------------------------
      ! Initialize MetaData STRING fields special for OMH2O same that for OMHCHO (not completely true)
      ! ----------------------------------------------------------------------------------------------
    CASE ( pge_h2o_idx)
      DO imd = 1, n_mdata_omhcho
        version = 1
        mdata_voc_values(imd) = str_missval

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_voc_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_voc_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "CLD" )
          md_stat = PGS_MET_getPCAttr_s( &
            voc_amf_luns(voc_omicld_idx), version, TRIM(ADJUSTL(cmd_str))//'.0', &
            TRIM(ADJUSTL(upper_case(mdata_voc_fields(1,imd)))), &
            mdata_voc_values(imd) )
          CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                            modulename//f_sep//TRIM(ADJUSTL(mdata_voc_fields(1,imd))), &
                            vb_lev_default, errstat )
        CASE DEFAULT
          ! No idea what to do here
        END SELECT
      END DO

    END SELECT

    ! ----------------------------------
    ! Initialize MetaData INTEGER fields
    ! ----------------------------------
    DO imd = 1, n_mdata_int
      IF ( TRIM(ADJUSTL(mdata_integer_fields(3,imd))) == "pcf" .OR. &
          TRIM(ADJUSTL(mdata_integer_fields(3,imd))) == "mcf" .OR. &
          TRIM(ADJUSTL(mdata_integer_fields(3,imd))) == "l1r"        ) THEN

        version = 1
        mdata_integer_values(imd) = int16_missval

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_integer_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_integer_fields(3,imd)))

        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
          mdata_index   = md_inventory_idx
          metadata_type = cmd_str//'.0'
        CASE ( "psa" )
          mdata_index   = md_inventory_idx
          metadata_type = cmd_str//'.0'
          ! Do we have any of these at all?
        CASE ( "arc" )
          mdata_index   = md_archive_idx
          metadata_type = amd_str
        CASE DEFAULT
          ! No idea what to do here
        END SELECT

        SELECT CASE ( mdata_loc )
        CASE ( "mcf" )
          md_stat = PGS_MET_getSetAttr_i( mcf_groups(mdata_index), &
                                         TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
                                         mdata_integer_values(imd) )
        CASE ( "pcf" )
          md_stat = PGS_MET_getSetAttr_i( mcf_groups(mdata_index), &
                                         TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
                                         mdata_integer_values(imd) )
          IF ( TRIM(ADJUSTL(mdata_integer_fields(1,imd))) == vid_str ) &
            mcf_versionid = mdata_integer_values(imd)
        CASE ( "l1r" )
          md_stat = PGS_MET_getPCAttr_i( &
            l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
            TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
            mdata_integer_values(imd) )
        CASE ( "l1R" )
          md_stat = PGS_MET_getPCAttr_i( &
            l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
            TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
            mdata_integer_values(imd) )
        CASE DEFAULT
          ! No idea what to do here
        END SELECT

        ! ------------------------------
        ! Error check for initialization
        ! ------------------------------
        CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_GETATTR, &
          modulename//f_sep//TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
          vb_lev_default, errstat )

      END IF
    END DO

    ! ----------------------------------
    ! Initialize MetaData DOUBLE fields
    ! ----------------------------------
    DO imd = 1, n_mdata_dbl
      IF ( TRIM(ADJUSTL(mdata_double_fields(3,imd))) == "mcf" .OR. &
          TRIM(ADJUSTL(mdata_double_fields(3,imd))) == "l1r"        ) THEN

        version = 1
        mdata_double_values(imd) = r8_missval

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_double_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_double_fields(3,imd)))

        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
          mdata_index   = md_inventory_idx
          metadata_type = cmd_str//'.0'
        CASE ( "psa" )
          mdata_index   = md_inventory_idx
          metadata_type = cmd_str//'.0'
          ! Do we have any of these at all?
        CASE ( "arc" )
          mdata_index   = md_archive_idx
          metadata_type = amd_str//'.0'
        CASE DEFAULT
          ! No idea what to do here
        END SELECT

        SELECT CASE ( mdata_loc )
        CASE ( "mcf" )
          md_stat = PGS_MET_getSetAttr_d( mcf_groups(mdata_index), &
                                         TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
                                         mdata_double_values(imd) )
        CASE ( "l1r" )
          md_stat = PGS_MET_getPCAttr_d( &
            l1b_radiance_lun, version, TRIM(ADjUSTL(metadata_type)), &
            TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
            mdata_double_values(imd) )
        CASE ( "l1R" )
          md_stat = PGS_MET_getPCAttr_d( &
            l1b_radiance_lun, version, TRIM(ADjUSTL(metadata_type)), &
            TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
            mdata_double_values(imd) )
        CASE DEFAULT
          ! No idea what to do here
        END SELECT

        ! ------------------------------
        ! Error check for initialization
        ! ------------------------------
        CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_GETATTR, &
          modulename//f_sep//TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
          vb_lev_default, errstat )

      END IF
    END DO

    ! -----------------------------
    ! Read L1B radiance OPF version
    ! -----------------------------
    ! (somehow this fails when called as subroutine)
    ! ----------------------------------------------
    l1br_inputp="" !JED
    md_stat = PGS_MET_getPCAttr_s ( &
      l1b_radiance_lun, version, TRIM(ADJUSTL(cmd_str//'.0')), 'INPUTPOINTER', l1br_inputp )
    l1br_opf_version = -1 ; version = 1
    get_opf: DO imd = 1, n_ipp
      tmp_string = l1br_inputp(imd)
      idx        = INDEX ( tmp_string, 'OML1BOPF' )
      IF ( idx > 0 ) THEN
        opf_string = tmp_string(idx+10:idx+14)
        READ ( opf_string, '(I4)') l1br_opf_version
      END IF
      IF ( l1br_opf_version > 0 ) EXIT get_opf
    END DO get_opf

    RETURN
  END SUBROUTINE DISABLED_init_metadata

end module metadata_tools_diabled
