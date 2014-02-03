MODULE he5_output_tools

  USE OMSAO_indices_module, ONLY: n_voc_amf_luns
  USE OMSAO_precision_module, ONLY: i2, i4, r4, r8
  USE OMSAO_parameters_module, ONLY: MAX_STR_LEN

  use ctrlvars, only: yn_diagnostic_run, yn_refseccor, yn_scat_weights

  USE OMSAO_he5_module, ONLY: he5_start_2d, he5_stride_2d, he5_edge_2d, &
    he5_start_3d, he5_stride_3d, he5_edge_3d, &
    HE5T_NATIVE_DOUBLE, HE5T_NATIVE_FLOAT, HE5T_NATIVE_INT, &
    HE5T_NATIVE_INT16, HE5T_NATIVE_CHAR, &
    pge_swath_id, vcoordinate_field, missval_attr, offset_attr, &
    rstemp_attr, scafac_attr, title_attr, ufd_attr, units_attr, &
    valids_attr, brofit_swath_file_id, brofit_swath_id, &
    brofit_swath_name, lqh2ofit_swath_file_id, lqh2ofit_swath_id, &
    o3fit_swath_file_id, o3fit_swath_id, vertical_coordinate, &
    pge_swath_file_id, lqh2ofit_swath_name, o3fit_swath_name, granule_day, &
    granule_month, granule_year, input_versions, l1b_orbitdata, &
    NrofBadOutputSamples, NrofConvergedSamples, NrofCrossTrackPixels, &
    NrofFailedConvergenceSamples, NrofExceededIterationsSamples, &
    NrofGoodInputSamples, NrofGoodOutputSamples, NrofOutofBoundsSamples, &
    NrofSuspectOutputSamples, PercentBadOutputSamples, NrofInputSamples, &
    NrofScanLines, PercentGoodOutputSamples, PercentSuspectOutputSamples, &
    TAI93At0zOfGranule, tai_attr, &
    HE5_SWrdfld, HE5_SWwrlattr, HE5_SWwrattr, HE5_SWsetfill, &
    HE5_SWdetach, HE5_SWclose, HE5_EHwrglatt, HE5_SWwrfld

  USE iso_c_binding, only: c_loc, c_ptr

  CHARACTER (LEN=MAX_STR_LEN), DIMENSION (n_voc_amf_luns) :: amf_swath_names    = 'undefined'
  INTEGER   (KIND=i4),      DIMENSION (n_voc_amf_luns) :: amf_swath_ids      = -1
  INTEGER   (KIND=i4),      DIMENSION (n_voc_amf_luns) :: amf_swath_file_ids = -1

  type, PRIVATE :: Array_Container_Type
    character (len=2) :: dtype
    INTEGER (kind=i2), DIMENSION(:), pointer :: array_i2
    INTEGER (kind=i4), DIMENSION(:), pointer :: array_i4
    REAL (kind=r4), DIMENSION(:), pointer :: array_r4
    REAL (kind=r8), DIMENSION(:), pointer :: array_r8
  end type Array_Container_Type

  PRIVATE
  PUBLIC he5_write_fitting_statistics, he5_write_omi_database, he5_write_radfit_output, &
    saopge_geofield_read, saopge_columninfo_read, he5_init_swath, &
    he5_define_fields, he5_close_output_file, he5_set_field_attributes, &
    he5_write_swath_attributes, he5_write_global_attributes, &
    he5_write_wavcal_output, he5_write_common_mode

CONTAINS

  FUNCTION he5_init_swath ( file_name, swath_name, nTimes, nXtrack, nSwLevels ) RESULT ( he5stat )

    !------------------------------------------------------------------------------
    ! This function initializes the HE5 output swath.
    !
    ! Input:
    !   file_name  - Name of HE5 output file
    !   swath_name - Name of swath to be created
    !   nXtrack    - Number of cross-track positions (mind spatial zoom mode!)
    !
    ! Return: he5stat
    !
    ! Variables passed through MODULE:
    !   pge_swath_file_id - id number for HE5 output file (required for closing it)
    !   pge_swath_id      - id number for swath (required for writing to swath)
    !
    !------------------------------------------------------------------------------

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: NWAVEL_MAX, nUTCdim
    USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
    USE OMSAO_omidata_module,   ONLY: nclenfit, n_comm_wvl
    USE OMSAO_he5_module
    USE OMSAO_errstat_module
    USE OMSAO_variables_module, ONLY: n_fitvar_rad

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=14), PARAMETER :: modulename = 'he5_init_swath'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4),      INTENT(IN) :: nTimes, nXtrack, nSwLevels
    CHARACTER (LEN=MAX_STR_LEN), INTENT(IN) :: file_name, swath_name

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: he5stat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=C_LONG), PARAMETER :: onecl = 1, twocl = 2, fourcl = 4
    INTEGER   (KIND=i4) :: errstat

    he5stat = pge_errstat_ok
    errstat = pge_errstat_ok

    ! ---------------------------------------------------------------
    ! Open HE5 output file and check PGE_SWATH_FILE_ID ( -1 if error)
    ! ---------------------------------------------------------------
    pge_swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_trunc )
    IF ( pge_swath_file_id == he5_stat_fail ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_fatal, OMSAO_F_HE5SWOPEN, modulename, vb_lev_default, he5stat )
      RETURN
    END IF

    ! ------------------------------------------------------
    ! Create HE5 swath and check PGE_SWATH_ID ( -1 if error)
    ! ------------------------------------------------------
    pge_swath_id = HE5_SWcreate ( pge_swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( pge_swath_id == he5_stat_fail ) THEN
      CALL error_check ( &
        0, 1, pge_errstat_fatal, OMSAO_F_HE5SWCREATE, modulename, vb_lev_default, he5stat )
      RETURN
    END IF

    ! ---------------------------------------------------------
    ! Define new dimensions in HE5 swath and check error status
    ! ---------------------------------------------------------
    errstat = HE5_SWdefdim  ( pge_swath_id, ntc,        INT(nTimes,    KIND=C_LONG) )
    errstat = HE5_SWdefdim  ( pge_swath_id, nxc,        INT(nXtrack,   KIND=C_LONG) )
    errstat = HE5_SWdefdim  ( pge_swath_id, nlc,        INT(nSwLevels, KIND=C_LONG) ) !gga
    errstat = HE5_SWdefdim  ( pge_swath_id, ntcp1,      INT(nTimes+1,  KIND=C_LONG) )
    errstat = HE5_SWdefdim  ( pge_swath_id, nxcp1,      INT(nXtrack+1, KIND=C_LONG) )
    errstat = HE5_SWDefdim  ( pge_swath_id, nutcd,      INT(nUTCdim,   KIND=C_LONG) )
    errstat = HE5_SWDefdim  ( pge_swath_id, "1",        onecl                       )
    errstat = HE5_SWDefdim  ( pge_swath_id, "2",        twocl                       )
    errstat = HE5_SWDefdim  ( pge_swath_id, "4",        fourcl                      )

    ! ---------------------------------------------------------
    ! Dimensions for Diagnostic Fields
    ! ---------------------------------------------------------
    IF ( yn_diagnostic_run ) THEN
      errstat = HE5_SWdefdim  ( pge_swath_id, nfv,   INT(n_fitvar_rad,   KIND=C_LONG) )
      errstat = HE5_SWdefdim  ( pge_swath_id, ncwvl, INT(n_comm_wvl,     KIND=C_LONG) )
      errstat = HE5_SWDefdim  ( pge_swath_id, ncv,   INT(nclenfit,       KIND=C_LONG) )
      errstat = HE5_SWdefdim  ( pge_swath_id, nwcp,  INT(max_calfit_idx, KIND=C_LONG) )

      ! CCM for refspec database
      errstat = HE5_SWdefdim  ( pge_swath_id, nrspc, INT(max_rs_idx,     KIND=C_LONG) )
      errstat = HE5_SWdefdim  ( pge_swath_id, nwalm, INT(nwavel_max,     KIND=C_LONG) )

    END IF

    CALL error_check ( &
      errstat, he5_stat_ok, pge_errstat_fatal, OMSAO_F_HE5SWDEFDIM, &
      modulename, vb_lev_default, he5stat )
    IF ( he5stat >= pge_errstat_error ) RETURN

    RETURN
  END FUNCTION he5_init_swath

  SUBROUTINE define_fields (swid, use_dfld, datafield_list, &
                            ntimes, nxtrack, nswlevels, errstat)

    USE OMSAO_precision_module, ONLY: C_LONG
    USE datafields, ONLY: DataField_HE5, Datafield_List_Type
    USE OMSAO_omidata_module, ONLY: n_field_maxdim
    USE OMSAO_he5_module, ONLY: he5_comp_type, he5_hdfe_nomerge, &
      he5_comp_par, he5_nocomp_par, he5_nocomp_type, &
      HE5_SWdefcomch, HE5_SWdefchunk, HE5_SWdefcomp, &
      HE5_SWdefgfld, HE5_SWdefdfld

    IMPLICIT NONE
    INTEGER (KIND=i4), INTENT(IN) :: swid, use_dfld, ntimes, &
      nxtrack, nswlevels
    INTEGER (KIND=i4), INTENT(INOUT) :: errstat
    INTEGER (KIND=i4), EXTERNAL :: deffunc
    type (Datafield_List_Type), intent(in) :: datafield_list
    LOGICAL :: compress_ok
    INTEGER (KIND=i4) :: n_chunk_dim
    INTEGER (KIND=C_LONG), DIMENSION(n_field_maxdim) :: chunk_dim
    TYPE (DataField_HE5), pointer :: ptr

    ptr => datafield_list%head
    do while (associated (ptr))
      if (.not. ptr%output) cycle

      CALL he5_check_for_compressibility ( &
        ntimes, nxtrack, nSwLevels, trim(ptr%dimensions), &
        compress_ok, n_chunk_dim, chunk_dim )

      IF (compress_ok) THEN
        errstat = HE5_SWdefcomch ( &
          swid, he5_comp_type, he5_comp_par, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
      ELSE
        errstat = HE5_SWdefchunk(swid, n_chunk_dim, chunk_dim(1:n_chunk_dim) )
        errstat = HE5_SWdefcomp (swid, he5_nocomp_type, he5_nocomp_par )
      END IF

      ptr%swath_id = swid
      CALL he5_set_fill_value (ptr, errstat )

      if (use_dfld /= 0) then
        errstat = HE5_SWdefdfld (swid, &
                                 trim(ptr%name), trim(ptr%dimensions), &
                                 " ", ptr%he5_datatype, &
                                 he5_hdfe_nomerge);
      else
        errstat = HE5_SWdefgfld (swid, &
                                 trim(ptr%name), trim(ptr%dimensions), &
                                 " ", ptr%he5_datatype, &
                                 he5_hdfe_nomerge);
      endif
      ptr => ptr%next
    end do
  END SUBROUTINE

  FUNCTION he5_define_fields ( pge_idx, swath_name, nTimes, nXtrack, nSwLevels ) RESULT ( he5stat )

    !------------------------------------------------------------------------------
    ! This function defines fields in the HE5 output swath.
    !
    ! Input:
    !   pge_idx    - Index of current PGE
    !
    ! Return: he5stat
    !
    ! Variables passed through MODULE:
    !   pge_swath_file_id - id number for HE5 output file (required for closing it)
    !   pge_swath_id      - id number for swath (required for writing to swath)
    !
    !------------------------------------------------------------------------------

    USE OMSAO_indices_module,    ONLY: &
      pge_hcho_idx, pge_gly_idx, pge_o3_idx, sao_molecule_names, &
      pge_h2o_idx
    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN
    USE OMSAO_he5_module
    USE datafields, ONLY: comdata_he5fields, diagnostic_he5fields, &
      geo_he5fields, o3_prefit_he5fields, o3_prefit_uncert_he5fields, &
      rad_calfit_he5fields, rad_reffit_he5fields, rs_he5fields, &
      sol_calfit_he5fields, sw_he5fields, voc_he5fields
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=17), PARAMETER :: modulename = 'he5_define_fields'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4),      INTENT(IN) :: pge_idx, nTimes, nXtrack, nSwLevels
    CHARACTER (LEN=MAX_STR_LEN), INTENT(IN) :: swath_name

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: he5stat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                             :: errstat
    CHARACTER (LEN=4)                               :: molstr

    ! ---------------------------------------------
    ! Variables related to compression and chunking
    ! ---------------------------------------------
    he5stat = pge_errstat_ok
    errstat = pge_errstat_ok

    ! ----------------------------------------------------
    ! Some fields have the PGE molecule name added to them
    ! ----------------------------------------------------
    molstr = TRIM(ADJUSTL(sao_molecule_names(pge_idx)))

    ! -------------------------------------------------------------
    ! Define geolocation fields in HE5 swath and check error status
    ! -------------------------------------------------------------
    call define_fields (pge_swath_id, 0, geo_he5fields, &
                        ntimes, nxtrack, nswlevels, errstat)

    ! ----------------------------------------------------------------------
    ! Create fields for solar and radiance wavelength calibration parameters
    ! ----------------------------------------------------------------------
    call define_fields (pge_swath_id, 1, sol_calfit_he5fields, &
                        ntimes, nxtrack, nswlevels, errstat)

    call define_fields (pge_swath_id, 1, rad_calfit_he5fields, &
                        ntimes, nxtrack, nswlevels, errstat)

    call define_fields (pge_swath_id, 1, rad_reffit_he5fields, &
                        ntimes, nxtrack, nswlevels, errstat)

    ! --------------------
    ! Radiance data fields
    ! --------------------

    call define_fields (pge_swath_id, 1, comdata_he5fields, &
                        ntimes, nxtrack, nswlevels, errstat)

    IF ( yn_diagnostic_run ) &
      call define_fields (pge_swath_id, 1, diagnostic_he5fields, &
                          ntimes, nxtrack, nswlevels, errstat)

    ! -------------------------------------------------------------------------
    ! (3) PGE-Specific Data Fields - Currently OMHCHO, OMCHOCHO and OMSAO3 only
    ! -------------------------------------------------------------------------
    SELECT CASE ( pge_idx )
    CASE ( pge_hcho_idx )
      ! ---------------------------------------------
      ! Special HCHO Data Fields - AMF Cloud Data
      ! ---------------------------------------------

      call define_fields (pge_swath_id, 1, voc_he5fields, &
                          ntimes, nxtrack, nswlevels, errstat)

      ! ------------------------------------------------------
      ! Reference Sector Correction 2 fields: Total Column and
      !   uncertainty but so far by limiting the number of loo
      !   ps only 1 is used Total Column gga
      ! ------------------------------------------------------
      IF (yn_refseccor) &
        call define_fields (pge_swath_id, 1, rs_he5fields, &
                            ntimes, nxtrack, nswlevels, errstat)

    CASE ( pge_gly_idx)
      ! ---------------------------------------------------------------------
      ! Special CHOCHO Data Fields - AMF Cloud Data
      ! ---------------------------------------------------------------------
      call define_fields (pge_swath_id, 1, voc_he5fields, &
                          ntimes, nxtrack, nswlevels, errstat)

    CASE ( pge_h2o_idx)
      ! ---------------------------------------------------------------------
      ! Special H2O Data Fields - AMF Cloud Data
      ! ---------------------------------------------------------------------
      call define_fields (pge_swath_id, 1, voc_he5fields, &
                          ntimes, nxtrack, nswlevels, errstat)

    CASE ( pge_o3_idx )
      ! ----------------------------------
      ! Special Data Fields for OMSAO3
      ! ----------------------------------
      call define_fields (pge_swath_id, 1, o3_prefit_he5fields, &
                          ntimes, nxtrack, nswlevels, errstat)
      call define_fields (pge_swath_id, 1, o3_prefit_uncert_he5fields, &
                          ntimes, nxtrack, nswlevels, errstat)
    CASE DEFAULT
      ! --------------------------------------------------------------------
      ! Nothing special for any other PGE yet
      ! (although OMCHOCHO will require the same AMF Cloud fields as OMHCHO)
      ! --------------------------------------------------------------------
    END SELECT

    ! -------------------------------------------------
    ! (4) Scattering weights, gas profile averaging kernels
    ! and albedo gga
    ! -------------------------------------------------
    IF (yn_scat_weights) &
      call define_fields (pge_swath_id, 1, sw_he5fields, &
                          ntimes, nxtrack, nswlevels, errstat)

    ! ------------------------------------------
    ! Check error status of swath initialization
    ! ------------------------------------------
    CALL error_check ( &
      errstat, he5_stat_ok, pge_errstat_error, OMSAO_E_HE5SWDEFFLD, &
      modulename, vb_lev_default, he5stat )

    ! -------------------------------------------------------------------------------
    ! Detach from and re-attach to created swath (recommended before adding to swath)
    ! -------------------------------------------------------------------------------
    errstat       = HE5_SWdetach ( pge_swath_id )
    pge_swath_id  = HE5_SWattach ( pge_swath_file_id, TRIM(ADJUSTL(swath_name)) )

    ! JED NOTE: Do we really want to ignore the error here?
    IF ( pge_swath_id == he5_stat_fail ) CALL error_check ( &
      0, 1, pge_errstat_fatal, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )

    RETURN

  END FUNCTION he5_define_fields

  function lookup_field (list, fieldname) result (ptr)

    use datafields, ONLY: Datafield_List_Type, DataField_HE5

    implicit none
    type (Datafield_List_Type), intent(in) :: list
    character (len=*), intent(in) :: fieldname
    type (DataField_HE5), pointer :: ptr

    ptr => list%head
    do while (associated (ptr))
      if (ptr%name == fieldname) return
      ptr => ptr%next
    end do
    return
  end function lookup_field

  ! ---------------------------------------------------------------------------
  subroutine do_he5_swwrfld (list, fieldname, start, stride, edge, &
                             array_struct, errstat)

    use OMSAO_precision_module, ONLY: C_LONG
    use datafields, ONLY: Datafield_List_Type, DataField_HE5
    use OMSAO_errstat_module, ONLY: pge_errstat_fatal
    use OMSAO_he5_module, ONLY: he5_swwrfld

    implicit none
    type (Datafield_List_Type), intent(in) :: list
    character (len=*), intent (in) :: fieldname
    integer (kind=C_LONG), dimension(:), intent(in) :: start, stride, edge
    integer (kind=i4), intent (inout) :: errstat
    type (Array_Container_Type), intent(in) :: array_struct

    ! locals
    type (DataField_HE5), pointer :: ptr
    integer (kind=i4) :: locerr

    ptr => lookup_field (list, fieldname)
    if (.not.associated (ptr)) then
      write (0, *) '***FATAL ERROR: in do_he5_swwrfld: field name "', &
        trim(fieldname), '" not found'
      errstat = pge_errstat_fatal
      return
    endif

    select case (array_struct%dtype)

    case ("i2")
      locerr = he5_swwrfld (ptr%swath_id, TRIM(fieldname), &
                            start, stride, edge, array_struct%array_i2)
    case ("i4")
      locerr = he5_swwrfld (ptr%swath_id, TRIM(fieldname), &
                            start, stride, edge, array_struct%array_i4)
    case ("r4")
      locerr = he5_swwrfld (ptr%swath_id, TRIM(fieldname), &
                            start, stride, edge, array_struct%array_r4)
    case ("r8")
      locerr = he5_swwrfld (ptr%swath_id, TRIM(fieldname), &
                            start, stride, edge, array_struct%array_r8)
    case default
      locerr = pge_errstat_fatal

    end select

    ! FIXME: there are 2 classes of errors, the pge ones and he5
    errstat = locerr

  end subroutine do_he5_swwrfld

  SUBROUTINE he5_write_wavcal_output ( nXtloc, fpix, lpix, errstat )

    USE OMSAO_precision_module, ONLY: i4, r8, i8, C_LONG
    USE OMSAO_indices_module,   ONLY: max_calfit_idx
    USE OMSAO_errstat_module
    USE OMSAO_omidata_module,   ONLY: &
      n_roff_dig,                                            &
      omi_solcal_xflag,  omi_radcal_xflag, omi_radref_xflag, &
      omi_solcal_pars,   omi_radcal_pars,  omi_radref_pars,  &
      omi_radref_col,    omi_radref_dcol,  omi_radref_rms,   &
      omi_radref_xtrcol
    USE datafields, ONLY: rad_reffit_he5fields, &
      rrcf_field, rrcol_field, rrdcol_field, rrlr_field, rrrms_field, &
      rrxcol_field, rwccf_field, rwclr_field, sol_calfit_he5fields, &
      swccf_field, rad_calfit_he5fields

    USE OMSAO_variables_module, ONLY: radref_latrange
    USE sao_pge_utils, ONLY: roundoff_2darr_r4, roundoff_1darr_r8
    !use addr_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_write_wavcal_output'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nXtloc, fpix, lpix

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                             :: locerrstat, i, j, npix
    REAL    (KIND=r8), DIMENSION (max_calfit_idx) :: tmpr8
    type (Array_Container_Type) :: array_struct

    locerrstat = pge_errstat_ok

    ! -----------------------------------------------
    ! Number of cross-track pixels actually processed
    ! -----------------------------------------------
    npix = lpix - fpix + 1

    ! -------------------------------------------------------------------------------------------
    ! Write results for solar and radiance wavelength calibration, and the radiance reference fit
    ! -------------------------------------------------------------------------------------------
    DO j = fpix, lpix
      tmpr8(1:max_calfit_idx) = omi_solcal_pars(1:max_calfit_idx,j)
      CALL roundoff_1darr_r8 ( n_roff_dig, max_calfit_idx, tmpr8(1:max_calfit_idx) )
      omi_solcal_pars(1:max_calfit_idx,j) = tmpr8(1:max_calfit_idx)

      tmpr8(1:max_calfit_idx) = omi_radcal_pars(1:max_calfit_idx,j)
      CALL roundoff_1darr_r8 ( n_roff_dig, max_calfit_idx, tmpr8(1:max_calfit_idx) )
      omi_radcal_pars(1:max_calfit_idx,j) = tmpr8(1:max_calfit_idx)

      tmpr8(1:max_calfit_idx) = omi_radref_pars(1:max_calfit_idx,j)
      CALL roundoff_1darr_r8 ( n_roff_dig, max_calfit_idx, tmpr8(1:max_calfit_idx) )
      omi_radref_pars(1:max_calfit_idx,j) = tmpr8(1:max_calfit_idx)
    END DO

    ! --------------------------------------------------
    ! There is only one variable in the data field array
    ! --------------------------------------------------
    i = 1

    ! ------------------------------------------
    ! The remaining parameters are 1-dimensional
    ! ------------------------------------------
    he5_start_2d  = (/                0, 0 /)
    he5_stride_2d = (/                1, 0 /)
    he5_edge_2d   = (/ INT(nXtloc,KIND=C_LONG), 0_i8 /)

    ! -----------------------------
    ! Exit value of Fitting Routine
    ! -----------------------------
    array_struct%array_i2 => omi_solcal_xflag; array_struct%dtype = "i2"
    call do_he5_swwrfld (sol_calfit_he5fields, swccf_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, &
                         locerrstat);

    array_struct%array_i2 => omi_radcal_xflag; array_struct%dtype = "i2"
    call do_he5_swwrfld (rad_calfit_he5fields, rwccf_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    array_struct%array_i2 => omi_radref_xflag; array_struct%dtype = "i2"
    call do_he5_swwrfld (rad_reffit_he5fields, rrcf_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    ! --------------------------------------------------------------------------
    ! Special output for radiance wavelength calibration, and the possibility of
    ! the usage of a radiance reference spectrum.
    ! --------------------------------------------------------------------------
    he5_start_2d  = (/ 0, 0 /)
    he5_stride_2d = (/ 1, 0 /)
    he5_edge_2d   = (/ 2, 0 /)

    ! FIXME: The original code used radref_latrange for both rwccf_field and rrlr_field
    ! Are they always the same??
    array_struct%array_r4 => radref_latrange; array_struct%dtype = "r4"
    call do_he5_swwrfld (rad_calfit_he5fields, rwclr_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    array_struct%array_r4 => radref_latrange; array_struct%dtype = "r4"
    call do_he5_swwrfld (rad_reffit_he5fields, rrlr_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    ! ---------------------------------------------------------------------------
    ! Write results for column and column uncertainty from radiance reference fit
    ! ---------------------------------------------------------------------------
    CALL roundoff_1darr_r8 ( n_roff_dig, npix, omi_radref_col   (fpix:lpix) )
    CALL roundoff_1darr_r8 ( n_roff_dig, npix, omi_radref_dcol  (fpix:lpix) )
    CALL roundoff_1darr_r8 ( n_roff_dig, npix, omi_radref_rms   (fpix:lpix) )
    CALL roundoff_1darr_r8 ( n_roff_dig, npix, omi_radref_xtrcol(fpix:lpix) )

    he5_start_2d  = (/ 0, 0 /)
    he5_stride_2d = (/ 1, 0 /)
    he5_edge_2d   = (/ INT(nXtloc,KIND=C_LONG), INT(0, KIND=C_LONG) /)

    array_struct%array_r8 => omi_radref_col; array_struct%dtype = "r8"
    call do_he5_swwrfld (rad_reffit_he5fields, rrcol_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    array_struct%array_r8 => omi_radref_dcol; array_struct%dtype = "r8"
    call do_he5_swwrfld (rad_reffit_he5fields, rrdcol_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    array_struct%array_r8 => omi_radref_xtrcol; array_struct%dtype = "r8"
    call do_he5_swwrfld (rad_reffit_he5fields, rrxcol_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    array_struct%array_r8 => omi_radref_rms; array_struct%dtype = "r8"
    call do_he5_swwrfld (rad_reffit_he5fields, rrrms_field, &
                         he5_start_2d, he5_stride_2d, he5_edge_2d, &
                         array_struct, locerrstat);

    ! ------------------
    ! Check error status
    ! ------------------

    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
                      modulename, vb_lev_default, errstat )

    RETURN
  END SUBROUTINE he5_write_wavcal_output

  !---------------------------------------------------------------------------
  SUBROUTINE he5_write_radfit_output ( &
      pge_idx, iline, nXtrack, nblock, fpix, lpix, &
      all_fitted_columns, all_fitted_errors, correlation_columns,&
      omi_fitspc, errstat )

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: nxtrack_max, nUTCdim
    USE OMSAO_variables_module,  ONLY: n_fitvar_rad, n_rad_wvl, n_rad_wvl_max
    USE OMSAO_indices_module,    ONLY: sao_molecule_names, &
      pge_bro_idx, pge_o3_idx, o3_t1_idx
!      pge_bro_idx, pge_o3_idx, o3_t1_idx,
!      corr_didx,  corrcol_didx, correrr_didx, itnum_didx,  &
!      fitwt_didx, posobs_didx,  spcobs_didx,  spcfit_didx, &
!      spcres_didx
    USE OMSAO_omidata_module
    USE OMSAO_errstat_module
    USE sao_pge_utils, ONLY: roundoff_2darr_r4, roundoff_1darr_r8, roundoff_2darr_r8
    USE datafields, ONLY: DataField_HE5, &
      o3_prefit_he5fields, o3_prefit_uncert_he5fields, &
      auraalt_field, col_field, corr_field, corrcol_field, correrr_field, &
      dcol_field, extr_field, lat_field, lon_field, posobs_field, &
      fitcon_field, fitrms_field, fitwt_field, itnum_field, &
      saa_field, spcfit_field, spcobs_field, spcres_field, &
      sza_field, thgt_field, time_field, utc_field, vaa_field, &
      vza_field, xtr_field

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_write_radfit_output'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, iline, nXtrack, nblock, fpix, lpix
    !REAL    (KIND=r8), INTENT (IN), DIMENSION (1:n_fitvar_rad,1:nXtrack,0:nblock-1) :: &
    !     all_fitted_columns, all_fitted_errors, correlation_columns

    ! CCM
    ! REAL (KIND=r8), INTENT (IN), DIMENSION(n_comm_wvl,nxtrack_max,4,0:nt-1) :: omi_fitspc
    REAL (KIND=r8), INTENT (IN), &
      DIMENSION(n_rad_wvl_max,nxtrack_max,4,0:nblock-1) :: omi_fitspc ! JED

    ! Temporary array
    !REAL (KIND=r8),DIMENSION(n_rad_wvl_max,nxtrack_max,0:nt-1) :: tmp_fitspc
    REAL (KIND=r8),DIMENSION(n_rad_wvl,nXtrack,0:nblock-1) :: tmp_fitspc

    ! ------------------
    ! Modified variables gga
    ! ------------------
    ! The reason these appear as INTENT (INOUT) is that the rounding
    ! will modify the values. Hence INTENT (IN) produces a conflict.
    ! --------------------------------------------------------------
    REAL (KIND=r8), INTENT (INOUT), DIMENSION (1:n_fitvar_rad,1:nXtrack,0:nblock-1) :: &
      all_fitted_columns, all_fitted_errors, correlation_columns

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    !INTEGER   (KIND=i4) :: locerrstat, i, npix
    !CHARACTER (LEN=4)   :: molstr

    ! ---------------
    ! Local variables gga
    ! ---------------

    INTEGER   (KIND=i4)                      :: locerrstat, i, npix
    CHARACTER (LEN=4)                        :: molstr

    ! CCM loop indices
    INTEGER (KIND=i4) :: ii !,jj,kk
    TYPE (DataField_HE5), pointer :: ptr

    locerrstat = pge_errstat_ok

    write (0,*) "he5_write_radfit_output: n_rad_wvl_max=", n_rad_wvl_max, &
      " n_comm_wvl=", n_comm_wvl, " n_rad_wvl=", n_rad_wvl, &
      " nxtrack_max=", nxtrack_max, " nxtrack=", nxtrack

    ! -------------------------------------------------------
    ! Some fields have the PGE molecule name appended to them
    ! -------------------------------------------------------
    molstr = sao_molecule_names(pge_idx)(1:len(molstr))

    ! ---------------------------------------------------
    ! Write current data block fitting output to HE5 file
    ! ---------------------------------------------------

    ! -----------------------------------------------
    ! Number of cross-track pixels actually processed
    ! -----------------------------------------------
    npix = lpix - fpix + 1

    ! ----------------------------------------------------------------------------
    ! Geolocation Fileds: Latitiude, Longitude, Solar Zenith, Viewing Zenith, Time
    ! ----------------------------------------------------------------------------
    ! NOTE: The Aura Altitude and TIME fields are one-dimensional, and so require
    !       a different stride for writing. We write those ones first, then set
    !       the strides for the rest of the fields.
    ! ----------------------------------------------------------------------------
    CALL roundoff_2darr_r4 ( n_roff_dig, npix, nblock, omi_latitude (fpix:lpix,0:nblock-1) )
    CALL roundoff_2darr_r4 ( n_roff_dig, npix, nblock, omi_longitude(fpix:lpix,0:nblock-1) )
    CALL roundoff_2darr_r4 ( n_roff_dig, npix, nblock, omi_sazimuth (fpix:lpix,0:nblock-1) )
    CALL roundoff_2darr_r4 ( n_roff_dig, npix, nblock, omi_szenith  (fpix:lpix,0:nblock-1) )
    CALL roundoff_2darr_r4 ( n_roff_dig, npix, nblock, omi_vazimuth (fpix:lpix,0:nblock-1) )
    CALL roundoff_2darr_r4 ( n_roff_dig, npix, nblock, omi_vzenith  (fpix:lpix,0:nblock-1) )

    he5_start_2d  = (/ iline, 0 /) ;  he5_stride_2d = (/ 1, 0 /) ; he5_edge_2d = (/ nblock, 0 /)
    locerrstat = HE5_SWWRFLD ( pge_swath_id, auraalt_field,   he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_auraalt(0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, time_field,   he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_time(0:nblock-1) )

    ! --------
    ! UTC time
    ! --------
    he5_start_2d = (/ 0, iline /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nUTCdim, nblock /)
    locerrstat = HE5_SWWRFLD ( pge_swath_id, utc_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_time_utc(1:nUTCdim, 0:nblock-1) )

    he5_start_2d = (/ 0, iline /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, nblock /)
    locerrstat = HE5_SWWRFLD ( pge_swath_id, lat_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_latitude(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, lon_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_longitude(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, saa_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_sazimuth(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, sza_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_szenith(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, vaa_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_vazimuth(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, vza_field,    he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_vzenith(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, thgt_field,   he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_height(1:nXtrack,0:nblock-1) )

    ! ----------------------------------------------------------------------------------
    ! Row-anomaly related cross-track quality flags in two flavors: The original format,
    ! copied "as is" from the L1b file, and an "expanded" version that doesn't require
    ! a bit-conversion. gga
    ! ----------------------------------------------------------------------------------
    locerrstat = HE5_SWWRFLD ( pge_swath_id, xtr_field,   he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_xtrflg_l1b(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, extr_field,  he5_start_2d, he5_stride_2d, he5_edge_2d, &
                              omi_xtrflg(1:nXtrack,0:nblock-1) )

    ! -----------------------------
    ! Column amount and uncertainty
    ! -----------------------------
    ! ------------------------------------------------------------------
    ! NOTE: We do NOT perform any rounding here since the columns and
    !       their uncertainties will be converted from slant to vertical
    !       by the molecular AMF further down (for some PGEs anyway). A
    !       double rounding can only lessen the precision.
    ! ------------------------------------------------------------------
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(col_field)),     &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, omi_column_amount(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(dcol_field)),    &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, omi_column_uncert(1:nXtrack,0:nblock-1) )

    ! -----------------------------------------------------
    ! Fitting RMS and EXIT values from fitting.
    ! The latter shows if convergence was achieved and how.
    ! -----------------------------------------------------
    CALL roundoff_2darr_r8 ( n_roff_dig, npix, nblock, omi_fit_rms(fpix:lpix,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(fitrms_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, omi_fit_rms(1:nXtrack,0:nblock-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(fitcon_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, omi_fitconv_flag(1:nXtrack,0:nblock-1) )

    ! ---------------------------------------------------------
    ! PGE specific output:
    !
    ! * Individual O3 slant column      (OMSAO3)
    ! ---------------------------------------------------------
    SELECT CASE ( pge_idx )
    CASE ( pge_bro_idx )
      ! Nothing here
    CASE ( pge_o3_idx )
      ! ------------------------------------------------------------------
      ! Individual O3 slant columns (possibly multiple temperature x-secs)
      ! ------------------------------------------------------------------
      i = o3_t1_idx
      ptr => o3_prefit_he5fields%head
      do while (associated (ptr))
        CALL roundoff_2darr_r8 (n_roff_dig, npix, nblock, omi_o3_amount(i,fpix:lpix,0:nblock-1))
        locerrstat = HE5_SWWRFLD (pge_swath_id, TRIM(ptr%name), &
                                  he5_start_2d, he5_stride_2d, he5_edge_2d, &
                                  omi_o3_amount(i,1:nXtrack,0:nblock-1))
        ptr => ptr % next
        i = i + 1
      end do

      i = o3_t1_idx
      ptr => o3_prefit_uncert_he5fields%head
      do while (associated (ptr))
        CALL roundoff_2darr_r8 (n_roff_dig, npix, nblock, omi_o3_uncert(i,fpix:lpix,0:nblock-1))
        locerrstat = HE5_SWWRFLD (pge_swath_id, TRIM(ptr%name), &
                                  he5_start_2d, he5_stride_2d, he5_edge_2d, &
                                  omi_o3_uncert(i,1:nXtrack,0:nblock-1))
        ptr => ptr % next
        i = i + 1
      end do
      !DO i = o3_t1_idx, o3_t3_idx
      !  CALL roundoff_2darr_r8 ( n_roff_dig, npix, nblock, omi_o3_amount(i,fpix:lpix,0:nblock-1) )
      !  CALL roundoff_2darr_r8 ( n_roff_dig, npix, nblock, omi_o3_uncert(i,fpix:lpix,0:nblock-1) )
      !  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(o3_prefit_fields(i,1))), &
      !                            he5_start_2d, he5_stride_2d, he5_edge_2d, omi_o3_amount(i,1:nXtrack,0:nblock-1) )
      !  locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(o3_prefit_fields(i,2))), &
      !                            he5_start_2d, he5_stride_2d, he5_edge_2d, omi_o3_uncert(i,1:nXtrack,0:nblock-1) )
      !END DO
    END SELECT

    ! --------------------------------------------------
    ! Correlation Information (requires additional rank)
    ! --------------------------------------------------
    IF ( yn_diagnostic_run ) THEN

      he5_start_2d = (/ 0, iline /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, nblock /)

      !FIXME: IF( yn_output_diag( itnum_didx ) ) &
        locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(itnum_field)), &
                                  he5_start_2d, he5_stride_2d, he5_edge_2d, omi_itnum_flag(1:nXtrack,0:nblock-1) )

      DO i = 1, n_fitvar_rad
        CALL roundoff_2darr_r8 ( n_roff_dig, npix, nblock, correlation_columns(i,fpix:lpix,0:nblock-1) )
        CALL roundoff_2darr_r8 ( n_roff_dig, npix, nblock, all_fitted_columns (i,fpix:lpix,0:nblock-1) )
        CALL roundoff_2darr_r8 ( n_roff_dig, npix, nblock, all_fitted_errors  (i,fpix:lpix,0:nblock-1) )
      END DO

      he5_start_3d  = (/ 0,            0,      iline       /)
      he5_stride_3d = (/ 1,            1,      1           /)
      he5_edge_3d   = (/ n_fitvar_rad, nXtrack, nblock /)

      !FIXME IF( yn_output_diag( corr_didx ) ) &
      if (yn_diagnostic_run) &
        locerrstat = HE5_SWWRFLD ( pge_swath_id, corr_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  correlation_columns(1:n_fitvar_rad,1:nXtrack,0:nblock-1) )

      !IF( yn_output_diag( corrcol_didx ) ) &
      if (yn_diagnostic_run) &
        locerrstat = HE5_SWWRFLD ( pge_swath_id, corrcol_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  all_fitted_columns(1:n_fitvar_rad,1:nXtrack,0:nblock-1) )

      !IF( yn_output_diag( correrr_didx ) ) THEN &
      if (yn_diagnostic_run) &
        locerrstat = HE5_SWWRFLD ( pge_swath_id, correrr_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  all_fitted_errors(1:n_fitvar_rad,1:nXtrack,0:nblock-1) )

    END IF

    ! -------------------------------
    ! CCM Write Fit residuals to disk
    ! -------------------------------
    IF( yn_diagnostic_run ) THEN

      ! Write to he5
      he5_start_3d  = (/ 0,            0,      iline       /)
      he5_stride_3d = (/ 1,            1,      1           /)
      !he5_edge_3d   = (/ n_comm_wvl, nXtrack, nblock /)
      he5_edge_3d   = (/ n_rad_wvl, nXtrack, nblock /)

      ! (1) Model Spectrum
      !FIXME IF( yn_output_diag( spcfit_didx ) ) THEN
      if (.true.) then
        !DO ii=1,n_comm_wvl
        !  DO jj=1,nXtrack
        !    DO kk=0,nblock-1
        !      tmp_fitspc(ii,jj,kk) = omi_fitspc(ii,jj,1,kk)
        !    ENDDO
        !  ENDDO
        !ENDDO
        !tmp_fitspc (1:n_comm_wvl, 1:nXtrack, 0:nblock-1) = omi_fitspc (1:n_comm_wvl, 1:nXtrack, 1, 0:nblock-1)
        tmp_fitspc = omi_fitspc (1:n_rad_wvl, 1:nXtrack, 1, 0:nblock-1)
        locerrstat = HE5_SWWRFLD ( pge_swath_id,spcfit_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  tmp_fitspc(1:n_rad_wvl,1:nXtrack,0:nblock-1))
      ENDIF

      ! (2) Measured Spectrum
      ! FIXME: IF( yn_output_diag( spcobs_didx ) ) THEN
      if (.true.) then
        !DO ii=1,n_comm_wvl
        !  DO jj=1,nXtrack
        !    DO kk=0,nblock-1
        !      tmp_fitspc(ii,jj,kk) = omi_fitspc(ii,jj,2,kk)
        !    ENDDO
        !  ENDDO
        !ENDDO
        !tmp_fitspc (1:n_comm_wvl, 1:nXtrack, 0:nblock-1) = omi_fitspc (1:n_comm_wvl, 1:nXtrack, 2, 0:nblock-1)
        tmp_fitspc = omi_fitspc (1:n_rad_wvl, 1:nXtrack, 2, 0:nblock-1)
        locerrstat = HE5_SWWRFLD ( pge_swath_id,spcobs_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  tmp_fitspc(1:n_rad_wvl,1:nXtrack,0:nblock-1) )
      ENDIF

      ! (3) Measured Position
      !IF( yn_output_diag( posobs_didx ) ) THEN
      if (.true.) then
        !DO ii=1,n_comm_wvl
        DO ii=1,n_rad_wvl
          !DO jj=1,nXtrack
          !  DO kk=0,nblock-1
          !    tmp_fitspc(ii,jj,kk) = omi_fitspc(ii,jj,3,kk)
          !  ENDDO
          !ENDDO
          !tmp_fitspc (ii, 1:nXtrack, 0:nblock-1) = omi_fitspc (ii, 1:nXtrack, 3, 0:nblock-1)
          tmp_fitspc(ii,:,:) = omi_fitspc (ii, 1:nXtrack, 3, 0:nblock-1)
          CALL roundoff_2darr_r8 ( n_roff_dig, nXtrack, nblock, tmp_fitspc(ii,1:nXtrack,0:nblock-1) )
        ENDDO

        locerrstat = HE5_SWWRFLD ( pge_swath_id,posobs_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  tmp_fitspc(1:n_rad_wvl,1:nXtrack,0:nblock-1) )
      ENDIF

      ! (4) Fit Weights
      ! FIXME: IF( yn_output_diag( fitwt_didx ) ) THEN
      if (.true.) then
        !DO ii=1,n_comm_wvl
        !  DO jj=1,nXtrack
        !    DO kk=0,nblock-1
        !      tmp_fitspc(ii,jj,kk) = omi_fitspc(ii,jj,4,kk)
        !    ENDDO
        !  ENDDO
        !ENDDO
        !tmp_fitspc (1:n_comm_wvl, 1:nXtrack, 0:nblock-1) = omi_fitspc (1:n_comm_wvl, 1:nXtrack, 4, 0:nblock-1)
        tmp_fitspc = omi_fitspc (1:n_rad_wvl, 1:nXtrack, 4, 0:nblock-1)
        locerrstat = HE5_SWWRFLD ( pge_swath_id,fitwt_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  tmp_fitspc(1:n_rad_wvl,1:nXtrack,0:nblock-1) )
      ENDIF

      ! Residual Spectrum
      ! FIXME: IF( yn_output_diag( spcres_didx ) ) THEN
      if (.true.) then
        !DO ii=1,n_comm_wvl
        !  DO jj=1,nXtrack
        !    DO kk=0,nblock-1
        !      tmp_fitspc(ii,jj,kk) = omi_fitspc(ii,jj,2,kk) - omi_fitspc(ii,jj,1,kk)
        !    ENDDO
        !  ENDDO
        !ENDDO
        !tmp_fitspc (1:n_comm_wvl, 1:nXtrack, 0:nblock-1) = &
        !  omi_fitspc (1:n_comm_wvl, 1:nXtrack, 2, 0:nblock-1) &
        !  - omi_fitspc (1:n_comm_wvl, 1:nXtrack, 1, 0:nblock-1)

        tmp_fitspc = &
          omi_fitspc (1:n_rad_wvl, 1:nXtrack, 2, 0:nblock-1) &
          - omi_fitspc (1:n_rad_wvl, 1:nXtrack, 1, 0:nblock-1)

        locerrstat = HE5_SWWRFLD ( pge_swath_id,spcres_field, he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                  tmp_fitspc(1:n_rad_wvl,1:nXtrack,0:nblock-1) )
      ENDIF

    END IF
    ! CCM

    ! ------------------
    ! Check error status
    ! ------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
                      modulename, vb_lev_default, errstat )
    RETURN

  END SUBROUTINE he5_write_radfit_output

  SUBROUTINE he5_write_common_mode ( nXtrack, npts, errstat )

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    !USE OMSAO_indices_module,   ONLY: commcnt_didx, commspc_didx, &
    !  commwvl_didx, ccdpix_didx
    USE OMSAO_variables_module, ONLY: common_mode_spec
    USE OMSAO_omidata_module,   ONLY: n_roff_dig
    USE sao_pge_utils, ONLY: roundoff_2darr_r4, roundoff_1darr_r8, roundoff_2darr_r8
    USE datafields

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_write_common_mode'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nXtrack, npts

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                            :: locerrstat, j1, j2
    INTEGER   (KIND=i2), DIMENSION (nXtrack,2)     :: locccd
    INTEGER   (KIND=i4), DIMENSION (nXtrack)       :: loccnt
    REAL      (KIND=r4), DIMENSION (nXtrack, npts) :: locwvl
    REAL      (KIND=r8), DIMENSION (nXtrack, npts) :: locspc

    locerrstat = pge_errstat_ok

    ! ----------------------------------------------------------------------------
    ! Geolocation Fields: Latitiude, Longitude, Solar Zenith, Viewing Zenith, Time
    ! ----------------------------------------------------------------------------
    ! NOTE: The Aura Altitude and TIME fields are one-dimensional, and so require
    !       a different stride for writing. We write those ones first, then set
    !       the strides for the rest of the fields.
    ! ----------------------------------------------------------------------------
    j1                       =       1
    j2                       = nXtrack
    locccd(1:nXtrack,1:2)    =        common_mode_spec%CCDPixel    (j1:j2,1:2)
    loccnt(1:nXtrack)        =        common_mode_spec%RefSpecCount(j1:j2)
    locspc(1:nXtrack,1:npts) =        common_mode_spec%RefSpecData (j1:j2,1:npts)
    locwvl(1:nXtrack,1:npts) = REAL ( common_mode_spec%RefSpecWavs (j1:j2,1:npts), KIND=r4 )

    !DO j1 = 1, npts
    !   WRITE (22,'(0PF10.4, 1PE15.5)') locwvl(10,j1), locspc(10,j1)
    !END DO

    CALL roundoff_2darr_r4 ( n_roff_dig, nXtrack, npts, locwvl (1:nXtrack,1:npts) )
    CALL roundoff_2darr_r8 ( n_roff_dig, nXtrack, npts, locspc (1:nXtrack,1:npts) )

    ! ------------------------------------------
    ! Common Mode Spectrum - Count for Averaging
    ! ------------------------------------------
    ! FIXME IF( yn_output_diag(commcnt_didx) ) THEN
    if (yn_diagnostic_run) then
      he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 0 /) ; he5_edge_2d = (/ nXtrack, 0 /)
      locerrstat = HE5_SWWRFLD ( &
        pge_swath_id, commcnt_field, he5_start_2d, he5_stride_2d, he5_edge_2d, loccnt(1:nXtrack) )
    ENDIF

    ! ----------------------------------------------
    ! Common Mode Spectrum - Wavelengths and Spectra
    ! ----------------------------------------------
    he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, npts /)

    !FIXME IF( yn_output_diag(commwvl_didx) ) THEN
    if (yn_diagnostic_run) &
      locerrstat = HE5_SWWRFLD ( &
        pge_swath_id, commwvl_field, he5_start_2d, he5_stride_2d, he5_edge_2d, locwvl(1:nXtrack,1:npts) )

    !FIXME: IF( yn_output_diag(commspc_didx) ) THEN
    if (yn_diagnostic_run) &
      locerrstat = HE5_SWWRFLD ( &
        pge_swath_id, commspc_field, he5_start_2d, he5_stride_2d, he5_edge_2d, locspc(1:nXtrack,1:npts) )

    ! --------------------------
    ! CCD Pixel - First and Last
    ! --------------------------
    ! FIXME: IF ( yn_diagnostic_run .AND. yn_output_diag(ccdpix_didx)) THEN
    if (yn_diagnostic_run) then
      he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nXtrack, 2 /)
      locerrstat = HE5_SWWRFLD ( &
        pge_swath_id, ccdpix_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
        locccd(1:nXtrack,1:2) )
    END IF
    ! ------------------
    ! Check error status
    ! ------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
                      modulename, vb_lev_default, errstat )

    RETURN
  END SUBROUTINE he5_write_common_mode

  SUBROUTINE he5_write_omi_database ( database_he5, database_he5_wvl, nRefSpec,  npts, nXtrack, errstat )

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    !USE OMSAO_omidata_module,   ONLY: n_roff_dig
    USE OMSAO_variables_module, ONLY: refspecs_original
    !USE OMSAO_indices_module,   ONLY: &
    !  spdata_didx, spnrmf_didx, spdatw_didx
    USE datafields
    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=23), PARAMETER :: modulename = 'he5_write_omi_database'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN)                                   :: nXtrack, npts, nRefSpec
    REAL    (KIND=r8), INTENT (IN), DIMENSION(nRefSpec,npts,nXtrack) :: database_he5
    REAL    (KIND=r8), INTENT (IN), DIMENSION(npts,nXtrack)          :: database_he5_wvl

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)                 :: locerrstat,ii

    REAL (KIND=r8), DIMENSION(nRefSpec) :: tmp_normfactor

    ! =============================================================================
    ! he5_write_omi_database starts here
    ! =============================================================================

    ! Set local error status
    locerrstat = pge_errstat_ok

    ! Datablock
    he5_start_3d  = (/ 0,            0,       0 /)
    he5_stride_3d = (/ 1,            1,       1 /)
    he5_edge_3d   = (/ nRefSpec,  npts, nXtrack /)

    ! Write refspec database
    ! FIXME IF( yn_output_diag(spdata_didx) ) THEN
    if (yn_diagnostic_run) &
      locerrstat = HE5_SWWRFLD (pge_swath_id, spdata_field, &
                                he5_start_3d, he5_stride_3d, he5_edge_3d, &
                                database_he5 )

    ! Datablock
    he5_start_2d  = (/    0,       0 /)
    he5_stride_2d = (/    1,       1 /)
    he5_edge_2d   = (/ npts, nXtrack /)
    ! Write refspec database wavelength
    !FIXME IF( yn_output_diag(spdatw_didx) ) THEN
    if (yn_diagnostic_run) &
      locerrstat = HE5_SWWRFLD ( pge_swath_id, spdatw_field, &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, &
                                database_he5_wvl )
    DO ii=1,nRefSpec
      tmp_normfactor(ii) = refspecs_original(ii)%NormFactor
    END DO

    ! Datablock
    he5_start_2d  = (/ 0, 0 /)
    he5_stride_2d = (/ 1, 0 /)
    he5_edge_2d = (/ nRefSpec, 0 /)

    ! Write Normalisation factors
    ! FIXME IF( yn_output_diag( spnrmf_didx ) ) THEN
    if (yn_diagnostic_run) &
      locerrstat = HE5_SWWRFLD ( pge_swath_id, spnrmf_field, &
                                he5_start_2d, he5_stride_2d, he5_edge_2d,&
                                tmp_normfactor)

    ! Write Refspec Names
    !he5_start_2d  = (/ zerocl, zerocl /)
    !he5_stride_2d = (/  onecl, zerocl /)
    !he5_edge_2d   = (/ INT(nclenfit,KIND=C_LONG), zerocl /)
    !locerrstat = HE5_SWWRFLD ( &
    !        pge_swath_id, correlm_field,  he5_start_2d, he5_stride_2d, he5_edge_2d, &
    !        TRIM(ADJUSTL(refspec_strings)) )
    ! ------------------
    ! Check error status
    ! ------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
                      modulename, vb_lev_default, errstat )

    RETURN

  END SUBROUTINE he5_write_omi_database

  SUBROUTINE he5_write_fitting_statistics ( &
      pge_idx, max_col, nx, nt, saomqf, avg_col, avg_dcol, avg_rms, errstat )

    USE OMSAO_precision_module, ONLY: C_LONG, r8, i2
    USE OMSAO_indices_module,    ONLY: sao_molecule_names !, correlm_didx
    USE OMSAO_errstat_module
    USE OMSAO_parameters_module, ONLY: NLINES_MAX
    USE OMSAO_omidata_module,    ONLY: correlation_names_concat, nclenfit
    USE datafields, ONLY: avgcol_field, avgdcol_field, avgrms_field, &
      correlm_field, mainqa_field, maxcol_field

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=28), PARAMETER :: modulename = 'he5_write_fitting_statistics'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nx, nt
    REAL    (KIND=r8),                          INTENT (IN) :: max_col, avg_col, avg_dcol, avg_rms
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: saomqf

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=C_LONG), PARAMETER :: zerocl = 0, onecl = 1
    CHARACTER (LEN=4)        :: molstr
    INTEGER   (KIND=i4)      :: locerrstat, iline, nt_loop

    ! ---------------------------------------------------
    ! Write current data block fitting output to HE5 file
    ! ---------------------------------------------------

    locerrstat = pge_errstat_ok

    ! ----------------------------------------------------------
    ! Average fields have the PGE molecule name appended to them
    ! ----------------------------------------------------------
    molstr = sao_molecule_names(pge_idx)(1:len(molstr))

    ! ----------------------------------------
    ! Write the Main Data Quality Flag to file
    ! ----------------------------------------
    DO iline = 0, nt-1, NLINES_MAX

      ! --------------------------------------------------------
      ! Check if loop ends before n_times_loop max is exhausted.
      ! --------------------------------------------------------
      nt_loop = MIN ( NLINES_MAX, nt-iline )

      ! --------------------------------------------
      ! Set Start, Stride, and Edge of writing block
      ! --------------------------------------------
      he5_start_2d  = (/            zerocl,   INT(iline,KIND=C_LONG) /)
      he5_stride_2d = (/            onecl,                 onecl /)
      he5_edge_2d   = (/ INT(nx,KIND=C_LONG), INT(nt_loop,KIND=C_LONG) /)

      ! ---------------
      ! Write MDQA flag
      ! ---------------
      locerrstat = HE5_SWwrfld ( pge_swath_id, TRIM(ADJUSTL(mainqa_field)),           &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, saomqf(1:nx,iline:iline+nt_loop-1) )

    END DO

    ! ----------------------------------------------
    ! All quantities written below are single numbers
    ! ----------------------------------------------
    he5_start_2d  = (/ zerocl, zerocl /)
    he5_stride_2d = (/  onecl, zerocl /)
    he5_edge_2d   = (/  onecl, zerocl /)

    ! ----------------------------------------
    ! Maximum column amount for QA flag "good"
    ! ----------------------------------------
    locerrstat = HE5_SWwrfld ( &
      pge_swath_id,  TRIM(ADJUSTL(maxcol_field)), &
      he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), max_col  )

    ! ---------------------------------------------------
    ! Average column amount, uncertainty, and fitting RMS
    ! ---------------------------------------------------
    locerrstat = HE5_SWwrfld ( &
      pge_swath_id,  TRIM(ADJUSTL(avgcol_field)), &
      he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), avg_col  )
    locerrstat = HE5_SWwrfld ( &
      pge_swath_id, TRIM(ADJUSTL(avgdcol_field)), &
      he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), avg_dcol )
    locerrstat = HE5_SWwrfld ( &
      pge_swath_id,  TRIM(ADJUSTL(avgrms_field)), &
      he5_start_2d(1), he5_stride_2d(1), he5_edge_2d(1), avg_rms  )

    ! -------------------------------------------
    ! Elements included in the Correlation Output
    ! -------------------------------------------
    ! This was added later, after the routine had been named to include AVERAGE.
    ! It is a bit confusing, but no matter where we places it at this point, it
    ! will be. And this is one of the few places outside the scan line loop, so
    ! the field won't be written over and over again.
    ! ----------------------------------------------------------------------------
    ! To make matters worse, STRINGS cannot be written to the swath (at least on
    ! non-TLCF implementations of HDF-EOS5). We have thus adopted the work-around
    ! solution of converting everything to INTEGERs first and write those to file.
    ! ----------------------------------------------------------------------------
    ! FIXME. Need to test output field of the linked list here
    !IF (yn_diagnostic_run & .AND. yn_output_diag(correlm_didx) ) THEN
    IF (yn_diagnostic_run) THEN
      he5_start_2d  = (/ zerocl, zerocl /)
      he5_stride_2d = (/  onecl, zerocl /)
      he5_edge_2d   = (/ INT(nclenfit,KIND=C_LONG), zerocl /)
      locerrstat = HE5_SWWRFLD ( &
        pge_swath_id, correlm_field,  he5_start_2d, he5_stride_2d, he5_edge_2d, &
        TRIM(ADJUSTL(correlation_names_concat)) )
    END IF
    ! ------------------
    ! Check error status
    ! ------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
                      modulename, vb_lev_default, errstat )

    RETURN
  END SUBROUTINE he5_write_fitting_statistics

!UNUSED!   SUBROUTINE he5_write_amf ( &
!UNUSED!       pge_idx, nx, nt, saocol, saodco, amfmol, amfgeo, amfdiag, &
!UNUSED!       amfcfr, amfctp, errstat )
!UNUSED!
!UNUSED!     USE OMSAO_precision_module, ONLY: i2, i4, r8
!UNUSED!     USE OMSAO_errstat_module
!UNUSED!     USE OMSAO_omidata_module,   ONLY: n_roff_dig
!UNUSED!     USE OMSAO_indices_module,   ONLY: pge_hcho_idx, pge_gly_idx, pge_h2o_idx
!UNUSED!     USE sao_pge_utils, ONLY: roundoff_2darr_r4, roundoff_1darr_r8, roundoff_2darr_r8
!UNUSED!     USE datafields
!UNUSED!
!UNUSED!     IMPLICIT NONE
!UNUSED!
!UNUSED!     ! ------------------------------
!UNUSED!     ! Name of this module/subroutine
!UNUSED!     ! ------------------------------
!UNUSED!     !CHARACTER (LEN=13), PARAMETER :: modulename = 'he5_write_amf'
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Input variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nx, nt
!UNUSED!     REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: saocol, saodco
!UNUSED!     REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfmol, amfgeo
!UNUSED!     REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfcfr, amfctp
!UNUSED!     INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfdiag
!UNUSED!
!UNUSED!     ! -----------------
!UNUSED!     ! Modified variable
!UNUSED!     ! -----------------
!UNUSED!     INTEGER (KIND=i4), INTENT (INOUT) :: errstat
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Local variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4)                          :: locerrstat
!UNUSED!     REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1) :: amfloc
!UNUSED!     REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc
!UNUSED!
!UNUSED!     ! -------------------------------------------------------------
!UNUSED!     ! Air mass factor plus diagnostic.
!UNUSED!     ! -------------------------------------------------------------
!UNUSED!     ! As yet, only OMBRO and OMHCHO have true, non-geometric AMFs.
!UNUSED!     ! But we try to have symmetric data fields as much as possible,
!UNUSED!     ! hence the presence of the "molecule specific" AMF and its
!UNUSED!     ! diagnostic for all PGEs. Non-OMBRO and -OMHCHO PGEs carry a
!UNUSED!     ! geometric AMF here.
!UNUSED!     !
!UNUSED!     ! For completeness, the geometric AMF is added.
!UNUSED!     ! -------------------------------------------------------------
!UNUSED!
!UNUSED!     locerrstat = pge_errstat_ok
!UNUSED!
!UNUSED!     he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)
!UNUSED!
!UNUSED!     ! ----------------------------------------
!UNUSED!     ! (1) AMF diagnostic. No rounding required
!UNUSED!     ! ----------------------------------------
!UNUSED!     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfdiag_field)), &
!UNUSED!                               he5_start_2d, he5_stride_2d, he5_edge_2d, amfdiag(1:nx,0:nt-1) )
!UNUSED!     errstat = MAX ( errstat, locerrstat )
!UNUSED!
!UNUSED!     ! -----------------
!UNUSED!     ! (2) Geometric AMF
!UNUSED!     ! -----------------
!UNUSED!     amfloc = REAL ( amfgeo, KIND=r4 )
!UNUSED!     CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
!UNUSED!     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfgeo_field)), &
!UNUSED!                               he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
!UNUSED!     errstat = MAX ( errstat, locerrstat )
!UNUSED!
!UNUSED!     ! -----------------
!UNUSED!     ! (3) Molecular AMF
!UNUSED!     ! -----------------
!UNUSED!     amfloc = REAL ( amfmol, KIND=r4 )
!UNUSED!     CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
!UNUSED!     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfmol_field)), &
!UNUSED!                               he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
!UNUSED!     errstat = MAX ( errstat, locerrstat )
!UNUSED!
!UNUSED!     ! ----------------------------------------------------------
!UNUSED!     ! (4) OMHCHO, OMCHOCHO only: AMF cloud fraction and pressure
!UNUSED!     ! ----------------------------------------------------------
!UNUSED!     IF ( pge_idx == pge_hcho_idx .OR. pge_idx == pge_gly_idx .OR. pge_idx == pge_h2o_idx) THEN
!UNUSED!       amfloc = REAL ( amfcfr, KIND=r4 )
!UNUSED!       CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
!UNUSED!       locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfcfr_field)), &
!UNUSED!                                 he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
!UNUSED!       errstat = MAX ( errstat, locerrstat )
!UNUSED!
!UNUSED!       amfloc = REAL ( amfctp, KIND=r4 )
!UNUSED!       CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
!UNUSED!       locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfctp_field)), &
!UNUSED!                                 he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
!UNUSED!       errstat = MAX ( errstat, locerrstat )
!UNUSED!     END IF
!UNUSED!
!UNUSED!     ! -----------------------------------------------------------------------
!UNUSED!     ! (5) All PGEs: Output of columns and column uncertainties. For some PGEs
!UNUSED!     !     (e.g., OMBRO, OMHCHO, OMCHOCHO) those have been adjusted by the AMF,
!UNUSED!     !     but we have as yet to perform the rounding for any of them.
!UNUSED!     ! -----------------------------------------------------------------------
!UNUSED!     colloc = saocol
!UNUSED!     CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
!UNUSED!
!UNUSED!     ! Why are *col_field being written again?  See he5_write_radfit_output. --JED
!UNUSED!     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(col_field)), &
!UNUSED!                               he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
!UNUSED!     errstat = MAX ( errstat, locerrstat )
!UNUSED!
!UNUSED!     colloc = saodco
!UNUSED!     CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
!UNUSED!     locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(dcol_field)), &
!UNUSED!                               he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
!UNUSED!     errstat = MAX ( errstat, locerrstat )
!UNUSED!
!UNUSED!     RETURN
!UNUSED!   END SUBROUTINE he5_write_amf

  FUNCTION he5_set_field_attributes ( pge_idx ) RESULT ( he5stat )

    !----------------------------------------------------------------------
    ! This function sets to FillValue all those entries in the swath fields
    ! that haven't been initialized properly at this point.
    !
    ! It is called after all processing is complete.
    ! ---------------------------------------------------------------------
    !
    ! Input:
    !   pge_idx    - Index for current PGE
    !
    ! Return: he5stat = OMI_E_SUCCESS if it is
    !
    ! Variables passed through MODULE:
    !   pge_swath_file_id - id number for HE5 output file (required for closing it)
    !   pge_swath_id      - id number for swath (required for writing to swath)
    !
    !------------------------------------------------------------------------------

    USE OMSAO_indices_module,    ONLY: &
      pge_hcho_idx, pge_gly_idx, pge_o3_idx, sao_molecule_names, &
      pge_h2o_idx
    !USE OMSAO_parameters_module, ONLY: r8_missval
    USE datafields
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=24), PARAMETER :: modulename = 'he5_set_field_attributes'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4),      INTENT(IN) :: pge_idx

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: he5stat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4) :: locerrstat
    CHARACTER (LEN=4)   :: molstr

    he5stat    = pge_errstat_ok
    locerrstat = pge_errstat_ok

    ! ---------------------------------------------------
    ! Some fields have PGE molecule name appended to them
    ! ---------------------------------------------------
    molstr = sao_molecule_names (pge_idx)(1:len(molstr))

    ! --------------------------------------------------------------------------
    ! Set the Title, Units, and MissingValue for all Geolocation and Data Fields.
    ! This is a non-critical operation, so for now we don't care whether the
    ! error status might be "not O.K."
    !
    ! This is a lot of repeated code, so eventually we want to come up with
    ! something more elegant. But for now it will have to do.
    ! --------------------------------------------------------------------------

    ! ------------------
    ! Geolocation Fields
    ! ------------------
    call write_attributes (geo_he5fields, locerrstat)

    ! ---------------------------------------------------------------------
    ! Data fields associated with solar and radiance wavelength calibration
    ! ---------------------------------------------------------------------
    CALL write_attributes (sol_calfit_he5fields, locerrstat)

    CALL write_attributes (rad_calfit_he5fields, locerrstat)

    CALL write_attributes (rad_reffit_he5fields, locerrstat)

    ! ----------------
    ! Main Data Fields
    ! ----------------
    CALL write_attributes (comdata_he5fields, locerrstat)

    ! ----------------------
    ! Diagnostic Data Fields
    ! ----------------------
    IF ( yn_diagnostic_run ) THEN
      CALL write_attributes (diagnostic_he5fields, locerrstat)
    END IF

    ! ----------------------------------
    ! Special fields for individual PGEs
    ! ----------------------------------
    SELECT CASE ( pge_idx )
    CASE ( pge_hcho_idx )          ! OMHCHO
      CALL write_attributes (voc_he5fields, locerrstat)

      ! ------------------------------------------------------
      ! Reference Sector Correction 2 fields: Total Column and
      !   uncertainty but so far by limiting the number of loo
      !   ps only 1 is used Total Column gga
      ! ------------------------------------------------------
      ! FIXME: This will do all elements of the list despite the above
      ! comment.
      IF (yn_refseccor) THEN
        CALL write_attributes (rs_he5fields, locerrstat)
      ENDIF

    CASE ( pge_gly_idx) ! OMCHOCHO
        CALL write_attributes (voc_he5fields, locerrstat)

    CASE ( pge_h2o_idx) ! OMH2O
      CALL write_attributes (voc_he5fields, locerrstat)

    CASE ( pge_o3_idx   )          ! OMSAO3
      ! FIXME: The "Special case" in write_an_attribute should be triggered here.
      ! The original code used "rstemp_attr" in these calls/
      CALL write_attributes (o3_prefit_he5fields, locerrstat)
      CALL write_attributes (o3_prefit_uncert_he5fields, locerrstat)

    CASE DEFAULT
      ! Nothing to be done here
    END SELECT

    ! --------------------------------------------------
    ! Scattering weights, gas profile, averaging kernels
    ! and albedo attributes. gga
    ! --------------------------------------------------
    IF (yn_scat_weights) THEN
      CALL write_attributes (sw_he5fields, locerrstat)
    ENDIF

    ! -----------------------------------------------------------------------
    ! Assign the error status returned from this routine. If we have
    ! reached here, then everything has gone well and we can exit with "O.K."
    ! -----------------------------------------------------------------------
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5SWRLATTR, &
                      modulename, vb_lev_default, he5stat )

    RETURN
  END FUNCTION he5_set_field_attributes

  FUNCTION he5_write_global_attributes ( ) RESULT ( he5stat )

    !----------------------------------------------------------------------
    ! This function writes HDF Global Attributes defined by the PGE.
    ! It uses the HE5_EHwrtglatt function, which writes the attribute,
    ! and defines it beforehand should it not exists.
    !
    ! It is called after all processing is complete.
    ! ---------------------------------------------------------------------
    !
    ! Input:
    !   none
    !
    ! Return: he5stat = OMI_S_SUCCESS if it is
    !
    ! Variables passed through MODULE:
    !   pge_swath_file_id - id number for HE5 output file (required for closing it)
    !   pge_swath_id      - id number for swath (required for writing to swath)
    !
    !------------------------------------------------------------------------------

    USE OMSAO_precision_module, ONLY: C_LONG, r4, r8
    USE metadata_tools, ONLY: n_mdata_dbl, n_mdata_str, mdata_string_fields, &
      mdata_string_values, mdata_double_fields, mdata_double_values
    USE OMSAO_errstat_module
    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, N_FIT_WINWAV
    USE OMSAO_indices_module,    ONLY: &
      n_config_luns, yn_config_lun_autocopy, config_lun_strings, config_lun_values
    USE OMSAO_variables_module,  ONLY: ctrl_fit_winwav_lim, ctrl_fit_winexc_lim

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=27), PARAMETER :: modulename = 'he5_write_global_attributes'

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: he5stat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=C_LONG), PARAMETER ::  onecl = 1
    INTEGER   (KIND=i4)      :: locerr, i, imd
    INTEGER   (KIND=C_LONG)     :: nlen48
    CHARACTER (LEN=MAX_STR_LEN) :: parname
    REAL      (KIND=r4), DIMENSION (N_FIT_WINWAV+2) :: fitwinlim

    locerr  = pge_errstat_ok
    he5stat = pge_errstat_ok

    ! -------------------------------------------------------------------
    ! First Global Attributes required by the Aura File Format Guidelines
    ! -------------------------------------------------------------------
    ! The first loop auto-copies some strings to the L2 output file
    ! -------------------------------------------------------------
    DO i = 1, n_config_luns
      IF ( yn_config_lun_autocopy(i) ) THEN
        parname = ADJUSTL(config_lun_strings(i))
        nlen48 = INT( LEN_TRIM(ADJUSTL(config_lun_values(i))), KIND=C_LONG )
        locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(parname), &
                                HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(config_lun_values(i))) )
        CALL error_check ( &
          locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
          modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
      END IF
    END DO

    parname = "InputVersions"
    nlen48 = INT ( LEN_TRIM(ADJUSTL(input_versions)), KIND=C_LONG )
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(input_versions)) )

    parname = "OrbitData"
    nlen48 = INT ( LEN_TRIM(ADJUSTL(l1b_orbitdata)), KIND=C_LONG )
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(l1b_orbitdata))           )

    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
    parname = "GranuleDay"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_INT, onecl, granule_day )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
    parname = "GranuleMonth"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_INT, onecl, granule_month )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )
    parname = "GranuleYear"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_INT, onecl, granule_year )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    ! -------------------------------------
    ! Finally the TAI time of 0z of Granule
    ! -------------------------------------
    parname = tai_attr
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), HE5T_NATIVE_DOUBLE, onecl, TAI93At0zOfGranule )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    ! ------------------------------------------------------------------------
    ! Now the Global Attributes that are set by the PGE: PGE statistics mostly
    ! ------------------------------------------------------------------------
    parname = "FittingWindowLimits"
    fitwinlim(1:N_FIT_WINWAV) = REAL( ctrl_fit_winwav_lim(1:N_FIT_WINWAV), KIND=r4 )
    fitwinlim(5:6)            = REAL( ctrl_fit_winexc_lim(1:2),            KIND=r4 )
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_FLOAT, INT( N_FIT_WINWAV+2, KIND=r8), fitwinlim(1:N_FIT_WINWAV+2) )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfCrossTrackPixels"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofCrossTrackPixels )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfScanLines"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofScanLines )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfInputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofInputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfGoodInputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofGoodInputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfGoodOutputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofGoodOutputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfSuspectOutputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofSuspectOutputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfBadOutputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofBadOutputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfConvergedSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofConvergedSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfFailedConvergenceSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofFailedConvergenceSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfExceededIterationsSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofExceededIterationsSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "NumberOfOutOfBoundsSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_INT, onecl, NrofOutofBoundsSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "PercentGoodOutputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_FLOAT, onecl, PercentGoodOutputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "PercentBadOutputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_FLOAT, onecl, PercentBadOutputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    parname = "PercentSuspectOutputSamples"
    locerr = HE5_EHwrglatt ( pge_swath_file_id, TRIM(ADJUSTL(parname)), &
                            HE5T_NATIVE_FLOAT, onecl, PercentSuspectOutputSamples )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5EHWRGLATT, &
      modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, he5stat )

    ! ------------------------------------------------------------------
    ! And here we have the ArchiveMetadata, which have to be written as
    ! Global Attributes. This section somewhat mirrors the one in
    ! SET_L2_METADATA, picking up the missing "arc" pieces. The Metadata
    ! section was written before it was known that the ArchivedMetadata
    ! could not be written as such.
    ! ------------------------------------------------------------------
    ! -----------------------------
    ! Write STRING MetaData to file
    ! -----------------------------
    DO imd = 1, n_mdata_str
      IF ( TRIM(ADJUSTL(mdata_string_fields(2,imd))) == "arc" ) THEN
        nlen48 = INT ( LEN_TRIM(ADJUSTL(mdata_string_values(imd))), KIND=C_LONG )
        locerr = HE5_EHwrglatt ( &
          pge_swath_file_id, TRIM(ADJUSTL(mdata_string_fields(1,imd))),  &
          HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(mdata_string_values(imd))) )
        CALL error_check ( &
          locerr, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2ARC, &
          modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,imd))), &
          vb_lev_default, he5stat )
      END IF
    END DO
    ! -------------------
    ! Set DOUBLE MetaData
    ! -------------------
    DO imd = 1, n_mdata_dbl
      IF ( TRIM(ADJUSTL(mdata_double_fields(2,imd))) == "arc" ) THEN
        locerr = HE5_EHwrglatt ( &
          pge_swath_file_id, TRIM(ADJUSTL(mdata_double_fields(1,imd))),  &
          HE5T_NATIVE_DOUBLE, onecl, mdata_double_values(imd) )
        CALL error_check ( &
          locerr, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2ARC, &
          modulename//f_sep//TRIM(ADJUSTL(mdata_double_fields(1,imd))), &
          vb_lev_default, he5stat )
      END IF
    END DO

    RETURN
  END FUNCTION he5_write_global_attributes

  FUNCTION he5_close_output_file ( pge_idx ) RESULT ( he5stat )

    !------------------------------------------------------------------------------
    ! This function detatches from the HE5 swath and closes the HE5 output file.
    !
    ! Input: NONE
    !
    !------------------------------------------------------------------------------

    USE OMSAO_indices_module,    ONLY: &
      pge_bro_idx, pge_hcho_idx, pge_oclo_idx, pge_gly_idx, &
      voc_isccp_idx, voc_isccp_idx, n_voc_amf_luns, n_voc_amf_luns
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------------------------------
    ! Name of this module/subroutine/function
    ! ---------------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'he5_close_output_file'

    ! --------------
    ! Input variable
    ! --------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: he5stat

    ! --------------
    ! Local variable
    ! --------------
    INTEGER (KIND=i4) :: locerr, i

    he5stat = pge_errstat_ok
    ! -----------------------------------------------
    ! Detach from HE5 swath and close HE5 output file
    ! -----------------------------------------------
    locerr = HE5_SWDETACH ( pge_swath_id )
    locerr = HE5_SWCLOSE  ( pge_swath_file_id )
    CALL error_check ( &
      locerr, HE5_STAT_OK, pge_errstat_warning, OMSAO_W_HE5SWCLOSE, &
      modulename, vb_lev_default, he5stat )

    ! -------------------------------------------------------------
    ! If we are running as OMHCHO and OMOCLO, we also need to close
    ! a few auxiliary files.
    ! -------------------------------------------------------------
    SELECT CASE ( pge_idx )
    CASE ( pge_bro_idx )
    CASE ( pge_hcho_idx )
      ! --------------------------------------------------------------------
      ! AMF table files; make sure we skip the ISCCP swath, because that one
      ! gets closed immediately after all the data are read from it.
      ! --------------------------------------------------------------------
      DO i = 1, n_voc_amf_luns
        IF ( i                                 /=  voc_isccp_idx .AND. &
            TRIM(ADJUSTL(amf_swath_names(i))) /= 'undefined'       .AND. &
            amf_swath_ids     (i)             /= -1                .AND. &
            amf_swath_file_ids(i)             /= -1                       ) THEN
          locerr = HE5_SWDETACH ( amf_swath_ids(i) )
          locerr = HE5_SWCLOSE  ( amf_swath_file_ids(i) )
        END IF
      END DO
      ! ------------------
      ! O3 and BrO prefits
      ! ------------------
      IF ( TRIM(ADJUSTL(o3fit_swath_name)) /= 'undefined' .AND. &
          o3fit_swath_id                  /= -1          .AND. &
          o3fit_swath_file_id             /= -1              ) THEN
        locerr = HE5_SWDETACH ( o3fit_swath_id )
        locerr = HE5_SWCLOSE  ( o3fit_swath_file_id )
      END IF
      IF ( TRIM(ADJUSTL(brofit_swath_name)) /= 'undefined' .AND. &
          brofit_swath_id                  /= -1          .AND. &
          brofit_swath_file_id             /= -1              ) THEN
        locerr = HE5_SWDETACH ( brofit_swath_id )
        locerr = HE5_SWCLOSE  ( brofit_swath_file_id )
      END IF
    CASE ( pge_gly_idx )
      ! --------------------------------------------------------------------
      ! AMF table files; make sure we skip the ISCCP swath, because that one
      ! gets closed immediately after all the data are read from it.
      ! --------------------------------------------------------------------
      DO i = 1, n_voc_amf_luns
        IF ( i                                 /=  voc_isccp_idx .AND. &
            TRIM(ADJUSTL(amf_swath_names(i))) /= 'undefined'       .AND. &
            amf_swath_ids     (i)             /= -1                .AND. &
            amf_swath_file_ids(i)             /= -1                       ) THEN
          locerr = HE5_SWDETACH ( amf_swath_ids(i) )
          locerr = HE5_SWCLOSE  ( amf_swath_file_ids(i) )
        END IF
      END DO
      ! --------------------
      ! Liquid Water prefits
      ! --------------------
      IF ( TRIM(ADJUSTL(lqh2ofit_swath_name)) /= 'undefined' .AND. &
          lqh2ofit_swath_id                  /= -1          .AND. &
          lqh2ofit_swath_file_id             /= -1              ) THEN
        locerr = HE5_SWDETACH ( lqh2ofit_swath_id )
        locerr = HE5_SWCLOSE  ( lqh2ofit_swath_file_id )
      END IF
    CASE ( pge_oclo_idx )
    END SELECT

    RETURN
  END FUNCTION he5_close_output_file

  SUBROUTINE he5_set_fill_value ( data_field, errstat )

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    USE OMSAO_parameters_module,     ONLY: str_missval
    USE datafields, ONLY: DataField_HE5

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    TYPE (DataField_HE5), INTENT (IN) :: data_field

    ! ----------------
    ! Output Variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    INTEGER (KIND=i2) :: locfil_i2
    INTEGER (KIND=i4) :: locfil_i4
    REAL    (KIND=r4) :: locfil_r4

    locerrstat = pge_errstat_ok

    ! MissingValue and ValidRange depend on data type
    SELECT CASE ( data_field%HE5_DataType )
    CASE ( HE5T_NATIVE_DOUBLE )
      errstat = HE5_SWsetfill (                   &
        data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
        data_field%HE5_DataType, data_field%FillValue )
    CASE ( HE5T_NATIVE_FLOAT )
      locfil_r4 = REAL (data_field%FillValue, KIND=r4)
      errstat = HE5_SWsetfill (                   &
        data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
        data_field%HE5_DataType, locfil_r4 )
    CASE ( HE5T_NATIVE_INT )
      locfil_i4 = INT (data_field%FillValue, KIND=i4)
      errstat = HE5_SWsetfill (                   &
        data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
        data_field%HE5_DataType, locfil_i4 )
    CASE ( HE5T_NATIVE_INT16 )
      locfil_i2 = INT (data_field%FillValue, KIND=i2)
      errstat = HE5_SWsetfill (                   &
        data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
        data_field%HE5_DataType, locfil_i2 )
    CASE ( HE5T_NATIVE_CHAR )
      errstat = HE5_SWsetfill (                   &
        data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
        data_field%HE5_DataType, str_missval )
    CASE DEFAULT
      ! We should never reach here
    END SELECT

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE he5_set_fill_value

  SUBROUTINE write_an_attribute ( addstr, data_field, errstat )

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    USE OMSAO_parameters_module,     ONLY: str_missval
    USE datafields, ONLY: DataField_HE5

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    CHARACTER (LEN=*),    INTENT (IN) :: addstr
    TYPE (DataField_HE5), INTENT (IN) :: data_field

    ! ----------------
    ! Output Variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER (KIND=C_LONG), PARAMETER :: onecl = 1 , twocl = 2
    INTEGER (KIND=i4) :: locerrstat
    INTEGER (KIND=i2)                :: locmis_i2
    INTEGER (KIND=i4)                :: locmis_i4
    REAL    (KIND=r4)                :: locmis_r4
    INTEGER (KIND=i2), DIMENSION (2) :: locrov_i2
    INTEGER (KIND=i4), DIMENSION (2) :: locrov_i4
    INTEGER (KIND=C_LONG)               :: nlen48
    REAL    (KIND=r4), DIMENSION (2) :: locrov_r4

    locerrstat = pge_errstat_ok

    ! Title
    nlen48 = data_field%LenTitle
    locerrstat = HE5_SWwrlattr ( &
      data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), title_attr, &
      HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(data_field%Title)) )
    ! Units
    nlen48 = data_field%LenUnits
    locerrstat = HE5_SWwrlattr ( &
      data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), units_attr, &
      HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(data_field%Units)) )
    ! UniqueFieldDefinition
    nlen48 = data_field%LenUniqueFD
    locerrstat = HE5_SWwrlattr ( &
      data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), ufd_attr, &
      HE5T_NATIVE_CHAR, nlen48, TRIM(ADJUSTL(data_field%UniqueFD)) )
    ! Offset and ScaleFactor
    ! (suppose to be always DoublePrecsion to avoid rounding errors and improve
    !  compressibility)
    locerrstat = HE5_SWwrlattr ( &
      data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), offset_attr, &
      HE5T_NATIVE_DOUBLE, onecl, data_field%Offset )
    locerrstat = HE5_SWwrlattr ( &
      data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), scafac_attr, &
      HE5T_NATIVE_DOUBLE, onecl, data_field%ScaleFactor )

    ! MissingValue and ValidRange depend on data type
    SELECT CASE ( data_field%HE5_DataType )
    CASE ( HE5T_NATIVE_DOUBLE )
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID,  &
                                  TRIM(ADJUSTL(data_field%Name)), missval_attr, &
                                  data_field%HE5_DataType, onecl, data_field%MissingValue )
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), valids_attr, &
                                  data_field%HE5_DataType, twocl, data_field%ValidRange )
    CASE ( HE5T_NATIVE_FLOAT )
      locmis_r4 = REAL (data_field%MissingValue, KIND=r4)
      locrov_r4 = REAL (data_field%ValidRange,   KIND=r4)
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), missval_attr, data_field%HE5_DataType, onecl, locmis_r4 )
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), valids_attr,  data_field%HE5_DataType, twocl, locrov_r4 )
    CASE ( HE5T_NATIVE_INT )
      locmis_i4 = INT (data_field%MissingValue, KIND=i4)
      locrov_i4 = INT (data_field%ValidRange,   KIND=i4)
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), missval_attr, data_field%HE5_DataType, onecl, locmis_i4 )
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), valids_attr,  data_field%HE5_DataType, twocl, locrov_i4 )
    CASE ( HE5T_NATIVE_INT16 )
      locmis_i2 = INT (data_field%MissingValue, KIND=i2)
      locrov_i2 = INT (data_field%ValidRange,   KIND=i2)
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), missval_attr, data_field%HE5_DataType, onecl, locmis_i2 )
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), valids_attr,  data_field%HE5_DataType, twocl, locrov_i2 )
    CASE ( HE5T_NATIVE_CHAR )
      nlen48 = LEN(str_missval)
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, &
                                  TRIM(ADJUSTL(data_field%Name)), missval_attr, &
                                  HE5T_NATIVE_CHAR, nlen48, str_missval )
    CASE DEFAULT
      ! We should never reach here
    END SELECT

    ! -------------
    ! Special cases
    ! -------------
    ! FIXME!! This special case is not triggered.
    IF ( TRIM(ADJUSTL(addstr)) == rstemp_attr ) THEN
      locerrstat = HE5_SWwrlattr ( data_field%Swath_ID, TRIM(ADJUSTL(data_field%Name)), &
                                  rstemp_attr, data_field%HE5_DataType, onecl, data_field%SpecTemp )
    END IF

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE write_an_attribute

  subroutine write_attributes (list, errstat)

    use datafields, ONLY: Datafield_List_Type, DataField_HE5
    implicit none
    type (Datafield_List_Type), intent(in) :: list
    integer (kind=i4), intent (inout) :: errstat
    ! local
    type (DataField_HE5), pointer :: ptr
    integer (kind=i4) :: locerr

    ptr => list % head
    do while (associated(ptr))
      locerr = 0
      CALL write_an_attribute ("", ptr, locerr)
      if (locerr /= 0) then
        write (*,*) "write_an_attribute err=", locerr, "name=", ptr%name
        errstat = locerr
      end if
      ptr => ptr % next
    end do
  end subroutine

  FUNCTION he5_write_swath_attributes ( pge_idx ) RESULT ( he5stat )

    !----------------------------------------------------------------------
    ! This function writes HE5 Swath Attributes defined by the PGE.
    ! It uses the HE5_SWwrtatt function, which writes the attribute,
    ! and defines it in case it does not exists.
    !
    ! It is called after all processing is complete.
    ! ---------------------------------------------------------------------
    !
    ! Input:  pge_idx
    !
    ! Return: he5stat = OMI_E_SUCCESS if it is
    !
    ! Variables passed through MODULE:
    !   pge_swath_id      - id number for swath (required for writing to swath)
    !
    !------------------------------------------------------------------------------

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    USE OMSAO_omidata_module, ONLY: EarthSunDistance
    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=26), PARAMETER :: modulename = 'he5_write_swath_attributes'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: he5stat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat
    INTEGER (KIND=C_LONG) :: count48

    ! ----------------------------------
    ! Write swath attributes to HE5 file
    ! ------------------------------------------------------------------
    ! So far there are only two: VerticalCoordinate and EarthSunDistance
    ! ------------------------------------------------------------------

    he5stat    = pge_errstat_ok
    locerrstat = pge_errstat_ok

    count48 = INT(LEN_TRIM(ADJUSTL(vertical_coordinate(pge_idx))), KIND=C_LONG)
    locerrstat = he5_swwrattr ( &
      pge_swath_id, vcoordinate_field, HE5T_NATIVE_CHAR, count48, &
      TRIM(ADJUSTL(vertical_coordinate(pge_idx))) )
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRATTR, &
                      modulename, vb_lev_default, he5stat )

    count48 = 1
    locerrstat = he5_swwrattr ( &
      pge_swath_id, 'EarthSunDistance', HE5T_NATIVE_FLOAT, count48, EarthSunDistance )
    CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRATTR, &
                      modulename, vb_lev_default, he5stat )

    RETURN
  END FUNCTION he5_write_swath_attributes

  SUBROUTINE he5_check_for_compressibility ( &
      nTimes, nXtrack, nSwLevels, field_dim, yn_compress_field, n_chunk_dim, chunk_dim )

    USE OMSAO_precision_module, ONLY: C_LONG
    USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
    USE OMSAO_parameters_module, ONLY: nUTCdim, NLINES_MAX, NWAVEL_MAX
    USE OMSAO_errstat_module
    USE OMSAO_variables_module, ONLY: n_fitvar_rad
    USE OMSAO_omidata_module,   ONLY: &
      nclenfit, n_field_maxdim, n_comm_wvl

    ! -----------------------------------------------------------------
    ! The purpose of this subroutine is to check whether compression
    ! should be enabled for a data field. The catch with compression
    ! is that
    !
    ! (1) it only works for 2- and higher dimensional arrays
    ! (2) it requires "chunking", i.e., the definition of the chunk
    !     of data that are being written at one time.
    !
    ! The second point causes some headache: For any fields that
    ! include the "nTimes" dimension,  we are limited by the number
    ! of swath lines that are read from the L1b file (usually 100),
    ! and we have either 2 or 3 dimensions; other fields may be
    ! written in one go. This defines the particular way of chunking.
    !
    ! The distinction can be made by checking the dimensions of the
    ! data fields and define chunking based on the strings that go
    ! into the definition of the data field.
    !
    ! The dimensions for which compression can be enabled are:
    !
    !    "nXtrack,2"
    !    "nXtrack,nCommonWavl"
    !    "nXtrack,nTimes"
    !    "nXtrack+1,nTimes+1"
    !    "nWavCalPars,nXtrack"
    !    "nUTCdim,nTimes"
    !    "nFitElements,nXtrack,nTimes"
    !    "nXtrack,nTimes,nSwLevels"
    !
    ! NO compression is possible for fields of dimension:
    !
    !    "1"
    !    "nTimes"
    !    "nXtrack"
    !    "nCharLenFitElements"
    !
    ! -----------------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4), INTENT (IN) :: nTimes, nXtrack, nSwLevels !gga
    CHARACTER (LEN=*),   INTENT (IN) :: field_dim

    ! ----------------
    ! Output variables
    ! ----------------
    LOGICAL,                                       INTENT (OUT) :: yn_compress_field
    INTEGER (KIND=i4),                             INTENT (OUT) :: n_chunk_dim
    INTEGER (KIND=C_LONG), DIMENSION (n_field_maxdim), INTENT (OUT) :: chunk_dim

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    yn_compress_field = .FALSE.
    n_chunk_dim       = 0
    chunk_dim(1:3)    = 0

    ! -----------------------------------------------
    ! Non-compression cases: Set only chunking
    ! -----------------------------------------------
    IF ( field_dim == '1' ) THEN
      n_chunk_dim  = 1
      chunk_dim(1) = 1
    END IF
    IF ( field_dim == '2' ) THEN
      n_chunk_dim  = 1
      chunk_dim(1) = 2
    END IF
    IF ( field_dim == '4' ) THEN
      n_chunk_dim  = 1
      chunk_dim(1) = 4
    END IF
    IF ( field_dim == "nTimes" ) THEN
      n_chunk_dim  = 1
      chunk_dim(1) = INT ( nTimes, KIND=C_LONG )
    END IF

    IF ( field_dim == "nXtrack" ) THEN
      n_chunk_dim  = 1
      chunk_dim(1) = INT ( nXtrack, KIND=C_LONG )
    END IF
    IF ( field_dim == "nCharLenFitElements" ) THEN
      n_chunk_dim  = 1
      chunk_dim(1) = INT ( nclenfit, KIND=C_LONG )
    END IF

    ! -------------------------------------------------------------------------
    ! Now the compressible cases, one by one. Since exact matches are required,
    ! there is no harm if we miss to RETURN on a non-compressible field.
    ! -------------------------------------------------------------------------

    IF ( field_dim == "nXtrack,2" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 2
      chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack,KIND=C_LONG), INT(2,KIND=C_LONG) /)
    END IF
    IF ( field_dim == "nXtrack,nCommonWavl" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 2
      chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack,KIND=C_LONG), INT(n_comm_wvl,KIND=C_LONG) /)
    END IF
    IF ( field_dim == "nXtrack,nTimes" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 2
      chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack,KIND=C_LONG), INT(NLINES_MAX,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nXtrack+1,nTimes+1" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 2
      chunk_dim(1:n_chunk_dim) = (/ INT(nXtrack+1,KIND=C_LONG), INT(NLINES_MAX,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nWavCalPars,nXtrack" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 2
      chunk_dim(1:n_chunk_dim) = (/ INT(max_calfit_idx,KIND=C_LONG), INT(nXtrack,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nUTCdim,nTimes" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 2
      chunk_dim(1:n_chunk_dim) = (/ INT(nUTCdim,KIND=C_LONG), INT(nTimes,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nFitElements,nXtrack,nTimes" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 3
      chunk_dim(1:n_chunk_dim) = &
        (/ INT(n_fitvar_rad,KIND=C_LONG), INT(nXtrack,KIND=C_LONG), INT(NLINES_MAX,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nXtrack,nTimes,nLevels" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 3
      chunk_dim(1:n_chunk_dim) = &
        (/ INT(nXtrack,KIND=C_LONG), INT(NLINES_MAX,KIND=C_LONG), INT(nSwLevels,KIND=C_LONG) /)
    END IF

    ! CCM New fields
    IF ( field_dim == "nCommonWavl,nXtrack,nTimes" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 3
      chunk_dim(1:n_chunk_dim) = &
        (/ INT(n_comm_wvl,KIND=C_LONG), INT(nXtrack,KIND=C_LONG), INT(NLINES_MAX,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nRfSpec,nwavel_max,nXtrack" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 3
      chunk_dim(1:n_chunk_dim) = &
        (/ INT(max_rs_idx,KIND=C_LONG), INT(nwavel_max,KIND=C_LONG), INT(nXtrack,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nwavel_max,nXtrack" ) THEN
      yn_compress_field = .TRUE.
      n_chunk_dim = 2
      chunk_dim(1:n_chunk_dim) = &
        (/ INT(nwavel_max,KIND=C_LONG), INT(nXtrack,KIND=C_LONG) /)
    END IF

    IF ( field_dim == "nRfSpec" ) THEN
      yn_compress_field = .FALSE.
      n_chunk_dim = 1
      chunk_dim(1:n_chunk_dim) = &
        (/ INT(max_rs_idx,KIND=C_LONG) /)
    END IF

    RETURN
  END SUBROUTINE he5_check_for_compressibility

!UNUSED!   FUNCTION he5_open_readwrite ( file_name, swath_name ) RESULT ( he5stat )
!UNUSED!
!UNUSED!     !------------------------------------------------------------------------------
!UNUSED!     ! This function initializes the HE5 output swath.
!UNUSED!     !
!UNUSED!     ! Input:
!UNUSED!     !   file_name  - Name of HE5 output file
!UNUSED!     !   swath_name - Name of swath to be created
!UNUSED!     !
!UNUSED!     ! Return: he5stat
!UNUSED!     !
!UNUSED!     ! Variables passed through MODULE:
!UNUSED!     !   pge_swath_file_id - id number for HE5 output file (required for closing it)
!UNUSED!     !   pge_swath_id      - id number for swath (required for writing to swath)
!UNUSED!     !
!UNUSED!     !------------------------------------------------------------------------------
!UNUSED!
!UNUSED!     USE OMSAO_errstat_module
!UNUSED!
!UNUSED!     IMPLICIT NONE
!UNUSED!
!UNUSED!     ! ------------------------------
!UNUSED!     ! Name of this module/subroutine
!UNUSED!     ! ------------------------------
!UNUSED!     CHARACTER (LEN=18), PARAMETER :: modulename = 'he5_open_readwrite'
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Input variables
!UNUSED!     ! ---------------
!UNUSED!     CHARACTER (LEN=MAX_STR_LEN), INTENT(IN) :: file_name, swath_name
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Result variable
!UNUSED!     ! ---------------
!UNUSED!     INTEGER (KIND=i4) :: he5stat
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Local variables
!UNUSED!     ! ---------------
!UNUSED!     INTEGER   (KIND=i4) :: errstat
!UNUSED!
!UNUSED!     he5stat = pge_errstat_ok
!UNUSED!     errstat = pge_errstat_ok
!UNUSED!
!UNUSED!     ! ---------------------------------------------------------------
!UNUSED!     ! Open HE5 output file and check PGE_SWATH_FILE_ID ( -1 if error)
!UNUSED!     ! ---------------------------------------------------------------
!UNUSED!     pge_swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(file_name)), he5f_acc_rdwr )
!UNUSED!     IF ( pge_swath_file_id == he5_stat_fail ) THEN
!UNUSED!       CALL error_check ( &
!UNUSED!         0, 1, pge_errstat_fatal, OMSAO_F_HE5SWOPEN, modulename, vb_lev_default, he5stat )
!UNUSED!       RETURN
!UNUSED!     END IF
!UNUSED!
!UNUSED!     ! ---------------
!UNUSED!     ! Attach to swath
!UNUSED!     ! ---------------
!UNUSED!     pge_swath_id  = HE5_SWattach ( pge_swath_file_id, TRIM(ADJUSTL(swath_name)) )
!UNUSED!     IF ( pge_swath_id == he5_stat_fail ) CALL error_check ( &
!UNUSED!       0, 1, pge_errstat_fatal, OMSAO_E_HE5SWATTACH, modulename, vb_lev_default, he5stat )
!UNUSED!
!UNUSED!     RETURN
!UNUSED!   END FUNCTION he5_open_readwrite

  SUBROUTINE saopge_geofield_read ( &
      ntimes, nxtrack, geodata_field, geodata, errstat )

    USE OMSAO_precision_module, ONLY: i2, i4, r4
    USE OMSAO_parameters_module,   ONLY: NLINES_MAX, r4_missval
    USE OMSAO_errstat_module
    USE datafields, ONLY: thgt_field
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
    CHARACTER (LEN=*), INTENT (IN) :: geodata_field

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Output variable
    ! ----------------
    REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: geodata

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: locerrstat, iline, nblock
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1) :: geodata_i2

    ! -------------------------------------
    ! Initialize output and local variables
    ! -------------------------------------
    geodata = r4_missval ; locerrstat = pge_errstat_ok

    ! -------------------------------------------------------
    ! Loop over all lines in the file in blocks of NLINES_MAX
    ! -------------------------------------------------------
    ScanLines: DO iline = 0, ntimes-1, NLINES_MAX

      ! --------------------------------------------------------
      ! Check if loop ends before n_times_loop max is exhausted.
      ! --------------------------------------------------------
      nblock = MIN ( NLINES_MAX, ntimes-iline )

      ! ----------------------------------------------------
      ! Read current data block fitting output from HE5 file
      ! ----------------------------------------------------
      he5_start_2d  = (/       0,  iline /)
      he5_stride_2d = (/       1,      1 /)
      he5_edge_2d   = (/ nxtrack, nblock /)

      ! --------------------------------------
      ! Geolocation field
      ! --------------------------------------
      IF (geodata_field .EQ. thgt_field) THEN
        locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(geodata_field)),         &
                                  he5_start_2d, he5_stride_2d, he5_edge_2d, geodata_i2(1:nxtrack,iline:iline+nblock-1) )
        geodata(1:nxtrack,iline:iline+nblock-1) = INT(geodata_i2(1:nxtrack,iline:iline+nblock-1), KIND = r4)
      ELSE
        locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(geodata_field)),         &
                                  he5_start_2d, he5_stride_2d, he5_edge_2d, geodata(1:nxtrack,iline:iline+nblock-1) )
      END IF
    END DO ScanLines

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE saopge_geofield_read

  SUBROUTINE saopge_columninfo_read (              &
      ntimes, nxtrack, saocol, saodco, saorms, saoamf, saofcf, errstat )

    USE OMSAO_precision_module,  ONLY: i2, i4, r8
    USE OMSAO_parameters_module, ONLY: r8_missval, i2_missval, NLINES_MAX
    USE OMSAO_errstat_module
    USE datafields

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saocol, saodco, saorms, saoamf
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saofcf

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: locerrstat, iline, nblock

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    saocol = r8_missval ; saodco = r8_missval
    saorms = r8_missval ; saoamf = r8_missval
    saofcf = i2_missval

    locerrstat = pge_errstat_ok

    ! -------------------------------------------------------
    ! Loop over all lines in the file in blocks of NLINES_MAX
    ! -------------------------------------------------------
    ScanLines: DO iline = 0, ntimes-1, NLINES_MAX

      ! --------------------------------------------------------
      ! Check if loop ends before n_times_loop max is exhausted.
      ! --------------------------------------------------------
      nblock = MINVAL ( (/NLINES_MAX, ntimes-iline /) )

      ! ----------------------------------------------------
      ! Read current data block fitting output from HE5 file
      ! ----------------------------------------------------
      he5_start_2d  = (/       0,  iline /)
      he5_stride_2d = (/       1,      1 /)
      he5_edge_2d   = (/ nxtrack, nblock /)

      ! ------------------------------------------
      ! Column amount, column uncertainty, and RMS
      ! ------------------------------------------
      locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(col_field)),         &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, saocol(1:nxtrack,iline:iline+nblock-1) )
      locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(dcol_field)),        &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, saodco(1:nxtrack,iline:iline+nblock-1) )
      locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(fitrms_field)),      &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, saorms(1:nxtrack,iline:iline+nblock-1) )

      ! ---------------
      ! Air Mass Factor
      ! ---------------
      locerrstat = HE5_SWrdfld (pge_swath_id, TRIM(ADJUSTL(amfmol_field)),      &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, saoamf(1:nxtrack,iline:iline+nblock-1) )

      ! ------------------------
      ! Fitting Convergence Flag
      ! ------------------------
      locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(fitcon_field)),      &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, saofcf(1:nxtrack,iline:iline+nblock-1) )

    END DO ScanLines

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE saopge_columninfo_read

END MODULE
