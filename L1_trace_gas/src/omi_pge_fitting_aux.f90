MODULE omi_pge_fitting_aux

  use tell_module
  use OMSAO_precision_module,  ONLY: i2, i4, r4, r8

  private
  public find_swathrange_by_latitude, read_latitude, &
    find_swathline_by_latitude, convert_tai_to_utc, &
    find_swathline_range, &
    compute_fitting_statistics, &! compute_fitting_statistics_nohe5, &
    omi_set_xtrpix_range, &
    omi_set_fitting_parameters, set_input_pointer_and_versions

  !> Radiance fit QA statistics
  type, public :: fitting_statistics_type
    ! quality_flag array is dimension(nxtrack,0:ntimes-1)
    integer (kind=i2), dimension(:,:), allocatable :: quality_flag
    real (kind=r8) :: col_avg = 0.0_r8
    real (kind=r8) :: dcol_avg = 0.0_r8
    real (kind=r8) :: rms_avg = 0.0_r8
    integer (kind=i4) :: num_col = 0_i4
    integer (kind=i4) :: num_scan_lines = 0_i4
    integer (kind=i4) :: num_crosstrack_pixels = 0_i4
    integer (kind=i4) :: num_input = 0_i4
    integer (kind=i4) :: num_good_input = 0_i4
    integer (kind=i4) :: num_good_output = 0_i4
    integer (kind=i4) :: num_missing = 0_i4
    integer (kind=i4) :: num_suspect_output = 0_i4
    integer (kind=i4) :: num_bad_output = 0_i4
    integer (kind=i4) :: num_converged = 0_i4
    integer (kind=i4) :: num_failed_convergence = 0_i4
    integer (kind=i4) :: num_exceeded_iterations = 0_i4
    integer (kind=i4) :: num_out_of_bounds = 0_i4
    real (kind=r4) :: percent_good_output = 0.0_r4
    real (kind=r4) :: percent_suspect_output = 0.0_r4
    real (kind=r4) :: percent_out_of_bounds = 0.0_r4
    real (kind=r4) :: percent_bad_output = 0.0_r4
    real (kind=r4) :: absolute_percent_missing = 0.0_r4
  end type

CONTAINS
  SUBROUTINE omi_set_fitting_parameters ( pge_idx, errstat )

    USE OMSAO_precision_module
    USE OMSAO_variables_module,  ONLY: l1b_channel
    USE OMSAO_he5_module,        ONLY: swath_base_name, pge_swath_name
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
    INTEGER (KIND=i4), INTENT (inout) :: errstat

    if (errstat /= 0) return

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
      ntimes, nxtrack, xtrange, saocol, saodco, saorms, saofcf, saoamf, amfdiag, &
      sza, vza, fit_stats, errstat )

    USE OMSAO_parameters_module, ONLY: &
      i2_missval, r8_missval, main_qa_good, main_qa_suspect, main_qa_bad, &
      yn_amf, deg2rad
    use optimizer_interface_module, only: &
      opt_convergence_failed, opt_convergence_maxiter_exceeded, opt_convergence_suspect, &
      opt_convergence_good
    USE metadata_tools,  ONLY:  set_automatic_quality_flag
    USE OMSAO_variables_module, ONLY: mdqf_max_good_col, mdqf_min_good_col, &
      mdqf_stddev_sus, mdqf_stddev_bad, mdqf_sza_sus, mdqf_sza_bad, &
      mdqf_amfgeo_sus, mdqf_amfgeo_bad, mdqf_amf_min

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
    INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),     INTENT (IN) :: xtrange
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saocol, saodco, saorms, saoamf
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saofcf, amfdiag
    REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: sza, vza

    ! ----------------
    ! Output variables
    ! ----------------
    type (fitting_statistics_type), intent(inout) :: fit_stats

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Local variables
    ! ----------------
    INTEGER (KIND=i4) :: ix, it, spix, epix
    REAL    (KIND=r8) :: col_avg, rms_avg, dcol_avg
    REAL    (KIND=r8) :: colsig_sus, colsig_bad, amfgeo
    integer (kind=i4) :: num_col, &
      num_good_input, num_good_output, num_missing, num_suspect_output, &
      num_bad_output, num_converged, num_failed_convergence, num_exceeded_iterations, &
      num_out_of_bounds
    character (len=256) :: out_string

    if (errstat /= 0) return

    ! ---------------------------------------------------------
    ! The total number of input samples is simply the number of
    ! pixels in the granule
    ! ---------------------------------------------------------
    fit_stats % num_input = nxtrack * ntimes
    fit_stats % num_crosstrack_pixels = nxtrack
    fit_stats % num_scan_lines = ntimes

    ! ------------------------------------------------------------------
    ! Compute all other fitting statistics variables over two nice loops
    ! ------------------------------------------------------------------
    fit_stats % quality_flag(:,:) = i2_missval

    num_good_input          = 0_i4
    num_good_output         = 0_i4
    num_suspect_output      = 0_i4
    num_bad_output          = 0_i4
    num_out_of_bounds       = 0_i4
    num_converged           = 0_i4
    num_failed_convergence  = 0_i4
    num_exceeded_iterations = 0_i4
    num_missing             = 0_i4

    num_col = 0_i4
    col_avg = 0.0_r8 ; rms_avg = 0.0_r8 ; dcol_avg = 0.0_r8
    DO it = 0, ntimes-1

      spix = xtrange(it,1) ; epix = xtrange(it,2)
      DO ix = spix, epix

        ! Calculate SCD uncertainty limit for filtering
        colsig_sus = (saocol(ix,it) + mdqf_stddev_sus * saodco(ix,it)) * saoamf(ix,it)
        colsig_bad = (saocol(ix,it) + mdqf_stddev_bad * saodco(ix,it)) * saoamf(ix,it)

        ! Calculate geometric AMF for filtering
        amfgeo = &
            1.0_r8 / cos ( real(sza(ix,it),KIND=r8)*deg2rad ) + &
            1.0_r8 / cos ( real(vza(ix,it),KIND=r8)*deg2rad )

        ! If negative inputs for SZA or geometric AMF limits are given in control file,
        ! set these to very large numbers so that the flag does not use them 
        ! as constraints. If mdqf_amf_min is negative, it will already be
        ! ignored.  
        IF (mdqf_amfgeo_sus < 0) mdqf_amfgeo_sus = 1.0e30
        IF (mdqf_amfgeo_bad < 0) mdqf_amfgeo_bad = 1.0e30
        IF (mdqf_sza_sus < 0)    mdqf_sza_sus    = 1.0e30
        IF (mdqf_sza_bad < 0)    mdqf_sza_bad    = 1.0e30

        ! ------------------------------------------------------
        ! The Good: Columns are postive within a defined sigma fitting
        !           uncertainty and the fitting has converged.
        !           For this "sweet spot" we compute the average
        !           fitting statistics.
        ! ------------------------------------------------------

        IF (saofcf(ix,it) == opt_convergence_good) THEN

          num_converged = num_converged + 1

          IF((saocol(ix,it) >= mdqf_min_good_col   ) .AND. &
             (saocol(ix,it) <= mdqf_max_good_col   ) .AND. &            
             (colsig_sus    >= 0.0_r8              ) .AND. &
             (sza(ix,it)    <= mdqf_sza_sus        ) .AND. &
             (amfgeo        <= mdqf_amfgeo_sus     ) .AND. &
             (saoamf(ix,it) >= mdqf_amf_min        ) )  THEN
 
            fit_stats % quality_flag(ix,it) = main_qa_good

            num_good_input = num_good_input + 1
            num_good_output = num_good_output + 1

            col_avg  = col_avg  + saocol(ix,it)
            dcol_avg = dcol_avg + saodco(ix,it)
            rms_avg  = rms_avg  + saorms(ix,it)
            num_col  = num_col + 1

            CYCLE
          END IF
        END IF

        ! ----------------------------------------------------------
        ! The Bad: Fitting hasn't converged or columns are negative
        !          within a defined sigma fitting uncertainty. Note that
        !          pixels can count towards both the number of out-
        !          of bounds and the failed convergence samples.
        ! ----------------------------------------------------------
        IF ((saofcf(ix,it) > i2_missval .AND. saofcf(ix,it) < 0_i2) .OR. &
            (saocol(ix,it) > r8_missval .AND. colsig_bad < 0.0_r8 ) .OR. &
            (btest(amfdiag(ix,it),yn_amf)                         ) .OR. &
            (sza(ix,it) > mdqf_sza_bad                            ) .OR. &
            (amfgeo > mdqf_amfgeo_bad                             ) ) THEN

          fit_stats % quality_flag(ix,it) = main_qa_bad

          num_good_input = num_good_input + 1
          num_bad_output = num_bad_output + 1

          IF (saocol(ix,it) > r8_missval .AND. (colsig_bad < 0.0_r8 .OR. &
             sza(ix,it) > mdqf_sza_bad .OR. amfgeo > mdqf_amfgeo_bad)) &
            num_out_of_bounds = num_out_of_bounds + 1
          IF (saofcf(ix,it) == opt_convergence_failed .or. saofcf(ix,it) == opt_convergence_maxiter_exceeded ) &
            num_failed_convergence = num_failed_convergence  + 1
          IF (saofcf(ix,it) == opt_convergence_maxiter_exceeded) &
            num_exceeded_iterations = num_exceeded_iterations + 1

          CYCLE
        END IF

        ! ----------------------------------------------------------
        ! The Ugly: Whatever is left (outside plain missing columns)
        ! ----------------------------------------------------------
        IF (saocol(ix,it) > r8_missval ) THEN

          IF ((saofcf(ix,it) == opt_convergence_suspect       ) .OR. &
              (colsig_sus < 0.0_r8 .AND. colsig_bad >= 0.0_r8 ) .OR. &
              (saocol(ix,it) > mdqf_max_good_col              ) .OR. &
              (saocol(ix,it) < mdqf_min_good_col              ) .OR. &
              (sza(ix,it) > mdqf_sza_sus .AND. sza(ix,it) <= mdqf_sza_bad) .OR. &   
              (amfgeo > mdqf_amfgeo_sus .AND. amfgeo <= mdqf_amfgeo_bad  ) .OR. &
              (saoamf(ix,it) <  mdqf_amf_min                  ) ) THEN

            fit_stats % quality_flag(ix,it) = main_qa_suspect

            num_good_input     = num_good_input     + 1
            num_suspect_output = num_suspect_output + 1

            CYCLE
          END IF

        ELSE

          ! ----------------------------------------------------------
          ! The Missing: Not processed because of either missing input
          !              or restrictions on lat, lon, sza, etc.
          ! ----------------------------------------------------------
          num_missing = num_missing + 1

        END IF

      END DO
    END DO

    ! --------------------------------------------
    ! Now we can compute averages and percentages,
    ! and write out the final statistics
    ! --------------------------------------------

    IF (num_col > 0) THEN
      col_avg  = col_avg  / num_col
      rms_avg  = rms_avg  / num_col
      dcol_avg = dcol_avg / num_col
    END IF

    fit_stats % percent_good_output = 100_r4 * &
      REAL(num_good_output, KIND=r4) / &
      MAX ( 1.0_r4, REAL(num_good_input,  KIND=r4) )

    fit_stats % percent_bad_output = 100_r4 * &
      REAL(num_bad_output, KIND=r4) / &
      MAX ( 1.0_r4, REAL(num_good_input, KIND=r4) )

    fit_stats % percent_suspect_output = 100.0_r4 * &
      REAL(num_suspect_output, KIND=r4) / &
      MAX ( 1.0_r4, REAL(num_good_input, KIND=r4) )

    fit_stats % percent_out_of_bounds = 100.0_r4 * &
      REAL(num_out_of_bounds, KIND=r4) / &
      MAX ( 1.0_r4, REAL(num_good_input, KIND=r4) )

    fit_stats % absolute_percent_missing = 100_r4 * &
      REAL(num_missing, KIND=r4) / &
      MAX ( 1.0_4, REAL(fit_stats % num_input, KIND=r4) )

    ! ------------------------------------------------------------------------
    ! With the above information we can easily determine the Automatic QA Flag
    ! ------------------------------------------------------------------------
    CALL set_automatic_quality_flag (fit_stats % percent_good_output)

    WRITE (out_string, '(A, 3(1PE15.5))')'Col-DCol-RMS: ', col_avg, dcol_avg, rms_avg
    call tell_log (0, out_string)

    WRITE (out_string, '(A, I7,A,I7,A,F7.1,A)')'Statistics:   ', &
      MAX(num_good_output,0), ' of ', MAX(num_good_input,0), ' converged and in bounds - ', &
      MAX(fit_stats % percent_good_output, 0.0), '%'
    call tell_log (0, out_string)

    WRITE (out_string, '(A, I7,A,I7,A,F7.1,A)')'Statistics:   ', &
      MAX(num_converged,0), ' of ', MAX(num_good_input,0), ' converged - ', &
      MAX(100_r4 * REAL(num_converged, KIND=r4) / &
      MAX ( 1.0_r4, REAL(num_good_input,  KIND=r4) ), 0.0), '%'
    call tell_log (0, out_string)

    WRITE (out_string, '(A, i10)') 'num_col =', num_col
    call tell_log (0, out_string)

    fit_stats % col_avg  = col_avg
    fit_stats % dcol_avg = dcol_avg
    fit_stats % rms_avg  = rms_avg
    fit_stats % num_col  = num_col

    fit_stats % num_good_input          = num_good_input
    fit_stats % num_good_output         = num_good_output
    fit_stats % num_suspect_output      = num_suspect_output
    fit_stats % num_bad_output          = num_bad_output
    fit_stats % num_out_of_bounds       = num_out_of_bounds
    fit_stats % num_converged           = num_converged
    fit_stats % num_failed_convergence  = num_failed_convergence
    fit_stats % num_exceeded_iterations = num_exceeded_iterations
    fit_stats % num_missing             = num_missing

    RETURN

  END SUBROUTINE compute_fitting_statistics

  SUBROUTINE set_input_pointer_and_versions ( )

    USE OMSAO_indices_module,      ONLY: &
      l1b_radiance_lun, l1b_irradiance_lun, &
      voc_amf_luns, voc_omicld_idx
    USE OMSAO_he5_module, ONLY: n_lun_inp, lun_input
    use ctrlvars, only: yn_I0

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
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
    ! (a) Whether a radiance reference from a granule other than the one
    !     being processed is being used.
    ! --------------------------------------------------------------------
    ! (0) Processed granule
    ! ---------------------
    n_lun_inp = 1
    lun_input(n_lun_inp) = l1b_radiance_lun
    ! --------------------
    ! (a) Solar Irradiance
    ! --------------------
    IF ( .NOT. yn_I0 ) THEN
      n_lun_inp = n_lun_inp + 1
      lun_input(n_lun_inp) = l1b_irradiance_lun
    END IF

    ! -----------------
    ! Cloud information
    ! -----------------
    n_lun_inp            = n_lun_inp + 1
    lun_input(n_lun_inp) = voc_amf_luns(voc_omicld_idx)

    RETURN
  END SUBROUTINE set_input_pointer_and_versions


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

  SUBROUTINE omi_set_xtrpix_range (                 &
      nTimes, nXtrack, pixnum_limits, omi_binfac,  &
      omi_xtrpix_range, first_wc_pix, last_wc_pix, &
      errstat )

    USE OMSAO_precision_module,  ONLY: i1, i4
    USE OMSAO_omidata_module,    ONLY: szoom_mode, global_mode, gzoom_spix, gzoom_epix
    !USE OMSAO_errstat_module

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
    INTEGER (KIND=i4)          :: first_pix, last_pix, i !, locerrstat

    ! -----------------------
    ! Name of this subroutine
    ! -----------------------
    !CHARACTER (LEN=20), PARAMETER :: modulename = 'omi_set_xtrpix_range'

    if (errstat /= 0) return

    ! ---------------------------
    ! Initialize return variables
    ! ---------------------------
    omi_xtrpix_range(0:nTimes-1,1:2) = -1
    first_wc_pix                     = -1
    last_wc_pix                      = -1

    !locerrstat = pge_errstat_ok

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

    !errstat = MAX ( errstat, locerrstat )
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
    !USE OMSAO_errstat_module
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
    INTEGER (KIND=i4)                          :: fpix, lpix, midnum !, locerrstat, estat
    REAL    (KIND=r4)                          :: midlat
    INTEGER (KIND=i4), DIMENSION (0:nt-1, 1:2) :: xtrange
    INTEGER (KIND=i1), DIMENSION (0:nt-1)      :: binfac
    LOGICAL,           DIMENSION (0:nt-1)      :: ynzoom

    if (errstat /= 0) return
    !locerrstat = pge_errstat_ok
    !estat      = pge_errstat_ok

    ! ----------------------------------------------------------------
    ! Read preparatory arrays for determining the range of swath lines
    ! that fall within the desired latitude interval.
    ! ----------------------------------------------------------------
    CALL omi_read_binning_factor ( &
      l1bfile, l1bswath, nt, binfac(0:nt-1), ynzoom(0:nt-1), errstat )

    CALL omi_set_xtrpix_range ( &
      nt, nx, pixnum_lim(3:4), binfac(0:nt-1), &
      xtrange(0:nt-1,1:2), fpix, lpix, errstat    )

    if (errstat /= 0) return

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

    !errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE find_swathline_range

  SUBROUTINE find_swathline_by_latitude ( &
      nxrr, sline, eline, latr4, lat, xtrange, lnum, was_found )

    USE OMSAO_precision_module, ONLY: i4, r4
    USE OMSAO_parameters_module, ONLY: r4_missval
    !USE OMSAO_errstat_module
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
    !INTEGER (KIND=i4), DIMENSION (sline:eline,2),      INTENT (IN) :: xtrange
    INTEGER (KIND=i4), DIMENSION (:,:),                INTENT (IN) :: xtrange
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
      j1, j2, icnt, iline, fpix, lpix !, locerr
    character (len=256) :: log_msg

    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    lnum = -1
    was_found = .FALSE.

    ! --------------------------------------------------------------------------
    ! First, start a bisection of the [0, NLINES-1] interval to find the closest
    ! match in latitude to the mipoint of the latitude regime to average.
    ! --------------------------------------------------------------------------
    ! FIXME (JCH) calling routine already provides eline=nlines-1, so
    ! setting j2=eline-1 seems excessive.
    j1 = sline ; j2 = eline-1  ;  icnt = 0
    FindLine: DO WHILE ( .NOT. was_found )
      icnt  = icnt + 1
      iline = (j1 + j2) / 2

      write(log_msg, *)'find_swathline_by_latitude: looking for lat=', &
        lat,' iline,j1,j2=',iline,j1,j2
      call tell_log (3, log_msg)

      ! -----------------------------------------------------------------------
      ! Get first and last pixel.
      ! -----------------------------------------------------------------------
      ! FIXME: (JCH) This can fail if iline=0. It works only if we have
      !        enough latitude points to put the midpoint at iline>=1.
      fpix = xtrange(iline,1)
      lpix = xtrange(iline,2)

      write(log_msg, *)'find_swathline_by_latitude: fpix,lpix=',fpix,lpix
      call tell_log (3, log_msg)

      IF ( iline < sline .OR. iline > eline ) THEN
        !unused! locerr = pge_errstat_error
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
    !USE OMSAO_errstat_module

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
    !use l1bread, only: l1bread_open_swath, l1bread_close, L1B_Object_Type, &
    !  l1bread_get2d_r4
    use tio_module
    use netcdf, only : nf90_nowrite

    implicit none
    CHARACTER (LEN=*),     INTENT (IN) :: l1bfile, l1bswath
    INTEGER (KIND=i4),     INTENT (IN) :: tstart, ntimes
    integer, intent(inout) :: errstat
    REAL    (KIND=r4), DIMENSION(:,:), INTENT (out) :: latr4

    !type (L1B_Object_Type) :: l1bobj
    type (tiof_file_type) :: tio_l1obj
    integer :: nxtrack

    if (errstat /= 0) return

    !call l1bread_open_swath (l1bfile, l1bswath, l1bobj, errstat)
    !if (errstat /= 0) return
    !if (size(latr4, 1) /= l1bobj%num_xtrack) then
    !  call tell_error (tell_runtime_error, "read_latitude: nxtrack dimension is not correct", errstat)
    !  call l1bread_close (l1bobj)
    !  return
    !endif
    !call l1bread_get2d_r4 (l1bobj, "latitude", tstart, ntimes, latr4, errstat)
    !call l1bread_close (l1bobj)
    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_l1obj, l1bswath, errstat)
    call tiof_inq_dimlen (tio_l1obj, "xtrack", nxtrack, errstat)
    if (errstat /= 0) return
    if (size(latr4, 1) /= nxtrack) then
      call tell_error (tell_io_read_error, &
                       "read_latitude: nxtrack dimension is not correct", errstat)
      call tiof_close (tio_l1obj, errstat)
      return
    endif
    call tiof_get2d_r4 (tio_l1obj, "latitude", [tstart,0], [ntimes,-1], latr4(1:nxtrack,1:ntimes), errstat)
    call tiof_close (tio_l1obj, errstat)

    return

  end subroutine read_latitude

!unused  SUBROUTINE compute_fitting_statistics_nohe5 ( &
!unused      pge_idx, ntimes, nxtrack, xtrange, saocol, saodco, saorms, &
!unused      saofcf, saomqf, errstat )
!unused
!unused    USE OMSAO_precision_module,  ONLY: i2, i4, r4, r8
!unused    USE OMSAO_parameters_module, ONLY: &
!unused      i2_missval, r8_missval, main_qa_good, main_qa_suspect, main_qa_bad
!unused    use optimizer_interface_module, only: &
!unused      opt_convergence_failed, opt_convergence_maxiter_exceeded, opt_convergence_suspect, &
!unused      opt_convergence_good
!unused    USE metadata_tools, ONLY:  QAPercentMissingData, QAPercentOutofBoundsData, &
!unused      set_automatic_quality_flag
!unused
!unused    USE OMSAO_he5_module,       ONLY:  &
!unused      NrOfInputSamples, NrofGoodOutputSamples, NrofSuspectOutputSamples,        &
!unused      NrofBadOutputSamples, NrofConvergedSamples, NrofFailedConvergenceSamples, &
!unused      NrofExceededIterationsSamples, NrofOutofBoundsSamples, NrofMissingSamples, &
!unused      NrofGoodInputSamples, NrofSuspectOutputSamples, NrofBadOutputSamples,      &
!unused      NrofConvergedSamples, NrofFailedConvergenceSamples, &
!unused      PercentGoodOutputSamples, PercentSuspectOutputSamples, &
!unused      PercentBadOutputSamples, &
!unused      AbsolutePercentMissingSamples
!unused    USE OMSAO_errstat_module,   ONLY: vb_lev_screen, pge_errstat_ok
!unused    USE OMSAO_variables_module, ONLY: verb_thresh_lev, max_good_col
!unused
!unused    IMPLICIT NONE
!unused
!unused    ! ---------------
!unused    ! Input variables
!unused    ! ---------------
!unused    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ntimes, nxtrack
!unused    INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),     INTENT (IN) :: xtrange
!unused    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saocol, saodco, saorms
!unused    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saofcf
!unused
!unused    ! ----------------
!unused    ! Output variables
!unused    ! ----------------
!unused    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: saomqf
!unused
!unused    ! -----------------
!unused    ! Modified variable
!unused    ! -----------------
!unused    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
!unused
!unused    ! ----------------
!unused    ! Local variables
!unused    ! ----------------
!unused    INTEGER (KIND=i4) :: locerrstat, ix, it, spix, epix
!unused    REAL    (KIND=r4) :: PercentOutofBoundsSamples
!unused    REAL    (KIND=r8) :: fitcol_avg, rms_avg, dfitcol_avg, nfitcol
!unused    REAL    (KIND=r8) :: col2sig, col3sig
!unused
!unused    locerrstat = pge_errstat_ok
!unused
!unused    ! ---------------------------------------------------------
!unused    ! The total number of input samples is simply the number of
!unused    ! pixels in the granule
!unused    ! ---------------------------------------------------------
!unused    NrofInputSamples = nxtrack*ntimes
!unused
!unused    ! ------------------------------------------------------------------
!unused    ! Compute all other fitting statistics variables over two nice loops
!unused    ! ------------------------------------------------------------------
!unused    saomqf                        = i2_missval
!unused    NrofGoodInputSamples          = 0_i4
!unused    NrofGoodOutputSamples         = 0_i4
!unused    NrofSuspectOutputSamples      = 0_i4
!unused    NrofBadOutputSamples          = 0_i4
!unused    NrofOutOfBoundsSamples        = 0_i4
!unused    NrofConvergedSamples          = 0_i4
!unused    NrofFailedConvergenceSamples  = 0_i4
!unused    NrofExceededIterationsSamples = 0_i4
!unused    NrofMissingSamples            = 0_i4
!unused
!unused    nfitcol    = 0.0_r8
!unused    fitcol_avg = 0.0_r8 ; rms_avg = 0.0_r8 ; dfitcol_avg = 0.0_r8
!unused    DO it = 0, ntimes-1
!unused
!unused      spix = xtrange(it,1) ; epix = xtrange(it,2)
!unused      DO ix = spix, epix
!unused
!unused        col2sig = saocol(ix,it)+2.0_r8*saodco(ix,it)
!unused        col3sig = saocol(ix,it)+3.0_r8*saodco(ix,it)
!unused
!unused        ! ------------------------------------------------------
!unused        ! The Good: Columns are postive within two sigma fitting
!unused        !           uncertainty and the fitting has converged.
!unused        !           For this "sweet spot" we compute the average
!unused        !           fitting statistics.
!unused        ! ------------------------------------------------------
!unused        IF ( (saofcf(ix,it) == opt_convergence_good) .AND. &
!unused            (saocol(ix,it)      >  r8_missval                       ) .AND. &
!unused            (ABS(saocol(ix,it)) <= max_good_col                     ) .AND. &
!unused            (col2sig            >= 0.0_r8                           ) ) THEN
!unused
!unused          saomqf(ix,it) = main_qa_good
!unused
!unused          NrofGoodInputSamples  = NrofGoodInputSamples   + 1
!unused          NrofConvergedSamples  = NrofConvergedSamples   + 1
!unused          NrofGoodOutputSamples = NrofGoodOutputSamples + 1
!unused
!unused          fitcol_avg  = fitcol_avg  + saocol(ix,it)
!unused          dfitcol_avg = dfitcol_avg + saodco(ix,it)
!unused          rms_avg     = rms_avg     + saorms(ix,it)
!unused          nfitcol     = nfitcol     + 1.0_r8
!unused          !write(*,'(1P3E15.5,2I5)') saocol(ix,it),saorms(ix,it), saodco(ix,it), ix, it
!unused
!unused          CYCLE
!unused        END IF
!unused
!unused        ! ----------------------------------------------------------
!unused        ! The Bad: Fitting hasn't converged or columns are negative
!unused        !          within three sigma fitting uncertainty. Note that
!unused        !          pixels can count towards both the number of out-
!unused        !          of bounds and the failed convergence samples.
!unused        ! ----------------------------------------------------------
!unused        IF ( (saofcf(ix,it)      > i2_missval .AND. saofcf(ix,it) < 0_i2) .OR. &
!unused            (saocol(ix,it)      > r8_missval .AND. col3sig < 0.0_r8    ) ) THEN
!unused
!unused          saomqf(ix,it) = main_qa_bad
!unused
!unused          NrofGoodInputSamples = NrofGoodInputSamples + 1
!unused          NrofBadOutputSamples = NrofBadOutputSamples + 1
!unused
!unused          IF ( saocol(ix,it) > r8_missval .AND. col3sig < 0.0_r8 ) &
!unused            NrofOutofBoundsSamples        = NrofOutofBoundsSamples        + 1
!unused          IF ( saofcf(ix,it) == opt_convergence_failed .or. saofcf(ix,it) == opt_convergence_maxiter_exceeded) &
!unused            NrofFailedConvergenceSamples  = NrofFailedConvergenceSamples  + 1
!unused          IF ( saofcf(ix,it) == opt_convergence_maxiter_exceeded)  &
!unused            NrofExceededIterationsSamples = NrofExceededIterationsSamples + 1
!unused
!unused          CYCLE
!unused        END IF
!unused
!unused        ! ----------------------------------------------------------
!unused        ! The Ugly: Whatever is left (outside plain missing columns)
!unused        ! ----------------------------------------------------------
!unused        IF ( saocol(ix,it) > r8_missval ) THEN
!unused
!unused          IF ( (saofcf(ix,it) == opt_convergence_suspect) .OR. &
!unused              (col2sig <  0.0_r8  .AND. col3sig >= 0.0_r8                      ) .OR. &
!unused              (ABS(saocol(ix,it)) > max_good_col                               ) ) THEN
!unused
!unused            saomqf(ix,it) = main_qa_suspect
!unused
!unused            NrofGoodInputSamples     = NrofGoodInputSamples     + 1
!unused            NrofSuspectOutputSamples = NrofSuspectOutputSamples + 1
!unused
!unused            CYCLE
!unused          END IF
!unused
!unused        ELSE
!unused
!unused          ! ----------------------------------------------------------
!unused          ! The Missing: Not processed because of either missing input
!unused          !              or restrictions on lat, lon, sza, etc.
!unused          ! ----------------------------------------------------------
!unused          NrofMissingSamples = NrofMissingSamples + 1
!unused
!unused        END IF
!unused
!unused      END DO
!unused    END DO
!unused
!unused    ! --------------------------------------------
!unused    ! Now we can compute averages and percentages,
!unused    ! and write out the final statistics
!unused    ! --------------------------------------------
!unused
!unused    IF ( nfitcol >= 1.0_r8 ) THEN
!unused      fitcol_avg  = fitcol_avg  / nfitcol
!unused      rms_avg     = rms_avg     / nfitcol
!unused      dfitcol_avg = dfitcol_avg / nfitcol
!unused    END IF
!unused
!unused    PercentGoodOutputSamples      = 100_r4    * &
!unused      REAL(NrofGoodOutputSamples, KIND=r4) / &
!unused      MAX ( 1.0_r4, REAL(NrofGoodInputSamples,  KIND=r4) )
!unused
!unused    PercentBadOutputSamples       = 100_r4        * &
!unused      REAL(NrofBadOutputSamples, KIND=r4) / &
!unused      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )
!unused
!unused    PercentSuspectOutputSamples   =  100.0_r4         * &
!unused      REAL(NrofSuspectOutputSamples, KIND=r4) / &
!unused      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )
!unused
!unused    PercentOutofBoundsSamples     =  100.0_r4         * &
!unused      REAL(NrofOutofBoundsSamples, KIND=r4) / &
!unused      MAX ( 1.0_r4, REAL(NrofGoodInputSamples, KIND=r4) )
!unused
!unused    AbsolutePercentMissingSamples = 100_r4 * &
!unused      REAL(NrofMissingSamples, KIND=r4) / &
!unused      MAX ( 1.0_4, REAL(NrofInputSamples, KIND=r4) )
!unused
!unused    QAPercentMissingData     = NINT ( AbsolutePercentMissingSamples, KIND=i4 )
!unused    QAPercentOutofBoundsData = NINT ( PercentOutofBoundsSamples,     KIND=i4 )
!unused
!unused    ! ------------------------------------------------------------------------
!unused    ! With the above information we can easily determine the Automatic QA Flag
!unused    ! ------------------------------------------------------------------------
!unused    CALL set_automatic_quality_flag ( PercentGoodOutputSamples )
!unused
!unused    IF ( verb_thresh_lev >= vb_lev_screen ) THEN
!unused      WRITE (*, '(A, 3(1PE15.5))')          'Col-DCol-RMS: ', fitcol_avg, dfitcol_avg, rms_avg
!unused      WRITE (*, '(A, I7,A,I7,A,F7.1,A)')  'Statistics:   ', &
!unused        MAX(NrofGoodOutputSamples,0), ' of ', MAX(NrofGoodInputSamples,0), ' converged - ', &
!unused        MAX(PercentGoodOutputSamples, 0.0), '%'
!unused      WRITE (*, '(A, F7.1)') 'Nfitcol =', nfitcol
!unused    END IF
!unused
!unused    errstat = locerrstat
!unused    RETURN
!unused
!unused  END SUBROUTINE compute_fitting_statistics_nohe5

END MODULE
