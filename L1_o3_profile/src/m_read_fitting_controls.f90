!
module m_read_fitting_controls

  public read_fitting_control_file
  private !get_mols_for_fitting

contains

  ! *********************** Modification History **********************
  ! Xiong Liu; July, 2003  (!xliu)
  ! 1. Read whether to do ozone profile retrieval and set a flag
  ! 2. Read additional fitting control variables if ozprof_flag is set
  ! 3. Read an option about variable slit width, add yn_varyslit variable
  !    n_slit_itnerval, slit_fname, slit_redo, wavcal_redo, wavcal_fname,
  !    in USE OMSAO_variables_module
  ! 4. Add 1 nm more for winwav_min, winwav_max to avoid interpolation
  !    out of bounds for solar spectrum calibration
  ! 5. Read option use_meas_sig
  ! 6. Read option use_pixel_bin
  ! *******************************************************************

  SUBROUTINE read_fitting_control_file( fit_ctrl_unit, pge_error_status)

    USE PGS_PC_class
    USE OMI_LUN_set
    ! ***********************************************************
    !
    !   Read fitting control parameters from input control file
    !
    ! ***********************************************************

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: &
         max_rs_idx, max_calfit_idx, mns_idx, mxs_idx,&
         calfit_strings, radfit_strings, refspec_strings, & 
         genline_str, socline_str, racline_str,rafline_str, eoi3str, & 
         solar_idx, us1_idx, us2_idx, shift_offset, &
         com_idx, com1_idx, com2_idx, com3_idx, &
         comfidx,  cm1fidx, cm2fidx,  cm3fidx, &
         comvidx,  cm1vidx,  cm2vidx, cm3vidx, & 
         hwe_idx, asy_idx, vgl_idx, hwr_idx, spk_idx, shi_idx, squ_idx,& 
         wr0_idx, wr7_idx , &
         so2_idx, o2_idx, o2o2_idx, h2o_idx, &
         which_instrument, max_instrument_idx, instrument_idx, &
         gome_idx, omi_idx, scia_idx, gome2_idx, tempo_idx
    USE OMSAO_parameters_module,   ONLY: mswath, maxchlen, maxwin, max_fit_pts, &
                                         zerospec_string

    USE OMSAO_variables_module,    ONLY: l1l2inp_unit,use_backup, use_solcomp,    &
         l1b_irrad_filename, l1b_rad_filename, l2_filename, l2_cld_filename, &
         avg_solcomp, avgsol_allorb, &
         n_fincol_idx, max_itnum_sol, &
         max_itnum_rad,  yn_smooth, yn_doas, weight_sun, fitvar_sol_saved, &
         fitvar_sol_init, fitvar_rad_init, fitvar_rad_saved, yn_varyslit,  &
         wavcal, which_slit,slit_name, n_slit_step, slit_fname, slit_redo, wavcal_redo, &
         wavcal_fname, swavcal_fname, use_meas_sig,  smooth_slit, &
         slit_fit_pts, wavcal_fit_pts, n_wavcal_step, wcal_bef_coadd,       &
         wavcal_sol, mask_fitvar_rad, mask_fitvar_sol, n_fitvar_sol, renorm,&
         weight_rad, szamax, zatmos, slit_rad, rslit_fname,  n_fitvar_rad,  &
         lo_sunbnd, up_sunbnd, lo_sunbnd_init,up_sunbnd_init, lo_radbnd,    &
         up_radbnd, refspec_fname, &
         radwavcal_freq, tol,  epsrel,  epsabs,  epsx, pm_one, phase, &
         linenum_lim,  pixnum_lim, coadd_uv2,       &
         fitvar_rad_str, winwav_min, winwav_max, &
         have_amftable, have_undersampling, winlim, sol_identifier, &
         rad_identifier, numwin, nviswin,scnwrt, calwrt,calscn, rtmdbg,& 
         band_selectors, do_bandavg, &
         n_band_avg, n_band_samp, outdir, atmdbdir, tabdir, refdbdir, &
         rmask_fitvar_rad, database_indices, slit_trunc_limit, &
         reduce_resolution, redsampr, redlam, reduce_slit, rm_mgline, &
         redfixwav,use_redfixwav, nredfixwav, redfixwav_fname, radnhtrunc, &
         refnhextra, l2_swathname, fitvar_rad_unit, l1b_rad_filename, &
         correct_merr, xbin_decerr, ybin_decerr, &
         do_xbin, do_ybin, nxbin, nybin, rmask_fitvar_sol,l2_hdf_flag, redslw, &
         upper_wvls, lower_wvls, upper_spec, lower_spec, retlbnd, retubnd,& 
         nswath, orbnum, orbnumsol,num_param, ncoadd, do_ch2reso
    USE OMSAO_errstat_module
    USE omi_read_l1b_data, only: find_scan_line_range
    USE OMSAO_omidata_module, ONLY:  &
        mswath_omi=>mswath, nxomi_max=>nxtrack_max, ntomi_max=>ntimes_max, &
        upper_wvls_omi=>upper_wvls, lower_wvls_omi=>lower_wvls, &
        upper_spec_omi=>upper_spec, lower_spec_omi=>lower_spec, omisol_version
    USE OMSAO_tmpodata_module, ONLY: &
        mswath_tmpo=>mswath, nxtmpo_max=>nxtrack_max,nttmpo_max=>ntimes_max, &
        upper_wvls_tmpo=>upper_wvls, lower_wvls_tmpo=>lower_wvls, &
        upper_spec_tmpo=>upper_spec, lower_spec_tmpo=>lower_spec
    USE ozprof_data_module, ONLY: radcalwrt,ozprof_flag, ozprof_input_fname,& 
        nlay, nos, nsh, nsl, do_simu, nfgas, ozfit_end_index, ozfit_start_index, &
        use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs 
    USE UTIL_tools_class
    USE m_utilities, only: get_substring, string2index,check_for_endofinput, &
                           skip_to_filemark
    USE m_read_ozprof_input

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,           INTENT (IN)    :: fit_ctrl_unit

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER,           INTENT (OUT)   :: pge_error_status
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: nxtrack_max, ntimes_max
    INTEGER :: i, j, k, file_read_stat, sidx, ridx, idx, cldorb, ntsh, ext
    INTEGER, DIMENSION(2)    :: pixlim, linelim
    CHARACTER (LEN=maxchlen) :: tmpchar, l1l2_files, fname
    CHARACTER (LEN=3)        :: idxchar, xbinchar, ybinchar
    CHARACTER (LEN=5)        :: idxchar1, cldorbc
    CHARACTER (LEN=4)        :: slinechar, elinechar
    CHARACTER (LEN=4)        :: sxchar, exchar
    LOGICAL                  :: yn_eoi, select_lonlat, &
         rw_l1l2_pcf
    REAL      (KIND=dp)      :: vartmp, lotmp, uptmp, slat, elat, slon, elon

    CHARACTER (LEN=100)      :: fit_ctrl_file
    CHARACTER (LEN=5)        :: orbc, orbcsol
    CHARACTER (LEN=3) :: sn
    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=30), PARAMETER :: modulename = 'read_fitting_control_file'

    ! ========================
    ! Error handling variables
    ! ========================
    INTEGER :: errstat, version!, ios

    ! -----------------------
    ! GOME specific additions
    ! -----------------------
    CHARACTER (LEN=26), PARAMETER :: lm_instrument = &
         '*Satellite instrument name'
    CHARACTER (LEN=17), PARAMETER :: lm_l2hdf      = '*Write HDF output'
    CHARACTER (LEN=31), PARAMETER :: lm_bandselect = &
         '*OMI radiance bands to be used'
    CHARACTER (LEN=28), PARAMETER :: lm_reduceres  = &
         '*Reduce spectral resolution'
    CHARACTER (LEN=6), PARAMETER :: lm_debug  = "*debug"

    CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
    INTEGER :: nc
    ! =================================
    ! External OMI and Toolkit routines
    ! =================================
    INTEGER :: OMI_SMF_setmsg

    pge_error_status = pge_errstat_ok
    coadd_uv2 = .FALSE.
    do_ch2reso= .FALSE.
    ! -----------------------------------------------------------
    ! Initialize array with reference spectrum names to Zero_Spec
    ! -----------------------------------------------------------
    DO j = solar_idx, max_rs_idx
      refspec_fname(j) = 'OMSAO_Zero_Spec.dat'
    END DO

    fit_ctrl_file(:) = ' '
    msg(:) = ' '
 
    ! -----------------------------------------------------------
    ! Open fitting control file
    ! -----------------------------------------------------------
    version = 1
    errstat = PGS_PC_getreference( fit_ctrl_unit, version, fit_ctrl_file )
    IF( errstat /= PGS_S_SUCCESS ) THEN
      WRITE(msg, '(A,I10,I4)') 'get file from lun=', fit_ctrl_unit, version
      errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    ELSE
      errstat = OMI_SMF_setmsg(OMI_S_SUCCESS, &
           'fit_ctrl_file ='//TRIM(fit_ctrl_file), modulename, 0)
    END IF

    OPEN ( UNIT=fit_ctrl_unit, FILE=TRIM(ADJUSTL(fit_ctrl_file)), &
         STATUS='OLD', IOSTAT=errstat)
    IF ( errstat /= pge_errstat_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
           TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
      WRITE(www_lun) 'failed to open:'//ADJUSTL(TRIM(fit_ctrl_file))
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    ! -----------------------------------------------------------
    ! Get database directories fro external
    ! -----------------------------------------------------------
    ! tabdir
    version = 1
    errstat = PGS_PC_getreference( TAB_DIR_LUN, version, tabdir )
    IF( errstat /= PGS_S_SUCCESS ) THEN
      WRITE(msg, '(A,I10,I4)') 'get file from lun=', TAB_DIR_LUN, version
      errstat = OMI_SMF_setmsg (OMI_E_INPUT, msg, modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    ELSE
      errstat = OMI_SMF_setmsg(OMI_S_SUCCESS, 'tabdir ='//TRIM(tabdir), &
           modulename, 0)
    END IF

    ! atmdbdir
    version = 1
    errstat = PGS_PC_getreference( ATMOSDB_DIR_LUN, version, atmdbdir )
    IF( errstat /= PGS_S_SUCCESS ) THEN
      WRITE(msg, '(A,I10,I4)') 'get file from lun=', ATMOSDB_DIR_LUN, version
      errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    ELSE
      errstat = OMI_SMF_setmsg(OMI_S_SUCCESS, TRIM(atmdbdir), modulename, 0)
    END IF
    IF( scnwrt ) WRITE (*, '(A)') TRIM(atmdbdir)

    ! refdbdir
    version = 1
    errstat = PGS_PC_getreference( REFDB_DIR_LUN, version, refdbdir )
    IF( errstat /= PGS_S_SUCCESS ) THEN
      WRITE(msg, '(A,I10,I4)') 'get file from lun=', REFDB_DIR_LUN, version
      errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    ELSE
      errstat = OMI_SMF_setmsg(OMI_S_SUCCESS, TRIM(refdbdir), modulename, 0)
    END IF

    ! -----------------------------------------------------------
    ! Position cursor to read instrument name
    ! -----------------------------------------------------------
    REWIND ( fit_ctrl_unit )
    CALL skip_to_filemark ( fit_ctrl_unit, lm_instrument, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_instrument, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF
    READ (fit_ctrl_unit, '(A)') tmpchar
    CALL string2index ( which_instrument, max_instrument_idx, tmpchar, &
         instrument_idx )

    ! Set up according to instrument 
    ! If Tempo synthetic data, we need to behave like OMI, but with an
    ! over-ride for some settings (e.g., uv2_coadd)
    IF (instrument_idx == omi_idx) THEN 
      nxtrack_max = nxomi_max ; ntimes_max = ntomi_max
      upper_wvls = upper_wvls_omi(1:mswath_omi)
      lower_wvls = lower_wvls_omi(1:mswath_omi)
      upper_spec = upper_spec_omi
      lower_spec = lower_spec_omi
    ELSE IF (instrument_idx == tempo_idx) THEN
      nxtrack_max = nxtmpo_max ; ntimes_max = nttmpo_max
      upper_wvls = upper_wvls_tmpo(1:mswath_tmpo)
      lower_wvls = lower_wvls_tmpo(1:mswath_tmpo)
      upper_spec = upper_spec_tmpo
      lower_spec = lower_spec_tmpo
      !instrument_idx = omi_idx
    ENDIF

    ! -----------------------------------------------------------
    ! Position cursor to debug options
    ! -----------------------------------------------------------
    REWIND ( fit_ctrl_unit )
    CALL skip_to_filemark ( fit_ctrl_unit, lm_debug, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_debug, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      WRITE(www_lun) 'failed to read:', lm_debug 
      RETURN
    ENDIF
    READ (fit_ctrl_unit, *) 
    READ (fit_ctrl_unit, *) scnwrt, calwrt, calscn, rtmdbg

    !----------------------------------------------------------------------
    !xliu, 01/03/2007, read options to degrade spectral resolution
    !------------------------------------------------------------------------
    REWIND ( fit_ctrl_unit )
    CALL skip_to_filemark ( fit_ctrl_unit, lm_reduceres, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_bandselect, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF
    READ (fit_ctrl_unit, *) reduce_resolution
    READ (fit_ctrl_unit, *) reduce_slit
    READ (fit_ctrl_unit, *) redslw(1:mswath)
    IF ( reduce_slit == 1 ) redslw(1:mswath) = &
         redslw(1:mswath) / 1.66511  ! convert from FWHM to hw1e
    READ (fit_ctrl_unit, *) use_redfixwav
    IF (.NOT. reduce_resolution) use_redfixwav = .FALSE.
    READ (fit_ctrl_unit, '(A)') redfixwav_fname
    IF (use_redfixwav) THEN
      OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(redfixwav_fname)), &
           STATUS='OLD', IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
             TRIM(ADJUSTL(redfixwav_fname)), modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      ELSE
        READ(l1l2inp_unit, *) nredfixwav
        IF (nredfixwav > max_fit_pts) THEN
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
        READ(l1l2inp_unit, *) (redfixwav(i), i = 1, nredfixwav)
        CLOSE(UNIT=l1l2inp_unit)
      END IF
    ENDIF
    READ (fit_ctrl_unit, *) redsampr
    READ (fit_ctrl_unit, *) redlam

    ! -----------------------------------------------------
    ! Position cursor to read OMI channels used for fitting
    ! -----------------------------------------------------
    REWIND ( fit_ctrl_unit )
    CALL skip_to_filemark ( fit_ctrl_unit, lm_bandselect, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_bandselect, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF
    READ (fit_ctrl_unit, *) do_xbin, nxbin
    READ (fit_ctrl_unit, *) do_ybin, nybin
    READ (fit_ctrl_unit, *) rm_mgline
    READ (fit_ctrl_unit, *) numwin, do_bandavg, wcal_bef_coadd !, winwav_min, winwav_max

    IF (nxbin == 1) do_xbin = .FALSE.
    IF (nybin == 1) do_ybin = .FALSE.
    IF (.not.do_xbin) nxbin = 1
    IF (.not.do_ybin) nybin = 1

    IF (reduce_resolution) THEN
      do_bandavg = .FALSE.
      rm_mgline = .FALSE.
    ENDIF
    IF (numwin > maxwin .OR. numwin < 1) THEN
      WRITE(www_lun, *) 'Number of windows exceeds maxwin or less than 1!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF

    retlbnd = 1000.0
    retubnd = 0.0
    nviswin = 0
    DO i = 1, numwin
      READ(fit_ctrl_unit, *) band_selectors(i), winlim(i, 1), winlim(i, 2), &
                            n_band_avg(i), n_band_samp(i)

      IF ((band_selectors(i) < 0) .OR. (band_selectors(i) > mswath)) THEN
        WRITE(www_lun, *) 'No such bands exist !!!'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF

      IF (numwin > 1) THEN
        IF (winlim(i, 1) < lower_wvls(band_selectors(i)) .OR. &
            winlim(i, 2) > upper_wvls(band_selectors(i))) THEN
          WRITE(www_lun, *) 'Specified fitting windows does not make sense!!!'
          WRITE(www_lun, *) ' winlim(i,:)', winlim(i,:)
          WRITE(www_lun, *) ' lower/upper wvls(i,:)',lower_wvls(band_selectors(i)), upper_wvls(band_selectors(i))
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
      ELSE
        ! Allow 2 extra nm for radiance calibration
        IF (winlim(i, 1) < lower_wvls(band_selectors(i)) - 1.0 .OR. &
            winlim(i, 2) > upper_wvls(band_selectors(i)) + 1.0) THEN
          WRITE(www_lun, *) 'Specified fitting windows does not make sense!!!'
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
      ENDIF

      IF (i > 1) THEN
        IF  (band_selectors(i) < band_selectors(i-1) .OR. &
             winlim(i, 1) < winlim(i-1, 2))  THEN
          WRITE(www_lun, *) &
          'Incorrect band selection (must be in increasing wavelength) !!'
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
      ENDIF

      IF (winlim(i, 1) < retlbnd(band_selectors(i))) &
           retlbnd(band_selectors(i)) = winlim(i, 1)
      IF (winlim(i, 2) > retubnd(band_selectors(i))) &
           retubnd(band_selectors(i)) = winlim(i, 2)
      IF (instrument_idx == tempo_idx .or. instrument_idx == gome2_idx) THEN 
         IF (winlim(i,1) > 400.0) nviswin = nviswin + 1
      ENDIF
    END DO
    IF (MAXVAL(n_band_avg(1:numwin)) == 1) do_bandavg = .FALSE.

    DO i = 1, mswath
      IF (retlbnd(i) == 1000.0) retlbnd(i) = lower_wvls(i)
      IF (retubnd(i) == 0.0)    retubnd(i) = upper_wvls(i)
    ENDDO


    IF (instrument_idx == omi_idx) THEN 
      IF (ANY(band_selectors(1:numwin) == 1) .AND. &
           ANY(band_selectors(1:numwin) == 2)) THEN
        coadd_uv2 = .TRUE.
        nswath = 2
      ELSE
        coadd_uv2 = .FALSE.
        nswath = 1
      ENDIF
      IF (nswath == 1 .AND. band_selectors(1) == 1) THEN
        nswath = 2 ! have to read measurements around 370 nm
        coadd_uv2 = .TRUE.
      ENDIF
    ELSE IF (instrument_idx == tempo_idx) THEN     
      IF (ANY(band_selectors(1:numwin) == 1) .AND. &
         ANY(band_selectors(1:numwin) == 2)) THEN
        nswath = 2
      ELSE
        nswath = 1
      ENDIF        
    ENDIF
  
    IF (coadd_uv2) THEN
      ncoadd = 2
    ELSE
      ncoadd = 1
    ENDIF

    IF (do_bandavg) THEN
      IF (ANY(n_band_avg(1:numwin) < 1)) THEN
        WRITE(www_lun, *) 'Number of points for averaging must >= 1!!!'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF

      IF (ANY(n_band_samp(1:numwin) < 1)) THEN
        WRITE(www_lun, *) 'Number of points for sampling must >= 1!!!'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF
    ENDIF
    winwav_min = winlim(1, 1) - 5.0
    winwav_max = winlim(numwin, 2) + 5.0

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Position cursor to read general input parameters
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    REWIND (fit_ctrl_unit)
    CALL skip_to_filemark ( fit_ctrl_unit, genline_str, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, genline_str, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    READ (fit_ctrl_unit, *) use_backup, use_solcomp, avg_solcomp, avgsol_allorb
    READ (fit_ctrl_unit, *) which_slit
    IF (reduce_resolution .AND. which_slit /= 0 .AND. which_slit /= 3) THEN
      WRITE(www_lun, *) 'Have to use consistent slit function!!!'
      IF (reduce_slit == 1) which_slit = 0
      IF (reduce_slit == 2) which_slit = 3
    ENDIF

    READ (fit_ctrl_unit, *) slit_trunc_limit
    READ (fit_ctrl_unit, *) yn_varyslit, wavcal, wavcal_sol, smooth_slit, &
         slit_rad
    IF (reduce_resolution .AND. yn_varyslit) THEN
      yn_varyslit = .FALSE.
      WRITE(www_lun, *) &
           'Could not use variable slit function when to reduce resolution!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF
    !IF (.NOT. wavcal) wavcal_sol = .FALSE.
    READ (fit_ctrl_unit, *) slit_fit_pts, n_slit_step, slit_redo
    READ (fit_ctrl_unit, *) wavcal_fit_pts, n_wavcal_step, wavcal_redo
    READ (fit_ctrl_unit, *) yn_smooth
    READ (fit_ctrl_unit, *) yn_doas
    READ (fit_ctrl_unit, *) use_meas_sig, correct_merr, xbin_decerr, ybin_decerr
    READ (fit_ctrl_unit, *) linenum_lim
    READ (fit_ctrl_unit, *) pixnum_lim
    READ (fit_ctrl_unit, *) tol
    READ (fit_ctrl_unit, *) epsrel
    READ (fit_ctrl_unit, *) epsabs
    READ (fit_ctrl_unit, *) epsx

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Position cursor to read solar calibration input parameters
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    REWIND (fit_ctrl_unit)
    CALL skip_to_filemark ( fit_ctrl_unit, socline_str, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, socline_str, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    ! ---------------------------------------------------
    ! First thing to read is the Solar Reference Spectrum
    ! ---------------------------------------------------
    READ (fit_ctrl_unit, '(A)') refspec_fname(solar_idx)
    refspec_fname(solar_idx) = &
         TRIM(ADJUSTL(refdbdir)) // refspec_fname(solar_idx)
    READ (fit_ctrl_unit, *) weight_sun
    READ (fit_ctrl_unit, *) max_itnum_sol

    n_fitvar_sol = 0
    fitvar_sol_init = 0.0
    mask_fitvar_sol = 0 
    rmask_fitvar_sol = 0
    solpars: DO i = 1, max_calfit_idx

      READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
      ! ---------------------------------------------------------
      ! Check for consitency of bounds and adjust where necessary
      ! ---------------------------------------------------------
      IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp
        uptmp = vartmp
      END IF
      IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp
        lotmp = vartmp
      END IF

      IF ( idxchar == eoi3str ) EXIT solpars

      CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
      IF (which_slit == 4) THEN  
         IF ((sidx == spk_idx .or. sidx == hwe_idx) .and. lotmp == uptmp ) THEN 
           WRITE(www_lun, '(a)') 'solar cali. pars of check main_control.f90'
           WRITE(www_lun, '(a)') 'k is fixed to 2 so super should be same as gauss' 
         ELSE IF (sidx == asy_idx .or. (sidx >= vgl_idx .and. sidx <= hwr_idx)) THEN 
           vartmp = 0.0 ; lotmp = 0.0 ; uptmp = 0.0  
         ENDIF
      ELSE
        IF (sidx == spk_idx) THEN 
           vartmp = 0.0 ; lotmp = 0.0 ; uptmp = 0.0  
        ENDIF
      ENDIF

      IF (which_slit == 0 ) THEN 
        IF (sidx == hwe_idx .and. lotmp == uptmp) THEN
            WRITE(www_lun, '(a)') 'sys gaussian does not have fit. variables w.r.t hwe'
            pge_error_status = pge_errstat_error; RETURN
         ELSE IF (sidx == asy_idx .or. (sidx >= vgl_idx .and. sidx <= hwr_idx)) THEN 
            vartmp = 0.0 ; lotmp = 0.0 ; uptmp = 0.0  
         ENDIF
      ELSE IF (which_slit == 1) THEN 
         IF ( (sidx == hwe_idx .or. sidx == asy_idx) .and. lotmp == uptmp ) THEN
           WRITE(www_lun, '(a)') 'aym gaussian does not have fit. variables w.r.t hwe or asy'
           pge_error_status = pge_errstat_error; RETURN
         ELSE IF (sidx >=vgl_idx .and. sidx <= hwr_idx) THEN
               vartmp = 0.0 ; lotmp = 0.0 ; uptmp = 0.0
         ENDIF
      ENDIF

      IF ( sidx > 0 ) THEN
        fitvar_sol_init(sidx) = vartmp
        lo_sunbnd(sidx) = lotmp
        up_sunbnd(sidx) = uptmp
        IF ( lotmp < uptmp ) THEN
          n_fitvar_sol = n_fitvar_sol + 1
          mask_fitvar_sol(n_fitvar_sol) = sidx
          rmask_fitvar_sol(sidx) = n_fitvar_sol
        ENDIF
      END IF
    END DO solpars
    fitvar_sol_saved = fitvar_sol_init
    lo_sunbnd_init = lo_sunbnd
    up_sunbnd_init = up_sunbnd

    IF (ANY(rmask_fitvar_sol(wr0_idx:wr7_idx) > 0) .AND. &
          (rmask_fitvar_sol(shi_idx) > 0 .OR. rmask_fitvar_sol(squ_idx) > 0 )) THEN
          WRITE(www_lun, '(A)') 'Wavelength shi/squ should not be used with wavelength registraion!!!'
          pge_error_status = pge_errstat_error; RETURN
    ENDIF

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Position cursor to read radiance calibration input parameters
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    REWIND (fit_ctrl_unit)
    CALL skip_to_filemark ( fit_ctrl_unit, racline_str, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, racline_str, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF
    READ (fit_ctrl_unit, *) renorm
    READ (fit_ctrl_unit, *) weight_rad
    READ (fit_ctrl_unit, *) max_itnum_rad
    READ (fit_ctrl_unit, *) radwavcal_freq
    READ (fit_ctrl_unit, *) szamax
    READ (fit_ctrl_unit, *) zatmos
    READ (fit_ctrl_unit, *) phase

    fitvar_rad_init = 0.0
    fitvar_rad_str = '      '
    radpars: DO i = 1, max_calfit_idx

      READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
      ! ---------------------------------------------------------
      ! Check for consitency of bounds and adjust where necessary
      ! ---------------------------------------------------------
      IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp
        uptmp = vartmp
      END IF
      IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp
        lotmp = vartmp
      END IF
      IF (lotmp /= uptmp ) THEN
         WRITE(*,*) 'Need to re-setup radiance calibration parameters'
         STOP
      ENDIF
      IF ( idxchar == eoi3str ) EXIT radpars
      CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
      IF ( sidx > 0 ) THEN
        fitvar_rad_init(sidx) = vartmp
        fitvar_rad_str (sidx) = TRIM(ADJUSTL(idxchar))
        lo_radbnd(sidx) = lotmp
        up_radbnd(sidx) = uptmp
      END IF
    END DO radpars

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Position cursor to read radiance fitting input parameters
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    REWIND (fit_ctrl_unit)
    CALL skip_to_filemark ( fit_ctrl_unit, rafline_str, tmpchar, &
         file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, rafline_str, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    END IF

    ! By default we set the undersampling spectrum to FALSE. Only if we select
    ! it to be included in the fitting does it become trues. This way we save
    ! computation time for cases where we don't include the undersampling.
    have_undersampling = .FALSE.

    ! -------------------------------------------------------------
    ! Now keep reading spectrum blocks until EOF. This, obviously,
    ! has to be the last READ action performed from the input file.
    ! -------------------------------------------------------------
    comvidx = 0
    cm1vidx = 0
    cm2vidx = 0
    cm3vidx = 0
    comfidx = 0
    cm1fidx = 0
    cm2fidx = 0
    cm3fidx = 0
    fitvar_rad_unit = 'NoUnits'

    getpars: DO j = 1, max_rs_idx
      ! Read the spectrum identification string (SIS)
      READ (UNIT=fit_ctrl_unit, FMT='(A)', IOSTAT=errstat) tmpchar
      IF ( errstat /= file_read_ok ) THEN
        errstat = OMI_SMF_setmsg ( omsao_e_read_fitctrl_file, &
             'radiance fitting parameters', modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF
      CALL check_for_endofinput ( TRIM(ADJUSTL(tmpchar)), yn_eoi )
      IF ( yn_eoi ) EXIT getpars

      ! Convert SIS to index
      CALL string2index ( refspec_strings, max_rs_idx, tmpchar, ridx )

      ! GOME specific: the first line in the "spectrum parameter block" is
      ! the name of the corresponding reference spectrum
      READ (UNIT=fit_ctrl_unit, FMT='(A)', IOSTAT=errstat) fname

      ! Read the block of fitting parameters for current reference spectrum
      DO k = 1, mxs_idx
        READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
        ! ---------------------------------------------------------
        ! Check for consitency of bounds and adjust where necessary
        ! ---------------------------------------------------------
        IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
          lotmp = vartmp
          uptmp = vartmp
        END IF
        IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
          uptmp = vartmp
          lotmp = vartmp
        END IF

        CALL string2index ( radfit_strings, mxs_idx, idxchar, sidx )
        IF ( sidx > 0 ) THEN
          i = max_calfit_idx + (ridx-1)*mxs_idx + sidx
          fitvar_rad_init(i) = vartmp
          fitvar_rad_str (i) = TRIM(ADJUSTL(tmpchar))
          lo_radbnd (i) = lotmp
          up_radbnd (i) = uptmp
          IF ( (ridx == us1_idx .OR. ridx == us2_idx) .AND. &
               ANY ( (/ vartmp,lotmp,uptmp /) /= 0.0 ) ) &
               have_undersampling = .TRUE.

          IF ( ridx == com_idx .AND. lotmp < uptmp ) THEN
             comvidx = i
          ENDIF

          IF ( ridx == com1_idx .AND. lotmp < uptmp ) THEN
             cm1vidx = i
          ENDIF

          IF ( ridx == com2_idx .AND. lotmp < uptmp ) THEN
             cm2vidx = i
          ENDIF

          IF ( ridx == com3_idx .AND. lotmp < uptmp ) THEN
             cm3vidx = i
          ENDIF

          IF (lotmp < uptmp) THEN
           refspec_fname(ridx) = TRIM(ADJUSTL(refdbdir)) // TRIM(ADJUSTL(fname))
          ENDIF
        END IF
      END DO
        
      ! read shift parameter
      READ (fit_ctrl_unit, *) idxchar1, vartmp, lotmp, uptmp
      IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp
        uptmp = vartmp
      END IF

      IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp
        lotmp = vartmp
      END IF

      i =  max_calfit_idx + (ridx-1)*mxs_idx + 1
      IF  (ALL(lo_radbnd(i:i+2) - up_radbnd(i:i+2) >= 0.0)) THEN
        vartmp = 0.0
        uptmp = 0.0
        lotmp = 0.0
      END IF

      i =  shift_offset + ridx
      fitvar_rad_init(i) = vartmp
      fitvar_rad_str(i) = TRIM(ADJUSTL(idxchar1))
      fitvar_rad_unit(i) = 'nm'

      lo_radbnd (i) = lotmp
      up_radbnd (i) = uptmp
    END DO getpars


    ! -------------------------------------------------------------
    ! Find the indices of those variables that are actually varied
    ! during the fitting, and save those in MASK_FITVAR_RAD. Save
    ! the number of varied parameters in N_FITVAR_RAD.
    !
    ! In addition, we need to determine the fitting indices that
    ! will make up the final fitted column of the molecule(s) in
    ! question. For this we have to jump through a double loop:
    ! Since we are compressing the fitting parameter array to
    ! include only the varied parameters, the final covariance
    ! matrix, which is crucial for determining the uncertainties,
    ! only knows the compressed indices. Therefore we have to have
    ! an "index of an index" type array, that remembers the index
    ! position of indices of the final molecule(s), AS THEY APPEAR
    ! IN THE COMPRESSED FITTING PARAMETER LIST.
    !
    ! For the latter task it is easier to split the loops into
    ! calibration parameters (unrelated to the final fitted column)
    ! and reference spectra parameters.
    ! -------------------------------------------------------------
    n_fitvar_rad = 0
    mask_fitvar_rad = 0
    rmask_fitvar_rad = 0
    database_indices = 0
    ! --------------------------------
    ! First the calibration parameters
    ! --------------------------------
    DO i = 1, max_calfit_idx
      IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad
      END IF
    END DO

    ! ------------------------------------
    ! Now the reference spectra parameters
    ! ------------------------------------
    n_fincol_idx = 0
    DO i = 1, max_rs_idx
      idx = max_calfit_idx + (i-1) * mxs_idx
      DO j = mns_idx, mxs_idx
        ! ----------------------------------------------
        ! Assign only entries that are varied in the fit
        ! ----------------------------------------------
        IF ( lo_radbnd(idx+j) < up_radbnd(idx+j) ) THEN
          n_fitvar_rad = n_fitvar_rad + 1
          mask_fitvar_rad(n_fitvar_rad) = idx+j
          rmask_fitvar_rad(idx+j) = n_fitvar_rad
          database_indices(n_fitvar_rad) = i

          ! -------------------------------------------------
          ! And here the loop over the final column molecules.
          ! We have to match the FITCOL_IDX with the current
          ! molecule index, and then remember the position of
          ! the fitting index in the MASK_FITVAR_RAD array.
          ! The second index remembers the reference spectrum
          ! that is associated with this molecule, so that we
          ! can easily access its normalization factor.
          ! -------------------------------------------------
          !IF (.NOT. ozprof_flag) THEN  !xliu
          !  getfincol: DO k = 1, n_mol_fit
          !    IF ( fitcol_idx(k) == i ) THEN
          !      n_fincol_idx = n_fincol_idx + 1
          !      fincol_idx (1,n_fincol_idx) = n_fitvar_rad
          !      fincol_idx (2,n_fincol_idx) = i
          !      EXIT getfincol
          !    END IF
          !  END DO getfincol
          !ENDIF   !xliu

          IF (idx + j == comvidx) comfidx = n_fitvar_rad
          IF (idx + j == cm1vidx) cm1fidx = n_fitvar_rad
          IF (idx + j == cm2vidx) cm2fidx = n_fitvar_rad
          IF (idx + j == cm3vidx) cm3fidx = n_fitvar_rad

        END IF
      END DO
    END DO

    ! For shift
    ntsh = 0
    DO i = 1, max_rs_idx
      j = shift_offset + i
      IF ( lo_radbnd(j) < up_radbnd(j) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = j
        rmask_fitvar_rad(j) = n_fitvar_rad
        ntsh = ntsh + 1
      END IF
    END DO


    !  ! --------------------------------------
    !  ! Position cursor to read AMF table file
    !  ! --------------------------------------
    !  REWIND ( fit_ctrl_unit )
    have_amftable = .FALSE.
    !  CALL skip_to_filemark ( fit_ctrl_unit, lm_amftable, tmpchar, file_read_stat )
    !  IF ( file_read_stat /= file_read_ok ) THEN
    !     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_amftable, modulename, 0)
    !     pge_error_status = pge_errstat_warning
    !  ELSE
    !     READ (fit_ctrl_unit, *) have_amftable
    !     IF ( have_amftable ) THEN
    !        READ (fit_ctrl_unit, '(A)') static_input_fnames(amf_idx)
    !        static_input_fnames(amf_idx) = TRIM(ADJUSTL(refdbdir)) // static_input_fnames(amf_idx)
    !     ENDIF
    !  END IF

    !  !xliu: add the following block
    !  ! ----------------------------------------------------------
    !  ! Position cursor to read whether to retrieve ozone profile
    !  ! ----------------------------------------------------------
    !  REWIND ( fit_ctrl_unit )
    !  CALL skip_to_filemark ( fit_ctrl_unit, ozprof_str, tmpchar, file_read_stat )
    !  IF ( file_read_stat /= file_read_ok ) THEN
    !     ozprof_flag = .FALSE.
    !     WRITE(*, *) 'This algorithm is only for ozone profile retrieval!!!'
    !     pge_error_status = pge_errstat_error; RETURN
    !  ELSE
    !     READ (fit_ctrl_unit, *) ozprof_flag
    !     IF (ozprof_flag)  THEN
    !!       READ (fit_ctrl_unit, '(A)') ozprof_input_fname !! commented Kai
    !     END IF
    !  END IF
    ozprof_flag = .TRUE.

    version = 1
    errstat = PGS_PC_getreference( OZPROF_CTRL_LUN, version, &
         ozprof_input_fname )
    IF( errstat /= PGS_S_SUCCESS ) THEN
      WRITE(msg, '(A,I10,I4)') 'get file from lun=', OZPROF_CTRL_LUN, version
      errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
           modulename, 0)
      pge_error_status = pge_errstat_error
      RETURN
    ELSE
      errstat = OMI_SMF_setmsg (OMI_S_SUCCESS,TRIM(ozprof_input_fname), & !! Kai
           modulename, 0)
    END IF


    ! -----------------------------------------------
    ! Position cursor to read molecule name(s) to fit
    ! -----------------------------------------------
    IF (.NOT. ozprof_flag) THEN  !xliu
    !  REWIND ( fit_ctrl_unit )
    !  CALL skip_to_filemark ( fit_ctrl_unit, molline_str, tmpchar, &
    !       file_read_stat )
    !  IF ( file_read_stat /= file_read_ok ) THEN
    !    errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, molline_str, &
    !         modulename, 0)
    !    pge_error_status = pge_errstat_error
    !    RETURN
    !  END IF
    !  READ (fit_ctrl_unit, '(A)') tmpchar
    !  CALL get_mols_for_fitting ( tmpchar, n_mol_fit, fitcol_idx, errstat )
    !  IF ( errstat /= pge_errstat_ok ) THEN
    !    errstat = OMI_SMF_setmsg (omsao_e_get_molfitname, '', modulename, 0)
    !    pge_error_status = pge_errstat_error
    !    RETURN
    !  END IF
    ENDIF

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Position cursor to read Level 1 input files
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    select_lonlat = .FALSE.
    rw_l1l2_pcf   = .TRUE.        ! read l1 (write l2) fname from PCF file
    pixlim = -9999
    linelim = -9999


    l1l2_files = TRIM(ADJUSTL(tabdir)) // '../control/L1L2_fnames.inp'

    IF (.NOT. rw_l1l2_pcf ) THEN
      OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(l1l2_files)), STATUS='OLD', &
           IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
             TRIM(ADJUSTL(l1l2_files)), modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      ELSE
        READ(l1l2inp_unit, '(A)') l1b_irrad_filename
        READ(l1l2inp_unit, '(A)') l1b_rad_filename
        READ(l1l2inp_unit, '(A)') l2_cld_filename
        READ(l1l2inp_unit, '(A)') l2_filename
        READ(l1l2inp_unit, *)     linelim, pixlim
        READ(l1l2inp_unit, *, IOSTAT=errstat)  select_lonlat
        IF ( errstat /= pge_errstat_ok ) select_lonlat = .FALSE.
        IF (select_lonlat) THEN
          READ(l1l2inp_unit, *)   slat, elat, slon, elon
          IF (slat >= elat .OR. slon >= elon) THEN
            WRITE(www_lun, *) 'Incorrect lat/lon range!!!'
            pge_error_status = pge_errstat_error
            RETURN
          ENDIF
          pixlim(1)  = -5
          pixlim(2) = -5
          linelim(1) = -5
          linelim(2) = -5
        ENDIF
        CLOSE(UNIT=l1l2inp_unit)
      END IF
    ELSE !! Start Read L1L2 file from PCF 
      ! read l1b_irrad_filename
      version = 1
      errstat = PGS_PC_getreference( L1B_IRR_FILE_LUN, version, &
           l1b_irrad_filename )
      IF( errstat /= PGS_S_SUCCESS ) THEN
        WRITE(msg, '(A,I10,I4)') 'get file from lun=', L1B_IRR_FILE_LUN, version             
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
             modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF
      ! read l1b_rad_filename
      version = 1
      errstat = PGS_PC_getreference( L1B_UV_FILE_LUN, version, &
           l1b_rad_filename )
      IF( errstat /= PGS_S_SUCCESS ) THEN
        WRITE(msg, '(A,I10,I4)') 'get file from lun=', L1B_UV_FILE_LUN, version
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
             modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      ELSE
        errstat = OMI_SMF_setmsg( OMI_S_SUCCESS, 'l1b_rad_filename ='// &
             TRIM(l1b_rad_filename), modulename, 0)
      END IF
      ! read l2_cld_filename
      version = 1
      errstat = PGS_PC_getreference( L2_CLD_FILE_LUN, version, l2_cld_filename )
      IF( errstat /= PGS_S_SUCCESS ) THEN
        WRITE(msg, '(A,I10,I4)') 'get file from lun=', L2_CLD_FILE_LUN, version
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
             modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      ELSE
        errstat = OMI_SMF_setmsg(OMI_S_SUCCESS,  &
             'l2_cld_filename ='//TRIM(l2_cld_filename), modulename, 0)
      END IF
      ! read l2_filename
      version = 1
      errstat = PGS_PC_getreference( L2_OUT_LUN, version, l2_filename )
      IF( errstat /= PGS_S_SUCCESS ) THEN
        WRITE(msg, '(A,I10,I4)') 'get file from lun=', L2_OUT_LUN, version
        errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, msg, &
             modulename, 0)
        pge_error_status = pge_errstat_error
        RETURN
      ELSE
        errstat = OMI_SMF_setmsg(OMI_S_SUCCESS,  &
             'l2_filename ='//TRIM(l2_filename), modulename, 0)
      END IF

      !! get run specific inputs
      errstat = PGS_PC_GetConfigData( LINE_SAMPLE_RANGE_LUN, msg )
      IF( errstat /= PGS_S_SUCCESS ) THEN
        WRITE( msg,'(A,I8)' ) "get line sample range failed at LUN = ", &
             LINE_SAMPLE_RANGE_LUN
        errstat = OMI_SMF_setmsg( errstat, msg, modulename, 0 )
        pge_error_status = pge_errstat_error
        RETURN
      ELSE
        READ( msg, *) linelim, pixlim
        errstat = OMI_SMF_setmsg( OMI_S_SUCCESS, &
             'PCF: UV2 Line sample ranges='//TRIM(msg), modulename, 0 )
      ENDIF
    END IF

    IF( scnwrt ) THEN
        WRITE(*, '(A)') TRIM(l1b_irrad_filename)
        WRITE(*, '(A)') TRIM(l1b_rad_filename)
        WRITE(*, '(A)') TRIM(l2_cld_filename)
    ENDIF

    IF (instrument_idx == omi_idx) THEN
       ! obtain orbit number from irradiance file
       i = INDEX(l1b_irrad_filename, '-o') + 2
       orbcsol = l1b_irrad_filename(i : i + 5)
       READ (orbcsol, *) orbnumsol
       READ (l1b_irrad_filename(i+7 : i + 9), *) omisol_version

       ! obtain orbit number from radiance file
       i = INDEX(l1b_rad_filename, '-o') + 2
       orbc = l1b_rad_filename(i : i + 5)
       READ (orbc, *) orbnum

       ! check cloud file
       i = INDEX(l2_cld_filename, '-o') + 2
       cldorbc = l2_cld_filename(i : i + 5)
       READ (cldorbc, *) cldorb

       IF (cldorb /= orbnum) THEN
          WRITE(www_lun, *) &
               'Inconsistent orbit number between radiance and cloud file!!!'
          pge_error_status = pge_errstat_error
          RETURN
       ENDIF
      ! generate identifer for irradiance and radiance spectrum
      i = INDEX(l2_filename, 'OMIO3PROF')
      outdir = l2_filename(1:i-1)
      rad_identifier = 'o' // orbc
      sol_identifier = 'o' // orbcsol
    ELSE IF (instrument_idx == tempo_idx) THEN 
      orbc = '00000'
      orbcsol = orbc
      !i = INDEX(l2_filename, 'OMIO3PROF')
      outdir = './' !l2_filename(1:i-1)
      rad_identifier = 'o' // orbc
      sol_identifier = 'o' // orbcsol
    ELSE IF (instrument_idx == gome_idx) THEN
      i = INDEX(l1b_irrad_filename, 'lv1_') + 11
      orbc = l1b_rad_filename(i-2:i)
      orbcsol = orbc
      i = INDEX(l1b_irrad_filename, 'lv1') + 4
      sol_identifier = l1b_irrad_filename(i:i+7)
      i = INDEX(l1b_rad_filename, 'lv1') + 4
      rad_identifier = l1b_rad_filename(i:i+7)
      j = INDEX(l2_filename, 'lv2')
      outdir = l2_filename(1:j-1)
        !ELSEIF (instrument_idx == scia_idx) THEN
        !   i = INDEX(l1b_rad_filename, 'Ch1orb') + 28
        !   orbc = l1b_rad_filename(i-2:i)
        !   i = INDEX(l1b_irrad_filename, 'Ch1orb') + 20
        !   sol_identifier = l1b_irrad_filename(i:i+4) // l1b_irrad_filename(i+6:i+8)
        !   rad_identifier = sol_identifier
        !   j = INDEX(l2_filename, 'lv2')        ; outdir = l2_filename(1:j-1)
        !   sciaorb_identifier =  l1b_irrad_filename(i-3:i+11)
    ELSE IF (instrument_idx == gome2_idx) THEN
      orbc = '0'
      i = INDEX(l1b_irrad_filename, 'GOME_xxx_1B') + 18
      sol_identifier = l1b_irrad_filename(i:i+7)
      i = INDEX(l1b_rad_filename, 'GOME_xxx_1B') + 18
      rad_identifier = l1b_rad_filename(i:i+7)
      j = INDEX(l2_filename, 'lv2')
      outdir = l2_filename(1:j-1)
    ENDIF
    IF (use_backup) THEN 
      sol_identifier = rad_identifier 
    ENDIF
    sn = slit_name(which_slit)
    slit_fname    = TRIM(ADJUSTL(outdir))// 'slit_'//sol_identifier//'.'//sn
    swavcal_fname = TRIM(ADJUSTL(outdir))// 'swavcal_'//sol_identifier// '.'//sn
    rslit_fname   = TRIM(ADJUSTL(outdir))// 'rslit_'// rad_identifier // '.'//sn
    wavcal_fname  = TRIM(ADJUSTL(outdir))//'wavcal_'// rad_identifier //'.'//sn

    ! ----------------------------------------------
    ! Position cursor to read HDF output flags (CRN)
    ! ----------------------------------------------
    REWIND ( fit_ctrl_unit )
    CALL skip_to_filemark ( fit_ctrl_unit, lm_l2hdf, tmpchar, file_read_stat )
    IF ( file_read_stat /= file_read_ok ) THEN
      errstat = OMI_SMF_setmsg (omsao_w_read_fitctrl_file, lm_l2hdf, &
           modulename, 0)
      pge_error_status = pge_errstat_warning
      l2_hdf_flag = 0  ! Write ASCII
    ELSE
      READ (fit_ctrl_unit, '(I8)') l2_hdf_flag

      ! Check output data format
      IF (instrument_idx == gome2_idx .AND. (l2_hdf_flag < 1 .OR. &
           l2_hdf_flag > 2)) THEN
        WRITE(www_lun, *) 'GOME-2 data can only be written in HDF'
        pge_error_status = pge_errstat_error
        RETURN
      ELSE IF ((instrument_idx == gome_idx .OR. &
           instrument_idx == scia_idx) .AND. l2_hdf_flag  > 0) THEN
        WRITE(www_lun, *) 'GOME-1/SCIA data can only be written in ASCII'
        pge_error_status = pge_errstat_error
        RETURN
      ELSE IF (instrument_idx == omi_idx .AND. (l2_hdf_flag /= 0 .AND. &
           l2_hdf_flag /= 3 .AND. l2_hdf_flag /=4)) THEN
        WRITE(www_lun, *) 'OMI data can only be written in ASCII or HDF-EOS5 or NC'
        pge_error_status = pge_errstat_error
        RETURN
      ELSE IF (instrument_idx == tempo_idx .AND. (l2_hdf_flag /= 0 .AND. &
           l2_hdf_flag /= 3 .AND. l2_hdf_flag /=4)) THEN
        WRITE(www_lun, *) 'TEMPO data can only be written in ASCII or HDF-EOS5 or NC'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF
    END IF


    ! ------------------------------------------------
    ! Check for consistency of pixel limits to process
    ! ------------------------------------------------
    IF (select_lonlat) THEN
      CALL find_scan_line_range(slat, elat, slon, elon, linelim(1), &
           linelim(2), pixlim(1), pixlim(2), pge_error_status )
      IF (pixlim(1) < 0 .OR. linelim(1) < 0 .OR. pge_error_status >= &
           pge_errstat_error) THEN
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF
      ! pixlim is based on UV1, if both channels are selected
      IF (coadd_uv2) THEN
        pixlim(1) = pixlim(1) * ncoadd - 1
        pixlim(2) = pixlim(2) * ncoadd
      ENDIF
    ENDIF

    ! check for boundaries
    IF (linelim(1) /= -9999) linenum_lim =  linelim
    IF (pixlim(1)  /= -9999) pixnum_lim  =  pixlim
    IF ( ALL ( linenum_lim < 0 ) )      linenum_lim(1:2) = (/ 1, ntimes_max /)
    IF ( linenum_lim(1) > linenum_lim(2) ) linenum_lim([1, 2]) = &
         linenum_lim([2, 1])
    IF ( linenum_lim(1) < 1 )              linenum_lim(1) = 1
    IF ( linenum_lim(2) > ntimes_max )     linenum_lim(2) = ntimes_max

    IF ( ALL ( pixnum_lim < 0 ) )       pixnum_lim(1:2) = (/ 1, nxtrack_max /)
    IF ( pixnum_lim(1) > pixnum_lim(2) )   pixnum_lim([1, 2]) = &
         pixnum_lim([2, 1])
    IF ( pixnum_lim(1) < 1 )               pixnum_lim(1) = 1
    IF ( pixnum_lim(2) > nxtrack_max )     pixnum_lim(2) = nxtrack_max

    ! check for selected across track position (must start from odd positions)
    IF (instrument_idx == omi_idx) THEN    
     IF (coadd_uv2)  THEN
      i = pixnum_lim(2)-pixnum_lim(1) + 1
      IF ( MOD(pixnum_lim(1), ncoadd) /= 1 .OR. MOD(i, ncoadd) /= 0 ) THEN
        WRITE(www_lun, '(A,2I4)') &
             'Incorrect across track positions to be coadded: ', &
             pixnum_lim(1:2)
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF
      ! NINT behaviour is not consistent between Intel and GNU builds
      ! when input is exactly integer+0.5
      ! pixnum_lim = NINT(1.0 * pixnum_lim / ncoadd)
      pixnum_lim = INT((1.0 * pixnum_lim / ncoadd)+0.5)
     ENDIF
     ! must divide and must start from odd coadded positions
     IF (do_xbin .AND. nxbin > 1) THEN
       i = pixnum_lim(2)-pixnum_lim(1) + 1
       IF ( MOD (i, nxbin) /= 0 .OR. MOD(pixnum_lim(1), nxbin) /= 1 ) THEN
         WRITE(www_lun, '(A,2I4)') &
             'Incorrect across track binning option: ', pixnum_lim(1:2)
              pge_error_status = pge_errstat_error
         RETURN
       ENDIF
     ELSE
       nxbin = 1
     ENDIF

     IF( coadd_uv2 ) THEN
       WRITE(msg, '(A,2I5,2I3)') 'Processing UV1 Line sample ranges=', &
           linenum_lim(1:2),pixnum_lim(1:2)
     ELSE
       WRITE(msg, '(A,2I5,2I3)') 'Processing UV2 Line sample ranges=', &
           linenum_lim(1:2),pixnum_lim(1:2)
     ENDIF
     errstat = OMI_SMF_setmsg( OMI_S_SUCCESS, TRIM(msg), modulename, 0 )
    ELSE IF (instrument_idx == tempo_idx) THEN 
      IF (do_xbin .and. nxbin > 1) THEN 
         IF (pixnum_lim(1) <= nybin) THEN 
             pixnum_lim(1) = 1
         ELSE
             pixnum_lim(1) = INT(( pixnum_lim(1)/nxbin)*nxbin   ) +1
         ENDIF
       
             i = pixnum_lim(2)-pixnum_lim(1) + 1
             IF (i < nxbin) THEN 
                 pixnum_lim(2) = pixnum_lim(1) + nxbin -1
             ENDIF
             IF (mod(i, nxbin) /= 0 ) THEN 
               pixnum_lim(2) = CEILING(1.0 * i / nxbin) * nxbin + pixnum_lim(1) - 1
             ENDIF     
             IF (pixnum_lim(2) > nxtrack_max) THEN
               pixnum_lim(2) = pixnum_lim(2) - nxbin 
             ENDIF       
      ENDIF
         ! print * , nint(0.4), ceiling(0.4), int(0.4)
         ! print * , nint(0.6), ceiling(0.6), int(0.6)
    ENDIF
    
    ! Could start from any positions, adjust line positions if necessary
    IF (do_ybin .AND. nybin > 1)  THEN
      IF( linenum_lim(1) <= nybin ) THEN
        linenum_lim(1) = 1
      ELSE
        linenum_lim(1) = INT((linenum_lim(1)-1)/nybin)*nybin + 1
      ENDIF
      IF ( linenum_lim(1) > linenum_lim(2) )   linenum_lim(2) = &
           linenum_lim(1)
      i = linenum_lim(2)-linenum_lim(1) + 1
      IF (MOD(i, nybin) /= 0) THEN
        linenum_lim(2) = CEILING(1.0 * i / nybin) * nybin + linenum_lim(1) - 1
        IF (linenum_lim(2) > ntimes_max) linenum_lim(2) = &
             linenum_lim(2) - nybin
      ENDIF
    ELSE
      nybin = 1
    ENDIF
    IF (instrument_idx /= tempo_idx ) THEN 
    WRITE(slinechar, '(I4.4)') linenum_lim(1)
    WRITE(elinechar, '(I4.4)') linenum_lim(2)
    WRITE(sxchar, '(I2.2)')    pixnum_lim(1)
    WRITE(exchar, '(I2.2)')    pixnum_lim(2)
    ELSE
    WRITE(slinechar, '(I3.3)') linenum_lim(1)
    WRITE(elinechar, '(I3.3)') linenum_lim(2)
    WRITE(sxchar, '(I4.4)')    pixnum_lim(1)
    WRITE(exchar, '(I4.4)')    pixnum_lim(2)
    ENDIF

    !IF( l2_hdf_flag == 0) THEN  !! Modify output file name only when it is L2 HE5 output, Kai
      IF ( .NOT. (linenum_lim(1) == 1 .AND. linenum_lim(2) == ntimes_max)) THEN
        l2_filename = TRIM(ADJUSTL(l2_filename)) &
         // '_L' // TRIM(ADJUSTL(slinechar)) // '-' //TRIM(ADJUSTL(elinechar))
      ENDIF
      IF ( .NOT. (pixnum_lim(1) == 1 .AND. (pixnum_lim(2) == nxtrack_max .OR. &
        (coadd_uv2 .AND. pixnum_lim(2) == nxtrack_max / ncoadd)))) THEN
        l2_filename = TRIM(ADJUSTL(l2_filename)) // '_X' &
        // ADJUSTL(TRIM(sxchar)) // '-' // ADJUSTL(TRIM(exchar))
      ENDIF

      IF (do_xbin) THEN
        WRITE(xbinchar, '(A2,I1)') 'BX', nxbin
        l2_filename = TRIM(ADJUSTL(l2_filename)) // '-' //ADJUSTL(TRIM(xbinchar))
      ENDIF
      IF (do_ybin) THEN
        WRITE(ybinchar, '(A2,I1)') 'BY', nybin
        l2_filename = TRIM(ADJUSTL(l2_filename)) // '-' //TRIM(ADJUSTL(ybinchar))
      ENDIF
    !ENDIF 

    IF (l2_hdf_flag == 4) THEN 
        l2_filename=TRIM(ADJUSTL(L2_filename))//'.nc'
    ELSE IF (l2_hdf_flag == 3) THEN 
        l2_filename=TRIM(ADJUSTL(L2_filename))//'.he5'
    ELSE IF (l2_hdf_flag == 0) THEN 
        l2_filename=TRIM(ADJUSTL(L2_filename))//'.out'
    ELSE IF (l2_hdf_flag == 1) THEN 
        l2_filename=TRIM(ADJUSTL(L2_filename))//'.h5'
    ELSE
          STOP 
    ENDIF
    ! -----------------------------------------------
    ! Close fitting control file, report SUCCESS read
    ! -----------------------------------------------
    CLOSE ( UNIT=fit_ctrl_unit )

    ! ------------------------------------------------------------------------
    ! Read fitting conrol parameters from input file for
    ! ozone profile variables
    IF (ozprof_flag) THEN 
      CALL read_ozprof_input ( &
           fit_ctrl_unit, ozprof_input_fname, pge_error_status )
      IF ( pge_error_status >= pge_errstat_error ) RETURN
    ENDIF
    ! refnhextra must >= 1 and radnhtrunc > refnhextra, if interpolation
    ! is performed radnhtrunc should be
    IF ( (ntsh == 0 .AND. nsh == 0 .AND. nos == 0 .AND. nsl == 0) &
         .OR. (do_simu .AND. .NOT. radcalwrt)) THEN
      radnhtrunc = 3
      refnhextra = 2
      IF (instrument_idx == tempo_idx) THEN 
        radnhtrunc = 5
        refnhextra = 2
      ENDIF
      !ELSE IF (reduce_resolution .AND. use_redfixwav) THEN
      !   radnhtrunc = 2; refnhextra = 1
    ELSE
      radnhtrunc = 3
      refnhextra = 2
      IF (instrument_idx == tempo_idx) THEN 
        radnhtrunc = 5
        refnhextra = 2
      ENDIF
    ENDIF
    IF (use_so2dtcrs) refspec_fname (so2_idx) = zerospec_string
    !IF (use_o4dtcrs) refspec_fname (o2o2_idx) = zerospec_string
    IF (use_o2dptcrs) refspec_fname (o2_idx) = zerospec_string
    IF (use_h2odptcrs) refspec_fname (h2o_idx) = zerospec_string
    ! ------------------------------------------------------------------------
    fitvar_rad_saved = fitvar_rad_init
    IF ( yn_doas ) THEN
      pm_one     = -1.D0
    ELSE
      pm_one     = 1.D0
    END IF

    errstat = OMI_SMF_setmsg(OMI_S_SUCCESS, "done reading.", modulename, 0) !! Kai

    num_param = n_fitvar_rad - nfgas - &
                (ozfit_end_index - ozfit_start_index + 1)
    IF (scnwrt) THEN
      WRITE(*, '(A, I8, A, I8)') ' n_fitvar_rad = ', n_fitvar_rad, '      nlayer =', nlay
      DO i = 1, n_fitvar_rad
         IF (i < ozfit_start_index .or. i > ozfit_end_index) THEN
          WRITE(*, '(2I5, A10,f5.2, A20)') i,mask_fitvar_rad(i),fitvar_rad_str(mask_fitvar_rad(i)),  &
          fitvar_rad_init(masK_fitvar_rad(i)),fitvar_rad_unit(mask_fitvar_rad(i))
         ENDIF
      ENDDO
    ENDIF
    RETURN
    DO i = 1, 192
      WRITE(*, '(2I5, A10,f5.2, A20)') i, rmask_fitvar_rad(i), fitvar_rad_str(i),&
      fitvar_rad_init(i),  fitvar_rad_unit(i)
    ENDDO
    !xliu, 09/23/05 Add direcotry, remove hard code directory
    ! ----------------------------------------------------------
    ! Position cursor to read database directory
    ! ----------------------------------------------------------
    !  REWIND ( fit_ctrl_unit )
    !  CALL skip_to_filemark ( fit_ctrl_unit, lm_atmdb, tmpchar, file_read_stat )
    !  IF ( file_read_stat /= file_read_ok ) THEN
    !     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_atmdb, modulename, 0)
    !     pge_error_status = pge_errstat_error; RETURN
    !  ELSE
    !!    READ (fit_ctrl_unit, '(A)') atmdbdir !! commented Kai
    !  ENDIF
    !  REWIND ( fit_ctrl_unit)
    !  CALL skip_to_filemark ( fit_ctrl_unit, lm_refdb, tmpchar, file_read_stat )
    !  IF ( file_read_stat /= file_read_ok ) THEN
    !     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, lm_refdb, modulename, 0)
    !     pge_error_status = pge_errstat_error; RETURN
    !  ELSE
    !!    READ (fit_ctrl_unit, '(A)') refdbdir !! commented Kai
    !  ENDIF
    RETURN
  END SUBROUTINE read_fitting_control_file


!  SUBROUTINE get_mols_for_fitting ( tmpchar, n_mol_fit, fitcol_idx, errstat )
!
!    USE OMSAO_indices_module,    ONLY: refspec_strings, max_rs_idx
!    USE OMSAO_parameters_module, ONLY: max_mol_fit
!    USE OMSAO_errstat_module,    ONLY: pge_errstat_error
!    use utilities, only: get_substring, string2index
!
!    IMPLICIT NONE
!
!    ! ===============
!    ! Input variables
!    ! ===============
!    CHARACTER (LEN=*), INTENT (INOUT) :: tmpchar
!
!    ! ================
!    ! Output variables
!    ! ================
!    INTEGER,                          INTENT (OUT) :: n_mol_fit, errstat
!    INTEGER, DIMENSION (max_mol_fit), INTENT (OUT) :: fitcol_idx
!
!    ! ===============
!    ! Local variables
!    ! ===============
!    INTEGER                      :: i, ncl, sstart, sidx
!    LOGICAL                      :: yn_eoc
!    CHARACTER (LEN=LEN(tmpchar)) :: tmpsub
!
!    ! ----------------------------
!    ! Initialize output quantities
!    ! ----------------------------
!    n_mol_fit = 0
!    fitcol_idx = 0
!
!    ! -----------------------------------------------
!    ! Get names and indices for main molecules to fit
!    ! -----------------------------------------------
!    sstart = 0
!    ncl = 0
!    yn_eoc = .FALSE.
!    getmolnames: DO i = 1, max_mol_fit
!      ! ---------------------------------------------------------
!      ! Extract index string, find index, then extract file name.
!      ! ---------------------------------------------------------
!      CALL get_substring ( tmpchar, sstart, tmpsub, ncl, yn_eoc )
!      IF ( ncl > 0 ) THEN
!        CALL string2index ( refspec_strings, max_rs_idx, tmpsub, sidx )
!        IF ( sidx > 0 ) THEN
!          n_mol_fit = n_mol_fit + 1
!          fitcol_idx(n_mol_fit) = sidx
!        END IF
!      END IF
!      IF ( yn_eoc ) EXIT getmolnames
!    END DO getmolnames
!    IF ( n_mol_fit == 0 .OR. ALL(fitcol_idx == 0) ) errstat = pge_errstat_error
!
!    RETURN
!  END SUBROUTINE get_mols_for_fitting


end module m_read_fitting_controls
