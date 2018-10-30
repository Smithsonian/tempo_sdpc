!
module omi_fitting_aux

  public omi_adj_solar_data, omi_adj_earthshine_data, &
       omi_set_fitting_parameters, corruv2wav, timestamp, get_doy
  private omi_adj_rad_sig, load_omi_comres

contains

  SUBROUTINE omi_adj_solar_data (pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx!, hwe_idx
    USE OMSAO_parameters_module, ONLY: normweight!, mrefl
    USE OMSAO_variables_module,  ONLY: curr_sol_spec, n_irrad_wvl, &
         use_meas_sig, numwin, nsol_ring, sol_spec_ring, nsolpix, &
         yn_varyslit, slit_rad, solwinfit, nslit, slitwav, slitfit, &
         sring_fidx, sring_lidx,  currpix, which_slit!, &
         !reduce_resolution, slitdis
    USE OMSAO_omidata_module,    ONLY: omi_nwav_irrad, omi_nsolpix, &
         omi_irradiance_wavl, omi_irradiance_spec, omi_irradiance_prec, &
         omi_nsolring, omi_solspec_ring, omi_solspecr, omi_nslit_rad, &
         omi_slitwav_rad, omi_radslitfit, omi_nslit_sol, &
         omi_slitwav_sol, omi_solslitfit, omi_solwinfit, solring_lin, &
         solring_uin, omi_solnorm!, orbnumsol
    USE ozprof_data_module,      ONLY: div_sun, sun_posr, sun_specr, nrefl
    USE OMSAO_errstat_module 

    IMPLICIT NONE

    ! =================
    ! Output variables
    ! =================
    INTEGER, INTENT (OUT) :: pge_error_status

    !INTEGER :: i

    pge_error_status = pge_errstat_ok

    ! Solar Spectrum
    n_irrad_wvl = omi_nwav_irrad(currpix) 
    div_sun     = omi_solnorm(currpix)

    curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = omi_irradiance_wavl(1:n_irrad_wvl, currpix) 
    curr_sol_spec(spc_idx, 1:n_irrad_wvl) = omi_irradiance_spec(1:n_irrad_wvl, currpix) 
    IF (use_meas_sig ) THEN
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = omi_irradiance_prec(1:n_irrad_wvl, currpix)
    ELSE
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
    ENDIF
    nsolpix(1:numwin) = omi_nsolpix(1:numwin, currpix) 

    ! Solar Spectrum for Ring Calculation
    nsol_ring = omi_nsolring(currpix)
    sol_spec_ring(1:2, 1:nsol_ring) = omi_solspec_ring(1:2, 1:nsol_ring, currpix)
    sring_fidx = solring_lin(currpix)
    sring_lidx = solring_uin(currpix)

    !WRITE(90, '(F12.6, D14.6)') (sol_spec_ring(1:2, i), i = 1, nsol_ring)

    ! Reflectance spectrum at ~370 nm
    sun_posr (1:nrefl) = omi_solspecr(wvl_idx, 1:nrefl, currpix)
    sun_specr(1:nrefl) = omi_solspecr(spc_idx, 1:nrefl, currpix)

    ! Load slit calibration parameters
    IF (which_slit < 5) THEN
      IF (yn_varyslit) THEN
        IF (slit_rad) THEN
          nslit = omi_nslit_rad(currpix)
          slitwav(1:nslit) = omi_slitwav_rad(1:nslit, currpix)
          slitfit = omi_radslitfit(:, :, :, currpix)
        ELSE
          nslit = omi_nslit_sol(currpix)
          slitwav(1:nslit) = omi_slitwav_sol(1:nslit, currpix)
          slitfit = omi_solslitfit(:, :, :, currpix) 
        ENDIF
      ELSE 
        solwinfit = omi_solwinfit(:, :, :, currpix)
      END IF
    ENDIF

    RETURN
  END SUBROUTINE omi_adj_solar_data


  SUBROUTINE omi_adj_earthshine_data (theline, pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, &
         n_max_fitpars, solar_idx!, rsl_idx, fsl_idx, comm_idx, com1_idx
    USE OMSAO_parameters_module, ONLY: normweight, max_fit_pts, &
         maxchlen!, mrefl, deg2rad
    USE OMSAO_variables_module,  ONLY: curr_rad_spec, curr_sol_spec, &
         n_rad_wvl, use_meas_sig, numwin, nradpix, the_sza_atm, the_vza_atm, &
         the_aza_atm, the_sca_atm, the_month, the_year, the_day, the_lon, &
         the_lat, the_lats, the_lons, edgelons, edgelats, the_surfalt, nview, &
         nloc, the_utc, n_radwvl_sav, radwvl_sav,  nradpix_sav, saa_minlat, &
         saa_maxlat, saa_minlon, saa_maxlon, saa_minlat1, saa_maxlat1, &
         saa_minlon1, saa_maxlon1, do_bandavg, refidx, fitvar_rad_saved, &
         n_fitvar_rad, radwavcal_freq, currpix, currloop, currline, &
         n_irrad_wvl, nsolpix, actspec_rad, database, band_selectors, &
         mask_fitvar_rad, radnhtrunc, refnhextra, &
         tempo_syn, GranuleYear, GranuleMonth, GranuleDay, &
         tabdir!, refwvl, refidx_sav, lo_radbnd, up_radbnd, i0sav
    USE OMSAO_omidata_module,    ONLY: omi_nwav_rad, omi_nradpix, &
         omi_radiance_wavl, omi_radiance_spec, omi_radiance_prec, &
         omi_specr, &!omi_szenith, omi_vzenith, omi_longitude, omi_latitude, &
         !omi_esca, omi_eaza, &
         omi_initval, omi_exitval, nlines_max, omi_fitvar, omi_radnorm, &
         nfxtrack, radwind, omi_irradiance_wavl, &
         omi_irradiance_spec, omi_irradiance_prec, orbnumsol, mswath, &
         nxtrack_max, glint_flg, land_water_flg, snow_ice_flg!, omi_saa_flag, &
         !omi_irrad_stray, omi_rad_stray, nxtrack, omi_height, omi_time, &
    USE OMSAO_pixelcorner_module, ONLY: omi_allclon, omi_allclat, &
         omi_allelon, omi_allelat, omi_allHeight, omi_alltime, &
         omi_alllat, omi_alllon, omi_allsza, omi_allvza, omi_allaza, &
         omi_allsca
    USE OMSAO_omicloud_module, ONLY: OMIL2_clouds 
    USE ozprof_data_module, ONLY: div_rad, div_sun, rad_posr, rad_specr, &
         nsaa_spike, saa_flag, the_cfrac, the_ctp, the_cld_flg, which_cld, &
         the_cod, the_orig_cfr, the_orig_ctp, scacld_initcod, the_orig_cod, &
         the_ai, radcalwrt, biasfname, biascorr, l1l2inp_unit, &
         which_biascorr, nrefl, aerosol, which_aerosol, scale_aod, &
         scaled_aod, do_simu, the_fixalb, do_lambcld, lambcld_refl, &
         has_glint, glintprob, sun_posr, sun_specr, pos_alb, &
         the_snowice!, ozprof_start_index, ozprof_end_index, nir
    USE OMSAO_errstat_module 
    use m_angle_sat2toa, only: sunglint_probability !, adjust_angle
    use avg_band, only: avg_band_radspec
    use m_ezspline_interpolation, only: bspline, interpol
    use get_cloud, only: get_tomsv8_ctp
    use m_prepare_databases
    use adj_measurement_data, only: rough_spike_detect

    IMPLICIT NONE

    ! =================
    ! In/Out variables
    ! =================
    INTEGER, INTENT (IN)        :: theline
    INTEGER, INTENT (OUT)       :: pge_error_status

    ! =================
    ! Local variables
    ! ================= 
    INTEGER :: hour, minute, fidx, lidx, i, j, west_idx, south_idx, idxoff, &
         nhtrunc, ntrunc, ntrunc1, errstat, ntempx, nch, ix, nord, ch, nw, &
         is, nsub, idum, iw
    INTEGER (KIND=i4)           :: estat
    REAL (KIND=r8)              :: second, finit
    REAL (KIND=dp), DIMENSION (n_max_fitpars) :: fitvar
    LOGICAL                     :: redo_database

    ! xliu (02/03/2007): variables for correcting across-track dependent biases
    INTEGER, PARAMETER  :: maxord = 12
    INTEGER, DIMENSION (mswath), SAVE :: corr_npars, nxcorr, nxwav
    LOGICAL, SAVE   :: first = .TRUE.
    REAL (KIND=dp), DIMENSION(mswath, nxtrack_max, 0:maxord), SAVE :: &
         corrpars, offset_pars, slope_pars
    REAL (KIND=dp), DIMENSION(mswath), SAVE    :: corr_woffset
    REAL (KIND=dp), DIMENSION(max_fit_pts)     :: corr, offset, slope
    REAL (KIND=dp), DIMENSION(nxtrack_max, max_fit_pts), SAVE :: allcorr, &
         alloffset
    REAL (KIND=dp), DIMENSION(mswath, nxtrack_max), SAVE  :: xcorr
    REAL (KIND=dp), DIMENSION(mswath,nxtrack_max,max_fit_pts),SAVE :: xwcorr, &
         xwslp, xwoff!, cldclrdf_xwcorr
    REAL (KIND=dp), DIMENSION(mswath, max_fit_pts), SAVE    :: xwavs
    INTEGER, SAVE           :: nxgascorr, nxw2corr
    REAL (KIND=dp), DIMENSION(nxtrack_max, max_fit_pts, 2), SAVE :: gascorr, &
         xw2corr
    INTEGER, DIMENSION(nxtrack_max, 3), SAVE  :: gascorr_npts, xw2corr_npts
    REAL (KIND=dp), DIMENSION(1:maxord, max_fit_pts)  :: del
    !REAL (KIND=dp), DIMENSION (2, max_fit_pts)        :: strayspec
    REAL (KIND=dp)                    :: woffset, rad347, irad347
    real (kind=dp), dimension(1) :: temp_pos_alb, temp_rad
    CHARACTER (LEN=maxchlen) :: gascorr_fname, xw2corr_fname!, &
         !cldclrdf_biasfname
    CHARACTER (LEN=255)      :: msg !! Kai

    ! ================================
    !   External functions
    ! ================================
    INTEGER (KIND=i4), EXTERNAL :: PGS_TD_TAItoUTC  
    INTEGER :: OMI_SMF_setmsg

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=23), PARAMETER :: modulename = 'omi_adj_earthshine_data'

    pge_error_status = pge_errstat_ok

    IF (first .AND. biascorr) THEN
      WRITE(msg, *) TRIM(ADJUSTL(biasfname))//',which_biascorr=',which_biascorr
      errstat = OMI_SMF_setmsg (OMI_W_GENERAL, TRIM(msg), modulename, 0)

      ! It is much better to use direct correction instead of parameterized correction
      IF ( which_biascorr == 1 ) THEN  ! Direct correction ( Y' = Y * c )
        ! Note that the nw has to be consistent with # of wavelength used in the retrieval
        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='unknown')
        READ(l1l2inp_unit, *) ntempx, nw
        DO i = 1, nw
          READ(l1l2inp_unit, *) allcorr(1:ntempx, i)
        ENDDO
        CLOSE(l1l2inp_unit)

      ELSE IF ( which_biascorr == 2 ) THEN  ! Same as 1 but parameterized as a function of wavelength

        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='OLD', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
               TRIM(ADJUSTL(biasfname)), modulename, 0)
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF

        READ(l1l2inp_unit, *) ntempx, nch
        READ(l1l2inp_unit, *) corr_npars(1:nch), corr_woffset(1:nch)

        DO ix = 1, ntempx
          READ(l1l2inp_unit, *)
          DO i = 1, nch
            READ(l1l2inp_unit, *) corrpars(i, ix, 0:corr_npars(i))
          ENDDO
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 

      ELSE IF (which_biascorr == 3 .OR. which_biascorr == 5 ) THEN  ! Direct correction ( Y' = Y * C + O)

        ! Note that the nw has to be consistent with # of wavelength used in the retrieval
        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='unknown')
        READ(l1l2inp_unit, *) ntempx, nw
        READ(l1l2inp_unit, *)
        DO i = 1, nw
          READ(l1l2inp_unit, *) allcorr(1:ntempx, i)
        ENDDO
        READ(l1l2inp_unit, *)
        DO i = 1, nw
          READ(l1l2inp_unit, *) alloffset(1:ntempx, i)
        ENDDO
        CLOSE(l1l2inp_unit)       
      ELSE IF ( which_biascorr == 4 ) THEN ! Same as 3 but parameterized as a function of wavelength

        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='OLD', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
               TRIM(ADJUSTL(biasfname)), modulename, 0)
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF

        READ(l1l2inp_unit, *) ntempx, nch
        READ(l1l2inp_unit, *) corr_npars(1:nch), corr_woffset(1:nch)

        DO ix = 1, ntempx
          READ(l1l2inp_unit, *)
          DO i = 1, nch
            READ(l1l2inp_unit, *) offset_pars(i, ix, 0:corr_npars(i)), slope_pars(i, ix, 0:corr_npars(i))
          ENDDO
        ENDDO

        CLOSE(UNIT=l1l2inp_unit) 
      ELSE IF ( which_biascorr == 6) THEN
        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='OLD', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
               TRIM(ADJUSTL(biasfname)), modulename, 0)
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF

        xcorr = 1.0  ! Initialize to 1.0
        DO is = 1, mswath
          READ (l1l2inp_unit, *) nxcorr(is)
          READ (l1l2inp_unit, *) xcorr(is, 1:nxcorr(is))

          IF (nxcorr(is) > nfxtrack) THEN
            nsub = nxcorr(is) / nfxtrack
            DO ix = 1, nfxtrack
              fidx = (ix - 1) * nsub + 1
              lidx = fidx + nsub - 1
              xcorr(is, ix) = SUM(xcorr(is, fidx:lidx)) / nsub
            ENDDO
          ENDIF
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 

      ELSE IF ( which_biascorr == 7) THEN
        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='OLD', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
               TRIM(ADJUSTL(biasfname)), modulename, 0)
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
        xwcorr = 1.0 ! Initialize to one
        DO is = 1, mswath
          READ (l1l2inp_unit, *) nxcorr(is), nxwav(is)
          DO iw = 1, nxwav(is)
            READ (l1l2inp_unit, *) xwavs(is, iw), xwcorr(is, 1:nxcorr(is), iw)
          ENDDO

          IF (nxcorr(is) > nfxtrack) THEN
            nsub = nxcorr(is) / nfxtrack
            DO ix = 1, nfxtrack
              fidx = (ix - 1) * nsub + 1
              lidx = fidx + nsub - 1
              DO iw = 1, nxwav(is)
                xwcorr(is, ix, iw) = SUM(xwcorr(is, fidx:lidx, iw)) / nsub
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 

        ! Correction for trace gases
        gascorr_fname = ADJUSTL(TRIM(tabdir)) // 'SoftCal/OMIO3PROF_corrgas_hres_o10582.dat'

        OPEN (UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(gascorr_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          WRITE(*, '(2A)') modulename, ': Cannot open trace gas correction file!!!'
          errstat = pge_errstat_error
          RETURN
        END IF

        READ (l1l2inp_unit, *)
        READ (l1l2inp_unit, *)
        READ(l1l2inp_unit, *) 
        READ (l1l2inp_unit, *) nxgascorr
        gascorr(:, :, 2) = 1.d0
        DO ix = 1, nxgascorr
          READ (l1l2inp_unit, *) idum, gascorr_npts(ix, 1:3)
          DO i = 1, gascorr_npts(ix, 1)
            READ (l1l2inp_unit, *) gascorr(ix, i, 1:2)
            gascorr(ix, i, 2) = EXP(gascorr(ix, i, 2))
          ENDDO
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 

        ! Additional correction for x-track dependent biases
        xw2corr_fname = ADJUSTL(TRIM(tabdir)) // 'SoftCal/OMIO3PROF_hres_xwcorr-2006m071116.dat'

        OPEN (UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(xw2corr_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          WRITE(*, '(2A)') modulename, ': Cannot open additional x-track dependent correction file!!!'
          errstat = pge_errstat_error
          RETURN
        END IF

        READ (l1l2inp_unit, *)
        READ (l1l2inp_unit, *) nxw2corr
        xw2corr(:, :, 2) = 0.0d0
        DO ix = 1, nxw2corr
          READ (l1l2inp_unit, *) idum, xw2corr_npts(ix, 1:3)
          DO i = 1, xw2corr_npts(ix, 1)
            READ (l1l2inp_unit, *) xw2corr(ix, i, 1:2)
          ENDDO
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 

      ELSE IF ( which_biascorr == 8 .OR. which_biascorr == 9) THEN
        OPEN(UNIT=l1l2inp_unit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='OLD', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
               TRIM(ADJUSTL(biasfname)), modulename, 0)
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
        xwslp = 0.0
        xwoff = 0.0  ! Initialize to zero
        DO is = 1, mswath
          READ (l1l2inp_unit, *) nxcorr(is), nxwav(is)
          DO iw = 1, nxwav(is)
            READ (l1l2inp_unit, *) xwavs(is, iw), xwoff(is, 1:nxcorr(is), iw), xwslp(is, 1:nxcorr(is), iw)
          ENDDO

          IF (nxcorr(is) > nfxtrack) THEN
            nsub = nxcorr(is) / nfxtrack
            DO ix = 1, nfxtrack
              fidx = (ix - 1) * nsub + 1
              lidx = fidx + nsub - 1
              DO iw = 1, nxwav(is)
                xwoff(is, ix, iw) = SUM(xwoff(is, fidx:lidx, iw)) / nsub
                xwslp(is, ix, iw) = SUM(xwslp(is, fidx:lidx, iw)) / nsub
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        CLOSE(UNIT=l1l2inp_unit) 
      ENDIF

      first = .FALSE.
    ENDIF

    ! Radiance Spectrum
    n_rad_wvl = omi_nwav_rad(currpix, currloop) 
    div_rad   = omi_radnorm (currpix, currloop)

    curr_rad_spec(wvl_idx, 1:n_rad_wvl) = omi_radiance_wavl(1:n_rad_wvl, currpix, currloop) 
    curr_rad_spec(spc_idx, 1:n_rad_wvl) = omi_radiance_spec(1:n_rad_wvl, currpix, currloop)   
    IF (use_meas_sig) THEN
      curr_rad_spec(sig_idx, 1:n_rad_wvl) = omi_radiance_prec(1:n_rad_wvl, currpix, currloop)
    ELSE
      curr_rad_spec(sig_idx, 1:n_rad_wvl) = normweight
    ENDIF
    nradpix(1:numwin) = omi_nradpix(1:numwin, currpix, currloop)     ! Solar Spectrum
    n_irrad_wvl = n_rad_wvl
    !div_sun     = omi_solnorm(currpix)

    IF ( biascorr .AND. which_biascorr == 7 ) THEN
      fidx = 1
      DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        ch = band_selectors(i) 
        !print *, ch, currpix, fidx, lidx, nradpix(i), nxwav(ch)
        CALL INTERPOL(xwavs(ch, 1:nxwav(ch)), xwcorr(ch, currpix, 1:nxwav(ch)), nxwav(ch), &
             curr_rad_spec(wvl_idx, fidx:lidx),  corr(1:nradpix(i)), nradpix(i), errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
          errstat = pge_errstat_error
          RETURN
        ENDIF

        !DO j = 1, nradpix(i)
        !   WRITE(77, *) curr_rad_spec(wvl_idx, fidx+j-1), corr(j)
        !ENDDO

        curr_rad_spec(spc_idx, fidx:lidx) = curr_rad_spec(spc_idx, fidx:lidx) / corr(1:nradpix(i))

        !! Get cldclrdf spectra
        !IF (nir > 0) THEN
        !   CALL INTERPOL(xwavs(ch, 1:nxwav(ch)), cldclrdf_xwcorr(ch, currpix, 1:nxwav(ch)), nxwav(ch), &
        !        curr_rad_spec(wvl_idx, fidx:lidx),  strayspec(2, fidx:lidx), nradpix(i), errstat)
        !   IF (errstat < 0) THEN
        !      WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
        !      errstat = pge_errstat_error; RETURN
        !   ENDIF
        !   
        !ENDIF
        fidx = lidx + 1
      ENDDO
    ENDIF

    curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = omi_irradiance_wavl(radwind(1:n_rad_wvl, currpix, currloop), currpix) 
    curr_sol_spec(spc_idx, 1:n_irrad_wvl) = omi_irradiance_spec(radwind(1:n_rad_wvl, currpix, currloop), currpix) 
    IF (use_meas_sig .AND. orbnumsol /= 99999) THEN
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = omi_irradiance_prec(radwind(1:n_rad_wvl, currpix, currloop), currpix)
    ELSE IF (orbnumsol == 99999) THEN
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = 0.0  ! Ignore error in solar irradiance
    ELSE
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
    ENDIF
    nsolpix(1:numwin) = nradpix(1:numwin)  
    !strayspec(1, 1:n_irrad_wvl) = omi_irrad_stray(radwind(1:n_rad_wvl, currpix, currloop), currpix) 
    !strayspec(2, 1:n_irrad_wvl) = omi_rad_stray(radwind(1:n_rad_wvl, currpix, currloop), currpix) !* div_sun / div_rad

    !DO i = 1, n_rad_wvl
    !   WRITE(www_lun, *) i, curr_sol_spec(1, i), curr_rad_spec(1, i)
    !ENDDO

    ! Reflectance spectrum at ~370 nm
    rad_posr (1:nrefl) = omi_specr(wvl_idx, 1:nrefl, currpix, currloop)
    rad_specr(1:nrefl) = omi_specr(spc_idx, 1:nrefl, currpix, currloop)

    IF (biascorr .AND. (which_biascorr == 8 .OR. which_biascorr == 9)) THEN
      IF (which_biascorr == 8) THEN   
        temp_pos_alb=pos_alb
        CALL BSPLINE(sun_posr(1:nrefl), sun_specr(1:nrefl), nrefl, &
             temp_pos_alb, temp_rad, 1, errstat)
        irad347=temp_rad(1)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
          STOP 1
        ENDIF
        CALL BSPLINE(rad_posr(1:nrefl), rad_specr(1:nrefl), nrefl, &
             temp_pos_alb, temp_rad, 1, errstat)
        rad347=temp_rad(1)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
          STOP 1
        ENDIF

        rad347 = rad347 / irad347
      ENDIF

      fidx = 1
      DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        ch = band_selectors(i) 
        CALL INTERPOL(xwavs(ch, 1:nxwav(ch)), xwoff(ch, currpix, 1:nxwav(ch)), nxwav(ch), &
             curr_rad_spec(wvl_idx, fidx:lidx),  offset(1:nradpix(i)), nradpix(i), errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
          errstat = pge_errstat_error
          RETURN
        ENDIF
        CALL INTERPOL(xwavs(ch, 1:nxwav(ch)), xwslp(ch, currpix, 1:nxwav(ch)), nxwav(ch), &
             curr_rad_spec(wvl_idx, fidx:lidx), slope(1:nradpix(i)), nradpix(i), errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
          errstat = pge_errstat_error
          RETURN
        ENDIF
        IF (which_biascorr == 8) THEN
          corr(1:nradpix(i)) = offset(1:nradpix(i)) + slope(1:nradpix(i)) * rad347
          corr(1:nradpix(i)) = corr(1:nradpix(i)) * curr_sol_spec(spc_idx, fidx:lidx) * div_sun / div_rad

          !DO j = 1, nradpix(i)
          !   WRITE(78, *) curr_rad_spec(wvl_idx, fidx+j-1), curr_rad_spec(spc_idx, fidx+j-1) &
          !        / (curr_rad_spec(spc_idx, fidx+j-1) - corr(j))
          !ENDDO          

          curr_rad_spec(spc_idx, fidx:lidx) = curr_rad_spec(spc_idx, fidx:lidx) - corr(1:nradpix(i))
        ELSE IF (which_biascorr == 9) THEN
          corr(1:nradpix(i)) = offset(1:nradpix(i)) + slope(1:nradpix(i)) * &
               curr_rad_spec(spc_idx, fidx:lidx) / curr_sol_spec(spc_idx, fidx:lidx) * div_rad / div_sun
          corr(1:nradpix(i)) = corr(1:nradpix(i)) / 100. + 1.0

          !DO j = 1, nradpix(i)
          !   WRITE(79, *) curr_rad_spec(wvl_idx, fidx+j-1), corr(j)
          !ENDDO
          curr_rad_spec(spc_idx, fidx:lidx) = curr_rad_spec(spc_idx, fidx:lidx) / corr(1:nradpix(i)) 
        ENDIF
        fidx = lidx + 1
      ENDDO
    ENDIF

    IF ( biascorr ) THEN
      IF ( which_biascorr == 6) THEN
        rad_specr(1:nrefl) = rad_specr(1:nrefl) / xcorr(mswath, currpix)
      ELSE IF ( which_biascorr == 7 ) THEN
        !print *, mswath, currpix, xwavs(mswath, nxwav(mswath)), rad_posr(1)
        IF ( xwavs(mswath, nxwav(mswath)) < rad_posr(1) ) THEN
          rad_specr(1:nrefl) = rad_specr(1:nrefl) / xwcorr(mswath, currpix, nxwav(mswath))
          !print *, xwcorr(mswath, currpix, nxwav(mswath))
        ELSE
          fidx = MINVAL( MINLOC( xwavs(mswath, 1:nxwav(mswath)), MASK = &
               (xwavs(mswath, 1:nxwav(mswath)) > rad_posr(1) )))
          lidx = MINVAL(MAXLOC( xwavs(mswath, 1:nxwav(mswath)), MASK = &
               (xwavs(mswath, 1:nxwav(mswath)) < rad_posr(nrefl) )))
          IF (fidx > lidx) THEN
            idum = fidx
            fidx = lidx
            lidx = idum
          ENDIF
          rad_specr(1:nrefl) = rad_specr(1:nrefl) / &
               (SUM(xwcorr(mswath, currpix, fidx:lidx)) / (lidx - fidx + 1))
        ENDIF
      ENDIF
    ENDIF

    nview       = 1
!    the_sza_atm = omi_szenith  (currpix, currloop)
    the_sza_atm = omi_allsza    (currpix, currloop)
!    the_vza_atm = omi_vzenith  (currpix, currloop)
    the_vza_atm = omi_allvza    (currpix, currloop)
!    the_aza_atm = omi_eaza     (currpix, currloop)
    the_aza_atm = omi_allaza    (currpix, currloop)
!    the_sca_atm = omi_esca     (currpix, currloop)
    the_sca_atm = omi_allsca    (currpix, currloop)
!    the_lon     = omi_longitude(currpix, currloop)
    the_lon     = omi_alllon    (currpix, currloop)
!    the_lat     = omi_latitude (currpix, currloop)
    the_lat     = omi_alllat    (currpix, currloop)
!    the_surfalt = omi_height   (currpix, currloop) / 1000.
    the_surfalt = omi_allHeight (currpix, currloop) / 1000.

    !print *, the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm
    !print *, the_surfalt
    !IF ( the_surfalt /= 0.0 ) CALL ADJUST_ANGLE(the_surfalt, the_sza_atm, &
    !     the_vza_atm, the_aza_atm, the_sca_atm)
    !print *, the_sza_atm, the_vza_atm, the_aza_atm, the_sca_atm

    nloc          = 5
    the_lons(1) = omi_allclon  (currpix-1, currline)
    the_lons(2) = omi_allclon  (currpix-1, currline + 1)
    the_lons(3) = omi_allclon  (currpix,   currline + 1)
    the_lons(4) = omi_allclon  (currpix,   currline)
    the_lons(5) = the_lon
    the_lats(1) = omi_allclat  (currpix-1, currline)
    the_lats(2) = omi_allclat  (currpix-1, currline + 1)
    the_lats(3) = omi_allclat  (currpix,   currline + 1)
    the_lats(4) = omi_allclat  (currpix,   currline)
    the_lats(5) = the_lat
    edgelons(1) = omi_allelon  (currpix-1, currline)
    edgelons(2) = omi_allelon  (currpix,   currline)
    edgelats(1) = omi_allelat  (currpix-1, currline)
    edgelats(2) = omi_allelat  (currpix,   currline)
    if (tempo_syn) then
      the_year  = GranuleYear
      the_month = GranuleMonth
      the_day   = GranuleDay
    else
      estat = PGS_TD_TAItoUTC(omi_alltime(currloop), the_utc)
      READ (the_utc, '(I4, 1x, I2, 1x, I2, 1x, I2, 1x, I2, 1x, F9.6)') &
         the_year, the_month, the_day, hour, minute, second
    endif

    ! Check cloud fractions
    IF (which_cld /= 2) THEN
      the_cld_flg = OMIL2_clouds%qflags(currpix, currline)
      the_ai      = OMIL2_clouds%ai    (currpix, currline)
      IF (the_cld_flg /= 10) THEN  ! Bad clouds for 10
        the_cfrac = OMIL2_clouds%cfr(currpix, currline)
        the_ctp   = OMIL2_clouds%ctp(currpix, currline)
      ELSE
        the_cfrac = 0.0
        the_ctp = 0.0
      ENDIF
    ENDIF

    IF (which_cld == 2 .OR. (the_cld_flg == 10 .AND. which_cld >= 3 ))  THEN
      the_cld_flg = 2  ! ISCCP
      CALL GET_TOMSV8_CTP(the_month, the_day, the_lon, the_lat, the_ctp, pge_error_status)
      the_cfrac = 0.5  ! will be updated anyway at longer wavelength
      the_ai = -999.0
    ENDIF
    print * , 'PLEASE REMOVE AFTER TESTING:ctp and cfrac is fixed to'
    the_cfrac =0.086421
    the_ctp   = 732.18
    ! Special treatments for sea glint
    has_glint = .FALSE.
    glintprob = 0.0
    ! Land-water flag=1: >=8 not used, else contain water
    IF  (land_water_flg(currpix, currloop) /= 1 .AND. land_water_flg(currpix, currloop) < 8) THEN 
      IF (glint_flg(currpix, currloop) == 1) THEN
        has_glint = .TRUE.
        CALL SUNGLINT_PROBABILITY (the_sza_atm, the_vza_atm, the_aza_atm, glintprob)
        !PRINT *, 'Glint Probability: ', glintprob, the_cfrac
        !IF (the_cfrac < 0.30 * glintprob) the_cfrac = 0.0
      ENDIF
    ENDIF

    ! Snow/ice flag
    ! 00: Snow-free land
    ! 1-100: sea ice concentration
    ! 101: permanent ice
    ! 102: not used
    ! 103: snow
    ! 104: ocean
    !105-123: Reserved
    !124: mixed pixels
    !125: suspect ice value
    !126: corners
    !17: Error
    the_snowice = snow_ice_flg(currpix, currloop)

    IF ( do_lambcld ) THEN
      the_cod = 0.0
    ELSE
      ! Pixel-independent approximation: cloudy scence with an effective COD 20.0 
      ! (cloud thickness 100 mb) and clear-sky scene. If cloud fraction is 20, 
      ! then rederive the effective COD.  Since CTP from OMI products are based on
      ! Lambertian cloud model, it is better to assume thin cloud layer (e.g., 100 mb)
      ! even for thick clouds
      the_cod = scacld_initcod
    ENDIF

    IF (do_simu .AND. .NOT. radcalwrt) THEN
      OPEN(UNIT=l1l2inp_unit, FILE='INP/sim.inp', STATUS='unknown')
      READ(l1l2inp_unit, *) the_sza_atm, the_vza_atm, the_aza_atm, the_fixalb, the_surfalt, &
           the_cfrac, the_ctp, the_cod, the_lon, the_lat, the_month, the_day, which_aerosol, &
           scaled_aod, do_lambcld, lambcld_refl
      IF (which_aerosol < 0 ) THEN
        aerosol = .FALSE.
        scale_aod = .FALSE.
        scaled_aod = 0.0
      ELSE
        aerosol = .TRUE.
        scale_aod = .TRUE.
      ENDIF

      IF (the_cfrac == 0.0 .OR. the_cod == 0) THEN
        the_ctp = 0.0
        the_cod = 0.0
        the_cfrac = 0
      ENDIF

      IF (do_lambcld ) THEN
        the_cod = 0.0
      ENDIF

      CLOSE (l1l2inp_unit)
    ENDIF

    ! These properties may be slightly modified later (save them)
    the_orig_cfr = the_cfrac
    the_orig_ctp = the_ctp
    the_orig_cod = the_cod

    ! Detect Spikes over the South Atlantic Anomaly region 
    IF (.NOT. do_simu .AND. .NOT. radcalwrt) THEN
      IF ( ((the_lat > saa_minlat  .AND. the_lat < saa_maxlat    .AND. &
           the_lon  > saa_minlon  .AND. the_lon < saa_maxlon)   .OR.  & 
           (the_lat > saa_minlat1 .AND. the_lat < saa_maxlat1   .AND. &
           the_lon  > saa_minlon1 .AND. the_lon < saa_maxlon1)) .AND. .NOT. do_simu )  THEN   
        CALL ROUGH_SPIKE_DETECT(n_rad_wvl, curr_rad_spec(wvl_idx, 1:n_rad_wvl), &
             curr_rad_spec(spc_idx, 1:n_rad_wvl), curr_sol_spec(spc_idx, 1:n_rad_wvl), nsaa_spike)
        saa_flag = .FALSE. !; nsaa_spike = 0
      ELSE
        saa_flag = .FALSE.
        nsaa_spike = 0
      ENDIF
    ELSE
      saa_flag = .FALSE.
      nsaa_spike = 0
    ENDIF

    ! Obtain measurement error in term sun-normalized radiance
    CALL omi_adj_rad_sig (curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl),  &
         curr_sol_spec(wvl_idx:sig_idx, 1:n_rad_wvl), n_rad_wvl)

    ! Make sure that reference spectra has  more wavelengths than
    ! irradiance interpolation and shifting      
    nhtrunc = radnhtrunc
    ntrunc = nhtrunc * 2
    ntrunc1 = ntrunc + 1

    fidx = 1
    DO i = 1, numwin
      lidx = fidx + nradpix(i) - ntrunc1
      curr_rad_spec(1:sig_idx, fidx:lidx) = curr_rad_spec(1:sig_idx, fidx + nhtrunc : lidx + nhtrunc)
      !strayspec(1:2, fidx:lidx) = strayspec(1:2, fidx + nhtrunc : lidx + nhtrunc)  
      !strayspec(2, fidx:lidx) = strayspec(2, fidx + nhtrunc : lidx + nhtrunc)    
      IF (lidx  < n_rad_wvl - ntrunc1 ) THEN
        curr_rad_spec(1:sig_idx, lidx+1:n_rad_wvl - ntrunc) =  &
             curr_rad_spec(1:sig_idx, lidx+ntrunc1:n_rad_wvl)
        !strayspec(1:2, lidx+1:n_rad_wvl - ntrunc) =  strayspec(1:2, lidx+ntrunc1:n_rad_wvl)
        !strayspec(2, lidx+1:n_rad_wvl - ntrunc) =  strayspec(2, lidx+ntrunc1:n_rad_wvl)
      ENDIF

      nradpix(i) = nradpix(i) - ntrunc
      fidx = lidx + 1
      n_rad_wvl = n_rad_wvl - ntrunc
    ENDDO

    ! save the original grids for later obtain ozone cross section
    n_radwvl_sav = n_rad_wvl
    nradpix_sav = nradpix
    radwvl_sav(1:n_rad_wvl) = curr_rad_spec(wvl_idx, 1:n_rad_wvl)  

    redo_database = .FALSE.
    IF (theline >= 1) THEN
      IF (omi_nwav_rad(currpix, currloop) /= omi_nwav_rad(currpix, currloop-1)) THEN
        redo_database = .TRUE.
      ELSE 
        IF (ANY(radwind(1:n_rad_wvl, currpix, currloop) - &
             radwind(1:n_rad_wvl, currpix, currloop-1) /= 0)) redo_database = .TRUE.
      ENDIF
    ENDIF

    IF ( MOD (theline, radwavcal_freq) == 0 .OR. redo_database) THEN 
      ! --------------------------------------------------------------
      ! Spline data bases, compute undersampling spectrum, and prepare
      ! reference spectra for fitting.
      ! --------------------------------------------------------------
      CALL prepare_databases ( n_rad_wvl, curr_rad_spec(wvl_idx,1:n_rad_wvl), pge_error_status )
      IF ( pge_error_status >= pge_errstat_error ) RETURN
    ENDIF

    ! average and subsampling for selected bands and update nradpix
    IF (do_bandavg) THEN 
      CALL avg_band_radspec (curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl), &
           n_rad_wvl, pge_error_status)
      IF ( pge_error_status >= pge_errstat_error ) RETURN
    ENDIF

    fidx = 1
    DO i = 1, numwin
      lidx = fidx + nradpix(i) - 1
      idxoff = refnhextra + (i - 1) * 2 * refnhextra
      refidx(fidx:lidx) = (/(j, j = fidx + idxoff, lidx + idxoff)/)
      fidx = lidx  + 1
    ENDDO
    IF (radcalwrt) THEN
      actspec_rad(1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) / &
           database(solar_idx, refidx(1:n_rad_wvl)) * div_rad / div_sun
    ENDIF

    ! load databases for common modes
    IF ( MOD (theline, radwavcal_freq) == 0 .OR. redo_database) THEN
      CALL load_omi_comres(pge_error_status)
      IF ( pge_error_status >= pge_errstat_error ) RETURN
    ENDIF

    !! restore straylight spectra to database
    !IF ( MOD (theline, radwavcal_freq) == 0 .OR. redo_database) THEN 
    !   !database(fsl_idx, refidx(1:n_rad_wvl)) = strayspec(1, 1:n_rad_wvl)
    !   database(rsl_idx, refidx(1:n_rad_wvl)) = strayspec(2, 1:n_rad_wvl)
    !
    !   !DO i = 1, n_rad_wvl
    !   !   WRITE(90, *) curr_rad_spec(_idx, i), strayspec(2, i)
    !   !ENDDO        
    !ENDIF

    IF (biascorr) THEN

      IF ( which_biascorr == 1 ) THEN

        curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) / allcorr(currpix, 1:n_rad_wvl)

      ELSE IF (which_biascorr == 2) THEN
        fidx = 1
        DO i = 1, numwin
          lidx = fidx + nradpix(i) - 1
          ch = band_selectors(i)
          nord = corr_npars(ch)
          woffset = corr_woffset(ch)

          del(1, fidx:lidx) = curr_rad_spec(wvl_idx, fidx:lidx) - woffset
          DO j = 2, nord
            del(j, fidx:lidx) = del(j-1, fidx:lidx) * del(1, fidx:lidx)
          ENDDO

          corr(fidx:lidx) = corrpars(ch, currpix, 0)
          DO j = 1, nord
            corr(fidx:lidx) =  corr(fidx:lidx) + corrpars(ch, currpix, j) * del(j, fidx:lidx)
          ENDDO

          fidx = lidx + 1
        ENDDO

        corr(1:n_rad_wvl) = 1.0 + corr(1:n_rad_wvl) / 100.0
        curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) / corr(1:n_rad_wvl)

      ELSE IF ( which_biascorr == 3 ) THEN

        curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) * allcorr(currpix, 1:n_rad_wvl)  &
             + alloffset(currpix, 1:n_rad_wvl) * database(solar_idx, refidx(1:n_rad_wvl) ) * div_sun / div_rad

      ELSE IF ( which_biascorr == 4 ) THEN

        fidx = 1
        DO i = 1, numwin
          lidx = fidx + nradpix(i) - 1
          ch = band_selectors(i)
          nord = corr_npars(ch)
          woffset = corr_woffset(ch)

          del(1, fidx:lidx) = curr_rad_spec(wvl_idx, fidx:lidx) - woffset
          DO j = 2, nord
            del(j, fidx:lidx) = del(j-1, fidx:lidx) * del(1, fidx:lidx)
          ENDDO

          offset(fidx:lidx) = offset_pars(ch, currpix, 0)
          slope(fidx:lidx)  = slope_pars (ch, currpix, 0)
          DO j = 1, nord
            offset(fidx:lidx) =  offset(fidx:lidx) + offset_pars(ch, currpix, j) * del(j, fidx:lidx)
            slope (fidx:lidx) =  slope (fidx:lidx) + slope_pars (ch, currpix, j) * del(j, fidx:lidx)
          ENDDO
          fidx = lidx + 1
        ENDDO

        curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) * slope(1:n_rad_wvl) &
             + offset(1:n_rad_wvl) * database(solar_idx, refidx(1:n_rad_wvl)) * div_sun / div_rad

      ELSE IF ( which_biascorr == 5 ) THEN

        curr_rad_spec(spc_idx, 1:n_rad_wvl) = EXP( allcorr(currpix, 1:n_rad_wvl) &
             * LOG( curr_rad_spec(spc_idx, 1:n_rad_wvl) ) + ( allcorr(currpix, 1:n_rad_wvl) - 1.0 )  &
             * LOG ( div_rad / ( database(solar_idx, refidx(1:n_rad_wvl)) * div_sun ) ) &
             + alloffset(currpix, 1:n_rad_wvl) )

      ELSE IF ( which_biascorr == 6 ) THEN
        fidx = 1
        DO i = 1, numwin
          lidx = fidx + nradpix(i) - 1
          ch = band_selectors(i)
          curr_rad_spec(spc_idx, fidx:lidx) = curr_rad_spec(spc_idx, fidx:lidx) / xcorr(ch, currpix)
          fidx = lidx + 1
        ENDDO
      ELSE IF ( which_biascorr == 7) THEN
        !! Kai need to comment out gascorr, i.e., the following two lines to make the code work for arbitrary window. 
        !curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) * gascorr(currpix, 1:n_rad_wvl, 2)
        !curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) * &
        !     (1.0d0 + xw2corr(currpix, 1:n_rad_wvl, 2) / 100.)
      ENDIF
    ENDIF

    !WRITE(95, *) currpix, n_rad_wvl
    !DO i = 1, n_rad_wvl
    !   WRITE(95, '(F10.4, 2D16.7)')  curr_rad_spec(wvl_idx, i), curr_rad_spec(spc_idx, i) * div_rad, &
    !        database(solar_idx, refidx(i)) * div_sun
    !ENDDO
    !OPEN(UNIT=95, FILE='corr.dat', STATUS='unknown')
    !READ(95, *) corr(1:n_rad_wvl)
    !curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) / corr(1:n_rad_wvl)
    !CLOSE(95)

    ! Initialized fitted variables from valid western and southern neighbors 
    ! since the retrievals are performed from west to east and from south to north
    west_idx  = currpix  - 1
    south_idx = currloop - 1
    IF (south_idx < 0) south_idx = nlines_max - 1

    omi_initval(currpix, currloop) = 0
    fitvar = 0.0
    finit = 0.0
    IF (west_idx > 0) THEN
      IF (omi_exitval(west_idx, currloop) > 0) THEN  ! Western pixel (success retrieval)
        fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) + &
             omi_fitvar(west_idx, currloop, 1:n_fitvar_rad)
        finit = finit + 1.0
      ENDIF

      IF (south_idx >= 0 .AND. south_idx /= nlines_max - 1) THEN
        IF (omi_exitval(west_idx, south_idx) > 0) THEN ! Southwestern pixel (success retrieval)
          fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) &
               + omi_fitvar(west_idx, south_idx, 1:n_fitvar_rad) * 0.5
          finit = finit + 0.5
        ENDIF
      ENDIF
    ENDIF

    IF ( south_idx >= 0 ) THEN
      IF (omi_exitval(currpix, south_idx) > 0) THEN     ! Southern pixel (success retrieval)
        fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) + &
             omi_fitvar(currpix, south_idx, 1:n_fitvar_rad) 
        finit = finit + 1.0
      ENDIF
    ENDIF

    IF (finit > 0) THEN
      fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad)) = fitvar(1:n_fitvar_rad) / finit
      omi_initval(currpix, currloop) = 1
    ENDIF

    RETURN
  END SUBROUTINE omi_adj_earthshine_data


  SUBROUTINE omi_adj_rad_sig (radspec, solspec, np)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
    !USE OMSAO_parameters_module, ONLY: maxwin
    !USE OMSAO_variables_module,  ONLY: numwin, nradpix, band_selectors
    USE ozprof_data_module,      ONLY: div_rad, div_sun, use_lograd, use_flns

    IMPLICIT NONE

    ! ====================
    ! In/Output variables
    ! ====================
    INTEGER, INTENT(IN)                             :: np
    REAL (KIND=dp), DIMENSION(3, np), INTENT(IN)    :: solspec
    REAL (KIND=dp), DIMENSION(3, np), INTENT(INOUT) :: radspec

    ! ====================
    ! Local variables
    ! ====================
    INTEGER, PARAMETER                :: nreg = 3
    INTEGER                           :: i, j, fidx, lidx!, iwin
    !REAL(KIND=dp)                     :: dw1, dw2
    REAL (KIND=dp), DIMENSION(np)     :: relsig, normrad, sig
    !REAL (KIND=dp), DIMENSION(maxwin) :: floor_noise =  &
    !     (/0.004, 0.002, 0.001, 0.001, 0.001/)
    REAL (KIND=dp), DIMENSION(nreg) :: reg_noise =  &
         (/0.004, 0.004, 0.002/)
    REAL (KIND=dp), DIMENSION(0:nreg) :: reg_waves = &
         (/260.0, 300.0, 310.0, 350./)

    ! I/F
    normrad = radspec(spc_idx,:) / solspec(spc_idx,:) * (div_rad / div_sun) 

    ! Relative photon noise (I)
    relsig = radspec(sig_idx,:) / radspec(spc_idx,:)
    relsig = SQRT( relsig ** 2.0 + (solspec(sig_idx, :) / solspec(spc_idx,:)) ** 2.0)

    !IF (use_flns) THEN
    !   fidx = 1
    !   DO iwin = 1, numwin
    !      lidx = fidx + nradpix(iwin) - 1
    !      IF (lidx - fidx > 0) THEN
    !         WHERE (relsig(fidx:lidx) < floor_noise(band_selectors(iwin)))               
    !            relsig(fidx:lidx) = floor_noise(band_selectors(iwin))
    !         END WHERE
    !      ENDIF
    !
    !      fidx = lidx + 1
    !   ENDDO
    !END IF

    IF (use_flns) THEN
      fidx = 1
      i = fidx
      DO j = 1, nreg
        DO WHILE (i <= np )
          IF (radspec(wvl_idx, i) < reg_waves(j)) THEN
            i = i + 1
          ELSE
            exit
          ENDIF
        ENDDO
        lidx = i - 1
        WHERE (relsig(fidx:lidx) < reg_noise(j))              
          relsig(fidx:lidx) = reg_noise(j)
        ENDWHERE
        fidx = lidx + 1
      ENDDO
    ENDIF

    IF (use_lograd) THEN        ! error in logarithmic radiance
      sig = LOG(1.0 + relsig)  ! ~relative error
    ELSE
      sig = relsig * normrad   ! absolute measurement error in I/F
    ENDIF

    radspec(sig_idx, :) = sig

    !DO i = 1, np - 1
    !   WRITE(92, *) radspec(1, i), sig(i)
    !ENDDO
    !STOP

    RETURN
  END SUBROUTINE omi_adj_rad_sig


  SUBROUTINE omi_set_fitting_parameters ( pge_error_status )

    USE OMSAO_precision_module
    USE OMSAO_omidata_module,   ONLY : nswath, omi_radiance_swathname, &
         omi_irradiance_swathname, omichs
    USE OMSAO_variables_module, ONLY : band_selectors
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER,           INTENT (OUT) :: pge_error_status

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=26), PARAMETER :: modulename = 'omi_set_fitting_parameters'

    ! --------------------------
    ! Initialize OUTPUT variable
    ! --------------------------
    pge_error_status = pge_errstat_ok

    ! ----------------------------------------------------------------
    ! Name of solar and earthshine swaths (normally obtained from PCF)
    ! ----------------------------------------------------------------
    IF (nswath == 2) THEN
      omi_irradiance_swathname(1) = 'Sun Volume UV-1 Swath' !'Sun Daily UV-1 Swath'
      omi_irradiance_swathname(2) = 'Sun Volume UV-2 Swath' !'Sun Daily UV-2 Swath'
      omi_radiance_swathname(1)   = 'Earth UV-1 Swath'      !'Earthshine Daily UV-1 Swath'
      omi_radiance_swathname(2)   = 'Earth UV-2 Swath'      !'Earthshine Daily UV-1 Swath' 
      omichs(1) = 1
      omichs(2) = 2
    ELSE IF (nswath == 1) THEN
      ! Ozone profile retrieval with channel 1 only (impossible due to always using uv2 for fc)
      IF (band_selectors(1) == 1) THEN
        omi_irradiance_swathname(1) = 'Sun Volume UV-1 Swath' !'Sun Daily UV-1 Swath'
        omi_radiance_swathname(1)   = 'Earth UV-1 Swath'      !'Earthshine Daily UV-1 Swath' 
        omichs(1) = 1
      ELSE
        ! Ozone profile retrieval with channel 2 only (total ozone retrieval)
        omi_irradiance_swathname(1) = 'Sun Volume UV-2 Swath' !'Sun Daily UV-2 Swath'
        omi_radiance_swathname(1)   = 'Earth UV-2 Swath'      !'Earthshine Daily UV-2 Swath' 
        omichs(1) = 2
      ENDIF
    ELSE
      WRITE(*, '(A)') 'Need and only need UV swathes for ozone profile retrieval!!!'
      pge_error_status = pge_errstat_error     
    ENDIF

    RETURN
  END SUBROUTINE omi_set_fitting_parameters

SUBROUTINE load_omi_comres(errstat)
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: max_fit_pts, maxchlen
    USE OMSAO_variables_module, ONLY: refdbdir, n_refwvl, database, &
         nradpix, refidx, currpix!, up_radbnd, refwvl, numwin, &
         !n_fitvar_rad, fitvar_rad_str, mask_fitvar_rad, lo_radbnd
    USE OMSAO_indices_module, ONLY: com1_idx, com_idx, comfidx, cm1fidx
    USE OMSAO_errstat_module, ONLY: pge_errstat_error, pge_errstat_ok, www_lun

    INTEGER, INTENT (OUT)      :: errstat

    CHARACTER (LEN=15), PARAMETER :: modulename = 'load_omi_comres'
    INTEGER, PARAMETER            :: nx = 30, lun = 12
    INTEGER                       :: ix, itemp, i, fidx, lidx
    CHARACTER (LEN=maxchlen)      :: comres_fname

    LOGICAL,                                        SAVE :: first = .TRUE.
    REAL (KIND=dp), DIMENSION (nx, max_fit_pts, 2), SAVE :: comres
    INTEGER, DIMENSION(nx, 3),                      SAVE :: npts

    errstat = pge_errstat_ok
    IF (first) THEN
      comres = 0.D0
      comres_fname = ADJUSTL(TRIM(refdbdir)) // 'OMI_hresjul11-30S-30N_comres.dat'

      OPEN (UNIT=lun, FILE=TRIM(ADJUSTL(comres_fname)), STATUS='UNKNOWN', IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        WRITE(www_lun, '(2A)') modulename, ': Cannot open common-mode residual file!!!'
        errstat = pge_errstat_error; RETURN
      END IF

      READ (lun, *)
      READ (lun, *)
      DO ix = 1, nx
        READ (lun, *) itemp, npts(ix, 1:3)
        DO i = 1, npts(ix, 1)
          READ (lun, *) comres(ix, i, 1:2)
        ENDDO
      ENDDO
      CLOSE (lun)

      first = .FALSE.
    ENDIF

    IF (comfidx > 0) THEN
      database(com_idx, 1:n_refwvl) = 0.D0
      fidx = 1; lidx = fidx + nradpix(1) - 1
      database(com_idx, refidx(fidx:lidx)) = comres(currpix, fidx:lidx, 2)
    ENDIF

    IF ( cm1fidx > 0 ) THEN
      database(com1_idx, 1:n_refwvl) = 0.D0
      fidx = nradpix(1) + 1; lidx = fidx + nradpix(2) - 1
      database(com1_idx, refidx(fidx:lidx)) = comres(currpix, fidx:lidx, 2)
    ENDIF

    RETURN

  END SUBROUTINE load_omi_comres

  ! Correction for wavelength registration at 1:67 and 498:557
  SUBROUTINE corruv2wav(nw, nx, waves)
    USE OMSAO_precision_module

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                  INTENT (IN)    :: nw, nx
    REAL (KIND=r4), DIMENSION (nw, nx), INTENT (INOUT) :: waves

    ! ----------------
    ! Local variables
    ! ---------------
    INTEGER        :: ix
    REAL (KIND=r4) :: del, ndel, delp1, delp2, delm1, delm2, ndelp1, ndelm1, sh1, sh2

    !RETURN

    !print *, nw, nx
    DO ix = 1, nx   
      ! At position 67
      del    = waves(68, ix) - waves(67, ix)
      delm1  = waves(67, ix) - waves(66, ix)
      delp1  = waves(69, ix) - waves(68, ix)
      delp2  = waves(70, ix) - waves(69, ix)
      ndel   = delp1 * 2 - delp2
      ndelm1 = ndel * 2  - delp1

      ! Shifts
      sh1 = (ndel - del)
      sh2 = sh1 + (ndelm1 - delm1)

      !print *, 67, sh1, sh2
      waves(67, ix)   = waves(67, ix) - sh1
      waves(1:66, ix) = waves(1:66, ix) - sh2

      ! At position 498
      delm2 = waves(496, ix) - waves(495, ix)
      delm1 = waves(497, ix) - waves(496, ix)
      del   = waves(498, ix) - waves(497, ix)
      delp1 = waves(499, ix) - waves(498, ix)
      ndel = delm1 * 2 - delm2
      ndelp1 = ndel * 2 - delm1

      sh1 = (ndel - del)
      sh2 = sh1 + (ndelp1 - delp1)
      !print *, 498, sh1, sh2

      waves(498, ix) = waves(498, ix) + sh1
      waves(499:nw, ix) = waves(499:nw, ix) + sh2
    ENDDO

    RETURN
  END SUBROUTINE corruv2wav

  subroutine timestamp (curtime )
    !
    !*******************************************************************************
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !
    !  Example:
    !
    !    May 31 2001   9:45:54.872 AM
    !
    !  Modified:
    !
    !    31 May 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none
    !
    character (len=24), INTENT(OUT) :: curtime
    integer d
    character ( len = 8 ) date
    integer h
    integer m
    integer mm
    character ( len = 3 ), parameter, dimension(12) :: month = (/ &
         'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)
    integer n
    integer s
    character ( len = 10 )  time
    integer values(8)
    integer y
    character ( len = 5 ) zone
    !
    call date_and_time ( date, time, zone, values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    !if ( h < 12 ) then
    !  ampm = 'AM'
    !else if ( h == 12 ) then
    !  if ( n == 0 .and. s == 0 ) then
    !    ampm = 'Noon'
    !  else
    !    ampm = 'PM'
    !  end if
    !else
    !  h = h - 12
    !  if ( h < 12 ) then
    !    ampm = 'PM'
    !  else if ( h == 12 ) then
    !    if ( n == 0 .and. s == 0 ) then
    !      ampm = 'Midnight'
    !    else
    !      ampm = 'AM'
    !    end if
    !  end if
    !end if

    write ( curtime, '(a3,1x,i2.2,1x,i4,1x,i2.2,a1,i2.2,a1,i2.2,a1,i3.3)' ) &
         trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm

    return
  end subroutine timestamp


  SUBROUTINE GET_DOY(theyear, themon,  theday, thedoy)
    USE OMSAO_precision_module
    IMPLICIT NONE

    ! --------------------------
    ! Input/Output variables
    ! --------------------------
    INTEGER, INTENT (IN)   :: theyear, themon, theday
    INTEGER, INTENT (OUT)  :: thedoy

    ! Local variables
    INTEGER, DIMENSION(12) :: ndays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    IF (MOD(theyear, 4) == 0) ndays(2) = 29
    IF (themon == 1) THEN
      thedoy = theday
    ELSE 
      thedoy = SUM(ndays(1:themon-1)) + theday
    ENDIF

    RETURN
  END SUBROUTINE GET_DOY
  

end module omi_fitting_aux
