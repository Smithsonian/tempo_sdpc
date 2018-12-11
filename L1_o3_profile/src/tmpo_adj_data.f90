!
MODULE tmpo_adj_data

  USE OMSAO_tmpodata_module, ONLY : nxtrack_max,nlines_max,  & 
  ring=>tmpo_ring, refl=>tmpo_refl, rad=>tmpo_rad, irrad=>tmpo_irrad,&
  cali=>tmpo_cali, geo1=>tmpo_geo1, geo2=>tmpo_geo2, o3p=>tmpo_o3p
  IMPLICIT NONE
  PUBLIC  adj_solar_data, adj_earthshine_data     
  PRIVATE adj_rad_sig, load_comres

CONTAINS

  SUBROUTINE adj_solar_data (pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, hwe_idx
    USE OMSAO_parameters_module, ONLY: normweight
    USE OMSAO_variables_module,  ONLY: curr_sol_spec, n_irrad_wvl, &
         use_meas_sig, numwin, nsol_ring, sol_spec_ring, nsolpix, &
         yn_varyslit, slit_rad, solwinfit, nslit, slitwav, slitfit, &
         sring_fidx, sring_lidx,  currpix, which_slit
    USE ozprof_data_module,      ONLY: div_sun, sun_posr, sun_specr, nrefl
    USE OMSAO_errstat_module 
    
    IMPLICIT NONE

    ! =================
    ! Output variables
    ! =================
    INTEGER, INTENT (OUT) :: pge_error_status


    pge_error_status = pge_errstat_ok

    ! Solar Spectrum
    n_irrad_wvl = irrad%nwav(currpix) 
    div_sun     = irrad%norm(currpix)

    curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = irrad%wavl(1:n_irrad_wvl, currpix) 
    curr_sol_spec(spc_idx, 1:n_irrad_wvl) = irrad%spec(1:n_irrad_wvl, currpix) 
    IF (use_meas_sig ) THEN
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = irrad%prec(1:n_irrad_wvl, currpix)
    ELSE
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
    ENDIF
    nsolpix(1:numwin) = irrad%npix(1:numwin, currpix) 

    ! Solar Spectrum for Ring Calculation
    nsol_ring = ring%nsol(currpix)
    sol_spec_ring(wvl_idx, 1:nsol_ring) = ring%wavl(1:nsol_ring, currpix)
    sol_spec_ring(spc_idx, 1:nsol_ring) = ring%spec(1:nsol_ring, currpix)
    sring_fidx = ring%winpix(currpix,1)
    sring_lidx = ring%winpix(currpix,2)

    ! Reflectance spectrum at ~370 nm
    sun_posr (1:nrefl) = refl%solwavl(1:nrefl, currpix)
    sun_specr(1:nrefl) = refl%solspec(1:nrefl, currpix)

    ! Load slit calibration parameters
    IF (which_slit < 5) THEN
      IF (yn_varyslit) THEN
        IF (slit_rad) THEN
          nslit = cali%nslit_rad(currpix)
          slitwav(1:nslit) = cali%slitwav_rad(1:nslit, currpix)
          slitfit(:,:,:) = cali%slitfit_rad(:, :, :, currpix)
        ELSE
          nslit = cali%nslit_sol(currpix)
          slitwav(1:nslit) = cali%slitwav_sol(1:nslit, currpix)
          slitfit = cali%slitfit_sol(:, :, :, currpix) 
        ENDIF
      ELSE 
        solwinfit(:,:, 1) = cali%solwinfit(:, :,  currpix)
        WRITE(www_lun,*) 'slit:', solwinfit(1:numwin, hwe_idx,1)
      END IF
    ENDIF

    RETURN
  END SUBROUTINE adj_solar_data

  SUBROUTINE adj_earthshine_data (theline, pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, &
         n_max_fitpars, solar_idx!, rsl_idx, fsl_idx, comm_idx, com1_idx
    USE OMSAO_parameters_module, ONLY: mswath, normweight, max_fit_pts, &
         maxchlen, calunit
    USE OMSAO_variables_module,  ONLY: curr_rad_spec, curr_sol_spec, &
         n_rad_wvl, use_meas_sig, numwin, nradpix, the_sza_atm, the_vza_atm, &
         the_aza_atm, the_sca_atm, the_month, the_year, the_day,the_jday, the_lon, &
         the_lat, the_lats, the_lons, edgelons, edgelats, the_surfalt, nview, &
         nloc, the_utc, n_radwvl_sav, radwvl_sav,  nradpix_sav, saa_minlat, &
         saa_maxlat, saa_minlon, saa_maxlon, saa_minlat1, saa_maxlat1, &
         saa_minlon1, saa_maxlon1, do_bandavg, refidx, fitvar_rad_saved, &
         n_fitvar_rad, radwavcal_freq, currpix, currloop, &
         n_irrad_wvl, nsolpix, actspec_rad, database, band_selectors, &
         mask_fitvar_rad, radnhtrunc, refnhextra, curr_rad_spec_save, &
         GranuleYear, GranuleMonth, GranuleDay,GranuleJDay,tabdir, currline
    USE OMSAO_omicloud_module, ONLY: OMIL2_clouds 
    USE ozprof_data_module, ONLY: div_rad, div_sun, rad_posr, rad_specr, &
         nsaa_spike, saa_flag, the_cfrac, the_ctp, the_cld_flg, which_cld, &
         the_cod, the_orig_cfr, the_orig_ctp, scacld_initcod, the_orig_cod, &
         the_ai, radcalwrt, biasfname, biascorr,  &
         which_biascorr, nrefl, aerosol, which_aerosol, scale_aod, &
         scaled_aod, do_simu, the_fixalb, do_lambcld, lambcld_refl, &
         has_glint, glintprob, sun_posr, sun_specr, pos_alb, &
         the_snowice, the_landwater_flg, the_glint_flg, the_geo1, the_geo2
    USE OMSAO_errstat_module 

    USE m_angle_sat2toa, only: sunglint_probability 
    USE m_avg_band, only: avg_band_radspec
    USE m_ezspline_interpolation, only: bspline, interpol
    USE m_get_cloud, only: get_tomsv8_ctp
    USE m_prepare_databases
    USE m_fitting_util, only: rough_spike_detect


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
    REAL (KIND=dp)              :: second, finit
    REAL (KIND=dp), DIMENSION (n_max_fitpars) :: fitvar
    LOGICAL                     :: redo_database
    REAL (KIND=dp), DIMENSION(max_fit_pts) :: corr
    ! xliu (02/03/2007): variables for correcting across-track dependent biases
    TYPE soft_group
      INTEGER, DIMENSION (mswath) ::  nxcorr, nxwav
      REAL (KIND=dp), DIMENSION(mswath, nxtrack_max) :: xcorr
      REAL (KIND=dp), DIMENSION(mswath, nxtrack_max,max_fit_pts) :: xwcorr
      REAL (KIND=dp), DIMENSION(mswath, max_fit_pts)     :: xwavs
    END TYPE soft_group
    TYPE (soft_group), SAVE :: soft
    LOGICAL, SAVE   :: first = .TRUE.
    CHARACTER (LEN=255)      :: msg !! Kai

    INTEGER, PARAMETER :: DBPRECISION = SELECTED_INT_KIND(PRECISION(1.0d0))
    INTEGER (DBPRECISION), PARAMETER :: NAN = Z"7FF8000000000000"
    INTEGER, PARAMETER :: DPSB = BIT_SIZE(NAN) - 1
    !   External functions
    ! ================================
    INTEGER (KIND=i4), EXTERNAL :: PGS_TD_TAItoUTC  
    INTEGER :: OMI_SMF_setmsg

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=23), PARAMETER :: modulename = 'tmpo_adj_earthshine_data'

    pge_error_status = pge_errstat_ok

    !--------------------------------------------------------------------------------------------
    ! geometry
    !---------------------------------------------------------------------------------------------
    nview       = 1
    the_sza_atm = geo1%sza    (currpix, currline)
    the_vza_atm = geo1%vza    (currpix, currline)
    the_aza_atm = geo1%aza    (currpix, currline)
    the_sca_atm = geo1%sca    (currpix, currline)
    the_lon     = geo1%lon    (currpix, currline)
    the_lat     = geo1%lat    (currpix, currline)
    the_surfalt = geo1%Height (currpix, currline) / 1000.
    nloc          = 5
    the_lons(1:4) = geo1%clon  (1:4,currpix, currline)
    the_lons(5) = the_lon
    the_lats(1:4) = geo1%clat  (1:4,currpix, currline)
    the_lats(5) = the_lat
    edgelons(1:2) = geo1%elon  (currpix-1:currpix, currline)
    edgelats(1:2) = geo1%elat  (currpix-1:currpix, currline)
    the_snowice = geo1%snow_ice_flg(currpix, currline)
    the_landwater_flg = geo1%land_water_flg(currpix,currline)
    the_glint_flg = geo1%glint_flg(currpix, currline)


!    estat = PGS_TD_TAItoUTC(geo1%time(currline), the_utc)
!    READ (the_utc, '(I4, 1x, I2, 1x, I2, 1x, I2, 1x, I2, 1x, F9.6)') &
     the_utc='1111-11-11-11-11-11-11.11'

    the_geo1%sza = geo1%sza (currpix, currline) 
    the_geo1%vza = geo1%vza (currpix, currline) 
    the_geo1%aza = geo1%aza (currpix, currline) 
    the_geo1%sca = geo1%sca (currpix, currline) 
    the_geo1%lon = geo1%lon (currpix, currline) 
    the_geo1%lat = geo1%lat (currpix, currline) 
    the_geo1%surfalt = geo1%height (currpix, currline) 
    the_geo1%lons(1:4) = geo1%clon(1:4,currpix, currline)
    the_geo1%lons(5) = geo1%lon(currpix, currline)
    the_geo1%lats(1:4) = geo1%clat(1:4,currpix, currline)
    the_geo1%lats(5) = geo1%lat(currpix, currline)
    the_geo1%elon(1:2) = geo1%elon(currpix-1:currpix, currline)
    the_geo1%elat(1:2) = geo1%elat(currpix-1:currpix, currline)

    the_geo2%sza = geo2%sza (currpix, currline) 
    the_geo2%vza = geo2%vza (currpix, currline) 
    the_geo2%aza = geo2%aza (currpix, currline) 
    the_geo2%sca = geo2%sca (currpix, currline) 
    the_geo2%lon = geo2%lon (currpix, currline) 
    the_geo2%lat = geo2%lat (currpix, currline) 
    the_geo2%surfalt = geo2%height (currpix, currline) 
    the_geo2%lons(1:4) = geo2%clon(1:4,currpix, currline)
    the_geo1%lons(5) = geo2%lon(currpix, currline)
    the_geo2%lats(1:4) = geo2%clat(1:4,currpix, currline)
    the_geo1%lats(5) = geo2%lat(currpix, currline)
    the_geo2%elon(1:2) = geo2%elon(currpix-1:currpix, currline)
    the_geo2%elat(1:2) = geo2%elat(currpix-1:currpix, currline)



    the_year  = GranuleYear
    the_month = GranuleMonth
    the_day   = GranuleDay
    the_jday  = GranuleJday
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

    ! Check cloud fractions
    IF (which_cld /= 2) THEN
      the_cld_flg = geo1%cloud_qflg(currpix,currline)
      the_ai      = geo1%ai    (currpix, currloop)
      IF (the_cld_flg /= 10) THEN  ! Bad clouds for 10
        the_cfrac = geo1%cfr(currpix, currline)
        the_ctp   = geo1%ctp(currpix, currline)
      ELSE
        !the_cfrac = 0.0
        the_ctp = 0.0
      ENDIF
    ENDIF
    !-----------------------------------------------------------------------
    ! load soft calibration spectra
    !-----------------------------------------------------------------------
    IF (first .AND. biascorr) THEN
      WRITE(msg, *) TRIM(ADJUSTL(biasfname))//',which_biascorr=',which_biascorr
      errstat = OMI_SMF_setmsg (OMI_W_GENERAL, TRIM(msg), modulename, 0)
      IF ( which_biascorr == 7) THEN
        OPEN(UNIT=calunit, FILE=TRIM(ADJUSTL(biasfname)), STATUS='OLD', IOSTAT=errstat)
        IF ( errstat /= pge_errstat_ok ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
               TRIM(ADJUSTL(biasfname)), modulename, 0)
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
        soft%xwcorr = 1.0 ! Initialize to one
        DO is = 1, mswath
          READ (calunit, *) soft%nxcorr(is), soft%nxwav(is)
          DO iw = 1, soft%nxwav(is)
            READ (calunit, *) soft%xwavs(is, iw), soft%xwcorr(is, 1:soft%nxcorr(is), iw)
          ENDDO
        ENDDO
        CLOSE(UNIT=calunit) 
      ENDIF

      first = .FALSE.
    ENDIF

    !------------------------------------------------------------------------------------------
    ! Radiance Spectrum
    !----------------------------------------------------------------------------------------
    n_rad_wvl = rad%nwav(currpix, currloop) 
    div_rad   = rad%norm(currpix, currloop)

    curr_rad_spec(wvl_idx, 1:n_rad_wvl) = rad%wavl(1:n_rad_wvl, currpix, currloop) 
    curr_rad_spec(spc_idx, 1:n_rad_wvl) = rad%spec(1:n_rad_wvl, currpix, currloop)   
    IF (use_meas_sig) THEN
      curr_rad_spec(sig_idx, 1:n_rad_wvl) = rad%prec(1:n_rad_wvl, currpix, currloop)
    ELSE
      curr_rad_spec(sig_idx, 1:n_rad_wvl) = normweight
    ENDIF
    nradpix(1:numwin) = rad%npix(1:numwin, currpix, currloop)     ! Solar Spectrum
    n_irrad_wvl = n_rad_wvl

    IF ( biascorr .AND. which_biascorr == 7 ) THEN
      fidx = 1
      DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        ch = band_selectors(i) 
        CALL INTERPOL(soft%xwavs(ch, 1:soft%nxwav(ch)), soft%xwcorr(ch, currpix, 1:soft%nxwav(ch)), soft%nxwav(ch), &
             curr_rad_spec(wvl_idx, fidx:lidx),  corr(1:nradpix(i)), nradpix(i), errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
          errstat = pge_errstat_error
          RETURN
        ENDIF
        curr_rad_spec(spc_idx, fidx:lidx) = curr_rad_spec(spc_idx, fidx:lidx) / corr(1:nradpix(i))

        fidx = lidx + 1
      ENDDO
    ENDIF

    curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = irrad%wavl(rad%wind(1:n_rad_wvl, currpix, currloop), currpix) 
    curr_sol_spec(spc_idx, 1:n_irrad_wvl) = irrad%spec(rad%wind(1:n_rad_wvl, currpix, currloop), currpix) 
    IF (use_meas_sig) THEN
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = irrad%prec(rad%wind(1:n_rad_wvl, currpix, currloop), currpix)
    ELSE
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = normweight
    ENDIF
    nsolpix(1:numwin) = nradpix(1:numwin)  
    
    !-------------------------------------------------------------------------------------------------
    ! Reflectance spectrum at ~370 nm
    !----------------------------------------------------------------------------------------
    rad_posr (1:nrefl) = refl%radwavl(1:nrefl, currpix, currloop)
    rad_specr(1:nrefl) = refl%radspec(1:nrefl, currpix, currloop)

    IF ( biascorr ) THEN

      IF ( which_biascorr == 7 ) THEN
        IF ( soft%xwavs(mswath, soft%nxwav(mswath)) < rad_posr(1) ) THEN
          rad_specr(1:nrefl) = rad_specr(1:nrefl) / soft%xwcorr(mswath, currpix, soft%nxwav(mswath))       
        ELSE
          fidx = MINVAL( MINLOC( soft%xwavs(mswath, 1:soft%nxwav(mswath)), MASK = &
               (soft%xwavs(mswath, 1:soft%nxwav(mswath)) > rad_posr(1) )))
          lidx = MINVAL(MAXLOC( soft%xwavs(mswath, 1:soft%nxwav(mswath)), MASK = &
               (soft%xwavs(mswath, 1:soft%nxwav(mswath)) < rad_posr(nrefl) )))
          IF (fidx > lidx) THEN
            idum = fidx
            fidx = lidx
            lidx = idum
          ENDIF
          rad_specr(1:nrefl) = rad_specr(1:nrefl) / &
               (SUM(soft%xwcorr(mswath, currpix, fidx:lidx)) / (lidx - fidx + 1))
        ENDIF
      ENDIF
    ENDIF


    IF (which_cld == 2 .OR. (the_cld_flg == 10 .AND. which_cld >= 3 ))  THEN
      the_cld_flg = 2  ! ISCCP
      CALL GET_TOMSV8_CTP(the_month, the_day, the_lon, the_lat, the_ctp, pge_error_status)
      the_ai = -999.0
    ENDIF
    ! check for NAN
    IF (IEOR(IBCLR(TRANSFER(the_cfrac, NAN), DPSB), NAN) == 0 .OR. &
        the_cfrac < 0 ) THEN
        the_cfrac = 0.5
    ENDIF
    has_glint = .FALSE.
    glintprob = 0.0
    ! Land-water flag=1: >=8 not used, else contain water
    IF  (the_landwater_flg /= 1 .AND. the_landwater_flg < 8) THEN
      IF (the_glint_flg == 1) THEN
        has_glint = .TRUE.
        CALL SUNGLINT_PROBABILITY (the_sza_atm, the_vza_atm, the_aza_atm, glintprob)
        !PRINT *, 'Glint Probability: ', glintprob, the_cfrac
        !IF (the_cfrac < 0.30 * glintprob) the_cfrac = 0.0
      ENDIF
    ENDIF

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
      OPEN(UNIT=calunit, FILE='INP/sim.inp', STATUS='unknown')
      READ(calunit, *) the_sza_atm, the_vza_atm, the_aza_atm, the_fixalb, the_surfalt, &
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

      CLOSE (calunit)
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
    CALL adj_rad_sig (curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl),  &
         curr_sol_spec(wvl_idx:sig_idx, 1:n_rad_wvl), n_rad_wvl)

    ! Make sure that reference spectra has  more wavelengths than
    ! irradiance interpolation and shifting      
    curr_rad_spec_save = curr_rad_spec ! used in prepare_database
    fidx = 1
    DO i = 1, numwin
      lidx = fidx + nradpix(i)
      WRITE(www_lun, '(A,i2,f8.4, A,f8.2,A,f8.2)') & 
      'shift rad/sol:',i, curr_rad_spec(1, fidx:fidx)-curr_sol_spec(1, fidx:fidx), &
      'sampleing rate:', curr_rad_spec(1, fidx+1:fidx+1)-curr_rad_spec(1,fidx:fidx), &
      'shift rad/sol/sample rate:',(curr_rad_spec(1,2)-curr_sol_spec(1,2))/(curr_rad_spec(1,2)-curr_rad_spec(1,1))
      fidx = lidx + 1
    ENDDO

    nhtrunc = radnhtrunc
    ntrunc = nhtrunc * 2
    ntrunc1 = ntrunc + 1
    fidx = 1
    DO i = 1, numwin
      lidx = fidx + nradpix(i) - ntrunc1
      curr_rad_spec(1:sig_idx, fidx:lidx) = curr_rad_spec(1:sig_idx, fidx + nhtrunc : lidx + nhtrunc)
      IF (lidx  < n_rad_wvl - ntrunc1 ) THEN
        curr_rad_spec(1:sig_idx, lidx+1:n_rad_wvl - ntrunc) =  &
        curr_rad_spec(1:sig_idx, lidx+ntrunc1:n_rad_wvl)
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
    IF (theline > 0) THEN
      IF (rad%nwav(currpix, currloop) /= rad%nwav(currpix, currloop-1)) THEN
        redo_database = .TRUE.
      ELSE 
        IF (ANY(rad%wind(1:n_rad_wvl, currpix, currloop) - &
             rad%wind(1:n_rad_wvl, currpix, currloop-1) /= 0)) redo_database = .TRUE.
      ENDIF
    ENDIF
    IF ( MOD (theline, radwavcal_freq) == 0 .OR. redo_database) THEN 
      ! --------------------------------------------------------------
      ! Spline data bases, compute undersampling spectrum, and prepare
      ! reference spectra for fitting.
      ! --------------------------------------------------------------
      CALL prepare_databases ( n_rad_wvl, curr_rad_spec(wvl_idx,1:n_rad_wvl), pge_error_status )
      IF ( pge_error_status >= pge_errstat_error ) THEN 
           redo_database = .TRUE.
           WRITE(*,'(A)') 'Errors in prepare_database'
           RETURN
      ENDIF
      WRITE(www_lun,'(A)')  'common mode load for future tempo, here' 
      ! CALL load_comres(pge_error_status)
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
      IF ( which_biascorr == 7) THEN
        !! Kai need to comment out gascorr, i.e., the following two lines to make the code work for arbitrary window. 
        !curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) * gascorr(currpix, 1:n_rad_wvl, 2)
        !curr_rad_spec(spc_idx, 1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl) * &
        !     (1.0d0 + xw2corr(currpix, 1:n_rad_wvl, 2) / 100.)
      ENDIF
    ENDIF


    ! Initialized fitted variables from valid western and southern neighbors 
    ! since the retrievals are performed from west to east and from south to north
    west_idx  = currpix  - 1
    south_idx = currloop - 1
    IF (south_idx < 0) south_idx = nlines_max - 1

    o3p%initval(currpix, currloop) = 0
    fitvar = 0.0
    finit = 0.0
    IF (west_idx > 0) THEN
      IF (o3p%exitval(west_idx, currloop) > 0) THEN  ! Western pixel (success retrieval)
        fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) + &
             o3p%fitvar(west_idx, currloop, 1:n_fitvar_rad)
        finit = finit + 1.0
      ENDIF

      IF (south_idx >= 0 .AND. south_idx /= nlines_max - 1) THEN
        IF (o3p%exitval(west_idx, south_idx) > 0) THEN ! Southwestern pixel (success retrieval)
          fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) &
               + o3p%fitvar(west_idx, south_idx, 1:n_fitvar_rad) * 0.5
          finit = finit + 0.5
        ENDIF
      ENDIF
    ENDIF

    IF ( south_idx >= 0 ) THEN
      IF (o3p%exitval(currpix, south_idx) > 0) THEN     ! Southern pixel (success retrieval)
        fitvar(1:n_fitvar_rad) = fitvar(1:n_fitvar_rad) + &
             o3p%fitvar(currpix, south_idx, 1:n_fitvar_rad) 
        finit = finit + 1.0
      ENDIF
    ENDIF

    IF (finit > 0) THEN
      fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad)) = fitvar(1:n_fitvar_rad) / finit
      o3p%initval(currpix, currloop) = 1
    ENDIF

    RETURN
  END SUBROUTINE adj_earthshine_data


  SUBROUTINE adj_rad_sig (radspec, solspec, np)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx
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


    RETURN
  END SUBROUTINE adj_rad_sig

  SUBROUTINE load_comres(errstat)
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: max_fit_pts, maxchlen
    USE OMSAO_variables_module, ONLY: refdbdir, n_refwvl, database, &
         nradpix, refidx, currpix
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

  END SUBROUTINE load_comres

END MODULE tmpo_adj_data
