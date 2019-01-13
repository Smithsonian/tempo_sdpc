!
module omi_adj_data
   
  USE OMSAO_omidata_module, ONLY: nxtrack_max,nlines_max, nfxtrack, & 
      ring=>omi_ring, refl=>omi_refl,irrad=>omi_irrad, &
      cali=>omi_cali, rad=>omi_rad, geo=>omi_geo, o3p=>omi_o3p
  USE OMSAO_pixelcorner_module, ONLY: omi_NSPC
  public adj_solar_data, adj_earthshine_data
  private adj_rad_sig, load_omi_comres
contains

  SUBROUTINE adj_solar_data (pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx!, hwe_idx
    USE OMSAO_parameters_module, ONLY: normweight!, mrefl
    USE OMSAO_variables_module,  ONLY: curr_sol_spec, n_irrad_wvl, &
         use_meas_sig, numwin, nsol_ring, sol_spec_ring, nsolpix, &
         yn_varyslit, slit_rad, solwinfit, nslit, slitwav, slitfit, &
         sring_fidx, sring_lidx,  currpix, which_slit!, &
         !reduce_resolution, slitdis
    USE ozprof_data_module,      ONLY: div_sun, sun_posr, sun_specr, nrefl
    USE OMSAO_errstat_module 

    IMPLICIT NONE

    ! =================
    ! Output variables
    ! =================
    INTEGER, INTENT (OUT) :: pge_error_status

    !INTEGER :: i

    pge_error_status = pge_errstat_ok

    ! Solar spectrum
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
      END IF
    ENDIF

    RETURN
  END SUBROUTINE adj_solar_data


  SUBROUTINE adj_earthshine_data (theline, pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: wvl_idx, spc_idx, sig_idx, &
         n_max_fitpars, solar_idx!, rsl_idx, fsl_idx, comm_idx, com1_idx
    USE OMSAO_parameters_module, ONLY: mswath, normweight, max_fit_pts,max_spec_pts, &
         maxchlen, l1l2inp_unit
    USE OMSAO_variables_module,  ONLY: curr_rad_spec,curr_sol_spec,curr_rad_spec_save, &
         n_rad_wvl, use_meas_sig, numwin, nradpix, the_sza_atm, the_vza_atm, &
         the_aza_atm, the_sca_atm, the_month, the_year, the_day, the_lon, &
         the_lat, the_lats, the_lons, edgelons, edgelats, the_surfalt, nview, &
         nloc, the_utc, n_radwvl_sav, radwvl_sav,  nradpix_sav, saa_minlat, &
         saa_maxlat, saa_minlon, saa_maxlon, saa_minlat1, saa_maxlat1, &
         saa_minlon1, saa_maxlon1, do_bandavg, refidx, fitvar_rad_saved, &
         n_fitvar_rad, radwavcal_freq, currpix, currloop, currline, &
         n_irrad_wvl, nsolpix, actspec_rad, database, band_selectors, &
         mask_fitvar_rad, radnhtrunc, refnhextra, &
         tabdir, orbnumsol, NSPC_omi!, refwvl, refidx_sav, lo_radbnd, up_radbnd, i0sav
    !USE OMSAO_pixelcorner_module, ONLY: & 
    !     omi_allclon, omi_allclat,  omi_allelon, omi_allelat,&
    !     omi_allHeight, omi_alltime,omi_alllat, omi_alllon, &
    !     omi_allsza, omi_allvza, omi_allaza, omi_allsca, omi_allglint_flg, omi_allland_water_flg, omi_allsnow_ice_flg
    USE OMSAO_omicloud_module, ONLY: OMIL2_clouds 
    USE ozprof_data_module, ONLY: div_rad, div_sun, rad_posr, rad_specr, &
         nsaa_spike, saa_flag, the_cfrac, the_ctp, the_cld_flg, which_cld, &
         the_cod, the_orig_cfr, the_orig_ctp, scacld_initcod, the_orig_cod, &
         the_ai, radcalwrt, biasfname, biascorr, &
         which_biascorr, nrefl, aerosol, which_aerosol, scale_aod, &
         scaled_aod, do_simu, the_fixalb, do_lambcld, lambcld_refl, &
         has_glint, glintprob, sun_posr, sun_specr, pos_alb, &
         the_snowice, the_landwater_flg, the_glint_flg !, ozprof_start_index, ozprof_end_index, nir
    USE OMSAO_errstat_module 
    use m_angle_sat2toa, only: sunglint_probability !, adjust_angle
    use m_avg_band, only: avg_band_radspec
    use m_ezspline_interpolation, only: bspline, interpol
    use m_get_cloud, only: get_tomsv8_ctp
    use m_prepare_databases
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
    LOGICAL                     :: redo_database
    INTEGER :: hour, minute, fidx, lidx, i, j, west_idx, south_idx, idxoff, &
         nhtrunc, ntrunc, ntrunc1, errstat, ntempx, nch, ix, nord, ch, nw, &
         is, nsub, idum, iw
    INTEGER (KIND=i4)           :: estat
    REAL (KIND=r8)              :: second, finit
    REAL (KIND=dp), DIMENSION (n_max_fitpars) :: fitvar

    ! xliu (02/03/2007): variables for correcting across-track dependent biases
    INTEGER, PARAMETER  :: maxord = 12
    INTEGER, DIMENSION (mswath), SAVE :: corr_npars, nxcorr, nxwav
    REAL (KIND=dp), DIMENSION(mswath), SAVE    :: corr_woffset
    REAL (KIND=dp), DIMENSION(:,:,:), SAVE, POINTER  :: & ! (mswath,nx,0:maxord)
                                      corrpars, offset_pars, slope_pars
    REAL (KIND=dp), DIMENSION(:), POINTER       :: corr, offset, slope !max_spec_pts
    REAL (KIND=dp), DIMENSION(:,:), SAVE, POINTER :: &
                                                  allcorr, alloffset ! nx, max_spec_pts
    REAL (KIND=dp), DIMENSION(:,:), SAVE, POINTER  :: xcorr
    REAL (KIND=dp), DIMENSION(:,:,:),SAVE,POINTER ::& ! mswath, nx, max_spec_pts 
                                                xwcorr,  xwslp, xwoff
    REAL (KIND=dp), DIMENSION(:,:), SAVE, POINTER    :: xwavs ! (mswath,max_spec_pts)
    INTEGER, SAVE                          :: nxgascorr, nxw2corr
    REAL (KIND=dp), DIMENSION(:,:,:),SAVE,POINTER  :: gascorr, xw2corr ! nx, max_spec_pts, 2
    INTEGER, DIMENSION(:,:),SAVE, POINTER  :: gascorr_npts, xw2corr_npts ! nx,3
    REAL (KIND=dp), DIMENSION (:, :), POINTER :: del ! ( 1:maxord, max_spec_pts) 

    REAL (KIND=dp)               :: woffset, rad347, irad347
    real (kind=dp), dimension(1) :: temp_pos_alb, temp_rad
    CHARACTER (LEN=maxchlen) :: gascorr_fname, xw2corr_fname
    CHARACTER (LEN=255)      :: msg !! Kai
    LOGICAL, SAVE   :: first = .TRUE.
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
    nview       = 1
    the_sza_atm = geo%sza    (currpix, currline)
    the_vza_atm = geo%vza    (currpix, currline)
    the_aza_atm = geo%aza    (currpix, currline)
    the_sca_atm = geo%sca    (currpix, currline)
    the_lon     = geo%lon    (currpix, currline)
    the_lat     = geo%lat    (currpix, currline)
    the_surfalt = geo%Height (currpix, currline) / 1000.
    nloc          = 5
    the_lons(1:4) = geo%clon  (1:4,currpix, currline)
    the_lons(5)   = the_lon
    the_lats(1:4) = geo%clat  (1:4,currpix, currline)
    the_lats(5)   = the_lat
    edgelons(1:2) = geo%elon  (currpix-1:currpix, currline)
    edgelats(1:2) = geo%elat  (currpix-1:currpix, currline)
    the_snowice   = geo%snow_ice_flg(currpix, currline)
    the_landwater_flg = geo%land_water_flg(currpix,currline)
    the_glint_flg = geo%glint_flg(currpix, currline)
    NSPC_omi = omi_NSPC(currline)
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
 
    estat = PGS_TD_TAItoUTC(geo%time(currline), the_utc)
    READ (the_utc, '(I4, 1x, I2, 1x, I2, 1x, I2, 1x, I2, 1x, F9.6)') &
         the_year, the_month, the_day, hour, minute, second

    IF (first .AND. biascorr) THEN

     allocate ( del(1:maxord, max_spec_pts))  
     allocate ( corr(max_spec_pts),offset(max_spec_pts),slope(max_spec_pts))
     IF (which_biascorr == 1) THEN
        allocate (allcorr(nxtrack_max, max_spec_pts))
        allocate (alloffset(nxtrack_max, max_spec_pts))
     ENDIF
     IF (which_biascorr == 2) THEN 
       allocate ( corrpars(mswath, nxtrack_max, 0:maxord))
     ENDIF
     IF (which_biascorr == 4) THEN 
       allocate ( offset_pars(mswath, nxtrack_max, 0:maxord))
       allocate ( slope_pars(mswath, nxtrack_max, 0:maxord))
     ENDIF
     IF (which_biascorr == 7) THEN
       allocate ( xwavs (mswath, max_spec_pts))
       allocate ( xcorr (mswath, nxtrack_max))
       allocate ( xwcorr(mswath, nxtrack_max, max_spec_pts))
       allocate ( gascorr_npts(nxtrack_max, 3), xw2corr_npts(nxtrack_max, 3))
       allocate ( gascorr(nxtrack_max, max_spec_pts, 2))
       allocate ( xw2corr (nxtrack_max, max_spec_pts, 2))
     ENDIF
     IF (which_biascorr == 8 .or. which_biascorr == 9) THEN 
       allocate ( xwavs (mswath, max_spec_pts))
       allocate ( xwslp(mswath, nxtrack_max, max_spec_pts))
       allocate ( xwoff(mswath, nxtrack_max, max_spec_pts))
     ENDIF

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
    n_rad_wvl = rad%nwav(currpix, currloop) 
    div_rad   = rad%norm (currpix, currloop)

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


    curr_sol_spec(wvl_idx, 1:n_irrad_wvl) = irrad%wavl(rad%wind(1:n_rad_wvl, currpix, currloop), currpix) 
    curr_sol_spec(spc_idx, 1:n_irrad_wvl) = irrad%spec(rad%wind(1:n_rad_wvl, currpix, currloop), currpix) 
    IF (use_meas_sig .AND. orbnumsol /= 99999) THEN
      curr_sol_spec(sig_idx, 1:n_irrad_wvl) = irrad%prec(rad%wind(1:n_rad_wvl, currpix, currloop), currpix)
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
    rad_posr (1:nrefl) = refl%radwavl(1:nrefl, currpix, currloop)
    rad_specr(1:nrefl) = refl%radspec(1:nrefl, currpix, currloop)
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
   ! print * , 'PLEASE REMOVE AFTER TESTING:ctp and cfrac is fixed to'
   ! the_cfrac =0.086421
   ! the_ctp   = 732.18
    ! Special treatments for sea glint
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
    CALL adj_rad_sig (curr_rad_spec(wvl_idx:sig_idx, 1:n_rad_wvl),  &
         curr_sol_spec(wvl_idx:sig_idx, 1:n_rad_wvl), n_rad_wvl)

    ! Make sure that reference spectra has  more wavelengths than
    ! irradiance interpolation and shifting      
    nhtrunc = radnhtrunc
    ntrunc = nhtrunc * 2
    ntrunc1 = ntrunc + 1
    fidx = 1
    curr_rad_spec_save = curr_rad_spec
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
  END SUBROUTINE adj_rad_sig

 
  SUBROUTINE load_omi_comres(errstat)
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: max_fit_pts, maxchlen
    USE OMSAO_variables_module, ONLY: refdbdir, n_refwvl, database, &
         nradpix, refidx, currpix!, up_radbnd, refwvl, numwin, &
         !n_fitvar_rad, fitvar_rad_str, mask_fitvar_rad, lo_radbnd
    USE OMSAO_indices_module, ONLY: com1_idx, com_idx, comfidx, cm1fidx
    USE OMSAO_errstat_module, ONLY: pge_errstat_error, pge_errstat_ok, www_lun

    INTEGER, INTENT (OUT)      :: errstat

    INTEGER, PARAMETER            :: nx = 30, lun = 12
    INTEGER                       :: ix, itemp, i, fidx, lidx
    CHARACTER (LEN=maxchlen)      :: comres_fname

    REAL (KIND=dp), DIMENSION (:,:,:), SAVE, POINTER :: comres
    INTEGER, DIMENSION(:,:),SAVE, POINTER :: npts
    LOGICAL,                            SAVE :: first = .TRUE.
    CHARACTER (LEN=15), PARAMETER :: modulename = 'load_omi_comres'

    errstat = pge_errstat_ok
    IF (first) THEN
      allocate (npts(nx, 3))
      allocate (comres(nx, max_fit_pts, 2))
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

end module omi_adj_data 
