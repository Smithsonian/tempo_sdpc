MODULE m_specfit_ozprof

  ! ***************************************************************
  ! Author:  Xiong Liu
  ! Date:    July 24, 2003
  ! Purpose: Ozone profile retrieval using PTR or OE
  ! 
  ! Modification History:
  ! 1.  xliu, Jan 19, 2004
  !     a. Remove using ELSUNC 
  !     b. Add a interface for preparing atmosphere using tomsv8 
  !        profile, EP total ozone, temperature profiles, surface 
  !        pressure
  !     c. Add a interface to prepare a priori covariance  
  !     d. Calculate error in total ozone here
  !     e. Move get_reg_matrix here
  !     f. Prepare measurment error and measurement vector
  ! jbak
  !     should be re-checked for setting up initial/apriori/covariance for
  !     spectral refpectance (03/28/2020)
  ! ***************************************************************

CONTAINS

SUBROUTINE specfit_ozprof (initval, fitcol, dfitcol, rms, exval)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: missing_value_dp, maxlay
  USE OMSAO_indices_module,      ONLY: instrument_idx, simu_idx, &
       n_max_fitpars, wvl_idx, spc_idx,         &
       sig_idx, maxalb, so2_idx, so2v_idx, us1_idx, &
       us2_idx, max_calfit_idx, max_rs_idx, mxs_idx, maxoth, shift_offset, maxwfc, &
       comvidx, cm1vidx, comfidx, cm1fidx, bro2_idx, o2o2_idx, o2_idx, o2t2_idx,   &
       h2o_idx, h2ot2_idx, cm2fidx, cm3fidx !bro_idx, hcho_idx, no2_t1_idx, no2_t2_idx
  USE OMSAO_variables_module,    ONLY: scnwrt, refspec_norm,&
       weight_rad, nradpix, npix_fitted, rad_wav_avg,&
       numwin, widx_vis, wcenter_uvvis, nviswin, &
       n_rad_wvl, curr_rad_spec,fitwavs, currspec, fitweights, &
       n_fitvar_rad, fitvar_rad_str, mask_fitvar_rad, lo_radbnd, up_radbnd, &
       fitvar_rad, fitvar_rad_saved, fitvar_rad_std, fitvar_rad_nstd, &
       fitvar_rad_apriori, fitvar_rad_aperror, fitvar_rad_init, fitvar_rad_init_saved, & 
       fitspec_rad, fitres_rad, fitspec_q, fitres_q, &
       currpix, currline, nloc,  the_lons, the_lats, the_surfalt, the_lon, the_lat, &
       the_year, the_day, the_month, the_jday
  USE ozprof_data_module,        ONLY: &
       use_oe, use_logstate, smooth_ozbc, ptr_order, ozwrtwf, & 
       use_large_so2_aperr, do_sy_diagonal, merr_corrlen, msyserr, merr_covar, &
       ozprof_start_index, ozprof_end_index, ozfit_start_index, ozfit_end_index, &
       covar, ncovar, contri, ozwrtcontri, weight_function, ndiv, &
       trace_profwf, trace_contri, trace_prof, & 
       ozprof_std, ozprof_nstd, ozprof_ap, ozprof_apstd, ozprof_init, ozprof, &
       start_layer, end_layer, nup2p, nlay, nflay, ntp, nsfc, &
       use_tropopause, which_toz, which_atm, ps0, pst, pst0, atmosprof, &
       nt_fit, tf_fidx, tf_lidx,t_fidx,t_lidx, &
       pos_alb, nalb, nfalb, eff_alb, eff_alb_init, albidx, albfidx,&
       nwfc, nfwfc, eff_wfc, eff_wfc_init, wfcidx, wfcfidx, &
       the_cfrac, the_cod, the_ctp, &
       maxawin, actawin, aerwavs,  &
       ngas, gasidxs, fgasidxs, tracegas, mgasprof, &
       do_subfit, fit_atanring, nos, nsh, & !nsl, nrn,nis, ndc
       osind, rnind, dcind, isind, irind, shind, shfind, osfind, & !slind
       np1, np2, np3, p1find, p2find, p3find, p1ind, p2ind, p3ind,ncm, cmind, &
       ozabs_convl, use_effcrs, radcalwrt, do_simu, &
       glintprob, the_snowice, the_landwater_flg, the_landfrac, &
       tropaod, tropsca, tropwaer, strataod, stratsca, taodind, taodfind, twaeind,   &
       saodfind, ecfrind, ecfrfind, ecodind, ecodfind, ectpind, ectpfind, has_glint, &
       twaefind, saodind, sprsind, sprsfind, so2zind, so2zfind,&
       albfpix, alblpix, is_albspcvar, use_albeofs, nalbspc, nactalbspc, &
       which_albspc,sfcalbs, albspcs, do_brdf, nalbwf, do_debug_o3p, which_alb
  USE OMSAO_errstat_module
  USE m_get_reg_matrix, ONLY: get_reg_matrix
  USE m_get_bclayer, ONLY: get_bc_layer
  USE m_get_toz,     ONLY:get_toz
  USE m_get_o3prof,  ONLY: get_apriori_covar
  USE m_set_tracegas,ONLY: set_tracegas
  USE m_make_atm,    ONLY: make_atm
  USE m_set_cldalb,  ONLY: set_cldalb, oceanflg, snowflg
  USE m_utilities,   ONLY: day_of_year
  USE m_set_brdf,    ONLY: Surface
  USE m_lidort_util, ONLY: get_hres_radcal_waves
  USE m_ozprof_inverse, ONLY: ozprof_inverse

  IMPLICIT NONE

  ! =============================
  !  Input / Output Variables
  ! =============================
  INTEGER,            INTENT(IN)  :: initval
  INTEGER,            INTENT(OUT) :: exval
  REAL     (KIND=dp), INTENT(OUT) :: rms
  REAL     (KIND=dp), DIMENSION(3), INTENT(OUT)    :: fitcol
  REAL     (KIND=dp), DIMENSION(3, 2), INTENT(OUT) :: dfitcol  ! smooth+noise, noise

  ! ===============
  ! Local variables
  ! ===============
  INTEGER  :: i, j, k, k1, is, nump, errstat, npoints, u1idx, u2idx, nsub, nord,fidx, lidx, unit_deg, off
  REAL (KIND=dp) :: asum, ssum, chisq, tmpsa, aodscl, waerscl, salbedo, wavavg
  REAL (KIND=dp) :: toz, albfc_aperr, albfc_aperr1, albfc_aperr2
  REAL (KIND=dp), DIMENSION(n_max_fitpars, n_max_fitpars) :: bb, sa
  REAL (KIND=dp), DIMENSION(nlay, nlay)                   :: sao3
  REAL (KIND=dp), DIMENSION (n_max_fitpars)               :: lowbond, upbond, fitvar, &
       stderr, fitvarap, stderr1
  REAL (KIND=dp), DIMENSION(maxlay)                       :: ozprof_std_sav
  REAL (KIND=dp), DIMENSION(n_max_fitpars, n_max_fitpars) :: covar_sav
  CHARACTER (LEN=6), DIMENSION(n_max_fitpars)             :: varname

  LOGICAL, SAVE :: first = .TRUE.
  INTEGER, SAVE :: ozp_fidx,  ozp_lidx, ozf_fidx, ozf_lidx, nf
  
  !xliu: 09/03/05, add sacldscal, scaling factor for scaling a priori covariance below clouds
  REAL (KIND=dp), DIMENSION(nlay)           :: sacldscl
  REAL (KIND=dp), ALLOCATABLE               :: tmpcovar(:,:), tmpspc(:)
  REAL (KIND=dp) :: init_alb
 ! systematic noise
  INTEGER, PARAMETER              :: nreg = 4 !should be updated for TEMPO/GOME2
  REAL (KIND=dp), DIMENSION(nreg) :: reg_noise =  &
       (/0.004, 0.004, 0.002, 0.002/)
  REAL (KIND=dp), DIMENSION(0:nreg) :: reg_waves = &
       (/260.0, 300.0, 310.0, 380., 800./)

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  LOGICAL :: do_debug_out = .false.
  CHARACTER (LEN=14), PARAMETER :: modulename = 'specfit_ozprof'

  ! Initialize variables for convenience
  toz = 0.0
  nalbwf  = 1 ! become more if do_brdf
  npoints = n_rad_wvl
  errstat = pge_errstat_ok
  IF (first) THEN  ! only need to be initialized once
    nf = n_fitvar_rad;   fitvar_rad = 0.0
    ozp_fidx = ozprof_start_index; ozp_lidx = ozprof_end_index
    ozf_fidx = ozfit_start_index; ozf_lidx = ozfit_end_index
    sa = 0.0; fitvar_rad_apriori = 0.0;   fitvar_rad_std = 0.0;  fitvar_rad_nstd = 0.0
    ozprof_std = 0.0;  ozprof_nstd = 0.0  

    fitvar = 0.0 ; lowbond = 0.0 ; upbond = 0.0 
    varname(1:nf) = fitvar_rad_str(mask_fitvar_rad(1:nf))
    the_jday = day_of_year(the_year, the_month, the_day)

    IF (.NOT. use_effcrs) THEN
      CALL get_hres_radcal_waves(errstat)
      IF (errstat == pge_errstat_error) THEN
        WRITE(*, *) modulename, ': Errors in getting fine & radiance calculation wavelength grids!!!'
        STOP 1
      ENDIF
    ENDIF
    CALL set_tracegas (ngas)
    first = .FALSE.
  ENDIF

  ! use previous fitting results except T, albedo, cloud will be updated
  ! use previous ozone will speed the convergence (could even double)
  fitvar_rad_init = fitvar_rad_saved

  ! ===================================================================
  !	         Set up measurement vector and measurement error
  ! ===================================================================
  fitwavs   (1:npoints) = curr_rad_spec(wvl_idx,1:npoints)
  currspec  (1:npoints) = curr_rad_spec(spc_idx,1:npoints)
  fitweights(1:npoints) = curr_rad_spec(sig_idx,1:npoints)

  widx_vis = npoints + 1
  IF (nviswin > 0) THEN
      widx_vis = MINVAL(MAXLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) <wcenter_uvvis)))+1
  ENDIF
 ! constract Sy, Sy = S_rnd + S_sys, S_rnd is diagonal, represented by
 ! measurements noise,
 ! S_sys is square, symetric, represented by sig^2*exp(-(lambda1-lambda2)/h)  
  do_sy_diagonal = .true.
  IF (.NOT.do_sy_diagonal) THEN
    IF (allocated(msyserr)) deallocate(msyserr, merr_covar)
    allocate (msyserr(npoints), merr_covar(npoints, npoints))
    msyserr = 0.002D0; merr_corrlen = 100.0D0
    fidx = 1; i = fidx
    DO j = 1, nreg
      DO WHILE (i <= npoints )
        IF (fitwavs(i) < reg_waves(j)) THEN
          i = i + 1
        ELSE
          EXIT
        ENDIF
      ENDDO
      lidx = i - 1
      msyserr(fidx:lidx) = reg_noise(j)
      fidx = lidx + 1
    ENDDO
    !msyserr = 0.01D0
    ALLOCATE(tmpcovar(npoints, npoints))
    tmpcovar = 0.0D0
    DO i = 1, npoints
      tmpcovar(i, i) = fitweights(i)**2 + msyserr(i)**2
    ENDDO
    DO i = 1, npoints
      DO j = 1, i-1
        tmpcovar(i, j) = msyserr(i)*msyserr(j) * EXP(-(ABS(fitwavs(i)-fitwavs(j))/merr_corrlen))
        !tmpcovar(i, j) =  sqrt(tmpcovar(i, i)* tmpcovar(j, j)) * &
        !                 EXP(-ABS((fitwavs(i)-fitwavs(j))/merr_corrlen))
        tmpcovar(j, i) = tmpcovar(i, j)
     ENDDO
    ENDDO

    ! UV1 and UV2 does not correlate
    merr_covar = 0.0D0
    fidx = 1
    DO i = 1, numwin
      lidx = fidx + nradpix(i) - 1
      merr_covar(fidx:lidx, fidx:lidx) = tmpcovar(fidx:lidx, fidx:lidx)
      fidx = lidx + 1
    ENDDO
    merr_covar(1:npoints,1:npoints)=tmpcovar(1:npoints,1:npoints)
    DEALLOCATE(tmpcovar)
  ENDIF

  IF (ozabs_convl) THEN
    ! For aerosol properties
    actawin = numwin + 2
    IF (actawin > maxawin) STOP 'Need to increase maxawin!!!'
    fidx = 1
    DO i = 1, numwin
       lidx = fidx + nradpix(i) - 1
       aerwavs(i) = fitwavs(fidx)
       fidx = lidx + 1
    ENDDO
    aerwavs(numwin + 1) = fitwavs(lidx)
    i = COUNT(mask = (aerwavs(1:numwin+1) < pos_alb))
    aerwavs(i+2:numwin+2) = aerwavs(i+1:numwin+1)
    aerwavs(i+1) = pos_alb 

    ! calculate approximate average wavelength for the window
    IF ( weight_rad ) THEN
       asum = SUM ( fitwavs(1:n_rad_wvl) / fitweights(1:n_rad_wvl)**2 )
       ssum = SUM ( 1.D0 / fitweights(1:n_rad_wvl)**2 )
       rad_wav_avg = asum / ssum
    ELSE
       rad_wav_avg = (fitwavs(n_rad_wvl) + fitwavs(1)) / 2.0
    END IF
  ENDIF
  
  ! =======================================================================
  !       Set up atmospheric cloud properties, albedo and atmosphere
  ! ======================================================================
  IF (instrument_idx /= simu_idx) THEN 
    CALL get_bc_layer (which_atm, nloc, the_lons, the_lats, ps0, pst, the_surfalt)
    IF (.NOT. use_tropopause) pst = pst0
  ENDIF
  
  ! ====================================================================
  !	                 Set up atmospheric profiles
  ! ====================================================================
  ! 1. Monthly mean total ozone (EP)
  ! 2. Surface and Tropopause Pressure (also used for getting albedo)
  ! 3. Aerosols (SAGE + GOCART)
  ! 4. Clouds   (GOMECAT)  
  ! 5. Albedo Database (If use clouds)
  ! 6. Temperature (ECMWF)
  ! 7. Ozone, BrO, SO2, NO2, HCHO
  ! 8. Atmos. Profiles
  !    a. TOMS V7 with TOMSV8 temperature profiles (T, P, h, O3) (deleted)
  !    b. GOME workgroup a priori (O3, P, h, T.)  (deleted)
  !    c. ECMWF Temperature and TOMS V8/McPeters Profiles (best conditions)
  ! 9. Get apriori covariance (ozone and non-ozone paramters)
  ! 10.Get initial albedo
  ! ======================================================================
  ! note: the returned ozprof,atmosprof and nup2p are counted bottom up
  ! set up tracegas index
  IF (which_toz /= 0) CALL get_toz(which_toz, toz)
  CALL make_atm(the_year, the_month, the_day, ndiv, &
       the_cod, the_cfrac, the_ctp, nlay, toz, ps0, pst, atmosprof(:,0:nlay),    &
       ozprof(1:nlay), nup2p(0:nlay), sacldscl, errstat)  
  IF (errstat == pge_errstat_error)  THEN
     exval = -2; RETURN
  ENDIF
  ! second step, use channel 1 retrieval results
  !IF (curr_rad_spec(1, 1) > 307.2) THEN
  !   CALL GET_FIRST_RETRIEVAL(the_month, the_lat, nlay, atmosprof(1, 0:nlay), &
  !        atmosprof(2, 0:nlay), atmosprof(2, 0:nlay), ozprof(1:nlay), sao3)
  !   STOP 'NO Two-Step Retrieval!!!'
  !ENDIF

  ! ======================================================================
  !	 Set up initial value
  ! ======================================================================
  ! Add ozone and trace gases into initialized array for the first retrieval
  ! Then use the previous retrieval as the initial
  ! NO2 and HCHO: always using GEOS-CHEM fields with 100% error
  ! For BrO:      a prioir from model fields but with 1.0E-14 error globally (enough information)
  ! For SO2:      a priori from model fields but with dynamic a priori error 
  !               to deal with volcanic eruption (implemented in ozone_reverse.f90)

  IF (initval == 0 .OR. ANY(fitvar_rad_init(ozp_fidx:ozp_lidx) <= 0.0))  THEN
     fitvar_rad_init(ozp_fidx:ozp_lidx) = ozprof(1:nlay)  
  ENDIF
  IF (nsfc < nlay) fitvar_rad_init(ozp_lidx + 1 - nlay + nsfc : ozp_lidx) = ozprof(nsfc+1:nlay)
  fitvar_rad_apriori(ozp_fidx:ozp_lidx) = ozprof(1:nlay)
   
  ! Always use climatological temperature profiles
  fitvar_rad_init(t_fidx:t_lidx) = &
       (atmosprof(3, 1:nlay) + atmosprof(3, 0:nlay-1)) / 2.0
  ! For aerosols, a priori is based on model, a priori error is assumed as 100%, 
  ! and initial value is from previous retrievals.
  fitvar_rad_apriori(taodind) = tropaod(actawin)
  fitvar_rad_apriori(twaeind) = tropwaer(actawin)
  fitvar_rad_apriori(saodind) = strataod(actawin)  
  IF (initval == 0 .OR. taodfind == 0) fitvar_rad_init(taodind) = tropaod(actawin)
  IF (initval == 0 .OR. twaefind == 0) fitvar_rad_init(twaeind) = tropwaer(actawin)
  IF (initval == 0 .OR. saodfind == 0) fitvar_rad_init(saodind) = strataod(actawin)

  ! Set up albedo and cloud fraction in the retrieval
  ! albedo and cloud fraction can be adjusted based on 370.2 nm reflectance
  CALL set_cldalb(npoints, fitwavs, the_cod, the_ctp, the_cfrac, salbedo, errstat)
  
  IF (do_debug_o3p) WRITE(www_lun, '(A)') '@ finish set up cldalb'
  IF (errstat == pge_errstat_error)  THEN
     WRITE(*,'(A)') modulename//': Errors in set_cldalb'
     exval = -1; RETURN
  ENDIF
  ! For clouds, initial ctp, cod is based on assumed input (e.g., 20/10) or from other products, 
  ! which maybe re-adjusted using longer wavelengths

  fitvar_rad_init(ecfrind) = the_cfrac; fitvar_rad_apriori(ecfrind) = the_cfrac
  fitvar_rad_init(ecodind) = the_cod;   fitvar_rad_apriori(ecodind) = the_cod
  fitvar_rad_init(ectpind) = the_ctp;   fitvar_rad_apriori(ectpind) = the_ctp
  IF (so2zfind > 0) THEN
     fitvar_rad_init(so2zind) = fitvar_rad_init_saved(so2zind)
     lo_radbnd(so2zind) = the_surfalt; up_radbnd(so2zind) = 30.
     fitvar_rad_apriori(so2zind)  = fitvar_rad_init(so2zind)
  ENDIF

  ! For surface pressure
  fitvar_rad_init(sprsind) = ps0;       fitvar_rad_apriori(sprsind) = ps0 
  IF (sprsfind > 0) THEN
     lo_radbnd(sprsind) = atmosprof(1, nlay) - (atmosprof(1, nlay)-atmosprof(1, nlay-1)) * 0.5
     up_radbnd(sprsind) = atmosprof(1, nlay) + (atmosprof(1, nlay)-atmosprof(1, nlay-1)) * 0.5
  ENDIF

  DO i = albidx, albidx + nalb - 1
     k = i - albidx + 1
     READ (fitvar_rad_str(i)(4:5), '(I2)') nord
     IF (nord /= 0) fitvar_rad_init(i) = 0.D0
     IF (nord == 0) fitvar_rad_apriori(i) = fitvar_rad_init(i)
     IF (is_albspcvar(k) .AND. use_albeofs .AND. nord == 1 .and. nactalbspc == 1) THEN
        fitvar_rad_init(i) = 1.0d0
        fitvar_rad_apriori(i) = 1.0d0
     ELSE IF (is_albspcvar(k) .AND. use_albeofs .AND. nactalbspc > 1 .and. nord > 0 .AND. &
          which_albspc == 2 .AND. oceanflg /= 1 .AND. snowflg /= 1) THEN 
        fitvar_rad_init(i) = Surface%Option4%FitCoeff(nord, 1)
        fitvar_rad_apriori(i) = Surface%Option4%Fitcoeff(nord,1)
     ELSE IF (is_albspcvar(k) .AND. do_brdf .and. nord > 0) THEN
        fitvar_rad_init(i) = 1.00
        fitvar_rad_apriori(i) = 1.00
     ENDIF
    ! IF (nord > 0 .and. which_albspc ==2) &
     !print *, nord, is_albspcvar(k), fitvar_rad_str(i), fitvar_rad_init(i),fitvar_rad_init(i)*100/sqrt(Surface%Option4%FactorUncertainty(nord,nord))
  ENDDO
  DO i = wfcidx, wfcidx + nwfc - 1
     IF (fitvar_rad_str(i)(4:4) /= '0') fitvar_rad_init(i) = 0.0D0
  ENDDO

  ! For trace gases 
  DO k = 1, ngas
     i = fgasidxs(k)
     IF (i > 0) THEN
        j = mask_fitvar_rad(i)          
        ! tracegas(k, 8) = 1.0 - the_cfrac * tracegas(k, 8)  
        ! fitvar_rad_apriori(j) = mgasprof(k, nflay + 1) * refspec_norm(gasidxs(k)) !* tracegas(k, 8)
        ! This is incorrect, since everything is taken into account implicitly through weighting function       
        fitvar_rad_apriori(j) = mgasprof(k, nflay + 1) * refspec_norm(gasidxs(k)) 
        !IF ( npix_fitted == 0 .OR. gasidxs(k) == so2_idx .OR. gasidxs(k) == so2v_idx &
        !     .OR. initval == 0 .OR. fitvar_rad_init(j) < 0.0) &

        fitvar_rad_init(j) =  fitvar_rad_apriori(j) 
        IF (gasidxs(k) == so2_idx .OR. gasidxs(k) == so2v_idx) fitvar_rad_init(j) = 0.0 
        IF (gasidxs(k) == bro2_idx) fitvar_rad_apriori(j) = fitvar_rad_apriori(j) * 5. / 2.
        
        ! initial values for trace gas shift parameters needs to be fixed for every pixel
        ! otherwise increasing from North to South to unreasonably large
        j = shift_offset + gasidxs(k)
        fitvar_rad_init(j) = 0.00
     ENDIF
  ENDDO

  ! For pseudo absorbers
  IF ( initval == 0 ) THEN
     DO i = 1, numwin
        fitvar_rad_init (osind(i, 1:maxoth)) = 0.0D0
        fitvar_rad_init (shind(i, 1:maxoth)) = 0.0D0
     ENDDO
  ENDIF
  DO i = 1, numwin
     fitvar_rad_init (osind(i, 1:maxoth)) = 0.0D0
     fitvar_rad_init (shind(i, 1:maxoth)) = 0.0D0
     IF (ncm > 0 ) THEN 
       fitvar_rad_init (cmind(i, 1)) = 1.0D0
       fitvar_rad_init (cmind(i, 2:maxoth)) = 0.0D0
     ENDIF
     IF (np1 > 0)  fitvar_rad_init (p1ind(i, 1:maxoth)) = 0.0D0
     IF (np2 > 0)  fitvar_rad_init (p2ind(i, 1:maxoth)) = 0.0D0
     IF (np3 > 0)  fitvar_rad_init (p3ind(i, 1:maxoth)) = 0.0D0
  ENDDO

  !fitvar_rad_init(irind(1:numwin, 1))  = -1.0E-5    ! non zero
  fitvar_rad_init(irind(1:numwin, 1))  = 1.0    ! non zero

  IF (do_subfit) THEN
     nsub = numwin
  ELSE
     nsub = 1
  ENDIF

  IF (scnwrt) THEN               
      WRITE(*, '(3(A,F8.2), 2(A,i4), f5.2)') ' spres =', ps0, ' tpres =', pst, &
     ' toz = ', toz, 'snow=', the_snowice, 'land/ocean:', the_landwater_flg,the_landfrac
  ENDIF
  ! ======================================================================
  !	 Set up state vector, a priori state vector and covariance matrix
  ! ======================================================================
  ! set up fitting variables
  fitvar_rad = fitvar_rad_init    
  IF (start_layer /= 1) THEN
     lo_radbnd(ozp_fidx:ozp_fidx+start_layer-2) = &
          fitvar_rad(ozp_fidx:ozp_fidx+start_layer-2)
     up_radbnd(ozp_fidx:ozp_fidx+start_layer-2) = &
          fitvar_rad(ozp_fidx:ozp_fidx+start_layer-2)     
  ENDIF
    
  IF (end_layer /= nlay) THEN
     i = nlay - end_layer - 1
     lo_radbnd(ozp_lidx-i:ozp_lidx) = fitvar_rad(ozp_lidx-i:ozp_lidx)
     up_radbnd(ozp_lidx-i:ozp_lidx) = fitvar_rad(ozp_lidx-i:ozp_lidx)    
  ENDIF

  ! Get a priori ozone covariance matrix
  IF (use_oe .OR. (.NOT. use_oe .AND. ptr_order == 6)) THEN
     CALL GET_APRIORI_COVAR(nlay, atmosprof(1,0:nlay), &
          atmosprof(2, 0:nlay), ozprof(1:nlay), toz, ntp, sao3)
     IF (nsfc < nlay) THEN
        sao3(nsfc+1:nlay, :) = 0.0; sao3(:, nsfc+1:nlay) = 0.0
     ENDIF     
  END IF

  !xliu, 08/29/05 scaling a priori for layers below clouds to avoid smoothing even
  !for full cloudy conditions
  IF (.NOT. smooth_ozbc) THEN
     sacldscl = sacldscl * the_cfrac + (1.0 - the_cfrac )
     DO i = 1, nsfc
        IF (sacldscl(i) < 1.0) THEN 
           tmpsa = sao3(i, i)
           sao3(i, 1:nlay) = sao3(i, 1:nlay) * sacldscl(i)
           sao3(1:nlay, i) = sao3(1:nlay, i) * sacldscl(i)
           sao3(i, i) = tmpsa
        ENDIF
     ENDDO
  ENDIF
  
  ! Set up a priori state vector and covariance matrix
  IF (.NOT. use_oe) THEN  ! PTR
     CALL get_reg_matrix(sao3(start_layer:end_layer, start_layer:end_layer), &
          nump, bb(1:nf,1:nf), errstat)
     IF (errstat == pge_errstat_error)  THEN
        exval = -2; RETURN
     ENDIF

     ! has to use apriori for T since no dfs is available
     ! for ozone, it is not good to use a priori
     fitvar_rad_apriori(t_fidx:t_lidx)= fitvar_rad(t_fidx:t_lidx)      
  ELSE                    ! OE   
     ! use a priori for O3, T, 0th albedo, trace gas, Ring effect
     ! Zero for others
     fitvar_rad_apriori(ozp_fidx:ozp_lidx) = ozprof(1:nlay)
     fitvar_rad_apriori(t_fidx:t_lidx)     = fitvar_rad(t_fidx:t_lidx)
     DO i = albidx, albidx + nalb - 1
        IF (fitvar_rad_str(i)(4:4) == '0') fitvar_rad_apriori(i) = fitvar_rad(i)
     ENDDO 
     DO i = wfcidx, wfcidx + nwfc - 1
        IF (fitvar_rad_str(i)(4:4) == '0') fitvar_rad_apriori(i) = fitvar_rad(i)
     ENDDO 

     IF (fit_atanring) THEN        
        fitvar_rad_apriori(rnind(1, 1) : rnind(1, 3) + nsub - 1) = &
             fitvar_rad(rnind(1, 1) : rnind(1, 3) + nsub - 1) 
     ELSE
        fitvar_rad_apriori(rnind(1, 1) : rnind(1, 1) + nsub - 1) = &
             fitvar_rad(rnind(1, 1) : rnind(1, 1) + nsub - 1) 
     ENDIF
     IF (comfidx > 0) fitvar_rad_apriori(comvidx) = 1.0
     IF (cm1fidx > 0) fitvar_rad_apriori(cm1vidx) = 1.0

     u1idx = max_calfit_idx + (us1_idx - 1) * mxs_idx 
     u2idx = max_calfit_idx + (us2_idx - 1) * mxs_idx 
    
     ! The following covariance matrix are fixed
     IF ( npix_fitted == 0) THEN       
        DO i = 1, nf 
           j = mask_fitvar_rad(i)
           k = j - shift_offset
           ! 2500% error for parameters unless specified  (xliu, 03/21/2006)
           IF (k > 0 .AND. k < max_rs_idx) THEN  ! shift parameters specified in BOREAS.inp
              sa(i, i) = 4.0E-1 !2.5E-3
           ELSE IF (j > u1idx .AND. j <= u1idx + mxs_idx) THEN
              sa(i, i) = 0.5
           ELSE IF (j > u2idx .AND. j <= u2idx + mxs_idx) THEN
              sa(i, i) = 0.5
           ELSE IF (j >= rnind(1, 1) .AND. j <= rnind(nsub, 1)) THEN
              IF (fit_atanring) THEN ! false is default
                 sa(i, i) = 0.5  ! changed from 1.0 to 0.5 for old code toupdates
              ELSE
                 sa(i, i) = 1.0
              ENDIF
           ELSE IF (j >= rnind(1, 2) .AND. j <= rnind(nsub, 2)) THEN
              IF (fit_atanring) THEN 
                 sa(i, i) = 9.
              ELSE
                 sa(i, i) = 2.0E-3
              ENDIF
           ELSE IF (j >= rnind(1, 3) .AND. j <= rnind(nsub, 3) ) THEN
              IF (fit_atanring) THEN 
                 sa(i, i) = 9. ! changed from 4 to 9 for update v2
              ELSE
                 sa(i, i) = 1.0E-4
              ENDIF
           ELSE IF (j >= rnind(1, 4) .AND. j <= rnind(nsub, maxoth) ) THEN
              sa(i, i) = 5.0E-5   
           ELSE IF (j >= isind(1, 1) .AND. j <= isind(nsub, 1)) THEN
              sa(i, i) = 1.0E-4
           ELSE IF (j >= isind(2, 1) .AND. j <= isind(nsub, maxoth)) THEN
              sa(i, i) = 1.0E-4
           ELSE IF (j >= irind(1, 1) .AND. j <= irind(nsub, 1)) THEN
              !sa(i, i) = 1.0E-4  !10.0
              sa(i, i) = 1.0E+5 
           ELSE IF (j >= irind(2, 1) .AND. j <= irind(nsub, maxoth)) THEN
              !sa(i, i) = 1.0E-3 !2.0
              sa(i, i) = 1.0E+3
           ELSE IF (j >= dcind(1, 1) .AND. j <= dcind(nsub, 1) ) THEN
              sa(i, i) = 0.05**2.0
           ELSE IF (j >= dcind(2, 1) .AND. j <= dcind(nsub, maxoth)) THEN
              sa(i, i) = 0.01**2
           ELSE IF (i == comfidx .OR. i == cm1fidx .OR. i == cm2fidx .OR. i == cm3fidx) THEN
              sa(i, i) = 1.0
          ELSE IF (j >= cmind(1,1) .AND. J <= cmind(nsub, maxoth)) THEN  !JBAK
              sa(i,i)  = 1.0
           ELSE IF (i < ozf_fidx .OR. i > ozf_lidx) THEN
              sa(i, i) = (fitvar_rad(j))**2.0 * 25.0 
           ENDIF
        ENDDO
     
        IF (nt_fit > 0) THEN   ! 5 K std. deviation
           DO i = tf_fidx, tf_lidx
              sa(i, i) =  25.0
           ENDDO
        ENDIF

        ! zero-order error +/-0.05 nm
        ! error decreases by a factor of 10 when the order increases by 1
        IF (nos > 0) THEN
           DO i = 1, nos
              DO j = 1, nsub
                 k = osfind(j, i)
                 IF (k > 0) sa(k, k) = 5.0E-4 * (10.0 ** (-(i - 1) * 2.0)) 
              ENDDO
           ENDDO
        ENDIF
        
        ! zero-order error +/-0.01 nm
        ! error decreases by a factor of 10 when the order increases by 1
        IF (nsh > 0) THEN
           DO i = 1, nsh
              DO j = 1, nsub
                 k = shfind(j, i)
                 IF (k > 0) sa(k, k) = 5.0E-4  * (10.0 ** (-(i - 1) * 2.0))
              ENDDO
           ENDDO
        ENDIF

        IF (np1 > 0) THEN
           !print *, ' ask to xiong what number should I put here ?'
           DO i = 1, np1
              DO j = 1, nsub
                 k = p1find(j, i)
                 iF (k > 0) sa(k, k) = 1*1.0E-1 * (10.0 ** (-(i - 1) * 2.0))
                 !iF (k > 0) sa(k, k) = 0.1 * (10.0 ** (-(2.*i - 2) * 2.0))
              ENDDO
           ENDDO
        ENDIF
        IF (np2 > 0) THEN
           DO i = 1, np2
              DO j = 1, nsub
                 k = p2find(j, i)
                 IF (k > 0) sa(k, k) = 1*1.0E-1  * (10.0 ** (-(i - 1) * 2.0))
                 !IF (k > 0) sa(k, k) = 0.1  * (10.0 ** (-(i*2. - 2) * 2.0))
              ENDDO
           ENDDO
        ENDIF
        IF (np3 > 0) THEN
           DO i = 1, np3
              DO j = 1, nsub
                 k = p3find(j, i)
                 IF (k > 0) sa(k, k) = 1*1.0E-1  * (10.0 ** (-(i - 1) * 2.0))
              ENDDO
           ENDDO
        ENDIF
     ENDIF

     ! need to update o3 a priori covariance matrix for each retrieval
     sa(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx) = &
          sao3(start_layer:end_layer,start_layer:end_layer) 

     IF (nfalb > 0 .OR. nfwfc > 0 .OR. ecfrfind > 0) THEN
        albfc_aperr = 0.05**2.0
        IF (tropaod(1) >= 0.25 .AND. taodfind == 0 .AND. twaefind == 0) THEN
           albfc_aperr1 = albfc_aperr * (1.0 + 8.0 * tropaod(1) - 2.0) 
        ELSE IF (has_glint) THEN  ! Assume a priori error of 0.2 instead of 0.05 for 100% sun glint
           albfc_aperr1 = albfc_aperr * ( 1.0 + 15.0 * glintprob) 
        ELSE
           albfc_aperr1 = albfc_aperr 
        ENDIF
        albfc_aperr2    = albfc_aperr * (1.0 + 15.0 * salbedo)
        albfc_aperr     = MAX(albfc_aperr1, albfc_aperr2)
     ENDIF

     IF (nfalb > 0) THEN
       off = 0
       DO i = albfidx, albfidx + nfalb - 1
         j = mask_fitvar_rad(i)
         k = j - albidx + 1
         ! The a priori std. for non-zero albedo terms are based on retrievals
         ! (1.6E-3, 1.0E-5, 1.0E-7, 1.0E-9)
         ! (3.0E-2, 4.0E-4, 1.6E-5, 6.4E-7)
         fidx = albfpix(mask_fitvar_rad(i)-albidx+1)
         lidx = alblpix(mask_fitvar_rad(i)-albidx+1)
         READ (fitvar_rad_str(mask_fitvar_rad(i))(4:5), '(I2.2)') nord
         IF (.NOT. is_albspcvar(k) ) THEN
           !sa(i, i) = albfc_aperr * (5.0 ** ( - nord * 2.0)) !IN OMI
           sa(i, i) = albfc_aperr * ( ((fitwavs(lidx)-fitwavs(fidx))*0.5) ** (-nord * 2.0) )  ! IN GOME
           off = off + 1
           IF (salbedo > 0.6 .AND. nord >= 1) sa(i, i) = 0.0
         ELSE ! BRDF or EOF or Spectrum
           IF (do_brdf) THEN
              sa(i,i) = 0.10
           ELSE IF (use_albeofs) THEN
             IF (which_albspc == 1 .OR. snowflg == 1 .OR. oceanflg == 1) THEN 
               IF (nactalbspc > 1 ) THEN  ! peter alb spectrum
                 sa(i,i) = 2.0**2.0
                 ! For water/snow, only 1 EOF (same as water/snow spectra),
                 ! A priori error for 2 to nalbspc EOF is zero, fixed to be zero
               ELSE IF (nactalbspc == 1) THEN
                 !sa(i, i) = (0.2d0 ** 2)
                 sa(i, i) = (0.02d0 ** 2) !* ( ((fitwavs(lidx)-fitwavs(fidx))*0.5)** (-(nord-1) * 2.0) )
                 IF (nord > 1) sa(i,i) = 0.0
               ENDIF
             ELSE IF (which_albspc == 2) THEN ! chris alb spectrum
                IF (nactalbspc < nalbspc ) THEN 
                   sa(i, i) = (0.02d0 *2)! 
                  ! IF (nord > 1) sa(i,i) = 0 
                ELSE
                   fidx = 1 + albfidx + off -1
                   lidx = nactalbspc + albfidx + off -1
                   sa(i,fidx:lidx) = Surface%Option4%FactorUncertainty(nord,1:nactalbspc)**2 ! more constraint with 1.5 scaling
                   !sa(i,i) = Surface%Option4%FactorUncertainty(nord,nord) ! more constraint with 1.5 scaling
                ENDIF
             ENDIF
           ELSE ! Scale albedo spectrum
             sa(i, i) = (0.2d0 ** 2) * ( ((fitwavs(lidx)-fitwavs(fidx))*0.5) ** (-nord * 2.0) )
           ENDIF 
         ENDIF
       ENDDO
     ENDIF
     IF (nfwfc > 0) THEN
        DO i = wfcfidx, wfcfidx + nfwfc - 1
           READ (fitvar_rad_str(mask_fitvar_rad(i))(4:4), '(I1.1)') nord
           sa(i, i) = albfc_aperr * (5.0 ** ( - nord * 2.0))
           IF (the_cfrac <= 1.0D-2 .OR. the_cfrac >= 0.99) sa(i, i) = 0.0 ! added for update v2
        ENDDO
     ENDIF
     
     ! Use 50% (NO2 and HCHO) or 100% (SO2 and BrO) error for other minor trace gases
     ! The apiori of these trace gases are determined by climatology (NO2, HCHO)
   
     DO i = 1, ngas
        j = fgasidxs(i)
        IF (j > 0) THEN
           sa(j, j) = (fitvar_rad_apriori(mask_fitvar_rad(j)))**2.0
           !IF (gasidxs(i) == hcho_idx .OR. gasidxs(i) == no2_t1_idx .OR. &
           !     gasidxs(i) == no2_t2_idx ) sa(j, j) = sa(j, j)* 0.25 !50% error
           IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. use_large_so2_aperr)  &
                sa(j, j) = (1.0E17 * refspec_norm(gasidxs(i))) ** 2.0
           IF (gasidxs(i) == o2o2_idx) sa(j, j) = sa(j, j) * 0.09 ! 30% error
           IF (gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx) sa(j, j) = sa(j, j) * 0.09 ! 30% error
           IF (gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx) sa(j, j) = sa(j, j) !100% error of AFGLUS
        ENDIF
     ENDDO

     ! A priori covariance matrix for aerosols and clouds
     IF (taodfind > 0) THEN  ! 100%
        sa(taodfind, taodfind) = fitvar_rad_apriori(taodind)**2   !* 0.25
     ENDIF
     IF (saodfind > 0) THEN  ! 50%
        sa(saodfind, saodfind) = fitvar_rad_apriori(saodind)**2 * 0.25
     ENDIF
     IF (twaefind > 0) THEN  ! aerosol single scattering albedo change by 0.05
        sa(twaefind, twaefind) = 2.5E-3
     ENDIF
     IF (ecfrfind > 0) THEN  
        sa(ecfrfind, ecfrfind) = albfc_aperr
     ENDIF
     IF (ecodfind > 0) THEN  ! 10%
        sa(ecodfind, ecodfind) = (fitvar_rad_apriori(ecodind) * 0.1)**2
     ENDIF
     IF (ectpfind > 0) THEN  ! 100 mb
        sa(ectpfind, ectpfind) = 100.**2
     ENDIF    

     IF (sprsfind > 0) THEN  ! 100 mb
        sa(sprsfind, sprsfind) = 30.**2
     ENDIF

     IF (so2zfind > 0) THEN  ! 1.0 km
        sa(so2zfind, so2zfind) = 1.0 ** 2.
     ENDIF
        
  ENDIF
    
  ! Initialize exit value to be zero (missing)
  exval = 0  
    
  ! Create a condensed array of fitting variables that are varied. 
  ! This considerably reduces the execution time of the fitting routine.
  IF (radcalwrt .AND. do_simu) fitvar_rad = fitvar_rad_apriori
  fitvar(1:nf)    = fitvar_rad(mask_fitvar_rad(1:nf))
  fitvarap(1:nf)  = fitvar_rad_apriori(mask_fitvar_rad(1:nf))
  lowbond(1:nf)   = lo_radbnd(mask_fitvar_rad(1:nf))
  upbond(1:nf)    = up_radbnd(mask_fitvar_rad(1:nf))

  !print*, fitvar_rad
  IF (do_debug_o3p) WRITE(www_lun, '(A)') 'start ozprof_invesre'
  CALL ozprof_inverse (nf, varname(1:nf), fitvar(1:nf), fitvarap(1:nf), &
       lowbond(1:nf), upbond(1:nf), npoints, nump, sa(1:nf,1:nf), bb(1:nf,1:nf), &
       !chisq, fitspec_rad(1:npoints), fitres_rad(1:npoints), exval)
       chisq, fitspec_rad(1:npoints), fitres_rad(1:npoints), fitspec_q(1:npoints), fitres_q(1:npoints), exval)

  fitvar_rad(mask_fitvar_rad(1:nf)) = fitvar(1:nf)  ! for safe 
  fitvar_rad_apriori(mask_fitvar_rad(1:nf)) = fitvarap(1:nf)  ! Some a priori values can be changed

  DO i = 1, nf
     fitvar_rad_aperror(i) = SQRT(sa(i, i))
  ENDDO
  
  ! reset the wavelength shifts to zero if retrievals are not successful
  ! because the failure of retrievals are due to too large wavelength shifts most of the time
  IF (exval <= 0) THEN
     DO i = 1, numwin
        fitvar_rad_init(osind(i, 1:maxoth)) = 0.0;   fitvar_rad_init(shind(i, 1:maxoth)) = 0.0
     ENDDO
  ENDIF
    
  IF (exval < 0) THEN     ! Terminate the whole retrieval
     fitcol  = missing_value_dp
     dfitcol = missing_value_dp    
     fitvar_rad_saved = fitvar_rad_init
     RETURN
  ELSE
     IF (exval == 0 .OR. the_cfrac > 0.4) THEN
        fitvar_rad_saved = fitvar_rad_init 
     ELSE
        fitvar_rad_saved = fitvar_rad 
     ENDIF
  ENDIF
  
  ! calculate rms difference between measurements and calculations
  ! If using measurement error, rms is ideally to be 1, if it > 1, suggesting
  ! the estiamted measurement error is too small, and vice versa too large
  rms = SQRT(chisq / REAL(npoints, KIND=dp))
  
  ! need to multiply rms to get the actual retrieval random noise error
  ! no matter whether use measurement error or not (sig=1.0)
  ! because measurement error is not reliable, a rms of > 1 indicates
  ! underestimated measurement error and vice versa
  ! but for smoothing error, do not do it
  ! covar(1:nf,1:nf)=rms**2.0*covar(1:nf, 1:nf) !*npoints/(npoints-nf)
  
  ozprof(1:nlay) = fitvar_rad(ozp_fidx:ozp_lidx)
  DO i = 1, nf
     stderr(i) = SQRT(covar(i, i)); stderr1(i) = SQRT(ncovar(i, i))
  END DO
  fitvar_rad_std(mask_fitvar_rad(1:nf))  = stderr(1:nf)
  fitvar_rad_nstd(mask_fitvar_rad(1:nf)) = stderr1(1:nf)
  
  ! compute error in total ozone, stratospheric ozone, tropospheric ozone
  ! -------------------------------------------------------------
  ! FITCOL:       Used for total o3, strat o3, and trop o3
  ! DFITCOL:      Uncertainty for corresponding columns
  ! NLAY_FIT:     Number of ozone layers that are varied
  ! STATRT_LAYER: first layer that are varied among all layers
  ! END_LAYER:    last layer that are varied among all layers
  ! -------------------------------------------------------------
  fitcol = 0.0  ;  dfitcol = 0.0
  fitcol(1) = SUM (ozprof(1:nsfc)) 
  j = start_layer 
  DO i = ozf_fidx, ozf_lidx
     ozprof_std(j) = stderr(i); ozprof_nstd(j) = stderr1(i); j = j + 1
  END DO

  k1 = ozf_fidx - start_layer
  DO is = 1, 2
     IF (is == 1 ) THEN
        ! xliu, 4/08/2016, bug, should be below
        !ozprof_std_sav(ozf_fidx:ozf_lidx) = ozprof_std(ozf_fidx:ozf_lidx)
        ozprof_std_sav(start_layer:nsfc) = ozprof_std(start_layer:nsfc)
        covar_sav(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx) = covar(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx)
     ELSE
        IF (.NOT. use_oe) EXIT
        ! xliu, 4/08/2016, bug, should be below
        !ozprof_std_sav(ozf_fidx:ozf_lidx) = ozprof_nstd(ozf_fidx:ozf_lidx)
        ozprof_std_sav(start_layer:nsfc) = ozprof_nstd(start_layer:nsfc)
        covar_sav(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx) = ncovar(ozf_fidx:ozf_lidx, ozf_fidx:ozf_lidx)
     ENDIF
            
     ! remove correlated error to get error in total ozone
     DO i = start_layer, nsfc
        IF (ozprof_std_sav(i) > 0.0) THEN
           dfitcol(1, is) = dfitcol(1, is) + ozprof_std_sav(i)  ** 2.0
        
           DO j = start_layer, i - 1
              dfitcol(1, is) = dfitcol(1, is) + 2.0 * ozprof_std_sav(i) * ozprof_std_sav(j) * &
                   covar_sav(k1+i, k1+j) / SQRT(covar_sav(k1+i, k1+i) * covar_sav(k1+j, k1+j))
           END DO
        ENDIF
     ENDDO
  
     IF (ntp > 0) THEN     
        ! stratospheric ozone
        DO i = start_layer, ntp
           IF (ozprof_std_sav(i) > 0.0) THEN
              dfitcol(2, is) = dfitcol(2, is) + ozprof_std_sav(i)**2.0
              DO j = start_layer, i - 1
                 dfitcol(2, is) = dfitcol(2, is) + 2.0 * ozprof_std_sav(i) * ozprof_std_sav(j) * &
                      covar_sav(k1+i, k1+j) / SQRT(covar_sav(k1+i, k1+i) * covar_sav(k1+j, k1+j))
              ENDDO
           ENDIF
        ENDDO

        ! tropospheric ozone
        DO i = ntp + 1, nsfc
           IF (ozprof_std_sav(i) > 0.0) THEN
              dfitcol(3, is) = dfitcol(3, is) + ozprof_std_sav(i)**2.0
              DO j = ntp + 1, i - 1
                 dfitcol(3, is) = dfitcol(3, is) + 2.0 * ozprof_std_sav(i) * ozprof_std_sav(j) *  &
                      covar_sav(k1+i, k1+j) / SQRT(covar_sav(k1+i, k1+i) * covar_sav(k1+j, k1+j))

                 !print *, ozprof_std_sav(i), ozprof_std_sav(j), sqrt(covar_sav(k1+i, k1+i)), sqrt(covar_sav(k1+j, k1+j))
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO

  dfitcol(1, :)    = SQRT (dfitcol(1, :))
  IF (ntp > 0) THEN
     fitcol(2)     = SUM (ozprof(1:ntp)) 
     fitcol(3)     = SUM (ozprof(ntp+1:nsfc))
     dfitcol(2, :) = SQRT (dfitcol(2, :))
     dfitcol(3, :) = SQRT (dfitcol(3, :))
  ENDIF

  IF (nwfc > 0) THEN
     DO i = wfcidx + nwfc - 1, wfcidx, -1
        IF (fitvar_rad_str(i)(4:4) == '0') EXIT
     ENDDO
     the_cfrac = fitvar_rad(i)   ! Use UV2/last channel cloud fraction
  ELSE IF ( ecfrfind > 0) THEN
     the_cfrac = fitvar_rad(ecfrind )
  ENDIF
  
  !IF (the_cfrac <= 1.0D-3 .AND. ecfrind > 0) THEN
  !   the_cod = 0.d0; the_ctp = 0.d0
  !ENDIF
  
  ozprof(1:nlay) = fitvar_rad(ozp_fidx:ozp_lidx)
  ozprof_init(1:nlay) = fitvar_rad_init(ozp_fidx:ozp_lidx)
  ozprof_ap(1:nlay) = fitvar_rad_apriori(ozp_fidx:ozp_lidx) 
  ozprof_apstd = 0.0
  j = start_layer
  IF (use_logstate) THEN
     DO i = ozf_fidx, ozf_lidx
        ozprof_apstd(j) = SQRT(sa(i, i))*ozprof_ap(j); j = j + 1
     ENDDO
  ELSE
     DO i = ozf_fidx, ozf_lidx
        ozprof_apstd(j) = SQRT(sa(i, i)); j = j + 1
     ENDDO
  ENDIF  
  
  IF (nsfc < nlay) THEN
     ozprof(nsfc+1:nlay) = -99.99;     ozprof_ap(nsfc+1:nlay) = -99.99
  ENDIF

  eff_alb_init = fitvar_rad_init(albidx:albidx+maxalb-1)
  eff_alb      = fitvar_rad(albidx:albidx+maxalb-1) 
  eff_wfc_init = fitvar_rad_init(wfcidx:wfcidx+maxwfc-1)
  eff_wfc      = fitvar_rad(wfcidx:wfcidx+maxwfc-1)
 !Junsung
 !print*, eff_alb_init 
 !print*, eff_alb
 !stop
  ! Derive Final Albedo spectraum
  !IF (wrtalbspc) THEN
    DO i = 1, nalb
      j = albidx - 1 + i
      READ(fitvar_rad_str(j)(4:4), '(I1)') nord
      fidx = albfpix(i); lidx = alblpix(i)
      wavavg = (fitwavs(fidx) + fitwavs(lidx)) / 2.0d0

      IF (.NOT. is_albspcvar(i)) THEN
        IF (nord == 0  ) THEN
          sfcalbs(fidx:lidx, 2) = fitvar_rad(j)
        ELSE
          sfcalbs(fidx:lidx, 2) = sfcalbs(fidx:lidx, 2) + fitvar_rad(j) * (fitwavs(fidx:lidx) - wavavg)**nord
        ENDIF
      ELSE
        IF (use_albeofs) THEN
          IF (nactalbspc == 1) THEN ! Snow/ice
            sfcalbs(fidx:lidx, 2) = albspcs(fidx:lidx, 0) * fitvar_rad(j)
          ELSE
            IF (nord == 1)  sfcalbs(fidx:lidx, 2) = albspcs(fidx:lidx, 0)
              sfcalbs(fidx:lidx, 2) = sfcalbs(fidx:lidx, 2) + albspcs(fidx:lidx, nord) * fitvar_rad(j)
          ENDIF
        ELSE
          IF (nord == 0) THEN
            sfcalbs(fidx:lidx, 2) = albspcs(fidx:lidx, 0) * fitvar_rad(j)
          ELSE
            sfcalbs(fidx:lidx, 2) =  sfcalbs(fidx:lidx, 2) + (fitwavs(fidx:lidx) - wavavg)**nord &
              * albspcs(fidx:lidx, 0) * fitvar_rad(j)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  !ENDIF

 IF (do_debug_out) THEN
 unit_deg = 2 ! 1=chris at 1, 2=peter at 2
 IF (which_albspc == 2) unit_deg = 1
 WRITE(unit_deg,*) '------------------------------------------------------'   
 WRITE(unit_deg, '(3i8,2f8.2,f8.5, 3i8)') exval,currpix, currline,the_lon, the_lat,the_landfrac,&
  the_landwater_flg, the_snowice, npoints-widx_vis+1
 allocate(tmpspc(npoints))
 tmpspc = exp(fitspec_rad(1:npoints) - fitres_rad(1:npoints))
 WRITE(unit_deg,'(800f8.2)') fitwavs(widx_vis:npoints)
 WRITE(unit_deg,'(800e15.7)') sfcalbs(widx_vis:npoints,1)
 WRITE(unit_deg,'(800e15.7)') sfcalbs(widx_vis:npoints,2)
 WRITE(unit_deg,'(800e15.7)') exp(fitspec_rad(widx_vis:npoints))
 WRITE(unit_deg,'(800e15.7)') tmpspc(widx_vis:npoints)
 deallocate(tmpspc)
 ENDIF

 ! trace gases
 DO k = 1, ngas
    i = fgasidxs(k)
    IF (i > 0) THEN
       j = mask_fitvar_rad(i)
       tracegas(k, 1) = fitvar_rad_init(j)
       tracegas(k, 2) = fitvar_rad_apriori(j)
       tracegas(k, 3) = SQRT(sa(i, i))
       tracegas(k, 4) = fitvar_rad(j)
       tracegas(k, 5) = fitvar_rad_std(j) 
       tracegas(k, 6) = fitvar_rad_nstd(j) 
       tracegas(k, 1:6) = tracegas(k, 1:6) / refspec_norm(gasidxs(k))  
       fitvar_rad(j) = tracegas(k, 4)
       fitvar_rad_std(j) = tracegas(k, 5); fitvar_rad_nstd(j) = tracegas(k, 6)
       ! xliu: 07/01/2010, 08/09/2010
       ! Change weighting function and contribution function for trace gas variables 
       ! wrt to the reported unit_deg instead of normalized quantities
       IF (ozwrtwf) THEN
          weight_function(1:npoints, i) = weight_function(1:npoints, i) * refspec_norm(gasidxs(k)) 
          trace_profwf(k, 1:npoints, 1:nlay) = trace_profwf(k, 1:npoints, 1:nlay) * refspec_norm(gasidxs(k))
       ENDIF
       IF (ozwrtcontri) THEN 
          contri(i, 1:npoints) = contri(i, 1:npoints) / refspec_norm(gasidxs(k)) 
          trace_contri(k, 1:npoints) = contri(i, 1:npoints)
       ENDIF

       ! Add trace gas profile (a priori profile/shape)
       DO j = 1, nlay
          fidx = nup2p(j - 1) + 1; lidx = nup2p(j)
          trace_prof(k, j) = SUM(mgasprof(k, fidx:lidx))
       ENDDO
    ENDIF
 ENDDO 
  
 IF ( taodfind > 0) THEN
    aodscl = fitvar_rad(taodind) / tropaod(actawin)
    tropaod(1:actawin) = tropaod(1:actawin) * aodscl
    IF ( twaefind == 0 ) THEN  ! Single scattering albedo does not change
       tropsca(1:actawin) = tropsca(1:actawin) * aodscl
    ENDIF
 ELSE
    aodscl = 1.0
 ENDIF

IF ( twaefind > 0 ) THEN
    waerscl = fitvar_rad(twaeind) / tropwaer(actawin)   ! Scale single scattering albedo
    tropwaer(1:actawin) = tropwaer(1:actawin) * waerscl
    tropsca(1:actawin)  = tropsca(1:actawin) * waerscl * aodscl
 ENDIF

 IF ( saodfind > 0 ) THEN
    aodscl = fitvar_rad(saodind) / strataod(actawin)
    strataod(1:actawin) = strataod(1:actawin) * aodscl    
    stratsca(1:actawin) = stratsca(1:actawin) * aodscl
 ENDIF

 IF ( ecfrfind > 0) THEN
    the_cfrac = fitvar_rad(ecfrind )
 ENDIF

 IF (sprsfind > 0) THEN
    atmosprof(1, nsfc) = fitvar_rad(sprsind)
 ENDIF
  
 RETURN
END SUBROUTINE specfit_ozprof

END MODULE m_specfit_ozprof
