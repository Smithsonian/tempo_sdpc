!
module m_set_cldalb

  public set_cldalb
  private !setcld

contains

  ! ===================================================================
  !	        Set up clouds and surface albedo  (pain in the neck)
  ! ===================================================================
  ! Search for cloud information and surface albedo
  ! If has clouds, use GOME derived albedo by Kolemeijer, then use fixed
  ! albedoes, need to adjust fitvar_rad and fitvar
  ! If there no clouds, then derive the surface albedo from 370.2 nm
  ! Need to override the specified albedo values
  ! Need to check for no albedo or multiple albedo specified at some
  ! wavelengths

  SUBROUTINE SET_CLDALB(npoints, fitwavs, ctau, ctp, cfrac, salbedo, &
       pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY:instrument_idx, omi_idx,tempo_idx
    USE OMSAO_variables_module, ONLY: fitvar_rad_init, fitvar_rad_str, &
         lo_radbnd, up_radbnd, nf=>n_fitvar_rad, &
         mask_fitvar_rad, rmask_fitvar_rad, &
         the_month, the_day, the_jday, edgelons, edgelats, & 
         the_sza_atm,the_vza_atm, the_aza_atm, numwin,nradpix, nviswin
    USE ozprof_data_module, ONLY: nrefl,albidx, albfidx, nalb, nfalb, albmin, &
         albmax, albfpix, alblpix, do_lambcld, lambcld_refl, &
         which_alb, the_fixalb, do_simu, radcalwrt, lambcld_initalb, &
         ecfrfind, wfcmax, wfcmin, nwfc, nfwfc, &
         wfcfpix, wfclpix, wfcfidx, wfcidx, the_snowice, pos_alb, &
         use_albspc, albspcs, malbspc, nactalbspc, is_albspcvar,sfcalbs
    !ps0, has_glint, glintprob, cloud
    !USE OMSAO_gome_data_module, ONLY: gome_pixdate, gome_orbc, gome_curpix, &
    !     gome_curscan
    USE OMSAO_errstat_module
    use m_get_initial_albedo, only: adj_albcfrac, get_gome_alb, &
       get_initial_albedo, get_omi_alb, get_omler_alb, get_toms_alb, &
       get_sciagm2_alb, get_surface_spectrum  
    IMPLICIT NONE

    ! ========================
    ! Input/Output Variables
    ! ========================
    INTEGER, INTENT(IN)                 :: npoints
    REAL(KIND=dp), DIMENSION(npoints)   :: fitwavs
    REAL (KIND=dp), INTENT(INOUT)       :: ctau, ctp, cfrac
    REAL (KIND=dp), INTENT(OUT)         :: salbedo

    INTEGER, INTENT(OUT)                :: pge_error_status

    ! ===============
    ! Local variables
    ! ===============
    INTEGER :: i, j, k, ninalb,nsub, region, fidx, lidx,nw,ntmp, snowflg, nalbw, &
               which_sciagm2
    INTEGER, DIMENSION(npoints)  :: hasalb, haswfc
    REAL (KIND=dp)               :: cfrac_old, albedo, wavg!, actp
    INTEGER, DIMENSION(2)        :: wavin
    REAL (KIND=dp), DIMENSION(2) :: wavfrac
    INTEGER, PARAMETER           :: nalbw0=29
    REAL (KIND=dp), DIMENSION (nalbw0):: albarr, albwave
    LOGICAL                      :: noalb, do_adjcfrac

    pge_error_status = pge_errstat_ok
    nalbw = 7
    albwave(1:nalbw) = (/335.0, 380.0, 440.0, 495., 555., 610., 670.0/)
    IF (which_alb == 1) THEN
      !ninalb = 1; region = 1
      ninalb=7 ; region = 3
      CALL GET_GOME_ALB(the_month, edgelons, edgelats, region, albarr, ninalb)
    ELSE IF (which_alb == 2) THEN
      CALL GET_TOMS_ALB(the_month, edgelons, edgelats, albedo)
      albarr(:) = albedo
    ELSE IF (which_alb == 3) THEN
      CALL GET_OMI_ALB(the_month, the_day, edgelons, edgelats, albedo)
      albarr(:) = albedo
    ELSE IF (which_alb == 4) THEN
     nalbw = 7
     albwave(1:nalbw) = (/328.0, 345.0, 354.0, 380.0, 442.0, 477.0, 499.0/)
      CALL GET_OMLER_ALB(the_month, the_day, edgelons, edgelats, albedo)
      albarr(:) = albedo
    ELSE IF (which_alb >=5 .AND. which_alb <=7) THEN
      which_sciagm2 = which_alb - 4
      CALL GET_SCIAGM2_ALB(which_sciagm2, the_month, the_day, edgelons, edgelats,albarr, albwave, nalbw)
    ELSE
      WRITE(www_lun, *) 'Albedo database: not implemented!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF
    
    IF (use_albspc) THEN
     IF (instrument_idx == omi_idx) THEN 
      IF (the_snowice>= 10 .and. the_snowice /= 104) THEN
          snowflg =1
      ELSE
          snowflg =0
      ENDIF
     ELSE IF (instrument_idx == tempo_idx) THEN 
       snowflg = 0
     ENDIF
      fidx = SUM(nradpix(1 : numwin - nviswin)) + 1
      lidx = npoints
      ntmp = lidx - fidx + 1

      CALL get_surface_spectrum (the_jday, the_sza_atm, the_aza_atm, the_vza_atm, edgelons,edgelats, &
          snowflg,  ntmp, fitwavs(fidx:lidx), albspcs(fidx:lidx, 0:malbspc-1),nactalbspc, pge_error_status)
      IF (snowflg == 1) THEN
        albspcs(fidx:lidx, 0) = albspcs(fidx:lidx, 0) * the_snowice / 100.0
      ENDIF
      IF (pge_error_status == pge_errstat_error) THEN
        WRITE(*, *) ': Error in getting surface albedo spectrum';
        RETURN
      ENDIF
    ENDIF
    ! get initial albedo
    IF (nrefl > 0 ) THEN
      IF (ctp > 0.0) THEN
        noalb = .TRUE.
      ELSE
        noalb = .FALSE.
      ENDIF  
      ! get effective surface albedo if noalb = .false.
      ! otherwise, compute I/F at albedo wavelength for further derivation of
      ! surface albedo and/or initial cloud fraction
      CALL GET_INITIAL_ALBEDO(noalb, albedo, pge_error_status)
      IF (pge_error_status == pge_errstat_error)  RETURN
    ELSE
      noalb = .TRUE.
    ENDIF


    IF (noalb) THEN 
      IF (which_alb == 2 .OR. which_alb == 3) THEN
        albedo = albarr(1)
      ELSE
        DO i = 1, nalbw
          IF (albwave(i) >= pos_alb) EXIT
        ENDDO
        IF (i == 1) THEN 
          nw = 1 ; wavin(1) = 1; wavfrac(1) = 1.0
        ELSE IF (i == nalbw  .AND. pos_alb >=albwave(nalbw)) THEN
          nw = 1 ; wavin(1) = nalbw; wavfrac(1) = 1.0
        ELSE
          nw = 2; wavin(1) = i - 1; wavin(2) = i
          wavfrac(2) = (pos_alb - albwave(i-1)) / (albwave(i)-albwave(i-1))
          wavfrac(1) = 1.0 - wavfrac(2)
        ENDIF
         albedo = 0.0
        DO i = 1, nw
           albedo = albedo + albarr(wavin(i)) * wavfrac(i)
        ENDDO
      ENDIF
    ELSE
     IF (albedo > 1.0) THEN
        WRITE(*, *) 'Thick clouds exist. No retrieval done or add cloud info!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF
    ENDIF

    IF (do_simu .AND. .NOT. radcalwrt) albedo = the_fixalb

    ! adjust the cloud fraction or surface albedo based on
    ! LIDORT-calculated radiance at 370.2 nm
    !CALL ADJ_ALBCFRAC(albedo, cfrac, pge_error_status)
    !IF (pge_error_status == pge_errstat_error)  RETURN

    ! Initialize Lambertian cloud albedo to be 80%
    IF ( .NOT. (do_simu .AND. .NOT. radcalwrt) ) THEN
      IF (do_lambcld ) lambcld_refl = lambcld_initalb
    ENDIF


    ! For simulation, if cfrac == 1.0, then need to adjust lambcld_refl.
    ! However, this already partly corrects for calibration offset by
    ! deriving a lambcld_refl that matches the measured radiance at
    ! cloud wavelength
    ! IF (.NOT. do_simu .OR. (cfrac == 1.0D0 .AND. radcalwrt) )
    !IF (the_snowice == 101) THEN
    !  albedo = 0.90
    !ELSE IF (the_snowice == 103) THEN
    !  albedo = 0.80
    !ELSE IF (the_snowice > 1 .AND. the_snowice <= 100) THEN
    !  albedo = MAX(albedo, 0.8d0 * DBLE(the_snowice) / 100.d0)
    !ENDIF
    ! 104 : OCEAN
    ! xliu, 12/07/2014, changed based on ASTER snow spectrum
    IF (instrument_idx == omi_idx) THEN 
     IF (the_snowice == 101) THEN
      albedo = 0.98
      albarr = albedo
     ELSE IF (the_snowice == 103) THEN
      albedo = 0.8
      albarr = albedo
     ELSE IF (the_snowice > 1 .AND. the_snowice <= 100) THEN
      albedo = MAX(albedo, 0.98 * the_snowice / 100.0)
      albarr = albedo
     ENDIF
    ELSE IF (instrument_idx == tempo_idx) THEN 
      
    ELSE
      STOP
    ENDIF
    do_adjcfrac = .TRUE.
    IF (nrefl < 1) do_adjcfrac = .FALSE.
    IF (do_adjcfrac) THEN 
      cfrac_old = cfrac  ! save cloud fraction from other products
      CALL ADJ_ALBCFRAC(albedo, cfrac, ctau, pge_error_status)
    ENDIF 
    salbedo = albedo   ! Surface albedo
    ! xliu: 08/16/2008, when surface albedo increases, it is more difficult
    ! to differentiate clouds/surfaces
    ! Assume a cloud fraction of 0 and increases a priori error for surface
    ! albedo and cloud fraction
    IF (albedo > 0.6 .AND. cfrac >= 0.6 .AND. .NOT. do_simu) cfrac = 0.5

    IF (ecfrfind == 0 .AND. nwfc == 0) THEN
      ! Final cloud fraction is computed here
      IF (ctp > 0.0) THEN
        IF (do_lambcld .AND. cfrac >= 1.0D0) albedo = lambcld_refl
      ENDIF
    ELSE
      ! The derived cloud fraction is derived as initial value
      ! Slightly change the value for clear-sky/cloud-sky, so weighting
      ! function for clouds are calculated
      IF ( cfrac == 0.0) cfrac = 0.01
      IF ( cfrac == 1.0) cfrac = 0.99
      !IF ( has_glint .AND. cfrac < 0.20 * glintprob ) cfrac = 0.0D0
    ENDIF

    DO i =  1, nalb
      j = albidx + i - 1
      IF  (fitvar_rad_str(j)(4:4) =='0' .AND. .NOT. is_albspcvar(i)) THEN
        IF (which_alb == 1 .OR. which_alb ==4) THEN
          wavg = (albmin(i) + albmax(i))/2.0
          DO k = 1, nalbw
            IF (albwave(k) >= wavg) EXIT
          ENDDO  
          IF (k == 1) THEN
            nw = 1; wavin(1) = 1; wavfrac(1) = 1.0
          ELSE IF ( k >= nalbw .AND. wavg >= albwave(nalbw) ) THEN
            nw = 1; wavin(1) = nalbw; wavfrac(1) = 1.0
          ELSE
            nw = 2; wavin(1) = k - 1; wavin(2) = k
            wavfrac(2) = (wavg - albwave(k-1)) / (albwave(k)-albwave(k-1))
            wavfrac(1) = 1.0 - wavfrac(2)
          ENDIF
          albedo = 0.0
          DO k = 1, nw
            albedo = albedo + albarr(wavin(k)) * wavfrac(k)
          ENDDO
        ENDIF
        
        fitvar_rad_init(j) = albedo
        IF (up_radbnd(j) == lo_radbnd(j)) THEN
          up_radbnd(j) = albedo
          lo_radbnd(j) = albedo
        ENDIF
      ENDIF
    ENDDO 

    ! Go thorugh albedo terms again to check for unused
   
    i= 1; k = albfidx - 1
    DO i = 1, nalb
      j = albidx + i - 1
      IF (lo_radbnd(j) < up_radbnd(j))  k = k + 1
      nsub = COUNT(MASK=(fitwavs(1:npoints) >= albmin(i) .AND.fitwavs(1:npoints) <= albmax(i)))
      IF (lo_radbnd(j) < up_radbnd(j) .AND. nsub <= 0)  THEN
        mask_fitvar_rad(k:nf-1) = mask_fitvar_rad(k+1:nf)
        rmask_fitvar_rad(mask_fitvar_rad(k:nf-1))=rmask_fitvar_rad(mask_fitvar_rad(k:nf-1)) - 1
        mask_fitvar_rad(nf) = 0
        nfalb = nfalb - 1; nf = nf - 1; k = k - 1

        lo_radbnd(j) = fitvar_rad_init(j)
        up_radbnd(j) = fitvar_rad_init(j)
        rmask_fitvar_rad(j) = 0
           !WRITE(*, '(4I5, 2F8.3, A10)') i,  j, k, nfalb, albmin(i), albmax(i)
           !, fitvar_rad_str(j)
      ENDIF
    ENDDO    
    IF (nalb < 1) THEN
      WRITE(www_lun, *) 'No valid albedo is specified!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF

    hasalb = 0
    sfcalbs(1:npoints, :) = 0.0d0
    DO i = 1, nalb
      j = albidx + i -1
      albfpix(i)= MINVAL(MINLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
           >= albmin(i) .AND. fitwavs(1:npoints) < albmax(i))))
      alblpix(i)= MINVAL(MAXLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
           >= albmin(i) .AND. fitwavs(1:npoints) < albmax(i))))
      IF (fitvar_rad_str(j)(4:4) == '0' .AND. .NOT. is_albspcvar(i) ) THEN
           hasalb(albfpix(i):alblpix(i)) = hasalb(albfpix(i):alblpix(i)) + 1
           sfcalbs(albfpix(i):alblpix(i), 1) = fitvar_rad_init(j)
      ELSE IF (is_albspcvar(i)) THEN 
          hasalb(albfpix(i):alblpix(i)) = 1
          sfcalbs(albfpix(i):alblpix(i), 1) = albspcs(albfpix(i):alblpix(i), 0) ! italize with MODIS-derived spectrum
      ENDIF
    ENDDO

    IF (ANY(hasalb(1:npoints) == 0)) THEN
      WRITE(www_lun, *) 'Albedo is not specified for all wavelengths!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ELSE IF (ANY(hasalb(1:npoints) > 1)) THEN
      WRITE(www_lun, *) &
           'Multiple albedos are specified for some wavelengths!!!'
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF


    DO i =  1, nwfc
      j = wfcidx + i - 1
      IF  (fitvar_rad_str(j)(4:4) =='0') THEN
        fitvar_rad_init(j) = cfrac
        IF (up_radbnd(j) == lo_radbnd(j)) THEN
          up_radbnd(j) = cfrac
          lo_radbnd(j) = cfrac
        ENDIF
      ENDIF
    ENDDO

  ! Go thorugh wfc terms again to check for unused
  IF (nwfc > 0) THEN
     i = 1; k = wfcfidx - 1
     DO i = 1, nwfc
        j = wfcidx + i - 1
        IF (lo_radbnd(j) < up_radbnd(j))  k = k + 1
        
        nsub = COUNT(MASK=(fitwavs(1:npoints) >= wfcmin(i) .AND. fitwavs(1:npoints) <= wfcmax(i)))
        IF (lo_radbnd(j) < up_radbnd(j) .AND. nsub <= 0)  THEN                      
           mask_fitvar_rad(k:nf-1) = mask_fitvar_rad(k+1:nf)
           rmask_fitvar_rad(mask_fitvar_rad(k:nf-1))= rmask_fitvar_rad(mask_fitvar_rad(k:nf-1)) - 1
           mask_fitvar_rad(nf) = 0          
           nfwfc = nfwfc - 1; nf = nf - 1; k = k - 1

           lo_radbnd(j) = fitvar_rad_init(j)
           up_radbnd(j) = fitvar_rad_init(j)
           rmask_fitvar_rad(j) = 0   
        ENDIF           
     ENDDO
  ENDIF


    IF ( nwfc > 0 ) THEN
      haswfc = 0
      DO i = 1, nwfc
        j = wfcidx + i -1
        wfcfpix(i)= MINVAL(MINLOC(fitwavs(1:npoints), &
             MASK=(fitwavs(1:npoints) &
             >= wfcmin(i) .AND. fitwavs(1:npoints) < wfcmax(i))))
        wfclpix(i)= MINVAL(MAXLOC(fitwavs(1:npoints), &
             MASK=(fitwavs(1:npoints) &
             >= wfcmin(i) .AND. fitwavs(1:npoints) < wfcmax(i))))
        IF (fitvar_rad_str(j)(4:4) == '0') &
             haswfc(wfcfpix(i):wfclpix(i)) = haswfc(wfcfpix(i):wfclpix(i)) + 1
      ENDDO

      IF (ANY(haswfc(1:npoints) == 0)) THEN
        WRITE(www_lun, *) &
             'Cloud fraction is not specified for all wavelengths!!!'
        pge_error_status = pge_errstat_error
        RETURN
      ELSE IF (ANY(haswfc(1:npoints) > 1)) THEN
        WRITE(www_lun, *) &
             'Multiple cloud fractions are specified for some wavelengths!!!'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF
    ENDIF


    RETURN

  END SUBROUTINE set_cldalb


  !  Unused?
  !
  !  SUBROUTINE SETCLD(npoints, fitwavs, nscene, ctaus, ctps, cfracs, errstat)
  !
  !    USE OMSAO_precision_module
  !    USE OMSAO_variables_module,    ONLY: the_month, the_year, the_day!, &
  !         !the_lon, the_lat, b1ab_div_wav
  !    USE ozprof_data_module,        ONLY: cloud!, do_ch2reso
  !    USE OMSAO_gome_data_module,    ONLY: gome_orbc, gome_curpix, &
  !         gome_curscan, gome_npix!, gome_pixdate,
  !    USE OMSAO_errstat_module
  !    use get_cloud, only: get_cloud_maprop
  !
  !    IMPLICIT NONE
  !
  !    ! ========================
  !    ! Input/Output Variables
  !    ! ========================
  !    INTEGER, INTENT(IN)                 :: npoints, nscene
  !    INTEGER, INTENT(OUT)                :: errstat
  !    REAL(KIND=dp), DIMENSION(npoints)   :: fitwavs
  !    REAL (KIND=dp), DIMENSION (nscene), INTENT(OUT) :: ctaus, ctps, cfracs
  !    ! 0.375s and 1.2s
  !
  !    ! ===============
  !    ! Local variables
  !    ! ===============
  !    INTEGER :: stpix, sndpix, cldflg!, i, j, k
  !
  !    errstat = pge_errstat_ok
  !    ctaus  = 0.0
  !    cfracs = 0.0
  !    ctps   = 0.0
  !
  !    IF (cloud)  THEN
  !      IF (gome_npix == 1 .AND. gome_curscan < 3) THEN
  !        stpix = gome_curpix
  !        cldflg = 1
  !      ELSE IF (gome_npix == 1 .AND. gome_curscan == 3) THEN
  !        stpix = gome_curpix - 3
  !        cldflg = 2
  !      ELSE IF (gome_npix == 8) THEN
  !        stpix = gome_curpix - 7
  !        cldflg = 0
  !      ELSE IF (gome_npix == 40) THEN
  !        stpix = gome_curpix - 39
  !        cldflg = 5
  !      ELSE
  !        WRITE(www_lun, *) 'SETCLD: Should not happen'
  !        errstat = pge_errstat_error
  !        RETURN
  !      ENDIF
  !
  !      ! always assume the same cloud info for channel 1 and channel 2
  !      ! current algorithm cannot handle different cloud information for
  !      ! chs 1 & 2
  !      sndpix = 0
  !
  !      !WRITE(www_lun,'(3I5,A4,3I5, 2F8.2)') the_year, the_month, the_day, &
  !      !     gome_orbc, stpix, sndpix, cldflg, the_lon, the_lat
  !
  !      CALL GET_CLOUD_MAPROP(the_year, the_month, the_day, gome_orbc, &
  !           stpix, sndpix, cldflg, ctaus, cfracs, ctps, errstat)
  !    ENDIF
  !
  !    RETURN
  !  END SUBROUTINE SETCLD

end module m_set_cldalb
