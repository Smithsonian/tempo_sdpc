!
module m_set_cldalb

  public set_cldalb
  private setcld

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
    USE OMSAO_variables_module, ONLY: fitvar_rad_init, fitvar_rad_str, &
         lo_radbnd, up_radbnd, nf=>n_fitvar_rad, mask_fitvar_rad, &
         the_month, the_day, edgelons, edgelats!, the_lon, the_lat, &
         !the_year, b1ab_div_wav
    USE ozprof_data_module, ONLY: albidx, albfidx, nalb, nfalb, albmin, &
         albmax, albfpix, alblpix, do_lambcld, lambcld_refl, &
         which_alb, the_fixalb, do_simu, radcalwrt, lambcld_initalb, &
         ecfrfind, wfcmax, wfcmin, nwfc, nfwfc, &
         wfcfpix, wfclpix, wfcfidx, wfcidx, the_snowice!, ps0, has_glint, glintprob, cloud
    !USE OMSAO_gome_data_module, ONLY: gome_pixdate, gome_orbc, gome_curpix, &
    !     gome_curscan
    USE OMSAO_errstat_module
    use m_get_initial_albedo, only: adj_albcfrac, get_gome_alb, &
         get_initial_albedo, get_omi_alb, get_omler_alb, get_toms_alb

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
    INTEGER :: i, j, k, ninalb, region
    INTEGER, DIMENSION(npoints)  :: hasalb, haswfc
    REAL (KIND=dp)               :: cfrac_old, albedo!, actp
    REAL (KIND=dp), DIMENSION (5):: albarr ! 335, 453, 525, 610, 670 nm
    LOGICAL                      :: noalb               

    pge_error_status = pge_errstat_ok

    IF (which_alb == 1) THEN
      ninalb = 1
 region = 1
      CALL GET_GOME_ALB(the_month, edgelons, edgelats, region, albarr, ninalb)
    ELSE IF (which_alb == 2) THEN
      CALL GET_TOMS_ALB(the_month, edgelons, edgelats, albedo)
      albarr(:) = albedo
    ELSE IF (which_alb == 3) THEN
      CALL GET_OMI_ALB(the_month, the_day, edgelons, edgelats, albedo)
      albarr(:) = albedo
    ELSE IF (which_alb == 4) THEN
      CALL GET_OMLER_ALB(the_month, the_day, edgelons, edgelats, albedo)
      albarr(:) = albedo
    ELSE
      WRITE(www_lun, *) 'Albedo database: not implemented!!!'
      pge_error_status = pge_errstat_error
 RETURN
    ENDIF

    ! get initial albedo
    IF (ctp > 0.0) THEN
      noalb = .TRUE.

      CALL GET_INITIAL_ALBEDO(noalb, albedo, pge_error_status) 
      albedo = albarr(1)

      IF (pge_error_status == pge_errstat_error)  RETURN
    ELSE
      ! get effective surface albedo
      noalb = .FALSE.
      CALL GET_INITIAL_ALBEDO(noalb, albedo, pge_error_status)     
      IF (pge_error_status == pge_errstat_error)  RETURN

      IF (albedo > 1.0) THEN
        WRITE(www_lun, *) &
             'Thick clouds exist. No retrieval done or add cloud info!!!'
        pge_error_status = pge_errstat_error
 RETURN
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
    IF (the_snowice == 101) THEN
      albedo = 0.90
    ELSE IF (the_snowice == 103) THEN
      albedo = 0.80
    ELSE IF (the_snowice > 1 .AND. the_snowice <= 100) THEN
      albedo = MAX(albedo, 0.8d0 * DBLE(the_snowice) / 100.d0)
    ENDIF
    cfrac_old = cfrac  ! save cloud fraction from other products
    CALL ADJ_ALBCFRAC(albedo, cfrac, ctau, pge_error_status) 
    salbedo = albedo   ! Surface albedo

    ! xliu: 08/16/2008, when surface albedo increases, it is more difficult 
    ! to differentiate clouds/surfaces
    ! Assume a cloud fraction of 0 and increases a priori error for surface 
    ! albedo and cloud fraction
    IF (albedo > 0.6 .AND. cfrac >= 0.6 .AND. .NOT. do_simu) cfrac = 0.5

    IF (ecfrfind == 0 .AND. nwfc == 0) THEN ! Final cloud fraction is computed here  
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
      IF  (fitvar_rad_str(j)(4:4) =='0') THEN
        fitvar_rad_init(j) = albedo 
        IF (up_radbnd(j) == lo_radbnd(j)) THEN
          up_radbnd(j) = albedo
 lo_radbnd(j) = albedo
        ENDIF
      ENDIF
    ENDDO

    ! Go thorugh albedo terms again to check for unused 
    i = 1
 k = albfidx - 1
    DO 
      j = albidx + i - 1
      IF (lo_radbnd(j) < up_radbnd(j))  k = k + 1

      IF (albmax(i) <= MINVAL(fitwavs(1:npoints)) .OR. &
           albmin(i) >= MAXVAL(fitwavs(1:npoints))) THEN          
        IF (i < nalb) THEN
          fitvar_rad_init(j) = fitvar_rad_init(j+1)
          fitvar_rad_str(j) = fitvar_rad_str(j+1)
          albmin(i) = albmin(i+1)
 albmax(i) = albmax(i+1)
        ENDIF
        nalb = nalb - 1

        IF (lo_radbnd(j) < up_radbnd(j)) THEN
          mask_fitvar_rad(k:nf-1) = mask_fitvar_rad(k+1:nf)
          nfalb = nfalb - 1
 nf = nf - 1
        ENDIF

      ELSE
        i = i + 1
      ENDIF

      IF (i >= nalb) EXIT       
    ENDDO

    IF (nalb < 1) THEN
      WRITE(www_lun, *) 'No valid albedo is specified!!!'
      pge_error_status = pge_errstat_error
 RETURN
    ENDIF

    hasalb = 0
    DO i = 1, nalb 
      j = albidx + i -1
      albfpix(i)= MINVAL(MINLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
           >= albmin(i) .AND. fitwavs(1:npoints) < albmax(i)))) 
      alblpix(i)= MINVAL(MAXLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
           >= albmin(i) .AND. fitwavs(1:npoints) < albmax(i))))  
      IF (fitvar_rad_str(j)(4:4) == '0') &
           hasalb(albfpix(i):alblpix(i)) = hasalb(albfpix(i):alblpix(i)) + 1
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
    i = 1
 k = wfcfidx - 1
    DO 
      j = wfcidx + i - 1
      IF (lo_radbnd(j) < up_radbnd(j))  k = k + 1

      IF (wfcmax(i) <= MINVAL(fitwavs(1:npoints)) .OR. &
           wfcmin(i) >= MAXVAL(fitwavs(1:npoints))) THEN          
        IF (i < nwfc) THEN
          fitvar_rad_init(j) = fitvar_rad_init(j+1)
          fitvar_rad_str (j) = fitvar_rad_str (j+1)
          wfcmin(i) = wfcmin(i+1)
 wfcmax(i) = wfcmax(i+1)
        ENDIF
        nwfc = nwfc - 1

        IF (lo_radbnd(j) < up_radbnd(j)) THEN
          mask_fitvar_rad(k:nf-1) = mask_fitvar_rad(k+1:nf)
          nfwfc = nfwfc - 1
 nf = nf - 1
        ENDIF

      ELSE
        i = i + 1
      ENDIF

      IF (i >= nwfc) EXIT       
    ENDDO


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
  SUBROUTINE SETCLD(npoints, fitwavs, nscene, ctaus, ctps, cfracs, errstat)

    USE OMSAO_precision_module
    USE OMSAO_variables_module,    ONLY: the_month, the_year, the_day!, &
         !the_lon, the_lat, b1ab_div_wav
    USE ozprof_data_module,        ONLY: cloud!, do_ch2reso
    USE OMSAO_gome_data_module,    ONLY: gome_orbc, gome_curpix, &
         gome_curscan, gome_npix!, gome_pixdate, 
    USE OMSAO_errstat_module
    use get_cloud, only: get_cloud_maprop

    IMPLICIT NONE

    ! ========================
    ! Input/Output Variables
    ! ========================
    INTEGER, INTENT(IN)                 :: npoints, nscene
    INTEGER, INTENT(OUT)                :: errstat
    REAL(KIND=dp), DIMENSION(npoints)   :: fitwavs
    REAL (KIND=dp), DIMENSION (nscene), INTENT(OUT) :: ctaus, ctps, cfracs 
    ! 0.375s and 1.2s

    ! ===============
    ! Local variables
    ! ===============
    INTEGER :: stpix, sndpix, cldflg!, i, j, k

    errstat = pge_errstat_ok
    ctaus  = 0.0
    cfracs = 0.0
    ctps   = 0.0

    IF (cloud)  THEN
      IF (gome_npix == 1 .AND. gome_curscan < 3) THEN  
        stpix = gome_curpix
        cldflg = 1      
      ELSE IF (gome_npix == 1 .AND. gome_curscan == 3) THEN
        stpix = gome_curpix - 3
        cldflg = 2
      ELSE IF (gome_npix == 8) THEN
        stpix = gome_curpix - 7
        cldflg = 0
      ELSE IF (gome_npix == 40) THEN
        stpix = gome_curpix - 39
        cldflg = 5
      ELSE
        WRITE(www_lun, *) 'SETCLD: Should not happen'
        errstat = pge_errstat_error
        RETURN
      ENDIF

      ! always assume the same cloud info for channel 1 and channel 2
      ! current algorithm cannot handle different cloud information for 
      ! chs 1 & 2
      sndpix = 0                 

      !WRITE(www_lun,'(3I5,A4,3I5, 2F8.2)') the_year, the_month, the_day, &
      !     gome_orbc, stpix, sndpix, cldflg, the_lon, the_lat

      CALL GET_CLOUD_MAPROP(the_year, the_month, the_day, gome_orbc, &
           stpix, sndpix, cldflg, ctaus, cfracs, ctps, errstat)
    ENDIF

    RETURN
  END SUBROUTINE SETCLD

end module m_set_cldalb
