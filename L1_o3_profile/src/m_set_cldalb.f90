! ===================================================================
!	        Set up clouds and surface albedo  (pain in the neck)
! ===================================================================     
! Search for cloud information and surface albedo
! If has clouds, use GOME derived albedo by Kolemeijer, then use fixed
! albedoes, need to adjust fitvar_rad and fitvar
! If there no clouds, then derive the surface albedo from 370.2 nm
! Need to override the specified albedo values 
! Need to check for no albedo or multiple albedo specified at some wavelengths
module m_set_cldalb
  USE m_get_initial_albedo, only: adj_albcfrac, get_gome_alb, &
        get_initial_albedo, get_omi_alb, get_omler_alb, get_omler_albs,get_toms_alb, &
        get_sciagm2_alb, get_surface_spectrum  

  PUBLIC set_cldalb
  INTEGER (kind=1) :: oceanflg, snowflg
  CONTAINS

  SUBROUTINE set_snowoceanflg (the_snowice, the_landocean)
  USE OMSAO_indices_module, ONLY: instrument_idx, omi_idx, tempo_idx, gome2_idx
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: the_snowice, the_landocean  ! given from l1b product

  snowflg = 0; oceanflg=0
  IF (instrument_idx == omi_idx .or. instrument_idx == tempo_idx) THEN 
   IF (the_snowice >= 10 .and. the_snowice <= 103 ) THEN
      snowflg = 1
   ELSE
      snowflg = 0
   ENDIF
   IF (the_landocean /=1 )  oceanflg = 1
  ELSE IF (instrument_idx == gome2_idx) THEN 
    IF (the_snowice >= 10) THEN 
      snowflg = 1 
      WRITE(*,*) 'snowflag should be checked for gome2' ; stop
    ENDIF
  ELSE
    WRITE(*,*) 'snowflag should be checked for this instrument' ; stop
  ENDIF
  RETURN
  END SUBROUTINE set_snowoceanflg

  SUBROUTINE set_cldalb (npoints, fitwavs, ctau, ctp, cfrac, salbedo, pge_error_status)
     
  USE OMSAO_precision_module
  USE OMSAO_variables_module,    ONLY: fitvar_rad_init, fitvar_rad_str, &
       lo_radbnd, up_radbnd, nf=>n_fitvar_rad, mask_fitvar_rad, &
       the_month, the_day, edgelons, edgelats, the_jday, &
       the_sza_atm, the_vza_atm, the_aza_atm, & 
       numwin, nviswin, widx_rvis,nradpix, rmask_fitvar_rad, wcenter_uvvis
  USE ozprof_data_module,        ONLY: albidx, albfidx, nalb, nfalb, albmin, &
       albmax, albfpix, alblpix, albfpix_r, alblpix_r, & 
        do_lambcld, lambcld_refl, which_alb, &
       the_fixalb, do_simu, radcalwrt, lambcld_initalb, ecfrfind, & !has_glint,glintprob
       wfcmax, wfcmin, nwfc, nfwfc, wfcfpix, wfclpix, wfcfidx, wfcidx, the_snowice, pos_alb, &
       do_brdf, use_albspc, use_albeofs, which_albspc ,the_landwater_flg,  & 
       albspcs, malbspc, nactalbspc, is_albspcvar, sfcalbs, nrefl, &
       albspcs_hres, use_effcrs, ncalcp, radcwav, the_landfrac !nalbspc
  USE OMSAO_errstat_module
  USE m_set_brdf, ONLY: set_brdf, Surface
  USE m_ezspline_interpolation, ONLY:bspline
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
  INTEGER                      :: i, j, k,  ninalb, region, nw, fidx, lidx, ntmp, & 
                                  nsub, which_sciagm2, nalbw
  INTEGER, DIMENSION(npoints)  :: hasalb, haswfc
  REAL (KIND=dp)               :: cfrac_old, albedo, wavg
  INTEGER, DIMENSION(2)        :: wavin
  REAL (KIND=dp), DIMENSION(2) :: wavfrac
  ! GOME albedo database: 335, 380, 440, 495, 555, 610, 670 nm
  ! OMI albedo database:  328, 345, 354, 380, 442, 477, 499 nm
  ! For GOME-2 MSC/PMD and SCIA, wavelengths are read from the data
  INTEGER, PARAMETER                 :: nalbw0 = 29
  REAL (KIND=dp), DIMENSION (nalbw0) :: albarr, albwave
  LOGICAL                            :: noalb, do_adjcfrac       
  INTEGER (KIND=i4), EXTERNAL       :: day_of_year       
  ! ===============
  ! module name
  ! ===============
  LOGICAL :: do_debug
  CHARACTER(LEN=10), PARAMETER      :: modulename = 'SET_CLDALB'

  pge_error_status = pge_errstat_ok

  nalbw = 7
  albwave(1:nalbw) = (/335.0, 380.0, 440.0, 495., 555., 610., 670.0/)
  IF (which_alb == 1) THEN
     !ninalb = 1; region = 1
     ninalb = 7; region = 3
     CALL GET_GOME_ALB(the_month, edgelons, edgelats, region, albarr, ninalb)
  ELSE IF (which_alb == 2) THEN
     CALL GET_TOMS_ALB(the_month, edgelons, edgelats, albedo)
     albarr(1:nalbw) = albedo
  ELSE IF (which_alb == 3) THEN
     CALL GET_OMI_ALB(the_month, the_day, edgelons, edgelats, albedo)
     albarr(1:nalbw) = albedo
  ELSE IF (which_alb == 4) THEN
     nalbw = 7
     albwave(1:nalbw) = (/328.0, 345.0, 354.0, 380.0, 442.0, 477.0, 499.0/)
     CALL GET_OMLER_ALBS(the_month, the_day, edgelons, edgelats, albarr)
     albarr(1:nalbw) = albarr
  ELSE IF (which_alb >= 5 .and. which_alb <=7) THEN
     which_sciagm2 = which_alb - 4 
     CALL GET_SCIAGM2_ALB(which_sciagm2, the_month, the_day, edgelons, edgelats, albarr, albwave, nalbw)
     DO i = 1, nalbw
        print *, albwave(i), albarr(i)
     ENDDO
     STOP
  ELSE
     WRITE(*, *) 'Albedo database: not implemented!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
 
  CALL set_snowoceanflg(the_snowice,the_landwater_flg) 
  IF (use_albspc .or. do_brdf) THEN
    fidx = SUM(nradpix(1 : numwin - nviswin)) + 1
    lidx = npoints
    ntmp = lidx - fidx + 1
    IF (allocated(albspcs)) deallocate (albspcs)
    allocate (albspcs(npoints, 0:malbspc-1))
    albspcs = 0.0
    IF (do_brdf) THEN 
       CALL set_brdf (npoints, fitwavs(1:npoints), nactalbspc,the_landfrac, pge_error_status)
    ELSE
      IF (which_albspc == 1 .or. snowflg == 1 .or. oceanflg == 1) then ! from peter  
        CALL get_surface_spectrum (the_jday, the_sza_atm, the_aza_atm, the_vza_atm, edgelons, edgelats, &
            snowflg,  ntmp, fitwavs(fidx:lidx), albspcs(fidx:lidx,0:malbspc-1), nactalbspc, the_landfrac, pge_error_status)
        IF (snowflg == 1) THEN
          albspcs(fidx:lidx, 0) = albspcs(fidx:lidx, 0) * the_snowice / 100.0
        ENDIF
        IF (the_landfrac == 0.0) nactalbspc = 1
        print * , 'get_surface_spectrum(nactalbspc,the_landfrac', nactalbspc,the_landfrac, snowflg
      ELSE  IF (which_albspc == 2) THEN 
        CALL SET_BRDF(ntmp, fitwavs(fidx:lidx), nactalbspc, the_landfrac, pge_error_status)
        IF (Surface%Option4%Wvl(1) > fitwavs(fidx) .or. & 
           Surface%Option4%Wvl(Surface%Option4%wmx) <  fitwavs(lidx) ) THEN
           WRITE(*,'(A)') modulename//':Check the boundaries of BRDF spectrum !!!' ; RETURN
        ENDIF 
        CALL BSPLINE(Surface%Option4%Wvl,Surface%Option4%Mu(:),Surface%Option4%wmx, fitwavs(fidx:lidx), & 
             albspcs(fidx:lidx,0),ntmp, pge_error_status)
        DO i=1,Surface%Option4%fmx
          CALL BSPLINE(Surface%Option4%Wvl, Surface%Option4%W(:,i),Surface%Option4%wmx, &
          fitwavs(fidx:lidx), albspcs(fidx:lidx,i),ntmp, pge_error_status)
        ENDDO
      ENDIF
    ENDIF

     !WRITE(999,*) ntmp
     !DO i  = fidx, lidx
     !WRITE(999, '(10e15.7)') fitwavs(i), albspcs(i, 0:malbspc-1)
     !ENDDO
     
    IF (.not. use_effcrs) THEN 
      IF (allocated(albspcs_hres)) deallocate(albspcs_hres)
      allocate (albspcs_hres(ncalcp, 0:malbspc-1))
      fidx=widx_rvis ; lidx=ncalcp
      ntmp = lidx-fidx+1
      IF (do_brdf) THEN 
        CALL set_brdf (ncalcp, radcwav(1:ncalcp), nactalbspc, the_landfrac, pge_error_status)
      ELSE
        albspcs_hres = 0.0
        IF (which_albspc == 1 .or. snowflg == 1 .or. oceanflg == 1) THEN 
          CALL get_surface_spectrum (the_jday, the_sza_atm, the_aza_atm, the_vza_atm,edgelons, edgelats, &
          snowflg,  ntmp, radcwav(fidx:lidx), albspcs_hres(fidx:lidx,0:malbspc-1), nactalbspc,the_landfrac, pge_error_status)
          IF (snowflg == 1) albspcs_hres(fidx:lidx, 0) = albspcs_hres(fidx:lidx,0) *the_snowice/100.0
          IF (the_landfrac == 0) nactalbspc =1
        ELSE IF (which_albspc == 2) THEN 
          IF (Surface%Option4%Wvl(1) > radcwav(fidx) .or. & 
             Surface%Option4%Wvl(Surface%Option4%wmx) <  radcwav(lidx) ) THEN
             WRITE(*,'(A)') modulename//':Check the boundaries of BRDF spectrum !!!' ; RETURN
          ENDIF 
          CALL BSPLINE(Surface%Option4%Wvl,Surface%Option4%Mu(:),Surface%Option4%wmx, radcwav(fidx:lidx), & 
               albspcs_hres(fidx:lidx,0),ntmp, pge_error_status)
          DO i=1,Surface%Option4%fmx
            CALL BSPLINE(Surface%Option4%Wvl, Surface%Option4%W(:,i),Surface%Option4%wmx, &
            radcwav(fidx:lidx), albspcs_hres(fidx:lidx,i),ntmp, pge_error_status)
          ENDDO
        ENDIF
      ENDIF
       !WRITE(999, *)  ntmp     
       !DO i  = fidx, lidx
       !   WRITE(999, *) radcwav(i),  albspcs_hres(i, 0:malbspc-1)
       !ENDDO
    ENDIF
    IF (pge_error_status == pge_errstat_error) THEN
      WRITE(*, *) modulename, ' Error in getting surface albedo spectrum'; RETURN
    ENDIF
  ENDIF

  ! get initial albedo
  IF (nrefl > 0) THEN
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
           nw = 1; wavin(1) = 1; wavfrac(1) = 1.0
        ELSE IF ( i == nalbw .AND. pos_alb >= albwave(nalbw) ) THEN
           nw = 1; wavin(1) = nalbw; wavfrac(1) = 1.0
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
  
  IF (do_simu .AND. .NOT. radcalwrt .AND. the_fixalb >= 0.0) albedo = the_fixalb

  ! adjust the cloud fraction or surface albedo based on 
  ! LIDORT-calculated radiance at 370.2 nm
  !CALL ADJ_ALBCFRAC(albedo, cfrac, pge_error_status)
  !IF (pge_error_status == pge_errstat_error)  RETURN

  ! Initialize Lambertian cloud albedo to be 80%
  IF ( .NOT. (do_simu .AND. .NOT. radcalwrt) ) THEN
     IF (do_lambcld ) lambcld_refl = lambcld_initalb
  ENDIF
 
   ! For simulation, if cfrac == 1.0, then need to adjust lambcld_refl.  However, this already partly corrects 
   ! for calibration offset by deriving a lambcld_refl that matches the measured radiance at cloud wavelength
   ! IF (.NOT. do_simu .OR. (cfrac == 1.0D0 .AND. radcalwrt) ) 
   !IF (the_snowice == 101) THEN
   !   albedo = 0.90
   !ELSE IF (the_snowice == 103) THEN
   !   albedo = 0.80
   !ELSE IF (the_snowice > 1 .AND. the_snowice <= 100) THEN
   !   albedo = MAX(albedo, 0.8 * the_snowice / 100.0)
   !ENDIF

  ! xliu, 12/07/2014, changed based on ASTER snow spectrum
   IF (the_snowice > 100 .and. snowflg == 1 ) THEN  ! permanently covered
      albedo = 0.98
      albarr = albedo
   ELSE IF (the_snowice > 1 .AND. the_snowice <= 100) THEN ! partialy covered
      albedo = MAX(albedo, 0.98 * the_snowice / 100.0)
      albarr = albedo
   ENDIF
  
   do_adjcfrac = .TRUE.
   IF (nrefl < 1) do_adjcfrac = .FALSE.
   IF ( do_adjcfrac ) THEN
     cfrac_old = cfrac  ! save cloud fraction from other products
     CALL ADJ_ALBCFRAC(albedo, cfrac, ctau, pge_error_status) 
   ENDIF
   salbedo = albedo   ! Surface albedo
   ! Test using FRECO cld fraction 
   !cfrac = cfrac_old 
 
   ! xliu: 08/16/2008, when surface albedo increases, it is more difficult to differentiate clouds/surfaces
   ! Assume a cloud fraction of 0 and increases a priori error for surface albedo and cloud fraction
   !IF (albedo > 0.6 .AND. cfrac >= 0.6 ) cfrac = 0.5

   IF (ecfrfind == 0 .AND. nwfc == 0) THEN ! Final cloud fraction is computed here  
     IF (ctp > 0.0) THEN  
        IF (do_lambcld .AND. cfrac >= 1.0D0) albedo = lambcld_refl
     ENDIF
  ELSE
     ! The derived cloud fraction is derived as initial value
     ! Slightly change the value for clear-sky/cloud-sky, so weighting function 
     ! for clouds are calculated
     IF ( cfrac == 0.0) cfrac = 0.01
     IF ( cfrac == 1.0) cfrac = 0.99
     !IF ( has_glint .AND. cfrac < 0.20 * glintprob ) cfrac = 0.0D0
  ENDIF
  
  ! copy initial albedo to fitvar_rad_init   
  DO i =  1, nalb
     j = albidx + i - 1
     IF  (fitvar_rad_str(j)(4:4) =='0' .AND. .NOT. is_albspcvar(i) ) THEN
        IF (which_alb == 1 .OR. which_alb == 4) THEN
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
           up_radbnd(j) = albedo; lo_radbnd(j) = albedo
        ENDIF
     ELSE  
        ! at specfit_ozprof 
     ENDIF
  ENDDO
 
  ! Go thorugh albedo terms again to check for unused 
  ! 04/21/2016, updated 

  !i = 1; k = albfidx - 1
  !DO i = 1, nalb
  !   j = albidx + i - 1
  !   IF (lo_radbnd(j) < up_radbnd(j))  k = k + 1
  !   write(*, '(6I5, F6.3, A10, 2F8.3)') i, j, k, nf, mask_fitvar_rad(k), rmask_fitvar_rad(j), &
  !        fitvar_rad_init(j), fitvar_rad_str(j), albmin(i), albmax(i)
  !ENDDO  

  IF (nalb > 0) THEN
     i = 1; k = albfidx - 1
     DO i = 1, nalb
        j = albidx + i - 1
        IF (lo_radbnd(j) < up_radbnd(j))  k = k + 1
        nsub = COUNT(MASK=(fitwavs(1:npoints) >= albmin(i) .AND. fitwavs(1:npoints) <= albmax(i)))
        IF (lo_radbnd(j) < up_radbnd(j) .AND. nsub <= 0)  THEN               
           mask_fitvar_rad(k:nf-1) = mask_fitvar_rad(k+1:nf)
           rmask_fitvar_rad(mask_fitvar_rad(k:nf-1))= rmask_fitvar_rad(mask_fitvar_rad(k:nf-1)) - 1
           mask_fitvar_rad(nf) = 0        
           nfalb = nfalb - 1; nf = nf - 1; k = k - 1

           lo_radbnd(j) = fitvar_rad_init(j)
           up_radbnd(j) = fitvar_rad_init(j)
           rmask_fitvar_rad(j) = 0          
           !WRITE(*, '(4I5, 2F8.3, A10)') i,  j, k, nfalb, albmin(i), albmax(i) , fitvar_rad_str(j)
        ENDIF
     ENDDO
  ENDIF
  
  IF (nalb < 1) THEN
     WRITE(*, *) 'No valid albedo is specified!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  
  hasalb = 0
  IF (allocated(sfcalbs)) deallocate(sfcalbs)
  allocate (sfcalbs(npoints, 2))
  DO i = 1, nalb 
     j = albidx + i -1
     albfpix(i)= MINVAL(MINLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
          >= albmin(i) .AND. fitwavs(1:npoints) < albmax(i)))) 
     alblpix(i)= MINVAL(MAXLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
          >= albmin(i) .AND. fitwavs(1:npoints) < albmax(i))))  
     IF (.NOT. use_effcrs) THEN 
      albfpix_r(i)= MINVAL(MINLOC(radcwav(1:ncalcp), MASK=(radcwav(1:ncalcp) &
          >= albmin(i) .AND. radcwav(1:ncalcp) < albmax(i)))) 
      alblpix_r(i)= MINVAL(MAXLOC(radcwav(1:ncalcp), MASK=(radcwav(1:ncalcp) &
          >= albmin(i) .AND. radcwav(1:ncalcp) < albmax(i))))  
     ENDIF
     IF (fitvar_rad_str(j)(4:4) == '0' .AND. .NOT. is_albspcvar(i)) THEN
        hasalb(albfpix(i):alblpix(i)) = hasalb(albfpix(i):alblpix(i)) + 1
     ELSE IF (is_albspcvar(i)) THEN
        hasalb(albfpix(i):alblpix(i)) = 1
     ENDIF        
  ENDDO
  IF (ANY(hasalb(1:npoints) == 0)) THEN
     WRITE(*, *) 'Albedo is not specified for all wavelengths!!!'
     pge_error_status = pge_errstat_error; RETURN
  ELSE IF (ANY(hasalb(1:npoints) > 1)) THEN
     WRITE(*, *) 'Multiple albedos are specified for some wavelengths!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  !DO i = 1, npoints
  !   WRITE(85, '(F8.4, 4F12.8)') fitwavs(i), sfcalbs(i, 1)
  !ENDDO
    
  !DO i = 1, nalb
  !   j = albidx +i - 1
  !   WRITE(*,'(A5,3D12.3,2F8.2,2I5)') fitvar_rad_str(j), fitvar_rad_init(j), &
  !        lo_radbnd(j), up_radbnd(j), albmin(i), albmax(i), &
  !        albfpix(i), alblpix(i)
  !ENDDO

  DO i =  1, nwfc
     j = wfcidx + i - 1
     IF  (fitvar_rad_str(j)(4:4) =='0') THEN
        fitvar_rad_init(j) = cfrac 
        IF (up_radbnd(j) == lo_radbnd(j)) THEN
           up_radbnd(j) = cfrac; lo_radbnd(j) = cfrac
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
  !IF (nwfc < 1) THEN
  !   WRITE(*, *) 'No valid wavelength-dependent cloud fraction is specified!!!'
  !   pge_error_status = pge_errstat_error; RETURN
  !ENDIF
  
  IF ( nwfc > 0 ) THEN 
     haswfc = 0
     DO i = 1, nwfc 
        j = wfcidx + i -1
        wfcfpix(i)= MINVAL(MINLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
              >= wfcmin(i) .AND. fitwavs(1:npoints) < wfcmax(i)))) 
        wfclpix(i)= MINVAL(MAXLOC(fitwavs(1:npoints), MASK=(fitwavs(1:npoints) &
             >= wfcmin(i) .AND. fitwavs(1:npoints) < wfcmax(i))))  
        IF (fitvar_rad_str(j)(4:4) == '0') &
             haswfc(wfcfpix(i):wfclpix(i)) = haswfc(wfcfpix(i):wfclpix(i)) + 1
     ENDDO
     
     IF (ANY(haswfc(1:npoints) == 0)) THEN
        WRITE(*, *) 'Cloud fraction is not specified for all wavelengths!!!'
        pge_error_status = pge_errstat_error; RETURN
     ELSE IF (ANY(haswfc(1:npoints) > 1)) THEN
        WRITE(*, *) 'Multiple cloud fractions are specified for some wavelengths!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ENDIF
  RETURN
  
END SUBROUTINE set_cldalb

end module m_set_cldalb
