module m_read_ozprof_input
  ! ***************************************************************************
  ! Author:  xiong liu
  ! Date  :  July 23, 2003
  ! Purpose: read input ozone profile variables and their bounds and add
  !          those variable to fitting control variables for radiance fit
  ! ***************************************************************************

CONTAINS
SUBROUTINE read_ozprof_input (fit_ctrl_unit, fit_ctrl_file, pge_error_status )     

  USE OMSAO_precision_module
  USE ozprof_data_module, ONLY: algorithm_name, algorithm_version, &
         do_multi_vza, do_radinter, do_simu, do_simu_rmring, radcalwrt, &
         which_caloz, caloz_fname, use_oe, use_lograd, use_logstate, use_flns, &
         ring_on_line, ring_convol, fit_atanring, ring_lut, & 
         do_twostep, do_bothstep, use_large_so2_aperr, do_tracewf, do_subfit,& 
         rtm_treatment, atmwrt, atmos_prof_fname, &
         gaswrt, ozwrtcorr, ozwrtcovar, ozwrtcontri, ozwrtres, &
         ozwrtvar, ozwrtwf, ozwrtsnr, ozwrtavgk,  ozwrtfavgk, &
         wrtring, wrtozcrs, wrtalbspc, ozwrtint, ozwrtint_fname, &
         lcurve_write, lcurve_fname,lcurve_gcv, ptr_order, ptr_w0, ptr_w1, ptr_w2, &
         biascorr, degcorr, which_biascorr, biasfname, degfname, &
         which_atm, which_clima, which_aperr, which_alb, which_cld,which_toz, norm_tropo3,&
         fnldir, loose_aperr, min_serr, min_terr, smooth_ozbc,  &
         aerosol, strat_aerosol, which_aerosol, scale_aod, scaled_aod, &
         cloud, do_lambcld, lambcld_initalb, scacld_initcod, useasy, nmom, &
         maxmom, ngksec, maxgksec, ngkmatc, maxgkmatc, &
         use_effcrs, hres_samprate, hres_vis_samprate, hres_slitwidth, hres_stok, ex_stok,&
         radc_msegsr, radc_nsegsr, radc_samprate, radc_lambnd, &
         nlay,  nlay_fit, ndiv, nt_fit, t_fidx, t_lidx, tf_fidx, tf_lidx, &
         ozprof_start_index, ozprof_end_index, ozabs_fname,& 
         ozfit_start_index, ozfit_end_index, start_layer, end_layer, &
         use_reg_presgrid, presgrid_fname, use_tropopause, &
         adjust_trop_layer, fixed_ptrop, pst0, ntp0, define_2km_layer, &
         use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
         ngas, nfgas, gasidxs, fgasidxs,fgassidxs, fgaspos, &
         pos_alb, toms_fwhm, ozcrs_alb_fname, vary_sfcalb, &
         nalb, nfalb, albmax, albmin, albidx, albfidx, thealbidx,&
         nwfc, nfwfc, wfcmax, wfcmin, wfcidx, wfcfidx, thewfcidx,&
         ncldaer,ecfrind, ecodind,   ectpind, taodind, twaeind, saodind, sprsind,so2zind,&
         ecfrfind, ecodfind, ectpfind,taodfind, twaefind, saodfind, sprsfind, so2zfind, &
         nothgrp, which_inr, nos, nsl, nsh, nrn, ndc, nis, nir, np1, np2, np3,ncm, &
         osind, osfind, slind, slfind, shind, shfind, rnind, rnfind, dcind, dcfind,&
         isind, isfind, irind, irfind, oswins, slwins, shwins, rnwins, dcwins, &
         iswins, irwins, p1ind, p1find, p1wins, p2ind, p2find, p2wins, &
         p3ind, p3find, p3wins, cmwins, cmind, cmfind, &
         polcorr, the_str,do_rtm_pca, do_brdf, & 
         do_alb_longwav, alb_swav,alb_ewav, use_prefitalb, use_albspc,which_albspc, &
         use_albeofs, nalbspc, nalbspcwin, nalbspcord, is_albspcvar
   USE  OMSAO_parameters_module, ONLY: maxlay,  maxchlen, maxwin, &
         ozstr, othgasstr, albstr,  othstr, cldaerstr, wfcstr, cntrstr
   USE OMSAO_indices_module, ONLY: max_rs_idx,  max_calfit_idx, mxs_idx, &
         maxalb, maxoth, maxgrp, shift_offset, maxcldaer, maxwfc, so2v_idx, &
         hwe_idx, asy_idx, spk_idx, & !vgr_idx, vgl_idx, hwr_idx, hwl_idx
         instrument_idx, gome_idx, gome2_idx
   USE OMSAO_variables_module, ONLY: &
         n_fitvar_rad, mask_fitvar_rad, rmask_fitvar_rad,fitvar_rad_str, &
         fitvar_rad_init, fitvar_rad_init_saved, fitvar_rad_saved, &
         fitvar_rad_unit, lo_radbnd, up_radbnd, numwin, nviswin,winlim, &
         rad_identifier, outdir,tabdir, refdbdir,  do_ch2reso,&
         fothvarpos,band_selectors, which_slit, yn_varyslit, &
         npsl, psl_fpos, scnwrt, instrument_sidx
  USE OMSAO_errstat_module
  USE m_utilities, only: skip_to_filemark
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,           INTENT (IN) :: fit_ctrl_unit
  CHARACTER (LEN=*), INTENT (IN) :: fit_ctrl_file

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER,           INTENT (OUT):: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                  :: i, j, k, iw, fidx, lidx, file_read_stat, idx, &
       swin, ewin, ntemp, nord, ntotp, thewin, theord, np, fstlay, lstlay, &
       fstlayT, lstlayT, numwin0, l
  INTEGER, DIMENSION (maxoth, 2)      :: tmpwins
  INTEGER, DIMENSION (maxwin, maxoth) :: tmpind, tmpfind
  CHARACTER (LEN=maxchlen)            :: tmpchar
  CHARACTER (LEN=6)                   :: idxchar1
  CHARACTER (LEN=1)                   :: winstr, albspc_typestr
  CHARACTER (LEN=2)                   :: albordstr
  REAL      (KIND=dp)                 :: vartmp, lotmp, uptmp, mnoz, minoz, maxoz, &
                                         mnT, minT, maxT
  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'read_ozprof_input'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER :: errstat

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  np1 = 0  ; np2 = 0; np3=0 ; ncm = 0
  ! -------------------------
  ! Open fitting control file
  ! -------------------------
  OPEN ( UNIT=fit_ctrl_unit, FILE=TRIM(ADJUSTL(fit_ctrl_file)), &
       STATUS='OLD', IOSTAT=errstat)
  IF ( errstat /= pge_errstat_ok ) THEN
     WRITE(www_lun, *) 'errors in reading ', fit_ctrl_file
     errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  WRITE(www_lun, *) 'O3CTR:'//TRIM(ADJUSTL(fit_ctrl_file)) 
  ! ------------------------------------------------
  ! Read ozone profile retrieval control variables
  ! ------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, cntrstr, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  WRITE(www_lun, *) 'O3CTR:'//TRIM(ADJUSTL(cntrstr)) 
  READ (fit_ctrl_unit, '(A)') algorithm_name
  READ (fit_ctrl_unit, '(A)') algorithm_version
  algorithm_name = TRIM(ADJUSTL(algorithm_name))
  algorithm_version = TRIM(ADJUSTL(algorithm_version))
 
 ! -------------------------------
  ! Read general control variables
  ! -------------------------------
  READ (fit_ctrl_unit, *) do_multi_vza
  IF (do_multi_vza) THEN
     WRITE(www_lun, *) 'Only effective viewing geometry is computed for a ground pixel!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (do_ch2reso) do_multi_vza = .FALSE.

  READ (fit_ctrl_unit, *) do_radinter
  READ (fit_ctrl_unit, *) do_simu, do_simu_rmring
  IF ( .NOT. do_simu ) do_simu_rmring = .FALSE.
  
  READ (fit_ctrl_unit, *) use_oe
  READ (fit_ctrl_unit, *) do_twostep, do_bothstep, use_large_so2_aperr
  READ (fit_ctrl_unit, *) do_tracewf 
  READ (fit_ctrl_unit, *) 
  READ (fit_ctrl_unit, *) atmwrt, scnwrt, ozwrtint, gaswrt, ozwrtvar, ozwrtcorr, &
       ozwrtcovar, ozwrtavgk, ozwrtfavgk, ozwrtcontri, ozwrtres, ozwrtwf, ozwrtsnr, &
       wrtring, wrtozcrs, wrtalbspc
  ozwrtint_fname = TRIM(ADJUSTL(outdir)) // 'inter_' // rad_identifier // '.dat'
  READ (fit_ctrl_unit, *) 
  READ (fit_ctrl_unit, '(A)') atmos_prof_fname
  atmos_prof_fname = TRIM(ADJUSTL(tabdir)) // TRIM(ADJUSTL(atmos_prof_fname))

  !  Calibration options
  READ (fit_ctrl_unit, *) 
  READ (fit_ctrl_unit, *) radcalwrt, which_caloz
  READ (fit_ctrl_unit, '(A)') caloz_fname
  caloz_fname = TRIM(ADJUSTL(tabdir)) // TRIM(ADJUSTL(caloz_fname))

  IF (radcalwrt .AND. .NOT. do_simu) ozwrtres = .FALSE.
  IF (radcalwrt .AND. do_simu) THEN
     ozwrtavgk = .FALSE.; ozwrtfavgk = .FALSE. ; ozwrtint = .FALSE.
     ozwrtcorr = .FALSE.; ozwrtcovar = .FALSE.;  ozwrtcontri = .FALSE.
     ozwrtwf   = .FALSE.; ozwrtsnr   = .FALSE.;  ozwrtres = .TRUE.
     ozwrtvar  = .FALSE.; wrtalbspc  = .FALSE.
  ENDIF

  READ (fit_ctrl_unit, *) use_lograd
  READ (fit_ctrl_unit, *) use_logstate
  READ (fit_ctrl_unit, *) use_flns
  READ (fit_ctrl_unit, *) ring_on_line, ring_convol, fit_atanring, ring_LUT
  READ (fit_ctrl_unit, *) biascorr, degcorr  
  READ (fit_ctrl_unit, *) which_biascorr
  READ (fit_ctrl_unit, '(A)') biasfname
  biasfname = TRIM(ADJUSTL(tabdir)) // TRIM(ADJUSTL(biasfname))
  READ (fit_ctrl_unit, '(A)') degfname
  degfname = TRIM(ADJUSTL(tabdir)) // TRIM(ADJUSTL(degfname))

  READ (fit_ctrl_unit, *) do_subfit
  READ (fit_ctrl_unit, *) lcurve_gcv
  READ (fit_ctrl_unit, *) lcurve_write
  lcurve_fname = TRIM(ADJUSTL(outdir)) // 'lcurve_' // rad_identifier // '.dat'
  READ (fit_ctrl_unit, *) ptr_order, ptr_w0, ptr_w1, ptr_w2
  READ (fit_ctrl_unit, *) which_atm
  READ (fit_ctrl_unit, *) which_clima
  READ (fit_ctrl_unit, *) which_aperr
  READ (fit_ctrl_unit, *) loose_aperr, min_serr, min_terr
  READ (fit_ctrl_unit, *) smooth_ozbc
  READ (fit_ctrl_unit, *) which_toz
  READ (fit_ctrl_unit, *) norm_tropo3
  READ (fit_ctrl_unit, *) which_alb
  READ (fit_ctrl_unit, *) which_cld
  
  fnldir='fnl13.75LST'
  IF (instrument_idx == gome_idx .or. instrument_idx == gome2_idx) THEN 
    fnldir = 'fnl9.5LST'
  ENDIF

  IF (which_atm > 3 .OR. which_atm < 0) THEN
      WRITE(www_lun, *) modulename, ' No such meterology '
      pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (which_clima > 13 .OR. which_clima <= 0) THEN
     WRITE(www_lun, *) modulename, ' No such ozone profile climatology!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (which_aperr > 13 .OR. which_aperr <= 0 .OR.  &
      which_aperr == 4 .or. which_aperr == 5) THEN
     WRITE(www_lun, *) modulename, ' No such ozone profile a priori error!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (which_toz > 3 .OR. which_toz < 0 ) THEN
     WRITE(www_lun, *) modulename, ' No such total ozone field!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF ( which_toz == 0) THEN 
     IF (which_clima == 1 ) THEN
        WRITE(www_lun, *) modulename, ' Total ozone is needed to use this climatology!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ENDIF
  IF (which_alb > 7) THEN
     WRITE(www_lun, *) modulename, ' No such albedo database!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (which_cld > 5) THEN
     WRITE(www_lun, *) modulename, ' No such cloud option!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (which_cld == 5 .OR. which_toz == 4) THEN
     !CALL read_cldalb_input(errstat)
     STOP
     IF (errstat /= pge_errstat_ok) THEN
        WRITE(*, *) modulename, ' Error in reading cldalb fit input!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ENDIF

  READ (fit_ctrl_unit, *) aerosol, strat_aerosol
  IF (.NOT. aerosol) strat_aerosol = .FALSE.
  READ (fit_ctrl_unit, *)
  READ (fit_ctrl_unit, *) which_aerosol, scale_aod, scaled_aod
  READ (fit_ctrl_unit, *) cloud, do_lambcld, lambcld_initalb, scacld_initcod
  IF (.NOT. cloud) do_lambcld = .FALSE.
  !IF (.NOT. do_lambcld) THEN
  !   WRITE(www_lun, *) modulename, ': Clouds must be assummed Lambertian!!!'
  !   pge_error_status = pge_errstat_error; RETURN
  !ENDIF
  READ (fit_ctrl_unit, *) useasy
  READ (fit_ctrl_unit, *) nmom
  IF (nmom > maxmom) THEN
     WRITE(www_lun, *) modulename, ': Input # of phase moments exceeds maxmom: ', maxmom
     WRITE(www_lun, *) modulename, ': Reduce it or increase maxmom!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (.NOT. aerosol .AND. (.NOT. cloud .OR. do_lambcld) ) THEN
     nmom = 2
  ENDIF

  READ (fit_ctrl_unit, *) use_effcrs
  IF (.NOT. use_effcrs) do_radinter = .FALSE.
  READ (fit_ctrl_unit, *) hres_samprate, hres_vis_samprate
  READ (fit_ctrl_unit, *) hres_slitwidth
  READ (fit_ctrl_unit, *) hres_stok, ex_stok
  IF (hres_samprate < 0.0005) THEN
     WRITE(www_lun, *) modulename, ':hres_samprate must be > 0.0005!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  
  IF ( MOD(hres_samprate * 100, 1.0) /= 0) THEN
     WRITE(www_lun, *) modulename, ': hres_samprate must be multiples of 0.01!!!'
     !hres_samprate = NINT(hres_samprate * 100) / 100.
     !WRITE(www_lun, *) modulename, ': hres_samprate is reset to: ', hres_samprate
     !stop
  ENDIF
  
  ! Read parameters for specifying sampling rate for specified # of spectral region
  ! Only valid if use_effcrs is set to false
  READ (fit_ctrl_unit, *) radc_nsegsr
  IF (radc_nsegsr > radc_msegsr) THEN
     WRITE(www_lun, *) modulename, ': Need to increase radc_msegsr!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  IF (radc_nsegsr <= 0) THEN
     WRITE(www_lun, *) modulename, ': radc_nsegsr must be > 0!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  READ (fit_ctrl_unit, *) (radc_lambnd(i), i = 1, radc_nsegsr+1)
  DO i = 2, radc_nsegsr
     IF (radc_lambnd(i) <= radc_lambnd(i-1)) THEN
        WRITE(www_lun, *) modulename, ': Need to be increasing order!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF
  ENDDO
  IF (radc_lambnd(1) > winlim(1, 1) - 5) THEN
     radc_lambnd(1) = winlim(1, 1) - 5
  ENDIF
  READ (fit_ctrl_unit, *) (radc_samprate(i), i = 1, radc_nsegsr)
  !IF ( MINVAL(radc_samprate(1:radc_nsegsr)) < hres_samprate ) THEN
  !   WRITE(www_lun, *) modulename, ': radc_samprate must be > hres_samprate!!!', radc_samprate(1:radc_nsegsr), hres_samprate
  !   pge_error_status = pge_errstat_error; RETURN
  !ENDIF
  READ (fit_ctrl_unit, *) the_str
  READ (fit_ctrl_unit, *) polcorr 
  READ (fit_ctrl_unit, *) do_rtm_pca
  READ (fit_ctrl_unit, *) do_brdf
  IF (do_rtm_pca) use_effcrs = .false.
  IF (polcorr > 7 ) THEN
     WRITE(www_lun, *) modulename, ': No such polarization correction option!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  IF (polcorr == 0 .OR.polcorr==2 .OR. polcorr >= 3 ) THEN
     useasy = .FALSE.
     ngksec = maxgksec; ngkmatc = maxgkmatc
  ELSE
     ngksec = 1; ngkmatc = 1
  ENDIF

  ! -------------------------------------
  ! Read ozone profile control variables
  ! -------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, ozstr, tmpchar, file_read_stat )
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  WRITE(www_lun, *) 'O3CTR:'//TRIM(ADJUSTL(ozstr)) 
  !xliu (02/08/2007): modify the way of reading input for atmospheric profiles
  !                   add more options for atmospheric layering scheme
  READ (fit_ctrl_unit, '(A)') ozabs_fname
  ozabs_fname = TRIM(ADJUSTL(refdbdir)) // ozabs_fname
  READ (fit_ctrl_unit, *) nlay
  IF (nlay > maxlay) THEN
     WRITE(www_lun, *) modulename, ' : # of layers exceed maximum # of layers allowed!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF 
  READ (fit_ctrl_unit, *) ndiv
  READ (fit_ctrl_unit, *) fstlay, lstlay
  READ (fit_ctrl_unit, *) mnoz, minoz, maxoz
  READ (fit_ctrl_unit, *) fstlayT, lstlayT
  READ (fit_ctrl_unit, *) mnT, minT, maxT
  READ (fit_ctrl_unit, *) use_reg_presgrid
  READ (fit_ctrl_unit, '(A)') presgrid_fname
  presgrid_fname = TRIM(ADJUSTL(tabdir)) // TRIM(ADJUSTL(presgrid_fname))
  READ (fit_ctrl_unit, *) use_tropopause
  READ (fit_ctrl_unit, *) fixed_ptrop, pst0, ntp0
  READ (fit_ctrl_unit, *) adjust_trop_layer  
  READ (fit_ctrl_unit, *) define_2km_layer
  
  ! -----------------------------------------------
  ! Read initial ozone profile variables
  ! -----------------------------------------------
  idx  = shift_offset + max_rs_idx
  ozprof_start_index = idx + 1
  ozprof_end_index = idx + nlay
  t_fidx = ozprof_start_index + maxlay
  t_lidx = ozprof_end_index + maxlay
  fitvar_rad_unit(ozprof_start_index:ozprof_end_index) = 'DU'
  fitvar_rad_unit(t_fidx:t_lidx) = 'K'
 
  fitvar_rad_init(ozprof_start_index:ozprof_end_index) = mnoz
  lo_radbnd(ozprof_start_index:ozprof_end_index)       = minoz
  up_radbnd(ozprof_start_index:ozprof_end_index)       = maxoz
  DO i = 1, fstlay-1
     lo_radbnd(idx+i) = mnoz
     up_radbnd(idx+i) = mnoz
  ENDDO
  DO i = lstlay+1, nlay
     lo_radbnd(idx+i) = mnoz
     up_radbnd(idx+i) = mnoz
  ENDDO

  fitvar_rad_init(t_fidx:t_lidx) = mnT
  lo_radbnd(t_fidx:t_lidx)       = minT
  up_radbnd(t_fidx:t_lidx)       = maxT
  DO i = 1, fstlayT-1
     lo_radbnd(t_fidx+i-1) = mnT
     up_radbnd(t_fidx+i-1) = mnT
  ENDDO
  DO i = lstlayT+1, nlay
     lo_radbnd(t_fidx+i-1) = mnT
     up_radbnd(t_fidx+i-1) = mnT
  ENDDO
  DO i = 1, nlay
     WRITE(fitvar_rad_str(idx+i), '(A2, I2.2)') 'oz', i
  ENDDO
  DO i = fstlayT, lstlayT
     WRITE(fitvar_rad_str(t_fidx+i-1), '(A2, I2.2)') 'tt', i
  ENDDO
!  ozpfpars: DO i = 1, nlay
!     READ (fit_ctrl_unit, *, IOSTAT=errstat) vartmp, lotmp, uptmp, vart
!     IF ( errstat /= pge_errstat_ok ) THEN
!        errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
!             TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
!        WRITE(www_lun, *) modulename, ' : Error in reading initial ozone variables!!!'
!        pge_error_status = pge_errstat_error; RETURN
!     END IF
!     ! ---------------------------------------------------------
!     ! Check for consitency of bounds and adjust where necessary
!     ! ---------------------------------------------------------
!     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
!        lotmp = vartmp ; uptmp = vartmp
!     END IF
!     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
!        uptmp = vartmp ; lotmp = vartmp
!     END IF
!     
!     fitvar_rad_init(idx + i) = vartmp
!     WRITE(fitvar_rad_str(idx+i), '(A2, I2.2)') 'oz', i
!     
!     lo_radbnd(idx + i) = lotmp
!     up_radbnd(idx + i) = uptmp  
!     
!     ! determine whether to vary temperature
!     IF (vart .AND. lotmp < uptmp) THEN
!        fitvar_rad_init(idx + i + maxlay) = 200.0  ! rechange after preparing atmos.
!        lo_radbnd(idx + i + maxlay) = 0.
!        up_radbnd(idx + i + maxlay) = 400.
!        WRITE(fitvar_rad_str(idx+i+maxlay), '(A2, I2)') 'tmp', i
!
!     ELSE
!        fitvar_rad_init(idx + i+ maxlay) = 0. 
!        lo_radbnd(idx + i + maxlay) = 0.
!        up_radbnd(idx + i + maxlay) = 0.
!     END IF
!
!  END DO ozpfpars
  ! -------------------------------------------------------------
  ! Add unfixed ozone profile variables to the variable list 
  ! -------------------------------------------------------------
  ntemp = n_fitvar_rad
  ozfit_start_index = n_fitvar_rad + 1
  DO i = idx + 1, idx + nlay
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad        
     END IF
     IF (n_fitvar_rad == ntemp + 1) THEN
        start_layer = i - idx
     END IF
  END DO
  ozfit_end_index = n_fitvar_rad
  nlay_fit = ozfit_end_index - ozfit_start_index + 1
  end_layer = start_layer + nlay_fit - 1
  ntemp = n_fitvar_rad 
  tf_fidx = 0; tf_lidx = 0; nt_fit = 0
  
  ! Add unfixed temperature variables to the variable list
  DO i = idx + 1 + maxlay, idx + nlay + maxlay
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad  
        IF (tf_fidx == 0) tf_fidx = n_fitvar_rad        
     END IF

  END DO
  IF (tf_fidx > 0) THEN
     tf_lidx = n_fitvar_rad
     nt_fit = tf_lidx - tf_fidx + 1
  ENDIF
  
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Read the control variable for T-dependent cross gases
   !------------------------------------------------
   REWIND ( fit_ctrl_unit )
    CALL skip_to_filemark ( fit_ctrl_unit, othgasstr, tmpchar, file_read_stat )
     IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
           TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
      pge_error_status = pge_errstat_error; RETURN
   END IF
   WRITE(www_lun,*) 'O3CTR'//ADJUSTL(TRIM(othgasstr))
   READ (fit_ctrl_unit, *) use_so2dtcrs
   READ (fit_ctrl_unit, *) use_o4dtcrs
   READ (fit_ctrl_unit, *) use_o2dptcrs
   READ (fit_ctrl_unit, *) use_h2odptcrs

  ! -----------------------------------------------
  ! Read wavelength dependent surface albedo terms
  ! -----------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, albstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  WRITE(www_lun,*) 'O3CTR'//ADJUSTL(TRIM(albstr))
  READ (fit_ctrl_unit, *) pos_alb, toms_fwhm
  READ (fit_ctrl_unit, '(A)') ozcrs_alb_fname
  READ (fit_ctrl_unit, *) vary_sfcalb
  READ (fit_ctrl_unit, *) do_alb_longwav, alb_swav, alb_ewav
  READ (fit_ctrl_unit, *) use_prefitalb
  READ (fit_ctrl_unit, *) use_albspc, use_albeofs, which_albspc

  IF (use_albspc .AND. nviswin == 0) THEN
    use_albspc = .FALSE.
    WRITE(www_lun,*) ' turn off use_albspc due to no visible fitting'
  ENDIF

  IF (.NOT. use_albspc) THEN 
    use_albeofs = .FALSE.
  ENDIF
  IF (use_albspc) vary_sfcalb = .TRUE.
  READ (fit_ctrl_unit, *) nalbspc
  READ (fit_ctrl_unit, *) nalbspcwin
  READ (fit_ctrl_unit, *) nalbspcord

  IF ( use_albspc ) THEN
    IF (nalbspcwin /= 1 .AND. nalbspcwin /= nviswin) THEN
      WRITE(*, *) modulename, ' : # windows for albedo spectrum must be 1 or # fitting windows!!!'
      pge_error_status = pge_errstat_error; RETURN
    ENDIF

    IF (use_albeofs) THEN
      IF (nalbspcord /= 1) nalbspcord = 1
       IF (which_albspc > 2) THEN 
         WRITE(*, *) modulename, ' : check which_albspc, it should be 1 or 2!!!'
         pge_error_status = pge_errstat_error; RETURN
       ENDIF
    ELSE
       IF (nalbspcord < 1) THEN
         WRITE(*, *) modulename, ' : Use >= 1 parameter for albedo spectrum in each window !!!'
         pge_error_status = pge_errstat_error; RETURN
       ENDIF
       nalbspc = 1
    ENDIF
  ELSE
    nalbspcwin = 0
    nalbspcord = 0
    nalbspc = 0
  ENDIF

  IF (do_alb_longwav .and. nviswin /= 0) THEN 
    WRITE(www_lun, *) TRIM(adjustl(modulename))//':longwav alb retrievals not completed!' ; stop
  ENDIF
  ! check the upper and lower range for albedo and cloud 
  ! Use waves < alb_swav and waves > alb_ewav to retrieve albedo
  ! We do some checking here
  ! Make sure that:
  ! (1) the start/end wavelength is within the fitting window.
  ! (2) if alb_swav >= alb_ewav then use waves > alb_ewav only
  ! (3) we have enough points for fitting, which depends on parameters.
  !     for maximum 6 pars, we need a window of at least 5 nm (fwhm=0.45) width
  IF (do_alb_longwav .AND. numwin == 2) THEN
     IF (alb_ewav >= winlim(numwin, 2) .OR. alb_ewav <= winlim(numwin, 1) ) THEN
          WRITE(*, *) modulename, ' : albswav is not within the fitting window.', alb_swav
     ENDIF
     IF (alb_swav >= alb_ewav) alb_swav = 0.0D0
  ENDIF
  numwin0 = numwin

  READ (fit_ctrl_unit, *) nalb
  IF (nalb > maxalb) THEN
     WRITE(www_lun, *) modulename, ' : # of albedo terms exceeds allowed  ', maxalb
     pge_error_status = pge_errstat_error; RETURN
  ELSE IF (nalb < 1) THEN
     WRITE(www_lun, *) modulename, ' : Need to specify at least 1 albedo term!!!'
     pge_error_status = pge_errstat_error; RETURN     
  ENDIF

  idx = shift_offset + max_rs_idx + maxlay * 2
  albidx = idx + 1
  albmin= 0.0; albmax = 0.0
  albpars: DO i = 1, nalb
     READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, &
          lotmp, uptmp, albmin(i), albmax(i)
     is_albspcvar(i) = .FALSE.
     IF ( errstat /= pge_errstat_ok ) THEN
       errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
       TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
       WRITE(www_lun, *) modulename, ' : Error in reading initial albedo variables!!!'
       pge_error_status = pge_errstat_error; RETURN
     END IF

     ! Check for spectral, must cover entire 1 or more fitting windows
     DO j = 1, numwin
       IF ((albmin(i) > winlim(j, 1) .AND. albmin(i) < winlim(j, 2)) .OR. &
           (albmax(i) > winlim(j, 1) .AND. albmax(i) < winlim(j, 2)) ) THEN
          WRITE(*, *) modulename, ' : This albedo only covers partial fitting window: ', TRIM(ADJUSTL(idxchar1))
          !pge_error_status = pge_errstat_error; RETURN
        ENDIF
     ENDDO

     IF (use_albspc) THEN
       IF (albmin(i) <= winlim(numwin-nviswin+1, 1) .AND. albmax(i) >=winlim(numwin-nviswin+1, 2)) THEN
         lotmp = uptmp  ! Disable albedo variables for this and afterward
         nalb = i - 1; EXIT
       ENDIF
     ENDIF
     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
          lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
          uptmp = vartmp ; lotmp = vartmp
     END IF
    
     fitvar_rad_init(idx+i) = vartmp
     fitvar_rad_str (idx+i) = TRIM(ADJUSTL(idxchar1))
     lo_radbnd(idx + i) = lotmp
     up_radbnd(idx + i) = uptmp
  END DO albpars
  
  ! Add albedo spectrum variables
  IF (use_albeofs) THEN
     albspc_typestr = 'e'
  ELSE
     albspc_typestr = 's'
  ENDIF

  ! Add albedo spectrum variables, index starting from i as
  ! there are i - 1 albedo vairables previously
  DO j = 1, nalbspcwin
    IF (nalbspcwin == nviswin ) THEN
      WRITE(winstr, '(I1)') band_selectors(numwin - nviswin + j )
    ELSE
      winstr = 'v'
    ENDIF
    DO k = 1, nalbspc
      DO l = 1, nalbspcord
        fitvar_rad_init(idx + 1) = 0.0d0
        IF (l == 1 .AND. .NOT. use_albeofs) THEN
          fitvar_rad_init(idx + i) = 1.0d0
        ENDIF
        IF (use_albeofs) THEN
          WRITE(albordstr, '(I2)') k
        ELSE
          WRITE(albordstr, '(I2)') l - 1
        ENDIF
        fitvar_rad_str(idx + i) = winstr // 'a' // albspc_typestr//ADJUSTL(TRIM(albordstr))
        IF (.NOT. use_albeofs) THEN
          IF ( l == 1 ) THEN
            lo_radbnd(idx + i) = 0.0d0
            up_radbnd(idx + i) = 5.0d0
          ELSE
            lo_radbnd(idx + i) = -1.0d0
            up_radbnd(idx + i) =  1.0d0
          ENDIF
        ELSE
            ! Need a differnet set for chris and peter ?
            lo_radbnd(idx + i) = -5.0d0 ! changed from = -10
            up_radbnd(idx + i) =  5.0d0 ! changed from = 10
        ENDIF

        is_albspcvar(i) = .TRUE.
        IF (winstr == 'v') THEN
          albmin(i) = winlim(numwin - nviswin + 1, 1)-2 ! more space need
          albmax(i) = winlim(numwin, 2)+2
        ELSE
          albmin(i) = winlim(numwin - nviswin + j, 1)-2
          albmax(i) = winlim(numwin - nviswin + j, 2)+2
        ENDIF
           
        i = i + 1
      ENDDO
   ENDDO
  ENDDO
  nalb = i - 1
  ! -------------------------------------------------------------
  ! Add albedo variable to the variable list
  ! -------------------------------------------------------------
  albfidx = 0; thealbidx = 0
  DO i = idx + 1, idx + nalb
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad  
        IF (albfidx == 0) albfidx = n_fitvar_rad
        IF ( fitvar_rad_str(i)(4:4)== '0') thealbidx = n_fitvar_rad
     ENDIF
  ENDDO
  IF (albfidx > 0) THEN
     nfalb = n_fitvar_rad - albfidx + 1
     thealbidx = thealbidx - albfidx + 1
  ENDIF
  ! ------------------------------------------------
  ! Read wavelength dependent cloud fraction terms
  ! ------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, wfcstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF

  READ (fit_ctrl_unit, *) nwfc
  IF (nwfc > maxwfc) THEN
     WRITE(www_lun, *) modulename, ' : # of wavelength dependent terms exceed allowed  ', maxwfc
     pge_error_status = pge_errstat_error; RETURN
  ELSE IF (nwfc < 1) THEN
     nwfc = 0
  ENDIF

  idx = idx + maxalb
  wfcidx = idx + 1
  wfcmin= 0.0; wfcmax = 0.0
  wfcpars: DO i = 1, nwfc
     READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, &
          lotmp, uptmp, wfcmin(i), wfcmax(i)
     IF ( errstat /= pge_errstat_ok ) THEN
       errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
       TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
        WRITE(www_lun, *) modulename, ' : Error in reading initial cloud fraction variables!!!'
        pge_error_status = pge_errstat_error; RETURN
     END IF

     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
          lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
          uptmp = vartmp ; lotmp = vartmp
     END IF
    
     fitvar_rad_init(idx + i) = vartmp
     fitvar_rad_str(idx+i) = TRIM(ADJUSTL(idxchar1))
     lo_radbnd(idx + i) = lotmp
     up_radbnd(idx + i) = uptmp
  END DO wfcpars

  ! ---------------------------------------------------------------------
  ! Add wavelength-dependent cloud fraction variables to the variable list
  ! ----------------------------------------------------------------------
  wfcfidx = 0; thewfcidx=0
  DO i = idx + 1, idx + nwfc
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad  
        IF (wfcfidx == 0) wfcfidx = n_fitvar_rad
        IF ( fitvar_rad_str(i)(4:4)== '0') thewfcidx = n_fitvar_rad
     END IF
  END DO
  IF (wfcfidx > 0) THEN
     nfwfc = n_fitvar_rad - wfcfidx + 1
     thewfcidx = thewfcidx - albfidx + 1
  ENDIF
 
  ! -------------------------------
  ! Read cloud/aerosol variables
  ! -------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, cldaerstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  WRITE(www_lun,*) 'O3CTR'//ADJUSTL(TRIM(cldaerstr))
  READ (fit_ctrl_unit, *) ncldaer
  IF (ncldaer > maxcldaer) THEN
     WRITE(www_lun, *) modulename, ' : Increase maxcldaer in OMSAO_indices...!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF

  idx = wfcidx + maxwfc - 1 
  DO i = 1, ncldaer
     READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, lotmp, uptmp
     IF ( errstat /= pge_errstat_ok ) THEN
       errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
                TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
       WRITE(www_lun, *) modulename, ' : Error in reading initial albedo variables!!!'
       pge_error_status = pge_errstat_error; RETURN
     END IF

     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
       lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
       uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( .NOT. cloud .AND. i > 0 .AND. i < 4) THEN   ! Cannot fit any cloud variables
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     ! Disalble this cloud fraction if wavelength dependent cloud fraction is selected
     IF (nwfc > 0 .AND. i == 1) THEN  
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF

     IF ( i == 8 ) THEN ! For SO2Z, SO2V needs to be selected
        j = max_calfit_idx + (so2v_idx - 1) * mxs_idx + 2  ! Index for SO2V
        IF (rmask_fitvar_rad(j) == 0) THEN
           vartmp = 0.0; lotmp = 0.0; uptmp = 0.0          ! Disable SO2Z if SO2V is not selected
        ELSE                                               
           IF (vartmp == 0.0 .AND. lotmp >= uptmp) THEN 
              vartmp = 5.0; lotmp = 5.0; uptmp = 5.0       ! Initialize to 5 km if not initialized or set to zero
           ENDIF
        ENDIF
     ENDIF

     IF ( do_lambcld .AND. i == 2 ) THEN              ! Cannot fit cloud optical thickness
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     IF ( .NOT. aerosol .AND. i > 3 .AND. i < 7) THEN ! Cannot fit aerosol variables
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     IF ( .NOT. strat_aerosol .AND. i == 6) THEN      ! Cannot fit stratospheric aerosols
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
     ! Disable fitting cloud top pressure (not implemented)
     ! Cloud optical thickness is directly fitted either from longer wavelengths, 
     ! which is enabled when using scattering clouds, or a cloud fraction is directly fitted 
     ! and when cloud fraction is greater than 1 (need to exchange fitting variable between cloud 
     ! fraction and cloud optical thickness)
     IF ( i == 2 .OR. i == 3) THEN                   
        vartmp = 0.0; lotmp = 0.0; uptmp = 0.0
     ENDIF
    
     fitvar_rad_init(idx + i) = vartmp
     fitvar_rad_str(idx+i) = TRIM(ADJUSTL(idxchar1))
     lo_radbnd(idx + i) = lotmp
     up_radbnd(idx + i) = uptmp          
  ENDDO
  ecfrind = idx + 1; ecodind = idx + 2; ectpind = idx + 3
  taodind = idx + 4; twaeind = idx + 5; saodind = idx + 6
  sprsind = idx + 7; so2zind = idx + 8
  
  ecfrfind = 0; ecodfind = 0; ectpfind = 0
  taodfind = 0; twaefind = 0; saodfind = 0
  sprsfind = 0; so2zfind = 0
 
  ! -------------------------------------------------------------
  ! Add Cloud/Aerosol variables to the variable list
  ! -------------------------------------------------------------
  DO i = idx + 1, idx + ncldaer
     j = i -idx
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        rmask_fitvar_rad(i) = n_fitvar_rad  
        IF      (j == 1 ) THEN
           ecfrfind = n_fitvar_rad
        ELSE IF (j == 2) THEN
           ecodfind = n_fitvar_rad
        ELSE IF (j == 3) THEN
           ectpfind = n_fitvar_rad
           fitvar_rad_unit(i) = 'mb'
        ELSE IF (j == 4) THEN
           taodfind = n_fitvar_rad
        ELSE IF (j == 5) THEN
           twaefind = n_fitvar_rad
        ELSE IF (j == 6) THEN
           saodfind = n_fitvar_rad
        ELSE IF (j == 7) THEN
           sprsfind = n_fitvar_rad
           fitvar_rad_unit(i) = 'mb'
        ELSE IF (j == 8) THEN
           so2zfind = n_fitvar_rad
           fitvar_rad_unit(i) = 'km'
        ENDIF
     END IF
  END DO
  
  ! Read other parameters (for multiple window)
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, othstr, tmpchar, file_read_stat)
  IF ( file_read_stat /= file_read_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF    
  WRITE(www_lun,*) 'O3CTR'//ADJUSTL(TRIM(othstr))
  READ (fit_ctrl_unit, *) which_inr
  READ (fit_ctrl_unit, *) nothgrp 
  IF (nothgrp > maxgrp) THEN
     WRITE(www_lun, *) modulename, ' : Increase maxgrp in OMSAO_indices...!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
   
  k = wfcidx + maxwfc + maxcldaer - 1 
  IF (do_subfit) THEN
     ntotp = numwin * maxoth; np = numwin
  ELSE
     ntotp = maxoth; np = 1
  ENDIF

  DO i = 1, nothgrp        ! for each group of parameters
     nord = 0; tmpind = 0; tmpfind = 0; tmpwins = 0
     DO j = 1, maxoth   ! for each order of parameters
       
        READ (fit_ctrl_unit, *, IOSTAT=errstat) idxchar1, vartmp, lotmp, uptmp, swin, ewin
        IF ( errstat /= pge_errstat_ok ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_read_fitctrl_file, &
                TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
           WRITE(www_lun, *) modulename, ' : Error in reading other variables!!!'
           pge_error_status = pge_errstat_error; RETURN
        END IF
        ! special conditions for fit_atanring
        IF (i == 4 .AND. fit_atanring) THEN
           swin = 1;   ewin = 1
           IF (j > 3) THEN
              lotmp = 0.0; uptmp = 0.0
           ELSE
              lotmp = -1.0D99; uptmp = 1.0D99
           ENDIF
           IF (j == 1) THEN
              vartmp = 0.5
           ELSE IF (j == 2) THEN
              vartmp = 300.
           ELSE IF (j == 3) THEN
              vartmp = 5.0
           ELSE
              vartmp = 0.0
           ENDIF
        ENDIF
        
        IF (swin > ewin) THEN
           ntemp = swin; swin = ewin; ewin = ntemp
        ENDIF
        IF (swin < 1 .AND. ewin < 1) THEN
           swin = 1; ewin = numwin
        ENDIF
        IF (swin < 1) swin = 1
        IF (ewin > numwin) ewin = numwin       
        tmpwins(j, 1) = swin; tmpwins(j, 2) = ewin
              
        ! ---------------------------------------------------------
        ! Check for consitency of bounds and adjust where necessary
        ! ---------------------------------------------------------
        IF (lotmp > vartmp .OR. uptmp < vartmp ) THEN
           lotmp = vartmp ; uptmp = vartmp
        END IF
        IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
           uptmp = vartmp ; lotmp = vartmp
        END IF
        
        IF (do_subfit) THEN
           lo_radbnd (k + 1 : k + numwin) = 0.0
           up_radbnd (k + 1 : k + numwin) = 0.0
           fitvar_rad_init (k + 1 : k + numwin) = 0.0
           IF (i == 1 .OR. i == 3) fitvar_rad_unit(k+1:k+numwin) = 'nm'
                      
           fitvar_rad_init(k + swin : k + ewin) = vartmp
           fitvar_rad_str (k + swin : k + ewin) = TRIM(ADJUSTL(idxchar1))
           lo_radbnd      (k + swin : k + ewin) = lotmp
           up_radbnd      (k + swin : k + ewin) = uptmp
           IF (i == 1 .OR. i == 3) fitvar_rad_unit(k+swin:k+ewin) = 'nm'
           tmpind(1:numwin, j) = k + (/(idx, idx = 1, numwin)/)
           k = k + numwin
        ELSE
           fitvar_rad_init(k + 1) = vartmp
           fitvar_rad_str (k + 1) = TRIM(ADJUSTL(idxchar1))
           IF (i == 1 .OR. i == 3) fitvar_rad_unit(k + 1) = 'nm'
           lo_radbnd(k + 1) = lotmp
           up_radbnd(k + 1) = uptmp
           k = k + 1; tmpind(1, j) = k

           ! The windows must be the same for all the orders
           IF (j > 1) THEN
              tmpwins(j, :) = tmpwins(1, :)
           ENDIF
        ENDIF

     END DO

     ! If lower order, then no higher orders
     DO j = 1, maxoth-1
        DO iw = 1, np
           fidx = tmpind(iw, j)
           IF (lo_radbnd(fidx) == up_radbnd(fidx)) THEN
              lo_radbnd(tmpind(iw, j+1:maxoth)) = 0.0
              up_radbnd(tmpind(iw, j+1:maxoth)) = 0.0
              fitvar_rad_init(tmpind(iw, j+1:maxoth)) = 0.0
           ENDIF
        ENDDO
     ENDDO

     ntemp = n_fitvar_rad; idx = tmpind(1, 1) - 1
     DO j = 1, ntotp
        IF (lo_radbnd(idx + j) < up_radbnd(idx + j)) THEN
           n_fitvar_rad = n_fitvar_rad + 1
           mask_fitvar_rad(n_fitvar_rad) = idx + j
           rmask_fitvar_rad(idx + j) = n_fitvar_rad  
           theord = CEILING(1.0 * j / np); thewin = j - (theord - 1) * np! a bug with using three windows for cloud fitting.
           ! jbak change
           tmpfind(thewin, theord) = n_fitvar_rad 
        ENDIF
     ENDDO
     IF (n_fitvar_rad > ntemp) nord =  theord
     IF (i == 1) THEN
        nos = nord; oswins = tmpwins; osind = tmpind; osfind = tmpfind
     ELSE IF ( i == 2) THEN
        nsl = nord; slwins = tmpwins; slind = tmpind; slfind = tmpfind
     ELSE IF (i == 3)  THEN
        nsh = nord; shwins = tmpwins; shind = tmpind; shfind = tmpfind
     ELSE IF (i == 4)  THEN
        nrn = nord; rnwins = tmpwins; rnind = tmpind; rnfind = tmpfind
     ELSE IF (i == 5 ) THEN 
        ndc = nord; dcwins = tmpwins; dcind = tmpind; dcfind = tmpfind
     ELSE IF (i == 6)  THEN
        nis = nord; iswins = tmpwins; isind = tmpind; isfind = tmpfind
     ELSE IF (i == 7) THEN
        nir = nord; irwins = tmpwins; irind = tmpind; irfind = tmpfind
     ELSE IF (i == 8) THEN 
        np1 = nord; p1wins = tmpwins; p1ind = tmpind; p1find = tmpfind
     ELSE IF (i == 9) THEN 
        np2 = nord; p2wins = tmpwins; p2ind = tmpind; p2find = tmpfind
     ELSE IF (i == 10 ) THEN    
        np3 = nord; p3wins = tmpwins; p3ind = tmpind; p3find = tmpfind
     ELSE IF (i == 11) THEN
        ncm = nord; cmwins = tmpwins; cmind = tmpind; cmfind = tmpfind
     ENDIF  
  ENDDO
 
  IF (np1 + np2 + np3 /= 0 ) THEN  
    IF (which_slit == 0 ) THEN   ! gaussian
        npsl = 1 ; psl_fpos(1:npsl) = [hwe_idx]
        np2 =0 ; np3 = 0
    ELSE IF (which_slit == 1) THEN  ! asym. gaussian
        npsl = 2 ; psl_fpos(1:npsl) = [hwe_idx, asy_idx]
        np2 = np3; p2wins = p3wins; p2ind = p3ind; p2find = p3find
        np3 = 0
    ELSE IF (which_slit == 2) THEN  ! triangle
         npsl = 1 ; psl_fpos(1:npsl) = [hwe_idx]
       print * , 'not yet implemented ' ; stop
    ELSE IF (which_slit == 3) THEN  ! viot
       !npsl = 4 ; psl_fpos(1:npsl) = [vgl_idx, vgr_idx, hwl_idx,hwr_idx ]
       print * , 'not yet implemented ' ; stop
    ELSE IF (which_slit == 4) THEN  ! super gaussian
        npsl = 2 ; psl_fpos(1:npsl) = [hwe_idx, spk_idx]
       !np2 = np3; p2wins = p3wins; p2ind = p3ind; p2find = p3find
       np3 = 0
       IF (np1 ==0 .or. np2 == 0) THEN 
          PRINT * , 'check for pw or pk'
          STOP
       ENDIF
    ELSE IF (which_slit == 5) THEN ! asym super gaussian
        npsl = 3 ; psl_fpos(1:npsl) = [hwe_idx, spk_idx, asy_idx]
       IF (np1 ==0 .or. np2 == 0 .or. np3 == 0) THEN 
          PRINT * , 'check for pw or pk'
          STOP
       ENDIF
    ELSE IF (which_slit == instrument_sidx) THEN  ! omi instrument
         np1 = 0 ; np2 = 0 ; np3 = 0
     ! current OMI instrument slit function is not fitted, but Kang Sun did it
     ! and jbak did it for OMPS
     ! so could be implemented later
    ENDIF
    IF (yn_varyslit) THEN
       np1 = 0 ; np2 = 0 ; np3 = 0
       WRITE(www_lun,*) modulename//&
         ':pseudo slit function is turn off if yn_varyslit'
    ENDIF
  ENDIF

  !WRITE(www_lun, '(7I5)') nos, nsl, nsh, nrn, ndc, nis, nir

  ! get indices for auxiliary variables in the final fitted array
  rtm_treatment(1:max_rs_idx) = .FALSE.
  rtm_treatment(gasidxs) = .TRUE.
  fgasidxs = 0; nfgas = 0
  DO i = 1, ngas
     fidx =  max_calfit_idx + (gasidxs(i) - 1) * mxs_idx + 1; lidx = fidx + 2
     fitvar_rad_unit(fidx:lidx) = 'molecumes cm^-2'
     fgasidxs(i) = MAXVAL(rmask_fitvar_rad(fidx:lidx))
     fgassidxs(i) = rmask_fitvar_rad(shift_offset + gasidxs(i))
     IF (fgasidxs(i) > 0) THEN
        nfgas = nfgas + 1; fgaspos(nfgas) = i
     ENDIF
  ENDDO

  ! Find indices of variables (other than trace gases and ozone)
  j = 1
  DO i = 1, n_fitvar_rad
     IF (nfgas > 0) THEN
        IF (i >= fgasidxs(fgaspos(1)) .AND. i <= fgasidxs(fgaspos(nfgas))) CYCLE
     ENDIF
     IF (i >= ozfit_start_index .AND. i <= ozfit_end_index) CYCLE
     fothvarpos (j) = i; j = j + 1
  ENDDO

  ! -----------------------------------------------
  ! Close fitting control file, report SUCCESS read
  ! -----------------------------------------------
  CLOSE ( UNIT=fit_ctrl_unit )

  fitvar_rad_saved = fitvar_rad_init
  fitvar_rad_init_saved = fitvar_rad_init
  numwin = numwin0
  !---------------------------
  ! debuging part
  !---------------------------
  IF (nfwfc > 0 .and. do_rtm_pca) THEN
      WRITE(*,'()') modulename//':Rethink, wcfrac is not implemented in PCA RTM !!!';STOP
  ENDIF
  IF (use_prefitalb) THEN 
      WRITE(*,'()') modulename//':use_prefitalb is not really implemented!!!';STOP
  ENDIF
  RETURN
END SUBROUTINE read_ozprof_input


SUBROUTINE read_simuret_input (fit_ctrl_unit, fit_ctrl_file, pge_error_status )     

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: rad2deg, deg2rad
  USE ozprof_data_module,        ONLY: simu_aerosol, simu_strataer, &
       simu_which_aerosol, simu_scale_aod, simu_scaled_aod, simu_cfrac, simu_ctp, simu_lambcld, simu_cod, &
       simu_effcrs, simu_ring, simu_polcorr, simu_snrtbl_fname, simu_nret, caloz_fname, ps0,  &
       simu_cloud, simu_lambcldalb, simu_snowice, &
       simu_usealbedo, simu_albedo, has_glint, saa_flag, nsaa_spike, simu_perturbapriori
  USE OMSAO_variables_module,    ONLY: the_year, the_month, the_day, the_utchr, the_utc, the_sza_atm, the_vza_atm, &
       the_aza_atm, the_sca_atm, the_lons, the_lats, the_lon, the_lat, the_surfalt, numwin, simu_hw1es, simu_easyms, &
       simu_deltalams, simu_doppler, simu_dopplershift, simu_dplshi_lamref, simu_shiftirrad, simu_irradshifts, &
       simu_shiftrad, simu_radshifts, simu_addirradsnr, simu_slitcal, simu_wavecal, edgelons, edgelats, nview, nloc
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,           INTENT (IN) :: fit_ctrl_unit
  CHARACTER (LEN=*), INTENT (IN) :: fit_ctrl_file

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER,           INTENT (OUT):: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                             :: i

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=18), PARAMETER :: modulename = 'read_simuret_input'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER :: errstat

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER :: OMI_SMF_setmsg

  ! -------------------------
  ! Open fitting control file
  ! -------------------------
  OPEN ( UNIT=fit_ctrl_unit, FILE=TRIM(ADJUSTL(fit_ctrl_file)), &
       STATUS='OLD', IOSTAT=errstat)

  IF ( errstat /= pge_errstat_ok ) THEN
     errstat = OMI_SMF_setmsg (omsao_e_open_fitctrl_file, &
          TRIM(ADJUSTL(fit_ctrl_file)), modulename, 0)
     pge_error_status = pge_errstat_error; RETURN
  END IF
  DO i = 1, numwin
    READ (fit_ctrl_unit, *) simu_hw1es(i), simu_easyms(i)
  ENDDO
  READ (fit_ctrl_unit, *) (simu_deltalams(i), i=1, numwin)
  READ (fit_ctrl_unit, *) simu_doppler, simu_dopplershift, simu_dplshi_lamref
  READ (fit_ctrl_unit, *) simu_shiftirrad,  simu_irradshifts
  READ (fit_ctrl_unit, *) simu_shiftrad,  simu_radshifts
  READ (fit_ctrl_unit, *) simu_addirradsnr, simu_slitcal, simu_wavecal
  READ (fit_ctrl_unit, *) the_year, the_month, the_day, the_utchr
  WRITE (the_utc, '(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A14)') the_year, '-', &
       the_month, '-', the_day, 'T', the_utchr, ':00:00.000000Z'
  DO i = 1, 5
    READ (fit_ctrl_unit, *) the_lats(i), the_lons(i)
  ENDDO
  the_lat = the_lats(5);  the_lon = the_lons(5)
  nloc = 5
  edgelons(1) = (the_lons(1) + the_lons(2))/2.0
  edgelons(2) = (the_lons(3) + the_lons(4))/2.0
  edgelats(1) = (the_lats(1) + the_lats(4))/2.0
  edgelats(2) = (the_lats(2) + the_lats(3))/2.0
  
  saa_flag = .FALSE.; nsaa_spike = 0

  READ (fit_ctrl_unit, *) the_sza_atm, the_vza_atm, the_aza_atm
  nview = 1

  the_sca_atm = ACOS(COS(the_sza_atm * deg2rad) * COS(the_vza_atm * deg2rad) &
  + SIN(the_sza_atm * deg2rad) * SIN(the_vza_atm * deg2rad) * COS(the_aza_atm * deg2rad))  * rad2deg
  the_sca_atm = 180.0 - the_sca_atm

  READ (fit_ctrl_unit, *) ps0, the_surfalt
  READ (fit_ctrl_unit, '(A)') caloz_fname
  i = INDEX(caloz_fname, ' '); caloz_fname = caloz_fname(1:i-1)
  caloz_fname = TRIM(ADJUSTL(caloz_fname))
  READ (fit_ctrl_unit, *) simu_perturbapriori

  READ (fit_ctrl_unit, *) simu_aerosol, simu_strataer
  IF (.NOT. simu_aerosol) simu_strataer = .FALSE.
  READ (fit_ctrl_unit, *) simu_which_aerosol, simu_scale_aod, simu_scaled_aod
  READ (fit_ctrl_unit, *) simu_cloud, simu_cfrac, simu_ctp, simu_lambcld, simu_lambcldalb, simu_cod
  IF (.NOT. simu_cloud) simu_lambcld = .FALSE.

  READ (fit_ctrl_unit, *) simu_usealbedo, simu_albedo
  has_glint = .FALSE.

  READ (fit_ctrl_unit, *) simu_snowice  
  READ (fit_ctrl_unit, *) simu_effcrs
  READ (fit_ctrl_unit, *) simu_ring
  READ (fit_ctrl_unit, *) simu_polcorr  
  IF (simu_polcorr > 3 .OR. simu_polcorr == 2 .OR. simu_polcorr < 0) THEN
     WRITE(*, *) modulename, ': No such polarization correction option!!!'
     pge_error_status = pge_errstat_error; RETURN
  ENDIF
  READ (fit_ctrl_unit, '(A)') simu_snrtbl_fname
  i = INDEX(simu_snrtbl_fname, ' '); simu_snrtbl_fname = simu_snrtbl_fname(1:i-1)
  simu_snrtbl_fname = TRIM(ADJUSTL( simu_snrtbl_fname))

  READ (fit_ctrl_unit, *) simu_nret

  CLOSE (UNIT = fit_ctrl_unit)

  RETURN
END SUBROUTINE read_simuret_input

END MODULE m_read_ozprof_input
