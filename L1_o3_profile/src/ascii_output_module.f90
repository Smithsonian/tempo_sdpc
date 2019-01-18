!
module ascii_output_module

  public omi_write_intermed, l2_ascii_create, l2_ascii_close
  private !gome_write_intermed

contains

  function filter_exp (x) result (exp_x)
    real (kind=8), dimension(:), intent(in) :: x
    real (kind=8), dimension(size(x)) :: exp_x
    real (kind=8), parameter :: maxexp = log(huge(1.0d0))
    where (x < maxexp)
      exp_x = exp(x)
    elsewhere
      exp_x = huge(1.0d0)
    endwhere
  end function

  ! ***************** Modification History ******************
  ! xiong liu, July 2003
  ! 1. Add an option for writing ozone retrieval results in 
  !    SUBROUTINE write_intermed
  ! *********************************************************

!  Unused?
!
!  SUBROUTINE gome_write_intermed (founit, npix, nsub, fitcol, dfitcol, &
!       rms, amf, amfgeo, sol_zen_eff, nang, sza_atm, vza_atm, ngeo, &
!       lat, lon, exval )
!
!    USE OMSAO_precision_module
!    !USE OMSAO_indices_module, ONLY: so2_idx, no2_t1_idx, hcho_idx, bro_idx
!    USE ozprof_data_module, ONLY: ozprof_flag, ozprof, ozprof_std, &
!         ozprof_ap, ozprof_apstd, eff_alb, eff_alb_init, &
!         nlay, ozdfs, ozinfo, num_iter, covar, contri, avg_kernel, &
!         use_lograd, nalb, atmosprof, ntp, nlay_fit, ozfit_start_index, &
!         ozfit_end_index, start_layer, end_layer, the_ctp, the_cfrac, &
!         the_orig_cfr, the_orig_ctp, the_cld_flg, &
!         the_ai, lambcld_refl, which_cld, which_alb, &
!         ozprof_nstd, strataod, stratsca, tropaod, tropsca, aerwavs, &
!         actawin, ozwrtcorr, ozwrtcovar, ozwrtcontri, ozwrtres, &
!         ozwrtavgk, ozwrtvar, fgasidxs, ngas, nfgas, gaswrt, tracegas, &
!         saa_flag, nsaa_spike, ozwrtfavgk, favg_kernel, nsfc!, the_cod, &
!         !the_orig_cod, ozprof_init, ozprof_start_index, ozprof_end_index, &
!         !maxawin, use_oe, do_lambcld
!    USE OMSAO_variables_module,  ONLY : the_sza_atm, the_vza_atm, &
!         the_aza_atm, fitvar_rad, mask_fitvar_rad, n_fitvar_rad, &
!         fitvar_rad_std, n_rad_wvl, fitspec_rad, fitres_rad, fitwavs, &
!         fitvar_rad_str, fitvar_rad_nstd!, chisq, fitvar_rad_apriori
!    USE OMSAO_gome_data_module,  ONLY : gome_orbc, gome_pixdate, &
!         orbnum, gome_stpix, gome_npix!, gome_angles_wrtn, gome_endpix
!
!    IMPLICIT NONE
!
!    ! ===============
!    ! Input Variables
!    ! ===============
!    INTEGER,        INTENT (IN) :: founit, npix, nsub, nang, ngeo, exval
!    REAL (KIND=dp), INTENT (IN) :: rms, amf, amfgeo, sol_zen_eff
!    REAL (KIND=dp), DIMENSION(3), INTENT (IN)    :: fitcol
!    REAL (KIND=dp), DIMENSION(3, 2), INTENT (IN) :: dfitcol
!    REAL (KIND=dp), DIMENSION (nang),INTENT (IN) :: sza_atm, vza_atm
!    REAL (KIND=dp), DIMENSION (ngeo),INTENT (IN) :: lat, lon
!
!    ! ===============
!    ! Local variables
!    ! ===============
!    INTEGER                                  :: i, j!, id
!    REAL (KIND=dp), DIMENSION (2*ngeo)       :: latlon
!    REAL (KIND=dp)                           :: avgres
!    REAL (KIND=dp), DIMENSION (n_fitvar_rad) :: correl
!    REAL (KIND=dp), DIMENSION (n_rad_wvl)    :: simrad
!
!    DO i = 1, ngeo
!      latlon(2*i-1) = lat(i)
!      latlon(2*i) = lon(i)
!    END DO
!
!    IF (.NOT. ozprof_flag) THEN
!      WRITE (founit, '(I5,I3,1P8E12.4,I8, 2X, 0P16F7.2)') &
!           npix, nsub, fitcol, dfitcol, fitcol/amf, dfitcol/amf, rms, &
!           amf, amfgeo, sol_zen_eff, exval, sza_atm(1:nang), &
!           vza_atm(1:nang), latlon(1:2*ngeo)
!    ELSE     
!      WRITE(founit, '(A12,I5,A5,2I5,1x,A24, 2I5)') 'GOME Pixel: ', orbnum, &
!           gome_orbc, npix, nsub, gome_pixdate, gome_stpix, gome_npix
!      WRITE(founit,'(3F7.2)') the_sza_atm,the_vza_atm,the_aza_atm
!      WRITE(founit,'(10F7.2)') latlon(1:2*ngeo)
!      WRITE(founit, '(A13,I5, A15, I5, A6, L5, I5)') 'Exit Status: ', &
!           exval, '# Iterations: ', num_iter, ' SAA: ', saa_flag, nsaa_spike
!
!      IF (exval > 0) THEN
!        simrad = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
!        IF (use_lograd) THEN
!          fitspec_rad(1:n_rad_wvl) = EXP(fitspec_rad(1:n_rad_wvl))
!          simrad(1:n_rad_wvl) = EXP(simrad(1:n_rad_wvl))
!          fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  simrad
!        END IF
!        avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
!             fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0
!
!        WRITE(founit, '(4(A9, f8.2))') 'rms = ', rms, ' avgres = ', &
!             avgres, ' dfs = ', ozdfs, ' info. = ', ozinfo
!        WRITE(founit, '(A8,2I3,1P28E12.4)') 'Albedo: ', which_alb, nalb, &
!             eff_alb_init(1:nalb), eff_alb(1:nalb)
!        WRITE(founit, '(A8,2I3, 1P5E12.4)') 'Cloud:  ', which_cld, &
!             the_cld_flg, the_cfrac, the_ctp, the_orig_cfr, the_orig_ctp, &
!             lambcld_refl
!        WRITE(founit, '(A8, I3, 1P36E10.3)') 'Aerosol:  ', actawin, the_ai, &
!             aerwavs(1:actawin), tropaod(1:actawin), tropsca(1:actawin), &
!             strataod(1:actawin), stratsca(1:actawin)
!
!        WRITE(founit, '(A31, 2I5)') 'Atmosphere and ozone profiles: ', nlay, &
!             ntp, nsfc
!        WRITE(founit, '(A)') '  #    P(mb)   Z(km)  T(K)     O3-ap O3-apstd   O3     STD     NSTD (DU)'
!        WRITE(founit, '(A2, 1x, F9.3, 2F8.3)') ' 0', atmosprof(1:3, 0)
!        DO i = 1, nlay 
!          WRITE(founit, '(I2, 1X, F9.3, 7F8.3)') i, atmosprof(1:3, i), &
!               ozprof_ap(i), ozprof_apstd(i), ozprof(i), ozprof_std(i), &
!               ozprof_nstd(i)
!        END DO
!        WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Total Ozone: ', &
!             SUM(ozprof_ap(1:nsfc)), fitcol(1), dfitcol(1, 1), dfitcol(1, 2)
!        IF (ntp > 0) THEN
!          WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Stratospheric Ozone: ',&
!               SUM(ozprof_ap(1:ntp)), fitcol(2), dfitcol(2, 1), dfitcol(2, 2)
!          WRITE(founit, '(A24,4x,F8.3,8X, 3F8.3)') ' Tropospheric Ozone: ',&
!               SUM(ozprof_ap(ntp+1:nsfc)), fitcol(3), dfitcol(3, 1), &
!               dfitcol(3, 2)          
!        ENDIF
!
!        IF (gaswrt) THEN
!          WRITE(founit, '(A, I5)') &
!               'Fitted trace gases and uncertainty: ', nfgas
!          WRITE(founit, '(A)')  &
!               '    Var    Initial   A Priori A Priori Std  VCD      STD       NSTD      AMF      ACFRAC    AVGK(1)   AVGK(2)'
!          DO i = 1, ngas
!            IF (fgasidxs(i) > 0) THEN
!              j = mask_fitvar_rad(fgasidxs(i))
!              WRITE(founit, '(I2,1x,A6,1P10d10.2)') fgasidxs(i), &
!                   fitvar_rad_str(j), tracegas(i, 1:10)
!            ENDIF
!          ENDDO
!        ENDIF
!
!        IF (ozwrtvar) THEN                      
!          WRITE(founit, '(A, I5)') &
!               'Fitted variables and uncertainty: ', n_fitvar_rad
!          DO i = 1, n_fitvar_rad
!            WRITE(founit, '(I2, 1x,A6,1P3d10.2)') i, &
!                 fitvar_rad_str(mask_fitvar_rad(i)), &
!                 fitvar_rad(mask_fitvar_rad(i)), &
!                 fitvar_rad_std(mask_fitvar_rad(i)), &
!                 fitvar_rad_nstd(mask_fitvar_rad(i))
!          END DO
!        ENDIF
!
!        IF (ozwrtcorr) THEN
!          WRITE(founit, '(A, I5)') 'Correlation matrix: ', n_fitvar_rad
!          DO i = 1, n_fitvar_rad
!            correl(i) = 1.0
!            DO j = 1, i - 1
!              correl(j) = covar(i, j) / SQRT(covar(i, i) * covar(j, j))
!            END DO
!            WRITE(founit, '(I2,1X,100f6.2)') i, correl(1:i)
!          END DO
!        ENDIF
!
!        IF (ozwrtcovar) THEN
!          WRITE(founit, '(A, 3I5)') 'Covariance matrix: ', nlay_fit, &
!               start_layer, end_layer
!          DO i = ozfit_start_index, ozfit_end_index
!            WRITE(founit, '(1p60d10.2)') &
!                 covar(i, ozfit_start_index:ozfit_end_index)
!          END DO
!        ENDIF
!
!        IF (ozwrtavgk) THEN
!          WRITE(founit, '(A, 3I5)') 'Average kernel: ', nlay_fit, &
!               start_layer, end_layer
!          DO i = ozfit_start_index, ozfit_end_index
!            WRITE(founit, '(1p60d10.2)') &
!                 avg_kernel(i, ozfit_start_index:ozfit_end_index)
!          END DO
!        ENDIF
!
!        IF (ozwrtfavgk) THEN
!          WRITE(founit, '(A, 3I5)') 'Average kernel1: ', nlay_fit, &
!               start_layer, end_layer
!          DO i = ozfit_start_index, ozfit_end_index
!            WRITE(founit, '(1p60d10.2)') &
!                 favg_kernel(i, ozfit_start_index:ozfit_end_index)
!          END DO
!        ENDIF
!
!        IF (ozwrtcontri) THEN
!          WRITE(founit, '(A, I5)') 'Contribution Function: ', n_rad_wvl
!          DO i = 1, n_rad_wvl
!            WRITE(founit, '(1p60d10.2)') &
!                 contri(ozfit_start_index:ozfit_end_index, i)
!          ENDDO
!        ENDIF
!
!        IF (ozwrtres) THEN
!          WRITE(founit, '(A, I5)') 'Fit residual: ', n_rad_wvl
!          DO i = 1, n_rad_wvl
!            WRITE(founit, '(f9.4, 1p2d12.4)') &
!                 fitwavs(i), fitspec_rad(i), simrad(i)
!          END DO
!        ENDIF
!
!      END IF
!      WRITE(founit, *)
!    END IF
!
!    RETURN
!  END SUBROUTINE gome_write_intermed

  SUBROUTINE l2_ascii_create (l2_filename,l2funit,errstat)
   USE OMSAO_errstat_module
   USE OMSAO_variables_module, ONLY:nxbin, nybin
   USE ozprof_data_module, ONLY: algorithm_name, algorithm_version
   USE m_utilities, ONLY:timestamp
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: l2funit
   CHARACTER (len=*), INTENT(IN)  :: l2_filename
   INTEGER, INTENT(OUT):: errstat
   CHARACTER (LEN=24) :: currtime

   
   open (UNIT=l2funit, FILE=trim(adjustl(l2_filename)), STATUS='UNKNOWN', &
        IOSTAT=errstat)
 
   call timestamp(currtime)
   write(l2funit, '(3A,1x,A27,A10,I5,A10,I5)') &
   trim(adjustl(algorithm_name)), ', ', &
   trim(adjustl(algorithm_version)), currtime, ' xbin = ', nxbin, ' ybin = ', nybin
 
   RETURN
  END SUBROUTINE

  SUBROUTINE l2_ascii_close (l2funit, errstat)
    USE m_utilities, ONLY:timestamp
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: l2funit
    INTEGER, INTENT(OUT):: errstat
    CHARACTER (LEN=24) :: currtime

    errstat = 0
    CALL timestamp(currtime)
    WRITE (l2funit, '(A27)') currtime
    CLOSE ( l2funit )
    RETURN
  END SUBROUTINE

  SUBROUTINE omi_write_intermed (founit, fitcol, dfitcol, exval )

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: maxloc
    !USE OMSAO_indices_module,    ONLY: ring_idx
    USE ozprof_data_module,      ONLY: ozprof, ozprof_std, &
         ozprof_ap, ozprof_apstd, eff_alb, eff_alb_init, nlay, ozdfs, &
         ozinfo, num_iter, ncovar, covar, contri, avg_kernel, use_lograd, &
         nalb, atmosprof, ntp, nlay_fit, ozfit_start_index, &
         ozfit_end_index, start_layer, end_layer, the_ctp, the_cfrac, &
         the_cod, lambcld_refl, do_lambcld, ozprof_nstd, strataod, stratsca, &
         tropaod, tropsca, aerwavs, actawin, ozwrtcorr, ozwrtcovar, &
         ozwrtcontri, the_ai, the_cld_flg, ozwrtres, ozwrtavgk, ozwrtvar, &
         fgasidxs, ngas, nfgas, gaswrt, tracegas, saa_flag, nsaa_spike, &
         ozwrtfavgk, favg_kernel, radcalwrt, nsfc, ozwrtwf, ozwrtsnr, &
         weight_function, do_simu, glintprob,the_snowice, the_landwater_flg, the_glint_flg, &
         allrms, allradrms
         !, ozprof_end_index, &
         !ozprof_start_index, which_alb, which_cld, &
         !use_oe, maxawin, ozprof_init
    USE OMSAO_variables_module, ONLY : & 
         the_sza_atm, the_vza_atm, the_aza_atm, &
         the_lons, the_lats, nloc, the_utc, fitvar_rad, &
         mask_fitvar_rad, n_fitvar_rad, fitvar_rad_std, n_rad_wvl, &
         fitspec_rad, fitres_rad, fitwavs, &
         fitvar_rad_str, fitvar_rad_nstd, &
         simspec_rad, clmspec_rad, actspec_rad, fitweights, &
         numwin, the_pix, the_line

    IMPLICIT NONE

    ! ===============
    ! Input Variables
    ! ===============
    INTEGER,        INTENT (IN) :: founit, exval
    REAL (KIND=dp), DIMENSION(3), INTENT (IN)    :: fitcol
    REAL (KIND=dp), DIMENSION(3, 2), INTENT (IN) :: dfitcol


    ! ===============
    ! Local variables
    ! ===============
    INTEGER                                  :: i, j
    REAL (KIND=dp), DIMENSION (2*maxloc)     :: latlon
    REAL (KIND=dp), DIMENSION (n_fitvar_rad) :: correl

    !REAL (KIND=dp), DIMENSION (n_rad_wvl)    :: simrad

    WRITE(founit, '(A11,2I5,1X,A28)') 'Line/XPix: ', the_line, the_pix, the_utc
    DO i = 1, nloc
      latlon(2*i-1) = the_lats(i)
      latlon(2*i) = the_lons(i)
    END DO
    WRITE(founit,'(13F8.2)') &
         latlon(1:2*nloc), the_sza_atm,the_vza_atm,the_aza_atm
    WRITE(founit, '(A)') &
         'Exit Status, # Iterations, SAA, # SAA Spike, Land/Water, Glint, Snow/Ice Glint Probability'
    WRITE(founit, '(2I4, L4, 4I4, F5.2)') &
         exval, num_iter, saa_flag, nsaa_spike,   &
         the_landwater_flg, the_glint_flg, the_snowice, glintprob

    ! rms = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
    !      fitweights(1:n_rad_wvl))**2.0)/n_rad_wvl)

     simspec_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) - fitres_rad(1:n_rad_wvl)
     IF (use_lograd) THEN
        fitspec_rad(1:n_rad_wvl) = filter_exp(fitspec_rad(1:n_rad_wvl))
        simspec_rad(1:n_rad_wvl) = filter_exp(simspec_rad(1:n_rad_wvl))
        fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) - simspec_rad(1:n_rad_wvl)
     END IF
     ! avgres = SQRT(SUM((ABS(fitres_rad(1:n_rad_wvl)) / &
     !      fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0

    !IF (exval >= 0) THEN
      WRITE(founit, '(A)')    'rms, avgres, dfs, info'
      WRITE(founit, '(14f8.3)') &
          allrms(0), allradrms(0), ozdfs, ozinfo, allrms(1:numwin), allradrms(1:numwin)

      WRITE(founit, '(A8, I3, 1P28E10.2)') &
           'Albedo: ', nalb, eff_alb_init(1:nalb), eff_alb(1:nalb)
      WRITE(founit, '(A8,3F10.3, I3, L3, F10.3)') &
           'Cloud :  ', the_cfrac, the_cod, the_ctp, &
           the_cld_flg, do_lambcld, lambcld_refl
      WRITE(founit, '(A8, I3, 1P36E10.2)') &
           'Aerosol:  ', actawin, aerwavs(1:actawin), tropaod(1:actawin), &
           tropsca(1:actawin), strataod(1:actawin), stratsca(1:actawin), the_ai

      WRITE(founit, '(A31, 3I5)') &
           'Atmosphere and ozone profiles: ', nlay, ntp, nsfc
      WRITE(founit, '(A)') '  #    P(mb)   Z(km)  T(K)     O3-ap O3-apstd   O3     STD     NSTD (DU)'
      WRITE(founit, '(A2, 1x, F9.3, 2F8.3)') ' 0', atmosprof(1:3, 0)
      DO i = 1, nlay 
        WRITE(founit, '(I2, 1X, F9.3, 7F8.3)') i, atmosprof(1:3, i), &
             ozprof_ap(i), ozprof_apstd(i), ozprof(i), ozprof_std(i), &
             ozprof_nstd(i)
      END DO
      WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Total Ozone: ', &
           SUM(ozprof_ap(1:nsfc)), fitcol(1), dfitcol(1, 1), dfitcol(1, 2)
      IF (ntp > 0) THEN
        WRITE(founit, '(A24,4x,F8.3,8X,3F8.3)') ' Stratospheric Ozone: ',&
             SUM(ozprof_ap(1:ntp)), fitcol(2), dfitcol(2, 1), dfitcol(2, 2)
        WRITE(founit, '(A24,4x,F8.3,8X, 3F8.3)') ' Tropospheric Ozone: ',&
             SUM(ozprof_ap(ntp+1:nsfc)), fitcol(3), dfitcol(3, 1), &
             dfitcol(3, 2)          
      ENDIF

      IF (gaswrt) THEN
        WRITE(founit, '(A, I5)') 'Fitted trace gases and uncertainty: ', nfgas
        WRITE(founit, '(A)')  &
             '    Var    Initial   A Priori A Priori Std  VCD      STD       NSTD      AMF      ACFRAC    AVGK(1)   AVGK(2)'
        DO i = 1, ngas
          IF (fgasidxs(i) > 0) THEN
            j = mask_fitvar_rad(fgasidxs(i))
            WRITE(founit, '(I2,1x,A6,1P10d10.2)') &
                 fgasidxs(i), fitvar_rad_str(j), tracegas(i, 1:10)
          ENDIF
        ENDDO
      ENDIF

      IF (ozwrtvar) THEN                      
        WRITE(founit, '(A, I5)') &
             'Fitted variables and uncertainty: ', n_fitvar_rad
        DO i = 1, n_fitvar_rad
          WRITE(founit, '(I2, 1x,A6,1P3d10.2)') i, &
               fitvar_rad_str(mask_fitvar_rad(i)), &
               fitvar_rad(mask_fitvar_rad(i)), &
               fitvar_rad_std(mask_fitvar_rad(i)), &
               fitvar_rad_nstd(mask_fitvar_rad(i))
        END DO
      ENDIF

      IF (ozwrtcorr) THEN
        WRITE(founit, '(A, I5)') &
             'Correlation matrix (Solution Error): ', n_fitvar_rad
        DO i = 1, n_fitvar_rad
          correl(i) = 1.0
          DO j = 1, i - 1
            correl(j) = covar(i, j) / SQRT(covar(i, i) * covar(j, j))
          END DO
          WRITE(founit, '(I2,1X,100f6.2)') i, correl(1:i)
        END DO
      ENDIF

      IF (ozwrtcovar) THEN
        WRITE(founit, '(A, 3I5)') &
             'Covariance matrix (Noise): ', nlay_fit, start_layer, end_layer
        DO i = ozfit_start_index, ozfit_end_index
          WRITE(founit, '(1p60d10.2)') &
               ncovar(i, ozfit_start_index:ozfit_end_index)
        END DO
      ENDIF

      IF (ozwrtavgk) THEN
        WRITE(founit, '(A, 3I5)') &
             'Average kernel: ', nlay_fit, start_layer, end_layer
        DO i = ozfit_start_index, ozfit_end_index
          WRITE(founit, '(1p60d10.2)') &
               avg_kernel(i, ozfit_start_index:ozfit_end_index)
        END DO
      ENDIF

      IF (ozwrtfavgk) THEN
        WRITE(founit, '(A, 3I5)') &
             'Average kernel1: ', nlay_fit, start_layer, end_layer
        DO i = ozfit_start_index, ozfit_end_index
          WRITE(founit, '(1p60d10.2)') &
               favg_kernel(i, ozfit_start_index:ozfit_end_index)
        END DO
      ENDIF

      IF (ozwrtcontri) THEN
        WRITE(founit, '(A, I5)') 'Contribution Function: ', n_rad_wvl
        DO i = 1, n_rad_wvl
          WRITE(founit, '(1p60d10.2)') &
               contri(ozfit_start_index:ozfit_end_index, i)
        ENDDO
      ENDIF

      IF (ozwrtwf) THEN
        WRITE(founit, '(A, 2I5)') &
             'Weighting Function: ', n_rad_wvl, n_fitvar_rad
        DO i = 1, n_rad_wvl
          WRITE(founit, '(1p100d11.3)') weight_function(i, 1:n_fitvar_Rad)
        ENDDO
      ENDIF

      IF (ozwrtsnr) THEN
        WRITE(founit, '(A, I5)') 'Measurement Error: ', n_rad_wvl
        WRITE(founit, '(1p100d11.3)') fitweights(1:n_rad_wvl)
      ENDIF

      IF (ozwrtres) THEN
        WRITE(founit, '(A, I5)') 'Fit residual: ', n_rad_wvl
        DO i = 1, n_rad_wvl
          WRITE(founit, '(f9.4, 1p2d12.4)') &
               fitwavs(i), fitspec_rad(i), simspec_rad(i)
        END DO
      ENDIF

      !IF (ozwrtres) THEN
      !   WRITE(founit, '(A, I5)') 'Ring+Fit residual: ', n_rad_wvl
      !   DO i = 1, n_rad_wvl
      !      WRITE(founit, '(f9.4, 1p3d12.4)') fitwavs(i), fitspec_rad(i), &
      !      simspec_rad(i), database(ring_idx, refidx(i))
      !   END DO
      !ENDIF

      IF (radcalwrt .AND. .NOT. do_simu) THEN
        WRITE(founit, '(A, I5)') 'Radiance Calibration: ', n_rad_wvl
        DO i = 1, n_rad_wvl
          WRITE(founit, '(f9.4, 1p4d12.4)') &
               fitwavs(i), fitspec_rad(i), simspec_rad(i), &
               clmspec_rad(i), actspec_rad(i)
          write(*,*) fitwavs(i), fitspec_rad(i), simspec_rad(i), &
               clmspec_rad(i), actspec_rad(i)
        END DO
      ENDIF

    !END IF

    RETURN
  END SUBROUTINE omi_write_intermed

end module ascii_output_module
