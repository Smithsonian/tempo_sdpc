!
module m_prepare_databases

  public prepare_databases
  private

contains

  ! *********************** Modification History *******************
  ! xiong liu, July 2003
  ! 1. Add vgr, vgl, hwl, hwr for voigt profile shape
  ! 2. Add one argument to subroutine undersample
  ! ****************************************************************

  SUBROUTINE prepare_databases (n_rad_wvl, curr_rad_wvl, pge_error_status )

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY: phase, n_refwvl, refwvl, database, &
         database_shiwf, do_bandavg, nradpix, numwin, lo_radbnd, up_radbnd, &
         n_refspec_pts, curr_sol_spec, i0sav, refidx_sav, database_save, &
         n_refwvl_sav, refwvl_sav, curr_sol_spec, nsolpix, refsol_idx, &
         radnhtrunc, refnhextra!, n_irrad_wvl, refspec_orig_data, &
         !fitvar_rad_str, have_undersampling
    USE OMSAO_indices_module,   ONLY: max_rs_idx, solar_idx, shift_offset, &
      wvl_idx, ring1_idx, ring_idx, sdc_idx, &
      com_idx, com1_idx, com2_idx, com3_idx
    USE OMSAO_errstat_module,   ONLY: pge_errstat_error
    USE ozprof_data_module,     ONLY: ring_on_line, ozprof_flag, do_tracewf!, &
         !radcalwrt, do_simu
    use avg_band, only: avg_band_refspec
    use m_dataspline, only: dataspline
    use m_undersample
    use m_prepare_refspecs, only: prepare_refspecs

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                               INTENT (IN) :: n_rad_wvl
    REAL (KIND=dp), DIMENSION (n_rad_wvl), INTENT (IN) :: curr_rad_wvl

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER, INTENT (OUT) :: pge_error_status

    ! Local variable
    INTEGER        :: ntemp, i, j, fidx, lidx, fidx1, lidx1, &
         nextra, nhsolextra, idxoff!, npts, k, choice
    REAL (KIND=dp) :: deltlam

    ! ensure reference spectra having larger wave range than required
    ! and avoid interpolation failure to get garbage values
    ! in gome_pge_fitting_process, special processing is made to enusre
    ! n_rad_wvl is smaller than n_irrad_wvl by radnhextra * 2 in each 
    ! fitting window
    nextra    = refnhextra * 2
    nhsolextra = radnhtrunc - refnhextra
    n_refwvl = n_rad_wvl + numwin * nextra

    ! Obtain reference position in solar spectra
    fidx = 1
    DO i = 1, numwin
      lidx = fidx + nradpix(i) + nextra - 1
      idxoff = nhsolextra * 2 * i - nhsolextra
      refsol_idx(fidx:lidx) = (/(j, j = fidx + idxoff, lidx + idxoff)/)
      fidx = lidx  + 1
    ENDDO

    fidx = 1; fidx1 = 1
    DO i = 1, numwin
      j = i - 1
      deltlam = &
           curr_sol_spec(wvl_idx, fidx1 + 1) - curr_sol_spec(wvl_idx, fidx1)
      lidx = fidx + nradpix(i) + nextra - 1 ; lidx1 = fidx + nsolpix(i) - 1

      IF (refnhextra > 0) refwvl(fidx:fidx+refnhextra-1) = &
           curr_sol_spec(wvl_idx, refsol_idx(fidx:fidx+refnhextra-1))
      IF (refnhextra > 0) refwvl(lidx-refnhextra+1:lidx) = &
           curr_sol_spec(wvl_idx, refsol_idx(lidx-refnhextra+1:lidx))
      refwvl(fidx+refnhextra:lidx-refnhextra) = &
           curr_rad_wvl(fidx - j * nextra : lidx - i * nextra)

      fidx = lidx + 1; fidx1 = lidx1 + 1
    ENDDO

    fidx = 1
    DO i = 1, numwin
      lidx = fidx + nradpix(i) - 1
      idxoff = refnhextra + (i - 1) * nextra
      refidx_sav(fidx:lidx) = (/(j, j = fidx + idxoff, lidx + idxoff)/)
      fidx = lidx  + 1
    ENDDO
    refwvl_sav(1:n_refwvl) = refwvl(1:n_refwvl); n_refwvl_sav = n_refwvl


    ! Initialize database
    database = 0.0; database_shiwf = 0.0

    ! --------------------------------------
    ! Calculate the splined fitting database
    ! --------------------------------------
    CALL prepare_refspecs (n_refwvl, refwvl(1:n_refwvl), pge_error_status)
   
    IF ( pge_error_status >= pge_errstat_error ) RETURN
    i0sav(1:n_refwvl) = database(solar_idx, 1:n_refwvl)


    ! ---------------------------------------------------------
    ! Spline external reference spectra to common radiance grid
    ! ---------------------------------------------------------
    CALL dataspline ( n_refwvl, refwvl(1:n_refwvl), pge_error_status)

    
    IF ( pge_error_status >= pge_errstat_error) RETURN

    ! ----------------------------------------------------------
    ! Calculate the undersampled spectrum
    ! -----------------------------------------------------------
    CALL undersample (n_refwvl, refwvl(1:n_refwvl), phase, pge_error_status)

   
    IF ( pge_error_status >= pge_errstat_error ) RETURN


    IF (do_bandavg) THEN
      DO i = 1, max_rs_idx 
        IF ( (i == ring_idx .OR. i == ring1_idx ) &
             .AND. ozprof_flag .AND. ring_on_line) CYCLE 

        IF ( i == com_idx .OR. i == com1_idx .OR. i == com2_idx .OR. i == com3_idx .OR. i == sdc_idx ) CYCLE

        IF (n_refspec_pts(i) > 0 ) THEN
          CALL avg_band_refspec(refwvl(1:n_refwvl), database(i,1:n_refwvl), &
               n_refwvl, ntemp, pge_error_status)
          IF ( pge_error_status >= pge_errstat_error ) RETURN
        ENDIF

        j = shift_offset + i
        IF (n_refspec_pts(i) > 0 .AND. lo_radbnd(j) < up_radbnd(j)) THEN
          CALL avg_band_refspec(refwvl(1:n_refwvl), &
               database_shiwf(i,1:n_refwvl), n_refwvl, ntemp, &
               pge_error_status)    

          IF ( pge_error_status >= pge_errstat_error ) RETURN
        ENDIF
      ENDDO

      CALL avg_band_refspec(refwvl(1:n_refwvl), refwvl(1:n_refwvl), &
           n_refwvl, ntemp, pge_error_status)
      IF ( pge_error_status >= pge_errstat_error ) RETURN
      n_refwvl = ntemp  
    ENDIF


    ! save it for fitting weighting function
    IF (do_tracewf) database_save = database


    RETURN
  END SUBROUTINE prepare_databases

end module m_prepare_databases
