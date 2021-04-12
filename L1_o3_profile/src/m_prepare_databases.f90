!
module m_prepare_databases

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,  ONLY: zerospec_string
    USE OMSAO_variables_module, ONLY: phase, have_undersampling, &
         do_bandavg, numwin, lo_radbnd, up_radbnd, &         
         refnhextra,radnhtrunc, refnhextra,  &       
         curr_sol_spec, curr_rad_spec_ori ,nsolpix, n_irrad_wvl,nradpix, & ! measured spectra
         refspec_orig_data, n_refspec_pts,  & ! refs
         n_refspec_pts,  refspec_norm, refspec_fname, &
         sring_fidx, sring_lidx, nsol_ring, sol_spec_ring, & ! ringspec
         yn_varyslit,slit_rad, which_slit, instrument_sidx, & ! slit variables
         solwinfit,nslit, nslit_rad, nslit_sol, slitwav, &
         slitwav_sol, slitwav_rad, slitfit, solslitfit, radslitfit, &
         i0sav, refidx_sav,refsol_idx, n_refwvl, refwvl, n_refwvl_sav, refwvl_sav, & ! output
         database, database_shiwf, database_save, curr_radresponse_spec ! output

    USE OMSAO_indices_module,   ONLY: max_rs_idx, solar_idx, shift_offset, &
         hwe_idx, hwr_idx, hwl_idx, asy_idx, spk_idx, &
         wvl_idx,spc_idx, ring1_idx, ring_idx, us1_idx, us2_idx, &        
         com_idx, com1_idx, com2_idx, com3_idx, vege_idx, &
         so2_idx, so2v_idx, bro_idx, bro2_idx, &
         hcho_idx, no2_t1_idx, no2_t2_idx,  o2o2_idx, chloro_idx, &
         o2_idx, o2t2_idx, h2o_idx, h2ot2_idx, lh2o_idx, rsl_idx

    USE OMSAO_errstat_module
    USE ozprof_data_module,     ONLY: ring_on_line, ozprof_flag, nsl, ring_convol, which_inr, &
                                      use_o4dtcrs, use_so2dtcrs, use_h2odptcrs, use_o2dptcrs

    USE m_ezspline_interpolation, ONLY:  interpolation, bspline, bspline1
    USE m_convol, ONLY: convol, convol_i0
    use m_avg_band, only: avg_band_refspec


    public prepare_databases,prepare_refspecs, undersample, dataspline, append_solring

contains

  ! 1) define refwvl
    ! ensure reference spectra having larger wave range than 
    ! required and avoid interpolation failure to get garbage values
    ! special processing is made to enusre
    !   n_rad_wvl is smaller than n_irrad_wvl by radnhextra * 2 in each 
    ! fitting window
  ! 2) interpolate irradiance or solor reference onto refwvl
  ! 3) interpolate/convolve reference spectrum onto refwvl
  ! 4) calculate undersample spectrum
  ! 5) if (avg_band) then re-sample above spectrum into avg wavel
  ! *********************** Modification History *******************
  ! xiong liu, July 2003
  ! 1. Add vgr, vgl, hwl, hwr for voigt profile shape
  ! 2. Add one argument to subroutine undersample
  ! ****************************************************************

  SUBROUTINE prepare_databases (n_rad_wvl, curr_rad_wvl, pge_error_status )

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

    ! --------------------------------------
    ! define refwvl
    ! --------------------------------------
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
      deltlam= curr_sol_spec(wvl_idx, fidx1 + 1) - curr_sol_spec(wvl_idx, fidx1)
      lidx = fidx + nradpix(i) + nextra - 1 
      lidx1 = fidx + nsolpix(i) - 1

      ! jbak , 'discontinuty occures if there is large shift btw rad/irrad wave'
      !IF (refnhextra > 0) refwvl(fidx:fidx+refnhextra-1) = &
      !     curr_rad_spec_ori(wvl_idx, refsol_idx(fidx:fidx+refnhextra-1))
      !IF (refnhextra > 0) refwvl(lidx-refnhextra+1:lidx) = &
      !     curr_rad_spec_ori(wvl_idx, refsol_idx(lidx-refnhextra+1:lidx))

      IF (refnhextra > 0) refwvl(fidx:fidx+refnhextra-1) = &
           curr_sol_spec(wvl_idx, refsol_idx(fidx:fidx+refnhextra-1))
      IF (refnhextra > 0) refwvl(lidx-refnhextra+1:lidx) = &
           curr_sol_spec(wvl_idx, refsol_idx(lidx-refnhextra+1:lidx))

      refwvl(fidx+refnhextra:lidx-refnhextra) = &
           curr_rad_wvl(fidx - j * nextra : lidx - i * nextra)

      fidx = lidx + 1; fidx1 = lidx1 + 1
    ENDDO

    IF (ANY( (refwvl(2:n_refwvl)-refwvl(1:n_refwvl-1)) < 0 )) THEN 
      WRITE(*,*) 'prepare_database: discontinuty in refwvl '
      !DO i = 1, n_refwvl-1
      !     write(*,'(i4, 5f8.2)') i,curr_rad_spec_save(1, i), curr_sol_spec(1,i), refwvl(i),refwvl(i+1) - refwvl(i)
      !ENDDO
      pge_error_status = pge_errstat_error 
      RETURN
    ENDIF


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
    IF (refwvl(1) < curr_sol_spec(1, 1) .OR. & 
      refwvl(n_refwvl) > curr_sol_spec(1,n_irrad_wvl) ) THEN 
      WRITE(*,*) 'prepare_database: refwvl out of curr_sol_wavl '
      WRITE(*,*) 'refwvl', refwvl(1), radnhtrunc, refnhextra
      WRITE(*,*) 'sol', curr_sol_spec(1, 1)
      pge_error_status = pge_errstat_error 
      RETURN
    ENDIF
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

        IF ( i == com_idx .OR. i == com1_idx .OR. i == com2_idx .OR. i == com3_idx  ) CYCLE

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
    database_save = database
    RETURN
  END SUBROUTINE prepare_databases

  SUBROUTINE dataspline ( n_radwvl, curr_rad_wvl, errstat)


    IMPLICIT NONE

    INTEGER,                              INTENT (IN)  :: n_radwvl
    REAL (KIND=dp), DIMENSION (n_radwvl), INTENT (IN)  :: curr_rad_wvl
    INTEGER,                              INTENT (OUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: j, fidx, lidx, stidx,  idx, npts
    REAL (KIND=dp), DIMENSION (:), POINTER :: specmod
    REAL (KIND=dp)                         :: frefw, lrefw, scalex

    ! ------------------
    ! External functions
    ! ------------------
    !INTEGER :: OMI_SMF_setmsg

    CHARACTER (LEN=11), PARAMETER :: modulename = 'dataspline'

    errstat = pge_errstat_ok

    ! ---------------------------------------------------------------------
    ! Load results into the database array. The order of the spectra is
    ! determined by the molecule indices in OMSAO_indices_module. Only 
    ! those spectra that are requested for read-in are interpolated.
    ! Note that this might be different from the actual fitting parameters
    ! to be varied in the fit - a non-zero but constant fitting parameter
    ! will still require a spectrum to make a contribution.
    !
    ! If the original reference spectrum does not cover the wavelength
    ! range of the current radiance wavelength, we interpolate only the
    ! part that is covered and set the rest to Zero.
    ! ---------------------------------------------------------------------
    ! ------------------------
    ! Spline Reference Spectra
    ! ------------------------
    IF (slit_rad) THEN   ! use derived slit from solar for ring effect
      nslit  = nslit_sol; slitwav = slitwav_sol; slitfit = solslitfit
    END IF

    ! 1: solar_idx, obtained latter in prepare_refspecs
    ! 2: ring_idx, done on line
    IF (ring_on_line) THEN
      stidx = ring_idx + 1
    ELSE 
      stidx = ring_idx
    ENDIF

    ! Perform solar i0 effect (no need to convolve) 
    DO idx = stidx, max_rs_idx
      IF ( n_refspec_pts(idx) >= 3 .AND. INDEX(TRIM(ADJUSTL(refspec_fname(idx))), &
           zerospec_string ) == 0) THEN
           npts = n_refspec_pts(idx)  ! Define short-hand
           allocate (specmod(npts))
           specmod(1:npts) = refspec_orig_data(idx,1:npts,spc_idx)
        !IF (idx == o2o2_idx) refspec_orig_data(idx,1:npts,wvl_idx) = refspec_orig_data(idx,1:npts,wvl_idx)
        IF (idx == bro_idx .OR. idx == bro2_idx .OR. idx == no2_t1_idx .OR. idx == no2_t2_idx .OR. &
            idx == so2_idx .OR. idx == so2v_idx .OR. idx == o2o2_idx .OR. &
            idx == o2_idx .OR. idx == o2t2_idx .OR. idx == h2o_idx .OR. idx == h2ot2_idx .OR. &
            idx == hcho_idx ) THEN
            IF ((idx == h2o_idx .or. idx == h2ot2_idx) .and. use_h2odptcrs) cycle  
            IF ((idx == o2_idx .or. idx == o2t2_idx) .and. use_o2dptcrs) cycle  
            IF (idx == o2o2_idx  .and. use_o4dtcrs) cycle  
            IF (idx == so2_idx  .and. use_so2dtcrs) cycle  
             
          IF (idx == bro_idx .OR. idx == bro2_idx) THEN
              scalex = 2.0E13 ! 1.0E-4
          ELSE IF (idx == o2o2_idx) THEN
              scalex = 2.6E33
          ELSE IF (idx == o2_idx .OR. idx == o2t2_idx) THEN
              scalex = 6.0E24
          ELSE IF (idx == h2o_idx .OR. idx == h2ot2_idx) THEN
              scalex = 1.0E23
          ELSE
              scalex = 5.0E15 ! 20E-4
          ENDIF

          scalex = scalex * refspec_norm(idx)
          CALL convol_i0(refspec_orig_data(idx,1:npts,wvl_idx),specmod(1:npts), npts,scalex)
          IF (errstat == pge_errstat_error) RETURN

          ! -----------------------------------------------------------------
          ! Call interpolation and check returned error status. 
          ! WARNING status indicates missing parts of the interpolated 
          ! spectrum, while ERROR status indicates a
          ! more serious condition that requires termination.
          ! -----------------------------------------------------------------  
        ELSE IF ((idx /= com_idx .AND. idx /= com1_idx .AND. idx /= com2_idx .AND. idx /= com3_idx  &
             .AND. idx /= ring_idx .AND. idx /= ring1_idx &
             .AND. idx /= vege_idx .AND. idx /= chloro_idx .AND. idx /= lh2o_idx ) .OR. &
             (idx == ring_idx .AND. ring_convol) .OR. (idx == ring1_idx .AND. ring_convol)) THEN

          CALL convol (refspec_orig_data(idx,1:npts, wvl_idx),specmod(1:npts),npts)

        ENDIF

        j = shift_offset + idx

        frefw = refspec_orig_data(idx,1, wvl_idx)
        lrefw = refspec_orig_data(idx,npts, wvl_idx)
        fidx = MINVAL(MINLOC(curr_rad_wvl, MASK=(curr_rad_wvl >= &
             frefw + 0.02 .AND. curr_rad_wvl <= lrefw - 0.02)))
        lidx = MINVAL(MAXLOC(curr_rad_wvl, MASK=(curr_rad_wvl >= &
             frefw + 0.02 .AND. curr_rad_wvl <= lrefw - 0.02)))

        IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN 
          IF (lo_radbnd(j) < up_radbnd(j)) THEN

            ! xliu: 02/28/2009
            ! Save convolved but at original reference wavelength grid for 
            ! further interpolation 
            ! Necessary when fitting a wavelength shift 
            refspec_orig_data(idx,1:npts,3) =  specmod(1:npts)

            CALL bspline1(refspec_orig_data(idx,1:npts,wvl_idx), &
                 specmod(1:npts),&
                 npts, curr_rad_wvl(fidx:lidx), database(idx, fidx:lidx), &
                 database_shiwf(idx, fidx:lidx), lidx - fidx + 1, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': BSPLINE1 error'
              errstat = pge_errstat_error
            ENDIF
          ELSE
            CALL bspline(refspec_orig_data(idx,1:npts,wvl_idx), &
                 specmod(1:npts), &
                 npts, curr_rad_wvl(fidx:lidx), database(idx, fidx:lidx), &
                 lidx - fidx + 1, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename,': BSPLINE error'
              errstat = pge_errstat_error
            ENDIF
          END IF
        ENDIF
        deallocate(specmod)
      END IF

!      IF (slit_rad .AND. i == 2 ) THEN   ! use derived slit for other reference spectra
      IF (slit_rad ) THEN   ! use derived slit for other reference spectra
        nslit  = nslit_rad
        slitwav = slitwav_rad
        slitfit = radslitfit
      END IF
    END DO
    RETURN
  END SUBROUTINE dataspline

  SUBROUTINE prepare_refspecs (n_radpts, curr_rad_wvl, pge_error_status )

    ! ***********************************************************
    !
    !   Calculate the splined fitting database.
    !   Note that the undersampled spectrum has just been done.
    !   Finish filling in database array.
    !
    ! ***********************************************************

    IMPLICIT NONE

    ! *******************************************************************
    ! CAREFUL: Assumes that radiance and solar wavelength arrays have the
    ! same number of points. That must not be the case if we read in a
    ! general EL1 file. Examine and adjust! (tpk, note to himself)
    ! *******************************************************************

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                              INTENT (IN)    :: n_radpts
    REAL (KIND=dp), DIMENSION (n_radpts), INTENT (IN)    :: curr_rad_wvl
    INTEGER,                              INTENT (INOUT) :: pge_error_status

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: errstat
    REAL  (KIND=dp), DIMENSION (n_radpts) :: spline_sun
    real (kind=dp), dimension(n_irrad_wvl) :: tmp_x_in, tmp_y_in

     errstat = pge_errstat_ok
     
    ! Spline irradiance spectrum onto radiance grid
    ! PROBLEM: first wavlength in POS is smaller than first wavelength in
    !          curr_sol_spec(wvl_idx,*).
    ! SOLUTION: Don't include POS(1) in the interpolation, and assign 
    !           SPLINE_SUN(1) = SPEC_SUN(1)

    

   ! FIXME - masking array temporaries
    tmp_x_in=curr_sol_spec(wvl_idx,1:n_irrad_wvl)
    tmp_y_in=curr_sol_spec(spc_idx,1:n_irrad_wvl)
    CALL interpolation (n_irrad_wvl,tmp_x_in, tmp_y_in, &
         n_radpts, curr_rad_wvl, spline_sun, errstat )  
    IF (which_inr == 1) THEN
        IF (.NOT. allocated(curr_radresponse_spec)) THEN 
           WRITE(*,*) 'curr_radresponse_spec is not allocated'; stop 1
        ENDIF
        CALL interpolation (n_irrad_wvl, curr_radresponse_spec(wvl_idx,1:n_irrad_wvl),        &
        curr_radresponse_spec(spc_idx,1:n_irrad_wvl), n_radpts,curr_rad_wvl(1:n_radpts), &
        database(rsl_idx, 1:n_radpts), errstat )
        deallocate(curr_radresponse_spec) 
    ENDIF 
    ! Save to database
    database(solar_idx, 1:n_radpts) = spline_sun(1:n_radpts) 

    pge_error_status = MAX ( errstat, pge_error_status )

    RETURN
  END SUBROUTINE prepare_refspecs
 ! *********************** Modification History ********
  ! xliu: 
  ! 1. Add call to asym_gauss_vary if vary slit width 
  !    option is set
  ! 2. Add variables hw1earr, slitwav, e_asymarr, n_slit_pts
  !    slitwav,  from OMSAO_variables_module
  ! 3. Add arguments vgl, vgr, hwl, hwr
  ! 4. Use voigt_gauss to replace asym_gauss
  ! *****************************************************
  SUBROUTINE undersample (n_rad_pts, curr_wvl, phase, pge_error_status )

    !  Convolves input spectrum with Gaussian slit function of specified
    !  HW1e, and samples at a particular input phase to give the OMI
    !  undersampling spectrum. This version calculates both phases of the
    !  undersampling spectrum, phase1 - i.e., underspec (1, i) - being the
    !  more common in OMI spectra.

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                               INTENT (IN) :: n_rad_pts
    REAL (KIND=dp),                        INTENT (IN) :: phase
    REAL (KIND=dp), DIMENSION (n_rad_pts), INTENT (IN) :: curr_wvl

    ! ---------------
    ! Output variable
    ! ---------------
    INTEGER, INTENT (OUT) :: pge_error_status

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=dp), DIMENSION (2,n_rad_pts+4) :: underspec
    REAL (KIND=dp), DIMENSION (:), POINTER    :: locwvl, locspec, specmod, specmod1
    REAL (KIND=dp), DIMENSION (n_rad_pts + 4) :: tmpwav, over, under, resample, & 
          resample1, subwav, tmpspec
    INTEGER :: npts, errstat, iwin, fidx, lidx, npoints

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=11), PARAMETER :: modulename = 'undersample'

    errstat = pge_errstat_ok

    IF (have_undersampling .OR. sring_fidx > 0 .OR. sring_lidx < nsol_ring) &
         THEN

      ! Use derived slit from irradiance
      IF (slit_rad) THEN
        nslit = nslit_sol
        slitwav = slitwav_sol
        slitfit = solslitfit
      ENDIF

      ! ==================================================
      ! Assign solar reference spectrum to local variables
      ! ==================================================
      npts = n_refspec_pts(solar_idx)
      allocate (locwvl(npts), locspec(npts))
      allocate (specmod(npts), specmod1(npts))
      locwvl (1:npts) = refspec_orig_data(solar_idx,1:npts,wvl_idx)
      locspec(1:npts) = refspec_orig_data(solar_idx,1:npts,spc_idx)

      CALL convol (locwvl, locspec,  npts)
      specmod = locspec
      IF (nsl > 0 ) THEN 
        IF (.NOT. yn_varyslit ) THEN
          IF (which_slit == 0) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
             solwinfit(1:numwin, hwe_idx,1) + 0.01
          ELSE IF (which_slit == 1) THEN 
             solwinfit(1:numwin, hwe_idx:asy_idx,1) = &
             solwinfit(1:numwin, hwe_idx:asy_idx,1) + 0.01
          ELSE IF (which_slit == 2) THEN
             solwinfit(1:numwin, hwl_idx:hwr_idx,1) = &
             solwinfit(1:numwin, hwl_idx:hwr_idx,1) + 0.01
          ELSE IF (which_slit == 3) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
             solwinfit(1:numwin, hwe_idx,1) + 0.01
          ELSE IF (which_slit == 4) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
                  solwinfit(1:numwin, hwe_idx,1) + 0.01
             solwinfit(1:numwin, spk_idx,1) = &
                  solwinfit(1:numwin, spk_idx,1) + 0.01
          ELSE IF (which_slit == 5) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
                  solwinfit(1:numwin, hwe_idx,1) + 0.01
             solwinfit(1:numwin, spk_idx,1) = &
                  solwinfit(1:numwin, spk_idx,1) + 0.01
             solwinfit(1:numwin, asy_idx,1) = &
                  solwinfit(1:numwin, asy_idx,1) + 0.01
          ELSE IF (which_slit == instrument_sidx) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
                  solwinfit(1:numwin, hwe_idx,1) + 0.01
             solwinfit(1:numwin, spk_idx,1) = &
                  solwinfit(1:numwin, spk_idx,1) + 0.01
             solwinfit(1:numwin, asy_idx,1) = &
                  solwinfit(1:numwin, asy_idx,1) + 0.01
          ENDIF
        ELSE
          IF (which_slit == 0) THEN
            slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
          ELSE IF (which_slit == 1) THEN 
            slitfit(1:nslit, hwe_idx:asy_idx, 1) = &
            slitfit(1:nslit, hwe_idx:asy_idx, 1) + 0.01
          ELSE IF (which_slit == 2) THEN
           slitfit(1:nslit, hwl_idx:hwr_idx, 1) = &
           slitfit(1:nslit, hwl_idx:hwr_idx, 1) + 0.01
          ELSE IF ( which_slit == 3) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
          ELSE IF ( which_slit == 4) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
           slitfit(1:nslit, spk_idx, 1) = slitfit(1:nslit, spk_idx, 1) + 0.01
          ELSE IF ( which_slit == 5) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
           slitfit(1:nslit, asy_idx, 1) = slitfit(1:nslit, asy_idx, 1) + 0.01
           slitfit(1:nslit, spk_idx, 1) = slitfit(1:nslit, spk_idx, 1) + 0.01
          ELSE IF ( which_slit == instrument_sidx) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) + 0.01
           slitfit(1:nslit, asy_idx, 1) = slitfit(1:nslit, asy_idx, 1) + 0.01
           slitfit(1:nslit, spk_idx, 1) = slitfit(1:nslit, spk_idx, 1) + 0.01
          ENDIF
        ENDIF
        CALL convol (locwvl, locspec,  npts)
        specmod1 = locspec
     
       IF (.NOT. yn_varyslit ) THEN
          IF (which_slit == 0) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
             solwinfit(1:numwin, hwe_idx,1) - 0.01
          ELSE IF (which_slit == 1) THEN 
             solwinfit(1:numwin, hwe_idx:asy_idx,1) = &
             solwinfit(1:numwin, hwe_idx:asy_idx,1) - 0.01
          ELSE IF (which_slit == 2) THEN
             solwinfit(1:numwin, hwl_idx:hwr_idx,1) = &
             solwinfit(1:numwin, hwl_idx:hwr_idx,1) - 0.01
          ELSE IF (which_slit == 3) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
             solwinfit(1:numwin, hwe_idx,1) - 0.01
          ELSE IF ( which_slit == 4) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
                  solwinfit(1:numwin, hwe_idx,1) - 0.01
             solwinfit(1:numwin, spk_idx,1) = &
                  solwinfit(1:numwin, spk_idx,1) - 0.01
          ELSE IF ( which_slit == 5) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
                  solwinfit(1:numwin, hwe_idx,1) - 0.01
             solwinfit(1:numwin, spk_idx,1) = &
                  solwinfit(1:numwin, spk_idx,1) - 0.01
             solwinfit(1:numwin, asy_idx,1) = &
                  solwinfit(1:numwin, asy_idx,1) - 0.01
          ELSE IF ( which_slit == instrument_sidx) THEN
             solwinfit(1:numwin, hwe_idx,1) = &
                  solwinfit(1:numwin, hwe_idx,1) - 0.01
             solwinfit(1:numwin, spk_idx,1) = &
                  solwinfit(1:numwin, spk_idx,1) - 0.01
             solwinfit(1:numwin, asy_idx,1) = &
                  solwinfit(1:numwin, asy_idx,1) - 0.01
          ENDIF
        ELSE
          IF (which_slit == 0) THEN
            slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
          ELSE IF (which_slit == 1) THEN 
            slitfit(1:nslit, hwe_idx:asy_idx, 1) = &
            slitfit(1:nslit, hwe_idx:asy_idx, 1) + 0.01
          ELSE IF (which_slit == 2) THEN
           slitfit(1:nslit, hwl_idx:hwr_idx, 1) = &
           slitfit(1:nslit, hwl_idx:hwr_idx, 1) + 0.01
          ELSE IF ( which_slit == 3) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
          ELSE IF ( which_slit == 4) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
           slitfit(1:nslit, spk_idx, 1) = slitfit(1:nslit, spk_idx, 1) - 0.01
          ELSE IF ( which_slit == 5) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
           slitfit(1:nslit, spk_idx, 1) = slitfit(1:nslit, spk_idx, 1) - 0.01
           slitfit(1:nslit, asy_idx, 1) = slitfit(1:nslit, asy_idx, 1) - 0.01
          ELSE IF ( which_slit == instrument_sidx) THEN
           slitfit(1:nslit, hwe_idx, 1) = slitfit(1:nslit, hwe_idx, 1) - 0.01
           slitfit(1:nslit, spk_idx, 1) = slitfit(1:nslit, spk_idx, 1) - 0.01
           slitfit(1:nslit, asy_idx, 1) = slitfit(1:nslit, asy_idx, 1) - 0.01
          ENDIF
        ENDIF
      ENDIF

      ! Append Ring Source Spectrum
      IF (sring_fidx > 0 .OR. sring_lidx < nsol_ring) THEN
        WRITE(www_lun,'(A)') 'undersample :Perform Append_Solring !!!'
        CALL append_solring(nsol_ring, sring_fidx, sring_lidx, &
             sol_spec_ring(wvl_idx, 1:nsol_ring), &
             sol_spec_ring(spc_idx, 1:nsol_ring), &
             npts, locwvl(1:npts), specmod(1:npts), pge_error_status)
        IF (pge_error_status > pge_errstat_warning) RETURN
       ENDIF
     ENDIF


    IF (.NOT. have_undersampling) RETURN

    ! do it separately for each window
    lidx = 0
    DO iwin = 1, numwin
!      npoints = nradpix(iwin) + 8  ! 
      npoints = nradpix(iwin) + refnhextra*2 + 4 ! JBAK reference has 4 more wavelengths and add 4 extra pixels 
      fidx =  lidx + 1            ! fidx:lidx refers to position in curr_wvl
      lidx =  fidx + npoints - 5

      IF (iwin == 2) THEN
        database(us1_idx, fidx:lidx) = 0.0
        database(us2_idx, fidx:lidx) = 0.0
        CYCLE
      ENDIF

      subwav(3:npoints-2) = curr_wvl(fidx:lidx)

      ! Add extra wavelengths
      subwav(2) = 2 * subwav(3) - subwav(4)
      subwav(1) = 2 * subwav(2) - subwav(3)
      subwav(npoints-1) = 2 * subwav(npoints-2) - subwav(npoints-3)
      subwav(npoints)   = 2 * subwav(npoints-1) - subwav(npoints-2)

      ! Phase1 calculation: Calculate spline derivatives for KPNO data
      !                     Calculate solar spectrum at OMI positions 
      CALL interpolation (npts, locwvl(1:npts), specmod(1:npts), npoints, &
           subwav(1:npoints), resample(1:npoints),  errstat)

      IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 

        pge_error_status = pge_errstat_error
        RETURN
      END IF

      ! weighting function for slit width
      IF (nsl > 0) THEN
        CALL interpolation (npts, locwvl(1:npts), specmod1(1:npts), npoints, &
             subwav(1:npoints),  resample1(1:npoints), errstat)

        IF ( errstat > pge_errstat_warning ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) 

          pge_error_status = pge_errstat_error
          RETURN
        END IF
        database_shiwf(1,fidx:lidx) = &
             (resample1(3:npoints-2) - resample(3:npoints-2)) / 0.01
      ENDIF

      ! Calculate solar spectrum at OMI + phase positions, 
      ! original and resampled.
      ! ----------------------------------------------------------------------
      ! The original ("modified K.C.) scheme to compute the UNDERSPEC 
      ! wavelength array
      ! ----------------------------------------------------------------------
      ! ( assumes ABS(PHASE) < 1.0 )
      ! ----------------------------
      tmpwav(1:npoints-1) = &
           (1.0-phase) * subwav(1:npoints-1) + phase * subwav(2:npoints)

      tmpwav(1) = subwav(1)
      tmpwav(npoints)   = subwav(npoints)
      IF ( tmpwav(2) <= tmpwav(1) ) tmpwav(2) = (tmpwav(1)+tmpwav(3))/2.0

      CALL interpolation (  npts, locwvl(1:npts), specmod(1:npts), npoints, &
           tmpwav(1:npoints), over(1:npoints), errstat )

      IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF

      CALL interpolation (npoints, subwav(1:npoints), resample(1:npoints),&
           npoints, tmpwav(1:npoints), under(1:npoints), errstat )

      IF ( errstat /= pge_errstat_ok ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF

      underspec(1,1:npoints) = over(1:npoints) - under(1:npoints)

      ! --------------------------------------------------------------
      ! Phase2 calculation: Calculate solar spectrum at OMI positions, 
      ! original and resampled.
      ! -------------------------------------------------------------- 
      tmpspec(1:npoints)   = resample(1:npoints)
      resample (1:npoints) = over(1:npoints)
      over(1:npoints)      = tmpspec(1:npoints)

      CALL interpolation (npoints, tmpwav(1:npoints), &
           resample(1:npoints), npoints, subwav(1:npoints), &
           under(1:npoints), errstat )
      IF ( errstat /= pge_errstat_ok ) THEN    
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
        pge_error_status = pge_errstat_error
        RETURN
      END IF

      ! ====================================
      ! Compute final undersampling spectrum
      ! ====================================
      underspec(2,1:npoints) = over(1:npoints) - under(1:npoints)

      ! check for missing F and thus spikes
      !dwavmax = subwav(2) - subwav(1)
      !DO i = 2, npoints
      !   IF (subwav(i) - subwav(i-1) > dwavmax * 1.2) THEN
      !      underspec(1, i-1:i) = 0.0
      !      underspec(2, i-1:i) = 0.0
      !   ENDIF
      !ENDDO

      refspec_orig_data(us1_idx, fidx:lidx, wvl_idx) = &
           tmpwav   (   3:npoints-2)
      refspec_orig_data(us1_idx, fidx:lidx, spc_idx) = &
           underspec(1, 3:npoints-2)
      refspec_orig_data(us2_idx, fidx:lidx, wvl_idx) = &
           subwav   (   3:npoints-2)
      refspec_orig_data(us2_idx, fidx:lidx, spc_idx) = &
           underspec(2, 3:npoints-2)

      database(us1_idx, fidx:lidx) = underspec(1, 3:npoints-2)
      database(us2_idx, fidx:lidx) = underspec(2, 3:npoints-2)    

    ENDDO ! end window loop

    n_refspec_pts (us1_idx)  = n_rad_pts
    n_refspec_pts (us2_idx)  = n_rad_pts

    !WRITE(90, *) n_rad_pts, nradpix(1:numwin) + 4
    !DO i = 1, n_rad_pts
    !   WRITE(91, *) curr_wvl(i), database(us1_idx, i),  database(us2_idx, i)
    !ENDDO
    !stop 1

    ! Use derived slit from radiance later on if it is available
    IF (slit_rad) THEN
      nslit = nslit_rad
      slitwav = slitwav_rad
      slitfit = radslitfit
    ENDIF

    deallocate(locwvl, locspec, specmod, specmod1)

    RETURN
  END SUBROUTINE undersample
  
  SUBROUTINE append_solring(nspec, n1, n2, wav, spec, nsol, solwav, solspec, errstat)

    IMPLICIT NONE  

    ! ================================
    ! Input and Output variables
    ! =================================
    INTEGER, INTENT(IN)                             :: nspec, n1, n2, nsol
    REAL (KIND=dp), DIMENSION(nsol),     INTENT(IN) :: solwav, solspec
    REAL (KIND=dp), DIMENSION(nspec), INTENT(INOUT) :: wav, spec
    INTEGER,                            INTENT(OUT) :: errstat

    ! ===============
    ! Local variables
    ! ===============
    REAL (KIND=dp)    :: delw, temp      
    INTEGER           :: i, fidx, lidx, nref

    INTEGER, EXTERNAL :: OMI_SMF_setmsg

    ! --------------------------------
    ! Name of this subroutine/module
    ! --------------------------------
    CHARACTER (LEN=14), PARAMETER :: modulename = 'append_solring'

    errstat = pge_errstat_ok

    IF (n1 > 1) THEN
      !delw = wav(n1 + 1) - wav(n1)
      delw = (wav(n2) - wav(n1))/(n2-n1) ! jbak, much better when there are large band pixels
      DO i = n1-1, 1, -1
        wav(i) = wav(i+1) - delw
      ENDDO
    ENDIF

    IF (n2 < nspec) THEN
      !delw = wav(n2) - wav(n2-1)
      delw = (wav(n2) - wav(n1))/(n2-n1) ! jbak, much better when there are large band pixels
      DO i = n2 + 1, nspec
        wav(i) = wav(i - 1) + delw
      ENDDO
    ENDIF

    IF (n1 > 1) THEN     
      ! Perform interpolatioon
      fidx = MINVAL(MAXLOC(solwav, MASK=(solwav < wav(1))))
      lidx = MINVAL(MINLOC(solwav, MASK=(solwav > wav(n1))))
      nref = lidx - fidx + 1
      IF (nref <= 4) THEN
        lidx = lidx + 2
        fidx = fidx - 2
      ENDIF
      IF (fidx < 1 .OR. lidx > nsol) THEN 
        WRITE(www_lun, *) modulename, ': Increase wavelength range of solar reference!!!'
        WRITE(www_lun, *)  fidx, lidx, nsol
        errstat = pge_errstat_error; RETURN
      ENDIF
      temp = spec(n1)

      CALL interpolation (nref, solwav(fidx:lidx), solspec(fidx:lidx), &
           n1,  wav(1:n1), spec(1:n1), errstat )

      IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); RETURN
      END IF

      spec(1:n1-1) = spec(1:n1-1) * temp / spec(n1)
    ENDIF

    IF (n2 < nspec) THEN  
      ! Perform interpolatioon
      fidx = MINVAL(MAXLOC(solwav, MASK=(solwav < wav(n2))))
      lidx = MINVAL(MINLOC(solwav, MASK=(solwav > wav(nspec))))
      nref = lidx - fidx + 1
      IF (nref <= 4) THEN
        lidx = lidx + 2; fidx = fidx -2
      ENDIF

      IF (fidx < 1 .OR. lidx > nsol) THEN 
        WRITE(www_lun, *) modulename, 'Increase wavelength range of solar reference!!!'
        WRITE(www_lun,*) fidx, nsol, solwav(1), wav(n2)
        errstat = pge_errstat_error; RETURN
      ENDIF

      temp = spec(n2)
      CALL interpolation (nref, solwav(fidx:lidx), solspec(fidx:lidx), &
           nspec-n2+1,  wav(n2:nspec), spec(n2:nspec), errstat )
      IF ( errstat > pge_errstat_warning ) THEN
        errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0); RETURN
      END IF
      spec(n2+1:nspec) = spec(n2+1:nspec) * temp / spec(n2)
    ENDIF

    RETURN

  END SUBROUTINE append_solring


end module m_prepare_databases

!   Unused?
!
!  SUBROUTINE CORRECT_COADDEFFECT(refwav, refspec, i0, nref, choice, nout, errstat)
!    USE OMSAO_precision_module     
!    USE OMSAO_errstat_module
!    IMPLICIT NONE
!
!    INTEGER,                          INTENT (IN)    :: nref, choice
!    INTEGER,                          INTENT (OUT)   :: errstat, nout
!    REAL (KIND=dp), DIMENSION (nref), INTENT (IN)    :: refwav, i0
!    REAL (KIND=dp), DIMENSION (nref), INTENT (INOUT) :: refspec
!
!    ! ---------------
!    ! Local variables
!    ! ---------------
!    !INTEGER                          :: i, j
!    REAL (KIND=dp), DIMENSION (nref) :: newi0, abspec
!    REAL (KIND=dp)                   :: scalex!, frefw, lrefw
!    !CHARACTER (LEN=19), PARAMETER    :: modulename = 'CORRECT_COADDEFFECT'
!
!    errstat = pge_errstat_ok
!    scalex = 0.08              !  Arbitarily assummed
!    newi0 = i0
!
!    IF (choice == 1 .OR. choice == 2) THEN 
!      abspec = i0 * EXP(-refspec * scalex)
!    ENDIF
!
!    IF (choice == 1) THEN 
!      CALL avg_band_refspec(refwav, abspec, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!      CALL avg_band_refspec(refwav, newi0, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!    ELSE IF (choice == 2) THEN     
!      CALL avg_band_ozcrs(refwav, abspec, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!      CALL avg_band_ozcrs(refwav, newi0, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!    ELSE IF (choice == 3) THEN
!      CALL avg_band_refspec(refwav, refspec, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!    ENDIF
!
!
!    IF (choice == 1 .OR. choice == 2) THEN
!      refspec(1:nout) = -LOG(abspec(1:nout) / newi0(1:nout)) / scalex
!    ENDIF
!
!    RETURN
!
!  END SUBROUTINE CORRECT_COADDEFFECT




!   Unused?
!   * here super_gauss is not implemented (jbak)
!  SUBROUTINE normalize_solar_refspec ( n_radwvl, curr_rad_wvl, solar_spec, errstat)
!
!    USE OMSAO_precision_module
!    USE OMSAO_indices_module,     ONLY: solar_idx, wvl_idx, spc_idx
!    USE OMSAO_variables_module,   ONLY: n_refspec_pts, refspec_orig_data, &
!         yn_varyslit, which_slit, refspec_norm, solar_refspec
!    USE OMSAO_parameters_module,  ONLY: max_spec_pts
!    USE ozprof_data_module,       ONLY: div_sun
!    USE OMSAO_slitfunction_module    
!    USE OMSAO_errstat_module
!    IMPLICIT NONE
!
!    INTEGER,                              INTENT (IN)  :: n_radwvl
!    REAL (KIND=dp), DIMENSION (n_radwvl), INTENT (IN)  :: curr_rad_wvl, solar_spec
!    INTEGER,                              INTENT (OUT) :: errstat
!
!    ! ---------------
!    ! Local variables
!    ! ---------------
!    INTEGER                                  :: fidx, lidx, ni0!, i
!    REAL (KIND=dp), DIMENSION (max_spec_pts) :: wave, specmod, ratio
!    REAL (KIND=dp), DIMENSION (n_radwvl)     :: solar_spec0, ratio0
!    REAL (KIND=dp)                           :: frefw, lrefw
!
!    ! ------------------
!    ! External functions
!    ! ------------------
!    !INTEGER :: OMI_SMF_setmsg
!
!    CHARACTER (LEN=23), PARAMETER :: modulename = 'normalize_solar_refspec'
!
!    errstat = pge_errstat_ok
!
!    ! Convole high-resolution solar reference spectrum
!    ni0  = n_refspec_pts(solar_idx)
!    wave = refspec_orig_data(solar_idx,1:ni0,wvl_idx)
!
!    IF (.NOT. yn_varyslit) THEN
!      IF (which_slit == 0) THEN
!        CALL gauss_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 1) THEN
!        CALL asym_gauss_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 2) THEN
!        CALL asym_voigt_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 3) THEN
!        CALL triangle_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE
!        CALL omislit_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      END IF
!    ELSE 
!      IF (which_slit == 0) THEN
!        CALL gauss_vary (wave(1:ni0), solar_refspec(1:ni0),specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 1) THEN
!        CALL asym_gauss_vary (wave(1:ni0), solar_refspec(1:ni0),specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 2) THEN
!        CALL asym_voigt_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 3) THEN
!        CALL triangle_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE
!        CALL omislit_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ENDIF
!    ENDIF
!
!    CALL bspline(wave(1:ni0), specmod(1:ni0), ni0, curr_rad_wvl, solar_spec0, n_radwvl, errstat)
!    ratio0 = solar_spec / solar_spec0 * div_sun / refspec_norm(solar_idx)
!
!    frefw = curr_rad_wvl(1)
!    lrefw = curr_rad_wvl(n_radwvl)
!    fidx = MINVAL(MINLOC(wave(1:ni0), MASK=(wave(1:ni0) >= frefw)))
!    lidx = MINVAL(MAXLOC(wave(1:ni0), MASK=(wave(1:ni0) <= lrefw)))
!
!    CALL bspline(curr_rad_wvl, ratio0, n_radwvl, wave(fidx:lidx), ratio(fidx:lidx), lidx-fidx + 1, errstat)
!    IF (errstat < 0) THEN
!      WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
!      errstat = pge_errstat_error
!    ENDIF
!    IF (fidx > 1) ratio(1:fidx-1) = ratio(fidx)
!    IF (fidx < ni0) ratio(lidx+1:ni0) = ratio(lidx)
!
!    refspec_orig_data(solar_idx, 1:ni0, spc_idx) = &
!         refspec_orig_data(solar_idx, 1:ni0, spc_idx) * ratio(1:ni0)
!
!    RETURN
!  END SUBROUTINE normalize_solar_refspec
