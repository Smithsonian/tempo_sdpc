!
MODULE tmpo_pge_process

  PUBLIC tmpo_fitting_process

CONTAINS
  SUBROUTINE tmpo_fitting_process  ( message, pge_error_status )
    USE OMSAO_precision_module
    USE OMSAO_indices_module, ONLY: max_calfit_idx, n_max_fitpars
    USE OMSAO_variables_module,  only:num_wav_max, num_param, & 
         nxtrack, ntimes, pixnum_lim, linenum_lim, &
         nxbin, nybin, currpix, currline, currloop,&
         n_fitvar_rad, mask_fitvar_rad,fitvar_rad_saved,& 
         npix_fitting, npix_fitted, numwin, radnhtrunc,&
         scnwrt, use_backup, use_solcomp,reduce_resolution, wavcal, which_slit, &
         l1b_irrad_filename, l1b_rad_filename,l2_cld_filename,l2_filename,  &
         l2_hdf_flag, pix_pos, line_pos, ntimes_loop, offset_line,&
         the_pix, the_line, calwrt

    USE OMSAO_parameters_module, ONLY: l2funit, lcurve_unit, ozwrtint_unit,calunit
    USE ozprof_data_module, only: lcurve_write, ozwrtint,&
         lcurve_fname, ozwrtint_fname,&
         ozabs_convl, so2crs_convl, o2crs_convl, o4crs_convl, h2ocrs_convl,allrms
    USE OMSAO_tmpodata_module, only: rad_swathname, & 
        tmpo_rad, tmpo_irrad,tmpo_refl, tmpo_ring, tmpo_cali, & 
        tmpo_geo1,tmpo_geo2,tmpo_o3p,nlines_max
    USE OMSAO_errstat_module
    USE m_specfit_ozprof
    USE m_allocate
    USE o3p_output_module
    USE m_read_geo_tio
    USE m_read_cloud_tio
    USE m_cross_calibrate, ONLY: calibrate_irrad_cross,calibrate_rad_cross
    USE tmpo_read_l1b_data
    USE tmpo_adj_data, ONLY:adj_solar_data, adj_earthshine_data
    USE tio_output_module
    USE ascii_output_module, only: write_final
    implicit none

    ! -----------------------
    ! Input/Output variables
    ! ----------------------
    INTEGER, INTENT (OUT) :: pge_error_status
    CHARACTER(LEN=100), INTENT(OUT) :: message
    ! -------------------------
    ! Local variables (for now)
    ! -------------------------
    LOGICAL :: reduce_resolution_save
    INTEGER :: exval, initval, errstat, curr_fitted_line,ix, iline
    INTEGER :: first_line, last_line, first_pix, last_pix
    INTEGER :: npix, nline ! actual number of pixel/line 
    INTEGER :: sline, eline,spix, epix ! actual position in L1 domain (both staring 1)
    REAL (kind=dp) :: rms, fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg
    REAL (kind=dp),DIMENSION(3)    :: fitcol
    REAL (kind=dp),DIMENSION(3, 2) :: dfitcol
    ! FIXME - should be input variable, not fixed value
    INTEGER :: processing_version = 1
    LOGICAL :: problems = .false.
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (len=23), parameter :: modulename = 'tmpo_fitting_process'

    !----------------------------------------------------------------------------
    ! @ Initial Setup
    !-----------------------------------------------------------------------------
    pge_error_status = pge_errstat_ok
    errstat = 0
    CALL tmpo_set_parameters (pge_error_status)
    ! insch, the_year, the_month_the_day, jday is set-up
    ! ntimes, nxtrack, nwavel, granuelyear/month/day
    IF ( pge_error_status >= pge_errstat_error ) THEN
        message = ':Failed in tmpo_set_parameters'
        RETURN
    ENDIF

    !  define the boundaries of along and across track domain 
    ! linenum_lim and pixnum_lim is actual location in TEMPO domain
    ! linelim and pixlim is location for current mpi process
    IF (linenum_lim(2) >= ntimes)  linenum_lim(2) = ntimes 
    IF (pixnum_lim(2) >= nxtrack)  pixnum_lim(2) = nxtrack
   
    first_pix  = ceiling(1.0 * pixnum_lim(1) / nxbin)
    last_pix   = nint(1.0 * pixnum_lim(2) / nxbin )
    offset_line = linenum_lim(1) - 1
    first_line  = 1
    last_line   = int ((linenum_lim(2) - linenum_lim(1) + 1.0) / nybin)

    npix  = last_pix - first_pix + 1
    nline = last_line - first_line + 1
    offset_line = linenum_lim(1) - 1

    spix = (first_pix-1)*nxbin +1 
    epix = (last_pix-1)*nxbin  +1
    sline = offset_line + 1
    eline = offset_line + (nline-1)*nybin + 1
    allocate (pix_pos(nxtrack), line_pos(ntimes))
    DO ix = first_pix, last_pix
       pix_pos(ix) = (ix-1)*nxbin + 1
    ENDDO
    DO iline = first_line, last_line
       line_pos(iline) = (iline-1)*nybin +1 + offset_line 
    ENDDO
    WRITE (*,*) '@ Define TEMPO Domain**'
    WRITE(*,'(A,2i5, A,i2)') '=>pixnum_lim :', pixnum_lim, "nxbin:", nxbin
    WRITE(*,'(A,2i5, A,i2)') '=>linenum lim:', linenum_lim, "nxbin:", nybin
    WRITE(*,'(A,2i5,A,2I5)') '=>first/last pix :',first_pix,  last_pix, 'In',spix,epix
    WRITE(*,'(A,2i5,A,2I5)') '=>first/last line:',first_line, last_line,'In',sline, eline
    !print * , pix_pos(first_pix), pix_pos(last_pix)
    !print * , line_pos(first_line), line_pos(first_line)

    !IF (nline  > nlines_max) THEN
    !  pge_error_status = pge_errstat_error
    !  message=": number of binned lines defined in PCF exceeds maximum of 100"
    !  RETURN
    !ENDIF
    !-----------------------------------------------------------------
    ! Allocate some memory - must be before cross_calibrate calls
    !-------------------------------------------------------------------
   
    CALL allocate_spec (numwin, nxtrack,nlines_max, & 
            tmpo_irrad, tmpo_rad, tmpo_ring, tmpo_refl, tmpo_cali, pge_error_status)
    IF ( pge_error_status >= pge_errstat_error ) THEN
        message = ':Failed to allocate variables'
        return
    ENDIF

    CALL allocate_geo (nxtrack, ntimes, tmpo_geo1, pge_error_status)
    call allocate_geo (nxtrack, ntimes, tmpo_geo2, pge_error_status)

    IF (pge_error_status /= pge_errstat_ok) THEN 
       message =": failed to allocate geolocation variables"
       go to 111
    ENDIF

    ALLOCATE(tmpo_o3p%exitval(nxtrack,0:nlines_max))
    ALLOCATE(tmpo_o3p%initval(nxtrack,0:nlines_max))
    ALLOCATE(tmpo_o3p%fitvar(nxtrack, 0:nlines_max, n_fitvar_rad))



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Read irradiance & slit/wavelength calibration
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    !------------------------------------------------------------------------    
    ! load instrument slit parameters
    !-------------------------------------------------------------------------
    IF (which_slit == 5) THEN
      PRINT *, 'Please update slit function for real TEMPO algorithm'
      !call load_slitpars (pge_error_status)
      IF (pge_error_status /= pge_errstat_ok) THEN 
        message=": failed to read instrument slit function"
        RETURN
      ENDIF
      if (scnwrt) write(*, '(A)') '@ Finish loading slit parameters !!!'
    ENDIF

    call tmpo_read_irradiance (first_pix, last_pix, pge_error_status)
    IF (calwrt) close(calunit)
    reduce_resolution_save = reduce_resolution; reduce_resolution = .false.         
    IF ( pge_error_status >= pge_errstat_error ) then
      use_backup = .true.
      call tmpo_read_irradiance (first_pix, last_pix, pge_error_status )
      IF ( pge_error_status >= pge_errstat_error ) THEN
         message =": failed to read irradiances"
         RETURN
      ENDIF
    endif

    if (scnwrt) write(*, '(A)') '@ Finish reading irradiances!!!'
    ! if (use_solcomp) then
    !   call replace_solar_irradiance(calunit, first_pix, last_pix, &
    !        pge_error_status)
    !   if ( pge_error_status >= pge_errstat_error ) return
    !   if (scnwrt) write(*, '(A)') &
    !        'Finish replacing irradiances with solar composite!!!'
    ! endif
    IF (reduce_resolution_save) reduce_resolution = .true.
    call calibrate_irrad_cross (tmpo_irrad, tmpo_ring,tmpo_cali, first_pix, last_pix, pge_error_status)
    IF (pge_error_status /= pge_errstat_ok) THEN 
       message =": failed to calibrate irradiance"
       go to 111
    ENDIF
    if (scnwrt) write(*, '(A)') '@ Finish calibrating irradiances!!!'

    IF (reduce_resolution_save) THEN
      !call omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, &
      !     pge_error_status )
      !if ( pge_error_status >= pge_errstat_error ) return
      !if (scnwrt) write(*, '(A)') &
      !     '@ Finish reading irradiances (with reduced resolution)!!!'
    ENDIF
    
    !----------------------------------------------------------------
    ! @ Reading geolocation variables
    !   with Computing spatial pixel corners and effective viewing geometry
    !------------------------------------------------------------------
    !FIXME
    ! At present, compute_pixel_corners sets values in arrays whose size
    ! is determined by nlines_max = 100. This suggests that it expects
    ! to operate within one OMI data black of 100 lines. When operated
    ! in it's current position, outside the OMIBlock loop, if last_line >100
    ! the array indices will exceed nlines_max, causing memory to be
    ! overwritten and causing the program to crash / give false results
    ! Consider moving compute_pixel_corners inside the OMIBlock loop.

    CALL read_geo_tio (rad_swathname(1),tmpo_geo1, ntimes, nxtrack, & 
        first_pix, last_pix, sline, eline, .false.,pge_error_status)
    CALL read_geo_tio (rad_swathname(2),tmpo_geo2, ntimes, nxtrack, & 
        first_pix, last_pix, sline, eline, .true., pge_error_status)
    IF (pge_error_status /= pge_errstat_ok) THEN
       message =": failed to read geo location"
       RETURN
    ENDIF
    IF (scnwrt) write(*, '(A)') '@ Finish reading geolocation data!!!'
    !-----------------------------------------------------------------------
    ! reading cloud product
    !----------------------------------------------------------------------
    CALL read_cloud_tio(ntimes, nxtrack,sline, eline,pge_error_status)
    IF (pge_error_status /= pge_errstat_ok) THEN
       message =": failed to read cloud"
       RETURN
    ENDIF
    tmpo_geo1%cfr(1:nxtrack,0:nline-1) = L2_cloud%cfr(1:nxtrack, 0:nline-1)
    tmpo_geo1%ctp(1:nxtrack,0:nline-1) = L2_cloud%ctp(1:nxtrack, 0:nline-1)
    tmpo_geo1%cloud_qflg(1:nxtrack, 0:nline-1) = L2_cloud%qflags(1:nxtrack,0:nline-1)
    if (scnwrt) write(*, '(A)') '@ Finish reading L2 clouds !!!'

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Open output files
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    num_wav_max = maxval(tmpo_irrad%nwav(first_pix:last_pix)) - &
             numwin * 2 * radnhtrunc
    IF (lcurve_write) THEN
      OPEN(UNIT=lcurve_unit, FILE=TRIM(ADJUSTL(lcurve_fname)), &
           STATUS='unknown', IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        message=': failed to open lcurve file!!!'
        pge_error_status = pge_errstat_error
        go to 111
      ENDIF
    ENDIF
  
    IF (ozwrtint) THEN
      OPEN(UNIT=ozwrtint_unit, FILE=TRIM(ADJUSTL(ozwrtint_fname)), &
           STATUS='unknown', IOSTAT=errstat)
       IF ( errstat /= pge_errstat_ok ) THEN
        message=': failed to open lcurve file!!!'
        pge_error_status = pge_errstat_error
        go to 111
      ENDIF
    ENDIF

    CALL L2_O3P_create (ntimes,first_pix, last_pix, first_line, last_line, errstat)
    IF (errstat /= 0) then
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF
    IF (scnwrt) WRITE(*,'(A)') '@ Create '//TRIM(ADJUSTL(l2_filename)) 
    CALL L2_O3P_WRITE_GEO (tmpo_geo1,first_pix, last_pix, first_line, last_line, errstat)
     IF (errstat /= 0) then
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    ! loop through each OMI data block
    ! 1. read each block
    ! 2. perform calibration for each block (middle line)
    ! 3. perform retrievals from first_pix to last_pix (all lines in block)     
    ! these pixels with exitval >= 0 will be used by subsequent retrievals
    !  as initial values
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Initialize fitting statistics
    npix_fitting = 0       ! number of pixels (failure + success)
    npix_fitted  = 0       ! number of successfully fitted pixels
    fitcol_avg   = 0.0
    rms_avg = 0.0
    dfitcol_avg = 0.0
    drel_fitcol_avg = 0.0
    tmpo_o3p%exitval = 10
    tmpo_o3p%fitvar = 0.0
    tmpo_o3p%initval = 0
    OMIBlock: do iline =  0, last_line-1, nlines_max
      
       ! Actually lines in L1B Data
      ntimes_loop = nlines_max
      if ( iline + ntimes_loop > last_line ) ntimes_loop = last_line - iline
      sline = offset_line + iline * nybin +1
      eline = sline + ntimes_loop * nybin -1
      call tmpo_read_radiance_lines (iline, first_pix, last_pix, sline, eline, pge_error_status)

      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A,I4,A,I4)') &
           '@ Finishing reading radiances for lines: ', sline, ' - ', eline
      ! Perform calibration for the middle scan line and apply to the other
      !  scan lines
      if (wavcal) then
        call calibrate_rad_cross (tmpo_rad, tmpo_cali,first_pix,last_pix,first_line, last_line, pge_error_status)
        if ( pge_error_status >= pge_errstat_error ) return
        if (scnwrt) write(*, '(A,I4,A,I4)') &
           '@ Finishing calibrating radiances for lines: ', sline, ' - ', eline
      endif
      if (calwrt) cycle
      ! loop through each xtrack position
      ! 1. Prepare databases, adjust radiances
      ! 2. Process all pixels at this poistion
      XtrackPix: do currpix = last_pix, first_pix, - 10        
        ! Need to convolve high-resolution ozone absorption cross section
        !  (for this position). Once the xsection is convolved, it will be
        !  set to false in ROUTINE getabs_crs
        ozabs_convl = .true.; so2crs_convl = .true. ; o4crs_convl = .true. 
        o2crs_convl = .true.;h2ocrs_convl = .true.

        IF (ALL(tmpo_rad%pix_errstat(currpix, 0:ntimes_loop-1) == pge_errstat_error) &
             .OR. tmpo_irrad%errstat(currpix) == pge_errstat_error) THEN
           tmpo_o3p%exitval(currpix, 0:ntimes_loop-1) = -10; CYCLE
        ENDIF
        ! Load/adjust irradiances and slit calibration parameters
        call adj_solar_data (pge_error_status) 
        if ( pge_error_status >= pge_errstat_error ) cycle

        curr_fitted_line = 0
        YfitLine: do currloop = 0, ntimes_loop - 1

          tmpo_o3p%exitval(currpix, currloop) = -10
          currline = iline + currloop
          the_line= currline * nybin + offset_line + 1
          the_pix = (currpix-1) * nxbin  + 1
          if (tmpo_rad%pix_errstat(currpix, currloop) == pge_errstat_error .or. &
               tmpo_irrad%errstat(currpix) == pge_errstat_error ) &
               tmpo_o3p%exitval(currpix, currloop) = -9

          ! Load/adjust radiances/geolocations fields for a particular pixel
          ! Prepare databases for the first pixel (ifitline == 1)
          if (tmpo_o3p%exitval(currpix, currloop) == -10) then
            call adj_earthshine_data (curr_fitted_line, pge_error_status) 
            if ( pge_error_status >= pge_errstat_error ) &
                 tmpo_o3p%exitval(currpix, currloop) = -9
          endif
        
          if (tmpo_o3p%exitval(currpix, currloop) == -10) then

            initval = tmpo_o3p%initval(currpix, currloop)
            initval = 0
            call specfit_ozprof (initval, fitcol, dfitcol,rms,  exval)
            ! Store exit status for current pixel
            tmpo_o3p%exitval(currpix, currloop) = exval
            tmpo_o3p%fitvar(currpix, currloop, 1:n_fitvar_rad) &
                 = fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad))
        
            !if (scnwrt) write(*, '(A,2I5,A,I5, A, I4, A, i3)') &
            ! '@ O3P Retrieval: Line =', the_line, iline, ' XPix= ', the_pix,' init =',initval, 'exval=', exval
          else
            exval = -9
          endif

         ! Write retrievals        
          CALL L2_O3P_write_data (currpix, first_pix, last_pix, currloop, currline, ntimes_loop,&
                      exval, fitcol, dfitcol, message, problems)
          IF (problems) THEN
           print *, message; stop
          ENDIF
          IF (exval >= 0 .and. fitcol(1) > 0.0 .and.dfitcol(1, 1) >= 0.0 ) THEN
            ! -----------------------------------------------------------------
            ! Some general statistics on the average fitted column and
            ! uncertainty. Again, we make sure that only "good" fits are
            ! included in the average.
            ! ----------------------------------------------------------------
            fitcol_avg       = fitcol_avg + fitcol(1)
            rms_avg          = rms_avg + rms
            dfitcol_avg      = dfitcol_avg + dfitcol(1, 1)
            drel_fitcol_avg  = drel_fitcol_avg + dfitcol(1, 1) / fitcol(1)
            npix_fitted      = npix_fitted + 1
            curr_fitted_line = curr_fitted_line + 1
          ENDIF
          npix_fitting        = npix_fitting + 1
        ENDDO YfitLine
      ENDDO XtrackPix
    ENDDO OMIBlock ! from first_line to last_line
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Final output
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IF (calwrt) THEN
       close(calunit)
       STOP
    ENDIF
    IF (scnwrt) THEN
      CALL write_final(fitcol_avg, rms_avg, dfitcol_avg,drel_fitcol_avg, &
           npix_fitted)
      WRITE(*, '(2(A,I5))') 'Number of pixels = ', &
      npix_fitting, '   Number of fitted pixels = ', npix_fitted
    ENDIF

    !----------------------------------------------------------------
    ! close l1b radiance file       
    !----------------------------------------------------------------
    ! If using netCDF inputs, copy critical metadata and label file
    IF (l2_hdf_flag == 4) THEN
      call copy_hdr_metadata (l1b_rad_filename, errstat)
      call label_output_file ("o3p", processing_version, errstat)
    ENDIF

    !----------------------------------------------------------------
    ! Close L2 output file
    !----------------------------------------------------------------
    IF (lcurve_write) CLOSE (lcurve_unit)
    IF (ozwrtint)     CLOSE (ozwrtint_unit)
    
     CALL L2_O3P_close (errstat)
    IF (errstat /=0) THEN
      pge_error_status = pge_errstat_error
      message =": failed to close "// ADJUSTL(TRIM((l2_filename)))
      RETURN
    ENDIF


    !-----------------------------------------------------------------
    ! Deallocate any remaining arrays
    !----------------------------------------------------------------
    call dealloc (errstat)
    RETURN
    !------------------------------------------------------------
    ! Error Handling
    !------------------------------------------------------------
    111 CONTINUE
    message = ADJUSTL(TRIM(modulename))//ADJUSTL(TRIM(message))
    WRITE(*,*) message
    call tell_error (tell_io_write_error, message, pge_error_status)
    STOP
  END SUBROUTINE tmpo_fitting_process

END MODULE tmpo_pge_process
