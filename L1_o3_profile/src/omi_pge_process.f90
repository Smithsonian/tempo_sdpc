!
MODULE omi_pge_process
  PUBLIC omi_fitting_process
  
CONTAINS

  SUBROUTINE omi_fitting_process ( message, pge_error_status )

    USE OMSAO_precision_module
    USE OMSAO_variables_module,  only: wavcal,  pixnum_lim, linenum_lim, &
         npix_fitting, npix_fitted, l2_filename, fitvar_rad_saved, &
         n_fitvar_rad, currpix, currline, currloop, which_slit, scnwrt, use_solcomp, use_backup, &
         mask_fitvar_rad, reduce_resolution, l2_cld_filename, &
         l1b_rad_filename,  numwin, radnhtrunc, &
         l1b_irrad_filename, nxbin, nybin, ncoadd,&
         nxtrack, ntimes, num_wav_max, ntimes_loop, offset_line, calwrt, the_pix, the_line
    USE OMSAO_parameters_module, ONLY: l2funit, lcurve_unit, ozwrtint_unit,calunit
    USE ozprof_data_module, only: lcurve_write, ozwrtint,  &
         lcurve_fname, ozwrtint_fname, &
         ozabs_convl, so2crs_convl, o2crs_convl, o4crs_convl,which_cld, allrms
    USE OMSAO_errstat_module
    USE ascii_output_module, only: write_final
    USE m_cross_calibrate, ONLY:calibrate_rad_cross, calibrate_irrad_cross
    USE m_specfit_ozprof
    use m_allocate
    USE O3P_output_module
    ! should be changed depending on instrument
    use OMSAO_pixelcorner_module
    use OMSAO_omicloud_module
    use OMSAO_slitfunction_module
    use OMSAO_omidata_module, only: nlines_max, &
         nfxtrack, zoom_mode, zoom_p1, zoom_p2,  &
         omi_cali, omi_irrad, omi_rad, omi_refl, omi_ring, omi_geo, omi_o3p
    USE omi_read_l1b_data, only: omi_read_irradiance_data, &
        omi_read_radiance_lines,  replace_solar_irradiance, omi_set_parameters
    USE omi_adj_data

    IMPLICIT NONE

    ! -----------------------
    ! Input/Output variables
    ! ----------------------
    integer, intent (OUT) :: pge_error_status
    character(len=100), INTENT(OUT) :: message
    ! -------------------------
    ! Local variables (for now)
    ! -------------------------
    integer :: i,first_line, last_line, iline, nxcoadd, first_pix, last_pix, &
         exval, initval, errstat, curr_fitted_line, sline, eline, ix, ext
    real (kind=dp), dimension(3)    :: fitcol
    real (kind=dp), dimension(3, 2) :: dfitcol
    real (kind=dp)     :: fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg,rms
    logical            :: reduce_resolution_save
    integer :: processing_version = 1
    !------------------------------------
    ! Error variables
    !------------------------------------
    LOGICAL :: problems=.false.
    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    character (len=23), parameter :: modulename = 'omi_fitting_process'

    pge_error_status = pge_errstat_ok
    errstat = 0
    !-------------------------------------------------------------------------------
    ! @ inital setup
    !-------------------------------------------------------------------------------
    call omi_set_parameters ( pge_error_status )
    IF ( pge_error_status >= pge_errstat_error ) THEN
        message = ':Failed in omi_set_parameters'
        return
    ENDIF

    ! xtrack positions to be processed (start from the actual binned position)
    IF (linenum_lim(2) >= ntimes)  linenum_lim(2) = ntimes
    IF (zoom_mode) then
      pixnum_lim(1) = max(nint(1.0 * (zoom_p1 + 1) / ncoadd), pixnum_lim(1))
      pixnum_lim(2) = min(zoom_p2 / ncoadd, pixnum_lim(2))
      if ( mod(pixnum_lim(2) - pixnum_lim(1) + 1, nxbin) /= 0) then
        write(www_lun, '(A,2I4)') 'Incorrect across track binning option: ', &
             pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; return
      endif

      if (pixnum_lim(1) > pixnum_lim(2)) then
        write(www_lun, *) 'This is a zoom in mode orbit!!!'
        write(www_lun, '(A, I2, A, I2)') &
             'Invalid pixel selection, must be between ', &
             (zoom_p1 / ncoadd) * ncoadd + 1, ' and ', &
             (zoom_p2 / ncoadd ) * ncoadd
        pge_error_status = pge_errstat_error; return
      endif
    endif

    first_pix  = ceiling(1.0 * pixnum_lim(1) / nxbin)
    last_pix  = nint(1.0 * pixnum_lim(2) / nxbin )
    nxcoadd = nxbin * ncoadd

    ! Line number, starting from zero and keep track of offset
    offset_line = linenum_lim(1) - 1; first_line = 1
    last_line   = int ((linenum_lim(2) - linenum_lim(1) + 1.0) / nybin)
    IF (last_line > 100) then
      WRITE( message, *)  ': number of binned lines (>100) = ', last_line
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Allocate some memory - must be before cross_calibrate calls
    !-------------------------------------------------------------------
    call allocate_geo(nxtrack, ntimes, omi_geo, pge_error_status)
    CALL allocate_spec (numwin, nxtrack,nlines_max, omi_irrad, omi_rad, omi_ring, omi_refl, omi_cali, pge_error_status)

    IF ( pge_error_status >= pge_errstat_error ) THEN
        message = ':Failed to allocate variables'
        return
    ENDIF

    allocate(omi_o3p%fitvar(nxtrack, 0:nlines_max-1, n_fitvar_rad))
    allocate(omi_o3p%initval(nxtrack, 0:nlines_max-1))
    allocate(omi_o3p%exitval(nxtrack, 0:nlines_max-1))
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Read irradiance & slit/wavelength calibration
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !------------------------------------------------------------------------    
    ! load instrument slit parameters
    !-------------------------------------------------------------------------
    IF (which_slit == 5) THEN
      call load_slitpars (pge_error_status)     
      IF ( pge_error_status >= pge_errstat_error ) THEN 
        message = ':Failed to read instrument slit function'
        RETURN
      ENDIF
      if (scnwrt) write(*, '(A)') 'Finish loading OMI slit parameters !!!'
    ENDIF
   
    reduce_resolution_save = reduce_resolution; reduce_resolution = .false.
    CALL omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, &
         pge_error_status )

    IF ( pge_error_status >= pge_errstat_error ) then
      use_backup = .true.
      CALL omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, &
           pge_error_status )
      IF ( pge_error_status >= pge_errstat_error ) return
    ENDIF
    IF (scnwrt) write(*, '(A)') 'Finish reading irradiances!!!'

    IF (use_solcomp) THEN
      print * ,'not implemented'
      STOP
      call replace_solar_irradiance(calunit, first_pix, last_pix, &
           pge_error_status)
      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A)') &
           'Finish replacing irradiances with solar composite!!!'
    ENDIF    
    ! Calibrate OMI irradiances for all cross track pixels and save results
    IF (reduce_resolution_save) reduce_resolution = .true.
    CALL calibrate_irrad_cross (omi_irrad, omi_ring, omi_cali, first_pix, last_pix, pge_error_status)
    IF ( pge_error_status >= pge_errstat_error ) return
    IF (scnwrt) write(*, '(A)') 'Finish calibrating irradiances!!!'
    
    IF (reduce_resolution_save) THEN
      CALL omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, &
           pge_error_status )
      IF ( pge_error_status >= pge_errstat_error ) return
      IF (scnwrt) write(*, '(A)') &
           'Finish reading irradiances (with reduced resolution)!!!'
    ENDIF
    !-----------------------------------------------------------------------
    ! reading geolocation
    !----------------------------------------------------------------------
    !FIXME
    ! At present, compute_pixel_corners sets values in arrays whose size
    ! is determined by nlines_max = 100. This suggests that it expects
    ! to operate within one OMI data black of 100 lines. When operated
    ! in it's current position, outside the OMIBlock loop, if last_line >100
    ! the array indices will exceed nlines_max, causing memory to be
    ! overwritten and causing the program to crash / give false results
    ! Consider moving compute_pixel_corners inside the OMIBlock loop.

    ! Compute spatial pixel corners and effective viewing geometry
    CALL compute_pixel_corners ( omi_geo, ntimes, nfxtrack, last_line, errstat)
    IF (errstat /= 0) then
      message=': failed to read geolocation data !!!'
      pge_error_status = pge_errstat_error
    ENDIF
    IF (scnwrt) write(*, '(A)') 'Finish reading geolocation data!!!'
    !-----------------------------------------------------------------------
    ! reading cloud product
    !----------------------------------------------------------------------
    IF (which_cld == 4 .or. which_cld == 1 ) THEN 
      CALL read_omicldo2_clouds (ntimes, nxtrack, last_line, errstat)
    ELSE IF (which_cld == 3 .or. which_cld == 0 ) THEN 
      CALL read_omicldrr_clouds (ntimes, nxtrack, last_line, errstat) 
    ENDIF
    IF (errstat /= 0) then
      message=': failed to read cloud product !!!'
      pge_error_status = pge_errstat_error
    ENDIF
    IF (scnwrt) write(*, '(A, i3)') 'Finish reading L2 clouds!!!', which_cld
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Open output files
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    num_wav_max = maxval(omi_irrad%nwav(first_pix:last_pix)) - &
                  numwin * 2 * radnhtrunc
    IF (lcurve_write) THEN
      OPEN(UNIT=lcurve_unit, FILE=TRIM(ADJUSTL(lcurve_fname)), &
           STATUS='unknown', IOSTAT=errstat)
      IF ( errstat /= pge_errstat_ok ) THEN
        message=': failed to open lcurve file!!!'
        pge_error_status = pge_errstat_error
      ENDIF
    ENDIF
  
    IF (ozwrtint) THEN
      OPEN(UNIT=ozwrtint_unit, FILE=TRIM(ADJUSTL(ozwrtint_fname)), &
           STATUS='unknown', IOSTAT=errstat)
       IF ( errstat /= pge_errstat_ok ) THEN
        message=': failed to open lcurve file!!!'
        pge_error_status = pge_errstat_error
      ENDIF
    ENDIF

    CALL L2_O3P_create (ntimes,first_pix, last_pix, first_line, last_line, errstat)
    IF (errstat /= 0) then
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF
    IF (scnwrt) WRITE(*,'(A)') '@ Create '//TRIM(ADJUSTL(l2_filename)) 

    CALL L2_O3P_WRITE_GEO (omi_geo,first_pix, last_pix, first_line, last_line, errstat)
     IF (errstat /= 0) then
      pge_error_status = pge_errstat_error
      RETURN
    ENDIF 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! loop through each OMI data block
    ! 1. read each block
    ! 2. perform calibration for each block (middle line)
    ! 3. perform retrievals from first_pix to last_pix (all lines in block)    
    ! these pixels with exitval >= 0 will be used by subsequent retrievals
    !  as initial values
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Initialize fitting statistics
    npix_fitting = 0       ! number of pixels (failure + success)
    npix_fitted  = 0       ! number of successfully fitted pixels
    fitcol_avg   = 0.0
    rms_avg = 0.0
    dfitcol_avg = 0.0
    drel_fitcol_avg = 0.0
    omi_o3p%exitval = -10 ! no retrievals yet
    omi_o3p%fitvar  = 0.0

    OMIBlock: DO iline = 0, last_line-1, nlines_max

      ! Actually lines in OMI data
      ntimes_loop = nlines_max
      if ( iline + ntimes_loop > last_line ) ntimes_loop = last_line - iline
      sline = offset_line + iline * nybin
      eline = sline + ntimes_loop * nybin - 1
 
      ! Get NTIMES_LOOP radiance lines (with effective viewing geometry)
      call omi_read_radiance_lines (iline, ntimes_loop, sline, &
           first_pix, last_pix, nxcoadd, pge_error_status)
      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A,I4,A,I4)') &
           'Finishing reading radiances for lines: ', sline + 1, ' - ', eline+1

      ! Perform calibration for the middle scan line and apply to the other
      !  scan lines
      if (wavcal) then
        call calibrate_rad_cross (omi_rad, omi_cali, first_pix,last_pix,first_line, last_line, pge_error_status)
        if ( pge_error_status >= pge_errstat_error ) return
        if (scnwrt) write(*, '(A,I4,A,I4)') &
             'Finishing calibrating radiances for lines: ', sline+1, ' - ', &
             eline+1
      endif

      ! loop through each xtrack position
      ! 1. Prepare databases, adjust radiances
      ! 2. Process all pixels at this poistion
      XtrackPix: do currpix = first_pix, last_pix
        !! Kai
        if( scnwrt ) write(*,'(A,I5,A,I3)') 'Doing Line=',iline,', iX=',currpix
        ! Need to convolve high-resolution ozone absorption cross section
        !  (for this position). Once the xsection is convolved, it will be
        !  set to false in ROUTINE getabs_crs
        ozabs_convl = .true.; so2crs_convl = .true. ; o4crs_convl = .true.
        IF (ALL(omi_rad%pix_errstat(currpix, 0:ntimes_loop-1) == pge_errstat_error) &
             .OR. omi_irrad%errstat(currpix) == pge_errstat_error) THEN
           omi_o3p%exitval(currpix, 0:ntimes_loop-1) = -10; CYCLE
        ENDIF
        ! Load/adjust irradiances and slit calibration parameters
        call adj_solar_data (pge_error_status)
        if ( pge_error_status >= pge_errstat_error ) cycle
        curr_fitted_line = 0
        YfitLine: do currloop = 0, ntimes_loop - 1

          omi_o3p%exitval(currpix, currloop) = -10
          currline = iline + currloop
          the_line= currline * nybin + offset_line + 1
          the_pix = (currpix-1) * nxbin  + 1
          if (omi_rad%pix_errstat(currpix, currloop) == pge_errstat_error .or. &
               omi_irrad%errstat(currpix) == pge_errstat_error ) &
               omi_o3p%exitval(currpix, currloop) = -9

          ! Load/adjust radiances/geolocations fields for a particular pixel
          ! Prepare databases for the first pixel (ifitline == 1)
          if (omi_o3p%exitval(currpix, currloop) == -10) then
            call adj_earthshine_data (curr_fitted_line, pge_error_status) 
            if ( pge_error_status >= pge_errstat_error ) &
                 omi_o3p%exitval(currpix, currloop) = -9
          endif
        
          if (omi_o3p%exitval(currpix, currloop) == -10) then
            if (scnwrt) write(*, '(A,I5,A10,I5, A10, I5)') &
                 'OMI Pixel: Line = ',the_line, ' XPix = ', the_pix, &
                 ' Loop = ', currloop

            initval = omi_o3p%initval(currpix, currloop)
            call specfit_ozprof (initval, fitcol, dfitcol,rms, exval)
            ! Store exit status for current pixel
            omi_o3p%exitval(currpix, currloop) = exval
            omi_o3p%fitvar(currpix, currloop, 1:n_fitvar_rad) &
                 = fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad))
          else
            exval = -9
          endif

          ! Write retrievals        
          CALL L2_O3P_write_data (currpix, first_pix, last_pix, currloop, currline, ntimes_loop,&
                      exval, fitcol, dfitcol, message, problems)

          IF (exval >= 0 .and. fitcol(1) > 0.0 .and. dfitcol(1, 1) >= 0.0 ) THEN
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
      ENDDO XtrackPix ! do currpix = first_pix, last_pix
    ENDDO OMIBlock ! do iline = 0, last_line-1, nlines_max

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Final 
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


    !-----------------------------------------------------------------
    ! Deallocate any remaining arrays
    !----------------------------------------------------------------
    call dealloc (errstat)

    !----------------------------------------------------------------
    ! Close L2 output file
    !----------------------------------------------------------------
    IF (lcurve_write) CLOSE (lcurve_unit)
    IF (ozwrtint)     CLOSE (ozwrtint_unit)
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

    RETURN
  END SUBROUTINE omi_fitting_process
END MODULE omi_pge_process
