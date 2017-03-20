!
module m_omi_fitting_process

  public omi_fitting_process
  private

contains

  subroutine omi_fitting_process ( l2_hdf_flag,  pge_error_status )

    use OMSAO_precision_module
    use OMSAO_variables_module,  only: wavcal,  pixnum_lim, linenum_lim, &
         npix_fitting, npix_fitted, l2_filename, fitvar_rad_saved, &
         n_fitvar_rad, currpix, currline, currloop, ozabs_convl, &
         so2crs_convl, which_slit, scnwrt, use_solcomp, use_backup, &
         mask_fitvar_rad, reduce_resolution, l2_cld_filename, &!, coadd_uv2
         use_he5_in, use_he5_out, use_tio_in, use_tio_out, l1b_rad_filename, &
         nc_rad_swathname, l1_rad_filename_nc, numwin, radnhtrunc, &
         l1_irrad_filename_nc, l1b_irrad_filename, &
         GranuleDay, GranuleMonth, GranuleYear!, GranuleJDay
    use ozprof_data_module, only: lcurve_write, ozwrtint, l2funit, &
         lcurve_fname, ozwrtint_fname, lcurve_unit, ozwrtint_unit, &
         algorithm_name, algorithm_version, calunit, nfgas, nlay, &
         ozfit_start_index, ozfit_end_index!, which_cld
    use OMSAO_pixelcorner_module
    use OMSAO_omicloud_module
    use OMSAO_slitfunction_module
    use OMSAO_omidata_module, only: nlines_max, ntimes, ntimes_loop, &
         nxtrack, nfxtrack, ncoadd, omi_radpix_errstat, omi_exitval, &
         omi_fitvar, omi_initval, omi_solpix_errstat, nxbin, nybin, &
         omi_irradiance_wavl, &
         offset_line, zoom_mode, zoom_p1, zoom_p2, omi_nwav_irrad
    use he5_output_module, only: he5_l2setgeofields, he5_l2setdatafields, &
         he5_l2wrtinit
    use OMSAO_errstat_module
    use OMI_metaData_class
    use PROFOZ_metaDef
    use OMI_LUN_set
    use omi_fitting_aux, only: omi_adj_earthshine_data, omi_adj_solar_data, &
         omi_set_fitting_parameters, timestamp
    use omi_cross_calibrate, only: omi_rad_cross_calibrate, &
         omi_irrad_cross_calibrate
    use omi_read_l1b_data, only: omi_read_irradiance_data, &
         omi_read_radiance_lines, omi_read_radiance_paras, &
         replace_solar_irradiance
    use m_write_final, only: omi_write_intermed, write_final
    use m_specfit_ozprof
    use tio_output_module
    use m_read_cloud_tio
    use m_read_l1_tio, only: read_l1_dims_tio
    use m_read_geo_tio
    use m_read_metadata_tio, only: read_date_tio



    implicit none

    ! -----------------------
    ! Input/Output variables
    ! ----------------------
    integer, intent (IN)  :: l2_hdf_flag
    integer, intent (OUT) :: pge_error_status

    ! -------------------------
    ! Local variables (for now)
    ! -------------------------
    integer :: first_line, last_line, iline, nxcoadd, first_pix, last_pix, &
         exval, initval, errstat, curr_fitted_line, sline, eline, ix, &!, i
         num_param, num_wav_max
    real (kind=dp), dimension(3)    :: fitcol
    real (kind=dp), dimension(3, 2) :: dfitcol
    real (kind=dp)     :: fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, rms
    character (len=24) :: currtime
    logical            :: reduce_resolution_save

    integer :: version, year, month, day, jday
    type (OMIECSMETA_T) :: L1BcoreMeta
    type (ECSMETA_ITEM_T), dimension(6) :: PROFOZ_metaItems 

    integer (kind=4), dimension(3) :: LUNinputPointer
    !CHARACTER(len=7) :: cldtype
    character(len=6) :: ShortName = 'PROFOZ'
    integer :: ext
    character(len=1024) :: nc_l2_filename, l2_cld_filename_nc

    ! FIXME - should be input variable, not fixed value
    integer :: processing_version = 1

    !INTEGER, parameter :: ntemp = 18001
    !REAL (kind=dp), DIMENSION(ntemp) :: waves, raycofs, depols

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    character (len=23), parameter :: modulename = 'omi_fitting_process'

    pge_error_status = pge_errstat_ok

    ! Initialize swath names 
    call omi_set_fitting_parameters ( pge_error_status )
    if ( pge_error_status >= pge_errstat_error ) return

    !FIXME - remove when errstat more widely used
    errstat = 0
    ! Initialize ntimes, nxtrack, nfxtrack (UV-1)
    if (use_tio_in) then ! read from netCDF file
      ext=index(l1b_rad_filename,'.he4')
      if (ext > 0) then ! input is .he4
        l1_rad_filename_nc=l1b_rad_filename(1:ext-1)//'.nc'
      else
        ext=index(l1b_rad_filename,'.nc')
        if (ext > 0) then ! input is .nc
          l1_rad_filename_nc = l1b_rad_filename
          l1b_rad_filename = l1b_rad_filename(1:ext-1)//'.he4'
        else ! input has no extension
          l1_rad_filename_nc=trim(l1b_rad_filename)//'.nc'
          l1b_rad_filename = trim(l1b_rad_filename)//'.he4'
        endif
      endif
      ext=index(l1b_irrad_filename,'.he4')
      if (ext > 0) then ! extension is .he4
        l1_irrad_filename_nc=l1b_irrad_filename(1:ext-1)//'.nc'
      else
        ext=index(l1b_irrad_filename,'.nc')
        if (ext > 0) then ! extension is .nc
          l1_irrad_filename_nc=l1b_irrad_filename
          l1b_irrad_filename=l1b_irrad_filename(1:ext-1)//'.he4'
        else
          l1_irrad_filename_nc=trim(l1b_irrad_filename)//'.he4'
          l1b_irrad_filename=trim(l1b_irrad_filename)//'.he4'
        endif
      endif
      call read_l1_dims_tio (l1_rad_filename_nc, nc_rad_swathname(1), &
           ntimes, nfxtrack, nxtrack, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_error, &
             "omi_fitting_process: failed to read dimensions", errstat)
        return
      endif
    else
      ! read from he4 OMI file
      call omi_read_radiance_paras (pge_error_status )
      if ( pge_error_status >= pge_errstat_error ) return
    endif

    ! xtrack positions to be processed (start from the actual binned position)
    if (linenum_lim(2) >= ntimes)  linenum_lim(2) = ntimes
    if (zoom_mode) then
      pixnum_lim(1) = max(nint(1.0 * (zoom_p1 + 1) / ncoadd), pixnum_lim(1))
      pixnum_lim(2) = min(zoom_p2 / ncoadd, pixnum_lim(2))
      if ( mod(pixnum_lim(2) - pixnum_lim(1) + 1, nxbin) /= 0) then
        write(www_lun, '(A,2I4)') 'Incorrect across track binning option: ', pixnum_lim(1:2)
        pge_error_status = pge_errstat_error; return
      endif

      if (pixnum_lim(1) > pixnum_lim(2)) then
        write(www_lun, *) 'This is a zoom in mode orbit!!!'
        write(www_lun, '(A, I2, A, I2)') 'Invalid pixel selection, must be between ', &
             (zoom_p1 / ncoadd) * ncoadd + 1, ' and ', (zoom_p2 / ncoadd ) * ncoadd 
        pge_error_status = pge_errstat_error; return
      endif
    endif

    first_pix  = ceiling(1.0 * pixnum_lim(1) / nxbin)
    last_pix  = nint(1.0 * pixnum_lim(2) / nxbin )
    nxcoadd = nxbin * ncoadd

    ! Line number, starting from zero and keep track of offset
    offset_line = linenum_lim(1) - 1; first_line = 1
    last_line   = int ((linenum_lim(2) - linenum_lim(1) + 1.0) / nybin) 
    if (last_line > 100) then
      call tell_error(tell_usage_error, &
           "omi_fitting_process: number of binned lines defined in PCF exceeds maximum of 100", errstat)
      if (scnwrt) print *, 'number of binned lines (>100) = ', last_line
      pge_error_status = pge_errstat_error
      return
    endif


    ! load omi slit parameters  (Need to use it while getting coadded irradiance data)
    if (which_slit == 4) then
      call load_slitpars (pge_error_status)
      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A)') 'Finish loading OMI slit parameters !!!'
    endif

    !If using nc inputs only, need to read the obs date a this point,
    ! as it will be needed if the backup solar spectrum is used
    if (use_tio_in) then
      errstat = 0 ! FIXME - remove when libtell more widely used
      call read_date_tio (l1_rad_filename_nc, year, month, day, jday, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_error, &
             "omi_fitting_process: failed to read date", &
             errstat)
        return
      else
        GranuleYear = year
        GranuleMonth = month
        GranuleDay = day
      endif
    endif

    ! Read OMI irradiances
    reduce_resolution_save = reduce_resolution; reduce_resolution = .false.
    call omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, pge_error_status ) 
    if ( pge_error_status >= pge_errstat_error ) then
      use_backup = .true.
      call omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, pge_error_status )  
      if ( pge_error_status >= pge_errstat_error ) return   
    endif
    if (scnwrt) write(*, '(A)') 'Finish reading irradiances!!!'

    if (use_solcomp) then
      call replace_solar_irradiance(calunit, nxcoadd, first_pix, last_pix, pge_error_status)
      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A)') 'Finish replacing irradiances with solar composite!!!'
    endif

    ! Calibrate OMI irradiances for all cross track pixels and save results
    if (reduce_resolution_save) reduce_resolution = .true.
    call omi_irrad_cross_calibrate (first_pix, last_pix, pge_error_status)
    if ( pge_error_status >= pge_errstat_error ) return
    if (scnwrt) write(*, '(A)') 'Finish calibrating irradiances!!!'

    if (reduce_resolution_save) then
      call omi_read_irradiance_data (calunit, nxcoadd, first_pix, last_pix, pge_error_status ) 
      if ( pge_error_status >= pge_errstat_error ) return  
      if (scnwrt) write(*, '(A)') 'Finish reading irradiances (with reduced resolution)!!!'      
    endif


    !FIXME
    ! At present, compute_pixel_corners sets values in arrays whose size
    ! is determined by nlines_max = 100. This suggests that it expects
    ! to operate within one OMI data black of 100 lines. When operated
    ! in it's current position, outside the OMIBlock loop, if last_line >100
    ! the array indices will exceed nlines_max, causing memory to be 
    ! overwritten and causing the program to crash / give false results
    ! Consider moving compute_pixel_corners inside the OMIBlock loop.

    ! Compute spatial pixel corners and effective viewing geometry
    if (use_he5_in) then
      call compute_pixel_corners ( ntimes, nfxtrack, last_line, &
           pge_error_status)
      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A)') 'Finish computing pixel corners!!!'
    endif
    if (use_tio_in) then
      call read_geo_tio (l1_rad_filename_nc, nc_rad_swathname(1), ntimes, &
           nfxtrack, last_line, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "omi_fitting_process: failed to read geolocation data", errstat)
        pge_error_status = pge_errstat_error
        return
      else if (scnwrt) then
        print *, 'Finished reading geolocation data'
      endif
    endif
    ! For TEMPO we will only have one cloud product to read from
    !    ! Determine OMCLDRR or OMCLDO2
    !    i = INDEX(l2_cld_filename, '-o') - 22
    !    cldtype = l2_cld_filename(i : i + 6)
    !    IF (cldtype == 'OMCLDO2') THEN
    !      IF (which_cld == 0) which_cld = 1
    !      IF (which_cld == 3) which_cld = 4
    !    ENDIF
    !
    !    IF (which_cld == 4 .OR. which_cld == 1) THEN
    !      CALL read_omicldo2_clouds (ntimes, nxtrack, last_line, pge_error_status)
    !      IF ( pge_error_status >= pge_errstat_error ) RETURN
    !      IF (scnwrt) WRITE(*, '(A)') 'Finish reading omi L2 clouds!!!'
    !    ELSE IF (which_cld == 3 .OR. which_cld == 0) THEN

    if (use_he5_in) then
      call read_omicldrr_clouds (ntimes, nxtrack, last_line, &
           pge_error_status)
      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A)') 'Finish reading L2 clouds from he5!!!'
    endif
    if (use_tio_in) then
      !NetCDF cloud file
      ext=index(l2_cld_filename,'.he5')
      if (ext > 0) then ! extension is .he5
        l2_cld_filename_nc=l2_cld_filename(1:ext-1)//'.nc'
      else
        ext=index(l2_cld_filename,'.nc')
        if (ext > 0) then ! extension is .nc
          l2_cld_filename_nc=l2_cld_filename
          l2_cld_filename=l2_cld_filename(1:ext-1)//'.he5'
        else !no extension
          l2_cld_filename=trim(l2_cld_filename)//'.he5'
          l2_cld_filename_nc=trim(l2_cld_filename)//'.nc'
        endif
      endif
      call read_cloud_tio(l2_cld_filename_nc, ntimes, nxtrack, last_line, &
           errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_write_error, &
             "omi_fitting_process: failed to read cloud data", errstat)
      else
        OMIL2_clouds%cfr = L2_cloud%cfr
        OMIL2_clouds%ctp = L2_cloud%ctp
        OMIL2_clouds%qflags = L2_cloud%qflags
      endif
      if (errstat /= 0) then
        !FIXME - remove when error handling becomes unified
        pge_error_status = pge_errstat_error
        return
      endif
      if (scnwrt) write(*, '(A)') 'Finish reading L2 clouds fron netCDF!!!'
    endif

    !    ENDIF

    ! Initialize fitting statistics
    npix_fitting = 0       ! number of pixels (failure + success)
    npix_fitted  = 0       ! number of successfully fitted pixels   
    fitcol_avg   = 0.0
    rms_avg = 0.0
    dfitcol_avg = 0.0
    drel_fitcol_avg = 0.0

    ! Open output files
    if (lcurve_write) then
      open(UNIT=lcurve_unit, file=trim(adjustl(lcurve_fname)), status='unknown', IOSTAT=errstat)
      if ( errstat /= pge_errstat_ok ) then
        write(www_lun, *) modulename, ': Cannot open lcurve file!!!'
        pge_error_status = pge_errstat_error; return
      end if
    endif
    if (ozwrtint) then
      open(UNIT=ozwrtint_unit, file=trim(adjustl(ozwrtint_fname)), status='unknown', IOSTAT=errstat)
      if ( errstat /= pge_errstat_ok ) then
        write(www_lun, *) modulename, ': Cannot open intermediate output file!!!'
        pge_error_status = pge_errstat_error; return
      end if
    endif


    version = 1
    if (use_he5_in) then
      errstat = OMI_getCoreMetaData( L1B_UV_FILE_LUN, version, &
           L1BcoreMeta, year, month, day, jday )
      if( errstat /= OMI_S_SUCCESS ) then
        write(www_lun, *)  modulename, "OMI_getCoreMetaData failed "
        pge_error_status = pge_errstat_error; return
      end if
    endif
    !! read_date_tio moved up from here



    if (l2_hdf_flag == 0) then
      ! text output
      open (UNIT=l2funit, FILE=trim(adjustl(l2_filename)), STATUS='UNKNOWN', IOSTAT=errstat)
      if ( errstat /= pge_errstat_ok ) then
        write(www_lun, *) modulename, ': Cannot open output file!!!'
        pge_error_status = pge_errstat_error; return
      end if
      call timestamp(currtime)
      write(l2funit, '(3A,1x,A27,A10,I5,A10,I5)') trim(adjustl(algorithm_name)), ', ', &
           trim(adjustl(algorithm_version)), currtime, ' xbin = ', nxbin, ' ybin = ', nybin

    else 

      ! he5 output
      if (use_he5_out) then
        call He5_L2WrtInit (first_pix, last_pix, first_line, last_line, errstat )
        if ( errstat /= pge_errstat_ok ) then
          write(www_lun, *) modulename, ' : Cannot create HE5 output file!!!'
          pge_error_status = pge_errstat_error; return
        end if

        call He5_L2SetGeoFields (first_pix, last_pix, errstat )
        if ( errstat /= pge_errstat_ok ) then
          write(www_lun, *) modulename, ' : Cannot write geolocation fields!!!'
          pge_error_status = pge_errstat_error; return
        endif
      endif

      ! netCDF output
      if (use_tio_out) then
        ext = index (l2_filename, '.he5', back=.true.)
        if (ext /= 0) then
          nc_l2_filename = l2_filename(:ext)//'nc'//c_null_char
        else
          nc_l2_filename = l2_filename !assume filename is ok as-is.
        endif
        num_param = n_fitvar_rad - nfgas - &
             (ozfit_end_index - ozfit_start_index + 1)
        num_wav_max = maxval(omi_nwav_irrad(first_pix:last_pix)) - &
             numwin * 2 * radnhtrunc
        call l2_tio_create(nc_l2_filename, first_pix, last_pix, first_line, &
             last_line, offset_line, nybin, nfgas, nlay, n_fitvar_rad, &
             numwin, num_param, num_wav_max, errstat)
        if (errstat < 0) then
          call tell_error (tell_io_write_error, &
               "omi_fitting_process: L2 file creation failed", errstat)
        endif
        ! write geolocation data to netCDF
        call l2_tio_write_geo(first_pix, last_pix, first_line, last_line, errstat)
      endif

    endif

    ! loop through each OMI data block
    ! 1. read each block
    ! 2. perform calibration for each block (middle line)
    ! 3. perform retrievals from first_pix (all lines within the block) to last_pix

    ! these pixels with exitval >= 0 will be used by subsequent retrievals as initial values
    omi_exitval = -10 ! no retrievals yet 
    omi_fitvar  = 0.0

    OMIBlock: do iline = 0, last_line-1, nlines_max
      ntimes_loop = nlines_max

      if ( iline + ntimes_loop > last_line )   ntimes_loop = last_line - iline

      ! Actually lines in OMI data
      sline = offset_line + iline * nybin
      eline = sline + ntimes_loop * nybin - 1     

      ! Get NTIMES_LOOP radiance lines (with effective viewing geometry)
      call omi_read_radiance_lines (iline, ntimes_loop, sline, eline, first_pix, last_pix, nxcoadd, pge_error_status)
      if ( pge_error_status >= pge_errstat_error ) return
      if (scnwrt) write(*, '(A,I4,A,I4)') 'Finishing reading radiances for lines: ', sline + 1, ' - ', eline+1

      ! Perform calibration for the middle scan line and apply to the other scan lines
      if (wavcal) then
        call omi_rad_cross_calibrate (first_pix, last_pix, pge_error_status)
        if ( pge_error_status >= pge_errstat_error ) return
        if (scnwrt) write(*, '(A,I4,A,I4)') 'Finishing calibrating radiances for lines: ', sline+1, ' - ', eline+1
      endif

      ! loop through each xtrack position
      ! 1. Prepare databases, adjust radiances
      ! 2. Process all pixels at this poistion
      XtrackPix: do currpix = first_pix, last_pix
        !! Kai
        if( scnwrt ) write(*,'(A,I5,A,I3)') 'Doing Line=',iline,', iX=',currpix
        ! Need to convolve high-resolution ozone absorption cross section (for this position)
        ! Once the xsection is convolved, it will be set to false in ROUTINE getabs_crs
        ozabs_convl = .true.; so2crs_convl = .true.

        !IF (ALL(omi_radpix_errstat(currpix, 0:ntimes_loop-1) == pge_errstat_error) &
        !     .OR. omi_solpix_errstat(currpix) == pge_errstat_error) THEN
        !   omi_exitval(currpix, 0:ntimes_loop-1) = -10; CYCLE
        !ENDIF

        ! Load/adjust irradiances and slit calibration parameters
        call omi_adj_solar_data (pge_error_status)
        if ( pge_error_status >= pge_errstat_error ) cycle

        curr_fitted_line = 0
        YfitLine: do currloop = 0, ntimes_loop - 1

          omi_exitval(currpix, currloop) = -10
          currline = iline + currloop  

          if (omi_radpix_errstat(currpix, currloop) == pge_errstat_error .or. &
               omi_solpix_errstat(currpix) == pge_errstat_error ) &
               omi_exitval(currpix, currloop) = -9

          ! Load/adjust radiances/geolocations fields for a particular pixel
          ! Prepare databases for the first pixel (ifitline == 1)   
          if (omi_exitval(currpix, currloop) == -10) then
            call omi_adj_earthshine_data (curr_fitted_line, pge_error_status)
            if ( pge_error_status >= pge_errstat_error ) omi_exitval(currpix, currloop) = -9
          endif

          if (omi_exitval(currpix, currloop) == -10) then
            if (scnwrt) write(*, '(A,I5,A10,I5, A10, I5)')       'OMI Pixel: Line = ', &
                 currline * nybin + offset_line + 1, ' XPix = ', (currpix-1) * nxbin  + 1, &
                 ' Loop = ', currloop

            initval = omi_initval(currpix, currloop)
            call specfit_ozprof (initval, fitcol, dfitcol, rms, exval)
            omi_exitval(currpix, currloop) = exval  ! Store exit status for current pixel
            omi_fitvar(currpix, currloop, 1:n_fitvar_rad) &
                 = fitvar_rad_saved(mask_fitvar_rad(1:n_fitvar_rad))
          else
            exval = -9
          endif

          ! Write retrievals
          if (l2_hdf_flag == 0) then
            ! text output
            if (exval > -9) call omi_write_intermed (l2funit, fitcol, &
                 dfitcol,  rms,  exval)
          else 
            ! write he5 output
            if (use_he5_out) then
              call He5_L2SetDataFields (currpix, first_pix, last_pix, &
                   currloop, currline, ntimes_loop, exval, fitcol, dfitcol, &
                   pge_error_status )
              if ( pge_error_status >= pge_errstat_error ) return
            endif
            ! write netCDF output
            if (use_tio_out) then
              ix = currpix - first_pix ! start from zero for l2_tio_write_data
              if (exval >= 0) then  ! Retrieval finished.
                call l2_tio_write_data (ix, currloop, &
                     exval, fitcol, dfitcol, nfgas, nlay, n_fitvar_rad, &
                     numwin, num_param, num_wav_max, ozfit_start_index, &
                     ozfit_end_index, errstat)
                if (errstat < 0) then
                  call tell_error (tell_io_write_error, &
                       "omi_fitting_process: L2 write failed", errstat)
                endif
              else  ! Retrieval failed. Fill in as missing values.
                call l2_tio_fill_data (ix, currloop, &
                     exval, nfgas, nlay, n_fitvar_rad, numwin, num_param, &
                     num_wav_max, ozfit_start_index, ozfit_end_index, errstat)
                if (errstat < 0) then
                  call tell_error (tell_io_write_error, &
                       "omi_fitting_process: L2 fill failed", errstat)
                endif
              endif
            endif
          endif

          if ( exval >= 0 .and. fitcol(1) > 0.0 .and. dfitcol(1, 1) >= 0.0 ) then                     
            ! ----------------------------------------------------------------------
            ! Some general statistics on the average fitted column and uncertainty.
            ! Again, we make sure that only "good" fits are included in the average.
            ! ----------------------------------------------------------------------              
            fitcol_avg       = fitcol_avg + fitcol(1)
            rms_avg          = rms_avg + rms
            dfitcol_avg      = dfitcol_avg + dfitcol(1, 1)
            drel_fitcol_avg  = drel_fitcol_avg + dfitcol(1, 1) / fitcol(1)
            npix_fitted      = npix_fitted + 1
            curr_fitted_line = curr_fitted_line + 1
          endif

          npix_fitting        = npix_fitting + 1          
        enddo YfitLine
      enddo XtrackPix
    enddo OMIBlock

    ! ------------
    ! Final output
    ! ------------
    PROFOZ_metaItems(1:6) = (/ it_ParameterName,      &
         it_QAPercentMissingData,    &
         it_QAPercentOutofBoundsData,&
         it_QAPercentCloudCover,     &
         it_AutoQaFlagExpl,          &
         it_AutomaticQualityFlag /)

    ! he5 metadata
    if (use_he5_out) then
      if( use_backup ) then
        LUNinputPointer(1:2) = (/ L1B_UV_FILE_LUN, L2_CLD_FILE_LUN /)
        errstat  = OMI_setCoreArchMetaData( L2_OUT_LUN , L1BcoreMeta, &
             LUNinputPointer(1:2), MCF_FILE_LUN,  &
             PROFOZ_metaItems(1:6), ShortName ) 
      else
        LUNinputPointer(1:3) = (/ L1B_IRR_FILE_LUN, L1B_UV_FILE_LUN, &
             L2_CLD_FILE_LUN /)
        errstat  = OMI_setCoreArchMetaData( L2_OUT_LUN , L1BcoreMeta, &
             LUNinputPointer(1:3), MCF_FILE_LUN,  &
             PROFOZ_metaItems(1:6), ShortName ) 
      endif

      if( errstat /= OMI_S_SUCCESS ) then
        write(*,*) trim(modulename)//': Set OMI_setCoreArchMetaData failed.'
        pge_error_status = pge_errstat_error; return
      end if
    endif

    if ( npix_fitted == 0) npix_fitted = 1
    if (scnwrt) call write_final(fitcol_avg, rms_avg, dfitcol_avg, drel_fitcol_avg, npix_fitted)
    if (scnwrt) write(*, '(2(A,I5))') 'Number of pixels = ', &
         npix_fitting, '   Number of fitted pixels = ', npix_fitted

    ! -----------------------------------------
    ! Close L1 radiance file and L2 output file
    ! -----------------------------------------
    if (l2_hdf_flag == 0) then
      call timestamp(currtime)
      write(l2funit, '(A27)') currtime 
      close ( l2funit )
    endif

    if (lcurve_write) close (lcurve_unit)
    if (ozwrtint)     close (ozwrtint_unit)

    ! If using netCDF inputs, copy critical metadata and label file
    if (use_tio_in .AND. use_tio_out) then
      call copy_hdr_metadata (l1_rad_filename_nc, errstat)
      call label_output_file ("o3p", processing_version, errstat)
    endif


    ! Close L2 netCDF output file
    if (use_tio_out) then
      call l2_tio_close(errstat)
      if (errstat < 0) then
        call tell_error (tell_io_write_error, &
             "omi_fitting_process: Failed to close L2 file", errstat)
        stop 1
      endif
    endif

    return
  end subroutine omi_fitting_process

end module m_omi_fitting_process
