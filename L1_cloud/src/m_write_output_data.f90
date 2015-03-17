module m_write_output_data

contains

  subroutine write_output_data(outfile, swathname)
    use m_write_swath_field
    use m_vars
    USE ISO_C_BINDING, ONLY: C_LONG
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  write_output_data
    ! 
    ! !DESCRIPTION: write_output_data writes OMI cloud level 2 data to HDF-EOS
    !		file
    !
    ! !CALLING SEQUENCE: 
    !
    !        call write_output_data
    !     
    ! !INPUT PARAMETERS:   
    !
    ! !OUTPUT PARAMETERS:  
    !
    ! !SEE ALSO:  
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner     original fortran 90
    !  15Jan03   Vasilkov   write GroundpixelQualityFlags
    !  30Aug03   Vasilkov   write MeasurementQualityFlags
    !  30Aug07   Joiner     turn off CloudPressure/Fraction
    !
    !EOP
    !-------------------------------------------------------------------------
    !
    !inputs
    !------
    character(len=*), intent(in) :: outfile, swathname

    !ouputs
    !------
    !       integer, intent(out) :: ierr

    ! Declare the HDF-EOS file and swath identification numbers, and
    ! the status of the HDF-EOS functions calls.
    !-----------------------------------------------------------------
    integer (kind = 4) swfid, swid

    ! Declare the HDF-EOS functions.
    !-------------------------------
    integer (kind = 4) :: he5_swcreate, he5_swdefdim
    integer, parameter :: HE5_ACC_RDWR=100
    integer, parameter :: HE5_ACC_TRUNC=102
    integer (kind = 4) :: he5_swopen, he5_swattach, he5_swdetach, he5_swclose

    !Local variables
    !---------------
    character(len=255) :: nTimesstr="nTimes"
    character(len=255) :: nXtrackstr="nXtrack"
    character(len=255) :: nWavestr="nWavel"
    character(len=255) :: dims3
    character(len=255) :: dims2
    character(len=255) :: dims1
    character(len=255) :: dimsw
    integer            :: nTime, iLine1, nWaveRes
    real (kind = 4) :: misval_r4 = fill_value
    real (kind = 8) :: misval_r8 = fill_value
    integer (kind = 2) :: misval_i2 = fill_value_int 
    integer (kind = 1) :: misval_i1 = fill_value_int1
    !       integer (kind = 4) :: misval_i2 = fill_value_int 
    !       integer (kind = 3) :: misval_i1 = fill_value_int1

    iLine1=iLine-start_line
    nTime=1


    !Create the file 
    !===============

    ! Open the OMI Level 2 HDF-EOS cloud output file.
    !-------------------------------------------------
    if (iprt >= 1) &
         write(6,*) 'write_output_file: creating ',trim(outfile)
    swfid = he5_swopen (outfile, HE5_ACC_TRUNC)
    if (iprt >= 1) &
         write (6, *) 'write_output_file: swfid ', swfid

    ! Create the swath.
    !-----------------
    swid = he5_swcreate (swfid, swathname)
    if (iprt >= 1) &
         write (6, *) 'write_output_file: swid  ', swid

    ! Write the dimensions
    !---------------------
    if (iprt >= 1) &
         write (6,*) 'write_output_file: dimensions, nTimes ',nTimes, &
         ' nXtrack ',nXtrack
    status = he5_swdefdim(swid,nTimesstr,INT(nTimes, KIND=C_LONG))
    if (iprt >= 1) &
         write(6, *) status, 'after define nTimesstr'
    status = he5_swdefdim(swid,nXtrackstr,INT(nXtrack, KIND=C_LONG))
    if (iprt >= 1) &
         write(6, *) status, 'after define nXtrackstr'
    if (write_resid) then
      nWaveRes=size(wave_resid)
      status = he5_swdefdim(swid,nWavestr,INT(nWaveRes, KIND=C_LONG))
      if (iprt >= 1) &
           write(6, *) status, 'after define nWavestr'
    endif
    dimsw = trim(nWavestr)
    dims1 = trim(nTimesstr)
    dims2 = trim(nXtrackstr)//','//trim(nTimesstr)
    dims3 = trim(nWavestr)//','//trim(nXtrackstr)//','//trim(nTimesstr)

    ! Detach from the swath
    !----------------------
    status = he5_swdetach(swid)

    ! Close the swath file
    !----------------------
    status = he5_swclose(swfid)

    !Now write to the file all at once!
    !==========================================

    ! Reopen the file
    !----------------------
    swfid = he5_swopen (outfile, HE5_ACC_RDWR)

    ! Attach to the swath
    !----------------------
    swid = he5_swattach(swfid, swathname)

    ! Write the fields
    !----------------------
    status = put_data (swid, "Latitude", &
         dims2, lat, misval_r4, "Latitude", "deg", &
         geo=.true.,offset=(/0,0/),iprt=iprt)
    where(lon(:,:) > 180.) lon(:,:) = lon(:,:) -360.
    status = put_data (swid, "Longitude", &
         dims2, lon, misval_r4, "Longitude", "deg", &
         geo=.true.,offset=(/0,0/),iprt=iprt)
    status = put_data (swid, "Time", &
         dims1, time, misval_r8, "Time at Start of Scan (s, TAI93)", "s", geo=.true.,iprt=iprt)
    if (write_geom) then
      status = put_data (swid, "SolarZenithAngle", &
           dims2, sza, misval_r4, "Solar Zenith Angle", "deg", &
           geo=.true.,offset=(/0,0/),iprt=iprt)
      status = put_data (swid, "ViewingZenithAngle", &
           dims2, sat_zen, misval_r4, "Viewing Zenith Angle","deg", &
           geo=.true.,offset=(/0,0/),iprt=iprt)
      status = put_data (swid, "RelativeAzimuthAngle", &
           dims2, azimuth, misval_r4, "Relative Azimuth Angle", "deg", &
           geo=.true.,offset=(/0,0/),iprt=iprt)
    endif
    status = put_data (swid, "TerrainHeight", &
         dims2, terr_height, misval_i2, "Terrain Height", "m", geo=.true.,offset=(/0,0/),iprt=iprt)
    !      if (.not. cloud_clear) cloud_pres=cloud_pres*1013.25
    !      status = put_data (swid, "CloudPressure", &
    !        dims2, cloud_pres, misval_r4, "Cloud Pressure", "hPa", &
    !        offset=(/0,0/),iprt=iprt)
    !      status = put_data (swid, "CloudFraction", &
    !        dims2, eff_cld_frac, misval_r4, "Cloud Fraction", "NoUnits", &
    !        offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "RadiativeCloudFraction", &
         dims2, rad_cld_frac, misval_r4, "Radiative Cloud Fraction", &
         "NoUnits", &
         offset=(/0,0/), iprt=iprt)
    if ((do_cloud_mask) .and. (allocated(cloud_mask))) then 
      status = put_data (swid, "CloudMask", &
           dims2, cloud_mask, misval_i2, "Cloud Mask", "NoUnits", &
           offset=(/0,0/), iprt=iprt)
      status = put_data (swid, "SmallPixelStddev", &
           dims2, smpx_stddev, misval_r4, "Small Pixel Standard Deviation", &
           "NoUnits", offset=(/0,0/), iprt=iprt)
      status = put_data (swid, "SmallPixelMean", &
           dims2, smpx_mean, misval_r4, "Small Pixel Mean", &
           "NoUnits", offset=(/0,0/), iprt=iprt)
      status = put_data (swid, "SmallPixelnPix", &
           dims2, smpx_nPix, misval_i2, "Small Pixel Number of Pixels", &
           "NoUnits", offset=(/0,0/), iprt=iprt)
      !status = put_data (swid, "SmallPixelWave", &
      ! dims2, smpx_wavel, misval_r4, "Small Pixel Wavelength", &
      ! "nm", offset=(/0,0/), iprt=iprt)
    endif
    if (write_ps) then
      ps=ps*1013.25
      status = put_data (swid, "TerrainPressure", &
           dims2, ps, misval_r4, "Terrain Pressure", "hPa", &
           offset=(/0,0/),iprt=iprt)
    endif
    if (wr_shift) then
      status = put_data (swid, "WavelengthShift", &
           dims2, shifts2, misval_r4, "Wavelength Shift", "nm", &
           offset=(/0,0/),iprt=iprt)
      if (squeeze) status = put_data (swid, "WavelengthSqueeze", &
           dims2, squeezes, misval_r4, "Wavelength Squeeze", "NoUnits", &
           offset=(/0,0/),iprt=iprt)
    endif
    status = put_data (swid, "SurfaceReflectivity", &
         dims2, ref_clr, misval_r4, &
         "Surface Reflectivity Climatology TOMS", "NoUnits", &
         offset=(/0,0/),iprt=iprt)
    if (write_fill) then
      status = put_data (swid, "Filling-In", &
           dims2, fill, misval_r4, "Effective filling-in", "NoUnits", &
           offset=(/0,0/),iprt=iprt)
    endif
    if (write_resid) then
      status = put_data (swid, "WavelengthResiduals", &
           dimsw, wave_resid, misval_r4, "Wavelengths", "nm", &
           offset=0,iprt=iprt)
      status = put_data (swid, "Residuals", &
           dims3, resid, misval_r4, "Observed minus retrieved", "percent", &
           offset=(/0,0,0/),iprt=iprt)
    endif
    if (cal_reflec) then
      status = put_data (swid, "dIdR", &
           dims2, dIdR, misval_r4, "Radiance (fractional) refl. sens.", "NoUnits", &
           offset=(/0,0/),iprt=iprt)
    endif
    if (.not. do_mler) then
      status = put_data (swid, "Cloud_reflectivity", &
           dims2, reflect_cld, misval_r4, "Cloud reflectivity", "NoUnits", &
           offset=(/0,0/), iprt=iprt)
    endif
    status = put_data (swid, "Reflectivity", &
         dims2, refl, misval_r4, "Reflectivity", "NoUnits", &
         offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "MeasurementQualityFlags", &
         dims1, meas_qual_flg, &
         misval_i2, "Measurement Quality Flags", "NoUnits")
    status = put_data (swid, "GroundPixelQualityFlags", &
         dims2, geoflg, misval_i2, "Ground Pixel Quality Flags", "NoUnits", &
         geo=.true., offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "XTrackQualityFlags", &
         dims2, anomflg, misval_i1, "Cross Track Quality Flags", "NoUnits", &
         geo=.true., offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "Residual_bias", &
         dims2, biases2, misval_r4, "Residual Bias", "NoUnits", &
         offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "Residual_stddev", &
         dims2, stds2, misval_r4, "Residual Standard Deviation", "NoUnits", &
         offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "Convergence_factor", &
         dims2, chi_sqr2, misval_r4, "Convergence Factor", "NoUnits", &
         offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "Chlorophyll", &
         dims2, chlorophyll, misval_r4, "Chlorophyll Concentration", "mg/m3", &
         offset=(/0,0/), iprt=iprt)
    !       status = put_data (swid, "ProcessingQualityFlags", &
    !         dims2, qc, misval_i2, "Processing Quality Flags", "NoUnits", &
    !         offset=(/0,0/), iprt=iprt)

    ! Detach from the swath interface.
    !-------------------------------------------------
    status = he5_swdetach (swid)
    if (iprt >= 2) &
         write (6, *) 'write_output_data: detached swath ',status

    ! Close the OMI Level2  HDF-EOS output file.
    !----------------------------------------------
    status = he5_swclose (swfid)
    if (iprt >= 2) &
         write (6, *) 'write_output_data: closed file ', status

  end subroutine write_output_data

end module m_write_output_data
