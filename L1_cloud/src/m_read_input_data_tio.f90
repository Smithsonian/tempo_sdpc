!read input data from L1B netCDF file
module m_read_input_data_tio

  private
  public read_input_data_tio

contains

  subroutine read_input_data_tio(l1bfile, errstat)
    use m_vars, only: iprt, input_data_path, iLine, wrt_solar, wmin, wmax, &
         wmin2, wmax2, set_wmin, set_wmax, wave_long, wave_short, nLines, &
         start_line, max_lines, nXtrack, nTimes, nWavel, nWavelCoef, &
         time, lat, lon, sza, sazimuth, sat_zen, vazimuth, terr_height, &
         geoflg, anomflg
    use m_strpos
    use m_instr_config
    use tio_module
    use tell_module
    use cld_names_module
    use netcdf, only : nf90_nowrite

    implicit none

    !input variables
    character (len=*), intent (in) :: l1bfile 

    !output variables
    integer (kind=4), intent (inout) :: errstat

    !local variables
    character (len=200) :: filenamepath, swathname
    integer (kind = 4) :: PGS_TD_TAItoUTC 
    logical :: uvswath

    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

 
    !-----------------------------------------------------------------
    !Note that unlike he5 version, this reads in all data in one pass,
    ! so no need for looping!
    !-----------------------------------------------------------------
    !First, perform some setup
    !set file path+name, swath name
    filenamepath=trim(input_data_path)//l1bfile
    uvswath = strpos (l1bfile, 'BRUG') > 0
    if (uvswath) then
      swathname = "band_540_740_nm"
    else
      call tell_error (tell_io_write_error, &
           "read_input_data_tio: input file is not OMI L1 UV swath", &
           errstat)
      return
    endif
      !set wavelength bounds
   !   if (wrt_solar) then
   !     wmin=355.0d0
   !     wmax=500.0d0
   !     wmin2=310.0d0
   !     wmax2=375.0d0
   !   endif
   !   wmin2 = 330.0d0
   !   wmax2 = 367.0d0
   !   if (.not. set_wmin) wmin = 345.5d0
   !   if (.not. set_wmax) wmax = 354.5d0
   !   wave_long=362.5d0 
   !   wave_short=345.4d0
      
    !open the file and get dimensions
    if (iprt > 0) print *,'read_input_data_tio: filename ', &
         trim(filenamepath), '   ', trim(swathname)
    call tiof_open (l1bfile, tio_l1obj, nf90_nowrite, errstat)
    call tiof_inq_group (tio_l1obj, swathname, errstat)
    call tiof_inq_dimlen (tio_l1obj, "xtrack", nXtrack, errstat)
    call tiof_inq_dimlen (tio_l1obj, "mirror_step", nTimes, errstat)
    call tiof_inq_dimlen (tio_l1obj, "wavelengths", nWavel, errstat)
    if(iprt > 0) then
      print *,'read_input_data: nTimes, nXtrack, nWavel, nWavelCoef '
      print *, nTimes,nXtrack,nWavel,nWavelCoef
    endif
    iLine=start_line
    if (max_lines > 0 .and. iprt > 0) then
      print *,'read_input_data_tio: changing nTimes to ',max_lines
      nTimes=max_lines+start_line
    endif
    
    if (errstat > 0) then
      call tell_error (tell_io_write_error, &
           "read_input_data_tio: failed to open L1B file", &
           errstat)
      return
    endif
    
    !------------------------------------------------------------------
    nLines=nTimes

    !Allocate arrays for variables to be read in
    call alloc_scan()

    !------------------------------------------------------------------
    ! Read in the Geolocation data
    call tiof_get1d_r8 (tio_l1obj, cld_var_time, [iLine], [nTimes], &
         time, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_latitude, [iLine,0], [nTimes,-1], &
         lat, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_longitude, [iLine,0], [nTimes,-1], &
         lon, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_sz_angle, [iLine,0], [nTimes,-1], &
         sza, errstat)
    call tiof_get2d_r4 (tio_l1obj, "solar_azimuth_angle", [iLine,0], &
         [nTimes,-1], sazimuth, errstat)
    call tiof_get2d_r4 (tio_l1obj, cld_var_vz_angle, [iLine,0], [nTimes,-1], &
         sat_zen, errstat)
    call tiof_get2d_r4 (tio_l1obj, "viewing_azimuth_angle", [iLine,0], &
         [nTimes,-1], vazimuth, errstat)
    call tiof_get2d_i2 (tio_l1obj, "ellipsoid_altitude", [iLine,0], &
         [nTimes,-1], terr_height, errstat)
    call tiof_get2d_i2 (tio_l1obj, "GroundPixelQualityFlags", [iLine,0], &
         [nTimes,-1], geoflg, errstat)
    call tiof_get2d_i1 (tio_l1obj, "XTrackQualityFlags", [iLine,0], &
         [nTimes,-1], anomflg, errstat)


  end subroutine read_input_data_tio





  !allocate memory for variables !!! NEEDS UPDATING!!!
  subroutine alloc_scan()

    use m_vars
    use m_pgs_include

    implicit none

    ! allocate memory for arrays
    if (allocated(lat)) deallocate (lat)   
    ALLOCATE( lat(nXtrack,nlines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "latitude allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    if (allocated(lon)) deallocate (lon)   
    ALLOCATE( lon(nXtrack,nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "longitude allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    if (allocated(sza)) deallocate (sza)
    ALLOCATE( sza(0:nXtrack-1,nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN 
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, & 
!           "szenith allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    if (allocated(sat_zen)) deallocate (sat_zen)
    ALLOCATE( sat_zen(0:nXtrack-1,nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN 
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, & 
!           "vzenith allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( sazimuth(nXtrack,nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN 
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "sazimuth allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( vazimuth(nXtrack,nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN 
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "vazimuth allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( terr_height(nXtrack,nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN 
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "terrain height allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( geoflg(nXtrack,nLines), STAT=ierr )
    geoflg=0 
!    IF( ierr .NE. zero ) THEN 
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "geoflg allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( anomflg(nXtrack,nLines), STAT=ierr )
    anomflg=0
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "anomflg allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( mflg(nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN 
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "measflg allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( quality_flagL(nWavel,nXtrack), STAT=ierr )
    quality_flagL=0
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "quality_flagL allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    if (allocated(w12d)) deallocate (w12d)
    ALLOCATE( w12d(0:nWavel-1,0:nXtrack-1), STAT=ierr )   ; w12d = 0.0
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "wavelengthL allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    if (allocated(f12d)) deallocate (f12d)
    ALLOCATE( f12d(0:nWavel-1,0:nXtrack-1), STAT=ierr ) ; f12d=0.
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "radianceL allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( time(nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "Time allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    ALLOCATE( meas_class(nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "MeasurementClass allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF
    ALLOCATE( config_rad(nLines), STAT=ierr )
!    IF( ierr .NE. zero ) THEN
!      ierr = OMI_SMF_setmsg( OMCLDRR_F_MEM_ALLOC, &
!           "InstrumentConfigurationID allocation failure, PGE aborting, exit code = 1", &
!           "read_input_data", 1 )
!      call exit(1)
!    END IF

    if (allocated(meas_qual_flg)) deallocate (meas_qual_flg)
    allocate  (meas_qual_flg(nLines)) ; meas_qual_flg = 0
    if (allocated(cloud_pres)) deallocate (cloud_pres)
    allocate (cloud_pres (0:nXtrack - 1,nLines)) ; cloud_pres=fill_value
    if (allocated(azimuth)) deallocate (azimuth)
    allocate (azimuth (0:nXtrack-1,nLines))      ; azimuth=fill_value
    if (allocated(refl)) deallocate (refl)  
    allocate (refl    (0:nXtrack-1,nLines))      ; refl=fill_value
    if (allocated(dIdR)) deallocate (dIdR)  
    allocate (dIdR    (0:nXtrack-1,nLines))      ; dIdR=fill_value
    if (allocated(ps)) deallocate (ps)     
    allocate (ps      (0:nXtrack-1,nLines))      ; ps=fill_value
    if (allocated(ref_clr)) deallocate (ref_clr)     
    allocate (ref_clr      (0:nXtrack-1,nLines))      ; ref_clr=fill_value
    if (allocated(reflect_cld)) deallocate (reflect_cld)      
    allocate (reflect_cld      (0:nXtrack-1,nLines)) ; reflect_cld=fill_value
    if (allocated(rad_cld_frac)) deallocate (rad_cld_frac)
    allocate (rad_cld_frac(0:nXtrack-1,nLines))  ; rad_cld_frac=fill_value
    if (allocated(eff_cld_frac)) deallocate (eff_cld_frac)
    allocate (eff_cld_frac(0:nXtrack-1,nLines))  ; eff_cld_frac=fill_value
    if (allocated(eff_cld_frac2)) deallocate (eff_cld_frac2)
    allocate (eff_cld_frac2(0:nXtrack-1,nLines))  ; eff_cld_frac2=fill_value
    if (allocated(cld_pres2)) deallocate (cld_pres2)
    allocate (cld_pres2(0:nXtrack-1,nLines))  ; cld_pres2=fill_value
    if (allocated(chlorophyll)) deallocate (chlorophyll)
    allocate (chlorophyll(0:nXtrack - 1,nLines)) ; chlorophyll=fill_value
    if (allocated(biases)) deallocate (biases)  
    allocate (biases  (0:nXtrack-1,nLines))      ; biases=fill_value
    if (allocated(biases2)) deallocate (biases2)  
    allocate (biases2  (0:nXtrack-1,nLines))      ; biases2=fill_value
    if (allocated(stds)) deallocate (stds) 
    allocate (stds    (0:nXtrack-1,nLines))      ; stds=fill_value
    if (allocated(stds2)) deallocate (stds2) 
    allocate (stds2    (0:nXtrack-1,nLines))      ; stds2=fill_value
    if (allocated(chi_sqr)) deallocate (chi_sqr)
    allocate (chi_sqr (0:nXtrack-1,nLines))      ; chi_sqr=fill_value
    if (allocated(chi_sqr2)) deallocate (chi_sqr2)
    allocate (chi_sqr2 (0:nXtrack-1,nLines))      ; chi_sqr2=fill_value
    if (allocated(land_flg)) deallocate (land_flg)
    allocate (land_flg(0:nXtrack-1)) ; land_flg=.FALSE.
    if (allocated(chlcl)) deallocate (chlcl)
    allocate (chlcl   (0:nXtrack-1)) ; chlcl=fill_value
    if (allocated(qc)) deallocate (qc)      
    allocate (qc      (0:nXtrack-1,nLines)) ; qc=0
    if (allocated(qc2)) deallocate (qc2)      
    allocate (qc2      (0:nXtrack-1,nLines)) ; qc2=0
    if (allocated(fill)) deallocate (fill)      
    allocate (fill    (0:nXtrack-1,nLines)) ; fill=fill_value
    if (allocated(shifts)) deallocate (shifts)      
    allocate (shifts  (0:nXtrack-1,nLines)) ; shifts=fill_value
    if (allocated(shifts2)) deallocate (shifts2)      
    allocate (shifts2  (0:nXtrack-1,nLines)) ; shifts2=fill_value
    if (allocated(squeezes)) deallocate (squeezes)      
    allocate (squeezes(0:nXtrack-1,nLines)) ; squeezes=1

  end subroutine alloc_scan





end module m_read_input_data_tio
