!***************
program OMCDO2N
!***************

  use iso_fortran_env, only: output_unit
  use hdfeos4_parameters
  use he5_swreader
  use m_vars
  use m_read_input_kleipool
  use m_read_GMI, only: read_GMI_TMP
  use m_read_lut
  use m_read_hdf5
  use tell_module
  use tio_module
  use m_read_input_tio, only: read_rad_tio, read_irr_tio, read_cldo4_tio, &
               get_tio_global_attr, get_tio_l1rad_glbattr
  use m_write_output_tio, only: create_output_file_tio, &
               update_output_file_tio, write_debug_processing_flags
  use m_read_input_clim, only: read_geoscf
  use m_read_input_gler, only: read_gler
  use m_cal_ecf, only: cal_ecf, allocate_ecf_arrays
  use m_cal_ocp, only: cal_ocp, allocate_ocp_arrays
  use m_cal_pscene, only : cal_pscene
  use HDF5

  implicit none

  !============================
  ! read control.txt 
  !============================
  include 'GetConfig.inc'

  character(len=CFG_VAL_LEN)::buf
  integer::gmonth
  integer::status

  character(len=255)::filename, l1radfnm
  character(len=255)::name_gmi_psfc
  character(len=255)::name_gmi_tmp
  character(len=255)::name_kleipool_rsfc

  integer(kind=4) :: errstat
  integer(kind=4) :: nx, nt, iter_ecfocp
  real :: fspecial

  character(len=80) :: logmsg
  character(len=15), parameter :: swathname = "band_290_490_nm"

  call tell_open ("L1_cloud_o4", 0)
  errstat = 0
  nx = 0
  nt = 0
  fspecial = fFillValue ! large negative value in m_vars

  !call unbufferstdout()

  call tell_log(0, 'Starting L1_cloud_o4')
  ! ----------------------------------------
  ! READ IN INPUT FILENAMES from control.txt
  ! ----------------------------------------

  buf(:) = ' '
  status=GetConfigString("E","Input Files TEMPOL1RAD_fnm",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading TEMPOL1RAD_fnm from control file", errstat)
    call exit(-1)
  endif
  name_rad_file=trim(buf)

  status=GetConfigString("E","Input Files TEMPOL1RAD_dir",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading TEMPOL1RAD_dir from control file", errstat)
    call exit(-1)
  endif
  name_rad_dir=trim(buf)

  status=GetConfigString("E","Input Files TEMPOO4SCD_fnm",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading TEMPOO4SCD_fnm from control file", errstat)
    call exit(-1)
  endif
  name_nasa_file=trim(buf)

  status=GetConfigString("E","Input Files TEMPOO4SCD_dir",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading TEMPO4SCD_dir from control file", errstat)
    call exit(-1)
  endif
  name_nasa_dir=trim(buf)

! TEMPO GLER is now handled by m_read_input_gler
! explicit GLER fnm is no longer needed for TEMPO, but is needed for OMI
  status=GetConfigString("E","Input Files OMGLER",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading OMGLER from control file", errstat)
    call exit(-1)
  endif
  name_brdf_file=trim(buf)

  status=GetConfigString("E","Input Files TEMPOL1IRR_fnm",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading TEMPOL1IRR fnm from control file", errstat)
    call exit(-1)
  endif
  name_irr_file=trim(buf)

  status=GetConfigString("E","Input Files TEMPOL1IRR_dir",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading TEMPOL1IRR dir from control file", errstat)
    call exit(-1)
  endif
  name_irr_dir=trim(buf)

  ! -----------------------------------------
  ! READ IN OUTPUT FILENAMES from control.txt
  ! -----------------------------------------
  status=GetConfigString("E","Output Files TEMPOCLDO4",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading Output Filename from control file", errstat)
    call exit(-1)
  endif
  name_out_ncdf=trim(buf)

  gmetadata%localgranid=name_out_ncdf

  ! -----------------------------------------
  ! READ NAME OPTIONS from control.txt
  ! -----------------------------------------
  status=GetConfigString("E","Name Options T-P profile", buf)
  if (status < 0) then
     call tell_error(tell_io_read_error,"Problem reading T-P profile option", errstat)
     call exit(-1)
  endif
  name_option_TemperaturePressure=trim(buf)
  write(*,*) 'name_option_TemperaturePressure=',trim(name_option_TemperaturePressure)

  status=GetConfigString("E","Name Options Surface Reflectivity",buf)
  if (status < 0) then
     call tell_error(tell_io_read_error,"Problem reading Surface Reflectivity option", errstat)
     call exit(-1)
  endif
  name_option_SurfaceReflectivity=trim(buf)
  write(*,*) 'name_option_SurfaceReflectivity=',trim(name_option_SurfaceReflectivity)

  ! -------------------------------------------
  ! READ RUNTIME PARAMETERS from control.txt
  ! -------------------------------------------
  status=GetConfigString("E","Runtime Parameters RunMode",buf)
  if(status <0) then
    call tell_error(tell_io_read_error,"Problem reading RunMode from control file", errstat)
    call exit(-1)
  endif
  run_mode = trim(buf)
  write(*,*)'run_mode=',trim(run_mode)

  status=GetConfigString("W","Runtime Parameters APPShortName",buf)
  if(status .NE. 0) then
    write(*,*)"no APPShortName from control file, use default"
  ! will use m_vars default
  else
     gmetadata%appshortname=trim(buf)
  endif

  status=GetConfigString("W","Runtime Parameters APPVersion",buf)
  if(status .NE. 0) then
    write(*,*) "no APPVersion from control file, use default"
  ! will use m_vars default
  else
     gmetadata%appversion=trim(buf)
  endif

  status=GetConfigString("W","Runtime Parameters AuthorName",buf)
  if(status .NE. 0) then
    write(*,*)"no AuthorName from control file, use default"
  !  will use default in m_vars
  else
     gmetadata%author_name=trim(buf)
  endif

   status=GetConfigString("W","Runtime Parameters AuthorAffiliation",buf)
   if (status .NE. 0) then
      write(*,*)"no AuthorAffiliation from control file, use default"
   !  will use default in m_vars
   else
      gmetadata%author_affiliation=trim(buf)
   endif

  status=GetConfigString("W","Runtime Parameters ProcessingCenter",buf)
  if(status .NE. 0) then
    write(*,*)"no ProcessingCenter from control file, use default"
  ! will use default in m_vars 
  else 
     gmetadata%processingcenter=trim(buf)
  endif

  status=GetConfigString("W","Runtime Parameters TEMPO Footprint",buf)
  if(status .NE. 0) then
    write(*,*) "no Footprint channel from control file, use default"
  ! will use default in m_vars 
  else
     gmetadata%omiwindow=trim(buf)
  endif

  status=GetConfigString("W","Runtime Parameters Collection",buf)
  if(status .NE. 0) then
    write(*,*)"no Collection Number from control file, use default"
  ! will use default in m_vars
  else
     gmetadata%omi_collection=trim(buf)
  endif

  status=GetConfigString("W","Runtime Parameters ixdebug",buf)
  if (status .NE. 0) then
     write(*,*) "No ixdebug, use default instead"
  else
     read(buf,*,iostat=status) ixdebug
  endif
  write(*,*) 'ixdebug=',ixdebug

  status=GetConfigString("W","Runtime Parameters itdebug",buf)
  if (status .NE. 0) then
     write(*,*) "No itdebug, use default instead"
  else 
     read(buf,*,iostat=status) itdebug
  endif
  write(*,*) 'itdebug=',itdebug
  
  status=GetConfigString("W","Runtime Parameters ecfocp_maxiter",buf)
  if (status .NE. 0) then 
     write(*,*) "No ecfocp_maxiter in control file, use default instead"
  else
     read(buf,*,iostat=status) ecfocp_maxiter
  endif
  write(*,*) 'ecfocp_maxiter=',ecfocp_maxiter

  status=GetConfigString("W","Runtime Parameters option_destripe_scd",buf)
  if (status .NE. 0) then
     write(*,*) "use default option_destripe_scd"
  else 
     read(buf,*,iostat=status) option_destripe_scd
  endif
  write(*,*) 'option_destripe_scd=',option_destripe_scd

  if (option_destripe_scd .eq. 1) then
    status=GetConfigString("W","Input Files desfac_dir",buf)
    if(status .ne. 0) then
      write(*,*) 'No desfac_dir in control file, use default'
    !call tell_error(tell_io_read_error,&
    !   "Problem reading desfac_dir from control file", errstat)
    !call exit(-1)
    else
       name_desfac_dir=trim(buf)
    endif

    status=GetConfigString("W","Input Files desfac_fnm",buf)
    if(status .ne. 0) then
       write(*,*) 'No desfac_fnm in control file, usedefault'
    !call tell_error(tell_io_read_error,&
    !   "Problem reading desfac_fnm from control file", errstat)
    !call exit(-1)
    else
       name_desfac_fnm=trim(buf)
    endif
    write(*,*) 'desfac_filename=',trim(name_desfac_dir),trim(name_desfac_fnm)
  endif
  
  flush (output_unit)
  call tell_log(0,'Read control file')

  !===========================================
  ! 0. get global attributes from TEMPO L1 RAD to gmetadata
  !===========================================
  ! gmetadata is needed for climatology
  l1radfnm = trim(adjustl(name_rad_dir))//trim(adjustl(name_rad_file))
  call get_tio_l1rad_glbattr(l1radfnm,errstat)
  if (errstat < 0) then
     call tell_error(tell_io_read_error,"Error getting time_coverage_start from TEMPO L1 RAD" ,errstat)
     call exit(-1)
  endif

  !==============================
  ! 1.1 read inputs from radiance file 
  !==============================
  ! read rad for 440nm and 466nm
  ! also read lat/lon needed for GEOS-CF T-P
  ! allocate and initialize processing quality flags to zero
  call read_rad_tio (l1radfnm, swathname, errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_rad_tio failed', errstat)
    call exit(-1);
  endif
  flush (output_unit)
  call tell_log(0,'Read RAD '//l1radfnm)

  !============================================
  ! 1.2 read inputs from irradiance 
  !============================================
  ! Ewan tested read_irr_tio, appears to work 
  filename = trim(adjustl(name_irr_dir))//trim(adjustl(name_irr_file))
  call read_irr_tio(filename,swathname,errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_irr_tio failed', errstat)
    call exit(-1);
  endif
  flush (output_unit)
  call tell_log(0,'Read IRR '//filename)

  ! --------------------
  ! 1.3 read TEMPO O4 SCD from intermediate L2 file
  ! --------------------
  filename = trim(adjustl(name_nasa_dir))//trim(adjustl(name_nasa_file))
  call read_cldo4_tio(filename,errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_cldo4_tio failed', errstat)
    call exit(-1);
  endif
  flush (output_unit)
  call tell_log(0, 'Read O4 SCD '//filename)

  ! ---------------------
  ! 2. read T/P/Psfc 
  ! ---------------------
  ! 2.1. GMI
  !gmonth is used to decide filename for GMI
  gmonth = gmetadata%granule_month
  flush (output_unit)

  if(name_option_TemperaturePressure.eq.'GMI') then
    status=GetConfigString("E","Input Files GMI_dir",name_gmi_dir)
    status=GetConfigString("E","Input Files "//lun_gmi_psfc(gmonth),buf)
    name_gmi_psfc=trim(name_gmi_dir)//trim(buf)

    status=GetConfigString("E","Input Files "//lun_gmi_tmp(gmonth),buf)
    name_gmi_tmp=trim(name_gmi_dir)//trim(buf)

    write(*,*)'   name_gmi_psfc=',trim(name_gmi_psfc)
    write(*,*)'   name_gmi_tmp=',trim(name_gmi_tmp)
    call read_GMI_TMP(name_gmi_psfc,name_gmi_tmp, errstat)
    call tell_log(0, "Read GMI T/P & GMI Ps for month ")
  endif

  flush (output_unit)

  ! 2.2. GEOS-CF
  ! TEMPO operational uses GEOS5 as the option for GEOS-CF
  if(name_option_TemperaturePressure.eq.'GEOS5') then
  ! T-P is not needed in cal_ecf, but surface pressure is
     call read_geoscf (errstat)
     if (errstat /= 0) then
       call tell_error (tell_runtime_error, 'read_geoscf failed', errstat)
       call exit(-1)
     endif
     flush (output_unit)
     call tell_log(0,'Read GEOS-CF profiles')
  endif

  ! ------------------------------
  ! 3. read Rsfc from Kleipool or BRDF
  ! ------------------------------
  ! 3.1. Kleipool
  if(name_option_SurfaceReflectivity.eq.'Kleipool') then
    status=GetConfigString("E","Input Files Kleipool_dir",buf)
    name_kleipool_dir = trim(adjustl(buf))
    status=GetConfigString("E","Input Files Kleipool_fnm",buf)
    name_kleipool_rsfc=trim(name_kleipool_dir)//trim(adjustl(buf))
    call read_Kleipool_Rsfc(name_kleipool_rsfc,gmonth,errstat)
    flush (output_unit)
    call tell_log(0,'Read Kleipool Rsfc climatology')
  endif

  ! 3.2. GLER
  ! use BRDF option for TEMPO GLER
  if(name_option_SurfaceReflectivity.eq.'BRDF') then
    call read_gler (errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, 'read_gler failed', errstat)
      call exit(-1)
    endif
    if (name_option_TemperaturePressure .eq. 'GMI') then
        write(*,*)'   WARNING: wind speed = 0 for GMI+BRDF combo.'
    endif
    flush (output_unit)
    call tell_log(0,'Read BRDF Rsfc from GLER')
  endif

  !================================
  ! 4. read radiance LUT at 466 nm and 440nm
  !================================
  status=GetConfigString("E","Input Files LUT_dir",buf)
  name_lut_dir = trim(adjustl(buf))

  status=GetConfigString("E","Input Files LUT_RAD_466",buf)
  name_lut_rad = trim(adjustl(buf))

  status=GetConfigString("E","Input Files LUT_RAD_440",buf)
  name_lut_rad440 = trim(adjustl(buf))

  call read_lut_rad (errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_lut_rad failed', errstat)
    call exit(-1)
  endif
  flush (output_unit)
  call tell_log(0,'Read radiance look-up table at 466 nm')

  call read_lut_rad440 (errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_lut_rad440 failed', errstat)
    call exit(-1)
  endif
  flush (output_unit)
  call tell_log(0,'Read radiance look-up table at 440 nm')

  !===========================
  ! 5. read AMF LUT at 477 nm
  !===========================
  call read_lut_amf_clr (errstat)
  call read_lut_amf_cld (errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'error reading AMF LUTs failed', errstat)
    call exit(-1)
  endif
  flush (output_unit)
  call tell_log(0,'Read AMF look-up table')

  !================================
  ! 6. Allocate arrays for ECF&OCP
  !================================
  nx=rad_nXtrack
  nt=rad_NumTimes
  call allocate_ecf_arrays(nx,nt,fspecial,errstat)
  call allocate_ocp_arrays(nx,nt,fspecial,errstat)

  if (errstat /= 0) then
     call tell_error(tell_runtime_error, 'error allocating ECFOCP failed',errstat)
     call exit(-1)
  endif
  flush (output_unit)
  call tell_log(0,'Allocate ECFOCP arrays')

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ecf-ocp iteration
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do iter_ecfocp = 1, ecfocp_maxiter
     write(*,*) '--------ECFOCP iteration#',iter_ecfocp
     !================================
     ! 6. calculate ECF/CRF at 466 nm
     !================================
     call cal_ecf(iter_ecfocp)
     flush (output_unit)
     call tell_log(0,'   Calculated effective cloud fraction')

     !==================
     ! 7. calculate OCP
     !==================
     call cal_ocp(iter_ecfocp)
     flush (output_unit)
     call tell_log(0,'   Calculated cloud pressure')
  enddo 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! end of ecf-ocp iteration
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
  !=============================
  ! 8. calculate Scene Pressure
  !=============================
  call read_lut_amf_ler (errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'error read_lut_amf_ler failed', errstat)
    call exit(-1)
  endif
  flush (output_unit)
  call tell_log(0,'Read AMF LER look-up table')

  call cal_pscene
  flush (output_unit)
  call tell_log(0,'Calculated LER and scene pressure')

  !===================
  ! 9. write outputs
  !===================
  logmsg = 'All calculation is done. Now writing '//trim(name_out_ncdf)
  call tell_log(0, logmsg)

  call write_debug_processing_flags

  if (run_mode .EQ. 'production') then
    ! JCH: in this mode, we add variables to an existing output file
    write(*,*) 'update file with output: '//trim(name_out_ncdf)
    call update_output_file_tio (name_out_ncdf, &
                                 rad_NumTimes, rad_nXtrack, errstat)
  else ! in development mode we create new file
    write(*,*) 'create file for output: '//trim(name_out_ncdf)
    call create_output_file_tio (name_out_ncdf, &
                                 l1radfnm, swathname, &
                                 rad_NumTimes, rad_nXtrack, errstat)
  endif

  if (errstat == 0) then
     call tell_log(0, 'Success! I have done everything.')
  else 
     call tell_log(0, 'Sorry, something is wrong.')
  endif

  call tell_close()
  !*******************
End program OMCDO2N
!*******************
