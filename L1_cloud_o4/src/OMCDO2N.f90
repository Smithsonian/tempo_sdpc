!***************
program OMCDO2N
!***************

  use iso_fortran_env, only: output_unit
  use hdfeos4_parameters
  use he5_swreader
  use m_vars
  use m_read_input_kleipool
  use m_read_GMI
  use m_read_lut
  use m_read_hdf5
  use tell_module
  use tio_module
  use m_read_input_tio, only: read_rad_tio, read_irr_tio, read_cldo4_tio, &
                  get_tio_global_attr, get_tio_l1rad_glbattr
  use m_write_output_tio, only: create_output_file_tio, update_output_file_tio
  use m_read_input_clim, only: read_geoscf
  use m_read_input_gler, only: read_gler
  use m_cal_ecf, only : cal_ecf
  use m_cal_ocp, only : cal_ocp
  use m_cal_pscene, only : cal_pscene
  use HDF5

  implicit none

  !============================
  ! read control.txt in ../utd
  !============================
  include 'GetConfig.inc'
  !include "Messages.inc"

  character(len=CFG_VAL_LEN)::buf
  integer::gmonth
  integer::status

  !character(len=255)::outfile
  character(len=255)::filename, l1radfnm
  character(len=255)::name_gmi_psfc
  character(len=255)::name_gmi_tmp
  integer::id_gmi_psfc,id_gmi_tmp
  character(len=255)::name_kleipool_rsfc

  !integer(kind=4)::ierr
  integer(kind=4) :: errstat
  character(len=80) :: logmsg
  character(len=15), parameter :: swathname = "band_290_490_nm"

  call tell_open ("L1_cloud_o4", 0)
  errstat = 0

  ! JCH: this doesn't seem to work.
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

  if(name_option_SlantColumnDensity.eq.'NASA') then
    status=GetConfigString("E","Input Files TEMPOO4SCD_fnm",buf)
    if(status < 0) then
      call tell_error(tell_io_read_error,"Problem reading OMO4SCD from control file", errstat)
      call exit(-1)
    endif
    name_nasa_file=trim(buf)
  endif

  status=GetConfigString("E","Input Files TEMPOO4SCD_dir",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading OML1BRVG from control file", errstat)
    call exit(-1)
  endif
  name_nasa_dir=trim(buf)

!hqw TEMPO GLER is now handled by m_read_input_gler
!explicit GLER fnm is no longer needed for TEMPO, but is needed for OMI
  if((name_option_SurfaceReflectivity.eq.'BRDF')) then 
    status=GetConfigString("E","Input Files OMGLER",buf)
    if(status < 0) then
      call tell_error(tell_io_read_error,"Problem reading GLER from control file", errstat)
      call exit(-1)
    endif
    name_brdf_file=trim(buf)
  endif

  status=GetConfigString("E","Input Files TEMPOL1IRR_fnm",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading IRR from control file", errstat)
    call exit(-1)
  endif
  name_irr_file=trim(buf)

  status=GetConfigString("E","Input Files TEMPOL1IRR_dir",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading IRR from control file", errstat)
    call exit(-1)
  endif
  name_irr_dir=trim(buf)

  ! -----------------------------------------
  ! READ IN OUTPUT FILENAMES from control.txt
  ! -----------------------------------------
  status=GetConfigString("E","Output Files TEMPOCLDO4",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading OMCDO2N Output Filename from control file", errstat)
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

  status=GetConfigString("E","Runtime Parameters APPShortName",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading APPName from control file", errstat)
    call exit(-1)
  endif
  gmetadata%appshortname=trim(buf)

  status=GetConfigString("E","Runtime Parameters APPVersion",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading parameter APPName from control file", errstat)
    call exit(-1)
  endif
  gmetadata%appversion=trim(buf)

  status=GetConfigString("E","Runtime Parameters AuthorName",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading parameter Author Name from control file", errstat)
    call exit(-1)
  endif
  gmetadata%author_name=trim(buf)

   status=GetConfigString("E","Runtime Parameters AuthorAffiliation",buf)
   if (status <0) then
      call tell_error(tell_io_read_error,"Problem reading Author Affiliation from control file", errstat)
      call exit(-1)
   endif
   gmetadata%author_affiliation=trim(buf)

  status=GetConfigString("E","Runtime Parameters ProcessingCenter",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading parameter Processing Center  from control file", errstat)
    call exit(-1)
  endif
  gmetadata%processingcenter=trim(buf)

  status=GetConfigString("E","Runtime Parameters TEMPO Footprint",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading parameter TEMPO footprint from control file", errstat)
    call exit(-1)
  endif
  gmetadata%omiwindow=trim(buf)

  status=GetConfigString("E","Runtime Parameters Collection",buf)
  if(status < 0) then
    call tell_error(tell_io_read_error,"Problem reading parameter Collection Number from control file", errstat)
    call exit(-1)
  endif
  gmetadata%omi_collection=trim(buf)

  flush (output_unit)
  call tell_log(0,'Read control file')

  !===========================================
  !hqw get global attributes from TEMPO L1 RAD to gmetadata
  !===========================================
  l1radfnm = trim(adjustl(name_rad_dir))//trim(adjustl(name_rad_file))
  call get_tio_l1rad_glbattr(l1radfnm,errstat)
  if (errstat < 0) then
     call tell_error(tell_io_read_error,"Error getting time_coverage_start from TEMPO L1 RAD" ,errstat)
     call exit(-1)
  endif

  !hqw moved 5 to 2.1 because GEOS-CF TP needs rad/lon
  !  and read_irr_tio need out_ProcessingQualityFlags
  !==============================
  ! 5. read inputs from radiance file -> 2.1
  !==============================
  !call read_rad
  ! Ewan Tested read_rad_tio, reports no errors
  call read_rad_tio (l1radfnm, swathname, errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_rad_tio failed', errstat)
    call exit(-1);
  endif
  flush (output_unit)
  call tell_log(0,'Read RAD '//l1radfnm)

  !============================================
  ! 1. read inputs from OML1BIRR on 12/22/2004
  !============================================
  !call read_irr ! replaced with TEMPO subroutine
  ! Ewan TIO irradiance input tested, appears to work OK
  filename = trim(adjustl(name_irr_dir))//trim(adjustl(name_irr_file))
  call read_irr_tio(filename,swathname,errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_irr_tio failed', errstat)
    call exit(-1);
  endif
  flush (output_unit)
  call tell_log(0,'Read IRR '//filename)

  ! --------------------
  ! 2.2 read TEMPO O4 SCD from intermediate L2 file
  ! --------------------
  ! GGA TIO cldo4 input tested, appears to work OK
  filename = trim(adjustl(name_nasa_dir))//trim(adjustl(name_nasa_file))
  call read_cldo4_tio(filename,errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'read_cldo4_tio failed', errstat)
    call exit(-1);
  endif
  flush (output_unit)
  call tell_log(0, 'Read O4 SCD '//filename)

  ! ---------------------
  ! 2.3 T/P/Psfc from GMI
  ! ---------------------
  !hqw added the following for gmonth
  !gmonth is used to decide filename for GMI
  gmonth = gmetadata%granule_month
  write(*,*) '   gmonth=',gmonth
  flush (output_unit)

  if(name_option_TemperaturePressure.eq.'GMI') then
    !hqw moved above read(gmeta_orbit_month,'(i2)') gmonth
    status=GetConfigString("E","Input Files "//lun_gmi_psfc(gmonth),buf)
    name_gmi_psfc=trim(name_gmi_dir)//trim(buf)
    id_gmi_psfc=ilun_gmi_psfc(gmonth)

    status=GetConfigString("E","Input Files "//lun_gmi_tmp(gmonth),buf)
    name_gmi_tmp=trim(name_gmi_dir)//trim(buf)
    id_gmi_tmp=ilun_gmi_tmp(gmonth)

    write(*,*)'   name_gmi_psfc=',trim(name_gmi_psfc)
    write(*,*)'   name_gmi_tmp=',trim(name_gmi_tmp)
    call read_GMI_TMP(name_gmi_psfc,name_gmi_tmp)
    call tell_log(0, "Read GMI T/P & GMI Ps for month ")
  endif

  flush (output_unit)

  !hqw TEMPO operational uses GEOS5 as the option for GEOS-CF
  if(name_option_TemperaturePressure.eq.'GEOS5') then
  !T-P is not needed in cal_ecf, but surface pressure is needed
     call read_geoscf (errstat)
     if (errstat /= 0) then
       call tell_error (tell_runtime_error, 'read_geoscf failed', errstat)
       call exit(-1)
     endif
     flush (output_unit)
     call tell_log(0,'Read GEOS-CF profiles')
  endif

  ! ------------------------------
  ! 2.4 Rsfc from Kleipool or BRDF
  ! ------------------------------
  if(name_option_SurfaceReflectivity.eq.'Kleipool') then
    status=GetConfigString("E","Input Files Kleipool_dir",buf)
    name_kleipool_dir = trim(adjustl(buf))
    status=GetConfigString("E","Input Files Kleipool_fnm",buf)
    name_kleipool_rsfc=trim(name_kleipool_dir)//trim(adjustl(buf))
    call read_Kleipool_Rsfc(name_kleipool_rsfc,gmonth)
    flush (output_unit)
    call tell_log(0,'Read Kleipool Rsfc climatology')
  endif

  !hqw use BRDF option for TEMPO GLER
  if(name_option_SurfaceReflectivity.eq.'BRDF') then
    !call read_BRDF_Rsfc_h5 ! this for OMI, change to
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
  ! 3. read radiance LUT at 466 nm and 440nm
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
  ! 4. read AMF LUT at 477 nm
  !===========================
  call read_lut_amf_clr (errstat)
  call read_lut_amf_cld (errstat)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, 'error reading AMF LUTs failed', errstat)
    call exit(-1)
  endif
  flush (output_unit)
  call tell_log(0,'Read AMF look-up table')

  !hqw moved out_ProcessingQualityFlags to read_rad_tio to do
  ! allocate and initialize quality flags for better memory usage
  !allocate(out_ProcessingQualityFlags(nx,nt),stat=ierr)
  !out_ProcessingQualityFlags=0

  !================================
  ! 6. calculate ECF/CRF at 466 nm
  !================================
  call cal_ecf
  flush (output_unit)
  call tell_log(0,'Calculated effective cloud fraction')

  !==================
  ! 7. calculate OCP
  !==================
  call cal_ocp
  flush (output_unit)
  call tell_log(0,'Calculated cloud pressure')

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
  logmsg = 'All job done. Now writing '//trim(name_out_ncdf)
  call tell_log(0, logmsg)

  !E. O'Sullivan TEMPO-format OUTPUT - only structure & geolocation so far
  !hqw added other components
  if (run_mode .EQ. 'production') then
    ! JCH: in this mode, we add variables to an existing output file
    write(*,*) 'update file with output'
    call update_output_file_tio (name_out_ncdf, &
                                 rad_NumTimes, rad_nXtrack, errstat)
  else
    write(*,*) 'create file for output'
    call create_output_file_tio (name_out_ncdf, &
                                 l1radfnm, swathname, &
                                 rad_NumTimes, rad_nXtrack, errstat)
  endif

  call tell_log(0, 'Success! I have done everything.')
  call tell_close()
  !*******************
End program OMCDO2N
!*******************
