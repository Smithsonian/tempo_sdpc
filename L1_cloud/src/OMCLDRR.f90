program OMCLDRR

  use m_vars
  use m_rd_chl
  use m_rd_terr
  use m_rd_toms_refl
  use m_initialize
  use m_read_input_data
  use m_read_input_data_tio
  use m_read_tables
  use m_read_ocean_table
  use m_write_output_data
  use m_write_output_data_2pres
  use m_write_output_data_tio
  use m_cloud_pres_ret
  use m_str_replace
  use m_cloud_mask
  use m_read_thresholds
  use m_read_resid
  use MetadataModule
  use m_write_HDFEOS_attr 
  use m_LUN_set
  use L1B_Reader_class
  use m_pgs_include
  use tell_module

  IMPLICIT NONE

  !-------------------------------------------------------------------------
  !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
  !-------------------------------------------------------------------------
  !BOP
  !
  ! !PROGRAM:  cloud_ret
  ! 
  ! !DESCRIPTION: main driver for OMI cloud product retrieval
  !		
  ! !SEE ALSO:  
  !
  ! !REVISION HISTORY: 
  !
  !  05Jan01   Joiner     original fortran 90
  !  28Mar02   Vasilkov   changes to read filename_out from PCF (marked ****)
  !  26Aug14   O'Sullivan tidying, simplification, annotation for TEMPO
  !
  !EOP

  !-------------------------------------------------------------------------
  !include 'PGS_PC.f'
  !include 'PGS_PC_9.f'
  !include 'PGS_OMI_1900.f'
  !include 'PGS_OMCLDRR_52251.f'
  !include 'PGS_SMF.f'

  integer (kind=4), external :: pgs_pc_getreference !, OMI_SMF_setmsg
  TYPE (L1B_block_type) :: blk

  integer (kind=4) :: ext_index, errstat
  character(len=255) :: filename_out_nc, filename_in_nc

  !************************************************************************

  iLine=0
  !Initialize (read resource file)
  !===============================
  if (iprt > 1) print *,'cloud_ret: initializing'
  call initialize(err_code)


  !Read in pre-computed Ring and radiance data
  !===========================================
  if (iprt > 1) print *,'cloud_ret: reading_tables'
  call read_tables(err_code)
  call read_ocean_table(err_code)
  call read_thresholds(err_code)
  if (using_resid) call read_resids
  if (do_o3) call read_o3

  !Assign name and open output file
  !======================================
  if (ex) then
    version = 1
    status = pgs_pc_getreference(L2_out,version,flnm_out)
    if(status.ne.0) then
      ierr=OMI_SMF_setmsg(OMCLDRR_F_FAILURE, &
           "error opening output L2 file, PGE aborting, exit code = 1", &
           "program cloud_ret",0)
      stop 1
    endif
    filename_out=trim(flnm_out)
    if (iprt > 1) print *,'cloud_ret: status',status,' output filename ',&
         trim(flnm_out)
  endif

  ! Main processing stage
  !============================================================


  !Get level 1 small pixel data and
  !compute cloud mask
  !===========================================================
  if (iprt >= 1) print *,'cloud_ret: calling cld_mask'
  call cld_mask()
  if (iprt >= 1) print *,'cloud_ret: finished calling cld_mask'

  !Get level 1 input data: 
  !radiance and solar irradiance spectra
  !===========================================================
  if (iprt > 1) print *,'cloud_ret: reading input data'
  !he5 version
!  call read_input_data(blk, err_code)
  !netCDF version
  ext_index=index(filename,'.he4')
 print *, filename
  filename_in_nc=filename(1:ext_index-1)//'.nc'
  call read_input_data_tio(filename_in_nc, errstat)
call exit(0)

  !loop over the # of lines
  !========================
  do while (iLine <= nTimes) 
    if (iprt >= 1) print *,'cloud_ret: processing iLine ',iLine,' of ', &
         nTimes


    if(iLine > start_line) then

      !Get level 1 input data: 
      !radiance and solar irradiance spectra
      !===========================================================
      if (iprt > 1) print *,'cloud_ret: reading input data'
      !he5 version
      call read_input_data(blk, err_code)
      if(err_code >= 1) goto 999
      !netCDF version
      call read_input_data_tio(filename_in_nc, errstat)
      if(errstat > 0) goto 999
      
      ext_index=index(filename_out,'.he5')
      filename_in_nc=filename_out(1:ext_index-1)//'.nc'

    endif ! start_line

    !get the climatological terrain pressure
    !=======================================
    call rd_terr ()

    !get the surface reflectivity climatology 
    !=======================================
    if (get_refl_clim) then
      call rd_toms_refl ()
    endif

    !get the climatological chlorophyll
    !=======================================
    call rd_chl ()
    land_flg=chlcl < 0.

    !do the retrieval
    !=================

    ! Note that the clear and cloudy reflectance values below (0.15/0.8)
    ! over-ride those in m_vars, and disagree with Joiner & Vasilkov (2006).
    ! However these values are required to produce results consistent 
    ! with the OMI pipeline. refl_clr will be replaced by TOMS reflect.
    ! climatology if get_refl_clim=.true. (as it usually is).
    refl_clr=0.15d0 
    refl_cld=0.80d0 


    if (.not. get_refl_clim) then
      ref_clr(:,iLine)=refl_clr
    endif

    if (iprt > 1) print *,'cloud_ret: retrieving cloud pressure'
    call cloud_pres_ret(refl_clr, refl_cld)

    cld_pres2(:,iLine)=cloud_pres(:,iLine)
    eff_cld_frac2(:,iLine)=eff_cld_frac(:,iLine)
    qc2(:,iLine)=qc(:,iLine)
    biases2(:,iLine)=biases(:,iLine)
    stds2(:,iLine)=stds(:,iLine)
    chi_sqr2(:,iLine)=chi_sqr(:,iLine)
    shifts2(:,iLine)=shifts(:,iLine)

    !write output 
    !============
    if (iprt >= 2) then
      print *, 'cloud_ret: pix, CP, R, f, ps, sza, land, biases, stds, chl. fg, chl, ref. clr.'
      do ispec=0, nXtrack-1
        write(6,"(i5,5f8.3,l3,5e11.3,f8.3)") ispec, cloud_pres(ispec,iLine), &
             refl(ispec,iLine), eff_cld_frac(ispec,iLine), ps(ispec,iLine), &
             sza(ispec,iLine), land_flg(ispec), biases(ispec,iLine), &
             stds(ispec,iLine), chlcl(ispec), chlorophyll(ispec,iLine), &
             ref_clr(ispec,iLine)
      enddo ! ispec
    endif ! iprt >= 2


999 continue

    iLine=iLine+1
  enddo ! iLine loop

  !Write an output .he5 file
  !=========================
  call write_output_data(filename_out,outswathname)
  call write_output_data_2pres(filename_out,outswathname)

  !Write out an output .nc file
  !============================
  ext_index=index(filename_out,'.he5')
  filename_out_nc=filename_out(1:ext_index-1)//'.nc'
  if (write_resid) then
    call create_output_file(filename_out_nc,nTimes,nXtrack,err_code, &
         size(wave_resid))
  else
    call create_output_file(filename_out_nc,nTimes,nXtrack,err_code)
  endif
  if (err_code < 0) then
    call tell_error (tell_io_write_error, &
           "create_output_file: failed", &
           err_code)
  endif
  call close_output_file(err_code)
  if (err_code < 0) then
    call tell_error (tell_io_write_error, &
           "close_output_file: failed", &
           err_code)
  endif

  !Writing Metadata including LocalGranuleId
  !=========================================
  call write_HDFEOS_attr(filename_out,outswathname,err_code)
  retstatus = RdWrMetadata(filename_out)

  ! close data block structure
  status = L1Br_close( blk )
  IF( status .NE. OMI_S_SUCCESS ) THEN
    ierr = OMI_SMF_setmsg( status, &
         "L1Br_close failed, PGE aborting, exit code = 1", &
         "OMCLDRR_2pres", 1 )
    call exit(1)
  END IF



  !exit with normal status
  !=======================
  status = omi_smf_setmsg(OMI_S_SUCCESS, &
       'PGE finishes normally, exit code = 0  ', 'cloud_ret',0)
  status = OMI_S_SUCCESS                                                      
  call exit(0)


END program OMCLDRR

