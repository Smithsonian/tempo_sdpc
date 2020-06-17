!>Retrieve cloud parameters via measurements of reflectance and 
!>Raman scattering
program OMCLDRR

  use m_vars
  use m_rd_chl
  use m_rd_terr
  use m_rd_toms_refl
  use m_initialize
  use m_read_input_data
  use m_read_input_data_tio
  use m_read_tables_tio
  use m_read_ocean_table_tio
  use m_write_output_data
  use m_write_output_data_2pres
  use m_write_output_data_tio
  use m_cloud_pres_ret
  use m_cloud_mask
  use m_read_thresholds
  use m_read_resid
  use MetadataModule
  use m_write_HDFEOS_attr 
  use m_LUN_set
  use L1B_Reader_class
  use m_pgs_include
  use tell_module
  use m_write_odl_metadata
  use md_module
  use tio_module, only : tempo_prod_type_cldrr

  IMPLICIT NONE

  !-------------------------------------------------------------------------
  ! !REVISION HISTORY:
  !
  !  05Jan01   Joiner     original fortran 90
  !  28Mar02   Vasilkov   changes to read filename_out from PCF (marked ****)
  !  26Aug14   O'Sullivan tidying, simplification, annotation for TEMPO
  !-------------------------------------------------------------------------

  integer (kind=4), external :: pgs_pc_getreference
  !>@param blk Defined type for passing L1B he4 input data from
  TYPE (L1B_block_type) :: blk

  !>@param errstat integer error code used throughout program
  integer (kind=4) :: errstat
  !>@param filename_out_nc netCDF output filename
  !>@param filename_in_nc netCDF input filename
  character(len=255) :: logmsg
  integer (kind=4) :: n, j, returnstatus, dummy_version
  integer (kind=4), parameter :: ninp=2 ! number of input files
  integer (kind=4), dimension(ninp) :: input_luns
  character (len=128), dimension(ninp) :: inputs
  character (len=256) :: buf

  type (boundary_polygon_type) :: bdry
  !************************************************************************

  errstat=0
  call tell_open ("L1_cloud", 0)

  iLine=0
  !Initialize (read resource file)
  !===============================
  call tell_log(2,'initializing')
  call initialize(errstat)


  !Read in pre-computed Ring and radiance data
  !===========================================
  call tell_log(2,'reading_tables')
  call read_tables_tio(errstat)
  call read_ocean_table_tio(errstat)
  if (errstat /= 0) then
    call tell_error (tell_io_error, &
         "Failed to read Ring-effect tables", &
         errstat)
    stop 1
  endif

  call read_thresholds(errstat)
  if (using_resid) call read_resids(errstat)
  if (do_o3) call read_o3(errstat)
  if (errstat /= 0) then
    call tell_error (tell_io_error, &
         "Failure in reading thresholds, residuals, or O3 cross-section", &
         errstat)
    stop 1
  endif


  ! Main processing stage
  !============================================================


  !Get level 1 small pixel data and
  !compute cloud mask
  !===========================================================
  if (.not. read_he4) do_cloud_mask = .false.
  if (do_cloud_mask) then
    call tell_log(1,'calling cld_mask')
    call cld_mask(errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "Failure to generate cloud mask", &
           errstat)
      stop 1
    else
      call tell_log(2,'finished calling cld_mask')
    endif
  else
    call tell_log(2,'skipping cloud mask')
  endif

  !Get level 1 input data: 
  !radiance and solar irradiance spectra
  !===========================================================
  call tell_log(2,'reading input data')
  !he5 version
  if (read_he4) call read_input_data(blk, errstat)
  !netCDF version
  if (read_nc) then
!    ext_index=index(filename,'.he4')
!    filename_in_nc=filename(1:ext_index-1)//'.nc'
    if (read_he4) iLine=iLine-1 
    call read_input_data_tio(filename_in_nc, errstat)
  endif
  if (errstat /= 0) then
    call tell_error (tell_io_error, &
         "Failed to input first line of data, exiting", &
         errstat)
    stop 1
  endif

  !loop over the # of lines
  !========================
  do while (iLine <= nTimes) 
    write(logmsg,"(A17,I4,A4,I4)") 'processing iLine ',iLine,' of ',nTimes
    call tell_log(1,logmsg)


    if(iLine > start_line) then

      !Get level 1 input data: 
      !radiance and solar irradiance spectra
      !===========================================================
      call tell_log(2,'reading input data')
      !he5 version
      if (read_he4) then
        call read_input_data(blk, errstat)
        if(errstat /= 0) goto 999
      endif
      !netCDF version
      if (read_nc) then
        call read_input_data_tio(filename_in_nc, errstat)
        if(errstat /= 0) goto 999
      endif

    endif ! start_line

    !get the climatological terrain pressure
    !=======================================
    ! Crude presssure climatology - higher resolution values could
    ! be read in from meteorological forecast
    call rd_terr (errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "Error opening terrain height file, exiting", &
           errstat)
      stop 1
    endif

    !get the surface reflectivity climatology 
    !=======================================
    if (get_refl_clim) then
      call rd_toms_refl (errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_error, &
             "Error opening surface reflectivity file, exiting", &
             errstat)
        stop 1
      endif
    endif

    !get the climatological chlorophyll
    !=======================================
    call rd_chl (errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "Error opening chlorophyll file, exiting", &
           errstat)
      stop 1
    endif
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
      ref_clr(:,iLine)=real(refl_clr, kind=4)
    endif

    call tell_log(2,'retrieving cloud pressure')
    call cloud_pres_ret(refl_clr, refl_cld, errstat)
    if (errstat /= 0) then
      call tell_error (tell_application_error, &
           "cloud_pres_ret: failed, exiting", &
           errstat)
      stop 1
    endif
 
    cld_pres2(:,iLine)=cloud_pres(:,iLine)
    eff_cld_frac2(:,iLine)=eff_cld_frac(:,iLine)
    qc2(:,iLine)=qc(:,iLine)
    biases2(:,iLine)=biases(:,iLine)
    stds2(:,iLine)=stds(:,iLine)
    chi_sqr2(:,iLine)=chi_sqr(:,iLine)
    shifts2(:,iLine)=shifts(:,iLine)

!    !write output 
!    !============
!    if (iprt >= 2) then
!      print *, 'cloud_ret: pix, CP, R, f, ps, sza, land, biases, stds, chl. fg, chl, ref. clr.'
!      do ispec=0, nXtrack-1
!        write(6,"(i5,5f8.3,l3,5e11.3,f8.3)") ispec, cloud_pres(ispec,iLine), &
!             refl(ispec,iLine), eff_cld_frac(ispec,iLine), ps(ispec,iLine), &
!             sza(ispec,iLine), land_flg(ispec), biases(ispec,iLine), &
!             stds(ispec,iLine), chlcl(ispec), chlorophyll(ispec,iLine), &
!             ref_clr(ispec,iLine)
!      enddo ! ispec
!    endif ! iprt >= 2


999 continue

    iLine=iLine+1
  enddo ! iLine loop

  !Write an output .he5 file
  !=========================
  if (write_he5) then
    call write_output_data(filename_out,outswathname)
    call write_output_data_2pres(filename_out,outswathname)
  endif

  !Write out an output .nc file
  !============================
  if (write_nc) then
!    ext_index=index(filename_out,'.he5')
!    filename_out_nc=filename_out(1:ext_index-1)//'.nc'
    if (write_resid) then
      call create_output_file(filename_out_nc,nTimes,nXtrack,errstat, &
           size(wave_resid))
    else
      call create_output_file(filename_out_nc,nTimes,nXtrack,errstat)
    endif
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "create_output_file: failed", &
           errstat)
    endif
    call copy_pixel_corners (filename_in_nc, nTimes, nXtrack, errstat)
    call copy_hdr_metadata (filename_in_nc, errstat)
    call label_output_file (tempo_prod_type_cldrr, processing_version, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, "failed writing metadata", &
           errstat)
    endif
    ! Even if create_output_file fails, try to close file to make sure we
    ! don't exit with a half-written file open.
    call close_output_file(errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
           "close_output_file: failed", &
           errstat)
      stop 1
    endif

    ! Proof-of-concept example of ODL ASCII and netCDF metadata
    if (wrt_odl) then
      !input_luns=(/ L1B_LUN, IRR1B_file, terr_prs_id, chl_id, oc_ram_id, &
      !     ring_id, thresh_id, resid_id_late, refl_id /)
      input_luns=(/ L1B_LUN, IRR1B_file /)
      do n=1,ninp
        dummy_version = 1
        returnstatus = PGS_PC_GetReference( input_luns(n), dummy_version, buf )
        if( returnstatus /= 0 ) then
          call tell_error(tell_io_read_error, &
               "write_odl_metadata: failed to read input file names", errstat)
          stop 1
        else
          j = index( buf, '/', BACK = .true. ) + 1
          inputs(n) = trim(buf(j:))
        endif
      enddo

      call make_bounding_polygon (bdry, errstat)

      call md_open (filename_out_nc, errstat)
      call md_write_geo_bounds (bdry % lons, bdry % lats, errstat)
      call md_write_inputs (ninp, inputs, errstat)
      call md_write_prodid (filename_out_nc, processing_version, errstat)
      call md_close (errstat)

      errstat = write_odl_metadata (filename_out_nc, bdry)
      if (errstat /= 0) then
        call tell_error(tell_io_write_error, "failed writing ODL metadata", &
             errstat)
        stop 1
      endif
    endif
    !
    call tell_log(1,'netCDF file output successfully')
  endif

  !Writing Metadata including LocalGranuleId
  !=========================================
  if (read_he4) then
    call write_HDFEOS_attr(filename_out,outswathname,err_code)
    retstatus = RdWrMetadata(filename_out)

    ! close data block structure
    status = L1Br_close( blk )
    IF( status .NE. OMI_S_SUCCESS ) THEN
      call tell_error(tell_io_error, &
           "L1Br_close failed, L1_cloud aborting, exit code = 1", &
           errstat)
      stop 1
    END IF
  endif


  !exit with normal status
  !=======================
  call tell_log(1,"L1_cloud finished normally")

  call tell_close()

  stop


END program OMCLDRR

