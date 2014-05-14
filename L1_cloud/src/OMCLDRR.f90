program OMCLDRR

  use m_vars
  use m_rd_chl
  use m_rd_terr
  use m_rd_toms_refl
  use m_initialize
  use m_read_input_data
  use m_read_tables
  use m_read_ocean_table
  use m_write_output_data
  use m_write_output_data_2pres
  use m_cloud_pres_ret
  !use m_cloud_clear_ret
  use m_str_replace
  !use m_num2string
  use m_cloud_mask
  use m_read_thresholds
  use m_read_references
  use m_read_resid
  use m_read_cal
  use MetadataModule
  use m_write_HDFEOS_attr 
  use m_LUN_set
  use L1B_Reader_class
  !use full_reader_class
  !use m_read_no2
  use m_pgs_include

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
  !  28Mar02   Vasilkov   changes to read filename_out from PCF (marked with ****)
  !
  !EOP

  !-------------------------------------------------------------------------
  !include 'PGS_PC.f'
  !include 'PGS_PC_9.f'
  !include 'PGS_OMI_1900.f'
  !include 'PGS_OMCLDRR_52251.f'
  !include 'PGS_SMF.f'

  integer :: i
  integer (kind=4), external :: pgs_pc_getreference !, OMI_SMF_setmsg
  TYPE (L1B_block_type) :: blk
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
  if (using_ref) call read_references(err_code)
  if (using_resid) call read_resids
  if (using_cal) call read_cals(err_code)
  if (do_o3) call read_o3
  !if (do_no2) call read_no2(err_code)

  !Assign name and open output file
  !======================================
  if (ex) then
    nfiles=1
    no_cl_filename=.true.
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

  !Loop over # of files
  !====================
  do ifile=1, nfiles
    if (iprt > 1) print *, 'cloud_ret: processing file ',ifile, ' of ',nfiles

    !if files given in command line, assign names
    !============================================ 
    if (.not. no_cl_filename) then
      filename=filenames(ifile)
      print *,' filename ',filename
      filename_out=trim(output_data_path)//filename
      filename_out = str_replace (filename_out, '.hdf', '')
      filename_out = trim(filename_out)//'_cloudret.hdf'
    endif

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
    call read_input_data(blk, err_code)

    !loop over the # of lines
    !========================
    do while (iLine+ny-1 <= nTimes) 
      if (iprt >= 1) print *,'cloud_ret: processing iLine ',iLine,' of ', &
           nTimes

      do ll=1, ny
        if (ll > 1) iLine=iLine+1

        if(iLine > start_line) then

          !Get level 1 input data: 
          !radiance and solar irradiance spectra
          !===========================================================
          if (iprt > 1) print *,'cloud_ret: reading input data'
          call read_input_data(blk, err_code)
          if(err_code >= 1) goto 999

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
      enddo ! loop over ny

      !do the retrieval
      !=================

      do i=1,n_products
        if(i == 1) then
          refl_clr=0.15d0 ! just for testing, set back later
          refl_cld=0.80d0
        else
          refl_clr=0.11d0
          refl_cld=0.40d0
        endif

        if (.not. get_refl_clim) then
          ref_clr(:,iLine)=refl_clr
        endif

        !if (cloud_clear) then
        ! if (iprt > 1) print *,'cloud_ret: cloud clearing'
        ! call cloud_clear_ret(refl_clr, refl_cld)
        !else
        if (iprt > 1) print *,'cloud_ret: retrieving cloud pressure'
        call cloud_pres_ret(refl_clr, refl_cld)
        !endif

        if (i == 1) then
          cld_pres2(:,iLine)=cloud_pres(:,iLine)
          eff_cld_frac2(:,iLine)=eff_cld_frac(:,iLine)
          qc2(:,iLine)=qc(:,iLine)
          biases2(:,iLine)=biases(:,iLine)
          stds2(:,iLine)=stds(:,iLine)
          chi_sqr2(:,iLine)=chi_sqr(:,iLine)
          shifts2(:,iLine)=shifts(:,iLine)
          !      qc=0
          !      where (btest(qc2,14)) qc=IBSET(qc,14) 
        endif
        if(i == 2 .and. iLine == start_line) then
          qc=0
          where (btest(qc2,14)) qc=IBSET(qc,14)
        endif

      enddo

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


999   continue

      iLine=iLine+1
    enddo ! iLine loop

    !Write the output 
    !=================
    if (write_he5) then
      call write_output_data(filename_out,outswathname)
      call write_output_data_2pres(filename_out,outswathname)
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

  enddo ! file loop


  !exit with normal status
  !=======================
  status = omi_smf_setmsg(OMI_S_SUCCESS, &
       'PGE finishes normally, exit code = 0  ', 'cloud_ret',0)
  status = OMI_S_SUCCESS                                                      
  call exit(0)


END program OMCLDRR

