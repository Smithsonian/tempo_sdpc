! Jbak
! l2_hdf_flag = 0 : ascii
!               3 : he5, only for OMI
!               4 : NC

MODULE O3P_output_module
  USE OMSAO_indices_module, ONLY:instrument_idx, omi_idx
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: l2_filename, l2funit,  & 
      num_param, num_wav_max, & 
      num_wav_max, n_fitvar_rad, numwin, &
      linenum_lim, l2_hdf_flag, use_backup
  USE ozprof_data_module, ONLY: nfgas, nlay,  ozfit_start_index, ozfit_end_index
  USE ascii_output_module, only: omi_write_intermed, &
                                l2_ascii_create, l2_ascii_close
  USE he5_output_module, only: he5_l2setgeofields, he5_l2setdatafields,he5_l2wrtinit
  USE tio_output_module
  USE OMI_metaData_class
  USE PROFOZ_metaDef
  USE OMI_LUN_SET
  CONTAINS

  SUBROUTINE L2_O3P_CREATE (ntimes,first_pix, last_pix, first_line, last_line, errstat)
   IMPLICIT NONE
   !--------------------
   ! IN/OUT variables
   !--------------------
   INTEGER, INTENT(IN)  :: ntimes, first_pix, last_pix, first_line, last_line
   INTEGER, INTENT(OUT) :: errstat
   !-------------------
   ! Local variables
   !-------------------
   ! variables for just OMI
   !TYPE(OMIECSMETA_T) :: L1BcoreMeta
   !TYPE(ECSMETA_ITEM_T), DIMENSION(6) :: PROFOZ_metaItems !
   !integer (kind=4), dimension(3) :: LUNinputPointer
   !INTEGER :: errstat, i, version, jday, the_year, the_month, the_day
   !character(len=6) :: ShortName = 'PROFOZ'
   !--------------------------------------------------------------
   INTEGER :: i
   INTEGER, DIMENSION(:), ALLOCATABLE :: step_idx

   ! Initialize
   allocate(step_idx(ntimes))
   step_idx(1:ntimes) = (/(i, i= 0, ntimes-1)/)   
   IF (l2_hdf_flag == 0) THEN ! TEXT output
    CALL l2_ascii_create(l2_filename, l2funit, errstat)
   ELSE IF (l2_hdf_flag == 3) THEN  
    CALL He5_L2WrtInit (first_pix, last_pix, first_line, last_line, errstat)
   ELSE IF (l2_hdf_flag == 4) then ! netCDF output
    CALL l2_tio_create(l2_filename, first_pix, last_pix, first_line, &
             last_line, nfgas, nlay, n_fitvar_rad, numwin, num_param, &
             num_wav_max, step_idx(linenum_lim(1):linenum_lim(2)), &
             errstat)
   ENDIF

  END SUBROUTINE L2_O3P_CREATE
   
 
  SUBROUTINE L2_O3P_WRITE_GEO (geo, first_pix, last_pix, first_line, last_line, errstat)
   USE OMSAO_variables_module, ONLY: geo_group
   IMPLICIT NONE
   !--------------------
   ! IN/OUT variables
   !--------------------
   TYPE(geo_group), INTENT(IN) :: geo
   INTEGER, INTENT(IN)  :: first_pix, last_pix, first_line, last_line
   INTEGER, INTENT(OUT) :: errstat
   
   ! Initialize
    errstat = 0
    IF (l2_hdf_flag == 0) THEN ! TEXT output
      RETURN
    ELSE IF (l2_hdf_flag == 3) THEN
      CALL He5_L2SetGeoFields (geo,first_pix, last_pix, errstat )          
    ELSE IF (l2_hdf_flag == 4) then ! netCDF output
      CALL l2_tio_write_geo(geo,first_pix, last_pix, first_line, last_line, &
             errstat)
    ENDIF
   RETURN  
  END SUBROUTINE L2_O3P_WRITE_GEO 

  SUBROUTINE L2_O3P_write_data (currpix, first_pix, last_pix, & 
                                currloop,currline, ntimes_loop, & 
                                exval, fitcol, dfitcol, message, problems)
   IMPLICIT NONE
   !--------------------
   ! IN/OUT variables
   !--------------------
   INTEGER, INTENT(IN)  :: currpix, first_pix, last_pix, currloop, currline, ntimes_loop, &
                           exval
   REAL (KIND=dp), DIMENSION(3) :: fitcol
   REAL (KIND=dp), DIMENSION(3, 2) :: dfitcol
   LOGICAL, INTENT(OUT) :: problems
   CHARACTER (100), INTENT(OUT) :: message
   !-------------------
   ! Local variables
   !-------------------
   INTEGER :: errstat, ix, iy
   errstat = 0
   IF (l2_hdf_flag == 0) THEN 
     IF (exval > 0) THEN 
         call omi_write_intermed (l2funit, fitcol, dfitcol, exval) 
     ENDIF
   ELSE IF (l2_hdf_flag == 3) THEN 
      call He5_L2SetDataFields (currpix, first_pix, last_pix, &
           currloop, currline, ntimes_loop, exval, fitcol, dfitcol, &
           errstat)
          IF ( errstat /= 0  ) THEN
              message = ':falied to write l2 o3p (he5)'
              problems=.true. ; return
          ENDIF                           
   ELSE IF (l2_hdf_flag == 4) THEN 
     ix = currpix - first_pix ! start from zero for l2_tio_write_data
     iy = currline
   
     IF (exval >= 0) then  ! Retrieval finished.
        call l2_tio_write_data (ix, iy, &
        exval, fitcol, dfitcol, nfgas, nlay, n_fitvar_rad, &
        numwin, num_param, num_wav_max, ozfit_start_index, &
        ozfit_end_index, errstat)
        if (errstat < 0) then
           message=": L2 write failed (nc)"
           problems=.true. ; return
        ENDIF
     ELSE  ! Retrieval failed. Fill in as missing values.
        call l2_tio_fill_data (ix, iy, &
        exval, nfgas, nlay, n_fitvar_rad, numwin, num_param, &
        num_wav_max, ozfit_start_index, ozfit_end_index, errstat)
        if (errstat < 0) then
           message=": L2 write failed (nc)"
           problems=.true. ; return
        endif
     ENDIF
   ENDIF
  END SUBROUTINE

  SUBROUTINE  L2_o3p_close (errstat)
  INTEGER, INTENT(OUT) :: errstat
  errstat = 0
  IF (l2_hdf_flag == 0) THEN
    CALL l2_ascii_close(l2funit, errstat)
  ELSE IF (l2_hdf_flag==3) THEN
  ELSE IF (l2_hdf_flag==4) THEN
    CALL l2_tio_close(errstat)
  ENDIF
  END SUBROUTINE L2_o3p_close 
END MODULE


