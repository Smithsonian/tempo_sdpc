MODULE swathline_loop
CONTAINS
SUBROUTINE swathline_loops (                               &
    pge_idx, rpt, n_max_rspec, do_process_line,                     &
    xtrange, do_remove_target, ntargpol,         &
    in_common_mode_loop, errstat, retrieval_opt)

  USE OMSAO_precision_module,  ONLY: i4, r8, i2, r4
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval, MAX_STR_LEN, &
    nlines_max, nUTCdim, NXTRACK_MAX
  USE OMSAO_indices_module,    ONLY: n_max_fitpars
  USE OMSAO_variables_module,  ONLY:  &
    n_fitvar_rad, l1b_rad_filename, verb_thresh_lev, n_fincol_idx, fincol_idx, &
    n_rad_wvl, n_rad_wvl_max, Radiance_Paras_Type, &
    fitvar_rad_init, fitvar_rad_saved
  use ctrlvars, only: yn_radiance_reference, yn_diagnostic_run
  USE OMSAO_omidata_module,    ONLY:  &
    omi_blockline_no,                  &
    omi_itnum_flag, omi_fitconv_flag, omi_column_amount,                     &
    omi_column_uncert, omi_time_utc, omi_time, omi_fit_rms,    &
    omi_radiance_errstat,  &
    omi_szenith, omi_vzenith, omi_latitude, omi_longitude, omi_xtrflg, omi_height, &
    retrieval_type
  USE OMSAO_prefitcol_module, ONLY: read_prefit_columns, init_prefit_files
  USE OMSAO_errstat_module
  USE OMSAO_radiance_ref_module, ONLY: remove_target_from_radiance
  USE omi_read_l1b_data, ONLY: omi_read_radiance_lines
  USE fitting_loops, ONLY: xtrack_radiance_fitting_loop
  USE omi_pge_fitting_aux, ONLY: convert_tai_to_utc
  USE he5_output_tools, ONLY: he5_write_radfit_output
  use errormodule
  use tell_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ntargpol, n_max_rspec
  TYPE (Radiance_Paras_Type), INTENT(IN) :: rpt
  INTEGER (KIND=i4), DIMENSION (0:rpt%ntimes-1,1:2),  INTENT (IN) :: xtrange
  LOGICAL,           DIMENSION (0:rpt%ntimes-1),      INTENT (IN) :: do_process_line
  LOGICAL,           INTENT (IN) :: in_common_mode_loop, do_remove_target
  type (retrieval_type), optional, intent(inout) ::retrieval_opt

  ! ------------------
  ! Modified variables
  ! ------------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)      :: iline, iloop, nblock, fpix, lpix, ipix, estat, locerrstat
  CHARACTER (LEN=MAX_STR_LEN) :: addmsg
  INTEGER (KIND=i4) :: nt, nx, nccd, scanline_no
  integer, parameter :: unit_column_amount=22

  ! ---------------------------------------------------------------
  ! Variables to remove target gas from radiance reference spectrum
  ! ---------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (n_fincol_idx,1:rpt%nxtrack) :: target_var, targsum, targcnt
  REAL (KIND=r8), DIMENSION (1:rpt%nxtrack)              :: target_fit, target_col

  ! ---------------------------------------------------------------------------------
  ! CCM Array to hold (1) Fitted Spec (2) Observed Spec (3) Spec Pos (4) Weight flags
  ! ---------------------------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (n_rad_wvl_max,nxtrack_max,4) :: fitspc_tmp
  REAL (KIND=r8), DIMENSION (:,:,:,:), allocatable :: omi_fitspc

  ! -------------------------------------
  ! Correlations with main output product
  ! -------------------------------------
  REAL (KIND=r8), DIMENSION (n_fitvar_rad,rpt%nxtrack,0:nlines_max-1) :: &
    all_fitted_columns, all_fitted_errors, correlation_columns

  if (errstat < 0) return

  locerrstat = pge_errstat_ok
  nt = rpt%ntimes
  nx = rpt%nxtrack
  nccd = rpt%nwavel_ccd

  ! --------------------------------
  ! Initialize fitting output arrays
  ! --------------------------------
  all_fitted_columns  = r8_missval
  all_fitted_errors   = r8_missval
  correlation_columns = r8_missval
  fitspc_tmp          = r8_missval

  IF ( yn_radiance_reference .AND. do_remove_target ) THEN
    target_var = 0.0_r8
    targsum    = 0.0_r8
    targcnt    = 0.0_r8
    target_fit = 0.0_r8
    target_col = 0.0_r8
  END IF

  if (.not.in_common_mode_loop) then
    allocate (omi_fitspc(n_rad_wvl_max,nxtrack_max,4,0:nlines_max-1), stat=locerrstat)
    if (locerrstat /= 0) then
      errstat = -1
      call err_message_error ("swathline_loops: allocate failed", &
                              errstat)
      return
    endif
    omi_fitspc = 0.0_r8
  endif

  if (yn_diagnostic_run) then
    open (unit=unit_column_amount, file='diag.column_amount', iostat=locerrstat)
      if (locerrstat /= 0) then
        call tell_error (tell_io_open_error, &
                         "error opening diag.column_amount", errstat)
        return
      endif
  endif

  ! ---------------------------------------------------------------------
  ! Loop over all scan lines, in multiples of NLINES_MAX (100 by default)
  ! ---------------------------------------------------------------------
  ScanLines: DO iline = 0, nt-1, nlines_max

    ! ---------------------------------------------------------
    ! Check if loop ends before n_times_loop max is exhausted.
    ! Not a serious problem but it saves a few bytes of memory.
    ! ---------------------------------------------------------
    nblock = nlines_max
    IF ( (iline+nblock) > nt ) nblock = nt - iline
    ! -----------------------------------------
    ! Skip if we don't have anything to process
    ! -----------------------------------------
    IF ( .NOT. ( ANY ( do_process_line(iline:iline+nblock-1) ) ) ) CYCLE

    ! ------------------------------
    ! Get NBLOCK radiance lines
    ! ------------------------------
    write(*,*)'swathline_loops calling omi_read_radiance_lines, iline=',iline
    CALL omi_read_radiance_lines (                   &
      l1b_rad_filename, iline, nx, nblock, nccd, locerrstat )
    ! -----------------------------------------------------------------------------------

    ! ------------------------------------------
    ! Initialize output fields with MissingValue
    ! ------------------------------------------
    omi_itnum_flag   (1:nx,     0:nblock-1) = i2_missval
    omi_fitconv_flag (1:nx,     0:nblock-1) = i2_missval
    omi_column_amount(1:nx,     0:nblock-1) = r8_missval
    omi_column_uncert(1:nx,     0:nblock-1) = r8_missval
    omi_fit_rms      (1:nx,     0:nblock-1) = r8_missval
    omi_time_utc     (1:nUTCdim,0:nblock-1) = i2_missval

    ! --------------------------------
    ! Read pre-fitted molecule columns
    ! --------------------------------
    IF (.NOT. yn_radiance_reference) then
      CALL read_prefit_columns ( pge_idx, nx, nblock, iline, locerrstat )
      errstat = MAX ( errstat, locerrstat )
      IF ( errstat >= pge_errstat_error ) RETURN
    END IF

    ! -------------------------------------
    ! Re-initialize saved fitting variables
    ! -------------------------------------
    fitvar_rad_saved(1:n_max_fitpars ) = fitvar_rad_init(1:n_max_fitpars)

    ! -----------------------------------------------
    ! Loops over all scan lines in current data block
    ! -----------------------------------------------
    ScanLineBlock: DO iloop = 0, nblock-1

      ! --------------------------------------------------------------------
      ! Further down, in deeper layers of the algorithm, we require both the
      ! current line in the data block and the absolute swath line number.
      ! Both values are initialized here.
      ! --------------------------------------------------------------------
      omi_blockline_no = iloop
      scanline_no  = iline+iloop

      IF (scanline_no > nt-1 ) EXIT ScanLines

      ! ----------------------------------------------------------
      ! Skip this line if it isn't in the list of those to process
      ! ----------------------------------------------------------
      IF ( .NOT. do_process_line(scanline_no) ) CYCLE

      ! ------------------
      ! Report on progress
      ! ------------------
      addmsg = ''
      WRITE (addmsg,'(A,I5)') 'Working on scan line', scanline_no
      estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), " ", vb_lev_omidebug )
      !IF ( verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

      IF ( omi_radiance_errstat(iloop) /= pge_errstat_error ) THEN

        fpix = xtrange(scanline_no,1)
        lpix = xtrange(scanline_no,2)

        ! One side effect of this routine is that the value of n_rad_wvl
        ! will change.
        CALL xtrack_radiance_fitting_loop ( &
          pge_idx, n_max_rspec, fpix, lpix, iloop,               &
          n_fitvar_rad,                              &
          all_fitted_columns (1:n_fitvar_rad,fpix:lpix,iloop),   & !gga (1:nx to fpix:lpix)
          all_fitted_errors  (1:n_fitvar_rad,fpix:lpix,iloop),   & !gga (1:nx to fpix:lpix)
          correlation_columns(1:n_fitvar_rad,fpix:lpix,iloop),   & !gga (1:nx to fpix:lpix)
          target_var(1:n_fincol_idx,fpix:lpix), locerrstat, &
          fitspc_tmp, n_rad_wvl_max)

        ipix = (fpix+lpix)/2
        addmsg = ''
        WRITE (addmsg,'(I5, 1x, I4, 3(1PE15.5),I5)') scanline_no, ipix, &
          omi_column_amount(ipix, iloop), omi_column_uncert(ipix, iloop), &
          omi_fit_rms   (ipix, iloop), MAX(-1,omi_itnum_flag(ipix, iloop))
        estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(addmsg), " ", vb_lev_omidebug )
        IF ( verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(addmsg)

        if (yn_diagnostic_run) then
          write (unit_column_amount, '(a)')trim(addmsg)
        endif

        ! CCM Add omi_fitspc - Assignment problem - do an inefficient loop for now
        !DO i=1,n_rad_wvl
        !  DO j=1,nxtrack_max
        !    DO k=1,4
        !      omi_fitspc(i,j,k,iloop) = fitspc_tmp(i,j,k)
        !    ENDDO
        !  ENDDO
        !ENDDO

        if (.not.in_common_mode_loop) &
          omi_fitspc(1:n_rad_wvl,:,:,iloop) = fitspc_tmp (1:n_rad_wvl,:,:)

        ! ---------------------------------------------------------------
        ! Add fitted columns for possible removal from radiance reference
        ! ---------------------------------------------------------------
        IF ( yn_radiance_reference .AND. do_remove_target ) THEN
          DO ipix = fpix, lpix
            IF ( &
              ( omi_fitconv_flag (ipix,iloop) > 0_i2       ) .AND. &
              ( omi_column_amount(ipix,iloop) > r8_missval ) .AND. &
              ( omi_column_amount(ipix,iloop) + &
              2.0_r8*omi_column_uncert(ipix,iloop) >= 0.0_r8 )  ) THEN
              targsum(1:n_fincol_idx,ipix) = &
                targsum(1:n_fincol_idx,ipix) + target_var(1:n_fincol_idx,ipix)
              targcnt(1:n_fincol_idx,ipix) = &
                targcnt(1:n_fincol_idx,ipix) + 1.0_r8

              target_col(ipix) = target_col(ipix) + omi_column_amount(ipix,iloop)
            END IF
          END DO
        END IF

        ! -----------------------------------------------------
        ! Optionally, keep the results of the fitting in memory
        ! -----------------------------------------------------
        if (present(retrieval_opt)) then
          retrieval_opt%column_amount(fpix:lpix,scanline_no)      = omi_column_amount(fpix:lpix,iloop)
          retrieval_opt%column_uncertainty(fpix:lpix,scanline_no) = omi_column_uncert(fpix:lpix,iloop)
          retrieval_opt%rms(fpix:lpix,scanline_no)                = omi_fit_rms(fpix:lpix,iloop)
          retrieval_opt%latitude(fpix:lpix,scanline_no)           = omi_latitude(fpix:lpix,iloop)
          retrieval_opt%longitude(fpix:lpix,scanline_no)          = omi_longitude(fpix:lpix,iloop)
          retrieval_opt%sza(fpix:lpix,scanline_no)                = omi_szenith(fpix:lpix,iloop)
          retrieval_opt%vza(fpix:lpix,scanline_no)                = omi_vzenith(fpix:lpix,iloop)
          retrieval_opt%fit_flag(fpix:lpix,scanline_no)           = omi_fitconv_flag(fpix:lpix,iloop)
          retrieval_opt%xtr_flag(fpix:lpix,scanline_no)           = omi_xtrflg(fpix:lpix,iloop)
          retrieval_opt%height(fpix:lpix,scanline_no)             = REAL(omi_height(fpix:lpix,iloop), KIND = r4)
        endif
      END IF

      ! -----------------------
      ! Convert TAI to UTC time
      ! -----------------------
      call tell_log (1,' *** skipping call to convert_tai_to_utc()')
      if (.false.) then
      CALL convert_tai_to_utc ( &
        nUTCdim, omi_time(iloop), omi_time_utc(1:nUTCdim,iloop) )
      endif

    END DO ScanLineBlock

    ! ----------------------------------------------------------------
    ! AMF calculation and update of fitting statistics only need to be
    ! done for the final round throught the common mode iteration loop
    ! ----------------------------------------------------------------
    IF ( .NOT. in_common_mode_loop ) THEN

      CALL he5_write_radfit_output (                            &
        pge_idx, iline, nx, nblock, fpix, lpix,              &
        all_fitted_columns (1:n_fitvar_rad,1:nx,0:nblock-1), &
        all_fitted_errors  (1:n_fitvar_rad,1:nx,0:nblock-1), &
        correlation_columns(1:n_fitvar_rad,1:nx,0:nblock-1), &
        omi_fitspc,locerrstat )
      errstat = MAX ( errstat, locerrstat )

    END IF

  END DO ScanLines

  if (yn_diagnostic_run) then
    close (unit_column_amount)
  endif

  ! -----------------------------------------
  ! Remove target gas from radiance reference
  ! -----------------------------------------
  IF ( yn_radiance_reference .AND. do_remove_target ) THEN

      WHERE ( targcnt > 0.0_r8 )
        targsum = targsum / targcnt
      ELSEWHERE
        targsum = r8_missval
      END WHERE

      ! ----------------------------------------------------------------
      ! Removing the target gas from the radiance reference will alter
      ! OMI_RADREF_SPEC (1:NWVL,FPIX:LPIX). This is being passed to the
      ! subroutine via MODULE use rather than through the argument list.
      ! ----------------------------------------------------------------
      CALL remove_target_from_radiance (                              &
        nccd, fpix, lpix, n_fincol_idx, fincol_idx(1:2,1:n_fincol_idx),  &
        ntargpol, targsum(1:n_fincol_idx,fpix:lpix), target_fit(fpix:lpix) )

  END IF

  RETURN

END SUBROUTINE swathline_loops
END MODULE

