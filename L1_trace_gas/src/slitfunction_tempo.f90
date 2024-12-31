! read and convolve with TEMPO preflight or irradiance file slit function
MODULE slitfunction_tempo
  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_parameters_module, ONLY: nxtrack_max
  IMPLICIT NONE

  public tempo_slitfunc_read, tempo_slitfunc_convolve, &
         solarcal_write_file, solarcal_read_file

  PRIVATE

  ! Mean supergaussian slit function parameters
  real (kind=8) :: mean_asym, mean_hw1e, mean_shape, mean_wl

CONTAINS

  !------------------------------------------------------------------------
  !
  ! @param  errstat  error-tracking variable, non-zero = failure
  !
  ! @author  E.O'Sullivan July 2020
  !------------------------------------------------------------------------
  SUBROUTINE tempo_slitfunc_read ( errstat )

    USE OMSAO_variables_module,  ONLY: omi_slitfunc_fname, l1b_channel, &
         ctrl_fit_winwav_lim
    use omsao_parameters_module, only: N_FIT_WINWAV
    use l1bread_utils, only: lookup_swathname
    use tio_module
    use tell_module
    use netcdf, only: nf90_nowrite, nf90_noerr, nf90_global, nf90_enotatt, &
         nf90_get_att, nf90_inquire_attribute

    IMPLICIT NONE

    !input variables
    integer (kind=4), intent(inout) :: errstat

    ! --------------
    ! Local variable
    ! --------------
    character (len=15) :: swathname
    character (len=16) :: prod_str
    integer (kind=4) :: ncerr
    integer (kind=4), parameter :: nxtrack = 2048, nwl = 1024
    logical :: preflight
    real (kind=8), dimension(nwl, nxtrack) :: sf_asym, sf_hw1e, sf_shape, &
         sf_wavelength
    real (kind=4), dimension(nwl, nxtrack) :: r4_asym, r4_hw1e, r4_shape, &
         r4_wavelength
    real (kind=4), dimension(nwl, nxtrack,1) :: tmp_asym, tmp_hw1e, tmp_shape, tmp_wl
    real (kind=4) :: minwl, maxwl
    logical, dimension(nwl, nxtrack) :: mask

    type (tiof_file_type) :: tio_l1obj

    if (errstat /= 0) return

    ! determine which band we're using
    call lookup_swathname (l1b_channel, swathname, errstat)

    ! determine whether we're reading from irradiance or pre-flight slit func.
    preflight = .true.
    call tiof_open (omi_slitfunc_fname, tio_l1obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "tempo_slitfunc_read: error opening slit function file", errstat)
      return
    endif
    prod_str=""
    ncerr = nf90_inquire_attribute (tio_l1obj%fileid, nf90_global, &
         "product_type")
    if (ncerr == nf90_enotatt) then
      preflight = .true.
    else if (ncerr == nf90_noerr) then
      ncerr = nf90_get_att (tio_l1obj%fileid, nf90_global, "product_type", &
           prod_str)
      if (ncerr .eq. nf90_noerr .and. prod_str(1:3) .eq. "IRR") then
        preflight = .false.
      else
        call tell_error (tell_io_open_error, &
             "tempo_slitfunc_read: undetermined slit function file type", &
             errstat)
        return
      endif
    endif

    if (preflight) then ! read from pre flight slit function
      call tiof_push_group (tio_l1obj, swathname, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_asym", [0,0], [nxtrack,nwl], &
           r4_asym, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_hw1e", [0,0], [nxtrack,nwl], &
           r4_hw1e, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_shape", [0,0], [nxtrack,nwl], &
           r4_shape, errstat)
      call tiof_get2d_r4 (tio_l1obj, "sf_wavelength", [0,0], [nxtrack,nwl], &
           r4_wavelength, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "tempo_slitfunc_read: failed to read pre-flight file", &
             errstat)
        return
      endif
      sf_asym=real(r4_asym, kind=8)
      sf_hw1e=real(r4_hw1e, kind=8)
      sf_shape=real(r4_shape, kind=8)
      sf_wavelength=real(r4_wavelength, kind=8)
      call tell_log (0, "read pre-flight slit function")
    else ! read from irradiance file
      call tiof_push_group (tio_l1obj, swathname, errstat)
      call tiof_get3d_r4 (tio_l1obj, "sf_asym", [0,0,0], [1,nxtrack,nwl], &
           tmp_asym, errstat)
      call tiof_get3d_r4 (tio_l1obj, "sf_hw1e", [0,0,0], [1,nxtrack,nwl], &
           tmp_hw1e, errstat)
      call tiof_get3d_r4 (tio_l1obj, "sf_shape", [0,0,0], [1,nxtrack,nwl], &
           tmp_shape, errstat)
      call tiof_get3d_r4 (tio_l1obj, "wavelength", [0,0,0], [1,nxtrack,nwl], &
           tmp_wl, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, &
             "tempo_slitfunc_read: failed to read irradiance file", &
             errstat)
        return
      endif
      call tell_log (0, "read slit function from irradiance file")
      sf_asym = real(tmp_asym(:,:,1), kind=8)
      sf_hw1e = real(tmp_hw1e(:,:,1), kind=8)
      sf_shape = real(tmp_shape(:,:,1), kind=8)
      sf_wavelength = real(tmp_wl(:,:,1), kind=8)
    endif

    ! determine the mean variable values in the wavelength window
    mask = .true.
    minwl = real(ctrl_fit_winwav_lim(2), kind=4)
    maxwl = real(ctrl_fit_winwav_lim(N_FIT_WINWAV-1), kind=4)
    where (sf_wavelength .lt. minwl .or. sf_wavelength .gt. maxwl)
      mask = .false.
    end where
    !also mask out any bad pixels in the slit function files
    where (sf_asym .gt. 1000.0 .or. sf_asym .lt. 0.0)
      mask = .false.
    end where
    mean_asym = real(sum(sf_asym,mask=mask)/count(mask), kind=8)
    mean_hw1e = real(sum(sf_hw1e,mask=mask)/count(mask), kind=8)
    mean_shape = real(sum(sf_shape,mask=mask)/count(mask), kind=8)
    mean_wl = real(sum(sf_wavelength,mask=mask)/count(mask), kind=8)

    write(*,*) mean_asym, mean_hw1e, mean_shape, mean_wl

  END SUBROUTINE tempo_slitfunc_read

  !--------------------------------------------------------------------------
  !
  ! @param[in]  xtrack_pix  current cross-track pixel number
  ! @param[in]  nwl         number of wavelengths in current spectrum
  ! @param[in]  wvl         wavelength array of current spectrum
  ! @param[in]  spec        current spectrum
  ! @param[out] spec_conv   convolved spectrum
  ! @param      errstat     error handling integer, non-zero = failure
  !
  ! @author  E. O'Sullivan, July 2020
  !--------------------------------------------------------------------------
  SUBROUTINE tempo_slitfunc_convolve ( xtrack_pix, nwvl, wvl, spec, &
       spec_conv, errstat )

    use slitfunction_super_gaussian, only: super_gaussian_sf
    use tell_module

    implicit none

    !input variables
    integer (kind=4), intent(in) :: xtrack_pix, nwvl
    real (kind=8), dimension(:), intent(in) :: wvl, spec

    !output variables
    real (kind=8), dimension(:), intent(out) :: spec_conv
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    call super_gaussian_sf (nwvl, mean_hw1e, mean_asym, mean_shape, 0.0_r8, &
         wvl, spec, spec_conv)

  END SUBROUTINE tempo_slitfunc_convolve

  !--------------------------------------------------------------------------
  ! @param[in]     l1b_irr_filename     irradiance measurement file name
  ! @param[in]     molname              name of retrived molecule (i.e., 'NO2')
  ! @param[in]     solcal_filename      name of file to write 
  ! @param[in]     save_solcal_wvl      wavelengths saved from solar calibration
  ! @param[in]     save_solcal_spec     measured spectrum from solar calibration
  ! @param[in]     save_solcal_resid    fitting residuals from solar calibration
  ! @param[inout]  errstat              error handling integer, non-zero = failure
  !
  ! @author  C. Nowlan, November 2024
  !--------------------------------------------------------------------------
  SUBROUTINE solarcal_write_file(l1b_irr_filename, molname, solcal_filename, save_solcal_wvl, &
        save_solcal_spec, save_solcal_resid, errstat)

    USE netcdf
    USE tell_module
    USE tio_module
    USE tg_names_module
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: nwavel_max, nxtrack_max
    USE OMSAO_variables_module, ONLY: n_fitvar_solcal, mask_fitvar_solcal
    USE OMSAO_omidata_module, ONLY: omi_solcal_pars, omi_solcal_xflag, omi_solcal_chisq, &
                                    omi_solcal_rms, omi_solcal_itnum, omi_irradiance_wght
    USE OMSAO_indices_module, ONLY: hwe_idx, asy_idx, sgk_idx, &
                                    shi_idx, squ_idx, calfit_strings

    IMPLICIT NONE

    ! Intput variables
    CHARACTER (LEN=*), INTENT(IN) :: l1b_irr_filename
    CHARACTER (LEN=*), INTENT(IN) :: molname, solcal_filename
    REAL (KIND=r8), DIMENSION(:,:), INTENT(IN) :: save_solcal_wvl, save_solcal_spec, &
                                                  save_solcal_resid
    INTEGER, INTENT(INOUT) :: errstat

    ! Local variables
    TYPE (tiof_file_type) :: output_file
    TYPE (tiof_file_type) :: l1b
    TYPE (tiof_dimlist_type) :: dimlist
    TYPE (tiof_varlist_type) :: varlist
    TYPE (tiof_attlist_type) :: attlist
    TYPE (tiof_attlist_type) :: att_convergence_flag
    TYPE (tiof_attlist_type) :: att_fit_params

    INTEGER, DIMENSION(1) :: dimid_xtrack
    INTEGER, DIMENSION(2) :: dimids_refwavl_xtrack, dimsizes_refwavl_xtrack
    INTEGER, DIMENSION(2) :: dimids_var_xtrack, dimsizes_var_xtrack
    INTEGER, DIMENSION(2) :: chunksizes
    INTEGER, PARAMETER :: deflate_level = 1
    LOGICAL, PARAMETER :: shuffle = .true.
    INTEGER (KIND=i4) :: i, sidx, eidx, nxtrack, nwvl, nfitvar
    INTEGER, DIMENSION(nxtrack_max) :: xtrack_indices
    REAL (KIND=8), DIMENSION(nxtrack_max) :: solcal_shi, solcal_squ, &
                                             solcal_hwe, solcal_asy, solcal_sgk
    REAL (KIND=8), DIMENSION(n_fitvar_solcal, nxtrack_max) :: solcal_fit_params
    CHARACTER (LEN=3), DIMENSION(n_fitvar_solcal) :: solcal_fit_str
    CHARACTER (LEN=n_fitvar_solcal*4) :: solcal_fit_str_all

    ! using fill values from the original code simplifies diffing output files
    REAL (KIND=8), PARAMETER :: &
      fill_short = -9999, &
      fill_double = -1.0e30_r8

    IF (errstat /= 0) RETURN

    ! Create the file
    CALL tiof_create (output_file, solcal_filename, nf90_clobber, errstat)
    IF (errstat /= 0) THEN
      CALL tell_error (tell_io_write_error, &
                       "solarcal_write_file: creating file " // TRIM(solcal_filename), &
                       errstat)
      RETURN
    ENDIF

    nxtrack = nxtrack_max
    nwvl = nwavel_max
    nfitvar = n_fitvar_solcal

    ! Extract fitted variables from indices
    solcal_fit_params = omi_solcal_pars(mask_fitvar_solcal(1:nfitvar),:)
    solcal_fit_str = calfit_strings(mask_fitvar_solcal(1:nfitvar))

    solcal_fit_str_all = ''
    DO i = 1, nfitvar
      sidx = (i - 1) * 4 + 1
      eidx = sidx + 3
      solcal_fit_str_all(sidx:eidx) = calfit_strings(mask_fitvar_solcal(i)) // " "
    ENDDO

    ! Define a dimension list.
    CALL tiof_dimlist_append (dimlist, tg_dim_xtrack, nxtrack, errstat)
    CALL tiof_dimlist_append (dimlist, tg_dim_refwavl, nwvl, errstat)
    CALL tiof_dimlist_append (dimlist, tg_dim_fitvar, nfitvar, errstat)
    CALL tiof_def_dims (output_file, dimlist, errstat)

    CALL tiof_dimlist_lookup (dimlist, [tg_dim_xtrack], dimid_xtrack, errstat)
    CALL tiof_dimlist_lookup (dimlist, &
                              [tg_dim_refwavl, tg_dim_xtrack], &
                              dimids_refwavl_xtrack, &
                              errstat, dimsizes = dimsizes_refwavl_xtrack)
    CALL tiof_dimlist_lookup (dimlist, &
                              [tg_dim_fitvar, tg_dim_xtrack], &
                              dimids_var_xtrack, &
                              errstat, dimsizes = dimsizes_var_xtrack)

    ! Global attributes
    CALL tiof_history_append_cmdline (output_file)
    CALL tiof_put_git_commit_hash (output_file, errstat)
    CALL tiof_attlist_append (attlist, errstat, "product_type", att_text = molname)
    CALL tiof_def_atts (output_file, attlist, nf90_global, errstat)
    CALL tiof_attlist_free (attlist)

    ! Copy irradiance timestamps to global attributes
    CALL tiof_open (l1b_irr_filename, l1b, nf90_nowrite, errstat)
    CALL tiof_copy_granule_ident (l1b, output_file, errstat)
    CALL tiof_close (l1b, errstat)
    CALL tiof_write_epoch_timestamp (output_file, errstat)
    IF (errstat /= 0) THEN
      CALL tell_error (tell_runtime_error, &
                       "copying metadata from " // TRIM(l1b_irr_filename), &
                       errstat)
      RETURN
    ENDIF
  
    ! Convergence flag attributes
    CALL tiof_attlist_append (att_convergence_flag, errstat, "flag_meanings", &
                              att_text = "failed maxiter_exceeded suspect good")
    CALL tiof_attlist_append (att_convergence_flag, errstat, "flag_values", &
                              att_i4 = [-2,-1,0,1])

    ! Fittied paramater attribute
    CALL tiof_attlist_append (att_fit_params, errstat, "fit_parameter_names", &
                              att_text = solcal_fit_str_all)

    ! Write variables
    chunksizes(1) = dimsizes_refwavl_xtrack(1)              ! wavelength dimension
    chunksizes(2) = MIN(dimsizes_refwavl_xtrack(2), 128)    ! xtrack dimension

    CALL tiof_varlist_append (varlist, errstat, tg_dim_xtrack, nf90_int, &
                             dimids = dimid_xtrack, &
                             long_name = "pixel index along slit")
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_convergence_flag, &
                              nf90_short, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar calibration convergence flag", &
                              valid_range = [-2.0_r8, 1.0_r8], &
                              fillvalue = fill_short, &
                              attlist = att_convergence_flag)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_shift, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              units = 'nm', &
                              long_name = "solar calibration wavelength shift", &
                              valid_range = [-10.0_r8, 10.0_r8], &
                              fillvalue = fill_double)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_squeeze, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar calibration wavelength squeeze", &
                              valid_range = [-10.0_r8, 10.0_r8], &
                              fillvalue = fill_double)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_chisq, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar calibration fit chi-squared", &
                              valid_range = [0.0_r8, 1000.0_r8], &
                              fillvalue = fill_double)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_itnum, &
                              nf90_short, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar calibration number of iterations", &
                              valid_range = [0.0_r8, 100.0_r8], &
                              fillvalue = fill_short)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_rms, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar calibration fitting RMS", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_double)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_sf_hw1e, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              units = 'nm', &
                              long_name = "solar calibration slit function Super-Gaussian half-width at 1/e", &
                              valid_range = [-10.0_r8, 10.0_r8], &
                              fillvalue = fill_double)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_sf_asym, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar calibration slit function Super-Gaussian asymmetry parameter", &
                              valid_range = [-10.0_r8, 10.0_r8], &
                              fillvalue = fill_double)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_sf_shape, &
                              nf90_double, &
                              dimids = dimid_xtrack,  &
                              long_name = "solar calibration slit function Super-Gaussian shape parameter", &
                              valid_range = [0.0_r8, 100.0_r8], &
                              fillvalue = fill_double)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_wavelengths, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              units = "nm", &
                              long_name = "calibrated solar spectrum wavelengths", &
                              valid_range = [100.0_r8, 1000.0_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_spectrum, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              long_name = "normalized solar spectrum", &
                              valid_range = [-100.0_r8, 100.0_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_weight, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              long_name = "solar spectrum weighting", &
                              valid_range = [0.0_r8, 1.0_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)
    CALL tiof_varlist_append (varlist, errstat, &
                              tg_var_solcal_residuals, &
                              nf90_double, &
                              dimids = dimids_refwavl_xtrack,  &
                              long_name = "solar spectrum residuals", &
                              comment = "fit residuals from solar spectrum wavelength and slit function calibration", &
                              valid_range = [-1e30_r8, 1e30_r8], &
                              fillvalue = fill_double, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              chunksizes = chunksizes)

   chunksizes(1) = dimsizes_var_xtrack(1)           ! var dimension
   chunksizes(2) = MIN(dimsizes_var_xtrack(2),128)  ! xtrack dimension

   CALL tiof_varlist_append (varlist, errstat, &
                             tg_var_solcal_fit_params, &
                             nf90_double, &
                             dimids = dimids_var_xtrack, &
                             long_name = "solar irradiance fit parameter", &
                             valid_range=[-1e30_r8, 1e30_r8], &
                             fillvalue = fill_double, &
                             attlist = att_fit_params, &
                             deflate_level = deflate_level, &
                             shuffle = shuffle, &
                             chunksizes = chunksizes)
    
    CALL tiof_def_vars (output_file, varlist, errstat)
    CALL tiof_varlist_free (varlist)
    CALL tiof_dimlist_free (dimlist)
    CALL tiof_attlist_free (att_convergence_flag)
    CALL tiof_attlist_free (att_fit_params)

    ! Extract information from 2D arrays
    solcal_shi = omi_solcal_pars(shi_idx,:)
    solcal_squ = omi_solcal_pars(squ_idx,:)
    solcal_hwe = omi_solcal_pars(hwe_idx,:)
    solcal_asy = omi_solcal_pars(asy_idx,:)    
    solcal_sgk = omi_solcal_pars(sgk_idx,:)
 
    ! Write data 
    xtrack_indices = [(i, i=0,nxtrack-1)]
    CALL tiof_put1d_i4 (output_file, tg_dim_xtrack, [0], [nxtrack], xtrack_indices, errstat) 
    CALL tiof_put1d_i2 (output_file, tg_var_solcal_convergence_flag, [0], [nxtrack], &
                        omi_solcal_xflag(1:nxtrack), errstat)
    CALL tiof_put1d_r8 (output_file, tg_var_solcal_shift, [0], [nxtrack], &
                        solcal_shi(1:nxtrack), errstat)
    CALL tiof_put1d_r8 (output_file, tg_var_solcal_squeeze, [0], [nxtrack], &
                        solcal_squ(1:nxtrack), errstat)
    CALL tiof_put1d_r8 (output_file, tg_var_solcal_chisq, [0], [nxtrack], &
                        omi_solcal_chisq(1:nxtrack), errstat)
    CALL tiof_put1d_i2 (output_file, tg_var_solcal_itnum, [0], [nxtrack], &
                        omi_solcal_itnum(1:nxtrack), errstat)
    CALL tiof_put1d_r8 (output_file, tg_var_solcal_rms, [0], [nxtrack], &
                        omi_solcal_rms(1:nxtrack), errstat)
    CALL tiof_put1d_r8 (output_file, tg_var_solcal_sf_hw1e, [0], [nxtrack], &
                        solcal_hwe(1:nxtrack), errstat)
    CALL tiof_put1d_r8 (output_file, tg_var_solcal_sf_asym, [0], [nxtrack], &
                        solcal_asy(1:nxtrack), errstat)
    CALL tiof_put1d_r8 (output_file, tg_var_solcal_sf_shape, [0], [nxtrack], &
                        solcal_sgk(1:nxtrack), errstat)
    CALL tiof_put2d_r8 (output_file, tg_var_solcal_wavelengths, [0,0], [nxtrack, nwvl], &
                        save_solcal_wvl(1:nwvl,1:nxtrack), errstat)
    CALL tiof_put2d_r8 (output_file, tg_var_solcal_spectrum, [0,0], [nxtrack, nwvl], &
                        save_solcal_spec(1:nwvl,1:nxtrack), errstat)
    CALL tiof_put2d_r8 (output_file, tg_var_solcal_weight, [0,0], [nxtrack, nwvl], &
                        omi_irradiance_wght(1:nwvl,1:nxtrack), errstat)
    CALL tiof_put2d_r8 (output_file, tg_var_solcal_residuals, [0,0], [nxtrack, nwvl], &
                        save_solcal_resid (1:nwvl,1:nxtrack), errstat)
    CALL tiof_put2d_r8 (output_file, tg_var_solcal_fit_params, [0,0], [nxtrack,nfitvar], &
                        solcal_fit_params(1:nfitvar,1:nxtrack), &
                        errstat)

    CALL tell_log (1, "Writing solar irradiance calibration to file: " // TRIM(solcal_filename))

    ! Close file
    CALL tiof_close (output_file, errstat)
    IF (errstat /= 0) THEN
      CALL tell_error (tell_io_error, "solarcal_write_file", errstat)
    ENDIF

  END SUBROUTINE solarcal_write_file

  !--------------------------------------------------------------------------
  !
  ! @param[in]     molname              name of retrived molecule (i.e., 'NO2')
  ! @param[in]     solcal_filename      name of file to write 
  ! @param[inout]  save_solcal_wvl      wavelengths saved from solar calibration
  ! @param[inout]  save_solcal_spec     measured spectrum from solar calibration
  ! @param[inout]  save_solcal_resid    fitting residuals from solar calibration
  ! @param         errstat              error handling integer, non-zero = failure
  !
  ! @author  C. Nowlan, November 2024
  !--------------------------------------------------------------------------
  SUBROUTINE solarcal_read_file(molname, solcal_filename, save_solcal_wvl, &
        save_solcal_spec, save_solcal_resid, errstat)

    USE netcdf
    USE tell_module
    USE tio_module
    USE tg_names_module
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: nwavel_max, nxtrack_max
    USE OMSAO_omidata_module, ONLY: omi_solcal_pars, omi_solcal_xflag, omi_solcal_chisq, &
                                    omi_solcal_rms, omi_solcal_shift, omi_irradiance_wght
    USE OMSAO_variables_module, ONLY: n_fitvar_solcal, mask_fitvar_solcal
    USE irradiance_data, ONLY: Irr_Data
    use ctrlvars, only: yn_I0, yn_spectrum_norm

    IMPLICIT NONE

    ! Input variables
    CHARACTER (LEN=*), INTENT(IN) :: molname, solcal_filename
    REAL (KIND=r8), DIMENSION(:,:), INTENT(INOUT) :: save_solcal_wvl, save_solcal_spec, &
                                                     save_solcal_resid
    INTEGER, INTENT(INOUT) :: errstat

    ! Local variables
    INTEGER(KIND=i4) :: i, nxtrack, nwvl, nfitvar
    REAL (KIND=8), DIMENSION(n_fitvar_solcal, nxtrack_max) :: solcal_fit_params
    CHARACTER (LEN=16) :: prod_str
    INTEGER (KIND=4) :: ncerr

    TYPE (tiof_file_type) :: input_file


    nxtrack = nxtrack_max
    nwvl = nwavel_max
    nfitvar = n_fitvar_solcal
  
    CALL tiof_open (solcal_filename, input_file, nf90_nowrite, errstat)
    IF (errstat /= 0) THEN
      CALL tell_error (tell_io_open_error, &
           "solarcal_read_file: error opening solar calibration file", errstat)
      RETURN
    ENDIF

    prod_str=""
    ncerr = nf90_inquire_attribute (input_file%fileid, nf90_global, &
         "product_type")
    IF (ncerr == nf90_enotatt) THEN
       CALL tell_error (tell_io_open_error, &
             "solarcal_read_file: undetermined product type", &
             errstat)
       RETURN
    ELSEIF (ncerr == nf90_noerr) THEN
       ncerr = nf90_get_att (input_file%fileid, nf90_global, "product_type", &
                             prod_str)
       IF (TRIM(prod_str) .NE. molname) THEN
          CALL tell_error (tell_io_open_error, &
                          "solarcal_read_file: wrong product type", &
                           errstat)
          RETURN
       ENDIF
    ENDIF

    CALL tell_log (1, "Reading solar irradiance calibration from file: " // TRIM(solcal_filename))

    ! Read data from file
    CALL tiof_get1d_i2 (input_file, tg_var_solcal_convergence_flag, [0], [nxtrack], &
                        omi_solcal_xflag(1:nxtrack), errstat)
    CALL tiof_get1d_r8 (input_file, tg_var_solcal_shift, [0], [nxtrack], &
                        omi_solcal_shift(1:nxtrack), errstat)
    CALL tiof_get1d_r8 (input_file, tg_var_solcal_chisq, [0], [nxtrack], &
                        omi_solcal_chisq(1:nxtrack), errstat)
    CALL tiof_get1d_r8 (input_file, tg_var_solcal_rms, [0], [nxtrack], &
                        omi_solcal_rms(1:nxtrack), errstat)
    CALL tiof_get2d_r8 (input_file, tg_var_solcal_wavelengths, [0,0], [nxtrack, nwvl], &
                        save_solcal_wvl(1:nwvl,1:nxtrack), errstat)
    CALL tiof_get2d_r8 (input_file, tg_var_solcal_spectrum, [0,0], [nxtrack, nwvl], &
                        save_solcal_spec (1:nwvl,1:nxtrack), errstat)
    CALL tiof_get2d_r8 (input_file, tg_var_solcal_weight, [0,0], [nxtrack, nwvl], &
                        omi_irradiance_wght (1:nwvl,1:nxtrack), errstat)
    CALL tiof_get2d_r8 (input_file, tg_var_solcal_residuals, [0,0], [nxtrack, nwvl], &
                        save_solcal_resid (1:nwvl,1:nxtrack), errstat)
    CALL tiof_get2d_r8 (input_file, tg_var_solcal_fit_params, [0,0], [nxtrack,nfitvar], &
                        solcal_fit_params(1:nfitvar,1:nxtrack), &
                        errstat)

    ! Close file
    CALL tiof_close (input_file, errstat)
    IF (errstat /= 0) THEN
      CALL tell_error (tell_io_error, "solarcal_read_file", errstat)
    ENDIF

    ! Store data in other variables that get used later
    if (.not. yn_I0) then
      Irr_Data%wavelengths = save_solcal_wvl
      Irr_Data%spectrum = save_solcal_spec
    else if (yn_spectrum_norm) then
      ! ---------------------------------------
      ! Normalize I0 spectrum if necessary.
      ! Usually is done during solar wavelength
      ! calibration but we are skipping it.
      ! ---------------------------------------
      do i = 1, nxtrack
        Irr_Data%spectrum(1:Irr_Data%nwaves(i),i) = &
            Irr_Data%spectrum(1:Irr_Data%nwaves(i),i) / &
            (SUM(Irr_Data%spectrum(1:Irr_Data%nwaves(i),i)) / REAL(Irr_Data%nwaves(i),kind=8))
      end do
    endif

    DO i = 1, n_fitvar_solcal
       omi_solcal_pars(mask_fitvar_solcal(i), :) = solcal_fit_params(i,:)
    END DO 
 
  END SUBROUTINE solarcal_read_file

END MODULE slitfunction_tempo
