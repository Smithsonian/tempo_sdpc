!>Subroutines to read L2 Ozone Profile netCDF files
module m_read_L2_o3p_tio
  use netcdf
  use tio_module
  use tell_module
  use o3p_names_module


  implicit none

  private
  public read_o3p_dims, read_o3p_dim_indices, read_o3p_product, &
       read_o3p_geolocation, read_o3p_support, read_o3p_diagnostic, &
       read_o3p_qastat, open_o3p, close_o3p



contains

  !> Open L2 netCDF O3 profile product and get dimensions
  !--------------------------------------------------------------------------
  !
  !> @param[in]  l2file       L2 O3 profile product filename
  !> @param      tio_l2obj    L2 file object
  !> @param[out] nstep        mirror step dimension size
  !> @param[out] nxtrack      cross-track dimension size
  !> @param[out] ncorner      corner dimension size
  !> @param[out] nfitvars     fit_variables dimension size
  !> @param[out] nfitwins     fitting_windows dimension size
  !> @param[out] ngas         gases dimension size (zero if unused)
  !> @param[out] nlayer       layer dimension size
  !> @param[out] nlayerp1     layer_plus1 dimension size
  !> @param[out] nmax_wavs    max_wavelengths dimension size (zero if unused)
  !> @param[out] nnoise_elems noise_elements dimension size (zero if unused)
  !> @param[out] nnongas      non_gas_variables dimension size
  !> @param[out] naeros_wavs  aerosol_wavelengths dim. size (zero if unused)
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine read_o3p_dims(l2file, tio_l2obj, nstep, nxtrack, ncorner, &
       nfitvars, nfitwins, ngas, nlayer, nlayerp1, nmax_wavs, nnoise_elems, &
       nnongas, naeros_wavs, errstat)

    implicit none

    !input variables
    character (len=*), intent (in) :: l2file

    !output variables
    integer (kind=4), intent (out) :: nstep, nxtrack, ncorner, &
         nfitvars, nfitwins, ngas, nlayer, nlayerp1, nmax_wavs, nnoise_elems, &
         nnongas, naeros_wavs
    integer (kind=4), intent (inout) :: errstat

    !local variables
    integer (kind=4) :: status, dimid

    type (tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_open (l2file, tio_l2obj, nf90_nowrite, errstat)
    ! dimension always present
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_step, nstep, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_xtrack, nxtrack, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_corner, ncorner, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_fitvar, nfitvars, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_windows, nfitwins, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_layer, nlayer, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_layerp1, nlayerp1, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_param, nnongas, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_o3p_dims: failed to read standard dimensions", errstat)
      call tiof_close(tio_l2obj, errstat)
      return
    endif

    ! optional dimensions - check if present in the file
    status = nf90_inq_dimid (tio_l2obj%groupid, o3p_dim_gas, dimid)
    if (status == nf90_noerr) then
      call tiof_inq_dimlen (tio_l2obj, o3p_dim_gas, ngas, errstat)
    else
      ngas = 0
    endif
    status = nf90_inq_dimid (tio_l2obj%groupid, o3p_dim_wavl_max, dimid)
    if (status == nf90_noerr) then
      call tiof_inq_dimlen (tio_l2obj, o3p_dim_wavl_max, nmax_wavs, errstat)
    else
      nmax_wavs = 0
    endif
    status = nf90_inq_dimid (tio_l2obj%groupid, o3p_dim_elms, dimid)
    if (status == nf90_noerr) then
      call tiof_inq_dimlen (tio_l2obj, o3p_dim_elms, nnoise_elems, errstat)
    else
      nnoise_elems = 0
    endif
    status = nf90_inq_dimid (tio_l2obj%groupid, o3p_dim_aeros_wavl, dimid)
    if (status == nf90_noerr) then
      call tiof_inq_dimlen (tio_l2obj, o3p_dim_aeros_wavl, naeros_wavs, &
           errstat)
    else
      naeros_wavs = 0
    endif

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_o3p_dims: failed to read optional dimensions", errstat)
    endif

    call tiof_close(tio_l2obj, errstat)

  end subroutine read_o3p_dims



  !> Open L2 netCDF O3 profile product and read step and xtrack index arrays
  !--------------------------------------------------------------------------
  !
  !> @param[in]  l2file       L2 O3 profile product filename
  !> @param      tio_l2obj    L2 file object
  !> @param[in]  nstep        mirror step dimension size
  !> @param[in]  nxtrack      cross-track dimension size
  !> @param[out] step         array of step indices
  !> @param[out] xtrack       array of xtrack indices
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine read_o3p_dim_indices (l2file, tio_l2obj, nstep, nxtrack, &
       step, xtrack, errstat)

    implicit none

    !input variables
    character (len=*), intent(in) :: l2file
    integer (kind=4), intent(in) :: nstep, nxtrack

    !output variables
    integer (kind=4), dimension(nstep), intent(out) :: step
    integer (kind=4), dimension(nxtrack), intent(out) :: xtrack
    integer (kind=4), intent(inout) :: errstat

    type(tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_open (l2file, tio_l2obj, nf90_nowrite, errstat)
    call tiof_get1d_i4 (tio_l2obj, o3p_dim_step, [0], [nstep], step, errstat)
    call tiof_get1d_i4 (tio_l2obj, o3p_dim_xtrack, [0], [nxtrack], &
         xtrack, errstat)
    call tiof_close(tio_l2obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_o3p_dim_indices: failed", errstat)
    endif

  end subroutine read_o3p_dim_indices



  !> Read geolocation parameters from an L2 netCDF O3 profile product
  !--------------------------------------------------------------------------
  !
  !> @param      tio_l2obj    L2 file object
  !> @param[in]  nstep        mirror step dimension size
  !> @param[in]  nxtrack      cross-track dimension size
  !> @param[in]  ncorner      corner dimension size
  !> @param[in]  min_xtrack   min value of xtrack array to be read into
  !> @param[in]  max_xtrack   max value of xtrack array to be read into
  !> @param[in]  min_step     min value of step array to be read into
  !> @param[in]  max_step     max value of step array to be read into
 !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine read_o3p_geolocation (tio_l2obj, nstep, nxtrack, ncorner, &
       min_xtrack, max_xtrack, min_step, max_step, errstat)

    use m_o3p_params, only: time, geoflg, lat, lon, aza, sza, vza, &
         corner_lat, corner_lon

    implicit none

    !input variables
    integer (kind=4), intent(in) :: nstep, nxtrack, ncorner, &
         min_xtrack, max_xtrack, min_step, max_step

    !output variables
    integer (kind=4), intent(inout) :: errstat

    type(tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_push_group (tio_l2obj, o3p_grp_geolocation, errstat)
    call tiof_get1d_r8 (tio_l2obj, o3p_var_time, [0], [nstep], &
         time(min_step:max_step), errstat)
    call tiof_get2d_ui2 (tio_l2obj, o3p_var_geoflg, [0, 0], &
         [nstep, nxtrack], geoflg(min_xtrack:max_xtrack, min_step:max_step), &
         errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_latitude, [0,0], &
         [nstep, nxtrack], lat(min_xtrack:max_xtrack, min_step:max_step), &
         errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_longitude, [0,0], &
         [nstep, nxtrack], lon(min_xtrack:max_xtrack, min_step:max_step), &
         errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_ra_angle, [0,0], &
         [nstep, nxtrack], aza(min_xtrack:max_xtrack, min_step:max_step), &
         errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_sz_angle, [0,0], &
         [nstep, nxtrack], sza(min_xtrack:max_xtrack, min_step:max_step), &
         errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_vz_angle, [0,0], &
         [nstep, nxtrack], vza(min_xtrack:max_xtrack, min_step:max_step), &
         errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_latitude_bounds, [0,0,0], &
         [nstep, nxtrack, ncorner], &
         corner_lat(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_longitude_bounds, [0,0,0], &
         [nstep, nxtrack, ncorner], &
         corner_lon(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_pop_group (tio_l2obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_o3p_geolocation failed", &
           errstat)
    endif

  end subroutine read_o3p_geolocation


  !> Read product parameters from an L2 netCDF O3 profile product
  !--------------------------------------------------------------------------
  !
  !> @param      tio_l2obj    L2 file object
  !> @param[in]  nstep        mirror step dimension size
  !> @param[in]  nxtrack      cross-track dimension size
  !> @param[in]  nlayer       layer dimension size
  !> @param[in]  ngas         gas dimension size
  !> @param[in]  nnongas      nongas dimension size
  !> @param[in]  min_xtrack   min value of xtrack array to be read into
  !> @param[in]  max_xtrack   max value of xtrack array to be read into
  !> @param[in]  min_step     min value of step array to be read into
  !> @param[in]  max_step     max value of step array to be read into
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine read_o3p_product (tio_l2obj, nstep, nxtrack, nlayer, ngas, &
       nnongas, min_xtrack, max_xtrack, min_step, max_step, errstat)

    use ozprof_data_module, only: gaswrt, ozwrtvar
    use m_o3p_params, only: o3tot, o3tot_prec, o3tot_err, o3strat, &
         o3strat_prec, o3strat_err, o3trop, o3trop_prec, o3trop_err, ozprof, &
         ozprof_prec, ozprof_err, gas, gas_prec, gas_err, nongas, &
         nongas_prec, nongas_err

    implicit none

    !input variables
    integer (kind=4), intent(in) :: nstep, nxtrack, ngas, nnongas, &
         nlayer, min_xtrack, max_xtrack, min_step, max_step

    !output variables
    integer (kind=4), intent(inout) :: errstat

    type(tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_push_group (tio_l2obj, o3p_grp_product, errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_o3_retrieve_prof, [0, 0, 0], &
         [nstep,nxtrack,nlayer], &
         ozprof(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_o3_retrieve_prof_prec, [0, 0, 0], &
         [nstep,nxtrack,nlayer], &
         ozprof_prec(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_o3_retrieve_prof_err, [0, 0, 0], &
         [nstep,nxtrack,nlayer], &
         ozprof_err(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_total_o3, [0, 0], &
         [nstep, nxtrack], &
         o3tot(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_total_o3_prec, [0, 0], &
         [nstep, nxtrack], &
         o3tot_prec(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_total_o3_err, [0, 0], &
         [nstep, nxtrack], &
         o3tot_err(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_strat_o3, [0, 0], &
         [nstep, nxtrack], &
         o3strat(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_strat_o3_prec, [0, 0], &
         [nstep, nxtrack], &
         o3strat_prec(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_strat_o3_err, [0, 0], &
         [nstep, nxtrack], &
         o3strat_err(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_tropo_o3, [0, 0], &
         [nstep, nxtrack], &
         o3trop(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_tropo_o3_prec, [0, 0], &
         [nstep, nxtrack], &
         o3trop_prec(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_tropo_o3_err, [0, 0], &
         [nstep, nxtrack], &
         o3trop_err(min_xtrack:max_xtrack, min_step:max_step), errstat)
    ! FIXME - add surface layer (0-2km) O3 column when UVIS code implemented
    ! Optional products - gas
    if (gaswrt .and. ngas > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_other_gas_retrieve, [0, 0, 0], &
           [nstep, nxtrack, ngas], &
           gas(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_other_gas_retrieve_prec, &
           [0, 0, 0], [nstep, nxtrack, ngas], &
           gas_prec(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_other_gas_retrieve_err, &
           [0, 0, 0], [nstep, nxtrack, ngas], &
           gas_err(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif
    ! Optional products - nongas
    if(ozwrtvar .and. nnongas > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_nongas_param_retrieve, &
           [0, 0, 0], [nstep, nxtrack, nnongas], &
           nongas(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_nongas_param_ret_prec, &
           [0, 0, 0], [nstep, nxtrack, nnongas], &
           nongas_prec(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_nongas_param_ret_err, &
           [0, 0, 0], [nstep, nxtrack, nnongas], &
           nongas_err(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif

    call tiof_pop_group (tio_l2obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_o3p_product failed", &
           errstat)
    endif

  end subroutine read_o3p_product


  !> Read support parameters from an L2 netCDF O3 profile product
  !--------------------------------------------------------------------------
  !
  !> @param      tio_l2obj     L2 file object
  !> @param[in]  nstep         mirror step dimension size
  !> @param[in]  nxtrack       cross-track dimension size
  !> @param[in]  nlayer        layer dimension size
  !> @param[in]  nlayerp1      layer+1 dimension size
  !> @param[in]  nwindow       fit windows dimension size
  !> @param[in]  ngas          gas dimension size
  !> @param[in]  nnongas       nongas dimension size
  !> @param[in]  nfitvars      fit variables dimension size
  !> @param[in]  nnoise_elms   noise elements dimension size
  !> @param[in]  nmax_wavs     maximum wavelength dimension size
  !> @param[in]  naeros_wavs   aerosol wavelength dimension size
  !> @param[in]  min_xtrack   min value of xtrack array to be read into
  !> @param[in]  max_xtrack   max value of xtrack array to be read into
  !> @param[in]  min_step     min value of step array to be read into
  !> @param[in]  max_step     max value of step array to be read into
  !> @param      errstat       error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine read_o3p_support (tio_l2obj, nstep, nxtrack, nlayer, nlayerp1, &
       nwindow, ngas, nnongas, nfitvars, nnoise_elms, nmax_wavs, naeros_wavs, &
       min_xtrack, max_xtrack, min_step, max_step, errstat)

    use ozprof_data_module, only: ozwrtavgk, ozwrtcorr, ozwrtcovar, &
         ozwrtcontri, ozwrtvar, gaswrt, aerosol, do_lambcld
    use m_o3p_params, only: aeros_idx, tropo_idx, cld_frac, cld_pres, &
         cld_flag, glintprob, eff_alb, o3apriori, o3apriori_err, &
         ozprof_pres, ozprof_alt, ozprof_temp, ozinfo, cld_opt_depth, &
         gas_apriori, gas_apriori_err, nongas_apriori, nongas_apriori_err, &
         correl_mtrx, contrib_mtrx, aeros_opt_thick, aeros_scatter_thick, &
         avg_kernel, noise_mtrx, n_fit_wvl, n_window_wvl, gas_names, &
         nongas_names, nongas_units

    implicit none

    !input variables
    integer (kind=4), intent(in) :: nstep, nxtrack, ngas, nnongas, &
         nlayer, nlayerp1, nwindow, nfitvars, nnoise_elms, nmax_wavs, &
         naeros_wavs, min_xtrack, max_xtrack, min_step, max_step

    !output variables
    integer (kind=4), intent(inout) :: errstat

    type(tiof_file_type) :: tio_l2obj

    integer :: i, j

    if (errstat /= 0) return


    call tiof_push_group (tio_l2obj, o3p_grp_support_data, errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_aeros_index, [0,0], &
         [nstep, nxtrack], &
         aeros_idx(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_profile_pres, [0, 0, 0], &
         [nstep, nxtrack, nlayerp1], &
         ozprof_pres(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_profile_alt, [0, 0, 0], &
         [nstep, nxtrack, nlayerp1], &
         ozprof_alt(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_profile_temp, [0, 0, 0], &
         [nstep, nxtrack, nlayerp1], &
         ozprof_temp(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_o3_apriori_prof, [0, 0, 0], &
         [nstep, nxtrack, nlayer], &
         o3apriori(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_o3_apriori_prof_err, [0, 0, 0], &
         [nstep, nxtrack, nlayer], &
         o3apriori_err(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_i4 (tio_l2obj, o3p_var_tropo_index, [0, 0], &
         [nstep, nxtrack], &
         tropo_idx(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_eff_cld_frac, [0, 0], &
         [nstep, nxtrack], &
         cld_frac(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_eff_cld_pres, [0,0], &
         [nstep, nxtrack], &
         cld_pres(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_i4 (tio_l2obj, o3p_var_cld_flag, [0, 0], &
         [nstep, nxtrack], &
         cld_flag(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_glint_prob, [0, 0], &
         [nstep, nxtrack], &
         glintprob(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_r4 (tio_l2obj, o3p_var_surf_albedo, [0, 0], &
         [nstep, nxtrack], &
         eff_alb(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get2d_ui4 (tio_l2obj, o3p_var_fit_wavel, [0, 0], &
         [nstep, nxtrack], &
         n_fit_wvl(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_ui4 (tio_l2obj, o3p_var_window_wavel, [0, 0, 0], &
         [nstep, nxtrack, nwindow], &
         n_window_wvl(:, min_xtrack:max_xtrack, min_step:max_step), errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_o3p_support: failed to read standar variables", errstat)
      return
    endif

    ! Optional products - gas
    if (gaswrt .and. ngas > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_other_gas_apriori, &
           [0, 0, 0], [nstep, nxtrack, ngas], &
           gas_apriori(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_other_gas_apriori_err, &
           [0, 0, 0], [nstep, nxtrack, ngas], &
           gas_apriori_err(:, min_xtrack:max_xtrack, min_step:max_step), &
           errstat)
      call tiof_get1d_string (tio_l2obj, o3p_var_other_gas_names, 0, &
           ngas, gas_names, errstat)
    endif
    ! Optional products - nongas
    if (ozwrtvar .and. nnongas > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_nongas_param_apriori, &
           [0, 0, 0], [nstep, nxtrack, nnongas], &
           nongas_apriori(:, min_xtrack:max_xtrack, min_step:max_step), &
           errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_nongas_param_apriori_err, &
           [0, 0, 0], [nstep, nxtrack, nnongas], &
           nongas_apriori_err(:, min_xtrack:max_xtrack, min_step:max_step), &
           errstat)
      call tiof_get1d_string (tio_l2obj, o3p_var_nongas_param_names, 0, &
           nnongas, nongas_names, errstat)
      call tiof_get1d_string (tio_l2obj, o3p_var_nongas_param_units, 0, &
           nnongas, nongas_units, errstat)
    endif
    ! Optional - averaging kernels
    if (ozwrtavgk .and. nlayer > 0) then
      call tiof_get4d_i2 (tio_l2obj, o3p_var_o3_avg_kernel, [0, 0, 0, 0], &
           [nstep, nxtrack, nlayer, nlayer], &
           avg_kernel(:,:,min_xtrack:max_xtrack,min_step:max_step), errstat)
    endif
    ! Optional - correlation matrix
    if (ozwrtcorr .and. nfitvars > 0) then
      call tiof_get4d_r4 (tio_l2obj, o3p_var_correl, [0, 0, 0, 0], &
           [nstep, nxtrack, nfitvars, nfitvars], &
           correl_mtrx(:, :, min_xtrack:max_xtrack, min_step:max_step), &
           errstat)
    endif
    ! Optional - noise matrix & ozone information content matrix
    if (ozwrtcovar .and. nnoise_elms > 0) then
      call tiof_get3d_i2 (tio_l2obj, o3p_var_o3_noise_matrix, [0, 0, 0], &
           [nstep, nxtrack, nnoise_elms], &
           noise_mtrx(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
      call tiof_get2d_r4 (tio_l2obj, o3p_var_o3_info_content, [0, 0], &
           [nstep, nxtrack], &
           ozinfo(min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif
    ! Optional - contribution matrix
    if (ozwrtcontri .and. nfitvars > 0 .and. nmax_wavs > 0) then
      call tiof_get4d_r4 (tio_l2obj, o3p_var_contrib_func, [0, 0, 0, 0], &
           [nstep, nxtrack, nfitvars, nmax_wavs], &
           contrib_mtrx(:, :, min_xtrack:max_xtrack, min_step:max_step), &
           errstat)
    endif
    ! Optional - cloud optical depth
    if (.not. do_lambcld) then
      call tiof_get2d_r4 (tio_l2obj, o3p_var_cld_opt_depth, [0, 0], &
           [nstep, nxtrack], &
           cld_opt_depth(min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif
    ! Optional - aerosols
    if (aerosol .and. naeros_wavs > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_aeros_opt_thick, [0, 0, 0], &
           [nstep, nxtrack, naeros_wavs], &
           aeros_opt_thick(:, min_xtrack:max_xtrack, min_step:max_step), &
           errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_aeros_scatter_thick, [0, 0, 0], &
           [nstep, nxtrack, naeros_wavs], &
           aeros_scatter_thick(:, min_xtrack:max_xtrack, min_step:max_step), &
           errstat)
    endif

    call tiof_pop_group (tio_l2obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_o3p_support: failed to read optional variables", errstat)
    endif

  end subroutine read_o3p_support



  !> Read diagnostic parameters from an L2 netCDF O3 profile product
  !--------------------------------------------------------------------------
  !
  !> @param      tio_l2obj    L2 file object
  !> @param[in]  nstep        mirror step dimension size
  !> @param[in]  nxtrack      cross-track dimension size
  !> @param[in]  nmax_wavs    maximum wavelengths dimension size
  !> @param[in]  nfitvars     fit variables dimension size
  !> @param[in]  min_xtrack   min value of xtrack array to be read into
  !> @param[in]  max_xtrack   max value of xtrack array to be read into
  !> @param[in]  min_step     min value of step array to be read into
  !> @param[in]  max_step     max value of step array to be read into
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine read_o3p_diagnostic (tio_l2obj, nstep, nxtrack, nmax_wavs, &
       nfitvars, min_xtrack, max_xtrack, min_step, max_step, errstat)

    use ozprof_data_module, only: ozwrtres, ozwrtwf
    use OMSAO_variables_module, only: reduce_resolution
    use m_o3p_params, only : wavelengths, norm_rad, sim_norm_rad, wgt_func

    implicit none

    !input variables
    integer (kind=4), intent(in) :: nstep, nxtrack, nmax_wavs, nfitvars, &
         min_xtrack, max_xtrack, min_step, max_step

    !output variables
    integer (kind=4), intent(inout) :: errstat

    type(tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_push_group (tio_l2obj, o3p_grp_diagnostic, errstat)
    ! All parameters are optional
    if (.not. reduce_resolution .and. nmax_wavs > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_wavel, [0, 0, 0], &
           [nstep, nxtrack, nmax_wavs], &
           wavelengths(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif
    if (ozwrtwf .and. nmax_wavs > 0 .and. nfitvars > 0) then
      call tiof_get4d_r4 (tio_l2obj, o3p_var_weight_func, [0, 0, 0, 0], &
           [nstep, nxtrack, nmax_wavs, nfitvars], &
           wgt_func(:, :, min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif
    if (ozwrtres .and. nmax_wavs > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_sim_norm_rad, [0, 0, 0], &
           [nstep, nxtrack, nmax_wavs], &
           sim_norm_rad(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
      call tiof_get3d_r4 (tio_l2obj, o3p_var_norm_radiance, [0, 0, 0], &
           [nstep, nxtrack, nmax_wavs], &
           norm_rad(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_o3p_diagnostic failed", &
           errstat)
    endif

  end subroutine read_o3p_diagnostic



  !> Read QA statistics parameters from an L2 netCDF O3 profile product
  !--------------------------------------------------------------------------
  !
  !> @param      tio_l2obj    L2 file object
  !> @param[in]  nstep        mirror step dimension size
  !> @param[in]  nxtrack      cross-track dimension size
  !> @param[in]  nwindow      fitting windows dimension size
  !> @param[in]  nmax_wavs    maximum wavelengths dimension size
  !> @param[in]  min_xtrack   min value of xtrack array to be read into
  !> @param[in]  max_xtrack   max value of xtrack array to be read into
  !> @param[in]  min_step     min value of step array to be read into
  !> @param[in]  max_step     max value of step array to be read into
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine read_o3p_qastat (tio_l2obj, nstep, nxtrack, nwindow, nmax_wavs, &
       min_xtrack, max_xtrack, min_step, max_step, errstat)

    use ozprof_data_module, only: ozwrtsnr
    use m_o3p_params, only: exval, iterations, rms, avg_resid, mqf, fit_wgt

    implicit none

    !input variables
    integer (kind=4), intent(in) :: nstep, nxtrack, nmax_wavs, nwindow, &
         min_xtrack, max_xtrack, min_step, max_step

    !output variables
    integer (kind=4), intent(inout) :: errstat

    type(tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_push_group (tio_l2obj, o3p_grp_qa_stats, errstat)
    call tiof_get2d_i4 (tio_l2obj, o3p_var_exit_stat, [0, 0], &
         [nstep, nxtrack], exval(min_xtrack:max_xtrack, min_step:max_step), &
         errstat)
    call tiof_get2d_i4 (tio_l2obj, o3p_var_iter, [0, 0], [nstep, nxtrack], &
         iterations(min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_fit_rms, [0, 0, 0], &
         [nstep, nxtrack, nwindow], &
         rms(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get3d_r4 (tio_l2obj, o3p_var_avg_resid, [0, 0, 0], &
         [nstep, nxtrack, nwindow], &
         avg_resid(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    call tiof_get1d_i2 (tio_l2obj, o3p_var_mqf, [0], [nstep], &
         mqf(min_step:max_step), errstat)
    ! Optional - fit weights
    if (ozwrtsnr .and. nmax_wavs > 0) then
      call tiof_get3d_r4 (tio_l2obj, o3p_var_fit_weight, [0, 0, 0], &
           [nstep, nxtrack, nmax_wavs], &
           fit_wgt(:, min_xtrack:max_xtrack, min_step:max_step), errstat)
    endif

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_o3p_qastat failed", &
           errstat)
    endif

  end subroutine read_o3p_qastat


  !> open an L2 netCDF O3 profile product file for reading
  !--------------------------------------------------------------------------
  !
  !> @param[in]  filename     name of input file
  !> @param      tio_l2obj    L2 file object
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan November 2016
  !--------------------------------------------------------------------------
  subroutine open_o3p(filename, tio_l2obj, errstat)

    implicit none

    character (len=*), intent(in) :: filename
    integer (kind=4), intent(inout) :: errstat
    type(tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_open(filename, tio_l2obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error(tell_io_open_error, &
           "failed to open input file for reading", errstat)
    endif

  end subroutine open_o3p


  !> close an L2 netCDF O3 profile product file after reading
  !--------------------------------------------------------------------------
  !
  !> @param      tio_l2obj    L2 file object
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan November 2016
  !--------------------------------------------------------------------------
  subroutine close_o3p(tio_l2obj, errstat)

    implicit none

    integer (kind=4), intent(inout) :: errstat
    type(tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_close(tio_l2obj, errstat)
    if (errstat /= 0) then
      call tell_error(tell_io_error, &
           "failed to close input file after reading", errstat)
    endif

  end subroutine close_o3p


end module m_read_L2_o3p_tio
