!> data parameters used in merge_o3p_files, and subroutine to allocate them
module m_o3p_params

  use tell_module

  implicit none


  !dimension index arrays
  integer (kind=4), dimension(:), allocatable :: nstep, nxtrack, ncorner, &
       nfitvars, nfitwins, ngas, nnongas, nlayer, nlayerp1, nmax_wavs, &
       nnoise_elems, naeros_wavs

  ! logical switches
  logical :: write_global_attr = .FALSE.

  ! output data parameter arrays
  ! allocatable to allow size to be set for full granule
  ! (nstep)
  real (kind=8), dimension(:), allocatable :: time
  integer (kind=2), dimension(:), allocatable :: mqf
  ! (ngas)
  character (len=64), dimension(:), allocatable :: gas_names
  ! (nnongas)
  character (len=64), dimension(:), allocatable :: nongas_names, nongas_units
  ! (nxtrack, nstep)
  integer (kind=2), dimension(:,:), allocatable :: geoflg
  real (kind=4), dimension(:,:), allocatable :: lat, lon, aza, sza, vza
  real (kind=4), dimension(:,:), allocatable :: o3tot, &
       o3tot_prec, o3tot_err, o3strat, o3strat_prec, o3strat_err, &
       o3trop, o3trop_prec, o3trop_err
  real (kind=4), dimension(:,:), allocatable :: aeros_idx, cld_frac, &
       cld_pres, glintprob, eff_alb, ozinfo, cld_opt_depth
  integer (kind=4), dimension(:,:), allocatable :: tropo_idx, &
         cld_flag, n_fit_wvl
  integer (kind=4), dimension(:,:), allocatable :: exval, iterations
  ! (ncorner, nxtrack, nstep)
  real (kind=4), dimension(:,:,:), allocatable :: corner_lat, corner_lon
  ! (nlayer, nxtrack, nstep)
  real (kind=4), dimension(:,:,:), allocatable :: ozprof, ozprof_prec, &
       ozprof_err
  real (kind=4), dimension(:,:,:), allocatable :: o3apriori, o3apriori_err
  ! (nlayer, nxtrack, nstep)
  real (kind=4), dimension(:,:,:), allocatable :: ozprof_pres, ozprof_alt, &
       ozprof_temp
  ! (2, nlayer, nxtrack, nstep)
  real (kind=4), dimension(:,:,:,:), allocatable :: ozprof_pres_bnds, ozprof_alt_bnds
  ! (ngas, nxtrack, nstep)
  real (kind=4), dimension(:,:,:), allocatable :: gas, gas_prec, gas_err
  real (kind=4), dimension(:,:,:), allocatable :: gas_apriori, &
       gas_apriori_err
  ! (nnongas, nxtrack, nstep)
  real (kind=4), dimension(:,:,:), allocatable :: nongas, nongas_prec, &
       nongas_err
  real (kind=4), dimension(:,:,:), allocatable :: nongas_apriori, &
       nongas_apriori_err
  ! (naeros_wavs, nxtrack, nstep)
  real (kind=4), dimension(:,:,:), allocatable :: aeros_opt_thick, &
       aeros_scatter_thick
  ! (nnoise_elems, nxtrack, nstep)
  integer (kind=2), dimension(:,:,:), allocatable :: noise_mtrx
  ! (nmax_wavs, nxtrack, nstep)
  real (kind=4), dimension(:,:,:), allocatable :: wavelengths, norm_rad, &
       sim_norm_rad
  real (kind=4), dimension(:,:,:), allocatable :: fit_wgt
  ! (nwindow, nxtrack, nstep)
  integer (kind=4), dimension(:,:,:), allocatable :: n_window_wvl
  real (kind=4), dimension(:,:,:), allocatable :: rms, avg_resid
  ! (nfitvars, nfitvars, nxtrack, nstep)
  real (kind=4), dimension(:,:,:,:), allocatable :: correl_mtrx
  ! (nfitvars, nmax_wavs, nxtrack, nstep)
  real (kind=4), dimension(:,:,:,:), allocatable :: wgt_func
  ! (nmax_wavs, nfitvars, nxtrack, nstep)
  real (kind=4), dimension(:,:,:,:), allocatable :: contrib_mtrx
  ! (nlayer, nlayer, nxtrack, nstep)
  integer (kind=2), dimension(:,:,:,:), allocatable :: avg_kernel


contains

  !> allocate data arrays for use in merging partial L2 ozone profile products
  !--------------------------------------------------------------------------
  !
  !> @param[out] nstep_tot    mirror step dimension size for full granule
  !! @param[out] nxtrack_tot  cross-track dimension size for full granule
  !! @param[out] ncorner      corner dimension size
  !! @param[out] nfitvars     fit_variables dimension size
  !! @param[out] nwindow      fitting_windows dimension size
  !! @param[out] ngas         gases dimension size (zero if unused)
  !! @param[out] nnongas      non_gas_variables dimension size
  !! @param[out] nlayer       layer dimension size
  !! @param[out] nmax_wavs    max_wavelengths dimension size (zero if unused)
  !! @param[out] nnoise_elems noise_elements dimension size (zero if unused)
  !! @param[out] naeros_wavs  aerosol_wavelengths dim. size (zero if unused)
  !! @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine o3p_param_alloc (min_step, max_step, min_xtrack, max_xtrack, ncorner, nfitvars, &
       nwindow, ngas, nnongas, nlayer, nmax_wavs, nnoise_elems, &
       naeros_wavs, errstat)

    use ozprof_data_module, only: ozwrtavgk, ozwrtcorr, ozwrtcovar, &
         ozwrtcontri, ozwrtres, ozwrtwf, ozwrtsnr, &
         do_lambcld
    use OMSAO_variables_module, only: reduce_resolution

    implicit none

    ! input variables
    integer (kind=4), intent(in) :: min_step, max_step, min_xtrack, &
         max_xtrack, ncorner, nfitvars, nwindow, ngas, nnongas, nlayer, &
         nmax_wavs, nnoise_elems, naeros_wavs
    ! output variables
    integer (kind=4), intent(inout) :: errstat


    if (errstat /= 0) return

    call tell_log(1,'allocating memory')

    ! allocate standard parameters
    allocate(time(min_step:max_step), &
         mqf(min_step:max_step), &
         geoflg(min_xtrack:max_xtrack, min_step:max_step), &
         lat(min_xtrack:max_xtrack, min_step:max_step), &
         lon(min_xtrack:max_xtrack, min_step:max_step), &
         aza(min_xtrack:max_xtrack, min_step:max_step), &
         sza(min_xtrack:max_xtrack, min_step:max_step), &
         vza(min_xtrack:max_xtrack, min_step:max_step), &
         o3tot(min_xtrack:max_xtrack, min_step:max_step), &
         o3tot_prec(min_xtrack:max_xtrack, min_step:max_step), &
         o3tot_err(min_xtrack:max_xtrack, min_step:max_step), &
         o3strat(min_xtrack:max_xtrack, min_step:max_step), &
         o3strat_prec(min_xtrack:max_xtrack, min_step:max_step), &
         o3strat_err(min_xtrack:max_xtrack, min_step:max_step), &
         o3trop(min_xtrack:max_xtrack, min_step:max_step), &
         o3trop_prec(min_xtrack:max_xtrack, min_step:max_step), &
         o3trop_err(min_xtrack:max_xtrack, min_step:max_step), &
         aeros_idx(min_xtrack:max_xtrack, min_step:max_step), &
         tropo_idx(min_xtrack:max_xtrack, min_step:max_step), &
         cld_pres(min_xtrack:max_xtrack, min_step:max_step), &
         cld_frac(min_xtrack:max_xtrack, min_step:max_step), &
         cld_flag(min_xtrack:max_xtrack, min_step:max_step), &
         glintprob(min_xtrack:max_xtrack, min_step:max_step), &
         eff_alb(min_xtrack:max_xtrack, min_step:max_step), &
         n_fit_wvl(min_xtrack:max_xtrack, min_step:max_step), &
         exval(min_xtrack:max_xtrack, min_step:max_step), &
         iterations(min_xtrack:max_xtrack, min_step:max_step), &
         corner_lat(ncorner, min_xtrack:max_xtrack, min_step:max_step), &
         corner_lon(ncorner, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof_prec(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof_err(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         o3apriori(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         o3apriori_err(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof_pres(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof_alt(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof_temp(nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof_pres_bnds(2,nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         ozprof_alt_bnds(2,nlayer, min_xtrack:max_xtrack, min_step:max_step), &
         n_window_wvl(nwindow, min_xtrack:max_xtrack, min_step:max_step), &
         rms(nwindow, min_xtrack:max_xtrack, min_step:max_step), &
         avg_resid(nwindow, min_xtrack:max_xtrack, min_step:max_step), &
         stat=errstat)

    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "o3p_param_alloc: failed to allocate standard parameters", errstat)
      return
    endif

    ! allocate optional parameters
    if (ngas > 0) then
      allocate(gas(ngas, min_xtrack:max_xtrack, min_step:max_step), &
           gas_prec(ngas, min_xtrack:max_xtrack, min_step:max_step), &
           gas_err(ngas, min_xtrack:max_xtrack, min_step:max_step), &
           gas_apriori(ngas, min_xtrack:max_xtrack, min_step:max_step), &
           gas_apriori_err(ngas, min_xtrack:max_xtrack, min_step:max_step), &
           gas_names(ngas), &
           stat = errstat)
    endif
    if (nnongas > 0) then
      allocate(nongas(nnongas, min_xtrack:max_xtrack, min_step:max_step), &
           nongas_prec(nnongas, min_xtrack:max_xtrack, min_step:max_step), &
           nongas_err(nnongas, min_xtrack:max_xtrack, min_step:max_step), &
           nongas_apriori(nnongas, min_xtrack:max_xtrack, min_step:max_step), &
           nongas_apriori_err(nnongas, min_xtrack:max_xtrack, min_step:max_step), &
           nongas_names(nnongas), &
           nongas_units(nnongas), &
           stat = errstat)
    endif
    if (nlayer > 0 .and. ozwrtavgk) then
      allocate(avg_kernel(nlayer, nlayer, min_xtrack:max_xtrack, min_step:max_step), &
           stat = errstat)
    endif
    if (nfitvars > 0 .and. ozwrtcorr) then
      allocate(correl_mtrx(nfitvars, nfitvars, min_xtrack:max_xtrack, min_step:max_step), &
           stat = errstat)
    endif
    if (nnoise_elems > 0 .and. ozwrtcovar) then
      allocate(noise_mtrx(nnoise_elems, min_xtrack:max_xtrack, min_step:max_step), &
           ozinfo(min_xtrack:max_xtrack, min_step:max_step), &
           stat = errstat)
    endif
    if (nmax_wavs > 0) then
      if (.not. reduce_resolution) &
           allocate(wavelengths(nmax_wavs, min_xtrack:max_xtrack, min_step:max_step), &
           stat=errstat)
      if (ozwrtres) allocate(norm_rad(nmax_wavs, min_xtrack:max_xtrack, min_step:max_step), &
           sim_norm_rad(nmax_wavs, min_xtrack:max_xtrack, min_step:max_step), &
           stat = errstat )
      if (ozwrtsnr) allocate(fit_wgt(nmax_wavs, min_xtrack:max_xtrack, min_step:max_step), &
           stat = errstat )
    endif
    if (naeros_wavs > 0) then
      allocate(aeros_opt_thick(naeros_wavs, min_xtrack:max_xtrack, min_step:max_step), &
           aeros_scatter_thick(naeros_wavs, min_xtrack:max_xtrack, min_step:max_step), &
           stat = errstat )
    endif
    if (nmax_wavs > 0 .and. nfitvars > 0) then
      if (ozwrtcontri) &
           allocate(contrib_mtrx(nmax_wavs, nfitvars, min_xtrack:max_xtrack, min_step:max_step),&
           stat = errstat )
      if (ozwrtwf) &
           allocate(wgt_func(nfitvars, nmax_wavs, min_xtrack:max_xtrack, min_step:max_step), &
           stat=errstat)
    endif
    !Cloud optical depth not goverened by whether optional dimension is used
    if (.not. do_lambcld) then
      allocate(cld_opt_depth(min_xtrack:max_xtrack, min_step:max_step), stat = errstat)
    endif

    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "o3p_param_alloc: failed to allocate optional parameters", errstat)
    endif

  end subroutine o3p_param_alloc





  !> deallocate data arrays used in merging partial L2 ozone profile products
  !--------------------------------------------------------------------------
  !
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan October 2016
  !--------------------------------------------------------------------------
  subroutine o3p_param_dealloc (errstat)

    implicit none

    ! output variables
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    call tell_log(1,'deallocating memory')

    ! allocate standard parameters
    if(allocated(time)) deallocate(time, stat = errstat)
    if(allocated(mqf)) deallocate(mqf , stat=errstat)
    if(allocated(geoflg)) deallocate(geoflg , stat=errstat)
    if(allocated(lat)) deallocate(lat , stat=errstat)
    if(allocated(lon)) deallocate(lon , stat=errstat)
    if(allocated(aza)) deallocate(aza , stat=errstat)
    if(allocated(sza)) deallocate(sza , stat=errstat)
    if(allocated(vza)) deallocate(vza , stat=errstat)
    if(allocated(o3tot)) deallocate(o3tot , stat=errstat)
    if(allocated(o3tot_prec)) deallocate(o3tot_prec , stat=errstat)
    if(allocated(o3tot_err)) deallocate(o3tot_err , stat=errstat)
    if(allocated(o3strat)) deallocate(o3strat , stat=errstat)
    if(allocated(o3strat_prec)) deallocate(o3strat_prec , stat=errstat)
    if(allocated(o3strat_err)) deallocate(o3strat_err , stat=errstat)
    if(allocated(o3trop)) deallocate(o3trop , stat=errstat)
    if(allocated(o3trop_prec)) deallocate(o3trop_prec , stat=errstat)
    if(allocated(o3trop_err)) deallocate(o3trop_err , stat=errstat)
    if(allocated(aeros_idx)) deallocate(aeros_idx , stat=errstat)
    if(allocated(tropo_idx)) deallocate(tropo_idx , stat=errstat)
    if(allocated(cld_pres)) deallocate(cld_pres , stat=errstat)
    if(allocated(cld_frac)) deallocate(cld_frac , stat=errstat)
    if(allocated(cld_flag)) deallocate(cld_flag , stat=errstat)
    if(allocated(glintprob)) deallocate(glintprob , stat=errstat)
    if(allocated(eff_alb)) deallocate(eff_alb , stat=errstat)
    if(allocated(n_fit_wvl)) deallocate(n_fit_wvl , stat=errstat)
    if(allocated(exval)) deallocate(exval , stat=errstat)
    if(allocated(iterations)) deallocate(iterations , stat=errstat)
    if(allocated(corner_lat)) deallocate(corner_lat , stat=errstat)
    if(allocated(corner_lon)) deallocate(corner_lon , stat=errstat)
    if(allocated(ozprof)) deallocate(ozprof , stat=errstat)
    if(allocated(ozprof_prec)) deallocate(ozprof_prec , stat=errstat)
    if(allocated(ozprof_err)) deallocate(ozprof_err , stat=errstat)
    if(allocated(o3apriori)) deallocate(o3apriori , stat=errstat)
    if(allocated(o3apriori_err)) deallocate(o3apriori_err , stat=errstat)
    if(allocated(ozprof_pres)) deallocate(ozprof_pres , stat=errstat)
    if(allocated(ozprof_alt)) deallocate(ozprof_alt , stat=errstat)
    if(allocated(ozprof_temp)) deallocate(ozprof_temp , stat=errstat)
    if(allocated(ozprof_pres_bnds)) deallocate(ozprof_pres_bnds , stat=errstat)
    if(allocated(ozprof_alt_bnds)) deallocate(ozprof_alt_bnds , stat=errstat)
    if(allocated(n_window_wvl)) deallocate(n_window_wvl , stat=errstat)
    if(allocated(rms)) deallocate(rms , stat=errstat)
    if(allocated(avg_resid)) deallocate(avg_resid , stat=errstat)

    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "o3p_param_dealloc: failed to deallocate standard parameters", &
           errstat)
      return
    endif

    ! deallocate optional parameters
    if(allocated(gas)) deallocate(gas, stat=errstat)
    if(allocated(gas_prec)) deallocate(gas_prec, stat=errstat)
    if(allocated(gas_err)) deallocate(gas_err , stat=errstat)
    if(allocated(gas_apriori)) deallocate(gas_apriori , stat=errstat)
    if(allocated(gas_apriori_err)) deallocate(gas_apriori_err , stat=errstat)
    if(allocated(gas_names)) deallocate(gas_names , stat=errstat)
    if(allocated(nongas)) deallocate(nongas , stat=errstat)
    if(allocated(nongas_prec)) deallocate(nongas_prec , stat=errstat)
    if(allocated(nongas_err)) deallocate(nongas_err , stat=errstat)
    if(allocated(nongas_apriori)) deallocate(nongas_apriori , stat=errstat)
    if(allocated(nongas_apriori_err)) deallocate(nongas_apriori_err , &
         stat=errstat)
    if(allocated(nongas_names)) deallocate(nongas_names , stat=errstat)
    if(allocated(nongas_units)) deallocate(nongas_units , stat=errstat)
    if(allocated(avg_kernel)) deallocate(avg_kernel , stat=errstat)
    if(allocated(correl_mtrx)) deallocate(correl_mtrx , stat=errstat)
    if(allocated(noise_mtrx)) deallocate(noise_mtrx , stat=errstat)
    if(allocated(ozinfo)) deallocate(ozinfo , stat=errstat)
    if(allocated(fit_wgt)) deallocate(fit_wgt , stat=errstat)
    if(allocated(wavelengths)) deallocate(wavelengths , stat=errstat)
    if(allocated(norm_rad)) deallocate(norm_rad , stat=errstat)
    if(allocated(sim_norm_rad)) deallocate(sim_norm_rad , stat=errstat)
    if(allocated(wgt_func)) deallocate(wgt_func , stat=errstat)
    if(allocated(aeros_opt_thick)) deallocate(aeros_opt_thick , stat=errstat)
    if(allocated(aeros_scatter_thick)) deallocate(aeros_scatter_thick , &
         stat=errstat)
    if(allocated(contrib_mtrx)) deallocate(contrib_mtrx , stat=errstat)
    if(allocated(cld_opt_depth)) deallocate(cld_opt_depth , stat=errstat)

    if (errstat /= 0) then
      call tell_error(tell_malloc_error, &
           "o3p_param_dealloc: failed to deallocate optional parameters", &
           errstat)
    endif

  end subroutine o3p_param_dealloc



  !> fill data arrays for merging partial L2 ozone profile products
  !--------------------------------------------------------------------------
  !
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan November 2016
  !--------------------------------------------------------------------------
  subroutine o3p_param_fill (errstat)

    implicit none

    ! output variables
    integer (kind=4), intent(inout) :: errstat
    ! internal variables
    real (kind=8) :: fill_double
    real (kind=4) :: fill_float
    real (kind=8), parameter :: fill_uint1 = 255
    !real (kind=8), parameter :: fill_int8 = -127
    !real (kind=8), parameter :: fill_uint8 = 65535 !-127
    real (kind=8), parameter :: fill_int16 = -32767
    !real (kind=8), parameter :: fill_uint16 = -32767
    integer (kind=4), external ::r8fill


    if (errstat /= 0) return

    if (r8fill(fill_double) /= 0) then
      call tell_error (tell_runtime_error, &
                       "o3p_param_fill: failed to define fill values", &
                       errstat)
      return
    endif
    fill_float = real(fill_double, kind=4)

    ! assign fill values to all basic variables
    time = fill_double
    geoflg = int(fill_uint1, kind=2)
    lat = fill_float
    lon = fill_float
    aza = fill_float
    sza = fill_float
    vza = fill_float
    corner_lat = fill_float
    corner_lon = fill_float
    mqf = 0   ! FIXME - should probably be flagged as bad if not included
    exval = int(fill_int16)
    ozprof = fill_float
    ozprof_prec = fill_float
    ozprof_err = fill_float
    ozprof_pres = fill_float
    ozprof_alt = fill_float
    ozprof_temp = fill_float
    ozprof_pres_bnds = fill_float
    ozprof_alt_bnds = fill_float
    o3tot = fill_float
    o3tot_prec = fill_float
    o3tot_err = fill_float
    o3strat = fill_float
    o3strat_prec = fill_float
    o3strat_err = fill_float
    o3trop = fill_float
    o3trop_prec = fill_float
    o3trop_err = fill_float
    o3apriori = fill_float
    o3apriori_err = fill_float
    aeros_idx = fill_float
    cld_frac = fill_float
    cld_pres = fill_float
    cld_flag = int(fill_int16)
    tropo_idx = int(fill_int16)
    glintprob = fill_float
    eff_alb = fill_float
    n_fit_wvl = int(fill_int16)
    n_window_wvl = int(fill_int16)
    iterations = int(fill_uint1)
    rms = fill_float
    avg_resid = fill_float

    ! assign fill values to optional variables
    if (allocated(gas)) then
      gas = fill_float
      gas_prec = fill_float
      gas_err = fill_float
      gas_apriori = fill_float
      gas_apriori_err = fill_float
    endif

    if (allocated(nongas)) then
      nongas = fill_float
      nongas_prec = fill_float
      nongas_err = fill_float
      nongas_apriori = fill_float
      nongas_apriori_err = fill_float
    endif

    if (allocated(avg_kernel)) avg_kernel = int(fill_int16, kind=2)

    if (allocated(correl_mtrx)) correl_mtrx = fill_float

    if (allocated(noise_mtrx)) then
      noise_mtrx = int(fill_int16, kind=2)
      ozinfo = fill_float
    endif

    if (allocated(contrib_mtrx)) contrib_mtrx = fill_float

    if (allocated(cld_opt_depth)) cld_opt_depth = fill_float

    if (allocated(aeros_opt_thick)) then
      aeros_opt_thick = fill_float
      aeros_scatter_thick = fill_float
    endif

    if (allocated(wavelengths)) wavelengths = fill_float

    if (allocated(wgt_func)) wgt_func = fill_float

    if (allocated(norm_rad)) then
      norm_rad = fill_float
      sim_norm_rad = fill_float
    endif

    if(allocated(fit_wgt)) fit_wgt = fill_float

  end subroutine o3p_param_fill



  !> allocate dimension size arrays used in merging L2 ozone profile products
  !--------------------------------------------------------------------------
  !
  !> @param[in]  ninputs      number of files being merged
  !! @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan November 2016
  !--------------------------------------------------------------------------
  subroutine o3p_dim_alloc (ninput, errstat)

    implicit none

    !input variables
    integer (kind = 4), intent(in) :: ninput
    ! output variables
    integer (kind = 4), intent(inout) :: errstat

    if (errstat /= 0) return

    allocate(nstep(ninput), nxtrack(ninput), ncorner(ninput), &
       nfitvars(ninput), nfitwins(ninput), ngas(ninput), &
       nnongas(ninput), nlayer(ninput), &
       nmax_wavs(ninput), nnoise_elems(ninput), naeros_wavs(ninput), &
       stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "o3p_dim_alloc: failed to allocate dimension size arrays", errstat)
    endif


  end subroutine o3p_dim_alloc



  !> deallocate dimension size arrays used in merging L2 ozone profile products
  !--------------------------------------------------------------------------
  !
  !> @param      errstat      error handling integer, non-zero = error
  !
  !> @author E. O'Sullivan November 2016
  !--------------------------------------------------------------------------
  subroutine o3p_dim_dealloc (errstat)

    implicit none

    ! output variables
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    if (allocated(nstep)) deallocate (nstep, stat=errstat)
    if (allocated(nxtrack)) deallocate (nxtrack, stat=errstat)
    if (allocated(ncorner)) deallocate (ncorner, stat=errstat)
    if (allocated(nfitvars)) deallocate (nfitvars, stat=errstat)
    if (allocated(nfitwins)) deallocate (nfitwins, stat=errstat)
    if (allocated(ngas)) deallocate (ngas, stat=errstat)
    if (allocated(nnongas)) deallocate (nnongas, stat=errstat)
    if (allocated(nlayer)) deallocate (nlayer, stat=errstat)
    if (allocated(nmax_wavs)) deallocate (nmax_wavs, stat=errstat)
    if (allocated(nnoise_elems)) deallocate (nnoise_elems, stat=errstat)
    if (allocated(naeros_wavs)) deallocate (naeros_wavs, stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "o3p_dim_alloc: failed to deallocate dimension size arrays", &
           errstat)
    endif


  end subroutine o3p_dim_dealloc


  end module m_o3p_params
