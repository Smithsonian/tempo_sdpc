!>Subroutines to write out L2 Ozone Profile netCDF file
module tio_output_module
  use netcdf
  use tio_module
  use tell_module
  use o3p_names_module
  use ozprof_data_module, only: ozwrtavgk, ozwrtcorr, ozwrtcovar, &
       ozwrtcontri, ozwrtres, ozwrtwf, ozwrtsnr, wrtring, &
       ozwrtvar, gaswrt, aerosol, do_lambcld
  use OMSAO_variables_module, only: reduce_resolution

  implicit none

  private write_coordinate_vars, append_geolocation_vars, &
       append_product_vars, append_support_vars, append_qa_vars, &
       append_diagnostic_vars
  public l2_tio_create, l2_tio_close, l2_tio_write_geo, l2_tio_write_data, &
       write_merged_geo, write_merged_data, copy_hdr_metadata, &
       copy_l2_metadata, label_output_file, set_production_date_time, &
       write_apriori_attr

  type (tiof_file_type), private, target :: primary_output_file

  !fill values to be filled in with r4fill, r8fill
  real (kind=8), private :: fill_float, fill_double
  !fill values taken from he5_l2writer_class
  real (kind=8), private, parameter :: fill_uint1 = 255
  real (kind=8), private, parameter :: fill_int16 = -32767
  real (kind=8), private, parameter :: fill_uint16 = 65535

  logical, private :: have_diagnostic_group = .false.

contains

  !> Create netCDF Level 2 product file
  !! @param[in] filename       output netCDF file name
  !! @param[in] first_pix      index of first xtrack pixel in output
  !! @param[in] last_pix       index of last xtrack pixel in output
  !! @param[in] first_line     index of first line in output
  !! @param[in] last_line      index of last line in output
  !! @param[in] ngas           number of other gas parameters (nfgas)
  !! @param[in] nlayer         number of layers in ozone profile (nlay)
  !! @param[in] nfitvar        number of fit variables (n_fitvar_rad)
  !! @param[in] nwindow        number of fit windows (numwin)
  !! @param[in] num_param      number of non-gas fit parameters
  !! @param[in] num_wav_max    maximum number of wavelengths in fit
  !! @param[in] step_indices   1D array of unbinned step indices
  !! @param[inout] errstat     error status variable
  subroutine l2_tio_create (filename, first_pix, last_pix, first_line, &
       last_line, ngas, nlayer, nfitvar, nwindow, &
       num_param, num_wav_max, step_indices, errstat)

    implicit none
    character (len=*), intent(in) :: filename
    integer (kind=4), intent(in) :: first_pix, last_pix, first_line, &
         last_line, ngas, nlayer, nfitvar, nwindow, num_param, num_wav_max
    integer (kind=4), dimension(:), intent(in) :: step_indices
    integer (kind=4), intent(inout) :: errstat
    integer (kind=4) :: num_steps, num_xtrack, num_elms, &
       num_aeros_wavl, num_corners, num_layer_bounds, bin_ms, n
    integer (kind=4), dimension(:), allocatable :: step_indices_bin

    integer (kind=4), external :: r8fill

    type (tiof_file_type), pointer :: obj
    type (tiof_dimlist_type) :: dimlist

    if (errstat < 0) return

    obj => primary_output_file

    if ((r8fill(fill_float) /= 0) .or. (r8fill(fill_double) /= 0)) then
      call tell_error (tell_runtime_error, &
                       "l2_tio_create: defining fill values", &
                       errstat)
      return
    endif

    num_corners = 4
    num_steps = last_line - first_line + 1
    num_xtrack = last_pix - first_pix  + 1
    num_aeros_wavl = nwindow + 2
    num_elms = (nlayer * (nlayer - 1))/2
    num_layer_bounds = 2

    ! get binned mirror step and xtrack indices
    allocate (step_indices_bin(num_steps), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "l2_tio_create: allocation error", &
           errstat)
      return
    endif
    bin_ms=size(step_indices)/num_steps
    do n=1,num_steps
      step_indices_bin(n)=INT(step_indices(n*bin_ms)/bin_ms)
    enddo

    ! Create a file.
    call tiof_create (obj, filename, nf90_clobber, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create: creating file "//trim(filename), &
                       errstat)
      return
    endif

    call tiof_put_git_commit_hash (obj, errstat)
    call tiof_write_epoch_timestamp (obj, errstat)

    ! Create default groups.
    call tiof_def_group (obj, o3p_grp_product, errstat)
    call tiof_def_group (obj, o3p_grp_geolocation, errstat)
    call tiof_def_group (obj, o3p_grp_support_data, errstat)
    call tiof_def_group (obj, o3p_grp_qa_stats, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_create:  defining groups in "//trim(filename), &
                       errstat)
      return
    endif

    !define dimension list
    call tiof_dimlist_append(dimlist, o3p_dim_xtrack, num_xtrack, errstat)
    call tiof_dimlist_append(dimlist, o3p_dim_step, num_steps, errstat)
    call tiof_dimlist_append(dimlist, o3p_dim_corner, num_corners, errstat)
    call tiof_dimlist_append(dimlist, o3p_dim_layer, nlayer, errstat)
    call tiof_dimlist_append(dimlist, o3p_dim_level, num_layer_bounds, errstat)
    call tiof_dimlist_append(dimlist, o3p_dim_fitvar, nfitvar, errstat)
    call tiof_dimlist_append(dimlist, o3p_dim_param, num_param, errstat)
    call tiof_dimlist_append(dimlist, o3p_dim_windows, nwindow, errstat)
    !dimensions for optional variables
    if (ozwrtcovar) then
      call tiof_dimlist_append(dimlist, o3p_dim_elms, num_elms, errstat)
    endif
    if (ozwrtres) then
      call tiof_dimlist_append(dimlist, o3p_dim_wavl_max, num_wav_max, &
           errstat)
    endif
    if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
      call tiof_dimlist_append(dimlist, o3p_dim_gas, ngas, errstat)
    endif
    if (aerosol) then
      call tiof_dimlist_append(dimlist, o3p_dim_aeros_wavl, num_aeros_wavl, &
           errstat)
    endif
    !define dimensions in output file
    call tiof_def_dims (obj, dimlist, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "l2_tio_create: defining dimensions in "//trim(filename), &
           errstat)
      return
    endif

    !write coordinate variables
    call write_coordinate_vars(obj, dimlist, first_pix, last_pix, num_steps, &
         num_xtrack, num_corners, nlayer, nfitvar, num_param, &
         nwindow, step_indices_bin, num_elms, num_wav_max, ngas, &
         num_aeros_wavl, errstat)

    !write geolocation group variables
    call tiof_push_group (obj, o3p_grp_geolocation, errstat)
    call append_geolocation_vars (obj, dimlist, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "l2_tio_create: creating geolocation group in "//trim(filename), &
           errstat)
      return
    endif

    !write product group variables
    call tiof_push_group (obj, o3p_grp_product, errstat)
    call append_product_vars (obj, dimlist, num_param, ngas, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "l2_tio_create: creating product group in "//trim(filename), &
           errstat)
      return
    endif

    !write support data group variables
    call tiof_push_group (obj, o3p_grp_support_data, errstat)
    call append_support_vars (obj, dimlist, num_param, ngas, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "l2_tio_create: creating support data group in "//trim(filename), &
           errstat)
      return
    endif

    !write diagnostic group variables

    call append_diagnostic_vars (obj, dimlist, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "l2_tio_create: creating diagnotic group in "//trim(filename), &
           errstat)
      return
    endif

    !write QA stats group variables
    call tiof_push_group (obj, o3p_grp_qa_stats, errstat)
    call append_qa_vars (obj, dimlist, errstat)
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "l2_tio_create: creating diagnotic group in "//trim(filename), &
           errstat)
      return
    endif

    !write metadata

    !before finishing, free the dimsension list, index arrays
    call tiof_dimlist_free(dimlist)
    deallocate (step_indices_bin, stat=errstat)

    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "l2_tio_create: creating file "//trim(filename), &
           errstat)
      return
    endif
  end subroutine l2_tio_create

  !> Close Level 2 product file
  !! @param[inout] errstat Error status variable
  subroutine l2_tio_close (errstat)
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    obj => primary_output_file
    call tiof_close (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_error, "l2_tio_close failed", errstat)
    endif

  end subroutine l2_tio_close

  !> Write coordinate variables to L2 netCDF file
  !! @param[inout]  obj        pointer to output file
  !! @param[in] dimlist        dimension list
  !! @param[in] first_pix      cross-track starting binned step number
  !! @param[in] last_pix       cross-track ending binned step number
  !! @param[in] num_steps      number of scan steps
  !! @param[in] num_xtrack     number of cross-track pixels
  !! @param[in] num_corners    number of pixel corners
  !! @param[in] num_layer      number of atmospheric layers
  !! @param[in] num_fitvar     number of fitted variables
  !! @param[in] num_param      number of non-gas fit variables
  !! @param[in] num_windows    number of fitting windows
  !! @param[in] step_indices_bin 1D arry of mirror step indices
  !! @param[in] num_elms       no. of elements in noise matrix (opt)
  !! @param[in] num_wav_max   number of max wavelength values (opt)
  !! @param[in] ngas        number of gasses in fit (opt)
  !! @param[in] num_aeros_wavl number of aerosol wavelengths (opt)
  !! @param[inout] errstat     error status variable
  subroutine write_coordinate_vars(obj, dimlist, first_pix, last_pix, &
       num_steps, num_xtrack, num_corners, num_layer, &
       num_fitvar, num_param,  num_windows, step_indices_bin, num_elms, &
       num_wav_max, ngas,  num_aeros_wavl, errstat)

    implicit none
    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer (kind=4), intent(in) :: num_steps, num_xtrack, num_corners, &
         num_layer, num_fitvar, num_param, num_windows, &
         num_elms, num_wav_max, ngas, num_aeros_wavl, first_pix, last_pix
    integer (kind=4), dimension(num_steps), intent(in) :: step_indices_bin
    integer (kind=4), intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    integer (kind=4), dimension(num_xtrack) :: xtrack_indices
    integer (kind=4), dimension(num_corners) :: corner_indices
    integer (kind=4), dimension(num_layer) :: layer_indices
    integer (kind=4), dimension(num_fitvar) :: fitvar_indices
    integer (kind=4), dimension(num_param) :: param_indices
    integer (kind=4), dimension(num_windows) :: window_indices
    integer (kind=4), dimension(:), allocatable :: elm_indices
    integer (kind=4), dimension(:), allocatable :: wav_max_indices
    integer (kind=4), dimension(:), allocatable :: gas_indices
    integer (kind=4), dimension(:), allocatable :: aeros_wavl_indices
    integer :: i, dimids(12), status

    if (errstat < 0) return

    if (.false.) then
      !If optional variables in use, allocate indices
      if (ozwrtcovar) then
        allocate(elm_indices(num_elms), stat=status)
        if (status /= 0) then
          call tell_error (tell_malloc_error, &
                           'write_coordinate_vars:  allocate failed', errstat)
          return
        endif
      endif

      if (ozwrtres) then
        allocate(wav_max_indices(num_wav_max), stat=status)
        if (status /= 0) then
          call tell_error (tell_malloc_error, &
                           'write_coordinate_vars:  allocate failed', errstat)
          return
        endif
      endif

      if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
        allocate(gas_indices(ngas), stat=status)
        if (status /= 0) then
          call tell_error (tell_malloc_error, &
                           'write_coordinate_vars:  allocate failed', errstat)
          return
        endif
      endif

      if (aerosol) then
        allocate(aeros_wavl_indices(num_aeros_wavl), stat=status)
        if (status /= 0) then
          call tell_error (tell_malloc_error, &
                           'write_coordinate_vars:  allocate failed', errstat)
          return
        endif
      endif
    endif

    ! Retrieve dimension info for baseline coordinates
    call tiof_dimlist_lookup(dimlist, [o3p_dim_xtrack, o3p_dim_step, &
         o3p_dim_corner, o3p_dim_layer, o3p_dim_level, o3p_dim_fitvar, &
         o3p_dim_param, o3p_dim_windows], dimids(1:8), errstat)
    ! Make the baseline variable list
    call tiof_varlist_append(varlist, errstat, o3p_dim_xtrack, nf90_int, &
                             dimids=[dimids(1)], &
                             long_name = "pixel index along slit")
    call tiof_varlist_append (varlist, errstat, o3p_dim_step, nf90_int, &
                              dimids=[dimids(2)], &
                              long_name = "scan mirror position index")
    if (.false.) then
      call tiof_varlist_append (varlist, errstat, o3p_dim_corner, nf90_int, &
                                dimids=[dimids(3)])
      call tiof_varlist_append (varlist, errstat, o3p_dim_layer, nf90_int, &
                                dimids=[dimids(4)])
      call tiof_varlist_append (varlist, errstat, o3p_dim_level, nf90_int, &
                                dimids=[dimids(5)])
      call tiof_varlist_append (varlist, errstat, o3p_dim_fitvar, nf90_int, &
                                dimids=[dimids(6)])
      call tiof_varlist_append (varlist, errstat, o3p_dim_param, nf90_int, &
                                dimids=[dimids(7)])
      call tiof_varlist_append (varlist, errstat, o3p_dim_windows, nf90_int, &
                                dimids=[dimids(8)])
      ! Add optional coordinates
      if (ozwrtcovar) then
        call tiof_dimlist_lookup(dimlist, [o3p_dim_elms], dimids(9:9), errstat)
        call tiof_varlist_append (varlist, errstat, o3p_dim_elms, nf90_int, &
                                  dimids=[dimids(9)])
      endif
      if (ozwrtres) then
        call tiof_dimlist_lookup(dimlist, [o3p_dim_wavl_max], dimids(10:10), &
           errstat)
        call tiof_varlist_append (varlist, errstat, o3p_dim_wavl_max, nf90_int, &
                                  dimids=[dimids(10)])
      endif
      if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
        call tiof_dimlist_lookup(dimlist, [o3p_dim_gas], dimids(11:11), errstat)
        call tiof_varlist_append (varlist, errstat, o3p_dim_gas, nf90_int, &
                                  dimids=[dimids(11)])
      endif
      if (aerosol) then
        call tiof_dimlist_lookup(dimlist, [o3p_dim_aeros_wavl], dimids(12:12), &
                                 errstat)
        call tiof_varlist_append (varlist, errstat, o3p_dim_aeros_wavl, &
                                  nf90_int, dimids=[dimids(12)])
      endif
    endif

    !Write variable list to L2 file
    call tiof_def_vars(obj, varlist, errstat)
    call tiof_varlist_free (varlist)

    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "write_coordinate_vars: error writing variable list", &
           errstat)
      return
    endif

    ! Use input mirror step indices, and set xtrack indices to binned pixel
    ! numbers so as to allow stitching together of output files when code
    ! is run on subsets of the granule
    !Write variable indices to L2 file
    call tiof_put1d_i4 (obj, o3p_dim_step, [0], [num_steps], &
         step_indices_bin, errstat)

    xtrack_indices = [(i, i=first_pix-1,last_pix-1)]
    call tiof_put1d_i4 (obj, o3p_dim_xtrack, [0], [num_xtrack], &
         xtrack_indices, errstat)

    if (.false.) then
      corner_indices = [(i, i=0,num_corners-1)]
      call tiof_put1d_i4 (obj, o3p_dim_corner, [0], [num_corners], &
                          corner_indices, errstat)

      layer_indices = [(i, i=0,num_layer-1)]
      call tiof_put1d_i4 (obj, o3p_dim_layer, [0], [num_layer], &
                          layer_indices, errstat)

      fitvar_indices = [(i, i=0,num_fitvar-1)]
      call tiof_put1d_i4 (obj, o3p_dim_fitvar, [0], [num_fitvar], &
                          fitvar_indices, errstat)

      param_indices = [(i, i=0,num_param-1)]
      call tiof_put1d_i4 (obj, o3p_dim_param, [0], [num_param], &
                          param_indices, errstat)

      window_indices = [(i, i=0,num_windows-1)]
      call tiof_put1d_i4 (obj, o3p_dim_windows, [0], [num_windows], &
                          window_indices, errstat)

      if (ozwrtcovar) then
        elm_indices = [(i, i=0,num_elms-1)]
        call tiof_put1d_i4 (obj, o3p_dim_elms, [0], [num_elms], &
                            elm_indices, errstat)
      endif

      if (ozwrtres) then
        wav_max_indices = [(i, i=0,num_wav_max-1)]
        call tiof_put1d_i4 (obj, o3p_dim_wavl_max, [0], [num_wav_max], &
                            wav_max_indices, errstat)
      endif

      if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
        gas_indices = [(i, i=0,ngas-1)]
        call tiof_put1d_i4 (obj, o3p_dim_gas, [0], [ngas], &
                            gas_indices, errstat)
      endif

      if (aerosol) then
        aeros_wavl_indices = [(i, i=0,num_aeros_wavl-1)]
        call tiof_put1d_i4 (obj, o3p_dim_aeros_wavl, [0], [num_aeros_wavl], &
                            aeros_wavl_indices, errstat)
      endif
    endif

    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
           "write_coordinate_vars: error writing variable indices", &
           errstat)
      return
    endif

    if (.false.) then
      !tidy up any allocated arrays
      if (ozwrtcovar) then
        deallocate(elm_indices, stat=errstat)
      endif

      if (ozwrtres) then
        deallocate(wav_max_indices, stat=errstat)
      endif

      if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
        deallocate(gas_indices, stat=errstat)
      endif

      if (aerosol) then
        deallocate(aeros_wavl_indices, stat=errstat)
      endif
    endif

  end subroutine write_coordinate_vars

  !> Define geolocation variables in L2 output file
  !! @param[in]    obj      pointer to output file
  !! @param[in]    dimlist  dimension list
  !! @param[inout] errstat  error status variable
  subroutine append_geolocation_vars(obj, dimlist, errstat)

    implicit none

    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord, att_lonbnd, att_latbnd, att_time
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_corner_xtrack_step
    character (len=32) :: epoch_buf

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_xtrack, o3p_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_corner, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_corner_xtrack_step, &
                              errstat)
    ! Construct a list of variables with their associated dimension ids
    ! and attributes:
    epoch_buf(:)=''
    call tiof_mktimestamp_str (0.0_8, epoch_buf, errstat)

    call tiof_attlist_append (att_time, errstat, "calendar", &
                              att_text = "gregorian")
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_time, &
                              nf90_double, &
                              dimids = [dimids_xtrack_step(2)],  &
                              standard_name = "time", &
                              comment = "exposure start time", &
                              units = "seconds since "//trim(epoch_buf), &
                              valid_range = [-5.0e9_8, 1.e10_8], &
                              fillvalue = fill_double, &
                              attlist = att_time)
    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude))
    call tiof_attlist_append (att_latbnd, errstat, "bounds", &
                              att_text = o3p_var_latitude_bounds)
    call tiof_attlist_append (att_lonbnd, errstat, "bounds", &
                              att_text = o3p_var_longitude_bounds)
    call tiof_attlist_append (att_latbnd, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude))
    call tiof_attlist_append (att_lonbnd, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude))
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_latitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              standard_name = "latitude", &
                              comment = "latitude at pixel center", &
                              units = "degrees_north", &
                              valid_range = [-90.0_8, 90.0_8], &
                              fillvalue = fill_float, &
                              attlist = att_latbnd)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_latitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              comment = "latitude at pixel corners", &
                              !units = "degrees_north", &
                              valid_range = [-90.0_8, 90.0_8], &
                              fillvalue = fill_float, &
                              attlist = att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_longitude, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              standard_name = "longitude", &
                              comment = "longitude at pixel center", &
                              units = "degrees_east", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist = att_lonbnd)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_longitude_bounds, &
                              nf90_float, &
                              dimids = dimids_corner_xtrack_step,  &
                              comment = "longitude at pixel corners", &
                              !units = "degrees_east", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist = att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_sz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "solar zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_vz_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "viewing zenith angle at pixel center", &
                              units = "degrees", &
                              valid_range = [0.0_8, 90.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_ra_angle, &
                              nf90_float, &
                              dimids = dimids_xtrack_step,  &
                              comment = "relative azimuth angle at pixel center", &
                              units = "degrees", &
                              valid_range = [-180.0_8, 180.0_8], &
                              fillvalue = fill_float, &
                              attlist=att_coord)
!    call tiof_varlist_append (varlist, errstat, &
!                              o3p_var_anomflg, &
!                              nf90_ubyte, &
!                              dimids = dimids_xtrack_step, &
!                              comment = "cross-track quality flag", &
!                              valid_range = [0.0_8, 254.0_8], &
!                              attlist=att_coord)

    call tiof_def_vars(obj,varlist,errstat)
    call tiof_varlist_free(varlist)
    call tiof_attlist_free(att_coord)
    call tiof_attlist_free(att_lonbnd)
    call tiof_attlist_free(att_latbnd)
    call tiof_attlist_free(att_time)

  end subroutine append_geolocation_vars

  !> Define product variables in L2 output file
  !! @param[in]    obj        pointer to output file
  !! @param[in]    dimlist    dimension list
  !! @param[in]    num_param  number of other fit parameters
  !! @param[in]    ngas    number of fitted other gas parameters
  !! @param[inout] errstat    error status variable
!  subroutine append_product_vars(obj, dimlist, num_param, othunitc, errstat)
  subroutine append_product_vars(obj, dimlist, num_param, ngas, errstat)

    implicit none

    ! input
    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(in) :: num_param, ngas
    ! output
    integer, intent(inout) :: errstat
    !local
    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord, att_coord_layer
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_layer_xtrack_step, &
         dimids_param_xtrack_step, dimids_gas_xtrack_step
    integer, parameter :: deflate_level = 1
    logical, parameter :: shuffle = .true.
    ! In case you want to write out the as-yet unused sub2km fields
    logical :: do_sub2km = .FALSE.

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_xtrack, o3p_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_layer, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_layer_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_param, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_param_xtrack_step, &
                              errstat)
    if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_gas, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_gas_xtrack_step, &
                              errstat)
    endif

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude))
    call tiof_attlist_append (att_coord_layer, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude) &
                              //' '//trim(o3p_var_profile_pres) &
                              //' '//trim(o3p_var_profile_alt))

    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_retrieve_prof, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "retrieved ozone profile", &
                              units = "DU", &
                              valid_range = [-100.0_8, 100.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_retrieve_prof_prec, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "retrieved ozone profile precision", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_retrieve_prof_err, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "retrieved ozone profile solution error", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_total_o3, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "total ozone column", &
                              units = "DU", &
                              valid_range = [0.0_8, 600.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_total_o3_prec, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "total ozone column precision", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_total_o3_err, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "total ozone column solution error", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_strat_o3, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "stratospheric ozone column", &
                              units = "DU", &
                              valid_range = [0.0_8, 600.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_strat_o3_prec, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "stratospheric ozone column precision", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_strat_o3_err, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "stratospheric ozone column solution error", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_tropo_o3, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "tropospheric ozone column", &
                              units = "DU", &
                              valid_range = [0.0_8, 100.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_tropo_o3_prec, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "tropospheric ozone column precision", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_tropo_o3_err, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "tropospheric ozone column solution error", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    ! Not in current algorithm, but expected for UVIS version
    ! Set do_sub2km = .TRUE. if you need to write (empty) fields
    ! to output file
    if (do_sub2km) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_surf_o3, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "0-2km ozone column", &
                              units = "DU", &
                              valid_range = [0.0_8, 600.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_surf_o3_prec, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "0-2km ozone column precision", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_surf_o3_err, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "0-2km ozone column solution error", &
                              units = "DU", &
                              valid_range = [0.0_8, 50.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif
    ! Optional products

    ! Other fitted gases
    if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_other_gas_retrieve, &
                              nf90_float, &
                              dimids = dimids_gas_xtrack_step, &
                              comment = "other gas retrieved vertical column density", &
                              units = "molecules cm^-2", &
                              valid_range = [-1.0E30_8, 1.0E30_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_other_gas_retrieve_prec, &
                              nf90_float, &
                              dimids = dimids_gas_xtrack_step, &
                              comment = "other gas retrieved vertical column density precision", &
                              units = "molecules cm^-2", &
                              valid_range = [0.0_8, 1.0E30_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_other_gas_retrieve_err, &
                              nf90_float, &
                              dimids = dimids_gas_xtrack_step, &
                              comment = "other gas retrieved vertical column density solution error", &
                              units = "molecules cm^-2", &
                              valid_range = [0.0_8, 1.0E30_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif ! other fitted gases

    ! other variables
    if (num_param > 0 .and. ozwrtvar) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_nongas_param_ret, &
                              nf90_float, &
                              dimids = dimids_param_xtrack_step, &
                              comment = "non-gas parameter retrieved values", &
                              units = "non-gas parameter units", &
                              valid_range = [-1.0E30_8, 1.0E30_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_nongas_param_ret_prec, &
                              nf90_float, &
                              dimids = dimids_param_xtrack_step, &
                              comment = "non-gas parameter precision", &
                              units = "non-gas parameter units", &
                              valid_range = [0.0_8, 1.0E30_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_nongas_param_ret_err, &
                              nf90_float, &
                              dimids = dimids_param_xtrack_step, &
                              comment = "non-gas parameter solution error", &
                              units = "non-gas parameter units", &
                              valid_range = [0.0_8, 1.0E30_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif ! other non-gas parameters

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_coord)
    call tiof_attlist_free (att_coord_layer)

  end subroutine append_product_vars

  !> Define support variables in L2 output file
  !! @param[in]    obj        pointer to output file
  !! @param[in]    dimlist    dimension list
  !! @param[in]    num_param  number of other fit parameters
  !! @param[in]    ngas    number of other gas parameters
  !! @param[inout] errstat    error status variable
!  subroutine append_support_vars(obj, dimlist, num_param, othunitc, errstat)
  subroutine append_support_vars(obj, dimlist, num_param, ngas, errstat)

    implicit none

    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(in) :: num_param, ngas
!    character (len=*), intent(in) :: othunitc
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord, att_coord_layer
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_layer_xtrack_step, &
         dimids_param_xtrack_step, dimids_gas_xtrack_step, &
         dimids_elms_xtrack_step, &
         dimids_wavmax_xtrack_step, dimids_aeros_xtrack_step, &
         dimids_window_xtrack_step
    integer, dimension(4) :: dimids_layer_layer_xtrack_step, &
         dimids_level_layer_xtrack_step, &
         dimids_wavmax_fitvar_xtrack_step, dimids_fitvar_fitvar_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_xtrack, o3p_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_layer, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_layer_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_level, o3p_dim_layer, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_level_layer_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_layer, o3p_dim_layer, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_layer_layer_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_param, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_param_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_fitvar, o3p_dim_fitvar, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_fitvar_fitvar_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_windows, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_window_xtrack_step, &
                              errstat)
    if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_gas, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_gas_xtrack_step, &
                              errstat)
    endif
    if (ozwrtcovar) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_elms, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_elms_xtrack_step, &
                              errstat)
    endif
    if (ozwrtres) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_wavl_max, o3p_dim_fitvar, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_wavmax_fitvar_xtrack_step, &
                              errstat)
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_wavl_max, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_wavmax_xtrack_step, &
                              errstat)
    endif
    if (aerosol) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_aeros_wavl, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_aeros_xtrack_step, &
                              errstat)
    endif

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude))
    call tiof_attlist_append (att_coord_layer, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude) &
                              //' '//trim(o3p_var_profile_pres) &
                              //' '//trim(o3p_var_profile_alt))
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_geoflg, &
                              nf90_int, &
                              dimids = dimids_xtrack_step, &
                              comment = "ground pixel quality flag", &
                              valid_range = [0.0_8, 65534.0_8], &
                              fillvalue = fill_uint16, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_aeros_index, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "UV aerosol index", &
                              valid_range = [-200.0_8, 200.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_cld_flag, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "cloud flag", &
                              valid_range = [0.0_8, 4.0_8], &
                              fillvalue = fill_int16, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_eff_cld_frac, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "effective cloud fraction", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_eff_cld_pres, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "effective cloud pressure", &
                              valid_range = [0.0_8, 1100.0_8], &
                              units = "hPa", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_glint_prob, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "glint probablility", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_apriori_prof, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "a priori ozone profile", &
                              valid_range = [0.0_8, 100.0_8], &
                              units = "DU", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_apriori_prof_err, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "a priori ozone profile error", &
                              valid_range = [0.0_8, 100.0_8], &
                              units = "DU", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_profile_alt, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "altitude of each retrieval layer", &
                              valid_range = [0.0_8, 100.0_8], &
                              units = "km", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_profile_alt_bnds, &
                              nf90_float, &
                              dimids = dimids_level_layer_xtrack_step, &
                              comment = "altitude of each retrieval layer boundary", &
                              valid_range = [0.0_8, 100.0_8], &
                              !units = "km", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_profile_temp, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "temperature of each retrieval layer", &
                              valid_range = [150.0_8, 350.0_8], &
                              units = "Kelvin", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_profile_pres, &
                              nf90_float, &
                              dimids = dimids_layer_xtrack_step, &
                              comment = "pressure of each retrieval layer", &
                              valid_range = [0.0_8, 1100.0_8], &
                              units = "hPa", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_profile_pres_bnds, &
                              nf90_float, &
                              dimids = dimids_level_layer_xtrack_step, &
                              comment = "pressure of each retrieval layer boundary", &
                              valid_range = [0.0_8, 1100.0_8], &
                              !units = "hPa", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_tropo_index, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "tropopause index", &
                              valid_range = [0.0_8, 100.0_8], &
                              fillvalue = fill_int16, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_surf_albedo, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "surface albedo", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_fit_wavel, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "number of wavelengths used in fitting", &
                              valid_range = [0.0_8, 2000.0_8], &
                              fillvalue = fill_int16, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_window_wavel, &
                              nf90_short, &
                              dimids = dimids_window_xtrack_step, &
                              comment = "number of wavelengths in each fitting window", &
                              valid_range = [0.0_8, 1000.0_8], &
                              fillvalue = fill_int16, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)

    ! Optional parameters
    ! Other fitted gases
    if (ngas > 0 .and. (gaswrt .or. ozwrtvar)) then
      ! OMI code writes variable names as a very long string in title field
      ! We will write them as a separate array of strings
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_other_gas_names, &
                              nf90_string, &
                              dimids = [dimids_gas_xtrack_step(1)], &
                              comment = "names of other fitted gases")
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_other_gas_apriori, &
                              nf90_float, &
                              dimids = dimids_gas_xtrack_step, &
                              comment = "other gas a priori column density", &
                              valid_range = [0.0_8, 1.0E30_8], &
                              units = "molecules cm^-2", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_other_gas_apriori_err, &
                              nf90_float, &
                              dimids = dimids_gas_xtrack_step, &
                              comment = "other gas a priori column density error", &
                              valid_range = [-1.0E30_8, 1.0E30_8], &
                              units = "molecules cm^-2", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif !other fitted gases

    ! non-gas fitted parameters
    if (num_param > 0 .and. ozwrtvar) then
      ! OMI code writes variable names as a very long string in title field
      ! We will write them as a separate array of strings
      ! Units also need their own array of strings
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_nongas_param_names, &
                              nf90_string, &
                              dimids = [dimids_param_xtrack_step(1)], &
                              comment = "names of non-gas fitted parameters")
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_nongas_param_units, &
                              nf90_string, &
                              dimids = [dimids_param_xtrack_step(1)], &
                              comment = "units of non-gas fitted parameters")
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_nongas_param_apriori, &
                              nf90_float, &
                              dimids = dimids_param_xtrack_step, &
                              comment = "non-gas parameter a priori", &
                              valid_range = [0.0_8, 1.0E30_8], &
                              units = "non-gas parameter units", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_nongas_param_apriori_err, &
                              nf90_float, &
                              dimids = dimids_param_xtrack_step, &
                              comment = "non-gas parameter a priori error", &
                              valid_range = [-1.0E30_8, 1.0E30_8], &
                              units = "non-gas parameter units", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif ! non-gas params

    ! Ozone averaging kernels
    if (ozwrtavgk) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_avg_kernel, &
                              nf90_short, &
                              dimids = dimids_layer_layer_xtrack_step, &
                              comment = "ozone profile averaging kernels", &
                              valid_range = [-32766.0_8, 32767.0_8], &
                              !units = "DU", &
                              fillvalue = fill_int16, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord_layer)

    endif

    ! Ozone correlation matrix
    if (ozwrtcorr) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_correl, &
                              nf90_float, &
                              dimids = dimids_fitvar_fitvar_xtrack_step, &
                              comment = "correlation matrix (upper matrix above diagonal)", &
                              valid_range = [-1.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    ! Ozone noise matrix
    if (ozwrtcovar) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_noise_matrix, &
                              nf90_short, &
                              dimids = dimids_elms_xtrack_step, &
                              comment = "O3 noise correlation matrix", &
                              valid_range = [-10000.0_8, 10000.0_8], &
                              !units = "DU", &
                              fillvalue = fill_int16, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    ! FIXME - information content should have its own switch rather than
    !   using the covariance matrix switch, but for now this is OK.
    if (ozwrtcovar) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_o3_info_content, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "ozone information content", &
                              valid_range = [0.0_8, 100.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif
    ! Contibution matrix
    if (ozwrtcontri) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_contrib_func, &
                              nf90_float, &
                              dimids = dimids_wavmax_fitvar_xtrack_step, &
                              comment = "contribution function matrix", &
                              valid_range = [-1.0E38_8, 1.0E38_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    ! Cloud optical depth
    if (.not. do_lambcld) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_cld_opt_depth, &
                              nf90_float, &
                              dimids = dimids_xtrack_step, &
                              comment = "effective cloud optical depth", &
                              valid_range = [0.0_8, 1000.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    ! Aerosols
    if (aerosol) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_aeros_opt_thick, &
                              nf90_float, &
                              dimids = dimids_aeros_xtrack_step, &
                              comment = "input aerosol optical thickness", &
                              valid_range = [0.0_8, 20.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_aeros_scatter_thick, &
                              nf90_float, &
                              dimids = dimids_aeros_xtrack_step, &
                              comment = "aerosol scattering optical thickness", &
                              valid_range = [0.0_8, 20.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_coord)
    call tiof_attlist_free (att_coord_layer)

  end subroutine append_support_vars

  !> Define support variables in L2 output file
  !! @param[in]    obj        pointer to output file
  !! @param[in]    dimlist    dimension list
  !! @param[inout] errstat    error status variable
  subroutine append_diagnostic_vars (obj, dimlist, errstat)

    implicit none

    type (tiof_file_type), intent(inout) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_wavmax_xtrack_step, &
         dimids_aeros_xtrack_step
    integer, dimension(4) :: dimids_fitvar_wavmax_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_xtrack, o3p_dim_step], &
                              dimids_xtrack_step, &
                              errstat)

    if (ozwrtres) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_wavl_max, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_wavmax_xtrack_step, &
                              errstat)

    endif
    if (aerosol) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_aeros_wavl, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_aeros_xtrack_step, &
                              errstat)
    endif
    if (ozwrtwf) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_fitvar, o3p_dim_wavl_max, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_fitvar_wavmax_xtrack_step, &
                              errstat)
    endif

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude))
    ! Optional variables
    ! fitted wavelengths
    if ( (.not. reduce_resolution) .and. &
         (ozwrtcontri .or. ozwrtres .or. ozwrtwf .or. ozwrtsnr .or. &
         wrtring)) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_wavel, &
                              nf90_float, &
                              dimids = dimids_wavmax_xtrack_step, &
                              comment = "wavelength used in the retrieval", &
                              valid_range = [100.0_8, 1000.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    ! Ring-effect spectrum
    if (wrtring) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_ring_spec, &
                              nf90_float, &
                              dimids = dimids_wavmax_xtrack_step, &
                              comment = "Ring effect spectrum, relative filling-in", &
                              valid_range = [-1.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      endif

    ! Weighting function
    if (ozwrtwf) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_weight_func, &
                              nf90_float, &
                              dimids = dimids_fitvar_wavmax_xtrack_step, &
                              comment = "weighting function matrix", &
                              valid_range = [-1.0E38_8, 1.0E38_8], &
                              units = "1/unit of respective variables", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif
    ! observed and simulated spectra
    if (ozwrtres) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_norm_radiance, &
                              nf90_float, &
                              dimids = dimids_wavmax_xtrack_step, &
                              comment = "observed normalized radiance, I/F", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_sim_norm_rad, &
                              nf90_float, &
                              dimids = dimids_wavmax_xtrack_step, &
                              comment = "simulated normalized radiance, I/F", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    ! aerosol wavelengths
    if (aerosol) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_aeros_wavel, &
                              nf90_float, &
                              dimids = dimids_aeros_xtrack_step, &
                              comment = "input aerosol wavelengths", &
                              valid_range = [100.0_8, 1000.0_8], &
                              units="nm", &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    if (varlist % num_items > 0) then
      call tiof_def_group (obj, o3p_grp_diagnostic, errstat)
      call tiof_push_group (obj, o3p_grp_diagnostic, errstat)
      call tiof_def_vars (obj, varlist, errstat)
      call tiof_pop_group (obj, errstat)
      have_diagnostic_group = (errstat == 0)
    endif
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_coord)
  end subroutine append_diagnostic_vars

  !> Define support variables in L2 output file
  !! @param[in]    obj        pointer to output file
  !! @param[in]    dimlist    dimension list
  !! @param[inout] errstat    error status variable
  subroutine append_qa_vars(obj, dimlist, errstat)

    use OMSAO_variables_module, only : max_itnum_rad
    implicit none

    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: dimlist
    integer, intent(inout) :: errstat

    type (tiof_varlist_type) :: varlist
    type (tiof_attlist_type) :: att_coord
    integer, dimension(2) :: dimids_xtrack_step
    integer, dimension(3) :: dimids_window_xtrack_step, &
         dimids_wavmax_xtrack_step
    integer, parameter :: deflate_level = 5
    logical, parameter :: shuffle = .true.

    if (errstat < 0) return

    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_xtrack, o3p_dim_step], &
                              dimids_xtrack_step, &
                              errstat)
    call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_windows, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_window_xtrack_step, &
                              errstat)
    if (ozwrtsnr) then
      call tiof_dimlist_lookup (dimlist, &
                              [o3p_dim_wavl_max, o3p_dim_xtrack, o3p_dim_step], &
                              dimids_wavmax_xtrack_step, &
                              errstat)
    endif

    ! Construct a list of variables with their associated dimension ids
    ! and attributes:

    call tiof_attlist_append (att_coord, errstat, "coordinates", &
                              att_text = trim(o3p_var_time) &
                              //' '//trim(o3p_var_longitude) &
                              //' '//trim(o3p_var_latitude))
    !call tiof_varlist_append (varlist, errstat, &
    !                          o3p_var_mqf, &
    !                          nf90_ushort, &
    !                          dimids = [dimids_xtrack_step(2)], &
    !                          comment = "measurement quality flags", &
    !                          valid_range = [0.0_8, 254.0_8], &
    !                          fillvalue = fill_uint1)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_exit_stat, &
                              nf90_short, &
                              dimids = dimids_xtrack_step, &
                              comment = "retrieval exit status (0=not converged, 1-99 converged, "// &
                              "100-110 converged but with negative ozone value)", &
                              valid_range = [-10.0_8, 110.0_8], &
                              fillvalue = fill_int16, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_iter, &
                              nf90_ushort, &
                              dimids = dimids_xtrack_step, &
                              comment = "number of iterations", &
                              valid_range = [0.0_8, max_itnum_rad * 1.0_8], &
                              fillvalue = fill_uint1, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_fit_rms, &
                              nf90_float, &
                              dimids = dimids_window_xtrack_step, &
                              comment = "ratio of fitting residual to measurement error", &
                              valid_range = [0.0_8, 100.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    call tiof_varlist_append (varlist, errstat, &
                              o3p_var_avg_resid, &
                              nf90_float, &
                              dimids = dimids_window_xtrack_step, &
                              comment = "average fitting residuals", &
                              units = "percent", &
                              valid_range = [0.0_8, 100.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)

    !optional parameter - relative measurement error
    if (ozwrtsnr) then
      call tiof_varlist_append (varlist, errstat, &
                              o3p_var_merr, &
                              nf90_float, &
                              dimids = dimids_wavmax_xtrack_step, &
                              comment = "relative measurement error", &
                              valid_range = [0.0_8, 1.0_8], &
                              fillvalue = fill_float, &
                              deflate_level = deflate_level, &
                              shuffle = shuffle, &
                              attlist=att_coord)
    endif

    call tiof_def_vars (obj, varlist, errstat)
    call tiof_varlist_free (varlist)
    call tiof_attlist_free (att_coord)

  end subroutine append_qa_vars

  !> Write all geolocation variables to Level 2 product file
  !! @param[in] first_pix  index of first xtrack pixel in use
  !! @param[in] last_pix   index of last xtrack pixel in use
  !! @param[in] ntimes     number of mirror steps in use
  !! @param[inout] errstat  Error status variable
  subroutine l2_tio_write_geo (geo, first_pix, last_pix, first_line, last_line, &
       errstat)

    USE OMSAO_variables_module, only: geo_group
    implicit none

    integer, intent(in) :: first_pix, last_pix, first_line, last_line
    integer, intent(inout) :: errstat
    TYPE(geo_group), INTENT(IN) :: geo

    integer :: num_lines, num_pix

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    num_lines = last_line - first_line + 1
    num_pix = last_pix - first_pix + 1

    ! root group
    call tiof_push_group (obj, "/", errstat)
    call tiof_put1d_i4 (obj, tempo_dim_step, [0], [num_lines], &
                        geo%step_idx(0:num_lines-1), errstat)
    call tiof_pop_group (obj, errstat)

    ! geolocation group
    call tiof_push_group (obj, o3p_grp_geolocation, errstat)
    call tiof_put1d_r8 (obj, o3p_var_time, [0], [num_lines], &
         geo%time(0:num_lines-1), errstat)
    call tiof_put2d_r4 (obj, o3p_var_latitude, [0,0], &
         [num_lines, num_pix], &
         real(geo%lat(first_pix:last_pix, 0:num_lines-1), kind=4), errstat)
    call tiof_put2d_r4 (obj, o3p_var_longitude, [0,0], &
         [num_lines, num_pix], &
         real(geo%lon(first_pix:last_pix, 0:num_lines-1), kind=4), errstat)
    call tiof_put2d_r4 (obj, o3p_var_ra_angle, [0,0], &
         [num_lines, num_pix], &
         real(geo%aza(first_pix:last_pix, 0:num_lines-1), kind=4), errstat)
    call tiof_put2d_r4 (obj, o3p_var_sz_angle, [0,0], &
         [num_lines, num_pix], &
         real(geo%sza(first_pix:last_pix, 0:num_lines-1), kind=4), errstat)
    call tiof_put2d_r4 (obj, o3p_var_vz_angle, [0,0], &
         [num_lines, num_pix], &
         real(geo%vza(first_pix:last_pix, 0:num_lines-1), kind=4), errstat)
    call tiof_put3d_r4 (obj, o3p_var_latitude_bounds, [0,0,0], &
         [num_lines, num_pix, 4], geo%clat(1:4,first_pix:last_pix, 0:num_lines-1), errstat)
    call tiof_put3d_r4 (obj, o3p_var_longitude_bounds, [0,0,0], &
         [num_lines, num_pix, 4], geo%clon(1:4,first_pix:last_pix, 0:num_lines-1), errstat)
    call tiof_pop_group (obj, errstat)

    ! support_data group
    call tiof_push_group (obj, o3p_grp_support_data, errstat)
    call tiof_put2d_i2 (obj, o3p_var_geoflg, &
         [0, 0], [num_lines, num_pix], &
         int(geo%GFlg(first_pix:last_pix, 0:num_lines-1), kind=2), &
         errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat < 0) then
      call tell_error (tell_io_write_error, "l2_tio_write_geo: failed", errstat)
      return
    endif

  end subroutine l2_tio_write_geo

  !> Write data for current pixel to L2 netCDF file
  !! @param[in]    ipix        current cross-track pixel index (0 based)
  !! @param[in]    iline       current mirror step index (0 based)
  !! @param[in]    exval       fitting exit value
  !! @param[in]    fitcol      ozone column parameters
  !! @param[in]    dfitcol     ozone column uncertainty parameters
  !! @param[in]    ngas        number of other gas parameters (nfgas)
  !! @param[in]    nlayer      number of layers in ozone profile (nlay)
  !! @param[in]    nfitvar     number of fit variables (n_fitvar_rad)
  !! @param[in]    nwindow     number of fit windows (numwin)
  !! @param[in]    num_param   number of non-gas fit parameters
  !! @param[in]    num_wav_max maximum number of wavelengths in fit
  !! @param[in]    ozfit_idxs  fit variable starting index (ozfit_start_index)
  !! @param[in]    ozfit_idxe  fit variable sending index (ozfit_start_index)
  !! @param[inout] errstat     error status variable
  subroutine l2_tio_write_data (ipix, iline, exval, fitcol, dfitcol, ngas, &
       nlayer, nfitvar, nwindow, num_param, num_wav_max, ozfit_idxs, &
       ozfit_idxe, errstat)

    use OMSAO_variables_module, only: n_rad_wvl, nradpix, fitres_rad, &
         fitweights, fitspec_rad, simspec_rad, mask_fitvar_rad, &
         fitvar_rad, fothvarpos, fitvar_rad_nstd, fitvar_rad_std, &
         fitvar_rad_str, fitvar_rad_apriori, fitvar_rad_aperror, fitwavs, &
         fitvar_rad_unit
    use OMSAO_pixelcorner_module, only: omi_Mflg
    use ozprof_data_module, only: use_lograd, the_ai, ozprof, ozprof_std, &
         ozprof_nstd, ozprof_ap, ozprof_apstd, tracegas, fgaspos, fgasidxs, &
         atmosprof, the_ctp, the_cfrac, the_cod, the_cld_flg, ozinfo, ntp, &
         eff_alb, thealbidx, glintprob, avg_kernel, covar, ncovar, contri, &
         weight_function, tropaod, tropsca, num_iter
    use ISO_C_BINDING, only: c_null_char

    implicit none

    integer, intent(in) :: ipix, iline, exval, ngas, nlayer, nfitvar, &
         nwindow, num_param, num_wav_max, ozfit_idxs, ozfit_idxe
    integer, intent(inout) ::errstat
    real (kind=8), dimension(3), intent(in) :: fitcol
    real (kind=8), dimension(3,2), intent(in) :: dfitcol

    real (kind=8), dimension (nwindow)           :: allrms, allavgres
    !real (kind=8), dimension (5, n_max_fitpars) :: tempvar
    real (kind=8), dimension (nfitvar, nfitvar) :: correl
    !real (kind=8), dimension (max_fit_pts)      :: tempring
    real (kind=8)                               :: ncorrl_foo
    real (kind=8) , dimension(n_rad_wvl) :: lfitres_rad
    real (kind=4), dimension(:), allocatable :: tmp_coord_var
    real (kind=4), dimension(:,:), allocatable :: tmp_coord_var_bnds
    integer (KIND= 2), dimension(nlayer,nlayer)  :: OzAvgK_I16
    integer (kind=2), dimension (:), allocatable :: ncorrl_1d
    character (len = 4), dimension(nfitvar)  :: varnames_nNum
    character (len = 20), dimension(nfitvar) :: units
    real (kind=8)  :: avgres
    integer (kind=4) :: i, j, fidx, lidx, ii, irow, jcol, nn
    integer (kind=4) :: num_elms, num_aeros_wavl
    logical, save :: first = .true.

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    num_elms = (nlayer * (nlayer - 1))/2
    num_aeros_wavl = nwindow + 2

    allocate (tmp_coord_var(nlayer), tmp_coord_var_bnds(2,nlayer))

    ! allocate ncorrl_1d if we need it
    if (ozwrtcovar) then
      allocate(ncorrl_1d(num_elms), stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             'l2_tio_write_data: allocate ncorrl_1d failed', errstat)
        return
      endif
    endif

    ! calculate some of the output values
    fidx = 1
    do i = 1, nwindow
      lidx = fidx + nradpix(i) - 1
      allrms(i) = &
           sqrt(sum((fitres_rad(fidx:lidx)/fitweights(fidx:lidx))**2) / &
           nradpix(i))
      fidx = lidx + 1
    enddo
      ! changes will already have been made if he5 output generated
      simspec_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) - &
           fitres_rad(1:n_rad_wvl)
      if (use_lograd) then
        fitspec_rad(1:n_rad_wvl) = exp(fitspec_rad(1:n_rad_wvl))
        simspec_rad(1:n_rad_wvl) = exp(simspec_rad(1:n_rad_wvl))
!        fitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  &
        lfitres_rad(1:n_rad_wvl) = fitspec_rad(1:n_rad_wvl) -  &
             simspec_rad(1:n_rad_wvl)
      else
        lfitres_rad(1:n_rad_wvl) = fitres_rad(1:n_rad_wvl)
      end if
!    avgres = sqrt(sum((abs(fitres_rad(1:n_rad_wvl)) / &
    avgres = sqrt(sum((abs(lfitres_rad(1:n_rad_wvl)) / &
         fitspec_rad(1:n_rad_wvl))**2.0)/n_rad_wvl)*100.0

    fidx = 1
    do i = 1, nwindow
      lidx = fidx + nradpix(i) - 1
!      allavgres(i) = sqrt(sum((abs(fitres_rad(fidx:lidx)) / &
      allavgres(i) = sqrt(sum((abs(lfitres_rad(fidx:lidx)) / &
           fitspec_rad(fidx:lidx))**2.0)/nradpix(i))*100.0
      fidx = lidx + 1
    enddo
    ! set up array of parameter names and units
    do i= 1, nfitvar
      write(varnames_nNum(i), '(A4)') &
           trim(fitvar_rad_str(mask_fitvar_rad(i)))
      write(units(i), '(A20)') &
           trim(fitvar_rad_unit(mask_fitvar_rad(i)))
    enddo

    ! Product group
    call tiof_push_group (obj, o3p_grp_product, errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_retrieve_prof, [iline, ipix, 0], &
         [1,1,nlayer], real(ozprof(1:nlayer), kind=4), errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_retrieve_prof_prec, [iline, ipix, 0], &
         [1,1,nlayer], real(ozprof_nstd(1:nlayer), kind=4), errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_retrieve_prof_err, [iline, ipix, 0], &
         [1,1,nlayer], real(ozprof_std(1:nlayer), kind=4), errstat)
    call tiof_put1d_r4 (obj, o3p_var_total_o3, [iline, ipix], &
         [1,1], [real(fitcol(1), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_total_o3_prec, [iline, ipix], &
         [1,1], [real(dfitcol(1, 2), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_total_o3_err, [iline, ipix], &
         [1,1], [real(dfitcol(1, 1), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_strat_o3, [iline, ipix], &
         [1,1], [real(fitcol(2), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_strat_o3_prec, [iline, ipix], &
         [1,1], [real(dfitcol(2, 2), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_strat_o3_err, [iline, ipix], &
         [1,1], [real(dfitcol(2, 1), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_tropo_o3, [iline, ipix], &
         [1,1], [real(fitcol(3), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_tropo_o3_prec, [iline, ipix], &
         [1,1], [real(dfitcol(3, 2), kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_tropo_o3_err, [iline, ipix], &
         [1,1], [real(dfitcol(3, 1), kind=4)], errstat)
    ! FIXME - Surface layer (0-2km) O3 column will go here in UVIS code
    ! Optional products
    ! other fitted gases
    if ( ngas > 0 .and. (gaswrt .or. ozwrtvar) ) then
      call tiof_put1d_r4 (obj, o3p_var_other_gas_retrieve, [iline, ipix, 0], &
         [1,1, ngas], real(tracegas(fgaspos(1:ngas), 4), kind=4), &
         errstat)
      call tiof_put1d_r4 (obj, o3p_var_other_gas_retrieve_prec, &
         [iline, ipix, 0], [1,1, ngas], &
         real(tracegas(fgaspos(1:ngas), 6), kind=4), errstat)
      call tiof_put1d_r4 (obj, o3p_var_other_gas_retrieve_err, &
         [iline, ipix, 0], [1,1, ngas], &
         real(tracegas(fgaspos(1:ngas), 5), kind=4), errstat)
    endif
    ! non-gas fitted params
    if (num_param > 0 .and. ozwrtvar) then
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_ret, &
         [iline, ipix, 0], [1,1,num_param], &
         real(fitvar_rad(mask_fitvar_rad(fothvarpos(1:num_param))), kind=4),&
         errstat)
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_ret_prec, &
         [iline, ipix, 0], [1,1, num_param], &
         real(fitvar_rad_nstd(mask_fitvar_rad(fothvarpos(1:num_param))), kind=4), &
         errstat)
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_ret_err, &
         [iline, ipix, 0], [1,1, num_param], &
         real(fitvar_rad_std(mask_fitvar_rad(fothvarpos(1:num_param))), kind=4), &
         errstat)
    endif

    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: writing to product group", &
                       errstat)
      return
    endif

    ! Support data group
    call tiof_push_group (obj, o3p_grp_support_data, errstat)
    call tiof_put1d_r4 (obj, o3p_var_aeros_index, [iline, ipix], [1,1], &
         [real(the_ai, kind=4)], errstat)

    tmp_coord_var(1:nlayer) = real(sqrt(atmosprof(1,0:nlayer-1) &
                                        * atmosprof(1,1:nlayer)), kind=4)
    call tiof_put1d_r4 (obj, o3p_var_profile_pres, [iline, ipix, 0], &
         [1,1, nlayer], tmp_coord_var(1:nlayer), errstat)

    tmp_coord_var(1:nlayer) = real(0.5*(atmosprof(2,0:nlayer-1) &
                                        + atmosprof(2,1:nlayer)), kind=4)
    call tiof_put1d_r4 (obj, o3p_var_profile_alt, [iline, ipix, 0], &
         [1,1, nlayer], tmp_coord_var(1:nlayer), errstat)

    tmp_coord_var(1:nlayer) = real(0.5*(atmosprof(3,0:nlayer-1) &
                                        + atmosprof(3,1:nlayer)), kind=4)
    call tiof_put1d_r4 (obj, o3p_var_profile_temp, [iline, ipix, 0], &
         [1,1, nlayer], tmp_coord_var(1:nlayer), errstat)

    tmp_coord_var_bnds(1,1:nlayer) = real(atmosprof(1,0:nlayer-1),kind=4)
    tmp_coord_var_bnds(2,1:nlayer) = real(atmosprof(1,1:nlayer  ),kind=4)
    call tiof_put2d_r4 (obj, o3p_var_profile_pres_bnds, [iline,ipix,0,0], &
         [1,1,nlayer,2], tmp_coord_var_bnds(1:2,1:nlayer), errstat)

    tmp_coord_var_bnds(1,1:nlayer) = real(atmosprof(2,0:nlayer-1),kind=4)
    tmp_coord_var_bnds(2,1:nlayer) = real(atmosprof(2,1:nlayer  ),kind=4)
    call tiof_put2d_r4 (obj, o3p_var_profile_alt_bnds, [iline,ipix,0,0], &
         [1,1,nlayer,2], tmp_coord_var_bnds(1:2,1:nlayer), errstat)

    call tiof_put1d_r4 (obj, o3p_var_o3_apriori_prof, [iline, ipix, 0], &
         [1,1,nlayer], real(ozprof_ap(1:nlayer), kind=4), errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_apriori_prof_err, [iline, ipix, 0], &
         [1,1,nlayer], real(ozprof_apstd(1:nlayer), kind=4), errstat)
    call tiof_put1d_i4 (obj, o3p_var_tropo_index, [iline, ipix], &
         [1,1], [ntp], errstat)
    call tiof_put1d_r4 (obj, o3p_var_eff_cld_frac, [iline, ipix], &
         [1,1], [real(the_cfrac, kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_eff_cld_pres, [iline, ipix], &
         [1,1], [real(the_ctp, kind=4)], errstat)
    call tiof_put1d_i4 (obj, o3p_var_cld_flag, [iline, ipix], &
         [1,1], [the_cld_flg], errstat)
    call tiof_put1d_r4 (obj, o3p_var_glint_prob, [iline, ipix], &
         [1,1], [real(glintprob, kind=4)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_surf_albedo, [iline, ipix], &
         [1,1], [real(eff_alb(thealbidx), kind=4)], errstat)
    call tiof_put1d_i4 (obj, o3p_var_fit_wavel, &
         [iline, ipix], [1,1], [n_rad_wvl], errstat)
    call tiof_put1d_i4 (obj, o3p_var_window_wavel, &
         [iline, ipix, 0], [1,1,nwindow], nradpix(1:nwindow), errstat)
    ! Optional parameters
    ! Other fitted gases
    if ( ngas > 0 .and. (gaswrt .or. ozwrtvar) ) then
      call tiof_put1d_r4 (obj, o3p_var_other_gas_apriori, &
         [iline, ipix, 0], [1,1, ngas], &
         real(tracegas(fgaspos(1:ngas), 2), kind=4), errstat)
      call tiof_put1d_r4 (obj, o3p_var_other_gas_apriori_err, &
         [iline, ipix, 0], [1,1, ngas], &
         real(tracegas(fgaspos(1:ngas), 3), kind=4), errstat)
      if (first) then
! FIXME
! see non-gas param names below
!        call tiof_put1d_string (obj, o3p_var_other_gas_names, 0, ngas, &
!             varnames_nNum(fgasidxs(fgaspos(1:ngas))), errstat)
        do i=1,ngas
          call tiof_put1d_string (obj, o3p_var_other_gas_names, i-1, 1, &
               [varnames_nNum(fgasidxs(fgaspos(i)))//c_null_char], errstat)
        enddo
      endif
    endif
    ! Non-gas fitted parameters
    if (num_param > 0 .and. ozwrtvar) then
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_apriori, &
         [iline, ipix, 0], [1,1,num_param], &
         real(fitvar_rad_apriori(mask_fitvar_rad(fothvarpos(1:num_param))), &
         kind=4), errstat)
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_apriori_err, &
         [iline, ipix, 0], [1,1,num_param], &
         real(fitvar_rad_aperror(fothvarpos(1:num_param)), &
         kind=4), errstat)
      if (first) then
        do i = 1, num_param
          call tiof_put1d_string (obj, o3p_var_nongas_param_names, i-1, &
               1, [varnames_nNum(fothvarpos(i))//c_null_char], errstat)
          call tiof_put1d_string (obj, o3p_var_nongas_param_units, i-1, &
               1, [units(fothvarpos(i))//c_null_char], errstat)
        enddo
      endif
    endif
    !averaging kernels
    if (ozwrtavgk) then
      i = ozfit_idxs; j = ozfit_idxe
      ! Transpose averaging kernels, so in he5, row x col (same as avg_kernel)
      OzAvgK_I16(1:nlayer,1:nlayer) = nint(avg_kernel(i:j, i:j)*1.0d4, KIND= 2)
      call tiof_put2d_i2 (obj, o3p_var_o3_avg_kernel, [iline, ipix, 0, 0], &
           [1,1,nlayer,nlayer], transpose(OzAvgK_I16(1:nlayer,1:nlayer)), errstat)
    endif
    !correlation matrix
    if (ozwrtcorr) then
      do i = 1, nfitvar
        correl(i, i) = 1.0
        do j = 1, i - 1
          correl(i, j) = covar(i, j) / sqrt(covar(i, i) * covar(j, j))
          correl(j, i) = correl(i, j)
        enddo
      enddo
      call tiof_put2d_r4 (obj, o3p_var_correl, &
         [iline, ipix, 0, 0], [1,1,nfitvar,nfitvar], &
         real(correl(1:nfitvar, 1:nfitvar), &
         kind=4), errstat)
    endif
    !noise matrix
    if (ozwrtcovar) then
      i = ozfit_idxs; j = ozfit_idxe
      nn = j-i+1
      if( nn /= nlayer ) then
        call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: nn not equal to nlayer", &
                       errstat)
        return
      endif
      ii = 0
      do irow=1,nn-1
        do jcol = irow+1, nn
          ii = ii + 1
          ncorrl_foo    = ncovar(irow,jcol) &
               / sqrt( ncovar(irow,irow)*ncovar(jcol,jcol))
          ncorrl_1d(ii) = nint( ncorrl_foo*1.0d04 , kind=2)
        enddo
      enddo
      if( ii /= num_elms) then
        call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: upper off-diagonal elements size mismatch", &
                       errstat)
        return
      endif
      call tiof_put1d_i2 (obj, o3p_var_o3_noise_matrix, &
         [iline, ipix, 0], [1,1,num_elms], &
         ncorrl_1d(1:num_elms), errstat)
      !
      ! FIXME - ozone information content should have its own switch
      call tiof_put1d_r4 (obj, o3p_var_o3_info_content, [iline, ipix], &
         [1,1], [real(ozinfo, kind=4)], errstat)
    endif
    !contribution matrix
    if (ozwrtcontri) then
      if (num_wav_max > n_rad_wvl) then
        contri(1:nfitvar, n_rad_wvl+1:num_wav_max) = 0.D0
      endif
      call tiof_put2d_r4 (obj, o3p_var_contrib_func, &
         [iline, ipix, 0, 0], [1,1,nfitvar,num_wav_max], &
         real(transpose(contri(1:nfitvar, 1:num_wav_max)), KIND=4), &
         errstat)
    endif
   !cloud optical depth
    if (.not. do_lambcld) then
      call tiof_put1d_r4 (obj, o3p_var_cld_opt_depth, &
         [iline, ipix], [1,1], &
         [real(the_cod, kind=4)], errstat)
    endif
    !aerosols
    if (aerosol) then
      call tiof_put1d_r4 (obj, o3p_var_aeros_opt_thick, &
         [iline, ipix, 0], [1,1,num_aeros_wavl], &
         real(tropaod(1:num_aeros_wavl), kind=4), errstat)
      call tiof_put1d_r4 (obj, o3p_var_aeros_scatter_thick, &
         [iline, ipix, 0], [1,1,num_aeros_wavl], &
         real(tropsca(1:num_aeros_wavl), kind=4), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: writing to support data group", &
                       errstat)
      return
    endif

    ! Diagnostic group
    if (have_diagnostic_group) then
      call tiof_push_group (obj, o3p_grp_diagnostic, errstat)
    endif
    ! Optional parameters
    ! fitted wavelengths
    if ( (.not. reduce_resolution) .and. &
         (ozwrtcontri .or. ozwrtres .or. ozwrtwf .or. ozwrtsnr &
         .or. wrtring)) then
      if (num_wav_max > n_rad_wvl) then
        fitwavs(n_rad_wvl + 1 : num_wav_max) = 0.D0
      endif
      call tiof_put1d_r4 (obj, o3p_var_wavel, &
           [iline, ipix, 0], [1,1,num_wav_max], &
           real(fitwavs(1:num_wav_max), kind=4), errstat)
    endif
    !weighting function
    if (ozwrtwf) then
      if (num_wav_max > n_rad_wvl) then
        weight_function(n_rad_wvl+1:num_wav_max, 1:nfitvar) = 0.D0
      endif
      call tiof_put2d_r4 (obj, o3p_var_weight_func, &
         [iline, ipix, 0, 0], [1,1,num_wav_max,nfitvar], &
         real(transpose(weight_function(1:num_wav_max,1:nfitvar)), &
         kind=4), errstat)
    endif
     ! normalized observed & simulated spectra
    if (ozwrtres) then
      if (num_wav_max > n_rad_wvl) then
        fitspec_rad(n_rad_wvl + 1 : num_wav_max) = 0.D0
        simspec_rad(n_rad_wvl + 1 : num_wav_max) = 0.D0
      endif

      call tiof_put1d_r4 (obj, o3p_var_sim_norm_rad, &
           [iline, ipix, 0], [1,1,num_wav_max], &
           real(simspec_rad(1:num_wav_max), kind=4), errstat)
      call tiof_put1d_r4 (obj, o3p_var_norm_radiance, &
           [iline, ipix, 0], [1,1,num_wav_max], &
           real(fitspec_rad(1:num_wav_max), kind=4), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: writing to diagnostic group", &
                       errstat)
      return
    endif

    ! QA stats group
    call tiof_push_group (obj, o3p_grp_qa_stats, errstat)
    call tiof_put1d_i4 (obj, o3p_var_exit_stat, &
         [iline, ipix], [1,1], [exval], errstat)
    call tiof_put1d_i4 (obj, o3p_var_iter, &
         [iline, ipix], [1,1], [num_iter], errstat)
    call tiof_put1d_r4 (obj, o3p_var_fit_rms, &
         [iline, ipix, 0], [1,1,nwindow], &
         real(allrms(1:nwindow), kind=4), errstat)
    call tiof_put1d_r4 (obj, o3p_var_avg_resid, &
         [iline, ipix, 0], [1,1,nwindow], &
         real(allavgres(1:nwindow), kind=4), errstat)
    !call tiof_put1d_i2 (obj, o3p_var_mqf, &
    !     [iline], [1], [omi_Mflg(iline)], errstat)
    !Optional param - relative measurement error
    if (ozwrtsnr) then
      call tiof_put1d_r4 (obj, o3p_var_merr, &
           [iline, ipix, 0], [1,1,num_wav_max], &
           real(fitweights(1:num_wav_max), kind=4), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: writing to qa_stats group", &
                       errstat)
      return
    endif

    !tidy up
    first = .false.
    if (allocated(ncorrl_1d)) deallocate(ncorrl_1d, stat=errstat)

  end subroutine l2_tio_write_data

  !> Write data for current pixel to L2 netCDF file
  !! @param[in]    ipix        current cross-track pixel index (0 based)
  !! @param[in]    iline       current mirror step index (0 based)
  !! @param[in]    exval       fitting exit value
  !! @param[in]    ngas        number of other gas parameters (nfgas)
  !! @param[in]    nlayer      number of layers in ozone profile (nfgas)
  !! @param[in]    nfitvar     number of fit variables (n_fitvar_rad)
  !! @param[in]    nwindow     number of fit windows (numwin)
  !! @param[in]    num_param   number of non-gas fit parameters
  !! @param[in]    num_wav_max maximum number of wavelengths in fit
  !! @param[in]    ozfit_idxs  fit variable starting index (ozfit_start_index)
  !! @param[in]    ozfit_idxe  fit variable sending index (ozfit_start_index)
  !! @param[inout] errstat     error status variable
  subroutine l2_tio_fill_data (ipix, iline, exval, ngas, nlayer, nfitvar, &
       nwindow, num_param, num_wav_max, ozfit_idxs, ozfit_idxe, errstat)

    use OMSAO_pixelcorner_module, only: omi_Mflg

    implicit none

    integer, intent(in) :: ipix, iline, exval, ngas, nlayer, nfitvar, &
         nwindow, ozfit_idxs, ozfit_idxe
    integer, intent(inout) ::errstat

    integer (kind=4) :: i, j
    integer (kind=4) :: num_elms, num_aeros_wavl, num_param, &
         num_wav_max
    real (KIND=4), dimension(nlayer)      :: tmp1D_layer
    real (KIND=4), dimension(2,nlayer) :: tmp1D_layer_bnds
    real (KIND=4), dimension(nfitvar) :: tmp1D_fitvar
    real (KIND=4), dimension(num_wav_max)   :: tmp1D_fitpts
    real (KIND=4), dimension(nwindow)        :: tmp1D_numwin
    !integer (KIND=2), dimension(maxwin)      :: tmp1D_num
    integer (KIND=4), dimension(nwindow)      :: tmp1D_num
    integer (KIND=2), dimension(nlayer*(nlayer-1)/2)       :: tmp1D_ncorrl
    integer(KIND=2), dimension(nfitvar, nfitvar) :: tmp2D_fitvarK16
    real (KIND=4), dimension(nwindow+2)      :: tmp1D_aer
    real (KIND=4), dimension(nfitvar, nfitvar) :: tmp2D_fitvar
    real (KIND=4), dimension(nfitvar, num_wav_max)   :: tmp2D_contri
    real (KIND=4), dimension(num_wav_max, nfitvar)   :: tmp2D_wf

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    num_elms = (nlayer * (nlayer - 1))/2
    num_aeros_wavl = nwindow + 2

    tmp1D_layer_bnds(1:2,1:nlayer)  = real(fill_double, kind=4)
    tmp1D_layer(1:nlayer)       = real(fill_double, kind=4)
    tmp1D_numwin(1:nwindow)     = real(fill_double, kind=4)
    tmp1D_num(1:nwindow)        = int(fill_int16, kind=4)
    tmp1D_fitvar(1:nfitvar)     = real(fill_double, kind=4)
    tmp1D_fitpts(1:num_wav_max) = real(fill_double, kind=4)

    ! Product group
    call tiof_push_group (obj, o3p_grp_product, errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_retrieve_prof, [iline, ipix, 0], &
         [1,1,nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_retrieve_prof_prec, [iline, ipix, 0], &
         [1,1,nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_retrieve_prof_err, [iline, ipix, 0], &
         [1,1,nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_total_o3, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_total_o3_prec, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_total_o3_err, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_strat_o3, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_strat_o3_prec, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_strat_o3_err, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_tropo_o3, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_tropo_o3_prec, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_tropo_o3_err, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    ! FIXME - Surface layer (0-2km) O3 column will go here in UVIS code
    ! Optional products
    ! other fitted gases
    if ( ngas > 0 .and. (gaswrt .or. ozwrtvar) ) then
      call tiof_put1d_r4 (obj, o3p_var_other_gas_retrieve, [iline, ipix, 0], &
         [1,1, ngas], tmp1D_fitvar(1:ngas), errstat)
      call tiof_put1d_r4 (obj, o3p_var_other_gas_retrieve_prec, &
         [iline, ipix, 0], [1,1, ngas], &
         tmp1D_fitvar(1:ngas), errstat)
      call tiof_put1d_r4 (obj, o3p_var_other_gas_retrieve_err, &
         [iline, ipix, 0], [1,1, ngas], &
         tmp1D_fitvar(1:ngas), errstat)
     endif
    ! non-gas fitted params
    if (num_param > 0 .and. ozwrtvar) then
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_ret, &
         [iline, ipix, 0], [1,1,num_param], &
         tmp1D_fitvar(1:num_param), errstat)
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_ret_prec, &
         [iline, ipix, 0], [1,1, num_param], &
         tmp1D_fitvar(1:num_param), errstat)
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_ret_err, &
         [iline, ipix, 0], [1,1, num_param], &
         tmp1D_fitvar(1:num_param), errstat)
    endif

    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_fill_data: filling in product group", &
                       errstat)
      return
    endif

    ! Support data group
    call tiof_push_group (obj, o3p_grp_support_data, errstat)
    call tiof_put1d_r4 (obj, o3p_var_aeros_index, [iline, ipix], [1,1], &
         [tmp1D_layer(1)], errstat)
    call tiof_put2d_r4 (obj, o3p_var_profile_pres_bnds, [iline, ipix, 0,0], &
         [1,1, nlayer,2], tmp1D_layer_bnds(1:2,1:nlayer), errstat)
    call tiof_put2d_r4 (obj, o3p_var_profile_alt_bnds, [iline, ipix, 0], &
         [1,1, nlayer,2], tmp1D_layer_bnds(1:2,1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_profile_pres, [iline, ipix, 0], &
         [1,1, nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_profile_alt, [iline, ipix, 0], &
         [1,1, nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_profile_temp, [iline, ipix, 0], &
         [1,1, nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_apriori_prof, [iline, ipix, 0], &
         [1,1,nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_r4 (obj, o3p_var_o3_apriori_prof_err, [iline, ipix, 0], &
         [1,1,nlayer], tmp1D_layer(1:nlayer), errstat)
    call tiof_put1d_i4(obj, o3p_var_tropo_index, [iline, ipix], &
         [1,1], [int(fill_int16)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_eff_cld_frac, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_eff_cld_pres, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_i4 (obj, o3p_var_cld_flag, [iline, ipix], &
         [1,1], [INT(fill_int16)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_glint_prob, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_surf_albedo, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    call tiof_put1d_i4 (obj, o3p_var_fit_wavel, &
         [iline, ipix], [1,1], [int(fill_int16)], errstat)
    call tiof_put1d_i4 (obj, o3p_var_window_wavel, &
         [iline, ipix, 0], [1,1,nwindow], tmp1D_num(1:nwindow), errstat)
    ! Optional parameters
    ! Other fitted gases
    if ( ngas > 0 .and. (gaswrt .or. ozwrtvar) ) then
      call tiof_put1d_r4 (obj, o3p_var_other_gas_apriori, &
         [iline, ipix, 0], [1,1, ngas], tmp1D_fitvar(1:ngas), errstat)
      call tiof_put1d_r4 (obj, o3p_var_other_gas_apriori_err, &
         [iline, ipix, 0], [1,1, ngas], tmp1D_fitvar(1:ngas), errstat)
    endif
    ! Non-gas fitted parameters
    if (num_param > 0 .and. ozwrtvar) then
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_apriori, &
         [iline, ipix, 0], [1,1,num_param], tmp1D_fitvar(1:num_param), &
         errstat)
      call tiof_put1d_r4 (obj, o3p_var_nongas_param_apriori_err, &
         [iline, ipix, 0], [1,1,num_param], tmp1D_fitvar(1:num_param), &
         errstat)
    endif
    !averaging kernels
    if (ozwrtavgk) then
      i = ozfit_idxs; j = ozfit_idxe
      tmp2D_fitvarK16(i:j, i:j) = int(fill_int16, kind=2)
      call tiof_put2d_i2 (obj, o3p_var_o3_avg_kernel, [iline, ipix, 0, 0], &
           [1,1,nlayer,nlayer], tmp2D_fitvarK16(i:j, i:j), errstat)
    endif
    !correlation matrix
    if (ozwrtcorr) then
      call tiof_put2d_r4 (obj, o3p_var_correl, &
         [iline, ipix, 0, 0], [1,1,nfitvar,nfitvar], &
         tmp2D_fitvar(1:nfitvar, 1:nfitvar), errstat)
    endif
    !noise matrix
    if (ozwrtcovar) then
      i = ozfit_idxs; j = ozfit_idxe
      tmp1D_ncorrl(1:num_elms) = int(fill_int16, kind=2)
      call tiof_put1d_i2 (obj, o3p_var_o3_noise_matrix, &
         [iline, ipix, 0], [1,1,num_elms], &
         tmp1D_ncorrl(1:num_elms), errstat)
      !
      ! FIXME - ozone info content should have its own switch eventually
      call tiof_put1d_r4 (obj, o3p_var_o3_info_content, [iline, ipix], &
         [1,1], [tmp1D_layer(1)], errstat)
    endif
    !contribution matrix
    if (ozwrtcontri) then
      tmp2D_contri(1:nfitvar, 1:num_wav_max) = real(fill_double, kind=4)
      call tiof_put2d_r4 (obj, o3p_var_contrib_func, &
         [iline, ipix, 0, 0], [1,1,nfitvar,num_wav_max], &
         tmp2D_contri(1:nfitvar, 1:num_wav_max), errstat)
    endif
    !cloud optical depth
    if (.not. do_lambcld) then
      call tiof_put1d_r4 (obj, o3p_var_cld_opt_depth, &
         [iline, ipix], [1,1], &
         [tmp1D_layer(1)], errstat)
    endif
    !aerosols
    if (aerosol) then
      call tiof_put1d_r4 (obj, o3p_var_aeros_opt_thick, &
         [iline, ipix, 0], [1,1,num_aeros_wavl], &
         tmp1D_aer(1:num_aeros_wavl), errstat)
      call tiof_put1d_r4 (obj, o3p_var_aeros_scatter_thick, &
         [iline, ipix, 0], [1,1,num_aeros_wavl], &
         tmp1D_aer(1:num_aeros_wavl), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_fill_data: filling in support data group", &
                       errstat)
      return
    endif

    ! Diagnostic group
    if (have_diagnostic_group) then
      call tiof_push_group (obj, o3p_grp_diagnostic, errstat)
    endif
    ! Optional parameters
    ! fitted wavelengths
    if ( (.not. reduce_resolution) .and. &
         (ozwrtcontri .or. ozwrtres .or. ozwrtwf .or. ozwrtsnr &
         .or. wrtring)) then
      call tiof_put1d_r4 (obj, o3p_var_wavel, [iline, ipix, 0], &
           [1,1,num_wav_max], tmp1D_fitpts(1:num_wav_max), errstat)
    endif
    !weighting function
    if (ozwrtwf) then
      tmp2D_wf(1:num_wav_max, 1:nfitvar) = real(fill_double, kind=4)
      call tiof_put2d_r4 (obj, o3p_var_weight_func, &
         [iline, ipix, 0, 0], [1,1,num_wav_max,nfitvar], &
         tmp2D_wf(1:num_wav_max,1:nfitvar), errstat)
    endif
    ! normalized observed & simulated spectra
    if (ozwrtres) then
      call tiof_put1d_r4 (obj, o3p_var_sim_norm_rad, [iline, ipix, 0], &
           [1,1,num_wav_max], tmp1D_fitpts(1:num_wav_max), errstat)
      call tiof_put1d_r4 (obj, o3p_var_norm_radiance, [iline, ipix, 0], &
           [1,1,num_wav_max], tmp1D_fitpts(1:num_wav_max), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_fill_data: filling in diagnostic group", &
                       errstat)
      return
    endif

    ! QA stats group
    call tiof_push_group (obj, o3p_grp_qa_stats, errstat)
    call tiof_put1d_i4 (obj, o3p_var_exit_stat, &
         [iline, ipix], [1,1], [exval], errstat) ! always write exit value
    call tiof_put1d_i4 (obj, o3p_var_iter, &
         [iline, ipix], [1,1], [int(fill_uint1)], errstat)
    call tiof_put1d_r4 (obj, o3p_var_fit_rms, &
         [iline, ipix, 0], [1,1,nwindow], &
         tmp1D_numwin(1:nwindow), errstat)
    call tiof_put1d_r4 (obj, o3p_var_avg_resid, &
         [iline, ipix, 0], [1,1,nwindow], &
         tmp1D_numwin(1:nwindow), errstat)
    !call tiof_put1d_i2 (obj, o3p_var_mqf, & ! there should always be an mflag
    !     [iline], [1], [omi_Mflg(iline)], errstat)
    !Optional param - relative measurement error
    if (ozwrtsnr) then
      call tiof_put1d_r4 (obj, o3p_var_merr, [iline, ipix, 0], &
           [1,1,num_wav_max], tmp1D_fitpts(1:num_wav_max), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_fill_data: filling in qa_stats group", &
                       errstat)
      return
    endif

  end subroutine l2_tio_fill_data

  !> Write all geolocation variables to merged Level 2 product file
  !-----------------------------------------------------------------
  !
  !> @param[in]    min_step   mirror step index at which to start writing
  !> @param[in]    max_step   mirror step index at which to stop writing
  !> @param[in]    min_xtrack cross-track index at which to start writing
  !> @param[in]    max_xtrack cross-track index at which to stop writing
  !> @param[in]    ncorner    size of corner dimension
  !> @param[inout] errstat    Error status variable
  !
  !> @author E. O'Sullivan November 2016
  !-----------------------------------------------------------------
  subroutine write_merged_geo (min_step, max_step, min_xtrack, max_xtrack, ncorner, errstat)

    use m_o3p_params, only: lon, lat, aza, sza, vza, corner_lon, corner_lat, &
         time, geoflg

    implicit none

    integer (kind=4), intent(in) :: min_step, max_step, min_xtrack, &
         max_xtrack, ncorner
    integer, intent(inout) :: errstat

    ! local variables
    integer (kind=4) :: nstep, nxtrack

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    nstep=max_step-min_step+1
    nxtrack=max_xtrack-min_xtrack+1

    ! geolocation group
    call tiof_push_group (obj, o3p_grp_geolocation, errstat)
    call tiof_put1d_r8 (obj, o3p_var_time, [min_step], [nstep], time, errstat)
    call tiof_put2d_r4 (obj, o3p_var_latitude, [min_step, min_xtrack], &
         [nstep, nxtrack], lat(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_longitude, [min_step, min_xtrack], &
         [nstep, nxtrack], lon(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_ra_angle, [min_step, min_xtrack], &
         [nstep, nxtrack], aza(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_sz_angle, [min_step, min_xtrack], &
         [nstep, nxtrack], sza(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_vz_angle, [min_step, min_xtrack], &
         [nstep, nxtrack], vza(1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_latitude_bounds, &
         [min_step, min_xtrack, 0], &
         [nstep, nxtrack, ncorner], corner_lat(1:ncorner,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_longitude_bounds, &
         [min_step, min_xtrack, 0], &
         [nstep, nxtrack, ncorner], corner_lon(1:ncorner,1:nxtrack,1:nstep), errstat)
    call tiof_pop_group (obj, errstat)

    ! support_data group
    call tiof_push_group (obj, o3p_grp_support_data, errstat)
    call tiof_put2d_i2 (obj, o3p_var_geoflg, [min_step, min_xtrack], &
         [nstep, nxtrack], geoflg(1:nxtrack,1:nstep), errstat)
    call tiof_pop_group (obj, errstat)

    if (errstat < 0) then
      call tell_error (tell_io_write_error, "write_merged_geo: failed", &
           errstat)
      return
    endif

  end subroutine write_merged_geo

  !> Write all geolocation variables to merged Level 2 product file
  !-----------------------------------------------------------------
  !
  !> @param[in]    min_step     mirror step index at which to start writing
  !> @param[in]    max_step     mirror step index at which to stop writing
  !> @param[in]    min_xtrack   cross-track index at which to start writing
  !> @param[in]    max_xtrack   cross-track index at which to stop writing
  !> @param[in]    ngas         size of gas variables dimension
  !> @param[in]    nnongas      size of nongas variables dimension
  !> @param[in]    nlayer       size of layer dimension
  !> @param[in]    nfitvars     size of fit variables dimension
  !> @param[in]    nfitwins     size of fit windows dimension
  !> @param[in]    nmax_wavs    size of wavelength dimension
  !> @param[in]    nnoise_elems size of noise elements dimension
  !> @param[in]    naeros_wavs  size of aerosol wavelengths dimension
  !> @param[inout] errstat      Error status variable
  !
  !> @author E. O'Sullivan November 2016
  !-----------------------------------------------------------------
  subroutine write_merged_data (min_step, max_step, min_xtrack, max_xtrack, &
       ngas, nnongas, nlayer, nfitvars, nfitwins, nmax_wavs, &
       nnoise_elems, naeros_wavs, errstat)

    use m_o3p_params, only: ozprof, ozprof_prec, ozprof_err, o3tot, &
         o3tot_prec, o3tot_err, o3strat, o3strat_prec, o3strat_err, &
         o3trop, o3trop_prec, o3trop_err, gas, gas_prec, gas_err, &
         nongas, nongas_prec, nongas_err, aeros_idx, ozprof_pres, &
         ozprof_alt, ozprof_temp, ozprof_pres_bnds, ozprof_alt_bnds, &
         o3apriori, o3apriori_err, tropo_idx, &
         cld_frac, cld_pres, cld_flag, glintprob, eff_alb, n_fit_wvl, &
         n_window_wvl, gas_apriori, gas_apriori_err, nongas_apriori, &
         nongas_apriori_err, gas_names, nongas_names, nongas_units, &
         avg_kernel, correl_mtrx, noise_mtrx, ozinfo, contrib_mtrx, &
         cld_opt_depth, aeros_opt_thick, aeros_scatter_thick, &
         exval,wavelengths, wgt_func, norm_rad, sim_norm_rad, exval, &
         iterations, mqf, rms, avg_resid, fit_wgt
    use ozprof_data_module, only: ozwrtavgk, ozwrtcorr, ozwrtcovar, &
         ozwrtcontri, ozwrtres, ozwrtwf, ozwrtsnr, &
         ozwrtvar, gaswrt, aerosol, do_lambcld
    use OMSAO_variables_module, only: reduce_resolution
    use ISO_C_BINDING, only: c_null_char

    implicit none

    integer, intent(in) :: min_step, max_step, min_xtrack, max_xtrack, &
         ngas, nnongas, nlayer, nfitvars, &
         nfitwins, nmax_wavs, nnoise_elems, naeros_wavs
    integer, intent(inout) ::errstat
    integer :: nstep, nxtrack

    type (tiof_file_type), pointer :: obj

    if (errstat < 0) return

    obj => primary_output_file

    nstep=max_step-min_step+1
    nxtrack=max_xtrack-min_xtrack+1

    ! Product group
    call tiof_push_group (obj, o3p_grp_product, errstat)
    call tiof_put3d_r4 (obj, o3p_var_o3_retrieve_prof, &
         [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nlayer], ozprof(1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_o3_retrieve_prof_prec, &
         [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nlayer], ozprof_prec(1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_o3_retrieve_prof_err, &
         [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nlayer], ozprof_err(1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_total_o3, [min_step, min_xtrack], &
         [nstep, nxtrack], o3tot(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_total_o3_prec, [min_step, min_xtrack], &
         [nstep, nxtrack], o3tot_prec(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_total_o3_err, [min_step, min_xtrack], &
         [nstep, nxtrack], o3tot_err(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_strat_o3, [min_step, min_xtrack], &
         [nstep, nxtrack], o3strat(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_strat_o3_prec, [min_step, min_xtrack], &
         [nstep, nxtrack], o3strat_prec(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_strat_o3_err, [min_step, min_xtrack], &
         [nstep, nxtrack], o3strat_err(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_tropo_o3, [min_step, min_xtrack], &
         [nstep, nxtrack], o3trop(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_tropo_o3_prec, [min_step, min_xtrack], &
         [nstep, nxtrack], o3trop_prec(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_tropo_o3_err, [min_step, min_xtrack], &
         [nstep, nxtrack], o3trop_err(1:nxtrack,1:nstep), errstat)
    ! FIXME - Surface layer (0-2km) O3 column will go here in UVIS code
    ! Optional products
    ! other fitted gases
    if ( ngas > 0 .and. gaswrt ) then
      call tiof_put3d_r4 (obj, o3p_var_other_gas_retrieve, &
           [min_step, min_xtrack, 0], [nstep, nxtrack, ngas], &
                          gas(1:ngas,1:nxtrack,1:nstep), errstat)
      call tiof_put3d_r4 (obj, o3p_var_other_gas_retrieve_prec, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, ngas], &
                          gas_prec(1:ngas,1:nxtrack,1:nstep), errstat)
      call tiof_put3d_r4 (obj, o3p_var_other_gas_retrieve_err, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, ngas], &
                          gas_err(1:ngas,1:nxtrack,1:nstep), errstat)
    endif
    ! non-gas fitted params
    if (nnongas > 0 .and. ozwrtvar) then
      call tiof_put3d_r4 (obj, o3p_var_nongas_param_ret, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nnongas], &
                          nongas(1:nnongas,1:nxtrack,1:nstep), errstat)
      call tiof_put3d_r4 (obj, o3p_var_nongas_param_ret_prec, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nnongas], &
                          nongas_prec(1:nnongas,1:nxtrack,1:nstep), &
         errstat)
      call tiof_put3d_r4 (obj, o3p_var_nongas_param_ret_err, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nnongas], &
                          nongas_err(1:nnongas,1:nxtrack,1:nstep), &
         errstat)
    endif

    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "write_merged_data: writing to product group", &
                       errstat)
      return
    endif

    ! Support data group
    call tiof_push_group (obj, o3p_grp_support_data, errstat)
    call tiof_put2d_r4 (obj, o3p_var_aeros_index, [min_step, min_xtrack], &
         [nstep, nxtrack], aeros_idx(1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_profile_pres, [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nlayer], ozprof_pres(1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_profile_alt, [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nlayer], ozprof_alt(1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_profile_temp, [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nlayer], ozprof_temp(1:nlayer,1:nxtrack,1:nstep), errstat)

    call tiof_put4d_r4 (obj, o3p_var_profile_pres_bnds, [min_step, min_xtrack,0,0], &
         [nstep, nxtrack, nlayer,2], ozprof_pres_bnds(1:2,1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put4d_r4 (obj, o3p_var_profile_alt_bnds, [min_step, min_xtrack,0,0], &
         [nstep, nxtrack, nlayer,2], ozprof_alt_bnds(1:2,1:nlayer,1:nxtrack,1:nstep), errstat)

    call tiof_put3d_r4 (obj, o3p_var_o3_apriori_prof, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nlayer], &
                        o3apriori(1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_o3_apriori_prof_err, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nlayer], &
                        o3apriori_err(1:nlayer,1:nxtrack,1:nstep), errstat)
    call tiof_put2d_i4 (obj, o3p_var_tropo_index, [min_step, min_xtrack], &
         [nstep, nxtrack], tropo_idx(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_eff_cld_frac, [min_step, min_xtrack], &
         [nstep, nxtrack], cld_frac(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_eff_cld_pres, [min_step, min_xtrack], &
         [nstep, nxtrack], cld_pres(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_i4 (obj, o3p_var_cld_flag, [min_step, min_xtrack], &
         [nstep, nxtrack], cld_flag(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_glint_prob, [min_step, min_xtrack], &
         [nstep, nxtrack], glintprob(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_r4 (obj, o3p_var_surf_albedo, [min_step, min_xtrack], &
         [nstep, nxtrack], eff_alb(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_i4 (obj, o3p_var_fit_wavel, &
         [min_step, min_xtrack], [nstep, nxtrack], n_fit_wvl(1:nxtrack,1:nstep), errstat)
    call tiof_put3d_i4 (obj, o3p_var_window_wavel, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nfitwins], &
                        n_window_wvl(1:nfitwins,1:nxtrack,1:nstep), errstat)
    ! Optional parameters
    ! Other fitted gases
    if ( ngas > 0 .and. gaswrt ) then
      call tiof_put3d_r4 (obj, o3p_var_other_gas_apriori, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, ngas], &
                          gas_apriori(1:ngas,1:nxtrack,1:nstep), errstat)
      call tiof_put3d_r4 (obj, o3p_var_other_gas_apriori_err, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, ngas], &
                          gas_apriori_err(1:ngas,1:nxtrack,1:nstep), errstat)
      call tiof_put1d_string (obj, o3p_var_other_gas_names, 0, ngas, &
             gas_names, errstat)
!        do i=1,ngas
!          call tiof_put1d_string (obj, o3p_var_other_gas_names, i-1, 1, &
!               [varnames_nNum(fgasidxs(fgaspos(i)))//c_null_char], errstat)
!        enddo
    endif
    ! Non-gas fitted parameters
    if (nnongas > 0 .and. ozwrtvar) then
      call tiof_put3d_r4 (obj, o3p_var_nongas_param_apriori, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nnongas], &
         nongas_apriori(1:nnongas,1:nxtrack,1:nstep), errstat)
      call tiof_put3d_r4 (obj, o3p_var_nongas_param_apriori_err, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nnongas], &
         nongas_apriori_err(1:nnongas,1:nxtrack,1:nstep), errstat)
      call tiof_put1d_string (obj, o3p_var_nongas_param_names, 0, &
           nnongas, nongas_names, errstat)
      call tiof_put1d_string (obj, o3p_var_nongas_param_units, 0, &
           nnongas, nongas_units, errstat)
!        do i = 1, nnongas
!          call tiof_put1d_string (obj, o3p_var_nongas_param_names, i-1, &
!               1, [varnames_nNum(fothvarpos(i))//c_null_char], errstat)
!          call tiof_put1d_string (obj, o3p_var_nongas_param_units, i-1, &
!               1, [units(fothvarpos(i))//c_null_char], errstat)
!        enddo
    endif
    !averaging kernels
    if (ozwrtavgk) then
      call tiof_put4d_i2 (obj, o3p_var_o3_avg_kernel, &
           [min_step, min_xtrack, 0, 0], &
           [nstep, nxtrack, nlayer, nlayer], &
                          avg_kernel(1:nlayer,1:nlayer,1:nxtrack,1:nstep), errstat)
    endif
    !correlation matrix
    if (ozwrtcorr) then
      call tiof_put4d_r4 (obj, o3p_var_correl, [min_step, min_xtrack, 0, 0], &
           [nstep, nxtrack, nfitvars, nfitvars], &
                          correl_mtrx(1:nfitvars,1:nfitvars,1:nxtrack,1:nstep), errstat)
    endif
    !noise matrix
    if (ozwrtcovar) then
      call tiof_put3d_i2 (obj, o3p_var_o3_noise_matrix, &
         [min_step, min_xtrack, 0], [nstep, nxtrack, nnoise_elems], &
         noise_mtrx(1:nnoise_elems,1:nxtrack,1:nstep), errstat)
      !
      ! FIXME - ozone information content should have its own switch
      call tiof_put2d_r4 (obj, o3p_var_o3_info_content, &
           [min_step, min_xtrack], [nstep, nxtrack], &
                          ozinfo(1:nxtrack,1:nstep), errstat)
    endif
    !contribution matrix
    if (ozwrtcontri) then
      call tiof_put4d_r4 (obj, o3p_var_contrib_func, &
         [min_step, min_xtrack, 0, 0], [nstep, nxtrack, nfitvars, nmax_wavs], &
         contrib_mtrx(1:nmax_wavs,1:nfitvars,1:nxtrack,1:nstep), errstat)
    endif
   !cloud optical depth
    if (.not. do_lambcld) then
      call tiof_put2d_r4 (obj, o3p_var_cld_opt_depth, [min_step, min_xtrack], &
           [nstep, nxtrack], cld_opt_depth(1:nxtrack,1:nstep), errstat)
    endif
    !aerosols
    if (aerosol) then
      call tiof_put3d_r4 (obj, o3p_var_aeros_opt_thick, &
           [min_step, min_xtrack, 0], [nstep, nxtrack, naeros_wavs], &
           aeros_opt_thick(1:naeros_wavs,1:nxtrack,1:nstep), errstat)
      call tiof_put3d_r4 (obj, o3p_var_aeros_scatter_thick, &
           [min_step, min_xtrack, 0], [nstep, nxtrack, naeros_wavs], &
           aeros_scatter_thick(1:naeros_wavs,1:nxtrack,1:nstep), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "write_merged_data: writing to support data group", &
                       errstat)
      return
    endif

    ! Diagnostic group
    if (have_diagnostic_group) then
      call tiof_push_group (obj, o3p_grp_diagnostic, errstat)
    endif
    ! Optional parameters
    ! fitted wavelengths
    if (.not. reduce_resolution) then
      call tiof_put3d_r4 (obj, o3p_var_wavel, [min_step, min_xtrack, 0], &
           [nstep, nxtrack, nmax_wavs], &
                          wavelengths(1:nmax_wavs,1:nxtrack,1:nstep), errstat)
    endif
    !weighting function
    if (ozwrtwf) then
      call tiof_put4d_r4 (obj, o3p_var_weight_func, &
           [min_step, min_xtrack, 0, 0], &
           [nstep, nxtrack, nmax_wavs, nfitvars], &
                          wgt_func(1:nfitvars,1:nmax_wavs,1:nxtrack,1:nstep), errstat)
    endif
     ! normalized observed & simulated spectra
    if (ozwrtres) then
      call tiof_put3d_r4 (obj, o3p_var_sim_norm_rad, &
           [min_step, min_xtrack, 0], &
           [nstep, nxtrack, nmax_wavs], &
                          sim_norm_rad(1:nmax_wavs,1:nxtrack,1:nstep), errstat)
      call tiof_put3d_r4 (obj, o3p_var_norm_radiance, &
           [min_step, min_xtrack, 0], &
           [nstep, nxtrack, nmax_wavs], norm_rad(1:nmax_wavs,1:nxtrack,1:nstep), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: writing to diagnostic group", &
                       errstat)
      return
    endif

    ! QA stats group
    call tiof_push_group (obj, o3p_grp_qa_stats, errstat)
    call tiof_put2d_i4 (obj, o3p_var_exit_stat, [min_step, min_xtrack], &
         [nstep, nxtrack], exval(1:nxtrack,1:nstep), errstat)
    call tiof_put2d_i4 (obj, o3p_var_iter, [min_step, min_xtrack], &
         [nstep, nxtrack], iterations(1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_fit_rms, [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nfitwins], rms(1:nfitwins,1:nxtrack,1:nstep), errstat)
    call tiof_put3d_r4 (obj, o3p_var_avg_resid, [min_step, min_xtrack, 0], &
         [nstep, nxtrack, nfitwins], avg_resid(1:nfitwins,1:nxtrack,1:nstep), errstat)
    !call tiof_put1d_i2 (obj, o3p_var_mqf, [min_step], [nstep], mqf, errstat)
    !Optional param - relative measurement error
    if (ozwrtsnr) then
      call tiof_put3d_r4 (obj, o3p_var_merr, [min_step, min_xtrack, 0], &
           [nstep, nxtrack, nmax_wavs], fit_wgt(1:nmax_wavs,1:nxtrack,1:nstep), errstat)
    endif
    call tiof_pop_group (obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_write_error, &
                       "l2_tio_write_data: writing to qa_stats group", &
                       errstat)
      return
    endif

  end subroutine write_merged_data

  !> Copy metadata required for processing from L1B input file
  !! @param[in]    l1bfile  Filename for input radiance file
  !! @param[inout] errstat  Error status variable
  subroutine copy_hdr_metadata (l1bfile, errstat)
    implicit none
    character (len=*), intent(in) :: l1bfile
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj
    type (tiof_file_type) :: l1b

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_open (l1bfile, l1b, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "copy_metadata: opening file "//trim(l1bfile), &
                       errstat)
      return
      endif

    call tiof_copy_granule_ident (l1b, obj, errstat)
    call tiof_close (l1b, errstat)

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "copy_metadata: copying from "//trim(l1bfile), &
                       errstat)
    endif
  end subroutine copy_hdr_metadata

  subroutine copy_l2_metadata (l2file_in, errstat)
    implicit none
    character (len=*), intent(in) :: l2file_in
    integer, intent(inout) :: errstat
    type (tiof_file_type), pointer :: obj
    type (tiof_file_type) :: l2_in

    character (len=*), parameter, dimension(*) :: &
      attr_names = (/'geospatial_bounds      ', &
                     'geospatial_bounds_crs  ', &
                     'geospatial_lon_min     ', &
                     'geospatial_lon_max     ', &
                     'geospatial_lat_min     ', &
                     'geospatial_lat_max     ', &
                     'input_files            ', &
                     'local_granule_id       ', &
                     'version_id             ', &
                     'apriori_source         '/)

    if (errstat /= 0) return

    obj => primary_output_file

    call tiof_open (l2file_in, l2_in, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, "copy_l2_metadata: opening file "//trim(l2file_in), &
                       errstat)
      return
    endif

    call tiof_copy_attr_id (l2_in, nf90_global, obj, nf90_global, attr_names, errstat)
    call tiof_close (l2_in, errstat)

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "copy_l2_metadata: copying from "//trim(l2file_in), &
                       errstat)
    endif
  end subroutine copy_l2_metadata

  !> Label product type in Level 2 product file
  !! @param[in]    label   product type label to apply
  !! @param{in]    processing_version processing version used in file creation
  !! @param[inout] errstat Error status variable
  subroutine label_output_file (label, processing_version, errstat)
    implicit none
    character (len=*), intent(in) :: label
    integer, intent(in) :: processing_version
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj

    obj => primary_output_file

    call tiof_label_product (obj, label, 2, processing_version, errstat)
    if (errstat < 0) then
      call tell_error (tell_io_error, "label_output_file failed", errstat)
    endif

  end subroutine label_output_file

  subroutine write_apriori_attr (errstat)
    use OMSAO_variables_module, only: apriori_source
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_attlist_type) :: attlist

    if (errstat /= 0) return

    obj => primary_output_file
    call tiof_attlist_append (attlist, errstat, "apriori_source", &
                              att_text=apriori_source)
    call tiof_def_atts (obj, attlist, nf90_global, errstat)
    call tiof_attlist_free (attlist)
    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "write_apriori_attr: failed", errstat)
      return
    endif

  end subroutine write_apriori_attr

  subroutine set_production_date_time (errstat)
    implicit none
    integer, intent(inout) :: errstat

    type (tiof_file_type), pointer :: obj
    type (tiof_attlist_type) :: attlist
    character (len=8) :: date
    character (len=10) :: time
    character (len=5) :: zone
    character (len=32) :: prod_datetime

    if (errstat /= 0) return

    obj => primary_output_file

    call date_and_time(DATE=date, TIME=time, ZONE=zone)
    prod_datetime=date(1:4)//'-'//date(5:6)//'-'//date(7:8) &
         //'T'//time(1:2)//':'//time(3:4)//':'//time(5:10)//' UTC'//zone

    call tiof_attlist_append (attlist, errstat, "production_date_time", &
         att_text=prod_datetime)
    call tiof_push_group (obj, "/", errstat)
    call tiof_def_atts (obj, attlist, nf90_global, errstat)
    call tiof_pop_group (obj, errstat)
    call tiof_attlist_free (attlist)

    if (errstat /= 0) then
      call tell_error (tell_io_write_error, &
                       "set_production_date_time: failed", errstat)
      return
    endif

  end subroutine set_production_date_time

end module tio_output_module
