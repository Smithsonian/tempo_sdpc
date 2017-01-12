program merge_o3p_files

  use netcdf
  use tio_module
  use tell_module
  use o3p_names_module
  use m_read_L2_o3p_tio
  use m_o3p_params
  use tio_output_module, only: l2_tio_create, l2_tio_close, write_merged_geo, &
       write_merged_data
  use ozprof_data_module, only: ozwrtavgk, ozwrtcorr, ozwrtcovar, &
       ozwrtcontri, ozwrtres, ozwrtwf, ozwrtsnr, &
       ozwrtvar, gaswrt, aerosol, do_lambcld
  use OMSAO_variables_module, only: reduce_resolution

  implicit none


  integer (kind=4) :: errstat
  ! input files
  ! note fixed max size of input_files array, required for namelist input
  integer (kind=4) :: ninput
  character (len=128), dimension(64) :: input_files
  ! output file
  character (len=128) :: outfile
  ! dimension indices
  integer (kind=4) :: min_step, max_step, min_xtrack, max_xtrack, &
       nstep_tot, nxtrack_tot, start_wstep, end_wstep, start_wxtrack, &
       end_wxtrack
  integer (kind=4) :: min_sf, max_sf, min_xf, max_xf
  integer (kind=4), dimension(:,:), allocatable :: step, xtrack, dup_check
  integer :: n, status, dummyid, i, j
  integer, parameter :: zero=0, one=1

  type (tiof_file_type) :: tio_l2in

  !input & output filenames entered via namelist
  namelist /merge_o3p_iolist/ ninput, input_files, outfile

  !--------------------------------------------------------------------

  errstat = 0
  call tell_open ("merge_o3p_files", 0)

  ! read filenames from namelist
  open (777, file='merge_o3p_iolist.nml', status='OLD', iostat=errstat)
  read (777, nml=merge_o3p_iolist, iostat=errstat)
  close(777, iostat=errstat)
  if (errstat /= 0) then
    call tell_error (tell_io_read_error, &
         "failed to read namelist merge_o3p_iolist.nml", errstat)
    stop 1
  endif

  min_step=0
  max_step=0
  min_xtrack=0
  max_xtrack=0
  ! allocate dmension size arrays to match number of input files
  call o3p_dim_alloc (ninput, errstat)
  if (errstat /= 0) stop 1

  ! Loop through input files, get dimension sizes
  do n=1,ninput
  ! From first file determine which dimensions/variables are used
    call read_o3p_dims(input_files(n), tio_l2in, nstep(n), nxtrack(n), &
         ncorner(n), nfitvars(n), nfitwins(n), ngas(n), nlayer(n), &
         nlayerp1(n), nmax_wavs(n), nnoise_elems(n), nnongas(n), &
         naeros_wavs(n), errstat)

    if (n == 1) then
      if (ngas(n) > 0) gaswrt = .true.
      if (nnongas(n) > 0) ozwrtvar = .true.
      if (nnoise_elems(n) > 0) ozwrtcovar = .true.
      if (naeros_wavs(n) > 0) aerosol = .true.
      call open_o3p(input_files(n), tio_l2in, errstat)
      call tiof_push_group (tio_l2in, o3p_grp_support_data, errstat)
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_o3_avg_kernel, dummyid)
      if (status == nf90_noerr) ozwrtavgk = .true.
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_correl, dummyid)
      if (status == nf90_noerr) ozwrtcorr = .true.
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_contrib_func, dummyid)
      if (status == nf90_noerr) ozwrtcontri = .true.
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_cld_opt_depth, dummyid)
      if (status == nf90_noerr) then
        do_lambcld = .false. !note false means include cld_opt_depth
      else
        do_lambcld = .true.
      endif
      call tiof_pop_group (tio_l2in, errstat)
      call tiof_push_group (tio_l2in, o3p_grp_qa_stats, errstat)
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_fit_weight, dummyid)
      if (status == nf90_noerr) ozwrtsnr = .true.
      call tiof_pop_group(tio_l2in, errstat)
      call tiof_push_group (tio_l2in, o3p_grp_diagnostic, errstat)
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_weight_func, dummyid)
      if (status == nf90_noerr) ozwrtwf = .true.
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_norm_radiance, dummyid)
      if (status == nf90_noerr) ozwrtres = .true.
      status = nf90_inq_varid(tio_l2in%groupid, o3p_var_wavel, dummyid)
      if (status == nf90_noerr) then
        reduce_resolution = .false. !note false means include wavelengths
      else
        reduce_resolution = .true.
      endif
      call close_o3p(tio_l2in, errstat)
      if (errstat /= 0) stop 1
    endif
  enddo

  do n=1,ninput ! loop over input files, get dimension indices
  ! For each file, read in xtrack and step values
    if (n == 1) then
      allocate(step(maxval(nstep), ninput), stat=errstat)
      allocate(xtrack(maxval(nxtrack), ninput), stat=errstat)
      if (errstat /= 0) then
        call tell_error(tell_malloc_error, &
             "failed to allocate step and xtrack arrays", errstat)
        stop 1
      endif
    endif
    call read_o3p_dim_indices (input_files(n), tio_l2in, nstep(n), &
         nxtrack(n), step(:,n), xtrack(:,n), errstat)

  ! check other dimensions have same sizes
    if (n > 1) then
      if (ncorner(n) /= ncorner(1)) errstat = -1
      if (nfitvars(n) /= nfitvars(1)) errstat = -2
      if (nfitwins(n) /= nfitwins(1)) errstat = -3
      if (ngas(n) /= ngas(1)) errstat = -4
      if (nlayer(n) /= nlayer(1)) errstat = -5
      if (nlayerp1(n) /= nlayerp1(1)) errstat = -6
      if (nmax_wavs(n) /= nmax_wavs(1)) errstat = -7
      if (nnoise_elems(n) /= nnoise_elems(1)) errstat = -8
      if (nnongas(n) /= nnongas(1)) errstat = -9
      if (naeros_wavs(n) /= naeros_wavs(1)) errstat = -10
      if (errstat /= 0) then
        call tell_error(tell_invalid_parm_error, &
             "mismatched dimension size in input files", errstat)
        stop 1
      endif
    endif
  enddo  ! loop over input files



  ! determine min max indices of step & xtrack
  max_step=maxval(step)
  min_step=minval(step)
  max_xtrack=maxval(xtrack)
  min_xtrack=minval(xtrack)
  nstep_tot=max_step-min_step+1
  nxtrack_tot=max_xtrack-min_xtrack+1

  ! check for overlaps - files should cover unique areas in step, xtrack
  allocate(dup_check(0:max_step, 0:max_xtrack), stat=errstat)
  if (errstat /= 0) then
    call tell_error(tell_malloc_error, "failed to allocate dup_check", errstat)
    stop 1
  endif
  dup_check=0
  do n=1,ninput
    do i=1,nstep(n)
      do j=1,nxtrack(n)
        dup_check(step(i,n),xtrack(j,n))=1
      enddo
    enddo
  enddo
  if (maxval(dup_check) > 1) then
    call tell_error(tell_runtime_error, &
         "input files have overlapping step, xtrack indices", errstat)
    stop 1
  endif

  ! allocate output data arrays and insert fill values
  call o3p_param_alloc (one, maxval(nstep), one, maxval(nxtrack), &
       ncorner(1), nfitvars(1), nfitwins(1), ngas(1), nnongas(1), nlayer(1), &
       nlayerp1(1), nmax_wavs(1), nnoise_elems(1), naeros_wavs(1), &
       errstat)
!  call o3p_param_fill (errstat)
  if (errstat /= 0) stop 1

  ! Create output file
  call l2_tio_create (outfile, min_xtrack+1, max_xtrack+1, min_step+1, &
       max_step+1, zero, one, ngas(1), nlayer(1), nfitvars(1), nfitwins(1), &
       nnongas(1), nmax_wavs(1), errstat)
  if (errstat /= 0) stop 1

  ! For each file, read in data to appropriate section of output arrays
  do n=1,ninput
    min_sf = minval(step(:,n))
    max_sf = maxval(step(:,n))
    min_xf = minval(xtrack(:,n))
    max_xf = maxval(xtrack(:,n))
    start_wstep=min_sf-min_step
    end_wstep=max_sf-min_step
    start_wxtrack=min_xf-min_xtrack
    end_wxtrack=max_xf-min_xtrack

    call o3p_param_fill (errstat)

    call open_o3p(input_files(n), tio_l2in, errstat)

    call read_o3p_geolocation(tio_l2in, nstep(n), nxtrack(n), ncorner(n), &
         one, nxtrack(n), one, nstep(n), errstat)

    call read_o3p_product(tio_l2in, nstep(n), nxtrack(n), nlayer(n), &
         ngas(n), nnongas(n), one, nxtrack(n), one, nstep(n), errstat)

    call read_o3p_support(tio_l2in, nstep(n), nxtrack(n), nlayer(n), &
         nlayerp1(n), nfitwins(n), ngas(n), nnongas(n), nfitvars(n), &
         nnoise_elems(n), nmax_wavs(n), naeros_wavs(n), &
         one, nxtrack(n), one, nstep(n), errstat)

    call read_o3p_diagnostic (tio_l2in, nstep(n), nxtrack(n), nmax_wavs(n), &
         nfitvars(n), one, nxtrack(n), one, nstep(n), errstat)

    call read_o3p_qastat (tio_l2in, nstep(n), nxtrack(n), nfitwins(n), &
         nmax_wavs(n), one, nxtrack(n), one, nstep(n), errstat)

    call close_o3p (tio_l2in, errstat)

    if (errstat /= 0) stop 1

    !write values
    call write_merged_geo (start_wstep, end_wstep, start_wxtrack, &
         end_wxtrack, ncorner(1), errstat)

    call write_merged_data (start_wstep, end_wstep, start_wxtrack, &
         end_wxtrack, ngas(1), &
         nnongas(1), nlayer(1), nfitvars(1), nfitwins(1), nmax_wavs(1), &
         nnoise_elems(1), naeros_wavs(1), errstat)
    
    if (errstat /= 0) stop 1

  enddo ! read/write data loop

  call l2_tio_close (errstat)
  if (errstat /= 0) stop 1

  ! deallocate memory
  call o3p_param_dealloc (errstat)
  call o3p_dim_dealloc (errstat)
  deallocate(dup_check, step, xtrack, stat=errstat)
  if (errstat /= 0) then
    call tell_error(tell_malloc_error, "deallocation failure", errstat)
    stop 1
  endif

end program merge_o3p_files
