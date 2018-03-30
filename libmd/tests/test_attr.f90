program test_attr

  use netcdf
  use tio_module
  use tell_module
  use md_module

  implicit none

  !Fixed inputs
  integer (kind=4), parameter :: nxtrack=21, nstep=11, ninp=5
  character (len=3), parameter :: version_str='(1)'
  character (len=32), dimension(ninp) :: inputs
  character (len=16), parameter :: l2file='test_attr.nc', &
       nlfile='boilerplate.nml'

  integer, parameter :: cld_max_name_len = 64
  character (len=cld_max_name_len), parameter :: grp_geo = "/geolocation"
  character (len=cld_max_name_len), parameter :: grp_md = "/metadata"

  !local variables read from namelist
  character (len=32) :: collection_shortname, platform
  character (len=5) :: collection_version
  character (len=2048) :: abstract
  character (len=128) :: access_description
  character (len=32) :: access_value
  character (len=1024) :: keywords

  integer (kind=4) :: n, m, errstat
  integer (kind=4), dimension(0:nstep-1) :: step_indices
  integer (kind=4), dimension (0:nxtrack-1) :: xtrack_indices
  real (kind=4), dimension (nxtrack,nstep) :: lat, lon
  integer(kind=4), dimension(2) :: dimid_2d

  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist
  type (tiof_file_type) :: obj

  namelist /metadata/  collection_shortname, collection_version, platform, &
         access_description, access_value, abstract, keywords


  !--------------------------------------------------------------

  errstat = 0

  ! read from namelist
  open (333, file=nlfile, status='OLD', iostat=errstat)
  if (errstat == 0) read (333, nml=metadata, iostat=errstat)
  if (errstat == 0) close (333, iostat=errstat)
  if (errstat /= 0) then
    call tell_error (tell_io_read_error, &
         "*** test_attr: failed to read from namelist", errstat)
    stop 1
  endif

  ! set input files
  inputs(1)='l1_radiance.nc'
  inputs(2)='l1_irradiance.nc'
  inputs(3)='lookup.txt        '
  inputs(4)='reference_table.nc'
  inputs(5)='   otherstuff.dat '

  ! set up lat, lon arrays
  do n=1,nstep
    do m=1,nxtrack
      lat(m,n)=-10+(m-1)
      lon(m,n)=-10+(2*n-2)
    enddo
  enddo

  ! build a netCDF file with some basic data and a  metadata group
  call tiof_dimlist_append (dimlist, "mirror_step", nstep, errstat)
  call tiof_dimlist_append (dimlist, "xtrack", nxtrack, errstat)
  if (errstat == 0) call tiof_create (obj, l2file, nf90_clobber, errstat)
  if (errstat == 0) call tiof_def_group (obj, grp_md, errstat)
  if (errstat == 0) call tiof_def_group (obj, grp_geo, errstat)
  if (errstat == 0) call tiof_def_dims (obj, dimlist, errstat)
  if (errstat /= 0) then
    write(*,*)'*** test_attr: tiof_create failed'
    stop 1
  endif

  step_indices=[(n, n=0, nstep-1)]
  xtrack_indices=[(m, m=0, nxtrack-1)]

  call tiof_dimlist_lookup (dimlist, ["     xtrack", "mirror_step"], dimid_2d,&
       errstat)
  if (errstat == 0) call tiof_varlist_append (varlist, errstat, "mirror_step",&
       nf90_int, dimids=[dimid_2d(2)])
  if (errstat == 0) call tiof_varlist_append (varlist, errstat, "xtrack", &
       nf90_int, dimids=[dimid_2d(1)])
  if (errstat == 0) call tiof_def_vars (obj, varlist, errstat)
  if (errstat == 0) call tiof_varlist_free (varlist)
  if (errstat == 0) call tiof_put1d_i4 (obj, "xtrack", [0], [nxtrack], &
         xtrack_indices, errstat)
  if (errstat == 0) call tiof_put1d_i4 (obj, "mirror_step", [0], [nstep], &
         step_indices, errstat)
  if (errstat /= 0) then
    write (*,*)'*** test_attr: failed to define coordinate vars'
    stop 1
  endif

  call tiof_varlist_append (varlist, errstat, "lon", nf90_float, &
       dimids=dimid_2d)
  if (errstat == 0) call tiof_varlist_append (varlist, errstat, "lat", &
       nf90_float, dimids=dimid_2d)
  if (errstat == 0) call tiof_push_group (obj, grp_geo, errstat)
  if (errstat == 0) call tiof_def_vars (obj, varlist, errstat)
  if (errstat == 0) call tiof_put2d_r4 (obj, "lon", [0,0], &
       [nstep, nxtrack], lon(1:nxtrack,1:nstep), errstat)
  if (errstat == 0) call tiof_put2d_r4 (obj, "lat", [0,0], &
       [nstep, nxtrack], lat(1:nxtrack,1:nstep), errstat)
  if (errstat == 0) call tiof_pop_group (obj, errstat)
  if (errstat == 0) call tiof_varlist_free (varlist)
  if (errstat /= 0) then
    write (*,*)'*** test_attr: failed to write data vars'
    stop 1
  endif
!
  call tiof_label_product (obj, "test", 1, errstat)
  if (errstat == 0) call tiof_close (obj, errstat)
  if (errstat /= 0) then
    write(*,*)'*** test_attr: failed to define metadata group'
    stop 1
  endif

  ! add the metadata
  call open_md(l2file, errstat)
  if (errstat == 0) call write_geo_bounds_md(nxtrack, nstep, lat, lon, errstat)
  if (errstat == 0) call write_inputs_md(ninp, inputs, errstat)
  if (errstat == 0) call write_fixed_md(nlfile,errstat)
  if (errstat == 0) call write_prodid_md(l2file,version_str,errstat)
  if (errstat == 0) call close_md(errstat)
  if (errstat /= 0) then
    write(*,*)'*** test_attr: failed to write metadata'
    stop 1
  endif

end program test_attr
