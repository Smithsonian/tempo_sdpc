program test_attr

  use netcdf
  use tio_module
  use tell_module
  use md_module

  implicit none

  !Fixed inputs
  integer (kind=4), parameter :: nxtrack=21, nstep=11, ninp=5, version=1
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
  real (kind=4), dimension (nxtrack,nstep) :: lat, lon, vza
  integer (kind=4), dimension (nxtrack,nstep) :: inrqf
  integer(kind=4), dimension(2) :: dimid_2d
  real (kind=4), dimension(:), allocatable :: bdry_lon, bdry_lat
  real (kind=4) :: centroid_lon, centroid_lat

  type (tiof_dimlist_type) :: dimlist
  type (tiof_varlist_type) :: varlist
  type (tiof_file_type) :: obj

  namelist /metadata/  collection_shortname, collection_version, platform, &
         access_description, access_value, abstract, keywords

  !--------------------------------------------------------------

  errstat = 0

  ! read from namelist
  open (333, file=nlfile, status='OLD', iostat=errstat)
  if (errstat /= 0) then
    write (*,*)'*** error opening file: '//trim(nlfile)
    stop 1
  endif

  read (333, nml=metadata, iostat=errstat)
  if (errstat /= 0) then
    write (*,*)'*** error reading file: '//trim(nlfile)
    stop 1
  endif

  close (333, iostat=errstat)
  if (errstat /= 0) then
    write (*,*)'*** error closing file: '//trim(nlfile)
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
      vza(m,n)=0.0
      inrqf(m,n) = 0
    enddo
  enddo

  ! build a netCDF file with some basic data and a  metadata group
  call tiof_dimlist_append (dimlist, "mirror_step", nstep, errstat)
  call tiof_dimlist_append (dimlist, "xtrack", nxtrack, errstat)
  call tiof_create (obj, l2file, nf90_clobber, errstat)
  call tiof_def_group (obj, grp_md, errstat)
  call tiof_def_group (obj, grp_geo, errstat)
  call tiof_def_dims (obj, dimlist, errstat)

  if (errstat /= 0) then
    write(*,*)'*** test_attr: tiof_create failed'
    stop 1
  endif

  step_indices=[(n, n=0, nstep-1)]
  xtrack_indices=[(m, m=0, nxtrack-1)]

  call tiof_dimlist_lookup (dimlist, ["     xtrack", "mirror_step"], dimid_2d,&
       errstat)
  call tiof_varlist_append (varlist, errstat, "mirror_step",&
       nf90_int, dimids=[dimid_2d(2)])
  call tiof_varlist_append (varlist, errstat, "xtrack", &
       nf90_int, dimids=[dimid_2d(1)])
  call tiof_def_vars (obj, varlist, errstat)
  call tiof_varlist_free (varlist)
  call tiof_put1d_i4 (obj, "xtrack", [0], [nxtrack], xtrack_indices, errstat)
  call tiof_put1d_i4 (obj, "mirror_step", [0], [nstep], step_indices, errstat)

  if (errstat /= 0) then
    write (*,*)'*** test_attr: failed to define coordinate vars'
    stop 1
  endif

  call tiof_varlist_append (varlist, errstat, tempo_var_longitude, nf90_float, &
       dimids=dimid_2d)
  call tiof_varlist_append (varlist, errstat, tempo_var_latitude, &
       nf90_float, dimids=dimid_2d)
  call tiof_varlist_append (varlist, errstat, tempo_var_vz_angle, &
       nf90_float, dimids=dimid_2d)
  call tiof_varlist_append (varlist, errstat, tempo_var_inrqf, &
       nf90_int, dimids=dimid_2d)

  call tiof_push_group (obj, grp_geo, errstat)
  call tiof_def_vars (obj, varlist, errstat)

  call tiof_put2d_r4 (obj, tempo_var_longitude, [0,0], &
       [nstep, nxtrack], lon(1:nxtrack,1:nstep), errstat)
  call tiof_put2d_r4 (obj, tempo_var_latitude, [0,0], &
       [nstep, nxtrack], lat(1:nxtrack,1:nstep), errstat)
  call tiof_put2d_r4 (obj, tempo_var_vz_angle, [0,0], &
       [nstep, nxtrack], vza(1:nxtrack,1:nstep), errstat)
  call tiof_put2d_i4 (obj, tempo_var_inrqf, [0,0], &
       [nstep, nxtrack], inrqf(1:nxtrack,1:nstep), errstat)

  call tiof_pop_group (obj, errstat)
  call tiof_varlist_free (varlist)

  if (errstat /= 0) then
    write (*,*)'*** test_attr: failed to write data vars'
    stop 1
  endif

  call tiof_label_product (obj, "test", 1, errstat)
  call tiof_close (obj, errstat)

  if (errstat /= 0) then
    write(*,*)'*** test_attr: failed to define metadata group'
    stop 1
  endif

  call tiof_open (l2file, obj, nf90_write, errstat)
  call tiof_push_group (obj, grp_geo, errstat)
  call tiof_make_lev1_bounding_polygon (obj, bdry_lon, bdry_lat, centroid_lon, centroid_lat, errstat)
  call tiof_close (obj, errstat)

  if (errstat /= 0) then
    call tell_error (tell_runtime_error, "generating bounding polygon", errstat)
    stop 1
  endif

  if (abs(centroid_lon) > 1.e-4) then
    write(*,*)'*** test_attr:  unexpected centroid_lon = ', centroid_lon
    stop 1
  endif
  if (abs(centroid_lat) > 1.e-4) then
    write(*,*)'*** test_attr:  unexpected centroid_lat = ', centroid_lat
    stop 1
  endif

  ! add the metadata
  call md_open (l2file, errstat)
  call md_write_geo_bounds (bdry_lon, bdry_lat, centroid_lon, centroid_lat, errstat)
  call md_write_inputs (ninp, inputs, errstat)
  call md_write_fixed (nlfile,errstat)
  call md_write_prodid (l2file,version,errstat)
  call md_close (errstat)

  deallocate (bdry_lon, bdry_lat)

  if (errstat /= 0) then
    write(*,*)'*** test_attr: failed to write metadata'
    stop 1
  endif

end program test_attr
