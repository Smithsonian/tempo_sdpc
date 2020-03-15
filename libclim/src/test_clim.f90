program test_clim
  use clim_module
  implicit none
  type (clim_pres_type) :: cpt
  type (clim_species_type) :: cst
  type (clim_pres_bounds_type) :: bounds
  type (clim_cloud_type) :: cct
  integer :: month, day, nz, i, errstat
  real (kind=4) :: hour, lon, lat, psurf, ptrop, cloud_pressure
  real (kind=4), dimension(:), allocatable :: pres_z, ap, bp
  real (kind=4), dimension(:), allocatable :: vmr_z, partial_column_z
  real (kind=4), dimension(:), allocatable :: vmr_stddev_z
  character (len=*), parameter :: species = 'NO2'
  logical, parameter :: select_vmr_stddev = .true.

  errstat = 0

  call clim_query_nz (nz, errstat)
  if (errstat /= 0) call exit(1)

  bounds % hour_beg = 18.0
  bounds % hour_end = 18.0 + 6.0/60
  bounds % lon_min = -90.0
  bounds % lon_max = -80.0
  bounds % lat_min = 15.0
  bounds % lat_max = 70.0

  month = 7
  day   = 1

  call clim_pres_init (cpt, month, day, bounds, errstat)
  if (errstat /= 0) call exit(1)

  nz = clim_pres_nz (cpt)

  allocate (pres_z(nz), ap(nz), bp(nz))

  hour = 18.0
  lon  = -85.0
  lat  = +36.0

  call clim_cloud_init (cct, errstat)
  if (errstat /= 0) call exit(1)

  call clim_cloud (cct, month, day, lon, lat, cloud_pressure, errstat)
  if (errstat /= 0) call exit(1)

  write (*,'(a,f10.4,a)')'cloud pressure = ',cloud_pressure,' hPa'

  write (*,*)'species = ',species
  write (*,'(a,i2,1x,i2,1x,f7.4)')' month, day, hour = ', month, day, hour
  write (*,'(a,f10.4,1x,f10.4)')' lon, lat = ', lon, lat

  call clim_pres_eta (cpt, ap, bp, errstat)
  call clim_pres (cpt, hour, lon, lat, pres_z, errstat, &
                  p_surf=psurf, p_trop=ptrop)
  if (errstat /= 0) call exit(1)

  write (*,*)'P(surface) = ',psurf,' hPa'
  write (*,*)'P(tropopause) = ',ptrop,' hPa'

  ! select_vmr_stddev is an optional, logical argument.
  ! If select_vmr_stddev=.true. the VMR standard deviation will be read
  ! in from the disk files and will be available in calls to clim_species_vmr
  call clim_species_init (cst, cpt, species, errstat, select_vmr_stddev)
  if (errstat /= 0) call exit(1)

  allocate (vmr_z(nz-1), partial_column_z(nz-1), vmr_stddev_z(nz-1))

  ! vmr_stddev_z is an optional argument.
  ! If clim_species_init was called with select_vmr_stddev=.true., then
  ! the VMR standard deviation, vmr_stddev_z, will be interpolated.
  ! Here, I'm calling clim_species_vmr both with and without the optional
  ! argument only for testing purposes.
  call clim_species_vmr (cst, cpt, hour, lon, lat, vmr_z, errstat)
  if (errstat /= 0) call exit(1)
  call clim_species_vmr (cst, cpt, hour, lon, lat, vmr_z, errstat, vmr_stddev_z)
  if (errstat /= 0) call exit(1)

  call clim_partial_column (pres_z, vmr_z, partial_column_z, errstat)
  if (errstat /= 0) call exit(1)

  write(*,*)'N(total) = ',sum(partial_column_z),' mol/cm^2'

  write(*,*)' i  pres_z        ap            bp            '// &
    'vmr_z         dN [cm^-2]    vmr_stddev_z'
  do i=1,nz-1
    write(*,'(i3,10(2x,1pe12.6))')i,pres_z(i),ap(i),bp(i), &
      vmr_z(i),partial_column_z(i),vmr_stddev_z(i)
  enddo
  i=nz
  write(*,'(i3,10(2x,1pe12.6))')i,pres_z(i),ap(i),bp(i)

  deallocate(pres_z, ap, bp, vmr_z, partial_column_z, vmr_stddev_z)

end program
