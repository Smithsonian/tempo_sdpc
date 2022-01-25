program test_clim
  use tell_module
  use clim_module
  implicit none
  type (clim_pres_type) :: cpt
  type (clim_val_type) :: cst, temp_cst, u2m_cst, v2m_cst
  type (clim_pres_bounds_type) :: bounds
  type (clim_cloud_type) :: cct
  integer :: year, month, day, nz, i, errstat
  real (kind=4) :: hour, lon, lat, psurf, ptrop, cloud_pressure
  real (kind=4), dimension(:), allocatable :: pres_z, ap, bp
  real (kind=4), dimension(:), allocatable :: vmr_z, partial_column_z, temp_z
  real (kind=4), dimension(1) :: u2m, v2m
  character (len=*), parameter :: species = 'NO2'

  errstat = 0

  call tell_open ("test_clim", 0)
  call tell_set_log_level (1)

  call clim_read_config ('clim_config.ini', errstat)
  if (errstat /= 0) call exit(1)

  call clim_query_nz (nz, errstat)
  if (errstat /= 0) call exit(1)

  bounds % hour_beg = 18.0
  bounds % hour_end = bounds % hour_beg + 6.0/60
  bounds % lon_min = -90.0
  bounds % lon_max = -80.0
  bounds % lat_min = 15.0
  bounds % lat_max = 70.0

  year = 2022
  month = 1
  day   = 25

  call clim_pres_init (cpt, year, month, day, bounds, errstat)
  if (errstat /= 0) call exit(1)

  allocate (pres_z(nz+1), ap(nz+1), bp(nz+1))

  hour = 18.0
  lon  = -85.0
  lat  = +36.0

  call clim_cloud_init (cct, errstat)
  if (errstat /= 0) call exit(1)

  call clim_cloud (cct, month, day, lon, lat, cloud_pressure, errstat)
  if (errstat /= 0) call exit(1)

  write (*,'(a,f10.4,a)')'cloud pressure = ',cloud_pressure,' hPa'

  write (*,*)'species = ',species
  write (*,'(a,i4)')' year = ',year
  write (*,'(a,i2,1x,i2,1x,f7.4)')' month, day, hour = ', month, day, hour
  write (*,'(a,f10.4,1x,f10.4)')' lon, lat = ', lon, lat

  call clim_pres_eta (ap, bp, errstat)
  call clim_pres (cpt, hour, lon, lat, pres_z, errstat, &
                  p_surf=psurf, p_trop=ptrop)
  if (errstat /= 0) call exit(1)

  write (*,*)'P(surface) = ',psurf,' hPa'
  write (*,*)'P(tropopause) = ',ptrop,' hPa'

  call clim_val_init (u2m_cst, cpt, 'U2M', errstat, single_layer=.true.)
  call clim_val_init (v2m_cst, cpt, 'V2M', errstat, single_layer=.true.)
  if (errstat /= 0) call exit(1)
  call clim_val_interp (u2m_cst, cpt, hour, lon, lat, u2m, errstat)
  call clim_val_interp (v2m_cst, cpt, hour, lon, lat, v2m, errstat)
  if (errstat /= 0) call exit(1)
  write (*,'(a,f7.3,a)')'U2M = ',u2m,' m/s'
  write (*,'(a,f7.3,a)')'V2M = ',v2m,' m/s'

  call clim_val_init (cst, cpt, species, errstat)
  if (errstat /= 0) call exit(1)

  allocate (vmr_z(nz), partial_column_z(nz), temp_z(nz))

  call clim_val_interp (cst, cpt, hour, lon, lat, vmr_z, errstat)
  if (errstat /= 0) call exit(1)

  call clim_val_init (temp_cst, cpt, 'T', errstat)
  if (errstat /= 0) call exit(1)
  call clim_val_interp (temp_cst, cpt, hour, lon, lat, temp_z, errstat)

  call clim_partial_column (pres_z, vmr_z, partial_column_z, errstat)
  if (errstat /= 0) call exit(1)

  write(*,*)'N(total) = ',sum(partial_column_z),' mol/cm^2'

  write(*,*)' i  pres_z        ap            bp            '// &
    'vmr_z         dN_z [cm^-2]  T_z[K]'
  do i=1,nz
    write(*,'(i3,10(2x,1pe12.6))')i,pres_z(i),ap(i),bp(i),vmr_z(i), &
      partial_column_z(i), temp_z(i)
  enddo
  i=nz+1
  write(*,'(i3,10(2x,1pe12.6))')i,pres_z(i),ap(i),bp(i)

  deallocate(pres_z, ap, bp, vmr_z, partial_column_z, temp_z)

end program
