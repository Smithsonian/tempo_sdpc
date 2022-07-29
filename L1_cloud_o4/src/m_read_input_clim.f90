! Read GEOS-CF input climatology
module m_read_input_clim
   use clim_module
   use tell_module

   public read_geoscf

   private prepare_geoscf
   real(kind=4), parameter :: r4_missval = -1.0e+30
   real(kind=8), parameter :: r8_missval = -1.0d+30
contains

   subroutine prepare_geoscf(cpt,temp_cst,u2m_cst,v2m_cst,&
              lon_min,lon_max,lat_min,lat_max, errstat)


   use m_vars, only: geos_np, nlayers, gmetadata, &
              rad_latitude, rad_longitude, rad_time

  implicit none

  type (clim_pres_type), intent(out) :: cpt
  type (clim_val_type), intent (out) :: temp_cst,u2m_cst,v2m_cst
  real(kind=4), intent(out) :: lon_min,lon_max,lat_min,lat_max
  integer, intent(inout) :: errstat

  type (clim_pres_bounds_type) :: bounds

  integer :: year(2), month(2), day(2)
  integer :: thisyear, thismonth, thisday

  integer :: nz
  logical :: have_forecast

  real (kind=8) :: t_beg, t_end, hour_beg, hour_end

  if (errstat /= 0) return

!  write(*,*) 'Preapre to read GEOS-CF for orbit'
  nlayers = geos_np

  t_beg = minval(rad_time, rad_time /= r8_missval)
  t_end = maxval(rad_time, rad_time /= r8_missval)

  call tio_f_taix_time_to_utc_caldate(t_beg, year(1),month(1),day(1),hour_beg)
  call tio_f_taix_time_to_utc_caldate(t_end, year(2),month(2),day(2),hour_end)

  if ((day(1) /= day(2)) .OR. (month(1) /= month(2))  &
      .OR. (year(1) /= year(2))) then
     call tell_error(tell_runtime_error,"time error",errstat)
     return
  endif
  thisyear = year(1)
  thismonth = month(1)
  thisday = day(1)
  bounds%hour_beg = real(hour_beg, kind=4)
  bounds%hour_end = real(hour_end, kind=4)

  errstat = 0

  ! The run-time environment should specify the config file location
  !call clim_read_config ('clim_config.ini', errstat)
  !if (errstat /= 0) call exit(1)

!  lon_min = minval(rad_longitude, rad_longitude /= r4_missval)
  lon_min = -180.0 !initialize
  lon_min = minval(rad_longitude, rad_longitude .gt. -180.0001)
  bounds%lon_min = lon_min
!  lon_max = maxval(rad_longitude, rad_longitude /= r4_missval)
  lon_max = 180.0
  lon_max = maxval(rad_longitude, rad_longitude .lt. 180.0001)
  bounds%lon_max = lon_max
!  lat_min = minval(rad_latitude, rad_latitude /= r4_missval)
  lat_min = -90.
  lat_min = minval(rad_latitude, rad_latitude .gt. -90.0001)
  bounds%lat_min = lat_min
!  lat_max = maxval(rad_latitude, rad_latitude /= r4_missval)
  lat_max = 90.
  lat_max = maxval(rad_latitude, rad_latitude .lt. 90.0001)
  bounds%lat_max = lat_max
  !hqw hour_start & hour_end = gmeta(hour+minute/60.+seconds/3600.)
  write(*,*)'   bounds year,month,day',year(1),month(1),day(1)
  write(*,*)'   bounds hour_beg,end = ',bounds%hour_beg, bounds%hour_end
  write(*,*)'   bounds lon:',bounds%lon_min,bounds%lon_max
  write(*,*)'   bounds lat:',bounds%lat_min,bounds%lat_max

  !hqw check consistency with gmatadata
  !there is slight difference with bounds for lon_min,lon_max,lat_min,lax_max
!  write(*,*) 'gmeta year,month,day:'
!  write(*,*)gmetadata%granule_year,gmetadata%granule_month, gmetadata%granule_day
!  write(*,*) 'gmeta begin hour,minute,seconds=',&
!     gmetadata%granule_hour_start,gmetadata%granule_minute_start,gmetadata%granule_seconds_start
!  write(*,*) 'gmeta end hour,minute,seconds=',&
!     gmetadata%granule_hour_end, gmetadata%granule_minute_end,gmetadata%granule_seconds_end
!  write(*,*)'meta lon:',gmetadata%geospatial_lon_min,gmetadata%geospatial_lon_max
!  write(*,*)'meta lat:',gmetadata%geospatial_lat_min,gmetadata%geospatial_lat_max

  call clim_pres_init(cpt,thisyear,thismonth,thisday,bounds,errstat)

  call clim_query_apriori_source (cpt, have_forecast, errstat)
  call clim_query_nz (nz, errstat)
  if (errstat /= 0) return
  if (nz .NE. geos_np) call exit(1)

  if (have_forecast) then
    gmetadata % apriori_source = 'GEOSCF:forecast'
  else
    gmetadata % apriori_source = 'GEOSCF:climatology'
  endif

  call clim_val_init (temp_cst, cpt, 'T', errstat)
  if (errstat /= 0) call exit(1)

  call clim_val_init (u2m_cst, cpt, 'U2M', errstat, single_layer=.true.)
  call clim_val_init (v2m_cst, cpt, 'V2M', errstat, single_layer=.true.)
  if (errstat /= 0) call exit(1)

  end subroutine prepare_geoscf

!
!------------
!
   subroutine read_geoscf (errstat)

   use m_vars, only: rad_NumTimes,rad_nXtrack, geos_np, &
       rad_latitude, rad_longitude, rad_time, &
       l2_TerrainPressure, windspeed2m, &
       ixdebug, itdebug

   use m_vars, only: geos_Temperature, geos_Pressure

   implicit none

   integer, intent(inout) :: errstat

   type(clim_pres_type) :: cpt
   type(clim_val_type) :: temp_cst, u2m_cst, v2m_cst

   real(kind=8) :: hour
   real(kind=4) :: thislon, thislat, thishour
   real(kind=4) :: lon_min,lon_max,lat_min,lat_max
   real(kind=4) :: pp(geos_np+1), tt(geos_np)

   real (kind=4) :: psurf, ptrop
   real (kind=4), dimension(1) :: u2m, v2m
   real(kind=4), dimension(:), allocatable :: pres_z, temp_z

   integer :: ix, it, iz, nx, nt, kk, nz
   integer :: thisyear, thismonth, thisday

   if (errstat /= 0) return

  ! initialize
  nx = rad_nXtrack
  nt = rad_NumTimes
  nz = geos_np
  write(*,*)'read_geoscf:nXtrack,nTimes,nz=',nx,nt,nz

  call prepare_geoscf(cpt,temp_cst,u2m_cst,v2m_cst,&
                   lon_min,lon_max,lat_min,lat_max, errstat)
  if (errstat /= 0) return

  !allocate geos_temperature and geos_pressure
  allocate(geos_Temperature(nx,nt,nz), stat=errstat)
  allocate(geos_Pressure(nx,nt,nz+1), stat=errstat)

  !allocate windspeed2m
  allocate(windspeed2m(nx,nt), stat=errstat)

  !allocate temporary variables
  allocate(temp_z(nz), pres_z(nz+1), stat=errstat)

  ! loop through pixels
   do it = 1, nt
     if (rad_time(it) == r8_missval) cycle

     call tio_f_taix_time_to_utc_caldate(rad_time(it), &
                    thisyear, thismonth, thisday, hour)
     thishour = real(hour)

     do ix = 1, nx
        thislon = rad_longitude(ix,it)
        thislat = rad_latitude(ix,it)
      ! make sure thislon,thislat are within bounds
      if (thislon .LT. lon_min) thislon = lon_min
      if (thislon .GT. lon_max) thislon = lon_max
      if (thislat .LT. lat_min) thislat = lat_min
      if (thislat .GT. lat_max) thislat = lat_max

      ! get pres_z at lon lat
        call clim_pres (cpt, thishour, thislon, thislat, pres_z, errstat, &
                  p_surf=psurf, p_trop=ptrop)
        if (errstat /= 0) call exit(1)

      ! assign surface pressure (hPa)
        geos_Pressure(ix,it,nz+1) = psurf
        l2_TerrainPressure(ix, it) = psurf

      ! get temperature at lon lat
        call clim_val_interp (temp_cst, cpt, thishour, thislon, thislat, &
                 temp_z, errstat)
        if (errstat /= 0) call exit(1)

      !reverse vertical order so that it is ordered TOA->BOA
      !pres_z in hPa on grid levels
      !temp_z in K between grid levels
      do iz=1,nz
         kk = nz - iz + 1
         pp(iz) = pres_z(kk+1)
         tt(iz) = temp_z(kk)
         geos_Pressure(ix,it,iz) = pp(iz)
         geos_Temperature(ix,it,iz) = tt(iz)
      enddo !iz
      pp(nz+1) = psurf

     ! get wind speed at 2m
     call clim_val_interp(u2m_cst, cpt, thishour, thislon, thislat,u2m, errstat)
     call clim_val_interp(v2m_cst, cpt, thishour, thislon, thislat,v2m, errstat)
     windspeed2m(ix,it) = sqrt(u2m(1)*u2m(1) + v2m(1)*v2m(1))

!hqw debug
      if ((ix .eq. ixdebug) .AND. (it .eq. itdebug)) then
         write(*,*) 'ix,it,u2m,v2m,windspeed2m:'
         write(*,*) ix, it, u2m, v2m, windspeed2m(ix,it)
         write(*,*) 'P-T profile:'
         do iz = 1, nz
            write(*,*) iz, pp(iz), tt(iz)
         enddo
         write(*,*) nz+1, pp(nz+1)
         write(*,*)'psurf=',psurf
      endif

    enddo ! ix
    !write(*,*) it, thisyear, thismonth, thisday, thishour
  enddo ! it

  deallocate(pres_z, temp_z)

  end subroutine read_geoscf

end module m_read_input_clim
