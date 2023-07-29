! Read GEOS-CF input climatology
module m_read_input_clim
   use clim_module
   use tell_module

   public read_geoscf

   private prepare_geoscf, adjust_surface_pressure, psfc_topo_adjust
   real(kind=8), parameter :: r8_missval = -1.0d+30

contains

   subroutine prepare_geoscf(cpt,temp_cst,u2m_cst,v2m_cst,q_cst,&
              phis_cst,lon_min,lon_max,lat_min,lat_max, errstat)

   use m_vars, only: geos_np, nlayers, gmetadata, &
              rad_latitude, rad_longitude, rad_Time

  implicit none

  type (clim_pres_type), intent(out) :: cpt
  type (clim_val_type), intent(out) :: temp_cst,u2m_cst,v2m_cst
  type (clim_val_type), intent(out) :: phis_cst, q_cst
  real(kind=4), intent(out) :: lon_min,lon_max,lat_min,lat_max
  integer, intent(inout) :: errstat

  type (clim_pres_bounds_type) :: bounds

  integer :: year(2), month(2), day(2)
  integer :: thisyear, thismonth, thisday

  integer :: nz
  logical :: have_forecast

  real (kind=8) :: t_beg, t_end, hour_beg, hour_end

  if (errstat /= 0) return

! write(*,*) 'Preapre to read GEOS-CF for orbit'
! rad_Time was read by m_read_input_tio
  t_beg = minval(rad_Time, rad_Time /= r8_missval)
  t_end = maxval(rad_Time, rad_Time /= r8_missval)

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

! The run-time environment should have specified config file location

  lon_min = -180.0 !initialize
  lon_max = 180.0
  lat_min = -90.
  lat_max = 90.

  lon_min = minval(rad_longitude, rad_longitude .gt. -180.0001)
  bounds%lon_min = lon_min
  lon_max = maxval(rad_longitude, rad_longitude .lt. 180.0001)
  bounds%lon_max = lon_max
  lat_min = minval(rad_latitude, rad_latitude .gt. -90.0001)
  bounds%lat_min = lat_min
  lat_max = maxval(rad_latitude, rad_latitude .lt. 90.0001)
  bounds%lat_max = lat_max
  !hqw hour_start & hour_end = gmeta(hour+minute/60.+seconds/3600.)
  write(*,*)'   bounds year,month,day',year(1),month(1),day(1)
  write(*,*)'   bounds hour_beg,end = ',bounds%hour_beg, bounds%hour_end
  write(*,*)'   bounds lon:',bounds%lon_min,bounds%lon_max
  write(*,*)'   bounds lat:',bounds%lat_min,bounds%lat_max

! check consistency with gmatadata
!there is slight difference with bounds for lon_min,lon_max,lat_min,lax_max
  write(*,*) 'gmeta year,month,day',&
     gmetadata%granule_year,gmetadata%granule_month, gmetadata%granule_day
  write(*,*) 'gmeta begin hour,minute,seconds=',&
     gmetadata%granule_hour_start,gmetadata%granule_minute_start,gmetadata%granule_seconds_start
  write(*,*) 'gmeta end hour,minute,seconds=',&
     gmetadata%granule_hour_end, gmetadata%granule_minute_end,gmetadata%granule_seconds_end
  write(*,*)'gmeta lon:',gmetadata%geospatial_lon_min,gmetadata%geospatial_lon_max
  write(*,*)'gmeta lat:',gmetadata%geospatial_lat_min,gmetadata%geospatial_lat_max

  call clim_pres_init(cpt,thisyear,thismonth,thisday,bounds,errstat)

  call clim_query_apriori_source (cpt, have_forecast, errstat)
  call clim_query_nz (nz, errstat)
  if (errstat /= 0) return

  ! initialize the global geos_np and nlayers
  ! to the number of layers in the GEOSCF file.
  geos_np = nz
  nlayers = nz

  if (have_forecast) then
    gmetadata % apriori_source = 'GEOSCF:forecast'
  else
    gmetadata % apriori_source = 'GEOSCF:climatology'
  endif

  ! temperature profile
  call clim_val_init (temp_cst, cpt, 'T', errstat)
  ! specific humidity profile
  call clim_val_init (q_cst, cpt, 'Q',errstat)
  if (errstat /= 0) call exit(1)

  call clim_val_init (u2m_cst, cpt, 'U2M', errstat, single_layer=.true.)
  call clim_val_init (v2m_cst, cpt, 'V2M', errstat, single_layer=.true.)
  call clim_val_init (phis_cst,cpt, 'PHIS',errstat, single_layer=.true.)

  if (errstat /= 0) call exit(1)

  end subroutine prepare_geoscf

!
!------------
!
   subroutine read_geoscf (errstat)

   use m_vars, only: rad_NumTimes,rad_nXtrack, geos_np, &
       rad_latitude, rad_longitude, rad_Time, out_TerrainHeight,&
       l2_TerrainPressure, windspeed2m, phisurf,&
       ixdebug, itdebug, run_mode, lun_debug_clim, lun_debug_psfc

   use m_vars, only: geos_Temperature, geos_Pressure, geos_Q

   implicit none

   integer, intent(inout) :: errstat

   type(clim_pres_type) :: cpt
   type(clim_val_type) :: temp_cst, u2m_cst, v2m_cst
   type(clim_val_type) :: q_cst, phis_cst

   real(kind=8) :: hour
   real(kind=4) :: thislon, thislat, thishour
   real(kind=4) :: lon_min,lon_max,lat_min,lat_max
   real(kind=4), dimension(:), allocatable :: pp(:), tt(:), qq(:)

   real (kind=4) :: psurf, ptrop
   real (kind=4), dimension(1) :: u2m, v2m, thisphis
   real (kind=4), dimension(:), allocatable :: pres_z, temp_z, q_z
   
   real (kind=4), parameter :: g_grav = 9.80665

   real (kind=4) :: model_tsurf, pixel_height, model_height, &
                    model_qsurf, adj_pressure

   real (kind=4), parameter :: psfc_recordhigh = 1085.0

   integer :: ix, it, iz, nx, nt, kk, nz
   integer :: thisyear, thismonth, thisday

   if (errstat /= 0) return

  ! initialize
  nx = rad_nXtrack
  nt = rad_NumTimes

  lon_min=-180.0
  lon_max=180.0
  lat_min=-90.0
  lat_max=90.0

  call prepare_geoscf(cpt,temp_cst,u2m_cst,v2m_cst,q_cst,&
       phis_cst,lon_min,lon_max,lat_min,lat_max, errstat)
  if (errstat /= 0) return

  nz = geos_np
  write(*,*)'read_geoscf:nXtrack,nTimes,nz=',nx,nt,nz

  !allocate & init geos_Temperature, geos_Q, geos_Pressure
  allocate(geos_Temperature(nx,nt,nz), stat=errstat)
  allocate(geos_Q(nx,nt,nz), stat=errstat)
  allocate(geos_Pressure(nx,nt,nz+1), stat=errstat)
  geos_Temperature = -999.
  geos_Q = 0.
  geos_Pressure = -999.

  !allocate & init windspeed2m
  allocate(windspeed2m(nx,nt), stat=errstat)
  ! initialization will be used by GMI, overwritten by GEOS-CF 
  windspeed2m = 5. 
  
  !allocate & init surface geopotential height
  allocate(phisurf(nx,nt), stat=errstat)
  phisurf = -999.

  !allocate temporary variables
  allocate(temp_z(nz), pres_z(nz+1), stat=errstat)
  allocate(q_z(nz), stat=errstat)
  allocate(tt(nz), pp(nz+1), stat=errstat)
  allocate(qq(nz), stat=errstat)

  ! loop through pixels
   do it = 1, nt
     if (rad_Time(it) == r8_missval) cycle

     call tio_f_taix_time_to_utc_caldate(rad_Time(it), &
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

      ! get temperature at lon lat
        call clim_val_interp (temp_cst, cpt, thishour, thislon, thislat, &
                 temp_z, errstat)
        if (errstat /= 0) call exit(1)

      ! get specific humidity q at lon lat
        call clim_val_interp (q_cst, cpt, thishour, thislon, thislat, &
                 q_z, errstat)

      !reverse vertical order so that it is ordered TOA->BOA
      !pres_z in hPa on grid levels
      !temp_z in K between grid levels
      pp(nz+1) = psurf
      do iz=1,nz
         kk = nz - iz + 1
         pp(iz) = pres_z(kk+1)
         tt(iz) = temp_z(kk)
         qq(iz) = q_z(kk)
         geos_Pressure(ix,it,iz) = pp(iz)
         geos_Temperature(ix,it,iz) = tt(iz)
         geos_Q(ix,it,iz) = qq(iz)
      enddo !iz

      ! assign surface pressure (hPa)
      geos_Pressure(ix,it,nz+1) = psurf

     ! get wind speed at 2m for GLER over water surface
     call clim_val_interp(u2m_cst, cpt, thishour, thislon, thislat,u2m, errstat)
     call clim_val_interp(v2m_cst, cpt, thishour, thislon, thislat,v2m, errstat)
     windspeed2m(ix,it) = sqrt(u2m(1)*u2m(1) + v2m(1)*v2m(1))

     ! get surface geopotential height for psfc topo correction
     call clim_val_interp(phis_cst,cpt,thishour,thislon,thislat,thisphis,errstat)
     phisurf(ix,it) = thisphis(1)

      ! l2_TerrainPressure is topography adjusted surface pressure
      ! can use geos_Pressure as an approximation, 
      ! when forecast or climatology is used for TEMPO
      ! it is better to adjust when reanalysis meteorology is used
      ! calculation should use l2_TerrainPressure 
      !  l2_TerrainPressure(ix, it) = psurf
      pixel_height = out_TerrainHeight(ix,it) ! meters

      !  do adjustment only if pixel_height is reasonable
      if (pixel_height .gt. -1000.) then
         model_height = thisphis(1) / g_grav ! meters
         ! tt,qq from TOA to BOA, thus, level closet to surface is nz
         model_tsurf = tt(nz)
         model_qsurf = qq(nz)
         ! adjust_surface_pressure is used in L1_trace_gas, assumes dry air
         ! adj_pressure has the same unit as psurf
         !call adjust_surface_pressure(pixel_height,model_height, &
         !     psurf,model_tsurf,adj_pressure,errstat)
         ! psfc_topo_adjust considers humidity
         call psfc_topo_adjust(pixel_height,model_height,psurf, &
            model_tsurf,model_qsurf,adj_pressure,errstat)
         if (errstat /= 0) adj_pressure = psurf
         ! safeguard surface pressure (wikipedia record is ~1085hPa)
         if (adj_pressure .GT. psfc_recordhigh) adj_pressure = psfc_recordhigh
         l2_TerrainPressure(ix,it) = adj_pressure
      else ! unrealistic pixel_height, typically a fill value 
         ! no adjustment in this case
         l2_TerrainPressure(ix,it) = psurf
      endif
     
    enddo ! ix
  enddo ! it

!----------------------
! debug
   if ((trim(run_mode) .eq. 'development') .and. &
       (ixdebug .ge. 0) .and. (itdebug .ge. 0)) then
      ix = ixdebug
      it = itdebug
      write(*,*) ' writing debug_geos.txt'
      open(unit=lun_debug_clim, file='debug_geos.txt')
         write(lun_debug_clim,*) 'ix,it=',ix,it
         write(lun_debug_clim,*) 'windspeed2m, phis ='
         write(lun_debug_clim,*) windspeed2m(ix,it),phisurf(ix,it)
         write(lun_debug_clim,*) 'iz    pp      tt      qq'
         do iz = 1, nz
            write(lun_debug_clim,*) iz, geos_Pressure(ix,it,iz), &
                 geos_Temperature(ix,it,iz), geos_Q(ix,it,iz)
         enddo
         write(lun_debug_clim,*) nz+1, geos_Pressure(ix,it,nz+1)
       close(lun_debug_clim)

      write(*,*) '   writing debug_psfc.txt'
      open(unit=lun_debug_psfc, file='debug_psfc.txt')
      write(lun_debug_psfc,*) 'itdebug = ',itdebug
      write(lun_debug_psfc,*) &
         'model_height   pixel_height   model_psfc   pixel_psfc'
      do ix = 1, nx
          model_height = phisurf(ix,it) / g_grav
          pixel_height = out_TerrainHeight(ix,it)
          psurf = geos_Pressure(ix,it,nz+1)
          adj_pressure = l2_TerrainPressure(ix,it)
          write(lun_debug_psfc,*) model_height,pixel_height,psurf,adj_pressure
      enddo
      close(lun_debug_psfc)
   endif

!------
  deallocate(pres_z, temp_z, q_z)
  deallocate(pp, tt, qq)

  end subroutine read_geoscf

!!!!!!!!!!!!
! hypsometric psfc topography adjustment
!!!!!!!!!!!!
  subroutine psfc_topo_adjust(pixel_height,model_height,model_pressure, &
             model_temperature, model_q, adj_pressure, errstat)
  ! hypsometric equation AMS glossary of meteorology:
  ! P2 = P1*exp((z1-z2)/(Rd/g*Tv))
  ! z1 & z1 in meters, Tv in K
  ! Tv = T * (1 + 0.608*q) ! Wikipedia
  ! q = specific humidity in unit of [g H2O/ g dry air] = [kg/kg]

  real(kind=4),intent(in):: pixel_height,model_height,model_pressure
  real(kind=4),intent(in):: model_temperature, model_q
  real(kind=4), intent(out):: adj_pressure
  integer(kind=4),intent(inout):: errstat
  
  real(kind=4):: Rdovg = 287.058 / 9.80665
  
  real(kind=4):: virtual_temperature, abc

  virtual_temperature = model_temperature*(1.+0.608*model_q)

  abc = (model_height - pixel_height)/virtual_temperature/Rdovg

  adj_pressure = model_pressure*exp(abc)

  end subroutine psfc_topo_adjust

!!!!!!!!!!!!
! borrow adjust_surface_pressure from L1_trace_gas/OMSAO_wramf_module
! this keeps consistency with trace gas AMF calculation, but for dry air
! slight modification is made for model_heigt to physical height (m)
!        instead of the original geopotential height (m^2/s^2)
!!!!!!!!!!!!

  subroutine adjust_surface_pressure (pixel_height, model_height, &
             model_pressure, model_temperature, adj_pressure, errstat)
   ! Hypsometric equation. Follows equation 3 of Boersma et al., 2011
   ! https://amt.copernicus.org/articles/4/1905/2011/
   ! Given the terrain height for a satellite pixel, model terrain height, model
   ! surface temperature and
   ! model surface pressure compute adjusted pressure
   ! P_pix = P_mod ( T_mod / (T_mod + lr (H_mod - H_pix) ) )^(-g / R / lr * 1000.0)
   ! Heights are in units of km for this equation
   ! input pixel_height should be in meters
   ! input model_height should also be in meters

   implicit none

   real (kind=4), intent(in) :: pixel_height, model_height, model_pressure, model_temperature
   integer (kind=4), intent(inout) :: errstat
   real (kind=4), intent(out) :: adj_pressure ! adjusted pressure

   real (kind=4), parameter :: lr = 6.5 ! lapse rate (K km^-1)
   real (kind=4), parameter :: g = 9.80665  ! Gravitational constant (m s^-2)
   real (kind=4), parameter :: R = 287.058  ! Gas constant for dry air (J kg^-1 K^-1)

   real (kind=4) :: a, b, mh, sh

   if (errstat /= 0) return

   ! model height from meter to km
   mh = model_height / 1000.0
   ! Satellite pixel height from meter to km
   sh = pixel_height / 1000.0

   ! Hypsometric equation
   a = - (g / R / lr * 1000.0)
   b = lr * (mh - sh)
   adj_pressure = model_pressure*(model_temperature/(model_temperature + b))**a

  end subroutine adjust_surface_pressure

end module m_read_input_clim
