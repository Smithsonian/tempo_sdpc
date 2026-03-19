MODULE  m_get_met_tempo
  USE OMSAO_precision_module, ONLY: i4,r4,sp, dp
  USE OMSAO_parameters_module, ONLY: maxchlen, accgrav
  USE OMSAO_variables_module, ONLY: l2_met_filenames, the_lon, the_lat, the_time, the_surfalt, &
    lat_min, lat_max, lon_min, lon_max, time_min, time_max
  USE OMI_LUN_SET, ONLY:num_met_luns
  TYPE :: met_type
     LOGICAL :: do_debug=.false., adj_spres = .true.
     INTEGER :: np
     REAL (kind=dp) :: ptrop, psurf, z0, phis, ztrop, zpbl2, ppbl2, zpbl, ppbl
     REAL (kind=dp), DIMENSION(:), allocatable :: ts, ps
  END TYPE 
  TYPE (met_type) :: thismet
  CONTAINS
  SUBROUTINE get_met_tempo (errstat)
  USE, intrinsic :: iso_c_binding, only: c_ptr, c_char,c_null_char, c_null_ptr, c_associated
  use tell_module
  USE met_module, only : synth_met_type, open_synth_met_data, read_synth_met_data
  USE m_ezspline_interpolation, ONLY: bspline
  use OMSAO_variables_module, only: apriori_source
  use clim_module
  IMPLICIT NONE
  !===========================
  ! input/output 
  !===========================
  INTEGER, INTENT(OUT) :: errstat
  !===========================
  ! local variables
  !===========================
  INTEGER :: np, ilay
  INTEGER, PARAMETER :: nlecm = 67 !-changed by junsung (JAN 2024)
  real(kind=dp), dimension(:), allocatable :: pprof, tprof, pmid
  real(kind=sp), dimension(:), allocatable :: pres_z, temp_z
  real(kind=sp) :: phi(1)!, zpbl(1) !-added by junsung (JAN 2024)
  real(kind=dp) :: phis !-added by junsung (JAN 2024)
  real(kind=sp) :: lon_f, lat_f, psurf, ptrop
  real(kind=dp) :: local_srf
  real(kind=dp) :: tempo_pres, model_pres, model_stemp, model_srf !added by junsung (JAN 2024)

  !-These temperatures calulcated based on the pressure (from GEMS-Chem vertical grid) values
  !-modified by junsung to change the number of layers (JAN 2024)
  real(kind=dp), dimension(0:nlecm), parameter :: clim = &
       (/288.20, 287.40, 286.55, 285.68, 284.81, 283.93, 283.03, 282.13, 281.21, 280.28, &
         279.34, 278.39, 277.42, 276.11, 274.43, 272.72, 270.95, 269.14, 267.28, 264.38, &
         261.34, 258.15, 254.78, 251.20, 247.40, 243.32, 238.92, 234.12, 227.04, 219.92, &
         218.20, 218.20, 218.20, 218.20, 218.20, 218.20, 218.20, 218.20, 218.20, 219.12, &
         220.20, 221.31, 222.44, 223.60, 224.80, 226.02, 227.29, 228.58, 229.92, 233.28, &
         237.31, 241.48, 245.91, 250.28, 254.91, 259.68, 264.60, 269.66, 272.20, 272.20, &
         269.99, 264.44, 258.78, 252.98, 247.16, 241.29, 235.39, 229.46/)

  type (synth_met_type), SAVE :: smt
  TYPE (clim_pres_type), SAVE :: cpt
  TYPE (clim_pres_bounds_type), SAVE :: bounds
  TYPE (clim_val_type), SAVE :: cst, cst_phis!, cst_zpbl !-added cst_phis by junsung (JAN 2024)
  INTEGER, SAVE :: nl0
  integer :: year(2), month(2), day(2), k, j
  real (kind=dp) :: hour
  real (kind=sp) :: hour_f
  logical, SAVE :: have_synthetic_met_data
  logical, SAVE :: first = .true.
  logical :: have_forecast
  
  ! Initialize dataset
  IF (first) THEN 
    !print * , l2_met_filenames(1)
    if (0 /= index (l2_met_filenames(1), '.nc', .true.)) then
      have_synthetic_met_data = .true.
      call open_synth_met_data (smt, trim(l2_met_filenames(1)),errstat)
      !print * , l2_met_filenames(1)
      if (errstat /= 0) then
        call tell_error (tell_runtime_error, "error opening synthetic met data", errstat)
        STOP 1
      endif
    else
      have_synthetic_met_data = .false.
      if (time_max - time_min > 86400.0) then
        call tell_error (tell_runtime_error, "get_met_tempo: granule duration exceeds 24 hours", errstat)
        return
      endif

      call tio_f_taix_time_to_utc_caldate(time_min, year(1), month(1), day(1), hour)
      bounds % hour_beg = real (hour, kind=r4)
      call tio_f_taix_time_to_utc_caldate(time_max, year(2), month(2), day(2), hour)
      bounds % hour_end = real (hour, kind=r4)
      bounds % lon_min = real(lon_min,kind=r4)
      bounds % lon_max = real(lon_max,kind=r4)
      bounds % lat_min = real(lat_min,kind=r4)
      bounds % lat_max = real(lat_max,kind=r4)

      !@ set bounds
      call clim_pres_init (cpt, year(1), month(1), day(1), bounds, errstat)
      call clim_query_nz (nl0, errstat)
      if (errstat /= 0) THEN
        call tell_error (tell_runtime_error, "get_met_tempo: errors in clim_pres_init", errstat)
        return
      endif
      call clim_query_apriori_source (cpt, have_forecast, errstat)
      if (have_forecast) then
        apriori_source = 'GEOSCF:forecast'
      else
        apriori_source = 'GEOSCF:climatology'
      endif
      !print *, nl0
      !print *, apriori_source

      call clim_val_init (cst, cpt, 'T'//c_null_char, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, "get_met_tempo: initializing air temperature", errstat)
        return
      endif

      !-added below 6 line by junsung (JAN 2024)-surface geopotential height (m^2 / s^2)
      !-and define the error state
      call clim_val_init (cst_phis, cpt, 'PHIS', errstat, single_layer=.true.)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, "get_met_tempo: initializing PHIS", errstat) !-added by junsung (JAN 2024)
        return
      endif

      ! xl: add below to get GEOS-CF planetary boundary layer height
      !call clim_val_init (cst_zpbl, cpt, 'ZPBL', errstat, single_layer=.true.)
      !if (errstat /= 0) then
      !  call tell_error (tell_io_read_error, "get_met_tempo: initializing ZPBL", errstat) 
      !  return
      !endif

    endif
    first = .false.
   endif

   call tio_f_taix_time_to_utc_caldate(the_time, year(1), month(1), day(1), hour)
   hour_f = real(hour, kind=r4)

   !xl, 1/2/2022 better to use the original pressure/T profiles,
   ! so most of the following block commented
   ! define user vertical grids 
   np = nlecm 
   thismet%np = np
   if (allocated(thismet%ts)) deallocate(thismet%ts, thismet%ps)
   allocate(thismet%ts(0:thismet%np),thismet%ps(0:thismet%np))
   allocate(tprof(0:np),pprof(0:np))
  
   !Interpolation dataset into current pixel
   lon_f = real(the_lon, kind=sp) ; lat_f = real(the_lat, kind=sp)
   IF (have_synthetic_met_data) then
    CALL read_synth_met_data(smt, lon_f, lat_f, ptrop, errstat, &
         pprof = pprof, tprof=tprof)
    thismet%psurf = pprof(0)
    thismet%ptrop = real(ptrop, kind=dp)
   ELSE
     allocate (pres_z(nl0+1), pmid(nl0), temp_z(nl0))
     call clim_pres (cpt, hour_f, lon_f, lat_f, pres_z, errstat, &
                     p_surf = psurf, p_trop = ptrop)
     call clim_val_interp (cst, cpt, hour_f, lon_f, lat_f, temp_z, errstat)
     pmid(1:nl0) = real(0.5 * (pres_z(1:nl0) + pres_z(2:nl0+1)), kind=dp)
     
     !-added below 3 lines to define the error state by junsung (JAN 2024)
     if (errstat /= 0) then
        call tell_error (tell_runtime_error, "get_met_tempo: calculating pressure", errstat)
     end if

     !-added below 5 lines to read geopotantial height and define error state by
     !junsung (JAN 2024)
     call clim_val_interp (cst_phis, cpt, hour_f, lon_f, lat_f, phi, errstat)
     phis = real(phi(1), kind=dp)
     if (errstat /= 0) then
       call tell_error (tell_runtime_error, "get_met_tempo: calculating surface geopotential height", errstat)
     end if     

     ! xl: added the following to read planetary boundary layer height
     !call clim_val_interp (cst_zpbl, cpt, hour_f, lon_f, lat_f, zpbl, errstat)
     !if (errstat /= 0) then
     !  call tell_error (tell_runtime_error, "get_met_tempo: calculating boundary layer height", errstat)
     !end if  

     ! Derive surface temperature from layer 1 & 2 temperature
     tprof(0) = temp_z(1) + (temp_z(1) - temp_z(2)) &
                / (LOG(pmid(1)) - LOG(pmid(2))) * (LOG(pres_z(1)) - LOG(pmid(1)))

     DO j = 1, np
        tprof(j) = 2.0 * temp_z(j) - tprof(j-1)
     ENDDO

     pprof(0:np) = pres_z(1:np+1)
     thismet%psurf = real(psurf, kind=dp)
     thismet%ptrop = real(ptrop, kind=dp)
     thismet%phis =  phis
     !thismet%zpbl = zpbl(1) / 1000.0_dp ! convert to km

     !print *, thismet%psurf, thismet%ptrop, phis
     !print *, pres_z(1:nl0+1)
     !print *, pmid(1:nl0)
     !print *, temp_z(1:nl0)
     !print *, pprof(0:np)
     !print *, tprof(0:np)

     !! May need to adjust pprof(:) near the surface
     !k = 0
     !j = 1
     !do while (k <= np .and. j <= nl0 .and. pprof(k) > pmid(j))
     !  pprof(k) = pmid(j)
     !  k = k + 1
     !  j = j + 1
     !enddo
     !
     !call BSPLINE(pmid(1:nl0), real(temp_z(1:nl0),kind=dp), nl0, &
     !             pprof(0:np), tprof(0:np), np+1, errstat)
     deallocate (pres_z, pmid, temp_z)
  ENDIF

  thismet%ps(0:thismet%np)=pprof
  thismet%ts(0:thismet%np)=tprof
  deallocate(tprof, pprof)
  
  ! heanding the pressure grids of out of range
  do ilay = 0, np
     if (thismet%ts(ilay) < 150. .or. thismet%ts(ilay) > 350.) THEN
        thismet%ts(ilay) = clim(ilay)
     endif
  enddo

  ! surface pressure
  local_srf = REAL(the_surfalt, kind=dp)

  thismet%z0 = the_surfalt
  IF (thismet%adj_spres) then
     model_pres = REAL(thismet%psurf, kind=dp)
     model_srf = REAL(phis / accgrav / 1000.0_dp, kind=dp)
     model_stemp = thismet%ts(0) 
     call adjust_surface_pressure(local_srf, model_srf, model_pres, model_stemp, tempo_pres)
     !print *, model_srf, model_pres, local_srf, tempo_pres
     
     !tempo_pres = model_pres * 10.0_dp ** ((model_srf - local_srf)/16.0_dp)
     !print *, model_srf, model_pres, local_srf, tempo_pres
     
     IF (tempo_pres > model_pres) THEN
        thismet%psurf = REAL(tempo_pres, kind=sp)
        thismet%ps(0) = tempo_pres
        
        ! Adjust surface temperature
        thismet%ts(0) = model_stemp + (model_stemp - thismet%ts(1)) &
             / (LOG(model_pres) - LOG(thismet%ps(1))) * (LOG(tempo_pres) - LOG(model_pres))
     ELSE IF(tempo_pres < model_pres) THEN
        DO j = 1, np
           IF (tempo_pres > thismet%ps(j)) EXIT
        ENDDO
        j = j - 1
        thismet%ts(j) = thismet%ts(j+1) + (thismet%ts(j) - thismet%ts(j+1)) &
             / (LOG(thismet%ps(j)) - LOG(thismet%ps(j+1))) &
             * (LOG(tempo_pres) - LOG(thismet%ps(j+1)))
        thismet%ps(j) = tempo_pres
        thismet%psurf = tempo_pres
        np = np -j
        thismet%np = np
        thismet%ps(0:np) = thismet%ps(j:j+np)
        thismet%ts(0:np) = thismet%ts(j:j+np)
     ENDIF
  ENDIF
  
  !print *, np
  !print *, thismet%ps(0:np)
  !print *, thismet%ts(0:np)
  !print *, thismet%psurf, thismet%ptrop, thismet%z0!, thismet%zpbl

  !local_srf = 1013.0_dp * (10.0_dp ** (local_srf / (-16.0_dp)))
  !thismet%psurf = real(local_srf, kind=sp)
  !thismet%z0 = real(the_surfalt, kind=sp)
  
  ! Debuging
  IF (thismet%do_debug) THEN 
     write(*,'(f8.1)')thismet%psurf
     write(*,'(f8.1)')thismet%ptrop
     write(*,'(f8.1)')thismet%ts
  ENDIF
  RETURN
END SUBROUTINE get_met_tempo


  ! added adjust_surface_pressure subroutine to adjust from the geopotential
  ! height to surface pressure (model --> tempo)
  ! this subroutine from the TEMPO tracegas algorithm
  SUBROUTINE adjust_surface_pressure (tempo_height, model_height, model_pressure, model_temperature, pressure)
   ! Hypsometric equation. Follows equation 3 of Boersma et al., 2011
   ! https://amt.copernicus.org/articles/4/1905/2011/
   ! Given the terrain height for a satellite pixel, model terrain height, model
   ! surface temperature and
   ! model surface pressure compute adjusted pressure
   ! P_pix = P_mod ( T_mod / (T_mod + lr (H_mod - H_pix) ) )^(-g / R / lr *
   ! 1000.0)
   ! Heights have to be in units of km.
   implicit none

   real (kind=dp), intent(in) :: tempo_height, model_height, model_pressure
   real (kind=dp), intent(in) :: model_temperature
   real (kind=dp), intent(out) :: pressure ! adjusted pressure

   real (kind=dp), parameter :: lr = 6.5 ! lapse rate (K km^-1)
   real (kind=dp), parameter :: g = 9.80665  ! Gravitational constant (m s^-2)
   real (kind=dp), parameter :: R = 287.058  ! Gas constant for dry air (J kg^-1 K^-1)
  
   real (kind=dp) :: a, b
 
   ! Hypsometric equation
   a = - (g / R / lr * 1000.0)
   b = lr * (model_height - tempo_height)
   pressure = REAL(model_pressure * ( model_temperature / (model_temperature + b) )**a, KIND=dp)
!print*, '11'
!print*, a
!print*, b
!print*, model_pressure
!print*, model_temperature
!print*, ( model_temperature / (model_temperature + b) )
!print*, ( model_temperature / (model_temperature + b) )**a
  END SUBROUTINE adjust_surface_pressure

END MODULE m_get_met_tempo
