MODULE  m_get_met_tempo
  USE OMSAO_precision_module, ONLY: i4,r4,sp, dp
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module, ONLY: l2_met_filenames, the_lon, the_lat, the_time, the_surfalt, &
    lat_min, lat_max, lon_min, lon_max, time_min, time_max
  USE OMI_LUN_SET, ONLY:num_met_luns
  TYPE :: met_type
     LOGICAL :: do_debug=.false.
     INTEGER :: np
     REAL (kind=r4) ::ptrop, psurf, z0
     REAL (kind=dp), DIMENSION(:), allocatable :: ts, ps
  END TYPE 
  TYPE (met_type) :: thismet
  CONTAINS
  SUBROUTINE get_met_tempo (errstat)
  USE, intrinsic :: iso_c_binding, only: c_ptr, c_char,c_null_char, c_null_ptr, c_associated
  use tell_module
  USE met_module, only : synth_met_type, open_synth_met_data, read_synth_met_data
  USE m_ezspline_interpolation, ONLY: bspline
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
  INTEGER, PARAMETER :: nlecm = 30
  real(kind=dp), dimension(:), allocatable :: pprof, tprof, pmid
  real(kind=sp), dimension(:), allocatable :: pres_z, temp_z
  real(kind=sp) :: lon_f, lat_f
  real(kind=dp) :: local_srf
  real(kind=dp), dimension(0:nlecm), parameter :: clim = & 
       (/288.2, 287.5,286.1, 284.7, 283.2, 281.8, &
        278.7, 275.5, 272.2, 268.6, 264.8, 260.8, & 
       256.6, 251.9, 246.9, 241.4,  235.4, 228.6, & 
       220.5,216.5, 216.7, 216.7,   216.7, 217.2,   &
       220.5, 223.1, 227.7, 232.4,  239.3,249.5, 257.9/)
  type (synth_met_type), SAVE :: smt
  TYPE (clim_pres_type), SAVE :: cpt
  TYPE (clim_pres_bounds_type), SAVE :: bounds
  TYPE (clim_val_type), SAVE :: cst
  INTEGER, SAVE :: nl0
  integer :: year, month, day, k, j
  real (kind=dp) :: hour
  real (kind=sp) :: hour_f
  logical, SAVE :: have_synthetic_met_data
  logical, SAVE :: first = .true.
  
  ! Initialize dataset
  IF (first) THEN 
    if (0 /= index (l2_met_filenames(1), '.nc', .true.)) then
      have_synthetic_met_data = .true.
      call open_synth_met_data (smt, trim(l2_met_filenames(1)),errstat)
      print * , l2_met_filenames(1)
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

      call tio_f_taix_time_to_utc_caldate(time_min, year, month, day,hour)
      bounds % hour_beg = real (hour, kind=r4)
      call tio_f_taix_time_to_utc_caldate(time_max, year, month, day,hour)
      bounds % hour_end = real (hour, kind=r4)
      bounds % lon_min = real(lon_min,kind=r4)
      bounds % lon_max = real(lon_max,kind=r4)
      bounds % lat_min = real(lat_min,kind=r4)
      bounds % lat_max = real(lat_max,kind=r4)

      !@ set bounds
      call clim_pres_init (cpt, year, month, day, bounds, errstat)
      call clim_query_nz (nl0, errstat)
      if (errstat /= 0) THEN
        call tell_error (tell_runtime_error, "get_met_tempo: errors in clim_pres_init", errstat)
        return
      endif

      call clim_val_init (cst, cpt, 'T'//c_null_char, errstat)
      if (errstat /= 0) then
        call tell_error (tell_io_read_error, "get_met_tempo: initializing air temperature", errstat)
        return
      endif
    endif
    first = .false.
   endif

   call tio_f_taix_time_to_utc_caldate(the_time, year, month, day, hour)
   hour_f = real(hour, kind=r4)

   !xl, 1/2/2022 better to use the original pressure/T profiles,
   ! so most of the following block commented
   ! define user vertical grids 
   np = nlecm 
   thismet%np = np
   if (allocated(thismet%ts)) deallocate(thismet%ts, thismet%ps)
   allocate(thismet%ts(0:thismet%np),thismet%ps(0:thismet%np))
   allocate(tprof(0:np),pprof(0:np))
   pprof(0:np) = (/1013.25, 1000., 975., 950., 925., 900., 850., &
                  800., 750., 700., 650., 600., 550., 500., 450., &
                  400., 350., 300.,250.,200., 150., 100., 70., 50., &
                  30., 20., 10., 7., 5., 3., 2./)

   !Interpolation dataset into current pixel
   lon_f = real(the_lon, kind=sp) ; lat_f = real(the_lat, kind=sp)
   IF (have_synthetic_met_data) then
    CALL read_synth_met_data(smt, lon_f, lat_f, thismet%ptrop, errstat, &
                             pprof = pprof, tprof=tprof)
   ELSE
     allocate (pres_z(nl0+1), pmid(nl0), temp_z(nl0))
     call clim_pres (cpt, hour_f, lon_f, lat_f, pres_z, errstat, &
                     p_surf = thismet%psurf, p_trop = thismet%ptrop)
     call clim_val_interp (cst, cpt, hour_f, lon_f, lat_f, temp_z, errstat)
     pmid(1:nl0) = real(0.5 * (pres_z(1:nl0) + pres_z(2:nl0+1)), kind=dp)
     ! May need to adjust pprof(:) near the surface
     k = 0
     j = 1
     do while (k <= np .and. j <= nl0 .and. pprof(k) > pmid(j))
       pprof(k) = pmid(j)
       k = k + 1
       j = j + 1
     enddo
    
     call BSPLINE(pmid(1:nl0), real(temp_z(1:nl0),kind=dp), nl0, &
                  pprof(0:np), tprof(0:np), np+1, errstat)
     deallocate (pres_z, pmid, temp_z)
     !write(*,'(a,5(2x,f10.4))')'get_met_tempo: hour,lon,lat=',hour_f,lon_f,lat_f
     !write(*,'(a,2x,f10.4,2x,a,f10.4)')'get_met_tempo: psurf=',thismet%psurf,' ptrop=',thismet%ptrop
     !write(*,'(a,5(2x,e12.4))')'tprof',tprof(0:np)
   ENDIF
   thismet%ps(0:thismet%np)=pprof
   thismet%ts(0:thismet%np)=tprof
   deallocate(tprof, pprof)
   IF (errstat /= 0) THEN 
      WRITE(*,*) 'Errors in get_met_tempo' ; STOP 1
   ENDIF
   
   ! heanding the pressure grids of out of range
   do ilay = 0, np
     if (thismet%ts(ilay) < 0 .or. thismet%ts(ilay) > 1000) THEN
         thismet%ts(ilay) = clim(ilay)
     endif
   enddo

   ! surface pressure 

   local_srf = REAL(the_surfalt, kind=dp)
   local_srf = 1013.0_dp * (10.0_dp ** (local_srf / 1000.0_dp / (-16.0_dp)))
   thismet%psurf = real(local_srf, kind=sp)
   thismet%z0 = real(the_surfalt, kind=sp)

   ! Debuging
   IF (thismet%do_debug) THEN 
    write(*,'(f8.1)')thismet%psurf
    write(*,'(f8.1)')thismet%ptrop
    write(*,'(f8.1)')thismet%ts
   ENDIF
   RETURN
  END SUBROUTINE get_met_tempo
END MODULE m_get_met_tempo
