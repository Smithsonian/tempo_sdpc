!Read GLER
module m_read_input_gler
   use gler_module
   use tell_module

contains
   subroutine read_gler (errstat)

   use m_vars,only: BRDF_SurfaceReflectivity466,&
       BRDF_SurfaceReflectivity440, rad_NumTimes, &
       rad_nXtrack, rad_longitude, rad_latitude, rad_time,&
       windspeed2m ,name_option_TemperaturePressure,&
       rad_SnowIceFraction

   use m_vars,only: itdebug, ixdebug, run_mode 

   implicit none
   integer, intent(inout) :: errstat
   type (gler_type) :: glt
   logical :: clip_opt ! wether to limit gler to [0.,1.] range
   real(kind=8) :: thistime
   real(kind=4) :: thislon, thislat, thisalb, wind_speed, thissnowice
   integer :: iwavelen, ix, it, nx, nt, nana
   real(kind=4) :: fspecial

   if (errstat /= 0) return

   fspecial = -999. ! make it negative
   clip_opt = .TRUE.
   write(*,*) '   GLER clip_opt=',clip_opt
 
   nx = rad_nXtrack
   nt = rad_NumTimes

   ! wind_speed set to zero for all options but 'GEOS5'
   ! GMI climatology does not have surface wind speed
   wind_speed = 0.

   !allocate m_vars arrays and initialize
   allocate(BRDF_SurfaceReflectivity466(nx,nt), stat=errstat)
   BRDF_SurfaceReflectivity466 = fspecial

   !------------------------------
   ! calculate GLER for 466nm
   !------------------------------
   iwavelen = 466
   write(*,*) '   Initializing GLER for ',iwavelen
   ! The run-time environment should specify the config file location
   call gler_open(glt, iwavelen, errstat) ! config_file='clim_config.ini')
   if (errstat /= 0) then
     call tell_error (tell_io_error, 'gler_open failed for 466', errstat)
     return
   endif

   !use the starting time of swath in gler_interp_time
   thistime = rad_time(1)
   call gler_interp_time(glt, thistime, errstat)
   if (errstat /= 0) then
     call tell_error (tell_runtime_error, 'gler_interp_time failed', errstat)
     return
   endif

   ! loop through pixels
   if (name_option_TemperaturePressure .eq. 'GEOS5') then
    nana = 0
    do it = 1, nt
      do ix = 1, nx
         thisalb = fspecial
         thislon = rad_longitude(ix,it)
         if ((thislon .lt. -360.) .or. (thislon .gt. 360.)) continue
         thislat = rad_latitude(ix,it)
         if ((thislat .lt. -90.) .or. (thislat .gt. 90.)) continue
         wind_speed = windspeed2m(ix,it)
         if (wind_speed .lt. 0.) continue
         thissnowice = rad_SnowIceFraction(ix,it)
         if ((thissnowice .lt. 0.) .or. (thissnowice .gt. 1.)) continue
         call gler_albedo(glt, thislon, thislat, wind_speed, &
                thissnowice, thisalb, errstat, clip_opt)
         if (isnan(thisalb)) then ! test NAN
                thisalb = fspecial
                nana = nana + 1
         endif
         BRDF_SurfaceReflectivity466(ix,it) = thisalb
      enddo
    enddo
   else
    wind_speed = 0.
    write(*,*)'note:  GMI climatology does not contain wind_speed'
    write(*,*)' GLER is thus calculated with wind_speed=0.'
    nana = 0
    do it = 1, nt
      do ix = 1, nx
         thisalb = fspecial
         thislon = rad_longitude(ix,it)
         if ((thislon .lt. -360.) .or. (thislon .gt. 360.)) continue
         thislat = rad_latitude(ix,it)
         if ((thislat .lt. -90.) .or. (thislat .gt. 90.)) continue
         thissnowice = rad_SnowIceFraction(ix,it)
         if ((thissnowice .lt. 0.) .or. (thissnowice .gt. 1.)) continue
         call gler_albedo(glt, thislon, thislat, wind_speed, &
              thissnowice, thisalb, errstat, clip_opt)
         if (isnan(thisalb)) then ! test NAN
              thisalb = fspecial
              nana = nana + 1
         endif 
         BRDF_SurfaceReflectivity466(ix,it) = thisalb
      enddo
    enddo
   endif

   call gler_close(glt)
   write(*,*) '   GLER n_NAN=',nana

   !------------------------------
   ! calculate GLER for 440nm
   ! 440nm was unavailable during initial development, thus
   ! 440nm is currently not used in calculation, only transfer to output
   !------------------------------
   iwavelen = 440
   allocate(BRDF_SurfaceReflectivity440(nx,nt), stat=errstat)
   BRDF_SurfaceReflectivity440 = fspecial

   write(*,*) '   GLER is not yet implemented for ',iwavelen

   ! debug
   if ((trim(run_mode) .eq. 'development').and. &
       (ixdebug .gt. 0).and.(itdebug .ge. 0)) then
      write(*,*) ' writing debug_gler.txt'
      open(unit=49,file='debug_gler.txt')
      do it = 1, nt
         write(49,*)(BRDF_SurfaceReflectivity466(ix,it),ix=1,nx)
      enddo
      close(49)
   endif

   end subroutine

end module m_read_input_gler
