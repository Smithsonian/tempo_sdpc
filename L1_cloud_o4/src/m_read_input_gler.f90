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
       rad_SnowIceFraction, rad_GroundPixelQualityFlags

   use m_vars,only: itdebug, ixdebug, run_mode, negative999

   use m_vars,only: PerturbAlb466, nord_Alb466pert, Alb466PertCoef

   implicit none
   integer, intent(inout) :: errstat
   type (gler_type) :: glt, glt440
   logical :: clip_opt ! wether to limit gler to [0.,1.] range
   real(kind=8) :: thistime
   real(kind=4) :: thislon, thislat, thisalb, wind_speed, thissnowice
   real(kind=4) :: xalb, yalb
   integer :: iwavelen, ix, it, nx, nt, nana, iord
   real(kind=4) :: fspecial
   integer(kind=2) :: thislandwater
   real (kind=8), parameter :: r8_missval=-1.0d+30

   if (errstat /= 0) return

   fspecial = negative999 ! make it negative
   clip_opt = .TRUE.
   write(*,*) '   GLER clip_opt=',clip_opt

   nx = rad_nXtrack
   nt = rad_NumTimes

   ! wind_speed set to 5m/s for all options but GEOS5
   ! GMI climatology does not have surface wind speed
   wind_speed = 5.

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

   !use the midpoint valid time of swath in gler_interp_time
   thistime = 0.5 * (minval (rad_time, rad_time /= r8_missval) &
                  + maxval(rad_time, rad_time /= r8_missval))
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
         thislandwater = int(ibits(rad_GroundPixelQualityFlags(ix,it), 0, 4), kind=2)
         call gler_albedo(glt, thislon, thislat, thislandwater, wind_speed, &
                thissnowice, thisalb, errstat, clip_opt)
         if (errstat /= 0) then
           errstat = 0
           call tell_set_error (0)
           write (*,*)'gler_albedo failed: lon=',thislon,'lat=',thislat
           thisalb = fspecial
         endif
         if (isnan(thisalb)) then ! test NAN
                thisalb = fspecial
                nana = nana + 1
         endif
         BRDF_SurfaceReflectivity466(ix,it) = thisalb
      enddo
    enddo
   else
    wind_speed = 5.
    write(*,*)'note:  GMI climatology does not contain wind_speed'
    write(*,*)' GLER is thus calculated with wind_speed=5.'
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
         thislandwater = int(ibits(rad_GroundPixelQualityFlags(ix,it), 0, 4), kind=2)
         call gler_albedo(glt, thislon, thislat, thislandwater, wind_speed, &
                thissnowice, thisalb, errstat, clip_opt)
         if (errstat /= 0) then
           errstat = 0
           call tell_set_error (0)
           write (*,*)'gler_albedo failed: lon=',thislon,'lat=',thislat
           thisalb = fspecial
         endif
         if (isnan(thisalb)) then ! test NAN
              thisalb = fspecial
              nana = nana + 1
         endif
         BRDF_SurfaceReflectivity466(ix,it) = thisalb
      enddo
    enddo
   endif

   call gler_close(glt)
   write(*,*) '   GLER 466nm n_NAN=',nana

   ! alb466 perturbation
   if (PerturbAlb466) then
      write(*,*) '   Apply Alb466PertCoef to alb466'
      do it = 1, nt
         do ix = 1, nx
            xalb = BRDF_SurfaceReflectivity466(ix,it)

          if ((xalb .ge. 0.) .and. (xalb .le. 1.)) then 
            yalb = 0.
            do iord = 0, nord_Alb466Pert
               ! fortran index from 1
               yalb = yalb + Alb466PertCoef(iord+1)*(xalb**iord)
            enddo ! iord
            ! clip out of range yalb
            if (yalb .gt. 1.) yalb = 1.0
            if (yalb .lt. 0.) yalb = 0.0
          else
               yalb = fspecial
          endif

           BRDF_SurfaceReflectivity466(ix,it) = yalb
         enddo ! ix
      enddo ! it
   endif  ! PerturbAlb466
     
   !------------------------------
   ! calculate GLER for 440nm
   ! 440nm was unavailable during initial development, activated apr2024
   ! 440nm is currently not used in calculation, only transfer to output
   !------------------------------

   iwavelen = 440
   allocate(BRDF_SurfaceReflectivity440(nx,nt), stat=errstat)
   BRDF_SurfaceReflectivity440 = fspecial

   call gler_open(glt440, iwavelen, errstat) ! config_file='clim_config.ini')
   if (errstat /= 0) then
     call tell_error (tell_io_error, 'gler_open failed for 440', errstat)
     return
   endif

   !use the midpoint valid time of swath in gler_interp_time
   thistime = 0.5 * (minval (rad_time, rad_time /= r8_missval) &
                  + maxval(rad_time, rad_time /= r8_missval))
   call gler_interp_time(glt440, thistime, errstat)
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
         thislandwater = int(ibits(rad_GroundPixelQualityFlags(ix,it), 0, 4), kind=2)
         call gler_albedo(glt440, thislon, thislat, thislandwater, wind_speed, &
                thissnowice, thisalb, errstat, clip_opt)
         if (errstat /= 0) then
           errstat = 0
           call tell_set_error (0)
           write (*,*)'gler_albedo failed: lon=',thislon,'lat=',thislat
           thisalb = fspecial
         endif
         if (isnan(thisalb)) then ! test NAN
                thisalb = fspecial
                nana = nana + 1
         endif
         BRDF_SurfaceReflectivity440(ix,it) = thisalb
      enddo
    enddo
   else
    wind_speed = 5.
    write(*,*)'note:  GMI climatology does not contain wind_speed'
    write(*,*)' GLER is thus calculated with wind_speed=5.'
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
         thislandwater = int(ibits(rad_GroundPixelQualityFlags(ix,it), 0, 4), kind=2)
         call gler_albedo(glt440, thislon, thislat, thislandwater, wind_speed, &
                thissnowice, thisalb, errstat, clip_opt)
         if (errstat /= 0) then
           errstat = 0
           call tell_set_error (0)
           write (*,*)'gler_albedo failed: lon=',thislon,'lat=',thislat
           thisalb = fspecial
         endif
         if (isnan(thisalb)) then ! test NAN
              thisalb = fspecial
              nana = nana + 1
         endif
         BRDF_SurfaceReflectivity440(ix,it) = thisalb
      enddo
    enddo
   endif

   call gler_close(glt440)
   write(*,*) '   GLER440 n_NAN=',nana

   !-------------------------------

888 continue

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
