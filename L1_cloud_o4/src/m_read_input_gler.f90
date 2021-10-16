!Read GLER 
module m_read_input_gler
   use gler_module
   use tell_module

contains
   subroutine read_gler

   use m_vars,only: BRDF_SurfaceReflectivity466,&
       BRDF_SurfaceReflectivity440, rad_NumTimes, &
       rad_nXtrack, rad_longitude, rad_latitude, rad_time,&
       windspeed2m ,name_option_TemperaturePressure,&
       rad_SnowIceFraction

   implicit none
   type (gler_type) :: glt
   logical :: clip_opt
   real(kind=8) :: thistime
   real(kind=4) :: thislon, thislat, thisalb, wind_speed, thissnowice
   integer :: iwavelen, ix, it, nx, nt, errstat
   real(kind=4) :: fspecial

   errstat = 0
   fspecial = -9999.
   clip_opt = .TRUE.

   nx = rad_nXtrack
   nt = rad_NumTimes

   !hqw wind_speed set to zero for all other options but 'GEOS5'
   wind_speed = 0.

   !allocate m_vars arrays and initialize
   allocate(BRDF_SurfaceReflectivity466(nx,nt), stat=errstat)
   BRDF_SurfaceReflectivity466 = -9.9

   !------------------------------
   ! calculate GLER for 466nm
   !------------------------------
   iwavelen = 466
!   write(*,*) 'Initializing GLER for 466'
   call gler_open(glt, iwavelen, errstat, config_file='clim_config.ini')
   if (errstat /=0) then
      write(*,*) 'gler_open failed for 466'
      stop 1
   endif

   !use the starting time of swath in gler_interp_time
   thistime = rad_time(1)
   call gler_interp_time(glt, thistime, errstat)
   if (errstat /= 0) then
       write(*,*)'gler_interp_time failed'
       stop 1
   endif

   ! loop through pixels
   if (name_option_TemperaturePressure .eq. 'GEOS5') then
    do it = 1, nt
      do ix = 1, nx
         thislon = rad_longitude(ix,it)
         if ((thislon .lt. -360.) .or. (thislon .gt. 360.)) continue
         thislat = rad_latitude(ix,it)
         if ((thislat .lt. -90.) .or. (thislat .gt. 90.)) continue
         wind_speed = windspeed2m(ix,it)
         if (wind_speed .lt. 0.) continue
         thissnowice = rad_SnowIceFraction(ix,it)
         if ((thissnowice .lt. 0.) .or. (thissnowice .gt. 1.)) continue
!hqw debug
!         if (ix .eq. 40) write(*,*) ix, it, wind_speed
         call gler_albedo(glt, thislon, thislat, wind_speed, &
                thissnowice, thisalb, errstat, clip_opt)
         BRDF_SurfaceReflectivity466(ix,it) = thisalb
      enddo
    enddo     
   else
    wind_speed = 0.
    write(*,*)'note: GLER calculate with wind_speed=0.'
    write(*,*)'    because GMI climatology does not contain wind_speed info'
    do it = 1, nt
      do ix = 1, nx
         thislon = rad_longitude(ix,it)
         if ((thislon .lt. -360.) .or. (thislon .gt. 360.)) continue
         thislat = rad_latitude(ix,it)
         if ((thislat .lt. -90.) .or. (thislat .gt. 90.)) continue
         thissnowice = rad_SnowIceFraction(ix,it)
         if ((thissnowice .lt. 0.) .or. (thissnowice .gt. 1.)) continue
         call gler_albedo(glt, thislon, thislat, wind_speed, &
                thissnowice, thisalb, errstat, clip_opt)
         BRDF_SurfaceReflectivity466(ix,it) = thisalb
      enddo
    enddo
   endif 

   !------------------------------
   ! calculate GLER for 440nm
   ! 440nm is currently not used in cloud calculation, only transfer to output
   !------------------------------
   iwavelen = 440
   allocate(BRDF_SurfaceReflectivity440(nx,nt), stat=errstat)
   BRDF_SurfaceReflectivity440 = fspecial

   call gler_close(glt)

   !hqw debug
   !write(*,*) 'write debug_gler.txt'
   !open(unit=19,file='debug_gler.txt')
   !do it = 1, nt
   !   write(19,*)(BRDF_SurfaceReflectivity466(ix,it),ix=1,nx) 
   !enddo
   !close(19)

   end subroutine

end module m_read_input_gler
