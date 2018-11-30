!> Routines to allocate and deallocate arrays
MODULE m_allocate

  USE OMSAO_indices_module, only:  max_calfit_idx, n_max_fitpars, & 
      omi_idx,instrument_idx
  USE OMSAO_parameters_module, only: max_fit_pts, max_ring_pts, mrefl
  USE OMSAO_variables_module, ONLY: &
      slit_rad, wavcal, wavcal_sol,yn_varyslit, &
      rad_group, irrad_group, refl_group, ring_group, cali_group, geo_group

  IMPLICIT NONE

  PUBLIC allocate_spec, allocate_geo, dealloc
  PRIVATE

CONTAINS

  SUBROUTINE allocate_spec (nwin, nx, ny, irrad,rad,  ring, refl, cali, status)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx, ny, nwin
  TYPE(irrad_group), INTENT(INOUT) :: irrad
  TYPE(rad_group), INTENT(INOUT) :: rad
  TYPE(ring_group), INTENT(INOUT) :: ring
  TYPE(refl_group), INTENT(INOUT) :: refl
  TYPE(cali_group), INTENT(INOUT) :: cali
  INTEGER, INTENT(OUT) :: status
  CHARACTER(10) :: varname

   IF (ASSOCIATED (irrad%nwav)) DEALLOCATE (irrad%nwav)
   IF (ASSOCIATED (irrad%errstat)) DEALLOCATE (irrad%errstat)
   IF (ASSOCIATED (irrad%npix)) DEALLOCATE (irrad%npix)
   IF (ASSOCIATED (irrad%winpix)) DEALLOCATE (irrad%winpix)
   IF (ASSOCIATED (irrad%wind)) DEALLOCATE (irrad%wind)
   IF (ASSOCIATED (irrad%wavl)) DEALLOCATE (irrad%wavl)
   IF (ASSOCIATED (irrad%spec)) DEALLOCATE (irrad%spec)
   IF (ASSOCIATED (irrad%prec)) DEALLOCATE (irrad%prec)
   IF (ASSOCIATED (irrad%qflg)) DEALLOCATE (irrad%qflg)
   IF (ASSOCIATED (irrad%norm)) DEALLOCATE (irrad%norm)

   varname='irrad%nwav'
   ALLOCATE (irrad%nwav(nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%errstat'
   ALLOCATE (irrad%errstat(nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%npix'
   ALLOCATE (irrad%npix(nwin, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%winpix'
   ALLOCATE (irrad%winpix(nwin, nx, 2), stat=status); IF (status .ne. 0 )GOTO 111

   varname='irrad%wind'
   ALLOCATE (irrad%wind(max_fit_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%wavl'
   ALLOCATE (irrad%wavl(max_fit_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%spec'
   ALLOCATE (irrad%spec(max_fit_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%prec'
   ALLOCATE (irrad%prec(max_fit_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%qflg'
   ALLOCATE (irrad%qflg(max_fit_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='irrad%norm'
   ALLOCATE (irrad%norm(nx), stat=status); IF (status .ne. 0 ) GOTO 111

   !--------------------------------------------------------------------------
   ! radiance
   !-------------------------------------------------------------------------
   IF (ASSOCIATED (rad%nwav)) DEALLOCATE (rad%nwav)
   IF (ASSOCIATED (rad%npix)) DEALLOCATE (rad%npix)
   IF (ASSOCIATED (rad%errstat)) DEALLOCATE (rad%errstat)
   IF (ASSOCIATED (rad%pix_errstat)) DEALLOCATE (rad%pix_errstat)
   IF (ASSOCIATED (rad%wind)) DEALLOCATE (rad%wind)
   IF (ASSOCIATED (rad%wavl)) DEALLOCATE (rad%wavl)
   IF (ASSOCIATED (rad%spec)) DEALLOCATE (rad%spec)
   IF (ASSOCIATED (rad%prec)) DEALLOCATE (rad%prec)
   IF (ASSOCIATED (rad%qflg)) DEALLOCATE (rad%qflg)
   IF (ASSOCIATED (rad%norm)) DEALLOCATE (rad%norm)

   varname='rad%nwav'
   ALLOCATE (rad%nwav(nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
  
   varname='rad%npix'
   ALLOCATE (rad%npix(nwin, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
  
   varname='rad%errstat'
   ALLOCATE (rad%errstat(0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   
   varname='rad%pix_errstat'
   ALLOCATE (rad%pix_errstat(nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   
   varname='rad%wind'

   ALLOCATE (rad%wind(max_fit_pts, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   
   varname='rad%wavl'
   ALLOCATE (rad%wavl(max_fit_pts, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   
   varname='rad%spec'
   ALLOCATE (rad%spec(max_fit_pts, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   
   varname='rad%prec'
   ALLOCATE (rad%prec(max_fit_pts, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   
   varname='rad%qflg'
   ALLOCATE (rad%qflg(max_fit_pts, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   
   varname='rad%norm'
   ALLOCATE (rad%norm(nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111

   !--------------------------------------------------------------------------
   ! RING 
   !-------------------------------------------------------------------------
   IF (ASSOCIATED (ring%spec)) DEALLOCATE (ring%spec)
   IF (ASSOCIATED (ring%wavl)) DEALLOCATE (ring%wavl)
   IF (ASSOCIATED (ring%nsol)) DEALLOCATE (ring%nsol)
   IF (ASSOCIATED (ring%ndiv)) DEALLOCATE (ring%ndiv)
   IF (ASSOCIATED (ring%winpix)) DEALLOCATE (ring%winpix)

   varname='ring%spec'
   ALLOCATE (ring%spec(max_ring_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='ring%wavl'

   ALLOCATE (ring%wavl(max_ring_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='ring%nsol'
   ALLOCATE (ring%nsol(nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='ring%ndiv'
   ALLOCATE (ring%ndiv(nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='ring%winpix'
   ALLOCATE (ring%winpix(nx,2), stat=status); IF (status .ne. 0 ) GOTO 111


   !--------------------------------------------------------------------------
   ! reflectance
   !--------------------------------------------------------------------------
   IF (ASSOCIATED (refl%winpix)) DEALLOCATE (refl%winpix)
   IF (ASSOCIATED (refl%solspec)) DEALLOCATE (refl%solspec)
   IF (ASSOCIATED (refl%solwavl)) DEALLOCATE (refl%solwavl)
   IF (ASSOCIATED (refl%radspec)) DEALLOCATE (refl%radspec)
   IF (ASSOCIATED (refl%radwavl)) DEALLOCATE (refl%radwavl)

   varname='refl%winpix'
   ALLOCATE (refl%winpix(nx,2), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='refl%solspec'
   ALLOCATE (refl%solspec(mrefl, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='refl%solwavl'
   ALLOCATE (refl%solwavl(mrefl, nx), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='refl%radspec'
   ALLOCATE (refl%radspec(mrefl, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111

   varname='refl%radwavl'
   ALLOCATE (refl%radwavl(mrefl, nx, 0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111

   !--------------------------------------------------------------------------
   ! calibration
   !--------------------------------------------------------------------------
   varname = 'cali%wincal_wav'
   IF (ASSOCIATED (cali%wincal_wav)) DEALLOCATE (cali%wincal_wav)    
   ALLOCATE (cali%wincal_wav(nwin, nx), stat=status); IF (status .ne. 0 ) GOTO 111    

   IF (wavcal) THEN      
     IF (ASSOCIATED (cali%radwinfit)) DEALLOCATE (cali%radwinfit)
     ALLOCATE (cali%radwinfit(nwin, max_calfit_idx, nx), stat=status)
   ENDIF
   IF (status .ne. 0 ) GOTO 111    

   varname='cali%solwinfit'
   IF (ASSOCIATED (cali%solwinfit)) DEALLOCATE (cali%solwinfit)
   ALLOCATE (cali%solwinfit(nwin, max_calfit_idx, nx), stat=status)
   IF (status .ne. 0 ) GOTO 111    


   IF (.NOT. slit_rad .AND. .NOT. yn_varyslit) return

     IF (slit_rad) THEN 
        IF (ASSOCIATED (cali%radwinfit)) DEALLOCATE (cali%radwinfit)
        ALLOCATE (cali%radwinfit(nwin, max_calfit_idx, nx), stat=status)
        IF (status .ne. 0 ) GOTO 111    
     ENDIF 
      
     IF (yn_varyslit) THEN
       IF (ASSOCIATED (cali%nslit_sol)) DEALLOCATE (cali%nslit_sol)
       ALLOCATE (cali%nslit_sol(nx), stat=status); IF (status .ne. 0 ) GOTO 111    

       IF (ASSOCIATED (cali%slitwav_sol)) DEALLOCATE (cali%slitwav_sol)
       ALLOCATE (cali%slitwav_sol(max_fit_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111    

       IF (ASSOCIATED (cali%slitfit_sol)) DEALLOCATE (cali%slitfit_sol)
       ALLOCATE (cali%slitfit_sol(max_fit_pts, max_calfit_idx, 2, nx), stat=status) 
       IF (status .ne. 0 ) GOTO 111    
       
       IF (wavcal .OR. wavcal_sol) THEN  
         IF (ASSOCIATED (cali%nwavcal_sol)) DEALLOCATE (cali%nwavcal_sol)
         ALLOCATE (cali%nwavcal_sol(nx), stat=status); IF (status .ne. 0 ) GOTO 111    

         IF (ASSOCIATED (cali%sswav_sol)) DEALLOCATE (cali%sswav_sol)
         ALLOCATE (cali%sswav_sol(max_fit_pts, nx), stat=status)
       ENDIF

       IF (status .ne. 0 ) GOTO 111    

         IF ( slit_rad) THEN 
            IF (ASSOCIATED (cali%nslit_rad)) DEALLOCATE (cali%nslit_rad)
            ALLOCATE (cali%nslit_rad(nx), stat=status)
            IF (status .ne. 0 ) GOTO 111    

            IF (ASSOCIATED (cali%slitwav_rad)) DEALLOCATE (cali%slitwav_rad)
            ALLOCATE (cali%slitwav_sol(max_fit_pts, nx), stat=status)
            IF (status .ne. 0 ) GOTO 111    
 
            IF (ASSOCIATED (cali%slitfit_rad)) DEALLOCATE (cali%slitfit_rad)
            ALLOCATE (cali%slitfit_rad(max_fit_pts, max_calfit_idx, 2, nx), stat=status)
            IF (status .ne. 0 ) GOTO 111    
         ENDIF

         IF (wavcal) THEN 
           IF (ASSOCIATED (cali%nwavcal_rad)) DEALLOCATE (cali%nwavcal_rad)
           ALLOCATE (cali%nwavcal_rad(nx), stat=status); IF (status .ne. 0 ) GOTO 111    

           IF (ASSOCIATED (cali%sswav_rad)) DEALLOCATE (cali%sswav_rad)
           ALLOCATE (cali%sswav_rad(max_fit_pts, nx), stat=status)
           IF (status .ne. 0 ) GOTO 111    
         ENDIF
     ENDIF

   RETURN

   111 continue
   WRITE(*,*) "allocation errors with "//TRIM(ADJUSTL(varname))
   status  = -1
   RETURN

  END SUBROUTINE allocate_spec

  SUBROUTINE allocate_geo (nx, ny, geo, status)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nx, ny
   TYPE(geo_group), INTENT(INOUT) :: geo
   INTEGER, INTENT(OUT) :: status   

   IF (ASSOCIATED (geo%time)) DEALLOCATE (geo%time)
   IF (ASSOCIATED (geo%lon)) DEALLOCATE (geo%lon)
   IF (ASSOCIATED (geo%lat)) DEALLOCATE (geo%lat)
   IF (ASSOCIATED (geo%sza)) DEALLOCATE (geo%sza)
   IF (ASSOCIATED (geo%vza)) DEALLOCATE (geo%vza)
   IF (ASSOCIATED (geo%aza)) DEALLOCATE (geo%aza)
   IF (ASSOCIATED (geo%sca)) DEALLOCATE (geo%sca)
   IF (ASSOCIATED (geo%clon)) DEALLOCATE (geo%clon)
   IF (ASSOCIATED (geo%clat)) DEALLOCATE (geo%clat)
   IF (ASSOCIATED (geo%elon)) DEALLOCATE (geo%elon)
   IF (ASSOCIATED (geo%elat)) DEALLOCATE (geo%elat)
   IF (ASSOCIATED (geo%height)) DEALLOCATE (geo%height)
   IF (ASSOCIATED (geo%cfr)) DEALLOCATE (geo%cfr)
   IF (ASSOCIATED (geo%ctp)) DEALLOCATE (geo%ctp)
   IF (ASSOCIATED (geo%cloud_qflg)) DEALLOCATE (geo%cloud_qflg)
   IF (ASSOCIATED (geo%xflg)) DEALLOCATE (geo%xflg)
   IF (ASSOCIATED (geo%gflg)) DEALLOCATE (geo%gflg)
   IF (ASSOCIATED (geo%land_water_flg)) DEALLOCATE (geo%land_water_flg)
   IF (ASSOCIATED (geo%snow_ice_flg)) DEALLOCATE (geo%snow_ice_flg)
   IF (ASSOCIATED (geo%glint_flg)) DEALLOCATE (geo%glint_flg)

   ALLOCATE (geo%time(0:ny-1), stat=status); IF (status .ne. 0 ) GOTO 111
   ALLOCATE (geo%lon(nx,0:ny-1),geo%lat(nx,0:ny-1),geo%sza(nx,0:ny-1), &
             geo%vza(nx,0:ny-1),geo%aza(nx,0:ny-1),geo%sca(nx,0:ny-1), &
             geo%cfr(nx,0:ny-1),geo%ctp(nx,0:ny-1),geo%ai(nx,0:ny-1), &            
             geo%height(nx,  0:ny-1), geo%gflg(nx, 0:ny-1), &
             geo%land_water_flg(nx,0:ny-1), geo%snow_ice_flg(nx,0:ny-1), geo%glint_flg(nx,0:ny-1), &
             geo%cloud_qflg(nx,0:ny-1), stat=status)              
   IF (instrument_idx == omi_idx) THEN 
      ALLOCATE(geo%xflg(nx, 0:ny-1))
   ENDIF
   IF (status .ne. 0 ) GOTO 111

   ALLOCATE (geo%clon(4, nx, 0:ny-1), geo%clat(4, nx, 0:ny-1), stat=status)
   ALLOCATE (geo%elon(0:nx, 0:ny-1), geo%elat(0:nx, 0:ny-1), stat=status)  
  
   RETURN
   111 continue
   status  = -1
   RETURN
  END SUBROUTINE allocate_geo

  SUBROUTINe dealloc (errstat)

    implicit none

    !output variables
    integer (kind=4), intent(inout) :: errstat

  END SUBROUTINE dealloc


end module m_allocate
