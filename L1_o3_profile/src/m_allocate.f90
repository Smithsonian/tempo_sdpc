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
  public deallocate_irrad, deallocate_rad, deallocate_ring, deallocate_refl
  public deallocate_cali
  PRIVATE

CONTAINS

  subroutine deallocate_irrad (irrad)
    implicit none
    type (irrad_group), intent(inout) :: irrad

    IF (ALLOCATED (irrad%nwav)) DEALLOCATE (irrad%nwav)
    IF (ALLOCATED (irrad%errstat)) DEALLOCATE (irrad%errstat)
    IF (ALLOCATED (irrad%npix)) DEALLOCATE (irrad%npix)
    IF (ALLOCATED (irrad%winpix)) DEALLOCATE (irrad%winpix)
    IF (ALLOCATED (irrad%wind)) DEALLOCATE (irrad%wind)
    IF (ALLOCATED (irrad%wavl)) DEALLOCATE (irrad%wavl)
    IF (ALLOCATED (irrad%spec)) DEALLOCATE (irrad%spec)
    IF (ALLOCATED (irrad%prec)) DEALLOCATE (irrad%prec)
    IF (ALLOCATED (irrad%qflg)) DEALLOCATE (irrad%qflg)
    IF (ALLOCATED (irrad%norm)) DEALLOCATE (irrad%norm)

   return
 end subroutine deallocate_irrad

 subroutine deallocate_rad (rad)
   implicit none
   type (rad_group), intent(inout) :: rad

   IF (ALLOCATED (rad%nwav)) DEALLOCATE (rad%nwav)
   IF (ALLOCATED (rad%npix)) DEALLOCATE (rad%npix)
   IF (ALLOCATED (rad%errstat)) DEALLOCATE (rad%errstat)
   IF (ALLOCATED (rad%pix_errstat)) DEALLOCATE (rad%pix_errstat)
   IF (ALLOCATED (rad%wind)) DEALLOCATE (rad%wind)
   IF (ALLOCATED (rad%wavl)) DEALLOCATE (rad%wavl)
   IF (ALLOCATED (rad%spec)) DEALLOCATE (rad%spec)
   IF (ALLOCATED (rad%prec)) DEALLOCATE (rad%prec)
   IF (ALLOCATED (rad%qflg)) DEALLOCATE (rad%qflg)
   IF (ALLOCATED (rad%norm)) DEALLOCATE (rad%norm)

   return
 end subroutine deallocate_rad

 subroutine deallocate_ring (ring)
   implicit none
   type (ring_group), intent(inout) :: ring

   IF (ALLOCATED (ring%spec)) DEALLOCATE (ring%spec)
   IF (ALLOCATED (ring%wavl)) DEALLOCATE (ring%wavl)
   IF (ALLOCATED (ring%nsol)) DEALLOCATE (ring%nsol)
   IF (ALLOCATED (ring%ndiv)) DEALLOCATE (ring%ndiv)
   IF (ALLOCATED (ring%winpix)) DEALLOCATE (ring%winpix)

   return
 end subroutine deallocate_ring

 subroutine deallocate_refl (refl)
   implicit none
   type (refl_group), intent(inout) :: refl

   IF (ALLOCATED (refl%winpix)) DEALLOCATE (refl%winpix)
   IF (ALLOCATED (refl%solspec)) DEALLOCATE (refl%solspec)
   IF (ALLOCATED (refl%solwavl)) DEALLOCATE (refl%solwavl)
   IF (ALLOCATED (refl%radspec)) DEALLOCATE (refl%radspec)
   IF (ALLOCATED (refl%radwavl)) DEALLOCATE (refl%radwavl)

   return
 end subroutine deallocate_refl

 subroutine deallocate_cali (cali)
   implicit none
   type (cali_group), intent(inout) :: cali

   if (allocated (cali%wincal_wav)) deallocate(cali%wincal_wav)
   if (allocated (cali%solwinfit)) deallocate(cali%solwinfit)
   if (allocated (cali%radwinfit)) deallocate(cali%radwinfit)
   if (allocated (cali%nslit_sol)) deallocate(cali%nslit_sol)
   if (allocated (cali%nslit_rad)) deallocate(cali%nslit_rad)
   if (allocated (cali%slitfit_sol)) deallocate(cali%slitfit_sol)
   if (allocated (cali%slitfit_rad)) deallocate(cali%slitfit_rad)
   if (allocated (cali%nwavcal_sol)) deallocate(cali%nwavcal_sol)
   if (allocated (cali%nwavcal_rad)) deallocate(cali%nwavcal_rad)
   if (allocated (cali%sswav_sol)) deallocate(cali%sswav_sol)
   if (allocated (cali%sswav_rad)) deallocate(cali%sswav_rad)

   return
 end subroutine deallocate_cali

  SUBROUTINE allocate_spec (nwin, nx, ny, irrad,rad,  ring, refl, cali, status)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx, ny, nwin
  TYPE(irrad_group), INTENT(INOUT) :: irrad
  TYPE(rad_group), INTENT(INOUT) :: rad
  TYPE(ring_group), INTENT(INOUT) :: ring
  TYPE(refl_group), INTENT(INOUT) :: refl
  TYPE(cali_group), INTENT(INOUT) :: cali
  INTEGER, INTENT(OUT) :: status
  CHARACTER(64) :: varname

  call deallocate_irrad (irrad)

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
   call deallocate_rad (rad)

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
   call deallocate_ring (ring)

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
   call deallocate_refl (refl)

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
   IF (ALLOCATED (cali%wincal_wav)) DEALLOCATE (cali%wincal_wav)
   ALLOCATE (cali%wincal_wav(nwin, nx), stat=status); IF (status .ne. 0 ) GOTO 111    

   IF (wavcal) THEN      
     IF (ALLOCATED (cali%radwinfit)) DEALLOCATE (cali%radwinfit)
     ALLOCATE (cali%radwinfit(nwin, max_calfit_idx, nx), stat=status)
   ENDIF
   IF (status .ne. 0 ) GOTO 111    

   varname='cali%solwinfit'
   IF (ALLOCATED (cali%solwinfit)) DEALLOCATE (cali%solwinfit)
   ALLOCATE (cali%solwinfit(nwin, max_calfit_idx, nx), stat=status)
   IF (status .ne. 0 ) GOTO 111    


   IF (.NOT. slit_rad .AND. .NOT. yn_varyslit) return

     IF (slit_rad) THEN 
        IF (ALLOCATED (cali%radwinfit)) DEALLOCATE (cali%radwinfit)
        ALLOCATE (cali%radwinfit(nwin, max_calfit_idx, nx), stat=status)
        IF (status .ne. 0 ) GOTO 111    
     ENDIF 
      
     IF (yn_varyslit) THEN
       IF (ALLOCATED (cali%nslit_sol)) DEALLOCATE (cali%nslit_sol)
       ALLOCATE (cali%nslit_sol(nx), stat=status); IF (status .ne. 0 ) GOTO 111    

       IF (ALLOCATED (cali%slitwav_sol)) DEALLOCATE (cali%slitwav_sol)
       ALLOCATE (cali%slitwav_sol(max_fit_pts, nx), stat=status); IF (status .ne. 0 ) GOTO 111    

       IF (ALLOCATED (cali%slitfit_sol)) DEALLOCATE (cali%slitfit_sol)
       ALLOCATE (cali%slitfit_sol(max_fit_pts, max_calfit_idx, 2, nx), stat=status) 
       IF (status .ne. 0 ) GOTO 111    
       
       IF (wavcal .OR. wavcal_sol) THEN  
         IF (ALLOCATED (cali%nwavcal_sol)) DEALLOCATE (cali%nwavcal_sol)
         ALLOCATE (cali%nwavcal_sol(nx), stat=status); IF (status .ne. 0 ) GOTO 111    

         IF (ALLOCATED (cali%sswav_sol)) DEALLOCATE (cali%sswav_sol)
         ALLOCATE (cali%sswav_sol(max_fit_pts, nx), stat=status)
       ENDIF

       IF (status .ne. 0 ) GOTO 111    

         IF ( slit_rad) THEN 
            IF (ALLOCATED (cali%nslit_rad)) DEALLOCATE (cali%nslit_rad)
            ALLOCATE (cali%nslit_rad(nx), stat=status)
            IF (status .ne. 0 ) GOTO 111    

            IF (ALLOCATED (cali%slitwav_rad)) DEALLOCATE (cali%slitwav_rad)
            ALLOCATE (cali%slitwav_sol(max_fit_pts, nx), stat=status)
            IF (status .ne. 0 ) GOTO 111    
 
            IF (ALLOCATED (cali%slitfit_rad)) DEALLOCATE (cali%slitfit_rad)
            ALLOCATE (cali%slitfit_rad(max_fit_pts, max_calfit_idx, 2, nx), stat=status)
            IF (status .ne. 0 ) GOTO 111    
         ENDIF

         IF (wavcal) THEN 
           IF (ALLOCATED (cali%nwavcal_rad)) DEALLOCATE (cali%nwavcal_rad)
           ALLOCATE (cali%nwavcal_rad(nx), stat=status); IF (status .ne. 0 ) GOTO 111    

           IF (ALLOCATED (cali%sswav_rad)) DEALLOCATE (cali%sswav_rad)
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

   IF (ALLOCATED (geo%time)) DEALLOCATE (geo%time)
   IF (ALLOCATED (geo%lon)) DEALLOCATE (geo%lon)
   IF (ALLOCATED (geo%lat)) DEALLOCATE (geo%lat)
   IF (ALLOCATED (geo%sza)) DEALLOCATE (geo%sza)
   IF (ALLOCATED (geo%vza)) DEALLOCATE (geo%vza)
   IF (ALLOCATED (geo%aza)) DEALLOCATE (geo%aza)
   IF (ALLOCATED (geo%sca)) DEALLOCATE (geo%sca)
   IF (ALLOCATED (geo%clon)) DEALLOCATE (geo%clon)
   IF (ALLOCATED (geo%clat)) DEALLOCATE (geo%clat)
   IF (ALLOCATED (geo%elon)) DEALLOCATE (geo%elon)
   IF (ALLOCATED (geo%elat)) DEALLOCATE (geo%elat)
   IF (ALLOCATED (geo%height)) DEALLOCATE (geo%height)
   IF (ALLOCATED (geo%cfr)) DEALLOCATE (geo%cfr)
   IF (ALLOCATED (geo%ctp)) DEALLOCATE (geo%ctp)
   IF (ALLOCATED (geo%cloud_qflg)) DEALLOCATE (geo%cloud_qflg)
   IF (ALLOCATED (geo%xflg)) DEALLOCATE (geo%xflg)
   IF (ALLOCATED (geo%gflg)) DEALLOCATE (geo%gflg)
   IF (ALLOCATED (geo%land_water_flg)) DEALLOCATE (geo%land_water_flg)
   IF (ALLOCATED (geo%snow_ice_flg)) DEALLOCATE (geo%snow_ice_flg)
   IF (ALLOCATED (geo%glint_flg)) DEALLOCATE (geo%glint_flg)

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

  SUBROUTINE dealloc (errstat)

    implicit none

    !output variables
    integer (kind=4), intent(inout) :: errstat

  END SUBROUTINE dealloc


end module m_allocate
