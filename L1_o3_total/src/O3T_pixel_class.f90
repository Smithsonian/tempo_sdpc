!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_pixel_class
! 
! This module defines the pixel geo-angle data type (O3T_pixgeo_type) 
! and pixel surface cover type (O3T_pixcover_type). And a set of functions
! to set the field values of these data types.
!
!!Revision History:
! Initial version 03/26/2002  Kai Yang/UMBC
!
!!Team-unique Header:
! This software was developed by the OMI Science Team Support
! Group for the National Aeronautics and Space Administration, Goddard
! Space Flight Center, under NASA Task 916-003-1
!
!!References and Credits
! Written by 
! Kai Yang 
! University of Maryland Baltimore County
! email: Kai.Yang-1@nasa.gov
! 
!!Design Notes
!
!!END
!!****************************************************************************
MODULE O3T_pixel_class
    IMPLICIT NONE
    REAL (KIND = 8),PARAMETER::DEGtoRAD=0.017453292519943295769236907684886
    REAL (KIND = 8),PARAMETER::RADtoDEG=57.29577951308232087679815481410517

    PUBLIC :: O3T_pixgeo_set
!    PUBLIC :: O3T_pixrefl_set
    
    TYPE, PUBLIC :: O3T_pixgeo_type
      ! values assocaitated with a set of solar and viwing zenith angles.
      INTEGER (KIND = 4) :: id
      REAL (KIND = 8) :: sza, vza, phi
      REAL (KIND = 8) :: p1, pr
      REAL (KIND = 8) ::    cos_sza, sin_sza
      REAL (KIND = 8) :: lg_cos_sza
      REAL (KIND = 8) ::    cos_vza, sin_vza
      REAL (KIND = 8) :: lg_cos_vza
      REAL (KIND = 8) :: pt, ptp, pc, pcp
      REAL (KIND = 8) :: log_pt, log_ptp, log_pc, log_pcp
      REAL (KIND = 8) :: cphi, c2phi
      REAL (KIND = 4) :: fteran
      REAL (KIND = 4) :: lat, lon
    END TYPE O3T_pixgeo_type

    TYPE, PUBLIC :: O3T_pixcover_type
      REAL (KIND=4) :: grref,  clref
      REAL (KIND=4) :: clfrac, pclfrac
      REAL (KIND=4) :: rcf1  !! Radiative cloud fraction = clfrac*(0.8/ref360)
      REAL (KIND=4) :: rcf2  !! Radiative cloud fraction = clfrac*(Ic331/Im331)
      REAL (KIND=4) :: ref,    pref
      REAL (KIND=4) :: ref360
      
      INTEGER (KIND=4) :: iwl_glint, iwl_refl_l, iwl_refl_h
      LOGICAL (KIND=4) :: glint_flag
      REAL (KIND=4) :: wl_glint, wl_refl_l, wl_refl_h
      REAL (KIND=4) :: nres_glnt_ub
      REAL (KIND=4) :: glnfrc

      INTEGER (KIND=4) :: isnow
      INTEGER (KIND=1) :: landsea_mask
      INTEGER (KIND=1) :: surface_category
      
    END TYPE O3T_pixcover_type

    CONTAINS

      FUNCTION O3T_pixgeo_set( pix_sza, pix_vza, pix_phi, &
                              pix_pt, pix_pc, lat, lon, &
                              pixel_id ) RESULT (this)
        REAL( KIND = 4), INTENT(IN) :: pix_sza, pix_vza, pix_phi, pix_pt, pix_pc
        REAL( KIND = 4), INTENT(IN) :: lat, lon
        INTEGER( KIND = 4), INTENT(IN), OPTIONAL :: pixel_id
        TYPE (O3T_pixgeo_type) :: this

        IF( PRESENT( pixel_id ) ) THEN
           this%id = pixel_id
        ELSE
           this%id = -32767
        ENDIF

        this%sza = pix_sza*DEGtoRAD
        this%vza = pix_vza*DEGtoRAD
        this%phi = pix_phi*DEGtoRAD
        this%pt  = pix_pt
        this%ptp = pix_pt*0.99
        this%pc  = pix_pc
        this%pcp = pix_pc*1.01
        IF( this%pcp > 1.0 ) this%pcp = pix_pc*0.99
        this%lat = lat
        this%lon = lon
        this%fteran = REAL((1.d0-this%pt)/0.5d0) ! calculate the terrain factor
         
        this%cos_sza = COS( this%sza )
        this%sin_sza = SIN( this%sza )
        this%cos_vza = COS( this%vza )
        this%sin_vza = SIN( this%vza )
        this%p1 = -0.375D0*this%cos_sza*this%sin_sza*this%sin_vza
        this%pr = 2.0D0*this%p1/(3.D0*this%cos_vza*this%cos_sza**2)
  
        this%lg_cos_sza = -LOG( this%cos_sza )
        this%lg_cos_vza = -LOG( this%cos_vza )
        this%log_pt = LOG10(this%pt)
        this%log_ptp= LOG10(this%ptp)
        this%log_pc = LOG10(this%pc)
        this%log_pcp= LOG10(this%pcp)
        this%cphi   = COS( this%phi )
        this%c2phi  = COS( 2.D0*this%phi )

      END FUNCTION O3T_pixgeo_set

!      FUNCTION O3T_pixrefl_set( grref, clref, clfrac ) RESULT(this)
!        REAL(KIND=4), INTENT(IN) :: grref, clref, clfrac
!        TYPE(O3T_pixcover_type) :: this
!        this%grref = grref
!        this%clref = clref
!        this%clfrac = clfrac
!        this%ref = (1.0-clfrac)*grref + clfrac*clref
!      END FUNCTION O3T_pixrefl_set

END MODULE O3T_pixel_class
