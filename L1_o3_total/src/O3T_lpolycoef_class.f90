!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_lpolycoef_class
! 
! This module contains a set of functions to calculate the coefficients
! used in Lagrangian interpolation.
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
MODULE O3T_lpolycoef_class
    USE OMI_SMF_class 
    USE O3T_pixel_class
    USE O3T_nval_class, ONLY : nval_read, nvzaT, nszaT, npresT, &
                                sza, vza, pressure 
    IMPLICIT NONE
    
    INTEGER (KIND = 4), PARAMETER :: NINTERP = 4
    INTEGER (KIND = 4), PARAMETER, PRIVATE :: zero = 0

    PUBLIC :: O3T_lpoly_cden
    PUBLIC :: O3T_lpolycden_dispose
    PUBLIC :: O3T_lpoly_coef
    PRIVATE:: get_index
    PRIVATE:: Lnk_denum

    TYPE, PUBLIC :: O3T_lpoly_cden_type
      ! values associated with with the grid points of LUT
      INTEGER (KIND = 4) :: nsza, nvza
      INTEGER (KIND = 4) :: npres
      LOGICAL :: setLUTcoefDen
  
      ! values needed for interpoaltion
      REAL (KIND = 8), DIMENSION(:), POINTER :: lg_csza, lg_cvza, log_pres
      REAL (KIND = 8), DIMENSION(:,:), POINTER :: densza, denvza
      REAL (KIND = 8), DIMENSION(:), POINTER :: denpres 
    END TYPE O3T_lpoly_cden_type

    TYPE, PUBLIC :: O3T_lpoly_coef_type
      ! values assocaitated with a set of solar and viwing zenith angles.
      INTEGER (KIND = 4) :: pixel_id
      INTEGER (KIND = 4) :: isza, ivza
      INTEGER (KIND = 4) :: isza_end, ivza_end
      REAL (KIND = 8), DIMENSION(NINTERP) :: cthet0, cthet
      REAL (KIND = 8), DIMENSION(NINTERP) :: pgr, pgrp, pcl, pclp
      REAL (KIND = 8), DIMENSION(NINTERP*NINTERP) :: Lnk
    END TYPE O3T_lpoly_coef_type

    CONTAINS

      FUNCTION O3T_lpoly_cden( this ) RESULT (status)
        TYPE (O3T_lpoly_cden_type), INTENT( OUT ) :: this
        INTEGER (KIND = 4) :: status
        INTEGER  :: np, j !,k, i,
        INTEGER  :: ierr
        !REAL (KIND = 8) :: xdenom
 
        IF( .NOT. nval_read ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, &
                                 "LUT not initialized", &
                                 "O3T_lpoly_cden", zero )
           RETURN
        ENDIF

        this%nsza  = nszaT
        this%nvza  = nvzaT
        this%npres = npresT

        ALLOCATE( this%lg_csza(this%nsza), this%lg_cvza(this%nvza), &
                  this%log_pres(this%npres), STAT = ierr )

        IF( ierr .NE. zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                 "this%.. allocation failure", &
                                 "O3T_lpoly_cden", zero )
           RETURN
        ENDIF

        this%lg_csza  = -LOG( COS(DBLE(sza*DEGtoRAD)) )
        this%lg_cvza  = -LOG( COS(DBLE(vza*DEGtoRAD)) )
        this%log_pres = LOG10( DBLE(pressure))
        
        np = MAX( 1, this%nsza - NINTERP + 1 )
        ALLOCATE( this%densza(NINTERP, np), STAT = ierr )
        DO j = 1, np
          this%densza(:, j) = Lnk_denum(this%lg_csza(j:j+NINTERP-1))
        END DO

        np = MAX( 1, this%nvza - NINTERP + 1 )
        ALLOCATE( this%denvza(NINTERP, np), STAT = ierr )
        DO j = 1, np
          this%denvza(:, j) = Lnk_denum(this%lg_cvza(j:j+NINTERP-1))
        END DO

        np = this%npres 
        IF( np .NE. NINTERP ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, &
                                 "This code assums npres_level == NINTERP", &
                                 "O3T_lpoly_cden", zero )
           RETURN
        ENDIF

        ALLOCATE( this%denpres(np), STAT = ierr )
        this%denpres = Lnk_denum( this%log_pres )

        this%setLUTcoefDen = .TRUE.
        status = OZT_S_SUCCESS
        RETURN
      END FUNCTION O3T_lpoly_cden

      FUNCTION O3T_lpolycden_dispose( this ) RESULT (status)
        TYPE (O3T_lpoly_cden_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4) :: status
        INTEGER  :: ierr
 
        IF( this%setLUTcoefDen ) THEN
        !  write(*,*) "lg_csza=", this%lg_csza
        !  write(*,*) "lg_cvza=", this%lg_cvza
        !  write(*,*) "log_pres=", this%log_pres
        !  write(*,*) "densza=", this%densza
        !  write(*,*) "denvza=", this%denvza
        !  write(*,*) "denpres=", this%denpres
           DEALLOCATE( this%lg_csza, this%lg_cvza, this%log_pres, STAT = ierr )
           DEALLOCATE( this%densza, STAT = ierr )
           DEALLOCATE( this%denvza, STAT = ierr )
           DEALLOCATE( this%denpres, STAT = ierr )
        ENDIF
        this%setLUTcoefDen = .FALSE.
        status = OZT_S_SUCCESS
      END FUNCTION O3T_lpolycden_dispose

      FUNCTION O3T_lpoly_coef( this, cden_blk, pixGEO ) &
               RESULT (status)
        TYPE (O3T_lpoly_coef_type), INTENT( OUT ) :: this
        TYPE (O3T_lpoly_cden_type), INTENT( IN ) :: cden_blk
        TYPE (O3T_pixgeo_type), INTENT( IN ) :: pixGEO
        INTEGER (KIND=4) :: isza, isza_end, ivza, ivza_end
        INTEGER (KIND=4) :: i,j,k
        INTEGER (KIND = 4) :: status
        INTEGER  :: ierr
  
        IF( .NOT. cden_blk%setLUTcoefDen ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, &
                                 "cden block not initialized", &
                                 "O3T_lpoly_coef", zero )
           RETURN
        ENDIF

        IF( pixGEO%lg_cos_sza >  cden_blk%lg_csza(cden_blk%nsza) ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "Input sza out of range", &
                                 "O3T_lpoly_coef", zero )
           RETURN
        ENDIF

        IF( pixGEO%lg_cos_vza >  cden_blk%lg_cvza(cden_blk%nvza) ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "Input vza out of range", &
                                 "O3T_lpoly_coef", zero )
           RETURN
        ENDIF
        
        this%pixel_id = pixGEO%id
        isza = get_index( pixGEO%lg_cos_sza, cden_blk%lg_csza )
        isza_end = isza + NINTERP-1
        this%isza = isza
        this%isza_end = isza_end
        this%cthet0 = Lnk_denum(cden_blk%lg_csza(isza:isza_end),&
                                pixGEO%lg_cos_sza)
        this%cthet0 = this%cthet0/cden_blk%densza(:,isza)

        ivza = get_index( pixGEO%lg_cos_vza, cden_blk%lg_cvza )
        ivza_end = ivza + NINTERP-1
        this%ivza = ivza
        this%ivza_end = ivza_end
        this%cthet  = Lnk_denum(cden_blk%lg_cvza(ivza:ivza_end),&
                                pixGEO%lg_cos_vza)
        this%cthet  = this%cthet/cden_blk%denvza(:,ivza)

        k = 0
        DO i = 1, NINTERP
        DO j = 1, NINTERP
            k = k + 1
            this%Lnk(k) = this%cthet(j)*this%cthet0(i)
        ENDDO 
        ENDDO 

        this%pgr = Lnk_denum(cden_blk%log_pres,pixGEO%log_pt)
        this%pgr = this%pgr / cden_blk%denpres
        this%pgrp = Lnk_denum(cden_blk%log_pres,pixGEO%log_ptp)
        this%pgrp = this%pgrp/ cden_blk%denpres

        this%pcl = Lnk_denum(cden_blk%log_pres,pixGEO%log_pc)
        this%pcl = this%pcl / cden_blk%denpres
        this%pclp= Lnk_denum(cden_blk%log_pres,pixGEO%log_pcp)
        this%pclp= this%pclp/ cden_blk%denpres
        ! write(*,*) "cofpgr=", this%pgr
        ! write(*,*) "cofpcl=", this%pcl
        ! write(*,*) "cofpclp=", this%pclp

        status = OZT_S_SUCCESS
      END FUNCTION O3T_lpoly_coef

      FUNCTION get_index( x0, xarray ) RESULT( index )
        REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: xarray
        REAL (KIND = 8), INTENT(IN) :: x0
        INTEGER (KIND=4) :: nn
        INTEGER (KIND=4) :: ifoo(1)
        INTEGER (KIND=4) :: index

        nn = SIZE( xarray )
        ifoo = MINLOC( xarray, MASK = xarray - x0 >= 0.0 )
        index = ifoo(1) - NINTERP/2
        IF( index < 1 ) THEN
           index = 1
        ELSE IF( index > nn - NINTERP+1 ) THEN
           index = nn-NINTERP+1
        ENDIF
      END FUNCTION get_index

      FUNCTION Lnk_denum( xarray, x0 ) RESULT( oarray )
        REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: xarray
        REAL (KIND = 8), INTENT(IN), OPTIONAL :: x0
        REAL (KIND = 8), DIMENSION(SIZE(xarray)) :: oarray
        INTEGER (KIND=4) :: k, i, nn
        REAL (KIND = 8) :: xdenom, xnumer
        nn = SIZE( xarray )

        IF( .NOT. PRESENT(x0) ) THEN
           DO k=1,nn
             xdenom=1.0D0
             DO i=1,nn
               IF (i == k) CYCLE
               xdenom = xdenom*(xarray(k)-xarray(i))
             END DO
             oarray(k)=xdenom
           END DO
        ELSE
           DO k=1,nn
             xnumer=1.0D0
             DO i=1,nn
               IF (i == k) CYCLE
               xnumer = xnumer*(x0-xarray(i))
             END DO
             oarray(k)=xnumer
           END DO
        ENDIF
      END FUNCTION Lnk_denum

END MODULE O3T_lpolycoef_class
