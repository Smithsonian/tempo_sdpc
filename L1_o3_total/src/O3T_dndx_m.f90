!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_dndx_m
! 
! computes dndx for selected lambdas for a given table pressure, profile, 
! albedos, zenith angle and scan angle.
! O3T_dndx: performs a table lookup on quantities i0, i1, i2,
! using precomputed lagrangian polynomial coefficients
! interpolates table values to evaluate quantities for
! current thet0 and theta
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
MODULE O3T_dndx_m
    USE OMI_SMF_class
    USE O3T_dndx_class, ONLY : npresT, nlyrT, nwl_sub, log_pres, c0, c1, c2, &
                               lyrDetO3 
    USE O3T_lpolyinterp_class, ONLY : O3T_lpoly_interpPLW
    USE O3T_stnprof_class, ONLY : NLYR, stdprf, O3T_ozfraction 
    USE O3T_pixel_class, ONLY : O3T_pixgeo_type, O3T_pixcover_type
    USE O3T_lpolycoef_class, ONLY : O3T_lpoly_coef_type
    use tell_module
 
    IMPLICIT NONE
    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0
 
    PUBLIC  :: O3T_dndx
    PUBLIC  :: O3T_delnbyT
 
    CONTAINS

      FUNCTION O3T_dndx( iplow, estozn, coefs, pixGEO, pixSURF, &
           dndx, dndomega ) RESULT( errstat )

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iplow
        REAL (KIND=4), INTENT(IN) :: estozn
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(IN) :: pixSURF
        REAL (KIND=4), DIMENSION(:,:), INTENT(OUT) :: dndx
        REAL (KIND=4), DIMENSION(:), INTENT(OUT), OPTIONAL :: dndomega
        INTEGER (KIND=4) :: iwl, nwl, ii
        INTEGER (KIND=4) :: ip, iprof
        INTEGER (KIND=4) :: errstat
        REAL (KIND=4), DIMENSION( nlyrT, npresT, nwl_sub ) ::  &
             ezero, tr, sb, Lgr_npres, Lcl_npres
        REAL (KIND=4), DIMENSION( nlyrT, nwl_sub ) :: Lgr, Lcl, Ltot
        REAL (KIND=4), DIMENSION(2,nwl_sub ) :: xnval
        REAL (KIND=4), DIMENSION(2,NLYR,nwl_sub) :: dxnv
        REAL (KIND=4) :: omeglo, omeghi, ozfrac, clfrac  
        REAL (KIND=4) :: grref, clref, ref
        !REAL (KIND=4) :: q_unper
        REAL (KIND=8) :: pc, fac_pcl

        errstat = 0
        fac_pcl = 0.0  ! silence compiler warning 'fac_pcl may be used uninitialized'

        nwl = SIZE(dndx, 1)
        IF( nwl  < nwl_sub  .OR.  SIZE(dndx, 2) < NLYR  ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_dndx: input array size too small", errstat)
          RETURN
        ENDIF


        clfrac = pixSURF%clfrac 
        grref  = pixSURF%grref
        clref  = pixSURF%clref
        ref    = pixSURF%ref

        ozfrac  = O3T_ozfraction( estozn, iplow, &
             omeglo_k=omeglo, omeghi_k=omeghi )
        pc = pixGEO%pc
        IF( pc < 0.25 ) THEN
          fac_pcl = (pixGEO%log_pc-log_pres(3))/(log_pres(4)-log_pres(3))
        ENDIF

        iprof = iplow
        DO ip = 1, 2
          errstat = O3T_lpoly_interpPLW( iprof, pixGEO, coefs, &
               ezero, tr, sb  )
          IF( errstat .NE. 0 ) THEN
            call tell_error(tell_runtime_error, &
                 "O3T_dndx: table interpolation failed", errstat)
            RETURN
          ENDIF

          IF( clfrac <= 0.0 ) THEN
            Lgr_npres = ezero + grref*tr/(1.0-grref*sb)
            DO iwl = 1, nwl_sub
              Ltot(:,iwl)= real( &
                   (/(SUM( Lgr_npres(ii,:,iwl)*coefs%pgr), ii=1,nlyrT)/) &
                   , kind=4)
            ENDDO
          ELSE IF( clfrac >= 1.0 ) THEN
            Lcl_npres = ezero + clref*tr/(1.0-clref*sb)
            IF( pc < 0.25 ) THEN
              DO iwl = 1, nwl_sub
                Ltot(:,iwl) = real(fac_pcl*(Lcl_npres(:,4,iwl) &
                     - Lcl_npres(:,3,iwl))+Lcl_npres(:,3,iwl) &
                     , kind=4)
              ENDDO
            ELSE
              DO iwl = 1, nwl_sub
                Ltot(:,iwl) = real((/ (SUM(Lcl_npres(ii,:,iwl)*coefs%pcl), &
                     ii=1,nlyrT) /), kind=4)
              ENDDO
            ENDIF
          ELSE
            Lgr_npres = ezero + grref*tr/(1.0-grref*sb)
            Lcl_npres = ezero + clref*tr/(1.0-clref*sb)
            DO iwl = 1, nwl_sub
              Lgr(:,iwl)=real((/(SUM(Lgr_npres(ii,:,iwl)*coefs%pgr), &
                   ii=1,nlyrT)/), kind=4)
            ENDDO
            IF( pc < 0.25 ) THEN
              DO iwl = 1, nwl_sub
                Lcl(:,iwl) = real(fac_pcl*(Lcl_npres(:,4,iwl) &
                     - Lcl_npres(:,3,iwl))+Lcl_npres(:,3,iwl), &
                     kind=4)
              ENDDO
            ELSE
              DO iwl = 1, nwl_sub
                Lcl(:,iwl) = real((/ (SUM(Lcl_npres(ii,:,iwl)*coefs%pcl), &
                     ii=1,nlyrT) /), kind=4)
              ENDDO
            ENDIF
            Ltot = clfrac*Lcl + (1.0-clfrac)*Lgr
          ENDIF
          IF( PRESENT( dndomega ) ) THEN
            xnval(ip,:)  = -100.0*LOG10(Ltot(1,:))
          ENDIF
          IF( lyrDetO3(1) > 0.0 ) THEN
            DO iwl = 1, nwl_sub
              dxnv(ip,:,iwl) = -100.0*LOG10(Ltot(2:nlyrT,iwl)/Ltot(1,iwl))
            ENDDO
          ELSE 
            DO iwl = 1, nwl_sub
              dxnv(ip,:,iwl) = -100.0*LOG10(Ltot(2:nlyrT,iwl)/Ltot(1,iwl)) &
                   / (0.1*stdprf(1:NLYR,iprof))
            ENDDO
          ENDIF
          iprof = iplow + ip
        ENDDO
        IF( PRESENT( dndomega ) ) THEN
          IF( SIZE(dndomega) < nwl ) THEN
            call tell_error(tell_runtime_error, &
                 "O3T_dndx: dndomega size too small", errstat)
            RETURN
          ENDIF
          dndomega(1:nwl_sub) = (xnval(2,:) - xnval(1,:))/(omeghi - omeglo)
        ENDIF
        DO iwl = 1, nwl_sub
          dndx(iwl,:) = dxnv(2,:,iwl)*ozfrac + dxnv(1,:,iwl)*(1.0-ozfrac)
        ENDDO

        ! set dndx for additional wavelength  to zero
        IF( nwl_sub < nwl ) THEN
          dndx(nwl_sub+1:nwl,:)   = 0.0
          IF( PRESENT( dndomega ) ) dndomega(nwl_sub+1:nwl) = 0.0
        ENDIF

      END FUNCTION O3T_dndx




!     FUNCTION O3T_delnbyT computes the Jacobian dN/dT for each of 11 Umkehr 
!        layers and all the selected wavelengths using a chain rule relating 
!        dndt to dndx.
!        the chain rule formula is: dN/dT = dN/dx * dx/da * da/dT,
!        where x is layer ozone amount and a is absorption coefficient.
!        this is used as:  
!        delN = dN/dx * x/a * ( c1*(T1 - T0)+c2*(T1-T0)*(T1+T0) 
!        to calculate the change in n-value at each wavelength due to
!        the temperature difference between first guess and a priori.
!        x is taken from the standard profile, dndx is from table
!        lookup (see O3T_dndx), and a and its temperature
!        dependence are taken from a = c0 + c1 * T + c2 * T ** 2.
!        the coefficients, c0, c1, and c2 are slit averages
!        calculated for each of the wavelengths less and 360nm and
!        each layer of the standard profile. they are stored by block
!        data for each spacecraft and transferred appropriately
!        by subroutine start to the variables c0, c1, and c2 for
!        use here by subroutine O3T_delnbyT

      FUNCTION O3T_delnbyT( dndx, fgprf, fgtmp, aprftm, delnT, dNdT ) &
                            RESULT( errstat )

        implicit none

        REAL (KIND=4), DIMENSION(:,:), INTENT(IN) :: dndx
        REAL( KIND=4 ), DIMENSION(:), INTENT(IN) :: fgprf, fgtmp, aprftm
        REAL (KIND=4), DIMENSION(:,:), INTENT(OUT) :: delnT
        REAL (KIND=4), DIMENSION(:), OPTIONAL, INTENT(OUT) :: dNdT !Oz weighted
        REAL( KIND=4 ), DIMENSION(SIZE(aprftm)) :: ctmp, delt, deltsq, t2
        INTEGER (KIND=4) :: iwl, nwl, errstat

        errstat = 0

        IF( ANY( SHAPE(delnT) .NE. SHAPE(dndx)) ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_delnbyT: Shapes dndx & delnT differ", errstat)
           RETURN
        ENDIF

        IF( SIZE( fgprf ) .NE. SIZE( fgtmp ) .OR. &
            SIZE( fgprf ) .NE. SIZE( aprftm )  ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_delnbyT: profile array sizes not equal", errstat)
           RETURN
        ENDIF

        nwl = SIZE( dndx, 1 )
        IF( SIZE( fgprf ) .NE. SIZE( dndx, 2 )  ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_delnbyT: layersizes not equal", errstat)
           RETURN
        ENDIF

        ctmp   = aprftm - 273.13
        delt   = aprftm - fgtmp
        deltsq = (aprftm - 273.13)**2 - (fgtmp - 273.13)**2

        DO iwl = 1, nwl_sub
          IF( c0(iwl) > 1.0e-6 ) THEN
             delnT(iwl,:) = (dndx(iwl,:)*fgprf &
                          / (c0(iwl) + c1(iwl)*ctmp + c2(iwl)*ctmp**2.0)) &
                          * (c1(iwl) * delt + c2(iwl) * deltsq)
          ELSE
             delnT(iwl,:) = 0.0
          ENDIF
        ENDDO

        ! set delnT for additional wavelength  to zero
        IF( nwl_sub < nwl ) THEN
           delnT(nwl_sub+1:nwl,:) = 0.0 
        ENDIF

        IF( PRESENT( dNdT ) ) THEN
           IF( SIZE( dNdT ) < nwl ) THEN 
          call tell_error(tell_runtime_error, &
               "O3T_delnbyT: dNdT size too small", errstat)
              RETURN
           ENDIF

           t2(:) = (aprftm(:) - 273.13) + (fgtmp(:) - 273.13)
           DO iwl = 1, nwl_sub
             IF( c0(iwl) > 1.0e-6 ) THEN !!Remember N value is -100*Log[rad/Irr]
             !! dNdT(iwl) = SUM( dndx(iwl,:) * fgprf(:)                 &
             !!   *( fgprf(:)/(c0(iwl)+c1(iwl)*ctmp(:)+c2(iwl)*ctmp(:)**2.0) ) &
             !!   *( c1(iwl) + c2(iwl)*t2(:) ) ) / SUM( fgprf(:) )
                dNdT(iwl) = SUM( dndx(iwl,:)                            &
                  *( fgprf(:)/(c0(iwl)+c1(iwl)*ctmp(:)+c2(iwl)*ctmp(:)**2.0) ) &
                  *( c1(iwl) + c2(iwl)*t2(:) ) ) 
             ELSE
                dNdT(iwl) = 0.0
             ENDIF
           ENDDO
           
           IF( nwl_sub < nwl ) THEN
              dNdT(nwl_sub+1:nwl) = 0.0 
           ENDIF
        ENDIF

      END FUNCTION O3T_delnbyT
END MODULE O3T_dndx_m
