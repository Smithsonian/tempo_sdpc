!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_class
!
! This module contains the functions that implement the algorithm steps in
! the process of successive improvement of ozone retrieval.
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
!
!!Return
! None
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
MODULE O3T_class
    USE OMI_SMF_class
    USE O3T_stnprof_class
    USE O3T_apriori_class
    USE O3T_lpolycoef_class, ONLY : O3T_lpoly_coef_type
    USE O3T_pixel_class, ONLY : O3T_pixgeo_type, O3T_pixcover_type, &
                                DEGtoRAD, RADtoDEG
    USE O3T_nval_class, ONLY: nwl_com, nvRRS
    USE O3T_dndx_class, ONLY: nwl_sub
    USE O3T_iztrsb_m, ONLY : O3T_iztrsb
    USE O3T_dndx_m, ONLY : O3T_dndx, O3T_delnbyT
    USE O3T_const, ONLY : swthrsh
    USE UTIL_tools_class

    IMPLICIT NONE
    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0
    INTEGER (KIND=4), PARAMETER, PRIVATE :: four = 4
    REAL (KIND=4), PARAMETER, PRIVATE :: EPSILON10 = 1.0E-5
    REAL (KIND=8), PARAMETER, PRIVATE :: PI = 3.1415926535897932384626433832795d0

    PUBLIC  :: O3T_nvalm
    PUBLIC  :: O3T_step1
    PUBLIC  :: O3T_step2!, O3T_step2_omps
    PUBLIC  :: O3T_step3
    PUBLIC  :: O3T_calcRefl
    PUBLIC  :: O3T_oznot
    PUBLIC  :: O3T_residue
    PUBLIC  :: O3T_stp2oz
    PUBLIC  :: O3T_getshp
    PUBLIC  :: O3T_resadj
    PUBLIC  :: O3T_blwcld
    PUBLIC  :: O3T_glnchk
    PUBLIC  :: O3T_rad
    PUBLIC  :: O3T_descendQ
    PUBLIC  :: O3T_getSatName
    PRIVATE :: O3T_dndomega
!    PRIVATE :: initRefl
    PRIVATE :: O3T_rcf1
    PRIVATE :: O3T_rcf2

    CONTAINS
!!
!!  Before xnvalm is calculated with the following way. Both radiance and
!!  irradiance are interpolated to the nval table wavelengths. The
!!  interpolations are done in the L1B reader function L1Bri_getLineWL.
!!  Then the rad/irr ratios are calculated and the negative log10 of the
!!  ratios are xnvalm. Now instead of calling L1Bri_getLineWL, L1Bri_getLine
!!  is called. This returns the orginal irradiance and radiance values
!!  stored in the L1B file. Then thses values are passed into O3T_xnvalm
!!  together with wl_com, the wavelength in obtained the nval LUT. Then first
!!  irradiance value are interpolated to the radiance wavelength. rad/irr
!!  can then be calculated because they are at the same wavelength. Next the
!!  ratio is interpolated to the tabulated wavelength (wl_com, the nval LUT
!!  wavelength). Then the negative log of the ratio at the table wavelength
!!  is calculated and returned at the xnvalm.
!!
      FUNCTION O3T_nvalm( irrWl, irradiance, &
           radWl, radiance,   &
           wl_com, xnvalm ) RESULT( errstat )

        use tell_module

        implicit none

        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: wl_com, &
             irrWl, irradiance, &
             radWl, radiance
        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: xnvalm
        REAL (KIND=4), DIMENSION(SIZE(radWl)) :: Rtoa
        REAL (KIND=4) :: irrInterp, RtoaInterp, frac, dwl
        INTEGER :: nwl_l, nwl_rad, nwl_irr
        INTEGER :: iwl, il, ih
        INTEGER :: errstat

        errstat = 0

        nwl_l = SIZE( wl_com )

        IF( SIZE( xnvalm ) /= nwl_l .OR. nwl_l /= nwl_com ) THEN
          call tell_error (tell_runtime_error, &
               "O3T_nvalm: xnvalm and wl_com array sizes not match", &
               errstat)
          RETURN
        ENDIF

        nwl_irr = SIZE( irrWl )
        IF( nwl_irr /= SIZE( irradiance) ) THEN
          call tell_error (tell_runtime_error, &
               "O3T_nvalm: irrWl and irradiance array sizes not match", &
               errstat)
          RETURN
        ENDIF

        nwl_rad = SIZE( radWl )
        IF( nwl_rad /= SIZE( radiance) ) THEN
          call tell_error (tell_runtime_error, &
               "O3T_nvalm: radWl and radiance array sizes not match", &
               errstat)
          RETURN
        ENDIF

        !! first get the irradiance value at radiance wavelength
        DO iwl = 1, nwl_rad
          il = iwl
          il = hunt( irrWl(:), radWl(iwl), il )
          IF( il == 0 .OR. il == nwl_irr ) THEN
            !! no irradiance data at this radiance wl, set it 0.
            Rtoa(iwl)  = 0.0
          ELSE
            ih = il+1
            dwl = irrWl(ih) - irrWl(il)
            IF( dwl < EPSILON10 ) THEN
              call tell_error (tell_runtime_error, &
                   "O3T_nvalm: irrWl grid spacing too small", errstat)
              RETURN
            ENDIF

            frac = ( radWl(iwl) - irrWl(il) ) / dwl
            IF( ABS( frac ) < EPSILON10 ) THEN
              irrInterp = irradiance(il)
            ELSE IF( ABS( frac -1.0 ) < EPSILON10 ) THEN
              irrInterp = irradiance(ih)
            ELSE
              irrInterp = (1.0-frac)*irradiance(il) + frac*irradiance(ih)
            ENDIF

            IF( irrInterp > EPSILON10 .AND. radiance(iwl) > EPSILON10 ) THEN
              Rtoa(iwl) = radiance(iwl)/irrInterp
            ELSE
              Rtoa(iwl) = 0.0
            ENDIF
          ENDIF
        ENDDO

        !! second get the N-value at fixed wavelength grid wl_com
        DO iwl = 1, nwl_l
          il = iwl
          il = hunt( radWl(:), wl_com(iwl), il )
          IF( il == 0 .OR. il == nwl_rad ) THEN
            !! no radiance data at this table wl, set it 0.
            xnvalm(iwl)  = 0.0
          ELSE
            ih = il+1
            dwl = radWl(ih) - radWl(il)

            IF( dwl < EPSILON10 ) THEN
              call tell_error (tell_runtime_error, &
                   "O3T_nvalm: wl_com grid spacing too small", errstat)
              RETURN
            ENDIF

            frac = ( wl_com(iwl) - radWl(il) ) / dwl
            IF( ABS( frac ) < EPSILON10 ) THEN
              RtoaInterp = Rtoa(il)
            ELSE IF( ABS( frac -1.0 ) < EPSILON10 ) THEN
              RtoaInterp = Rtoa(ih)
            ELSE
              RtoaInterp = (1.0-frac)*Rtoa(il) + frac*Rtoa(ih)
            ENDIF
            IF( RtoaInterp > EPSILON10 ) THEN
              xnvalm(iwl) = -ALOG10( RtoaInterp )
            ELSE
              xnvalm(iwl) = 0.0
            ENDIF
          ENDIF
        ENDDO

        RETURN
      END FUNCTION O3T_nvalm



!!Step 1:
!!the algorithm uses a pair of wavelengths by deriving reflectivity from 331 nm
!!and ozone from 318nm. This is an iterative process which starts with an
!!initial guess ozone, and converges when the reflectivity and the ozone
!!become consistent, meaning that the measured radiances for the pair are equal
!!(within limit) to the calculated radiances based on the reflectivity and the
!!a standard profile of the total ozone loading. The different pair (331&360)
!!may be used under high ozone and high solar zenith angles.
!!
      FUNCTION O3T_step1( iwl_oz, iwl_refl_l, iwl_refl_h, xnvalm, pixGEO, &
           pixSURF, coefs, iwl_refl, iplow, guesoz, stp1oz, &
           res_stp1, dndomega_t, dndr, algflg, absrfl, &
           maxitr, stp1oz_valid, skipit, doz_limit )  &
           RESULT( errstat )

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iwl_refl_l, iwl_refl_h
        INTEGER (KIND=4), INTENT(INOUT) :: iwl_oz
        TYPE (O3T_pixgeo_type), INTENT(INOUT) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: xnvalm
        REAL (KIND=4), OPTIONAL, INTENT(IN) :: doz_limit
        INTEGER (KIND=4), INTENT(OUT) :: iplow
        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: res_stp1
        REAL (KIND=4), INTENT(INOUT) :: guesoz
        REAL (KIND=4), INTENT(OUT) :: stp1oz
        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: dndomega_t, dndr
        INTEGER (KIND=4), INTENT(OUT) :: iwl_refl
        INTEGER (KIND=1), INTENT(INOUT) :: algflg
        LOGICAL, INTENT(INOUT) :: skipit
        LOGICAL, INTENT(OUT) :: absrfl, maxitr, stp1oz_valid
        INTEGER (KIND=4) :: iphigh !, iprof
        INTEGER (KIND=4) :: ilat
        INTEGER (KIND=4) :: iter, itermax
        LOGICAL (KIND=4) :: first_call
        REAL (KIND=4) :: ozonin, estozn
        REAL (KIND=4) :: dndoz_r, dndoz_o, doz_limit_l! , dr
        INTEGER :: errstat

        errstat = 0

        IF( SIZE( xnvalm ) < nwl_com ) THEN
          call tell_error(tell_runtime_error, "array size too small", errstat)
          RETURN
        ENDIF

        IF( SIZE( xnvalm ) .NE. SIZE( dndomega_t ) .OR. &
             SIZE( xnvalm ) .NE. SIZE( dndr ) .OR. &
             SIZE( xnvalm ) .NE. SIZE( res_stp1 ) ) THEN
          call tell_error(tell_runtime_error, "array size not equal", errstat)
          RETURN
        ENDIF

        !! determine latitude band based in pixel center latitude.
        IF( ABS(pixGEO%lat) <  30.0 )  THEN
          ilat = 1
        ELSE IF( ABS(pixGEO%lat) <=  60.0 ) THEN
          ilat = 2
        ELSE
          ilat = 3
        ENDIF

        !! set initial ozone guess, if input guess ozone out of valid range
        !! use guess that is set for the latitude band.
        IF( guesoz < 50.0 .OR. guesoz > 600.0 ) THEN
          IF( ABS(pixGEO%lat) <=  45.0 ) THEN
            guesoz = 260
          ELSE IF(ABS(pixGEO%lat) <= 75.0 ) THEN
            guesoz = 340
          ELSE
            guesoz = 360
          ENDIF
        ENDIF

        itermax = 8
        ozonin = guesoz
        iplow = O3T_stnprof_idxf( ozonin, ilat )
        absrfl = .TRUE.
        maxitr = .FALSE.
        stp1oz_valid = .TRUE.

        first_call = .TRUE.
        iwl_refl   = iwl_refl_l
        IF( PRESENT( doz_limit ) ) THEN !! doz_limit_l: convergent threshold
          doz_limit_l = doz_limit
        ELSE
          doz_limit_l = 1.0
        ENDIF

        pixSURF%iwl_refl_l = iwl_refl_l
        pixSURF%iwl_refl_h = iwl_refl_h
        pixSURF%iwl_glint  = iwl_refl_l
        !! Iteration to find the step 1 ozone.
        algflg = 1
        DO iter = 1, itermax
          errstat = O3T_calcRefl( iplow, ozonin, iwl_refl, xnvalm, coefs, &
               pixGEO, pixSURF, dndoz_r, first_call )
          IF( errstat .NE. 0 ) THEN
            call tell_error(tell_runtime_error, "O3T_calcRefl error", errstat)
            skipit = .TRUE.
            RETURN
          ENDIF

          errstat = O3T_oznot( iwl_oz, xnvalm, coefs, pixGEO, pixSURF, &
               iplow, iphigh, estozn, dndoz_o, skipit )
          IF( errstat .NE. 0 ) THEN
            call tell_error(tell_runtime_error, "O3T_oznot error", errstat)
            RETURN
          ENDIF
          !!        IF( estozn <= 0.0 .OR. estozn > 900.0 .OR. &
          !!            pixSURF%ref < -0.25 .OR.  pixSURF%ref > 1.25 )  THEN
          IF( estozn <= 0.0 .OR. estozn > 900.0  ) THEN
            stp1oz = estozn
            skipit = .TRUE.                !! check to see if step 1 ozone
            stp1oz_valid = .FALSE.         !! is in the valid range.
          ENDIF
          IF( skipit ) RETURN

          IF( iter == 1 ) THEN                              !!This part decides
            ozonin = estozn                                !!whether C pair is
            IF( (dndoz_o-dndoz_r)/dndoz_r < swthrsh ) THEN !!used. Note that
              ozonin = guesoz                             !!dndoz_o & dndoz_r
              iwl_refl = iwl_refl_h                       !!actually depends
              iwl_oz   = iwl_refl_l                       !!on ozonin,
              first_call = .TRUE.                         !!consequently that
              iplow = O3T_stnprof_idxf( ozonin, ilat )    !!the switch may
              absrfl = .FALSE.                            !!depends on the
              algflg = 3                                  !!initial O3 guess.
            ENDIF
            CYCLE                                          !!next iteration
          ENDIF

          IF( ABS( estozn - ozonin ) > doz_limit_l ) THEN
            ozonin = estozn              !! not converged yet
          ELSE
            EXIT                         !! iteration converged,
          ENDIF                           !! exit the do-while loop.
        ENDDO

        IF( iter > itermax ) THEN         !! check to if maximum iteration
          skipit = .TRUE.                !! is exceeded. If yes, bad ozone,
          maxitr = .TRUE.                !! return to main code.
          stp1oz = estozn
          RETURN
        ENDIF

        !! Calculate step 1 residue and other useful quantities. Note that
        !! dndomega_t (dN/dOz) is calculated in O3T_residue, which uses the
        !! the N values from the NVAL table. This dndomega_t is used
        !! through out the program. Note that a similar quantity dndomega
        !! is calculated in O3T_step2 using the function O3T_dndx, which
        !! uses the N values from the DNDX table. In principle dndomega_t
        !! and dndomega are the same. However if NVAL table includes the LRRS
        !! ring, while DNDX table does not, so dndomega_t and dndomega will
        !! have different values. if a no-ring NVAL table is used, the
        !! two quanities will have very close values (but not idential values)
        !! due to numerical round-off and different N values in the NVAL and
        !! DNDX tables.

        stp1oz = estozn
        errstat = O3T_residue( iplow, stp1oz, xnvalm, coefs, pixGEO, &
             pixSURF, res_stp1, dndomega_t, dndr )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step1: determine step 1 residue failed", errstat)
          RETURN
        ENDIF

        !! calculate reflectivity at 360 nm.
        pixSURF%ref360 = pixSURF%ref &
             + (res_stp1(iwl_refl_h)-res_stp1(iwl_refl_l)) &
             / dndr(iwl_refl_h)
        !! Compute Radiative Cloud Fraction based on reflectivity 360
      END FUNCTION O3T_step1




!! Step 2:
!! Ozone and temperature climatologies are applied at all levels to account
!! for seasonal and latitudinal variations in profile shape.
!!
      FUNCTION O3T_step2( iwl_oz, iwl_refl, iplow, jday, pixGEO, pixSURF, &
           coefs, stp1oz, res_stp1, dndomega_t, dndr, absrfl, &
           stp2oz, res_stp2, dr, stp2prf, eff, aprfoz, &
           stp2oz_valid, skipit, dNdT )  RESULT( errstat )

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iwl_oz, iwl_refl
        TYPE (O3T_pixgeo_type), INTENT(INOUT) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        INTEGER (KIND=4), INTENT(IN) :: iplow, jday
        REAL (KIND=4), INTENT(IN) :: stp1oz
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: res_stp1, dndomega_t, dndr
        REAL (KIND=4), DIMENSION(:), OPTIONAL, INTENT(OUT) :: dNdT
        LOGICAL, INTENT(IN) :: absrfl
        REAL (KIND=4), INTENT(OUT) :: stp2oz, dr
        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: res_stp2
        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: stp2prf, eff, aprfoz
        LOGICAL, INTENT(OUT) :: stp2oz_valid
        LOGICAL, INTENT(INOUT) :: skipit

        !INTEGER (KIND=4) :: iprof
        INTEGER :: errstat

        REAL (KIND=4), DIMENSION(nwl_com) :: dndomega
        REAL (KIND=4), DIMENSION(NLYR) :: fgprf, dxdomega
        REAL (KIND=4), DIMENSION(NLYR) :: aprftm, delshp
        REAL (KIND=4), DIMENSION(nwl_com, NLYR) :: dndx, delnT

        errstat = 0

        IF( SIZE( res_stp1 ) .NE. SIZE( res_stp2 ) ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: array size not equal", errstat)
          RETURN
        ENDIF

        IF( SIZE( stp2prf ) < NLYR .OR. SIZE( eff ) < NLYR )  THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: array size too small", errstat)
          RETURN
        ENDIF

        !! The quantity dndomega calculated from O3T_dndx is not
        !! used at all in the function. It is discarded. While
        !! dndomega_t passed into function is used instead.

        errstat = O3T_dndx( iplow, stp1oz, coefs, pixGEO, pixSURF, &
             dndx, dndomega )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: calculate DNDX failed", errstat)
          RETURN
        ENDIF

        ! determine first guess profile and dxdomega
        errstat = O3T_stnprof_1stG( stp1oz, iplow, fgprf, dxdomega )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: determine first guess profile failed", errstat)
          RETURN
        ENDIF

        ! determine apriori temperature profile
        errstat = O3T_apriori_prf(pixGEO%lat, jday, aprftm )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: determine apriori temperature profile failed", &
               errstat)
          RETURN
        ENDIF

        ! calculate temperature adjustments
        errstat = O3T_delnbyT( dndx, fgprf, fgtmp, aprftm, delnT, dNdT )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: calculate temperature adjustment failed", errstat)
          RETURN
        ENDIF

        ! determine apriori profile shape
        errstat =  O3T_getshp( pixGEO%lat, jday, stp1oz, fgprf, aprfoz, &
             delshp )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: determine apriori profile shape failed", errstat)
          RETURN
        ENDIF

        ! calculate step 2 ozone
        stp2oz = stp1oz    !initialized to step 1 ozone
        stp2oz_valid = .FALSE.

        errstat = O3T_stp2oz( iwl_oz, iwl_refl, stp1oz, fgprf, &
             dndx, delnT, dndomega_t, dndr, dxdomega, &
             delshp, absrfl, dr, stp2prf, eff, stp2oz )
!             fgtmp, dndx, delnT, dndomega_t, dndr, dxdomega, &
!             aprftm, delshp, absrfl, dr, stp2prf, eff, stp2oz )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: determine step 2 ozone failed", errstat)
          RETURN
        ENDIF
        IF( stp2oz > 0.0 .AND. stp2oz < 900.0 )  THEN
          stp2oz_valid = .TRUE.
        ELSE
          skipit = .TRUE.  !! temporay disabled for comparison wtih v806
        ENDIF

        errstat = O3T_resadj( iwl_refl, res_stp1, dndr, dndx, delnT, &
             fgprf, stp2prf, dr, res_stp2 )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step2: determine step 2 residue failed", errstat)
          RETURN
        ENDIF
        pixSURF%ref    =  pixSURF%ref + dr

      END FUNCTION O3T_step2




!      FUNCTION O3T_step2_omps(iwl_oz, iwl_refl, iplow, &
!                              pixGEO, pixSURF, aprftm, aprfoz, &
!                              coefs, stp1oz, res_stp1, dndomega_t, dndr,  &
!                              absrfl, stp2oz, res_stp2, dr, stp2prf, eff, &
!                              stp2oz_valid, skipit, dNdT)  RESULT( status )
!
!        ! Description:
!        !     It is a clone of O3T_step2(), except a caller now must supply
!        !     apriori temperature profile and apiori ozone profile (namely
!        !     aprftm(:) and aprfoz(:) respectively) to this function as
!        !     input arguments. They are for a pixel of interest.
!        !
!        !     As a result, we no longer have use for the jday input, nor the
!        !     need to call the private function O3T_getshp(). The shape is
!        !     calculated simply by:  delshp(:) = aprfoz(:) - fgprf(:)
!        !
!        ! [JYLI @ Sun Feb 13 09:37:59 EST 2011]
!
!        INTEGER (KIND=4), INTENT(IN) :: iwl_oz, iwl_refl
!        TYPE (O3T_pixgeo_type), INTENT(INOUT) :: pixGEO
!        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
!        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: aprftm
!        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: aprfoz
!        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
!        INTEGER (KIND=4), INTENT(IN) :: iplow
!        REAL (KIND=4), INTENT(IN) :: stp1oz
!        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: res_stp1, dndomega_t, dndr
!        REAL (KIND=4), DIMENSION(:), OPTIONAL, INTENT(OUT) :: dNdT
!        LOGICAL, INTENT(IN) :: absrfl
!        REAL (KIND=4), INTENT(OUT) :: stp2oz, dr
!        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: res_stp2
!        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: stp2prf, eff
!        LOGICAL, INTENT(OUT) :: stp2oz_valid
!        LOGICAL, INTENT(INOUT) :: skipit
!
!        INTEGER (KIND=4) :: iprof
!        INTEGER :: ierr, status
!
!        REAL (KIND=4), DIMENSION(nwl_com) :: dndomega
!        REAL (KIND=4), DIMENSION(NLYR) :: fgprf, dxdomega
!        REAL (KIND=4), DIMENSION(NLYR) :: delshp
!        REAL (KIND=4), DIMENSION(nwl_com, NLYR) :: dndx, delnT
!
!
!        IF( SIZE( res_stp1 ) .NE. SIZE( res_stp2 ) ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
!                                "res_stp1 and res_stp2 array size not equal", &
!                                "O3T_step2", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!
!        IF( SIZE( stp2prf ) < NLYR .OR. SIZE( eff ) < NLYR )  THEN
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
!                                 "stp2prf[] and eff array size too small", &
!                                 "O3T_step2", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!
!        IF( SIZE( aprftm ) < NLYR )  THEN
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
!                                 "aprftm[] array size too small", &
!                                 "O3T_step2", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!
!        IF( SIZE( aprfoz ) < NLYR )  THEN
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
!                                 "aprfoz[] array size too small", &
!                                 "O3T_step2", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!
!        !! The quantity dndomega calculated from O3T_dndx is not
!        !! used at all in the function. It is discarded. While
!        !! dndomega_t passed into function is used instead.
!
!        status = O3T_dndx( iplow, stp1oz, coefs, pixGEO, pixSURF, &
!                           dndx, dndomega )
!        IF( status .NE. OZT_S_SUCCESS ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
!                                 "Calculated DNDX failed", &
!                                 "O3T_step2", zero )
!           RETURN
!        ENDIF
!
!        ! determine first guess profile and dxdomega
!        status = O3T_stnprof_1stG( stp1oz, iplow, fgprf, dxdomega )
!        IF( status .NE. OZT_S_SUCCESS ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
!                    "O3T_stnprof_1stG:determine first guess profile failed", &
!                    "O3T_step2", zero )
!           RETURN
!        ENDIF
!
!        ! calculate temperature adjustments
!        status = O3T_delnbyT( dndx, fgprf, fgtmp, aprftm, delnT, dNdT )
!        IF( status .NE. OZT_S_SUCCESS ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
!                     "O3T_delnbyT: calculate temperature adjustment failed", &
!                     "O3T_step2", zero )
!           RETURN
!        ENDIF
!
!        ! determine apriori profile shape
!        IF( ANY(aprfoz(1:NLYR) < 0.0) ) THEN ! trap negative O3 value
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
!                     "Unable to determine apriori profile shape", &
!                     "O3T_step2", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!        delshp(1:NLYR) = aprfoz(1:NLYR) - fgprf(1:NLYR)
!
!        ! calculate step 2 ozone
!        stp2oz = stp1oz    !initialized to step 1 ozone
!        stp2oz_valid = .FALSE.
!
!        status = O3T_stp2oz( iwl_oz, iwl_refl, stp1oz, fgprf, fgtmp, &
!                            dndx, delnT, dndomega_t, dndr, dxdomega, &
!                            aprftm, delshp, absrfl, dr, stp2prf, eff, stp2oz )
!        IF( status .NE. OZT_S_SUCCESS ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
!                     "O3T_stp2oz: determine step 2 ozone failed", &
!                     "O3T_step2", zero )
!           RETURN
!        ENDIF
!        IF( stp2oz > 0.0 .AND. stp2oz < 900.0 )  THEN
!           stp2oz_valid = .TRUE.
!        ELSE
!           skipit = .TRUE.  !! temporay disabled for comparison wtih v806
!        ENDIF
!
!        status = O3T_resadj( iwl_refl, res_stp1, dndr, dndx, delnT, &
!                             fgprf, stp2prf, dr, res_stp2 )
!        IF( status .NE. OZT_S_SUCCESS ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
!                     "O3T_resadj: determine step 2 residue failed", &
!                     "O3T_step2", zero )
!           RETURN
!        ENDIF
!        pixSURF%ref =  pixSURF%ref + dr
!
!      END FUNCTION O3T_step2_omps

!! Step 3:
!! Correct step 2 ozone for wavelength dependence effects, such as
!! tropospheric aerosol, sun glint, and local upper level profile
!! shape effects.
!!
      FUNCTION O3T_step3( iglnt, irflo, imixr, pixGEO, stp2oz, f313, &
           f360, res_stp2, dndomega_t, &
           stp3oz, res_stp3, pathl, algflg ) &
           RESULT( errstat )

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iglnt, irflo, imixr
        TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
        REAL (KIND=4), INTENT(IN) :: stp2oz, f313
        REAL (KIND=4), DIMENSION(3), INTENT(IN) :: f360
        REAL (KIND=4), DIMENSION(:), INTENT(IN) ::  res_stp2, dndomega_t
        REAL (KIND=4), INTENT(OUT) :: stp3oz, pathl
        REAL (KIND=4), DIMENSION(:), INTENT(OUT) ::  res_stp3
        INTEGER (KIND=1), INTENT(INOUT) :: algflg

        INTEGER (KIND=4) :: irfhi
        REAL (KIND=4) :: ozpath, aiadj
        INTEGER (KIND=4) :: errstat

        errstat = 0

        IF( SIZE( res_stp2 ) .NE. SIZE( res_stp3 ) .OR. &
             SIZE( res_stp2 ) .NE. SIZE(dndomega_t) ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_step3: array size not equal", errstat)
          RETURN
        ENDIF

        pathl  = real(1.0/pixGEO%cos_sza + 1.0/pixGEO%cos_vza , kind=4)
        ozpath = pathl * stp2oz

        stp3oz = stp2oz
        IF( ozpath < 2000.0 .AND. algflg == 1 ) THEN
          aiadj=f360(1)+f360(2)*pathl+f360(3)*pathl**2.0
          stp3oz = stp2oz + aiadj* (res_stp2(iglnt)/100.0 )*stp2oz
        ELSE IF( ozpath >= 2000.0 .AND. algflg == 1 ) THEN
          stp3oz = stp2oz + f313 * (res_stp2(imixr)-res_stp2(irflo))
          algflg = 2
        ENDIF
        IF( algflg == 3 ) THEN
          irfhi = iglnt
          stp3oz = stp2oz + (  res_stp2(irflo) - res_stp2(irfhi)  ) / &
               (dndomega_t(irflo) - dndomega_t(irfhi))
        ENDIF

        !   adjust residues for step 3 ozone change
        res_stp3 = res_stp2 - dndomega_t*(stp3oz-stp2oz)

      END FUNCTION O3T_step3

!!Description:
! FUNCTION O3T_oznot
!  This function brackets measured n-value with table n-values, computes
!  total ozone by linear interpolation using the sensitivity dndomega.
!!Input Parameters:
!   iwl_oz   : the wavelength index used for ozone retrieval
!   xnvalm   : the measured nvalues for the all the wavelength
!   coefs   : Lagrangian interpolation coefficients associated with
!             pixGEO
!   pixGEO  : the geometry inforamtion for the pixel to the output
!             parameters are associated with. pixGEO has angular and
!             terrain and cloud pressure info.
!   pixSURF : the surface, cloud reflectance, and cloud fraction
!   iplow    : the index of lower bound ozone profile
!   iphigh   : the index of upper bound ozone profile
!!Output Parameters:
!   iplow    : the index of lower bound ozone profile
!   iphigh   : the index of upper bound ozone profile (equal to iplow+1)
!   estozn   : the output ozone amount in D.U. (normally between iplow, iplow+1)
!!Return
!   satus : The SMF stauts
!!References and Credits
! Written by
! Kai Yang
! University of Maryland Baltimore County
!!END Description:

      FUNCTION O3T_oznot( iwl_oz, xnvalm, coefs, pixGEO, pixSURF, &
                          iplow, iphigh, estozn, dndoz_o, skipit ) &
                          RESULT( errstat )

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iwl_oz
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: xnvalm
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(IN) :: pixSURF
        INTEGER (KIND=4), INTENT(INOUT) :: iplow, iphigh
        LOGICAL, INTENT(OUT) :: skipit
        REAL (KIND=4), INTENT(OUT) :: estozn, dndoz_o
        REAL (KIND=4) :: xlown, xhighn
        INTEGER (KIND=4) :: errstat
        INTEGER (KIND=4) :: iprof
        REAL (KIND=4) :: ezgr, tgr, sbgr, knbgr, ezcl, tcl, sbcl, knbcl
        REAL (KIND=4) :: radgr, rsbgr, radcl, rsbcl,  xnvalc
        REAL (KIND=4) :: omeglo, omeghi
        LOGICAL :: lb_found, ub_found
        CHARACTER (LEN =255) :: msg

        errstat = 0

        iprof = iplow
        lb_found = .FALSE.
        ub_found = .FALSE.
        radcl = 0.0
        radgr = 0.0

        DO
          errstat = O3T_iztrsb( nvRRS, iwl_oz, iprof, pixGEO, coefs,  &
                               ezgr, tgr, sbgr, knbgr, ezcl, tcl, sbcl, knbcl )
          IF( errstat .NE. 0 ) THEN
            call tell_error(tell_runtime_error, &
                 "O3T_oznot: O3T_iztrsb error", errstat)
             skipit = .TRUE.
             RETURN
          ENDIF

          IF( pixSURF%clfrac <= 0.0 ) THEN
             rsbgr = pixSURF%grref/(1.0-pixSURF%grref*sbgr )
             radgr = ezgr + tgr*rsbgr*(1.0+knbgr*rsbgr)
             xnvalc = -LOG10( radgr )
          ELSE IF( pixSURF%clfrac >= 1.0 ) THEN
             rsbcl = pixSURF%clref/(1.0-pixSURF%clref*sbcl )
             radcl = ezcl + rsbcl*tcl*(1.0+knbcl*rsbcl)
             xnvalc = -LOG10( radcl )
          ELSE
             rsbgr = pixSURF%grref/(1.0-pixSURF%grref*sbgr )
             radgr = ezgr + tgr*rsbgr*(1.0+knbgr*rsbgr)
             rsbcl = pixSURF%clref/(1.0-pixSURF%clref*sbcl )
             radcl = ezcl + rsbcl*tcl*(1.0+knbcl*rsbcl)
             xnvalc = -LOG10( radcl*pixSURF%clfrac  &
                            + radgr*(1.0-pixSURF%clfrac) )
          ENDIF
          IF( radgr <=  0.0 .AND. pixSURF%clfrac < 1.0 ) THEN
             WRITE(msg, *) "radgr =", radgr, ", refl =", pixSURF%ref, &
                           ", clfrac =", pixSURF%clfrac, &
                           ", pixID =", pixGEO%id
             call tell_log(0,msg)
             skipit=.TRUE.
             RETURN
          ENDIF

          IF( radcl <=  0.0 .AND. pixSURF%clfrac > 0.0 ) THEN
             WRITE(msg, *) "radcl =", radcl, ", refl =", pixSURF%ref, &
                           ", clfrac =", pixSURF%clfrac, &
                           ", pixID =", pixGEO%id
             call tell_log(0,msg)
             skipit=.TRUE.
             RETURN
          ENDIF

          IF( xnvalm(iwl_oz) <= xnvalc ) THEN
             ! computed value high
             ! check if at low profile bound
             IF( ANY( iprof .EQ. idx_lb ) ) THEN
                ! at low bound
                xlown = xnvalc
                iplow = iprof
                iprof = iplow + 1
                lb_found = .TRUE.
             ELSE
                ! not a low bound
                iphigh = iprof
                xhighn = xnvalc
                iprof = iphigh - 1
                ub_found = .TRUE.
             ENDIF
          ELSE
             ! computed value low
             ! check if at high bound
             IF( ANY( iprof .EQ. idx_ub ) ) THEN
                ! at high bound
                iphigh = iprof
                iprof = iphigh -1
                xhighn = xnvalc
                ub_found = .TRUE.
             ELSE
                ! not a high bound
                xlown = xnvalc
                iplow = iprof
                iprof = iplow + 1
                lb_found = .TRUE.
             ENDIF
          ENDIF
          IF( lb_found .AND. ub_found ) EXIT
        ENDDO

        ! table lookups complete, compute ozone

        omeghi  = prof_toz(iphigh) - terroz(iphigh)*pixGEO%fteran
        omeglo  = prof_toz(iplow) - terroz(iplow)*pixGEO%fteran

        dndoz_o = (xhighn-xlown) / (omeghi-omeglo)
        estozn  = omeglo + (xnvalm(iwl_oz)-xlown)/dndoz_o
        dndoz_o = 100.0*dndoz_o
        skipit = .FALSE.

      END FUNCTION O3T_oznot

!!Description:
! FUNCTION O3T_calcRefl
!   This function use the interpolation indices for table lookup, and
!   perform table interpolations in pressure and ozone, assumed ground
!   and cloud reflectivities,  then calculate cloud fraction.
!   check for the prsence of snow or sunglint, for which clear sky is
!   assumed. calculate reflectivity by inverting the surface radiance
!   formula for clear or completly cloudy, and use the partial cloud model
!   for partly clouded scenes.
!!Input Parameters:
!   iplow    : the index of lower bound ozone profile
!   ozonin   : the input ozone amount in D.U. (normally between iplow, iplow+1)
!   iwl_refl : the wavelength index used for reflectance calcualtion
!   xnvalm   : the measured nvalues for the all the wavelength
!   coefs   : Lagrangian interpolation coefficients associated with
!             pixGEO
!   pixGEO  : the geometry inforamtion for the pixel to the output
!             parameters are associated with. pixGEO has angular and
!             terrain and cloud pressure info.
!!Output Parameters:
!   pixSURF : the surface, cloud reflectance, and cloud fraction
!!Return
!   satus : The SMF stauts
!!References and Credits
! Written by
! Kai Yang
! University of Maryland Baltimore County
!!END Description:

      FUNCTION O3T_calcRefl( iplow, ozonin, iwl_refl, xnvalm, coefs, &
                             pixGEO, pixSURF, dndoz_r, first_call ) &
                           RESULT(errstat)
        INTEGER (KIND=4), INTENT(IN) :: iplow, iwl_refl
        REAL (KIND=4), INTENT(IN) :: ozonin
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: xnvalm
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        TYPE (O3T_pixgeo_type), INTENT(INOUT) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
        LOGICAL (KIND=4), INTENT(INOUT) :: first_call
        REAL (KIND=4), INTENT(OUT) :: dndoz_r
        REAL (KIND=4) :: grref,  clref,  clfrac
        REAL (KIND=4) :: pgrref, pclref, pclfrac
        REAL (KIND=4), DIMENSION(2), SAVE :: ezgr, tgr, sbgr, knbgr, &
                                             ezcl, tcl, sbcl, knbcl
        REAL (KIND=4) :: alb,  den,  ref
        REAL (KIND=4) :: palb, pden, pref, pxnvalm
        REAL (KIND=4) :: grad, crad
        REAL (KIND=4) :: ozfrac, psi
        REAL (KIND=4) :: omeglo, omeghi
        REAL (KIND=4) :: ezgcor, trgcor, sbgro, knbgro, &
                         ezccor, trccor, sbclo, knbclo
        INTEGER (KIND=4), SAVE :: ip_l = -1, pixel_id = -32767, &
                                  iwl_s = -1 !, ip_u = -1
        INTEGER (KIND=4) :: ioz, iprof
        INTEGER (KIND=4) :: errstat
        CHARACTER (LEN =255) :: msg

        errstat = 0
        alb  = 10.0**(-xnvalm(iwl_refl))

        ozfrac  = O3T_ozfraction( ozonin, iplow, pixGEO%fteran, omeglo, omeghi )
        IF( ozfrac <= -4.0 ) THEN
           WRITE( msg, '(A,I8,A,I3,2(A,F10.4))') 'pixID='  , pixel_id,     &
                                            ',iplow=' , iplow,        &
                                            ',ozonin=', ozonin,       &
                                            ',fteran=', pixGEO%fteran
           call tell_error(tell_runtime_error, "O3T_calcRefl: ozfrac <= -4", &
                errstat)
           RETURN
        ENDIF

        iprof = iplow
        IF( iplow .NE. ip_l .OR. pixGEO%id .NE. pixel_id &
                            .OR. iwl_refl  .NE. iwl_s ) THEN
           DO ioz = 1, 2
             errstat = O3T_iztrsb( nvRRS, iwl_refl, iprof, pixGEO, coefs, &
                                  ezgr(ioz), tgr(ioz), sbgr(ioz), knbgr(ioz), &
                                  ezcl(ioz), tcl(ioz), sbcl(ioz), knbcl(ioz)  )
             IF( errstat .NE. 0 ) THEN
               call tell_error(tell_runtime_error, &
                    "O3T_calcRefl: O3T_iztrsb error" , errstat)
                RETURN
             ENDIF
             iprof = iprof + 1
           ENDDO
           ip_l = iplow
           pixel_id = pixGEO%id
           iwl_s = iwl_refl
        ENDIF

        clref = 0.80
        IF( first_call ) THEN
           grref = 0.15
        ELSE
           grref = pixSURF%grref;
           IF( grref < 0.15 ) grref = 0.15
        ENDIF

        grad = O3T_rad( ezgr, tgr, sbgr, knbgr, &
                        grref, ozfrac, ezgcor, &
                        trgcor, sbgro, knbgro )

        crad = O3T_rad( ezcl, tcl, sbcl, knbcl, &
                        clref, ozfrac, ezccor, &
                        trccor, sbclo, knbclo )

        ! perturb measured reflectivity channel n value by 0.2%
        ! all varaibles beginning with p are from perturbed n value
        ! (needed for dN/dR calculation)
        pxnvalm = 1.002*xnvalm(iwl_refl)
        palb = 10.0**(-pxnvalm)
        pgrref = grref
        pclref = clref

        IF( ABS( grad - crad ) > EPSILON10 ) THEN
           clfrac  = ( alb-grad) / (crad-grad)
           pclfrac = (palb-grad) / (crad-grad)
        ELSE
           clfrac  = 0.0
           pclfrac = 0.0
        ENDIF

        !  assume clear sky if snow is present (v7 used ge 5 total.f)
        IF( pixSURF%isnow == 10 ) clfrac = 0

        IF( clfrac > 0.0 .AND. pixSURF%surface_category == 0 &  !! cloudy ocean
           .AND. first_call ) THEN
           psi = real( ACOS(pixGEO%cos_sza*pixGEO%cos_vza &
                    + pixGEO%sin_sza*pixGEO%sin_vza*pixGEO%cphi ) &
                    *RADtoDEG , kind=4)

           IF( psi < 20.0 ) THEN
              den =  alb  - ezgcor                 !! compute Refl331 assuming
              IF( ABS(den) > EPSILON10 ) THEN       ! clear scene so that
                 ref = 1. / (trgcor/den + sbgro )   ! O3T_glnchk can compute
              ELSE                                  ! corresponding residue at
                 ref = 0.                           ! 360 nm., leaving out RRS
              ENDIF                                !! correction for this calc.
              errstat = O3T_glnchk( iplow, xnvalm, ref, &
                                   coefs, pixGEO, pixSURF )
              IF( errstat .NE. 0 ) THEN
                call tell_error (tell_runtime_error, &
                     "O3T_calcRefl: O3T_glnchk error", errstat)
                 RETURN
              ENDIF
           ENDIF
        ENDIF

        !!When glint cloud fraction less than the calculated cloud fraction,
        !!assume cloud free.
        !!IF( pixSURF%glint_flag .AND. pixSURF%glnfrc < clfrac ) clfrac = 0.0

        IF( clfrac <= 0.0 ) THEN
        ! calculate reflectivity for cloud free
           den  =  alb - ezgcor
           pden = palb - ezgcor
           IF( ABS(den) > EPSILON10 ) THEN
              ref = 1. / (trgcor/den + sbgro)
           ELSE
              ref = 0.0
           ENDIF
           IF( ABS(pden) > EPSILON10 ) THEN
              pref = 1. / (trgcor/pden + sbgro)
           ELSE
              pref = 0.0
           ENDIF
           clfrac  = 0.0
           pclfrac = 0.0
           pixSURF%grref = ref
           pixSURF%clref = clref
           grref         = ref
        ELSE IF( clfrac >= 1.0 ) THEN
        !  calculate reflectivity for cloudy case
           den  =  alb - ezccor
           pden = palb - ezccor
           IF( ABS(den) > EPSILON10 ) THEN
              ref = 1. / (trccor/den + sbclo)
           ELSE
              ref = 0.0
           ENDIF
           IF( ABS(pden) > EPSILON10 ) THEN
              pref = 1. / (trccor/pden + sbclo)
           ELSE
              pref = 0.0
           ENDIF
           clfrac  = 1.0
           pclfrac = 1.0
           pixSURF%grref = grref
           pixSURF%clref = ref
                   clref = ref
        ELSE
        ! calculate reflectivity using cloud fraction and nominal refls
           ref =   grref +  clfrac*( clref- grref)
           pref = pgrref + pclfrac*(pclref-pgrref)
           pixSURF%grref = grref
           pixSURF%clref = clref
        ENDIF
        pixSURF%clfrac  = clfrac
        pixSURF%pclfrac = pclfrac
        pixSURF%ref     = ref
        pixSURF%pref    = pref
        IF( first_call ) dndoz_r = O3T_dndomega( ezgr, tgr, sbgr, knbgr, &
                                                 ezcl, tcl, sbcl, knbcl, &
                                                 grref, clref, clfrac, &
                                                 omeghi, omeglo )
        first_call = .FALSE.

      END FUNCTION O3T_calcRefl




      FUNCTION O3T_residue( iplow, estozn, xnvalm, &
           coefs, pixGEO, pixSURF, &
           n_residue, dndomega_t, dndr, nlambda ) &
           RESULT( errstat )

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iplow
        INTEGER (KIND=4), OPTIONAL, INTENT(IN) :: nlambda
        REAL (KIND=4), INTENT(IN) :: estozn
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: xnvalm
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        TYPE (O3T_pixgeo_type), INTENT(INOUT) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
        REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: n_residue, dndomega_t, dndr
        REAL (KIND=4), DIMENSION(2) :: ezgr, tgr, sbgr, knbgr, &
             ezcl, tcl, sbcl, knbcl
        REAL (KIND=4), DIMENSION(2) :: radgr, radcl, xnval, pxnval, &
             rsbgr, rsbcl
        REAL (KIND=4), DIMENSION(SIZE(xnvalm)) :: xnvalc, pxnvalc
        INTEGER (KIND=4) :: ierr, errstat
        INTEGER (KIND=4) :: iwl, nwl, ioz, iprof
        REAL (KIND=4) :: omeglo, omeghi, ozfrac, clfrac, pclfrac
        REAL (KIND=4) :: grref, clref, ref, pref
        REAL (KIND=4) :: Ic331, Im331, Rm360
        REAL (KIND=4), SAVE :: Rc360 = 0.8

        errstat = 0

        ozfrac  = O3T_ozfraction( estozn, iplow, pixGEO%fteran, omeglo, omeghi )
        clfrac = pixSURF%clfrac
        grref  = pixSURF%grref
        clref  = pixSURF%clref
        pclfrac= pixSURF%pclfrac
        ref    = pixSURF%ref
        pref   = pixSURF%pref

        IF( SIZE( xnvalm ) .NE. SIZE( n_residue ) .OR. &
             SIZE( xnvalm ) .NE. SIZE( dndomega_t ) .OR. &
             SIZE( xnvalm ) .NE. SIZE( dndr ) ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_residue: array size not equal", errstat)
          RETURN
        ENDIF

        IF( PRESENT( nlambda ) ) THEN
          nwl = nlambda
          IF( nwl > SIZE(xnvalm) ) THEN
          call tell_error(tell_runtime_error, &
               "O3T_residue: nlambda too large", errstat)
            RETURN
          ENDIF
        ELSE
          nwl = SIZE( xnvalm )
          IF( nwl > nwl_com ) nwl = nwl_com
        ENDIF

        DO iwl = 1, nwl
          iprof = iplow
          DO ioz = 1, 2
            errstat = O3T_iztrsb( nvRRS, iwl, iprof, pixGEO, coefs, &
                 ezgr(ioz), tgr(ioz), sbgr(ioz), knbgr(ioz), &
                 ezcl(ioz), tcl(ioz), sbcl(ioz), knbcl(ioz) )
            iprof = iprof + 1
            IF( errstat .NE. 0 ) THEN
              call tell_error(tell_runtime_error, &
                   "O3T_residue: O3T_iztrsb error", errstat)
              RETURN
            ENDIF
          ENDDO

          IF( clfrac <= 0.0 ) THEN
            rsbgr(:) = ref/( 1.0-ref*sbgr(:) )
            radgr(:) = ezgr(:)+tgr(:)*rsbgr(:)*(1.0+knbgr(:)*rsbgr(:) )
            xnval(:) = -LOG10(radgr(:))

            rsbgr(:) = pref/( 1.0-pref*sbgr(:) )
            radgr(:) = ezgr(:)+tgr(:)*rsbgr(:)*(1.0+knbgr(:)*rsbgr(:) )
            pxnval(:)= -LOG10(radgr(:))
            pixSURF%rcf1 = 0.0
            pixSURF%rcf2 = 0.0
          ELSE IF( clfrac >= 1.0 ) THEN
            rsbcl(:) = ref/( 1.0-ref*sbcl(:) )
            radcl(:) = ezcl(:)+tcl(:)*rsbcl(:)*(1.0+knbcl(:)*rsbcl(:) )
            xnval(:) = -LOG10(radcl(:))

            rsbcl(:) = pref/( 1.0-pref*sbcl(:) )
            radcl(:) = ezcl(:)+tcl(:)*rsbcl(:)*(1.0+knbcl(:)*rsbcl(:) )
            pxnval(:)= -LOG10(radcl(:))
            pixSURF%rcf1 = 1.0
            pixSURF%rcf2 = 1.0
          ELSE
            rsbgr(:) = grref/( 1.0-grref*sbgr(:) )
            rsbcl(:) = clref/( 1.0-clref*sbcl(:) )

            radgr(:) = ezgr(:)+tgr(:)*rsbgr(:)*(1.0+knbgr(:)*rsbgr(:) )
            radcl(:) = ezcl(:)+tcl(:)*rsbcl(:)*(1.0+knbcl(:)*rsbcl(:) )

            IF( iwl == pixSURF%iwl_refl_l ) THEN
              xnval(:) = -LOG10(radcl(:))
              Ic331    =  xnval(1) + (  xnval(2)- xnval(1) )*ozfrac
              Ic331    = 10.0**(-Ic331)
              Im331    = 10.0**(-xnvalm(iwl))
              !! Compute Radiative Cloud Fraction based on Ic331 & Im331
              !NB - use of ierr indicates output status unimportant 
              ierr     = O3T_rcf1( Ic331, Im331, pixSURF )
            ELSE IF( iwl == pixSURF%iwl_refl_h ) THEN
              Rm360    = real(PI*10.0**(-xnvalm(iwl))/pixGEO%cos_sza , kind=4)
              ierr     = O3T_rcf2( Rc360, Rm360, pixSURF )
            ENDIF

            xnval(:) = -LOG10(  clfrac*radcl(:) + (1.0- clfrac)*radgr(:) )
            pxnval(:)= -LOG10( pclfrac*radcl(:) + (1.0-pclfrac)*radgr(:) )
          ENDIF
          xnvalc(iwl)    =  xnval(1) + (  xnval(2)- xnval(1) )*ozfrac
          pxnvalc(iwl)   = pxnval(1) + ( pxnval(2)-pxnval(1) )*ozfrac
          dndomega_t(iwl)= 100.0*(xnval(2) - xnval(1))/(omeghi-omeglo)

        ENDDO
        n_residue(1:nwl) = 100.0*(xnvalm(1:nwl)  - xnvalc(1:nwl))
        dndr(1:nwl)      = 100.0*(pxnvalc(1:nwl)-xnvalc(1:nwl))/(pref-ref)

      END FUNCTION O3T_residue



      FUNCTION O3T_stp2oz( iwl_oz, iwl_refl, stp1oz, fgprf,&
                           dndx, delnT, dndomega_t, dndr, dxdomega, &
                           delshp, absrfl, dr, stp2prf, eff, stp2oz ) &
                           RESULT( errstat )
!      FUNCTION O3T_stp2oz( iwl_oz, iwl_refl, stp1oz, fgprf, fgtmp, &
!                           dndx, delnT, dndomega_t, dndr, dxdomega, &
!                           aprftm, delshp, absrfl, dr, stp2prf, eff, stp2oz ) &
!                           RESULT( errstat )

        implicit none

        INTEGER(KIND=4), INTENT(IN) :: iwl_oz, iwl_refl
        REAL(KIND=4), DIMENSION(:,:), INTENT(IN) :: dndx, delnT
        REAL(KIND=4), DIMENSION(:), INTENT(IN) :: dndomega_t, dndr
        REAL(KIND=4), INTENT(IN) :: stp1oz
        LOGICAL, INTENT(IN) :: absrfl
        REAL(KIND=4), INTENT(OUT) :: stp2oz, dr
        REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: stp2prf, eff
        REAL(KIND=4), DIMENSION(:), INTENT(IN) :: fgprf, delshp, dxdomega!, &
                                                  !fgtmp, aprftm
        REAL(KIND=4) :: dn_at_wloz, dn_at_wlrefl, crrcto, domega
        REAL(KIND=4),DIMENSION(SIZE(eff)) :: crrctx
        INTEGER(KIND=4) :: errstat

        errstat = 0

        IF( absrfl ) THEN
           crrctx(:) = (dndr( iwl_oz )/dndr( iwl_refl ))*dndx(iwl_refl,:)
           crrcto    = (dndr( iwl_oz )/dndr( iwl_refl ))*dndomega_t(iwl_refl)
        ELSE
           crrctx(:) = 0.0
           crrcto    = 0.0
        ENDIF
        eff(:)    = (dndx(iwl_oz,:)-crrctx(:))/(dndomega_t(iwl_oz)-crrcto)

        WHERE( eff < 0.0 ) eff = 0.0

        IF( absrfl ) THEN
           dn_at_wlrefl = SUM(delshp(:)*dndx(iwl_refl,:) + delnT(iwl_refl,:))
        ELSE
           dn_at_wlrefl = 0.0
        ENDIF

        dn_at_wloz   = SUM(delshp(:)*dndx(iwl_oz,  :) + delnT(iwl_oz,  :))
        dr         = dn_at_wlrefl / dndr(iwl_refl)
        dn_at_wloz = dn_at_wloz - dr * dndr(iwl_oz)
        domega     = dn_at_wloz / dndomega_t(iwl_oz)
        stp2oz     = stp1oz - domega
        IF( stp2oz < 0.0 .OR. stp2oz > 900.0 ) RETURN

        stp2prf(:) = fgprf(:) + delshp(:) - domega * dxdomega(:)
        errstat =  O3T_prof_check( stp2prf, fgprf + delshp )
      END FUNCTION O3T_stp2oz




      FUNCTION O3T_resadj( iwl_refl, res_stp1, dndr, dndx, delnT, &
                           fgprf, stp2prf, dr, res_stp2 ) &
                           RESULT( errstat )

        implicit none

        REAL(KIND=4), DIMENSION(:), INTENT(IN) :: res_stp1, dndr
        REAL(KIND=4), DIMENSION(:,:), INTENT(IN) :: dndx, delnT
        REAL(KIND=4), DIMENSION(:), INTENT(IN) :: fgprf, stp2prf
        REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: res_stp2
        REAL(KIND=4), INTENT(OUT) :: dr
        REAL(KIND=4) :: dn
        INTEGER(KIND=4), INTENT(IN) :: iwl_refl
        INTEGER(KIND=4) :: iwl, nwl
        INTEGER(KIND=4) :: errstat

        errstat = 0

        nwl = SIZE(res_stp1)

        dn = SUM( (stp2prf - fgprf)*dndx(iwl_refl,:) + delnT(iwl_refl,:) )
        dr = (res_stp1(iwl_refl) - dn)/dndr(iwl_refl)

        DO iwl = 1, nwl_sub
          dn = SUM( (stp2prf - fgprf)*dndx(iwl,:) + delnT(iwl,:) )
          res_stp2(iwl) = res_stp1(iwl)-dn-dr*dndr(iwl)
        ENDDO

        DO iwl = nwl_sub+1, nwl
          res_stp2(iwl) = res_stp1(iwl) - dr*dndr(iwl)
        ENDDO

      END FUNCTION O3T_resadj

      FUNCTION O3T_blwcld( stp2prf, pixGEO, pixSURF ) RESULT( oz_cld )

        implicit none

        REAL(KIND=4), DIMENSION(:), INTENT(IN) :: stp2prf
        TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(IN) :: pixSURF
        REAL(KIND=4) :: oz_cld
        REAL(KIND=4) :: pc, clfrac, fteran
        REAL(KIND=4) :: fcloud

        pc = real(pixGEO%pc, kind=4)
        fteran = pixGEO%fteran
        clfrac = pixSURF%clfrac

       ! calculate ozone amount beneath cloud for fully clouded scene
       ! if pcloud higher than 0.5 atm, we are in umkehr layer 0

       IF( pc >= 0.5 ) THEN
         fcloud = LOG10(1.00/pc)/LOG10(1.00/0.500)
         oz_cld = fcloud*stp2prf(1)
       ELSE IF( pc >= 0.250 ) THEN
       ! if pcloud lower than 0.5 atm and higher than 0.25, we are in layer 1
         fcloud = LOG10(0.500/pc)/LOG10(0.500/0.250)
         oz_cld = fcloud*stp2prf(2) + stp2prf(1)
       ELSE
       ! if pcloud lower than 0.25 atm, we are in layer 2
         fcloud = LOG10(0.250/pc)/LOG10(0.250/0.125)
         oz_cld = fcloud*stp2prf(3) + stp2prf(2) + stp2prf(1)
       ENDIF

       ! include correction for ozone below terrain pressure
       oz_cld = oz_cld - fteran*stp2prf(1)

       ! calculate amount beneath cloud for partially clouded situation

       IF( clfrac <= 0.) THEN
         clfrac = 0.0
       ELSE IF( clfrac >= 1.) THEN
         clfrac = 1.0
       ENDIF
       oz_cld = clfrac*oz_cld
      END FUNCTION O3T_blwcld




      FUNCTION O3T_getshp( latitude, jday, estozn, fgprf, apprf, delshp ) &
                           RESULT( errstat )

        implicit none

        REAL( KIND=4 ), INTENT(IN) :: latitude, estozn
        INTEGER( KIND=4 ), INTENT(IN) :: jday
        REAL( KIND=4 ), DIMENSION(:), INTENT(IN) :: fgprf
        REAL( KIND=4 ), DIMENSION(:), INTENT(OUT) :: apprf, delshp
        REAL( KIND=4 ),DIMENSION(SIZE(fgprf)) :: tmpprf
        INTEGER( KIND=4 ) :: i,j,k,l
        REAL( KIND=4 ) :: foz, flt, fmn
        INTEGER( KIND=4 ) :: errstat!, ierr

        errstat = 0

        ! determine interpolation indices
        j = int(( estozn - 25.0 )/50.0 - 1)
        IF( j <    1 ) j = 1
        IF( j >=  10 ) j = 9
        IF( ABS( latitude ) <=  65 .AND. j <  3 ) j = 3
        IF( ABS( latitude ) <=  35 .AND. j >= 5 ) j = 4
        foz = ( estozn - toz(j) )/50.0

        k = int(( latitude + 85.0)/10.0 + 1)
        IF( k <   1 ) k =  1
        IF( k >= 18 ) k = 17
        flt = ( latitude - rlats(k) )/10.0

        l = int(jday/30.5 + 0.51)
        IF( l < 1 .OR. l > 12 ) l = 12
        If( jday >= 15 ) THEN
           fmn = (jday-((l-1)*30.5+15.0))/((l*30.5+15.0)-((l-1)*30.5+15.0))
        ELSE
           fmn = ((jday+366.0)-((l-1)*30.5+15.0))/ &
                               ((l*30.5+15.0)-((l-1)*30.5+15.0))
        ENDIF

        DO i = 1, NLYR

          apprf(i) = (1.0-foz)*(1.0-flt)*(1.0-fmn)*tzaprf(i,j,k,l) &
                   + foz*(1.0-flt)*(1.0-fmn)*tzaprf(i,j+1,k,l) &
                   + foz*flt*(1.0-fmn)*tzaprf(i,j+1,k+1,l) &
                   + (1.0-foz)*flt*(1.0-fmn)*tzaprf(i,j,k+1,l) &
                   + (1.0-foz)*(1.0-flt)*fmn*tzaprf(i,j,k,l+1) &
                   + foz*(1.0-flt)*fmn*tzaprf(i,j+1,k,l+1) &
                   + foz*flt*fmn*tzaprf(i,j+1,k+1,l+1) &
                   + (1.0-foz)*flt*fmn*tzaprf(i,j,k+1,l+1)
        ENDDO

        IF( estozn < toz(j) ) THEN
           IF( ANY(apprf < 0.0)  ) THEN
              DO i = 1, NLYR
                tmpprf(i) = (1.0-flt)*(1.0-fmn)*tzaprf(i,j,k  ,l  ) &
                          +      flt *(1.0-fmn)*tzaprf(i,j,k+1,l  ) &
                          + (1.0-flt)*     fmn *tzaprf(i,j,k  ,l+1) &
                          +      flt *     fmn *tzaprf(i,j,k+1,l+1)
              ENDDO
              errstat =  O3T_prof_check( apprf, tmpprf )
           ENDIF
        ENDIF

        delshp(:) = apprf(:)-fgprf(:)

      END FUNCTION O3T_getshp

!     FUNCTION O3T_getshp_omps(fgprf, apprf, delshp) RESULT( status )
!     ! COMMENT-JYLI:
!     ! code logic has been folded into O3T_step2_omps().  Currently, this
!     ! routine is not being used.
!
!       REAL( KIND=4 ), DIMENSION(:), INTENT(IN)  :: fgprf
!       REAL( KIND=4 ), DIMENSION(:), INTENT(IN)  :: apprf
!       REAL( KIND=4 ), DIMENSION(:), INTENT(OUT) :: delshp
!       INTEGER( KIND=4 ) :: status, ierr
!
!       IF( ANY(apprf < 0.0) ) THEN
!            status = OZT_E_FAILURE
!            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
!                                  "negative O3 value(s) in apprf", &
!                                  "O3T_getshp_omps", zero )
!       ELSE
!          status = OZT_S_SUCCESS
!       ENDIF
!
!       delshp(1:NLYR) = apprf(1:NLYR) - fgprf(1:NLYR)
!
!     END FUNCTION O3T_getshp_omps

      FUNCTION O3T_glnchk( iplow, xnvalm, refl, coefs, pixGEO, pixSURF ) &
                           RESULT(errstat)

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iplow
        REAL (KIND=4), INTENT(IN) :: refl
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: xnvalm
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        TYPE (O3T_pixgeo_type), INTENT(INOUT) :: pixGEO
        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
        REAL (KIND=4), SAVE :: ezgr, tgr, sbgr, knbgr, rsbgr, &
                               ezcl, tcl, sbcl, knbcl
        REAL (KIND=4) :: grrad, xncalc, nres_glint
        INTEGER( KIND=4 ) :: iwl_glint
        INTEGER( KIND=4 ) :: errstat

        errstat = 0

        iwl_glint = pixSURF%iwl_glint
        errstat = O3T_iztrsb( nvRRS, iwl_glint, iplow, pixGEO, coefs, &
                             ezgr, tgr, sbgr, knbgr, ezcl, tcl, sbcl, knbcl )
        IF( errstat .NE. 0 ) THEN
          call tell_error(tell_runtime_error, "O3T_glnchk: O3T_iztrsb error", &
               errstat)
           RETURN
        ENDIF

        ! assume clear scene for glint check
        rsbgr = refl/(1.0-refl*sbgr)
        grrad = ezgr+tgr*rsbgr*(1.0+knbgr*rsbgr)
        xncalc= -LOG10(grrad)

        ! compute residue, logical glint indicator, and glint cloud frac

        nres_glint = 100.0 * (xnvalm(iwl_glint) - xncalc)

        IF( nres_glint <  pixSURF%nres_glnt_ub) pixSURF%glint_flag = .TRUE.
        IF( pixSURF%glint_flag ) THEN
           pixSURF%glnfrc=(1.0 - nres_glint/(refl*(-12.0)))
           IF(pixSURF%glnfrc > 1.0) pixSURF%glnfrc = 1.0
           IF(pixSURF%glnfrc < 0.0) pixSURF%glnfrc = 0.0
           IF(nres_glint < -8.0) pixSURF%glnfrc = 0.0
        ENDIF

      END FUNCTION O3T_glnchk

!      SUBROUTINE initRefl( grref, clref )
!        REAL (KIND=4), INTENT(OUT) :: grref, clref
!        grref = 0.15
!        clref = 0.80
!      END SUBROUTINE initRefl

!!Description:
! PRIVATE FUNCTION O3T_rad
!   This function performs log linear interpolation or extrapolation
!   for the input parameters of two bracketing ozone profile and produce
!   the normalized radiance.
!!Input Parameters:
!   ez(2), tr(2), sb(2) : ez, tg, ans sb already calculated for
!                         two pofiles at a specific set of angles
!                         and surface (terrain or cloud) height
!                         specified by pressure.
!   refl   : the surface reflectance
!   ozfrac : the linear interpolation fraction
!!Output Parameters:
!   ezo, to, sbo : optional output, intermediate quantities after liner
!                  interpolation.
!!Return
!   satus : the normalized radiance
!!References and Credits
! Written by
! Kai Yang
! University of Maryland Baltimore County
!!END Description:

      FUNCTION O3T_rad( ez, t, sb, knb, refl, ozfrac,  &
                        ezo,to,sbo,knbo  ) RESULT( Lrad )
        REAL (KIND=4), DIMENSION(2), INTENT(IN) :: ez, t, sb, knb
        REAL (KIND=4), INTENT(IN) :: ozfrac, refl
        REAL (KIND=4), INTENT(OUT), OPTIONAL :: ezo, to, sbo, knbo
        REAL (KIND=4) :: ezl, tl, rsb
        REAL (KIND=4) :: ezlog, tlog, sblin, knblin
        REAL (KIND=4) :: Lrad

        ezlog = LOG10(ez(1)) + LOG10(ez(2) / ez(1))*ozfrac
        tlog  = LOG10( t(1)) + LOG10( t(2) /  t(1))*ozfrac
        sblin =       sb(1)  +      (sb(2) - sb(1))*ozfrac
        knblin=      knb(1)  +     (knb(2) -knb(1))*ozfrac

        ezl   = 10.0**ezlog
        tl    = 10.0**tlog
        rsb   = refl/( 1.0 - refl*sblin )
        Lrad  = ezl + rsb*tl*( 1.0 + knblin*rsb )

        IF( PRESENT(ezo ) ) ezo  = ezl
        IF( PRESENT( to ) )  to  = tl
        IF( PRESENT(sbo ) ) sbo  = sblin
        IF( PRESENT(knbo) ) knbo = knblin

      END FUNCTION O3T_rad

      FUNCTION O3T_dndomega( ezgr, tgr, sbgr, knbgr, &
                             ezcl, tcl, sbcl, knbcl, &
                             grref, clref, clfrac, omeghi, omeglo ) &
                             RESULT( dndomega )
        REAL (KIND=4), DIMENSION(2), INTENT(IN) :: ezgr, tgr, sbgr, knbgr, &
                                                   ezcl, tcl, sbcl, knbcl

        REAL (KIND=4), DIMENSION(2) :: radgr, radcl, rsbgr, rsbcl
        REAL (KIND=4), INTENT(IN) :: grref, clref, clfrac, omeghi, omeglo
        REAL (KIND=4) :: dndomega

        IF( clfrac <= 0.0 ) THEN
           rsbgr(:) = grref/(1.0-grref*sbgr(:))
           radgr(:) = ezgr(:)+rsbgr(:)*tgr(:)*(1.0+rsbgr(:)*knbgr(:))
           dndomega = (-100.0*LOG10(radgr(2)/radgr(1)))/(omeghi-omeglo)
        ELSE IF( clfrac >= 1.0 ) THEN
           rsbcl(:) = clref/(1.0-clref*sbcl(:))
           radcl(:) = ezcl(:)+rsbcl(:)*tcl(:)*(1.0+rsbcl(:)*knbcl(:))
           dndomega = (-100.0*LOG10(radcl(2)/radcl(1)))/(omeghi-omeglo)
        ELSE
           rsbgr(:) = grref/(1.0-grref*sbgr(:))
           radgr(:) = ezgr(:)+rsbgr(:)*tgr(:)*(1.0+rsbgr(:)*knbgr(:))
           rsbcl(:) = clref/(1.0-clref*sbcl(:))
           radcl(:) = ezcl(:)+rsbcl(:)*tcl(:)*(1.0+rsbcl(:)*knbcl(:))
           dndomega = ( -100.0*LOG10( &
                        ( radgr(2)+clfrac*(radcl(2)-radgr(2)) ) &
                       /( radgr(1)+clfrac*(radcl(1)-radgr(1)) ) &
                      )             ) /( omeghi-omeglo )
        ENDIF
      END FUNCTION O3T_dndomega

      FUNCTION O3T_descendQ( lat1, lat2 ) RESULT( descendQ )
        REAL (KIND=4), INTENT(IN) :: lat1, lat2
        LOGICAL :: descendQ
        IF( lat1 > lat2 ) THEN
           descendQ = .TRUE.        !! descending part of the orbit
        ELSE
           descendQ = .FALSE.       !! ascending part of the orbit
        ENDIF
      END FUNCTION O3T_descendQ


      !FIXME - will need to be updated for TEMPO
      FUNCTION O3T_getSatName( ShortName, nval_lun ) RESULT( satname )

        use tell_module

        implicit none

        CHARACTER(LEN=*), INTENT(IN) :: ShortName
        INTEGER (KIND = 4), INTENT(IN) :: nval_lun
        CHARACTER(LEN=2) :: satname
        CHARACTER (LEN =255) :: msg, fn_nval
        INTEGER (KIND = 4) :: version, status, di, errstat

        !! get the nval table filename and path
        errstat = 0
        version = 1
        status = PGS_PC_GetReference( nval_lun, version, fn_nval )
        IF( status /= PGS_S_SUCCESS ) THEN
          WRITE( msg,'(A,I9)' ) &
               "O3T_getSatName: get LUT name failed at LUN = ", nval_lun
          call tell_error(tell_io_read_error, msg, errstat)
          fn_nval = ""
        ENDIF

        !! remove path
        di = INDEX( fn_nval, "/", BACK=.TRUE. )
        fn_nval = fn_nval( di+1: )

        IF( INDEX( ShortName, "TOMS" )>0 .AND. INDEX( fn_nval, "TOMS" )>0 ) THEN
          IF( INDEX( fn_nval, "N7" ) > 0 ) THEN
            satname = "N7";
          ELSE IF( INDEX( fn_nval, "EP" ) > 0 ) THEN
            satname = "EP";
          ELSE IF( INDEX( fn_nval, "M3" ) > 0 ) THEN
            satname = "M3";
          ELSE IF( INDEX( fn_nval, "AD" ) > 0 ) THEN
            satname = "AD";
          ELSE
            WRITE( msg,'(A)' ) "O3T_getSatName: ShortName:"// &
                 TRIM( ShortName ) // &
                 "and NVAL file name"// TRIM(fn_nval) // &
                 "incompatiable."
            call tell_error(tell_runtime_error, msg, errstat)
            satname = "";
          ENDIF
        ELSE IF ( TRIM(ShortName) /= "OMTO3" ) THEN
          satname = "OM";
        ELSE IF ( TRIM(ShortName) /= "GMTO3" ) THEN
          satname = "OM";
        ELSE
          WRITE( msg,'(A)' ) "O3T_getSatName: ShortName:"// &
               TRIM( ShortName ) // &
               "and NVAL file name"// TRIM(fn_nval) // &
               "incompatiable."
          call tell_error(tell_runtime_error, msg, errstat)
          satname = "";
        ENDIF
        WRITE( msg,'(A)' ) "used satname:"// satname
        call tell_log(4, msg)
        RETURN
      END FUNCTION O3T_getSatName



      FUNCTION O3T_rcf2( Rc360, Rm360, pixSURF ) RESULT( errstat )

        implicit none

        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
        REAL (KIND=4), INTENT(IN) :: Rc360, Rm360
        INTEGER (KIND=4) :: ierr, errstat
        REAL (KIND=4) :: fc
        REAL (KIND=4), SAVE :: fill_float32
        LOGICAL (KIND=4), SAVE :: firsttime = .TRUE.
        INTEGER (KIND=4), EXTERNAL :: r4Fill

        errstat = 0

        IF( firsttime ) THEN
          ierr = r4Fill( fill_float32 )
          firsttime = .FALSE.
        ENDIF

        fc = pixSURF%clfrac
        IF( fc < 0.0 .OR. fc > 1.0 ) THEN
          call tell_error(tell_invalid_parm_error, &
               "O3T_rcf2: Effective cloud fraction out of range" , errstat)
          pixSURF%rcf2 = fill_float32
          RETURN
        ELSE IF( fc == 0.0 .OR. fc == 1.0 ) THEN
          pixSURF%rcf2 = fc
          RETURN
        ELSE
          IF( Rm360 <= 0.0 ) THEN
            call tell_error(tell_invalid_parm_error, &
                 "O3T_rcf2: Rm360 out of range" , errstat)
            pixSURF%rcf2 = fill_float32
            RETURN
          ELSE
            pixSURF%rcf2 = fc*Rc360/Rm360
          ENDIF
        ENDIF
        RETURN
      END FUNCTION O3T_rcf2

      FUNCTION O3T_rcf1( Ic331, Im331, pixSURF ) RESULT( errstat )

        implicit none

        TYPE (O3T_pixcover_type), INTENT(INOUT) :: pixSURF
        REAL (KIND=4), INTENT(IN) :: Ic331, Im331
        INTEGER (KIND=4) :: ierr, errstat
        REAL (KIND=4) :: fc
        REAL (KIND=4), SAVE :: fill_float32
        LOGICAL (KIND=4), SAVE :: firsttime = .TRUE.
        INTEGER (KIND=4), EXTERNAL :: r4Fill

        errstat = 0

        IF( firsttime ) THEN
          ierr = r4Fill( fill_float32 )
          firsttime = .FALSE.
        ENDIF

        fc = pixSURF%clfrac
        IF( fc < 0.0 .OR. fc > 1.0 ) THEN
          call tell_error(tell_invalid_parm_error, &
               "O3T_rcf1: Effective cloud fraction out of range", errstat)
          pixSURF%rcf2 = fill_float32
          RETURN
        ELSE IF( fc == 0.0 .OR. fc == 1.0 ) THEN
          pixSURF%rcf1 = fc
          RETURN
        ELSE
          IF( Im331 <= 0.0 ) THEN
            call tell_error(tell_invalid_parm_error, &
                 "O3T_rcf1: Im331 out of range", errstat)
            pixSURF%rcf1 = fill_float32
            RETURN
          ELSE
            pixSURF%rcf1 = fc*Ic331/Im331
          ENDIF
        ENDIF
        RETURN
      END FUNCTION O3T_rcf1

END MODULE O3T_class
