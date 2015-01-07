!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_stnprof_class
! this module contains the standard O3 profiles used in the master
! Look-Up Table calculations. 
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
! Adopted blockdata.f from TOMS V8 code.
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
MODULE O3T_stnprof_class
    USE OMI_SMF_class
    IMPLICIT NONE
    INTEGER (KIND = 4), PARAMETER, PRIVATE :: zero = 0

    INTEGER (KIND=4), PARAMETER :: NPROF = 21, NLYR = 11, NBlat = 3
    REAL (KIND=4), PARAMETER :: oz_step = 50.0
    REAL (KIND=4), DIMENSION(NPROF) :: prof_toz = (/225.,275.,325.,225.,275.,&
                                                    325.,375.,425.,475.,525.,&
                                                    575.,125.,175.,225.,275.,&
                                                    325.,375.,425.,475.,525.,&
                                                    575./)

    REAL (KIND=4), DIMENSION(NBlat) :: oz_lb = (/ 225.,225.,125. /)
    REAL (KIND=4), DIMENSION(NBlat) :: oz_ub = (/ 325.,575.,575. /)
    INTEGER (KIND=4), DIMENSION(NBlat) :: idx_lb = (/ 1,  4, 12 /)
    INTEGER (KIND=4), DIMENSION(NBlat) :: idx_ub = (/ 3, 11, 21 /)

    REAL (KIND=4), DIMENSION(NPROF) :: terroz =  &
          (/15.0,15.0,15.0,15.09,15.09,15.09,15.09,15.09,15.09,15.09,15.09, &
            14.11,14.11,14.11,14.11,14.11,14.11,14.11,14.11,14.11,14.11 /)

    REAL (KIND=4), DIMENSION(NLYR) :: fgtmp = (/ 269.7,236.4,216.7,216.7, &
                                                 217.9,222.4,226.7,237.0, &
                                                 250.8,265.9,270.6 /)
    REAL (KIND=4), DIMENSION(NLYR,NPROF) :: stdprf  =  RESHAPE( &
   (/ 15.00, 9.00,  4.94,  5.91, 25.75,62.60,57.00,29.40,10.90, 3.20, 1.30,&
      15.00, 9.00,  5.61, 14.66, 50.68,78.25,57.00,29.40,10.90, 3.20, 1.30,&
      15.00, 9.00, 10.98, 30.51, 71.05,86.66,57.00,29.40,10.90, 3.20, 1.30,&
      15.09, 9.03,  8.68, 18.55, 45.11,53.84,36.26,22.56,10.75, 3.70, 1.43,&
      15.09,11.46, 15.41, 29.67, 61.57,63.64,39.72,22.56,10.75, 3.70, 1.43,&
      15.09,14.08, 25.51, 43.64, 75.49,70.32,42.43,22.56,10.75, 3.70, 1.43,&
      15.09,16.91, 39.73, 61.73, 87.09,73.51,42.50,22.56,10.75, 3.70, 1.43,&
      15.09,19.93, 56.13, 82.00, 96.26,74.25,42.90,22.56,10.75, 3.70, 1.43,&
      15.09,23.15, 73.24,103.71,104.09,74.28,43.00,22.56,10.75, 3.70, 1.43,&
      15.09,26.57, 91.54,126.67,110.15,73.54,43.00,22.56,10.75, 3.70, 1.43,&
      15.09,30.19,111.11,150.11,114.67,72.39,43.00,22.56,10.75, 3.70, 1.43,&
      14.11, 8.47, 18.56,  5.91,  5.00,23.28,22.94,15.87, 7.37, 2.53, 0.96,&
      14.11, 8.73, 18.53, 16.91, 23.61,33.26,30.08,17.91, 8.06, 2.73, 1.07,&
      14.11, 9.89, 21.62, 28.94, 41.27,41.10,35.61,19.65, 8.71, 2.93, 1.17,&
      14.11,11.85, 27.72, 42.14, 57.09,47.89,39.54,21.02, 9.28, 3.11, 1.25,&
      14.11,14.61, 36.84, 56.58, 71.10,53.60,41.88,21.98, 9.72, 3.27, 1.31,&
      14.11,17.80, 49.19, 72.51, 83.51,58.34,42.33,22.35,10.09, 3.41, 1.36,&
      14.11,22.50, 64.50, 89.46, 93.69,61.65,41.70,22.40,10.08, 3.51, 1.40,&
      14.11,29.41, 82.32,106.83,101.07,63.21,40.80,22.34, 9.95, 3.54, 1.42,&
      14.11,38.08,102.60,124.66,105.87,63.22,39.45,22.17, 9.90, 3.52, 1.42,&
      14.11,48.45,125.21,142.90,108.23,61.79,37.65,21.90, 9.88, 3.47, 1.41/), &
      (/ NLYR, NPROF /) )

    PUBLIC  :: O3T_stnprof_idxf
    PUBLIC  :: O3T_ozfraction
!    PUBLIC  :: O3T_stnprof
    PUBLIC  :: O3T_stnprof_1stG
    PUBLIC  :: O3T_prof_check

    CONTAINS

!!Description:
! FUNCTION FUNCTION  O3T_stnprof_idxf
!   This function calculates the profile index based on the input
!   latitude band and the ozone amount. This the ozone amount is in 
!   between the return profile index iplow and iplow+1 for most
!   realistic ozone amount in the latitude band.
!!Input Parameters:
!   oz     : the ozone amount (in D.U.)
!   ilat   : the latitude band (1, 2, or 3)
!!Output Parameters:
!   None
!!Return
!   iplow : the lower bound profile index (the input ozone amount is
!           between column ozone of iplow and iplow+1) or -1 if 
!           input ilat is out of range.
!!References and Credits
! Written by 
! Kai Yang 
! University of Maryland Baltimore County
!!END Description:

     FUNCTION  O3T_stnprof_idxf( oz, ilat ) RESULT( iplow )
       REAL (KIND=4), INTENT(IN) :: oz
       INTEGER (KIND=4), INTENT(IN) :: ilat
       INTEGER (KIND=4) :: iplow, ierr

       IF( ilat < 1 .OR. ilat > 3 ) THEN
          iplow = -1
          ierr = OMI_SMF_setmsg( OZT_E_INPUT, "ilat out of range", &
                                 "O3T_stnprof_idxf", zero )
          RETURN
       ENDIF

       iplow = FLOOR((oz-oz_lb(ilat))/oz_step) + idx_lb(ilat)
       IF (iplow >= idx_ub(ilat)) THEN
         iplow = idx_ub(ilat) - 1
       ELSE IF (iplow < idx_lb(ilat)) THEN
         iplow = idx_lb(ilat)
       ENDIF
     END FUNCTION  O3T_stnprof_idxf

!!Description:
! FUNCTION O3T_ozfraction
!   This function calculates the linear interpolation (if input oz
!   lies between the input lower loading ozone profile iplow) or 
!   linear extrapolation (if input oz lie outside the two profiles
!   iplow and iplow+1) fraction. Note that the ozone amounts are 
!   pressure corrected for the profiles.
!!Input Parameters:
!   oz     : the ozone amount
!   iplow: the lower limit ozone profile index
!   fteran : the pressure correction factor (the higher the terrain result
!            in lesser the amount of column ozone for the profile
!!Output Parameters:
!   None
!!Return
!   ozfrac : the linear interpolation (or linear extrapolation) fraction
!            when input iplow is within the range of profile index, other
!            -HUGE(1.0) is returned.
!!References and Credits
! Written by 
! Kai Yang 
! University of Maryland Baltimore County
!!END Description:

     FUNCTION O3T_ozfraction( oz, iplow, fteran, omeglo_k, omeghi_k ) &
                              RESULT( ozfrac )
       REAL (KIND=4), INTENT(IN) :: oz 
       REAL (KIND=4), INTENT(IN), OPTIONAL :: fteran
       INTEGER (KIND=4), INTENT(IN) :: iplow
       INTEGER (KIND=4) :: iphigh, ierr
       REAL (KIND=4), INTENT(OUT), OPTIONAL ::omeglo_k, omeghi_k 
       REAL (KIND=4) :: ozfrac, omeglo, omeghi
       IF( iplow < 1 .OR. ANY(iplow .EQ. idx_ub ) .OR. iplow > 21 ) THEN
          ozfrac = -HUGE(1.0)
          ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iplow out of range", &
                                 "O3T_ozfraction", zero )
          RETURN
       ENDIF
       iphigh = iplow + 1 
       IF( PRESENT( fteran ) )THEN
          omeglo = prof_toz(iplow) - terroz(iplow)*fteran
          omeghi = prof_toz(iphigh) - terroz(iphigh)*fteran
       ELSE
          omeglo = prof_toz(iplow) 
          omeghi = prof_toz(iphigh)
       ENDIF
       ozfrac = (oz-omeglo)/(omeghi-omeglo)
       IF( PRESENT( omeglo_k ) ) omeglo_k = omeglo
       IF( PRESENT( omeghi_k ) ) omeghi_k = omeghi
     END FUNCTION O3T_ozfraction

!!Description:
! FUNCTION O3T_stnprof
!   This function calculates standard ozone profile for a specific total
!   ozone at a specific latitude using linear interpolation of the TOMS V8
!   low, mid, and high standard profiles
!!Input Parameters:
!   latitude: latitude in degree
!   ozone  : the total ozone amount in D.U.
!!Output Parameters:
!   stnprof : the inepolated standard ozone profile
!!Return
!   status  : SMF status
!!References and Credits
! Written by
! Kai Yang
! University of Maryland Baltimore County
!!END Description:

!     FUNCTION O3T_stnprof( latitude, ozone, stnprof ) RESULT( status )
!       REAL (KIND=4), INTENT(IN) :: latitude, ozone
!       REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: stnprof
!       REAL (KIND=4), DIMENSION(SIZE(stnprof)) :: stnprm
!       INTEGER (KIND=4) :: nlayer
!       INTEGER (KIND=4) :: status, ierr
!       REAL (KIND=4) :: abslat, latfrac, ozfrac
!       INTEGER (KIND=4) :: ilat
!       INTEGER (KIND=4) :: indm1, indm2
!       INTEGER (KIND=4) :: indl1, indl2
!       INTEGER (KIND=4) :: indh1, indh2
!       CHARACTER (LEN =255) :: msg
!
!       status = OZT_S_SUCCESS
!       nlayer = SIZE( stnprof )
!
!       IF( nlayer < NLYR ) THEN
!          status = OZT_E_FAILURE
!          ierr = OMI_SMF_setmsg( OZT_E_INPUT, "array stnprof size too small", &
!                                 "O3T_stnprof", zero )
!          RETURN
!       ENDIF
!    
!       abslat = ABS( latitude )
!       ! Estimate standard profile for given ozone and latitude
!       ! by interpolating a mid latitude standard profile
!       IF( abslat >= 20.0 .AND. abslat <= 60.0 ) THEN
!          ilat = 2
!          indm1 = O3T_stnprof_idxf( ozone, ilat ) 
!          indm2 = indm1 + 1
!          ozfrac = O3T_ozfraction( ozone, indm1 )
!          stnprm = stdprf(:,indm1)+ozfrac*(stdprf(:,indm2)-stdprf(:,indm1))
!       ENDIF
!
!       ! low to mid latitude interpolation
!       IF( abslat < 45.0 ) THEN
!          ilat = 1
!          indl1 = O3T_stnprof_idxf( ozone, ilat ) 
!          indl2 = indl1 + 1
!          ozfrac = O3T_ozfraction( ozone, indl1 )
!          stnprof = stdprf(:,indl1)+ozfrac*(stdprf(:,indl2)-stdprf(:,indl1))
!
!          IF( abslat > 20.0 ) THEN
!             latfrac = (abslat - 20.0)/25.0    ! 25.0 = (45.0-20.0)
!             stnprof = stnprof + latfrac*(stnprm - stnprof)
!          ENDIF
!       ELSE
!          ilat = 3
!          indh1 = O3T_stnprof_idxf( ozone, ilat ) 
!          indh2 = indh1 + 1
!          ozfrac = O3T_ozfraction( ozone, indh1 )
!          stnprof = stdprf(:,indh1)+ozfrac*(stdprf(:,indh2)-stdprf(:,indh1))
!!         write(*,*) 'indh1, indh2 =', indh1, indh2
!!         write(*,*) 'stnprof =', stnprof
!
!          IF( abslat < 60.0 ) THEN
!             latfrac = (abslat - 45.0)/15.0    ! 15.0 = (60.0-45.0)
!             stnprof = stnprm + latfrac*(stnprof - stnprm)
!          ENDIF
!       ENDIF
!!      write(*,*) 'stnprof =', stnprof
!
!       IF( ABS( SUM(stnprof(1:NLYR)) - ozone ) > 0.1 ) THEN
!          WRITE( msg, * ) 'non-linearity error in stnprf:', &
!                          'SUM(stnprof) =',  SUM(stnprof(1:NLYR)), &
!                          'Input Ozone  =',  ozone, &
!                          'difference   =',   SUM(stnprof(1:NLYR)) - ozone   
!          ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, "O3T_stnprof", zero )
!       ENDIF
!
!     END FUNCTION O3T_stnprof

     FUNCTION O3T_stnprof_1stG( estozn, iplow, fgprf, dxdomega ) &
                                RESULT( status )
       REAL (KIND=4), INTENT(IN) :: estozn
       REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: dxdomega, fgprf
       INTEGER (KIND=4), INTENT(IN) :: iplow
       REAL (KIND=4) :: omeghi, omeglo, ozfrac, domega
       INTEGER (KIND=4) :: iphigh, status, ierr

       IF( SIZE( dxdomega ) < NLYR .OR. SIZE( fgprf ) < NLYR ) THEN
          status = OZT_E_FAILURE
          ierr = OMI_SMF_setmsg( OZT_E_INPUT, "input array size too small", &
                                 "O3T_stnprof_1stG", zero )
          RETURN
       ENDIF

       iphigh = iplow+1
       omeghi = prof_toz(iphigh)  
       omeglo = prof_toz(iplow)
       domega = omeghi - omeglo
       ozfrac = (estozn-omeglo)/domega
       dxdomega = (stdprf(:,iphigh)- stdprf(:,iplow))/domega
       fgprf    = stdprf(:,iplow) + (stdprf(:,iphigh)-stdprf(:,iplow))*ozfrac 

       IF( estozn < omeglo ) THEN
         status = O3T_prof_check( fgprf, stdprf(:,iplow) )
       ELSE
         status = OZT_S_SUCCESS
       ENDIF
       
     END FUNCTION O3T_stnprof_1stG

     FUNCTION O3T_prof_check( o3prof, o3prof_base ) RESULT( status )
       REAL (KIND=4), DIMENSION(:), INTENT(INOUT) :: o3prof
       REAL (KIND=4), DIMENSION(:), INTENT(IN)    :: o3prof_base
       REAL (KIND=4), DIMENSION(SIZE(o3prof_base)):: ratio
       INTEGER (KIND=4), DIMENSION(1) :: ml
       INTEGER (KIND=4) :: i, j, l, nL 
       INTEGER (KIND=4) :: status, ierr
       CHARACTER (LEN =255) :: msg

       !! make sure array sizes are the same.
       nL = SIZE( o3prof_base ) 
       IF( SIZE( o3prof ) /= nL ) THEN
          status = OZT_E_FAILURE
          ierr = OMI_SMF_setmsg( OZT_E_INPUT, "input array size not equal", &
                                 "O3T_prof_check", zero )
          RETURN
       ENDIF
       
       !! check to see if there is any value profile value is negative
       !! If yes, find the array index range where the negative values
       !! occur.
       i = 0
       DO l = 1, nL
         IF( o3prof(l) < 0.0 ) THEN
            IF( i == 0 ) THEN
               i = l
               j = l
            ELSE
               j = l
            ENDIF
         ENDIF
       ENDDO

       IF( i > 0 ) THEN          !! TRUE, further checking and correction 
          IF( SUM( o3prof ) < 0.0 ) THEN
             status = OZT_E_FAILURE
             WRITE( msg, * ) "input total O3 =", SUM( o3prof ), o3prof 
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_prof_check", zero )
             RETURN
          ENDIF

          !! make sure reference profile is all positive before continuing.
          IF( ANY( o3prof_base < 0.0 ) ) THEN
             status = OZT_E_FAILURE
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "negative O3 value(s) "//&
                                    "in input reference profile", &
                                    "O3T_prof_check", zero )
             RETURN
          ENDIF

          !! Extend the low and upper bound by 1 if not yet at the start or
          !! the end of the array
          IF( i > 1  ) i = i - 1
          IF( j < nL ) j = j + 1

          !! Change the boundary until within the range, the sum of the
          !! the ozone is positive.
          DO WHILE( SUM( o3prof(i:j) ) <= 0.0 ) 
            IF( i > 1 ) THEN
               i = i - 1        !! first go down 
            ELSE
               j = j + 1        !! there if can't go further, go up.
            ENDIF
          ENDDO
          !! find the ratio in each of the layer in the reference prfoile
          !! redistribute the sum according to the ratio.
           ratio(i:j) = o3prof_base(i:j)/SUM( o3prof_base(i:j))
          o3prof(i:j) =       ratio(i:j)*SUM( o3prof(i:j)     )
       ENDIF
       status = OZT_S_SUCCESS
       RETURN
     END FUNCTION O3T_prof_check

END MODULE O3T_stnprof_class
