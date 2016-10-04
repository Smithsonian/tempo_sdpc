!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE O3T_iztrsb_m
!  
!  This module contain function for forward calculation of the components that
!  made up the nvalue.
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
MODULE O3T_iztrsb_m
    USE OMI_SMF_class 
    USE O3T_nval_class, ONLY :  nvalLUT_t
    USE O3T_pixel_class, ONLY : O3T_pixgeo_type
    USE O3T_lpolycoef_class, ONLY: O3T_lpoly_coef_type, NINTERP
    use tell_module

    IMPLICIT NONE
    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0

    PUBLIC  :: O3T_iztrsb
!    PUBLIC  :: O3T_iztrsbp

    CONTAINS

!!Description:
! FUNCTION O3T_iztrsb
!   This function performs angles and pressure interpolations for one 
!   wavelength index one one ozone profil index of the table parameters 
!   ezero tr and sb, which are used by the calling routines to calculate
!   table radiances or n-values. uses lagrangian interpolation for angles
!   and pressures except for low pressure (high clouds) for which a 
!   linear extrapolation is used.
!!Input Parameters:
!   iwl     : the wavelength index (index into the LUT)
!   iprof   : the index of the ozone profile (1, 21)
!   pixGEO  : the geometry inforamtion for the pixel to the output
!             parameters are associated with. pixGEO has angular and
!             terrain and cloud pressure info.
!   coefs   : Lagrangian interpolation coefficients associated with
!             pixGEO
!!Output Parameters:
!   
!    ezgr, tgr, sbgr : ez, tg, ans sb for the ground terrain pressure
!    ezcl, tcl, sbcl : ez, tg, ans sb for the cloud top pressure
!!Return
!   satus : The SMF stauts 
!!References and Credits
! Written by
! Kai Yang
! University of Maryland Baltimore County
!!END Description:

      FUNCTION O3T_iztrsb( nv, iwl, iprof, pixGEO, coefs, &
                           ezgr, tgr, sbgr, knbgr, &
                           ezcl, tcl, sbcl, knbcl ) RESULT( errstat)

        implicit none

        INTEGER (KIND=4), INTENT(IN) :: iwl, iprof
        TYPE (nvalLUT_t), INTENT(IN) :: nv
        TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
        REAL (KIND=4), INTENT(OUT) :: ezgr, tgr, sbgr, knbgr, &
                                      ezcl, tcl, sbcl, knbcl 
        INTEGER (KIND=4) :: errstat
        INTEGER :: ipres
        REAL (KIND=4), DIMENSION(nv%npres) :: ezofp, tofp, sbofp, knbofp
        REAL (KIND=8) :: fac

        errstat = 0

        DO ipres = 1, nv%npres 
          errstat = O3T_lpoly_interp1( nv, ipres, iwl, iprof, pixGEO, &
                              coefs, ezofp(ipres), tofp(ipres), sbofp(ipres), &
                              knbofp(ipres) )
          IF( errstat .NE. 0 ) THEN
            call tell_error(tell_runtime_error, &
                 "O3T_iztrsb: interpolation error" , errstat)
             RETURN
          ENDIF
        ENDDO 

        ezgr  = real(SUM(  ezofp*coefs%pgr ) , kind=4)
         tgr  = real(SUM(   tofp*coefs%pgr ) , kind=4)
        sbgr  = real(SUM(  sbofp*coefs%pgr ) , kind=4)
        knbgr = real(SUM( knbofp*coefs%pgr ) , kind=4)
        IF( pixGEO%pc >= 0.25 ) THEN 
           ezcl = real(SUM(  ezofp*coefs%pcl ) , kind=4)
            tcl = real(SUM(   tofp*coefs%pcl ) , kind=4)
           sbcl = real(SUM(  sbofp*coefs%pcl ) , kind=4)
          knbcl = real(SUM( knbofp*coefs%pcl ) , kind=4)
        ELSE                                            ! linear extrapolation
           fac  = (pixGEO%log_pc-nv%presLog(3)) &
                 /(nv%presLog(4)-nv%presLog(3))
           ezcl = real(fac*( ezofp(4)- ezofp(3))+ ezofp(3) , kind=4)
            tcl = real(fac*(  tofp(4)-  tofp(3))+  tofp(3) , kind=4)
           sbcl = real(fac*( sbofp(4)- sbofp(3))+ sbofp(3) , kind=4)
          knbcl = real(fac*(knbofp(4)-knbofp(3))+knbofp(3) , kind=4)
        ENDIF

      END FUNCTION O3T_iztrsb

!      FUNCTION O3T_iztrsbp( nv, iwl, iprof, pixGEO, coefs, &
!                            ezgrp, tgrp, sbgrp, knbgrp, &
!                            ezclp, tclp, sbclp, knbclp ) RESULT( status)
!        INTEGER (KIND=4), INTENT(IN) :: iwl, iprof
!        TYPE (nvalLUT_t), INTENT(IN) :: nv
!        TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
!        TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
!        REAL (KIND=4), INTENT(OUT) :: ezgrp, tgrp, sbgrp, knbgrp, &
!                                      ezclp, tclp, sbclp, knbclp 
!        INTEGER (KIND=4) :: status, ierr
!        INTEGER :: ipres
!        REAL (KIND=4), DIMENSION(nv%npres) :: ezofp, tofp, sbofp, knbofp
!        REAL (KIND=8) :: fac
!
!        DO ipres = 1, nv%npres 
!          status = O3T_lpoly_interp1( nv, ipres, iwl, iprof, pixGEO, &
!                              coefs, ezofp(ipres), tofp(ipres), sbofp(ipres), &
!                              knbofp(ipres) )
!          IF( status .NE. OZT_S_SUCCESS ) THEN
!             ierr = OMI_SMF_setmsg( OZT_E_INPUT, "interpolation error", &
!                                   "O3T_iztrsb", zero )
!             RETURN
!          ENDIF
!        ENDDO 
!
!        ezgrp  = SUM(  ezofp*coefs%pgrp ) 
!         tgrp  = SUM(   tofp*coefs%pgrp ) 
!        sbgrp  = SUM(  sbofp*coefs%pgrp ) 
!        knbgrp = SUM( knbofp*coefs%pgrp )
!        IF( pixGEO%pcp >= 0.25 ) THEN 
!           ezclp = SUM(  ezofp*coefs%pclp ) 
!            tclp = SUM(   tofp*coefs%pclp ) 
!           sbclp = SUM(  sbofp*coefs%pclp ) 
!          knbclp = SUM( knbofp*coefs%pclp )
!        ELSE                                            ! linear extrapolation
!           fac  = (pixGEO%log_pcp-nv%presLog(3)) &
!                 /(nv%presLog(4) -nv%presLog(3))
!           ezclp = fac*( ezofp(4)- ezofp(3))+ ezofp(3)
!            tclp = fac*(  tofp(4)-  tofp(3))+  tofp(3)
!           sbclp = fac*( sbofp(4)- sbofp(3))+ sbofp(3)
!          knbclp = fac*(knbofp(4)-knbofp(3))+knbofp(3)
!        ENDIF
!        status = OZT_S_SUCCESS
!      END FUNCTION O3T_iztrsbp

!!Description:
!   FUNCTION O3T_lpoly_interp1
!   This function performs angles( Log of 1/COS(angles) ) Lagrangian
!   interpolations for one wavelength index, one ozone profile
!   index, and one pressure level index of the table parameters
!   ezero, tr and sb.
!!Input Parameters:
!   ipres   : the pressure index
!   iwl     : the wavelength index (index into the LUT)
!   iprof   : the index of the ozone profile [1:21] 
!   pixGEO  : the geometry inforamtion for the pixel to the output
!             parameters are associated with. pixGEO has angular and
!             terrain and cloud pressure info.
!   coefs   : Lagrangian interpolation coefficients associated with
!             pixGEO
!!Output Parameters:
!
!    ezero, tr, sb : at the wavelength, profile, and  pressure level indices
!!Return
!   satus : The SMF stauts
!!References and Credits
! Written by
! Kai Yang
! University of Maryland Baltimore County
!!END Description:

    FUNCTION O3T_lpoly_interp1( nv, ipres, iwl, iprof, pixGEO, &
                                coefs, ezero, tr, sb, knb ) RESULT( errstat )

      implicit none

      INTEGER (KIND=4), INTENT(IN) :: ipres, iwl, iprof
      TYPE (nvalLUT_t), INTENT(IN) :: nv
      TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
      TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
      REAL (KIND=8) :: inot, ione, itwo, tran, knr2, y
      REAL (KIND=8) :: zone, ztwo
      REAL (KIND=4), INTENT(OUT) :: ezero, tr, sb, knb
      INTEGER (KIND=4) :: isza_b, ivza_b, LL0, LL, iv, is, m 
      INTEGER (KIND=4) :: errstat
      CHARACTER (LEN =255) :: msg

      errstat = 0

      IF( .NOT. nv%nval_read ) THEN
         WRITE(msg,*) "O3T_lpoly_interp1: NVAL LUT has not been read in yet"
         call tell_error(tell_application_error, msg, errstat)
         RETURN
      ENDIF

      IF( iprof < 1 .OR. iprof > nv%nprof ) THEN
         WRITE(msg,*) "O3T_lpoly_interp1: iprof = ", iprof, &
              "out of range[1,",nv%nprof,"]"
         call tell_error(tell_application_error, msg, errstat)
         RETURN
      ENDIF
      IF( iwl < 1 .OR. iwl > nv%nwl ) THEN
         WRITE(msg,*) "O3T_lpoly_interp1: iwl = ", iwl, &
              "out of range[1,",nv%nwl,"]"
         call tell_error(tell_application_error, msg, errstat)
         RETURN
      ENDIF
      IF( ipres < 1 .OR. ipres > nv%npres ) THEN
         WRITE(msg,*) "O3T_lpoly_interp1: ipres = ", ipres, &
              "out of range[1,",nv%npres,"]"
         call tell_error(tell_application_error, msg, errstat)
         RETURN
      ENDIF

      isza_b = coefs%isza
      ivza_b = coefs%ivza

      LL0 = (ipres-1)*nv%n5size(2)+(iwl-1)*nv%n5size(3) &
           +(iprof-1)*nv%n5size(4)+(isza_b-1)*nv%n5size(5)+ivza_b
      
      m = 0
      inot = 0.0
      zone = 0.0
      ztwo = 0.0
      tran = 0.0
      knr2 = 0.0
      DO is = 1, NINTERP
        LL = LL0+(is-1)*nv%n5size(5)
        DO iv = 1, NINTERP
          m = m + 1
          y = coefs%Lnk(m)
          inot = inot +  nv%lgi0(LL)*y
          zone = zone +    nv%z1(LL)*y
          ztwo = ztwo +    nv%z2(LL)*y
          tran = tran +    nv%tr(LL)*y
          knr2 = knr2 +   nv%knb(LL)*y
          LL   = LL + 1
        ENDDO
      ENDDO
 
      inot = 10.0** inot
      ione = zone * pixGEO%p1 * inot
      itwo = ztwo * pixGEO%pr * pixGEO%p1 * inot
      tr   = real((tran * inot), kind=4) 
      knb  = real(knr2, kind=4)
 
      ezero = real((inot + ione*pixGEO%cphi + itwo*pixGEO%c2phi), kind=4)
      LL = (ipres-1)*nv%n3size(2)+(iwl-1)*nv%n3size(3)+iprof
      sb    = nv%sb(LL)

    END FUNCTION O3T_lpoly_interp1

END MODULE O3T_iztrsb_m
