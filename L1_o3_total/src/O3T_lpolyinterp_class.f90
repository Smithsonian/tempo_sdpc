!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE O3T_lpolyinterp_class
!  
!  This module contains functions that perform the Lagrange interpolation
!  of the N-value components.
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
MODULE O3T_lpolyinterp_class
  USE OMI_SMF_class
  USE O3T_pixel_class
  USE O3T_nval_class
  USE O3T_dndx_class
  USE O3T_lpolycoef_class

  IMPLICIT NONE
  INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0

  PUBLIC :: O3T_lpoly_interpPLW

  CONTAINS  

!   This function performs solar and view angles( Log of 1/COS(angles) ) 
!   Lagrangian interpolations for all (4) the pressure levels (P), all (12) 
!   the layers (L), and all the wavelength (W) index range from [iwl_s, iwl_e], 
!   of the table parameters ezero, tr and sb.
!!Input Parameters:
!   pixGEO  : the geometry inforamtion for the pixel to the output
!             parameters are associated with. pixGEO has angular and
!             terrain and cloud pressure info.
!   coefs   : Lagrangian interpolation coefficients associated with
!             pixGEO
!!Output Parameters:
!
!   ezero(nlyrT,npresT,nwl_sub), tr(nlyrT,npresT,nwl_sub), sb(nlyrT,npresT,nwl_sub) 
!

    FUNCTION O3T_lpoly_interpPLW( iprof, pixGEO, coefs, &
                                  ezero, tr, sb  ) RESULT( status )
      INTEGER, INTENT(IN) :: iprof
      TYPE (O3T_pixgeo_type), INTENT(IN) :: pixGEO
      TYPE (O3T_lpoly_coef_type), INTENT(IN) :: coefs
      REAL (KIND=4), DIMENSION(:,:,:), INTENT(OUT) :: ezero, tr, sb
      REAL (KIND=4), DIMENSION(SIZE(sb,1), SIZE(sb,2), SIZE(sb,3)):: &
                                 zone, ztwo, inot, ione, itwo, tran 
      REAL (KIND=8) :: i0_s, z1_s, z2_s, tr_s, y 
      INTEGER (KIND=4) :: isza_b, ivza_b 
      INTEGER (KIND=4) :: ilyr, ipres, iwl
      INTEGER (KIND=4) :: LL0, LLp0, LLpw0, LLpwl0, LL, is, iv, m
      INTEGER (KIND=4) :: L20, L2p0, L2pw0, L2
      INTEGER (KIND=4) :: status, ierr
      CHARACTER (LEN =255) :: msg
 
      IF( .NOT. dndx_read ) THEN
         status = OZT_E_FAILURE
         WRITE(msg,*) "DNDX LUT has not been read in yet"
         ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_lpoly_interpPLW", zero )
         RETURN
      ENDIF

      isza_b = coefs%isza
      ivza_b = coefs%ivza

      LL0 = (iprof-1)*d6size(2)+(isza_b-1)*d6size(6) + ivza_b
      L20 = (iprof-1)*d4size(2)
      DO ipres = 1, npresT
        LLp0 = LL0+(ipres-1)*d6size(3)  
        L2p0 = L20+(ipres-1)*d4size(3)
        DO iwl = 1, nwl_sub
          LLpw0 = LLp0 + (iwl-1)*d6size(4)
          L2pw0 = L2p0 + (iwl-1)*d4size(4)
          DO ilyr = 1, nlyrT
            LLpwl0 = LLpw0+(ilyr-1)*d6size(5)
            m = 0
            i0_s = 0.0
            z1_s = 0.0
            z2_s = 0.0
            tr_s = 0.0
            DO is = 1, NINTERP
              LL = LLpwl0+(is-1)*d6size(6)
              DO iv = 1, NINTERP
                m = m + 1
                y = coefs%Lnk(m)
                i0_s = i0_s +  lgi0_d(LL)*y
                z1_s = z1_s +    z1_d(LL)*y
                z2_s = z2_s +    z2_d(LL)*y
                tr_s = tr_s +    tr_d(LL)*y
                LL   = LL + 1
              ENDDO
            ENDDO
            inot(ilyr,ipres,iwl) = real(i0_s, kind=4)
            zone(ilyr,ipres,iwl) = real(z1_s, kind=4)
            ztwo(ilyr,ipres,iwl) = real(z2_s, kind=4)
            tran(ilyr,ipres,iwl) = real(tr_s, kind=4)
            L2 = L2pw0 + ilyr
              sb(ilyr,ipres,iwl) = sb_d(L2)
          ENDDO
        ENDDO
      ENDDO 

       inot = 10.0** inot
       ione = real(zone * inot * pixGEO%p1 , kind=4)
       itwo = real(ztwo * inot * pixGEO%pr * pixGEO%p1 , kind=4)
       tr   = tran * inot
       ezero = real(inot + ione*pixGEO%cphi + itwo*pixGEO%c2phi , kind=4)
       status = OZT_S_SUCCESS
    END FUNCTION O3T_lpoly_interpPLW

END MODULE O3T_lpolyinterp_class
