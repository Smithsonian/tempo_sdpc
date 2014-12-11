!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_apriori_class
! 
! determine a priori ozone profile for a given day and latitude
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
MODULE O3T_apriori_class
    USE O3T_stnprof_class
    USE OMI_LUN_set
    USE PGS_PC_class
    IMPLICIT NONE
    INTEGER (KIND = 4), PARAMETER, PRIVATE :: zero = 0
    INTEGER (KIND = 4), PARAMETER, PRIVATE :: ntoz = 10
    INTEGER (KIND = 4), PARAMETER, PRIVATE :: nlat_ctrs = 18
    
    REAL (KIND=4), DIMENSION(ntoz) :: toz
    REAL (KIND=4), DIMENSION(nlat_ctrs) :: lat_ctrs=(/-85.0,-75.0,-65.0,-55.0,&
    -45.0,-35.0,-25.0,-15.0,-5.0,5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0/)
    REAL (KIND=4), DIMENSION(nlat_ctrs) :: rlats
    REAL (KIND=4), DIMENSION(NLYR,ntoz,nlat_ctrs,13) :: tzaprf
    REAL (KIND=4), DIMENSION(NLYR,nlat_ctrs,12) :: climtm
    
   CONTAINS

     FUNCTION O3T_apriori_rd( ) RESULT( status )
       INCLUDE 'PGS_IO.f'
       INCLUDE 'PGS_IO_1.f'
       INTEGER (KIND=4), EXTERNAL :: pgs_io_gen_openf, pgs_io_gen_closef
       INTEGER (KIND=4) :: file_version, record_length, &
                           climoz_handle, climtm_handle
       INTEGER (KIND=4) :: status, ierr, ios
       INTEGER (KIND=4) :: kmonth, jlat, ioz, ilyr, mnth
       CHARACTER (LEN =255) :: msg

       record_length= 0

       file_version = 1
       status = PGS_IO_Gen_OpenF( O3_CLIM_LUN, PGSd_IO_Gen_RSeqFrm, &
                                  record_length, climoz_handle, &
                                  file_version )
       IF( status /= PGS_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( status,  "open clim oz file for read failed", &
                                 "O3T_apriori_rd", zero )
          status = OZT_E_FAILURE
          RETURN
       ELSE
          DO kmonth = 1, 12
            DO jlat = 1, nlat_ctrs
              READ( climoz_handle, "(1X,'Month = ',I2,' Latitude = ',F6.1)", &
                                     IOSTAT = ios ) mnth, rlats(jlat)
              IF( ios /= zero ) THEN
                 status = OZT_E_FAILURE
                 WRITE(msg,*) 'error: reading climoz file '
                 ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                       "O3T_apriori_rd", zero )
                 RETURN
              ENDIF

              DO ioz = 1, ntoz
                READ( climoz_handle, "(1X,F6.1,6F7.2,5F6.2)", IOSTAT = ios ) &
                      toz(ioz), ( tzaprf(ilyr,ioz,jlat,kmonth), ilyr = 1, NLYR )
                IF( ios /= zero ) THEN
                   status = OZT_E_FAILURE
                   WRITE(msg,*) 'error: reading climoz file '
                   ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                         "O3T_apriori_rd", zero )
                   RETURN
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          ! put month 1 in position 13 to handle wrap-around
          tzaprf(:,:,:,13) = tzaprf(:,:,:,1)
          status = pgs_io_gen_closef( climoz_handle )
       ENDIF
 
       file_version = 1
       status = PGS_IO_Gen_OpenF( TM_CLIM_LUN, PGSd_IO_Gen_RSeqFrm, &
                                  record_length, climtm_handle, &
                                  file_version )
       IF( status /= PGS_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( status,  "open clim tm file for read failed", &
                                 "O3T_apriori_rd", zero )
          status = OZT_E_FAILURE
          RETURN
       ELSE
          READ (climtm_handle,'(11F7.1)', IOSTAT = ios ) climtm

          IF( ios /= zero ) THEN
             status = OZT_E_FAILURE
             WRITE(msg,*) 'error: reading climtm file '
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_apriori_rd", zero )
             RETURN
          ENDIF
     
          status = pgs_io_gen_closef( climtm_handle )
       ENDIF
       status = OZT_S_SUCCESS
     END FUNCTION O3T_apriori_rd

     FUNCTION O3T_apriori_prf( latitude, jday, aprftm_k ) &
                               RESULT( status )
       REAL (KIND=4), INTENT(IN) :: latitude
       INTEGER (KIND=4), INTENT(IN) :: jday
       REAL (KIND=4), DIMENSION(:), INTENT(OUT) :: aprftm_k
       INTEGER (KIND=4) :: status, ierr, ios
       CHARACTER (LEN =255) :: msg
       INTEGER :: l1, l2, m1, m2
       REAL (KIND=4) :: fracl, xmon, fracm

       !
       ! -- compute indeces for bracketing months and latitudes
       !
       l1 = (latitude + 85.0) / 10.0 + 1.0
       IF(l1  <=  0) l1 = 1
       IF(l1  >= 18) l1 = 17
       l2 = l1 + 1
       fracl = (latitude - lat_ctrs(l1)) / (lat_ctrs(l2) - lat_ctrs(l1))
 
       xmon = (jday + 15.25) / 30.5
       m1 = xmon
       m2 = m1 + 1
       fracm = xmon - m1
       if(m1 == 0) m1 = 12
       if(m2 == 13) m2 = 1

       IF( SIZE(aprftm_k) < NLYR ) THEN
          status = OZT_E_FAILURE
          ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input array aprftm size too small", &
                                 "O3T_apriori_prf", zero )
          RETURN
       ENDIF
       aprftm_k=(1.0-fracl)*(1.0-fracm)*climtm(:,l1,m1) + &
                     fracl *(1.0-fracm)*climtm(:,l2,m1) + &
                     fracl *     fracm *climtm(:,l2,m2) + &
                (1.0-fracl)*     fracm *climtm(:,l1,m2)

       status = OZT_S_SUCCESS
     END FUNCTION O3T_apriori_prf
END MODULE O3T_apriori_class
