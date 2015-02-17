!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_so2_class 
! 
! This module compute SO2 index from residues and pre-stored coefficients
! for OMI, GOMI, and for the different TOMS instruments on different 
! satellite platforms.
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
MODULE O3T_so2_class 
    USE OMI_SMF_class
    USE O3T_const
    IMPLICIT NONE
    INTEGER, PARAMETER, PRIVATE :: zero = 0
    INTEGER (KIND=4), DIMENSION(4), PRIVATE :: &
                                      iso2wn7 = (/1, 2, 3, 5/), &
                                      iso2wm3 = (/1, 2, 3, 5/), &
                                      iso2wad = (/3, 4, 5, 6/), &
                                      iso2wep = (/3, 4, 5, 6/), &
                                      iso2wqt = (/3, 4, 5, 6/)
    INTEGER (KIND=4), DIMENSION(6), PRIVATE :: iwl = (/1,2,3,4,5,6/)
    REAL (KIND=4), DIMENSION(6), PRIVATE :: wlep = &
             (/308.7, 312.61, 317.62, 322.42, 331.34, 360.15/)
    REAL (KIND=4), DIMENSION(5), PRIVATE :: &
               o3an7  = (/1.684, 0.877, 0.145, 0.022, 0.000000001/), &
               so2an7 = (/3.942, 2.230, 0.055, 0.006, 0.0        /), &
               o3am3  = (/1.679, 0.867, 0.143, 0.020, 0.000000001/), &
               so2am3 = (/3.880, 2.220, 0.025, 0.011, 0.0        /), &
               o3aad  = (/2.956, 1.630, 0.877, 0.474, 0.137      /), &
               so2aad = (/9.090, 4.485, 2.384, 0.643, 0.050      /), &     
               o3aep  = (/2.970, 1.640, 0.880, 0.470, 0.140      /), &
               so2aep = (/9.620, 4.360, 2.210, 0.650, 0.050      /), &
               o3aqt  = (/2.987, 1.669, 0.883, 0.471, 0.142      /), &
               so2aqt = (/9.214, 4.238, 2.267, 0.727, 0.042      /)

    PUBLIC :: O3_so2_setCoef
    PUBLIC :: O3_so2_index

    CONTAINS
      !! this function set the coefficients depending on the input
      !! satellite name.
      FUNCTION O3_so2_setCoef( satname, o3abs, so2abs, iso2w, wlen, soilim ) &
                               RESULT( status )
        CHARACTER (LEN=*), INTENT(IN) :: satname 
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: wlen
        INTEGER (KIND=4), DIMENSION(4), INTENT(OUT) :: iso2w
        REAL (KIND=4), DIMENSION(5), INTENT(OUT) :: o3abs, so2abs
        REAL (KIND=4), INTENT(OUT) :: soilim
        INTEGER :: di, i_foo(1)
        INTEGER (KIND=4) :: ierr, status
        
        status = OZT_S_SUCCESS
        SELECT CASE( satname ) 
          CASE ('ep', 'EP' )
            o3abs  = o3aep 
            so2abs = so2aep 
            iso2w  = iso2wep
            soilim = soilimEP
            RETURN 
          CASE ('n7', 'N7' )
            o3abs  = o3an7 
            so2abs = so2an7 
            iso2w  = iso2wn7
            soilim = soilimN7
            RETURN 
          CASE ('m3', 'M3' )
            o3abs  = o3am3 
            so2abs = so2am3 
            iso2w  = iso2wm3
            soilim = soilimM3
            RETURN 
          CASE ('ad', 'AD' )
            o3abs  = o3aad 
            so2abs = so2aad 
            iso2w  = iso2wad
            soilim = soilimAD
            RETURN 
          CASE ('OM', 'om' )
            o3abs  = o3aep 
            so2abs = so2aep 
            iso2w  = iso2wep
            soilim = soilimOM
            IF( wlep(1) > MAXVAL( wlen ) .OR. &
                wlep(4) < MINVAL( wlen ) ) THEN
               status = OZT_E_FAILURE
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, "wlT_min too large", &
                                     "O3_so2_setCoef", zero )
               RETURN
            ENDIF
            DO di = 1, 6
              i_foo = MINLOC( ABS(wlen-wlep(di)) )
              iwl(di) = i_foo(1)
            ENDDO  
          CASE DEFAULT
            status = OZT_E_FAILURE
            RETURN 
        END SELECT 

      END FUNCTION O3_so2_setCoef

      ! compute soi from residues and pre-stored coefficients
      FUNCTION O3_so2_index( wlen, dndomega, residue, iso2w, &
                             o3abs, so2abs ) RESULT( so2ind ) 
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: wlen   !wavelengths
        INTEGER (KIND=4), DIMENSION(4), INTENT(IN) :: iso2w !soi channels
        REAL (KIND=4), DIMENSION(5), INTENT(IN) :: o3abs, & !o3 absorp. coef.
                                                   so2abs   !so absorp. coef.

        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: dndomega, &  !sensitivities
                                                   residue    !input v8 residues
        REAL (KIND=4) :: so2ind

        !- internal parameters
        REAL (KIND=4), DIMENSION(5) :: abs_rat, dwav, snso
        REAL (KIND=4), DIMENSION(6) :: sens, resn
        REAL (KIND=4) :: sens13, sens23, snso13, snso23, resn13, resn23
        !INTEGER (KIND=4) :: di

        abs_rat(1:5) = so2abs(1:5)/o3abs(1:5)
        ! calculate separation between wavls and reflectivity wavl.
        ! also calculate so2 sensitivities
        dwav(1:5)    = wlen(iwl(1:5))-wlen(iwl(iso2w(4)))
        snso(1:5)    = dndomega(iwl(1:5))*abs_rat(1:5) 

        resn(1:6)  = residue( iwl(1:6) )
        sens(1:6)  = dndomega( iwl(1:6) )
      
        ! calculate "cross" sensitivities and residues for wavelengths
        ! used in the retrieval 

        sens13 = sens(iso2w(1))*dwav(iso2w(3))-sens(iso2w(3))*dwav(iso2w(1))
        sens23 = sens(iso2w(2))*dwav(iso2w(3))-sens(iso2w(3))*dwav(iso2w(2))
        snso13 = snso(iso2w(1))*dwav(iso2w(3))-snso(iso2w(3))*dwav(iso2w(1))
        snso23 = snso(iso2w(2))*dwav(iso2w(3))-snso(iso2w(3))*dwav(iso2w(2))
        resn13 = resn(iso2w(1))*dwav(iso2w(3))-resn(iso2w(3))*dwav(iso2w(1))
        resn23 = resn(iso2w(2))*dwav(iso2w(3))-resn(iso2w(3))*dwav(iso2w(2))

        ! calculate so2 index
        so2ind = (resn23*sens13-resn13*sens23)/(snso23*sens13-snso13*sens23)
        IF( so2ind > 800.0 ) so2ind = 800.0 

      END FUNCTION O3_so2_index

END MODULE O3T_so2_class
