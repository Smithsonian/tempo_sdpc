!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_const
!
! store of a number of constants used to be in data file
! CONST.EP, CONST.N7, CONST.M3 of the TOMSV8 codes.
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

MODULE O3T_const
    IMPLICIT NONE
    REAL (KIND=4), DIMENSION(3) :: f360 = &         ! aerosol adjustment factor
                                           (/ 4.223, -1.460, 0.145 /) 
    REAL (KIND=4) :: f313    =-8.0  !profile shape adjustment factor
    REAL (KIND=4) :: swthrsh = 1.25 !threshold ratio for hi sza !Junsung: need to be checked (1.25 or other values)
    REAL (KIND=4) :: soilimEP= 12.5 !so2 flag limit
    REAL (KIND=4) :: soilimN7= 25.0 !so2 flag limit
    REAL (KIND=4) :: soilimM3= 24.0 !so2 flag limit
    REAL (KIND=4) :: soilimAD= 12.5 !so2 flag limit
    REAL (KIND=4) :: soilimOM= 12.5 !so2 flag limit
    REAL (KIND=4) :: glntlm  =-1.5  !R360-R331 glint flag limit 
    REAL(KIND=4), DIMENSION(3) :: flg3lm = (/ 5.0, 2.5, 99.0 /) !flag 3 limits
    REAL(KIND=4), DIMENSION(3) :: flg4lm = (/  3.5, 2.0,  5.0 /) !flag 4 limits
END MODULE O3T_const
