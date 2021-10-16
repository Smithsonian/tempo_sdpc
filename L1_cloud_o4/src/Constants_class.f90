!*********************
MODULE Constants_class
!*********************

use DataTypeDef
implicit none

!...............................................................................
! return status codes: more or less aligned with HDF5 status codes
!...............................................................................

INTEGER,PARAMETER::FAILURE_STATE=-1
INTEGER,PARAMETER::SUCCESS_STATE=0
INTEGER,PARAMETER::WARNING_STATE=1

!...............................................................................
! unsigned integers (as defined in OMI):
!...............................................................................

integer(I4B),public,parameter::FILLVALUE_UI1B=255    ! = 2^8  - 1
integer(I4B),public,parameter::FILLVALUE_UI2B=65535  ! = 2^16 - 1

!...............................................................................
! signed integers:
!...............................................................................

integer(I1B),public,parameter::FILLVALUE_I1B = - HUGE(1_I1B)
integer(I2B),public,parameter::FILLVALUE_I2B = - HUGE(1_I2B)
integer(I4B),public,parameter::FILLVALUE_I4B = - HUGE(1_I4B) 

!...............................................................................
! float type (as defined in OMI):
!...............................................................................

real(SP),public,parameter::FILLVALUE_SP=-2.0**100    
real(DP),public,parameter::FILLVALUE_DP=-2.0d0**100  

!...............................................................................
! Math constants:
!...............................................................................

REAL(SP),PARAMETER::PI     =3.141592653589793238462643383279502884197_sp
REAL(SP),PARAMETER::PIO2   =1.57079632679489661923132169163975144209858_sp
REAL(SP),PARAMETER::TWOPI  =6.283185307179586476925286766559005768394_sp
REAL(SP),PARAMETER::SQRT2  =1.41421356237309504880168872420969807856967_sp
REAL(SP),PARAMETER::EULER  =0.5772156649015328606065120900824024310422_sp
REAL(DP),PARAMETER::PI_D   =3.141592653589793238462643383279502884197_dp
REAL(DP),PARAMETER::PIO2_D =1.57079632679489661923132169163975144209858_dp
REAL(DP),PARAMETER::TWOPI_D=6.283185307179586476925286766559005768394_dp

!*************************
END MODULE Constants_class
!*************************
