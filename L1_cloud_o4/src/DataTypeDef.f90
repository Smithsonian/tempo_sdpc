!*****************
MODULE DataTypeDef
!*****************

! PURPOSE: 
!   Define common data types, some useful math constants and status codes.
!   It is borrowed from "Nuberical Recipes" and is used by all of my Fortran
!   routines to ensure consistency across the board.  
!
!   Note: treat these as reserved keywords, avoid using them as variable names!

IMPLICIT NONE
SAVE

! Integer data types:

INTEGER,PARAMETER::I8B=SELECTED_INT_KIND(18) ! (+/-) 9223372036854775807
INTEGER,PARAMETER::I4B=SELECTED_INT_KIND(9)  ! (+/-)          2147483647
INTEGER,PARAMETER::I2B=SELECTED_INT_KIND(4)  ! (+/-)               32767
INTEGER,PARAMETER::I1B=SELECTED_INT_KIND(2)  ! (+/-)                 127

! Floating point data types: 

! SP=SELECTED_REAL_KIND(6): (+/-) 3.4028235E+38
!   =SELECTED_REAL_KIND(6,37) 
!   =KIND(1.0)
INTEGER,PARAMETER::SP=SELECTED_REAL_KIND(6,37)

! DP=SELECTED_REAL_KIND(15): (+/-) 1.7976931348623167E+308
!   =SELECTED_REAL_KIND(15,307)
!   =KIND(1.0D0)  
INTEGER,PARAMETER::DP=SELECTED_REAL_KIND(15,307)

! Complex data type:

INTEGER,PARAMETER::SPC=KIND((1.0,1.0))
INTEGER,PARAMETER::DPC=KIND((1.0D0,1.0D0))

!*********************
END MODULE DataTypeDef
!*********************
