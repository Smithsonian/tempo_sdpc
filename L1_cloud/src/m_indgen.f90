module m_indgen

public indgen

contains

function indgen(length) result (vector)
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  indgen
! 
! !DESCRIPTION: similar to IDL's "indgen" function
!
! !CALLING SEQUENCE: 
!
!        vector = indgen(length)
!     
! !INPUT PARAMETERS:   
integer, intent(in)               :: length 
!                        length : length of vector to create
!
! !OUTPUT PARAMETERS:  
integer, dimension(length) :: vector
!                        vector : vector filled with 0,1,2,...length-1
!
! !SEE ALSO:  IDL documentation, findgen.f90
!
! !REVISION HISTORY: 
!
!  13Aug97   Joiner     original code
!
!EOP
!-------------------------------------------------------------------------

!local variables
!---------------
integer                                :: i

vector = (/ (i,i=0,length-1) /)

end function indgen

end module m_indgen
