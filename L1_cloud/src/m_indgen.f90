!>Create integer vector containing values 0,1,2,...length-1
!
!-------------------------------------------------------------------------
!
! !ROUTINE:  indgen
! 
! !DESCRIPTION: create integer vector containing values
!               0,1,2,...length-1
!
! !CALLING SEQUENCE: 
!
!        vector = indgen(length)
!     
! !INPUT PARAMETERS:   
!> @param  length[in]  length of vector to create
!
! !OUTPUT PARAMETERS:  
!> @param vector[out] vector filled with 0,1,2,...length-1
!
! !SEE ALSO:  IDL documentation, findgen.f90
!
! !REVISION HISTORY: 
!
!> @author  13Aug97   Joiner     original code
!
!-------------------------------------------------------------------------
module m_indgen

  private
  public indgen

contains

  function indgen(length) result (vector)
    implicit none
    ! !INPUT PARAMETERS:   
    integer, intent(in)               :: length 
    ! !OUTPUT PARAMETERS:  
    integer, dimension(length) :: vector
    !local variables
    !---------------
    integer                                :: i

    vector = (/ (i,i=0,length-1) /)

  end function indgen

end module m_indgen
