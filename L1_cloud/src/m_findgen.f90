!>Create double precision vector containing integer values 0,1,2...length-1
!
!-------------------------------------------------------------------------
!
! !ROUTINE:  findgen
! 
! !DESCRIPTION: 
!> Create double precision vector containing
!>               integer values 0,1,2,...length-1
!
! !CALLING SEQUENCE: 
!
!        vector = findgen(length)
!     
! !INPUT PARAMETERS:   
!> @param length[in] length of vector to create
!
! !OUTPUT PARAMETERS:  
!> @param vector[out] vector filled with 0,1,2,...length-1
!
! !SEE ALSO:  IDL documentation, indgen.f90
!
! !REVISION HISTORY: 
!
!> @author  13Aug97   Joiner     original code
!
!-------------------------------------------------------------------------
module m_findgen

  private
  public findgen

contains 

  function findgen(length) result (vector)

    implicit none
    ! !INPUT PARAMETERS:   
    integer, intent(in)            :: length 
    ! !OUTPUT PARAMETERS:  
    real (KIND=8), dimension(length) :: vector
    !local variables
    !---------------
    integer                                :: i

    vector = (/ (i,i=0,length-1) /)
  end function findgen

end module m_findgen

