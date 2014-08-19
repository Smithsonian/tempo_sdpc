module m_findgen

  public 

contains 

  function findgen(length) result (vector)
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  findgen
    ! 
    ! !DESCRIPTION: create double precision vector containing
    !               values 0,1,2,...length-1
    !
    ! !CALLING SEQUENCE: 
    !
    !        vector = findgen(length)
    !     
    ! !INPUT PARAMETERS:   
    integer, intent(in)            :: length 
    !                        length : length of vector to create
    !
    ! !OUTPUT PARAMETERS:  
    real (KIND=8), dimension(length) :: vector
    !                        vector : vector filled with 0,1,2,...length-1
    !
    ! !SEE ALSO:  IDL documentation, indgen.f90
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
  end function findgen

end module m_findgen

