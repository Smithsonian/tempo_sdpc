!>Vector averaging routines
! 
!-------------------------------------------------------------------------
!
! !ROUTINE:  i_avg, r_avg
! 
! !DESCRIPTION: 
!> similar to IDL avg function (from code 916 userlib)
!! finds average of a vector (not a matrix)
!
! !CALLING SEQUENCE: 
!
!        avg = i_avg(vectorin)
!     
! !INPUT PARAMETERS:   
!> @param vectorin  input vector to average
!
! !OUTPUT PARAMETERS:  
!> @return avg average of vectorin
!
! !SEE ALSO:  IDL documentation 
!
! !REVISION HISTORY: 
!
! @author  13Aug96   Joiner     Original code.
!----------------------------------------------
module m_avg

  interface avg
    module procedure r_avg
    module procedure r4_avg
    module procedure i_avg
    module procedure r_avg2D
  end interface

contains

  !>Average integer vector
  function i_avg(vectorin) result (avg)

    implicit none
    ! !INPUT PARAMETERS:   
    integer, dimension(:), intent(in)   :: vectorin
    ! !OUTPUT PARAMETERS:  
    real (KIND=8)                       :: avg

    avg=sum(vectorin)/real(size(vectorin))

  end function i_avg

  !>Average real_8 vector
  function r_avg(vectorin) result (avg)
    implicit none

    real (KIND=8), dimension(:), intent(in)  :: vectorin
    real (KIND=8)                            :: avg

    avg=sum(vectorin)/size(vectorin)

  end function r_avg

  !>Average real_4 vector
  function r4_avg(vectorin) result (avg)
    implicit none

    real (kind=4), dimension(:), intent(in)  :: vectorin
    real (kind=4)                            :: avg

    avg=sum(vectorin)/size(vectorin)

  end function r4_avg

  !>Average real_8 2D matrix
  function r_avg2D(vectorin) result (avg)
    implicit none

    real (KIND=8), dimension(:,:), intent(in) :: vectorin
    real (KIND=8)                             :: avg

    avg=sum(reshape(vectorin,(/size(vectorin,1)*size(vectorin,2) /) )) &
         / (size(vectorin,1) * size(vectorin,2))

  end function r_avg2D

end module m_avg
