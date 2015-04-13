!>Routines to find standard deviation of vectors
!
!-------------------------------------------------------------------------
!
! !DESCRIPTION:  finds std dev of a vector (not a matrix)
!
! !CALLING SEQUENCE:
!
!        std = sigma(vectorin)
!
!> @param   ar1[in]  input vector to find std deviation
!> @param   avg1[out] optional output of vector average
!> @param    std[out]        std deviation

! !SEE ALSO:  m_avg
!
! !REVISION HISTORY:
!
!> @author  21Oct97   M. Karki     Original code.
!
!-------------------------------------------------------------------------
module m_sigma

  interface sigma
    module procedure r_sigma
    module procedure r4_sigma
  end interface

contains

  !>Standard dev of a real_8 vector
  function  r_sigma(ar1, avg1) result(std)

    use m_avg
    implicit none

    !-------------------------------------------------------------------------
    ! !INPUT PARAMETERS:
    real (KIND=8), dimension(:), intent(in)  :: ar1
    ! !OUTPUT PARAMETERS:
    real (KIND=8), optional, intent(out)     :: avg1
    real (KIND=8)               :: std
    !-------------------------------------------------------------------------

    !! LOCAL PARAMETERS:
    real (KIND=8) :: avgarg

    ! note this is the correct standard deviation when estimating
    ! the mean
    !    std = sqrt(sum( (ar1-avg(ar1))**2) /float(size(ar1)-1) )
    ! this is the function evaluated in the IDL sigma function

    avgarg=avg(ar1)

    !     std = sqrt(sum( (ar1-avgarg)**2) /float(size(ar1)) )
    std = sqrt(sum( (ar1-avgarg)**2) /real(size(ar1), kind=8) )

    if (present(avg1)) avg1=avgarg

  end function r_sigma



  !>Standard dev of a real_4 vector
  function  r4_sigma(ar1, avg1) result(std)

    use m_avg
    implicit none

    !-------------------------------------------------------------------------
    ! !INPUT PARAMETERS:
    real (kind=4), dimension(:), intent(in)  :: ar1
    ! !OUTPUT PARAMETERS:
    real (kind=4), optional, intent(out)     :: avg1
    real (kind=4)                :: std
    !-------------------------------------------------------------------------

    !! LOCAL PARAMETERS: none
    real (KIND=4) :: avgarg

    ! note this is the correct standard deviation when estimating
    ! the mean
    !    std = sqrt(sum( (ar1-avg(ar1))**2) /float(size(ar1)-1) )
    ! this is the function evaluated in the IDL sigma function

    avgarg=avg(ar1)

    !     std = sqrt(sum( (ar1-avgarg)**2) /float(size(ar1)) )
    std = sqrt(sum( (ar1-avgarg)**2) /real(size(ar1), kind=4) )

    if (present(avg1)) avg1=avgarg

  end function r4_sigma

end module m_sigma
