module m_sigma

 interface sigma
        module procedure r_sigma
        module procedure r4_sigma
 end interface

contains

     function  r_sigma(ar1, avg1, iprt) result(std)

     use m_avg
     implicit none

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sigma 
!
! !DESCRIPTION:  finds std dev of a vector (not a matrix)
!
! !CALLING SEQUENCE:
!
!        std = sigma(vectorin)
!
! !INPUT PARAMETERS:
     real (KIND=8), dimension(:)  :: ar1
!                       ar1 : input vector to find std deviation
     real (KIND=8),    optional   :: avg1
     integer, optional   :: iprt
!
! !OUTPUT PARAMETERS:
     real (KIND=8)                :: std
!                       std      :  std deviation
!
! !SEE ALSO:  m_avg
!
! !REVISION HISTORY:
!
!  21Oct97   M. Karki     Original code.
!EOP
!-------------------------------------------------------------------------

!! LOCAL PARAMETERS: none

! note this is the correct standard deviation when estimating
! the mean
!    std = sqrt(sum( (ar1-avg(ar1))**2) /float(size(ar1)-1) )
! this is the function evaluated in the IDL sigma function
    real (KIND=8) :: avgarg

     avgarg=avg(ar1)

!     std = sqrt(sum( (ar1-avgarg)**2) /float(size(ar1)) )
     std = sqrt(sum( (ar1-avgarg)**2) /real(size(ar1)) )

     if (present(avg1)) avg1=avgarg

     end function r_sigma

     function  r4_sigma(ar1, avg1, iprt) result(std)

     use m_avg
     implicit none

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sigma 
!
! !DESCRIPTION:  finds std dev of a vector (not a matrix)
!
! !CALLING SEQUENCE:
!
!        std = sigma(vectorin)
!
! !INPUT PARAMETERS:
     real (kind=4), dimension(:)  :: ar1
!                       ar1 : input vector to find std deviation
     real (kind=4),    optional   :: avg1
     integer, optional   :: iprt
!
! !OUTPUT PARAMETERS:
     real  (kind=4)               :: std
!                       std      :  std deviation
!
! !SEE ALSO:  m_avg
!
! !REVISION HISTORY:
!
!  21Oct97   M. Karki     Original code.
!EOP
!-------------------------------------------------------------------------

!! LOCAL PARAMETERS: none

! note this is the correct standard deviation when estimating
! the mean
!    std = sqrt(sum( (ar1-avg(ar1))**2) /float(size(ar1)-1) )
! this is the function evaluated in the IDL sigma function
    real (KIND=8) :: avgarg

     avgarg=avg(ar1)

!     std = sqrt(sum( (ar1-avgarg)**2) /float(size(ar1)) )
     std = sqrt(sum( (ar1-avgarg)**2) /real(size(ar1)) )

     if (present(avg1)) avg1=avgarg

     end function r4_sigma

end module m_sigma
