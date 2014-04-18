module m_alloc2

contains

subroutine alloc2()

use m_cloud_pres_mod
   implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  alloc2
! 
! !DESCRIPTION: alloc2 allocates/deallocates memory for 
!               retrievals		
!
! !CALLING SEQUENCE: 
!
!        call alloc2
!     
! !INPUT PARAMETERS:   
!
! !OUTPUT PARAMETERS:  
!
! !SEE ALSO:  
!
! !REVISION HISTORY: 
!
!  05Jan01   Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------
!
!character(len=50) :: myname='alloc2: '

!**************************************************************************

if (allocated(x)) then
  !deallocate memory
  !=================
  deallocate(x)
  deallocate(x_fg)
  deallocate(h)
  deallocate(htr)
  deallocate(err_cov)
  deallocate(corr)
  deallocate(b_i)
endif
if(allocated(y_back)) deallocate(y_back)
allocate(x(0:nst-1,1))
allocate(x_fg(0:nst-1,1))
allocate(h(0:nobs-1,0:nst-1))   
allocate(htr(0:nst-1,0:nobs-1))   
allocate(err_cov(0:nst-1,0:nst-1))   
allocate(corr(0:nst-1,0:nst-1))   
allocate(b_i(0:nst-1)) 
allocate(y_back(0:nst-1,1)) 

end subroutine alloc2

end module m_alloc2
