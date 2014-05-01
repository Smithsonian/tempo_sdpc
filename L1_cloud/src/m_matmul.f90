module m_matmul

public 

interface operator (.mm.)
  module procedure matmul1 
  module procedure matmul2 
  module procedure matmul3 
  module procedure matmul4 
end interface

contains

   function matmul1 (a,b) result (c)
implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  matmul1
! 
! !DESCRIPTION: calls fortran matmul, designed to work with .mm. operator
!                c = a .mm. b (from c = a \# b in IDL)
!
! !CALLING SEQUENCE: 
!
!        c = matmul1(a,b)
!     
! !INPUT PARAMETERS:   
     real (KIND=8), dimension(:,:), intent(in) :: a
!                        a : 1st matrix operand
     real (KIND=8), dimension(:,:), intent(in) :: b
!                        b : 2nd matrix operand
!
! !OUTPUT PARAMETERS:  
     real (KIND=8), dimension(size(a,1),size(b,2))   :: c
!                        c : result matrix
!
! !SEE ALSO:  IDL documentation
!
! !REVISION HISTORY: 
!
!  13Aug97   Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------

     if (size(a,2) /= size(b,1) ) then
       print *,'incorrect matrix dimensions in .mm.'
       return
     endif
     c = matmul(a, b)

   end function matmul1

   function matmul2 (a,b) result (c)
implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  matmul2
! 
! !DESCRIPTION: calls fortran matmul, designed to work with .mm. operator
!                c = a .mm. b (from c = a \# b in IDL), same as matmul1,
!                but with matrix and vector
!
! !CALLING SEQUENCE: 
!
!        c = matmul1(a,b)
!     
! !INPUT PARAMETERS:   
     real (KIND=8), dimension(:,:), intent(in)         :: a
!                        a : 1st matrix operand
     real (KIND=8), dimension(:),   intent(in)         :: b
!                        b : 2nd vector operand
!
! !OUTPUT PARAMETERS:  
     real (KIND=8), dimension(size(a,1))           :: c
!                        c : result matrix
!
! !SEE ALSO:  IDL documentation
!
! !REVISION HISTORY: 
!
!  13Aug97   Joiner     original fortran 90
!
!EOP
!----------------------------------------------


     if (size(a,2) /= size(b,1) ) then
       print *,'incorrect matrix dimensions in .mm.'
       return
     endif
     c = matmul(a, b)

   end function matmul2

   function matmul3 (a,b) result (c)
   implicit none

     real (KIND=8), dimension(:),   intent(in)         :: a
     real (KIND=8), dimension(:,:), intent(in)         :: b
     real (KIND=8), dimension(size(b,2))           :: c

     if (size(a,1) /= size(b,1) ) then
       print *,'incorrect matrix dimensions in .mm.'
       return
     endif
     c = matmul(a, b)

   end function matmul3

   function matmul4 (a,b) result (c)
    implicit none

     real (KIND=8), dimension(:),   intent(in)         :: a
     real (KIND=8), dimension(:), intent(in)         :: b
     real (KIND=8), dimension(1)                   :: c

     if (size(a,1) /= size(b,1) ) then
       print *,'incorrect matrix dimensions in .mm.'
       return
     endif
     c = dot_product(a, b)

   end function matmul4

end module m_matmul
