module m_interp_pres

contains

!subroutine interp_pres(ix1, ix2, rad, rads, pres, pres_int)
subroutine interp_pres(ix1, ix2, rad, pres, pres_int, jacob)

use m_cloud_pres_mod, ONLY: temp2D
   implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  interp_pres
! 
! !DESCRIPTION: interp_pres allocates/deallocates memory for 
!               retrievals		
!
! !CALLING SEQUENCE: 
!
!        call interp_pres
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
integer, intent(in) :: ix1, ix2
real,    intent(in) :: pres_int
real, dimension(:), intent(out) :: rad, jacob
real, dimension(:), intent(in) :: pres
!real, dimension(:,:), intent(in) :: rads
real :: temp1, temp2
!real, dimension(:,:), pointer :: rads

character(len=50) :: myname='interp_pres: '
integer :: i

!**************************************************************************

  !temp=(pres(ix2)-pres_int)/(pres(ix2)-pres(ix1))
  temp1=(pres_int-pres(ix1))
  temp2=(pres(ix2)-pres(ix1))
  do i=1,size(rad)
   jacob(i)=(temp2D(2,i)-temp2D(1,i))/temp2 
  enddo
  if (temp1 /= 0) then
   do i=1,size(rad)
     rad(i)=temp2D(1,i) + (jacob(i))*temp1 
    ! rad(i)=temp2D(2,i) - (temp2D(2,i)-temp2D(1,i))*temp 
    !rad(i)=rads(ix2,i) - (rads(ix2,i)-rads(ix1,i))*temp &
    !print *, i, rads(ix2,i), rads(ix1,i), (rads(ix2,i)-rads(ix1,i)), &
    !    (pres(ix2)-pres(ix1))*(pres(ix2)-pres_int), rad(i)
   enddo ! nobs
  else
   rad=temp2D(1,:)
  endif

end subroutine interp_pres

subroutine interp_rads(ix1, ix2, pres, pres_int,  i0_1, i0_2, sb_1, sb_2, &
     tr_1, tr_2, i0, sb, tr )

   implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  interp_pres
! 
! !DESCRIPTION: interp_pres allocates/deallocates memory for 
!               retrievals		
!
! !CALLING SEQUENCE: 
!
!        call interp_pres
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
integer, intent(in)  :: ix1, ix2
real,    intent(in)  :: pres_int, i0_1, i0_2, sb_1, sb_2, tr_1, tr_2
real,    intent(out) :: i0, sb, tr 
real, dimension(:), intent(in) :: pres

character(len=50) :: myname='interp_rads: '
real :: temp
!integer :: i

!**************************************************************************

  temp=(pres(ix2)-pres_int)/(pres(ix2)-pres(ix1))
!  do i=1,size(rad)
    i0=i0_2 - (i0_2-i0_1)*temp 
    sb=sb_2 - (sb_2-sb_1)*temp 
    tr=tr_2 - (tr_2-tr_1)*temp 
    !print *, i, &
!  enddo ! nobs

end subroutine interp_rads

end module m_interp_pres
