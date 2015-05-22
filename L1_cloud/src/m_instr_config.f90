module m_instr_config

private
public instr_config 

contains

subroutine instr_config(ierr, izoom)

   use m_vars, ONLY: iLine, iprt, config_rad, config_irr, mflg

   implicit none

   integer, intent(out) :: ierr, izoom
   integer :: i 
   integer (kind=1), dimension(4) :: id1=(/0_1,1_1,2_1,7_1/)
!   integer (kind=1), dimension(7) :: id2=(/42,43,44,49,56,57,58/)  

izoom=1
if(iprt>1) then 
  print *,'config_irr',config_irr,' config_rad',config_rad(iLine),'rebinning flag',btest(mflg(iLine),7)
endif

! simplified testing the zoom radiance measurements
    ierr=0
    do i=1,4
      if(config_rad(iLine)==id1(i)) izoom=0
    enddo

! reserved for future use
!ierr=1
!if(config_irr==8) then
!testing the rebinning flag
!  if(btest(mflg(iLine),7)) then
!    do i=1,7
!      if(config_rad(iLine)==id2(i)) ierr=0
!    enddo
!  else
!    do i=1,4
!      if(config_rad(iLine)==id1(i)) ierr=0
!    enddo
!  endif
!endif

!if(config_irr==50) then
!    do i=1,4
!      if(config_rad(iLine)==id2(i)) ierr=0
!    enddo
!endif

!if(config_irr==62) then
!    do i=5,7
!      if(config_rad(iLine)==id2(i)) ierr=0
!    enddo
!endif

!set instrument settings error flag
!if(ierr==1) meas_qual_flg(iLine)=IBSET(meas_qual_flg(iLine),7)

end subroutine instr_config

end module m_instr_config
       
    

