module m_str_replace

public str_replace

contains

function str_replace (instring, searchstr, newstr) result (repl_str)

   use m_strpos
   implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  str_replace
! 
! !DESCRIPTION: replaces a substring within a string to another string
!
! !CALLING SEQUENCE: 
!
!   Result = str_replace (instring, searchstr, newstr)
!     
! !INPUT PARAMETERS:   
character(len=*), intent(in) :: instring
character(len=*), intent(in) :: searchstr
character(len=*), intent(in) :: newstr
!
! !OUTPUT PARAMETERS:  
character(len=len(instring)-len(searchstr)+len(newstr)) :: repl_str
!
! !SEE ALSO:  
!
! !REVISION HISTORY: 
!
!  13Dec99   Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------
integer :: ipos

ipos = strpos(instring,searchstr)
repl_str = 'x'
if (ipos > 0 .and. ipos <= len(instring)) then
  repl_str(1:ipos-1)=instring(1:ipos-1)
  repl_str(ipos:ipos+len(newstr)-1)=newstr
  repl_str(ipos+len(newstr):len(instring)-len(searchstr)+len(newstr)) &
    = instring(ipos+len(searchstr):len(instring))
endif

end function str_replace

end module m_str_replace
