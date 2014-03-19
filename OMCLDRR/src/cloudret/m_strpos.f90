module m_strpos

public strpos

contains

function strpos (string, searchstr) result (position)
   implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  strpos
! 
! !DESCRIPTION: similar to IDL function
!
! !CALLING SEQUENCE: 
!
!   Result = STRUPCASE(String)
!     
! !INPUT PARAMETERS:   
character(len=*), intent(in)		:: string
character(len=*), intent(in)		:: searchstr
!
! !OUTPUT PARAMETERS:  
integer :: position
!
! !SEE ALSO:  IDL documentation
!
! !REVISION HISTORY: 
!
!  13Dec99   Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------
integer            :: i, j, temp, tot
logical            :: goodsearch, done

  position = -1
   i=1
   tot=1
   done=.false.
   goodsearch=.false.
   do while (i <= len(string) .and. .not. done)
    ! find the first occurrence
    !==========================
    goodsearch=.false.
    temp = scan(string(i:len(string)),searchstr(1:1))
    if (temp /= 0) then
      goodsearch=.true.
      j=2
      do while (j <= len(searchstr) .and. goodsearch)
        goodsearch=string(i+temp+j-2:i+temp+j-2) == searchstr(j:j)
        j=j+1
      enddo
      done = goodsearch .or. i == len(string)! done if found a match
      if (.not. done) i=temp+i ! set to next substring to search
    else
      done=.true.
    endif ! first character found
   enddo
   if (goodsearch) position=temp+i-1

end function strpos

end module m_strpos
