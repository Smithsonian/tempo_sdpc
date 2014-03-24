module m_num2string

private
public num2string

 interface num2string
   module procedure i_string
   module procedure i_string_s
   module procedure r_string
   module procedure r_string_s
 end interface 

contains

function i_string (number, formatstr) result (line)
   implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  i_string, etc
! 
! !DESCRIPTION: similar to IDL "string" function, designed to work
!		with string interface
!
! !CALLING SEQUENCE: 
!
!        line = i_string(number, formatstr)
!     
! !INPUT PARAMETERS:   
integer, dimension(:), intent(in) 		:: number
!			number    : input number (or array)
character(len=*), optional, intent(in) 		:: formatstr
!			formatstr : optional format string
!
! !OUTPUT PARAMETERS:  
character(len=72) 				:: line
!			line      : formatted character string
!
! !SEE ALSO:  IDL documentation
!
! !REVISION HISTORY: 
!
!  13Aug97   Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------

if ( present(formatstr) ) then
  write(line, fmt=formatstr )  number
else
  write(line, *) number
  line=trim(adjustl(line))
endif

end function i_string

function i_string_s (number, formatstr) result (line)
   implicit none

integer, intent(in) 				:: number
character(len=*), optional, intent(in) 		:: formatstr
character(len=72) 				:: line

if ( present(formatstr) ) then
  write(line,formatstr) number
else
  write(line,*) number
  line=trim(adjustl(line))
endif

end function i_string_s


function r_string (number, formatstr) result (line)
   implicit none

real, dimension(:), intent(in) 			:: number
character(len=*), optional, intent(in)		:: formatstr
character(len=72) 				:: line

if ( present(formatstr) ) then
  write(line,formatstr) number
else
  write(line,*) number
  line=trim(adjustl(line))
endif

end function r_string

function r_string_s (number, formatstr) result (line)
   implicit none

real, intent(in)	 			:: number
character(len=*), optional, intent(in)		:: formatstr
character(len=72) 				:: line

if ( present(formatstr) ) then
  write(line,formatstr) number
else
  write(line,*) number
  line=trim(adjustl(line))
endif

end function r_string_s

end module m_num2string
