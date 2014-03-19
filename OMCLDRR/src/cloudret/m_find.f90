module m_find

!interface find
!  module procedure find2
!  module procedure find1
!end interface 

public

contains

function find2(mask,count) result (index)
implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  find
! 
! !DESCRIPTION: similar to IDL's "where" function, but needs to know
!               the dimension of the output (can get from "count")
!
! !CALLING SEQUENCE: 
!
!        index = find(mask, count)
!     
! !INPUT PARAMETERS:   
logical, dimension(:), 	intent(in)     	:: mask
!			mask : mask vector
integer, 		intent(in)   	:: count
!			count: size of output vector
!
! !OUTPUT PARAMETERS:  
integer, dimension(count)		:: index
!			index : vector of indices where mask is true
!
! !SEE ALSO:  IDL documentation
!
! !REVISION HISTORY: 
!
!  13Aug97   Joiner     original code
!
!EOP
!-------------------------------------------------------------------------

!local variables
!---------------
integer					:: i
integer					:: j

!----------------------------------------------
j = 1 
if (count > 0) then
do i=1, size(mask,1) 
  if (mask(i)) then
    index(j) = i
    j = j + 1
  endif
enddo
j = j - 1
endif
if (j /= count .or. count <= 0) then
  print *, 'ERROR!!! in find (dimensions not correct)', j, count
endif

end function find2

function findp(mask,counts) result(index)
implicit none
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  findp
! 
! !DESCRIPTION: similar to IDL's "where" function, but needs to know
!               the dimension of the output (can get from "count")
!
! !CALLING SEQUENCE: 
!
!        index = find(mask, counts)
!     
! !INPUT PARAMETERS:   
logical, dimension(:), 	intent(in)     	:: mask
!			mask : mask vector
 integer, 		intent(out), optional :: counts
!			counts: size of output vector
!
! !OUTPUT PARAMETERS:  
integer, dimension(:), pointer		:: index
!			index : vector of indices where mask is true
!
! !SEE ALSO:  IDL documentation
!
! !REVISION HISTORY: 
!
!  13Aug97   Joiner     original code
!
!EOP
!-------------------------------------------------------------------------

!local variables
!---------------
integer					:: i
integer					:: j
integer                                 :: tempsize
integer                                 :: error
integer                                 :: ncounts

!----------------------------------------------
j = 1 
ncounts = count(mask)
deallocate(index)
if (ncounts > 0) then
 allocate(index(ncounts))
 do i=1, size(mask,1) 
  if (mask(i)) then
    index(j) = i
    j = j + 1
  endif
 enddo
else
  allocate(index(1))
  index=-1
endif
j = j - 1
if (present(counts)) counts=ncounts

end function findp


function find1(mask,iprt) result (index)

implicit none

logical, dimension(:),  intent(in)      :: mask
integer, optional,      intent(in)      :: iprt
integer                                 :: index, iprt1

integer, dimension(1) :: temp
integer :: n
integer, dimension(:), allocatable :: temp2

iprt1=1
if (present(iprt)) iprt1=iprt
n = count(mask)
if (n == 1) then
  temp(1:1)=find2(mask,1)
  index=temp(1)
else
  if (iprt1 >= 1) then
    print *,'find: WARNING, found more than 1 value'
    print *,'taking the first'
  endif 
  allocate(temp2(n))
  temp2 = find2(mask,n)
  index=temp2(1)
  deallocate(temp2)
endif

end function find1

end module m_find
