!>Return vector of indices for which input expression is true/non-zero
!
!-------------------------------------------------------------------------
!
! !ROUTINE:  find
! 
! !DESCRIPTION: 
!> similar to IDL's "where" function, but needs to know
!> the dimension of the output (can get from "count")
!
! !CALLING SEQUENCE: 
!
!        index = find(mask, count)
!     
! !INPUT PARAMETERS:   
!> @param mask[in]  mask vector
!> @param count[in] size of output vector
!
! !OUTPUT PARAMETERS:  
!> @param index[out] vector of indices where mask is true
!
! !SEE ALSO:  IDL documentation
!
! !REVISION HISTORY: 
!
!> @author  13Aug97   Joiner     original code
!
!-------------------------------------------------------------------------
module m_find

  public

contains

  function find2(mask,count) result (index)

    implicit none
    ! !INPUT PARAMETERS:   
    logical, dimension(:),         intent(in)             :: mask
    integer,                 intent(in)           :: count
    ! !OUTPUT PARAMETERS:  
    integer, dimension(count)                :: index
    !local variables
    !---------------
    integer                                        :: i
    integer                                        :: j

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


  function find1(mask) result (index)

    use tell_module

    implicit none

    logical, dimension(:),  intent(in)      :: mask
    integer                                 :: index

    integer, dimension(1) :: temp
    integer :: n, errstat
    integer, dimension(:), allocatable :: temp2

    n = count(mask)
    if (n == 1) then
      temp(1:1)=find2(mask,1)
      index=temp(1)
    else
      call tell_log(1,'find: WARNING, found more than 1 value')
      call tell_log(1,'find: taking the first')

      if (allocated(temp2))  deallocate(temp2, stat=errstat)
      allocate(temp2(n), stat=errstat)

      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "find1: allocation failed", &
             errstat)
        stop 1
      endif

      temp2 = find2(mask,n)
      index=temp2(1)
    endif

  end function find1

end module m_find
