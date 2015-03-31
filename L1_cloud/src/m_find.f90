module m_find

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
    logical, dimension(:),         intent(in)             :: mask
    !                        mask : mask vector
    integer,                 intent(in)           :: count
    !                        count: size of output vector
    !
    ! !OUTPUT PARAMETERS:  
    integer, dimension(count)                :: index
    !                        index : vector of indices where mask is true
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


  function find1(mask,iprt) result (index)

    use tell_module

    implicit none

    logical, dimension(:),  intent(in)      :: mask
    integer, optional,      intent(in)      :: iprt
    integer                                 :: index, iprt1

    integer, dimension(1) :: temp
    integer :: n, errstat
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

      if (allocated(temp2))  deallocate(temp2, stat=errstat)
      allocate(temp2(n), stat=errstat)

      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "find1: allocation failed", &
             errstat)
        call exit(-1)
      endif

      temp2 = find2(mask,n)
      index=temp2(1)
    endif

  end function find1

end module m_find
