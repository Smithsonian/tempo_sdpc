module m_bilinear

interface bilinear
  module procedure bilineara
  module procedure bilinear1
end interface 

contains

  FUNCTION BILINEARa (P,IXin,JYin, x, y ) result (bilinear_res)

    use m_interpol
    use m_findgen
    use tell_module
    implicit none

    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  bilinear
    ! 
    ! !DESCRIPTION: bilinear interpolation, similar to IDL routine
    !
    ! !CALLING SEQUENCE: 
    !
    !        bilinear_res = bilinear(p, ixin, iyin, x, y)
    !     
    ! !INPUT PARAMETERS:   
    real (KIND=8), dimension(:),    intent(in)           :: ixin, jyin
    real (KIND=8), dimension(:),    intent(in), optional :: x, y
    real (KIND=8), dimension(:,:),  intent(in)           :: p
    !
    ! !OUTPUT PARAMETERS:  
    real (KIND=8), dimension(size(ixin))           ::  bilinear_res
    !   bilinear_res : interpolated point(s)
    !
    ! !SEE ALSO:  IDL documentation
    !
    ! !REVISION HISTORY: 
    !
    !  13Aug96   Joiner       Original code converted from IDL code
    !  31Jul14   O'Sullivan   Tidied up for TEMPO. Comments below suggest code 
    !                         should allow array P to have indexed from 0, but
    !                         J.Joiner insists it should be indexed from 1.
    !                         Code added to allow index=0, but commented out
    !                         for now.
    !
    ! !COMMENTS FROM IDL VERSION
    !     IX must satisfy the expression,
    !      0 <= MIN(IX) < N0  and 0 < MAX(IX) <= N0
    !     where N0 is the total number of subscripts in the first dimension
    !     of P.
    !
    !     JY must satisfy the expression,
    !      0 <= MIN(JY) < M0  and 0 < MAX(JY) <= M0
    !     where M0 is the total number of subscripts in the second dimension
    !     of P.
    !
    ! !SIDE EFFECTS:
    !     Note: If x and y are not specified, the grid coordinates are
    !     taken to be 1, 2, 3, .... 
    !
    !EOP
    !---------------------------------------------------------------------

    ! local variables
    real (KIND=8),     dimension(:),  allocatable :: ix, jy
    real (KIND=8),     dimension(:),  allocatable :: dx, dy, dx1, dy1
    integer,  dimension(:),  allocatable :: i, j 
    integer,  dimension(:),  allocatable :: ip, jp
    integer                              :: icnt, sizex
    integer                              :: errstat
    !if p indexed from zero
    !        real (KIND=8), dimension(:,:),allocatable :: p0

    !===============================================================

    sizex=size(ixin)

    allocate(ix(sizex), stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "bilineara: failured to allocate ix", &
           errstat)
      stop 1
    endif

    allocate(jy(size(ix)), &
         i(size(ix)), &
         j(size(ix)), &
         ip(size(ix)), &
         jp(size(ix)), &
         dx(size(ix)), &
         dy(size(ix)), &
         dx1(size(ix)), &
         dy1(size(ix)), &
         stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "bilineara: allocation failure", &
           errstat)
      stop 1
    endif


    !if p indexed from zero
    !   allocate(p0(0:size(p,1)-1,0:size(p,2)-1))
    !   p0=p

    if (present(x)) then
      ix=interpol(findgen(size(x))+1, x, ixin)
    else
      !print *,size(ix), size(ixin)
      !print *, ixin
      ix=ixin
    endif
    if (present(y)) then
      jy=interpol(findgen(size(y))+1, y, jyin)
    else
      jy=jyin
    endif
    I=int(IX)
    J=int(JY)
    IP=I+1  
    JP=J+1
    where (i .ge. size(p,1)) ip=ip-1
    where (j .ge. size(p,2)) jp=jp-1
    !if p indexed from zero, swap in following 4 lines
    !   where (i .ge. size(p0,1)-1) ip=ip-1
    !   where (j .ge. size(p0,2)-1) jp=jp-1
    !   DX=IX-FLOAT(I) 
    !   DY=JY-FLOAT(J)
    DX=IX-real(I) 
    DY=JY-real(J)
    DX1=(1.-DX) 
    DY1=(1.-DY) 
    do icnt=1, size(ix)
      bilinear_res(icnt) = ( P(I(icnt),J(icnt))*DX1(icnt)*DY1(icnt) & 
           + P(I(icnt),JP(icnt))*DX1(icnt)*DY(icnt)  &
           + P(IP(icnt),J(icnt))*DX(icnt)*DY1(icnt)  &
           + P(IP(icnt),JP(icnt))*DX(icnt)*DY(icnt))
      !if p indexed from zero
      !    bilinear_res(icnt) = ( P0(I(icnt),J(icnt))*DX1(icnt)*DY1(icnt) & 
      !        + P0(I(icnt),JP(icnt))*DX1(icnt)*DY(icnt)  &
      !        + P0(IP(icnt),J(icnt))*DX(icnt)*DY1(icnt)  &
      !        + P0(IP(icnt),JP(icnt))*DX(icnt)*DY(icnt))
    enddo

    deallocate(i, j, ip, jp, dx, dy, dx1, dy1, ix, jy, stat=errstat)
    !if p indexed from zero
    !   deallocate(p0)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "bilineara: deallocation failure", &
           errstat)
      stop 1
    endif



  END  function bilineara

FUNCTION BILINEAR1 (P,IXin,JYin, x, y ) result (bilinear_res)
! !INPUT PARAMETERS:   
        real (KIND=8),                  intent(in)           :: ixin, jyin
        real (KIND=8), dimension(:),    intent(in), optional :: x, y
        real (KIND=8), dimension(:,:),  intent(in)           :: p
!
! !OUTPUT PARAMETERS:  
       real (KIND=8) ::  bilinear_res

       real (KIND=8), dimension(1) :: bilin_interp, ixin_arr, jyin_arr

       ixin_arr=ixin
       jyin_arr=jyin
       if (present(x) .and. present(y)) then
         bilin_interp=bilinear(p,ixin_arr,jyin_arr,x=x,y=y)
       else
         bilin_interp=bilinear(p,ixin_arr,jyin_arr)
       endif
       bilinear_res=bilin_interp(1)

END  function bilinear1 

end module m_bilinear
