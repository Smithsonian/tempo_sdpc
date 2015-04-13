!>Linear interpolation routines
!
!-------------------------------------------------------------------------
! 
! !DESCRIPTION: 
!> @brief Linear interpolation, similar to IDL routine
!>               extrapolate outside input range
!
! !INPUT PARAMETERS:   
!> @param v[in] input ordinate vector
!> @param x[in] input abcissa vector
!> @param u[in] abcissa values for output vector r
! !OUTPUT PARAMETERS:  
!> @param r[out] output ordinate vector
!
!> @author  02Mar96   Joiner     Original code.
!> @author  13Aug97   Joiner     Update to fortran 90
!
!-------------------------------------------------------------------------
module m_interpol 

  interface interpol
    module procedure interpol1
    module procedure interpol2
    module procedure interpold
  end interface

contains

  function interpol1(v, x, u) result (r)

    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  interpol
    ! 
    ! !DESCRIPTION: linear interpolation, similar to IDL routine
    !               extrapolates outside input range
    !
    ! !CALLING SEQUENCE: 
    !
    !        r = interpol(v, x, u)
    !     
    ! !INPUT PARAMETERS:   
    real (kind = 4), dimension(:), intent(in) :: v
    real (kind = 4), dimension(:), intent(in) :: x
    real (kind = 4), dimension(:), intent(in) :: u
    !
    ! !OUTPUT PARAMETERS:  
    real (kind = 4), dimension(lbound(u,1):ubound(u,1)) :: r
    !
    ! !SEE ALSO:  IDL documentation
    !
    ! !REVISION HISTORY: 
    !
    !  02Mar96   Joiner     Original code.
    !  13Aug97   Joiner     Update to fortran 90
    !
    !EOP
    !-------------------------------------------------------------------------

    !Local variables
    !---------------
    real (KIND=8)           d
    real (KIND=8)           s1
    integer        i, ix, m2
    integer        n, m

    n = size(u,1)
    m = size(v,1)

    m2=m-1
    r(:)=1 !v(lbound(v,1))

    if (x(2) - x(1) >= 0) then
      s1 = 1
    else
      s1=-1
    endif

    ix=1 !lbound(v,1)
    do i= 1, ubound(u,1)
      d = s1 * (u(i)-x(ix))
      if (d == 0.d0) then
        r(i)=v(ix)
      else
        if (d > 0) then
          do while ((s1*(u(i)-x(ix+1))) > 0 .and. (ix < m2))
            ix=ix+1
          enddo
        else
          do while ((s1*(u(i)-x(ix))) < 0 .and. (ix > 1))
            ix=ix-1
          enddo
        endif
        r(i) = v(ix) + (u(i)-x(ix))*(v(ix+1)-v(ix))/(x(ix+1)-x(ix))
      endif
    enddo

  end function interpol1
  !***********************************************************************



  function interpol2(v, x, u) result (r)

    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  interpol2 
    ! 
    ! !DESCRIPTION: linear interpolation, similar to IDL routine
    !               extrapolates outside input range, calls 
    !          interpol, output is scalar rather than a vector
    !
    ! !CALLING SEQUENCE: 
    !
    !        r = interpol(v, x, u)
    !     
    ! !INPUT PARAMETERS:   
    real (kind = 8), dimension(:), intent(in) :: v
    real (kind = 8), dimension(:), intent(in) :: x
    real (kind = 8),               intent(in) :: u
    !
    ! !OUTPUT PARAMETERS:  
    real (kind = 8)                    :: r
    !
    ! !SEE ALSO:  IDL documentation, interpol, idlmod.f90 (interface)
    !
    ! !REVISION HISTORY: 
    !
    !  02Mar96   Joiner     Original code.
    !  13Aug97   Joiner     Update to fortran 90
    !
    !EOP
    !-------------------------------------------------------------------------

    real (kind = 8), dimension(1)          :: dumu, dumr

    dumu(1) = u
    dumr = interpold(v, x, dumu)
    r = dumr(1)
  end function interpol2



  function interpold(v, x, u) result (r)

    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  interpol
    ! 
    ! !DESCRIPTION: linear interpolation, similar to IDL routine
    !               extrapolates outside input range
    !
    ! !CALLING SEQUENCE: 
    !
    !        r = interpol(v, x, u)
    !     
    ! !INPUT PARAMETERS:   
    real (kind = 8), dimension(:), intent(in) :: v
    real (kind = 8), dimension(:), intent(in) :: x
    real (kind = 8), dimension(:), intent(in) :: u
    !
    ! !OUTPUT PARAMETERS:  
    real (kind = 8), dimension(lbound(u,1):ubound(u,1)) :: r
    !
    ! !SEE ALSO:  IDL documentation
    !
    ! !REVISION HISTORY: 
    !
    !  02Mar96   Joiner     Original code.
    !  13Aug97   Joiner     Update to fortran 90
    !
    !EOP
    !-------------------------------------------------------------------------

    !Local variables
    !---------------
    real (KIND=8)           d
    real (KIND=8)           s1
    integer        i, ix, m2
    integer        n, m

    n = size(u,1)
    m = size(v,1)

    m2=m-1
    r(:)=1 !v(lbound(v,1))

    if (x(2) - x(1) >= 0) then
      s1 = 1
    else
      s1=-1
    endif

    ix=1 !lbound(v,1)
    do i= 1, ubound(u,1)
      d = s1 * (u(i)-x(ix))
      if (d == 0.d0) then
        r(i)=v(ix)
      else
        if (d > 0) then
          do while ((s1*(u(i)-x(ix+1))) > 0 .and. (ix < m2))
            ix=ix+1
          enddo
        else
          do while ((s1*(u(i)-x(ix))) < 0 .and. (ix > 1))
            ix=ix-1
          enddo
        endif
        r(i) = v(ix) + (u(i)-x(ix))*(v(ix+1)-v(ix))/(x(ix+1)-x(ix))
      endif
    enddo

  end function interpold


end module m_interpol
