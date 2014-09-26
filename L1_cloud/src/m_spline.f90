module m_spline

  interface spline
    module procedure spline1
    module procedure spline2
  end interface

contains

  function spline2(x,y,t,sigma_in) result(spl)   

    implicit NONE          ! *** IDL2F9O ***

    !+
    ! NAME:
    !        SPLINE
    !
    ! PURPOSE:
    !        This function performs cubic spline interpolation.
    !
    ! CATEGORY:
    !        Interpolation - E1.
    !
    ! CALLING SEQUENCE:
    !        Result = SPLINE(X, Y, T [, Sigma])
    !
    ! INPUTS:
    real (KIND=8), dimension(:), intent(in) :: x, y, t
    real (KIND=8), intent(in), optional :: sigma_in
    !        X:        The abcissa vector. Values MUST be monotonically increasing.
    !
    !        Y:        The vector of ordinate values corresponding to X.
    !
    !        T:        The vector of abcissae values for which the ordinate is 
    !                desired. The values of T MUST be monotonically increasing.
    !
    ! OPTIONAL INPUT PARAMETERS:
    real (KIND=8), dimension(size(t)) :: spl
    !        Sigma:        The amount of "tension" that is applied to the curve. The 
    !                default value is 1.0. If sigma is close to 0, (e.g., .01),
    !                then effectively there is a cubic spline fit. If sigma
    !                is large, (e.g., greater than 10), then the fit will be like
    !                a polynomial interpolation.
    !
    ! OUTPUTS:
    !        SPLINE returns a vector of interpolated ordinates.
    !        Result(i) = value of the function at T(i).
    !
    ! RESTRICTIONS:
    !        Abcissa values must be monotonically increasing.
    !
    ! EXAMPLE:
    !        The commands below show a typical use of SPLINE:
    !
    !                X = [2.,3.,4.]          ;X values of original function
    !                Y = (X-3)^2             ;Make a quadratic
    !                T = FINDGEN(20)/10.+2         ;Values for interpolated points.
    !                                        ;twenty values from 2 to 3.9.
    !                Z = SPLINE(X,Y,T)         ;Do the interpolation.
    !
    !
    !
    ! MODIFICATION HISTORY:
    !        Author:        Walter W. Jones, Naval Research Laboratory, Sept 26, 1976.
    !        Reviewer: Sidney Prahl, Texas Instruments.
    !        Adapted for IDL: DMS, Research Systems, March, 1983.
    !        Converted to f90: Joanna Joiner, March 2001
    !
    !-
    !
    !real (KIND=8), dimension(:), allocatable :: xx, yp
    !integer, dimension(:), allocatable :: subs, subs1
    real (KIND=8), dimension(0:size(x)*2-1) :: yp
    real (KIND=8), dimension(0:size(x)-1) :: xx
    integer, dimension(0:size(t)-1) :: subs, subs1
    real (KIND=8) :: sigma, sigmap, dels, sinhs, sinhin, diag1, diagin
    real (KIND=8) :: spdiag, delx2, dx2, diag2, s, exps
    !real (KIND=8), dimension(:), allocatable :: sinhs2, del1, del2, dels2, exps1, exps2, sinhd1, sinhd2
    real (KIND=8), dimension(0:size(t)-1) :: sinhs2, del1, del2, dels2, &
         exps1, exps2, sinhd1, sinhd2
    real (KIND=8) :: delx1, dx1, c1, c2, c3, slpp1, deln, delnm1, delnn
    real (KIND=8) :: delx12, slppn
    integer :: nm1, np1, i, m, j, n

    !call pzeitbeg('nspline')
    n = size(x) ! < size(y)   
    !print *,size(x), size(t)
    !print *,lbound(xx), ubound(xx)
    !print *,lbound(del1), ubound(del1)
    !

    !to avoid uninitialized variable warnings, set some params to zero
    slpp1=0.
    slppn=0.
    dx2=0.

    sigma=1.
    if (present(sigma_in)) sigma=sigma_in
    !allocate(xx(0:n-1))
    xx = x * 1.                           !Make X values floating if not.
    !allocate(yp(0:n*2-1))                   !temp storage
    delx1 = xx(1)-xx(0)                   !1st incr
    dx1=(y(2)-y(1))/delx1   
    !
    nm1 = n-1
    np1 = n+1
    if (n == 2) then   
      yp(0)=0.   
      yp(1)=0.   
    else   
      delx2 = xx(2)-xx(1)   
      delx12 = xx(2)-xx(0)   
      c1 = -(delx12+delx1)/delx12/delx1   
      c2 = delx12/delx1/delx2   
      c3 = -delx1/delx12/delx2   
      !
      slpp1 = c1*y(1)+c2*y(2)+c3*y(3)   
      deln = xx(nm1)-xx(nm1-1)   
      delnm1 = xx(nm1-1)-xx(nm1-2)   
      delnn = xx(nm1)-xx(nm1-2)   
      c1=(delnn+deln)/delnn/deln   
      c2=-delnn/deln/delnm1   
      c3=deln/delnn/delnm1   
      slppn = c3*y(nm1-1)+c2*y(nm1)+c1*y(nm1+1)   
    endif
    !
    sigmap = sigma*nm1/(xx(nm1)-xx(0))   
    dels = sigmap*delx1   
    exps = exp(dels)   
    sinhs = .5d0*(exps-1./exps)   
    sinhin=1./(delx1*sinhs)   
    diag1 = sinhin*(dels*0.5d0*(exps+1./exps)-sinhs)   
    diagin = 1./diag1   
    yp(0)=diagin*(dx1-slpp1)   
    spdiag = sinhin*(sinhs-dels)   
    yp(n)=diagin*spdiag   
    !
    if  (n > 2) then 
      do i=1,nm1-1    
        delx2 = xx(i+1)-xx(i)   
        dx2=(y(i+2)-y(i+1))/delx2   
        dels = sigmap*delx2   
        exps = exp(dels)   
        sinhs = .5d00 *(exps-1./exps)   
        sinhin=1./(delx2*sinhs)   
        diag2 = sinhin*(dels*(.5*(exps+1./exps))-sinhs)   
        diagin = 1./(diag1+diag2-spdiag*yp(n+i-1))   
        yp(i)=diagin*(dx2-dx1-spdiag*yp(i-1))   
        spdiag=sinhin*(sinhs-dels)   
        yp(i+n)=diagin*spdiag   
        dx1=dx2   
        diag1=diag2   
      enddo
    endif
    !

    diagin=1./(diag1-spdiag*yp(n+nm1-1))   
    yp(nm1)=diagin*(slppn-dx2-spdiag*yp(nm1-1))   
    do i=n-2,0,-1 
      yp(i)=yp(i)-yp(i+n)*yp(i+1)                   
    enddo
    !
    !
    m = size(t)   
    !allocate(subs(0:m-1))
    !allocate(subs1(0:m-1))
    subs=nm1
    s = xx(nm1)-xx(0)   
    sigmap = sigma*nm1/s   
    j=0
    do i=1,nm1 !find subscript where xx(subs) > t(j+1) > xx(subs-1)
      do while (xx(i) > t(j+1))     
        subs(j)=i   
        j=j+1   
        if (j == m) goto 999
      enddo
    enddo

999 continue
    subs1 = subs-1   
    !allocate(del1(0:m-1))
    !allocate(del2(0:m-1))
    !allocate(dels2(0:m-1))
    !allocate(exps1(0:m-1))
    !allocate(exps2(0:m-1))
    !allocate(sinhd1(0:m-1))
    !allocate(sinhd2(0:m-1))
    !allocate(sinhs2(0:m-1))
    del1 = t-xx(subs1)   
    del2 = xx(subs)-t   
    dels2 = xx(subs)-xx(subs1)   
    exps1=exp(sigmap*del1)   
    sinhd1 = .5*(exps1-1./exps1)   
    exps2=exp(sigmap*del2)   
    sinhd2=.5*(exps2-1./exps2)   
    exps2 = exps1*exps2   
    sinhs2=.5*(exps2-1./exps2)   
    spl=(yp(subs)*sinhd1+yp(subs1)*sinhd2)/sinhs2+ &   
         ((y(subs+1)-yp(subs))*del1+(y(subs1+1)-yp(subs1))*del2)/dels2   
    !deallocate(del1)
    !deallocate(del2)
    !deallocate(dels2)
    !deallocate(exps1)
    !deallocate(exps2)
    !deallocate(sinhd1)
    !deallocate(sinhd2)
    !deallocate(sinhs2)
    !deallocate(subs)
    !deallocate(subs1)
    !deallocate(xx)
    !deallocate(yp)
    !       if (m == 1) then return,spl(0) else return,spl   
    !call pzeitend
    ! ***********************************************************
  end  function spline2

  function spline1(x,y,t,sigma_in) result(spl)
    real (KIND=8), dimension(:), intent(in) :: x, y
    real (KIND=8), intent(in) :: t
    real (KIND=8), intent(in), optional :: sigma_in
    real (KIND=8) :: spl
    real (KIND=8), dimension(1)          :: dumu, dumr

    dumu(1) = t
    if (present(sigma_in)) then
      dumr = spline2(x, y, dumu, sigma_in=sigma_in)
    else 
      dumr = spline2(x, y, dumu)
    endif
    spl = dumr(1)
  end function spline1

end module m_spline
