module m_poly_fit   

private
public poly_fit

contains
! $Id$
!+
! NAME:
!        POLY_FIT
!
! PURPOSE:
!        Perform a least-square polynomial fit with optional error estimates.
!
!        This routine uses matrix inversion.  A newer version of this routine,
!        SVDFIT, uses Singular Value Decomposition.  The SVD technique is more
!        flexible, but slower.
!
! CATEGORY:
!        Curve fitting.
!
! CALLING SEQUENCE:
!   Result = POLY_FIT(X, Y, Degree)
!
! INPUTS:
!   X:  The independent variable vector.
!
!   Y:  The dependent variable vector, should be same length as x.
!
!   Degree: The degree of the polynomial to fit.
!
! OUTPUTS:
!        POLY_FIT returns a vector of coefficients with a length of NDegree+1.
!
! KEYWORDS:
!   CHISQ:   Sum of squared errors divided by MEASURE_ERRORS if specified.
!
!   COVAR:   Covariance matrix of the coefficients.
!
!   MEASURE_ERRORS: Set this keyword to a vector containing standard
!       measurement errors for each point Y[i].  This vector must be the same
!       length as X and Y.
!
!     Note - For Gaussian errors (e.g. instrumental uncertainties),
!        MEASURE_ERRORS should be set to the standard
!              deviations of each point in Y. For Poisson or statistical weighting
!              MEASURE_ERRORS should be set to sqrt(Y).
!
!   SIGMA:   The 1-sigma error estimates of the returned parameters.
!
!     Note: if MEASURE_ERRORS is omitted, then you are assuming that
!           your model is correct. In this case,
!           SIGMA is multiplied by SQRT(CHISQ/(N-M)), where N is the
!           number of points in X. See section 15.2 of Numerical Recipes
!           in C (Second Edition) for details.
!
!   STATUS = Set this keyword to a named variable to receive the status
!          of the operation. Possible status values are:
!          0 for successful completion, 1 for a singular array (which
!          indicates that the inversion is invalid), and 2 which is a
!          warning that a small pivot element was used and that significant
!          accuracy was probably lost.
!
!    Note: if STATUS is not specified then any error messages will be output
!          to the screen.
!
!   YBAND:        1 standard deviation error estimate for each point.
!
!   YERROR: The standard error between YFIT and Y.
!
!   YFIT:   Vector of calculated Y's. These values have an error
!           of + or - YBAND.
!
! COMMON BLOCKS:
!        None.
!
! SIDE EFFECTS:
!        None.
!
! MODIFICATION HISTORY:
!        Written by: George Lawrence, LASP, University of Colorado,
!                December, 1981.
!
!        Adapted to VAX IDL by: David Stern, Jan, 1982.
!       Modified:    GGS, RSI, March 1996
!                    Corrected a condition which explicitly converted all
!                    internal variables to single-precision float.
!                    Added a check for singular array inversion.
!                     SVP, RSI, June 1996
!                     Changed A to Corrm to match IDL5.0 docs.
!                    S. Lett, RSI, December 1997
!                     Changed inversion status check to check only for
!                     numerically singular matrix.
!                    S. Lett, RSI, March 1998
!                     Initialize local copy of the independent variable
!                     to be of type DOUBLE when working in double precision.
!       CT, RSI, March 2000: Changed to call POLYFITW.
!       CT, RSI, July-Aug 2000: Removed call to POLYFITW,
!                   added MEASURE_ERRORS keyword,
!                   added all other keywords (except DOUBLE),
!                   made output arguments obsolete.
!-
   
function poly_fit( x, y, ndegree, &   
        yfit, measure_errors ) result(res)   
   
use m_invert
use m_matmul
       implicit NONE          ! *** IDL2F9O ***

integer, intent(in) :: ndegree
real (KIND=8), intent(in), dimension(:) :: x,y
real (KIND=8), intent(out), dimension(:), optional :: yfit
real (KIND=8), intent(out), dimension(:), optional :: measure_errors
!real (KIND=8), dimension(:), allocatable :: sdev, sdev2
real (KIND=8), dimension(size(x)) :: sdev, sdev2
!real (KIND=8), dimension(:,:), allocatable :: covar
real (KIND=8), dimension(0:ndegree,0:ndegree) :: covar
!real (KIND=8), dimension(:), allocatable :: b, z, wy
real (KIND=8), dimension(0:ndegree) :: b
real (KIND=8), dimension(size(x)) :: z, wy
real (KIND=8), dimension(ndegree+1) :: res
real (KIND=8) :: sum1
integer :: status
integer :: n,m,p,k,j
logical :: no_weight

        n = size(x)   
        if(n /= size(y)) then 
           print *, 'x and y must have same number of elements.'   
           return
        endif
        m = ndegree + 1           ! # of elements in coeff vec
   
        no_weight =.not. present(measure_errors)
        !allocate(sdev(n)) ; 
        sdev = 1.0   
        if (.not. no_weight) sdev = sdev*measure_errors   
        !allocate(sdev2(n)) 
        sdev2 = sdev**2   
   
           ! construct work arrays
           !allocate(covar(0:m-1,0:m-1)) 
           covar=0 
                    ! least square matrix, weighted matrix
           !allocate(b(0:m-1)) ; 
           b=0  ! will contain sum weights*y*x^j
           !allocate(z(n)) ; 
           z=1. ! basis vector for constant term
           !allocate(wy(n)) 
           wy = y/sdev2   
 !print *, 'y',y
 !print *, 'sdev',sdev
 !print *, 'sdev2',sdev2
 !print *, 'wy',wy
   
  !     covar(1,1) = no_weight ? n : sum(1/sdev2)   
        if (no_weight) then
          covar(0,0)=n
        else
          covar(0,0)=sum(1/sdev2)
        endif
            
        b(0) = sum(wy)   
   
        do p = 1,2*ndegree        ! power loop
                z = z*x           ! z is now x^p
                if (p < m) b(p) = sum(wy*z)   ! b is sum weights*y*x^j
                sum1 = sum(z/sdev2)   
                do j = maxval((/0,p-ndegree/)), minval((/ndegree,p/))
                   covar(j,p-j) = sum1   
                enddo
        enddo    ! end of p loop, construction of covar and b
   
   
        covar = invert(covar,status)   
   
        !endif   
   
        !print *, 'b' ,b
        !print *, 'covar' ,covar
        res = b .mm. covar  ! construct coefficients
   
        ! compute optional output parameters.
   
      if (present(yfit)) then
        yfit = res(1) !res(ndegree+1)   
        do k = 2, ndegree+1
          yfit=yfit+res(k)*x**(k-1)
        enddo
      endif
        !do k = ndegree-1, 0, -1 
        !  yfit = res(k+1) + yfit*x  ! sum basis vectors
        !enddo
   
 ! deallocate(b)
 ! deallocate(z)
 ! deallocate(covar)
 ! deallocate(wy)
   
! ***********************************************************
end function poly_fit

end module m_poly_fit   
