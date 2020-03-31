MODULE chebyshev_module
  
  USE parameters_module
  USE error_module,     ONLY : ErrorType, RaisePixelError, CheckError

  !********************************************************
  !*    Collection of Chebyshev approximation routines    *
  !* ---------------------------------------------------- *
  !* REFERENCE: "Numerical Recipes, The Art of Scientific *
  !*             Computing By W.H. Press, B.P. Flannery,  *
  !*             S.A. Teukolsky and W.T. Vetterling,      *
  !*             Cambridge University Press, 1986"        *
  !*             [BIBLI 08].                              *
  !*                                                      *
  !*                F90 Release 1.1 By J-P Moreau, Paris. *
  !*                        (www.jpmoreau.fr)             *
  !* ---------------------------------------------------- *
  !* Release 1.1: added subroutines CHINT and CHDER.      *
  !********************************************************
  
  IMPLICIT NONE
  
  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'chebyshev_module'
  REAL(KIND=8),     PARAMETER :: ZERO = 0.0d0
  REAL(KIND=8),     PARAMETER :: ONE  = 1.0d0
  REAL(KIND=8),     PARAMETER :: TWO  = 2.0d0
  REAL(KIND=8),     PARAMETER :: HALF = 0.5d0
  REAL(KIND=8),     PARAMETER :: QUART= 0.25d0
  INTEGER,          PARAMETER :: NMAX = 50
  
  PRIVATE :: ModuleName, ZERO, ONE, TWO, HALF, NMAX, QUART

  CONTAINS
  
  SUBROUTINE ConstructChebyshevBasis(nx,xmin,xmax,xgrid,nT,T,Error)

      ! --------------------
      ! subroutine arguments
      ! --------------------
      INTEGER,        INTENT(IN)    :: nx
      REAL(KIND=8),   INTENT(IN)    :: xmin, xmax
      REAL(KIND=8),   INTENT(IN)    :: xgrid(nx)
      INTEGER,        INTENT(IN)    :: nT
      REAL(KIND=8),   INTENT(OUT)   :: T(nx,nT)
      TYPE(ErrorType),INTENT(INOUT) :: Error

      ! ---------------
      ! local variables
      ! ---------------
      INTEGER      :: i, b
      REAL(KIND=8) :: x_new

      CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ConstructChebyshevBasis'

      ! ============================================================
      ! ConstructChebyshevBasis starts here
      ! ============================================================

      ! Error checking goes here
      IF(nT .LE. 0) STOP 'Chebyshev basis must be > 0'
      IF(MINVAL(xgrid) .LT. xmin-TINY(0.0d0) ) THEN
        CALL RaisePixelError(Error, ErrorCode_SpecFit,            &
                             ModuleName, SubroutineName,          &
                             'xgrid is outside chebyshev interval')
      ENDIF
      IF(MAXVAL(xgrid) .GT. xmax+TINY(0.0d0) ) THEN
        CALL RaisePixelError(Error, ErrorCode_SpecFit,            &
                             ModuleName, SubroutineName,          &
                             'xgrid is outside chebyshev interval')
      ENDIF

      ! T_0
      T(:,1) = 1.0d0

      IF(nT .GT. 1) THEN

        DO i=1,nx

          ! Transform value for grid [-1,1]
          x_new = 2.0d0*(xgrid(i)-xmin)/(xmax-xmin)-1.0d0

          ! T_1
          T(i,2) = x_new

          ! T_b-1
          DO b=3,nT
            T(i,b) = 2*x_new*T(i,b-1)-T(i,b-2)
          ENDDO
        ENDDO

      ENDIF

  END SUBROUTINE ConstructChebyshevBasis

  SUBROUTINE EvaluateChebyshev(nx,xmin,xmax,xgrid,nT,c,y,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,        INTENT(IN)    :: nx
    REAL(KIND=8),   INTENT(IN)    :: xmin, xmax
    REAL(KIND=8),   INTENT(IN)    :: xgrid(nx)
    INTEGER,        INTENT(IN)    :: nT
    REAL(KIND=8),   INTENT(IN)    :: c(nT)
    REAL(KIND=8),   INTENT(OUT)   :: y(nx)
    TYPE(ErrorType),INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: i, b
    REAL(KIND=8) :: x_new
    REAL(KIND=8) :: T(nx,nT)

    ! ============================================================
    ! EvaluateChebyshev starts here
    ! ============================================================

    CALL ConstructChebyshevBasis(nx,xmin,xmax,xgrid,nT,T,Error)
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    y(:) = 0.0d0
    DO b=1,nT
      y(:) = y(:) + c(b)*T(:,b)
    ENDDO

  END SUBROUTINE EvaluateChebyshev
  
  SUBROUTINE EvaluateChebyshevPoint(xmin,xmax,x_out,nT,c,y_out,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    REAL(KIND=8),   INTENT(IN)    :: xmin, xmax
    REAL(KIND=8),   INTENT(IN)    :: x_out
    INTEGER,        INTENT(IN)    :: nT
    REAL(KIND=8),   INTENT(IN)    :: c(nT)
    REAL(KIND=8),   INTENT(OUT)   :: y_out
    TYPE(ErrorType),INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: i, b
    REAL(KIND=8) :: x_out_arr(1)
    REAL(KIND=8) :: T(1,nT)
    ! ============================================================
    ! EvaluateChebyshevPoint starts here
    ! ============================================================
    
    x_out_arr(1) = x_out
    
    CALL ConstructChebyshevBasis(1,xmin,xmax,x_out_arr,nT,T,Error)

    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    y_out = 0.0d0
    DO b=1,nT
      y_out = y_out + c(b)*T(1,b)
    ENDDO

  END SUBROUTINE EvaluateChebyshevPoint
  
  subroutine CHEBFT(A,B,C,N,FUNC)
    real*8 A,B,C(N)
    integer :: n
    external FUNC
    !********************************************************
    !* Chebyshev fit: Given a real function FUNC(X), lower  *
    !* and upper limits of the interval [A,B] for X, and a  *
    !* maximum degree N, this routine computes the N Cheby- *
    !* shev coefficients Ck, such that FUNC(X) is approxima-*
    !* ted by:  N                                           *
    !*         [Sum Ck Tk-1(Y)] - C1/2, where X and Y are   *
    !*         k=1                                          *
    !* related by:     Y = (X - 1/2(A+B)) / (1/2(B-A))      *
    !* This routine is to be used with moderately large N   *
    !* (e.g. 30 or 50), the array of C's subsequently to be *
    !* truncated at the smaller value m such that Cm+1 and  *
    !* subsequent elements are negligible.                  *
    !********************************************************
    real*8 PI,SUM,F(NMAX)
    real*8 BMA,BPA,FAC, Y
    real*8 FUNC
    integer :: j, k

    PI=4.d0*datan(1.d0)
    BMA=HALF*(B-A); BPA=HALF*(B+A)
    do K=1,N
      Y=DCOS(PI*(K-HALF)/REAL(N,KIND=8))
      F(K)=FUNC(Y*BMA+BPA)
    end do
    FAC=TWO/REAL(N,KIND=8)
    do J=1,N
      SUM=ZERO
      do K=1,N
        SUM=SUM+F(K)*DCOS((PI*REAL(J-1,KIND=8))*((REAL(K,KIND=8)-HALF)/REAL(N,KIND=8)))
      end do
      C(J)=FAC*SUM
    end do
    return
  end subroutine CHEBFT

  real*8 function CHEBEV(A,B,C,M,X)
  !**********************************************************
  !* Chebyshev evaluation: All arguments are input. C is an *
  !* array of Chebyshev coefficients, of length M, the first*
  !* M elements of Coutput from subroutine CHEBFT (which    *
  !* must have been called with the same A and B). The Che- *
  !* byshev polynomial is evaluated at a point Y determined *
  !* from X, A and B, and the result FUNC(X) is returned as *
  !* the function value.                                    *
  !**********************************************************
  !parameter(HALF=0.5d0,TWO=2.d0,ZERO=0.d0)
  real*8 A,B,C(M),X
  real*8 D,DD,SV,Y,Y2
  integer :: j, k, m
  if ((X-A)*(X-B).gt.ZERO) STOP ' X not in range.'
  D=ZERO; DD=ZERO
  Y=(TWO*X-A-B)/(B-A)  !change of variable
  Y2=TWO*Y
  do J=M,2,-1
    SV=D
    D=Y2*D-DD+C(J)
    DD=SV
  end do
  CHEBEV=Y*D-DD+HALF*C(1)
  return
  end function CHEBEV

  subroutine CHINT(A,B,C,CINT,N)
  !**********************************************************
  !* Given A,B,C, as output from routine CHEBFT, and given  *
  !* N, the desired degree of approximation (length of C to *
  !* be used), this routine returns the array CINT, the Che-*
  !* byshev coefficients of the integral of the function    *
  !* whose coefficients are C. The constant of integration  *
  !* is set so that the integral vanishes at A.             *
  !**********************************************************
  real*8 A,B,C(N),CINT(N)
  !parameter(ONE=1.d0,QUART=0.25d0,TWO=2.d0,ZERO=0.d0)
  real*8 CON,FAC,SUM
  integer :: j, n

  CON=QUART*(B-A)
  SUM=ZERO
  FAC=ONE
  do J=2, N-1
    CINT(J)=CON*(C(J-1)-C(J+1))/REAL(J-1,KIND=8)
    SUM=SUM+FAC*CINT(J)
    FAC=-FAC
  end do
  CINT(N)=CON*C(N-1)/REAL(N-1,KIND=8)
  SUM=SUM+FAC*CINT(N)
  CINT(1)=TWO*SUM  !set the constant of integration
  return
  end subroutine CHINT
  
  subroutine CHDER(A,B,C,CDER,N)
  !**********************************************************
  !* Given A,B,C, as output from routine CHEBFT, and given  *
  !* N, the desired degree of approximation (length of C to *
  !* be used), this routine returns the array CDER, the Che-*
  !* byshev coefficients of the derivative of the function  *
  !* whose coefficients are C.                              *
  !********************************************************** 
  real*8 A,B,C(N),CDER(N)
  !parameter(TWO=2.d0,ZERO=0.d0)
  real*8 CON
  integer :: j, n
  CDER(N)=ZERO
  CDER(N-1)=TWO*REAL(N-1,KIND=8)*C(N)
  if (N.ge.3) then
    do J=N-2, 1, -1
      CDER(J)=CDER(J+2)+TWO*REAL(J,KIND=8)*C(J+1)
    end do
  end if
  CON=TWO/(B-A)
  do J=1,N                !normalize to interval B - A
    CDER(J)=CDER(J)*CON
  end do
  return
  end subroutine CHDER
  
END MODULE chebyshev_module