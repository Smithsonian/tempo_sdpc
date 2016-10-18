!>Perform matrix inversion by LU decomposition, using LAPACK routines
module m_invert2

use tell_module

private !dgetri, dgetrf, dtrtri, dswap, dtrmm, dtrmv, dscal, dgemm, dgemv, &
     !dgetf2, dlaswp, dtrsm, dger, ilaenv, lsame, dlamch, idamax

public invert2

contains

  function invert2 (amat, error) result (amatinv)

    implicit none

    !--------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  invert2
    !
    ! !DESCRIPTION:
    !>  Matrix inversion interface to LAPACK routines,
    !>  performing inversion by LU decomposition.
    !>  Replaces m_invert, which used Numerical Recipes code
    !
    ! !CALLING SEQUENCE:
    !
    !        amatinv = invert(amat)
    ! 
    ! !INPUT PARAMETERS:
    real (KIND=8), dimension(:,:), intent(in) :: amat
    !> @param[in]   amat   2D matrix to invert
    !
    ! !OUTPUT PARAMETERS:
    integer,              intent(out):: error
    !> @param[out]  error  non-zero value indicates matrix inversion failed
    !>/verbatim
    !>   Note that error value is determined by the LAPACK routines
    !>    DGETRF and DGETRI.
    !>    error <1 => illegal element value in matrix
    !>    error >1 => singular matrix
    !>/endverbatim
    !
    real (KIND=8), dimension(lbound(amat,1):ubound(amat,1), &
         lbound(amat,2):ubound(amat,2)) :: amatinv
    !> @param   amatinv   2D inverted matrix
    !
    ! !REVISION HISTORY:
    !
    !> @author  31Jul14   O'Sullivan    Initial version
    !>   4Aug15   O'Sullivan    Brought LAPACK and BLAS routines into module
    !
    !EOP
    !-------------------------------------------------------------------------

    !local variables
    integer, dimension (size(amat,1)) :: indx
    real (KIND=8), dimension (size(amat,1)**2) :: work
    real (KIND=8), dimension (size(amat,1),size(amat,2)) :: temp
    integer                           :: nsampl, lwork


    nsampl=size(amat,1)
    lwork=nsampl**2

    amatinv=0.d0
    temp=amat

    if (nsampl /= size(amat,2)) then
      call tell_error(tell_runtime_error, &
           "invert2 can't invert a non-square matrix", error)
      return
    endif

    call dgetrf(nsampl, nsampl, temp, nsampl, indx, error)
    if (error == 0) then
      call dgetri(nsampl, temp, nsampl, indx, work, lwork, error)
    endif
    if (error == 0) then
      amatinv = temp
    else if (error > 0) then
      call tell_error(tell_runtime_error, "invert2: singular matrix", error)
      return
    else
      call tell_error(tell_runtime_error, &
           "invert2: illegal value in matrix", error)
      return
    endif

  end function invert2


!! everything below this taken from LAPACK and BLAS libraries,
!! with minor editing
!!------------------------------------------------------------------------
!
!!> \brief \b DGETRI
!!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!!
!!       .. Scalar Arguments ..
!!       INTEGER            INFO, LDA, LWORK, N
!!       ..
!!       .. Array Arguments ..
!!       INTEGER            IPIV( * )
!!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DGETRI computes the inverse of a matrix using the LU factorization
!!> computed by DGETRF.
!!>
!!> This method inverts U and then computes inv(A) by solving the system
!!> inv(A)*L = inv(U) for inv(A).
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>          The order of the matrix A.  N >= 0.
!!> \endverbatim
!!>
!!> \param[in,out] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!>          On entry, the factors L and U from the factorization
!!>          A = P*L*U as computed by DGETRF.
!!>          On exit, if INFO = 0, the inverse of the original matrix A.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>          The leading dimension of the array A.  LDA >= max(1,N).
!!> \endverbatim
!!>
!!> \param[in] IPIV
!!> \verbatim
!!>          IPIV is INTEGER array, dimension (N)
!!>          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!!>          matrix was interchanged with row IPIV(i).
!!> \endverbatim
!!>
!!> \param[out] WORK
!!> \verbatim
!!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!>          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!!> \endverbatim
!!>
!!> \param[in] LWORK
!!> \verbatim
!!>          LWORK is INTEGER
!!>          The dimension of the array WORK.  LWORK >= max(1,N).
!!>          For optimal performance LWORK >= N*NB, where NB is
!!>          the optimal blocksize returned by ILAENV.
!!>
!!>          If LWORK = -1, then a workspace query is assumed; the routine
!!>          only calculates the optimal size of the WORK array, returns
!!>          this value as the first entry of the WORK array, and no error
!!>          message related to LWORK is issued by XERBLA.
!!> \endverbatim
!!>
!!> \param[out] INFO
!!> \verbatim
!!>          INFO is INTEGER
!!>          = 0:  successful exit
!!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!!>          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!!>                singular and its inverse could not be computed.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup doubleGEcomputational
!!
!!  =====================================================================
!      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!!
!!  -- LAPACK computational routine (version 3.4.0) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      INTEGER            INFO, LDA, LWORK, N
!!     ..
!!     .. Array Arguments ..
!      INTEGER            IPIV( * )
!      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION   ZERO, ONE
!      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!!     ..
!!     .. Local Scalars ..
!      LOGICAL            LQUERY
!      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, &
!                        NBMIN, NN, errstat
!!     ..
!!     .. External Functions ..
!!      INTEGER            ILAENV
!!      EXTERNAL           ILAENV
!!     ..
!!     .. External Subroutines ..
!!      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI!, XERBLA
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!!     ..
!!     .. Executable Statements ..
!!
!!     Test the input parameters.
!!
!      INFO = 0
!      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
!      LWKOPT = N*NB
!      WORK( 1 ) = LWKOPT
!      LQUERY = ( LWORK.EQ.-1 )
!      IF( N.LT.0 ) THEN
!         INFO = -1
!      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!         INFO = -3
!      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
!         INFO = -6
!      END IF
!      IF( INFO.NE.0 ) THEN
!!        CALL XERBLA( 'DGETRI', -INFO )
!         errstat=INFO
!         call tell_error(tell_invalid_parm,'DGETRI', errstat)
!         RETURN
!      ELSE IF( LQUERY ) THEN
!         RETURN
!      END IF
!!
!!     Quick return if possible
!!
!      IF( N.EQ.0 ) RETURN
!!
!!     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
!!     and the inverse is not computed.
!!
!      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
!      IF( INFO.GT.0 ) RETURN
!!
!      NBMIN = 2
!      LDWORK = N
!      IF( NB.GT.1 .AND. NB.LT.N ) THEN
!         IWS = MAX( LDWORK*NB, 1 )
!         IF( LWORK.LT.IWS ) THEN
!            NB = LWORK / LDWORK
!            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
!         END IF
!      ELSE
!         IWS = N
!      END IF
!!
!!     Solve the equation inv(A)*L = inv(U) for inv(A).
!!
!      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
!!
!!        Use unblocked code.
!!
!         DO 20 J = N, 1, -1
!!
!!           Copy current column of L to WORK and replace with zeros.
!!
!            DO 10 I = J + 1, N
!               WORK( I ) = A( I, J )
!               A( I, J ) = ZERO
!   10       CONTINUE
!!
!!           Compute current column of inv(A).
!!
!            IF( J.LT.N ) CALL DGEMV( 'No transpose', N, N-J, -ONE, &
!                 A( 1, J+1 ), LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
!   20    CONTINUE
!      ELSE
!!
!!        Use blocked code.
!!
!         NN = ( ( N-1 ) / NB )*NB + 1
!         DO 50 J = NN, 1, -NB
!            JB = MIN( NB, N-J+1 )
!!
!!           Copy current block column of L to WORK and replace with
!!           zeros.
!!
!            DO 40 JJ = J, J + JB - 1
!               DO 30 I = JJ + 1, N
!                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
!                  A( I, JJ ) = ZERO
!   30          CONTINUE
!   40       CONTINUE
!!
!!           Compute current block column of inv(A).
!!
!            IF( J+JB.LE.N ) CALL DGEMM( 'No transpose', 'No transpose', N, &
!                 JB, N-J-JB+1, -ONE, A( 1, J+JB ), LDA, WORK( J+JB ), &
!                 LDWORK, ONE, A( 1, J ), LDA )
!            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB, &
!                 ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
!   50    CONTINUE
!      END IF
!!
!!     Apply column interchanges.
!!
!      DO 60 J = N - 1, 1, -1
!         JP = IPIV( J )
!         IF( JP.NE.J ) CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
!   60 CONTINUE
!!
!      WORK( 1 ) = IWS
!      RETURN
!!
!!     End of DGETRI
!!
!      END subroutine dgetri
!
!
!
!!> \brief \b DTRTRI
!!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!!
!!       .. Scalar Arguments ..
!!       CHARACTER          DIAG, UPLO
!!       INTEGER            INFO, LDA, N
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION   A( LDA, * )
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DTRTRI computes the inverse of a real upper or lower triangular
!!> matrix A.
!!>
!!> This is the Level 3 BLAS version of the algorithm.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] UPLO
!!> \verbatim
!!>          UPLO is CHARACTER*1
!!>          = 'U':  A is upper triangular;
!!>          = 'L':  A is lower triangular.
!!> \endverbatim
!!>
!!> \param[in] DIAG
!!> \verbatim
!!>          DIAG is CHARACTER*1
!!>          = 'N':  A is non-unit triangular;
!!>          = 'U':  A is unit triangular.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>          The order of the matrix A.  N >= 0.
!!> \endverbatim
!!>
!!> \param[in,out] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!>          On entry, the triangular matrix A.  If UPLO = 'U', the
!!>          leading N-by-N upper triangular part of the array A contains
!!>          the upper triangular matrix, and the strictly lower
!!>          triangular part of A is not referenced.  If UPLO = 'L', the
!!>          leading N-by-N lower triangular part of the array A contains
!!>          the lower triangular matrix, and the strictly upper
!!>          triangular part of A is not referenced.  If DIAG = 'U', the
!!>          diagonal elements of A are also not referenced and are
!!>          assumed to be 1.
!!>          On exit, the (triangular) inverse of the original matrix, in
!!>          the same storage format.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>          The leading dimension of the array A.  LDA >= max(1,N).
!!> \endverbatim
!!>
!!> \param[out] INFO
!!> \verbatim
!!>          INFO is INTEGER
!!>          = 0: successful exit
!!>          < 0: if INFO = -i, the i-th argument had an illegal value
!!>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!!>               matrix is singular and its inverse can not be computed.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup doubleOTHERcomputational
!!
!!  =====================================================================
!      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!!
!!  -- LAPACK computational routine (version 3.4.0) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      CHARACTER          DIAG, UPLO
!      INTEGER            INFO, LDA, N
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION   A( LDA, * )
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION   ONE, ZERO
!      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!!     ..
!!     .. Local Scalars ..
!      LOGICAL            NOUNIT, UPPER
!      INTEGER            J, JB, NB, NN, errstat
!!     ..
!!     .. External Functions ..
!!     LOGICAL            LSAME
!!     INTEGER            ILAENV
!!     EXTERNAL           LSAME, ILAENV
!!     ..
!!     .. External Subroutines ..
!!     EXTERNAL           DTRMM, DTRSM, DTRTI2 !, XERBLA
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!!     ..
!!     .. Executable Statements ..
!!
!!     Test the input parameters.
!!
!      INFO = 0
!      UPPER = LSAME( UPLO, 'U' )
!      NOUNIT = LSAME( DIAG, 'N' )
!      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
!         INFO = -1
!      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
!         INFO = -2
!      ELSE IF( N.LT.0 ) THEN
!         INFO = -3
!      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!         INFO = -5
!      END IF
!      IF( INFO.NE.0 ) THEN
!!         CALL XERBLA( 'DTRTRI', -INFO )
!         errstat=INFO
!         call tell_error(tell_invalid_parm,'DTRTRI', errstat)
!         RETURN
!      END IF
!!
!!     Quick return if possible
!!
!      IF( N.EQ.0 ) RETURN
!!
!!     Check for singularity if non-unit.
!!
!      IF( NOUNIT ) THEN
!         DO 10 INFO = 1, N
!            IF( A( INFO, INFO ).EQ.ZERO ) RETURN
!   10    CONTINUE
!         INFO = 0
!      END IF
!!
!!     Determine the block size for this environment.
!!
!      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!!
!!        Use unblocked code
!!
!         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!      ELSE
!!
!!        Use blocked code
!!
!         IF( UPPER ) THEN
!!
!!           Compute inverse of upper triangular matrix
!!
!            DO 20 J = 1, N, NB
!               JB = MIN( NB, N-J+1 )
!!
!!              Compute rows 1:j-1 of current block column
!!
!               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1, &
!                    JB, ONE, A, LDA, A( 1, J ), LDA )
!               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1, &
!                    JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
!!
!!              Compute inverse of current diagonal block
!!
!               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
!   20       CONTINUE
!         ELSE
!!
!!           Compute inverse of lower triangular matrix
!!
!            NN = ( ( N-1 ) / NB )*NB + 1
!            DO 30 J = NN, 1, -NB
!               JB = MIN( NB, N-J+1 )
!               IF( J+JB.LE.N ) THEN
!!
!!                 Compute rows j+jb:n of current block column
!!
!                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG, &
!                       N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, &
!                       A( J+JB, J ), LDA )
!                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG, &
!                       N-J-JB+1, JB, -ONE, A( J, J ), LDA, &
!                       A( J+JB, J ), LDA )
!               END IF
!!
!!              Compute inverse of current diagonal block
!!
!               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
!   30       CONTINUE
!         END IF
!      END IF
!!
!      RETURN
!!
!!     End of DTRTRI
!!
!      END subroutine dtrtri
!
!
!!> \brief \b DTRTI2 computes the inverse of a triangular matrix (unblocked algorithm).
!!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!!
!!       .. Scalar Arguments ..
!!       CHARACTER          DIAG, UPLO
!!       INTEGER            INFO, LDA, N
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION   A( LDA, * )
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DTRTI2 computes the inverse of a real upper or lower triangular
!!> matrix.
!!>
!!> This is the Level 2 BLAS version of the algorithm.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] UPLO
!!> \verbatim
!!>          UPLO is CHARACTER*1
!!>          Specifies whether the matrix A is upper or lower triangular.
!!>          = 'U':  Upper triangular
!!>          = 'L':  Lower triangular
!!> \endverbatim
!!>
!!> \param[in] DIAG
!!> \verbatim
!!>          DIAG is CHARACTER*1
!!>          Specifies whether or not the matrix A is unit triangular.
!!>          = 'N':  Non-unit triangular
!!>          = 'U':  Unit triangular
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>          The order of the matrix A.  N >= 0.
!!> \endverbatim
!!>
!!> \param[in,out] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!>          On entry, the triangular matrix A.  If UPLO = 'U', the
!!>          leading n by n upper triangular part of the array A contains
!!>          the upper triangular matrix, and the strictly lower
!!>          triangular part of A is not referenced.  If UPLO = 'L', the
!!>          leading n by n lower triangular part of the array A contains
!!>          the lower triangular matrix, and the strictly upper
!!>          triangular part of A is not referenced.  If DIAG = 'U', the
!!>          diagonal elements of A are also not referenced and are
!!>          assumed to be 1.
!!>
!!>          On exit, the (triangular) inverse of the original matrix, in
!!>          the same storage format.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>          The leading dimension of the array A.  LDA >= max(1,N).
!!> \endverbatim
!!>
!!> \param[out] INFO
!!> \verbatim
!!>          INFO is INTEGER
!!>          = 0: successful exit
!!>          < 0: if INFO = -k, the k-th argument had an illegal value
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date September 2012
!!
!!> \ingroup doubleOTHERcomputational
!!
!!  =====================================================================
!      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!!
!!  -- LAPACK computational routine (version 3.4.2) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     September 2012
!!
!!     .. Scalar Arguments ..
!      CHARACTER          DIAG, UPLO
!      INTEGER            INFO, LDA, N
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION   A( LDA, * )
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION   ONE
!      PARAMETER          ( ONE = 1.0D+0 )
!!     ..
!!     .. Local Scalars ..
!      LOGICAL            NOUNIT, UPPER
!      INTEGER            J, errstat
!      DOUBLE PRECISION   AJJ
!!     ..
!!     .. External Functions ..
!!     LOGICAL            LSAME
!!     EXTERNAL           LSAME
!!     ..
!!     .. External Subroutines ..
!!     EXTERNAL           DSCAL, DTRMV !, XERBLA
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          MAX
!!     ..
!!     .. Executable Statements ..
!!
!!     Test the input parameters.
!!
!      INFO = 0
!      UPPER = LSAME( UPLO, 'U' )
!      NOUNIT = LSAME( DIAG, 'N' )
!      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
!         INFO = -1
!      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
!         INFO = -2
!      ELSE IF( N.LT.0 ) THEN
!         INFO = -3
!      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!         INFO = -5
!      END IF
!      IF( INFO.NE.0 ) THEN
!!         CALL XERBLA( 'DTRTI2', -INFO )
!         errstat=INFO
!         call tell_error(tell_invalid_parm,'DTRTI2', errstat)
!         RETURN
!      END IF
!!
!      IF( UPPER ) THEN
!!
!!        Compute inverse of upper triangular matrix.
!!
!         DO 10 J = 1, N
!            IF( NOUNIT ) THEN
!               A( J, J ) = ONE / A( J, J )
!               AJJ = -A( J, J )
!            ELSE
!               AJJ = -ONE
!            END IF
!!
!!           Compute elements 1:j-1 of j-th column.
!!
!            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA, &
!                 A( 1, J ), 1 )
!            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
!   10    CONTINUE
!      ELSE
!!
!!        Compute inverse of lower triangular matrix.
!!
!         DO 20 J = N, 1, -1
!            IF( NOUNIT ) THEN
!               A( J, J ) = ONE / A( J, J )
!               AJJ = -A( J, J )
!            ELSE
!               AJJ = -ONE
!            END IF
!            IF( J.LT.N ) THEN
!!
!!              Compute elements j+1:n of j-th column.
!!
!               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J, &
!                    A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
!               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
!            END IF
!   20    CONTINUE
!      END IF
!!
!      RETURN
!!
!!     End of DTRTI2
!!
!      END subroutine dtrti2
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!!
!!       .. Scalar Arguments ..
!!       INTEGER INCX,INCY,N
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION DX(*),DY(*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!>    interchanges two vectors.
!!>    uses unrolled loops for increments equal one.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level1
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>     jack dongarra, linpack, 3/11/78.
!!>     modified 12/3/93, array(1) declarations changed to array(*)
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE dswap(N,DX,INCX,DY,INCY)
!!
!!  -- Reference BLAS level1 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      INTEGER incx,incy,n
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION dx(*),dy(*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. Local Scalars ..
!      DOUBLE PRECISION dtemp
!      INTEGER i,ix,iy,m,mp1
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC mod
!!     ..
!      IF (n.LE.0) RETURN
!      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!!
!!       code for both increments equal to 1
!!
!!
!!       clean-up loop
!!
!         m = mod(n,3)
!         IF (m.NE.0) THEN
!            DO i = 1,m
!               dtemp = dx(i)
!               dx(i) = dy(i)
!               dy(i) = dtemp
!            END DO
!            IF (n.LT.3) RETURN
!         END IF
!         mp1 = m + 1
!         DO i = mp1,n,3
!            dtemp = dx(i)
!            dx(i) = dy(i)
!            dy(i) = dtemp
!            dtemp = dx(i+1)
!            dx(i+1) = dy(i+1)
!            dy(i+1) = dtemp
!            dtemp = dx(i+2)
!            dx(i+2) = dy(i+2)
!            dy(i+2) = dtemp
!         END DO
!      ELSE
!!
!!       code for unequal increments or equal increments not equal
!!         to 1
!!
!         ix = 1
!         iy = 1
!         IF (incx.LT.0) ix = (-n+1)*incx + 1
!         IF (incy.LT.0) iy = (-n+1)*incy + 1
!         DO i = 1,n
!            dtemp = dx(ix)
!            dx(ix) = dy(iy)
!            dy(iy) = dtemp
!            ix = ix + incx
!            iy = iy + incy
!         END DO
!      END IF
!      RETURN
!
!      END subroutine dswap
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!!
!!       .. Scalar Arguments ..
!!       DOUBLE PRECISION ALPHA
!!       INTEGER LDA,LDB,M,N
!!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION A(LDA,*),B(LDB,*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DTRMM  performs one of the matrix-matrix operations
!!>
!!>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!!>
!!> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!!>
!!>    op( A ) = A   or   op( A ) = A**T.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] SIDE
!!> \verbatim
!!>          SIDE is CHARACTER*1
!!>           On entry,  SIDE specifies whether  op( A ) multiplies B from
!!>           the left or right as follows:
!!>
!!>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!!>
!!>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!!> \endverbatim
!!>
!!> \param[in] UPLO
!!> \verbatim
!!>          UPLO is CHARACTER*1
!!>           On entry, UPLO specifies whether the matrix A is an upper or
!!>           lower triangular matrix as follows:
!!>
!!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!!>
!!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!!> \endverbatim
!!>
!!> \param[in] TRANSA
!!> \verbatim
!!>          TRANSA is CHARACTER*1
!!>           On entry, TRANSA specifies the form of op( A ) to be used in
!!>           the matrix multiplication as follows:
!!>
!!>              TRANSA = 'N' or 'n'   op( A ) = A.
!!>
!!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!!>
!!>              TRANSA = 'C' or 'c'   op( A ) = A**T.
!!> \endverbatim
!!>
!!> \param[in] DIAG
!!> \verbatim
!!>          DIAG is CHARACTER*1
!!>           On entry, DIAG specifies whether or not A is unit triangular
!!>           as follows:
!!>
!!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!!>
!!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!!>                                  triangular.
!!> \endverbatim
!!>
!!> \param[in] M
!!> \verbatim
!!>          M is INTEGER
!!>           On entry, M specifies the number of rows of B. M must be at
!!>           least zero.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>           On entry, N specifies the number of columns of B.  N must be
!!>           at least zero.
!!> \endverbatim
!!>
!!> \param[in] ALPHA
!!> \verbatim
!!>          ALPHA is DOUBLE PRECISION.
!!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!!>           zero then  A is not referenced and  B need not be set before
!!>           entry.
!!> \endverbatim
!!>
!!> \param[in] A
!!> \verbatim
!!>           A is DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!!>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!!>           upper triangular part of the array  A must contain the upper
!!>           triangular matrix  and the strictly lower triangular part of
!!>           A is not referenced.
!!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!!>           lower triangular part of the array  A must contain the lower
!!>           triangular matrix  and the strictly upper triangular part of
!!>           A is not referenced.
!!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!!>           A  are not referenced either,  but are assumed to be  unity.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>           On entry, LDA specifies the first dimension of A as declared
!!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!!>           then LDA must be at least max( 1, n ).
!!> \endverbatim
!!>
!!> \param[in,out] B
!!> \verbatim
!!>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!!>           Before entry,  the leading  m by n part of the array  B must
!!>           contain the matrix  B,  and  on exit  is overwritten  by the
!!>           transformed matrix.
!!> \endverbatim
!!>
!!> \param[in] LDB
!!> \verbatim
!!>          LDB is INTEGER
!!>           On entry, LDB specifies the first dimension of B as declared
!!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!!>           max( 1, m ).
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level3
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  Level 3 Blas routine.
!!>
!!>  -- Written on 8-February-1989.
!!>     Jack Dongarra, Argonne National Laboratory.
!!>     Iain Duff, AERE Harwell.
!!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE dtrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!!
!!  -- Reference BLAS level3 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      DOUBLE PRECISION alpha
!      INTEGER lda,ldb,m,n
!      CHARACTER diag,side,transa,uplo
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION a(lda,*),b(ldb,*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. External Functions ..
!!     LOGICAL lsame
!!     EXTERNAL lsame
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC max
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION temp
!      INTEGER i,info,j,k,nrowa
!      LOGICAL lside,nounit,upper
!!     ..
!!     .. Parameters ..
!      DOUBLE PRECISION one,zero
!      parameter(one=1.0d+0,zero=0.0d+0)
!!     ..
!!
!!     Test the input parameters.
!!
!      lside = lsame(side,'L')
!      IF (lside) THEN
!          nrowa = m
!      ELSE
!          nrowa = n
!      END IF
!      nounit = lsame(diag,'N')
!      upper = lsame(uplo,'U')
!!
!      info = 0
!      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
!          info = 1
!      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
!          info = 2
!      ELSE IF ((.NOT.lsame(transa,'N')) .AND. &
!           (.NOT.lsame(transa,'T')) .AND. &
!           (.NOT.lsame(transa,'C'))) THEN
!          info = 3
!      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
!          info = 4
!      ELSE IF (m.LT.0) THEN
!          info = 5
!      ELSE IF (n.LT.0) THEN
!          info = 6
!      ELSE IF (lda.LT.max(1,nrowa)) THEN
!          info = 9
!      ELSE IF (ldb.LT.max(1,m)) THEN
!          info = 11
!      END IF
!      IF (info.NE.0) THEN
!!          CALL xerbla('DTRMM ',info)
!          call tell_error(tell_invalid_parm, 'DTRMM ',info)
!          RETURN
!      END IF
!!
!!     Quick return if possible.
!!
!      IF (m.EQ.0 .OR. n.EQ.0) RETURN
!!
!!     And when  alpha.eq.zero.
!!
!      IF (alpha.EQ.zero) THEN
!          DO 20 j = 1,n
!              DO 10 i = 1,m
!                  b(i,j) = zero
!   10         CONTINUE
!   20     CONTINUE
!          RETURN
!      END IF
!!
!!     Start the operations.
!!
!      IF (lside) THEN
!          IF (lsame(transa,'N')) THEN
!!
!!           Form  B := alpha*A*B.
!!
!              IF (upper) THEN
!                  DO 50 j = 1,n
!                      DO 40 k = 1,m
!                          IF (b(k,j).NE.zero) THEN
!                              temp = alpha*b(k,j)
!                              DO 30 i = 1,k - 1
!                                  b(i,j) = b(i,j) + temp*a(i,k)
!   30                         CONTINUE
!                              IF (nounit) temp = temp*a(k,k)
!                              b(k,j) = temp
!                          END IF
!   40                 CONTINUE
!   50             CONTINUE
!              ELSE
!                  DO 80 j = 1,n
!                      DO 70 k = m,1,-1
!                          IF (b(k,j).NE.zero) THEN
!                              temp = alpha*b(k,j)
!                              b(k,j) = temp
!                              IF (nounit) b(k,j) = b(k,j)*a(k,k)
!                              DO 60 i = k + 1,m
!                                  b(i,j) = b(i,j) + temp*a(i,k)
!   60                         CONTINUE
!                          END IF
!   70                 CONTINUE
!   80             CONTINUE
!              END IF
!          ELSE
!!
!!           Form  B := alpha*A**T*B.
!!
!              IF (upper) THEN
!                  DO 110 j = 1,n
!                      DO 100 i = m,1,-1
!                          temp = b(i,j)
!                          IF (nounit) temp = temp*a(i,i)
!                          DO 90 k = 1,i - 1
!                              temp = temp + a(k,i)*b(k,j)
!   90                     CONTINUE
!                          b(i,j) = alpha*temp
!  100                 CONTINUE
!  110             CONTINUE
!              ELSE
!                  DO 140 j = 1,n
!                      DO 130 i = 1,m
!                          temp = b(i,j)
!                          IF (nounit) temp = temp*a(i,i)
!                          DO 120 k = i + 1,m
!                              temp = temp + a(k,i)*b(k,j)
!  120                     CONTINUE
!                          b(i,j) = alpha*temp
!  130                 CONTINUE
!  140             CONTINUE
!              END IF
!          END IF
!      ELSE
!          IF (lsame(transa,'N')) THEN
!!
!!           Form  B := alpha*B*A.
!!
!              IF (upper) THEN
!                  DO 180 j = n,1,-1
!                      temp = alpha
!                      IF (nounit) temp = temp*a(j,j)
!                      DO 150 i = 1,m
!                          b(i,j) = temp*b(i,j)
!  150                 CONTINUE
!                      DO 170 k = 1,j - 1
!                          IF (a(k,j).NE.zero) THEN
!                              temp = alpha*a(k,j)
!                              DO 160 i = 1,m
!                                  b(i,j) = b(i,j) + temp*b(i,k)
!  160                         CONTINUE
!                          END IF
!  170                 CONTINUE
!  180             CONTINUE
!              ELSE
!                  DO 220 j = 1,n
!                      temp = alpha
!                      IF (nounit) temp = temp*a(j,j)
!                      DO 190 i = 1,m
!                          b(i,j) = temp*b(i,j)
!  190                 CONTINUE
!                      DO 210 k = j + 1,n
!                          IF (a(k,j).NE.zero) THEN
!                              temp = alpha*a(k,j)
!                              DO 200 i = 1,m
!                                  b(i,j) = b(i,j) + temp*b(i,k)
!  200                         CONTINUE
!                          END IF
!  210                 CONTINUE
!  220             CONTINUE
!              END IF
!          ELSE
!!
!!           Form  B := alpha*B*A**T.
!!
!              IF (upper) THEN
!                  DO 260 k = 1,n
!                      DO 240 j = 1,k - 1
!                          IF (a(j,k).NE.zero) THEN
!                              temp = alpha*a(j,k)
!                              DO 230 i = 1,m
!                                  b(i,j) = b(i,j) + temp*b(i,k)
!  230                         CONTINUE
!                          END IF
!  240                 CONTINUE
!                      temp = alpha
!                      IF (nounit) temp = temp*a(k,k)
!                      IF (temp.NE.one) THEN
!                          DO 250 i = 1,m
!                              b(i,k) = temp*b(i,k)
!  250                     CONTINUE
!                      END IF
!  260             CONTINUE
!              ELSE
!                  DO 300 k = n,1,-1
!                      DO 280 j = k + 1,n
!                          IF (a(j,k).NE.zero) THEN
!                              temp = alpha*a(j,k)
!                              DO 270 i = 1,m
!                                  b(i,j) = b(i,j) + temp*b(i,k)
!  270                         CONTINUE
!                          END IF
!  280                 CONTINUE
!                      temp = alpha
!                      IF (nounit) temp = temp*a(k,k)
!                      IF (temp.NE.one) THEN
!                          DO 290 i = 1,m
!                              b(i,k) = temp*b(i,k)
!  290                     CONTINUE
!                      END IF
!  300             CONTINUE
!              END IF
!          END IF
!      END IF
!!
!      RETURN
!!
!!     End of DTRMM .
!!
!      END subroutine dtrmm
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DSCAL(N,DA,DX,INCX)
!!
!!       .. Scalar Arguments ..
!!       DOUBLE PRECISION DA
!!       INTEGER INCX,N
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION DX(*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!>    DSCAL scales a vector by a constant.
!!>    uses unrolled loops for increment equal to one.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level1
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>     jack dongarra, linpack, 3/11/78.
!!>     modified 3/93 to return if incx .le. 0.
!!>     modified 12/3/93, array(1) declarations changed to array(*)
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE dscal(N,DA,DX,INCX)
!!
!!  -- Reference BLAS level1 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      DOUBLE PRECISION da
!      INTEGER incx,n
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION dx(*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. Local Scalars ..
!      INTEGER i,m,mp1,nincx
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC mod
!!     ..
!      IF (n.LE.0 .OR. incx.LE.0) RETURN
!      IF (incx.EQ.1) THEN
!!
!!        code for increment equal to 1
!!
!!
!!        clean-up loop
!!
!         m = mod(n,5)
!         IF (m.NE.0) THEN
!            DO i = 1,m
!               dx(i) = da*dx(i)
!            END DO
!            IF (n.LT.5) RETURN
!         END IF
!         mp1 = m + 1
!         DO i = mp1,n,5
!            dx(i) = da*dx(i)
!            dx(i+1) = da*dx(i+1)
!            dx(i+2) = da*dx(i+2)
!            dx(i+3) = da*dx(i+3)
!            dx(i+4) = da*dx(i+4)
!         END DO
!      ELSE
!!
!!        code for increment not equal to 1
!!
!         nincx = n*incx
!         DO i = 1,nincx,incx
!            dx(i) = da*dx(i)
!         END DO
!      END IF
!      RETURN
!      END subroutine dscal
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!!
!!       .. Scalar Arguments ..
!!       INTEGER INCX,LDA,N
!!       CHARACTER DIAG,TRANS,UPLO
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION A(LDA,*),X(*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DTRMV  performs one of the matrix-vector operations
!!>
!!>    x := A*x,   or   x := A**T*x,
!!>
!!> where x is an n element vector and  A is an n by n unit, or non-unit,
!!> upper or lower triangular matrix.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] UPLO
!!> \verbatim
!!>          UPLO is CHARACTER*1
!!>           On entry, UPLO specifies whether the matrix is an upper or
!!>           lower triangular matrix as follows:
!!>
!!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!!>
!!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!!> \endverbatim
!!>
!!> \param[in] TRANS
!!> \verbatim
!!>          TRANS is CHARACTER*1
!!>           On entry, TRANS specifies the operation to be performed as
!!>           follows:
!!>
!!>              TRANS = 'N' or 'n'   x := A*x.
!!>
!!>              TRANS = 'T' or 't'   x := A**T*x.
!!>
!!>              TRANS = 'C' or 'c'   x := A**T*x.
!!> \endverbatim
!!>
!!> \param[in] DIAG
!!> \verbatim
!!>          DIAG is CHARACTER*1
!!>           On entry, DIAG specifies whether or not A is unit
!!>           triangular as follows:
!!>
!!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!!>
!!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!!>                                  triangular.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>           On entry, N specifies the order of the matrix A.
!!>           N must be at least zero.
!!> \endverbatim
!!>
!!> \param[in] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!!>           upper triangular part of the array A must contain the upper
!!>           triangular matrix and the strictly lower triangular part of
!!>           A is not referenced.
!!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!!>           lower triangular part of the array A must contain the lower
!!>           triangular matrix and the strictly upper triangular part of
!!>           A is not referenced.
!!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!!>           A are not referenced either, but are assumed to be unity.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>           On entry, LDA specifies the first dimension of A as declared
!!>           in the calling (sub) program. LDA must be at least
!!>           max( 1, n ).
!!> \endverbatim
!!>
!!> \param[in,out] X
!!> \verbatim
!!>          X is DOUBLE PRECISION array of dimension at least
!!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!!>           Before entry, the incremented array X must contain the n
!!>           element vector x. On exit, X is overwritten with the
!!>           tranformed vector x.
!!> \endverbatim
!!>
!!> \param[in] INCX
!!> \verbatim
!!>          INCX is INTEGER
!!>           On entry, INCX specifies the increment for the elements of
!!>           X. INCX must not be zero.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level2
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  Level 2 Blas routine.
!!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!!>
!!>  -- Written on 22-October-1986.
!!>     Jack Dongarra, Argonne National Lab.
!!>     Jeremy Du Croz, Nag Central Office.
!!>     Sven Hammarling, Nag Central Office.
!!>     Richard Hanson, Sandia National Labs.
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE dtrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!!
!!  -- Reference BLAS level2 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      INTEGER incx,lda,n
!      CHARACTER diag,trans,uplo
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION a(lda,*),x(*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION zero
!      parameter(zero=0.0d+0)
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION temp
!      INTEGER i,info,ix,j,jx,kx
!      LOGICAL nounit
!!     ..
!!     .. External Functions ..
!!     LOGICAL lsame
!!     EXTERNAL lsame
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC max
!!     ..
!!
!!     Test the input parameters.
!!
!      info = 0
!      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
!          info = 1
!      ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND. &
!           .NOT.lsame(trans,'C')) THEN
!          info = 2
!      ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
!          info = 3
!      ELSE IF (n.LT.0) THEN
!          info = 4
!      ELSE IF (lda.LT.max(1,n)) THEN
!          info = 6
!      ELSE IF (incx.EQ.0) THEN
!          info = 8
!      END IF
!      IF (info.NE.0) THEN
!!          CALL xerbla('DTRMV ',info)
!          call tell_error(tell_invalid_parm,'DTRMV ',info)
!          RETURN
!      END IF
!!
!!     Quick return if possible.
!!
!      IF (n.EQ.0) RETURN
!!
!      nounit = lsame(diag,'N')
!!
!!     Set up the start point in X if the increment is not unity. This
!!     will be  ( N - 1 )*INCX  too small for descending loops.
!!
!      IF (incx.LE.0) THEN
!          kx = 1 - (n-1)*incx
!      ELSE IF (incx.NE.1) THEN
!          kx = 1
!      END IF
!!
!!     Start the operations. In this version the elements of A are
!!     accessed sequentially with one pass through A.
!!
!      IF (lsame(trans,'N')) THEN
!!
!!        Form  x := A*x.
!!
!          IF (lsame(uplo,'U')) THEN
!              IF (incx.EQ.1) THEN
!                  DO 20 j = 1,n
!                      IF (x(j).NE.zero) THEN
!                          temp = x(j)
!                          DO 10 i = 1,j - 1
!                              x(i) = x(i) + temp*a(i,j)
!   10                     CONTINUE
!                          IF (nounit) x(j) = x(j)*a(j,j)
!                      END IF
!   20             CONTINUE
!              ELSE
!                  jx = kx
!                  DO 40 j = 1,n
!                      IF (x(jx).NE.zero) THEN
!                          temp = x(jx)
!                          ix = kx
!                          DO 30 i = 1,j - 1
!                              x(ix) = x(ix) + temp*a(i,j)
!                              ix = ix + incx
!   30                     CONTINUE
!                          IF (nounit) x(jx) = x(jx)*a(j,j)
!                      END IF
!                      jx = jx + incx
!   40             CONTINUE
!              END IF
!          ELSE
!              IF (incx.EQ.1) THEN
!                  DO 60 j = n,1,-1
!                      IF (x(j).NE.zero) THEN
!                          temp = x(j)
!                          DO 50 i = n,j + 1,-1
!                              x(i) = x(i) + temp*a(i,j)
!   50                     CONTINUE
!                          IF (nounit) x(j) = x(j)*a(j,j)
!                      END IF
!   60             CONTINUE
!              ELSE
!                  kx = kx + (n-1)*incx
!                  jx = kx
!                  DO 80 j = n,1,-1
!                      IF (x(jx).NE.zero) THEN
!                          temp = x(jx)
!                          ix = kx
!                          DO 70 i = n,j + 1,-1
!                              x(ix) = x(ix) + temp*a(i,j)
!                              ix = ix - incx
!   70                     CONTINUE
!                          IF (nounit) x(jx) = x(jx)*a(j,j)
!                      END IF
!                      jx = jx - incx
!   80             CONTINUE
!              END IF
!          END IF
!      ELSE
!!
!!        Form  x := A**T*x.
!!
!          IF (lsame(uplo,'U')) THEN
!              IF (incx.EQ.1) THEN
!                  DO 100 j = n,1,-1
!                      temp = x(j)
!                      IF (nounit) temp = temp*a(j,j)
!                      DO 90 i = j - 1,1,-1
!                          temp = temp + a(i,j)*x(i)
!   90                 CONTINUE
!                      x(j) = temp
!  100             CONTINUE
!              ELSE
!                  jx = kx + (n-1)*incx
!                  DO 120 j = n,1,-1
!                      temp = x(jx)
!                      ix = jx
!                      IF (nounit) temp = temp*a(j,j)
!                      DO 110 i = j - 1,1,-1
!                          ix = ix - incx
!                          temp = temp + a(i,j)*x(ix)
!  110                 CONTINUE
!                      x(jx) = temp
!                      jx = jx - incx
!  120             CONTINUE
!              END IF
!          ELSE
!              IF (incx.EQ.1) THEN
!                  DO 140 j = 1,n
!                      temp = x(j)
!                      IF (nounit) temp = temp*a(j,j)
!                      DO 130 i = j + 1,n
!                          temp = temp + a(i,j)*x(i)
!  130                 CONTINUE
!                      x(j) = temp
!  140             CONTINUE
!              ELSE
!                  jx = kx
!                  DO 160 j = 1,n
!                      temp = x(jx)
!                      ix = jx
!                      IF (nounit) temp = temp*a(j,j)
!                      DO 150 i = j + 1,n
!                          ix = ix + incx
!                          temp = temp + a(i,j)*x(ix)
!  150                 CONTINUE
!                      x(jx) = temp
!                      jx = jx + incx
!  160             CONTINUE
!              END IF
!          END IF
!      END IF
!!
!      RETURN
!!
!!     End of DTRMV .
!!
!      END subroutine dtrmv
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!
!!       .. Scalar Arguments ..
!!       DOUBLE PRECISION ALPHA,BETA
!!       INTEGER K,LDA,LDB,LDC,M,N
!!       CHARACTER TRANSA,TRANSB
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DGEMM  performs one of the matrix-matrix operations
!!>
!!>    C := alpha*op( A )*op( B ) + beta*C,
!!>
!!> where  op( X ) is one of
!!>
!!>    op( X ) = X   or   op( X ) = X**T,
!!>
!!> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!!> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] TRANSA
!!> \verbatim
!!>          TRANSA is CHARACTER*1
!!>           On entry, TRANSA specifies the form of op( A ) to be used in
!!>           the matrix multiplication as follows:
!!>
!!>              TRANSA = 'N' or 'n',  op( A ) = A.
!!>
!!>              TRANSA = 'T' or 't',  op( A ) = A**T.
!!>
!!>              TRANSA = 'C' or 'c',  op( A ) = A**T.
!!> \endverbatim
!!>
!!> \param[in] TRANSB
!!> \verbatim
!!>          TRANSB is CHARACTER*1
!!>           On entry, TRANSB specifies the form of op( B ) to be used in
!!>           the matrix multiplication as follows:
!!>
!!>              TRANSB = 'N' or 'n',  op( B ) = B.
!!>
!!>              TRANSB = 'T' or 't',  op( B ) = B**T.
!!>
!!>              TRANSB = 'C' or 'c',  op( B ) = B**T.
!!> \endverbatim
!!>
!!> \param[in] M
!!> \verbatim
!!>          M is INTEGER
!!>           On entry,  M  specifies  the number  of rows  of the  matrix
!!>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>           On entry,  N  specifies the number  of columns of the matrix
!!>           op( B ) and the number of columns of the matrix C. N must be
!!>           at least zero.
!!> \endverbatim
!!>
!!> \param[in] K
!!> \verbatim
!!>          K is INTEGER
!!>           On entry,  K  specifies  the number of columns of the matrix
!!>           op( A ) and the number of rows of the matrix op( B ). K must
!!>           be at least  zero.
!!> \endverbatim
!!>
!!> \param[in] ALPHA
!!> \verbatim
!!>          ALPHA is DOUBLE PRECISION.
!!>           On entry, ALPHA specifies the scalar alpha.
!!> \endverbatim
!!>
!!> \param[in] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!!>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!!>           part of the array  A  must contain the matrix  A,  otherwise
!!>           the leading  k by m  part of the array  A  must contain  the
!!>           matrix A.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>           On entry, LDA specifies the first dimension of A as declared
!!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!!>           least  max( 1, k ).
!!> \endverbatim
!!>
!!> \param[in] B
!!> \verbatim
!!>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!!>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!!>           part of the array  B  must contain the matrix  B,  otherwise
!!>           the leading  n by k  part of the array  B  must contain  the
!!>           matrix B.
!!> \endverbatim
!!>
!!> \param[in] LDB
!!> \verbatim
!!>          LDB is INTEGER
!!>           On entry, LDB specifies the first dimension of B as declared
!!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!!>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!!>           least  max( 1, n ).
!!> \endverbatim
!!>
!!> \param[in] BETA
!!> \verbatim
!!>          BETA is DOUBLE PRECISION.
!!>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!!>           supplied as zero then C need not be set on input.
!!> \endverbatim
!!>
!!> \param[in,out] C
!!> \verbatim
!!>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!!>           Before entry, the leading  m by n  part of the array  C must
!!>           contain the matrix  C,  except when  beta  is zero, in which
!!>           case C need not be set on entry.
!!>           On exit, the array  C  is overwritten by the  m by n  matrix
!!>           ( alpha*op( A )*op( B ) + beta*C ).
!!> \endverbatim
!!>
!!> \param[in] LDC
!!> \verbatim
!!>          LDC is INTEGER
!!>           On entry, LDC specifies the first dimension of C as declared
!!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!!>           max( 1, m ).
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level3
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  Level 3 Blas routine.
!!>
!!>  -- Written on 8-February-1989.
!!>     Jack Dongarra, Argonne National Laboratory.
!!>     Iain Duff, AERE Harwell.
!!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!!
!!  -- Reference BLAS level3 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      DOUBLE PRECISION alpha,beta
!      INTEGER k,lda,ldb,ldc,m,n
!      CHARACTER transa,transb
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION a(lda,*),b(ldb,*),c(ldc,*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. External Functions ..
!!     LOGICAL lsame
!!     EXTERNAL lsame
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC max
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION temp
!      INTEGER i,info,j,l,ncola,nrowa,nrowb
!      LOGICAL nota,notb
!!     ..
!!     .. Parameters ..
!      DOUBLE PRECISION one,zero
!      parameter(one=1.0d+0,zero=0.0d+0)
!!     ..
!!
!!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!!     and  columns of  A  and the  number of  rows  of  B  respectively.
!!
!      nota = lsame(transa,'N')
!      notb = lsame(transb,'N')
!      IF (nota) THEN
!          nrowa = m
!          ncola = k
!      ELSE
!          nrowa = k
!          ncola = m
!      END IF
!      IF (notb) THEN
!          nrowb = k
!      ELSE
!          nrowb = n
!      END IF
!!
!!     Test the input parameters.
!!
!      info = 0
!      IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND. &
!           (.NOT.lsame(transa,'T'))) THEN
!          info = 1
!      ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND. &
!           (.NOT.lsame(transb,'T'))) THEN
!          info = 2
!      ELSE IF (m.LT.0) THEN
!          info = 3
!      ELSE IF (n.LT.0) THEN
!          info = 4
!      ELSE IF (k.LT.0) THEN
!          info = 5
!      ELSE IF (lda.LT.max(1,nrowa)) THEN
!          info = 8
!      ELSE IF (ldb.LT.max(1,nrowb)) THEN
!          info = 10
!      ELSE IF (ldc.LT.max(1,m)) THEN
!          info = 13
!      END IF
!      IF (info.NE.0) THEN
!!          CALL xerbla('DGEMM ',info)
!          call tell_error(tell_invalid_parm,'DGEMM ',info)
!          RETURN
!      END IF
!!
!!     Quick return if possible.
!!
!      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. &
!           (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
!!
!!     And if  alpha.eq.zero.
!!
!      IF (alpha.EQ.zero) THEN
!          IF (beta.EQ.zero) THEN
!              DO 20 j = 1,n
!                  DO 10 i = 1,m
!                      c(i,j) = zero
!   10             CONTINUE
!   20         CONTINUE
!          ELSE
!              DO 40 j = 1,n
!                  DO 30 i = 1,m
!                      c(i,j) = beta*c(i,j)
!   30             CONTINUE
!   40         CONTINUE
!          END IF
!          RETURN
!      END IF
!!
!!     Start the operations.
!!
!      IF (notb) THEN
!          IF (nota) THEN
!!
!!           Form  C := alpha*A*B + beta*C.
!!
!              DO 90 j = 1,n
!                  IF (beta.EQ.zero) THEN
!                      DO 50 i = 1,m
!                          c(i,j) = zero
!   50                 CONTINUE
!                  ELSE IF (beta.NE.one) THEN
!                      DO 60 i = 1,m
!                          c(i,j) = beta*c(i,j)
!   60                 CONTINUE
!                  END IF
!                  DO 80 l = 1,k
!                      IF (b(l,j).NE.zero) THEN
!                          temp = alpha*b(l,j)
!                          DO 70 i = 1,m
!                              c(i,j) = c(i,j) + temp*a(i,l)
!   70                     CONTINUE
!                      END IF
!   80             CONTINUE
!   90         CONTINUE
!          ELSE
!!
!!           Form  C := alpha*A**T*B + beta*C
!!
!              DO 120 j = 1,n
!                  DO 110 i = 1,m
!                      temp = zero
!                      DO 100 l = 1,k
!                          temp = temp + a(l,i)*b(l,j)
!  100                 CONTINUE
!                      IF (beta.EQ.zero) THEN
!                          c(i,j) = alpha*temp
!                      ELSE
!                          c(i,j) = alpha*temp + beta*c(i,j)
!                      END IF
!  110             CONTINUE
!  120         CONTINUE
!          END IF
!      ELSE
!          IF (nota) THEN
!!
!!           Form  C := alpha*A*B**T + beta*C
!!
!              DO 170 j = 1,n
!                  IF (beta.EQ.zero) THEN
!                      DO 130 i = 1,m
!                          c(i,j) = zero
!  130                 CONTINUE
!                  ELSE IF (beta.NE.one) THEN
!                      DO 140 i = 1,m
!                          c(i,j) = beta*c(i,j)
!  140                 CONTINUE
!                  END IF
!                  DO 160 l = 1,k
!                      IF (b(j,l).NE.zero) THEN
!                          temp = alpha*b(j,l)
!                          DO 150 i = 1,m
!                              c(i,j) = c(i,j) + temp*a(i,l)
!  150                     CONTINUE
!                      END IF
!  160             CONTINUE
!  170         CONTINUE
!          ELSE
!!
!!           Form  C := alpha*A**T*B**T + beta*C
!!
!              DO 200 j = 1,n
!                  DO 190 i = 1,m
!                      temp = zero
!                      DO 180 l = 1,k
!                          temp = temp + a(l,i)*b(j,l)
!  180                 CONTINUE
!                      IF (beta.EQ.zero) THEN
!                          c(i,j) = alpha*temp
!                      ELSE
!                          c(i,j) = alpha*temp + beta*c(i,j)
!                      END IF
!  190             CONTINUE
!  200         CONTINUE
!          END IF
!      END IF
!!
!      RETURN
!!
!!     End of DGEMM .
!!
!      END subroutine dgemm
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!!
!!       .. Scalar Arguments ..
!!       DOUBLE PRECISION ALPHA,BETA
!!       INTEGER INCX,INCY,LDA,M,N
!!       CHARACTER TRANS
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DGEMV  performs one of the matrix-vector operations
!!>
!!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!!>
!!> where alpha and beta are scalars, x and y are vectors and A is an
!!> m by n matrix.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] TRANS
!!> \verbatim
!!>          TRANS is CHARACTER*1
!!>           On entry, TRANS specifies the operation to be performed as
!!>           follows:
!!>
!!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!!>
!!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!!>
!!>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!!> \endverbatim
!!>
!!> \param[in] M
!!> \verbatim
!!>          M is INTEGER
!!>           On entry, M specifies the number of rows of the matrix A.
!!>           M must be at least zero.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>           On entry, N specifies the number of columns of the matrix A.
!!>           N must be at least zero.
!!> \endverbatim
!!>
!!> \param[in] ALPHA
!!> \verbatim
!!>          ALPHA is DOUBLE PRECISION.
!!>           On entry, ALPHA specifies the scalar alpha.
!!> \endverbatim
!!>
!!> \param[in] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!!>           Before entry, the leading m by n part of the array A must
!!>           contain the matrix of coefficients.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>           On entry, LDA specifies the first dimension of A as declared
!!>           in the calling (sub) program. LDA must be at least
!!>           max( 1, m ).
!!> \endverbatim
!!>
!!> \param[in] X
!!> \verbatim
!!>          X is DOUBLE PRECISION array of DIMENSION at least
!!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!!>           and at least
!!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!!>           Before entry, the incremented array X must contain the
!!>           vector x.
!!> \endverbatim
!!>
!!> \param[in] INCX
!!> \verbatim
!!>          INCX is INTEGER
!!>           On entry, INCX specifies the increment for the elements of
!!>           X. INCX must not be zero.
!!> \endverbatim
!!>
!!> \param[in] BETA
!!> \verbatim
!!>          BETA is DOUBLE PRECISION.
!!>           On entry, BETA specifies the scalar beta. When BETA is
!!>           supplied as zero then Y need not be set on input.
!!> \endverbatim
!!>
!!> \param[in,out] Y
!!> \verbatim
!!>          Y is DOUBLE PRECISION array of DIMENSION at least
!!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!!>           and at least
!!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!!>           Before entry with BETA non-zero, the incremented array Y
!!>           must contain the vector y. On exit, Y is overwritten by the
!!>           updated vector y.
!!> \endverbatim
!!>
!!> \param[in] INCY
!!> \verbatim
!!>          INCY is INTEGER
!!>           On entry, INCY specifies the increment for the elements of
!!>           Y. INCY must not be zero.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level2
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  Level 2 Blas routine.
!!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!!>
!!>  -- Written on 22-October-1986.
!!>     Jack Dongarra, Argonne National Lab.
!!>     Jeremy Du Croz, Nag Central Office.
!!>     Sven Hammarling, Nag Central Office.
!!>     Richard Hanson, Sandia National Labs.
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!!
!!  -- Reference BLAS level2 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      DOUBLE PRECISION alpha,beta
!      INTEGER incx,incy,lda,m,n
!      CHARACTER trans
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION a(lda,*),x(*),y(*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION one,zero
!      parameter(one=1.0d+0,zero=0.0d+0)
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION temp
!      INTEGER i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
!!     ..
!!     .. External Functions ..
!!     LOGICAL lsame
!!     EXTERNAL lsame
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC max
!!     ..
!!
!!     Test the input parameters.
!!
!      info = 0
!      IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') &
!           .AND. .NOT.lsame(trans,'C')) THEN
!          info = 1
!      ELSE IF (m.LT.0) THEN
!          info = 2
!      ELSE IF (n.LT.0) THEN
!          info = 3
!      ELSE IF (lda.LT.max(1,m)) THEN
!          info = 6
!      ELSE IF (incx.EQ.0) THEN
!          info = 8
!      ELSE IF (incy.EQ.0) THEN
!          info = 11
!      END IF
!      IF (info.NE.0) THEN
!!          CALL xerbla('DGEMV ',info)
!          call tell_error(tell_invalid_parm,'DGEMV',info)
!          RETURN
!      END IF
!!
!!     Quick return if possible.
!!
!      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. &
!         ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
!!
!!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!!     up the start points in  X  and  Y.
!!
!      IF (lsame(trans,'N')) THEN
!          lenx = n
!          leny = m
!      ELSE
!          lenx = m
!          leny = n
!      END IF
!      IF (incx.GT.0) THEN
!          kx = 1
!      ELSE
!          kx = 1 - (lenx-1)*incx
!      END IF
!      IF (incy.GT.0) THEN
!          ky = 1
!      ELSE
!          ky = 1 - (leny-1)*incy
!      END IF
!!
!!     Start the operations. In this version the elements of A are
!!     accessed sequentially with one pass through A.
!!
!!     First form  y := beta*y.
!!
!      IF (beta.NE.one) THEN
!          IF (incy.EQ.1) THEN
!              IF (beta.EQ.zero) THEN
!                  DO 10 i = 1,leny
!                      y(i) = zero
!   10             CONTINUE
!              ELSE
!                  DO 20 i = 1,leny
!                      y(i) = beta*y(i)
!   20             CONTINUE
!              END IF
!          ELSE
!              iy = ky
!              IF (beta.EQ.zero) THEN
!                  DO 30 i = 1,leny
!                      y(iy) = zero
!                      iy = iy + incy
!   30             CONTINUE
!              ELSE
!                  DO 40 i = 1,leny
!                      y(iy) = beta*y(iy)
!                      iy = iy + incy
!   40             CONTINUE
!              END IF
!          END IF
!      END IF
!      IF (alpha.EQ.zero) RETURN
!      IF (lsame(trans,'N')) THEN
!!
!!        Form  y := alpha*A*x + y.
!!
!          jx = kx
!          IF (incy.EQ.1) THEN
!              DO 60 j = 1,n
!                  IF (x(jx).NE.zero) THEN
!                      temp = alpha*x(jx)
!                      DO 50 i = 1,m
!                          y(i) = y(i) + temp*a(i,j)
!   50                 CONTINUE
!                  END IF
!                  jx = jx + incx
!   60         CONTINUE
!          ELSE
!              DO 80 j = 1,n
!                  IF (x(jx).NE.zero) THEN
!                      temp = alpha*x(jx)
!                      iy = ky
!                      DO 70 i = 1,m
!                          y(iy) = y(iy) + temp*a(i,j)
!                          iy = iy + incy
!   70                 CONTINUE
!                  END IF
!                  jx = jx + incx
!   80         CONTINUE
!          END IF
!      ELSE
!!
!!        Form  y := alpha*A**T*x + y.
!!
!          jy = ky
!          IF (incx.EQ.1) THEN
!              DO 100 j = 1,n
!                  temp = zero
!                  DO 90 i = 1,m
!                      temp = temp + a(i,j)*x(i)
!   90             CONTINUE
!                  y(jy) = y(jy) + alpha*temp
!                  jy = jy + incy
!  100         CONTINUE
!          ELSE
!              DO 120 j = 1,n
!                  temp = zero
!                  ix = kx
!                  DO 110 i = 1,m
!                      temp = temp + a(i,j)*x(ix)
!                      ix = ix + incx
!  110             CONTINUE
!                  y(jy) = y(jy) + alpha*temp
!                  jy = jy + incy
!  120         CONTINUE
!          END IF
!      END IF
!!
!      RETURN
!!
!!     End of DGEMV .
!!
!      END subroutine dgemv
!
!
!
!!  Definition:
!!  ===========
!!
!!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!!
!!       .. Scalar Arguments ..
!!       CHARACTER*( * )    NAME, OPTS
!!       INTEGER            ISPEC, N1, N2, N3, N4
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> ILAENV is called from the LAPACK routines to choose problem-dependent
!!> parameters for the local environment.  See ISPEC for a description of
!!> the parameters.
!!>
!!> ILAENV returns an INTEGER
!!> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!!> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!!>
!!> This version provides a set of parameters which should give good,
!!> but not optimal, performance on many of the currently available
!!> computers.  Users are encouraged to modify this subroutine to set
!!> the tuning parameters for their particular machine using the option
!!> and problem size information in the arguments.
!!>
!!> This routine will not function correctly if it is converted to all
!!> lower case.  Converting it to all upper case is allowed.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] ISPEC
!!> \verbatim
!!>          ISPEC is INTEGER
!!>          Specifies the parameter to be returned as the value of
!!>          ILAENV.
!!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!!>               algorithm will give the best performance.
!!>          = 2: the minimum block size for which the block routine
!!>               should be used; if the usable block size is less than
!!>               this value, an unblocked routine should be used.
!!>          = 3: the crossover point (in a block routine, for N less
!!>               than this value, an unblocked routine should be used)
!!>          = 4: the number of shifts, used in the nonsymmetric
!!>               eigenvalue routines (DEPRECATED)
!!>          = 5: the minimum column dimension for blocking to be used;
!!>               rectangular blocks must have dimension at least k by m,
!!>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!!>          = 6: the crossover point for the SVD (when reducing an m by n
!!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!!>               this value, a QR factorization is used first to reduce
!!>               the matrix to a triangular form.)
!!>          = 7: the number of processors
!!>          = 8: the crossover point for the multishift QR method
!!>               for nonsymmetric eigenvalue problems (DEPRECATED)
!!>          = 9: maximum size of the subproblems at the bottom of the
!!>               computation tree in the divide-and-conquer algorithm
!!>               (used by xGELSD and xGESDD)
!!>          =10: ieee NaN arithmetic can be trusted not to trap
!!>          =11: infinity arithmetic can be trusted not to trap
!!>          12 <= ISPEC <= 16:
!!>               xHSEQR or one of its subroutines,
!!>               see IPARMQ for detailed explanation
!!> \endverbatim
!!>
!!> \param[in] NAME
!!> \verbatim
!!>          NAME is CHARACTER*(*)
!!>          The name of the calling subroutine, in either upper case or
!!>          lower case.
!!> \endverbatim
!!>
!!> \param[in] OPTS
!!> \verbatim
!!>          OPTS is CHARACTER*(*)
!!>          The character options to the subroutine NAME, concatenated
!!>          into a single character string.  For example, UPLO = 'U',
!!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!!>          be specified as OPTS = 'UTN'.
!!> \endverbatim
!!>
!!> \param[in] N1
!!> \verbatim
!!>          N1 is INTEGER
!!> \endverbatim
!!>
!!> \param[in] N2
!!> \verbatim
!!>          N2 is INTEGER
!!> \endverbatim
!!>
!!> \param[in] N3
!!> \verbatim
!!>          N3 is INTEGER
!!> \endverbatim
!!>
!!> \param[in] N4
!!> \verbatim
!!>          N4 is INTEGER
!!>          Problem dimensions for the subroutine NAME; these may not all
!!>          be required.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup auxOTHERauxiliary
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  The following conventions have been used when calling ILAENV from the
!!>  LAPACK routines:
!!>  1)  OPTS is a concatenation of all of the character options to
!!>      subroutine NAME, in the same order that they appear in the
!!>      argument list for NAME, even if they are not used in determining
!!>      the value of the parameter specified by ISPEC.
!!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!!>      that they appear in the argument list for NAME.  N1 is used
!!>      first, N2 second, and so on, and unused problem dimensions are
!!>      passed a value of -1.
!!>  3)  The parameter value returned by ILAENV is checked for validity in
!!>      the calling subroutine.  For example, ILAENV is used to retrieve
!!>      the optimal blocksize for STRTRI as follows:
!!>
!!>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!!> \endverbatim
!!>
!!  =====================================================================
!      FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!!
!!  -- LAPACK auxiliary routine (version 3.4.0) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      CHARACTER*( * )    NAME, OPTS
!      INTEGER            ISPEC, N1, N2, N3, N4
!!     ..
!!
!!  =====================================================================
!!
!!     .. Local Scalars ..
!      INTEGER            I, IC, IZ, NB !, NBMIN, NX
!      LOGICAL            CNAME, SNAME
!      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!!     ..
!!     .. External Functions ..
!!      INTEGER            IEEECK !, IPARMQ
!!      EXTERNAL           IEEECK, IPARMQ
!!
!!
!      if (ISPEC > 16) then
!        print *,'inputs:',ISPEC, NAME, OPTS, N1, N2, N3, N4
!      endif
!      if (ILAENV.NE.1) then
!        ILAENV = -1
!        RETURN
!      endif
!!
!!     Convert NAME to upper case if the first character is lower case.
!!
!      ILAENV = 1
!      SUBNAM = NAME
!      IC = ICHAR( SUBNAM( 1: 1 ) )
!      IZ = ICHAR( 'Z' )
!      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!!
!!        ASCII character set
!!
!         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
!            SUBNAM( 1: 1 ) = CHAR( IC-32 )
!            DO 20 I = 2, 6
!               IC = ICHAR( SUBNAM( I: I ) )
!               IF( IC.GE.97 .AND. IC.LE.122 ) SUBNAM( I: I ) = CHAR( IC-32 )
!   20       CONTINUE
!         END IF
!!
!      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!!
!!        EBCDIC character set
!!
!         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
!              ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
!              ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
!            SUBNAM( 1: 1 ) = CHAR( IC+64 )
!            DO 30 I = 2, 6
!               IC = ICHAR( SUBNAM( I: I ) )
!               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
!                    ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
!                    ( IC.GE.162 .AND. IC.LE.169 ) ) &
!                    SUBNAM( I:I ) = CHAR( IC+64 )
!   30       CONTINUE
!         END IF
!!
!      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!!
!!        Prime machines:  ASCII+128
!!
!         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
!            SUBNAM( 1: 1 ) = CHAR( IC-32 )
!            DO 40 I = 2, 6
!               IC = ICHAR( SUBNAM( I: I ) )
!               IF( IC.GE.225 .AND. IC.LE.250 ) &
!                    SUBNAM( I: I ) = CHAR( IC-32 )
!   40       CONTINUE
!         END IF
!      END IF
!!
!      C1 = SUBNAM( 1: 1 )
!      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
!      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
!      IF( .NOT.( CNAME .OR. SNAME ) )   RETURN
!      C2 = SUBNAM( 2: 3 )
!      C3 = SUBNAM( 4: 6 )
!      C4 = C3( 2: 3 )
!!
!!
!!     ISPEC = 1:  block size
!!
!!     In these examples, separate code is provided for setting NB for
!!     real and complex.  We assume that NB will take the same value in
!!     single or double precision.
!!
!      NB = 1
!!
!      IF( C2.EQ.'GE' ) THEN
!         IF( C3.EQ.'TRF' ) THEN
!            IF( SNAME ) THEN
!               NB = 64
!            ELSE
!               NB = 64
!            END IF
!         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
!                 C3.EQ.'QLF' ) THEN
!            IF( SNAME ) THEN
!               NB = 32
!            ELSE
!               NB = 32
!            END IF
!         ELSE IF( C3.EQ.'HRD' ) THEN
!            IF( SNAME ) THEN
!               NB = 32
!            ELSE
!               NB = 32
!            END IF
!         ELSE IF( C3.EQ.'BRD' ) THEN
!            IF( SNAME ) THEN
!               NB = 32
!            ELSE
!               NB = 32
!            END IF
!         ELSE IF( C3.EQ.'TRI' ) THEN
!            IF( SNAME ) THEN
!               NB = 64
!            ELSE
!               NB = 64
!            END IF
!         END IF
!      ELSE IF( C2.EQ.'PO' ) THEN
!         IF( C3.EQ.'TRF' ) THEN
!            IF( SNAME ) THEN
!               NB = 64
!            ELSE
!               NB = 64
!            END IF
!         END IF
!      ELSE IF( C2.EQ.'SY' ) THEN
!         IF( C3.EQ.'TRF' ) THEN
!            IF( SNAME ) THEN
!               NB = 64
!            ELSE
!               NB = 64
!            END IF
!         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
!            NB = 32
!         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
!            NB = 64
!         END IF
!      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
!         IF( C3.EQ.'TRF' ) THEN
!            NB = 64
!         ELSE IF( C3.EQ.'TRD' ) THEN
!            NB = 32
!         ELSE IF( C3.EQ.'GST' ) THEN
!            NB = 64
!         END IF
!      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
!         IF( C3( 1: 1 ).EQ.'G' ) THEN
!            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
!               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
!                THEN
!               NB = 32
!            END IF
!         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
!            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
!               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
!                THEN
!               NB = 32
!            END IF
!         END IF
!      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
!         IF( C3( 1: 1 ).EQ.'G' ) THEN
!            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
!               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
!                THEN
!               NB = 32
!            END IF
!         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
!            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
!               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
!                THEN
!               NB = 32
!            END IF
!         END IF
!      ELSE IF( C2.EQ.'GB' ) THEN
!         IF( C3.EQ.'TRF' ) THEN
!            IF( SNAME ) THEN
!               IF( N4.LE.64 ) THEN
!                  NB = 1
!               ELSE
!                  NB = 32
!               END IF
!            ELSE
!               IF( N4.LE.64 ) THEN
!                  NB = 1
!               ELSE
!                  NB = 32
!               END IF
!            END IF
!         END IF
!      ELSE IF( C2.EQ.'PB' ) THEN
!         IF( C3.EQ.'TRF' ) THEN
!            IF( SNAME ) THEN
!               IF( N2.LE.64 ) THEN
!                  NB = 1
!               ELSE
!                  NB = 32
!               END IF
!            ELSE
!               IF( N2.LE.64 ) THEN
!                  NB = 1
!               ELSE
!                  NB = 32
!               END IF
!            END IF
!         END IF
!      ELSE IF( C2.EQ.'TR' ) THEN
!         IF( C3.EQ.'TRI' ) THEN
!            IF( SNAME ) THEN
!               NB = 64
!            ELSE
!               NB = 64
!            END IF
!         END IF
!      ELSE IF( C2.EQ.'LA' ) THEN
!         IF( C3.EQ.'UUM' ) THEN
!            IF( SNAME ) THEN
!               NB = 64
!            ELSE
!               NB = 64
!            END IF
!         END IF
!      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
!         IF( C3.EQ.'EBZ' ) THEN
!            NB = 1
!         END IF
!      END IF
!      ILAENV = NB
!      RETURN
!!
!
!!
!!     End of ILAENV
!!
!      END function ilaenv
!
!
!!  =========== DOCUMENTATION ===========
!!
!! Online html documentation available at
!!            http://www.netlib.org/lapack/explore-html/
!!
!!  Definition:
!!  ===========
!!
!!      LOGICAL FUNCTION LSAME( CA, CB )
!!
!!     .. Scalar Arguments ..
!!      CHARACTER          CA, CB
!!     ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!!> case.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] CA
!!> \verbatim
!!> \endverbatim
!!>
!!> \param[in] CB
!!> \verbatim
!!>          CA and CB specify the single characters to be compared.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup auxOTHERauxiliary
!!
!!  =====================================================================
!      LOGICAL FUNCTION LSAME( CA, CB )
!!
!!  -- LAPACK auxiliary routine (version 3.4.0) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      CHARACTER          CA, CB
!!     ..
!!
!! =====================================================================
!!
!!     .. Intrinsic Functions ..
!      INTRINSIC          ICHAR
!!     ..
!!     .. Local Scalars ..
!      INTEGER            INTA, INTB, ZCODE
!!     ..
!!     .. Executable Statements ..
!!
!!     Test if the characters are equal
!!
!      LSAME = CA.EQ.CB
!      IF( LSAME )   RETURN
!!
!!     Now test for equivalence if both characters are alphabetic.
!!
!      ZCODE = ICHAR( 'Z' )
!!
!!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!!     machines, on which ICHAR returns a value with bit 8 set.
!!     ICHAR('A') on Prime machines returns 193 which is the same as
!!     ICHAR('A') on an EBCDIC machine.
!!
!      INTA = ICHAR( CA )
!      INTB = ICHAR( CB )
!!
!      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!!
!!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!!        upper case 'Z'.
!!
!         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
!         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!!
!      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!!
!!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!!        upper case 'Z'.
!!
!         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
!            INTA.GE.145 .AND. INTA.LE.153 .OR. &
!            INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
!         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
!            INTB.GE.145 .AND. INTB.LE.153 .OR. &
!            INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!!
!      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!!
!!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!!        plus 128 of either lower or upper case 'Z'.
!!
!         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
!         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
!      END IF
!      LSAME = INTA.EQ.INTB
!!
!!     RETURN
!!
!!     End of LSAME
!!
!      END function lsame
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!!
!!       .. Scalar Arguments ..
!!       INTEGER            INFO, LDA, M, N
!!       ..
!!       .. Array Arguments ..
!!       INTEGER            IPIV( * )
!!       DOUBLE PRECISION   A( LDA, * )
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DGETRF computes an LU factorization of a general M-by-N matrix A
!!> using partial pivoting with row interchanges.
!!>
!!> The factorization has the form
!!>    A = P * L * U
!!> where P is a permutation matrix, L is lower triangular with unit
!!> diagonal elements (lower trapezoidal if m > n), and U is upper
!!> triangular (upper trapezoidal if m < n).
!!>
!!> This is the right-looking Level 3 BLAS version of the algorithm.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] M
!!> \verbatim
!!>          M is INTEGER
!!>          The number of rows of the matrix A.  M >= 0.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>          The number of columns of the matrix A.  N >= 0.
!!> \endverbatim
!!>
!!> \param[in,out] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!>          On entry, the M-by-N matrix to be factored.
!!>          On exit, the factors L and U from the factorization
!!>          A = P*L*U; the unit diagonal elements of L are not stored.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>          The leading dimension of the array A.  LDA >= max(1,M).
!!> \endverbatim
!!>
!!> \param[out] IPIV
!!> \verbatim
!!>          IPIV is INTEGER array, dimension (min(M,N))
!!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!!>          matrix was interchanged with row IPIV(i).
!!> \endverbatim
!!>
!!> \param[out] INFO
!!> \verbatim
!!>          INFO is INTEGER
!!>          = 0:  successful exit
!!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!!>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!!>                has been completed, but the factor U is exactly
!!>                singular, and division by zero will occur if it is used
!!>                to solve a system of equations.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup doubleGEcomputational
!!
!!  =====================================================================
!      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!!
!!  -- LAPACK computational routine (version 3.4.0) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      INTEGER            INFO, LDA, M, N
!!     ..
!!     .. Array Arguments ..
!      INTEGER            IPIV( * )
!      DOUBLE PRECISION   A( LDA, * )
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION   ONE
!      PARAMETER          ( ONE = 1.0D+0 )
!!     ..
!!     .. Local Scalars ..
!      INTEGER            I, IINFO, J, JB, NB, errstat
!!     ..
!!     .. External Subroutines ..
!!     EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM!, XERBLA
!!     ..
!!     .. External Functions ..
!!     INTEGER            ILAENV
!!     EXTERNAL           ILAENV
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!!     ..
!!     .. Executable Statements ..
!!
!!     Test the input parameters.
!!
!      INFO = 0
!      IF( M.LT.0 ) THEN
!         INFO = -1
!      ELSE IF( N.LT.0 ) THEN
!         INFO = -2
!      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
!         INFO = -4
!      END IF
!      IF( INFO.NE.0 ) THEN
!!         CALL XERBLA( 'DGETRF', -INFO )
!        errstat=INFO
!        call tell_error(tell_invalid_parm,'DGETRF', errstat)
!        RETURN
!      END IF
!!
!!     Quick return if possible
!!
!      IF( M.EQ.0 .OR. N.EQ.0 )   RETURN
!!
!!     Determine the block size for this environment.
!!
!      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
!      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!!
!!        Use unblocked code.
!!
!         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
!      ELSE
!!
!!        Use blocked code.
!!
!         DO 20 J = 1, MIN( M, N ), NB
!            JB = MIN( MIN( M, N )-J+1, NB )
!!
!!           Factor diagonal and subdiagonal blocks and test for exact
!!           singularity.
!!
!            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!!
!!           Adjust INFO and the pivot indices.
!!
!            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) &
!                 INFO = IINFO + J - 1
!            DO 10 I = J, MIN( M, J+JB-1 )
!               IPIV( I ) = J - 1 + IPIV( I )
!   10       CONTINUE
!!
!!           Apply interchanges to columns 1:J-1.
!!
!            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!!
!            IF( J+JB.LE.N ) THEN
!!
!!              Apply interchanges to columns J+JB:N.
!!
!               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, IPIV, 1 )
!!
!!              Compute block row of U.
!!
!               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
!                          N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
!                          LDA )
!               IF( J+JB.LE.M ) THEN
!!
!!                 Update trailing submatrix.
!!
!                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
!                             N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
!                             A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
!                             LDA )
!               END IF
!            END IF
!   20    CONTINUE
!      END IF
!      RETURN
!!
!!     End of DGETRF
!!
!      END subroutine dgetrf
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!!
!!       .. Scalar Arguments ..
!!       INTEGER            INFO, LDA, M, N
!!       ..
!!       .. Array Arguments ..
!!       INTEGER            IPIV( * )
!!       DOUBLE PRECISION   A( LDA, * )
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DGETF2 computes an LU factorization of a general m-by-n matrix A
!!> using partial pivoting with row interchanges.
!!>
!!> The factorization has the form
!!>    A = P * L * U
!!> where P is a permutation matrix, L is lower triangular with unit
!!> diagonal elements (lower trapezoidal if m > n), and U is upper
!!> triangular (upper trapezoidal if m < n).
!!>
!!> This is the right-looking Level 2 BLAS version of the algorithm.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] M
!!> \verbatim
!!>          M is INTEGER
!!>          The number of rows of the matrix A.  M >= 0.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>          The number of columns of the matrix A.  N >= 0.
!!> \endverbatim
!!>
!!> \param[in,out] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!>          On entry, the m by n matrix to be factored.
!!>          On exit, the factors L and U from the factorization
!!>          A = P*L*U; the unit diagonal elements of L are not stored.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>          The leading dimension of the array A.  LDA >= max(1,M).
!!> \endverbatim
!!>
!!> \param[out] IPIV
!!> \verbatim
!!>          IPIV is INTEGER array, dimension (min(M,N))
!!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!!>          matrix was interchanged with row IPIV(i).
!!> \endverbatim
!!>
!!> \param[out] INFO
!!> \verbatim
!!>          INFO is INTEGER
!!>          = 0: successful exit
!!>          < 0: if INFO = -k, the k-th argument had an illegal value
!!>          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!!>               has been completed, but the factor U is exactly
!!>               singular, and division by zero will occur if it is used
!!>               to solve a system of equations.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date September 2012
!!
!!> \ingroup doubleGEcomputational
!!
!!  =====================================================================
!      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!!
!!  -- LAPACK computational routine (version 3.4.2) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     September 2012
!!
!!     .. Scalar Arguments ..
!      INTEGER            INFO, LDA, M, N
!!     ..
!!     .. Array Arguments ..
!      INTEGER            IPIV( * )
!      DOUBLE PRECISION   A( LDA, * )
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION   ONE, ZERO
!      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION   SFMIN
!      INTEGER            I, J, JP, errstat
!!     ..
!!     .. External Functions ..
!!     DOUBLE PRECISION   DLAMCH
!!     INTEGER            IDAMAX
!!     EXTERNAL           DLAMCH, IDAMAX
!!     ..
!!     .. External Subroutines ..
!!     EXTERNAL           DGER, DSCAL, DSWAP !, XERBLA
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          MAX, MIN
!!     ..
!!     .. Executable Statements ..
!!
!!     Test the input parameters.
!!
!      INFO = 0
!      IF( M.LT.0 ) THEN
!         INFO = -1
!      ELSE IF( N.LT.0 ) THEN
!         INFO = -2
!      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
!         INFO = -4
!      END IF
!      IF( INFO.NE.0 ) THEN
!!         CALL XERBLA( 'DGETF2', -INFO )
!         errstat=INFO
!         call tell_error(tell_invalid_parm,'DGETF2', errstat)
!         RETURN
!      END IF
!!
!!     Quick return if possible
!!
!      IF( M.EQ.0 .OR. N.EQ.0 )   RETURN
!!
!!     Compute machine safe minimum
!!
!      SFMIN = DLAMCH('S')
!!
!      DO 10 J = 1, MIN( M, N )
!!
!!        Find pivot and test for singularity.
!!
!         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
!         IPIV( J ) = JP
!         IF( A( JP, J ).NE.ZERO ) THEN
!!
!!           Apply the interchange to columns 1:N.
!!
!            IF( JP.NE.J ) CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!!
!!           Compute elements J+1:M of J-th column.
!!
!            IF( J.LT.M ) THEN
!               IF( ABS(A( J, J )) .GE. SFMIN ) THEN
!                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!               ELSE
!                 DO 20 I = 1, M-J
!                    A( J+I, J ) = A( J+I, J ) / A( J, J )
!   20            CONTINUE
!               END IF
!            END IF
!!
!         ELSE IF( INFO.EQ.0 ) THEN
!!
!            INFO = J
!         END IF
!!
!         IF( J.LT.MIN( M, N ) ) THEN
!!
!!           Update trailing submatrix.
!!
!            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
!                 A( J+1, J+1 ), LDA )
!         END IF
!   10 CONTINUE
!      RETURN
!!
!!     End of DGETF2
!!
!      END subroutine dgetf2
!
!
!
!!  Definition:
!!  ===========
!!
!!      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DLAMCH determines double precision machine parameters.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] CMACH
!!> \verbatim
!!>          Specifies the value to be returned by DLAMCH:
!!>          = 'E' or 'e',   DLAMCH := eps
!!>          = 'S' or 's' ,   DLAMCH := sfmin
!!>          = 'B' or 'b',   DLAMCH := base
!!>          = 'P' or 'p',   DLAMCH := eps*base
!!>          = 'N' or 'n',   DLAMCH := t
!!>          = 'R' or 'r',   DLAMCH := rnd
!!>          = 'M' or 'm',   DLAMCH := emin
!!>          = 'U' or 'u',   DLAMCH := rmin
!!>          = 'L' or 'l',   DLAMCH := emax
!!>          = 'O' or 'o',   DLAMCH := rmax
!!>          where
!!>          eps   = relative machine precision
!!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!!>          base  = base of the machine
!!>          prec  = eps*base
!!>          t     = number of (base) digits in the mantissa
!!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!!>          emin  = minimum exponent before (gradual) underflow
!!>          rmin  = underflow threshold - base**(emin-1)
!!>          emax  = largest exponent before overflow
!!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup auxOTHERauxiliary
!!
!!  =====================================================================
!      FUNCTION DLAMCH( CMACH )
!!
!!  -- LAPACK auxiliary routine (version 3.4.0) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      CHARACTER          CMACH
!!     ..
!!
!!     .. Scalar Arguments ..
!!      DOUBLE PRECISION   A, B
!!     ..
!!
!! =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION   ONE, ZERO
!      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH, DLAMCH
!!     ..
!!     .. External Functions ..
!!     LOGICAL            LSAME
!!     EXTERNAL           LSAME
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, &
!                         MINEXPONENT, RADIX, TINY
!!     ..
!!     .. Executable Statements ..
!!
!!
!!     Assume rounding, not chopping. Always.
!!
!      RND = ONE
!!
!      IF( ONE.EQ.RND ) THEN
!         EPS = EPSILON(ZERO) * 0.5
!      ELSE
!         EPS = EPSILON(ZERO)
!      END IF
!!
!      IF( LSAME( CMACH, 'E' ) ) THEN
!         RMACH = EPS
!      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
!         SFMIN = TINY(ZERO)
!         SMALL = ONE / HUGE(ZERO)
!         IF( SMALL.GE.SFMIN ) THEN
!!
!!           Use SMALL plus a bit, to avoid the possibility of rounding
!!           causing overflow when computing  1/sfmin.
!!
!            SFMIN = SMALL*( ONE+EPS )
!         END IF
!         RMACH = SFMIN
!      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
!         RMACH = RADIX(ZERO)
!      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
!         RMACH = EPS * RADIX(ZERO)
!      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
!         RMACH = DIGITS(ZERO)
!      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
!         RMACH = RND
!      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
!         RMACH = MINEXPONENT(ZERO)
!      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
!         RMACH = tiny(zero)
!      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
!         RMACH = MAXEXPONENT(ZERO)
!      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
!         RMACH = HUGE(ZERO)
!      ELSE
!         RMACH = ZERO
!      END IF
!!
!      DLAMCH = RMACH
!      RETURN
!!
!!     End of DLAMCH
!!
!    end function dlamch
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!!
!!       .. Scalar Arguments ..
!!       INTEGER            INCX, K1, K2, LDA, N
!!       ..
!!       .. Array Arguments ..
!!       INTEGER            IPIV( * )
!!       DOUBLE PRECISION   A( LDA, * )
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DLASWP performs a series of row interchanges on the matrix A.
!!> One row interchange is initiated for each of rows K1 through K2 of A.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>          The number of columns of the matrix A.
!!> \endverbatim
!!>
!!> \param[in,out] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!>          On entry, the matrix of column dimension N to which the row
!!>          interchanges will be applied.
!!>          On exit, the permuted matrix.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>          The leading dimension of the array A.
!!> \endverbatim
!!>
!!> \param[in] K1
!!> \verbatim
!!>          K1 is INTEGER
!!>          The first element of IPIV for which a row interchange will
!!>          be done.
!!> \endverbatim
!!>
!!> \param[in] K2
!!> \verbatim
!!>          K2 is INTEGER
!!>          The last element of IPIV for which a row interchange will
!!>          be done.
!!> \endverbatim
!!>
!!> \param[in] IPIV
!!> \verbatim
!!>          IPIV is INTEGER array, dimension (K2*abs(INCX))
!!>          The vector of pivot indices.  Only the elements in positions
!!>          K1 through K2 of IPIV are accessed.
!!>          IPIV(K) = L implies rows K and L are to be interchanged.
!!> \endverbatim
!!>
!!> \param[in] INCX
!!> \verbatim
!!>          INCX is INTEGER
!!>          The increment between successive values of IPIV.  If IPIV
!!>          is negative, the pivots are applied in reverse order.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date September 2012
!!
!!> \ingroup doubleOTHERauxiliary
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  Modified by
!!>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!!
!!  -- LAPACK auxiliary routine (version 3.4.2) --
!!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     September 2012
!!
!!     .. Scalar Arguments ..
!      INTEGER            INCX, K1, K2, LDA, N
!!     ..
!!     .. Array Arguments ..
!      INTEGER            IPIV( * )
!      DOUBLE PRECISION   A( LDA, * )
!!     ..
!!
!! =====================================================================
!!
!!     .. Local Scalars ..
!      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
!      DOUBLE PRECISION   TEMP
!!     ..
!!     .. Executable Statements ..
!!
!!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!!
!      IF( INCX.GT.0 ) THEN
!         IX0 = K1
!         I1 = K1
!         I2 = K2
!         INC = 1
!      ELSE IF( INCX.LT.0 ) THEN
!         IX0 = 1 + ( 1-K2 )*INCX
!         I1 = K2
!         I2 = K1
!         INC = -1
!      ELSE
!         RETURN
!      END IF
!!
!      N32 = ( N / 32 )*32
!      IF( N32.NE.0 ) THEN
!         DO 30 J = 1, N32, 32
!            IX = IX0
!            DO 20 I = I1, I2, INC
!               IP = IPIV( IX )
!               IF( IP.NE.I ) THEN
!                  DO 10 K = J, J + 31
!                     TEMP = A( I, K )
!                     A( I, K ) = A( IP, K )
!                     A( IP, K ) = TEMP
!   10             CONTINUE
!               END IF
!               IX = IX + INCX
!   20       CONTINUE
!   30    CONTINUE
!      END IF
!      IF( N32.NE.N ) THEN
!         N32 = N32 + 1
!         IX = IX0
!         DO 50 I = I1, I2, INC
!            IP = IPIV( IX )
!            IF( IP.NE.I ) THEN
!               DO 40 K = N32, N
!                  TEMP = A( I, K )
!                  A( I, K ) = A( IP, K )
!                  A( IP, K ) = TEMP
!   40          CONTINUE
!            END IF
!            IX = IX + INCX
!   50    CONTINUE
!      END IF
!!
!      RETURN
!!
!!     End of DLASWP
!!
!      END subroutine dlaswp
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!!
!!       .. Scalar Arguments ..
!!       DOUBLE PRECISION ALPHA
!!       INTEGER LDA,LDB,M,N
!!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION A(LDA,*),B(LDB,*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DTRSM  solves one of the matrix equations
!!>
!!>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!!>
!!> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!!>
!!>    op( A ) = A   or   op( A ) = A**T.
!!>
!!> The matrix X is overwritten on B.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] SIDE
!!> \verbatim
!!>          SIDE is CHARACTER*1
!!>           On entry, SIDE specifies whether op( A ) appears on the left
!!>           or right of X as follows:
!!>
!!>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!!>
!!>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!!> \endverbatim
!!>
!!> \param[in] UPLO
!!> \verbatim
!!>          UPLO is CHARACTER*1
!!>           On entry, UPLO specifies whether the matrix A is an upper or
!!>           lower triangular matrix as follows:
!!>
!!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!!>
!!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!!> \endverbatim
!!>
!!> \param[in] TRANSA
!!> \verbatim
!!>          TRANSA is CHARACTER*1
!!>           On entry, TRANSA specifies the form of op( A ) to be used in
!!>           the matrix multiplication as follows:
!!>
!!>              TRANSA = 'N' or 'n'   op( A ) = A.
!!>
!!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!!>
!!>              TRANSA = 'C' or 'c'   op( A ) = A**T.
!!> \endverbatim
!!>
!!> \param[in] DIAG
!!> \verbatim
!!>          DIAG is CHARACTER*1
!!>           On entry, DIAG specifies whether or not A is unit triangular
!!>           as follows:
!!>
!!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!!>
!!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!!>                                  triangular.
!!> \endverbatim
!!>
!!> \param[in] M
!!> \verbatim
!!>          M is INTEGER
!!>           On entry, M specifies the number of rows of B. M must be at
!!>           least zero.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>           On entry, N specifies the number of columns of B.  N must be
!!>           at least zero.
!!> \endverbatim
!!>
!!> \param[in] ALPHA
!!> \verbatim
!!>          ALPHA is DOUBLE PRECISION.
!!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!!>           zero then  A is not referenced and  B need not be set before
!!>           entry.
!!> \endverbatim
!!>
!!> \param[in] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, k ),
!!>           where k is m when SIDE = 'L' or 'l'
!!>             and k is n when SIDE = 'R' or 'r'.
!!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!!>           upper triangular part of the array  A must contain the upper
!!>           triangular matrix  and the strictly lower triangular part of
!!>           A is not referenced.
!!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!!>           lower triangular part of the array  A must contain the lower
!!>           triangular matrix  and the strictly upper triangular part of
!!>           A is not referenced.
!!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!!>           A  are not referenced either,  but are assumed to be  unity.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>           On entry, LDA specifies the first dimension of A as declared
!!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!!>           then LDA must be at least max( 1, n ).
!!> \endverbatim
!!>
!!> \param[in,out] B
!!> \verbatim
!!>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!!>           Before entry,  the leading  m by n part of the array  B must
!!>           contain  the  right-hand  side  matrix  B,  and  on exit  is
!!>           overwritten by the solution matrix  X.
!!> \endverbatim
!!>
!!> \param[in] LDB
!!> \verbatim
!!>          LDB is INTEGER
!!>           On entry, LDB specifies the first dimension of B as declared
!!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!!>           max( 1, m ).
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level3
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  Level 3 Blas routine.
!!>
!!>
!!>  -- Written on 8-February-1989.
!!>     Jack Dongarra, Argonne National Laboratory.
!!>     Iain Duff, AERE Harwell.
!!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!!> \endverbatim
!!>
!!  =====================================================================
!      SUBROUTINE dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!!
!!  -- Reference BLAS level3 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      DOUBLE PRECISION alpha
!      INTEGER lda,ldb,m,n
!      CHARACTER diag,side,transa,uplo
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION a(lda,*),b(ldb,*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. External Functions ..
!!     LOGICAL lsame
!!     EXTERNAL lsame
!!     ..
!!     .. External Subroutines ..
!!      EXTERNAL xerbla
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC max
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION temp
!      INTEGER i,info,j,k,nrowa
!      LOGICAL lside,nounit,upper
!!     ..
!!     .. Parameters ..
!      DOUBLE PRECISION one,zero
!      parameter(one=1.0d+0,zero=0.0d+0)
!!     ..
!!
!!     Test the input parameters.
!!
!      lside = lsame(side,'L')
!      IF (lside) THEN
!          nrowa = m
!      ELSE
!          nrowa = n
!      END IF
!      nounit = lsame(diag,'N')
!      upper = lsame(uplo,'U')
!!
!      info = 0
!      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
!          info = 1
!      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
!          info = 2
!      ELSE IF ((.NOT.lsame(transa,'N')) .AND. &
!               (.NOT.lsame(transa,'T')) .AND. &
!               (.NOT.lsame(transa,'C'))) THEN
!          info = 3
!      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
!          info = 4
!      ELSE IF (m.LT.0) THEN
!          info = 5
!      ELSE IF (n.LT.0) THEN
!          info = 6
!      ELSE IF (lda.LT.max(1,nrowa)) THEN
!          info = 9
!      ELSE IF (ldb.LT.max(1,m)) THEN
!          info = 11
!      END IF
!      IF (info.NE.0) THEN
!!          !CALL xerbla('DTRSM ',info)
!          call tell_error(tell_invalid_parm,'DTRSM ',info)
!          RETURN
!      END IF
!!
!!     Quick return if possible.
!!
!      IF (m.EQ.0 .OR. n.EQ.0) RETURN
!!
!!     And when  alpha.eq.zero.
!!
!      IF (alpha.EQ.zero) THEN
!          DO 20 j = 1,n
!              DO 10 i = 1,m
!                  b(i,j) = zero
!   10         CONTINUE
!   20     CONTINUE
!          RETURN
!      END IF
!!
!!     Start the operations.
!!
!      IF (lside) THEN
!          IF (lsame(transa,'N')) THEN
!!
!!           Form  B := alpha*inv( A )*B.
!!
!              IF (upper) THEN
!                  DO 60 j = 1,n
!                      IF (alpha.NE.one) THEN
!                          DO 30 i = 1,m
!                              b(i,j) = alpha*b(i,j)
!   30                     CONTINUE
!                      END IF
!                      DO 50 k = m,1,-1
!                          IF (b(k,j).NE.zero) THEN
!                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
!                              DO 40 i = 1,k - 1
!                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
!   40                         CONTINUE
!                          END IF
!   50                 CONTINUE
!   60             CONTINUE
!              ELSE
!                  DO 100 j = 1,n
!                      IF (alpha.NE.one) THEN
!                          DO 70 i = 1,m
!                              b(i,j) = alpha*b(i,j)
!   70                     CONTINUE
!                      END IF
!                      DO 90 k = 1,m
!                          IF (b(k,j).NE.zero) THEN
!                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
!                              DO 80 i = k + 1,m
!                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
!   80                         CONTINUE
!                          END IF
!   90                 CONTINUE
!  100             CONTINUE
!              END IF
!          ELSE
!!
!!           Form  B := alpha*inv( A**T )*B.
!!
!              IF (upper) THEN
!                  DO 130 j = 1,n
!                      DO 120 i = 1,m
!                          temp = alpha*b(i,j)
!                          DO 110 k = 1,i - 1
!                              temp = temp - a(k,i)*b(k,j)
!  110                     CONTINUE
!                          IF (nounit) temp = temp/a(i,i)
!                          b(i,j) = temp
!  120                 CONTINUE
!  130             CONTINUE
!              ELSE
!                  DO 160 j = 1,n
!                      DO 150 i = m,1,-1
!                          temp = alpha*b(i,j)
!                          DO 140 k = i + 1,m
!                              temp = temp - a(k,i)*b(k,j)
!  140                     CONTINUE
!                          IF (nounit) temp = temp/a(i,i)
!                          b(i,j) = temp
!  150                 CONTINUE
!  160             CONTINUE
!              END IF
!          END IF
!      ELSE
!          IF (lsame(transa,'N')) THEN
!!
!!           Form  B := alpha*B*inv( A ).
!!
!              IF (upper) THEN
!                  DO 210 j = 1,n
!                      IF (alpha.NE.one) THEN
!                          DO 170 i = 1,m
!                              b(i,j) = alpha*b(i,j)
!  170                     CONTINUE
!                      END IF
!                      DO 190 k = 1,j - 1
!                          IF (a(k,j).NE.zero) THEN
!                              DO 180 i = 1,m
!                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
!  180                         CONTINUE
!                          END IF
!  190                 CONTINUE
!                      IF (nounit) THEN
!                          temp = one/a(j,j)
!                          DO 200 i = 1,m
!                              b(i,j) = temp*b(i,j)
!  200                     CONTINUE
!                      END IF
!  210             CONTINUE
!              ELSE
!                  DO 260 j = n,1,-1
!                      IF (alpha.NE.one) THEN
!                          DO 220 i = 1,m
!                              b(i,j) = alpha*b(i,j)
!  220                     CONTINUE
!                      END IF
!                      DO 240 k = j + 1,n
!                          IF (a(k,j).NE.zero) THEN
!                              DO 230 i = 1,m
!                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
!  230                         CONTINUE
!                          END IF
!  240                 CONTINUE
!                      IF (nounit) THEN
!                          temp = one/a(j,j)
!                          DO 250 i = 1,m
!                              b(i,j) = temp*b(i,j)
!  250                     CONTINUE
!                      END IF
!  260             CONTINUE
!              END IF
!          ELSE
!!
!!           Form  B := alpha*B*inv( A**T ).
!!
!              IF (upper) THEN
!                  DO 310 k = n,1,-1
!                      IF (nounit) THEN
!                          temp = one/a(k,k)
!                          DO 270 i = 1,m
!                              b(i,k) = temp*b(i,k)
!  270                     CONTINUE
!                      END IF
!                      DO 290 j = 1,k - 1
!                          IF (a(j,k).NE.zero) THEN
!                              temp = a(j,k)
!                              DO 280 i = 1,m
!                                  b(i,j) = b(i,j) - temp*b(i,k)
!  280                         CONTINUE
!                          END IF
!  290                 CONTINUE
!                      IF (alpha.NE.one) THEN
!                          DO 300 i = 1,m
!                              b(i,k) = alpha*b(i,k)
!  300                     CONTINUE
!                      END IF
!  310             CONTINUE
!              ELSE
!                  DO 360 k = 1,n
!                      IF (nounit) THEN
!                          temp = one/a(k,k)
!                          DO 320 i = 1,m
!                              b(i,k) = temp*b(i,k)
!  320                     CONTINUE
!                      END IF
!                      DO 340 j = k + 1,n
!                          IF (a(j,k).NE.zero) THEN
!                              temp = a(j,k)
!                              DO 330 i = 1,m
!                                  b(i,j) = b(i,j) - temp*b(i,k)
!  330                         CONTINUE
!                          END IF
!  340                 CONTINUE
!                      IF (alpha.NE.one) THEN
!                          DO 350 i = 1,m
!                              b(i,k) = alpha*b(i,k)
!  350                     CONTINUE
!                      END IF
!  360             CONTINUE
!              END IF
!          END IF
!      END IF
!!
!      RETURN
!!
!!     End of DTRSM .
!!
!      END subroutine dtrsm
!
!
!
!!  Definition:
!!  ===========
!!
!!       INTEGER FUNCTION IDAMAX(N,DX,INCX)
!!
!!       .. Scalar Arguments ..
!!       INTEGER INCX,N
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION DX(*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!>    IDAMAX finds the index of element having max. absolute value.
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup aux_blas
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>     jack dongarra, linpack, 3/11/78.
!!>     modified 3/93 to return if incx .le. 0.
!!>     modified 12/3/93, array(1) declarations changed to array(*)
!!> \endverbatim
!!>
!!  =====================================================================
!      FUNCTION idamax(N,DX,INCX)
!!
!!  -- Reference BLAS level1 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      INTEGER incx,n
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION dx(*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. Local Scalars ..
!      DOUBLE PRECISION dmax
!      INTEGER i,ix
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC dabs
!!     ..
!      idamax = 0
!      IF (n.LT.1 .OR. incx.LE.0) RETURN
!      idamax = 1
!      IF (n.EQ.1) RETURN
!      IF (incx.EQ.1) THEN
!!
!!        code for increment equal to 1
!!
!         dmax = dabs(dx(1))
!         DO i = 2,n
!            IF (dabs(dx(i)).GT.dmax) THEN
!               idamax = i
!               dmax = dabs(dx(i))
!            END IF
!         END DO
!      ELSE
!!
!!        code for increment not equal to 1
!!
!         ix = 1
!         dmax = dabs(dx(1))
!         ix = ix + incx
!         DO i = 2,n
!            IF (dabs(dx(ix)).GT.dmax) THEN
!               idamax = i
!               dmax = dabs(dx(ix))
!            END IF
!            ix = ix + incx
!         END DO
!      END IF
!      RETURN
!      END function idamax
!
!
!
!!  Definition:
!!  ===========
!!
!!       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!!
!!       .. Scalar Arguments ..
!!       DOUBLE PRECISION ALPHA
!!       INTEGER INCX,INCY,LDA,M,N
!!       ..
!!       .. Array Arguments ..
!!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!!       ..
!!
!!
!!> \par Purpose:
!!  =============
!!>
!!> \verbatim
!!>
!!> DGER   performs the rank 1 operation
!!>
!!>    A := alpha*x*y**T + A,
!!>
!!> where alpha is a scalar, x is an m element vector, y is an n element
!!> vector and A is an m by n matrix.
!!> \endverbatim
!!
!!  Arguments:
!!  ==========
!!
!!> \param[in] M
!!> \verbatim
!!>          M is INTEGER
!!>           On entry, M specifies the number of rows of the matrix A.
!!>           M must be at least zero.
!!> \endverbatim
!!>
!!> \param[in] N
!!> \verbatim
!!>          N is INTEGER
!!>           On entry, N specifies the number of columns of the matrix A.
!!>           N must be at least zero.
!!> \endverbatim
!!>
!!> \param[in] ALPHA
!!> \verbatim
!!>          ALPHA is DOUBLE PRECISION.
!!>           On entry, ALPHA specifies the scalar alpha.
!!> \endverbatim
!!>
!!> \param[in] X
!!> \verbatim
!!>          X is DOUBLE PRECISION array of dimension at least
!!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!!>           Before entry, the incremented array X must contain the m
!!>           element vector x.
!!> \endverbatim
!!>
!!> \param[in] INCX
!!> \verbatim
!!>          INCX is INTEGER
!!>           On entry, INCX specifies the increment for the elements of
!!>           X. INCX must not be zero.
!!> \endverbatim
!!>
!!> \param[in] Y
!!> \verbatim
!!>          Y is DOUBLE PRECISION array of dimension at least
!!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!!>           Before entry, the incremented array Y must contain the n
!!>           element vector y.
!!> \endverbatim
!!>
!!> \param[in] INCY
!!> \verbatim
!!>          INCY is INTEGER
!!>           On entry, INCY specifies the increment for the elements of
!!>           Y. INCY must not be zero.
!!> \endverbatim
!!>
!!> \param[in,out] A
!!> \verbatim
!!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!!>           Before entry, the leading m by n part of the array A must
!!>           contain the matrix of coefficients. On exit, A is
!!>           overwritten by the updated matrix.
!!> \endverbatim
!!>
!!> \param[in] LDA
!!> \verbatim
!!>          LDA is INTEGER
!!>           On entry, LDA specifies the first dimension of A as declared
!!>           in the calling (sub) program. LDA must be at least
!!>           max( 1, m ).
!!> \endverbatim
!!
!!  Authors:
!!  ========
!!
!!> \author Univ. of Tennessee
!!> \author Univ. of California Berkeley
!!> \author Univ. of Colorado Denver
!!> \author NAG Ltd.
!!
!!> \date November 2011
!!
!!> \ingroup double_blas_level2
!!
!!> \par Further Details:
!!  =====================
!!>
!!> \verbatim
!!>
!!>  Level 2 Blas routine.
!!>
!!>  -- Written on 22-October-1986.
!!>     Jack Dongarra, Argonne National Lab.
!!>     Jeremy Du Croz, Nag Central Office.
!!>     Sven Hammarling, Nag Central Office.
!!>     Richard Hanson, Sandia National Labs.
!!> \endverbatim
!!>
!!  =====================================================================
!     SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!!
!!  -- Reference BLAS level2 routine (version 3.4.0) --
!!  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
!!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!!     November 2011
!!
!!     .. Scalar Arguments ..
!      DOUBLE PRECISION alpha
!      INTEGER incx,incy,lda,m,n
!!     ..
!!     .. Array Arguments ..
!      DOUBLE PRECISION a(lda,*),x(*),y(*)
!!     ..
!!
!!  =====================================================================
!!
!!     .. Parameters ..
!      DOUBLE PRECISION zero
!      parameter(zero=0.0d+0)
!!     ..
!!     .. Local Scalars ..
!      DOUBLE PRECISION temp
!      INTEGER i,info,ix,j,jy,kx
!!     ..
!!     .. External Subroutines ..
!!      EXTERNAL xerbla
!!     ..
!!     .. Intrinsic Functions ..
!      INTRINSIC max
!!     ..
!!
!!     Test the input parameters.
!!
!      info = 0
!      IF (m.LT.0) THEN
!          info = 1
!      ELSE IF (n.LT.0) THEN
!          info = 2
!      ELSE IF (incx.EQ.0) THEN
!          info = 5
!      ELSE IF (incy.EQ.0) THEN
!          info = 7
!      ELSE IF (lda.LT.max(1,m)) THEN
!          info = 9
!      END IF
!      IF (info.NE.0) THEN
!!          CALL xerbla('DGER  ',info)
!          call tell_error(tell_invalid_parm,'DGER  ',info)
!          RETURN
!      END IF
!!
!!     Quick return if possible.
!!
!      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
!!
!!     Start the operations. In this version the elements of A are
!!     accessed sequentially with one pass through A.
!!
!      IF (incy.GT.0) THEN
!          jy = 1
!      ELSE
!          jy = 1 - (n-1)*incy
!      END IF
!      IF (incx.EQ.1) THEN
!          DO 20 j = 1,n
!              IF (y(jy).NE.zero) THEN
!                  temp = alpha*y(jy)
!                  DO 10 i = 1,m
!                      a(i,j) = a(i,j) + x(i)*temp
!   10             CONTINUE
!              END IF
!              jy = jy + incy
!   20     CONTINUE
!      ELSE
!          IF (incx.GT.0) THEN
!              kx = 1
!          ELSE
!              kx = 1 - (m-1)*incx
!          END IF
!          DO 40 j = 1,n
!              IF (y(jy).NE.zero) THEN
!                  temp = alpha*y(jy)
!                  ix = kx
!                  DO 30 i = 1,m
!                      a(i,j) = a(i,j) + x(ix)*temp
!                      ix = ix + incx
!   30             CONTINUE
!              END IF
!              jy = jy + incy
!   40     CONTINUE
!      END IF
!!
!      RETURN
!!
!!     End of DGER  .
!!
!      END subroutine dger
!
!


end module m_invert2
