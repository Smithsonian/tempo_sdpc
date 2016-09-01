!> LAPACK routines for complex numbers
module LAPACK_Z16_tools

  use LAPACK_BLAS_tools
  use LAPACK_tools, only: lsame, ilaenv
  use m_xerbla


  !! Note many of the routines in this module appear to be unused!

  public zswap

  private zgetf2, zgetrf, zlaswp, zgetri, ztrti2, ztrtri, zgeru, ztrsm, &
       zgemm, ztrmv, ztrmm, zgemv, izamax, dcabs1

contains

  SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
    !
    !  -- LAPACK routine (version 3.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     September 30, 1994
    !
    !     .. Scalar Arguments ..
    INTEGER            INFO, LDA, M, N
    !     ..
    !     .. Array Arguments ..
    INTEGER            IPIV( * )
    COMPLEX*16         A( LDA, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZGETF2 computes an LU factorization of a general m-by-n matrix A
    !  using partial pivoting with row interchanges.
    !
    !  The factorization has the form
    !     A = P * L * U
    !  where P is a permutation matrix, L is lower triangular with unit
    !  diagonal elements (lower trapezoidal if m > n), and U is upper
    !  triangular (upper trapezoidal if m < n).
    !
    !  This is the right-looking Level 2 BLAS version of the algorithm.
    !
    !  Arguments
    !  =========
    !
    !  M       (input) INTEGER
    !          The number of rows of the matrix A.  M >= 0.
    !
    !  N       (input) INTEGER
    !          The number of columns of the matrix A.  N >= 0.
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          On entry, the m by n matrix to be factored.
    !          On exit, the factors L and U from the factorization
    !          A = P*L*U; the unit diagonal elements of L are not stored.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,M).
    !
    !  IPIV    (output) INTEGER array, dimension (min(M,N))
    !          The pivot indices; for 1 <= i <= min(M,N), row i of the
    !          matrix was interchanged with row IPIV(i).
    !
    !  INFO    (output) INTEGER
    !          = 0: successful exit
    !          < 0: if INFO = -k, the k-th argument had an illegal value
    !          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
    !               has been completed, but the factor U is exactly
    !               singular, and division by zero will occur if it is used
    !               to solve a system of equations.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    COMPLEX*16         ONE, ZERO
    PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
         ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Local Scalars ..
    INTEGER            J, JP
    !     ..
    !     .. External Functions ..
    !INTEGER            IZAMAX
    !EXTERNAL           IZAMAX
    !     ..
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC          MAX, MIN
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF( M.LT.0 ) THEN
      INFO = -1
    ELSE IF( N.LT.0 ) THEN
      INFO = -2
    ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
      INFO = -4
    END IF
    IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'ZGETF2', -INFO )
      RETURN
    END IF
    !
    !     Quick return if possible
    !
    IF( M.EQ.0 .OR. N.EQ.0 )   RETURN
    !
    DO J = 1, MIN( M, N )
      !
      !        Find pivot and test for singularity.
      !
      JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
      IPIV( J ) = JP
      IF( A( JP, J ).NE.ZERO ) THEN
        !
        !           Apply the interchange to columns 1:N.
        !
        IF( JP.NE.J ) &
             CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
        !
        !           Compute elements J+1:M of J-th column.
        !
        IF( J.LT.M ) &
             CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
        !
      ELSE IF( INFO.EQ.0 ) THEN
        !
        INFO = J
      END IF
      !
      IF( J.LT.MIN( M, N ) ) THEN
        !
        !           Update trailing submatrix.
        !
        CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), &
             LDA, A( J+1, J+1 ), LDA )
      END IF
    end do
    RETURN
    !
    !     End of ZGETF2
    !
  END SUBROUTINE ZGETF2



  SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
    !
    !  -- LAPACK routine (version 3.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     September 30, 1994
    !
    !     .. Scalar Arguments ..
    INTEGER            INFO, LDA, M, N
    !     ..
    !     .. Array Arguments ..
    INTEGER            IPIV( * )
    COMPLEX*16         A( LDA, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZGETRF computes an LU factorization of a general M-by-N matrix A
    !  using partial pivoting with row interchanges.
    !
    !  The factorization has the form
    !     A = P * L * U
    !  where P is a permutation matrix, L is lower triangular with unit
    !  diagonal elements (lower trapezoidal if m > n), and U is upper
    !  triangular (upper trapezoidal if m < n).
    !
    !  This is the right-looking Level 3 BLAS version of the algorithm.
    !
    !  Arguments
    !  =========
    !
    !  M       (input) INTEGER
    !          The number of rows of the matrix A.  M >= 0.
    !
    !  N       (input) INTEGER
    !          The number of columns of the matrix A.  N >= 0.
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          On entry, the M-by-N matrix to be factored.
    !          On exit, the factors L and U from the factorization
    !          A = P*L*U; the unit diagonal elements of L are not stored.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,M).
    !
    !  IPIV    (output) INTEGER array, dimension (min(M,N))
    !          The pivot indices; for 1 <= i <= min(M,N), row i of the
    !          matrix was interchanged with row IPIV(i).
    !
    !  INFO    (output) INTEGER
    !          = 0:  successful exit
    !          < 0:  if INFO = -i, the i-th argument had an illegal value
    !          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
    !                has been completed, but the factor U is exactly
    !                singular, and division by zero will occur if it is used
    !                to solve a system of equations.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    COMPLEX*16         ONE
    PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Local Scalars ..
    INTEGER            I, IINFO, J, JB, NB
    !     ..
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
    !     ..
    !     .. External Functions ..
    !INTEGER            ILAENV
    !EXTERNAL           ILAENV
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC          MAX, MIN
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF( M.LT.0 ) THEN
      INFO = -1
    ELSE IF( N.LT.0 ) THEN
      INFO = -2
    ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
      INFO = -4
    END IF
    IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'ZGETRF', -INFO )
      RETURN
    END IF
    !
    !     Quick return if possible
    !
    IF( M.EQ.0 .OR. N.EQ.0 )   RETURN
    !
    !     Determine the block size for this environment.
    !
    NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
    IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
      !
      !        Use unblocked code.
      !
      CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
    ELSE
      !
      !        Use blocked code.
      !
      DO J = 1, MIN( M, N ), NB
        JB = MIN( MIN( M, N )-J+1, NB )
        !
        !           Factor diagonal and subdiagonal blocks and test for exact
        !           singularity.
        !
        CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
        !
        !           Adjust INFO and the pivot indices.
        !
        IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
        DO  I = J, MIN( M, J+JB-1 )
          IPIV( I ) = J - 1 + IPIV( I )
        end do
        !
        !           Apply interchanges to columns 1:J-1.
        !
        CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
        !
        IF( J+JB.LE.N ) THEN
          !
          !              Apply interchanges to columns J+JB:N.
          !
          CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
               IPIV, 1 )
          !
          !              Compute block row of U.
          !
          CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
               N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
               LDA )
          IF( J+JB.LE.M ) THEN
            !
            !                 Update trailing submatrix.
            !
            CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                 N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                 A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                 LDA )
          END IF
        END IF
      end do
    END IF
    RETURN
    !
    !     End of ZGETRF
    !
  END SUBROUTINE ZGETRF



  SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
    !
    !  -- LAPACK auxiliary routine (version 3.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     June 30, 1999
    !
    !     .. Scalar Arguments ..
    INTEGER            INCX, K1, K2, LDA, N
    !     ..
    !     .. Array Arguments ..
    INTEGER            IPIV( * )
    COMPLEX*16         A( LDA, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZLASWP performs a series of row interchanges on the matrix A.
    !  One row interchange is initiated for each of rows K1 through K2 of A.
    !
    !  Arguments
    !  =========
    !
    !  N       (input) INTEGER
    !          The number of columns of the matrix A.
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          On entry, the matrix of column dimension N to which the row
    !          interchanges will be applied.
    !          On exit, the permuted matrix.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.
    !
    !  K1      (input) INTEGER
    !          The first element of IPIV for which a row interchange will
    !          be done.
    !
    !  K2      (input) INTEGER
    !          The last element of IPIV for which a row interchange will
    !          be done.
    !
    !  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
    !          The vector of pivot indices.  Only the elements in positions
    !          K1 through K2 of IPIV are accessed.
    !          IPIV(K) = L implies rows K and L are to be interchanged.
    !
    !  INCX    (input) INTEGER
    !          The increment between successive values of IPIV.  If IPIV
    !          is negative, the pivots are applied in reverse order.
    !
    !  Further Details
    !  ===============
    !
    !  Modified by
    !   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
    !
    ! =====================================================================
    !
    !     .. Local Scalars ..
    INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
    COMPLEX*16         TEMP
    !     ..
    !     .. Executable Statements ..
    !
    !     Interchange row I with row IPIV(I) for each of rows K1 through K2.
    !
    IF( INCX.GT.0 ) THEN
      IX0 = K1
      I1 = K1
      I2 = K2
      INC = 1
    ELSE IF( INCX.LT.0 ) THEN
      IX0 = 1 + ( 1-K2 )*INCX
      I1 = K2
      I2 = K1
      INC = -1
    ELSE
      RETURN
    END IF
    !
    N32 = ( N / 32 )*32
    IF( N32.NE.0 ) THEN
      DO J = 1, N32, 32
        IX = IX0
        DO I = I1, I2, INC
          IP = IPIV( IX )
          IF( IP.NE.I ) THEN
            DO K = J, J + 31
              TEMP = A( I, K )
              A( I, K ) = A( IP, K )
              A( IP, K ) = TEMP
            end do
          END IF
          IX = IX + INCX
        end do
      end do
    END IF
    IF( N32.NE.N ) THEN
      N32 = N32 + 1
      IX = IX0
      DO I = I1, I2, INC
        IP = IPIV( IX )
        IF( IP.NE.I ) THEN
          DO K = N32, N
            TEMP = A( I, K )
            A( I, K ) = A( IP, K )
            A( IP, K ) = TEMP
          end do
        END IF
        IX = IX + INCX
      end do
    END IF
    !
    RETURN
    !
    !     End of ZLASWP
    !
  END SUBROUTINE ZLASWP



  SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    !
    !  -- LAPACK routine (version 3.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     June 30, 1999
    !
    !     .. Scalar Arguments ..
    INTEGER            INFO, LDA, LWORK, N
    !     ..
    !     .. Array Arguments ..
    INTEGER            IPIV( * )
    COMPLEX*16         A( LDA, * ), WORK( * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZGETRI computes the inverse of a matrix using the LU factorization
    !  computed by ZGETRF.
    !
    !  This method inverts U and then computes inv(A) by solving the system
    !  inv(A)*L = inv(U) for inv(A).
    !
    !  Arguments
    !  =========
    !
    !  N       (input) INTEGER
    !          The order of the matrix A.  N >= 0.
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          On entry, the factors L and U from the factorization
    !          A = P*L*U as computed by ZGETRF.
    !          On exit, if INFO = 0, the inverse of the original matrix A.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,N).
    !
    !  IPIV    (input) INTEGER array, dimension (N)
    !          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
    !          matrix was interchanged with row IPIV(i).
    !
    !  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
    !          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
    !
    !  LWORK   (input) INTEGER
    !          The dimension of the array WORK.  LWORK >= max(1,N).
    !          For optimal performance LWORK >= N*NB, where NB is
    !          the optimal blocksize returned by ILAENV.
    !
    !          If LWORK = -1, then a workspace query is assumed; the routine
    !          only calculates the optimal size of the WORK array, returns
    !          this value as the first entry of the WORK array, and no error
    !          message related to LWORK is issued by XERBLA.
    !
    !  INFO    (output) INTEGER
    !          = 0:  successful exit
    !          < 0:  if INFO = -i, the i-th argument had an illegal value
    !          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
    !                singular and its inverse could not be computed.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    COMPLEX*16         ZERO, ONE
    PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), &
         ONE = ( 1.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Local Scalars ..
    LOGICAL            LQUERY
    INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, &
         NBMIN, NN
    !     ..
    !     .. External Functions ..
    !INTEGER            ILAENV
    !EXTERNAL           ILAENV
    !     ..
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA, ZGEMM, ZGEMV, ZSWAP, ZTRSM, ZTRTRI
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC          MAX, MIN
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    NB = ILAENV( 1, 'ZGETRI', ' ', N, -1, -1, -1 )
    LWKOPT = N*NB
    WORK( 1 ) = LWKOPT
    LQUERY = ( LWORK.EQ.-1 )
    IF( N.LT.0 ) THEN
      INFO = -1
    ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
      INFO = -3
    ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
      INFO = -6
    END IF
    IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'ZGETRI', -INFO )
      RETURN
    ELSE IF( LQUERY ) THEN
      RETURN
    END IF
    !
    !     Quick return if possible
    !
    IF( N.EQ.0 )   RETURN
    !
    !     Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
    !     and the inverse is not computed.
    !
    CALL ZTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
    IF( INFO.GT.0 )   RETURN
    !
    NBMIN = 2
    LDWORK = N
    IF( NB.GT.1 .AND. NB.LT.N ) THEN
      IWS = MAX( LDWORK*NB, 1 )
      IF( LWORK.LT.IWS ) THEN
        NB = LWORK / LDWORK
        NBMIN = MAX( 2, ILAENV( 2, 'ZGETRI', ' ', N, -1, -1, -1 ) )
      END IF
    ELSE
      IWS = N
    END IF
    !
    !     Solve the equation inv(A)*L = inv(U) for inv(A).
    !
    IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
      !
      !        Use unblocked code.
      !
      DO J = N, 1, -1
        !
        !           Copy current column of L to WORK and replace with zeros.
        !
        DO I = J + 1, N
          WORK( I ) = A( I, J )
          A( I, J ) = ZERO
        end do
        !
        !           Compute current column of inv(A).
        !
        IF( J.LT.N ) &
             CALL ZGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ), &
             LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
      end do
    ELSE
      !
      !        Use blocked code.
      !
      NN = ( ( N-1 ) / NB )*NB + 1
      DO J = NN, 1, -NB
        JB = MIN( NB, N-J+1 )
        !
        !           Copy current block column of L to WORK and replace with
        !           zeros.
        !
        DO JJ = J, J + JB - 1
          DO I = JJ + 1, N
            WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
            A( I, JJ ) = ZERO
          end do
        end do
        !
        !           Compute current block column of inv(A).
        !
        IF( J+JB.LE.N ) &
             CALL ZGEMM( 'No transpose', 'No transpose', N, JB, &
             N-J-JB+1, -ONE, A( 1, J+JB ), LDA, &
             WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
        CALL ZTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB, &
             ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
      end do
    END IF
    !
    !     Apply column interchanges.
    !
    DO  J = N - 1, 1, -1
      JP = IPIV( J )
      IF( JP.NE.J ) CALL ZSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
    end do
    !
    WORK( 1 ) = IWS
    RETURN
    !
    !     End of ZGETRI
    !
  END SUBROUTINE ZGETRI



  SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
    !
    !  -- LAPACK routine (version 3.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     September 30, 1994
    !
    !     .. Scalar Arguments ..
    CHARACTER          DIAG, UPLO
    INTEGER            INFO, LDA, N
    !     ..
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZTRTI2 computes the inverse of a complex upper or lower triangular
    !  matrix.
    !
    !  This is the Level 2 BLAS version of the algorithm.
    !
    !  Arguments
    !  =========
    !
    !  UPLO    (input) CHARACTER*1
    !          Specifies whether the matrix A is upper or lower triangular.
    !          = 'U':  Upper triangular
    !          = 'L':  Lower triangular
    !
    !  DIAG    (input) CHARACTER*1
    !          Specifies whether or not the matrix A is unit triangular.
    !          = 'N':  Non-unit triangular
    !          = 'U':  Unit triangular
    !
    !  N       (input) INTEGER
    !          The order of the matrix A.  N >= 0.
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          On entry, the triangular matrix A.  If UPLO = 'U', the
    !          leading n by n upper triangular part of the array A contains
    !          the upper triangular matrix, and the strictly lower
    !          triangular part of A is not referenced.  If UPLO = 'L', the
    !          leading n by n lower triangular part of the array A contains
    !          the lower triangular matrix, and the strictly upper
    !          triangular part of A is not referenced.  If DIAG = 'U', the
    !          diagonal elements of A are also not referenced and are
    !          assumed to be 1.
    !
    !          On exit, the (triangular) inverse of the original matrix, in
    !          the same storage format.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,N).
    !
    !  INFO    (output) INTEGER
    !          = 0: successful exit
    !          < 0: if INFO = -k, the k-th argument had an illegal value
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    COMPLEX*16         ONE
    PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Local Scalars ..
    LOGICAL            NOUNIT, UPPER
    INTEGER            J
    COMPLEX*16         AJJ
    !     ..
    !     .. External Functions ..
    !LOGICAL            LSAME
    !EXTERNAL           LSAME
    !     ..
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA, ZSCAL, ZTRMV
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC          MAX
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    UPPER = LSAME( UPLO, 'U' )
    NOUNIT = LSAME( DIAG, 'N' )
    IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
    ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
      INFO = -2
    ELSE IF( N.LT.0 ) THEN
      INFO = -3
    ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
      INFO = -5
    END IF
    IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'ZTRTI2', -INFO )
      RETURN
    END IF
    !
    IF( UPPER ) THEN
      !
      !        Compute inverse of upper triangular matrix.
      !
      DO J = 1, N
        IF( NOUNIT ) THEN
          A( J, J ) = ONE / A( J, J )
          AJJ = -A( J, J )
        ELSE
          AJJ = -ONE
        END IF
        !
        !           Compute elements 1:j-1 of j-th column.
        !
        CALL ZTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA, &
             A( 1, J ), 1 )
        CALL ZSCAL( J-1, AJJ, A( 1, J ), 1 )
      end do
    ELSE
      !
      !        Compute inverse of lower triangular matrix.
      !
      DO J = N, 1, -1
        IF( NOUNIT ) THEN
          A( J, J ) = ONE / A( J, J )
          AJJ = -A( J, J )
        ELSE
          AJJ = -ONE
        END IF
        IF( J.LT.N ) THEN
          !
          !              Compute elements j+1:n of j-th column.
          !
          CALL ZTRMV( 'Lower', 'No transpose', DIAG, N-J, &
               A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
          CALL ZSCAL( N-J, AJJ, A( J+1, J ), 1 )
        END IF
      end do
    END IF
    !
    RETURN
    !
    !     End of ZTRTI2
    !
  END SUBROUTINE ZTRTI2



  SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
    !
    !  -- LAPACK routine (version 3.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     September 30, 1994
    !
    !     .. Scalar Arguments ..
    CHARACTER          DIAG, UPLO
    INTEGER            INFO, LDA, N
    !     ..
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZTRTRI computes the inverse of a complex upper or lower triangular
    !  matrix A.
    !
    !  This is the Level 3 BLAS version of the algorithm.
    !
    !  Arguments
    !  =========
    !
    !  UPLO    (input) CHARACTER*1
    !          = 'U':  A is upper triangular;
    !          = 'L':  A is lower triangular.
    !
    !  DIAG    (input) CHARACTER*1
    !          = 'N':  A is non-unit triangular;
    !          = 'U':  A is unit triangular.
    !
    !  N       (input) INTEGER
    !          The order of the matrix A.  N >= 0.
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          On entry, the triangular matrix A.  If UPLO = 'U', the
    !          leading N-by-N upper triangular part of the array A contains
    !          the upper triangular matrix, and the strictly lower
    !          triangular part of A is not referenced.  If UPLO = 'L', the
    !          leading N-by-N lower triangular part of the array A contains
    !          the lower triangular matrix, and the strictly upper
    !          triangular part of A is not referenced.  If DIAG = 'U', the
    !          diagonal elements of A are also not referenced and are
    !          assumed to be 1.
    !          On exit, the (triangular) inverse of the original matrix, in
    !          the same storage format.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,N).
    !
    !  INFO    (output) INTEGER
    !          = 0: successful exit
    !          < 0: if INFO = -i, the i-th argument had an illegal value
    !          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
    !               matrix is singular and its inverse can not be computed.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    COMPLEX*16         ONE, ZERO
    PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
         ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Local Scalars ..
    LOGICAL            NOUNIT, UPPER
    INTEGER            J, JB, NB, NN
    !     ..
    !     .. External Functions ..
    !LOGICAL            LSAME
    !INTEGER            ILAENV
    !EXTERNAL           LSAME, ILAENV
    !     ..
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA, ZTRMM, ZTRSM, ZTRTI2
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC          MAX, MIN
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    UPPER = LSAME( UPLO, 'U' )
    NOUNIT = LSAME( DIAG, 'N' )
    IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
      INFO = -1
    ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
      INFO = -2
    ELSE IF( N.LT.0 ) THEN
      INFO = -3
    ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
      INFO = -5
    END IF
    IF( INFO.NE.0 ) THEN
      CALL XERBLA( 'ZTRTRI', -INFO )
      RETURN
    END IF
    !
    !     Quick return if possible
    !
    IF( N.EQ.0 )   RETURN
    !
    !     Check for singularity if non-unit.
    !
    IF( NOUNIT ) THEN
      DO INFO = 1, N
        IF( A( INFO, INFO ).EQ.ZERO )   RETURN
      end do
      INFO = 0
    END IF
    !
    !     Determine the block size for this environment.
    !
    NB = ILAENV( 1, 'ZTRTRI', UPLO // DIAG, N, -1, -1, -1 )
    IF( NB.LE.1 .OR. NB.GE.N ) THEN
      !
      !        Use unblocked code
      !
      CALL ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
    ELSE
      !
      !        Use blocked code
      !
      IF( UPPER ) THEN
        !
        !           Compute inverse of upper triangular matrix
        !
        DO J = 1, N, NB
          JB = MIN( NB, N-J+1 )
          !
          !              Compute rows 1:j-1 of current block column
          !
          CALL ZTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1, &
               JB, ONE, A, LDA, A( 1, J ), LDA )
          CALL ZTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1, &
               JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
          !
          !              Compute inverse of current diagonal block
          !
          CALL ZTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
        end do
      ELSE
        !
        !           Compute inverse of lower triangular matrix
        !
        NN = ( ( N-1 ) / NB )*NB + 1
        DO J = NN, 1, -NB
          JB = MIN( NB, N-J+1 )
          IF( J+JB.LE.N ) THEN
            !
            !                 Compute rows j+jb:n of current block column
            !
            CALL ZTRMM( 'Left', 'Lower', 'No transpose', DIAG, &
                 N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, &
                 A( J+JB, J ), LDA )
            CALL ZTRSM( 'Right', 'Lower', 'No transpose', DIAG, &
                 N-J-JB+1, JB, -ONE, A( J, J ), LDA, &
                 A( J+JB, J ), LDA )
          END IF
          !
          !              Compute inverse of current diagonal block
          !
          CALL ZTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
        end do
      END IF
    END IF
    !
    RETURN
    !
    !     End of ZTRTRI
    !
  END SUBROUTINE ZTRTRI


  !================================================================
  !================================================================
  !================================================================

  integer function izamax(n,zx,incx)
    !
    !     finds the index of element having max. absolute value.
    !     jack dongarra, 1/15/85.
    !     modified 3/93 to return if incx .le. 0.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    double complex zx(*)
    double precision smax
    integer i,incx,ix,n
    !double precision dcabs1
    !
    izamax = 0
    if( n.lt.1 .or. incx.le.0 )return
    izamax = 1
    if(n.eq.1)return
    if(incx.eq.1)go to 20
    !
    !        code for increment not equal to 1
    !
    ix = 1
    smax = dcabs1(zx(1))
    ix = ix + incx
    do i = 2,n
      if(dcabs1(zx(ix)).le.smax) go to 5
      izamax = i
      smax = dcabs1(zx(ix))
5     ix = ix + incx
    end do
    return
    !
    !        code for increment equal to 1
    !
20  smax = dcabs1(zx(1))
    do i = 2,n
      if(dcabs1(zx(i)).le.smax) cycle
      izamax = i
      smax = dcabs1(zx(i))
    end do
    return
  end function izamax



  double precision function dcabs1(z)
    double complex z,zz
    double precision t(2)
    equivalence (zz,t(1))
    zz = z
    dcabs1 = dabs(t(1)) + dabs(t(2))
    return
  end function dcabs1



  subroutine  zswap (n,zx,incx,zy,incy)
    !
    !     interchanges two vectors.
    !     jack dongarra, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    double complex zx(*),zy(*),ztemp
    integer i,incx,incy,ix,iy,n
    !
    if(n.le.0)return
    if(incx.eq.1.and.incy.eq.1)go to 20
    !
    !       code for unequal increments or equal increments not equal
    !         to 1
    !
    ix = 1
    iy = 1
    if(incx.lt.0)ix = (-n+1)*incx + 1
    if(incy.lt.0)iy = (-n+1)*incy + 1
    do i = 1,n
      ztemp = zx(ix)
      zx(ix) = zy(iy)
      zy(iy) = ztemp
      ix = ix + incx
      iy = iy + incy
    end do
    return
    !
    !       code for both increments equal to 1
20  do i = 1,n
      ztemp = zx(i)
      zx(i) = zy(i)
      zy(i) = ztemp
    end do
    return
  end subroutine  zswap



  subroutine  zscal(n,za,zx,incx)
    !
    !     scales a vector by a constant.
    !     jack dongarra, 3/11/78.
    !     modified 3/93 to return if incx .le. 0.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    double complex za,zx(*)
    integer i,incx,ix,n
    !
    if( n.le.0 .or. incx.le.0 )return
    if(incx.eq.1)go to 20
    !
    !        code for increment not equal to 1
    !
    ix = 1
    do i = 1,n
      zx(ix) = za*zx(ix)
      ix = ix + incx
    end do
    return
    !
    !        code for increment equal to 1
    !
20  do i = 1,n
      zx(i) = za*zx(i)
    end do
    return
  end  subroutine  zscal



  SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
    !     .. Scalar Arguments ..
    COMPLEX*16         ALPHA
    INTEGER            INCX, INCY, LDA, M, N
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * ), X( * ), Y( * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZGERU  performs the rank 1 operation
    !
    !     A := alpha*x*y' + A,
    !
    !  where alpha is a scalar, x is an m element vector, y is an n element
    !  vector and A is an m by n matrix.
    !
    !  Parameters
    !  ==========
    !
    !  M      - INTEGER.
    !           On entry, M specifies the number of rows of the matrix A.
    !           M must be at least zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry, N specifies the number of columns of the matrix A.
    !           N must be at least zero.
    !           Unchanged on exit.
    !
    !  ALPHA  - COMPLEX*16      .
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  X      - COMPLEX*16       array of dimension at least
    !           ( 1 + ( m - 1 )*abs( INCX ) ).
    !           Before entry, the incremented array X must contain the m
    !           element vector x.
    !           Unchanged on exit.
    !
    !  INCX   - INTEGER.
    !           On entry, INCX specifies the increment for the elements of
    !           X. INCX must not be zero.
    !           Unchanged on exit.
    !
    !  Y      - COMPLEX*16       array of dimension at least
    !           ( 1 + ( n - 1 )*abs( INCY ) ).
    !           Before entry, the incremented array Y must contain the n
    !           element vector y.
    !           Unchanged on exit.
    !
    !  INCY   - INTEGER.
    !           On entry, INCY specifies the increment for the elements of
    !           Y. INCY must not be zero.
    !           Unchanged on exit.
    !
    !  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
    !           Before entry, the leading m by n part of the array A must
    !           contain the matrix of coefficients. On exit, A is
    !           overwritten by the updated matrix.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. LDA must be at least
    !           max( 1, m ).
    !           Unchanged on exit.
    !
    !
    !  Level 2 Blas routine.
    !
    !  -- Written on 22-October-1986.
    !     Jack Dongarra, Argonne National Lab.
    !     Jeremy Du Croz, Nag Central Office.
    !     Sven Hammarling, Nag Central Office.
    !     Richard Hanson, Sandia National Labs.
    !
    !
    !     .. Parameters ..
    COMPLEX*16         ZERO
    PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     .. Local Scalars ..
    COMPLEX*16         TEMP
    INTEGER            I, INFO, IX, J, JY, KX
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA
    !     .. Intrinsic Functions ..
    INTRINSIC          MAX
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF     ( M.LT.0 )THEN
      INFO = 1
    ELSE IF( N.LT.0 )THEN
      INFO = 2
    ELSE IF( INCX.EQ.0 )THEN
      INFO = 5
    ELSE IF( INCY.EQ.0 )THEN
      INFO = 7
    ELSE IF( LDA.LT.MAX( 1, M ) )THEN
      INFO = 9
    END IF
    IF( INFO.NE.0 )THEN
      CALL XERBLA( 'ZGERU ', INFO )
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )   RETURN
    !
    !     Start the operations. In this version the elements of A are
    !     accessed sequentially with one pass through A.
    !
    IF( INCY.GT.0 )THEN
      JY = 1
    ELSE
      JY = 1 - ( N - 1 )*INCY
    END IF
    IF( INCX.EQ.1 )THEN
      DO J = 1, N
        IF( Y( JY ).NE.ZERO )THEN
          TEMP = ALPHA*Y( JY )
          DO I = 1, M
            A( I, J ) = A( I, J ) + X( I )*TEMP
          end do
        END IF
        JY = JY + INCY
      end do
    ELSE
      IF( INCX.GT.0 )THEN
        KX = 1
      ELSE
        KX = 1 - ( M - 1 )*INCX
      END IF
      DO  J = 1, N
        IF( Y( JY ).NE.ZERO )THEN
          TEMP = ALPHA*Y( JY )
          IX   = KX
          DO I = 1, M
            A( I, J ) = A( I, J ) + X( IX )*TEMP
            IX        = IX        + INCX
          end do
        END IF
        JY = JY + INCY
      end do
    END IF
    !
    RETURN
    !
    !     End of ZGERU .
    !
  END SUBROUTINE ZGERU



  SUBROUTINE ZTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
       B, LDB )
    !     .. Scalar Arguments ..
    CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
    INTEGER            M, N, LDA, LDB
    COMPLEX*16         ALPHA
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * ), B( LDB, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZTRSM  solves one of the matrix equations
    !
    !     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
    !
    !  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
    !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
    !
    !     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
    !
    !  The matrix X is overwritten on B.
    !
    !  Parameters
    !  ==========
    !
    !  SIDE   - CHARACTER*1.
    !           On entry, SIDE specifies whether op( A ) appears on the left
    !           or right of X as follows:
    !
    !              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
    !
    !              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
    !
    !           Unchanged on exit.
    !
    !  UPLO   - CHARACTER*1.
    !           On entry, UPLO specifies whether the matrix A is an upper or
    !           lower triangular matrix as follows:
    !
    !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
    !
    !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
    !
    !           Unchanged on exit.
    !
    !  TRANSA - CHARACTER*1.
    !           On entry, TRANSA specifies the form of op( A ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSA = 'N' or 'n'   op( A ) = A.
    !
    !              TRANSA = 'T' or 't'   op( A ) = A'.
    !
    !              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
    !
    !           Unchanged on exit.
    !
    !  DIAG   - CHARACTER*1.
    !           On entry, DIAG specifies whether or not A is unit triangular
    !           as follows:
    !
    !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
    !
    !              DIAG = 'N' or 'n'   A is not assumed to be unit
    !                                  triangular.
    !
    !           Unchanged on exit.
    !
    !  M      - INTEGER.
    !           On entry, M specifies the number of rows of B. M must be at
    !           least zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry, N specifies the number of columns of B.  N must be
    !           at least zero.
    !           Unchanged on exit.
    !
    !  ALPHA  - COMPLEX*16      .
    !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
    !           zero then  A is not referenced and  B need not be set before
    !           entry.
    !           Unchanged on exit.
    !
    !  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
    !           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
    !           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
    !           upper triangular part of the array  A must contain the upper
    !           triangular matrix  and the strictly lower triangular part of
    !           A is not referenced.
    !           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
    !           lower triangular part of the array  A must contain the lower
    !           triangular matrix  and the strictly upper triangular part of
    !           A is not referenced.
    !           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
    !           A  are not referenced either,  but are assumed to be  unity.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
    !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
    !           then LDA must be at least max( 1, n ).
    !           Unchanged on exit.
    !
    !  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
    !           Before entry,  the leading  m by n part of the array  B must
    !           contain  the  right-hand  side  matrix  B,  and  on exit  is
    !           overwritten by the solution matrix  X.
    !
    !  LDB    - INTEGER.
    !           On entry, LDB specifies the first dimension of B as declared
    !           in  the  calling  (sub)  program.   LDB  must  be  at  least
    !           max( 1, m ).
    !           Unchanged on exit.
    !
    !
    !  Level 3 Blas routine.
    !
    !  -- Written on 8-February-1989.
    !     Jack Dongarra, Argonne National Laboratory.
    !     Iain Duff, AERE Harwell.
    !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !     Sven Hammarling, Numerical Algorithms Group Ltd.
    !
    !
    !     .. External Functions ..
    !LOGICAL            LSAME
    !EXTERNAL           LSAME
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA
    !     .. Intrinsic Functions ..
    INTRINSIC          DCONJG, MAX
    !     .. Local Scalars ..
    LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
    INTEGER            I, INFO, J, K, NROWA
    COMPLEX*16         TEMP
    !     .. Parameters ..
    COMPLEX*16         ONE
    PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
    COMPLEX*16         ZERO
    PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    LSIDE  = LSAME( SIDE  , 'L' )
    IF( LSIDE )THEN
      NROWA = M
    ELSE
      NROWA = N
    END IF
    NOCONJ = LSAME( TRANSA, 'T' )
    NOUNIT = LSAME( DIAG  , 'N' )
    UPPER  = LSAME( UPLO  , 'U' )
    !
    INFO   = 0
    IF(      ( .NOT.LSIDE                ).AND. &
         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
      INFO = 1
    ELSE IF( ( .NOT.UPPER                ).AND. &
         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
      INFO = 2
    ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
         ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
      INFO = 3
    ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
      INFO = 4
    ELSE IF( M  .LT.0               )THEN
      INFO = 5
    ELSE IF( N  .LT.0               )THEN
      INFO = 6
    ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
      INFO = 9
    ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
      INFO = 11
    END IF
    IF( INFO.NE.0 )THEN
      CALL XERBLA( 'ZTRSM ', INFO )
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF( N.EQ.0 )   RETURN
    !
    !     And when  alpha.eq.zero.
    !
    IF( ALPHA.EQ.ZERO )THEN
      DO J = 1, N
        DO I = 1, M
          B( I, J ) = ZERO
        end do
      end do
      RETURN
    END IF
    !
    !     Start the operations.
    !
    IF( LSIDE )THEN
      IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*inv( A )*B.
        !
        IF( UPPER )THEN
          DO J = 1, N
            IF( ALPHA.NE.ONE )THEN
              DO I = 1, M
                B( I, J ) = ALPHA*B( I, J )
              end do
            END IF
            DO K = M, 1, -1
              IF( B( K, J ).NE.ZERO )THEN
                IF( NOUNIT )  B( K, J ) = B( K, J )/A( K, K )
                DO I = 1, K - 1
                  B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                end do
              END IF
            end do
          end do
        ELSE
          DO  J = 1, N
            IF( ALPHA.NE.ONE )THEN
              DO I = 1, M
                B( I, J ) = ALPHA*B( I, J )
              end do
            END IF
            DO K = 1, M
              IF( B( K, J ).NE.ZERO )THEN
                IF( NOUNIT ) B( K, J ) = B( K, J )/A( K, K )
                DO I = K + 1, M
                  B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                end do
              END IF
            end do
          end do
        END IF
      ELSE
        !
        !           Form  B := alpha*inv( A' )*B
        !           or    B := alpha*inv( conjg( A' ) )*B.
        !
        IF( UPPER )THEN
          DO J = 1, N
            DO I = 1, M
              TEMP = ALPHA*B( I, J )
              IF( NOCONJ )THEN
                DO K = 1, I - 1
                  TEMP = TEMP - A( K, I )*B( K, J )
                end do
                IF( NOUNIT ) TEMP = TEMP/A( I, I )
              ELSE
                DO K = 1, I - 1
                  TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
                end do
                IF( NOUNIT ) TEMP = TEMP/DCONJG( A( I, I ) )
              END IF
              B( I, J ) = TEMP
            end do
          end do
        ELSE
          DO J = 1, N
            DO I = M, 1, -1
              TEMP = ALPHA*B( I, J )
              IF( NOCONJ )THEN
                DO K = I + 1, M
                  TEMP = TEMP - A( K, I )*B( K, J )
                end do
                IF( NOUNIT ) TEMP = TEMP/A( I, I )
              ELSE
                DO K = I + 1, M
                  TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
                end do
                IF( NOUNIT ) TEMP = TEMP/DCONJG( A( I, I ) )
              END IF
              B( I, J ) = TEMP
            end do
          end do
        END IF
      END IF
    ELSE
      IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*B*inv( A ).
        !
        IF( UPPER )THEN
          DO J = 1, N
            IF( ALPHA.NE.ONE )THEN
              DO I = 1, M
                B( I, J ) = ALPHA*B( I, J )
              end do
            END IF
            DO K = 1, J - 1
              IF( A( K, J ).NE.ZERO )THEN
                DO I = 1, M
                  B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                end do
              END IF
            end do

            IF( NOUNIT )THEN
              TEMP = ONE/A( J, J )
              DO  I = 1, M
                B( I, J ) = TEMP*B( I, J )
              end do
            END IF
          end do

        ELSE
          DO J = N, 1, -1
            IF( ALPHA.NE.ONE )THEN
              DO I = 1, M
                B( I, J ) = ALPHA*B( I, J )
              end do
            END IF
            DO K = J + 1, N
              IF( A( K, J ).NE.ZERO )THEN
                DO I = 1, M
                  B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                end do
              END IF
            end do
            IF( NOUNIT )THEN
              TEMP = ONE/A( J, J )
              DO I = 1, M
                B( I, J ) = TEMP*B( I, J )
              end do
            END IF
          end do
        END IF
      ELSE
        !
        !           Form  B := alpha*B*inv( A' )
        !           or    B := alpha*B*inv( conjg( A' ) ).
        !
        IF( UPPER )THEN
          DO K = N, 1, -1
            IF( NOUNIT )THEN
              IF( NOCONJ )THEN
                TEMP = ONE/A( K, K )
              ELSE
                TEMP = ONE/DCONJG( A( K, K ) )
              END IF
              DO  I = 1, M
                B( I, K ) = TEMP*B( I, K )
              end do
            END IF
            DO J = 1, K - 1
              IF( A( J, K ).NE.ZERO )THEN
                IF( NOCONJ )THEN
                  TEMP = A( J, K )
                ELSE
                  TEMP = DCONJG( A( J, K ) )
                END IF
                DO I = 1, M
                  B( I, J ) = B( I, J ) - TEMP*B( I, K )
                end do
              END IF
            end do
            IF( ALPHA.NE.ONE )THEN
              DO  I = 1, M
                B( I, K ) = ALPHA*B( I, K )
              end do
            END IF
          end do
        ELSE
          DO K = 1, N
            IF( NOUNIT )THEN
              IF( NOCONJ )THEN
                TEMP = ONE/A( K, K )
              ELSE
                TEMP = ONE/DCONJG( A( K, K ) )
              END IF
              DO I = 1, M
                B( I, K ) = TEMP*B( I, K )
              end do
            END IF
            DO J = K + 1, N
              IF( A( J, K ).NE.ZERO )THEN
                IF( NOCONJ )THEN
                  TEMP = A( J, K )
                ELSE
                  TEMP = DCONJG( A( J, K ) )
                END IF
                DO I = 1, M
                  B( I, J ) = B( I, J ) - TEMP*B( I, K )
                end do
              END IF
            end do
            IF( ALPHA.NE.ONE )THEN
              DO I = 1, M
                B( I, K ) = ALPHA*B( I, K )
              end do
            END IF
          end do
        END IF
      END IF
    END IF
    !
    RETURN
    !
    !     End of ZTRSM .
    !
  END SUBROUTINE ZTRSM



  SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
       BETA, C, LDC )
    !     .. Scalar Arguments ..
    CHARACTER*1        TRANSA, TRANSB
    INTEGER            M, N, K, LDA, LDB, LDC
    COMPLEX*16         ALPHA, BETA
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZGEMM  performs one of the matrix-matrix operations
    !
    !     C := alpha*op( A )*op( B ) + beta*C,
    !
    !  where  op( X ) is one of
    !
    !     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
    !
    !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
    !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
    !
    !  Parameters
    !  ==========
    !
    !  TRANSA - CHARACTER*1.
    !           On entry, TRANSA specifies the form of op( A ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSA = 'N' or 'n',  op( A ) = A.
    !
    !              TRANSA = 'T' or 't',  op( A ) = A'.
    !
    !              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
    !
    !           Unchanged on exit.
    !
    !  TRANSB - CHARACTER*1.
    !           On entry, TRANSB specifies the form of op( B ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSB = 'N' or 'n',  op( B ) = B.
    !
    !              TRANSB = 'T' or 't',  op( B ) = B'.
    !
    !              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
    !
    !           Unchanged on exit.
    !
    !  M      - INTEGER.
    !           On entry,  M  specifies  the number  of rows  of the  matrix
    !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry,  N  specifies the number  of columns of the matrix
    !           op( B ) and the number of columns of the matrix C. N must be
    !           at least zero.
    !           Unchanged on exit.
    !
    !  K      - INTEGER.
    !           On entry,  K  specifies  the number of columns of the matrix
    !           op( A ) and the number of rows of the matrix op( B ). K must
    !           be at least  zero.
    !           Unchanged on exit.
    !
    !  ALPHA  - COMPLEX*16      .
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
    !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
    !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
    !           part of the array  A  must contain the matrix  A,  otherwise
    !           the leading  k by m  part of the array  A  must contain  the
    !           matrix A.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
    !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
    !           least  max( 1, k ).
    !           Unchanged on exit.
    !
    !  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
    !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
    !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
    !           part of the array  B  must contain the matrix  B,  otherwise
    !           the leading  n by k  part of the array  B  must contain  the
    !           matrix B.
    !           Unchanged on exit.
    !
    !  LDB    - INTEGER.
    !           On entry, LDB specifies the first dimension of B as declared
    !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
    !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
    !           least  max( 1, n ).
    !           Unchanged on exit.
    !
    !  BETA   - COMPLEX*16      .
    !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
    !           supplied as zero then C need not be set on input.
    !           Unchanged on exit.
    !
    !  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
    !           Before entry, the leading  m by n  part of the array  C must
    !           contain the matrix  C,  except when  beta  is zero, in which
    !           case C need not be set on entry.
    !           On exit, the array  C  is overwritten by the  m by n  matrix
    !           ( alpha*op( A )*op( B ) + beta*C ).
    !
    !  LDC    - INTEGER.
    !           On entry, LDC specifies the first dimension of C as declared
    !           in  the  calling  (sub)  program.   LDC  must  be  at  least
    !           max( 1, m ).
    !           Unchanged on exit.
    !
    !
    !  Level 3 Blas routine.
    !
    !  -- Written on 8-February-1989.
    !     Jack Dongarra, Argonne National Laboratory.
    !     Iain Duff, AERE Harwell.
    !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !     Sven Hammarling, Numerical Algorithms Group Ltd.
    !
    !
    !     .. External Functions ..
    !LOGICAL            LSAME
    !EXTERNAL           LSAME
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA
    !     .. Intrinsic Functions ..
    INTRINSIC          DCONJG, MAX
    !     .. Local Scalars ..
    LOGICAL            CONJA, CONJB, NOTA, NOTB
    INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
    COMPLEX*16         TEMP
    !     .. Parameters ..
    COMPLEX*16         ONE
    PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
    COMPLEX*16         ZERO
    PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Executable Statements ..
    !
    !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    !     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
    !     B  respectively are to be  transposed but  not conjugated  and set
    !     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
    !     and the number of rows of  B  respectively.
    !
    NOTA  = LSAME( TRANSA, 'N' )
    NOTB  = LSAME( TRANSB, 'N' )
    CONJA = LSAME( TRANSA, 'C' )
    CONJB = LSAME( TRANSB, 'C' )
    IF( NOTA )THEN
      NROWA = M
      NCOLA = K
    ELSE
      NROWA = K
      NCOLA = M
    END IF
    IF( NOTB )THEN
      NROWB = K
    ELSE
      NROWB = N
    END IF
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF(      ( .NOT.NOTA                 ).AND. &
         ( .NOT.CONJA                ).AND. &
         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
      INFO = 1
    ELSE IF( ( .NOT.NOTB                 ).AND. &
         ( .NOT.CONJB                ).AND. &
         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
      INFO = 2
    ELSE IF( M  .LT.0               )THEN
      INFO = 3
    ELSE IF( N  .LT.0               )THEN
      INFO = 4
    ELSE IF( K  .LT.0               )THEN
      INFO = 5
    ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
      INFO = 8
    ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
      INFO = 10
    ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
      INFO = 13
    END IF
    IF( INFO.NE.0 )THEN
      CALL XERBLA( 'ZGEMM ', INFO )
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
         ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
         RETURN
    !
    !     And when  alpha.eq.zero.
    !
    IF( ALPHA.EQ.ZERO )THEN
      IF( BETA.EQ.ZERO )THEN
        DO J = 1, N
          DO I = 1, M
            C( I, J ) = ZERO
          end do
        end do
      ELSE
        DO J = 1, N
          DO I = 1, M
            C( I, J ) = BETA*C( I, J )
          end do
        end do
      END IF
      RETURN
    END IF
    !
    !     Start the operations.
    !
    IF( NOTB )THEN
      IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B + beta*C.
        !
        DO J = 1, N
          IF( BETA.EQ.ZERO )THEN
            DO I = 1, M
              C( I, J ) = ZERO
            end do

          ELSE IF( BETA.NE.ONE )THEN
            DO I = 1, M
              C( I, J ) = BETA*C( I, J )
            end do
          END IF
          DO L = 1, K
            IF( B( L, J ).NE.ZERO )THEN
              TEMP = ALPHA*B( L, J )
              DO I = 1, M
                C( I, J ) = C( I, J ) + TEMP*A( I, L )
              end do
            END IF
          end do
        end do
      ELSE IF( CONJA )THEN
        !
        !           Form  C := alpha*conjg( A' )*B + beta*C.
        !
        DO J = 1, N
          DO I = 1, M
            TEMP = ZERO
            DO L = 1, K
              TEMP = TEMP + DCONJG( A( L, I ) )*B( L, J )
            end do
            IF( BETA.EQ.ZERO )THEN
              C( I, J ) = ALPHA*TEMP
            ELSE
              C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
            END IF
          end do
        end do
      ELSE
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        DO J = 1, N
          DO I = 1, M
            TEMP = ZERO
            DO L = 1, K
              TEMP = TEMP + A( L, I )*B( L, J )
            end do
            IF( BETA.EQ.ZERO )THEN
              C( I, J ) = ALPHA*TEMP
            ELSE
              C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
            END IF
          end do
        end do
      END IF
    ELSE IF( NOTA )THEN
      IF( CONJB )THEN
        !
        !           Form  C := alpha*A*conjg( B' ) + beta*C.
        !
        DO J = 1, N
          IF( BETA.EQ.ZERO )THEN
            DO I = 1, M
              C( I, J ) = ZERO
            end do
          ELSE IF( BETA.NE.ONE )THEN
            DO I = 1, M
              C( I, J ) = BETA*C( I, J )
            end do
          END IF
          DO L = 1, K
            IF( B( J, L ).NE.ZERO )THEN
              TEMP = ALPHA*DCONJG( B( J, L ) )
              DO I = 1, M
                C( I, J ) = C( I, J ) + TEMP*A( I, L )
              end do
            END IF
          end do
        end do
      ELSE
        !
        !           Form  C := alpha*A*B'          + beta*C
        !
        DO J = 1, N
          IF( BETA.EQ.ZERO )THEN
            DO I = 1, M
              C( I, J ) = ZERO
            end do
          ELSE IF( BETA.NE.ONE )THEN
            DO I = 1, M
              C( I, J ) = BETA*C( I, J )
            end do
          END IF
          DO L = 1, K
            IF( B( J, L ).NE.ZERO )THEN
              TEMP = ALPHA*B( J, L )
              DO I = 1, M
                C( I, J ) = C( I, J ) + TEMP*A( I, L )
              end do
            END IF
          end do
        end do
      END IF
    ELSE IF( CONJA )THEN
      IF( CONJB )THEN
        !
        !           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
        !
        DO J = 1, N
          DO I = 1, M
            TEMP = ZERO
            DO L = 1, K
              TEMP = TEMP + DCONJG( A( L, I ) )*DCONJG( B( J, L ) )
            end do
            IF( BETA.EQ.ZERO )THEN
              C( I, J ) = ALPHA*TEMP
            ELSE
              C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
            END IF
          end do
        end do
      ELSE
        !
        !           Form  C := alpha*conjg( A' )*B' + beta*C
        !
        DO J = 1, N
          DO I = 1, M
            TEMP = ZERO
            DO L = 1, K
              TEMP = TEMP + DCONJG( A( L, I ) )*B( J, L )
            end do
            IF( BETA.EQ.ZERO )THEN
              C( I, J ) = ALPHA*TEMP
            ELSE
              C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
            END IF
          end do
        end do
      END IF
    ELSE
      IF( CONJB )THEN
        !
        !           Form  C := alpha*A'*conjg( B' ) + beta*C
        !
        DO J = 1, N
          DO I = 1, M
            TEMP = ZERO
            DO L = 1, K
              TEMP = TEMP + A( L, I )*DCONJG( B( J, L ) )
            end do
            IF( BETA.EQ.ZERO )THEN
              C( I, J ) = ALPHA*TEMP
            ELSE
              C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
            END IF
          end do
        end do
      ELSE
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        DO J = 1, N
          DO I = 1, M
            TEMP = ZERO
            DO L = 1, K
              TEMP = TEMP + A( L, I )*B( J, L )
            end do
            IF( BETA.EQ.ZERO )THEN
              C( I, J ) = ALPHA*TEMP
            ELSE
              C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
            END IF
          end do
        end do
      END IF
    END IF
    !
    RETURN
    !
    !     End of ZGEMM .
    !
  END SUBROUTINE ZGEMM




  SUBROUTINE ZTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
    !     .. Scalar Arguments ..
    INTEGER            INCX, LDA, N
    CHARACTER*1        DIAG, TRANS, UPLO
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * ), X( * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZTRMV  performs one of the matrix-vector operations
    !
    !     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
    !
    !  where x is an n element vector and  A is an n by n unit, or non-unit,
    !  upper or lower triangular matrix.
    !
    !  Parameters
    !  ==========
    !
    !  UPLO   - CHARACTER*1.
    !           On entry, UPLO specifies whether the matrix is an upper or
    !           lower triangular matrix as follows:
    !
    !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
    !
    !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
    !
    !           Unchanged on exit.
    !
    !  TRANS  - CHARACTER*1.
    !           On entry, TRANS specifies the operation to be performed as
    !           follows:
    !
    !              TRANS = 'N' or 'n'   x := A*x.
    !
    !              TRANS = 'T' or 't'   x := A'*x.
    !
    !              TRANS = 'C' or 'c'   x := conjg( A' )*x.
    !
    !           Unchanged on exit.
    !
    !  DIAG   - CHARACTER*1.
    !           On entry, DIAG specifies whether or not A is unit
    !           triangular as follows:
    !
    !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
    !
    !              DIAG = 'N' or 'n'   A is not assumed to be unit
    !                                  triangular.
    !
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry, N specifies the order of the matrix A.
    !           N must be at least zero.
    !           Unchanged on exit.
    !
    !  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
    !           Before entry with  UPLO = 'U' or 'u', the leading n by n
    !           upper triangular part of the array A must contain the upper
    !           triangular matrix and the strictly lower triangular part of
    !           A is not referenced.
    !           Before entry with UPLO = 'L' or 'l', the leading n by n
    !           lower triangular part of the array A must contain the lower
    !           triangular matrix and the strictly upper triangular part of
    !           A is not referenced.
    !           Note that when  DIAG = 'U' or 'u', the diagonal elements of
    !           A are not referenced either, but are assumed to be unity.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. LDA must be at least
    !           max( 1, n ).
    !           Unchanged on exit.
    !
    !  X      - COMPLEX*16       array of dimension at least
    !           ( 1 + ( n - 1 )*abs( INCX ) ).
    !           Before entry, the incremented array X must contain the n
    !           element vector x. On exit, X is overwritten with the
    !           tranformed vector x.
    !
    !  INCX   - INTEGER.
    !           On entry, INCX specifies the increment for the elements of
    !           X. INCX must not be zero.
    !           Unchanged on exit.
    !
    !
    !  Level 2 Blas routine.
    !
    !  -- Written on 22-October-1986.
    !     Jack Dongarra, Argonne National Lab.
    !     Jeremy Du Croz, Nag Central Office.
    !     Sven Hammarling, Nag Central Office.
    !     Richard Hanson, Sandia National Labs.
    !
    !
    !     .. Parameters ..
    COMPLEX*16         ZERO
    PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     .. Local Scalars ..
    COMPLEX*16         TEMP
    INTEGER            I, INFO, IX, J, JX, KX
    LOGICAL            NOCONJ, NOUNIT
    !     .. External Functions ..
    !LOGICAL            LSAME
    !EXTERNAL           LSAME
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA
    !     .. Intrinsic Functions ..
    INTRINSIC          DCONJG, MAX
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF     ( .NOT.LSAME( UPLO , 'U' ).AND. &
         .NOT.LSAME( UPLO , 'L' )      )THEN
      INFO = 1
    ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND. &
         .NOT.LSAME( TRANS, 'T' ).AND. &
         .NOT.LSAME( TRANS, 'C' )      )THEN
      INFO = 2
    ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND. &
         .NOT.LSAME( DIAG , 'N' )      )THEN
      INFO = 3
    ELSE IF( N.LT.0 )THEN
      INFO = 4
    ELSE IF( LDA.LT.MAX( 1, N ) )THEN
      INFO = 6
    ELSE IF( INCX.EQ.0 )THEN
      INFO = 8
    END IF
    IF( INFO.NE.0 )THEN
      CALL XERBLA( 'ZTRMV ', INFO )
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF( N.EQ.0 )   RETURN
    !
    NOCONJ = LSAME( TRANS, 'T' )
    NOUNIT = LSAME( DIAG , 'N' )
    !
    !     Set up the start point in X if the increment is not unity. This
    !     will be  ( N - 1 )*INCX  too small for descending loops.
    !
    IF( INCX.LE.0 )THEN
      KX = 1 - ( N - 1 )*INCX
    ELSE IF( INCX.NE.1 )THEN
      KX = 1
    END IF
    !
    !     Start the operations. In this version the elements of A are
    !     accessed sequentially with one pass through A.
    !
    IF( LSAME( TRANS, 'N' ) )THEN
      !
      !        Form  x := A*x.
      !
      IF( LSAME( UPLO, 'U' ) )THEN
        IF( INCX.EQ.1 )THEN
          DO J = 1, N
            IF( X( J ).NE.ZERO )THEN
              TEMP = X( J )
              DO I = 1, J - 1
                X( I ) = X( I ) + TEMP*A( I, J )
              end do
              IF( NOUNIT ) X( J ) = X( J )*A( J, J )
            END IF
          end do
        ELSE
          JX = KX
          DO J = 1, N
            IF( X( JX ).NE.ZERO )THEN
              TEMP = X( JX )
              IX   = KX
              DO I = 1, J - 1
                X( IX ) = X( IX ) + TEMP*A( I, J )
                IX      = IX      + INCX
              end do
              IF( NOUNIT ) X( JX ) = X( JX )*A( J, J )
            END IF
            JX = JX + INCX
          end do
        END IF
      ELSE
        IF( INCX.EQ.1 )THEN
          DO J = N, 1, -1
            IF( X( J ).NE.ZERO )THEN
              TEMP = X( J )
              DO I = N, J + 1, -1
                X( I ) = X( I ) + TEMP*A( I, J )
              end do
              IF( NOUNIT ) X( J ) = X( J )*A( J, J )
            END IF
          end do
        ELSE
          KX = KX + ( N - 1 )*INCX
          JX = KX
          DO J = N, 1, -1
            IF( X( JX ).NE.ZERO )THEN
              TEMP = X( JX )
              IX   = KX
              DO I = N, J + 1, -1
                X( IX ) = X( IX ) + TEMP*A( I, J )
                IX      = IX      - INCX
              end do
              IF( NOUNIT )  X( JX ) = X( JX )*A( J, J )
            END IF
            JX = JX - INCX
          end do
        END IF
      END IF
    ELSE
      !
      !        Form  x := A'*x  or  x := conjg( A' )*x.
      !
      IF( LSAME( UPLO, 'U' ) )THEN
        IF( INCX.EQ.1 )THEN
          DO J = N, 1, -1
            TEMP = X( J )
            IF( NOCONJ )THEN
              IF( NOUNIT )  TEMP = TEMP*A( J, J )
              DO I = J - 1, 1, -1
                TEMP = TEMP + A( I, J )*X( I )
              end do
            ELSE
              IF( NOUNIT )  TEMP = TEMP*DCONJG( A( J, J ) )
              DO I = J - 1, 1, -1
                TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
              end do
            END IF
            X( J ) = TEMP
          end do
        ELSE
          JX = KX + ( N - 1 )*INCX
          DO J = N, 1, -1
            TEMP = X( JX )
            IX   = JX
            IF( NOCONJ )THEN
              IF( NOUNIT )    TEMP = TEMP*A( J, J )
              DO I = J - 1, 1, -1
                IX   = IX   - INCX
                TEMP = TEMP + A( I, J )*X( IX )
              end do
            ELSE
              IF( NOUNIT )  TEMP = TEMP*DCONJG( A( J, J ) )
              DO I = J - 1, 1, -1
                IX   = IX   - INCX
                TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
              end do
            END IF
            X( JX ) = TEMP
            JX      = JX   - INCX
          end do
        END IF
      ELSE
        IF( INCX.EQ.1 )THEN
          DO J = 1, N
            TEMP = X( J )
            IF( NOCONJ )THEN
              IF( NOUNIT ) TEMP = TEMP*A( J, J )
              DO I = J + 1, N
                TEMP = TEMP + A( I, J )*X( I )
              end do
            ELSE
              IF( NOUNIT )  TEMP = TEMP*DCONJG( A( J, J ) )
              DO I = J + 1, N
                TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
              end do
            END IF
            X( J ) = TEMP
          end do
        ELSE
          JX = KX
          DO J = 1, N
            TEMP = X( JX )
            IX   = JX
            IF( NOCONJ )THEN
              IF( NOUNIT )  TEMP = TEMP*A( J, J )
              DO I = J + 1, N
                IX   = IX   + INCX
                TEMP = TEMP + A( I, J )*X( IX )
              end do
            ELSE
              IF( NOUNIT ) TEMP = TEMP*DCONJG( A( J, J ) )
              DO I = J + 1, N
                IX   = IX   + INCX
                TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
              end do
            END IF
            X( JX ) = TEMP
            JX      = JX   + INCX
          end do
        END IF
      END IF
    END IF
    !
    RETURN
    !
    !     End of ZTRMV .
    !
  END SUBROUTINE ZTRMV



  SUBROUTINE ZTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
       B, LDB )
    !     .. Scalar Arguments ..
    CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
    INTEGER            M, N, LDA, LDB
    COMPLEX*16         ALPHA
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * ), B( LDB, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZTRMM  performs one of the matrix-matrix operations
    !
    !     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
    !
    !  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
    !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
    !
    !     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
    !
    !  Parameters
    !  ==========
    !
    !  SIDE   - CHARACTER*1.
    !           On entry,  SIDE specifies whether  op( A ) multiplies B from
    !           the left or right as follows:
    !
    !              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
    !
    !              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
    !
    !           Unchanged on exit.
    !
    !  UPLO   - CHARACTER*1.
    !           On entry, UPLO specifies whether the matrix A is an upper or
    !           lower triangular matrix as follows:
    !
    !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
    !
    !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
    !
    !           Unchanged on exit.
    !
    !  TRANSA - CHARACTER*1.
    !           On entry, TRANSA specifies the form of op( A ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSA = 'N' or 'n'   op( A ) = A.
    !
    !              TRANSA = 'T' or 't'   op( A ) = A'.
    !
    !              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
    !
    !           Unchanged on exit.
    !
    !  DIAG   - CHARACTER*1.
    !           On entry, DIAG specifies whether or not A is unit triangular
    !           as follows:
    !
    !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
    !
    !              DIAG = 'N' or 'n'   A is not assumed to be unit
    !                                  triangular.
    !
    !           Unchanged on exit.
    !
    !  M      - INTEGER.
    !           On entry, M specifies the number of rows of B. M must be at
    !           least zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry, N specifies the number of columns of B.  N must be
    !           at least zero.
    !           Unchanged on exit.
    !
    !  ALPHA  - COMPLEX*16      .
    !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
    !           zero then  A is not referenced and  B need not be set before
    !           entry.
    !           Unchanged on exit.
    !
    !  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
    !           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
    !           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
    !           upper triangular part of the array  A must contain the upper
    !           triangular matrix  and the strictly lower triangular part of
    !           A is not referenced.
    !           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
    !           lower triangular part of the array  A must contain the lower
    !           triangular matrix  and the strictly upper triangular part of
    !           A is not referenced.
    !           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
    !           A  are not referenced either,  but are assumed to be  unity.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
    !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
    !           then LDA must be at least max( 1, n ).
    !           Unchanged on exit.
    !
    !  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
    !           Before entry,  the leading  m by n part of the array  B must
    !           contain the matrix  B,  and  on exit  is overwritten  by the
    !           transformed matrix.
    !
    !  LDB    - INTEGER.
    !           On entry, LDB specifies the first dimension of B as declared
    !           in  the  calling  (sub)  program.   LDB  must  be  at  least
    !           max( 1, m ).
    !           Unchanged on exit.
    !
    !
    !  Level 3 Blas routine.
    !
    !  -- Written on 8-February-1989.
    !     Jack Dongarra, Argonne National Laboratory.
    !     Iain Duff, AERE Harwell.
    !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !     Sven Hammarling, Numerical Algorithms Group Ltd.
    !
    !
    !     .. External Functions ..
    !LOGICAL            LSAME
    !EXTERNAL           LSAME
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA
    !     .. Intrinsic Functions ..
    INTRINSIC          DCONJG, MAX
    !     .. Local Scalars ..
    LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
    INTEGER            I, INFO, J, K, NROWA
    COMPLEX*16         TEMP
    !     .. Parameters ..
    COMPLEX*16         ONE
    PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
    COMPLEX*16         ZERO
    PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    LSIDE  = LSAME( SIDE  , 'L' )
    IF( LSIDE )THEN
      NROWA = M
    ELSE
      NROWA = N
    END IF
    NOCONJ = LSAME( TRANSA, 'T' )
    NOUNIT = LSAME( DIAG  , 'N' )
    UPPER  = LSAME( UPLO  , 'U' )
    !
    INFO   = 0
    IF(      ( .NOT.LSIDE                ).AND. &
         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
      INFO = 1
    ELSE IF( ( .NOT.UPPER                ).AND. &
         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
      INFO = 2
    ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
         ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
      INFO = 3
    ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
      INFO = 4
    ELSE IF( M  .LT.0               )THEN
      INFO = 5
    ELSE IF( N  .LT.0               )THEN
      INFO = 6
    ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
      INFO = 9
    ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
      INFO = 11
    END IF
    IF( INFO.NE.0 )THEN
      CALL XERBLA( 'ZTRMM ', INFO )
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF( N.EQ.0 )   RETURN
    !
    !     And when  alpha.eq.zero.
    !
    IF( ALPHA.EQ.ZERO )THEN
      DO J = 1, N
        DO I = 1, M
          B( I, J ) = ZERO
        end do
      end do
      RETURN
    END IF
    !
    !     Start the operations.
    !
    IF( LSIDE )THEN
      IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*A*B.
        !
        IF( UPPER )THEN
          DO J = 1, N
            DO K = 1, M
              IF( B( K, J ).NE.ZERO )THEN
                TEMP = ALPHA*B( K, J )
                DO I = 1, K - 1
                  B( I, J ) = B( I, J ) + TEMP*A( I, K )
                end do
                IF( NOUNIT ) TEMP = TEMP*A( K, K )
                B( K, J ) = TEMP
              END IF
            end do
          end do
        ELSE
          DO J = 1, N
            DO K = M, 1, -1
              IF( B( K, J ).NE.ZERO )THEN
                TEMP      = ALPHA*B( K, J )
                B( K, J ) = TEMP
                IF( NOUNIT ) B( K, J ) = B( K, J )*A( K, K )
                DO I = K + 1, M
                  B( I, J ) = B( I, J ) + TEMP*A( I, K )
                end do
              END IF
            end do
          end do
        END IF
      ELSE
        !
        !           Form  B := alpha*A'*B   or   B := alpha*conjg( A' )*B.
        !
        IF( UPPER )THEN
          DO J = 1, N
            DO I = M, 1, -1
              TEMP = B( I, J )
              IF( NOCONJ )THEN
                IF( NOUNIT ) TEMP = TEMP*A( I, I )
                DO K = 1, I - 1
                  TEMP = TEMP + A( K, I )*B( K, J )
                end do
              ELSE
                IF( NOUNIT ) TEMP = TEMP*DCONJG( A( I, I ) )
                DO K = 1, I - 1
                  TEMP = TEMP + DCONJG( A( K, I ) )*B( K, J )
                end do
              END IF
              B( I, J ) = ALPHA*TEMP
            end do
          end do
        ELSE
          DO J = 1, N
            DO I = 1, M
              TEMP = B( I, J )
              IF( NOCONJ )THEN
                IF( NOUNIT )  TEMP = TEMP*A( I, I )
                DO K = I + 1, M
                  TEMP = TEMP + A( K, I )*B( K, J )
                end do
              ELSE
                IF( NOUNIT ) TEMP = TEMP*DCONJG( A( I, I ) )
                DO K = I + 1, M
                  TEMP = TEMP + DCONJG( A( K, I ) )*B( K, J )
                end do
              END IF
              B( I, J ) = ALPHA*TEMP
            end do
          end do
        END IF
      END IF
    ELSE
      IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*B*A.
        !
        IF( UPPER )THEN
          DO J = N, 1, -1
            TEMP = ALPHA
            IF( NOUNIT )    TEMP = TEMP*A( J, J )
            DO I = 1, M
              B( I, J ) = TEMP*B( I, J )
            end do
            DO K = 1, J - 1
              IF( A( K, J ).NE.ZERO )THEN
                TEMP = ALPHA*A( K, J )
                DO I = 1, M
                  B( I, J ) = B( I, J ) + TEMP*B( I, K )
                end do
              END IF
            end do
          end do
        ELSE
          DO J = 1, N
            TEMP = ALPHA
            IF( NOUNIT )   TEMP = TEMP*A( J, J )
            DO I = 1, M
              B( I, J ) = TEMP*B( I, J )
            end do
            DO  K = J + 1, N
              IF( A( K, J ).NE.ZERO )THEN
                TEMP = ALPHA*A( K, J )
                DO I = 1, M
                  B( I, J ) = B( I, J ) + TEMP*B( I, K )
                end do
              END IF
            end do

          end do
        END IF
      ELSE
        !
        !           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
        !
        IF( UPPER )THEN
          DO K = 1, N
            DO J = 1, K - 1
              IF( A( J, K ).NE.ZERO )THEN
                IF( NOCONJ )THEN
                  TEMP = ALPHA*A( J, K )
                ELSE
                  TEMP = ALPHA*DCONJG( A( J, K ) )
                END IF
                DO I = 1, M
                  B( I, J ) = B( I, J ) + TEMP*B( I, K )
                end do
              END IF
            end do
            TEMP = ALPHA
            IF( NOUNIT )THEN
              IF( NOCONJ )THEN
                TEMP = TEMP*A( K, K )
              ELSE
                TEMP = TEMP*DCONJG( A( K, K ) )
              END IF
            END IF
            IF( TEMP.NE.ONE )THEN
              DO I = 1, M
                B( I, K ) = TEMP*B( I, K )
              end do
            END IF
          end do
        ELSE
          DO K = N, 1, -1
            DO J = K + 1, N
              IF( A( J, K ).NE.ZERO )THEN
                IF( NOCONJ )THEN
                  TEMP = ALPHA*A( J, K )
                ELSE
                  TEMP = ALPHA*DCONJG( A( J, K ) )
                END IF
                DO I = 1, M
                  B( I, J ) = B( I, J ) + TEMP*B( I, K )
                end do
              END IF
            end do
            TEMP = ALPHA
            IF( NOUNIT )THEN
              IF( NOCONJ )THEN
                TEMP = TEMP*A( K, K )
              ELSE
                TEMP = TEMP*DCONJG( A( K, K ) )
              END IF
            END IF
            IF( TEMP.NE.ONE )THEN
              DO I = 1, M
                B( I, K ) = TEMP*B( I, K )
              end do
            END IF
          end do
        END IF
      END IF
    END IF
    !
    RETURN
    !
    !     End of ZTRMM .
    !
  END SUBROUTINE ZTRMM



  SUBROUTINE ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
    !     .. Scalar Arguments ..
    COMPLEX*16         ALPHA, BETA
    INTEGER            INCX, INCY, LDA, M, N
    CHARACTER*1        TRANS
    !     .. Array Arguments ..
    COMPLEX*16         A( LDA, * ), X( * ), Y( * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZGEMV  performs one of the matrix-vector operations
    !
    !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
    !
    !     y := alpha*conjg( A' )*x + beta*y,
    !
    !  where alpha and beta are scalars, x and y are vectors and A is an
    !  m by n matrix.
    !
    !  Parameters
    !  ==========
    !
    !  TRANS  - CHARACTER*1.
    !           On entry, TRANS specifies the operation to be performed as
    !           follows:
    !
    !              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
    !
    !              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
    !
    !              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
    !
    !           Unchanged on exit.
    !
    !  M      - INTEGER.
    !           On entry, M specifies the number of rows of the matrix A.
    !           M must be at least zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry, N specifies the number of columns of the matrix A.
    !           N must be at least zero.
    !           Unchanged on exit.
    !
    !  ALPHA  - COMPLEX*16      .
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
    !           Before entry, the leading m by n part of the array A must
    !           contain the matrix of coefficients.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. LDA must be at least
    !           max( 1, m ).
    !           Unchanged on exit.
    !
    !  X      - COMPLEX*16       array of DIMENSION at least
    !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
    !           and at least
    !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
    !           Before entry, the incremented array X must contain the
    !           vector x.
    !           Unchanged on exit.
    !
    !  INCX   - INTEGER.
    !           On entry, INCX specifies the increment for the elements of
    !           X. INCX must not be zero.
    !           Unchanged on exit.
    !
    !  BETA   - COMPLEX*16      .
    !           On entry, BETA specifies the scalar beta. When BETA is
    !           supplied as zero then Y need not be set on input.
    !           Unchanged on exit.
    !
    !  Y      - COMPLEX*16       array of DIMENSION at least
    !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
    !           and at least
    !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
    !           Before entry with BETA non-zero, the incremented array Y
    !           must contain the vector y. On exit, Y is overwritten by the
    !           updated vector y.
    !
    !  INCY   - INTEGER.
    !           On entry, INCY specifies the increment for the elements of
    !           Y. INCY must not be zero.
    !           Unchanged on exit.
    !
    !
    !  Level 2 Blas routine.
    !
    !  -- Written on 22-October-1986.
    !     Jack Dongarra, Argonne National Lab.
    !     Jeremy Du Croz, Nag Central Office.
    !     Sven Hammarling, Nag Central Office.
    !     Richard Hanson, Sandia National Labs.
    !
    !
    !     .. Parameters ..
    COMPLEX*16         ONE
    PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
    COMPLEX*16         ZERO
    PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
    !     .. Local Scalars ..
    COMPLEX*16         TEMP
    INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
    LOGICAL            NOCONJ
    !     .. External Functions ..
    !LOGICAL            LSAME
    !EXTERNAL           LSAME
    !     .. External Subroutines ..
    !EXTERNAL           XERBLA
    !     .. Intrinsic Functions ..
    INTRINSIC          DCONJG, MAX
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
         .NOT.LSAME( TRANS, 'T' ).AND. &
         .NOT.LSAME( TRANS, 'C' )      )THEN
      INFO = 1
    ELSE IF( M.LT.0 )THEN
      INFO = 2
    ELSE IF( N.LT.0 )THEN
      INFO = 3
    ELSE IF( LDA.LT.MAX( 1, M ) )THEN
      INFO = 6
    ELSE IF( INCX.EQ.0 )THEN
      INFO = 8
    ELSE IF( INCY.EQ.0 )THEN
      INFO = 11
    END IF
    IF( INFO.NE.0 )THEN
      CALL XERBLA( 'ZGEMV ', INFO )
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
         ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
         RETURN
    !
    NOCONJ = LSAME( TRANS, 'T' )
    !
    !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    !     up the start points in  X  and  Y.
    !
    IF( LSAME( TRANS, 'N' ) )THEN
      LENX = N
      LENY = M
    ELSE
      LENX = M
      LENY = N
    END IF
    IF( INCX.GT.0 )THEN
      KX = 1
    ELSE
      KX = 1 - ( LENX - 1 )*INCX
    END IF
    IF( INCY.GT.0 )THEN
      KY = 1
    ELSE
      KY = 1 - ( LENY - 1 )*INCY
    END IF
    !
    !     Start the operations. In this version the elements of A are
    !     accessed sequentially with one pass through A.
    !
    !     First form  y := beta*y.
    !
    IF( BETA.NE.ONE )THEN
      IF( INCY.EQ.1 )THEN
        IF( BETA.EQ.ZERO )THEN
          DO I = 1, LENY
            Y( I ) = ZERO
          end do
        ELSE
          DO I = 1, LENY
            Y( I ) = BETA*Y( I )
          end do
        END IF
      ELSE
        IY = KY
        IF( BETA.EQ.ZERO )THEN
          DO I = 1, LENY
            Y( IY ) = ZERO
            IY      = IY   + INCY
          end do
        ELSE
          DO I = 1, LENY
            Y( IY ) = BETA*Y( IY )
            IY      = IY           + INCY
          end do
        END IF
      END IF
    END IF
    IF( ALPHA.EQ.ZERO )   RETURN
    IF( LSAME( TRANS, 'N' ) )THEN
      !
      !        Form  y := alpha*A*x + y.
      !
      JX = KX
      IF( INCY.EQ.1 )THEN
        DO J = 1, N
          IF( X( JX ).NE.ZERO )THEN
            TEMP = ALPHA*X( JX )
            DO I = 1, M
              Y( I ) = Y( I ) + TEMP*A( I, J )
            end do
          END IF
          JX = JX + INCX
        end do
      ELSE
        DO J = 1, N
          IF( X( JX ).NE.ZERO )THEN
            TEMP = ALPHA*X( JX )
            IY   = KY
            DO I = 1, M
              Y( IY ) = Y( IY ) + TEMP*A( I, J )
              IY      = IY      + INCY
            end do
          END IF
          JX = JX + INCX
        end do
      END IF
    ELSE
      !
      !        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
      !
      JY = KY
      IF( INCX.EQ.1 )THEN
        DO J = 1, N
          TEMP = ZERO
          IF( NOCONJ )THEN
            DO I = 1, M
              TEMP = TEMP + A( I, J )*X( I )
            end do
          ELSE
            DO I = 1, M
              TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
            end do
          END IF
          Y( JY ) = Y( JY ) + ALPHA*TEMP
          JY      = JY      + INCY
        end do
      ELSE
        DO J = 1, N
          TEMP = ZERO
          IX   = KX
          IF( NOCONJ )THEN
            DO I = 1, M
              TEMP = TEMP + A( I, J )*X( IX )
              IX   = IX   + INCX
            end do
          ELSE
            DO I = 1, M
              TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
              IX   = IX   + INCX
            end do
          END IF
          Y( JY ) = Y( JY ) + ALPHA*TEMP
          JY      = JY      + INCY
        end do
      END IF
    END IF
    !
    RETURN
    !
    !     End of ZGEMV .
    !
  END SUBROUTINE ZGEMV


end module LAPACK_Z16_tools
