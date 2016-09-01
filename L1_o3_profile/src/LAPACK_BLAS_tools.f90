!> BLAS routines used by LAPACK modules
module lapack_blas_tools

  use LAPACK_tools, only: ilaenv, lsame
  use gsvd_o3prof_utilities, only: dtrti2, dtrtri
  use m_xerbla

  public daxpy, dcopy, dgemm, dgemv, dger, dscal, dswap, dtbsv, dtrsm, &
       dtrsv, drot, dasum, ddot, idamax
  private

contains

  SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
    !     .. Scalar Arguments ..
    DOUBLE PRECISION DA
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     constant times a vector plus a vector.
    !     uses unrolled loops for increments equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    INTEGER I,IX,IY,M,MP1
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !     ..
    IF (N.LE.0) RETURN
    IF (DA.EQ.0.0d0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    !
    !        code for unequal increments or equal increments
    !          not equal to 1
    !
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DY(IY) = DY(IY) + DA*DX(IX)
      IX = IX + INCX
      IY = IY + INCY
    end do
    RETURN
    !
    !        code for both increments equal to 1
    !
    !
    !        clean-up loop
    !
20  M = MOD(N,4)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
      DY(I) = DY(I) + DA*DX(I)
    end do
    IF (N.LT.4) RETURN
40  MP1 = M + 1
    DO I = MP1,N,4
      DY(I) = DY(I) + DA*DX(I)
      DY(I+1) = DY(I+1) + DA*DX(I+1)
      DY(I+2) = DY(I+2) + DA*DX(I+2)
      DY(I+3) = DY(I+3) + DA*DX(I+3)
    end do
    RETURN
  END SUBROUTINE DAXPY




  SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
    !     .. Scalar Arguments ..
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     copies a vector, x, to a vector, y.
    !     uses unrolled loops for increments equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    INTEGER I,IX,IY,M,MP1
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !     ..
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    !
    !        code for unequal increments or equal increments
    !          not equal to 1
    !
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DY(IY) = DX(IX)
      IX = IX + INCX
      IY = IY + INCY
    end do
    RETURN
    !
    !        code for both increments equal to 1
    !
    !
    !        clean-up loop
    !
20  M = MOD(N,7)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
      DY(I) = DX(I)
    end do
    IF (N.LT.7) RETURN
40  MP1 = M + 1
    DO I = MP1,N,7
      DY(I) = DX(I)
      DY(I+1) = DX(I+1)
      DY(I+2) = DX(I+2)
      DY(I+3) = DX(I+3)
      DY(I+4) = DX(I+4)
      DY(I+5) = DX(I+5)
      DY(I+6) = DX(I+6)
    end do
    RETURN
  END SUBROUTINE DCOPY




  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    use m_xerbla, only: xerbla

    implicit none

    !     .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA,BETA
    INTEGER K,LDA,LDB,LDC,M,N
    CHARACTER TRANSA,TRANSB
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DGEMM  performs one of the matrix-matrix operations
    !
    !     C := alpha*op( A )*op( B ) + beta*C,
    !
    !  where  op( X ) is one of
    !
    !     op( X ) = X   or   op( X ) = X',
    !
    !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
    !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
    !
    !  Arguments
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
    !              TRANSA = 'C' or 'c',  op( A ) = A'.
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
    !              TRANSB = 'C' or 'c',  op( B ) = B'.
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
    !  ALPHA  - DOUBLE PRECISION.
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
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
    !  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
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
    !  BETA   - DOUBLE PRECISION.
    !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
    !           supplied as zero then C need not be set on input.
    !           Unchanged on exit.
    !
    !  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
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
    !LOGICAL LSAME
    !EXTERNAL LSAME
    !     ..
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MAX
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
    LOGICAL NOTA,NOTB
    !     ..
    !     .. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    !     ..
    !
    !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    !     and  columns of  A  and the  number of  rows  of  B  respectively.
    !
    NOTA = LSAME(TRANSA,'N')
    NOTB = LSAME(TRANSB,'N')
    IF (NOTA) THEN
      NROWA = M
      NCOLA = K
    ELSE
      NROWA = K
      NCOLA = M
    END IF
    IF (NOTB) THEN
      NROWB = K
    ELSE
      NROWB = N
    END IF
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
         (.NOT.LSAME(TRANSA,'T'))) THEN
      INFO = 1
    ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
         (.NOT.LSAME(TRANSB,'T'))) THEN
      INFO = 2
    ELSE IF (M.LT.0) THEN
      INFO = 3
    ELSE IF (N.LT.0) THEN
      INFO = 4
    ELSE IF (K.LT.0) THEN
      INFO = 5
    ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
      INFO = 8
    ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
      INFO = 10
    ELSE IF (LDC.LT.MAX(1,M)) THEN
      INFO = 13
    END IF
    IF (INFO.NE.0) THEN
      CALL XERBLA('DGEMM ',INFO)
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
         (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
    !
    !     And if  alpha.eq.zero.
    !
    IF (ALPHA.EQ.ZERO) THEN
      IF (BETA.EQ.ZERO) THEN
        DO J = 1,N
          DO I = 1,M
            C(I,J) = ZERO
          end do
        end do
      ELSE
        DO J = 1,N
          DO I = 1,M
            C(I,J) = BETA*C(I,J)
          end do
        end do
      END IF
      RETURN
    END IF
    !
    !     Start the operations.
    !
    IF (NOTB) THEN
      IF (NOTA) THEN
        !
        !           Form  C := alpha*A*B + beta*C.
        !
        DO J = 1,N
          IF (BETA.EQ.ZERO) THEN
            DO I = 1,M
              C(I,J) = ZERO
            end do
          ELSE IF (BETA.NE.ONE) THEN
            DO I = 1,M
              C(I,J) = BETA*C(I,J)
            end do
          END IF
          DO L = 1,K
            IF (B(L,J).NE.ZERO) THEN
              TEMP = ALPHA*B(L,J)
              DO I = 1,M
                C(I,J) = C(I,J) + TEMP*A(I,L)
              end do
            END IF
          end do
        end do
      ELSE
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        DO J = 1,N
          DO I = 1,M
            TEMP = ZERO
            DO L = 1,K
              TEMP = TEMP + A(L,I)*B(L,J)
            end do
            IF (BETA.EQ.ZERO) THEN
              C(I,J) = ALPHA*TEMP
            ELSE
              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
            END IF
          end do
        end do
      END IF
    ELSE
      IF (NOTA) THEN
        !
        !           Form  C := alpha*A*B' + beta*C
        !
        DO J = 1,N
          IF (BETA.EQ.ZERO) THEN
            DO I = 1,M
              C(I,J) = ZERO
            end do
          ELSE IF (BETA.NE.ONE) THEN
            DO I = 1,M
              C(I,J) = BETA*C(I,J)
            end do
          END IF
          DO L = 1,K
            IF (B(J,L).NE.ZERO) THEN
              TEMP = ALPHA*B(J,L)
              DO I = 1,M
                C(I,J) = C(I,J) + TEMP*A(I,L)
              end do
            END IF
          end do
        end do
      ELSE
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        DO J = 1,N
          DO I = 1,M
            TEMP = ZERO
            DO L = 1,K
              TEMP = TEMP + A(L,I)*B(J,L)
            end do
            IF (BETA.EQ.ZERO) THEN
              C(I,J) = ALPHA*TEMP
            ELSE
              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
            END IF
          end do
        end do

      END IF
    END IF
    !
    RETURN
    !
    !     End of DGEMM .
    !
  END SUBROUTINE DGEMM


  SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

    use m_xerbla, only: xerbla

    implicit none
    !     .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA,BETA
    INTEGER INCX,INCY,LDA,M,N
    CHARACTER TRANS
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*),Y(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DGEMV  performs one of the matrix-vector operations
    !
    !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    !
    !  where alpha and beta are scalars, x and y are vectors and A is an
    !  m by n matrix.
    !
    !  Arguments
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
    !              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
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
    !  ALPHA  - DOUBLE PRECISION.
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
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
    !  X      - DOUBLE PRECISION array of DIMENSION at least
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
    !  BETA   - DOUBLE PRECISION.
    !           On entry, BETA specifies the scalar beta. When BETA is
    !           supplied as zero then Y need not be set on input.
    !           Unchanged on exit.
    !
    !  Y      - DOUBLE PRECISION array of DIMENSION at least
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
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
    !     ..
    !     .. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    !     ..
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MAX
    !     ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
         .NOT.LSAME(TRANS,'C')) THEN
      INFO = 1
    ELSE IF (M.LT.0) THEN
      INFO = 2
    ELSE IF (N.LT.0) THEN
      INFO = 3
    ELSE IF (LDA.LT.MAX(1,M)) THEN
      INFO = 6
    ELSE IF (INCX.EQ.0) THEN
      INFO = 8
    ELSE IF (INCY.EQ.0) THEN
      INFO = 11
    END IF
    IF (INFO.NE.0) THEN
      CALL XERBLA('DGEMV ',INFO)
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
         ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
    !
    !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    !     up the start points in  X  and  Y.
    !
    IF (LSAME(TRANS,'N')) THEN
      LENX = N
      LENY = M
    ELSE
      LENX = M
      LENY = N
    END IF
    IF (INCX.GT.0) THEN
      KX = 1
    ELSE
      KX = 1 - (LENX-1)*INCX
    END IF
    IF (INCY.GT.0) THEN
      KY = 1
    ELSE
      KY = 1 - (LENY-1)*INCY
    END IF
    !
    !     Start the operations. In this version the elements of A are
    !     accessed sequentially with one pass through A.
    !
    !     First form  y := beta*y.
    !
    IF (BETA.NE.ONE) THEN
      IF (INCY.EQ.1) THEN
        IF (BETA.EQ.ZERO) THEN
          DO I = 1,LENY
            Y(I) = ZERO
          end do
        ELSE
          DO I = 1,LENY
            Y(I) = BETA*Y(I)
          end do
        END IF
      ELSE
        IY = KY
        IF (BETA.EQ.ZERO) THEN
          DO I = 1,LENY
            Y(IY) = ZERO
            IY = IY + INCY
          end do
        ELSE
          DO I = 1,LENY
            Y(IY) = BETA*Y(IY)
            IY = IY + INCY
          end do
        END IF
      END IF
    END IF
    IF (ALPHA.EQ.ZERO) RETURN
    IF (LSAME(TRANS,'N')) THEN
      !
      !        Form  y := alpha*A*x + y.
      !
      JX = KX
      IF (INCY.EQ.1) THEN
        DO J = 1,N
          IF (X(JX).NE.ZERO) THEN
            TEMP = ALPHA*X(JX)
            DO I = 1,M
              Y(I) = Y(I) + TEMP*A(I,J)
            end do
          END IF
          JX = JX + INCX
        end do
      ELSE
        DO J = 1,N
          IF (X(JX).NE.ZERO) THEN
            TEMP = ALPHA*X(JX)
            IY = KY
            DO I = 1,M
              Y(IY) = Y(IY) + TEMP*A(I,J)
              IY = IY + INCY
            end do
          END IF
          JX = JX + INCX
        end do
      END IF
    ELSE
      !
      !        Form  y := alpha*A'*x + y.
      !
      JY = KY
      IF (INCX.EQ.1) THEN
        DO J = 1,N
          TEMP = ZERO
          DO I = 1,M
            TEMP = TEMP + A(I,J)*X(I)
          end do
          Y(JY) = Y(JY) + ALPHA*TEMP
          JY = JY + INCY
        end do
      ELSE
        DO J = 1,N
          TEMP = ZERO
          IX = KX
          DO I = 1,M
            TEMP = TEMP + A(I,J)*X(IX)
            IX = IX + INCX
          end do
          Y(JY) = Y(JY) + ALPHA*TEMP
          JY = JY + INCY
        end do
      END IF
    END IF
    !
    RETURN
    !
    !     End of DGEMV .
    !
  END SUBROUTINE DGEMV


  SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)

    use m_xerbla, only: xerbla

    implicit none
    !     .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA
    INTEGER INCX,INCY,LDA,M,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*),Y(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DGER   performs the rank 1 operation
    !
    !     A := alpha*x*y' + A,
    !
    !  where alpha is a scalar, x is an m element vector, y is an n element
    !  vector and A is an m by n matrix.
    !
    !  Arguments
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
    !  ALPHA  - DOUBLE PRECISION.
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  X      - DOUBLE PRECISION array of dimension at least
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
    !  Y      - DOUBLE PRECISION array of dimension at least
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
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
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
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D+0)
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,IX,J,JY,KX
    !     ..
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MAX
    !     ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF (M.LT.0) THEN
      INFO = 1
    ELSE IF (N.LT.0) THEN
      INFO = 2
    ELSE IF (INCX.EQ.0) THEN
      INFO = 5
    ELSE IF (INCY.EQ.0) THEN
      INFO = 7
    ELSE IF (LDA.LT.MAX(1,M)) THEN
      INFO = 9
    END IF
    IF (INFO.NE.0) THEN
      CALL XERBLA('DGER  ',INFO)
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
    !
    !     Start the operations. In this version the elements of A are
    !     accessed sequentially with one pass through A.
    !
    IF (INCY.GT.0) THEN
      JY = 1
    ELSE
      JY = 1 - (N-1)*INCY
    END IF
    IF (INCX.EQ.1) THEN
      DO J = 1,N
        IF (Y(JY).NE.ZERO) THEN
          TEMP = ALPHA*Y(JY)
          DO I = 1,M
            A(I,J) = A(I,J) + X(I)*TEMP
          end do
        END IF
        JY = JY + INCY
      end do
    ELSE
      IF (INCX.GT.0) THEN
        KX = 1
      ELSE
        KX = 1 - (M-1)*INCX
      END IF
      DO J = 1,N
        IF (Y(JY).NE.ZERO) THEN
          TEMP = ALPHA*Y(JY)
          IX = KX
          DO I = 1,M
            A(I,J) = A(I,J) + X(IX)*TEMP
            IX = IX + INCX
          end do
        END IF
        JY = JY + INCY
      end do
    END IF
    !
    RETURN
    !
    !     End of DGER  .
    !
  END SUBROUTINE DGER


  SUBROUTINE DSCAL(N,DA,DX,INCX)
    !     .. Scalar Arguments ..
    DOUBLE PRECISION DA
    INTEGER INCX,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*)
    !     ..
    !
    !  Purpose
    !  =======
    !*
    !     scales a vector by a constant.
    !     uses unrolled loops for increment equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 3/93 to return if incx .le. 0.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    INTEGER I,M,MP1,NINCX
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !     ..
    IF (N.LE.0 .OR. INCX.LE.0) RETURN
    IF (INCX.EQ.1) GO TO 20
    !
    !        code for increment not equal to 1
    !
    NINCX = N*INCX
    DO I = 1,NINCX,INCX
      DX(I) = DA*DX(I)
    end do
    RETURN
    !
    !        code for increment equal to 1
    !
    !
    !        clean-up loop
    !
20  M = MOD(N,5)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
      DX(I) = DA*DX(I)
    end do
    IF (N.LT.5) RETURN
40  MP1 = M + 1
    DO I = MP1,N,5
      DX(I) = DA*DX(I)
      DX(I+1) = DA*DX(I+1)
      DX(I+2) = DA*DX(I+2)
      DX(I+3) = DA*DX(I+3)
      DX(I+4) = DA*DX(I+4)
    end do
    RETURN
  END SUBROUTINE DSCAL




  SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
    !     .. Scalar Arguments ..
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     interchanges two vectors.
    !     uses unrolled loops for increments equal one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY,M,MP1
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !     ..
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    !
    !       code for unequal increments or equal increments not equal
    !         to 1
    !
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DTEMP = DX(IX)
      DX(IX) = DY(IY)
      DY(IY) = DTEMP
      IX = IX + INCX
      IY = IY + INCY
    end do
    RETURN
    !
    !       code for both increments equal to 1
    !
    !
    !       clean-up loop
    !
20  M = MOD(N,3)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
      DTEMP = DX(I)
      DX(I) = DY(I)
      DY(I) = DTEMP
    end do
    IF (N.LT.3) RETURN
40  MP1 = M + 1
    DO I = MP1,N,3
      DTEMP = DX(I)
      DX(I) = DY(I)
      DY(I) = DTEMP
      DTEMP = DX(I+1)
      DX(I+1) = DY(I+1)
      DY(I+1) = DTEMP
      DTEMP = DX(I+2)
      DX(I+2) = DY(I+2)
      DY(I+2) = DTEMP
    end do
    RETURN
  END SUBROUTINE DSWAP




  SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)

    use m_xerbla, only: xerbla

    implicit none
    !     .. Scalar Arguments ..
    INTEGER INCX,K,LDA,N
    CHARACTER DIAG,TRANS,UPLO
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DTBSV  solves one of the systems of equations
    !
    !     A*x = b,   or   A'*x = b,
    !
    !  where b and x are n element vectors and A is an n by n unit, or
    !  non-unit, upper or lower triangular band matrix, with ( k + 1 )
    !  diagonals.
    !
    !  No test for singularity or near-singularity is included in this
    !  routine. Such tests must be performed before calling this routine.
    !
    !  Arguments
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
    !           On entry, TRANS specifies the equations to be solved as
    !           follows:
    !
    !              TRANS = 'N' or 'n'   A*x = b.
    !
    !              TRANS = 'T' or 't'   A'*x = b.
    !
    !              TRANS = 'C' or 'c'   A'*x = b.
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
    !  K      - INTEGER.
    !           On entry with UPLO = 'U' or 'u', K specifies the number of
    !           super-diagonals of the matrix A.
    !           On entry with UPLO = 'L' or 'l', K specifies the number of
    !           sub-diagonals of the matrix A.
    !           K must satisfy  0 .le. K.
    !           Unchanged on exit.
    !
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
    !           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
    !           by n part of the array A must contain the upper triangular
    !           band part of the matrix of coefficients, supplied column by
    !           column, with the leading diagonal of the matrix in row
    !           ( k + 1 ) of the array, the first super-diagonal starting at
    !           position 2 in row k, and so on. The top left k by k triangle
    !           of the array A is not referenced.
    !           The following program segment will transfer an upper
    !           triangular band matrix from conventional full matrix storage
    !           to band storage:
    !
    !                 DO 20, J = 1, N
    !                    M = K + 1 - J
    !                    DO 10, I = MAX( 1, J - K ), J
    !                       A( M + I, J ) = matrix( I, J )
    !              10    CONTINUE
    !              20 CONTINUE
    !
    !           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
    !           by n part of the array A must contain the lower triangular
    !           band part of the matrix of coefficients, supplied column by
    !           column, with the leading diagonal of the matrix in row 1 of
    !           the array, the first sub-diagonal starting at position 1 in
    !           row 2, and so on. The bottom right k by k triangle of the
    !           array A is not referenced.
    !           The following program segment will transfer a lower
    !           triangular band matrix from conventional full matrix storage
    !           to band storage:
    !
    !                 DO 20, J = 1, N
    !                    M = 1 - J
    !                    DO 10, I = J, MIN( N, J + K )
    !                       A( M + I, J ) = matrix( I, J )
    !              10    CONTINUE
    !              20 CONTINUE
    !
    !           Note that when DIAG = 'U' or 'u' the elements of the array A
    !           corresponding to the diagonal elements of the matrix are not
    !           referenced, but are assumed to be unity.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. LDA must be at least
    !           ( k + 1 ).
    !           Unchanged on exit.
    !
    !  X      - DOUBLE PRECISION array of dimension at least
    !           ( 1 + ( n - 1 )*abs( INCX ) ).
    !           Before entry, the incremented array X must contain the n
    !           element right-hand side vector b. On exit, X is overwritten
    !           with the solution vector x.
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
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D+0)
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
    LOGICAL NOUNIT
    !     ..
    !     .. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    !     ..
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MAX,MIN
    !     ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
      INFO = 1
    ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
         .NOT.LSAME(TRANS,'C')) THEN
      INFO = 2
    ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
      INFO = 3
    ELSE IF (N.LT.0) THEN
      INFO = 4
    ELSE IF (K.LT.0) THEN
      INFO = 5
    ELSE IF (LDA.LT. (K+1)) THEN
      INFO = 7
    ELSE IF (INCX.EQ.0) THEN
      INFO = 9
    END IF
    IF (INFO.NE.0) THEN
      CALL XERBLA('DTBSV ',INFO)
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF (N.EQ.0) RETURN
    !
    NOUNIT = LSAME(DIAG,'N')
    !
    !     Set up the start point in X if the increment is not unity. This
    !     will be  ( N - 1 )*INCX  too small for descending loops.
    !
    IF (INCX.LE.0) THEN
      KX = 1 - (N-1)*INCX
    ELSE IF (INCX.NE.1) THEN
      KX = 1
    END IF
    !
    !     Start the operations. In this version the elements of A are
    !     accessed by sequentially with one pass through A.
    !
    IF (LSAME(TRANS,'N')) THEN
      !
      !        Form  x := inv( A )*x.
      !
      IF (LSAME(UPLO,'U')) THEN
        KPLUS1 = K + 1
        IF (INCX.EQ.1) THEN
          DO J = N,1,-1
            IF (X(J).NE.ZERO) THEN
              L = KPLUS1 - J
              IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
              TEMP = X(J)
              DO I = J - 1,MAX(1,J-K),-1
                X(I) = X(I) - TEMP*A(L+I,J)
              end do
            END IF
          end do
        ELSE
          KX = KX + (N-1)*INCX
          JX = KX
          DO J = N,1,-1
            KX = KX - INCX
            IF (X(JX).NE.ZERO) THEN
              IX = KX
              L = KPLUS1 - J
              IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
              TEMP = X(JX)
              DO I = J - 1,MAX(1,J-K),-1
                X(IX) = X(IX) - TEMP*A(L+I,J)
                IX = IX - INCX
              end do
            END IF
            JX = JX - INCX
          end do
        END IF
      ELSE
        IF (INCX.EQ.1) THEN
          DO J = 1,N
            IF (X(J).NE.ZERO) THEN
              L = 1 - J
              IF (NOUNIT) X(J) = X(J)/A(1,J)
              TEMP = X(J)
              DO I = J + 1,MIN(N,J+K)
                X(I) = X(I) - TEMP*A(L+I,J)
              end do
            END IF
          end do
        ELSE
          JX = KX
          DO J = 1,N
            KX = KX + INCX
            IF (X(JX).NE.ZERO) THEN
              IX = KX
              L = 1 - J
              IF (NOUNIT) X(JX) = X(JX)/A(1,J)
              TEMP = X(JX)
              DO I = J + 1,MIN(N,J+K)
                X(IX) = X(IX) - TEMP*A(L+I,J)
                IX = IX + INCX
              end do
            END IF
            JX = JX + INCX
          end do
        END IF
      END IF
    ELSE
      !
      !        Form  x := inv( A')*x.
      !
      IF (LSAME(UPLO,'U')) THEN
        KPLUS1 = K + 1
        IF (INCX.EQ.1) THEN
          DO J = 1,N
            TEMP = X(J)
            L = KPLUS1 - J
            DO I = MAX(1,J-K),J - 1
              TEMP = TEMP - A(L+I,J)*X(I)
            end do
            IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
            X(J) = TEMP
          end do
        ELSE
          JX = KX
          DO J = 1,N
            TEMP = X(JX)
            IX = KX
            L = KPLUS1 - J
            DO I = MAX(1,J-K),J - 1
              TEMP = TEMP - A(L+I,J)*X(IX)
              IX = IX + INCX
            end do
            IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
            X(JX) = TEMP
            JX = JX + INCX
            IF (J.GT.K) KX = KX + INCX
          end do
        END IF
      ELSE
        IF (INCX.EQ.1) THEN
          DO J = N,1,-1
            TEMP = X(J)
            L = 1 - J
            DO I = MIN(N,J+K),J + 1,-1
              TEMP = TEMP - A(L+I,J)*X(I)
            end do
            IF (NOUNIT) TEMP = TEMP/A(1,J)
            X(J) = TEMP
          end do
        ELSE
          KX = KX + (N-1)*INCX
          JX = KX
          DO J = N,1,-1
            TEMP = X(JX)
            IX = KX
            L = 1 - J
            DO I = MIN(N,J+K),J + 1,-1
              TEMP = TEMP - A(L+I,J)*X(IX)
              IX = IX - INCX
            end do
            IF (NOUNIT) TEMP = TEMP/A(1,J)
            X(JX) = TEMP
            JX = JX - INCX
            IF ((N-J).GE.K) KX = KX - INCX
          end do
        END IF
      END IF
    END IF
    !
    RETURN
    !
    !     End of DTBSV .
    !
  END SUBROUTINE DTBSV


  SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

    use m_xerbla, only: xerbla

    implicit none

    !     .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA
    INTEGER LDA,LDB,M,N
    CHARACTER DIAG,SIDE,TRANSA,UPLO
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),B(LDB,*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DTRSM  solves one of the matrix equations
    !
    !     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
    !
    !  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
    !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
    !
    !     op( A ) = A   or   op( A ) = A'.
    !
    !  The matrix X is overwritten on B.
    !
    !  Arguments
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
    !              TRANSA = 'C' or 'c'   op( A ) = A'.
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
    !  ALPHA  - DOUBLE PRECISION.
    !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
    !           zero then  A is not referenced and  B need not be set before
    !           entry.
    !           Unchanged on exit.
    !
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
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
    !  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
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
    !
    !  -- Written on 8-February-1989.
    !     Jack Dongarra, Argonne National Laboratory.
    !     Iain Duff, AERE Harwell.
    !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !     Sven Hammarling, Numerical Algorithms Group Ltd.
    !
    !
    !     .. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    !     ..
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MAX
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,J,K,NROWA
    LOGICAL LSIDE,NOUNIT,UPPER
    !     ..
    !     .. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    !     ..
    !
    !     Test the input parameters.
    !
    LSIDE = LSAME(SIDE,'L')
    IF (LSIDE) THEN
      NROWA = M
    ELSE
      NROWA = N
    END IF
    NOUNIT = LSAME(DIAG,'N')
    UPPER = LSAME(UPLO,'U')
    !
    INFO = 0
    IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
      INFO = 1
    ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
      INFO = 2
    ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
         (.NOT.LSAME(TRANSA,'T')) .AND. &
         (.NOT.LSAME(TRANSA,'C'))) THEN
      INFO = 3
    ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
      INFO = 4
    ELSE IF (M.LT.0) THEN
      INFO = 5
    ELSE IF (N.LT.0) THEN
      INFO = 6
    ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
      INFO = 9
    ELSE IF (LDB.LT.MAX(1,M)) THEN
      INFO = 11
    END IF
    IF (INFO.NE.0) THEN
      CALL XERBLA('DTRSM ',INFO)
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF (N.EQ.0) RETURN
    !
    !     And when  alpha.eq.zero.
    !
    IF (ALPHA.EQ.ZERO) THEN
      DO J = 1,N
        DO I = 1,M
          B(I,J) = ZERO
        end do
      end do
      RETURN
    END IF
    !
    !     Start the operations.
    !
    IF (LSIDE) THEN
      IF (LSAME(TRANSA,'N')) THEN
        !
        !           Form  B := alpha*inv( A )*B.
        !
        IF (UPPER) THEN
          DO J = 1,N
            IF (ALPHA.NE.ONE) THEN
              DO I = 1,M
                B(I,J) = ALPHA*B(I,J)
              end do
            END IF
            DO K = M,1,-1
              IF (B(K,J).NE.ZERO) THEN
                IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                DO I = 1,K - 1
                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
                end do
              END IF
            end do
          end do
        ELSE
          DO J = 1,N
            IF (ALPHA.NE.ONE) THEN
              DO I = 1,M
                B(I,J) = ALPHA*B(I,J)
              end do
            END IF
            DO K = 1,M
              IF (B(K,J).NE.ZERO) THEN
                IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                DO I = K + 1,M
                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
                end do
              END IF
            end do
          end do
        END IF
      ELSE
        !
        !           Form  B := alpha*inv( A' )*B.
        !
        IF (UPPER) THEN
          DO J = 1,N
            DO I = 1,M
              TEMP = ALPHA*B(I,J)
              DO K = 1,I - 1
                TEMP = TEMP - A(K,I)*B(K,J)
              end do
              IF (NOUNIT) TEMP = TEMP/A(I,I)
              B(I,J) = TEMP
            end do
          end do
        ELSE
          DO J = 1,N
            DO I = M,1,-1
              TEMP = ALPHA*B(I,J)
              DO K = I + 1,M
                TEMP = TEMP - A(K,I)*B(K,J)
              end do
              IF (NOUNIT) TEMP = TEMP/A(I,I)
              B(I,J) = TEMP
            end do
          end do
        END IF
      END IF
    ELSE
      IF (LSAME(TRANSA,'N')) THEN
        !
        !           Form  B := alpha*B*inv( A ).
        !
        IF (UPPER) THEN
          DO J = 1,N
            IF (ALPHA.NE.ONE) THEN
              DO I = 1,M
                B(I,J) = ALPHA*B(I,J)
              end do
            END IF
            DO K = 1,J - 1
              IF (A(K,J).NE.ZERO) THEN
                DO I = 1,M
                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
                end do
              END IF
            end do
            IF (NOUNIT) THEN
              TEMP = ONE/A(J,J)
              DO I = 1,M
                B(I,J) = TEMP*B(I,J)
              end do
            END IF
          end do
        ELSE
          DO J = N,1,-1
            IF (ALPHA.NE.ONE) THEN
              DO I = 1,M
                B(I,J) = ALPHA*B(I,J)
              end do
            END IF
            DO K = J + 1,N
              IF (A(K,J).NE.ZERO) THEN
                DO I = 1,M
                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
                end do
              END IF
            end do
            IF (NOUNIT) THEN
              TEMP = ONE/A(J,J)
              DO I = 1,M
                B(I,J) = TEMP*B(I,J)
              end do
            END IF
          end do
        END IF
      ELSE
        !
        !           Form  B := alpha*B*inv( A' ).
        !
        IF (UPPER) THEN
          DO K = N,1,-1
            IF (NOUNIT) THEN
              TEMP = ONE/A(K,K)
              DO I = 1,M
                B(I,K) = TEMP*B(I,K)
              end do
            END IF
            DO J = 1,K - 1
              IF (A(J,K).NE.ZERO) THEN
                TEMP = A(J,K)
                DO I = 1,M
                  B(I,J) = B(I,J) - TEMP*B(I,K)
                end do
              END IF
            end do
            IF (ALPHA.NE.ONE) THEN
              DO I = 1,M
                B(I,K) = ALPHA*B(I,K)
              end do
            END IF
          end do
        ELSE
          DO K = 1,N
            IF (NOUNIT) THEN
              TEMP = ONE/A(K,K)
              DO I = 1,M
                B(I,K) = TEMP*B(I,K)
              end do
            END IF
            DO J = K + 1,N
              IF (A(J,K).NE.ZERO) THEN
                TEMP = A(J,K)
                DO I = 1,M
                  B(I,J) = B(I,J) - TEMP*B(I,K)
                end do
              END IF
            end do
            IF (ALPHA.NE.ONE) THEN
              DO I = 1,M
                B(I,K) = ALPHA*B(I,K)
              end do
            END IF
          end do
        END IF
      END IF
    END IF
    !
    RETURN
    !
    !     End of DTRSM .
    !
  END SUBROUTINE DTRSM


  SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)

    use m_xerbla, only: xerbla

    implicit none

    !     .. Scalar Arguments ..
    INTEGER INCX,LDA,N
    CHARACTER DIAG,TRANS,UPLO
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),X(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DTRSV  solves one of the systems of equations
    !
    !     A*x = b,   or   A'*x = b,
    !
    !  where b and x are n element vectors and A is an n by n unit, or
    !  non-unit, upper or lower triangular matrix.
    !
    !  No test for singularity or near-singularity is included in this
    !  routine. Such tests must be performed before calling this routine.
    !
    !  Arguments
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
    !           On entry, TRANS specifies the equations to be solved as
    !           follows:
    !
    !              TRANS = 'N' or 'n'   A*x = b.
    !
    !              TRANS = 'T' or 't'   A'*x = b.
    !
    !              TRANS = 'C' or 'c'   A'*x = b.
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
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
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
    !  X      - DOUBLE PRECISION array of dimension at least
    !           ( 1 + ( n - 1 )*abs( INCX ) ).
    !           Before entry, the incremented array X must contain the n
    !           element right-hand side vector b. On exit, X is overwritten
    !           with the solution vector x.
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
    DOUBLE PRECISION ZERO
    PARAMETER (ZERO=0.0D+0)
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,IX,J,JX,KX
    LOGICAL NOUNIT
    !     ..
    !     .. External Functions ..
    !LOGICAL LSAME
    !EXTERNAL LSAME
    !     ..
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MAX
    !     ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
      INFO = 1
    ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
         .NOT.LSAME(TRANS,'C')) THEN
      INFO = 2
    ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
      INFO = 3
    ELSE IF (N.LT.0) THEN
      INFO = 4
    ELSE IF (LDA.LT.MAX(1,N)) THEN
      INFO = 6
    ELSE IF (INCX.EQ.0) THEN
      INFO = 8
    END IF
    IF (INFO.NE.0) THEN
      CALL XERBLA('DTRSV ',INFO)
      RETURN
    END IF
    !
    !     Quick return if possible.
    !
    IF (N.EQ.0) RETURN
    !
    NOUNIT = LSAME(DIAG,'N')
    !
    !     Set up the start point in X if the increment is not unity. This
    !     will be  ( N - 1 )*INCX  too small for descending loops.
    !
    IF (INCX.LE.0) THEN
      KX = 1 - (N-1)*INCX
    ELSE IF (INCX.NE.1) THEN
      KX = 1
    END IF
    !
    !     Start the operations. In this version the elements of A are
    !     accessed sequentially with one pass through A.
    !
    IF (LSAME(TRANS,'N')) THEN
      !
      !        Form  x := inv( A )*x.
      !
      IF (LSAME(UPLO,'U')) THEN
        IF (INCX.EQ.1) THEN
          DO J = N,1,-1
            IF (X(J).NE.ZERO) THEN
              IF (NOUNIT) X(J) = X(J)/A(J,J)
              TEMP = X(J)
              DO I = J - 1,1,-1
                X(I) = X(I) - TEMP*A(I,J)
              end do
            END IF
          end do
        ELSE
          JX = KX + (N-1)*INCX
          DO J = N,1,-1
            IF (X(JX).NE.ZERO) THEN
              IF (NOUNIT) X(JX) = X(JX)/A(J,J)
              TEMP = X(JX)
              IX = JX
              DO I = J - 1,1,-1
                IX = IX - INCX
                X(IX) = X(IX) - TEMP*A(I,J)
              end do
            END IF
            JX = JX - INCX
          end do
        END IF
      ELSE
        IF (INCX.EQ.1) THEN
          DO  J = 1,N
            IF (X(J).NE.ZERO) THEN
              IF (NOUNIT) X(J) = X(J)/A(J,J)
              TEMP = X(J)
              DO I = J + 1,N
                X(I) = X(I) - TEMP*A(I,J)
              end do
            END IF
          end do
        ELSE
          JX = KX
          DO J = 1,N
            IF (X(JX).NE.ZERO) THEN
              IF (NOUNIT) X(JX) = X(JX)/A(J,J)
              TEMP = X(JX)
              IX = JX
              DO I = J + 1,N
                IX = IX + INCX
                X(IX) = X(IX) - TEMP*A(I,J)
              end do
            END IF
            JX = JX + INCX
          end do
        END IF
      END IF
    ELSE
      !
      !        Form  x := inv( A' )*x.
      !
      IF (LSAME(UPLO,'U')) THEN
        IF (INCX.EQ.1) THEN
          DO J = 1,N
            TEMP = X(J)
            DO I = 1,J - 1
              TEMP = TEMP - A(I,J)*X(I)
            end do
            IF (NOUNIT) TEMP = TEMP/A(J,J)
            X(J) = TEMP
          end do
        ELSE
          JX = KX
          DO J = 1,N
            TEMP = X(JX)
            IX = KX
            DO I = 1,J - 1
              TEMP = TEMP - A(I,J)*X(IX)
              IX = IX + INCX
            end do
            IF (NOUNIT) TEMP = TEMP/A(J,J)
            X(JX) = TEMP
            JX = JX + INCX
          end do
        END IF
      ELSE
        IF (INCX.EQ.1) THEN
          DO J = N,1,-1
            TEMP = X(J)
            DO I = N,J + 1,-1
              TEMP = TEMP - A(I,J)*X(I)
            end do
            IF (NOUNIT) TEMP = TEMP/A(J,J)
            X(J) = TEMP
          end do
        ELSE
          KX = KX + (N-1)*INCX
          JX = KX
          DO J = N,1,-1
            TEMP = X(JX)
            IX = KX
            DO I = N,J + 1,-1
              TEMP = TEMP - A(I,J)*X(IX)
              IX = IX - INCX
            end do
            IF (NOUNIT) TEMP = TEMP/A(J,J)
            X(JX) = TEMP
            JX = JX - INCX
          end do
        END IF
      END IF
    END IF
    !
    RETURN
    !
    !     End of DTRSV .
    !
  END SUBROUTINE DTRSV


  DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
    !     .. Scalar Arguments ..
    INTEGER INCX,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     takes the sum of the absolute values.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 3/93 to return if incx .le. 0.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,M,MP1,NINCX
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC DABS,MOD
    !     ..
    DASUM = 0.0d0
    DTEMP = 0.0d0
    IF (N.LE.0 .OR. INCX.LE.0) RETURN
    IF (INCX.EQ.1) GO TO 20
    !
    !        code for increment not equal to 1
    !
    NINCX = N*INCX
    DO  I = 1,NINCX,INCX
      DTEMP = DTEMP + DABS(DX(I))
    end do
    DASUM = DTEMP
    RETURN
    !
    !        code for increment equal to 1
    !
    !
    !        clean-up loop
    !
20  M = MOD(N,6)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
      DTEMP = DTEMP + DABS(DX(I))
    end do
    IF (N.LT.6) GO TO 60
40  MP1 = M + 1
    DO I = MP1,N,6
      DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) + &
           DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
    end do
60  DASUM = DTEMP
    RETURN
  END FUNCTION DASUM




  DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
    !     .. Scalar Arguments ..
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     forms the dot product of two vectors.
    !     uses unrolled loops for increments equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY,M,MP1
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !     ..
    DDOT = 0.0d0
    DTEMP = 0.0d0
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    !
    !        code for unequal increments or equal increments
    !          not equal to 1
    !
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DTEMP = DTEMP + DX(IX)*DY(IY)
      IX = IX + INCX
      IY = IY + INCY
    end do
    DDOT = DTEMP
    RETURN
    !
    !        code for both increments equal to 1
    !
    !
    !        clean-up loop
    !
20  M = MOD(N,5)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
      DTEMP = DTEMP + DX(I)*DY(I)
    end do
    IF (N.LT.5) GO TO 60
40  MP1 = M + 1
    DO  I = MP1,N,5
      DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
           DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
    end do
60  DDOT = DTEMP
    RETURN
  END FUNCTION DDOT




  INTEGER FUNCTION IDAMAX(N,DX,INCX)
    !     .. Scalar Arguments ..
    INTEGER INCX,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     finds the index of element having max. absolute value.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 3/93 to return if incx .le. 0.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    DOUBLE PRECISION DMAX
    INTEGER I,IX
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC DABS
    !     ..
    IDAMAX = 0
    IF (N.LT.1 .OR. INCX.LE.0) RETURN
    IDAMAX = 1
    IF (N.EQ.1) RETURN
    IF (INCX.EQ.1) GO TO 20
    !
    !        code for increment not equal to 1
    !
    IX = 1
    DMAX = DABS(DX(1))
    IX = IX + INCX
    DO I = 2,N
      IF (DABS(DX(IX)).LE.DMAX) GO TO 5
      IDAMAX = I
      DMAX = DABS(DX(IX))
5     IX = IX + INCX
    end do
    RETURN
    !
    !        code for increment equal to 1
    !
20  DMAX = DABS(DX(1))
    DO I = 2,N
      IF (DABS(DX(I)).LE.DMAX) CYCLE
      IDAMAX = I
      DMAX = DABS(DX(I))
    end do
    RETURN
  END FUNCTION IDAMAX


  SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
    !     .. Scalar Arguments ..
    DOUBLE PRECISION C,S
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(*),DY(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     applies a plane rotation.
    !     jack dongarra, linpack, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY
    !     ..
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    !
    !       code for unequal increments or equal increments not equal
    !         to 1
    !
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DTEMP = C*DX(IX) + S*DY(IY)
      DY(IY) = C*DY(IY) - S*DX(IX)
      DX(IX) = DTEMP
      IX = IX + INCX
      IY = IY + INCY
    end do
    RETURN
    !
    !       code for both increments equal to 1
    !
20  DO I = 1,N
      DTEMP = C*DX(I) + S*DY(I)
      DY(I) = C*DY(I) - S*DX(I)
      DX(I) = DTEMP
    end do
    RETURN
  END SUBROUTINE DROT


end module lapack_blas_tools
