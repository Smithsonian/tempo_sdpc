module m_invert

public invert

contains

function invert (amat, error) result (amatinv)

 !use m_die
 implicit none

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  invert
! 
! !DESCRIPTION: matrix inversion interface to look like IDL interface
!               calls routines from numerical recipes (may not be
!  fastest)
!
! !CALLING SEQUENCE: 
!
!        amatinv = invert(amat)
!     
! !INPUT PARAMETERS:   
       real (KIND=8), dimension(:,:), intent(in) :: amat
       integer,              intent(out):: error
!   amat : 2D matrix to invert
!
! !OUTPUT PARAMETERS:  
       real (KIND=8), dimension(lbound(amat,1):ubound(amat,1), &
                      lbound(amat,2):ubound(amat,2)) :: amatinv
!   amatinv : 2D inverted matrix
!
! !SEE ALSO:  numerical recipes
!
! !REVISION HISTORY: 
!
!  13Aug97   Joiner     fortran 90 version from AIRS code
!
!EOP
!-------------------------------------------------------------------------

       !integer, dimension (:), allocatable :: indx
       !real (KIND=8), dimension (:,:), allocatable :: work
       integer, dimension (size(amat,1)) :: indx
       real (KIND=8), dimension (size(amat,1),size(amat,2)) :: work
       real (KIND=8)              :: d
       integer           :: j
       integer           :: nsampl
!       integer           :: ierr

       nsampl=size(amat,1)
       if (nsampl /= size(amat,2)) then
         !call die('invert','can''t invert a non-square matrix')
         print *,'invert','can''t invert a non-square matrix'
         return
       endif
       amatinv=0.
       !allocate(indx(nsampl))
       !allocate(work(nsampl,nsampl),stat=ierr)
!       if (ierr /= 0) then
!         print *,'invert: tried to create ',nsampl,' square matrix'
         !call die('invert','can''t allocate enough memory')
!       endif
       work=amat
       do j = 0, nsampl-1
         amatinv(j+lbound(amat,1),j+lbound(amat,2))=1.
       end do
       call ludcmp(work, nsampl, nsampl, indx, d, error)
       if (error == 0) then
        do j = 0, nsampl-1
         call lubksb(work, nsampl, nsampl, indx, & 
             amatinv(:,j+lbound(amat,2)))
        end do
       endif
       !deallocate(indx)
       !deallocate(work)
       end function invert

      SUBROUTINE LUDCMP(A,N,NP,INDX,D, error)
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ludcmp
! 
! !DESCRIPTION: LU decomposition
!
! !CALLING SEQUENCE: 
!
!        call ludcmp(a,n,np,indx,d)
!     
! !INPUT PARAMETERS:   
!			a    : input matrix to decompose
!			n    : size of 2D square matrix a
!		 	np   : physical dimensions of a
!
! !OUTPUT PARAMETERS:  
!			indx : records row permutation
!			d    : +/-1 depending on whether # row 
!				interchanges even or odd
!
! !SEE ALSO:  numerical recipes, lubksb.f
!
! !REVISION HISTORY: 
!
!  02Mar96   Joiner     no changes to numerical recipes code
!
!EOP
!-------------------------------------------------------------------------
!

!      PARAMETER (NMAX=600,TINY=1.0E-20)
      integer n, np
      real (KIND=8) A(NP,NP)
      integer INDX(NP), error
      !real (KIND=8), dimension(:), allocatable :: VV
      real (KIND=8), dimension(n) :: VV
      real (KIND=8) d

      !allocate(vv(n))
      error=0
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) then
!          print *,'ludcmp: Singular matrix'
          error=1
          return
        endif
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      !deallocate(vv)
      RETURN
      END subroutine ludcmp

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  lubksb
! 
! !DESCRIPTION: solves set of linear equations Ax=b 
!
! !CALLING SEQUENCE: 
!
!        call lubksb(a,n,np,indx,b)
!     
! !INPUT PARAMETERS:   
!			a    : LU decomposition of A
!			indx : permutation vector from ludcmp
!			b    : RHS vector 
!			n    : size of b
!		 	np   : dimension of square matrix A		
!
! !OUTPUT PARAMETERS:  
!			b    : returns solution vector x 
!
! !SEE ALSO:  numerical recipes, ludcmp.f
!
! !REVISION HISTORY: 
!
!  02Mar96   Joiner     no changes to numerical recipes code
!
!EOP
!-------------------------------------------------------------------------
!

!!Added following parameter declarations since OMI code lacked them
      integer N, NP, I, II, J
      integer INDX(NP)
      real (KIND=8) A(NP,NP), B(NP)
!      DIMENSION A(NP,NP),INDX(NP),B(NP)


      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END subroutine lubksb

end module m_invert
