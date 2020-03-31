MODULE string_utils_module

USE parameters_module

IMPLICIT NONE

PUBLIC :: READ_ONE_LINE,  &
          SPLIT_ONE_LINE, &
          StrSplit,       &
          StrRepl,        &
          TxtExt,         &
          CntMat,         &
          CopyTxt,        &
          TRANUC,         &
          MatchListName
          
CONTAINS
    
    
      !###################################################################
      !#                              SPLAT                              #
      !###################################################################
      
      ! SUBROUTINE: SkipFileLines
      ! 
      ! DESCRIPTION: Skips the inputed number of lines from an open file
      
      SUBROUTINE SkipFileLines(funit,nLines)
            
            ! --------------------
            ! subroutine arguments
            ! --------------------
            INTEGER,  INTENT(IN)    :: funit
            INTEGER,  INTENT(IN)    :: nlines
            
            ! ---------------
            ! Local variables
            ! ---------------
            INTEGER                :: F
            CHARACTER(LEN=maxChar) :: tmpStr
            ! ============================================================
            ! SkipFileLines starts here
            ! ============================================================
            
            DO F=1,nLines
            READ(funit,'(A)') tmpStr
            ENDDO
            
      END SUBROUTINE SkipFileLines

      SUBROUTINE MatchListName( List, nList, Targ, nTarg, Idx )
        
        ! --------------------
        ! subroutine arguments
        ! --------------------
        CHARACTER(LEN=*), INTENT(IN)  :: List(nList)
        INTEGER,          INTENT(IN)  :: nList
        CHARACTER(LEN=*), INTENT(IN)  :: Targ(nTarg)
        INTEGER,          INTENT(IN)  :: nTarg
        INTEGER,          INTENT(OUT) :: Idx(nTarg)
        
        ! ---------------
        ! local variables
        ! ---------------
        INTEGER :: N_FOUND, N, L
        
        ! =====================================================================
        ! MatchGasProfileIndex starts here
        ! =====================================================================
        
        ! Initialize Idx
        Idx(:) = -1
        
        N_FOUND = 0
        
        DO N=1,nTarg
          DO L=1,nList
            IF(TRIM(ADJUSTL(List(L))) .EQ. TRIM(ADJUSTL(Targ(N)))) THEN
              Idx(N) = L
              N_FOUND=N_FOUND+1
            ENDIF
          ENDDO

          IF(N_FOUND .EQ. nTarg) EXIT
        ENDDO
        
      END SUBROUTINE MatchlistName

      SUBROUTINE MatchListOneName( List, nList, Targ, Idx )

        ! --------------------
        ! subroutine arguments
        ! --------------------
        CHARACTER(LEN=*), INTENT(IN)  :: List(nList)
        INTEGER,          INTENT(IN)  :: nList
        CHARACTER(LEN=*), INTENT(IN)  :: Targ
        INTEGER,          INTENT(OUT) :: Idx
        
        ! ---------------
        ! local variables
        ! ---------------
        CHARACTER(LEN=maxChar) :: Targ_in(1)
        INTEGER                :: Idx_out(1)

        ! =====================================================================
        ! MatchListOneName starts here
        ! =====================================================================
        
        Targ_in(1) = Targ
        CALL MatchListName(List, nList, Targ_in, 1, Idx_out)
        Idx = Idx_out(1)

      END SUBROUTINE MatchListOneName
    
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_one_line
!
! !DESCRIPTION: Subroutine READ\_ONE\_LINE reads a line from the input file.  
!  If the global variable VERBOSE is set, the line will be printed to stdout.  
!  READ\_ONE\_LINE can trap an unexpected EOF if LOCATION is passed.  
!  Otherwise, it will pass a logical flag back to the calling routine, 
!  where the error trapping will be done.
!\\
!\\
! !INTERFACE:
!
      FUNCTION READ_ONE_LINE( IU_GEOS, EOF, LOCATION ) RESULT( LINE )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!     
      INTEGER,          INTENT(IN)           :: IU_GEOS
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: LOCATION    ! Msg to display
!
! !OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(OUT)          :: EOF         ! Denotes EOF 
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Aug 2012 - R. Yantosca - Now make IU_GEOS a global module variable
!                              so that we can define it with findFreeLun
!  17 Sep 2013 - R. Yantosca - Extend line length to read in more tracers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IOS
!------------------------------------------------------------------------------
! Prior to 9/17/13:
! Need to extend the line length for many more tracers (bmy, 9/17/13)
!      CHARACTER(LEN=255) :: LINE, MSG
!------------------------------------------------------------------------------
      CHARACTER(LEN=maxChar) :: LINE, MSG

      !=================================================================
      ! READ_ONE_LINE begins here!
      !=================================================================

      ! Initialize
      EOF = .FALSE.

      ! Read a line from the file
      READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE

      ! IO Status < 0: EOF condition
      IF ( IOS < 0 ) THEN
         EOF = .TRUE.

         ! Trap unexpected EOF -- stop w/ error msg if LOCATION is passed
         ! Otherwise, return EOF to the calling program
         IF ( PRESENT( LOCATION ) ) THEN
            MSG = 'READ_ONE_LINE: error at: ' // TRIM( LOCATION )
            WRITE( 6, '(a)' ) MSG
            WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
            WRITE( 6, '(a)' ) 'STOP in READ_ONE_LINE (input_mod.f)'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            STOP
         ELSE
            RETURN
         ENDIF
      ENDIF

      ! IO Status > 0: true I/O error condition
!       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_one_line:1' )

      ! Print the line (if necessary)
!       IF ( VERBOSE ) WRITE( 6, '(a)' ) TRIM( LINE )

      END FUNCTION READ_ONE_LINE

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: split_one_line
!
! !DESCRIPTION: Subroutine SPLIT\_ONE\_LINE reads a line from the input file 
!  (via routine READ\_ONE\_LINE), and separates it into substrings.
!\\
!\\
!  SPLIT\_ONE\_LINE also checks to see if the number of substrings found is 
!  equal to the number of substrings that we expected to find.  However, if
!  you don't know a-priori how many substrings to expect a-priori, 
!  you can skip the error check.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SPLIT_ONE_LINE( funit, SUBSTRS, N_SUBSTRS, N_EXP, LOCATION ) 
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      ! Number of substrings we expect to find
      INTEGER,            INTENT(IN)  :: funit, N_EXP

      ! Name of routine that called SPLIT_ONE_LINE
      CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION 
!
! !OUTPUT PARAMETERS:
!
      ! Array of substrings (separated by " ")
      CHARACTER(LEN=maxChar), INTENT(OUT) :: SUBSTRS(maxChar)

      ! Number of substrings actually found
      INTEGER,            INTENT(OUT) :: N_SUBSTRS
      
      INTEGER, PARAMETER :: FIRSTCOL = 34

! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  17 Sep 2013 - R. Yantosca - Extend LINE to 500 chars to allow more tracers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                         :: EOF
!----------------------------------------------------------------------------
! Prior to 9/17/13:
! Extend LINE to 500 chars to allow more tracers (bmy, 9/17/13)
!      CHARACTER(LEN=255)              :: LINE, MSG
!----------------------------------------------------------------------------
      CHARACTER(LEN=maxChar)              :: LINE
      CHARACTER(LEN=maxChar)              :: MSG

      !=================================================================
      ! SPLIT_ONE_LINE begins here!
      !=================================================================      

      ! Create error msg
      MSG = 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )

      !=================================================================
      ! Read a line from disk
      !=================================================================
      LINE = READ_ONE_LINE( funit, EOF )

      ! STOP on End-of-File w/ error msg
      IF ( EOF ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, '(a)' ) 'End of file encountered!' 
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      !=================================================================
      ! Split the lines between spaces -- start at column FIRSTCOL
      !=================================================================
      CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

      ! Sometimes we don't know how many substrings to expect,
      ! if N_EXP is greater than MAXDIM, then skip the error check
      IF ( N_EXP < 0 ) RETURN

      ! Stop if we found the wrong 
      IF ( N_EXP /= N_SUBSTRS ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, 100   ) N_EXP, N_SUBSTRS
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
 100     FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
      ENDIF
       
      END SUBROUTINE SPLIT_ONE_LINE

!------------------------------------------------------------------------------

  SUBROUTINE StrSplit( STR, SEP, RESULT, N_SUBSTRS )
!
!******************************************************************************
!  Subroutine STRSPLIT returns substrings in a string, separated by a 
!  separator character (similar to IDL's StrSplit function).  This is mainly
!  a convenience wrapper for CHARPAK routine TxtExt. (bmy, 7/11/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STR       (CHARACTER*(*)) : String to be searched (variable length)  
!  (2 ) SEP       (CHARACTER*1  ) : Separator character
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) RESULT    (CHARACTER*255) : Array containing substrings (255 elements)
!  (4 ) N_SUBSTRS (INTEGER      ) : Number of substrings returned (optional)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)            :: STR
      CHARACTER(LEN=1), INTENT(IN)            :: SEP
      CHARACTER(LEN=*), INTENT(OUT)           :: RESULT(maxChar)
      INTEGER,          INTENT(OUT), OPTIONAL :: N_SUBSTRS

      ! Local variables
      INTEGER                                 :: I, IFLAG, COL
      CHARACTER (LEN=maxChar)                 :: WORD

      !=================================================================
      ! STRSPLIT begins here!
      !=================================================================

      ! Initialize
      I         = 0
      COL       = 1 
      IFLAG     = 0
      RESULT(:) = ''
      
      ! Loop until all matches found, or end of string
      DO WHILE ( IFLAG == 0 )

         ! Look for strings beteeen separator string
         CALL TXTEXT ( SEP, TRIM( STR ), COL, WORD, IFLAG )

         ! Store substrings in RESULT array
         I         = I + 1
         RESULT(I) = TRIM( WORD )

      ENDDO

      ! Optional argument: return # of substrings found
      IF ( PRESENT( N_SUBSTRS ) ) N_SUBSTRS = I

      ! Return to calling program
    END SUBROUTINE StrSplit

!------------------------------------------------------------------------------

      SUBROUTINE StrRepl( STR, PATTERN, REPLTXT )

      !=================================================================
      ! Subroutine STRREPL replaces all instances of PATTERN within
      ! a string STR with replacement text REPLTXT. 
      ! (bmy, 6/25/02, 7/20/04)
      !
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) STR     : String to be searched
      ! (2 ) PATTERN : Pattern of characters to replace w/in STR
      ! (3 ) REPLTXT : Replacement text for PATTERN
      !
      ! Arguments as Output:
      ! ----------------------------------------------------------------
      ! (1 ) STR     : String with new replacement text 
      !
      ! NOTES
      ! (1 ) REPLTXT must have the same # of characters as PATTERN.
      ! (2 ) Replace LEN_TRIM with LEN (bmy, 7/20/04)
      !=================================================================

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: STR
      CHARACTER(LEN=*), INTENT(IN)    :: PATTERN, REPLTXT
      
      ! Local variables
      INTEGER                         :: I1, I2

      !=================================================================
      ! STRREPL begins here!
      !=================================================================

      ! Error check: make sure PATTERN and REPLTXT have the same # of chars
      IF ( LEN( PATTERN ) /= LEN( REPLTXT ) ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'STRREPL: PATTERN and REPLTXT must have same # of characters!'
         WRITE( 6, '(a)' ) 'STOP in STRREPL (charpak_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      ! Loop over all instances of PATTERN in STR
      DO 

         ! I1 is the starting location of PATTERN w/in STR  
         I1 = INDEX( STR, PATTERN )

         ! If pattern is not found, then return to calling program
         IF ( I1 < 1 ) RETURN

         ! I2 is the ending location of PATTERN w/in STR
         I2 = I1 + LEN_TRIM( PATTERN ) - 1
      
         ! Replace text
         STR(I1:I2) = REPLTXT

      ENDDO
         
      ! Return to calling program
      END SUBROUTINE StrRepl

!------------------------------------------------------------------------------


      SUBROUTINE TxtExt(ch,text,col,word,iflg)
!
!     PURPOSE: TxtExt extracts a sequence of characters from
!              text and transfers them to word.  The extraction
!              procedure uses a set of character "delimiters"
!              to denote the desired sequence of characters.
!              For example if ch=' ', the first character sequence
!              bracketed by blank spaces will be returned in word.
!              The extraction procedure begins in column, col,
!              of TEXT.  If text(col:col) = ch (any character in
!              the character string), the text is returned beginning
!              with col+1 in text (i.e., the first match with ch
!              is ignored).
!
!              After completing the extraction, col is incremented to
!              the location of the first character following the
!              end of the extracted text.
!
!              A status flag is also returned with the following
!              meaning(s)
!
!                 IF iflg = -1, found a text block, but no more characters
!                               are available in TEXT
!                    iflg = 0, task completed sucessfully (normal term)
!                    iflg = 1, ran out of text before finding a block of
!                              text.
!
!       COMMENTS: TxtExt is short for Text Extraction.  This routine
!                 provides a set of powerful line-by-line
!                 text search and extraction capabilities in
!                 standard FORTRAN.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        CntMat                    CHARPAK.FOR
!        TxtExt                    CHARPAK.FOR
!        FillStr                   CHARPAK.FOR
!        CopyTxt                   CHARPAK.FOR
!
!        other routines are indirectly called.
!      AUTHOR: Robert D. Stewart
!        DATE: Jan. 1st, 1995
!
!      REVISIONS: FEB 22, 1996.  Slight bug fix (introduced by a
!        (recent = FLIB 1.04) change in the CntMat routine)
!        so that TxtExt correctlyhandles groups of characters
!        delimited by blanks).
!
!      MODIFICATIONS by Bob Yantosca (6/25/02)
!        (1) Replace call to FILLSTR with F90 intrinsic REPEAT
!
      CHARACTER*(*) ch,text,word
      INTEGER col,iflg
      INTEGER Tmax,T1,T2,imat
      LOGICAL again,prev

!     Length of text
      Tmax = LEN(text)

!     Fill Word with blanks
      WORD = REPEAT( ' ', LEN( WORD ) )
      
      IF (col.GT.Tmax) THEN
!       Text does not contain any characters past Tmax.
!       Reset col to one and return flag = {error condition}
        iflg = 1
        col = 1
      ELSEIF (col.EQ.Tmax) THEN
!       End of TEXT reached
        CALL CntMat(ch,text(Tmax:Tmax),imat)
        IF (imat.EQ.0) THEN
!         Copy character into Word and set col=1
          CALL CopyTxt(1,Text(Tmax:Tmax),Word)
          col = 1
          iflg = -1
        ELSE
!         Same error condition as if col.GT.Tmax
          iflg = 1
        ENDIF
      ELSE
!       Make sure column is not less than 1
        IF (col.LT.1) col=1
        CALL CntMat(ch,text(col:col),imat)
        IF (imat.GT.0) THEN
          prev=.true.
        ELSE
          prev=.false.
        ENDIF
        T1=col
        T2 = T1

        again = .true.
        DO WHILE (again)
!         Check for a match with a character in ch
          CALL CntMat(ch,text(T2:T2),imat)
          IF (imat.GT.0) THEN
!           Current character in TEXT matches one (or more) of the
!           characters in ch.
            IF (prev) THEN
              IF (T2.LT.Tmax) THEN
!               Keep searching for a block of text
                T2=T2+1
                T1=T2
              ELSE
!               Did not find any text blocks before running
!               out of characters in TEXT.
                again=.false.
                iflg=1
              ENDIF
            ELSE
!             Previous character did not match ch, so terminate.
!             NOTE: This is "NORMAL" termination of the loop
              again=.false.
              T2=T2-1
              iflg = 0
            ENDIF
          ELSEIF (T2.LT.Tmax) THEN
!           Add a letter to the current block of text
            prev = .false.
            T2=T2+1
          ELSE
!           Reached the end of the characters in TEXT before reaching
!           another delimiting character.  A text block was identified
!           however.
            again=.false.
            iflg=-1
          ENDIF
        ENDDO

        IF (iflg.EQ.0) THEN
!         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = T2+1
        ELSE
!         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = 1
        ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE TxtExt

!------------------------------------------------------------------------------

      SUBROUTINE CntMat(str1,str2,imat)
!
!     Count the number of characters in str1 that match
!     a character in str2.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!          LENTRIM                CharPak
!
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!     COMMENTS: Revised slightly (2-5-1996) so that trailing
!               blanks in str1 are ignored.  Revised again
!               on 3-6-1996.
!
      CHARACTER*(*) str1,str2
      INTEGER imat
      INTEGER L1,L2,i,j
      LOGICAL again

      L1 = MAX(1,LEN_TRIM(str1))
      L2 = LEN(str2)
      imat = 0
      DO i=1,L1
        again = .true.
        j = 1
        DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
            imat = imat+1
            again = .false.
          ELSEIF (j.LT.L2) THEN
            j=j+1
          ELSE
            again = .false.
          ENDIF
        ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CntMat

!------------------------------------------------------------------------------

      SUBROUTINE CopyTxt(col,str1,str2)
!
!     PURPOSE: Write all of the characters in str1 into variable
!              str2 beginning at column, col.  If the length of str1
!              + col is longer than the number of characters str2
!              can store, some characters will not be transfered to
!              str2.  Any characters already existing in str2 will
!              will be overwritten.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   DEC. 24, 1993
!     AUTHOR: R.D. STEWART
!
      CHARACTER*(*) str2,str1
      INTEGER col,ilt1,i1,i,j,ic

      i1 = LEN(str2)
      IF (i1.GT.0) THEN
        ilt1 = LEN(str1)
        IF (ilt1.GT.0) THEN
          ic = MAX0(col,1)
          i = 1
          j = ic
          DO WHILE ((i.LE.ilt1).and.(j.LE.i1))
            str2(j:j) = str1(i:i)
            i = i + 1
            j = ic + (i-1)
          ENDDO
        ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE CopyTxt

!------------------------------------------------------------------------------

      SUBROUTINE TRANUC(text)
!
!     PURPOSE: Tranlate a character variable to all upper case letters.
!              Non-alphabetic characters are not affected.
!
!    COMMENTS: The original "text" is destroyed.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
      CHARACTER*(*) text
      INTEGER iasc,i,ilen

      ilen = LEN(text)
      DO i=1,ilen
        iasc = ICHAR(text(i:i))
        IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
          text(i:i) = CHAR(iasc-32)
        ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE TRANUC

END MODULE string_utils_module