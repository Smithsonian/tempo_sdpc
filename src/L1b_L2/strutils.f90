MODULE strutils

CONTAINS

  !---------------------------------------------------------------------------
  SUBROUTINE find_endstring ( slen, cstring, istart, iend)
    IMPLICIT NONE

    INTEGER,              INTENT (IN) :: slen, istart
    CHARACTER (LEN=slen), INTENT (IN) :: cstring

    INTEGER, INTENT (OUT) :: iend

    iend = INDEX ( cstring(istart:slen), "," )
    RETURN
  END SUBROUTINE find_endstring

  !---------------------------------------------------------------------------
  SUBROUTINE find_dimension ( slen, cstring, istart, iend, dimstring )
    IMPLICIT NONE

    INTEGER,              INTENT (IN ) :: slen, istart, iend
    CHARACTER (LEN=slen), INTENT (IN) :: cstring

    CHARACTER (LEN=*),    INTENT (OUT) :: dimstring

    dimstring = cstring(istart:iend)

    RETURN
  END SUBROUTINE find_dimension

  !---------------------------------------------------------------------------
  SUBROUTINE slice_string ( slen, instring, istart, ostring )
    IMPLICIT NONE

    INTEGER,              INTENT (IN) :: slen, istart
    CHARACTER (LEN=slen), INTENT (IN) :: instring

    CHARACTER (LEN=*),    INTENT (OUT) :: ostring

    ostring = TRIM(ADJUSTL( instring(istart:slen) ))

    RETURN
  END SUBROUTINE slice_string

  !---------------------------------------------------------------------------

  SUBROUTINE remove_quotes ( quotestring )

    ! ==========================================================
    ! Subroutine to remove any quotation marks (") from a string
    ! ==========================================================

    USE OMSAO_precision_module, ONLY: i4
    IMPLICIT NONE

    ! ---------------------
    ! Input/modified string
    ! ---------------------
    CHARACTER (LEN=*), INTENT (INOUT) :: quotestring

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: i, j
    CHARACTER (LEN=LEN(quotestring)) :: tmpstring

    IF ( INDEX(quotestring, '"') == 0 .OR. LEN(quotestring) == 0 ) RETURN

    tmpstring = ""; j = 1
    DO i = 1, LEN(quotestring)
      SELECT CASE ( quotestring(i:i) )
      CASE ( '"' )
        ! skip; nothing else to do
      CASE DEFAULT
        tmpstring(j:j) = quotestring(i:i)
        j = j+1
      END SELECT
    END DO

    quotestring = tmpstring

    RETURN
  END SUBROUTINE remove_quotes

  !---------------------------------------------------------------------------

  SUBROUTINE get_substring ( string, sstart, substring, nsubstring, eostring )

    USE OMSAO_precision_module, ONLY: i4
    IMPLICIT NONE

    ! ==========================================================================
    ! Given string STRING, extract first space- or comma-delimited substring
    ! beginning on or after position SSTART, and return its value, SUBSTRING and
    ! its length, NSUBSTRING; update STRING and SSTART to remove this substring.
    !
    ! F90 version of the original GET_TOKEN by J. Lavanigno
    ! ==========================================================================

    ! NSUBSTRING represents the token's length without trailing blanks. It's
    ! set to 0 if no substring was found in STRING: this can be used to determine
    ! when to stop looking.
    !
    ! This routine makes no attempt to detect or handle the case in which
    ! the token is bigger than the token buffer.  It's assumed that the
    ! caller will declare line and token to be the same size so that this
    ! won't ever happen.

    CHARACTER (LEN = *), INTENT (INOUT) :: string
    INTEGER   (KIND=i4), INTENT (INOUT) :: sstart

    CHARACTER (LEN =LEN(string)), INTENT (OUT) :: substring
    INTEGER   (KIND=i4),          INTENT (OUT) :: nsubstring
    LOGICAL,                      INTENT (OUT) :: eostring

    ! local variables
    CHARACTER (LEN=1)   :: char
    INTEGER   (KIND=i4) :: lstart, lend, nline

    ! ------------------
    ! Initialize outputs
    ! ------------------
    substring  = ' '  ;  nsubstring = 0  ;  eostring = .FALSE.

    nline = LEN(TRIM(ADJUSTL(string)) )

    ! ------------------------------------------------------
    ! We are working on character variables, i.e., positions
    ! are always larger than 0
    ! ------------------------------------------------------
    IF ( sstart <= 0 ) sstart = 1

    ! ----------------------------------------------
    ! Check first whether we have to any work at all
    ! ----------------------------------------------
    IF ( sstart >= nline ) THEN
      eostring = .TRUE.  ;  RETURN
    END IF

    ! --------------------------------------------------------
    ! Find first character in line that's not a blank or comma
    ! --------------------------------------------------------
    findchar: DO lstart = sstart, nline
      char = string(lstart:lstart)
      SELECT CASE ( char )
      CASE ( ' ' )
        IF ( lstart == nline ) THEN
          sstart = nline + 1; RETURN  ! there are no further substrings in STRING
        END IF
      CASE ( ',' )
        IF ( lstart == nline ) THEN
          sstart = nline + 1; RETURN  ! there are no further substrings in STRING
        END IF
      CASE DEFAULT
        EXIT findchar
      END SELECT
    END DO findchar

    ! --------------------------------------------------
    ! Start of the next substring is at position LSTART.
    ! Now find separator that ends the substring.
    ! --------------------------------------------------
    findsep: DO lend = lstart + 1, nline
      char = string (lend:lend)
      If ( char == ' ' .OR. char == ',' ) EXIT findsep
    END DO findsep
    IF ( (lend == nline) .AND. (char /= ' ' .AND. char /= ',') ) lend = lend + 1
    lend = lend - 1

    ! --------------------
    ! Output the substring
    ! --------------------
    nsubstring = lend - lstart + 1

    substring(1:nsubstring) = string (lstart:lend)

    ! -----------------------------------------------------------------
    ! The next substring, if any, starts at least two characters beyond
    ! the end of thelast (we have to skip over the comma or space that
    ! marks the substring's end).
    ! -----------------------------------------------------------------
    sstart = lend + 2

    ! ---------------------------------------------
    ! Final check, whether we have to any more work
    ! ---------------------------------------------
    IF ( sstart >= nline ) eostring = .TRUE.

    RETURN
  END SUBROUTINE get_substring

END MODULE
