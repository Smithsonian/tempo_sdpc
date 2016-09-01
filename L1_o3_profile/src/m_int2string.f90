!>Convert integer to string. 
module m_int2string

  public int2string
  private

contains

  !Had to be moved to separate module to avoid circular dependency in build

  FUNCTION int2string ( int_num, ni ) RESULT ( int_str )

    ! ===============================================================
    ! Converts INTEGER number INT_NUM to STRING INT_STR of length NI
    ! or the number of digits in INT_NUM, whatever is larger.
    ! ===============================================================

    USE OMSAO_precision_module, ONLY: i4
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),  INTENT (IN)  :: int_num, ni

    ! ---------------
    ! Result variable
    ! ---------------
    CHARACTER (LEN=4) :: int_str

    ! ------------------------------
    ! Local variables and parameters
    ! ------------------------------
    ! * Arrays containing indices for number strings in ASCII table
    INTEGER (KIND=i4),                   PARAMETER :: n = 10
    INTEGER (KIND=i4), DIMENSION(0:n-1)            :: aidx
    CHARACTER (LEN=1), DIMENSION(0:n-1), PARAMETER :: astr = &
         (/ "0","1","2","3","4","5","6","7","8","9" /)
    ! * Temporary and loop variables
    INTEGER (KIND=i4)                              :: i, k, nd, tmpint, ld

    ! ----------------------------------------------------
    ! Compute the index entries of ASTR in the ASCII table
    ! ----------------------------------------------------
    aidx = IACHAR(astr)

    ! ----------------------------------------------------------
    ! Find the number of digits in INT_NUM. This is equal to the 
    ! truncated integer of LOG10(INT_NUM) plus 1.
    ! ----------------------------------------------------------
    SELECT CASE ( int_num )
    CASE ( 0:9 ) 
      nd = 1
    CASE ( 10: )
      nd = INT ( LOG10(REAL(int_num)) ) + 1
    CASE DEFAULT
      int_str = '?' 
      RETURN
    END SELECT

    ! -------------------------------------------------------------
    ! We may want to create a string that is longer than the number
    ! of digits in INT_NUM. This will create leading "0"s.
    ! -------------------------------------------------------------
    nd = MAX ( nd, ni )

    ! ----------------------------------------------
    ! Convert the integer to string "digit by digit"
    ! ----------------------------------------------
    int_str = "" 
    tmpint = int_num
    DO i = 1, nd
      ld = MOD ( tmpint, 10 )        ! Current last digit
      tmpint = ( tmpint - ld ) / 10  ! Remove current last digit from INT_STR
      k = nd - i + 1                 ! Position of current digit in INT_STR
      int_str(k:k) = ACHAR(aidx(ld)) ! Convert INTEGER digit to CHAR
    END DO

    RETURN
  END FUNCTION int2string



end module m_int2string
