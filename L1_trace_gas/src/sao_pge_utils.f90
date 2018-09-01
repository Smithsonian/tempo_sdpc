MODULE sao_pge_utils

  use tell_module

  implicit none

  public day_of_year, utc_julian_date_and_time, skip_to_filemark, signdp, &
       interpolation, print_array, get_pge_ident
  private year_month_day, check_read_status, fill_nonoverlap

CONTAINS
  ! ===========================================================================
  !
  ! Collection of subroutines that relate to PGE identification and processing 
  ! aids.
  !
  ! ===========================================================================

  FUNCTION day_of_year ( year, month, day ) RESULT ( jday )

    USE OMSAO_precision_module, ONLY: i4
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: year, month, day

    ! ---------------
    ! Result variable
    ! ---------------
    INTEGER (KIND=i4) :: jday

    ! ------------------------------
    ! Local variables and parameters
    ! ------------------------------
    INTEGER (KIND=i4), PARAMETER :: n_month = 12
    INTEGER (KIND=i4), DIMENSION (n_month), PARAMETER :: &
         days_per_month = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    ! -------------------------------------------------------
    ! First find day of the year in a regular (non leap) year
    ! -------------------------------------------------------
    SELECT CASE ( month )
    CASE ( 1 )
      jday = day
    CASE (2:12)
      jday = SUM ( days_per_month(1:month-1) ) + day
    CASE DEFAULT
      jday = -9999
    END SELECT

    ! -------------------------
    ! Now apply leap year rules
    ! ------------------------------------
    ! * Divisible by 4:   leap year
    ! * Divisible by 100: not a leap year
    ! * Divisible by 400: leap year
    ! ------------------------------------
    IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. MOD(year,400) == 0 ) ) jday = jday+1

    RETURN
  END FUNCTION day_of_year

  SUBROUTINE get_pge_ident (in_name, out_idx, out_name, errstat)

    ! ====================================================
    ! Find name and index of current PGE from input string
    ! ====================================================

    USE OMSAO_precision_module, ONLY: i4
    USE OMSAO_indices_module, ONLY: sao_pge_names, sao_pge_min_idx, sao_pge_max_idx
    !USE OMSAO_errstat_module

    IMPLICIT NONE

    ! --------------
    ! Input variable
    ! --------------
    CHARACTER (LEN=*), INTENT (IN) :: in_name

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT(OUT) :: out_idx
    CHARACTER (LEN=*), INTENT(OUT) :: out_name
    INTEGER   (KIND=i4), INTENT (INOUT) :: errstat

    ! --------------
    ! Local variable
    ! --------------
    INTEGER (KIND=i4) :: i !, locerrstat

    if (errstat /= 0) return

    ! ---------------------------
    ! Initialize variables
    ! ---------------------------
    out_idx = -1
    out_name(1:1) = '?'
    !locerrstat = pge_errstat_ok

    ! -------------------------------------------------
    ! Find name and index by looping over all SAO PGEs.
    ! Not very elegant but simple and effective.
    ! -------------------------------------------------
    getpge: DO i = sao_pge_min_idx, sao_pge_max_idx
      IF ( TRIM(ADJUSTL(in_name)) == TRIM(ADJUSTL(sao_pge_names(i))) ) THEN
        out_name = sao_pge_names(i)
        out_idx = i
        EXIT getpge
      END IF
    END DO getpge

    IF ( (out_idx == -1) .OR. (out_name(1:1) == '?' )) then
      call tell_set_error (tell_runtime_error) !locerrstat = pge_errstat_error
    endif
    !errstat = MAX ( errstat, locerrstat )
    RETURN

  END SUBROUTINE get_pge_ident

  !unused SUBROUTINE pge_error_message ( errstat, ok_msg, warn_msg, err_msg )
  !unused 
  !unused   ! ====================================================================
  !unused   ! Check PGE error status and report appropriate message. STOP on error
  !unused   ! ====================================================================
  !unused   USE OMSAO_precision_module, ONLY: i4
  !unused   USE OMSAO_errstat_module
  !unused 
  !unused   IMPLICIT NONE
  !unused 
  !unused   ! ---------------
  !unused   ! Input variables
  !unused   ! ---------------
  !unused   INTEGER   (KIND=i4), INTENT (IN) :: errstat
  !unused   CHARACTER (LEN=*),   INTENT (IN) :: ok_msg, warn_msg, err_msg
  !unused 
  !unused   SELECT CASE ( errstat )
  !unused   CASE ( pge_errstat_ok )
  !unused     WRITE (*, '(A)') TRIM(ADJUSTL(ok_msg))
  !unused   CASE ( pge_errstat_warning )
  !unused     WRITE (*, '(A,A)') 'WARNING: ', TRIM(ADJUSTL(warn_msg))
  !unused   CASE ( pge_errstat_error )
  !unused     WRITE (*, '(A,A)') 'ERROR: ', TRIM(ADJUSTL(err_msg))
  !unused     STOP 1
  !unused   END SELECT
  !unused 
  !unused   RETURN
  !unused END SUBROUTINE pge_error_message

  !JCH - A character function must not be declared to return a value with LEN=*.
  !JCH - On the other hand, this function is equivalent to
  !JCH -  write (int_str, '(i0.${ni})') int_num
  !JCH - where ni is the field width (necessary only if ni>1),
  !JCH - so I don't think we need this function at all.
  !JCH - For that reason, I'm commenting it out.
  !JCH FUNCTION int2string ( int_num, ni ) RESULT ( int_str )
  !JCH 
  !JCH   ! ===============================================================
  !JCH   ! Converts INTEGER number INT_NUM to STRING INT_STR of length NI
  !JCH   ! or the number of digits in INT_NUM, whatever is larger.
  !JCH   ! ===============================================================
  !JCH 
  !JCH   USE OMSAO_precision_module, ONLY: i4
  !JCH   IMPLICIT NONE
  !JCH 
  !JCH   ! ---------------
  !JCH   ! Input variables
  !JCH   ! ---------------
  !JCH   INTEGER (KIND=i4),  INTENT (IN)  :: int_num, ni
  !JCH 
  !JCH   ! ---------------
  !JCH   ! Result variable
  !JCH   ! ---------------
  !JCH   CHARACTER (LEN=*) :: int_str
  !JCH 
  !JCH   ! ------------------------------
  !JCH   ! Local variables and parameters
  !JCH   ! ------------------------------
  !JCH   ! * Arrays containing indices for number strings in ASCII table
  !JCH   INTEGER (KIND=i4),                   PARAMETER :: n = 10
  !JCH   INTEGER (KIND=i4), DIMENSION(0:n-1)            :: aidx
  !JCH   CHARACTER (LEN=1), DIMENSION(0:n-1), PARAMETER :: astr = (/ "0","1","2","3","4","5","6","7","8","9" /)
  !JCH   ! * Temporary and loop variables
  !JCH   INTEGER (KIND=i4)                              :: i, k, nd, tmpint, ld
  !JCH 
  !JCH   ! ----------------------------------------------------
  !JCH   ! Compute the index entries of ASTR in the ASCII table
  !JCH   ! ----------------------------------------------------
  !JCH   aidx = IACHAR(astr)
  !JCH 
  !JCH   ! ----------------------------------------------------------
  !JCH   ! Find the number of digits in INT_NUM. This is equal to the
  !JCH   ! truncated integer of LOG10(INT_NUM) plus 1.
  !JCH   ! ----------------------------------------------------------
  !JCH   SELECT CASE ( int_num )
  !JCH   CASE ( 0:9 )
  !JCH     nd = 1
  !JCH   CASE ( 10: )
  !JCH     nd = INT ( LOG10(REAL(int_num)) ) + 1
  !JCH   CASE DEFAULT
  !JCH     int_str = '?' ; RETURN
  !JCH   END SELECT
  !JCH 
  !JCH   ! -------------------------------------------------------------
  !JCH   ! We may want to create a string that is longer than the number
  !JCH   ! of digits in INT_NUM. This will create leading "0"s.
  !JCH   ! -------------------------------------------------------------
  !JCH   nd = MAX ( nd, ni )
  !JCH 
  !JCH   ! ----------------------------------------------
  !JCH   ! Convert the integer to string "digit by digit"
  !JCH   ! ----------------------------------------------
  !JCH   int_str = "" ; tmpint = int_num
  !JCH   DO i = 1, nd
  !JCH     ld = MOD ( tmpint, 10 )        ! Current last digit
  !JCH     tmpint = ( tmpint - ld ) / 10  ! Remove current last digit from INT_STR
  !JCH     k = nd - i + 1                 ! Position of current digit in INT_STR
  !JCH     int_str(k:k) = ACHAR(aidx(ld)) ! Convert INTEGER digit to CHAR
  !JCH   END DO
  !JCH 
  !JCH   RETURN
  !JCH END FUNCTION int2string

  SUBROUTINE utc_julian_date_and_time ( year, month, day, julday, hour, minute, second )

    ! -------------------------------------------------------------
    ! Converts current date and time to UTC time and "Julian" date,
    ! i.e., day of the year.
    ! -------------------------------------------------------------

    USE OMSAO_precision_module, ONLY: i4
    !USE sao_pge_utils, ONLY: day_of_year
    IMPLICIT NONE

    ! ---------------
    ! RESULT
    ! ---------------
    INTEGER (KIND=i4), INTENT (OUT) :: year, month, day, julday, hour, minute, second

    !INTEGER (KIND=i4), DIMENSION(8), INTENT (OUT) :: utc_juldate

    ! ---------------
    ! Local variables
    ! ---------------
    ! * Arguments for DATE_AND_TIME, and position indices for Year, Month, Day, Zone, Hour, Minute, Seconds
    CHARACTER (LEN= 8)                 :: date
    CHARACTER (LEN=10)                 :: time
    CHARACTER (LEN= 5)                 :: zone
    INTEGER   (KIND=i4), DIMENSION (8) :: date_vector
    INTEGER   (KIND=i4), PARAMETER     :: y_idx=1, m_idx=2, d_idx=3, z_idx=4, hh_idx=5, mm_idx=6, ss_idx=7
    ! * Some MAX values for Minutes and Hours (not that we would expect those to change :-)
    INTEGER   (KIND=i4), PARAMETER     :: max_mi = 60, max_hr = 24, max_dy = 31, max_mo = 12
    ! * Other local variables
    INTEGER   (KIND=i4)                :: del_mm, del_hh, max_julday, yyy, mmm, ddd

    ! ------------------
    ! External functions
    ! ------------------
    !INTEGER (KIND=i4), EXTERNAL :: day_of_year

    ! ----------------------------------------------------------
    ! First copy the input variable to the output variable. That
    ! saves work on premature return.
    ! ----------------------------------------------------------
    CALL DATE_AND_TIME ( date, time, zone, date_vector )

    ! -----------------------------------------------
    ! Initialize all output variables (except JULDAY)
    ! -----------------------------------------------
    year   = date_vector(y_idx )
    month  = date_vector(m_idx )
    day    = date_vector(d_idx )
    hour   = date_vector(hh_idx)
    minute = date_vector(mm_idx)
    second = date_vector(ss_idx)

    ! ----------------------------------------------------------------------
    ! Next, find the day of the year. At this point, we also compute the
    ! maximum number of days per year. This will come in handy further down.
    ! ----------------------------------------------------------------------
    julday     = day_of_year ( year, month,  day    )
    max_julday = day_of_year ( year, max_mo, max_dy )

    ! -----------------------------------------------
    ! Check if there is any difference to UTC at all.
    ! -----------------------------------------------
    IF ( ABS(date_vector(z_idx)) == 0 ) RETURN

    ! ------------------------------------------------------------
    ! Apply difference in UTC time to local time vector. This is
    ! given in minutes, thus apply to minutes index.
    !
    ! REMEMBER: This is the DIFFERENCE to UTC, so we have to
    !           SUBTRACT whatever difference from the current time.
    ! -------------------------------------------------------------
    del_mm = MOD ( date_vector(z_idx), max_mi )
    minute = minute - del_mm

    ! ----------------------------------------
    ! Compute hour difference, return if none.
    ! ----------------------------------------
    del_hh =   ( date_vector(z_idx) - del_mm ) / max_mi
    IF ( del_hh == 0 ) RETURN

    ! --------------------------------------------------------
    ! Apply hour difference and check for/apply day difference
    ! --------------------------------------------------------
    hour = hour - del_hh

    ! -------------------------------------------------------------
    ! If we reach here, we may have had a day difference. We could
    ! make our lives miserable by checking for all the months of
    ! the year. But what we ultimately want is the day of the year
    ! (a proxy for Julian day), and that makes life a lot easier.
    ! -------------------------------------------------------------
    SELECT CASE ( hour )
    CASE ( 0:23 )
      RETURN
    CASE ( :-1 )
      hour = max_hr + hour
      day  = day - 1
    CASE ( 24: )
      hour = hour - max_hr
      day  = day  + 1
    END SELECT

    ! --------------------------------------------------------------
    ! Finally, if we ever reach here, we've had a year difference.
    ! Apply correction, compute the new day of the year, and return.
    ! --------------------------------------------------------------
    yyy = year ; mmm = month ; ddd = day                     ! save old values
    CALL year_month_day ( yyy, mmm, ddd, year, month, day )  ! overwrite with new ones
    julday = day_of_year ( year, month, day )                ! recompute day-of-year

    RETURN
  END SUBROUTINE utc_julian_date_and_time

  SUBROUTINE year_month_day ( year, month, day, newyear, newmonth, newday )

    USE OMSAO_precision_module, ONLY: i4
    IMPLICIT NONE

    ! ----------------------
    ! Input/Output variables
    ! ----------------------
    INTEGER (KIND=i4), INTENT (IN)    :: year, month, day
    INTEGER (KIND=i4), INTENT (INOUT) :: newyear, newmonth, newday

    ! ------------------------------
    ! Local variables and parameters
    ! ------------------------------
    INTEGER (KIND=i4), PARAMETER :: n_month = 12
    INTEGER (KIND=i4), DIMENSION (n_month), PARAMETER :: &
         days_per_month    = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
         ly_days_per_month = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    INTEGER (KIND=i4), DIMENSION (n_month) :: dpm

    ! -------------------------
    ! Check for leap year
    ! ------------------------------------
    ! * Divisible by 4:   leap year
    ! * Divisible by 100: not a leap year
    ! * Divisible by 400: leap year
    ! ------------------------------------
    IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. MOD(year,400) == 0 ) ) THEN
      dpm = days_per_month
    ELSE
      dpm = ly_days_per_month
    END IF

    newyear = year  ;  newmonth = month  ;  newday = day

    ! -------------------------------------------------------
    ! If DAY is within bounds there's nothing to do
    ! (MONTH hasn't been adjusted yet in the calling program)
    ! -------------------------------------------------------
    IF ( day > 0 .AND. day <= dpm(month) ) RETURN

    ! ---------------------------------------------------------------------------
    ! If we reach here, DAY must be either too large or <= 0. The crucial
    ! variable now is MONTH, since it might lead to an additional change in YEAR.
    ! ---------------------------------------------------------------------------
    SELECT CASE ( month )
    CASE ( 2:11 )   ! The easiest case: change only MONTH and DAY
      IF ( day <= 0 ) THEN
        newmonth = month - 1
        newday   = dpm(newmonth) + day  ! DAY is negative
      ELSE
        newmonth = month + 1
        newday   = day - dpm(month)
      END IF
    CASE ( 1 )      ! January requires a YEAR change if DAY <= 0
      IF ( day <= 0 ) THEN
        newmonth = 12
        newday   = dpm(newmonth) + day  ! DAY is negative
        newyear  = year - 1
      ELSE
        newmonth = month + 1
        newday   = day - dpm(month)
      END IF
    CASE ( 12 )      ! December requires a YEAR change if DAY > 31
      IF ( day <= 0 ) THEN
        newmonth = month - 1
        newday   = dpm(newmonth) + day  ! DAY is negative
      ELSE
        newmonth = 1
        newday   = day - dpm(month)
        newyear  = year + 1
      END IF
    END SELECT

    RETURN
  END SUBROUTINE year_month_day

  !FUNCTION upper_case ( mixstring ) RESULT ( upcase )
  !
  !  ! =============================================
  !  ! Function to convert strings to all upper case
  !  ! =============================================
  !
  !  USE OMSAO_precision_module, ONLY: i4
  !  IMPLICIT NONE
  !
  !  ! ------------------------------------------------------------
  !  CHARACTER (LEN=*), INTENT (IN) :: mixstring   ! Input string
  !  ! ------------------------------------------------------------
  !
  !  ! -------------------------------------------------------------
  !  CHARACTER (LEN=len(mixstring)) :: upcase      ! Result string
  !  ! -------------------------------------------------------------
  !
  !  ! -------------------------------------------------------------
  !  INTEGER (KIND=i4) :: i                         ! Local variables
  !  ! -------------------------------------------------------------
  !
  !  DO i = 1, LEN(mixstring)
  !    SELECT CASE ( ICHAR(mixstring(i:i)) )
  !    CASE ( 97:122 )
  !      upcase(i:i) = ACHAR(ICHAR(mixstring(i:i))-32)
  !    CASE DEFAULT
  !      upcase(i:i) = mixstring(i:i)
  !    END SELECT
  !  END DO
  !
  !  RETURN
  !END FUNCTION upper_case

  SUBROUTINE check_read_status ( ios, file_read_stat )

    USE OMSAO_precision_module, ONLY: i4
    USE OMSAO_errstat_module,   ONLY: file_read_ok, file_read_failed, file_read_eof
    IMPLICIT NONE

    INTEGER (KIND=i4), INTENT (IN)  :: ios
    INTEGER (KIND=i4), INTENT (OUT) :: file_read_stat

    SELECT CASE ( ios )
    CASE ( :-1 )  ;  file_read_stat = file_read_eof
    CASE (   0 )  ;  file_read_stat = file_read_ok
    CASE DEFAULT  ;  file_read_stat = file_read_failed
    END SELECT

    RETURN
  END SUBROUTINE check_read_status

  SUBROUTINE skip_to_filemark ( funit, lm_string, lastline, file_read_stat )

    USE OMSAO_precision_module, ONLY: i4
    USE OMSAO_parameters_module, ONLY: MAX_STR_LEN
    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER   (KIND=i4), INTENT (IN) :: funit
    CHARACTER (LEN=*),   INTENT (IN) :: lm_string

    ! ================
    ! Output variables
    ! ================
    INTEGER   (KIND=i4), INTENT (OUT) :: file_read_stat
    CHARACTER (LEN=*),   INTENT (OUT) :: lastline

    ! ===============
    ! Local variables
    ! ===============
    INTEGER   (KIND=i4)      :: lmlen, ios
    CHARACTER (LEN=MAX_STR_LEN) :: tmpline

    ! -------------------------------------------
    ! Determine the length of the string landmark
    ! -------------------------------------------
    lmlen = LEN(TRIM(ADJUSTL(lm_string)))

    ! -------------------------------------------------------
    ! Read lines in the file until we either find the string,
    ! reach the end of the file, or reading fails otherwise.
    ! ----------------------------------------------------
    ios = 0
    getlm: DO WHILE ( ios == 0 )
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) tmpline
      tmpline = TRIM(ADJUSTL(tmpline))
      IF ( ios /= 0 .OR. tmpline(1:lmlen) == lm_string ) EXIT getlm
    END DO getlm

    ! ---------------------------------------------------
    ! Return the last line read for the case that we need
    ! to extract further information from it
    ! ---------------------------------------------------
    lastline = TRIM(ADJUSTL(tmpline))

    CALL check_read_status ( ios, file_read_stat )

    RETURN
  END SUBROUTINE skip_to_filemark


  REAL (KIND=KIND(1.0D0)) FUNCTION signdp ( x )

    USE OMSAO_precision_module, ONLY: r8
    IMPLICIT NONE

    REAL (KIND=r8), INTENT (IN) :: x

    signdp = 0.0_r8
    IF ( x < 0.0_r8 ) THEN
      signdp = -1.0_r8
    ELSE
      signdp = +1.0_r8
    END IF

    RETURN
  END FUNCTION signdp

  SUBROUTINE interpolation ( &
       n1, x1, y1, n2, x2, y2, filltype, fillval, did_full_range, err )

    USE OMSAO_precision_module
    !USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_error
    USE ezspline_interpolation, ONLY: ezspline_1d_interpolation
    use arrayutils, only: array_locate_r8
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*),                 INTENT (IN) :: filltype
    INTEGER (KIND=i4),                 INTENT (IN) :: n1, n2
    REAL    (KIND=r8),                 INTENT (IN) :: fillval
    REAL    (KIND=r8), DIMENSION (n1), INTENT (IN) :: x1, y1
    REAL    (KIND=r8), DIMENSION (n2), INTENT (IN) :: x2

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                 INTENT (INOUT) :: err
    REAL    (KIND=r8), DIMENSION (n2), INTENT (OUT)   :: y2
    LOGICAL,                           INTENT (OUT)   :: did_full_range

    ! --------------
    ! Local variable
    ! --------------
    INTEGER (KIND=i4)                 :: imin, imax, nloc !, locerrstat
    REAL    (KIND=r8), DIMENSION (n2) :: xtmp, ytmp

    if (err /= 0) return
    !locerrstat = 0 ! pge_errstat_ok

    ! -------------------------------
    ! Initialize interpolation output
    ! -------------------------------
    y2(1:n2) = 0.0_r8

    ! --------------------------------------------------------------------------
    ! Find indices in radiance wavelength spectrum that cover reference spectrum
    ! --------------------------------------------------------------------------
    imin = -1 ; imax = -1
    !imin = MAXVAL ( MINLOC ( x2(1:n2), MASK=(x2(1:n2) >= x1( 1)) ) )
    !imax = MAXVAL ( MAXLOC ( x2(1:n2), MASK=(x2(1:n2) <= x1(n1)) ) )
    CALL array_locate_r8 ( n2, x2(1:n2), x1( 1), 'GE', imin )
    CALL array_locate_r8 ( n2, x2(1:n2), x1(n1), 'LE', imax )

    ! -------------------------------------------------------------------------------
    ! Check whether we have the whole wavelength range. We don't set the error status
    ! variable for this case, because it is too insignificant to set non-Zero exit at
    ! PGE termination. Instead, we check in the calling subroutine for did_full_range
    ! and report this at higher verbosity thresholds.
    ! ------------------------------------------------------------------------------
    did_full_range = .TRUE.
    IF ( imin /= 1 .OR. imax /= n2 ) did_full_range = .FALSE.

    ! --------------------------------------------------------------------------
    ! Now that we know the first and last index to cover with the interpolation,
    ! we can treat all cases alike. Only we need make sure that the number of
    ! interpolation points is consistent. And if we don't have *any* points,
    ! then we must return without calling the interpolation routine.
    ! --------------------------------------------------------------------------
    nloc = imax - imin + 1

    ! ------------------------------------------------------------------
    ! Anything less than 4 data points will make the interpolation fail.
    ! So before we choose to interpolate or not, we set the interpolated
    ! array to Zero. That is easier than a IF_THEN_ELSEing based on the
    ! IMIN and IMAX indices.
    ! ------------------------------------------------------------------
    SELECT CASE ( nloc )
    CASE ( :3 ) ! Less than 4 data points available for interpolation
      !errstat = -1 ! pge_errstat_error
      call tell_log (0, "interpolation: less than 4 data points available for interpolation")
      err = -1
      RETURN
    CASE DEFAULT
      xtmp(1:nloc) = x2(imin:imax)
      CALL ezspline_1d_interpolation (                &
           n1,   x1  (1:n1),   y1  (1:n1),            &
           nloc, xtmp(1:nloc), ytmp(1:nloc), err) ! locerrstat )
      if (err /= 0) return
      y2(imin:imax) = ytmp(1:nloc)
      CALL fill_nonoverlap ( n2, y2(1:n2), imin, imax, filltype, fillval )
      ! -----------------------------------------------------------------
      ! If we have non-zero exit status, something must have gone wrong
      ! in the interpolation. Set PGE_ERROR_STATUS to ERROR in this case.
      ! -----------------------------------------------------------------
      !errstat = MAX ( errstat, locerrstat )
    END SELECT

    RETURN
  END SUBROUTINE interpolation


!  SUBROUTINE find_overlap ( n1, arr1, n2, arr2, i1, i2 )
!
!    ! ----------------------------------------------------------------
!    !
!    ! Given arrays ARR1 of dimension N1 and ARR2 of dimension N2, find
!    ! array positions I1 and I2 in ARR2 such that ARR2(I1:i2) fully
!    ! overlaps with ARR1(1:N1):
!    !
!    !        ARR1(1) <= ARR2(I1) <= ARR2(I2) <= ARR1(N1)
!    !
!    ! ----------------------------------------------------------------
!
!    USE OMSAO_precision_module, ONLY: i4, r8
!    IMPLICIT NONE
!
!    ! ---------------
!    ! Input variables
!    ! ---------------
!    INTEGER (KIND=i4), INTENT (IN) :: n1, n2
!    REAL    (KIND=r8), DIMENSION (n1), INTENT (IN) :: arr1
!    REAL    (KIND=r8), DIMENSION (n2), INTENT (IN) :: arr2
!
!    ! ----------------
!    ! Output variables
!    ! ----------------
!    INTEGER (KIND=i4), INTENT (OUT) :: i1, i2
!
!    ! -----------------------------------------------------------------
!    ! ARRAY_LOCATE is the preferred way of determining the indices, but
!    ! this branch of the code has not yet been tested. tpk 09 Feb 2007
!    ! -----------------------------------------------------------------
!    i1 = -1 ; i2 = -1
!    i1 = MAXVAL ( MINLOC ( arr2(1:n2), MASK=(arr2(1:n2) >= arr1(1 )) ) )
!    i2 = MAXVAL ( MAXLOC ( arr2(1:n2), MASK=(arr2(1:n2) <= arr1(n1)) ) )
!    !CALL array_locate_r8 ( n2, arr2(1:n2), arr1( 1), 'GE', i1 )
!    !CALL array_locate_r8 ( n2, arr2(1:n2), arr1(n1), 'LE', i2 )
!
!    RETURN
!  END SUBROUTINE find_overlap

  SUBROUTINE fill_nonoverlap ( n, arr, i1, i2, filltype, fillval )

    ! ----------------------------------------------------------------
    !
    ! Given array ARR of dimension N and array positions I1 and I2 in
    ! ARR as returned from FIND_OVERLAP, fill in any non-overlapping
    ! segments either by constant extrapolation of endpoints or by
    ! using the value of FILLVAL.
    !
    ! ----------------------------------------------------------------

    USE OMSAO_precision_module, ONLY: i4, r8

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=9),                INTENT (IN)    :: filltype
    INTEGER (KIND=i4),                INTENT (IN)    :: n, i1, i2
    REAL    (KIND=r8),                INTENT (IN)    :: fillval

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL    (KIND=r8), DIMENSION (n), INTENT (INOUT) :: arr

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: low_end, hig_end

    low_end = 0.0_r8 ; hig_end = 0.0_r8

    SELECT CASE ( filltype )
    CASE ( 'endpoints' )
      low_end = arr(i1) ; hig_end = arr(i2)
    CASE ( 'fillvalue' )
      low_end = fillval ; hig_end = fillval
    END SELECT

    IF ( i1 > 1 ) arr (   1:i1-1) = low_end
    IF ( i2 < n ) arr (i2+1:n   ) = hig_end

    RETURN
  END SUBROUTINE fill_nonoverlap

!FUNCTION roundoff_r8 ( ndecim, r8value ) RESULT ( r8rounded )
!
!  USE OMSAO_precision_module, ONLY: i4, r8
!  IMPLICIT NONE
!
!  ! ---------------------------------------------------------------------
!  ! Explanation of FUNCTION arguments:
!  !
!  !    ndecim ........... number of significant digits to keep
!  !    r8value .......... DOUBLE PRECISION number to be truncated
!  !    r8rounded ........ truncated value
!  ! ---------------------------------------------------------------------
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4), INTENT (IN)    :: ndecim
!  REAL    (KIND=r8), INTENT (INOUT) :: r8value
!
!  ! --------------
!  ! Returned value
!  ! --------------
!  REAL    (KIND=r8) :: r8rounded
!
!  ! ---------------
!  ! Local variables
!  ! ---------------
!  INTEGER (KIND=i4) :: pow
!  REAL    (KIND=r8) :: tmpval
!
!  ! ------------------------------------------
!  ! Nothing to be done if we have a Zero value
!  ! ------------------------------------------
!  IF (  r8value == 0.0_r8 ) THEN
!    r8rounded = 0.0_r8 ; RETURN
!  END IF
!
!  ! ---------------------------------
!  ! Power of 10 of the original value
!  ! ---------------------------------
!  pow = INT ( LOG10 (ABS(r8value)), KIND=i4 )
!
!  ! ------------------------------------------------------
!  ! Remove original power of 10 and shift by NDECIM powers
!  ! ------------------------------------------------------
!  tmpval = r8value * 10.0_r8**(ndecim-pow)
!
!  ! ---------------------------------------------
!  ! Find the nearest INTEGER and undo power-shift
!  ! ---------------------------------------------
!  r8rounded = ANINT ( tmpval )  * 10.0_r8**(pow-ndecim)
!
!  RETURN
!END FUNCTION roundoff_r8
!
!FUNCTION roundoff_r4 ( ndecim, r4value ) RESULT ( r4rounded )
!
!  USE OMSAO_precision_module, ONLY: i4, r4
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4), INTENT (IN)    :: ndecim
!  REAL    (KIND=r4), INTENT (INOUT) :: r4value
!
!  INTEGER (KIND=i4) :: pow
!  REAL    (KIND=r4) :: r4rounded, tmpval
!
!  ! ------------------------------------------
!  ! Nothing to be done if we have a Zero value
!  ! ------------------------------------------
!  IF (  r4value == 0.0_r4 ) THEN
!    r4rounded = 0.0_r4 ; RETURN
!  END IF
!
!  ! ---------------------------------
!  ! Power of 10 of the original value
!  ! ---------------------------------
!  pow = INT ( LOG10 (ABS(r4value)), KIND=i4 )
!
!  ! ------------------------------------------------------
!  ! Remove original power of 10 and shift by NDECIM powers
!  ! ------------------------------------------------------
!  tmpval = r4value * 10.0_r4**(ndecim-pow)
!
!  ! ---------------------------------------------
!  ! Find the nearest INTEGER and undo power-shift
!  ! ---------------------------------------------
!  r4rounded = ANINT ( tmpval ) * 10.0_r4**(pow-ndecim)
!
!  RETURN
!END FUNCTION roundoff_r4
!
!SUBROUTINE roundoff_1darr_r8 ( ndecim, ndim, r8value )
!
!  USE OMSAO_precision_module, ONLY: i4, r8
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4),                   INTENT (IN)    :: ndecim, ndim
!  !REAL    (KIND=r8), DIMENSION (ndim), INTENT (INOUT) :: r8value
!  REAL    (KIND=r8), DIMENSION (:), INTENT (INOUT) :: r8value
!
!  INTEGER (KIND=i4) :: i
!  !REAL    (KIND=r8) :: roundoff_r8
!
!  DO i = 1, ndim
!    r8value(i) = ROUNDOFF_R8 ( ndecim, r8value(i) )
!  END DO
!
!  RETURN
!END SUBROUTINE roundoff_1darr_r8
!
!SUBROUTINE roundoff_1darr_r4 ( ndecim, ndim, r4value )
!
!  USE OMSAO_precision_module, ONLY: i4, r4
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4),                   INTENT (IN)    :: ndecim, ndim
!  !REAL    (KIND=r4), DIMENSION (ndim), INTENT (INOUT) :: r4value
!  REAL    (KIND=r4), DIMENSION (:), INTENT (INOUT) :: r4value
!
!  INTEGER (KIND=i4) :: i
!  !REAL    (KIND=r4) :: roundoff_r4
!
!  DO i = 1, ndim
!    r4value(i) = ROUNDOFF_R4 ( ndecim, r4value(i) )
!  END DO
!
!  RETURN
!END SUBROUTINE roundoff_1darr_r4
!
!SUBROUTINE roundoff_2darr_r8 ( ndecim, n1, n2, r8value )
!
!  USE OMSAO_precision_module, ONLY: i4, r8
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4),                    INTENT (IN)    :: ndecim, n1, n2
!  !REAL    (KIND=r8), DIMENSION (n1,n2), INTENT (INOUT) :: r8value
!  REAL    (KIND=r8), DIMENSION (:,:), INTENT (INOUT) :: r8value
!
!  INTEGER (KIND=i4) :: i, j
!
!  DO i = 1, n1
!    DO j = 1, n2
!      r8value(i,j) = roundoff_r8 ( ndecim, r8value(i,j) )
!    END DO
!  END DO
!
!  RETURN
!END SUBROUTINE roundoff_2darr_r8
!
!SUBROUTINE roundoff_2darr_r4 ( ndecim, n1, n2, r4value )
!
!  USE OMSAO_precision_module, ONLY: i4, r4
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4),                    INTENT (IN)    :: ndecim, n1, n2
!  !REAL    (KIND=r4), DIMENSION (n1,n2), INTENT (INOUT) :: r4value
!  REAL    (KIND=r4), DIMENSION (:,:), INTENT (INOUT) :: r4value
!  INTEGER (KIND=i4) :: i, j
!
!  DO i = 1, n1
!    DO j = 1, n2
!      r4value(i,j) = roundoff_r4 ( ndecim, r4value(i,j) )
!    END DO
!  END DO
!
!  RETURN
!END SUBROUTINE roundoff_2darr_r4
!
!SUBROUTINE roundoff_3darr_r8 ( ndecim, n1, n2, n3, r8value )
!
!  USE OMSAO_precision_module, ONLY: i4, r8
!  IMPLICIT NONE
!
!  ! ---------------
!  ! Input variables
!  ! ---------------
!  INTEGER (KIND=i4),                       INTENT (IN)    :: ndecim, n1, n2, n3
!  !REAL    (KIND=r8), DIMENSION (n1,n2,n3), INTENT (INOUT) :: r8value
!  REAL    (KIND=r8), DIMENSION (:,:,:), INTENT (INOUT) :: r8value
!
!  INTEGER (KIND=i4) :: i, j, k
!  !REAL    (KIND=r8) :: roundoff_r8
!
!  DO i = 1, n1
!    DO j = 1, n2
!      DO k = 1, n3
!        r8value(i,j,k) = ROUNDOFF_R8 ( ndecim, r8value(i,j,k) )
!      END DO
!    END DO
!  END DO
!
!  RETURN
!END SUBROUTINE roundoff_3darr_r8

  subroutine print_array (a, n)
    use OMSAO_precision_module, only: i4, r8
    implicit none
    integer (kind=i4), intent(in) :: n
    real (kind=r8), intent(in) :: a(n)

    integer(kind=i4) :: i
    do i=1, n
      write(*,*) "a[", i, "]=", a(i)
    enddo
  end subroutine print_array

!  SUBROUTINE get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid, lon0, lat0, mon0, & 
!       lon, lat, mon, nblon, nblat, nbmon, lonfrac, latfrac, monfrac, lonin, latin, monin)
!
!    USE OMSAO_precision_module, only: r8
!    IMPLICIT NONE
!
!    ! ======================
!    ! Input/Output variables
!    ! ======================
!    INTEGER, INTENT(IN) :: nlon, nlat, nmon
!    REAL (KIND=r8), INTENT(IN) :: lon0, lat0, mon0, lat, lon, mon, longrid, latgrid, mongrid 
!    INTEGER, INTENT(OUT) :: nblon, nblat, nbmon
!    INTEGER, DIMENSION(2), INTENT(OUT) :: latin, lonin, monin
!    REAL (KIND=r8), DIMENSION(2), INTENT(OUT) :: latfrac, lonfrac, monfrac 
!
!    ! ======================
!    ! Local variables
!    ! ======================  
!    REAL (KIND=r8) :: frac, lat_offset, lon_offset, mon_offset
!
!    lat_offset   = lat0   + latgrid / 2.0
!    lon_offset   = lon0   + longrid / 2.0
!    mon_offset   = mon0   + mongrid / 2.0
!
!    nblat = 2; frac = (lat - lat_offset) / latgrid + 1
!    latin(1) = INT(frac); latin(2) = latin(1) + 1
!    latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
!    IF (latin(1) == 0)   THEN    
!      latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
!    ENDIF
!
!    IF (latin(2) > nlat) THEN
!      latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
!    ENDIF
!
!    ! Circular in longitude direction
!    nblon = 2; frac = (lon - lon_offset) / longrid + 1
!    lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
!    lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
!    IF (lonin(1) == 0)   lonin(1) = nlon
!    IF (lonin(2) > nlon) lonin(2) = 1
!
!    ! Circular in year
!    nbmon = 2; frac = (mon - mon_offset) / mongrid + 1
!    monin(1) = INT(frac); monin(2) = monin(1) + 1
!    monfrac(1) = monin(2) - frac; monfrac(2) = 1.0 - monfrac(1)
!    IF (monin(1) == 0)   monin(1) = nmon
!    IF (monin(2) > nmon) monin(2) = 1
!
!    RETURN  
!
!  END SUBROUTINE get_gridfrac1

  !>Calculate relative azimuth angle
  !------------------------------------------------------------------------
  !
  !> @param[in]   saa  Solar Azimuth Angle (degrees)
  !> @param[in]   vaa  Viewing Azimuth Angle (degrees)
  !> @param[out]  raa  Relative Azimuth Angle = SAA + 180 - VAA
  !
  !> @author E. O'Sullivan  Decmber 2017
  !------------------------------------------------------------------------
  function calc_relaz_angle (saa,vaa) result (raa)
    use OMSAO_parameters_module, only: r4_missval

    implicit none

    real (kind=4), intent(in) :: saa,vaa
    real (kind=4) :: raa

    if (saa == r4_missval .or. vaa == r4_missval) then
      raa = r4_missval
      return
    endif

    raa = 180.0 + saa - vaa
    if (raa .gt. 180.0) then
      raa = raa - 360.0
    else if (raa .lt. -180.0) then
      raa = raa + 360.0
    endif

    return

  end function calc_relaz_angle

END MODULE sao_pge_utils
