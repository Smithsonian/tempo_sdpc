!
module m_utilities

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY : eoi_str
  USE OMSAO_errstat_module, ONLY : file_read_ok, file_read_failed, &
         file_read_eof
  USE OMSAO_variables_module, ONLY : maxchlen

  public check_for_endofinput, skip_to_filemark, gome_check_read_status, &
         get_substring, string2index, upper_case, & 
         get_doy, day_of_year,find_pos, signdp, &
         get_monfrac, get_latfrac, get_gridfrac, get_gridfrac1

  !private convert_16bits_to_2bytes,  &
         ! convert_8bits_to_byte
         !, year_month_day, &
         ! upper_case, utc_julian_date_and_time, 

contains

  FUNCTION day_of_year ( year, month, day ) RESULT ( jday )

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
    IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. &
         MOD(year,400) == 0 ) .and. month .gt. 2) jday = jday+1

    RETURN
  END FUNCTION day_of_year

  SUBROUTINE get_doy (theyear, themon,  theday, thedoy)

    IMPLICIT NONE

    ! --------------------------
    ! Input/Output variables
    ! --------------------------
    INTEGER, INTENT (IN)   :: theyear, themon, theday
    INTEGER, INTENT (OUT)  :: thedoy

    ! Local variables
    INTEGER, DIMENSION(12) :: ndays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    IF (MOD(theyear, 4) == 0) ndays(2) = 29
    IF (themon == 1) THEN
      thedoy = theday
    ELSE 
      thedoy = SUM(ndays(1:themon-1)) + theday
    ENDIF

    RETURN
  END SUBROUTINE get_doy

  SUBROUTINE timestamp (curtime )
    !
    !*******************************************************************************
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !
    !  Example:
    !
    !    May 31 2001   9:45:54.872 AM
    !
    !  Modified:
    !
    !    31 May 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none
    !
    character (len=24), INTENT(OUT) :: curtime
    integer d
    character ( len = 8 ) date
    integer h
    integer m
    integer mm
    character ( len = 3 ), parameter, dimension(12) :: month = (/ &
         'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)
    integer n
    integer s
    character ( len = 10 )  time
    integer values(8)
    integer y
    character ( len = 5 ) zone
    !
    call date_and_time ( date, time, zone, values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    !if ( h < 12 ) then
    !  ampm = 'AM'
    !else if ( h == 12 ) then
    !  if ( n == 0 .and. s == 0 ) then
    !    ampm = 'Noon'
    !  else
    !    ampm = 'PM'
    !  end if
    !else
    !  h = h - 12
    !  if ( h < 12 ) then
    !    ampm = 'PM'
    !  else if ( h == 12 ) then
    !    if ( n == 0 .and. s == 0 ) then
    !      ampm = 'Midnight'
    !    else
    !      ampm = 'AM'
    !    end if
    !  end if
    !end if

    WRITE ( curtime, '(a3,1x,i2.2,1x,i4,1x,i2.2,a1,i2.2,a1,i2.2,a1,i3.3)' ) &
         TRIM ( month(m) ), d, y, h, ':', n, ':', s, '.', mm

    RETURN
  END SUBROUTINE timestamp

  SUBROUTINE check_for_endofinput ( iostring, yn_eoi )

    IMPLICIT NONE

    ! ==============
    ! Input variable
    ! ==============
    CHARACTER (LEN=*), INTENT (IN) :: iostring

    ! ===============
    ! Output variable
    ! ===============
    LOGICAL, INTENT (OUT) :: yn_eoi

    yn_eoi = .FALSE.
    IF ( TRIM(ADJUSTL(iostring)) == eoi_str ) yn_eoi = .TRUE.

    RETURN
  END SUBROUTINE check_for_endofinput

  SUBROUTINE gome_check_read_status ( ios, file_read_stat )

    IMPLICIT NONE

    INTEGER, INTENT (IN)  :: ios
    INTEGER, INTENT (OUT) :: file_read_stat

    SELECT CASE ( ios )
    CASE ( :-1 )  
  file_read_stat = file_read_eof
    CASE (   0 )  
  file_read_stat = file_read_ok
    CASE DEFAULT  
  file_read_stat = file_read_failed
    END SELECT

    RETURN
  END SUBROUTINE gome_check_read_status


  SUBROUTINE skip_to_filemark ( funit, lm_string, lastline, file_read_stat )

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER,           INTENT (IN) :: funit
    CHARACTER (LEN=*), INTENT (IN) :: lm_string

    ! ================
    ! Output variables
    ! ================
    INTEGER,           INTENT (OUT) :: file_read_stat
    CHARACTER (LEN=*), INTENT (OUT) :: lastline

    ! ===============
    ! Local variables
    ! ===============
    INTEGER                  :: lmlen, ios!, iline
    CHARACTER (LEN=maxchlen) :: tmpline

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

    CALL gome_check_read_status ( ios, file_read_stat )

    RETURN
  END SUBROUTINE skip_to_filemark


  SUBROUTINE get_substring ( string, sstart, substring, nsubstring, eostring )

    IMPLICIT NONE

    ! =========================================================================
    ! Given string STRING, extract first space- or comma-delimited substring
    ! beginning on or after position SSTART, and return its value, SUBSTRING 
    ! and length, NSUBSTRING; update STRING and SSTART to remove the substring.
    ! 
    ! F90 version of the original GET_TOKEN by J. Lavanigno
    ! =========================================================================

    ! NSUBSTRING represents the token's length without trailing blanks. It's
    ! set to 0 if no substring was found in STRING: this can be used to 
    ! determine when to stop looking.
    !
    ! This routine makes no attempt to detect or handle the case in which
    ! the token is bigger than the token buffer.  It's assumed that the
    ! caller will declare line and token to be the same size so that this
    ! won't ever happen.

    ! ==================
    ! Modified arguments.
    ! ==================
    CHARACTER (LEN = *), INTENT (INOUT) :: string
    INTEGER,             INTENT (INOUT) :: sstart

    ! =================
    ! Output arguments.
    ! =================
    CHARACTER (LEN =LEN(string)), INTENT (OUT) :: substring
    INTEGER,                      INTENT (OUT) :: nsubstring
    LOGICAL,                      INTENT (OUT) :: eostring

    ! ================
    ! Local arguments.
    ! ================
    CHARACTER :: char
    INTEGER   :: lstart, lend, nline


    ! ------------------
    ! Initialize outputs
    ! ------------------
    substring  = ' '  
    nsubstring = 0  
    eostring = .FALSE.

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
      eostring = .TRUE.  
      RETURN
    END IF

    ! --------------------------------------------------------
    ! Find first character in line that's not a blank or comma
    ! --------------------------------------------------------
    findchar: DO lstart = sstart, nline
      char = string(lstart:lstart)
      SELECT CASE ( char )
      CASE ( ' ' )
        IF ( lstart == nline ) THEN
          sstart = nline + 1
         RETURN  ! there are no further substrings in STRING
        END IF
      CASE ( ',' )
        IF ( lstart == nline ) THEN
          sstart = nline + 1
         RETURN  ! there are no further substrings in STRING
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

  SUBROUTINE string2index ( table, ntable, string, stridx )

    ! =====================================================
    ! Looks up STRING in character table TABLE of dimension
    ! NTABLE, and returns position STRIDX. Defaults to
    ! STRIDX = -1 if STRING is not found in TABLE.
    ! =====================================================

    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER,                               INTENT (IN) :: ntable
    CHARACTER (LEN=*), DIMENSION (ntable), INTENT (IN) :: table
    CHARACTER (LEN=*),                     INTENT (IN) :: string

    ! ================
    ! Output variables
    ! ================
    INTEGER, INTENT (OUT) :: stridx

    ! ===============
    ! Local variables
    ! ===============
    INTEGER :: i

    stridx = -1

    getidx: DO i = 1,  ntable
      IF ( TRIM(ADJUSTL(string)) == TRIM(ADJUSTL(table(i))) ) THEN
        stridx = i
        EXIT getidx
      END IF
    END DO getidx

    RETURN
  END SUBROUTINE string2index

  CHARACTER (LEN=maxchlen) FUNCTION int2string ( int_num, ni ) RESULT ( int_str )

  ! ===============================================================
  ! Converts INTEGER number INT_NUM to STRING INT_STR of length NI
  ! or the number of digits in INT_NUM, whatever is larger.
  ! ===============================================================

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),  INTENT (IN)  :: int_num, ni

  !! ---------------
  !! Result variable
  !! ---------------
  !CHARACTER (LEN=*) :: int_str

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  ! * Arrays containing indices for number strings in ASCII table
  INTEGER (KIND=i4),                   PARAMETER :: n = 10
  INTEGER (KIND=i4), DIMENSION(0:n-1)            :: aidx
  CHARACTER (LEN=1), DIMENSION(0:n-1), PARAMETER :: astr = (/"0","1","2","3","4","5","6","7","8","9" /)
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
     int_str = '?' ; RETURN
  END SELECT
 ! -------------------------------------------------------------
  ! We may want to create a string that is longer than the number
  ! of digits in INT_NUM. This will create leading "0"s.
  ! -------------------------------------------------------------
  nd = MAX ( nd, ni )

  ! ----------------------------------------------
  ! Convert the integer to string "digit by digit"
  ! ----------------------------------------------
  int_str = "" ; tmpint = int_num
  DO i = 1, nd
     ld = MOD ( tmpint, 10 )        ! Current last digit
     tmpint = ( tmpint - ld ) / 10  ! Remove current last digit from INT_STR
     k = nd - i + 1                 ! Position of current digit in INT_STR
     int_str(k:k) = ACHAR(aidx(ld)) ! Convert INTEGER digit to CHAR
  END DO

  RETURN
  END FUNCTION int2string


  REAL (KIND=KIND(1.0D0)) FUNCTION signdp ( x )

    IMPLICIT NONE

    REAL (KIND=dp), INTENT (IN) :: x

    signdp = 0.0_dp
    IF ( x < 0.0_dp ) THEN
      signdp = -1.0_dp
    ELSE
      signdp = +1.0_dp
    END IF

    RETURN
  END FUNCTION signdp

  SUBROUTINE reverse ( inarr, num )
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)

    INTEGER, INTENT(IN) :: num
    INTEGER             :: i
    REAL (KIND=dp), DIMENSION(1: num), INTENT(INOUT) :: inarr
    REAL (KIND=dp), DIMENSION(1: num)                :: temp

    DO i = 1, num
      temp(i) = inarr(num - i + 1)
    ENDDO
    inarr = temp

    RETURN
  END SUBROUTINE reverse

  SUBROUTINE find_pos (fwave,nf, cwave,nc, pos)

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY : winlim, numwin
    IMPLICIT NONE

    ! ===============
    ! Input variables
    ! ===============
    INTEGER,                        INTENT (IN)         :: nc, nf
    REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
    REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
    
    ! ===============
    ! Output  variables
    ! ===============
    INTEGER, DIMENSION (nc) :: pos
    ! ===============
    ! Local variables
    ! ===============
    REAL (KIND=dp) :: dfw, dcw, temp
    INTEGER        :: i, iwin, npos, fidx, fidxc, lidx, lidxc, midx

    fidx = 1; fidxc = 1
    dcw  = cwave(2) - cwave(1)
    dfw  = fwave(2) - fwave(1) !xliu, 10/22/2009
    npos = 0
    DO iwin = 1, numwin
      IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
      ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
      ENDIF
      IF (fwave(1) > winlim(iwin, 2)) CYCLE         
      DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx
        pos(i) = midx
     !   print *, iwin, i, cwave(i), fwave(midx)
      ENDDO

      fidx = lidx + 1; fidxc = lidxc + 1
    ENDDO
    RETURN
  END SUBROUTINE find_pos

 SUBROUTINE get_monfrac(nmon, mon, day, nbmon, monfrac, monin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nmon, mon, day
  INTEGER, INTENT(OUT)                      :: nbmon
  INTEGER, DIMENSION(2), INTENT(OUT)        :: monin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: monfrac

  IF (day <= 15) THEN
     monin(1) = mon - 1
     IF (monin(1) == 0) monin(1) = 12
     monin(2) = mon
     monfrac(1) = (15.0 - day) / 30.0
     monfrac(2) = 1.0 - monfrac(1)
  ELSE
     monin(2) = mon + 1
     IF (monin(2) == 13) monin(2) = 1
     monin(1) = mon
     monfrac(2) = (day - 15) / 30.0
     monfrac(1) = 1.0 - monfrac(2)
  ENDIF
     nbmon=2
  END SUBROUTINE get_monfrac

  SUBROUTINE get_latfrac( nlat, latgrid, lat0, lat,  nblat, latfrac,latin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlat
  REAL (KIND=dp), INTENT(IN)                :: lat0, lat, latgrid
  INTEGER, INTENT(OUT)                      :: nblat
  INTEGER, DIMENSION(2), INTENT(OUT)        :: latin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: latfrac

  ! ======================
  ! Local variables
  ! ======================  
  REAL (KIND=dp) :: frac, lat_offset

  lat_offset   = lat0 + latgrid / 2.0
  nblat = 2; frac = (lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)

  IF (latin(1) == 0)   THEN
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF

  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF
  RETURN
  END SUBROUTINE get_latfrac

  SUBROUTINE get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
  lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlon, nlat
  REAL (KIND=dp), INTENT(IN)                :: lon0, lat0, lat, lon,longrid, latgrid
  INTEGER, INTENT(OUT)                      :: nblon, nblat
  INTEGER, DIMENSION(2), INTENT(OUT)        :: latin, lonin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: latfrac, lonfrac

  ! ======================
  ! Local variables
  ! ======================  
  REAL (KIND=dp) :: frac, lat_offset, lon_offset

  lat_offset   = lat0 + latgrid / 2.0
  lon_offset   = lon0 + longrid  / 2.0

  nblat = 2; frac = (lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  IF (latin(1) == 0)   THEN
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF

  ! Circular in longitude direction
  nblon = 2; frac = (lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1

  RETURN

  END SUBROUTINE get_gridfrac

  SUBROUTINE get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid,lon0, lat0, mon0, &
     lon, lat, mon, nblon, nblat, nbmon, lonfrac, latfrac, monfrac,lonin, latin, monin)

  USE OMSAO_precision_module
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlon, nlat, nmon
  REAL (KIND=dp), INTENT(IN)                :: lon0, lat0, mon0, lat,lon, mon, longrid, latgrid, mongrid
  INTEGER, INTENT(OUT)                      :: nblon, nblat, nbmon
  INTEGER, DIMENSION(2), INTENT(OUT)        :: latin, lonin, monin
  REAL (KIND=dp), DIMENSION(2), INTENT(OUT) :: latfrac, lonfrac, monfrac

  ! ======================
  ! Local variables
  ! ======================  
  REAL (KIND=dp) :: frac, lat_offset, lon_offset, mon_offset

  lat_offset   = lat0   + latgrid / 2.0
  lon_offset   = lon0   + longrid / 2.0
  mon_offset   = mon0   + mongrid / 2.0

  nblat = 2; frac = (lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  IF (latin(1) == 0)   THEN
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF

  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF

  ! Circular in longitude direction
  nblon = 2; frac = (lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1

  ! Circular in year
  nbmon = 2; frac = (mon - mon_offset) / mongrid + 1
  monin(1) = INT(frac); monin(2) = monin(1) + 1
  monfrac(1) = monin(2) - frac; monfrac(2) = 1.0 - monfrac(1)
  IF (monin(1) == 0)   monin(1) = nmon
  IF (monin(2) > nmon) monin(2) = 1

  RETURN

  END SUBROUTINE get_gridfrac1
  !  Unused?
!
!  SUBROUTINE utc_julian_date_and_time ( year, month, day, julday, hour, &
!       minute, second )
!
!    ! -------------------------------------------------------------
!    ! Converts current date and time to UTC time and "Julian" date,
!    ! i.e., day of the year.
!    ! -------------------------------------------------------------
!
!    USE OMSAO_precision_module, ONLY: i4!, r4
!    IMPLICIT NONE
!
!    ! ---------------
!    ! RESULT
!    ! ---------------
!    INTEGER (KIND=i4), INTENT (OUT) :: year, month, day, julday, hour, &
!         minute, second
!
!    !INTEGER (KIND=i4), DIMENSION(8), INTENT (OUT) :: utc_juldate
!
!    ! ---------------
!    ! Local variables
!    ! ---------------
!    ! * Arguments for DATE_AND_TIME, and position indices for Year, Month, Day, Zone, Hour, Minute, Seconds
!    CHARACTER (LEN= 8)                 :: date
!    CHARACTER (LEN=10)                 :: time
!    CHARACTER (LEN= 5)                 :: zone
!    INTEGER   (KIND=i4), DIMENSION (8) :: date_vector
!    INTEGER   (KIND=i4), PARAMETER     :: y_idx=1, m_idx=2, d_idx=3, &
!         z_idx=4, hh_idx=5, mm_idx=6, ss_idx=7
!    ! * Some MAX values for Minutes and Hours 
!    INTEGER   (KIND=i4), PARAMETER     :: max_mi = 60, max_hr = 24, &
!         max_dy = 31, max_mo = 12
!    ! * Other local variables
!    INTEGER   (KIND=i4)    :: del_mm, del_hh, max_julday, yyy, mmm, ddd
!
!
!    ! ----------------------------------------------------------
!    ! First copy the input variable to the output variable. That
!    ! saves work on premature return.
!    ! ----------------------------------------------------------
!    CALL DATE_AND_TIME ( date, time, zone, date_vector )
!
!    ! -----------------------------------------------
!    ! Initialize all output variables (except JULDAY)
!    ! -----------------------------------------------
!    year   = date_vector(y_idx )
!    month  = date_vector(m_idx )
!    day    = date_vector(d_idx )
!    hour   = date_vector(hh_idx)
!    minute = date_vector(mm_idx)
!    second = date_vector(ss_idx)
!
!    ! ----------------------------------------------------------------------
!    ! Next, find the day of the year. At this point, we also compute the
!    ! maximum number of days per year. This will come in handy further down.
!    ! ----------------------------------------------------------------------
!    julday     = day_of_year ( year, month,  day    )
!    max_julday = day_of_year ( year, max_mo, max_dy )
!
!    ! -----------------------------------------------
!    ! Check if there is any difference to UTC at all.
!    ! -----------------------------------------------
!    IF ( ABS(date_vector(z_idx)) == 0 ) RETURN
!
!    ! ------------------------------------------------------------
!    ! Apply difference in UTC time to local time vector. This is
!    ! given in minutes, thus apply to minutes index. 
!    !
!    ! REMEMBER: This is the DIFFERENCE to UTC, so we have to
!    !           SUBTRACT whatever difference from the current time.
!    ! -------------------------------------------------------------
!    del_mm = MOD ( date_vector(z_idx), max_mi )
!    minute = minute - del_mm
!
!    ! ----------------------------------------
!    ! Compute hour difference, return if none.
!    ! ----------------------------------------
!    del_hh =   ( date_vector(z_idx) - del_mm ) / max_mi
!    IF ( del_hh == 0 ) RETURN
!
!    ! --------------------------------------------------------
!    ! Apply hour difference and check for/apply day difference
!    ! --------------------------------------------------------
!    hour = hour - del_hh
!
!    ! -------------------------------------------------------------
!    ! If we reach here, we may have had a day difference. We could
!    ! make our lives miserable by checking for all the months of 
!    ! the year. But what we ultimately want is the day of the year
!    ! (a proxy for Julian day), and that makes life a lot easier.
!    ! -------------------------------------------------------------
!    SELECT CASE ( hour )
!    CASE ( 0:23 )
!      RETURN
!    CASE ( :-1 )
!      hour = max_hr + hour
!      day  = day - 1
!    CASE ( 24: )
!      hour = hour - max_hr
!      day  = day  + 1
!    END SELECT
!
!    ! --------------------------------------------------------------
!    ! Finally, if we ever reach here, we've had a year difference.
!    ! Apply correction, compute the new day of the year, and return.
!    ! --------------------------------------------------------------
!    yyy = year 
! mmm = month 
! ddd = day                     ! save old values
!    CALL year_month_day ( yyy, mmm, ddd, year, month, day )  ! overwrite with new ones
!    julday = day_of_year ( year, month, day )                ! recompute day-of-year
!
!    RETURN
!  END SUBROUTINE utc_julian_date_and_time
!  Unused
!
!  SUBROUTINE year_month_day ( year, month, day, newyear, newmonth, newday )
!
!    USE OMSAO_precision_module, ONLY: i4
!    IMPLICIT NONE
!
!    ! ----------------------
!    ! Input/Output variables
!    ! ----------------------
!    INTEGER (KIND=i4), INTENT (IN)    :: year, month, day
!    INTEGER (KIND=i4), INTENT (INOUT) :: newyear, newmonth, newday
!
!    ! ------------------------------
!    ! Local variables and parameters
!    ! ------------------------------
!    INTEGER (KIND=i4), PARAMETER :: n_month = 12
!    INTEGER (KIND=i4), DIMENSION (n_month), PARAMETER :: &
!         days_per_month    = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
!         ly_days_per_month = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
!
!    INTEGER (KIND=i4), DIMENSION (n_month) :: dpm
!
!    ! -------------------------
!    ! Check for leap year
!    ! ------------------------------------
!    ! * Divisible by 4:   leap year
!    ! * Divisible by 100: not a leap year
!    ! * Divisible by 400: leap year
!    ! ------------------------------------
!    IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. &
!         MOD(year,400) == 0 ) ) THEN
!      dpm = days_per_month
!    ELSE
!      dpm = ly_days_per_month
!    END IF
!
!    newyear = year  
!  newmonth = month  
!  newday = day
!
!    ! -------------------------------------------------------
!    ! If DAY is within bounds there's nothing to do
!    ! (MONTH hasn't been adjusted yet in the calling program)
!    ! -------------------------------------------------------
!    IF ( day > 0 .AND. day <= dpm(month) ) RETURN
!
!    ! ------------------------------------------------------------------------
!    ! If we reach here, DAY must be either too large or <= 0. The crucial
!    ! variable now is MONTH, since it might lead to an additional change 
!    ! in YEAR.
!    ! ------------------------------------------------------------------------
!    SELECT CASE ( month )
!    CASE ( 2:11 )   ! The easiest case: change only MONTH and DAY
!      IF ( day <= 0 ) THEN
!        newmonth = month - 1
!        newday   = dpm(newmonth) + day  ! DAY is negative
!      ELSE
!        newmonth = month + 1
!        newday   = day - dpm(month)
!      END IF
!    CASE ( 1 )      ! January requires a YEAR change if DAY <= 0
!      IF ( day <= 0 ) THEN
!        newmonth = 12
!        newday   = dpm(newmonth) + day  ! DAY is negative
!        newyear  = year - 1
!      ELSE
!        newmonth = month + 1
!        newday   = day - dpm(month)
!      END IF
!    CASE ( 12 )      ! December requires a YEAR change if DAY > 31
!      IF ( day <= 0 ) THEN
!        newmonth = month - 1
!        newday   = dpm(newmonth) + day  ! DAY is negative
!      ELSE
!        newmonth = 1
!        newday   = day - dpm(month)
!        newyear  = year + 1
!      END IF
!    END SELECT
!
!    RETURN
!  END SUBROUTINE year_month_day


!   Unused
!
  FUNCTION upper_case ( mixstring ) RESULT ( upcase )

    ! =============================================
    ! Function to convert strings to all upper case
    ! =============================================

    IMPLICIT NONE

    ! ------------------------------------------------------------
    CHARACTER (LEN=*), INTENT (IN) :: mixstring   ! Input string
    ! ------------------------------------------------------------

    ! -------------------------------------------------------------
    CHARACTER (LEN=len(mixstring))              :: upcase      ! Result string
    ! -------------------------------------------------------------

    ! -------------------------------------------------------------
    INTEGER :: i!, slen                            ! Local variables
    ! -------------------------------------------------------------

    DO i = 1, LEN(mixstring)
      SELECT CASE ( ICHAR(mixstring(i:i)) )
      CASE ( 97:122 )
        upcase(i:i) = ACHAR(ICHAR(mixstring(i:i))-32)
      CASE DEFAULT
        upcase(i:i) = mixstring(i:i)
      END SELECT
    END DO

    RETURN
  END FUNCTION upper_case
   
end module m_utilities
