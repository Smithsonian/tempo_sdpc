!
module prepare_atmosphere

  public get_gridfrac1, get_geoschem_o3mean, get_geoschem_o3std, get_spres, &
       get_tpres, get_toz, get_ncept, get_v8prof, get_v8temp, get_mipasig2t, &
       get_surfalt, get_ncepreso_surfalt, get_gridfrac, &
       get_finereso_surfalt, get_geoschem_o31, get_logan_clima, &
       get_geoschem_hcho, get_bro, get_no2, get_geoschem_so2, &
       get_mlso3prof, get_mlso3prof_single, get_tomsv8_clima, &
       get_normtoz, get_mipasig2o3
  private !get_ecmwft

  integer, parameter, private :: max_pathlen = 1024

contains
  ! ************************************************************************
  ! Author:  xiong liu
  ! Date  :  July 24, 2003
  ! Purpose: Routine to read atmospheric surface and tropopause pressure, 
  !          temperature, ozone, trace gases, surface altitdue and so on.
  ! ************************************************************************
  SUBROUTINE get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)           :: nlon, nlat
    REAL (KIND=dp), INTENT(IN)    :: lon0, lat0, lat, lon, longrid, latgrid
    INTEGER, INTENT(OUT)          :: nblon, nblat
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


  SUBROUTINE get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid, lon0, lat0, mon0, &
       lon, lat, mon, nblon, nblat, nbmon, lonfrac, latfrac, monfrac, lonin, latin, monin)

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                       :: nlon, nlat, nmon
    REAL (KIND=dp), INTENT(IN)                :: lon0, lat0, mon0, lat, lon, mon, longrid, latgrid, mongrid
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

  SUBROUTINE get_geoschem_o3mean(month, lon, lat, ps, ozprof, nz, ntp)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)         :: month, nz, ntp
    REAL (KIND=dp), INTENT(IN)  :: lat, lon
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER               :: nlat=18, nlon=12, nalt=19
    INTEGER                          :: errstat, i, j, k, nblat, nblon, nalt0, ntp0
    REAL (KIND=dp), PARAMETER        :: latgrid=10.0, longrid=30.0, lon0=-180.0, lat0=-90.0

    INTEGER, DIMENSION(2)            :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)     :: latfrac, lonfrac
    REAL (KIND=dp), DIMENSION(nalt)  :: gprof
    REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz

    ! Saved variables
    REAL (KIND=dp), DIMENSION(nlon, nlat, nalt), SAVE :: geosoz
    LOGICAL                                    , SAVE :: first = .TRUE.

    REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
    CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb',&
         'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=max_pathlen)                         :: geosfile

    ! Correct coordinates
    REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
         .987871d0, .954730d0, .905120d0, .845000d0, .78d0,    .710000d0,   &
         .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
         .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
         .102098d0/)

    IF (first) THEN
      geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(month) // '_o3_mean.dat'

      OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
      READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

    gprof = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO

    geospres = pres * ps(0)

    cumoz = 0.0
    DO i = 1, nalt
      cumoz(i) = cumoz(i-1) + gprof(i) * 1000.0 / 1.25 * &
           (geospres(i-1) - geospres(i)) / 1013.25 
      IF (ANY(geosoz(lonin(1:nblon), latin(1:nblat), i) <= 0.0)) THEN
        j = i - 1; EXIT
      ELSE
        j = i
      ENDIF
    ENDDO
    nalt0 = j

    DO i = 1, ntp
      IF (ps(i) < geospres(nalt0)) THEN
        ntp0 = i - 1; EXIT
      ENDIF
    ENDDO

    CALL BSPLINE(geospres, cumoz, nalt0+1, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
    tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)     
    ozprof(1:ntp0) =  tempoz(1:ntp0) !* SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) *

    RETURN  
  END SUBROUTINE get_geoschem_o3mean

  SUBROUTINE get_geoschem_o3std(month, lon, lat, ps, ozprof, nz, ntp)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)         :: month, nz, ntp
    REAL (KIND=dp), INTENT(IN)  :: lat, lon
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER               :: nlat=18, nlon=12, nalt=19
    INTEGER                          :: errstat, i, j, k, nblat, nblon, nalt0, ntp0
    REAL (KIND=dp), PARAMETER        :: latgrid=10.0, longrid=30.0, lon0=-180.0, lat0=-90.0

    INTEGER, DIMENSION(2)            :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)     :: latfrac, lonfrac
    REAL (KIND=dp), DIMENSION(nalt)  :: gprof
    REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz

    ! Saved variables
    REAL (KIND=dp), DIMENSION(nlon, nlat, nalt), SAVE :: geosoz
    LOGICAL                                    , SAVE :: first = .TRUE.

    REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
    CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb',&
         'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=max_pathlen)                         :: geosfile

    ! Correct coordinates
    REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
         .987871d0, .954730d0, .905120d0, .845000d0, .78d0,    .710000d0,   &
         .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
         .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
         .102098d0/)

    IF (first) THEN
      geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(month) // '_o3_std.dat'

      OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
      READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

    gprof = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO
    geospres = pres * ps(0)

    cumoz = 0.0
    DO i = 1, nalt
      cumoz(i) = cumoz(i-1) + gprof(i) * 1000.0 / 1.25 * &
           (geospres(i-1) - geospres(i)) / 1013.25 
      IF (ANY(geosoz(lonin(1:nblon), latin(1:nblat), i) <= 0.0)) THEN
        j = i - 1; EXIT
      ELSE
        j = i
      ENDIF
    ENDDO
    nalt0 = j

    DO i = 1, ntp
      IF (ps(i) < geospres(nalt0)) THEN
        ntp0 = i - 1; EXIT
      ENDIF
    ENDDO

    CALL BSPLINE(geospres, cumoz, nalt+1, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
    tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)    
    ozprof(1:ntp0) =  tempoz(1:ntp0) 

    RETURN  
  END SUBROUTINE get_geoschem_o3std

  ! ====================================================
  ! Obtain NCAR/NCEP 12pm surface pressure (mb) 
  !   for  each 2.5 by 2.5 region
  ! If no data is available, then use mean surface
  !    pressure from all years
  ! =====================================================
  SUBROUTINE get_spres(year, month, day, lon, lat, spres)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)           :: month, year, day
    REAL (KIND=dp), INTENT(IN)    :: lon, lat
    REAL (KIND=dp), INTENT(OUT)   :: spres

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER             :: nlat=72, nlon=144
    REAL (KIND=dp), PARAMETER      :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
    INTEGER                        :: i, j,  nblat, nblon
    LOGICAL                        :: file_exist
    CHARACTER (LEN=2)              :: monc, yrc, dayc
    CHARACTER (LEN=max_pathlen)            :: spres_fname
    INTEGER, DIMENSION(2)          :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

    INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbspres
    LOGICAL, SAVE                        :: first = .TRUE.

    IF (first) THEN
      WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
      WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
      WRITE(yrc, '(I2.2)')  MOD(year, 100) ! from 1997 to '97'

      spres_fname =TRIM(ADJUSTL(atmdbdir)) // 'nspres/spres' // yrc // monc // dayc // '.dat'

      ! Determine if file exists or not
      INQUIRE (FILE= spres_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'Warning: no surface pressure file found, use monthly mean!!!'
        spres_fname = TRIM(ADJUSTL(atmdbdir)) // 'nspres/spresavg' // monc // '.dat'
      ENDIF

      OPEN (UNIT = atmos_unit, file = spres_fname, status = 'unknown')
      READ (atmos_unit, '(144I4)') ((glbspres(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    spres = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        spres = spres + glbspres(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_spres

  ! ====================================================
  ! Obtain NCAR/NCEP 12pm tropopause pressure (mb) 
  !   for  each 2.5 by 2.5 region
  ! If no data is available, then use mean surface
  !    pressure from all years
  ! =====================================================
  SUBROUTINE get_tpres(year, month, day, lon, lat, tpres)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)           :: month, year, day
    REAL (KIND=dp), INTENT(IN)    :: lon, lat
    REAL (KIND=dp), INTENT(OUT)   :: tpres

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER             :: nlat=72, nlon=144
    REAL (KIND=dp), PARAMETER      :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
    INTEGER                        :: i, j,  nblat, nblon
    LOGICAL                        :: file_exist
    CHARACTER (LEN=2)              :: monc, yrc, dayc
    CHARACTER (LEN=max_pathlen)            :: tpres_fname
    INTEGER, DIMENSION(2)          :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac


    INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbtpres
    LOGICAL, SAVE                        :: first = .TRUE.

    IF (first) THEN
      WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
      WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
      WRITE(yrc, '(I2.2)')  MOD(year, 100) ! from 1997 to '97'

      tpres_fname =TRIM(ADJUSTL(atmdbdir)) // 'ntpres/tpres' // yrc // monc // dayc // '.dat'

      ! Determine if file exists or not
      INQUIRE (FILE= tpres_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'Warning: no tropopause pressure file found, use monthly mean!!!'
        tpres_fname = TRIM(ADJUSTL(atmdbdir)) // 'ntpres/tpresavg' // monc // '.dat'
      ENDIF

      OPEN (UNIT = atmos_unit, file = tpres_fname, status = 'unknown')
      READ (atmos_unit, '(144I3)') ((glbtpres(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    tpres = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        tpres = tpres + glbtpres(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_tpres


  ! ===================================================
  ! Obtain EP TOMS monthly mean total ozone (DU) for
  !    each 1.25 by 1 region
  ! If no data is available, then use mean total ozone 
  !    from all years
  ! ====================================================
  SUBROUTINE get_toz(year, month, day, lon, lat, toz)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)         :: month, year, day
    REAL (KIND=dp),INTENT(IN)   :: lon, lat
    REAL (KIND=dp), INTENT(OUT) :: toz

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER           :: nlat=180, nlon=288
    REAL (KIND=dp), PARAMETER    :: longrid = 1.25, latgrid = 1.0, lon0=-180.0, lat0=-90.0
    CHARACTER (LEN=2)            :: monc, yrc, dayc
    CHARACTER (LEN=max_pathlen)          :: toz_fname

    INTEGER                      :: i, j, nblat, nblon
    LOGICAL                      :: file_exist
    INTEGER, DIMENSION(2)        :: latin, lonin
    REAL (KIND=dp), DIMENSION(2) :: latfrac, lonfrac

    INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbtoz
    LOGICAL, SAVE                        :: first = .TRUE.

    IF (first) THEN
      WRITE(dayc, '(I2.2)') day             ! from 9 to '09'  
      WRITE(monc, '(I2.2)') month           ! from 9 to '09'  
      WRITE(yrc,  '(I2.2)') MOD(year, 100)  ! from 1997 to '97'

      toz_fname = TRIM(ADJUSTL(atmdbdir)) // 'eptoz/ep' // yrc // monc // '.dat'

      ! Determine if file exists or not
      INQUIRE (FILE= toz_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'Warning: no EP O3 file found, use monthly mean!!!'
        toz_fname = TRIM(ADJUSTL(atmdbdir)) // 'eptoz/avgep' // monc // '.dat'
      ENDIF

      OPEN (UNIT = atmos_unit, file=toz_fname, status = 'unknown')
      DO i = 1, 3
        READ (atmos_unit, '(A)')
      END DO
      READ (atmos_unit, *) ((glbtoz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    toz = 0.0
    DO i = 1, nblon
      DO j = 1, nblat
        IF (glbtoz(lonin(i), latin(j)) > 0) &
             toz = toz + glbtoz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_toz

!  Unused?
!
!  ! Obtain ECMWF temperature profile
!  SUBROUTINE get_ecmwft(year, month, day, lon, lat, ecmwft)
!
!    USE OMSAO_precision_module
!    USE OMSAO_variables_module, ONLY: atmdbdir
!    USE ozprof_data_module,     ONLY: atmos_unit
!    USE OMSAO_errstat_module
!    IMPLICIT NONE
!
!    ! ======================
!    ! Input/Output variables
!    ! ======================
!    INTEGER, PARAMETER          :: nlecm = 23
!    INTEGER, INTENT(IN)         :: month, year, day
!    REAL (KIND=dp), INTENT(IN)  :: lon, lat
!    REAL (KIND=dp), DIMENSION(nlecm), INTENT(OUT) :: ecmwft
!
!    ! ======================
!    ! Local variables
!    ! ======================
!    INTEGER, PARAMETER           :: nlat=72, nlon=144, nalt=23
!    REAL (KIND=dp), PARAMETER    :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
!    CHARACTER (LEN=2)            :: yrc, monc, dayc
!    CHARACTER (LEN=max_pathlen)          :: ecmwft_fname, ncep_fname
!    INTEGER                      :: i, j, k, nblat, nblon
!    INTEGER, DIMENSION(2)        :: latin, lonin
!    REAL (KIND=dp), DIMENSION(2) :: latfrac, lonfrac
!    LOGICAL                      :: file_exist
!
!    INTEGER, SAVE, DIMENSION(nlon, nlat, nalt) :: glbecmwft
!    LOGICAL, SAVE                              :: first = .TRUE.
!
!    IF (first) THEN
!      WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
!      WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
!      WRITE(yrc,  '(I2.2)') MOD(year, 100) ! from 1997 to '97'
!
!      ! use ECMWF
!      IF (year <= 2001) THEN      
!        ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwft' // yrc // monc // dayc // '.dat'      
!        ! Determine if file exists or not
!        INQUIRE (FILE= ecmwft_fname, EXIST= file_exist)
!        IF (.NOT. file_exist) THEN
!          WRITE(www_lun, *) 'Warning: no T profile file found, use monthly mean!!!'
!          ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
!        ENDIF
!        OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
!        READ (atmos_unit, '(144i3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
!        CLOSE(atmos_unit)
!      ELSE  ! Use NCEP for up to 10 mb and ECMWFT average for up between 10 and 1 mb
!        ! ECMWFT average between 10mb and 1mb (7, 5, 3, 2, 1)', other layers will be overlapped if no more data
!        ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
!        OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
!        READ (atmos_unit, '(144i3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
!        CLOSE(atmos_unit) 
!
!        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ncep' // yrc // monc // dayc // '.dat'      
!        ! Determine if file exists or not
!        INQUIRE (FILE= ncep_fname, EXIST= file_exist)
!        IF (.NOT. file_exist) THEN
!          WRITE(www_lun, *) &
!               'Warning: no T profile file found, use monthly mean!!!'
!          ! already read the monthly mean above
!        ELSE
!          OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
!          ! NCEP misses the 775 level, which is shown in ECMWFT, other levels are the same
!          READ(atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i = 1, nlon), &
!               j = 1, nlat), k = 1, 3)
!          READ(atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i = 1, nlon), &
!               j = 1, nlat), k = 5, 18)
!!          glbecmwft(:, :, 4) = (glbecmwft(:, :, 3) + glbecmwft(:, :, 5)) / 2.0
!! Note glbecmwft is an integer
!          glbecmwft(:, :, 4) = (glbecmwft(:, :, 3) + glbecmwft(:, :, 5)) / 2
!        ENDIF
!      ENDIF
!
!      first = .FALSE.
!    ENDIF
!
!    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
!         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
!    ecmwft = 0.0
!    DO i = 1, nblon
!      DO j = 1, nblat
!        ecmwft = ecmwft + glbecmwft(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
!      ENDDO
!    ENDDO
!
!    RETURN
!  END SUBROUTINE get_ecmwft

  ! Obtain ECMWF temperature profile
  SUBROUTINE get_ncept(year, month, day, lon, lat, ncept)

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, PARAMETER                            :: nlecm=22
    INTEGER, INTENT(IN)                           :: month, year, day
    REAL (KIND=dp), INTENT(IN)                    :: lon, lat
    REAL (KIND=dp), DIMENSION(nlecm), INTENT(OUT) :: ncept

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER              :: nlat=72, nlon=144, nalt=22
    REAL (KIND=dp), PARAMETER       :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
    CHARACTER (LEN=2)               :: yrc, monc, dayc
    CHARACTER (LEN=max_pathlen)             :: ecmwft_fname, ncep_fname
    INTEGER                         :: i, j, k, nblat, nblon
    INTEGER, DIMENSION(2)           :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)    :: latfrac, lonfrac
    LOGICAL                         :: file_exist


    INTEGER, SAVE, DIMENSION(nlon, nlat, nalt) :: glbncept
    LOGICAL, SAVE                              :: first = .TRUE.

    IF (first) THEN
      WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
      WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
      WRITE(yrc,  '(I2.2)') MOD(year, 100) ! from 1997 to '97'

      ! Use NCEP for up to 10 mb and ECMWFT average for up between 10 and 1 mb
      ! ECMWFT average between 10mb and 1mb (7, 5, 3, 2, 1)', 
      ! NCEP: 17 layers ECMWFT: 23 layers (including 7, 5, 3, 2, 1, 775 mb)
      ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
      OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
      ! nalt + 1 = 23
      READ (atmos_unit, '(144I3)') (((glbncept(i, j, k), i=1, nlon), j=1, nlat), k=1, 1)
      READ (atmos_unit, '(144I3)') (((glbncept(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
      CLOSE(atmos_unit) 

      ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'nncept/ncep' // yrc // monc // dayc // '.dat'      
      ! Determine if file exists or not
      INQUIRE (FILE= ncep_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'Warning: no T profile file found, use monthly mean!!!'
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'nncept/ncepavg' // monc // '.dat'
      ENDIF
      OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
      READ(atmos_unit, '(144I3)') (((glbncept(i, j, k), i = 1, nlon), j = 1, nlat), k = 1, 17)

      CLOSE(atmos_unit) 
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    ncept = 0.0
    DO i = 1, nblon
      DO j = 1, nblat
        ncept = ncept + glbncept(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_ncept


  ! ===============================================================
  ! Obtain TOMS V8 ozone profiles (12 month, 18 latitude bands,
  !   3-10 profiles with total ozone at a step of 50 DU
  ! ===============================================================
  SUBROUTINE get_v8prof(month, day, lat, toz, which_clima, oz, ozref)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    INTEGER, PARAMETER                          :: nl = 11, nref = 60
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: month, day, which_clima
    REAL (KIND=dp),INTENT(IN)                    :: lat
    REAL (KIND=dp), INTENT(INOUT)                :: toz
    REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)   :: oz
    REAL (KIND=dp), DIMENSION(nref), INTENT(OUT) :: ozref

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER :: nlat=18, maxprof=10, nmon=12
    !CHARACTER (LEN=3), DIMENSION(12) :: months = (/'jan', 'feb','mar', &
    !'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=max_pathlen)                                :: ozprof_fname
    CHARACTER (LEN=200)                                :: line


    ! saved variables
    REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, maxprof, nl) :: ozprofs
    INTEGER,        SAVE, DIMENSION(nmon, nlat)              :: nprofs
    REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, nref)        :: ozrefs
    LOGICAL,        SAVE                                     :: first = .TRUE.

    REAL (KIND=dp)                                           :: frac, fdum, maxoz, minoz
    REAL (KIND=dp), DIMENSION(2)                             :: latfrac, monfrac
    INTEGER,        DIMENSION(2)                             :: latin, monin
    INTEGER :: i, j, ib, profin, nprof, nband, nm, im

    IF (first) THEN
      ! read the reference profile for climatology
      ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'mpclima/llmclima_prof.dat'
      OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')
      DO im = 1, nmon
        READ (atmos_unit, *)
        DO i = nref, 1, -1
          READ (atmos_unit, *) ozrefs(im, :, i)
        ENDDO
      ENDDO
      CLOSE (atmos_unit)

      ! read the TOMS V8 profiles
      IF (which_clima == 1) THEN 
        ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tomsv8_ozone_clima.dat'
        OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')

        ! Read until the target month        
        DO im = 1, nmon
          DO i = 1, nlat 
            READ(atmos_unit, *) 
            nprof = 1
            DO j = 1, maxprof
              READ (atmos_unit, '(A)') line;  READ (line, *) fdum

              IF (fdum < 999.0) THEN
                READ (line, *) fdum, ozprofs(im, i, nprof, :)
                nprof = nprof + 1
              ENDIF
            ENDDO
            nprofs(im, i) = nprof - 1              
          ENDDO
        ENDDO
        CLOSE (atmos_unit)
      ENDIF

      first = .FALSE.
    ENDIF

    IF (day <= 15) THEN
      monin(1) = month - 1
      IF (monin(1) == 0) monin(1) = 12
      monin(2) = month
      monfrac(1) = (15.0 - day) / 30.0
      monfrac(2) = 1.0 - monfrac(1)
    ELSE 
      monin(2) = month + 1
      IF (monin(2) == 13) monin(2) = 1
      monin(1) = month
      monfrac(2) = (day - 15) / 30.0
      monfrac(1) = 1.0 - monfrac(2)
    ENDIF
    nm = 2

    IF (lat <= -85.0) THEN
      nband = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= 85.0) THEN
      nband = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      nband = 2     ; frac = (lat + 85.0) / 10.0 + 1
      latin(1) = INT(frac); latin(2) = latin(1) + 1
      latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
    ENDIF

    ozref = 0.0
    DO im = 1, nm
      DO ib = 1, nband
        ozref = ozref + latfrac(ib) * monfrac(im) * ozrefs(monin(im), latin(ib), :) 
      ENDDO
    ENDDO

    oz = 0.0
    IF (which_clima == 1) THEN
      DO im = 1, nm
        DO ib = 1, nband   
          nprof = nprofs(monin(im), latin(ib))
          minoz = SUM(ozprofs(monin(im), latin(ib), 1, :))
          maxoz = SUM(ozprofs(monin(im), latin(ib), nprof, :))

          IF (toz < minoz) THEN
            WRITE(www_lun,*) 'Warning: no a priori profile available!!!'
            oz  = oz + ozprofs(monin(im), latin(ib), 1, :) * latfrac(ib) * toz / minoz * latfrac(ib)
          ELSE IF (toz > maxoz) THEN
            WRITE(www_lun,*) 'Warning: no a priori profile available!!!'
            oz = oz + ozprofs(monin(im), latin(ib), nprof, :) * latfrac(ib) * toz / minoz * latfrac(ib)
          ELSE
            profin = INT ((toz - minoz ) / 50.0) + 1
            IF (profin == 0) THEN 
              profin = 1
            ELSE IF (profin == nprof) THEN
              profin = profin - 1
            ENDIF

            frac = 1.0 - (toz - (minoz + (profin-1) * 50.0)) / 50.0
            oz = oz + latfrac(ib) * monfrac(im) * (frac * ozprofs(monin(im), latin(ib), profin, :) &
                 + (1.0 - frac) * ozprofs(monin(im), latin(ib), profin+1, :))
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    RETURN
  END SUBROUTINE get_v8prof

  ! ===============================================================================
  ! Obtain TOMS V8 temperatire profiles (11, levels, 12 months, 18 latitude bands)
  ! ===============================================================================
  SUBROUTINE get_v8temp(month, day, lat, v8temp)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    INTEGER, PARAMETER                              :: nl = 11
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                             :: month, day
    REAL (KIND=dp), INTENT(IN)                      :: lat
    REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: v8temp

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER                              :: nlat=18, nmon=12
    CHARACTER (LEN=max_pathlen)                             :: tfname

    ! saved variables
    REAL (KIND=dp), SAVE, DIMENSION(nl, nlat, nmon) :: tprofs
    LOGICAL,        SAVE                            :: first = .TRUE.

    REAL (KIND=dp)                                  :: frac
    REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
    INTEGER,        DIMENSION(2)                    :: latin, monin
    INTEGER                                         :: ib, nb, nm, im

    IF (first) THEN
      ! read the reference profile for climatology
      tfname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tv8_temp.dat'
      OPEN (UNIT = atmos_unit, file= tfname, status = 'unknown')
      READ  (atmos_unit, *) tprofs
      CLOSE (atmos_unit)     
      first = .FALSE.
    ENDIF

    IF (day <= 15) THEN
      monin(1) = month - 1
      IF (monin(1) == 0) monin(1) = 12
      monin(2) = month
      monfrac(1) = (15.0 - day) / 30.0
      monfrac(2) = 1.0 - monfrac(1)
    ELSE 
      monin(2) = month + 1
      IF (monin(2) == 13) monin(2) = 1
      monin(1) = month
      monfrac(2) = (day - 15) / 30.0
      monfrac(1) = 1.0 - monfrac(2)
    ENDIF
    nm = 2

    IF (lat <= -85.0) THEN
      nb = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= 85.0) THEN
      nb = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      nb = 2     ; frac = (lat + 85.0) / 10.0 + 1
      latin(1) = INT(frac); latin(2) = latin(1) + 1
      latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
    ENDIF

    v8temp = 0.0
    DO im = 1, nm
      DO ib = 1, nb
        v8temp = v8temp + latfrac(ib) * monfrac(im) * tprofs(:, latin(ib), monin(im)) 
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_v8temp


  ! Use MIPAS IG2 Temperature Profile cimatology 
  ! 121 levels (pressre altitude from 120 km to 0 km), 4 months (1,4,7,10)
  ! and 6 latitude bands (-75, -45, -10, 10, 45, 75)
  SUBROUTINE GET_MIPASIG2T(month, day, lat, xx, yy)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    INTEGER, PARAMETER                              :: nl = 121
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                             :: month, day
    REAL (KIND=dp), INTENT(IN)                      :: lat
    REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: xx, yy

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER                              :: nlat=6, nmon=4
    REAL (KIND=dp), DIMENSION(1:nmon), PARAMETER    :: mons = (/0.5, 3.5, 6.5, 9.5/)
    REAL (KIND=dp), DIMENSION(1:nlat), PARAMETER    :: lats = (/-75.0, -45.0, -10.0, 10.0, 45.0, 75.0/)

    CHARACTER (LEN=max_pathlen)                             :: fname

    ! saved variables
    REAL (KIND=dp), SAVE, DIMENSION(nl, nlat, nmon) :: profs
    REAL (KIND=dp), SAVE, DIMENSION(nl)             :: pres0
    LOGICAL,        SAVE                            :: first = .TRUE.

    REAL (KIND=dp), DIMENSION(0:nlat)               :: temp
    REAL (KIND=dp)                                  :: frac, fmon
    REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
    INTEGER,        DIMENSION(2)                    :: latin, monin
    INTEGER                                         :: ib, nb, nm, im, i, nheader

    IF (first) THEN
      fname = TRIM(ADJUSTL(atmdbdir)) // 'mipasprof/MIPAS_IG2_Tclima.dat'
      nheader = 8

      OPEN (UNIT = atmos_unit, file= fname, status = 'unknown')
      DO i = 1, nheader
        READ (atmos_unit, *) 
      ENDDO

      DO im = 1, nmon
        READ (atmos_unit, *)
        DO i = nl, 1, -1
          READ (atmos_unit, *) temp(0:nlat)
          profs(i, 1:nlat, im) = temp(1:nlat)
        ENDDO
        READ (atmos_unit, *)        
      ENDDO

      DO i = 1, nl
        pres0(i) = 1013.25 * 10. ** (- (i - 1.0) / 16. )
      ENDDO

      CLOSE (atmos_unit)     
      first = .FALSE.
    ENDIF

    fmon = month - 1.0 + 1.0 * day / 31.
    IF (fmon <= mons(1)) THEN
      monin(1) = nmon; monin(2) = 1
      frac = 1.0 - (fmon + 2.5) / 3.0
      monfrac(1) = frac; monfrac(2) = 1.0 - frac
    ELSE IF (fmon >= mons(nmon) ) THEN
      monin(1) = nmon; monin(2) = 1
      frac = 1.0 - (fmon - mons(nmon)) / 3.0
      monfrac(1) = frac; monfrac(2) = 1.0 - frac
    ELSE
      DO i = 2, nmon
        IF (fmon < mons(i)) EXIT
      ENDDO

      monin(1) = i - 1; monin(2) = i
      frac     = 1.0 - (fmon - mons(i-1)) / (mons(i) - mons(i-1))
      monfrac(1) = frac; monfrac(2) = 1.0 - frac
    ENDIF
    nm = 2

    IF (lat <= lats(1)) THEN
      nb = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= lats(nlat)) THEN
      nb = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      DO i = 2, nlat
        IF (lat < lats(i)) EXIT
      ENDDO

      nb = 2; latin(1) = i - 1; latin(2) = i
      frac     = 1.0 - (lat - lats(i-1)) / (lats(i) - lats(i-1))
      latfrac(1) = frac; latfrac(2) = 1.0 - frac
    ENDIF

    xx = pres0; yy = 0.0
    DO im = 1, nm
      DO ib = 1, nb
        yy = yy + latfrac(ib) * monfrac(im) * profs(:, latin(ib), monin(im)) 
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE GET_MIPASIG2T



  ! Use MIPAS IG2 Temperature Profile cimatology 
  ! 121 levels (pressre altitude from 120 km to 0 km), 4 months (1,4,7,10)
  ! and 6 latitude bands (-75, -45, -10, 10, 45, 75)
  SUBROUTINE GET_MIPASIG2O3(month, day, lat, xx, yy)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    INTEGER, PARAMETER                              :: nl = 121
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                             :: month, day
    REAL (KIND=dp), INTENT(IN)                      :: lat
    REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: xx, yy

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER                              :: nlat=6, nmon=4
    REAL (KIND=dp), DIMENSION(1:nmon), PARAMETER    :: mons = (/0.5, 3.5, 6.5, 9.5/)
    REAL (KIND=dp), DIMENSION(1:nlat), PARAMETER    :: lats = (/-75.0, -45.0, -10.0, 10.0, 45.0, 75.0/)

    CHARACTER (LEN=max_pathlen)                             :: fname

    ! saved variables
    REAL (KIND=dp), SAVE, DIMENSION(nl, nlat, nmon) :: profs
    REAL (KIND=dp), SAVE, DIMENSION(nl)             :: pres0
    LOGICAL,        SAVE                            :: first = .TRUE.

    REAL (KIND=dp), DIMENSION(0:nlat)               :: temp
    REAL (KIND=dp)                                  :: frac, fmon
    REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
    INTEGER,        DIMENSION(2)                    :: latin, monin
    INTEGER                                         :: ib, nb, nm, im, i, nheader

    IF (first) THEN
      fname = TRIM(ADJUSTL(atmdbdir)) // 'mipasprof/MIPAS_IG2_O3clima.dat'
      nheader = 9

      OPEN (UNIT = atmos_unit, file= fname, status = 'unknown')
      DO i = 1, nheader
        READ (atmos_unit, *) 
      ENDDO

      DO im = 1, nmon
        READ (atmos_unit, *)
        DO i = nl, 1, -1
          READ (atmos_unit, *) temp(0:nlat)
          profs(i, 1:nlat, im) = temp(1:nlat)
        ENDDO
        READ (atmos_unit, *)        
      ENDDO

      DO i = 1, nl
        pres0(i) = 1013.25 * 10. ** (- (i - 1.0) / 16. )
      ENDDO

      CLOSE (atmos_unit)     
      first = .FALSE.
    ENDIF

    fmon = month - 1.0 + 1.0 * day / 31.
    IF (fmon <= mons(1)) THEN
      monin(1) = nmon; monin(2) = 1
      frac = 1.0 - (fmon + 2.5) / 3.0
      monfrac(1) = frac; monfrac(2) = 1.0 - frac
    ELSE IF (fmon >= mons(nmon) ) THEN
      monin(1) = nmon; monin(2) = 1
      frac = 1.0 - (fmon - mons(nmon)) / 3.0
      monfrac(1) = frac; monfrac(2) = 1.0 - frac
    ELSE
      DO i = 2, nmon
        IF (fmon < mons(i)) EXIT
      ENDDO

      monin(1) = i - 1; monin(2) = i
      frac     = 1.0 - (fmon - mons(i-1)) / (mons(i) - mons(i-1))
      monfrac(1) = frac; monfrac(2) = 1.0 - frac
    ENDIF
    nm = 2

    IF (lat <= lats(1)) THEN
      nb = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= lats(nlat)) THEN
      nb = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      DO i = 2, nlat
        IF (lat < lats(i)) EXIT
      ENDDO

      nb = 2; latin(1) = i - 1; latin(2) = i
      frac     = 1.0 - (lat - lats(i-1)) / (lats(i) - lats(i-1))
      latfrac(1) = frac; latfrac(2) = 1.0 - frac
    ENDIF

    xx = pres0; yy = 0.0
    DO im = 1, nm
      DO ib = 1, nb
        yy = yy + latfrac(ib) * monfrac(im) * profs(:, latin(ib), monin(im)) 
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE GET_MIPASIG2O3

  SUBROUTINE get_surfalt(lon, lat, z0)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    REAL (KIND=dp), INTENT(IN)     :: lon, lat
    REAL (KIND=dp), INTENT(OUT)    :: z0

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER             :: nlat=360, nlon=720
    REAL (KIND=dp), PARAMETER      :: longrid = 0.5, latgrid = 0.5, lon0=-180.0, lat0=-90.0

    INTEGER                        :: i, j, nblat, nblon
    LOGICAL                        :: file_exist
    CHARACTER (LEN=max_pathlen)            :: surfalt_fname
    INTEGER, DIMENSION(2)          :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

    INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbz
    LOGICAL, SAVE                        :: first = .TRUE.

    IF (first) THEN
      surfalt_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/tomsv7_terrain.dat'

      ! Determine if file exists or not
      INQUIRE (FILE= surfalt_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        STOP 'No Terrain Elevation datafile found!!!'
      ENDIF

      OPEN (UNIT = atmos_unit, file = surfalt_fname, status = 'unknown')
      DO i = 1, 4
        READ(atmos_unit, *)
      ENDDO

      READ (atmos_unit, '(720I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    z0 = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO
    z0 = z0 / 1000.0  ! convert tp km

    RETURN
  END SUBROUTINE get_surfalt

  SUBROUTINE get_ncepreso_surfalt(lon, lat, z0)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    REAL (KIND=dp), INTENT(IN)     :: lon, lat
    REAL (KIND=dp), INTENT(OUT)    :: z0

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER             :: nlat=72, nlon=144
    REAL (KIND=dp), PARAMETER      :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0

    INTEGER                        :: i, j, nblat, nblon
    LOGICAL                        :: file_exist
    CHARACTER (LEN=max_pathlen)            :: surfalt_fname
    INTEGER, DIMENSION(2)          :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

    INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbz
    LOGICAL, SAVE                        :: first = .TRUE.

    IF (first) THEN
      surfalt_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/dem2.5x2.5.dat'

      ! Determine if file exists or not
      INQUIRE (FILE= surfalt_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        STOP 'No Terrain Elevation datafile found!!!'
      ENDIF

      OPEN (UNIT = atmos_unit, file = surfalt_fname, status = 'unknown')
      DO i = 1, 4
        READ(atmos_unit, *)
      ENDDO

      READ (atmos_unit, '(144I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    z0 = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO
    z0 = z0 / 1000.0  ! convert tp km

    RETURN
  END SUBROUTINE get_ncepreso_surfalt

  SUBROUTINE get_finereso_surfalt(lon, lat, z0)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    REAL (KIND=dp), INTENT(IN)     :: lon, lat
    REAL (KIND=dp), INTENT(OUT)    :: z0

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER             :: nlat=1800, nlon=3600
    REAL (KIND=dp), PARAMETER      :: longrid = 0.1, latgrid = 0.1, lon0=-180.0, lat0=-90.0

    INTEGER                        :: i, j, nblat, nblon
    LOGICAL                        :: file_exist
    CHARACTER (LEN=max_pathlen)            :: surfalt_fname
    INTEGER, DIMENSION(2)          :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

    INTEGER, SAVE, DIMENSION(nlon, nlat) :: glbz
    LOGICAL, SAVE                        :: first = .TRUE.

    IF (first) THEN
      surfalt_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/dem0.1x0.1.dat'

      ! Determine if file exists or not
      INQUIRE (FILE= surfalt_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        STOP 'No Terrain Elevation datafile found!!!'
      ENDIF

      OPEN (UNIT = atmos_unit, file = surfalt_fname, status = 'unknown')
      DO i = 1, 4
        READ(atmos_unit, *)
      ENDDO

      READ (atmos_unit, '(3600I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    z0 = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO
    z0 = z0 / 1000.0  ! convert tp km

    RETURN
  END SUBROUTINE get_finereso_surfalt


  SUBROUTINE get_geoschem_o31(month, lon, lat, ps, ozprof, nz, ntp)  
    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY: atmdbdir 
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)         :: month, nz, ntp
    REAL (KIND=dp), INTENT(IN)  :: lon, lat
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER               :: nlat=91, nlon=144, nalt=19
    REAL (KIND=dp), PARAMETER        :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
    INTEGER                          :: errstat, i, j, k, nblat, nblon, ntp0, nalt0

    REAL (KIND=dp), DIMENSION(nalt)  :: gprof
    REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz
    INTEGER, DIMENSION(2)            :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)     :: latfrac, lonfrac

    REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat, nalt) :: geosoz
    LOGICAL, SAVE                                     :: first = .TRUE.

    ! Correct coordinates
    REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
         .987871d0, .954730d0, .905120d0, .845000d0, .78d0, .710000d0,      &
         .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
         .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
         .102098d0/)

    REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
    CHARACTER (LEN=3), DIMENSION(12) :: months = (/'jan', 'feb','mar', 'apr', &
         'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)

    CHARACTER (LEN=max_pathlen)      :: geosfile

    IF (first) THEN
      geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(month) // '_o3_avg.dat'
      OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
      READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    gprof = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO
    geospres = pres * ps(0)

    ! Integrate from ppm to DU  
    cumoz = 0.0
    DO i = 1, nalt  ! 1000.0 / 1.25 / 1013.25 = 1.0 / 1.2665625
      cumoz(i) = cumoz(i-1) + gprof(i) * (geospres(i-1) - geospres(i)) / 1.266525
      IF (ANY(geosoz(lonin(1:nblon), latin(1:nblat), i) <= 0.0)) THEN
        j = i - 1; EXIT
      ELSE
        j = i
      ENDIF
    ENDDO
    nalt0 = j

    DO i = 1, ntp
      IF (ps(i) < geospres(nalt0)) THEN
        ntp0 = i - 1; EXIT
      ENDIF
    ENDDO

    CALL BSPLINE(geospres, cumoz, nalt0+1, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
    tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)     
    ozprof(1:ntp0) =  tempoz(1:ntp0) !* SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) *
    ! use profile shape only
    !ozprof(1:ntp) =  tempoz(1:ntp) * SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) 

    RETURN  
  END SUBROUTINE GET_GEOSCHEM_O31

  SUBROUTINE get_logan_clima(month, lon, lat, ps, ozprof, nz, ntp)  
    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY: atmdbdir 
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                             :: month, nz, ntp
    REAL (KIND=dp), INTENT(IN)                      :: lon, lat
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER        :: nlat=46, nlon=72, nalt=13
    REAL (KIND=dp), PARAMETER :: longrid = 5.0, latgrid = 4.0, lon0=-180.0, lat0=-92.0
    INTEGER                   :: errstat, i, j, nblat, nblon, ntp0!, k

    REAL (KIND=dp), DIMENSION(nalt)             :: gprof
    REAL (KIND=dp), DIMENSION(0:nz)             :: tempoz
    INTEGER, DIMENSION(2)                       :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)                :: latfrac, lonfrac

    REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat, nalt) :: geosoz
    LOGICAL, SAVE                                     :: first = .TRUE.

    ! Correct coordinates
    REAL (KIND=DP), DIMENSION(1:nalt), PARAMETER:: pres = (/1000., 900., &
         800., 700., 600., 500., 400., 300., 250., 200., 150., 125., 100./)

    REAL (KIND=dp), DIMENSION(1:nalt)           :: cumoz, presmod
!    CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb', &
!       'mar', 'apr',  'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=2)                           :: monc
    CHARACTER (LEN=max_pathlen)                         :: geosfile

    IF (first) THEN
      WRITE(monc, '(I2.2)') month
      geosfile = TRIM(ADJUSTL(atmdbdir)) // 'logan_clima/ozone.13.4x5.' // monc

      OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
      READ(atmos_unit, '(9E10.3)') geosoz
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    gprof = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        gprof = gprof + geosoz(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO

    ! Integrate from ppb to DU
    cumoz = 0.0
    DO i = 2, nalt  ! 2533.125 = 2 * 1.25 * 1013.25
      cumoz(i) = cumoz(i-1) + (gprof(i-1) + gprof(i)) * (pres(i-1) - pres(i)) / 2533.125 
    ENDDO

    presmod = pres
    DO i = 1, ntp
      IF (ps(i) < presmod(nalt)) THEN
        ntp0 = i - 1; EXIT
      ENDIF
    ENDDO
    IF (presmod(1) < ps(0))  presmod(1) = ps(0)

    CALL BSPLINE(presmod, cumoz, nalt, ps(0:ntp0), tempoz(0:ntp0), ntp0+1, errstat)
    tempoz(1:ntp0) = tempoz(1:ntp0) - tempoz(0:ntp0-1)    
    ozprof(1:ntp0) =  tempoz(1:ntp0)  ! use actual profile shape

    ! use profile shape only
    ! ozprof(1:ntp) =  tempoz(1:ntp) * SUM(ozprof(1:ntp)) / SUM(tempoz(1:ntp)) 

    RETURN  
  END SUBROUTINE GET_LOGAN_CLIMA

  ! ps: bottom up
  ! read GEOS-STRAT (V6.13) from May Fu
  SUBROUTINE GET_GEOSCHEM_HCHO(month, lon, lat, ps, hcho, nz)  
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: du2mol
    USE OMSAO_variables_module, ONLY: atmdbdir 
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: interpol

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: month, nz
    REAL (KIND=dp), INTENT(IN)                   :: lon, lat
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: hcho  ! in Dobson Units

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER        :: nlat=91, nlon=144, nalt=19
    REAL (KIND=dp), PARAMETER :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
    INTEGER                   :: errstat, i, j, k, nblat, nblon, ntp

    REAL (KIND=dp), DIMENSION(nalt)             :: gprof
    REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumhcho
    REAL (KIND=dp), DIMENSION(0:nz)             :: temphcho
    INTEGER, DIMENSION(2)                       :: latin, lonin
    REAL (KIND=dp), DIMENSION(2)                :: latfrac, lonfrac

    REAL (KIND=dp), SAVE, DIMENSION(nlon, nlat, nalt) :: geoshcho
    LOGICAL, SAVE                                     :: first = .TRUE.

    ! Correct coordinates
    REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
         .987871d0, .954730d0, .905120d0, .845000d0, .78d0, .710000d0,      &
         .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
         .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
         .102098d0/)
    CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr', &
         'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=max_pathlen)               :: geosfile

    IF (first) THEN
      geosfile = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_hcho/' // months(month) // '_hcho_avg.dat'
      OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
      DO i = 1, 5
        READ(atmos_unit, *) 
      ENDDO
      READ(atmos_unit, *) (((geoshcho(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
      CLOSE (atmos_unit)
      first = .FALSE.
    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    gprof = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        gprof = gprof + geoshcho(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO
    geospres = pres * ps(0)

    ! Integrate from ppb to DU  
    cumhcho = 0.0
    DO i = 1, nalt  ! 1266.5625 = 1.25 * 1013.25
      cumhcho(i) = cumhcho(i-1) + gprof(i) * (geospres(i-1) - geospres(i)) / 1266.5625
    ENDDO

    ! MAXLOC (maximum value), the index of maxloc starts from 1
    ntp = MINVAL(MINLOC(ps(0:nz), MASK = (ps(0:nz) >= geospres(nalt)))) - 1

    hcho = 0.0
    CALL INTERPOL(geospres, cumhcho, nalt+1, ps(0:ntp), temphcho(0:ntp), ntp+1, errstat)
    hcho(1:ntp) = temphcho(1:ntp) - temphcho(0:ntp-1)   ! DU at each layer 

    hcho(1:ntp) = hcho(1:ntp) * du2mol

    RETURN  
  END SUBROUTINE GET_GEOSCHEM_HCHO


  ! zs, ps: bottom up (BOS -> TOA)
  ! read stratospheric BrO from Chris + 0.2 ppbv in the troposphere (too small)
  SUBROUTINE GET_BRO(month, lat, zs, ps, sza, bro, nz)  

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: deg2rad, du2mol
    USE OMSAO_variables_module,  ONLY: atmdbdir 
    USE ozprof_data_module,      ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: month, nz
    REAL (KIND=dp),                  INTENT(IN)  :: lat
    REAL (KIND=dp), INTENT(IN)                   :: sza
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: zs, ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: bro     ! In Dobson Units

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER          :: nlat = 18, nalt = 25, maxsza = 17, nmax=601
    REAL (KIND=dp), PARAMETER   :: latgrid = 10.0
    INTEGER                     :: errstat, i, j, nsza, fidx, lidx, nband
    INTEGER, DIMENSION(2)       :: latin
    REAL(KIND=dp), DIMENSION(2) :: latfrac

    INTEGER,        SAVE, DIMENSION(nlat)               :: nszas
    REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt, maxsza) :: allbro
    REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt)         :: alts
    REAL (KIND=dp), SAVE, DIMENSION(nlat, maxsza)       :: szas
    REAL (KIND=dp), SAVE, DIMENSION(nmax)               :: ptmp0, brotmp, ptmp
    INTEGER, SAVE                                       :: ntmp
    LOGICAL, SAVE                                       :: first = .TRUE.

    REAL (KIND=dp), DIMENSION(maxsza)           :: tempszas
    REAL (KIND=dp), DIMENSION(nalt)             :: bprof
    REAL (KIND=dp), DIMENSION(nalt+5)           :: temprof, tempalt      ! Stuff to the surface
    REAL (KIND=dp), DIMENSION(0:nz)             :: tempbro, temp
    REAL (KIND=dp)                              :: frac, csza, airdens

    INTEGER :: which_bro = 0  ! 0: from PRATMO, 1: from SAO  2: from GEOS-5


    CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr', &
         'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=max_pathlen)               :: fname

    IF (first) THEN
      fname = TRIM(ADJUSTL(atmdbdir)) // 'BRO/' // months(month) // '_atm_BrO.dat' 
      OPEN(UNIT = atmos_unit, FILE = fname, status='old')

      DO i = 1, nlat 
        READ(atmos_unit, '(5X,2I5)') nszas(i), nszas(i)
        ! Only use AM data (from mid-night to noon)
        nszas(i) = nszas(i) / 2 
        READ(atmos_unit, '(5X, 34E13.4)') szas(i, 1:nszas(i)), szas(i, 1:nszas(i))
        DO j = nalt, 1, -1  ! Read top-down (j = 1, lowest level) 
          READ(atmos_unit, *) alts(i, j), allbro(i, j, 1:nszas(i)), &
               allbro(i, j, 1:nszas(i)), airdens
          allbro(i, j, 1:nszas(i)) = allbro(i, j, 1:nszas(i)) / airdens * 1.0E9 ! to ppbv           
        ENDDO
        READ(atmos_unit, *)
      ENDDO
      CLOSE (atmos_unit)

      IF (which_bro == 1) THEN
        fname = TRIM(ADJUSTL(atmdbdir)) // 'BRO/' // 'omsao_bro_oper.dat' 
        OPEN(UNIT = atmos_unit, FILE = fname, status='old')     
        READ(atmos_unit, *) ntmp
        DO i = 1, ntmp
          READ(atmos_unit, *) ptmp0(i), brotmp(i)
        ENDDO
        CLOSE(atmos_unit)
      ELSE IF (which_bro == 2) THEN
        fname = TRIM(ADJUSTL(atmdbdir)) // 'BRO/' // 'geos5_bro_barr_2008m0406.dat' 
        OPEN(UNIT = atmos_unit, FILE = fname, status='old')     
        READ(atmos_unit, *) ntmp
        DO i = 1, ntmp
          READ(atmos_unit, *) ptmp0(i), brotmp(i)
        ENDDO
        CLOSE(atmos_unit)
      ENDIF

      first = .FALSE.
    ENDIF

    ! Interpolate over latitude (avoid discontinuity) 
    IF (lat <= -85.0) THEN
      nband = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= 85.0) THEN
      nband = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      nband = 2; frac = (lat + 85.0) / latgrid + 1
      latin(1) = INT(frac); latin(2) = latin(1) + 1
      latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
    ENDIF

    tempalt(1:5) = (/0.0, 2.0, 4.0, 6.0, 8.0/)
    ! temprof(1:5) = 2.0E-4   ! Assume 0.2 pptv in the troposphere (from S.E. Chris)
    ! Now assume constant mixing ratio with lowest level mixing ratio
    ! This will better capture seasonal and latitudinal variation of the tropospheric BrO

    tempbro = 0.0
    csza = COS(sza * deg2rad)
    DO i = 1, nband
      nsza = nszas(latin(i))
      tempszas(1:nsza) = COS(szas(latin(i), 1:nsza) * deg2rad)

      fidx = MINVAL(MAXLOC( tempszas(1:nsza), MASK=(tempszas(1:nsza) <= csza)))
      IF (fidx == nsza) THEN
        bprof = allbro(latin(i), :, nsza)
      ELSE IF (fidx == 0) THEN
        bprof = allbro(latin(i), :, 1)
      ELSE
        lidx = fidx + 1
        frac = 1.0 - (csza - tempszas(fidx)) / (tempszas(lidx) - tempszas(fidx))
        bprof = frac * allbro(latin(i), :, fidx) + (1 - frac) * allbro(latin(i), :, lidx)
      ENDIF

      tempalt(6:nalt+5) = alts(latin(i), :);     temprof(6:nalt+5) = bprof
      temprof(1:5) = 0. !temprof(6)  ! Extratoplate to the whole troposphere
      IF (tempalt(nalt+5) < zs(nz)) tempalt(nalt+5) = zs(nz)
      IF (tempalt(1) > zs(0)) tempalt(1) = zs(0)  ! Avoid extrapolation

      ! Interpolate to GOME retrieval altitudes
      CALL BSPLINE(tempalt, temprof, nalt+5, zs(0:nz), temp(0:nz), nz+1, errstat)
      IF (temprof(5) == 0.0) THEN
        WHERE (zs(0:nz) < tempalt(6))
          temp(0:nz) = 0.0
        ENDWHERE
      ENDIF
      tempbro = tempbro + temp * latfrac(i)     
    ENDDO

    IF (which_bro /= 0) THEN
      ptmp(1:ntmp) = ptmp0(1:ntmp)
      IF (ptmp0(1) < ps(0)) ptmp(1) = ps(0)
      IF (ptmp0(ntmp) > ps(nz)) ptmp(ntmp) = ps(nz)  
      CALL BSPLINE(ptmp(1:ntmp), brotmp(1:ntmp), ntmp, ps(0:nz), tempbro(0:nz), nz+1, errstat)
    ENDIF

    !WRITE(90, *) nz + 1
    !WRITE(90, '(2D14.5)') ((ps(i), tempbro(i)), i=0, nz)
    !STOP

    ! Integrate from ppbv to DU  
    DO i = 1, nz
      ! accurate to within 1% (2533.125 = 2 * 1.25 * 1013.25
      bro(i) = (tempbro(i) + tempbro(i-1)) * (ps(i-1) - ps(i)) / 2533.125
    ENDDO
    bro(1:nz) = bro(1:nz) * du2mol

    RETURN  
  END SUBROUTINE GET_BRO


  ! zs, ps: bottom up (BOS -> TOA)
  ! read GEOS-STRAT (V6.13) from May Fu + stratospheric NO2 from Chris
  SUBROUTINE GET_NO2(month, lon, lat, zs, ps, sza, no2, nz)  

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: deg2rad, du2mol
    USE OMSAO_variables_module,  ONLY: atmdbdir 
    USE ozprof_data_module,      ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline, interpol

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: month, nz
    REAL (KIND=dp), INTENT(IN)                   :: lat, lon
    REAL (KIND=dp), INTENT(IN)                   :: sza
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: zs, ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: no2     ! In Dobson Units

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER          :: nlat = 18, nalt = 25, maxsza = 17
    INTEGER, PARAMETER          :: nglat=91, nglon=144, ngalt=19
    REAL (KIND=dp), PARAMETER   :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
    INTEGER                     :: errstat, i, j, k, nsza, fidx, lidx, nband, nblat, nblon, ntp
    INTEGER, DIMENSION(2)       :: latin,   lonin
    REAL(KIND=dp), DIMENSION(2) :: latfrac, lonfrac

    INTEGER,        SAVE, DIMENSION(nlat)                :: nszas
    REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt, maxsza)  :: allno2
    REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt)          :: alts
    REAL (KIND=dp), SAVE, DIMENSION(nlat, maxsza)        :: szas
    REAL (KIND=dp), SAVE, DIMENSION(nglon, nglat, ngalt) :: geosno2
    LOGICAL, SAVE                                        :: first = .TRUE.

    REAL (KIND=dp), DIMENSION(ngalt)            :: gprof
    REAL (KIND=dp), DIMENSION(0:ngalt)          :: geospres, geosalt
    REAL (KIND=dp), DIMENSION(maxsza)           :: tempszas
    REAL (KIND=dp), DIMENSION(nalt)             :: bprof, bpres
    REAL (KIND=dp), DIMENSION(0:nalt+ngalt)     :: cumno2, tempalt      
    REAL (KIND=dp), DIMENSION(0:nz)             :: temp
    REAL (KIND=dp)                              :: frac, csza, airdens

    ! Correct coordinates
    REAL (KIND=DP), DIMENSION(0:ngalt), PARAMETER:: pres = (/1.0d0,              &
         0.987871d0, 0.954730d0, 0.905120d0, 0.845000d0, 0.7800d0,   0.710000d0, &
         0.639000d0, 0.570000d0, 0.503000d0, 0.440000d0, 0.3800d0,   0.325000d0, &
         0.278000d0, 0.237954d0, 0.202593d0, 0.171495d0, 0.144267d0, 0.121347d0, &
         0.102098d0/)
    CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr',    &
         'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=max_pathlen)               :: fname

    IF (first) THEN
      fname = TRIM(ADJUSTL(atmdbdir)) // 'NO2/' // months(month) // '_strat_NO2.dat' 
      OPEN(UNIT = atmos_unit, FILE = fname, status='old')

      DO i = 1, nlat 
        READ(atmos_unit, '(5X,2I5)') nszas(i), nszas(i)
        ! Only use AM data (from mid-night to noon)
        nszas(i) = nszas(i) / 2 
        READ(atmos_unit, '(5X, 34E13.4)') szas(i, 1:nszas(i)), szas(i, 1:nszas(i))
        DO j = nalt, 1, -1  ! Read top-down (j = 1, lowest level) 
          READ(atmos_unit, *) alts(i, j), allno2(i, j, 1:nszas(i)), &
               allno2(i, j, 1:nszas(i)), airdens
          allno2(i, j, 1:nszas(i)) = allno2(i, j, 1:nszas(i)) / airdens * 1.0E9 ! to ppbv           
        ENDDO
        READ(atmos_unit, *)
      ENDDO
      CLOSE (atmos_unit)

      fname = TRIM(ADJUSTL(atmdbdir)) // 'NO2/' // months(month) // '_no2_avg.dat'
      OPEN(UNIT = atmos_unit, FILE = fname, status='old')
      DO i = 1, 5
        READ(atmos_unit, *) 
      ENDDO
      READ(atmos_unit, *) (((geosno2(i, j, k), k = 1, ngalt), j = 1, nglat), i = 1, nglon)
      CLOSE (atmos_unit)

      first = .FALSE.
    ENDIF

    !  Get GEOS-CHEM profiles (ppbv)  
    CALL get_gridfrac(nglon, nglat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
    gprof = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        gprof = gprof + geosno2(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO

    geospres = pres * ps(0)  
    CALL BSPLINE(ps, zs, nz+1, geospres,  geosalt, ngalt+1, errstat)

    ! Interpolate over latitude (avoid discontinuity) 
    IF (lat <= -85.0) THEN
      nband = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= 85.0) THEN
      nband = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      nband = 2; frac = (lat + 85.0) / 10.0 + 1
      latin(1) = INT(frac); latin(2) = latin(1) + 1
      latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
    ENDIF

    ! Get Profile in the stratosphere
    no2 = 0.0
    csza = COS(sza * deg2rad)
    DO i = 1, nband
      nsza = nszas(latin(i))
      tempszas(1:nsza) = COS(szas(latin(i), 1:nsza) * deg2rad)

      fidx = MINVAL(MAXLOC( tempszas(1:nsza), MASK=(tempszas(1:nsza) <= csza)))
      IF (fidx == nsza) THEN
        bprof = allno2(latin(i), :, nsza)
      ELSE IF (fidx == 0) THEN
        bprof = allno2(latin(i), :, 1)
      ELSE
        lidx = fidx + 1
        frac = 1.0 - (csza - tempszas(fidx)) / (tempszas(lidx) - tempszas(fidx))
        bprof = frac * allno2(latin(i), :, fidx) + (1 - frac) * allno2(latin(i), :, lidx)
      ENDIF

      ! Approximation: the layer-averaged profile is used as level (bottom) profile
      ntp = MINVAL(MAXLOC(geosalt(1:ngalt), MASK=(gprof > 0 .AND. geosalt(1:ngalt) < alts(latin(i), 1))))

      tempalt(0:ntp)          = geosalt(0:ntp)
      tempalt(ntp+1:ntp+nalt) = alts(latin(i), :)
      IF (tempalt(0) > zs(0)) tempalt(0) = zs(0)
      IF (tempalt(ntp+nalt) < zs(nz)) tempalt(ntp+nalt) = zs(nz)

      CALL BSPLINE(zs, ps, nz+1, alts(latin(i), :), bpres, nalt, errstat)
      cumno2 = 0.0
      DO j = 1, ntp
        cumno2(j) = cumno2(j-1) + gprof(i) * (geospres(j-1) - geospres(j)) / 1266.5625
      ENDDO
      cumno2(ntp + 1) = cumno2(ntp) + (gprof(ntp) + bprof(1)) * (geospres(ntp) - bpres(1)) / 2533.125
      DO j = 2, nalt
        k = ntp + j
        cumno2(k) = cumno2(k-1) + (bprof(j-1) + bprof(j)) * (bpres(j-1)-bpres(j)) / 2533.125
      ENDDO

      ! Interpolate to GOME retrieval altitudes
      CALL INTERPOL(tempalt(0:nalt+ntp), cumno2(0:nalt+ntp), nalt + ntp + 1, zs(0:nz), temp(0:nz), nz+1, errstat)
      temp(1:nz) = temp(1:nz) - temp(0:nz-1)
      no2 = no2 + temp(1:nz) * latfrac(i)     
    ENDDO

    no2(1:nz) = no2(1:nz) * du2mol

    RETURN  
  END SUBROUTINE GET_NO2

  ! ps: bottom up
  ! read GEOS-4 fields (96-97, 99-00) and GEOS-3 fields (98) of SO2 from Randall Martin and Neil Moore
  ! use average for other years
  ! 2006 is provide by Chulkyu Lee on 30 altitude grids (reduced GEOS-4) in HDF5 format
  ! will be used for after 2004
  SUBROUTINE GET_GEOSCHEM_SO2(year, month, lon, lat, ps, so2, ntp, nz)  
    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY: du2mol
    USE OMSAO_variables_module, ONLY: atmdbdir 
    USE ozprof_data_module,     ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: interpol

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: month, nz, year, ntp
    REAL (KIND=dp), INTENT(IN)                   :: lon, lat
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: so2  ! in Dobson Units

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER        :: nlat = 91, nlon = 144, nalt_pre2006 = 21, &
         nalt_post2006 = 30
    INTEGER, SAVE             :: nalt
    REAL (KIND=dp), PARAMETER :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, &
         lat0=-91.0
    INTEGER                   :: errstat, i, j, k, nblat, nblon!, error
    REAL (KIND=dp), DIMENSION(0:nalt_post2006)          :: geospres, cumso2
    REAL (KIND=dp), DIMENSION(nalt_post2006)            :: gprof
    REAL (KIND=dp), DIMENSION(0:nz)                     :: tempso2
    LOGICAL, SAVE                                       :: file_exist
    INTEGER, DIMENSION(2)                               :: latin,   lonin
    REAL(KIND=dp), DIMENSION(2)                         :: latfrac, lonfrac

    REAL (KIND=dp), SAVE, DIMENSION(nlon,nlat,nalt_post2006) :: geosso2
    LOGICAL, SAVE                                       :: first = .TRUE.

    ! Correct coordinates (for geos3 fields)
    REAL (KIND=DP), DIMENSION(0:nalt_pre2006), PARAMETER:: pres3 = (/ &
         1.00000D+00, 9.97095D-01, 9.91200D-01, 9.81500D-01, 9.67100D-01, 9.46800D-01, &
         9.19500D-01, 8.84000D-01, 8.39000D-01, 7.83000D-01, 7.18200D-01, 6.47600D-01, &
         5.74100D-01, 5.00000D-01, 4.27800D-01, 3.59500D-01, 2.97050D-01, 2.41950D-01, &
         1.94640D-01, 1.55000D-01, 1.22680D-01, 9.69000D-02/)

    ! for geos4 fields
    REAL (KIND=DP), DIMENSION(0:nalt_pre2006), PARAMETER:: ap4    = (/  &
         0.000000D0,   0.000000D0,   12.704939D0,  35.465965D0,  66.098427D0,  101.671654D0,  &
         138.744400D0, 173.403183D0, 198.737839D0, 215.417526D0, 223.884689D0, 224.362869D0, &
         216.864929D0, 201.192093D0, 176.929993D0, 150.393005D0, 127.837006D0, 108.663429D0, &
         92.365662D0,  78.512299D0,  56.387939D0,  40.175419D0/)
    REAL (KIND=DP), DIMENSION(0:nalt_pre2006), PARAMETER:: bp4    = (/  &
         1.000000D0,   0.985110D0,   0.943290D0,   0.867830D0, 0.764920D0,  0.642710D0,  &
         0.510460D0,   0.378440D0,   0.270330D0,   0.183300D0, 0.115030D0,  0.063720D0,  &
         0.028010D0,   0.006960D0,   0.000000D0,   0.000000D0, 0.000000D0,   0.000000D0, &
         0.000000D0,   0.000000D0,   0.000000D0,   0.000000D0/)
    !p(L) = AP4(L) + BP4(L) * ps

    ! for GEOS-4, 30 levels
!!$  REAL (KIND=DP), DIMENSION(nalt_post2006), PARAMETER:: etac    = (/  &
!!$       0.9926, 0.9707, 0.9300, 0.8680, 0.7891, 0.6987, 0.6031, 0.5135, 0.4373, 0.3724, &
!!$       0.3171, 0.2701, 0.2299, 0.1956, 0.1663, 0.1414, 0.1202, 0.1021, 0.0868, 0.0685, &
!!$       0.0491, 0.0348, 0.0245, 0.0148, 0.0068, 0.0029, 0.0011, 0.0004, 0.0001, 0.0000/)
    REAL (KIND=DP), DIMENSION(0:nalt_post2006), PARAMETER:: etae    = (/  &
         1.0000,  0.9851, 0.9562, 0.9039, 0.8321, 0.7460, 0.6515, 0.5547, 0.4723, 0.4022, &
         0.3425,  0.2917, 0.2484, 0.2114, 0.1798, 0.1528, 0.1299, 0.1104, 0.0939, 0.0798, &
         0.0573,  0.0408, 0.0288, 0.0201, 0.0094, 0.0041, 0.0017, 0.0006, 0.0002, 0.0001, 0.0/)



    CHARACTER (LEN=max_pathlen)               :: geosfile
    CHARACTER (LEN=2)                 :: yearc, monc


    IF (first) THEN
      WRITE(yearc, '(I2.2)') MOD(year, 100)
      WRITE(monc,  '(I2.2)') month  

      IF (year >= 2004) THEN
        ! Currently only have 2006 from Chulkyu, we will use for all years since
        geosfile = TRIM(ADJUSTL(atmdbdir)) // 'gcso2/so2_06' // monc // '.dat'
      ELSE
        geosfile = TRIM(ADJUSTL(atmdbdir)) // 'gcso2/so2_' // yearc // monc // '.dat'
      ENDIF

      ! Determine if file exists or not
      INQUIRE (FILE= geosfile, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'Warning: no SO2 profile file found, use 96-97, 99-00 average!!!'
        geosfile = TRIM(ADJUSTL(atmdbdir)) // 'gcso2/so2_avg' // monc // '.dat'
        nalt = nalt_pre2006
      ELSE 
        IF (year >= 2004) THEN 
          nalt = nalt_post2006
        ELSE
          nalt = nalt_pre2006
        ENDIF
      ENDIF

      OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
      DO i = 1, 5
        READ(atmos_unit, *) 
      ENDDO
      IF (year >= 2004 .AND. file_exist) THEN
        READ(atmos_unit, '(30E8.2)') (((geosso2(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
      ELSE
        READ(atmos_unit, '(21E8.2)') (((geosso2(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
      ENDIF
      CLOSE (atmos_unit)
      first = .FALSE.

    ENDIF

    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
         lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)


    gprof = 0.0
    DO i = 1, nblon
      DO j = 1, nblat 
        gprof = gprof + geosso2(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
      ENDDO
    ENDDO


    ! Need special processing for getting pressure profile for GEOS-4
    IF (year == 1998) THEN
      DO i = 0, nalt
        geospres(i) = pres3(i+1) * ps(0)
      ENDDO
    ELSEIF (year >= 2004 .AND. file_exist) THEN
      geospres = etae * ps(0)
    ELSE
      DO i = 0, nalt
        !geospres(i) = ap4(i+1) + bp4(i+1) * ps(0)
        geospres(i) = ap4(i) + bp4(i) * ps(0)
      ENDDO
    ENDIF


    ! Integrate from ppb to DU  
    cumso2 = 0.0
    DO i = 1, nalt  ! 1266.5625 = 1.25 * 1013.25
      cumso2(i) = cumso2(i-1) + gprof(i) * (geospres(i-1) - geospres(i)) / 1266.5625
    ENDDO
    !WRITE(90, *) nalt
    !WRITE(90, '(30D14.6)') (geospres(i), i=0, nalt)
    !WRITE(90, '(30D14.6)') (gprof(i), i=1, nalt)

    ! MAXLOC (maximum value), the index of maxloc starts from 1
    ! Note geos4 fields include some stratospheric part, which should not be used here
    !ntp = MINVAL(MINLOC(ps(0:nz), MASK = (ps(0:nz) >= geospres(nalt)))) - 1

    so2 = 0.0
    CALL INTERPOL(geospres, cumso2, nalt+1, ps(0:ntp), tempso2(0:ntp), ntp+1, errstat)
    so2(1:ntp) = tempso2(1:ntp) - tempso2(0:ntp-1)   ! DU at each layer 

    ! Assume 0.015 ppbv for stratospheric SO2
    ! 0.015 / 1.25 / 1013.25 = 1.1843E-5
    !DO i = ntp+1, nz
    !   so2(i) = (ps(i-1) - ps(i)) * 1.1843E-5
    !ENDDO


    so2(1:nz) = so2(1:nz) * du2mol


    RETURN  
  END SUBROUTINE GET_GEOSCHEM_SO2

  ! =====================================================================
  ! Obtain AURA MLS zonal mean ozone profiles and its standard deviations
  ! (quality flags applied) 0.1-215 mb (i.e., 10-64 km), 36 latitude bins
  ! =====================================================================
  SUBROUTINE get_mlso3prof(year, month, day, lat, nz, mnorstd, ps, zs, oz, ntp, errstat)
    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit!, which_clima
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline

    IMPLICIT NONE

    INTEGER, PARAMETER                           :: ml = 37, mlat=36
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: year, month, day, nz, mnorstd
    INTEGER, INTENT(OUT)                         :: errstat, ntp
    REAL (KIND=dp),INTENT(IN)                    :: lat
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
    REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz

    ! ======================
    ! Local variables
    ! ======================
    CHARACTER (LEN=max_pathlen)              :: mlsfname
    CHARACTER (LEN=2)                :: monc, dayc
    CHARACTER (LEN=4)                :: yrc
    LOGICAL                          :: file_exist
    INTEGER                          :: i, j, ib, nband, fidx, lidx, ntmpl, sl, el
    REAL (KIND=dp), DIMENSION (ml)   :: tmpoz, tmpozstd, ratio
    REAL (KIND=dp), DIMENSION (0:ml) :: cumoz
    REAL (KIND=dp), DIMENSION (0:nz) :: tmpcumoz, tmps
    REAL (KIND=dp), DIMENSION (2)    :: latfrac
    !REAL (KIND=dp)                   :: sumfrac
    INTEGER,        DIMENSION (2)    :: latin

    ! Saved variables
    INTEGER, SAVE                             :: nlat, nl
    REAL (KIND=dp), SAVE, DIMENSION(mlat, ml) :: mlsprofs, mlstds
    REAL (KIND=dp), SAVE, DIMENSION(0:ml)     :: mlsps
    REAL (KIND=dp), SAVE, DIMENSION(mlat)     :: mlslats
    LOGICAL,        SAVE                      :: first = .TRUE.

    errstat = 0
    IF (first) THEN
      WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
      WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
      WRITE(yrc,  '(I4.4)') year           ! from 1997 to '1997'

      ! Check the availablity of MLS ozone profiles
      mlsfname =TRIM(ADJUSTL(atmdbdir)) // 'MLSO3/zm_v02_' // yrc // monc // dayc // '.dat'

      ! Determine if file exists or not
      INQUIRE (FILE= mlsfname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'No MLS ozone profile found!!!'; errstat = -1; RETURN
      ENDIF

      OPEN (UNIT = atmos_unit, file = mlsfname, status = 'unknown')
      READ (atmos_unit, *) nl
      ! nl = nl-1      ! Only use above 147 mb
      IF (nl > ml) THEN
        WRITE(www_lun, *) 'Need to increase ml from ', ml, ' to ', nl
        errstat = -1; RETURN
      ENDIF
      ! Reading pressure bottom up
      READ (atmos_unit, *); READ (atmos_unit, *) (mlsps(i), i = nl, 0, -1)
      mlsps(0:nl) = LOG(mlsps(0:nl)) 

      READ (atmos_unit, *) nlat
      IF (nlat > mlat) THEN
        WRITE(www_lun, *) 'Need to increase mlat from ', mlat, ' to ', nlat
        errstat = -1; RETURN
      ENDIF
      READ (atmos_unit, *); READ (atmos_unit, *) mlslats(1:nlat)

      ! Reading bottom up
      READ (atmos_unit, *)
      DO i = 1, nlat
        READ (atmos_unit, *) (mlsprofs(i, j), j = nl, 1, -1)
      ENDDO
      READ (atmos_unit, *)
      DO i = 1, nlat
        READ (atmos_unit, *) (mlstds(i, j), j = nl, 1, -1)
      ENDDO
      CLOSE (atmos_unit)

      first = .FALSE.
    ENDIF

    IF (lat <= mlslats(1)) THEN
      nband = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= mlslats(nlat)) THEN
      nband = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      nband = 2 
      DO i = 2, nlat
        IF ( lat <= mlslats(i) ) THEN
          latin(1) = i - 1; latin(2) = i
          latfrac(2) = (lat - mlslats(i - 1)) / (mlslats(i) - mlslats(i-1))
          latfrac(1) = 1.0d0 - latfrac(2)
          EXIT
        ENDIF
      ENDDO
    ENDIF

    tmpoz(1:nl) = 0.0; tmpozstd(1:nl) = 0.0
    DO ib = 1, nband
      tmpoz(1:nl) = tmpoz(1:nl) + mlsprofs(latin(ib), 1:nl) * latfrac(ib)
      tmpozstd(1:nl) = tmpozstd(1:nl) + mlstds(latin(ib), 1:nl) * latfrac(ib)
    ENDDO
    ratio(1:nl) = tmpozstd(1:nl) / tmpoz(1:nl) * 100.0

    ! Only use MLS altitude range where reltative variability is < 50%
    ! Find first MLS layer to be used
    DO i = 1, nl
      IF (ratio(i) <= 50.0) THEN
        sl = i; EXIT
      ENDIF
    ENDDO

    ! Find last MLS layer to be used
    DO i = nl, 1, -1
      IF (ratio(i) <= 50.0) THEN
        el = i; EXIT
      ENDIF
    ENDDO

    IF (mnorstd == 2) tmpoz(1:nl) = tmpozstd(1:nz)

    ! Get cumulative ozone profile from (215 mb to 0.1 mb)
    cumoz(0) = 0.0
    DO i = 1, nl
      cumoz(i) = cumoz(i-1) + tmpoz(i)
    ENDDO
    tmps = LOG(ps(0:nz))

    fidx = MINVAL(MAXLOC(tmps(0:nz), MASK = (tmps(0:nz) <= mlsps(sl-1))) - 1)
    lidx = MINVAL(MINLOC(tmps(0:nz), MASK = (tmps(0:nz) >= mlsps(el))) - 1)
    ntmpl = lidx - fidx + 1
    !print *, sl, el, fidx, lidx, ntmpl
    !print *, EXP(mlsps(sl-1)), EXP(mlsps(el))
    !print *, EXP(tmps(fidx)), EXP(tmps(lidx))

    !print *, ' ozone before: ', SUM(oz)
    !print *, oz
    CALL BSPLINE(mlsps(0:nl), cumoz(0:nl), nl+1, tmps(fidx:lidx), &
         tmpcumoz(0:ntmpl-1), ntmpl, errstat)
    IF (errstat < 0) THEN
      WRITE(www_lun, *) 'GET_MLSO3PROF: BSPLINE error, errstat = ', errstat; RETURN
    ENDIF
    oz(fidx+1:lidx) = tmpcumoz(1:ntmpl-1) - tmpcumoz(0:ntmpl-2)
    !print *, ' ozone after: ', SUM(oz)
    !print *, fidx+1, lidx, SUM(oz(fidx+1:lidx))
    !print *, oz

    ntp = fidx

    RETURN
  END SUBROUTINE get_mlso3prof

  ! =====================================================================
  ! Obtain AURA MLS zonal mean ozone profiles and its standard deviations
  ! (quality flags applied) 0.1-215 mb (i.e., 10-64 km), 36 latitude bins
  ! =====================================================================
  SUBROUTINE get_mlso3prof_single(year, month, day, lat, nz, mnorstd, ps, &
       zs, oz, ntp, errstat)
    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir, tabdir
    USE ozprof_data_module,     ONLY: atmos_unit!, which_clima
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline

    IMPLICIT NONE

    INTEGER, PARAMETER                           :: ml = 37
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: year, month, day, nz, mnorstd
    INTEGER, INTENT(OUT)                         :: errstat, ntp
    REAL (KIND=dp),INTENT(IN)                    :: lat
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
    REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz

    ! ======================
    ! Local variables
    ! ======================
    CHARACTER (LEN=max_pathlen)              :: mlsfname
    CHARACTER (LEN=2)                :: monc, dayc
    CHARACTER (LEN=4)                :: yrc
    LOGICAL                          :: file_exist
    INTEGER                          :: i, j, fidx, lidx, ntmpl, sl, el, nl, theprof, ios, nm
    REAL (KIND=dp), DIMENSION (ml)   :: tmpoz, ratio, mlsprof, mlstd
    REAL (KIND=dp), DIMENSION (0:ml) :: cumoz, mlsps
    REAL (KIND=dp), DIMENSION (0:nz) :: tmpcumoz, tmps
    REAL (KIND=dp)                   :: tmplon, tmplat, tmpsza, tmptime

    errstat = 0

    WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
    WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
    WRITE(yrc,  '(I4.4)') year           ! from 1997 to '1997'

    ! Check the availablity of MLS ozone profiles
    mlsfname =TRIM(ADJUSTL(atmdbdir)) // 'MLSO3/mlso3_v02_' // yrc // monc // dayc // '.dat'

    ! Determine if file exists or not
    INQUIRE (FILE= mlsfname, EXIST= file_exist)
    IF (.NOT. file_exist) THEN
      WRITE(www_lun, *) 'No MLS ozone profile found!!!'; errstat = -1; RETURN
    ENDIF

    OPEN (UNIT = atmos_unit, file = TRIM(tabdir)//'INP/mlsprof_index.inp', status = 'unknown', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(www_lun, *) 'Do not know which profile to choose!!!'; errstat = -1; RETURN
    ELSE
      READ (atmos_unit, *) theprof; CLOSE (atmos_unit)
    ENDIF

    OPEN (UNIT = atmos_unit, file = mlsfname, status = 'unknown')
    READ (atmos_unit, *) nm, nl
    IF (nl > ml) THEN
      WRITE(www_lun, *) 'Need to increase ml from ', ml, ' to ', nl
      errstat = -1; CLOSE(atmos_unit); RETURN
    ENDIF
    IF (theprof > nm - 1) THEN
      WRITE(www_lun, *) 'Do not have this profile!!!'
      errstat = -1; CLOSE(atmos_unit); RETURN
    ENDIF

    ! Reading pressure bottom up
    READ (atmos_unit, *) (mlsps(i), i = nl, 0, -1)
    mlsps(0:nl) = LOG(mlsps(0:nl)) 

    ! Skip profiles until the one we want
    DO i = 1, theprof
      READ (atmos_unit, *); READ (atmos_unit, *); READ (atmos_unit, *)
    ENDDO

    READ (atmos_unit, *) i, tmplon, tmplat, tmpsza, tmptime
    !WRITE(www_lun, '(4F10.4)') tmplon, tmplat, tmpsza, tmptime

    ! Reading bottom up
    READ (atmos_unit, *) (mlsprof(j), j = nl, 1, -1)
    READ (atmos_unit, *) (mlstd(j),   j = nl, 1, -1)
    CLOSE (atmos_unit)
    ratio(1:nl) = mlstd(1:nl) / mlsprof(1:nl) * 100.0

    ! Only use MLS altitude range where reltative variability is < 50%
    ! Find first MLS layer to be used
    DO i = 1, nl
      IF (ratio(i) <= 50.0) THEN
        sl = i; EXIT
      ENDIF
    ENDDO

    ! Find last MLS layer to be used
    DO i = nl, 1, -1
      IF (ratio(i) <= 50.0) THEN
        el = i; EXIT
      ENDIF
    ENDDO

    IF (mnorstd == 1) THEN
      tmpoz(1:nl) = mlsprof(1:nz)
    ELSE IF (mnorstd == 2) THEN
      tmpoz(1:nl) = mlstd(1:nz)
    ENDIF

    ! Get cumulative ozone profile from 
    cumoz(0) = 0.0
    DO i = 1, nl
      cumoz(i) = cumoz(i-1) + tmpoz(i)
    ENDDO
    tmps = LOG(ps(0:nz))

    fidx = MINVAL(MAXLOC(tmps(0:nz), MASK = (tmps(0:nz) <= mlsps(sl-1))) - 1)
    lidx = MINVAL(MINLOC(tmps(0:nz), MASK = (tmps(0:nz) >= mlsps(el))) - 1)
    ntmpl = lidx - fidx + 1
    !print *, sl, el, fidx, lidx, ntmpl
    !print *, EXP(mlsps(sl-1)), EXP(mlsps(el))
    !print *, EXP(tmps(fidx)), EXP(tmps(lidx))

    !print *, ' ozone before: ', SUM(oz)
    !print *, oz
    CALL BSPLINE(mlsps(0:nl), cumoz(0:nl), nl+1, tmps(fidx:lidx), &
         tmpcumoz(0:ntmpl-1), ntmpl, errstat)
    IF (errstat < 0) THEN
      WRITE(www_lun, *) 'GET_MLSO3PROF: BSPLINE error, errstat = ', errstat; RETURN
    ENDIF
    oz(fidx+1:lidx) = tmpcumoz(1:ntmpl-1) - tmpcumoz(0:ntmpl-2)
    !print *, ' ozone after: ', SUM(oz)
    !print *, fidx+1, lidx, SUM(oz(fidx+1:lidx))
    !print *, oz

    ntp = fidx

    RETURN
  END SUBROUTINE get_mlso3prof_single



  SUBROUTINE GET_NORMTOZ(year, month, day, lat, toz, nz, ntp, ps, zs, oz, errstat)

    USE OMSAO_precision_module 
    USE OMSAO_variables_module, ONLY: atmdbdir
    USE ozprof_data_module,     ONLY: atmos_unit, norm_tropo3, which_toz
    USE OMSAO_errstat_module
    IMPLICIT NONE

    INTEGER, PARAMETER                           :: ntlat = 180

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: year, month, day, nz, ntp
    INTEGER, INTENT(OUT)                         :: errstat
    REAL (KIND=dp),INTENT(IN)                    :: lat
    REAL (KIND=dp),INTENT(INOUT)                 :: toz
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
    REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz

    ! ======================
    ! Local variables
    ! ======================
    CHARACTER (LEN=max_pathlen)              :: omto3fname
    CHARACTER (LEN=2)                :: monc, dayc
    CHARACTER (LEN=4)                :: yrc
    LOGICAL                          :: file_exist
    INTEGER                          :: i, ib, nband!, j
    REAL (KIND=dp), DIMENSION (2)    :: latfrac
    REAL (KIND=dp)                   :: sumfrac, mnalt, do3
    INTEGER,        DIMENSION (2)    :: latin

    ! Saved variables
    REAL (KIND=dp), SAVE, DIMENSION(ntlat) :: zmto3, tlats, zmalt
    LOGICAL,        SAVE                   :: first = .TRUE.

    errstat = 0
    IF (which_toz == 2) THEN 
      IF (first) THEN
        WRITE(monc, '(I2.2)') month          ! from 9 to '09' 
        WRITE(dayc, '(I2.2)') day            ! from 9 to '09'     
        WRITE(yrc,  '(I4.4)') year           ! from 1997 to '1997'

        ! Check the availablity of MLS ozone profiles
        omto3fname =TRIM(ADJUSTL(atmdbdir)) // 'OMTO3/zm_v003_' // yrc // 'm' // monc // dayc // '.dat'

        ! Determine if file exists or not
        INQUIRE (FILE= omto3fname, EXIST= file_exist)
        IF (.NOT. file_exist) THEN
          WRITE(www_lun, *) 'No Zonal Mean OMTO3 found!!!'; errstat = -1; RETURN
        ENDIF
        OPEN (UNIT = atmos_unit, file = omto3fname, status = 'unknown')
        DO i = 1, ntlat
          tlats(i) = REAL(i, KIND=dp) - 89.5
        ENDDO
        READ (atmos_unit, *)
!        READ (atmos_unit, *) ((zmto3(i), zmalt(i)), i = 1, ntlat)
        do i = 1, ntlat
          READ (atmos_unit, *) zmto3(i), zmalt(i)
        enddo
        CLOSE(atmos_unit)

        first = .FALSE.
      ENDIF

      IF (lat <= -89.5) THEN
        nband = 1; latin(1) = 1; latfrac(1) = 1.0
      ELSE IF (lat >= 89.5) THEN
        nband = 1; latin(1) = ntlat; latfrac(1) = 1.0
      ELSE
        nband = 2 
        DO i = 2, ntlat
          IF ( lat <= tlats(i)) THEN
            latin(1) = i - 1; latin(2) = i
            latfrac(2) = (lat - tlats(i - 1)) / (tlats(i) - tlats(i-1))
            latfrac(1) = 1.0d0 - latfrac(2)
            EXIT
          ENDIF
        ENDDO
      ENDIF

      toz = 0.0; sumfrac = 0.0; mnalt = 0.0
      DO ib = 1, nband
        IF ( zmto3(latin(ib)) > 0.0 ) THEN
          toz = toz + zmto3(latin(ib)) * latfrac(ib)
          mnalt = mnalt + zmalt(latin(ib)) * latfrac(ib)
          sumfrac = sumfrac + latfrac(ib)
        ENDIF
      ENDDO
      toz  = toz / sumfrac
      mnalt = mnalt / sumfrac / 1000.0 

      ! Accounting for different terrain height using approximate pressure conversion
      do3 = ( 1013.25 * (10.0**(-mnalt / 16.0)) - ps(0) ) / (ps(0) - ps(1)) * oz(1)
      toz = toz - do3
    ENDIF

    IF (norm_tropo3) THEN
      oz(1:ntp) = oz(1:ntp) * (toz - SUM(oz(ntp+1:nz))) / SUM(oz(1:ntp))
    ELSE
      oz(1:nz) = oz(1:nz) * toz / SUM(oz(1:nz))
    ENDIF

    RETURN
  END SUBROUTINE GET_NORMTOZ




  ! ===============================================================
  ! Obtain TOMS V8 ozone profiles (12 month, 18 latitude bands,
  !   3-10 profiles with total ozone at a step of 50 DU
  ! ===============================================================
  SUBROUTINE get_tomsv8_clima(month, day, lat, toz, nl, ps, apoz, oz, errstat)

    USE OMSAO_parameters_module, ONLY: p0
    USE OMSAO_precision_module 
    USE OMSAO_variables_module,  ONLY: atmdbdir
    USE ozprof_data_module,      ONLY: atmos_unit
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline, reverse

    IMPLICIT NONE

    INTEGER, PARAMETER                          :: nl0 = 11
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: month, day, nl
    INTEGER, INTENT(OUT)                         :: errstat
    REAL (KIND=dp),INTENT(IN)                    :: lat
    REAL (KIND=dp), INTENT(IN)                   :: toz
    REAL (KIND=dp), DIMENSION(0:nl), INTENT(IN)  :: ps
    REAL (KIND=dp), DIMENSION(nl), INTENT(IN)    :: apoz
    REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)   :: oz

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER :: nlat=18, maxprof=10, nmon=12
    !CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', &
    !'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=max_pathlen)                                :: ozprof_fname
    CHARACTER (LEN=200)                                :: line

    ! saved variables
    REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, maxprof, nl0) :: ozprofs
    INTEGER,        SAVE, DIMENSION(nmon, nlat)               :: nprofs
    REAL (KIND=dp), SAVE, DIMENSION(0:nl0)                    :: pv80
    LOGICAL,        SAVE                                      :: first = .TRUE.

    REAL (KIND=dp)                                            :: frac, fdum, maxoz, minoz
    REAL (KIND=dp), DIMENSION(nl0)                            :: oz0
    REAL (KIND=dp), DIMENSION(0:nl0)                          :: cum0
    REAL (KIND=dp), DIMENSION(0:nl)                           :: logps, cum
    REAL (KIND=dp), DIMENSION(2)                              :: latfrac, monfrac
    INTEGER,        DIMENSION(2)                              :: latin, monin
    INTEGER :: i, j, ib, profin, nprof, nband, nm, im

    CHARACTER (LEN=16), PARAMETER :: modulename = 'get_tomsv8_clima'

    IF (first) THEN
      ! read the TOMS V8 profiles
      ozprof_fname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tomsv8_ozone_clima.dat'
      OPEN (UNIT = atmos_unit, file= ozprof_fname, status = 'unknown')

      ! Read until the target month        
      DO im = 1, nmon
        DO i = 1, nlat 
          READ(atmos_unit, *) 
          nprof = 1
          DO j = 1, maxprof
            READ (atmos_unit, '(A)') line;  READ (line, *) fdum

            IF (fdum < 999.0) THEN
              READ (line, *) fdum, ozprofs(im, i, nprof, :)
              nprof = nprof + 1
            ENDIF
          ENDDO
          nprofs(im, i) = nprof - 1              
        ENDDO
      ENDDO
      CLOSE (atmos_unit)

      pv80(0) = ps(0)
      DO i = 1, nl0
        pv80(i) = p0 * 2.0D0 ** (+ i - nl0)
      ENDDO
      pv80(0:nl0) = LOG(pv80(0:nl0))

      first = .FALSE.
    ENDIF

    IF (day <= 15) THEN
      monin(1) = month - 1
      IF (monin(1) == 0) monin(1) = 12
      monin(2) = month
      monfrac(1) = (15.0 - day) / 30.0
      monfrac(2) = 1.0 - monfrac(1)
    ELSE 
      monin(2) = month + 1
      IF (monin(2) == 13) monin(2) = 1
      monin(1) = month
      monfrac(2) = (day - 15) / 30.0
      monfrac(1) = 1.0 - monfrac(2)
    ENDIF
    nm = 2

    IF (lat <= -85.0) THEN
      nband = 1; latin(1) = 1; latfrac(1) = 1.0
    ELSE IF (lat >= 85.0) THEN
      nband = 1; latin(1) = nlat; latfrac(1) = 1.0
    ELSE
      nband = 2     ; frac = (lat + 85.0) / 10.0 + 1
      latin(1) = INT(frac); latin(2) = latin(1) + 1
      latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
    ENDIF

    oz0 = 0.0
    DO im = 1, nm
      DO ib = 1, nband   
        nprof = nprofs(monin(im), latin(ib))
        minoz = SUM(ozprofs(monin(im), latin(ib), 1, :))
        maxoz = SUM(ozprofs(monin(im), latin(ib), nprof, :))

        IF (toz < minoz) THEN
          !WRITE(*,*), 'Warning: no a priori profile available!!!'
          oz0  = oz0 + ozprofs(monin(im), latin(ib), 1, :) * toz / minoz * latfrac(ib)
        ELSE IF (toz > maxoz) THEN
          !WRITE(*,*), 'Warning: no a priori profile available!!!'
          oz0 = oz0 + ozprofs(monin(im), latin(ib), nprof, :) * toz / maxoz * latfrac(ib)
        ELSE
          profin = INT ((toz - minoz ) / 50.0) + 1
          IF (profin == 0) THEN 
            profin = 1
          ELSE IF (profin == nprof) THEN
            profin = profin - 1
          ENDIF

          frac = 1.0 - (toz - (minoz + (profin-1) * 50.0)) / 50.0
          oz0 = oz0 + latfrac(ib) * monfrac(im) * (frac * ozprofs(monin(im), latin(ib), profin, :) &
               + (1.0 - frac) * ozprofs(monin(im), latin(ib), profin+1, :))
        ENDIF
      ENDDO
    ENDDO
    CALL REVERSE(oz0(1:nl0), nl0)

    ! Interpolate ozone profile to the input pressure grid
    cum0(0) = 0.0
    DO i = 1, nl0
      cum0(i) = cum0(i-1) + oz0(i)
    ENDDO
    logps = LOG(ps)
    pv80(0) = logps(0)

    errstat = pge_errstat_ok
    CALL BSPLINE(pv80, cum0, nl0+1, logps, cum, nl+1, errstat)
    IF (errstat < 0) THEN
      WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
      errstat = pge_errstat_error; RETURN
    ENDIF
    oz = cum(1:nl) - cum(0:nl-1)

    ! Correct for top few layers using original xap (based on McPeters Clima)
    DO i = 0, nl
      IF (logps(i) >= pv80(1)) EXIT
    ENDDO
    oz(1:i) = apoz(1:i) * SUM(oz(1:i)) / SUM(apoz(1:i))

    RETURN
  END SUBROUTINE get_tomsv8_clima


end module prepare_atmosphere
