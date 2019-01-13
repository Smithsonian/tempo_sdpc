MODULE prepare_atmosphere
  ! ************************************************************************
  ! Author:  xiong liu
  ! Date  :  July 24, 2003
  ! Purpose: Routine to read atmospheric surface and tropopause pressure, 
  !          temperature, ozone, trace gases, surface altitdue and so on.
  ! ************************************************************************
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_year, the_day, the_lon, the_lat
  USE OMSAO_parameters_module, ONLY: atmos_unit
  USE OMSAO_errstat_module
CONTAINS
SUBROUTINE get_ecmwft( ecmwft) ! not used anywhere

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER          :: nlecm = 23
  REAL (KIND=dp), DIMENSION(nlecm), INTENT(OUT) :: ecmwft

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER           :: nlat=72, nlon=144, nalt=23
  REAL (KIND=dp), PARAMETER    :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
  CHARACTER (LEN=2)            :: yrc, monc, dayc
  CHARACTER (LEN=130)          :: ecmwft_fname, ncep_fname
  INTEGER                      :: i, j, k, nblat, nblon
  INTEGER, DIMENSION(2)        :: latin, lonin
  REAL (KIND=dp), DIMENSION(2) :: latfrac, lonfrac
  LOGICAL                      :: file_exist

  INTEGER, SAVE, DIMENSION(:,:,:), POINTER :: glbecmwft
  LOGICAL, SAVE  :: first = .TRUE.

  IF (first) THEN
    allocate (glbecmwft(nlon, nlat, nalt))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc,  '(I2.2)') MOD(the_year, 100) ! from 1997 to '97'

     ! use ECMWF
     IF (the_year <= 2001) THEN      
        ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwft' // yrc // monc // dayc // '.dat'      
        ! Determine if file exists or not
        INQUIRE (FILE= ecmwft_fname, EXIST= file_exist)
        IF (.NOT. file_exist) THEN
           WRITE(www_lun, *) 'Warning: no T profile file found, use monthly mean!!!'
           ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
        ENDIF
        OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
        READ (atmos_unit, '(144i3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
        CLOSE(atmos_unit)
     ELSE  ! Use NCEP for up to 10 mb and ECMWFT average for up between 10 and 1 mb
        ! ECMWFT average between 10mb and 1mb (7, 5, 3, 2, 1)', other layers will be overlapped if no more data
        ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
        OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
        READ (atmos_unit, '(144i3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nalt)
        CLOSE(atmos_unit) 
        
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ncep' // yrc // monc // dayc // '.dat'      
        ! Determine if file exists or not
        INQUIRE (FILE= ncep_fname, EXIST= file_exist)
        IF (.NOT. file_exist) THEN
           WRITE(www_lun, *) 'Warning: no T profile file found, use monthly mean!!!'
           ! already read the monthly mean above
        ELSE
           OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
           ! NCEP misses the 775 level, which is shown in ECMWFT, other levels are the same
           READ(atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i = 1, nlon), j = 1, nlat), k = 1, 3)
           READ(atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i = 1, nlon), j = 1, nlat), k = 5, 18)
           glbecmwft(:, :, 4) = (glbecmwft(:, :, 3) + glbecmwft(:, :, 5)) / 2.0
        ENDIF
     ENDIF

     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  ecmwft = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        ecmwft = ecmwft + glbecmwft(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE get_ecmwft

! Obtain monthly mean ECMWF temperature for 7, 5, 3, 2, 1 mb 
SUBROUTINE get_ecmwfavgt(ecmwft)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                            :: nl=5
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)    :: ecmwft

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER              :: nlat=72, nlon=144
  REAL (KIND=dp), PARAMETER       :: longrid = 2.5, latgrid = 2.5, lon0=-180.0,lat0=-90.0
  CHARACTER (LEN=2)               :: yrc, monc, dayc
  CHARACTER (LEN=130)             :: ecmwft_fname
  INTEGER                         :: il, i, j, k, nblat, nblon
  INTEGER, DIMENSION(2)           :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)    :: latfrac, lonfrac
  LOGICAL                         :: file_exist


  INTEGER, SAVE, DIMENSION(:,:,:), POINTER :: glbecmwft
  LOGICAL, SAVE                            :: first = .TRUE.

  IF (first) THEN
     allocate( glbecmwft(nlon, nlat, nl))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     

     ! There are 23 layers in the data, but only read the last 
     ! five layers from ECMWF, i.e., at 7, 5, 3, 2, 1 mb, respectively
     ecmwft_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc //'.dat'
     OPEN (UNIT = atmos_unit, file = ecmwft_fname, status = 'unknown')
     DO il = 1, 18
        READ (atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, 1)
     ENDDO
     READ (atmos_unit, '(144I3)') (((glbecmwft(i, j, k), i=1, nlon), j=1, nlat), k=1, nl)
     CLOSE(atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  ecmwft = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        ecmwft = ecmwft + glbecmwft(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE get_ecmwfavgt

! Obtain TOMS V8 temperatire profiles (11, levels, 12 months, 18 latitude bands)
! ===============================================================================
SUBROUTINE get_v8temp(v8temp)

  IMPLICIT NONE
  
  INTEGER, PARAMETER                              :: nl = 11
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: v8temp
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER                              :: nlat=18, nmon=12
  REAL (KIND=dp), PARAMETER                       :: latgrid=10, lat0=-90 
  CHARACTER (LEN=130)                             :: tfname

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), POINTER :: tprofs
  LOGICAL,        SAVE                            :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                    :: latin, monin
  INTEGER                                         :: ib, nb, nm, im

  IF (first) THEN
     allocate (tprofs(nl, nlat, nmon))
     ! read the reference profile for climatology
     tfname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tv8_temp.dat'
     OPEN (UNIT = atmos_unit, file= tfname, status = 'unknown')
     READ  (atmos_unit, *) tprofs
     CLOSE (atmos_unit)     
     first = .FALSE.
  ENDIF

  CALL get_monfrac(nmon, the_month, the_day, nm, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0, the_lat, nb, latfrac, latin)

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
SUBROUTINE GET_MIPASIG2T(xx, yy)
  
  IMPLICIT NONE
  
  INTEGER, PARAMETER                              :: nl = 121
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)      :: xx, yy
  
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER                              :: nlat=6, nmon=4
  REAL (KIND=dp), DIMENSION(1:nmon), PARAMETER    :: mons = (/0.5, 3.5, 6.5, 9.5/)
  REAL (KIND=dp), DIMENSION(1:nlat), PARAMETER    :: lats = (/-75.0, -45.0, -10.0, 10.0, 45.0, 75.0/)

  CHARACTER (LEN=130)                             :: fname

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), POINTER :: profs
  REAL (KIND=dp), SAVE, DIMENSION(nl)             :: pres0
  LOGICAL,        SAVE                            :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nlat)               :: temp
  REAL (KIND=dp)                                  :: frac, fmon
  REAL (KIND=dp), DIMENSION(2)                    :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                    :: latin, monin
  INTEGER                                         :: ib, nb, nm, im, i, nheader

  IF (first) THEN
     allocate (profs(nl, nlat, nmon))
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

  fmon = the_month - 1.0 + 1.0 * the_day / 31.
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

  IF (the_lat <= lats(1)) THEN
     nb = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (the_lat >= lats(nlat)) THEN
     nb = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     DO i = 2, nlat
        IF (the_lat < lats(i)) EXIT
     ENDDO

     nb = 2; latin(1) = i - 1; latin(2) = i
     frac     = 1.0 - (the_lat - lats(i-1)) / (lats(i) - lats(i-1))
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

SUBROUTINE get_monfrac(nmon, mon, day, nbmon, monfrac, monin)

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

SUBROUTINE get_latfrac( nlat, latgrid, lat0, lat,  nblat, latfrac, latin)

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

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlon, nlat
  REAL (KIND=dp), INTENT(IN)                :: lon0, lat0, lat, lon, longrid, latgrid
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

SUBROUTINE get_gridfrac1(nlon, nlat, nmon, longrid, latgrid, mongrid, lon0,lat0, mon0, &
       lon, lat, mon, nblon, nblat, nbmon, lonfrac, latfrac, monfrac, lonin,latin, monin)

    USE OMSAO_precision_module
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                       :: nlon, nlat, nmon
    REAL (KIND=dp), INTENT(IN)                :: lon0, lat0, mon0, lat, lon,mon, longrid, latgrid, mongrid
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
  
END MODULE prepare_atmosphere


