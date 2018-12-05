MODULE m_get_ncep
  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir, the_month, the_year, the_day, the_lon, the_lat
  USE OMSAO_parameters_module,     ONLY: atmos_unit
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! private variables
  INTEGER, PARAMETER         :: nlat=72, nlon=144,  nlecm=22
  REAL (KIND=dp), PARAMETER  :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
  INTEGER                    :: nblat, nblon, i,j, k
  INTEGER, DIMENSION(2)      :: latin, lonin
  REAL (KIND=dp), DIMENSION(2), PRIVATE   :: latfrac, lonfrac
  CHARACTER (LEN=2)          :: monc, yrc, dayc
  CHARACTER (LEN=130)        :: ncep_fname
  LOGICAL                    :: file_exist

  public  get_spres, get_ncepreso_surfalt, get_tpres, get_ncep_temp
  private get_gridfrac
  
  ! public variables
  CONTAINS
SUBROUTINE get_spres(spres)

  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(OUT)   :: spres
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, SAVE, DIMENSION(:,:), POINTER :: glbspres
  LOGICAL, SAVE :: first = .TRUE.

  IF (first) THEN
    allocate(glbspres(nlon, nlat))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc, '(I2.2)')  MOD(the_year, 100) ! from 1997 to '97'
     
     ncep_fname =TRIM(ADJUSTL(atmdbdir)) // 'nspres/spres' // yrc // monc // dayc // '.dat'
    
     ! Determine if file exists or not
     INQUIRE (FILE= ncep_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'Warning: no surface pressure file found, use monthly mean!!!'
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'nspres/spresavg' // monc // '.dat'
     ENDIF
     
     OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
     READ (atmos_unit, '(144I4)') ((glbspres(i, j), i=1, nlon), j=1, nlat)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
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
SUBROUTINE get_tpres(tpres)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(OUT)   :: tpres
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, SAVE, DIMENSION(:,:),POINTER :: glbtpres
  LOGICAL, SAVE :: first = .TRUE.

  IF (first) THEN
    allocate(glbtpres(nlon, nlat))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc, '(I2.2)')  MOD(the_year, 100) ! from 1997 to '97'
     
     ncep_fname =TRIM(ADJUSTL(atmdbdir)) // 'ntpres/tpres' // yrc // monc // dayc // '.dat'
    
     ! Determine if file exists or not
     INQUIRE (FILE= ncep_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(www_lun, *) 'Warning: no tropopause pressure file found, use monthly mean!!!'
        ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'ntpres/tpresavg' // monc // '.dat'
     ENDIF
     
     OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
     READ (atmos_unit, '(144I3)') ((glbtpres(i, j), i=1, nlon), j=1, nlat)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  tpres = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        tpres = tpres + glbtpres(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
      
  RETURN
END SUBROUTINE get_tpres


! Obtain ECMWF temperature profile
SUBROUTINE get_ncep_temp (ncept)
  
  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), DIMENSION(nlecm), INTENT(OUT) :: ncept

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, SAVE, DIMENSION(:,:,:),POINTER :: glbncept
  LOGICAL, SAVE  :: first = .TRUE.

  IF (first) THEN
     allocate(glbncept(nlon, nlat, nlecm))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc,  '(I2.2)') MOD(the_year, 100) ! from 1997 to '97'

     ! Use NCEP for up to 10 mb and ECMWFT average for up between 10 and 1 mb
     ! ECMWFT average between 10mb and 1mb (7, 5, 3, 2, 1)', 
     ! NCEP: 17 layers ECMWFT: 23 layers (including 7, 5, 3, 2, 1, 775 mb)
     ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'ecmwft/ecmwftavg' // monc // '.dat'
     OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
      
     ! nalt + 1 = 23
     READ (atmos_unit, '(144I3)') (((glbncept(i, j, k), i=1, nlon), j=1, nlat), k=1, 1)
     READ (atmos_unit, '(144I3)') (((glbncept(i, j, k), i=1, nlon), j=1, nlat), k=1, nlecm)
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
     first = .FALSE.
  ENDIF
 
  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  ncept = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        ncept = ncept + glbncept(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE get_ncep_temp

SUBROUTINE get_ncepreso_surfalt(z0)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(OUT)    :: z0
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, SAVE, DIMENSION(:,:), POINTER :: glbz
  LOGICAL, SAVE :: first = .TRUE.

  IF (first) THEN
    allocate(glbz(nlon, nlat))
      ncep_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/dem2.5x2.5.dat'
      
      ! Determine if file exists or not
      INQUIRE (FILE= ncep_fname, EXIST= file_exist)
      IF (.NOT. file_exist) THEN
         STOP 'No Terrain Elevation datafile found!!!'
      ENDIF
      
      OPEN (UNIT = atmos_unit, file = ncep_fname, status = 'unknown')
      DO i = 1, 4
         READ(atmos_unit, *)
      ENDDO
      
      READ (atmos_unit, '(144I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
      CLOSE (atmos_unit)
      first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  z0 = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  z0 = z0 / 1000.0  ! convert tp km
  
  RETURN
END SUBROUTINE get_ncepreso_surfalt

SUBROUTINE get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
  lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  USE OMSAO_precision_module
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

END MODULE m_get_ncep
