MODULE m_get_ncep

  USE OMSAO_precision_module 
  USE OMSAO_variables_module,  ONLY: atmdbdir, atmos_unit, &
                                the_month, the_year, the_day, the_lon, the_lat
  USE OMSAO_errstat_module
  USE m_utilities, ONLY: get_gridfrac
  IMPLICIT NONE

  ! dimension
  INTEGER, PARAMETER, PRIVATE :: nlat=72, nlon=144,  nlecm=22
  ! variables used for interpolation
  REAL (KIND=dp), PARAMETER, PRIVATE     :: longrid = 2.5, latgrid = 2.5, lon0=-180.0, lat0=-90.0
  INTEGER ,PRIVATE                       :: nblat, nblon
  INTEGER, DIMENSION(2),PRIVATE          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2), PRIVATE  :: latfrac, lonfrac
  ! others
  LOGICAL                     :: file_exist
  INTEGER                     :: i, j, k
  CHARACTER (LEN=2), PRIVATE  :: monc, yrc, dayc
  CHARACTER (LEN=130),PRIVATE :: ncep_fname

  PUBLIC  get_spres, get_ncepreso_surfalt, get_tpres, get_ncep_temp
  PRIVATE 
  
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

END MODULE m_get_ncep
