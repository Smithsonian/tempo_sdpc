MODULE m_get_fnl
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir, &
                              the_month, the_year, the_day, &
                              the_lon, the_lat
  USE OMSAO_parameters_module,     ONLY: atmos_unit

  ! used for interpolation
  INTEGER, PARAMETER :: nlfnl = 26
  INTEGER, PARAMETER, PRIVATE  :: nlat=180, nlon=360 
  REAL (KIND=dp), PARAMETER,PRIVATE  :: longrid = 1.0, latgrid = 1.0, lon0=-180.0, lat0=-90
  REAL (KIND=dp), PARAMETER,PRIVATE  :: lat_offset   = lat0 + latgrid / 2.0
  REAL (KIND=dp), PARAMETER,PRIVATE  :: lon_offset   = lon0 + longrid  / 2.0
  REAL (KIND=dp),PRIVATE                :: frac
  INTEGER,PRIVATE                       :: nblon, nblat
  INTEGER, DIMENSION(2),PRIVATE         :: latin, lonin
  REAL (KIND=dp), DIMENSION(2),PRIVATE  :: latfrac, lonfrac
  LOGICAL,PRIVATE                    :: file_exist
! others
  INTEGER, PRIVATE                   :: i, j, k
  CHARACTER (LEN=2),PRIVATE          :: monc, dayc
  CHARACTER (LEN=4),PRIVATE          :: yrc
  CHARACTER (LEN=130),PRIVATE        :: fnl_fname
  
  public get_fnl_temp, get_fnl_spres, & 
         get_fnl_sfct, get_fnl_tpres, get_fnl_surfalt

  private  
  CONTAINS

  SUBROUTINE get_fnl_temp(tprof)
  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
! FNL data modules
  REAL (KIND=DP), DIMENSION (nlfnl) :: tprof
  INTEGER (KIND=DP), DIMENSION (:,:,:), SAVE, POINTER :: glbtemp
  LOGICAL, SAVE :: first=.TRUE.

  IF (first) THEN
     allocate (glbtemp(nlon, nlat, nlfnl))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc,  '(I4.4)') the_year

     fnl_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltemp/fnltemp_' //yrc // monc // dayc // '.dat'    

     ! Determine if file exists or not
     INQUIRE (FILE= fnl_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no T profile file found, use monthly mean!!!'
          print * , fnl_fname
         fnl_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltemp/fnltempavg'// monc // '.dat'
     ENDIF
     
     ! NCEP FNL: 26 layers (top down from 10 to 1000 mb), but data will be
     ! bottom up after being read     
     OPEN (UNIT = atmos_unit, file = fnl_fname, status = 'unknown')
     READ(atmos_unit, '(360I3)') (((glbtemp(i, j, k), i = 1, nlon), j = 1,nlat), k = nlfnl, 1, -1)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF
 ! ================================================================
  nblat = 2; frac = (the_lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)

  IF (latin(1) == 0)   THEN
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF

  ! Circular in longitude direction
  nblon = 2; frac = (the_lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1

  tprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        tprof = tprof + glbtemp(lonin(i), latin(j), :) * lonfrac(i) *latfrac(j)
     ENDDO
  ENDDO
RETURN
END SUBROUTINE get_fnl_temp
 
SUBROUTINE get_fnl_spres (spres)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT(OUT)  :: spres
  INTEGER (KIND=DP), DIMENSION (:,:),SAVE, POINTER ::glbspres
  LOGICAL, SAVE                         :: first=.TRUE.
  IF (first) THEN
     allocate (glbspres(nlon, nlat))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc,  '(I4.4)') the_year
     fnl_fname =TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlsp/fnlsp_' // yrc// monc // dayc // '.dat'
     ! Determine if file exists or not
     INQUIRE (FILE= fnl_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no surface pressure file found, use monthlymean!!!'
       fnl_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlsp/fnlspavg' //monc // '.dat'
     ENDIF

     OPEN (UNIT = atmos_unit, file = fnl_fname, status = 'unknown')
     READ (atmos_unit, '(360I4)') ((glbspres(i, j), i=1, nlon), j=1, nlat)
     first = .FALSE.
  ENDIF

  nblat = 2; frac = (the_lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)

  IF (latin(1) == 0)   THEN
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF

  ! Circular in longitude direction
  nblon = 2; frac = (the_lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1

  spres = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        spres = spres + glbspres(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
RETURN
END SUBROUTINE get_fnl_spres

SUBROUTINE get_fnl_sfct (sfct)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT(OUT)  :: sfct
  INTEGER (KIND=DP), DIMENSION (:,:),SAVE, POINTER ::glbsfct
  LOGICAL, SAVE                :: first=.TRUE.

  IF (first) THEN
     allocate (glbsfct(nlon, nlat))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc,  '(I4.4)') the_year

     fnl_fname =TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlst/fnlst_' // yrc //monc // dayc // '.dat'

     ! Determine if file exists or not
     INQUIRE (FILE= fnl_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no surface temperature file found, use monthlymean!!!'
        fnl_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnlst/fnlstavg' //monc // '.dat'
     ENDIF

     OPEN (UNIT = atmos_unit, file = fnl_fname, status = 'unknown')
     READ (atmos_unit, '(360I3)') ((glbsfct(i, j), i=1, nlon), j=1, nlat)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  nblat = 2; frac = (the_lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)

  IF (latin(1) == 0)   THEN
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF

  ! Circular in longitude direction
  nblon = 2; frac = (the_lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1

  sfct = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        sfct = sfct + glbsfct(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
RETURN
RETURN
END SUBROUTINE get_fnl_sfct


SUBROUTINE get_fnl_tpres (tpres)
  IMPLICIT NONE
  REAL (KIND=dp), INTENT(OUT)  :: tpres
  INTEGER (KIND=dp), DIMENSION(:,:), SAVE, POINTER :: glbtpres
  LOGICAL, SAVE                         :: first=.TRUE.
  IF (first) THEN
    allocate (glbtpres(nlon, nlat))
    WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
    WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
    WRITE(yrc,  '(I4.4)') the_year

    fnl_fname =TRIM(ADJUSTL(atmdbdir)) // &
               'fnl13.75LST/fnltp/fnltp_' //yrc // monc // dayc // '.dat'
     ! Determine if file exists or not
     INQUIRE (FILE= fnl_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no tropopause pressure file found, use monthlymean!!!'
        fnl_fname = TRIM(ADJUSTL(atmdbdir)) // 'fnl13.75LST/fnltp/fnltpavg' //monc // '.dat'
     ENDIF

     OPEN (UNIT = atmos_unit, file = fnl_fname, status = 'unknown')
     READ (atmos_unit, '(360I3)') ((glbtpres(i, j), i=1, nlon), j=1, nlat)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  nblat = 2; frac = (the_lat - lat_offset) / latgrid + 1
  latin(1) = INT(frac); latin(2) = latin(1) + 1
  latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)

  IF (latin(1) == 0)   THEN
     latin(1) = 1;    latfrac(1) = 1.0; nblat = 1
  ENDIF
  IF (latin(2) > nlat) THEN
     latin(1) = nlat; latfrac(1) = 1.0; nblat = 1
  ENDIF

  ! Circular in longitude direction
  nblon = 2; frac = (the_lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1

  tpres = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        tpres = tpres + glbtpres(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
RETURN
END SUBROUTINE get_fnl_tpres

SUBROUTINE get_fnl_surfalt(z0)

  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(OUT)    :: z0
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, SAVE, DIMENSION(:,:), POINTER :: glbz
  LOGICAL, SAVE  :: first = .TRUE.

  IF (first) THEN
    allocate(glbz(nlon, nlat))
     fnl_fname = TRIM(ADJUSTL(atmdbdir)) // 'terrain_height/fnlsh1x1.dat'
     ! Determine if file exists or not
     INQUIRE (FILE= fnl_fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        STOP 'No Terrain Elevation datafile found!!!'
     ENDIF
     OPEN (UNIT = atmos_unit, file = fnl_fname, status = 'unknown')
     DO i = 1, 4
        READ(atmos_unit, *)
     ENDDO
     READ (atmos_unit, '(360I4)') ((glbz(i, j), i=1, nlon), j=1, nlat)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  ! Circular in longitude direction
  nblon = 2; frac = (the_lon - lon_offset) / longrid + 1
  lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
  lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  IF (lonin(1) == 0)   lonin(1) = nlon
  IF (lonin(2) > nlon) lonin(2) = 1
  z0 = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        z0 = z0 + glbz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  z0 = z0 / 1000.0  ! convert tp km

  RETURN
END SUBROUTINE get_fnl_surfalt
END MODULE m_get_fnl
