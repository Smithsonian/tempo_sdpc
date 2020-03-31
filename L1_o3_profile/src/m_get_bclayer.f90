!
module m_get_bclayer

  public get_bc_layer
  private

  contains
  ! To get the spres for OMI pixel:
  ! 1. get spres at ncep/ncar reso
  ! 2. get z0 at ncep/ncar reso
  ! 3. get z0 at omi spatial resolution
  ! 4. get spres at omi spatial resolution
SUBROUTINE get_bc_layer (which_atm, nloc, the_lons, the_lats, &
                         ps0, pst, the_surfalt)
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: currpix, currline
  USE m_get_omigeos5, ONLY: geos5
  USE m_get_fnl, ONLY:get_fnl_spres,get_fnl_tpres, get_fnl_surfalt
  USE m_get_ncep, only: get_spres, get_ncepreso_surfalt, get_tpres
  USE m_get_met_tempo, ONLY: get_met_tempo, tempo=>thismet
  IMPLICIT NONE
  ! input/output variables
  INTEGER, INTENT(IN) :: which_atm, nloc
  REAL (KIND=dp), DIMENSION(nloc), INTENT(IN) :: the_lons, the_lats
  REAL (KIND=dp), INTENT(OUT) :: the_surfalt, ps0, pst
  ! local variables
  INTEGER :: i, errstat
  REAL (KIND=dp) :: ncepreso_z0, omi_z0
  REAL (KIND=dp), DIMENSION(nloc) :: fine_z0

  IF (which_atm == 0 ) THEN
     CALL get_spres(ps0)
     CALL get_ncepreso_surfalt(ncepreso_z0)
     CALL get_tpres(pst)
  ELSE IF (which_atm == 1) THEN
     CALL get_fnl_spres(ps0)
     CALL get_fnl_surfalt (ncepreso_z0)
     CALL get_fnl_tpres(pst)
  ELSE IF (which_atm == 2 ) THEN
     ps0 = geos5%spres(currpix, currline)
     ncepreso_z0 = geos5%phis(currpix, currline)
     pst = geos5%ptrp(currpix, currline)
  ELSE IF (which_atm == 3) THEN 
    ! defined in make_atm exept that the_surfalb is from TEMPO level 1b
    ! CALL get_met_tempo(errstat)
    ! ps0 = tempo%psurf
    ! ncepreso_z0 = tempo%z0
    ! pst = tempo%ptrop
  ENDIF
  IF (which_atm /= 2 .and. which_atm /=3) THEN 
    DO i = 1, nloc
     CALL get_finereso_surfalt(the_lons(i), the_lats(i), fine_z0(i))
    ENDDO
    omi_z0 = (SUM(fine_z0(1:4)) + fine_z0(5) * 4.) / 8.
    ! Adjust surface pressure
    ps0 = ps0 + 1013.25 * (10.**(-omi_z0/16.) - 10.**(-ncepreso_z0/16.))
    the_surfalt = omi_z0
  ELSE
    IF (which_atm /= 3) THEN 
      the_surfalt = ncepreso_z0
    ENDIF
  ENDIF
  RETURN
END SUBROUTINE get_bc_layer

SUBROUTINE get_finereso_surfalt(lon, lat, z0)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir, atmos_unit
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
  CHARACTER (LEN=130)            :: surfalt_fname
  INTEGER, DIMENSION(2)          :: latin, lonin
  REAL (KIND=dp), DIMENSION(2)   :: latfrac, lonfrac

  INTEGER, SAVE, DIMENSION(:,:), POINTER :: glbz
  LOGICAL, SAVE                          :: first = .TRUE.

  IF (first) THEN
      allocate (glbz(nlon, nlat))
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
END MODULE m_get_bclayer
