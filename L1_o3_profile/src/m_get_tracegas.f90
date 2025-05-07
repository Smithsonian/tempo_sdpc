!
! Re-structered by jbak
MODULE m_get_tracegas

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: deg2rad, du2mol
  USE OMSAO_variables_module, ONLY: atmdbdir, atmos_unit
  USE OMSAO_errstat_module
  USE m_ezspline_interpolation, ONLY: interpol, bspline
  USE m_utilities, ONLY: get_gridfrac, get_latfrac
  public get_geoschem_hcho, & ! unit is converted from ppb to DU
         get_geoschem_h2o, &
         get_geoschem_so2, &
         get_bro, &
         get_no2, &
         get_afglus_h2o
  private

  contains  
  ! ps: bottom up
  ! read GEOS-STRAT (V6.13) from May Fu
  SUBROUTINE GET_GEOSCHEM_HCHO(month, lon, lat, ps, hcho, nz)  
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

    REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: geoshcho
    LOGICAL, SAVE  :: first = .TRUE.

    ! Correct coordinates
    REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0, .710000d0,      &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)
    CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr', &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
    CHARACTER (LEN=130)               :: geosfile
  
    IF (first) THEN
       allocate(geoshcho(nlon, nlat, nalt))
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
    REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: allbro
    REAL (KIND=dp), SAVE, DIMENSION(:,:),   ALLOCATABLE :: alts
    REAL (KIND=dp), SAVE, DIMENSION(:,:),   ALLOCATABLE :: szas
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
    CHARACTER (LEN=130)               :: fname
  
    IF (first) THEN
       allocate (allbro(nlat, nalt, maxsza), alts(nlat, nalt), szas(nlat, maxsza))

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
  ! CRN: After 2006, use NO2 from Lok Lamsal produced for SCIA overpass times in 2006 (Originally provided
  ! in HDF on variable pressure grid, changed to ASCII format and interpolated in log(p) to constant grid 
  ! on GEOS-4 reduced levels
  SUBROUTINE GET_NO2(year, month, lon, lat, zs, ps, sza, no2, nz)  

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: year, month, nz
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
    REAL (KIND=dp), SAVE, DIMENSION(:,:,:),ALLOCATABLE  :: allno2
    REAL (KIND=dp), SAVE, DIMENSION(nlat, nalt)          :: alts
    REAL (KIND=dp), SAVE, DIMENSION(nlat, maxsza)        :: szas
    REAL (KIND=dp), SAVE, DIMENSION(:,:,:),ALLOCATABLE :: geosno2
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
    CHARACTER (LEN=130)               :: fname
    CHARACTER (LEN=2)                 :: monc
  
    IF (first) THEN
       allocate(geosno2(nglon, nglat, ngalt))
       allocate(allno2(nlat, nalt, maxsza))
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

       IF (year >= 2004) THEN
          WRITE(monc, '(I2.2)') month
          fname = TRIM(ADJUSTL(atmdbdir)) // 'NO2/gcno2_06' //monc// '.dat'
       ELSE
          fname = TRIM(ADJUSTL(atmdbdir)) // 'NO2/' // months(month) // '_no2_avg.dat'
       ENDIF

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

       IF (ntp == 0) ntp =1 
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
    INTEGER, PARAMETER        :: nlat = 91, nlon = 144, nalt_pre2006 = 21, nalt_post2006 = 30
    INTEGER, SAVE             :: nalt
    REAL (KIND=dp), PARAMETER :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
    INTEGER                   :: errstat, i, j, k, nblat, nblon
    REAL (KIND=dp), DIMENSION(0:nalt_post2006)          :: geospres, cumso2
    REAL (KIND=dp), DIMENSION(nalt_post2006)            :: gprof
    REAL (KIND=dp), DIMENSION(0:nz)                     :: tempso2
    LOGICAL, sAVE                                       :: file_exist
    INTEGER, DIMENSION(2)                               :: latin,   lonin
    REAL(KIND=dp), DIMENSION(2)                         :: latfrac, lonfrac

    REAL (KIND=dp), SAVE, DIMENSION(:,:,:),ALLOCATABLE :: geosso2
    LOGICAL, SAVE   :: first = .TRUE.

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

    CHARACTER (LEN=130)               :: geosfile
    CHARACTER (LEN=2)                 :: yearc, monc

    IF (first) THEN
     allocate(geosso2(nlon,nlat,nalt_post2006))
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

  ! xliu, 11/04/2011
  ! Use a single US standard atmospheric water vapor profile
  SUBROUTINE GET_AFGLUS_H2O (ps, h2o, nz)
    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                          :: nz
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps
    REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: h2o
    ! Local variables
    ! AFGL US standard atmosphere
    INTEGER, PARAMETER              :: nz0 = 51
    REAL (KIND=dp), DIMENSION(0:nz0), PARAMETER :: ps0 = (/ &
       1.0130E+03,1.0000E+03,9.8500E+02,9.7000E+02,9.5500E+02,9.4000E+02,9.2000E+02,9.0000E+02,8.7500E+02,8.5000E+02, &
       8.2500E+02,8.0000E+02,7.7000E+02,7.4000E+02,7.1000E+02,6.8000E+02,6.4000E+02,6.0000E+02,5.5000E+02,5.0000E+02, &
       4.5000E+02,4.1110E+02,3.5650E+02,3.0800E+02,2.6500E+02,2.2700E+02,1.9400E+02,1.6580E+02,1.4170E+02,1.0350E+02, &
       7.5650E+01,5.5290E+01,4.0470E+01,2.9720E+01,2.1867E+01,1.6186E+01,1.1970E+01,8.6499E+00,6.5195E+00,5.0669E+00, &
       3.8613E+00,2.8710E+00,2.1972E+00,1.6958E+00,1.3141E+00,1.0242E+00,7.9780E-01,6.2130E-01,4.8300E-01,3.7330E-01, &
       2.8670E-01,8.5000E-02/)
    REAL (KIND=dp), DIMENSION(0:nz0), PARAMETER :: cumh2o0 = (/ &  ! cumulative profile in molecules cm^-2
       0.0000000E+00, 2.1127570E+21, 4.4882070E+21, 6.7953580E+21, 9.0331740E+21, 1.1201484E+22, &
       1.3983768E+22, 1.6641685E+22, 1.9791958E+22, 2.2755888E+22, 2.5536794E+22, 2.8133193E+22, &
       3.0992833E+22, 3.3556958E+22, 3.5822291E+22, 3.7808143E+22, 4.0076334E+22, 4.1956598E+22, &
       4.3810544E+22, 4.5194388E+22, 4.6207524E+22, 4.6776434E+22, 4.7328851E+22, 4.7603807E+22, &
       4.7709864E+22, 4.7753148E+22, 4.7772574E+22, 4.7781574E+22, 4.7785887E+22, 4.7789937E+22, &
       4.7792267E+22, 4.7793955E+22, 4.7795220E+22, 4.7796186E+22, 4.7796926E+22, 4.7797486E+22, &
       4.7797906E+22, 4.7798216E+22, 4.7798445E+22, 4.7798616E+22, 4.7798748E+22, 4.7798845E+22, &
       4.7798918E+22, 4.7798974E+22, 4.7799017E+22, 4.7799050E+22, 4.7799075E+22, 4.7799095E+22, &
       4.7799111E+22, 4.7799123E+22, 4.7799133E+22, 4.7799154E+22/)
    INTEGER                         :: errstat
    REAL (KIND=dp), DIMENSION(0:nz) :: cumh2o

    h2o = 0.0
    CALL INTERPOL(ps0, cumh2o0, nz0+1, ps(0:nz), cumh2o(0:nz), nz+1, errstat)
    h2o(1:nz) = (cumh2o(1:nz) - cumh2o(0:nz-1))*3.0     ! 3 times
    RETURN
  END SUBROUTINE GET_AFGLUS_H2O

 !======================================================================
 ! xliu, 1/2/2015
 ! GEOS-Chem H2O (used in Helen Wang's OMI H2O product
 ! H2O actually comes from GEOS MERRA product (1x1 downgraded to 2.5x2)
 ! Data are generated from full chemistry run with 54 tracers
 !    with daily GFED3 biomass burning emissions
 ! GEOS-CHEM SIMULATION v9-01-03: GEOS-5 NOx-Ox-HC-Aerosol simulation
 !======================================================================

  SUBROUTINE GET_GEOSCHEM_H2O(month, lon, lat, ps, h2o, ntp, nz)  
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY: du2mol
  USE OMSAO_variables_module, ONLY: atmdbdir , atmos_unit

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: month, nz, ntp
  REAL (KIND=dp), INTENT(IN)                   :: lon, lat
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(OUT) :: h2o  ! in Dobson Units

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER        :: nlat = 91, nlon = 144, nalt = 35
  REAL (KIND=dp), PARAMETER :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
  INTEGER                   :: errstat, i, j, k, nblat, nblon
  REAL (KIND=dp), DIMENSION(0:nalt)          :: geospres, cumh2o
  REAL (KIND=dp), DIMENSION(nalt)            :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)            :: temph2o
  INTEGER, DIMENSION(2)                      :: latin,   lonin
  REAL(KIND=dp), DIMENSION(2)                :: latfrac, lonfrac

  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: geosh2o
  LOGICAL, SAVE                                       :: first = .TRUE.

  ! for MERRA fields (only first 35 layers up to 92 mb)
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: ap4    = (/  &
       0.000000,    0.048048,    6.593752,   13.134800,  19.613110,  26.092010, &
       32.570808,   38.982010,   45.339008,   51.696110,  58.053211,  64.362640, &
       70.621979,   78.834221,   89.099922,   99.365211, 109.181702, 118.958603, &
       128.695908,  142.910004,  156.259995,  169.608994, 181.619003, 193.097000, &
       203.259003,  212.149994,  218.776001,  223.897995, 224.363007, 216.865005, &
       201.192001,  176.929993,  150.393005,  127.836998, 108.663002,  92.365723/)
  
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: bp4    = (/  &
       1.000000,   0.984952,   0.963406,   0.941865,  0.920387,   0.898908, &
       0.877429,   0.856018,   0.834661,   0.813304,  0.791947,   0.770638, &
       0.749378,   0.721166,   0.685900,   0.650635,  0.615818,   0.581042, &
       0.546304,   0.494590,   0.443740,   0.392891,  0.343381,   0.294403, &
       0.246741,   0.200350,   0.156224,   0.113602,  0.063720,   0.028010, &
       0.006960,   0.000000,   0.000000,   0.000000,  0.000000,   0.000000/)  
  !p(L) = AP4(L) + BP4(L) * ps  
  CHARACTER (LEN=130)               :: geosfile
  CHARACTER (LEN=2)                 :: monc

  IF (first) THEN
     allocate (geosh2o(nlon, nlat, nalt))
     WRITE(monc,  '(I2.2)') month  
     geosfile = TRIM(ADJUSTL(atmdbdir)) // 'gch2o/gc_merra_h2o_2007' // monc // '.dat'

     OPEN(UNIT = atmos_unit, FILE = geosfile, status='old')
     DO i = 1, 5
        READ(atmos_unit, *) 
     ENDDO
     READ(atmos_unit, '(35E12.4)') (((geosh2o(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  gprof = 0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        gprof = gprof + geosh2o(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO

  geospres = ap4 + bp4 * ps(0)

  ! Integrate from ppb to DU  
  cumh2o = 0.0
  DO i = 1, nalt  ! 1266.5625 = 1.25 * 1013.25
     cumh2o(i) = cumh2o(i-1) + gprof(i) * (geospres(i-1) - geospres(i)) / 1266.5625
  ENDDO
  !WRITE(90, *) nalt
  !WRITE(90, '(30D14.6)') (geospres(i), i=0, nalt)
  !WRITE(90, '(30D14.6)') (gprof(i), i=1, nalt)

  h2o = 0.0
  CALL INTERPOL(geospres, cumh2o, nalt+1, ps(0:ntp), temph2o(0:ntp), ntp+1, errstat)
  h2o(1:ntp) = temph2o(1:ntp) - temph2o(0:ntp-1)   ! DU at each layer 
  h2o(1:nz) = h2o(1:nz) * du2mol

  RETURN  
  END SUBROUTINE GET_GEOSCHEM_H2O

END MODULE m_get_tracegas

