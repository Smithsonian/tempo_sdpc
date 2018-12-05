MODULE m_get_toz

  USE OMSAO_precision_module 
  USE OMSAO_variables_module, ONLY: atmdbdir,&
      the_month, the_year, the_day, the_lon, the_lat, l3_toc_filename
  USE OMSAO_parameters_module,     ONLY: atmos_unit
  USE OMSAO_errstat_module
  USE OMSAO_he5_module

  IMPLICIT NONE
  INTEGER, PRIVATE :: i, j, nblat, nblon
  INTEGER, DIMENSION(2), PRIVATE :: latin, lonin
  REAL (KIND=dp), PRIVATE :: sumfrac
  REAL (KIND=dp), DIMENSION(2), PRIVATE :: latfrac, lonfrac
  CHARACTER(LEN=130), PRIVATE :: toz_fname
  LOGICAL, PRIVATE :: file_exist

  public get_toz
  private

  CONTAINS
  SUBROUTINE get_toz (which_toz, toz)
  
    IMPLICIT NONE
    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER , INTENT(IN) :: which_toz
    REAL (KIND=dp), INTENT(OUT) :: toz
    CHARACTER (len=4) :: ctoz
    IF (which_toz ==  1 ) THEN
       CALL get_eptoz (toz) 
    ELSE IF (which_toz == 2) THEN 
       CALL get_omtoz_zm (toz)
    ELSE IF (which_toz == 3) THEN 
       CALL get_omtoz (toz)
     IF (toz <= 0.0) THEN 
          CALL get_omtoz_zm(toz)
     ENDIF
    ELSE 
     toz = 0.0
    ENDIF
  
  END SUBROUTINE get_toz

  ! ===================================================
  ! Obtain EP TOMS monthly mean total ozone (DU) for
  !    each 1.25 by 1 region
  ! If no data is available, then use mean total ozone 
  !    from all years
  ! ====================================================

  SUBROUTINE get_eptoz(toz)
  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(OUT) :: toz

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER           :: nlat=180, nlon=288
  REAL (KIND=dp), PARAMETER    :: longrid = 1.25, latgrid = 1.0, lon0=-180.0, lat0=-90.0
  CHARACTER (LEN=2)            :: monc, yrc, dayc
  INTEGER, SAVE, DIMENSION(:,:), POINTER:: glbtoz
  LOGICAL, SAVE                        :: first = .TRUE.

  IF (first) THEN
     allocate(glbtoz(nlon, nlat))
     WRITE(dayc, '(I2.2)') the_day             ! from 9 to '09'  
     WRITE(monc, '(I2.2)') the_month           ! from 9 to '09'  
     WRITE(yrc,  '(I2.2)') MOD(the_year, 100)  ! from 1997 to '97'

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
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  toz = 0.0 ; sumfrac = 0.0
  DO i = 1, nblon
     DO j = 1, nblat
        IF (glbtoz(lonin(i), latin(j)) > 0) then 
             toz = toz + glbtoz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
             sumfrac = sumfrac + lonfrac(i)*latfrac(j)
        ENDIF
     ENDDO
  ENDDO
  IF (toz > 0 ) THEN 
      toz = toz/sumfrac
  ENDIF
  RETURN
  END SUBROUTINE get_eptoz

  SUBROUTINE get_omtoz (toz)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), INTENT(OUT) :: toz
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER           :: nlat=180, nlon=360
  REAL (KIND=dp), PARAMETER    :: longrid = 1, latgrid = 1.0, lon0=-180.0, lat0=-90.0
  INTEGER                      :: i, j, k,  fidx, lidx, sidx, eidx
  REAL                         :: dis, frac, sumfrac, toz0
  INTEGER :: fid, grid_id, status
  INTEGER (KIND=8), DIMENSION(2)       :: start, stride
  INTEGER (KIND=8), DIMENSION(2)       :: edge =(/nlon, nlat/)
  REAL, SAVE, DIMENSION(:,:),POINTER   :: glbtoz
  LOGICAL, PARAMETER                   :: do_fillin =.TRUE.
  LOGICAL, SAVE                        :: first = .TRUE.
   
  IF (first) THEN 
     ! Determine if file exists or not
     allocate(glbtoz(nlon, nlat))
     INQUIRE (FILE = TRIM(ADJUSTL(l3_toc_filename)), EXIST = file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*,*) 'please prepare OMI TO3 L3 climatological data'
        WRITE(*, *) 'GET_OMTOZ: TOC file does not exist!!!', l3_toc_filename ; STOP
     ENDIF
     start(:) = 0 ; stride(:) = 1
     fid = HE5_GDopen(TRIM(ADJUSTL(l3_toc_filename)), he5f_acc_rdonly)
     grid_id = HE5_GDattach(fid,TRIM(ADJUSTL('OMI Column Amount O3')))
     status = HE5_GDrdfld(grid_id, 'ColumnAmountO3',  start, stride, edge, glbtoz)
    
     ! fill in the bad data
     DO j = 2, nlat -1
         DO i = 1, nlon
            IF ( glbtoz(i, j) <= 0 .and. glbtoz( i, j-1) >0 .and. glbtoz(i, j+1) > 0 ) THEN 
                 glbtoz(i,j) = ( glbtoz(i,j-1) + glbtoz(i,j+1)) /2.0
            ENDIF
         ENDDO
     ENDDO    
     DO j = 1, nlat
     DO i = 2, nlon-1
            IF ( glbtoz(i, j) <= 0 .and. glbtoz( i-1, j) >0 .and. glbtoz(i+1, j) > 0 ) THEN 
                 glbtoz(i,j) = ( glbtoz(i-1,j) + glbtoz(i+1,j)) /2.0
            ENDIF
         ENDDO
     ENDDO
     IF (do_fillin == .true. ) THEN
     ! linear interpolation along the track
     DO i = 1, nlon 
      DO j = 1, nlat 
        IF (glbtoz(i,j) > 0.0 ) EXIT
      ENDDO
      fidx = j
      DO j = nlat, 1, -1 
         IF (glbtoz(i,j) > 0.0) EXIT
      ENDDO
      lidx = j
      !glbtoz(i, 1:fidx) = glbtoz(i,fidx)
      !glbtoz(i, lidx:nlat) = glbtoz(i,lidx)
      IF (fidx >= lidx ) CYCLE
      j = fidx + 1
      DO while (j <=lidx)
         IF (glbtoz(i,j) <= 0.0 ) THEN 
             sidx = j -1 ; j = j + 1
             eidx = sidx -1 
             DO WHILE ( j <=lidx)
               IF (glbtoz(i,j) > 0.0 ) THEN 
                   eidx = j ; j = j + 1 ;EXIT
               ELSE 
                   j = j + 1
               ENDIF
             ENDDO
             dis = real((eidx - sidx), KIND = dp )
             IF (dis <= 10) THEN 
                   DO k= sidx +1, eidx -1
                      frac = 1.0 - REAL((k - sidx), KIND=dp) /dis
                      glbtoz(i,k) = frac*glbtoz(i, sidx) + (1.0 - frac)*glbtoz(i,eidx)
                   ENDDO
             ENDIF
         ELSE
           j = j + 1
         ENDIF
      ENDDO ! loop of endwhile
     ENDDO  ! loop of i
     
     ! linear interpolation across the track
     Do j = 1, nlat
       DO i = 1, nlon
          IF (glbtoz(i,j) > 0.0 ) EXIT
       ENDDO
       fidx = i 
       DO i = nlon, 1, -1
          IF (glbtoz(i,j) > 0.0 ) EXIT
       ENDDO 
       lidx = i
       IF (fidx >= lidx) CYCLE
       
       i = fidx + 1
       DO WHILE ( i <= lidx)
         IF (glbtoz(i,j) <= 0.0 ) THEN 
             sidx = i -1 ; i = i + 1
             eidx = sidx -1 
             DO WHILE ( i <=lidx)
               IF (glbtoz(i,j) > 0.0 ) THEN 
                   eidx = i ; i = i + 1 ;EXIT
               ELSE 
                   i = i + 1
               ENDIF
             ENDDO
             dis = real((eidx - sidx), KIND = dp )
             IF (dis <= 10) THEN 
                   DO k= sidx +1, eidx -1
                      frac = 1.0 - REAL((k - sidx), KIND=dp) /dis
                      glbtoz(k,j) = frac*glbtoz(sidx, j) + (1.0 - frac)*glbtoz(eidx,j)
                   ENDDO
             ENDIF
         ELSE
           i = i + 1
         ENDIF
      ENDDO ! loop of endwhile
     ENDDO ! loop of latitude
     ENDIF
     first = .FALSE.
  ENDIF
  
  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
  toz = 0.0 ; sumfrac=0.0
  DO i = 1, nblon
     DO j = 1, nblat
        IF (glbtoz(lonin(i), latin(j)) > 0) THEN 
             toz = toz + glbtoz(lonin(i), latin(j)) * lonfrac(i) * latfrac(j)
             sumfrac = sumfrac + lonfrac(i)*latfrac(j)
        ENDIF
     ENDDO
  ENDDO
  toz0=toz 
  IF (toz > 0) toz = toz/sumfrac
  IF (toz > 0) toz = toz +3  
  RETURN
  END SUBROUTINE get_omtoz

  SUBROUTINE get_omtoz_zm (toz)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
   REAL (KIND=dp),INTENT(OUT)                   :: toz
   
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: ntlat=180
  REAL (KIND=dp), PARAMETER        :: latgrid=1.0, lat0=-90.
  CHARACTER (LEN=130)              :: omto3fname
  CHARACTER (LEN=2)                :: monc, dayc
  CHARACTER (LEN=4)                :: yrc
  INTEGER                          :: i, j, ib, nband
  REAL (KIND=dp)                   :: mnalt, do3

  ! Saved variables
  REAL (KIND=dp), SAVE, DIMENSION(ntlat) :: zmto3, tlats, zmalt
  LOGICAL,        SAVE                   :: first = .TRUE.

  IF (first) THEN
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc,  '(I4.4)') the_year           ! from 1997 to '1997'
        
     ! Check the availablity of MLS ozone profiles
     omto3fname =TRIM(ADJUSTL(atmdbdir)) // 'OMTO3/zm_v003_' // yrc // 'm' // monc // dayc // '.dat'
        
     ! Determine if file exists or not
     INQUIRE (FILE= omto3fname, EXIST= file_exist)
     IF (.NOT. file_exist) THEN
        WRITE(*, *) 'No Zonal Mean OMTO3 found!!!' ; stop
     ENDIF
     OPEN (UNIT = atmos_unit, file = omto3fname, status = 'unknown')
     DO i = 1, ntlat
        tlats(i) = REAL(i, KIND=dp) - 89.5
     ENDDO
     READ (atmos_unit, *)
     READ (atmos_unit, *) ((zmto3(i), zmalt(i)), i = 1, ntlat)
     CLOSE(atmos_unit)
     first = .FALSE.
   ENDIF
     
   CALL get_latfrac(ntlat,latgrid, lat0, the_lat, nband, latfrac, latin)

   toz = 0.0; sumfrac = 0.0; mnalt = 0.0
   DO ib = 1, nband
     IF ( zmto3(latin(ib)) > 5.0 ) THEN
        toz = toz + zmto3(latin(ib)) * latfrac(ib)
        mnalt = mnalt + zmalt(latin(ib)) * latfrac(ib)
        sumfrac = sumfrac + latfrac(ib)
     ENDIF
   ENDDO
   IF (toz > 0 ) THEN 
     toz  = toz / sumfrac
     mnalt = mnalt / sumfrac / 1000.0 
   ENDIF
     ! Accounting for different terrain height using approximate pressure conversion
!     IF (ps(0) < ps(1)) THEN 
!         WRITE(*, *) 'check the order of pressure , should be down to top here'
!         stop
!     ENDIF
!     do3 = 0.0
!     IF (mnalt > 0.0 ) THEN 
!     do3 = ( 1013.25 * (10.0**(-mnalt / 16.0)) - ps(0) ) / (ps(0) - ps(1)) * oz(1)
!     ENDIF
!     toz = toz - do3
  END SUBROUTINE get_omtoz_zm


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
END MODULE m_get_toz
