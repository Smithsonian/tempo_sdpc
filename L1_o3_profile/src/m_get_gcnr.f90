!JBAK 2019-10-02: for TEMPO hourly meterological and gas climatology from Chris's modeling product
! spatial interpolation :: linear interpolation depending on the distance based on lon/lat (degree)
!                          * index, fraction is determined from get_gridfrac
! time interpolation    :: linear interpolation depending on UTC between 0 and 23
!                          * interpolation is not done for month or date
! list of subroutines
!  1. gcnr_init    : lon/lat/EtaA/EtaB/ is saved internally.
!                    interpolated spres and tpres and calculated vertical layer
!                    is saved globally to gcnr%
!  2. get_gcnr_gas : interpolated gas profile is saved to gcnr%
!  3. get_gridfrac : get variables for geolocation interpolation

MODULE m_get_gcnr
  USE OMSAO_precision_module, ONLY: sp, dp
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE NETCDF
  INCLUDE 'netcdf.inc'

  !INTEGER, PARAMETER :: sp = KIND(1.0)
  !INTEGER, PARAMETER :: dp = KIND(1.0D0) !dp = KIND(1.0)
  !CHARACTER (LEN=100), PARAMETER :: atmdbdir= './gcnr/'
  ! dimension variables
  INTEGER, PARAMETER :: nlayer_gc=47, nlon_gc=1440, nlat_gc=720, nh_gc=24
  ! geolocation interpolation variables
  INTEGER :: nblon, nblat, lonin(2), latin(2)
  REAL(KIND=sp) :: lonfrac(2), latfrac(2)
  ! time interpolation variables
  INTEGER, DIMENSION (nh_gc) :: utc_gc=(/0,1,2,3,4,5,6,7,8,9,10, &
              11, 12, 13, 14, 15, 16, 17, 18, 19, 20,21,22,23/)
  INTEGER :: nbt, timein(2)
  REAL (KIND=sp) :: timefrac(2)
  ! gncr output variables
  TYPE :: gcnr_type
    REAL (KIND=sp) :: spres, tpres
    ! dimension (1:nlayer_gc)
    REAL (KIND=sp), DIMENSION (:), ALLOCATABLE :: o3, bro,ch2o, glyx,no2, so2
    REAL (KIND=sp), DIMENSION (:), ALLOCATABLE :: o3std, brostd,ch2ostd, glyxstd,no2std, so2std
    ! dimension (0:nlayer_gc)
    REAL (KIND=sp), DIMENSION (:), ALLOCATABLE :: ps
  END TYPE gcnr_type
  TYPE (gcnr_type) :: gcnr

  PUBLIC gcnr_init, get_gcnr_gas, gcnr
  PRIVATE

  Contains
  SUBROUTINE gcnr_init (the_month, the_utc, the_lon, the_lat)
    IMPLICIT NONE
    ! INPUT VARIABLES
    INTEGER, INTENT(IN) :: the_month
    REAL (KIND=sp), INTENT(IN) :: the_utc, the_lon, the_lat
    ! LOCAL VARAIBELS
    INTEGER :: i,j,k, rcode, vid, nl, nlon, nlat, idx_h
    REAL (KIND=sp):: spres,tpres
    REAL (KIND=sp),  DIMENSION (:,:,:), ALLOCATABLE:: tmp
    CHARACTER (LEN=6)   :: date
    CHARACTER (LEN=100) :: fname
    ! SAVE variables
    LOGICAL, SAVE :: first = .true.
    INTEGER, SAVE :: ncid
    REAL (KIND=sp), SAVE  :: lon0, lat0, longrid, latgrid
    REAL (KIND=sp), DIMENSION (:), ALLOCATABLE, SAVE :: lon, lat, EtaA, EtaB

    ! define dimension, simply
    nl = nlayer_gc 
    nlon = nlon_gc
    nlat = nlat_gc
    ! check out of bounds 
    IF (the_utc < 0 .or. the_utc >= 24 ) THEN
       WRITE (*,*) 'gcnr: it is out of bounds' ; stop
    ENDIF
 
    IF (first) THEN 
      IF (the_month >= 7 ) THEN 
        WRITE(date, '(A4,i2.2)' ) '2013', the_month
      ELSE  
        WRITE(date, '(A4,i2.2)' ) '2014', the_month
      ENDIF
      fname = adjustl(trim(atmdbdir))//'gcnr/PRES/gcnr_pressure_'//date//'.nc'   
      ncid = ncopn(trim(adjustl(fname)), nf_Nowrite, rcode)
      IF (rcode /= 0 ) THEN 
       WRITE(*,*) 'Error in opening ',trim(adjustl(fname)) 
      ENDIF
      
      ! read lon/lat
      ALLOCATE(lon(nlon), lat(nlat))
      CALL check ( nf_inq_varid(ncid, 'Longitude', vid))
      CALL check ( nf_get_var_real(ncid, vid, lon))
      CALL check ( nf_inq_varid(ncid, 'Latitude', vid))
      CALL check ( nf_get_var_real(ncid, vid, lat))
      ! read level-related variables
      ALLOCATE(EtaA(0:nl), EtaB(0:nl))
      CALL check ( nf_inq_varid(ncid, 'EtaA', vid))
      CALL check ( nf_get_var_real(ncid, vid, EtaA)) 
      CALL check ( nf_inq_varid(ncid, 'EtaB', vid))
      CALL check ( nf_get_var_real(ncid, vid, EtaB)) 
      allocate (gcnr%ps(0:nl))
      longrid = abs(lon(2)-lon(1))
      latgrid = abs(lat(2)-lat(1))
      lat0=lat(1)-latgrid/2.0
      lon0=lon(1)-longrid/2.0
      first = .false.
    ENDIF   

    ! check out of bounds 
    IF (the_lon  < lon(1) .or. the_lon >  lon(nlon) .or. &
       the_lat < lat(1) .or. the_lat > lat(nlat) ) THEN
       WRITE (*,*) 'gcnr: it is out of bounds ; stop'
    ENDIF
    print * , the_lon, the_lat, the_utc
    ! define hourly interpolation variables    
    nbt = 2
    timein(1) = INT(the_utc) + 1 ; timein(2) = timein(1) + 1
    if (timein(2) == 25) timein(2) = 1
    timefrac(2)= the_utc - INT(the_utc)
    timefrac(1) = 1.0 - timefrac(2)
    if (timefrac(1) == 1.0) nbt =1
    ! define geolocation interpolation variables
    CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, the_lon, the_lat, &
         nblon, nblat, lonfrac, latfrac, lonin, latin)
    gcnr%spres = 0.0 ; gcnr%tpres = 0.0
    ! getting surface pressure at given location and utc
    CALL check (nf_inq_varid(ncid, 'SurfacePressure', vid))
    DO i = 1, nblon 
     DO j = 1, nblat
      DO k = 1, nbt
     WRITE(*,'(2f10.4, i4, 3f8.2)')  lon(lonin(i)), lat(latin(j)), utc_gc(timein(k)), lonfrac(i),latfrac(j), timefrac(k)
       CALL check( nf_get_vara_real(ncid, vid, & 
           (/lonin(i),latin(j),timein(k)/),(/1,1,1/),spres)) 
       gcnr%spres = gcnr%spres + spres*lonfrac(i)*latfrac(j)*timefrac(k) 
      ENDDO
     ENDDO
    ENDDO
    ! getting surface pressure at given location and utc
    CALL check( nf_inq_varid(ncid, 'TropopausePressure', vid))
    DO i = 1, nblon 
     DO j = 1, nblat
      DO k = 1, nbt
       CALL check( nf_get_vara_real(ncid, vid, & 
           (/lonin(i),latin(j),timein(k)/),(/1,1,1/),tpres)) 
       gcnr%tpres = gcnr%tpres + tpres*lonfrac(i)*latfrac(j)*timefrac(k) 
      ENDDO
     ENDDO
    ENDDO
    ! calculate vertical levels 
    gcnr%ps = 0.0
    DO i = 0, nl
      gcnr%ps(i) = EtaA(i) + EtaB(i)*spres
    ENDDO
    ! finish     
    CALL ncclos(ncid, rcode)
    RETURN    
  END SUBROUTINE gcnr_init

  ! Trace Gas VMR (v/v) 
  ! gcnr_init should be called before reading gas profile

  SUBROUTINE get_gcnr_gas (the_month, the_utc, the_lon, the_lat, which_gas)
  IMPLICIT NONE 
  ! Input variables
  INTEGER, INTENT(IN) :: the_month
  REAL (Kind=sp), INTENT(IN) :: the_utc, the_lon, the_lat
  CHARACTER(LEN=*), INTENT(IN) :: which_gas
  ! Local variables
  INTEGER :: ncid, rcode, vid, nl, i,j, k
  REAL (kind=sp), DIMENSION (:), ALLOCATABLE :: tmp, gas,std
  CHARACTER (LEN=11) :: date
  CHARACTER (LEN=20)  :: varname1, varname2
  CHARACTER (LEN=100):: fname, subfname
  
  ! initialize
   CALL gcnr_init (the_month, the_utc, the_lon, the_lat)
   nl = nlayer_gc
   IF (adjustl(trim(which_gas)) == 'O3' ) then 
     varname1 = 'TRCO3'
     subfname=adjustl(trim(atmdbdir))//'gcnr/O3_nostd/gcnr_TRCO3_'  
   ELSE IF (adjustl(trim(which_gas)) == 'BrO' ) then 
     varname1 = 'TRCBrO'
     subfname=adjustl(trim(atmdbdir))//'gcnr/BrO/gcnr_TRCBrO_'  
   ELSE IF (adjustl(trim(which_gas)) == 'CH2O' ) then 
     varname1 = 'TRCCH2O'
     subfname=adjustl(trim(atmdbdir))//'gcnr/CH2O/gcnr_TRCCH2O_'  
   ELSE IF (adjustl(trim(which_gas)) == 'GLYX' ) then 
     varname1 = 'SPCGLYX' 
     subfname=adjustl(trim(atmdbdir))//'gcnr/GLYX/gcnr_SPCGLYX_'  
   ELSE IF (adjustl(trim(which_gas)) == 'NO2' ) then 
     varname1 = 'TRCNO2'
     subfname=adjustl(trim(atmdbdir))//'gcnr/NO2/gcnr_TRCNO2_' 
   ELSE IF (adjustl(trim(which_gas)) == 'SO2' ) then 
     varname1 = 'TRCSO2'
     subfname=adjustl(trim(atmdbdir))//'gcnr/SO2/gcnr_TRCSO2_' 
   ELSE
     WRITE (*,*) 'no exist for '//which_gas//' in gcnr climatology'
   ENDIF
   varname2=adjustl(trim(varname1))//'_variance'  

   ALLOCATE ( tmp(nl), gas(nl),  std(nl))
   gas = 0.0 ; std = 0.0
   DO k = 1, nbt 
     IF (the_month >= 7) THEN 
       WRITE(date, '(A4, i2.2, A1, I2.2, A2)') '2013', the_month,'_',utc_gc(timein(k)), '00'
     ELSE
       WRITE(date, '(A4, i2.2, A1, I2.2, A2)') '2014', the_month,'_',utc_gc(timein(k)), '00'
     ENDIF
     fname=adjustl(trim(subfname))//date//'.nc'   

     ! open file after testing  if file exists or not 
     ncid = ncopn(trim(adjustl(fname)), nf_Nowrite, rcode)
     IF (rcode /= 0) THEN 
        WRITE(*,*) 'Error in operning ',trim(adjustl(fname)) 
     ENDIF
  
     ! read/interpolate gas profile
     CALL check(nf_inq_varid(ncid,ADJUSTL(TRIM(varname1)), vid))
     DO i = 1, nblon
      DO j = 1, nblat
        CALL check(nf_get_vara_real(ncid, vid, &
             (/lonin(i),latin(j),1/), (/1, 1, nl/), tmp(:)))
        gas = gas + tmp*lonfrac(i)*latfrac(i)*timefrac(k)
      ENDDO
     ENDDO
     CALL check(nf_inq_varid(ncid,ADJUSTL(TRIM(varname2)), vid))
     DO i = 1, nblon
      DO j = 1, nblat
        CALL check(nf_get_vara_real(ncid, vid, &
             (/lonin(i),latin(j),1/), (/1, 1, nl/), tmp(:)))
        std = std + tmp*lonfrac(i)*latfrac(i)*timefrac(k)
      ENDDO
     ENDDO
     CALL ncclos(ncid, rcode)
   ENDDO
   ! copy local variable to global variable  
   IF (adjustl(trim(which_gas)) == 'O3' ) then 
     if (.not. allocated (gcnr%o3)) allocate(gcnr%o3(nl), gcnr%o3std(nl))
     gcnr%o3 = gas !; gcnr%o3std = std 
   ELSE IF (adjustl(trim(which_gas)) == 'BrO' ) then 
     if (.not. allocated (gcnr%BrO)) allocate(gcnr%BrO(nl),gcnr%BrOstd(nl))
     gcnr%BrO = gas !; gcnr%BrOstd = std 
   ELSE IF (adjustl(trim(which_gas)) == 'CH2O' ) then 
     if (.not. allocated (gcnr%CH2O)) allocate(gcnr%CH2O(nl),gcnr%CH2Ostd(nl))
     gcnr%CH2O = gas !; gcnr%CH2Ostd = std 
   ELSE IF (adjustl(trim(which_gas)) == 'GLYX' ) then 
     if (.not. allocated (gcnr%GLYX)) allocate(gcnr%GLYX(nl),gcnr%GLYXstd(nl))
     gcnr%GLYX = gas !; gcnr%GLYXstd = std 
   ELSE IF (adjustl(trim(which_gas)) == 'NO2' ) then 
     if (.not. allocated (gcnr%NO2)) allocate(gcnr%NO2(nl),gcnr%NO2std(nl))
     gcnr%NO2 = gas !; gcnr%NO2std = std 
   ELSE IF (adjustl(trim(which_gas)) == 'SO2' ) then 
     if (.not. allocated (gcnr%SO2)) allocate(gcnr%SO2(nl),gcnr%SO2std(nl))
     gcnr%SO2 = gas !; gcnr%SO2std = std 
   ENDIF
   
   print * , gcnr%spres, gcnr%tpres, gas
   ! finish 
   deallocate (gas, tmp, std)
   RETURN
  END SUBROUTINE

  SUBROUTINE check(status)
  INTEGER, intent ( in) :: status
  
  IF (status /= nf90_noerr) THEN
     print*,NF_STRERROR(status)
     WRITE(*,*)  'Errors in reading NETCDF'
     stop
  ENDIF
 END SUBROUTINE check

 SUBROUTINE get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
  lon, lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                       :: nlon, nlat
  REAL (KIND=sp), INTENT(IN)                :: lon0, lat0, lat, lon, longrid,latgrid
  INTEGER, INTENT(OUT)                      :: nblon, nblat
  INTEGER, DIMENSION(2), INTENT(OUT)        :: latin, lonin
  REAL (KIND=sp), DIMENSION(2), INTENT(OUT) :: latfrac, lonfrac

  ! ======================
  ! Local variables
  ! ======================
  REAL (KIND=sp) :: frac, lat_offset, lon_offset

  ! the starting value of center lon/lat in the grid box
  lat_offset   = lat0 + latgrid / 2.0 
  lon_offset   = lon0 + longrid  / 2.0

  ! latitude direction
  IF (lat <=lat_offset) THEN 
     nblat =1 ; latin(1) =1 ; latfrac(1) = 1.0
  ELSE IF (lat >= lat_offset + (nlat-1)*latgrid) THEN 
     nblat =1 ; latin(1)=nlat; latfrac(1) = 1.0
  ELSE
   nblat = 2; frac = (lat - lat_offset) / latgrid + 1
   latin(1) = INT(frac); latin(2) = latin(1) + 1
   latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
  ENDIF

  ! longitude direction
  IF (lon <= lon_offset) then 
     nblon =1 ; lonin(1) =1 ; lonfrac(1) = 1.0
  ELSE IF (lon >= lon_offset + (nlon-1)*longrid) THEN 
     nblon =1 ; lonin(1)=nlon; lonfrac(1) = 1.0
  ELSE
   nblon = 2; frac = (lon - lon_offset) / longrid + 1
   lonin(1) = INT(frac); lonin(2) = lonin(1) + 1
   lonfrac(1) = lonin(2) - frac; lonfrac(2) = 1.0 - lonfrac(1)
  ENDIF
  RETURN

END SUBROUTINE get_gridfrac

END MODULE m_get_gcnr
