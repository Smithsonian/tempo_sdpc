!
module m_get_apriori_covar

  public get_apriori_covar
  private

  integer, parameter, private :: max_pathlen = 1024

contains

  ! ==============================================================
  ! Construct a priori covariance for ozone based on
  ! ozone standard deviation of Fortuin Climatology 
  ! Diagonal elements are directly from this climatology
  ! Non-diagonal elements are calculated by assuming a 
  ! correlation length (5 km for now)
  ! ==============================================================
  SUBROUTINE get_apriori_covar(month, day, lat, lon, ps, zs, ts, nz, ntp, &
       ozprof, sao3)

    USE OMSAO_precision_module 
    USE ozprof_data_module,     ONLY: use_logstate, atmos_unit, &
         min_serr, min_terr, loose_aperr, which_aperr!, which_clima
    USE OMSAO_variables_module, ONLY: atmdbdir, the_year
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline
    use prepare_atmosphere, only: get_geoschem_o3std, get_mlso3prof, &
         get_mlso3prof_single

    IMPLICIT NONE

    ! ======================
    ! Input/Output variables
    ! ======================
    INTEGER, INTENT(IN)                            :: month, day, nz, ntp
    REAL (KIND=dp), INTENT(IN)                     :: lat, lon
    REAL (KIND=dp), DIMENSION(0:nz),   INTENT(INOUT) :: zs, ps, ts
    REAL (KIND=dp), DIMENSION(nz, nz), INTENT(OUT) :: sao3
    REAL (KIND=dp), DIMENSION(nz),      INTENT(IN) :: ozprof

    ! ======================
    ! Local variables
    ! ======================
    INTEGER, PARAMETER :: mlat=18, nm=12, mlay=62
    ! 1. McPeters Clima 2. Fortuin Clima  3. McPeters + GEOS-CHEM  4. McPeters + ZM MLS
    INTEGER            :: errstat
    REAL (KIND=dp), PARAMETER       :: corrlen=6.0     ! changed from 6 km to 8 km (more uniform in a priori influence)

    REAL (KIND=dp), DIMENSION(mlay)       :: preslg
    REAL (KIND=dp), SAVE, DIMENSION(mlay) :: pres
    INTEGER, SAVE                         :: nlay, nlat, nmon

    CHARACTER (LEN=max_pathlen)           :: apfname
    !CHARACTER (LEN=2)             :: monc
    REAL (KIND=dp), SAVE, DIMENSION(nm, mlat, mlay) :: stds 
    REAL (KIND=dp), DIMENSION(mlay)           :: astd, cumastd
    REAL (KIND=dp), DIMENSION(nz)             :: zmid
    REAL (KIND=dp), DIMENSION(0:nz)           :: pslg, nstd, nstd1, ps1, zs1
    INTEGER                       :: i, j, k,  idum, nband, mnorstd, tmpntp
    INTEGER, DIMENSION(2)         :: latin, monin
    REAL (KIND=dp)                :: frac
    REAL (KIND=dp), DIMENSION(2)  :: latfrac, monfrac
    LOGICAL, SAVE                 :: first = .TRUE.

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=17), PARAMETER :: modulename = 'get_apriori_covar'  

    !IF (whichap == 1 .OR. whichap == 2) THEN
    !   IF (whichap == 1) THEN
    !      apfname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tomsv7_apcovar.dat'
    !   ELSE 
    !      WRITE(monc, '(I2.2)') month
    !      apfname = TRIM(ADJUSTL(atmdbdir)) // 'hocovar/hoh99_mean_covar_' // monc // '.dat'
    !   ENDIF 
    !
    !   IF (nz == 11) THEN   !use TOMS V7 apriori covariance 
    !      OPEN (UNIT = atmos_unit, file=apfname, status = 'old')
    !      IF (whichap == 2) READ(atmos_unit, '(A)')  
    !      READ(atmos_unit, *) ((sao3(i, j), j=1, nz), i=1, nz)
    !      CLOSE (atmos_unit) 
    !      
    !      DO i = 1,  nz
    !         IF (sao3(i, i) <= 1.d0) sao3(i, i) = 1.0
    !      ENDDO       
    !      RETURN
    !   ELSE
    !      whichap = 3
    !   ENDIF
    !ENDIF

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
    nmon = 2

    sao3 = 0.d0; astd = 0.d0
    IF (which_aperr == 1 .OR. which_aperr == 3 .OR. which_aperr == 4) THEN     
      IF (first) THEN
        nlay = 62; nlat = 18
        pres(2:62) = (/(1013.25d0*10.d0**(-1.d0*DBLE(i)/16.d0), i = 60, 0, -1)/)
        pres(1) = 0.05d0  ! about 70 km
        apfname = TRIM(ADJUSTL(atmdbdir)) // 'mpclima/llmclima_std.dat'
        OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
        READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)') 
        DO i = 1, nm 
          READ(atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')  ! read month label
          DO k = nlay, 2, -1
            READ(atmos_unit, *) idum, (stds(i, j, k), j=1, nlat) ! ppmv
          ENDDO
        ENDDO
        CLOSE(atmos_unit)
        first = .FALSE.
      ENDIF

      IF (lat <= -85.0) THEN
        nband = 1; latin(1) = 1; latfrac(1) = 1.0
      ELSE IF (lat >= 85.0) THEN
        nband = 1; latin(1) = nlat; latfrac(1) = 1.0
      ELSE
        nband = 2     ; frac = (lat + 85.0) / 10.0 + 1
        latin(1) = INT(frac); latin(2) = latin(1) + 1
        latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
      ENDIF

      DO i = 1, nband
        DO j = 1, nmon
          astd(2:nlay) =  astd(2:nlay) + stds(monin(j), latin(i), 2:nlay) * monfrac(j) * latfrac(i)
        ENDDO
      ENDDO
      astd(1) =  astd(2)
    ELSE 
      IF (first) THEN
        nlay = 20; nlat = 17
        apfname = TRIM(ADJUSTL(atmdbdir)) // 'fkclima/fortuin_o3_sdev.dat'
        pres(1:nlay) = (/0.05, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20., &
             30., 50., 70., 100., 150., 200., 300., 500., 700., 1000.0/)

        OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
        DO i = 1, nm 
          READ(atmos_unit, '(A)')  ! read month label
          READ(atmos_unit, *) ((stds(i, j, k), j=1, nlat), k=nlay, 2, -1) ! ppmv
        ENDDO
        CLOSE(atmos_unit)
        first = .FALSE.
      ENDIF

      IF (lat <= -80.0) THEN
        nband = 1; latin(1) = 1; latfrac(1) = 1.0
      ELSE IF (lat >= 80.0) THEN
        nband = 1; latin(1) = nlat; latfrac(1) = 1.0
      ELSE
        nband = 2     ; frac = (lat + 80.0) / 10.0 + 1
        latin(1) = INT(frac); latin(2) = latin(1) + 1
        latfrac(1) = latin(2) - frac; latfrac(2) = 1.0 - latfrac(1)
      ENDIF

      DO i = 1, nband
        DO j = 1, nmon
          astd(2:nlay) =  astd(2:nlay) + stds(monin(j), latin(i), 2:nlay) * monfrac(j) * latfrac(i)
        ENDDO
      ENDDO
      astd(1) =  astd(2)
    ENDIF

    cumastd(1) = 0.0
    ! convert from ppmv to partial column and accumulate
    DO i = 2, nlay
      cumastd(i) = cumastd(i-1) + (astd(i) + astd(i-1)) * &
           0.5 * (pres(i)-pres(i-1)) / 1.267
    ENDDO

    ! Interpolate to LIDORT grid
    IF (ps(nz) > pres(nlay)) pres(nlay) = ps(nz)
    preslg = LOG(pres); pslg = LOG(ps)

    CALL BSPLINE(preslg(1:nlay), cumastd(1:nlay), nlay, pslg(0:nz),&
         nstd(0:nz), nz+1, errstat)
    IF (errstat < 0) THEN
      WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat; STOP
    ENDIF

    ! Contruct the full covariance matrix for ozone (in Dobson units)
    nstd(1:nz) = nstd(1:nz) - nstd(0:nz-1)
    !print *, SUM(nstd(ntp+1:nz)) / SUM(ozprof(ntp+1:nz))
    !nstd(1:nz) =  ozprof(1:nz) * 0.5

    IF (which_aperr == 3) THEN
      ps1(0) = ps(nz)
      DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
      ENDDO
      CALL GET_GEOSCHEM_O3STD(month, lon, lat, ps1, nstd1(1:nz), nz, nz-ntp)  
      DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
      ENDDO
    ELSE IF (which_aperr == 4) THEN
      ps1(0) = ps(nz)
      DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
      ENDDO

      mnorstd = 2
      CALL get_mlso3prof(the_year, month, day, lat, nz, mnorstd, ps1(0:nz), zs1(0:nz), nstd1(1:nz), tmpntp, errstat)
      IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename, ': Error in getting MLS ozone variabilities!!!'; STOP
      ENDIF
      DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
      ENDDO
    ELSE IF (which_aperr == 5) THEN
      ps1(0) = ps(nz)
      DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
      ENDDO

      mnorstd = 2
      CALL get_mlso3prof_single(the_year, month, day, lat, nz, mnorstd, ps1(0:nz), zs1(0:nz), nstd1(1:nz), tmpntp, errstat)
      IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename, ': Error in getting MLS ozone variabilities!!!'; STOP
      ENDIF
      DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
      ENDDO
    ENDIF

    ! Loose a priori constraint (because those from climatology are sometimes too small)
    IF (loose_aperr) THEN
      DO i = 1, ntp-1 
        IF (nstd(i) / ozprof(i) < min_serr) THEN
          nstd(i) = ozprof(i) * min_serr
        ENDIF
      ENDDO

      DO i = ntp, nz 
        IF (nstd(i) / ozprof(i) < min_terr) THEN
          nstd(i) = ozprof(i) * min_terr
        ENDIF
      ENDDO

    ENDIF

    IF (use_logstate) nstd(1:nz) = nstd(1:nz)/ozprof(1:nz)

    DO i = 1, nz
      sao3(i, i)= nstd(i) ** 2.0 
    ENDDO

    ! This is based on retrieval stastistics 
    zmid = (zs(0:nz-1) + zs(1:nz)) / 2.0  
    DO i = 1, nz
      DO j = 1, i - 1
        sao3(i, j) = SQRT(sao3(i,i) * sao3(j, j)) * &
             EXP(- ABS((zmid(i)-zmid(j)) / corrlen)**2 )
        sao3(j, i) = sao3(i, j) 
      ENDDO
    ENDDO

    !WRITE(77, *) 'Se  (ozone parameters only)'
    !WRITE(77, '(11d12.4)') sao3

    RETURN
  END SUBROUTINE get_apriori_covar



end module m_get_apriori_covar
