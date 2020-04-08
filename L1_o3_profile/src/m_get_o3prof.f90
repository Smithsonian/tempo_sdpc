
MODULE m_get_o3prof

  USE OMSAO_parameters_module, ONLY: p0, du2mol
  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: atmdbdir, tabdir, atmos_unit, & 
  the_month, the_year, the_day, the_utc, the_lon, the_lat, the_time,&
  lat_min, lat_max, lon_min, lon_max, time_min, time_max
  USE ozprof_data_module,     ONLY: which_clima, which_aperr, which_toz, & 
                                    trpz, pst, ozone_above60km, use_logstate, & 
                                    min_serr, min_terr, loose_aperr, norm_tropo3
  USE OMSAO_errstat_module
  USE NETCDF
  USE m_ezspline_interpolation, ONLY: bspline
  USE m_utilities, ONLY:reverse, get_monfrac, get_latfrac, get_gridfrac

  !USE m_get_m2prof, ONLY:m2du, get_m2prof
  IMPLICIT NONE
  ! common variables used in this module
  INTEGER, PARAMETER :: which_m2 = 2, neof=72
  INTEGER, PRIVATE   :: nblat, nblon , nbmon
  INTEGER, DIMENSION(2), PRIVATE        :: latin, lonin, monin
  REAL (KIND=dp), DIMENSION(2), PRIVATE :: latfrac, lonfrac, monfrac
  CHARACTER (LEN=130), PRIVATE          :: apfname

  public get_o3prof, get_apriori_covar, get_tomsv8_clima, get_normtoz, test
  private
  
  CONTAINS
  ! main subroutine
  ! get_o3prof, get_apriori_covar 
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  ! get_normtoz
  ! get_tbprof (ozref, out_prof)
  ! get_tb
  ! get_ab
  ! get_mcprof(ozprof)
  ! get_mlprof(out, index_out)
  ! get_logan_clima(month, lon, lat, ps, ozprof, nz, ntp)
  ! GET_MIPASIG2O3(month, day, lat, xx, yy) 
  ! get_v8prof(toz, oz)
  ! get_tomsv8_clima
  ! get_geoschem_o3p
  ! get_geoschem_o3mean/std
  ! get_geoschem_o31
  ! get_mlso3prof(year, month, day, lat, nz, mnorstd, ps, oz, oz, ntp, errstat)
  ! get_mlso3prof_singple
  ! get_fortstd
  ! get_gcnr
  ! *external subroutines : BSPLINE, REVERSE in EZSPLINE_INTERPOLATION
  ! NOte : o3prof is from down to top , std is from to to down

SUBROUTINE test (error)
  LOGICAL :: norm_o3p = .false.  
  REAL (KIND=dp) :: toz
  INTEGER, PARAMETER :: numk =  65
  INTEGER :: i, ntp, error
  REAL (KIND=dp) , DIMENSION (0:numk) :: umkp, umkz 
  REAL (KIND=dp) , DIMENSION (1:numk) :: ozprof
  REAL (KIND=dp) , DIMENSION (1:numk, 1:numk) :: sao3
  
  which_clima = 12 
  which_aperr = 12
  the_lon =  -97.46540832519531
  the_lat =  36.584938049316406
  the_year  = 2018
  the_month = 10
  the_day   = 29
  umkp (0:numk)=(/(p0 * 10.0D0 ** (-REAL(i, KIND=dp)/16.D0), i=0, numk )/)
  umkz (0:numk)=(/(i, i=0, numk )/)

  call get_o3prof (numk, umkp, umkz, ntp, norm_o3p,toz, ozprof)
  call reverse(ozprof, numk)
  call reverse(umkp, numk + 1) 
  call reverse(umkz, numk + 1) 
  call get_apriori_covar(numk, umkp, umkz, ozprof,toz,  ntp, sao3)  

  write(123, *) numk, ": Number of layers", trpz, ":tropopause (km)"
  WRITE(123, *) 'ps(hPa)  zs (km)  o3 (du) SD(%) COVAR (du)'
  WRITE(123, '(f5.1, f10.4)') umkz(0), umkp(0)
  DO i = 1, numk
    WRITE(123, '(f5.1, f10.4, 70e17.5)') umkz(i), umkp(i), ozprof(i),sqrt(sao3(i, i))*100/ozprof(i), sao3(i, 1:numk)
     print * , umkz(i), sqrt(sao3(i,i))*100./ozprof(i)
  ENDDO
  stop 1
END SUBROUTINE

SUBROUTINE get_o3prof (numk, umkp, umkz, ntp, norm_o3p,toz, ozprof)

  IMPLICIT NONE

  INTEGER, PARAMETER         :: nmpref=61, nref=71, nmipas=121, nv8=11
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN) :: numk, ntp
  LOGICAL, INTENT(IN) :: norm_o3p
  REAL (KIND=dp), INTENT(IN) :: toz
  REAL (KIND=dp), DIMENSION(0:numk), INTENT(IN) :: umkp, umkz
  REAL (KIND=dp), DIMENSION(1:numk), INTENT(OUT):: ozprof ! down-top

  !============================
  ! local variables
  !============================
  INTEGER :: i, errstat, mnorstd, fidx, lidx, tmpntp
  REAL (KIND=dp) :: tmp
  REAL (KIND=dp), DIMENSION(:), allocatable :: oztmp
  REAL (KIND=dp), DIMENSION(:), allocatable :: mipasp, mipaso3 !(nmipas)
  REAL (KIND=dp), DIMENSION(0:nref)    :: ozref, refp
  REAL (KIND=DP), DIMENSION(:), allocatable :: pv8, v8oz !(0:nv8)
  REAL (KIND=dp), DIMENSION(0:numk)    :: umkoz, umkpg
  CHARACTER(10), PARAMETER :: modulename='get_o3prof'

  ! the pressure grids to call the climatology
  refp(0:nref) = (/(p0 * 10.0D0 ** (-REAL(i, KIND=dp)/16.D0), i=0, nref )/)
  refp(0) = p0

  ! Get a priori
  ozref(0:nref) = 0.0
  ! a. a priori climatology for 60-70 km from MIPAS climatology
  allocate(mipasp(nmipas), mipaso3(nmipas))
  CALL GET_MIPASIG2O3(mipasp, mipaso3)
  ozref(nmpref:nref-1) = mipaso3(nmpref:nref-1)
  ozref(nref)          = sum(mipaso3(nref:nmipas))
  ozone_above60km      = SUM( ozref(nmpref:nref)) ! used to correct to-dependent profiles 
  deallocate(mipasp, mipaso3)

  ! b. Get a priori climatology for 0-60 km (pressure altitude)
  IF (which_clima == 1) THEN 
     allocate(pv8(0:nv8), v8oz(0:nv8))
     DO i = 0, nv8 - 1
           pv8(i) = p0*2.0D0 ** (-i) 
     ENDDO
     pv8(nv8)   = refp(nref) ! reconsider before pv8(nv8) = umkp(num) : 2^(-13.5)
     pv8(0:nv8) = LOG(pv8(0:nv8) )
    !CALL get_v8prof(toz, v8oz(1:nv8))
     PRINT * , 'not well implemented' ; stop 1
  ELSE IF (which_clima == 2) THEN
     CALL GET_MCPROF (ozref(1:nmpref-1), 1) 
  ELSE IF (which_clima >= 8 .AND. which_clima <=9) THEN
     CALL GET_TBPROF (ozref(1:nmpref-1), 1)
  ELSE IF (which_clima == 10) THEN
     CALL GET_MLprof(ozref(1:nmpref-1), 1) 
  ELSE IF (which_clima == 11) THEN
    
  ELSE IF (which_clima == 12) THEN 
     IF (which_m2 == 1) THEN
        !CALL get_m2prof ('TO3',toz,nmpref-1,refp(0:nmpref-1), neof)
     ELSE IF (which_m2 == 2) THEN 
        !CALL get_m2prof ('TPP',trpz,nmpref-1,refp(0:nmpref-1), neof)
     ELSE IF (which_m2 == 3) THEN 
        !CALL get_m2prof ('LAZ',the_lat,nmpref-1,refp(0:nmpref-1), neof)
     ENDIF
     !ozref(1:nmpref-1) = m2du%o3p(1:nmpref-1)
     !print * , ozprof(1:nmpref-1), sum(ozprof(1:nmpref-1))
     !CALL GET_MLprof(ozref(1:nmpref-1), 1) 
     !print * , ozprof(1:nmpref-1), sum(ozprof(1:nmpref-1))
  ELSE IF (which_clima == 13) THEN 
     allocate(oztmp(nmpref-1))
     tmp = exp((log(pst) + log(p0))*0.5)
     tmpntp = MINVAL(MINLOC(refp(0:nmpref-1), MASK=(refp(0:nmpref-1) >  tmp)))-1
     CALL GET_tempoprof(tmpntp,refp(0:tmpntp),oztmp(1:tmpntp),1)
     ozref(1:tmpntp) = oztmp(1:tmpntp)
     CALL GET_TBPROF (oztmp(1:nmpref-1), 1)
     ozref(tmpntp) = ozref(tmpntp)*0.5 + oztmp(tmpntp)*0.5
     ozref(tmpntp+1:nmpref-1) = oztmp(tmpntp+1:nmpref-1)
     deallocate(oztmp)
  ELSE 
     CALL GET_MLprof(ozref(1:nmpref-1), 1)
  ENDIF
 
  DO i = 1, nref
     ozref(i) = ozref(i-1) + ozref(i)
  ENDDO

  ! Bondary layer correction
  IF (umkp(0) > p0) THEN !sfc
     tmp = ( umkp(0) - p0)/(refp(0)-refp(1))
     ozref(1) = ozref(1)*(1+tmp)
     refp(0) = umkp(0)
  ENDIF

  refp = LOG(refp)
  umkpg = LOG(umkp)
  ! @ Interpolate Ozone to Retrieved Grid
  CALL BSPLINE(refp, ozref, nref+1, umkpg(0:numk), umkoz(0:numk), numk+1,errstat)
  IF (errstat < 0) THEN
     WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat ; stop 1
  ENDIF
  umkoz(1:numk) = umkoz(1:numk) - umkoz(0:numk-1)
  ozprof (1:numk) = umkoz(1:numk)
  IF (which_clima == 11) THEN
     CALL get_geoschem_o3p (umkp, umkz, ozprof, numk, ntp)
  ELSE IF (which_clima == 7 ) THEN
     mnorstd = 1
     CALL get_mlso3prof_single(numk, mnorstd, umkp, umkz, ozprof, tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename, ': Error in getting MLS ozone profiles!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
  ELSE IF (which_clima == 6 .or. which_clima == 11) THEN
     mnorstd = 1
     CALL get_mlso3prof(numk, mnorstd, umkp, umkz, ozprof, tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename, ': Error in getting MLS ozone profiles!!!'
        errstat = pge_errstat_error; RETURN
     ENDIF
  ELSE IF (which_clima == 5) THEN      ! 72 x 46, Logan clima only
     CALL get_logan_clima(umkp, ozprof, numk, ntp)
  ELSE IF (which_clima == 4) THEN      ! 144 x 91, profile only
     CALL get_geoschem_o31(umkp, ozprof, numk, ntp)
  ELSE IF (which_clima == 3) THEN      ! 18 x 12, profile
     CALL get_geoschem_o3mean(umkp, ozprof, numk, ntp)
  ENDIF

  !WRITE(*,*) norm_o3p, sum(ozprof(1:numk)), toz
  IF (norm_o3p) THEN 
     IF (which_clima /=6 .AND. which_clima /= 7 ) tmpntp = ntp
     ! jbak the limitation of the vertical range to the layers below the ozone layer is better 
      IF (norm_tropo3) THEN 
       tmpntp  =  MINVAL(MAXLOC ( ozprof (1:numk) , MASK =(ozprof(1:numk) <= maxval(ozprof))))     
       fidx=1 ; lidx=tmpntp
      ELSE
       fidx=1 ; lidx=numk
      ENDIF
      CALL get_normtoz (toz, numk, ozprof, fidx, lidx)
  ENDIF
  return
END SUBROUTINE get_o3prof

SUBROUTINE get_normtoz (toz, nz,  oz, fidx, lidx)

 IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: nz, fidx, lidx
  REAL (KIND=dp),INTENT(IN)                    :: toz
  REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz

  ! ======================
  ! Local variables
  ! ======================
  REAL (KIND=dp):: res_to3
  ! 1:sfc, nz:top
   IF (fidx /= 1 .or. lidx  /= nz) THEN
     !oz(1:ntp) = oz(1:ntp) * (toz - SUM(oz(ntp+1:nz))) / SUM(oz(1:ntp))
     !oz(1:tmpntp) = oz(1:tmpntp) * (toz - SUM(oz(1:tmpntp)))
     !/SUM(oz(tmpntp+1:nz))
     !fidx = INT(ntp/2.0)
     !lidx = tmpntp
     res_to3 = 0.0
     IF (lidx < nz) res_to3 = SUM(oz(lidx+1:nz))
     IF (fidx > 1 ) res_to3 = res_to3 + SUM(oz(1:fidx-1))
     oz(fidx:lidx) = oz(fidx:lidx) * (toz - res_to3) /SUM(oz(fidx:lidx))
   ELSE
     oz(1:nz) = oz(1:nz) * toz /SUM(oz(1:nz))
   ENDIF
  RETURN
END SUBROUTINE get_normtoz


! ==============================================================
! Construct a priori covariance for ozone based on
! ozone standard deviation of Fortuin Climatology
! Diagonal elements are directly from this climatology
! Non-diagonal elements are calculated by assuming a
! correlation length (5 km for now)
! ==============================================================

SUBROUTINE get_apriori_covar( nz, ps, zs, ozprof, toz, ntp,  sao3)
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN) :: nz, ntp
  REAL (KIND=dp) :: toz
  REAL (KIND=dp), DIMENSION(0:nz),    INTENT(IN) :: ps, zs
  REAL (KIND=dp), DIMENSION(nz),      INTENT(IN) :: ozprof
  REAL (KIND=dp), DIMENSION(nz, nz),  INTENT(OUT) :: sao3 ! top-down

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: mref=60
  INTEGER            :: errstat
  REAL (KIND=dp), PARAMETER       :: corrlen=6.0     ! changed from 6 km to 8 km (more uniform in a priori influence)

  REAL (KIND=dp), DIMENSION(nz)       :: zmid
  REAL (KIND=dp), DIMENSION(0:nz)     :: pslg, nstd, nstd1, ps1, zs1
  INTEGER                             :: i, j, k,mnorstd, tmpntp, nref
  REAL (KIND=dp) :: tmp
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: oztmp, ozavg
  REAL (KIND=dp), DIMENSION(mref)     :: astd, a1, a2, a3
  REAL (KIND=dp), DIMENSION(0: mref)  :: cumastd, preslg, pres

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'get_apriori_covar'

  ! ==============================
  ! get astd
  ! ==============================
  sao3 = 0.0; astd = 0.0
  IF (which_aperr == 12 ) THEN 
    IF (which_m2 == 1) THEN 
    !  CALL get_m2prof ('TO3',toz,nz,ps(0:nz), neof)
    ELSE IF (which_m2 == 2) THEN 
    !  CALL get_m2prof ('TPP',trpz,nz,ps(0:nz), neof)
    ELSE IF (which_m2 == 3) THEN 
    !  CALL get_m2prof ('LAZ',the_lat,nz,ps(0:nz),neof)
    ENDIF
    !sao3 = m2du%sa
  ELSE ! sa is calculated from aperr
  nref = 60
  pres(1:60) = (/(1013.25*10.0**(-1.0*i/16.0), i = 59, 0, -1)/)
  pres(0) = 0.05  ! about 70 km
  IF ( which_aperr == 1) then  ! Fourtuine
     nref = 19 ! need more '
     pres(0:nref) = (/0.05, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20., &
             30., 50., 70., 100., 150., 200., 300., 500., 700., 1000.0/)
     call get_fortstd (astd(1:nref )) 
  ELSE IF (which_aperr == 2 ) THEN
     call get_mcprof(astd(1:nref), 2) 
  ELSE IF (which_aperr >= 8 .and. which_aperr <=9) THEN
     call get_tbprof (astd(1:nref),2) 
  ELSE IF (which_aperr == 10) THEN 
     call get_mlprof (astd(1:nref),2) 
  ELSE IF (which_aperr == 13) THEN
     allocate(oztmp(nref), ozavg(nref))
     tmp = exp((log(pst) + log(p0))*0.5)
     tmpntp = MAXVAL(MAXLOC(pres(0:nref), MASK=(pres(0:nref) <  tmp)))+1
     CALL GET_TBPROF (astd(1:nref), 2)
     CALL GET_tempoprof(nref-tmpntp+1,pres(tmpntp-1:nref),oztmp(tmpntp:nref),2)
     astd(tmpntp) = astd(tmpntp)*0.5 + oztmp(tmpntp)*0.5
     astd(tmpntp+1:nref) = oztmp(tmpntp+1:nref)
     deallocate(oztmp, ozavg)
  ELSE 
     call get_mlprof (astd(1:nref),2) 
  ENDIF
 ! ==============================
  ! get nstd
  ! ==============================
   !IF (return_v1) THEN
   ! IF (ps(nz) > pres(nref)) pres(nref) = ps(nz)
   !ELSE
  ! Bondary layer correction
    IF (ps(nz) > p0) then !sfc
        tmp = ( ps(nz) - p0)/(pres(nref)-pres(nref-1))
        astd(nref) = astd(nref)*(1+tmp)
        pres(nref) = ps(nz)
    ENDIF
    IF ( ps(0) < pres(0) ) pres(0) = ps(0) !top
   !ENDIF
  ! convert  partial column to accumulate
    cumastd(0) = 0.0
    DO i = 1, nref
       cumastd(i) = cumastd(i-1) + astd(i)
    ENDDO

    preslg = LOG(pres); pslg = LOG(ps)
    CALL BSPLINE(preslg(0:nref), cumastd(0:nref),nref+1, pslg(0:nz),&
         nstd(0:nz), nz+1, errstat)
    IF (errstat < 0) THEN
       WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat; stop 1
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
     CALL GET_GEOSCHEM_O3STD(ps1, nstd1(1:nz), nz, nz-ntp)
     DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
     ENDDO
  ELSE IF (which_aperr == 6) THEN 
     ps1(0) = ps(nz)
     DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
     ENDDO

     mnorstd = 2
     CALL get_mlso3prof(nz, mnorstd, ps1(0:nz), zs1(0:nz), nstd1(1:nz), tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': Error in getting MLS ozone variabilities!!!'; stop 1
     ENDIF
     DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
     ENDDO
  ELSE IF (which_aperr == 7) THEN 
     ps1(0) = ps(nz)
     DO i = 1, nz
        ps1(i) = ps(nz-i); nstd1(i) = nstd(nz-i+1)
     ENDDO

     mnorstd = 2
     CALL get_mlso3prof_single(nz, mnorstd, ps1(0:nz), zs1(0:nz), nstd1(1:nz), tmpntp, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': Error in getting MLS ozone variabilities!!!'; stop 1
     ENDIF
     DO i = 1, nz
        nstd(i) = nstd1(nz-i+1)
     ENDDO
  ENDIF

 ! Loose a priori constraint (because those from climatology are sometimes too  small)
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
        sao3(i, j) = SQRT(sao3(i,i) * sao3(j, j)) * EXP(- ABS((zmid(i)-zmid(j)) / corrlen)**2 )
        sao3(j, i) = sao3(i, j)
     ENDDO
  ENDDO
  ENDIF
  RETURN
  END SUBROUTINE get_apriori_covar

  SUBROUTINE get_geoschem_o3p  (refps, refzs, refo3, nz, ntp )

  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: nz, ntp
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: refps, refzs
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: refo3

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: nlat=91, nlon=144, nps=47
  REAL (KIND=dp), PARAMETER        :: latgrid=2, longrid=2.5, lon0=-180, lat0=-90.0

  INTEGER :: ncid, varid, i, j, fidx, lidx, errstat
  REAL, SAVE  :: lon(nlon), lat(nlat), gsps(nlon, nlat, nps)
  REAL, SAVE, DIMENSION (:,:,:), ALLOCATABLE:: gso3 !(nlon, nlat, nps)
  LOGICAL, SAVE :: first=.true.

  REAL (KIND=dp), DIMENSION(nps)   :: o3
  REAL (KIND=dp), DIMENSION(0:nps) :: ps, cumo3
  REAL (KIND=dp), DIMENSION(0:nz)  :: tempo3
  LOGICAL :: file_exist
  CHARACTER (LEN=8) :: cdate

  IF (first) THEN 
   allocate(gso3(nlon, nlat, nps))
   write(cdate, '(i4.4, i2.2, i2.2)') the_year, the_month, the_day
   apfname = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_o3p/ts_satellite.'//cdate//'.nc'
  ! Determine if file exists or not
  INQUIRE (FILE= apfname, EXIST= file_exist)
  IF (.NOT. file_exist) THEN
        WRITE(*, *) 'Warning: no geoschem o3p file!!!'
        stop 1
  ENDIF
  ! OPEN
  CALL check( nf90_open(trim(ADJUSTL(apfname)), NF90_NOWRITE, ncid))

  CALL check( nf90_inq_varid(ncid, "lon", varid) ) ! nx, ny
  CALL check( nf90_get_var(ncid, varid, lon) )

  CALL check( nf90_inq_varid(ncid, "lat", varid) ) ! nx, ny
  CALL check( nf90_get_var(ncid, varid, lat) )

  CALL check( nf90_inq_varid(ncid, "O3", varid) ) ! nx, ny
  CALL check( nf90_get_var(ncid, varid, gso3) )

  CALL check( nf90_inq_varid(ncid, "PSURF", varid) ) ! nx, ny
  CALL check( nf90_get_var(ncid, varid, gsps) )
  first = .false.
  ENDIF
  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

  o3 = 0.0 ; ps=0.0
  DO i = 1, nblon
     DO j = 1, nblat 
        o3(1:nps) = o3(1:nps) + gso3(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
        ps(0:nps-1) = ps(0:nps-1) + gsps(lonin(i), latin(j), :) * lonfrac(i) * latfrac(j)
     ENDDO
  ENDDO
  ps(nps) = 0.01

  cumo3 = 0.0
  DO i = 1, nps
     cumo3(i) = cumo3(i-1) + (o3(i)  ) * 1.0E9 / 1.25 * &
                abs(ps(i-1) - ps(i)) / 1013.25 
  ENDDO
  if (refps(0) > ps(0)) ps(0) = refps(0)
  fidx = MAXVAL(MAXLOC(refps(0:nz), MASK= (refps(0:nz) <= ps(0)))) 
  lidx = MINVAL(MINLOC(refps(0:nz), MASK= (refps(0:nz) >= refps(ntp))))-1 
  CALL BSPLINE ( ps(0:nps), cumo3(0:nps), nps+1, refps(fidx:lidx), tempo3(fidx:lidx),lidx-fidx+1 , errstat)
  refo3 (fidx+1:lidx) = tempo3(fidx+1:lidx) - tempo3(fidx:lidx-1)
  IF (any (refo3 <= 0 ) ) THEN
     print * , cumo3(2:nps)-cumo3(1:nps-1), fidx, lidx
  ENDIF
  RETURN
 END SUBROUTINE get_geoschem_o3p

 SUBROUTINE check(status)
 INTEGER, intent ( in) :: status
 IF (status /= nf90_noerr) THEN
     print *, trim(nf90_strerror(status))
     stop 1
 ENDIF
 END SUBROUTINE check

 SUBROUTINE get_geoschem_o3mean( ps, ozprof, nz, ntp)
  
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: nz, ntp
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: nlat=18, nlon=12, nalt=19
  INTEGER                          :: errstat, i, j, k,  nalt0, ntp0
  REAL (KIND=dp), PARAMETER        :: latgrid=10.0, longrid=30.0, lon0=-180.0, lat0=-90.0
  REAL (KIND=dp), DIMENSION(nalt)  :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz
  
  ! Saved variables
  REAL (KIND=dp), DIMENSION(:,:,:),ALLOCATABLE, SAVE :: geosoz
  LOGICAL                                 , SAVE :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
  CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb',&
       'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  
  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0,    .710000d0,   &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)
  
  IF (first) THEN
     allocate(geosoz(nlon, nlat, nalt))
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(the_month) // '_o3_mean.dat'
     
     OPEN(UNIT = atmos_unit, FILE = apfname, status='old')
     READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF 

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

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

  SUBROUTINE get_geoschem_o3std(ps, ozprof, nz, ntp)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: nz, ntp
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: nlat=18, nlon=12, nalt=19
  INTEGER                          :: errstat, i, j, k,  nalt0, ntp0
  REAL (KIND=dp), PARAMETER        :: latgrid=10.0, longrid=30.0, lon0=-180.0, lat0=-90.0
  REAL (KIND=dp), DIMENSION(nalt)  :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz
  
  ! Saved variables
  REAL (KIND=dp), DIMENSION(:,:,:),ALLOCATABLE, SAVE :: geosoz
  LOGICAL                                 , SAVE :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
  CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb',&
       'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  
  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0,    .710000d0,   &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)
  
  IF (first) THEN
     allocate(geosoz(nlon, nlat, nalt))
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(the_month) // '_o3_std.dat'
     OPEN(UNIT = atmos_unit, FILE = apfname, status='old')
     READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)

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

! Obtain TB hybrid oz profiles
! 2011.6.15 Jbak
! ======================================================================
  SUBROUTINE get_tbprof (ozref, out_prof)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             ::  nref = 60
  INTEGER, INTENT(IN)                            ::  out_prof ! 1 = ozref 2 =std 
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   ::  ozref

  ! ======================
  ! Local variables
  ! ======================
  INTEGER :: i, which_tb
  REAL (KIND=dp), DIMENSION(nref):: llm,ab, tmp, refz
  REAL (KIND=dp)                 :: weight1, weight2, del1, del2, meg1, meg2

  ! sketch for vertical merging
  !------------------------------- TOP
  ! 100 % ML
  !------------------------------- at trpz + del1 + meg1
  ! merging btw
  !------------------------------- at trpz + del1
  ! 100 % TB
  !------------------------------- at trpz - del2
  ! merging btw
  !------------------------------- at trpz - del2 - meg2
  ! 100 % AB
  !------------------------------- at bottom

  ! set up vertical smoothing parameters
  del1   = 1 ; del2   = 5
  meg1   = 1 ; meg2   = 1

  ! decide which TB is selected depending on the tropopause height                           
  IF ( trpz <= 14 ) THEN
     which_tb = 1 ! extra-tropicsl TB
  ELSE
     which_tb = 2 ! tropical TB
  ENDIF

  ! load climatology
  IF (out_prof == 1 ) THEN  ! mean ozone profile : bottom-up
    CALL get_mlprof(llm, 1)
    CALL get_tb (ozref, tmp,which_tb)
    CALL get_ab (AB, tmp)
    refz(1:nref) = (/(i*1.0+0.5, i = 0, nref-1 )/)
  ELSE IF (out_prof == 2) THEN  ! error profile : top-down
    CALL get_mlprof(llm,2)
    CALL get_tb (tmp,ozref,which_tb)
    CALL get_ab (tmp, AB)
    refz(1:nref) = (/(i*1.0+0.5, i= nref-1, 0,-1 )/)
  ENDIF
    IF (which_clima == 9) ozref = ab
  ! Vertical mixing. 
  DO i = 1, nref
     weight1 = 1-(abs( refz(i) - trpz )-del1)/(meg1)
     weight2 = 1-(abs( refz(i) - trpz )-del2)/(meg2)
     IF ( weight1 < 0 )  weight1 = 0
     IF ( weight1 > 1 )  weight1 = 1
     IF ( weight2 < 0 )  weight2 = 0
     IF ( weight2 > 1 )  weight2 = 1

     IF ( refz(i) >= trpz ) then ! stratosphere
       ozref(i) = ozref(i)*weight1 +LLM(i)*(1-weight1)
     ELSE ! troposphere
       ozref(i) = ozref(i)*weight2 +AB(i)*(1-weight2)
     ENDIF
  ENDDO

  IF (any(ozref(:) < 0)) then ; print * , 'error at get_tbprof' ; stop 1 ; ENDIF
  RETURN

  END SUBROUTINE get_tbprof

! ===============================================================
! Obtain TB-based oz profiles (12 month, 18 latitude bands, 80 layers : ppb)
! Variable shifht
! 2011.6.2 Jbak
! ===============================================================

  SUBROUTINE get_tb(ozref,std, which_tb)

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             ::  nref = 60
  INTEGER, INTENT(IN)                            ::  which_tb
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   ::  ozref,std

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, nmon=12, nlay=80
  REAL (KIND=dp), PARAMETER         :: lat0=-90., latgrid=10.
  REAL (KIND=dp), DIMENSION(nlay)   :: ozref0,std0 ! orignal profile
  REAL (KIND=dp), DIMENSION(0:nlay) :: cum0,cums0, refz0, zstar, tb0
  REAL (KIND=dp), DIMENSION(0:nref) :: cum,cums,refz, offset, tb

  INTEGER                           :: i, j, k,fidx, lidx, errstat
  REAL (KIND=dp)                    :: frac,fdum
  REAL (KIND=dp)                    :: meg
  REAL (KIND=dp)                    :: gravity_correct ! used for convertingunit

  LOGICAL, SAVE                     :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:),ALLOCATABLE ::ozrefs,ozrefs1, ozrefs2
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:),ALLOCATABLE ::stds,stds1, stds2
  REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat)    ::mtropz, mtropz1, mtropz2
  REAL (KIND=dp), SAVE, DIMENSION(nlat)          ::lats
  REAL (KIND=dp), SAVE, DIMENSION(nlay) :: z0

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'get_tb'

! ** load oz profiles ** !
  IF (first) THEN
     allocate(ozrefs(nmon, nlat, nlay), &
             ozrefs1(nmon, nlat, nlay), ozrefs2(nmon,nlat, nlay))

     allocate(stds(nmon, nlat, nlay), &
              stds1(nmon, nlat, nlay), stds2(nmon,nlat, nlay))

     apfname = TRIM(ADJUSTL(atmdbdir)) // 'tbclima/TB14L-5.vs' ! jbak^M
     OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
     READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')
     DO i = 1, nmon
       DO j = nlat, 1, -1
         READ(atmos_unit, *) fdum, lats(j), mtropz1(i,j) ! nsample, lat, meanztrop 
         READ(atmos_unit, *) (z0(k), ozrefs1(i, j, k),stds1(i, j, k), k = nlay,1, -1) ! ppb 
       ENDDO
     ENDDO
     CLOSE(atmos_unit)

     apfname = TRIM(ADJUSTL(atmdbdir)) // 'tbclima/TB14H-5.vs' ! jbak^M
     OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
     READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')
     DO i = 1, nmon
       DO j = nlat, 1, -1
          READ(atmos_unit, *) fdum, lats(j),mtropz2(i,j)
          READ(atmos_unit, *) (z0(k), ozrefs2(i, j, k),stds2(i, j, k), k = nlay,1, -1) ! ppb 
        ENDDO
     ENDDO
     CLOSE(atmos_unit)
    ! extratropical TB: fill 5, 15, 25 
      fidx=minval(minloc(lats,mask = (lats(1:nlat) < 35 .and. lats(1:nlat) >0)))
      lidx=minval(maxloc(lats,mask = (lats(1:nlat) < 35 .and. lats(1:nlat) >0)))

     DO i = fidx, lidx
        ozrefs1(:,i, :) =ozrefs1(:, lidx+1, :)
        stds1(:, i, :)  =stds1(:, lidx+1, :)
        mtropz1(:,i)    =mtropz1(:,lidx+1)
     ENDDO

     fidx=minval(minloc(lats,mask = (lats(1:nlat) < 0 .and. lats(1:nlat) > -35)))
     lidx=minval(maxloc(lats,mask = (lats(1:nlat) < 0 .and. lats(1:nlat) >-35)))

     DO i = fidx, lidx
        ozrefs1(:,i, :) =ozrefs1(:, fidx-1, :)
        stds1(:, i, :)  =stds1(:, fidx-1, :)
        mtropz1(:,i)    =mtropz1(:,fidx-1)
     ENDDO
     ! tropical TB: fill -85~-35 with -25, fill 35~85 with 35
     DO j = 1, 6
        ozrefs2(:, j, :) =ozrefs2(:, 7, :) !at -25
        stds2(:, j, :)   =stds2(:, 7, :)
        mtropz2(:,j)     =mtropz2(:,7)
        ozrefs2(:, j+12, :)   =ozrefs1(:, 13, :) ! at 35
        stds2(:, j+12, :)     =stds2(:, 13, :)
        mtropz2(:,j+12)       =mtropz2(:,13)
     ENDDO
     IF (any(ozrefs1 < 0) .or. any(ozrefs2 <0) ) then
         print *, 'TB clima contain -999'  ; stop 1
         stop 1
     ENDIF
    first = .FALSE.
  ENDIF
  IF (which_tb == 1 ) THEN
   ozrefs (:,:,:) = ozrefs1(:,:,:)
   stds (:,:,:)   = stds1(:,:,:)
   mtropz(:,:)    = mtropz1(:,:)
  ELSE
   ozrefs (:,:,:) = ozrefs2(:,:,:)
   stds (:,:,:)   = stds2(:,:,:)
   mtropz(:,:)    = mtropz2(:,:)
  ENDIF
    ! check
  ! ** interpolation for lat, mon** 

  CALL get_monfrac(nmon, the_month, the_day, nbmon, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nblat, latfrac, latin)

  ozref0 =0.0 ; std0 = 0.0 
  DO i = 1, nblat
        DO j = 1, nbmon
            ozref0 =  ozref0+ ozrefs(monin(j), latin(i), :) * monfrac(j) *latfrac(i)
            std0   =  std0+ stds(monin(j), latin(i), :) * monfrac(j) *latfrac(i)
            !mzt    =  mzt+mtropz(monin(j), latin(i)) * monfrac(j) * latfrac(i)
        ENDDO
  ENDDO

  !** convert tb reference into regular reference^M
  !   zs0 = [ -20, 60] is data grid^M
  !   tb  = reg - offset ^M
  !   reg = tb + offset ^M
  !   offset = (tropz - mzt)*(1-|reg(i)-tropz|/5)+mzt^M
  !   reg(i) = tb(i) +tropz+(mtz-tropz)|reg(i)-tropz|/5   ^M
  !   convert PPB to DU [ here, just use constant shift ]^M
  tb0(0:nlay-1)   = z0-0.5 ; tb0(nlay)= 60 
  refz0(0:nlay)   = tb0(0:nlay)+trpz
  zstar(0:nlay)   = 1.0/(10**((refz0)/16.0))
  cum0(:) = 0.0 ; cums0(:)=0.0

  DO i = 1, nlay
    gravity_correct = (6367. / (6367. + refz0(i)+0.5 ))**2.
    ozref0(i) = ozref0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
    cum0(i)   = cum0(i-1) + ozref0(i)
    std0(i)   = std0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
    cums0(i)   = cums0(i-1) + std0(i)
  ENDDO

  ! derive variable shifht from LLM grid algitude covering 0 to 60 km^M
  refz(0:nref)   = (/(i*1.0, i= 0, 60 )/)
  tb(0:nref)     = refz(0:nref)- trpz

  !  DO i = 0, nref
  !     IF ( abs(refz(i)-trpz ) <= 6. ) then ^M
  !        offset(i) = (trpz-mzt)*(1-abs(refz(i)-trpz)/6.)+mzt^M
  !     ELSE
  !        offset(i) = mzt
  !     ENDIF^M
  !  ENDDO
  !  tb(0:nref) = refz(0:nref)-offset(0:nref) ^M
  !ENDIF

 IF (tb(0) < tb0(0) .or. tb(nref) > tb0(nlay) ) then
      print * , 'check boundary condition in TB clim'
      print * , TB(0), tb0(0), tb(nref), tb0(nlay), trpz ; stop 1
  ENDIF
  CALL BSPLINE(tb0, cum0, nlay+1, tb, cum, nref+1, errstat)
  CALL BSPLINE(tb0, cums0, nlay+1, tb, cums, nref+1, errstat)
  IF (errstat < 0) THEN
    WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat ; stop 1
  ENDIF

  ozref(1:nref) = cum(1:nref)-cum(0:nref-1)
  std(1:nref)   = cums(1:nref)-cums(0:nref-1)
  IF (any(ozref(:) < 0)) then
     ozref(:) = -999 ; std(:) = -999 ; print *, 'TB <0' ; return
  endif
  CALL REVERSE(STD(1:nref), nref)

  RETURN
  END SUBROUTINE get_tb

  SUBROUTINE get_ab (ozref,std)
  ! remove the option of selecting tropical, extratropical AB
  ! just use AB all by Jbak 2017-07-11
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, PARAMETER                             ::  nref = 60
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)   ::  ozref,std
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, nmon=12, nlay=60
  REAL (KIND=dp), parameter::  latgrid=10., lat0=-90
  REAL (KIND=dp), DIMENSION(nlay)   :: ozref0,std0
  REAL (KIND=dp), DIMENSION(0:nlay) :: refz, zstar
  REAL (KIND=dp)                    :: gravity_correct

  REAL (KIND=dp)                :: frac,fdum
  INTEGER                       :: i, j, k, errstat

  LOGICAL, SAVE                 :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE::ozrefs
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE::stds
  REAL (KIND=dp), SAVE, DIMENSION(nlay) :: zs0

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'get_ab'      

  ! ** load std profiles ** !
  IF (first) THEN
     allocate( ozrefs(nmon, nlat, nlay), stds(nmon, nlat, nlay))
     ! AB clima. with all sonde profiles
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'tbclima/ABall.ns' ! jbak^M
     OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
     READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')

     DO i = 1, nmon
     DO j = nlat, 1, -1
         READ(atmos_unit, *) fdum, fdum
         READ(atmos_unit, *) (zs0(k), ozrefs(i, j, k),stds(i, j, k), k=nlay, 1,-1)! ppb
     ENDDO
     ENDDO
     CLOSE(atmos_unit)
     first = .FALSE.
  ENDIF

  ! ** interpolation for lat, mon** ! 
  CALL get_monfrac(nmon, the_month, the_day, nbmon, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nblat, latfrac, latin)

  ozref0(:) =0.0 ; std0(:)=0.0
  DO i = 1, nblat
        DO j = 1, nbmon
           ozref0=  ozref0+ ozrefs(monin(j), latin(i), :) * monfrac(j) *latfrac(i)
           std0=  std0+ stds(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
        ENDDO
  ENDDO

  ! ** convert ppb into DU ** ! 
  refz(0:nref) = (/(i*1.0, i= 0, 60 )/)
  zstar(0:nref) = 1.0/(10**((refz)/16.0))
  DO i = 1, nlay
   gravity_correct = (6367. / (6367. + refz(i)+0.5 ))**2.
   ozref0(i) = ozref0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
   std0(i) = std0(i)*( zstar(i-1)-zstar(i)) / ( 1.25 * gravity_correct)
  ENDDO

  ozref(:) = ozref0(:)
  std(:)   = std0(:)
  CALL REVERSE(STD(1:nref), nref)

  RETURN
END SUBROUTINE get_ab

! ===============================================================
! Obtain TOMS V8 ozone profiles (12 month, 18 latitude bands,
!   3-10 profiles with total ozone at a step of 50 DU
! ===============================================================
SUBROUTINE get_v8prof(toz, oz)

  IMPLICIT NONE

  INTEGER, PARAMETER                           :: nl = 11
  ! ======================
  ! Input/Output variables
  ! ======================

  REAL (KIND=dp), INTENT(INOUT)                :: toz
  REAL (KIND=dp), DIMENSION(nl), INTENT(OUT)   :: oz

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, maxprof=10, nmon=12
  REAL (KIND=dp), parameter :: latgrid=10, lat0=-90
  CHARACTER (LEN=200)                                :: line

  REAL (KIND=dp) :: frac, fdum, maxoz,minoz
  INTEGER        :: i, j, ib, profin, nprof, im

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: ozprofs
  INTEGER,        SAVE, DIMENSION(:,:), ALLOCATABLE:: nprofs
  LOGICAL,        SAVE  :: first = .TRUE.

! ** load oz profiles ** !

  IF (first) THEN
   allocate(ozprofs(nmon, nlat, maxprof, nl), nprofs(nmon, nlat))
  apfname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tomsv8_ozone_clima.dat'
  OPEN (UNIT = atmos_unit, file= apfname, status = 'unknown')
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
  first = .FALSE.
  ENDIF

  CALL get_monfrac(nmon, the_month, the_day, nbmon, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0, the_lat, nblat, latfrac, latin)

  oz = 0.0
  DO im = 1, nbmon
        DO ib = 1, nblat
           nprof = nprofs(monin(im), latin(ib))
           minoz = SUM(ozprofs(monin(im), latin(ib), 1, :))
           maxoz = SUM(ozprofs(monin(im), latin(ib), nprof, :))

           IF (toz < minoz) THEN
              WRITE(*,*) 'Warning: no a priori profile available!!!'
              oz  = oz + ozprofs(monin(im), latin(ib), 1, :) * toz / minoz *latfrac(ib)
           ELSE IF (toz > maxoz) THEN
              WRITE(*,*) 'Warning: no a priori profile available!!!'
              oz = oz + ozprofs(monin(im), latin(ib), nprof, :) * toz / maxoz *latfrac(ib)
           ELSE
              profin = INT ((toz - minoz ) / 50.0)+1
              IF (profin == 0) THEN
                 profin = 1
              ELSE IF (profin == nprof) THEN
                 profin = profin - 1
              ENDIF

              frac = 1.0 - (toz - (minoz + (profin-1) * 50.0)) / 50.0
              oz = oz + latfrac(ib) * monfrac(im) * (frac * ozprofs(monin(im),latin(ib), profin, :) &
                   + (1.0 - frac) * ozprofs(monin(im), latin(ib), profin+1, :))           
           ENDIF
        ENDDO
  ENDDO
  RETURN
  END SUBROUTINE get_v8prof

  SUBROUTINE get_mcprof(ozref, which_out)

  IMPLICIT NONE

  INTEGER, PARAMETER                           :: nref = 60
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER , INTENT (IN)               :: which_out ! 1=o3p, 2=std
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT) :: ozref

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, nmon=12, nlevel=61 ! o3 (DU) for 60 layer, std (mr) for 61 level
  INTEGER :: i, j, im,ib
  REAL (KIND=dp), parameter :: latgrid=10, lat0=-90
  REAL (KIND=dp)                                  :: frac, idum
  REAL (KIND=dp), DIMENSION(nlevel)               :: std0, pres
  REAL (KIND=dp), DIMENSION(nref)                 :: std
  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: ozrefs
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: stds
  LOGICAL,        SAVE                            :: first = .TRUE.

! ** load oz profiles ** !
  IF (first) THEN
     allocate(ozrefs(nmon, nlat, nref), stds(nmon, nlat, nlevel))
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'mpclima/llmclima_prof.dat'
     OPEN (UNIT = atmos_unit, file= apfname, status = 'unknown')
     DO im = 1, nmon
        READ (atmos_unit, *)
        DO i = nref, 1, -1
           READ (atmos_unit, *) ozrefs(im, :, i) ! down to top
        ENDDO
     ENDDO
     CLOSE (atmos_unit)
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'mpclima/llmclima_std.dat'
     OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
     READ (atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')
     DO im = 1, nmon
        READ(atmos_unit, '(A)') ;  READ(atmos_unit, '(A)')  ! read month label
        DO i = nlevel, 1, -1
              READ(atmos_unit, *) idum, (stds(im, j, i), j=1, nlat) ! ppm top to down
        ENDDO
     ENDDO
     CLOSE(atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_monfrac(nmon, the_month, the_day, nbmon, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nblat, latfrac, latin)
  
  std0  = 0.0
  ozref = 0.0
  DO im = 1, nbmon
     DO ib = 1, nblat
        ozref = ozref + latfrac(ib) * monfrac(im) * ozrefs(monin(im), latin(ib),:)
        std0(1:nlevel) = std0(1:nlevel) + latfrac(ib)*monfrac(im)*stds(monin(im), latin(ib),1:nlevel)
     ENDDO
  ENDDO
  ! convert ppm into DU
  
  pres(1) = 0.05  ! about 70 km
  pres(2:61) = (/(1013.25*10.0**(-1.0*i/16.0), i = 59, 0, -1)/)
  DO i = 1, nref
     std(i) = (std0(i+1) + std0(i))*0.5 *(pres(i+1) - pres(i))/ 1.267
  ENDDO
  IF (which_out == 2) ozref(1:nref) = std(1:nref)  

  RETURN
  END SUBROUTINE get_mcprof

!  DU table 
!  lat [-85, 85]
!  mon [1, 12]
!  lat [0-1, 64-65, 66-90]
  SUBROUTINE get_mlprof(out, index_out)

  IMPLICIT NONE

  INTEGER, PARAMETER                           :: nref = 60
  ! ======================
  ! Input/Output variables
  ! ======================

  INTEGER, INTENT(IN) :: index_out
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT) :: out
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=18, nmon=12, nlay=66
  REAL (KIND=dp), parameter :: latgrid=10, lat0=-90
  CHARACTER (LEN=10) :: cdum
  REAL (KIND=dp)     :: frac
  INTEGER :: i, j, im,ib
  REAL (KIND=dp),DIMENSION(nlay) :: ozref0,std0, pres
  REAL (KIND=dp),DIMENSION(nref) :: std, ozref
  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: ozrefs, stds
  LOGICAL,        SAVE                            :: first = .TRUE.
! ** load oz profiles ** !
  IF (first) THEN
     ! LOAD ozone DU table
     allocate( ozrefs(nmon, nlat, nlay), stds(nmon, nlat, nlay))
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'MLclima/ML_du_table.dat'
     OPEN (UNIT = atmos_unit, file= apfname, status = 'unknown')

     DO ib = 1, nlat
        READ (atmos_unit, *) ;READ (atmos_unit, *)
        READ (atmos_unit, *) ;READ (atmos_unit, *)
        DO i =  1,nlay
           READ (atmos_unit, '(a10, 12f7.3)') cdum, ozrefs(:, ib, i) !bottom-top
        ENDDO
     ENDDO
     CLOSE (atmos_unit)

     ! LOAD STD ppmv table
     apfname = TRIM(ADJUSTL(atmdbdir)) //'MLclima/ML_ppmv_stats.dat'
     OPEN (UNIT = atmos_unit, file= apfname, status = 'unknown')
      READ (atmos_unit, *) ;READ (atmos_unit, *)
     DO ib = 1, nlat
        READ (atmos_unit, *) ;READ (atmos_unit, *)
        READ (atmos_unit, *)
        DO i =  nlay,1, -1
           READ (atmos_unit, '(a4, 12f7.3)') cdum, stds(:, ib, i) ! top-bottom
        ENDDO
     ENDDO
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_monfrac(nmon, the_month, the_day, nbmon, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nblat, latfrac, latin)

  ozref0 = 0.0; std0=0.0
  DO im = 1, nbmon
     DO ib = 1, nblat
        ozref0 = ozref0 + latfrac(ib) * monfrac(im) * ozrefs(monin(im),latin(ib), :)
        std0 = std0 + latfrac(ib) * monfrac(im) * stds(monin(im), latin(ib), :)
     ENDDO
  ENDDO
  ozref(:) = ozref0(1:nref)

! convert ppmb to DU for std profile
  pres(1) = 0.05  ! about 70 km
  pres(2:nlay) = (/(1013.25*10.0**(-1.0*i/16.0), i = nlay-2, 0, -1)/)

  DO i = nlay-nref, nlay-1
    std(i-5)= (std0(i+1) + std0(i))*0.5 *(pres(i+1) - pres(i))/ 1.267
  ENDDO

  IF ( index_out == 1 ) out(:) = ozref(:)
  IF ( index_out == 2 ) out(:) = std(:)
  RETURN
  END SUBROUTINE get_mlprof

! Use MIPAS IG2 Temperature Profile cimatology 
! 121 levels (pressre altitude from 120 km to 0 km), 4 months (1,4,7,10)
! and 6 latitude bands (-75, -45, -10, 10, 45, 75)
  SUBROUTINE GET_MIPASIG2O3(xx, yy)
  
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

  ! saved variables
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: profs
  REAL (KIND=dp), SAVE, DIMENSION(nl)             :: pres0
  LOGICAL,        SAVE                            :: first = .TRUE.

  REAL (KIND=dp), DIMENSION(0:nlat)               :: temp
  REAL (KIND=dp)                                  :: frac, fmon
  INTEGER                                         :: ib, nb, nm, im, i, nheader

  IF (first) THEN
     allocate(profs(nl, nlat, nmon))
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'mipasprof/MIPAS_IG2_O3clima.dat'
     nheader = 9
     
     OPEN (UNIT = atmos_unit, file= apfname, status = 'unknown')
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

  END SUBROUTINE GET_MIPASIG2O3

  SUBROUTINE get_geoschem_o31(ps, ozprof, nz, ntp)  

  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)         :: nz, ntp
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER               :: nlat=91, nlon=144, nalt=19
  REAL (KIND=dp), PARAMETER        :: longrid = 2.5, latgrid = 2.0, lon0=-181.25, lat0=-91.0
  INTEGER                          :: errstat, i, j, k, ntp0, nalt0

  REAL (KIND=dp), DIMENSION(nalt)  :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)  :: tempoz

  REAL (KIND=dp), SAVE, DIMENSION(:,:,:),ALLOCATABLE :: geosoz
  LOGICAL, SAVE                                  :: first = .TRUE.

  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(0:nalt), PARAMETER:: pres = (/1.0d0,          &
       .987871d0, .954730d0, .905120d0, .845000d0, .78d0, .710000d0,      &
       .639000d0, .570000d0, .503000d0, .440000d0,.380000d0, .325000d0,   &
       .278000d0, .237954d0, .202593d0, .171495d0, .144267d0, .121347d0,  &
       .102098d0/)

  REAL (KIND=dp), DIMENSION(0:nalt)           :: geospres, cumoz
  CHARACTER (LEN=3), DIMENSION(12) :: months = (/'jan', 'feb','mar', 'apr', &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  
  IF (first) THEN
     allocate (geosoz(nlon, nlat, nalt))
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'geoschem_tropclima/' // months(the_month) // '_o3_avg.dat'
     OPEN(UNIT = atmos_unit, FILE = apfname, status='old')
     READ(atmos_unit, *) (((geosoz(i, j, k), k = 1, nalt), j = 1, nlat), i = 1, nlon)
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
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


  SUBROUTINE get_logan_clima( ps, ozprof, nz, ntp)  
  IMPLICIT NONE

  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                             :: nz, ntp
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)     :: ps
  REAL (KIND=dp), DIMENSION(nz),   INTENT(INOUT)  :: ozprof

  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER        :: nlat=46, nlon=72, nalt=13
  REAL (KIND=dp), PARAMETER :: longrid = 5.0, latgrid = 4.0, lon0=-180.0, lat0=-92.0
  INTEGER                   :: errstat, i, j, k, ntp0

  REAL (KIND=dp), DIMENSION(nalt)             :: gprof
  REAL (KIND=dp), DIMENSION(0:nz)             :: tempoz

  REAL (KIND=dp), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: geosoz
  LOGICAL, SAVE                                     :: first = .TRUE.

  ! Correct coordinates
  REAL (KIND=DP), DIMENSION(1:nalt), PARAMETER:: pres = (/1000., 900., &
       800., 700., 600., 500., 400., 300., 250., 200., 150., 125., 100./)

  REAL (KIND=dp), DIMENSION(1:nalt)           :: cumoz, presmod
  CHARACTER (LEN=3), DIMENSION(12)            :: months = (/'jan', 'feb', &
       'mar', 'apr',  'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=2)                           :: monc
  
  IF (first) THEN
     allocate(geosoz(nlon, nlat, nalt))
     WRITE(monc, '(I2.2)') the_month
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'logan_clima/ozone.13.4x5.' // monc

     OPEN(UNIT = atmos_unit, FILE = apfname, status='old')
     READ(atmos_unit, '(9E10.3)') geosoz
     CLOSE (atmos_unit)
     first = .FALSE.
  ENDIF

  CALL get_gridfrac(nlon, nlat, longrid, latgrid, lon0, lat0, &
       the_lon, the_lat, nblon, nblat, lonfrac, latfrac, lonin, latin)
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


! =====================================================================
! Obtain AURA MLS zonal mean ozone profiles and its standard deviations
! (quality flags applied) 0.1-215 mb (i.e., 10-64 km), 36 latitude bins
! updated for MLSv4.2 by jbak on 2017-08-30
! =====================================================================
  SUBROUTINE get_mlso3prof(nz, mnorstd, ps, zs, oz, ntp, errstat)
  IMPLICIT NONE

  INTEGER, PARAMETER                           :: ml = 37, mlat=36
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          :: nz, mnorstd
  INTEGER, INTENT(OUT)                         :: errstat, ntp
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
  REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz
   
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=2)                :: monc, dayc
  CHARACTER (LEN=4)                :: yrc
  LOGICAL                          :: file_exist
  INTEGER                          :: i, j, ib, nband, fidx, lidx, ntmpl, sl, el
  REAL (KIND=dp), DIMENSION (ml)   :: tmpoz, tmpozstd, ratio
  REAL (KIND=dp), DIMENSION (0:ml) :: cumoz
  REAL (KIND=dp), DIMENSION (0:nz) :: tmpcumoz, tmps
  REAL (KIND=dp), DIMENSION (2)    :: latfrac
  REAL (KIND=dp)                   :: sumfrac
  INTEGER,        DIMENSION (2)    :: latin

  ! Saved variables
  INTEGER, SAVE                             :: nlat, nl
  REAL (KIND=dp), SAVE, DIMENSION(:,:), ALLOCATABLE :: mlsprofs, mlstds
  REAL (KIND=dp), SAVE, DIMENSION(0:ml)     :: mlsps
  REAL (KIND=dp), SAVE, DIMENSION(mlat)     :: mlslats
  LOGICAL,        SAVE                      :: first = .TRUE.

  errstat = 0
  IF (first) THEN
     allocate (mlsprofs(mlat, ml), mlstds(mlat, ml))
     WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
     WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
     WRITE(yrc,  '(I4.4)') the_year           ! from 1997 to '1997'
     
     ! Check the availablity of MLS ozone profiles
     !apfname =TRIM(ADJUSTL(atmdbdir)) // 'MLSO3/zm_v02_' // yrc // monc // dayc // '.dat'
     apfname =TRIM(ADJUSTL(atmdbdir)) // 'MLSO3V4/zm_v04_' // yrc // monc // dayc // '.dat'
      
     ! Determine if file exists or not
     INQUIRE (FILE= apfname, EXIST= file_exist) 
     IF (.NOT. file_exist) THEN
        WRITE(*,*) apfname
        WRITE(www_lun, *) 'No MLS ozone profile found!!!'; errstat = -1; stop 1
     ENDIF
     
     OPEN (UNIT = atmos_unit, file = apfname, status = 'unknown')
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

  IF (the_lat <= mlslats(1)) THEN
     nband = 1; latin(1) = 1; latfrac(1) = 1.0
  ELSE IF (the_lat >= mlslats(nlat)) THEN
     nband = 1; latin(1) = nlat; latfrac(1) = 1.0
  ELSE
     nband = 2 
     DO i = 2, nlat
        IF ( the_lat <= mlslats(i) ) THEN
           latin(1) = i - 1; latin(2) = i
           latfrac(2) = (the_lat - mlslats(i - 1)) / (mlslats(i) - mlslats(i-1))
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
  !jbak
  DO i = 1, nl
     IF (exp(mlsps(i)) <= 215.0) THEN
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
  !print *,exp(tmps(fidx-1)), EXP(tmps(fidx)), EXP(tmps(lidx))
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
  SUBROUTINE get_mlso3prof_single(nz, mnorstd, ps, zs, oz, ntp, errstat)
  IMPLICIT NONE

  INTEGER, PARAMETER                           :: ml = 37
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)                          ::  nz, mnorstd
  INTEGER, INTENT(OUT)                         :: errstat, ntp
  REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: ps, zs
  REAL (KIND=dp), DIMENSION(nz), INTENT(INOUT) :: oz
   
  ! ======================
  ! Local variables
  ! ======================
  CHARACTER (LEN=2)                :: monc, dayc
  CHARACTER (LEN=4)                :: yrc
  LOGICAL                          :: file_exist
  INTEGER                          :: i, j, fidx, lidx, ntmpl, sl, el, nl, theprof, ios, nm
  REAL (KIND=dp), DIMENSION (ml)   :: tmpoz, ratio, mlsprof, mlstd
  REAL (KIND=dp), DIMENSION (0:ml) :: cumoz, mlsps
  REAL (KIND=dp), DIMENSION (0:nz) :: tmpcumoz, tmps
  REAL (KIND=dp)                   :: tmplon, tmplat, tmpsza, tmptime

  errstat = 0

  WRITE(monc, '(I2.2)') the_month          ! from 9 to '09' 
  WRITE(dayc, '(I2.2)') the_day            ! from 9 to '09'     
  WRITE(yrc,  '(I4.4)') the_year           ! from 1997 to '1997'
  
  ! Check the availablity of MLS ozone profiles
  apfname =TRIM(ADJUSTL(atmdbdir)) // 'MLSO3/mlso3_v02_' // yrc // monc // dayc // '.dat'
    
  ! Determine if file exists or not
  INQUIRE (FILE= apfname, EXIST= file_exist)
  IF (.NOT. file_exist) THEN
     WRITE(www_lun, *) 'No MLS ozone profile found!!!'; errstat = -1; RETURN
  ENDIF

  OPEN (UNIT = atmos_unit, file = TRIM(tabdir)//'INP/mlsprof_index.inp', status = 'unknown', IOSTAT=ios)
  IF (ios /= 0) THEN
     WRITE(www_lun, *) 'Do not know which profile to choose!!!'; errstat = -1; RETURN
  ELSE
     READ (atmos_unit, *) theprof; CLOSE (atmos_unit)
  ENDIF
   
  OPEN (UNIT = atmos_unit, file = apfname, status = 'unknown')
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

  CALL BSPLINE(mlsps(0:nl), cumoz(0:nl), nl+1, tmps(fidx:lidx), &
       tmpcumoz(0:ntmpl-1), ntmpl, errstat)
  IF (errstat < 0) THEN
     WRITE(www_lun, *) 'GET_MLSO3PROF: BSPLINE error, errstat = ', errstat; RETURN
  ENDIF
  oz(fidx+1:lidx) = tmpcumoz(1:ntmpl-1) - tmpcumoz(0:ntmpl-2)
  ntp = fidx
  
  RETURN
  END SUBROUTINE get_mlso3prof_single

! ===============================================================
! Obtain TOMS V8 ozone profiles (12 month, 18 latitude bands,
!   3-10 profiles with total ozone at a step of 50 DU
! ===============================================================
  SUBROUTINE get_tomsv8_clima(month, day, lat, toz, nl, ps, apoz, oz, errstat)

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
  INTEGER, PARAMETER :: nmon=12, nlat=18, maxprof=10
  REAL (KIND=dp), PARAMETER :: lat0=-90.0, latgrid=10.0
  CHARACTER (LEN=3), DIMENSION(12)  :: months = (/'jan', 'feb','mar', 'apr', &
       'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/)
  CHARACTER (LEN=200)                                :: line

  ! saved variables
  !REAL (KIND=dp), SAVE, DIMENSION(nmon, nlat, maxprof, nl0) :: ozprofs
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: ozprofs
  INTEGER,        SAVE, DIMENSION(nmon, nlat)               :: nprofs
  REAL (KIND=dp), SAVE, DIMENSION(0:nl0)                    :: pv80
  LOGICAL,        SAVE                                      :: first = .TRUE.

  REAL (KIND=dp)                                            :: frac, fdum, maxoz, minoz
  REAL (KIND=dp), DIMENSION(nl0)                            :: oz0
  REAL (KIND=dp), DIMENSION(0:nl0)                          :: cum0
  REAL (KIND=dp), DIMENSION(0:nl)                           :: logps, cum
  REAL (KIND=dp), DIMENSION(2)                              :: latfrac, monfrac
  INTEGER,        DIMENSION(2)                              :: latin, monin
  INTEGER :: i, j, ib, profin, nprof, im

  CHARACTER (LEN=16), PARAMETER :: modulename = 'get_tomsv8_clima'

  IF (first) THEN
     allocate (ozprofs(nmon, nlat, maxprof, nl0))
     ! read the TOMS V8 profiles
     apfname = TRIM(ADJUSTL(atmdbdir)) // 'v8clima/tomsv8_ozone_clima.dat'
     OPEN (UNIT = atmos_unit, file= apfname, status = 'unknown')
     
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

  CALL get_monfrac(nmon, month, day, nbmon, monfrac, monin)
  CALL get_latfrac(nlat, latgrid, lat0, lat, nblat, latfrac, latin)

  oz0 = 0.0
  DO im = 1, nbmon
     DO ib = 1, nblat   
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

! ===============================================================
! Obtain fortstd profiles (12 month, 17 latitude bands, 19 levels)
! ===============================================================
  SUBROUTINE get_fortstd(std)
  IMPLICIT NONE
  INTEGER, PARAMETER                            ::  nref = 19
  ! ======================
  ! Input/Output variables
  ! ======================
  REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)  :: std
  ! ======================
  ! Local variables
  ! ======================
  INTEGER, PARAMETER :: nlat=17, nmon=12, nlay=20
  REAL (KIND=dp), PARAMETER :: lat0=-90.0, latgrid=10.0
  REAL (KIND=dp), DIMENSION(nlay) :: std0, pres
  INTEGER                         :: i, j, k
  LOGICAL, SAVE                   :: first = .TRUE.
  REAL (KIND=dp), SAVE, DIMENSION(:,:,:),ALLOCATABLE :: stds 

  
! ** load std profiles ** !
  IF (first) THEN
        allocate(stds(nmon, nlat, nlay))
        apfname = TRIM(ADJUSTL(atmdbdir)) // 'fkclima/fortuin_o3_sdev.dat'
        OPEN (UNIT = atmos_unit, file=apfname, status = 'unknown')
        DO i = 1, nmon
           READ(atmos_unit, '(A)')  ! read month label
           READ(atmos_unit, *) ((stds(i, j, k), j=1, nlat), k=nlay, 1, -1) ! ppmv
        ENDDO
        CLOSE(atmos_unit)
        first = .FALSE.
 ENDIF

! ** interpolation for lat, mon** ! 

  CALL get_monfrac(nmon, the_month, the_day, nbmon, monfrac, monin)
  CALL get_latfrac(nlat,latgrid, lat0,the_lat, nblat, latfrac, latin)

  std0 =0.0
  DO i = 1, nblat
     DO j = 1, nbmon
           std0 =  std0 + stds(monin(j), latin(i), :) * monfrac(j) * latfrac(i)
     ENDDO
  ENDDO

  pres(1:nlay) = (/0.05, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20., &
             30., 50., 70., 100., 150., 200., 300., 500., 700., 1000.0/) 
  ! convert to Du
  DO i = 1, nref
     std(i) = (std0(i+1) + std0(i))*0.5 *(pres(i+1) - pres(i))/ 1.267
  ENDDO
  RETURN

  END SUBROUTINE get_fortstd 

  SUBROUTINE get_tempoprof(nref,ps, ozprof, which_out)
  use clim_module
  use tell_module
  use OMSAO_variables_module, ONLY: do_geoloc_init
  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
   INTEGER, INTENT(IN) :: nref, which_out
   REAL (KIND=dp), DIMENSION(0:nref), INTENT(IN) :: ps
   REAL (KIND=dp), DIMENSION(nref), INTENT(OUT)  :: ozprof
  ! ======================
  ! local variables
  ! ======================
   TYPE (clim_pres_type), SAVE :: cpt
   TYPE (clim_pres_bounds_type), SAVE :: bounds
   TYPE (clim_species_type), SAVE :: cst
   INTEGER, SAVE  :: nl0
   INTEGER :: errstat, nl
   INTEGER :: year, month, day
   REAL (KIND=r8) :: hour
   REAL (KIND=r4) :: lon, lat
   REAL (KIND=r4) :: hour_f, lon_f, lat_f
   REAL (kind=r4), dimension(:), allocatable, SAVE :: pres, vmr,vmr_stddev
   REAL (KIND=r4), dimension(:), allocatable ::  partial_column, tmp
   character (len=6), PARAMETER :: clim_db_molecule_name  ='O3    '
   logical :: is_reord
   INTEGER :: i
   REAL(KIND=dp), DIMENSION(:), allocatable :: pstmp
  ! ======================
  ! module name
  ! ======================
   character (len=13), parameter :: modulename = 'get_tempoprof'
   
   errstat=0
   IF (do_geoloc_init) THEN  
   !--------------------------------------
   ! initialize climatology
   !-------------------------------------
   !@ set bounds
   if (time_max - time_min > 86400.0) then
      call tell_error (tell_runtime_error, "libclim_climatology: granule duration exceeds 24 hours", errstat)
      return
   endif

   call tio_f_taix_time_to_utc_caldate(time_min, year, month, day,hour)
   bounds % hour_beg = real (hour, kind=r4)
   call tio_f_taix_time_to_utc_caldate(time_max, year, month, day,hour)
   bounds % hour_end = real (hour, kind=r4)
   bounds % lon_min = real(lon_min,kind=r4)
   bounds % lon_max = real(lon_max,kind=r4)
   bounds % lat_min = real(lat_min,kind=r4)
   bounds % lat_max = real(lat_max,kind=r4)

   !@ set bounds
   call clim_pres_init (cpt, month, day, bounds, errstat)
   if (errstat /= 0) THEN 
      call tell_error (tell_runtime_error, TRIM(ADJUSTL(modulename))//": errors in clim_pres_init", errstat)
      return
   endif
   nl0 = clim_pres_nz (cpt)

   call clim_species_init (cst, cpt, trim(clim_db_molecule_name), errstat,.true.)
   if (errstat /= 0) then
     call tell_error ( tell_io_read_error, "libclim_climatology: initializing "//trim(clim_db_molecule_name), errstat)
     return
   end if
   do_geoloc_init = .false.
   ENDIF
  !--------------------------------------
  ! get climatology for this pixel
  !---------------------------------
  lon_f = real(the_lon, kind=r4)
  lat_f = real(the_lat, kind=r4)
  call tio_f_taix_time_to_utc_caldate(the_time, year, month, day,hour)
  hour_f = real(hour, kind=r4)

  nl = nl0 ! number of level
  allocate (pres(nl))
  call clim_pres (cpt, hour_f, lon_f, lat_f, pres, errstat)
   
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, TRIM(ADJUSTL(modulename))//": errors in clim_pres", errstat)
    deallocate(pres)
    return
  end if

  allocate (vmr(nl), vmr_stddev(nl))
  ! Get vmr profile
  call clim_species_vmr (cst, cpt, hour_f, lon_f, lat_f, vmr, errstat, vmr_stddev)
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, TRIM(ADJUSTL(modulename))//": errors in clim_species_vmr", errstat)
    deallocate(pres, vmr, vmr_stddev)
    return
  end if
  
 
  ! Compute partical columns
  nl  = nl -1 ! number of layer
  allocate (partial_column(nl))
  IF (which_out == 1) THEN 
    call clim_partial_column (pres, vmr, partial_column, errstat)
  ELSE
    call clim_partial_column (pres, vmr_stddev, partial_column, errstat)
  ENDIF
  if (errstat /= 0) then
    call tell_error (tell_runtime_error, "libclim_climatology:calculating partiacl column", errstat)
    deallocate(partial_column, pres, vmr, vmr_stddev)
    return
  else
    partial_column = partial_column/du2mol
    ! Fix non-physical partial columns
    where (partial_column < 0.0_r8)
     partial_column = 0.0_r8
    end where
    deallocate(vmr, vmr_stddev)
  endif

  ! interpolation to user grids
  is_reord = .false.
  allocate(pstmp(0:nref))
  pstmp(0:nref) =ps(0:nref)
  IF (ps(0) < ps(1)) THEN
     is_reord = .true.
     pstmp(0:nref) =(/(ps(i), i = nref,0,-1)/)
   ENDIF
  CALL bspline_partial_column (nl, pres(1:nl+1)*1.D0, partial_column(1:nl)*1.D0, &
       nref, pstmp(0:nref), ozprof(1:nref), errstat)
  if (errstat /= 0) then
     call tell_error (tell_runtime_error, TRIM(ADJUSTL(modulename))//": errors in bspline_partial_column", errstat)
  else 
    IF (is_reord) THEN 
     ozprof(1:nref) = (/(ozprof(i), i = nref, 1, -1)/)
    ENDIF
  endif
  IF (errstat /= 0) THEN 
      print * , pres(1), pres(nl+1), pstmp(0), pstmp(nref); stop 1
  ENDIF
  deallocate (pstmp,pres,partial_column)
  RETURN
  END SUBROUTINE get_tempoprof
   
  SUBROUTINE bspline_partial_column(nl1,ps1,prof1, nl2,ps2,prof2, errstat)
  ! bottom to top
  IMPLICIT NONE
  ! ======================
  ! Input/Output variables
  ! ======================
  INTEGER, INTENT(IN)            :: nl1, nl2
  REAL(kind=r8), DIMENSION(nl1+1),INTENT(IN)  :: ps1
  REAL(kind=r8), DIMENSION(nl1), INTENT(IN)   :: prof1
  REAL(kind=r8), DIMENSION(nl2+1), INTENT(IN) :: ps2
  REAL(kind=r8), DIMENSION(nl2), INTENT(OUT)  :: prof2
  INTEGER, INTENT(OUT)            :: errstat
  ! ======================
  ! local variables
  ! ======================
  INTEGER :: i, fidx, lidx
  REAL (kind=r8) :: tmp
  REAL(kind=r8), ALLOCATABLE, DIMENSION(:) :: psg1, psg2, cum1, cum2
  CHARACTER(LEN=100), PARAMETER :: modulename = 'bspline_partial_column'
  errstat = 0
  IF ((ps1(1) > ps1(nl1+1) .and. ps2(1) <ps2(nl2+1)) .or. &
   (ps1(1) < ps1(nl1+1) .and. ps2(1) > ps2(nl2+1))) THEN 
    WRITE(*,*) modulename//'vertical grids are inconsistent'
    errstat = -1 ; return
  ENDIF

  fidx = 1; lidx = nl2+1
  IF (ps2(1) > ps1(1) ) THEN 
      do i = 1, nl2+1
         IF ( ps2(i) <=  ps1(1) ) exit
      enddo
      fidx = i
  ENDIF

  IF (ps2(nl2+1) < ps1(nl1+1) ) THEN 
      do i = nl2+1, 1, -1
         IF ( ps2(i) >=  ps1(nl1+1) ) exit
      enddo
      lidx = i
  ENDIF

  IF (allocated(psg1)) THEN 
     deallocate(psg1, cum1, psg2, cum2)
  ENDIF

  allocate (psg1(nl1+1), cum1(nl1+1), psg2(nl2+1), cum2(nl2+1))
  psg1 = log(ps1) ;   psg2 = log(ps2)
  cum1 = 0.D0
  DO i = 2, nl1+1
    cum1(i) = cum1(i-1) + prof1(i-1)
  ENDDO
  
  CALL BSPLINE(psg1(1:nl1+1),cum1(1:nl1+1),nl1+1,psg2(fidx:lidx),cum2(fidx:lidx),lidx-fidx + 1,errstat)
  IF (errstat /= 0 ) THEN 
      WRITE(*,*) modulename//'errors in bspline' ;stop 1
  ENDIF
  DO i = fidx-1, 1, -1
     tmp = (cum2(i+1)-cum2(i+2))/(psg2(i+1) - psg2(i+2))
     cum2(i) = cum2(i+1) + tmp*(psg2(i) - psg2(i+1))
  !   print * ,'(a)', i , tmp, cum2(i:i+2), ps2(i:i+2)
  ENDDO
  DO i = lidx+1, nl2+1
     tmp = (cum2(i-1)-cum2(i-2))/(psg2(i-1) - psg2(i-2))
     cum2(i) = cum2(i-1) + tmp*(psg2(i) - psg2(i-1))
  !   print * , '(b)',i , cum2(i+1), cum2(i), ps2(i+1), ps2(i)
  ENDDO
  prof2(1:nl2) = cum2(2:nl2+1) - cum2(1:nl2)
  deallocate(psg1, cum1, psg2, cum2)
  RETURN
  END SUBROUTINE bspline_partial_column

END MODULE m_get_o3prof
