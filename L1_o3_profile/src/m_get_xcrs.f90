! Author : jbak 
! * A part of preparing cross section dataset in lidort_prof_utilites.f90 is modified
! * O4, H2O, O2 cross section is added
! * Version of hitran LUT : 2016
! * add functionality to calculate pslwf for any highres spectrum
MODULE m_get_xcrs
   USE OMSAO_precision_module  
   USE OMSAO_parameters_module, ONLY :  maxchlen, zerok
   USE OMSAO_indices_module, ONLY: so2_idx, so2v_idx, o2o2_idx, &
       o2_idx, o2t2_idx, h2o_idx, h2ot2_idx, &
       solar_idx, wvl_idx, spc_idx, hwe_idx, spk_idx, o3_t1_idx
   USE OMSAO_variables_module,  ONLY : refdbdir, ozabs_unit,do_bandavg, &
       do_dsdw, do_dsdk, do_dsda, npsl, winwav_min,winwav_max, &
       n_refspec_pts, refspec_orig_data,  refspec_norm, database,refidx
   USE OMSAO_errstat_module
   USE ozprof_data_module,      ONLY : num_iter,ozabs_fname, &
       ozcrs_alb_fname,pos_alb, toms_fwhm, hres_slitwidth, &
       ozabs_convl, so2crs_convl,o4crs_convl,o2crs_convl, h2ocrs_convl, &
       do_so2shi, do_o4shi, do_o2shi, do_h2oshi,&
       do_so2tmp, do_o4tmp, do_o2tmp, do_h2otmp,&
       do_so2psl, do_o4psl, do_o2psl, do_h2opsl

   USE m_ezspline_interpolation, only: &
       bspline, bspline2, bspline1,interpol, interpol2, interpolation, linearinterpolationweights
   USE m_convol, ONLY: convol_f2c, convol_i0, convol_i0f2c, convol, get_i0
   USE m_avg_band, ONLY: avg_band_effozcrs
   USE m_utilities, ONLY: find_pos   
   IMPLICIT NONE

   ! Name of T-dependent tracegases cross section file  
   INTEGER, PARAMETER :: which_hitran = 2
   CHARACTER(maxchlen), PARAMETER :: so2abs_fname='OMSAO_SO2_scia_fm.dat'
!   CHARACTER(maxchlen), PARAMETER :: o4abs_fname='OMSAO_Thalman_O4quad_337-654nm.dat'
   CHARACTER(maxchlen), PARAMETER :: o4abs_fname= 'OMSAO_Thalman_O4quad_extended654nm.dat'
   !CHARACTER(maxchlen), PARAMETER :: o4abs_fname='OMSAO_Thalman_O4ts_extended654nm.dat'
   CHARACTER(maxchlen), PARAMETER :: o2abs1_fname='hitran_lut/HITRAN2016_O2_530-660nm_0p00_reduced.nc'
   CHARACTER(maxchlen), PARAMETER :: h2oabs1_fname='hitran_lut/HITRAN2016_H2O_530-660nm_0p00_reduced.nc'
   CHARACTER(maxchlen), PARAMETER :: o2abs2_fname='hitran_lut/HITRAN2016_O2_530-660nm_0p01_reduced.nc'
   CHARACTER(maxchlen), PARAMETER :: h2oabs2_fname='hitran_lut/HITRAN2016_H2O_530-660nm_0p01_reduced.nc'

   !----------------------------------------------------------------
   ! Typed cross section for O3, O4, So2, T-dependent cross section
   !------------------------------------------------------------------
   TYPE  txcrs_set
   INTEGER         :: nw, nt
   LOGICAL         :: tdepend, slitconv
   REAL (KIND=dp)  :: normc
   REAL (KIND=dp), DIMENSION(:),  ALLOCATABLE :: ts
   REAL (KIND=dp), DIMENSION(:),  ALLOCATABLE :: wvl 
   REAL (KIND=dp), DIMENSION(:,:),ALLOCATABLE :: crs  ! after convolution
   END TYPE  txcrs_set

   !-----------------------------------------------------------------
   ! T-P dependent hitran LUT
   !-----------------------------------------------------------------
   TYPE  hitran16_set
   INTEGER         :: ncid ! id of file
   INTEGER         :: wmx, bmx, pmx, tmx ! dimension
   INTEGER         :: widx0, widxf ! index range coveraing window
   REAL (KIND=dp)  :: normc
   REAL (KIND=dp), DIMENSION(:),ALLOCATABLE:: br
   REAL (KIND=dp), DIMENSION(:),ALLOCATABLE:: ps
   REAL (KIND=dp), DIMENSION(:,:),ALLOCATABLE :: ts
   REAL (KIND=dp), DIMENSION(:),  ALLOCATABLE :: wvl 
   REAL (KIND=dp), DIMENSION(:,:,:,:),ALLOCATABLE :: crs  
   END TYPE  hitran16_set

   !-----------------------------------------------------------------
   ! help variables
   !-----------------------------------------------------------------
   REAL (KIND=dp), PARAMETER :: rm_uv = 370, rm_vis = 530   

   !-----------------------------------------------------------------
   ! output variables
   !-----------------------------------------------------------------
   TYPE crsz_set ! used getgas_crs 
     LOGICAL :: do_shiwf=.false., do_tmpwf=.false., do_pslwf=.false.
     REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: crs, dads, dadt
     REAL (KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: dadp
   END TYPE crsz_set


   PUBLIC  crsz_set, geto3_crs , &! called in raman
           get_all_raycof,  &     ! called in raman
           get_alb_ozcrs_ray, &   ! nw = 1 for 347 nm
           get_hres_gascrs_ray, & ! nw > 1 without convolution
           get_effres_gascrs_ray, calc_pslwf  ! nw > 1 with convolution

   PRIVATE
   !PRIVATE read_txcrs, calc_crsz, & 
   !        read_hitran_lut,calc_hitran_crsz, & 
   !        get_all_raycof_depol,get_all_raycof_depol1
   !        geto4_crs, getso2_crs
   !        geth2o_crs_hitran16, geto2_crs_hitran16

CONTAINS 

  SUBROUTINE geto3_crs  (lamda, nlsav, nlamda, nlayers, tsgrid, crsz,problems)

  IMPLICIT NONE
  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN) :: nlamda, nlsav, nlayers
  ! #ofwave could vary for in/output
  
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)   :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers) :: tsgrid
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  TYPE (crsz_set), INTENT(INOUT) :: crsz
  LOGICAL, INTENT(OUT) :: problems
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  LOGICAL            :: dods, dodt, dodp!, do_i0convol=.true.
  INTEGER            :: i, errstat, ntemp
  REAL (KIND=dp)     :: scalex
  REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: savabs,savabs_d1 !(maxt, nlsav)
  REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: tmpabs,tmpabs_d1 !(maxt, nlamda)
  REAL (KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: dadp
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  INTEGER, SAVE :: nline, nt
  REAL (KIND=dp), DIMENSION(:,:), SAVE, ALLOCATABLE :: crs
  LOGICAL, SAVE :: first = .true.
  TYPE (txcrs_set), SAVE :: lut
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'geto3_crs'
  problems = .FALSE.
  dods = crsz%do_shiwf
  dodt = crsz%do_tmpwf
  dodp = crsz%do_pslwf
 
  ! load origianl spectrum
  IF (first) THEN
    ! read cross section 
    CALL read_txcrs (ozabs_fname,winwav_min, winwav_max,lut)  
    IF (lut%wvl(1)  > lamda(1) .OR. lut%wvl(lut%nw)  < lamda(nlsav)) THEN
       WRITE(www_lun, *) modulename//': O3abs should cover the whole fit wavelenth!!!'
       WRITE(www_lun, *)'crswav::', lut%wvl(1)  , lut%wvl(lut%nw),'omiwav:',lamda(1), lamda(nlsav)
       problems = .TRUE. ; return
    ENDIF
    IF (.NOT. lut%slitconv ) THEN
       WRITE(www_lun, *) modulename//': Need to use high-resolution cross section for O3!!!'
       problems = .TRUE. ; return
    ENDIF
    nline = lut%nw ; nt = lut%nt 
    refspec_norm (o3_t1_idx) = lut%normc
    allocate (crs(nline, nt))
    crs = lut%crs
    ozabs_convl = .true.
    first = .FALSE.
  ENDIF
  allocate (savabs(nt, nlsav), savabs_d1(nt, nlsav) )
  allocate (tmpabs(nt, nlamda),tmpabs_d1(nt, nlamda))
  !-----------------------------------------------------------------------------------------
  ! convolution     
  !-----------------------------------------------------------------------------------------
  IF (ozabs_convl .AND. lut%slitconv ) THEN  
    do_dsdw = .false.; do_dsdk = .false. ; do_dsda = .false. 
    crs(1:nline, 1:nt) = lut%crs(1:nline, 1:nt) 
    scalex = 0.2 ! ~600 DU SUM(fozs(1:nflay)) * 2.69E16 * normc, now a dummy number, not used
    ! Perform solar i0 effect on ozone cross-section (no need to convolve)
    DO i = 1, nt 
       CALL convol_i0(lut%wvl(1:nline),crs(1:nline, i), nline, scalex)
       !CALL convol(lut%wvl(1:nline),crs(1:nline, i), nline) 
    ENDDO
    ozabs_convl = .FALSE.
  ENDIF
  !------------------------------------------------------------------------
  ! interpolation
  !-------------------------------------------------------------------------
  savabs = 0.0 ; savabs_d1 = 0.0  ! onto refwvl 
  tmpabs = 0.0 ; tmpabs_d1 = 0.0  ! onto re-sampled wavelengths 

  DO i = 1, nt
    CALL BSPLINE2(lut%wvl(1:nline), crs( 1:nline, i), nline, dods, lamda, &
         savabs(i, :), savabs_d1(i, :), nlsav, errstat)

    IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename//': BSPLINE2 error, errstat = ', errstat
        problems = .TRUE.; RETURN
    ENDIF

    IF (do_bandavg) THEN
      CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)
      IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
        WRITE(www_lun, *) modulename//'O3 Spectra Averaging Error: ', nlsav, nlamda, ntemp
        problems = .TRUE.; RETURN
      ENDIF
      tmpabs(i, :) = savabs(i, 1:nlamda)

      IF (dods) THEN      
        CALL avg_band_effozcrs(lamda, savabs_d1(i, :), nlsav, ntemp, errstat)
        IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
          WRITE(www_lun, *) modulename//'O3 Spectra Averaging Error: ', nlsav, nlamda, ntemp
          problems = .TRUE.; RETURN 
        ENDIF
      tmpabs_d1(i, :) = savabs_d1(i, 1:nlamda)
      ENDIF
    ELSE
      tmpabs(i, :) = savabs(i, 1:nlamda)
      tmpabs_d1(i, :) = savabs_d1(i, 1:nlamda)
    ENDIF
  ENDDO

  IF (allocated(crsz%crs))  deallocate (crsz%crs) 
  IF (allocated(crsz%dads)) deallocate (crsz%dads) 
  IF (allocated(crsz%dadt)) deallocate (crsz%dadt) 
  IF (allocated(crsz%dadp)) deallocate (crsz%dadp) 

  allocate (crsz%crs( nlamda, nlayers)) ; crsz%crs =0.0
  crsz%crs = calc_crsz (tmpabs(1:nt,1:nlamda), &
             nt, nlamda,lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers)

  IF (dods) THEN
     allocate (crsz%dads(nlamda, nlayers))
     crsz%dads = calc_crsz( tmpabs_d1(1:nt,1:nlamda),nt,nlamda,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     crsz%dads = crsz%dads / crsz%crs   ! get relative sensitivty to shift
  ENDIF 
  IF (dodt) THEN 
     allocate (crsz%dadt(nlamda, nlayers))
     crsz%dadt = calc_tmpwf(tmpabs(1:nt,1:nlamda),nt,nlamda,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     crsz%dadt = crsz%dadt / crsz%crs   ! get relative sensitivity to T
  ENDIF
  IF (dodp) THEN 
     allocate ( crsz%dadp(nlamda, nlayers, npsl))
     allocate ( dadp(nt, nlamda, npsl))
     DO i = 1, nt 
       dadp(i, 1:nlamda, 1:npsl) =calc_pslwf ( lut%wvl, lut%crs(1:nline, i), nline, npsl,&
                                   .false., scalex, lamda, nlamda)
     ENDDO
     DO i = 1, npsl
        crsz%dadp(:,:,i)=  calc_crsz( dadp(1:nt,1:nlamda, i),nt,nlamda,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     ENDDO
     deallocate (dadp)
  ENDIF
  crsz%crs = crsz%crs*lut%normc 
  deallocate (savabs,savabs_d1)
  deallocate (tmpabs,tmpabs_d1)
  RETURN  
  END SUBROUTINE geto3_crs

  SUBROUTINE getso2_crs  (lamda, nlsav, nlamda, nlayers, tsgrid, & 
                          crsz,problems)
  IMPLICIT NONE

  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN) :: nlamda, nlsav, nlayers
       ! #ofwave could vary for in/output
  
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)   :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers) :: tsgrid
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  TYPE (crsz_set), INTENT(INOUT) :: crsz
  LOGICAL, INTENT(OUT) :: problems
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  LOGICAL            :: dods, dodt, dodp!, do_i0convol=.true.
  INTEGER            :: i, errstat, ntemp, fidx, lidx
  REAL (KIND=dp)     :: scalex
  REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: savabs,savabs_d1 !(maxt, nlsav)
  REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: tmpabs,tmpabs_d1 !(maxt, nlamda)
  REAL (KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: dadp
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  INTEGER, SAVE :: nline, nt
  REAL (KIND=dp), DIMENSION(:,:), SAVE, ALLOCATABLE :: crs
  LOGICAL, SAVE :: first = .true.
  TYPE (txcrs_set), SAVE :: lut
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'getso2_crs'
  problems = .FALSE.
  dods = crsz%do_shiwf
  dodt = crsz%do_tmpwf
  dodp = crsz%do_pslwf
 
  ! load origianl spectrum
  IF (first) THEN
    ! read cross section 
    CALL read_txcrs (ADJUSTL(TRIM(refdbdir))//so2abs_fname,winwav_min, winwav_max,lut)  
    nline = lut%nw ; nt = lut%nt 
    refspec_norm (so2_idx) = lut%normc
    allocate (crs(nline, nt))
    crs = lut%crs
    so2crs_convl = .true.
    first = .FALSE.
  ENDIF
  allocate (savabs(nt, nlsav), savabs_d1(nt, nlsav) )
  allocate (tmpabs(nt, nlamda),tmpabs_d1(nt, nlamda))
  !-----------------------------------------------------------------------------------------
  ! convolution     
  !-----------------------------------------------------------------------------------------
  IF (so2crs_convl .AND. lut%slitconv ) THEN  
    do_dsdw = .false.; do_dsdk = .false. ; do_dsda = .false. 
    crs(1:nline, 1:nt) = lut%crs(1:nline, 1:nt) 
    scalex = 0.1 ! ~600 DU SUM(fozs(1:nflay)) * 2.69E16 * normc, now a dummy number, not used
    ! Perform solar i0 effect on ozone cross-section (no need to convolve)
    DO i = 1, nt 
       CALL convol_i0(lut%wvl(1:nline),crs(1:nline, i), nline, scalex)
       !CALL convol(lut%wvl(1:nline),crs(1:nline, i), nline) 
    ENDDO
    ozabs_convl = .FALSE.
  ENDIF
   
  !------------------------------------------------------------------------
  ! interpolation
  !-------------------------------------------------------------------------
  savabs = 0.0 ; savabs_d1 = 0.0  ! onto refwvl 
  tmpabs = 0.0 ; tmpabs_d1 = 0.0  ! onto re-sampled wavelengths 

  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= lut%wvl(1))))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= lut%wvl(nline))))
  DO i = 1, nt
    CALL BSPLINE2(lut%wvl(1:nline), crs( 1:nline, i), nline, dods, lamda(fidx:lidx), &
         savabs(i, fidx:lidx), savabs_d1(i, fidx:lidx), lidx-fidx+1, errstat)

    IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename//': BSPLINE2 error, errstat = ', errstat
        problems = .TRUE.; RETURN
    ENDIF

    IF (do_bandavg) THEN
      CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)
      IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
        WRITE(www_lun, *) modulename//'so2 Spectra Averaging Error: ', nlsav, nlamda, ntemp
        problems = .TRUE.; RETURN
      ENDIF
      tmpabs(i, :) = savabs(i, 1:nlamda)

      IF (dods) THEN      
        CALL avg_band_effozcrs(lamda, savabs_d1(i, :), nlsav, ntemp, errstat)
        IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
          WRITE(www_lun, *) modulename//'so2 Spectra Averaging Error: ', nlsav, nlamda, ntemp
          problems = .TRUE.; RETURN 
        ENDIF
      tmpabs_d1(i, :) = savabs_d1(i, 1:nlamda)
      ENDIF
    ELSE
      tmpabs(i, :) = savabs(i, 1:nlamda)
      tmpabs_d1(i, :) = savabs_d1(i, 1:nlamda)
    ENDIF
  ENDDO

  IF (allocated(crsz%crs))  deallocate (crsz%crs) 
  IF (allocated(crsz%dads)) deallocate (crsz%dads) 
  IF (allocated(crsz%dadt)) deallocate (crsz%dadt) 
  IF (allocated(crsz%dadp)) deallocate (crsz%dadp) 

  allocate (crsz%crs( nlamda, nlayers)) ; crsz%crs =0.0
  crsz%crs(fidx:lidx,:) = calc_crsz (tmpabs(1:nt,fidx:lidx),nt, lidx-fidx+1, & 
                          lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers)
  IF (dods) THEN
     allocate (crsz%dads(nlamda, nlayers)) ; crsz%dads = 0.0
     crsz%dads(fidx:lidx,:) = calc_crsz( tmpabs_d1(1:nt,fidx:lidx),nt,lidx-fidx+1,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     crsz%dads(fidx:lidx,:) = crsz%dads(fidx:lidx,:) / crsz%crs(fidx:lidx,:)   ! get relative sensitivty to shift
  ENDIF 
  IF (dodt) THEN 
     allocate (crsz%dadt(nlamda, nlayers)) ; crsz%dadt = 0.0
     crsz%dadt(fidx:lidx,:) = calc_tmpwf(tmpabs(1:nt,fidx:lidx),nt,lidx-fidx+1,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     crsz%dadt(fidx:lidx,:) = crsz%dadt(fidx:lidx,:) / crsz%crs(fidx:lidx,:)   ! get relative sensitivity to T
  ENDIF
  IF (dodp) THEN 
     allocate ( crsz%dadp(nlamda, nlayers, npsl)) ; crsz%dadp = 0.0
     allocate ( dadp(nt, nlamda, npsl))
     DO i = 1, nt 
       dadp(i, fidx:lidx, 1:npsl) =calc_pslwf ( lut%wvl, lut%crs(1:nline, i), nline, npsl,&
                                   .false., scalex, lamda(fidx:lidx), lidx-fidx+1)
     ENDDO
     DO i = 1, npsl
        crsz%dadp(fidx:lidx,:,i)=  calc_crsz( dadp(1:nt,fidx:lidx, i),nt,lidx-fidx+1,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     ENDDO
     deallocate (dadp)
  ENDIF

  crsz%crs = crsz%crs*lut%normc 
  deallocate (savabs,savabs_d1)
  deallocate (tmpabs,tmpabs_d1)
  RETURN  
  END SUBROUTINE getso2_crs

  SUBROUTINE geto4_crs  (lamda, nlsav, nlamda, nlayers, tsgrid, crsz, problems)

  IMPLICIT NONE

  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN) :: nlamda, nlsav, nlayers
       ! #ofwave could vary for in/output
  
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)   :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers) :: tsgrid
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  TYPE (crsz_set), INTENT(INOUT) :: crsz
  LOGICAL, INTENT(OUT) :: problems
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  LOGICAL            :: dods, dodt, dodp!, do_i0convol=.true.
  INTEGER            :: i, errstat, ntemp, fidx, lidx
  REAL (KIND=dp)     :: scalex
  REAL (KIND=dp), DIMENSION (:,:),   ALLOCATABLE :: savabs,savabs_d1 !(maxt, nlsav)
  REAL (KIND=dp), DIMENSION (:,:),   ALLOCATABLE :: tmpabs,tmpabs_d1 !(maxt, nlamda)
  REAL (KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: dadp
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  INTEGER, SAVE :: nline, nt
  REAL (KIND=dp), DIMENSION(:,:), SAVE, ALLOCATABLE :: crs
  LOGICAL, SAVE :: first = .true.
  TYPE (txcrs_set), SAVE :: lut
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'geto4_crs'
  problems = .FALSE.
  dods = crsz%do_shiwf
  dodt = crsz%do_tmpwf
  dodp = crsz%do_pslwf
 
  ! load origianl spectrum
  IF (first) THEN
    ! read cross section 
    CALL read_txcrs (TRIM(ADJUSTL(refdbdir))//o4abs_fname,winwav_min, winwav_max,lut)  
    nline = lut%nw ; nt = lut%nt 
    refspec_norm (o2o2_idx) = lut%normc
    allocate (crs(nline, nt))
    crs = lut%crs
    o4crs_convl = .true.
    first = .FALSE. 
    WRITE(www_lun, *) dods, dodt, dodp
  ENDIF
  allocate (savabs(nt, nlsav), savabs_d1(nt, nlsav) )
  allocate (tmpabs(nt, nlamda),tmpabs_d1(nt, nlamda))
  !-----------------------------------------------------------------------------------------
  ! convolution     
  !-----------------------------------------------------------------------------------------
  IF (o4crs_convl .AND. lut%slitconv ) THEN  
    do_dsdw = .false.; do_dsdk = .false. ; do_dsda = .false. 
    crs(1:nline, 1:nt) = lut%crs(1:nline, 1:nt) 
    scalex = 0.01 ! ~600 DU SUM(fozs(1:nflay)) * 2.69E16 * normc, now a dummy number, not used
    ! Perform solar i0 effect on ozone cross-section (no need to convolve)
    DO i = 1, nt 
       CALL convol_i0(lut%wvl(1:nline),crs(1:nline, i), nline, scalex)
       !CALL convol(lut%wvl(1:nline),crs(1:nline, i), nline) 
    ENDDO
    o4crs_convl = .FALSE.
  ENDIF
  !------------------------------------------------------------------------
  ! interpolation
  !-------------------------------------------------------------------------
  savabs = 0.0 ; savabs_d1 = 0.0  ! onto refwvl 
  tmpabs = 0.0 ; tmpabs_d1 = 0.0  ! onto re-sampled wavelengths 

  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= lut%wvl(1))))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= lut%wvl(nline))))

  DO i = 1, nt 

    CALL BSPLINE2(lut%wvl(1:nline), crs( 1:nline, i), nline, dods, lamda(fidx:lidx), &
         savabs(i, fidx:lidx), savabs_d1(i, fidx:lidx), lidx-fidx+1, errstat)

    IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename//': BSPLINE2 error, errstat = ', errstat
        problems = .TRUE.; RETURN
    ENDIF

    IF (do_bandavg) THEN
      CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)
      IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
        WRITE(www_lun, *) modulename//'o4 Spectra Averaging Error: ', nlsav, nlamda, ntemp
        problems = .TRUE.; RETURN
      ENDIF
      tmpabs(i, :) = savabs(i, 1:nlamda)

      IF (dods) THEN      
        CALL avg_band_effozcrs(lamda, savabs_d1(i, :), nlsav, ntemp, errstat)
        IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
          WRITE(www_lun, *) modulename//'o4 Spectra Averaging Error: ', nlsav, nlamda, ntemp
          problems = .TRUE.; RETURN 
        ENDIF
      tmpabs_d1(i, :) = savabs_d1(i, 1:nlamda)
      ENDIF
    ELSE
      tmpabs(i, 1:nlamda) = savabs(i, 1:nlamda)
      tmpabs_d1(i, 1:nlamda) = savabs_d1(i, 1:nlamda)
    ENDIF
  ENDDO

  IF (allocated(crsz%crs))  deallocate (crsz%crs) 
  IF (allocated(crsz%dads)) deallocate (crsz%dads) 
  IF (allocated(crsz%dadt)) deallocate (crsz%dadt) 
  IF (allocated(crsz%dadp)) deallocate (crsz%dadp) 

  allocate (crsz%crs( nlamda, nlayers)) ; crsz%crs = 0.0
  crsz%crs(fidx:lidx,:) = calc_crsz (tmpabs(1:nt,fidx:lidx),nt, lidx-fidx+1, & 
                          lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers)

  IF (dods) THEN
     allocate (crsz%dads(nlamda, nlayers)) ;crsz%dads = 0.0
     crsz%dads(fidx:lidx,:) = calc_crsz( tmpabs_d1(1:nt,fidx:lidx),nt,lidx-fidx+1,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     crsz%dads(fidx:lidx,:) = crsz%dads(fidx:lidx,:) / crsz%crs(fidx:lidx,:)   ! get relative sensitivty to shift
  ENDIF 
  IF (dodt) THEN 
     allocate (crsz%dadt(nlamda, nlayers)) ; crsz%dadt = 0.0
     crsz%dadt(fidx:lidx,:) = calc_tmpwf(tmpabs(1:nt,fidx:lidx),nt,lidx-fidx+1,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     crsz%dadt(fidx:lidx,:) = crsz%dadt(fidx:lidx,:) / crsz%crs(fidx:lidx,:)   ! get relative sensitivity to T
  ENDIF
  IF (dodp) THEN 
     allocate ( crsz%dadp(nlamda, nlayers, npsl)) ; crsz%dadp = 0.0
     allocate ( dadp(nt, nlamda, npsl))
     DO i = 1, nt 
       dadp(i, fidx:lidx, 1:npsl) =calc_pslwf ( lut%wvl, lut%crs(1:nline, i), nline, npsl,&
                                   .false., scalex, lamda(fidx:lidx), lidx-fidx+1)
     ENDDO
     DO i = 1, npsl
        crsz%dadp(fidx:lidx,:,i)=  calc_crsz( dadp(1:nt,fidx:lidx, i),nt,lidx-fidx+1,  & 
        lut%tdepend,lut%ts(1:nt),tsgrid(1:nlayers), nlayers )
     ENDDO
     deallocate (dadp)
  ENDIF

  crsz%crs(fidx:lidx,:) = crsz%crs(fidx:lidx,:)*lut%normc 
  deallocate (savabs,savabs_d1)
  deallocate (tmpabs,tmpabs_d1)
  RETURN  
  END SUBROUTINE geto4_crs    

  SUBROUTINE geto2_crs_hitran16(lamda, nlsav, nlamda, nz, tsgrid, psgrid,crsz, do_convl, vs)
  USE m_convol
  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN)                                :: nlamda, nlsav, nz, vs
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)       :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nz)          :: tsgrid, psgrid
  LOGICAL, INTENT(IN) :: do_convl
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  TYPE (crsz_set), INTENT(INOUT) :: crsz
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  LOGICAL, PARAMETER :: DoSpeciesBroadening = .false.
  LOGICAL, PARAMETER :: do_i0corr=.false.
  LOGICAL            :: dods, dodt, dodp
  INTEGER            :: fidx, lidx, i, errstat, npts, rcode, vid
  INTEGER            :: idp0, idt0_p0, idt0_p1, idb0
  REAL (KIND=dp)     :: pwt0, pwt1, bwt0, bwt1,twt0_p0,twt0_p1,twt1_p0,twt1_p1, dt_p0, dt_p1
  REAL (KIND=dp)     :: scalex
  REAL (KIND=dp), DIMENSION ( :, :, :, :), ALLOCATABLE :: xs_lut
  REAL (KIND=dp), DIMENSION ( :, :, :), ALLOCATABLE :: dadp
  REAL (KIND=dp), DIMENSION ( :), ALLOCATABLE :: xs_hr, dxdt
  REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: ht_temp
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  INTEGER, SAVE :: nwvl_lut
  LOGICAL, SAVE :: first = .true.
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: crs, dxdtz
  TYPE(hitran16_set),SAVE  :: lut
  CHARACTER (LEN=256), SAVE :: filename
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  !CHARACTER (LEN=9), PARAMETER    :: modulename = 'geto2_crs'
 
  IF (first) THEN
     !---------------------------------------------------------------------
     ! read look-up table
     !----------------------------------------------------------------------
     filename = TRIM(ADJUSTL(refdbdir))//TRIM(ADJUSTL(o2abs1_fname))
     if (vs == 1) filename = TRIM(ADJUSTL(refdbdir))//TRIM(ADJUSTL(o2abs2_fname))
     CALL read_hitran16_lut(filename, winwav_min, winwav_max, lut)
     refspec_norm(o2_idx) = 1.0E-25
     lut%ps = lut%ps*0.01 ! convert to hPa
     !lut%crs = lut%crs/refspec_norm(o2_idx) 
     o2crs_convl = .TRUE.
     nwvl_lut = lut%wmx
  ENDIF
  dods = crsz%do_shiwf
  dodt = crsz%do_tmpwf
  dodp = crsz%do_pslwf
 !---------------------------------------------------------------------
 ! calculate 
 !----------------------------------------------------------------------
 !IF (num_iter == 0 ) THEN 
   allocate (xs_lut(nwvl_lut,2,2, 2))
   allocate (xs_hr(nwvl_lut), dxdt(nwvl_lut))
   allocate (ht_temp(lut%tmx))

   IF (allocated (crs)) deallocate (crs)
   IF (allocated (dxdtz)) deallocate (dxdtz)
   allocate (crs(nwvl_lut, nz))
   allocate (dxdtz(nwvl_lut, nz))
   
   lut%ncid = ncopn(trim(adjustl(filename)), nf_Nowrite, rcode)
   rcode = nf_inq_varid(lut%ncid, 'CrossSection', vid)

   DO i = 1, nz 
      ! Get pressure interpolation coefficients

      CALL LinearInterpolationWeights(lut%pmx,lut%ps, psgrid(i),idp0,pwt0,pwt1)
      ! Temperature
      ht_temp = lut%ts(:, idp0)
      CALL linearinterpolationweights(lut%tmx, ht_temp, tsgrid(i), idt0_p0, twt0_p0, twt1_p0)
      ht_temp = lut%ts(:, idp0+1)
      CALL linearinterpolationweights(lut%tmx, ht_temp, tsgrid(i), idt0_p1, twt0_p1, twt1_p1)

      ! Broaderner
      IF(DoSpeciesBroadening) THEN
          !CALL LinearInterpolationWeights(lut%bmx,lut%Broadener1Grid,MixingRatio(i, bidx1), idb0,bwt0,bwt1)
      ELSE
          idb0 = 1
          bwt0 = 1.0d0
          bwt1 = 0.0d0
      ENDIF

      ! Load the cross section LUT

      rcode = nf_get_vara_double( lut%ncid, vid,                &
                                  (/lut%widx0,idb0,idT0_p0,idp0/),&
                                  (/nwvl_lut,2,2,1/), xs_lut(:,:,:,1)      )
      rcode = nf_get_vara_double( lut%ncid, vid,                  &
                                  (/lut%widx0,idb0,idT0_p1,idp0+1/),&
                                  (/nwvl_lut,2,2,1/), xs_lut(:,:,:,2)        )

     !print *, tsgrid(i), lut%ts(idt0_p0, idp0), lut%ts(idt0_p0+1, idp0)
     ! Perform trilinear interpolation to spline grid
      xs_hr(:) =  pwt0*(bwt0*(Twt0_p0*xs_lut(:,1,1,1) + Twt1_p0*xs_lut(:,1,2,1)) +  &
                         bwt1*(Twt0_p0*xs_lut(:,2,1,1) + Twt1_p0*xs_lut(:,2,2,1))  ) &
                 + pwt1*(bwt0*(Twt0_p1*xs_lut(:,1,1,2) + Twt1_p1*xs_lut(:,1,2,2)) +  &
                         bwt1*(Twt0_p1*xs_lut(:,2,1,2) + Twt1_p1*xs_lut(:,2,2,2))  )
       
       
     ! Grid delta T
     dT_p0 = lut%ts(idT0_p0+1,idp0  )-lut%ts(idT0_p0,idp0)
     dT_p1 = lut%ts(idT0_p1+1,idp0+1)-lut%ts(idT0_p1,idp0+1)

     dxdt(:)  = pwt0*( bwt0*( xs_lut(:,1,2,1) - xs_lut(:,1,1,1) ) +        &
                bwt1*( xs_lut(:,2,2,1) - xs_lut(:,2,1,1))  )/dT_p0  &
              + pwt1*( bwt0*( xs_lut(:,1,2,2) - xs_lut(:,1,1,2)) +         &
                bwt1*( xs_lut(:,2,2,2) - xs_lut(:,2,1,2))  )/dT_p1       
     crs(1:nwvl_lut, i) = xs_hr(1:nwvl_lut) 
     dxdtz(1:nwvl_lut, i) = dxdt(1:nwvl_lut) 
    ENDDO
    call ncclos(lut%ncid, rcode)
    deallocate (xs_lut, xs_hr, dxdt, ht_temp )
  !ENDIF
  ! Convolution 

  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= lut%wvl(1) )))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= lut%wvl(nwvl_lut))))
  npts = lidx - fidx + 1

  IF (allocated(crsz%crs))  deallocate (crsz%crs) 
  IF (allocated(crsz%dads)) deallocate (crsz%dads) 
  IF (allocated(crsz%dadt)) deallocate (crsz%dadt) 
  IF (allocated(crsz%dadp)) deallocate (crsz%dadp) 

  allocate (crsz%crs( nlamda, nz)) ; crsz%crs = 0.0

  IF (fidx >= 1 .and. lidx > fidx) THEN
    IF (do_convl .and. o2crs_convl ) THEN
       scalex = 0.1 !1.0E23
       IF (do_i0corr) THEN
         CALL convol_i0f2c(lut%wvl(1:nwvl_lut), crs(1:nwvl_lut, 1:nz), nwvl_lut, nz,scalex, &
            lamda(fidx:lidx), crsz%crs(fidx:lidx, 1:nz), lidx-fidx+1)
       ELSE
         CALL convol_f2c(lut%wvl(1:nwvl_lut), crs(1:nwvl_lut, 1:nz), nwvl_lut, nz,& 
           lamda(fidx:lidx), crsz%crs(fidx:lidx, 1:nz), lidx-fidx+1)
       ENDIF
    ELSE
      DO i = 1, nz 
        CALL INTERPOL(lut%wvl(1:nwvl_lut), crs(1:nwvl_lut, i), nwvl_lut, &
             lamda(fidx:lidx), crsz%crs(fidx:lidx, i), lidx-fidx + 1, errstat)   
      ENDDO        
    ENDIF
  ENDIF

  IF (dodp) THEN 
     allocate( crsz%dadp(nlamda, nz, npsl))
     allocate ( dadp(nz, nlamda, npsl))
     DO i = 1, nz
       dadp(i, 1:nlamda, 1:npsl) =calc_pslwf ( lut%wvl, crs(1:nwvl_lut, i),nwvl_lut, npsl,&
                                   .false., scalex, lamda, nwvl_lut)
     ENDDO
     crsz%dadp  = dadp(:,:,:)
     deallocate (dadp)
  ENDIF
  !crsz%crs(fidx:lidx, 1:nz) = crsz%crs(fidx:lidx, 1:nz)*refspec_norm(o2_idx)
  IF (do_bandavg) THEN 
    WRITE(*,*) 'geto2_hitran: not implemented for do bandavg'
    stop 1
  ENDIF
  IF (first) THEN 
     first = .FALSE.
     WRITE(www_lun,*) ADJUSTL(TRIM(filename)), refspec_norm(o2_idx)
     WRITE(www_lun,*) lamda(fidx), lamda(lidx), fidx, lidx
     WRITE(www_lun,*) 'do_conv/do_io',do_convl, do_i0corr
  ENDIF
  RETURN  
  END SUBROUTINE geto2_crs_hitran16

  SUBROUTINE geth2o_crs_hitran16(lamda, nlsav, nlamda, nz, tsgrid, psgrid, crsz, do_convl, vs)
  USE m_convol
  USE OMSAO_variables_module, ONLY: npsl
  USE ozprof_data_module, ONLY: mgasprof, h2oidx, h2ot2idx
  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'
  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN)                                :: nlamda, nlsav, nz, vs
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)       :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nz)          :: tsgrid, psgrid
  LOGICAL, INTENT(IN) :: do_convl
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  TYPE (crsz_set), INTENT(INOUT) :: crsz
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  LOGICAL, PARAMETER :: DoSpeciesBroadening = .false.
  LOGICAL, PARAMETER :: do_i0corr=.false.
  LOGICAL            :: dods, dodt, dodp
  INTEGER            :: fidx, lidx, i,  errstat, npts, rcode, vid
  INTEGER            :: idp0, idt0_p0, idt0_p1, idb0
  REAL (KIND=dp)     :: pwt0, pwt1, bwt0, bwt1,twt0_p0,twt0_p1,twt1_p0,twt1_p1, dt_p0, dt_p1
  REAL (KIND=dp)     :: scalex
  REAL (KIND=dp), DIMENSION(:),ALLOCATABLE :: h2oprof
  REAL (KIND=dp), DIMENSION ( :, :, :, :), ALLOCATABLE :: xs_lut
  REAL (KIND=dp), DIMENSION ( :, :, :), ALLOCATABLE :: dadp
  REAL (KIND=dp), DIMENSION ( :), ALLOCATABLE :: xs_hr, dxdt
  REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: ht_temp
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  INTEGER, SAVE :: nwvl_lut, z1, z2
  LOGICAL, SAVE :: first = .true.
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: crs, dxdtz
  TYPE(hitran16_set), save  :: lut
  CHARACTER (LEN=256), SAVE :: filename
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  !CHARACTER (LEN=10), PARAMETER    :: modulename = 'geth2o_crs'
 
  IF (first) THEN
     !---------------------------------------------------------------------
     ! read look-up table
     !----------------------------------------------------------------------
     filename = TRIM(ADJUSTL(refdbdir))//TRIM(ADJUSTL(h2oabs1_fname))
     if (vs == 1) filename = TRIM(ADJUSTL(refdbdir))//TRIM(ADJUSTL(h2oabs2_fname))
     CALL read_hitran16_lut(filename, winwav_min, winwav_max, lut)
     lut%ps = lut%ps*0.01 ! convert to hPa
     refspec_norm(h2o_idx) = 1.0E-25
     h2ocrs_convl = .TRUE.
     nwvl_lut = lut%wmx
  ENDIF
  dods = crsz%do_shiwf
  dodt = crsz%do_tmpwf
  dodp = crsz%do_pslwf
 !---------------------------------------------------------------------
 ! calculate 
 !----------------------------------------------------------------------
 IF (num_iter == 0 ) THEN 
   allocate (xs_lut(nwvl_lut,2,2, 2))
   allocate (xs_hr(nwvl_lut), dxdt(nwvl_lut))
   allocate (ht_temp(lut%tmx))

   IF (allocated (crs)) deallocate (crs)
   IF (allocated (dxdtz)) deallocate (dxdtz)
   allocate (crs(nwvl_lut, nz))
   allocate (dxdtz(nwvl_lut, nz))
   
   lut%ncid = ncopn(trim(adjustl(filename)), nf_Nowrite, rcode)
   rcode = nf_inq_varid(lut%ncid, 'CrossSection', vid)

   allocate (h2oprof(nz))
   h2oprof(1:nz) = mgasprof(h2oidx, 1:nz) + mgasprof(h2ot2idx, 1:nz)
   DO i =  1, nz
     IF (h2oprof(i) > 0 ) EXIT
   ENDDO
   deallocate (h2oprof)
   z1 = i
   z2 = nz
   DO i = z1, z2 
      ! Get pressure interpolation coefficients

      CALL LinearInterpolationWeights(lut%pmx,lut%ps, psgrid(i),idp0,pwt0,pwt1)

      ! Temperature
      ht_temp = lut%ts(:, idp0)
      CALL linearinterpolationweights(lut%tmx, ht_temp, tsgrid(i), idt0_p0, twt0_p0, twt1_p0)
      ht_temp = lut%ts(:, idp0+1)
      CALL linearinterpolationweights(lut%tmx, ht_temp, tsgrid(i), idt0_p1, twt0_p1, twt1_p1)

      ! Broaderner
      IF(DoSpeciesBroadening) THEN
          !CALL LinearInterpolationWeights(lut%bmx,lut%Broadener1Grid,MixingRatio(i, bidx1), idb0,bwt0,bwt1)
      ELSE
          idb0 = 1
          bwt0 = 1.0d0
          bwt1 = 0.0d0
      ENDIF

      ! Load the cross section LUT
      rcode = nf_get_vara_double( lut%ncid, vid,                &
                                  (/lut%widx0,idb0,idT0_p0,idp0/),&
                                  (/nwvl_lut,2,2,1/), xs_lut(:,:,:,1)      )
      rcode = nf_get_vara_double( lut%ncid, vid,                  &
                                  (/lut%widx0,idb0,idT0_p1,idp0+1/),&
                                  (/nwvl_lut,2,2,1/), xs_lut(:,:,:,2)        )

     !print *, tsgrid(i), lut%ts(idt0_p0, idp0), lut%ts(idt0_p0+1, idp0)
     ! Perform trilinear interpolation to spline grid
      xs_hr(:) =  pwt0*(bwt0*(Twt0_p0*xs_lut(:,1,1,1) + Twt1_p0*xs_lut(:,1,2,1)) +  &
                         bwt1*(Twt0_p0*xs_lut(:,2,1,1) + Twt1_p0*xs_lut(:,2,2,1))  ) &
                 + pwt1*(bwt0*(Twt0_p1*xs_lut(:,1,1,2) + Twt1_p1*xs_lut(:,1,2,2)) +  &
                         bwt1*(Twt0_p1*xs_lut(:,2,1,2) + Twt1_p1*xs_lut(:,2,2,2))  )
       
       
     ! Grid delta T
     dT_p0 = lut%ts(idT0_p0+1,idp0  )-lut%ts(idT0_p0,idp0)
     dT_p1 = lut%ts(idT0_p1+1,idp0+1)-lut%ts(idT0_p1,idp0+1)

     dxdt(:)  = pwt0*( bwt0*( xs_lut(:,1,2,1) - xs_lut(:,1,1,1) ) +        &
                bwt1*( xs_lut(:,2,2,1) - xs_lut(:,2,1,1))  )/dT_p0  &
              + pwt1*( bwt0*( xs_lut(:,1,2,2) - xs_lut(:,1,1,2)) +         &
                bwt1*( xs_lut(:,2,2,2) - xs_lut(:,2,1,2))  )/dT_p1       
     crs(1:nwvl_lut, i) = xs_hr(1:nwvl_lut) 
     dxdtz(1:nwvl_lut, i) = dxdt(1:nwvl_lut) 
    ENDDO
    call ncclos(lut%ncid, rcode)
    deallocate (xs_lut, xs_hr, dxdt, ht_temp )
  ENDIF

  IF (allocated(crsz%crs))  deallocate (crsz%crs) 
  IF (allocated(crsz%dads)) deallocate (crsz%dads) 
  IF (allocated(crsz%dadt)) deallocate (crsz%dadt) 
  IF (allocated(crsz%dadp)) deallocate (crsz%dadp) 

  allocate (crsz%crs( nlamda, nz)) ; crsz%crs = 0.0

  ! Convolution 
  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= lut%wvl(1) )))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= lut%wvl(nwvl_lut))))
  npts = lidx - fidx + 1
  crsz%crs = 0.0
  IF (fidx >= 1 .and. lidx > fidx) THEN
    !PRINT * , 'geth2o_jbak', do_convl, do_i0corr
    IF (do_convl .and. h2ocrs_convl ) THEN
       scalex = 0.1 !1.0E23
       IF (do_i0corr) THEN
         CALL convol_i0f2c(lut%wvl(1:nwvl_lut), crs(1:nwvl_lut, z1:z2), nwvl_lut, z2-z1+1,scalex, &
            lamda(fidx:lidx), crsz%crs(fidx:lidx, z1:z2), lidx-fidx+1)
       ELSE
         CALL convol_f2c(lut%wvl(1:nwvl_lut), crs(1:nwvl_lut, z1:z2), nwvl_lut, z2-z1+1,& 
           lamda(fidx:lidx), crsz%crs(fidx:lidx, z1:z2), lidx-fidx+1)
       ENDIF
    ELSE
      DO i = z1, z2
        CALL INTERPOL(lut%wvl(1:nwvl_lut), crs(1:nwvl_lut, i), nwvl_lut, &
             lamda(fidx:lidx), crsz%crs(fidx:lidx, i), lidx-fidx + 1, errstat)  
        !CALL INTERPOL(lut%wvl(1:nwvl_lut), dxdtz(1:nwvl_lut, i), nwvl_lut, &
        !     lamda(fidx:lidx), crsz%crs(fidx:lidx, i), lidx-fidx + 1, errstat)    
      ENDDO        
    ENDIF
  ENDIF
  IF (dodp) THEN 
     allocate( crsz%dadp(nlamda, nz, npsl))
     allocate ( dadp(nz, nlamda, npsl))
     dadp(:,:,:) = 0.0
     DO i = z1, z2
       dadp(i, 1:nlamda, 1:npsl) =calc_pslwf ( lut%wvl, crs(1:nwvl_lut, i),nwvl_lut, npsl,&
                                   .false., scalex, lamda, nwvl_lut)
     ENDDO
     crsz%dadp  = dadp(:,:,:)
     deallocate (dadp)
  ENDIF
  crsz%crs(fidx:lidx, 1:nz) = crsz%crs(fidx:lidx, 1:nz) !*refspec_norm(h2o_idx)
  IF (do_bandavg) THEN 
    WRITE(*,*) 'geth2o_hitran: not implemented for do bandavg'
    stop 1
  ENDIF
  IF (first) THEN 
     first = .FALSE.
     WRITE(www_lun,*) ADJUSTL(TRIM(filename)), refspec_norm(h2o_idx)
     WRITE(www_lun,*) lamda(fidx), lamda(lidx), fidx, lidx
     WRITE(www_lun,*) 'do_conv/do_io',do_convl, do_i0corr
  ENDIF
  RETURN  
  END SUBROUTINE geth2o_crs_hitran16

  SUBROUTINE get_all_raycof_depol(nw, waves, raycof, depol)
  IMPLICIT NONE
  !------------------
  !Input/Output
  !------------------
  INTEGER, INTENT(IN)                        :: nw
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(:),  INTENT(OUT), ALLOCATABLE :: raycof, depol
  !-----------------
  !Local variables
  !-----------------
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: & 
       sig, sig2, sig2p, sig4, fk_n2, fk_o2, fking ! nw
  REAL (KIND=dp), PARAMETER     :: abod = 1.0455996d0, bbod = -341.29061d0, &
       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
  !--------------------------------------------------------------------
  !     Rayleigh coefficient
  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
  !--------------------------------------------------------------------
  IF (allocated (raycof)) deallocate (raycof, depol)
  allocate (raycof(nw), depol(nw))

  allocate ( sig(nw), sig2(nw), sig2p(nw), sig4(nw), fk_n2(nw),fk_o2(nw), fking(nw))
  sig =    1.0d3 / waves
  sig2 =   sig * sig
  sig2p =  1.d0 / sig2
  sig4 =   sig2 * sig2
  raycof = (abod + bbod * sig2 + cbod * sig2p) &
       / (1.d0 + dbod * sig2 + ebod * sig2p) * 1.d-28

  !     Derivation of depolarization factor d from king factors for air,
  !     fking = (6 + 3.depol) / (6 - 7.depol)
  !     bodhaine et al., 370 ppmv co2
  fk_n2 = 1.034d0 + 3.17d-4 * sig2
  fk_o2 = 1.096d0 + 1.385d-3 * sig2 + 1.448d-4 * sig4
  fking = (78.084d0 * fk_n2 + 20.946d0 * fk_o2 + 0.97655d0) / 100.001d0
  depol = 6.d0 * (fking - 1.d0) / (3.d0 + 7.d0 * fking)

  deallocate ( sig, sig2, sig2p, sig4, fk_n2,fk_o2, fking)

  RETURN

  END SUBROUTINE get_all_raycof_depol

  SUBROUTINE get_all_raycof(nw, waves, raycof)

  IMPLICIT NONE
  !-------------------
  !Input/Output
  !------------------
  INTEGER, INTENT(IN)                        :: nw
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(nw), INTENT(OUT) :: raycof
  !------------------
  !Local variables
  !------------------
  REAL (KIND=dp), PARAMETER     :: abod = 1.0455996d0, bbod = -341.29061d0, &
       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: sig, sig2, sig2p, sig4

  !---------------------------
  ! Rayleigh coefficient
  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
  !--------------------------
  
  allocate (sig(nw), sig2(nw), sig2p(nw), sig4(nw))
  sig =    1.0d3 / waves
  sig2 =   sig * sig
  sig2p =  1.d0 / sig2
  sig4 =   sig2 * sig2
  raycof = (abod + bbod * sig2 + cbod * sig2p) &
       / (1.d0 + dbod * sig2 + ebod * sig2p) * 1.d-28
  deallocate (sig, sig2, sig2p, sig4)
  RETURN

  END SUBROUTINE get_all_raycof 

  SUBROUTINE get_all_raycof_depol1 (nw, waves, nw1, raycof, depol, problems)
  ! it is called for every pixels, so allocation/deaollocation is not applied.
  IMPLICIT NONE

  !     Input/Output
  INTEGER, INTENT(IN)                        :: nw, nw1
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  !REAL (KIND=dp), DIMENSION(nw1), INTENT(OUT):: raycof, depol
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT):: raycof, depol
  LOGICAL, INTENT(OUT)                       :: problems

  ! Local variables
  INTEGER                               :: errstat, ntemp
  REAL (KIND=dp)                        :: scalex
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE      :: raycof1, depol1

  INTEGER, SAVE                                  :: nref
  REAL (KIND=dp), DIMENSION(:),ALLOCATABLE, SAVE :: ray, dep, refwavs
  REAL (KIND=dp),                           SAVE :: rnorm, dnorm
  LOGICAL, SAVE                                  :: first = .TRUE.

  problems = .FALSE.
  IF (first) THEN
     nref = n_refspec_pts(1)
     allocate (refwavs(nref))
     refwavs(1:nref) = refspec_orig_data(1, 1:nref, 1)

     CALL get_all_raycof_depol(nref, refwavs, ray, dep)
     rnorm = 1.0E-25; dnorm = 1.0E-2
     ray(1:nref) = ray(1:nref) / rnorm; dep(1:nref) = dep(1:nref) / dnorm

     scalex = 1.0  ! dummy variable here
     CALL convol_i0(refwavs(1:nref), ray(1:nref), nref, scalex)
     CALL convol_i0(refwavs(1:nref), dep(1:nref), nref, scalex)
     first = .FALSE.
  ENDIF 
  allocate ( raycof1(nw), depol1(nw))
  ! interpolation
  CALL BSPLINE(refwavs(1:nref), ray(1:nref), nref, waves, raycof1, nw, errstat)
  IF (errstat < 0) THEN
     WRITE(www_lun, *) 'Error in interpolating raycof!!!'
     problems = .TRUE.; RETURN
  ENDIF

  CALL BSPLINE(refwavs(1:nref), dep(1:nref), nref, waves, depol1,  nw, errstat)
  IF (errstat < 0) THEN
     WRITE(www_lun, *) 'Error in interpolating depol!!!'
     problems = .TRUE.; RETURN
  ENDIF

  IF (do_bandavg) THEN
     CALL avg_band_effozcrs(waves, raycof1, nw, ntemp, errstat)
     IF ( errstat /= 0 .OR. ntemp /= nw1) THEN
        WRITE(www_lun, *) 'Raycof Spectra Averaging Error: ', nw, nw1, ntemp
        problems = .TRUE.; RETURN
     ENDIF

     CALL avg_band_effozcrs(waves, depol1, nw, ntemp, errstat)
     IF ( errstat /= 0 .OR. ntemp /= nw1) THEN
        WRITE(www_lun, *) 'Depol Spectra Averaging Error: ', nw, nw1, ntemp
        problems = .TRUE.; RETURN
     ENDIF
  ENDIF

  IF (allocated (raycof))  deallocate(raycof, depol) 
  allocate (raycof(nw1), depol(nw1))
  raycof = raycof1(1:nw1) * rnorm
  depol  = depol1(1:nw1)  * dnorm
  deallocate (raycof1, depol1)
  RETURN
  END SUBROUTINE get_all_raycof_depol1

  SUBROUTINE get_alb_ozcrs_ray (nz, ts, ngas, abscrs, raycof, depol, problems)

  IMPLICIT NONE
  
  ! Input variables
  INTEGER, INTENT(IN)                              :: nz, ngas
  REAL (KIND=dp), INTENT(IN),  DIMENSION(nz)       :: ts
  REAL (KIND=dp), INTENT(OUT)                      :: raycof, depol
  REAL (KIND=dp), INTENT(OUT), DIMENSION(ngas, nz) :: abscrs
  LOGICAL, INTENT(OUT)                             :: problems
  
  ! Local variable
  INTEGER                                          :: nline, nw, i
  REAL (KIND=dp)                                   :: fwav, lwav, swav, ewav
  REAL (KIND=dp), DIMENSION(11)                    :: temp
  REAL (KIND=dp), DIMENSION(:),  ALLOCATABLE  :: waves, sol, weights, rays, dpols 
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: ozcrs  
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: gcrs   
  CHARACTER (len=maxchlen)                   :: crs_fname
  
  ! Saved variables
  LOGICAL,                      SAVE :: first = .TRUE.
  REAL (KIND=dp), DIMENSION(3), SAVE :: cozcrs
  REAL (KIND=dp), DIMENSION(6), SAVE :: cgcrs
  REAL (KIND=dp),               SAVE :: craycof, cdepol

  problems = .FALSE.
  IF (first) THEN
     crs_fname = TRIM(ADJUSTL(refdbdir)) // '/' // TRIM(ADJUSTL(ozcrs_alb_fname))
     OPEN(UNIT = ozabs_unit, file=crs_fname, status='old')     
     DO i = 1, 5
        READ(ozabs_unit, *)
     ENDDO
     READ(ozabs_unit, *) nline, fwav, lwav
     READ(ozabs_unit, *)
     swav = pos_alb - toms_fwhm;  ewav = pos_alb + toms_fwhm
     IF (swav < fwav .OR. ewav > lwav) THEN
        WRITE(*, *) 'Ozone cross section does not cover region for determining fcld!!!'
        problems = .TRUE.; CLOSE(ozabs_unit); RETURN
     ENDIF

     allocate (waves(nline), sol(nline), weights(nline))
     allocate (ozcrs(3, nline), gcrs(6, nline))

     nw = 0
     DO i = 1, nline
        READ (ozabs_unit, *) temp
        IF (temp(1) >= swav .AND. temp(1) <= ewav) THEN
           nw = nw + 1; waves(nw) = temp(1)
           ozcrs(1:3, nw) = temp(2:4); gcrs(1:6, nw) = temp(5:10); sol(nw) = temp(11)
        ELSE IF (temp(1) > ewav) THEN
           EXIT
        ENDIF
     ENDDO
     CLOSE (ozabs_unit)

     ! Compute weights
     weights(1:nw) = (1.0 - ABS(waves(1:nw) - pos_alb) / toms_fwhm) * sol(1:nw)
     weights(1:nw) = weights(1:nw) /  SUM(weights(1:nw))

     DO i = 1, 3
        cozcrs(i) = SUM(ozcrs(i, 1:nw) *  weights(1:nw))
     ENDDO

     DO i = 1, 6
        cgcrs(i) = SUM(gcrs(i, 1:nw) *  weights(1:nw))
     ENDDO

     CALL GET_ALL_RAYCOF_DEPOL(nw, waves(1:nw), rays, dpols)
     craycof = SUM(rays(1:nw)   * weights(1:nw)) 
     cdepol  = SUM(dpols(1:nw)  * weights(1:nw)) 

     first = .FALSE.
     deallocate (waves, sol, weights, ozcrs, gcrs)
  ENDIF
     
  abscrs(1, :) = cozcrs(1) + (ts - zerok) * cozcrs(2) + (ts - zerok) ** 2.0 * cozcrs(3)
  DO i = 1, 6
     abscrs(i+1, :) = cgcrs(i)
  ENDDO
  ! O3, NO2, SO2, BrO, HCO, O4, OCLO, F0
  raycof = craycof; depol = cdepol

  RETURN
  END SUBROUTINE get_alb_ozcrs_ray

  SUBROUTINE get_effres_gascrs_ray (& 
   num_iter, nlsav, lamda, nlamda, nz, ts, ps,nfgas, allcol, rhos, &
   do_o3shi, o3shi,do_tmpwf, do_pslwf, allcrs, raycof, depol,problems)

  USE OMSAO_indices_module, ONLY: so2_idx, so2v_idx, o2o2_idx, o2_idx, h2o_idx
  USE OMSAO_variables_module, ONLY: numwin, winlim, & 
      database_save, database_shiwf,  &
      n_refspec_pts,  refspec_orig_data, refidx, &
      fitvar_rad, rmask_fitvar_rad, database, npsl
  USE ozprof_data_module, ONLY: do_subfit, nos, oswins, osfind,  &
      use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, & 
      ngas, gasidxs,fgasidxs, fgassidxs, ccrs, dadp, dads, dadt

  IMPLICIT NONE
  ! Input/output variables
  INTEGER, INTENT (IN)                                  :: nlsav,nlamda, nz, nfgas, num_iter
  LOGICAL, INTENT (IN)                                  :: do_o3shi, do_tmpwf,do_pslwf
  REAL (KIND=dp), DIMENSION(nlsav), INTENT (IN )        :: lamda
  REAL (KIND=dp), DIMENSION(nz), INTENT (IN )           :: ts, ps, rhos
  REAL (KIND=dp), DIMENSION(nfgas, nz), INTENT(IN)          :: allcol
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)        :: o3shi
  REAL (KIND=dp), DIMENSION(nlamda), INTENT (OUT)           :: raycof, depol
  REAL (KIND=dp), DIMENSION(nlamda, nfgas, nz), INTENT(OUT) :: allcrs
  LOGICAL, INTENT (OUT)                                     :: problems
  ! Local variables
  LOGICAL :: do_convl
  INTEGER :: fidx, lidx, i, j,k, npts, nfgas1, errstat
  REAL (kind=dp) :: tmp, temp, normc
  REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: delshi, delpos, gshiwf !(nlsav)
  REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: refspec, refwav
  !--------------------------------
  ! Save variables (max_fit_pts)
  !--------------------------------
  TYPE (crsz_set), SAVE :: o3, so2, o4, o2, h2o
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE,SAVE :: raycof0, depol0
  !-------------------------------------------------
  ! Name of this subroutine/module 
  !-------------------------------------------------
  CHARACTER(10), PARAMETER :: modulename ='get_effcrs'

  allcrs = 0.0

  IF (num_iter == 0) THEN
     CALL get_all_raycof_depol1(nlsav, lamda, nlamda,& 
          raycof0,depol0,problems)
        IF (allocated (ccrs%o3)) deallocate(ccrs%o3)
        allocate (ccrs%o3(nlamda, nz))
     IF (use_so2dtcrs) then
        IF (allocated (ccrs%so2)) deallocate(ccrs%so2)
        allocate (ccrs%so2(nlamda, nz))
     ENDIF
     IF (use_o4dtcrs)  then
        IF (allocated (ccrs%o4)) deallocate(ccrs%o4)
        allocate (ccrs%o4(nlamda, nz)) 
     ENDIF
     IF (use_h2odptcrs) then 
        IF (allocated (ccrs%h2o)) deallocate(ccrs%h2o)
        allocate (ccrs%h2o(nlamda, nz)) 
     ENDIF
     IF (use_o2dptcrs) then 
        IF (allocated (ccrs%o2)) deallocate(ccrs%o2)
        allocate (ccrs%o2(nlamda, nz))
     ENDIF
     IF (do_o3shi) THEN 
         IF (allocated(dads%o3)) deallocate(dads%o3)
         allocate (dads%o3(nlamda, nz))
     ENDIF
     IF (do_pslwf) THEN 
         IF (allocated(dadp%o3)) deallocate(dadp%o3)
         allocate (dadp%o3(nlamda, nz, npsl))
         IF (use_so2dtcrs .and. do_so2psl) THEN 
            IF (allocated(dadp%so2)) deallocate(dadp%so2)
            allocate (dadp%so2(nlamda, nz, npsl))
         ENDIF
         IF (use_o4dtcrs .and. do_o4psl) THEN 
            IF (allocated(dadp%o4)) deallocate(dadp%o4)
            allocate (dadp%o4(nlamda, nz, npsl))
         ENDIF
         IF (use_h2odptcrs .and. do_h2opsl) THEN 
             IF (allocated(dadp%h2o)) deallocate(dadp%h2o)
            allocate (dadp%h2o(nlamda, nz, npsl))
         ENDIF
         IF (use_o2dptcrs .and. do_o2psl) THEN 
            IF (allocated(dadp%o2)) deallocate(dadp%o2)
            allocate (dadp%o2(nlamda, nz, npsl))
         ENDIF
     ENDIF
     IF (do_tmpwf) THEN 
         IF (allocated(dadt%o3)) deallocate(dadt%o3)
         allocate (dadt%o3(nlamda, nz))
         IF (use_so2dtcrs .and. do_so2tmp) THEN 
            IF (allocated(dadt%so2)) deallocate(dadt%so2)
            allocate (dadt%so2(nlamda, nz))
         ENDIF
         IF (use_o4dtcrs .and. do_o4tmp) THEN 
            IF (allocated(dadt%o4)) deallocate(dadt%o4)
            allocate (dadt%o4(nlamda, nz))
         ENDIF
         IF (use_h2odptcrs .and. do_h2otmp) THEN 
             IF (allocated(dadt%h2o)) deallocate(dadt%h2o)
            allocate (dadt%h2o(nlamda, nz))
         ENDIF
         IF (use_o2dptcrs .and. do_o2tmp) THEN 
            IF (allocated(dadt%o2)) deallocate(dadt%o2)
            allocate (dadt%o2(nlamda, nz))
         ENDIF
     ENDIF
  ENDIF 

  raycof(1:nlamda) = raycof0(1:nlamda) 
  depol(1:nlamda) =  depol0(1:nlamda)

  allocate(delshi(nlsav), delpos(nlsav), gshiwf(nlsav))
  delshi = 0.0
  IF (do_o3shi .AND. nos > 0 ) THEN 
    IF (do_subfit) THEN 
      fidx = 1
      DO j = 1, numwin
        IF (j == numwin) THEN 
           lidx = nlsav
        ELSE
           tmp = (winlim(j, 2) + winlim(j+1, 1)) / 2.0
           lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <=tmp)))
        ENDIF
        delpos(fidx:lidx) =  lamda(fidx:lidx) - (lamda(fidx) +lamda(lidx)) / 2.0
        IF (osfind(j, 1) > 0) delshi(fidx:lidx) =  o3shi(j, 1)
          DO i = 2, nos
            IF (osfind(j, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx) + &
                         o3shi(j, i) * delpos(fidx:lidx) ** (i-1)
          ENDDO
          fidx = lidx + 1
      ENDDO    
    ELSE
      IF (oswins(1, 1) == 1) THEN
        fidx = 1
      ELSE
        tmp = (winlim(oswins(1, 1), 1) + winlim(oswins(1, 1) - 1, 2))/2.
        fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >=tmp)))
      ENDIF
      IF (oswins(1, 2) == numwin) THEN
         lidx = nlsav
      ELSE
         tmp = (winlim(oswins(1, 2), 2) + winlim(oswins(1, 2) + 1, 1))/2.
         lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <=tmp)))
      ENDIF
      delpos(fidx:lidx) =  lamda(fidx:lidx) - (lamda(fidx) +lamda(lidx)) / 2.0
      IF (osfind(1, 1) > 0) delshi(fidx:lidx) =  + o3shi(1, 1)
      DO i = 2, nos
          IF (osfind(1, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx) + &
             o3shi(1, i) * delpos(fidx:lidx) ** (i-1)
      ENDDO
    ENDIF
  ENDIF
  
  IF (num_iter == 0. .OR.  do_o3shi .OR. do_pslwf) THEN 
     !IF (num_iter == 0) WRITE(*, *) modulename,' : Set up O3 absorption !!!'
      o3%do_shiwf = do_o3shi ; o3%do_tmpwf = do_tmpwf ; o3%do_pslwf = do_pslwf
      CALL geto3_crs(lamda - delshi, nlsav, nlamda, nz, ts, o3, problems)
      IF (problems) THEN
        WRITE(www_lun, *) modulename,' : Problems in reading O3 absorption !!!'
        RETURN
      ENDIF
  ENDIF
  ! GET T-dependent So2 cross section at instrument spectral resolution
  IF ( num_iter == 0 .AND.  use_so2dtcrs) THEN
    so2%do_shiwf = do_so2shi ; so2%do_tmpwf = do_so2tmp ; so2%do_pslwf = do_so2psl
    !WRITE(*, *) modulename,' : Set up SO2 absorption !!!'
    CALL getso2_crs(lamda, nlsav, nlamda, nz, ts, so2, problems)  
    IF (problems) THEN
       WRITE(*, *) modulename, ' : Problems in reading SO2 absorption !!!'
      RETURN
    ENDIF
  ENDIF
 
  ! GET T-dependent O4 cross section at instrument spectral resolution
  IF (num_iter == 0 .AND. use_o4dtcrs) THEN
    !WRITE(*, *) modulename,' : Set up O4 absorption !!!'
    o4%do_shiwf = do_o4shi ; o4%do_tmpwf = do_o4tmp ; o4%do_pslwf = do_o4psl
   ! CALL geto4_crs_old(lamda, nlsav, nlamda,  nz, ts, ccrs%o4, problems)
    CALL geto4_crs(lamda, nlsav, nlamda,  nz, ts, o4, problems)
    IF (problems) THEN
       WRITE(*, *) modulename, ' : Problems in reading O4 absorption !!!'
      RETURN
    ENDIF
  ENDIF
  
  do_convl = .true.
  ! GET T-dependent o2 cross section at instrument spectral resolution
  IF ( num_iter == 0 .AND. use_o2dptcrs) THEN
    !WRITE(*, *) modulename,' : Set up O2 absorption !!!'
    o2%do_shiwf = do_o2shi ; o2%do_tmpwf = do_o2tmp ; o2%do_pslwf = do_o2psl
     !CALL GET_O2_CiRS(lamda, nlsav, nlamda, nz, ps(1:nz), ts(1:nz),& 
     ! crsz%o2(1:nlamda, 1:nz), problems)
    CALL geto2_crs_hitran16(lamda, nlsav, nlamda, nz,ts, ps,o2, do_convl,0)    
    IF (problems) THEN
      WRITE(*, *) modulename, ' : Problems in getting O2 cross section!!!'
      RETURN
    ENDIF
  ENDIF
  
  ! GET T-dependent h2o cross section at instrument spectral resolution
  IF ( num_iter == 0 .AND. use_h2odptcrs) THEN
    !WRITE(*, *) modulename,' : Set up H2O absorption !!!'
     !h2oprof(1:nz) = mgasprof(h2oidx, 1:nz) + mgasprof(h2ot2idx, 1:nz)
     !DO k =  1, nz
     !   IF (h2oprof(k) > 0 ) EXIT
     !ENDDO
     !fidx = k
     !lidx = nz
     !CALL GET_H2O_CRS(lamda, nlsav, nlamda, nz, ps(1:nz)/1013.25, ts(1:nz), &
     !     h2oprof(1:nz), ccrs%h2o(1:nlamda, 1:nz), problems) ! XLIU
     h2o%do_shiwf = do_h2oshi ; h2o%do_tmpwf = do_h2otmp ; h2o%do_pslwf = do_h2opsl
     CALL geth2o_crs_hitran16(lamda, nlsav, nlamda, nz, ts, ps, h2o, do_convl,0)
    IF (problems) THEN
      WRITE(*, *) modulename, ' : Problems in getting H2O cross section!!!'
      RETURN
    ENDIF
  ENDIF

  allcrs(1:nlamda, 1, 1:nz) = o3%crs(1:nlamda, 1:nz)
  IF (do_o3shi) dads%o3(1:nlamda, 1:nz) = o3%dads(1:nlamda, 1:nz)
  IF (do_tmpwf) dadt%o3(1:nlamda, 1:nz) = o3%dadt(1:nlamda, 1:nz)
  
  nfgas1 = 1
  ccrs%o3 = o3%crs(1:nlamda, 1:nz)
  DO i = 1, ngas
     IF (fgasidxs(i) > 0 ) THEN 
         nfgas1 = nfgas1 + 1
         normc = refspec_norm(gasidxs(i))
         IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. use_so2dtcrs) THEN 
           ccrs%so2 = so2%crs
           allcrs(1:nlamda, nfgas1, 1:nz) = ccrs%so2(1:nlamda,1:nz)/normc
           IF (do_so2psl) dadp%so2(:,:,:) = so2%dadp 
         ELSE IF (gasidxs(i) == o2o2_idx .AND. use_o4dtcrs) THEN 
           ccrs%o4 = o4%crs
           allcrs(1:nlamda, nfgas1, 1:nz) = ccrs%o4(1:nlamda,1:nz)/normc
           IF (do_o4psl) dadp%o4(:,:,:) = o4%dadp 
         ELSE IF ((gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx) .AND. use_o2dptcrs) THEN 
           ccrs%o2 = o2%crs
           allcrs(1:nlamda, nfgas1, 1:nz) = ccrs%o2(1:nlamda,1:nz)/normc
           IF (do_o2psl) dadp%o2(:,:,:) = o2%dadp 
         ELSE IF ((gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx) .AND. use_h2odptcrs) THEN 
           ccrs%h2o = h2o%crs
           allcrs(1:nlamda, nfgas1, 1:nz) = ccrs%h2o(1:nlamda,1:nz)/normc
           IF (do_h2opsl) dadp%h2o(:,:,:) = h2o%dadp 
         ELSE 
           IF (fgassidxs(i) > 0 ) THEN 
              npts = n_refspec_pts(gasidxs(k))
              allocate (refspec(npts), refwav(npts)) 
              refwav  = refspec_orig_data(gasidxs(k), 1:npts, 1)
              refspec = refspec_orig_data(gasidxs(k),1:npts, 3)
    
              fidx = MINVAL(MINLOC(lamda(1:nlamda), MASK=(lamda(1:nlamda) >= &
                         refwav(1) + 0.1 .AND. lamda(1:nlamda) <= refwav(npts)- 0.1)))
              lidx = MINVAL(MAXLOC(lamda(1:nlamda), MASK=(lamda(1:nlamda) >= &
                          refwav(1) + 0.1 .AND. lamda(1:nlamda) <= refwav(npts)- 0.1)))
              IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN
                temp = fitvar_rad(rmask_fitvar_rad(fgassidxs(k)))
                CALL BSPLINE1(  refwav - temp,refspec, npts, &
                     lamda(fidx:lidx), allcrs(fidx:lidx, nfgas1, 1),gshiwf(fidx:lidx), lidx-fidx+1, errstat)
                database_shiwf(gasidxs(k), refidx(fidx:lidx)) = gshiwf(fidx:lidx)
                database(gasidxs(k), refidx(fidx:lidx)) =allcrs(fidx:lidx, nfgas, 1)
                IF (errstat < 0) THEN
                  WRITE(*, *) modulename, ' : BSPLINE error, errstat =', errstat; RETURN
                ENDIF
              ENDIF
              deallocate (refspec, refwav) 
           ELSE
              allcrs(1:nlamda, nfgas1, 1) = database_save(gasidxs(i),refidx(1:nlamda))
           ENDIF
           DO k = 2, nz
              allcrs(1:nlamda, nfgas1, k) = allcrs(1:nlamda, nfgas1, 1)
           ENDDO
         ENDIF
     ENDIF
  ENDDO   
  deallocate(delshi, delpos, gshiwf)

  RETURN
  END SUBROUTINE get_effres_gascrs_ray

!******************************************************************************************
! Prepare high resolution spectra at fine grid: solar reference, trace gas cross sections
!   (o3 shift and o3 temperature), Raleigh scattering coefficient, depolarization factor
! just need to get once for each retrieval
! Other trace gas cross section: just need to get it once for all the retrievals if no shifts
! 1) read origianl cross section ==> o3%
! 2) interpolate onto hreswav  and wf  ==> hres_o3, hres_o3shi (mxsect, max_spec_pts)
! 3) calculate T or P dependent cross section and wf ==> o3crsz, o3dadsz, o3dadtz  (max_spec_pts, mflay)
!******************************************************************************************
  SUBROUTINE get_hres_gascrs_ray ( &
             num_iter, nw, waves, nz, ts, ps, nfgas,allcol, rhos, &
             do_o3shi, o3shi, do_tmpwf, & 
             allcrs, raycof, depol, problems)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY  : numwin, winlim, &
       n_refspec_pts, refspec_orig_data, refspec_norm, & 
       fitvar_rad, rmask_fitvar_rad, solwinfit,slitfit, nslit,yn_varyslit
  USE OMSAO_indices_module,   ONLY : &
       so2_idx, so2v_idx, o2o2_idx, o2_idx, h2o_idx, o2t2_idx, h2ot2_idx, &
       hwe_idx, asy_idx, spk_idx, max_calfit_idx
  USE ozprof_data_module,     ONLY : do_subfit, nos, oswins, osfind, &
       ncalcp, radcidxs, hcrs,&
       ngas, gasidxs, fgasidxs,  fgassidxs,&
       so2sfidx, so2vsfidx, o4sfidx, &
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
       nhresp, hreswav, &
       hres_gas, hres_gasshi, hresgabs, hresgabs0, hresray, &
       hres_raycof, hres_depol, &
       hres_o3, hres_so2, hres_o4, &
       hres_o3shi, hres_so2shi, hres_o4shi, &
       o3crsz, o3dadsz, o3dadtz, & 
       so2crsz,o4crsz, o2crsz, h2ocrsz, o2crsz0, h2ocrsz0, &
       so2dads,o4dads, o2dads, h2odads

  USE OMSAO_errstat_module
  IMPLICIT NONE

  !-----------------------------
  ! Input/output variables
  !--------------------------------
  INTEGER, INTENT (IN)                                  :: nw, nz, nfgas, num_iter ! nw = ncalcp
  LOGICAL, INTENT (IN)                                  :: do_o3shi, do_tmpwf
  REAL (KIND=dp), DIMENSION(nw), INTENT (IN )           :: waves
  REAL (KIND=dp), DIMENSION(nz), INTENT (IN )           :: ts, ps, rhos
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)    :: o3shi
  REAL (KIND=dp), DIMENSION(nfgas, nz), INTENT(IN)      :: allcol
  REAL (KIND=dp), DIMENSION(nw),    INTENT (OUT)        :: raycof, depol
  REAL (KIND=dp), DIMENSION(nw,nfgas,nz),INTENT (OUT)        :: allcrs
  LOGICAL, INTENT (OUT)                                  :: problems
  !-----------------------------
  ! Local variables
  !--------------------------------
  INTEGER :: i, j, fidx, lidx,  npts, idum, nfgas1, errstat
  LOGICAL                              :: do_shi, do_convl
  REAL (KIND=dp)                       :: tmp, so2sum, o4sum, o2sum,h2osum
  REAL (KIND=dp), DIMENSION(:),   ALLOCATABLE :: delshi, tmpwav, delpos ! (nhresp)
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: so2dadsz, o4dadsz, o2dadsz, h2odadsz
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: so2dadtz, o4dadtz!, o2dadtz, h2odadtz
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: slitfit_save, solwinfit_save
  ! Save variables
  LOGICAL, SAVE :: first = .TRUE.
  TYPE(txcrs_set), SAVE  :: o3, so2, o4
  TYPE(crsz_set), SAVE :: o2, h2o
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=19), PARAMETER :: modulename = 'GET_HRES_GASCRS_RAY'

   problems = .false. 
   ! allocating output variables
   !=====================================================
   ! Load cross section at hres
   !=====================================================
   IF (first) THEN
     allocate (hres_gas (ngas, nhresp), hres_gasshi(ngas, nhresp))
     ! Obtain hres solar reference spectra

     ! Obtain hres rayleigh scattering coefficients and depolarization factor
     CALL GET_ALL_RAYCOF_DEPOL(nhresp, hreswav(1:nhresp), hres_raycof,hres_depol)

     ! Obtainhres cross sections of other trace gases (except for  O3,  SO2, O4, O2, H2O)
     hres_gas(1:ngas, 1:nhresp) = 0.0D0
     DO i = 1, ngas
        IF (fgasidxs(i) > 0 ) THEN
           ! find indices for shift
           IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx)  .AND. use_so2dtcrs) CYCLE
           IF ( gasidxs(i) == o2o2_idx                              .AND. use_o4dtcrs) CYCLE
           IF ((gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx)   .AND. use_o2dptcrs) CYCLE
           IF ((gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx) .AND. use_h2odptcrs) CYCLE
           IF (fgassidxs(i) > 0 ) CYCLE

           
           npts = n_refspec_pts(gasidxs(i))
           fidx = MINVAL(MINLOC(hreswav(1:nhresp), MASK = (hreswav(1:nhresp) >= &
                refspec_orig_data(gasidxs(i), 1, 1) .AND. hreswav(1:nhresp)  &
                <= refspec_orig_data(gasidxs(i), npts, 1))))
           lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK = (hreswav(1:nhresp) >= &
                refspec_orig_data(gasidxs(i), 1, 1) .AND. hreswav(1:nhresp) &
                <= refspec_orig_data(gasidxs(i), npts, 1))))

           IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN
              !print * , hreswav(fidx), hreswav(lidx), refspec_orig_data(gasidxs(i), 1, 1),refspec_orig_data(gasidxs(i), npts, 1)
              CALL BSPLINE(refspec_orig_data(gasidxs(i), 1:npts, 1), &
                   refspec_orig_data(gasidxs(i), 1:npts, 2), npts, hreswav(fidx:lidx), &
                   hres_gas(i, fidx:lidx), lidx - fidx + 1, errstat)

              IF (errstat < 0) THEN
                !    print * , hreswav(fidx:lidx)
                 WRITE(*, *) modulename, ' : BSPLINE2 error, errstat = ',errstat
                  stop 1
                 problems=.true. ; return
              ENDIF
           ENDIF
        ENDIF
     ENDDO

     ! Obtain original O3 cross section (quadratic or individual T, shift, T-depen)
     CALL read_txcrs(ozabs_fname,winwav_min,winwav_max,o3)
     refspec_norm (o3_t1_idx) = 1.0 ! o3%normc
     allocate (hres_o3 (o3%nt, nhresp))
     IF (do_o3shi)  THEN
        allocate (hres_o3shi (o3%nt, nhresp))
     ELSE IF (.NOT. do_o3shi) THEN
        DO i = 1, o3%nt
           CALL BSPLINE(o3%wvl(1:o3%nw), o3%crs(1:o3%nw, i),o3%nw,& 
                hreswav(1:nhresp), hres_o3(i, 1:nhresp), nhresp, errstat)
           IF (errstat < 0) THEN
              WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
                 problems=.true. ; return
           ENDIF
        ENDDO
     ENDIF
     ! Obtain original SO2 cross section (quadratic or individual T, shift)
     IF (use_so2dtcrs) THEN
        CALL read_txcrs(so2abs_fname,winwav_min, winwav_max,so2)
        refspec_norm (so2_idx) = so2%normc
        allocate (hres_so2 (so2%nt, nhresp))
        IF (do_so2shi) THEN 
          allocate (hres_so2shi (so2%nt, nhresp))
        ELSE
          fidx=MINVAL(MINLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) > so2%wvl(1) )))
          lidx=MINVAL(MAXLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) < so2%wvl(so2%nw) )))
          DO i = 1, so2%nt
             CALL BSPLINE(so2%wvl(1:so2%nw), so2%crs(1:so2%nw,i),so2%nw, & 
                hreswav(fidx:lidx),hres_so2(i, fidx:lidx), lidx-fidx+1, errstat)
             IF (errstat < 0) THEN
                WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
                 problems=.true. ; return
             ENDIF
          ENDDO
        ENDIF
     ENDIF
     IF (use_o4dtcrs) THEN
        CALL read_txcrs(ADJUSTL(TRIM(refdbdir))//o4abs_fname,winwav_min, winwav_max,o4)
        refspec_norm (o2o2_idx) = o4%normc
        allocate (hres_o4 (o4%nt, nhresp))
        IF (do_o4shi) THEN 
          allocate (hres_o4shi (o4%nt, nhresp))
        ELSE
          fidx=MINVAL(MINLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) > o4%wvl(1) )))
          lidx=MINVAL(MAXLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) < o4%wvl(o4%nw) )))
          DO i = 1, o4%nt
             CALL BSPLINE(o4%wvl(1:o4%nw), o4%crs(1:o4%nw,i),o4%nw, & 
                hreswav(fidx:lidx),hres_o4(i, fidx:lidx), lidx-fidx+1, errstat)
             IF (errstat < 0) THEN
                WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
                 problems=.true. ; return
             ENDIF
          ENDDO
        ENDIF
     ENDIF
     first = .false.
   ENDIF

   ! allocating saving variables (==> no dellocating)
   IF (num_iter == 0) THEN 
     IF (allocated(hresgabs0)) deallocate (hresgabs0)
     IF (allocated(hresgabs)) deallocate (hresgabs)
     IF (allocated(hresray)) deallocate (hresray)
     allocate (hresgabs0(nhresp,nz), hresgabs (nhresp, nz), hresray(nhresp, nz))
     IF (allocated(o3crsz))  deallocate (o3crsz)
     IF (allocated(o3dadsz)) deallocate (o3dadsz)
     IF (allocated(o3dadtz)) deallocate (o3dadtz)
     allocate (o3crsz (nhresp, nz),o3dadsz(nhresp, nz))
     IF (do_o3shi) allocate( o3dadtz(nhresp, nz)) 
     IF (use_so2dtcrs) THEN 
       IF (allocated(so2crsz)) deallocate(so2crsz)
       IF (allocated(so2dads) )deallocate(so2dads)
       IF (allocated(so2dadsz))deallocate(so2dadsz)
       allocate (so2crsz(nhresp, nz))
       IF (do_so2shi) THEN
         allocate(so2dads(nhresp), so2dadsz(nhresp, nz))
       ENDIF
     ENDIF
     IF (use_o4dtcrs) THEN 
       IF (allocated(o4crsz)) deallocate(o4crsz)
       IF (allocated(o4dads) )deallocate(o4dads)
       IF (allocated(o4dadsz))deallocate(o4dadsz)
       allocate (o4crsz(nhresp, nz))
       IF (do_o4shi) THEN
         allocate(o4dads(nhresp), o4dadsz(nhresp, nz))
       ENDIF
     ENDIF
     IF (use_o2dptcrs) THEN 
       IF (allocated(o2crsz)) deallocate(o2crsz)
       IF (allocated(o2crsz0)) deallocate(o2crsz0)
       IF (allocated(o2dads) )deallocate(o2dads)
       IF (allocated(o2dadsz))deallocate(o2dadsz)
       allocate (o2crsz(nhresp, nz), o2crsz0(nhresp, nz))
       IF (do_o2shi) THEN
         allocate(o2dads(nhresp), o2dadsz(nhresp, nz))
       ENDIF
     ENDIF
     IF (use_h2odptcrs) THEN 
       IF (allocated(h2ocrsz)) deallocate(h2ocrsz)
       IF (allocated(h2ocrsz0)) deallocate(h2ocrsz0)
       IF (allocated(h2odads) )deallocate(h2odads)
       IF (allocated(h2odadsz))deallocate(h2odadsz)
       allocate (h2ocrsz(nhresp, nz), h2ocrsz0(nhresp,nz))
       IF (do_h2oshi) THEN
         allocate(h2odads(nhresp), h2odadsz(nhresp, nz))
       ENDIF
     ENDIF  
   ENDIF
   
   ! allocating temporary local variables 
   allocate (delshi(nhresp), tmpwav(nhresp), delpos(nhresp))

   ! Interpolate onto instrument resolution
   !======================================================================
   ! Interpolate ozone cross section
   IF (do_o3shi .AND. nos > 0) THEN
      ! determine ozone wavelength shifts
      delshi = 0.0
      IF (do_subfit) THEN
        fidx = 1
        DO j = 1, numwin
          IF (j == numwin) THEN
            lidx = nhresp
          ELSE
            tmp = (winlim(j, 2) + winlim(j+1, 1)) / 2.0
            lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) <=tmp)))
          ENDIF
          delpos(fidx:lidx) =  hreswav(fidx:lidx) - (hreswav(fidx) + hreswav(lidx)) / 2.0
          IF (osfind(j, 1) > 0) delshi(fidx:lidx) =  o3shi(j, 1)
          DO i = 2, nos
           IF (osfind(j, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx)  + &
                            o3shi(j, i) * delpos(fidx:lidx) ** (i-1)
          ENDDO
          fidx = lidx + 1
        ENDDO
      ELSE
        IF (oswins(1, 1) == 1) THEN
          fidx = 1
        ELSE
          tmp = (winlim(oswins(1, 1), 1) + winlim(oswins(1, 1) - 1, 2))/2.
          fidx = MINVAL(MINLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) >=tmp)))
        ENDIF

        IF (oswins(1, 2) == numwin) THEN
          lidx = nhresp
        ELSE
          tmp = (winlim(oswins(1, 2), 2) + winlim(oswins(1, 2) + 1, 1))/2.
          lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) <=tmp)))
        ENDIF

        delpos(fidx:lidx) =  hreswav(fidx:lidx) - (hreswav(fidx) + hreswav(lidx)) / 2.0
        IF (osfind(1, 1) > 0) delshi(fidx:lidx) =  + o3shi(1, 1)
        DO i = 2, nos
           IF (osfind(1, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx) + &
                o3shi(1, i) * delpos(fidx:lidx) ** (i-1)
        ENDDO
      ENDIF
        DO i = 1, o3%nt
          CALL BSPLINE2(o3%wvl(1:o3%nw), o3%crs(1:o3%nw, i),o3%nw, do_o3shi,& 
           hreswav(1:nhresp)-delshi, hres_o3(i, 1:nhresp), hres_o3shi(i, 1:nhresp), nhresp, errstat)
           IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
                 problems=.true. ; return
          ENDIF
        ENDDO
   ENDIF
   
   ! Get ozone cross section at each layer
    IF (num_iter == 0 .OR. do_o3shi .OR. do_tmpwf) THEN
      o3crsz(1:nhresp, 1:nz) = calc_crsz(hres_o3(1:o3%nt, 1:nhresp),& 
          o3%nt, nhresp, o3%tdepend, o3%ts(1:o3%nt), ts(1:nz), nz)
      IF (do_o3shi) THEN 
        o3dadsz(1:nhresp, 1:nz) = calc_crsz(hres_o3shi(1:o3%nt, 1:nhresp), &
          o3%nt, nhresp, o3%tdepend, o3%ts(1:o3%nt), ts(1:nz), nz)
        o3dadsz(1:nhresp, 1:nz) = o3dadsz(1:nhresp, 1:nz) / o3crsz(1:nhresp,1:nz)
      ENDIF
      IF (do_tmpwf) THEN
        o3dadtz(1:nhresp, 1:nz) = calc_tmpwf(hres_o3(1:o3%nt, 1:nhresp), & 
         o3%nt, nhresp,o3%tdepend, o3%ts(1:o3%nt),ts(1:nz), nz)
        o3dadtz(1:nhresp, 1:nz) = o3dadtz(1:nhresp, 1:nz) / o3crsz(1:nhresp,1:nz)
      ENDIF
      o3crsz(1:nhresp, 1:nz) = o3crsz(1:nhresp, 1:nz) * o3%normc
   ENDIF

   ! Saved O3 crs variables
   allcrs(1:ncalcp, 1, 1:nz) = o3crsz(radcidxs(1:ncalcp), 1:nz)
   DO i = 1, nz
       hresgabs(1:nhresp, i) = o3crsz(1:nhresp, i) * allcol(1, i)
       hresgabs0(1:nhresp, i) = hresgabs(1:nhresp,i)
   ENDDO

   ! Get SO2 cross sections
   IF (use_so2dtcrs) THEN
     fidx=MINVAL(MINLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) >so2%wvl(1) )))
     lidx=MINVAL(MAXLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) <so2%wvl(so2%nw) )))
     IF (do_so2shi) THEN
       DO i = 1, so2%nt
         idum = MAX(so2sfidx, so2vsfidx)
         tmp = fitvar_rad(rmask_fitvar_rad(idum))
         CALL BSPLINE2(so2%wvl(1:so2%nw), so2%crs(1:so2%nw, i),so2%nw,do_so2shi, &
                 hreswav(fidx:lidx) - tmp, hres_so2(i, fidx:lidx), hres_so2shi(i,fidx:lidx), &
                 lidx-fidx+1, errstat)
         IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
                 problems=.true. ; return
         ENDIF
       ENDDO
     ENDIF
     IF (num_iter == 0 .OR. do_so2shi .OR. do_tmpwf) THEN
       so2crsz = 0.0 !; so2dadsz=0.0
       so2crsz (fidx:lidx, 1:nz) = calc_crsz (hres_so2(1:so2%nt, fidx:lidx), &
       so2%nt,lidx-fidx+1, so2%tdepend, so2%ts(1:so2%nt), ts(1:nz), nz ) 
       IF (do_so2shi) THEN  
        so2dadsz(fidx:lidx,1:nz)= calc_crsz(hres_so2shi(1:so2%nt,fidx:lidx),& 
         so2%nt,lidx-fidx+1, so2%tdepend, so2%ts(1:so2%nt), ts(1:nz), nz ) 
        so2dadsz(1:nhresp, 1:nz) = so2dadsz(1:nhresp, 1:nz) /so2crsz(1:nhresp, 1:nz)
       ENDIF
       IF (do_tmpwf) THEN 
        so2dadtz(1:nhresp, 1:nz) = calc_tmpwf(hres_so2(1:so2%nt, 1:nhresp), & 
         so2%nt, nhresp,so2%tdepend, so2%ts(1:so2%nt),ts(1:nz), nz)
        so2dadtz(1:nhresp, 1:nz) = so2dadtz(1:nhresp, 1:nz) / so2crsz(1:nhresp,1:nz)
       ENDIF
        so2crsz(1:nhresp, 1:nz) = so2crsz(1:nhresp, 1:nz) * so2%normc
     ENDIF
   ENDIF
  
   ! Get O4 cross sections
   IF (use_o4dtcrs) THEN
     fidx=MINVAL(MINLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) >o4%wvl(1))))
     lidx=MINVAL(MAXLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) <o4%wvl(o4%nw) )))
     IF (do_o4shi) THEN
       DO i = 1, o4%nt
         tmp = fitvar_rad(rmask_fitvar_rad(o4sfidx))
         CALL BSPLINE2(o4%wvl(1:o4%nw), o4%crs(1:o4%nw, i),o4%nw,do_o4shi, &
                 hreswav(fidx:lidx) - tmp, hres_o4(i, fidx:lidx),hres_o4shi(i,fidx:lidx), &
                 lidx-fidx+1, errstat)
         IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
           problems=.true. ; return
         ENDIF
       ENDDO
     ENDIF
     IF (num_iter == 0 .OR. do_o4shi .OR. do_tmpwf) THEN
       o4crsz = 0.0 ! ;o4dadsz=0.0
       o4crsz (fidx:lidx, 1:nz)  = calc_crsz (hres_o4(1:o4%nt, fidx:lidx), &
       o4%nt,lidx-fidx+1, o4%tdepend, o4%ts(1:o4%nt), ts(1:nz), nz )
       IF (do_o4shi) THEN
           o4dadsz (fidx:lidx, 1:nz)  = calc_crsz (hres_o4shi(1:o4%nt,fidx:lidx), &
           o4%nt,lidx-fidx+1, o4%tdepend, o4%ts(1:o4%nt), ts(1:nz), nz )
           o4dadsz(1:nhresp, 1:nz) = o4dadsz(1:nhresp, 1:nz)/o4crsz(1:nhresp, 1:nz)
       ENDIF
       IF (do_tmpwf) THEN
        o4dadtz(1:nhresp, 1:nz) = calc_tmpwf(hres_o4(1:o4%nt, 1:nhresp), & 
         o4%nt, nhresp,o4%tdepend,o4%ts(1:o4%nt),ts(1:nz), nz)
        o4dadtz(1:nhresp, 1:nz) = o4dadtz(1:nhresp, 1:nz) /o4crsz(1:nhresp,1:nz)
       ENDIF
       o4crsz(1:nhresp, 1:nz) = o4crsz(1:nhresp, 1:nz) * o4%normc
     ENDIF
   ENDIF 
   
   ! Get O2 cross section
   IF (use_o2dptcrs) THEN 
     IF (num_iter == 0) THEN  
       do_convl = .false.       
       IF (hres_slitwidth /= 0.0D0) do_convl = .true.
       o2crsz(1:nhresp, 1:nz)=0.0 ! ;o2dadsz(1:nhresp, 1:nz)=0.0
       o2crsz0(1:nhresp, 1:nz)=0.0 
       IF (do_convl) THEN 
        IF (yn_varyslit) THEN 
         allocate(slitfit_save(nslit, max_calfit_idx))
         slitfit_save(1:nslit, :)  = slitfit(1:nslit, :, 1)
         slitfit(:, hwe_idx, 1) = hres_slitwidth
         slitfit(:, asy_idx, 1) = 0.00
         slitfit(:, spk_idx, 1) = 2.0
        ELSE
         allocate(solwinfit_save(numwin, max_calfit_idx))
         solwinfit_save(1:numwin, :)  = solwinfit(1:numwin, :, 1)
         solwinfit(:, hwe_idx, 1) = hres_slitwidth
         solwinfit(:, asy_idx, 1) = 0.00
         solwinfit(:, spk_idx, 1) = 2.0
        ENDIF
       ENDIF
       CALL geto2_crs_hitran16 (hreswav(1:nhresp), nhresp, nhresp, nz,ts, ps, o2, do_convl, 1)
       !CALL geto2_crs_hitran16 (hreswav(1:nhresp), nhresp, nhresp, nz,ts, ps, o20, do_convl,0)
       IF (do_convl) THEN 
        IF (yn_varyslit) THEN
          slitfit(1:nslit, :,1)  = slitfit_save(1:nslit, :)
          deallocate (slitfit_save)
        ELSE
          solwinfit(1:numwin, :,1)  = solwinfit_save(1:numwin, :)
          deallocate (solwinfit_save)
        ENDIF
       ENDIF
       o2crsz(1:nhresp, 1:nz) = o2%crs
       o2crsz0(1:nhresp, 1:nz) = o2%crs
       !o2crsz(1:nhresp, 1:nz) = o2crsz(1:nhresp, 1:nz) !/refspec_norm(o2_idx)         
     ENDIF     
   ENDIF  
   
   ! Get h2o cross section
   IF (use_h2odptcrs) THEN 
     IF (num_iter == 0) THEN  
       do_convl = .false.       
       IF (hres_slitwidth /= 0.0D0) do_convl = .true.
       IF (do_convl) THEN 
         IF (yn_varyslit) THEN   
          allocate(slitfit_save(nslit, max_calfit_idx))
          slitfit_save(1:nslit, :)  = slitfit(1:nslit, :, 1)
          slitfit(:, hwe_idx, 1) = hres_slitwidth
          slitfit(:, asy_idx, 1) = 0.000
          slitfit(:, spk_idx, 1) = 2.0
         ELSE
          allocate(solwinfit_save(numwin, max_calfit_idx))
          solwinfit_save(1:numwin, :)  = solwinfit(1:numwin, :, 1)
          solwinfit(:, hwe_idx, 1) = hres_slitwidth
          solwinfit(:, asy_idx, 1) = 0.00
          solwinfit(:, spk_idx, 1) = 2.0
         ENDIF
       ENDIF
       h2ocrsz(1:nhresp, 1:nz)=0.0 !;h2odadsz(1:nhresp, 1:nz)=0.0
       h2ocrsz0(1:nhresp, 1:nz)=0.0 !;h2odadsz(1:nhresp, 1:nz)=0.0
       CALL geth2o_crs_hitran16 (hreswav(1:nhresp), nhresp, nhresp, nz,ts, ps,h2o, do_convl,1)       
       !CALL geth2o_crs_hitran16 (hreswav(1:nhresp), nhresp, nhresp, nz,ts, ps,h2o0, do_convl,0)       
       h2ocrsz = h2o%crs
       h2ocrsz0 = h2o%crs
       IF (do_convl)  THEN
        IF (yn_varyslit) THEN
          slitfit(1:nslit, :,1)  = slitfit_save(1:nslit, :)
          deallocate (slitfit_save)
        ELSE
          solwinfit(1:numwin, :,1)  = solwinfit_save(1:numwin, :)
          deallocate (solwinfit_save)
        ENDIF
       ENDIF
       !h2ocrsz(1:nhresp, 1:nz) = h2ocrsz(1:nhresp, 1:nz) !*refspec_norm(h2o_idx)         
     ENDIF
   ENDIF  

   ! Obtain high resolution cross sections of other trace gases
   IF (do_so2shi) THEN
      so2dads(1:nhresp) = 0.D0; so2sum = 0.D0
   ENDIF
   IF (do_o4shi) THEN
      o4dads(1:nhresp) = 0.D0; o4sum =0.D0
   ENDIF  
   IF (do_o2shi) THEN
     o2dads(1:nhresp) = 0.D0; o2sum =0.D0
   ENDIF  
   IF (do_h2oshi) THEN
      h2odads(1:nhresp) = 0.D0; h2osum =0.D0
   ENDIF  
   IF (allocated(hcrs%o3)) THEN 
       deallocate (hcrs%o3)
   ENDIF
   allocate(hcrs%o3(ncalcp, nz))
   hcrs%o3(:,:) = o3crsz(radcidxs(1:ncalcp), 1:nz)
   
   nfgas1 = 1
   DO i = 1, ngas
     IF (fgasidxs(i) > 0 ) THEN
       nfgas1 = nfgas1 + 1
    
       IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. use_so2dtcrs) THEN
         allcrs(1:ncalcp, nfgas1, 1:nz) = so2crsz(radcidxs(1:ncalcp),1:nz)/refspec_norm(gasidxs(i))         
         IF (do_so2shi) THEN
           DO j = 1, nz
             so2dads(1:nhresp) = so2dads(1:nhresp) + so2dadsz(1:nhresp, j) *allcol(nfgas1, j)
             so2sum = so2sum + allcol(nfgas1, j)
           ENDDO
         ENDIF
           DO j = 1, nz
             hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + so2crsz(1:nhresp,j) &
                                    * allcol(nfgas1,j)/refspec_norm(gasidxs(i))
             hresgabs0(1:nhresp, j) = hresgabs0(1:nhresp, j) + so2crsz(1:nhresp,j) &
                                    * allcol(nfgas1,j)/refspec_norm(gasidxs(i))
           ENDDO
       ELSE IF (gasidxs(i) == o2o2_idx .AND. use_o4dtcrs) THEN
         allcrs(1:ncalcp, nfgas1, 1:nz) =o4crsz(radcidxs(1:ncalcp),1:nz)/refspec_norm(gasidxs(i))
         IF (do_o4shi) THEN
           DO j = 1, nz
             o4dads(1:nhresp) = o4dads(1:nhresp) + o4dadsz(1:nhresp, j) *allcol(nfgas1, j) 
             o4sum = o4sum + allcol(nfgas1, j)
           ENDDO
         ENDIF
           DO j = 1, nz
             hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + & 
            o4crsz(1:nhresp,j)* allcol(nfgas1, j)/refspec_norm(gasidxs(i))
             hresgabs0(1:nhresp, j) = hresgabs0(1:nhresp, j) + & 
            o4crsz(1:nhresp,j)* allcol(nfgas1, j)/refspec_norm(gasidxs(i))
           ENDDO
       ELSE IF (gasidxs(i) == o2_idx .AND. use_o2dptcrs) THEN
         allcrs(1:ncalcp, nfgas1, 1:nz) = o2crsz(radcidxs(1:ncalcp),1:nz)/refspec_norm(gasidxs(i))
         IF (do_o2shi) THEN
           DO j = 1, nz
             o2dads(1:nhresp) = o2dads(1:nhresp) + o2dadsz(1:nhresp, j) *allcol(nfgas1, j)
             o2sum = o2sum + allcol(nfgas1, j)
           ENDDO
         ENDIF
           DO j = 1, nz
             hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + & 
                  o2crsz(1:nhresp,j)* allcol(nfgas1, j)/refspec_norm(gasidxs(i))
             hresgabs0(1:nhresp, j) = hresgabs0(1:nhresp, j) + & 
                  o2crsz0(1:nhresp,j)* allcol(nfgas1, j)/refspec_norm(gasidxs(i))
           ENDDO
         
       ELSE IF (gasidxs(i) == h2o_idx .AND. use_h2odptcrs) THEN
         allcrs(1:ncalcp, nfgas1, 1:nz) = h2ocrsz(radcidxs(1:ncalcp), 1:nz)/refspec_norm(gasidxs(i))
         IF (do_h2oshi) THEN
           DO j = 1, nz
             h2odads(1:nhresp) = h2odads(1:nhresp) + h2odadsz(1:nhresp, j) *allcol(nfgas1, j)
             h2osum = h2osum + allcol(nfgas1, j)
           ENDDO
         ENDIF
         DO j = 1, nz
           hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + & 
                                   h2ocrsz(1:nhresp,j)* allcol(nfgas1, j)/refspec_norm(gasidxs(i))

           hresgabs0(1:nhresp, j) = hresgabs0(1:nhresp, j) + & 
                                   h2ocrsz0(1:nhresp,j)* allcol(nfgas1, j)/refspec_norm(gasidxs(i))
         ENDDO
       ELSE
         IF (fgassidxs(i) > 0) THEN
           do_shi = .TRUE.
           idum = rmask_fitvar_rad(fgassidxs(i))
           tmpwav = hreswav(1:nhresp) - fitvar_rad(idum)  ! wavelength shifts
           hres_gas(i, 1:nhresp) = 0.0D0
           npts = n_refspec_pts(gasidxs(i))
           fidx = MINVAL(MINLOC(tmpwav(1:nhresp), MASK = (tmpwav(1:nhresp) >= &
                   refspec_orig_data(gasidxs(i), 1, 1) .AND. tmpwav(1:nhresp)  &
                   <= refspec_orig_data(gasidxs(i), npts, 1))))
           lidx = MINVAL(MAXLOC(tmpwav(1:nhresp), MASK = (tmpwav(1:nhresp) >=&
                   refspec_orig_data(gasidxs(i), 1, 1) .AND. tmpwav(1:nhresp) &
                   <= refspec_orig_data(gasidxs(i), npts, 1))))

           IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN
             CALL BSPLINE2(refspec_orig_data(gasidxs(i), 1:npts, 1), &
                 refspec_orig_data(gasidxs(i), 1:npts, 2), npts, do_shi,tmpwav(fidx:lidx), &
                 hres_gas(i, fidx:lidx), hres_gasshi(i, fidx:lidx), lidx -fidx + 1, errstat)
             IF (errstat < 0) THEN
               WRITE(*, *) modulename, ' : BSPLINE2 error, errstat = ',errstat
               problems = .true. ; return
             ENDIF
           ENDIF
         ENDIF
         DO j = 1, nz
           allcrs(1:ncalcp, nfgas1, j) = hres_gas(i, radcidxs(1:ncalcp))
           hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + hres_gas(i,1:nhresp) * allcol(nfgas1, j)
           hresgabs0(1:nhresp, j) = hresgabs0(1:nhresp, j) + hres_gas(i,1:nhresp) * allcol(nfgas1, j)
         ENDDO
       ENDIF ! which gas ?
     ENDIF ! fgas > 0
   ENDDO ! loop of ngas
   IF (do_so2shi .AND. use_so2dtcrs) so2dads(1:nhresp) = so2dads(1:nhresp) / so2sum
   IF (do_o4shi  .AND. use_o4dtcrs)   o4dads(1:nhresp) = o4dads(1:nhresp) / o4sum
   IF (do_o2shi  .AND. use_o2dptcrs)  o2dads(1:nhresp) = o2dads(1:nhresp) / o2sum
   IF (do_h2oshi .AND. use_h2odptcrs)h2odads(1:nhresp) = h2odads(1:nhresp) / h2osum

    ! Rayleigh scattering and depolarization factor
    raycof(1:ncalcp) = hres_raycof(radcidxs(1:ncalcp))
    depol(1:ncalcp)  = hres_depol(radcidxs(1:ncalcp))
    ! compute total rayleigh and absorption optical thickness (o3 + other gas) for
    ! later radiance interpolaiton
    DO i = 1, nz
       hresray(1:nhresp, i) = rhos(i) * hres_raycof(1:nhresp)
    ENDDO  
   ! Finish with deallocating local variables
   deallocate (delshi, tmpwav, delpos)
   RETURN
  END SUBROUTINE get_hres_gascrs_ray

  ! Read T-dependent cross sections (O3, SO2, O4) at original grids
  !----------------------------------------------------------------------------------
  ! Note : temperature dependent cross section data should have the follwoing form:
  !----------------------------------------------------------------------------------
  ! START OF TABLE
  ! 1402       238.95810       395.02670       1.0E-20
  !  F 5 F ; 1: Temperature dependent? 2. Number of coefficients/temperature 3. Convolve?
  ! 203.00  223.00  243.00  273.00  293.00 !nm, cm2/molec
  ! abscrs X normc = cm2/molec
  
  SUBROUTINE read_txcrs (absfname,minw, maxw, txcrs)
    IMPLICIT NONE
  
    ! Input variables
    CHARACTER (LEN=maxchlen)        :: absfname
    REAL (KIND=dp), INTENT(IN)      :: maxw, minw
    ! Output variables
    TYPE(txcrs_set), INTENT(OUT)    :: txcrs
    ! Local variables    
    LOGICAL                         :: file_exist
    INTEGER                         :: nline, i, j
    REAL (KIND=dp)                  :: tmp
    REAL (KIND=dp), DIMENSION (:),   ALLOCATABLE :: ts
    REAL (KIND=dp), DIMENSION (:),   ALLOCATABLE :: wvl
    REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: crs
    CHARACTER (LEN=14)              :: tmpchar
    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=10), PARAMETER    :: modulename = 'READ_TXCRS'
  
    !----------------------------
    ! initalize & allocation
    !----------------------------
    tmpchar = ' ' 
    INQUIRE (FILE= TRIM(ADJUSTL(absfname)), EXIST= file_exist)
    IF (file_exist) THEN
        WRITE(www_lun,'(A)') TRIM(ADJUSTL(absfname))
        OPEN(UNIT = ozabs_unit, file=TRIM(ADJUSTL(absfname)), status='old')
    ELSE
        WRITE(*,*) modulename//': No file of '//TRIM(ADJUSTL(absfname)); stop 1
    ENDIF

    DO WHILE (tmpchar /= 'START OF TABLE') 
       READ (ozabs_unit, '(A14)') tmpchar
    ENDDO
     
    READ (ozabs_unit, *) nline, tmp, tmp, txcrs%normc
    READ (ozabs_unit, *) txcrs%tdepend, txcrs%nt, txcrs%slitconv
    allocate (ts(txcrs%nt),wvl(nline), crs(nline, txcrs%nt))
    READ (ozabs_unit, *) ts(1:txcrs%nt)
    IF ((txcrs%nt > 1 .AND. .NOT. txcrs%tdepend ) .AND. (MINVAL(ts(1:txcrs%nt)) > 220. &
        .OR.  MAXVAL(ts(1:txcrs%nt)) < 280. )) THEN
       WRITE(*, *) modulename//': Temperature range for X-section not enough!!!';stop 1
    ENDIF

    j = 1
    DO i = 1, nline
       READ(ozabs_unit, *) wvl(j), crs(j, 1:txcrs%nt)
       !IF (txcrs%crs(1, j) < 0 ) txcrs%crs(:,j) = 0.0
       IF (wvl(j) > minw  .AND. wvl(j) <  maxw .AND. &
          .NOT. (wvl(j) > rm_uv .AND. wvl(j) < rm_vis) )  j = j + 1
    ENDDO
    CLOSE (ozabs_unit)
    txcrs%nw = j - 1    
    allocate (txcrs%ts(txcrs%nt), txcrs%wvl(txcrs%nw), & 
              txcrs%crs(txcrs%nw, txcrs%nt))
    txcrs%ts  = ts(1:txcrs%nt)
    txcrs%wvl = wvl(1:txcrs%nw)
    txcrs%crs = crs(1:txcrs%nw, 1:txcrs%nt)
    deallocate (wvl, crs, ts)
    WRITE(www_lun,*) 'Fname:'//TRIM(ADJUSTL(absfname))
    WRITE(www_lun,*) 'Nw:', txcrs%nw ,'Nt:', 'normc', txcrs%nt, txcrs%normc
   
  RETURN
  END SUBROUTINE read_txcrs

  FUNCTION calc_crsz (crs, nt, ncrs, tdepend, ts, tz, nz) RESULT (crsz)                       
     IMPLICIT NONE
     !------------------------------
     !INPUT 
     !------------------------------
     LOGICAL, INTENT(IN) :: tdepend
     INTEGER, INTENT(IN) :: nt, ncrs, nz
     REAL(KIND=dp), INTENT(IN)   :: ts(nt), tz(nz)
     REAL(KIND=dp), INTENT(IN)   :: crs(nt, ncrs)
     !-----------------------------
     !OUTPUT VARIABLES
     !-----------------------------
     REAL (KIND=dp), DIMENSION(ncrs, nz) :: crsz
     !local variables
     INTEGER :: i , errstat
     REAL (KIND=dp) :: thet, frac

     crsz = 0.0 
     IF (tdepend .AND. nt == 3) THEN     ! qudratic T dependent coefficients
        DO i = 1, nz
           thet = tz(i) - zerok
           crsz(1:ncrs, i) = crs(1, 1:ncrs) + crs(2, 1:ncrs) *thet &
                              + crs(3, 1:ncrs) * thet * thet 
        ENDDO
     ELSE IF (tdepend .AND. nt == 2) THEN   ! linear T dependent coefficients
        DO i = 1, nz
           thet = tz(i) - zerok
           crsz(1:ncrs, i) = crs(1, 1:ncrs) + crs(2, 1:ncrs) *thet
        ENDDO
     ELSE IF (tdepend .AND. nt == 1)  THEN           ! only 1 T
        DO i = 1, nz
           crsz(1:ncrs, i) = crs(1, 1:ncrs)
        ENDDO
     ELSE IF (.NOT. tdepend .AND. nt == 2) THEN            ! have 2 T values
        DO i = 1, nz
           frac = 1.0 - (tz(i) - ts(1)) / (ts(2) - ts(1))
           crsz(1:ncrs, i) = (frac * crs(1, 1:ncrs) + (1.0 - frac) *crs(2, 1:ncrs))
        END DO
     ELSE  IF (.NOT. tdepend .AND. nt >= 3) THEN  ! have more than n T
        DO i = 1, ncrs                           ! Interpolate/extrapolate over T  
           CALL INTERPOL(ts(1:nt), crs(1:nt,i), nt, tz(1:nz), &
                crsz(i, 1:nz), nz, errstat)
           IF (errstat < 0) THEN
              WRITE(*, *) 'calc_crsz:  INTERPOL2 error, errstat = ', errstat ;stop 1
           ENDIF
        ENDDO
     ELSE
        WRITE(*, *) 'calc_crsz:Such type of ozone cross sections not implemented'; stop 1
     ENDIF
     RETURN
  END FUNCTION calc_crsz

  FUNCTION calc_tmpwf (crs, nt, ncrs, tdepend, ts, tz, nz) RESULT (dadtz)                       
     IMPLICIT NONE
     !INPUT VARIABLES
     LOGICAL, INTENT(IN) :: tdepend
     INTEGER, INTENT(IN) :: nt, ncrs, nz
     REAL(KIND=dp), INTENT(IN)  :: ts(nt), tz(nz)
     REAL(KIND=dp), INTENT(IN)  ::  crs(nt, ncrs)
     !OUTPUT VARIABLES
     REAL (KIND=dp), DIMENSION(ncrs, nz) ::  dadtz
     !local variables
     INTEGER :: i , errstat
     REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: tmpcrs 
     REAL (KIND=dp) :: thet, frac
     dadtz = 0.0
     IF (tdepend .AND. nt == 3) THEN     ! qudratic T dependent coefficients
        DO i = 1, nz
           thet = tz(i) - zerok
           dadtz(1:ncrs, i) = crs(2, 1:ncrs)  + 2.0 * crs(3,1:ncrs) * thet
        ENDDO
     ELSE IF (tdepend .AND. nt == 2) THEN   ! linear T dependent coefficients
        DO i = 1, nz
           thet = tz(i) - zerok
           dadtz(1:ncrs, i) = crs(2, 1:ncrs)
        ENDDO
     ELSE IF (tdepend .AND. nt == 1)  THEN           ! only 1 T
        DO i = 1, nz
           dadtz(1:ncrs, i) = 0.0
        ENDDO
     ELSE IF (.NOT. tdepend .AND. nt == 2) THEN            ! have 2 T values
        DO i = 1, nz
           frac = 1.0 - (tz(i) - ts(1)) / (ts(2) - ts(1))
           dadtz(1:ncrs, i) = (crs(1, 1:ncrs) -crs(2, 1:ncrs)) / (ts(1) - ts(2))
        END DO
     ELSE  IF (.NOT. tdepend .AND. nt > 3) THEN  ! have more than n T
        allocate(tmpcrs(nz))
        DO i = 1, ncrs  ! Interpolate/extrapolate over T  
           CALL INTERPOL2(ts(1:nt), crs(1:nt,i), nt, .TRUE., tz(1:nz), &
                tmpcrs(1:nz), dadtz(i, 1:nz), nz, errstat)
           IF (errstat < 0) THEN
              WRITE(*, *) ' INTERPOL2 error, errstat = ', errstat;stop 1
           ENDIF
        ENDDO
        deallocate(tmpcrs)
     ELSE
        WRITE(*, *)  'Such type of ozone cross sections not implemented'
        stop 1
     ENDIF
     RETURN
  END FUNCTION calc_tmpwf

  ! the unit of wf = dI/dp/I
  FUNCTION calc_pslwf (hwav, habs,  nhres, npsl, do_i0convol, scalex, cwav, ncres) RESULT (dadp) 
     USE OMSAO_indices_module,   ONLY: hwe_idx, spk_idx, asy_idx
     USE OMSAO_variables_module, ONLY: numwin, psl_fpos, nradpix_sav,solwinfit, solwinfit_save, & 
         do_dsdw, do_dsdk, do_dsda
     IMPLICIT NONE
     !INPUT VARS
     LOGICAL, INTENT(IN) :: do_i0convol
     INTEGER, INTENT(IN) :: nhres, npsl, ncres
     REAL (KIND=DP), INTENT(IN) :: scalex
     REAL (KIND=DP), DIMENSION (nhres), INTENT(IN)     :: hwav
     REAL (KIND=DP), DIMENSION (nhres), INTENT(IN)     :: habs
     REAL (KIND=DP), DIMENSION (ncres), INTENT(IN)     :: cwav
     !OUTPUT VARS
     REAL (KIND=DP), DIMENSION (ncres, npsl) :: dadp
     !LOCAL VARIABLES
     INTEGER, PARAMETER :: which_pslwf = 2
     INTEGER :: idx , j, k, fidx, lidx
     REAL (KIND=DP), DIMENSION(:), ALLOCATABLE :: hi0, hspec
     REAL (KIND=DP), DIMENSION(:), ALLOCATABLE :: ci0, cspec,  absc, dfabs
    
     allocate (hi0(nhres), hspec(nhres))
     allocate (ci0(ncres), cspec(ncres), absc(ncres), dfabs(ncres))
     call get_i0(nhres, hwav, hi0)
     IF (which_pslwf == 1) THEN       
         hspec = habs
         ! convolve hres without perturbation
         IF (do_i0convol) THEN 
           CALL convol_i0f2c(hwav, hspec, nhres, 1,scalex, cwav, cspec,ncres)
         ELSE
           CALL convol_f2c(hwav, hspec, nhres, 1, cwav, cspec,ncres)
         ENDIF
         absc = cspec 

         solwinfit_save = solwinfit
         DO k = 1, npsl
            hspec = habs                
            idx = psl_fpos(k)
            solwinfit(1:numwin, idx, 1) = solwinfit_save(1:numwin, idx,1) *1.001
            IF (do_i0convol) THEN 
               hspec = hspec*hi0
               CALL convol_f2c(hwav, hspec, nhres, 1, cwav, cspec,ncres)
               CALL convol_f2c(hwav, hi0, nhres, 1, cwav, ci0,ncres)
               dfabs =cspec/ci0 - absc
            ELSE
               CALL convol_f2c(hwav, hspec, nhres, 1, cwav, cspec,ncres)
               dfabs =(cspec - absc)/absc
            ENDIF
            solwinfit = solwinfit_save
            fidx = 1
            DO j = 1, numwin
              lidx = fidx + nradpix_sav(j) - 1
              dadp(fidx:lidx, k) =  dfabs(fidx:lidx)/(solwinfit(j, idx, 1)*0.001)
              fidx = lidx + 1
            ENDDO
         ENDDO
     ELSE
         hspec = habs
         ! convolve hres without perturbation
         do_dsdw = .false. ; do_dsdk=.false. ; do_dsda = .false.
         IF (do_i0convol) THEN 
           CALL convol_i0f2c(hwav, hspec, nhres, 1,scalex, cwav, cspec,ncres)
         ELSE
           CALL convol_f2c(hwav, hspec, nhres, 1, cwav, cspec,ncres)
         ENDIF
         absc = cspec 
         
         DO k = 1, npsl
           idx = psl_fpos(k)
           IF ( idx == hwe_idx) THEN
             do_dsdw = .true.  ; do_dsdk = .false. ; do_dsda = .false.
           ELSE IF (idx == spk_idx) THEN
             do_dsdw = .false. ; do_dsdk = .true. ; do_dsda = .false.
           ELSE IF (idx == asy_idx) THEN 
             do_dsdw = .false. ; do_dsdk = .false. ; do_dsda = .true.
           ENDIF 
           hspec=habs
           IF (do_i0convol) THEN 
             hspec = hspec*hi0
             CALL convol_f2c(hwav, hspec, nhres, 1, cwav, cspec,ncres)
             CALL convol_f2c(hwav, hi0, nhres, 1, cwav, ci0,ncres)
             cspec = cspec/ci0/absc
           ELSE
             CALL convol_f2c(hwav, hspec, nhres, 1, cwav, cspec,ncres)
           ENDIF
             dadp(1:ncres, k) = cspec/absc
         ENDDO
           do_dsdw = .false. ; do_dsdk = .false. ; do_dsda = .false.          
    ENDIF     
    deallocate(hi0, hspec, ci0, cspec, absc, dfabs) 
  RETURN
    
  END FUNCTION calc_pslwf

    SUBROUTINE read_hitran16_lut(abs_fname, win_min, win_max,lut)

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    !INTEGER, PARAMETER :: maxh_ht = 52002
    !----------------------------------
    !Input
    !---------------------------------
    CHARACTER(LEN=*), INTENT (IN) :: abs_fname
    REAL (KIND=dp), INTENT(IN)    :: win_min, win_max
    !---------------------------------
    !Output 
    !---------------------------------
    TYPE(hitran16_set), INTENT(OUT) :: lut
    !---------------------------------
    !Local variables
    !---------------------------------
    LOGICAL :: file_exist, fail
    INTEGER :: rcode,ncid,  vid, dimid(4), wmx, pmx, bmx, tmx, fidx, lidx
    CHARACTER (LEN=100) :: tmpchar
    !-------------------------------------------------------------------
    ! LUT variables
    !--------------------------------------------------------------------
    !INTEGER :: nwvl_lut, np_lut, nT_lut
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: p_lut, wvl_lut, br_lut
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:):: T_lut
    !REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:)  :: xs_lut
    ! =====================================================================
    ! InitType_pTLookupTable starts here
    ! =====================================================================

    ! Open file in read mode
    INQUIRE (FILE= TRIM(ADJUSTL(abs_fname)), EXIST= file_exist)
    IF (.not. file_exist) THEN
        write(*,*) "No exsit:"//abs_fname ; stop 1
    ENDIF 
    ncid = ncopn(trim(adjustl(abs_fname)), nf_Nowrite, rcode)
    WRITE(*,*) trim(adjustl(abs_fname))
    ! Read the grid dimensions
    rcode = nf_inq_varid(ncid, 'CrossSection', vid)
    rcode = nf_inq_vardimid(ncid, vid, dimid)
    rcode = nf_inq_dim(ncid, dimid(1),tmpchar, wmx)
    rcode = nf_inq_dim(ncid, dimid(2),tmpchar, bmx)
    rcode = nf_inq_dim(ncid, dimid(3),tmpchar, tmx)
    rcode = nf_inq_dim(ncid, dimid(4),tmpchar, pmx)

    ! Allocate grid arrays
    ALLOCATE(wvl_lut(wmx))
    ALLOCATE(br_lut(bmx))
    ALLOCATE(p_lut(pmx))
    ALLOCATE(t_lut(tmx, pmx))

    ! Read Grids
    rcode = nf_inq_varid(ncid, 'Wavelength', vid)
    rcode = nf_get_var_double(ncid, vid,wvl_lut)
    rcode = nf_inq_varid(ncid, 'Pressure', vid)
    rcode = nf_get_var_double(ncid, vid,p_lut)
    rcode = nf_inq_varid(ncid, 'Temperature', vid)
    rcode = nf_get_var_double(ncid,vid,t_lut)
    rcode = nf_inq_varid(ncid, 'Broadener_01_VMR', vid)
    rcode = nf_get_var_double(ncid, vid,br_lut)
    ! ----------
    ! Close file
    ! ----------

    call ncclos(ncid, rcode)
    if (rcode .eq. -1) then
       !message =  ' error in netcdf_rd_dim: ncclos'
       fail = .true.; return
    endif
    ! ================================================================
    ! subtrac for ozone fitting window
    ! ================================================================
     ! Get index range covering window
    fidx = MINVAL(MINLOC(wvl_lut, MASK= wvl_lut .GE. win_min))
    lidx = MINVAL(MAXLOC(wvl_lut, MASK= wvl_lut .LE. win_max))

    lut%widx0 = fidx
    lut%widxf = lidx

    lut%wmx = lidx - fidx + 1 
    lut%tmx = tmx
    lut%pmx = pmx
    lut%bmx = bmx
    lut%ncid = ncid
    allocate (lut%wvl(lut%wmx), lut%ps(pmx), lut%ts(tmx, pmx), &
              lut%crs(lut%wmx, bmx, tmx, pmx))

    lut%wvl = wvl_lut(fidx:lidx)
    lut%ps  = p_lut
    lut%ts  = t_lut
    lut%br  = br_lut
                
    deallocate (wvl_lut, br_lut, t_lut, p_lut)     
  END SUBROUTINE read_hitran16_lut

  SUBROUTINE calc_hitran_crsz (lut,nw_lut, nt_lut, np_lut, &
                               t_lut, p_lut,  & 
                               nz, temp,press, crsz)
   IMPLICIT NONE
  !------------------------------
  ! input variables
  !------------------------------
   INTEGER, INTENT(IN) :: nz, nw_lut, nt_lut, np_lut
   REAL (KIND=dp), DIMENSION (nw_lut,nt_lut, np_lut), INTENT(IN) :: lut
   REAL (KIND=dp), DIMENSION (nt_lut), INTENT(IN) :: t_lut
   REAL (KIND=dp), DIMENSION (np_lut), INTENT(IN) :: p_lut
   REAL (KIND=dp), DIMENSION (nz), INTENT(IN) :: press, temp
  !------------------------------
  ! output variables
  !------------------------------
   REAL (KIND=dp), DIMENSION (nw_lut, nz), INTENT(OUT) :: crsz
  !-----------------------------
  ! local variables
  !-----------------------------
   INTEGER :: i, pidx, tidx, np, nt
   REAL (KIND=dp), DIMENSION (nw_lut, nt_lut) :: tmpcrs
   REAL (KIND=dp) :: pmin, pmax, tmin, tmax
   REAL (KIND=sp) :: pfrac, tfrac

  ! Get max/min p and T
    pmin = minval(p_lut)
    pmax = maxval(p_lut)
    Tmin = minval(T_lut)
    Tmax = maxval(T_lut)  
   
    tmpcrs = 0.0D0
    ! interpolation
    DO i = 1, nz
    pidx = minval(maxloc(p_lut(1:np_lut), MASK=(p_lut(1:np_lut) < press(i))))
    np = 1
    IF (press(i) <= pmin) THEN
        pidx = 1 ; pfrac = 1.0 
    ELSE IF (press(i) >= pmax) THEN
        pidx = np_lut ; pfrac = 1.0 
    ELSE
        np = 2
        pfrac = 1.0- (press(i)-p_lut(pidx))/(p_lut(pidx+1) - p_lut(pidx)) 
    ENDIF
    IF (np == 1) THEN 
       tmpcrs(1:nw_lut, 1:nt_lut) = lut(1:nw_lut,1:nt_lut, pidx)
    ELSE
       tmpcrs(1:nw_lut, 1:nt_lut) = lut(1:nw_lut,1:nt_lut, pidx)*pfrac+ & 
                                    lut(1:nw_lut,1:nt_lut, pidx+1)*(1.0-pfrac)
    ENDIF
    tidx = minval(maxloc(t_lut(1:nt_lut), MASK=(t_lut(1:nt_lut) < temp(i))))
    nt = 1
    IF (temp(i) <= tmin) THEN
       tidx = 1 ; tfrac = 1.0 
    ELSE IF (temp(i) >= tmax) THEN
       tidx = nt_lut ; tfrac = 1.0 
    ELSE
    nt = 2
    tfrac = 1.0- (temp(i)-t_lut(tidx))/(t_lut(tidx+1) - t_lut(tidx))
    ENDIF
    !WRITE(*,'(i3, 2f8.3, 2f15.2, 2f15.2)') i, pfrac, tfrac,  press(i), p_lut(pidx), temp(i),  t_lut(tidx)
    IF (nt == 1) THEN 
       crsz(1:nw_lut, i) = tmpcrs(1:nw_lut, tidx)
    ELSE
       crsz(1:nw_lut, i) = tmpcrs(1:nw_lut, tidx)*tfrac + tmpcrs(1:nw_lut, tidx+1)*(1.0-tfrac)
    ENDIF
  ENDDO
  RETURN
  END SUBROUTINE calc_hitran_crsz

  END MODULE m_get_xcrs

!SUBROUTINE GET_RAYCOF_DEPOL(lamda, raycof, depol)
!
!  USE OMSAO_precision_module
!  IMPLICIT NONE
!
!  !    Input/Output
!  REAL (KIND=dp), INTENT(IN)  :: lamda
!  REAL (KIND=dp), INTENT(OUT) :: raycof, depol
!
!  !    Local variables
!  REAL (KIND=dp) :: sig, sig2, sig2p, sig4, fk_n2, fk_o2, fking
!  REAL (KIND=dp), PARAMETER :: abod = 1.0455996d0, bbod = -341.29061d0, &
!       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
!
!  !    Rayleigh coefficient
!  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
!  sig =    1.0d3 / lamda
!  sig2 =   sig * sig
!  sig2p =  1.d0 / sig2
!  sig4 =   sig2 * sig2
!  raycof = (abod + bbod * sig2 + cbod * sig2p) &
!       / (1.d0 + dbod * sig2 + ebod * sig2p) * 1.d-28
!
!  !    Derivation of depolarization factor d from king factors for air,
!  !    fking = (6 + 3.depol) / (6 - 7.depol)
!  !    bodhaine et al., 370 ppmv co2
!  fk_n2 = 1.034d0 + 3.17d-4 * sig2
!  fk_o2 = 1.096d0 + 1.385d-3 * sig2 + 1.448d-4 * sig4
!  !fk_ar = 1.00
!  !fk_co2 = 1.15d0
!  !fking = (78.084d0 * fk_n2 + 20.946d0 * fk_o2 + 0.934d0 * fk_ar + &
!  !     0.037d0 * fk_co2) / (78.084d0 + 20.946d0 + 0.934d0 + 0.037d0)
!  fking = (78.084d0 * fk_n2 + 20.946d0 * fk_o2 + 0.97655d0) / 100.001d0
!  depol = 6.d0 * (fking - 1.d0) / (3.d0 + 7.d0 * fking)
!
!  RETURN
!
!END SUBROUTINE GET_RAYCOF_DEPOL


