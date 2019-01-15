! temperature dependent cross section should "crs*norm" and normc=1

MODULE m_get_xcrs
   
   USE OMSAO_precision_module  
   USE OMSAO_parameters_module, ONLY : max_spec_pts, maxchlen, zerok
   USE OMSAO_indices_module, ONLY: so2_idx, so2v_idx, o2o2_idx, &
       o2_idx, o2t2_idx, h2o_idx, h2ot2_idx, &
       solar_idx, wvl_idx, spc_idx, hwe_idx, spk_idx, &
                                   o3_t1_idx
   USE OMSAO_variables_module,  ONLY : refdbdir, ozabs_unit,do_bandavg, &
       do_dsdw, do_dsdk,winwav_min,winwav_max, &
       n_refspec_pts, refspec_orig_data,  refspec_norm, database,refidx
   USE OMSAO_errstat_module
   USE ozprof_data_module,      ONLY : mxsect,mflay,ozabs_fname, &
       ozcrs_alb_fname,pos_alb, toms_fwhm, &
       ozabs_convl, so2crs_convl, o4crs_convl,o2crs_convl, h2ocrs_convl

   USE m_ezspline_interpolation, only: &
       bspline, bspline2, bspline1,interpol, interpol2, interpolation
   USE m_convol, ONLY: convol_f2c, convol_i0effect, convolf2c_i0effect, convol
   USE m_avg_band, ONLY: avg_band_effozcrs

   IMPLICIT NONE
   ! Name of T-dependent tracegases cross section file  
   CHARACTER(maxchlen), PARAMETER :: so2abs_fname='OMSAO_SO2_scia_fm.dat'
   CHARACTER(maxchlen), PARAMETER :: o4abs_fname='OMSAO_Thalman_O4quad_extended654nm.dat'
   !CHARACTER(maxchlen), PARAMETER :: o4abs_fname='OMSAO_Thalman_O4ts_extended654nm.dat'
   CHARACTER(maxchlen), PARAMETER :: h2oabs_fname='hitran_lut/h2o_lut_280-800nm_0p04fwhm.nc'
   CHARACTER(maxchlen), PARAMETER :: o2abs_fname='hitran_lut/o2_lut_280-800nm_0p04fwhm.nc'

   !----------------------------------------------------------------
   ! Typed cross section for O3, O4, So2, T-dependent cross section
   !------------------------------------------------------------------
   TYPE  txcrs_set
   INTEGER         :: nw, nt
   LOGICAL         :: tdepend, slitconv
   REAL (KIND=dp)  :: normc
   REAL (KIND=dp), DIMENSION(mxsect)      :: ts
   REAL (KIND=dp), DIMENSION(:),  POINTER :: wvl ! before convolution
   REAL (KIND=dp), DIMENSION(:,:),POINTER :: crs0 ! before convolution
   REAL (KIND=dp), DIMENSION(:,:),POINTER :: crs  ! after convolution
   END TYPE  txcrs_set

   !-----------------------------------------------------------------
   ! T-P dependent hitran LUT
   !-----------------------------------------------------------------
   INTEGER, PARAMETER :: maxw_ht=max_spec_pts, maxt_ht=19, maxp_ht=9
   REAL (KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: o2_lut, h2o_lut, & !
     o2_lut0,h2o_lut0

   PUBLIC  geto3_crs, &           ! called in raman
           get_all_raycof,  &     ! called in raman
           get_alb_ozcrs_ray, &   ! nw = 1 for 347 nm
           get_hres_gascrs_ray, & ! nw > 1 without convolution
           get_effres_gascrs_ray  ! nw > 1 with convolution

   PRIVATE
   !PRIVATE read_txcrs, calc_crsz, & 
   !        read_hitran_lut,calc_hitran_crsz, & 
   !        get_all_raycof_depol,get_all_raycof_depol1
   !        getabs_crs_hitran, geto4_crs, getso2_crs

CONTAINS 

  SUBROUTINE allocate_txcrs (crs)
    TYPE (txcrs_set), INTENT(INOUT) :: crs
    allocate (crs%wvl  (max_spec_pts))
    allocate (crs%crs0 (mxsect, max_spec_pts))
    allocate (crs%crs  (mxsect, max_spec_pts))
  END SUBROUTINE  

  SUBROUTINE geto3_crs  (lamda, nlsav, nlamda, nlayers, tsgrid, & 
                        abscrs, dods, dodt, dads, dadt,problems)

  IMPLICIT NONE

  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN) :: nlamda, nlsav, nlayers
       ! #ofwave could vary for in/output
  LOGICAL, INTENT(IN) :: dods, dodt
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)   :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers) :: tsgrid
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: abscrs, dads, dadt
  LOGICAL, INTENT(OUT) :: problems
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  INTEGER, PARAMETER :: maxline  = max_spec_pts      ! # of wavelengths
  INTEGER, PARAMETER :: maxt = 3  ! most 3 except for o4
  INTEGER            :: i, j, errstat, ntemp,nline, nt
  REAL (KIND=dp)     :: scalex, normc
  REAL (KIND=dp), DIMENSION(maxt, nlsav)   :: savabs, savabs_d1
  REAL (KIND=dp), DIMENSION(maxt, nlamda)  :: tmpabs, tmpabs_d1
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  TYPE(txcrs_set), SAVE :: o3
  LOGICAL, SAVE :: first = .true.
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'geto3_crs'
  
  problems = .FALSE.

  ! load origianl spectrum
  IF (first) THEN
    call allocate_txcrs (o3)
    CALL read_txcrs (winwav_min, winwav_max, o3_t1_idx, maxt, o3, problems)
    o3%crs = o3%crs0
    IF (o3%wvl(1)  > lamda(1) .OR. o3%wvl(o3%nw)  < lamda(nlsav)) THEN
       WRITE(www_lun, *) modulename//': O3abs should cover the whole fit wavelenth!!!'
       problems = .TRUE. ; return
    ENDIF
    IF (.NOT. o3%slitconv ) THEN
       WRITE(www_lun, *) modulename//': Need to use high-resolution cross section for O3!!!'
       problems = .TRUE. ; return
    ENDIF
    IF (o3%slitconv) ozabs_convl = .true.     
    first = .FALSE.
  ENDIF
  nline = o3%nw ; nt = o3%nt ; normc=o3%normc
  !-----------------------------------------------------------------------------------------
  ! convolution     
  !-----------------------------------------------------------------------------------------
  IF (ozabs_convl .AND. o3%slitconv ) THEN  
    do_dsdw = .false.; do_dsdk = .false. 
    scalex = 0.2 ! ~600 DU SUM(fozs(1:nflay)) * 2.69E16 * normc, now a dummy number, not used
    ! Perform solar i0 effect on ozone cross-section (no need to convolve)
    o3%crs=o3%crs0
    DO i = 1, nt 
         CALL convol_i0effect(o3%wvl(1:nline),o3%crs(i, 1:nline), nline, scalex, errstat)
         IF ( errstat /= 0 ) THEN
           WRITE(*, *) modulename//': Error in Correct I0 Effect!!!'
           problems = .TRUE.; RETURN
         ENDIF
    ENDDO
    ozabs_convl = .FALSE.
  ENDIF

  !------------------------------------------------------------------------
  ! interpolation
  !-------------------------------------------------------------------------
  savabs = 0.0 ; savabs_d1 = 0.0  ! onto refwvl 
  tmpabs = 0.0 ; tmpabs_d1 = 0.0  ! onto re-sampled wavelengths 
  dads = 0.0; dadt = 0.0  

  DO i = 1, nt
    CALL BSPLINE2(o3%wvl(1:nline),o3%crs(i, 1:nline), nline, dods, lamda, &
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
  !database (o3_t1_idx,refidx(1:nlamda)) = tmpabs(1, 1:nlamda)
  abscrs = 0.0
  problems = calc_crsz (tmpabs(1:nt,1:nlamda),tmpabs_d1(1:nt, 1:nlamda), &
             nt, nlamda,o3%tdepend,o3%ts(1:nt),tsgrid(1:nlayers), nlayers, &
             dods,dodt ,abscrs(1:nlamda,1:nlayers), &
             dadsz=dads(1:nlamda, 1:nlayers), dadtz=dadt(1:nlamda, 1:nlayers) )

  IF (dods) dads = dads / abscrs   ! get relative sensitivty to shift
  IF (dodt) dadt = dadt / abscrs   ! get relative sensitivity to T
  abscrs = abscrs*normc

  RETURN  
  END SUBROUTINE geto3_crs

  SUBROUTINE getso2_crs  (lamda, nlsav, nlamda, nlayers, tsgrid, & 
                         abscrs, problems)

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
  REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: abscrs
  LOGICAL, INTENT(OUT) :: problems
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  INTEGER, PARAMETER :: maxline  = max_spec_pts      ! # of wavelengths
  INTEGER, PARAMETER :: maxt = 3  ! most 3 except for o4
  INTEGER            :: fidx, lidx, i, j, errstat, nline, nt, ntemp
  REAL (KIND=dp)     :: scalex, normc
  LOGICAL            :: dods, dodt
  REAL (KIND=dp), DIMENSION(maxt, nlsav)   :: savabs, savabs_d1
  REAL (KIND=dp), DIMENSION(maxt, nlamda)  :: tmpabs, tmpabs_d1
  REAL (KIND=dp), DIMENSION(nlamda, nlayers) :: dads, dadt
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  TYPE(txcrs_set), SAVE :: so2
  LOGICAL, SAVE :: first = .true.
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'getso2_crs'
  dods=.false. ; dodt = .false.
  problems = .FALSE.
  ! load origianl spectrum
  IF (first) THEN
    call allocate_txcrs (so2)
    CALL read_txcrs (winwav_min, winwav_max, so2_idx, maxt, so2, problems)    
    so2%crs = so2%crs0   
    IF (so2%slitconv) so2crs_convl = .true.     
    first = .FALSE.
  ENDIF
  nline = so2%nw ; nt = so2%nt ; normc=so2%normc
  !-----------------------------------------------------------------------------------------
  ! convolution     
  !-----------------------------------------------------------------------------------------
  IF (so2crs_convl .AND. so2%slitconv) THEN  
    scalex = 0.1! ~600 DU SUM(fozs(1:nflay)) * 2.69E16 * normc, now a dummy number, not used
    so2%crs(1:nt, 1:nline) = so2%crs0(1:nt, 1:nline)
    ! Perform solar i0 effect on ozone cross-section (no need to convolve)
    DO i = 1, nt 
         CALL convol_i0effect(so2%wvl(1:nline), so2%crs(i, 1:nline), nline, &
              scalex, errstat)
         IF ( errstat /= 0 ) THEN
           WRITE(*, *) modulename//': Error in Correct I0 Effect!!!'
           problems = .TRUE.; RETURN
         ENDIF
    ENDDO
    so2crs_convl = .FALSE.
  ENDIF

  !------------------------------------------------------------------------
  ! interpolation
  !-------------------------------------------------------------------------
  savabs = 0.0 ; savabs_d1 = 0.0  ! onto refwvl 
  tmpabs = 0.0 ; tmpabs_d1 = 0.0  ! onto re-sampled wavelengths 

  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= so2%wvl(1))))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= so2%wvl(so2%nw))))

  DO i = 1, nt
    CALL BSPLINE(so2%wvl(1:nline), so2%crs(i, 1:nline), nline, &
          lamda(fidx:lidx), savabs(i, fidx:lidx),lidx-fidx+1, errstat)
    IF (errstat < 0) THEN
          WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
          problems = .TRUE.; RETURN
    ENDIF
    IF (do_bandavg) THEN
      CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)
      IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
         WRITE(*, *) modulename//':SO2 crs Averaging Error: ', nlsav, nlamda, ntemp
         problems = .TRUE.; RETURN
      ENDIF
      tmpabs(i, :) = savabs(i, 1:nlamda)
    ELSE
      tmpabs(i, :) = savabs(i, 1:nlamda)
    ENDIF
  ENDDO

  abscrs = 0.0
  ! second tmpabs should be changed to crsshi if do_shi is true.
  problems = calc_crsz (tmpabs(1:nt,fidx:lidx),tmpabs(1:nt, fidx:lidx), &
               nt,lidx-fidx+1,so2%tdepend,so2%ts(1:nt),tsgrid(1:nlayers), nlayers, &
               dods,dodt,abscrs(fidx:lidx,1:nlayers), & 
               dadsz=dads(1:nlamda, 1:nlayers), dadtz=dadt(1:nlamda, 1:nlayers) )
  IF (dods) dads = dads / abscrs   ! get relative sensitivty to shift
  IF (dodt) dadt = dadt / abscrs   ! get relative sensitivity to T
  abscrs = abscrs * normc
  RETURN  
  END SUBROUTINE getso2_crs

  SUBROUTINE getabs_pslwf(lamda, nlsav, nlamda, nlayers, tsgrid, & 
                         slit_idx, dadp,problems)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : zerok
  USE OMSAO_variables_module, ONLY : n_refspec_pts, refspec_orig_data, &
      solwinfit, solwinfit_save, numwin, nradpix_sav, do_dsdw, do_dsdk
  USE OMSAO_indices_module,   ONLY : solar_idx, wvl_idx, spc_idx, hwe_idx, spk_idx
  USE OMSAO_slitfunction_module
  IMPLICIT NONE

  ! Input variables
  INTEGER, INTENT(IN) :: nlamda, nlsav, nlayers, slit_idx
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)   :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers) :: tsgrid

  ! Output variables
  REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: dadp
  LOGICAL, INTENT(OUT) :: problems

  ! Local variables
  INTEGER, PARAMETER :: which_pslwf = 2 ! 1 finite 2 analytic
  INTEGER, PARAMETER :: maxline  = max_spec_pts     ! # of wavelengths
  INTEGER, PARAMETER :: maxt = 3 

  REAL (KIND=dp), DIMENSION(maxline) :: specmod
  REAL (KIND=dp), DIMENSION(maxt, nlamda)  :: tmpabs
  REAL (KIND=dp), DIMENSION(nlamda, nlayers) :: abscrs

  INTEGER            :: i, j, k, errstat, npts, idx, lidx, fidx, nline
  REAL (KIND=dp)     :: thet
  REAL (KIND=dp), DIMENSION (1) :: scalex
  !----------------------------------------------
  ! Save variables
  !------------------------------------------------
  TYPE(txcrs_set), SAVE :: o3
  LOGICAL, SAVE :: first = .true.
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'getabs_pslwf'

  problems = .FALSE.
  IF (first) THEN
    CALL allocate_txcrs(o3) 
    CALL read_txcrs (winwav_min, winwav_max, o3_t1_idx, maxt, o3, problems)         
     o3%crs  = o3%crs0
     ! Obtain high resolution solar reference spectra
     IF ( problems ) THEN
        print * , modulename//'errors in read_txcrs'
        problems = .TRUE.; RETURN
     ENDIF
     ozabs_convl = .true.
     first = .FALSE.
  ENDIF
   
  nline = o3%nw
  scalex = 0.2
  IF (ozabs_convl) THEN  
    do_dsdw = .false.  ; do_dsdk = .false. 
    o3%crs=o3%crs0
    DO i = 1, o3%nt
      CALL convol_i0effect ( o3%wvl(1:nline), o3%crs(i,1:nline),nline,scalex(1), &
                             errstat)    
    ENDDO
      ozabs_convl = .false.
  ENDIF

  !---------------------------------------------------------
  ! interpolation
  !---------------------------------------------------------
  DO i = 1, o3%nt
    CALL BSPLINE(o3%wvl(1:nline), o3%crs(i, 1:nline), nline, &
          lamda(1:nlamda), tmpabs(i, 1:nlamda),nlamda, errstat)
    IF (errstat < 0) THEN
          WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
          problems = .TRUE.; RETURN
    ENDIF
    IF (do_bandavg) THEN
      print * , 'not implemented'; STOP 1
    ENDIF
  ENDDO

  abscrs = 0.0
  DO i = 1, nlayers
       thet = tsgrid(i) - zerok
       abscrs(1:nlamda,i) = tmpabs(1,1:nlamda) + tmpabs(2,1:nlamda) * thet + tmpabs(3,1:nlamda)*thet *thet
  ENDDO

  ! with perterbated slit function
  IF (which_pslwf  == 1) THEN 
     solwinfit_save = solwinfit
     solwinfit(1:numwin, slit_idx, 1) = solwinfit_save(1:numwin,slit_idx,1) *1.001
     
     DO i = 1, o3%nt
      CALL convolf2c_i0effect(o3%wvl(1:nline),o3%crs0(i,1:nline),nline,1,scalex(1), &
                                lamda(1:nlamda), tmpabs(i,1:nlamda), nlamda)    
     ENDDO

     dadp = 0.0
     DO i = 1, nlayers
       thet = tsgrid(i) - zerok
       dadp(:,i) = tmpabs(1,:) + tmpabs(2,:) * thet + tmpabs(3,:)*thet *thet
     ENDDO
     dadp = (dadp - abscrs)
     solwinfit = solwinfit_save
     fidx =1
     Do i = 1, numwin
       lidx = fidx + nradpix_sav(i) -1
       dadp(fidx:lidx, 1:nlayers) = dadp(fidx:lidx, 1:nlayers)/(solwinfit(i,slit_idx,1)*0.001)
       fidx = lidx + 1               
     ENDDO
  ELSE IF (which_pslwf == 2) THEN 
     IF (slit_idx == hwe_idx) THEN 
         do_dsdw = .true.  ; do_dsdk = .false. 
     ELSE IF (slit_idx == spk_idx) THEN 
         do_dsdw = .false.  ; do_dsdk = .true. 
     ENDIF
     DO i = 1, o3%nt
      CALL convolf2c_i0effect ( o3%wvl(1:nline), o3%crs0(i,1:nline),nline,1,scalex, &
                                lamda(1:nlamda), tmpabs(i,1:nlamda), nlamda)    
     ENDDO
     DO i = 1, nlayers
       thet = tsgrid(i) - zerok
       dadp(:,i) = tmpabs(1,:) + tmpabs(2,:) * thet + tmpabs(3,:)*thet *thet
     ENDDO
     do_dsdw = .false. ; do_dsdk = .false.
  ENDIF 

  dadp(1:nlamda, 1:nlayers) =dadp(1:nlamda, 1:nlayers)/abscrs(1:nlamda, 1:nlayers) 
  RETURN
 END SUBROUTINE getabs_pslwf


  SUBROUTINE geto4_crs  (lamda, nlsav, nlamda, nlayers, tsgrid, & 
                         abscrs, problems)

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
  REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: abscrs
  LOGICAL, INTENT(OUT) :: problems
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  INTEGER, PARAMETER :: maxline  = max_spec_pts      ! # of wavelengths
  INTEGER, PARAMETER :: maxt = 3  ! most 3 except for o4
  INTEGER            :: fidx, lidx, i, j, errstat, nline, nt, ntemp
  REAL (KIND=dp)     :: scalex, normc
  LOGICAL            :: dods, dodt
  REAL (KIND=dp), DIMENSION(maxt, nlsav)   :: savabs, savabs_d1
  REAL (KIND=dp), DIMENSION(maxt, nlamda)  :: tmpabs, tmpabs_d1
  REAL (KIND=dp), DIMENSION(nlamda, nlayers):: dads, dadt
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  TYPE(txcrs_set), SAVE :: o4
  LOGICAL,SAVE :: first=.true.
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'geto4_crs'
  dods=.false. ; dodt = .false.
  problems = .FALSE.
  ! load origianl spectrum
  IF (first) THEN
    CALL allocate_txcrs(o4)
    CALL read_txcrs (winwav_min, winwav_max, o2o2_idx, maxt,o4, problems)    
    o4%crs = o4%crs0   
    IF (o4%slitconv) o4crs_convl = .true.     
    first = .FALSE.
  ENDIF
  nline = o4%nw ; nt = o4%nt ; normc = o4%normc
  !-----------------------------------------------------------------------------------------
  ! convolution     
  !-----------------------------------------------------------------------------------------
  IF (o4crs_convl .AND. o4%slitconv) THEN  
    scalex = 1 
    o4%crs(1:nt, 1:nline) = o4%crs0(1:nt, 1:nline)
    DO i = 1, nt 
         CALL convol_i0effect(o4%wvl(1:nline), o4%crs(i, 1:nline), nline, &
              scalex,  errstat)
         IF ( errstat /= 0 ) THEN
           WRITE(*, *) modulename//': Error in Correct I0 Effect!!!'
           problems = .TRUE.; RETURN
         ENDIF
    ENDDO
    o4crs_convl = .FALSE.
  ENDIF

  !------------------------------------------------------------------------
  ! interpolation
  !-------------------------------------------------------------------------
  savabs = 0.0 ; savabs_d1 = 0.0  ! onto refwvl 
  tmpabs = 0.0 ; tmpabs_d1 = 0.0  ! onto re-sampled wavelengths 

  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= o4%wvl(1))))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= o4%wvl(o4%nw))))

  DO i = 1, nt
    CALL BSPLINE(o4%wvl(1:nline), o4%crs(i, 1:nline), nline, &
          lamda(fidx:lidx), savabs(i, fidx:lidx),lidx-fidx+1, errstat)
    IF (errstat < 0) THEN
          WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
          problems = .TRUE.; RETURN
    ENDIF
    IF (do_bandavg) THEN
      CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)
      IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
         WRITE(*, *) modulename//':SO2 crs Averaging Error: ', nlsav, nlamda, ntemp
         problems = .TRUE.; RETURN
      ENDIF
      tmpabs(i, :) = savabs(i, 1:nlamda)
    ELSE
      tmpabs(i, :) = savabs(i, 1:nlamda)
    ENDIF
  ENDDO

  abscrs = 0.0
  ! second tmpabs should be changed to crsshi if do_shi is true.
  problems = calc_crsz (tmpabs(1:nt,fidx:lidx),tmpabs(1:nt, fidx:lidx), &
               nt,lidx-fidx+1,o4%tdepend,o4%ts(1:nt),tsgrid(1:nlayers), nlayers, &
               dods,dodt,abscrs(fidx:lidx,1:nlayers), & 
               dadsz=dads(1:nlamda, 1:nlayers), dadtz=dadt(1:nlamda, 1:nlayers) )
  IF (dods) dads = dads / abscrs   ! get relative sensitivty to shift
  IF (dodt) dadt = dadt / abscrs   ! get relative sensitivity to T
  abscrs = abscrs * normc
  RETURN  
  END SUBROUTINE geto4_crs

  SUBROUTINE geto2_crs_hitran(lamda, nlsav, nlamda, nlayers, tsgrid, psgrid, & 
                              abscrs, do_convl)
  IMPLICIT NONE
  
  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN)                                     :: nlamda, nlsav, nlayers
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)            :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers)          :: tsgrid, psgrid
  LOGICAL, INTENT(IN) :: do_convl
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: abscrs
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  INTEGER, PARAMETER :: maxline  = max_spec_pts      ! # of wavelengths
  INTEGER            :: fidx, lidx, i, j, errstat, nline, nz, ntemp
  REAL (KIND=dp)     :: scalex
  LOGICAL            :: do_i0corr, problems
  REAL (KIND=dp), DIMENSION(nlayers, nlsav)   :: savabs
  REAL (KIND=dp), DIMENSION(nlayers, nlamda)  :: tmpabs
  REAL (KIND=dp), DIMENSION(maxline, nlayers) :: crsz
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  INTEGER, SAVE :: nw_lut, nt_lut, np_lut
  REAL (KIND=dp), SAVE ::wvl_lut(maxline), t_lut (maxt_ht), p_lut(maxp_ht)
  LOGICAL, SAVE :: first = .TRUE.
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'geto2_crs'
 
  !------------ 
  ! initialize
  !-----------
  do_i0corr = .FALSE.
  
  WRITE(www_lun, '(A)')  'scalex should be reconsidered'
  ! load origianl spectrum
  IF (first) THEN
     CALL read_hitran_lut(o2_idx, winwav_min, winwav_max, & 
          nw_lut, nt_lut, np_lut, wvl_lut, t_lut, p_lut)
     first = .FALSE.
     o2crs_convl = .true.
     o2_lut = o2_lut0
  ENDIF
    
  nline = nw_lut
  !-----------------------------------------------------------
  ! convolution     
  !-----------------------------------------------------------
  IF (do_convl .and. o2crs_convl) THEN  
    o2_lut = o2_lut0
    scalex = 6.0E24! ~600 DU
    ! Perform solar i0 effect on ozone cross-section (no need to convolve)
    DO i = 1, nt_lut
    DO j = 1, np_lut
        IF (do_i0corr) THEN 
          CALL convol_i0effect(wvl_lut(1:nline),o2_lut(1:nline,i,j), nline, &
              scalex,  errstat)
        ELSE
          CALL convol(wvl_lut(1:nline),o2_lut(1:nline,i,j), nline, errstat) 
        ENDIF
        IF ( errstat /= 0 ) THEN
          WRITE(*, *) modulename//': Error in Correct I0 Effect!!!'
          problems = .TRUE.; RETURN
        ENDIF
    ENDDO
    ENDDO
    o2crs_convl = .false.
  ENDIF

  !---------------------------------------------------------------------
  ! interpolate onto T and P
  !----------------------------------------------------------------------
  nz = nlayers
  CALL calc_hitran_crsz (o2_idx,nw_lut, nt_lut, np_lut, & 
               t_lut, p_lut,& 
               nz, tsgrid, psgrid, crsz(1:nw_lut,1:nz),problems)

  !------------------------------------------------------------------------
  ! interpolation onto W
  !-------------------------------------------------------------------------
  savabs = 0.0 
  tmpabs = 0.0 

  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= wvl_lut(1))))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= wvl_lut(nline))))

  DO i = 1, nz
    CALL BSPLINE(wvl_lut(1:nline), crsz(1:nline,i), nline, &
          lamda(fidx:lidx), savabs(i, fidx:lidx),lidx-fidx+1, errstat)
    IF (errstat < 0) THEN
          WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
          problems = .TRUE.; RETURN
    ENDIF
    IF (do_bandavg) THEN
      CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)
      IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
         WRITE(*, *) modulename//':o2 crs Averaging Error: ', nlsav, nlamda, ntemp
         problems = .TRUE.; RETURN
      ENDIF
      tmpabs(i, :) = savabs(i, 1:nlamda)
    ELSE
      tmpabs(i, :) = savabs(i, 1:nlamda)
    ENDIF
    abscrs(:, i) = tmpabs(i,:)
  ENDDO
    
  RETURN  
  END SUBROUTINE geto2_crs_hitran

  SUBROUTINE geth2o_crs_hitran(lamda, nlsav, nlamda, nlayers, tsgrid, psgrid, & 
                              abscrs, do_convl)
  IMPLICIT NONE
  
  !----------------------------------------------
  ! Input variables
  !------------------------------------------------
  INTEGER, INTENT(IN)                                     :: nlamda, nlsav, nlayers
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)            :: lamda
  REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers)          :: tsgrid, psgrid
  LOGICAL, INTENT(IN) :: do_convl
  !------------------------------------------------- 
  ! Output variables
  !-------------------------------------------------
  REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: abscrs
  !-------------------------------------------------
  ! Local variables
  !-------------------------------------------------
  INTEGER, PARAMETER :: maxline  = max_spec_pts      ! # of wavelengths
  INTEGER            :: fidx, lidx, i, j, errstat, nline, nz, ntemp
  REAL (KIND=dp)     :: scalex
  LOGICAL            :: do_i0corr, problems
  REAL (KIND=dp), DIMENSION(nlayers, nlsav)   :: savabs
  REAL (KIND=dp), DIMENSION(nlayers, nlamda)  :: tmpabs
  REAL (KIND=dp), DIMENSION(maxline, nlayers) :: crsz
  !-------------------------------------------------
  ! Save variables
  !-------------------------------------------------
  INTEGER, SAVE :: nw_lut, nt_lut, np_lut
  REAL (KIND=dp), SAVE ::wvl_lut(maxline), t_lut (maxt_ht), p_lut(maxp_ht)
  LOGICAL, SAVE :: first = .true.
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=10), PARAMETER    :: modulename = 'geth2o_crs'
 

  !------------
  ! initialize
  !-----------
  do_i0corr = .FALSE.

  WRITE(www_lun, '(A)')  'scalex should be reconsidered'
  ! load origianl spectrum
  IF (first) THEN
     CALL read_hitran_lut(h2o_idx, winwav_min, winwav_max, &
          nw_lut, nt_lut, np_lut, wvl_lut, t_lut, p_lut)
     first = .FALSE.
     h2ocrs_convl = .true.
     h2o_lut = h2o_lut0
  ENDIF

  nline = nw_lut
  !-----------------------------------------------------------
  ! convolution
  !-----------------------------------------------------------
  IF (do_convl .and. h2ocrs_convl) THEN
    h2o_lut = h2o_lut0
    scalex = 6.0E24! ~600 DU
    ! Perform solar i0 effect on ozone cross-section (no need to convolve)
    DO i = 1, nt_lut
    DO j = 1, np_lut
        IF (do_i0corr) THEN
          CALL convol_i0effect(wvl_lut(1:nline),h2o_lut(1:nline,i,j), nline, &
              scalex,  errstat)
        ELSE
          CALL convol(wvl_lut(1:nline),h2o_lut(1:nline,i,j), nline, errstat)
        ENDIF
        IF ( errstat /= 0 ) THEN
          WRITE(*, *) modulename//': Error in Correct I0 Effect!!!'
          problems = .TRUE.; RETURN
        ENDIF
    ENDDO
    ENDDO
    h2ocrs_convl = .false.
  ENDIF
 
  !---------------------------------------------------------------------
  ! interpolate onto T and P
  !----------------------------------------------------------------------
  nz = nlayers
  CALL calc_hitran_crsz (h2o_idx,nw_lut, nt_lut, np_lut, &
               t_lut, p_lut,&
               nz, tsgrid, psgrid, crsz(1:nw_lut,1:nz),problems)

  !------------------------------------------------------------------------
  ! interpolation onto W
  !------------------------------------------------------------------------- 
  savabs = 0.0 
  tmpabs = 0.0 

  fidx = MINVAL(MINLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) >= wvl_lut(1))))
  lidx = MINVAL(MAXLOC(lamda(1:nlsav), MASK=(lamda(1:nlsav) <= wvl_lut(nline))))

  DO i = 1, nz
    CALL BSPLINE(wvl_lut(1:nline), crsz(1:nline,i), nline, &
          lamda(fidx:lidx), savabs(i, fidx:lidx),lidx-fidx+1, errstat)
    IF (errstat < 0) THEN
          WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
          problems = .TRUE.; RETURN
    ENDIF
    IF (do_bandavg) THEN
      CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)
      IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
         WRITE(*, *) modulename//':h2o crs Averaging Error: ', nlsav, nlamda, ntemp
         problems = .TRUE.; RETURN
      ENDIF
      tmpabs(i, :) = savabs(i, 1:nlamda)
    ENDIF
    abscrs(:, i) = tmpabs(i,:)
  ENDDO
    
  RETURN  
  END SUBROUTINE geth2o_crs_hitran

  SUBROUTINE get_all_raycof_depol(nw, waves, raycof, depol)
  IMPLICIT NONE
  !------------------
  !Input/Output
  !------------------
  INTEGER, INTENT(IN)                        :: nw
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(nw), INTENT(OUT) :: raycof, depol
  !-----------------
  !Local variables
  !-----------------
  REAL (KIND=dp), DIMENSION(:), POINTER :: & 
       sig, sig2, sig2p, sig4, fk_n2, fk_o2, fking ! nw
  REAL (KIND=dp), PARAMETER     :: abod = 1.0455996d0, bbod = -341.29061d0, &
       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
  !--------------------------------------------------------------------
  !     Rayleigh coefficient
  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
  !--------------------------------------------------------------------
  
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
  REAL (KIND=dp), DIMENSION(nw) :: sig, sig2, sig2p, sig4
  REAL (KIND=dp), PARAMETER     :: abod = 1.0455996d0, bbod = -341.29061d0, &
       cbod = -0.90230850d0, dbod = 0.0027059889d0, ebod = -85.968563d0
  !---------------------------
  !     Rayleigh coefficient
  ! Using bodhaine et al, j. atm. oceanic tech. 16, 1854-1861, 1999.
  !--------------------------
  sig =    1.0d3 / waves
  sig2 =   sig * sig
  sig2p =  1.d0 / sig2
  sig4 =   sig2 * sig2
  raycof = (abod + bbod * sig2 + cbod * sig2p) &
       / (1.d0 + dbod * sig2 + ebod * sig2p) * 1.d-28

  RETURN

  END SUBROUTINE get_all_raycof 

  SUBROUTINE get_all_raycof_depol1 (nw, waves, nw1, raycof, depol, do_first,problems)
  ! it is called for every pixels, so allocation/deaollocation is not applied.
  IMPLICIT NONE

  !     Input/Output
  INTEGER, INTENT(IN)                        :: nw, nw1
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)  :: waves
  REAL (KIND=dp), DIMENSION(nw1), INTENT(OUT):: raycof, depol
  LOGICAL, INTENT(IN)                        :: do_first
  LOGICAL, INTENT(OUT)                       :: problems

  ! Local variables
  INTEGER                               :: i, errstat, ni0, ntemp
  REAL (KIND=dp)                        :: scalex
  REAL (KIND=dp), DIMENSION(nw)         :: raycof1, depol1

  INTEGER, SAVE                              :: nref
  REAL (KIND=dp), DIMENSION(:),POINTER, SAVE :: ray, dep, refwavs
  REAL (KIND=dp),                       SAVE :: rnorm, dnorm
  LOGICAL, SAVE                              :: first = .TRUE.

  problems = .FALSE.
  IF (do_first) first = .true.
  IF (first) THEN

     allocate (ray(max_spec_pts), dep(max_spec_pts), refwavs(max_spec_pts))
     ni0 = n_refspec_pts(1); nref = ni0
     refwavs(1:nref) = refspec_orig_data(1, 1:ni0, 1)

     CALL get_all_raycof_depol(nref, refwavs, ray(1:nref), dep(1:nref))
     rnorm = 1.0E-25; dnorm = 1.0E-2
     ray(1:nref) = ray(1:nref) / rnorm; dep(1:nref) = dep(1:nref) / dnorm

     scalex = 1.0  ! dummy variable here
     CALL convol_i0effect(refwavs(1:nref), ray(1:nref), nref, &
          scalex,  errstat)
     IF ( errstat /= 0 ) THEN
        WRITE(www_lun, *) 'Error in correcting I0 effect for raycof!!!'
        problems = .TRUE.; RETURN
     ENDIF
     CALL convol_i0effect(refwavs(1:nref), dep(1:nref), nref, scalex, errstat )
     IF ( errstat /= 0 ) THEN
        WRITE(www_lun, *) 'Error in correcting I0 effect for depol!!!'
        problems = .TRUE.; RETURN
     ENDIF

     first = .FALSE.
  ENDIF

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

  raycof = raycof1(1:nw1) * rnorm
  depol  = depol1(1:nw1)  * dnorm

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
  REAL (KIND=dp), DIMENSION(:),  POINTER  :: waves, sol, weights, rays, dpols !( max_spec_pts)
  REAL (KIND=dp), DIMENSION(:,:), POINTER :: ozcrs  !(3, max_spec_pts)
  REAL (KIND=dp), DIMENSION(:,:), POINTER :: gcrs   !(6, max_spec_pts)
  CHARACTER (len=maxchlen)                   :: crs_fname
  
  ! Saved variables
  LOGICAL,                      SAVE :: first = .TRUE.
  REAL (KIND=dp), DIMENSION(3), SAVE :: cozcrs
  REAL (KIND=dp), DIMENSION(6), SAVE :: cgcrs
  REAL (KIND=dp),               SAVE :: craycof, cdepol

  problems = .FALSE.
  IF (first) THEN
     allocate (waves(max_spec_pts), sol(max_spec_pts), weights(max_spec_pts))
     allocate (rays(max_spec_pts), dpols(max_spec_pts))
     allocate (ozcrs(3, max_spec_pts), gcrs(6, max_spec_pts))

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

     CALL GET_ALL_RAYCOF_DEPOL(nw, waves(1:nw), rays(1:nw), dpols(1:nw))
     craycof = SUM(rays(1:nw)   * weights(1:nw)) 
     cdepol  = SUM(dpols(1:nw)  * weights(1:nw)) 

     first = .FALSE.
     deallocate (waves, sol, weights, rays, dpols, ozcrs, gcrs)
  ENDIF
     
  abscrs(1, :) = cozcrs(1) + (ts - zerok) * cozcrs(2) + (ts - zerok) ** 2.0 * cozcrs(3)
  DO i = 1, 6
     abscrs(i+1, :) = cgcrs(i)
  ENDDO
  raycof = craycof; depol = cdepol
  

  RETURN
  END SUBROUTINE get_alb_ozcrs_ray

  SUBROUTINE get_effres_gascrs_ray (& 
   num_iter, nlsav, lamda, nlamda, nz, ts, ps,nfgas, allcol, rhos, &
   do_o3shi, o3shi,do_tmpwf, do_o3hwe,do_o3spk, &
   allcrs, dadsz, dadtz, dadp1, dadp2,raycof, depol,problems)

  USE OMSAO_indices_module, ONLY: so2_idx, so2v_idx, o2o2_idx, o2_idx, h2o_idx
  USE OMSAO_variables_module, ONLY: numwin, winlim, & 
      nradpix_sav, radwvl_sav, nradpix, database_save, database_shiwf,  &
      n_refspec_pts,  refspec_orig_data, refidx, &
      fitvar_rad, rmask_fitvar_rad, database
  USE OMSAO_parameters_module, ONLY: max_fit_pts, max_spec_pts
  USE ozprof_data_module, ONLY: mflay, do_subfit, nos, oswins, osfind, &
      use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, & 
      ngas, gasidxs,fgasidxs, fgassidxs

  IMPLICIT NONE
  ! Input/output variables
  INTEGER, INTENT (IN)                                  :: nlsav,nlamda, nz, nfgas, num_iter 
  LOGICAL, INTENT (IN)                                  :: do_o3shi, do_tmpwf,do_o3hwe, do_o3spk
  REAL (KIND=dp), DIMENSION(nlsav), INTENT (IN )        :: lamda
  REAL (KIND=dp), DIMENSION(nz), INTENT (IN )           :: ts, ps, rhos
  REAL (KIND=dp), DIMENSION(nfgas, nz), INTENT(IN)      :: allcol
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)    :: o3shi
  REAL (KIND=dp), DIMENSION(nlamda), INTENT (OUT)       :: raycof, depol
  REAL (KIND=dp), DIMENSION(nlamda, nz), INTENT(OUT) :: dadsz, dadtz, dadp1, dadp2
  REAL (KIND=dp), DIMENSION(nlamda, nfgas, nz), INTENT(OUT) :: allcrs
  LOGICAL, INTENT (OUT)                                 :: problems
  ! Local variables
  LOGICAL :: do_convl
  INTEGER :: fidx, lidx, i, j,k, npts, nfgas1, errstat
  REAL (kind=dp) :: tmp, temp, normc
  REAL (KIND=dp), DIMENSION (nlsav) :: delshi, tmpwav, delpos, gshiwf
  REAL (KIND=dp), DIMENSION (:), POINTER :: tmprefspec, tmprefwav
  !--------------------------------
  ! Save variables (max_fit_pts)
  !--------------------------------
  LOGICAL :: first=.TRUE.
  TYPE crsz_group
    REAL (KIND=dp), DIMENSION(:), POINTER :: raycof, depol
    REAL (KIND=dp), DIMENSION(:,:), POINTER  :: o3, so2, o4, o2,h2o, & 
                                           do3ds, do3dt, do3dhwe, do3dspk
  END TYPE crsz_group
  TYPE (crsz_group), SAVE :: crsz
  !-------------------------------------------------
  ! Name of this subroutine/module 
  !-------------------------------------------------
  CHARACTER(10), PARAMETER :: modulename ='get_effcrs'
  
  IF (first) THEN
     allocate (tmprefspec(max_spec_pts), tmprefwav(max_spec_pts))
     allocate (crsz%raycof(max_fit_pts), crsz%depol(max_fit_pts)) 
     allocate (crsz%o3(max_fit_pts, mflay),crsz%do3ds(max_fit_pts, mflay))
     allocate (crsz%do3dt(max_fit_pts, mflay))
     allocate (crsz%do3dhwe(max_fit_pts, mflay))
     allocate (crsz%do3dspk(max_fit_pts, mflay))
     IF (use_so2dtcrs) allocate (crsz%so2(max_fit_pts, mflay))
     IF (use_o4dtcrs)  allocate (crsz%o4(max_fit_pts, mflay))
     IF (use_o2dptcrs) allocate (crsz%o2(max_fit_pts, mflay))
     IF (use_h2odptcrs)allocate (crsz%h2o(max_fit_pts, mflay))
     first = .false. 
  ENDIF

  allcrs = 0.0
  IF (num_iter == 0) THEN
     CALL GET_ALL_RAYCOF_DEPOL1(nlsav, lamda, nlamda,& 
          crsz%raycof(1:nlamda),crsz%depol(1:nlamda), .false.,problems)
  ENDIF 
  raycof(1:nlamda) = crsz%raycof(1:nlamda) 
  depol(1:nlamda) = crsz%depol(1:nlamda)
  

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
  
  IF (num_iter == 0. .OR.  do_o3shi .OR. do_o3hwe .OR. do_o3spk) THEN 
     !IF (num_iter == 0) WRITE(*, *) modulename,' : Set up O3 absorption !!!'
     CALL geto3_crs(lamda - delshi, nlsav, nlamda, nz, ts, crsz%o3(1:nlamda,1:nz), &
      do_o3shi, do_tmpwf,crsz%do3ds(1:nlamda, 1:nz),crsz%do3dt(1:nlamda, 1:nz), problems)
      IF (problems) THEN
        WRITE(www_lun, *) modulename,' : Problems in reading O3 absorption !!!'
        RETURN
     ENDIF
     IF (do_o3hwe) THEN 
        CALL getabs_pslwf(lamda-delshi, nlsav, nlamda,nz, ts,hwe_idx,& 
                         crsz%do3dhwe(1:nlamda,1:nz), problems)
     ENDIF
     IF (do_o3spk ) THEN
        ! w.r.t slit shape factor
        CALL getabs_pslwf(lamda-delshi, nlsav, nlamda,nz, ts,spk_idx, &
                         crsz%do3dhwe(1:nlamda,1:nz), problems)
     ENDIF
  ENDIF

  ! GET T-dependent So2 cross section at instrument spectral resolution
  IF ( num_iter == 0 .AND.  use_so2dtcrs) THEN
    !WRITE(*, *) modulename,' : Set up SO2 absorption !!!'
    CALL getso2_crs(lamda, nlsav, nlamda, nz, ts, &
                    crsz%so2(1:nlamda, 1:nz), problems)
    IF (problems) THEN
       WRITE(*, *) modulename, ' : Problems in reading SO2 absorption !!!'
      RETURN
    ENDIF
  ENDIF

  ! GET T-dependent O4 cross section at instrument spectral resolution
  IF (num_iter == 0 .AND. use_o4dtcrs) THEN
    !WRITE(*, *) modulename,' : Set up O4 absorption !!!'
    CALL geto4_crs(lamda, nlsav, nlamda,  nz, ts, &
                  crsz%o4(1:nlamda, 1:nz), problems)
    IF (problems) THEN
       WRITE(*, *) modulename, ' : Problems in reading O4 absorption !!!'
      RETURN
    ENDIF
  ENDIF
  
  do_convl = .true.
  ! GET T-dependent o2 cross section at instrument spectral resolution
  IF ( num_iter == 0 .AND. use_o2dptcrs) THEN
    !WRITE(*, *) modulename,' : Set up O2 absorption !!!'
    CALL geto2_crs_hitran(lamda, nlsav, nlamda, nz,ts, ps,&
         crsz%o2(1:nlamda, 1:nz), do_convl)
    IF (problems) THEN
      WRITE(*, *) modulename, ' : Problems in getting H2O cross section!!!'
      RETURN
    ENDIF
  ENDIF

  ! GET T-dependent h2o cross section at instrument spectral resolution
  IF ( num_iter == 0 .AND. use_h2odptcrs) THEN
    !WRITE(*, *) modulename,' : Set up H2O absorption !!!'
    CALL geth2o_crs_hitran(lamda, nlsav, nlamda, nz, ts, ps, &
                           crsz%h2o(1:nlamda, 1:nz), do_convl)
    IF (problems) THEN
      WRITE(*, *) modulename, ' : Problems in getting H2O cross section!!!'
      RETURN
    ENDIF
  ENDIF

  dadsz(1:nlamda, 1:nz) = crsz%do3ds(1:nlamda, 1:nz)
  dadtz(1:nlamda, 1:nz) = crsz%do3dt(1:nlamda, 1:nz)
  dadp1(1:nlamda, 1:nz) = crsz%do3dhwe(1:nlamda, 1:nz)
  dadp2(1:nlamda, 1:nz) = crsz%do3dspk(1:nlamda, 1:nz)
  allcrs(1:nlamda, 1, 1:nz) = crsz%o3(1:nlamda, 1:nz)
  nfgas1 = 1
  DO i = 1, ngas
     IF (fgasidxs(i) > 0 ) THEN 
         nfgas1 = nfgas1 + 1
         normc = refspec_norm(gasidxs(i))
         IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. use_so2dtcrs) THEN 
           allcrs(1:nlamda, nfgas1, 1:nz) = crsz%so2(1:nlamda,1:nz)/normc
         ELSE IF (gasidxs(i) == o2o2_idx .AND. use_o4dtcrs) THEN 
           allcrs(1:nlamda, nfgas1, 1:nz) = crsz%o4(1:nlamda,1:nz)/normc
         ELSE IF (gasidxs(i) == o2_idx .AND. use_o2dptcrs) THEN 
           allcrs(1:nlamda, nfgas1, 1:nz) = crsz%o2(1:nlamda,1:nz)/normc
         ELSE IF (gasidxs(i) == h2o_idx .AND. use_h2odptcrs) THEN 
           allcrs(1:nlamda, nfgas1, 1:nz) = crsz%h2o(1:nlamda,1:nz)/normc
         ELSE 
           IF (fgassidxs(i) > 0 ) THEN 
              npts = n_refspec_pts(gasidxs(k))
              tmprefwav(1:npts) = refspec_orig_data(gasidxs(k), 1:npts, 1)
              tmprefspec(1:npts) = refspec_orig_data(gasidxs(k),1:npts, 3)
              fidx = MINVAL(MINLOC(lamda(1:nlamda), MASK=(lamda(1:nlamda) >= &
                         tmprefwav(1) + 0.1 .AND. lamda(1:nlamda) <= tmprefwav(npts)- 0.1)))
              lidx = MINVAL(MAXLOC(lamda(1:nlamda), MASK=(lamda(1:nlamda) >= &
                          tmprefwav(1) + 0.1 .AND. lamda(1:nlamda) <= tmprefwav(npts)- 0.1)))
              IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN
                temp = fitvar_rad(rmask_fitvar_rad(fgassidxs(k)))
                CALL BSPLINE1(  tmprefwav(1:npts) - temp,tmprefspec(1:npts), npts, &
                     lamda(fidx:lidx), allcrs(fidx:lidx, nfgas1, 1),gshiwf(fidx:lidx), lidx-fidx+1, errstat)
                database_shiwf(gasidxs(k), refidx(fidx:lidx)) = gshiwf(fidx:lidx)
                database(gasidxs(k), refidx(fidx:lidx)) =allcrs(fidx:lidx, nfgas, 1)
                IF (errstat < 0) THEN
                  WRITE(*, *) modulename, ' : BSPLINE error, errstat =', errstat; RETURN
                ENDIF
              ENDIF
           ELSE
              allcrs(1:nlamda, nfgas1, 1) = database_save(gasidxs(i),refidx(1:nlamda))
           ENDIF
           DO k = 2, nz
              allcrs(1:nlamda, nfgas1, k) = allcrs(1:nlamda, nfgas1, 1)
           ENDDO
         ENDIF
     ENDIF
  ENDDO
  RETURN
  END SUBROUTINE get_effres_gascrs_ray

! Prepare high resolution spectra at fine grid: solar reference, trace gas cross
! sections
!   (o3 shift and o3 temperature), Raleigh scattering coefficient,
!   depolarization factor
! O3/SO2 (use_so2dtcrs=.TRUE.) cross section: if do_tmpwf = .FALSE. and do_o3shi
! is false,
! just need to get once for each retrieval
! Other trace gas cross section: just need to get it once for all the retrievals
! if no shifts
  SUBROUTINE get_hres_gascrs_ray ( &
             num_iter, nw, waves, nz, ts, ps, nfgas,allcol, rhos, &
             do_o3shi, o3shi, do_tmpwf, & 
             allcrs, raycof, depol, pge_error_status)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY  : max_spec_pts, zerok
  USE OMSAO_variables_module, ONLY  : numwin, nradpix, refspec_orig_data,    &
       n_refspec_pts, solwinfit, which_slit, curr_rad_spec, use_redfixwav,   &
       winlim, fitvar_rad, rmask_fitvar_rad, refspec_norm
  USE OMSAO_indices_module,   ONLY : hwe_idx, wvl_idx, solar_idx, spc_idx,  &
       so2_idx, so2v_idx, o2o2_idx, o2_idx, h2o_idx, o2t2_idx, h2ot2_idx
  USE ozprof_data_module,     ONLY : mxsect,mflay, do_subfit, nos, oswins, osfind, &
       ngas, gasidxs, fgasidxs,  fgassidxs, &
       nhresp, nhresp0, hreswav0, hres_samprate, hreswav, radcidxs, ncalcp, &
       hres_gas, hres_gasshi,hres_o3, hresgabs, hresray, &
       hres_i0, hres_raycof, hres_depol, &
       hres_o3shi, hres_so2, hres_so2shi, hres_o4, hres_o4shi,  & ! could be local variables
       hres_o2, hres_h2o, hres_o2shi, hres_h2oshi, &
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
       o3crsz, o3dadsz, o3dadtz, so2crsz, so2dads, o4crsz, o4dads, &
       h2ocrsz,h2odads, o2crsz, o2dads  ! could be local variables
  USE OMSAO_errstat_module
  IMPLICIT NONE

  !-----------------------------
  ! Input/output variables
  !--------------------------------
  INTEGER, INTENT (IN)                                  :: nw, nz, nfgas, num_iter ! nw = ncalcp
  INTEGER, INTENT (OUT)                                 :: pge_error_status
  LOGICAL, INTENT (IN)                                  :: do_o3shi, do_tmpwf
  REAL (KIND=dp), DIMENSION(nw), INTENT (IN )           :: waves
  REAL (KIND=dp), DIMENSION(nz), INTENT (IN )           :: ts, ps, rhos
  REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)    :: o3shi
  REAL (KIND=dp), DIMENSION(nfgas, nz), INTENT(IN)      :: allcol
  REAL (KIND=dp), DIMENSION(nw), INTENT (OUT)           :: raycof, depol
  REAL (KIND=dp), DIMENSION(nw, nfgas, nz), INTENT(OUT) :: allcrs

  !-----------------------------
  ! Local variables
  !--------------------------------
  INTEGER :: i, j, k, fidx, lidx,  npts, idum, nfgas1, errstat, nratio, nhalf
  LOGICAL                              :: problems, do_shi, do_convl
  REAL (KIND=dp)                       :: tmp, so2sum, o4sum, o2sum,h2osum
  REAL (KIND=dp), DIMENSION(nhresp)    :: delshi, tmpwav, delpos ! (nhresp)
!  REAL (KIND=dp), DIMENSION(:, :), ALLOCATABLE ::so2dadsz, o4dadsz, o2dadsz, h2odadsz
!  REAL (KIND=dp), DIMENSION(:, :), ALLOCATABLE ::so2dadtz, o4dadtz, o2dadtz, h2odadtz
  REAL (KIND=dp), DIMENSION(nhresp, nz) ::so2dadsz, o4dadsz, o2dadsz, h2odadsz
  REAL (KIND=dp), DIMENSION(nhresp, nz) ::so2dadtz, o4dadtz, o2dadtz, h2odadtz

  ! Save original O3/SO2 cross sections
  LOGICAL, SAVE :: do_so2shi, do_o4shi, do_o2shi, do_h2oshi, first = .TRUE.
  LOGICAL, SAVE :: do_so2tmp=.false., do_o4tmp=.false., do_o2tmp=.false.,do_h2otmp=.false.
  INTEGER, SAVE :: so2sfidx, so2vsfidx, o4sfidx, o2sfidx, h2osfidx
  TYPE(txcrs_set), SAVE  :: o3, so2, o4
 ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=19), PARAMETER :: modulename = 'GET_HRES_GASCRS_RAY'

  pge_error_status = pge_errstat_ok
   !=====================================================
   ! Load cross section at original grids
   !=====================================================
   IF (first) THEN
     ! Obtain high resolution solar reference spectra
     npts = n_refspec_pts(solar_idx)
     CALL interpolation (npts, refspec_orig_data(solar_idx,1:npts,wvl_idx), &
          refspec_orig_data(solar_idx,1:npts,spc_idx), nhresp0,hreswav0(1:nhresp0),  &
          hres_i0(1:nhresp0), errstat)
     IF ( errstat > pge_errstat_warning ) THEN
        pge_error_status = pge_errstat_error; RETURN
     ENDIF

     nratio = hres_samprate / 0.01
     nhalf =  nratio / 2

     j = 1
     DO i = nhalf + 1, nhresp0 - nhalf, nratio
        fidx = i - nhalf; lidx = i - nhalf + nratio - 1
        hres_i0(j) = SUM(hres_i0(fidx:lidx)) / REAL(nratio)
        j = j + 1
     ENDDO

     ! Obtain high resolution rayleigh scattering coefficients and
     ! depolarization factor
     CALL GET_ALL_RAYCOF_DEPOL(nhresp, hreswav(1:nhresp), hres_raycof(1:nhresp),hres_depol(1:nhresp))

     ! Obtain high resolution cross sections of other trace gases (except for  O3,  SO2, O4)
     do_so2shi = .FALSE. ; do_o4shi = .FALSE. ; do_o2shi = .FALSE.; do_h2oshi = .FALSE.
     so2sfidx = 0; so2vsfidx = 0 ; o4sfidx=0 ; o2sfidx = 0; h2osfidx = 0

     hres_gas(1:ngas, 1:nhresp) = 0.0D0
     DO i = 1, ngas
        IF (fgasidxs(i) > 0 ) THEN
           ! find indices for shift
           IF (gasidxs(i) == so2_idx)  so2sfidx  = fgassidxs(i)
           IF (gasidxs(i) == so2v_idx) so2vsfidx = fgassidxs(i)
           IF (gasidxs(i) == o2o2_idx) o4sfidx   = fgassidxs(i)
           IF (gasidxs(i) == o2_idx)   o2sfidx   = fgassidxs(i)
           IF (gasidxs(i) == h2o_idx)  h2osfidx  = fgassidxs(i)
           
           IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. &
               fgassidxs(i) > 0) do_so2shi = .TRUE.
           IF (gasidxs(i) == o2o2_idx .AND. fgassidxs(i) > 0) do_o4shi = .TRUE.
           IF ((gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx) .AND. fgassidxs(i) > 0) do_o2shi = .TRUE.
           IF ((gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx) .AND. fgassidxs(i) > 0) do_h2oshi = .TRUE.

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
              CALL BSPLINE(refspec_orig_data(gasidxs(i), 1:npts, 1), &
                   refspec_orig_data(gasidxs(i), 1:npts, 2), npts, hreswav(fidx:lidx), &
                   hres_gas(i, fidx:lidx), lidx - fidx + 1, errstat)

              IF (errstat < 0) THEN
                 WRITE(*, *) modulename, ' : BSPLINE2 error, errstat = ',errstat
                 pge_error_status = pge_errstat_error; RETURN
              ENDIF
           ENDIF
        ENDIF
     ENDDO

     CALL allocate_txcrs (o3)
     IF (use_so2dtcrs) THEN 
        CALL allocate_txcrs (so2)
        !allocate (so2dadsz(nhresp, mflay))
        !allocate (so2dadtz(nhresp, mflay))
     ENDIF
     IF (use_o4dtcrs) THEN 
        CALL allocate_txcrs (o4)
        !allocate (o4dadsz(nhresp, mflay))
        !allocate (o4dadtz(nhresp, mflay))
     ENDIF
     IF (use_o2dptcrs) THEN 
        !allocate (o2dadsz(nhresp, mflay))
        !allocate (o2dadtz(nhresp, mflay))
     ENDIF
     IF (use_h2odptcrs) THEN 
        !allocate (h2odadsz(nhresp, mflay))
        !allocate (h2odadtz(nhresp, mflay))
     ENDIF 

     ! Obtain original O3 cross section (quadratic or individual T, shift, T-depen)
     CALL read_txcrs(winwav_min, winwav_max,o3_t1_idx,mxsect,o3,problems)
     IF (problems) THEN
        WRITE(*, *) modulename, ' : Error in reading ozone cross sections!!!'
        pge_error_status = pge_errstat_error; RETURN
     ENDIF

     IF (.NOT. do_o3shi) THEN
        DO i = 1, o3%nt
           CALL BSPLINE(o3%wvl(1:o3%nw), o3%crs0(i,1:o3%nw),o3%nw,& 
                hreswav(1:nhresp), hres_o3(i, 1:nhresp), nhresp, errstat)
           IF (errstat < 0) THEN
              WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
              pge_error_status = pge_errstat_error; RETURN
           ENDIF
        ENDDO
     ENDIF

     ! Obtain original SO2 cross section (quadratic or individual T, shift)
     IF (use_so2dtcrs) THEN
        CALL read_txcrs(winwav_min, winwav_max,so2_idx, mxsect,so2,problems)
        IF (problems) THEN
           WRITE(*, *) modulename, ' : Error in reading SO2 cross sections!!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF

        IF (.NOT. do_so2shi) THEN
           DO i = 1, so2%nt
                   fidx=MINVAL(MINLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) > so2%wvl(1) )))
                   lidx=MINVAL(MAXLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) < so2%wvl(so2%nw) )))
              CALL BSPLINE(so2%wvl(1:so2%nw), so2%crs0(i, 1:so2%nw),so2%nw, & 
                   hreswav(fidx:lidx),hres_so2(i, fidx:lidx), lidx-fidx+1, errstat)
              IF (errstat < 0) THEN
                 WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
                 pge_error_status = pge_errstat_error; RETURN
              ENDIF
           ENDDO
        ENDIF
     ENDIF
     IF (use_o4dtcrs) THEN
        CALL read_txcrs(winwav_min, winwav_max,o2o2_idx, mxsect,o4,problems)
        IF (problems) THEN
           WRITE(*, *) modulename, ' : Error in reading SO2 cross sections!!!'
           pge_error_status = pge_errstat_error; RETURN
        ENDIF

        IF (.NOT. do_o4shi) THEN
           DO i = 1, o4%nt
              fidx=MINVAL(MINLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) > o4%wvl(1) )))
              lidx=MINVAL(MAXLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) < o4%wvl(o4%nw) )))
              CALL BSPLINE(o4%wvl(1:o4%nw),o4%crs0(i,1:o4%nw),o4%nw, & 
                   hreswav(fidx:lidx),hres_o4(i, fidx:lidx), lidx-fidx+1, errstat)
              IF (errstat < 0) THEN
                 WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
                 pge_error_status = pge_errstat_error; RETURN
              ENDIF
           ENDDO
        ENDIF
     ENDIF

    first = .FALSE.
    ENDIF

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
          CALL BSPLINE2(o3%wvl(1:o3%nw), o3%crs0(i,1:o3%nw),o3%nw, do_o3shi,& 
           hreswav(1:nhresp)-delshi, hres_o3(i, 1:nhresp), hres_o3shi(i, 1:nhresp), nhresp, errstat)
          IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
           pge_error_status = pge_errstat_error; RETURN
          ENDIF
        ENDDO
      ENDIF

    ! Get ozone cross section at each layer
    IF (num_iter == 0 .OR. do_o3shi .OR. do_tmpwf) THEN
      problems = calc_crsz (hres_o3(1:o3%nt, 1:nhresp), hres_o3shi(1:o3%nt, 1:nhresp), &
                 o3%nt, nhresp, o3%tdepend, o3%ts(1:o3%nt), ts(1:nz), nz, &
                 do_o3shi, do_tmpwf, o3crsz(1:nhresp, 1:nz), &
                 dadsz=o3dadsz(1:nhresp, 1:nz),dadtz=o3dadtz(1:nhresp, 1:nz)) 
       IF (do_tmpwf) THEN
           o3dadtz(1:nhresp, 1:nz) = o3dadtz(1:nhresp, 1:nz) / o3crsz(1:nhresp,1:nz)
       ENDIF
       IF (do_o3shi) THEN
          o3dadsz(1:nhresp, 1:nz) = o3dadsz(1:nhresp, 1:nz) / o3crsz(1:nhresp,1:nz)
       ENDIF
       o3crsz(1:nhresp, 1:nz) = o3crsz(1:nhresp, 1:nz) * o3%normc
   ENDIF

   ! Saved O3 crs variables
   allcrs(1:ncalcp, 1, 1:nz) = o3crsz(radcidxs(1:ncalcp), 1:nz)
   DO i = 1, nz
       hresgabs(1:nhresp, i) = o3crsz(1:nhresp, i) * allcol(1, i)
   ENDDO

   ! Get SO2 cross sections
   IF (use_so2dtcrs) THEN
     fidx=MINVAL(MINLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) >so2%wvl(1) )))
     lidx=MINVAL(MAXLOC( hreswav(1:nhresp),MASK=(hreswav(1:nhresp) <so2%wvl(so2%nw) )))
     IF (do_so2shi) THEN
       DO i = 1, so2%nt
         idum = MAX(so2sfidx, so2vsfidx)
         tmp = fitvar_rad(rmask_fitvar_rad(idum))
         CALL BSPLINE2(so2%wvl(1:so2%nw), so2%crs0(i,1:so2%nw),so2%nw,do_so2shi, &
                 hreswav(fidx:lidx) - tmp, hres_so2(i, fidx:lidx), hres_so2shi(i,fidx:lidx), &
                 lidx-fidx+1, errstat)
         IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
           pge_error_status = pge_errstat_error; RETURN
         ENDIF
       ENDDO
     ENDIF
     IF (num_iter == 0 .OR. do_so2shi .OR. do_tmpwf) THEN
       so2crsz = 0.0; so2dadsz=0.0
       problems = calc_crsz (hres_so2(1:so2%nt, fidx:lidx), hres_so2shi(1:so2%nt,fidx:lidx), &
       so2%nt,lidx-fidx+1, so2%tdepend, so2%ts(1:so2%nt), ts(1:nz), nz, do_so2shi,do_so2tmp, &
       so2crsz(fidx:lidx, 1:nz), dadsz=so2dadsz(fidx:lidx, 1:nz),dadtz=so2dadtz(fidx:lidx, 1:nz)) 
       IF (do_so2shi) THEN  
         so2dadsz(1:nhresp, 1:nz) = so2dadsz(1:nhresp, 1:nz) /so2crsz(1:nhresp, 1:nz)
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
         idum = o4sfidx
         tmp = fitvar_rad(rmask_fitvar_rad(idum))
         CALL BSPLINE2(o4%wvl(1:o4%nw), o4%crs0(i,1:o4%nw),o4%nw,do_o4shi, &
                 hreswav(fidx:lidx) - tmp, hres_o4(i, fidx:lidx),hres_o4shi(i,fidx:lidx), &
                 lidx-fidx+1, errstat)
         IF (errstat < 0) THEN
           WRITE(*, *) modulename, ': BSPLINE2 error, errstat = ', errstat
           pge_error_status = pge_errstat_error; RETURN
         ENDIF
       ENDDO
     ENDIF
     IF (num_iter == 0 .OR. do_o4shi .OR. do_tmpwf) THEN
       o4crsz=0.0 ; o4dadsz = 0.0
       problems = calc_crsz (hres_o4(1:o4%nt, fidx:lidx),hres_o4shi(1:o4%nt,fidx:lidx), &
           o4%nt,lidx-fidx+1, o4%tdepend, o4%ts(1:o4%nt), ts(1:nz), nz,do_o4shi,do_o4tmp, &
           o4crsz(fidx:lidx, 1:nz), dadsz=o4dadsz(fidx:lidx, 1:nz),dadtz=o4dadtz(fidx:lidx, 1:nz))
       IF (do_o4shi) THEN              
         o4dadsz(fidx:lidx, 1:nz) =o4dadsz(fidx:lidx, 1:nz) /o4crsz(fidx:lidx, 1:nz)
       ENDIF
       o4crsz(1:nhresp, 1:nz) = o4crsz(1:nhresp, 1:nz) * o4%normc
     ENDIF
   ENDIF 

   ! Get O2 cross section
   IF (use_o2dptcrs) THEN 
     IF (num_iter == 0) THEN  
       do_convl = .false.
       o2crsz=0.0 ;o2dadsz=0.0
       CALL geto2_crs_hitran (hreswav(1:nhresp), nhresp, nhresp, & 
             nz,ts, ps, o2crsz(1:nhresp, 1:nz), do_convl)
     ENDIF
   ENDIF  

   ! Get h2o cross section
   IF (use_h2odptcrs) THEN 
     IF (num_iter == 0) THEN  
       do_convl = .false.
       h2ocrsz=0.0 ;h2odadsz=0.0
       CALL geth2o_crs_hitran (hreswav(1:nhresp), nhresp, nhresp, & 
             nz,ts, ps, h2ocrsz(1:nhresp, 1:nz), do_convl)
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
               pge_error_status = pge_errstat_error; RETURN
             ENDIF
           ENDIF
         ENDIF
         DO j = 1, nz
           allcrs(1:ncalcp, nfgas1, j) = hres_gas(i, radcidxs(1:ncalcp))
           hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + hres_gas(i,1:nhresp) * allcol(nfgas1, j)
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
  
  SUBROUTINE read_txcrs (minw, maxw, gas_idx, maxt, txcrs, problems)
    ! max_spec_pts, refdbdir, ozabs_unit
    IMPLICIT NONE
  
    ! Input variables
    INTEGER, INTENT(IN)             :: maxt, gas_idx
    REAL (KIND=dp), INTENT(IN)      :: maxw, minw
    ! Output variables
    TYPE(txcrs_set), INTENT(OUT)    :: txcrs
    LOGICAL, INTENT(OUT)            :: problems
    ! Local variables
    CHARACTER (LEN=maxchlen)        :: absfname
    INTEGER                         :: nline, i, j, errstat
    LOGICAL                         :: file_exist
    REAL (KIND=dp)                  :: tmp
    CHARACTER (LEN=14)              :: tmpchar
    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=10), PARAMETER    :: modulename = 'READ_TXCRS'
  
    CALL allocate_txcrs(txcrs)

    problems = .FALSE.; tmpchar = ' ' 
    IF (gas_idx == o3_t1_idx ) THEN   
        absfname = TRIM(ADJUSTL(ozabs_fname))
    ELSE IF (gas_idx == so2_idx .or. gas_idx == so2v_idx) THEN 
        absfname = TRIM(ADJUSTL(refdbdir)) //TRIM(ADJUSTL(so2abs_fname))
    ELSE IF (gas_idx == o2o2_idx) THEN 
        absfname = TRIM(ADJUSTL(refdbdir)) //TRIM(ADJUSTL(o4abs_fname))
    ENDIF

    INQUIRE (FILE= TRIM(ADJUSTL(absfname)), EXIST= file_exist)
    IF (file_exist) THEN
        OPEN(UNIT = ozabs_unit, file=TRIM(ADJUSTL(absfname)), status='old')
    ELSE
        write(*,*) 'read_txcrs: No file of '//TRIM(ADJUSTL(absfname)); STOP
    ENDIF

    DO WHILE (tmpchar /= 'START OF TABLE') 
       READ (ozabs_unit, '(A14)') tmpchar
    ENDDO
     
    READ (ozabs_unit, *) nline, tmp, tmp, txcrs%normc
    READ (ozabs_unit, *) txcrs%tdepend, txcrs%nt, txcrs%slitconv
    IF (txcrs%nt > maxt) THEN
       WRITE(*, *) 'Need to increase parameter mxsect!!!'
       problems = .TRUE.; CLOSE (ozabs_unit); RETURN
    ENDIF

    READ (ozabs_unit, *) txcrs%ts(1:txcrs%nt)
    IF ((txcrs%nt > 1 .AND. .NOT. txcrs%tdepend ) .AND. (MINVAL(txcrs%ts(1:txcrs%nt)) > 220. &
        .OR.  MAXVAL(txcrs%ts(1:txcrs%nt)) < 280. )) THEN
       WRITE(*, *) 'Temperature range for X-section not enough!!!'
       problems = .TRUE.; CLOSE (ozabs_unit); RETURN
    ENDIF

    j = 1
    DO i = 1, nline
       READ(ozabs_unit, *) txcrs%wvl(j), txcrs%crs0(1:txcrs%nt, j)
       IF (txcrs%wvl(j) > minw  .AND. txcrs%wvl(j) < &
           maxw)  j = j + 1
    ENDDO
    txcrs%nw = j - 1
    CLOSE (ozabs_unit)
  
    IF (txcrs%nw > max_spec_pts) THEN
       WRITE(*, *) 'Need to increase parameter max_spec_pts!!!', txcrs%nw,max_spec_pts
       problems = .TRUE.; CLOSE (ozabs_unit); RETURN
    ENDIF
    
    !txcrs%crs0(1:txcrs%nt, 1:txcrs%nw) = txcrs%crs0(1:txcrs%nt, 1:txcrs%nw)*txcrs%normc
    refspec_norm (gas_idx) = txcrs%normc

    WRITE(www_lun,*) 'N of refspectrum:', txcrs%nw ,'normc', txcrs%normc
  RETURN
  END SUBROUTINE read_txcrs

  FUNCTION calc_crsz (crs, crsshi, nt, ncrs, tdepend, ts, tz, nz, &
                      do_shi,do_tmpwf ,crsz, dadsz, dadtz) RESULT (problem)
                       
     IMPLICIT NONE
     !INPUT VARIABLES
     LOGICAL, INTENT(IN) :: tdepend, do_shi, do_tmpwf
     INTEGER, INTENT(IN) :: nt, ncrs, nz
     REAL(KIND=dp), INTENT(IN)   :: ts(nt), tz(nz)
     REAL(KIND=dp), INTENT(IN)  ::  crs(nt, ncrs), crsshi(nt, ncrs)
     !OUTPUT VARIABLES
     REAL (KIND=dp), DIMENSION(ncrs, nz) :: crsz
     REAL (KIND=dp), DIMENSION(ncrs, nz), OPTIONAL ::  dadsz, dadtz
     LOGICAL :: problem
     !local variables
     INTEGER :: i , errstat
     REAL (KIND=dp) :: thet, frac
     problem = .false.
     crsz = 0.0 ; dadsz = 0.0; dadtz = 0.0
     IF (tdepend .AND. nt == 3) THEN     ! qudratic T dependent coefficients
        DO i = 1, nz
           thet = tz(i) - zerok
           crsz(1:ncrs, i) = crs(1, 1:ncrs) + crs(2, 1:ncrs) *thet &
                              + crs(3, 1:ncrs) * thet * thet 
           IF (do_shi) THEN
             dadsz(1:ncrs, i) = crsshi(1, 1:ncrs) + crsshi(2,1:ncrs) * thet &
                                     + crsshi(3, 1:ncrs) * thet * thet 
           ENDIF
           IF (do_tmpwf) THEN
              dadtz(1:ncrs, i) = crs(2, 1:ncrs)  + 2.0 * crs(3,1:ncrs) * thet
           ENDIF
        ENDDO
     ELSE IF (tdepend .AND. nt == 2) THEN   ! linear T dependent coefficients
        DO i = 1, nz
           thet = tz(i) - zerok
           crsz(1:ncrs, i) = crs(1, 1:ncrs) + crs(2, 1:ncrs) *thet
           IF (do_shi) THEN
             dadsz(1:ncrs, i) = crsshi(1, 1:ncrs) + crsshi(2,1:ncrs) * thet 
           ENDIF
           IF (do_tmpwf) THEN
              dadtz(1:ncrs, i) = crs(2, 1:ncrs)
           ENDIF
        ENDDO
     ELSE IF (tdepend .AND. nt == 1)  THEN           ! only 1 T
        DO i = 1, nz
           crsz(1:ncrs, i) = crs(1, 1:ncrs)
           IF (do_shi) dadsz(1:ncrs, i) = crsshi(1, 1:ncrs)
           IF (do_tmpwf) dadtz(1:ncrs, i) = 0.0
        ENDDO
     ELSE IF (.NOT. tdepend .AND. nt == 2) THEN            ! have 2 T values
        DO i = 1, nz
           frac = 1.0 - (tz(i) - ts(1)) / (ts(2) - ts(1))
           crsz(1:ncrs, i) = (frac * crs(1, 1:ncrs) + (1.0 - frac) *crs(2, 1:ncrs))
           IF (do_shi) dadsz(1:ncrs, i) = (frac * crsshi(1, 1:ncrs)+ &
                (1.0 - frac) * crsshi(2, 1:ncrs))
           IF (do_tmpwf) dadtz(1:ncrs, i) = (crs(1, 1:ncrs) -crs(2, 1:ncrs)) / (ts(1) - ts(2))
        END DO
     ELSE  IF (.NOT. tdepend .AND. nt > 3) THEN  ! have more than n T
        DO i = 1, ncrs                    ! Interpolate/extrapolate over T  
           CALL INTERPOL2(ts(1:nt), crs(1:nt,i), nt, do_tmpwf, tz(1:nz), &
                crsz(i, 1:nz), dadtz(i, 1:nz), nz, errstat)
           IF (errstat < 0) THEN
              WRITE(*, *) ' INTERPOL2 error, errstat = ', errstat
              problem=.TRUE.
           ENDIF
           IF (do_shi) THEN
              CALL INTERPOL(ts(1:nt), crsshi(1:nt,i), nt, tz(1:nz), &
                   dadsz(i, 1:nz), nz, errstat)
              IF (errstat < 0) THEN
                 WRITE(*, *) 'INTERPOL error, errstat = ', errstat
               problem=.TRUE.
              ENDIF
           ENDIF
        ENDDO
     ELSE
        WRITE(*, *)  'Such type of ozone cross sections not implemented'
        problem = .TRUE.
     ENDIF
     RETURN
  END FUNCTION calc_crsz

  SUBROUTINE read_hitran_lut (gas_idx, win_min, win_max, & 
                              nw, nt, np, wvl, ts, ps)

    ! ORIGINAL CODE : GC_xsections_module.f90 (M. Chris ?)
    USE OMSAO_indices_module, ONLY: h2o_idx, h2ot2_idx, o2_idx, o2t2_idx
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    !----------------------------------
    !Input
    !---------------------------------
    INTEGER, INTENT(IN) :: gas_idx
    REAL (KIND=dp) :: win_min, win_max
    !---------------------------------
    !Output 
    !---------------------------------
    INTEGER, INTENT(OUT) :: nw, nt, np
    REAL (KIND=dp), INTENT(OUT), DIMENSION (max_spec_pts) :: wvl
    REAL (KIND=dp), INTENT(OUT), DIMENSION (maxt_ht) :: ts
    REAL (KIND=dp), INTENT(OUT), DIMENSION (maxp_ht) :: ps
    !---------------------------------
    !Local variables
    !---------------------------------
    ! helper for reading nc file
    CHARACTER (LEN=100) :: filename
    INTEGER :: ncid, rcode, var_id
    INTEGER :: posdim, pdim, Tdim
    CHARACTER(len=maxchlen) :: message
    CHARACTER(len=31) :: dimname
    LOGICAL :: file_exist, fail
    ! helper for LUT interpolation
    INTEGER :: i, j, errstat
    INTEGER :: tidx, pidx, fidx, lidx
    REAL (KIND=dp) :: tfrac, pfrac
    !-------------------------------------------------------------------
    ! LUT variables
    !--------------------------------------------------------------------
    INTEGER :: nwvl_lut, np_lut, nT_lut
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)  :: T_lut, p_lut, wvl_lut
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)  :: xs_lut
    !------------------------------------------------------------------

    fail = .false.
    IF (gas_idx == h2o_idx .or. gas_idx == h2ot2_idx) THEN 
        filename = ADJUSTL(TRIM(refdbdir))//ADJUSTL(TRIM(h2oabs_fname))
    ELSE IF (gas_idx == o2_idx .or. gas_idx == o2t2_idx) THEN 
        filename = ADJUSTL(TRIM(refdbdir))//ADJUSTL(TRIM(o2abs_fname))
    ELSE
        write(*,*) 'this gas cross section is not provided from hitran'
    ENDIF

    INQUIRE (FILE= TRIM(ADJUSTL(filename)), EXIST= file_exist)
    IF (.not. file_exist) THEN
        write(*,*) "No exsit:"//filename ; stop
    ENDIF 

    ! ================================================================
    ! Open netCDF file and allocate arrays
    ! ================================================================
    ! Open file in read mode
    ncid = ncopn(trim(adjustl(filename)), nf_Nowrite, rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in read_xy_nc_prof: ncopn failed' ; stop
    endif
   ! Get the wavelength dimension
    posdim = ncdid(ncid, 'npos', RCODE)
    if (rcode  .eq. -1 ) then
       message =  ' error in read_xy_nc_prof: ncdid failed(npos)' ;stop
    endif

    ! Read wavelength dimension
    call ncdinq(ncid, posdim, dimname, nwvl_lut, rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in read_xy_nc_prof: ncdinq failed(npos)' ;stop
    endif

    ! Get the temperature dimension
    Tdim = ncdid(ncid, 'nT', RCODE)
    if (rcode  .eq. -1 ) then
       message =  ' error in read_xy_nc_prof: ncdid failed(nT)' ;stop
    endif
   ! Read temperature  dimension
    call ncdinq(ncid, Tdim, dimname, nT_lut, rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in read_xy_nc_prof: ncdinq failed(nT)' ; stop
    endif

    ! Get the pressure dimension
    pdim = ncdid(ncid, 'np', RCODE)
    if (rcode  .eq. -1 ) then
       message =  ' error in read_xy_nc_prof: ncdid failed(np)' ; stop
    endif

    ! Read the pressure  dimension
    call ncdinq(ncid, pdim, dimname, np_lut, rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in read_xy_nc_prof: ncdinq failed(np)' ; stop
    endif

    ! Allocate arrays
    ALLOCATE(wvl_lut(nwvl_lut))
    ALLOCATE(t_lut(nt_lut))
    ALLOCATE(p_lut(np_lut))
    ALLOCATE(xs_lut(nwvl_lut,nt_lut,np_lut))
    ! ================================================================
    ! Read the variables
    ! ================================================================

    ! ---------------------
    ! Read wavelength grid
    ! ---------------------

    var_id = ncvid(ncid, 'Wavelength', rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvid failed(Wavelength)'
       fail = .true.; return
    endif
    call ncvgt(ncid, var_id, (/1/), (/nwvl_lut/), wvl_lut, rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvgt failed(Longitude)'
       fail = .true.; return
    endif

    ! ----------------------
    ! Read temperature grid
    ! ----------------------

    var_id = ncvid(ncid, 'Temperature', rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvid failed(Temperature)'
       fail = .true.; return
    endif
    call ncvgt(ncid, var_id, (/1/), (/nT_lut/), T_lut, rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvgt failed(Longitude)'
       fail = .true.; return
    endif

    ! ----------------------
    ! Read pressure grid
    ! ----------------------

    var_id = ncvid(ncid, 'Pressure', rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvid failed(Pressure)'
       fail = .true.; return
    endif
    call ncvgt(ncid, var_id, (/1/), (/np_lut/), p_lut, rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvgt failed(Pressure)'
       fail = .true.; return
    endif

    ! ----------------------
    ! Read cross sections
    ! ----------------------

    var_id = ncvid(ncid, 'CrossSection', rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvid failed(Pressure)'
       fail = .true.; return
    endif
    call ncvgt(ncid, var_id, (/1,1,1/), (/nwvl_lut,nT_lut,np_lut/), xs_lut,rcode)
    if (rcode  .eq. -1 ) then
       message =  ' error in netcdf_rd_dim: ncvgt failed(Pressure)'
       fail = .true.; return
    endif

    ! ----------
    ! Close file
    ! ----------

    call ncclos(ncid, rcode)
    if (rcode .eq. -1) then
       message =  ' error in netcdf_rd_dim: ncclos'
       fail = .true.; return
    endif
    

    ! ================================================================
    ! Interpolate spectra
    ! ================================================================
    ! this spectra has the dependence on temperature and thence
    ! it should be updated if temperature is included as state vector     
    fidx=MINVAL(MINLOC( wvl_lut(1:nwvl_lut),MASK=(wvl_lut(1:nwvl_lut) > win_min )))
    lidx=MINVAL(MAXLOC( wvl_lut(1:nwvl_lut),MASK=(wvl_lut(1:nwvl_lut) < win_max )))
    nw = lidx - fidx + 1
    wvl (1:nw) = wvl_lut(fidx:lidx)
    nt = nt_lut 
    np = np_lut
    ps (1:np) = p_lut(1:np)
    ts (1:nt) = t_lut(1:nt)
    IF (gas_idx == o2_idx) THEN 
        allocate (o2_lut(nw, nt, np), o2_lut0(nw, nt, np))
        o2_lut0 (1:nw, 1:nt, 1:np) = xs_lut(fidx:lidx, 1:nt,1:np)
    ELSE IF (gas_idx == h2o_idx) THEN 
        allocate (h2o_lut(nw, nt, np), h2o_lut0(nw, nt, np))
        h2o_lut0 (1:nw, 1:nt, 1:np) = xs_lut(fidx:lidx, 1:nt,1:np)
    ENDIF
    DEALLOCATE (wvl_lut, t_lut, p_lut, xs_lut)
    RETURN 
  END SUBROUTINE read_hitran_lut 
 
  SUBROUTINE calc_hitran_crsz (gas_idx,nw_lut, nt_lut, np_lut, &
                               t_lut, p_lut,  & 
                              nz, temp,press, crsz,fail)
   IMPLICIT NONE
  !------------------------------
  ! input variables
  !------------------------------
   INTEGER, INTENT(IN) :: gas_idx,nz, nw_lut, nt_lut, np_lut
   REAL (KIND=dp), DIMENSION (nt_lut), INTENT(IN) :: t_lut
   REAL (KIND=dp), DIMENSION (np_lut), INTENT(IN) :: p_lut
   REAL (KIND=dp), DIMENSION (nz), INTENT(IN) :: press, temp
  !------------------------------
  ! output variables
  !------------------------------
   REAL (KIND=dp), DIMENSION (nw_lut, nz), INTENT(OUT) :: crsz
   LOGICAL, INTENT(OUT) :: fail
  !-----------------------------
  ! local variables
  !-----------------------------
   INTEGER :: i, pidx, tidx
   REAL (KIND=dp) :: pmin, pmax, tmin, tmax, pfrac, tfrac
   fail = .false.
  ! Get max/min p and T
    pmin = minval(p_lut)
    pmax = maxval(p_lut)
    Tmin = minval(T_lut)
    Tmax = maxval(T_lut)  
  
  ! interpolation
  DO i = 1, nz
    pidx = minval(maxloc(p_lut(1:np_lut), MASK=(p_lut(1:np_lut) < press(i))))
    IF (press(i) <= pmin) THEN
        pidx = 1 ; pfrac = 1.0
    ENDIF
    IF (press(i) >= pmax) THEN
        pidx = np_lut-1 ; pfrac = 0.0
    ENDIF
    pfrac = 1.0- (press(i)-p_lut(pidx))/(p_lut(pidx+1) - p_lut(pidx))

    tidx = minval(maxloc(t_lut(1:nt_lut), MASK=(t_lut(1:nt_lut) < temp(i))))
    IF (temp(i) <= tmin) THEN
       tidx = 1 ; tfrac = 1.0
    ENDIF
    IF (temp(i) >= tmax) THEN
       tidx = nt_lut-1 ; tfrac = 0.0
    ENDIF
    tfrac = 1.0- (temp(i)-t_lut(tidx))/(t_lut(tidx+1) - t_lut(tidx))
    IF (gas_idx == o2_idx) THEN 
    crsz(1:nw_lut, i) = o2_lut(1:nw_lut, tidx, pidx)*tfrac*pfrac + &
                        o2_lut(1:nw_lut, tidx+1,pidx+1)*(1.0-tfrac)*(1.0-pfrac) + &
                        o2_lut(1:nw_lut, tidx, pidx+1)*tfrac*(1.0-pfrac) + &
                        o2_lut(1:nw_lut, tidx+1, pidx)*(1.0-tfrac)*pfrac
    ELSE IF (gas_idx == h2o_idx) THEN 
    crsz(1:nw_lut, i) = h2o_lut(1:nw_lut, tidx, pidx)*tfrac*pfrac + &
                        h2o_lut(1:nw_lut, tidx+1,pidx+1)*(1.0-tfrac)*(1.0-pfrac) + &
                        h2o_lut(1:nw_lut, tidx, pidx+1)*tfrac*(1.0-pfrac) + &
                        h2o_lut(1:nw_lut, tidx+1, pidx)*(1.0-tfrac)*pfrac
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

