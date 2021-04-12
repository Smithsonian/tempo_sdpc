!> Subroutines to apply slit functions to high-resoution spectra (convolution)
!> author J.Bak 2018.12
!> which_slit = 0 (sym gauss) 1 (asym gauss) 2 (voigt) 3 (triangle) 
!>              4 (sym super) 5 (asym super) 6 (instrment slit function)
!> convol (hwave, hspec, nwave) : convolution is applied without interpolation
!> convol_f2c (hwave, hspec, nhwave, nspec, cwave, cspec, ncwave) : both convolution and interpolation is applied

MODULE m_convol
  
  USE OMSAO_precision_module     
  USE OMSAO_indices_module,     ONLY: solar_idx, wvl_idx, spc_idx, &
      spk_idx, hwr_idx, hwl_idx, vgl_idx, vgr_idx, asy_idx, hwe_idx, &
      instrument_idx, omi_idx, gome_idx, gome2_idx, asy_idx, hwe_idx, &
      spk_idx, tempo_idx
  USE OMSAO_variables_module,   ONLY: yn_varyslit, which_slit, &
       fixslitcal,instrument_sidx, n_refspec_pts, refspec_orig_data, &
       fitvar_sol, solwinfit, mean_hw1e, mean_asym, mean_shape, numwin

  USE OMSAO_errstat_module
  USE m_ezspline_interpolation, ONLY: interpolation

  use m_gauss, only: asym_gauss, gauss, gauss_tophat_f2c, &
                     asym_gauss_multi, asym_gauss_vary, gauss_multi, gauss_vary, &
                     asym_gauss_f2c, asym_gauss_vary_f2c, gauss_f2c,gauss_vary_f2c
  use m_voigt, only: asym_voigt, asym_voigt_multi, asym_voigt_vary, &
                     asym_voigt_f2c, asym_voigt_vary_f2c
  use m_triangle, only: triangle, triangle_multi, triangle_vary, &
                        triangle_f2c, triangle_vary_f2c
  USE m_super_gauss, ONLY: super_gauss, super_agauss, super_gauss_multi, super_gauss_vary, & 
                           super_gauss_f2c, super_gauss_vary_f2c, & 
                           super_agauss_multi, super_agauss_vary, & 
                           super_agauss_f2c, super_agauss_vary_f2c
  USE OMSAO_slitfunction_module
  use tell_module

  PUBLIC simple_convol, convol, convol_i0, convol_f2c,convol_i0f2c, convol_f2c_stk, get_i0
  PRIVATE !get_i0 !correct_coaddeffect, normalize_solar_refspec convol_f2c is much faster than convol

  CONTAINS 

  SUBROUTINE  simple_convol (kppos, kpspec, kpspec_gauss, npts) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npts
  REAL (KIND=dp), DIMENSION (npts), INTENT(IN) :: kppos, kpspec
  REAL (KIND=dp), DIMENSION (npts), INTENT(OUT) :: kpspec_gauss
  

    IF (.NOT. yn_varyslit .OR. fixslitcal) THEN 
      IF (which_slit == 0) THEN 
         CALL gauss (kppos, kpspec, kpspec_gauss, npts, &
              fitvar_sol(hwe_idx))
      ELSE IF (which_slit == 1) THEN 
         CALL asym_gauss (kppos, kpspec, kpspec_gauss, npts, &
              fitvar_sol(hwe_idx), fitvar_sol(asy_idx))
      ELSE IF (which_slit == 2) THEN 
         CALL asym_voigt (kppos, kpspec, kpspec_gauss, npts, &
              fitvar_sol(vgl_idx), fitvar_sol(vgr_idx), &
              fitvar_sol(hwl_idx), fitvar_sol(hwr_idx) )
      ELSE IF (which_slit == 3) THEN 
         CALL triangle (kppos, kpspec, kpspec_gauss, npts, &
              fitvar_sol(hwe_idx))
      ELSE IF (which_slit == 4) THEN 
         CALL super_gauss (kppos, kpspec, kpspec_gauss, npts, &
              fitvar_sol(hwe_idx), fitvar_sol(spk_idx))
      ELSE IF (which_slit == 5) THEN 
         CALL super_agauss (kppos, kpspec, kpspec_gauss, npts, &
              fitvar_sol(hwe_idx), fitvar_sol(asy_idx), fitvar_sol(spk_idx))
      ELSE IF (which_slit >= instrument_sidx) THEN 
        IF (instrument_idx == omi_idx) THEN 
          CALL omislit_multi (kppos, kpspec, kpspec_gauss, npts)
        else if (instrument_idx == tempo_idx) then
          solwinfit(1:numwin,asy_idx,1) = mean_asym
          solwinfit(1:numwin,hwe_idx,1) = mean_hw1e
          solwinfit(1:numwin,spk_idx,1) = mean_shape
          CALL super_agauss_multi (kppos, kpspec, kpspec_gauss, npts)
        ENDIF
      ENDIF
    ELSE
      IF (which_slit == 0) THEN 
         CALL gauss_vary (kppos, kpspec, kpspec_gauss, npts)
      ELSE IF (which_slit == 1) THEN 
         CALL asym_gauss_vary (kppos, kpspec, kpspec_gauss, npts)
      ELSE IF (which_slit == 2) THEN 
         CALL asym_voigt_vary (kppos, kpspec, kpspec_gauss, npts)
      ELSE IF (which_slit == 3) THEN 
         CALL triangle_vary (kppos, kpspec, kpspec_gauss, npts)
      ELSE IF (which_slit == 4) THEN 
         CALL super_gauss_vary  (kppos, kpspec, kpspec_gauss, npts)
      ELSE IF (which_slit == 5) THEN 
         CALL super_agauss_vary  (kppos, kpspec, kpspec_gauss, npts)
      ELSE IF (which_slit >= instrument_sidx) THEN 
        IF (instrument_idx == omi_idx) THEN 
         CALL omislit_vary  (kppos, kpspec, kpspec_gauss, npts)
        else if (instrument_idx == tempo_idx) then
          print *, "Simple convolve should not be used with TEMPO"
          stop 1
        ENDIF
      ENDIF
    ENDIF  

  END SUBROUTINE simple_convol

  SUBROUTINE convol (refwav, refspec, nref)

    IMPLICIT NONE

    INTEGER,                          INTENT (IN)    :: nref
    REAL (KIND=dp), DIMENSION (nref), INTENT (IN)    :: refwav
    REAL (KIND=dp), DIMENSION (nref), INTENT (INOUT) :: refspec

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                          :: fidx, lidx, npts
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: abspecmod, abspec

    allocate (abspecmod(nref), abspec(nref))
    fidx = 1
    lidx = nref
    npts = lidx - fidx + 1  
    abspec(fidx:lidx)  = refspec(fidx:lidx)
    
    IF (.NOT. yn_varyslit) THEN
      IF (which_slit == 0) THEN
        CALL gauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 3) THEN
        CALL triangle_multi   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_multi   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 5) THEN 
        CALL super_agauss_multi   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == instrument_sidx) THEN
        IF (instrument_idx == omi_idx) THEN  
          CALL omislit_multi    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        else if (instrument_idx == tempo_idx) then
          solwinfit(1:numwin,asy_idx,1) = mean_asym
          solwinfit(1:numwin,hwe_idx,1) = mean_hw1e
          solwinfit(1:numwin,spk_idx,1) = mean_shape
          CALL super_agauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), &
               abspecmod(fidx:lidx), npts)
        ENDIF
      ENDIF
    ELSE 
      IF (which_slit == 0) THEN
        CALL gauss_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 3) THEN
        CALL triangle_vary   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_vary   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == 5) THEN 
        CALL super_agauss_vary   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ELSE IF (which_slit == instrument_sidx) THEN 
        IF (instrument_idx == omi_idx) THEN  
          CALL omislit_vary    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        else if (instrument_idx == tempo_idx) then
          print *, "Variable instrument slit convolution not implemented"
          stop 1
        ENDIF
      ENDIF
    ENDIF

    refspec(fidx:lidx) = abspecmod(fidx:lidx)
    deallocate (abspec, abspecmod)
    RETURN

  END SUBROUTINE convol

  SUBROUTINE convol_i0(refwav, refspec, npts, scalex)

    IMPLICIT NONE

    INTEGER,                          INTENT (IN)    :: npts
    REAL (KIND=dp),                   INTENT (IN)    :: scalex
    REAL (KIND=dp), DIMENSION (npts), INTENT (IN)    :: refwav
    REAL (KIND=dp), DIMENSION (npts), INTENT (INOUT) :: refspec
    ! ---------------
    ! Local variables
    ! ---------------
    LOGICAL, PARAMETER               :: weight_irrad = .true.
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: i0, abspec

    allocate (i0(npts), abspec(npts))

    ! Interpolate i0 to refwav positions
    CALL get_i0( npts, refwav, i0)

    IF (weight_irrad) THEN   
       abspec  = i0 * refspec
    ELSE
       abspec  = i0*EXP(-refspec*scalex) 
    ENDIF

    CALL convol (refwav, i0, npts)
    CALL convol (refwav, abspec, npts)

    IF (weight_irrad) THEN
       refspec = abspec / i0
    ELSE
       refspec = - LOG(abspec / i0) /  scalex
    ENDIF
    deallocate (i0, abspec)
    RETURN

  END SUBROUTINE convol_i0

  SUBROUTINE convol_i0f2c (fwave, fspec0, nf, nspec,scalex,cwave, cspec, nc)
    IMPLICIT NONE
    ! ===============
    ! Input variables
    ! ===============
    INTEGER,                        INTENT (IN)         :: nc, nf, nspec
    REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
    REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec0
    REAL (KIND=dp), INTENT (IN)                         :: scalex
    REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
    REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   :: i
    REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: fspec !(nf, nspec)
    REAL (KIND=dp), DIMENSION (:)  , ALLOCATABLE :: i0 !(nf)
    REAL (KIND=dp), DIMENSION (:)  , ALLOCATABLE :: newi0 !(nc)
    !CHARACTER (LEN=16), PARAMETER    :: modulename = 'CORRECT_I0EFFECT'
    LOGICAL, PARAMETER :: weight_irrad = .true.
  
    allocate (fspec(nf, nspec), i0(nf), newi0(nc))

    ! Interpolate i0 to refwav positions
    CALL get_i0( nf, fwave,i0)

    DO i = 1, nspec
      IF (weight_irrad) THEN   
         fspec(1:nf,i)  = i0 * fspec0(1:nf, i)
      ELSE
         fspec(1:nf,i)  = i0*EXP(-fspec0(1:nf, i)*scalex) 
      ENDIF
    ENDDO

    CALL convol_f2c (fwave, i0, nf, 1, cwave, newi0, nc)
    CALL convol_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)

    DO i = 1, nspec
      IF (weight_irrad) THEN   
         cspec(1:nc, i)  = cspec(1:nc, i)/newi0 
      ELSE
         cspec(1:nc, i)  = -LOG(cspec(1:nc, i)/newi0) /scalex
      ENDIF
    ENDDO

    deallocate (fspec, i0, newi0)

    RETURN
  END SUBROUTINE convol_i0f2c

  SUBROUTINE convol_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
    IMPLICIT NONE
    ! ===============
    ! Input variables
    ! ===============
    INTEGER,                        INTENT (IN)         :: nc, nf, nspec
    REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
    REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
    REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
    REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

    IF (.NOT. yn_varyslit) THEN
      IF (which_slit == 0) THEN
        CALL gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 3) THEN
        CALL triangle_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 4) THEN
        CALL super_gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 5) THEN
        CALL super_agauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == instrument_sidx) THEN
        IF (instrument_idx == omi_idx) THEN  
             CALL omislit_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
          solwinfit(1:numwin,asy_idx,1) = mean_asym
          solwinfit(1:numwin,hwe_idx,1) = mean_hw1e
          solwinfit(1:numwin,spk_idx,1) = mean_shape
          CALL super_agauss_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
        ENDIF
      ENDIF
    ELSE
      IF (which_slit == 0) THEN
        CALL gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 3) THEN
        CALL triangle_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 4) THEN
        CALL super_gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 5) THEN
        CALL super_agauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == instrument_sidx) THEN
        IF (instrument_idx == omi_idx) THEN  
          CALL omislit_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        else if (instrument_idx == tempo_idx) then
          print *, "Variable instrument slit function not implemented"
          stop 1
        ENDIF
      ENDIF
    ENDIF

    RETURN
  END SUBROUTINE convol_f2c

  SUBROUTINE convol_f2c_stk(fwave, fspec, nf, cwave, cspec, nc)

  USE OMSAO_precision_module
  
  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)  :: nc, nf
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)  :: fwave, fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)  :: cwave
  REAL (KIND=dp), DIMENSION (nc), INTENT (OUT) :: cspec
 
  CALL gauss_tophat_f2c (fwave, fspec, nf, cwave, cspec, nc)

  END SUBROUTINE convol_f2c_stk

  SUBROUTINE get_i0 (ni0, wvl, i0)
    !USE ozprof_data_module, ONLY: hw1e=>hres_slitwidth
    IMPLICIT NONE
    ! INPUT/OUTPUT variables
    INTEGER, INTENT(IN) :: ni0
    REAL (kind=dp), DIMENSION(ni0), INTENT(IN) :: wvl
    REAL (kind=dp), DIMENSION(ni0), INTENT(OUT):: i0
    ! local variables
    INTEGER :: errstat, npts
    REAL (KIND=dp), ALLOCATABLE, DIMENSION (:) :: wave, spec
    CHARACTER (len=20), PARAMETER :: modulename='get_i0'

    npts = n_refspec_pts(solar_idx)
    allocate (wave(npts), spec(npts))
    wave = refspec_orig_data(solar_idx, 1:npts, wvl_idx)
    spec = refspec_orig_data(solar_idx, 1:npts, spc_idx)
    !IF (hw1e .ge. 0.03D0 ) THEN
    !   CALL gauss(wave, spec, spec, npts, 0.05D0)
    !ENDIF
    CALL interpolation (npts, wave, spec, ni0, wvl, i0, errstat)

    IF (errstat /=0 ) THEN
       WRITE(*,*) ADJUSTL(TRIM(modulename))//': interpolation error', errstat
       WRITE(*,*) 'out', wvl(1), wvl(ni0), ni0
       WRITE(*,*) 'in', refspec_orig_data(solar_idx, 1, wvl_idx)
       WRITE(*,*) 'in', refspec_orig_data(solar_idx, n_refspec_pts(solar_idx), wvl_idx)
       STOP 1
    ENDIF
    deallocate (wave, spec)
    RETURN
  END SUBROUTINE get_i0

!  SUBROUTINE CORRECT_COADDEFFECT(refwav, refspec, i0, nref, choice, nout, errstat)
!    USE OMSAO_precision_module     
!    USE OMSAO_errstat_module
!    IMPLICIT NONE
!
!    INTEGER,                          INTENT (IN)    :: nref, choice
!    INTEGER,                          INTENT (OUT)   :: errstat, nout
!    REAL (KIND=dp), DIMENSION (nref), INTENT (IN)    :: refwav, i0
!    REAL (KIND=dp), DIMENSION (nref), INTENT (INOUT) :: refspec
!
!    ! ---------------
!    ! Local variables
!    ! ---------------
!    !INTEGER                          :: i, j
!    REAL (KIND=dp), DIMENSION (nref) :: newi0, abspec
!    REAL (KIND=dp)                   :: scalex!, frefw, lrefw
!    !CHARACTER (LEN=19), PARAMETER    :: modulename = 'CORRECT_COADDEFFECT'
!
!    errstat = pge_errstat_ok
!    scalex = 0.08              !  Arbitarily assummed
!    newi0 = i0
!
!    IF (choice == 1 .OR. choice == 2) THEN 
!      abspec = i0 * EXP(-refspec * scalex)
!    ENDIF
!
!    IF (choice == 1) THEN 
!      CALL avg_band_refspec(refwav, abspec, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!      CALL avg_band_refspec(refwav, newi0, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!    ELSE IF (choice == 2) THEN     
!      CALL avg_band_ozcrs(refwav, abspec, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!      CALL avg_band_ozcrs(refwav, newi0, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!    ELSE IF (choice == 3) THEN
!      CALL avg_band_refspec(refwav, refspec, nref, nout, errstat)
!      IF ( errstat >= pge_errstat_error ) RETURN
!    ENDIF
!
!
!    IF (choice == 1 .OR. choice == 2) THEN
!      refspec(1:nout) = -LOG(abspec(1:nout) / newi0(1:nout)) / scalex
!    ENDIF
!
!    RETURN
!
!  END SUBROUTINE CORRECT_COADDEFFECT


!   * here super_gauss is not implemented (jbak)
!  SUBROUTINE normalize_solar_refspec ( n_radwvl, curr_rad_wvl, solar_spec, errstat)
!
!    USE OMSAO_precision_module
!    USE OMSAO_indices_module,     ONLY: solar_idx, wvl_idx, spc_idx
!    USE OMSAO_variables_module,   ONLY: n_refspec_pts, refspec_orig_data, &
!         yn_varyslit, which_slit, refspec_norm, solar_refspec
!    USE OMSAO_parameters_module,  ONLY: max_spec_pts
!    USE ozprof_data_module,       ONLY: div_sun
!    USE OMSAO_slitfunction_module    
!    USE OMSAO_errstat_module
!    IMPLICIT NONE
!
!    INTEGER,                              INTENT (IN)  :: n_radwvl
!    REAL (KIND=dp), DIMENSION (n_radwvl), INTENT (IN)  :: curr_rad_wvl, solar_spec
!    INTEGER,                              INTENT (OUT) :: errstat
!
!    ! ---------------
!    ! Local variables
!    ! ---------------
!    INTEGER                                  :: fidx, lidx, ni0!, i
!    REAL (KIND=dp), DIMENSION (max_spec_pts) :: wave, specmod, ratio
!    REAL (KIND=dp), DIMENSION (n_radwvl)     :: solar_spec0, ratio0
!    REAL (KIND=dp)                           :: frefw, lrefw
!
!    ! ------------------
!    ! External functions
!    ! ------------------
!    !INTEGER :: OMI_SMF_setmsg
!
!    CHARACTER (LEN=23), PARAMETER :: modulename = 'normalize_solar_refspec'
!
!    errstat = pge_errstat_ok
!
!    ! Convole high-resolution solar reference spectrum
!    ni0  = n_refspec_pts(solar_idx)
!    wave = refspec_orig_data(solar_idx,1:ni0,wvl_idx)
!
!    IF (.NOT. yn_varyslit) THEN
!      IF (which_slit == 0) THEN
!        CALL gauss_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 1) THEN
!        CALL asym_gauss_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 2) THEN
!        CALL asym_voigt_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 3) THEN
!        CALL triangle_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE
!        CALL omislit_multi (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      END IF
!    ELSE 
!      IF (which_slit == 0) THEN
!        CALL gauss_vary (wave(1:ni0), solar_refspec(1:ni0),specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 1) THEN
!        CALL asym_gauss_vary (wave(1:ni0), solar_refspec(1:ni0),specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 2) THEN
!        CALL asym_voigt_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE IF (which_slit == 3) THEN
!        CALL triangle_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ELSE
!        CALL omislit_vary (wave(1:ni0), solar_refspec(1:ni0), specmod(1:ni0), ni0)
!      ENDIF
!    ENDIF
!
!    CALL bspline(wave(1:ni0), specmod(1:ni0), ni0, curr_rad_wvl, solar_spec0, n_radwvl, errstat)
!    ratio0 = solar_spec / solar_spec0 * div_sun / refspec_norm(solar_idx)
!
!    frefw = curr_rad_wvl(1)
!    lrefw = curr_rad_wvl(n_radwvl)
!    fidx = MINVAL(MINLOC(wave(1:ni0), MASK=(wave(1:ni0) >= frefw)))
!    lidx = MINVAL(MAXLOC(wave(1:ni0), MASK=(wave(1:ni0) <= lrefw)))
!
!    CALL bspline(curr_rad_wvl, ratio0, n_radwvl, wave(fidx:lidx), ratio(fidx:lidx), lidx-fidx + 1, errstat)
!    IF (errstat < 0) THEN
!      WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
!      errstat = pge_errstat_error
!    ENDIF
!    IF (fidx > 1) ratio(1:fidx-1) = ratio(fidx)
!    IF (fidx < ni0) ratio(lidx+1:ni0) = ratio(lidx)
!
!    refspec_orig_data(solar_idx, 1:ni0, spc_idx) = &
!         refspec_orig_data(solar_idx, 1:ni0, spc_idx) * ratio(1:ni0)
!
!    RETURN
!  END SUBROUTINE normalize_solar_refspec


end module m_convol
