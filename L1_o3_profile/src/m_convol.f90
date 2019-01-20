!
module m_convol
  
  use m_gauss, only: asym_gauss_multi, asym_gauss_vary, gauss_multi, gauss_vary, &
                     asym_gauss_f2c, asym_gauss_vary_f2c, gauss_f2c,gauss_vary_f2c
  use m_voigt, only: asym_voigt_multi, asym_voigt_vary, &
                     asym_voigt_f2c, asym_voigt_vary_f2c
  use m_triangle, only: triangle_multi, triangle_vary, &
                        triangle_f2c, triangle_vary_f2c
  USE m_super_gauss, ONLY: super_gauss_multi, super_gauss_vary, & 
                           super_gauss_f2c, super_gauss_vary_f2c
  USE OMSAO_slitfunction_module
  USE m_ezspline_interpolation, ONLY: interpolation
  USE OMSAO_precision_module     
  USE OMSAO_indices_module, ONLY: solar_idx, wvl_idx, spc_idx
  USE OMSAO_errstat_module
  USE OMSAO_variables_module,   ONLY: yn_varyslit, which_slit, & 
      n_refspec_pts,refspec_orig_data
  public convol, convol_i0effect, convol_f2c,convolf2c_i0effect
  private get_refspec !correct_coaddeffect, normalize_solar_refspec
  ! convol_f2c is much faster than convol
contains

  SUBROUTINE convol (refwav, refspec, nref,errstat)

    IMPLICIT NONE

    INTEGER,                          INTENT (IN)    :: nref
    INTEGER,                          INTENT (OUT)   :: errstat
    REAL (KIND=dp), DIMENSION (nref), INTENT (IN)    :: refwav
    REAL (KIND=dp), DIMENSION (nref), INTENT (INOUT) :: refspec

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                          :: fidx, lidx, npts
    REAL (KIND=dp), DIMENSION (nref) :: abspecmod, abspec
    !CHARACTER (LEN=*), PARAMETER    :: modulename = 'CORRECT_I0EFFECT'
    errstat = pge_errstat_ok

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
      ELSE
        CALL omislit_multi    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      END IF
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
      ELSE
        CALL omislit_vary    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
      ENDIF
    ENDIF

    refspec(fidx:lidx) = abspecmod(fidx:lidx)

    RETURN

  END SUBROUTINE convol

  SUBROUTINE convol_i0effect(refwav, refspec, nref, scalex, errstat)

    IMPLICIT NONE

    INTEGER,                       INTENT (IN)    :: nref
    REAL (KIND=dp),                INTENT (IN)    :: scalex
    REAL (KIND=dp), DIMENSION (:), INTENT (IN)    :: refwav  ! (nref)
    REAL (KIND=dp), DIMENSION (:), INTENT (INOUT) :: refspec ! (nref)
    INTEGER,                       INTENT (OUT)   :: errstat
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                          :: fidx, lidx, npts
    REAL (KIND=dp), DIMENSION (nref) :: specmod, newi0, abspecmod, abspec
    !CHARACTER (LEN=*), PARAMETER    :: modulename = 'CORRECT_I0EFFECT'
    LOGICAL, PARAMETER :: weight_irrad = .true.
    errstat = pge_errstat_ok

   fidx = 1
   lidx = nref
   npts = lidx - fidx + 1

    ! Interpolate i0 to refwav positions
    CALL get_refspec( npts, refwav(fidx:lidx), newi0(fidx:lidx))

    IF (weight_irrad) THEN   
       abspec(fidx:lidx)  = newi0(fidx:lidx) * refspec(fidx:lidx)
    ELSE
       abspec(fidx:lidx) = newi0(fidx:lidx)*EXP(-refspec(fidx:lidx)*scalex) 
    ENDIF

    IF (.NOT. yn_varyslit) THEN
      IF (which_slit == 0) THEN
        CALL gauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL gauss_multi (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL asym_gauss_multi (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_multi (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL asym_voigt_multi (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 3) THEN
        CALL triangle_multi   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL triangle_multi   (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_multi   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL super_gauss_multi   (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE
        CALL omislit_multi    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL omislit_multi    (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      END IF
    ELSE 
      IF (which_slit == 0) THEN
        CALL gauss_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL gauss_vary (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL asym_gauss_vary (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL asym_voigt_vary (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 3) THEN
        CALL triangle_vary   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL triangle_vary   (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE IF (which_slit == 4) THEN 
        CALL super_gauss_vary   (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL super_gauss_vary   (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ELSE
        CALL omislit_vary    (refwav(fidx:lidx), abspec(fidx:lidx), abspecmod(fidx:lidx), npts)
        CALL omislit_vary    (refwav(fidx:lidx), newi0(fidx:lidx),  specmod(fidx:lidx),   npts)
      ENDIF
    ENDIF

    IF (weight_irrad) THEN
       refspec(fidx:lidx) = abspecmod(fidx:lidx) / specmod(fidx:lidx)
    ELSE
       refspec(fidx:lidx) = - LOG(abspecmod(fidx:lidx) / specmod(fidx:lidx)) /  scalex
    ENDIF

    RETURN

  END SUBROUTINE convol_i0effect 

  SUBROUTINE convolf2c_i0effect (fwave, fspec0, nf, nspec,scalex,cwave, cspec, nc)
    IMPLICIT NONE
    ! ===============
    ! Input variables
    ! ===============
    INTEGER,                        INTENT (IN)         :: nc, nf, nspec
    REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
    REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec0
    REAL (KIND=dp), DIMENSION (nspec), INTENT (IN)      :: scalex
    REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
    REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=dp), DIMENSION (nf, nspec) :: fspec
    INTEGER                          :: i
    REAL (KIND=dp), DIMENSION (nf) :: i0
    REAL (KIND=dp), DIMENSION (nc) :: newi0
    !CHARACTER (LEN=*), PARAMETER    :: modulename = 'CORRECT_I0EFFECT'
    LOGICAL, PARAMETER :: weight_irrad = .true.
  

    ! Interpolate i0 to refwav positions
    CALL get_refspec( nf, fwave,i0)

    DO i = 1, nspec
      IF (weight_irrad) THEN   
         fspec(1:nf,i)  = i0 * fspec0(1:nf, i)
      ELSE
         fspec(1:nf, i)  = i0*EXP(-fspec0(1:nf, i)*scalex(i)) 
      ENDIF
    ENDDO

    IF (.NOT. yn_varyslit) THEN
      IF (which_slit == 0) THEN
        CALL gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL gauss_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL asym_gauss_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL asym_voigt_f2c(fwave, i0, nf, 1, cwave, newi0, nc) 
      ELSE IF (which_slit == 3) THEN
        CALL triangle_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL triangle_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 4) THEN
        CALL super_gauss_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL super_gauss_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 5) THEN
        CALL omislit_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL omislit_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      END IF
    ELSE
      IF (which_slit == 0) THEN
        CALL gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL gauss_vary_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 1) THEN
        CALL asym_gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL asym_gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ELSE IF (which_slit == 2) THEN
        CALL asym_voigt_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL asym_voigt_vary_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 3) THEN
        CALL triangle_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL triangle_vary_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 4) THEN
        CALL super_gauss_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL super_gauss_vary_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ELSE IF (which_slit == 5) THEN
        CALL omislit_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
        CALL omislit_vary_f2c(fwave, i0, nf, 1, cwave, newi0, nc)
      ENDIF
    ENDIF
 
    DO i = 1, nspec
      IF (weight_irrad) THEN   
         cspec(1:nc, i)  = cspec(1:nc, i)/newi0 
      ELSE
         cspec(1:nc, i)  = -LOG(cspec(1:nc, i)/newi0) /scalex(i)
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE convolf2c_i0effect

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
        CALL omislit_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      END IF
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
        CALL omislit_vary_f2c(fwave, fspec, nf, nspec, cwave, cspec, nc)
      ENDIF
    ENDIF

    RETURN
  END SUBROUTINE convol_f2c

  SUBROUTINE get_refspec (nwvl, wvl, refsol)
    IMPLICIT NONE
    ! INPUT
    INTEGER, INTENT(IN) :: nwvl
    REAL (kind=dp), DIMENSION(:), INTENT(IN) :: wvl  ! nwvl
    ! outpit
    REAL (kind=dp), DIMENSION(nwvl), INTENT(OUT) :: refsol
    ! local
    INTEGER :: errstat
    CALL interpolation (n_refspec_pts(solar_idx), &
         refspec_orig_data (solar_idx, 1:n_refspec_pts(solar_idx), wvl_idx), &
         refspec_orig_data (solar_idx, 1:n_refspec_pts(solar_idx), spc_idx), &
         nwvl, wvl, refsol, errstat)
    IF (errstat /=0 ) THEN 
       WRITE(*,*) 'get_refspec: interpolation error'
       STOP 1
    ENDIF     
    RETURN  
  END SUBROUTINE get_refspec

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
