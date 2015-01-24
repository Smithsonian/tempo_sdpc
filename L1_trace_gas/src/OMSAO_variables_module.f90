MODULE OMSAO_variables_module

  USE OMSAO_precision_module, ONLY: i4, r8, r4, i2
  USE OMSAO_indices_module, ONLY: &
    max_rs_idx, max_calfit_idx, n_max_fitpars, mxs_idx, sig_idx, icf_idx, &
    o3_t1_idx, o3_t2_idx, o3_t3_idx,                                      &
    us1_idx, us2_idx, n_voc_amf_luns, ccd_idx, radfit_idx

  USE OMSAO_parameters_module,   ONLY: MAX_STR_LEN, max_spec_pts, N_FIT_WINWAV, &
    max_mol_fit, nwavel_max, nxtrack_max

  IMPLICIT NONE

  PRIVATE i4, r8, r4, i2, &
    max_rs_idx, max_calfit_idx, n_max_fitpars, mxs_idx, sig_idx, icf_idx, &
    o3_t1_idx, o3_t2_idx, o3_t3_idx, us1_idx, us2_idx, n_voc_amf_luns, &
    ccd_idx, radfit_idx, &
    MAX_STR_LEN, max_spec_pts, N_FIT_WINWAV, max_mol_fit, nwavel_max, nxtrack_max

  ! -----------------------------
  ! Variables read from PCF file
  ! -----------------------------

  ! ---------------------
  ! * Verbosity threshold
  ! ---------------------
  INTEGER (KIND=I4) :: verb_thresh_lev
  ! ----------------------------------------------------
  ! * Orbit number, Version ID, L1B radiance OPF version
  ! ----------------------------------------------------
  INTEGER (KIND=I4) :: orbit_number, ecs_version_id, l1br_opf_version

  ! -------------------------------------------------
  ! Variables defined in preamble of original program
  ! -------------------------------------------------

  INTEGER (KIND=I4), DIMENSION (n_max_fitpars)  :: mask_fitvar_rad, mask_fitvar_cal, all_radfit_idx

  INTEGER (KIND=I4)                             :: n_fitvar_rad, n_fitvar_cal

  REAL    (KIND=r8), DIMENSION (max_calfit_idx) :: &
    fitvar_cal, fitvar_cal_saved, fitvar_sol_init, lo_sunbnd, up_sunbnd

  REAL    (KIND=r8), DIMENSION (n_max_fitpars)  :: fitvar_rad, fitvar_rad_init
  REAL    (KIND=r8), DIMENSION (n_max_fitpars)  :: fitvar_rad_saved
  REAL    (KIND=r8), DIMENSION (n_max_fitpars)  :: lo_radbnd, up_radbnd
  CHARACTER (LEN=6), DIMENSION (n_max_fitpars)  :: fitvar_rad_str, fitvar_sol_str

  REAL    (KIND=r8), DIMENSION (nwavel_max, max_rs_idx) :: database

  ! -------------------------------------
  ! Variables related to Air Mass Factors
  ! -------------------------------------
  !INTEGER (KIND=I4)  :: n_amftab_dim, n_amftab_ang
  LOGICAL            :: Have_AMF_Table
  !!REAL    (KIND=r8)  :: amf_esza_min, amf_esza_max
  !REAL    (KIND=r8)  :: amf_tab_wvl, amf_tab_alb
  !REAL    (KIND=r8), DIMENSION (n_amftab_ang_max, n_amftab_dim_max) :: amf_table_bro

  !! ---------------------------------------------------------------
  !! Some BrO specific AMF variables: Fitted polynomial coefficients
  !! that are used in lieu of interpolation.
  !! ---------------------------------------------------------------
  !INTEGER (KIND=i4)                                  :: n_amftab_coef
  !REAL    (KIND=r8)                                  :: amftab_coef_norm
  !REAL    (KIND=r8), DIMENSION (0:n_amftab_coef_max) :: amftab_coef

  ! -----------------------------
  ! Previously IMPLICIT variables
  ! -----------------------------
  REAL (KIND=r8) :: Undersample_Phase, szamax, chisq, sol_wav_avg
  REAL (KIND=r8) :: Slit_Half_Width_1e, Slit_Asym_Factor
  REAL (KIND=r4) :: zatmos

  ! -----------------------------------------------------------
  ! Variables related to reference spectra
  !  * number of reference spectra:       N_REfSPEC
  !  * indentification strings:           FITPAR_IDXNAME
  !  * file names with reference spectra: REFSPEC_FNAME
  !  * original (uniterpolated) data:     REFSPEC_ORIG_DATA
  !  * number of spectral points:         N_REFSPEC_PTS
  !  * first and last wavelenghts:        REFSPEC_FIRSTLAST_WAV
  ! -----------------------------------------------------------
  INTEGER   (KIND=I4) :: n_refspec

  ! -------------------------------------
  ! TYPE declaration for Reference Specta
  ! -------------------------------------
  TYPE, PUBLIC :: reference_spectrum_type
    CHARACTER (LEN=MAX_STR_LEN)                      :: Title, Units
    CHARACTER (LEN=MAX_STR_LEN)                      :: FileName
    CHARACTER (LEN=MAX_STR_LEN)                      :: FittingIdxName
    INTEGER   (KIND=I4)                           :: nPoints
    REAL      (KIND=r8)                           :: NormFactor, Temperature
    REAL      (KIND=r8), DIMENSION (2)            :: FirstLastWav
    REAL      (KIND=r8), DIMENSION (max_spec_pts) :: RefSpecWavs
    REAL      (KIND=r8), DIMENSION (max_spec_pts) :: RefSpecData
  END TYPE reference_spectrum_type

  ! -----------------------------------------
  ! TYPE declaration for Common Mode Spectrum
  ! -----------------------------------------
  TYPE, PUBLIC :: common_mode_spectrum_type
    CHARACTER (LEN=MAX_STR_LEN)                                  :: Title, Units
    CHARACTER (LEN=MAX_STR_LEN)                                  :: FileName
    CHARACTER (LEN=MAX_STR_LEN)                                  :: FittingIdxName
    INTEGER   (KIND=I4)                                       :: nPoints
    REAL      (KIND=r8)                                       :: NormFactor, Temperature
    REAL      (KIND=r8), DIMENSION (2)                        :: FirstLastWav
    !INTEGER   (KIND=I2), DIMENSION (nxtrack_max,2)            :: CCDPixel
    INTEGER   (KIND=I2), DIMENSION (:,:), allocatable          :: CCDPixel
    !INTEGER   (KIND=I4), DIMENSION (nxtrack_max)              :: RefSpecCount
    INTEGER   (KIND=I4), DIMENSION (:), allocatable            :: RefSpecCount
    !REAL      (KIND=r8), DIMENSION (max_spec_pts,nxtrack_max) :: RefSpecWavs
    REAL      (KIND=r8), DIMENSION (:,:), allocatable          :: RefSpecWavs
    !REAL      (KIND=r8), DIMENSION (max_spec_pts,nxtrack_max) :: RefSpecData
    REAL      (KIND=r8), DIMENSION (:,:), allocatable          :: RefSpecData
    !integer   (kind=i4), dimension (nxtrack_max)              :: num_wavelengths = 0
    integer   (kind=i4), dimension (:), allocatable            :: num_wavelengths
  END TYPE common_mode_spectrum_type

  ! -------------------------------
  ! Array for all Reference Spectra
  ! -------------------------------
  !TYPE (reference_spectrum_type),  DIMENSION (max_rs_idx) :: refspecs_original
  TYPE (reference_spectrum_type),  DIMENSION (:), allocatable :: refspecs_original
  TYPE (common_mode_spectrum_type)                            :: common_mode_spec

  ! -------------------------------------------
  ! A special beast: The undersampling spectrum
  ! -------------------------------------------
  LOGICAL, DIMENSION (us1_idx:us2_idx) :: have_undersampling

  CHARACTER (LEN=MAX_STR_LEN), DIMENSION (icf_idx:max_rs_idx) :: static_input_fnames

  ! --------------------------------------
  ! Solar and Earth shine wavlength limits
  ! --------------------------------------
  REAL (KIND=r8), DIMENSION (N_FIT_WINWAV) :: ctrl_fit_winwav_lim
  REAL (KIND=r8), DIMENSION (2)            :: ctrl_fit_winexc_lim
  REAL (KIND=r8)                           :: winwav_min, winwav_max

  ! ------------------------------------------------------------------
  ! Indices of fitting window defining wavelengths in current spectrum
  ! ------------------------------------------------------------------
  INTEGER (KIND=I4), DIMENSION (N_FIT_WINWAV) :: fit_winwav_idx

  ! --------------------------------------------------------------------------
  ! The current solar and radiance spectrum, including wavelengths and weights
  ! --------------------------------------------------------------------------
  REAL    (KIND=r8), DIMENSION (nwavel_max, ccd_idx)    :: curr_sol_spec

  ! ----------------------------------------
  ! Pixel number limits:
  !  * 1,2: First and last scan line number
  !  * 3,4: First and last cross-track pixel
  ! ----------------------------------------
  INTEGER (KIND=I4), DIMENSION (4) :: pixnum_lim

  ! ----------------------------------------------------
  ! Latitude limits: Lower and upper latitude to process
  ! ----------------------------------------------------
  REAL (KIND=r4), DIMENSION (2) :: radfit_latrange

  ! ---------------------------------------------------------------------
  ! Contstraints on the fitting residual: Window and Number of Iterations
  ! ---------------------------------------------------------------------
  INTEGER (KIND=i4), DIMENSION (radfit_idx)  :: &
    ctrl_n_fitres_loop, ctrl_fitres_range
  REAL    (KIND=r8), DIMENSION (nxtrack_max) :: xtrack_fitres_limit

  ! --------------------------------------------
  ! Frequency of radiance wavelength calibration
  ! --------------------------------------------
  INTEGER (KIND=I4) :: radwavcal_freq

  ! -----------------------------------------------------------------
  ! Variables connected with numerical precision/convergence criteria
  ! -----------------------------------------------------------------
  REAL (KIND=r8) :: tol,  epsrel,  epsabs,  epsx

  ! ----------------------------------------------------------------------
  ! Index for the fitting parameters carrying the fitted column value.
  !
  ! N_MOL_FIT:    Number of "molecules" that carry the final column; this
  !               can be one molecule at different temperatures.
  ! FITCOL_IDX:   The main molecule indices, corresponding to the list of
  !               reference spectra.
  ! FINCOL_IDX:   For the final summation of the fitted column: The total
  !               number is the number of different molecules times the
  !               allowed sub-indices. The second dimension is for the
  !               reference spectrum index - this eases the final sum over
  !               the fitted columns (see RADIANCE_FIT subroutine).
  !                includes
  !               the subindices, hence the dimension.
  ! N_FINCOL_IDX: Number of final column indices.
  ! ----------------------------------------------------------------------
  INTEGER (KIND=I4)                                    :: n_mol_fit, n_fincol_idx
  INTEGER (KIND=I4), DIMENSION (max_mol_fit)           :: fitcol_idx
  INTEGER (KIND=I4), DIMENSION (2,max_mol_fit*mxs_idx) :: fincol_idx

  ! ------------------------------------
  ! Maximum number of fitting iterations
  ! ------------------------------------
  INTEGER (KIND=I4) :: max_itnum_sol, max_itnum_rad

  ! ------------------------------------
  ! Maximum good column amount
  ! ------------------------------------
  REAL (KIND=r8) :: max_good_col

  ! ---------------------
  ! L1B and L2 file names
  ! ---------------------
  CHARACTER (LEN=MAX_STR_LEN) :: l1b_rad_filename, l1b_irrad_filename, l2_filename

  ! -----------------------------------------------------------------
  ! Generic dimension variables (initialized from either GOME or OMI)
  ! -----------------------------------------------------------------
  INTEGER (KIND=I4) :: n_rad_wvl, n_database_wvl
  INTEGER (KIND=I4) :: n_rad_wvl_max

  ! --------------------------------------------
  ! Name of the tabulated OMI slit function data
  ! --------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: omi_slitfunc_fname

  ! ---------------------------------------------------------
  ! Filenames specific for the AMF scheme in OMBRO and OMHCHO
  ! ---------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN)                               :: OMBro_AMF_Filename
  CHARACTER (LEN=MAX_STR_LEN), DIMENSION (N_VOC_AMF_LUNS)   :: VOC_AMF_Filenames

  ! ---------------------------------------------------------------
  ! Filename, logical and type indices for composite Solar Spectrum
  ! ---------------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: OMSAO_solcomp_filename
  !INTEGER (KIND=i4)        :: solar_comp_orb

  ! -------------------------------------------------------
  ! Filename and logical for solar monthly average spectrum
  ! -------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: OMSAO_solmonthave_filename

  ! --------------------------------------------
  ! * Logical for Common Mode Iteration
  ! * Index position of Common Mode add-on
  ! * Array for initial fitting variables
  ! * Latitude range for common mode computation
  ! --------------------------------------------
  INTEGER (KIND=i4)                :: common_fitpos
  REAL    (KIND=r8), DIMENSION (3) :: common_fitvar
  REAL    (KIND=r4), DIMENSION (2) :: common_latrange
  INTEGER (KIND=i4), DIMENSION (2) :: common_latlines

  ! ---------------------------------------------------------------------
  ! 3-letter string to identify the OMI channel (UV2 or VIS) we are using
  ! ---------------------------------------------------------------------
  CHARACTER (LEN=3) :: l1b_channel

  ! --------------------------------------------------------
  ! Filename and logical for Reference Sector Correction gga
  ! --------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: OMSAO_refseccor_filename
  CHARACTER (LEN=MAX_STR_LEN) :: OMSAO_refseccor_cld_filename

  ! -----------------------------------------------------------------
  ! Logical for Scattering Weights, Gas Profile and Averaging Kernels
  ! Also filename
  ! -----------------------------------------------------------------
  CHARACTER (LEN=MAX_STR_LEN) :: OMSAO_OMLER_filename

  ! -------------------------------------------------------
  ! Variables connected with  a radiance reference spectrum
  ! -------------------------------------------------------
  INTEGER (KIND=i4)                :: target_npol
  INTEGER (KIND=i4), DIMENSION (2) :: radiance_reference_lnums
  REAL    (KIND=r4), DIMENSION (2), TARGET :: radref_latrange
  CHARACTER (LEN=MAX_STR_LEN) :: l1b_radref_filename

  ! This type gets filled in by omi_read_radiance_paras
  TYPE, PUBLIC :: Radiance_Paras_Type
    character (len=MAX_STR_LEN) :: l1bfilename
    integer (kind=i4)        :: ntimes, nxtrack, nwavel_ccd
    character (len=MAX_STR_LEN) :: swathname
    character (len=MAX_STR_LEN) :: l1bchannel
  END TYPE Radiance_Paras_Type

  ! --------------------------------
  ! Current cross-track pixel number
  ! --------------------------------
  INTEGER (KIND=i4) :: curr_xtrack_pixnum

contains

  subroutine allocate_refspec_storage (errstat)
    use tell_module
    implicit none
    integer, intent(inout) :: errstat
    if (errstat < 0) return
    allocate (refspecs_original(max_rs_idx), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "allocate_refspec_storage: allocate failed", &
                       errstat)
      return
    endif
  end subroutine allocate_refspec_storage

  subroutine allocate_common_mode_storage (cms, errstat)
    use tell_module
    implicit none
    type (common_mode_spectrum_type), intent(inout) :: cms
    integer, intent(inout) :: errstat
    if (errstat < 0) return

    allocate (cms % RefSpecWavs(nwavel_max, nxtrack_max), &
              cms % RefSpecData(nwavel_max, nxtrack_max), &
              cms % CCDPixel (nxtrack_max, 2), &
              cms % RefSpecCount (nxtrack_max), &
              cms % num_wavelengths (nxtrack_max), &
              stat = errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "allocate_common_mode_storage: allocate failed", &
                       errstat)
      return      
    endif
    cms % num_wavelengths(:) = 0
    
  end subroutine allocate_common_mode_storage

END MODULE OMSAO_variables_module
