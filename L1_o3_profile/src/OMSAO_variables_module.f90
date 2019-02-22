 ! *********************** Modification History *******************
  ! xiong liu, July 2003
  ! 1. Add three variables: the_sza_atm, the_vza_atm, the_aza_atm 
  !    Effective siewing geometry at TOA averaged on A, B, C for
  !    inputting into LIDORT
  ! 2. Add the_lat, the_lon (used in preparing atmospheric profiles)
  ! 3. Add yn_varyslit, hwlarr, hwrarr, vglarr, vgrarr, slitwav, n_slit_pts,
  !    n_slit_interval, slit_fname, slit_redo, wavcal_redo, 
  !    wavcal_fname, shiarr, squarr, sswav, for implementing 
  !    variable slit width ( for voigt profile shape only)
  ! 5. Add varaible use_meas_sig (use measurement error, otherwise
  !    use normal weight = 1 for all pixels except edge pixels)
  ! 6. Add variable for using pixel_bin, else use interpoaltion
  ! 7. Add # of wavelengths in Channel 1 and 2 respectively
  ! ****************************************************************

MODULE OMSAO_variables_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, max_calfit_idx, n_max_fitpars, mxs_idx, sig_idx, icf_idx, &
       amf_idx, n_amftab_ang_max,  n_amftab_dim_max
  USE OMSAO_parameters_module,   ONLY: &
       maxchlen, max_spec_pts, max_fit_pts, n_sol_winwav, n_rad_winwav,      &
       max_mol_fit, maxwin, maxview, maxloc, max_ref_pts, mswath, max_ring_pts,mrefl
  IMPLICIT NONE

  ! Time of a Granule
  REAL (KIND=8) :: TAI93At0ZOfGranule, TAI93StartOfGranule, GranuleSecond
  INTEGER       :: GranuleYear, GranuleMonth, GranuleDay, GranuleHour,&
                   GranuleMinute, GranuleJDay

  ! -----------------------------
  ! Variables read from PCF file
  ! -----------------------------
  ! * Current PGE name and index
  ! -----------------------------
  INTEGER           :: pge_idx
  CHARACTER (LEN=6) :: pge_name

  ! ---------------------
  ! * Verbosity threshold
  ! ---------------------
  CHARACTER (LEN=1) :: verb_thresh_char
  INTEGER           :: verb_thresh_lev
  ! -------------------------------
  ! * Swath name (really from PCF?)
  ! -------------------------------
  CHARACTER (LEN=maxchlen) :: swath_name
  ! --------------
  ! * Orbit number
  ! --------------
  INTEGER          :: pcf_orbit_number

  ! ---------------------------------------------
  ! Ground Scan Lines and xtrack pixels Limits
  ! --------------------------------------------
  INTEGER :: nswath, nxtrack, ntimes, nwavel
  INTEGER :: ntimes_loop, offset_line ! subset of fitted lines
  INTEGER, DIMENSION (2) :: linenum_lim  ! scan lines or pixel numbers
  INTEGER, DIMENSION (2) :: pixnum_lim   ! across the track pixels
  INTEGER :: the_pix, the_line           ! current pix/line num in actualy array
  INTEGER :: currloop, currline, currpix ! current pix/line num in coadded/offset array 
  CHARACTER (LEN=3) :: currpixchar
  !* coadding
  LOGICAL           :: do_xbin, do_ybin ! binning across and along the track
  INTEGER           :: nxbin, nybin, ncoadd
  !* orbit number
  INTEGER (KIND=I4) :: orbnum, orbnumsol
  !* Identifier string for irradiacne and radiance input file
  CHARACTER(LEN=6)  :: sol_identifier
  CHARACTER(LEN=6)  :: rad_identifier
  
  ! array of mirror_step indices for TEMPO synthetic data
  integer (kind=4), dimension(:), allocatable :: step_idx

  !------------------------------------------------
  ! variables used to select fitting windows
  !------------------------------------------------
  INTEGER                                :: numwin, nviswin
  REAL (KIND=dp), DIMENSION(maxwin, 2)   :: winlim
  INTEGER, DIMENSION (maxwin, 2)         :: winpix ! find locatio of winlim in the original spectrum
  LOGICAL                                :: do_bandavg
  INTEGER, DIMENSION(maxwin)             :: band_selectors, n_band_avg, n_band_samp
  !------------------------------------------------
  ! Variables related to subseting l1b spectrum 
  !------------------------------------------------
  REAL (KIND=dp)                         :: upper_spec, lower_spec
  INTEGER,        DIMENSION (mswath)     :: inschs
  REAL (KIND=dp), DIMENSION (mswath)     :: upper_wvls, lower_wvls
  REAL (KIND=dp), DIMENSION (mswath)     :: retubnd,retlbnd
  REAL (KIND=dp), DIMENSION (maxwin)     :: redslw
  REAL (KIND=dp), DIMENSION (mswath)     :: reduce_ubnd, reduce_lbnd

  !------------------------------------------------
  ! variables related to reduceing spectral resolution
  !------------------------------------------------
  LOGICAL                                :: reduce_resolution, use_redfixwav
  INTEGER                                :: reduce_slit, nredfixwav
  REAL (KIND=dp)                         :: redsampr, redlam
  REAL (KIND=dp), DIMENSION(max_fit_pts) :: redfixwav
  CHARACTER (LEN=maxchlen)               :: redfixwav_fname

  ! -------------------------------------------------
  ! control variables for spectral fitting
  ! -------------------------------------------------
  ! * Maximum number of fitting iterations
  INTEGER  :: max_itnum_sol, max_itnum_rad

  ! * which wavelength correction 
  INTEGER, PARAMETER :: correct_lambda = 2
  !1 :  lamda =( fitwavs - shi)/(1.0 + squ) ; IN OLD VERSION
  !2 :  lamda =( fitwavs - shi + sol_wav_avg*squ) /(1.0+squ) IN UPDATED

  ! * used for sol_wav_avg, rad_wav_avg
  LOGICAL :: weight_sun, weight_rad

  ! * used in determineing refwav and  radwav to avoid extrapolation 
  ! radnhtrunc: number of unused radiance pixels at each end of a spectralregion
  ! refnhextra: number of extra pixels added to the reference spectra
  INTEGER                                :: radnhtrunc, refnhextra

  ! * Whehter to filter spectral pixels around 280 and 285 nm
  LOGICAL  :: rm_mgline

  ! * Whether to do undersampling correction
  LOGICAL :: have_undersampling

  ! * whether to coadd UV2 specta to the UV-1 resolution when both UV1/UV exist
  LOGICAL  :: coadd_uv2 ! for OMI

  ! * Whether to perform wavelength calibration (due to spatial smile) before spatial coadding 
  LOGICAL  :: wcal_bef_coadd

  ! * irradiance
  LOGICAL::use_backup, use_solcomp, avg_solcomp, avgsol_allorb  

  ! * measurement error 
  LOGICAL :: use_meas_sig              ! 
  LOGICAL :: correct_merr              ! Correct OMI COL3 measurement error,xliu: 09/25/12
  LOGICAL :: xbin_decerr, ybin_decerr  ! Reduce meas error when coadding in x/ydirection, xliu: 09/25/12

  ! * slit function calibration
  INTEGER, PARAMETER        :: max_slit = 5
  INTEGER                   :: which_slit   ! 1. Gauss 2. Voigt 3. Triangle 4 super 5 instrument
  CHARACTER (3), PARAMETER, DIMENSION(0:max_slit) :: slit_name=(/'sga','aga','voi','tri','spg','ins'/)

  LOGICAL :: yn_smooth, yn_doas, yn_varyslit 

  ! * Frequency of radiance wavelength calibration
  INTEGER :: radwavcal_freq

  ! * Define SAA boundary 
  REAL (KIND=sp) :: saa_minlon = -75.0, saa_maxlon = 0.00, &
                    saa_minlat = -50.0, saa_maxlat = -5.00, &
                    saa_minlon1= 0.00 , saa_maxlon1= 30.0, &
                    saa_minlat1= -35.0, saa_maxlat1= -15.0

  LOGICAL        :: renorm ! not used anywhere
  ! * debug variable
  LOGICAL        :: debug_boreas ! not used anywhere
  ! * Full orbit processing vs. pixel processing
  LOGICAL        :: fullorb=.TRUE., do_ch2reso ! not used anywhere

  !----------------------------------------------------------------
  ! Variables related to refernce spectra
  !---------------------------------------------------------------
  ! * Directory for reference spectra and atmospheric databases
  CHARACTER (LEN=maxchlen) :: atmdbdir, refdbdir, tabdir

  ! * reference spectra at origianl resolution
  REAL (KIND=dp) :: winwav_min, winwav_max ! winwav of reference spectrum

  !CHARACTER(LEN=6),       DIMENSION (max_rs_idx):: fitpar_idxname
  !!identification strings
  INTEGER,        DIMENSION (max_rs_idx)        :: n_refspec_pts  ! number ofreferepnce spectrum
  CHARACTER (LEN=maxchlen),DIMENSION (max_rs_idx):: refspec_fname
  !REAL (KIND=dp),DIMENSION (max_rs_idx,2) :: refspec_firstlast_wav
  REAL (KIND=dp), DIMENSION (max_rs_idx)               :: refspec_norm
  REAL (KIND=dp), DIMENSION (max_rs_idx,max_spec_pts,3):: refspec_orig_data
  ! * reference spectra at instrument resolution
  INTEGER                                             :: n_refwvl, n_refwvl_sav
  INTEGER, DIMENSION (max_ref_pts)                    :: refidx, refsol_idx
  INTEGER, DIMENSION (max_ref_pts)                    :: refidx_sav
  REAL (KIND=dp), DIMENSION (max_ref_pts)             :: refwvl, refwvl_sav
  REAL (KIND=dp), DIMENSION (max_ref_pts)             :: i0sav
  REAL (KIND=dp), DIMENSION (max_ref_pts)             :: slwf, dfdsl
  REAL (KIND=dp), DIMENSION (max_rs_idx, max_ref_pts) :: database,&
                                         database_shiwf, database_save
  ! * common mode
  REAL (KIND=dp), DIMENSION (max_ref_pts) :: database_cmwf
  ! pseudo slit fitting variables
  INTEGER, PARAMETER :: max_psl=2
  INTEGER  :: npsl
  INTEGER, DIMENSION (max_psl) :: psl_fpos=0
  LOGICAL                      :: do_dsdk, do_dsdw
  REAL (KIND=dp), DIMENSION (max_psl, max_ref_pts) :: database_pslwf = 0.0
 !* Special for ozone, use Tdependent coefficients or at several T
  !INTEGER, PARAMETER :: maxozabs = 5
  !INTEGER            :: numozabs = 0, n_ozref_pts = 0
  !LOGICAL            :: oztdepend, ozconv
  !REAL (KIND=dp)            :: ozrefspec_norm
  !REAL (KIND=dp), DIMENSION (maxozabs)                  :: ozrefts
  !REAL (KIND=dp), DIMENSION (maxozabs, max_spec_pts)    :: ozrefspec
  !REAL (KIND=dp), DIMENSION (max_spec_pts) :: ozrefpos
  !REAL (KIND=dp), DIMENSION (maxozabs, max_fit_pts+4)   :: ozdb, ozdb_shiwf

  ! --------------------------------------------------------------------------
  ! The current solar and radiance spectrum, including wavelengths and weights
  ! --------------------------------------------------------------------------
  REAL (KIND=dp) :: EarthSundistance
  REAL (KIND=dp) :: curr_exposuretime ! not used anywhere
  INTEGER :: NSPC_omi
  ! * measurement spectrum in calibration process
  INTEGER ::                                n_irrad_wvl, n_rad_wvl, n_radwvl_sav
  INTEGER ::                                n_irrad_wvl_ori,n_rad_wvl_ori
  INTEGER, DIMENSION(maxwin)             :: nradpix, nsolpix,nradpix_sav
  REAL (KIND=dp), DIMENSION (sig_idx, max_fit_pts) :: curr_rad_spec,curr_rad_spec_ori
  REAL (KIND=dp), DIMENSION (sig_idx, max_fit_pts) :: curr_sol_spec,curr_sol_spec_ori
  ! * measurement spectrum used in ozone fitting process
  REAL (KIND=dp), DIMENSION (max_fit_pts):: fitwavs, fitweights, currspec
  REAL (KIND=dp), DIMENSION (max_fit_pts):: radwvl_sav
  !  REAL (KIND=dp), DIMENSION (2, max_fit_pts) :: curr_radresponse_spec
  ! * ring spectrum
  INTEGER                                       :: nsol_ring, sring_fidx,sring_lidx
  REAL (KIND=dp), DIMENSION (2, max_ring_pts)   :: sol_spec_ring
  ! ----------------------
  ! variable slit width
  ! ----------------------
  LOGICAL        :: slit_redo, wavcal_redo, wavcal_sol, wavcal, slit_rad
  LOGICAL        :: fixslitcal, smooth_slit, slitcal
  REAL (KIND=dp) :: slit_trunc_limit
  REAL (KIND=dp) :: dwavmax
  REAL (KIND=dp) :: sol_wav_avg, rad_wav_avg
  INTEGER  :: slit_fit_pts, n_slit_step, nslit,nslit_rad, nslit_sol, &
              wavcal_fit_pts, n_wavcal_step, nwavcal_sol, nwavcal_rad
  REAL (KIND=dp), DIMENSION (max_fit_pts) :: slitwav, slitwav_rad=0., &
       slitwav_sol=0.0, sswav_rad=0.0, sswav_sol=0.0, slitdis=0.0
  REAL (KIND=dp), DIMENSION (max_fit_pts, max_calfit_idx, 2) :: solslitfit=0.0,&
       radslitfit=0.0, solwavfit=0.0, radwavfit=0.0, slitfit=0.0
  CHARACTER (LEN=maxchlen)  :: slit_fname, rslit_fname, swavcal_fname,wavcal_fname
  ! hw1e, e_asym, shi, squ, hwl, hwr, vgl, vgr at each window
  !(value, standard deviation)
  REAL (KIND=dp), DIMENSION(maxwin)                   :: wincal_wav
  REAL (KIND=dp), DIMENSION(maxwin,max_calfit_idx, 2) :: solwinfit, radwinfit,solwinfit_save

  REAL (Kind=dp), DIMENSION(max_fit_pts) :: tmp_rad =0.0

  REAL (KIND=dp), DIMENSION (max_spec_pts) :: cubic_x, cubic_y, cubic_w
  REAL (KIND=dp), DIMENSION (max_spec_pts) :: poly_x, poly_y, poly_w
  INTEGER                                  :: poly_order

  !REAL (KIND=dp), DIMENSION (max_spec_pts) :: step2_y
  !REAL (KIND=dp), DIMENSION (max_spec_pts, max_fit_pts) :: step2_dyda
  !---------------------------------------------------
  ! fitting variables
  !---------------------------------------------------
  ! * Number of pixels read and number of pixels with successful retrieval
  INTEGER :: npix_fitting, npix_fitted
  ! * actual number of fitted variables and wavelengths
  INTEGER :: num_param    ! n_fitvar_rad - nfgas - (ozfit_end_index - ozfit_start_index + 1)
  INTEGER :: num_wav_max  ! maxval(irrad%nwav(:) - numwin*2*radnhtrunc  
  ! * calibration fit variables
  INTEGER                                       :: n_fitvar_rad, n_fitvar_sol
  INTEGER,           DIMENSION (max_calfit_idx) :: mask_fitvar_sol, rmask_fitvar_sol
  REAL (KIND=dp),    DIMENSION (max_calfit_idx) :: & 
       fitvar_sol, fitvar_sol_saved,fitvar_sol_init, & 
       lo_sunbnd, up_sunbnd, lo_sunbnd_init, up_sunbnd_init
  CHARACTER (LEN=6), DIMENSION (max_calfit_idx) :: fitvar_sol_str
  ! * geophysical fit variables 
  INTEGER,           DIMENSION (n_max_fitpars)  :: mask_fitvar_rad, rmask_fitvar_rad, &
                                                   database_indices, fothvarpos
  REAL (KIND=dp),    DIMENSION (n_max_fitpars)  :: & 
                                                   fitvar_rad, fitvar_rad_saved, &
                                                   fitvar_rad_init,fitvar_rad_init_saved, & 
                                                   fitvar_rad_apriori, fitvar_rad_aperror, &
                                                   fitvar_rad_std, fitvar_rad_nstd, &
                                                   lo_radbnd, up_radbnd, lo_radbnd_init, up_radbnd_init
  CHARACTER (LEN=6),  DIMENSION (n_max_fitpars)  :: fitvar_rad_str
  CHARACTER (LEN=15), DIMENSION (n_max_fitpars)  :: fitvar_rad_unit

  ! * spectrum input variables in ozone fitting process 
  REAL (KIND=dp),    DIMENSION (max_fit_pts)    :: & 
      fitspec_rad, &! I/F after removing non-ozone and albedo components
      simspec_rad, &! 
      clmspec_rad, &! simulated spectrum, ozone and albedo terms only, but with a priori climatology 
      fitres_rad, & ! fitspec_rad - simspec_rad
      actspec_rad   ! I/F without removing anything

  ! * geophysical input variables in ozone fitting process
  CHARACTER (LEN = 28)              :: the_utc
  INTEGER                           :: the_month, the_year, the_day, the_jday
  !REAL (KIND=dp), DIMENSION (3)     :: sza_atm, vza_atm, aza_atm
  !REAL (KIND=dp) :: avgsza, avgvza,avgaza, avgsca
  REAL (KIND=dp)                    :: the_sza_atm, the_vza_atm, the_aza_atm, &
       the_sca_atm, the_lat, the_lon, the_surfalt
  INTEGER                           :: nview, nloc
  REAL (KIND=dp), DIMENSION(maxloc) :: the_lons, the_lats
  REAL (KIND=dp), DIMENSION(2)      :: edgelons, edgelats 
 
  ! -------------------------------------
  ! Variables related to Air Mass Factors
  ! -------------------------------------
  INTEGER            :: n_amftab_dim, n_amftab_ang
  LOGICAL            :: have_amf, have_amftable
  REAL    (KIND=dp)  :: amfgeo, amf, sol_zen_eff
  REAL    (KIND=dp)  :: amf_esza_min, amf_esza_max
  REAL    (KIND=dp), DIMENSION (n_amftab_ang_max, n_amftab_dim_max) :: amf_table

  ! -----------------------------
  ! Previously IMPLICIT variables
  ! -----------------------------
  REAL (KIND=dp) :: phase, szamax, zatmos, chisq
  ! ------------------------------------------------------------------------
  ! Variables connected with ELSUNC numerical precision/convergence criteria
  ! ------------------------------------------------------------------------
  REAL (KIND=dp) :: tol,  epsrel,  epsabs,  epsx

  ! ----------------------------------------
  ! Variable for +1.0 or -1.0 multiplication
  ! ----------------------------------------
  REAL (KIND=dp) :: pm_one
  ! --------------------------------------
  ! Solar and Earth shine wavlength limits
  ! --------------------------------------
  !REAL (KIND=dp), DIMENSION (n_sol_winwav) :: sol_winwav_lim ! not used anywhere
  !REAL (KIND=dp), DIMENSION (n_rad_winwav) :: rad_winwav_lim ! not used anywhere

  ! ------------------------------------------------------------------
  ! Indices of fitting window defining wavelengths in current spectrum
  ! ------------------------------------------------------------------
  !INTEGER, DIMENSION (n_rad_winwav)  :: rad_winwav_idx ! not used anywhere


  !----------------------
  ! Write HDF output ?  0:ascii, 1:(GOME2), 2:(GOME2),  3:he5(OMI) 4:nc (TEMPO)
  !---------------------
  INTEGER :: l2_hdf_flag ! determined in main_control.inp
  CHARACTER(LEN=maxchlen)  :: outdir
  ! ---------------------
  ! L1B and L2 file names
  ! ---------------------
  CHARACTER (LEN=maxchlen) :: l1b_rad_filename, l1b_irrad_filename, l2_filename, &
       l2_cld_filename, l2_geos5_filename, l3_toc_filename
  CHARACTER (LEN=maxchlen), PARAMETER :: l2_swathname = 'O3 Profile'

  INTEGER, PARAMETER :: lcurve_unit = 50, ozwrtint_unit = 51, atmos_unit = 52, &
       tprofunit1 = 53, tprofunit2 = 54, profunit = 55,  &  ! two toms std profiles, ozprof
       cldunit = 56, corr_unit = 57, l1l2inp_unit = 58, ozabs_unit = 59, &
       calunit = 22, l2funit = 23, l2swathunit = 24, rtm_unit=25

  ! Writing limited output to screen for debugging and testing
  LOGICAL     :: scnwrt
  LOGICAL     :: calwrt ! for slit/wavelength calibration
  LOGICAL     :: calscn ! for slit/wavelength calibration
  LOGICAL     :: rtmdbg 
  
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
  INTEGER                                    :: n_mol_fit, n_fincol_idx
  INTEGER, DIMENSION (max_mol_fit)           :: fitcol_idx
  INTEGER, DIMENSION (2,max_mol_fit*mxs_idx) :: fincol_idx

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! satellite variables as a dimension of nxtrack, nytrack
  !  - bad pixel is filtered
  !  - sub-wave is taken in defined fitting window
  !  - co-adding is applied
  ! 1. irrad_group, rad_group, ring_group, refl_group 
  ! 2. cali_group
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ------------------------------------------------------------
  ! ozone fitting spectra
  ! ------------------------------------------------------------
  TYPE irrad_group
  INTEGER, ALLOCATABLE, DIMENSION (:)     :: errstat! nxtrack
  INTEGER, ALLOCATABLE, DIMENSION (:)     :: nwav   ! (nxtrack)
  INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: npix   ! (numwin, nxtrack)
  INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: winpix ! (numwin, nxtrack, 2)
  INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: wind   ! (max_fit_pts, nxtrack)
  INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: qflg   ! (max_fit_pts, nxtrack)
  REAL (kind=dp), ALLOCATABLE, DIMENSION (:,:) :: wavl,spec,prec !(max_fit_pts,nxtrack)
  REAL (kind=dp), ALLOCATABLE, DIMENSION (:) :: norm !(nxtrack)
  ENDTYPE irrad_group

  TYPE rad_group
  INTEGER, ALLOCATABLE, DIMENSION (:)     :: errstat
  INTEGER, ALLOCATABLE, DIMENSION (:,:)     :: pix_errstat
  INTEGER, ALLOCATABLE, DIMENSION (:,:)     :: nwav   ! (nxtrack, nytrack)
  INTEGER, ALLOCATABLE, DIMENSION (:,:,:)   :: npix   ! (numwin, nxtrack,nytrack)
  INTEGER, ALLOCATABLE, DIMENSION (:,:,:,:) :: winpix ! (numwin, nxtrack,nytrack 2)
  INTEGER, ALLOCATABLE, DIMENSION (:,:,:)   :: wind   ! (max_fit_pts, nxtrack)
  INTEGER, ALLOCATABLE, DIMENSION (:,:,:)   :: qflg   ! (max_fit_pts, nxtrack, nytrack)
  REAL (kind=dp), ALLOCATABLE, DIMENSION (:,:,:)      :: wavl,spec,prec !(max_fit_pts,nxtrack,nytrack)
  REAL (kind=dp), ALLOCATABLE, DIMENSION (:,:) :: norm !(nxtrack, nytrack)
  END TYPE rad_group
  ! ------------------------------------------------------------
  ! tempo ring fitting spectra
  ! ------------------------------------------------------------
  TYPE ring_group
    REAL    (KIND=dp), ALLOCATABLE, DIMENSION (:,:) ::  spec, wavl !(max_ring_pts, nxtrack_max)
    INTEGER, ALLOCATABLE, DIMENSION (:)   :: nsol,ndiv ! (nxtrack_max)
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: winpix ! (nxtrack_max,2)
  END TYPE ring_group

  !-----------------------
  ! tempo surface albedo fitting spectra
  !-----------------------
  TYPE refl_group
    REAL (KIND=r4), ALLOCATABLE, DIMENSION (:,:)   :: solspec, solwavl ! (mrefl, nxtrack)
    REAL (KIND=r4), ALLOCATABLE, DIMENSION (:,:,:) :: radspec, radwavl  !(mrefl, nxtrack, ntimes)
    INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: winpix ! (nxtrack, 2)
  END TYPE refl_group

  !-----------------------------------------
  ! geolocation 
  !-----------------------------------------
  TYPE geo_group
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)       :: time
    REAL (KIND=r4), ALLOCATABLE, DIMENSION (:,:)    :: lon, lat, sza, vza, aza, sca  !  (mrefl, nxtrack)
    REAL (KIND=r4), ALLOCATABLE, DIMENSION (:,:,:)    :: clon, clat  !  (mrefl, nxtrack)
    REAL (KIND=r4), ALLOCATABLE, DIMENSION (:,:)    :: elon, elat  !  (mrefl, nxtrack)
    REAL (KIND=r4), ALLOCATABLE, DIMENSION(:,:)     :: cfr, ctp, ai
    integer (kind=4), allocatable, dimension(:)     :: step_idx    ! (ntimes)
    INTEGER (KIND=1), ALLOCATABLE, DIMENSION(:,:)   :: cloud_qflg
    INTEGER (KIND=4), ALLOCATABLE, DIMENSION(:,:)   :: gflg
    INTEGER (KIND=1), ALLOCATABLE, DIMENSION(:,:)   :: xflg
    INTEGER (KIND=2), ALLOCATABLE, DIMENSION(:,:) :: height
    INTEGER (KIND=2), ALLOCATABLE, DIMENSION(:,:) :: land_water_flg
    INTEGER (KIND=2), ALLOCATABLE, DIMENSION(:,:) :: snow_ice_flg
    INTEGER (KIND=2), ALLOCATABLE, DIMENSION(:,:) :: glint_flg
  END TYPE geo_group

  !-------------------------------------------
  ! calibration variables for solar spectra
  !-------------------------------------------
  TYPE cali_group
     ! calibration results for each fitting window
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:)      :: wincal_wav     !maxwin, nxtrack
     REAL (KIND=dp), ALLOCATABLE, DIMENSION (:,:,:)   :: solwinfit, radwinfit ! maxwin, max_calfit_idx, 2, nxtrack
     ! slit fit results for sub-wavelengths
     INTEGER, ALLOCATABLE, DIMENSION (:)              :: nslit_sol, nslit_rad
     REAL (KIND=dp), ALLOCATABLE, DIMENSION (:,:,:,:) :: slitfit_sol, slitfit_rad
     REAL (KIND=dp), ALLOCATABLE, DIMENSION (:,:)     :: slitwav_sol, slitwav_rad
     ! wave fit results for sub
     INTEGER, ALLOCATABLE, DIMENSION (:)              :: nwavcal_sol, nwavcal_rad ! nxtrack_max
     REAL (KIND=dp), ALLOCATABLE, DIMENSION (:,:)     ::sswav_sol, sswav_rad !max_fit_pts, nxtrack
  END TYPE cali_group

  !------------------------------------------------------------
  ! retrieval status for a block
  !------------------------------------------------------------
  TYPE o3p_group
  INTEGER, DIMENSION (:, :), ALLOCATABLE :: exitval, initval
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: fitvar
  END TYPE  o3p_group
  
END MODULE OMSAO_variables_module


