! ***********************************************************
! Author:  xiong liu
! Date  :  July 23, 2003
! Purpose: Data module for ozone profile retrieval
!          Used for sharing data across different subroutines
! ***********************************************************

MODULE ozprof_data_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxchlen, maxlay, max_fit_pts, &
       max_iter, mrefl, maxwin, mflay, max_spec_pts , mreflcld
  USE OMSAO_indices_module,    ONLY : n_max_fitpars,  max_rs_idx, &
       maxoth, maxalb, maxwfc,  shift_offset, &
       no2_t1_idx, no2_t2_idx, o2o2_idx, so2_idx, bro_idx, oclo_idx,     &
       hcho_idx, o2_idx, h2o_idx, so2v_idx, bro2_idx,o2t2_idx,h2ot2_idx, &
       hwe_idx, spk_idx, asy_idx       
 
  IMPLICIT NONE

  ! ---------------------------------------------------
  ! filename for input/output profile variables
  ! ---------------------------------------------------
  CHARACTER (LEN=maxchlen) :: algorithm_name
  CHARACTER (LEN=maxchlen) :: algorithm_version
  CHARACTER (LEN=maxchlen) :: ozprof_input_fname   ! profile input file (unit=11)
  CHARACTER (LEN=maxchlen) :: atmos_prof_fname     ! atmos. profiles for LIDORT (unit=55)
  CHARACTER (LEN=maxchlen) :: ozabs_fname          ! o3 T-depend absorption coefficents
   
  ! ---------------------------------------
  ! whether to write retrieval at each step
  ! very useful during developing process
  ! ---------------------------------------
  LOGICAL                       :: ozwrtint
  CHARACTER (LEN=maxchlen)      :: ozwrtint_fname    

  ! Whether to write correlation maxtrix, covariance matrix, contribution function, 
  ! and residual, atmospheric condition used in radiative transfer, and trace gases
  LOGICAL :: ozwrtvar, ozwrtcorr, ozwrtcovar, ozwrtavgk, ozwrtfavgk, ozwrtcontri, &
       ozwrtres, atmwrt, gaswrt, radcalwrt, ozwrtwf, ozwrtsnr, wrtring,wrtozcrs, wrtalbspc

  ! Options and ozone profiles used in performing radiometric calibration
  INTEGER                       :: which_caloz ! 1 climatology 2 true ozone profile
  CHARACTER (LEN=maxchlen)      :: caloz_fname    
  ! ----------------------------------------------
  ! Variables for regularization
  ! ----------------------------------------------
  LOGICAL                   :: lcurve_write         ! write Lcurve output
  CHARACTER (LEN=maxchlen)  :: lcurve_fname         

  ! 1: lcurve, 2: GCV, 3: use GSV with lcurve validation (GCV sometimes fails)
  INTEGER                   :: lcurve_gcv   
  ! --------------------------------------------------------
  ! Use logarithmic radiance / logarithmic state vectir
  ! --------------------------------------------------------
  LOGICAL                       :: use_lograd 
  LOGICAL                       :: use_logstate  

  ! --------------------------------------------------------
  ! Use optimal estimation / phillips-tikhonov regularization
  ! --------------------------------------------------------
  LOGICAL                       :: use_oe   

  ! --------------------------------------------------------
  ! For trace gases like SO2, which has both natural and 
  ! anthrologenic sources that varies a lot, using a priori
  ! to consrain the retrievals will significantly underestimate
  ! volcanic SO2. One way to avoid this is use a second step 
  ! (lienar retrieval) after each iteration of ozone retrieval.
  ! --------------------------------------------------------
  LOGICAL                      :: do_twostep
  LOGICAL                      :: do_bothstep
  LOGICAL                      :: use_large_so2_aperr

  ! ------------------
  ! Use Floor Noise
  ! ------------------
  LOGICAL                       :: use_flns

  ! -------------------------------------
  ! whether to retrieve ozone profile
  ! -------------------------------------
  LOGICAL                       :: ozprof_flag

  ! -------------------------------------
  ! South Atlantic Amomaly Flag
  ! -------------------------------------
  LOGICAL                       :: saa_flag
  INTEGER                       :: nsaa_spike

  ! Need to convolve high-resolution ozone absorption cross section 
  ! at the beginning of processing each scan position of each block
  ! Once the xsection is convolved, it will be set to false in ROUTINE getabs_crs
  LOGICAL                       :: ozabs_convl, so2crs_convl, &
                                   o4crs_convl,  o2crs_convl, h2ocrs_convl

  !----------------------------------------------------------
  ! Temperatrue dependent cross section for other tracegases
  !----------------------------------------------------------
  LOGICAL                       :: use_so2dtcrs
  LOGICAL                       :: use_o4dtcrs
  LOGICAL                       :: use_o2dptcrs
  LOGICAL                       :: use_h2odptcrs

  !-----------------------------------------------------------
  ! indices for set up atmosphere
  ! --------------------------------------------------------
  !INTEGER                       :: which_atmos      
  ! use what atmospheric profiles
  ! 1: use interpolated TOMS V7 profiles
  ! 2: use GOME ozone profile workgroup a priori dataset
  ! 3. Best auxilary information
  ! --------------------------------------------------------
  INTEGER                       :: which_atm 
  ! 0 NCEP 1 FNL 2 GEOS5 (OMI)
  INTEGER                       :: which_clima         
  !1. V8+EP
  !2. McPeters Clima 
  !3. McPeters+GEOS-CHEM   18 x 12
  !4. McPeters+GEOS-CHEM   (144 x 91)
  !5. McPeters+LOGAN CLIMA (72 x 46)
  !6. McPeters+zonal mean MLS (0.1 mb - 215 mb)
  !7. McPeters+individual MLS / other profiles
  !8. TB
  !9. AB
  !10 ML

  INTEGER                       :: which_aperr
  !1: McPeters Clima  
  !2: Fortuin
  !3. McPter Clima+GEOS-CHEM  
  !4. MLS zonal mean  
  !5. MLS/other individual profile  
  !8. TB
  !9  AB
  !10  ML

  LOGICAL                       :: loose_aperr         
  ! Whehter to loose a priori constraint
  REAL (KIND=dp)                :: min_serr, min_terr  
  ! Minimum relative a priori error in the statosphere and troposphere     

  INTEGER                       :: which_toz         
  ! 0: None 
  ! 1: gridded (e.g. TOMS 1.25 x 1)      
  ! 2. zonal mean      
  ! 3. Individual (not implemented)
  LOGICAL                       :: norm_tropo3  
  !T: normalize trop O3 only  F: normalize whole profile
  INTEGER                      :: which_alb
  ! 1: GOME, 2: TOMS 3 4 KNMI-OMI
  INTEGER                      :: which_cld
  ! 1: , 2: climatology 5 self-retrievals

  ! variables to control atmospheric layering scheme
  LOGICAL :: use_reg_presgrid    ! T: automatically generated F: read from an file
  LOGICAL :: use_tropopause      ! T: varying tropopause F: fixed (e.g., 200 mb)
  LOGICAL :: fixed_ptrop         ! T: fixed tropopopause based on fixed pressure, F: based on pressure layer (from bottom)
  LOGICAL :: adjust_trop_layer   ! T: redistribute layers between tropopause and surface F: not except for surface
  LOGICAL :: insert_sfc_layer=.false. ! 
  !T: insert a small surface layer of 10 hPa, not used anywhere
  CHARACTER (LEN=maxchlen) :: presgrid_fname  ! File contains the pressure levels

  LOGICAL        :: smooth_ozbc !for processing ozone below clouds
  ! T: update ozone below clouds through smoothing 
  ! F: Don't update ozone below clouds (need to modify a priori ozone variance matrix

  !-----------------------------------------------------------------
  ! state vector set up variables
  !------------------------------------------------------------------
  INTEGER :: nlay, nlay_fit, nt_fit 
  ! # of layers, # of layers in fine-grid LIDORT calculation, 
  ! # of layers varied for ozone, and # of layers varied for T

  ! Starting and ending indices for fitting variables for convenient acess
  INTEGER :: ozfit_start_index, ozfit_end_index, tf_fidx, tf_lidx

  ! Indices in actual fitting array
  INTEGER :: ozprof_start_index, ozprof_end_index, t_fidx, t_lidx
  ! assumed that varied layers are continuous
  INTEGER :: start_layer, end_layer     ! First and last varied layers in o3 array  

  ! ---------------------------------------------
  ! Variables for surface and tropopause pressure
  ! ---------------------------------------------
  REAL (KIND=dp) :: ps0, pst, pst0, trpz 

  ! Number of group of auxiliar parameters
  INTEGER                             :: nothgrp  
  ! Number of fitted variabled for o3 crs shift, sol/rad shi, sol/rad slit (0-4)
  ! first-order ring (1/0)         
  INTEGER :: nos, nsh, nsl, nrn, ndc, nis, nir, np1, np2, ncm
 
  ! additional variables to set up inr
  INTEGER                        :: which_inr
  !-------------------------------------------
  ! inr spectrum 
  ! 0) from mean fitting residuals
  ! 1) from OPT radiance sensitivity function
  !--------------------------------------------

  INTEGER, DIMENSION (maxwin, maxoth) :: &
       osind, osfind, slind, slfind, shind, shfind, &
       rnind, rnfind, dcind, dcfind, isind, isfind, &
       irind, irfind, p1ind, p1find, p2ind, p2find, & 
       cmind, cmfind
  INTEGER, DIMENSION (maxoth, 2)      :: &
       oswins, slwins, shwins, rnwins, dcwins, iswins, irwins, &
       p1wins, p2wins, cmwins

  INTEGER :: ncldaer, ecfrind, ecodind, ectpind, taodind, twaeind, saodind, &
       ecfrfind, ecodfind, ectpfind, taodfind, twaefind, saodfind, &
       sprsind, sprsfind, so2zind, so2zfind

  ! atmospheric profiles (P, T, Z, O3, and other trace gases) 
  REAL (KIND=dp), DIMENSION (3, 0:maxlay)   :: atmosprof   ! 1: P, 2: Z  3: T 

  ! Ozone, target gases, all information ( initital, a priori) will be saved
  REAL (KIND=dp), DIMENSION (maxlay)        :: ozprof, ozprof_init, &
       ozprof_std, ozprof_ap, ozprof_apstd, ozprof_nstd
  INTEGER, DIMENSION(0:maxlay)  :: nup2p
  INTEGER                       :: ntp, ntp0    ! tropopause layer
  INTEGER                       :: nsfc, nfsfc  ! surface layer in retrieval and find altitude grids
  ! ---------------------------------------------------------
  ! Variables for updating ozone and sao3 at first iteration
  ! ---------------------------------------------------------
  REAL(KIND=dp)  :: ozone_above60km
  !REAL(KIND=dp), DIMENSION (0:11) :: pv811

  ! Other traces gases (: NO2, SO2, BrO, HCHO)   
  INTEGER, PARAMETER :: ngas = 13, nallgas = 14
  INTEGER            :: nfgas
  INTEGER, DIMENSION(ngas), PARAMETER         :: & 
   gasidxs = (/no2_t1_idx, no2_t2_idx, &
              o2_idx, o2t2_idx, o2o2_idx,  bro_idx, bro2_idx, h2o_idx, h2ot2_idx, &
              so2_idx, so2v_idx, hcho_idx, oclo_idx/)
  INTEGER :: so2idx, so2vidx, o4idx, o2idx, o2t2idx, h2oidx, h2ot2idx, broidx,hchoidx, no2idx
  INTEGER :: so2fidx,so2vfidx,o4fidx,o2fidx,o2t2fidx,h2ofidx,h2ot2fidx,brofidx, hchofidx,no2fidx
  INTEGER :: o3crsidx, so2crsidx, o4crsidx, o2crsidx, h2ocrsidx ! pos in allcol or allcrs (1:o3, 2~~)
  LOGICAL, DIMENSION(max_rs_idx)              :: rtm_treatment   ! included in RTM calculation
  INTEGER, DIMENSION(ngas), PARAMETER         :: gassidxs = gasidxs + shift_offset
  REAL (KIND=dp), DIMENSION (ngas, mflay)     :: mgasprof = 0.0D0
  INTEGER, DIMENSION(ngas)                    :: fgasidxs, fgassidxs, fgaspos
  ! initial, a priori, a priori std, retrieved, uncertainty (S+N), uncertainty (N), air mass factor, 
  ! above cloud fraction, average kernel (a priori influence), avgk (consider influence from others)
  LOGICAL                       :: do_tracewf
  REAL (KIND=dp), DIMENSION(ngas, 10)          :: tracegas = 1.0  ! initial AMF to 1.0

  ! SO2V profiles if the cental altitude is decreased/increased by 1 km
  REAL (KIND=dp), DIMENSION (mflay, 2)         :: so2vprofn1p1 = 0.0D0
  REAL (KIND=dp), DIMENSION (-1:1)             :: so2valts

  ! ------------------------------------------
  ! Variables for aerosols and clouds
  ! ------------------------------------------
  ! cloud: true, need to consider cloud but it may be clear 
  ! and cloud info is not available
  ! has_clouds: true, all necessary clouds info. is there
  LOGICAL        :: aerosol, strat_aerosol, cloud, has_clouds, useasy, &
       do_lambcld, do_aerosols
  INTEGER        :: which_aerosol
  LOGICAL        :: scale_aod
  REAL (KIND=dp) :: lambcld_refl, scacld_cod, scaled_aod, lambcld_initalb, scacld_initcod

  ! glint variables: in the presence of glint, increase the a priori error for surface albedo
  LOGICAL        :: has_glint
  REAL (KIND=dp) :: glintprob
  INTEGER        :: the_snowice, the_landwater_flg, the_glint_flg


  ! ------------------------------------------------------------------------------
  ! variables for albedo: albedo can be at most quadratically wavelength-dependent
  ! ------------------------------------------------------------------------------
  LOGICAL                   :: vary_sfcalb
  REAL (KIND=dp)            :: pf2ba0, pf2ba1, pf2fc0, pf2fc1
  LOGICAL                   :: do_alb_longwav=.false., use_prefitalb=.false. ! used for derving alb from non-ozoneband
  REAL (KIND=dp)            :: alb_swav, alb_ewav!, cld_swav, cld_ewav
  REAL (KIND=dp)            :: pos_alb, toms_fwhm 
  REAL (KIND=dp)            :: measref
  CHARACTER(LEN=maxchlen)   :: alb_tbl_fname, ozcrs_alb_fname       ! reflectance table
  ! @ variables for EOF/BRDF albedo spectrum
  LOGICAL                   :: use_albspc=.false.   ! Use albedo spectra/EOFs (bands 3&4only)
  LOGICAL                   :: use_albeofs=.false.  ! Use albedo EOFs (use_albspec must be set)
  INTEGER, PARAMETER        :: malbspc = 5  ! maximum number of albedo spectra/EOFs
  INTEGER                   :: nalbspc      ! # albedo spectra/# EOFs (does not include mean spectrum for albedo EOFs)
  INTEGER                   :: nactalbspc   ! # actually used albedo spectra/EOFs (include mean)
  INTEGER                   :: nalbspcwin   ! number of windows for usign albedo spectra
  INTEGER                   :: nalbspcord   ! # of parameters for scaling albedo spectra/EOFs
  REAL (KIND=dp), DIMENSION(max_fit_pts, 0:malbspc-1)  :: albspcs ! 0 store spectrum or mean EOF
  REAL (KIND=dp), DIMENSION(max_fit_pts, 2)            :: sfcalbs ! Initial andderived surface albedo values
  LOGICAL, DIMENSION (maxalb)       :: is_albspcvar

  ! number of wavelengths around center wavelength for calculating cloud fraction
  INTEGER :: nrefl
  ! wavelengths and spectra in solar/earthshine for computing reflectance
  REAL (KIND=dp), DIMENSION(mrefl) :: rad_posr, rad_specr, sun_posr, sun_specr 

 ! Number of wavelength for cloud retrieval cai
  INTEGER                          :: nreflcld
  ! wavelengths and spectra in solar/earthshine for cloud retrieval
  REAL(KIND=dp), DIMENSION(mreflcld) :: rad_posc, rad_spec, sun_posc, sun_specc

  ! wave-variant albedo
  INTEGER :: nalb, nfalb  ! total # of albedoes and fixed albedoes
  REAL (KIND=dp), DIMENSION(maxalb) :: eff_alb, eff_alb_init ! albedoes
  REAL (KIND=dp), DIMENSION(maxalb) :: albmax, albmin   ! max and min lamda 
  INTEGER, DIMENSION (maxalb)       :: albfpix, alblpix ! first and last pixel
  INTEGER :: albidx, albfidx, thealbidx ! star index for albedo in whole and varied array

  ! Wavelength-dependent cloud fraction (wfc) terms
  INTEGER :: nwfc, nfwfc                                     ! total # of used and varied wfc
  REAL (KIND=dp), DIMENSION(maxwfc) :: eff_wfc, eff_wfc_init ! wfc terms
  INTEGER, DIMENSION (maxwfc)       :: wfcfpix, wfclpix      ! first and last pixel
  INTEGER :: wfcidx, wfcfidx, thewfcidx ! star index for wfc in whole and varied array
  REAL (KIND=dp), DIMENSION(maxwfc) :: wfcmax, wfcmin   ! max and min wfc

  ! Aerosol information
  INTEGER, PARAMETER                 :: maxawin = 6
  INTEGER                            :: actawin
  REAL (KIND=dp), DIMENSION(maxawin) :: aerwavs
  REAL (KIND=dp), DIMENSION(maxawin) :: strataod, stratsca, tropaod, &
       tropsca, stratwaer, tropwaer

  ! Number of legendre phase moments and number of greek scattering expansion coefficients
  INTEGER, PARAMETER                        :: maxmom = 64, maxgksec = 6, maxgkmatc = 8
  INTEGER                                   :: nmom, ngksec, ngkmatc

  ! atmospheric data on fine grids (save it, to avoid repetitive reading)
  INTEGER                                   :: nflay, ncbp, nctp
  REAL (KIND=dp), DIMENSION(0:mflay)        :: fts, fps, fzs, fozs, frhos
  REAL (KIND=dp), DIMENSION(maxawin, mflay) :: gaext, gasca, gaasy
  REAL (KIND=dp), DIMENSION(maxawin)        :: gcq, gcw, gcasy
  REAL (KIND=dp), DIMENSION(maxawin, mflay, 0:maxmom, maxgksec) :: gamoms
  REAL (KIND=dp), DIMENSION(maxawin, 0:maxmom, maxgksec)        :: gcmoms

  REAL (KIND=dp) :: the_cfrac, the_ctp, the_cod, the_cbeta, &
       the_orig_cfr, the_orig_ctp, the_orig_cod, the_ai 
  INTEGER        :: the_cld_flg

  INTEGER, PARAMETER :: maxloc =5
  TYPE the_geo
    REAL(KIND=dp) :: sza, vza, aza, sca, lat, lon,surfalt, &
                   lons(maxloc), lats(maxloc),elon(2), elat(2), &
                   cfrac, ctp
  END TYPE  the_geo  
  TYPE (the_geo) :: the_geo1, the_geo2

  ! ---------------------------------------------------------
  ! Variables for whehter performing radiance calculation
  ! using effective cross section or high resolution 
  ! ---------------------------------------------------------
  ! Parameters for T-dependent cross sections
  INTEGER, PARAMETER :: mxsect    =  5
  LOGICAL            :: use_effcrs
  INTEGER, PARAMETER :: radc_msegsr = 5
  INTEGER            :: radc_nsegsr 
  INTEGER            :: nhresp, ncalcp, nhresp0
  REAL(KIND=dp), DIMENSION(radc_msegsr)          :: radc_samprate, radc_lambnd
  REAL(KIND=dp)                                  :: hres_samprate
  INTEGER, DIMENSION(max_fit_pts)                :: radcidxs
  REAL(KIND=dp), DIMENSION(max_fit_pts)          :: radcwav
 ! Rayleigh scattering coefficients and molecular depolarization correction factor
  REAL (KIND=dp), DIMENSION(max_fit_pts) :: raycof, depol 

  REAL(KIND=dp), DIMENSION(max_spec_pts)         :: hreswav, hreswav0, hres_i0, hres_raycof, hres_depol
  REAL(KIND=dp), DIMENSION(mxsect, max_spec_pts) :: hres_o3,  hres_o3shi, hres_so2, hres_so2shi, &
                                                    hres_o4, hres_o4shi
  REAL(KIND=dp), DIMENSION(mflay, max_spec_pts)  :: hres_o2, hres_h2o, hres_o2shi, hres_h2oshi
  REAL(KIND=dp), DIMENSION(ngas, max_spec_pts)   :: hres_gas, hres_gasshi
  REAL(KIND=dp), DIMENSION(max_spec_pts, mflay)  :: o3crsz,so2crsz,o4crsz,o2crsz,h2ocrsz
  REAL(KIND=dp), DIMENSION(max_spec_pts, mflay)  :: o3dadsz,o3dadtz
  REAL(KIND=dp), DIMENSION(max_spec_pts, mflay)  :: hresgabs, hresray      
  REAL(KIND=dp), DIMENSION(max_spec_pts)         :: so2dads, o4dads, o2dads, h2odads ! Weighted by profiles

  ! --------------------------------------
  ! Variables for ring effect
  ! --------------------------------------
  LOGICAL        :: ring_on_line, ring_convol, fit_atanring , ring_LUT

  ! 0: zero order, 1: first order, 2: second order, 3: x-avg(x) 
  ! 4: 0/1/2 combination, 5: both T and O3 (zero order) 6: use oe like
  INTEGER         :: ptr_order 

  ! only useful when ptr_order = 4
  REAL (KIND=dp)  :: ptr_w0, ptr_w1, ptr_w2

  ! ---------------------------------------------
  ! Variables for degradation correction
  ! Variables for bias correction (in GOME data)
  ! ---------------------------------------------
  LOGICAL                 :: degcorr, biascorr
  CHARACTER(LEN=maxchlen) :: degfname, biasfname     ! reflectance table
  INTEGER                 :: which_biascorr

  ! --------------------------------------
  ! Variables for polarization correction
  ! --------------------------------------
  INTEGER        :: VlidortNstream, polcorr

  ! -----------------------------------------
  ! Do claculation at three VZAs (A, B, C), but
  ! still use the effective SZA
  ! -----------------------------------------
  LOGICAL        :: do_multi_vza
  LOGICAL        :: do_simu, do_simu_rmring ! simulation, and remove ring effect in simulation
  REAL           :: the_fixalb

  ! ------------------------------------------------------
  ! Do radiance calculation at selected wavelengths and 
  ! perform interpolations for the others
  ! -----------------------------------------------------
  LOGICAL        :: do_radinter

  ! -----------------------------------------------------
  ! same parameters are derived for each different window
  ! -----------------------------------------------------
  LOGICAL        :: do_subfit

  ! ----------------------------------------------------------------------
  ! Whehter to coadd after b1a/b boundary change from 307 to 282 nm if 
  ! the starting wavelength is larger than the new boundary 282 nm.  If 
  ! false, the retrieval is done at 320 x 40 km2 and backscan is not used.
  ! OtherwISE, it is consistent with before the b1a/b change
  ! ---------------------------------------------------------------------
  LOGICAL    :: coadd_after_b1ab ! not used for OMI, maybe it is only for gome
  LOGICAL    :: b1ab_change ! not used for OMI
  LOGICAL    :: scia_coadd ! not used for OMI

  ! # of division from retrieval grid to find grid LIDORT calculation
  INTEGER :: ndiv   

  ! --------------------------------------
  ! Variables for reflectance calibration
  ! --------------------------------------
  ! spectra normalization constant, used for getting absolute reflectance
  REAL(KIND=dp)  :: div_sun, div_rad  
  REAL(KIND=dp), DIMENSION(maxwin)  :: div_suns, div_rads  

  ! ---------------------------------------------------------
  ! Measurement error matrix
  ! ---------------------------------------------------------
  REAL(KIND=dp)                          :: merr_corrlen   ! systematic error correlation lenth (nm)
  REAL(KIND=dp), dimension(max_fit_pts)  :: msyserr        ! systematic error
  LOGICAL                                :: do_sy_diagonal
  REAL(KIND=dp), DIMENSION(max_fit_pts, max_fit_pts) :: merr_covar  

  !------------------------------------------------------------------------------
  ! Fit resuls variables
  !------------------------------------------------------------------------------ 
  ! Covariance Matrix, Conbtribution Function, Average kernel
  REAL (KIND=dp), DIMENSION(n_max_fitpars, n_max_fitpars) :: avg_kernel, favg_kernel
  REAL (KIND=dp), DIMENSION(n_max_fitpars, max_fit_pts)   :: contri
  REAL (KIND=dp), DIMENSION(max_fit_pts, n_max_fitpars)   :: weight_function
  REAL (KIND=dp), DIMENSION(n_max_fitpars, n_max_fitpars) :: covar, ncovar
  REAL (KIND=dp) :: ozdfs, ozinfo ! DFS and information content for ozone
  INTEGER        :: num_iter      ! Current numer of iterations
  ! xliu, 08/06/2010: Add trace gas averaging kernels, profile weighting
  ! function, contribution function (VCD)
  ! trace gas averaging kernel (dC/dx): sensitivity of retrieved trace gas
  ! vertical column
  ! to actual trace gas at each individual layer
  ! profile weighting function: dy/dx
  ! contribution function: dC/dy
  ! dy = dlnI if use logarithm of measurements else dy = dI
  REAL (KIND=dp), DIMENSION(ngas, maxlay)              :: trace_avgk
  REAL (KIND=dp), DIMENSION(ngas, max_fit_pts, maxlay) :: trace_profwf
  REAL (KIND=dp), DIMENSION(ngas, max_fit_pts)         :: trace_contri
  REAL (KIND=dp), DIMENSION(ngas, maxlay)              :: trace_prof   ! Retrieval grid
  REAL(KIND=dp), DIMENSION(0:maxwin) :: allrms, allradrms


END MODULE ozprof_data_module
