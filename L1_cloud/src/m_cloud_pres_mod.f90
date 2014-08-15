module m_cloud_pres_mod

  !!-------------------------------------------------------------------
  !
  ! Variables and parameters used by m_cloud_pres_ret and related code
  !
  ! all variable descriptions added by E. O'Sullivan, 13-15Aug14,
  ! some guesswork involved.
  !
  !!-------------------------------------------------------------------

  INTEGER(KIND=4), EXTERNAL :: OMI_SMF_setmsg
  integer :: status, NISE 
  real (KIND=8), parameter :: sz_max=86d0 !max solar zenith
  real (KIND=8), parameter :: bad_thresh=0.02d0 !maximum acceptable
  !fractional error on observed flux
  real (KIND=8), parameter :: refl_cld_mask=0.4d0 !reflectance threshold
  !above which pixels are added to cloud mask 
  real (KIND=8), parameter :: cld_mask_press_diff=-0.2d0 !threshold for
  !identifying regions of ice & snow from unphysical pressure results,
  !though this check is commented out of m_cloud_pres_mod
  integer :: nbad !number of pixels with bad fluxes (large errors)
  integer, allocatable, dimension(:) :: bad_obs, good_obs !vectors 
  !identifying pixels with good or bad flux measurements
  integer :: nterms=2 !2
  integer :: nst !number of terms in state vector (=nterms +2,3 or 4)
  integer :: ntm !iteration index used with nterms
  real (KIND=8), dimension(:,:), allocatable :: x, x_fg !state vector
  ! and first guess value of state vector
  integer :: iter,nwave !iter=iteration counter for chi^2 fit
  ! nwave=number of wavelengths in Ring table, 
  integer :: nobs=0 !number of data points in usable wavelength range
  integer :: ix1, ix2, ix1_old, ic1_old, ic1, ic2, ierr !ix=cloud pressure
  ! index, ic=chlorophyll index
  integer :: i0x1, i0x2, i_np0 !i0x=surface pressure index(?), i_np0=
  ! integer value of lowest pressure value np0
  integer :: ip, i, ii, jj !counters, ii & jj are general purpose, 
  ! ip is cross-track pixel, i is for usable wavelength range 
  integer :: i1_ler,i2_ler !indices for lambert-equivalent relectance data
  real (KIND=8) :: int1_ler,int2_ler !interpolated reflectances 
  real (KIND=8) :: noise=0.005d0 !used in constructing covariance
  real (KIND=8) :: noise2=0.10d0 !used in constructing covariance
  real (KIND=8) :: delx=0.00001d0  !used in bracketing cloud press. index 
  real (KIND=8) :: delc=0.00001d0  !used in bracketing chlorophyll index
  real (KIND=8) :: diff_chi !fractional change in ch^2 from last iteration
  real (KIND=8) :: chisq_old, chisq, bias, std, bias_old !chi^2 values, 
  ! radian bias values, radiance standard deviation
  real (KIND=8) :: ixd, icd, np, np0, nc
  ! np0 = interpolated index of lowest pressure 
  ! np = interpolated index of current pressure estimate
  ! ixd = difference between interpolated index and current index value
  ! icd = as for ixd, but in chlorophyll retrieval
  ! nc = interpolated index of current chlorophyll concentration est.
  real (KIND=8) :: nt, j, l, nt_o, j_o, l_o ! interpolated values of 
  ! viewing geometry, nt=solar zenith, j=sat. zenith, l=azimuth, o=ocean
  real (KIND=8) :: refl_oc !ocean reflectivity
  integer :: ngood !number of good solar irradiances
  integer, dimension(:), allocatable :: good !vector identifying wavelengths
  ! for which solar irradiances are good
  real (KIND=8), dimension(3) :: res3 !polynomial fit to spectrum
  real (KIND=8), dimension(3) :: res !dummy result variable
  real (KIND=8), dimension(:), allocatable :: chls !chlorophyll concentration
  ! axis of ocean raman scattering table
  integer, dimension(:), allocatable :: ind ! vector of cross-track indices 
  ! for pixels in valid wavelength range
  real (KIND=8), dimension(:), allocatable :: y_obs, y_obs1, y_calc, &
       y_frac, y_calc_sh, y_calc_squeeze, sflx, wave_diff
  ! y_obs = flux in pixel     y_obs1 = pi*flux in pixel
  ! y_calc = model flux in pixel     y_frac = fractional residual
  ! y_calc_sh, y_calc_squeeze = flux after spectral shift or squeeze
  ! sflx = solar flux    
  ! wave_diff = change in wavelength from spectral shift or squeeze
  real (KIND=8), dimension(:), allocatable :: waves, wavesd  ! waves=vector 
  ! of valid wavelengths, wavesd=waves-waves(0)
  real (KIND=8), dimension(:,:), allocatable :: wavesp !wavesd^ntm, so 
  ! vectors of wavelengths squared, cubed, etc, depending on number of
  ! terms in state vector
  real (KIND=8), dimension(:,:), allocatable :: y_resid, y_back
  ! y_resid = residual flux in pixel  y_back = background flux 
  real (KIND=8), dimension(:), allocatable :: ycalc, rad_tot
  ! ycalc = model flux in pixel before squeeze/shift corrections
  ! rad_tot = total radiance of scene
  real (KIND=8), dimension(:), allocatable :: rad_clr, rad_cld, &
       ring_clr, ring_cld ! radiance & Ring effect fill-in in cloudy and
  !                         clear scenes (?)
  real (KIND=8), dimension(:), allocatable :: jacob_dummy, jacob_rad, fit_rad
  ! dummy jacobian, jacobian for computed radiance, linear radiance gradient 
  ! for subtraction.
  real (KIND=8), dimension(:), allocatable :: r_i, b_i
  ! r_i=observation errors,      b_i=background errors
  real (KIND=8), dimension(:,:), allocatable :: h, htr, err_cov, corr
  ! h=jacobian,    htr=jacobian transpose, 
  ! err_cov=covariance matrix,     corr=error correlatons
  real (KIND=8), dimension(:), allocatable :: ring_oc, rad_clr_oc
  ! ring effect fill-in and clear sky radiance in ocean scenes
  real (KIND=8), dimension(:,:), pointer :: ring_clds, rad_clds, ring_ocs, &
       rad_clrs, ring_clrs
  ! radiance or ring effect fill in for clear (clrs), cloudy (cld), ocean (oc)
  real (KIND=8), dimension(:), allocatable :: w1p, f1p
  ! initail wavelength and flux vectors from dataset
  real (KIND=8), dimension(:), allocatable :: o3_xsect !Ozone cross-section
  real (KIND=8), dimension(:,:,:,:), pointer :: table !temporary 4D array
  real (KIND=8), dimension(:,:,:), pointer :: temp3D !temporary 3D array
  real (KIND=8), dimension(:,:), pointer :: temp2D !temporary 2D array
  real (KIND=8) :: sz, satz, az !solar zenith, satelite zenith, and 
  !satellite azimuth angles

  real (KIND=8) :: wavetol=0.25d0 !0.22! nm tolerance for marking observations as bad due to bad qc flag
  real (KIND=8) :: var_inv_chl=1./0.5d0**2 !1./5.0**2 ! inverse variance for 
  !b_i, used if retrieving chlorophyll
  real (KIND=8) :: var_inv_big=1./1e8**2 !inverse variance for b_i
  real (KIND=8) :: var_inv_small=(1.d0/1e-15)**2 !inverse variance for b_i
  real (KIND=8) :: wave_min = 300.d0 !minimum acceptable wavelength
  real (KIND=8) :: cld_frac, cld_frac_oc !cloud fraction, ocean cloud fraction
  real (KIND=8) :: diff_chi_max=0.03d0 !threshold for judging convergence,
  ! iteration stops if diff_chi is below this value.
  real (KIND=8) :: diff_chi_iter_max=-0.05d0 !threshold to check if chi^2
  ! has increased in current iteration
  real (KIND=8) :: i_obs_s, i_obs_l !normalized flux for longest (l) or 
  ! shortest (s) good wavelength
  logical, dimension(:), allocatable :: computed, comp_clr, comp_all 
  !  computed = vector tracking whether ring and rad have been computed
  !  comp_clr = vector tracking whether cloud fraction has been comp'd?
  !  comp_all = switch governing computation of radiance params in 
  !             m_interp_ring_rads
  logical, dimension(:), allocatable :: computed_oc, comp_oc_clr, comp_all_ring
  !  computed_oc = swicth governing whether to calculate ocean Ring effect
  !  computed = vector tracking whether ring and rad have been comp'd in 
  !  case where ocean effect is included. 
  !  comp_all_ring = switch governing computation of radiance params 
  !             in m_interp_ring_rads when ring effect is calculated.
  logical :: new_h, new_hcl, new_r
  !  switches indicating whether to update main jacobian (new_h) or
  !  chlorophyll jacobian (new_hcl). new_r indicates r_i has updated, 
  !  which can trigger the need to update the jacobians
  logical :: check_solar=.true. !check for bad solar irradiances
  logical :: check_rad=.true.  !check for bad radiances
  logical :: ret_chl  ! retrieve chlorophyll?
  logical :: add_oc   ! include ocean Ring effect?
  logical :: comp_clear  !only compute cloud-free pixels (?)
  logical :: bad_pix  !switch used in counting excluded pixels
  real (KIND=8) :: chloro  !chlorophyll concentration in single pixel
  real (KIND=8) :: adj  !spectral squeeze adjustment factor
  logical :: add_background_error = .true.  !use bg error covar. matrix?
  integer :: nsh  !number of state vector elements arising from spectral shift
  real (KIND=8)   :: reflec_cld  !cloud reflectivity
  real (KIND=8)   :: psurf, reflec  !surface pressure, reflectivity
  real (KIND=8), parameter :: refl_ice = 0.7d0  !ice reflectivity
  integer :: niters  !maximum number of iterations in chi^2 fit
  integer :: counts  !used in checking for bad radiances & irradiances
  integer, dimension(:), allocatable  :: indw  !wavelength index
  logical :: set_cld_frac  ! calculate cloud fraction?
  real (KIND=8) :: rad_tot_oc  ! total radiance of ocean scene
  integer :: i_r  ! index of good wavelength grid closest to wavelength
  !  of interest
  logical :: do_o2_jacob=.false. ! subtract a gradient across waveband?
  integer :: nchlr  !number of state vector elements arising from
  !  chlorophyll retrieval
  integer :: ind1  ! index of good wavelength grid closest to fill-in
  !  wavelength 
  real (KIND=8) :: elastic  !elastic scattering component of fill-in factor(?)

  ! Radiance params in Ring effect calculations
    !   i0_*: backscattered intensity at top of atmosphere?
    !   sb_*: fraction of intensity reflected by surface that is then
    !         scattered back to surface by atmosphere?
    !   tr_*: fractional transmittance of atmosphere?
  real (KIND=8), dimension(:,:), allocatable :: i0a_clds, tra_clds, &
       nra_clds, nia_clds, i01a_clds
  real (KIND=8), dimension(:,:), allocatable :: z1_clds, z2_clds
  real (KIND=8), dimension(:,:), allocatable :: ntot
  real (KIND=8) :: i0_l1, i0_l2, sb_l1, sb_l2, tr_l1, tr_l2
  real (KIND=8) :: i0_s1, i0_s2, sb_s1, sb_s2, tr_s1, tr_s2
  real (KIND=8) :: i0_ls, sb_ls, tr_ls
  real (KIND=8) :: i0_ss, sb_ss, tr_ss
  real (KIND=8) :: I0_l, I0_s, sb_l, sb_s, tr_l, tr_s
  integer :: ind0, indt  !indices
end module m_cloud_pres_mod

