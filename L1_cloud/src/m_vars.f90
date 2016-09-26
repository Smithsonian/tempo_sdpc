module m_vars 

  implicit none
  !-------------------------------------------------------------------------
  !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
  !-------------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE:  m_mvars
  ! 
  ! !DESCRIPTION: m_vars contains all the variables global to all routines used
  !
  ! !CALLING SEQUENCE: 
  !
  !        
  !     
  ! !INPUT PARAMETERS:   
  !
  ! !OUTPUT PARAMETERS:  
  !
  ! !SEE ALSO:  
  !
  ! !REVISION HISTORY: 
  !
  !  09Jan01   Murrin      original fortran 90
  !  09Aug02   Vasilkov    to exclude the hard wired filenames
  !  02Jan03   Vasilkov    updates for V1.2 L1B Reader
  !  19Aug14   O'Sullivan  added variable descriptions, removed those
  !                        which were unused.
  !
  !EOP
  !-------------------------------------------------------------------

  !LOGICAL CONTROL SWITCHES
  logical :: cal_reflec=.true. ! calculate dIdR (sensitivity of radiance 
  ! to reflectance)?
  logical :: get_cloud_frac=.true. ! calculate cloud fraction?
  logical :: write_resid=.false. ! write out residuals in output file?
  ! Note that param is set in PCF file
  logical :: write_obs=.false. ! write out observation fluxes instead of 
  ! residuals. Note: Set in PCF file, requires write_resid=true as well
  logical :: no_ret_ps=.false. ! do not retrieve pressure, use surface
  ! pressure instead. Note: set in PCF file
  logical :: noret=.false. ! do not perform retrievals, just test 
  ! input stages of code. Note: set in PCF file
  logical :: wrt_solar=.false. ! write out solar data and quit
  logical :: do_chl=.true. ! do chlorophyll correction only based on climatology
  logical :: cloud_fr_corr=.true. ! do cloud fraction corrections?
  logical :: done_read_terr=.false. ! has terrain pressure dataset been read?
  logical :: done_read_refl=.false. ! has TOMS reflectance been read?
  logical :: done_read_chl=.false. ! has chlorophyll clim. been read in?
  logical :: do_alloc=.true. ! call m_alloc1 in m_cloud_pres_ret
  logical :: do_alloc2=.true. ! call m_alloc2 in m_read_input_data
  logical :: squeeze=.false. ! include spectral squeeze in fit?
  logical :: shift=.true. ! include wavelength shift in fit?
  logical :: wr_shift=.true. ! write out fitted wavelength shift?
  logical :: do_short_wave=.true. ! include short wavlength bound in calc?
  logical :: ex= .true. ! read from L1B, write to L2 hdf files
  logical :: write_geom= .true. ! write out viewing geometry parameters
  logical :: write_ps =.true. ! write out surface (terrain) pressure
  logical :: using_resid=.true. ! read in residual (soft cal) data?
  logical :: using_spline=.false. ! use spline rather than linear 
  ! interpolation in m_cloud_pres_ret?
  logical :: transient_check=.false. ! check for transients in input 
  ! data (temporary hot pixels)? Set in PCF file.
  logical :: set_wmin=.false. ! has wmin been set?
  logical :: set_wmax=.false. ! has wmax been set?
  logical :: do_o3=.false. ! include Ozone correction? set by PCF
  logical :: write_fill=.true. ! write out filling-in factor?
  logical :: do_zoom=.true. !include measurements in spectral zoom mode?
  logical :: get_refl_clim=.true. ! read & use TOMS reflectance clim?
  logical :: do_LER=.false. ! do Lambertian Equiv. Reflector version of calc.
  logical :: do_mler=.true. ! do Mixed LER version of calculation
  logical :: test_solar=.false. ! use solar spectrum as input data?
  logical :: read_he4=.true. ! read/write he4 input & he5 output files?
  logical :: read_nc=.true. ! read/write netCDF input/output files?
  logical :: do_cloud_mask=.true. ! create cloud mask product?

  ! FILENAMES and PATHS
  character(len=255) :: solar_path='../out/' ! path used by write_solar
  character(len=255) :: sfile='first_omi_solar.dat' !filename for write_solar
  character(len=255) :: input_data_path='' ! used by read_input_data, 
  ! read_thresholds, cloud_mask
  character(len=255) :: output_data_path='' ! used by OMCLDRR
  character(len=255) :: out_path='' ! used by read_tables
  character(len=255) :: filename='' ! widely used placeholder
  character(len=255) :: thresh_file='' ! threshold reference file name
  character(len=255) :: filename_out='' ! output .he5 filename
  character(len=255) :: outswathname='Cloud Product' ! in output file
  character(len=255) :: ring_file_pre='ring_tab_omi_p' !Ring effect table
  character(len=255) :: ring_file_suf='.dat' !Rifng effect table suffix
  character(len=255) :: nc_swathname='band_540_740_nm' !netCDF input swathname

  ! FLAGS
  logical, allocatable, dimension (:) :: land_flg ! pixel on land?
  INTEGER (KIND=2), DIMENSION(:,:), allocatable :: geoflg ! geolocation flags
  INTEGER (KIND=1), DIMENSION(:,:), allocatable :: anomflg ! cross-track
  ! quality flags
  INTEGER (KIND=2), dimension(:), allocatable :: mflg ! measurement 
  ! quality flag from L1B data
  INTEGER (KIND=2), dimension(:), allocatable :: meas_qual_flg ! measurement
  ! quality flags determined internally
  integer (kind=2), allocatable, dimension (:,:) :: qc, qc2 ! quality flags
  INTEGER (KIND=2), DIMENSION(:,:), allocatable :: quality_flagL, irr_quality_flagL
  ! radiance and irradiance quality flags
  integer, parameter :: min_refl_flag=6 ! flag for violations of min_refl
  integer, parameter :: bad_obs_flag=7 ! flag for spectr in which more
  ! than 50% of fluxes are bad


  ! M_READ_INPUT_DATA related
  integer :: start_line=1 ! for read_input_data
  integer :: max_lines=0 ! for read_input_data


  ! CLOUD MASK related
  real (KIND=8), allocatable, dimension (:,:) :: thresholds ! thresholds on
  ! standard deviation for use in cloud mask code
  integer, pointer, dimension (:) :: npixels !valid numbers of small pixels
  ! for which thresholds exist, used in cloud mask code
  real (KIND=8) :: stddev_thresh=0.001d0 ! initial standard deviation
  ! threshold for cloud_maks, value based on OMI small pix.
  integer (kind=2), allocatable, dimension (:,:) :: cloud_mask ! mask value 


  ! SOLAR related
  real (KIND=8), allocatable, dimension (:,:) :: ws ! solar wavelengths (nm)
  real (KIND=8), allocatable, dimension (:,:) :: fs ! solar irradiances
  real (kind=4), allocatable, dimension (:,:) :: sazimuth, vazimuth ! solar
  ! and viewing azimuth angles
  integer (kind=1) :: config_irr ! instrument config for solar obs.


  ! TOMS reflectance climatology related
  integer :: ref_nlat=180 !latitude range for TOMS reflectances
  integer :: ref_nlon=360 !longitude range for TOMS reflectances
  integer :: ref_nmon=12 !month range for TOMS reflectances
  integer,parameter :: ler_nsz=16, ler_nth=16, ler_nph=16 ! number of theta,
  ! phi, solar zenith points in TOMS reflectance climatology
  real (KIND=8), dimension(ler_nsz) :: ler_sz ! TOMS solar zenith values
  real (KIND=8), dimension(ler_nth) :: ler_th ! TOMS theta values
  real (KIND=8), dimension(ler_nph) :: ler_ph ! TOMS phi values
  real (KIND=8), dimension(ler_nsz,ler_nth,ler_nph) :: ler354 ! TOMS 
  ! reflectance climatology table
  real (kind=4), dimension(:), allocatable :: ref_lats, ref_lons !latitude
  ! and longitude values for TOMS reflectance climatology
  REAL (KIND=8) :: startlat,startlon,deltlat,deltlon ! starting latitude
  ! and longitude, fractional angles for matching to observed grid.
  real (kind=4), dimension(:,:,:), allocatable :: toms_refl ! TOMS reflectances


  ! CHLOROPHYLL related
  integer, parameter :: nlat=361 ! latitude range for chlorophyll clim.
  integer, parameter :: nlon=720 ! longitude range for chlorophyll clim.
  real (kind=4), dimension(nlon,nlat) :: chl2d ! chlorophyll climatology
  real (kind=4), allocatable, dimension(:,:) :: chlorophyll
  !chlorophyll concentration


  ! OCEAN table related
  integer :: nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl ! for ocean
  ! tables, numbers of wavelengths, theta, scan, phi, reflectances, 
  ! chlorophyll concentrations.
  real (kind=8), dimension(:), allocatable :: theta_oc, scan_oc, phi_oc
  ! ocean theta, scan, phi values
  real (kind=8), dimension(:), allocatable :: ocrefl, chl
  ! ocean reflectance, chlorophyll concentration values
  real (KIND=8),   dimension(:,:,:,:), pointer :: oc_table ! ocean data
  ! table (thet_oc,scan_oc,nchl,nwave2)
  real (kind=4), dimension(:,:,:,:,:,:), allocatable :: oc_perms ! ocean
  ! data table (oc_perms,wgrid_out_oc,theta_oc,scan_oc,phi_oc,ocrefl,chl)
  real (KIND=8), dimension(:), allocatable :: wgrid_oc !ocean wavelength grid
  real (kind=4), dimension(:), allocatable :: wgrid_out_oc !wavlength grid
  ! used when writing out ocean values


  ! TERRAIN heigh/pressure related
  integer (kind=2), pointer, dimension(:,:) :: terr_height ! terrain height
  integer, parameter :: l1b_nlat=360 !latitude range in deg, used by pterr
  integer, parameter :: l1b_nlon=720 !longitude range in deg, used by pterr
  real (kind=4), dimension(l1b_nlat,l1b_nlon) :: p_terr ! terrain pressure



  !integer :: iprt=0  ! verbosity level
  integer :: err_code ! error code in OMCLDRR
  integer :: retstatus ! status of metadata write-out
  character (len=200) :: flnm_out ! L2 output filename
  integer (kind=1), dimension(:), allocatable :: meas_class, config_rad
  ! measurement class (Earth, Dark, Sun, etc), instrument config.
  integer :: ispec ! cross-track iteration counter in OMCLDRR output
  integer :: nsolwave ! number of wavelength bins in solar spectrum
  integer :: nwave2 ! nuumber of wavelengths, used in read_ocean_table
  integer :: niter= 10 !6 ! max number of iterations in chi^2 fit
  INTEGER (KIND=4) :: version, status, ierr ! version is used in L1B access
  ! status = result of PGS operations, ierr = result of SMF operations
  INTEGER (KIND=4) :: nTimes=-1, nXtrack, nWavel, nWavelCoef
  ! nTimes = number of pix along swath,    nXtranck = no. of pix across track
  ! nWavel = no. of wavelengths in spectrum, nWavelCoef = parameter used 
  ! in calculation of wavelength scale
  INTEGER (KIND=4) :: iLine=0, nwl, nLines ! nwl = no. of valid wavelengths
  ! nLines = number of pixels along swath, iLine = index for nLines counter
  INTEGER (KIND=4), PARAMETER :: zero=0 ! the additive identity
  integer :: npres=5 ! number of pressures
  integer :: ntheta, nscan, nphi ! number of thetas, scans, phis

  real (KIND=8), allocatable, dimension (:,:) :: resid_spec !residual
  ! spectrum used for "soft calibration"

  real (KIND=4), dimension (:),     pointer :: wave_resid ! wavelength
  ! grid for output residual spectrum
  real (KIND=4), dimension (:,:,:), pointer :: resid ! output residuals
  real (KIND=4), dimension (:,:),   allocatable :: w12d
  real (KIND=4), dimension (:,:),   allocatable :: f12d

  real (KIND=8), pointer, dimension (:) :: k1bar
  real (KIND=8), pointer, dimension (:,:) :: sba, nba
  real (KIND=8), pointer, dimension (:,:,:,:,:) :: i01a, nia
  real (KIND=8), pointer, dimension (:,:,:,:) :: tra, nra, i0a, z1, z2
  real (kind=4), allocatable, dimension (:,:) :: refl, azimuth, sza, sat_zen
  ! reflectance, satellite azimuth, solar zenith angle, satellite zenith
  real (kind=4), allocatable, dimension (:,:) :: reflect_cld ! Cloud reflectivity
  real (KIND=8), allocatable, dimension (:) :: w_grid, chlcl ! wavelength
  ! grid, chlorophyll concentration
  real (kind=4), allocatable, dimension (:,:) :: ps ! output terrain pressure 
  real (kind=4), allocatable, dimension (:,:) :: shifts, shifts2, squeezes
  ! wavelength shift and squeeze values
  real (kind=4), allocatable, dimension (:,:) :: ref_clr ! reflectance
  ! with clear sky 
  real (kind=4), allocatable, dimension (:,:) :: dIdR ! sensitivity of
  ! radiance to reflectance
  real (kind=4), allocatable, dimension (:,:) :: fill ! filling-in factor
  real (kind=4), allocatable, dimension (:,:) :: lat, lon, & 
       eff_cld_frac, eff_cld_frac2, rad_cld_frac
  ! latitude, longitude, effective and radiative cloud fractions
  real (kind=4), allocatable, dimension (:,:) :: chi_sqr, chi_sqr2
  ! chi squared in each pixel
  real (kind=4), allocatable, dimension (:,:) :: biases, biases2
  ! residual bias in each pixel
  real (kind=4), allocatable, dimension (:,:) :: stds, stds2
  ! residual standard deviation in each pixel, 
  integer (kind=2), pointer, dimension (:,:) :: smpx_nPix ! number of
  ! small pixels in each full size pixel
  real (kind=4), pointer, dimension (:,:) :: smpx_mean, smpx_stddev
  ! mean and standard deviation of small pixel values in each full pixel
  real (kind=4), pointer, dimension (:,:) :: smpx_wavel ! mean of small
  ! pixel wavelengths in each ful pixel
  real (kind=8), dimension(:), pointer :: time ! time at start of scan
  real (KIND=8) , dimension(:), allocatable :: theta, scan, phi
  ! theta = solar zenith angle, scan = satellite zenith angle
  ! phi = azimuth angle      - all from Ring effect table
  real (KIND=8) , dimension(:), allocatable :: wgrid_out, sflx, pres
  ! sflx = solar flux, pres = pressure       - all from Ring effect table
  ! wgrid_out = output wavelength grid
  integer :: nwave, nwav ! number of wavelengths in use

  integer :: npixs, nscanpos ! npixels, n scan positions, see read_thresholds

  real (kind=4), allocatable, dimension(:,:) :: cloud_pres, cld_pres2
  !cloud pressure for O3

  real (KIND=8) :: refl_clr=0.15d0!0.11d0 !clear sky reflectance
  real (KIND=8) :: refl_cld=0.80d0!0.40d0 !cloud reflectance
  ! Note that refl_clr and refl_cld were 0.11 and 0.4 in pipeline version
  ! of m_vars, but were overwritten to the current values by OMCLDRR
  real (KIND=8) :: refl_l, refl_s ! reflectances at long and short wavel.
  real (KIND=8) :: refl_clr_oc=0.10d0 !0.08 !ocean reflectance
  real (KIND=8) :: cld_frac_min=0.05d0 !threshold for clear scene (?)
  real (KIND=8), parameter :: max_refl=1.00d0 ! max allowed reflectance
  real (KIND=8), parameter :: min_refl=0.00d0 ! min allowed reflectance

  real (KIND=4) :: wmin= 391 ! min wavelength in nm. Set by PCF
  real (KIND=4) :: wmax= 398.5 ! max wavelength in nm. Set by PCF
  real (KIND=4) :: wmin2= 350 ! min avelength in nm
  real (KIND=4) :: wmax2= 405 ! max wavelength in nm
  real (KIND=8) :: wdelt=0.5d0 ! accuracy to match wavelengths in nm
  real (KIND=8) :: wave_short= 376.4d0 !346.8!340.4
  real (KIND=8) :: wave_long=394.1 !376.4 ! 390 where Ring close to 0, 


  real (KIND=4), parameter :: fill_value=-9999.0 ! fill value for null data
  integer (KIND=2), parameter :: fill_value_int=-1 ! fill value for nulls
  integer (KIND=1), parameter :: fill_value_int1=-1 ! fill value for nulls

  integer :: n_good_input=0, n_good_output=0 ! no of good input and output
  ! values, used to get good percentage for dataset
  integer :: n_input=0, n_missing=0 ! no of cross-track data points, 
  ! number missing after qualityfalgging
  integer :: highqual=80, badqual=20 ! used in MetadataModule to class data
  ! as "passed" (>80% good), "suspect (20-80% good), "failed" (<20% good)
  integer, parameter :: min_wl=10 ! minimum number of valid wavelengths to 
  ! continue using spectrum.
  real (KIND=8) :: cld_frac=0.02d0 ! cloud_fraction

  real(kind=4) :: dist_rad, dist_irrad ! radiance and irradiance distances
  ! from m_earth_sun_dist

  real (KIND=8)    :: wave_fill=352.6d0 ! nm wavelength at which to output 
  ! filling-in factor
  real (KIND=8), allocatable, dimension(:) :: wave_o3, xsect_o3 ! Ozone
  ! cross-section table wavelength and relative flux 
  integer :: year, month, day ! observation date 

  real (kind=4) :: add_shift=0.0 ! Optional wavelength shift, to be added
                                 ! to solar spectrum if used as input

  !A(I,J)= corresponding albedo spectrum J at wavelength index I
  !          ie, A = (I/F)
  !=====================================================================
  !real (KIND=8), allocatable, dimension (:,:) :: A


end module m_vars
