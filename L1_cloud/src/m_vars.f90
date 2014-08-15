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
  !  09Jan01   Murrin     original fortran 90
  !  09Aug02   Vasilkov   to exclude the hard wired filenames
  !  02Jan03   Vasilkov   updates for V1.2 L1B Reader
  !
  !EOP
  !-------------------------------------------------------------------
  integer :: retstatus
  character (len=200) :: flnm_out

  real (KIND=8), allocatable, dimension (:,:) :: reference_spec
  real (KIND=8), allocatable, dimension (:,:) :: reference_wave
  real (KIND=8), allocatable, dimension (:,:) :: resid_spec
  real (KIND=8), allocatable, dimension (:)   :: resid_wave
  real (KIND=8), allocatable, dimension (:)   :: cal_fact
  real (KIND=8), allocatable, dimension (:,:) :: reference_ring
  real (KIND=8), allocatable, dimension (:,:) :: reference_rad
  real (KIND=8), allocatable, dimension (:,:) :: ws ! solar wavelengths (nm)
  real (KIND=8), allocatable, dimension (:,:) :: fs ! solar irradiances

  real (KIND = 4), dimension (:),     pointer :: wave_resid
  real (KIND = 4), dimension (:,:,:), pointer :: resid
  real (KIND = 4), dimension (:,:,:), allocatable :: w1
  real (KIND = 4), dimension (:,:,:), allocatable :: f1
  real (KIND = 4), dimension (:,:),   allocatable :: w12d
  real (KIND = 4), dimension (:,:),   allocatable :: f12d

  logical :: swap_geo = .false. !.true.
  real (KIND=8), allocatable, dimension (:,:) :: thresholds
  integer, pointer, dimension (:) :: npixels
  real (KIND=8), pointer, dimension (:) :: k1bar
  real (KIND=8), pointer, dimension (:,:) :: sba, nba
  real (KIND=8), pointer, dimension (:,:,:,:,:) :: i01a, nia
  real (KIND=8), pointer, dimension (:,:,:,:) :: tra, nra, i0a, z1, z2
  real (kind = 4), allocatable, dimension (:,:) :: refl, azimuth, sza, sat_zen
  real (kind = 4), allocatable, dimension (:,:) :: ai, reflect_cld
  real (KIND=8), allocatable, dimension (:) :: w_grid, chlcl
  real (kind = 4), allocatable, dimension (:,:) :: ps
  real (kind = 4), allocatable, dimension (:,:) :: shifts, shifts2, squeezes
  real (kind = 4), allocatable, dimension (:,:) :: ref_clr
  real (kind = 4), allocatable, dimension (:,:) :: dIdR
  real (kind = 4), allocatable, dimension (:,:) :: fill
  integer (kind = 2), pointer,dimension(:,:) :: terr_height
  integer (kind = 2), allocatable, dimension (:,:) :: qc, qc2 ! quality control flag
  !!integer (kind = 2), pointer, dimension (:,:) :: cloud_mask 
  integer (kind = 2), allocatable, dimension (:,:) :: cloud_mask 
  real (kind = 4), allocatable, dimension (:,:) :: sazimuth, vazimuth 
  real (kind = 4), allocatable, dimension (:,:) :: lat, lon, &
       eff_cld_frac, eff_cld_frac2, rad_cld_frac
  real (kind = 4), allocatable, dimension (:,:) :: chi_sqr, chi_sqr2, biases, biases2, stds, stds2, rms
  integer (kind = 2), pointer, dimension (:,:) :: smpx_nPix
  real (kind = 4), pointer, dimension (:,:) :: smpx_mean, smpx_stddev
  real (kind = 4), pointer, dimension (:,:) :: smpx_wavel
  real (KIND=8), pointer, dimension (:,:) ::  info
  logical, allocatable, dimension (:) :: land_flg
  INTEGER (KIND = 2), DIMENSION(:,:), allocatable :: geoflg
  INTEGER (KIND = 1), DIMENSION(:,:), allocatable :: anomflg
  INTEGER (KIND = 2), dimension(:), allocatable :: mflg
  INTEGER (KIND = 2), dimension(:), allocatable :: meas_qual_flg
  integer (kind = 1), dimension(:), allocatable :: meas_class, config_rad
  integer (kind = 1) :: config_irr
  real (kind = 8), dimension(:), pointer :: time
  REAL (KIND = 4), DIMENSION(:,:), allocatable :: rad_precisionL
  INTEGER (KIND = 2), DIMENSION(:,:), allocatable :: quality_flagL, irr_quality_flagL
  real (KIND=8) , dimension(:), allocatable :: theta, scan, phi, wgrid_out, sflx, pres
  real (KIND=8) , dimension(:), allocatable :: scan_cal, wgrid_cal, time_cal
  real (KIND=8) , dimension(:), allocatable :: refl_ref, psurf_ref, sza_ref, &
       satz_ref, az_ref, lat_ref, lon_ref
  real (KIND=8) , dimension(:,:,:), allocatable :: cal_table
  real (kind = 4), dimension(:), allocatable :: wgrid_out2
  integer :: nwave, nwav

  real (KIND=8), dimension(:), allocatable, target :: w0
  real (KIND=8), dimension(:), allocatable, target :: f0

  integer :: npixs, nscanpos
  integer :: n_sol_spec=1
  integer, parameter :: lun_out_resid=5200008
  integer, parameter :: lun_out_cal=5200007
  character(len=2) :: nstr=''
  character(len=255) :: solar_path='../out/'
  character(len=255) :: sfile='first_omi_solar.dat'
  character(len=255) :: input_data_path=''
  character(len=255) :: output_data_path=''
  character(len=255) :: out_path=' /omi/testdata/new/OMCTPo/coef/'
  character(len=255) :: out_solar_filename=''
  !! I don't understand why these are set to 200 when all others are 255,
  !! and as 255-length names can get written into 200-length, altering
  !character(len=200) :: filename=''
  !character(len=200) :: filename_cm=''
  !character(len=200) :: filename_gome=''
  character(len=255) :: filename=''
  character(len=255) :: filename_cm=''
  character(len=255) :: filename_gome=''
  character(len=255) :: filename_solar=''
  character(len=255) :: thresh_file='cld_mask_thresh.txt'
  character(len=255) :: filename_out='cld_test_output.hdf'
  character(len=255) :: filename_resid='cld_test_resid.txt'
  character(len=255) :: filename_cal='cld_test_cal.txt'
  character(len=255) :: filename_no2='no2.asc'
  character(len=255) :: outswathname='Cloud Product'
  character(len=255) :: resource_file='cloudret.rc'
  integer :: nfiles=0
  integer, parameter :: maxfiles=255
  character(len=255), dimension(maxfiles) :: filenames

  !A(I,J)= corresponding albedo spectrum J at wavelength index I
  !          ie, A = (I/F)
  !=====================================================================
  !real (KIND=8), allocatable, dimension (:,:) :: A

  real (kind=4), allocatable, dimension(:,:) :: cloud_pres, cld_pres2
  real (kind=4), allocatable, dimension(:,:) :: eta
  real (kind = 4), allocatable, dimension(:,:) :: chlorophyll

  integer :: ispec
  real (KIND=8)    :: dx
  integer ::  nsolwave
  integer :: err_code

  real (KIND=8), dimension(:,:), pointer :: sol2, sol3
  real (KIND=8), dimension(:,:,:), pointer :: rad2,  rad3
  real (KIND=8), pointer, dimension(:,:) :: s2, s3
  real (kind = 4), allocatable, dimension(:,:) :: s2b, s3b
  real (kind = 4), pointer, dimension(:,:,:) :: rad2b, rad3b
  real (KIND=8) :: wbeg, wend
  integer :: iprt=0
  integer :: date
  character(len=255) :: ring_file_pre='ring_tab_omi_p'
  character(len=255) :: ocean_file_pre
  character(len=255) :: ring_file_suf='.dat'
  character(len=255) :: ocean_file_suf
  logical :: cal_reflec=.true.
  logical :: get_cloud_frac=.true.
  logical :: write_resid=.false.
  logical :: write_resid_all =.false. !.true. !.false.
  !logical :: write_resid=.false.
  logical :: write_obs=.false.
  logical :: opened_write_resid=.false.
  logical :: no_ret_ps=.false.
  !logical :: no_ret_ps=.true.
  integer :: niter= 10 !6

  real (KIND=8) :: cal_const
  real (KIND=8) :: refl_chl_max = 0.40d0
  real (KIND=8) :: refl_clr_oc = 0.10d0 !0.08
  real (KIND=8) :: refl_clr = 0.11d0 !0.08 !0.15
  real (KIND=8) :: refl_cld = 0.40d0 !0.80
  real (KIND=8) :: refl_clr2 = 0.11d0 ! 0.15
  real (KIND=8) :: refl_cld2 = 0.80d0
  real (KIND=8) :: cld_frac_min = 0.05d0 ! 0.20 !0.15
  real (KIND=8) :: cld_frac_max = 1.0d0
  real (KIND=8), parameter :: max_refl = 1.00d0
  real (KIND=8), parameter :: min_refl = 0.00d0
  real (KIND=8), parameter :: max_ai   = 1.00d0
  integer, parameter :: min_refl_flag = 6
  integer, parameter :: bad_obs_flag = 7
  integer, parameter :: ai_flag = 8
  INTEGER (KIND = 4) :: version, status, ierr, &
       nTimes=-1, nXtrack, nWavel, nWavelCoef, &
       iLine=0, nwl, nLines
  INTEGER (KIND = 4), PARAMETER :: zero = 0
  logical :: noret = .false.
  logical :: wrt_solar = .false.
  integer :: ifile
  logical :: no_cl_filename = .false.
  logical :: wave_adjust = .true.
  integer :: interp = 1
  integer :: nwave2 
  integer :: nwave_cal, ntime_cal, nscan_cal
  real (KIND = 4) :: wmin= 391 !365.! 356.
  real (KIND = 4) :: wmax= 398.5
  real (KIND = 4) :: wmin2= 350 !365.! 356.
  real (KIND = 4) :: wmax2=405 !398.1
  real (KIND=8)            :: wdelt=0.5d0

  integer :: npres=5
  integer :: ntheta, nscan, nphi
  integer :: nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl
  real (KIND=8),   dimension(:,:,:,:), pointer :: oc_table
  real (kind = 4), dimension(:,:,:,:,:,:), allocatable :: oc_perms
  real (kind = 4), dimension(:), allocatable :: wgrid_out_oc
  real (KIND=8),   dimension(:), allocatable :: wgrid_oc
  real (kind = 8), dimension(:), allocatable :: theta_oc,scan_oc,phi_oc,ocrefl,chl
  logical :: read_resource_file
  !logical :: retrieve_chl_pres = .true. ! retrieve chl and pres when clear
  logical :: retrieve_chl_pres = .false. ! retrieve chl and pres when clear
  logical :: retrieve_chl_clr = .false. ! retrieve only chl when clear (not cld pres)
  !logical :: retrieve_chl_clr = .true.  ! retrieve only chl when clear (not cld pres)
  logical :: retrieve_chl_cld = .false. ! proxy for ret. cld liq water when cloudy
  logical :: do_chl = .true. ! do chlorophyll correction only based on climatology
  !logical :: do_chl = .false.
  logical :: cloud_fr_corr = .true.
  !logical :: squeeze = .true.
  logical :: squeeze = .false.
  integer, parameter :: l1b_nlat=360
  integer, parameter :: l1b_nlon=720
  real (kind = 4), dimension(l1b_nlat,l1b_nlon) :: p_terr
  integer :: ref_nlat=180
  integer :: ref_nlon=360
  integer :: ref_nmon=12
  integer,parameter :: ler_nsz=16, ler_nth=16, ler_nph=16
  real (KIND=8), dimension(ler_nsz) :: ler_sz
  real (KIND=8), dimension(ler_nth) :: ler_th
  real (KIND=8), dimension(ler_nph) :: ler_ph
  real (KIND=8), dimension(ler_nsz,ler_nth,ler_nph) :: ler354
  logical :: done_read_terr = .false.
  logical :: done_read_refl = .false.
  integer, parameter :: nlat=361
  integer, parameter :: nlon=720
  integer            :: start_line=1
  integer            :: max_lines=0
  real (kind = 4), dimension(nlon,nlat) :: chl2d
  logical :: done_read_chl = .false.
!  !logical :: simul=.true. 
!  logical :: simul=.false.
  logical :: do_alloc=.true., do_alloc2=.true.
  logical :: write_he5=.true.
  logical :: do_no2 =.false.
  !logical :: shift=.false. 
  logical :: wr_shift=.true.
  logical :: shift=.true.
  logical :: do_short_wave=.true.
  real (KIND=8), dimension(:,:), allocatable :: no2
  real (KIND=8) :: refl_l, refl_s

  real (KIND=8)    :: wave_short= 376.4d0 !346.8!340.4
  !real (KIND=8)    :: wave_long = 395.5d0 !376.4 ! 390 where Ring close to 0, 
  !real (KIND=8)    :: wave_long = 386.3d0 !376.4 ! 390 where Ring close to 0, 
  real (KIND=8)    :: wave_long = 394.1 !376.4 ! 390 where Ring close to 0, 
  ! old 386.3 !373.2 ! for GOMI 386.3! 
  !logical :: ex= .true
  logical :: ex= .true.
  logical :: write_geom= .true.
  logical :: write_ps  = .true. !.false.
  logical :: using_ref = .false. !.true.
  logical :: using_resid = .true. !.true.
  logical :: using_cal = .true.
  !logical :: using_resid = .false. !.true.
  logical :: using_spline = .false.
  !logical :: using_spline = .true.
  logical :: gomi=.false.
  !real (KIND=8)    :: stddev_thresh=0.01d0 ! initial value based on GOME PMD
  real (KIND=8)    :: stddev_thresh=0.001d0 ! initial value based on OMI small pix.
  real (KIND=4), parameter :: fill_value = -9999.0
!!!Altered since m_write_ouput_data seems to expect these to be int*2 and int*1
  !integer, parameter :: fill_value_int = 65535
  !integer, parameter :: fill_value_int1 = 255
  integer (KIND=2), parameter :: fill_value_int = -1
  integer (KIND=1), parameter :: fill_value_int1 = -1
  integer :: n_good_input = 0, n_good_output = 0, n_input = 0, n_missing = 0
  integer :: highqual = 80, badqual = 20
  integer, parameter :: min_wl = 10
  real (KIND=8) :: cloud_pres_max = 1100.0, cld_frac = 0.02d0

  real(kind=4) :: dist_rad, dist_irrad 
  logical :: cloud_clear=.false.!.true.
  integer :: ny=1, nx=1
  integer :: ll
  integer :: nfov=2
  real (KIND=8)    :: wave_dpdf = 393.6d0 !nm
  logical :: transient_check = .false. !.true.
  logical :: set_wmin = .false. !.true.
  logical :: set_wmax = .false. !.true.
  integer :: n_products = 1! 2
  logical :: do_o3=.false.
  logical :: write_fill=.true.
  logical :: do_zoom=.true.
  logical :: get_refl_clim=.true.
  real (KIND=8)    :: wave_fill=352.6d0 ! nm wavelength to output filling in
  real (KIND=8), allocatable, dimension(:) :: wave_o3, xsect_o3
  integer :: year, month, day

  real (kind = 4), dimension(:,:,:), allocatable :: toms_refl
  real (kind = 4), dimension(:), allocatable :: ref_lats, ref_lons
  !logical :: do_LER=.true.
  logical :: do_LER=.false.
  logical :: do_mler=.true.
  !logical :: write_ai=.true.
  logical :: write_ai=.false.

  !used in m_rd_toms_refl, need values to persist over multiple
  !iterations of cloud_ret
  REAL (KIND=8) :: startlat,startlon,deltlat,deltlon

end module m_vars
