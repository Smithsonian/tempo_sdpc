!*************
module m_vars
!*************

!---------------------------------------------------------------------72
! ROUTINE: m_mvars
! 
! DESCRIPTION: m_vars contains all filenames and variables
!
! REVISION HISTORY: 
!
!  04/23/15 Yang original fortran 90
!  2021 Wang adaption to TEMPO
!  2023 Wang modification for TEMPO
!---------------------------------------------------------------------72
  implicit none

!------------
! misc
!-----------
! operation mode ='production' or 'development', set through control
   character(len=255):: run_mode='development'
   ! number of bytes used in processing_quality_flag
   ! integer, parameter:: pflag_nbyte=4 
   integer, parameter:: pflag_nbyte=2
!----------------
! L1B irradiance
!----------------
  character(len=255)::name_irr_dir='./'
  character(len=255)::name_irr_file='empty'
  character(len=255)::name_irr_swath='Sun Volume VIS Swath'
  real(kind=8)::irr_Time
  real(kind=4)::irr_SecondsInDay
  integer(kind=4)::irr_NumTimes
  integer(kind=4)::irr_nXtrack
  integer(kind=4)::irr_nWavel
  integer(kind=4)::irr_nWavelCoef
  real::irr_EarthSunDist 
  real(kind=4),dimension(:),pointer::irr_out_irradiance_440nm
  real(kind=4),dimension(:),pointer::irr_out_irradiance_466nm
  real(kind=4),dimension(:),pointer::irr_out_irradiance_477nm

!--------------
! L1B radiance
!--------------
  character(len=255)::name_rad_dir='./'
  character(len=255)::name_rad_file='empty'
  character(len=255)::name_rad_swath='Earth VIS Swath'
  real(kind=8),   dimension(:),    pointer::rad_Time
  real(kind=4),   dimension(:),    pointer::rad_SecondsInDay
  real(kind=4),   dimension(:,:),  pointer::rad_Longitude
  real(kind=4),   dimension(:,:),  pointer::rad_Latitude
  real(kind=4),   dimension(:,:),  pointer::rad_SolarZenithAngle
  real(kind=4),   dimension(:,:),  pointer::rad_SolarAzimuthAngle
  real(kind=4),   dimension(:,:),  pointer::rad_ViewingZenithAngle
  real(kind=4),   dimension(:,:),  pointer::rad_ViewingAzimuthAngle
  real(kind=4),   dimension(:,:),  pointer::rad_RelativeAzimuthAngle

  integer(kind=4)::rad_NumTimes
  integer(kind=4)::rad_nXtrack
  integer(kind=4)::rad_nWavel
  !integer(kind=4)::rad_nWavelCoef !not used
  real::rad_EarthSunDist

! read in and transfer to out_TerrainHeight
  integer(kind=2),dimension(:,:),  pointer::rad_TerrainHeight
!OMI uses GroundPixelQualityFlags to decide snow/ice
!TEMPO uses L1 snow_ice_fraction instead
!  integer(kind=4),dimension(:,:),  pointer::rad_GroundPixelQualityFlags
  real(kind=4),   dimension(:,:),   pointer::rad_SnowIceFraction

! added rad_466nm,rad_477nm,rad_440nm which is what needed
  real(kind=4), dimension(:,:), pointer::rad_466nm,rad_477nm,rad_440nm

!--------------
! Lookup table
!--------------
! these initial names can be changed through control file
  character(len=255)::name_lut_dir='./refdata/'
  character(len=255)::name_lut_rad440='LUT_4400_RAD.h5'
  character(len=255)::name_lut_rad='LUT_4660_RAD.h5'
  character(len=255)::name_lut_ler='LUT_4770_CLOUD_LER_6D.h5'
  character(len=255)::name_lut_amf_clr='LUT_4770_CLEAR.h5'
  character(len=255)::name_lut_amf_cld='LUT_4770_CLOUD_MLER.h5'
  real(kind=4),dimension(:,:,:,:,:,:),pointer::lut_amf_ler
  real(kind=4),dimension(:,:,:,:,:),pointer::lut_rad_ler
  real(kind=4),dimension(:,:,:,:,:),pointer::lut_rad_clr440
  real(kind=4),dimension(:,:,:,:,:),pointer::lut_rad_clr
  real(kind=4),dimension(:,:,:,:,:),pointer::lut_amf_clr
  real(kind=4),dimension(:,:,:,:,:),pointer::lut_amf_cld
  real(kind=4),dimension(:),pointer::lut_alb
  real(kind=4),dimension(:),pointer::lut_sza
  real(kind=4),dimension(:),pointer::lut_vza
  real(kind=4),dimension(:),pointer::lut_raa
  real(kind=4),dimension(:),pointer::lut_psfc
  real(kind=4),dimension(:),pointer::lut_pcld

  ! LUT cloud node
  ! ALB(cloud) = 0.8 (ind=18), Pcld = 700hPa (ind=18)
  integer, parameter:: LUTrad_cloud_albid = 18
  integer, parameter:: LUT_pcld_700hPa = 18
  ! during ecfocp iteration LUTrad_cloud_psfcid changes for each pix
  integer:: LUTrad_cloud_psfcid = 18 !initial

  ! 1-based lut_pcld index for ocp whose height is just above cloud
  ! i.e. pressure is just smaller than cloud pressure
  integer,dimension(:,:),allocatable:: lut_pcld_indarr
  ! max number of ecfocp iterations 
  integer(kind=4):: ecfocp_maxiter = 10
  ! results from previous iteration
  real(kind=4),dimension(:,:),allocatable:: prev_ecf,prev_ecf_notclipped
  real(kind=4),dimension(:,:),allocatable:: prev_ocp,prev_ocp_notclipped
  integer(kind=4),dimension(:,:),allocatable:: ecf_niter,ocp_niter

  ! ecfocp stopping threshold
  real(kind=4), parameter:: delta_ocp_abs = 1.0 ! hPa
  real(kind=4), parameter:: delta_ocp_pct = 0.01 ! relative percent change
  real(kind=4), parameter:: delta_ecf = 0.005

  ! LUT albedo node used for solving sbar, trans 
  ! ALB=0.0, 0.1, 0.2 used to calc tran & sbar in pscene
  integer, parameter:: LUT_ALBID_0p0 = 1
  integer, parameter:: LUT_ALBID_0p1 = 7
  integer, parameter:: LUT_ALBID_0p2 = 12

  !OMI LUT dimension, not used for TEMPO
  !integer,parameter::nalb=20, nsza=30, nvza=19, nraa=37
  !integer,parameter::npsfc=23, npcld=23, nrsfc=23
  !TEMPO LUT dimension
  integer,parameter:: nalb=20, nsza=30, nvza=25, nraa=37
  integer,parameter:: npsfc=23, npcld=23

  ! added the following to remove hardcoded 
  ! these are determined by LUTs
  !real(kind=4),parameter::max_SZA=89.,max_VZA=72. ! OMI
  real(kind=4),parameter:: max_SZA=89., max_VZA=89. ! TEMPO

!------------------------
! O4 SCD temperature correction coefficients 
! y(T2) = a * y(T1) + b; T1 = 223K, T2=263, 293K (Finkenzeller)
!------------------------
! coefs are derived using ops3_4p3_livetest 20240216 with
! new O2O2_template_feb2024.pcf & control.O2O2_feb2024.in  
   real, parameter:: TrefO4 = 223. 
   real, parameter:: a263 = 1.0043, b263 = 0.0066
   real, parameter:: a293 = 1.0219, b293 = 0.1117

   ! maximum number of iteration for SCD temperature adjustment
   integer, parameter :: max_scd_iter = 20 
   ! if dT < dt_threshold, then stop iteration
   real, parameter :: dt_threshold = 1.0 !K 

!-------------------------
! vertical column density
!-------------------------
  integer,parameter::nvcd=npcld
! vvcd will be replaced with actual gmi_vcd//geos_vcd 
  real,dimension(nvcd):: vvcd &
         =(/0.00472129,0.00648191,0.00889265,0.0121931,0.0166994, &
           0.0228845, 0.0313553, 0.0429509, 0.0588269,0.0805602, &
           0.109367,  0.146245,  0.193142,  0.252468, 0.326837,  &
           0.419688,  0.534917,  0.677226,  0.852079, 1.06600,  &
           1.32557,1.41490,1.53974/)

 ! added the multiplicative conversion factor for calculating O4 VCD
 !    this removes hardcoded constant in many routines
 ! EY suggests change to 6.733e-4 to be more accurate
  real,parameter::vcd_convfac = 6.733e-4 !previously 6.765e-4

 ! add the fraction used for cpp during scd temperature correction
  real,parameter:: frac4cpp = 0.7937 !previously 0.5

!-----------
! input LUN
!-----------
character(len=6)::lun_lut_amf_clear='477000'
character(len=6)::lun_lut_amf_cloud='477001'
character(len=6)::lun_lut_amf_ler6d='477010'
! month dependent GMI files should keep the lun_ stuff
character(len=6),dimension(12):: &
  lun_gmi_psfc=(/'400201','400202','400203','400204','400205','400206', &
                 '400207','400208','400209','400210','400211','400212'/)
character(len=6),dimension(12):: &
  lun_gmi_tmp=(/'400301','400302','400303','400304','400305','400306', &
                '400307','400308','400309','400310','400311','400312'/)

integer:: ilun_gmi_psfc = 4002
integer:: ilun_gmi_tmp = 4003

integer::ilun_lut_amf_clear=477000
integer::ilun_lut_amf_cloud=477001
integer::ilun_lut_amf_ler6d=477010


! =============
! input options
! =============

! ----------------------------
! option 1: SlantColumnDensity
! ----------------------------
! name_option_SlantColumnDensity:
!hqw TEMPO always use NASA, thus changed to a parameter
  character(len=255), parameter::name_option_SlantColumnDensity='NASA'

! -----------------------------
! option 2: TemperaturePressure
! -----------------------------
! name_option_TemperaturePressure:
!   This option will change VCD and Psfc values, but will not affect AMFcalculated
!    GMI: GMI monthly T/P/Psfc
!  GEOS5: GEOS-5 T/P/Psfc

  character(len=255)::name_gmi_dir='./refdata/'
!  character(len=255)::name_geos5_dir='refdata/'
!  character(len=255)::name_geos5_file

! character(len=255)::name_option_TemperaturePressure='GMI'
! TEMPO uses GEOS5 option to read GEOS-CF
 character(len=255)::name_option_TemperaturePressure='GEOS5'

  integer :: nlayers

! debug
  ! ixdebug & itdebug can be overwrite through control file
  integer :: ixdebug=-1800 !set to negative to prevent debug output
  integer :: itdebug=-60 !set to negative to prevent debug output

  integer :: lun_debug_irr=19
  integer :: lun_debug_rad=29
  integer :: lun_debug_clim=39
  integer :: lun_debug_psfc=49
  integer :: lun_debug_ecf=59
  integer :: lun_debug_scdadj=69
  integer :: lun_debug_ocp=79
  integer :: lun_debug_pcldind=89
  integer :: lun_debug_shift=109
  integer :: lun_debug_pflags=119

! GMI as a backup and testing only, TEMPO usually uses GEOS-CF
! change gmi variables to allocatable
!   so that they won't be allocated if not needed
  integer,parameter::gmi_np=72,gmi_nx=288,gmi_ny=181
  real(kind=4),dimension(:), allocatable :: gmi_lon
  real(kind=4),dimension(:), allocatable :: gmi_lat
  real(kind=4),dimension(:,:,:), allocatable :: gmi_Temperature
  real(kind=4),dimension(:,:,:), allocatable :: gmi_Pressure !include psfc
  real(kind=4),dimension(:,:), allocatable :: gmi_TerrainPressure

! GEOS-CF
  ! geos_np = the number of layers, initialized when reading GEOS-CF
  integer :: geos_np=72
  real(kind=4),dimension(:,:,:),pointer::geos_Temperature
  real(kind=4),dimension(:,:,:),pointer::geos_Q
  real(kind=4),dimension(:,:,:),pointer::geos_Pressure

! calculated VCD at the LUT pressure levels
! these will replace vvcd
  real,dimension(npcld)::gmi_vcd        
  real,dimension(npcld)::geos_vcd

! -----------------------------
! option 3: SurfaceReflectivity
! -----------------------------
! name_option_SurfaceReflectivity:
!   Rsfc(Kleipool) vs. Rsfc(BRDF)

!  character(len=255)::name_option_SurfaceReflectivity='Kleipool'
  character(len=255)::name_option_SurfaceReflectivity='BRDF'

! name_kleipool_dir can be changed by control.txt
  character(len=255)::name_kleipool_dir='./refdata/'
  integer,parameter::kleipool_nx=720,kleipool_ny=360

! changed these to pointer so that they won't allocate if not needed
! seems that only 466 is actually used
! changed pointer to allocatable which can be tested with allocated function
  real,dimension(:),pointer :: kleipool_lon, kleipool_lat
  real(kind=4),dimension(:,:),allocatable :: kleipool_SurfaceReflectivity440
  real(kind=4),dimension(:,:),allocatable :: kleipool_SurfaceReflectivity466
  real(kind=4),dimension(:,:),allocatable :: kleipool_SurfaceReflectivity477

  character(len=255)::name_brdf_dir="./"
  character(len=255)::name_brdf_file='empty'
  real(kind=4),dimension(:,:),allocatable::BRDF_SurfaceReflectivity440
  real(kind=4),dimension(:,:),allocatable::BRDF_SurfaceReflectivity466
! SurfaceReflectivity477 is not used
!  real(kind=4),dimension(:,:),pointer::BRDF_SurfaceReflectivity477
! surface windspeed is needed for GLER
   real(kind=4),dimension(:,:),allocatable:: windspeed2m
! surface geopotential is needed for psfc topo correction
   real(kind=4),dimension(:,:),allocatable:: phisurf

! -----------------
! option 4: SnowIce
! -----------------
! name_option_skipSnowocp 
! 1: skip ocp calculation when snowice fraction > min_snowice
! 0: calculte ocp regardless of snowice fraction
  integer:: name_option_skipSnowocp=0

! name_option_SnowIce:
!  Pcld calculations over SnowIce  vs. Pscene calculations over SnowIce
! changed default to Pcld
  character(len=255)::name_option_SnowIce='Pcld'
!  character(len=255)::name_option_SnowIce='Pscene'

  real:: min_snowice=0.05
! ----------------
! option 5: skipECF005ocp
! ----------------
! name_option_skipECF005ocp:
!   Pcld calculations ECF >= min_ecf  vs. ECF >= 0.00
! 1: skip calculate ocp for ecf = (0., min_ecf]
!   i.e., ocp is for ecf>=min_ecf
! 0: calucate ocp for ecf = (0.0, min_ecf]
!   i.e., ocp is for ecf> 0.0
! Note, for small ecf, ocp is highly uncertainty
!   ocp for ecf <= 0.0 is always skipped
  integer::name_option_skipECFminocp=0

  real:: min_ecf=0.05

! whether to replace pcld with pscene for ecf<min_ecf 
  character(len=255)::name_option_MinECF='Pscene' ! 'Pcld'
! when name_option_MinECF='Pscene',set name_option_skipECFminocp=1 saves time
! as pscene will replace pcld for ecf = (0.,min_ecf)
! when name_option_MinECF='Pcld', name_option_skipECFminocp controls
! whether ocp is calculated or skipped for ecf=(0.,min_ecf)

! ----------------
! option 6: adjust total to dry pressure for vvcd
! ---------------
! 1 : will adjust total to dry for vvcd
! 0 : will not adjust, but assume total=dry 
  integer, parameter::name_option_adjdry=1

! -----------------------------------
! option 7: SceneAlbedo/ScenePressure
! -----------------------------------
! name_option_SceneAlbedo:
!   SceneAlbedoAtTerrain: both, yes or no

 character(len=255)::name_option_SceneAlbedoAtTerrain='both'
! character(len=255)::name_option_SceneAlbedoAtTerrain='yes'
! character(len=255)::name_option_SceneAlbedoAtTerrain='no'
! in production mode, m_cal_pscene force this option is forced to 'no'
!-----------------------------------------
! option 8: option_psfc_clear
! use clear or cloud for high-P interp
! to find pressure for AMF*VCD
! 0: Pclr=Psfc for Pcld>Psfc (original default)
! 1: Pclr=Pcld if Pcld>Psfc
!-------------------------------------
 integer :: option_psfc_clear = 0

!-----------------------------------------
! option 9: option_clip_pcld
! whether to clip pcld within [lut_pcld[1], psfc0]
!-----------------------------------------
 character(len=255):: option_clip_pcld='no'

! ===== end of input options =====
!------------
! input data 
!------------
! inp_ variables are from OMCLDO2 product 
! not needed for TEMPO, thus deleted

!----------------
! input NASA SCD
!----------------
  character(len=255)::name_nasa_dir='./'
  character(len=255)::name_nasa_file='empty'
  integer(kind=4)::nasa_NumTimes
  integer(kind=4)::nasa_nXtrack
  real(kind=4),dimension(:,:),pointer::nasa_SlantColumnAmountO2O2
  real(kind=4),dimension(:,:),allocatable::nasa_scduncertainty
  real(kind=4),dimension(:,:),allocatable::nasa_scdrms
  ! adds l2_TerrainPressure 
  real(kind=4),dimension(:,:),pointer::l2_TerrainPressure
  integer(kind=2), dimension(:,:),  pointer::scd_mdqfl

!----------------
! input tracegas diagnostic & log file
! for wavelength shifts
!----------------
  integer:: option_apply_solshift = 1
  integer:: option_apply_radshift = 1
  ! option_apply_waveshift = 0 will not apply wavelength shift
  !                        = 1 will attempt to apply shift
  character(len=255):: name_diaglog_dir = './'
  character(len=255):: name_diaglog_fnm='empty'
  real(kind=4),dimension(:,:),allocatable:: rad_waveshift
  real(kind=4),dimension(:),allocatable:: irr_waveshift
  real(kind=4):: maxradshift = 0.5 ! larger rad_waveshift will be set to 0.
  real(kind=4):: maxirrshift = 0.5 ! larger irr_waveshift will be set to 0.

!------------
! input extra
!------------

! wavelength
  real::w440=440.0 ! nm for cloud fraction calculation
  real::w466=466.0 ! nm for cloud fraction calculation
  real::w477=477.0 ! nm for cloud pressure calculation

! clip threshold
  real, parameter:: ecf_lowclip=-0.2, ecf_highclip=1.5

  real(kind=4),dimension(:,:),pointer::rad_of_irr440  ! radiance/irradiance at 440 nm calculated by "cal_ecf.f90"
  real(kind=4),dimension(:,:),pointer::rad_of_irr466  ! radiance/irradiance at 466 nm calculated by "cal_ecf.f90"
  real(kind=4),dimension(:,:),pointer::rad_of_irr477  ! radiance/irradiance at 477 nm calculated by "cal_ecf.f90"

!-----------
! frequently used variables in ecf, ocp, pscene calculation
!-----------
  ! moved the following into modules 
  !real::alb0,sza0,vza0,raa0,psfc0 ! input values

  real(kind=4),dimension(:,:),pointer::cal_rad_clr,cal_rad_cld,cal_rad_cld440
  real(kind=8),dimension(npsfc)::cal_amf_clr
  real(kind=8),dimension(npcld)::cal_amf_cld
  real(kind=8),dimension(npcld)::cal_ler_amf
  real(kind=8),dimension(nalb)::cal_ler_r466,cal_ler_r440

!---------------
! write outputs
!---------------
  character(len=255)::name_out_dir='./'
  character(len=255)::name_out_ncdf='empty'
  character(len=255)::name_out_txt='OMCDO2N.out'
  real(kind=8),   dimension(:),    pointer::out_Time
!  real(kind=4),   dimension(:),    pointer::out_SecondsInDay
  real(kind=4),   dimension(:,:),  pointer::out_Longitude
  real(kind=4),   dimension(:,:),  pointer::out_Latitude
  real(kind=4),   dimension(:,:),  pointer::out_SolarZenithAngle
  real(kind=4),   dimension(:,:),  pointer::out_ViewingZenithAngle
  real(kind=4),   dimension(:,:),  pointer::out_RelativeAzimuthAngle
!  integer(kind=2),dimension(:,:),  pointer::out_GroundPixelQualityFlags
!  integer(kind=1),dimension(:,:),  pointer::out_XTrackQualityFlags
!  integer(kind=2),dimension(:),    pointer::out_MeasurementQualityFlags
! 2- or 4-byte out_ProcessingQualityFlags determined by pflag_nbyte 
!  only bit00-15 are currently used
! as the leftest bit for an integer indicate sign and
! all 16 bits are used, apps should not be interpreted as an integer
! but to check individual bits instead,
! to resolve this, it was changed to 4-byte, however,
! total ozone and ozone profile uses 2-byte,to avoid problems there, 
! pflag_nbyte parameter is used as a temporary workaround
  integer(kind=pflag_nbyte),dimension(:,:), pointer::out_ProcessingQualityFlags
  integer(kind=pflag_nbyte),dimension(:,:), pointer::prev_processingflags
  real(kind=4),dimension(:,:),pointer::out_SlantColumnAmountO2O2
  real(kind=4),dimension(:,:),pointer::out_SlantColumnSceneO2O2
  real(kind=4),dimension(:,:),pointer::out_SlantColumnTerrainO2O2
! out_TerrainPressure now holds calculated cpp using LER466 in pscene
  real(kind=4),dimension(:,:),pointer::out_TerrainPressure
  real(kind=4),dimension(:,:),pointer::out_TerrainHeight
!  RAD terrain height is short integer, transfer to real out_TerrainHeight
  real(kind=4),dimension(:,:),pointer::out_SurfaceReflectivity440
  real(kind=4),dimension(:,:),pointer::out_SurfaceReflectivity466
!  integer(kind=2),dimension(:,:),pointer::out_LandAreaFraction

  integer(kind=4)::out_NumTimes
  integer(kind=4)::out_nXtrack
! changed cloud fraction and cloud pressure to real 
  real(kind=4),dimension(:,:),pointer::out_EffectiveCloudFraction
  real(kind=4),dimension(:,:),pointer::out_EffectiveCloudFractionNotClipped
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFraction440
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFractionNotClipped440
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFraction466
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFractionNotClipped466

  real(kind=4),dimension(:,:),pointer::out_CloudPressure
  real(kind=4),dimension(:,:),pointer::out_CloudPressureNotClipped

  real(kind=4),dimension(:,:),pointer::out_SurfaceLER466
  real(kind=4),dimension(:,:),pointer::out_SurfaceLER440
  real(kind=4),dimension(:,:),pointer::out_SceneLER466
  real(kind=4),dimension(:,:),pointer::out_SceneLER440
  real(kind=4),dimension(:,:),pointer::out_ScenePressure
  real(kind=4),dimension(:,:),pointer::out_ReflectanceFactor
! added temperature variables for SCD T correction
  real(kind=4),dimension(:,:),pointer::out_O2O2CloudTemperature
  real(kind=4),dimension(:,:),pointer::out_O2O2SceneTemperature
  real(kind=4),dimension(:,:),pointer::out_O2O2TerrainTemperature

  real, parameter ::fFillValue=-1.2676506E30
  integer, parameter ::iFillValue=-32767

!-------------
! scd filtering and destriping
!-------------
  integer :: option_scdfullfilter = 0
  integer :: option_destripe_scd = 0
  character(len=255)::name_desfac_dir='./'
  character(len=255)::name_desfac_fnm='test_desfac.txt'
  integer:: lun_desfac_fnm=177

  real(kind=4),dimension(:),allocatable:: scddes_hour

!-------------
! gmeta 
!-------------
type gmeta
  character(len=255)::author_affiliation='SAO'
  character(len=255)::author_name='TEMPO STM'
  character(len= 3)::DayNightFlag='Day' 
  character(len= 12)::platformShortName='Intelsat 40e'
  character(len=12)::omiwindow='VIS' 
  character(len=19)::ProcessingCenter='SAO'
  character(len= 1)::ProcessingLevel='2'
  character(len= 9)::InstrumentName='TEMPO'
  character(len= 7)::APPShortName='TEMPOCLDO4'
  character(len= 7)::APPVersion='1.0.0.0'
  character(len=255)::localgranID='empty'
  character(len=48)::APPLongName='TEMPO Cloud O4 Product 1-Orbit L2 Swath'
  character(len=11)::HDFVersion='empty'
  character(len=50)::parameterdescription='Geophysical Cloud Parameters'
  character(len= 3)::omi_collection='empty'
  character(len= 8)::starttime='00:00:00'
  character(len= 8)::endtime='00:00:00'
  character(len=10)::startdate='2023_08_02'
  character(len=10)::enddate='2023_08_02'
  real(kind=4) :: geospatial_lon_min=-180.
  real(kind=4) :: geospatial_lon_max=180.
  real(kind=4) :: geospatial_lat_min=-90.
  real(kind=4) :: geospatial_lat_max=90.
  character(len=13)::leadscientist='hwang'
  character(len=23)::Swathname = 'Cloud O4 Product'
  character(len=32) :: apriori_source = 'empty'
  integer :: granule_year=0, granule_month=0,granule_day=0
  integer:: granule_hour_start=0,granule_minute_start=0,granule_seconds_start=0
  integer:: granule_hour_end=0,granule_minute_end=0,granule_seconds_end=0
  integer(kind=4) :: scan_num = 0
  integer(kind=4) :: granule_num = 0
  real(kind=8)::tai 
end type gmeta

! moved gmetadata def from OMCDO2N.f90 here
type (gmeta) :: gmetadata

!TEMPO time stamp
real(kind=8)::gmeta_tai

!*****************
end module m_vars
!*****************
