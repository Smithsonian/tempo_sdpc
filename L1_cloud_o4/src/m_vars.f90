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
!---------------------------------------------------------------------72

  implicit none

!------------
! misc
!-----------
! operation mode ='production' or 'development', set through control
   character(len=255):: run_mode='development'

!-------------
! allocatable
!-------------
  integer(kind=8),dimension(:),allocatable::edge
  integer(kind=8),dimension(:),allocatable::start
  integer(kind=8),dimension(:),allocatable::stride

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
! reference Earth-Sun distance on 12/22/2014 !1.4715342E11m
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
  integer(kind=4)::rad_nWavelCoef
  real::rad_EarthSunDist

!hqw use out_TerrainHeight in reading for each option
!  integer(kind=2),dimension(:,:),  pointer::rad_TerrainHeight
!OMI uses GroundPixelQualityFlags to decide snow/ice
!TEMPO uses L1 snow_ice_fraction instead
!  integer(kind=4),dimension(:,:),  pointer::rad_GroundPixelQualityFlags
! hqw adds rad_SnowIceFraction
  real(kind=4),   dimension(:,:),   pointer::rad_SnowIceFraction

!hqw added rad_466nm,rad_477nm,rad_440nm which is what needed
  real(kind=4), dimension(:,:), pointer::rad_466nm,rad_477nm,rad_440nm

!--------------
! Lookup table
!--------------
!hqw these initial names can be changed through control file
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

  !hqw LUT cloud node
  ! ALB(cloud) = 0.8, Psfc(cloud) = 700hPa
  integer, parameter:: LUT466rad_cloud_albid = 18
  integer, parameter:: LUT466rad_cloud_psfcid = 18
  integer, parameter:: LUT440rad_cloud_albid = 18
  integer, parameter:: LUT440rad_cloud_psfcid = 18

  !hqw LUT albedo node
  ! ALB=0.0, 0.1, 0.2 used to calc tran & sbar in pscene
  integer, parameter:: LUT_ALBID_0p0 = 1
  integer, parameter:: LUT_ALBID_0p1 = 7
  integer, parameter:: LUT_ALBID_0p2 = 12

  !old LUT dimension
  !integer,parameter::nalb=20, nsza=30, nvza=19, nraa=37
  !integer,parameter::npsfc=23, npcld=23, nrsfc=23
  !new LUT dimension
  integer,parameter:: nalb=20, nsza=30, nvza=25, nraa=37
  integer,parameter:: npsfc=23, npcld=23

  !hqw added the following to remove hard_code in various places
  ! these are limited by oldLUTs
  !real(kind=4),parameter::max_SZA=89.,max_VZA=72.
  ! these are limited by newLUTs
  real(kind=4),parameter:: max_SZA=89., max_VZA=89.

!------------------------
! hqw O4 SCD temperature correction coefficients
! y(T2) = a * y(T1) + b; T1 = 273K, T2=203, 233, 253, 293K
!------------------------
!   real, parameter:: a203=0.8423, b203=-2.0170e-2
!   real, parameter:: a233=0.9318, b233=-1.3336e-2
!   real, parameter:: a253=0.9680, b253=-5.8256e-3
!   real, parameter:: a293=1.0270, b293=-2.4064e-3
! coefs updates to account for changes associated with RJH HITRAN2020 H2O
!  the coefs are derived using Thalman O4 and OMC4 20050701,20060101,20060715
   real, parameter:: TrefO4 = 273. 
   real, parameter:: a203 = 0.90, b203 = -0.03
   real, parameter:: a233 = 0.96, b233 = -0.02
   real, parameter:: a253 = 0.97, b253 = 0.00
   real, parameter:: a293 = 1.01, b293 = 0.00
   ! maximum number of iteration for SCD temperature adjustment
   integer, parameter :: max_scd_iter = 20
   ! if dT < dt_threshold, then stop iteration
   real, parameter :: dt_threshold = 1.0 !K 

!-------------------------
! vertical column density
!-------------------------
  integer,parameter::nvcd=npcld
!hqw vvcd here is replaced with actual gmi_vcd//geos_vcd 
  real,dimension(nvcd):: &
    vvcd=(/0.00472129,0.00648191,0.00889265,0.0121931,0.0166994, &
           0.0228845, 0.0313553, 0.0429509, 0.0588269,0.0805602, &
           0.109367,  0.146245,  0.193142,  0.252468, 0.326837,  &
           0.419688,  0.534917,  0.677226,  0.852079, 1.06600,  &
           1.32557,1.41490,1.53974/)

 !hqw added the multiplicative conversion factor for calculating O4 VCD
 !    this removes hardcoded constant in many routines
  real,parameter::vcd_convfac = 6.733e-4 !previously = 6.765e-4
 !E. Yang suggested change to 6.733e-4 to be more accurate

!-----------
! input LUN
!-----------
character(len=6)::lun_lut_amf_clear='477000'
character(len=6)::lun_lut_amf_cloud='477001'
character(len=6)::lun_lut_amf_ler6d='477010'
!hqw month dependent GMI files should keep the lun_ stuff
character(len=6),dimension(12):: &
  lun_gmi_psfc=(/'400201','400202','400203','400204','400205','400206', &
                 '400207','400208','400209','400210','400211','400212'/)
character(len=6),dimension(12):: &
  lun_gmi_tmp=(/'400301','400302','400303','400304','400305','400306', &
                '400307','400308','400309','400310','400311','400312'/)

integer::ilun_lut_amf_clear=477000
integer::ilun_lut_amf_cloud=477001
integer::ilun_lut_amf_ler6d=477010
integer,dimension(12):: &
  ilun_gmi_psfc=(/400201,400202,400203,400204,400205,400206, &
                  400207,400208,400209,400210,400211,400212/)
integer,dimension(12):: &
  ilun_gmi_tmp=(/400301,400302,400303,400304,400305,400306, &
                 400307,400308,400309,400310,400311,400312/)

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

  character(len=255)::name_gmi_dir='refdata/'
! geos5 is replaced with geoscf using libclim
!  character(len=255)::name_geos5_dir='refdata/'
!  character(len=255)::name_geos5_file

! character(len=255)::name_option_TemperaturePressure='GMI'
 character(len=255)::name_option_TemperaturePressure='GEOS5'

  integer :: nlayers
  integer :: ixdebug=-1800 !set to negative to prevent writing debug output
  integer :: itdebug=-40

  integer,parameter::gmi_np=72,gmi_nx=288,gmi_ny=181
!hqw change these to allocatable
!   so that they won't be allocated if not needed
!  real,dimension(gmi_nx)::gmi_lon
!  real,dimension(gmi_ny)::gmi_lat
!  real(kind=4),dimension(gmi_nx,gmi_ny,gmi_np)::gmi_Temperature
!  real(kind=4),dimension(gmi_nx,gmi_ny,gmi_np+1)::gmi_Pressure !include Psfc
!  real(kind=4),dimension(gmi_nx,gmi_ny)::gmi_TerrainPressure
  real(kind=4),dimension(:), allocatable :: gmi_lon
  real(kind=4),dimension(:), allocatable :: gmi_lat
  real(kind=4),dimension(:,:,:), allocatable :: gmi_Temperature
  real(kind=4),dimension(:,:,:), allocatable :: gmi_Pressure
  real(kind=4),dimension(:,:), allocatable :: gmi_TerrainPressure
  real::gmi_psfc

  ! geos_np = the number of layers is initialized when reading the GEOS-CF forecast
  integer ::geos_np=0
  real(kind=4),dimension(:,:,:),pointer::geos_Temperature
  real(kind=4),dimension(:,:,:),pointer::geos_Pressure

! calculate VCD at the LUT pressure level
  real,dimension(npcld)::gmi_vcd        
  real,dimension(npcld)::geos_vcd

! -----------------------------
! option 3: SurfaceReflectivity
! -----------------------------
! name_option_SurfaceReflectivity:
!   Rsfc(Kleipool) vs. Rsfc(BRDF)

!name_kleipool_dir can be changed by control.txt
  character(len=255)::name_option_SurfaceReflectivity='Kleipool'
!  character(len=255)::name_option_SurfaceReflectivity='BRDF'

  character(len=255)::name_kleipool_dir='./refdata/'
  integer,parameter::kleipool_nx=720,kleipool_ny=360
!hqw changed these to pointer so that they won't allocate if not needed
!  real,dimension(kleipool_nx)::kleipool_lon
!  real,dimension(kleipool_ny)::kleipool_lat
!  real(kind=4),dimension(kleipool_nx,kleipool_ny)::kleipool_SurfaceReflectivity440
!  real(kind=4),dimension(kleipool_nx,kleipool_ny)::kleipool_SurfaceReflectivity466
!  real(kind=4),dimension(kleipool_nx,kleipool_ny)::kleipool_SurfaceReflectivity477
!hqw seems that only 466 is actually used
!hqw changed pointer to allocatable which can be tested with allocated function
  real,dimension(:),pointer :: kleipool_lon, kleipool_lat
  real(kind=4),dimension(:,:),allocatable :: kleipool_SurfaceReflectivity440
  real(kind=4),dimension(:,:),allocatable :: kleipool_SurfaceReflectivity466
  real(kind=4),dimension(:,:),allocatable :: kleipool_SurfaceReflectivity477

  character(len=255)::name_brdf_dir="./"
  character(len=255)::name_brdf_file='empty'
  real(kind=4),dimension(:,:),allocatable::BRDF_SurfaceReflectivity440
  real(kind=4),dimension(:,:),allocatable::BRDF_SurfaceReflectivity466
!hqw SurfaceReflectivity477 is not used
!  real(kind=4),dimension(:,:),pointer::BRDF_SurfaceReflectivity477
!hqw 2m windspeed is needed for GLER
   real(kind=4),dimension(:,:),allocatable:: windspeed2m

! -----------------
! option 4: SnowIce
! -----------------
! name_option_SnowIce:
!   Pcld calculations over SnowIce  vs. Pscene calculations over SnowIce
!  character(len=255)::name_option_SnowIce='Pcld'
  character(len=255)::name_option_SnowIce='Pscene'

! ----------------
! option 5: ECF005
! ----------------
! name_option_ECF005:
!   Pcld calculations ECF >= 0.05  vs. ECF >= 0.00

  real,parameter::min_ecf=0.05, min_snowice=0.05
!hqw changed default name_option_MinECF from yes to no
!  character(len=255)::name_option_MinECF='yes'
  character(len=255)::name_option_MinECF='no'

! -----------------------------------
! option 7: SceneAlbedo/ScenePressure
! -----------------------------------
! name_option_SceneAlbedo:
!   SceneAlbedoAtTerrain: both, yes or no

 character(len=255)::name_option_SceneAlbedoAtTerrain='both'
! character(len=255)::name_option_SceneAlbedoAtTerrain='yes'
! character(len=255)::name_option_SceneAlbedoAtTerrain='no'
!hqw in production mode, m_cal_pscene force this option to 'no'
!-----------------------------------------
! option 8: option_psfc_clear
!hqw moved option_psfc_clear from m_cal_ocp here
! clear/cloud for high-P interp
! find pressure for AMF*VCD
! 0: Pclr=Psfc (fixed) & Pcld>Psfc (original default)
! 1: Pclr=Pcld if Pcld>Psfc
!-------------------------------------
 integer :: option_psfc_clear = 0

!-----------------------------------------
! option 9: option_clip_pcld
! hqw adds this option for whether to clip pcld within [lut_pcld[1], psfc0]
!-----------------------------------------
 character(len=255):: option_clip_pcld='no'

! ===== end of input options =====
!------------
! input data 
!------------
!hqw inp_ variables are from OMCLDO2 product, not needed for TEMPO
!    which is no longer needed, thus deleted

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
  !hqw adds l2_TerrainPressure 
  real(kind=4),dimension(:,:),pointer::l2_TerrainPressure
  integer(kind=2), dimension(:,:),  pointer::scd_mdqfl

!------------
! input extra
!------------
! layer
!46 standard layers + 2 bottom layers to extend to 1100 hPa
!  integer,parameter::nlay=48 

! wavelength
  real::w440=440.0 ! nm for cloud fraction calculation
  real::w466=466.0 ! nm for cloud fraction calculation
  real::w477=477.0 ! nm for cloud pressure calculation

  real(kind=4),dimension(:,:),pointer::rad_of_irr440  ! radiance/irradiance at 440 nm calculated by "cal_ecf.f90"
  real(kind=4),dimension(:,:),pointer::rad_of_irr466  ! radiance/irradiance at 466 nm calculated by "cal_ecf.f90"
  real(kind=4),dimension(:,:),pointer::rad_of_irr477  ! radiance/irradiance at 477 nm calculated by "cal_ecf.f90"

!-----------
! calculate
!-----------
  real::alb0,sza0,vza0,raa0,psfc0,rsfc0 ! input values
  real::alb1,sza1,vza1,raa1,psfc1,rsfc1 ! LUT node1 for interpolation
  real::alb2,sza2,vza2,raa2,psfc2,rsfc2 ! LUT node2 for interpolation

  real(kind=4),dimension(:,:),pointer::cal_rad_clr,cal_rad_cld,cal_rad_cld440
  real(kind=8),dimension(npsfc)::cal_amf_clr
  real(kind=8),dimension(npcld)::cal_amf_cld
  real(kind=8),dimension(npcld)::cal_ler_amf
  real(kind=8),dimension(nalb)::cal_ler_r466,cal_ler_r440

!---------------
! write outputs
!---------------
  character(len=255)::name_out_dir='./'
!  character(len=255)::name_out_he5
  character(len=255)::name_out_ncdf='empty'
  character(len=255)::name_out_txt='OMCDO2N.out'
!  character(len=255)::name_out_swath='Cloud Product'
  real(kind=8),   dimension(:),    pointer::out_Time
!  real(kind=4),   dimension(:),    pointer::out_SecondsInDay
  real(kind=4),   dimension(:,:),  pointer::out_Longitude
  real(kind=4),   dimension(:,:),  pointer::out_Latitude
  real(kind=4),   dimension(:,:),  pointer::out_SolarZenithAngle
  real(kind=4),   dimension(:,:),  pointer::out_ViewingZenithAngle
  real(kind=4),   dimension(:,:),  pointer::out_RelativeAzimuthAngle
  integer(kind=2),dimension(:,:),  pointer::out_GroundPixelQualityFlags
!  integer(kind=1),dimension(:,:),  pointer::out_XTrackQualityFlags
!  integer(kind=2),dimension(:),    pointer::out_MeasurementQualityFlags
  integer(kind=2),dimension(:,:),  pointer::out_ProcessingQualityFlags
  real(kind=4),dimension(:,:),pointer::out_SlantColumnAmountO2O2
  real(kind=4),dimension(:,:),pointer::out_SlantColumnSceneO2O2
  real(kind=4),dimension(:,:),pointer::out_SlantColumnTerrainO2O2
  real(kind=4),dimension(:,:),pointer::out_TerrainPressure
!  real(kind=4),dimension(:,:),pointer::out_TerrainPressureStdDev
  real(kind=4),dimension(:,:),pointer::out_TerrainHeight
!  real(kind=4),dimension(:,:),pointer::out_TerrainHeightStdDev
  real(kind=4),dimension(:,:),pointer::out_SurfaceReflectivity440
  real(kind=4),dimension(:,:),pointer::out_SurfaceReflectivity466
!  integer(kind=2),dimension(:,:),pointer::out_LandAreaFraction
  integer(kind=4)::out_NumTimes
  integer(kind=4)::out_nXtrack
!hqw changed to cloudfraction to real and removed STDs
!  integer(kind=2),dimension(:,:),pointer::out_EffectiveCloudFraction
!  integer(kind=2),dimension(:,:),pointer::out_EffectiveCloudFractionNotClipped
!  integer(kind=2),dimension(:,:),pointer::out_EffectiveCloudFractionSTD
!  integer(kind=2),dimension(:,:),pointer::out_CloudRadianceFraction440
!  integer(kind=2),dimension(:,:),pointer::out_CloudRadianceFractionNotClipped440
!  integer(kind=2),dimension(:,:),pointer::out_CloudRadianceFractionSTD440
!  integer(kind=2),dimension(:,:),pointer::out_CloudRadianceFraction466
!  integer(kind=2),dimension(:,:),pointer::out_CloudRadianceFractionNotClipped466
!  integer(kind=2),dimension(:,:),pointer::out_CloudRadianceFractionSTD466
  real(kind=4),dimension(:,:),pointer::out_EffectiveCloudFraction
  real(kind=4),dimension(:,:),pointer::out_EffectiveCloudFractionNotClipped
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFraction440
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFractionNotClipped440
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFraction466
  real(kind=4),dimension(:,:),pointer::out_CloudRadianceFractionNotClipped466

!hqw changed out_CloudPressure to real(kind=4)  
!  integer(kind=2),dimension(:,:),pointer::out_CloudPressure
!  integer(kind=2),dimension(:,:),pointer::out_CloudPressureNotClipped
!  integer(kind=2),dimension(:,:),pointer::out_CloudPressureSTD
  real(kind=4),dimension(:,:),pointer::out_CloudPressure
  real(kind=4),dimension(:,:),pointer::out_CloudPressureNotClipped

  real(kind=4),dimension(:,:),pointer::out_SurfaceLER466
  real(kind=4),dimension(:,:),pointer::out_SurfaceLER440
  real(kind=4),dimension(:,:),pointer::out_SceneLER466
  real(kind=4),dimension(:,:),pointer::out_SceneLER440
  real(kind=4),dimension(:,:),pointer::out_ScenePressure
  real(kind=4),dimension(:,:),pointer::out_ReflectanceFactor
!hqw added temperature variables for SCD T correction
  real(kind=4),dimension(:,:),pointer::out_O2O2CloudTemperature
  real(kind=4),dimension(:,:),pointer::out_O2O2SceneTemperature
  real(kind=4),dimension(:,:),pointer::out_O2O2TerrainTemperature

  real, parameter ::fFillValue=-1.2676506E30
  integer, parameter ::iFillValue=-32767

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
  character(len=48)::APPLongName='TEMPO Cloud Product 1-Orbit L2 Swath'
  character(len=11)::HDFVersion='empty'
  character(len=50)::parameterdescription='Geophysical Cloud Parameters'
  character(len= 3)::omi_collection='empty'
  character(len= 8)::starttime='00:00:00'
  character(len= 8)::endtime='00:00:00'
  character(len=10)::startdate='2023_04_07'
  character(len=10)::enddate='2023_04_07'
  real(kind=4) :: geospatial_lon_min=-180.
  real(kind=4) :: geospatial_lon_max=180.
  real(kind=4) :: geospatial_lat_min=-90.
  real(kind=4) :: geospatial_lat_max=90.
  character(len=13)::leadscientist='Science Team'
  character(len=23)::Swathname = 'Cloud Product'
  character(len=32) :: apriori_source = 'empty'
  integer :: granule_year=0, granule_month=0,granule_day=0
  integer:: granule_hour_start=0,granule_minute_start=0,granule_seconds_start=0
  integer:: granule_hour_end=0,granule_minute_end=0,granule_seconds_end=0
  integer(kind=4) :: scan_num = 0
  integer(kind=4) :: granule_num = 0
  real(kind=8)::tai 
end type gmeta

!hqw moved gmetadata def from OMCDO2N.f90 here
type (gmeta) :: gmetadata

!TEMPO time stamp
real(kind=8)::gmeta_tai

!*****************
end module m_vars
!*****************
