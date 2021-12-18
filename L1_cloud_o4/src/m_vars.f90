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
  character(len=255)::name_irr_file
  character(len=255)::name_irr_swath='Sun Volume VIS Swath'
  real(kind=8)::irr_Time
  real(kind=4)::irr_SecondsInDay
!hqw comments out OMI specific variables 
!  integer(kind=2),dimension(:,:),pointer::irr_IrradianceMantissa
!  integer(kind=2),dimension(:,:),pointer::irr_IrradiancePrecision
!  integer(kind=1),dimension(:,:),pointer::irr_IrradianceExponent
!  integer(kind=2),dimension(:,:),pointer::irr_PixelQualityFlags
!  real(kind=4),   dimension(:,:),pointer::irr_WavelengthCoefficient
!  real(kind=4),   dimension(:,:),pointer::irr_WavelengthCoefficientPrecision
!  integer(kind=2)::irr_WavelengthReferenceColumn
!  integer(kind=2)::irr_MeasurementQualityFlags
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
  character(len=255)::name_rad_file
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
!hqw to save memory, use out_TerrainHeight in reading for each option
!  integer(kind=2),dimension(:,:),  pointer::rad_TerrainHeight
!OMI uses GroundPixelQualityFlags to decide snow/ice
!TEMPO uses L1 snow_ice_fraction instead
!  integer(kind=4),dimension(:,:),  pointer::rad_GroundPixelQualityFlags
!hqw comments out OMI specific variables
!  integer(kind=2),dimension(:,:,:),pointer::rad_RadianceMantissa
!  integer(kind=2),dimension(:,:,:),pointer::rad_RadiancePrecision
!  integer(kind=1),dimension(:,:,:),pointer::rad_RadianceExponent
!  integer(kind=1),dimension(:,:),  pointer::rad_XTrackQualityFlags
!hqw moved rad_Radiance and rad_Wavelength inside read_rad_tio
!  real(kind=4), dimension(:,:,:), pointer :: rad_Radiance, rad_Wavelength
!hqw added rad_466nm,rad_477nm,rad_440nm which is what really needed
  real(kind=4), dimension(:,:), pointer::rad_466nm,rad_477nm,rad_440nm
!  integer(kind=2),dimension(:,:,:),pointer::rad_PixelQualityFlags
  real(kind=4),   dimension(:,:,:),pointer::rad_WavelengthCoefficient
  real(kind=4),   dimension(:,:,:),pointer::rad_WavelengthCoefficientPrecision
  integer(kind=2),dimension(:),    pointer::rad_WavelengthReferenceColumn
  integer(kind=2),dimension(:),    pointer::rad_MeasurementQualityFlags
  integer(kind=4)::rad_NumTimes
  integer(kind=4)::rad_nXtrack
  integer(kind=4)::rad_nWavel
  integer(kind=4)::rad_nWavelCoef
  real::rad_EarthSunDist

  ! hqw adds rad_SnowIceFraction
  real(kind=4),   dimension(:,:),   pointer::rad_SnowIceFraction

!--------------
! Lookup table
!--------------
!
! ALB=(/0.00,0.01,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00/)
! SZA=(/0,5,10,15,20,25,30,34,38,42,46,50,54,57,60,63,66,69,72,75,78,80,82,84,85,86,87,88,88.5,89/)
! VZA=(/0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72/)
! RAA=(/0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95, &
!       100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180/)
! name_LCLD=(/'026','027','028','029','030','031','032','033','034','035',&
!             '036','037','038','039','040','041','042','043','044','045','046','047','048'/)
! name_PSFC=(/'055','065','076','089','104','121','142','166','194','227',&
!             '265','308','357','411','472','541','617','701','795','899','1013','1050','1100'/) 

!hqw these initialized names can be changed through control file
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

  !LUT cloud node
  ! ALB(cloud) = 0.8, Psfc(cloud) = 700hPa
  integer, parameter:: LUT466rad_cloud_albid = 18
  integer, parameter:: LUT466rad_cloud_psfcid = 18
  integer, parameter:: LUT440rad_cloud_albid = 18
  integer, parameter:: LUT440rad_cloud_psfcid = 18

  integer,parameter::nalb=20, nsza=30, nvza=19, nraa=37
  integer,parameter::npsfc=23, npcld=23, nrsfc=23

  !hqw added the following to remove hard_code in various places
  ! these are limited by LUTs
  real(kind=4),parameter::max_SZA=89.,max_VZA=72.

!------------------------
! hqw O4 SCD temperature correction coefficients
! y(T2) = a * y(T1) + b; T1 = 273K, T2=203, 233, 253, 293K
!------------------------
   real, parameter:: a203=0.8423, b203=-2.0170e-2
   real, parameter:: a233=0.9318, b233=-1.3336e-2
   real, parameter:: a253=0.9680, b253=-5.8256e-3
   real, parameter:: a293=1.0270, b293=-2.4064e-3
   ! maximum number of iteration for SCD temperature adjustment
   integer, parameter :: max_scd_iter = 10
   ! if dT < dt_threshold, then stop iteration
   real, parameter :: dt_threshold = 0.5 

!-------------------------
! vertical column density
!-------------------------
  integer,parameter::nvcd=npcld
!hqw comment out lvcd, pvcd, tvcd  which are not used in program
!  integer,dimension(nvcd):: &
!    lvcd=(/ 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, &
!            36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48/)
!  real,dimension(nvcd):: &
!    pvcd=(/ 55., 65., 76., 89.,104.,121.,142.,166.,194.,227., &
!           265.,308.,357.,411.,472.,541.,617.,701.,795.,899.,1013.,1050.,1100./)
!  real,dimension(nvcd):: &
!    tvcd=(/216.70,216.70,216.70,216.70,216.70,216.70,216.70,216.70,216.70,216.80, &
!           223.30,229.70,236.20,242.70,249.20,255.70,262.20,268.70,275.20,281.70, &
!    288.20,290.16,292.68/)
!hqw vvcd initialized here is replaced with actual gmi_vcd//dem_vcd//geos_vcd 
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
!character(len=6)::lun_irr_file='400100' !hqw no longer needed
!character(len=6)::lun_kleipool_file='400110'
!character(len=6)::lun_lut_rad440='440000'
!character(len=6)::lun_lut_rad='466000'
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

!integer::ilun_irr_file=400100 !hqw no longer needed
!integer::ilun_kleipool_file=400110
!integer::ilun_lut_rad440=440000
!integer::ilun_lut_rad=466000
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
!   SCD: KNMI vs. NASA 	!KNMI(OMCLDO2)   vs. NASA direct fit

!  character(len=255)::name_option_SlantColumnDensity='KNMI'
!hqw TEMPO always use NASA, thus changed to a parameter
  character(len=255), parameter::name_option_SlantColumnDensity='NASA'

! -----------------------------
! option 2: TemperaturePressure
! -----------------------------
! name_option_TemperaturePressure:
!   This option will change VCD and Psfc values, but will not affect AMFcalculated
!    GMI: GMI monthly T/P/Psfc
!    DEM: GMI monthly T/P and DEM Psfc
!   BDEM: GMI monthly T/P and (GLER) DEM Psfc
!  GEOS5: GEOS-5 T/P/Psfc
!   KNMI: US Standard T/P and DEM Psfc

  character(len=255)::name_gmi_dir='refdata/'
! geos5 is replaced with geoscf using libclim
!  character(len=255)::name_geos5_dir='refdata/'
!  character(len=255)::name_geos5_file

! character(len=255)::name_option_TemperaturePressure='GMI'
! character(len=255)::name_option_TemperaturePressure='DEM'
! character(len=255)::name_option_TemperaturePressure='BDEM'
 character(len=255)::name_option_TemperaturePressure='GEOS5'
! TEMPO does not use KNMI
! character(len=255)::name_option_TemperaturePressure='KNMI'

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
  
  real(kind=4),dimension(:,:),pointer::BDEM_TerrainPressure
  real(kind=4),dimension(:,:),pointer::BDEM_TerrainPressureStdDev
  real(kind=4),dimension(:,:),pointer::BDEM_TerrainHeight
  real(kind=4),dimension(:,:),pointer::BDEM_TerrainHeightStdDev
!  integer,dimension(:,:),pointer::BDEM_LandAreaFraction
!  integer(kind=2),dimension(:,:),pointer::BDEM_LandAreaFraction

  integer,parameter::geos_np=72
  real(kind=4),dimension(:,:,:),pointer::geos_Temperature
  real(kind=4),dimension(:,:,:),pointer::geos_Pressure

! calculate VCD at the LUT pressure level
  real,dimension(npcld)::gmi_vcd        
  real,dimension(npcld)::dem_vcd
  real,dimension(npcld)::geos_vcd
  real,dimension(npcld)::knmi_vcd

! -----------------------------
! option 3: SurfaceReflectivity
! -----------------------------
! name_option_SurfaceReflectivity:
!   Rsfc(OMCLDO2)  vs. Rsfc(Kleipool) vs. Rsfc(BRDF)

!name_kleipool_dir can be changed by control.txt
  character(len=255)::name_kleipool_dir='./refdata/'
!hqw OMCLDO2 is not an option for TEMPO
!  character(len=255)::name_option_SurfaceReflectivity='OMCLDO2'
  character(len=255)::name_option_SurfaceReflectivity='Kleipool'
!  character(len=255)::name_option_SurfaceReflectivity='BRDF'
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
  character(len=255)::name_brdf_file
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

  real,parameter::min_ecf=0.05
  character(len=255)::name_option_MinECF='yes'
!  character(len=255)::name_option_MinECF='no'

!hqw newKNMI is not an option for TEMPO
! -----------------
! option 6: NewKNMI
! -----------------
! name_option_NewKNMI:
!   read NewKNMI?: yes or no

!character(len=255)::name_option_NewKNMI='yes'
! character(len=255)::name_option_NewKNMI='no'

! -----------------------------------
! option 7: SceneAlbedo/ScenePressure
! -----------------------------------
! name_option_SceneAlbedo:
!   SceneAlbedoAtTerrain: both, yes or no

 character(len=255)::name_option_SceneAlbedoAtTerrain='both'
! character(len=255)::name_option_SceneAlbedoAtTerrain='yes'
! character(len=255)::name_option_SceneAlbedoAtTerrain='no'

! -----------------------------------
! Test: zoom-mode
! -----------------------------------
! 0 (global mode); 15 (zoom mode)
!hqw TEMPO does not use izoom, keep it 0
!  integer(kind=4), parameter::izoom=0
 
! ===== end of input options =====
!hqw inp_ variables are from OMCLDO2 product, not needed for TEMPO
!   they are allocated in m_read_input.f90/read_input
!    which is no longer called
!------------
! input data 
!------------
!  character(len=255)::name_he5_dir='./'
!  character(len=255)::name_he5_file
!  character(len=255)::name_inp_swath='CloudFractionAndPressure'
!  integer(kind=4)::inp_NumTimes
!  integer(kind=4)::inp_nXtrack
!  real(kind=8),dimension(:),  pointer::inp_Time
!  real(kind=4),dimension(:,:),pointer::inp_Longitude
!  real(kind=4),dimension(:,:),pointer::inp_Latitude
!  real(kind=4),dimension(:,:),pointer::inp_CloudFraction
!  real(kind=4),dimension(:,:),pointer::inp_CloudFractionNotClipped
!  real(kind=4),dimension(:,:),pointer::inp_CloudFractionSTD
!  real(kind=4),dimension(:,:),pointer::inp_SlantColumnAmountO2O2
!  real(kind=4),dimension(:,:),pointer::inp_SlantColumnAmountO2O2cf
!  real(kind=4),dimension(:,:),pointer::inp_SceneAlbedo
!  real(kind=4),dimension(:,:),pointer::inp_ScenePressure
!  integer(kind=2),dimension(:,:),pointer::inp_CloudPressure
!  integer(kind=2),dimension(:,:),pointer::inp_CloudPressureNotClipped
!  integer(kind=2),dimension(:,:),pointer::inp_CloudPressureSTD
!  integer(kind=1),dimension(:,:),pointer::inp_TerrainReflectivity
!  integer(kind=2),dimension(:,:),pointer::inp_TerrainPressure

!----------------
! input NASA SCD
!----------------
  character(len=255)::name_nasa_dir='./'
  character(len=255)::name_nasa_file
!  character(len=255)::name_nasa_swath='ColumnAmountO4'
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
  integer,parameter::nlay=48 !46 standard layers + 2 bottom layers to extend to 1100 hPa

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
!  real::tsfc0                           ! surface temperature
  real::cal_ler_rad
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
  character(len=255)::name_out_ncdf
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
  
  integer(kind=2),dimension(:,:),pointer::out_CloudPressure
  integer(kind=2),dimension(:,:),pointer::out_CloudPressureNotClipped
!  integer(kind=2),dimension(:,:),pointer::out_CloudPressureSTD
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
! write gmeta 
!-------------
type gmeta
  character(len=255)::author_affiliation
  character(len=255)::author_name
  character(len= 3)::DayNightFlag 
  character(len= 12)::platformShortName='Intelsat 40e'
  character(len=12)::omiwindow='VIS' 
  character(len=19)::ProcessingCenter='SAO'
  character(len= 1)::ProcessingLevel='2'
  character(len= 9)::InstrumentName='TEMPO'
  character(len= 7)::APPShortName
  character(len= 7)::APPVersion
  character(len=255)::localgranID
  character(len=48)::APPLongName='TEMPO Cloud Product 1-Orbit L2 Swath'
  character(len=11)::HDFVersion
  character(len=50)::parameterdescription='Geophysical Cloud Parameters'
  character(len= 3)::omi_collection
  character(len= 8)::starttime
  character(len= 8)::endtime
  character(len=10)::startdate
  character(len=10)::enddate
  real(kind=4) :: geospatial_lon_min, geospatial_lon_max
  real(kind=4) :: geospatial_lat_min, geospatial_lat_max
  character(len=13)::leadscientist='TEMPO'
  character(len=23)::Swathname = 'Cloud Product'
  integer :: granule_year, granule_month,granule_day
  integer:: granule_hour_start,granule_minute_start,granule_seconds_start
  integer:: granule_hour_end, granule_minute_end,granule_seconds_end
  integer(kind=4) :: scan_num
  integer(kind=4) :: granule_num
  real(kind=8)::tai
end type gmeta

!hqw moved gmetadata def from OMCDO2N.f90 here
type (gmeta) :: gmetadata

real(kind=8)::gmeta_tai

!*****************
end module m_vars
!*****************
