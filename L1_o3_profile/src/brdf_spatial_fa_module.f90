MODULE spatial_fa_module

  USE parameters_module
  USE error_module,         ONLY : ErrorType, CheckError, RaisePixelError, RaiseWarning
  USE interpolation_module, ONLY : SPLINT1, SPLINE1, BSPLINE, BSPLINE_EdgeFill
  USE level1_def,           ONLY : GeolocationType
  USE profile_module,       ONLY : SurfProfType
  USE time_module,          ONLY : TimeType, SET_TIME_FROM_TAU
  USE netcdf_module,        ONLY : CheckNetCDFErrorStatus,ncdf_var_exists

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'
  
  TYPE SpatialFAType
    CHARACTER(LEN=maxChar)       :: Infile
    CHARACTER(LEN=maxChar)       :: RootBRDFDir ! BRDF Climatology 
    CHARACTER(LEN=maxChar)       :: RootDataDir ! Actually the root data directory
    LOGICAL                      :: DoIsotropic
    INTEGER                      :: WhichAlbedo
    LOGICAL                      :: DoOceanGlint
    INTEGER                      :: nkern
    INTEGER                      :: maxpar
    INTEGER,         ALLOCATABLE :: npar(:)
    INTEGER,         ALLOCATABLE :: KernIdx(:)
    REAL(KIND=8),    ALLOCATABLE :: KernPar(:,:)
    INTEGER                      :: ncid
    INTEGER                      :: wmx
    INTEGER                      :: imx
    INTEGER                      :: jmx
    INTEGER                      :: fmx
    INTEGER                      :: bmx ! # bands
    INTEGER                      :: I
    INTEGER                      :: J
    REAL(KIND=8)                 :: PixelLongitude
    REAL(KIND=8)                 :: PixelLatitude
    REAL(KIND=8)                 :: PixelCornerLongitudes(4)
    REAL(KIND=8)                 :: PixelCornerlatitudes(4)
    INTEGER                      :: maxdim
    REAL(KIND=8),    ALLOCATABLE :: WaterAlbedo(:)
    REAL(KIND=8),    ALLOCATABLE :: SnowFineAlbedo(:)
    REAL(KIND=8),    ALLOCATABLE :: SnowCoarseAlbedo(:)
    REAL(KIND=8),    ALLOCATABLE :: SnowAlbedo(:) ! For the given Scene
    REAL(KIND=8),    ALLOCATABLE :: Mu(:)
    REAL(KIND=8),    ALLOCATABLE :: MuSP(:)
    REAL(KIND=8),    ALLOCATABLE :: MuBands(:)
    REAL(KIND=8),    ALLOCATABLE :: Psi(:)
    REAL(KIND=8),    ALLOCATABLE :: Longitude(:)
    REAL(KIND=8),    ALLOCATABLE :: Latitude(:)
    REAL(KIND=8),    ALLOCATABLE :: Wvl(:)
    REAL(KIND=8),    ALLOCATABLE :: W(:,:)
    REAL(KIND=8),    ALLOCATABLE :: FitCoeff(:,:)
    REAL(KIND=8),    ALLOCATABLE :: WSP(:,:)
    REAL(KIND=8),    ALLOCATABLE :: G(:,:)
    REAL(KIND=8),    ALLOCATABLE :: RSR(:,:)
    REAL(KIND=8),    ALLOCATABLE :: BandUncertainty(:)
    REAL(KIND=8),    ALLOCATABLE :: FactorUncertainty(:,:)
    INTEGER(KIND=2), ALLOCATABLE :: FAExists(:,:)
    REAL(KIND=8),    ALLOCATABLE :: KernelAmplitudes(:,:)
    REAL(KIND=8),    ALLOCATABLE :: KernelAmplitudesSP(:,:)
    REAL(KIND=8)                 :: EffectiveSnowFraction
    REAL(KIND=8)                 :: SeaIceFraction
    REAL(KIND=8)                 :: SnowAgeWeight
    REAL(KIND=8)                 :: LandFraction
    INTEGER                      :: DayOfYear
    REAL(KIND=8)                 :: TauTime

    ! Climatology Fields
    INTEGER                      :: BRDF_imx
    INTEGER                      :: BRDF_jmx
    INTEGER                      :: BRDF_tmx
    REAL(KIND=8),    ALLOCATABLE :: BRDF_Lon(:)
    REAL(KIND=8),    ALLOCATABLE :: BRDF_Lat(:)
    REAL(KIND=8),    ALLOCATABLE :: BRDF_Tau(:)
    REAL(KIND=8),    ALLOCATABLE :: BRDF_MinLat(:)
    REAL(KIND=8),    ALLOCATABLE :: BRDF_MaxLat(:)

    ! For climatology times
    LOGICAL                      :: IsClimatology
    REAL(KIND=8),    ALLOCATABLE :: BRDF_DayDbl(:)
    INTEGER,         ALLOCATABLE :: BRDF_DayIdx(:)

    REAL(KIND=8),    ALLOCATABLE :: PixelBRDF(:,:)
    REAL(KIND=8),    ALLOCATABLE :: SubPixBRDF(:,:,:,:)
    INTEGER                      :: SubPixImx
    INTEGER                      :: SubPixJmx
    REAL(KIND=8),    ALLOCATABLE :: SubPixLon(:)
    REAL(KIND=8),    ALLOCATABLE :: SubPixLat(:)
    REAL(KIND=8),    ALLOCATABLE :: SubPixLandFrac(:,:)
    
    ! Constant open/closing files lead to errors - need to keep 
    ! current file open and only update when a new one is needed
    INTEGER              :: current_tidx
    INTEGER, ALLOCATABLE :: obs_ncid(:)
    INTEGER, ALLOCATABLE :: obs_vid(:,:)


    ! Land-Water Flag Database


    ! INTEGER                      :: MODIS_imx
    ! INTEGER                      :: MODIS_jmx
    ! REAL(KIND=8)                 :: PixelBRDF(4,3)
    ! REAL(KIND=8),    ALLOCATABLE :: MODIS_BRDF(:,:,:,:)
    ! REAL(KIND=8),    ALLOCATABLE :: MODIS_Lon(:)
    ! REAL(KIND=8),    ALLOCATABLE :: MODIS_Lat(:)
    ! REAL(KIND=8),    ALLOCATABLE :: MODIS_LandFrac(:,:)
    ! INTEGER(KIND=1), ALLOCATABLE :: MODIS_fld_i1(:,:)
    ! INTEGER(KIND=2), ALLOCATABLE :: MODIS_fld_i2(:,:)
  ENDTYPE SpatialFAType

  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'spatial_fa_module'
  PRIVATE :: ModuleName

  CONTAINS 

  SUBROUTINE InitSpatialFA(SpatialFA,Error)

!#if defined( USE_VLIDORT2p7) || defined(USE_VLIDORTpca )
    USE VLIDORT_PARS,         ONLY : NewCMGLINT_IDX
!#elif defined( USE_VLIDORT2p8 )
!    USE VLIDORT_PARS_M,       ONLY : NewCMGLINT_IDX
!#endif

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(SpatialFAType) :: SpatialFA
    TYPE(ErrorType)     :: Error

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                :: rcode, vid, dimid, dimid2(2), usgs_ncid, clim_ncid
    CHARACTER(LEN=maxChar) :: tmpchar
    CHARACTER(LEN=maxChar) :: usgs_infile
    
    ! For reading Archived float data
    REAL(KIND=4), ALLOCATABLE :: tmparr_r4(:)
    REAL(KIND=4)              :: usgs_wvl(480), usgs_water(480)
    REAL(KIND=4)              :: usgs_snow_fine(480), usgs_snow_coarse(480)
    INTEGER                   :: maxdim
    INTEGER(KIND=2)           :: is_clim(1)
    TYPE(TimeType)            :: TimeStruct
    ! Wavelength Grid for relative spectral response functions
    INTEGER                   :: n_rsr, b, i
    REAL(KIND=8), ALLOCATABLE :: wvl_rsr(:), rsr(:,:)

    ! Weights for regridding RSR
    REAL(KIND=8), ALLOCATABLE :: rsr_wt(:)
    REAL(KIND=8)              :: rsr_sum
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'InitSpatialFA'
    LOGICAL                     :: oob_low, oob_hi
    INTEGER                     :: errstat
    
    ! =====================================================================
    ! InitialSpatialFA starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Figure out the number of kernels
    IF( SpatialFA%DoIsotropic ) THEN
      SpatialFA%nkern = 1
      SpatialFA%maxpar = 1 ! Actually zero but 
      ALLOCATE(SpatialFA%KernIdx(SpatialFA%nkern)) ; SpatialFA%KernIdx(1) = 1
      ALLOCATE(SpatialFA%npar(SpatialFA%nkern))    ; SpatialFA%npar(1) = 0
    ELSE
      IF( SpatialFA%DoOceanGlint ) THEN
        SpatialFA%nkern = 4
        SpatialFA%maxpar = 2
        ALLOCATE(SpatialFA%KernIdx(SpatialFA%nkern)) ; SpatialFA%KernIdx(:) = (/1,3,4,NewCMGLINT_IDX/)
        ALLOCATE(SpatialFA%npar(SpatialFA%nkern))    ; SpatialFA%npar(:) = (/0,0,2,2/)
        ALLOCATE(SpatialFA%kernPar(SpatialFA%maxpar,SpatialFA%nkern)) ; SpatialFA%kernPar(:,:) = 0.0d0
        SpatialFA%kernPar(1,3) = 2.0d0 ; SpatialFA%kernPar(2,3) = 1.0d0
        SpatialFA%kernPar(1,4) = 0.0d0 ; SpatialFA%kernPar(2,4) = 0.0d0 ! Must be supplied later by windspeed
      ELSE
        SpatialFA%nkern = 3
        SpatialFA%maxpar = 2
        ALLOCATE(SpatialFA%KernIdx(SpatialFA%nkern)) ; SpatialFA%KernIdx(:) = (/1,3,4/)
        ALLOCATE(SpatialFA%npar(SpatialFA%nkern))    ; SpatialFA%npar(:) = (/0,0,2/)
        ALLOCATE(SpatialFA%kernPar(SpatialFA%maxpar,SpatialFA%nkern)) ; SpatialFA%kernPar(:,:) = 0.0d0
        SpatialFA%kernPar(1,3) = 2.0d0 ; SpatialFA%kernPar(2,3) = 1.0d0
      ENDIF
    ENDIF
    
    ! Initialize the input geolocation to some nonsense
    SpatialFA%PixelLongitude = -1e30
    SpatialFA%PixelLatitude = -1e30
    SpatialFA%PixelCornerLongitudes(:) = -1e30
    SpatialFA%PixelCornerlatitudes(:)  = -1e30
    SpatialFA%maxdim = 1

    ! --------------------------------------------------------------
    ! USGS Water Spectra
    ! --------------------------------------------------------------
    
    ! Path to data
    usgs_infile = trim(adjustl(SpatialFA%RootDataDir)) // '/BRDF_EOF/AlbSpec/' // &
                  'usgs_snow_water.nc'
    !print * , usgs_infile
    ! Attach file
    rcode = nf_open( trim(adjustl(usgs_infile)), nf_Share, usgs_ncid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs:open')

    ! Read Wavelength
    rcode = nf_inq_varid( usgs_ncid, 'Wavelength', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_wvl:inq_varid')
    rcode = nf_get_var_real( usgs_ncid, vid, usgs_wvl )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_wvl:get_var')

    ! Read Water Spectrum
    rcode = nf_inq_varid( usgs_ncid, 'Water', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_water:inq_varid')
    rcode = nf_get_var_real( usgs_ncid, vid, usgs_water )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_water:get_var')
    
    ! Read Fine Snow Spectrum
    rcode = nf_inq_varid( usgs_ncid, 'SnowFine', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_snowf:inq_varid')
    rcode = nf_get_var_real( usgs_ncid, vid, usgs_snow_fine )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_snowf:get_var')
    
    ! Read Coarse Snow Spectrum
    rcode = nf_inq_varid( usgs_ncid, 'SnowCoarse', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_snowc:inq_varid')
    rcode = nf_get_var_real( usgs_ncid, vid, usgs_snow_coarse )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs_snowc:get_var')
    
    ! Close file
    rcode = nf_close( usgs_ncid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'usgs:close')

    ! --------------------------------------------------------------
    ! Reflectance BRDF database
    ! --------------------------------------------------------------

    ! Attach file
    rcode = nf_open(TRIM(ADJUSTL(SpatialFA%RootBRDFDir))//'/brdf_grid_info.nc', nf_Share, clim_ncid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim:open')

    ! Dimension of the relative spectral response functions
    rcode = nf_inq_varid(clim_ncid, 'wvl',    vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_wvl:inq_varid')
    rcode = nf_inq_vardimid(clim_ncid,   vid,  dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_wvl:inq_vardimid')
    rcode = nf_inq_dim(clim_ncid, dimid,tmpchar,n_rsr)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_wvl:inq_dim')

    ! Read RSR wavelength grid
    ALLOCATE(wvl_rsr(n_rsr)) ; rcode = nf_get_var_double( clim_ncid, vid, wvl_rsr )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_wvl:get_var')

    ! Number of observed wavelength bands
    rcode = nf_inq_varid(clim_ncid, 'rsr',    vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_rsr:inq_varid')
    rcode = nf_inq_vardimid(clim_ncid,   vid,  dimid2)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_rsr:inq_vardimid')
    rcode = nf_inq_dim(clim_ncid, dimid2(2),tmpchar,SpatialFA%bmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_rsr:inq_dim')

    ! Read RSR
    ALLOCATE(rsr(n_rsr,SpatialFA%bmx)) ; rcode = nf_get_var_double( clim_ncid, vid, rsr )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_rsr:get_var')

    ! Band Uncertainty
    ALLOCATE(SpatialFA%BandUncertainty(SpatialFA%bmx)) ; SpatialFA%BandUncertainty(SpatialFA%bmx) = 0.0d0
    IF(ncdf_var_exists(clim_ncid,'band_uncertainty')) THEN
      rcode = nf_inq_varid(clim_ncid, 'band_uncertainty',vid)
      rcode = nf_get_var_double( clim_ncid, vid,SpatialFA%BandUncertainty)
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_bnduncert:get_var')
    ENDIF
    
    ! Read Longitude
    rcode = nf_inq_varid(clim_ncid, 'lon',    vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lon:inq_varid')
    rcode = nf_inq_vardimid(clim_ncid,   vid,  dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lon:inq_vardimid')
    rcode = nf_inq_dim(clim_ncid, dimid,tmpchar,SpatialFA%BRDF_imx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lon:inq_dim')
    ALLOCATE(SpatialFA%BRDF_Lon(SpatialFA%BRDF_imx))
    rcode = nf_get_var_double( clim_ncid, vid, SpatialFA%BRDF_Lon )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lon:get_var')

    ! Read Latitude 
    rcode = nf_inq_varid(clim_ncid, 'lat', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat:inq_varid')
    rcode = nf_inq_vardimid(clim_ncid,   vid,  dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat:inq_vardimid')
    rcode = nf_inq_dim(clim_ncid, dimid,tmpchar,SpatialFA%BRDF_jmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat:inq_dim')
    ALLOCATE(SpatialFA%BRDF_Lat(SpatialFA%BRDF_jmx))
    rcode = nf_get_var_double( clim_ncid, vid, SpatialFA%BRDF_Lat )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat:get_var')

    ! Read Time
    rcode = nf_inq_varid(clim_ncid, 'time', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_time:inq_varid')
    rcode = nf_inq_vardimid(clim_ncid,   vid,  dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_time:inq_vardimid')
    rcode = nf_inq_dim(clim_ncid, dimid,tmpchar,SpatialFA%BRDF_tmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_time:inq_dim')
    ALLOCATE(SpatialFA%BRDF_Tau(SpatialFA%BRDF_tmx))
    rcode = nf_get_var_double( clim_ncid, vid, SpatialFA%BRDF_Tau )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_time:get_var')

    ! Read Latitude coverage of fields
    rcode = nf_inq_varid(clim_ncid, 'lat_min', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat_min:inq_varid')
    ALLOCATE(SpatialFA%BRDF_MinLat(SpatialFA%BRDF_tmx))
    rcode = nf_get_var_double( clim_ncid, vid, SpatialFA%BRDF_MinLat )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat_min:get_var')
    rcode = nf_inq_varid(clim_ncid, 'lat_max', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat_max:inq_varid')
    ALLOCATE(SpatialFA%BRDF_MaxLat(SpatialFA%BRDF_tmx))
    rcode = nf_get_var_double( clim_ncid, vid, SpatialFA%BRDF_MaxLat )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_lat_max:get_var')

    ! Check if timestamps are for climatology or exact
    rcode = nf_inq_varid(clim_ncid,'is_clim', vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_is_clim:inq_varid')
    rcode = nf_get_var_int2(clim_ncid,vid,is_clim)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim_is_clim:get_var')
    IF(is_clim(1) .GT. 0) THEN
      SpatialFA%IsClimatology = .TRUE.

      ! Populate the day index arrays
      ALLOCATE(SpatialFA%BRDF_DayDbl(SpatialFA%BRDF_tmx+4))
      ALLOCATE(SpatialFA%BRDF_DayIdx(SpatialFA%BRDF_tmx+4))

      DO I=1,SpatialFA%BRDF_tmx

        ! Compute DOY
        CALL SET_TIME_FROM_TAU(SpatialFA%BRDF_Tau(I),TimeStruct,Error)

        ! Set Value
        SpatialFA%BRDF_DayDbl(I+2) = REAL(TimeStruct%DayOfYear,KIND=8)
        SpatialFA%BRDF_DayIdx(I+2) = I

      ENDDO

      ! Ensure periodicity
      SpatialFA%BRDF_DayIdx(1) = SpatialFA%BRDF_tmx-1
      SpatialFA%BRDF_DayIdx(2) = SpatialFA%BRDF_tmx
      SpatialFA%BRDF_DayDbl(1) = SpatialFA%BRDF_DayDbl(SpatialFA%BRDF_tmx+1)-365.0d0
      SpatialFA%BRDF_DayDbl(2) = SpatialFA%BRDF_DayDbl(SpatialFA%BRDF_tmx+2)-365.0d0
      SpatialFA%BRDF_DayIdx(SpatialFA%BRDF_tmx+3) = 1
      SpatialFA%BRDF_DayIdx(SpatialFA%BRDF_tmx+4) = 2
      SpatialFA%BRDF_DayIdx(SpatialFA%BRDF_tmx+3) = SpatialFA%BRDF_DayDbl(3)+365.0d0
      SpatialFA%BRDF_DayIdx(SpatialFA%BRDF_tmx+4) = SpatialFA%BRDF_DayDbl(4)+365.0d0

    ELSE
      SpatialFA%IsClimatology = .FALSE.
    ENDIF

    ! Close file
    rcode = nf_close( clim_ncid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'brdf_clim:close')
    
    ! Allocate Arrays for Pixel BRDF
    ALLOCATE(SpatialFA%PixelBRDF(SpatialFA%bmx,3))

    ! --------------------------------------------------------------
    ! Factor analysis results
    ! --------------------------------------------------------------
 
    ! Attach file
    rcode = nf_open(trim(adjustl(SpatialFA%RootDataDir))//trim(adjustl(SpatialFA%Infile)), &
            nf_Share, SpatialFA%ncid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa:open')

    ! ---------------------------------------------------------------
    ! Get the dimensions to allocate arrays
    ! ---------------------------------------------------------------
    
    ! Find the dimension of the spectrum
    rcode = nf_inq_varid(SpatialFA%ncid, 'wvl',    vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_wvl:inq_varid')
    rcode = nf_inq_vardimid(SpatialFA%ncid,   vid,  dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_wvl:inq_vardimid')
    rcode = nf_inq_dim(SpatialFA%ncid, dimid,tmpchar,SpatialFA%wmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_wvl:inq_dim')
    maxdim = SpatialFA%wmx
    
    ! Find the dimension of longitudes
    rcode = nf_inq_varid(SpatialFA%ncid, 'lon',    vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lon:inq_varid')
    rcode = nf_inq_vardimid(SpatialFA%ncid,   vid,  dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lon:inq_vardimid')
    rcode = nf_inq_dim(SpatialFA%ncid, dimid,tmpchar,SpatialFA%imx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lon:inq_dim')
    maxdim = MAX(maxdim,SpatialFA%imx)
    
    ! Find the dimension of latitudes
    rcode = nf_inq_varid(SpatialFA%ncid, 'lat',    vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lat:inq_varid')
    rcode = nf_inq_vardimid(SpatialFA%ncid,   vid,  dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lat:inq_vardimid')
    rcode = nf_inq_dim(SpatialFA%ncid, dimid,tmpchar,SpatialFA%jmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lat:inq_dim')
    maxdim = MAX(maxdim,SpatialFA%jmx)
    
    ! Find the dimension for the factors
    rcode = nf_inq_varid(SpatialFA%ncid, 'W_global',    vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_W:inq_varid')
    rcode = nf_inq_vardimid(SpatialFA%ncid,   vid,  dimid2)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_W:inq_vardimid')
    rcode = nf_inq_dim(SpatialFA%ncid, dimid2(1),tmpchar,SpatialFA%fmx)
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_W:inq_dim')
    maxdim = MAX(maxdim,SpatialFA%fmx)
    
    ! Set max dimension
    SpatialFA%maxdim = maxdim

    ! Allocate arrays
    ALLOCATE(SpatialFA%Mu(SpatialFA%wmx))
    ALLOCATE(SpatialFA%Psi(SpatialFA%wmx))
    ALLOCATE(SpatialFA%MuSP(SpatialFA%wmx))
    ALLOCATE(SpatialFA%MuBands(SpatialFA%bmx))
    ALLOCATE(SpatialFA%Wvl(SpatialFA%wmx))
    ALLOCATE(SpatialFA%Longitude(SpatialFA%imx))
    ALLOCATE(SpatialFA%Latitude(SpatialFA%jmx))
    ALLOCATE(SpatialFA%W(SpatialFA%wmx,SpatialFA%fmx))
    ALLOCATE(SpatialFA%WSP(SpatialFA%wmx,SpatialFA%fmx))
    ALLOCATE(SpatialFA%RSR(SpatialFA%wmx,SpatialFA%bmx))
    ALLOCATE(SpatialFA%G(SpatialFA%fmx,SpatialFA%bmx))
    ALLOCATE(SpatialFA%FAExists(SpatialFA%jmx,SpatialFA%imx))
    ALLOCATE(SpatialFA%WaterAlbedo(SpatialFA%wmx))
    ALLOCATE(SpatialFA%SnowFineAlbedo(SpatialFA%wmx))
    ALLOCATE(SpatialFA%SnowAlbedo(SpatialFA%wmx))
    ALLOCATE(SpatialFA%SnowCoarseAlbedo(SpatialFA%wmx))
    ALLOCATE(SpatialFA%KernelAmplitudes(SpatialFA%wmx,SpatialFA%nkern))
    ALLOCATE(SpatialFA%KernelAmplitudesSP(SpatialFA%wmx,SpatialFA%nkern))
    ALLOCATE(SpatialFA%FactorUncertainty(SpatialFA%fmx,SpatialFA%fmx))
    ALLOCATE(SpatialFA%FitCoeff(SpatialFA%fmx,3))
    
    ! Allocate temporary float array
    ALLOCATE(tmparr_r4(maxdim))
    
    ! Initialize x/y indices
    SpatialFA%I = -1
    SpatialFA%J = -1
    
    ! Initialize other variables
    SpatialFA%EffectiveSnowFraction = 0.0d0
    SpatialFA%SnowAlbedo(:) = 0.0d0

    ! Initialize MODIS read array indices
    SpatialFA%SubPixImx = -1
    SpatialFA%SubPixJmx = -1
    
    ! --------------------------------------------------------------
    ! Read data
    ! --------------------------------------------------------------
    
    ! (1) Wavelength
    rcode = nf_inq_varid( SpatialFA%ncid, 'wvl', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_wvl:inq_varid')
    rcode = nf_get_var_real( SpatialFA%ncid, vid, tmparr_r4(1:SpatialFA%wmx) )
    SpatialFA%Wvl(:) = REAL(tmparr_r4(1:SpatialFA%wmx),KIND=8)
    
    ! (2) Longitudes
    rcode = nf_inq_varid( SpatialFA%ncid, 'lon', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lon:inq_varid')
    rcode = nf_get_var_real( SpatialFA%ncid, vid, tmparr_r4(1:SpatialFA%imx) )
    SpatialFA%Longitude(:) = REAL(tmparr_r4(1:SpatialFA%imx),KIND=8)
    
    ! (3) Latitudes
    rcode = nf_inq_varid( SpatialFA%ncid, 'lat', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_lat:inq_varid')
    rcode = nf_get_var_real( SpatialFA%ncid, vid, tmparr_r4(1:SpatialFA%jmx) )
    SpatialFA%Latitude(:) = REAL(tmparr_r4(1:SpatialFA%jmx),KIND=8)
    
    ! (4) Array to check if FA exists for given lon/lat
    rcode = nf_inq_varid( SpatialFA%ncid, 'yn_fa', vid )
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'fa_yn_fa:inq_varid')
    rcode = nf_get_var_int2( SpatialFA%ncid, vid, SpatialFA%FAExists )
    
    ! Spline USGS Water/Snow albedos to FA grid
    CALL BSPLINE( REAL(usgs_wvl,KIND=8), REAL(usgs_water,KIND=8), 480,&
                  SpatialFA%Wvl, SpatialFA%WaterAlbedo, SpatialFA%wmx,&
                  rcode                                               )
    CALL BSPLINE( REAL(usgs_wvl,KIND=8), REAL(usgs_snow_fine,KIND=8), 480,&
                  SpatialFA%Wvl, SpatialFA%SnowFineAlbedo, SpatialFA%wmx, &
                  rcode                                                   )
    CALL BSPLINE( REAL(usgs_wvl,KIND=8), REAL(usgs_snow_coarse,KIND=8), 480,&
                  SpatialFA%Wvl, SpatialFA%SnowCoarseAlbedo, SpatialFA%wmx, &
                  rcode                                                     )

    ! --------------------------------------------------------------
    ! Compute RSR on wavelength grid
    ! --------------------------------------------------------------

    ! For now just interpolate to FA grid and renormalize
    ! Need to check exact definition of RSR 
    SpatialFA%RSR(:,:)= 0.0d0

    DO b=1,SpatialFA%bmx

      ! Spline Value
      CALL BSPLINE_EdgeFill(wvl_rsr, rsr(:,b), n_rsr,                        & 
                            SpatialFA%Wvl, SpatialFA%RSR(:,b), SpatialFA%wmx,&
                            errstat,oob_low,oob_hi,0.0d0                     )

      ! Renormalize
      rsr_sum = SUM(SpatialFA%RSR(:,b))
      IF(rsr_sum .GT. TINY(0.0d0)) THEN
        SpatialFA%RSR(:,b) = SpatialFA%RSR(:,b) / rsr_sum
      ENDIF

    ENDDO
    
    ! Warn if out of bounds flags have been raised
    ! IF(oob_low) THEN
    !   CALL RaiseWarning( Error, ErrorCode_OptProp, ModuleName, SubroutineName, &
    !                      Message_in='Wavelength Grid in FA > Min. Wavelength in Climatology RSR Grid',&
    !                      Action_in='Check Band overlap with climatology and FA' )
    ! ENDIF

    ! IF(oob_hi) THEN
    !   CALL RaiseWarning( Error, ErrorCode_OptProp, ModuleName, SubroutineName, &
    !                      Message_in='Wavelength Grid in FA < Max. Wavelength in Climatology RSR Grid',&
    !                      Action_in='Check Band overlap with climatology and FA' )
    ! ENDIF
    

    ! Initialize tracking of NetCDF indices
    ! -------------------------------------
    SpatialFA%current_tidx = -1
    ALLOCATE(SpatialFA%obs_ncid(SpatialFA%bmx))  ; SpatialFA%obs_ncid(:) = 0
    ALLOCATE(SpatialFA%obs_vid(SpatialFA%bmx,3)) ; SpatialFA%obs_vid(:,:) = 0

  END SUBROUTINE InitSpatialFA

  SUBROUTINE SampleSpatialFA( Geolocation, SurfProf, SpatialFA, Error )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(GeolocationType), INTENT(IN)      :: Geolocation
    TYPE(SurfProfType),    INTENT(IN)      :: SurfProf
    TYPE(SpatialFAType),   INTENT(INOUT)   :: SpatialFA
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: i, j, n, w, f, rcode, vid
    INTEGER :: tmp(1)
    REAL(KIND=4), ALLOCATABLE :: tmparr_r4(:), tmpW(:,:)
    

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'SampleSpatialFA'

    ! =====================================================================
    ! SampleSpatialFA starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Set windspeed [m/s] for ocean kernel if present
    IF(SpatialFA%DoOceanGlint .AND. .NOT. SpatialFA%DoIsotropic) THEN
      SpatialFA%kernPar(1,4) = SurfProf%WindSpeed
      SpatialFA%kernPar(2,4) = SurfProf%OceanSalinity
    ENDIF

    ! Check if center is out of range
    IF( ABS(Geolocation%Longitude) .GT. 180.0 .OR. ABS(Geolocation%Latitude) .GT. 90.0 ) THEN
      
      ! Load global factors
      i = 0
      j = 0
      
    ELSE
      
      ! Find the nearest Index
      tmp = MINLOC(ABS(SpatialFA%Longitude - Geolocation%Longitude)) ; i = tmp(1)
      tmp = MINLOC(ABS(SpatialFA%Latitude - Geolocation%Latitude))   ; j = tmp(1)
      
      ! Check if local FA results exist
      IF( SpatialFA%FAExists(j,i) .EQ. 0 ) THEN
        i = 0
        j = 0
      ENDIF
      
    ENDIF
    
    ! Check if case is already loaded
    IF( SpatialFA%I .NE. i .AND. SpatialFA%J .NE. j ) THEN
      
      ! Allocate temporary read array
      ALLOCATE(tmparr_r4(SpatialFA%maxdim))
      ALLOCATE(tmpW(SpatialFA%fmx,SpatialFA%wmx))
      
      ! Store new coords
      SpatialFA%I = i ; SpatialFA%J = j

      ! Global Case
      IF( i .EQ. 0 .OR. j .EQ. 0 ) THEN
        
        ! (1) Mean spectrum
        rcode = nf_inq_varid( SpatialFA%ncid, 'mu_global', vid )
        rcode = nf_get_var_real( SpatialFA%ncid, vid, tmparr_r4(1:SpatialFA%wmx) )
        SpatialFA%Mu = REAL(tmparr_r4(1:SpatialFA%wmx) ,KIND=8)
        
        ! (2) Uncertainty in FA
        rcode = nf_inq_varid( SpatialFA%ncid, 'Psi_global', vid )
        rcode = nf_get_var_real( SpatialFA%ncid, vid, tmparr_r4(1:SpatialFA%wmx) )
        SpatialFA%Psi = REAL(tmparr_r4(1:SpatialFA%wmx) ,KIND=8)

        ! (3) Factor Loading Matrix
        rcode = nf_inq_varid( SpatialFA%ncid, 'W_global', vid )
        rcode = nf_get_var_real( SpatialFA%ncid, vid, tmpW )
        DO f=1,SpatialFA%fmx
          SpatialFA%W(:,f) = REAL(tmpW(f,:) ,KIND=8)
        ENDDO

      ! Regional Case
      ELSE
        
        ! (1) Mean spectrum
        rcode = nf_inq_varid( SpatialFA%ncid, 'mu', vid )
        rcode = nf_get_vara_real( SpatialFA%ncid, vid,             &
                                  (/j,i,1/),(/1,1,SpatialFA%wmx/), &
                                  tmparr_r4(1:SpatialFA%wmx)       )
        SpatialFA%Mu = REAL(tmparr_r4(1:SpatialFA%wmx) ,KIND=8)
        
        ! (2) Uncertainty in FA
        rcode = nf_inq_varid( SpatialFA%ncid, 'Psi', vid )
        rcode = nf_get_vara_real( SpatialFA%ncid, vid,             &
                                  (/j,i,1/),(/1,1,SpatialFA%wmx/), &
                                  tmparr_r4(1:SpatialFA%wmx)       )
        SpatialFA%Psi = REAL(tmparr_r4(1:SpatialFA%wmx) ,KIND=8)

        ! (3) Factor Loading Matrix
        rcode = nf_inq_varid( SpatialFA%ncid, 'W', vid )
        CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'W:inq_varid')
        rcode = nf_get_vara_real( SpatialFA%ncid, vid,                            &
                                  (/j,i,1,1/),(/1,1,SpatialFA%fmx,SpatialFA%wmx/),&
                                  tmpW                                            )
        CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'W:get_vara')
        DO f=1,SpatialFA%fmx
          SpatialFA%W(:,f) = REAL(tmpW(f,:) ,KIND=8)
        ENDDO

      ENDIF

      ! Basis spline expansion coefficients
      CALL SPLINE1(SpatialFA%Wvl,SpatialFA%Mu,SpatialFA%wmx,SpatialFA%MuSP)
      DO N=1,SpatialFA%fmx
        CALL SPLINE1(SpatialFA%Wvl,SpatialFA%W(:,N),SpatialFA%wmx,SpatialFA%WSP(:,N))
      ENDDO
      
      ! Deallocate temporary float arrays
      DEALLOCATE(tmparr_r4,tmpW)
      
    ENDIF
    
    ! Load MODIS Pixel
    CALL LoadBRDFPixel(Geolocation, SurfProf, SpatialFA, Error)
    
  END SUBROUTINE SampleSpatialFA

  SUBROUTINE UpdateFaGain(nWvl, nBand, nFa, Mu, W, Psi, S_brdf, RSR, G, MuBand, Sx_hat, Error)

    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    INTEGER,         INTENT(IN)    :: nWvl, nBand, nFa
    REAL(KIND=8),    INTENT(IN)    :: Mu(nWvl)
    REAL(KIND=8),    INTENT(IN)    :: W(nWvl,nFa)
    REAL(KIND=8),    INTENT(IN)    :: Psi(nWvl)
    REAL(KIND=8),    INTENT(IN)    :: S_brdf(nBand)
    REAL(KIND=8),    INTENT(IN)    :: RSR(nWvl,nBand)
    REAL(KIND=8),    INTENT(OUT)   :: G(nFa,nBand)
    REAL(KIND=8),    INTENT(OUT)   :: MuBand(nBand)
    REAL(KIND=8),    INTENT(OUT)   :: Sx_hat(nFa,nFa)
    TYPE(ErrorType), INTENT(INOUT) :: Error

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: Alpha, Beta
    REAL(KIND=8) :: MuVec(nWvl,1), MuBandVec(nBand,1)
    REAL(KIND=8) :: W_b(nBand,nFa)
    REAL(KIND=8) :: iPiv(nBand), WorkX(nBand)
    REAL(KIND=8) :: SoRSR(nWvl,nBand), So_b(nBand,nBand), So_b_inv(nBand,nBand)
    REAL(KIND=8) :: tmp1(nBand,nBand), tmp2(nFa,nBand)
    REAL(KIND=8) :: iPiv_f(nFa),WorkX_f(nFa)
    INTEGER      :: I, ReturnStatus

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'UpdateFaGain'
    
    ! =====================================================================
    ! UpdateFaGain starts here
    ! =====================================================================
    
    ! Compute the band version of the mean reflectance
    Alpha = 1.0d0 ; Beta = 0.0d0 ; MuVec(:,1) = Mu(:) ; MuBandVec(:,1) = 0.0d0
    CALL DGEMM('T','N',nBand,1,nWvl,Alpha,RSR,nWvl,MuVec,nWvl,Beta,MuBandVec,nBand)
    MuBand(:) = MuBandVec(:,1)
    
    ! Compute the Factors on the Band grid
    Alpha = 1.0d0 ; Beta = 0.0d0 ; W_b(:,:) = 0.0d0
    CALL DGEMM('T','N',nBand,nFA,nWvl,Alpha,RSR,nWvl,W,nWvl,Beta,W_b,nBand)

    ! Compute the observational error
    SoRSR(:,:) = 0.0d0
    DO I=1,nWvl
      SoRSR(I,:) = Psi(I)*RSR(I,:)
    ENDDO
    Alpha = 1.0d0 ; Beta = 1.0d0 ; So_b(:,:) = 0.0d0
    DO I=1,nBand
      So_b(I,I) = S_brdf(I)
    ENDDO
    CALL DGEMM('T','N',nBand,nBand,nWvl,Alpha,RSR,nWvl,SoRSR,nWvl,Beta,So_b,nBand)
    
    ! Invert WW^T+So
    Alpha = 1.0d0 ; Beta = 1.0d0 ; tmp1 = So_b
    CALL DGEMM('N','T',nBand,nBand,nFa,Alpha,W_b,nBand,W_b,nBand,Beta,tmp1,nBand)

    CALL DGETRF(nBand,nBand, tmp1, nBand, iPiv, ReturnStatus) ! LU Fac tmp1
    IF(ReturnStatus .NE. 0) THEN
      CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,&
                           Message_in='WW^T + So is singular!!'                 )
      RETURN
    ENDIF
    CALL DGETRI(nBand, tmp1, nBand, iPiv, WorkX, nBand, ReturnStatus)
    IF(ReturnStatus .NE. 0) THEN
      CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,&
                           Message_in='Inversion of WW^T + So failed!!'         )
      RETURN
    ENDIF

    ! Compute Gain
    Alpha = 1.0d0 ; Beta = 0.0d0 ; G(:,:) = 0.0d0
    CALL DGEMM('T','N',nFa,nBand,nBand,Alpha,W_b,nBand,tmp1,nBand,Beta,G,nFa)
    
    ! Invert observational error
    So_b_inv(:,:) = So_b(:,:)
    CALL DGETRF(nBand,nBand, So_b_inv, nBand, iPiv, ReturnStatus) ! LU Fac tmp1
    IF(ReturnStatus .NE. 0) THEN
      CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,&
                           Message_in='So_b is singular!!'                      )
      RETURN
    ENDIF
    CALL DGETRI(nBand, So_b_inv, nBand, iPiv, WorkX, nBand, ReturnStatus)
    IF(ReturnStatus .NE. 0) THEN
      CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,&
                           Message_in='Inversion of So_b failed!!'              )
      RETURN
    ENDIF

    ! Compute the uncertainty in the 
    Alpha = 1.0d0 ; Beta = 0.0d0 ; tmp2(:,:) = 0.0d0
    CALL DGEMM('T','N',nFa,nBand,nBand,Alpha,W_b,nBand,So_b_inv,nBand,Beta,tmp2,nFa)
    Alpha = 1.0d0 ; Beta = 1.0d0 ; Sx_hat(:,:) = 0.0d0
    DO I=1,nFa
      Sx_hat(I,I) = 1.0d0
    ENDDO
    
    CALL DGEMM('N','N',nFa,nFa,nBand,Alpha,tmp2,nFa,W_b,nBand,Beta,Sx_hat,nFa)
    CALL DGETRF(nFa,nFa, Sx_hat, nFa, iPiv_f, ReturnStatus) ! LU Fac tmp1
    IF(ReturnStatus .NE. 0) THEN
      CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,&
                           Message_in='Sx_hat is singular!!'                    )
      RETURN
    ENDIF
    CALL DGETRI(nFa, Sx_hat, nFa, iPiv_f, WorkX_f, nFa, ReturnStatus)
    IF(ReturnStatus .NE. 0) THEN
      CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,&
                           Message_in='Inversion of Sx_hat failed!!'            )
      RETURN
    ENDIF
    
  END SUBROUTINE UpdateFaGain

  SUBROUTINE LoadBRDFPixel( Geolocation, SurfProf, SpatialFA, Error )
    
    ! --------------------
    ! Subroutine Arguments
    ! --------------------
    TYPE(GeolocationType), INTENT(IN)      :: Geolocation
    TYPE(SurfProfType),    INTENT(IN)      :: SurfProf
    TYPE(SpatialFAType),   INTENT(INOUT)   :: SpatialFA
    TYPE(ErrorType),       INTENT(INOUT)   :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER      :: I,J,N,B, idt0, idtf
    REAL(KIND=8) :: max_xdiff, max_ydiff
    REAL(KIND=8) :: elons(2), elats(2), fracs(2), tmp
    REAL(KIND=8) :: DayOfYear
    REAL(KIND=8) :: DayDist(SpatialFA%BRDF_tmx+4)
    REAL(KIND=8) :: MinLat(SpatialFA%BRDF_tmx+4),MaxLat(SpatialFA%BRDF_tmx+4)
    INTEGER      :: idxs(2)
    LOGICAL      :: valid_brdf, IsNewTime
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadBRDFPixel'

    ! Masking depth values[m] for IGBP land types
    REAL(KIND=8), DIMENSION(17) :: ds_star = (/ 1e9,& !  0 - Water
                                               12.0,& !  1 - Evergreen Needleleaf
                                                8.0,& !  2 - Evergreen Broadleaf
                                                8.0,& !  3 - Deciduous Needleleaf
                                                8.0,& !  4 - Deciduous Broadleaf
                                                8.0,& !  5 - Mixed Forest (Use evergreen broadleaf factor)
                                                1.0,& !  6 - Closed shrublands
                                                1.0,& !  7 - Open shrublands
                                               0.15,& !  8 - Woody savannas
                                               0.15,& !  9 - Savannas
                                                0.1,& ! 10 - Grasslands
                                               0.01,& ! 11 - Permanent Wetlands
                                                0.1,& ! 12 - Croplands
                                                4.0,& ! 13 - Urban and buit up
                                                0.1,& ! 14 - Cropland/ Natural vegetation mosaic
                                               0.01,& ! 15 - Snow and Ice
                                               0.01 /)! 16 - Barren or Sparsely vegetated
    
    ! =====================================================================
    ! LoadBRDFPIxel starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! ------------------------------------------
    ! Snow Spectrum
    ! ------------------------------------------
    
    ! Get Snow spectrum
    IF( SurfProf%SnowFraction .GT. TINY(SurfProf%SnowFraction) ) THEN

      ! Compute effective snow fraction (accounting for masking depths of different land types)
      DO i=1,17
        SpatialFA%EffectiveSnowFraction = SpatialFA%EffectiveSnowFraction                   &
        + SurfProf%LandCoverFraction(i)*( 1.0d0-exp(-1.0*SurfProf%SnowDepth/ds_star(i)))
      ENDDO
      SpatialFA%EffectiveSnowFraction = SpatialFA%EffectiveSnowFraction * SurfProf%SnowFraction

      ! Compute Snowage weight
      SpatialFA%SnowAgeWeight = EXP( -0.2d0*SurfProf%SnowAge )

      ! Compute snow spectrum
      SpatialFA%SnowAlbedo = SpatialFA%SnowAgeWeight*SpatialFA%SnowFineAlbedo &
                           + (1.0d0-SpatialFA%SnowAgeWeight)*SpatialFA%SnowCoarseAlbedo
    ENDIF
    
    ! Sea Ice fraction
    SpatialFA%SeaIceFraction = SurfProf%SeaIceFraction
    
    ! Check if we need to load the pixel
    max_xdiff = MAXVAL(ABS(Geolocation%CornerLongitudes - SpatialFA%PixelCornerLongitudes))
    max_ydiff = MAXVAL(ABS(Geolocation%CornerLatitudes - SpatialFA%PixelCornerLatitudes))

    ! ------------------------------------------
    ! Load the pixel
    ! ------------------------------------------

    ! Check if its a new time
    IF(SpatialFA%IsClimatology) THEN
      IsNewTime = ABS(SpatialFA%DayOfYear - Geolocation%Time%DayOfYear) .GT. TINY(0.0d0)
    ELSE
      IsNewTime = ABS(SpatialFA%TauTime-Geolocation%Time%Tau) .GT. TINY(0.0d0)
    ENDIF
    
    ! If the pixel coordinates are significantly different reload pixel
    IF( MAX(max_xdiff,max_ydiff) .GT. TINY(max_xdiff) .OR. IsNewTime ) THEN

      ! Get the longitude/latitude bounds of the pixel
      elons(1) = minval(Geolocation%CornerLongitudes)
      elons(2) = maxval(Geolocation%CornerLongitudes)
      elats(1) = minval(Geolocation%CornerLatitudes)
      elats(2) = maxval(Geolocation%CornerLatitudes)

      ! If pixel overlaps +/- 180 lon, then make points > -180 < 180
      ! IF(elons(2)-elons(1) .GT. 180.0d0) THEN
      !   tmp = elons(2)-180.0d0
      !   elons(2) = elons(1) ; elons(1) = tmp
      ! ENDIF

      ! -----------------------------------------------------------------
      ! Determine which MODIS files to read and their weighting fractions
      ! -----------------------------------------------------------------
      
      ! Archive DOY 
      SpatialFA%DayOfYear = Geolocation%Time%DayOfYear
      SpatialFA%TauTime   = Geolocation%Time%Tau

      ! Archive new coordinates
      SpatialFA%PixelCornerLongitudes = Geolocation%CornerLongitudes
      SpatialFA%PixelCornerLatitudes  = Geolocation%CornerLatitudes
      SpatialFA%PixelLongitude        = Geolocation%Longitude
      SpatialFA%PixelLatitude         = Geolocation%Latitude
      
      ! Find the nearest time
      IF(SpatialFA%IsClimatology) THEN

        ! Day of year as double
        DayOfYear = REAL(SpatialFA%DayOfYear,KIND=8)

        ! Get indices and time-weighted fractions
        idxs(1) = MINVAL(MAXLOC(SpatialFA%BRDF_DayDbl, MASK=(SpatialFA%BRDF_DayDbl < DayOfYear )))
        idxs(2) = idxs(1) + 1
        fracs(1) = (DayOfYear - SpatialFA%BRDF_DayDbl(idxs(1))) &
                 / (SpatialFA%BRDF_DayDbl(idxs(2)) - SpatialFA%BRDF_DayDbl(idxs(1)))
        fracs(2) = 1.0 - fracs(2)
        
        ! Load indices
        idt0 = SpatialFA%BRDF_DayIdx(idxs(1)) ; idtf = SpatialFA%BRDF_DayIdx(idxs(2)) 

        ! Check that there are observations for given MODIS fraction
        ! ----------------------------------------------------------
        IF( elats(2) > SpatialFA%BRDF_MaxLat(idt0) .OR. &
            elats(2) > SpatialFA%BRDF_MaxLat(idtf) .OR. &
            elats(1) < SpatialFA%BRDF_MinLat(idt0) .OR. &
            elats(1) < SpatialFA%BRDF_MinLat(idtf)      ) THEN

            ! Compute day distances
            DO i=1,SpatialFA%BRDF_tmx+4
              DayDist(i) = MIN( ABS(SpatialFA%BRDF_DayDbl(i)-DayOfYear),        &
                                ABS(365.0d0-DayOfYear)+SpatialFA%BRDF_DayDbl(i) )
              MaxLat(i) = SpatialFA%BRDF_MaxLat(SpatialFA%BRDF_DayIdx(i))
              MinLat(i) = SpatialFA%BRDF_MinLat(SpatialFA%BRDF_DayIdx(i))
            ENDDO

            ! Find the nearest valid latitude
            idxs(1) = MINVAL( MINLOC( DayDist, MASK=( MaxLat > elats(2) .AND. MinLat < elats(1)  ) ) )
            idt0 = SpatialFA%BRDF_DayIdx(idxs(1)) ; idtf = idt0
            fracs(1) = 1.0d0
            fracs(2) = 0.0d0
            
            ! Check if there are any valid indices
            IF(idxs(1) .LE. 0) THEN

              ! Raise Pixel Error
              CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,      &
                                   Message_in='BRDF Climatology does not cover pixel location')

              ! Return to calling program
              RETURN

            ENDIF
        ENDIF

      ! Use Tau time
      ELSE

        IF( Geolocation%Time%Tau .LE. SpatialFA%BRDF_Tau(1) ) THEN
          idt0     = 1     ; idtf     = 1
          fracs(1) = 1.0d0 ; fracs(2) = 0.0d0

        ELSEIF( Geolocation%Time%Tau .GE. SpatialFA%BRDF_Tau(SpatialFA%BRDF_tmx) ) THEN
          idt0     = SpatialFA%BRDF_tmx ; idtf     = idt0
          fracs(1) = 1.0d0              ; fracs(2) = 0.0d0

        ELSE
          idt0 = MINVAL(MAXLOC(SpatialFA%BRDF_Tau, MASK=(SpatialFA%BRDF_Tau < Geolocation%Time%Tau )))
          idtf = idt0 + 1
          fracs(1) = (Geolocation%Time%Tau - SpatialFA%BRDF_Tau(idt0)   ) &
                   / (SpatialFA%BRDF_Tau(idtf) - SpatialFA%BRDF_Tau(idt0))
          fracs(2) = 1.0 - fracs(1)
        ENDIF

        ! Check that there are observations for given MODIS fraction
        ! ----------------------------------------------------------
        IF( elats(2) > SpatialFA%BRDF_MaxLat(idt0) .OR. &
            elats(2) > SpatialFA%BRDF_MaxLat(idtf) .OR. &
            elats(1) < SpatialFA%BRDF_MinLat(idt0) .OR. &
            elats(1) < SpatialFA%BRDF_MinLat(idtf)      ) THEN


            ! Find closest valid index
            idt0 = MINVAL( MINLOC( ABS(SpatialFA%BRDF_Tau-Geolocation%Time%Tau), &
                           MASK=( SpatialFA%BRDF_MaxLat > elats(2) .AND.         & 
                                  SpatialFA%BRDF_MinLat < elats(1)     )     )   )
            idtf = idt0 ; fracs(1) = 1.0d0 ; fracs(2) = 0.0d0

            ! Check if there are any valid indices
            IF(idxs(1) .LE. 0) THEN

              ! Raise Pixel Error
              CALL RaisePixelError(Error, ErrorCode_OptProp, ModuleName, SubroutineName,      &
                                   Message_in='BRDF Climatology does not cover pixel location')

              ! Return to calling program
              RETURN

            ENDIF
            
        ENDIF

      ENDIF
      
      ! Reload BRDF
      CALL ReadClimatologyBRDF(elons,elats,(/idt0,idtf/),fracs,SpatialFA,Error)
      
      ! Now Compute the gain Matrix
      CALL UpdateFaGain(SpatialFA%wmx,SpatialFA%bmx,SpatialFA%fmx,    &
                        SpatialFA%Mu,SpatialFA%W,SpatialFA%Psi,       &
                        SpatialFA%BandUncertainty,                    &
                        SpatialFA%RSR,SpatialFA%G,SpatialFA%MuBands,  &
                        SpatialFA%FactorUncertainty, Error            )

      ! Compute Kernel Amplitudes
      IF( SpatialFA%DoIsotropic ) THEN
        CALL compute_RTLS_albedo(Geolocation, SpatialFA, Error)
      ELSE
        CALL compute_RTLS_brdf(Geolocation, SpatialFA, Error)
      ENDIF

    ENDIF
    

  END SUBROUTINE LoadBRDFPixel
  
  SUBROUTINE ReadClimatologyBRDF(elons, elats, idt, dayfrac, SpatialFA, Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    REAL(KIND=8),        INTENT(IN)    :: elons(2)
    REAL(KIND=8),        INTENT(IN)    :: elats(2)
    INTEGER,             INTENT(IN)    :: idt(2)
    REAL(KIND=8),        INTENT(IN)    :: dayfrac(2)
    TYPE(SpatialFAType), INTENT(INOUT) :: SpatialFA
    TYPE(ErrorType),     INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER :: ix0_m, ixf_m, ix0_p, ixf_p, nx_m, nx_p, nx, nx_aloc, ny_aloc
    INTEGER :: iy0, iyf, ny, ct, i, j, n, b, t, nt
    REAL(KIND=8) :: LonCorner(4),LatCorner(4),TimeFrac(2)
    LOGICAL :: valid_brdf

    INTEGER(KIND=2), ALLOCATABLE :: BRDF_int(:,:,:,:)
    LOGICAL,         ALLOCATABLE :: PixIsLand(:,:)
    REAL(KIND=8),    ALLOCATABLE :: SubPixLon_NoAdjust(:)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ReadClimatologyBRDF'
    
    ! =====================================================================
    ! ReadClimatologyBRDF starts here
    ! =====================================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    ! Copy Corners
    LonCorner = SpatialFA%PixelCornerLongitudes
    LatCorner = SpatialFA%PixelCornerLatitudes
    TimeFrac  = dayfrac


    ! This checks for longitudes overlapping +/- 180
    IF(elons(2)-elons(1) .GT. 180.0d0) THEN

      ! Load block between [elons(1)+360.0,180]
      ix0_m = MINVAL(MAXLOC(SpatialFA%BRDF_Lon, MASK=(SpatialFA%BRDF_lon .LT. elons(1)+360.0d0 )))
      ixf_m = SpatialFA%BRDF_imx
      nx_m  = ixf_m - ix0_m + 1

      ! Block between [elons(2),360]
      ix0_p = MINVAL(MINLOC(SpatialFA%BRDF_lon,MASK=(SpatialFA%BRDF_lon .GE. elons(2))))
      ixf_p = SpatialFA%BRDF_imx
      nx_p  = ixf_p - ix0_p + 1

      ! If we are overlapping then convert corner longitudes > 0 to < -180
      DO n=1,4
        IF( LonCorner(n) .GT. 0.0 ) LonCorner(n) = LonCorner(n) - 360.0d0
      ENDDO
    
    ELSE

      ! No points < -180
      nx_m = 0

      ! Region to load
      ix0_p = MINVAL(MAXLOC(SpatialFA%BRDF_Lon, MASK=(SpatialFA%BRDF_lon .LE. elons(1))) )
      ixf_p = MAXVAL(MINLOC(SpatialFA%BRDF_lon, MASK=(SpatialFA%BRDF_lon .GE. elons(2))) )
      nx_p  = ixf_p - ix0_p + 1

    ENDIF
    
    ! Total x points
    nx = nx_p + nx_m

    ! Get y limits
    iy0 = MINVAL(MAXLOC(SpatialFA%BRDF_Lat, MASK=(SpatialFA%BRDF_lat .LE. elats(1))) )
    iyf = MAXVAL(MINLOC(SpatialFA%BRDF_lat, MASK=(SpatialFA%BRDF_lat .GE. elats(2))) )
    ny  = iyf - iy0 + 1

    ! Allocate temporary read arrays
    ALLOCATE(BRDF_int(nx,ny,SpatialFA%bmx,3)) ; ALLOCATE(PixIsLand(nx,ny)) 

    ! Expand array dimension if necessary
    IF(nx > SpatialFA%SubPixImx .OR. ny > SpatialFA%SubPixJmx) THEN
        
      ! Dimensions for array reallocation
      nx_aloc = max( nx, SpatialFA%SubPixImx )
      ny_aloc = max( ny, SpatialFA%SubPixJmx )
      
      IF(ALLOCATED(SpatialFA%SubPixBRDF))     DEALLOCATE(SpatialFA%SubPixBRDF)
      IF(ALLOCATED(SpatialFA%SubPixLandFrac)) DEALLOCATE(SpatialFA%SubPixLandFrac)
      IF(ALLOCATED(SpatialFA%SubPixLon))      DEALLOCATE(SpatialFA%SubPixLon)
      IF(ALLOCATED(SpatialFA%SubPixLat))      DEALLOCATE(SpatialFA%SubPixLat)
      
      ALLOCATE(SpatialFA%SubPixBRDF(nx_aloc, ny_aloc, SpatialFA%bmx, 3))
      ALLOCATE(SpatialFA%SubPixLandFrac(nx_aloc, ny_aloc))
      ALLOCATE(SpatialFA%SubPixLon(nx_aloc) )
      ALLOCATE(SpatialFA%SubPixLat(ny_aloc) )
      
      SpatialFA%SubPixImx = nx
      SpatialFA%SubPixJmx = ny
      
    ENDIF

    ALLOCATE(SubPixLon_NoAdjust(nx))

    ! Compute lon/lat coordinates
    ct = 0
    DO i=1,nx_m
      ct = ct + 1
      SpatialFA%SubPixLon(ct) = SpatialFA%BRDF_Lon(ix0_m+i-1)-360.0d0
      SubPixLon_NoAdjust(ct) = SpatialFA%BRDF_Lon(ix0_m+i-1)
    ENDDO
    DO i=1,nx_p
      ct = ct + 1
      SpatialFA%SubPixLon(ct) = SpatialFA%BRDF_Lon(ix0_p+i-1)
      SubPixLon_NoAdjust(ct) = SpatialFA%BRDF_Lon(ix0_p+i-1)
    ENDDO
    DO j=1,ny
      SpatialFA%SubPixLat(j) = SpatialFA%BRDF_Lat(iy0+j-1) 
    ENDDO
    
    ! Determine Number of files to load
    nt = 2 
    IF(idt(1) .EQ. idt(2))  THEN
      nt = 1
      TimeFrac(1) = 1.0d0
    ENDIF
    
    ! Now we need to load the fields
    SpatialFA%SubPixBRDF(:,:,:,:) = 0.0d0 ; SpatialFA%SubPixLandFrac(:,:) = 1.0d0
    
    ! --------------------------------------------
    ! Load subregion for pixels < -180 if required
    ! --------------------------------------------
    IF(nx_m .GT. 0) THEN
      
      DO t=1,nt
        CALL LoadSubRegionBRDF(SpatialFA,ix0_m,iy0,nx_m,ny,idt(t),             &
                               BRDF_int(1:nx_m,1:ny,:,:),PixIsLand(1:nx_m,:),Error)

        ! Set BRDF For sub region
        SpatialFA%SubPixBRDF(1:nx_m,1:ny,:,:) = SpatialFA%SubPixBRDF(1:nx_m,1:ny,:,:) &
           + REAL(BRDF_int(1:nx_m,1:ny,:,:),KIND=8)*1.0d-3*TimeFrac(t)

        ! Check for ocean
        DO i=1,nx_m
        DO j=1,ny
          IF(.NOT. PixIsLand(i,j)) THEN
            SpatialFA%SubPixLandFrac(i,j) = 0.0d0
          ENDIF
        ENDDO
        ENDDO
        
      ENDDO

    ENDIF
    
    ! --------------------------------------------
    ! Load subregion for pixels > -180 if required
    ! --------------------------------------------
    IF(nx_p .GT. 0) THEN

      DO t=1,nt
        CALL LoadSubRegionBRDF(SpatialFA,ix0_p,iy0,nx_p,ny,idt(t),                &
                            BRDF_int(1+nx_m:nx,:,:,:),PixIsLand(1+nx_m:nx,:),Error)

        ! Set BRDF For sub region
        SpatialFA%SubPixBRDF(1+nx_m:nx,1:ny,:,:) = SpatialFA%SubPixBRDF(1+nx_m:nx,1:ny,:,:) &
           + REAL(BRDF_int(1+nx_m:nx,1:ny,:,:),KIND=8)*1.0d-3*TimeFrac(t)

        ! Check for ocean
        DO i=1+nx_m,nx
        DO j=1,ny
          IF(.NOT. PixIsLand(i,j)) THEN
            SpatialFA%SubPixLandFrac(i,j) = 0.0d0
          ENDIF
        ENDDO
        ENDDO
        
      ENDDO

    ENDIF
    
    ! Zero pixel averaged quantities
    N = 0
    SpatialFA%PixelBRDF(:,:) = 0.0d0
    SpatialFA%LandFraction = 0.0d0
      
    ! Average BRDF Within Pixel
    DO I=1,SpatialFA%SubPixImx
    DO J=1,SpatialFA%SubPixJmx
      
      ! Add points within pixel
      IF( PNPOLY_r8(4,LonCorner,LatCorner, &
                    SpatialFA%SubPixLon(I),&
                    SpatialFA%SubPixLat(J))) THEN
          
        ! Increment count
        n = n + 1
          
        ! Check isotropic kernel value for each band
        valid_brdf = .TRUE.
        DO b=1,SpatialFA%bmx
          IF( SpatialFA%SubPixBRDF(i,j,b,1)  .LT. 0.0 .OR. &
               SpatialFA%SubPixBRDF(i,j,b,1) .GT. 10.0     ) THEN
            valid_brdf = .FALSE.
          ENDIF
        ENDDO

        ! Add to average  
        IF( SpatialFA%SubPixLandFrac(i,j) > 0.0d0 .AND. valid_brdf ) THEN ! Its either 0 or 1
          SpatialFA%PixelBRDF = SpatialFA%PixelBRDF + SpatialFA%SubPixBRDF(i, j, :, :)
          SpatialFA%LandFraction = SpatialFA%LandFraction + 1.0d0
        ENDIF
          
      ENDIF
        
    ENDDO
    ENDDO
    
    ! Complete Average
    IF( n .GT. 0 ) THEN

      SpatialFA%LandFraction = SpatialFA%LandFraction / REAL(n,KIND=8)
      SpatialFA%PixelBRDF = SpatialFA%PixelBRDF / REAL(n,KIND=8)
        
    ! Nearest neighbour sample if no points within pixel
    ELSE

      ! get indices
      i = MINVAL( MINLOC( ABS(SubPixLon_NoAdjust(1:nx)-SpatialFA%PixelLongitude) ) )
      j = MINVAL( MINLOC( ABS(SpatialFA%SubPixLat(1:ny)-SpatialFA%PixelLatitude ) ) )
      
      ! Check valid
      valid_brdf = .TRUE.
      DO b=1,SpatialFA%bmx
        IF( SpatialFA%SubPixBRDF(i,j,b,1) < 0.0 .OR. &
            SpatialFA%SubPixBRDF(i,j,b,1) > 10.0     ) THEN
          valid_brdf = .FALSE.
        ENDIF
      ENDDO
      
      IF( SpatialFA%SubPixLandFrac(i,j) > 0.0 .AND. valid_brdf ) THEN
        SpatialFA%LandFraction = 1.0d0
        SpatialFA%PixelBRDF    = SpatialFA%SubPixBRDF(i,j,:,:)
      ENDIF

   ENDIF
    
    DEALLOCATE(SubPixLon_NoAdjust)

  END SUBROUTINE ReadClimatologyBRDF
  
  SUBROUTINE LoadSubRegionBRDF(SpatialFA,ix0,iy0,nx,ny,it0,BRDF,PixIsLand,Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(SpatialFAType), INTENT(INOUT) :: SpatialFA
    INTEGER,             INTENT(IN)    :: ix0, iy0, nx, ny, it0
    INTEGER(KIND=2),     INTENT(OUT)   :: BRDF(nx,ny,SpatialFA%bmx,3)
    LOGICAL,             INTENT(OUT)   :: PixIsLand(nx,ny)
    TYPE(ErrorType),     INTENT(INOUT) :: Error
    
    ! ---------------
    ! local variables
    ! ---------------
    INTEGER            :: i, j, b, k, rcode, ncid, vid
    CHARACTER(LEN=100) :: bstr,tstr
    INTEGER            :: vmin_vec(2), nv_vec(2)
    CHARACTER(LEN=5), PARAMETER :: varname(3) = (/'f_iso','f_vol','f_geo'/)
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'LoadSubRegionBRDF'
    
    ! =====================================================================
    ! LoadSubRegionBRDF starts here
    ! =====================================================================

    ! Set load indices
    
    ! Time string
    WRITE(tstr,'(I100)') it0

    ! Set starts and strides
    vmin_vec(1) = ix0 ; vmin_vec(2) = iy0 
    nv_vec(1)   = nx  ; nv_vec(2)   = ny

    ! Initialize flag for land
    PixIsLand(:,:) = .TRUE.

    ! Open file and attach file indices if needed
    IF(SpatialFA%current_tidx .NE. it0) THEN
      
      ! Shut previous files
      IF(SpatialFA%current_tidx .GT. 0) THEN
        DO b=1,SpatialFA%bmx
          rcode = nf_close(SpatialFA%obs_ncid(b))
          CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'close')
        ENDDO
      ENDIF

      ! Open datasets for it0
      DO b=1,SpatialFA%bmx

        ! Band Index string
        WRITE(bstr,'(I100)') b
        
        ! Open File
        rcode = nf_open(TRIM(ADJUSTL(SpatialFA%RootBRDFDir)) // '/BRDF_Band'  //      &
                        TRIM(ADJUSTL(bstr)) // '_Time' //TRIM(ADJUSTL(tstr)) // '.nc',&
                        nf_Share,SpatialFA%obs_ncid(b)                                )
        CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,'open')

        ! Get variable indices
        DO k=1,3
          rcode = nf_inq_varid(SpatialFA%obs_ncid(b) , varname(k), SpatialFA%obs_vid(b,k))
          CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
                                      varname(k)//':inq_varid')
        ENDDO
        
      ENDDO

      ! Update the current index
      SpatialFA%current_tidx = it0

    ENDIF

    ! Read Data
    DO b=1,SpatialFA%bmx
      DO k=1,3

        ! Load Field
        rcode = nf_get_vara_int2(SpatialFA%obs_ncid(b), SpatialFA%obs_vid(b,k),&
                                 vmin_vec, nv_vec, BRDF(:,:,b,k)               )
        CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName,&
                                    varname(k)//':get_vara'              )
        
        ! Set Logical for water
        DO i=1,nx
        DO j=1,ny
          IF(BRDF(i,j,b,k) .LT. 0) THEN
            PixIsLand(i,j) = .FALSE.
          ENDIF
        ENDDO
        ENDDO

      ENDDO

    ENDDO
    
  END SUBROUTINE LoadSubRegionBRDF



  LOGICAL FUNCTION PNPOLY_r8( nvert, vertx, verty, testx, testy)
    
    IMPLICIT NONE
    
    ! Return variable
    ! ---------------
    !LOGICAL :: PNPOLY ! True if point is in polygon
    
    ! Input
    INTEGER :: nvert ! # Vertices in polygon
    REAL(KIND=8), DIMENSION(nvert) :: vertx ! X coordinates of the polygon vertices
    REAL(KIND=8), DIMENSION(nvert) :: verty ! Y coordinates of the polygon vertices
    REAL(KIND=8)                   :: testx ! X coordinate of test point
    REAL(KIND=8)                   :: testy ! Y coordinate of test point
    
    ! Local variables
    INTEGER :: i, j
    
    ! =====================================================================
    ! PNPOLY_r8 starts here
    ! =====================================================================

    ! Initialize PNPOLY
    PNPOLY_r8 = .FALSE.
    
    ! Initialize j
    j = nvert
    
    ! Perform crossings test (Jordan curve thm.) to see if point is in polygon
    DO i=1,nvert
      
      IF ( ((verty(i) .gt. testy) .neqv. (verty(j) .gt. testy)) .and. &
          (testx .lt. (vertx(j)-vertx(i)) * &
          (testy-verty(i)) / (verty(j)-verty(i)) + vertx(i)) ) THEN
         
         PNPOLY_r8 = .NOT. PNPOLY_r8
         
      ENDIF
      
      j = i
      
    ENDDO
    
    RETURN
    
  END FUNCTION PNPOLY_r8

  SUBROUTINE BRDF_kernels(sza, aza, vza, kvol, kgeo)
  
    ! --------------------
    ! subroutine arguments
    ! --------------------
    REAL (KIND=8), INTENT(IN)  :: sza, aza, vza
    REAL (KIND=8), INTENT(OUT) :: kvol, kgeo

    ! ---------------
    ! local variables
    ! ---------------
    !Ratio of elevation to height of (spherical) tree crowns
    REAL (KIND=8), PARAMETER   :: hb = 2.0d0 
    
    REAL (KIND=8) :: sza1, aza1, vza1, cossza, cosaza, cosvza, &
        sinsza, sinaza, sinvza, tansza, tanvza, secsza, secvza, &
        cosxi, xi, sinxi, d, o, t, cost
    
    ! =====================================================================
    ! BRDF_kernels starts here
    ! =====================================================================

    sza1 = sza * deg2rad
    !aza1 = aza * deg2rad  
    aza1 = (180.0 - aza) * deg2rad  ! VLIDORT convention is different from MODIS BRDF convention 
    vza1 = vza * deg2rad

    cossza = COS(sza1)
    cosaza = COS(aza1)
    cosvza = COS(vza1)

    sinsza = SIN(sza1)
    sinaza = SIN(aza1)
    sinvza = SIN(vza1)

    tansza = TAN(sza1)
    tanvza = TAN(vza1)

    secsza = 1.d0 / cossza
    secvza = 1.d0 / cosvza
    
    cosxi = cossza * cosvza + sinsza * sinvza * cosaza 

    xi = ACOS(cosxi)

    sinxi = SIN(xi)
    
    ! RossThick
    kvol = ( ( (Constants_pi /2.0d0 - xi) * cosxi + sinxi ) / (cossza + cosvza)) &
         - (Constants_pi / 4.0d0)

    d = SQRT( tansza**2 + tanvza**2 - 2 * tansza * tanvza * cosaza)

    cost = hb * SQRT(d**2 + (tansza * tanvza * sinaza)**2) / (secsza + secvza)

    IF (cost > 1) THEN
      cost = 1.0d0
    ENDIF

    t = ACOS(cost)

    o = (1.0d0 / Constants_pi) * (t - SIN(t) * cost) * (secsza + secvza)
    
    ! LiSparse
    kgeo = o - secsza - secvza + (1 + cosxi) * secsza * secvza / 2.d0

  END SUBROUTINE BRDF_kernels

  SUBROUTINE compute_RTLS_albedo(Geolocation, SpatialFA, Error)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(GeolocationType), INTENT(IN)     :: Geolocation
    TYPE(SpatialFAType),   INTENT(INOUT)  :: SpatialFA
    TYPE(ErrorType),       INTENT(INOUT)  :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER, PARAMETER       :: nband=4, nkernel=3, neof = 4
    INTEGER, PARAMETER       :: nsza0 = 10, nsnowpos = 480
    REAL (KIND=8), DIMENSION(nkernel), PARAMETER    :: wsg = (/1.0d0, 0.18984d0, -1.377622d0/)
    REAL (KIND=8), DIMENSION(3, nkernel), PARAMETER :: bsg = reshape ( (/ &
          1.0d0,       0.d0,        0.d0,      &
         -0.007574d0, -0.070987d0, 0.307588d0, &
         -1.284909d0, -0.166314d0, 0.041840d0 /), (/3,nkernel/) )
    
    INTEGER                                   :: i
    REAL(KIND=8)                              :: kgeo, kvol, sza2, sza3, wtemp,&
                                                 secsza, szafrac
    REAL(KIND=8), DIMENSION(SpatialFA%bmx)    :: mBRDF
    REAL(KIND=8), DIMENSION(nband)            ::  dirfrac, fracslp
    
    REAL(KIND=8), DIMENSION(nsza0), SAVE        :: secsza0
    REAL(KIND=8), DIMENSION(nband, nsza0), SAVE :: dirfrac0, fracslp0
    
    LOGICAL, SAVE :: FIRST = .TRUE.
    
    ! ==============================================================
    ! compute_RTLS_albedo starts here
    ! ==============================================================
    
    ! Check error status before computation
    IF(CheckError(Error)) RETURN

    
    IF( FIRST ) THEN
    
      ! Read LUT of ratio of direct irradiance to total downwelling irradiance
      OPEN(UNIT=eofunit, file=trim(SpatialFA%RootDataDir)//"/BRDF_EOF/AlbSpec/" // &
           "directflux_to_totdnflux_modis_band1-4.txt", status='OLD')
      READ (eofunit, *)
      READ (eofunit, *) secsza0
      secsza0 = 1.0/COS(secsza0 * deg2rad)
      DO i = 1, 4
        READ(eofunit, *)
      ENDDO
      DO i = 1, nband
        READ(eofunit, *) wtemp, dirfrac0(i, 1:nsza0), fracslp0(i, 1:nsza0)
      ENDDO
      CLOSE(UNIT = eofunit)
      
    ENDIF
    
    ! Reflectance/Albedo
    mBRDF(:) = 0.0
    IF (SpatialFA%LandFraction .GT. 0.0) THEN ! At least partial land
      
      IF (SpatialFA%WhichAlbedo == 3) THEN
        
        ! #########################################################
        ! Here the weightings need to be interpolated to the given 
        ! bands 
        ! calculate_blueskyalb needs to be updated to account for 
        ! flexible band implementation
        ! ########################################################

        STOP 'Blue sky albedo not yet implemented'
        ! secsza = 1.0/COS(Geolocation%SZA * deg2rad)
          
        !   IF (secsza <= secsza0(1)) THEN
        !     dirfrac = dirfrac0(:, 1)
        !     fracslp = fracslp0(:, 1)
        !   ELSE IF (secsza >= secsza0(nsza0)) THEN
        !     dirfrac = dirfrac0(:, nsza0)
        !     fracslp = fracslp0(:, nsza0)
        !   ELSE
        !     DO i = 2, nsza0
        !         IF (secsza < secsza0(i)) EXIT
        !     ENDDO
        !     szafrac = (secsza - secsza0(i-1)) / (secsza0(i) - secsza0(i-1))
        !     dirfrac = dirfrac0(:, i-1) * (1.0 - szafrac) + dirfrac0(:, i) * szafrac
        !     fracslp = fracslp0(:, i-1) * (1.0 - szafrac) + fracslp0(:, i) * szafrac
        !   ENDIF
      ENDIF
      
      IF (SpatialFA%WhichAlbedo == 0) THEN   ! Directional albedo
          !Get BRDF kernels at MODIS points for viewing geometry
          CALL BRDF_kernels(Geolocation%SZA, Geolocation%AZA, Geolocation%VZA, kvol, kgeo)
      ELSE IF (SpatialFA%WhichAlbedo == 1 .OR. SpatialFA%WhichAlbedo == 3) THEN  ! For Black/Blue albedo
          sza2 = (Geolocation%SZA * deg2rad) ** 2
          sza3 = sza2 * Geolocation%SZA * deg2rad
          kvol = bsg(1, 2) + sza2 * bsg(2, 2) + sza3 * bsg(3, 2)
          kgeo = bsg(1, 3) + sza2 * bsg(2, 3) + sza3 * bsg(3, 3)
      ELSE IF (SpatialFA%WhichAlbedo == 2) THEN  ! White albedo
          kvol = wsg(2); kgeo = wsg(3)
      ENDIF
      
      !Calculate BRDF/Directional Reflectance at MODIS points
      IF (SpatialFA%WhichAlbedo <= 2) THEN
          DO i = 1, SpatialFA%bmx
            mBRDF(i) = SpatialFA%PixelBRDF(i,1)      &
                     + SpatialFA%PixelBRDF(i,2)*kvol &
                     + SpatialFA%PixelBRDF(i,3)*kgeo
          ENDDO
      ELSE IF (SpatialFA%WhichAlbedo == 3) THEN
          STOP 'Blue Sky albedo needs re-implementation in spatial FA'
          ! CALL calculate_blueskyalb(nband, &
          !                           SpatialFA%PixelBRDF(1:SpatialFA%bmx,1),&
          !                           SpatialFA%PixelBRDF(1:SpatialFA%bmx,2),&
          !                           SpatialFA%PixelBRDF(1:SpatialFA%bmx,3),&
          !                           kvol, kgeo, wsg(2), wsg(3), &
          !                           dirfrac, fracslp, mBRDF)
      ENDIF
      
      ! STOP 'Testing SpatialFA'
      ! Compute the land-only Kernel amplitudes
      CALL ComputeFACoeff(mBRDF, SpatialFA%MuBands, SpatialFA%Mu, SpatialFA%G, &
                          SpatialFA%W, SpatialFA%FitCoeff(:,1),                &
                          SpatialFA%KernelAmplitudes(:,1),                     &
                          SpatialFA%bmx, SpatialFA%wmx, SpatialFA%fmx          )

    ELSE
    
      ! Set land amplitudes to zero
      SpatialFA%KernelAmplitudes(:,1) = 0.0d0
      
    ENDIF
    
    ! Weight the land spectrum by ocean
    SpatialFA%KernelAmplitudes(:,1) = SpatialFA%LandFraction*(                 &
       SpatialFA%KernelAmplitudes(:,1)*(1.0d0-SpatialFA%EffectiveSnowFraction) &
     + SpatialFA%SnowAlbedo(:)*SpatialFA%EffectiveSnowFraction  )              &
                                    + (1.0d0-SpatialFA%LandFraction)*(         &
       SpatialFA%WaterAlbedo(:)*(1.0d0-SpatialFA%SeaIceFraction)               &
     + SpatialFA%SnowAlbedo(:)*SpatialFA%SeaIceFraction              )
    
    ! Compute basis spline expansion coefficients
    CALL SPLINE1(SpatialFA%Wvl,SpatialFA%KernelAmplitudes(:,1),&
                 SpatialFA%wmx,SpatialFA%KernelAmplitudesSP(:,1))
    
    
  END SUBROUTINE compute_RTLS_albedo

  SUBROUTINE compute_RTLS_brdf(Geolocation, SpatialFA, Error)

    ! --------------------
    ! subroutine arguments
    ! --------------------
    TYPE(GeolocationType), INTENT(IN)     :: Geolocation
    TYPE(SpatialFAType),   INTENT(INOUT)  :: SpatialFA
    TYPE(ErrorType),       INTENT(INOUT)  :: Error

    ! ---------------
    ! local variables
    ! ---------------
    REAL(KIND=8) :: f_land, f_landice, f_sea, f_seaice
    INTEGER      :: n

    ! =====================================================================
    ! compute_RTLS_brdf starts here
    ! =====================================================================

    ! Check error status before computation
    IF(CheckError(Error)) RETURN
    
    ! Compute wavelength dependent kernel amplitudes ComputeFACoeff(modis_refl, MuMODIS, Mu, G, W, coeff, hyprefl, wmx_modis, wmx, fmx)
    DO n=1,3
      CALL ComputeFACoeff(SpatialFA%PixelBRDF(:,n),                        &
                          SpatialFA%MuBands, SpatialFA%Mu, SpatialFA%G,    &
                          SpatialFA%W, SpatialFA%FitCoeff(:,n),            &
                          SpatialFA%KernelAmplitudes(:,n),                 &
                          SpatialFA%bmx, SpatialFA%wmx, SpatialFA%fmx      )
    ENDDO

    ! Compute land cover fractions
    f_land    = SpatialFA%LandFraction*(1.0-SpatialFA%EffectiveSnowFraction)
    f_landice = SpatialFA%LandFraction*SpatialFA%EffectiveSnowFraction
    f_sea     = (1.0-SpatialFA%LandFraction)*(1.0-SpatialFA%SeaIceFraction)
    f_seaice  = (1.0-SpatialFA%LandFraction)*SpatialFA%SeaIceFraction
    
    IF( SpatialFA%DoOceanGlint ) THEN

      ! Isotropic Kernel
      SpatialFA%KernelAmplitudes(:,1)  = SpatialFA%KernelAmplitudes(:,1)*f_land &
                                       + SpatialFA%SnowAlbedo(:)*(f_landice+f_seaice)
      
      ! Volumetric (RossThick)
      SpatialFA%KernelAmplitudes(:,2)  = SpatialFA%KernelAmplitudes(:,2)*f_land

      ! Geometric (LiSparse)
      SpatialFA%KernelAmplitudes(:,3)  = SpatialFA%KernelAmplitudes(:,3)*f_land

      ! Ocean Kernel (New CoxMunk + Glint)
      SpatialFA%KernelAmplitudes(:,4)  = f_sea

    ELSE

      ! Isotropic Kernel
      SpatialFA%KernelAmplitudes(:,1)  = SpatialFA%KernelAmplitudes(:,1)*f_land       &
                                       + SpatialFA%SnowAlbedo(:)*(f_landice+f_seaice) &
                                       + SpatialFA%WaterAlbedo(:)*f_sea
      
      ! Volumetric (RossThick)
      SpatialFA%KernelAmplitudes(:,2)  = SpatialFA%KernelAmplitudes(:,2)*f_land

      ! Geometric (LiSparse)
      SpatialFA%KernelAmplitudes(:,3)  = SpatialFA%KernelAmplitudes(:,3)*f_land

    ENDIF

    ! Compute basis spline expansion coefficients
    DO n=1,SpatialFA%nkern
      CALL SPLINE1(SpatialFA%Wvl,SpatialFA%KernelAmplitudes(:,n),&
                   SpatialFA%wmx,SpatialFA%KernelAmplitudesSP(:,n))
    ENDDO


  END SUBROUTINE compute_RTLS_brdf 

  SUBROUTINE calculate_blueskyalb(nw, fiso,fvol,fgeo, kvol1, kgeo1, kvol2, kgeo2, dirfrac, fracslp, alb) 

    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER (KIND=4), INTENT(IN)                      :: nw
    REAL (KIND=8), INTENT(IN)                         :: kgeo1, kvol1, kgeo2, kvol2
    REAL (KIND=8), DIMENSION(nw), INTENT(IN)          :: fiso,fvol,fgeo
    REAL (KIND=8), DIMENSION(nw), INTENT(IN)          :: dirfrac, fracslp
    REAL (KIND=8), DIMENSION(nw), INTENT(OUT)         :: alb

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                       :: i
    REAL (KIND=8), DIMENSION(nw) :: balb, walb, dalb, frac

    ! =====================================================================
    ! calculate_blueskyalb starts here
    ! =====================================================================

    DO i = 1, nw
      balb(i) = fiso(i) + kvol1 * fvol(i) + kgeo1 * fgeo(i)
      walb(i) = fiso(i) + kvol2 * fvol(i) + kgeo2 * fgeo(i)
    enddo

    alb = balb * dirfrac + walb * ( 1.0 - dirfrac) ! use dirfrac (i.e., albedo=0.001)
    dalb = 1.0
    ! Iterative Derivation of blue-sky albedo as weighting depends on surface albedo
    DO WHILE (ANY(dalb > 1.0E-4) ) 
      frac = dirfrac + fracslp * (alb - 0.001) ! Derive actual fraction of direct irradiance for given alb
      dalb = balb * frac + walb * (1.0 - frac) - alb
      alb = alb + dalb
    ENDDO
    
    RETURN
    
  END SUBROUTINE calculate_blueskyalb
  
  SUBROUTINE ComputeFACoeff(modis_refl, MuMODIS, Mu, G, W, coeff, hyprefl, wmx_modis, wmx, fmx)
    
    ! --------------------
    ! subroutine arguments
    ! --------------------
    INTEGER,      INTENT(IN)                  :: wmx
    INTEGER,      INTENT(IN)                  :: fmx
    INTEGER,      INTENT(IN)                  :: wmx_modis
    REAL(KIND=8), INTENT(IN)                  :: modis_refl(wmx_modis)
    REAL(KIND=8), INTENT(IN)                  :: MuMODIS(wmx_modis)
    REAL(KIND=8), INTENT(IN)                  :: Mu(wmx)
    REAL(KIND=8), INTENT(IN)                  :: G(fmx,wmx_modis)
    REAL(KIND=8), INTENT(IN)                  :: W(wmx,fmx)
    REAL(KIND=8), INTENT(OUT)                 :: coeff(fmx)
    REAL(KIND=8), INTENT(OUT)                 :: hyprefl(wmx)
    
    ! ---------------
    ! Local variables
    ! ---------------
    REAL(KIND=8), DIMENSION(wmx_modis)  :: refl_diff
    
    
    INTEGER :: i,j

    ! ==============================================================
    ! computeFACoeff starts here
    ! ==============================================================
    
    ! Subtract mean
    refl_diff = modis_refl - MuMODIS
      
    ! Compute coefficients
    coeff(:) = 0.0
    DO i=1,fmx
    DO j=1,wmx_modis
      coeff(i) = coeff(i) + G(i,j)*refl_diff(j)
    ENDDO
    ENDDO
      
    ! Initialize with mean
    hyprefl = Mu
    
    ! Add factor loadings
    DO i=1,fmx
      hyprefl = hyprefl + coeff(i)*W(:,i)
    ENDDO
    
  END SUBROUTINE computeFACoeff

END MODULE spatial_fa_module

