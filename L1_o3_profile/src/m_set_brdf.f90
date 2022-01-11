! ===================================================================     
! Search for cloud information and surface albedo
! If has clouds, use GOME derived albedo by Kolemeijer, then use fixed
! albedoes, need to adjust fitvar_rad and fitvar
! If there no clouds, then derive the surface albedo from 370.2 nm
! Need to override the specified albedo values 
! Need to check for no albedo or multiple albedo specified at some wavelengths
module m_set_brdf
   
  USE window_module, ONLY: WinType
  USE level1_def, ONLY: GeolocationType, AuxSurfType
  USE Input_module, ONLY: readinputfile, InpOptType
  USE surface_module, ONLY: SurfOptType, SurfProfType, SurfaceType
  !USE surface_module, ONLY: SurfOptType, SurfProfType, SurfaceType, InitSurface,SampleSurfaceProperties, SetSurfaceLinearization
  USE error_module, ONLY: ErrorType
  USE OMSAO_precision_module, ONLY:sp,dp
  TYPE (InpOptType)  :: InpOpt
  TYPE (SurfOptType) :: SurfOpt
  TYPE (WinType) :: window
  TYPE (GeolocationType) :: geolocation
  TYPE (SurfaceType) :: Surface
  TYPE (SurfProfType) :: SurfProf
  TYPE (AuxSurfType)  :: L2Surface
  TYPE (ErrorType), PUBLIC :: error

  public set_brdf, geolocation
  contains

SUBROUTINE set_brdf (nw, wavs,nactalbspc,landfrac, errstat)
  USE OMSAO_precision_module, ONLY:dp
  USE OMSAO_variables_module, ONLY: numwin, winlim, atmdbdir, ctrdbdir,&
      the_lon, the_lat, the_lons, the_lats, the_jday, the_sza_atm, the_vza_atm,the_aza_atm
  USE ozprof_data_module, ONLY: nalbwf
  USE m_ezspline_interpolation, ONLY:bspline
  ! INPUT/OUTPUT variables
  INTEGER, INTENT(IN) :: nw
  INTEGER, INTENT (INOUT) :: errstat
  REAL (KIND=dp), DIMENSION (nw), INTENT(IN) :: wavs
  REAL (KIND=sp), INTENT(OUT) :: landfrac
  ! Local Variables
  INTEGER :: i
  INTEGER,DIMENSION(4), PARAMETER :: ord = (/1, 3 ,4, 2/)
  CHARACTER(LEN=100) :: ctrfname
  ! Save varaibles
  LOGICAL, SAVE :: first = .true.
    ctrfname = ADJUSTL(TRIM(ctrdbdir))//'sfc.inp'
    nalbwf = 1
    Window%Settings%StartWvl = winlim(1, 1)
    Window%Settings%EndWvl   = winlim(numwin, 2) 
    IF (allocated(Window%RTM_WVL)) deallocate (Window%RTM_Wvl)
    ALLOCATE (Window%RTM_Wvl(nw))
    Window%nRTM_wvl = nw 
    Window%RTM_Wvl  = wavs(1:nw)
    Window%RTMSettings%RTMName =  'VL-LBL'
    Geolocation%Longitude= the_lon
    Geolocation%Latitude = the_lat
    Geolocation%CornerLongitudes = the_lons(ord)
    Geolocation%CornerLatitudes = the_lats(ord)
    Geolocation%Time%DayofYear = the_jday
    Geolocation%sza = the_sza_atm 
    Geolocation%vza = the_vza_atm
    Geolocation%aza = the_aza_atm
  IF (first) THEN 
    CALL readInputFile(ctrfname, SurfOpt, Error) 
    SurfOpt%RootDataDir = ADJUSTL(TRIM(atmdbdir))
    SurfOpt%DoAmplitudeLinearization = .True.
    SurfOpt%DoParameterLinearization = .False.
    SurfProf%WindSpeed = 5.0
    SurfProf%WindDirection = 0.01  
    SurfProf%Chlorophyll = 0.01
    SurfProf%OceanSalinity = 0.01
    SurfProf%SnowDepth = 0.00
    SurfProf%SnowFraction = 0.0
    SurfProf%SnowAge = 0.00
    SurfProf%SeaIceFraction = 0.0
    SurfProf%LandCoverFraction(:) = 0.0
    SurfProf%LandCoverFraction(17) = 1.0 ! Barren Land
    Surface%DoPlantFluorescence = .FALSE.
    CALL Surface%Init(SurfOpt, Window,L2Surface, Error)
    first = .false.
  ENDIF
  !Surface%Option4%W(Wav,EOF), Surface%Option4%WSP(Wav,EOF)
  !Surface%Option4%Mu ! Mean
  !Surface%Option4%FitCoeff(EOF,kernel) ! apriori scaling factors
  !Surface%Option4%FactorUncertainty(EOF,EOF) ! error covariance for FitCoeff
  CALL Surface%SamplePixel (Window, Geolocation, SurfProf, L2Surface, Error)
  nalbwf = 1
  nactalbspc = Surface%option4%fmx
  landfrac=Surface%Option4%LandFraction
  ! interpolation is done at set_cldalb_spc
  !RETURN

  IF (SurfOpt%Option4_DoIsotropic) THEN 
    DO i=1,Surface%Option4%fmx
      ! print*,i,Surface%Option4%FitCoeff(i,1),SQRT(Surface%Option4%FactorUncertainty(i,i))
    ENDDO
    CALL BSPLINE(Surface%Option4%Wvl, Surface%Option4%Mu(:),Surface%Option4%wmx, & 
       wavs(1:nw), Surface%Option4%MuSP(1:nw),nw, errstat)
    DO i=1,Surface%Option4%fmx
       CALL BSPLINE(Surface%Option4%Wvl, Surface%Option4%W(:, i),Surface%Option4%wmx, & 
       wavs(1:nw), Surface%Option4%WSP(1:nw,i),nw, errstat)
    !   print * , i, wavs(nw), Surface%Option4%WSP(nw,i)
    ENDDO
    !print * , 'finish chris surface spectrum loading !!!'
  ELSE
    !PRINT *, Surface%DoLambertian, Surface%fixed_par, surface%fixed_amp
    !DO i = 1 , surface%nkern
    !WRITE(*,'(i2, a16 i2, e15.7, i4, 10e15.7)')  surface%npar(i), ADJUSTL(TRIM(Surface%kern_name(i))), & 
    !          Surface%kern_idx(i), Surface%kern_amp(i,1:surface%npar(i)), & 
    !          Surface%npar(i),Surface%kern_par(i, 1, 1)
    !ENDDO
  ENDIF
    nalbwf = Surface%njac
    WRITE(321,*) nw, surface%nkern
    WRITE(321,'(4A15)') surface%kern_name
    DO i = 1, nw 
      WRITE(321,'(i5, f8.3, 4e15.7)') i,  wavs(i), surface%kern_amp(i, 1:surface%nkern)
    ENDDO 
  !ENDIF
  RETURN
END SUBROUTINE SET_BRDF

end module m_set_brdf
