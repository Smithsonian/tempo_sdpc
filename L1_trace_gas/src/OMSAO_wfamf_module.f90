MODULE OMSAO_wfamf_module

  ! ====================================================================
  ! This module defines variables associated with the wavelength depende
  ! nt AMF calculations and contains necessary subroutines to read files
  ! and calculate them
  ! ====================================================================
  USE OMSAO_precision_module, ONLY: i2, i4, r8, C_LONG, r4
  USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, i2_missval, i4_missval, r4_missval, r8_missval
  use tell_module
  use tio_module
  use ctrlvars, only: yn_gems
  USE OMSAO_he5_module, ONLY: pge_swath_id, &
    he5_start_4d, he5_edge_4d, he5_stride_4d, &
    he5_start_3d, he5_edge_3d, he5_stride_3d, &
    he5_start_2d, he5_edge_2d, he5_stride_2d, &
    he5_start_1d, he5_edge_1d, he5_stride_1d
  USE HDF5, ONLY: HSIZE_T
  IMPLICIT NONE
  private

  public read_climatology_dimensions, amf_calculation, &
    wfamf_deallocate

  ! ---------
  ! PCF stuff
  ! ---------
  INTEGER(KIND=i4), PARAMETER, public :: wfamf_table_lun = 700250
  INTEGER(KIND=i4), PARAMETER, public :: climatology_lun = 700270
  CHARACTER(LEN=MAX_STR_LEN), public  :: OMSAO_wfamf_table_filename
  CHARACTER(LEN=MAX_STR_LEN), public  :: OMSAO_climatology_filename
  integer(kind=i4), parameter, public, dimension(*) :: &
    meteorology_lun = (/700290, 700291/)
  integer(kind=i4), parameter, public :: num_met_luns = size(meteorology_lun)
  CHARACTER(LEN=MAX_STR_LEN), public, dimension(num_met_luns)  :: &
    OMSAO_meteorology_filename

  ! -----------------------------
  ! Dimensions of the climatology
  ! -----------------------------
  INTEGER (KIND=i4), public :: CmETA

  ! =============================
  ! AMF factor specific variables
  ! =============================
  REAL(KIND=r8), public :: amf_wvl, amf_alb_lnd, amf_alb_sno, amf_alb_cld

  ! ------------------------------
  ! amfdiag bit meaning parameters
  ! ------------------------------
  integer(kind=i2), parameter :: yn_amf_geo=0, yn_glint=1, yn_snow=2, &
       yn_cld_cli=3, yn_adj_srf_pre=4, yn_adj_cld_pre=5, yn_albedo=11, yn_cld=12, &
       yn_gas_cli=13, yn_sca=14, yn_amf_cor=15 

  ! ------------------------------
  ! Vlidort lookup table variables
  ! ------------------------------
  ! --------------
  ! Grid variables
  ! --------------
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: lut_alb, lut_sza, lut_vza, lut_srf, lut_wav, &
       lut_alt_lay, lut_alt_lev, lut_pre_lay, lut_pre_lev
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: lut_ozo

  ! -----------------
  ! Profile variables
  ! -----------------
  REAL(KIND=r4), DIMENSION(:,:), ALLOCATABLE :: lut_air_col, lut_ozo_col, lut_tmp

  ! --------------------------
  ! Parameterization variables
  ! --------------------------
  REAL(KIND=r4), DIMENSION(:,:), ALLOCATABLE :: lut_Sb
  REAL(KIND=r4), DIMENSION(:,:,:,:), ALLOCATABLE :: lut_I0, lut_I1, lut_I2, lut_Ir
  REAL(KIND=r4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: lut_dI0, lut_dI1, lut_dI2

  ! -------------------
  ! Dimension variables
  ! -------------------
  INTEGER(HSIZE_T), DIMENSION(1) :: alb_dim, alb_maxdim, ozo_dim, ozo_maxdim, &
       sza_dim, sza_maxdim, vza_dim, vza_maxdim, srf_dim, srf_maxdim, wav_dim, wav_maxdim, &
       lay_dim, lay_maxdim, lev_dim, lev_maxdim

  ! ---------------------------
  ! 32bit/64bit C_LONG integers
  ! ---------------------------
  INTEGER (KIND=C_LONG), PARAMETER :: zerocl = 0, onecl = 1, twocl = 2
  ! Integer parameters
  INTEGER (KIND=i4), PARAMETER :: one = 1, two = 2

CONTAINS

  subroutine wfamf_deallocate (errstat)
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return
    call vlidort_deallocate(errstat)
    if (errstat /= 0) return
  end subroutine wfamf_deallocate

  SUBROUTINE amf_calculation (            &
      pge_idx, nt, nx, lat, lon, sza, vza, saa, vaa, time,  &
      snow, glint, xtrange,    &
      saocol, saodco, saoamf, terrain_height,&
      do_write, errstat                                )

    ! =================================================================
    ! This subroutine computes the AMF factor using the following eleme
    ! nts:
    !     - Kleipool OMLER database
    !     - GEOS Chem climatology
    !     - VLIDORT calculated scattering weights
    ! =================================================================
    USE OMSAO_errstat_module, only: pge_errstat_ok!, pge_errstat_error
    use OMSAO_indices_module, only: voc_omicld_idx
    use OMSAO_omidata_module, only : amf_correction_type
    use output_tools, only : write_albedo, write_gas_profile, &
      write_scattering_weights, write_amf_correction
    USE OMSAO_variables_module,  ONLY: voc_amf_filenames
    use output_tools, only: read_cloud_params
    use ctrlvars, only : yn_stratrop, yn_gems
    use clim_module
    use m_read_gems, only: gems_read_cld
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx, pge_idx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, &
         vza, saa, vaa, terrain_height
    REAL    (KIND=r8), DIMENSION (0:nt-1),      INTENT (IN) :: time
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    LOGICAL                                                 :: do_write
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: saocol, saodco
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT), target :: saoamf
    INTEGER (KIND=i4),                          INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), target :: amfgeo, tropospheric_amf, &
         stratospheric_amf
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), target :: l2cfr, l2ctp
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)       :: albedo
    REAL    (KIND=r8), DIMENSION (CmETA,1:nx,0:nt-1) :: climatology
    REAL    (KIND=r8), DIMENSION (CmETA,1:nx,0:nt-1) :: scattw
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,2) :: cli_wgh_ozo_pro
    INTEGER (KIND=i4), DIMENSION (1:nx,0:nt-1,2) :: cli_idx_ozo_pro
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), target :: surface_pressure, tropopause_pressure

    real    (kind=r4), dimension (:), allocatable, target :: eta_a, eta_b
    integer :: nz

    type (clim_pres_type) :: cpt
    type (clim_cloud_type) :: cct
    type (amf_correction_type) :: amf_corr
    logical :: yn_write_cloud_variables
    character (len=256) :: cloud_file

    ! bitwise like amf calculation flags
    integer (kind=i2), dimension (1:nx,0:nt-1), target :: amfdiag

    if (errstat /= 0) return
    locerrstat  = pge_errstat_ok

    ! ------------------------------------
    ! Initialize variables that are output
    ! ------------------------------------
    albedo       = r8_missval
    climatology  = r8_missval
    cli_wgh_ozo_pro = r8_missval ! Not output
    cli_idx_ozo_pro = i4_missval ! Not output
    scattw       = r8_missval
    saoamf       = r8_missval
    amfgeo       = r8_missval
    amfdiag      = 0
    surface_pressure = r4_missval
    IF (yn_stratrop) then
       tropopause_pressure = r4_missval
       stratospheric_amf = r8_missval
       tropospheric_amf = r8_missval
    ENDIF

    ! -----------------------------------------
    ! If amf_wvl < 0.0 then the slant column is
    ! reported and AMFs equal to 1 with scattw,
    ! akernels and climatology set to missval
    ! -----------------------------------------
    IF (amf_wvl .LT. 0.0) THEN
       saoamf = 1.0_r8
       amfdiag=ibset(amfdiag,yn_amf_cor)
    ELSE

       ! -----------------
       ! Geolocation check
       ! -----------------
       call tell_log (1, 'amf_calculation: check geolocation information')
       call check_geolocation ( nt, nx, lat, lon, sza, vza, terrain_height, amfdiag )
      
       ! -------------------------
       ! Compute the geometric AMF
       ! -------------------------
       call tell_log (1, 'amf_calculation: compute geometric amf')
       CALL compute_geometric_amf ( nt, nx, sza, vza, amfgeo, amfdiag )

       ! -------------------------------------------------------
       ! Initialize molecular AMF with geometric AMF. Subsequent
       ! subroutines will replace any entries where the true
       ! molecular AMF can be computed.
       ! -------------------------------------------------------
       saoamf = amfgeo

       ! ------------------------------------------------------------------
       ! Read VLIDORT look up table. Variables are declared at module level
       ! (Input is read only on the first pass. Subsequent passes use
       ! cached values)
       ! ------------------------------------------------------------------
       call tell_log (1, 'amf_calculation: read scattering weights LUT')
       CALL read_vlidort (errstat)
       if (errstat /= 0) then
          call tell_error (tell_io_read_error, 'reading scattering weights LUT', errstat) 
          call vlidort_deallocate(errstat)
          return
       endif

       ! --------------------------------------------------------------
       ! Read and interpolate albedo database. If no albedo information
       ! set amfdiag bit 2. Set amfdiag bit 3 for glint.
       ! ------------------------------------------------------------
       call tell_log (1, 'amf_calculation: read and prepare albedo')
       call read_albedo ( nt, nx, lat, lon, glint, amfdiag, &
            albedo, errstat)
       if (errstat /= 0) then
          call tell_error (tell_io_read_error, "reading albedo", errstat)
          return
       endif

       ! -------------------------------
       ! Apply snow correction to albedo
       ! -------------------------------
       call tell_log (1, 'amf_calculation: snow correction')
       call snow_correction ( nt, nx, snow, albedo, amfdiag )

       ! ---------------------------------------
       ! Write the albedo to the output file he5
       ! ---------------------------------------
       IF (do_write) then
          call tell_log (1, 'amf_calculation: write albedo to L2 file')
          call write_albedo (albedo, nx, nt, errstat)
          if (errstat /= 0) return
       endif

       ! ---------------------
       ! Read L2 cloud product
       ! ---------------------
       cloud_file = voc_amf_filenames(voc_omicld_idx)
       call tell_log (1, 'amf_calculation: read cloud-top pressure, cloud fraction from '//trim(cloud_file))
       if (.not. yn_gems) then
         call read_cloud_params (cloud_file, nt, nx, l2cfr, l2ctp, errstat)
       else
         call gems_read_cld (cloud_file, nt, nx, l2cfr, l2ctp, errstat)
       endif
       if (errstat /= 0) then
          call tell_error (tell_io_read_error, "reading cloud file: "//trim(cloud_file), errstat)
          return
       endif


       ! Read cloud climatology
       call tell_log (1, 'amf_calculation: initialize cloud climatology')
       call clim_cloud_init (cct, errstat)
       if (errstat /= 0) then
          call tell_error (tell_io_read_error, "reading cloud pressure climatology", errstat)
          return
       end if
       ! ---------------------------------------------------
       ! Use cloud climatology if needed and set cloud flags
       ! ---------------------------------------------------              
       call tell_log (1, 'amf_calculation: read cloud climatology')
       call read_cloud_climatology (nt, nx, lat, lon, time, cct, &
            l2cfr, l2ctp, amfdiag, errstat)

       ! ------------------------------------------------
       ! Read climatology and interpolate to lon/lat/time
       ! ------------------------------------------------
       call tell_log (1, 'amf_calculation: read gas profile climatology')
       CALL get_climatology (cpt, pge_idx, climatology, cli_wgh_ozo_pro, &
            cli_idx_ozo_pro, lat, lon, time, nt, nx, errstat, amfdiag)
       if (errstat /= 0) then
          call tell_error (tell_io_read_error, 'reading gas profile climatology', errstat)
          return
       endif

       ! -------------------------------------
       ! Write the climatology to the he5 file
       ! -------------------------------------
       IF (do_write) then
          call tell_log (1, 'amf_calculation: write gas profile climatology to L2 file')
          call write_gas_profile (climatology, nx, nt, CmETA, errstat)
          if (errstat /= 0) return
       endif

       ! --------------------------------------------------------
       ! Compute Scattering weights in the look up table grid but
       ! with the correct albedo. amfdiag is used to skip pixel
       ! --------------------------------------------------------
       call tell_log (1, 'amf_calculation: compute scattering weights')
       CALL compute_scatt (cpt, nt, nx, time, albedo, sza, vza, saa, vaa, l2ctp, l2cfr, &
            terrain_height, surface_pressure, cli_wgh_ozo_pro, cli_idx_ozo_pro, &
            lat, lon, xtrange, amfdiag, scattw)

       ! -----------------------------------------------------------------
       ! Work out the AMF using the scattering weights and the climatology
       ! Work out Averaging Kernels
       ! -----------------------------------------------------------------
       call tell_log (1, 'amf_calculation: compute amfs')
       CALL compute_amf (cpt,  nt, nx, CmETA, climatology, &
                         scattw, saoamf, stratospheric_amf, tropospheric_amf, &
                         surface_pressure, tropopause_pressure, lat, lon, amfdiag, &
                         locerrstat)

       ! -----------------------------------------------------------------
       ! Write out scattering weights, altitude grid and averaging kernels
       ! -----------------------------------------------------------------
       IF (do_write) then
          call tell_log (1, 'amf_calculation: write scattering weights to L2 file')
          call write_scattering_weights (scattw, nx, nt, CmETA, errstat)
          if (errstat /= 0) then
             call tell_error (tell_io_read_error, 'writting scattering weights to L2 file', errstat)
             return
          endif
       endif

    END IF

    ! --------------------------
    ! Apply the air mass factors
    ! --------------------------
    WHERE ( saoamf > 0.0_r8 .AND. saocol > r8_missval .AND. saodco > r8_missval )
      saocol = saocol / saoamf
      saodco = saodco / saoamf
    END WHERE

    ! -----------------------------------------------
    ! Write AMFs, AMF diagnosting, and AMF-adjusted
    ! columns and column uncertainties to output file
    ! -----------------------------------------------
    IF (do_write) then
      call tell_log (1, 'amf_calculation: write amf correction to L2 file')
      nz = clim_pres_nz (cpt)
      allocate (eta_a(nz), eta_b(nz))
      call clim_pres_eta (cpt, eta_a, eta_b, errstat)
      if (errstat /= 0) return
      amf_corr % amf_molecule_specific => saoamf
      amf_corr % amf_molecule_stratospheric => stratospheric_amf
      amf_corr % amf_molecule_tropospheric => tropospheric_amf
      amf_corr % amf_geometric => amfgeo
      amf_corr % diagnostic_flag => amfdiag
      amf_corr % cloud_fraction => l2cfr
      amf_corr % cloud_pressure => l2ctp
      amf_corr % surface_pressure => surface_pressure
      amf_corr % tropopause_pressure => tropopause_pressure
      amf_corr % eta_a => eta_a
      amf_corr % eta_b => eta_b
      yn_write_cloud_variables = .TRUE.
      call write_amf_correction (nx, nt, amf_corr, saocol, saodco, &
                                 yn_write_cloud_variables, errstat)
      if (errstat /= 0) then
         call tell_error (tell_io_read_error, 'writting amf correction to L2 file', errstat)
         return
      endif
   endif

  END SUBROUTINE amf_calculation

  subroutine check_geolocation ( nt, nx, lat, lon, sza, vza, terrain_height, amfdiag )
    implicit none
    integer (kind=i4), intent (in) :: nt, nx
    real (kind=r4), dimension (1:nx,0:nt-1), intent (in) :: lat, lon, terrain_height
    real (kind=r4), dimension (1:nx,0:nt-1), intent (in) :: sza, vza
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (out) :: amfdiag
    real (kind=r4), parameter :: terrain_height_missval=-32737.0

    ! Check that a complete set of geolocation information is available to
    ! complete the AMF calculation. SZA and VZA have to be between 0 and 90
    ! and latitude and longitude have to be not equal to r4_missval
    ! Pixels without complete geolocation information get amfdiag bit 0 set
    ! FIXME, some pixels have no terrain height information. Those are given
    ! the value -32767. Will be good to know 
    where ( &
         sza(1:nx,0:nt-1) == r4_missval .or. &
         sza(1:nx,0:nt-1) < 0.0_r4 .or. &
         sza(1:nx,0:nt-1) > 90.0_r4 .or. &
         vza(1:nx,0:nt-1) == r4_missval .or. &
         vza(1:nx,0:nt-1) < 0.0_r4 .or. &
         vza(1:nx,0:nt-1) > 90.0_r4 .or. &
         lat(1:nx,0:nt-1) == r4_missval .or. &
         lon(1:nx,0:nt-1) == r4_missval .or. &
         terrain_height(1:nx,0:nt-1) == terrain_height_missval )
       amfdiag(1:nx,0:nt-1) = ibset(amfdiag(1:nx,0:nt-1),yn_amf_cor)
    end where
  end subroutine check_geolocation

  subroutine compute_geometric_amf ( nt, nx, sza, vza, amfgeo, amfdiag )
    use OMSAO_parameters_module, only: deg2rad

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    integer (kind=i4),                         intent (IN) :: nx, nt
    real    (kind=r4), dimension (nx,0:nt-1),  intent (IN) :: sza, vza

    ! ----------------
    ! Output variables
    ! ----------------
    real (kind=r8), dimension (1:nx,0:nt-1), intent (OUT) :: amfgeo
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (OUT) :: amfdiag

    ! ---------------------------------------------------
    ! Compute geometric AMF and set diagnostic flag bit 1
    ! ---------------------------------------------------
    where (.not. btest(amfdiag(1:nx,0:nt-1),yn_amf_cor))
       amfgeo(1:nx,0:nt-1) = &
            1.0_r8 / cos ( real(sza(1:nx,0:nt-1),KIND=r8)*deg2rad ) + &
            1.0_r8 / cos ( real(vza(1:nx,0:nt-1),KIND=r8)*deg2rad )
       amfdiag(1:nx,0:nt-1) = ibset(amfdiag(1:nx,0:nt-1),yn_amf_geo)
    end where
    return
  end subroutine compute_geometric_amf

  subroutine read_albedo ( nt, nx, lat, lon, glint, amfdiag, &
            albedo, errstat)

    ! ==================================================================
    ! This subroutine reads the OMLER albedo data base for the month of
    ! the orbit to processed. Then it interpolates the values for each
    ! one of the pixels of the orbit to be analyzed
    ! ==================================================================
    use OMSAO_linterpolation_module, only: lininterpol, GetNode
    use OMSAO_variables_module, only: OMSAO_OMLER_filename, &
      winwav_min, winwav_max
    use ezspline_interpolation, only: ezspline_1d_interpolation, &
      ezspline_2d_interpolation
    use OMSAO_errstat_module, only : he5_stat_fail, pge_errstat_ok
    use OMSAO_he5_module, only: HE5_GDOPEN, HE5_GDattach, HE5_GDRDFLD, &
         HE5_GDRDLATTR, HE5_GDDETACH, HE5_GDclose, he5f_acc_rdonly, &
         granule_month

    implicit none

    ! ---------------
    ! Input variables
    ! ---------------
    integer (kind=i4), intent (in) :: nt, nx
    real (kind=r4), dimension (1:nx,0:nt-1), intent (in) :: lat, lon
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (in) :: glint

    ! ------------------
    ! Modified variables
    ! ------------------
    integer (kind=i4), intent (inout) :: errstat
    real (kind=r8), dimension (1:nx,0:nt-1), intent (inout) :: albedo

    ! ----------------
    ! Output variables
    ! ----------------
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (OUT) :: amfdiag

    ! ------------------------------------------------------------------
    ! Local variables, the variables to hold the OMLER data are going to
    ! be allocated and deallocated within this subroutine.
    ! ------------------------------------------------------------------
    REAL    (KIND=r4), ALLOCATABLE, DIMENSION(:) :: OMLER_longitude,   &
      OMLER_latitude,    &
      OMLER_wvl
    INTEGER (KIND=i2), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
      OMLER_monthly_albedo
    REAL    (KIND=r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
      OMLER_wvl_albedo
    REAL    (KIND=r8), ALLOCATABLE, DIMENSION(:,:) :: &
      OMLER_albedo

    ! --------------------
    ! More Local variables
    ! --------------------
    CHARACTER (LEN=34), PARAMETER :: grid_name = 'EarthSurfaceReflectanceClimatology'
    CHARACTER (LEN=MAX_STR_LEN) :: grid_file
    INTEGER (KIND=i4), PARAMETER :: OMLER_n_latitudes = 360, &
      OMLER_n_longitudes = 720, OMLER_n_wavelenghts =  23, one = 1
    INTEGER (KIND=i4) :: itimes, ixtrack, ilon, ilat, nlon, &
      nlat, OMnwvl, grid_id, grid_file_id, month, minwvl, maxwvl
    INTEGER (KIND=i4), DIMENSION(2) :: lon_idx, lat_idx
    REAL (KIND=r4) :: scale_factor, offset
    REAL (KIND=r8) :: lonp, latp

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: locerrstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=16), PARAMETER :: modulename = 'omi_omler_albedo'

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    locerrstat = pge_errstat_ok
    if (errstat /= 0) return

    grid_file = TRIM(ADJUSTL(OMSAO_OMLER_filename))

    ! -------------------------------------------------------------------------------
    ! Open he5 OMI OMLER grid file and check GRID_FILE_ID (-1 if error)
    ! -------------------------------------------------------------------------------
    grid_file_id = HE5_GDOPEN (grid_file, he5f_acc_rdonly)
    IF (grid_file_id == he5_stat_fail) THEN
      call tell_error (tell_io_open_error, "omi_omler_albedo: opening Omler grid file"//trim(grid_file), &
                       errstat)
      RETURN
    END IF

    ! ----------------------------------------------
    ! Attach to grid and check GRID_ID (-1 if error)
    ! ----------------------------------------------
    grid_id = HE5_GDattach (grid_file_id, grid_name)
    IF (grid_id == he5_stat_fail) THEN
      call tell_error (tell_io_read_error, "omi_omler_albedo: attaching to Omler grid file"// &
                       trim(grid_file), errstat)
      return
    END IF

    ALLOCATE (OMLER_longitude(OMLER_n_longitudes), &
              OMLER_latitude(OMLER_n_latitudes), &
              OMLER_wvl(OMLER_n_wavelenghts), &
              OMLER_albedo(OMLER_n_longitudes, OMLER_n_latitudes), &
              stat=locerrstat)
    if (locerrstat /= 0) then
      call tell_error (tell_malloc_error, "omi_omler_albedo:  allocate failed", &
                       errstat)
      return
    endif

    ! -------------------------
    ! Read longitude data field
    ! -------------------------
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_longitudes
    locerrstat = HE5_GDRDFLD(grid_id, "Longitude", he5_start_1d, he5_stride_1d, &
      he5_edge_1d, OMLER_longitude)
    errstat = MAX (errstat, locerrstat)

    ! ------------------------
    ! Read latitude data field
    ! ------------------------
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_latitudes
    locerrstat = HE5_GDRDFLD(grid_id, "Latitude", he5_start_1d, he5_stride_1d, &
      he5_edge_1d, OMLER_latitude)
    errstat = MAX (errstat, locerrstat)

    ! ---------------------------
    ! Read wavelenghts data field
    ! ---------------------------
    he5_start_1d = 0; he5_stride_1d = 1; he5_edge_1d = OMLER_n_wavelenghts
    locerrstat = HE5_GDRDFLD(grid_id, "Wavelength", he5_start_1d, he5_stride_1d, &
      he5_edge_1d, OMLER_wvl)
    errstat = MAX (errstat, locerrstat)

    ! -----------------------------------------------
    ! Select the wavelenghts; finding array positions
    ! -----------------------------------------------
    minwvl = MINLOC(OMLER_wvl, 1, OMLER_wvl .GE. REAL(winwav_min,KIND=r4))
    maxwvl = MAXLOC(OMLER_wvl, 1, OMLER_wvl .LE. REAL(winwav_max,KIND=r4))
    OMnwvl = maxwvl-minwvl+1

    allocate (OMLER_monthly_albedo(OMLER_n_longitudes, OMLER_n_latitudes, OMnwvl,1), &
              OMLER_wvl_albedo(OMLER_n_longitudes, OMLER_n_latitudes, OMnwvl,1), &
              stat=locerrstat)
    if (locerrstat /= 0) then
      call tell_error (tell_malloc_error, "omi_omler_albedo:  allocate failed", &
                       errstat)
      return
    endif

    ! ---------------------------
    ! Read the albedo data field:
    ! -Month
    ! -Selected wavelenghts
    ! ---------------------------
    month = granule_month - 1

    he5_start_4d  = (/  0,  0, minwvl-1, month/)
    he5_stride_4d = (/  1,  1,        1,     1/)
    he5_edge_4d   = (/OMLER_n_longitudes,OMLER_n_latitudes, &
      OMnwvl, 1/)
    locerrstat = HE5_GDRDFLD(grid_id, "MonthlyMinimumSurfaceReflectance",  &
      he5_start_4d, he5_stride_4d, he5_edge_4d, OMLER_monthly_albedo)
    errstat = MAX (errstat, locerrstat)

    ! ----------------------------
    ! Read the albedo scale factor
    ! ----------------------------
    locerrstat = HE5_GDRDLATTR(grid_id, "MonthlyMinimumSurfaceReflectance", &
      "ScaleFactor", scale_factor)
    locerrstat = HE5_GDRDLATTR(grid_id, "MonthlyMinimumSurfaceReflectance", &
      "Offset", offset)

    ! --------------------
    ! Deattached from grid
    ! --------------------
    locerrstat = HE5_GDDETACH(grid_id)

    ! -------------------------------------------
    ! Close he5 OMI OMLER grid file (-1 if error)
    ! -------------------------------------------
    locerrstat = HE5_GDclose ( grid_file_id)
    IF ( locerrstat == he5_stat_fail) THEN
      call tell_error (tell_io_error, "omi_omler_albedo: closing Omler grid file"// &
                       trim(grid_file), errstat)
      return
    END IF

    OMLER_wvl_albedo = real(offset, KIND = r8) +          &
      real(scale_factor, KIND = r8)*      &
      real(OMLER_monthly_albedo, KIND=r8)

    ! ------------------------------------------------
    ! Interpolate for each pixel to amf_wvl wavelenght
    ! ------------------------------------------------
    DO ilon = 1, OMLER_n_longitudes
      DO ilat = 1, OMLER_n_latitudes
        CALL ezspline_1d_interpolation ( &
          OMnwvl, REAL(OMLER_wvl(minwvl:maxwvl), KIND=r8), &
          OMLER_wvl_albedo(ilon,ilat,1:OMnwvl,1), &
          one, [amf_wvl], OMLER_albedo(ilon,ilat), locerrstat )
      END DO
    END DO

    ! --------------------------------------------------
    ! Interpolate to the lat and longitude of each pixel
    ! --------------------------------------------------
    do itimes = 0, nt-1
       do ixtrack = 1, nx

          ! Convert longitude/latitude to real 8 for interpolation
          lonp = real(lon(ixtrack,itimes), kind=r8)
          latp = real(lat(ixtrack,itimes), kind=r8)

          ! Only work out surface reflectance if the albedo database contains information for
          ! that pixel. The interpolation method will have to change if we move to a higher
          ! resolution database
          if ( lonp < minval(OMLER_longitude) .or. lonp > maxval(OMLER_longitude) .or. &
               latp < minval(OMLER_latitude) .or. latp > maxval(OMLER_latitude) ) then
             amfdiag(ixtrack,itimes) = ibset(amfdiag(ixtrack,itimes),yn_albedo)
             cycle
          end if

          ! -----------------------------------------
          ! Locate two closest indices to lon and lat
          ! in OMLER_longitude and OMLER_latitudes.
          ! If result out of bounds bring it to the
          ! closest boundary.
          ! ------------------------------------------
          CALL GetNode(REAL(OMLER_longitude,KIND=r8),lonp, &
               lon_idx(1), 'Lower')
          IF (lon_idx(1) .EQ. -2) lon_idx(1) = 1
          IF (lon_idx(1) .EQ. -3) lon_idx(1) = OMLER_n_longitudes-1
          CALL GetNode(REAL(OMLER_longitude,KIND=r8),lonp, &
               lon_idx(2), 'Upper')
          IF (lon_idx(2) .EQ. -2) lon_idx(2) = 2
          IF (lon_idx(2) .EQ. -3) lon_idx(2) = OMLER_n_longitudes
          nlon = lon_idx(2)-lon_idx(1)+1
          IF (nlon == 1) THEN
             IF (lon_idx(1) <= OMLER_n_longitudes-1) THEN
                lon_idx(2) = lon_idx(1) + 1
             ELSE IF (lon_idx(1) == OMLER_n_longitudes) THEN
                lon_idx(1) = OMLER_n_longitudes-1
                lon_idx(2) = OMLER_n_longitudes
             ENDIF
             nlon = lon_idx(2)-lon_idx(1)+1
          ENDIF

          CALL GetNode(REAL(OMLER_latitude,KIND=r8),latp, &
               lat_idx(1), 'Lower')
          IF (lat_idx(1) .EQ. -2) lat_idx(1) = 1
          IF (lat_idx(1) .EQ. -3) lat_idx(1) = OMLER_n_latitudes-1
          CALL GetNode(REAL(OMLER_latitude,KIND=r8),latp, &
               lat_idx(2), 'Upper')
          IF (lat_idx(2) .EQ. -2) lat_idx(2) = 2
          IF (lat_idx(2) .EQ. -3) lat_idx(2) = OMLER_n_latitudes
          nlat = lat_idx(2)-lat_idx(1)+1
          IF (nlat == 1) THEN
             IF (lat_idx(1) <= OMLER_n_latitudes-1) THEN
                lat_idx(2) = lat_idx(1) + 1
             ELSE IF (lat_idx(1) == OMLER_n_latitudes) THEN
                lat_idx(1) = OMLER_n_latitudes-1
                lat_idx(2) = OMLER_n_latitudes
             ENDIF
             nlat = lat_idx(2)-lat_idx(1)+1
          ENDIF

          albedo(ixtrack,itimes) = linInterpol(nlon,nlat,&
               real(OMLER_longitude(lon_idx(1):lon_idx(2)),KIND=r8), &
               real(OMLER_latitude(lat_idx(1):lat_idx(2)),KIND=r8), &
               OMLER_albedo(lon_idx(1):lon_idx(2),lat_idx(1):lat_idx(2)), &
               lonp, latp, status=locerrstat)
          if (locerrstat /= 0) then
             call tell_error (tell_runtime_error, &
                  "omi_omler_albedo: lon/lat interpolation failed", errstat)
             return
          endif
       end do
    end do

    ! --------------------
    ! Deallocate variables
    ! --------------------
    DEALLOCATE (OMLER_monthly_albedo, OMLER_wvl_albedo, OMLER_albedo, &
         OMLER_longitude, OMLER_latitude, OMLER_wvl, stat=locerrstat)
    if (locerrstat /= 0) then
       call tell_error (tell_malloc_error, &
            "omi_omler_albedo: deallocate failed", errstat)
       return
    endif
    
    ! ---------------------------
    ! Set amfdiag bit 3 for glint
    ! ---------------------------
    where (glint(1:nx,0:nt-1) /= 0)
       amfdiag(1:nx,0:nt-1) = ibset(amfdiag(1:nx,0:nt-1),yn_glint)
    end where

    errstat = max(errstat, locerrstat)

  end subroutine read_albedo

  subroutine snow_correction ( nt, nx, snow, albedo, amfdiag )

    use OMSAO_omidata_module, only: &
         NISE_snowfree, NISE_permice

    implicit none
    ! ---------------
    ! Input variables
    ! ---------------
    integer (kind=i4), intent (in) :: nt, nx
    integer (KIND=i2), dimension (1:nx,0:nt-1), intent (in) :: snow

    ! ------------------
    ! Modified variables
    ! ------------------
    real (kind=r8), dimension (1:nx,0:nt-1), intent (inout) :: albedo

    ! ----------------
    ! Output variables
    ! ----------------
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (out) :: amfdiag

    ! Local variables
    integer (kind=i4) :: ix, it
    real (kind=r8) :: frc

    ! This are the values we are using for the correction:
    ! snowfree =   0, all snow = 100, permanent ice = 101, dry snow = 103
    ! ocean    = 104, suspect  = 125, error         = 127
    ! However, it may be worth to update to more recent NISE products or IMS.
    ! High temporal resolution GLER or BRDF products may provide better reflectance
    ! estimates that the current amf_alb_sno constant value currently applied.
    ! The correction uses independent pixel approximation for partially covered
    ! snow pixels and leaves un-altered the albedo for pixels identified as
    ! permanent_ice
    do ix = 1, nx
       do it = 0, nt-1
          if (snow(ix,it) == NISE_snowfree .or. snow (ix,it) >= NISE_permice) then
             cycle
          else if (snow(ix,it) > NISE_snowfree .and. snow(ix,it) < NISE_permice) then
             frc = real(snow(ix,it),kind=r8)/100.0_r8
             albedo(ix,it) = (1.0_r8-frc) * albedo (ix,it) + frc * amf_alb_sno
             amfdiag(ix,it) = ibset(amfdiag(ix,it),yn_snow)
          end if
       end do
    end do

  end subroutine snow_correction

  subroutine read_cloud_climatology (nt, nx, lat, lon, time, cct, &
       l2cfr, l2ctp, amfdiag, errstat)
    use clim_module
    implicit none

    ! ---------------
    ! Input variables
    ! ---------------
    integer (kind=i4), intent (in) :: nt, nx
    real (kind=r8), dimension (1:nx,0:nt-1), intent (in) :: l2cfr
    real (KIND=r4), dimension (1:nx,0:nt-1), intent (in) :: lat, lon
    real (kind=r8), dimension (0:nt-1), intent (in) :: time
    type (clim_cloud_type), intent(in) :: cct
    ! ------------------
    ! Modified variables
    ! ------------------
    integer (kind=i4), intent (inout) :: errstat
    ! -----------------
    ! Ouptput variables
    ! -----------------
    real (kind=r8), dimension (1:nx,0:nt-1), intent(out) :: l2ctp
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (out) :: amfdiag

    ! Local variables
    integer :: year, month, day
    integer (kind=i4) :: ix, it
    real (kind=r8) :: hour, tai93_offset
    real (kind=r4) :: pressure
    ! If gems, need to offset for different epoch
    tai93_offset=0.0d0
    if (yn_gems) then
      call tiof_time_set_taix_epoch ("1993-01-01T00:00:00Z", errstat)
      if (errstat == 0) then
        call tiof_utcstr_to_taix_time ("2000-01-01T00:00:00Z", tai93_offset, &
             errstat)
      endif
    endif
    ! ---------------------------------------------------------------
    ! Check if we have cloud information from satellite retrievals.
    ! If pressure is incomplete (or below 10 hPa assuming no clouds
    ! exist above that pressure then correct/complete information
    ! with cloud climatology.
    ! ---------------------------------------------------------------
    where ( l2cfr(1:nx,0:nt-1) < 0.0_r8 .or. l2cfr(1:nx,0:nt-1) > 1.0_r8 )
       amfdiag(1:nx,0:nt-1) = ibset(amfdiag(1:nx,0:nt-1),yn_cld)
    elsewhere ( l2ctp(1:nx,0:nt-1) < real(minval(lut_srf),kind=r8) .or. &
         l2ctp(1:nx,0:nt-1) > real(maxval(lut_srf),kind=r8) )
       amfdiag(1:nx,0:nt-1) = ibset(amfdiag(1:nx,0:nt-1),yn_cld_cli)
    endwhere
    !reset epoch if gems data is in use
    do ix=1,nx
       do it=0,nt-1
          if ( btest(amfdiag(ix,it),yn_cld_cli) ) then
             call tio_f_taix_time_to_utc_caldate (time(it)-tai93_offset, year, month, day, hour)
             call clim_cloud (cct, month, day, lon(ix,it), lat(ix,it), pressure, errstat)
             if (errstat /= 0) then
                call tell_error (tell_io_read_error, 'reading cloud pressure climatology', errstat)
                return
             end if
             l2ctp(ix,it) = real(pressure,kind=r8)
          else
             cycle
          endif
       end do
    end do
  end subroutine read_cloud_climatology

  subroutine clim_get_climatology (cpt, pge_idx, climatology, cli_wgh_ozo_pro, &
                                   cli_idx_ozo_pro, lat, lon, time, nt, nx, &
                                   errstat, amfdiag)
    use clim_module
    use omsao_indices_module, only: sao_molecule_names
    use ctrlvars, only: yn_gems
    implicit none

    type (clim_pres_type), intent(inout) :: cpt
    integer (kind=i4), intent(in) :: pge_idx
    real (kind=r8), dimension(cmeta,1:nx,0:nt-1), intent (inout) :: climatology
    real (kind=r8), dimension(1:nx,0:nt-1, 2), intent (inout) :: cli_wgh_ozo_pro
    integer (kind=i4), dimension(1:nx,0:nt-1, 2), intent (inout) :: cli_idx_ozo_pro
    real (kind=r4), dimension (1:nx,0:nt-1), intent (in) :: lat, lon
    real (kind=r8), dimension (0:nt-1), intent (in) :: time
    integer (kind=i4), intent (in) :: nt, nx
    integer (kind=i4), intent (inout) :: errstat
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (out) :: amfdiag

    type (clim_pres_bounds_type) :: bounds
    type (clim_species_type) :: cst
    integer :: year(2), month(2), day(2)
    integer :: nz, nlayers, itimes, ixtrack
    real (kind=r8) :: t_beg, t_end, hour, hour_beg, hour_end, tai93_offset
    real (kind=r4), dimension(:), allocatable :: pres, vmr, partial_column
    real (kind=r4) :: hour_f, lon_f, lat_f
    character (len=6) :: clim_db_molecule_name
    real (kind=r4), dimension (1:nx,0:nt-1) :: fudge_lon, fudge_lat


    if (errstat /= 0) return

    ! If gems, need to offset for different epoch
    tai93_offset=0.0d0
    if (yn_gems) then
      call tiof_time_set_taix_epoch ("1993-01-01T00:00:00Z", errstat)
      if (errstat == 0) then
        call tiof_utcstr_to_taix_time ("2000-01-01T00:00:00Z", tai93_offset, &
             errstat)
      endif
    endif

    t_beg = minval(time-tai93_offset, time /= r8_missval)
    t_end = maxval(time-tai93_offset, time /= r8_missval)

    if (t_end - t_beg > 86400.0) then
      call tell_error (tell_runtime_error, "libclim_climatology: granule duration exceeds 24 hours", errstat)
      return
    endif

    call tio_f_taix_time_to_utc_caldate (t_beg, year(1), month(1), day(1), hour_beg)
    call tio_f_taix_time_to_utc_caldate (t_end, year(2), month(2), day(2), hour_end)

    if (yn_gems) then
      call tell_log (1,"WARNING: CALIBRATION NOT IN PLACE FOR GEMS, FUDGING LAT & LON LIMITS")
      fudge_lon = 0.0-lon
      fudge_lat = lat+10.0
    else
      fudge_lon=lon
      fudge_lat=lat
    endif

    bounds % hour_beg = real (hour_beg, kind=r4)
    bounds % hour_end = real (hour_end, kind=r4)
    bounds % lon_min = minval(fudge_lon, fudge_lon /= r4_missval)
    bounds % lon_max = maxval(fudge_lon, fudge_lon /= r4_missval)
    bounds % lat_min = minval(fudge_lat, fudge_lat /= r4_missval)
    bounds % lat_max = maxval(fudge_lat, fudge_lat /= r4_missval)

    call clim_pres_init (cpt, month(1), day(1), bounds, errstat)
    if (errstat /= 0) return
    nz = clim_pres_nz (cpt)
    nlayers = nz - 1
    allocate (pres(nz), vmr(nlayers), partial_column(nlayers))

    clim_db_molecule_name = sao_molecule_names(pge_idx)

    ! We can't agree on how to spell the names of molecules.
    if (clim_db_molecule_name == 'HCHO') then
      clim_db_molecule_name = 'CH2O  '
    endif

    call clim_species_init (cst, cpt, trim(clim_db_molecule_name), errstat)
    if (errstat /= 0) then
       call tell_error ( tell_io_read_error, "libclim_climatology: initializing "//trim(clim_db_molecule_name), errstat)
       return
    end if
    ! FIXME. Instead of using fix values for cli_wgh_ozo_pro and
    ! cli_idx_ozo_pro we can use the lat/lon and ozone total column
    ! to decide at runtime. However, that requires to have O3 clima
    ! tologies not available now
!!$    call clim_species_init (cst_o3, cpt, 'O3', errstat)
!!$    if (errstat /= 0) then
!!$       call tell_error ( tell_io_read_error, "libclim_climatology: initializing O3", errstat)
!!$    end if
    
    do itimes = 0, nt-1
       do ixtrack = 1, nx
          ! Skip this pixel if geolocation information is not available
          if (btest(amfdiag(ixtrack,itimes),yn_amf_cor)) then
             amfdiag(ixtrack,itimes) = ibset(amfdiag(ixtrack,itimes),yn_gas_cli)
             cycle
          end if

          ! Work out hour of interest
          call tio_f_taix_time_to_utc_caldate (time(itimes)-tai93_offset, year(1), month(1), day(1), hour)
          hour_f = real (hour, kind=r4)
          
          lon_f = fudge_lon(ixtrack,itimes)
          lat_f = fudge_lat(ixtrack,itimes)

          ! FIXME - this is just temporary
          cli_wgh_ozo_pro(ixtrack,itimes,1:2) = 0.5
          cli_idx_ozo_pro(ixtrack,itimes,1:2) = 8
          ! Get pressure grid
          call clim_pres (cpt, hour_f, lon_f, lat_f, pres, errstat)
          if (errstat /= 0) then
             call tell_error (tell_runtime_error, "libclim_climatology: calculating pressure", errstat)
             amfdiag(ixtrack,itimes) = ibset(amfdiag(ixtrack,itimes),yn_gas_cli)
             cycle
          end if
          ! Get vmr profile
          call clim_species_vmr (cst, cpt, hour_f, lon_f, lat_f, vmr, errstat)
          if (errstat /= 0) then
             call tell_error (tell_runtime_error, "libclim_climatology: calculating vmr", errstat)
             amfdiag(ixtrack,itimes) = ibset(amfdiag(ixtrack,itimes),yn_gas_cli)
             cycle
          end if
          ! Compute partical columns
          call clim_partial_column (pres, vmr, partial_column, errstat)
          if (errstat /= 0) then
             call tell_error (tell_runtime_error, "libclim_climatology: calculating partial column", errstat)
             amfdiag(ixtrack,itimes) = ibset(amfdiag(ixtrack,itimes),yn_gas_cli)
             cycle
          end if
          ! Fix non-physical partial columns
          where (partial_column < 0.0_r8)
             partial_column = 0.0_r8
          end where
          ! Assign climatology values
          climatology(1:nlayers,ixtrack,itimes) = real (partial_column(1:nlayers), kind=r8)
       enddo
    enddo
  end subroutine clim_get_climatology

  subroutine get_climatology (cpt, pge_idx, climatology, cli_wgh_ozo_pro, &
                                   cli_idx_ozo_pro, lat, lon, time, nt, nx, &
                                   errstat, amfdiag)
    use clim_module
    implicit none
    type (clim_pres_type), intent(inout) :: cpt
    integer (kind=i4), intent(in) :: pge_idx
    real (kind=r8), dimension(cmeta,1:nx,0:nt-1), intent (inout) :: climatology
    real (kind=r8), dimension(1:nx,0:nt-1, 2), intent (inout) :: cli_wgh_ozo_pro
    integer (kind=i4), dimension(1:nx,0:nt-1, 2), intent (inout) :: cli_idx_ozo_pro
    real (kind=r4), dimension (1:nx,0:nt-1), intent (in) :: lat, lon
    real (kind=r8), dimension (0:nt-1), intent (in) :: time
    integer (kind=i4), intent (in) :: nt, nx
    integer (kind=i4), intent (inout) :: errstat
    integer (kind=i2), dimension (1:nx,0:nt-1), intent (out) :: amfdiag

    call clim_get_climatology (cpt, pge_idx, climatology, cli_wgh_ozo_pro, &
                               cli_idx_ozo_pro, lat, lon, time, nt, nx, &
                               errstat, amfdiag)
  end subroutine

  SUBROUTINE vlidort_deallocate (errstat)
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return
    if (allocated(lut_alb)) then
      deallocate (lut_alb, lut_sza, lut_vza, lut_srf, lut_wav, lut_ozo, &
           lut_alt_lay, lut_alt_lev, lut_pre_lay, lut_pre_lev, lut_air_col, &
           lut_ozo_col, lut_tmp, lut_Sb, lut_I0, lut_I1, lut_I2, lut_Ir, lut_dI0, &
           lut_dI1, lut_dI2, stat=errstat)
    endif
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, "vlidort_deallocate failed", errstat)
      return
    endif
  END SUBROUTINE

  SUBROUTINE vlidort_allocate (errstat)
    IMPLICIT NONE
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    INTEGER   (KIND=i4) :: estat

    if (errstat /= 0) return

    ALLOCATE (lut_alb(alb_dim(1)), &
              lut_ozo(ozo_dim(1)), &
              lut_sza(sza_dim(1)), &
              lut_vza(vza_dim(1)), &
              lut_srf(srf_dim(1)), &
              lut_wav(wav_dim(1)), &
              lut_alt_lay(lay_dim(1)), &
              lut_alt_lev(lev_dim(1)), &
              lut_pre_lay(lay_dim(1)), &
              lut_pre_lev(lev_dim(1)), &
              lut_air_col(ozo_dim(1),lay_dim(1)), &
              lut_ozo_col(ozo_dim(1),lay_dim(1)), &
              lut_tmp(ozo_dim(1),lev_dim(1)), &
              lut_I0(ozo_dim(1),srf_dim(1),vza_dim(1),sza_dim(1)), &
              lut_I1(ozo_dim(1),srf_dim(1),vza_dim(1),sza_dim(1)), &
              lut_I2(ozo_dim(1),srf_dim(1),vza_dim(1),sza_dim(1)), &
              lut_Ir(ozo_dim(1),srf_dim(1),vza_dim(1),sza_dim(1)), &
              lut_Sb(ozo_dim(1),srf_dim(1)), &
              lut_dI0(ozo_dim(1),srf_dim(1),lay_dim(1),alb_dim(1),vza_dim(1),sza_dim(1)), &
              lut_dI1(ozo_dim(1),srf_dim(1),lay_dim(1),alb_dim(1),vza_dim(1),sza_dim(1)), &
              lut_dI2(ozo_dim(1),srf_dim(1),lay_dim(1),alb_dim(1),vza_dim(1),sza_dim(1)), &
              STAT=estat )
    if (estat /= 0) then
      call tell_error (tell_malloc_error, "vlidort_allocate: failed", errstat)
      return
    endif

  END SUBROUTINE vlidort_allocate


  SUBROUTINE read_vlidort (errstat)

    ! ====================================================
    ! This subroutine reads in the VLIDORT calculations to
    ! compute the Scattering Weights.
    ! ====================================================

    USE HDF5, ONLY: HID_T, SIZE_T, h5dopen_f, h5dget_space_f, &
      h5dread_f, h5sget_simple_extent_dims_f, h5open_f, h5tcopy_f, &
      h5tset_size_f, h5dclose_f, h5fopen_f, h5fclose_f, &
      H5F_ACC_RDONLY_F, H5T_NATIVE_CHARACTER, H5T_NATIVE_REAL

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: hdferr

    INTEGER(HID_T) :: fid, alb_did, ozo_did, sza_did, vza_did, srf_did, &
         wav_did, i0_did, i1_did, i2_did, ir_did, sb_did, air_lay_did, alt_lay_did, &
         alt_lev_did, ozo_lay_did, pre_lay_did, pre_lev_did, tmp_lev_did, &
         di0_did, di1_did, di2_did, dspace, ozo_datatype_id

    INTEGER(SIZE_T) :: size
    INTEGER(HSIZE_T), DIMENSION(2) :: tmp_2d_dim
    INTEGER(HSIZE_T), DIMENSION(4) :: tmp_4d_dim
    INTEGER(HSIZE_T), DIMENSION(6) :: tmp_6d_dim

    LOGICAL, SAVE :: h5inited = .FALSE.

    CHARACTER(LEN=MAX_STR_LEN)       :: filename
    CHARACTER(LEN=MAX_STR_LEN), SAVE :: cached_filename = ""

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    if (errstat /= 0) return

    filename = TRIM(ADJUSTL(OMSAO_wfamf_table_filename))
    if (filename .eq. cached_filename) then
      call tell_log (1, "read_vlidort: using cached data from "//trim(filename))
      return
    endif

    ! ---------------------------------
    ! Initialize hdf5 FORTRAN Interface
    ! ---------------------------------
    if (.NOT.h5inited) then
      CALL h5open_f(hdferr)
      h5inited = .TRUE.
    endif

    ! ------------------
    ! Dataset data types
    ! ------------------
    size = 4
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, ozo_datatype_id, hdferr)
    CALL h5tset_size_f(ozo_datatype_id, size, hdferr)

    ! ******************************************
    ! Find out the dimensions of the input file:
    ! ******************************************
    ! -------------------
    ! Opening input TABLE
    ! -------------------
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, hdferr)
    IF (hdferr .eq. -1) THEN
      call tell_error (tell_io_open_error, 'opening '//trim(filename), &
                       errstat)
      return
    END IF

    ! -------------
    ! Open datasets
    ! -------------
    CALL h5dopen_f(fid,'/Grid/Albedo', alb_did, hdferr)
    CALL h5dopen_f(fid,'/Grid/OZO', ozo_did, hdferr)
    CALL h5dopen_f(fid,'/Grid/Surface Pressure', srf_did, hdferr)
    CALL h5dopen_f(fid,'/Grid/SZA', sza_did, hdferr)
    CALL h5dopen_f(fid,'/Grid/VZA', vza_did, hdferr)
    CALL h5dopen_f(fid,'/Grid/Wavelength', wav_did,hdferr)

    CALL h5dopen_f(fid,'Profiles/Air Column Layer',air_lay_did,hdferr)
    CALL h5dopen_f(fid,'Profiles/Altitude Layer',alt_lay_did,hdferr)
    CALL h5dopen_f(fid,'Profiles/Altitude Level',alt_lev_did,hdferr)
    CALL h5dopen_f(fid,'Profiles/Ozone Column Layer',ozo_lay_did,hdferr)
    CALL h5dopen_f(fid,'Profiles/Pressure Layer',pre_lay_did,hdferr)
    CALL h5dopen_f(fid,'Profiles/Pressure Level',pre_lev_did,hdferr)
    CALL h5dopen_f(fid,'Profiles/Temperature Level',tmp_lev_did,hdferr)

    CALL h5dopen_f(fid,'/Intensity/I0', I0_did, hdferr)
    CALL h5dopen_f(fid,'/Intensity/I1', I1_did, hdferr)
    CALL h5dopen_f(fid,'/Intensity/I2', I2_did, hdferr)
    CALL h5dopen_f(fid,'/Intensity/Ir', Ir_did, hdferr)
    CALL h5dopen_f(fid,'/Intensity/Sb', Sb_did, hdferr)

    CALL h5dopen_f(fid,'/Scattering Weights/dI0', dI0_did, hdferr)
    CALL h5dopen_f(fid,'/Scattering Weights/dI1', dI1_did, hdferr)
    CALL h5dopen_f(fid,'/Scattering Weights/dI2', dI2_did, hdferr)

    ! -----------------------
    ! Find out the dimensions
    ! -----------------------
    CALL h5dget_space_f(alb_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, alb_dim, alb_maxdim, hdferr)
    CALL h5dget_space_f(ozo_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, ozo_dim, ozo_maxdim, hdferr)
    CALL h5dget_space_f(sza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, sza_dim, sza_maxdim, hdferr)
    CALL h5dget_space_f(srf_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, srf_dim, srf_maxdim, hdferr)
    CALL h5dget_space_f(vza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, vza_dim, vza_maxdim, hdferr)
    CALL h5dget_space_f(wav_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, wav_dim, wav_maxdim, hdferr)
    CALL h5dget_space_f(alt_lay_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, lay_dim, lay_maxdim, hdferr)
    CALL h5dget_space_f(alt_lev_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, lev_dim, lev_maxdim, hdferr)

    ! ---------------------------------------------------------------
    ! Allocate & initialize variables now that we have the dimensions
    ! ---------------------------------------------------------------
    CALL vlidort_allocate (errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_vlidort: allocate failed", &
                       errstat)
      return
    endif

    ! ----------------------------------------------------
    ! Read from the h5 file all these small size variables
    ! ----------------------------------------------------
    ! Grids: albedo, ozone profiles, sza, vza, surface pressure, wavelength
    CALL h5dread_f(alb_did, H5T_NATIVE_REAL, lut_alb, alb_dim, hdferr)
    CALL h5dread_f(ozo_did, ozo_datatype_id, lut_ozo, ozo_dim, hdferr)
    CALL h5dread_f(sza_did, H5T_NATIVE_REAL, lut_sza, sza_dim, hdferr)
    CALL h5dread_f(vza_did, H5T_NATIVE_REAL, lut_vza, vza_dim, hdferr)
    CALL h5dread_f(srf_did, H5T_NATIVE_REAL, lut_srf, srf_dim, hdferr)
    CALL h5dread_f(wav_did, H5T_NATIVE_REAL, lut_wav, wav_dim, hdferr)
    ! Altitude and pressure layer
    CALL h5dread_f(alt_lay_did, H5T_NATIVE_REAL, lut_alt_lay, lay_dim, hdferr)
    CALL h5dread_f(pre_lay_did, H5T_NATIVE_REAL, lut_pre_lay, lay_dim, hdferr)
    ! Altitude and pressure level
    CALL h5dread_f(alt_lev_did, H5T_NATIVE_REAL, lut_alt_lev, lev_dim, hdferr)
    CALL h5dread_f(pre_lev_did, H5T_NATIVE_REAL, lut_pre_lev, lev_dim, hdferr)
    ! Air Column, Ozone Column
    tmp_2d_dim(1) = ozo_dim(1); tmp_2d_dim(2) = lay_dim(1)
    CALL h5dread_f(air_lay_did, H5T_NATIVE_REAL, lut_air_col, tmp_2d_dim, hdferr)
    CALL h5dread_f(ozo_lay_did, H5T_NATIVE_REAL, lut_ozo_col, tmp_2d_dim, hdferr)
    ! Temperature
    tmp_2d_dim(1) = ozo_dim(1); tmp_2d_dim(2) = lev_dim(1)
    CALL h5dread_f(tmp_lev_did, H5T_NATIVE_REAL, lut_tmp, tmp_2d_dim, hdferr)
    ! Sb
    tmp_2d_dim(1) = ozo_dim(1); tmp_2d_dim(2) = srf_dim(1)
    CALL h5dread_f(sb_did, H5T_NATIVE_REAL, lut_Sb, tmp_2d_dim, hdferr)
    ! I0, I1, I2 & Ir
    tmp_4d_dim(1) = ozo_dim(1); tmp_4d_dim(2) = srf_dim(1)
    tmp_4d_dim(3) = vza_dim(1); tmp_4d_dim(4) = sza_dim(1)
    CALL h5dread_f(i0_did, H5T_NATIVE_REAL, lut_i0, tmp_4d_dim, hdferr)
    CALL h5dread_f(i1_did, H5T_NATIVE_REAL, lut_i1, tmp_4d_dim, hdferr)
    CALL h5dread_f(i2_did, H5T_NATIVE_REAL, lut_i2, tmp_4d_dim, hdferr)
    CALL h5dread_f(ir_did, H5T_NATIVE_REAL, lut_ir, tmp_4d_dim, hdferr)
    ! dI0, dI1 & dI2
    tmp_6d_dim(1) = ozo_dim(1); tmp_6d_dim(2) = srf_dim(1)
    tmp_6d_dim(3) = lay_dim(1); tmp_6d_dim(4) = alb_dim(1)
    tmp_6d_dim(5) = vza_dim(1); tmp_6d_dim(6) = sza_dim(1)
    CALL h5dread_f(di0_did, H5T_NATIVE_REAL, lut_di0, tmp_6d_dim, hdferr)
    CALL h5dread_f(di1_did, H5T_NATIVE_REAL, lut_di1, tmp_6d_dim, hdferr)
    CALL h5dread_f(di2_did, H5T_NATIVE_REAL, lut_di2, tmp_6d_dim, hdferr)

    ! --------------
    ! Close datasets
    ! --------------
    CALL h5dclose_f (alb_did, hdferr)
    CALL h5dclose_f (ozo_did, hdferr)
    CALL h5dclose_f (sza_did, hdferr)
    CALL h5dclose_f (vza_did, hdferr)
    CALL h5dclose_f (srf_did, hdferr)
    CALL h5dclose_f (wav_did, hdferr)

    CALL h5dclose_f (air_lay_did, hdferr)
    CALL h5dclose_f (alt_lay_did, hdferr)
    CALL h5dclose_f (alt_lev_did, hdferr)
    CALL h5dclose_f (ozo_lay_did, hdferr)
    CALL h5dclose_f (pre_lay_did, hdferr)
    CALL h5dclose_f (pre_lev_did, hdferr)
    CALL h5dclose_f (tmp_lev_did, hdferr)

    CALL h5dclose_f(I0_did, hdferr)
    CALL h5dclose_f(I1_did, hdferr)
    CALL h5dclose_f(I2_did, hdferr)
    CALL h5dclose_f(Ir_did, hdferr)
    CALL h5dclose_f(Sb_did, hdferr)

    CALL h5dclose_f(dI0_did, hdferr)
    CALL h5dclose_f(dI1_did, hdferr)
    CALL h5dclose_f(dI2_did, hdferr)

    ! ----------
    ! Close file
    ! ----------
    CALL h5fclose_f(fid, hdferr)

    cached_filename = filename

    errstat = hdferr

  END SUBROUTINE read_vlidort

  SUBROUTINE compute_scatt (cpt, nt, nx, time, albedo, sza, vza, saa, vaa, l2ctp, l2cfr, &
                            terrain_height, surface_pressure, cli_wgh_ozo_pro, cli_idx_ozo_pro, &
                            lat, lon, xtrange, amfdiag, scattw)

    USE OMSAO_linterpolation_module, ONLY: lininterpol, GetNode
    USE ezspline_interpolation, ONLY: ezspline_2d_interpolation
    use sao_pge_utils, only: calc_relaz_angle
    use clim_module

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    type (clim_pres_type), intent(in) :: cpt
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    REAL    (KIND=r8), DIMENSION (0:nt-1),      INTENT (IN) :: time
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza, saa, vaa, terrain_height, lat, lon
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: albedo, l2cfr
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,1:2), INTENT (IN) :: cli_wgh_ozo_pro
    INTEGER (KIND=i4), DIMENSION (1:nx,0:nt-1,1:2), INTENT (IN) :: cli_idx_ozo_pro
    integer (KIND=i4), dimension (0:nt-1,1:2),  intent (IN) :: xtrange

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (CmETA,1:nx,0:nt-1), INTENT (INOUT) :: scattw
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: l2ctp
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: surface_pressure

    ! ----------------
    ! output variables
    ! ----------------
    integer (KIND=i2), dimension (1:nx,0:nt-1), intent (out) :: amfdiag

    ! ---------------
    ! Local variables
    ! ---------------
    real (kind=r4), dimension(CmETA+1) :: eta_a, eta_b
    INTEGER (KIND=i4) :: ialb, ictp, ilay, isrf, isza, ivza, itime, ixtrack, &
         nsza, nvza, nalb, ncld_alb, nsrf, nctp, status
    character (len=72) :: logmsg

    ! LUT ozone profile variables
    INTEGER (KIND=i4), PARAMETER :: nozo = 2
    REAL (KIND=r8), DIMENSION(nozo) :: local_ozo_wgh
    INTEGER (KIND=i4), DIMENSION(nozo) :: local_ozo_idx

    ! Interpolation variables
    INTEGER (KIND=i4), DIMENSION(2) :: idx_sza, idx_vza, idx_alb, idx_cld_alb, idx_srf, idx_ctp
    REAL (KIND=r8) :: Radiance_cld, Radiance_clr, delta1, delta2, delta3
    REAL (KIND=r8), DIMENSION(:), ALLOCATABLE :: Sca_1D, Sca_1D_cloud
    REAL (KIND=r8), DIMENSION(:,:), ALLOCATABLE :: Sca_2D, Sca_2D_cloud
    REAL (KIND=r8), DIMENSION(:,:,:), ALLOCATABLE :: Rad_3D_clear, Rad_3D_cloud
    REAL (KIND=r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sca_5D_clear, Sca_5D_cloud
    REAL (KIND=r8) :: local_alb, local_sza, local_vza, local_srf, local_ctp, &
         local_cfr, local_raa, out_pre_lay
    real (kind=4) :: local_saa, local_vaa

    ! Error variables
    INTEGER (KIND=i4) :: locerrstat

    REAL (kind=8), PARAMETER :: d2r = 3.141592653589793d0/180.0  !! JED fix

    write(logmsg, '(a)')'Computing scattering weights...'
    call tell_log (1, logmsg)

    locerrstat = 0
    call clim_pres_eta (cpt, eta_a, eta_b, locerrstat)
    if (locerrstat /= 0) return
  
    ! ---------------
    ! Loop over lines
    ! ---------------
    DO itime = 0, nt-1
       ! --------------------------
       ! Loop over xtrack positions
       ! --------------------------
       DO ixtrack = 1, nx
          ! Only calculate scattering weights if we have geolocation,
          ! albedo, and cloud information. 
          if ( btest(amfdiag(ixtrack,itime),yn_amf_cor) .or. &
               btest(amfdiag(ixtrack,itime),yn_albedo) .or. &
               btest(amfdiag(ixtrack,itime),yn_cld) .or. &
               ixtrack < xtrange(itime,1) .or. ixtrack > xtrange(itime,2) ) then
             amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
             cycle
          end if

          ! ----------------------------------------------
          ! If this point is reached then scattw should be
          ! different from r8_missval and it needs to be
          ! initialized to 0.0 to work out the average
          ! ----------------------------------------------
          scattw(:,ixtrack,itime) = 0.0_r8

          ! -----------------------------------
          ! Fill up local values for this pixel
          ! -----------------------------------
          local_alb = albedo(ixtrack,itime)
          local_ctp = l2ctp(ixtrack,itime)
          local_cfr = l2cfr(ixtrack,itime)
          local_sza = REAL(sza(ixtrack,itime), KIND = r8)  ! JED fix
          local_vza = REAL(vza(ixtrack,itime), KIND = r8)  ! JED fix
          local_saa = saa(ixtrack,itime)
          local_vaa = vaa(ixtrack,itime)
          local_srf = REAL(terrain_height(ixtrack,itime), KIND = r8)
          local_ozo_wgh(1:2) = cli_wgh_ozo_pro(ixtrack,itime,1:2)
          local_ozo_idx(1:2) = cli_idx_ozo_pro(ixtrack,itime,1:2)

          ! ----------------------
          ! Relative azimuth angle
          ! ----------------------
          local_raa = calc_relaz_angle(local_saa,local_vaa)

          ! ----------------------------------------------
          ! Convert pixel terrain height to pressure using
          ! Xiong suggested to use pressure altitude:
          !  Z = -16 alog10 (P / Po) Z in km and P in hPa.
          ! ----------------------------------------------
          local_srf = 1013.0_r8 * (10.0_r8 ** (local_srf / 1000.0_r8 / (-16.0_r8)))


          ! Make sure that clouds are above or at the surface
          if ( local_ctp > local_srf ) then
             amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_adj_cld_pre)
             local_ctp = local_srf
             l2ctp(ixtrack,itime) = local_ctp
          end if
          ! FIXME???
          ! Make sure surface and local
          ! pressures are within LUT limits
          if ( local_srf > maxval(lut_srf) ) then
             local_srf = maxval(lut_srf)
             amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_adj_srf_pre)
          else if (local_srf < minval(lut_srf) ) then
             amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_adj_srf_pre)
             local_srf = minval(lut_srf)
          end if
          if ( local_ctp > maxval(lut_srf) ) then
             amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_adj_cld_pre)
             local_ctp = maxval(lut_srf)
             l2ctp(ixtrack,itime) = local_ctp
          else if (local_ctp < minval(lut_srf) ) then
             amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_adj_cld_pre)
             local_ctp = minval(lut_srf)
             l2ctp(ixtrack,itime) = local_ctp
          end if             

          ! FIXME!!! Assign surface pressure values. They should be comming from the climatology
          surface_pressure(ixtrack,itime) = real(local_srf,kind=r4)

          ! -----------------------------------------------
          ! Don't compute scattering weight if local_sza or
          ! local_vza are outside LUT limits
          ! -----------------------------------------------
          if ( local_sza > maxval(lut_sza) .or. local_vza > maxval(lut_vza) ) then
             amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
             cycle
          end if

          ! ----------
          ! Find nodes
          ! ----------
          ! SZA
          CALL GetNode(REAL(lut_sza,KIND=8),local_sza,idx_sza(1), 'Lower')
          IF (idx_sza(1) .EQ. -2) idx_sza(1) = 1
          IF (idx_sza(1) .EQ. -3) idx_sza(1) = INT(sza_dim(1),KIND=i4)-1
          CALL GetNode(REAL(lut_sza,KIND=8),local_sza,idx_sza(2), 'Upper')
          IF (idx_sza(2) .EQ. -2) idx_sza(2) = 2
          IF (idx_sza(2) .EQ. -3) idx_sza(2) = INT(sza_dim(1),KIND=i4)
          nsza = idx_sza(2)-idx_sza(1) + 1
          IF (nsza == 1) THEN
             IF (idx_sza(1) <= sza_dim(1)-1) THEN
                idx_sza(2) = idx_sza(1) + 1
             ELSE IF (idx_sza(1) == sza_dim(1)) THEN
                idx_sza(1) = INT(sza_dim(1),KIND=i4)-1
                idx_sza(2) = INT(sza_dim(1),KIND=i4)
             ENDIF
             nsza = idx_sza(2)-idx_sza(1)+1
          ENDIF
          ! VZA
          CALL GetNode(REAL(lut_vza,KIND=8),local_vza,idx_vza(1), 'Lower')
          IF (idx_vza(1) .EQ. -2) idx_vza(1) = 1
          IF (idx_vza(1) .EQ. -3) idx_vza(1) = INT(vza_dim(1),KIND=i4)-1
          CALL GetNode(REAL(lut_vza,KIND=8),local_vza,idx_vza(2), 'Upper')
          IF (idx_vza(2) .EQ. -2) idx_vza(2) = 2
          IF (idx_vza(2) .EQ. -3) idx_vza(2) = INT(vza_dim(1),KIND=i4)
          nvza = idx_vza(2)-idx_vza(1) + 1
          IF (nvza == 1) THEN
             IF (idx_vza(1) <= vza_dim(1)-1) THEN
                idx_vza(2) = idx_vza(1) + 1
             ELSE IF (idx_vza(1) == vza_dim(1)) THEN
                idx_vza(1) = INT(vza_dim(1),KIND=i4)-1
                idx_vza(2) = INT(vza_dim(1),KIND=i4)
             ENDIF
             nvza = idx_vza(2)-idx_vza(1)+1
          ENDIF
          ! ALBEDO
          CALL GetNode(REAL(lut_alb,KIND=8),local_alb,idx_alb(1), 'Lower')
          IF (idx_alb(1) .EQ. -2) idx_alb(1) = 1
          IF (idx_alb(1) .EQ. -3) idx_alb(1) = INT(alb_dim(1),KIND=i4)-1
          CALL GetNode(REAL(lut_alb,KIND=8),local_alb,idx_alb(2), 'Upper')
          IF (idx_alb(2) .EQ. -2) idx_alb(2) = 2
          IF (idx_alb(2) .EQ. -3) idx_alb(2) = INT(alb_dim(1),KIND=i4)
          nalb = idx_alb(2)-idx_alb(1) + 1
          IF (nalb == 1) THEN
             IF (idx_alb(1) <= alb_dim(1)-1) THEN
                idx_alb(2) = idx_alb(1) + 1
             ELSE IF (idx_alb(1) == alb_dim(1)) THEN
                idx_alb(1) = INT(alb_dim(1),KIND=i4)-1
                idx_alb(2) = INT(alb_dim(1),KIND=i4)
             ENDIF
             nalb = idx_alb(2)-idx_alb(1)+1
          ENDIF
          ! CLOUD ALBEDO
          CALL GetNode(REAL(lut_alb,KIND=8),amf_alb_cld,idx_cld_alb(1), 'Lower')
          IF (idx_cld_alb(1) .EQ. -2) idx_cld_alb(1) = 1
          IF (idx_cld_alb(1) .EQ. -3) idx_cld_alb(1) = INT(alb_dim(1),KIND=i4)-1
          CALL GetNode(REAL(lut_alb,KIND=8),amf_alb_cld,idx_cld_alb(2), 'Upper')
          IF (idx_cld_alb(2) .EQ. -2) idx_cld_alb(2) = 2
          IF (idx_cld_alb(2) .EQ. -3) idx_cld_alb(2) = INT(alb_dim(1),KIND=i4)
          ncld_alb = idx_cld_alb(2)-idx_cld_alb(1) + 1
          IF (ncld_alb == 1) THEN
             IF (idx_cld_alb(1) <= alb_dim(1)-1) THEN
                idx_cld_alb(2) = idx_cld_alb(1) + 1
             ELSE IF (idx_cld_alb(1) == alb_dim(1)) THEN
                idx_cld_alb(1) = INT(alb_dim(1),KIND=i4)-1
                idx_cld_alb(2) = INT(alb_dim(1),KIND=i4)
             ENDIF
             ncld_alb = idx_cld_alb(2)-idx_cld_alb(1)+1
          ENDIF
          ! SURFACE PRESSURE
          CALL GetNode(REAL(lut_srf,KIND=8),local_srf,idx_srf(1), 'Lower')
          IF (idx_srf(1) .EQ. -2) idx_srf(1) = 1
          IF (idx_srf(1) .EQ. -3) idx_srf(1) = INT(srf_dim(1),KIND=i4)-1
          CALL GetNode(REAL(lut_srf,KIND=8),local_srf,idx_srf(2), 'Upper')
          IF (idx_srf(2) .EQ. -2) idx_srf(2) = 2
          IF (idx_srf(2) .EQ. -3) idx_srf(2) = INT(srf_dim(1),KIND=i4)
          nsrf = idx_srf(2)-idx_srf(1) + 1
          IF (nsrf == 1) THEN
             IF (idx_srf(1) <= srf_dim(1)-1) THEN
                idx_srf(2) = idx_srf(1) + 1
             ELSE IF (idx_srf(1) == srf_dim(1)) THEN
                idx_srf(1) = INT(srf_dim(1),KIND=i4)-1
                idx_srf(2) = INT(srf_dim(1),KIND=i4)
             ENDIF
             nsrf = idx_srf(2)-idx_srf(1)+1
          ENDIF
          ! CLOUD PRESSURE
          CALL GetNode(REAL(lut_srf,KIND=8),local_ctp,idx_ctp(1), 'Lower')
          IF (idx_ctp(1) .EQ. -2) idx_ctp(1) = 1
          IF (idx_ctp(1) .EQ. -3) idx_ctp(1) = INT(srf_dim(1),KIND=i4)-1
          CALL GetNode(REAL(lut_srf,KIND=8),local_ctp,idx_ctp(2), 'Upper')
          IF (idx_ctp(2) .EQ. -2) idx_ctp(2) = 2
          IF (idx_ctp(2) .EQ. -3) idx_ctp(2) = INT(srf_dim(1),KIND=i4)
          nctp = idx_ctp(2)-idx_ctp(1) + 1
          IF (nctp == 1) THEN
             IF (idx_ctp(1) <= srf_dim(1)-1) THEN
                idx_ctp(2) = idx_ctp(1) + 1
             ELSE IF (idx_ctp(1) == srf_dim(1)) THEN
                idx_ctp(1) = INT(srf_dim(1),KIND=i4)-1
                idx_ctp(2) = INT(srf_dim(1),KIND=i4)
             ENDIF
             nctp = idx_ctp(2)-idx_ctp(1)+1
          ENDIF

          ! -----------------------------------------------
          ! Allocate and initialize interpolation variables
          ! -----------------------------------------------
          ALLOCATE(Sca_1D(lay_dim(1)), &
               Sca_1D_cloud(lay_dim(1)), &
               Sca_2D(nsrf,lay_dim(1)), &
               Sca_2D_cloud(nctp,lay_dim(1)), &
               Rad_3D_clear(nsrf,nvza,nsza), &
               Rad_3D_cloud(nctp,nvza,nsza), &
               Sca_5D_clear(nsrf,lay_dim(1),nalb,nvza,nsza), &
               Sca_5D_cloud(nctp,lay_dim(1),ncld_alb,nvza,nsza), STAT=locerrstat)
          if (locerrstat /= 0) then
             call tell_error (tell_malloc_error, "compute_scatt:  allocate failed", &
                  status)
             return
          endif
          Rad_3D_clear = 0.0; Rad_3D_cloud = 0.0
          Sca_1D = 0.0; Sca_1D_cloud = 0.0; Sca_2D = 0.0; Sca_2D_cloud = 0.0
          Sca_5D_clear = 0.0; Sca_5D_cloud = 0.

          ! -------------------------------------------------------------
          ! First compute back from the parametrization on RAA and albedo
          ! -------------------------------------------------------------
          DO isza = idx_sza(1), idx_sza(2)
             DO ivza = idx_vza(1), idx_vza(2)
                DO isrf = idx_srf(1), idx_srf(2)
                   Rad_3D_clear(isrf-idx_srf(1)+1,&
                        ivza-idx_vza(1)+1,&
                        isza-idx_sza(1)+1) = REAL(lut_I0(local_ozo_idx(1),isrf,ivza,isza) + &
                        lut_I1(local_ozo_idx(1),isrf,ivza,isza) * COS(local_raa*d2r) + &
                        lut_I2(local_ozo_idx(1),isrf,ivza,isza) * COS(2.0*local_raa*d2r) + &
                        ( lut_Ir(local_ozo_idx(1),isrf,ivza,isza) * local_alb / &
                        (1.0 - local_alb * lut_Sb(local_ozo_idx(1),isrf))),KIND=8) &
                        * local_ozo_wgh(1) + REAL(lut_I0(local_ozo_idx(2),isrf,ivza,isza) + &
                        lut_I1(local_ozo_idx(2),isrf,ivza,isza) * COS(local_raa*d2r) + &
                        lut_I2(local_ozo_idx(2),isrf,ivza,isza) * COS(2.0*local_raa*d2r) + &
                        ( lut_Ir(local_ozo_idx(2),isrf,ivza,isza) * local_alb / &
                        (1.0 - local_alb * lut_Sb(local_ozo_idx(2),isrf))),KIND=8) &
                        * local_ozo_wgh(2)
                   DO ialb = idx_alb(1), idx_alb(2)
                      DO ilay = 1, INT(lay_dim(1),KIND=4)
                         Sca_5D_clear(isrf-idx_srf(1)+1, ilay, ialb-idx_alb(1)+1, &
                              ivza-idx_vza(1)+1, isza-idx_sza(1)+1) = REAL(&
                              lut_dI0(local_ozo_idx(1),isrf,ilay,ialb,ivza,isza) + &
                              lut_dI1(local_ozo_idx(1),isrf,ilay,ialb,ivza,isza) * COS(local_raa*d2r) + &
                              lut_dI2(local_ozo_idx(1),isrf,ilay,ialb,ivza,isza) * COS(2.0*local_raa*d2r), KIND=8) &
                              * local_ozo_wgh(1) + REAL(lut_dI0(local_ozo_idx(2),isrf,ilay,ialb,ivza,isza) + &
                              lut_dI1(local_ozo_idx(2),isrf,ilay,ialb,ivza,isza) * COS(local_raa*d2r) + &
                              lut_dI2(local_ozo_idx(2),isrf,ilay,ialb,ivza,isza) * COS(2.0*local_raa*d2r), KIND=8) &
                              * local_ozo_wgh(2)

                      END DO
                   END DO
                END DO
                DO ictp = idx_ctp(1), idx_ctp(2)
                   Rad_3D_cloud(ictp-idx_ctp(1)+1,&
                        ivza-idx_vza(1)+1,&
                        isza-idx_sza(1)+1) = REAL(lut_I0(local_ozo_idx(1),ictp,ivza,isza) + &
                        lut_I1(local_ozo_idx(1),ictp,ivza,isza) * COS(local_raa*d2r) + &
                        lut_I2(local_ozo_idx(1),ictp,ivza,isza) * COS(2.0*local_raa*d2r) + &
                        ( lut_Ir(local_ozo_idx(1),ictp,ivza,isza) * amf_alb_cld / &
                        (1.0 - amf_alb_cld * lut_Sb(local_ozo_idx(1),ictp))),KIND=8) * &
                        local_ozo_wgh(1) +  REAL(lut_I0(local_ozo_idx(2),ictp,ivza,isza) + &
                        lut_I1(local_ozo_idx(2),ictp,ivza,isza) * COS(local_raa*d2r) + &
                        lut_I2(local_ozo_idx(2),ictp,ivza,isza) * COS(2.0*local_raa*d2r) + &
                        ( lut_Ir(local_ozo_idx(2),ictp,ivza,isza) * amf_alb_cld / &
                        (1.0 - amf_alb_cld * lut_Sb(local_ozo_idx(2),ictp))),KIND=8) * &
                        local_ozo_wgh(2)
                   DO ialb = idx_cld_alb(1), idx_cld_alb(2)
                      DO ilay = 1, INT(lay_dim(1),KIND=4)
                         Sca_5D_cloud(ictp-idx_ctp(1)+1, ilay, ialb-idx_cld_alb(1)+1, &
                              ivza-idx_vza(1)+1, isza-idx_sza(1)+1) = REAL(&
                              lut_dI0(local_ozo_idx(1),ictp,ilay,ialb,ivza,isza) + &
                              lut_dI1(local_ozo_idx(1),ictp,ilay,ialb,ivza,isza) * COS(local_raa*d2r) + &
                              lut_dI2(local_ozo_idx(1),ictp,ilay,ialb,ivza,isza) * COS(2.0*local_raa*d2r),KIND=8) &
                              * local_ozo_wgh(1) + REAL(lut_dI0(local_ozo_idx(2),ictp,ilay,ialb,ivza,isza) + &
                              lut_dI1(local_ozo_idx(2),ictp,ilay,ialb,ivza,isza) * COS(local_raa*d2r) + &
                              lut_dI2(local_ozo_idx(2),ictp,ilay,ialb,ivza,isza) * COS(2.0*local_raa*d2r),KIND=8) &
                              * local_ozo_wgh(2)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          ! Radiance (perform linear interpolation on srf, vza, and sza)
          Radiance_clr = linInterpol(nsrf,nvza,nsza, &
               REAL(lut_srf(idx_srf(1):idx_srf(2)),KIND=8), &
               REAL(SIN(lut_vza(idx_vza(1):idx_vza(2))*d2r),KIND=8), &
               REAL(SIN(lut_sza(idx_sza(1):idx_sza(2))*d2r),KIND=8), &
               Rad_3D_clear, local_srf, SIN(local_vza*d2r), SIN(local_sza*d2r), status=status)
          IF ( status /= 0 ) THEN
            amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
            write(logmsg, '(a55,i4,i4)') &
                 "compute_scatt: Radiance_clr interpol failed at ", &
                 ixtrack,itime
             call tell_log (1,logmsg)
             goto 999
          END IF
          Radiance_cld = linInterpol(nctp,nvza,nsza, &
               REAL(lut_srf(idx_srf(1):idx_srf(2)),KIND=8), &
               REAL(SIN(lut_vza(idx_vza(1):idx_vza(2))*d2r),KIND=8), &
               REAL(SIN(lut_sza(idx_sza(1):idx_sza(2))*d2r),KIND=8), &
               Rad_3D_cloud, local_srf, SIN(local_vza*d2r), SIN(local_sza*d2r), status=status)
          IF ( status /= 0 ) THEN
            amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
            write(logmsg, '(a55,i4,i4)') &
                 "compute_scatt: Radiance_cld interpol failed at ", &
                 ixtrack,itime
             call tell_log (1,logmsg)
             goto 999
          END IF

          ! Scattering Weights linear interpolation on alb, vza, sza
          DO isrf = 1, nsrf
             DO ilay = 1, INT(lay_dim(1),KIND=4)
                Sca_2D(isrf,ilay) = linInterpol(nalb,nvza,nsza, &
                     REAL(lut_alb(idx_alb(1):idx_alb(2)),KIND=8), &
                     REAL(SIN(lut_vza(idx_vza(1):idx_vza(2))*d2r),KIND=8), &
                     REAL(SIN(lut_sza(idx_sza(1):idx_sza(2))*d2r),KIND=8), &
                     Sca_5D_clear(isrf,ilay,1:nalb,1:nvza,1:nsza), &
                     local_alb, SIN(local_vza*d2r), SIN(local_sza*d2r), status=status)
                IF ( status /= 0 ) THEN
                   amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
                   write(logmsg, '(a55,i4,i4)') &
                        "compute_scatt: Sca_2D interpol failed at ", &
                        ixtrack,itime
                   call tell_log (1,logmsg)
                   goto 999
                END IF
             END DO
          END DO
          ! Scattering Weights, linear interpolation on surface pressure
          DO ilay = 1, INT(lay_dim(1),KIND=4)
             IF (REAL(lut_pre_lay(ilay),KIND=8) .GT. local_srf) THEN
                Sca_1D(ilay) = 0.0
                cycle
             ENDIF
             IF (nsrf .EQ. 2 .AND. Sca_2D(1,ilay) .LE. 0 .AND. Sca_2D(2,ilay) .GT. 0) THEN
                delta1=(Sca_2D(2,ilay-2)-Sca_2D(2,ilay-1)) / &
                     (LOG(lut_pre_lay(ilay-2)) - LOG(lut_pre_lay(ilay-1)))
                delta2=(Sca_1D(ilay-2)-Sca_1D(ilay-1)) / &
                     (LOG(lut_pre_lay(ilay-2)) - LOG(lut_pre_lay(ilay-1)))
                delta3=(Sca_2D(2,ilay-1)-Sca_2D(2,ilay)) / &
                     (LOG(lut_pre_lay(ilay-1)) - LOG(lut_pre_lay(ilay)))
                Sca_1D(ilay) = &
                     Sca_1d(ilay-1) - delta3*delta2/delta1 * &
                     (LOG(lut_pre_lay(ilay-1)) - LOG(lut_pre_lay(ilay)))
                cycle
             END IF
             Sca_1D(ilay) = linInterpol(nsrf, REAL(lut_srf(idx_srf(1):idx_srf(2)),KIND=8), &
                  Sca_2D(1:nsrf,ilay), REAL(local_srf,KIND=8), status=status)
             IF ( status /= 0 ) THEN
               amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
               write(logmsg, '(a55,i4,i4)') &
                    "compute_scatt: Sca_1D interpol failed at ", &
                    ixtrack,itime
               call tell_log (1,logmsg)
               goto 999
             END IF
          END DO

          ! Cloudy scattering weights linear interpolation on amf_alb_cld, vza, sza
          DO ictp = 1, nctp
             DO ilay = 1, INT(lay_dim(1),KIND=4)
                Sca_2D_cloud(ictp,ilay) = linInterpol(ncld_alb,nvza,nsza, &
                     REAL(lut_alb(idx_cld_alb(1):idx_cld_alb(2)),KIND=8), &
                     REAL(SIN(lut_vza(idx_vza(1):idx_vza(2))*d2r),KIND=8), &
                     REAL(SIN(lut_sza(idx_sza(1):idx_sza(2))*d2r),KIND=8), &
                     Sca_5D_cloud(ictp,ilay,1:ncld_alb,1:nvza,1:nsza), &
                     amf_alb_cld, SIN(local_vza*d2r), SIN(local_sza*d2r), status=status)
                IF ( status /= 0 ) THEN
                  amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
                  write(logmsg, '(a55,i4,i4)') &
                       "compute_scatt: Sca_2D_cloud interpol failed at ", &
                       ixtrack,itime
                  call tell_log (1,logmsg)
                  goto 999
                END IF
             END DO
          END DO
          ! Cloudy scattering weights linear interpolation on cloud pressure
          DO ilay = 1, INT(lay_dim(1),KIND=4)
             IF (REAL(lut_pre_lay(ilay),KIND=8) .GT. local_ctp) THEN
                Sca_1D_cloud(ilay) = 0.0
                cycle
             ENDIF
             IF (nctp .EQ. 2 .AND. Sca_2D_cloud(1,ilay) .LE. 0 .AND. Sca_2D_cloud(2,ilay) .GT. 0) THEN
                delta1=(Sca_2D_cloud(2,ilay-2)-Sca_2D_cloud(2,ilay-1)) / &
                     (LOG(lut_pre_lay(ilay-2)) - LOG(lut_pre_lay(ilay-1)))
                delta2=(Sca_1D_cloud(ilay-2)-Sca_1D_cloud(ilay-1)) / &
                     (LOG(lut_pre_lay(ilay-2)) - LOG(lut_pre_lay(ilay-1)))
                delta3=(Sca_2D_cloud(2,ilay-1)-Sca_2D_cloud(2,ilay)) / &
                     (LOG(lut_pre_lay(ilay-1)) - LOG(lut_pre_lay(ilay)))
                Sca_1D_cloud(ilay) = &
                     Sca_1d_cloud(ilay-1) - delta3*delta2/delta1 * &
                     (LOG(lut_pre_lay(ilay-1)) - LOG(lut_pre_lay(ilay)))
                cycle
             END IF
             Sca_1D_cloud(ilay) = linInterpol(nctp, REAL(lut_srf(idx_ctp(1):idx_ctp(2)),KIND=8), &
                  Sca_2D_cloud(1:nctp,ilay), REAL(local_ctp,KIND=8), status=status)
             IF ( status /= 0 ) THEN
               amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
               write(logmsg, '(a55,i4,i4)') &
                    "compute_scatt: Sca_1D_cloud interpol failed at ", &
                    ixtrack,itime
               call tell_log (1,logmsg)
               goto 999
             END IF
          END DO

          ! Convert effective cloud fraction to radiance cloud fraction
          local_cfr = local_cfr * Radiance_cld / ( local_cfr * Radiance_cld + (1.0 - local_cfr) * Radiance_clr)

          ! ---------------------------------------------------------------------------------
          ! Boersma et al. 2011 AMT, 4, 2011
          ! Cloud radiance fraction: Crf= Cfr * Icl / Ir
          !  We define Icl = From the Vlidort calculation, see above
          !            Icr = From the Vlidort calculation, see above
          !             Ir = Cfr * Icl + (1 - Cfr) * Icr (Total pixel radiance)
          !
          ! Now the scattering weights become w = crf * scatt_cloud + (1 - crf) * scatt_clear
          !  We add the scattweights calculated in the previous wavelengths.
          ! --------------------------------------------------------------------------------
          DO ilay = 1, INT(lay_dim(1),KIND=4)
             Sca_1D(ilay) = ( Sca_1D_cloud(ilay) * local_cfr + Sca_1D(ilay) * (1.0 - local_cfr) )
             IF (Sca_1D(ilay) .LT. 0.0) Sca_1D(ilay) = 0.0
          END DO

          ! --------------------------------------------------------------
          ! Interpolate Sca_1D from lut_pre_lay grid to the one defined by
          ! local_srf (from climatology or L1 file), eta_a and eta_b
          ! --------------------------------------------------------------
          DO ilay = 1, CmETA
             out_pre_lay = (( real(eta_a(ilay),kind=r8) + &
                              local_srf * real(eta_b(ilay),kind=r8)  ) + &
                            ( real(eta_a(ilay+1),kind=r8) + &
                              local_srf * real(eta_b(ilay+1),kind=r8) )) / 2.0
             IF ( (out_pre_lay > MAXVAL(lut_pre_lay)) .OR. (out_pre_lay < MINVAL(lut_pre_lay)) ) THEN
                scattw(ilay,ixtrack,itime) = 0.0
                cycle
             ENDIF
             scattw(ilay,ixtrack,itime) = linInterpol( (INT(lay_dim(1),KIND=i4)), REAL(LOG(lut_pre_lay),KIND=r8), &
                  Sca_1D, LOG(out_pre_lay), status=status)
             IF ( status /= 0 ) THEN
               amfdiag(ixtrack,itime) = ibset(amfdiag(ixtrack,itime),yn_sca)
               write(logmsg, '(a55,i4,i4)') &
                    "compute_scatt: output press grid interpol failed at ", &
                    ixtrack,itime
               call tell_log (1,logmsg)
               goto 999
             END IF
          END DO

 999      continue
          if (allocated(Sca_1D)) then
            DEALLOCATE(Rad_3D_clear, Rad_3D_cloud, Sca_5D_clear, &
                 Sca_5D_cloud, Sca_2D, Sca_2D_cloud, &
                 Sca_1D, Sca_1D_cloud, STAT=locerrstat)
            if (locerrstat /= 0) then
              call tell_error (tell_malloc_error, "compute_scatt:  de-allocate failed", &
                   status)
              return
            endif
          endif

          !  Set non-physical entries to zero.
          WHERE ( scattw(1:CmETA,ixtrack,itime) < 0.0_r8 )
             scattw(1:CmETA,ixtrack,itime) = 0.0_r8
          END WHERE

       END DO ! End loop xtrack

    END DO ! End loop lines

  END SUBROUTINE COMPUTE_SCATT

  SUBROUTINE compute_amf (cpt, nt, nx, CmETA, climatology, &
      scattw, saoamf, stratospheric_amf, tropospheric_amf, surface_pressure, &
      tropopause_pressure, lat, lon, amfdiag, errstat)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_char, c_null_ptr, c_associated
    use ctrlvars, only: yn_stratrop
    use met_module
    use clim_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    type (clim_pres_type), intent(in) :: cpt
    INTEGER (KIND=i4), INTENT(IN) :: nt, nx, CmETA
    REAL (KIND=r8), DIMENSION(CmETA,1:nx,0:nt-1), INTENT(IN) :: climatology, scattw
    INTEGER (KIND=i2), DIMENSION(1:nx,0:nt-1), INTENT(out) :: amfdiag
    REAL (KIND=r4), DIMENSION(1:nx,0:nt-1), INTENT(IN) :: surface_pressure
    real (kind=r4), dimension(1:nx,0:nt-1), intent(in) :: lat, lon

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT(INOUT) :: saoamf, stratospheric_amf, tropospheric_amf
    REAL (KIND=r4), DIMENSION(1:nx,0:nt-1), INTENT(INOUT) :: tropopause_pressure
    INTEGER (KIND=i4), INTENT(INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: ixtrack, itimes, ilay, tropopause_idx, imet, status
    REAL (KIND=r8), DIMENSION(:), ALLOCATABLE :: pressure_grid, temperature_profile, alpha
    real (kind=r4), dimension(CmETA+1) :: eta_a, eta_b
    real (kind=r4) :: lon_f, lat_f, ptrop
    real (kind=r4), dimension(CmETA) :: isobar_f, temp_on_isobar_f
    logical :: have_synthetic_met_data
    character (len=256) :: errmsg
    type (synth_met_type) :: smt

    real (kind=r8), parameter :: amf_magic_temperature_bucsela = 220.0_r8

    ! met_flags bit definitions:
    ! bit 0 set => read surface pressure
    ! bit 1 set => read tropopause pressure
    ! bit 2 set => read temperature vs height
    integer, parameter :: met_flags = 7

    type (c_ptr) :: met = c_null_ptr

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=11), PARAMETER :: modulename = 'compute_amf'

    if (errstat /= 0) return

    call clim_pres_eta (cpt, eta_a, eta_b, errstat)
    if (errstat /= 0) return

    if (0 /= index (OMSAO_meteorology_filename(1), '.nc', .true.)) then
      have_synthetic_met_data = .true.
      met = c_null_ptr
      call open_synth_met_data (smt, trim(OMSAO_meteorology_filename(1)), errstat)
      if (errstat /= 0) then
        call tell_error (tell_runtime_error, "error opening synthetic met data", &
                         errstat)
        return
      endif
    else
      have_synthetic_met_data = .false.
      met = met_list_new (met_flags)
      if (.not. c_associated(met)) then
        call tell_error (tell_runtime_error, "met_list_new returned NULL", &
                         errstat)
        return
      endif
      do imet=1, num_met_luns

        if (0 /= index(OMSAO_meteorology_filename(imet), 'grib2', .true.)) then
          status = met_list_add_file (met, trim(OMSAO_meteorology_filename(imet))//c_null_char)
          if (status /= 0) then
            call met_list_free (met)
            call tell_error (tell_runtime_error, &
                             "reading: "//trim(OMSAO_meteorology_filename(imet)), &
                             errstat)
            return
          endif
        endif

      enddo
    endif

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    DO itimes  = 0, nt-1 ! Swath lines loop
      DO ixtrack = 1, nx ! Xtrack pixel loop

        ! -------------------------------------------
        ! We can only calculate AMFs if we have both,
        ! scattering weights and gas climatology
        ! -------------------------------------------
        IF ( btest(amfdiag(ixtrack,itimes),yn_gas_cli) .or. btest(amfdiag(ixtrack,itimes),yn_sca) ) cycle

        ! ---------------------------------------------------
        ! Read tropopause pressure from met forecast file
        ! ---------------------------------------------------
        IF (yn_stratrop) THEN
           ! Allocate pressure_grid and temperature vertical profile
           ALLOCATE(pressure_grid(1:CmETA),temperature_profile(1:CmETA), &
                alpha(1:CmETA))
           ! Work out pressure_grid
           DO ilay = 1, CmETA
              pressure_grid(ilay) = (( real(eta_a(ilay),kind=r8) + &
                                       surface_pressure(ixtrack,itimes) * real(eta_b(ilay),kind=r8)  ) + &
                                     ( real(eta_a(ilay+1),kind=r8) + & 
                                       surface_pressure(ixtrack,itimes) * real(eta_b(ilay+1),kind=r8) )) / 2.0
           END DO

           ! Get tropopause pressure and temperature profile
           if (have_synthetic_met_data) then
             call read_synth_met_data(smt, lat(ixtrack,itimes), lon(ixtrack,itimes), &
                                      tropopause_pressure(ixtrack,itimes), errstat, &
                                      pprof = pressure_grid, tprof = temperature_profile)
             ! If any pressure grid values are out of range, the interpolated temperature
             ! will be NaN.  Replace such temperatures with the magic buscela temp.
             do ilay=1,CmETA
               if (isnan(temperature_profile(ilay))) then
                 temperature_profile(ilay) = amf_magic_temperature_bucsela
               endif
             enddo
           else
             lon_f = real (lon(ixtrack,itimes), kind=r4)
             lat_f = real (lat(ixtrack,itimes), kind=r4)
             isobar_f(1:CmETA) = real (pressure_grid(1:CmETA), kind=r4)

             call met_list_interp_f (met, lon_f, lat_f, errstat, &
                                     ptrop=ptrop, isobars=isobar_f, &
                                     temp_on_isobar=temp_on_isobar_f)
             if (errstat /= 0) then
               write(errmsg, *)'interpolating forecast for lon=',lon_f,' lat=',lat_f
               call tell_error (tell_runtime_error, errmsg, errstat)
               call met_list_free (met)
               return
             endif

             ! FIXME? Where the climatology pressure grid is not covered by the
             ! forecast pressure grid, the interpolated temperature will be NaN.
             ! For pressures with unknown temperature, we'll just assume the
             ! magic/fiducial temperature.
             do ilay = 1, CmETA
               if (isnan(temp_on_isobar_f(ilay))) then
                 temp_on_isobar_f(ilay) = real (amf_magic_temperature_bucsela, kind=r4)
               endif
             enddo

             tropopause_pressure(ixtrack,itimes) = ptrop
             temperature_profile(1:CmETA) = real (temp_on_isobar_f, kind=r8)
           endif

           ! Find which layer is closer to the tropopause.
           tropopause_idx = MINLOC(ABS(pressure_grid-REAL(tropopause_pressure(ixtrack,itimes),KIND=r4)),1)

           ! Compute stratospheric and tropospheric AMFs following Bucsela et al., 2013
           ! DOI:10.5194/amt-6-2607-2013
           ! Apply temperature correction factor alpha(p) = 1-0.003 [T(p)-T0] with T0 .EQ. 220K
           ! EJOS adding a test for zero in climatology to avoid NaN AMFs
           alpha = 1.0_r8-0.003_r8*(temperature_profile-amf_magic_temperature_bucsela)
           if (SUM(climatology(1:tropopause_idx,ixtrack,itimes)).eq.0) then
             tropospheric_amf(ixtrack,itimes) = 0.0d0
           else
             tropospheric_amf(ixtrack,itimes) = SUM(scattw(1:tropopause_idx,ixtrack, itimes) * &
                  climatology(1:tropopause_idx,ixtrack,itimes) * alpha(1:tropopause_idx))     / &
                  SUM(climatology(1:tropopause_idx,ixtrack,itimes))
           endif
           if (SUM(climatology(tropopause_idx+1:CmETA,ixtrack,itimes)).eq.0) then
             stratospheric_amf(ixtrack,itimes) = 0.0d0
           else
             stratospheric_amf(ixtrack,itimes) = SUM(scattw(tropopause_idx+1:CmETA,ixtrack, itimes) * &
                  climatology(tropopause_idx+1:CmETA,ixtrack,itimes) * alpha(tropopause_idx+1:CmETA) ) / &
                  SUM(climatology(tropopause_idx+1:CmETA,ixtrack,itimes))
           endif
           DEALLOCATE(pressure_grid,temperature_profile,alpha)
        END IF

        ! -------------------------
        ! Finally work out the AMFs
        ! -------------------------
        saoamf(ixtrack,itimes) = SUM(scattw(1:CmETA,ixtrack, itimes) * &
             climatology(1:CmETA,ixtrack,itimes))     / &
             SUM(climatology(1:CmETA,ixtrack,itimes))

        ! --------------------------------------------------------------------
        ! Unset amfdiag pixel 1 to indicate molecular instead of geometric amf
        ! --------------------------------------------------------------------
        amfdiag(ixtrack,itimes) = ibclr(amfdiag(ixtrack,itimes),yn_amf_geo)

      END DO ! Finish xtrack pixel loop
    END DO ! Finish

    if (c_associated(met)) then
      call met_list_free (met)
    else if (have_synthetic_met_data) then
      call close_synth_met_data (smt, errstat)
    endif

  END SUBROUTINE compute_amf

  subroutine read_climatology_dimensions (errstat)
    use clim_module
    implicit none
    integer, intent(inout) :: errstat

    integer :: nz

    if (errstat /= 0) return

    call clim_query_nz (nz, errstat)
    if ( errstat /= 0) return
    CmETA = nz - 1

  end subroutine read_climatology_dimensions

END MODULE OMSAO_wfamf_module

