MODULE OMSAO_wfamf_module

  ! ====================================================================
  ! This module defines variables associated with the wavelength depende
  ! nt AMF calculations and contains necessary subroutines to read files
  ! and calculate them
  ! ====================================================================
  USE OMSAO_precision_module, ONLY: i2, i4, r8, C_LONG, r4
  USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, i2_missval, i4_missval, r4_missval, r8_missval
  use tell_module
  USE OMSAO_he5_module, ONLY: pge_swath_id, &
    he5_start_4d, he5_edge_4d, he5_stride_4d, &
    he5_start_3d, he5_edge_3d, he5_stride_3d, &
    he5_start_2d, he5_edge_2d, he5_stride_2d, &
    he5_start_1d, he5_edge_1d, he5_stride_1d
  USE HDF5, ONLY: HSIZE_T
  IMPLICIT NONE
  private

  public read_climatology_dimensions, amf_calculation_bis, &
    wfamf_deallocate

  ! ---------
  ! PCF stuff
  ! ---------
  INTEGER(KIND=i4), PARAMETER, public :: wfamf_table_lun = 700250
  INTEGER(KIND=i4), PARAMETER, public :: climatology_lun = 700270
  integer(kind=i4), parameter, public :: meteorology_lun = 700290
  CHARACTER(LEN=MAX_STR_LEN), public  :: OMSAO_wfamf_table_filename
  CHARACTER(LEN=MAX_STR_LEN), public  :: OMSAO_climatology_filename
  CHARACTER(LEN=MAX_STR_LEN), public  :: OMSAO_meteorology_filename

  ! -----------------------------
  ! Dimensions of the climatology
  ! -----------------------------
  INTEGER (KIND=i4), public :: Cmlat, Cmlon, CmETA, CmHRS

  ! ====================================================================
  ! Wavelength dependent AMF factor specific variables
  ! ====================================================================
  LOGICAL, public       :: yn_amf_wfmod
  INTEGER, public       :: amf_wfmod_idx
  REAL(KIND=r8), public :: amf_alb_lnd, amf_alb_sno, amf_wvl, amf_wvl2, amf_alb_cld, amf_max_sza

  ! ---------------------------------
  ! GMAO GEOS-5 hybrid grid Ap and Bp
  ! ---------------------------------
  REAL(KIND=r8), DIMENSION(48), PARAMETER :: Ap=(/0.000000E+00, 4.804826E-02, 6.593752E+00, 1.313480E+01, &
       1.961311E+01, 2.609201E+01, 3.257081E+01, 3.898201E+01, &
       4.533901E+01, 5.169611E+01, 5.805321E+01, 6.436264E+01, &
       7.062198E+01, 7.883422E+01, 8.909992E+01, 9.936521E+01, &
       1.091817E+02, 1.189586E+02, 1.286959E+02, 1.429100E+02, &
       1.562600E+02, 1.696090E+02, 1.816190E+02, 1.930970E+02, &
       2.032590E+02, 2.121500E+02, 2.187760E+02, 2.238980E+02, &
       2.243630E+02, 2.168650E+02, 2.011920E+02, 1.769300E+02, &
       1.503930E+02, 1.278370E+02, 1.086630E+02, 9.236572E+01, &
       7.851231E+01, 5.638791E+01, 4.017541E+01, 2.836781E+01, &
       1.979160E+01, 9.292942E+00, 4.076571E+00, 1.650790E+00, &
       6.167791E-01, 2.113490E-01, 6.600001E-02, 1.000000E-02/)
  REAL(KIND=r8), DIMENSION(48), PARAMETER :: Bp=(/1.000000E+00, 9.849520E-01, 9.634060E-01, 9.418650E-01, &
       9.203870E-01, 8.989080E-01, 8.774290E-01, 8.560180E-01, &
       8.346609E-01, 8.133039E-01, 7.919469E-01, 7.706375E-01, &
       7.493782E-01, 7.211660E-01, 6.858999E-01, 6.506349E-01, &
       6.158184E-01, 5.810415E-01, 5.463042E-01, 4.945902E-01, &
       4.437402E-01, 3.928911E-01, 3.433811E-01, 2.944031E-01, &
       2.467411E-01, 2.003501E-01, 1.562241E-01, 1.136021E-01, &
       6.372006E-02, 2.801004E-02, 6.960025E-03, 8.175413E-09, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00/)
  INTEGER(KIND=i4), PARAMETER :: ngeos5 = 48

  ! ---------------------------------------
  ! Data obtained from the climatology file
  ! ---------------------------------------
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: latvals, lonvals, timevals
  REAL(KIND=r4), DIMENSION(:,:,:,:), ALLOCATABLE :: Gas_profiles, wgh_ozo_pro
  INTEGER(KIND=i4), DIMENSION(:,:,:,:), ALLOCATABLE :: idx_ozo_pro

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

  ! -----------------------
  ! To find the right swath
  ! -----------------------
  INTEGER   (KIND=i4),                    PARAMETER :: nmonths = 12
  CHARACTER (LEN=3), DIMENSION (nmonths), PARAMETER :: &
    months = (/ &
    'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

  ! ---------------------------
  ! 32bit/64bit C_LONG integers
  ! ---------------------------
  INTEGER (KIND=C_LONG), PARAMETER :: zerocl = 0, onecl = 1, twocl = 2
  ! Integer parameters
  INTEGER (KIND=i4), PARAMETER :: one = 1, two = 2

  ! -----------
  ! ISCCP stuff
  ! -----------
  ! --------------------------------------------
  ! TYPE declaration for ISCCP cloud climatology
  ! --------------------------------------------
  INTEGER (KIND=i4), PARAMETER, PRIVATE :: nlat_isccp=72, nlon_isccp=6596
  TYPE :: CloudClimatology
     REAL    (KIND=r8)                         :: &
          scale_ctp, scale_cfr, delta_lat, missval_cfr, missval_ctp
     REAL    (KIND=r4), DIMENSION (nlat_isccp) :: latvals, delta_lon
     INTEGER (KIND=i4), DIMENSION (nlat_isccp) :: n_lonvals
     REAL    (KIND=r4), DIMENSION (nlon_isccp) :: lonvals
     REAL    (KIND=r8), DIMENSION (nlon_isccp) :: cfr, ctp
  END TYPE CloudClimatology
  ! ----------------------------------------------
  ! Composite variable for ISCCP Cloud Climatology
  ! ----------------------------------------------
  TYPE (CloudClimatology) :: ISCCP_CloudClim

  ! --------------------------
  !(3) ISCCP Cloud Climatology
  ! --------------------------
  CHARACTER (LEN=10), PARAMETER :: isccp_lat_field  = 'ISCCP_Lats'
  CHARACTER (LEN=15), PARAMETER :: isccp_dlon_field = 'ISCCP_DeltaLons'
  CHARACTER (LEN=13), PARAMETER :: isccp_nlon_field = 'ISCCP_NumLons'
  CHARACTER (LEN=10), PARAMETER :: isccp_lon_field  = 'ISCCP_Lons'
  CHARACTER (LEN=30), PARAMETER :: isccp_mcfr_field = 'ISCCP_MonthlyAVG_CloudFraction'
  CHARACTER (LEN=30), PARAMETER :: isccp_mctp_field = 'ISCCP_MonthlyAVG_CloudPressure'

CONTAINS

  subroutine wfamf_deallocate (errstat)
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return
    call climatology_deallocate(errstat)
    call vlidort_deallocate(errstat)
    if (errstat /= 0) return
  end subroutine wfamf_deallocate

  SUBROUTINE amf_calculation_bis (            &
      pge_idx, nt, nx, lat, lon, sza, vza, saa, vaa, time,  &
      snow, glint, xtrange, do_szoom,        &
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
    use OMSAO_indices_module, only: voc_omicld_idx, &
         voc_isccp_idx
    use OMSAO_omidata_module, only : amf_correction_type
    use output_tools, only : write_albedo, write_gas_profile, &
      write_scattering_weights, write_amf_correction
    USE OMSAO_variables_module,  ONLY: voc_amf_filenames
    use output_tools, only: read_cloud_params
    use ctrlvars, only : yn_do_he5_output, yn_stratrop
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx, pge_idx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, &
         vza, saa, vaa, terrain_height
    REAL    (KIND=r8), DIMENSION (0:nt-1),      INTENT (IN) :: time
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    LOGICAL,           DIMENSION (     0:nt-1), INTENT (IN) :: do_szoom
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
    INTEGER (KIND=i4)                                :: locerrstat, itt, spixx, epixx
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), target :: amfdiag
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), target :: amfgeo, tropospheric_amf, &
         stratospheric_amf
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), target :: l2cfr, l2ctp
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)       :: albedo
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: climatology
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: scattw
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,2) :: cli_wgh_ozo_pro
    INTEGER (KIND=i4), DIMENSION (1:nx,0:nt-1,2) :: cli_idx_ozo_pro
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), target :: surface_pressure, tropopause_pressure

    type (amf_correction_type) :: amf_corr
    logical :: yn_write_cloud_variables
    character (len=256) :: cloud_file

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
    amfdiag      = i2_missval
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

       DO itt = 0, nt-1
          spixx = xtrange(itt,1) ; epixx = xtrange(itt,2)
          saoamf(spixx:epixx,itt) = 1.0_r8
       END DO
       
    ELSE

       ! -------------------------
       ! Compute the geometric AMF
       ! -------------------------
       CALL compute_geometric_wfamf ( nt, nx, sza, vza, xtrange, amfgeo, amfdiag )

       ! -------------------------------------------------------
       ! Initialize molecular AMF with geometric AMF. Subsequent
       ! subroutines will replace any entries where the true
       ! molecular AMF can be computed.
       ! -------------------------------------------------------
       saoamf = amfgeo
      
       ! ----------------------------------------------------
       ! Read OMLER albedo database stored in variable albedo
       ! ----------------------------------------------------
       CALL omi_omler_albedo ( lat, lon, albedo, nt, nx, xtrange, errstat)

       ! ---------------------------------------
       ! Write the albedo to the output file he5
       ! ---------------------------------------
       IF (do_write) then
          if (yn_do_he5_output) then
             CALL write_albedo_he5 ( albedo, nt, nx, locerrstat)
          endif
          call write_albedo (albedo, nx, nt, errstat)
          if (errstat /= 0) return
       endif

       ! -----------------------------
       ! Read the OMI L2 cloud product
       ! -----------------------------
       cloud_file = voc_amf_filenames(voc_omicld_idx)
       if (0 /= index (cloud_file, ".he5", .true.)) then
          ! FIXME: amf_read_omiclouds to be removed
          CALL amf_read_omiclouds ( nt, nx, do_szoom, l2cfr, l2ctp, errstat )
       else if (0 /= index (cloud_file, ".nc", .true.)) then
          call read_cloud_params (cloud_file, nt, nx, l2cfr, l2ctp, errstat)
       else
          call tell_error (tell_runtime_error, "unexpected cloud file extension: "//trim(cloud_file), errstat)
          return
       endif
       if (errstat /= 0) then
          call tell_error (tell_io_read_error, "reading cloud file: "//trim(cloud_file), errstat)
          return
       endif
       call tell_log (1, 'Read cloud-top pressure, cloud fraction from: '//trim(cloud_file))

       ! ----------------------------
       ! Read ISCCP cloud climatology
       ! ----------------------------
       cloud_file = voc_amf_filenames(voc_isccp_idx)
       CALL voc_amf_readisccp  ( errstat )
       if (errstat /= 0) then
          call tell_error (tell_io_read_error, "reading ISCCP cloud file: "//trim(cloud_file), errstat)
          return
       endif
       call tell_log (1, 'Read ISCCP climatology from: '//trim(cloud_file))

       ! ------------------------------------------------
       ! Read climatology and interpolate to lon/lat/time
       ! ------------------------------------------------
       CALL omi_climatology (pge_idx, climatology, cli_wgh_ozo_pro, &
            cli_idx_ozo_pro, lat, lon, time, nt, nx, xtrange, errstat, amfdiag)
       ! -------------------------------------
       ! Write the climatology to the he5 file
       ! -------------------------------------
       IF (do_write) then
          if (yn_do_he5_output) then
             CALL write_climatology_he5 (climatology, nt, nx, CmETA, &
                  locerrstat)
          endif
          call write_gas_profile (climatology, nx, nt, CmETA, errstat)
          if (errstat /= 0) return
       endif

       ! ------------------------------------------------------------------
       ! Read VLIDORT look up table. Variables are declared at module level
       ! (Input is read only on the first pass.  Subsequent passes use
       ! cached values)
       ! ------------------------------------------------------------------
       CALL read_vlidort (errstat)
       if (errstat /= 0) then
          call vlidort_deallocate(errstat)
          return
       endif

       ! ----------------------------------------------------------------------
       ! amfdiag is used to keep track of the pixels were enough information is
       ! available to carry on the AMFs calculation.
       ! ----------------------------------------------------------------------
       CALL amf_diagnostic (nt, nx, lat, lon, &
            sza, vza, snow, glint, xtrange, &
            l2cfr, l2ctp, &
            amfdiag  )

       WHERE ((saocol <= r8_missval).or.(saodco<=r8_missval))
          amfdiag = i2_missval
       END WHERE

       ! --------------------------------------------------------
       ! Compute Scattering weights in the look up table grid but
       ! with the correct albedo. amfdiag is used to skip pixel
       ! ---------------------------------------------------------
       CALL compute_scatt ( nt, nx, albedo, sza, vza, saa, vaa, l2ctp, l2cfr, &
            terrain_height, surface_pressure, cli_wgh_ozo_pro, cli_idx_ozo_pro, &
            lat, lon, amfdiag, scattw)

       ! -----------------------------------------------------------------
       ! Work out the AMF using the scattering weights and the climatology
       ! Work out Averaging Kernels
       ! -----------------------------------------------------------------
       CALL compute_amf ( nt, nx, CmETA, climatology, &
            scattw, saoamf, stratospheric_amf, tropospheric_amf, surface_pressure, tropopause_pressure, lat, lon, amfdiag, &
            locerrstat)

       ! -----------------------------------------------------------------
       ! Write out scattering weights, altitude grid and averaging kernels
       ! -----------------------------------------------------------------
       IF (do_write) then
          if (yn_do_he5_output) then
             CALL write_scatt_he5 (scattw, nt, nx, CmETA, locerrstat)
          endif
          call write_scattering_weights (scattw, nx, nt, CmETA, errstat)
          if (errstat /= 0) return
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
      if (yn_do_he5_output) then
        CALL he5_amf_write ( nx, nt, saocol, saodco, saoamf, &
                            amfgeo, amfdiag, l2cfr, l2ctp, surface_pressure, &
                            tropopause_pressure, tropospheric_amf, &
                            stratospheric_amf, locerrstat )
      endif
      amf_corr % amf_molecule_specific => saoamf
      amf_corr % amf_molecule_stratospheric => stratospheric_amf
      amf_corr % amf_molecule_tropospheric => tropospheric_amf
      amf_corr % amf_geometric => amfgeo
      amf_corr % diagnostic_flag => amfdiag
      amf_corr % cloud_fraction => l2cfr
      amf_corr % cloud_pressure => l2ctp
      amf_corr % surface_pressure => surface_pressure
      amf_corr % tropopause_pressure => tropopause_pressure
      yn_write_cloud_variables = .TRUE.
      call write_amf_correction (nx, nt, amf_corr, saocol, saodco, &
                                 yn_write_cloud_variables, errstat)
      if (errstat /= 0) return
    endif

  END SUBROUTINE amf_calculation_bis

  SUBROUTINE omi_climatology (pge_idx, climatology, cli_wgh_ozo_pro, &
    cli_idx_ozo_pro, lat, lon, time, nt, nx, xtrange, errstat, amfdiag)
    
    ! =========================================
    ! Extract Gas climatology to granule pixels
    ! No interpolation or something like that,
    ! Just pick the closest model grid
    ! =========================================
    use OMSAO_omidata_module, only: omi_scattfail_amf
    USE OMSAO_indices_module, ONLY: sao_molecule_names, pge_h2o_idx
    USE omi_pge_fitting_aux, ONLY: convert_tai_to_utc
    USE OMSAO_parameters_module, ONLY: nUTCdim
    USE OMSAO_linterpolation_module, ONLY: lininterpol, GetNode
    USE OMSAO_errstat_module, only : he5_stat_fail
    USE OMSAO_he5_module, ONLY: he5_swclose, he5_swopen, &
         he5f_acc_rdonly, granule_month, he5_swattach, &
         he5_swinqswath, he5_swinqdflds

    IMPLICIT NONE
    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx, pge_idx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon
    REAL    (KIND=r8), DIMENSION (0:nt-1), INTENT (IN) :: time
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
    
    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1, CmETA), INTENT (INOUT) :: climatology
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1, 2), INTENT (INOUT) :: cli_wgh_ozo_pro
    INTEGER (KIND=i4), DIMENSION(1:nx,0:nt-1, 2), INTENT (INOUT) :: cli_idx_ozo_pro
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: amfdiag

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: itimes, ixtrack, spix, epix, n, status , &
         nlon, nlat, ntim, locerrstat
    INTEGER (KIND=i4), DIMENSION(2) :: idx_lat, idx_lon, idx_tim
    REAL (KIND=r8) :: llon, llat, ltime

    INTEGER (KIND=i4) :: swath_file_id, swath_id, nswath, ndatafields
    INTEGER   (KIND=i4), DIMENSION(10) :: datafield_rank, datafield_type
    CHARACTER (LEN=MAX_STR_LEN) :: swath_file, swath_name, locswathname, datafield_name, &
         gasdatafieldname, h2odatafieldname
    INTEGER (KIND=C_LONG) :: nswathcl, swlen
    INTEGER (KIND=i2), DIMENSION(nUTCdim) :: time_utc
    character (len=72) :: logmsg    

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    if (errstat /= 0) return

    locerrstat = 0
    
    ! ----------------------------------------------------
    ! Open climatology file (done here to do it only once)
    ! and look for the swath and data fields we are need.
    ! ----------------------------------------------------
    if (granule_month < 1 .OR. granule_month > nmonths) then
       call tell_error (tell_runtime_error, "omi_climatology: invalid month", &
            errstat)
       return
    endif
    swath_file = TRIM(ADJUSTL(OMSAO_climatology_filename))
  
    ! --------------------------------------------------------------
    ! Open he5 OMI climatology and check SWATH_FILE_ID (-1 if error)
    ! --------------------------------------------------------------
    swath_file_id = HE5_SWOPEN (swath_file, he5f_acc_rdonly)
    IF (swath_file_id == he5_stat_fail) THEN
       call tell_error (tell_io_open_error, "omi_climatology: opening climatology file"//trim(swath_file), &
            errstat)
       RETURN
    END IF
  
    ! -----------------------------------------------------------
    ! Check for existing HE5 swathw and attach to the one we need
    ! -----------------------------------------------------------
    swath_name = "" !JED
    nswathcl = HE5_SWinqswath(TRIM(ADJUSTL(swath_file)), swath_name, swlen )
    nswath   = INT(nswathcl, KIND=i4 )
    
    ! ----------------------------------------------------------------
    ! If there is only one swath in the file, we can attach to it but
    ! if there are more (NSWATH > 1), then we must find the swath that
    ! corresponds to the current month.
    ! ----------------------------------------------------------------
    IF (nswath > 1) THEN
       CALL extract_swathname(nswath, TRIM(ADJUSTL(swath_name)), &
            TRIM(ADJUSTL(months(granule_month))), locswathname)
       ! ---------------------------------------------------------------------------
       ! Check if we found the correct swath name. If not, report an error and exit.
       ! ---------------------------------------------------------------------------
       IF ( INDEX (TRIM(ADJUSTL(locswathname)),TRIM(ADJUSTL(months(granule_month)))) == 0 ) THEN
          call tell_error (tell_runtime_error, "omi_climatology: finding month in "// &
               trim(locswathname), errstat)
          RETURN
       END IF
    ELSE
       locswathname = TRIM(ADJUSTL(swath_name))
    END IF
    
    ! -----------------------------
    ! Attach to current month swath
    ! -----------------------------
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(locswathname)) )
    IF ( swath_id == he5_stat_fail ) THEN
       call tell_error (tell_io_error, "omi_climatology: attaching to swath "// &
            trim(locswathname), errstat)
       RETURN
    END IF
    
    ! -----------------------------------------------------------------------
    ! Finding out the data field for the gas of interest (.eq. to target gas)
    ! -----------------------------------------------------------------------
    datafield_name=""
    ndatafields = HE5_swinqdflds(swath_id, datafield_name, datafield_rank, datafield_type)
    CALL extract_swathname(ndatafields, TRIM(ADJUSTL(datafield_name)), &
      TRIM(ADJUSTL(sao_molecule_names(pge_idx))), gasdatafieldname)

    ! ---------------------------------------------------------------------------
    ! Check if we found the correct swath name. If not, report an error and exit.
    ! ---------------------------------------------------------------------------
    IF ( INDEX (TRIM(ADJUSTL(gasdatafieldname)),TRIM(ADJUSTL(sao_molecule_names(pge_idx)))) == 0 ) THEN
      call tell_error (tell_runtime_error, "omi_climatology: climatology file "// &
                       "does not contain data for "// trim(sao_molecule_names(pge_idx)), &
                       errstat)
      RETURN
    END IF

    ! -------------------------------------------------------------------
    ! Finding out the data field for water vapor (to compute air density)
    ! -------------------------------------------------------------------
    h2odatafieldname=""
    CALL extract_swathname(ndatafields, TRIM(ADJUSTL(datafield_name)), &
         TRIM(ADJUSTL(sao_molecule_names(pge_h2o_idx))), h2odatafieldname)
    
    ! ---------------------------------------------------------------------------
    ! Check if we found the correct swath name. If not, report an error and exit.
    ! ---------------------------------------------------------------------------
    IF ( INDEX (TRIM(ADJUSTL(h2odatafieldname)),TRIM(ADJUSTL(sao_molecule_names(pge_h2o_idx)))) == 0 ) THEN
      call tell_error (tell_runtime_error, "omi_climatology: climatology file "// &
                       "does not contain data for "// trim(sao_molecule_names(pge_h2o_idx)), &
                       errstat)
      RETURN
    END IF

    ! ------------------------------------------------------------
    ! Find the Climatology corresponding to each lat and lon pixel
    ! ------------------------------------------------------------
    DO itimes = 0, nt-1

       ! Only work out climatology if we have complete geolocation information
       IF (time(itimes) /= r8_missval) THEN ! FIXME
          CALL convert_tai_to_utc(nUTCdim, time(itimes), time_utc(1:nUTCdim))
          ltime = REAL(time_utc(4),KIND=r8)+REAL(time_utc(5),KIND=r8)/60.0
          IF (ltime .LT. MINVAL(timevals)) ltime = REAL(MINVAL(timevals),KIND=r8)
          IF (ltime .GT. MAXVAL(timevals)) ltime = REAL(MAXVAL(timevals),KIND=r8)
       ELSE
         write(logmsg, '(a,1x,i4)')'Failed to read time for step', itimes
         call tell_log (1, logmsg)
         amfdiag(:,itimes) = omi_scattfail_amf
         CYCLE
       END IF

       spix = xtrange(itimes,1); epix = xtrange(itimes,2)
       
       DO ixtrack = spix, epix
          
          llon = REAL(lon(ixtrack,itimes),KIND=r8)
          llat = REAL(lat(ixtrack,itimes),KIND=r8)

          ! Only complete climatology calculation if we have complete geolocation information
          IF ( (lon(ixtrack,itimes) /= r4_missval) .AND. (lat(ixtrack,itimes) /= r4_missval) ) THEN
             IF (llon .LT. MINVAL(lonvals)) llon = REAL(MINVAL(lonvals),KIND=r8)
             IF (llon .GT. MAXVAL(lonvals)) llon = REAL(MAXVAL(lonvals),KIND=r8)
             IF (llat .LT. MINVAL(latvals)) llat = REAL(MINVAL(latvals),KIND=r8)
             IF (llat .GT. MAXVAL(latvals)) llat = REAL(MAXVAL(latvals),KIND=r8)
          ELSE
             CYCLE
          END IF

          ! ------------------------------------------------------------------------
          ! Given the values of lonvals, latvals, timevals and llon, llat, and ltime
          ! determine the indices of climatolgoy values to be read.
          ! Using linear interpolation only 2 nodes needed in each dimension or if 
          ! outbounds, closest node is selected
          ! ------------------------------------------------------------------------
          CALL GetNode(REAL(lonvals,KIND=r8),llon, &
               idx_lon(1), 'Lower')
          IF (idx_lon(1) .EQ. -2) idx_lon(1) = 1
          IF (idx_lon(1) .EQ. -3) idx_lon(1) = Cmlon-1
          CALL GetNode(REAL(lonvals,KIND=r8),llon, &
               idx_lon(2), 'Upper')
          IF (idx_lon(2) .EQ. -2) idx_lon(2) = 2
          IF (idx_lon(2) .EQ. -3) idx_lon(2) = Cmlon
          nlon = idx_lon(2)-idx_lon(1) + 1
          IF (nlon == 1) THEN
             IF (idx_lon(1) <= Cmlon-1) THEN
                idx_lon(2) = idx_lon(1) + 1
             ELSE IF (idx_lon(1) == Cmlon) THEN
                idx_lon(1) = Cmlon-1
                idx_lon(2) = Cmlon
             ENDIF
             nlon = idx_lon(2)-idx_lon(1)+1
          ENDIF

          CALL GetNode(REAL(latvals,KIND=r8),llat, &
               idx_lat(1), 'Lower')
          IF (idx_lat(1) .EQ. -2) idx_lat(1) = 1
          IF (idx_lat(1) .EQ. -3) idx_lat(1) = Cmlat-1
          CALL GetNode(REAL(latvals,KIND=r8),llat, &
               idx_lat(2), 'Upper')
          IF (idx_lat(2) .EQ. -2) idx_lat(2) = 2
          IF (idx_lat(2) .EQ. -3) idx_lat(2) = Cmlat
          nlat = idx_lat(2)-idx_lat(1) + 1
          IF (nlat == 1) THEN
             IF (idx_lat(1) <= Cmlat-1) THEN
                idx_lat(2) = idx_lat(1) + 1
             ELSE IF (idx_lat(1) == Cmlat) THEN
                idx_lat(1) = Cmlat-1
                idx_lat(2) = Cmlat
             ENDIF
             nlat = idx_lat(2)-idx_lat(1)+1
          ENDIF

          CALL GetNode(REAL(timevals,KIND=r8),ltime, &
               idx_tim(1), 'Lower')
          IF (idx_tim(1) .EQ. -2) idx_tim(1) = 1
          IF (idx_tim(1) .EQ. -3) idx_tim(1) = CmHRS-1
          CALL GetNode(REAL(timevals,KIND=r8),ltime, &
               idx_tim(2), 'Upper')
          IF (idx_tim(2) .EQ. -2) idx_tim(2) = 2
          IF (idx_tim(2) .EQ. -3) idx_tim(2) = CmHRS
          ntim = idx_tim(2)-idx_tim(1) + 1
          IF (ntim == 1) THEN
             IF (idx_tim(1) <= CmHRS-1) THEN
                idx_tim(2) = idx_tim(1) + 1
             ELSE IF (idx_tim(1) == CmHRS) THEN
                idx_tim(1) = CmHRS-1
                idx_tim(2) = CmHRS
             ENDIF
             ntim = idx_tim(2)-idx_tim(1)+1
          ENDIF

          ! -----------------------------
          ! Allocate Climatology arrays
          ! nlon, nlat, CmETA, ntim since
          ! linear interpolation is used
          ! -----------------------------
          CALL climatology_allocate (nlat,nlon,CmETA,ntim,locerrstat)
          IF ( locerrstat /= 0 ) THEN
             call tell_error (tell_runtime_error, "omi_climatology: Climatology allocate failed", errstat)
             RETURN
          END IF

          ! ----------------
          ! Read climatology
          ! ----------------
          CALL read_climatology(idx_lon, idx_lat, idx_tim, &
               swath_id, gasdatafieldname, locerrstat)
          IF ( locerrstat /= 0 ) THEN
             call tell_error (tell_runtime_error, "omi_climatology: read climatology failed", errstat)
             RETURN
          END IF

          ! ----------------------------------------------------------
          ! Save ozone profile index and weight.
          ! Choose climatology pixel closest to satellite pixel center
          ! ----------------------------------------------------------
          cli_wgh_ozo_pro(ixtrack,itimes,1:2) = REAL(wgh_ozo_pro(MINLOC(ABS(lonvals(idx_lon)-llon),1), &
               MINLOC(ABS(latvals(idx_lat)-llat),1), &
               MINLOC(ABS(timevals(idx_tim)-ltime),1),1:2), KIND=r8)
          cli_idx_ozo_pro(ixtrack,itimes,1:2) = idx_ozo_pro(MINLOC(ABS(lonvals(idx_lon)-llon),1), &
               MINLOC(ABS(latvals(idx_lat)-llat),1), &
               MINLOC(ABS(timevals(idx_tim)-ltime),1),1:2)

          DO n = 1, CmETA
             ! Interpolate trace gas profile to lon,lat,hrs [GAS]/cm^2
             climatology(ixtrack,itimes,CmETA-n+1) = linInterpol(nlon,nlat,ntim, &
                  REAL(lonvals(idx_lon(1):idx_lon(2)),KIND=r8), &
                  REAL(latvals(idx_lat(1):idx_lat(2)),KIND=r8), &
                  REAL(timevals(idx_tim(1):idx_tim(2)),KIND=r8), &
                  REAL(Gas_profiles(1:nlon,1:nlat,n,1:ntim),KIND=r8), &
                  llon, llat, ltime, status=status)
             IF ( locerrstat /= 0 ) THEN
                call tell_error (tell_runtime_error, &
                     "omi_climatology: gas interpolation failed", errstat)
                RETURN
             END IF
          END DO

          !  Set non-physical entries to zero.
          WHERE ( climatology(ixtrack,itimes,1:CmETA) < 0.0_r8 )
             climatology(ixtrack,itimes,1:CmETA) = 0.0_r8
          END WHERE

          ! De-allocate climatology arrays
          CALL climatology_deallocate_arrays(locerrstat)
          IF ( locerrstat /= 0 ) THEN
             call tell_error (tell_runtime_error, &
                  "omi_climatology: climatology_deallocate_arrays failed", errstat)
             RETURN
          END IF
          
       END DO
       write(logmsg, '(a,1x,i5)')'Preparing climatology line', itimes
       call tell_log (1, logmsg)
    END DO
   
    ! Close climatology file (done here to do it only once)
    errstat = HE5_SWCLOSE(swath_file_id)
    IF ( errstat == he5_stat_fail ) THEN
       call tell_error (tell_io_error, "read_climatology_dimensions: closing climatology file "// &
            trim(swath_file), errstat)
       RETURN
    END IF

  END SUBROUTINE omi_climatology
  
  SUBROUTINE read_climatology (idx_lon, idx_lat, idx_tim, swath_id, &
       gasdatafieldname, errstat)
    ! ==========================================================
    ! This subroutine reads in the climatology from GEOS-Chem or
    ! other source. The climatology file needs to be conform to
    ! the format assumed here.
    ! ==========================================================    
    USE OMSAO_he5_module, ONLY: HE5_SWrdfld, HE5_SWrdlattr
    IMPLICIT NONE

    INTEGER (KIND=i4), INTENT (IN) :: swath_id
    INTEGER (KIND=i4), DIMENSION(2), INTENT (IN) :: idx_lon, idx_lat, idx_tim
    CHARACTER (LEN=MAX_STR_LEN) :: gasdatafieldname

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4) :: he5stat
    INTEGER   (KIND=C_LONG) :: nlon, nlat, nlev, ntim
    CHARACTER (LEN=16), PARAMETER :: wgh_ozo_pro_field = 'LUT Ozone Weight'
    CHARACTER (LEN=17), PARAMETER :: idx_ozo_pro_field = 'LUT Ozone Profile'
    REAL (KIND=r4) :: scale_gas, scale_wgh, scale_idx

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    if (errstat /= 0) return

    ! -----------------------------------------------------------------------
    ! Create KIND=4/KIND=8 variables. We have to use the vertical dimension a
    ! few times in this subroutine, so it saves some typing if we do the 
    ! conversion once and save them in new variables.
    ! -----------------------------------------------------------------------
    ! ---------------------------------------------------------
    ! Work out the number of nodes to be read in each dimension
    ! ---------------------------------------------------------
    nlon = INT(idx_lon(2)-idx_lon(1)+1,KIND=C_LONG)
    nlat = INT(idx_lat(2)-idx_lat(1)+1,KIND=C_LONG)
    nlev = INT(CmETA, KIND=C_LONG )
    ntim = INT(idx_tim(2)-idx_tim(1)+1,KIND=C_LONG)
    
    if (errstat /= 0) return

    ! ----------------------------
    ! Read data from gas datafield
    ! ----------------------------
    he5_start_4d  = (/ INT(idx_lon(1)-1,8),INT(idx_lat(1)-1,8),zerocl,INT(idx_tim(1)-1,8) /)
    he5_stride_4d = (/ onecl, onecl, onecl, onecl /)
    he5_edge_4d   = (/  nlon,  nlat, nlev, ntim /)

    he5stat = HE5_SWrdfld ( &
      swath_id, TRIM(ADJUSTL(gasdatafieldname)), &
      he5_start_4d, he5_stride_4d, he5_edge_4d, &
      Gas_profiles(1:nlon,1:nlat,1:nlev,1:ntim) )
    
    ! -----------------------------------------
    ! Read gas datafield scale factor attribute
    ! -----------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(gasdatafieldname)),&
      "ScaleFactor", scale_gas       )

    ! ---------------------------------------
    ! Read data from ozone profile data field
    ! ---------------------------------------
    he5_start_4d  = (/ INT(idx_lon(1)-1,8),INT(idx_lat(1)-1,8),INT(idx_tim(1)-1,8),zerocl /)
    he5_stride_4d = (/ onecl, onecl, onecl, onecl /)
    he5_edge_4d   = (/  nlon,  nlat, ntim, INT(2,8) /)
    he5stat = HE5_SWrdfld ( &
      swath_id, TRIM(ADJUSTL(idx_ozo_pro_field)), &
      he5_start_4d, he5_stride_4d, he5_edge_4d, &
      idx_ozo_pro(1:nlon,1:nlat,1:ntim,1:2) )
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(idx_ozo_pro_field)),&
      "ScaleFactor", scale_idx       )
    ! ----------------------------------------------
    ! Read data from ozone profile weight data field
    ! ----------------------------------------------
    he5_start_4d  = (/ INT(idx_lon(1)-1,8),INT(idx_lat(1)-1,8),INT(idx_tim(1)-1,8),zerocl /)
    he5_stride_4d = (/ onecl, onecl, onecl, onecl /)
    he5_edge_4d   = (/  nlon,  nlat, ntim, INT(2,8) /)
    he5stat = HE5_SWrdfld ( &
      swath_id, TRIM(ADJUSTL(wgh_ozo_pro_field)), &
      he5_start_4d, he5_stride_4d, he5_edge_4d, &
      wgh_ozo_pro(1:nlon,1:nlat,1:ntim,1:2) )
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(wgh_ozo_pro_field)),&
      "ScaleFactor", scale_wgh       )
    
    ! ------------------------------------
    ! Apply scaling factors to data fields
    ! ------------------------------------
    Gas_profiles = Gas_profiles * scale_gas
    idx_ozo_pro  = idx_ozo_pro  * INT(scale_idx,4)
    wgh_ozo_pro  = wgh_ozo_pro  * scale_wgh

  END SUBROUTINE read_climatology

  SUBROUTINE omi_omler_albedo( lat, lon, albedo, nt, nx, xtrange, &
      errstat)

    ! ==================================================================
    ! This subroutine reads the OMLER albedo data base for the month of
    ! the orbit to processed. Then it interpolates the values for each
    ! one of the pixels of the orbit to be analyzed
    ! ==================================================================
    USE OMSAO_linterpolation_module, ONLY: lininterpol, GetNode
    USE OMSAO_variables_module, ONLY: OMSAO_OMLER_filename, &
      winwav_min, winwav_max
    USE ezspline_interpolation, ONLY: ezspline_1d_interpolation, &
      ezspline_2d_interpolation
    USE OMSAO_errstat_module, only : he5_stat_fail, pge_errstat_ok
    USE OMSAO_he5_module, ONLY: HE5_GDOPEN, HE5_GDattach, HE5_GDRDFLD, &
         HE5_GDRDLATTR, HE5_GDDETACH, HE5_GDclose, he5f_acc_rdonly, &
         granule_month

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (INOUT) :: albedo

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
    INTEGER (KIND=i4) :: itimes, ixtrack, spix, epix, ilon, ilat, nlon, &
      nlat, OMnwvl, grid_id, grid_file_id, month, minwvl, maxwvl
    INTEGER (KIND=i4), DIMENSION(2) :: lon_idx, lat_idx
    REAL (KIND=r4) :: scale_factor, offset
    REAL (KIND=r8) :: lonp, latp
    REAL (KIND=r8), DIMENSION(1) :: midwvl

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
    midwvl(1) = REAL((amf_wvl + amf_wvl2) / 2.0_r8,KIND=r8)
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

    OMLER_wvl_albedo = REAL(offset, KIND = r8) +          &
      REAL(scale_factor, KIND = r8)*      &
      REAL(OMLER_monthly_albedo, KIND=r8)

    ! ----------------------------------------------------------------
    ! Interpolate for each pixel to one single wavelenght, just at the
    ! mid point of the fitting window: midwvl.
    ! ----------------------------------------------------------------
    DO ilon = 1, OMLER_n_longitudes
      DO ilat = 1, OMLER_n_latitudes
        CALL ezspline_1d_interpolation ( &
          OMnwvl, REAL(OMLER_wvl(minwvl:maxwvl), KIND=r8), &
          OMLER_wvl_albedo(ilon,ilat,1:OMnwvl,1), &
          one, midwvl, OMLER_albedo(ilon,ilat), locerrstat )
      END DO
    END DO

    ! --------------------------------------------------
    ! Interpolate to the lat and longitude of each pixel
    ! --------------------------------------------------
    DO itimes = 0, nt-1

       spix = xtrange(itimes,1); epix = xtrange(itimes,2)
       DO ixtrack = spix, epix
  
          lonp = REAL(lon(ixtrack,itimes), KIND=r8)
          latp = REAL(lat(ixtrack,itimes), KIND=r8)

          ! Only work out surface reflectance if we have geolocation information
          IF ( (lon(ixtrack,itimes) /= r4_missval) .AND. (lat(ixtrack,itimes) /= r4_missval) ) THEN
             ! Be sure that lonp and latp are within surface albedo boundaries
             IF (lonp .LT. MINVAL(OMLER_longitude)) lonp = REAL(MINVAL(OMLER_longitude),KIND=r8)
             IF (lonp .GT. MAXVAL(OMLER_longitude)) lonp = REAL(MAXVAL(OMLER_longitude),KIND=r8)
             IF (latp .LT. MINVAL(OMLER_latitude))  latp = REAL(MINVAL(OMLER_latitude),KIND=r8)
             IF (latp .GT. MAXVAL(OMLER_latitude))  latp = REAL(MAXVAL(OMLER_latitude),KIND=r8)
          ELSE
             CYCLE
          END IF

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
               REAL(OMLER_longitude(lon_idx(1):lon_idx(2)),KIND=r8), &
               REAL(OMLER_latitude(lat_idx(1):lat_idx(2)),KIND=r8), &
               OMLER_albedo(lon_idx(1):lon_idx(2),lat_idx(1):lat_idx(2)), &
               lonp, latp, status=locerrstat)
          if (locerrstat /= 0) then
             call tell_error (tell_runtime_error, &
                  "omi_omler_albedo: lon/lat interpolation failed", errstat)
             return
          endif
  
       END DO
    END DO

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
   
    errstat = MAX(errstat, locerrstat)

  END SUBROUTINE omi_omler_albedo

  SUBROUTINE extract_swathname ( nswath, multi_swath, swathstr, single_swath )

    ! ---------------------------------------------------------------------
    ! Extracts SINGLE_SWATH from MULTI_SWATH, based on presence of SWATHSTR
    ! ---------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4), INTENT (IN) :: nswath
    CHARACTER (LEN=*),   INTENT (IN) :: multi_swath, swathstr

    ! ----------------
    ! Output variables
    ! ----------------
    CHARACTER (LEN=*), INTENT (OUT) :: single_swath

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4), DIMENSION(0:nswath) :: swsep
    INTEGER   (KIND=i4)                      :: mslen, k, j1, j2, nsep
    CHARACTER (LEN=LEN(multi_swath))         :: tmpstr

    call tell_log (1, 'extract_swathname:  multi_swath='//trim(multi_swath))
    call tell_log (1, 'extract_swathname:  swathstr='//swathstr)

    ! --------------------------
    ! Initialize output variable
    ! --------------------------
    single_swath = ''

    ! ---------------------------------
    ! Find length of MULTI_SWATH string
    ! ---------------------------------
    tmpstr = TRIM(ADJUSTL(multi_swath))
    mslen  = LEN_TRIM(ADJUSTL(tmpstr))

    ! --------------------------------------
    ! First find the number of "," in TMPSTR
    ! --------------------------------------
    nsep = 0 ; swsep(0:nswath) = 0
    DO k = 1, mslen
      IF ( tmpstr(k:k) == ',' ) THEN
        nsep = nsep + 1
        swsep(nsep) = k
      END IF
    END DO
    IF ( nsep == nswath-1 ) THEN
      nsep = nsep + 1 ; swsep(nsep) = mslen + 1
    END IF

    ! ----------------------------------------------------------------------
    ! Hangle along the positions of separators (commas, ",") between the
    ! concatinated Swath Name entries. The first Swath Name to contain
    ! SWATHSTR is taken as the match - not a perfect rationale but simple
    ! enough if we have set up the AMF table file correctly.
    ! ----------------------------------------------------------------------
    getswath: DO k = 1, nswath
      j1 = swsep(k-1)+1  ;  j2 = swsep(k)-1
      IF ( INDEX ( tmpstr(j1:j2), TRIM(ADJUSTL(swathstr)) ) > 0 ) THEN
        single_swath = TRIM(ADJUSTL(tmpstr(j1:j2)))
        EXIT getswath
      END IF
    END DO getswath

    RETURN
  END SUBROUTINE extract_swathname

  SUBROUTINE climatology_getdim ( &
      swath_id, errstat )

    ! --------------------------------
    ! Return dimensions of Climatology
    ! --------------------------------
    USE OMSAO_he5_module, ONLY: HE5_SWinqdims, HE5_SWrdfld, HE5_SWrdlattr

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: swath_id

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4), PARAMETER :: maxdim = 100
    INTEGER   (KIND=i4) :: ndim, nsep, he5stat, i, j, swlen, iend, istart
    INTEGER   (KIND=i4), DIMENSION(0:maxdim) :: dim_array, dim_seps
    INTEGER   (KIND=C_LONG) :: ndimcl
    INTEGER   (KIND=C_LONG), DIMENSION(0:maxdim) :: dim_arraycl

    CHARACTER (LEN=8), PARAMETER :: cli_lat_field = 'Latitude'
    CHARACTER (LEN=9), PARAMETER :: cli_lon_field = 'Longitude'
    CHARACTER (LEN=8), PARAMETER :: cli_time_field = 'UTC_time'
    CHARACTER (LEN=10*maxdim) :: dim_chars

    REAL (KIND=r4) :: scale_lat, scale_lon, scale_time

    if (errstat /= 0) return
    ! ---------------------------
    ! Initialize output variables
    ! ---------------------------
    Cmlat = -1 ; Cmlon = -1 ; CmETA = -1; CmHRS = -1

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    dim_chars="" !JED
    ndimcl = HE5_SWinqdims  ( swath_id, dim_chars, dim_arraycl(0:maxdim) )
    ndim   = INT ( ndimcl, KIND=i4 )
    IF ( ndim <= 0 ) THEN
      call tell_error (tell_runtime_error, "in climatology_getdim", errstat)
      RETURN
    END IF
    dim_array(0:maxdim) = INT ( dim_arraycl(0:maxdim), KIND=i4 )

    dim_chars =     TRIM(ADJUSTL(dim_chars))
    swlen     = LEN_TRIM(ADJUSTL(dim_chars))

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    nsep = 0 ; dim_seps(0:ndim) = 0
    getseps: DO i = 1, swlen
      IF ( dim_chars(i:i) == ',' ) THEN
        nsep           = nsep + 1
        dim_seps(nsep) = i
      END IF
    END DO getseps
    nsep = nsep + 1 ; dim_seps(nsep) = swlen+1

    ! --------------------------------------------------------------------
    ! Hangle along the NSEP indices until we have found the two dimensions
    ! we are interested in.
    ! --------------------------------------------------------------------
    getdims:DO j = 0, nsep-1
      istart = dim_seps(j)+1 ; iend = dim_seps(j+1)-1

      SELECT CASE ( dim_chars(istart:iend) )
      CASE ( "nLat" )
        Cmlat = dim_array(j)
      CASE ( "nLon" )
        Cmlon = dim_array(j)
      CASE ( "nLev" )
        CmETA = dim_array(j)
      CASE ( "nHrs")
        CmHRS = dim_array(j)
      CASE DEFAULT
        ! Whatever. Nothing to be done here.
      END SELECT

    END DO getdims

    ALLOCATE (latvals (Cmlat),  &
              lonvals (Cmlon),  &
              timevals (CmHRS), STAT=errstat);
    if (errstat /= 0) then
       call tell_error (tell_malloc_error, "in climatology_getdim: allocating climalogy grid arrays failed", errstat)
       return
    endif

    ! -------------------------------
    ! Read dimension-defining arrays
    ! -------------------------------
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = INT(Cmlat,KIND=C_LONG)
    he5stat = HE5_SWrdfld ( swath_id, cli_lat_field, &
      he5_start_1d, he5_stride_1d, he5_edge_1d, latvals(1:Cmlat) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = INT(Cmlon,KIND=C_LONG)
    he5stat = HE5_SWrdfld ( swath_id, cli_lon_field, &
      he5_start_1d, he5_stride_1d, he5_edge_1d, lonvals(1:Cmlon) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = INT(CmHRS,KIND=C_LONG)
    he5stat = HE5_SWrdfld ( swath_id, cli_time_field, &
      he5_start_1d, he5_stride_1d, he5_edge_1d, timevals(1:CmHRS) )

    IF ( he5stat /= 0 ) then
      call tell_error (tell_io_read_error, "in climatology_getdim: reading climatology grid failed", &
                       errstat)
      return
    endif

    ! -----------------------------------------------
    ! Read dimension-defining scale factor attributes
    ! -----------------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, cli_lat_field, "ScaleFactor", scale_lat )
    he5stat = HE5_SWrdlattr ( swath_id, cli_lon_field, "ScaleFactor", scale_lon )
    he5stat = HE5_SWrdlattr ( swath_id, cli_time_field, "ScaleFactor", scale_time )
    IF ( he5stat /= 0 ) then
      call tell_error (tell_io_read_error, "in climatology_getdim: reading climatology attributes failed", &
                       errstat)
      return
    endif

    ! -----------------------------------
    ! Apply scaling factors to geo fields
    ! -----------------------------------
    lonvals = lonvals * scale_lon
    latvals = latvals * scale_lat
    timevals = timevals * scale_time

    RETURN
  END SUBROUTINE climatology_getdim

  subroutine climatology_deallocate (errstat)
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return
    deallocate (latvals, lonvals, timevals, stat=errstat)
    if (allocated(Gas_profiles)) DEALLOCATE(Gas_profiles, stat=errstat)
    if (allocated(wgh_ozo_pro)) DEALLOCATE(wgh_ozo_pro, stat=errstat)
    if (allocated(idx_ozo_pro)) DEALLOCATE(idx_ozo_pro, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, "climatology_deallocate failed", &
           errstat)
      return
    endif
  end subroutine climatology_deallocate

  subroutine climatology_deallocate_arrays (errstat)
    implicit none
    integer, intent(inout) :: errstat

    if (errstat /= 0) return
    deallocate ( Gas_profiles, wgh_ozo_pro, idx_ozo_pro, stat=errstat)
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, "climatology_deallocate failed", &
           errstat)
      return
    endif
  end subroutine climatology_deallocate_arrays

  SUBROUTINE climatology_allocate (Cmlat, Cmlon, CmETA, CmHRS, errstat)

    !USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: Cmlat, Cmlon, CmETA, CmHRS

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4) :: estat

    if (errstat /= 0) return

    ALLOCATE (Gas_profiles(Cmlon, Cmlat, CmETA, CmHRS), &
              wgh_ozo_pro(Cmlon, Cmlat, CmHRS, 2), &
              idx_ozo_pro(Cmlon, Cmlat, CmHRS, 2), STAT=estat ) ;
    if (estat /= 0) then
      call tell_error (tell_malloc_error, "climatology_allocate: failed", errstat)
      return
    endif

    RETURN
  END SUBROUTINE climatology_allocate

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

  SUBROUTINE compute_geometric_wfamf ( nt, nx, sza, vza, xtrange, amfgeo, amfdiag )

    USE OMSAO_parameters_module, ONLY: deg2rad
    USE OMSAO_omidata_module,    ONLY: omi_geo_amf, omi_oobview_amf

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                         INTENT (IN) :: nx, nt
    REAL    (KIND=r4), DIMENSION (nx,0:nt-1),  INTENT (IN) :: sza, vza
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2), INTENT (IN) :: xtrange

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: amfgeo
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: amfdiag

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: it, spix, epix

    ! ---------------------------------------------
    ! Compute geometric AMF and set diagnostic flag
    ! ---------------------------------------------
    ! ----------------------------------------------------------------
    ! Checking is done within a loop over NT to assure that we have
    ! "missing" values in all the right places. A single comprehensive
    ! WHERE statement over "1:nx,0:nt-1" would be more efficient but
    ! would also overwrite missing values with OMI_OOBVIEW_AMF.
    ! ----------------------------------------------------------------
    DO it = 0, nt-1
      spix = xtrange(it,1) ; epix = xtrange(it,2)
      WHERE ( &
          sza(spix:epix,it) /= r4_missval .AND. &
          sza(spix:epix,it) >=     0.0_r4 .AND. &
          sza(spix:epix,it) <     90.0_r4 .AND. &
          vza(spix:epix,it) /= r4_missval .AND. &
          vza(spix:epix,it) >=     0.0_r4 .AND. &
          vza(spix:epix,it) <     90.0_r4         )
        amfgeo(spix:epix,it) = &
          1.0_r8 / COS ( REAL(sza(spix:epix,it),KIND=r8)*deg2rad ) + &
          1.0_r8 / COS ( REAL(vza(spix:epix,it),KIND=r8)*deg2rad )
        amfdiag(spix:epix,it) = omi_geo_amf
      ELSEWHERE
        amfdiag(spix:epix,it) = omi_oobview_amf
      ENDWHERE
    END DO

    RETURN
  END SUBROUTINE compute_geometric_wfamf

  SUBROUTINE read_vlidort (errstat)

    ! ====================================================
    ! This subroutine reads in the VLIDORT calculations to
    ! compute the Scattering Weights.
    ! It should check if the fitting window is included in
    ! the file, if not a warning should be printed and all
    ! the AMF diagnostic set to non computed.
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

  SUBROUTINE amf_read_omiclouds ( nt, nx, do_szoom, l2cfr, l2ctp, errstat )

    USE OMSAO_variables_module,  ONLY: voc_amf_filenames
    USE OMSAO_indices_module,    ONLY: voc_omicld_idx
    USE OMSAO_omidata_module,    ONLY: gzoom_spix, gzoom_epix, gzoom_npix
    USE OMSAO_he5_module, ONLY: HE5_SWrdlattr, HE5_SWrdfld, HE5_SWinqdims, &
      he5_init_input_file
    USE OMSAO_errstat_module, only : pge_errstat_ok

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),           INTENT (IN) :: nt, nx
    LOGICAL, DIMENSION (0:nt-1), INTENT (IN) :: do_szoom
    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: l2cfr, l2ctp

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)        :: it, nt_loc, nx_loc, locerrstat
    REAL    (KIND=r4)        :: scale_cfr, offset_cfr, missval_cfr, scale_ctp, offset_ctp
    !INTEGER (KIND=i2)       :: missval_ctp
    real (kind=r4)           :: missval_ctp
    CHARACTER (LEN=5)        :: addstr
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1) :: o4ctp
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1) :: cfr, ctp
    LOGICAL                  :: do_raman_clouds

    ! ---------------------------------------
    ! For accesing the file (local variables)
    ! ---------------------------------------
    CHARACTER (LEN=MAX_STR_LEN) :: omicloud_swath_name    = 'undefined'
    INTEGER   (KIND=i4)      :: omicloud_swath_id      = -1
    INTEGER   (KIND=i4)      :: omicloud_swath_file_id = -1

    ! ----------------------------------------------
    ! Names of various HE5 fields to read from files
    ! ----------------------------------------------
    !(1) OMI Cloud Fields
    ! -------------------
    CHARACTER (LEN=13), PARAMETER :: omicld_cfrac_field      = 'CloudFraction'
    CHARACTER (LEN=13), PARAMETER :: omicld_cpres_field      = 'CloudPressure'

    ! ----------------------
    ! Name of the subroutine
    ! ----------------------
    !CHARACTER (LEN=22), PARAMETER :: modulename = 'amf_read_omiclouds'

    locerrstat = pge_errstat_ok

    ! -------------------------
    ! Attach to OMI Cloud Swath
    ! -------------------------
    CALL he5_init_input_file ( &
      voc_amf_filenames(voc_omicld_idx), omicloud_swath_name, &
      omicloud_swath_id, omicloud_swath_file_id,              &
      nt_loc, nx_loc, locerrstat )

    ! ------------------------------------------------------------------------
    ! Usually we would compare nx and nt here, but OMIL2 contains a
    ! different value for nt (=1). Thus we skip the test.
    ! ------------------------------------------------------------------------
    IF ( locerrstat /= pge_errstat_ok ) THEN
      call tell_error (tell_io_error, "amf_read_omiclouds:  reading cloud file"// &
                       trim(voc_amf_filenames(voc_omicld_idx)), errstat)
      RETURN
    END IF

    IF ( INDEX( voc_amf_filenames(voc_omicld_idx), 'CLDRR' ) /= 0 ) THEN
      do_raman_clouds = .TRUE.
      addstr          = "forO3"
    ELSE
      do_raman_clouds = .FALSE.
      addstr          = ""
    END IF

    ! --------------------------------------------
    ! Read scaling of cloud data fields (working?)
    ! --------------------------------------------
    scale_cfr = 1.0_r4 ; offset_cfr = 0.0_r4 ; missval_cfr = 0.0_r4
    locerrstat = HE5_SWrdlattr ( omicloud_swath_id, omicld_cfrac_field//TRIM(ADJUSTL(addstr)), 'MissingValue', missval_cfr )

    scale_ctp = 1.0_r4 ; offset_ctp = 0.0_r4 ; missval_ctp = 0.0_r4 ! 0
    locerrstat = HE5_SWrdlattr ( omicloud_swath_id, omicld_cpres_field//TRIM(ADJUSTL(addstr)), 'MissingValue', missval_ctp )

    ! -----------------------
    ! Read current data block
    ! -----------------------
    he5_start_2d = (/ 0, 0 /) ; he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)

    ! ----------------------------------------------------------------
    ! Read cloud fraction and cloud top pressure, and check for error.
    ! Eventually we may read the cloud uncertainties also, but for the
    ! first version we stick with just the basic cloud products.
    ! ----------------------------------------------------------------

    ! ---------------------------------------------------------------------------
    ! (1) Cloud Fraction is of type REAL*4 in both Raman and O2-O2 cloud products
    ! ---------------------------------------------------------------------------
    locerrstat = HE5_SWrdfld ( &
      omicloud_swath_id, omicld_cfrac_field//TRIM(ADJUSTL(addstr)),   &
      he5_start_2d, he5_stride_2d, he5_edge_2d, cfr(1:nx,0:nt-1) )
    IF ( locerrstat /= pge_errstat_ok ) then
      call tell_error (tell_io_read_error, "amf_read_omiclouds: reading cloud fraction", &
                       errstat)
      return
    endif

    ! ---------------------------------------------------------------
    ! Check for rebinned zoom data swath storage ("1-30" vs. "16-45")
    ! ---------------------------------------------------------------
    DO it = 0, nt-1
      IF ( do_szoom(it) .AND. &
        ALL ( cfr(gzoom_epix:nx,it) <= missval_cfr ) ) THEN
        cfr(gzoom_spix:gzoom_epix,it) = cfr(1:gzoom_npix,it)
        cfr(1:gzoom_spix-1,       it) = missval_cfr
      END IF
    END DO

    ! -----------------------------------------------------------
    ! Assign the cloud fraction array used in the AMF calculation
    ! -----------------------------------------------------------
    l2cfr = REAL(cfr, KIND=r8)
    WHERE ( cfr > r4_missval )
      l2cfr = l2cfr * scale_cfr + offset_cfr
    END WHERE

    ! ---------------------------------------------------------------------------
    ! (2) Cloud Pressure of type REAL*4 in Raman but INT*2 in O2-O2
    ! ---------------------------------------------------------------------------
    IF ( do_raman_clouds ) THEN
      locerrstat = HE5_SWrdfld ( &
        omicloud_swath_id, omicld_cpres_field//TRIM(ADJUSTL(addstr)),   &
        he5_start_2d, he5_stride_2d, he5_edge_2d, ctp(1:nx,0:nt-1) )
    ELSE
      locerrstat = HE5_SWrdfld ( &
        omicloud_swath_id, omicld_cpres_field,   &
        he5_start_2d, he5_stride_2d, he5_edge_2d, o4ctp(1:nx,0:nt-1) )

      ! -------------------------------------------
      ! Temporary copy of O4CTP to an R4 type array
      ! -------------------------------------------
      ctp(1:nx,0:nt-1) = REAL ( o4ctp(1:nx,0:nt-1), KIND=r4 )
    END IF

    IF ( locerrstat /= pge_errstat_ok ) then
      call tell_error (tell_io_read_error, "amf_read_omiclouds: reading cloud top pressure", &
                       errstat)
      return
    endif

    ! ---------------------------------------------------------------
    ! Check for rebinned zoom data swath storage ("1-30" vs. "16-45")
    ! ---------------------------------------------------------------
    DO it = 0, nt-1
      IF ( do_szoom(it) .AND. &
        ALL ( ctp(gzoom_epix:nx,it) <= REAL(missval_ctp, KIND=r4) ) ) THEN
        ctp(gzoom_spix:gzoom_epix,it)  = ctp(1:gzoom_npix,it)
        ctp(1:gzoom_spix-1,       it)  = REAL(missval_ctp, KIND=r4)
      END IF
    END DO

    ! ---------------------------------------------------------------
    ! Assign the cloud top pressure array used in the AMF calculation
    ! ---------------------------------------------------------------
    l2ctp  = REAL(ctp, KIND=r8)
    WHERE ( ctp > REAL(missval_ctp, KIND=r4) )
      l2ctp = l2ctp * scale_ctp + offset_ctp
    END WHERE

    ! ------------------------------------------------
    ! Force the cloud parameters into physical bounds.
    ! But make sure not to remove MissingValues.
    ! ------------------------------------------------
    WHERE ( l2cfr > REAL(missval_cfr, KIND=r8) .AND. l2cfr < 0.0_r8 )
      l2cfr = 0.0_r8
    ENDWHERE
    WHERE ( l2cfr > 1.0_r8 )
      l2cfr = 1.0_r8
    ENDWHERE
    WHERE ( l2ctp > REAL(missval_ctp, KIND=r8) .AND. l2ctp < 0.0_r8 )
      l2ctp = 0.0_r8
    ENDWHERE

    ! ------------------------------------------------------
    ! Replace cloud missing values by SAO PGE missing values
    ! ------------------------------------------------------
    WHERE ( l2cfr <= REAL(missval_cfr, KIND=r8) )
      l2cfr = r8_missval
    ENDWHERE
    WHERE ( l2ctp <= REAL(missval_ctp, KIND=r8) )
      l2ctp = r8_missval
    ENDWHERE

    RETURN
  END SUBROUTINE amf_read_omiclouds

  SUBROUTINE amf_diagnostic ( nt, nx, lat, lon, &
                             sza, vza, snow, glint, xtrange, &
                             l2cfr, l2ctp, amfdiag )

    USE OMSAO_omidata_module, ONLY: omi_geo_amf, omi_cld_addmiss, &
         omi_ooblut_amf, omi_glint_add, omi_bigsza_amf
    !USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza, lat, lon
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: snow, glint
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange

    ! ----------------
    ! Modified variabe
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: amfdiag
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: l2cfr, l2ctp

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: j1, j2, ix, it, ilat, ilon, spix, epix
    REAL    (KIND=r8) :: latdp, londp

    ! -------------------------------------------------------------------
    ! AMFDIAG has already been set to "geometric" AMF where SZA and VZA
    ! information was available, and "missing"/"out of bounds" otherwise.
    ! Here we need to check for cloud information as well as for snow
    ! and glint.
    ! -------------------------------------------------------------------

    ! ----------------------------------------------------------------
    ! Checking is done within a loop over NT to assure that we have
    ! "missing" values in all the right places. A single comprehensive
    ! WHERE statement over "1:nx,0:nt-1" would be more efficient but
    ! would also overwrite missing values with real diagnostic flags.
    ! ----------------------------------------------------------------

    ! -------------------------------------------------------------------
    ! If the file privided for Scattering weights has no information with
    ! in the fitting window, then nothing is to be done here.
    ! The geometric AMFs flag will remain in amfdiag and no further calcu
    ! lation will be performed inside the calculate_amf and calculate_sca
    ! subroutines.
    ! -------------------------------------------------------------------
    ! ----------------------------------------------------------
    ! Check that the value in amf_wvl is within the range of the
    ! suplied scattering weights file.
    ! ----------------------------------------------------------
    IF ((amf_wvl  .LT. MINVAL(lut_wav) ) .OR. &
        (amf_wvl2 .GT. MAXVAL(lut_wav) ) ) THEN
      RETURN
    END IF

    DO it = 0, nt-1
      spix = xtrange(it,1) ; epix = xtrange(it,2)

      ! ----------------------------------
      ! Check for ice/snow/ocean flag
      ! ----------------------------------
      WHERE ( &
           amfdiag(spix:epix,it) >= omi_geo_amf .AND. &
           snow(spix:epix,it) >= 0_i2         )
         amfdiag(spix:epix,it) = snow(spix:epix,it)
      END WHERE

      ! ---------------------------
      ! Check for glint possibility
      ! ---------------------------
      WHERE (                                    &
           amfdiag  (spix:epix,it) >= 0_i2 .AND. &
           glint(spix:epix,it) > 0_i2           )
         amfdiag(spix:epix,it) = amfdiag(spix:epix,it) + omi_glint_add
      END WHERE

      ! ---------------------------------------
      ! Check if we have cloud information from
      ! satellite retrievals. If not complete
      ! with cloud climatology.
      ! ---------------------------------------
       IF ( ( ANY( l2cfr(spix:epix,it) < 0.0_r8 ) ) .OR. &
            ( ANY( l2ctp(spix:epix,it) < 0.0_r8 ) ) ) THEN
          DO ix = spix, epix
             IF ((l2cfr(ix,it) < 0.0_r8 .or. l2ctp(ix,it) < 0.0_r8) .and. amfdiag(ix,it) >= 0_i2) THEN
                latdp = REAL ( lat(ix,it), KIND=r8 ) ; londp = REAL ( lon(ix,it), KIND=r8 ) ;
                ilat = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%latvals-latdp) ))
                j1 = SUM(ISCCP_CloudClim%n_lonvals(1:ilat-1)) + 1
                j2 = ISCCP_CloudClim%n_lonvals(ilat) + j1
                ilon = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%lonvals(j1:j2)-londp) ))
                l2ctp(ix,it) = ISCCP_CloudClim%ctp(ilon)
                l2cfr(ix,it) = ISCCP_CloudClim%cfr(ilon)                   
                amfdiag(ix,it) = omi_cld_addmiss + amfdiag(ix,it)
             ENDIF
          END DO
       END IF

      ! --------------------------------------------------
      ! Angles above the top value set on the control file
      ! are calculated "using this maximum value".
      ! --------------------------------------------------
      WHERE ( &
          ( sza(spix:epix,it) >= amf_max_sza) .AND. &
          ( amfdiag(spix:epix,it) >= 0_i2 ) )
        amfdiag(spix:epix,it) = omi_bigsza_amf + amfdiag(spix:epix,it)
      END WHERE

      ! ----------------------------------------------------------
      ! Out-of-Bound SZA,VZA,and cloud pressure (but not missing!)
      ! ----------------------------------------------------------
      WHERE( ( sza(spix:epix,it)   < MINVAL(lut_sza) ) .OR. &
             ( sza(spix:epix,it)   > MAXVAL(lut_sza) ) .OR. &
             ( vza(spix:epix,it)   < MINVAL(lut_vza) ) .OR. &
             ( vza(spix:epix,it)   > MAXVAL(lut_vza) ) .OR. &
             ( l2ctp(spix:epix,it) < MINVAL(lut_srf) ) .OR. &
             ( l2ctp(spix:epix,it) > MAXVAL(lut_srf) ) .AND. &
             ( amfdiag(spix:epix,it) >= 0_i2 ) )
         amfdiag(spix:epix,it) = omi_ooblut_amf + amfdiag(spix:epix,it)
      END WHERE

    END DO

    RETURN
  END SUBROUTINE amf_diagnostic

  SUBROUTINE compute_scatt ( nt, nx, albedo, sza, vza, saa, vaa, l2ctp, l2cfr, terrain_height, &
       surface_pressure, cli_wgh_ozo_pro, cli_idx_ozo_pro, lat, lon, amfdiag, scattw)

    use OMSAO_omidata_module, only: omi_scattfail_amf
    USE OMSAO_linterpolation_module, ONLY: lininterpol, GetNode
    USE ezspline_interpolation, ONLY: ezspline_2d_interpolation
    use sao_pge_utils, only: calc_relaz_angle

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nt, nx
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (inout) :: amfdiag
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza, saa, vaa, terrain_height, lat, lon
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: albedo, l2cfr
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,1:2), INTENT (IN) :: cli_wgh_ozo_pro
    INTEGER (KIND=i4), DIMENSION (1:nx,0:nt-1,1:2), INTENT (IN) :: cli_idx_ozo_pro

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT (INOUT) :: scattw
    REAL (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: l2ctp
    REAL (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (INOUT) :: surface_pressure

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: ialb, ictp, ilay, isrf, isza, ivza, itime, ixtrack, &
         iwavs, iwavf, nsza, nvza, nalb, ncld_alb, nsrf, nctp, nwav, &
         j1, j2, ilat, ilon, status
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

    ! -----------------------------------
    ! Find look up table wavelength index
    ! No interpolation, closest available
    ! is selected.
    ! -----------------------------------
    iwavs = MINLOC(ABS(lut_wav - REAL(amf_wvl,  KIND = r4) ),1)
    iwavf = MINLOC(ABS(lut_wav - REAL(amf_wvl2, KIND = r4) ),1)
    nwav = iwavf - iwavs + 1

    write(logmsg, '(a)')'Computing scattering weights...'
    call tell_log (1, logmsg)

    ! ---------------
    ! Loop over lines
    ! ---------------
    DO itime = 0, nt-1
       ! --------------------------
       ! Loop over xtrack positions
       ! --------------------------
       DO ixtrack = 1, nx

          IF (amfdiag(ixtrack,itime) .LT. 0) CYCLE

          ! ----------------------------------------------
          ! If this point is reached then scattw should be
          ! different from r8_missval and it needs to be
          ! initialized to 0.0 to work out the average
          ! ----------------------------------------------
          scattw(ixtrack,itime,:) = 0.0_r8

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
          ! ----------------------------------------------
          ! If sza > amf_max_sza set it for calculation to
          ! amf_max_sza
          ! ----------------------------------------------
          IF (local_sza .GT. amf_max_sza) local_sza = amf_max_sza
          !------------------------------------------------------
          ! We also need to check vza for values >max(vza_lut)
          ! or else flag pixel as bad and move on...
          ! -----------------------------------------------------
          ! reset value is arbitrary, but approved by Gonzalo
          IF (local_vza .GT. maxval(lut_vza)) local_vza = maxval(lut_vza)-0.1

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

          ! ------------------------------------------------------
          ! Bringing it to the lowest available pressure if needed
          ! ------------------------------------------------------
          IF (local_srf .GT. 1030.0) local_srf = 1030.0_r8
          ! Save it in ouptut variable
          surface_pressure(ixtrack,itime) = REAL(local_srf, KIND=r4)

          ! ---------------------------------------------------------------------------
          ! If cloud pressure is greater than surface pressure (cloud below surface!!!)
          ! then use climatology to correct cloud pressure.
          ! ---------------------------------------------------------------------------
          IF (local_ctp .GT. local_srf) THEN 
             ilat = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%latvals-REAL(lat(ixtrack,itime),KIND=r8))))
             j1 = SUM(ISCCP_CloudClim%n_lonvals(1:ilat-1)) + 1
             j2 = ISCCP_CloudClim%n_lonvals(ilat) + j1
             ilon = MAXVAL(MINLOC( ABS(ISCCP_CloudClim%lonvals(j1:j2)-REAL(lon(ixtrack,itime),KIND=r8))))
             IF (ISCCP_CloudClim%ctp(ilon) < local_srf) THEN
                local_ctp = ISCCP_CloudClim%ctp(ilon)
                l2ctp(ixtrack,itime) = ISCCP_CloudClim%ctp(ilon)
             ELSE
                local_ctp = local_srf
                l2ctp(ixtrack,itime) = local_srf
             END IF
          END IF

          ! -----------------------------------------------------------------
          ! If cloud pressure is lower than the minimum available in LUT then
          ! force it to match the lowest available pressure
          ! -----------------------------------------------------------------
          IF (local_ctp .LT. MINVAL(lut_srf)) local_ctp = REAL(MINVAL(lut_srf),KIND=r8)

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
            amfdiag(ixtrack,itime) = omi_scattfail_amf
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
            amfdiag(ixtrack,itime) = omi_scattfail_amf
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
                  amfdiag(ixtrack,itime) = omi_scattfail_amf
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
               amfdiag(ixtrack,itime) = omi_scattfail_amf
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
                  amfdiag(ixtrack,itime) = omi_scattfail_amf
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
               amfdiag(ixtrack,itime) = omi_scattfail_amf
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
          ! local_srf (from climatology or L1 file), Ap and Bp.
          ! --------------------------------------------------------------
          DO ilay = 1, CmETA
             out_pre_lay = (( Ap(ilay) + local_srf * Bp(ilay)  ) + &
                  ( Ap(ilay+1) + local_srf * Bp(ilay+1) )) / 2.0
             IF ( (out_pre_lay > MAXVAL(lut_pre_lay)) .OR. (out_pre_lay < MINVAL(lut_pre_lay)) ) THEN
                scattw(ixtrack,itime,CmETA-ilay+1) = 0.0
                cycle
             ENDIF
             scattw(ixtrack,itime,CmETA-ilay+1) = linInterpol( (INT(lay_dim(1),KIND=i4)), REAL(LOG(lut_pre_lay),KIND=r8), &
                  Sca_1D, LOG(out_pre_lay), status=status)
             IF ( status /= 0 ) THEN
               amfdiag(ixtrack,itime) = omi_scattfail_amf
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
          WHERE ( scattw(ixtrack,itime,1:CmETA) < 0.0_r8 )
             scattw(ixtrack,itime,1:CmETA) = 0.0_r8
          END WHERE
          
       END DO ! End loop xtrack

    END DO ! End loop lines
 
  END SUBROUTINE COMPUTE_SCATT

  SUBROUTINE compute_amf ( nt, nx, CmETA, climatology, &
      scattw, saoamf, stratospheric_amf, tropospheric_amf, surface_pressure, &
      tropopause_pressure, lat, lon, amfdiag, errstat)

    use ctrlvars, only: yn_stratrop
    use met_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT(IN) :: nt, nx, CmETA
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1,CmETA), INTENT(IN) :: climatology, scattw
    INTEGER (KIND=i2), DIMENSION(1:nx,0:nt-1), INTENT(IN) :: amfdiag
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
    INTEGER (KIND=i4) :: ixtrack, itimes, ilay, tropopause_idx
    REAL (KIND=r8), DIMENSION(:), ALLOCATABLE :: pressure_grid, temperature_profile, alpha

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=11), PARAMETER :: modulename = 'compute_amf'

    if (errstat /= 0) return

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    DO itimes  = 0, nt-1 ! Swath lines loop
      DO ixtrack = 1, nx ! Xtrack pixel loop

        ! ------------------------------------------------
        ! Set all amf with amfdiag < 0 to r8_missval
        ! ------------------------------------------------
        IF ( amfdiag(ixtrack,itimes) .LT. 0 ) THEN
          saoamf(ixtrack,itimes) = r8_missval
          CYCLE
        ENDIF

        ! ---------------------------------------------------
        ! Read tropopause pressure from met forecast file
        ! ---------------------------------------------------
        IF (yn_stratrop) THEN
           ! Allocate pressure_grid and temperature vertical profile
           ALLOCATE(pressure_grid(1:CmETA),temperature_profile(1:CmETA), &
                alpha(1:CmETA))
           ! Read in tropopause pressure
           call read_synth_met_data(trim(OMSAO_meteorology_filename), &
                lat(ixtrack,itimes), lon(ixtrack,itimes), &
                tropopause_pressure(ixtrack,itimes), errstat=errstat)

           ! Set temperature profile constant and equal to 220K
           temperature_profile = 220.0_r8
           ! Work out pressure_grid
           DO ilay = 1, CmETA
              pressure_grid(CmETA-ilay+1) = (( Ap(ilay) + surface_pressure(ixtrack,itimes) * Bp(ilay)  ) + &
                   ( Ap(ilay+1) + surface_pressure(ixtrack,itimes) * Bp(ilay+1) )) / 2.0
           END DO

           ! Find which layer is closer to the tropopause.
           tropopause_idx = MINLOC(ABS(pressure_grid-REAL(tropopause_pressure(ixtrack,itimes),KIND=r4)),1)

           ! Compute stratospheric and tropospheric AMFs following Bucsela et al., 2013
           ! DOI:10.5194/amt-6-2607-2013
           ! Apply temperature correction factor alpha(p) = 1-0.003 [T(p)-T0] with T0 .EQ. 220K
           ! EJOS adding a test for zero in climatology to avoid NaN AMFs
           alpha = 1.0_r8-0.003_r8*(temperature_profile-220.0_r8)
           if (SUM(climatology(ixtrack,itimes,1:tropopause_idx)).eq.0) then
             stratospheric_amf(ixtrack,itimes) = 0.0d0
           else
             stratospheric_amf(ixtrack,itimes) = SUM(scattw(ixtrack, itimes, 1:tropopause_idx) * &
                  climatology(ixtrack,itimes,1:tropopause_idx) * alpha(1:tropopause_idx))     / &
                  SUM(climatology(ixtrack,itimes,1:tropopause_idx))
           endif
           if (SUM(climatology(ixtrack,itimes,tropopause_idx+1:CmETA)).eq.0) then
             tropospheric_amf(ixtrack,itimes) = 0.0d0
           else
             tropospheric_amf(ixtrack,itimes) = SUM(scattw(ixtrack, itimes, tropopause_idx+1:CmETA) * &
                  climatology(ixtrack,itimes,tropopause_idx+1:CmETA) * alpha(tropopause_idx+1:CmETA) ) / &
                  SUM(climatology(ixtrack,itimes,tropopause_idx+1:CmETA))
           endif
           DEALLOCATE(pressure_grid,temperature_profile,alpha)
        END IF

        ! -------------------------
        ! Finally work out the AMFs
        ! -------------------------
        saoamf(ixtrack,itimes) = SUM(scattw(ixtrack, itimes, 1:CmETA) * &
             climatology(ixtrack,itimes,1:CmETA))     / &
             SUM(climatology(ixtrack,itimes,1:CmETA))

      END DO ! Finish xtrack pixel loop
    END DO ! Finish

  END SUBROUTINE compute_amf

  SUBROUTINE write_albedo_he5(albedo, nt, nx, errstat)

    ! ==================================================================
    ! This routines writes the albedos obtained from the OMLER climatolo
    ! gy to the output file.
    ! ==================================================================
    use datafields, only: albedo_field
    use OMSAO_he5_module, ONLY: HE5_SWWRFLD
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                         INTENT (IN) :: nt, nx
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (IN) :: albedo

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=16), PARAMETER :: modulename = 'write_albedo_he5'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                          :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc

    !locerrstat = pge_errstat_ok

    he5_start_2d  = (/ 0, 0 /)
    he5_stride_2d = (/ 1, 1 /)
    he5_edge_2d   = (/ nx, nt /)

    colloc = albedo
    !CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                              TRIM(ADJUSTL(albedo_field)),             &
                              he5_start_2d, he5_stride_2d, he5_edge_2d,&
                              colloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

  END SUBROUTINE write_albedo_he5

  SUBROUTINE write_climatology_he5(climatology, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the Target Gas Profiles from the GEOS-Chem
    ! climatology to the output file.
    ! ===============================================================
    use datafields, only: gasprofile_field
    use OMSAO_he5_module, ONLY: HE5_SWWRFLD
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                              INTENT (IN) :: nt, nx, nl
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1,1:nl), INTENT (IN) :: climatology

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=24), PARAMETER :: modulename = 'write_climatology_he5'  ! JED fix

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                               :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,1:nl) :: colloc

    !locerrstat = pge_errstat_ok

    he5_start_3d  = (/ 0, 0, 0 /)
    he5_stride_3d = (/ 1, 1, 1 /)
    he5_edge_3d   = (/ nx, nt, nl /)

    colloc = climatology
    !CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                              TRIM(ADJUSTL(gasprofile_field)),         &
                              he5_start_3d, he5_stride_3d, he5_edge_3d,&
                              colloc(1:nx,0:nt-1,1:nl) )
    errstat = MAX ( errstat, locerrstat )

  END SUBROUTINE write_climatology_he5

  SUBROUTINE write_scatt_he5(scattw, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the scattering weigths to the output file.
    ! ===============================================================
    use datafields, only: scaweights_field
    use OMSAO_he5_module, ONLY: HE5_SWWRFLD

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                              INTENT (IN) :: nt, nx, nl
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1,1:nl), INTENT (IN) :: scattw !, akernels

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4),                         INTENT (INOUT) :: errstat

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=16), PARAMETER :: modulename = 'write_scatt_he5'

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                               :: locerrstat
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,1:nl) :: colloc

    he5_start_3d  = (/ 0, 0, 0 /)
    he5_stride_3d = (/ 1, 1, 1 /)
    he5_edge_3d   = (/ nx, nt, nl /)

    colloc = scattw
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                              TRIM(ADJUSTL(scaweights_field)),         &
                              he5_start_3d, he5_stride_3d, he5_edge_3d,&
                              colloc(1:nx,0:nt-1,1:nl) )

    errstat = MAX ( errstat, locerrstat )

  END SUBROUTINE write_scatt_he5

  SUBROUTINE he5_amf_write ( &
      nx, nt, saocol, saodco, amfmol, amfgeo, amfdiag, &
      amfcfr, amfctp, surface_pressure, tropopause_pressure, &
      tropospheric_amf, stratospheric_amf, errstat )

    USE OMSAO_precision_module, ONLY: i2, i4, r8
    USE OMSAO_he5_module, ONLY: HE5_SWWRFLD, he5_start_2d, he5_stride_2d, &
      he5_edge_2d
    USE OMSAO_errstat_module, only : pge_errstat_ok
    use datafields, only: amfcfr_field, amfctp_field, amfdiag_field, &
      amfgeo_field, amfmol_field, col_field, dcol_field, surpre_field, &
      amfstr_field, amftro_field, tropre_field
    use ctrlvars, only: yn_stratrop

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=13), PARAMETER :: modulename = 'he5_write_amf'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nx, nt
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: saocol, saodco
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfmol, amfgeo, tropospheric_amf, &
         stratospheric_amf
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfcfr, amfctp
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: surface_pressure, tropopause_pressure
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfdiag

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                          :: locerrstat
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1) :: amfloc
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc

    ! -------------------------------------------------------------
    ! Air mass factor plus diagnostic.
    ! -------------------------------------------------------------
    ! As yet, only OMBRO and OMHCHO have true, non-geometric AMFs.
    ! But we try to have symmetric data fields as much as possible,
    ! hence the presence of the "molecule specific" AMF and its
    ! diagnostic for all PGEs. Non-OMBRO and -OMHCHO PGEs carry a
    ! geometric AMF here.
    !
    ! For completeness, the geometric AMF is added.
    ! -------------------------------------------------------------

    locerrstat = pge_errstat_ok

    he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)

    ! ----------------------------------------
    ! (1) AMF diagnostic. No rounding required
    ! ----------------------------------------
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfdiag_field)), &
      he5_start_2d, he5_stride_2d, he5_edge_2d, amfdiag(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    ! -----------------
    ! (2) Geometric AMF
    ! -----------------
    amfloc = REAL ( amfgeo, KIND=r4 )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfgeo_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    ! -----------------
    ! (3) Molecular AMF
    ! -----------------
    amfloc = REAL ( amfmol, KIND=r4 )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfmol_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    IF (yn_stratrop) then
       amfloc = REAL ( tropospheric_amf, KIND=r4 )
       locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amftro_field)), &
                                  he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
       errstat = MAX ( errstat, locerrstat )
       amfloc = REAL ( stratospheric_amf, KIND=r4 )
       locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfstr_field)), &
                                  he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
       errstat = MAX ( errstat, locerrstat )
    END IF

    ! ----------------------------------------------------------
    ! (4) AMF cloud fraction and pressure
    ! ----------------------------------------------------------
    amfloc = REAL ( amfcfr, KIND=r4 )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfcfr_field)), &
         he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )
    
    amfloc = REAL ( amfctp, KIND=r4 )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfctp_field)), &
         he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    ! -----------------------------------------------------------------------
    ! (5) All PGEs: Output of columns and column uncertainties. For some PGEs
    !     (e.g., OMBRO, OMHCHO, OMCHOCHO) those have been adjusted by the AMF,
    !     but we have as yet to perform the rounding for any of them.
    ! -----------------------------------------------------------------------
    colloc = saocol
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(col_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    colloc = saodco
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(dcol_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    ! ----------------
    ! Surface pressure
    ! ----------------
    amfloc = surface_pressure
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(surpre_field)), &
         he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    ! -------------------
    ! Tropopause pressure
    ! -------------------
    IF (yn_stratrop) then
       amfloc = tropopause_pressure
       locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(tropre_field)), &
                                 he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
       errstat = MAX ( errstat, locerrstat )
    END IF

    RETURN
  END SUBROUTINE he5_amf_write

  SUBROUTINE read_climatology_dimensions(errstat)

    USE OMSAO_he5_module, ONLY: he5_swopen, he5_swattach, &
      he5_SWrdlattr, he5_swclose, he5f_acc_rdonly, &
      he5_swinqswath, granule_month
    USE OMSAO_errstat_module, ONLY: he5_stat_fail

    IMPLICIT NONE
    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4) :: nswath, swath_id, swath_file_id
    CHARACTER (LEN=   MAX_STR_LEN) :: swath_file, locswathname
    INTEGER   (KIND=C_LONG) :: nswathcl, swlen
    CHARACTER (LEN=10*MAX_STR_LEN) :: swath_name

    swath_file = TRIM(ADJUSTL(OMSAO_climatology_filename))

    ! --------------------------------------------------------------
    ! Open he5 OMI climatology and check SWATH_FILE_ID (-1 if error)
    ! --------------------------------------------------------------
    swath_file_id = HE5_SWOPEN (swath_file, he5f_acc_rdonly)
    IF (swath_file_id == he5_stat_fail) THEN
      call tell_error (tell_io_open_error, &
           "read_climatology_dimensions: opening climatology file"//trim(swath_file), &
           errstat)
      RETURN
    END IF

    ! -----------------------------------------------------------
    ! Check for existing HE5 swathw and attach to the one we need
    ! -----------------------------------------------------------
    swath_name = "" !JED
    nswathcl = HE5_SWinqswath(TRIM(ADJUSTL(swath_file)), swath_name, swlen )
    nswath   = INT(nswathcl, KIND=i4 )

    ! ----------------------------------------------------------------
    ! If there is only one swath in the file, we can attach to it but
    ! if there are more (NSWATH > 1), then we must find the swath that
    ! corresponds to the current month.
    ! ----------------------------------------------------------------
    IF (nswath > 1) THEN
      CALL extract_swathname(nswath, TRIM(ADJUSTL(swath_name)), &
        TRIM(ADJUSTL(months(granule_month))), locswathname)
      ! ---------------------------------------------------------------------------
      ! Check if we found the correct swath name. If not, report an error and exit.
      ! ---------------------------------------------------------------------------
      IF ( INDEX (TRIM(ADJUSTL(locswathname)),TRIM(ADJUSTL(months(granule_month)))) == 0 ) THEN
        call tell_error (tell_runtime_error, "read_climatology_dimensions: finding month in "// &
                         trim(locswathname), errstat)
        RETURN
      END IF
    ELSE
      locswathname = TRIM(ADJUSTL(swath_name))
    END IF

    ! -----------------------------
    ! Attach to current month swath
    ! -----------------------------
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(locswathname)) )
    IF ( swath_id == he5_stat_fail ) THEN
      call tell_error (tell_io_error, "read_climatology_dimensions: attaching to swath "// &
                       trim(locswathname), errstat)
      RETURN
    END IF

    ! ------------------------------------
    ! Read dimensions of Climatology swath
    ! ------------------------------------
    CALL climatology_getdim ( swath_id, errstat )
    if (errstat /= 0) return

    errstat = HE5_SWCLOSE(swath_file_id)
    IF ( errstat == he5_stat_fail ) THEN
      call tell_error (tell_io_error, "read_climatology_dimensions: closing climatology file "// &
                       trim(swath_file), errstat)
      RETURN
    END IF

  END SUBROUTINE read_climatology_dimensions

  SUBROUTINE voc_amf_readisccp ( errstat )

    USE OMSAO_variables_module, ONLY: voc_amf_filenames
    USE OMSAO_indices_module, ONLY: voc_isccp_idx
    USE OMSAO_errstat_module, ONLY: he5_stat_fail
    USE OMSAO_he5_module, ONLY: HE5_SWclose, HE5_SWopen, &
         he5f_acc_rdonly, granule_month, HE5_SWattach, &
         HE5_SWinqswath, HE5_SWrdfld, HE5_SWrdlattr, &
         HE5_SWdetach

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=C_LONG) :: swath_id, swath_file_id, swlen, he5stat
    INTEGER   (KIND=C_LONG) :: locerrstat
    CHARACTER (LEN=MAX_STR_LEN) :: swath_file, swath_name
    INTEGER   (KIND=i4), DIMENSION (nlon_isccp) :: tmparr

    swath_file = TRIM(ADJUSTL(voc_amf_filenames(voc_isccp_idx)))
    
    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    swath_file_id = HE5_SWopen ( swath_file, he5f_acc_rdonly )
    IF ( swath_file_id == he5_stat_fail ) THEN
      call tell_error (tell_io_error, "voc_amf_readisccp: opening climatology file "// &
                       trim(swath_file), errstat)
      RETURN
    END IF
    
    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    swath_name = ''
    locerrstat  = HE5_SWinqswath  ( TRIM(ADJUSTL(swath_file)), swath_name, swlen )
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id == he5_stat_fail ) THEN
       call tell_error (tell_io_error, "voc_amf_readisccp: attaching to swath "// &
            trim(swath_name), errstat)
       RETURN
    END IF

    ! ----------------------------
    ! Read ISCCP Cloud Climatoloty
    ! ----------------------------
    ! * Latitudes
    he5_start_1d = 0 ; he5_stride_1d = 1 ; he5_edge_1d = nlat_isccp
    he5stat = HE5_SWrdfld ( swath_id, isccp_lat_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%latvals(1:nlat_isccp) )
    ! * Number of longitudes per latitude
    he5stat = HE5_SWrdfld ( swath_id, isccp_nlon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%n_lonvals(1:nlat_isccp) )
    ! * Delta-Longitudes    
    he5stat = HE5_SWrdfld ( swath_id, isccp_dlon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%delta_lon(1:nlat_isccp) )
    ! * Longitudes
    he5_start_1d = 0 ; he5_stride_1d = 1 ; he5_edge_1d = nlon_isccp
    he5stat = HE5_SWrdfld ( swath_id, isccp_lon_field, &
         he5_start_1d, he5_stride_1d, he5_edge_1d, ISCCP_CloudClim%lonvals(1:nlon_isccp) )
    
    ! * Cloud fraction
    he5_start_2d  = (/ granule_month-1, 0 /) ; he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ 1, nlon_isccp /)
    he5stat = HE5_SWrdfld ( swath_id, isccp_mcfr_field,                 &
         he5_start_2d, he5_stride_2d, he5_edge_2d, tmparr(1:nlon_isccp) )
    ISCCP_CloudClim%cfr(1:nlon_isccp) = REAL (tmparr(1:nlon_isccp), KIND=r8)
    ! * Cloud top pressure
    he5stat = HE5_SWrdfld ( swath_id, isccp_mctp_field,                 &
         he5_start_2d, he5_stride_2d, he5_edge_2d, tmparr(1:nlon_isccp) )
    ISCCP_CloudClim%ctp(1:nlon_isccp) = REAL (tmparr(1:nlon_isccp), KIND=r8)
    
    ! --------------------
    ! Read some attributes
    ! --------------------
    he5stat = HE5_SWrdlattr ( swath_id, isccp_lat_field,  "DeltaGrid",    ISCCP_CloudClim%delta_lat   )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mcfr_field, "ScaleFactor",  ISCCP_CloudClim%scale_cfr   )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mctp_field, "ScaleFactor",  ISCCP_CloudClim%scale_ctp   )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mcfr_field, "MissingValue", ISCCP_CloudClim%missval_cfr )
    he5stat = HE5_SWrdlattr ( swath_id, isccp_mctp_field, "MissingValue", ISCCP_CloudClim%missval_ctp )

    ! -------------------------------------------------------
    ! Scale the ISCCP cloud values with their scaling factors
    ! -------------------------------------------------------
    WHERE ( ISCCP_CloudClim%cfr /= ISCCP_CloudClim%missval_ctp )
       ISCCP_CloudClim%cfr = ISCCP_CloudClim%cfr * ISCCP_CloudClim%scale_cfr
    END WHERE
    WHERE ( ISCCP_CloudClim%ctp /= ISCCP_CloudClim%missval_ctp )
       ISCCP_CloudClim%ctp = ISCCP_CloudClim%ctp * ISCCP_CloudClim%scale_ctp
    END WHERE
   
    ! -----------------------------------------------
    ! Detach from HE5 swath and close HE5 output file
    ! -----------------------------------------------
    he5stat = HE5_SWdetach ( swath_id )
    he5stat = HE5_SWclose  ( swath_file_id )

    RETURN
  END SUBROUTINE voc_amf_readisccp

END MODULE OMSAO_wfamf_module

