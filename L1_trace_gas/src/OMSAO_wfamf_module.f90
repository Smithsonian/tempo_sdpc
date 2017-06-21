MODULE OMSAO_wfamf_module

  ! ====================================================================
  ! This module defines variables associated with the wavelength depende
  ! nt AMF calculations and contains necessary subroutines to read files
  ! and calculate them
  ! ====================================================================
  USE OMSAO_precision_module, ONLY: i2, i4, r8, C_LONG, r4
  USE OMSAO_parameters_module, ONLY: MAX_STR_LEN, i2_missval, r4_missval, r8_missval
  use tell_module
  USE OMSAO_he5_module, ONLY: pge_swath_id, &
    he5_start_4d, he5_edge_4d, he5_stride_4d, &
    he5_start_3d, he5_edge_3d, he5_stride_3d, &
    he5_start_2d, he5_edge_2d, he5_stride_2d, &
    he5_start_1d, he5_edge_1d, he5_stride_1d

  IMPLICIT NONE
  private

  public read_climatology_dimensions, amf_calculation_bis, &
    wfamf_deallocate

  ! ---------
  ! PCF stuff
  ! ---------
  INTEGER(KIND=i4), PARAMETER, public :: wfamf_table_lun = 700250
  INTEGER(KIND=i4), PARAMETER, public :: climatology_lun = 700270
  CHARACTER(LEN=MAX_STR_LEN), public  :: OMSAO_wfamf_table_filename
  CHARACTER(LEN=MAX_STR_LEN), public  :: OMSAO_climatology_filename

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
  REAL(KIND=r4), DIMENSION(:,:,:), ALLOCATABLE :: Psurface
  REAL(KIND=r4), DIMENSION(:,:,:,:), ALLOCATABLE :: Temperature, Gas_profiles, H2O_profiles

  ! ---------------------------------------
  ! Data obtained from Vlidort lookup table
  ! ---------------------------------------
  ! --------------------------------------------------
  ! Parameter for the definition of the vlidort arrays
  ! --------------------------------------------------
  ! ------------------------
  ! Cross sections variables
  ! ------------------------
  REAL(KIND=r4), DIMENSION(:), ALLOCATABLE :: vl_OzC0, vl_OzC1, vl_OzC2

  ! --------------
  ! Grid variables
  ! --------------
  REAL(KIND=r4),    DIMENSION(:), ALLOCATABLE :: vl_pre
  REAL(KIND=r4),    DIMENSION(:), ALLOCATABLE :: vl_sza
  REAL(KIND=r4),    DIMENSION(:), ALLOCATABLE :: vl_vza
  REAL(KIND=r4),    DIMENSION(:), ALLOCATABLE :: vl_wav
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: vl_toms

  ! ------------------
  ! Profiles variables
  ! ------------------
  REAL(KIND=r4), DIMENSION(:,:,:), ALLOCATABLE :: vl_air, vl_alt, vl_ozo, vl_tem

  ! --------------------------
  ! Parameterization variables
  ! --------------------------
  REAL(KIND=r4), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: vl_I0, vl_I1, vl_I2, vl_Ir
  REAL(KIND=r4), DIMENSION(:,:,:),       ALLOCATABLE :: vl_Sb
  REAL(KIND=r4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: vl_dI0, vl_dI1, vl_dI2, vl_dIr
  REAL(KIND=r4)                                      :: vl_Factor

  !!$  INTEGER(KIND=i4), PARAMETER :: vl_wavmax = 100, &
  !!$                                 vl_premax =   6, &
  !!$                                 vl_szamax =  12, &
  !!$                                 vl_vzamax =   8, &
  !!$                                 vl_altmax =  73, &
  !!$                                 vl_ozomax =  26
  !!$  ! ------------------------
  !!$  ! Cross sections variables
  !!$  ! ------------------------
  !!$  REAL(KIND=r4), DIMENSION(vl_wavmax) :: vl_OzC0, vl_OzC1, vl_OzC2
  !!$
  !!$  ! --------------
  !!$  ! Grid variables
  !!$  ! --------------
  !!$  REAL(KIND=r4),    DIMENSION(vl_premax) :: vl_pre
  !!$  REAL(KIND=r4),    DIMENSION(vl_szamax) :: vl_sza
  !!$  REAL(KIND=r4),    DIMENSION(vl_vzamax) :: vl_vza
  !!$  REAL(KIND=r4),    DIMENSION(vl_wavmax) :: vl_wav
  !!$  CHARACTER(LEN=4), DIMENSION(vl_ozomax) :: vl_toms
  !!$
  !!$  ! ------------------
  !!$  ! Profiles variables
  !!$  ! ------------------
  !!$  REAL(KIND=r4), DIMENSION(vl_ozomax,vl_premax,vl_altmax) :: vl_air, vl_alt, vl_ozo, vl_tem
  !!$
  !!$  ! --------------------------
  !!$  ! Parameterization variables
  !!$  ! --------------------------
  !!$  REAL(KIND=r4), DIMENSION(vl_ozomax,vl_premax,vl_szamax,vl_vzamax,vl_wavmax)           :: vl_I0, &
  !!$                                                                                           vl_I1, &
  !!$                                                                                           vl_I2, &
  !!$                                                                                           vl_Ir
  !!$  REAL(KIND=r4), DIMENSION(vl_ozomax,vl_premax,vl_wavmax)                               :: vl_Sb
  !!$  REAL(KIND=r4), DIMENSION(vl_ozomax,vl_premax,vl_szamax,vl_vzamax,vl_wavmax,vl_altmax) :: vl_dI0,&
  !!$                                                                                           vl_dI1,&
  !!$                                                                                           vl_dI2,&
  !!$                                                                                           vl_dIr
  !!$  REAL(KIND=r4)                                                                         :: vl_Factor

  ! -------------------
  ! Dimension variables
  ! -------------------
  INTEGER(KIND=i4) :: vl_nozo, vl_ncld, vl_nsza, vl_nvza, vl_nwav, vl_nalt

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
      pge_idx, nt, nx, lat, lon, sza, vza, time,  &
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
    USE OMSAO_errstat_module, only : pge_errstat_ok!, pge_errstat_error
    use OMSAO_indices_module, only: pge_hcho_idx, pge_gly_idx, voc_omicld_idx
    use OMSAO_omidata_module, only : amf_correction_type
    use output_tools, only : write_albedo, write_gas_profile, &
      write_scattering_weights, write_amf_correction
    USE OMSAO_variables_module,  ONLY: voc_amf_filenames
    use output_tools, only: read_cloud_params
    use ctrlvars, only : yn_do_he5_output
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx, pge_idx
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon, sza, vza, terrain_height
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
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), target :: amfgeo
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), target :: l2cfr, l2ctp
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1)       :: albedo, cli_psurface
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: climatology, cli_temperature, cli_heights
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA) :: scattw
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
    cli_heights  = r8_missval
    cli_psurface = r8_missval
    scattw       = r8_missval
    saoamf       = r8_missval
    amfgeo       = r8_missval
    amfdiag      = i2_missval

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
          CALL write_albedo_he5 ( albedo, nt, nx, locerrstat)! <-- FIXME: (to be removed)
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

      ! ------------------------------------------------
      ! Read climatology and interpolate to lon/lat/time
      ! ------------------------------------------------
      CALL omi_climatology (pge_idx, climatology, cli_heights, cli_psurface, cli_temperature, lat, lon, &
        nt, nx, xtrange, errstat)

      ! -------------------------------------
      ! Write the climatology to the he5 file
      ! -------------------------------------
      IF (do_write) then
        if (yn_do_he5_output) then
          CALL write_climatology_he5 (climatology, cli_heights, nt, nx, CmETA, locerrstat) ! <-- FIXME: (to be removed)
        endif
        call write_gas_profile (climatology, cli_heights, nx, nt, CmETA, errstat)
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
      CALL amf_diagnostic (nt, nx, &! lat, lon,
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
      CALL compute_scatt ( nt, nx, albedo, sza, vza, l2ctp, l2cfr, &
                          terrain_height, cli_heights, amfdiag, &
                          scattw)

      ! -----------------------------------------------------------------
      ! Work out the AMF using the scattering weights and the climatology
      ! Work out Averaging Kernels
      ! -----------------------------------------------------------------
      CALL compute_amf ( nt, nx, CmETA, climatology, &
           scattw, saoamf, amfdiag, locerrstat)

      ! -----------------------------------------------------------------
      ! Write out scattering weights, altitude grid and averaging kernels
      ! -----------------------------------------------------------------
      IF (do_write) then
        if (yn_do_he5_output) then
          CALL write_scatt_he5 (scattw, nt, nx, CmETA, locerrstat) ! FIXME <-- (to be removed)
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
        CALL he5_amf_write ( pge_idx, nx, nt, saocol, saodco, saoamf, &
                            amfgeo, amfdiag, l2cfr, l2ctp, locerrstat ) ! FIXME <-- (to be removed)
      endif
      amf_corr % amf_molecule_specific => saoamf
      amf_corr % amf_geometric => amfgeo
      amf_corr % diagnostic_flag => amfdiag
      amf_corr % cloud_fraction => l2cfr
      amf_corr % cloud_pressure => l2ctp
      yn_write_cloud_variables = (pge_idx == pge_hcho_idx) .or. (pge_idx == pge_gly_idx)
      call write_amf_correction (nx, nt, amf_corr, saocol, saodco, &
                                 yn_write_cloud_variables, errstat)
      if (errstat /= 0) return
    endif

  END SUBROUTINE amf_calculation_bis

  SUBROUTINE omi_climatology (pge_idx, climatology, local_heights, local_psurf, local_temperature, &
       lat, lon, nt, nx, xtrange, errstat)
    
    ! =========================================
    ! Extract Gas climatology to granule pixels
    ! No interpolation or something like that,
    ! Just pick the closest model grid
    ! =========================================
    USE OMSAO_indices_module, ONLY: sao_molecule_names, pge_h2o_idx
    USE OMSAO_omidata_module, ONLY: omi_oob_cli, omi_time, omi_time_utc
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
    INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
    
    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1, CmETA), INTENT (INOUT) :: climatology, &
         local_temperature, local_heights
    REAL (KIND=r8), DIMENSION(1:nx,0:nt-1), INTENT (INOUT) :: local_psurf
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: itimes, ixtrack, spix, epix, n, status , &
         nlon, nlat, ntim, locerrstat
    INTEGER (KIND=i4), DIMENSION(2) :: idx_lat, idx_lon, idx_tim
    REAL    (KIND=r8) :: rho, lhgt, aircolumn
    REAL    (KIND=r8), DIMENSION (1:CmETA) :: lh2o, lgas
    REAL    (KIND=r8) :: thish2omxr, Mwet, Rwet, detlnp, llon, llat, ltime

    INTEGER (KIND=i4) :: swath_file_id, swath_id, nswath, ndatafields
    INTEGER   (KIND=i4), DIMENSION(10) :: datafield_rank, datafield_type
    CHARACTER (LEN=MAX_STR_LEN) :: swath_file, swath_name, locswathname, datafield_name, &
         gasdatafieldname, h2odatafieldname
    INTEGER (KIND=C_LONG) :: nswathcl, swlen

    ! -----------------------
    ! Some physical constants
    ! -----------------------
    REAL (KIND=r8), PARAMETER ::  &
         Mdry = 0.02896,   & !kg mol-1
         Mh2o = 0.018,     & !kg mol-1
         Rstar = 8.314,    & !N m mol-1 K-1
         Navogadro = 6.02214e+23, & ! mol-1
         gplanet = 9.806,  &  ! m s-1
         m2tocm2 = 1.0e+4

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
       
       spix = xtrange(itimes,1); epix = xtrange(itimes,2)
       CALL convert_tai_to_utc(nUTCdim, omi_time(itimes), omi_time_utc(1:nUTCdim,itimes))
       
       DO ixtrack = spix, epix
          
          llon = REAL(lon(ixtrack,itimes),KIND=r8)
          llat = REAL(lat(ixtrack,itimes),KIND=r8)
          ltime = REAL(omi_time_utc(4,itimes),KIND=r8)
          
          IF (llon .LT. MINVAL(lonvals)) llon = REAL(MINVAL(lonvals),KIND=r8)
          IF (llon .GT. MAXVAL(lonvals)) llon = REAL(MAXVAL(lonvals),KIND=r8)
          IF (llat .LT. MINVAL(latvals)) llat = REAL(MINVAL(latvals),KIND=r8)
          IF (llat .GT. MAXVAL(latvals)) llat = REAL(MAXVAL(latvals),KIND=r8)

          ! ------------------------------------------------------------------------
          ! Given the values of lonvals, latvals, timevals and llon, llat, and ltime
          ! determine the indices of climatolgoy values to be read.
          ! Using linear interpolation only 2 nodes needed in each dimension or if 
          ! outbounds, closest node is selected
          ! ------------------------------------------------------------------------
          CALL GetNode(REAL(lonvals,KIND=r8),llon, &
               idx_lon(1), 'Lower')
          IF (idx_lon(1) .EQ. -2) idx_lon(1) = 1
          IF (idx_lon(1) .EQ. -3) idx_lon(1) = Cmlon
          CALL GetNode(REAL(lonvals,KIND=r8),llon, &
               idx_lon(2), 'Upper')
          IF (idx_lon(2) .EQ. -2) idx_lon(2) = 1
          IF (idx_lon(2) .EQ. -3) idx_lon(2) = Cmlon
          nlon = idx_lon(2)-idx_lon(1) + 1

          CALL GetNode(REAL(latvals,KIND=r8),llat, &
               idx_lat(1), 'Lower')
          IF (idx_lat(1) .EQ. -2) idx_lat(1) = 1
          IF (idx_lat(1) .EQ. -3) idx_lat(1) = Cmlat
          CALL GetNode(REAL(latvals,KIND=r8),llat, &
               idx_lat(2), 'Upper')
          IF (idx_lat(2) .EQ. -2) idx_lat(2) = 1
          IF (idx_lat(2) .EQ. -3) idx_lat(2) = Cmlat
          nlat = idx_lat(2)-idx_lat(1) + 1

          CALL GetNode(REAL(timevals,KIND=r8),ltime, &
               idx_tim(1), 'Lower')
          IF (idx_tim(1) .EQ. -2) idx_tim(1) = 1
          IF (idx_tim(1) .EQ. -3) idx_tim(1) = CmHRS
          CALL GetNode(REAL(timevals,KIND=r8),ltime, &
               idx_tim(2), 'Upper')
          IF (idx_tim(2) .EQ. -2) idx_tim(2) = 1
          IF (idx_tim(2) .EQ. -3) idx_tim(2) = CmHRS
          ntim = idx_tim(2)-idx_tim(1) + 1

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
               swath_id, gasdatafieldname, h2odatafieldname, locerrstat)
          IF ( locerrstat /= 0 ) THEN
             call tell_error (tell_runtime_error, "omi_climatology: read climatology failed", errstat)
             RETURN
          END IF

          ! ----------------------------
          ! Interpolate Surface pressure
          ! ----------------------------
          local_psurf(ixtrack,itimes) = linInterpol(nlon,nlat,ntim, &
               REAL(lonvals(idx_lon(1):idx_lon(2)),KIND=r8), &
               REAL(latvals(idx_lat(1):idx_lat(2)),KIND=r8), &
               REAL(timevals(idx_tim(1):idx_tim(2)),KIND=r8), &
               REAL(Psurface(1:nlon,1:nlat,1:ntim),KIND=r8), &
               llon, llat, ltime, status=locerrstat)
          IF ( locerrstat /= 0 ) THEN
             call tell_error (tell_runtime_error, &
                  "omi_climatology: surface pressure interpolation failed", errstat)
             RETURN
          END IF

          DO n = 1, CmETA
             
             local_heights(ixtrack,itimes,n) = (( Ap(n) + local_psurf(ixtrack,itimes) * Bp(n)  ) + &
                  ( Ap(n+1) + local_psurf(ixtrack,itimes) * Bp(n+1) )) / 2.0 * 1.D2
             
             ! Interpolate temperature to lon,lat,hrs
             local_temperature(ixtrack,itimes,n) =  linInterpol(nlon,nlat,ntim, &
                  REAL(lonvals(idx_lon(1):idx_lon(2)),KIND=r8), &
                  REAL(latvals(idx_lat(1):idx_lat(2)),KIND=r8), &
                  REAL(timevals(idx_tim(1):idx_tim(2)),KIND=r8), &
                  REAL(Temperature(1:nlon,1:nlat,n,1:ntim),KIND=r8), &
                  llon, llat, ltime, status=status)
             IF ( locerrstat /= 0 ) THEN
                call tell_error (tell_runtime_error, &
                     "omi_climatology: temperature interpolation failed", errstat)
                RETURN
             END IF

             
             ! Interpolate water vapor profile to lon,lat,hrs
             lh2o(n) =  linInterpol(nlon,nlat,ntim, &
                  REAL(lonvals(idx_lon(1):idx_lon(2)),KIND=r8), &
                  REAL(latvals(idx_lat(1):idx_lat(2)),KIND=r8), &
                  REAL(timevals(idx_tim(1):idx_tim(2)),KIND=r8), &
                  REAL(H2O_profiles(1:nlon,1:nlat,n,1:ntim),KIND=r8), &
                  llon, llat, ltime, status=status)
             IF ( locerrstat /= 0 ) THEN
                call tell_error (tell_runtime_error, &
                     "omi_climatology: water vapor interpolation failed", errstat)
                RETURN
             END IF

             ! Interpolate trace gas profile to lon,lat,hrs
             lgas(n) =  linInterpol(nlon,nlat,ntim, &
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

             rho = 0.0_r8
             aircolumn = 0.0_r8
             lhgt = 0.0_r8
             
             ! Convert input water vapor mixing ratio from PPB to unitless
             thish2omxr = lh2o(n) / 1.0E9
            
             ! Calculate mean molecular weight of wet air 
             Mwet = (1.0_r8 - thish2omxr)*Mdry + thish2omxr*Mh2o
             
             ! Calculate gas constant for wet air
             Rwet = Rstar / Mwet
             
             ! Calculate layer thickness using the following
             ! dz = -(R*T/g) * dlnP
             detlnp = LOG( ( Ap(n+1) + local_psurf(ixtrack,itimes) * Bp(n+1) ) ) - &
                  LOG( ( Ap(n) + local_psurf(ixtrack,itimes) * Bp(n)   ) ) !(in Pa)
             
             lhgt = -Rwet * local_temperature(ixtrack,itimes,n) * detlnp / gplanet ! meter
             rho  = local_heights(ixtrack,itimes,n) / local_temperature(ixtrack,itimes,n) / Rwet ! Kg m-3
             aircolumn = rho*lhgt*Navogadro/m2tocm2 ! # air/cm^2
             
             ! -------------------------------------------------------------
             climatology(ixtrack,itimes,n) = aircolumn * lgas(n) / 1.0E9 ! [GAS]/cm^2
             
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
       gasdatafieldname, h2odatafieldname, errstat)
    ! ==========================================================
    ! This subroutine reads in the climatology from GEOS-Chem or
    ! other source. The climatology file needs to be conform to
    ! the format assumed here.
    ! ==========================================================    
    USE OMSAO_he5_module, ONLY: HE5_SWrdfld, HE5_SWrdlattr
    USE OMSAO_errstat_module, only : pge_errstat_ok
    IMPLICIT NONE

    INTEGER (KIND=i4), INTENT (IN) :: swath_id
    INTEGER (KIND=i4), DIMENSION(2), INTENT (IN) :: idx_lon, idx_lat, idx_tim
    CHARACTER (LEN=MAX_STR_LEN) :: gasdatafieldname, h2odatafieldname

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4) :: he5stat
    INTEGER   (KIND=C_LONG) :: nlon, nlat, nlev, ntim
    CHARACTER (LEN=15), PARAMETER :: cli_Psurf_field       = 'SurfacePressure'
    CHARACTER (LEN=18), PARAMETER :: cli_Temperature_field = 'TemperatureProfile'
    REAL (KIND=r4) :: scale_gas, scale_Psurf, scale_temperature, scale_H2O

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

    ! -------------------------
    ! Read the tables: Psurface
    ! -------------------------
    he5_start_3d  = (/ INT(idx_lon(1)-1,8), INT(idx_lat(1)-1,8), INT(idx_tim(1)-1,8) /)
    he5_stride_3d = (/ onecl, onecl, onecl /)
    he5_edge_3d   = (/  nlon,  nlat,  ntim /)
    he5stat = HE5_SWrdfld ( &
         swath_id, cli_Psurf_field, &
         he5_start_3d, he5_stride_3d, he5_edge_3d, &
         Psurface(1:nlon,1:nlat,1:ntim) )

    ! ----------------------------
    ! Read the tables: Temperature
    ! ----------------------------
    he5_start_4d  = (/ INT(idx_lon(1)-1,8),INT(idx_lat(1)-1,8),zerocl,INT(idx_tim(1)-1,8) /)
    he5_stride_4d = (/ onecl, onecl, onecl, onecl /)
    he5_edge_4d   = (/  nlon,  nlat, nlev, ntim /)
    he5stat = HE5_SWrdfld ( &
      swath_id, cli_Temperature_field, &
      he5_start_4d, he5_stride_4d, he5_edge_4d, &
      Temperature(1:nlon,1:nlat,1:nlev,1:ntim) )

    ! --------------------------------------
    ! Read datafields scale factor attribute
    ! --------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, cli_Psurf_field, "ScaleFactor", scale_Psurf       )
    he5stat = HE5_SWrdlattr ( swath_id, cli_Temperature_field, "ScaleFactor", scale_Temperature )

    IF ( he5stat /= pge_errstat_ok ) then
      call tell_error (tell_io_read_error, "read_climatology: reading climatology data fields", &
                       errstat)
      return
    endif

    ! ----------------------------
    ! Read data from gas datafield
    ! ----------------------------
    he5stat = HE5_SWrdfld ( &
      swath_id, TRIM(ADJUSTL(gasdatafieldname)), &
      he5_start_4d, he5_stride_4d, he5_edge_4d, &
      Gas_profiles(1:nlon,1:nlat,1:nlev,1:ntim) )
    
    ! -----------------------------------------
    ! Read gas datafield scale factor attribute
    ! -----------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(gasdatafieldname)),&
      "ScaleFactor", scale_gas       )

    ! ----------------------------
    ! Read data from H2O datafield
    ! ----------------------------
    he5stat = HE5_SWrdfld (                                &
         swath_id, TRIM(ADJUSTL(h2odatafieldname)),        &
         he5_start_4d, he5_stride_4d, he5_edge_4d,         &
         H2O_profiles(1:nlon,1:nlat,1:nlev,1:ntim) )

    ! -----------------------------------------
    ! Read gas datafield scale factor attribute
    ! -----------------------------------------
    he5stat = HE5_SWrdlattr ( swath_id, TRIM(ADJUSTL(h2odatafieldname)),&
              "ScaleFactor", scale_H2O       )
    
    ! ------------------------------------
    ! Apply scaling factors to data fields
    ! ------------------------------------
    Temperature  = Temperature  * scale_Temperature
    Psurface     = Psurface     * scale_Psurf
    Gas_profiles = Gas_profiles * scale_gas
    H2O_profiles = H2O_profiles * scale_H2O

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

        ! -----------------------------------------
        ! Locate two closest indices to lon and lat
        ! in OMLER_longitude and OMLER_latitudes.
        ! If result out of bounds bring it to the
        ! closest boundary.
        ! ------------------------------------------
        CALL GetNode(REAL(OMLER_longitude,KIND=r8),lonp, &
             lon_idx(1), 'Lower')
        IF (lon_idx(1) .EQ. -2) lon_idx(1) = 1
        IF (lon_idx(1) .EQ. -3) lon_idx(1) = OMLER_n_longitudes
        CALL GetNode(REAL(OMLER_longitude,KIND=r8),lonp, &
             lon_idx(2), 'Upper')
        IF (lon_idx(2) .EQ. -2) lon_idx(2) = 1
        IF (lon_idx(2) .EQ. -3) lon_idx(2) = OMLER_n_longitudes
        nlon = lon_idx(2)-lon_idx(1)+1

        CALL GetNode(REAL(OMLER_latitude,KIND=r8),latp, &
             lat_idx(1), 'Lower')
        IF (lat_idx(1) .EQ. -2) lat_idx(1) = 1
        IF (lat_idx(1) .EQ. -3) lat_idx(1) = OMLER_n_latitudes
        CALL GetNode(REAL(OMLER_latitude,KIND=r8),latp, &
             lat_idx(2), 'Upper')
        IF (lat_idx(2) .EQ. -2) lat_idx(2) = 1
        IF (lat_idx(2) .EQ. -3) lat_idx(2) = OMLER_n_latitudes
        nlat = lat_idx(2)-lat_idx(1)+1

        IF (nlon .EQ. 1) lonp = REAL(OMLER_longitude(lon_idx(1)),KIND=r8)
        IF (nlon .EQ. 1) latp = REAL(OMLER_latitude(lat_idx(1)),KIND=r8)
        
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
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = REAL(Cmlat,KIND=r8)
    he5stat = HE5_SWrdfld ( swath_id, cli_lat_field, &
      he5_start_1d, he5_stride_1d, he5_edge_1d, latvals(1:Cmlat) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = REAL(Cmlon,KIND=r8)
    he5stat = HE5_SWrdfld ( swath_id, cli_lon_field, &
      he5_start_1d, he5_stride_1d, he5_edge_1d, lonvals(1:Cmlon) )
    he5_start_1d = zerocl ; he5_stride_1d = onecl ; he5_edge_1d = REAL(CmHRS,KIND=r8)
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
    if (allocated(Temperature)) DEALLOCATE(Temperature, stat=errstat)
    if (allocated(Gas_profiles)) DEALLOCATE(Gas_profiles, stat=errstat)
    if (allocated(H2O_profiles)) DEALLOCATE(H2O_profiles, stat=errstat)
    if (allocated(Psurface)) DEALLOCATE(Psurface, stat=errstat)
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
    deallocate ( Temperature, &
                 Gas_profiles, H2O_profiles, Psurface, stat=errstat)
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

    ALLOCATE (Temperature (Cmlon, Cmlat, CmETA, CmHRS), &
              Gas_profiles(Cmlon, Cmlat, CmETA, CmHRS), &
              H2O_profiles(Cmlon, Cmlat, CmETA, CmHRS), &
              Psurface(Cmlon,Cmlat,CmHRS), STAT=estat ) ;
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
    if (allocated(vl_OzC0)) then
      deallocate (vl_OzC0, vl_OzC1, vl_OzC2, vl_pre, vl_sza, vl_vza, &
                  vl_wav, vl_toms, vl_air, vl_alt, vl_ozo, vl_tem, &
                  vl_I0, vl_I1, vl_I2, vl_Ir, vl_Sb, &
                  vl_dI0, vl_dI1, vl_dI2, vl_dIr, stat=errstat)
    endif
    if (errstat /= 0) then
      call tell_error(tell_malloc_error, "vlidort_deallocate failed", errstat)
      return
    endif
  END SUBROUTINE

  SUBROUTINE vlidort_allocate (anozo, ancld, ansza, anvza, anwav, analt, errstat)
    IMPLICIT NONE
    INTEGER (KIND=i4), INTENT (IN) :: anozo, ancld, ansza, anvza, anwav, analt
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    INTEGER   (KIND=i4) :: estat

    if (errstat /= 0) return

    ALLOCATE (vl_OzC0(anwav), &
              vl_OzC1(anwav), &
              vl_OzC2(anwav), &
              vl_pre(ancld), &
              vl_sza(ansza), &
              vl_vza(anvza), &
              vl_wav(anwav), &
              vl_toms(anozo), &
              vl_air(anozo,ancld,analt), &
              vl_alt(anozo,ancld,analt), &
              vl_ozo(anozo,ancld,analt), &
              vl_tem(anozo,ancld,analt), &
              !vl_I0(anozo,ancld,ansza,anvza,anwav), &
              !vl_I1(anozo,ancld,ansza,anvza,anwav), &
              !vl_I2(anozo,ancld,ansza,anvza,anwav), &
              !vl_Ir(anozo,ancld,ansza,anvza,anwav), &
              vl_I0(anvza,ansza,ancld,anwav,anozo), &
              vl_I1(anvza,ansza,ancld,anwav,anozo), &
              vl_I2(anvza,ansza,ancld,anwav,anozo), &
              vl_Ir(anvza,ansza,ancld,anwav,anozo), &
              vl_Sb(anozo,ancld,anwav), &
              !vl_dI0(anozo,ancld,ansza,anvza,anwav,analt), &
              !vl_dI1(anozo,ancld,ansza,anvza,anwav,analt), &
              !vl_dI2(anozo,ancld,ansza,anvza,anwav,analt), &
              !vl_dIr(anozo,ancld,ansza,anvza,anwav,analt), &
              vl_dI0(analt,anvza,ansza,ancld,anwav,anozo), &
              vl_dI1(analt,anvza,ansza,ancld,anwav,anozo), &
              vl_dI2(analt,anvza,ansza,ancld,anwav,anozo), &
              vl_dIr(analt,anvza,ansza,ancld,anwav,anozo), STAT=estat )
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

    USE HDF5, ONLY: HID_T, HSIZE_T, SIZE_T, h5dopen_f, h5dget_space_f, &
      h5dread_f, h5sget_simple_extent_dims_f, h5open_f, h5tcopy_f, &
      h5tset_size_f, h5dclose_f, h5fopen_f, h5fclose_f, &
      H5F_ACC_RDONLY_F, H5T_NATIVE_CHARACTER, H5T_NATIVE_REAL
    !USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: hdferr

    INTEGER(HID_T) :: input_file_id                                  ! File identifier
    INTEGER(HID_T) :: OzC0_did, OzC1_did, OzC2_did,                & ! Dataset identifiers
      llp_did, sza_did, toz_did, vza_did, wav_did, &
      I0_did, I1_did, I2_did, Ir_did, Sb_did,      &
      dI0_did, dI1_did, dI2_did, dIr_did, Fac_did, &
      air_did, alt_did, ozo_did, tem_did, dspace,  &
      Tozo_datatype_id

    INTEGER(HSIZE_T), DIMENSION(1) :: hllp_dim, hllp_maxdim, &
      hsza_dim, hsza_maxdim, &
      htoz_dim, htoz_maxdim, &
      hvza_dim, hvza_maxdim, &
      hwav_dim, hwav_maxdim, &
      hfac_dim, hfac_maxdim
    INTEGER(HSIZE_T), DIMENSION(3) :: halt_dim, halt_maxdim, &
      hSb_dim,  hSb_maxdim
    INTEGER(HSIZE_T), DIMENSION(5) :: hI0_dim,  hI0_maxdim, perm5, dims5
    INTEGER(HSIZE_T), DIMENSION(6) :: hdI0_dim, hdI0_maxdim, perm6, dims6

    INTEGER(SIZE_T)                :: size
    LOGICAL, SAVE :: h5inited = .FALSE.

    INTEGER   (KIND=i4) :: estat
    REAL(KIND=r4), DIMENSION(:,:,:,:,:), ALLOCATABLE :: tmpspc_r4d5
    REAL(KIND=r4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: tmpspc_r4d6
    CHARACTER(LEN=MAX_STR_LEN)       :: filename
    CHARACTER(LEN=MAX_STR_LEN), SAVE :: cached_filename = ""

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=26), PARAMETER :: modulename = 'read_vlidort'

    ! ----------------------
    ! Subroutine starts here
    ! ----------------------
    if (errstat /= 0) return
    !errstat = pge_errstat_ok

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
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, Tozo_datatype_id, hdferr)
    CALL h5tset_size_f(Tozo_datatype_id, size, hdferr)

    ! ******************************************
    ! Find out the dimensions of the input file:
    !  # of pressure levels
    !  # of SZA
    !  # of Ozone profiles
    !  # of VZA
    !  # of wavelenghts
    !  # of altitude levels
    ! ******************************************
    ! -------------------
    ! Opening input TABLE
    ! -------------------
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, input_file_id, hdferr)
    IF (hdferr .eq. -1) THEN
      call tell_error (tell_io_open_error, 'opening '//trim(filename), &
                       errstat)
      return
    END IF

    ! --------------------------------------------------------------------------
    ! Open ozone cross sections, grid, intensity, jacobians and profile datasets
    ! --------------------------------------------------------------------------
    CALL h5dopen_f(input_file_id,'/Cross sections/Ozone C0', OzC0_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Cross sections/Ozone C1', OzC1_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Cross sections/Ozone C2', OzC2_did, hdferr)

    CALL h5dopen_f(input_file_id,'/Grid/Lower level pressure', llp_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/SZA', sza_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/TOMS ozone', toz_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/VZA', vza_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Grid/Wavelength', wav_did,hdferr)

    CALL h5dopen_f(input_file_id,'/Intensity/I0', I0_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/I1', I1_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/I2', I2_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Ir', Ir_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Intensity/Sb', Sb_did, hdferr)

    CALL h5dopen_f(input_file_id,'/Jacobians/Factor', Fac_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dI0', dI0_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dI1', dI1_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dI2', dI2_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Jacobians/dIr', dIr_did, hdferr)

    CALL h5dopen_f(input_file_id,'/Profiles/Air profile', air_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Altitude', alt_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Ozone profile', ozo_did, hdferr)
    CALL h5dopen_f(input_file_id,'/Profiles/Temperature profile', tem_did, hdferr)

    ! -----------------------
    ! Find out the dimensions
    ! -----------------------
    CALL h5dget_space_f(llp_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hllp_dim, hllp_maxdim, hdferr)
    CALL h5dget_space_f(sza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hsza_dim, hsza_maxdim, hdferr)
    CALL h5dget_space_f(toz_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, htoz_dim, htoz_maxdim, hdferr)
    CALL h5dget_space_f(vza_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hvza_dim, hvza_maxdim, hdferr)
    CALL h5dget_space_f(wav_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hwav_dim, hwav_maxdim, hdferr)
    CALL h5dget_space_f(I0_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hI0_dim, hI0_maxdim, hdferr)
    CALL h5dget_space_f(Sb_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hSb_dim, hSb_maxdim, hdferr)
    CALL h5dget_space_f(Fac_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hfac_dim, hfac_maxdim, hdferr)
    CALL h5dget_space_f(dI0_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, hdI0_dim, hdI0_maxdim, hdferr)
    CALL h5dget_space_f(alt_did,dspace,hdferr)
    CALL h5sget_simple_extent_dims_f (dspace, halt_dim, halt_maxdim, hdferr)

    ! ---------------------------------------------------------------
    ! Allocate & initialize variables now that we have the dimensions
    ! ---------------------------------------------------------------
    vl_nozo = INT(hdI0_dim(1), KIND=i4)
    vl_ncld = INT(hdI0_dim(2), KIND=i4)
    vl_nsza = INT(hdI0_dim(3), KIND=i4)
    vl_nvza = INT(hdI0_dim(4), KIND=i4)
    vl_nwav = INT(hdI0_dim(5), KIND=i4)
    vl_nalt = INT(hdI0_dim(6), KIND=i4)

    CALL vlidort_allocate (vl_nozo, vl_ncld, vl_nsza, vl_nvza, vl_nwav, vl_nalt, errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "read_vlidort: allocate failed", &
                       errstat)
      return
    endif

    ! ----------------------------------------------------
    ! Read from the h5 file all these small size variables
    ! ----------------------------------------------------
    CALL h5dread_f(OzC0_did, H5T_NATIVE_REAL, vl_OzC0(1:vl_nwav), hwav_dim, hdferr)
    CALL h5dread_f(OzC1_did, H5T_NATIVE_REAL, vl_OzC1(1:vl_nwav), hwav_dim, hdferr)
    CALL h5dread_f(OzC2_did, H5T_NATIVE_REAL, vl_OzC2(1:vl_nwav), hwav_dim, hdferr)

    CALL h5dread_f(llp_did, H5T_NATIVE_REAL,  vl_pre(1:vl_ncld),  hllp_dim, hdferr)
    CALL h5dread_f(sza_did, H5T_NATIVE_REAL,  vl_sza(1:vl_nsza),  hsza_dim, hdferr)
    CALL h5dread_f(toz_did, Tozo_datatype_id, vl_toms(1:vl_nozo), htoz_dim, hdferr)
    CALL h5dread_f(vza_did, H5T_NATIVE_REAL,  vl_vza(1:vl_nvza),  hvza_dim, hdferr)
    CALL h5dread_f(wav_did, H5T_NATIVE_REAL,  vl_wav(1:vl_nwav),  hwav_dim, hdferr)

    ALLOCATE (tmpspc_r4d5(vl_nozo,vl_ncld,vl_nsza,vl_nvza,vl_nwav), STAT=estat)
    if (estat /= 0) then
      call tell_error (tell_malloc_error, "read_vlidort: allocate failed", &
                       errstat)
      return
    endif

    perm5 = (/ 5, 3, 2, 1, 4 /)
    dims5 = (/ vl_nvza, vl_nsza, vl_ncld, vl_nwav, vl_nozo /)
    CALL h5dread_f(I0_did, H5T_NATIVE_REAL, tmpspc_r4d5, hI0_dim, hdferr)
    vl_I0 = reshape (tmpspc_r4d5, dims5, order=perm5)
    CALL h5dread_f(I1_did, H5T_NATIVE_REAL, tmpspc_r4d5, hI0_dim, hdferr)
    vl_I1 = reshape (tmpspc_r4d5, dims5, order=perm5)
    CALL h5dread_f(I2_did, H5T_NATIVE_REAL, tmpspc_r4d5, hI0_dim, hdferr)
    vl_I2 = reshape (tmpspc_r4d5, dims5, order=perm5)
    CALL h5dread_f(Ir_did, H5T_NATIVE_REAL, tmpspc_r4d5, hI0_dim, hdferr)
    vl_Ir = reshape (tmpspc_r4d5, dims5, order=perm5)
    where (vl_Ir /= vl_Ir)    ! isnan is preferable, but non-standard
      vl_Ir = 0.0
    end where

    DEALLOCATE(tmpspc_r4d5, stat=estat)
    if (estat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_vlidort: deallocate tmpspc_r4d5 failed", errstat)
      return
    endif

    CALL h5dread_f(Sb_did, H5T_NATIVE_REAL,  &
      vl_Sb(1:vl_nozo,1:vl_ncld,1:vl_nwav),  &
      hSb_dim, hdferr)

    CALL h5dread_f(Fac_did, H5T_NATIVE_REAL,  vl_Factor, hfac_dim, hdferr)

    ALLOCATE (tmpspc_r4d6(vl_nozo,vl_ncld,vl_nsza,vl_nvza,vl_nwav,vl_nalt), STAT=estat)
    if (estat /= 0) then
      call tell_error (tell_malloc_error, "read_vlidort: allocate failed", &
                       errstat)
      return
    endif

    perm6 = (/ 6, 4, 3, 2, 5, 1 /)
    dims6 = (/ vl_nalt, vl_nvza, vl_nsza, vl_ncld, vl_nwav, vl_nozo /)
    CALL h5dread_f(dI0_did, H5T_NATIVE_REAL, tmpspc_r4d6, hdI0_dim, hdferr)
    vl_dI0 = reshape (tmpspc_r4d6, dims6, order=perm6)
    CALL h5dread_f(dI1_did, H5T_NATIVE_REAL, tmpspc_r4d6, hdI0_dim, hdferr)
    vl_dI1 = reshape (tmpspc_r4d6, dims6, order=perm6)
    CALL h5dread_f(dI2_did, H5T_NATIVE_REAL, tmpspc_r4d6, hdI0_dim, hdferr)
    vl_dI2 = reshape (tmpspc_r4d6, dims6, order=perm6)
    CALL h5dread_f(dIr_did, H5T_NATIVE_REAL, tmpspc_r4d6, hdI0_dim, hdferr)
    vl_dIr = reshape (tmpspc_r4d6, dims6, order=perm6)

    DEALLOCATE (tmpspc_r4d6, stat=estat)
    if (estat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_vlidort: deallocate tmpspc_r4d6 failed", errstat)
      return
    endif

    CALL h5dread_f(air_did, H5T_NATIVE_REAL,  vl_air(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)
    CALL h5dread_f(alt_did, H5T_NATIVE_REAL,  vl_alt(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)
    CALL h5dread_f(ozo_did, H5T_NATIVE_REAL,  vl_ozo(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)
    CALL h5dread_f(tem_did, H5T_NATIVE_REAL,  vl_tem(1:vl_nozo,1:vl_ncld,1:vl_nalt),  halt_dim, hdferr)

    ! --------------
    ! Close datasets
    ! --------------
    CALL h5dclose_f (OzC0_did, hdferr)
    CALL h5dclose_f (OzC1_did, hdferr)
    CALL h5dclose_f (OzC2_did, hdferr)

    CALL h5dclose_f (llp_did, hdferr)
    CALL h5dclose_f (sza_did, hdferr)
    CALL h5dclose_f (toz_did, hdferr)
    CALL h5dclose_f (vza_did, hdferr)
    CALL h5dclose_f (wav_did, hdferr)

    CALL h5dclose_f(I0_did, hdferr)
    CALL h5dclose_f(I1_did, hdferr)
    CALL h5dclose_f(I2_did, hdferr)
    CALL h5dclose_f(Ir_did, hdferr)
    CALL h5dclose_f(Sb_did, hdferr)

    CALL h5dclose_f(Fac_did, hdferr)
    CALL h5dclose_f(dI0_did, hdferr)
    CALL h5dclose_f(dI1_did, hdferr)
    CALL h5dclose_f(dI2_did, hdferr)
    CALL h5dclose_f(dIr_did, hdferr)

    CALL h5dclose_f (air_did, hdferr)
    CALL h5dclose_f (alt_did, hdferr)
    CALL h5dclose_f (ozo_did, hdferr)
    CALL h5dclose_f (tem_did, hdferr)

    ! ----------
    ! Close file
    ! ----------
    CALL h5fclose_f(input_file_id, hdferr)

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
      !CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
      !  modulename//f_sep//'OMIL2 access failed.', vb_lev_default, errstat )
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
      !CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
      !modulename//f_sep//'OMIL2 CFR access failed.', vb_lev_default, errstat )
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
      !CALL error_check ( locerrstat, OMI_S_SUCCESS, pge_errstat_error, OMSAO_E_PREFITCOL, &
      !modulename//f_sep//'OMIL2 CTP access failed.', vb_lev_default, errstat )
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

  SUBROUTINE amf_diagnostic ( nt, nx, & !lat, lon,
                             sza, vza, snow, glint, xtrange, &
                             l2cfr, l2ctp, amfdiag )

    USE OMSAO_parameters_module, only: i2_missval
    USE OMSAO_omidata_module,   ONLY: omi_oobview_amf, omi_glint_add, omi_bigsza_amf
    !USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    !REAL    (KIND=r4),                          INTENT (IN) :: ctpmin, ctpmax
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: sza, vza
    !REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: lat, lon
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
    INTEGER (KIND=i4) :: it, spix, epix

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
    IF ((amf_wvl  .LT. MINVAL(vl_wav) ) .OR. &
        (amf_wvl2 .GT. MAXVAL(vl_wav) ) ) THEN
      RETURN
    END IF

    DO it = 0, nt-1
      spix = xtrange(it,1) ; epix = xtrange(it,2)

      ! ----------------
      ! Missing SZA, VZA
      ! ----------------
      WHERE (                                   &
          sza(spix:epix,it) <= r8_missval .OR. &
          vza(spix:epix,it) <= r8_missval        )
        amfdiag(spix:epix,it) = i2_missval
      END WHERE

      ! ----------------------------------------
      ! Out-of-Bound SZA, VZA (but not missing!)
      ! ----------------------------------------
      WHERE ( &
          ( sza(spix:epix,it) < MINVAL(vl_sza) ) .OR. &
          ( sza(spix:epix,it) > MAXVAL(vl_sza) ) .OR. &
          ( vza(spix:epix,it) < MINVAL(vl_vza) ) .OR. &
          ( vza(spix:epix,it) > MAXVAL(vl_vza) )      )
        amfdiag(spix:epix,it) = omi_oobview_amf
      END WHERE

      ! -------------------------------------------------
      ! Out of bounds clouds (to high over land), make it
      ! the highest possible value in the look up table.
      ! Ask Xiong...
      ! -------------------------------------------------
      WHERE ( &
          ( l2ctp(spix:epix,it) < MINVAL(vl_pre(1:vl_ncld))*1013.0_r8 ) )
        l2ctp(spix:epix,it) = MINVAL(vl_pre(1:vl_ncld))*1013.0_r8
      END WHERE

      ! ------------------------------------------------------
      ! For pixel without cloud information set amf to missval
      ! and flag to missval
      ! ------------------------------------------------------
      WHERE ( &
          ( l2cfr(spix:epix,it) .EQ. r8_missval ) .OR. &
          ( l2ctp(spix:epix,it) .EQ. r8_missval)       )
        amfdiag(spix:epix,it) = omi_oobview_amf
      END WHERE

      ! ------------------------------------------------------
      ! And AMFDIAG values > OOB must be good and are set to 0
      ! if we have "good" clouds
      ! ------------------------------------------------------
      WHERE ( &
          ( amfdiag(spix:epix,it) > omi_oobview_amf ) .AND. &
          (l2cfr(spix:epix,it) >= 0.0_r8            ) .AND. &
          (l2ctp(spix:epix,it) >= 0.0_r8            )       )
        amfdiag(spix:epix,it) = 0_i2
      END WHERE

      ! --------------------------------------------------
      ! Angles above the top value set on the control file
      ! are calculated "using this maximum value".
      ! --------------------------------------------------
      WHERE ( &
          ( sza(spix:epix,it)     .GE. amf_max_sza       ) .AND. &
          ( amfdiag(spix:epix,it) .GT. omi_oobview_amf ) )
        amfdiag(spix:epix,it) = omi_bigsza_amf + amfdiag(spix:epix,it)
      END WHERE

      ! -----------------------
      ! Start with the ice flag
      ! -----------------------
      WHERE (                                       &
          amfdiag     (spix:epix,it) >= 0_i2 .AND. &
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

    END DO

    RETURN
  END SUBROUTINE amf_diagnostic

  SUBROUTINE compute_scatt ( nt, nx, albedo, sza, vza, l2ctp, l2cfr, terrain_height, cli_heights, amfdiag, &
      scattw)

    USE OMSAO_linterpolation_module, ONLY: lininterpol
    !USE OMSAO_variables_module, ONLY: verb_thresh_lev
    USE ezspline_interpolation, ONLY: ezspline_2d_interpolation
    !USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                INTENT (IN) :: nt, nx
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1),       INTENT (IN) :: amfdiag
    REAL    (KIND=r4), DIMENSION (1:nx,0:nt-1),       INTENT (IN) :: sza, vza, terrain_height
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1),       INTENT (IN) :: albedo, l2ctp, l2cfr
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT (IN) :: cli_heights
    ! ------------------
    ! Modified variables
    ! ------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT (INOUT) :: scattw

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: itime, ixtrack, ispre, isozo, isalt, iswav, issza, isvza, status, one
    INTEGER (KIND=i4), DIMENSION(1) :: iwavs, iwavf, index_thg, index_cld
    REAL    (KIND=r8) :: temp, tempsquare, grad
    REAL    (KIND=r8) :: ozo_abs, Intensity, Jacobian, Oz_xs, Intensity_cld, Jacobian_cld
    REAL    (KIND=r8) :: crf, nwavs!, Icr, Icl ,cloud_scattw, clear_scattw
    REAL    (KIND=r8), DIMENSION(vl_ncld, vl_nsza, vl_nvza, vl_nalt) :: scattwe, scattwe_cld
    REAL    (KIND=r8), DIMENSION(vl_ncld, vl_nsza, vl_nvza)          :: Inte_clear, Inte_cloud
    REAL    (KIND=r8), DIMENSION(vl_nalt) :: re_alt
    REAL    (KIND=r8), DIMENSION(vl_ncld) :: re_pre
    REAL    (KIND=r8), DIMENSION(vl_nsza) :: re_sza
    REAL    (KIND=r8), DIMENSION(vl_nvza) :: re_vza
    REAL    (KIND=r8)                     :: local_alb, local_sza, local_vza, &
      local_thg, local_cld, local_cfr, vl_Isum
    REAL    (KIND=r8), DIMENSION(1)       :: ezlocal_sza, ezlocal_vza
    REAL    (KIND=r8), DIMENSION(1,1,1)   :: cloud_scattw, clear_scattw
    REAL    (KIND=r8), DIMENSION(1,1)     :: Icr, Icl
    REAL    (KIND=r8), DIMENSION(CmETA)   :: local_chg
    REAL    (kind=8),  PARAMETER :: d2r = 3.141592653589793d0/180.0  !! JED fix
    character (len=72) :: logmsg
    ! -----------------------------------
    ! Find look up table wavelength index
    ! No interpolation, closest available
    ! is selected.
    ! -----------------------------------
    one = 1_i4
    iwavs = MINLOC(ABS(vl_wav - REAL(amf_wvl,  KIND = r4) ))
    iwavf = MINLOC(ABS(vl_wav - REAL(amf_wvl2, KIND = r4) ))
    nwavs = REAL(iwavf(1), KIND=r8) - REAL(iwavs(1), KIND=r8) + 1.0_r8

    ! --------------------------------------------------------
    ! Re-order the pressure, altitude, sza, and vza dimensions
    ! from (look up table). Needed for interpolation.
    ! This should be moved to the program that creates the
    ! look up tables
    ! ------------------------------- ------------------------
    DO issza = 1, vl_nsza
      !re_sza(issza) = cosd(REAL(vl_sza(vl_nsza+1-issza), KIND = r8))
      re_sza(issza) = cos(d2r*REAL(vl_sza(vl_nsza+1-issza), KIND = r8))  ! JED fix
    END DO
    DO isvza = 1, vl_nvza
      !re_vza(isvza) = cosd(REAL(vl_vza(vl_nvza+1-isvza), KIND = r8))
      re_vza(isvza) = cos(d2r*REAL(vl_vza(vl_nvza+1-isvza), KIND = r8)) ! JED fix
    END DO
    DO ispre = 1, vl_ncld
      re_pre(ispre) = REAL(vl_pre(ispre), KIND = r8) * 1013.0_r8
    END DO
    DO isalt = 1, vl_nalt
      !!$       re_alt(isalt) = 1013.0_r8 * (10.0_r8 ** ( REAL(vl_alt(1,1,vl_nalt+1-isalt), KIND = r8) / (-16.0_r8)))
      re_alt(isalt) = 1013.0_r8 * (10.0_r8 ** ( REAL(vl_alt(1,1,isalt), KIND = r8) / (-16.0_r8)))
    END DO

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

        ! ----------------------------------------------
        ! If sza > amf_max_sza set it for calculation to
        ! amf_max_sza
        ! ----------------------------------------------
        local_sza = REAL(sza(ixtrack,itime), KIND = r8)
        IF (local_sza .GT. amf_max_sza) local_sza = amf_max_sza

        ! ---------------------
        ! Albedo for this pixel
        ! ---------------------
        local_alb    = albedo(ixtrack,itime)
        local_cld    = l2ctp(ixtrack,itime)
        local_cfr    = l2cfr(ixtrack,itime)
        !local_sza    = cosd(REAL(sza(ixtrack,itime), KIND = r8))
        local_sza    = cos(d2r*REAL(sza(ixtrack,itime), KIND = r8))  ! JED fix
        !local_vza    = cosd(REAL(vza(ixtrack,itime), KIND = r8))
        local_vza    = cos(d2r*REAL(vza(ixtrack,itime), KIND = r8))  ! JED fix
        local_thg    = REAL(terrain_height(ixtrack,itime), KIND = r8)
        local_chg(:) = cli_heights(ixtrack,itime,:)

        ! ----------------------------------------------
        ! Convert pixel terrain height to pressure using
        ! Xiong suggested to use pressure altitude:
        !  Z = -16 alog10 (P / Po) Z in km and P in hPa.
        ! ----------------------------------------------
        local_thg = 1013.0_r8 * (10.0_r8 ** (local_thg / 1000.0_r8 / (-16.0_r8)))

        !Bringing it to the lowest available pressure if needed
        IF (local_thg .GT. 1013.0) local_thg = 1013.0_r8
        !Bringing clouds heights to lowest available pressure if needed. Weird yes, but just in case
        IF (local_cld .GT. 1013.0) local_cld = 1013.0_r8

        ! -----------------------------------------------
        ! Find the two closest cloud top levels from the
        ! scattering weights table suitable for local_thg
        ! and local_cld
        ! -----------------------------------------------
        index_thg = MINLOC(ABS(re_pre - local_thg))
        ! Be sure that we are below the level of the surface in the table
        IF ( (local_thg .GT. re_pre(index_thg(1))) .AND. (index_thg(1) .GT. 1) ) index_thg = index_thg-1
        index_cld = MINLOC(ABS(re_pre - local_cld))
        ! Be sure that we are below the level of the cloud in the table
        IF ( (local_cld .GT. re_pre(index_cld(1))) .AND. (index_cld(1) .GT. 1) ) index_cld = index_cld-1

        ! ----------------------------------------------------------
        ! First compute back from the parametrization the scattering
        ! weights for the given wavelength and albedo.
        ! ----------------------------------------------------------
        DO isozo = 1, 1
          DO iswav = iwavs(1), iwavf(1) !vl_nwav

            ! ---------------------------
            ! Initialize scattwe to zeros
            ! ---------------------------
            scattwe      = 0.0_r8
            scattwe_cld  = 0.0_r8

            DO ispre = 1, vl_ncld
              DO issza = 1, vl_nsza
                DO isvza = 1, vl_nvza

                  vl_Isum = REAL(vl_I0(isvza,issza,ispre,iswav,isozo), KIND = r8) &
                    + REAL(vl_I1(isvza,issza,ispre,iswav,isozo), KIND = r8) &
                    + REAL(vl_I2(isvza,issza,ispre,iswav,isozo), KIND = r8)

                  Intensity = vl_Isum &
                    + ((local_alb * REAL(vl_Ir(isvza,issza,ispre,iswav,isozo), KIND = r8)) &
                       / (1.0_r8 - local_alb * REAL(vl_Sb(isozo,ispre,iswav), KIND = r8)))

                  Intensity_cld = vl_Isum &
                    + ((amf_alb_cld * REAL(vl_Ir(isvza,issza,ispre,iswav,isozo), KIND = r8)) &
                       / (1.0_r8 - (amf_alb_cld * REAL(vl_Sb(isozo,ispre,iswav), KIND = r8))))

                  ! -----------------------------------------
                  ! vl_ncld+1-ispre & vl_nalt+1-isalt to have
                  ! ascending orden for the interpolation
                  ! ----------------------------------------------------
                  ! Intensities for the calculation of the cloudy pixels
                  ! ----------------------------------------------------
                  Inte_clear(ispre, vl_nsza+1-issza, vl_nvza+1-isvza) = Intensity
                  Inte_cloud(ispre, vl_nsza+1-issza, vl_nvza+1-isvza) = Intensity_cld

                  DO isalt = 1, vl_nalt

                    Temp       = REAL(vl_tem(isozo,ispre,isalt), KIND = r8) - 273.15_r8
                    TempSquare = Temp * Temp

                    Oz_xs = REAL(vl_OzC0(iswav), KIND = r8) + &
                      REAL(vl_OzC1(iswav), KIND = r8) * Temp + &
                      REAL(vl_OzC2(iswav), KIND = r8) * TempSquare

                    ozo_abs = Oz_xs * REAL(vl_ozo(isozo,ispre,isalt), KIND = r8)

                    IF (ozo_abs.eq.0.0_r8) CYCLE

                    IF (intensity.ne.0.0_r8) THEN
                      Jacobian = &
                        ( REAL(vl_dI0(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8) &
                         + REAL(vl_dI1(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8) &
                         + REAL(vl_dI2(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8) &
                         + ((local_alb * REAL(vl_dIr(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8)) &
                            / (1.0_r8 - (local_alb * REAL(vl_Sb(isozo,ispre,iswav), KIND = r8)))) &
                        ) / REAL(vl_Factor, KIND = r8)

                      IF (Jacobian.ne.0_r8) THEN
                        scattwe(ispre, vl_nsza+1-issza, vl_nvza+1-isvza, isalt) &
                          = -Jacobian / Intensity / ozo_abs
                      ENDIF
                    ENDIF

                    IF (intensity_cld.ne.0.0_r8) THEN
                      Jacobian_cld = &
                        ( REAL(vl_dI0(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8) &
                         + REAL(vl_dI1(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8) &
                         + REAL(vl_dI2(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8) &
                         + ((amf_alb_cld &
                             * REAL(vl_dIr(isalt,isvza,issza,ispre,iswav,isozo), KIND = r8) &
                            ) / (1.0_r8 - (amf_alb_cld &
                                           * REAL(vl_Sb(isozo,ispre,iswav), KIND = r8)))) &
                        ) / REAL(vl_Factor, KIND = r8)

                      IF (Jacobian_cld .NE. 0.0_r8) THEN
                        scattwe_cld(ispre, vl_nsza+1-issza, vl_nvza+1-isvza, isalt) &
                          = -Jacobian_cld / Intensity_cld / ozo_abs
                      ENDIF
                    ENDIF
                  END DO ! End altitudes (scattering look up tables)
                END DO ! End vza loop
              END DO ! End sza loop
            END DO ! End ispre (pressure) loop

            ! ---------------------------------------
            ! Working out the cloud radiance fraction
            ! See note below, Boersma et al. 2011 and
            ! Martin et al. 2003 (crf used below)
            ! ---------------------------------------
            ezlocal_sza = local_sza
            ezlocal_vza = local_vza

            Icr = 0.0_r8
            !!$                Icr = linInterpol (                       &
            !!$                           vl_ncld,   vl_nsza,   vl_nvza, &
            !!$                            re_pre,    re_sza,    re_vza, &
            !!$                        Inte_clear,                       &
            !!$                         local_thg, local_sza, local_vza, &
            !!$                           status=status)
            CALL ezspline_2d_interpolation (vl_nsza,vl_nvza,re_sza,re_vza,Inte_clear(index_thg(1),:,:), &
              one,one,ezlocal_sza(1:1),ezlocal_vza(1:1),Icr(one,one), &
              status)
            Icl = 0.0_r8
            !!$                Icl = linInterpol (                       &
            !!$                           vl_ncld,   vl_nsza,   vl_nvza, &
            !!$                            re_pre,    re_sza,    re_vza, &
            !!$                        Inte_cloud,                       &
            !!$                         local_cld, local_sza, local_vza, &
            !!$                           status=status)
            CALL ezspline_2d_interpolation (vl_nsza,vl_nvza,re_sza,re_vza,Inte_cloud(index_cld(1),:,:), &
              one,one,ezlocal_sza(1:1),ezlocal_vza(1:1),Icl(one,one), &
              status)
            crf = 0.0_r8
            crf = local_cfr * Icl(one,one) / &
              (local_cfr * Icl(one,one) + (1 - local_cfr) * Icr(one,one) )

            ! ----------------------------------
            ! Interpolate for each altitude
            ! to the given sza, vza and pressure
            ! ----------------------------------
            DO isalt = 1, CmETA ! Loop over altitudes, climatology

              cloud_scattw = 0.0_r8
              clear_scattw = 0.0_r8
              ! --------------------------------------------------
              ! If we are below the level of the clouds, the cloud
              ! scattering weights are not needed. Cloud fraction
              ! must be GT than 0.0. Only interpolate for values
              ! above the cloud top.
              ! --------------------------------------------------
              IF ( local_chg(isalt) .LE. local_cld .AND. local_cfr .GT. 0.0 ) THEN
                !!$                      cloud_scattw =  linInterpol (                           &
                !!$                           vl_ncld, vl_nsza, vl_nvza, vl_nalt,                &
                !!$                           re_pre,  re_sza,  re_vza,  re_alt,                 &
                !!$                           scattwe_cld,                                       &
                !!$                           local_cld, local_sza, local_vza, local_chg(isalt), &
                !!$                           status=status)
                cloud_scattw(one,one,one) =  linInterpol (              &
                  vl_nsza, vl_nvza, vl_nalt,                         &
                  re_sza,  re_vza,  re_alt,                          &
                  scattwe_cld(index_cld(1),:,:,:),                   &
                  local_sza, local_vza, local_chg(isalt),            &
                  status=status)
                !!$                      CALL ezspline_3d_interpolation (vl_nsza,vl_nvza,vl_nalt,re_sza,re_vza,re_alt,scattwe_cld(index_cld(1),:,:,:), &
                !!$                           one,one,one,ezlocal_sza(one),ezlocal_vza(one),local_chg(isalt),cloud_scattw(one,one,one), &
                !!$                           status)
              END IF

              ! --------------------------------------------------
              ! If we are below the level of the land, no need to
              ! work out those scattering weights.
              ! --------------------------------------------------
              IF ( local_chg(isalt) .LE. local_thg ) THEN
                !!$                      clear_scattw =  linInterpol (                           &
                !!$                           vl_ncld, vl_nsza, vl_nvza, vl_nalt,                &
                !!$                           re_pre,  re_sza,  re_vza,  re_alt,                 &
                !!$                           scattwe,                                           &
                !!$                           local_thg, local_sza, local_vza, local_chg(isalt), &
                !!$                           status=status)
                clear_scattw(one,one,one) = linInterpol (           &
                  vl_nsza, vl_nvza, vl_nalt,                     &
                  re_sza,  re_vza,  re_alt,                     &
                  scattwe(index_thg(1),:,:,:),                   &
                  local_sza, local_vza, local_chg(isalt),        &
                  status=status)
                !!$                      CALL ezspline_3d_interpolation (vl_nsza,vl_nvza,vl_nalt,re_sza,re_vza,re_alt,scattwe(index_thg(1),:,:,:), &
                !!$                           one,one,one,ezlocal_sza(one),ezlocal_vza(one),local_chg(isalt),clear_scattw(one,one,one), &
                !!$                           status)
              END IF

              ! ---------------------------------------------------------------------------------
              ! Boersma et al. 2011 AMT, 4, 2011
              ! Cloud radiance fraction: Crf= Cfr * Icl / Ir
              !  We define Icl = From the Vlidort calculation, see above
              !            Icr = From the Vlidort calculation, see above
              !             Ir = Cfr * Icl + (1 - Cfr) * Icr (Total pixel radiance)
              !
              ! Now the scattering weights become w = crf * scatt_cloud + (1 - crf) * scatt_clear
              !  We add the scattweights calculated in the previous wavelengths.
              ! ---------------------------------------------------------------------------------
              scattw(ixtrack,itime,isalt) = (crf * cloud_scattw(one,one,one) + (1.0_r8 - crf) * clear_scattw(one,one,one)) &
                + scattw(ixtrack,itime,isalt)

            END DO ! End looop over altitudes
            ! -------------------------------------------------------------------------------------------------
            ! If there is a negative value on the lower layer of the scattw it has to do with the interpolation
            ! Quick fix, work out a gradient from two layers above and apply it to this first layer assuming
            ! linear behaviour. Next version of the lookup tables should have at least two levels below surface
            ! level to prevent the need for this. Even worst the 0.7 scaling for extreme cases
            ! -------------------------------------------------------------------------------------------------
            IF (scattw(ixtrack, itime, 1) .LT. 0.0) THEN
              grad = (scattw(ixtrack, itime,3) - scattw(ixtrack, itime, 2)) / log(local_chg(3)/local_chg(2))
              scattw(ixtrack, itime, 1) = scattw(ixtrack, itime, 2) - grad * log(local_chg(2)/local_chg(1))
              IF (scattw(ixtrack,itime,1) .LT. 0.0) scattw(ixtrack,itime,1) = scattw(ixtrack,itime,2) * 0.7
            END IF

          END DO ! End wavelength loop

        END DO ! End ozone profile loop
        scattw(ixtrack,itime,1:CmETA) = scattw(ixtrack,itime,1:CmETA) / nwavs

        !  Set non-physical entries to zero.
        WHERE ( scattw(ixtrack,itime,1:CmETA) < 0.0_r8 )
          scattw(ixtrack,itime,1:CmETA) = 0.0_r8
        END WHERE

      END DO ! End loop xtrack
      write(logmsg, '(a,1x,i5)')'Scattering weights line', itime
      call tell_log (1, logmsg)
      !IF ( verb_thresh_lev .GE. vb_lev_screen ) WRITE(*,*) 'Scattering weights line', itime

    END DO ! End loop lines

  END SUBROUTINE COMPUTE_SCATT

  SUBROUTINE compute_amf ( nt, nx, CmETA, climatology, &
      scattw, saoamf, amfdiag, errstat)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                INTENT(IN) :: nt, nx, CmETA
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1,CmETA), INTENT(IN) :: climatology, scattw
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1),       INTENT(IN) :: amfdiag

    ! -----------------------------
    ! Output and modified variables
    ! -----------------------------
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1),       INTENT (INOUT) :: saoamf
    INTEGER (KIND=i4),                                INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                      :: ixtrack, itimes

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
    !USE OMSAO_omidata_module,   ONLY: n_roff_dig
    !USE sao_pge_utils, ONLY: roundoff_2darr_r4, roundoff_2darr_r8
    use datafields, only: albedo_field
    use OMSAO_he5_module, ONLY: HE5_SWWRFLD
    !USE OMSAO_errstat_module
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

  SUBROUTINE write_climatology_he5(climatology, cli_heights, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the Target Gas Profiles from the GEOS-Chem
    ! climatology to the output file.
    ! ===============================================================
    !USE OMSAO_omidata_module,   ONLY: n_roff_dig
    !USE sao_pge_utils, ONLY: roundoff_3darr_r8
    use datafields, only: clialtgrid_field, gasprofile_field
    use OMSAO_he5_module, ONLY: HE5_SWWRFLD
    !USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                              INTENT (IN) :: nt, nx, nl
    REAL    (KIND=r8), DIMENSION(1:nx,0:nt-1,1:nl), INTENT (IN) :: climatology, cli_heights

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

    colloc = cli_heights
    !CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                              TRIM(ADJUSTL(clialtgrid_field)),         &
                              he5_start_3d, he5_stride_3d, he5_edge_3d,&
                              colloc(1:nx,0:nt-1,1:nl) )
    errstat = MAX ( errstat, locerrstat )

  END SUBROUTINE write_climatology_he5

  SUBROUTINE write_scatt_he5(scattw, nt, nx, nl, errstat)

    ! ===============================================================
    ! This routines writes the scattering weigths to the output file.
    ! ===============================================================
    !USE OMSAO_omidata_module,   ONLY: n_roff_dig
    !USE sao_pge_utils, ONLY: roundoff_3darr_r8
    use datafields, only: scaweights_field
    use OMSAO_he5_module, ONLY: HE5_SWWRFLD
    !USE OMSAO_errstat_module

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
    !CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
                              TRIM(ADJUSTL(scaweights_field)),         &
                              he5_start_3d, he5_stride_3d, he5_edge_3d,&
                              colloc(1:nx,0:nt-1,1:nl) )

    !!$    colloc = akernels
    !!$    CALL roundoff_3darr_r8 ( n_roff_dig, nx, nt, nl, colloc(1:nx,0:nt-1,1:nl) )
    !!$    locerrstat = HE5_SWWRFLD ( pge_swath_id,                            &
    !!$                               TRIM(ADJUSTL(avekernels_field)),         &
    !!$                               he5_start_3d, he5_stride_3d, he5_edge_3d,&
    !!$                               colloc(1:nx,0:nt-1,1:nl) )

    errstat = MAX ( errstat, locerrstat )

  END SUBROUTINE write_scatt_he5

  SUBROUTINE he5_amf_write ( &
      pge_idx, nx, nt, saocol, saodco, amfmol, amfgeo, amfdiag, &
      amfcfr, amfctp, errstat )

    USE OMSAO_precision_module, ONLY: i2, i4, r8
    USE OMSAO_he5_module, ONLY: HE5_SWWRFLD, he5_start_2d, he5_stride_2d, &
      he5_edge_2d
    USE OMSAO_errstat_module, only : pge_errstat_ok
    !USE OMSAO_omidata_module,   ONLY: n_roff_dig
    USE OMSAO_indices_module,   ONLY: pge_hcho_idx, pge_gly_idx
    !USE sao_pge_utils, ONLY: roundoff_2darr_r4, roundoff_2darr_r8
    use datafields, only: amfcfr_field, amfctp_field, amfdiag_field, &
      amfgeo_field, amfmol_field, col_field, dcol_field

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    !CHARACTER (LEN=13), PARAMETER :: modulename = 'he5_write_amf'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nx, nt
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: saocol, saodco
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfmol, amfgeo
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: amfcfr, amfctp
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
    !CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfgeo_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    ! -----------------
    ! (3) Molecular AMF
    ! -----------------
    amfloc = REAL ( amfmol, KIND=r4 )
    !CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfmol_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    ! ----------------------------------------------------------
    ! (4) OMHCHO, OMCHOCHO only: AMF cloud fraction and pressure
    ! ----------------------------------------------------------
    IF ( pge_idx == pge_hcho_idx .OR. pge_idx == pge_gly_idx ) THEN
      amfloc = REAL ( amfcfr, KIND=r4 )
      !CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
      locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfcfr_field)), &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
      errstat = MAX ( errstat, locerrstat )

      amfloc = REAL ( amfctp, KIND=r4 )
      !CALL roundoff_2darr_r4 ( n_roff_dig, nx, nt, amfloc(1:nx,0:nt-1) )
      locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(amfctp_field)), &
                                he5_start_2d, he5_stride_2d, he5_edge_2d, amfloc(1:nx,0:nt-1) )
      errstat = MAX ( errstat, locerrstat )
    END IF

    ! -----------------------------------------------------------------------
    ! (5) All PGEs: Output of columns and column uncertainties. For some PGEs
    !     (e.g., OMBRO, OMHCHO, OMCHOCHO) those have been adjusted by the AMF,
    !     but we have as yet to perform the rounding for any of them.
    ! -----------------------------------------------------------------------
    colloc = saocol
    !CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(col_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

    colloc = saodco
    !CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
    locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(dcol_field)), &
                              he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
    errstat = MAX ( errstat, locerrstat )

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

END MODULE OMSAO_wfamf_module

