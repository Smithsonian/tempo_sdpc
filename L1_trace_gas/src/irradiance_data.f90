module irradiance_data

  use OMSAO_precision_module, only: i2, i4, r4, r8
  use tell_module
  USE sao_pge_utils, ONLY: print_array
  implicit none

  type, public :: Irradiance_Data_Type
    integer (kind=i4) :: max_nwave, nxtrack
    ! Each 2d array A has size A[max_nwave, nxtrack], where only
    ! A[1:nwaves[j], j], j=1:nxtrack are meaningful.
    ! The 1d arrays have size [nxtrack].
    integer(kind=i2), dimension(:,:), allocatable :: qflags
    real (kind=r8), dimension(:,:), allocatable :: &
      wavelengths, spectrum
    integer (kind=i4), dimension(:), allocatable :: nwaves
    real (kind=r8), dimension(:), allocatable :: avg_wavelengths
    integer (kind=i4), dimension(:,:), allocatable :: ccdpix_selection
    integer (kind=i4), dimension(:,:), allocatable :: ccdpix_exclusion
  end type

  type (Irradiance_Data_Type) :: Irr_Data

  private
  public irradiance_data_init, Irr_Data

contains

  subroutine irradiance_data_init (rpt_rad, errstat)

    use OMSAO_variables_module, only: Radiance_Paras_Type
    use ctrlvars, only: yn_solar_comp, yn_solmonthave
    implicit none
    type(Radiance_Paras_Type), intent(in) :: rpt_rad
    integer, intent(inout) :: errstat

    if (errstat < 0) return

    ! --------------------------------------------------------------------
    ! Solar Irradiance Processing: If we don't do a solar composite, we can
    ! use a solar monthly average, if not we have to read the irradiance
    ! data.
    ! Otherwise we need to compute them from the solar composite
    ! parameterization on a equidistant grid.
    ! -------------------------------------------------------------------

    if ( yn_solar_comp ) then
      call omi_create_solcomp_irradiance (rpt_rad%nxtrack, errstat)
    else if (yn_solmonthave) then
      call omi_read_monthly_average_irradiance (errstat)
    else
      call read_irradiance_data (errstat)
    end if

  end subroutine irradiance_data_init

  ! =========================================================================

  subroutine allocate_irr_data_type (idt, max_nwave, nxtrack, errstat)

    implicit none
    type (Irradiance_Data_Type), intent(inout) :: idt
    integer (kind=i4), intent(in) :: max_nwave, nxtrack
    integer, intent(inout) :: errstat

    integer :: locerr

    if (errstat < 0) return

    allocate (idt%qflags(max_nwave, nxtrack), &
              idt%wavelengths(max_nwave, nxtrack), &
              idt%spectrum(max_nwave, nxtrack), &
              idt%nwaves(nxtrack), &
              idt%avg_wavelengths(nxtrack), &
              idt%ccdpix_selection(4, nxtrack), &
              idt%ccdpix_exclusion(2, nxtrack), &
              stat=locerr)
    if (locerr /= 0) then
      call tell_error (tell_malloc_error, &
                       "allocate_irr_data_type: allocate failed", &
                       errstat)
      return
    endif

    idt%max_nwave = max_nwave
    idt%nxtrack = nxtrack

  end subroutine allocate_irr_data_type

  ! =========================================================================

  subroutine package_irradiance_data (nwl, nxtrack, &
                                      wavelengths, spectrum, qflags, &
                                      errstat)
    ! In this routine, only the first nwl wavelengths are meaningful.
    use arrayutils, only: array_find_inner_bounding_indices_r8, &
      array_find_bounding_indices_r8
    use OMSAO_variables_module, only: &
      ctrl_fit_winexc_lim, ctrl_fit_winwav_lim

    implicit none
    integer (kind=i4), intent(in) :: nwl, nxtrack
    real (kind=r8), dimension(:, :), intent(in) :: wavelengths, spectrum
    integer (kind=i2), dimension(:,:), intent(in) :: qflags
    integer, intent(inout) :: errstat

    integer (kind=i4), dimension (:,:), allocatable :: ccdpix_sel
    integer (kind=i4) :: ix, icnt, imin, imax, max_nwavel
    integer :: locerrstat

    if (errstat < 0) return

    allocate (ccdpix_sel (4, nxtrack), stat = locerrstat)
    if (locerrstat < 0) then
      call tell_error (tell_malloc_error, "allocate_irr_data_type: allocate failed", errstat)
      return
    endif
    ccdpix_sel = -1

    ! ----------------------------------------------------
    ! Limit irradiance arrays to fitting window. Check for
    ! strictly ascending wavelengths in the process.
    ! ----------------------------------------------------
    DO ix = 1, nxtrack

      ! Restrict wavelengths to selected wavelength fitting window
      call array_find_bounding_indices_r8 (nwl, wavelengths(1:nwl, ix), &
                                           ctrl_fit_winwav_lim(1), &
                                           ctrl_fit_winwav_lim(2), &
                                           ccdpix_sel (1,ix), ccdpix_sel(2,ix))

      call array_find_bounding_indices_r8 (nwl, wavelengths(1:nwl, ix), &
                                           ctrl_fit_winwav_lim(3), &
                                           ctrl_fit_winwav_lim(4), &
                                           ccdpix_sel (3,ix), ccdpix_sel(4,ix))
    end do

    max_nwavel = 1 + MAXVAL (ccdpix_sel(4, :) - ccdpix_sel(1,:))

    call allocate_irr_data_type (Irr_Data, max_nwavel, nxtrack, errstat)
    if (errstat < 0) return

    Irr_Data%ccdpix_selection = ccdpix_sel

    do ix=1, nxtrack

      imin = ccdpix_sel(1,ix)
      imax = ccdpix_sel(4,ix)
      icnt = imax - imin + 1

      Irr_Data%wavelengths(1:icnt,ix) = wavelengths(imin:imax,ix)

      Irr_Data%spectrum(1:icnt, ix) = spectrum(imin:imax,ix)

      Irr_Data%qflags(1:icnt, ix) = qflags(imin:imax,ix)

      Irr_Data%avg_wavelengths(ix) &
        = SUM(wavelengths(imin:imax, ix))/REAL(icnt, kind=r8)

      Irr_Data%nwaves(ix) = icnt

      ! ------------------------------------------------------------------------------
      ! If any window is excluded, find the corresponding indices. This has to be done
      ! after the array assignements above because we need to know which indices to
      ! exclude from the final arrays, not the complete ones read from the HE4 file.
      ! ------------------------------------------------------------------------------
      Irr_Data%ccdpix_exclusion(1:2, ix) = -1

      if ( minval(ctrl_fit_winexc_lim(1:2)) > 0.0_r8 ) then

        call array_find_inner_bounding_indices_r8 ( &
          nwl, wavelengths(1:nwl, ix), &
          ctrl_fit_winexc_lim(1), &
          ctrl_fit_winexc_lim(2), &
          Irr_Data%ccdpix_exclusion(1, ix), &
          Irr_Data%ccdpix_exclusion(2, ix))
      end if

    end do

  end subroutine

  ! =========================================================================

  SUBROUTINE read_irradiance_data (errstat)

    USE OMSAO_precision_module
    USE OMSAO_variables_module,  ONLY: &
      l1b_irrad_filename, l1b_channel
    USE arrayutils, only: array_locate_r4
    use ctrlvars, only: yn_disable_omi_features

    !use l1bread
    use l1bread_utils
    use tio_module
    use netcdf, only : nf90_nowrite

    implicit none
    integer (kind=i4), intent (inout) :: errstat
    !
    integer :: locerrstat
    integer (kind=i4) :: nwavel, ix, nxtrack
    real (kind=r4), dimension(:,:,:), allocatable :: &
      tmp_wavelengths, tmp_spectrum
    real (kind=r8), dimension(:,:), allocatable :: &
      wavelengths, spectrum
    integer (kind=i2), dimension (:,:,:), allocatable :: tmp_qflags
    character (len=64) :: swathname
    !type (L1B_Object_Type) :: l1bobj
    type (tiof_file_type) :: tio_l1obj

    if (errstat < 0) return

    ! Allow errstat to flow

    call tell_log (1, 'reading irradiances = '//trim(l1b_irrad_filename))

    !call l1bread_swathname (l1b_irrad_filename, l1b_channel, swathname, errstat)
    !call l1bread_open_swath (l1b_irrad_filename, swathname, l1bobj, errstat)
    call tiof_open (l1b_irrad_filename, tio_l1obj, nf90_nowrite, errstat)
    call lookup_swathname (l1b_channel, swathname, errstat)
    call tiof_inq_group (tio_l1obj, swathname, errstat)
    call tiof_inq_dimlen (tio_l1obj, "xtrack", nxtrack, errstat)
    call tiof_inq_dimlen (tio_l1obj, "spectral_channel", nwavel, errstat)
    if (errstat < 0) return

    !nwavel = l1bobj%num_wavelengths
    !nxtrack = l1bobj%num_xtrack

    allocate (tmp_wavelengths(nwavel, nxtrack,1), &
              tmp_spectrum(nwavel, nxtrack,1), &
              tmp_qflags (nwavel, nxtrack,1), &
              wavelengths (nwavel, nxtrack), &
              spectrum (nwavel, nxtrack), &
              stat=locerrstat)
    if (locerrstat /= 0) then
      call tell_error (tell_malloc_error, "read_irradiance_data: allocate failed", errstat)
      return
    endif

    ! Allow errstat to flow through
    !call l1bread_get2d_r4 (l1bobj, "Irradiance", 0, 1, tmp_spectrum, errstat)
    !call l1bread_get2d_i2 (l1bobj, "PixelQualityFlags", 0, 1, tmp_qflags, errstat)
    !call l1bread_get2d_r4 (l1bobj, "Wavelength", 0, 1, tmp_wavelengths, errstat)
    !call l1bread_close (l1bobj)
    call tiof_get3d_r4 (tio_l1obj, "irradiance", [0,0,0], [1,nxtrack,nwavel], &
                        tmp_spectrum(:,1:nxtrack,1:1), errstat)
    call tiof_get3d_i2 (tio_l1obj, "pixel_quality_flag", [0,0,0], [1,nxtrack,-1], &
                        tmp_qflags(:,1:nxtrack,1:1), errstat)
    call tiof_get3d_r4 (tio_l1obj, "wavelength", [0,0,0], [1,nxtrack,nwavel], &
                        tmp_wavelengths(:,1:nxtrack,1:1), errstat)
    call tiof_close (tio_l1obj, errstat)
    if (errstat < 0) then
      call tell_error (tell_runtime_error, "read_irradiance_data:  failed reading irradiance data", &
                       errstat)
      return
    endif

    ! -------------------------------
    ! Reverse arrays for UV-1 channel.
    !   Why?  Are they stored in descending order?  --JED
    ! -------------------------------
    if (.not.yn_disable_omi_features) then
    IF ( l1b_channel == 'UV1' ) THEN
      DO ix = 1, nxtrack
        tmp_wavelengths(nwavel:1:-1, ix,1) = tmp_wavelengths(1:nwavel,ix,1)
        tmp_spectrum(nwavel:1:-1, ix,1) = tmp_spectrum(1:nwavel,ix,1)
        tmp_qflags(nwavel:1:-1, ix,1) = tmp_qflags(1:nwavel,ix,1)
      END DO
    END IF
    endif

    wavelengths = real (tmp_wavelengths(:,:,1), kind=r8)
    spectrum = real (tmp_spectrum(:,:,1), kind=r8)
    deallocate (tmp_wavelengths)
    deallocate (tmp_spectrum)

    call package_irradiance_data (nwavel, nxtrack, &
                                  wavelengths, spectrum, tmp_qflags(:,:,1), &
                                  errstat)

    return
  end subroutine read_irradiance_data

  ! ========================================================================

  SUBROUTINE omi_create_solcomp_irradiance (nxtrack, errstat)

    ! ------------------------------------------------------------------
    ! Compute an initial spectrum for an equidistant wavelength array.
    ! This will be used in the solar wavelength calibration to determine
    ! the shift of the composite solar spectrum.
    ! ------------------------------------------------------------------

    USE OMSAO_parameters_module, ONLY: i2, i4, r8, N_FIT_WINWAV
    USE OMSAO_solcomp_module,    ONLY: soco_compute
    USE arrayutils, only: array_locate_r8
    use OMSAO_variables_module, only: ctrl_fit_winwav_lim
    IMPLICIT NONE
    INTEGER (KIND=i4), INTENT (IN) :: nxtrack
    integer, intent(inout) :: errstat

    ! --------------------------------------------
    ! Spacing of the wavelength array to be set up
    ! --------------------------------------------
    REAL    (KIND=r8), PARAMETER :: dwvl = 0.1_r8

    INTEGER (KIND=i4)                         :: j, ix, nwvl
    REAL    (KIND=r8)                         :: swvl, ewvl
    real (kind=r8), dimension (:), allocatable :: tmpwvl
    real (kind=r8), dimension(:,:), allocatable :: &
      tmp_wavelengths, tmp_spectrum
    integer (kind=i2), dimension (:,:), allocatable :: tmp_qflags
    integer :: locerrstat

    if (errstat < 0) return;

    ! -----------------------------------------------------------
    ! Compute number of wavelengths and assign to temporary array
    ! -----------------------------------------------------------
    swvl = ctrl_fit_winwav_lim(1) ; ewvl = ctrl_fit_winwav_lim(N_FIT_WINWAV)
    nwvl = INT ( (ewvl-swvl) / dwvl, KIND=i4 ) + 1

    allocate (tmp_wavelengths(nwvl, nxtrack), &
              tmp_spectrum(nwvl, nxtrack), &
              tmp_qflags (nwvl, nxtrack), &
              tmpwvl(0:nwvl-1), stat=locerrstat)
    if (locerrstat /= 0) then
      call tell_error (tell_malloc_error, "omi_create_solcomp_irradiance: allocate failed", errstat)
      return
    endif

    tmpwvl = swvl + (/ (REAL(j, KIND=r8), j = 0, nwvl-1) /) * dwvl

    tmp_qflags = 0

    DO ix = 1, nxtrack

      tmp_wavelengths(:,ix) = tmpwvl

      ! ---------------------------------------------------------------
      ! Compute the solar spectrum. Note that we are not requesting the
      ! normalized spectrum here, even in cases where we DO want to use
      ! one. Rather, we are keeping the Solar Composite branch as close
      ! as possible to the regular L1b irradiance branch, which at this
      ! point is not normalized. This will be done in a later routine.
      ! ---------------------------------------------------------------
      CALL soco_compute (.FALSE., ix, nwvl, tmpwvl, tmp_spectrum(:,ix))
      !call print_array (tmp_spectrum(:,ix), nwvl)
    end do

    call package_irradiance_data (nwvl, nxtrack, &
                                  tmp_wavelengths, tmp_spectrum, tmp_qflags, &
                                  errstat)

    RETURN
  END SUBROUTINE omi_create_solcomp_irradiance

  SUBROUTINE omi_read_monthly_average_irradiance (errstat)

    ! ---------------------------------------------------------------------
    ! At this point the subroutine is ready to read the solar monthly mean
    ! irradiance from the ASCII files supplied by X. Liu, maybe in the futu
    ! re we will move to an hdf file. However, having to produce a new spec
    ! tra for each day makes that option no too charming
    ! ---------------------------------------------------------------------

    USE OMSAO_precision_module, ONLY: r4, i2
    USE OMSAO_parameters_module, ONLY: nwavel_max, nxtrack_max
    USE OMSAO_variables_module, ONLY: &
      OMSAO_solmonthave_filename, l1b_channel
    USE OMSAO_omidata_module, only: EarthSunDistance
    USE OMSAO_indices_module,   ONLY: &
      OMSAO_solmonthave_lun
    USE arrayutils, only: array_locate_r4
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: funit, ios, dummy, ix, jw, nwavel, iw
    !REAL    (KIND=r8), DIMENSION (nwavel_max,nxtrack_max) :: tmp_spc, tmp_wvl, tmp_prc
    REAL    (KIND=r8), DIMENSION (:,:), allocatable :: tmp_spc, tmp_wvl, tmp_prc
    !INTEGER (KIND=i2), DIMENSION (nwavel_max,nxtrack_max) :: tmp_flg, tmp_n
    INTEGER (KIND=i2), DIMENSION (:,:), allocatable :: tmp_flg, tmp_n
    ! ------------------------------
    ! Astronomical unit AU in meters
    ! ------------------------------
    REAL (KIND=r8), PARAMETER :: AU_m = 1.495978707d11
    REAL (KIND=r8)            :: sun_earth_distance
    ! --------------------------------------------------------------------
    ! Variable that holds the solar monthly averaged spectra
    ! --------------------------------------------------------------------
    INTEGER   (KIND = i4)     :: nxUV1
    INTEGER   (KIND = i4)     :: nxUV2
    INTEGER   (KIND = i4)     :: nxVIS
    INTEGER   (KIND = i4)     :: nwUV1
    INTEGER   (KIND = i4)     :: nwUV2
    INTEGER   (KIND = i4)     :: nwVIS
    REAL      (KIND = r4), DIMENSION(:,:,:), ALLOCATABLE :: comUV1
    REAL      (KIND = r4), DIMENSION(:,:,:), ALLOCATABLE :: comUV2
    REAL      (KIND = r4), DIMENSION(:,:,:), ALLOCATABLE :: comVIS
    INTEGER   (KIND = i4), DIMENSION(:,:),   ALLOCATABLE :: ncomUV1
    INTEGER   (KIND = i4), DIMENSION(:,:),   ALLOCATABLE :: ncomUV2
    INTEGER   (KIND = i4), DIMENSION(:,:),   ALLOCATABLE :: ncomVIS
    INTEGER (KIND=i4) :: nxtrack

    ! ------------------------
    ! Error handling variables
    ! ------------------------
    INTEGER (KIND=i4) :: version, locerrstat

    ! ---------------------------------
    ! External OMI and Toolkit routines
    ! ---------------------------------
    INTEGER (KIND=i4), EXTERNAL :: &
      pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=35), PARAMETER :: modulename = 'omi_read_monthly_average_irradiance'

    if (errstat < 0) return

    locerrstat = pge_errstat_ok

    nxtrack = 0 ; nwavel = 0 ;

    ! -------------------------
    ! Open monthly average file
    ! -------------------------
    version = 1
    locerrstat = PGS_IO_GEN_OPENF ( OMSAO_solmonthave_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
    CALL error_check ( &
      locerrstat, pgs_smf_mask_lev_s, pge_errstat_error, OMSAO_E_OPEN_SOLMONAVE_FILE, &
      modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )
    IF (  errstat /= pge_errstat_ok ) RETURN

    ! --------------------------------
    ! Reading the monthly average file
    ! --------------------------------
    ! UV1 channel
    ! -----------
    READ(UNIT=funit, FMT=*, IOSTAT=ios) nxUV1, nwUV1
    IF ( ios /= 0 ) THEN
      CALL error_check ( &
        ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
        modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )
      IF (  errstat /= pge_errstat_ok ) RETURN
    END IF

    ALLOCATE (comUV1(nwUV1,3,nxUV1))
    ALLOCATE (ncomUV1(nwUV1,nxUV1))

    DO ix = 1, nxUV1

      READ(UNIT=funit, FMT=*, IOSTAT=ios) dummy
      IF ( ios /= 0 ) THEN
        CALL error_check ( &
          ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
          modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
      END IF

      DO jw = 1, nwUV1

        READ(UNIT=funit, FMT=*, IOSTAT=ios) comUV1(jw,1,ix), &
          comUV1(jw,2,ix), comUV1(jw,3,ix), dummy, ncomUV1(jw,ix)
        IF ( ios /= 0 ) THEN
          CALL error_check ( &
            ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
            modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)),      &
            vb_lev_default, errstat )
          IF (  errstat /= pge_errstat_ok ) RETURN
        END IF

      END DO

    END DO

    ! -----------
    ! UV2 channel
    ! -----------
    READ(UNIT=funit, FMT=*, IOSTAT=ios) nxUV2, nwUV2
    IF ( ios /= 0 ) THEN
      CALL error_check ( &
        ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
        modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )
      IF (  errstat /= pge_errstat_ok ) RETURN
    END IF

    ALLOCATE (comUV2(nwUV2,3,nxUV2))
    ALLOCATE (ncomUV2(nwUV2,nxUV2))

    DO ix = 1, nxUV2

      READ(UNIT=funit, FMT=*, IOSTAT=ios) dummy
      IF ( ios /= 0 ) THEN
        CALL error_check ( &
          ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
          modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
      END IF

      DO jw = 1, nwUV2

        READ(UNIT=funit, FMT=*, IOSTAT=ios) comUV2(jw,1,ix), comUV2(jw,2,ix), comUV2(jw,3,ix), &
          dummy, ncomUV2(jw,ix)
        IF ( ios /= 0 ) THEN
          CALL error_check ( &
            ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
            modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)),      &
            vb_lev_default, errstat )
          IF (  errstat /= pge_errstat_ok ) RETURN
        END IF

      END DO

    END DO

    ! -----------
    ! VIS channel
    ! -----------
    READ(UNIT=funit, FMT=*, IOSTAT=ios) nxVIS, nwVIS
    IF ( ios /= 0 ) THEN
      CALL error_check ( &
        ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
        modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )
      IF (  errstat /= pge_errstat_ok ) RETURN
    END IF

    ALLOCATE (comVIS(nwVIS,3,nxVIS))
    ALLOCATE (ncomVIS(nwVIS,nxVIS))

    DO ix = 1, nxVIS

      READ(UNIT=funit, FMT=*, IOSTAT=ios) dummy
      IF ( ios /= 0 ) THEN
        CALL error_check ( &
          ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
          modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )
        IF (  errstat /= pge_errstat_ok ) RETURN
      END IF

      DO jw = 1, nwVIS

        READ(UNIT=funit, FMT=*, IOSTAT=ios) comVIS(jw,1,ix), comVIS(jw,2,ix), comVIS(jw,3,ix), &
          dummy, ncomVIS(jw,ix)
        IF ( ios /= 0 ) THEN
          CALL error_check ( &
            ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_SOLMONAVE_FILE, &
            modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)),      &
            vb_lev_default, errstat )
          IF (  errstat /= pge_errstat_ok ) RETURN
        END IF

      END DO

    END DO

    ! -----------------------------------------------
    ! Close monthly average file, report SUCCESS read
    ! -----------------------------------------------
    locerrstat = PGS_IO_GEN_CLOSEF ( funit )
    locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
    CALL error_check ( &
      locerrstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_SOLMONAVE_FILE, &
      modulename//f_sep//TRIM(ADJUSTL(OMSAO_solmonthave_filename)), vb_lev_default, errstat )

    IF ( errstat >= pge_errstat_error ) RETURN

    ! ---------------------------------------------
    ! Now is time to convert EarthSunDistance to AU
    ! The averaged irradiance is normalized to 1 AU
    ! so below I have to apply the scaling factor:
    ! 1.0 / sun_earth_distance / sun_earth_distance
    ! ---------------------------------------------
    sun_earth_distance = EarthSunDistance / AU_m

    allocate (tmp_spc(nwavel_max,nxtrack_max), &
              tmp_wvl(nwavel_max,nxtrack_max), &
              tmp_prc(nwavel_max,nxtrack_max), &
              tmp_flg(nwavel_max,nxtrack_max), &
              tmp_n(nwavel_max,nxtrack_max), stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, "omi_read_monthly_average_irradiance:  allocate failed", &
                       errstat)
      return
    endif
    ! ---------------------------------------------------------------
    ! Work out which channel we are interested on, UV1, UV2 or VIS
    ! Find out number of wavelengths and number of cross track pixels
    ! ---------------------------------------------------------------
    SELECT CASE (l1b_channel)
    CASE ('UV1')
      nwavel = nwUV1 ; nxtrack = nxUV1
      tmp_wvl(1:nwavel,1:nxtrack) = comUV1(1:nwavel,1,1:nxtrack)
      tmp_spc(1:nwavel,1:nxtrack) = REAL(comUV1(1:nwavel,2,1:nxtrack) / sun_earth_distance / sun_earth_distance, KIND=r4)
      DO iw = 1, nwavel
        DO ix = 1, nxtrack
          tmp_prc(iw,ix) = REAL(comUV1(iw,3,ix) / SQRT(REAL(ncomUV1(iw,ix), KIND=r8)), KIND=r4)
        END DO
      END DO
      tmp_flg(1:nwavel,1:nxtrack) = 128_i2 ! No problem with the pixel set for all of them
      tmp_n(1:nwavel, 1:nxtrack)  = INT(ncomUV1(1:nwavel,1:nxtrack),KIND=i2) !Number of irradiance averaged
    CASE ('UV2')
      nwavel = nwUV2 ; nxtrack = nxUV2
      tmp_wvl(1:nwavel,1:nxtrack) = comUV2(1:nwavel,1,1:nxtrack)
      tmp_spc(1:nwavel,1:nxtrack) = REAL( comUV2(1:nwavel,2,1:nxtrack) / sun_earth_distance / sun_earth_distance, KIND=r4)
      DO iw = 1, nwavel
        DO ix = 1, nxtrack
          tmp_prc(iw,ix) = REAL(comUV2(iw,3,ix) / SQRT(REAL(ncomUV2(iw,ix), KIND=r8)), KIND=r4)
        END DO
      END DO
      tmp_flg(1:nwavel,1:nxtrack) = 128_i2 ! No problem with the pixel set for all of them
      tmp_n(1:nwavel, 1:nxtrack)  = INT(ncomUV2(1:nwavel,1:nxtrack),KIND=i2) !Number of irradiance averaged
    CASE ('VIS')
      nwavel = nwVIS ; nxtrack = nxVIS
      tmp_wvl(1:nwavel,1:nxtrack) = comVIS(1:nwavel,1,1:nxtrack)
      tmp_spc(1:nwavel,1:nxtrack) = REAL(comVIS(1:nwavel,2,1:nxtrack) / sun_earth_distance / sun_earth_distance, KIND=r4)
      DO iw = 1, nwavel
        DO ix = 1, nxtrack
          tmp_prc(iw,ix) = REAL(comVIS(iw,3,ix) / SQRT(REAL(ncomVIS(iw,ix), KIND=r8)), KIND=r4)
        END DO
      END DO
      tmp_flg(1:nwavel,1:nxtrack) = 128_i2 ! No problem with the pixel set for all of them
      tmp_n(1:nwavel, 1:nxtrack)  = INT(ncomVIS(1:nwavel,1:nxtrack),KIND=i2) !Number of irradiance averaged
    CASE DEFAULT
      !Nothing to do here except to fold
    END SELECT

    DEALLOCATE(comUV1) ; DEALLOCATE(comUV2) ; DEALLOCATE(comVIS)
    DEALLOCATE(ncomUV1) ;  DEALLOCATE(ncomUV2) ; DEALLOCATE(ncomVIS)

    call package_irradiance_data (nwavel, nxtrack, &
                                  tmp_wvl, tmp_spc, tmp_flg, &
                                  errstat)

    ! FIXME (JCH) tmp_n unused??
    deallocate (tmp_spc,tmp_wvl,tmp_prc,tmp_flg, tmp_n)

  END SUBROUTINE omi_read_monthly_average_irradiance

end module
