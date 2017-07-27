!
module lidort_prof_utilities

  use m_ezspline_interpolation, only: bspline, bspline2, interpol, interpol2

  public hres_radcalc_env, get_hres_radcal_waves, get_hres_gascrs_ray, &
       hres_radwf_inter_convol, getabs_crs, getso2_crs, get_alb_ozcrs_ray, &
       get_slant_tau, get_tracegas_wf, radwf_interpol, get_efft, &
       lidort_prof_env
  private read_txcrs 

  integer, parameter, private :: max_pathlen = 1024

contains

  SUBROUTINE HRES_RADCALC_ENV (nw0, do_ozwf, do_albwf, do_tmpwf, do_o3shi, &
       ozvary, do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, &
       do_codwf, do_sprswf, do_so2zwf, nw, waves, nos, o3shi, sza, vza, aza, &
       nl, ozprof, tprof, n0alb, albarr, albpmin, albpmax, n0wfc, wfcarr, &
       wfcpmin, wfcpmax, nostk, albwf, ozwf, tmpwf, o3shiwf, cfracwf, codwf, &
       ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf, rad, errstat)

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY : fitvar_rad_str, numwin!, fitwavs
    USE ozprof_data_module,     ONLY : radcwav, ncalcp, &
         nalb, albidx, albmin, albmax, nwfc, wfcidx, wfcmin, wfcmax!, &
    !use_effcrs, wfcfidx, wfcfpix, wfclpix, nfalb, nfwfc, albfidx, &
    !albfpix, alblpix

    USE OMSAO_errstat_module


    IMPLICIT NONE

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN) :: nw0, nw, nl, nos, n0alb, nostk, n0wfc
    LOGICAL, INTENT(IN) :: do_ozwf, do_albwf, do_tmpwf, do_o3shi, do_taodwf, &
         do_twaewf, do_saodwf, do_cfracwf, do_codwf, do_ctpwf, do_sprswf, &
         do_so2zwf
    INTEGER, INTENT(OUT)                                 :: errstat
    INTEGER, DIMENSION(n0alb), INTENT(IN)                :: albpmax, albpmin
    INTEGER, DIMENSION(n0wfc), INTENT(IN)                :: wfcpmax, wfcpmin
    LOGICAL, DIMENSION(nl), INTENT(IN)                   :: ozvary
    REAL (KIND=dp), DIMENSION(nw),  INTENT(IN)           :: waves
    REAL (KIND=dp), DIMENSION(nw, nostk),  INTENT(OUT)   :: rad, albwf, &
         cfracwf, o3shiwf, codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf
    REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)   :: o3shi
    REAL (KIND=dp), DIMENSION(nl),  INTENT(IN)           :: ozprof, tprof
    REAL (KIND=dp), DIMENSION(nw, nl, nostk),INTENT(OUT) :: ozwf, tmpwf
    REAL (KIND=dp), DIMENSION(n0alb), INTENT(IN)         :: albarr 
    REAL (KIND=dp), DIMENSION(n0wfc), INTENT(IN)         :: wfcarr 
    REAL (KIND=dp), INTENT(IN)                           :: sza, vza, aza

    ! =======================
    ! Local variables
    ! =======================
    INTEGER :: i, j, k 
    REAL (KIND=dp), DIMENSION(nw0)            :: waves0
    REAL (KIND=dp), DIMENSION(nw0, nostk)     :: rad0, albwf0, cfracwf0, &
         o3shiwf0, codwf0, ctpwf0, taodwf0, twaewf0, saodwf0, sprswf0, so2zwf0
    REAL (KIND=dp), DIMENSION(nw0, nl, nostk) :: ozwf0, tmpwf0
    INTEGER, DIMENSION(n0alb)                 :: albpmax0, albpmin0
    INTEGER, DIMENSION(n0wfc)                 :: wfcpmax0, wfcpmin0

    errstat = pge_errstat_ok

    ! Note nw0 is MAX(ncalcp, nw)
    waves0 = 0.0
    waves0(1:ncalcp) = radcwav(1:ncalcp)

    ! Need to find indices of boundaries for using different 
    ! surface albedos/cloud fractions
    k = 0
    albpmin0(1) = 1
    albpmax0(n0alb) = nw0
    DO i = 1, nalb
      j = albidx - 1 + i
      IF (fitvar_rad_str(j)(4:4) == '0') THEN
        k = k + 1

        IF (k > 1) albpmin0(k)= albpmax0(k-1) + 1 
        IF (k < n0alb) albpmax0(k)= MINVAL(MAXLOC(waves0(1:ncalcp), MASK=(waves0(1:ncalcp) &
             >= albmin(i) .AND. waves0(1:ncalcp) < albmax(i))))  
      ENDIF
    ENDDO

    k = 0
    DO i = 1, nwfc
      j = wfcidx - 1 + i
      IF (fitvar_rad_str(j)(4:4) == '0') THEN
        k = k + 1

        IF (k > 1) wfcpmin0(k)= wfcpmax0(k-1) + 1 
        IF (k < n0wfc) wfcpmax0(k)= MINVAL(MAXLOC(waves0(1:ncalcp), MASK=(waves0(1:ncalcp) &
             >= wfcmin(i) .AND. waves0(1:ncalcp) < wfcmax(i))))  
      ENDIF
    ENDDO

    ! print *, n0alb, albpmin(1:n0alb), albpmax(1:n0alb)
    ! print *, waves(albpmin(1:n0alb)), waves(albpmax(1:n0alb))
    ! print *, nw0, ncalcp, nw
    ! print *, albpmin0(1:n0alb), albpmax0(1:n0alb)
    ! print *, hreswav(albpmin0(1:n0alb)), hreswav(albpmax0(1:n0alb))

    ! Call LIDORT_PROF_ENV on fine wavelength grids
    ! Return raidances and weighting functions on required resolution wavelength grid
    CALL LIDORT_PROF_ENV (do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary, &
         do_taodwf,do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
         do_sprswf, do_so2zwf, nw0, waves0, nos, o3shi, sza, vza, aza, nl, &
         ozprof, tprof, n0alb, albarr, albpmin0, albpmax0, n0wfc, wfcarr, &
         wfcpmin0, wfcpmax0, nostk, albwf0, ozwf0, tmpwf0, o3shiwf0, &
         cfracwf0, codwf0, ctpwf0, taodwf0, twaewf0, saodwf0, sprswf0, &
         so2zwf0, rad0, errstat)

    rad(1:nw, 1:nostk)         = rad0(1:nw, 1:nostk)
    albwf(1:nw, 1:nostk)       = albwf0(1:nw, 1:nostk)
    o3shiwf(1:nw, 1:nostk)     = o3shiwf0(1:nw, 1:nostk)
    cfracwf(1:nw, 1:nostk)     = cfracwf0(1:nw, 1:nostk)
    codwf(1:nw, 1:nostk)       = codwf0(1:nw, 1:nostk)
    ctpwf(1:nw, 1:nostk)       = ctpwf0(1:nw, 1:nostk)
    taodwf(1:nw, 1:nostk)      = taodwf0(1:nw, 1:nostk)
    twaewf(1:nw, 1:nostk)      = twaewf0(1:nw, 1:nostk)
    saodwf(1:nw, 1:nostk)      = saodwf0(1:nw, 1:nostk)
    sprswf(1:nw, 1:nostk)      = sprswf0(1:nw, 1:nostk)
    so2zwf(1:nw, 1:nostk)      = so2zwf0(1:nw, 1:nostk)
    ozwf(1:nw, 1:nl, 1:nostk)  = ozwf0(1:nw, 1:nl, 1:nostk)
    tmpwf(1:nw, 1:nl, 1:nostk) = tmpwf0(1:nw, 1:nl, 1:nostk)

    RETURN
  END SUBROUTINE HRES_RADCALC_ENV

  !1.	Establish fine wavelength grid: 0.01 nm now, may change to 0.05 nm later
  !2.	Establish radiance calculation grid, based on spectral sampling rate for 
  !   different specral regions
  !3.	Find indices of radiance calculation grid in fine wavelength grid

  SUBROUTINE get_hres_radcal_waves(errstat)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,ONLY  : max_spec_pts, max_fit_pts
    USE OMSAO_variables_module, ONLY  : numwin, nradpix, solwinfit, &
         which_slit, curr_rad_spec, use_redfixwav, winlim
    USE OMSAO_indices_module, ONLY  : hwe_idx, wvl_idx
    USE ozprof_data_module,  ONLY  : radc_nsegsr, &
         radc_samprate, radc_lambnd, hreswav, radcwav, nhresp, ncalcp, &
         radcidxs, nhresp0, hreswav0, hres_samprate!, radc_msegsr, ozabs_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! Output variables
    INTEGER, INTENT(OUT) :: errstat

    ! Local variables
    REAL (KIND=dp), PARAMETER :: dhw0 = 0.01  ! at 0.01 nm
    !INTEGER, PARAMETER        :: mextraw = 100

    INTEGER              :: i, j, fidx, lidx, nsub, nratio, nhalf!, &
    !n0, nextra, k
    REAL (KIND=dp)       :: tmp, swav, ewav, slw, samprate, invdhw!, ds1, ds2
    !REAL (KIND=dp), DIMENSION (mextraw) :: extrawave

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=21), PARAMETER :: modulename = 'get_hres_radcal_waves'

    errstat = pge_errstat_ok
    invdhw = 1.0 / dhw0

    ! Establish high resolution grid (0.01 nm grid)
    ! Do it for each fitting window
    nhresp = 0
    IF (.NOT. use_redfixwav) THEN
      DO i = 1, numwin
        swav = FLOOR(invdhw * (winlim(i, 1))) / invdhw
        IF (i > 1) THEN
          IF (swav < ewav) swav = ewav + dhw0
        ENDIF
        ewav = CEILING(invdhw * (winlim(i, 2))) / invdhw 

        nsub = NINT((ewav - swav) / dhw0) + 1  
        hreswav(nhresp+1:nhresp+nsub) = swav + dhw0 * (/(j, j = 0, nsub-1)/)      
        nhresp = nhresp + nsub
      ENDDO
    ELSE 
      ! If use fixed wavelengths, then go through it one by one
      fidx = 1
      DO i = 1, numwin
        lidx = fidx + nradpix(i) - 1
        slw = solwinfit(i, hwe_idx, 1)
        IF (which_slit == 0) THEN
          !slw = slw * 2.63  ! to truncate those < 0.001
          slw = slw * 2.2  ! to truncate those  < 0.01
        ELSE
          slw = slw * 2.0
        ENDIF

        DO j = fidx, lidx
          swav = FLOOR(invdhw * (curr_rad_spec(wvl_idx, j) - slw)) / invdhw
          IF (j > 1) THEN
            IF (swav < ewav) swav = ewav + invdhw
          ENDIF
          ewav = CEILING(invdhw * (curr_rad_spec(wvl_idx, j) + slw)) / invdhw

          nsub = NINT((ewav - swav) / dhw0) + 1  
          hreswav(nhresp+1:nhresp+nsub) = swav + dhw0 * (/(j, j = 0, nsub-1)/)
          nhresp = nhresp + nsub
        ENDDO

        fidx = lidx + 1
      ENDDO
    ENDIF
    nhresp0 = nhresp
    hreswav0(1:nhresp) = hreswav(1:nhresp)
    invdhw = 1.0 / hres_samprate
    nratio = INT(hres_samprate / dhw0)
    nhalf =  nratio / 2

    j = 1
    DO i = nhalf + 1, nhresp0 - nhalf, nratio
      fidx = i - nhalf
      lidx = i - nhalf + nratio - 1
      hreswav(j) = SUM(hreswav(fidx:lidx)) / REAL(nratio)
      j = j + 1
    ENDDO
    nhresp = j - 1

    IF (nhresp0 > max_spec_pts .OR. nhresp > max_spec_pts) THEN
      WRITE(www_lun, *) modulename, ': Need to increase max_spec_pts!!!'
      errstat = pge_errstat_error
      RETURN
    ENDIF

    ! Establish radiance calculation grid, based on spectral sampling
    ! Sampling rate are specified in several segments (e.g., < 295 nm, 295-308 nm, > 308 nm)
    ! Make sure that radiance will be done for the first and last points. 
    radc_samprate = FLOOR(radc_samprate * invdhw) / invdhw  ! multiples of dhw
    ncalcp = 1
    radcwav(1) = hreswav(1)
    fidx = 2
    DO i = 1, radc_nsegsr
      samprate = radc_samprate(i)

      IF (i == radc_nsegsr) THEN
        lidx = nhresp - 1
      ELSE
        lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) < radc_lambnd(i+1))))
      ENDIF

      DO j = fidx, lidx
        tmp = ABS(hreswav(j) - radcwav(ncalcp) - samprate)
        IF (tmp < hres_samprate * 0.1) THEN
          ncalcp = ncalcp + 1
          radcwav(ncalcp) = hreswav(j)
        ELSE IF (hreswav(j + 1) >= radcwav(ncalcp) + samprate * 2 &
             .AND. hreswav(j)   >  radcwav(ncalcp) + samprate * 0.5) THEN
          ncalcp = ncalcp + 1
          radcwav(ncalcp) = hreswav(j)
        ENDIF
      ENDDO

      fidx = lidx + 1
    ENDDO
    IF (hreswav(nhresp) > radcwav(ncalcp) + samprate * 0.5) THEN
      ncalcp = ncalcp + 1
      radcwav(ncalcp) = hreswav(nhresp)
    ELSE
      radcwav(ncalcp) = hreswav(nhresp)
    ENDIF

    !DO i = 1, nhresp
    !   WRITE(90, *) hreswav(i)
    !ENDDO
    !
    !DO i = 1, ncalcp
    !   WRITE(91, *) radcwav(i)
    !ENDDO

    !! **** May add some spectral points later where large errors can occur **
    !OPEN(UNIT = ozabs_unit, file='INP/wave_o3abs_minmax.dat', status='old')
    !READ(ozabs_unit, *) nextra
    !READ(ozabs_unit, *) extrawave(1:nextra)
    !CLOSE(ozabs_unit)
    !
    !IF (nextra > mextraw) THEN
    !   WRITE(www_lun, *) modulename, ': Need to increase mextraw!!!'
    !   errstat = pge_errstat_error
    !ENDIF
    !
    !! Insert these wavelengths
    !j = 1; n0 = ncalcp
    !DO i = 1, n0
    !   IF (radcwav(i) > extrawave(j)) THEN
    !      ds2 = radcwav(i) - extrawave(j)
    !      ds1 = extrawave(j) - radcwav(i-1) 
    !      IF (ds1 > hres_samprate * 2.0 .AND. ds2 > hres_samprate * 2.0 ) THEN ! Add this wavelength
    !         radcwav(i + 1 : ncalcp + 1) = radcwav(i : ncalcp)
    !         radcwav(i) = radcwav(i-1) + NINT( ds1 / hres_samprate) * hres_samprate 
    !         ncalcp = ncalcp + 1
    !      ENDIF
    !      
    !      j = j + 1
    !   ENDIF
    !ENDDO

    IF (ncalcp > max_fit_pts) THEN
      WRITE(www_lun, *) modulename, ': Need to increase max_fit_pts!!!'
      errstat = pge_errstat_error
      RETURN
    ENDIF

    ! Find the indices of radiance calc. wavelength in high resolution
    radcidxs(1:ncalcp) = 0
    radcidxs(1) = 1
    j = 2
    DO i = 2, nhresp
      IF ( ABS(hreswav(i) - radcwav(j)) <= hres_samprate * 0.2 ) THEN
        radcidxs(j) = i
        j = j + 1
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE get_hres_radcal_waves



  ! Prepare high resolution spectra at fine grid: solar reference, trace 
  ! gas cross sections (o3 shift and o3 temperature), Raleigh scattering 
  ! coefficient, depolarization factor O3/SO2 (use_so2dtcrs=.TRUE.) cross 
  ! section: if do_tmpwf = .FALSE. and do_o3shi is false, just need to get 
  ! once for each retrieval Other trace gas cross section: just need to 
  ! get it once for all the retrievals if no shifts 
  SUBROUTINE GET_HRES_GASCRS_RAY(nw, waves, nz, ts, do_o3shi, o3shi, &
       do_tmpwf, nfgas, use_so2dtcrs, num_iter, allcrs, allcol, rhos, &
       abscrs_qtdepen, raycof, depol, pge_error_status)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,ONLY  : max_spec_pts, zerok
    USE OMSAO_variables_module, ONLY  : numwin, refspec_orig_data, &
         n_refspec_pts, winlim, fitvar_rad, rmask_fitvar_rad, refspec_norm!, &
    !which_slit, use_redfixwav, solwinfit, nradpix, curr_rad_spec
    USE OMSAO_indices_module,   ONLY : wvl_idx, solar_idx, spc_idx, &
         so2_idx, so2v_idx!, hwe_idx
    USE ozprof_data_module,     ONLY : mxsect, nos, ngas, gasidxs, fgasidxs, &
         fgassidxs, nos, hreswav, radcidxs, ncalcp, hres_gas, hres_gasshi, &
         hres_o3, hres_o3shi, hres_so2, hres_so2shi, o3crsz, o3dadsz, &
         o3dadtz, so2crsz, so2dads, oswins, do_subfit, hres_i0, hres_raycof, &
         hres_depol, wrtozcrs, osfind, hresgabs, hresray, nhresp, nhresp0, &
         hreswav0, hres_samprate
    USE OMSAO_errstat_module
    use m_lidort_prof_prep, only: get_all_raycof_depol
    use utilities, only: interpolation

    IMPLICIT NONE

    ! Input/output variables
    INTEGER, INTENT (IN)       :: nw, nz, nfgas, num_iter ! nw = ncalcp
    INTEGER, INTENT (OUT)      :: pge_error_status
    LOGICAL, INTENT (IN)       :: do_o3shi, do_tmpwf, use_so2dtcrs
    REAL (KIND=dp), DIMENSION(nw), INTENT (IN )           :: waves
    REAL (KIND=dp), DIMENSION(nz), INTENT (IN )           :: ts, rhos
    REAL (KIND=dp), DIMENSION(nw), INTENT (OUT)           :: raycof, depol
    REAL (KIND=dp), DIMENSION(3, nw), INTENT (OUT)        :: abscrs_qtdepen
    REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)    :: o3shi
    REAL (KIND=dp), DIMENSION(nfgas, nz), INTENT(IN)      :: allcol
    REAL (KIND=dp), DIMENSION(nw, nfgas, nz), INTENT(OUT) :: allcrs

    ! Local variables
    INTEGER :: i, j, fidx, lidx, npts, idum, nfgas1, errstat, nratio, nhalf!, k
    LOGICAL                              :: problems, do_shi
    REAL (KIND=dp)                       :: tmp, frac, thet, so2sum
    REAL (KIND=dp), DIMENSION (nhresp)   :: delshi, tmpwav, delpos
    REAL (KIND=dp), DIMENSION (nhresp,nz):: so2dadsz   

    ! Save original O3/SO2 cross sections
    LOGICAL, SAVE :: o3tdepend, so2tdepend, do_so2shi, first = .TRUE.
    INTEGER, SAVE :: no3t, nso2t, no3, nso2, &
         so2sfidx, so2vsfidx!, so2idx, so2vidx
    REAL (KIND=dp), SAVE :: o3normc, so2normc
    REAL (KIND=dp), DIMENSION(mxsect),                 SAVE :: o3ts, so2ts
    REAL (KIND=dp), DIMENSION(0:mxsect, max_spec_pts), SAVE :: o3crs0, so2crs0 ! Original

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=19), PARAMETER :: modulename = 'GET_HRES_GASCRS_RAY'

    pge_error_status = pge_errstat_ok

    IF (first) THEN

      ! Obtain high resolution solar reference spectra
      npts = n_refspec_pts(solar_idx)
      CALL interpolation (npts, refspec_orig_data(solar_idx,1:npts,wvl_idx), &
           refspec_orig_data(solar_idx,1:npts,spc_idx), nhresp0, &
           hreswav0(1:nhresp0), hres_i0(1:nhresp0), errstat)   
      IF ( errstat > pge_errstat_warning ) THEN
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF

      nratio = INT(hres_samprate / 0.01)
      nhalf =  nratio / 2

      j = 1
      DO i = nhalf + 1, nhresp0 - nhalf, nratio
        fidx = i - nhalf
        lidx = i - nhalf + nratio - 1
        hres_i0(j) = SUM(hres_i0(fidx:lidx)) / REAL(nratio)
        j = j + 1
      ENDDO


      ! Obtain high resolution rayleigh scattering coefficients and depolarization factor
      CALL GET_ALL_RAYCOF_DEPOL(nhresp, hreswav(1:nhresp), hres_raycof(1:nhresp), hres_depol(1:nhresp))

      ! Obtain high resolution cross sections of other trace gases (except for O3 and SO2)
      do_so2shi = .FALSE.
      so2sfidx = 0
      so2vsfidx = 0
      DO i = 1, ngas
        IF (fgasidxs(i) > 0 ) THEN
          IF (gasidxs(i) == so2_idx) so2sfidx = fgassidxs(i)
          IF (gasidxs(i) == so2v_idx) so2vsfidx = fgassidxs(i)

          IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. fgassidxs(i) > 0) do_so2shi = .TRUE.
          IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. use_so2dtcrs) CYCLE
          IF (fgassidxs(i) > 0 ) CYCLE

          hres_gas(i, 1:nhresp) = 0.0D0
          npts = n_refspec_pts(gasidxs(i))
          fidx = MINVAL(MINLOC(hreswav(1:nhresp), MASK = (hreswav(1:nhresp) >= &
               refspec_orig_data(gasidxs(i), 1, 1) .AND. hreswav(1:nhresp)  &
               <= refspec_orig_data(gasidxs(i), npts, 1))))
          lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK = (hreswav(1:nhresp) >= &
               refspec_orig_data(gasidxs(i), 1, 1) .AND. hreswav(1:nhresp) &
               <= refspec_orig_data(gasidxs(i), npts, 1))))

          IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN 
            CALL BSPLINE(refspec_orig_data(gasidxs(i), 1:npts, 1), &
                 refspec_orig_data(gasidxs(i), 1:npts, 2), npts, hreswav(fidx:lidx), &
                 hres_gas(i, fidx:lidx), lidx - fidx + 1, errstat) 

            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ' : BSPLINE2 error, errstat = ', errstat
              pge_error_status = pge_errstat_error
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      ! Obtain original O3 cross section (quadratic or individual T, shift, T-depen)
      CALL READ_TXCRS(hreswav(1:nhresp), nhresp, 1, mxsect, o3crs0, no3, o3ts, &
           o3tdepend, no3t, o3normc, problems)
      IF (problems) THEN
        WRITE(www_lun, *) modulename, ' : Error in reading ozone cross sections!!!'
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF

      IF (.NOT. do_o3shi) THEN
        DO i = 1, no3t
          CALL BSPLINE(o3crs0(0, 1:no3), o3crs0(i, 1:no3), no3, hreswav(1:nhresp), &
               hres_o3(i, 1:nhresp), nhresp, errstat)
          IF (errstat < 0) THEN
            WRITE(www_lun, *) modulename, ': BSPLINE2 error, errstat = ', errstat
            pge_error_status = pge_errstat_error
            RETURN
          ENDIF
        ENDDO
      ENDIF

      ! Obtain original SO2 cross section (quadratic or individual T, shift)
      IF (use_so2dtcrs) THEN
        CALL READ_TXCRS(hreswav(1:nhresp), nhresp, 2, mxsect, so2crs0, nso2, so2ts, &
             so2tdepend, nso2t, so2normc, problems)
        IF (problems) THEN
          WRITE(www_lun, *) modulename, ' : Error in reading SO2 cross sections!!!'
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF

        IF (.NOT. do_so2shi) THEN
          DO i = 1, nso2t
            CALL BSPLINE(so2crs0(0, 1:nso2), so2crs0(i, 1:nso2), nso2, hreswav(1:nhresp), &
                 hres_so2(i, 1:nhresp), nhresp, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': BSPLINE2 error, errstat = ', errstat
              pge_error_status = pge_errstat_error
              RETURN
            ENDIF
          ENDDO
        ENDIF

      ENDIF

      first = .FALSE.  
    ENDIF

    ! Interpolate ozone cross section
    IF (do_o3shi .AND. nos > 0) THEN

      ! determine ozone wavelength shifts
      delshi = 0.0
      IF (do_subfit) THEN
        fidx = 1
        DO j = 1, numwin
          IF (j == numwin) THEN
            lidx = nhresp
          ELSE
            tmp = (winlim(j, 2) + winlim(j+1, 1)) / 2.0
            lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) <= tmp)))
          ENDIF
          delpos(fidx:lidx) =  hreswav(fidx:lidx) - (hreswav(fidx) + hreswav(lidx)) / 2.0
          IF (osfind(j, 1) > 0) delshi(fidx:lidx) =  o3shi(j, 1) 

          DO i = 2, nos
            IF (osfind(j, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx)  + &
                 o3shi(j, i) * delpos(fidx:lidx) ** (i-1)
          ENDDO
          fidx = lidx + 1
        ENDDO
      ELSE
        IF (oswins(1, 1) == 1) THEN
          fidx = 1
        ELSE
          tmp = (winlim(oswins(1, 1), 1) + winlim(oswins(1, 1) - 1, 2))/2.
          fidx = MINVAL(MINLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) >= tmp)))
        ENDIF

        IF (oswins(1, 2) == numwin) THEN
          lidx = nhresp
        ELSE
          tmp = (winlim(oswins(1, 2), 2) + winlim(oswins(1, 2) + 1, 1))/2.
          lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) <= tmp)))
        ENDIF

        delpos(fidx:lidx) =  hreswav(fidx:lidx) - (hreswav(fidx) + hreswav(lidx)) / 2.0
        IF (osfind(1, 1) > 0) delshi(fidx:lidx) =  + o3shi(1, 1) 

        DO i = 2, nos  
          IF (osfind(1, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx) + &
               o3shi(1, i) * delpos(fidx:lidx) ** (i-1)
        ENDDO
      ENDIF

      DO i = 1, no3t
        CALL BSPLINE2(o3crs0(0, 1:no3), o3crs0(i, 1:no3), no3, do_o3shi, hreswav(1:nhresp)-delshi, &
             hres_o3(i, 1:nhresp), hres_o3shi(i, 1:nhresp), nhresp, errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': BSPLINE2 error, errstat = ', errstat
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
      ENDDO
    ENDIF

    ! Get ozone cross section at each layer
    IF (num_iter == 0 .OR. do_o3shi .OR. do_tmpwf) THEN
      IF (o3tdepend .AND. no3t == 3) THEN     ! qudratic T dependent coefficients
        DO i = 1, nz
          thet = ts(i) - zerok
          o3crsz(1:nhresp, i) = (hres_o3(1, 1:nhresp) + hres_o3(2, 1:nhresp) * thet &
               + hres_o3(3, 1:nhresp) * thet * thet )
          IF (do_o3shi) THEN
            o3dadsz(1:nhresp, i) = (hres_o3shi(1, 1:nhresp) + hres_o3shi(2, 1:nhresp) * thet &
                 + hres_o3shi(3, 1:nhresp) * thet * thet )
          ENDIF
          IF (do_tmpwf) THEN
            o3dadtz(1:nhresp, i) = (hres_o3(2, 1:nhresp)  + 2.0 * hres_o3(3, 1:nhresp) * thet)
          ENDIF
        ENDDO
      ELSE IF (o3tdepend .AND. no3t == 2) THEN   ! linear T dependent coefficients
        DO i = 1, nz
          thet = ts(i) - zerok
          o3crsz(1:nhresp, i) = hres_o3(1, 1:nhresp) + hres_o3(2, 1:nhresp) * thet 
          IF (do_o3shi) o3dadsz(1:nhresp, i) = (hres_o3shi(1, 1:nhresp) + hres_o3shi(2, 1:nhresp) * thet)
          IF (do_tmpwf) o3dadtz(1:nhresp, i) = hres_o3(2, 1:nhresp)
        ENDDO
      ELSE IF (no3t == 1)  THEN           ! only 1 T
        DO i = 1, nz
          o3crsz(1:nhresp, i) = hres_o3(1, 1:nhresp)
          IF (do_o3shi) o3dadsz(1:nhresp, i) = hres_o3shi(1, 1:nhresp)
          IF (do_tmpwf) o3dadtz(1:nhresp, i) = 0.0
        ENDDO
      ELSE IF (no3t == 2) THEN            ! have 2 T values
        DO i = 1, nz
          frac = 1.0 - (ts(i) - o3ts(1)) / (o3ts(2) - o3ts(1))
          o3crsz(1:nhresp, i) = (frac * hres_o3(1, 1:nhresp) + (1.0 - frac) * hres_o3(2, 1:nhresp))
          IF (do_o3shi) o3dadsz(1:nhresp, i) = (frac * hres_o3shi(1, 1:nhresp) + &
               (1.0 - frac) * hres_o3shi(2, 1:nhresp))
          IF (do_tmpwf) o3dadtz(1:nhresp, i) = (hres_o3(1, 1:nhresp) - hres_o3(2, 1:nw)) / (o3ts(1) - o3ts(2))
        END DO
      ELSE  IF (.NOT. o3tdepend .AND. no3t > 3) THEN  ! have more than n T     
        DO i = 1, nhresp                     ! Interpolate/extrapolate over T   
          CALL INTERPOL2(o3ts(1:no3t), hres_o3(1:no3t,i), no3t, do_tmpwf, ts, &
               o3crsz(i, 1:nz), o3dadtz(i, 1:nz),nz, errstat)
          IF (errstat < 0) THEN
            WRITE(www_lun, *) modulename, ': INTERPOL2 error, errstat = ', errstat
            pge_error_status = pge_errstat_error
            RETURN
          ENDIF

          IF (do_o3shi) THEN 
            CALL INTERPOL(o3ts(1:no3t), hres_o3shi(1:no3t,i), no3t, ts, &
                 o3dadsz(i, 1:nz), nz, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
              pge_error_status = pge_errstat_error
              RETURN
            ENDIF
          ENDIF
        ENDDO
      ELSE 
        WRITE(www_lun, *) modulename, ': Such type of ozone cross sections not implemented' 
        pge_error_status = pge_errstat_error
        RETURN
      ENDIF
      IF (do_tmpwf) THEN
        o3dadtz(1:nhresp, 1:nz) = o3dadtz(1:nhresp, 1:nz) / o3crsz(1:nhresp, 1:nz)
      ENDIF
      IF (do_o3shi) THEN
        o3dadsz(1:nhresp, 1:nz) = o3dadsz(1:nhresp, 1:nz) / o3crsz(1:nhresp, 1:nz)
      ENDIF
      o3crsz(1:nhresp, 1:nz) = o3crsz(1:nhresp, 1:nz) * o3normc
      IF (wrtozcrs .AND. o3tdepend .AND. no3t <= 3) &
           abscrs_qtdepen(1:3, 1:ncalcp) = hres_o3(1:3, radcidxs(1:ncalcp)) * o3normc
    ENDIF
    allcrs(1:ncalcp, 1, 1:nz) = o3crsz(radcidxs(1:ncalcp), 1:nz)
    DO i = 1, nz
      hresgabs(1:nhresp, i) = o3crsz(1:nhresp, i) * allcol(1, i)
    ENDDO

    ! Get SO2 cross sections
    IF (use_so2dtcrs) THEN
      IF (do_so2shi) THEN
        DO i = 1, nso2t
          idum = MAX(so2sfidx, so2vsfidx)
          tmp = fitvar_rad(rmask_fitvar_rad(idum))
          CALL BSPLINE2(so2crs0(0, 1:nso2), so2crs0(i, 1:nso2), nso2, do_so2shi, &
               hreswav(1:nhresp) - tmp, hres_so2(i, 1:nhresp), hres_so2shi(i, 1:nhresp), &
               nhresp, errstat)
          IF (errstat < 0) THEN
            WRITE(www_lun, *) modulename, ': BSPLINE2 error, errstat = ', errstat
            pge_error_status = pge_errstat_error
            RETURN
          ENDIF
        ENDDO
      ENDIF

      IF (num_iter == 0 .OR. do_so2shi .OR. do_tmpwf) THEN
        IF (so2tdepend .AND. nso2t == 3) THEN     ! qudratic T dependent coefficients
          DO i = 1, nz
            thet = ts(i) - zerok
            so2crsz(1:nhresp, i) = (hres_so2(1, 1:nhresp) + hres_so2(2, 1:nhresp) * thet &
                 + hres_so2(3, 1:nhresp) * thet * thet )
            IF (do_so2shi) THEN
              so2dadsz(1:nhresp, i) = (hres_so2shi(1, 1:nhresp) + hres_so2shi(2, 1:nhresp) * thet &
                   + hres_so2shi(3, 1:nhresp) * thet * thet )
            ENDIF
          ENDDO
        ELSE IF (so2tdepend .AND. nso2t == 2) THEN   ! linear T dependent coefficients
          DO i = 1, nz
            thet = ts(i) - zerok
            so2crsz(1:nhresp, i) = hres_so2(1, 1:nhresp) + hres_so2(2, 1:nhresp) * thet 
            IF (do_so2shi) so2dadsz(1:nhresp, i) = (hres_so2shi(1, 1:nhresp) + hres_so2shi(2, 1:nhresp) * thet)
          ENDDO
        ELSE IF (nso2t == 1)  THEN           ! only 1 T
          DO i = 1, nz
            so2crsz(1:nhresp, i) = hres_so2(1, 1:nhresp)
            IF (do_so2shi) so2dadsz(1:nhresp, i) = hres_so2shi(1, 1:nhresp)
          ENDDO
        ELSE IF (nso2t == 2) THEN            ! have 2 T values
          DO i = 1, nz
            frac = 1.0 - (ts(i) - so2ts(1)) / (so2ts(2) - so2ts(1))
            so2crsz(1:nhresp, i) = (frac * hres_so2(1, 1:nhresp) + (1.0 - frac) * hres_so2(2, 1:nhresp))
            IF (do_so2shi) so2dadsz(1:nhresp, i) = (frac * hres_so2shi(1, 1:nhresp) + &
                 (1.0 - frac) * hres_so2shi(2, 1:nhresp))
          END DO
        ELSE IF (.NOT. so2tdepend .AND. nso2t >= 3) THEN  ! have more than n T     
          DO i = 1, nhresp         ! interpolate/extrapolate over T   
            CALL INTERPOL(so2ts(1:nso2t), hres_so2(1:nso2t,i), nso2t, ts,  &
                 so2crsz(i, 1:nz), nz, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
              pge_error_status = pge_errstat_error
              RETURN
            ENDIF

            IF (do_so2shi) THEN 
              CALL INTERPOL(so2ts(1:nso2t), hres_so2shi(1:nso2t,i), nso2t, ts, &
                   so2dadsz(i, 1:nz), nz, errstat)
              IF (errstat < 0) THEN
                WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
                pge_error_status = pge_errstat_error
                RETURN
              ENDIF
            ENDIF
          ENDDO
        ELSE 
          WRITE(www_lun, *) modulename, ': Such type of SO2 cross sections not implemented' 
          pge_error_status = pge_errstat_error
          RETURN
        ENDIF
        IF (do_so2shi) so2dadsz(1:nhresp, 1:nz) = so2dadsz(1:nhresp, 1:nz) / so2crsz(1:nhresp, 1:nz)
        so2crsz(1:nhresp, 1:nz) = so2crsz(1:nhresp, 1:nz) * so2normc
      ENDIF
    ENDIF

    ! Obtain high resolution cross sections of other trace gases 
    nfgas1 = 1
    IF (do_so2shi) THEN
      so2dads(1:nhresp) = 0.D0
      so2sum = 0.D0
    ENDIF
    DO i = 1, ngas
      IF (fgasidxs(i) > 0 ) THEN
        nfgas1 = nfgas1 + 1

        IF ((gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) .AND. use_so2dtcrs) THEN
          allcrs(1:ncalcp, nfgas1, 1:nz) = so2crsz(radcidxs(1:ncalcp), 1:nz) / refspec_norm(gasidxs(i))
          IF (do_so2shi) THEN
            DO j = 1, nz
              so2dads(1:nhresp) = so2dads(1:nhresp) + so2dadsz(1:nhresp, j) * allcol(nfgas1, j)
              so2sum = so2sum + allcol(nfgas1, j)
            ENDDO
          ENDIF
          DO j = 1, nz
            hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + so2crsz(1:nhresp, j) &
                 / refspec_norm(gasidxs(i)) * allcol(nfgas1, j)
          ENDDO
        ELSE
          IF (fgassidxs(i) > 0) THEN
            do_shi = .TRUE.
            idum = rmask_fitvar_rad(fgassidxs(i))
            tmpwav = hreswav(1:nhresp) - fitvar_rad(idum)  ! wavelength shifts
            hres_gas(i, 1:nhresp) = 0.0D0
            npts = n_refspec_pts(gasidxs(i))
            fidx = MINVAL(MINLOC(tmpwav(1:nhresp), MASK = (tmpwav(1:nhresp) >= &
                 refspec_orig_data(gasidxs(i), 1, 1) .AND. tmpwav(1:nhresp)  &
                 <= refspec_orig_data(gasidxs(i), npts, 1))))
            lidx = MINVAL(MAXLOC(tmpwav(1:nhresp), MASK = (tmpwav(1:nhresp) >= &
                 refspec_orig_data(gasidxs(i), 1, 1) .AND. tmpwav(1:nhresp) &
                 <= refspec_orig_data(gasidxs(i), npts, 1))))

            IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN 
              CALL BSPLINE2(refspec_orig_data(gasidxs(i), 1:npts, 1), &
                   refspec_orig_data(gasidxs(i), 1:npts, 2), npts, do_shi, tmpwav(fidx:lidx), &
                   hres_gas(i, fidx:lidx), hres_gasshi(i, fidx:lidx), lidx - fidx + 1, errstat) 

              IF (errstat < 0) THEN
                WRITE(www_lun, *) modulename, ' : BSPLINE2 error, errstat = ', errstat
                pge_error_status = pge_errstat_error
                RETURN
              ENDIF
            ENDIF
          ENDIF

          DO j = 1, nz
            allcrs(1:ncalcp, nfgas1, j) = hres_gas(i, radcidxs(1:ncalcp))
            hresgabs(1:nhresp, j) = hresgabs(1:nhresp, j) + hres_gas(i, 1:nhresp) * allcol(nfgas1, j)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    IF (do_so2shi .AND. use_so2dtcrs) so2dads(1:nhresp) = so2dads(1:nhresp) / so2sum

    ! Rayleigh scattering and depolarization factor
    raycof(1:ncalcp) = hres_raycof(radcidxs(1:ncalcp))
    depol(1:ncalcp)  = hres_depol(radcidxs(1:ncalcp))

    ! compute total rayleigh and absorption optical thickness (o3 + other gas) for 
    ! later radiance interpolaiton
    DO i = 1, nz 
      hresray(1:nhresp, i) = rhos(i) * hres_raycof(1:nhresp)
    ENDDO

    RETURN
  END SUBROUTINE GET_HRES_GASCRS_RAY

  SUBROUTINE hres_radwf_inter_convol(nw, nz, nctp, ncbp, nsprs, faerlvl,  &
       do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf, &
       do_tracewf, use_so2dtcrs, do_o3shi, do_tmpwf, wave, ozs, rad, &
       fozwf, albwf, cfracwf, faerwf, faerswf, fcodwf, fsprswf, fraywf, &
       dads, dadt, abscrs, so2crs, errstat)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,   ONLY  : so2_idx, so2v_idx
    USE OMSAO_variables_module, ONLY  : numwin, &
         winlim, owave=>radwvl_sav, now=>n_radwvl_sav, i0sav, &
         nrad=>n_rad_wvl,do_bandavg, refidx_sav!, nradpix, fitwavs, &
    !refidx, database, database_shiwf, curr_rad_spec, band_selectors
    USE ozprof_data_module,     ONLY  : hwave=>hreswav, &
         radcidxs, hres_i0, nhw=>nhresp, hresgabs, hresray, nw0=>ncalcp, &
         o3crsz, o3dadtz, o3dadsz, so2crsz, ngas, gasidxs, fgasidxs, &
         fgassidxs!, so2dads, nup2p, radcwav, hres_gasshi, hres_gas
    USE OMSAO_errstat_module
    use m_gauss, only: convol_f2c
    use avg_band, only: avg_band_spec

    IMPLICIT NONE

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN) :: nw, nz, nctp, ncbp, faerlvl, nsprs
    INTEGER, INTENT(OUT) :: errstat
    LOGICAL, INTENT(IN) :: do_albwf, do_faerwf, do_faerswf, do_codwf, &
         do_sprswf, do_cfracwf, do_o3shi, do_tmpwf, do_tracewf, use_so2dtcrs

    REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
    REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fozwf, faerwf, &
         faerswf, fcodwf, fsprswf, fraywf
    REAL (KIND=dp), DIMENSION(nw, nz), INTENT(OUT)   :: dads, dadt, abscrs, &
         so2crs
    REAL (KIND=dp), DIMENSION(nw),     INTENT(IN)    :: wave
    REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: rad, albwf, cfracwf

    ! Local variables
    INTEGER :: i, iwin, fidx, lidx, fidxc, lidxc, idx, iw, ntemp, &
         sidx, eidx!, nspec, j
    LOGICAL :: do_so2shi
    REAL (KIND=dp)                      :: temp
    INTEGER, DIMENSION (nw)             :: c2hfidx, c2hlidx

    REAL (KIND=dp), DIMENSION (nhw)     :: hrad, halbwf, dtau, dray, &
         hcfracwf!, tmparr
    REAL (KIND=dp), DIMENSION (now)     :: oi0, otmp, tmpi0!, so2dads1
    REAL (KIND=dp), DIMENSION (nw, nz)  :: tauwf
    REAL (KIND=dp), DIMENSION (now, nz) :: dads1, dadt1, abscrs1, so2crs1
    !REAL (KIND=dp), DIMENSION (ngas,now):: tmp_gas, tmp_gasshi
    REAL (KIND=dp), DIMENSION (nhw, nz) :: hozwf, haerwf, haerswf, hcodwf, &
         hsprswf!, hraywf 
    REAL (KIND=dp), DIMENSION (nhw, nz*8) :: inarr
    REAL (KIND=dp), DIMENSION (now, nz*8) :: outarr

    !INTEGER :: ntime = 1

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    !CHARACTER (LEN=17), PARAMETER :: modulename = 'hres_radwf_convol'

    errstat = pge_errstat_ok

    ! get weighting function in dlnI/dx and take the logarithm of radiances
    DO i = 1, nz
      fozwf(1:nw0, i) = fozwf(1:nw0, i) / rad(1:nw0)
    ENDDO

    IF (do_albwf) THEN
      albwf(1:nw0) = albwf(1:nw0) / rad(1:nw0)
    ENDIF

    IF (do_cfracwf) THEN
      cfracwf(1:nw0) = cfracwf(1:nw0) / rad(1:nw0)
    ENDIF

    IF (do_faerwf) THEN
      DO i = faerlvl, nz
        faerwf(1:nw0, i) = faerwf(1:nw0, i) / rad(1:nw0)
      ENDDO
    ENDIF

    IF (do_faerswf) THEN
      DO i = faerlvl, nz
        faerswf(1:nw0, i) = faerswf(1:nw0, i) / rad(1:nw0)
      ENDDO
    ENDIF

    IF (do_codwf) THEN
      DO i = nctp, ncbp
        fcodwf(1:nw0, i) = fcodwf(1:nw0, i) / rad(1:nw0)
      ENDDO
    ENDIF

    IF (do_sprswf) THEN
      DO i = nsprs, nz
        fsprswf(1:nw0, i) = fsprswf(1:nw0, i) / rad(1:nw0)
      ENDDO
    ENDIF

    DO i = 1, nz
      fraywf(1:nw0, i) = fraywf(1:nw0, i) / rad(1:nw0)
    ENDDO
    rad(1:nw0) = LOG(rad(1:nw0))

    ! convert ozone weighting function to gas absorption weighting function
    DO i = 1, nz
      tauwf(1:nw0, i) = fozwf(1:nw0, i) * ozs(i) / hresgabs(radcidxs(1:nw0), i)
    ENDDO

    c2hfidx(1) = 1
    fidxc = 1
    DO iwin = 1, numwin
      IF (iwin == numwin) THEN
        lidx = nhw
        lidxc = nw0
      ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(hwave(1:nhw), MASK=(hwave(1:nhw) <= temp)))
        lidxc = MINVAL(MAXLOC(wave(1:nw0), MASK=(wave(1:nw0) <= temp)))
      ENDIF

      ! Find range of indices that map coarse-grid to fine grid
      DO i = fidxc, lidxc
        IF (i < lidxc) THEN
          temp = (wave(i) + wave(i+1)) / 2.0
          idx = MINVAL(MAXLOC(hwave(1:nhw), MASK=(hwave(1:nhw) <= temp)))
          c2hlidx(i) = idx 
        ELSE
          c2hlidx(i) = lidx
        ENDIF

        IF (i > 1) THEN
          c2hfidx(i) = c2hlidx(i-1) + 1
        ENDIF
      ENDDO

      fidxc = lidxc + 1
    ENDDO

    ! Perform correction
    ! Radiance: use o3/tau weighting function
    ! O3 weighting function: same scaled by o3 absorption cross section
    ! Rayleigh/surface pressure/other weighting function: cublic interpolatin
    DO iw = 1, nw0
      ! Correction for radiance
      fidx = c2hfidx(iw)
      lidx = c2hlidx(iw)
      hrad(fidx:lidx) = rad(iw)
      DO i = 1, nz
        dtau(fidx:lidx) = hresgabs(fidx:lidx, i) - hresgabs(radcidxs(iw), i)
        dray(fidx:lidx) = hresray(fidx:lidx, i)  - hresray(radcidxs(iw), i)
        hrad(fidx:lidx) = hrad(fidx:lidx) + tauwf(iw, i) * dtau(fidx:lidx) + &
             fraywf(iw, i) * dray(fidx:lidx)

        ! Is it better to do cubic interpolation or assume same wf (after normalized by o3 xsec)
        hozwf(fidx:lidx, i) = fozwf(iw, i) * o3crsz(fidx:lidx, i) / o3crsz(radcidxs(iw), i)
      ENDDO
    ENDDO

    !DO i = 1, nz 
    !   CALL BSPLINE(wave(1:nw0), fozwf(1:nw0, i), nw0, hwave(1:nhw), hozwf(1:nhw, i), nhw, errstat)
    !   IF (errstat < 0) THEN
    !      WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
    !      errstat = pge_errstat_error; RETURN
    !   ENDIF
    !ENDDO

    IF (do_albwf) THEN
      CALL BSPLINE(wave(1:nw0), albwf(1:nw0), nw0, hwave(1:nhw), halbwf(1:nhw), nhw, errstat)
    ENDIF

    IF (do_cfracwf) THEN
      CALL BSPLINE(wave(1:nw0), cfracwf(1:nw0), nw0, hwave(1:nhw), hcfracwf(1:nhw), nhw, errstat)
    ENDIF

    IF (do_faerwf) THEN
      DO i = faerlvl, nz
        CALL BSPLINE(wave(1:nw0), faerwf(1:nw0, i), nw0, hwave(1:nhw), haerwf(1:nhw, i), nhw, errstat)
      ENDDO
    ENDIF

    IF (do_faerswf) THEN
      DO i = faerlvl, nz
        CALL BSPLINE(wave(1:nw0), faerswf(1:nw0, i), nw0, hwave(1:nhw), haerswf(1:nhw, i), nhw, errstat)
      ENDDO
    ENDIF

    IF (do_codwf) THEN
      DO i = nctp, ncbp
        CALL BSPLINE(wave(1:nw0), fcodwf(1:nw0, i), nw0, hwave(1:nhw), hcodwf(1:nhw, i), nhw, errstat)
      ENDDO
    ENDIF

    IF (do_sprswf) THEN
      DO i = nsprs, nz
        CALL BSPLINE(wave(1:nw0), fsprswf(1:nw0, i), nw0, hwave(1:nhw), hsprswf(1:nhw, i), nhw, errstat)
      ENDDO
    ENDIF

    ! Convert radiances back
    hrad = EXP(hrad)

    ! convert radiance/weighting function to dlnI/dx from dI/dx
    DO i = 1, nz
      hozwf(:, i) = hozwf(:, i) * hrad
    ENDDO

    IF (do_albwf) THEN
      halbwf = halbwf * hrad
    ENDIF

    IF (do_cfracwf) THEN
      hcfracwf = hcfracwf * hrad
    ENDIF

    IF (do_faerwf) THEN
      DO i = faerlvl, nz
        haerwf(:, i) = haerwf(:, i) * hrad
      ENDDO
    ENDIF

    IF (do_faerswf) THEN
      DO i = faerlvl, nz
        haerswf(:, i) = haerswf(:, i) * hrad
      ENDDO
    ENDIF

    IF (do_codwf) THEN
      DO i = nctp, ncbp
        hcodwf(:, i) = hcodwf(:, i) * hrad
      ENDDO
    ENDIF

    IF (do_sprswf) THEN
      DO i = nsprs, nz
        hsprswf(:, i) = hsprswf(:, i) * hrad
      ENDDO
    ENDIF


    ! Convolve radiance/weighting functions/solar reference at fine grids into measurement grid 
    ! convolve all spectra at once to speed up computation
    ! *** First, transfer all spectra to inarr ***
    inarr(1:nhw, 1) = hres_i0(1:nhw)
    inarr(1:nhw, 2) = hres_i0(1:nhw) * hrad(1:nhw)
    DO i = 1, nz
      sidx = 2 + i
      inarr(1:nhw, sidx) = hres_i0(1:nhw) * hozwf(1:nhw, i)
    ENDDO
    IF (do_albwf) THEN
      sidx = sidx + 1
      inarr(1:nhw, sidx) = hres_i0(1:nhw) * halbwf(1:nhw)
    ENDIF
    IF (do_cfracwf) THEN
      sidx = sidx + 1
      inarr(1:nhw, sidx) = hres_i0(1:nhw) * hcfracwf(1:nhw)
    ENDIF
    IF (do_faerwf) THEN
      DO i = faerlvl, nz
        sidx = sidx + 1
        inarr(1:nhw, sidx) = hres_i0(1:nhw) * haerwf(1:nhw, i)
      ENDDO
    ENDIF
    IF (do_faerswf) THEN
      DO i = faerlvl, nz
        sidx = sidx + 1
        inarr(1:nhw, sidx) = hres_i0(1:nhw) * haerswf(1:nhw, i)
      ENDDO
    ENDIF
    IF (do_codwf) THEN
      DO i = nctp, ncbp
        sidx = sidx + 1
        inarr(1:nhw, sidx) = hres_i0(1:nhw) * hcodwf(1:nhw, i)
      ENDDO
    ENDIF
    IF (do_sprswf) THEN
      DO i = nsprs, nz
        sidx = sidx + 1
        inarr(1:nhw, sidx) = hres_i0(1:nhw) * hsprswf(1:nhw, i)
      ENDDO
    ENDIF

    ! convolve ozone shift/temperature 
    IF (do_tracewf) THEN
      DO i = 1, nz
        sidx = sidx + 1
        inarr(1:nhw, sidx) = o3crsz(1:nhw, i) * hres_i0(1:nhw) 
      ENDDO

      do_so2shi = .FALSE.
      DO i = 1, ngas
        IF (fgasidxs(i) > 0) THEN
          IF (gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN
            IF (fgassidxs(i) > 0 ) do_so2shi = .TRUE.
            IF (use_so2dtcrs) CYCLE
          ENDIF

          ! This is not necessary: could still use those effective cross sections
          ! sidx = sidx + 1;  inarr(1:nhw, sidx) = hres_gas(i, 1:nhw) * hres_i0(1:nhw) 
          !
          !IF (fgassidxs(i) > 0) THEN
          ! sidx = sidx + 1;  inarr(1:nw, sidx) = hres_gasshi(i, 1:nhw)
          !ENDIF
        ENDIF
      ENDDO

      IF (use_so2dtcrs) THEN
        DO i = 1, nz
          sidx = sidx + 1
          inarr(1:nhw, sidx) = so2crsz(1:nhw, i) * hres_i0(1:nhw) 
        ENDDO

        !IF (do_so2shi) THEN
        !   sidx = sidx + 1;  inarr(1:nhw, sidx) = so2dads(1:nhw) 
        !ENDIF
      ENDIF
    ENDIF

    IF (do_o3shi) THEN
      DO i = 1, nz
        sidx = sidx + 1
        inarr(1:nhw, sidx) = o3dadsz(1:nhw, i) !* hres_i0(1:nhw)
      ENDDO
    ENDIF
    IF (do_tmpwf) THEN
      DO i = 1, nz
        sidx = sidx + 1
        inarr(1:nhw, sidx) = o3dadtz(1:nhw, i) !* hres_i0(1:nhw)
      ENDDO
    ENDIF
    eidx = sidx

    ! *** second, convole all spectra at once ****
    CALL convol_f2c(hwave(1:nhw), inarr(1:nhw, 1:eidx), nhw, eidx, owave(1:now), outarr(1:now, 1:eidx), now)

    ! *** third, transfer all convolved spectra back ***
    oi0(1:now)  = outarr(1:now, 1)
    hrad(1:now) = outarr(1:now, 2) / oi0(1:now)
    DO i = 1, nz
      sidx = 2 + i
      hozwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
    ENDDO
    IF (do_albwf) THEN
      sidx = sidx + 1
      halbwf(1:now) = outarr(1:now, sidx) / oi0(1:now)
    ENDIF
    IF (do_cfracwf) THEN
      sidx = sidx + 1
      hcfracwf(1:now) = outarr(1:now, sidx) / oi0(1:now)
    ENDIF
    IF (do_faerwf) THEN
      DO i = faerlvl, nz
        sidx = sidx + 1
        haerwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
      ENDDO
    ENDIF
    IF (do_faerswf) THEN
      DO i = faerlvl, nz
        sidx = sidx + 1
        haerswf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
      ENDDO
    ENDIF
    IF (do_codwf) THEN
      DO i = nctp, ncbp
        sidx = sidx + 1
        hcodwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
      ENDDO
    ENDIF
    IF (do_sprswf) THEN
      DO i = nsprs, nz
        sidx = sidx + 1
        hsprswf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
      ENDDO
    ENDIF

    ! convolve ozone shift/temperature 
    IF (do_tracewf) THEN
      DO i = 1, nz
        sidx = sidx + 1
        abscrs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
      ENDDO

      do_so2shi = .FALSE.
      DO i = 1, ngas
        IF (fgasidxs(i) > 0) THEN
          IF (gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN
            IF (fgassidxs(i) > 0 ) do_so2shi = .TRUE.
            IF (use_so2dtcrs) CYCLE
          ENDIF

          ! This is not necessary: could still use those effective cross sections
          ! sidx = sidx + 1;  tmp_gas(i, 1:now) = outarr(1:now, sidx) / oi0(1:now) 
          !
          !IF (fgassidxs(i) > 0) THEN
          ! sidx = sidx + 1;  tmp_gasshi(i, 1:now) = outarr(1:now, sidx) 
          !ENDIF
        ENDIF
      ENDDO

      IF (use_so2dtcrs) THEN
        DO i = 1, nz
          sidx = sidx + 1
          so2crs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO

        !IF (do_so2shi) THEN
        ! sidx = sidx + 1;  so2dads1(1:now) = outarr(1:now, sidx) 
        !ENDIF
      ENDIF
    ENDIF

    IF (do_o3shi) THEN
      DO i = 1, nz
        sidx = sidx + 1
        dads1(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
      ENDDO
    ENDIF
    IF (do_tmpwf) THEN
      DO i = 1, nz
        sidx = sidx + 1
        dadt1(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
      ENDDO
    ENDIF

    ! Perform additional coadding (if necesary)
    ! Use OMI solar spectra here
    IF (do_bandavg) THEN
      tmpi0(1:now) = i0sav(refidx_sav(1:now))  ! better to use OMI solar spectra     
      oi0(1:now) = tmpi0
      CALL avg_band_spec(owave(1:now), oi0(1:now), now, ntemp, errstat)
      IF (ntemp /= nrad .OR. errstat /= 0) THEN
        WRITE(www_lun, *) 'Spectra Averaging Error: ', now, ntemp, nrad 
        errstat = pge_errstat_error
        RETURN
      ENDIF

      otmp(1:now) = hrad(1:now) * tmpi0(1:now)
      CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
      rad(1:nrad) = otmp(1:nrad) / oi0(1:nrad)

      DO i = 1, nz
        otmp(1:now) = hozwf(1:now, i) * tmpi0(1:now)
        CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        fozwf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
      ENDDO

      !DO i = 1, nz
      !   otmp(1:now) = hraywf(1:now, i) * tmpi0(1:now)
      !   CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
      !   fraywf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
      !ENDDO

      IF (do_albwf) THEN
        otmp(1:now) = halbwf(1:now) * tmpi0(1:now)
        CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        albwf(1:nrad) = otmp(1:nrad) / oi0(1:nrad)
      ENDIF

      IF (do_cfracwf) THEN
        otmp(1:now) = hcfracwf(1:now) * tmpi0(1:now)
        CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        cfracwf(1:nrad) = otmp(1:nrad) / oi0(1:nrad)
      ENDIF

      IF (do_faerwf) THEN
        DO i = faerlvl, nz
          otmp(1:now) = haerwf(1:now, i) * tmpi0(1:now)
          CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          faerwf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO
      ENDIF

      IF (do_faerswf) THEN
        DO i = faerlvl, nz
          otmp(1:now) = haerswf(1:now, i) * tmpi0(1:now)
          CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          faerswf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO
      ENDIF

      IF (do_codwf) THEN
        DO i = nctp, ncbp
          otmp(1:now) = hcodwf(1:now, i) * tmpi0(1:now)
          CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          fcodwf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO
      ENDIF

      IF (do_sprswf) THEN
        DO i = nsprs, nz
          otmp(1:now) = hsprswf(1:now, i) * tmpi0(1:now)
          CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          fsprswf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO
      ENDIF

      IF (do_o3shi) THEN
        DO i = 1, nz
          otmp(1:now) = dads1(1:now, i)  !* tmpi0(1:now)
          CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          dads(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
        ENDDO
      ENDIF

      IF (do_tmpwf) THEN
        DO i = 1, nz
          otmp(1:now) = dadt1(1:now, i) !* tmpi0(1:now)
          CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          dadt(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
        ENDDO
      ENDIF

      IF (do_tracewf) THEN
        DO i = 1, nz
          otmp(1:now) = abscrs1(1:now, i) * tmpi0(1:now)
          CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          abscrs(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO

        !DO i = 1, ngas
        !   IF (fgasidxs(i) > 0) THEN
        !      IF (gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN
        !         IF (use_so2dtcrs) CYCLE
        !      ENDIF
        !      
        !      otmp(1:now) = tmp_gas(i, 1:now) * tmpi0(1:now)
        !      CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        !      database(gasidxs(i), refidx(1:nrad)) = otmp(1:nrad) / oi0(1:nrad)
        !      
        !      IF (fgassidxs(i) > 0) THEN
        !         otmp(1:now) = tmp_gasshi(i, 1:now) 
        !         CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        !         database_shiwf(gasidxs(i), refidx(1:nrad)) = otmp(1:nrad) 
        !      ENDIF
        !   ENDIF
        !ENDDO

        IF (use_so2dtcrs) THEN
          DO i = 1, nz
            otmp(1:now) = so2crs1(1:now, i) * tmpi0(1:now) 
            CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
            so2crs(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
          ENDDO

          !IF (do_so2shi) THEN
          !   otmp(1:now) = so2dads1(1:now) 
          !   CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
          !   database_shiwf(so2_idx, refidx(1:nrad))  = otmp(1:nrad)  
          !   database_shiwf(so2v_idx, refidx(1:nrad)) = otmp(1:nrad) 
          !ENDIF
        ENDIF
      ENDIF

    ELSE
      rad(1:now) = hrad(1:now)
      fozwf(1:now, 1:nz) = hozwf(1:now, 1:nz)
      !fraywf(1:now, 1:nz) = hraywf(1:now, 1:nz)
      IF (do_albwf) albwf(1:now) = halbwf(1:now)
      IF (do_cfracwf) cfracwf(1:now) = hcfracwf(1:now)
      IF (do_faerwf) faerwf(1:now, faerlvl:nz) = haerwf(1:now, faerlvl:nz)
      IF (do_faerwf) faerswf(1:now, faerlvl:nz) = haerswf(1:now, faerlvl:nz)
      IF (do_codwf) fcodwf(1:now, nctp:ncbp) = hcodwf(1:now, nctp:ncbp)
      IF (do_sprswf) fsprswf(1:now, nsprs:nz) = hcodwf(1:now, nsprs:nz)
      IF (do_o3shi) dads(1:now, 1:nz) = dads1(1:now, 1:nz)
      IF (do_tmpwf) dadt(1:now, 1:nz) = dadt1(1:now, 1:nz)
      IF (do_tracewf) THEN
        abscrs(1:now, 1:nz) = abscrs1(1:now, 1:nz)
        !DO i = 1, ngas
        !   IF (fgasidxs(i) > 0) THEN
        !      IF (gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN
        !         IF (use_so2dtcrs) CYCLE
        !      ENDIF
        !      database(gasidxs(i), refidx(1:now)) = tmp_gas(i, 1:now)
        !      IF (fgassidxs(i) > 0) database(gasidxs(i), refidx(1:now)) = tmp_gasshi(i, 1:now)
        !   ENDIF
        !ENDDO
        IF (use_so2dtcrs) THEN
          so2crs(1:now, 1:nz) = so2crs1(1:now, 1:nz)
          !IF (do_so2shi) THEN
          !   database_shiwf(so2_idx, refidx(1:now)) = so2dads1(1:now)  
          !   database_shiwf(so2v_idx, refidx(1:now)) = so2dads1(1:now) 
          !ENDIF
        ENDIF
      ENDIF
    ENDIF

    RETURN

  END SUBROUTINE hres_radwf_inter_convol


  ! Read T-dependent cross sections (O3 and SO2)
  ! Inteprolate original cross section to input wavelength grid
  SUBROUTINE READ_TXCRS(waves, nw0, which_gas, maxt, abscrs, nw, ts, tdepend, nt, normc, problems)

    USE OMSAO_precision_module  
    USE OMSAO_parameters_module, ONLY : max_spec_pts
    USE OMSAO_variables_module,  ONLY : refdbdir
    USE ozprof_data_module,      ONLY : ozabs_fname, ozabs_unit
    USE OMSAO_errstat_module
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN)                                :: which_gas, maxt, nw0
    REAL (KIND=dp), INTENT(IN), DIMENSION(nw0)         :: waves

    ! Output variables
    INTEGER, INTENT(OUT)                               :: nt, nw
    REAL (KIND=dp), INTENT(OUT)                        :: normc
    REAL (KIND=dp), DIMENSION(0:maxt, max_spec_pts), INTENT(OUT) :: abscrs
    REAL (KIND=dp), DIMENSION(maxt), INTENT(OUT)       :: ts
    LOGICAL, INTENT(OUT)                               :: problems, tdepend

    ! Local variables
    INTEGER                                            :: nline, i, j!, errstat
    LOGICAL                                            :: slitconv
    REAL (KIND=dp)                                     :: maxw, minw
    CHARACTER (LEN=14)                                 :: tmpchar
    CHARACTER (LEN=max_pathlen)                                :: absfname

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    !CHARACTER (LEN=10), PARAMETER    :: modulename = 'READ_TXCRS'

    problems = .FALSE.
    tmpchar = ' ' 

    IF (which_gas == 1) THEN
      absfname = ozabs_fname
    ELSE IF (which_gas == 2) THEN
      absfname = TRIM(ADJUSTL(refdbdir)) // 'OMSAO_SO2_scia_fm.dat'
    ELSE
      WRITE(www_lun, *) 'This cross section not implemented!!!'
      problems = .TRUE.
      RETURN
    ENDIF

    OPEN(UNIT = ozabs_unit, file=TRIM(ADJUSTL(absfname)), status='old')

    DO WHILE (tmpchar /= 'START OF TABLE') 
      READ (ozabs_unit, '(A14)') tmpchar
    ENDDO

    READ (ozabs_unit, *) nline, minw, maxw, normc  
    IF (minw > waves(1) .OR. maxw < waves(nw0)) THEN
      !print *, waves(1), waves(nw0), minw, maxw
      WRITE(www_lun, *) 'Need to increase reference cross section range!!!'
      problems = .TRUE.
      CLOSE (ozabs_unit)
      RETURN
    ENDIF

    READ (ozabs_unit, *) tdepend, nt, slitconv
    IF (nt > maxt) THEN
      WRITE(www_lun, *) 'Need to increase parameter mxsect!!!'
      problems = .TRUE.
      CLOSE (ozabs_unit)
      RETURN
    ENDIF

    IF (.NOT. slitconv .AND. which_gas == 1) THEN
      WRITE(www_lun, *) 'Need to use high-resolution cross section for O3!!!'
      problems = .TRUE.
      CLOSE (ozabs_unit)
      RETURN
    ENDIF

    IF (tdepend .AND. nt > 3) THEN
      WRITE(www_lun, *) 'High-order (> 2) T-dependent x-section not implemented!!!'
      problems = .TRUE.
      CLOSE (ozabs_unit)
      RETURN
    ENDIF

    READ (ozabs_unit, *) ts(1:nt)
    IF ((nt > 1 .AND. .NOT. tdepend ) .AND. (MINVAL(ts(1:nt)) > 220. &
         .OR.  MAXVAL(ts(1:nt)) < 280. )) THEN
      WRITE(www_lun, *) 'Temperature range for X-section not enough!!!'
      problems = .TRUE.
      CLOSE (ozabs_unit)
      RETURN
    ENDIF

    j = 1
    DO i = 1, nline
      READ(ozabs_unit, *) abscrs(0, j), abscrs(1:nt, j)
      IF (abscrs(0, j) > waves(1) - 1.0 .AND. abscrs(0, j) < &
           waves(nw0) + 1.0)  j = j + 1
    ENDDO
    nw = j - 1
    CLOSE (ozabs_unit)

    IF (nw > max_spec_pts) THEN
      WRITE(www_lun, *) 'Need to increase parameter max_spec_pts!!!'
      problems = .TRUE.
      CLOSE (ozabs_unit)
      RETURN
    ENDIF

    !abscrs(1:nt, 1:nw) = abscrs(1:nt, 1:nw) * normc

    RETURN

  END SUBROUTINE READ_TXCRS

  ! ===========================================================
  !  Author: Xiong Liu
  !  Date: Jan. 30, 2004
  !  Purpose: Get absorption cross section for several species 
  ! ===========================================================

  SUBROUTINE GETABS_CRS(lamda, nlsav, nlamda, whichgas,  &
       nlayers, tsgrid, abscrs, dods, dodt, dads, dadt, problems, crsqtd_k)

    USE OMSAO_precision_module  
    USE OMSAO_parameters_module,ONLY : zerok
    USE OMSAO_variables_module, ONLY : do_bandavg, &
         n_refspec_pts, refspec_orig_data, ozabs_convl!, yn_varyslit, &
    !which_slit, refwvl_sav, refdbdir, refidx_sav, n_refwvl_sav, i0sav
    USE ozprof_data_module, ONLY : ozabs_fname, ozabs_unit, wrtozcrs!, &
    !nflay, fozs
    USE OMSAO_slitfunction_module
    use avg_band, only: avg_band_effozcrs
    use m_dataspline, only: correct_i0effect

    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: nlamda, nlsav, nlayers, whichgas
    LOGICAL, INTENT(IN) :: dods, dodt
    REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)   :: lamda
    REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers) :: tsgrid

    ! Output variables
    REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: abscrs, dads, dadt
    REAL (KIND=dp), DIMENSION(3, nlamda), INTENT(OUT)       :: crsqtd_k
    LOGICAL, INTENT(OUT) :: problems

    ! Local variables
    INTEGER, PARAMETER :: ngas = 1   ! O3, NO2,O4,BrO,So2, HCHO, CLO2
    INTEGER, PARAMETER :: maxline  = 10001      ! # of wavelengths
    INTEGER, PARAMETER :: maxt = 5              ! # of Ts or # of coeff.

    INTEGER, SAVE        :: nline, nt
    LOGICAL, SAVE        :: first = .TRUE., tdepend, slitconv
    REAL (KIND=dp), SAVE :: maxw, minw, normc
    REAL (KIND=dp), SAVE, DIMENSION(maxline)       :: refwavs, lowresi0
    REAL (KIND=dp), SAVE, DIMENSION(maxt)          :: ts
    REAL (KIND=dp), SAVE, DIMENSION(maxt, maxline) :: refabs0, refabs

    REAL (KIND=dp), DIMENSION(maxt, nlsav)   :: savabs, savabs_d1
    REAL (KIND=dp), DIMENSION(maxt, nlamda)  :: tmpabs, tmpabs_d1

    !REAL (KIND=dp), DIMENSION(nlsav)         :: tmpi0

    INTEGER            :: i, j, errstat, ntemp, ni0!, nline1
    REAL (KIND=dp)     :: thet, frac, scalex
    LOGICAL            :: get_lresi0

    CHARACTER (LEN=max_pathlen):: absfname
    CHARACTER (LEN=14) :: tmpchar

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=10), PARAMETER    :: modulename = 'getabs_crs'

    problems = .FALSE.
    IF ( whichgas > ngas ) THEN
      WRITE(www_lun, *) 'Absorption cross section for this gas not available!!!'
      problems = .TRUE.
      RETURN
    ENDIF

    ! Read data
    IF (whichgas == 1)  THEN ! O3
      absfname = ozabs_fname
    ENDIF


    IF (first) THEN
      tmpchar = ' ' 
      OPEN(UNIT = ozabs_unit, file=absfname, status='old')
      DO WHILE (tmpchar /= 'START OF TABLE') 
        READ (ozabs_unit, '(A14)') tmpchar
      ENDDO

      READ (ozabs_unit, *) nline, minw, maxw, normc

      IF (minw > lamda(1) .OR. maxw < lamda(nlsav)) THEN
        WRITE(www_lun, *) 'Need to increase reference cross section range!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      READ (ozabs_unit, *) tdepend, nt, slitconv
      IF (nt > maxt) THEN
        WRITE(www_lun, *) 'Need to increase parameter maxt!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      IF (tdepend .AND. nt > 3) THEN
        WRITE(www_lun, *) 'High-order (>2) T-dependent cross section not implemented!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      READ (ozabs_unit, *) ts(1:nt)
      IF ((nt > 1 .AND. .NOT. tdepend ) .AND. (MINVAL(ts(1:nt)) > &
           MINVAL(tsgrid) .OR.  MAXVAL(ts(1:nt)) < MAXVAL(tsgrid))) THEN
        WRITE(www_lun, *) 'Temperature range for cross section not enough!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      j = 1
      DO i = 1, nline
        READ(ozabs_unit, *) refwavs(j), refabs0(1:nt, j)
        IF (refwavs(j) > lamda(1) - 5.0 .AND. refwavs(j) < &
             lamda(nlsav) + 5.0)  j = j + 1
      ENDDO
      nline = j - 1
      CLOSE (ozabs_unit)

      IF (nline > maxline) THEN
        WRITE(www_lun, *) 'Need to increase parameter maxline!!!'
        problems = .TRUE.
        RETURN
      ENDIF

      first = .FALSE.
    ENDIF

    IF (ozabs_convl) THEN   

      refabs(1:nt, 1:nline) = refabs0(1:nt, 1:nline)

      ! Perform solar i0 effect on ozone cross-section (no need to convolve)
      IF (slitconv) THEN
        ni0 = n_refspec_pts(1)
        scalex = 0.2 ! ~600 DU SUM(fozs(1:nflay)) * 2.69E16 * normc, now a dummy number, not used
        DO i = 1, nt 
          IF (i == 1) get_lresi0 = .TRUE.
          CALL CORRECT_I0EFFECT(refwavs(1:nline), refabs(i, 1:nline), nline, refspec_orig_data(1, 1:ni0, 1), &
               refspec_orig_data(1, 1:ni0, 2), ni0, scalex, get_lresi0, errstat, lowresi0(1:nline))
          IF (i == 1) get_lresi0 = .FALSE.

          IF ( errstat /= 0 ) THEN
            WRITE(www_lun, *) 'Error in Correct I0 Effect!!!'
            problems = .TRUE.
            RETURN
          ENDIF
        ENDDO
      ENDIF

      ozabs_convl = .FALSE.
    ENDIF

    dads = 0.0
    dadt = 0.0  

    ! Use high-resolution solar reference for weighting 
    ! (difference is very small between using different i0)
    !IF (do_bandavg) THEN
    !   CALL BSPLINE(refwavs(1:nline), lowresi0(1:nline), nline, refwvl_sav, &
    !        i0sav(1:n_refwvl_sav), n_refwvl_sav, errstat)
    !
    !   IF (errstat < 0) THEN
    !      WRITE(www_lun, *) modulename, ': BSPLINE2 error, errstat = ', errstat
    !      problems = .TRUE.; RETURN
    !   ENDIF
    !ENDIF

    ! Convolution and interpolation
    DO i = 1, nt
      CALL BSPLINE2(refwavs(1:nline), refabs(i, 1:nline), nline, dods, lamda, &
           savabs(i, :), savabs_d1(i, :), nlsav, errstat)
      IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename, ': BSPLINE2 error, errstat = ', errstat
        problems = .TRUE.
        RETURN
      ENDIF

      IF (do_bandavg) THEN
        CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)

        IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
          WRITE(www_lun, *) 'O3 Spectra Averaging Error: ', nlsav, nlamda, ntemp
          problems = .TRUE.
          RETURN
        ENDIF
        tmpabs(i, :) = savabs(i, 1:nlamda)

        CALL avg_band_effozcrs(lamda, savabs_d1(i, :), nlsav, ntemp, errstat)
        IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
          WRITE(www_lun, *) 'O3 Spectra Averaging Error: ', nlsav, nlamda, ntemp
          problems = .TRUE.
          RETURN 
        ENDIF
        tmpabs_d1(i, :) = savabs_d1(i, 1:nlamda)

      ELSE
        tmpabs(i, :) = savabs(i, 1:nlamda)
        tmpabs_d1(i, :) = savabs_d1(i, 1:nlamda)
      ENDIF

    ENDDO

    IF (wrtozcrs .AND. tdepend .AND. nt <= 3) THEN
      crsqtd_k(1:3, 1:nlamda) = tmpabs(1:3, 1:nlamda) * normc
    ENDIF

    ! interpolate with respect to temperature
    IF (tdepend .AND. nt == 3) THEN        ! qudratic T dependent coefficients
      DO i = 1, nlayers
        thet = tsgrid(i) - zerok
        abscrs(:, i) = (tmpabs(1, :) + tmpabs(2, :) * thet &
             + tmpabs(3, :) * thet * thet )
        IF (dods) THEN
          dads(:, i) = (tmpabs_d1(1, :) + tmpabs_d1(2, :) * thet &
               + tmpabs_d1(3, :) * (thet**2))
        ENDIF
        IF (dodt) THEN
          dadt(:, i) = (tmpabs(2, :)  + 2.0 * tmpabs(3, :) * thet)
        ENDIF
      ENDDO
    ELSE IF (tdepend .AND. nt == 2) THEN   ! linear T dependent coefficients
      DO i = 1, nlayers
        thet = tsgrid(i) - zerok
        abscrs(:, i) = (tmpabs(1, :) + tmpabs(2, :) * thet)
        IF (dods) dads(:, i) = (tmpabs_d1(1, :) + tmpabs_d1(2, :) * thet)
        IF (dodt) dadt(:, i) = tmpabs(2, :)
      ENDDO
    ELSE IF (nt == 1)  THEN           ! only 1 T
      DO i = 1, nlayers
        abscrs(:, i) = tmpabs(1, :)
        IF (dods) dads(:, i) = tmpabs_d1(1, :)
        IF (dodt) dadt(:, i) = 0.0
      ENDDO
    ELSE IF (nt == 2) THEN            ! have 2 T values
      DO i = 1, nlayers
        frac = 1.0 - (tsgrid(i) - ts(1)) / (ts(2) - ts(1))
        abscrs(:, i) = (frac * tmpabs(1, :) + (1.0 - frac) * tmpabs(2, :))
        IF (dods) dads(:, i) = (frac * tmpabs_d1(1, :) + &
             (1.0 - frac) * tmpabs_d1(2, :))
        IF (dodt) dadt(:, i) = (tmpabs(1, :) - tmpabs(2, :)) / (ts(1) - ts(2))
      END DO
    ELSE IF (.NOT. tdepend .AND. nt >= 3) THEN    ! have more than n T     
      DO i = 1, nlamda           ! interpolate over T   
        CALL INTERPOL2(ts(1:nt), tmpabs(1:nt,i), nt, dodt, tsgrid,&
             abscrs(i, :), dadt(i, :), nlayers, errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': INTEPROL2 error, errstat = ', errstat
          problems = .TRUE.
          RETURN
        ENDIF

        IF (dods) THEN 
          CALL INTERPOL(ts(1:nt), tmpabs_d1(1:nt,i), nt, tsgrid, &
               dads(i, :), nlayers, errstat)
          IF (errstat < 0) THEN
            WRITE(www_lun, *) modulename, ': INTERPOL error, errstat = ', errstat
            problems = .TRUE.
            RETURN
          ENDIF
        ENDIF
      ENDDO
    ELSE 
      WRITE(www_lun, *) modulename, ': Such type of ozone cross sections not implemented' 
      problems = .TRUE.
      RETURN
    ENDIF

    IF (dods) dads = dads / abscrs   ! get relative sensitivty to shift
    IF (dodt) dadt = dadt / abscrs   ! get relative sensitivity to T
    abscrs = abscrs * normc

    RETURN  
  END SUBROUTINE GETABS_CRS

  SUBROUTINE GETSO2_CRS(lamda, nlsav, nlamda, nlayers, tsgrid, abscrs, &
       problems)

    USE OMSAO_precision_module  
    USE OMSAO_parameters_module,ONLY : zerok
    USE OMSAO_variables_module, ONLY : do_bandavg, refdbdir, &
         n_refspec_pts, refspec_orig_data, so2crs_convl!, &
    !yn_varyslit, which_slit, refwvl_sav, n_refwvl_sav, refidx_sav
    USE ozprof_data_module,     ONLY : ozabs_unit!, mxsect
    USE OMSAO_slitfunction_module
    use avg_band, only: avg_band_effozcrs
    use m_dataspline, only: correct_i0effect

    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN)  :: nlamda, nlsav, nlayers
    REAL (KIND=dp), INTENT(IN), DIMENSION(nlsav)            :: lamda
    REAL (KIND=dp), INTENT(IN), DIMENSION(nlayers)          :: tsgrid

    ! Output variables
    REAL (KIND=dp), DIMENSION(nlamda, nlayers), INTENT(OUT) :: abscrs
    LOGICAL, INTENT(OUT)                                    :: problems

    ! Local variables
    INTEGER, PARAMETER   :: maxline  = 2000      ! # of wavelengths
    INTEGER, PARAMETER   :: maxt = 5             ! # of Ts or # of coeff.

    INTEGER, SAVE        :: nline, nt
    LOGICAL, SAVE        :: first = .TRUE., tdepend, slitconv
    REAL (KIND=dp), SAVE :: maxw, minw, normc
    REAL (KIND=dp), SAVE, DIMENSION(maxline)       :: refwavs, lowresi0
    REAL (KIND=dp), SAVE, DIMENSION(maxt)          :: ts
    REAL (KIND=dp), SAVE, DIMENSION(maxt, maxline) :: refabs0, refabs

    REAL (KIND=dp), DIMENSION(maxt, nlsav)   :: savabs
    REAL (KIND=dp), DIMENSION(maxt, nlamda)  :: tmpabs

    !REAL (KIND=dp), DIMENSION(nlsav)         :: tmpi0

    INTEGER            :: i, j, errstat, ntemp, ni0!, nline1
    REAL (KIND=dp)     :: thet, frac, scalex
    LOGICAL            :: get_lresi0

    CHARACTER (LEN=max_pathlen):: absfname
    CHARACTER (LEN=14) :: tmpchar

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=10), PARAMETER    :: modulename = 'getso2_crs'

    problems = .FALSE.

    absfname = TRIM(ADJUSTL(refdbdir)) // 'OMSAO_SO2_scia_fm.dat'

    IF (first) THEN
      tmpchar = ' ' 
      OPEN(UNIT = ozabs_unit, file=absfname, status='old')
      DO WHILE (tmpchar /= 'START OF TABLE') 
        READ (ozabs_unit, '(A14)') tmpchar
      ENDDO

      READ (ozabs_unit, *) nline, minw, maxw, normc
      IF (nline > maxline) THEN
        WRITE(www_lun, *) 'Need to increase parameter maxline!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      IF (minw > lamda(1) .OR. maxw < lamda(nlsav)) THEN
        WRITE(www_lun, *) 'Need to increase reference cross section range!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      READ (ozabs_unit, *) tdepend, nt, slitconv
      IF (nt > maxt) THEN
        WRITE(www_lun, *) 'Need to increase parameter maxt!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      IF (tdepend .AND. nt > 3) THEN
        WRITE(www_lun, *) 'High-order (>2) T-dependent cross section not implemented!!!'
        problems = .TRUE.
        CLOSE (ozabs_unit)
        RETURN
      ENDIF

      READ (ozabs_unit, *) ts(1:nt)
      !Allow extrapolation over temperature here
      !IF ((nt > 1 .AND. .NOT. tdepend ) .AND. (MINVAL(ts(1:nt)) > &
      !     MINVAL(tsgrid) .OR.  MAXVAL(ts(1:nt)) < MAXVAL(tsgrid))) THEN
      !   WRITE(www_lun, *) 'Temperature range for cross section not enough!!!'
      !   problems = .TRUE.; CLOSE (ozabs_unit); RETURN
      !ENDIF

      j = 1
      DO i = 1, nline
        READ(ozabs_unit, *) refwavs(j), refabs0(1:nt, j)
        IF (refwavs(j) > lamda(1) - 5.0 .AND. refwavs(j) < &
             lamda(nlsav) + 5.0)  j = j + 1
      ENDDO
      nline = j - 1
      CLOSE (ozabs_unit)
      first = .FALSE.
    ENDIF

    IF (so2crs_convl) THEN   
      refabs(1:nt, 1:nline) = refabs0(1:nt, 1:nline)

      IF (slitconv) THEN
        ni0 = n_refspec_pts(1)
        scalex = 0.1  ! dummy number, not used any more
        DO i = 1, nt 
          CALL CORRECT_I0EFFECT(refwavs(1:nline), refabs(i, 1:nline), nline, refspec_orig_data(1, 1:ni0, 1), &
               refspec_orig_data(1, 1:ni0, 2), ni0, scalex, get_lresi0, errstat, lowresi0(1:nline))
          IF ( errstat /= 0 ) THEN
            WRITE(www_lun, *) 'Error in Correct I0 Effect!!!'
            problems = .TRUE.
            RETURN
          ENDIF
        ENDDO
      ENDIF
      so2crs_convl = .FALSE.
    ENDIF

    ! Convolution and interpolation
    DO i = 1, nt
      CALL BSPLINE(refwavs(1:nline), refabs(i, 1:nline), nline, lamda, &
           savabs(i, :), nlsav, errstat)
      IF (errstat < 0) THEN
        WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
        problems = .TRUE.
        RETURN
      ENDIF

      IF (do_bandavg) THEN
        CALL avg_band_effozcrs(lamda, savabs(i, :), nlsav, ntemp, errstat)

        IF ( errstat /= 0 .OR. ntemp /= nlamda) THEN
          WRITE(www_lun, *) 'SO2 Spectra Averaging Error: ', nlsav, nlamda, ntemp
          problems = .TRUE.
          RETURN
        ENDIF
        tmpabs(i, :) = savabs(i, 1:nlamda)
      ELSE
        tmpabs(i, :) = savabs(i, 1:nlamda)
      ENDIF

    ENDDO

    ! Interpolate with respect to temperature
    IF (tdepend .AND. nt == 3) THEN       ! qudratic T dependent coefficients
      DO i = 1, nlayers
        thet = tsgrid(i) - zerok
        abscrs(:, i) = (tmpabs(1, :) + tmpabs(2, :) * thet &
             + tmpabs(3, :) * thet * thet )
      ENDDO
    ELSE IF (tdepend .AND. nt == 2) THEN   ! linear T dependent coefficients
      DO i = 1, nlayers
        thet = tsgrid(i) - zerok
        abscrs(:, i) = (tmpabs(1, :) + tmpabs(2, :) * thet)
      ENDDO
    ELSE IF (nt == 1)  THEN           ! only 1 T
      DO i = 1, nlayers
        abscrs(:, i) = tmpabs(1, :)
      ENDDO
    ELSE IF (nt == 2) THEN            ! have 2 T values
      DO i = 1, nlayers
        frac = 1.0 - (tsgrid(i) - ts(1)) / (ts(2) - ts(1))
        abscrs(:, i) = (frac * tmpabs(1, :) + (1.0 - frac) * tmpabs(2, :))
      END DO
    ELSE IF (.NOT. tdepend .AND. nt >= 3) THEN ! have more than n T     
      DO i = 1, nlamda           ! interpolate over T   
        CALL INTERPOL(ts(1:nt), tmpabs(1:nt,i), nt, tsgrid, abscrs(i, :), nlayers, errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
          problems = .TRUE.
          RETURN
        ENDIF
      ENDDO
    ELSE 
      WRITE(www_lun, *) modulename, ': Such type of SO2 cross sections not implemented' 
      problems = .TRUE.
      RETURN
    ENDIF

    abscrs = abscrs * normc

    RETURN  
  END SUBROUTINE GETSO2_CRS





  SUBROUTINE GET_ALB_OZCRS_RAY (nz, ts, ngas, abscrs, raycof, depol, problems)

    USE OMSAO_precision_module  
    USE OMSAO_parameters_module,ONLY : zerok, maxchlen, max_spec_pts
    USE OMSAO_variables_module, ONLY : refdbdir
    USE ozprof_data_module,     ONLY : ozcrs_alb_fname, ozabs_unit, pos_alb, &
         toms_fwhm
    USE OMSAO_errstat_module
    use m_lidort_prof_prep, only: get_all_raycof_depol

    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN)                              :: nz, ngas
    REAL (KIND=dp), INTENT(IN),  DIMENSION(nz)       :: ts
    REAL (KIND=dp), INTENT(OUT)                      :: raycof, depol
    REAL (KIND=dp), INTENT(OUT), DIMENSION(ngas, nz) :: abscrs
    LOGICAL, INTENT(OUT)                             :: problems

    ! Local variable
    INTEGER                                          :: nline, nw, i
    REAL (KIND=dp)                                   :: fwav, lwav, swav, ewav
    REAL (KIND=dp), DIMENSION(11)                    :: temp
    REAL (KIND=dp), DIMENSION(max_spec_pts)          :: waves, sol, &
         rays, dpols!, weights
    !Make weights allocatable to avoid ridiculously large array
    real (kind=dp), dimension(:), allocatable :: weights
    REAL (KIND=dp), DIMENSION(3, max_spec_pts)       :: ozcrs  
    REAL (KIND=dp), DIMENSION(6, max_spec_pts)       :: gcrs 
    CHARACTER (len=maxchlen)                         :: ozcrs_fname

    ! Saved variables
    LOGICAL,                      SAVE :: first = .TRUE.
    REAL (KIND=dp), DIMENSION(3), SAVE :: cozcrs
    REAL (KIND=dp), DIMENSION(6), SAVE :: cgcrs
    REAL (KIND=dp),               SAVE :: craycof, cdepol

    problems = .FALSE.
    IF (first) THEN
      ozcrs_fname = TRIM(ADJUSTL(refdbdir)) // '/' // TRIM(ADJUSTL(ozcrs_alb_fname))
      OPEN(UNIT = ozabs_unit, file=ozcrs_fname, status='old')     
      DO i = 1, 5
        READ(ozabs_unit, *)
      ENDDO
      READ(ozabs_unit, *) nline, fwav, lwav
      READ(ozabs_unit, *)
      swav = pos_alb - toms_fwhm
      ewav = pos_alb + toms_fwhm
      IF (swav < fwav .OR. ewav > lwav) THEN
        WRITE(www_lun, *) 'Ozone cross section does not cover region for determining fcld!!!'
        problems = .TRUE.
        CLOSE(ozabs_unit)
        RETURN
      ENDIF

      nw = 0
      DO i = 1, nline
        READ (ozabs_unit, *) temp
        IF (temp(1) >= swav .AND. temp(1) <= ewav) THEN
          nw = nw + 1
          waves(nw) = temp(1)
          ozcrs(1:3, nw) = temp(2:4)
          gcrs(1:6, nw) = temp(5:10)
          sol(nw) = temp(11)
        ELSE IF (temp(1) > ewav) THEN
          EXIT
        ENDIF
      ENDDO
      CLOSE (ozabs_unit)

      ! Compute weights
      allocate(weights(nw))
      weights(1:nw) = (1.0 - ABS(waves(1:nw) - pos_alb) / &
           toms_fwhm) * sol(1:nw)
      weights = weights(1:nw) /  SUM(weights(1:nw))

      DO i = 1, 3
        cozcrs(i) = SUM(ozcrs(i, 1:nw) *  weights(1:nw))
      ENDDO

      DO i = 1, 6
        cgcrs(i) = SUM(gcrs(i, 1:nw) *  weights(1:nw))
      ENDDO

      CALL GET_ALL_RAYCOF_DEPOL(nw, waves(1:nw), rays(1:nw), dpols(1:nw))
      craycof = SUM(rays(1:nw)   * weights(1:nw)) 
      cdepol  = SUM(dpols(1:nw)  * weights(1:nw)) 

      first = .FALSE.
    ENDIF

    abscrs(1, :) = cozcrs(1) + (ts - zerok) * cozcrs(2) + &
         (ts - zerok) ** 2.0 * cozcrs(3)
    DO i = 1, 6
      abscrs(i+1, :) = cgcrs(i)
    ENDDO
    raycof = craycof
    depol = cdepol

    if(allocated(weights)) deallocate(weights)

    RETURN
  END SUBROUTINE GET_ALB_OZCRS_RAY


  !   Unused ?
  !
  ! Compute slant optical thickness using Chapman function in LIDORT
  SUBROUTINE GET_SLANT_TAU(nz, zs, tauin, sza, tauout)
    USE OMSAO_parameters_module, ONLY  : maxchlen, rearth
    USE OMSAO_precision_module
    use m_vlidort90_include
    !    use vlidort_miscsetups_src, only: vlidort_chapman

    IMPLICIT NONE

    !===============================  Define Variables ===========================
    ! Include files of dimensions and numbers
!    INCLUDE 'VLIDORT90.PARS'
    INCLUDE 'VLIDORT_INPUTS90.VARS'
    INCLUDE 'VLIDORT_SETUPS90.VARS'

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN)                          :: nz
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(IN)  :: zs, tauin
    REAL (KIND=dp), INTENT(IN)                   :: sza
    REAL (KIND=dp), DIMENSION(0:nz), INTENT(out) :: tauout

    ! =======================
    ! Local variables
    ! ======================= 
    INTEGER                   :: i, j
    LOGICAL                   :: fail
    CHARACTER (len=maxchlen)  :: message, trace

    n_szangles = 1
    szangles(1) = sza
    IF (sza >= 90.0 .OR. sza < 0) THEN
      print *, 'GET_SLANT_TAU: SZA is >= 90 or < 0!!!'
      STOP 1
    ENDIF

    nlayers = nz
    taugrid_input(0:nz) = tauin
    height_grid(0:nz) = zs
    earth_radius = rearth
    IF (nz > maxlayers) THEN
      print *, 'LIDORT_PROF_ENV: # of layers exceeded allowed !!!'
      STOP 1
    ENDIF

    CALL VLIDORT_CHAPMAN(fail, message, trace)
    tauout = 0.0
    DO i = 1, nz
      DO j = 1, i
        tauout(i) = tauout(i) + deltau_slant(i, j, 1)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE GET_SLANT_TAU

  ! Obtain minor trace gas weighting functions from ozwf
  ! Note: it is not the exact wf but negative of the WF divided by radiances
  SUBROUTINE GET_TRACEGAS_WF (ozwf, ozabs, so2crs, use_so2dtcrs, rad, nw, &
       nz, nz1, &
       ozs, waves, do_so2zwf, so2zwf)
    USE OMSAO_precision_module
    USE OMSAO_parameters_module,ONLY : du2mol
    USE OMSAO_indices_module,   ONLY : so2_idx, so2v_idx!, &
    !ring_idx, ring1_idx, no2_t1_idx, hcho_idx, bro_idx
    USE OMSAO_variables_module, ONLY : refidx, database, &
         database_save, refspec_norm!, mask_fitvar_rad
    USE ozprof_data_module,     ONLY : mgasprof, fgasidxs, tracegas, &
         ngas, gasidxs, fgassidxs, so2valts, so2vprofn1p1!, nlay, fps, fzs
    IMPLICIT NONE

    ! Input/Output variables
    INTEGER, INTENT(IN)                           :: nw, nz, nz1
    REAL (KIND=dp), DIMENSION(nw, nz), INTENT(IN) :: ozwf, ozabs, so2crs
    REAL (KIND=dp), DIMENSION(nw),     INTENT(IN) :: rad, waves
    REAL (KIND=dp), DIMENSION(nw),    INTENT(OUT) :: so2zwf
    REAL (KIND=dp), DIMENSION(nz),     INTENT(IN) :: ozs
    LOGICAL,                           INTENT(IN) :: use_so2dtcrs, do_so2zwf

    ! Local variables
    INTEGER                                       :: i, j, nk!, k, lidx, fidx
    REAL (KIND=dp)                                :: tmp
    REAL (KIND=dp), DIMENSION(nw, nz)             :: amf
    REAL (KIND=dp), DIMENSION(ngas, nw)           :: tamf
    REAL (KIND=dp)                                :: avcd

    ! Obtain AMF  each wavelength and at each layer
    IF (ANY(fgasidxs > 0)) THEN
      DO i = 1, nz1
        amf(:, i) = -ozwf(:, i) / rad / ozabs(:, i) / du2mol
      ENDDO
    ENDIF

    ! Replace cross sections with weighting functions
    DO i = 1, ngas 
      IF (fgasidxs(i) > 0) THEN
        avcd = mgasprof(i, nz+1)

        IF ((gasidxs(i) /= so2_idx .AND. gasidxs(i) /= so2v_idx) .OR. .NOT. use_so2dtcrs) THEN
          DO j = 1, nw
            tamf(i, j) = SUM(amf(j, 1:nz1) * mgasprof(i, 1:nz1)) / avcd          
          ENDDO

          IF (fgassidxs(i) > 0) THEN
            database(gasidxs(i), refidx(1:nw)) = database(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
          ELSE
            database(gasidxs(i), refidx(1:nw)) = database_save(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
          ENDIF

          tracegas(i, 7) = 0.0
          nk = 0
          DO j = 1, nw
            IF (database_save(gasidxs(i), refidx(j)) > 0.) THEN
              tracegas(i, 7) = tracegas(i, 7) + tamf(i, j)
              nk = nk + 1
            ENDIF
          ENDDO
          tracegas(i, 7) = tracegas(i, 7) / nk
        ELSE
          tmp = avcd * refspec_norm(gasidxs(i))
          DO j = 1, nw
            tamf(i, j) = SUM(amf(j, 1:nz1) * mgasprof(i, 1:nz1) * so2crs(j, 1:nz1) ) / tmp
          ENDDO
          database(gasidxs(i), refidx(1:nw)) = tamf(i, 1:nw) 

          tracegas(i, 7) = 0.0
          nk = 0
          DO j = 1, nw
            IF (database_save(gasidxs(i), refidx(j)) > 0.) THEN
              tracegas(i, 7) = tracegas(i, 7) + tamf(i, j) / database_save(gasidxs(i), refidx(j)) 
              nk = nk + 1
            ENDIF
          ENDDO
          tracegas(i, 7) = tracegas(i, 7) / nk
        ENDIF

        ! Calculate weighting function for SO2V plume height, verified with finite difference
        IF (do_so2zwf .AND. gasidxs(i) == so2v_idx) THEN
          tmp = du2mol * refspec_norm(gasidxs(i)) * (so2valts(1)-so2valts(-1)) * avcd / tracegas(i, 4)

          IF (use_so2dtcrs) THEN
            DO j = 1, nw
              so2zwf(j) = SUM(ozwf(j, 1:nz1) / ozabs(j, 1:nz1) * so2crs(j, 1:nz1) &
                   * (so2vprofn1p1(1:nz1, 2)-so2vprofn1p1(1:nz1, 1)) ) / tmp
            ENDDO
          ELSE
            tmp =tmp / refspec_norm(gasidxs(i))
            DO j = 1, nw
              so2zwf(j) = SUM(ozwf(j, 1:nz1) / ozabs(j, 1:nz1) * database_save(gasidxs(i), refidx(j)) &
                   * (so2vprofn1p1(1:nz1, 2)-so2vprofn1p1(1:nz1, 1)) ) / tmp
            ENDDO
          ENDIF
        ENDIF

      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE GET_TRACEGAS_WF

  SUBROUTINE radwf_interpol(nw, nz, nctp, ncbp, nsprs, faerlvl, do_radcals, &
       do_fozwf, do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, &
       do_cfracwf, wave, abscrs, ozs, rad, fozwf, albwf, cfracwf, faerwf, &
       faerswf, fcodwf, fsprswf, errstat)

    USE OMSAO_precision_module
    USE OMSAO_variables_module, ONLY  : numwin, nradpix, band_selectors 
    !USE ozprof_data_module,     ONLY  : nup2p 
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN)                              :: nw, nz, nctp, ncbp, faerlvl, nsprs
    INTEGER, INTENT(OUT)                             :: errstat                                                   
    LOGICAL, INTENT(IN)                              :: do_fozwf, do_albwf, &
         do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf
    LOGICAL, DIMENSION(nw), INTENT(IN)               :: do_radcals

    REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
    REAL (KIND=dp), DIMENSION(nw, nz), INTENT(IN)    :: abscrs
    REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fozwf, faerwf, faerswf, fcodwf, fsprswf
    REAL (KIND=dp), DIMENSION(nw),     INTENT(IN)    :: wave
    REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: rad, albwf, cfracwf

    ! Local variables
    INTEGER :: i, j, k,  iw, nd, nd1, nud, p1, p2, dp1, dp2, fidx, lidx
    INTEGER, DIMENSION(nw)                           :: didxs, uidxs!, radcals
    REAL (KIND=dp)                                   :: toz, df
    REAL (KIND=dp), DIMENSION(nw)                    :: effcrs, a, b, crs1, crs2, tmp1, tmp2, wav1, wav2

    !REAL (KIND=dp), DIMENSION(nw, nz) :: fozwf1
    !REAL (KIND=dp), DIMENSION(nw, nl) :: ozwf1, tmpwf1
    !REAL (KIND=dp), DIMENSION(nw)     :: rad1, albwf1, shiwf1

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=14), PARAMETER :: modulename = 'radwf_interpol'


    errstat = pge_errstat_ok

    !fozwf1 = fozwf
    !ozwf1 = ozwf
    !albwf1 = albwf
    !rad1 = rad
    !shiwf1 = shiwf
    !
    !radcals = 0
    !WHERE (do_radcals)
    !   radcals = 1
    !ENDWHERE

    ! Calculate effective cross section or first derivative wrt wavelength
    effcrs = 0.0
    toz = SUM(ozs)
    DO i = 1, nw - 1
      effcrs(i) = SUM(abscrs(i, :) * ozs) / toz
    ENDDO

    ! For ozone weighting functions (fine grids), divide by radiance and change its negative sign
    IF (do_fozwf) THEN
      DO i = 1, nz
        WHERE(do_radcals)
          fozwf(:, i) = -fozwf(:, i) / rad
        ENDWHERE
      ENDDO
    ENDIF

    ! For albedo weighting functions, divide by radiance 
    IF (do_albwf) THEN
      WHERE(do_radcals)
        albwf = albwf / rad
      ENDWHERE
    ENDIF

    IF (do_cfracwf) THEN
      WHERE(do_radcals)
        cfracwf = cfracwf / rad
      ENDWHERE
    ENDIF

    IF (do_faerwf) THEN
      DO i = faerlvl, nz
        WHERE(do_radcals)
          faerwf(:, i) = faerwf(:, i) / rad
        ENDWHERE
      ENDDO
    ENDIF

    IF (do_faerswf) THEN
      DO i = faerlvl, nz
        WHERE(do_radcals)
          faerswf(:, i) = faerswf(:, i) / rad
        ENDWHERE
      ENDDO
    ENDIF

    IF (do_codwf) THEN
      DO i = nctp, ncbp
        WHERE(do_radcals)
          fcodwf(:, i) = fcodwf(:, i) / rad
        ENDWHERE
      ENDDO
    ENDIF

    IF (do_sprswf) THEN
      DO i = nsprs, nz
        WHERE(do_radcals)
          fsprswf(:, i) = fsprswf(:, i) / rad
        ENDWHERE
      ENDDO
    ENDIF

    ! Take the logarithm of radiances
    WHERE (do_radcals)
      rad = LOG(rad)
    ENDWHERE

    fidx = 1
    DO iw = 1, numwin
      lidx = fidx + nradpix(iw) - 1

      IF (band_selectors(iw) == 1) THEN ! channel 1, use cubic-spline, ozone cross is decreasing

        nd = 0
        nud = 0
        DO i = fidx, lidx
          IF ( do_radcals(i) ) THEN
            nd = nd + 1
            didxs(nd)  = i
          ELSE
            nud = nud + 1
            uidxs(nud) = i
          ENDIF
        ENDDO

        ! for radiances (lnI vs. ozcrs)
        k = MIN(uidxs(1)-2, 1)
        nd1 = nd - k + 1
        crs1(1:nd1)  = effcrs(didxs(k:nd))
        tmp1(1:nd1)  = rad(didxs(k:nd))
        crs2(1:nud) = effcrs(uidxs(1:nud))

        CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
          errstat = pge_errstat_error
          RETURN
        ENDIF

        rad(uidxs(1:nud)) =  tmp2(1:nud)

        ! for albedo weighting funtions: ln(albwf/rad) vs ozcrs
        IF (do_albwf .AND. MAXVAL(albwf) > 0.0) THEN
          tmp1(1:nd1) = LOG(albwf(didxs(k:nd)))

          CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
          IF (errstat < 0) THEN
            WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
            errstat = pge_errstat_error
            RETURN
          ENDIF

          albwf(uidxs(1:nud)) =  EXP(tmp2(1:nud))
        ENDIF

        IF (do_cfracwf .AND. MAXVAL(cfracwf) > 0.0) THEN
          tmp1(1:nd1) = LOG(cfracwf(didxs(k:nd)))

          CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
          IF (errstat < 0) THEN
            WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
            errstat = pge_errstat_error
            RETURN
          ENDIF

          cfracwf(uidxs(1:nud)) =  EXP(tmp2(1:nud))
        ENDIF

        ! for ozone weighting function (fine grids) ln(ozwf/rad) vs ozcrs
        IF (do_fozwf) THEN
          DO j = 1, nz 
            crs1(1:nd1)  = abscrs(didxs(k:nd), j)
            tmp1(1:nd1) = LOG(fozwf(didxs(k:nd), j))
            crs2(1:nud) = abscrs(uidxs(1:nud), j)

            CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
              errstat = pge_errstat_error
              RETURN
            ENDIF

            fozwf(uidxs(1:nud), j) = EXP(tmp2(1:nud))
          ENDDO
        ENDIF

        ! Linear interpolation over wavelength
        IF (do_faerwf) THEN
          wav1(1:nd1) = wave(didxs(k:nd))
          wav2(1:nud) = wave(uidxs(1:nud))
          DO j = faerlvl, nz 
            tmp1(1:nd1) = faerwf(didxs(k:nd), j)

            CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
              errstat = pge_errstat_error
              RETURN
            ENDIF

            faerwf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
        ENDIF

        IF (do_faerswf) THEN
          wav1(1:nd1) = wave(didxs(k:nd))
          wav2(1:nud) = wave(uidxs(1:nud))
          DO j = faerlvl, nz 
            tmp1(1:nd1) = faerswf(didxs(k:nd), j)

            CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
              errstat = pge_errstat_error
              RETURN
            ENDIF

            faerswf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
        ENDIF

        IF (do_codwf) THEN
          wav1(1:nd1) = wave(didxs(k:nd))
          wav2(1:nud) = wave(uidxs(1:nud))
          DO j = nctp, ncbp 
            tmp1(1:nd1) = fcodwf(didxs(k:nd), j)

            CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
              errstat = pge_errstat_error
              RETURN
            ENDIF

            fcodwf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
        ENDIF

        IF (do_sprswf) THEN
          wav1(1:nd1) = wave(didxs(k:nd))
          wav2(1:nud) = wave(uidxs(1:nud))
          DO j = nsprs, nz 
            tmp1(1:nd1) = fsprswf(didxs(k:nd), j)

            CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
            IF (errstat < 0) THEN
              WRITE(www_lun, *) modulename, ': BSPLINE error, errstat = ', errstat
              errstat = pge_errstat_error
              RETURN
            ENDIF

            fsprswf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
        ENDIF

      ELSE IF (band_selectors(iw) == 2) THEN ! Channel 2, linear interpolation vs cross section

        ! p1,  p2:  points without radiative calculation
        ! dp1, dp2: points with radiative calculation
        dp1 = fidx
        i = dp1 + 1
        DO WHILE (i <= lidx)
          IF ( do_radcals(i) ) THEN
            dp2 = i

            IF (dp2 - dp1 > 1) THEN  ! there are points with no rad/wf in between
              ! get rad/wf in between using rad/wf at dp1:dp2
              df = effcrs(dp2) - effcrs(dp1)
              p1 = dp1 + 1
              p2 = dp2 - 1
              a(p1 : p2) = ( effcrs(dp2) - effcrs(p1 : p2) ) / df
              b(p1 : p2) = 1.0 - a(p1 : p2)
              rad(p1 : p2) = a(p1 : p2) * rad(dp1) + b(p1 : p2) * rad(dp2)
              IF (do_albwf) albwf(p1 : p2) = a(p1 : p2) * albwf(dp1) + b(p1 : p2) * albwf(dp2)
              IF (do_cfracwf) cfracwf(p1 : p2) = a(p1 : p2) * cfracwf(dp1) + b(p1 : p2) * cfracwf(dp2)

              IF (do_fozwf) THEN
                DO j = 1, nz
                  df = abscrs(dp2, j) - abscrs(dp1, j)
                  a(p1 : p2) = ( abscrs(dp2, j) - abscrs(p1 : p2, j) ) / df
                  b(p1 : p2) = 1.0 - a(p1 : p2)
                  fozwf(p1 : p2, j) = a(p1 : p2) * fozwf(dp1, j) + b(p1 : p2) * fozwf(dp2, j)
                ENDDO
              ENDIF

              IF (do_faerwf) THEN
                DO j = faerlvl, nz
                  df = wave(dp2) - wave(dp1)
                  a(p1 : p2) = ( wave(dp2) - wave(p1 : p2) ) / df
                  b(p1 : p2) = 1.0 - a(p1 : p2)
                  faerwf(p1 : p2, j) = a(p1 : p2) * faerwf(dp1, j) + b(p1 : p2) * faerwf(dp2, j)
                ENDDO
              ENDIF

              IF (do_faerswf) THEN
                DO j = faerlvl, nz
                  df = wave(dp2) - wave(dp1)
                  a(p1 : p2) = ( wave(dp2) - wave(p1 : p2) ) / df
                  b(p1 : p2) = 1.0 - a(p1 : p2)
                  faerswf(p1 : p2, j) = a(p1 : p2) * faerswf(dp1, j) + b(p1 : p2) * faerswf(dp2, j)
                ENDDO
              ENDIF

              IF (do_codwf) THEN
                DO j = nctp, ncbp
                  df = wave(dp2) - wave(dp1)
                  a(p1 : p2) = ( wave(dp2) - wave(p1 : p2) ) / df
                  b(p1 : p2) = 1.0 - a(p1 : p2)
                  fcodwf(p1 : p2, j) = a(p1 : p2) * fcodwf(dp1, j) + b(p1 : p2) * fcodwf(dp2, j)
                ENDDO
              ENDIF

              IF (do_sprswf) THEN
                DO j = nsprs, nz
                  df = wave(dp2) - wave(dp1)
                  a(p1 : p2) = ( wave(dp2) - wave(p1 : p2) ) / df
                  b(p1 : p2) = 1.0 - a(p1 : p2)
                  fsprswf(p1 : p2, j) = a(p1 : p2) * fsprswf(dp1, j) + b(p1 : p2) * fsprswf(dp2, j)
                ENDDO
              ENDIF

            ENDIF
            dp1 = dp2
          ENDIF
          i = i + 1
        ENDDO
      ELSE                                 ! Not implemented
        WRITE(www_lun, *) 'Interpolation for this channel is not implemented!!!'
        errstat = pge_errstat_error
      ENDIF
      fidx = lidx + 1
    ENDDO

    ! Convert radiances back
    rad = EXP(rad)

    !DO i = 1, nw
    !   WRITE(90, '(F10.4, I5, 2D14.6)') wave(i), radcals(i), rad1(i), rad(i)
    !ENDDO

    ! Convert ozone weighting functions back
    IF (do_fozwf) THEN
      DO i = 1, nz
        fozwf(:, i) = - fozwf(:, i) * rad
      ENDDO
    ENDIF

    IF (do_albwf) THEN
      albwf = albwf * rad
    ENDIF

    IF (do_cfracwf) THEN
      cfracwf = cfracwf * rad
    ENDIF

    IF (do_faerwf) THEN
      DO i = faerlvl, nz
        faerwf(:, i) = faerwf(:, i) * rad
      ENDDO
    ENDIF

    IF (do_faerswf) THEN
      DO i = faerlvl, nz
        faerswf(:, i) = faerswf(:, i) * rad
      ENDDO
    ENDIF

    IF (do_codwf) THEN
      DO i = nctp, ncbp
        fcodwf(:, i) = fcodwf(:, i) * rad
      ENDDO
    ENDIF

    IF (do_sprswf) THEN
      DO i = nsprs, nz
        fsprswf(:, i) = fsprswf(:, i) * rad
      ENDDO
    ENDIF

    !DO i = 1, nw
    !   WRITE(91, '(F10.4, I5, 2D14.6)') wave(i),  radcals(i), albwf1(i), albwf(i)
    !ENDDO

    !DO i = 1, nw
    !   WRITE(92, '(F10.4, I5, 200D14.6)') wave(i), radcals(i), fozwf1(i, 1:nz), &
    !        fozwf(i, 1:nz)
    !ENDDO  
    !print *, nz, nl, nw
    !!STOP

    RETURN

  END SUBROUTINE radwf_interpol

  SUBROUTINE polcorr_online(niter, which_polcorr, nw,  nz, nctp, ncbp, &
       nsprs, faerlvl, ncorr, polidxs, do_fozwf, do_albwf, do_faerwf, &
       do_faerswf, do_codwf, do_sprswf, do_fraywf, do_cfracwf, wave, &
       rad, prad, tabs, ozs, albwf, palbwf, fozwf, pfozwf, faerwf, &
       pfaerwf, faerswf, pfaerswf, fcodwf, &
       pfcodwf, fsprswf, pfsprswf, fraywf, pfraywf, cfracwf, pcfracwf)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module, ONLY : mflay!, maxlay
    USE OMSAO_variables_module,  ONLY : currloop
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN) :: nw, nz, ncorr, niter, which_polcorr, nctp, ncbp, &
         faerlvl, nsprs
    LOGICAL, INTENT(IN) :: do_fozwf, do_albwf, do_faerwf, do_faerswf, &
         do_codwf, do_sprswf, do_fraywf, do_cfracwf

    INTEGER, DIMENSION (ncorr),          INTENT (IN) :: polidxs
    REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fozwf, faerwf, &
         faerswf, fcodwf, fsprswf, fraywf
    REAL (KIND=dp), DIMENSION(nw, nz),    INTENT(IN) :: pfozwf, pfaerwf, &
         pfaerswf, pfcodwf, pfsprswf, tabs, pfraywf
    REAL (KIND=dp), DIMENSION(nw),        INTENT(IN) :: wave
    REAL (KIND=dp), DIMENSION(nz),        INTENT(IN) :: ozs
    REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: rad, albwf, cfracwf
    REAL (KIND=dp), DIMENSION(nw),        INTENT(IN) :: prad, palbwf, pcfracwf

    ! Local variables
    INTEGER, PARAMETER      :: mcorr = 20, mw = 500
    INTEGER                 :: i, j, fidx, lidx
    REAL (KIND=dp)          :: frac, frac1, frac2, tmprad
    REAL (KIND=dp), DIMENSION(nz)                    :: tmpcorr
    REAL (KIND=dp), DIMENSION(mcorr, mflay), SAVE    :: tauwf, ptauwf
    REAL (KIND=dp), DIMENSION(mcorr),      SAVE  :: dfrad, dfalbwf, dfcfracwf
    REAL (KIND=dp), DIMENSION(mcorr, mflay), SAVE    :: dftauwf
    REAL (KIND=dp), DIMENSION(mcorr, mflay), SAVE :: dffozwf, dffaerwf, &
         dffaerswf, dffcodwf, dffsprswf, dffraywf
    REAL (KIND=dp), DIMENSION(mw, mflay), SAVE       :: fozwf_sav
    REAL (KIND=dp), DIMENSION(mw), SAVE              :: albwf_sav, cfracwf_sav
    LOGICAL                                          :: newcorr = .TRUE.

    ! Compute differences at positions where corrections are explicitly calcualted
    IF (niter == 0 .OR. which_polcorr == 3 .OR. which_polcorr == 5 ) THEN
      dfrad(1:ncorr) = 0.0D0

      ! Corrected in different way using DlnI/Dtau: 08/28/2008
      IF (newcorr) THEN
        WHERE ( prad(polidxs) /= 0.0 )
          dfrad(1:ncorr) = LOG(prad(polidxs))- LOG(rad(polidxs))
        ENDWHERE
      ELSE
        WHERE ( prad(polidxs) /= 0.0 )
          dfrad(1:ncorr) = (rad(polidxs) - prad(polidxs)) / prad(polidxs)
        ENDWHERE
      ENDIF

      IF ( which_polcorr /= 5 .OR. (which_polcorr == 5 .AND. (niter == 1 .OR. (niter == 0 .AND. currloop == 0))) ) THEN

        IF (do_albwf)   dfalbwf(1:ncorr) = 0.D0
        dftauwf(1:ncorr, 1:nz) = 0.D0
        IF (do_fozwf)   THEN
          dffozwf(1:ncorr, 1:nz) = 0.D0
          tauwf  (1:ncorr, 1:nz) = 0.D0
          ptauwf (1:ncorr, 1:nz) = 0.D0
        ENDIF
        IF (do_cfracwf) dfcfracwf(1:ncorr) = 0.D0
        IF (do_faerwf)  dffaerwf(1:ncorr, faerlvl:nz)  = 0.D0
        IF (do_faerswf) dffaerswf(1:ncorr, faerlvl:nz) = 0.D0
        IF (do_codwf)   dffcodwf(1:ncorr, nctp:ncbp)   = 0.D0
        IF (do_sprswf)  dffsprswf(1:ncorr, nsprs:nz)   = 0.D0
        IF (do_fraywf)  dffraywf(1:ncorr, 1:nz) = 0.D0

        !print *, polidxs
        !print *, wave(polidxs)
        !print *, dfrad(1:ncorr)

        IF ( do_albwf) THEN
          WHERE (palbwf(polidxs) /= 0.0)
            dfalbwf(1:ncorr) = (albwf(polidxs) - palbwf(polidxs)) / palbwf(polidxs)
          ENDWHERE
        ENDIF

        IF ( do_cfracwf) THEN
          WHERE (pcfracwf(polidxs) /= 0.0)
            dfcfracwf(1:ncorr) = (cfracwf(polidxs) - pcfracwf(polidxs)) / pcfracwf(polidxs)
          ENDWHERE
        ENDIF

        IF ( do_fozwf ) THEN
          WHERE( pfozwf(polidxs, 1:nz) /= 0.0)
            dffozwf(1:ncorr, 1:nz) = (fozwf(polidxs, 1:nz) - pfozwf(polidxs, 1:nz)) / pfozwf(polidxs, 1:nz)
          ENDWHERE

          DO i = 1, ncorr
            ptauwf(i, 1:nz) = pfozwf(polidxs(i), 1:nz) * ozs / tabs(polidxs(i), 1:nz) / prad(polidxs(i))
            tauwf (i, 1:nz) = fozwf (polidxs(i), 1:nz) * ozs / tabs(polidxs(i), 1:nz) / rad (polidxs(i))
          ENDDO
          dftauwf(1:ncorr, 1:nz) = ptauwf(1:ncorr, 1:nz) - tauwf(1:ncorr, 1:nz)
          !        WRITE(www_lun, *) ' ***'
          !        print *, prad(polidxs(4)), rad(polidxs(4))
          !        WRITE(www_lun, '(10D12.4)') dftauwf(4, 1:nz)      
          !        WRITE(www_lun, '(10D12.4)') ptauwf(4, 1:nz)
          !        WRITE(www_lun, '(10D12.4)') tauwf(4, 1:nz)
          !        WRITE(www_lun, '(10D12.4)') pfozwf(polidxs(4), 1:nz)
          !        WRITE(www_lun, '(10D12.4)') fozwf(polidxs(4), 1:nz)
          !        WRITE(www_lun, *) ' ***'
        ENDIF

        IF ( do_faerwf ) THEN
          WHERE( pfaerwf(polidxs, faerlvl:nz) /= 0.0 )
            dffaerwf(1:ncorr, faerlvl:nz) = (faerwf(polidxs, faerlvl:nz) &
                 - pfaerwf(polidxs, faerlvl:nz)) / pfaerwf(polidxs, faerlvl:nz)
          ENDWHERE
        ENDIF

        IF ( do_faerswf ) THEN
          WHERE ( pfaerswf(polidxs, faerlvl:nz) /= 0.0 )
            dffaerswf(1:ncorr, faerlvl:nz) = (faerswf(polidxs, faerlvl:nz) &
                 - pfaerswf(polidxs, faerlvl:nz)) / pfaerswf(polidxs, faerlvl:nz)
          ENDWHERE
        ENDIF

        IF ( do_codwf ) THEN
          WHERE ( pfcodwf( polidxs, nctp:ncbp ) /= 0.0 )
            dffcodwf(1:ncorr, nctp:ncbp) = (fcodwf(polidxs, nctp:ncbp) - pfcodwf(polidxs, nctp:ncbp)) &
                 / pfcodwf(polidxs, nctp:ncbp)
          ENDWHERE
        ENDIF

        IF ( do_sprswf ) THEN
          WHERE ( pfsprswf( polidxs, nsprs:nz ) /= 0.0 )
            dffsprswf(1:ncorr, nsprs:nz) = (fsprswf(polidxs, nsprs:nz) - pfsprswf(polidxs, nsprs:nz)) &
                 / pfsprswf(polidxs, nsprs:nz)
          ENDWHERE
        ENDIF

        IF ( do_fraywf ) THEN
          WHERE( pfraywf(polidxs, 1:nz) /= 0.0)
            dffraywf(1:ncorr, 1:nz) = (fraywf(polidxs, 1:nz) - pfraywf(polidxs, 1:nz)) / pfraywf(polidxs, 1:nz)
          ENDWHERE
        ENDIF

      ENDIF
    ENDIF

    ! Special correction for wavelengths before polidxs(1)
    fidx = 1
    lidx=polidxs(1)
    DO j = fidx, lidx
      IF (lidx /= 1) THEN
        frac = (wave(j) - wave(1) ) / ( wave(lidx) - wave(1) )
      ELSE
        frac = 1.0
      ENDIF

      IF (newcorr) THEN
        tmprad = frac * (SUM(dftauwf(1, 1:nz) * (tabs(j, 1:nz)-tabs(polidxs(1), 1:nz))) + dfrad(1))
        rad(j) = rad(j) * EXP(tmprad)
      ELSE
        rad(j) = rad(j) / (1.0 + frac * dfrad(1))
      ENDIF

      IF ( which_polcorr /= 5 .OR. (which_polcorr == 5 .AND. (niter == 1 .OR. (niter == 0 .AND. currloop == 0))) ) THEN

        IF ( do_albwf)  THEN
          tmpcorr(1) = (1.0 + frac * dfalbwf(1)) 
          IF (tmpcorr(1) /= 0.0) albwf(j) = albwf(j) / tmpcorr(1)
        ENDIF
        IF ( do_cfracwf)  THEN
          tmpcorr(1) = (1.0 + frac * dfcfracwf(1)) 
          IF (tmpcorr(1) /= 0.0) cfracwf(j) = cfracwf(j) / tmpcorr(1)
        ENDIF
        IF ( do_fozwf  ) THEN
          tmpcorr(1:nz) = (1.0 + frac * dffozwf(1, 1:nz)) 
          WHERE (tmpcorr(1:nz) /= 0.0)
            fozwf(j, :) = fozwf(j, :) / tmpcorr(1:nz)
          ENDWHERE
        ENDIF
        IF ( do_faerwf ) THEN
          tmpcorr(faerlvl:nz) = (1.0 + frac * dffaerwf(1, faerlvl:nz))
          WHERE (tmpcorr(faerlvl:nz) /= 0.0) 
            faerwf(j, faerlvl:nz)   = faerwf(j, faerlvl:nz)  / tmpcorr(faerlvl:nz)
          ENDWHERE
        ENDIF
        IF ( do_faerswf ) THEN
          tmpcorr(faerlvl:nz) = (1.0 + frac * dffaerswf(1, faerlvl:nz))
          WHERE (tmpcorr(faerlvl:nz) /= 0.0) 
            faerswf(j, faerlvl:nz)   = faerswf(j, faerlvl:nz)  / tmpcorr(faerlvl:nz)
          ENDWHERE
        ENDIF
        IF ( do_codwf )  THEN
          tmpcorr(nctp:ncbp) = (1.0 + frac * dffcodwf(1, nctp:ncbp))
          WHERE (tmpcorr(nctp:ncbp) /= 0.0)
            fcodwf(j, nctp:ncbp) = fcodwf(j, nctp:ncbp) / tmpcorr(nctp:ncbp)
          ENDWHERE
        ENDIF

        IF ( do_sprswf )  THEN
          tmpcorr(nsprs:nz) = (1.0 + frac * dffsprswf(1, nsprs:nz))
          WHERE (tmpcorr(nsprs:nz) /= 0.0)
            fsprswf(j, nsprs:nz) = fsprswf(j, nsprs:nz) / tmpcorr(nsprs:nz)
          ENDWHERE
        ENDIF

        IF ( do_fraywf  ) THEN
          tmpcorr(1:nz) = (1.0 + frac * dffraywf(1, 1:nz)) 
          WHERE (tmpcorr(1:nz) /= 0.0)
            fraywf(j, :) = fraywf(j, :) / tmpcorr(1:nz)
          ENDWHERE
        ENDIF

      ENDIF
    ENDDO
    !print *, nw
    !print *, fidx, lidx

    fidx = lidx + 1
    DO i = 2, ncorr
      lidx = polidxs(i)
      DO j = fidx, lidx
        frac2 = (wave(j) - wave(fidx-1)) / (wave(lidx) - wave(fidx-1))
        frac1 = 1.0 - frac2

        IF (newcorr) THEN
          !IF (j > 102 .AND. j < 104) THEN
          !   print *, j, rad(j), frac1, frac2, dfrad(i-1), dfrad(i), &
          !   SUM(dftauwf(i-1, 1:nz) * (tabs(j, 1:nz)-tabs(polidxs(i-1), 1:nz))), &
          !   SUM(dftauwf(i, 1:nz) * (tabs(j, 1:nz)-tabs(polidxs(i), 1:nz)))
          !   write(*, '(10D12.4)') dftauwf(i-1, 1:nz)
          !   write(*, '(10D12.4)') tabs(j, 1:nz)
          !   write(*, '(10D12.4)') tabs(polidxs(i-1), 1:nz)
          !   write(*, '(10D12.4)') dftauwf(i, 1:nz)
          !   write(*, '(10D12.4)') tabs(polidxs(i), 1:nz)
          !ENDIF
          tmprad = frac1 * (SUM(dftauwf(i-1, 1:nz) * (tabs(j, 1:nz)-tabs(polidxs(i-1), 1:nz))) + dfrad(i-1)) + &
               frac2 * (SUM(dftauwf(i, 1:nz) * (tabs(j, 1:nz)-tabs(polidxs(i), 1:nz))) + dfrad(i))
          rad(j) = rad(j) * EXP(tmprad)
        ELSE
          rad(j) = rad(j) / (1.0 + frac1 * dfrad(i-1) + frac2 * dfrad(i))
        ENDIF

        IF ( which_polcorr /= 5 .OR. (which_polcorr == 5 .AND. (niter == 1 .OR. (niter == 0 .AND. currloop == 0))) ) THEN

          IF ( do_albwf)  THEN
            tmpcorr(1) = (1.0 + frac1 * dfalbwf(i-1) + frac2 * dfalbwf(i))
            IF (tmpcorr(1) /= 0.0) albwf(j) = albwf(j) / tmpcorr(1)
          ENDIF
          IF ( do_cfracwf)  THEN
            tmpcorr(1) = (1.0 + frac1 * dfcfracwf(i-1) + frac2 * dfcfracwf(i))
            IF (tmpcorr(1) /= 0.0) cfracwf(j) = cfracwf(j) / tmpcorr(1)
          ENDIF
          IF ( do_fozwf  ) THEN
            tmpcorr(1:nz) = (1.0 + frac1 * dffozwf(i-1, 1:nz) + frac2 * dffozwf(i, 1:nz))
            WHERE (tmpcorr(1:nz) /= 0.0)
              fozwf(j, :) = fozwf(j, :) / tmpcorr(1:nz)
            ENDWHERE
          ENDIF
          IF ( do_faerwf ) THEN
            tmpcorr(faerlvl:nz) = (1.0 + frac1 * dffaerwf(i-1, faerlvl:nz) + frac2 * dffaerwf(i, faerlvl:nz))
            WHERE (tmpcorr(faerlvl:nz) /= 0.0) 
              faerwf(j, faerlvl:nz)   = faerwf(j, faerlvl:nz)  / tmpcorr(faerlvl:nz)
            ENDWHERE
          ENDIF
          IF ( do_faerswf ) THEN
            tmpcorr(faerlvl:nz) = (1.0 + frac1 * dffaerswf(i-1, faerlvl:nz) + frac2 * dffaerswf(i, faerlvl:nz))
            WHERE (tmpcorr(faerlvl:nz) /= 0.0) 
              faerswf(j, faerlvl:nz)   = faerswf(j, faerlvl:nz)  / tmpcorr(faerlvl:nz)
            ENDWHERE
          ENDIF
          IF ( do_codwf )  THEN
            tmpcorr(nctp:ncbp) = (1.0 + frac1 * dffcodwf(i-1, nctp:ncbp) + frac2 * dffcodwf(i, nctp:ncbp))
            WHERE (tmpcorr(nctp:ncbp) /= 0.0)
              fcodwf(j, nctp:ncbp) = fcodwf(j, nctp:ncbp) / tmpcorr(nctp:ncbp)
            ENDWHERE
          ENDIF

          IF ( do_sprswf )  THEN
            tmpcorr(nsprs:nz) = (1.0 + frac1 * dffsprswf(i-1, nsprs:nz) + frac2 * dffsprswf(i, nsprs:nz))
            WHERE (tmpcorr(nsprs:nz) /= 0.0)
              fsprswf(j, nsprs:nz) = fsprswf(j, nsprs:nz) / tmpcorr(nsprs:nz)
            ENDWHERE
          ENDIF

          IF ( do_fraywf  ) THEN
            tmpcorr(1:nz) = (1.0 + frac1 * dffraywf(i-1, 1:nz) + frac2 * dffraywf(i, 1:nz))
            WHERE (tmpcorr(1:nz) /= 0.0)
              fraywf(j, :) = fraywf(j, :) / tmpcorr(1:nz)
            ENDWHERE
          ENDIF

        ENDIF
      ENDDO

      fidx = lidx + 1
    ENDDO

    IF (which_polcorr == 5) THEN
      IF (niter == 1) THEN
        fozwf_sav(1:nw, 1:nz) = fozwf(1:nw, 1:nz)
        albwf_sav(1:nw) = albwf(1:nw)
        cfracwf_sav(1:nw) = cfracwf(1:nw)
      ELSE IF (niter > 1 .OR. (niter == 0 .AND. currloop /= 0) ) THEN
        fozwf(1:nw, 1:nz) = fozwf_sav(1:nw, 1:nz)
        albwf(1:nw) = albwf_sav(1:nw)
        cfracwf(1:nw) = cfracwf_sav(1:nw)
      ENDIF
    ENDIF

    RETURN

  END SUBROUTINE polcorr_online


  SUBROUTINE get_efft(nz, zs, ozs, fts, ts, errstat)
    USE OMSAO_precision_module
    USE OMSAO_errstat_module

    INTEGER, INTENT(IN)                              :: nz
    INTEGER, INTENT(OUT)                             :: errstat
    REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
    REAL (KIND=dp), DIMENSION(0:nz),   INTENT(IN)    :: zs, fts
    REAL (KIND=dp), DIMENSION(1:nz),   INTENT(OUT)   :: ts

    INTEGER :: i, j, fidx, lidx, nz1

    REAL (KIND=dp), DIMENSION(0:nz) :: cumo3  
    REAL (KIND=dp), DIMENSION(0:nz * 10) :: zs1, ts1, cumo31
    REAL (KIND=dp), DIMENSION(nz * 10)   :: ozs1
    real (kind=dp) :: dz

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=8), PARAMETER :: modulename = 'get_efft'

    nz1 = nz * 10

    cumo3(0) = 0
    DO i = 1, nz
      cumo3(i) = cumo3(i-1) + ozs(i)
    ENDDO

    zs1(0) = zs(0)
    fidx = 1
    DO i = 1, nz 
      lidx = fidx + 9
      dz = (zs(i) - zs(i-1)) / 10.
      DO j = fidx, lidx
        zs1(j) = zs(i-1) + dz * (j - fidx + 1)
      ENDDO
      fidx  = lidx + 1
    ENDDO

    !WRITE(www_lun, '(12F8.3)') zs(0:nz)
    !WRITE(www_lun, *)
    !WRITE(www_lun, '(12F8.3)') zs1(0:nz1)

    CALL BSPLINE(zs, fts, nz+1, zs1(0:nz1), ts1(0:nz1), nz1+1, errstat)
    IF (errstat < 0) THEN
      WRITE(www_lun, *) modulename, ' : BSPLINE error, errstat = ', errstat
      RETURN
    ENDIF

    CALL BSPLINE(zs, cumo3, nz+1, zs1(0:nz1), cumo31(0:nz1), nz1+1, errstat)
    IF (errstat < 0) THEN
      WRITE(www_lun, *) modulename, ' : BSPLINE error, errstat = ', errstat
      RETURN
    ENDIF
    !print *, cumo3(nz), cumo31(nz1)

    ozs1 = cumo31(1:nz1) - cumo31(0:nz1-1)
    ts1(1:nz1) = (ts1(0:nz1-1) + ts1(1:nz1)) / 2.0

    fidx = 1
    DO i = 1, nz
      lidx = fidx + 9
      ts(i) = SUM(ts1(fidx:lidx) * ozs1(fidx:lidx)) / SUM(ozs1(fidx:lidx))
      fidx = lidx + 1    
    ENDDO

    RETURN

  END SUBROUTINE get_efft


  ! **************************************************************************
  ! When use_effcrs is set to true, outputs (radiance and weighting function) 
  ! correspond to input wavelengths. 
  ! When use_effcrs is set to false, i.e., use high resolution cross 
  ! section in the radiative transfer calculation and convolve with slot 
  ! functions after radiative transfer calculation, waves do not correspond 
  ! to the output rad and weighting function, waves(1:nw) store the 
  ! wavelengths (waves(1:ncalcp) where radiance calculation are done 
  ! (waves(ncalp+1:nw) = 0.0. And outputs corresponds to measurents (1:ns).
  ! Here nw is max(ncalcp, ns).
  ! **************************************************************************

  SUBROUTINE lidort_prof_env (do_ozwf, do_albwf, do_tmpwf, do_o3shi, ozvary, &
       do_taodwf, do_twaewf, do_saodwf, do_cfracwf, do_ctpwf, do_codwf, &
       do_sprswf,  do_so2zwf, nw, waves, nos, o3shi, sza, vza, aza, nl, &
       ozprof, tprof, n0alb, albarr, albpmin, albpmax, n0wfc, wfcarr, &
       wfcpmin, wfcpmax, nostk, albwf, ozwf, tmpwf, o3shiwf, cfracwf, &
       codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf, rad, errstat)

    USE OMSAO_precision_module
    USE OMSAO_parameters_module,ONLY  : max_fit_pts, du2mol, &
         rearth, max_spec_pts!, maxlay
    USE OMSAO_variables_module, ONLY  : radwvl_sav, n_radwvl_sav, &
         fitvar_rad, rmask_fitvar_rad, numwin, nradpix, nradpix_sav, &
         database, database_shiwf, database_save, refidx, band_selectors, &
         refspec_norm, &
         the_surfalt, currloop, refspec_orig_data, n_refspec_pts, winlim, &
         n_rad_wvl, tabdir!, vza_atm, refwvl, rad_wav_avg, n_refwvl, &
    !fitvar_rad_init, b1ab_div_wav, aza_atm, avgaza, avgsza, avgvza
    USE OMSAO_indices_module,   ONLY : so2_idx, so2v_idx, o2o2_idx
    USE ozprof_data_module,     ONLY  : nt_fit, atmosprof, nup2p, nflay, &
         do_lambcld, lambcld_refl, ntp, &
         do_subfit, polcorr, the_cod, the_cfrac, &
         nmom, fts, fps, fzs, fozs, frhos, gaext, gasca, gaasy, gamoms, gcq,  &
         gcasy, gcmoms, the_cbeta, aerwavs, has_clouds, aerosol, useasy, &
         aerwavs, actawin, ncbp, nctp, num_iter, do_tracewf, mgasprof, ngas,  &
         tracegas, fgasidxs, gasidxs, nallgas, raycof, depol, do_radinter,    &
         fgassidxs, strat_aerosol, osfind, oswins, mflay, nfsfc, &
         nsfc, do_simu, maxgksec, ngksec, radcalwrt, &
         tropsca, tropaod, tropwaer, strataod, stratsca, taodind,  &
         saodind, twaeind, sprsind, nwfc, use_effcrs, ncalcp, VlidortNstream!,&
    !radcidxs, saodfind, start_layer, twaefind, taodfind, wrtozcrs, &
    !gcw, gassidxs, do_multi_vza, do_ch2reso, end_layer, maxmom, maxawin
    USE OMSAO_errstat_module
    use m_ezspline_interpolation, only: bspline, bspline1
    use m_lidort_prof_prep, only: get_all_raycof_depol, &
         get_all_raycof_depol1, lidort_prof_prep
    use m_vlidort90_include
    !use lib_vlidort, only: vlidort_l_input_master
    !use vlidort_l_masters_lambertian, only: vlidort_l_master_lambertian
    !use vlidort_aux, only: vlidort_status

    IMPLICIT NONE

    !===========================  Define Variables ===========================
    ! Include files of dimensions and numbers
!    INCLUDE 'VLIDORT90.PARS'

    ! Include files of input variables
    INCLUDE 'VLIDORT_INPUTS90.VARS'
    INCLUDE 'VLIDORT_SETUPS90.VARS'
    INCLUDE 'VLIDORT_L_INPUTS90.VARS'
    INCLUDE 'VLIDORT_BOOKKEEP90.VARS'

    ! Include files of result variables
    INCLUDE 'VLIDORT_RESULTS90.VARS'
    INCLUDE 'VLIDORT_L_RESULTS90.VARS'

    !INTEGER, PARAMETER :: ngas=1 

    ! =======================
    ! Input/Output variables
    ! =======================
    INTEGER, INTENT(IN)  :: nw, nl, nos, n0alb, nostk, n0wfc
    LOGICAL, INTENT(IN)  :: do_ozwf, do_albwf, do_tmpwf, do_o3shi, do_taodwf, &
         do_twaewf, do_saodwf, do_cfracwf, do_codwf, do_ctpwf, do_sprswf, &
         do_so2zwf
    INTEGER, INTENT(OUT) :: errstat
    INTEGER, DIMENSION(n0alb), INTENT(IN)        :: albpmax, albpmin
    INTEGER, DIMENSION(n0wfc), INTENT(IN)        :: wfcpmax, wfcpmin
    LOGICAL, DIMENSION(nl), INTENT(IN)           :: ozvary
    REAL (KIND=dp), DIMENSION(nw),  INTENT(IN)   :: waves
    REAL (KIND=dp), DIMENSION(nw, nostk),  INTENT(OUT)   :: rad, albwf, &
         cfracwf, o3shiwf, codwf, ctpwf, taodwf, twaewf, saodwf, sprswf, so2zwf
    REAL (KIND=dp), DIMENSION(numwin, nos), INTENT(IN)   :: o3shi
    REAL (KIND=dp), DIMENSION(nl),  INTENT(IN)           :: ozprof, tprof
    REAL (KIND=dp), DIMENSION(nw, nl, nostk),INTENT(OUT) :: ozwf, tmpwf
    REAL (KIND=dp), DIMENSION(n0alb), INTENT(IN)         :: albarr 
    REAL (KIND=dp), DIMENSION(n0wfc), INTENT(IN)         :: wfcarr 
    REAL (KIND=dp), INTENT(IN)                           :: sza, vza, aza

    ! =======================
    ! Local variables
    ! =======================
    LOGICAL :: problems, do_clouds, do_fozwf, do_faerwf, do_fraywf!, ex
    INTEGER :: status_inputread, status_inputcheck, status_calculation, &
         ic, iw, i, j, k, kk, jj, jk, idum, low, hgh, fidx, lidx, &
         nz1, nfgas, nstep, istk, npolmod, ipol, nsprs, npts, nw0!, &
    !tempi, nvza, naza, ii
    INTEGER :: ozwfidx, aodwfidx, twaewfidx, codwfidx, sprswfidx, raywfidx
    LOGICAL, DIMENSION(nflay)                    :: cldmsk, varyprof!, aermsk
    REAL (KIND=dp)                               :: lamda, xg, frac, toz, &
         temp, aodscl, waerscl
    REAL (KIND=dp), DIMENSION(0:nflay)           :: ozs, delps
    REAL (KIND=dp), DIMENSION(nflay)             :: cldsca, cldext, cldasy, &
         cldext0, aersca, aerext, aerasy!, eta
    REAL (KIND=dp), DIMENSION(0:nmom,maxgksec,nflay) :: cldmoms, aermoms
    !REAL (KIND=dp), DIMENSION(0:nl)              :: cumoz
    !REAL (KIND=dp), DIMENSION(nl)                :: tmptwf, oztmpwf
    REAL (KIND=dp), DIMENSION(nw)                :: albs, tmpalbs, gshiwf, wfcs
    REAL (KIND=dp), DIMENSION(nw, 2, nostk)      :: polerr
    REAL (KIND=dp), DIMENSION(nw, nostk)         :: radclr, radcld
    REAL (KIND=dp), DIMENSION(n_radwvl_sav)      :: delpos, swaves, delshi
    !INTEGER,        DIMENSION(3)                 :: vind
    REAL (KIND=dp), DIMENSION(2, nostk)          :: radclrcld
    REAL (KIND=dp), DIMENSION(nw, nflay, nostk)  :: fozwf, faerwf, faerswf, &
         fcodwf, fsprswf, fraywf
    REAL (KIND=dp), DIMENSION(nw, nflay)         :: pfozwf, pfaerwf, &
         pfaerswf, pfcodwf, pfsprswf, pfraywf
    !REAL (KIND=dp), DIMENSION(nw, nl)            :: ptmpwf, pozwf
    REAL (KIND=dp), DIMENSION(nw)                :: prad, palbwf, pctpwf, &
         pcfracwf
    REAL (KIND=dp), DIMENSION (max_spec_pts)     :: tmprefwav, tmprefspec

    INTEGER, DIMENSION (nallgas)                  :: gasin   
    INTEGER, DIMENSION (5)          :: tmp_gasidxs = (/1, 4, 5, 7, 3/)
    REAL (KIND=dp), DIMENSION(nw, nallgas, nflay) :: allcrs
    REAL (KIND=dp), DIMENSION(nallgas, nflay)     :: alleta, allcol  
    REAL (KIND=dp), DIMENSION(nw, nflay) :: deltau, delsca, delray, &
         delo3abs, delabs

    LOGICAL, SAVE    :: first = .TRUE.
    INTEGER, SAVE    :: nz, faer_lvl, npolcorr, nradcal
    REAL (KIND=dp), DIMENSION(0:mflay),  SAVE :: ts
    !REAL (KIND=dp), DIMENSION(max_fit_pts, mflay),            SAVE :: aersca, aerext, aerasy 
    !REAL (KIND=dp), DIMENSION(max_fit_pts,0:maxmom,maxgksec,mflay), SAVE :: aermoms
    REAL (KIND=dp), DIMENSION(3, max_fit_pts),     SAVE :: abscrs_qtdepen
    LOGICAL, DIMENSION(mflay),                     SAVE :: aermsk
    LOGICAL, DIMENSION(max_fit_pts),       SAVE :: do_radcals, do_polcorrs
    INTEGER, DIMENSION(max_fit_pts),       SAVE :: polcorr_idxs, radcal_idxs
    REAL (KIND=dp), DIMENSION(max_fit_pts, mflay), SAVE :: abscrs, dads, dadt
    REAL (KIND=dp), DIMENSION(max_fit_pts, mflay), SAVE :: so2crs
    INTEGER, SAVE :: so2idx, so2vidx, o4idx, so2crsidx
    !real (kind=dp), dimension(nw,mflay) :: tmp_dads, tmp_dadt, tmp_abscrs, &
    !     tmp_so2crs

    ! Will become an input parameter later in ozprof.inp
    LOGICAL :: do_ssfullb295 = .FALSE.
    LOGICAL :: use_so2dtcrs  = .FALSE.  ! Use SO2 cross section at different temperatures

    ! ==============================
    ! Name of this module/subroutine
    ! ==============================
    CHARACTER (LEN=15), PARAMETER :: modulename = 'LIDORT_PROF_ENV'

    errstat = pge_errstat_ok

      status_inputcheck  = vlidort_success
      status_calculation = vlidort_success
      status_inputread   = vlidort_success
    IF (first) THEN
!      status_inputcheck  = vlidort_success
!      status_calculation = vlidort_success
!      status_inputread   = vlidort_success

      ! ======================= Read LIDORT Control Input ==========================
!!! FIXME 
!!! file path below should NOT be hard coded!
      CALL VLIDORT_L_INPUT_MASTER ( TRIM(tabdir)//'../control/vlidort_control.inp', &
           'o3prof_lidort_error', status_inputread)
      IF (status_inputread .NE. vlidort_success) THEN
        WRITE(www_lun, *) modulename, ': Problems encountered with input read!!!'
        errstat = pge_errstat_error
        RETURN
      ENDIF
      ngreek_moments_input = nmom
      earth_radius = rearth

      IF (.NOT. aerosol .AND. (.NOT. has_clouds .OR. do_lambcld)) THEN
        DO_SSCORR_TRUNCATION = .FALSE.
        DO_DELTAM_SCALING = .FALSE.
        DO_RAYLEIGH_ONLY = .TRUE.
        DO_SOLUTION_SAVING = .FALSE.
        DO_BVP_TELESCOPING = .FALSE.
      ENDIF

      DO i = 1, ngas
        IF (gasidxs(i) == so2v_idx)   so2vidx = i  
        IF (gasidxs(i) == so2_idx)    so2idx  = i    
        IF (gasidxs(i) == o2o2_idx)    o4idx  = i
      ENDDO

      IF (fgasidxs(so2vidx) <= 0 .AND. fgasidxs(so2idx) <= 0) use_so2dtcrs = .FALSE.

      ! Could not handle two different wavelength shifts for the same cross sections
      IF (fgassidxs(so2vidx) > 0 .AND. fgassidxs(so2idx) > 0 .AND. use_so2dtcrs) use_so2dtcrs = .FALSE.

      first = .FALSE.
    ENDIF

    ! ============= Overridden some control and atmospheric variables ============== 
    IF (num_iter == 0 ) THEN
      n_szangles = 1
      szangles(1) = sza 
      IF (sza >= 90.0 .OR. sza < 0) THEN
        WRITE(www_lun, *) modulename, ' : SZA is >= 90 or < 0 !!!'
        errstat = pge_errstat_error
        RETURN
      ENDIF

      n_user_vzangles = 1
      n_user_relazms = 1
      user_vzangles(1) = vza
      user_relazms(1) = aza

      !IF (.NOT. do_multi_vza) THEN
      !ELSE
      !   n_user_streams = 3; n_user_relazms = 3
      !   user_angles_input(1:3) = ABS(vza_atm); user_relazms(1:3) = aza_atm
      !   ! vza will be sorted automatically in LIDORT, need to get order index
      !   vind(1) = MINVAL(MINLOC(ABS(vza_atm)))
      !   vind(3) = MINVAL(MAXLOC(ABS(vza_atm)))
      !   vind(2) = 6 - vind(1) - vind(3)
      !   IF (vind(2) < 1 .OR. vind(2) > 3) THEN
      !      WRITE(www_lun, *) modulename, ' : Sth. wrong in view geometry!!!'
      !      errstat = pge_errstat_error; RETURN
      !   ENDIF
      !
      !   vind((/vind(1), vind(2), vind(3)/)) = (/1, 2, 3/)
      !   naza =3; nvza = 3
      !ENDIF
      nz = nflay
      nlayers = nz

      IF (nz > maxlayers) THEN
        WRITE(www_lun, *) modulename, ' : # of layers exceeded allowed !!!'
        errstat = pge_errstat_error
        RETURN
      ENDIF
      IF (nl > nlayers) THEN
        WRITE(www_lun, *) modulename, ' : Coarse grids cannot be finer than fine grids!!!'
        errstat = pge_errstat_error
        RETURN
      ENDIF

      ! Set height grid for doing Chapman Function Calculation 
      height_grid(0:nlayers) = fzs(0:nz)
      geometry_specheight  = the_surfalt
      !print *, the_surfalt, height_grid(nlayers), fzs(nz)

      ! set the first aerosol layer (from TOA down) 1:faer_lvl-1 (without aerosols)
      faer_lvl = 1
      IF (.NOT. strat_aerosol) faer_lvl = nup2p(ntp) + 1

      IF (aerosol) THEN 
        aermsk = .TRUE.
        IF (.NOT. strat_aerosol) aermsk(1:faer_lvl-1) = .FALSE.
      ELSE
        aermsk = .FALSE.
      ENDIF

      ts(1:nz) = (fts(1:nz) + fts(0:nz-1)) / 2.0
      !WRITE(www_lun, '(12F8.2)') ts(1:nz)
      !CALL get_efft(nz, fzs(0:nz), fozs(1:nz), fts(0:nz), ts(1:nz), errstat)
      !WRITE(www_lun, *)
      !WRITE(www_lun, '(12F8.2)') ts(1:nz)
    ENDIF

    ! =================== Interpolate Ozone, T to fine grids =======================
    IF (nz /= nl) THEN
      DO i = 1, nl
        varyprof(nup2p(i-1)+1:nup2p(i)) = ozvary(i)
        ozs(nup2p(i-1)+1:nup2p(i)) = fozs(nup2p(i-1)+1:nup2p(i)) * ozprof(i) / &
             SUM(fozs(nup2p(i-1)+1:nup2p(i)))
      ENDDO

      IF (nt_fit > 0) THEN

        DO i = 1, nl
          atmosprof(3, i) = tprof(i) * 2.0 - atmosprof(3,  i-1)
        ENDDO

        CALL BSPLINE(atmosprof(2, 0:nl), atmosprof(3, 0:nl), &
             nl+1, fzs(0:nz), ts(0:nz), nz+1, errstat)
        IF (errstat < 0) THEN
          WRITE(www_lun, *) modulename, ' : BSPLINE error, errstat = ', errstat
          RETURN
        ENDIF
        ts(1:nz) = (ts(1:nz) + ts(0:nz-1)) / 2.0
      ENDIF
    ELSE
      ozs(1:nz) = ozprof(1:nl)
      varyprof(1:nz) = ozvary(1:nl)
    ENDIF

    nz1 = nfsfc - 1 
    ! Update aerosol fields: first AOD
    IF (do_taodwf .AND. num_iter > 0) THEN
      aodscl = fitvar_rad(taodind) / tropaod(actawin)
      tropaod(1:actawin) = tropaod(1:actawin) * aodscl
      gaext(1:actawin, nup2p(ntp)+1:nz1) = gaext(1:actawin, nup2p(ntp)+1:nz1) * aodscl

      IF (.NOT. do_twaewf) THEN ! Single scattering albedo does not change
        tropsca(1:actawin) = tropsca(1:actawin) * aodscl
        gasca(1:actawin, nup2p(ntp)+1:nz1) = gasca(1:actawin, nup2p(ntp)+1:nz1) * aodscl
      ENDIF
    ELSE
      aodscl = 1.0
    ENDIF

    IF (do_twaewf .AND. num_iter > 0) THEN
      waerscl = fitvar_rad(twaeind) / tropwaer(actawin)   ! Scale single scattering albedo
      tropwaer(1:actawin) = tropwaer(1:actawin) * waerscl
      tropsca(1:actawin)  = tropsca(1:actawin) * waerscl * aodscl
      gasca(1:actawin, nup2p(ntp)+1:nz1) = gasca(1:actawin, nup2p(ntp)+1:nz1) * waerscl * aodscl
    ENDIF

    IF (do_saodwf .AND. num_iter > 0) THEN
      aodscl = fitvar_rad(saodind) / strataod(actawin)
      strataod(1:actawin) = strataod(1:actawin) * aodscl
      stratsca(1:actawin) = stratsca(1:actawin) * aodscl
      gaext(1:actawin, 1:nup2p(ntp)) = gaext(1:actawin, 1:nup2p(ntp)) * aodscl
      gasca(1:actawin, 1:nup2p(ntp)) = gasca(1:actawin, 1:nup2p(ntp)) * aodscl
    ENDIF

    IF (do_sprswf .AND. num_iter > 0) THEN
      temp = (fitvar_rad(sprsind) - fps(nup2p(nsfc-1)))/ (fps(nz1) - fps(nup2p(nsfc-1)))
      frhos(nup2p(nsfc-1)+1:nz1) = frhos(nup2p(nsfc-1)+1:nz1) * temp
      delps(nup2p(nsfc-1)+1:nz1) = (fps(nup2p(nsfc-1)+1:nz1)-fps(nup2p(nsfc-1):nz1-1))  * temp
      DO i = nup2p(nsfc-1)+1, nz1
        fps(i) = fps(i-1) + delps(i)
      ENDDO
    ENDIF

    ! albedo array
    DO i = 1, n0alb
      albs(albpmin(i):albpmax(i)) = albarr(i)
    ENDDO
    DO i = 1, n0wfc
      wfcs(wfcpmin(i):wfcpmax(i)) = wfcarr(i)
    ENDDO

    ! Determinine atmospheric weighting functions to be calculated
    IF (do_ozwf .OR.  do_tmpwf .OR. do_o3shi .OR. do_taodwf .OR. &
         do_twaewf .OR. do_saodwf .OR. do_codwf .OR. do_sprswf )  do_atmos_linearization = .TRUE.
    do_surface_linearization = do_albwf    
    IF (do_atmos_linearization .OR. do_albwf) THEN
      do_simulation_only =    .FALSE.
      do_linearization = .TRUE.
    ELSE 
      do_simulation_only =    .TRUE.
      do_linearization = .FALSE.
    ENDIF

    i = 0
    ozwfidx = 0
    aodwfidx = 0
    twaewfidx = 0
    codwfidx = 0
    sprswfidx = 0
    layer_vary_flag(1:nz) = .FALSE.
    layer_vary_number(1:nz) = 0
    IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi) THEN
      i = i + 1
      ozwfidx = i
      profilewf_names(i) = 'ozone volume mixing ratio------'

      layer_vary_flag(1:nz) = varyprof(1:nz)
      WHERE (varyprof(1:nz) )
        layer_vary_number(1:nz) = layer_vary_number(1:nz) + 1
      ENDWHERE
    ENDIF
    IF ( do_taodwf .OR. do_saodwf) THEN
      i = i + 1
      aodwfidx = i
      profilewf_names(i) = 'aerosol extinction coefficient-'

      IF (do_taodwf) THEN
        layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(ntp)+1:nz1) = layer_vary_number(nup2p(ntp)+1:nz1) + 1
      ENDIF

      IF (do_saodwf) THEN
        layer_vary_flag(1:nup2p(ntp)) = .TRUE.
        layer_vary_number(1:nup2p(ntp)) = layer_vary_number(1:nup2p(ntp)) + 1
      ENDIF
    ENDIF
    IF ( do_twaewf) THEN
      i = i + 1
      twaewfidx = i
      profilewf_names(i) = 'aerosol scattering coefficient-'
      layer_vary_flag(nup2p(ntp)+1:nz1) = .TRUE.
      layer_vary_number(nup2p(ntp)+1:nz1) = layer_vary_number(nup2p(ntp)+1:nz1) + 1
    ENDIF
    IF ( do_codwf) THEN
      i = i + 1
      codwfidx = i
      profilewf_names(i) = 'cloud extinction coefficient---'
      layer_vary_flag(nctp:ncbp) = .TRUE.
      layer_vary_number(nctp:ncbp) = layer_vary_number(nctp:ncbp) + 1
    ENDIF
    do_fraywf = .FALSE.
    IF ( do_sprswf  .OR. (.NOT. use_effcrs .AND. nw > 1) ) THEN
      ! Need to use jacobians wrt rayleigh OD to perform interpolation

      i = i + 1
      sprswfidx = i
      raywfidx = i
      profilewf_names(i) = 'rayleigh optical thickness-----'
      IF (.NOT. use_effcrs ) THEN
        do_fraywf = .TRUE.
        layer_vary_flag(1:nz1) = .TRUE.
        layer_vary_number(1:nz1) = layer_vary_number(1:nz1) + 1
      ELSE
        layer_vary_flag(nup2p(nsfc-1)+1:nz1) = .TRUE.
        layer_vary_number(nup2p(nsfc-1)+1:nz1) = layer_vary_number(nup2p(nsfc-1)+1:nz1) + 1
      ENDIF
    ENDIF

    n_totalatmos_wfs = i

    IF (n_totalatmos_wfs > 1) THEN
      WHERE(layer_vary_number(1:nz) > 0) 
        layer_vary_number(1:nz) = n_totalatmos_wfs
      ENDWHERE
    ENDIF

    do_fozwf = .FALSE.
    IF (do_ozwf .OR. do_tmpwf .OR. do_o3shi .OR. (.NOT. use_effcrs .AND. nw > 1)) do_fozwf = .TRUE.
    do_faerwf = .FALSE.
    IF (do_taodwf .OR. do_saodwf) do_faerwf = .TRUE.

    !print *, n_totalatmos_wfs, nz
    !print *, layer_vary_flag(1:nz)
    !print *, layer_vary_number(1:nz)

    ! ==================== Get Ozone Absorption Cross Section ====================   
    IF (nw > 1) THEN
      allcol(1, 1:nz1) = ozs(1:nz1) * du2mol
      nfgas = 1
      gasin(1) = 1

      DO k = 1, ngas
        IF (fgasidxs(k) > 0) THEN
          nfgas = nfgas + 1
          gasin(nfgas) = nfgas
          IF (use_so2dtcrs) so2crsidx = nfgas

          ! molecules cm^-2, but normalized by refspec_norm(gasidxs(k))
          ! allcol / refspec_norm will be molecules cm^-2
          allcol(nfgas, 1:nz1) = mgasprof(k, 1:nz1) * tracegas(k, 4) / mgasprof(k, nz+1)
        ENDIF
      ENDDO

      IF (use_effcrs) THEN
        IF (nos > 0) THEN
          delshi = 0.0

          IF (do_subfit) THEN
            fidx = 1
            DO j = 1, numwin
              lidx = fidx + nradpix_sav(j) - 1
              delpos(fidx:lidx) =  radwvl_sav(fidx:lidx) - (radwvl_sav(fidx) + radwvl_sav(lidx)) / 2.0
              IF (osfind(j, 1) > 0) delshi(fidx:lidx) =  o3shi(j, 1) 

              DO i = 2, nos
                IF (osfind(j, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx)  + &
                     o3shi(j, i) * delpos(fidx:lidx) ** (i-1)
              ENDDO
              fidx = lidx + 1
            ENDDO
          ELSE
            IF (oswins(1, 1) == 1) THEN
              fidx = 1
            ELSE
              fidx = SUM(nradpix(1: oswins(1, 1)-1)) + 1
            ENDIF
            lidx = SUM(nradpix(1: oswins(1, 2)))

            delpos(fidx:lidx) =  radwvl_sav(fidx:lidx) - (radwvl_sav(fidx) + radwvl_sav(lidx)) / 2.0
            IF (osfind(1, 1) > 0) delshi(fidx:lidx) =  + o3shi(1, 1) 

            DO i = 2, nos  
              IF (osfind(1, i) > 0) delshi(fidx:lidx) = delshi(fidx:lidx)  + &
                   o3shi(1, i) * delpos(fidx:lidx) ** (i-1)
            ENDDO
          ENDIF

          swaves = radwvl_sav - delshi        
        ELSE
          swaves = radwvl_sav
        ENDIF

        ! Get ozone absorption coefficients
        IF (ANY(ABS(delshi) >= 2.0)) THEN
          WRITE(www_lun, *) modulename, ' : Ozone wavelength shifts are too large!!!'
          errstat = pge_errstat_error
          RETURN
        ENDIF

        ! Get temperature-dependent ozone cross section at instrumental spectral resolution
        IF (num_iter == 0 .OR. (do_o3shi .AND. nw > 1) ) THEN
          CALL GETABS_CRS(swaves, n_radwvl_sav, nw, 1, nz1, ts(1:nz1), abscrs(1:nw, 1:nz1), &
               do_o3shi, do_tmpwf, dads(1:nw, 1:nz1), dadt(1:nw, 1:nz1), problems, &
               abscrs_qtdepen(1:3, 1:nw))

          IF (problems) THEN
            WRITE(www_lun, *) modulename, ' : Problems in reading ozone absorption !!!'
            errstat = pge_errstat_error
            RETURN
          ENDIF
        ENDIF
        allcrs(1:nw, 1, 1:nz1) = abscrs(1:nw, 1:nz1)

        ! Get Rayleigh scattering coefficients and deplolarizaiton factor
        !IF (num_iter == 0) THEN 
        !   CALL GET_ALL_RAYCOF_DEPOL(nw, waves, raycof(1:nw), depol(1:nw))
        IF (num_iter == 0) THEN 
          CALL GET_ALL_RAYCOF_DEPOL1(n_radwvl_sav, radwvl_sav, nw, raycof(1:nw), depol(1:nw), problems)
        ENDIF

        IF ( num_iter == 0 .AND. (fgasidxs(so2idx) > 0 .OR. fgasidxs(so2vidx) > 0) .AND. use_so2dtcrs) THEN
          CALL GETSO2_CRS(radwvl_sav, n_radwvl_sav, nw, nz1, ts(1:nz1), so2crs(1:nw, 1:nz1), problems)
          IF (problems) THEN
            WRITE(www_lun, *) modulename, ' : Problems in reading SO2 absorption !!!'
            errstat = pge_errstat_error
            RETURN
          ENDIF
        ENDIF

        nfgas = 1
        DO k = 1, ngas
          IF (fgasidxs(k) > 0) THEN
            nfgas = nfgas + 1
            gasin(nfgas) = nfgas

            IF ((k /= so2idx .AND. k /= so2vidx) .OR. .NOT. use_so2dtcrs) THEN
              IF (fgassidxs(k) > 0 ) THEN
                npts = n_refspec_pts(gasidxs(k))
                tmprefwav(1:npts) = refspec_orig_data(gasidxs(k), 1:npts, 1)
                tmprefspec(1:npts) = refspec_orig_data(gasidxs(k), 1:npts, 3)
                fidx = MINVAL(MINLOC(waves(1:nw), MASK=(waves(1:nw) >= &
                     tmprefwav(1) + 0.1 .AND. waves(1:nw) <= tmprefwav(npts) - 0.1)))
                lidx = MINVAL(MAXLOC(waves(1:nw), MASK=(waves(1:nw) >= &
                     tmprefwav(1) + 0.1 .AND. waves(1:nw) <= tmprefwav(npts) - 0.1)))

                IF (lidx > fidx .AND. lidx > 0 .AND. fidx > 0) THEN 
                  temp = fitvar_rad(rmask_fitvar_rad(fgassidxs(k)))
                  CALL BSPLINE1(  tmprefwav(1:npts) - temp, tmprefspec(1:npts), npts, &
                       waves(fidx:lidx), allcrs(fidx:lidx, nfgas, 1),  gshiwf(fidx:lidx), lidx-fidx+1, errstat) 

                  !print *, k, gasidxs(k), fgasidxs(k), gassidxs(k), fgassidxs(k), fitvar_rad(gassidxs(k)), npts
                  !WRITE(90, '(F10.4, 2D14.6)') ((waves(i), database_shiwf(gasidxs(k), refidx(i)), gshiwf(i)), i=fidx, lidx)

                  database_shiwf(gasidxs(k), refidx(fidx:lidx)) = gshiwf(fidx:lidx)
                  database(gasidxs(k), refidx(fidx:lidx)) = allcrs(fidx:lidx, nfgas, 1)               
                  IF (errstat < 0) THEN
                    WRITE(www_lun, *) modulename, ' : BSPLINE error, errstat = ', errstat
                    RETURN
                  ENDIF
                ENDIF
              ELSE
                allcrs(1:nw, nfgas, 1) = database_save(gasidxs(k), refidx(1:nw))
              ENDIF

              DO i = 2, nz1
                allcrs(1:nw, nfgas, i) = allcrs(1:nw, nfgas, 1)
              ENDDO
            ELSE
              allcrs(1:nw, nfgas, 1:nz1) = so2crs(1:nw, 1:nz1) / refspec_norm(gasidxs(k))
            ENDIF
          ENDIF
        ENDDO

        IF (do_radinter) THEN
          IF (num_iter == 0) THEN
            do_radcals(1:nw) = .FALSE.
            fidx = 1
            k = nz1 / 2

            DO iw = 1, numwin
              lidx = fidx + nradpix(iw) - 1

              ! Do radiative transfer calculations at local maxima/minima
              ! always do radiative transfer calculation at end points
              IF (band_selectors(iw) == 2) THEN
                DO i = fidx + 1, lidx - 1
                  IF (abscrs(i, k) > abscrs(i-1, k) .AND. abscrs(i, k) > abscrs(i+1, k)) &
                       do_radcals(i) = .TRUE.
                  IF (abscrs(i, k) < abscrs(i-1, k) .AND. abscrs(i, k) < abscrs(i+1, k)) &
                       do_radcals(i) = .TRUE.
                  !IF (abscrs(i, nz1) > abscrs(i-1, nz1) .AND. abscrs(i, nz1) > abscrs(i+1, nz1)) &
                  !    do_radcals(i) = .TRUE.
                  !IF (abscrs(i, nz1) < abscrs(i-1, nz1) .AND. abscrs(i, nz1) < abscrs(i+1, nz1)) &
                  !    do_radcals(i) = .TRUE.
                ENDDO

                nstep = 5
                do_radcals(fidx) = .TRUE.
                do_radcals(lidx-1:lidx) = .TRUE.                
              ELSE IF (band_selectors(iw) == 1) THEN
                nstep = 3
                do_radcals(fidx) = .TRUE.
                do_radcals(lidx) = .TRUE.
              ELSE
                nstep = 1
              ENDIF

              DO i = fidx + 1, lidx - 1, nstep
                do_radcals(i) = .TRUE.
              ENDDO

              fidx = lidx + 1        
            ENDDO

            nradcal = 0
            DO i = 1, nw
              IF (do_radcals(i)) THEN
                nradcal = nradcal + 1
                radcal_idxs(nradcal) = i
              ENDIF
            ENDDO
          ENDIF
        ELSE
          IF (num_iter == 0) THEN
            do_radcals(1:nw) = .TRUE.
            nradcal = nw
            DO i = 1, nw
              radcal_idxs(i) = i
            ENDDO
          ENDIF
        ENDIF

        !WRITE(77, *) nw, nz1
        !WRITE(77, '(30F8.3)') ts(1:nz1)
        !WRITE(77, '(30F8.3)') ozs(1:nz1)
        !DO i = 1, nw
        !   WRITE(77, '(60D14.6)') waves(i), abscrs(i, 12) !, dads(i, 1:nz1)
        !ENDDO  
        !CLOSE (77)

        !WRITE(91, *) nz1
        !WRITE(91, '(30D14.6)')  frhos(1:nz1)
        !WRITE(91, '(30D14.6)')  fzs(0:nz1)
        !WRITE(91, *) nw
        !DO i = 1, nw
        !   WRITE(91, '(F10.4, D14.6)') waves(i), raycof(i)
        !ENDDO
      ELSE   ! .NOT. use_effcrs
        ! O3/SO2 (use_so2dtcrs=.TRUE.) cross section: if do_tmpwf = .FALSE. and do_o3shi is false, 
        ! just need to get once for each retrieval
        ! Other trace gas cross section: just need to get it once for all the retrievals if no shifts 
        CALL GET_HRES_GASCRS_RAY(nw, waves, nz1, ts(1:nz1), do_o3shi, o3shi, do_tmpwf, nfgas, &
             use_so2dtcrs, num_iter, allcrs(1:nw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), frhos(1:nz1), &
             abscrs_qtdepen(1:3, 1:nw), raycof(1:nw), depol(1:nw), errstat)
        abscrs(1:ncalcp, 1:nz1) = allcrs(1:ncalcp, 1, 1:nz1)
        IF (use_so2dtcrs) so2crs(1:ncalcp, 1:nz1) = allcrs(1:ncalcp, so2crsidx, 1:nz1)

        IF (num_iter == 0) THEN
          do_radcals(1:nw) = .FALSE.
          nradcal = ncalcp
          do_radcals(1:ncalcp) = .TRUE.
          radcal_idxs(1:ncalcp) = (/(i, i=1, ncalcp)/)
        ENDIF
      ENDIF

      !IF (wrtozcrs) THEN
      !   WRITE(91, *) nw
      !   DO iw = 1, nw 
      !      WRITE(91, '(F10.4, 6D14.6)') waves(iw), abscrs_qtdepen(1:3, iw), raycof(iw), depol(iw) !,  &
      !           !database(6, refidx(iw))*refspec_norm(6), database(9, refidx(iw))*refspec_norm(9), &
      !           !database(10, refidx(iw))*refspec_norm(10), database(12, refidx(iw))*refspec_norm(12)
      !   ENDDO
      !   
      !   WRITE(91, *) nz1
      !   WRITE(91, '(2D14.6)') fps(0), fzs(0)
      !   DO i = 1, nz1 
      !      WRITE(91, '(10D14.6)') fps(i), fzs(i), ts(i), frhos(i), allcol(1, i), allcol(2, i)/refspec_norm(6), &
      !           allcol(3, i)/refspec_norm(9), allcol(4, i)/refspec_norm(10), allcol(5, i)/refspec_norm(12) 
      !   ENDDO
      !ENDIF

    ELSE
      do_radcals(1) = .TRUE.
      nradcal = 1
      radcal_idxs(1) = 1

      ! o3 absorption coefficient at 370.2 nm with TOMS FWHM
      dads(1, 1:nz1) = 0.0
      dadt(1, 1:nz1) = 0.0

      ! Weighted by solar flux
      !abscrs(1, 1:nz) = 9.1231787D-24 + (ts(1:nz) - 273.15) * 1.9005502D-25 + &
      !     (ts(1:nz) - 273.15)**2.0 * 1.2275286D-27     
      !raycof(1) = 2.3184501D-26; depol(1) = 0.030247913D0

      nfgas = 7
      CALL GET_ALB_OZCRS_RAY(nz1, ts(1:nz1), nfgas, allcrs(1, 1:nfgas, 1:nz1), raycof(1), depol(1), problems)
      IF (problems) THEN
        WRITE(www_lun, *) modulename, ' : Problems in reading O3 XSec for determining Fc!!!'
        errstat = pge_errstat_error
        RETURN
      ENDIF

      !CALL GET_ALL_RAYCOF_DEPOL(1, 360.D0, raycof(1), depol(1))
      !print *, raycof(1), depol(1)

      ! Ozone
      nfgas = 6
      gasin(1) = 1
      allcol(1, 1:nz1) = ozs(1:nz1) * du2mol

      DO i = 1, 4
        allcol(i+1, 1:nz1) = mgasprof(tmp_gasidxs(i), 1:nz1)
        gasin(i+1) = i + 1
      ENDDO

      ! O2-O2 Concentration (Oxygen: 20.95%)
      allcol(6, 1:nz1) = ( frhos(1:nz1) * 0.2095 ) ** 2.0 / (fzs(0:nz1-1) - fzs(1:nz1)) / 1.0D5
      gasin(6) = 6
    ENDIF

    ! Initialize aerosol/cloud property profiles
    IF (num_iter == 0) THEN
      aersca(1:nz1) = 0.0
      aerext(1:nz1) = 0.0
      aerasy(1:nz1) = 0.0
      aermoms(0:nmom, 1:ngksec, 1:nz1) = 0.0
    ENDIF

    cldmsk = .FALSE.
    IF (nctp /= 0) cldmsk(nctp:ncbp)=.TRUE.
    cldsca=0.0
    cldext=0.0
    cldasy=0.0
    cldmoms = 0.0

    ! Initialize output variables  
    rad = 0.0
    radclr = 0.0
    radcld = 0.0
    cfracwf = 0.0
    IF ( do_fozwf )  fozwf   = 0.0
    IF ( do_albwf)   albwf   = 0.0
    IF ( do_codwf )  fcodwf  = 0.0
    IF ( do_sprswf ) fsprswf  = 0.0
    IF ( do_fraywf ) fraywf  = 0.0
    IF ( do_ctpwf )  ctpwf   = 0.0
    IF ( do_faerwf)  faerwf  = 0.0
    IF ( do_twaewf ) faerswf = 0.0

    ! senstitivity of absorption cross section to temperature, used 
    ! for calculating temperature wf directly with LIDORT
    !eta = 0.0             ! dummy variable here
    alleta = 0.0

    polerr = 1.0
    IF ( polcorr == 2) THEN 
      IF (nwfc > 0) THEN
        print *, 'Polarization correction + variable fc: not implemented!!!'
        STOP 1
      ENDIF
      tmpalbs = albs 
      ! clear-sky
      IF (the_cfrac /= 1.0) THEN
        toz = SUM(ozs(1:nz1))

        WHERE(tmpalbs >= 1.0)
          tmpalbs = 0.999
        ENDWHERE

        WHERE(tmpalbs < 0.0)
          tmpalbs = 0.001
        ENDWHERE

        !        CALL lup_polerror(toz, fps(nz1), nw, waves, tmpalbs, polerr(1:nw, 1, 1), errstat)
        IF (errstat == pge_errstat_error) RETURN
        polerr(1:nw, 1, 1) = 1.0 + polerr(1:nw, 1, 1) / 100.0
      ENDIF

      ! cloudy part
      IF (the_cfrac /= 0.0) THEN
        toz = SUM(ozs(1:nctp-1))
        IF (do_lambcld) THEN
          tmpalbs = lambcld_refl
        ELSE
          ! approxmation of cod as cloud albedo for polarization correction
          ! Follows Kokhanovsky et al., 2003
          tmpalbs = 1.0 - 1.0 / (1.072 + 0.1125 * the_cod)
        ENDIF

        WHERE(tmpalbs >= 1.0)
          tmpalbs = 0.999
        ENDWHERE

        WHERE(tmpalbs < 0.0)
          tmpalbs = 0.001
        ENDWHERE

        !        CALL lup_polerror(toz, fps(nctp-1), nw, waves, tmpalbs, polerr(1:nw, 2, 1), errstat)
        IF (errstat == pge_errstat_error) RETURN
        polerr(1:nw, 2, 1) = 1.0 + polerr(1:nw, 2, 1) / 100.0
      ENDIF
    ENDIF

    ! Determine wavelengths where exact polarization correction (NSTOKES: 4 vs 1) is calculated
    ! In UV1 (or between 270 and 310 nm): ~292 nm, ~298 nm, ~300 nm, ~302 nm, ~304 nm, ~306 nm, last wavelength 
    ! In UV2 (or between 310 and 340 nm): first, 1/4, middle and last wavelength
    ! So exact vector LIDORT calculation is done at 11 wavelengths.
    ! This option works when radiance interpolation option is turned on
    IF (num_iter == 0) THEN
      do_polcorrs(1:nw) = .FALSE.
      npolcorr = 0
      IF ( (polcorr >= 3 .AND. polcorr <= 5) .AND. nw > 1 ) THEN
        fidx = 1
        idum = 0
        DO iw = 1, numwin
          IF (iw == numwin) THEN
            temp = winlim(iw, 2)
          ELSE
            temp = (winlim(iw, 2) + winlim(iw + 1, 1)) / 2.
          ENDIF
          lidx = MINVAL(MAXLOC(waves(1:nw), MASK=(waves(1:nw) < temp .AND. waves(1:nw) > 0))) 
          !lidx = fidx + nradpix(iw) - 1
          IF ( waves(lidx) <= 312.0 ) THEN
            ! Error from using single scattering is about 0.2% at 270 nm, it needs o be corrected
            !IF (do_ssfullb295) do_polcorrs(fidx) = .TRUE. 
            IF (idum == 0) idum = 1
            DO i = fidx + 1, lidx - 1
              IF (do_radcals(i) ) THEN
                IF ( (waves(idum) < 290. .AND. waves(i) >= 290.0) .OR.  &
                     (waves(idum) < 295. .AND. waves(i) >= 295.0) .OR.  &
                     (waves(idum) < 299. .AND. waves(i) >= 299.0) .OR.  &
                     (waves(idum) < 301. .AND. waves(i) >= 301.0) .OR.  &
                     (waves(idum) < 303. .AND. waves(i) >= 303.) .OR.  &
                                !(waves(idum) < 304. .AND. waves(i) >= 304.0) .OR.  &
                     (waves(idum) < 305. .AND. waves(i) >= 305.0) )  do_polcorrs(i) = .TRUE.
                idum = i
              ENDIF
            ENDDO
            do_polcorrs(lidx) = .TRUE.
          ELSE
            IF (idum == 0) idum = 1
            do_polcorrs(fidx) = .TRUE.
            do_polcorrs(lidx) = .TRUE.     
            j = (fidx + lidx ) / 2
            k = fidx + (j - fidx) / 3  

            jj = (j + lidx )  / 2
            kk = fidx + (k - fidx) / 3
            jk = (j + k ) / 2

            DO i = fidx + 1, lidx - 1
              IF (do_radcals(i) ) THEN
                !IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                !     (waves(idum) < waves(k) .AND. waves(i) >= waves(k)) .OR. &
                !     (waves(idum) < waves(jj) .AND. waves(i) >= waves(jj)) .OR. &
                !     (waves(idum) < waves(kk) .AND. waves(i) >= waves(kk)) .OR. &
                !     (waves(idum) < waves(jk) .AND. waves(i) >= waves(jk)))  do_polcorrs(i) = .TRUE.
                IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                     (waves(idum) < waves(k) .AND. waves(i) >= waves(k)))  do_polcorrs(i) = .TRUE.
                idum = i
              ENDIF
            ENDDO
          ENDIF

          fidx = lidx + 1
        ENDDO
        DO i = 1, nw
          IF ( do_polcorrs(i) ) THEN
            !IF (do_ssfullb295 .AND. npolcorr == 2 .AND. waves(1) < 290.0) THEN
            !   ! Use single scattering (ss), scalar below 295 nm
            !   ! Use ss + multiple scattering, scalar above 295 nm
            !   ! Correction: 270 (ss + scalr), 290  (ss + scalar),
            !   !             first lamda before 295 (ss + scalar)
            !   !             first lamda after  295 (scalar)
            !   do_polcorrs(i-1) = .TRUE.; npolcorr = npolcorr + 2
            !   polcorr_idxs(3) = i-1; polcorr_idxs(4) = i
            !ELSE
            npolcorr = npolcorr + 1
            polcorr_idxs(npolcorr) = i
            !ENDIF
          ENDIF
        ENDDO

        !print *, nradcal, npolcorr
        !print *, polcorr_idxs(1:npolcorr)
        !print *, do_radcals(polcorr_idxs(1:npolcorr))
        !print *, waves(polcorr_idxs(1:npolcorr))
        !STOP
      ENDIF
    ENDIF
    ! ====================== Call LIDORT and Do Post Processing =====================
    DO iw = 1, nw 

      IF (.NOT. do_radcals(iw) ) CYCLE   ! Radiances and weighting functions will be interpolated 
      lamda = waves(iw) 
      IF (nwfc > 0) the_cfrac = wfcs(iw)

      ! NSTOKES = 4 when
      ! (1) Vector LIDORT or 
      ! (2) on-line polarization correction (polcorr=3) at cloud wavelength or wavelengths 
      !     where exact radiances are calcualted or
      ! (3) on-line polarization correction (polcorr=4) at cloud wavelength or wavelengths 
      !     where exact radiances are calcualted in the first iteration
      IF (polcorr == 0 .OR. ((polcorr == 3 .OR. polcorr == 5) .AND. (nw == 1 .OR. do_polcorrs(iw)))  &
           .OR. (polcorr == 4 .AND. (nw == 1 .OR. (do_polcorrs(iw) .AND. num_iter == 0))) ) THEN
        NSTOKES = 3 
        NSTREAMS = 4
      ELSE
        NSTOKES = 1 
        NSTREAMS = 4
      ENDIF
      VlidortNstream = NSTREAMS

      IF ( lamda < 295.0 .AND. NSTOKES == 1 .AND. nw > 1 .AND. do_ssfullb295) THEN
        DO_SSFULL = .TRUE.
        DO_SSCORR_TRUNCATION = .FALSE.
        DO_DELTAM_SCALING = .FALSE.
      ELSE
        DO_SSFULL = .FALSE.
        IF (aerosol .OR. (has_clouds .AND. .NOT. do_lambcld)) THEN
          DO_SSCORR_TRUNCATION = .TRUE.
          DO_DELTAM_SCALING = .TRUE.
        ENDIF
      ENDIF

      IF ( aerosol ) THEN  !.AND. num_iter == 0 ) THEN
        ! Interpolate/extrapoalte for aerosol properties
        hgh = actawin
        DO i = 1, actawin
          IF (lamda < aerwavs(i)) THEN
            hgh = i
            EXIT
          ENDIF
        ENDDO
        IF (hgh==1) hgh = 2    ! extrapolation
        low = hgh - 1
        xg = (lamda - aerwavs(low)) / (aerwavs(hgh) - aerwavs(low))

        aersca(faer_lvl:nfsfc-1) = gasca(low, faer_lvl:nfsfc-1) * (1.0 - xg) +  gasca(hgh, faer_lvl:nfsfc-1) * xg
        aerext(faer_lvl:nfsfc-1) = gaext(low, faer_lvl:nfsfc-1) * (1.0 - xg) +  gaext(hgh, faer_lvl:nfsfc-1) * xg

        !IF (num_iter == 0) THEN  ! Don't need to be updated
        IF (useasy) THEN
          aerasy(faer_lvl:nfsfc-1) = gaasy(low, faer_lvl:nfsfc-1) * (1.0 - xg) + &
               gaasy(hgh, faer_lvl:nfsfc-1) * xg
        ELSE        
          DO i = 0, nmom
            DO j = 1, ngksec
              aermoms(i, j, faer_lvl:nfsfc-1) = &
                   gamoms(low, faer_lvl:nfsfc-1, i, j) * (1.0 - xg) +  &
                   gamoms(hgh, faer_lvl:nfsfc-1, i, j) * xg
            ENDDO
          ENDDO
        ENDIF
        !ENDIF
      ENDIF

      IF (has_clouds .AND. .NOT. do_lambcld) THEN
        cldext(nctp:ncbp) = (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) * &
             (fzs(nctp-1 : ncbp-1) - fzs(nctp:ncbp)) * the_cbeta
        cldsca(nctp:ncbp) = cldext(nctp:ncbp)  ! w0 = 1.0
        IF (iw == 1) cldext0(nctp:ncbp) = cldext(nctp:ncbp) * gcq(actawin) / (gcq(low) * (1.0 - xg) + gcq(hgh) * xg) 

        IF (useasy) THEN
          cldasy(nctp:ncbp) = gcasy(low) * (1.0 - xg) + gcasy(hgh) * xg
        ELSE
          DO i = 0, nmom
            DO j = 1, ngksec
              cldmoms(i, j, nctp:ncbp) = gcmoms(low, i, j) &
                   * (1.0 - xg) + gcmoms(hgh, i, j) * xg 
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF ( do_polcorrs(iw) .AND. ((polcorr == 3 .OR. polcorr == 5) &
           .OR. (polcorr == 4 .AND. num_iter == 0) ) ) THEN
        npolmod = 2   ! Twice, one vector and one scalar          
      ELSE 
        npolmod = 1   ! Only once either scalar or vector
      ENDIF

      ! When polcorr = 5 is selected, only calculate weighting function for iteration 1
      ! iteration 0 if 1st pixel of a x-track position is being retrieved
      IF (polcorr == 5 .AND. nw > 1 .AND. (num_iter > 1 .OR. (num_iter == 0 .AND. currloop /= 0))) THEN
        do_simulation_only = .TRUE.
        do_linearization   = .FALSE.
      ENDIF

      DO ipol = 1, npolmod
        IF (ipol == 2) THEN
          NSTOKES = 1 
          NSTREAMS = 4 ! Always scalar for second mode

          IF ( lamda < 295.0 .AND. nw > 1 .AND. do_ssfullb295) THEN
            DO_SSFULL = .TRUE.
            DO_SSCORR_TRUNCATION = .FALSE.
            DO_DELTAM_SCALING = .FALSE.
          ELSE
            DO_SSFULL = .FALSE.
            DO_SSCORR_TRUNCATION = .TRUE.
            DO_DELTAM_SCALING = .TRUE.
            !DO_SSCORR_TRUNCATION = .FALSE.; DO_DELTAM_SCALING = .FALSE.
          ENDIF
          IF (.NOT. aerosol) THEN
            DO_SSCORR_TRUNCATION = .FALSE.
            DO_DELTAM_SCALING = .FALSE.
          ENDIF

          ! Save VECTOR LIDORT results
          prad(iw) = rad(iw, 1)
          rad(iw, 1:nostk) = 0.d0
          IF (do_cfracwf) pcfracwf(iw) = cfracwf(iw, 1)
          IF (do_cfracwf) cfracwf(iw, 1:nostk) =  0.d0

          IF ( do_linearization ) THEN
            IF (do_albwf)  palbwf(iw)             = albwf(iw, 1)
            IF (do_fozwf)  pfozwf(iw, :)          = fozwf(iw, :, 1)
            IF (do_codwf)  pfcodwf(iw, nctp:ncbp) = fcodwf(iw, nctp:ncbp, 1)
            IF (do_sprswf) pfsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1) = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1)
            IF (do_fraywf) pfraywf(iw, :)         = fraywf(iw, :, 1)
            IF (do_ctpwf)  pctpwf(iw)             = ctpwf(iw, 1)
            IF (do_faerwf) pfaerwf(iw, faer_lvl:nfsfc-1)  = faerwf(iw, faer_lvl:nfsfc-1, 1)           
            IF (do_twaewf) pfaerswf(iw, faer_lvl:nfsfc-1) = faerwf(iw, faer_lvl:nfsfc-1, 1)

            ! Initalize those variables to zero again          
            IF (do_albwf) albwf(iw, 1:nostk) = 0.d0
            IF (do_fozwf) fozwf(iw, :, 1:nostk) = 0.d0 
            IF (do_codwf) fcodwf(iw, nctp:ncbp, 1:nostk) = 0.d0
            IF (do_sprswf) fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, 1:nostk) = 0.d0
            IF (do_fraywf) fraywf(iw, :, 1:nostk) = 0.d0
            IF (do_ctpwf)  ctpwf(iw, 1:nostk)       = 0.d0
            IF (do_faerwf) faerwf(iw, faer_lvl:nfsfc-1, 1:nostk)  = 0.d0
            IF (do_twaewf) faerswf(iw, faer_lvl:nfsfc-1, 1:nostk) = 0.d0

          ENDIF
        ENDIF

        !print *, iw, ipol, do_simulation_only, do_linearization, do_atmos_linearization, do_surface_linearization

        radclrcld = 0.0
        DO ic = 1, 2               ! for clear and cloud

          IF (ic == 1) THEN
            do_clouds = .FALSE.
            frac = 1.0 - the_cfrac
          ELSE
            do_clouds = .TRUE. 
            frac = the_cfrac
          ENDIF
          !WRITE(www_lun, *) ic, frac, do_clouds, do_lambcld, lambcld_refl, do_multi_vza
          IF (frac == 0.0) CYCLE  ! No clear/cloudy part

          ! Note convert ozone from DU to molecule/cm^2  here
          IF ((ic == 1) .OR. (.NOT. do_lambcld))  THEN
            lambertian_albedo = albs(iw)

            ! Reset up number of layers since Lambertian 
            ! cloudy scene got different layers
            nlayers = nfsfc - 1
            nz1 = nlayers   

            ! zero O3 weighting function below surface (This is necessary)           
            profilewf(1:n_totalatmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0

            IF (ipol == 1) THEN
              !WRITE(91, '(4F10.4,3I5)') sza, vza, aza, lambertian_albedo, 1, nlayers, ngreek_moments_input
              !WRITE(91, '(40D14.6)') height_grid(0:nlayers)
              !WRITE(91, '(2D14.6)') lamda, depol(iw)

              CALL LIDORT_PROF_PREP(raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
                   varyprof(1:nz1), nfgas, gasin(1:nfgas), allcrs(iw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), &
                   alleta(1:nfgas, 1:nz1), useasy, nmom, aerosol, aersca(1:nz1),      &
                   aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
                   do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
                   cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), problems, &
                   deltau(iw, 1:nz1), delsca(iw, 1:nz1), delo3abs(iw, 1:nz1), delray(iw, 1:nz1))
              IF (problems) then
                WRITE(*, *) modulename, ' : Problems encountered in lidort preparation!!!'
                errstat = pge_errstat_error
                RETURN
              END IF
            ENDIF
          ELSE                 ! lambertian clouds
            nlayers = nctp-1  ! from cloud top to TOA
            nz1 = nlayers
            do_clouds = .FALSE.

            IF (the_cfrac == 1.0 .AND. nw /= 1) THEN
              lambertian_albedo = albs(iw)
              lambcld_refl = lambertian_albedo
            ELSE
              lambertian_albedo = lambcld_refl ! use 80% (could be adjusted when the_cfrac gt 0.90)
            ENDIF

            ! zero O3 weighting function (This is necessary)           
            profilewf(1:n_totalatmos_wfs, nz1+1:nz, 1, 1, 1:NSTOKES, 1) = 0.0

            IF (frac == 1.0 .AND. ipol == 1) THEN
              CALL LIDORT_PROF_PREP(raycof(iw), depol(iw), fzs(0:nz1), frhos(1:nz1), &
                   varyprof(1:nz1), nfgas, gasin(1:nfgas), allcrs(iw, 1:nfgas, 1:nz1), allcol(1:nfgas, 1:nz1), &
                   alleta(1:nfgas, 1:nz1), useasy, nmom, aerosol, aersca(1:nz1),      &
                   aerext(1:nz1), aerasy(1:nz1), aermoms(0:nmom, 1:maxgksec, 1:nz1), aermsk(1:nz1), &
                   do_clouds, cldsca(1:nz1), cldext(1:nz1), cldasy(1:nz1), &
                   cldmoms(0:nmom, 1:maxgksec, 1:nz1), cldmsk(1:nz1), problems , &
                   deltau(iw, 1:nz1), delsca(iw, 1:nz1), delo3abs(iw, 1:nz1), delray(iw, 1:nz1))

              IF (problems) THEN
                WRITE(www_lun, *) modulename, ' : Problems encountered in lidort preparation!!!'
                errstat = pge_errstat_error
                RETURN
              END IF
            ENDIF
          ENDIF

          !print *, iw, lamda, ic, frac, nlayers 
          ! No Need to reset viewing geometry for each LIDORT call
          !geometry_specheight = height_grid(nlayers)
          !szangles(1) = sza; user_vzangles(1) = vza; user_relazms(1) = aza  
          !print *, nlayers 
          !print *, geometry_specheight, height_grid(0:nlayers)
          CALL VLIDORT_L_MASTER_LAMBERTIAN (status_inputcheck, status_calculation)   
          !IF (iw >= 102 .AND. iw <= 103) THEN
          !   WRITE(www_lun, '(3I5,F8.3,30D14.6)') iw, ic, ipol, lamda, &
          !    stokes(1, 1, 1:nostk, 1)!, atmoswf(ozwfidx, 1:nfsfc-1, 1, 1, 1, 1)
          !ENDIF

          !CALL VLIDORT_STATUS ( status_inputcheck, status_calculation ) 
          IF (status_inputcheck   /= vlidort_success) THEN
            WRITE(www_lun, *) modulename, ' : Problems encountered in lidort input check!!!'
            errstat = pge_errstat_error
            RETURN
          ENDIF
          IF (status_calculation  /= vlidort_success) THEN
            WRITE(www_lun, *) modulename, ' : Problems encountered in lidort calculation!!!'
            errstat = pge_errstat_error
            RETURN
          ENDIF

          ! Pixel-independent approximation
          radclrcld(ic, 1:nostk) = stokes(1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk)
          rad(iw, 1:nostk)       = rad(iw, 1:nostk) + radclrcld(ic, 1:nostk) * frac

          IF (do_linearization ) THEN

            ! weighting function per Dobson Unit
            IF ( do_ozwf ) THEN
              DO istk = 1, nostk
                fozwf(iw, 1:nfsfc-1, istk) = fozwf(iw, 1:nfsfc-1, istk) + profilewf(ozwfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                     / ozs(1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
              ENDDO
            ENDIF

            ! Weighting function with respect to aerosol/cloud optical depth at the last wavelength
            ! so as to keep the scaling since we are fitting the aod at that wavelength
            IF ( do_taodwf .OR. do_saodwf ) THEN
              !print *, ipol, faer_lvl, nfsfc-1, frac
              DO istk = 1, nostk
                faerwf(iw, faer_lvl:nfsfc-1, istk)  = faerwf(iw, faer_lvl:nfsfc-1, istk) + &
                     profilewf(aodwfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                     gaext(actawin, faer_lvl:nfsfc-1) * polerr(iw, ic, istk) *  frac
                !print *, faerwf(iw, faer_lvl:nfsfc-1, 1)
              ENDDO
            ENDIF
            IF ( do_twaewf ) THEN
              DO istk = 1, nostk
                faerswf(iw, faer_lvl:nfsfc-1, istk)  = faerswf(iw, faer_lvl:nfsfc-1, istk) + &
                     profilewf(twaewfidx, faer_lvl:nfsfc-1, 1, 1, istk, 1) / &
                     gasca(actawin, faer_lvl:nfsfc-1) * gaext(actawin, faer_lvl:nfsfc-1) * &
                     polerr(iw, ic, istk) *  frac
              ENDDO
            ENDIF
            IF ( do_codwf ) THEN
              DO istk = 1, nostk
                fcodwf(iw, nctp:ncbp, istk)  = fcodwf(iw, nctp:ncbp, istk) + &
                     profilewf(codwfidx, nctp:ncbp, 1, 1, istk, 1) / &
                     cldext0(nctp:ncbp) * polerr(iw, ic, istk) *  frac
              ENDDO
            ENDIF
            IF ( do_sprswf ) THEN
              DO istk = 1, nostk 
                fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk)  = fsprswf(iw, nup2p(nsfc-1)+1:nfsfc-1, istk) + &
                     profilewf(sprswfidx, nup2p(nsfc-1)+1:nfsfc-1, 1, 1, istk, 1) / &
                     (fps(nup2p(nsfc-1)+1:nfsfc-1) - fps(nup2p(nsfc-1):nfsfc-2)) * polerr(iw, ic, istk) *  frac
              ENDDO
            ENDIF

            IF ( do_fraywf ) THEN
              DO istk = 1, nostk
                fraywf(iw, 1:nfsfc-1, istk) = fraywf(iw, 1:nfsfc-1, istk) + profilewf(raywfidx, 1:nfsfc-1, 1, 1, istk, 1) &
                     / delray(iw, 1:nfsfc-1)  * polerr(iw, ic, istk) * frac 
              ENDDO
            ENDIF

            ! Non-lambertian clouds: albedo wf is from both clear and cloud
            ! Lamberitan clouds:     if the_cfrac  < 1, albwf from clear only
            !                        if the_cfrac == 1, albwf from cloud only
            IF (do_albwf) THEN 
              IF (.NOT. do_lambcld) THEN
                albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + surfacewf(1, 1, 1, 1:nostk, 1) &
                     * polerr(iw, ic, 1:nostk) * frac / lambertian_albedo
              ELSE IF (the_cfrac < 1.0) THEN
                IF (ic == 1) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
                     surfacewf(1, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac / lambertian_albedo
              ELSE IF (the_cfrac == 1.0) THEN
                IF (ic == 2) albwf(iw, 1:nostk) = albwf(iw, 1:nostk) + &
                     surfacewf(1, 1, 1, 1:nostk, 1) * polerr(iw, ic, 1:nostk) * frac / lambertian_albedo
              ENDIF
            ENDIF

          ENDIF

        ENDDO ! end clear/cloudy scene loop
        radclr(iw, 1:nostk) = radclrcld(1, 1:nostk)
        radcld(iw, 1:nostk) = radclrcld(2, 1:nostk)
        IF (do_cfracwf) cfracwf(iw, 1:nostk) = radcld(iw, 1:nostk) - radclr(iw, 1:nostk)
      ENDDO    ! end scalar/vector modes
    ENDDO       ! end wavelength loop 

    ! On-line polarization correction
    IF ( the_cfrac == 1.0 .AND. do_lambcld) THEN
      nz1 = nctp - 1
    ELSE
      nz1 = nfsfc - 1
    ENDIF

    IF ( (polcorr >= 3 .AND. polcorr <= 5 ) .AND. nw > 1) THEN
      IF (use_effcrs) THEN
        nw0 = nw
      ELSE
        nw0 = ncalcp
      ENDIF

      nsprs = nup2p(nsfc-1)+1
      delabs(1:nw0, 1:nz1) = deltau(1:nw0, 1:nz1) - delsca(1:nw0, 1:nz1)
      CALL polcorr_online(num_iter, polcorr, nw0, nz1, nctp, ncbp, nsprs, &
           faer_lvl, npolcorr, polcorr_idxs(1:npolcorr), do_fozwf, do_albwf, &
           do_faerwf, do_twaewf, do_codwf, do_sprswf, do_fraywf, do_cfracwf, &
           waves(1:nw0), rad(1:nw0, 1), prad(1:nw0), delabs(1:nw0, 1:nz1), &
           ozs(1:nz1), albwf(1:nw0, 1), palbwf(1:nw0), &
!           fozwf(1:nw0, 1:nz1, 1), pfozwf(1:nw0, 1:nz1), &
!           faerwf(1:nw0, 1:nz1, 1), pfaerwf(1:nw0, 1:nz1), &
!           faerswf(1:nw0, 1:nz1, 1), pfaerswf(1:nw0, 1:nz1), &
!           fcodwf(1:nw0, 1:nz1, 1), pfcodwf(1:nw0, 1:nz1), &
!           fsprswf(1:nw0, 1:nz1, 1), pfsprswf(1:nw0, 1:nz1), &
!           fraywf(1:nw0, 1:nz1, 1), pfraywf(1:nw0, 1:nz1), &
           fozwf(1:nw0, 1:nz1, 1), pfozwf(1:nw0, 1:nz1), &
           faerwf(1:nw0, 1:nz1, 1), pfaerwf(1:nw0, 1:nz1), &
           faerswf(1:nw0, 1:nz1, 1), pfaerswf(1:nw0, 1:nz1), &
           fcodwf(1:nw0, 1:nz1, 1), pfcodwf(1:nw0, 1:nz1), &
           fsprswf(1:nw0, 1:nz1, 1), pfsprswf(1:nw0, 1:nz1), &
           fraywf(1:nw0, 1:nz1, 1), pfraywf(1:nw0, 1:nz1), &
           cfracwf(1:nw0, 1), pcfracwf(1:nw0))
    ENDIF

    ! Radiance Interpolation
    IF (nw > 1 .AND. do_radinter ) THEN
      nsprs = nup2p(nsfc-1)+1
      CALL radwf_interpol(nw, nz1, nctp, ncbp, nsprs, faer_lvl, &
           do_radcals(1:nw), do_fozwf, do_albwf, do_faerwf, do_twaewf, &
           do_codwf, do_sprswf, do_cfracwf, waves, abscrs(1:nw, 1:nz1), &
           ozs(1:nz1), rad(1:nw, 1), fozwf(1:nw, 1:nz1, 1), albwf(1:nw, 1), &
           cfracwf(1:nw, 1), faerwf(1:nw, 1:nz1, 1), &
           faerswf(1:nw, 1:nz1, 1), fcodwf(1:nw, 1:nz1, 1), &
           fsprswf(1:nw, 1:nz1, 1), errstat)
      IF (errstat == pge_errstat_error) RETURN
    ENDIF

    IF (nw > 1 .AND. .NOT. use_effcrs) THEN
      IF (do_fraywf .AND. do_fozwf) THEN 
        nsprs = nup2p(nsfc-1)+1
        CALL hres_radwf_inter_convol(nw, nz1, nctp, ncbp, nsprs, faer_lvl, &
             do_albwf, do_faerwf, do_twaewf, do_codwf, do_sprswf, do_cfracwf, &
             do_tracewf, use_so2dtcrs, do_o3shi, do_tmpwf, waves, &
             ozs(1:nz1), rad(1:nw, 1), fozwf(1:nw, 1:nz1, 1), &
             albwf(1:nw, 1), cfracwf(1:nw, 1), faerwf(1:nw, 1:nz1, 1), &
             faerswf(1:nw, 1:nz1, 1), fcodwf(1:nw, 1:nz1, 1), &
             fsprswf(1:nw, 1:nz1, 1), fraywf(1:nw, 1:nz1, 1), &
             dads(1:nw, 1:nz1), dadt(1:nw, 1:nz1), abscrs(1:nw, 1:nz), &
             so2crs(1:nw, 1:nz), errstat)
        IF (errstat == pge_errstat_error) RETURN
      ELSE
        WRITE(www_lun, *) 'Must have O3/rayleigh weighting function to do radiance interpolation!!!'
        errstat = pge_errstat_error
        RETURN
      ENDIF
    ENDIF


    !WRITE(71, *) nw, nflay
    !DO iw = 1, nw
    !   IF (do_radcals(iw) ) THEN
    !      write(71, *) '1 '
    !   ELSE
    !      write(71, *) '0 '
    !   ENDIF
    !ENDDO
    !    
    !DO iw = 1, nw 
    !   WRITE(71, '(1000D16.7)') waves(iw), rad(iw, 1), fozwf(iw, 1:nflay, 1), albwf(iw, 1), &
    !        deltau(iw, 1:nflay), delsca(iw, 1:nflay), delo3abs(iw, 1:nflay)
    !ENDDO
    !WRITE(71, '(10D16.7)') ozs(1:nflay)
    !print *, fozwf(100, 1:nflay, 1)
    !
    !STOP

    IF (nw > 1 .AND. .NOT. use_effcrs) THEN
      nw0 = n_rad_wvl
    ELSE
      nw0 = nw
    ENDIF

    ! Calculate desired weighting functions at the end after applying all the correction
    IF ( do_ozwf ) THEN  
      DO i = 1, nl
        fidx = nup2p(i - 1) + 1
        lidx = nup2p(i)         
        DO iw = 1, nw0
          DO istk = 1, nostk
            ozwf(iw, i, istk) = SUM(fozwf(iw, fidx:lidx, istk) * ozs(fidx:lidx)) / SUM(ozs(fidx:lidx))
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF ( do_tmpwf ) THEN  
      DO i = 1, nl
        fidx = nup2p(i - 1) + 1
        lidx = nup2p(i)         
        DO iw = 1, nw0
          DO istk = 1, nostk
            tmpwf(iw, i, istk) = SUM(fozwf(iw, fidx:lidx, istk) * ozs(fidx:lidx) &
                 * dadt(iw, fidx:lidx) * (fzs(fidx-1:lidx-1)-fzs(fidx:lidx))) / (fzs(fidx-1) - fzs(lidx))
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (do_o3shi) THEN
      DO iw = 1, nw0
        DO istk = 1, nostk
          o3shiwf(iw, istk) = SUM(fozwf(iw, 1:nz1, istk) * ozs(1:nz1) * dads(iw, 1:nz1))
        ENDDO
      ENDDO
    ENDIF

    IF (do_taodwf) THEN
      DO iw = 1, nw0
        DO istk = 1, nostk
          taodwf(iw, istk) = SUM( faerwf(iw, nup2p(ntp)+1:nz1, istk) * &
               gaext(actawin, nup2p(ntp)+1:nz1) ) / SUM (gaext(actawin, nup2p(ntp)+1:nz1))
        ENDDO
      ENDDO
    ENDIF

    IF (do_twaewf) THEN
      DO iw = 1, nw0
        DO istk = 1, nostk
          twaewf(iw, istk) = SUM( faerswf(iw, nup2p(ntp)+1:nz1, istk) * &
               gasca(actawin, nup2p(ntp)+1:nz1) / gaext(actawin, nup2p(ntp)+1:nz1) ) &
               / SUM (gasca(actawin, nup2p(ntp)+1:nz1) / gaext(actawin, nup2p(ntp)+1:nz1))
        ENDDO
      ENDDO
    ENDIF

    IF (do_saodwf) THEN
      DO iw = 1, nw0
        DO istk = 1, nostk
          saodwf(iw, istk) = SUM( faerwf(iw, 1:nup2p(ntp), istk) * &
               gaext(actawin, 1:nup2p(ntp)) ) / SUM (gaext(actawin, 1:nup2p(ntp)))
        ENDDO
      ENDDO
    ENDIF

    IF (do_codwf) THEN
      DO iw = 1, nw0
        DO istk = 1, nostk
          codwf(iw, istk) = SUM( fcodwf(iw, nctp:ncbp, istk) * &
               cldext0(nctp:ncbp) ) / SUM (cldext0(nctp:ncbp))
        ENDDO
      ENDDO
    ENDIF

    IF (do_sprswf) THEN
      DO iw = 1, nw0
        DO istk = 1, nostk
          sprswf(iw, istk) = SUM( fsprswf(iw, nup2p(nsfc-1)+1:nz1, istk) * &
               (fps(nup2p(nsfc-1)+1:nz1) - fps(nup2p(nsfc-1):nz1-1))) / (fps(nz1)-fps(nup2p(nsfc-1)))
        ENDDO
      ENDDO
    ENDIF

    !IF (nw > 1) THEN
    !   DO iw = 1, nw
    !      WRITE(78, '(4D16.7)') waves(iw), albwf(iw, 1), sprswf(iw, 1), cfracwf(iw, 1)
    !   ENDDO
    !   STOP
    !ENDIF

    IF (nw > 1 .AND. do_ozwf .AND. do_tracewf ) &
         CALL GET_TRACEGAS_WF (fozwf(1:nw0, 1:nz, 1), abscrs(1:nw0, 1:nz), so2crs(1:nw0, 1:nz), use_so2dtcrs, &
         rad(1:nw0, 1), nw0, nz, nz1, ozs(1:nz), waves(1:nw0), do_so2zwf, so2zwf(1:nw0, 1))

    IF (nw > 1 .AND. do_simu .AND. .NOT. radcalwrt) THEN
      WRITE(78, *) 2, nz1
      WRITE(78, *) 'Profile: Z, P, T, O3, SO2'
      DO i = 1, nz1
        WRITE(78, '(F10.4, 4D16.7)') fzs(i), fps(i), fts(i), fozs(i), allcol(3, i)/refspec_norm(9)
      ENDDO
      WRITE(78, '(A)') 'Wavelength, radiance, albedo wf, ozone wf, aerosol wf, ozcrs, so2crs'
      DO i = 1, nw !160, 369, 209
        WRITE(78, '(F10.4, 500D16.7)') waves(i), rad(i, 1), albwf(i, 1), fozwf(i, 1:nz1, 1), & !, faerwf(i, 1:nz1, 1), &
             allcrs(i, 1, 1:nz1), allcrs(i, 3, nz1)*refspec_norm(9)
      ENDDO
      print *, nw, nz
      STOP 1
    ENDIF

    RETURN
  END SUBROUTINE lidort_prof_env



end module lidort_prof_utilities
