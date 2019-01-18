!
MODULE m_lidort_util

  USE m_ezspline_interpolation, only: bspline, bspline2, interpol, interpol2
  USE m_avg_band, ONLY: avg_band_spec
  USE m_convol, ONLY: convol_f2c
  PUBLIC get_hres_radcal_waves, get_hres_gascrs_ray, &
       hres_radwf_inter_convol, &
       get_slant_tau, get_tracegas_wf, radwf_interpol, get_efft

  INTEGER, PARAMETER, PRIVATE :: max_pathlen = 1024

CONTAINS


  !1.	Establish fine wavelength grid: 0.01 nm now, may change to 0.05 nm later
  !2.	Establish radiance calculation grid, based on spectral sampling rate for 
  !   different specral regions
  !3.	Find indices of radiance calculation grid in fine wavelength grid

SUBROUTINE get_hres_radcal_waves(errstat)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY  : max_spec_pts, max_fit_pts
  USE OMSAO_variables_module, ONLY  : numwin, nradpix, solwinfit, which_slit, &
       curr_rad_spec, use_redfixwav, winlim
  USE OMSAO_indices_module,   ONLY  : hwe_idx, wvl_idx
  USE ozprof_data_module,     ONLY  : radc_msegsr, radc_nsegsr, radc_samprate, radc_lambnd,  &
       hreswav, radcwav, nhresp, ncalcp, radcidxs, nhresp0, hreswav0, hres_samprate
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! Output variables
  INTEGER, INTENT(OUT) :: errstat

  ! Local variables
  REAL (KIND=dp), PARAMETER :: dhw0 = 0.01  ! at 0.01 nm
  INTEGER, PARAMETER        :: mextraw = 100

  INTEGER              :: i, j, k, fidx, lidx, nsub, nratio, nhalf, nextra, n0
  REAL (KIND=dp)       :: tmp, swav, ewav, slw, samprate, invdhw, ds1, ds2
  REAL (KIND=dp), DIMENSION (mextraw) :: extrawave

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
  nhresp0 = nhresp; hreswav0(1:nhresp) = hreswav(1:nhresp)
  invdhw = 1.0 / hres_samprate; nratio = hres_samprate / dhw0
  nhalf =  nratio / 2
  
  j = 1
  DO i = nhalf + 1, nhresp0 - nhalf, nratio
     fidx = i - nhalf; lidx = i - nhalf + nratio - 1
     hreswav(j) = SUM(hreswav(fidx:lidx)) / REAL(nratio)
     j = j + 1
  ENDDO
  nhresp = j - 1
  
  IF (nhresp0 > max_spec_pts .OR. nhresp > max_spec_pts) THEN
     WRITE(*, *) modulename, ': Need to increase max_spec_pts!!!'
     errstat = pge_errstat_error; RETURN
  ENDIF
  
  ! Establish radiance calculation grid, based on spectral sampling
  ! Sampling rate are specified in several segments (e.g., < 295 nm, 295-308 nm, > 308 nm)
  ! Make sure that radiance will be done for the first and last points. 
  radc_samprate = FLOOR(radc_samprate * invdhw) / invdhw  ! multiples of dhw
  ncalcp = 1; radcwav(1) = hreswav(1)
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
  WRITE (www_lun,*) 'N of hreswav:', nhresp
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
  !   WRITE(*, *) modulename, ': Need to increase mextraw!!!'
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
     WRITE(*, *) modulename, ': Need to increase max_fit_pts!!!'
     errstat = pge_errstat_error; RETURN
  ENDIF    

  ! Find the indices of radiance calc. wavelength in high resolution
  radcidxs(1:ncalcp) = 0
  radcidxs(1) = 1; j = 2
  DO i = 2, nhresp
     IF ( ABS(hreswav(i) - radcwav(j)) <= hres_samprate * 0.2 ) THEN
        radcidxs(j) = i; j = j + 1
     ENDIF
  ENDDO
  
  RETURN
END SUBROUTINE get_hres_radcal_waves

SUBROUTINE hres_radwf_inter_convol(nw, nz, nctp, ncbp, nsprs, faerlvl,  &
     do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf, do_tracewf, &
     do_o3shi, do_tmpwf, do_pslwf, wave, ozs, rad, fozwf, albwf, cfracwf, faerwf, &
     faerswf, fcodwf, fsprswf, fraywf, dads, dadt, abscrs, &
     so2crs, o4crs, o2crs, h2ocrs, errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,   ONLY  : hwe_idx, spk_idx, &
      so2_idx, so2v_idx, o2o2_idx,o2_idx, o2t2_idx, h2o_idx, h2ot2_idx
  USE OMSAO_parameters_module,ONLY  : du2mol
  USE OMSAO_variables_module, ONLY  : numwin, nradpix, band_selectors, winlim, &
       owave=>radwvl_sav, now=>n_radwvl_sav,onpix=>nradpix_sav, i0sav, refidx, fitwavs, & 
       nrad=>n_rad_wvl, &
       do_bandavg, curr_rad_spec, refidx_sav, database, database_shiwf, &
       database_pslwf, solwinfit, solwinfit_save, npsl,psl_fpos,max_psl,do_dsdw, do_dsdk
  USE ozprof_data_module,     ONLY  : nup2p, hwave=>hreswav, radcwav, &
       radcidxs, hres_i0, nhw=>nhresp, hresgabs, hresray, nw0=>ncalcp, &
       hres_gas, hres_gasshi, ngas, gasidxs, fgasidxs, fgassidxs, &
       o3crsz, o3dadtz, o3dadsz, so2crsz, o4crsz, o2crsz, h2ocrsz, &
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                              :: nw, nz, nctp, ncbp, faerlvl, nsprs
  INTEGER, INTENT(OUT)                             :: errstat                                                   
  LOGICAL, INTENT(IN)                              :: do_albwf, do_faerwf, do_faerswf, &
       do_codwf, do_sprswf, do_cfracwf, do_o3shi, do_tmpwf, do_tracewf, do_pslwf

  REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fozwf, faerwf, faerswf, fcodwf, &
       fsprswf, fraywf
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(OUT)   :: dads, dadt, abscrs, so2crs, o4crs, o2crs, h2ocrs
  REAL (KIND=dp), DIMENSION(nw),     INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: rad, albwf, cfracwf

  ! Local variables
  INTEGER, PARAMETER :: which_pslwf = 2
  INTEGER :: i, j, iwin, fidx, lidx, fidxc, lidxc, idx, iw, ntemp, nspec, sidx, eidx
  LOGICAL :: do_so2shi, do_o4shi, do_o2shi, do_h2oshi
  REAL (KIND=dp)                      :: temp
  INTEGER, DIMENSION (nw)             :: c2hfidx, c2hlidx

  REAL (KIND=dp), DIMENSION (nhw)     :: hrad,hrad1,halbwf, tmparr, dtau, dray, hcfracwf
  REAL (KIND=dp), DIMENSION (now)     :: oi0, otmp, tmpi0 !, so2dads1, o4dads1
  REAL (KIND=dp), DIMENSION (nw,  nz) :: tauwf
  REAL (KIND=dp), DIMENSION (now, nz) :: dads1, dadt1, abscrs1, so2crs1, o4crs1,o2crs1, h2ocrs1
  !REAL (KIND=dp), DIMENSION (ngas,now):: tmp_gas , tmp_gasshi
  REAL (KIND=dp), DIMENSION (nhw, nz) :: hozwf, haerwf, haerswf, hcodwf, hsprswf, hraywf 
  REAL (KIND=dp), DIMENSION (nhw, nz*8) :: inarr
  REAL (KIND=dp), DIMENSION (now, nz*8) :: outarr, outarr1
  REAL (KIND=dp), DIMENSION (now) :: dpabs
  !INTEGER :: ntime = 1

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'hres_radwf_convol'

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
  ! it is alsoe minus of altitude-depedent air mass factor
  !DO i = 1, nz
  !   tauwf(1:nw0, i) = fozwf(1:nw0, i) * ozs(i) / hresgabs(radcidxs(1:nw0), i)
  !ENDDO

  ! xliu, 11/02/2011, the above is incorrect as hresgabs is the total absorption
  ! It should be as follows by using the ozone absorption
  DO i = 1, nz
     tauwf(1:nw0, i) = fozwf(1:nw0, i) / o3crsz(radcidxs(1:nw0), i) / du2mol
  ENDDO

  c2hfidx(1) = 1; fidxc = 1
  DO iwin = 1, numwin
     IF (iwin == numwin) THEN
        lidx = nhw; lidxc = nw0
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
     fidx = c2hfidx(iw); lidx = c2hlidx(iw)
     hrad(fidx:lidx) = rad(iw)
     DO i = 1, nz
        dtau(fidx:lidx) = hresgabs(fidx:lidx, i) - hresgabs(radcidxs(iw), i)
        dray(fidx:lidx) = hresray(fidx:lidx, i)  - hresray(radcidxs(iw), i)
        hrad(fidx:lidx) = hrad(fidx:lidx) + tauwf(iw, i) * dtau(fidx:lidx) + &
             fraywf(iw, i) * dray(fidx:lidx)

        ! Is it better to assume same wf (after normalized by o3 xsec)
        hozwf(fidx:lidx, i) = fozwf(iw, i) * o3crsz(fidx:lidx, i) / o3crsz(radcidxs(iw), i)
     ENDDO
  ENDDO
  
  !DO i = 1, nz 
  !   CALL BSPLINE(wave(1:nw0), fozwf(1:nw0, i), nw0, hwave(1:nhw), hozwf(1:nhw, i), nhw, errstat)
  !   IF (errstat < 0) THEN
  !      WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
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
     sidx = 2 + i;        inarr(1:nhw, sidx) = hres_i0(1:nhw) * hozwf(1:nhw, i)
  ENDDO
  IF (do_albwf) THEN
     sidx = sidx + 1;     inarr(1:nhw, sidx) = hres_i0(1:nhw) * halbwf(1:nhw)
  ENDIF
  IF (do_cfracwf) THEN
     sidx = sidx + 1;     inarr(1:nhw, sidx) = hres_i0(1:nhw) * hcfracwf(1:nhw)
  ENDIF
  IF (do_faerwf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hres_i0(1:nhw) * haerwf(1:nhw, i)
     ENDDO
  ENDIF
  IF (do_faerswf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hres_i0(1:nhw) * haerswf(1:nhw, i)
     ENDDO
  ENDIF
  IF (do_codwf) THEN
     DO i = nctp, ncbp
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hres_i0(1:nhw) * hcodwf(1:nhw, i)
     ENDDO
  ENDIF
  IF (do_sprswf) THEN
     DO i = nsprs, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hres_i0(1:nhw) * hsprswf(1:nhw, i)
     ENDDO
  ENDIF

  ! convolve ozone shift/temperature 
  IF (do_tracewf) THEN
     DO i = 1, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = o3crsz(1:nhw, i) * hres_i0(1:nhw) 
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
           
           !IF (fgassidxs(i) > 0) THEN
           ! sidx = sidx + 1;  inarr(1:nw, sidx) = hres_gasshi(i, 1:nhw)
           !ENDIF
        ENDIF
     ENDDO

     IF (use_so2dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = so2crsz(1:nhw, i) * hres_i0(1:nhw) 
        ENDDO
        !IF (do_so2shi) THEN
        !   sidx = sidx + 1;  inarr(1:nhw, sidx) = so2dads(1:nhw) 
        !ENDIF
     ENDIF

     do_o4shi = .FALSE. ; do_o2shi = .FALSE. ; do_h2oshi = .FALSE.
     DO i = 1, ngas
        IF (fgasidxs(i) > 0) THEN
           IF (gasidxs(i) == o2o2_idx ) THEN
              IF (fgassidxs(i) > 0 ) do_o4shi = .TRUE.
              IF (use_o4dtcrs) CYCLE
           ENDIF
           IF (gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx ) THEN
              IF (fgassidxs(i) > 0 ) do_o2shi = .TRUE.
           ENDIF
           IF (gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx ) THEN
              IF (fgassidxs(i) > 0 ) do_o2shi = .TRUE.
           ENDIF
        ENDIF       
     ENDDO

     IF (use_o4dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = o4crsz(1:nhw, i) * hres_i0(1:nhw) 
        ENDDO
        !DO i = 1, nhw 
        !   print * ,hwave(i),  o4crsz(i, 20)
        !ENDDO
         
     ENDIF
     IF (use_o2dptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = o2crsz(1:nhw, i)  * hres_i0(1:nhw) 
        ENDDO
        !DO i = 1, nhw 
        !   print * ,hwave(i),  o4crsz(i, 20)
        !ENDDO
     ENDIF
     IF (use_h2odptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = h2ocrsz(1:nhw, i) !* hres_i0(1:nhw) 
        ENDDO
     ENDIF
  ENDIF

  IF (do_o3shi) THEN
     DO i = 1, nz
        sidx = sidx + 1;   inarr(1:nhw, sidx) = o3dadsz(1:nhw, i) !* hres_i0(1:nhw)
     ENDDO
  ENDIF
  IF (do_tmpwf) THEN
     DO i = 1, nz
        sidx = sidx + 1;   inarr(1:nhw, sidx) = o3dadtz(1:nhw, i) !* hres_i0(1:nhw)
     ENDDO
  ENDIF
  eidx = sidx

  ! *** second, convole all spectra at once ****
  CALL convol_f2c(hwave(1:nhw), inarr(1:nhw, 1:eidx), nhw, eidx, owave(1:now), outarr(1:now, 1:eidx), now)
   
  ! *** third, transfer all convolved spectra back ***
  oi0(1:now)  = outarr(1:now, 1)
  hrad(1:now) = outarr(1:now, 2) / oi0(1:now)

  IF (do_pslwf .and. which_pslwf ==1  ) THEN 
    solwinfit_save = solwinfit
    do_dsdw = .false. ; do_dsdk = .false.
    DO i = 1, npsl
     idx = psl_fpos(i)
     solwinfit(1:numwin, idx, 1) = solwinfit_save(1:numwin, idx,1) *1.001
     !CALL convol_f2c(hwave(1:nhw), inarr(1:nhw,2)/inarr(1:nhw,1), nhw, 1, owave(1:now), outarr1(1:now,2), now)
     !hrad1(1:now) = outarr1(1:now, 2) !/ outarr1(1:now,1) 
     CALL convol_f2c(hwave(1:nhw), inarr(1:nhw,1:2), nhw, 2, owave(1:now), outarr1(1:now, 1:2), now)
     hrad1(1:now) = outarr1(1:now, 2)/ outarr1(1:now,1) 
     dpabs(1:now) = 0.0
     fidx = 1
     DO iwin = 1, numwin
       lidx = fidx + onpix(iwin) -1
       dpabs(fidx:lidx) =  solwinfit_save(iwin,idx,1)*0.001
       fidx = lidx+1
     ENDDO
     database_pslwf(i, refidx(1:now)) =(hrad1(1:now)-hrad(1:now))/(dpabs(1:now)) /hrad1(1:now)
     solwinfit = solwinfit_save 
    ENDDO
  ENDIF 
  
  IF (do_pslwf .and. which_pslwf == 2) THEN 

    ! this show a better fitting acurrcy
    DO i = 1, npsl 
      IF ( psl_fpos(i) == hwe_idx) THEN  
         do_dsdw = .true. ; do_dsdk = .false.
      ELSE IF (psl_fpos(i) == spk_idx) THEN 
         do_dsdw = .false. ; do_dsdk = .true.
      ENDIF
      !CALL convol_f2c(hwave(1:nhw), inarr(1:nhw, 2)/inarr(1:nhw,1), nhw, 1, &
      !    owave(1:now), outarr1(1:now, 2), now)
      !hrad1(1:now) = outarr1(1:now, 2)
      CALL convol_f2c(hwave(1:nhw), inarr(1:nhw,1:2), nhw, 2, owave(1:now), outarr1(1:now,1:2), now)
      hrad1(1:now) = outarr1(1:now, 2)/ outarr1(1:now,1) 
      database_pslwf(i,refidx(1:now)) = hrad1(1:now)/hrad(1:now)
    ENDDO
    do_dsdw = .false. ; do_dsdk = .false.
  ENDIF


  DO i = 1, nz
     sidx = 2 + i;        hozwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
  ENDDO
  IF (do_albwf) THEN
     sidx = sidx + 1;     halbwf(1:now) = outarr(1:now, sidx) / oi0(1:now)
  ENDIF
  IF (do_cfracwf) THEN
     sidx = sidx + 1;     hcfracwf(1:now) = outarr(1:now, sidx) / oi0(1:now)
  ENDIF
  IF (do_faerwf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  haerwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF
  IF (do_faerswf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  haerswf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF
  IF (do_codwf) THEN
     DO i = nctp, ncbp
        sidx = sidx + 1;  hcodwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF
  IF (do_sprswf) THEN
     DO i = nsprs, nz
        sidx = sidx + 1;  hsprswf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF

  ! convolve ozone shift/temperature 
  IF (do_tracewf) THEN
     DO i = 1, nz
        sidx = sidx + 1;  abscrs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
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
           
           !IF (fgassidxs(i) > 0) THEN
           !   sidx = sidx + 1;  tmp_gasshi(i, 1:now) = outarr(1:now, sidx) 
           !ENDIF
        ENDIF
     ENDDO

     IF (use_so2dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  so2crs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO

        !IF (do_so2shi) THEN
        !   sidx = sidx + 1;  so2dads1(1:now) = outarr(1:now, sidx) 
        !ENDIF
     ENDIF

     do_o4shi = .FALSE.
     DO i = 1, ngas
        IF (fgasidxs(i) > 0) THEN
           IF (gasidxs(i) == o2o2_idx) THEN
              IF (fgassidxs(i) > 0 ) do_o4shi = .TRUE.
           ELSE IF (gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx) THEN 
              IF (fgassidxs(i) > 0 ) do_o2shi = .TRUE.
           ELSE IF (gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx) THEN 
              IF (fgassidxs(i) > 0 ) do_h2oshi = .TRUE.
           ENDIF
        ENDIF
     ENDDO

     IF (use_o4dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  o4crs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO 
        !print * , o4crs1(now-10:now, 1),'now' ; stop
     ENDIF
     IF (use_o2dptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  o2crs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO
     ENDIF
     IF (use_h2odptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  h2ocrs1(1:now, i)  = outarr(1:now, sidx)  !/ oi0(1:now)
        ENDDO
     ENDIF

  ENDIF

  IF (do_o3shi) THEN
     DO i = 1, nz
        sidx = sidx + 1;  dads1(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
     ENDDO
  ENDIF
  IF (do_tmpwf) THEN
     DO i = 1, nz
        sidx = sidx + 1;  dadt1(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
     ENDDO
  ENDIF
   
  ! Perform additional coadding (if necesary)
  ! Use OMI solar spectra here
  IF (do_bandavg) THEN
     tmpi0(1:now) = i0sav(refidx_sav(1:now))  ! better to use OMI solar spectra     
     oi0(1:now) = tmpi0
     CALL avg_band_spec(owave(1:now), oi0(1:now), now, ntemp, errstat)
     IF (ntemp /= nrad .OR. errstat /= 0) THEN
        WRITE(*, *) 'Spectra Averaging Error: ', now, ntemp, nrad 
        errstat = pge_errstat_error; RETURN
     ENDIF

     otmp(1:now) = hrad(1:now) * tmpi0(1:now)
     CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
     rad(1:nrad) = otmp(1:nrad) / oi0(1:nrad)

     IF (do_pslwf) THEN
        print * , 'not implemented '; stop
     ENDIF

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
        IF (use_o4dtcrs) THEN 
           DO i = 1, nz
              otmp(1:now) = o4crs1(1:now, i) * tmpi0(1:now) 
              CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
              o4crs(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
           ENDDO
        ENDIF
        IF (use_o2dptcrs) THEN 
           DO i = 1, nz
              otmp(1:now) = o2crs1(1:now, i) * tmpi0(1:now) 
              CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
              o2crs(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
           ENDDO
        ENDIF
        IF (use_h2odptcrs) THEN 
           DO i = 1, nz
              otmp(1:now) = h2ocrs1(1:now, i) !* tmpi0(1:now) 
              CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
              h2ocrs(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
           ENDDO
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
        IF (use_o4dtcrs) THEN
           o4crs(1:now, 1:nz) = o4crs1(1:now, 1:nz)
        ENDIF
        IF (use_o2dptcrs) THEN
           o2crs(1:now, 1:nz) = o2crs1(1:now, 1:nz)
        ENDIF
        IF (use_h2odptcrs) THEN
           h2ocrs(1:now, 1:nz) = h2ocrs1(1:now, 1:nz)
        ENDIF
     ENDIF
  ENDIF
  RETURN

END SUBROUTINE hres_radwf_inter_convol

SUBROUTINE radwf_interpol(nw, nz, nctp, ncbp, nsprs, faerlvl, do_radcals, &
     do_fozwf, do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf, wave, abscrs, &
     ozs, rad, fozwf, albwf, cfracwf, faerwf, faerswf, fcodwf, fsprswf, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY  : numwin, nradpix, band_selectors 
  USE ozprof_data_module,     ONLY  : nup2p 
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
  effcrs = 0.0; toz = SUM(ozs)
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

        nd = 0; nud = 0
        DO i = fidx, lidx
           IF ( do_radcals(i) ) THEN
              nd = nd + 1;   didxs(nd)  = i
           ELSE
              nud = nud + 1; uidxs(nud) = i
           ENDIF
        ENDDO
        
       ! for radiances (lnI vs. ozcrs)
       k = MIN(uidxs(1)-2, 1); nd1 = nd - k + 1
       crs1(1:nd1)  = effcrs(didxs(k:nd));   tmp1(1:nd1)  = rad(didxs(k:nd))
       crs2(1:nud) = effcrs(uidxs(1:nud))

       CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
       IF (errstat < 0) THEN
          WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
          errstat = pge_errstat_error; RETURN
       ENDIF

       rad(uidxs(1:nud)) =  tmp2(1:nud)

       ! for albedo weighting funtions: ln(albwf/rad) vs ozcrs
       IF (do_albwf .AND. MAXVAL(albwf) > 0.0) THEN
          tmp1(1:nd1) = LOG(albwf(didxs(k:nd)))
 
          CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
          IF (errstat < 0) THEN
             WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
             errstat = pge_errstat_error; RETURN
          ENDIF

          albwf(uidxs(1:nud)) =  EXP(tmp2(1:nud))
       ENDIF

       IF (do_cfracwf .AND. MAXVAL(cfracwf) > 0.0) THEN
          tmp1(1:nd1) = LOG(cfracwf(didxs(k:nd)))
 
          CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
          IF (errstat < 0) THEN
             WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
             errstat = pge_errstat_error; RETURN
          ENDIF

          cfracwf(uidxs(1:nud)) =  EXP(tmp2(1:nud))
       ENDIF

       ! for ozone weighting function (fine grids) ln(ozwf/rad) vs ozcrs
       IF (do_fozwf) THEN
          DO j = 1, nz 
             crs1(1:nd1)  = abscrs(didxs(k:nd), j);   tmp1(1:nd1) = LOG(fozwf(didxs(k:nd), j))
             crs2(1:nud) = abscrs(uidxs(1:nud), j)

             CALL BSPLINE(crs1(1:nd1), tmp1(1:nd1), nd1, crs2(1:nud), tmp2(1:nud), nud, errstat)
             IF (errstat < 0) THEN
                WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
                errstat = pge_errstat_error; RETURN
             ENDIF

             fozwf(uidxs(1:nud), j) = EXP(tmp2(1:nud))
          ENDDO
       ENDIF

       ! Linear interpolation over wavelength
       IF (do_faerwf) THEN
          wav1(1:nd1) = wave(didxs(k:nd));  wav2(1:nud) = wave(uidxs(1:nud))
          DO j = faerlvl, nz 
             tmp1(1:nd1) = faerwf(didxs(k:nd), j)

             CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
             IF (errstat < 0) THEN
                WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
                errstat = pge_errstat_error; RETURN
             ENDIF

             faerwf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
       ENDIF

       IF (do_faerswf) THEN
          wav1(1:nd1) = wave(didxs(k:nd));  wav2(1:nud) = wave(uidxs(1:nud))
          DO j = faerlvl, nz 
             tmp1(1:nd1) = faerswf(didxs(k:nd), j)

             CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
             IF (errstat < 0) THEN
                WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
                errstat = pge_errstat_error; RETURN
             ENDIF

             faerswf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
       ENDIF

       IF (do_codwf) THEN
          wav1(1:nd1) = wave(didxs(k:nd));  wav2(1:nud) = wave(uidxs(1:nud))
          DO j = nctp, ncbp 
             tmp1(1:nd1) = fcodwf(didxs(k:nd), j)

             CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
             IF (errstat < 0) THEN
                WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
                errstat = pge_errstat_error; RETURN
             ENDIF

             fcodwf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
       ENDIF

       IF (do_sprswf) THEN
          wav1(1:nd1) = wave(didxs(k:nd));  wav2(1:nud) = wave(uidxs(1:nud))
          DO j = nsprs, nz 
             tmp1(1:nd1) = fsprswf(didxs(k:nd), j)

             CALL BSPLINE(wav1(1:nd1), tmp1(1:nd1), nd1, wav2(1:nud), tmp2(1:nud), nud, errstat)
             IF (errstat < 0) THEN
                WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
                errstat = pge_errstat_error; RETURN
             ENDIF

             fsprswf(uidxs(1:nud), j) = tmp2(1:nud)
          ENDDO
       ENDIF

    ELSE IF (band_selectors(iw) == 2) THEN ! Channel 2, linear interpolation vs cross section

       ! p1,  p2:  points without radiative calculation
       ! dp1, dp2: points with radiative calculation
       dp1 = fidx; i = dp1 + 1
       DO WHILE (i <= lidx)
          IF ( do_radcals(i) ) THEN
             dp2 = i
             
             IF (dp2 - dp1 > 1) THEN  ! there are points with no rad/wf in between
                ! get rad/wf in between using rad/wf at dp1:dp2
                df = effcrs(dp2) - effcrs(dp1)
                p1 = dp1 + 1; p2 = dp2 - 1
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
       WRITE(*, *) 'Interpolation for this channel is not implemented!!!'
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


SUBROUTINE polcorr_online(niter, which_polcorr, nw,  nz, nctp, ncbp, nsprs, &
           faerlvl,   ncorr, polidxs, &
   do_fozwf, do_albwf, do_faerwf, do_faerswf,do_codwf, do_sprswf, do_fraywf, do_cfracwf, &
   wave,  rad, prad, tabs, o3crs, ozs, &
     albwf, palbwf,   fozwf, pfozwf,    faerwf, pfaerwf, faerswf, pfaerswf, &
     fcodwf,pfcodwf, fsprswf, pfsprswf, fraywf, pfraywf, cfracwf, pcfracwf)

  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : maxlay, mflay, du2mol
  USE OMSAO_variables_module,  ONLY : currloop
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN) :: nw, nz, ncorr, niter, which_polcorr, nctp, ncbp, faerlvl, nsprs
  LOGICAL, INTENT(IN) :: do_fozwf, do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_fraywf, do_cfracwf

  INTEGER, DIMENSION (ncorr),          INTENT (IN) :: polidxs
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fozwf, faerwf,   faerswf, fcodwf, fsprswf, fraywf
  REAL (KIND=dp), DIMENSION(nw, nz),    INTENT(IN) :: pfozwf, pfaerwf, pfaerswf, pfcodwf, pfsprswf, tabs, o3crs, pfraywf
  REAL (KIND=dp), DIMENSION(nw),        INTENT(IN) :: wave
  REAL (KIND=dp), DIMENSION(nz),        INTENT(IN) :: ozs
  REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: rad, albwf, cfracwf
  REAL (KIND=dp), DIMENSION(nw),        INTENT(IN) :: prad, palbwf, pcfracwf
  


  ! Local variables
  INTEGER, PARAMETER                               :: mcorr = 20, mw = 500
  INTEGER                                          :: i, j, fidx, lidx
  REAL (KIND=dp)                                   :: frac, frac1, frac2, tmprad
  REAL (KIND=dp), DIMENSION(nz)                    :: tmpcorr
  REAL (KIND=dp), DIMENSION(mcorr, mflay), SAVE    :: tauwf, ptauwf
  REAL (KIND=dp), DIMENSION(mcorr),        SAVE    :: dfrad, dfalbwf, dfcfracwf
  REAL (KIND=dp), DIMENSION(mcorr, mflay), SAVE    :: dftauwf
  REAL (KIND=dp), DIMENSION(mcorr, mflay), SAVE    :: dffozwf, dffaerwf, dffaerswf, dffcodwf, dffsprswf, dffraywf
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
!           ptauwf(i, 1:nz) = pfozwf(polidxs(i), 1:nz) * ozs / tabs(polidxs(i), 1:nz) / prad(polidxs(i))
!           tauwf (i, 1:nz) = fozwf (polidxs(i), 1:nz) * ozs / tabs(polidxs(i), 1:nz) / rad (polidxs(i))
           ! xliu, 11/02/2011, the above is incorrect as tabs is the total absorption
           ! It should be as follows by using ozone absorption
           ptauwf(i, 1:nz) = pfozwf(polidxs(i), 1:nz) / o3crs(polidxs(i), 1:nz) / prad(polidxs(i)) / du2mol
           tauwf (i, 1:nz) = fozwf (polidxs(i), 1:nz) / o3crs(polidxs(i), 1:nz) / rad (polidxs(i)) / du2mol
        ENDDO

        dftauwf(1:ncorr, 1:nz) = ptauwf(1:ncorr, 1:nz) - tauwf(1:ncorr, 1:nz)
!        WRITE(*, *) ' ***'
!        print *, prad(polidxs(4)), rad(polidxs(4))
!        WRITE(*, '(10D12.4)') dftauwf(4, 1:nz)      
!        WRITE(*, '(10D12.4)') ptauwf(4, 1:nz)
!        WRITE(*, '(10D12.4)') tauwf(4, 1:nz)
!        WRITE(*, '(10D12.4)') pfozwf(polidxs(4), 1:nz)
!        WRITE(*, '(10D12.4)') fozwf(polidxs(4), 1:nz)
!        WRITE(*, *) ' ***'
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
  fidx = 1; lidx=polidxs(1)
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

! Obtain minor trace gas weighting functions from ozwf
! Note: it is not the exact wf but negative of the WF divided by radiances
!SUBROUTINE GET_TRACEGAS_WF (ozwf, ozabs, so2crs, use_so2dtcrs, rad, nw, nz, nz1, &
!     ozs, waves, do_so2zwf, so2zwf)
SUBROUTINE get_tracegas_wf (nw, nz, nz1, rad, ozwf, ozabs, & 
  use_so2dtcrs, so2crs, use_o4dtcrs, o4crs, use_o2dptcrs, o2crs, use_h2odptcrs,h2ocrs,do_so2zwf, so2zwf)
  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : du2mol
  USE OMSAO_indices_module,   ONLY : so2_idx, so2v_idx, bro_idx, o2o2_idx, o2_idx, h2o_idx,  &
                                o2t2_idx, h2ot2_idx, hcho_idx, no2_t1_idx, ring_idx, ring1_idx
  USE OMSAO_variables_module, ONLY : mask_fitvar_rad, refidx, database, database_save, refspec_norm
  USE ozprof_data_module,     ONLY : fps, fzs, nlay, mgasprof, fgasidxs, &
       tracegas, ngas, gasidxs, fgassidxs, so2valts, so2vprofn1p1, trace_profwf, nup2p, use_lograd
  IMPLICIT NONE

  ! Input/Output variables
  INTEGER, INTENT(IN)                           :: nw, nz, nz1
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(IN) :: ozwf, ozabs, so2crs, o4crs,o2crs, h2ocrs
  REAL (KIND=dp), DIMENSION(nw),     INTENT(IN) :: rad !, waves
!  REAL (KIND=dp), DIMENSION(nz),     INTENT(IN) :: ozs
  LOGICAL,                           INTENT(IN) :: use_so2dtcrs, do_so2zwf, use_o4dtcrs, &
  use_o2dptcrs, use_h2odptcrs
  REAL (KIND=dp), DIMENSION(nw),    INTENT(OUT) :: so2zwf
  ! Local variables
  INTEGER                                       :: i, j, k, fidx, lidx, nk
  REAL (KIND=dp)                                :: tmp
  REAL (KIND=dp), DIMENSION(nw, nz)             :: amf, tmpcrs
  REAL (KIND=dp), DIMENSION(ngas, nw)           :: tamf
  REAL (KIND=dp)                                :: avcd, svcd
 
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
        IF ( ((gasidxs(i) /= so2_idx .AND. gasidxs(i) /= so2v_idx)  .OR. .NOT. use_so2dtcrs) .AND. &
             (gasidxs(i) /= o2o2_idx  .OR. .NOT. use_o4dtcrs) .AND. &
             (gasidxs(i) /= h2o_idx   .OR. .NOT. use_h2odptcrs) .AND. &
             (gasidxs(i) /= o2_idx    .OR. .NOT. use_o2dptcrs) ) THEN
           DO j = 1, nw
              tamf(i, j) = SUM(amf(j, 1:nz1) * mgasprof(i, 1:nz1)) / avcd          
           ENDDO
           
           IF (fgassidxs(i) > 0) THEN
              database(gasidxs(i), refidx(1:nw)) = database(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
           ELSE
              database(gasidxs(i), refidx(1:nw)) = database_save(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
              !database(gasidxs(i), refidx(1:nw)) = database(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
           ENDIF
          
           tracegas(i, 7) = 0.0; nk = 0
           DO j = 1, nw
              IF (database_save(gasidxs(i), refidx(j)) > 0.) THEN
                 tracegas(i, 7) = tracegas(i, 7) + tamf(i, j)
                 nk = nk + 1
              ENDIF
           ENDDO
           tracegas(i, 7) = tracegas(i, 7) / nk

           ! xliu: 08/06/2010, Add trace gas profile weighting function (dY/dx)
           ! Note dy = dI, instead of DlnI
           ! x is the normalized (by refspec_norm) quantity at each ozone retrieval layer 
           ! instead of VLIDORT computation layer
           tmp = 0.0
           DO j = 1, nlay
              fidx = nup2p(j - 1) + 1; lidx = nup2p(j)
              svcd = SUM(mgasprof(i, fidx:lidx))
              IF (svcd > 0.0d0) THEN
                 trace_profwf(i, 1:nw, j) = amf(:, fidx) * mgasprof(i, fidx)
                 DO k = fidx + 1, lidx
                    trace_profwf(i, 1:nw, j) = trace_profwf(i, 1:nw, j) + amf(:, k) * mgasprof(i, k)
                 ENDDO
                 IF (.NOT. use_lograd) trace_profwf(i, 1:nw, j) = trace_profwf(i, 1:nw, j) * rad(1:nw) 
                 IF (fgassidxs(i) > 0) THEN
                    trace_profwf(i, 1:nw, j) = trace_profwf(i, 1:nw, j) * database(gasidxs(i), refidx(1:nw)) 
                 ELSE
                    trace_profwf(i, 1:nw, j) = trace_profwf(i, 1:nw, j) * database_save(gasidxs(i), refidx(1:nw))
                 ENDIF
                 trace_profwf(i, 1:nw, j) = -trace_profwf(i, 1:nw, j) / svcd
               ELSE
                  trace_profwf(i, 1:nw, j) = 0.0d0
               ENDIF
               tmp = tmp + trace_profwf(i, 50, j) * svcd / avcd

               !IF ( i == 5) then
               !   WRITE(*, '(3I5, 2D16.5)') j, fidx, lidx, avcd, svcd
               !   print *, trace_profwf(i, 1, j), trace_profwf(i, nw, j)
               !   print *, amf(1, j), amf(nw, j)
               !   print *, database(gasidxs(i), refidx(1)), database(gasidxs(i), refidx(nw))
               !ENDIF
           ENDDO           
        ELSE
           tmp = avcd * refspec_norm(gasidxs(i))
           IF ( gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN 
               tmpcrs = so2crs
           ELSE IF (gasidxs(i) == o2o2_idx) THEN 
               tmpcrs = o4crs
           ELSE IF (gasidxs(i) == o2_idx ) THEN 
               tmpcrs = o2crs
           ELSE IF (gasidxs(i) == h2o_idx) THEN
               tmpcrs = h2ocrs
           ENDIF
           DO j = 1, nw
              tamf(i, j) = SUM(amf(j, 1:nz1) * mgasprof(i, 1:nz1) * tmpcrs(j, 1:nz1) ) / tmp
           ENDDO
           database(gasidxs(i), refidx(1:nw)) = tamf(i, 1:nw) 
           tracegas(i, 7) = 0.0; nk = 0
           !if (gasidxs(i) == o2o2_idx )  then 
           !        print * , database(gasidxs(i),refidx(1:nw))
           !ENDIF
           DO j = 1, nw
              IF (database_save(gasidxs(i), refidx(j)) > 0.) THEN
                 tracegas(i, 7) = tracegas(i, 7) + tamf(i, j) / database_save(gasidxs(i), refidx(j)) 
                 nk = nk + 1
              ENDIF
           ENDDO
           IF (nk > 0 ) THEN 
              tracegas(i, 7) = tracegas(i, 7) / nk
           ENDIF
           ! xliu: 08/06/2010, Add trace gas profile weighting function (dY/dx)
           ! Note dy = dI, instead of DlnI
           ! x is the normalized (by refspec_norm) quantity at each ozone retrieval layer 
           ! instead of VLIDORT computation layer
           DO j = 1, nlay
              fidx = nup2p(j - 1) + 1; lidx = nup2p(j)
              svcd = SUM(mgasprof(i, fidx:lidx))
              IF (svcd > 0.0d0) THEN
                 trace_profwf(i, 1:nw, j) = mgasprof(i, fidx) * amf(:, fidx) * tmpcrs(1:nw, fidx)
                 DO k = fidx + 1, lidx
                    trace_profwf(i, 1:nw, j) = trace_profwf(i, 1:nw, j) + mgasprof(i, k) * amf(:, k) * tmpcrs(1:nw, k)
                 ENDDO
                 IF (.NOT. use_lograd) trace_profwf(i, 1:nw, j) = trace_profwf(i, 1:nw, j) * rad(1:nw) 
                 trace_profwf(i, 1:nw, j) = -trace_profwf(i, 1:nw, j) / svcd
              ELSE
                 trace_profwf(i, 1:nw, j) = 0.0d0
              ENDIF
           ENDDO
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
END SUBROUTINE get_tracegas_wf

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
     lidx = fidx + 9; dz = (zs(i) - zs(i-1)) / 10.
     DO j = fidx, lidx
        zs1(j) = zs(i-1) + dz * (j - fidx + 1)
     ENDDO
     fidx  = lidx + 1
  ENDDO

  !WRITE(*, '(12F8.3)') zs(0:nz)
  !WRITE(*, *)
  !WRITE(*, '(12F8.3)') zs1(0:nz1)

  CALL BSPLINE(zs, fts, nz+1, zs1(0:nz1), ts1(0:nz1), nz1+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ' : BSPLINE error, errstat = ', errstat; RETURN
  ENDIF

  CALL BSPLINE(zs, cumo3, nz+1, zs1(0:nz1), cumo31(0:nz1), nz1+1, errstat)
  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ' : BSPLINE error, errstat = ', errstat; RETURN
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
END MODULE m_lidort_util
