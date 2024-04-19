!
MODULE m_lidort_util

  USE m_ezspline_interpolation, only: bspline, bspline2, interpol, interpol2, interpolation
  USE m_avg_band, ONLY: avg_band_spec
  USE m_convol, ONLY: convol_f2c, get_i0
  !USE OMSAO_gome_data_module, ONLY  : n_gome_q, gome_q 

  PUBLIC get_hres_radcal_waves, &
         hres_radwf_inter_convol, radwf_inter_convol, radwf_interpol, & 
         polcorr_online,polcorr_online_with_lut, set_polcorr, &
         get_slant_tau, get_efft, get_tracegas_wf, debug_rtm, debug_taug

  INTEGER, PARAMETER, PRIVATE :: max_pathlen = 1024
  INTEGER, PARAMETER, PRIVATE :: n_gome_q = 1
  REAL, DIMENSION (3,n_gome_q), PRIVATE :: gome_q

CONTAINS

  !1.	Establish fine wavelength grid: 0.01 nm now, may change to 0.05 nm later
  !2.	Establish radiance calculation grid, based on spectral sampling rate for 
  !   different specral regions
  !3.	Find indices of radiance calculation grid in fine wavelength grid

SUBROUTINE get_hres_radcal_waves(errstat)

  !1.	Establish fine wavelength grid: 0.01 nm now, may change to 0.05 nm later
  !2.	Establish radiance calculation grid, based on spectral sampling rate for 
  !   different specral regions
  !3.	Find indices of radiance calculation grid in fine wavelength grid

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY  : numwin, nradpix, solwinfit, which_slit, &
      curr_rad_spec, use_redfixwav, winlim, nuvwin, nviswin, widx_hvis,widx_rvis, & 
      wcenter_uvvis, numwin, winwav_max, winwav_min
  USE OMSAO_indices_module,   ONLY  : hwe_idx, wvl_idx
  USE ozprof_data_module,     ONLY  : radc_nsegsr, radc_samprate, radc_lambnd,  &
        hres_samprate, hres_vis_samprate, hreswav, hres_i0, radcwav, nhresp, ncalcp, radcidxs, & 
        npca, npcapix, which_pcabin, winpca, &
        nalb, albmin, albmax
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! Output variables
  INTEGER, INTENT(OUT) :: errstat

  ! Local variables
  REAL (KIND=dp), PARAMETER     :: dhw0_uv = 0.01, dhw0_vis=0.001 ! I0 grids
  INTEGER                       :: i, j, k, fidx, lidx, nsub, nratio, nhalf, ntmp, maxline, nhresp0  !nextra, n0
  INTEGER, DIMENSION (numwin,2) :: nhpix
  REAL (KIND=dp)       :: tmp, swav, ewav, slw, samprate, hsamprate,invdhw, dhw0
  REAL (KIND=dp), DIMENSION (:), allocatable :: tmp_i0,tmp_wav
  REAL (KIND=dp), DIMENSION (:), allocatable :: hres0_i0, hres0_wav
  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=21), PARAMETER :: modulename = 'get_hres_radcal_waves'

  errstat = pge_errstat_ok
  WRITE(www_lun,*) '******'//modulename//'*****'
  !**  Establish high resolution grid 
  maxline = CEILING((winwav_max - winwav_min)/dhw0_vis) + 1
  allocate ( tmp_wav(maxline), tmp_i0(maxline))
  nhresp = 0   
  IF (.NOT. use_redfixwav) THEN
     fidx = 1
     DO i = 1, numwin
        dhw0 = dhw0_uv; IF (i > nuvwin) dhw0 = dhw0_vis
        invdhw = 1.0 / dhw0
        swav = FLOOR(invdhw * (winlim(i, 1))) / invdhw
        IF (i > 1) THEN
           IF (swav < ewav) swav = ewav + dhw0
        ENDIF
        ewav = CEILING(invdhw * (winlim(i, 2))) / invdhw         
        nsub = NINT((ewav - swav) / dhw0) + 1  
        tmp_wav(nhresp+1:nhresp+nsub) = swav + dhw0 * (/(j, j = 0, nsub-1)/)      
        nhresp = nhresp + nsub
        nhpix(i,1) = fidx ; nhpix(i,2) = nhresp
        fidx = nhresp + 1  
     ENDDO
  ELSE 
     ! If use fixed wavelengths, then go through it one by one
     fidx = 1
     DO i = 1, numwin
        dhw0 = dhw0_uv; IF (i > nuvwin) dhw0 = dhw0_vis
        invdhw = 1.0 / dhw0
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
           tmp_wav(nhresp+1:nhresp+nsub) = swav + dhw0 * (/(j, j = 0, nsub-1)/)
           nhresp = nhresp + nsub
        ENDDO
        fidx = lidx + 1
     ENDDO
  ENDIF
  
  ! ** define wav & i0 @ original highresolution (0.01 nm in UV)
  nhresp0 = nhresp
  allocate (hres0_wav(nhresp0), hres0_i0(nhresp0))
  hres0_wav(1:nhresp0) = tmp_wav(1:nhresp0)
  call get_i0(nhresp0, hres0_wav, hres0_i0)
  ! ** re-grided wavelengths (.e.g 0.05 nm in UV)
  J =1
  DO k = 1, numwin
    dhw0 = dhw0_uv ; hsamprate = hres_samprate
    IF (k > nuvwin) THEN 
      dhw0 = dhw0_vis ; hsamprate = hres_vis_samprate
    ENDIF
    invdhw = 1.0 /hsamprate
    nratio = NINT(hsamprate / dhw0) !RECHECK
    nhalf =  nratio / 2
    DO i = nhpix(k,1) +nhalf , nhpix(k,2) -nhalf, nratio
       fidx = i - nhalf; lidx = i - nhalf + nratio - 1
       IF (lidx > nhpix(k,2)) THEN
           lidx = nhpix(k,2)
       ENDIF
       ntmp = lidx - fidx + 1
       tmp_wav(j)= SUM(hres0_wav(fidx:lidx))/ REAL(ntmp)
       tmp_i0(j) = SUM(hres0_i0(fidx:lidx)) / REAL(ntmp)
       IF (ntmp /= nratio) THEN
           j = j -1
           WRITE(*,*) 'undersampleing' , ntmp, 'in get_pcahres_wave'
       ENDIF
      j = j + 1
    ENDDO
  ENDDO
  nhresp = j -1
  allocate (hreswav(nhresp), hres_i0(nhresp))
  hreswav(1:nhresp) = tmp_wav(1:nhresp)
  hres_i0(1:nhresp) = tmp_i0(1:nhresp)

  ! position bwt UV and VIS in hreswav
  widx_hvis = nhresp + 1
  IF ( nviswin > 0) THEN
     widx_hvis = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) <  wcenter_uvvis))) +1
  ELSE
     widx_hvis = 0
  ENDIF

  ! Establish radiance calculation grid, based on spectral sampling
  ! Sampling rate are specified in several segments (e.g., < 295 nm, 295-308 nm, > 308 nm)
  ! Make sure that radiance will be done for the first and last points.  
  ncalcp = 1; tmp_wav(1) = hreswav(1)
  fidx = 2
  DO i = 1, radc_nsegsr
     IF (i == radc_nsegsr) THEN
        lidx = nhresp - 1
     ELSE
        lidx = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) < radc_lambnd(i+1))))
     IF (lidx == 0) CYCLE
     ENDIF
     hsamprate = hres_samprate
     IF ( fidx >= widx_hvis ) hsamprate = hres_vis_samprate
     invdhw = 1.0 /hsamprate
     samprate = FLOOR(radc_samprate(i) * invdhw) / invdhw  ! multiples of dhw
     DO j = fidx, lidx
        IF (j >= nhresp) cycle
        tmp = ABS(hreswav(j) - tmp_wav(ncalcp) - samprate)
        IF (tmp < hsamprate * 0.1) THEN
           ncalcp = ncalcp + 1
           tmp_wav(ncalcp) = hreswav(j)
        ELSE IF (hreswav(j + 1) >= tmp_wav(ncalcp) + samprate * 2 &
             .AND. hreswav(j)   >  tmp_wav(ncalcp) + samprate * 0.5) THEN
           ncalcp = ncalcp + 1
           tmp_wav(ncalcp) = hreswav(j)
        ENDIF
     ENDDO
     WRITE(www_lun,'(A,i3,10f8.3)') 'radcwav:',i, samprate, hsamprate, hreswav(fidx), hreswav(lidx) 
     IF (lidx >= nhresp) exit
     fidx = lidx + 1
  ENDDO

  IF (hreswav(nhresp) > tmp_wav(ncalcp) + samprate * 0.5) THEN
     ncalcp = ncalcp + 1
     tmp_wav(ncalcp) = hreswav(nhresp)
  ELSE
     tmp_wav(ncalcp) = hreswav(nhresp)
  ENDIF

  ! Find the indices of radiance calc. wavelength in high resolution
  allocate (radcwav(ncalcp), radcidxs(ncalcp))
  radcwav = tmp_wav(1:ncalcp) 
  radcidxs(1:ncalcp) = 0
  radcidxs(1) = 1; j = 2
  DO i = 2, nhresp
     IF ( ABS(hreswav(i) - radcwav(j)) <= hres_samprate * 0.2 ) THEN
        radcidxs(j) = i; j = j + 1
     ENDIF
  ENDDO

  ! position bwt UV and VIS in hreswav
  ! Definition of widx_hvis moved to above, as it is used before it is defined.
  !widx_hvis = nhresp + 1
  IF ( nviswin > 0) THEN 
  !   widx_hvis = MINVAL(MAXLOC(hreswav(1:nhresp), MASK=(hreswav(1:nhresp) < wcenter_uvvis))) +1 
     widx_rvis = MINVAL(MAXLOC(radcwav(1:ncalcp), MASK=(radcwav(1:ncalcp) < wcenter_uvvis))) +1 
  ELSE
    !widx_hvis = 0  
    widx_rvis = 0
  ENDIF

  ! set PCA variables
  npca = 0
  DO i = 1, 1
    npca = npca + 1
    winpca(npca,1) = winlim(i,1) 
    winpca(npca,2) = winlim(nuvwin,2) 
    which_pcabin(npca) = 1
  ENDDO
  IF (nviswin > 0) THEN 
   DO i = 1,nviswin
    npca = npca + 1
    which_pcabin(npca) = 2
    !winpca(npca,1) = winlim(nuvwin+1,1) ; winpca(npca,2) = winlim(numwin,2)
    winpca(npca,1) = winlim(nuvwin+i,1) ; winpca(npca,2) = winlim(nuvwin+i,2)
   ENDDO
  ENDIF
  fidx = 1
  DO i = 1, npca
    lidx = MINVAL(MAXLOC(radcwav(1:ncalcp), MASK=(radcwav(1:ncalcp) <= winpca(i,2) &
            .and. radcwav(1:ncalcp) > 0 )))
    if (i == npca) lidx = ncalcp
    npcapix(i, 1) = fidx
    npcapix(i, 2) = lidx
    fidx = lidx + 1 
    WRITE(www_lun,'(A,3i5, 4f8.3)') 'PCAPIX:',i, npcapix(i,1), npcapix(i,2), &
    radcwav(npcapix(i,1)),radcwav(npcapix(i,2))
  ENDDO

  ! check boundaries
  IF (radcwav(1) < albmin(1) .OR. radcwav(ncalcp) > albmax(nalb)) THEN
     WRITE(*, *) modulename, ': Check albedo coverange for hreswav !!!'
     print * , albmin(1), albmax(nalb), radcwav(1), radcwav(ncalcp)
     errstat = pge_errstat_error; RETURN
  ENDIF
  WRITE (www_lun,*) 'N of hreswav/radcwav:', nhresp, ncalcp
  deallocate (tmp_i0, tmp_wav, hres0_wav, hres0_i0)
  RETURN
END SUBROUTINE get_hres_radcal_waves

SUBROUTINE hres_radwf_inter_convol(nw, nz, nctp, ncbp, nsprs, nalb, faerlvl,  &
     do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf, do_tracewf, &
     do_o3shi, do_tmpwf, do_pslwf, wave, ozs, do_abs, delabs, & 
     rad, fozwf, albwf, cfracwf, faerwf, &
     faerswf, fcodwf, fsprswf, fraywf, errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,   ONLY  :  &
      so2_idx, so2v_idx, o2o2_idx,o2_idx, o2t2_idx, h2o_idx, h2ot2_idx
  USE OMSAO_parameters_module,ONLY  : du2mol
  USE OMSAO_variables_module, ONLY  : numwin, winlim, &
       owave=>radwvl_sav, now=>n_radwvl_sav, i0sav, refidx,  & 
       nrad=>n_rad_wvl, &
       do_bandavg, refidx_sav, database_pslwf,npsl
  USE ozprof_data_module,     ONLY  : num_iter, hwave=>hreswav, &
       radcidxs, hres_i0, nhw=>nhresp, hresgabs, hresray, nw0=>ncalcp, &
       ngas, gasidxs, fgasidxs, fgassidxs, &
       o3crsz, o3dadtz, o3dadsz, so2crsz, o4crsz, o2crsz, h2ocrsz, &
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, & 
       ccrs, dads, dadt
  USE OMSAO_errstat_module
  USE m_get_xcrs, ONLY: calc_pslwf
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                              :: nw, nz, nctp, ncbp, faerlvl, nsprs, nalb
  INTEGER, INTENT(OUT)                             :: errstat
  LOGICAL, INTENT(IN)                              :: do_albwf, do_faerwf, do_faerswf, &
       do_codwf, do_sprswf, do_cfracwf, do_o3shi, do_tmpwf, do_tracewf, do_pslwf, do_abs

  REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
  REAL (KIND=dp), DIMENSION(nw),     INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fozwf, faerwf, faerswf, fcodwf, &
       fsprswf, fraywf, delabs
  REAL (KIND=dp), DIMENSION(nw, nalb),INTENT(INOUT) :: albwf
  REAL (KIND=dp), DIMENSION(nw),      INTENT(INOUT) :: rad, cfracwf

  ! Local variables
  !INTEGER, PARAMETER :: which_pslwf = 2
  INTEGER :: i, iwin, fidx, lidx, fidxc, lidxc, idx, iw, ntemp, sidx, eidx
  LOGICAL :: do_so2shi, do_o4shi, do_o2shi, do_h2oshi
  ! tracegases shift is not asscounted for T/P dependent cross section
  REAL (KIND=dp)                      :: temp
  INTEGER, DIMENSION (nw)             :: c2hfidx, c2hlidx

  REAL (KIND=dp), DIMENSION (nhw)     :: hrad, dtau, dray, hcfracwf
  REAL (KIND=dp), DIMENSION (nhw, nalb):: halbwf
  REAL (KIND=dp), DIMENSION (now)     :: oi0, otmp, tmpi0 
  REAL (KIND=dp), DIMENSION (nw,  nz) :: tauwf
  REAL (KIND=dp), DIMENSION (now, nz) :: dads1, dadt1, abscrs1, so2crs1, o4crs1,o2crs1, h2ocrs1
  REAL (KIND=dp), DIMENSION (nhw, nz) :: hozwf, haerwf, haerswf, hcodwf, hsprswf!, hraywf 
  REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: inarr
  REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: outarr

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'hres_radwf_convol'

  errstat = pge_errstat_ok
  allocate ( inarr(nhw, nz*8), outarr(now, nz*8))

  IF (num_iter == 0) THEN 
     IF (allocated (ccrs%o3)) deallocate(ccrs%o3)
     allocate (ccrs%o3(now, nz))
   IF (use_so2dtcrs) then
     IF (allocated (ccrs%so2)) deallocate(ccrs%so2)
     allocate (ccrs%so2(now, nz))
   ENDIF
   IF (use_o4dtcrs)  then
     IF (allocated (ccrs%o4)) deallocate(ccrs%o4)
     allocate (ccrs%o4(now, nz))
   ENDIF
   IF (use_h2odptcrs) then
     IF (allocated (ccrs%h2o)) deallocate(ccrs%h2o)
     allocate (ccrs%h2o(now, nz))
   ENDIF
   IF (use_o2dptcrs) then
     IF (allocated (ccrs%o2)) deallocate(ccrs%o2)
     allocate (ccrs%o2(now, nz))
   ENDIF
   IF (do_o3shi) THEN 
     IF (allocated(dads%o3)) deallocate(dads%o3)
     allocate (dads%o3(now, nz))
   ENDIF
   IF (do_tmpwf) THEN 
     IF (allocated(dadt%o3)) deallocate(dadt%o3)
     allocate (dadt%o3(now, nz))
   ENDIF
  ENDIF

  ! get weighting function in dlnI/dx and take the logarithm of radiances
  DO i = 1, nz
     fozwf(1:nw0, i) = fozwf(1:nw0, i) / rad(1:nw0)
  ENDDO

  IF (do_albwf) THEN
     DO i = 1, nalb
       albwf(1:nw0, i) = albwf(1:nw0, i) / rad(1:nw0)
     ENDDO
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
     DO i = 1, nalb
       CALL BSPLINE(wave(1:nw0), albwf(1:nw0,i), nw0, hwave(1:nhw), halbwf(1:nhw, i), nhw, errstat)
     ENDDO
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
    DO i = 1, nalb
      halbwf(:, i) = halbwf(:,i) * hrad
    ENDDO
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
    DO i = 1, nalb
     sidx = sidx + 1;     inarr(1:nhw, sidx) = hres_i0(1:nhw) * halbwf(1:nhw, i)
    ENDDO
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
     do_o4shi  = .FALSE. ; do_o2shi = .FALSE. ; do_h2oshi = .FALSE.
     DO i = 1, ngas
        IF (fgasidxs(i) > 0) THEN
           IF (gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN
              IF (fgassidxs(i) > 0 ) do_so2shi = .TRUE.
              !IF (use_so2dtcrs) CYCLE
           ENDIF
           IF (gasidxs(i) == o2o2_idx) THEN 
              IF (fgassidxs(i) > 0 ) do_o4shi = .TRUE.
           ENDIF     
           IF (gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx ) THEN
              IF (fgassidxs(i) > 0 ) do_o2shi = .TRUE.
           ENDIF
           IF (gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx ) THEN
              IF (fgassidxs(i) > 0 ) do_h2oshi = .TRUE.
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
           sidx = sidx + 1;   inarr(1:nhw, sidx) = o2crsz(1:nhw, i)  !* hres_i0(1:nhw) 
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

  IF (do_abs) THEN
     DO i = 1, nz
        sidx = sidx + 1;   inarr(1:nhw, sidx) = delabs(1:nhw, i) !* hi0(1:nhw)
     ENDDO
  ENDIF

  eidx = sidx
  
  IF (do_pslwf ) THEN 
     database_pslwf(refidx(1:now), 1:npsl) = calc_pslwf(hwave(1:nhw),hrad(1:nhw), nhw, &
              npsl, .false., 1.0D0, owave(1:now), now)
  ENDIF

  ! *** second, convole all spectra at once ****
  CALL convol_f2c(hwave(1:nhw), inarr(1:nhw, 1:eidx), nhw, eidx, owave(1:now), outarr(1:now, 1:eidx), now)
  ! *** third, transfer all convolved spectra back ***
  oi0(1:now)  = outarr(1:now, 1)
  hrad(1:now) = outarr(1:now, 2)  / oi0(1:now)

  DO i = 1, nz
     sidx = 2 + i;        hozwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
  ENDDO
  IF (do_albwf) THEN
    DO i = 1, nalb
     sidx = sidx + 1;     halbwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
    ENDDO
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

     IF (use_so2dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  so2crs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO

        !IF (do_so2shi) THEN
        !   sidx = sidx + 1;  so2dads1(1:now) = outarr(1:now, sidx) 
        !ENDIF
     ENDIF

     IF (use_o4dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  o4crs1(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO 
     ENDIF
     IF (use_o2dptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  o2crs1(1:now, i)  = outarr(1:now, sidx)  !/ oi0(1:now)
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
        print * , 'not implemented '; stop 1
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
       DO i = 1, nalb
        otmp(1:now) = halbwf(1:now, i) * tmpi0(1:now)
        CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        albwf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
      ENDDO
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
           dads%o3(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
        ENDDO
     ENDIF

     IF (do_tmpwf) THEN
        DO i = 1, nz
           otmp(1:now) = dadt1(1:now, i) !* tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           dadt%o3(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
        ENDDO
     ENDIF

     IF (do_tracewf) THEN
        DO i = 1, nz
           otmp(1:now) = abscrs1(1:now, i) * tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           ccrs%o3(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
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
              ccrs%so2(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
           ENDDO
        ENDIF
        IF (use_o4dtcrs) THEN 
           DO i = 1, nz
              otmp(1:now) = o4crs1(1:now, i) * tmpi0(1:now) 
              CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
              ccrs%o4(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
           ENDDO
        ENDIF
        IF (use_o2dptcrs) THEN 
           DO i = 1, nz
              otmp(1:now) = o2crs1(1:now, i) * tmpi0(1:now) 
              CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
              ccrs%o2(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
           ENDDO
        ENDIF
        IF (use_h2odptcrs) THEN 
           DO i = 1, nz
              otmp(1:now) = h2ocrs1(1:now, i) !* tmpi0(1:now) 
              CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
              ccrs%h2o(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
           ENDDO
        ENDIF
        
     ENDIF
     
  ELSE
     rad(1:now) = hrad(1:now)
     fozwf(1:now, 1:nz) = hozwf(1:now, 1:nz)
     !fraywf(1:now, 1:nz) = hraywf(1:now, 1:nz)
     IF (do_albwf)  albwf(1:now, 1:nalb) = halbwf(1:now, 1:nalb)
     IF (do_cfracwf)cfracwf(1:now) = hcfracwf(1:now)
     IF (do_faerwf) faerwf(1:now, faerlvl:nz) = haerwf(1:now, faerlvl:nz)
     IF (do_faerwf) faerswf(1:now, faerlvl:nz) = haerswf(1:now, faerlvl:nz)
     IF (do_codwf)  fcodwf(1:now, nctp:ncbp) = hcodwf(1:now, nctp:ncbp)
     IF (do_sprswf) fsprswf(1:now, nsprs:nz) = hcodwf(1:now, nsprs:nz)
     IF (do_o3shi)  dads%o3(1:now, 1:nz) = dads1(1:now, 1:nz)
     IF (do_tmpwf)  dadt%o3(1:now, 1:nz) = dadt1(1:now, 1:nz)
     IF (do_tracewf) THEN
           ccrs%o3(1:now, 1:nz)  = abscrs1(1:now, 1:nz)
        IF (use_so2dtcrs) THEN
           ccrs%so2(1:now, 1:nz) = so2crs1(1:now, 1:nz)
        ENDIF
        IF (use_o4dtcrs) THEN
           ccrs%o4(1:now, 1:nz) = o4crs1(1:now, 1:nz)
        ENDIF
        IF (use_o2dptcrs) THEN
           ccrs%o2(1:now, 1:nz) = o2crs1(1:now, 1:nz)
        ENDIF
        IF (use_h2odptcrs) THEN
           ccrs%h2o(1:now, 1:nz) = h2ocrs1(1:now, 1:nz)
        ENDIF
     ENDIF
  ENDIF

  IF (do_abs) THEN
     DO i = 1, nz
        sidx = sidx + 1;  delabs(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
     ENDDO
  ENDIF

  deallocate (inarr, outarr)
  RETURN

END SUBROUTINE hres_radwf_inter_convol

!1.	Establish fine wavelength grid: 0.01 nm now, may change to 0.05 nm later
!2.	Establish radiance calculation grid, based on spectral sampling rate for 
!   different specral regions
!3.	Find indices of radiance calculation grid in fine wavelength grid

SUBROUTINE hres_stkwf_inter_convol(nw, nz, nctp, ncbp, nsprs, faerlvl,  &
     do_faerwf, &
     wave, radq, fqozwf, fqaerwf, npol, pol_idx, errstat)

  USE OMSAO_precision_module   
  USE OMSAO_errstat_module
  USE m_convol, ONLY: convol_f2c_stk
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                              :: nw, nz, nctp, ncbp, faerlvl, nsprs
  INTEGER, INTENT(IN)                              :: npol
  INTEGER, DIMENSION(npol), INTENT(IN)             :: pol_idx                                                 
  LOGICAL, INTENT(IN)                              :: do_faerwf
  INTEGER, INTENT(OUT)                             :: errstat  

  !REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)   :: ozs
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fqozwf, fqaerwf 
       !faerswf, fcodwf, fsprswf, fraywf
  
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw), INTENT(INOUT) :: radq!, albwf, cfracwf

  ! Local variables
  INTEGER :: i
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: hwave, tmparr !(npol)  
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: owave, otmp !(n_gome_q)   
  
  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'hres_stkwf_convol'

  errstat = pge_errstat_ok

  allocate (hwave(npol), owave(n_gome_q), otmp(n_gome_q), tmparr(npol))

  owave = gome_q(1,1:n_gome_q) 

  hwave = wave(pol_idx)
  ! Convolve stokes fraction/weighting functions/solar reference at fine grids into PMD grid 
  ! (before coadding) solar reference
  tmparr = radq(pol_idx)
  !do i =1 , npol
  !   write(96,*) hwave(i), tmparr(i)
  !enddo
  !write(*,*) pol_idx
  !pause 
  CALL convol_f2c_stk(hwave(1:npol), tmparr, npol, owave(1:n_gome_q), otmp(1:n_gome_q), n_gome_q)
  radq(1:n_gome_q) = otmp(1:n_gome_q) 

  DO i =1, nz
    tmparr = fqozwf(pol_idx, i)
    CALL convol_f2c_stk(hwave(1:npol), tmparr, npol, owave(1:n_gome_q), otmp(1:n_gome_q), n_gome_q)
    fqozwf(1:n_gome_q, i) = otmp(1:n_gome_q) 
  ENDDO
  
  IF (do_faerwf) THEN
     DO i = faerlvl, nz
        tmparr = fqaerwf(pol_idx, i)
        CALL convol_f2c_stk(hwave(1:npol), tmparr, npol, owave(1:n_gome_q), otmp(1:n_gome_q), n_gome_q)
        fqaerwf(1:n_gome_q, i) = otmp(1:n_gome_q) 
     ENDDO
  ENDIF
 
  deallocate (hwave, owave, otmp, tmparr)
  RETURN

END SUBROUTINE hres_stkwf_inter_convol

SUBROUTINE hres_stkwf_inter_convol2(nw, wave, radq, errstat)

  USE OMSAO_precision_module   
  USE OMSAO_errstat_module
  USE m_convol, ONLY: convol_f2c_stk
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                              :: nw  
  INTEGER, INTENT(OUT)                             :: errstat                                                   
  
  REAL (KIND=dp), DIMENSION(nw),     INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: radq

  ! Local variables
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: hwave, tmparr !(nw)
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE  :: owave, otmp  !(n_gome_q)
  
  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=18), PARAMETER :: modulename = 'hres_stkwf_convol2'

  errstat = pge_errstat_ok  
  allocate (hwave(nw), owave(n_gome_q), otmp(n_gome_q), tmparr(nw))
  owave = gome_q(1,1:n_gome_q) 
  hwave = wave
  tmparr = radq

  CALL convol_f2c_stk(hwave(1:nw), tmparr, nw, owave(1:n_gome_q), otmp(1:n_gome_q), n_gome_q)
  radq(1:n_gome_q) = otmp(1:n_gome_q) 
  deallocate (hwave, owave, otmp, tmparr)
  RETURN

END SUBROUTINE hres_stkwf_inter_convol2

!-------------------------------------------------------------
!The same as hres_radwf_inter_convol but convol all stokes elements
!and wfs     (zcai, 2009-09-09 to be checked......)
SUBROUTINE hres_radwf_inter_convol_all(nostk, nw, nz, nctp, ncbp, nsprs, faerlvl,  &
     do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf, do_tracewf, &
     do_o3shi, do_tmpwf, wave, ozs, rad, radq,fqozwf, fozwf, albwf, cfracwf, faerwf, &
     faerswf, fcodwf, fsprswf, fraywf, errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,   ONLY  : so2_idx, so2v_idx
  USE OMSAO_variables_module, ONLY  : numwin,  winlim, &
       owave=>radwvl_sav, now=>n_radwvl_sav, i0sav, nrad=>n_rad_wvl, &
       do_bandavg, refidx_sav 
  USE ozprof_data_module,     ONLY  : num_iter, hwave=>hreswav, &
       radcidxs, hres_i0, nhw=>nhresp, hresgabs, hresray, nw0=>ncalcp, o3crsz, &
       o3dadsz, so2crsz, so2dads, ngas, &
       gasidxs, fgasidxs, fgassidxs, & 
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, & 
       ccrs, dads, dadt
  USE OMSAO_errstat_module
  USE m_avg_band, ONLY: avg_band_spec
  USE m_convol, ONLY: convol_f2c
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                              :: nostk, nw, nz, nctp, ncbp, faerlvl, nsprs
  INTEGER, INTENT(OUT)                             :: errstat                                                   
  LOGICAL, INTENT(IN)                              :: do_albwf, do_faerwf, do_faerswf, &
       do_codwf, do_sprswf, do_cfracwf, do_o3shi, do_tmpwf, do_tracewf

  REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
  REAL (KIND=dp), DIMENSION(nw, nz, nostk), INTENT(INOUT) :: fozwf, faerwf, faerswf, fcodwf, &
       fsprswf, fraywf
  REAL (KIND=dp), DIMENSION(nw),     INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw, nostk),     INTENT(INOUT) :: rad, albwf, cfracwf
  REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: radq    !-->
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fqozwf  !-->

  ! Local variables
  INTEGER :: i, j, iwin, fidx, lidx, fidxc, lidxc, idx, iw, ntemp
  LOGICAL :: do_so2shi
  REAL (KIND=dp)                      :: temp
  INTEGER, DIMENSION (nw)             :: c2hfidx, c2hlidx
  REAL (KIND=dp), DIMENSION (nw)      :: direc

  REAL (KIND=dp), DIMENSION (:,:),   ALLOCATABLE :: hrad, halbwf, hcfracwf !(nhw, nostk)
  REAL (KIND=dp), DIMENSION (:),     ALLOCATABLE :: dtau, dray, tmparr, hdirec,hradq !(nhw)
  REAL (KIND=dp), DIMENSION (:),     ALLOCATABLE :: oi0, otmp, tmpi0, so2dads1 !(now)
  REAL (KIND=dp), DIMENSION (:,:),   ALLOCATABLE :: dads1, dadt1, abscrs1, so2crs1 !(now, nz)
  REAL (KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: tauwf !(nw, nz, nostk)
  REAL (KIND=dp), DIMENSION (:,:),   ALLOCATABLE :: tmp_gas, tmp_gasshi !(ngas,now)
  REAL (KIND=dp), DIMENSION (:,:),   ALLOCATABLE :: hqozwf  !-->(nhw, nz)
  REAL (KIND=dp), DIMENSION (:,:,:), ALLOCATABLE :: hozwf, haerwf, haerswf, hcodwf, hsprswf, hraywf !(nhw, nz, nostk)
  !INTEGER :: ntime = 1

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=17), PARAMETER :: modulename = 'hres_radwf_convol'

  errstat = pge_errstat_ok

  ! local variables
    allocate (hrad(nhw, nostk), halbwf(nhw, nostk), hcfracwf(nhw, nostk))
    allocate (dtau(nhw), dray(nhw), tmparr(nhw), hdirec(nhw),hradq(nhw))
    allocate (oi0(now), otmp(now), tmpi0(now), so2dads1(now))    
    allocate (dads1(now, nz), dadt1(now, nz), abscrs1(now, nz), so2crs1(now, nz))
    allocate (tauwf(nw, nz, nostk))
    allocate (tmp_gas(ngas,now), tmp_gasshi(ngas,now))
    allocate (hozwf(nhw, nz, nostk), haerwf(nhw, nz, nostk), haerswf(nhw, nz, nostk), & 
             hcodwf(nhw, nz, nostk), hsprswf(nhw, nz, nostk), hraywf(nhw, nz, nostk))
    allocate (hqozwf(nhw, nz))
     
  ! output variables
  IF (num_iter == 0) THEN 
     IF (allocated (ccrs%o3)) deallocate(ccrs%o3)
     allocate (ccrs%o3(now, nz))
   IF (use_so2dtcrs) then
     IF (allocated (ccrs%so2)) deallocate(ccrs%so2)
     allocate (ccrs%so2(now, nz))
   ENDIF
   IF (use_o4dtcrs)  then
     IF (allocated (ccrs%o4)) deallocate(ccrs%o4)
     allocate (ccrs%o4(now, nz))
   ENDIF
   IF (use_h2odptcrs) then
     IF (allocated (ccrs%h2o)) deallocate(ccrs%h2o)
     allocate (ccrs%h2o(now, nz))
   ENDIF
   IF (use_o2dptcrs) then
     IF (allocated (ccrs%o2)) deallocate(ccrs%o2)
     allocate (ccrs%o2(now, nz))
   ENDIF
   IF (do_o3shi) THEN 
     IF (allocated(dads%o3)) deallocate(dads%o3)
     allocate (dads%o3(now, nz))
   ENDIF
   IF (do_tmpwf) THEN 
     IF (allocated(dadt%o3)) deallocate(dadt%o3)
     allocate (dadt%o3(now, nz))
   ENDIF
  ENDIF
  

  !write(*,*) 'nw, nw0, nhw', nw, nw0, nhw        !-->
  ! get weighting function in dlnI/dx and take the logarithm of radiances
  DO j = 1, nostk
    DO i = 1, nz
       fozwf(1:nw0, i, j) = fozwf(1:nw0, i, j) / rad(1:nw0, j)
    ENDDO
  
    IF (do_albwf) THEN
       albwf(1:nw0, j) = albwf(1:nw0, j) / rad(1:nw0, j)
    ENDIF

    IF (do_cfracwf) THEN
       cfracwf(1:nw0, j) = cfracwf(1:nw0, j) / rad(1:nw0, j)
    ENDIF

    IF (do_faerwf) THEN
       DO i = faerlvl, nz
          faerwf(1:nw0, i, j) = faerwf(1:nw0, i, j) / rad(1:nw0, j)
       ENDDO
    ENDIF

    IF (do_faerswf) THEN
       DO i = faerlvl, nz
          faerswf(1:nw0, i, j) = faerswf(1:nw0, i, j) / rad(1:nw0, j)
       ENDDO
    ENDIF

    IF (do_codwf) THEN
       DO i = nctp, ncbp
          fcodwf(1:nw0, i, j) = fcodwf(1:nw0, i, j) / rad(1:nw0, j)
       ENDDO
    ENDIF

    IF (do_sprswf) THEN
       DO i = nsprs, nz
          fsprswf(1:nw0, i, j) = fsprswf(1:nw0, i, j) / rad(1:nw0, j)
       ENDDO
    ENDIF

    DO i = 1, nz
       fraywf(1:nw0, i, j) = fraywf(1:nw0, i, j) / rad(1:nw0, j)
    ENDDO 
   ! never convert Q or V to ln(Q or V) directly
    IF (j>1) THEN
      direc=1.0     !remember the direction of Q,V
      WHERE(rad(1:nw0, j) < -1e-16)
         direc=-1.0
      END WHERE
      rad(1:nw0, j) = LOG(ABS(rad(1:nw0, j)))
    ELSE
      rad(1:nw0, j) = LOG(rad(1:nw0, j))
    ENDIF


  ! convert ozone weighting function to gas absorption weighting function
    DO i = 1, nz
       tauwf(1:nw0, i, j) = fozwf(1:nw0, i, j) * ozs(i) / hresgabs(radcidxs(1:nw0), i)
    ENDDO
  ENDDO ! end nostk loop


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
     DO j = 1, nostk ! Note: I, Q, U corrections are the same. ---> to be checked! 
        hrad(fidx:lidx, j) = rad(iw, j)   !Note: 
        IF (j>1) hdirec(fidx:lidx) = direc(iw) ! get high res direction--> to be checked!(*better to simulate high res radiance)
        DO i = 1, nz
           dtau(fidx:lidx) = hresgabs(fidx:lidx, i) - hresgabs(radcidxs(iw), i)
           dray(fidx:lidx) = hresray(fidx:lidx, i)  - hresray(radcidxs(iw), i)
           hrad(fidx:lidx, j) = hrad(fidx:lidx, j) + tauwf(iw, i, j) * dtau(fidx:lidx) + &
                                                   fraywf(iw, i, j) * dray(fidx:lidx)
        
        ! Is it better to do cubic interpolation ???
           hozwf(fidx:lidx, i, j) = fozwf(iw, i, j) * o3crsz(fidx:lidx, i) / o3crsz(radcidxs(iw), i)
           hqozwf(fidx:lidx, i) = fqozwf(iw, i) * o3crsz(fidx:lidx, i) / o3crsz(radcidxs(iw), i) !-->
        ENDDO
     ENDDO ! end nstok loop  
  ENDDO
  !DO i = 1, nz 
  !   CALL BSPLINE(wave(1:nw0), fraywf(1:nw0, i), nw0, hwave(1:nhw), hraywf(1:nhw, i), nhw, errstat)
  !   IF (errstat < 0) THEN
  !      WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
  !      errstat = pge_errstat_error; RETURN
  !   ENDIF
  !ENDDO
  !
  !DO i = 1, nz 
  !   CALL BSPLINE(wave(1:nw0), fozwf(1:nw0, i), nw0, hwave(1:nhw), hozwf(1:nhw, i), nhw, errstat)
  !   IF (errstat < 0) THEN
  !      WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
  !      errstat = pge_errstat_error; RETURN
  !   ENDIF
  !ENDDO
  DO j = 1, nostk
     IF (do_albwf) THEN
        CALL BSPLINE(wave(1:nw0), albwf(1:nw0, j), nw0, hwave(1:nhw), halbwf(1:nhw, j), nhw, errstat)
     ENDIF

     IF (do_cfracwf) THEN
        CALL BSPLINE(wave(1:nw0), cfracwf(1:nw0, j), nw0, hwave(1:nhw), hcfracwf(1:nhw, j), nhw, errstat)
     ENDIF

     IF (do_faerwf) THEN
        DO i = faerlvl, nz
           CALL BSPLINE(wave(1:nw0), faerwf(1:nw0, i, j), nw0, hwave(1:nhw), haerwf(1:nhw, i, j), nhw, errstat)
        ENDDO
     ENDIF
  
     IF (do_faerswf) THEN
        DO i = faerlvl, nz
           CALL BSPLINE(wave(1:nw0), faerswf(1:nw0, i, j), nw0, hwave(1:nhw), haerswf(1:nhw, i, j), nhw, errstat)
        ENDDO
     ENDIF

     IF (do_codwf) THEN
        DO i = nctp, ncbp
           CALL BSPLINE(wave(1:nw0), fcodwf(1:nw0, i, j), nw0, hwave(1:nhw), hcodwf(1:nhw, i, j), nhw, errstat)
        ENDDO
     ENDIF

     IF (do_sprswf) THEN
        DO i = nsprs, nz
           CALL BSPLINE(wave(1:nw0), fsprswf(1:nw0, i, j), nw0, hwave(1:nhw), hsprswf(1:nhw, i, j), nhw, errstat)
        ENDDO
     ENDIF
 
    ! Convert radiances back
    !hrad(:, j) = EXP(hrad(:, j))
    IF (j>1) THEN 
       hrad(:, j) = EXP(hrad(:, j))*hdirec(:) 
    ELSE 
       hrad(:, j) = EXP(hrad(:, j))
    ENDIF
    ! convert radiance/weighting function to dlnI/dx from dI/dx
    DO i = 1, nz
       hozwf(:, i, j) = hozwf(:, i, j) * hrad(:, j)
    ENDDO

    !DO i = 1, nz
    !   hraywf(:, i, j) = hraywf(:, i, j) * hrad(:, j)
    !ENDDO

    IF (do_albwf) THEN
       halbwf(:, j) = halbwf(:, j) * hrad(:, j)
    ENDIF

    IF (do_cfracwf) THEN
       hcfracwf(:, j) = hcfracwf(:, j) * hrad(:, j)
    ENDIF

    IF (do_faerwf) THEN
       DO i = faerlvl, nz
          haerwf(:, i, j) = haerwf(:, i, j) * hrad(:, j)
       ENDDO
    ENDIF

    IF (do_faerswf) THEN
       DO i = faerlvl, nz
          haerswf(:, i, j) = haerswf(:, i, j) * hrad(:, j)
       ENDDO
    ENDIF

    IF (do_codwf) THEN
       DO i = nctp, ncbp
          hcodwf(:, i, j) = hcodwf(:, i, j) * hrad(:, j)
       ENDDO
    ENDIF

    IF (do_sprswf) THEN
       DO i = nsprs, nz
          hsprswf(:, i, j) = hsprswf(:, i, j) * hrad(:, j)
       ENDDO
    ENDIF
  ENDDO ! end nstok loop
  !get q
  CALL BSPLINE(wave(1:nw0), radq(1:nw0), nw0, hwave(1:nhw), hradq(1:nhw), nhw, errstat)
  

  ! Convolve radiance/weighting functions/solar reference at fine grids into measurement grid 
  ! (before coadding) solar reference
  CALL convol_f2c(hwave(1:nhw), hres_i0(1:nhw), nhw, 1, owave(1:now), oi0(1:now), now)

  tmparr = hradq(1:nhw) !* hres_i0(1:nhw)							!-->
  CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)    !-->
  hradq(1:now) = otmp(1:now)  !/ oi0(1:now)         


DO j = 1, nostk

  tmparr = hres_i0(1:nhw) * hrad(1:nhw, j)
  CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
  hrad(1:now, j) = otmp(1:now) / oi0(1:now) 

  DO i = 1, nz
     tmparr = hres_i0(1:nhw) * hozwf(1:nhw, i, j)
     CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
     hozwf(1:now, i, j) = otmp(1:now) / oi0(1:now)

    IF (j==nostk) THEN
     tmparr =  hres_i0(1:nhw) * hqozwf(1:nhw, i)                                                !-->
     CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now) !-->
     hqozwf(1:now, i) = otmp(1:now)  / oi0(1:now)                                           !-->
    ENDIF
  ENDDO
  !DO i = 1, nz
  !   tmparr = hres_i0(1:nhw) * hraywf(1:nhw, i)
  !   CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
  !   hraywf(1:now, i) = otmp(1:now) / oi0(1:now)
  !ENDDO
  IF (do_albwf) THEN
     tmparr = hres_i0(1:nhw) * halbwf(1:nhw, j)
     CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
     halbwf(1:now, j) = otmp(1:now) / oi0(1:now)
  ENDIF
  IF (do_cfracwf) THEN
     tmparr = hres_i0(1:nhw) * hcfracwf(1:nhw, j)
     CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
     hcfracwf(1:now, j) = otmp(1:now) / oi0(1:now)
  ENDIF
  IF (do_faerwf) THEN
     DO i = faerlvl, nz
        tmparr = hres_i0(1:nhw) * haerwf(1:nhw, i, j)
        CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        haerwf(1:now, i, j) = otmp(1:now) / oi0(1:now) 
     ENDDO
  ENDIF
  IF (do_faerswf) THEN
     DO i = faerlvl, nz
        tmparr = hres_i0(1:nhw) * haerswf(1:nhw, i, j)
        CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        haerswf(1:now, i, j) = otmp(1:now) / oi0(1:now) 
     ENDDO
  ENDIF
  IF (do_codwf) THEN
     DO i = nctp, ncbp
        tmparr = hres_i0(1:nhw) * hcodwf(1:nhw, i, j)
        CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        hcodwf(1:now, i, j) = otmp(1:now) / oi0(1:now) 
     ENDDO
  ENDIF
  IF (do_sprswf) THEN
     DO i = nsprs, nz
        tmparr = hres_i0(1:nhw) * hsprswf(1:nhw, i, j)
        CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        hsprswf(1:now, i, j) = otmp(1:now) / oi0(1:now) 
     ENDDO
  ENDIF
ENDDO ! end nstok loop
  ! convolve ozone shift/temperature 
  IF (do_tracewf) THEN
     DO i = 1, nz
        tmparr = o3crsz(1:nhw, i) * hres_i0(1:nhw) 
        CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        abscrs1(1:now, i) = otmp(1:now) / oi0(1:now)
     ENDDO

     do_so2shi = .FALSE.
     DO i = 1, ngas
        IF (fgasidxs(i) > 0) THEN
           IF (gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN
              IF (fgassidxs(i) > 0 ) do_so2shi = .TRUE.
              IF (use_so2dtcrs) CYCLE
           ENDIF
                
           ! This is not necessary: could still use those effective cross sections
           !tmparr = hres_gas(i, 1:nhw)* hres_i0(1:nhw) 
           !CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
           !tmp_gas(i, 1:now) = otmp(1:now) / oi0(1:now)
           !
           !IF (fgassidxs(i) > 0) THEN
           !   tmparr = hres_gasshi(i, 1:nhw)
           !   CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
           !   tmp_gasshi(i, 1:now) = otmp(1:now)
           !ENDIF
        ENDIF
     ENDDO

     IF (use_so2dtcrs) THEN
        DO i = 1, nz
           tmparr = so2crsz(1:nhw, i) * hres_i0(1:nhw) 
           CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
           so2crs1(1:now, i) = otmp(1:now) / oi0(1:now)
        ENDDO

        !IF (do_so2shi) THEN
        !   tmparr = so2dads(1:nhw) 
        !   CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        !   so2dads1(1:now) = otmp(1:now) 
        !ENDIF
     ENDIF
  ENDIF

  IF (do_o3shi) THEN
     DO i = 1, nz
        tmparr = o3dadsz(1:nhw, i) * hres_i0(1:nhw) 
        CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        dads1(1:now, i) = otmp(1:now) / oi0(1:now)
     ENDDO
  ENDIF

  IF (do_tmpwf) THEN
     DO i = 1, nz
        tmparr = hres_i0(1:nhw) !* o3dadtz(1:nhw, i)
        CALL convol_f2c(hwave(1:nhw), tmparr, nhw, 1, owave(1:now), otmp(1:now), now)
        dadt1(1:now, i) = otmp(1:now) ! / oi0(1:now)
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
    DO j = 1, nostk
     otmp(1:now) = hrad(1:now, j) * tmpi0(1:now)
     CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
     rad(1:nrad, j) = otmp(1:nrad) / oi0(1:nrad)

     DO i = 1, nz
        otmp(1:now) = hozwf(1:now, i, j) * tmpi0(1:now)
        CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        fozwf(1:nrad, i, j) = otmp(1:nrad) / oi0(1:nrad)
     ENDDO

     !DO i = 1, nz
     !   otmp(1:now) = hraywf(1:now, i) * tmpi0(1:now)
     !   CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
     !   fraywf(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
     !ENDDO

     IF (do_albwf) THEN
        otmp(1:now) = halbwf(1:now, j) * tmpi0(1:now)
        CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        albwf(1:nrad, j) = otmp(1:nrad) / oi0(1:nrad)
     ENDIF

     IF (do_cfracwf) THEN
        otmp(1:now) = hcfracwf(1:now, j) * tmpi0(1:now)
        CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
        cfracwf(1:nrad, j) = otmp(1:nrad) / oi0(1:nrad)
     ENDIF

     IF (do_faerwf) THEN
        DO i = faerlvl, nz
           otmp(1:now) = haerwf(1:now, i, j) * tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           faerwf(1:nrad, i, j) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO 
     ENDIF

     IF (do_faerswf) THEN
        DO i = faerlvl, nz
           otmp(1:now) = haerswf(1:now, i, j) * tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           faerswf(1:nrad, i, j) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO 
     ENDIF
     
     IF (do_codwf) THEN
        DO i = nctp, ncbp
           otmp(1:now) = hcodwf(1:now, i, j) * tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           fcodwf(1:nrad, i, j) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO 
     ENDIF

     IF (do_sprswf) THEN
        DO i = nsprs, nz
           otmp(1:now) = hsprswf(1:now, i, j) * tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           fsprswf(1:nrad, i, j) = otmp(1:nrad) / oi0(1:nrad)
        ENDDO 
     ENDIF
   ENDDO ! end stokes loop
     IF (do_o3shi) THEN
        DO i = 1, nz
           otmp(1:now) = dads1(1:now, i)  !* tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           dads%o3(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
        ENDDO
     ENDIF

     IF (do_tmpwf) THEN
        DO i = 1, nz
           otmp(1:now) = dadt1(1:now, i) !* tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           dadt%o3(1:nrad, i) = otmp(1:nrad) !/ oi0(1:nrad)
        ENDDO
     ENDIF

     IF (do_tracewf) THEN
        DO i = 1, nz
           otmp(1:now) = abscrs1(1:now, i) * tmpi0(1:now)
           CALL avg_band_spec(owave(1:now), otmp(1:now), now, ntemp, errstat)
           ccrs%o3(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
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
              ccrs%so2(1:nrad, i) = otmp(1:nrad) / oi0(1:nrad)
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
     radq(1:now)= hradq(1:now)               !-->
     fqozwf(1:now, 1:nz)= hqozwf(1:now, 1:nz)!-->
     rad(1:now, 1:nostk) = hrad(1:now, 1:nostk)    
     fozwf(1:now, 1:nz, 1:nostk) = hozwf(1:now, 1:nz, 1:nostk)
     !fraywf(1:now, 1:nz) = hraywf(1:now, 1:nz)
     IF (do_albwf) albwf(1:now, 1:nostk) = halbwf(1:now, 1:nostk)
     IF (do_cfracwf) cfracwf(1:now, 1:nostk) = hcfracwf(1:now, 1:nostk)
     IF (do_faerwf) faerwf(1:now, faerlvl:nz, 1:nostk) = haerwf(1:now, faerlvl:nz, 1:nostk)
     IF (do_faerwf) faerswf(1:now, faerlvl:nz, 1:nostk) = haerswf(1:now, faerlvl:nz, 1:nostk)
     IF (do_codwf) fcodwf(1:now, nctp:ncbp, 1:nostk) = hcodwf(1:now, nctp:ncbp, 1:nostk)
     IF (do_sprswf) fsprswf(1:now, nsprs:nz, 1:nostk) = hcodwf(1:now, nsprs:nz, 1:nostk)
     IF (do_o3shi) dads%o3(1:now, 1:nz) = dads1(1:now, 1:nz)
     IF (do_tmpwf) dadt%o3(1:now, 1:nz) = dadt1(1:now, 1:nz)
     IF (do_tracewf) THEN
        ccrs%o3(1:now, 1:nz) = abscrs1(1:now, 1:nz)
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
           ccrs%so2(1:now, 1:nz) = so2crs1(1:now, 1:nz)
           !IF (do_so2shi) THEN
           !   database_shiwf(so2_idx, refidx(1:now)) = so2dads1(1:now)  
           !   database_shiwf(so2v_idx, refidx(1:now)) = so2dads1(1:now) 
           !ENDIF
        ENDIF
     ENDIF
  ENDIF

  deallocate (hrad, halbwf, hcfracwf)
  deallocate (dtau, dray, tmparr, hdirec,hradq)
  deallocate (oi0, otmp, tmpi0, so2dads)
  deallocate (dads1, dadt1, abscrs1, so2crs1, tauwf)
  deallocate (tmp_gas, tmp_gasshi)
  deallocate (hozwf, haerwf, haerswf, hcodwf, hsprswf, hraywf)
  deallocate (hqozwf)
  RETURN

END SUBROUTINE hres_radwf_inter_convol_all

SUBROUTINE radwf_inter_convol(nw, nz, nctp, ncbp, nsprs, nalbwf, faerlvl,  &
     do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf, do_tracewf, &
     do_o3shi, do_tmpwf, do_pslwf, wave, ozs, do_abs, delabs, rad, & 
     fozwf, albwf, cfracwf, faerwf, faerswf, fcodwf, fsprswf, fraywf, & 
     errstat)

  USE OMSAO_precision_module
  USE OMSAO_indices_module,   ONLY  :  &
      so2_idx, so2v_idx, o2o2_idx,o2_idx, o2t2_idx, h2o_idx, h2ot2_idx
  USE OMSAO_variables_module, ONLY  : nviswin, &
       owave=>radwvl_sav, now=>n_radwvl_sav, refidx, & 
       database_pslwf,npsl, widx_hvis
  USE ozprof_data_module,     ONLY  : num_iter,& 
       hwave=>hreswav,hres_i0, nhw=>nhresp, &
       ngas, gasidxs, fgasidxs, fgassidxs, &
       o3crsz, o3dadtz, o3dadsz, so2crsz,o4crsz,o2crsz, h2ocrsz, &
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
       ccrs, dadt, dads
  USE OMSAO_errstat_module
  USE m_get_xcrs, ONLY: calc_pslwf
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN)                              :: nw, nz, nctp, ncbp, faerlvl, nsprs, nalbwf
  INTEGER, INTENT(OUT)                             :: errstat                                                   
  LOGICAL, INTENT(IN)                              :: do_albwf, do_faerwf, do_faerswf, &
   do_codwf, do_sprswf, do_cfracwf, do_o3shi, do_tmpwf, do_tracewf,do_pslwf, do_abs

  REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(INOUT) :: fozwf, faerwf, faerswf, fcodwf, &
       fsprswf, fraywf, delabs
  REAL (KIND=dp), DIMENSION(nw),     INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw),     INTENT(INOUT) :: rad, cfracwf
  REAL (KIND=dp), DIMENSION(nw, nalbwf),  INTENT(INOUT) :: albwf

  ! Local variables
  LOGICAL, PARAMETER :: do_rmvis_i0corr = .false.
  LOGICAL :: do_so2shi, do_o4shi, do_o2shi, do_h2oshi
  INTEGER :: i, sidx, eidx, nvar
  REAL (KIND=dp), DIMENSION (:), ALLOCATABLE  :: hi0, oi0
  REAL (KIND=dp), DIMENSION (:,:),ALLOCATABLE :: inarr
  REAL (KIND=dp), DIMENSION (:,:),ALLOCATABLE :: outarr
  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=18), PARAMETER :: modulename = 'radwf_inter_convol'
  errstat = pge_errstat_ok
  nvar = nz*10
  allocate (inarr(nhw, nvar), oi0(now), outarr(now, nvar))
  IF (num_iter == 0) THEN
     IF (allocated (ccrs%o3)) deallocate(ccrs%o3)
     allocate (ccrs%o3(now, nz))
   IF (use_so2dtcrs) then
     IF (allocated (ccrs%so2)) deallocate(ccrs%so2)
     allocate (ccrs%so2(now, nz))
   ENDIF
   IF (use_o4dtcrs)  then
     IF (allocated (ccrs%o4)) deallocate(ccrs%o4)
     allocate (ccrs%o4(now, nz))
   ENDIF
   IF (use_h2odptcrs) then
     IF (allocated (ccrs%h2o)) deallocate(ccrs%h2o)
     allocate (ccrs%h2o(now, nz))
   ENDIF
   IF (use_o2dptcrs) then
     IF (allocated (ccrs%o2)) deallocate(ccrs%o2)
     allocate (ccrs%o2(now, nz))
   ENDIF
   IF (do_o3shi) THEN
     IF (allocated (dads%o3)) deallocate(dads%o3)
     allocate (dads%o3(now, nz))
   ENDIF
   IF (do_tmpwf) THEN 
     IF (allocated (dadt%o3)) deallocate(dadt%o3)
     allocate (dadt%o3(now, nz))
   ENDIF 
 ENDIF 
   ! Convolve radiance/weighting functions/solar reference at fine grids into measurement grid 
  ! convolve all spectra at once to speed up computation
  ! *** First, transfer all spectra to inarr ***
  allocate (hi0(nhw))
  hi0 = hres_i0
  IF (do_rmvis_i0corr .and. nviswin > 0) THEN 
    hi0(widx_hvis:nhw) = 1.0
  ENDIF
  inarr(1:nhw, 1) = hi0(1:nhw)
  inarr(1:nhw, 2) = hi0(1:nhw) * rad(1:nhw)

  IF (do_pslwf) THEN 
     database_pslwf(refidx(1:now), 1:npsl) = calc_pslwf(hwave(1:nhw),rad(1:nhw), nhw, &
              npsl, .false., 1.0D0, owave(1:now), now)
  ENDIF 

  DO i = 1, nz
     sidx = 2 + i;        inarr(1:nhw, sidx) = hi0(1:nhw) * fozwf(1:nhw, i)
  ENDDO
  IF (do_albwf) THEN
    DO i = 1, nalbwf 
     sidx = sidx + 1;     inarr(1:nhw, sidx) = hi0(1:nhw) * albwf(1:nhw, i)
    ENDDO
  ENDIF
  IF (do_cfracwf) THEN
     sidx = sidx + 1;     inarr(1:nhw, sidx) = hi0(1:nhw) * cfracwf(1:nhw)
  ENDIF
  IF (do_faerwf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hi0(1:nhw) * faerwf(1:nhw, i)
     ENDDO
  ENDIF
  IF (do_faerswf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hi0(1:nhw) * faerswf(1:nhw, i)
     ENDDO
  ENDIF
  IF (do_codwf) THEN
     DO i = nctp, ncbp
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hi0(1:nhw) * fcodwf(1:nhw, i)
     ENDDO
  ENDIF
  IF (do_sprswf) THEN
     DO i = nsprs, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = hi0(1:nhw) * fsprswf(1:nhw, i)
     ENDDO
  ENDIF
  ! convolve ozone shift/temperature 
  IF (do_tracewf) THEN
     DO i = 1, nz
        sidx = sidx + 1;  inarr(1:nhw, sidx) = o3crsz(1:nhw, i) * hi0(1:nhw) 
     ENDDO

     do_so2shi = .FALSE.
     do_o4shi  = .FALSE. ; do_o2shi = .FALSE. ; do_h2oshi = .FALSE.
     DO i = 1, ngas
        IF (fgasidxs(i) > 0) THEN
           IF (gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN
              IF (fgassidxs(i) > 0 ) do_so2shi = .TRUE.
              !IF (use_so2dtcrs) CYCLE
           ENDIF
           IF (gasidxs(i) == o2o2_idx) THEN 
              IF (fgassidxs(i) > 0 ) do_o4shi = .TRUE.
           ENDIF     
           IF (gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx ) THEN
              IF (fgassidxs(i) > 0 ) do_o2shi = .TRUE.
           ENDIF
           IF (gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx ) THEN
              IF (fgassidxs(i) > 0 ) do_h2oshi = .TRUE.
           ENDIF
           ! This is not necessary: could still use those effective cross sections
           ! sidx = sidx + 1;  inarr(1:nhw, sidx) = hres_gas(i, 1:nhw) * hi0(1:nhw) 
           
           !IF (fgassidxs(i) > 0) THEN
           ! sidx = sidx + 1;  inarr(1:nw, sidx) = hres_gasshi(i, 1:nhw)
           !ENDIF
        ENDIF
     ENDDO

     IF (use_so2dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = so2crsz(1:nhw, i) * hi0(1:nhw) 
        ENDDO
        !IF (do_so2shi) THEN
        !   sidx = sidx + 1;  inarr(1:nhw, sidx) = so2dads(1:nhw) 
        !ENDIF
     ENDIF

     IF (use_o4dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = o4crsz(1:nhw, i) * hi0(1:nhw) 
        ENDDO
        !DO i = 1, nhw 
        !   print * ,hwave(i),  o4crsz(i, 20)
        !ENDDO
     ENDIF

     IF (use_o2dptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = o2crsz(1:nhw, i)  !* hi0(1:nhw) 
        ENDDO
       !DO i = 1, nhw 
       !   print * ,hwave(i),  o4crsz(i, 20)
       !ENDDO
     ENDIF

     IF (use_h2odptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;   inarr(1:nhw, sidx) = h2ocrsz(1:nhw, i) !* hi0(1:nhw) 
        ENDDO
     ENDIF
  ENDIF

  IF (do_o3shi) THEN
     DO i = 1, nz
        sidx = sidx + 1;   inarr(1:nhw, sidx) = o3dadsz(1:nhw, i) !* hi0(1:nhw)
     ENDDO
  ENDIF
  IF (do_tmpwf) THEN
     DO i = 1, nz
        sidx = sidx + 1;   inarr(1:nhw, sidx) = o3dadtz(1:nhw, i) !* hi0(1:nhw)
     ENDDO
  ENDIF
  IF (do_abs) THEN
     DO i = 1, nz
        sidx = sidx + 1;   inarr(1:nhw, sidx) = delabs(1:nhw, i) !* hi0(1:nhw)
     ENDDO
  ENDIF
  eidx = sidx
  IF (nvar < eidx) THEN 
     WRITE(*,*) ADJUSTL(TRIM(modulename))//'nvar < eidx'
     stop 1
  ENDIF
  ! *** second, convole all spectra at once ****
  CALL convol_f2c(hwave(1:nhw), inarr(1:nhw, 1:eidx), nhw, eidx, owave(1:now), outarr(1:now, 1:eidx), now)
   
  ! *** third, transfer all convolved spectra back ***
  oi0(1:now)  = outarr(1:now, 1)
  rad(1:now)  = outarr(1:now, 2)  / oi0(1:now)

  DO i = 1, nz
     sidx = 2 + i;        fozwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
  ENDDO
  IF (do_albwf) THEN
    DO i = 1, nalbwf
     sidx = sidx + 1;     albwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
    ENDDO
  ENDIF
  IF (do_cfracwf) THEN
     sidx = sidx + 1;     cfracwf(1:now) = outarr(1:now, sidx) / oi0(1:now)
  ENDIF
  IF (do_faerwf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  faerwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF
  IF (do_faerswf) THEN
     DO i = faerlvl, nz
        sidx = sidx + 1;  faerswf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF
  IF (do_codwf) THEN
     DO i = nctp, ncbp
        sidx = sidx + 1;  fcodwf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF
  IF (do_sprswf) THEN
     DO i = nsprs, nz
        sidx = sidx + 1;  fsprswf(1:now, i) = outarr(1:now, sidx) / oi0(1:now)
     ENDDO
  ENDIF

  ! convolve ozone shift/temperature 
  IF (do_tracewf) THEN
     DO i = 1, nz
        sidx = sidx + 1;  ccrs%o3(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)

     ENDDO
     
     IF (use_so2dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  ccrs%so2(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO

        !IF (do_so2shi) THEN
        !   sidx = sidx + 1;  so2dads1(1:now) = outarr(1:now, sidx) 
        !ENDIF
     ENDIF

     IF (use_o4dtcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  ccrs%o4(1:now, i)  = outarr(1:now, sidx) / oi0(1:now)
        ENDDO 
     ENDIF
     IF (use_o2dptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  ccrs%o2(1:now, i)  = outarr(1:now, sidx)  !/ oi0(1:now)
        ENDDO
     ENDIF
     IF (use_h2odptcrs) THEN
        DO i = 1, nz
           sidx = sidx + 1;  ccrs%h2o(1:now, i)  = outarr(1:now, sidx)  !/ oi0(1:now)
        ENDDO
     ENDIF
  ENDIF

  IF (do_o3shi) THEN
     DO i = 1, nz
        sidx = sidx + 1;  dads%o3(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
     ENDDO
  ENDIF
  IF (do_tmpwf) THEN
     DO i = 1, nz
        sidx = sidx + 1;  dadt%o3(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
     ENDDO
  ENDIF
  IF (do_abs) THEN
     DO i = 1, nz
        sidx = sidx + 1;  delabs(1:now, i)  = outarr(1:now, sidx) !/ oi0(1:now)
     ENDDO
  ENDIF
  deallocate (inarr, oi0, outarr, hi0)
  RETURN

END SUBROUTINE radwf_inter_convol

SUBROUTINE radwf_interpol(nw, nz, nctp, ncbp, nsprs, faerlvl, do_radcals, &
     do_fozwf, do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_cfracwf, wave, abscrs, &
     ozs, rad, fozwf, albwf, cfracwf, faerwf, faerswf, fcodwf, fsprswf, errstat)

  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY  : numwin, nradpix, band_selectors 
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
  REAL (KIND=dp)                            :: toz, df
  INTEGER, DIMENSION(:), ALLOCATABLE        :: didxs, uidxs!, radcals
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: effcrs, a, b, &
       crs1, crs2, tmp1, tmp2, wav1, wav2

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=14), PARAMETER :: modulename = 'radwf_interpol'


  errstat = pge_errstat_ok
  allocate ( didxs(nw), uidxs(nw))
  allocate (effcrs(nw), a(nw), b(nw), & 
            crs1(nw), crs2(nw), tmp1(nw), tmp2(nw), wav1(nw), wav2(nw))

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
  !!stop 1

  deallocate (didxs, uidxs, effcrs, a, b, crs1, crs2, tmp1, tmp2, wav1, wav2)
  RETURN

END SUBROUTINE radwf_interpol

! Determine wavelengths where exact polarization correction (NSTOKES: 4 vs 1) is
! calculated
! In UV1 (or between 270 and 310 nm): ~292 nm, ~298 nm, ~300 nm, ~302 nm, ~304
! nm, ~306 nm, last wavelength
! In UV2 (or between 310 and 340 nm): first, 1/4, middle and last wavelength
! So exact vector LIDORT calculation is done at 11 wavelengths.
! This option works when radiance interpolation option is turned on
SUBROUTINE set_polcorr (numwin, winlim, nw, waves, do_radcals, npolcorr, do_polcorrs, polidx, polcorr_idxs)
  USE OMSAO_precision_module, ONLY:dp
  USE ozprof_data_module, ONLY:mpolcorr
  IMPLICIT NONE
  !=============================================
  !INPUT VARIABLES
  INTEGER, INTENT(IN) :: nw, numwin
  LOGICAL, INTENT(IN) :: do_radcals(nw)
  REAL (KIND=dp), INTENT(IN), DIMENSION(numwin, 2) :: winlim
  REAL (KIND=dp), INTENT(IN) :: waves (nw)
  !=============================================
  !OUTOUT VARIABLES
  INTEGER, INTENT(OUT) :: npolcorr
  INTEGER, INTENT(OUT), DIMENSION (mpolcorr) :: polcorr_idxs
  INTEGER, INTENT(OUT), DIMENSION (nw)       :: polidx
  LOGICAL, INTENT(OUT), DIMENSION (nw)       :: do_polcorrs
  !LOCAL VARIABLES
  INTEGER :: fidx, idum, iw, lidx, i, j, k, jj, kk, jk
  REAL (KIND=dp) :: temp
  
    do_polcorrs(1:nw) = .FALSE. ; npolcorr = 0
    fidx = 1; idum = 0
     DO iw = 1, numwin
       IF (iw == numwin) THEN
           temp = winlim(iw, 2)
       ELSE
           temp = (winlim(iw, 2) + winlim(iw + 1, 1)) / 2.
       ENDIF
       lidx = MINVAL(MAXLOC(waves(1:nw), MASK=(waves(1:nw) < temp .AND. waves(1:nw) > 0)))
          IF ( waves(lidx) <= 312.0 ) THEN
          ! Error from using single scattering is about 0.2% at 270 nm, it
          ! needs o be corrected
          !IF (do_ssfullb295) do_polcorrs(fidx) = .TRUE.
            IF (idum == 0) idum = 1
             DO i = fidx + 1, lidx - 1
                IF (do_radcals(i) ) THEN
                  IF ( (waves(idum) < 290. .AND. waves(i) >= 290.0) .OR.  &
                       (waves(idum) < 295. .AND. waves(i) >= 295.0) .OR.  &
                       (waves(idum) < 299. .AND. waves(i) >= 299.0) .OR.  &
                       (waves(idum) < 301. .AND. waves(i) >= 301.0) .OR.  &
                       (waves(idum) < 303. .AND. waves(i) >= 303.) .OR.  &
                       !(waves(idum) < 304. .AND. waves(i) >= 304.0).OR.  &
                       (waves(idum) < 305. .AND. waves(i) >= 305.0) ) &
                       do_polcorrs(i) = .TRUE.
         !          print * , do_polcorrs(i), i, waves(idum)             
                       idum = i
                ENDIF
             ENDDO
             do_polcorrs(lidx) = .TRUE.
          ELSE
            IF (idum == 0) idum = 1
              do_polcorrs(fidx) = .TRUE.; do_polcorrs(lidx) = .TRUE.
              j = (fidx + lidx ) / 2; k = fidx + (j - fidx) / 3

              jj = (j + lidx )  / 2; kk = fidx + (k - fidx) / 3
              jk = (j + k ) / 2

              DO i = fidx + 1, lidx - 1
                 IF (do_radcals(i) ) THEN
                    IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                         (waves(idum) < waves(k) .AND. waves(i) >= waves(k)) .OR. &
                         (waves(idum) < waves(jj) .AND. waves(i) >= waves(jj)) .OR. &
                         (waves(idum) < waves(kk) .AND. waves(i) >= waves(kk)) .OR. &
                         (waves(idum) < waves(jk) .AND. waves(i) >= waves(jk))) do_polcorrs(i) = .TRUE.
                    !IF ( (waves(idum) < waves(j) .AND. waves(i) >= waves(j)) .OR. &
                    !     (waves(idum) < waves(k) .AND. waves(i) >= waves(k))) & 
                    !      do_polcorrs(i) = .TRUE.
                    idum = i
                 ENDIF
              ENDDO
          ENDIF
          fidx = lidx + 1
     ENDDO

     DO i = 1, nw
       IF ( do_polcorrs(i) ) THEN
          npolcorr = npolcorr + 1
          polcorr_idxs(npolcorr) = i
          polidx(i) = npolcorr
       ENDIF
     ENDDO
 RETURN
END SUBROUTINE

SUBROUTINE polcorr_online_with_lut(niter, VLDLUTdir, nw,nz, nctp,nsprs,nalb, &
           do_albwf, do_cfracwf, the_cfrac, albs, &
           sza, vza, phi, lat, ozdu,&
           ps, pmid, &
           wave,  rad, tauwf, tabs, tray,&
           albwf, fozwf,cfracwf)
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : mflay, max_spec_pts
  USE OMSAO_variables_module,  ONLY : currloop
  USE ozprof_data_module,      ONLY : & 
      vza_min, vza_max, sza_min, sza_max,do_polut_init
  USE OMSAO_errstat_module
  USE m_ezspline_interpolation, ONLY: interpol
  USE LaGrangePolynomialCoefficient_m !, ONLY : MAX_LG_INTERP, get_index
  USE LUT_Stage_m
  USE LUT_Itoa_m
  USE LUT_util_m, ONLY: getLwt
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  CHARACTER (LEN=*), INTENT(IN) :: VLDLUTdir !='/home/jbak/data/GEMSTOOL/lutdatav2.8-r/LUT-48/'
  LOGICAL, INTENT(IN) :: do_cfracwf, do_albwf
  INTEGER, INTENT(IN) :: nw, nz,  niter, nctp,  nsprs, nalb
  REAL (KIND=dp),                       INTENT(IN)    :: the_cfrac
  REAL (KIND=dp),                       INTENT(IN)    :: & 
      sza, vza,phi, lat, ozdu
  REAL (KIND=dp),                       INTENT(IN)    :: ps(0:nz),pmid(1:nz) 
  REAL (KIND=dp), DIMENSION(nw,2),      INTENT(IN)    :: albs
  REAL (KIND=dp), DIMENSION(nw, nz),    INTENT(IN)    :: tabs, tray
  REAL (KIND=dp), DIMENSION(nw, nz),    INTENT(INOUT) :: fozwf, tauwf
  REAL (KIND=dp), DIMENSION(nw),        INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw),        INTENT(INOUT) :: rad, cfracwf
  REAL (KIND=dp), DIMENSION(nw, nalb),  INTENT(INOUT)   :: albwf
  ! Local variables
  INTEGER    :: i, j, ic,fidx, lidx, errstat
  REAL (KIND=dp):: frac,frac1, frac2, pt,ctp,sfcp, adj_o3, adj_ray
  REAL (KIND=dp), DIMENSION(nz) :: dtau
  ! LUT variables
  LOGICAL, PARAMETER :: do_plan = .true., do_debug=.false.
  CHARACTER (LEN=255) :: msg
  !CHARACTER (LEN=255), parameter :: VLDLUTdir1='/home/jbak/data/GEMSTOOL/lutdatav2.8-r/LUT-48/'
  !CHARACTER (LEN=255), parameter :: VLDLUTdir1='/home/jbak/data/GEMSTOOL/lutdatav2.8-o3-UV/LUT-24/'
  CHARACTER (LEN=255), parameter :: VLDLUTdir1='/home/jbak/OzoneFit/tbl/vldlut/'
  CHARACTER(LEN=12), DIMENSION(2)  :: LUT_type = (/"vec06st72nl_", "sca02st26nl_"  /)
  LOGICAL :: log_ret = .FALSE., L_stageJ, LFAIL
  REAL (KIND=4), DIMENSION(1:13) :: vza_grid = & 
       (/0.0, 15.0, 30.0, 43.0, 53.0, 61.0, 67.0, 72.0, 76.0, 80.0, 84.0, 86.0, 88.0/)
  REAL (KIND=4), DIMENSION(1:12) :: sza_grid = & 
       (/0.0, 16.0, 31.0, 44.0, 55.0, 64.0, 71.0, 76.5, 80.5, 83.5, 86.0, 88.0/)
  INTEGER :: ib,ilat,ipt, iv_min, iv_max
  INTEGER, SAVE :: npJ, nwLUT
  REAL (KIND=4), DIMENSION(1:2) :: szam, vzam, wlm
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:), SAVE :: pJ,pL,wl
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:), SAVE :: & 
    Itoa_ic, i0, tr, sb, wta_ic, wtc_ic
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:),  SAVE :: wto_ic
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:),  SAVE :: i0V,i0S,trV,trS,sbS,sbV
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:),SAVE :: di0dtV,di0dtS,dtrdtV,dtrdtS,dsbdtS,dsbdtV
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:),SAVE :: di0dtV2,di0dtS2,dtrdtV2,dtrdtS2,dsbdtS2,dsbdtV2
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)     ::  alb, alb1, alb2, qr1, qr2
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:)   ::  Jvert0, Jvert1,  Jvert2, Jvert3,Jvert4
  ! I and wts @ LUT grids used to derive correction spectrum
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:),  SAVE :: wta, wtc,itoa,taucum_lut, taucum2_lut
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:),SAVE :: wto, wto2
  ! correction spectrum @ LUT grids
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)    :: drad_lut, dcfracwf_lut,dalbwf_lut
  REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:)  :: dozwf_lut, dtauwf_lut,draywf_lut
  ! correction spectrum @ user girds
  REAL (KIND=dp), ALLOCATABLE,DIMENSION(:), SAVE:: dfrad, dfcfracwf, dfalbwf
  REAL (KIND=dp), ALLOCATABLE,DIMENSION(:,:),SAVE:: dffozwf, dftauwf, dfraywf
  REAL (KIND=dp), ALLOCATABLE,DIMENSION(:,:), SAVE:: taucum, taucum2
  LOGICAL, SAVE :: first=.TRUE.
  do_raywf=.false.
  IF (do_polut_init) THEN 
    CALL Free_Ti0trsb
    L_stageJ = .TRUE.
    ! define the range of LUT to be used
    wlm(1:2)=REAL((/wave(1), wave(nw)/), KIND=4)
    !vza_min = vza ; vza_max=vza
    iv_min = get_index(DBLE(vza_min), DBLE(vza_grid(1:8)), MAX_LG_INTERP)
    iv_max = get_index(DBLE(vza_max), DBLE(vza_grid(1:8)), MAX_LG_INTERP)
    iv_max = MIN( 8, iv_max + MAX_LG_INTERP-1 )
    vzam= vza_grid((/iv_min, iv_max/))
    print * , 'LUT:vza', iv_min, iv_max
    iv_min = get_index(DBLE(sza_min), DBLE(sza_grid(1:12)), MAX_LG_INTERP)
    iv_max = get_index(DBLE(sza_max), DBLE(sza_grid(1:12)), MAX_LG_INTERP)
    iv_max = MIN( 12, iv_max + MAX_LG_INTERP-1 )
    szam= sza_grid((/iv_min, iv_max/))
    print * , 'LUT:sza', iv_min, iv_max
    DO ib = 1, 2 !! ib = 1, vecLUT; ib = 2, scaLUT
      CALL Init__VLDLUT( VLDLUTdir1, ib, LUT_type(ib), wlm(:), szam(:), &
                         vzam(:), L_stageJ, LFAIL, msg  )
      IF( LFAIL ) THEN
         WRITE(*,*) TRIM(msg) ;STOP 'Error'
      ENDIF
    ENDDO ! ib = 1, 2
    ib = 1; nwLUT=nwl__B(ib); npJ=nlyrsMAX !! from LUT_Stage_m
    IF (first) then
      ! Save variables @ user_grids
      allocate(dfrad(max_spec_pts), dfcfracwf(max_spec_pts),dfalbwf(max_spec_pts))
      allocate(dffozwf(max_spec_pts, mflay), dftauwf(max_spec_pts, mflay))
      allocate(taucum(max_spec_pts, 0:mflay))
      print * , ADJUSTL(TRIM(VLDLUTDir1)), do_raywf 
      allocate (wl(1:nwLUT),pJ(npJ),pL(0:npJ),taucum_lut(nwLUT,0:npJ), &
               Itoa(nwLUT,2), wta(nwLUT, 2), wtc(nwLUT,2), wto(nwLUT, npJ, 2))
      IF (do_raywf) THEN 
         allocate(dfraywf(max_spec_pts, mflay))
         allocate(taucum2(max_spec_pts, 0:mflay))
         allocate(taucum2_lut(nwLUT,0:npJ), wto2(nwLUT, npJ, 2))
      ENDIF 
      IF (do_plan) THEN 
       print *, "I am here 1"
       allocate (i0V(1:nwLUT,2), trV(1:nwLUT,2), sbV(1:nwLUT,2), & 
                i0S(1:nwLUT,2), trS(1:nwLUT,2), sbS(1:nwLUT,2))
       allocate (di0dtV(nwLUT,npJ,2), dtrdtV(nwLUT,npJ,2),dsbdtV(nwLUT,npJ,2), & 
                di0dtS(nwLUT,npJ,2), dtrdtS(nwLUT,npJ,2),dsbdtS(nwLUT,npJ,2))
       IF (do_raywf) THEN 
       allocate (di0dtV2(nwLUT,npJ,2), dtrdtV2(nwLUT,npJ,2),dsbdtV2(nwLUT,npJ,2), & 
                di0dtS2(nwLUT,npJ,2), dtrdtS2(nwLUT,npJ,2),dsbdtS2(nwLUT,npJ,2))
       ENDIF
      ELSE
        allocate(i0(1:nwLUT), tr(1:nwLUT), sb(1:nwLUT), &
             Itoa_ic(1:nwLUT), wta_ic(1:nwLUT), & 
             wtc_ic(1:nwLUT),  wto_ic(1:nwLUT, npJ))
        IF (do_raywf) THEN 
          PRINT * , 'not implemented in polcorr_lut'
          stop 1
        ENDIF
      ENDIF
    ENDIF
    wl(1:nwLUT) =  wlLUT_B(ib,1:nwLUT)
    do_polut_init = .false.
  ENDIF

  sfcp=ps(nsprs-1) ; ctp = ps(nctp-1)
  ilat = INT(ABS(lat/30.0)) + 1
  allocate (alb(nwLUT), alb1(nwLUT),alb2(nwLUT))
  CALL INTERPOL (wave(1:nw), albs(1:nw, 1), nw, wl(1:nwLUT), alb1(1:nwLUT),nwLUT, errstat)
  CALL INTERPOL (wave(1:nw), albs(1:nw, 2), nw, wl(1:nwLUT), alb2(1:nwLUT),nwLUT, errstat)

  ! calculate I and wfs from LUT
  ! do_plan = T calculated here using i0,tr,sb
  !           F calculated in Itoa_rpro
  IF (do_plan) THEN 
    IF (niter >= 0) THEN
     log_ret = .false. 
     IF (do_debug) WRITE (*,'(A, 3L,10f8.2)') 'call itoa_rpro',do_plan,do_cfracwf,do_albwf, vza, sza, lat, ozdu, sfcp, ctp
     print *, "I am here 2"
     i0v = 0.0 ; trv=0.0; sbv = 0.0 ; di0dtv=0.0; dtrdtv=0.0; dsbdtv=0.0
     i0s = 0.0 ; trs=0.0; sbs = 0.0 ; di0dts=0.0; dtrdts=0.0; dsbdts=0.0
     IF (do_raywf) THEN 
        di0dtv2=0.0 ; dtrdtv2=0.0; dsbdtv2=0.0
        di0dts2=0.0 ; dtrdts2=0.0; dsbdts2=0.0
     ENDIF
     DO ic = 2,1,-1
       IF (ic == 1) pt = sfcp
       IF (ic == 2) pt = ctp
       IF (do_raywf) THEN 
       ib = 1 ! vector
       CALL Itoa_rpro( ib, ilat, pt, ozdu, sza, vza, phi, &
                log_ret,  nwLUT, wl(:), i0V(:,ic), trV(:,ic), sbV(:,ic), &
                LFAIL, msg,npJ_k=npJ, pJ_k=pJ(:),& 
                di0dt_k=di0dtV(:,:,ic),dtrdt_k=dtrdtV(:,:,ic),dsbdt_k=dsbdtV(:,:,ic), &
                di0dt2_k=di0dtV2(:,:,ic),dtrdt2_k=dtrdtV2(:,:,ic),dsbdt2_k=dsbdtV2(:,:,ic))
       ib = 2 ! scalar
       CALL Itoa_rpro( ib, ilat, pt, ozdu, sza, vza, phi, &
                log_ret,  nwLUT, wl(:), i0S(:,ic), trS(:,ic), sbS(:,ic), &
                LFAIL, msg,npJ_k=npJ, pJ_k=pJ(:),&
                di0dt_k=di0dtS(:,:,ic),dtrdt_k=dtrdtS(:,:,ic),dsbdt_k=dsbdtS(:,:,ic), &
                di0dt2_k=di0dtS2(:,:,ic),dtrdt2_k=dtrdtS2(:,:,ic),dsbdt2_k=dsbdtS2(:,:,ic))
       ELSE
       ib = 1 ! vector
       CALL Itoa_rpro( ib, ilat, pt, ozdu, sza, vza, phi, &
                log_ret,  nwLUT, wl(:), i0V(:,ic), trV(:,ic), sbV(:,ic), &
                LFAIL, msg,npJ_k=npJ, pJ_k=pJ(:),& 
                di0dt_k=di0dtV(:,:,ic),dtrdt_k=dtrdtV(:,:,ic),dsbdt_k=dsbdtV(:,:,ic))
       ib = 2 ! scalar
       CALL Itoa_rpro( ib, ilat, pt, ozdu, sza, vza, phi, &
                log_ret,  nwLUT, wl(:), i0S(:,ic), trS(:,ic), sbS(:,ic), &
                LFAIL, msg,npJ_k=npJ, pJ_k=pJ(:),&
                di0dt_k=di0dtS(:,:,ic),dtrdt_k=dtrdtS(:,:,ic),dsbdt_k=dsbdtS(:,:,ic))
       ENDIF
     ENDDO
     ib = 1 ; ipt = tsi%ipt
     pL(0:npJ) = plevLUT(0:npJ, ipt)
     taucum_lut(1:nwLUT,0:npJ)= Tsi%wts(1)*Ti0trsb(ib,ipt,Tsi%io3p(1))%taucum(1:nwLUT,0:npJ) &
           +Tsi%wts(2)*Ti0trsb(ib,ipt,Tsi%io3p(2))%taucum(1:nwLUT,0:npJ)
     IF (do_raywf) THEN 
     taucum2_lut(1:nwLUT,0:npJ)= Tsi%wts(1)*Ti0trsb(ib,ipt,Tsi%io3p(1))%taucum2(1:nwLUT,0:npJ) &
           +Tsi%wts(2)*Ti0trsb(ib,ipt,Tsi%io3p(2))%taucum2(1:nwLUT,0:npJ)
     ENDIF
    ENDIF
     allocate (qr1(nwLUT), qr2(nwLUT))
     itoa=0.0; wta=0.0; wto=0.0;wtc=0.0
     IF (do_raywf) wto2 = 0.0  
     frac = 1.0- the_cfrac 
     IF (.NOT. do_cfracwf) frac = 1.0
     IF (.NOT. do_cfracwf .and. the_cfrac > 0.9) frac = 0.0
     qr1 = alb1/(1-alb1*sbV(:,1))  
     qr2 = alb2/(1-alb2*sbV(:,2))
     Itoa(:,1) = [i0V(:,1) +trV(:,1)*qr1(:)]*frac+ & 
                 [i0V(:,2) +trV(:,2)*qr2(:)]*(1-frac) 
     IF (do_cfracwf) wtc(:,1) = ([i0V(:,1) + trV(:,1)*qr1(:)]-[i0V(:,2)+trV(:,2)*qr2(:)])/Itoa(:,1)
     IF (do_albwf) wta(:,1) = ([trV(:,1)*qr1(:)/alb1(:)]*frac + & 
               [trV(:,2)*qr2(:)/alb2(:)]*(1-frac))/itoa(:,1)
     DO i  = 1, npJ
      wto(:,i,1) = [di0dtV(:,i,1) + qr1*dtrdtV(:,i,1) +trV(:,1)*qr1*qr1*dsbdtV(:,i,1)]*frac + &
                   [di0dtV(:,i,2) + qr2*dtrdtV(:,i,2) +trV(:,2)*qr2*qr2*dsbdtV(:,i,2)]*(1.-frac)
      wto(:,i,1) = wto(:,i,1)/itoa(:,1)
     ENDDO
     IF (do_raywf) THEN 
     DO i  = 1, npJ
      wto2(:,i,1) = [di0dtV2(:,i,1) + qr1*dtrdtV2(:,i,1) +trV(:,1)*qr1*qr1*dsbdtV2(:,i,1)]*frac + &
                    [di0dtV2(:,i,2) + qr2*dtrdtV2(:,i,2) +trV(:,2)*qr2*qr2*dsbdtV2(:,i,2)]*(1.-frac)
      wto2(:,i,1) = wto2(:,i,1)/itoa(:,1)
     ENDDO
     ENDIF

     qr1 = alb1/(1-alb1*sbS(:,1)) ;  qr2 = alb2/(1-alb2*sbS(:,2))
     Itoa(:,2) = [i0S(:,1) +trS(:,1)*qr1(:)]*frac+[i0S(:,2) +trS(:,2)*qr2(:)]*(1-frac) 
     IF (do_cfracwf) wtc(:,2) = ([i0S(:,1) + trS(:,1)*qr1(:)]-[i0S(:,2)+trS(:,2)*qr2(:)])/itoa(:,2)
     IF (do_albwf) wta(:,2) = ([trS(:,1)*qr1(:)/alb1(:)]*frac + & 
                               [trS(:,2)*qr2(:)/alb2(:)]*(1-frac))/Itoa(:,2)
     DO i = 1, npJ
      wto(:,i,2) = [di0dtS(:,i,1) + qr1*dtrdtS(:,i,1) +trS(:,1)*qr1**2*dsbdtS(:,i,1)]*frac+&
                   [di0dtS(:,i,2) + qr2*dtrdtS(:,i,2) +trS(:,2)*qr2**2*dsbdtS(:,i,2)]*(1.-frac)
      wto(:,i,2) = wto(:,i,2)/itoa(:,2)
     ENDDO
     IF (do_raywf) THEN 
     DO i = 1, npJ
      wto2(:,i,2) = [di0dtS2(:,i,1) + qr1*dtrdtS2(:,i,1) +trS(:,1)*qr1**2*dsbdtS2(:,i,1)]*frac+&
                    [di0dtS2(:,i,2) + qr2*dtrdtS2(2:,i,2)+trS(:,2)*qr2**2*dsbdtS2(:,i,2)]*(1.-frac)
      wto2(:,i,2) = wto2(:,i,2)/itoa(:,2)
     ENDDO
     ENDIF
    deallocate(qr1,qr2)
  ELSE
    IF (niter == 0) THEN 
    itoa=0.0; wta=0.0; wto=0.0;wtc=0.0
    DO ib = 1, 2 ! loop for vector and scalar
      DO ic =  2, 1, -1  ! loop for ic=2 cloud, ic = 1 clear
       IF (ic == 1) THEN 
        frac = 1.0 - the_cfrac ;  pt = sfcp ; alb=alb1
        IF (.NOT. do_cfracwf) frac = 1.0 
       ELSE
        frac = the_cfrac ; pt = ctp ; alb=alb2
       ENDIF
       IF (frac == 0.0) CYCLE
       CALL Itoa_rpro( ib, ilat, pt, ozdu, sza, vza, phi, &
                log_ret,  nwLUT, wl(:), i0(:), tr(:), sb(:),     &
                LFAIL, msg,  alb_k = alb(:), npJ_k=npJ, pJ_k=pJ(:),&
                Itoa_k = Itoa_ic(:),dlnIdR_k = wta_ic(:),  dlnIdt_k = wto_ic(:,:) )
       Itoa(:,ib) = Itoa(:,ib) + itoa_ic*frac
       wta(:,ib)  = wta(:,ib) + wta_ic(:)*frac
       wto(:,1:npJ,ib)  = wto(:,1:npJ,ib) + wto_ic(:, 1:npJ)*frac
      IF (ib==1 .and. ic==1 .and. do_debug) & 
        WRITE (*,'(A, 3L,10f8.2)') 'call itoa_rpro',do_plan,do_cfracwf,do_albwf, & 
                                    vza, sza, lat, ozdu
       IF( LFAIL ) THEN
            WRITE(*,*) TRIM(msg) ; stop 1
       ENDIF
       IF (frac /= 0.0) THEN 
           IF (ic == 2) THEN 
             wtc (:, ib) = Itoa_ic(:)  
           ELSE IF (ic == 1) THEN 
             wtc (:, ib) = (wtc(:,ib) - Itoa_ic(:))/Itoa(:,ib)
           ENDIF
       ENDIF
      ENDDO ! ic
    ENDDO ! ib
    ib = 1 ; ipt = tsi%ipt
    pL(0:npJ) = plevLUT(0:npJ, ipt)
    taucum_lut(1:nwLUT,0:npJ)= Tsi%wts(1)*Ti0trsb(ib,ipt,Tsi%io3p(1))%taucum(1:nwLUT,0:npJ) &
           +Tsi%wts(2)*Ti0trsb(ib,ipt,Tsi%io3p(2))%taucum(1:nwLUT,0:npJ)
    ENDIF
  ENDIF
  deallocate (alb, alb1,alb2)

  IF (niter == 0 .or. do_plan ) THEN 
    allocate(drad_lut(nwLUT), dalbwf_lut(nwLUT), dcfracwf_lut(nwLUT),&
               dtauwf_lut(nwLUT, npJ), dozwf_lut(nwLUT, npJ))
    ! calculate difference 
    drad_lut = 1.0 ;  dalbwf_lut = 1.0 ; dcfracwf_lut = 1.0 ; dozwf_lut = 1.0
    dtauwf_lut = 0.0
    IF (do_raywf) then  
       allocate(draywf_lut(nwLUT,npJ)) ; draywf_lut = 0.0
    ENDIF
    WHERE (Itoa(:,2) /= 0.0)
       drad_lut(:) =  Itoa(:,1)/Itoa(:,2) !ItoaV /ItoaS
    END WHERE
    IF (do_albwf) THEN 
     WHERE (wta(:,2) /= 0.0)
       dalbwf_lut(:) =wta(:,1)/wta(:,2)*drad_lut(:)
     END WHERE
    ENDIF
    IF (do_cfracwf) THEN 
     WHERE (wtc(:,2) /= 0.0)
       dcfracwf_lut(:) =wtc(:,1)/wtc(:,2)*drad_lut(:)
     END WHERE
    ENDIF
    DO i = 1, npJ
     WHERE (wto(:,i,2) /= 0.0  )
       dozwf_lut(:,i) = wto(:,i,1)/wto(:,i,2)*drad_lut(:)   !didtV(:,i)/didtS(:,i)
       dtauwf_lut(:,i)= wto(:,i,1)-wto(:,i,2) !didtV(:,i)/itoaV-didtS(:,i)/itoaS
     END WHERE
     IF (do_raywf) draywf_lut(:,i)= wto2(:,i,1)-wto2(:,i,2) !didtV(:,i)/itoaV-didtS(:,i)/itoaS
     !WRITE(*,'(i5,10e15.6)') i,  wto(nwlut, i, :), Itoa(nwlut,:), wta(nwlut,:), wtc(nwlut,:)
    ENDDO

  !-------------------------------------------------------------------
  ! LUT interpolation to user grids
  !------------------------------------------------------------------
    ! derive ratio or difference onto on-line RTM grids
    allocate (Jvert1(1:nw,1:npJ), Jvert2(1:nw, 1:npJ),Jvert3(1:nw, 0:npJ))
    IF (do_raywf) allocate(Jvert0(1:nw, 1:npJ),Jvert4(nw, 0:npJ))
    ! interpolation I w.r.t wavelength
    dffozwf = 1.0 ; dfalbwf=1.0; dfcfracwf=1.0
    DO i = 1, nw
      CALL getLwt(nwLUT, wl(1:nwLUT)*1.D0, wave(i), fidx, lidx, frac1, frac2) 
      ! log(Iv/Is)
      dfrad(i) = log(drad_lut(fidx)*frac1+drad_lut(lidx)*frac2)
      ! dI/da_v / dI/da_s
      dfalbwf(i) = dalbwf_lut(fidx)*frac1+dalbwf_lut(lidx)*frac2
      ! dlnI/dc_v - dlnI/dc_s
      dfcfracwf(i) = dcfracwf_lut(fidx)*frac1+dcfracwf_lut(lidx)*frac2
      ! dlnI/dtau_v - dlnI/dtau_s
      IF (do_raywf) Jvert0(i,1:npJ) = draywf_lut(fidx, 1:npJ)*frac1 + draywf_lut(lidx, 1:npJ)*frac2
      Jvert1(i,1:npJ) = dtauwf_lut(fidx, 1:npJ)*frac1 + dtauwf_lut(lidx, 1:npJ)*frac2
      ! dI/do3_v / dI/do3_s
      Jvert2(i,1:npJ) = dozwf_lut(fidx, 1:npJ)*frac1 + dozwf_lut(lidx, 1:npJ)*frac2
      Jvert3(i,0:npJ) = taucum_lut(fidx, 0:npJ)*frac1 + taucum_lut(lidx, 0:npJ)*frac2
      IF (do_raywf) &
      Jvert4(i,0:npJ) = taucum2_lut(fidx, 0:npJ)*frac1 + taucum2_lut(lidx, 0:npJ)*frac2
    ENDDO
    ! interpolation II w.r.t layers
    DO i = 1, nz
        CALL getLwt( npJ, LOG(pJ(1:npJ)), LOG(pmid(i)), fidx, lidx, frac1, frac2)
        dftauwf(1:nw, i) = Jvert1(1:nw, fidx)*frac1 + Jvert1(1:nw, lidx)*frac2
        IF (do_raywf) dfraywf(1:nw, i) = Jvert0(1:nw, fidx)*frac1 + Jvert0(1:nw, lidx)*frac2
        dffozwf(1:nw, i) = Jvert2(1:nw, fidx)*frac1 + Jvert2(1:nw, lidx)*frac2
    ENDDO
    DO i = 0, nz
        CALL getLwt( npJ+1, LOG(pL(0:npJ)), LOG(ps(i)),fidx, lidx, frac1, frac2)
        taucum(1:nw, i) = Jvert3(1:nw, fidx-1)*frac1 + Jvert3(1:nw, lidx-1)*frac2
        IF (do_raywf) & 
        taucum2(1:nw, i) = Jvert4(1:nw, fidx-1)*frac1 + Jvert4(1:nw, lidx-1)*frac2
    ENDDO
    deallocate(drad_lut, dalbwf_lut, dcfracwf_lut, dtauwf_lut, dozwf_lut)
    deallocate (Jvert1, Jvert2, Jvert3)
    IF (do_raywf) deallocate(Jvert0, Jvert4)
  ENDIF 
  !======================================================
  ! applying correction w.r.t rad, wf wrs alb,cfrac,foz
  !=======================================================
  adj_o3 = 0.0 ; adj_ray = 0.0
  DO i = 1, nw
      dtau (1:nz) = tabs(i, 1:nz) - (taucum(i,1:nz) - taucum(i,0:nz-1))
      adj_o3 = sum(dftauwf(i,1:nz)*dtau(1:nz))

      IF (do_raywf) THEN    
        dtau (1:nz) = tray(i, 1:nz) - (taucum2(i,1:nz) - taucum2(i,0:nz-1))
        adj_ray = sum(dfraywf(i,1:nz)*dtau(1:nz))
      ENDIF

      rad(i)       = rad(i)*exp(dfrad(i)  + adj_o3 + adj_ray)
      fozwf(i, 1:nz) = fozwf(i, 1:nz)*dffozwf(i,1:nz)
      if (do_albwf) albwf(i,1)   = albwf(i,1)*dfalbwf(i)
      if (do_cfracwf) cfracwf(i)   = cfracwf(i)*dfcfracwf(i)
      !IF (i == nw .and.  niter==0 ) & 
      !WRITE(*,'(A,L,i5,A,2f8.2,A,3e15.7, A,3e15.7, A, f6.2)') & 
      !'LUT:',do_plan, i,'dfrad=', exp(dfrad(i)),exp(dfrad(i)+adj_o3), & 
      !        'cf/alb/o3=',dfcfracwf(i),dfalbwf(i),dffozwf(i,nz), & 
      !        'tau=',sum(tabs(i,1:nz)), taucum(i, nz), dftauwf(i,nz),'ozdu=',ozdu
  ENDDO 
    
  IF(first) first=.false.
END SUBROUTINE polcorr_online_with_lut

SUBROUTINE polcorr_online(niter, which_polcorr, nw,  nz, nctp, ncbp, nsprs,nalb, faerlvl,   ncorr, polidxs, &
           do_fozwf, do_albwf, do_faerwf, do_faerswf,do_codwf, do_sprswf,do_fraywf, do_cfracwf, &
           wave,  rad, pwave, prad, tauwf, ptauwf, tabs, &
           albwf, palbwf,  fozwf,   pfozwf,   faerwf, pfaerwf, faerswf,pfaerswf, &
           fcodwf,pfcodwf, fsprswf, pfsprswf, fraywf, pfraywf, cfracwf,pcfracwf)
  USE OMSAO_precision_module
  USE OMSAO_parameters_module, ONLY : mflay
  USE OMSAO_variables_module,  ONLY : currloop
  USE ozprof_data_module,      ONLY : mpolcorr, malbwf
  USE OMSAO_errstat_module
  
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER, INTENT(IN) :: nw, nz, ncorr, niter, which_polcorr, nctp, ncbp, faerlvl, nsprs, nalb
  LOGICAL, INTENT(IN) :: do_fozwf, do_albwf, do_faerwf, do_faerswf, do_codwf, do_sprswf, do_fraywf, do_cfracwf

  INTEGER, DIMENSION (ncorr),           INTENT(IN)    :: polidxs
  REAL (KIND=dp), DIMENSION(nw, nz),    INTENT(IN)    :: tabs
  REAL (KIND=dp), DIMENSION(nw, nz),    INTENT(INOUT) :: fozwf, faerwf,   faerswf, fcodwf, fsprswf, fraywf, tauwf
  REAL (KIND=dp), DIMENSION(ncorr, nz), INTENT(IN)    :: pfozwf, pfaerwf, pfaerswf, pfcodwf, pfsprswf, pfraywf, ptauwf
  REAL (KIND=dp), DIMENSION(nw),        INTENT(IN)    :: wave
  REAL (KIND=dp), DIMENSION(nw),        INTENT(INOUT) :: rad, cfracwf
  REAL (KIND=dp), DIMENSION(nw, nalb),INTENT(INOUT)   :: albwf
  REAL (KIND=dp), DIMENSION(ncorr),       INTENT(IN)  ::  prad, pwave, pcfracwf
  REAL (KIND=dp), DIMENSION(ncorr,nalb),       INTENT(IN) :: palbwf
  ! Local variables
  INTEGER                                          :: i, j, fidx, lidx, errstat
  REAL (KIND=dp)                                   :: frac1, frac2, tmprad
  REAL (KIND=dp), DIMENSION(nw)                    :: tmpdf
  REAL (KIND=dp), DIMENSION(mpolcorr),         SAVE:: dfrad, dfcfracwf
  REAL (KIND=dp), DIMENSION(mpolcorr, malbwf), SAVE:: dfalbwf
  REAL (KIND=dp), DIMENSION(mpolcorr, mflay),  SAVE:: dftauwf
  REAL (KIND=dp), DIMENSION(mpolcorr, mflay),  SAVE:: dffozwf, dffaerwf, dffaerswf, dffcodwf, dffsprswf, dffraywf
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE:: fozwf_sav !(mw, mflay)
  REAL (KIND=dp), DIMENSION(:),   ALLOCATABLE, SAVE:: cfracwf_sav !(mw)
  REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE:: albwf_sav
  LOGICAL                                          :: newcorr = .true.
  LOGICAL, SAVE :: first=.TRUE.
  
  IF (first) THEN 
    WRITE (www_lun, *) 'npolcorr', ncorr
    WRITE (www_lun, *) 'polwave:', wave(polidxs)
    IF (ncorr > mpolcorr) THEN
     WRITE (*, *) 'polcorr_online : ncorr > mpolcorr'
    ENDIF
    first = .false.
  ENDIF
  ! Compute differences at positions where corrections are explicitly calcualted
  IF (niter == 0 .OR. which_polcorr == 3 .OR. which_polcorr == 5 ) THEN
     dfrad(1:ncorr) = 0.0D0

     ! Corrected in different way using DlnI/Dtau: 08/28/2008
     IF (newcorr) THEN
        WHERE ( prad(1:ncorr) /= 0.0 )
           dfrad(1:ncorr) = LOG(prad(1:ncorr))- LOG(rad(polidxs))
        ENDWHERE
     ELSE
        WHERE ( prad(1:ncorr) /= 0.0 )
           dfrad(1:ncorr) = (rad(polidxs) - prad(1:ncorr)) / prad(1:ncorr)
        ENDWHERE
     ENDIF

     IF ( which_polcorr /= 5 .OR. (which_polcorr == 5 .AND. (niter == 1 .OR. (niter == 0 .AND. currloop == 0))) ) THEN
     
     IF (do_albwf)   dfalbwf(1:ncorr, 1:nalb) = 0.D0
     dftauwf(1:ncorr, 1:nz) = 0.D0
     IF (do_fozwf)   THEN
        dffozwf(1:ncorr, 1:nz) = 0.D0
     ENDIF    
     IF (do_cfracwf) dfcfracwf(1:ncorr) = 0.D0
     IF (do_faerwf)  dffaerwf(1:ncorr, faerlvl:nz)  = 0.D0
     IF (do_faerswf) dffaerswf(1:ncorr, faerlvl:nz) = 0.D0
     IF (do_codwf)   dffcodwf(1:ncorr, nctp:ncbp)   = 0.D0
     IF (do_sprswf)  dffsprswf(1:ncorr, nsprs:nz)   = 0.D0
     IF (do_fraywf)  dffraywf(1:ncorr, 1:nz) = 0.D0
    
     
     IF ( do_albwf) THEN
        WHERE( palbwf(1:ncorr, 1:nalb) /= 0.0)
           dfalbwf(1:ncorr, 1:nalb) = (albwf(polidxs, 1:nalb) - palbwf(1:ncorr, 1:nalb)) /palbwf(1:ncorr, 1:nalb)
        ENDWHERE
     ENDIF

     IF ( do_cfracwf) THEN
        WHERE (pcfracwf(1:ncorr) /= 0.0)
           dfcfracwf(1:ncorr) = (cfracwf(polidxs) - pcfracwf(1:ncorr)) / pcfracwf(1:ncorr)
        ENDWHERE
     ENDIF
     
     IF ( do_fozwf ) THEN
        WHERE( pfozwf(1:ncorr, 1:nz) /= 0.0)
           dffozwf(1:ncorr, 1:nz) = (fozwf(polidxs, 1:nz) - pfozwf(1:ncorr, 1:nz)) / pfozwf(1:ncorr, 1:nz)
        ENDWHERE
        DO i = 1 , nz 
           !dftauwf(1:ncorr, i) = ptauwf(1:ncorr, i)/prad(1:ncorr) -tauwf(polidxs,i)/rad(polidxs)
           dftauwf(1:ncorr, i) = ptauwf(1:ncorr, i) -tauwf(polidxs,i)
        ENDDO
    ENDIF
     
     IF ( do_faerwf ) THEN
        WHERE( pfaerwf(1:ncorr, faerlvl:nz) /= 0.0 )
           dffaerwf(1:ncorr, faerlvl:nz) = (faerwf(polidxs, faerlvl:nz) &
                - pfaerwf(1:ncorr, faerlvl:nz)) / pfaerwf(1:ncorr, faerlvl:nz)
        ENDWHERE
     ENDIF
     
     IF ( do_faerswf ) THEN
        WHERE ( pfaerswf(1:ncorr, faerlvl:nz) /= 0.0 )
           dffaerswf(1:ncorr, faerlvl:nz) = (faerswf(polidxs, faerlvl:nz) &
                - pfaerswf(1:ncorr, faerlvl:nz)) / pfaerswf(1:ncorr, faerlvl:nz)
        ENDWHERE
     ENDIF
     
     IF ( do_codwf ) THEN
        WHERE ( pfcodwf( 1:ncorr, nctp:ncbp ) /= 0.0 )
           dffcodwf(1:ncorr, nctp:ncbp) = (fcodwf(polidxs, nctp:ncbp) - pfcodwf(1:ncorr, nctp:ncbp)) &
                / pfcodwf(1:ncorr, nctp:ncbp)
        ENDWHERE
     ENDIF
     
     IF ( do_sprswf ) THEN
        WHERE ( pfsprswf( 1:ncorr, nsprs:nz ) /= 0.0 )
           dffsprswf(1:ncorr, nsprs:nz) = (fsprswf(polidxs, nsprs:nz) - pfsprswf(1:ncorr, nsprs:nz)) &
                / pfsprswf(1:ncorr, nsprs:nz)
        ENDWHERE
     ENDIF  

     IF ( do_fraywf ) THEN
        WHERE( pfraywf(1:ncorr, 1:nz) /= 0.0)
           dffraywf(1:ncorr, 1:nz) = (fraywf(polidxs, 1:nz) - pfraywf(1:ncorr, 1:nz)) / pfraywf(1:ncorr, 1:nz)
        ENDWHERE
     ENDIF
     
     ENDIF
  ENDIF

    
  IF ( which_polcorr /= 5 .OR. (which_polcorr == 5 .AND. (niter == 1 .OR. (niter == 0 .AND. currloop == 0))) ) THEN

    fidx = polidxs(1)
    lidx = polidxs(ncorr)

    IF (do_albwf) THEN 
     DO i = 1 ,nalb
      CALL interpol (pwave, dfalbwf(1:ncorr, i), ncorr, wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
      albwf(fidx:lidx, i) = albwf(fidx:lidx, i) / (1.0 + tmpdf(fidx:lidx)) 
     ENDDO
    ENDIF

    IF (do_cfracwf) THEN 
      CALL interpol (pwave, dfcfracwf(1:ncorr), ncorr,wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
      cfracwf(fidx:lidx) = cfracwf(fidx:lidx) / (1.0 + tmpdf(fidx:lidx)) 
    ENDIF

    IF (do_fozwf) THEN
       DO i = 1, nz 
         CALL interpol (pwave, dffozwf(1:ncorr, i), ncorr,wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
        fozwf(fidx:lidx,i) = fozwf(fidx:lidx, i) / (1.0 + tmpdf(fidx:lidx)) 
       ENDDO
    ENDIF
    IF ( do_faerwf ) THEN
      DO i = faerlvl, nz
         CALL interpol (pwave, dffaerwf(1:ncorr, i), ncorr,wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
        faerwf(fidx:lidx,i) = faerwf(fidx:lidx, i) / (1.0 + tmpdf(fidx:lidx)) 
      ENDDO
    ENDIF
    IF ( do_faerswf ) THEN
      DO i = faerlvl, nz
         CALL interpol (pwave, dffaerswf(1:ncorr, i), ncorr,wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
        faerswf(fidx:lidx,i) = faerswf(fidx:lidx, i) / (1.0 + tmpdf(fidx:lidx)) 
      ENDDO
    ENDIF
     
    IF ( do_codwf )  THEN
      DO i = nctp, ncbp
         CALL interpol (pwave, dffcodwf(1:ncorr, i), ncorr,wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
        fcodwf(fidx:lidx,i) = fcodwf(fidx:lidx, i) / (1.0 + tmpdf(fidx:lidx)) 
      ENDDO
    ENDIF

    IF ( do_sprswf )  THEN
       DO i = nsprs, nz
         CALL interpol (pwave, dffsprswf(1:ncorr, i), ncorr,wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
        fsprswf(fidx:lidx,i) = fsprswf(fidx:lidx, i) / (1.0 + tmpdf(fidx:lidx)) 
       ENDDO
    ENDIF

    IF ( do_fraywf  ) THEN
       DO i = 1, nz
         CALL interpol (pwave, dffraywf(1:ncorr, i), ncorr,wave(fidx:lidx), tmpdf(fidx:lidx),lidx-fidx+1, errstat)
        fraywf(fidx:lidx,i) = fraywf(fidx:lidx, i) / (1.0 + tmpdf(fidx:lidx)) 
       ENDDO
    ENDIF
  ENDIF

  fidx = polidxs(1)
  DO i = 2, ncorr
     lidx = polidxs(i)
     !print * , fidx, lidx
     DO j = fidx, lidx
        frac2 = (wave(j) - pwave(i-1)) / (pwave(i) - pwave(i-1))
        frac1 = 1.0 - frac2
        IF (newcorr) THEN
           tmprad = frac1 * (SUM(dftauwf(i-1, 1:nz) * (tabs(j, 1:nz)-tabs(polidxs(i-1), 1:nz))) + dfrad(i-1)) + &
                frac2 * (SUM(dftauwf(i, 1:nz) * (tabs(j, 1:nz)-tabs(polidxs(i), 1:nz))) + dfrad(i))
           rad(j) = rad(j) * EXP(tmprad)
           !print * , i, exp(tmprad), exp(dfrad(i))
           !write(www_lun,*) i,j, tmprad, wave(j), dfrad(i-1), dftauwf(i-1,10),tabs(j, 10), tabs(polidxs(i-1), 10)
        ELSE
           rad(j) = rad(j) / (1.0 + frac1 * dfrad(i-1) + frac2 * dfrad(i))
        ENDIF
     ENDDO
     fidx = lidx + 1
  ENDDO

  !print * , 1/(1+dffozwf(1:ncorr, nz))
  !print * , 1/(1+dfcfracwf(1:ncorr))
  !print * , 1/(1+dfalbwf(1:ncorr,1))
  IF (which_polcorr == 5) THEN
     IF (niter == 1) THEN
        IF (allocated(fozwf_sav)) deallocate (fozwf_sav, albwf_sav,cfracwf_sav)
        allocate (fozwf_sav(nw, nz), albwf_sav(nw, malbwf), cfracwf_sav(nw))
        fozwf_sav(1:nw, 1:nz) = fozwf(1:nw, 1:nz)
        albwf_sav(1:nw, 1:nalb) = albwf(1:nw, 1:nalb)
        cfracwf_sav(1:nw) = cfracwf(1:nw)
     ELSE IF (niter > 1 .OR. (niter == 0 .AND. currloop /= 0) ) THEN
        fozwf(1:nw, 1:nz) = fozwf_sav(1:nw, 1:nz)
        albwf(1:nw, 1:nalb) = albwf_sav(1:nw, 1:nalb)
        cfracwf(1:nw) = cfracwf_sav(1:nw)
     ENDIF
  ENDIF
  RETURN
  
  END SUBROUTINE polcorr_online

  ! Obtain minor trace gas weighting functions from ozwf
  ! Note: it is not the exact wf but negative of the WF divided by radiances
  ! SUBROUTINE GET_TRACEGAS_WF (ozwf, ozabs, so2crs, use_so2dtcrs, rad, nw, nz, nz1, &
  !     ozs, waves, do_so2zwf, so2zwf)
  ! Note : both nz and nz1 should be there
  SUBROUTINE get_tracegas_wf (nw, nz, nz1, rad, ozwf, do_so2zwf, so2zwf)
  USE OMSAO_precision_module
  USE OMSAO_parameters_module,ONLY : du2mol
  USE OMSAO_indices_module,   ONLY : so2_idx, so2v_idx, o2o2_idx, o2_idx, h2o_idx,  &
                                     o2t2_idx, h2ot2_idx
  USE OMSAO_variables_module, ONLY : refidx, database, database_save, refspec_norm
  USE ozprof_data_module,     ONLY : nlay, mgasprof, fgasidxs, &
      tracegas, ngas, gasidxs, fgassidxs, so2valts, so2vprofn1p1, trace_profwf,& 
      nup2p, use_lograd, use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs,ccrs
  IMPLICIT NONE

  ! Input/Output variables
  INTEGER, INTENT(IN)                           :: nw, nz, nz1
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(IN) :: ozwf
  REAL (KIND=dp), DIMENSION(nw),     INTENT(IN) :: rad !, waves
  LOGICAL,                           INTENT(IN) :: do_so2zwf
  REAL (KIND=dp), DIMENSION(nw),    INTENT(OUT) :: so2zwf
  !Local variables
  LOGICAL                                        :: do_wf = .false.
  INTEGER                                        :: i, j, k, fidx, lidx, nk
  REAL (KIND=dp)                                 :: tmp, tmpsum
  REAL (KIND=dp), ALLOCATABLE, DIMENSION (:,:)   :: amf, tmpcrs
  REAL (KIND=dp), ALLOCATABLE, DIMENSION (:,:)   :: tamf
  REAL (KIND=dp)                                 :: avcd, svcd
 
  ! Obtain AMF  each wavelength and at each layer
  IF (ANY(fgasidxs > 0)) THEN ! dlnI/dtau
     allocate (amf(nw, nz), tmpcrs(nw, nz), tamf(ngas, nw))
     DO i = 1, nz1
        amf(:, i) = -ozwf(:, i) / rad(:) / ccrs%o3(:, i) / du2mol
     ENDDO
     amf (:,nz1+1:nz) = 0.0D0 ; tmpcrs(:,nz1+1:nz) = 0.0
     do_wf = .true.
  ENDIF 
    
  IF (do_wf ) THEN
  ! Replace cross sections with weighting functions
  DO i = 1, ngas 
     IF (fgasidxs(i) > 0) THEN
        avcd = mgasprof(i, nz+1)
        IF ( ((gasidxs(i) /= so2_idx .AND. gasidxs(i) /= so2v_idx)  .OR. .NOT. use_so2dtcrs) .AND. &
             (gasidxs(i) /= o2o2_idx  .OR. .NOT. use_o4dtcrs) .AND. &
             ((gasidxs(i) /= h2o_idx .AND. gasidxs(i) /= h2ot2_idx) .OR. .NOT. use_h2odptcrs) .AND. &
             ((gasidxs(i) /= o2_idx .AND. gasidxs(i) /=o2t2_idx)    .OR. .NOT. use_o2dptcrs) ) THEN
           DO j = 1, nw
              tamf(i, j) = SUM(amf(j, 1:nz1) * mgasprof(i, 1:nz1)) / avcd
           ENDDO
           
           IF (fgassidxs(i) > 0) THEN
              database(gasidxs(i), refidx(1:nw)) = database(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
           ELSE
              database(gasidxs(i), refidx(1:nw)) = database_save(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
              !database(gasidxs(i), refidx(1:nw)) = database(gasidxs(i), refidx(1:nw)) * tamf(i, 1:nw)
           ENDIF
          
           !xliu, 3/31/2015, weighted air mass factor by cross sections
           !tracegas(i, 7) = 0.0; nk = 0
           !DO j = 1, nw
           !   IF (database_save(gasidxs(i), refidx(j)) > 0.) THEN
           !      tracegas(i, 7) = tracegas(i, 7) + tamf(i, j)
           !      nk = nk + 1
           !   ENDIF
           !ENDDO
           !tracegas(i, 7) = tracegas(i, 7) / nk
           !print * , tracegas(i, 7)

           tracegas(i, 7) = 0.0; tmpsum = 0.0d0
           DO j = 1, nw
              IF (database_save(gasidxs(i), refidx(j)) > 0.) THEN
                 tracegas(i, 7) = tracegas(i, 7) + database(gasidxs(i),refidx(j))
                 tmpsum = tmpsum + database_save(gasidxs(i), refidx(j))
              ENDIF
           ENDDO
           tracegas(i, 7) = tracegas(i, 7) / tmpsum
           !print * , tracegas(i,7)
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
           ENDDO           
        ELSE
           tmp = avcd * refspec_norm(gasidxs(i))
           IF ( gasidxs(i) == so2_idx .OR. gasidxs(i) == so2v_idx) THEN 
               tmpcrs(1:nw,1:nz1) = ccrs%so2(1:nw, 1:nz1)
           ELSE IF (gasidxs(i) == o2o2_idx) THEN
               tmpcrs(1:nw,1:nz1) = ccrs%o4(1:nw, 1:nz1)
           ELSE IF (gasidxs(i) == o2_idx .OR. gasidxs(i) == o2t2_idx) THEN 
               tmpcrs(1:nw,1:nz1) = ccrs%o2(1:nw, 1:nz1)
           ELSE IF (gasidxs(i) == h2o_idx .OR. gasidxs(i) == h2ot2_idx) THEN
               tmpcrs(1:nw,1:nz1) = ccrs%h2o(1:nw, 1:nz1)
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
                 so2zwf(j) = SUM(ozwf(j, 1:nz1) / ccrs%o3(j, 1:nz1) * ccrs%so2(j, 1:nz1) &
                      * (so2vprofn1p1(1:nz1, 2)-so2vprofn1p1(1:nz1, 1)) ) / tmp
              ENDDO
           ELSE
              tmp =tmp / refspec_norm(gasidxs(i))
              DO j = 1, nw
                 so2zwf(j) = SUM(ozwf(j, 1:nz1) / ccrs%o3(j, 1:nz1) * database_save(gasidxs(i), refidx(j)) &
                      * (so2vprofn1p1(1:nz1, 2)-so2vprofn1p1(1:nz1, 1)) ) / tmp
              ENDDO
           ENDIF
        ENDIF  
     ENDIF
  ENDDO
   deallocate(amf, tmpcrs, tamf)
  ENDIF ! do_wf
  RETURN
  END SUBROUTINE get_tracegas_wf

  SUBROUTINE debug_rtm (funit, rtm, ccrs)
  !USE OMSAO_indices_module, ONLY: o2o2_idx, o2_idx, h2o_idx
  USE ozprof_data_module, ONLY: rtm_outputs, ccrs_set !,o4crsidx, o2crsidx, h2ocrsidx
  IMPLICIT NONE
  INTEGEr, INTENT(IN) :: funit
  TYPE (rtm_outputs), INTENT(IN) :: rtm
  TYPE (ccrs_set),    INTENT(IN) :: ccrs
  INTEGER :: i,nw,nz, ngas, nalb
  nw = rtm%nw
  nz = rtm%nl
  nalb = rtm%nalbwf
  ngas = 1
  WRITE(funit, *) nw, ngas, nalb
  WRITE(funit, '(A)') '@wave,  alb, rad, albwf, ozwf, o3crs'
  DO i = 1, nw
     WRITE(funit, '(F10.4, 1000D16.7)') rtm%wav(i), rtm%alb(i), rtm%rad(i, 1), &
     rtm%albwf(i, 1:nalb,1), rtm%ozwf(i, nz, 1), ccrs%o3(i,nz), ccrs%o4(i,nz),ccrs%o2(i, nz), ccrs%h2o(i, nz)
  ENDDO
  END SUBROUTINE debug_rtm
  
  SUBROUTINE debug_taug (funit, nw, nz, ngas, wav, allcol, allcrs)
  USE OMSAO_precision_module, ONLY: dp
  USE OMSAO_indices_module, ONLY: refspec_strings
  USE ozprof_data_module, ONLY: num_iter, gasidxs, fgaspos
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: funit, nw, nz, ngas
  REAL(KIND=dp), INTENT(IN) :: wav(nw), allcol(ngas, nz), allcrs(nw, ngas, nz)
  INTEGER :: i, j
  REAL (KIND=dp) :: tmp (ngas, nw )
  WRITE(funit, *) num_iter, nw, nz, ngas
  WRITE(funit, '(10a10)') 'o3  ',refspec_strings(gasidxs(fgaspos(1:ngas-1)))
  tmp = 0.0
  DO i = 1, ngas
    DO j = 1, nz 
      tmp(i, 1:nw) = tmp(i, 1:nw) + allcol(i, j)*allcrs(1:nw, i, j)
    ENDDO
  ENDDO
  DO i = 1, nw
     WRITE(funit, '(i5,F10.4, 1000D16.7)') i, wav(i), tmp(1:ngas, i),allcrs(i,:, 10)
  ENDDO
  WRITE(funit, *) 'gas profile'
  DO i = 1, nz 
     WRITE(funit, *) allcol(:, i)
  ENDDO
  END SUBROUTINE

  SUBROUTINE get_efft(nz, zs, ozs, fts, ts, errstat)
  USE OMSAO_precision_module
  USE OMSAO_errstat_module

  INTEGER, INTENT(IN)                              :: nz
  INTEGER, INTENT(OUT)                             :: errstat                                                   

  REAL (KIND=dp), DIMENSION(nz),     INTENT(IN)    :: ozs
  REAL (KIND=dp), DIMENSION(0:nz),   INTENT(IN)    :: zs, fts
  REAL (KIND=dp), DIMENSION(1:nz),   INTENT(OUT)   :: ts

  INTEGER :: i, j, fidx, lidx, nz1
  REAL (KIND=dp) :: dz
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
  FUNCTION get_index( x0, xarray, NLG ) RESULT( idx )
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xarray
        REAL(KIND=8), INTENT(IN) :: x0
        INTEGER, INTENT(IN) :: NLG
        INTEGER :: nn
        INTEGER :: ifoo
        INTEGER :: idx

        nn = SIZE( xarray )
        IF( x0 <= xarray(1) ) THEN
           idx = 1
           RETURN
        ELSE IF( x0 >= xarray(nn) ) THEN
           idx = nn -NLG+1
           RETURN
        ENDIF
        ifoo = MINLOC( xarray, DIM = 1, MASK = xarray - x0 >= 0.d0 )
        idx = ifoo - NLG/2
        IF( idx < 1 ) THEN
           idx = 1
        ELSE IF( idx > nn-NLG+1 ) THEN
           idx = nn-NLG+1
        ENDIF
      END FUNCTION get_index
END MODULE m_lidort_util
