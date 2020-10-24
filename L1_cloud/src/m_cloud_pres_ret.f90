!>Main processing routine for cloud pressure and cloud fraction calculations
module m_cloud_pres_ret

  private

  public cloud_pres_ret

contains

  !!==================================================================
  !
  ! Subroutine cloud_pres_ret:
  !
  !> @brief
  !> This is the main code performing the retrieval (via Chi^2 fit) 
  !> of the primary cloud code data products:
  !>
  !> Cloud Pressure for O3\n
  !> Cloud Fraction for O3\n
  !>
  !> It also calculates a number of other output variables 
  !> (in some cases via calls to other modules) including:
  !>
  !> Convergence Factor\n
  !> dIdR (Fractional Sensitivity of Radiance to Reflectivity)\n
  !> Effective Filling-in\n
  !> Radiative Cloud Fraction\n
  !> Reflectivity\n
  !> Residual Bias and Residual Standard Deviation\n
  !> Wavelength Shift\n
  !>
  !> It can also update:
  !>
  !> Measurement Quality (flags)\n
  !> Cloud Mask\n
  !>
  !> There are a number of options which can generate additional outputs,
  !> including:
  !>
  !> Wavelength Squeeze\n
  !> 
  !> INPUTS/OUTPUTS:
  !>
  !> The two nominal input parameters refl_clr and refl_cld represent the 
  !> reflectivity of the earth in clear or fully overcast conditions.
  !> 
  !> For documentation of all other parameters used, SEE:\n
  !> m_cloud_pres_mod\n
  !> m_vars \n
  !> and other modules called herein (m_get_f, m_get_ai_refl, etc.)\n
  !>
  !> REFERENCES;
  !>
  !> Primary (most useful):\n
  !> Joiner et al. (1995) Applied Optica, 34, 4513\n
  !> Joiner & Vasilkov (2006) IEEE Trans. Geo. and Rem. Sens., 44, 1272\n
  !>
  !> Secondary:\n
  !> Joiner & Bhartia (1995) J. of Geophys. Res., 100, 23109\n
  !> Joiner et al. (2011) Atmos. Meas. Teach. Dscuss., 4, 6186\n
  !> Sneep et al. (2008) J. of Geophys. Res., 113, D15S23\n
  !> Optional Chlorophyll retrieval:\n
  !> Joiner et al. (2004) J. of Geophy. Res., 109, D01109\n
  !
  ! AUTHORS:
  ! 
  !> @author Joiner & Vasilkov - original OMI code (probably c2001)
  !> @author O'Sullivan (2014) - updates and documentation for TEMPO
  !
  ! VARIABLES:
  !
  !> @param refl_clr surface reflectivity of cloud free pixels
  !> @param refl_cld reflectivity of completely cloud filled pixels
  !> @param errstat Error return code, non-zero indicates problem
  !
  !!====================================================================
  subroutine cloud_pres_ret(refl_clr, refl_cld, errstat)

    use mathcons
    use m_get_ai_refl
    use m_get_f
!    use m_invert
    use m_invert2
    use m_find
    use m_findgen
    use m_indgen
    use m_poly_fit
    use m_trilin
    use m_matmul
    use m_interpol
    use m_bilinear
    use m_spline
    use m_sigma
    use m_interp_ring_rad
    use m_interp_pres
    use m_alloc1
    use m_alloc2
    use m_lambda_qual
    use m_pgs_include
    use tell_module
    use m_cloud_pres_mod
    use m_read_input_data_gems, only: bad_rad_lambda_gems
    use m_vars, ONLY: azimuth, bad_obs_flag, biases, & 
         chi_sqr, chlcl, chlorophyll, cld_frac_min, cld_pres2, &
         cloud_fr_corr, cloud_mask, cloud_pres, do_alloc, do_chl, &
         do_LER, do_mler, do_o3, do_short_wave, eff_cld_frac, &
         eff_cld_frac2, f12d, fill, fill_value, fs, geoflg, &
         get_cloud_frac, get_refl_clim, iLine, land_flg, lat, &
         ler354, ler_nph, ler_nsz, ler_nth, ler_ph, ler_sz, ler_th, lon, &
         meas_qual_flg, nchl, n_good_input, n_good_output, niter, &
         n_missing, noret, no_ret_ps, nphi, npres, nscan, nscan_oc, &
         ntheta, nthet_oc, nwl, nXtrack, oc_table, phi, pres, ps, qc, &
         rad_cld_frac, ref_clr, refl, refl_clr_oc, reflect_cld, resid, &
         resid_spec, sat_zen, scan, scan_oc, shift, shifts, squeeze, &
         squeezes, stds, sza, theta, theta_oc, using_resid, & 
         using_spline, w12d, wave_fill, wave_long, wave_o3, wave_resid, &
         wave_short, wdelt, w_grid, wmax, wmin, write_fill, write_obs, & 
         write_resid, ws, xsect_o3, test_solar, add_shift, do_cloud_mask, &
         read_gems
    implicit none

    real (KIND=8), intent(inout) :: refl_clr
    real (KIND=8), intent(inout) :: refl_cld
    integer, intent(inout) :: errstat

    !local vraiables
    character (len=128) :: logmsg
    integer :: polyfit_status

    if (errstat /= 0) return

    !*************************************************************************
    ! m_cloud_pres_ret is called within a loop in OMCLDRR. If this is the 
    ! first iteration, allocate memory and set up wavelength arrays
    if (do_alloc) then 

      ! get wavelengths for retrieval
      !==============================
      nobs=count(w_grid >= wmin+wdelt .and. w_grid <= wmax-wdelt)   

      !allocate memory
      !===============
      call alloc1(errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "alloc1: failure", &
             errstat)
        return
      endif
      if (nobs > 0) then
        ind=find2(w_grid >= wmin+wdelt .and. w_grid <= wmax-wdelt,nobs)   
      else
        write(logmsg,"(A40,2X,F7.3,2X,F7.3)") &
             'cloud_pres_ret: min and max wavelengths ',wmin, wmax
        call tell_log(1,logmsg)
!        if (iprt >= 1) write(6,'(7f12.4)') w_grid
        errstat = -1
        call tell_error(tell_unknown_error, &
             "cloud_pres_ret: no valid wavelengths", errstat)
        return
      endif ! nobs > 0

      !get wavelengths to use in retrieval and initialize arrays
      !retrieval done on subset of table wavelengths!
      !=========================================================
      waves=w_grid(ind)   
      wavesd=waves - waves(0)
      do ntm=2,nterms
        wavesp(ntm-1,:)=wavesd**ntm ! wavelength^n
      enddo

      if (do_o3) then
        o3_xsect=interpol(xsect_o3,wave_o3,waves)  
      endif

      do_alloc=.false.
    endif ! do_alloc

    ! If we're actually doing the retrievals, rather than testing IO..
    ! ================================================================
    if (.not. noret) then

      !f1p and w1p need to be allocated in each iteration to allow for
      !variation in number of valid wavelengths (nwl)
      if(allocated(f1p)) deallocate(f1p, stat=errstat)
      if(allocated(w1p)) deallocate(w1p, stat=errstat)
      allocate(f1p(nwl), w1p(nwl), stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "cloud_pres_ret: allocation failure", &
             errstat)
        return
      endif


      !profile loop
      !============
      do ip=0, nXtrack-1 

        !do some initialization
        !=======================
        sz=sza(ip,iLine)
        !switch to allow use of solar spectrum as input, for testing
        !===========================================================
        if (.not. test_solar) then
          f1p = f12d(0:nwl-1,ip)
        else
          f1p = fs(:,ip)
        endif

        !check for missing or bad geolocation data
        !=========================================
        if (lat(ip+1,iLine) < -90. .or. lat(ip+1,iLine) > 90.  &
             .or. lon(ip+1,iLine) < -180. .or. lon(ip+1,iLine) > 180.  &
             .or. sz > sz_max .or. sz < theta(1) .or. maxval(f1p) < 0. &
             ) qc(ip,iLine)=IBSET(qc(ip,iLine),1)

        ! Add check for extreme viewing zenith angles
        if (sat_zen(ip,iLine) > sat_zen_max) qc(ip,iLine)=IBSET(qc(ip,iLine),1)

        if (BTEST(geoflg(ip+1,iLine),6)) qc(ip,iLine)=IBSET(qc(ip,iLine),15)

        bad_pix = BTEST(qc(ip,iLine),15) .or. &
             BTEST(meas_qual_flg(iLine),1) .or. &
             BTEST(qc(ip,iLine),14) .or. &
             BTEST(qc(ip,iLine),1)

        if (bad_pix) then
          n_missing = n_missing+1
          cycle
        endif

        !preliminary estimate of number of good input pixels
        !===================================================
        n_good_input = n_good_input + 1

        !do more initialization
        !=======================
        r_i=0. ! set observation errors to zero
        !if testing with solar spectrum allow articfical wavelength shift
        !================================================================
        if (.not. test_solar) then
          w1p = w12d(0:nwl-1,ip)
        else
          w1p = ws(:,ip)-add_shift
        endif

        !check for bad radiances
        !=======================
        if (.not. read_gems) then
          call bad_rad_lambda(ip, iLine, errstat)
        else
          call bad_rad_lambda_gems (ip, iLine, errstat)
        endif
        if (check_solar) then
          if (btest(qc(ip,iLine),11)) then
            do i=0, nobs-1
              counts=count(abs(ws(:,ip) - waves(i)) < wavetol)
              if (counts > 0) then

                if (allocated(indw)) deallocate(indw, stat=errstat)
                allocate(indw(counts), stat=errstat)
                if (errstat /= 0) then
                  call tell_error (tell_malloc_error, &
                       "cloud_pres_ret: allocation failure", &
                       errstat)
                  return
                endif

                indw=find2(abs(ws(:,ip) - waves(i)) < wavetol,counts)
                if (any(fs(indw-1,ip) == 0.)) then
                  r_i(i)=var_inv_big
                endif ! if bad data found

              endif! if wavelengths found
            enddo ! loop over wavelengths
          endif ! if bad solar data found
        endif ! check_solar

        if (check_rad) then
          if (btest(qc(ip,iLine),9)) then
            do i=0, nobs-1
              counts=count(abs(w1p - waves(i)) < wavetol)
              if (counts > 0) then

                if (allocated(indw)) deallocate(indw, stat=errstat)
                allocate(indw(counts), stat=errstat)
                if (errstat /= 0) then
                  call tell_error (tell_malloc_error, &
                       "cloud_pres_ret: allocation failure", &
                       errstat)
                  return
                endif

                indw=find2(abs(w1p - waves(i)) < wavetol,counts)
                if (any(f1p(indw) == 0.)) then
                  r_i(i)=var_inv_big
                endif ! bad data found

              endif ! wavelengths found
            enddo ! loop over wavelengths
          endif ! check_rad
        endif ! check_rad

        !interpolate solar flux to table wavelengths 
        !===========================================
        if (using_spline) then
          ngood=count(fs(:,ip) > 0)
          if (ngood /= size(fs(:,ip))) then
            if (ngood > 0) then
              if (allocated(good)) deallocate(good, stat=errstat)
              allocate(good(ngood), stat=errstat)
              if (errstat /= 0) then
                call tell_error (tell_malloc_error, &
                     "cloud_pres_ret: allocation failure", &
                     errstat)
                return
              endif
              good=find2(fs(:,ip) > 0,ngood)-1
              sflx=spline(ws(good,ip),fs(good,ip),waves)
            else
              qc(ip,iLine)=IBSET(qc(ip,iLine),1)
            endif
          else
            sflx=spline(ws(:,ip),fs(:,ip),waves)
          endif

        else ! not spline, use linear interpolation instead
          ngood=count(fs(:,ip) > 0)
          if (ngood /= size(fs(:,ip))) then
            if (ngood > 0) then
              if (allocated(good)) deallocate(good, stat=errstat)
              allocate(good(ngood), stat=errstat)
              if (errstat /= 0) then
                call tell_error (tell_malloc_error, &
                     "cloud_pres_ret: allocation failure", &
                     errstat)
                return
              endif
              good=find2(fs(:,ip) > 0,ngood)-1
              sflx=interpol(fs(good,ip),ws(good,ip),waves)
            else
              qc(ip,iLine)=IBSET(qc(ip,iLine),1)
            endif
          else
            sflx=interpol(fs(:,ip),ws(:,ip),waves)
          endif

        endif ! not spline

        bad_pix = BTEST(qc(ip,iLine),1) 
        if (bad_pix) then
          n_missing = n_missing+1
          cycle
        endif

        ! get number of elements in state vector
        !========================================
        chloro=chlcl(ip)
        add_oc = .false.
        nst=nterms+2
        nchlr=0
        if (shift) then
          nst=nst+1
          nsh=1
          if (squeeze) nst=nst+1
        else
          nsh=0
        endif
        if (do_o3) then
          nst=nst+1
        endif

        !allocate memory
        !===============
        call alloc2(errstat)
        if (errstat /= 0) then
          call tell_error (tell_malloc_error, &
               "alloc2: failed, exiting", &
               errstat)
          return
        endif

        !profile independent Jacobian
        !============================
        h(:,nst-1-nsh)=1.0 ! constant   
        h(:,nst-2-nsh)=wavesd ! wavelength
        do ntm=2,nterms
          h(:,nst-1-nsh-ntm)=wavesp(ntm-1,:) ! wavelength^n
        enddo
        if (do_o3) then
          h(:,nst-nsh-nterms-nchlr-2)=o3_xsect
        endif

        ! define background error covariance matrix
        !==========================================
        if (add_background_error) then 
          b_i(nst-1-nsh)=var_inv_big ! const
          do ntm=1,nterms
            b_i(nst-1-nsh-ntm)=var_inv_big ! wavelength^n
          enddo
          if (shift) then
            b_i(nst-1)=var_inv_big ! shift
            if (squeeze) then
              b_i(1)=var_inv_big
            endif
          endif
          if (do_o3) then
            b_i(nst-nterms-nsh-2-nchlr)=var_inv_big
          endif
        endif

        satz=sat_zen(ip,iLine)
        az=azimuth(ip,iLine)
        psurf=ps(ip,iLine)
        if (get_refl_clim) refl_clr=ref_clr(ip,iLine)

        if(.not.land_flg(ip) .and. refl_clr .lt. refl_ice) then
          nt=interpol(findgen(ler_nsz)+1,ler_sz,sz)
          if(nt < 0) call tell_log(1,'negative nt in interpolation of ler sza')
          i1_ler=int(nt)
          i2_ler=i1_ler+1
          j=interpol(findgen(ler_nth)+1,ler_th,satz)
          l=interpol(findgen(ler_nph)+1,ler_ph,az)
          if(j < 0 .or. l < 0) call tell_log(1, &
               'negative input to bilinear interpolation of ler')
          !per Joanna Joiner, set a minium of j=1
          if(j < 1) j=1.0
          int1_ler=bilinear(ler354(i1_ler,:,:),j,l)
          int2_ler=bilinear(ler354(i2_ler,:,:),j,l)
          refl_clr=refl_clr+(nt-i1_ler)*(int2_ler-int1_ler)+int1_ler
        endif

        ! get interpolate values for viewing geometry
        !============================================
        nt=interpol(findgen(ntheta)+1,theta,sz)
        j=interpol(findgen(nscan)+1,scan,satz)
        l=interpol(findgen(nphi)+1,phi,az)
        np0=interpol(findgen(npres)+1,pres,psurf)
        !  np0=minval((/np0,float(npres)/))
        np0=minval((/np0,real(npres, KIND=8)/))
        nt_o=interpol(findgen(nthet_oc)+1,theta_oc,sz)
        j_o=interpol(findgen(nscan_oc)+1,scan_oc,satz)
        i_np0=nint(np0)

        !bracket the surface pressure
        !=================================
        i0x1=int(np0)
        ixd=abs(i0x1-np0)   
        if (ixd < delx) then 
          if (i0x1 < npres-1 .and. (i0x1 < np0 .or. i0x1 == 1)) then 
            i0x1=i0x1+1 
          else 
            i0x1=i0x1-1   
          endif
        endif
        i0x2=i0x1+1

        !construct observation vector
        !============================
        ngood=count(w1p > 0. .and. f1p > 0.)
        if (ngood > 0) then
          if (allocated(good)) deallocate(good, stat=errstat)
          allocate(good(ngood), stat=errstat)
          if (errstat /= 0) then
            call tell_error (tell_malloc_error, &
                 "cloud_pres_ret: allocation failure", &
                 errstat)
            return
          endif
          good=find2(w1p > 0. .and. f1p > 0.,ngood)
          if (using_spline) then
            y_obs=spline(w1p(good),f1p(good),waves)/sflx 
          else ! not spline
            y_obs=interpol(f1p(good),w1p(good),waves)/sflx !*pi
          endif ! not spline
          y_obs1=y_obs*pi
        else ! no good data found
          qc(ip,iLine)=IBSET(qc(ip,iLine),1)
        endif ! no good data found

        bad_pix = BTEST(qc(ip,iLine),1) 
        if (bad_pix) then
          n_missing = n_missing+1
          cycle
        endif

        !check for good residuals and apply correction if required
        !=========================================================
        if (using_resid) then
          if (size(resid_spec,dim=1) /= nobs) then
            call tell_error(tell_io_error, &
                 "cloud_pres_ret: incompatible resid table, not using corrections", &
                 errstat)
          else ! good residuals found
            y_obs=y_obs-resid_spec(:,ip+1)*y_obs
            y_obs1=y_obs1-resid_spec(:,ip+1)*y_obs1
          endif ! good residuals found
        endif ! using_resid

        !initialize computed arrays to false
        !===================================
        comp_clear=.false.
        comp_all=.false.
        comp_all_ring=.false.
        computed=.false.
        comp_clr=.false.
        comp_oc_clr=.false.
        bias=99999.
        nbad=0

        !construct first guess
        !=====================
        !res3(1:nterms+1)=poly_fit(wavesd,y_obs,nterms)
        polyfit_status = 0
        call poly_fit (wavesd, y_obs, nterms, res3(1:nterms+1), polyfit_status)
        if (polyfit_status /= 0) then
          qc(ip,iLine)=IBSET(qc(ip,iLine),4)
          cloud_pres(ip,iLine) = fill_value
          cld_pres2(ip,iLine) = fill_value
          cycle
        endif
        if (shift) then
          x(nst-1,1)=0.01 ! wavelength shift
          if (squeeze)   then
            x(1,1)=1.001 !0.99998 GOME ! wavelength squeeze factor (prof 574)
          endif
        endif
        x(nst-1-nsh,1)=res3(1) ! constant
        do ntm=1, nterms
          x(nst-1-nsh-ntm,1)=res3(ntm+1) ! wavelength scale factor
        enddo
        x(0,1)=0.7 ! Cloud pressure in atmospheres
        if (no_ret_ps) x(0,1)=psurf

        !set initial values for iterative loop
        !=====================================
        iter=0   
        chisq_old=9999.   
        diff_chi=9999.   
        ix1_old=9999
        ic1_old=9999

        !Construct observation error covariance
        !======================================
        do i=0, nobs-1    
          if (r_i(i) == 0) then
            if (waves(i) < wave_min) then   
              r_i(i)=1./(y_obs(i)*noise2)**2   
            else   
              r_i(i)=1./(y_obs(i)*noise)**2   
            endif
          endif
        enddo

        !set qc if snow/ice
        !====================
        NISE = IBITS( geoflg(ip+1, iLine), 8, 7)
        if( NISE >= 50 .and. NISE <= 103) qc(ip, iLine) = IBSET(qc(ip, iLine),5)

        if (no_ret_ps .and. .not. shift .and. .not. squeeze) then
          niters = 1
        else
          niters = niter
        endif

        !iteration loop
        !==============
        iteration_loop: &
          do while ( iter < niters .and. diff_chi > diff_chi_max .and. .not. &
             (btest(qc(ip,iLine),0) .or. btest(qc(ip,iLine),4) &
             .or. btest(qc(ip,iLine),bad_obs_flag)) )

          !bracket the cloud pressure between min and max
          !==============================================
          if(x(0,1) < pres(1)) then
            x(0,1)=pres(1)   
          elseif(x(0,1) > pres(npres)) then
            x(0,1)=pres(npres)   
          endif

          !bracket the cloud pressure index
          !=================================
          np=interpol(findgen(npres)+1,pres,x(0,1))   
          ix1=int(np)
          ixd=abs(ix1-np)   
          if (ixd < delx) then 
            if (ix1 < npres-1 .and. (ix1 < np .or. ix1 == 1)) then 
              ix1=ix1+1 
            else 
              ix1=ix1-1   
            endif
          endif
          if (ix1 >= npres) ix1=npres-1
          if ((np < ix1 .and. ix1 > 1) .or. ix1 == npres) then
            ix1=ix1-1
          endif
          ix2=ix1+1

          !interpolate observation to wavelengths for reflectivity calculation
          !put shift in opposite direction
          !===================================================================
          if (iter == 0 .or. cld_frac == 1 .or. cld_frac == 0) then
            ! at iteration 0 compute reflectivity at 500 hPa pressure
            ! unless no_ret_ps = true, then use surface pressure
            ! interpolate table values
            !=====================================================
            !JJ need to make sure this is not a bad pixel!!!!!!!
            !===================================================

            ! Assumes the reflectivity wavelength is within the fitting window
            !==================================================================
            ngood=count(r_i /= var_inv_big)
            if (ngood > 0) then
              good_obs(1:ngood)= &
                   find2(r_i /= var_inv_big,ngood)-1
              i_r=find1(abs(waves(good_obs(1:ngood))-wave_long) == &
                   minval(abs(waves(good_obs(1:ngood))-wave_long)) )
              indt=good_obs(i_r)
              ind0=ind(indt)
              call interp_rad(ix1,ind0, i0_l1, sb_l1, tr_l1)
              call interp_rad(ix2,ind0, i0_l2, sb_l2, tr_l2)
              call interp_rads(ix1, ix2, pres, x(0,1), i0_l1, i0_l2, sb_l1, &
                   sb_l2, tr_l1, tr_l2, i0_l, sb_l, tr_l)
              if (iter == 0) then
                call interp_rad(i0x1,ind0, i0_l1, sb_l1, tr_l1)
                call interp_rad(i0x2,ind0, i0_l2, sb_l2, tr_l2)
                call interp_rads(i0x1, i0x2, pres, psurf, i0_l1, i0_l2, sb_l1,&
                     sb_l2, tr_l1, tr_l2, i0_ls, sb_ls, tr_ls)
              endif
              i_obs_l=y_obs1(indt)

              if (do_short_wave) then 
                i_r=find1(abs(waves(good_obs(1:ngood))-wave_short) == &
                     minval(abs(waves(good_obs(1:ngood))-wave_short)) )
                indt=good_obs(i_r)
                ind0=ind(indt)
                call interp_rad(ix1,ind0, i0_s1, sb_s1, tr_s1)
                call interp_rad(ix2,ind0, i0_s2, sb_s2, tr_s2)
                call interp_rads(ix1, ix2, pres, x(0,1), i0_s1, i0_s2, sb_s1,&
                     sb_s2, tr_s1, tr_s2, i0_s, sb_s, tr_s)
                if (iter == 0) then
                  call interp_rad(i0x1,ind0, i0_l1, sb_l1, tr_l1)
                  call interp_rad(i0x2,ind0, i0_l2, sb_l2, tr_l2)
                  call interp_rads(i0x1, i0x2, pres, psurf, i0_l1, i0_l2, &
                       sb_l1, sb_l2, tr_l1, tr_l2, i0_ss, sb_ss, tr_ss)
                endif
                i_obs_s=y_obs1(indt)
              endif
              set_cld_frac=iter == 0
              if (do_mler) then 
                call get_ai_refl(refl_clr, refl_cld, I_obs_l, I_obs_s, ip, &
                     i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, sb_ls, &
                     tr_ls, set_cld_frac, i0_ss,sb_ss,tr_ss)
              else
                call get_f(refl_clr, refl_cld, I_obs_l, I_obs_s, ip, &
                     i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, sb_ls, &
                     tr_ls, set_cld_frac, i0_ss,sb_ss,tr_ss)
                reflect_cld(ip,iLine)=real(refl_cld, kind=4)
              endif
              reflec=refl(ip,iLine)
              if (reflec > 1) then
                reflec = 0.99999
              elseif (reflec < 0) then
                reflec = 0.00001
              endif

              if (iter == 0) then
                !set effective cloud fraction to 1 if retrieving 
                !effective scene pressure
                !===============================================
                if (.not. get_cloud_frac) then
                  cld_frac=1.
                else
                  cld_frac=eff_cld_frac(ip,iLine)
                endif

                !set cloud fraction, if clear or overcast, set fraction to 1
                !revert to LER
                !============================================================
                if (cld_frac <= cld_frac_min) then
                  cld_frac=1.
                  qc(ip,iLine)=IBSET(qc(ip,iLine),13)
                  x(0,1)=psurf
                endif

                !Snow Ice, revert to LER approximation
                !=====================================
                if (btest(qc(ip,iLine),5) .or. do_LER) then
                  cld_frac=1.
                  x(0,1)=psurf
                  if (.not. do_LER) &
                       eff_cld_frac(ip,iLine)=1.
                endif

                !set first guess to current values
                !=================================
                x_fg=x

                !Do correction based on reflectivity (cloud fraction?)
                !====================================================
                cld_frac_oc = eff_cld_frac(ip,iLine)
                add_oc = add_oc .or. (.not. land_flg(ip) .and. do_chl .and. &
                     cld_frac_oc < 1.) 
                if (add_oc) then
                  computed_oc=.false.
                else
                  computed_oc=.true.
                endif
              endif ! iter == 0

              !set the cloud reflectivity
              !==========================
              if (cld_frac == 1.) then
                reflec_cld=reflec
              else
                reflec_cld=real(refl_cld, kind=4)
              endif

            else  ! no good pixels
              qc(ip,iLine)=IBSET(qc(ip,iLine),bad_obs_flag)
              cycle
            endif ! good pixels

          endif ! iter == 0 or cld_frac == 1

          !set the ocean reflectivity. If clear, set to retrieved reflectivity,
          !else, set to predefined ocean clear reflectivity
          !====================================================================
          if (add_oc .and. cld_frac_oc == 0.) then
            refl_oc=reflec
          else
            refl_oc=refl_clr_oc
          endif


          b_i(0)=var_inv_big


          !bracket the chlorophyll index
          !=================================
          if (add_oc) then
            nc=interpol(findgen(nchl)+1,chls,chloro)   
            ic1=int(nc)
            if (ic1 > nchl-1) then
              ic1 = nchl-1
            elseif (ic1 < 1) then
              ic1=1
            endif
            icd=abs(ic1-nc)   
            if (icd < delc) then 
              if (ic1 < nchl-1 .and. (ic1 < nc .or. ic1 == 1)) then 
                ic1=ic1+1 
              else 
                ic1=ic1-1   
              endif
            endif
          else
            ic1=1
            ic2=1
          endif

          !compute radiance and atmospheric Ring at lower cloud press. bound
          !=================================================================
          if (cld_frac == 1) computed=.false. ! recompute at new reflectivity
          call interp_ring_rad(ix1,reflec_cld,computed,rad_clds(ix1,:),&
               ring=ring_clds(ix1,:), errstat=errstat)
          call interp_ring_rad(ix2,reflec_cld,computed,rad_clds(ix2,:),&
               ring=ring_clds(ix2,:), errstat=errstat)

          if (add_oc .and. iter == 0) then
            call interp_ring_rad(i_np0,refl_oc,comp_oc_clr,rad_clr_oc, &
                 errstat=errstat)
          endif

          if ( (eff_cld_frac2(ip,iLine) < 1.0 .or. & !cld_frac > cld_frac_min .and. &
               cld_frac < 1.0) .and. get_cloud_frac) then
            call interp_ring_rad(i0x1,refl_clr, comp_clr, rad_clrs(i0x1,:), &
                 ring=ring_clrs(i0x1,:), errstat=errstat)
            call interp_ring_rad(i0x2,refl_clr, comp_clr, rad_clrs(i0x2,:), &
                 ring=ring_clrs(i0x2,:), errstat=errstat)
          endif ! get_cloud_frac

          if (errstat /= 0) then
            call tell_error (tell_malloc_error, &
                 "interp_ring_rad: failed, exiting", &
                 errstat)
            return
          endif

          !compute ocean filling lower chlorophyll boundary
          !================================================
          if (add_oc) then
            if (.not. computed_oc(ic1)) then
              temp3D => oc_table (:,:,ic1,ind(0):ind(nobs-1))
              ring_ocs(ic1,:) = bilin(nt_o,j_o)
              computed_oc(ic1) = .true.
            endif ! computed_oc(ic1)
          endif ! add_oc

          !compute ocean filling at upper chlorophyll boundary
          !===================================================
          ic2=ic1+1
          if (.not. computed_oc(ic2) .and. add_oc) then
            temp3D => oc_table (:,:,ic2,ind(0):ind(nobs-1))
            ring_ocs(ic2,:) = bilin(nt_o,j_o)
            computed_oc(ic2) = .true.
          endif ! computed(ic2) and add_oc

          !interpolate at each wavelength
          !==============================
          temp2D=>ring_clds(ix1:ix2,:)
          call interp_pres(ix1, ix2, ring_cld, pres, x(0,1), h(:,0))
          temp2D=>rad_clds(ix1:ix2,:)
          call interp_pres(ix1, ix2, rad_cld, pres, x(0,1),jacob_rad)
          if ( (eff_cld_frac2(ip,iLine) < 1.0 .or. & 
                                !cld_frac > cld_frac_min .and. &
               cld_frac < 1.0) .and. get_cloud_frac .and. .not. comp_clear) then
            temp2D=>ring_clrs(i0x1:i0x2,:)
            call interp_pres(i0x1, i0x2, ring_clr, pres, psurf, jacob_dummy)
            temp2D=>rad_clrs(i0x1:i0x2,:)
            call interp_pres(i0x1, i0x2, rad_clr, pres, psurf, jacob_dummy)
            comp_clear=.true.
          endif ! get_cloud_frac
          if (add_oc .and. (iter == 0)) then
            temp2D=>ring_ocs(ic1:ic2,:)
            call interp_pres(ic1, ic2, ring_oc, chls, chloro, jacob_dummy)
            new_hcl=ic1_old /= ic1 
          endif

          !update cloud pressure jacobian if necessary
          !===========================================
          new_h=ix1_old /= ix1 
          !print *,'contributions to jacobian from ring '
          !print *,h(:,0)*rad_cld
          h(:,0)=(h(:,0)*rad_cld)*cld_frac
          if (do_o2_jacob) then
            !fit a line to first and last points of computed radiance
            !=========================================================
            fit_rad=((jacob_rad(nobs-1)-jacob_rad(0))/(nobs-1))*findgen(nobs)+jacob_rad(0)

            !subtract the line from the computed Jacobian 
            !(accounted for in polynomial)
            !===============================================
            h(:,0)=h(:,0)+(jacob_rad-fit_rad)*(1+ring_cld)*cld_frac
          endif

          !compute cloudy radiance as function of cloud fraction
          !-----------------------------------------------------
          rad_tot=cld_frac*rad_cld*(1+ring_cld)
          if ( cld_frac < 1.0 .and. get_cloud_frac) then
            rad_tot= rad_tot + (1-cld_frac)*rad_clr*(1+ring_clr)
          endif
          if (add_oc) then
            if (cloud_fr_corr .or. cld_frac_oc == 0.) then
              if (cld_frac_oc < 1.) then
                rad_tot_oc=(1-cld_frac_oc)!*rad_clr_oc
              endif
            else
              rad_tot_oc=(1-reflec)**2!*rad_clr_oc
            endif ! cloud_fr_corr
            rad_tot= rad_tot + rad_tot_oc*rad_clr_oc*ring_oc
          endif ! add_oc

          !fit polynomial to computed radiance and subtract 
          !------------------------------------------------
          !res=poly_fit(wavesd,rad_tot,nterms,yfit=ycalc)
          call poly_fit (wavesd,rad_tot,nterms,res,polyfit_status,yfit=ycalc)
          if (polyfit_status /= 0) then
            qc(ip,iLine)=IBSET(qc(ip,iLine),4)
            diff_chi = 0.0
            cloud_pres(ip,iLine) = fill_value
            cld_pres2(ip,iLine) = fill_value
            exit iteration_loop
          endif
          y_calc_sh=x(nst-1-nsh,1)+x(nst-2-nsh,1)*wavesd
          do ntm=2,nterms
            y_calc_sh=y_calc_sh+x(nst-1-nsh-ntm,1)*wavesp(ntm-1,:)
          enddo
          y_calc_sh=y_calc_sh*(((rad_tot-ycalc)/ycalc)+1.)

          if (shift) then
            !apply wavelength shift
            !----------------------
            if (using_spline) then
              y_calc=spline(waves,y_calc_sh*sflx,waves+x(nst-1,1))/sflx   
            else
              y_calc=interpol(y_calc_sh*sflx,waves,waves+x(nst-1,1))/sflx   
            endif

            !compute the Jacobian for wavelength shift
            !-----------------------------------------
            h(:,nst-1)=(y_calc-y_calc_sh)/(x(nst-1,1))   
            wave_diff=waves+x(nst-1,1)

            !compute the Jacobian for wavelength squeeze
            !Assumes shift = .true.
            !-------------------------------------------
            if (squeeze) then
              adj=waves(0)*x(1,1)-waves(0)

              ! this one shifts then squeezes shifted wavelengths
              !==================================================
              if (using_spline) then
                y_calc_squeeze=spline(waves,y_calc_sh*sflx, &
                     (waves*x(1,1) - adj + x(nst-1,1)) ) /sflx   
              else
                y_calc_squeeze=interpol(y_calc_sh*sflx, waves, &
                     (waves*x(1,1) - adj + x(nst-1,1)) ) /sflx   
              endif
              h(:,1)=(y_calc_squeeze-y_calc)/(x(1,1)-1.0)   
              y_calc=y_calc_squeeze
              wave_diff=waves*x(1,1) - adj + x(nst-1,1)
            endif
          else
            wave_diff=waves
            y_calc=y_calc_sh
          endif

          y_resid(:,1)=y_obs-y_calc   
          y_frac=y_resid(:,1)/y_obs
          new_r=.false.
          if (iter > 0 .and. nbad == 0) then
            nbad=count(abs(y_frac) > bad_thresh .and. r_i /= var_inv_big)
            if (nbad > 0) then
              bad_obs(1:nbad)= &
                   find2(abs(y_frac) > bad_thresh .and. r_i /= var_inv_big,nbad)-1
              r_i(bad_obs(1:nbad))= var_inv_big
              new_r=.true.
              chisq_old=9999.   
            endif ! nbad > 0
          endif ! iter > 0

          !load transpose of Jacobian and multiply by
          !observation error covariance (diagonal)
          !speed up by only recomputing for columns of
          !H that are changing
          !===========================================
          if (iter == 0 .or. new_r) then
            do ii=0, nobs-1
              htr(:,ii)=h(ii,:)*r_i(ii)
            enddo ! ii
          else
            if (new_h) then
              htr(0,:)=h(:,0)*r_i   
            endif
            if (shift) then
              htr(nst-1,:)=h(:,nst-1)*r_i
              if (squeeze) htr(1,:)=h(:,1)*r_i
            endif
          endif

          !compute retrieval error covariance and next
          !iteration state estimate
          !can speed this up by only computing on diagonal
          !===============================================
          err_cov=(htr .mm. h) 
          if (add_background_error) then
            do ii=0, nst-1
              err_cov(ii,ii) = err_cov(ii,ii) + b_i(ii)
            enddo ! ii
          endif
          err_cov = invert2(err_cov,ierr)

          if (ierr == 0) then ! good retrieval
            y_back(:,1)=b_i * (x(:,1) - x_fg(:,1))   
            x = x + (err_cov .mm. ((htr .mm. y_resid) + y_back))

            !compute chisq estimate
            !-----------------------
            chisq=sum(y_resid(:,1)**2*r_i)
            chisq=chisq**0.5/(nobs)   

            !compute radiance bias and standard deviation
            !============================================
            bias_old=bias
            std=sigma(y_resid(:,1)*r_i**0.5,avg1=bias)   

            !compute residual in terms of fractional amount
            !==============================================
            y_resid(:,1)=y_frac !y_resid(:,1)/y_obs

            !update for next iteration
            !=========================
            diff_chi=(chisq_old-chisq)/chisq_old
            if (diff_chi < diff_chi_iter_max .and. &
                 abs(bias) > abs(bias_old) .and. .not. new_r) then
              qc(ip,iLine)=IBSET(qc(ip,iLine),0)
            endif
            chisq_old=chisq   
            ix1_old=ix1
            ic1_old=ic1
            iter=iter+1   

          else
            qc(ip,iLine)=IBSET(qc(ip,iLine),4)
            diff_chi=0.
            cloud_pres(ip,iLine)=fill_value
            cld_pres2(ip,iLine)=fill_value
            !eff_cld_frac(ip,iLine)=fill_value
            !eff_cld_frac2(ip,iLine)=fill_value
            !rad_cld_frac(ip,iLine)=fill_value
          endif ! matrix inversion check
          !print out results
          !=================
!          if (iprt >= 3) call print_ret()

        enddo iteration_loop    ! iter loop

        !print out retrieval error covariance
        !====================================
        ! print *, 'standard deviations'   
!        if (iprt >= 4) then
!          do i=0, size(err_cov(:,0))-1 
!            write(6,'(7e12.4)') err_cov(i,i)**0.5   
!          enddo ! i
!
!          corr=err_cov
!          do ii=0, nst-1    
!            do jj=0, nst-1
!              corr(ii,jj)=corr(ii,jj)/err_cov(ii,ii)**0.5/err_cov(jj,jj)**0.5
!            enddo ! jj   
!          enddo ! ii   
!          print *, 'error correlations'   
!          do jj=0, nst-1
!            write(6,'(7f9.3)') corr(:,jj)
!          enddo ! jj   
!        endif ! iprt >= 4

        !store cloud pressure
        !====================
        if (.not. do_LER)  then
          cloud_pres(ip,iLine)=real(x(0,1), kind=4)
        else
          cloud_pres(ip,iLine)=real( &
               (x(0,1)-psurf*(1-rad_cld_frac(ip,iLine)))/rad_cld_frac(ip,iLine) &
               , kind=4)
        endif
        if (shift) then
          shifts(ip,iLine)=real(x(nst-1,1), kind=4)
        endif
        if (squeeze) then
          squeezes(ip,iLine)=real(x(1,1), kind=4)
        endif
        chlorophyll(ip,iLine)=real(chlcl(ip), kind=4)
        biases(ip,iLine)=real(bias, kind=4)
        stds(ip,iLine)=real(std, kind=4)
        chi_sqr(ip,iLine)=real(chisq, kind=4)

        if (count(r_i == var_inv_big) >= nobs/2) &
             qc(ip,iLine)=IBSET(qc(ip,iLine),bad_obs_flag)

        !convergence check
        !=================
        if (iter == niters .and. diff_chi > diff_chi_max .and. ierr == 0 &
             .and. abs(bias) > abs(bias_old) .and. .not. new_r) then
          qc(ip,iLine)=IBSET(qc(ip,iLine),0)
        endif

        if(cloud_pres(ip,iLine) >    psurf) qc(ip,iLine)=IBSET(qc(ip,iLine),3)
        if(cloud_pres(ip,iLine) <  pres(1)) qc(ip,iLine)=IBSET(qc(ip,iLine),2)
        if(cloud_pres(ip,iLine) > fill_value) cloud_pres(ip,iLine)=cloud_pres(ip,iLine)*1013.25

        !update the cloud mask, if no contrast and high reflectivity
        !============================================================
        !to do, add more QC checks
        if ((do_cloud_mask) .and. (allocated(cloud_mask))) then
          if (cloud_mask(ip+1,iLine) == 0 .and. reflec > refl_cld_mask) then
            if (.not.(btest(qc(ip,iLine),5))) then
              cloud_mask(ip+1,iLine) = 1
            else
              !JJ comment out check over snow/ice for now
              !if (cloud_pres(ip,iLine) - psurf .lt. cld_mask_press_diff) then
              !  cloud_mask(ip+1,iLine) = 1
              !else
              !  cloud_mask(ip+1,iLine) = 3
              !endif
            endif
          endif
        endif

        !preliminary estimate of number of good retrievals
        if( .not. ( btest(qc(ip,iLine),0) .or. btest(qc(ip,iLine),2) &
             .or. btest(qc(ip,iLine),3) .or. btest(qc(ip,iLine),4) &
             .or. btest(qc(ip,iLine),6))) &
                                !.or. btest(qc(ip,iLine),7) .or. btest(qc(ip,iLine),8))) &
             n_good_output = n_good_output + 1

        !print final result
        !==================
        write(logmsg,'(i6,i3,f8.3,4e11.3,3f8.2,i12)') ip, iter, &
             cloud_pres(ip,iLine), x(nst-1,1), bias, std, &
             chlorophyll(ip,iLine), chlcl(ip), reflec
        call tell_log(2,logmsg)
        !write file for residuals if write_resid
        !======================================
        if (write_resid) then
          if (wave_resid(1) == 0) then
            wave_resid=real(waves, kind=4)
          endif
          if (write_obs) then
            resid(:,ip+1,iLine)=real(y_obs, kind=4)
          else
            resid(:,ip+1,iLine)=real(y_resid(:,1), kind=4)
          endif
        endif

        !Write effective filling in a given wavelength
        !=============================================

        if (write_fill) then
          ind1 = find1 (abs(waves - wave_fill) == minval(abs(waves - wave_fill)) ) - 1
          fill(ip,iLine) = &
               real(cld_frac*rad_cld(ind1)*(ring_cld(ind1)), kind=4)
          elastic=cld_frac*rad_cld(ind1)
          if ( cld_frac < 1.0 .and. get_cloud_frac) then
            fill(ip,iLine) = fill(ip,iLine) + &
                 real((1-cld_frac)*rad_clr(ind1)*(ring_clr(ind1)), kind=4)
            elastic=elastic +  (1-cld_frac)*rad_clr(ind1)
          endif
          fill(ip,iLine) = fill(ip,iLine)/real(elastic, kind=4)
          !print *,'filling-in wavelength ',waves(ind1), fill(ip,iLine)
        endif

      enddo    ! profile loop
    else ! noret=.true. so we're only testing IO
      fill(:,iLine)=0.
      cloud_pres(:,iLine)=0.5
      cld_pres2(:,iLine)=0.5
      do ip=0, nXtrack-1 
        chlorophyll(ip,iLine)=real(chlcl(ip), kind=4)
      enddo
      refl(:,iLine)=0.1
      rad_cld_frac(:,iLine)=1.
      eff_cld_frac(:,iLine)=1.
      eff_cld_frac2(:,iLine)=1.
      biases(:,iLine)=9999.
      stds(:,iLine)=9999.
      chi_sqr(:,iLine)=9999.
    endif ! if (.not. noret)


    ! ***********************************************************
  end subroutine cloud_pres_ret
end module m_cloud_pres_ret
