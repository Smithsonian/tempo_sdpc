!>Main memory allocation subroutine
module m_alloc1

contains

  !-------------------------------------------------------------------------
  ! !ROUTINE:  alloc1
  ! 
  !  DESCRIPTION:
  !> alloc1 allocates/deallocates memory for retrievals
  !
  ! !CALLING SEQUENCE: 
  !
  !        call alloc1
  !     
  ! !INPUT PARAMETERS:   
  !> @param errstat error reporting integer, non-zero = failure
  !
  ! !REVISION HISTORY: 
  !
  !> @author  05Jan01   Joiner      original fortran 90
  !> @author  26Mar15   O'Sullivan  updated for TEMPO
  !
  !-------------------------------------------------------------------------
  subroutine alloc1(errstat)

    use m_vars, ONLY: nchl, chl, npres, & 
         nscanpos, nTimes, wave_resid, resid, write_resid, nwl
    use m_cloud_pres_mod
    use tell_module
    implicit none

    integer, intent(inout) :: errstat

    if (errstat /= 0) return


    !deallocate memory
    !=================
    if (allocated(ind)) then
      deallocate(ind, waves, wavesp, w1p, f1p, wave_diff, sflx, &
           bad_obs, good_obs, r_i, y_obs, y_obs1, y_frac, y_calc_squeeze, &
           y_calc, y_calc_sh, y_resid, y_back, jacob_dummy, jacob_rad, &
           fit_rad, rad_cld, rad_clr, rad_clr_oc, ring_cld, ring_clr, &
           ring_oc, comp_all, comp_all_ring, computed, computed_oc, &
           comp_clr, comp_oc_clr, rad_tot, ycalc, wavesd, o3_xsect, &
           chls, ntot, nia_clds, nra_clds, tra_clds, i0a_clds, &
           z1_clds, z2_clds, i01a_clds, ring_clds, rad_clds, rad_clrs, &
           ring_clrs, ring_ocs, &
           stat=errstat)

      if (write_resid) then
        deallocate(wave_resid, resid, stat=errstat)
      endif

      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "alloc1: deallocation failure", &
             errstat)
        return
      endif

    endif ! allocated


    allocate(ind(0:nobs-1), &
         waves(0:nobs-1), &
         wavesp(nterms-1,0:nobs-1), &
         w1p(nwl), &
         f1p(nwl), &
         wave_diff(0:nobs-1), &
         sflx(0:nobs-1), &
         bad_obs(nobs), &
         good_obs(nobs), &
         r_i(0:nobs-1), &
         y_obs(0:nobs-1)  , &
         y_obs1(0:nobs-1)  , &
         y_frac(0:nobs-1)  , &
         y_calc(0:nobs-1)  , &
         y_calc_sh(0:nobs-1)  , &
         y_calc_squeeze(0:nobs-1)  , &
         y_resid(0:nobs-1,1)  , &
         y_back(0:nobs-1,1)  , &
         jacob_dummy(0:nobs-1)  , &
         jacob_rad(0:nobs-1)  , &
         fit_rad(0:nobs-1)  , &
         rad_cld(0:nobs-1)  , &
         rad_clr(0:nobs-1)  , &
         rad_clr_oc(0:nobs-1)  , &
         ring_cld(0:nobs-1)  , &
         ring_clr(0:nobs-1)  , &
         ring_oc(0:nobs-1)  , &
         computed_oc(nchl), &
         computed(npres), &
         comp_all(npres), &
         comp_all_ring(npres), &
         chls(0:nchl-1), &
         rad_tot(0:nobs-1)  , &
         ycalc(0:nobs-1)  , &
         wavesd(0:nobs-1), &
         o3_xsect(0:nobs-1), &
         ntot(npres,0:nobs-1), &
         nia_clds(npres,0:nobs-1), &
         nra_clds(npres,0:nobs-1), &
         tra_clds(npres,0:nobs-1), &
         i0a_clds(npres,0:nobs-1), &
         z1_clds(npres,0:nobs-1), &
         z2_clds(npres,0:nobs-1), &
         i01a_clds(npres,0:nobs-1), &
         ring_clds(npres,0:nobs-1), &
         rad_clds(npres,0:nobs-1), &
         rad_clrs(npres,0:nobs-1), &
         ring_clrs(npres,0:nobs-1), &
         comp_oc_clr(npres), &
         comp_clr(npres), &
         ring_ocs(nchl,0:nobs-1), &
         stat=errstat)

    chls=chl

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "alloc1: allocation failure", &
           errstat)
      return
    endif

    if (write_resid) then
      allocate(wave_resid(nobs), resid(nobs,nscanpos,nTimes), stat=errstat) 
      resid=-9999.
      wave_resid = 0.
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "alloc1: allocation failure for resid or wave_resid", &
             errstat)
        return
      endif
    endif


  end subroutine alloc1

end module m_alloc1
