module m_alloc1

contains

  subroutine alloc1()

    use m_vars, ONLY: interp, nchl, chl, npres, w12d, f12d, &
         reference_ring, reference_rad, nscanpos, nTimes, wave_resid, resid, &
         write_resid, nwl
    use m_cloud_pres_mod
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  alloc1
    ! 
    ! !DESCRIPTION: alloc1 allocates/deallocates memory for retrievals
    !
    ! !CALLING SEQUENCE: 
    !
    !        call alloc1
    !     
    ! !INPUT PARAMETERS:   
    !
    ! !OUTPUT PARAMETERS:  
    !
    ! !SEE ALSO:  
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner     original fortran 90
    !
    !EOP
    !-------------------------------------------------------------------------
    !
    !character(len=50) :: myname='alloc1: '

    !**************************************************************************

    if (allocated(ind)) then
      !deallocate memory
      !=================
      deallocate(ind)
      deallocate(waves)
      deallocate(wavesp)
      deallocate(w1p)
      deallocate(f1p)
      deallocate(wave_diff)
      deallocate(sflx)
      deallocate(bad_obs)
      deallocate(good_obs)
      deallocate(r_i)
      deallocate(y_obs)
      deallocate(y_obs1)
      deallocate(y_frac)
      deallocate(y_calc_squeeze)
      deallocate(y_calc)
      deallocate(y_calc_sh)
      deallocate(y_resid)
      deallocate(y_back)
      deallocate(jacob_dummy)
      deallocate(jacob_rad)
      deallocate(fit_rad)
      deallocate(rad_cld)
      deallocate(rad_clr)
      deallocate(rad_clr_oc)
      deallocate(ring_cld)
      deallocate(ring_clr)
      deallocate(ring_oc)
      deallocate(comp_all)
      deallocate(comp_all_ring)
      deallocate(computed)
      deallocate(computed_oc)
      deallocate(comp_clr)
      deallocate(comp_oc_clr)
      deallocate(rad_tot)
      deallocate(ycalc)
      deallocate(wavesd)
      deallocate(o3_xsect)
      deallocate(chls)
      deallocate(ntot)
      deallocate(nia_clds)
      deallocate(nra_clds)
      deallocate(tra_clds)
      deallocate(i0a_clds)
      deallocate(z1_clds)
      deallocate(z2_clds)
      deallocate(i01a_clds)
      deallocate(ring_clds)
      deallocate(rad_clds)
      deallocate(rad_clrs)
      deallocate(ring_clrs)
      deallocate(ring_ocs)
      deallocate(reference_ring)
      deallocate(reference_rad)
      if (write_resid) then
        deallocate(wave_resid)
        deallocate(resid)
      endif
    endif ! allocated

    allocate(ind(0:nobs-1))
    allocate(waves(0:nobs-1))
    allocate(wavesp(nterms-1,0:nobs-1))
    allocate(w1p(nwl))
    allocate(f1p(nwl))
    allocate(wave_diff(0:nobs-1))
    allocate(sflx(0:nobs-1))
    allocate(bad_obs(nobs)) 
    allocate(good_obs(nobs)) 
    allocate(r_i(0:nobs-1)) 
    allocate(y_obs(0:nobs-1))   
    allocate(y_obs1(0:nobs-1))   
    allocate(y_frac(0:nobs-1))   
    allocate(y_calc(0:nobs-1))   
    allocate(y_calc_sh(0:nobs-1))   
    allocate(y_calc_squeeze(0:nobs-1))   
    allocate(y_resid(0:nobs-1,1))   
    allocate(y_back(0:nobs-1,1))   
    allocate(jacob_dummy(0:nobs-1))   
    allocate(jacob_rad(0:nobs-1))   
    allocate(fit_rad(0:nobs-1))   
    allocate(rad_cld(0:nobs-1))   
    allocate(rad_clr(0:nobs-1))   
    allocate(rad_clr_oc(0:nobs-1))   
    allocate(ring_cld(0:nobs-1))   
    allocate(ring_clr(0:nobs-1))   
    allocate(ring_oc(0:nobs-1))   
    allocate(computed_oc(nchl))
    allocate(computed(npres))
    allocate(comp_all(npres))
    allocate(comp_all_ring(npres))
    allocate(chls(0:nchl-1)) ; chls=chl
    allocate(rad_tot(0:nobs-1))   
    allocate(ycalc(0:nobs-1))   
    allocate(wavesd(0:nobs-1))
    allocate(o3_xsect(0:nobs-1))
    allocate(ntot(npres,0:nobs-1))
    allocate(nia_clds(npres,0:nobs-1))
    allocate(nra_clds(npres,0:nobs-1))
    allocate(tra_clds(npres,0:nobs-1))
    allocate(i0a_clds(npres,0:nobs-1))
    allocate(z1_clds(npres,0:nobs-1))
    allocate(z2_clds(npres,0:nobs-1))
    allocate(i01a_clds(npres,0:nobs-1))
    allocate(ring_clds(npres,0:nobs-1))
    allocate(rad_clds(npres,0:nobs-1))
    allocate(rad_clrs(npres,0:nobs-1))
    allocate(ring_clrs(npres,0:nobs-1))
    allocate(comp_oc_clr(npres))
    allocate(comp_clr(npres))
    allocate(ring_ocs(nchl,0:nobs-1))
    allocate(reference_ring(nobs,nscanpos)) ; reference_ring=0.
    allocate(reference_rad(nobs,nscanpos)) ; reference_rad=0.
    if (write_resid) then
      allocate(wave_resid(nobs)) ; wave_resid = 0.
      allocate(resid(nobs,nscanpos,nTimes)) ; resid=9999.
    endif

  end subroutine alloc1

end module m_alloc1
