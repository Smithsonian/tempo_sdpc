module m_interp_ring_rad

  public interp_ring_rad, interp_rad

contains

  subroutine interp_ring_rad(ix1,reflec,computed, rad, ring, drad_tot_dr, &
       errstat)

    use m_cloud_pres_mod, ONLY: l,j,nt,i01a_clds, i0a_clds, tra_clds, &
         nia_clds, nra_clds, ind, nobs, ntot, table, temp3D, comp_all, &
         z1_clds, z2_clds, sz, satz, az, comp_all_ring
    use m_vars, ONLY: i01a, i0a, tra, nia, nra, sba, nba, k1bar, z1, z2
    use m_trilin
    use mathcons, ONLY: deg2rad
    use tell_module
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  interp_ring_rad
    ! 
    ! !DESCRIPTION: interp_ring_rad computes radiance and atmospheric ring
    !   for current pixel and pressure estimate
    !
    ! !INPUT PARAMETERS:   
    !   ix1: cloud pressure index
    !   reflec: reflectivity 
    !
    ! !OUTPUT PARAMETERS:  
    !   computed: vector tracking whether ring and rad have been computed
    !   rad: total reflected radiance at top of atmosphere, I think...
    !   ring: ring effect factor
    !   drad_tot_dr: change in reflected radiance with reflectance.
    !
    ! !SEE ALSO:  
    !   Joiner et al 1995, Appl. Optics, 34, 4513
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner     original fortran 90
    !  12Aug14  O'Sullivan  added documentation, updated for TEMPO
    !
    !EOP
    !-------------------------------------------------------------------------
    !
    logical, dimension(:), intent(out) :: computed
    integer, intent(in) :: ix1
    integer, intent(inout) :: errstat

    real (KIND=8),    intent(in) :: reflec
    real (KIND=8), dimension(:), intent(inout), optional :: ring
    real (KIND=8), dimension(:), intent(inout) :: rad
    real (KIND=8), dimension(:), intent(inout), optional :: drad_tot_dr
    real (KIND=8) :: refl2, r1, r2
    real (KIND=8), allocatable, dimension(:) :: den, den2, dring_dr, drad_dr

    !**************************************************************************

    if (errstat /= 0) return

    if (.not. comp_all(ix1)) then
      temp3D => i0a (ix1,:,:,ind(0):ind(nobs-1))
      i0a_clds(ix1,:) = bilin(j,nt)
      temp3D => z1 (ix1,:,:,ind(0):ind(nobs-1))
      z1_clds(ix1,:) = bilin(j,nt)
      temp3D => z2 (ix1,:,:,ind(0):ind(nobs-1))
      z2_clds(ix1,:) = bilin(j,nt)
      temp3D => tra (ix1,:,:,ind(0):ind(nobs-1))
      tra_clds(ix1,:) = bilin(j,nt)
    endif
    if (.not. computed(ix1)) then
      allocate(den(nobs), den2(nobs), stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "interp_ring_rad: allocation failure", &
             errstat)
        return
      endif

      den=(1-reflec*sba(ix1,ind))
      den2=den**2
      r1=-3./8.*cos(sz*deg2rad)*sin(sz*deg2rad)*sin(satz*deg2rad)*cos(az*deg2rad)
      r2=3./32.*sin(sz*deg2rad)**2*sin(satz*deg2rad)**2/cos(satz*deg2rad)*cos(2.0*az*deg2rad)
      rad = exp(i0a_clds(ix1,:)) + z1_clds(ix1,:)*r1 + &
           z2_clds(ix1,:)*r2 + &
           reflec*tra_clds(ix1,:) / den

    endif

    if (present(ring)) then
      if (.not. comp_all_ring(ix1)) then
        table => i01a(ix1,:,:,:,ind(0):ind(nobs-1))
        i01a_clds(ix1,:) = trilinear(l,j,nt)
        table => nia (ix1,:,:,:,ind(0):ind(nobs-1))
        nia_clds(ix1,:) = trilinear(l,j,nt)
        temp3D => nra (ix1,:,:,ind(0):ind(nobs-1))
        nra_clds(ix1,:) = bilin(j,nt)
        comp_all_ring(ix1) = .true.
      endif ! if not comp_all_ring
      if (.not. computed(ix1)) then
        refl2=reflec**2
        ntot(ix1,:) = (nia_clds(ix1,:)+(reflec*nra_clds(ix1,:)/den) + &
             (refl2*nba(ix1,ind)*tra_clds(ix1,:)/den2) ) &
             / rad
        ring=i01a_clds(ix1,:)/rad+k1bar(ind)*ntot(ix1,:)
        if (present(drad_tot_dr)) then 
          if (allocated(dring_dr)) deallocate(dring_dr, drad_dr, stat=errstat)
          allocate(dring_dr(nobs), drad_dr(nobs), stat=errstat)
          if (errstat /= 0) then
            call tell_error (tell_malloc_error, &
                 "interp_ring_rad: allocation failure", &
                 errstat)
            return
          endif

          ! see Joiner et al (1995) eqn 31
          dring_dr=k1bar(ind)/rad*( (nra_clds(ix1,:)/den) + &
               reflec*nra_clds(ix1,:)*sba(ix1,ind)/den2 + &
               2*reflec*nba(ix1,ind)*tra_clds(ix1,:)/den2 + &
               2*refl2*nba(ix1,ind)*tra_clds(ix1,:)*sba(ix1,ind)/(den2*den) )
          drad_dr=tra_clds(ix1,:)/den+reflec*tra_clds(ix1,:)*sba(ix1,ind)/den2
          drad_tot_dr=drad_dr*(1+ring) + rad*dring_dr
        endif
      endif ! if not computed
    endif ! if present(ring)
    computed(ix1) = .true.
    comp_all(ix1) = .true.
    if (allocated(den)) deallocate(den, den2, stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "interp_ring_rad: deallocation failure", &
           errstat)
      return
    endif


  end subroutine interp_ring_rad

  subroutine interp_rad (ix1,ind, i0, sb, tr)

    !!------------------------------------------------------------------
    !
    ! interpolates various radiance-related parameters to their 
    ! expected values for the satelite and solar zenith angle.
    !
    ! Some guesswork involved in meaning of parameters.
    ! See also m_read_tables.f90
    ! INPUT PARAMETERS
    !  ix1: cloud pressure index
    !  ind: wavelength index 
    !  j: solar zenith angle index (theta)
    !  nt: satellite zenith angle index (scan)
    !  sz: solar zenith angle
    !  satz: satellite zenith angle
    !  az: satellite azimuth angle
    !  i0a: backscattered intensity at top of atmosphere 
    !  tra: transmittance factor ?
    !  sba: fraction of intensity reflected by surface that is then
    !       scattered back to surface by atmosphere?
    !  z1: prop. to intensity from single scatterings?
    !  z2: prop. to intensity from double scatterings?
    ! OUTPUT PARAMETERS
    !  i0: backscattered intensity in pixel
    !  tr: transmittance factor in pixel?
    !  sb: surface light lost to scattering in pixel?
    !
    !!------------------------------------------------------------------

    use m_trilin
    use m_triquad
    use m_vars, ONLY: i0a, tra, sba, z1, z2
!    use m_cloud_pres_mod, ONLY: l,j,nt, table, temp3D, sz, satz, az
    use m_cloud_pres_mod, ONLY: j,nt, temp3D, sz, satz, az
    use mathcons, ONLY: deg2rad
    implicit none

    integer, intent(in) :: ix1
    integer, intent(in) :: ind
    real (KIND=8),    intent(out) :: i0, sb, tr

    real (KIND=8),    dimension(1) :: i01, tr1, z11, z21
    real (KIND=8)                  :: R1, R2

    sb=sba(ix1,ind)
    temp3D => i0a (ix1,:,:,ind:ind)
    i01 = biquad(j,nt)
    temp3D => z1 (ix1,:,:,ind:ind)
    z11 = biquad(j,nt)
    temp3D => z2 (ix1,:,:,ind:ind)
    z21 = biquad(j,nt)
    temp3D => tra (ix1,:,:,ind:ind)
    tr1 = biquad(j,nt)
    r1=-3./8.*cos(sz*deg2rad)*sin(sz*deg2rad)*sin(satz*deg2rad)
    r2=3./32.*sin(sz*deg2rad)**2*sin(satz*deg2rad)**2/cos(satz*deg2rad)
    i0=exp(i01(1))+z11(1)*cos(az*deg2rad)*R1+z21(1)*cos(2.*az*deg2rad)*R2
    tr=tr1(1)

  end subroutine interp_rad

end module m_interp_ring_rad
