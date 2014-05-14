module m_interp_ring_rad

  public interp_ring_rad, interp_rad

contains

  subroutine interp_ring_rad(ix1,reflec,computed, rad, ring, drad_tot_dr)

    use m_cloud_pres_mod, ONLY: l,j,nt,i01a_clds, i0a_clds, tra_clds, &
         nia_clds, nra_clds, ind, nobs, ntot, table, temp3D, comp_all, &
         z1_clds, z2_clds, sz, satz, az, comp_all_ring
    use m_vars, ONLY: i01a, i0a, tra, nia, nra, sba, nba, k1bar, z1, z2
    use m_trilin
    use mathcons, ONLY: deg2rad
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  interp_ring_rad
    ! 
    ! !DESCRIPTION: interp_ring_rad allocates/deallocates memory for 
    !               retrievals		
    !
    ! !CALLING SEQUENCE: 
    !
    !        call interp_ring_rad
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
    logical, dimension(:), intent(out) :: computed
    integer, intent(in) :: ix1
    real (KIND=8),    intent(in) :: reflec
    real (KIND=8), dimension(:), intent(inout), optional :: ring
    real (KIND=8), dimension(:), intent(inout) :: rad
    real (KIND=8), dimension(:), intent(inout), optional :: drad_tot_dr
    real (KIND=8) :: refl2, r1, r2
    real (KIND=8), allocatable, dimension(:) :: den, den2, dring_dr, drad_dr

    !**************************************************************************
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
      allocate(den(nobs))
      allocate(den2(nobs))
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
          allocate(dring_dr(nobs))
          allocate(drad_dr(nobs))
          dring_dr=k1bar(ind)/rad*( (nra_clds(ix1,:)/den) + &
               reflec*nra_clds(ix1,:)*sba(ix1,ind)/den2 + &
               2*reflec*nba(ix1,ind)*tra_clds(ix1,:)/den2 + &
               2*refl2*nba(ix1,ind)*tra_clds(ix1,:)*sba(ix1,ind)/(den2*den) )
          drad_dr=tra_clds(ix1,:)/den+reflec*tra_clds(ix1,:)*sba(ix1,ind)/den2
          drad_tot_dr=drad_dr*(1+ring) + rad*dring_dr
          deallocate(dring_dr)
          deallocate(drad_dr)
        endif
      endif ! if not computed
    endif ! if present(ring)
    computed(ix1) = .true.
    comp_all(ix1) = .true.
    if (allocated(den)) deallocate(den)
    if (allocated(den2)) deallocate(den2)


  end subroutine interp_ring_rad

  subroutine interp_rad (ix1,ind, i0, sb, tr)

    use m_trilin
    use m_triquad
    use m_vars, ONLY: i0a, tra, sba, z1, z2
    use m_cloud_pres_mod, ONLY: l,j,nt, table, temp3D, sz, satz, az
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
