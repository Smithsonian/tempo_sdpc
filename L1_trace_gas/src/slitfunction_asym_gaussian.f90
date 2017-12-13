module slitfunction_asym_gaussian
  use OMSAO_precision_module, only: i4, r8
contains
  SUBROUTINE asymmetric_gaussian_sf ( npoints, hw1e, e_asym, wvlarr, specarr, specmod)

    ! =========================================================================
    !
    ! Convolves input spectrum with an asymmetric Gaussian slit function of
    ! specified HW1E (half-width at 1/e intensity) and asymmetry factor E_ASYM.
    !
    ! The asymetric Gaussian g(x) is defined as
    !                   _                                   _
    !                  |               x^2                   |
    !      g(x) =  EXP | - --------------------------------- |
    !                  |_   (hw1e * (1 + SIGN(x)*e_asym))^2 _|
    !
    ! g(x) becomes symmetric for E_ASYM = 0.
    !
    ! =========================================================================

    USE sao_pge_utils, ONLY: signdp
    USE integration_routines, ONLY: cubint
    use slatec_davint, only : davint
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: npoints
    REAL    (KIND=r8),                      INTENT (IN) :: hw1e, e_asym
    REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (npoints), INTENT (OUT) :: specmod

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                        :: i, j, nslit, sslit, eslit, davint_err
    REAL    (KIND=r8)                        :: slitsum, sliterr, cwvl, lwvl, rwvl
    REAL    (KIND=r8), DIMENSION (3*npoints) :: spc_temp, wvl_temp, sf_val, xtmp, ytmp

    !REAL (KIND=r8) :: signdp
    !EXTERNAL signdp

    sslit = 1; eslit = 1   ! silence compiler warning

    ! --------------------------------------------------------
    ! Initialize output variable (default for "no convolution"
    ! --------------------------------------------------------
    specmod(1:npoints) = specarr(1:npoints)

    ! -----------------------------------------------
    ! No Gaussian convolution if Halfwidth @ 1/e is 0
    ! -----------------------------------------------
    IF ( hw1e == 0.0_r8 ) RETURN

    ! ------------------------------------------------------------------------
    ! One temporary variable is SPC_TEMP, which is three times the size of
    ! SPEC. For the convolution routine to work (hopefully) in each and
    ! every case, we reflect the spectrum at its end points to always have a
    ! fully filled slit function. But this causes some real index headaches
    ! when the slit function wraps around at the ends. Performing the mirror
    ! imaging before we get to the convolution helps to keep things a little
    ! more simple.
    !
    ! Note that this approach is the same as for the pre-tabulated OMI lab
    ! slit function. It is adopted here because now we not only convolve the
    ! solar spectrum, but also any higher resolution reference cross sections,
    ! and these may not necessarily be equidistant in wavelength.
    ! ------------------------------------------------------------------------
    spc_temp(npoints+1:2*npoints) = specarr(1:npoints)
    wvl_temp(npoints+1:2*npoints) = wvlarr (1:npoints)
    DO i = 1, npoints
      spc_temp(npoints+1-i) = specarr(i)
      wvl_temp(npoints+1-i) = 2.0_r8*wvlarr(1)-wvlarr(i) -0.001_r8
      spc_temp(2*npoints+i) = specarr(npoints+1-i)
      wvl_temp(2*npoints+i) = 2.0_r8*wvlarr(npoints)-wvlarr(npoints+1-i) +0.001_r8
    END DO

    ! ------------------------------------------------------------------------
    ! We now compute the asymmetric Gaussian for every point in the spectrum.
    ! Starting from the center point, we go outwards and stop accumulating
    ! points when both sides are less than 0.001 of the maximum slit function.
    ! Since we are starting at the center wavelength, this can be set to 1.0.
    ! Remember that the original wavelength array is now located at indices
    ! NPOINTS+1:2*NPOINTS
    ! ------------------------------------------------------------------------
    DO i = 1, npoints
      sf_val = 0.0_r8
      cwvl = wvl_temp(npoints+i)

      sf_val(npoints+i) = 1.0_r8
      getslit: DO j = 1, npoints
        sslit = npoints+i-j ; lwvl = wvl_temp(sslit) - cwvl
        eslit = npoints+i+j ; rwvl = wvl_temp(eslit) - cwvl
        sf_val(sslit) = EXP(-lwvl**2 / ( hw1e * (1.0_r8 + signdp(lwvl)*e_asym) )**2)
        sf_val(eslit) = EXP(-rwvl**2 / ( hw1e * (1.0_r8 + signdp(rwvl)*e_asym) )**2)
        IF ( sf_val(sslit) < 0.0005_r8 .AND. sf_val(sslit) < 0.0005_r8 ) EXIT getslit
      END DO getslit

      ! ----------------------------------
      ! The number of slit function points
      ! ----------------------------------
      nslit = eslit - sslit + 1
      ! ----------------------------------------------------------------
      ! Compute the norm of the slitfunction. It should be close to 1
      ! already, but making sure doesn't hurt.
      ! ----------------------------------------------------------------
      xtmp(1:nslit) = wvl_temp(sslit:eslit)-cwvl
      ytmp(1:nslit) = sf_val  (sslit:eslit)
!      CALL cubint ( &
!        nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, slitsum, sliterr)
!      !!CALL DAVINT ( &
!      !!     xtmp(1:nslit), sf_val(1:nslit), nslit, xtmp(1), xtmp(nslit), &
!      !!     slitsum, locerrstat )
      if (nslit > 3) then
        ! jch - it's an error to call cubint with nslit < 4
        CALL cubint ( &
          nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, slitsum, sliterr)
      else
        CALL DAVINT ( &
          xtmp(1:nslit), sf_val(1:nslit), nslit, xtmp(1), xtmp(nslit), &
          slitsum, davint_err )
        if (davint_err /= 1) then
          write (*,*)'*** asymmetric_gaussian_sf::davint failed, davint_err=', &
            davint_err
        endif
      endif

      IF ( slitsum > 0.0_r8 ) sf_val(sslit:eslit) = sf_val(sslit:eslit) / slitsum

      ! ---------------------------------------------------------------------
      ! Prepare array for integration: Multiply slit function values with the
      ! spectrum array to be convolved.
      ! ---------------------------------------------------------------------
      ytmp(1:nslit) = sf_val(sslit:eslit) * spc_temp(sslit:eslit)

      ! ----------------------------------------------------------
      ! Folding (a.k.a. integration) of spectrum and slit function
      ! ----------------------------------------------------------
!      !!CALL DAVINT ( &
!      !!     xtmp(1:nslit), sf_val(1:nslit), nslit, xtmp(1), xtmp(nslit), &
!      !!     specmod(i), locerrstat )
!      CALL cubint ( &
!        nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, specmod(i), sliterr)
      if (nslit > 3) then
        ! jch - it's an error to call cubint with nslit < 4
        CALL cubint ( &
          nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, specmod(i), sliterr)
      else
        CALL DAVINT ( &
           xtmp(1:nslit), sf_val(1:nslit), nslit, xtmp(1), xtmp(nslit), &
           specmod(i), davint_err )
        if (davint_err /= 1) then
          write (*,*)'*** asymmetric_gaussian_sf::davint failed, davint_err=', &
            davint_err
        endif
      endif
    END DO

    RETURN
  END SUBROUTINE asymmetric_gaussian_sf  
end module slitfunction_asym_gaussian
