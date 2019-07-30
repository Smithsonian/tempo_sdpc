module slitfunction_super_gaussian
  use OMSAO_precision_module, only: i4, r8
contains
  SUBROUTINE super_gaussian_sf ( npoints, hwem, aw, ak,k, wvl, spec, conv_spec)

    ! =========================================================================
    !
    ! Convolves input spectrum with an assymetric super gaussian slit function of
    ! specified HWEM (half-width at 1/e intensity) and a shape variable (k).
    !
    ! The symmetric super gaussian is defined by:
    !
    !      S(x) = A(w,k) x exp (-(x/w)^k)
    ! where
    !                        k
    !     A(w,k) =   ==================
    !                 2 x w x Gamma(1/k)
    !
    !    w = HWEM (or FWEM = 2w)
    !    k is the shape factor
    !    k=2 makes it a normal gaussian
    !    Gamma function: https://en.wikipedia.org/wiki/Gamma_function
    !    
    !  For having an asymmetric super gaussian, we need to add two parameters:
    !  
    !     S(x) = A(w,k) x exp (- | x/(w - aw) | ^ (k - ak)) for x<=0
    !  	  S(x) = A(w,k)	x exp (- | x/(w + aw) | ^ (k + ak)) for	x>0
    ! 
    !     where aw and ak are asymmetry parameters
    !
    !     For now on, ak should be set zero until we find an appropirate normalized
    !     function. 
    !     The function gets shifted based on aw to make the center of mass = 0
    !     
    !     The asymmetric super gaussion equals to a symmetric one under (ak=aw=0.0_r8)
    ! =========================================================================


    USE sao_pge_utils, ONLY: signdp
    USE integration_routines, ONLY: cubint
    use slatec_davint, only : davint
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: npoints
    REAL    (KIND=r8),                      INTENT (IN) :: hwem, k
    REAL    (KIND=r8), DIMENSION (npoints), INTENT (IN) :: wvl, spec
    REAL    (KIND=r8),                      INTENT (IN) :: aw,ak
    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (npoints), INTENT (OUT) :: conv_spec

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                        :: i, j, nslit, sslit, eslit, davint_err
    REAL    (KIND=r8)                        :: slitsum, sliterr, cwvl, lwvl, rwvl
    REAL    (KIND=r8)                        :: amp_sg ! A(w,k)
    REAL    (KIND=r8), DIMENSION (3*npoints) :: spc_temp, wvl_temp, sf_val, xtmp, ytmp

    sslit = 1; eslit = 1   ! silence compiler warning

    ! --------------------------------------------------------
    ! Initialize output variable (default for "no convolution"
    ! --------------------------------------------------------
    conv_spec(1:npoints) = spec(1:npoints)

    ! -----------------------------------------------
    ! No Super Gaussian convolution if Halfwidth @ 1/e is 0
    ! -----------------------------------------------
    IF ( hwem == 0.0_r8 ) RETURN
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
    spc_temp(npoints+1:2*npoints) = spec(1:npoints)
    wvl_temp(npoints+1:2*npoints) = wvl (1:npoints)
    DO i = 1, npoints
      spc_temp(npoints+1-i) = spec(i)
      wvl_temp(npoints+1-i) = 2.0_r8*wvl(1)-wvl(i) -0.001_r8
      spc_temp(2*npoints+i) = spec(npoints+1-i)
      wvl_temp(2*npoints+i) = 2.0_r8*wvl(npoints)-wvl(npoints+1-i) +0.001_r8
    END DO

    ! ------------------------------------------------------------------------
    ! We now compute the super Gaussian for every point in the spectrum.
    ! Starting from the center point, we go outwards and stop accumulating
    ! points when both sides are less than 0.001 of the maximum slit function.
    ! Since we are starting at the center wavelength, this can be set to 1.0.
    ! Remember that the original wavelength array is now located at indices
    ! NPOINTS+1:2*NPOINTS
    ! ------------------------------------------------------------------------

    DO i = 1, npoints
      amp_sg = k / ( 2.0_r8 * hwem * GAMMA(1.0_r8/k) )
      sf_val = 0.0_r8
      cwvl = wvl_temp(npoints+i)
      sf_val(npoints+i) = amp_sg * EXP(-(ABS((aw)/(hwem+aw)))**(k-ak)) ! value at center
      getslit: DO j = 1, npoints
        sslit = npoints+i-j ; lwvl = wvl_temp(sslit) - cwvl + aw 
        eslit = npoints+i+j ; rwvl = wvl_temp(eslit) - cwvl + aw
        IF ( lwvl .LE. 0) THEN
        sf_val(sslit) = amp_sg * EXP(-(ABS((lwvl)/(hwem-aw)))**(k-ak))
	ELSE
        sf_val(sslit) = amp_sg * EXP(-(ABS((lwvl)/(hwem+aw)))**(k+ak))
	ENDIF
        IF ( rwvl .LE. 0) THEN
        sf_val(eslit) = amp_sg * EXP(-(ABS((rwvl)/(hwem-aw)))**(k-ak))
        ELSE
	sf_val(eslit) = amp_sg * EXP(-(ABS((rwvl)/(hwem+aw)))**(k+ak))
        ENDIF
        IF ( sf_val(sslit) < 0.0005_r8 .AND. sf_val(eslit) < 0.0005_r8 ) EXIT getslit
      END DO getslit

!	OPEN(UNIT=15, FILE="yslit_10.dat", ACTION="write", STATUS="replace")
!	WRITE(15,*) sf_val(sslit:eslit)
!	CLOSE(UNIT=15)
!	OPEN(UNIT=14, FILE="xslit_10.dat", ACTION="write", STATUS="replace")
!	WRITE(14,*) wvl_temp(sslit:eslit)-cwvl
!	CLOSE(UNIT=14)
!	STOP
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
          write (*,*)'*** super_gaussian_sf::davint failed, davint_err=', &
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

      if (nslit > 3) then
        ! jch - it's an error to call cubint with nslit < 4
        CALL cubint ( &
          nslit, xtmp(1:nslit), ytmp(1:nslit), 1, nslit, conv_spec(i), sliterr)
      else
	CALL DAVINT ( &
           xtmp(1:nslit), sf_val(1:nslit), nslit, xtmp(1), xtmp(nslit), &
           conv_spec(i), davint_err )
        if (davint_err /= 1) then
          write (*,*)'*** super_gaussian_sf::davint failed, davint_err=', &
            davint_err
        endif
      endif
    END DO

    RETURN
  END SUBROUTINE super_gaussian_sf
end module slitfunction_super_gaussian


