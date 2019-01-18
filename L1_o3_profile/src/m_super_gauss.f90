! =========================================================================
! Author : Juseon Bak
! Date : Sep. 2017
! it is modified from gauss.f90
! Purpose :
! Convolves input spectrum with an super Gaussian slit function 
! where the exponent of the Gaussian is a variable k and will give a Gaussian for k=2. 
! For big N the function will describe a more rectangular distribution, 
! while for small k it fits to a distribution with long tails on both sides. 
!
! The symetric super Gaussian g(x) is defined as
!                               _               _
!              K               |      abs(x)^k   |
!      g(x) =-----         EXP | - ------------- |
!            2hw1e gam(1/k)    |_     hw1e ^k   _|
! Reminder
!    slit function derivatives are implemented in super_gauss_multi, super_gauss_fc
! =========================================================================

MODULE m_super_gauss

 USE OMSAO_precision_module
 USE OMSAO_indices_module,   ONLY : hwe_idx, spk_idx
 USE OMSAO_variables_module, ONLY : slit_trunc_limit,  winlim, numwin,&
     solwinfit, slitwav, slitfit, nslit, & ! slit function variables
     do_dsdw, do_dsdk ! logical variables to implement sf derivaties
 USE m_ezspline_interpolation, ONLY: interpolation 
 USE OMSAO_errstat_module

 IMPLICIT NONE

 PUBLIC :: super_gauss, super_gauss_multi, super_gauss_vary, &
           super_gauss_f2c, super_gauss_vary_f2c
 PRIVATE

CONTAINS
SUBROUTINE super_gauss (wvlarr, specarr, specmod, npoints, hw1e, power)

  IMPLICIT NONE
  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN)    :: npoints
  REAL (KIND=dp),                      INTENT (IN)    :: hw1e, power
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT)   :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: nhi, nlo, i, j, num_slit
  REAL (KIND=dp)                      :: emult, delwvl, slitsum, slit0, coef
  REAL (KIND=dp), DIMENSION (npoints) :: slit

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  ! --------------------------------------
  ! No convolution if halfwidth @ 1/e is 0
  ! --------------------------------------
  IF ( hw1e == 0.0 .or. power ==0.0) THEN 
     WRITE(*,*) 'super gaussian error'
     STOP
  ENDIF
  
  delwvl = wvlarr(2) - wvlarr(1)
  !  Calculate slit function values out to 0.001 times x0 value,
  !     normalize so that sum = 1.
  coef   = power / ( 2*hw1e  *gamma(power)) 
  coef   = 1
  emult  = -1.0/(hw1e)**power
  slitsum = coef ! this initial value equal to center of slit
  slit0   = coef
  i = 1  ;  num_slit = 0
  DO WHILE ( num_slit <= npoints )
     slit (i) = coef*EXP (emult*(abs(delwvl*i))**power)
     slitsum = slitsum + 2.0 * slit (i)
     IF (slit (i) <= slit_trunc_limit ) EXIT 
     i = i + 1
  ENDDO
  num_slit = i

  slit0 = slit0 / slitsum
  slit(1:num_slit) = slit(1:num_slit) / slitsum
 
  ! Convolve spectrum. reflect at endpoints.
  ! Doesn't look right
  specmod(1:npoints) = slit0 * specarr(1:npoints)
  DO i = 1, npoints
     DO j = 1, num_slit
        nlo = i - j 
        IF (nlo < 1) nlo = -nlo + 2 
        nhi = i + j
        IF ( nhi > npoints ) nhi = npoints - MOD(nhi, npoints)
        specmod(i) = specmod(i) + slit(j) * ( specarr(nlo) + specarr(nhi) )
     END DO
  END DO
  RETURN
END SUBROUTINE super_gauss


SUBROUTINE super_gauss_multi (wvlarr, specarr, specmod, npoints)


  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN) :: npoints
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: i, j, ii, j1, j2, num_slit, mslit, fidx, lidx, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum, hw1e, power, coef, pert, ssum1
  REAL (KIND=dp), DIMENSION (npoints) :: slit, locwvl, upbnd, slit1

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  fidx = 1
  DO iwin = 1, numwin

     hw1e  = solwinfit(iwin, hwe_idx, 1)
     power = solwinfit(iwin, spk_idx, 1)
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF

     ! -----------------------------------------------
     ! No Gaussian convolution if Halfwidth @ 1/e is 0
     ! -----------------------------------------------
     IF ( hw1e == 0.0 .or. power == 0.0 ) STOP

     ! --------------------------------------------------------------
     ! Find the number of spectral points that fall within a Gaussian
     ! slit function with values >= 0.001. Remember that we have an
     ! asymmetric Gaussian, so we create a wavelength array symmetric
     ! around 0. The spacing is provided by the equidistant WVLARR.
     ! --------------------------------------------------------------
    ! sw = - hw1e ** 2.0
     IF (fidx > npoints .or. lidx < 1 ) return
     delwvl = wvlarr(fidx+1) - wvlarr(fidx)
    ! mslit = NINT( SQRT( LOG(slit_trunc_limit) * sw) / delwvl) 
     mslit = NINT((hw1e*((-LOG(slit_trunc_limit))**(1/power)))/delwvl)
     num_slit = mslit * 2 + 1
     
     IF (num_slit > npoints) THEN
        mslit = (npoints-1) / 2 ;    num_slit = mslit * 2 + 1
     ENDIF
         
     ! only compute half of the slit shape (say left side)
     j = -mslit
     DO i = 1, mslit
        locwvl(i) = delwvl * j
        j = j + 1
     ENDDO
     coef   = power / ( 2*hw1e  *gamma(power))
     coef   = 1.0
     slit (1:mslit) = coef*EXP (-1.*(abs(locwvl(1:mslit))/hw1e)**power)
     ! Normalization
     slitsum = SUM(slit(1:mslit)) * 2.0 + coef   ! (center value is 1.0)
     slit(1:mslit) = slit(1:mslit) / slitsum

      IF (do_dsdw) THEN
        slit1(1:mslit) =coef* EXP (-1.*(abs(locwvl(1:mslit))/(hw1e*1.001))**power)
        ssum1 = sum(slit1(1:mslit))*2.0 + coef
        pert = hw1e*0.001
        slit (1:mslit) = (slit1(1:mslit)/ssum1 - slit(1:mslit))/pert
      ELSE IF (do_dsdk) THEN
        slit1(1:mslit) =coef* EXP(-1.*(abs(locwvl(1:mslit))/hw1e)**(power*1.001))
        ssum1 = sum(slit1(1:mslit))*2.0 + coef
        pert = power*0.001
        slit (1:mslit) = (slit1(1:mslit)/ssum1 - slit(1:mslit))/pert
      ENDIF
          
     ! Make a local copy of the NSLIT spectrum points to be convolved
     ! with the slit function. The spectrum points to be convolved are
     ! arranged such that the updated index corresponds to the maximum
     ! of the slit function (MSLIT). For simplicity we reflect the
     ! spectrum at the array end points.
      
     ! ----------------------------------------------------
     ! Loop over all points of the spectrum to be convolved
     ! ----------------------------------------------------
     DO i = fidx, lidx
        ! Center value
        specmod(i) = specarr(i)* coef/ slitsum  ! (i.e., 1.0 / slitsum for center)

        ! Add both left and side (symmetric)
        DO j = 1, mslit
           ii = mslit + 1 - j
           j1 = i + ii; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
           j2 = i - ii; IF ( j2 < 1 ) j2 = ABS(j2) + 2 
           specmod(i) = specmod(i) + slit(j) * ( specarr(j1) + specarr(j2))
        ENDDO
     END DO

     fidx = lidx + 1
  END DO

  RETURN
END SUBROUTINE super_gauss_multi

SUBROUTINE super_gauss_vary (wvlarr, specarr, specmod, npoints) ! need to check in detail

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                             INTENT (IN) :: npoints
  REAL (KIND=dp), DIMENSION (npoints), INTENT (IN) :: wvlarr, specarr

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=dp), DIMENSION (npoints), INTENT (OUT) :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: i, ii, j, j1, j2, num_slit, mslit, &
       errstat,  fidx, lidx, fslit, lslit, finter, linter, iwin
  REAL (KIND=dp)                      :: delwvl, slitsum,  hw1e, sw, upbnd, POWER, coeff
  REAL (KIND=dp), DIMENSION (npoints) :: slit, lochwe, locspk

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=*), PARAMETER :: modulename = 'super_gauss_vary'


  !WRITE(*, *) 'nslit = ', nslit
  !WRITE(*, '(10f8.3)') slitwav(1:nslit), slitfit(1:nslit, hwe_idx, 1)

  ! --------------------------------------------------------
  ! Initialize output variable (default for "no convolution"
  ! --------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  

  ! -----------------------------------------------
  ! No Gaussian convolution if Halfwidth @ 1/e is 0
  ! -----------------------------------------------
  hw1e   = MAXVAL(slitfit(1:nslit, hwe_idx, 1))
  power  = MAXVAL(slitfit(1:nslit,spk_idx, 1))

  IF ( hw1e == 0.0 .or. power == 0.0 ) RETURN

  ! --------------------------------------------------------------
  ! Find the number of spectral points that fall within a Gaussian
  ! slit function with values >= 0.01. Remember that we have an
  ! asymmetric Gaussian, so we create a wavelength array symmetric
  ! around 0. The spacing is provided by the equidistant WVLARR.
  ! --------------------------------------------------------------  
  sw = -hw1e ** 2.0
  delwvl = wvlarr(2) - wvlarr(1)
  ! mslit = NINT( SQRT( LOG(slit_trunc_limit) * sw) / delwvl) 
  mslit = NINT((hw1e*((-LOG(slit_trunc_limit))**(1/power)))/delwvl)
  num_slit = mslit * 2 + 1

  IF (num_slit > npoints) THEN
     mslit = (npoints-1) / 2 + 1
  ENDIF
 
  ! ---------------------------------------------------------------
  ! Convolve spectrum. First do the middle part, where we have full
  ! overlap coverage of the slit function. Again, remember the
  ! asymmetry of the Gaussian, which makes impossible a simple
  ! 50-50 division of the summation interval.
  ! ---------------------------------------------------------------

  ! Make a local copy of the NSLIT spectrum points to be convolved
  ! with the slit function. The spectrum points to be convolved are
  ! arranged such that the updated index corresponds to the maximum
  ! of the slit function (MSLIT). For simplicity we reflect the
  ! spectrum at the array end points.

  errstat = pge_errstat_ok

  ! get an array of slit variables lochwe for each wavelength position
  lochwe = 0.0

  fidx = 1; lidx = npoints
  DO iwin = 1, numwin
     IF (iwin < numwin) THEN
        upbnd = (winlim(iwin, 2) + winlim(iwin+1, 1))/2.0 
        lidx = MINVAL(MAXLOC(wvlarr, MASK=(wvlarr <= upbnd )))
     ELSE 
        lidx = npoints
     END IF
               
     fslit = MINVAL(MINLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) >= wvlarr(fidx))))
     lslit = MINVAL(MAXLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) <= wvlarr(lidx))))
 
     ! no slit between fidx:lidx, should never happen
     IF (fslit <= 0) THEN
        fslit = lslit
     ELSE IF (lslit <= 0) THEN
        lslit = fslit
     ENDIF

     IF (fslit > nslit .OR. lslit > nslit) THEN  
        WRITE(*, *) fslit, lslit, wvlarr(fidx), wvlarr(lidx), slitwav(1), &
             slitwav(nslit), nslit, wvlarr(1), wvlarr(npoints)
        WRITE(*, *) modulename, ': Not slit available for this window!!!'
        STOP
     ENDIF
     
     IF (lslit < fslit + 3) THEN  ! extrapolate, use the nearest value
        lochwe(fidx:lidx) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
        locspk(fidx:lidx) = SUM(slitfit(fslit:lslit, spk_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(wvlarr, MASK = (wvlarr >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(wvlarr, MASK = (wvlarr <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, wvlarr(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, spk_idx, 1), &
             linter-finter+1, wvlarr(finter:linter), locspk(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        IF (finter > fidx) THEN
           lochwe(fidx:finter-1)=slitfit(fslit, hwe_idx, 1)
           locspk(fidx:finter-1)=slitfit(fslit, spk_idx, 1)
        END IF
        
        IF (linter < lidx)  THEN
           lochwe(linter+1:lidx)=slitfit(lslit, hwe_idx, 1)
           locspk(linter+1:lidx)=slitfit(lslit, spk_idx, 1)
        END IF


     ENDIF
     fidx = lidx + 1
  ENDDO
    
  DO i = 1, npoints
     ! Get center value
     coeff = locspk(i) / ( 2 *lochwe(i) *gamma(locspk(i)))
     coeff = 1.0
     specmod(i) = specarr(i)
     slitsum = coeff

     ! Compute slit contribution from both sides
     DO j = 1, mslit
        ii = mslit + 1 - j
        j1 = i + ii; IF ( j1 > npoints ) j1 = npoints - MOD(j1, npoints)
        j2 = i - ii; IF ( j2 < 1 ) j2 = ABS(j2) + 2 
        slit(j) = coeff*EXP (-1.*(abs(wvlarr(j1)-wvlarr(i))/lochwe(i))**locspk(i))
        slitsum = slitsum + 2.0 * slit(j)
        specmod(i) = specmod(i) + slit(j) * ( specarr(j1) + specarr(j2))
     ENDDO
     
     specmod(i) = specmod(i)*coeff / slitsum !normalized to center value
  END DO
  !stop ! not tested
  RETURN
END SUBROUTINE super_gauss_vary


SUBROUTINE super_gauss_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
  
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, iwin, fidx, fidxc, lidx, lidxc, midx, sidx, eidx, nhalf
  REAL (KIND=dp) :: power, temp, hw1e, dfw, ssum, ssum1, coeff, pert
  REAL (KIND=dp), DIMENSION (nf) :: slit, slit1
  
  !slit(:) = 0.0; cspec(:, :)= 0.0
  fidx = 1; fidxc = 1
  dfw  = fwave(2) - fwave(1) 
  DO iwin = 1, numwin

     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF    
     power  = solwinfit(iwin, spk_idx, 1)
     hw1e   = solwinfit(iwin, hwe_idx, 1)
     nhalf  = CEILING(hw1e/dfw*(-LOG(slit_trunc_limit))**(1/power))
     coeff  = power / ( 2 * hw1e *gamma(power))
     coeff = 1.0
     DO i = fidxc, lidxc
        ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx
        sidx = MAX(midx - nhalf, 1)
        eidx = MIN(nf, midx + nhalf)
        slit (sidx:eidx) =coeff* EXP (-1.*(abs(cwave(i) - fwave(sidx:eidx))/hw1e)**power)         
        ssum = SUM(slit(sidx:eidx))
        slit (sidx:eidx) = slit(sidx:eidx)/ssum
        IF (do_dsdw) THEN 
          slit1(sidx:eidx) =coeff* EXP (-1.*(abs(cwave(i) - fwave(sidx:eidx))/(hw1e*1.001))**power)         
          ssum1 = sum(slit1(sidx:eidx))
          pert = hw1e*0.001
          slit (sidx:eidx) = (slit1(sidx:eidx)/ssum1 - slit(sidx:eidx))/pert
        ELSE IF (do_dsdk) THEN 
          slit1 (sidx:eidx) =coeff* EXP (-1.*(abs(cwave(i) - fwave(sidx:eidx))/hw1e)**(power*1.001))         
          ssum1 = sum(slit1(sidx:eidx))
          pert = power*0.001
          slit (sidx:eidx) =(slit1(sidx:eidx)/ssum1 - slit(sidx:eidx))/pert
        ENDIF  
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) !/pert
        ENDDO
     ENDDO
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO 
  RETURN

END SUBROUTINE super_gauss_f2c


SUBROUTINE super_gauss_vary_f2c (fwave, fspec, nf, nspec, cwave, cspec, nc)
 
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER,                        INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=dp), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=dp), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=dp), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=dp), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER        :: i, j, iwin, fidx, fidxc, lidx, lidxc, &
       midx, sidx, eidx, fslit,lslit, finter, linter, errstat, nhalf
  REAL (KIND=dp) :: temp, hw1e, power, dfw, ssum, coeff
  REAL (KIND=dp), DIMENSION (nf) :: slit
  REAL (KIND=dp), DIMENSION (nc) :: lochwe, locspk

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER :: OMI_SMF_setmsg

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=*), PARAMETER :: modulename = 'super_gauss_vary_f2c'

  errstat = pge_errstat_ok 
  dfw  = fwave(2) - fwave(1) 
  fidx = 1; fidxc = 1

  DO iwin = 1, numwin
     IF (iwin == numwin) THEN
        lidx = nf; lidxc = nc
     ELSE
        temp = (winlim(iwin, 2) + winlim(iwin + 1, 1)) / 2.0
        lidx =  MINVAL(MAXLOC(fwave, MASK=(fwave <= temp)))
        lidxc = MINVAL(MAXLOC(cwave, MASK=(cwave <= temp)))
     ENDIF 
               
     fslit = MINVAL(MINLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) >= cwave(fidxc))))
     lslit = MINVAL(MAXLOC(slitwav(1:nslit), &
          MASK=(slitwav(1:nslit) <= cwave(lidxc))))

     ! no slit between fidx:lidx, should never happen
     IF (fslit <= 0) THEN
        fslit = lslit
     ELSE IF (lslit <= 0) THEN
        lslit = fslit
     ENDIF

     IF (fslit > nslit .OR. lslit > nslit) THEN  
        WRITE(*, *) fslit, lslit, cwave(fidxc), cwave(lidxc), slitwav(1), &
                    slitwav(nslit), nslit, cwave(1), cwave(nc)
        WRITE(*, *) modulename, ': Not slit available for this window!!!'
        STOP
     ENDIF

     IF (lslit < fslit + 3) THEN  ! extrapolate, use the nearest value
        lochwe(fidxc:lidxc) = SUM(slitfit(fslit:lslit, hwe_idx, 1))/(lslit-fslit+1)
        locspk(fidxc:lidxc) = SUM(slitfit(fslit:lslit, spk_idx, 1))/(lslit-fslit+1)
     ELSE
        ! finter <= linter here
        finter = MINVAL(MINLOC(cwave, MASK = (cwave >= slitwav(fslit))))
        linter = MINVAL(MAXLOC(cwave, MASK = (cwave <= slitwav(lslit))))
        
        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, hwe_idx, 1), &
             linter-finter+1, cwave(finter:linter), lochwe(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF

        CALL interpolation ( &
             lslit-fslit+1, slitwav(fslit:lslit), slitfit(fslit:lslit, spk_idx, 1), &
             linter-finter+1, cwave(finter:linter), locspk(finter:linter), errstat )
        IF ( errstat > pge_errstat_warning ) THEN
           errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0) ; STOP 1
        END IF
        
        IF (finter > fidxc) THEN
           lochwe(fidxc:finter-1)=slitfit(fslit, hwe_idx, 1)
           locspk(fidxc:finter-1)=slitfit(fslit, spk_idx, 1)
        END IF
        
        IF (linter < lidxc)  THEN
           lochwe(linter+1:lidxc)=slitfit(lslit, hwe_idx, 1)
           locspk(linter+1:lidxc)=slitfit(lslit, spk_idx, 1)
        END IF
     ENDIF
   
     DO i = fidxc, lidxc

        hw1e   = lochwe(i)
        power  = locspk(i)
        coeff  = power / ( 2 * hw1e *gamma(power))
        coeff  = 1.0
        nhalf  = CEILING(hw1e / dfw * (-LOG(slit_trunc_limit))**(1/power)) !xliu, 10/22/2009

       ! Find the closest pixel
        midx = MINVAL(MAXLOC(fwave(fidx:lidx), MASK=(fwave(fidx:lidx) <= cwave(i)))) + fidx

        !xliu, 10/22/2009, replace above with following
        sidx = MAX(midx - nhalf, 1)
        eidx = MIN(nf, midx + nhalf)
        slit(sidx:eidx) = coeff*EXP (-1.*(abs(cwave(i) - fwave(sidx:eidx))/hw1e)**power)         
     
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ENDDO
     
     fidx = lidx + 1; fidxc = lidxc + 1
  ENDDO

  RETURN
END SUBROUTINE super_gauss_vary_f2c

END MODULE m_super_gauss
