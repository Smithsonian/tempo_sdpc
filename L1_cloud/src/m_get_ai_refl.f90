module m_get_ai_refl

contains

  subroutine get_ai_refl(refl_clr, refl_cld, i_obs_l, i_obs_s, i, &
       i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
       sb_ls, tr_ls, set_cld_frac, i0_ss, sb_ss, tr_ss)

    use m_vars, ONLY: &
         refl, iprt, refl_l, dIdR, min_refl, max_refl, min_refl_flag, &
         qc, eff_cld_frac, do_short_wave, iLine, &
         cal_reflec, rad_cld_frac 

    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  get_ai_refl
    ! 
    ! !DESCRIPTION: get_ai_refl calculates radiative cloud fraction, 
    !   cloud reflectivity, aerosol index, etc., using the 
    !   Mixed Lambert-Equivalent Reflectivity (MLER) concept and
    !   Independent Pixel Approximation (IPA). Assumes that a partly 
    !   cloudy pixel is the sum of a mix of cloudy and clear sub-pixels.
    !
    ! !INPUT PARAMETERS:   
    !   refl_clr: reflectance of a clear pixel (fixed) 
    !   refl_cld: reflectance of clouds (fixed)
    !   i: cross-track pixel index (ip in m_cloud_pres_ret)
    !   iLine: along-swath scan row index
    !   i_obs_l, i_obs_s: normalized flux for longest (l) or shortest (s)
    !                     "good" wavelength
    !   i0_*: backscattered intensity at top of atmosphere?
    !   sb_*: fraction of intensity reflected by surface that is then
    !         scattered back to surface by atmosphere?
    !   tr_*: fractional transmittance of atmosphere?
    !   set_cld_frac: if true, calculate cloud fraction.
    !   min_refl, max_refl: minimum & maximum allowed reflectance (fixed)
    !   min_refl_flag: flag value for violations of min_refl (fixed)
    !   do_short_wave: include shortest wavelength bound in calculations?
    !   cal_reflec: calculate dIdR (radiance reflectance sensitivity)?
    !   iprt: verbosity level
    !
    ! !OUTPUT PARAMETERS:  
    !   refl: retreived reflectivity (array over whole swath)
    !   refl_l: retreived reflectivity in current pixel
    !   dIdR: Radiance (fractional) reflectance sensitivity
    !   qc: quality control array containing flags
    !   eff_cld_frac: effective cloud fraction
    !   rad_cld_frac: radiative cloud fraction
    !
    ! !SEE ALSO:  
    !   Joiner & Vasilkov (2006), IEEE transactions on geoscience and 
    !     remote sensing, 44, 1272, section III B.
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner     original fortran 90
    !  13Aug14  O'Sullivan  added documentation, some guesswork involved
    !
    !EOP
    !-------------------------------------------------------------------------
    !
    ! input/ouput variables
    real (KIND=8), intent(in) :: refl_clr, refl_cld
    integer, intent(in) :: i
    real (KIND=8), intent(in) :: i_obs_l, i_obs_s
    real (KIND=8), intent(in) :: i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
         sb_ls, tr_ls, i0_ss, sb_ss, tr_ss
    logical, intent(in) :: set_cld_frac
    ! local variables
    real (KIND=8) ::  I_clr_l, I_cld_l, I_clr_s, I_cld_s
    real (KIND=8) ::  eff_cld_frac_l

    !**************************************************************************

    !retrieve the reflectivity
    !=========================
    refl_l = (i_obs_l - i0_l)/(tr_l+(i_obs_l-i0_l)*Sb_l)
    if (cal_reflec) &
         dIdR(i,iLine)= &
         (tr_l/(1-refl_l*Sb_l)+refl_l*tr_l*Sb_l/ &
         (1-refl_l*Sb_l)**2)/i_obs_l
    refl(i,iLine)=refl_l

    !get effective cloud fraction from IPA method
    !============================================
    if (set_cld_frac) then
      I_clr_l=i0_ls + (refl_clr*tr_ls)/(1-refl_clr*Sb_ls)
      I_cld_l=i0_l + (refl_cld*tr_l)/(1-refl_cld*Sb_l)
      eff_cld_frac_l=(i_obs_l-I_clr_l)/(I_cld_l-I_clr_l)
      eff_cld_frac(i,iLine)= eff_cld_frac_l
      if (eff_cld_frac(i,iLine) > 1. .or. i_obs_l > I_cld_l) then
        eff_cld_frac(i,iLine) = 1.
      elseif (eff_cld_frac(i,iLine) < 0.) then
        eff_cld_frac(i,iLine) = 0.
        !! JJ take out for now
        ! else
        !refl(i,iLine)=refl_clr*(1-eff_cld_frac(i,iLine)) + refl_cld*eff_cld_frac(i,iLine)
      endif
      rad_cld_frac(i,iLine) = eff_cld_frac(i,iLine)*I_cld_l/i_obs_l
      if (rad_cld_frac(i,iLine) > 1. .or. i_obs_l > I_cld_l) then
        rad_cld_frac(i,iLine) = 1.
      elseif (rad_cld_frac(i,iLine) < 0.) then
        rad_cld_frac(i,iLine) = 0.
      endif

    endif ! set_cld_frac


    if (refl(i,iLine) < min_refl) then
      qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
    elseif (refl(i,iLine) > max_refl) then
      qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
    endif


    if (iprt >= 3) then
      print *, refl_clr, refl_cld
      print *, I_clr_l, I_cld_l, I_obs_l!, eff_cld_frac_l
      print *, (i_obs_l-I_clr_l), (I_cld_l-I_clr_l), &
           (i_obs_l-I_clr_l)/(I_cld_l-I_clr_l)
      print *,i0_ss + (refl_clr*tr_ss)/(1-refl_clr*Sb_ss),i0_s + &
           (refl_cld*tr_s)/(1-refl_cld*Sb_s)
      print *,'get_ai_refl: rad_cld_frac, eff_cld_frac, refl'
      print *, rad_cld_frac(i,iLine), eff_cld_frac(i,iLine), refl(i,iLine)
    endif


  end subroutine get_ai_refl

end module m_get_ai_refl
