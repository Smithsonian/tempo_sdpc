module m_get_ai_refl

contains

  subroutine get_ai_refl(refl_clr, refl_cld, iter, i_obs_l, i_obs_s, i, &
       i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
       sb_ls, tr_ls, set_cld_frac, i0_ss, sb_ss, tr_ss)

    use m_vars, ONLY: &
         refl, ai, iprt, refl_l, dIdR, & 
         min_refl, max_refl, min_refl_flag, ai_flag, qc, &
         max_ai, eff_cld_frac, do_short_wave, iLine, cal_reflec, rad_cld_frac

    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  get_ai_refl
    ! 
    ! !DESCRIPTION: get_ai_refl reads precomputed tables needed for cloud
    !               parameter retrievals		
    !
    ! !CALLING SEQUENCE: 
    !
    !        call get_ai_refl
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
    real (KIND=8), intent(in) :: refl_clr, refl_cld

    integer, intent(in) :: iter, i
    real (KIND=8), intent(in) :: i_obs_l, i_obs_s
    real (KIND=8), intent(in) :: i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
         sb_ls, tr_ls, i0_ss, sb_ss, tr_ss
    logical, intent(in) :: set_cld_frac

    real (KIND=8) :: i_ray_l, i_ray_s, I_clr_l, I_cld_l, I_clr_s, I_cld_s
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

    !retrieve the aerosol index
    !===========================
    if (do_short_wave) then
      i_ray_l=i_obs_l
      if (eff_cld_frac(i,iLine) == 0. .or. eff_cld_frac(i,iLine) == 1.) then
        ! calculated short wavelength using retrieved reflectivity
        i_ray_s=i0_s + (refl_l*tr_s)/(1-refl_l*Sb_s)
      elseif (eff_cld_frac(i,iLine) < 1.) then
        ! calculated short wavelength using retrieved effective cloud fraction
        I_clr_s=i0_ss + (refl_clr*tr_ss)/(1-refl_clr*Sb_ss)
        I_cld_s=i0_s + (refl_cld*tr_s)/(1-refl_cld*Sb_s)
        i_ray_s=i_clr_s*(1-eff_cld_frac(i,iLine))+i_cld_s*eff_cld_frac(i,iLine)
      endif
      ai(i,iLine)=-100.*(log10(i_obs_s/i_obs_l) - log10(i_ray_s/i_ray_l))
      if (ai(i,iLine) > max_ai) then
        qc(i,iLine) = IBSET(qc(i,iLine),ai_flag)
      endif
    endif

    if (iprt >= 3) then
      print *, refl_clr, refl_cld
      print *, I_clr_l, I_cld_l, I_obs_l!, eff_cld_frac_l
      print *, (i_obs_l-I_clr_l), (I_cld_l-I_clr_l), &
           (i_obs_l-I_clr_l)/(I_cld_l-I_clr_l)
      print *, I_clr_s, I_cld_s, I_obs_s!, eff_cld_frac_l
      print *,i0_ss + (refl_clr*tr_ss)/(1-refl_clr*Sb_ss),i0_s + &
           (refl_cld*tr_s)/(1-refl_cld*Sb_s)
      print *,'get_ai_refl: rad_cld_frac, eff_cld_frac, refl'
      print *, rad_cld_frac(i,iLine), eff_cld_frac(i,iLine), refl(i,iLine)
    endif


  end subroutine get_ai_refl

end module m_get_ai_refl
