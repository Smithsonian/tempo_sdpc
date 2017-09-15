!>Calculate relative cloud fraction using MLER method
module m_get_ai_refl

  private
  public get_ai_refl

contains

  !--------------------------------------------------------------------
  !
  !> DESCRIPTION: get_ai_refl calculates radiative cloud fraction, 
  !>   cloud reflectivity, aerosol index, etc., using the 
  !>   Mixed Lambert-Equivalent Reflectivity (MLER) concept and
  !>   Independent Pixel Approximation (IPA). Assumes that a partly 
  !>   cloudy pixel is the sum of a mix of cloudy and clear sub-pixels.\n
  !
  !> REFERENCE:\n  
  !>   Joiner & Vasilkov (2006), IEEE transactions on geoscience and 
  !>     remote sensing, 44, 1272, section III B.
  !-------------------------------------------------------------------------
  ! !INPUT PARAMETERS:   
  !> @param refl_clr  reflectance of a clear pixel (fixed) 
  !> @param refl_cld  reflectance of clouds (fixed)
  !> @param i  cross-track pixel index (ip in m_cloud_pres_ret)
  !  iLine  along-swath scan row index
  !> @param i_obs_l  normalized flux for longest (l) "good" wavelength
  !> @param i_obs_s  normalized flux for shortest (s) "good" wavelength
  !> @param i0_l  backscattered intensity at top of atmosphere?
  !> @param i0_s  backscattered intensity at top of atmosphere?
  !> @param sb_l  fraction of intensity reflected by surface that is then
  !>         scattered back to surface by atmosphere?
  !> @param sb_s  fraction of intensity reflected by surface that is then
  !>         scattered back to surface by atmosphere?
  !> @param tr_l  fractional transmittance of atmosphere?
  !> @param tr_s  fractional transmittance of atmosphere?
  !> @param set_cld_frac  if true, calculate cloud fraction.
  !  min_refl, max_refl  minimum & maximum allowed reflectance (fixed)
  !  min_refl_flag  flag value for violations of min_refl (fixed)
  !  do_short_wave  include shortest wavelength bound in calculations?
  !  cal_reflec  calculate dIdR (radiance reflectance sensitivity)?
  !
  ! !OUTPUT PARAMETERS:  
  !  refl  retreived reflectivity (array over whole swath)
  !  refl_l  retreived reflectivity in current pixel
  !  dIdR  Radiance (fractional) reflectance sensitivity
  !  qc  quality control array containing flags
  !  eff_cld_frac  effective cloud fraction
  !  rad_cld_frac  radiative cloud fraction
  !
  !
  ! !REVISION HISTORY: 
  !
  !> @author  05Jan01   Joiner     original fortran 90
  !> @author 13Aug14  O'Sullivan  added documentation, 
  !>                              some guesswork involved
  !
  !-------------------------------------------------------------------------

  subroutine get_ai_refl(refl_clr, refl_cld, i_obs_l, i_obs_s, i, &
       i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
       sb_ls, tr_ls, set_cld_frac, i0_ss, sb_ss, tr_ss)

    use m_vars, ONLY: refl, refl_l, dIdR, min_refl, max_refl, min_refl_flag, &
         qc, eff_cld_frac, iLine, cal_reflec, rad_cld_frac
    use tell_module

    implicit none

    ! input/ouput variables
    real (KIND=8), intent(in) :: refl_clr, refl_cld
    integer, intent(in) :: i
    real (KIND=8), intent(in) :: i_obs_l, i_obs_s
    real (KIND=8), intent(in) :: i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
         sb_ls, tr_ls, i0_ss, sb_ss, tr_ss
    logical, intent(in) :: set_cld_frac
    ! local variables
    real (KIND=8) ::  I_clr_l, I_cld_l
    real (KIND=8) ::  eff_cld_frac_l
    character (len=128) :: logmsg

    !**************************************************************************


    !retrieve the reflectivity
    !=========================
    refl_l = (i_obs_l - i0_l)/(tr_l+(i_obs_l-i0_l)*Sb_l)
    if (cal_reflec) &
         dIdR(i,iLine)= real(&
         (tr_l/(1-refl_l*Sb_l)+refl_l*tr_l*Sb_l/ &
         (1-refl_l*Sb_l)**2)/i_obs_l &
         , kind=4)
    refl(i,iLine)=real(refl_l, kind=4)

    !get effective cloud fraction from IPA method
    !============================================
    if (set_cld_frac) then
      I_clr_l=i0_ls + (refl_clr*tr_ls)/(1-refl_clr*Sb_ls)
      I_cld_l=i0_l + (refl_cld*tr_l)/(1-refl_cld*Sb_l)
      eff_cld_frac_l=(i_obs_l-I_clr_l)/(I_cld_l-I_clr_l)
      eff_cld_frac(i,iLine)= real(eff_cld_frac_l, kind=4)
      if (eff_cld_frac(i,iLine) > 1. .or. i_obs_l > I_cld_l) then
        eff_cld_frac(i,iLine) = 1.
      elseif (eff_cld_frac(i,iLine) < 0.) then
        eff_cld_frac(i,iLine) = 0.
        !! JJ take out for now
        ! else
        !refl(i,iLine)=refl_clr*(1-eff_cld_frac(i,iLine)) + refl_cld*eff_cld_frac(i,iLine)
      endif
      rad_cld_frac(i,iLine) = &
           real(eff_cld_frac(i,iLine)*I_cld_l/i_obs_l , kind=4) 
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


    write(logmsg,*) refl_clr, refl_cld
    call tell_log(4,logmsg)
    if (set_cld_frac) then
      write(logmsg,*) I_clr_l, I_cld_l, I_obs_l
      call tell_log(4,logmsg)
      write(logmsg,*) (i_obs_l-I_clr_l), (I_cld_l-I_clr_l), &
           (i_obs_l-I_clr_l)/(I_cld_l-I_clr_l)
      call tell_log(4,logmsg)
      write(logmsg,*) i0_ss + (refl_clr*tr_ss)/(1-refl_clr*Sb_ss),i0_s + &
           (refl_cld*tr_s)/(1-refl_cld*Sb_s)
      call tell_log(4,logmsg)
      call tell_log(4,' get_ai_refl: rad_cld_frac, eff_cld_frac, refl')
      write(logmsg,*) rad_cld_frac(i,iLine), eff_cld_frac(i,iLine), &
           refl(i,iLine)
      call tell_log(4,logmsg)
    endif


  end subroutine get_ai_refl

end module m_get_ai_refl
