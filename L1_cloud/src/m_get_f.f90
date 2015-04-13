!>Calculate radiative cloud fraction using LER method
module m_get_f

  private 
  public get_f

contains

  !-------------------------------------------------------------------------
  !> DESCRIPTION: get_f calculates radiative cloud fraction, 
  !>    cloud reflectivity, aerosol index, etc., using the 
  !>    Lambert-Equivalent Reflectivity (LER) concept. Assumes
  !>    cloud entirely fills pixel, so cloud fraction relates
  !>    to fraction of light that manages to pass through that
  !>    cloud layer.\n
  !
  !> REFERENCE:  \n
  !>   Joiner & Vasilkov (2006), IEEE transactions on geoscience and 
  !>     remote sensing, 44, 1272, section III A.
  !
  ! !INPUT PARAMETERS:   
  !> @param   refl_clr  reflectance of a clear pixel (fixed) 
  !> @param   refl_cld  reflectance of clouds (fixed)
  !  i  cross-track pixel index (ip in m_cloud_pres_ret)
  !  iLine  along-swath scan row index
  !> @param   i_obs_l, i_obs_s  normalized flux for longest (l) or shortest (s)
  !>                     "good" wavelength
  !> @param   i0_*  backscattered intensity at top of atmosphere?
  !> @param   sb_*  fraction of intensity reflected by surface that is then
  !>         scattered back to surface by atmosphere?
  !> @param   tr_*  fractional transmittance of atmosphere?
  !  set_cld_frac  if true, calculate cloud fraction.
  !  min_refl, max_refl  minimum & maximum allowed reflectance (fixed)
  !  min_refl_flag  flag value for violations of min_refl (fixed)
  !  do_short_wave  include shortest wavelength bound in calculations?
  !  cal_reflec  calculate dIdR (radiance reflectance sensitivity)?
  !  iprt  verbosity level
  !
  ! !OUTPUT PARAMETERS:  
  !  refl  retreived reflectivity (array over whole swath)
  !  refl_l  retreived reflectivity in current pixel
  !  dIdR  Radiance (fractional) reflectance sensitivity
  !  qc  quality control array containing flags
  !  eff_cld_frac  effective cloud fraction
  !  rad_cld_frac  radiative cloud fraction
  !
  ! !REVISION HISTORY: 
  !
  !> @author  05Jan01   Joiner     original fortran 90
  !> @author  13Aug14  O'Sullivan  added documentation
  !
  !-------------------------------------------------------------------------

  subroutine get_f(refl_clr, refl_cld, i_obs_l, i_obs_s, i, &
       i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
       sb_ls, tr_ls, set_cld_frac, i0_ss, sb_ss, tr_ss)

    use m_vars, ONLY: &
         refl, iprt, refl_l, dIdR, min_refl, max_refl, min_refl_flag, &
         qc, eff_cld_frac, iLine, &
         cal_reflec, rad_cld_frac 

    implicit none

    !input/output variables
    real (KIND=8), intent(in) :: refl_clr
    real (KIND=8), intent(out):: refl_cld
    integer, intent(in) :: i
    real (KIND=8), intent(in) :: i_obs_l, i_obs_s
    real (KIND=8), intent(in) :: i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, &
         sb_ls, tr_ls, i0_ss, sb_ss, tr_ss
    logical, intent(in) :: set_cld_frac
    !local variables
    real (KIND=8) :: I_clr_l, I_cld_l, I_clr_s
    real (KIND=8) ::  ratio_obs, ratio_clr

    !**************************************************************************

    refl_l = (i_obs_l - i0_l)/(tr_l+(i_obs_l-i0_l)*Sb_l)
    if (cal_reflec) &
         dIdR(i,iLine)= real( &
         (tr_l/(1-refl_l*Sb_l)+refl_l*tr_l*Sb_l/ &
         (1-refl_l*Sb_l)**2)/i_obs_l &
         , kind=4)
    refl(i,iLine)=real(refl_l, kind=4)

    I_clr_l=i0_ls + (refl_clr*tr_ls)/(1-refl_clr*Sb_ls)
    I_clr_s=i0_ss + (refl_clr*tr_ss)/(1-refl_clr*Sb_ss)
    !to avoid an unitialized variable error
    I_cld_l=I_clr_l

    if (set_cld_frac) then

      !get observed ratio
      !==================
      ratio_obs=i_obs_s/i_obs_l

      !compute clear scene ratio
      !=========================
      ratio_clr= I_clr_s/I_clr_l

      !retrieve the opacity fraction
      !=============================
      eff_cld_frac(i,iLine)=real((ratio_obs-ratio_clr)/(1-ratio_clr), kind=4)
      if (eff_cld_frac(i,iLine) > 1.) then
        eff_cld_frac(i,iLine) = 1.
      elseif (eff_cld_frac(i,iLine) < 0.) then
        eff_cld_frac(i,iLine) = 0.
      endif

      !get radiative cloud fraction
      !============================
      I_cld_l=(i_obs_l-i_clr_l*(1-eff_cld_frac(i,iLine)))/eff_cld_frac(i,iLine)
      rad_cld_frac(i,iLine) = &
           real(eff_cld_frac(i,iLine)*I_cld_l/i_obs_l, kind=4)

      !Alternatives
      !============
      !rad_cld_frac(i,iLine) = (i_obs_l-i_clr_l*(1-eff_cld_frac(i,iLine)))/i_obs_l
      !I_cld_l=rad_cld_frac(i,iLine)*i_obs_l/eff_cld_frac(i,iLine)

      !I_cld_l=i0_l + (refl_cld*tr_l)/(1-refl_cld*Sb_l) !old LER method

      if (rad_cld_frac(i,iLine) > 1. .or. i_obs_l > I_cld_l) then
        rad_cld_frac(i,iLine) = 1.
        qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
      elseif (rad_cld_frac(i,iLine) < 0.) then
        rad_cld_frac(i,iLine) = 0.
        qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
      endif

    endif ! set_cld_frac

    !get the cloud reflectivity
    !==========================
    refl_cld = (i_cld_l - i0_l)/(tr_l+(i_cld_l-i0_l)*Sb_l)

    if (refl(i,iLine) < min_refl) then
      qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
    elseif (refl(i,iLine) > max_refl) then
      qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
    endif


    if (iprt >= 3) then
      print *,'get_f: refl_clr, refl_cld'
      print *, refl_clr, refl_cld
      print *,'get_f: ratio_clr, ratio_obs'
      print *, ratio_clr, ratio_obs
      print *,'get_f: refl, ai, eff_cld_frac'
      print *, refl(i,iLine), eff_cld_frac(i,iLine)
    endif


  end subroutine get_f

end module m_get_f
