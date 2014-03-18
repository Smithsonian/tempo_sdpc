module m_get_f

contains

subroutine get_f(refl_clr, refl_cld, iter, i_obs_l, i_obs_s, i, &
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
real, intent(in) :: refl_clr
real, intent(out):: refl_cld
character(len=50) :: myname='get_f: '

integer, intent(in) :: iter, i
real, intent(in) :: i_obs_l, i_obs_s
real, intent(in) :: i0_l, i0_s, sb_l, sb_s, tr_l, tr_s, i0_ls, sb_ls, tr_ls, &
  i0_ss, sb_ss, tr_ss
logical, intent(in) :: set_cld_frac

real :: i_ray_l, i_ray_s, I_clr_l, I_cld_l, I_clr_s, I_cld_s
real :: I_clr_l2, I_cld_l2
real ::  eff_cld_frac_l, eff_cld_frac_l2
real ::  ratio_obs, ratio_clr
logical :: do_short

!**************************************************************************

  refl_l = (i_obs_l - i0_l)/(tr_l+(i_obs_l-i0_l)*Sb_l)
  if (cal_reflec) &
   dIdR(i,iLine)= &
    (tr_l/(1-refl_l*Sb_l)+refl_l*tr_l*Sb_l/ &
    (1-refl_l*Sb_l)**2)/i_obs_l
  refl(i,iLine)=refl_l

  I_clr_l=i0_ls + (refl_clr*tr_ls)/(1-refl_clr*Sb_ls)
  I_clr_s=i0_ss + (refl_clr*tr_ss)/(1-refl_clr*Sb_ss)

if (set_cld_frac) then

  !get observed ratio
  !==================
  ratio_obs=i_obs_s/i_obs_l

  !compute clear scene ratio
  !=========================
  ratio_clr= I_clr_s/I_clr_l

  !retrieve the opacity fraction
  !=============================
  eff_cld_frac(i,iLine)=(ratio_obs-ratio_clr)/(1-ratio_clr) 
  if (eff_cld_frac(i,iLine) > 1.) then
     eff_cld_frac(i,iLine) = 1.
  elseif (eff_cld_frac(i,iLine) < 0.) then
     eff_cld_frac(i,iLine) = 0.
  endif

  !get radiative cloud fraction
  !============================
  I_cld_l=(i_obs_l-i_clr_l*(1-eff_cld_frac(i,iLine)))/eff_cld_frac(i,iLine)
  rad_cld_frac(i,iLine) = eff_cld_frac(i,iLine)*I_cld_l/i_obs_l

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

! use refl_cld as QC check for now
! if (refl_cld > 1.) then
!    refl_cld = 1.
!    qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
! elseif (refl_cld < 0.) then
!    refl_cld = 0.
!    qc(i,iLine) = IBSET(qc(i,iLine),min_refl_flag)
! endif

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
     I_cld_s=i0_s + (refl_cld*tr_s)/(1-refl_cld*Sb_s)
     i_ray_s=i_clr_s*(1-eff_cld_frac(i,iLine))+i_cld_s*eff_cld_frac(i,iLine)
   endif
   ai(i,iLine)=-100.*(alog10(i_obs_s/i_obs_l) - alog10(i_ray_s/i_ray_l))
   if (ai(i,iLine) > max_ai) then
     qc(i,iLine) = IBSET(qc(i,iLine),ai_flag)
   endif
  endif

  if (iprt >= 3) then
   !print *, i0_l, i0_s, sb_l, sb_s, tr_l, tr_s
   !print *, i_obs_l, i_obs_s
    write(6,104) 
    write(6,101) refl_clr, refl_cld
    write(6,105) 
    write(6,101) ratio_clr, ratio_obs
   !write(6,101) I_clr_l, I_cld_l, I_obs_l!, eff_cld_frac_l
   !write(6,101) i0_ls, tr_ls, Sb_ls
   !write(6,101) i0_l, tr_l, Sb_l
   !write(6,101) (i_obs_l-I_clr_l), (I_cld_l-I_clr_l), (i_obs_l-I_clr_l)/(I_cld_l-I_clr_l)
   !write(6,101) I_clr_s, I_cld_s, I_obs_s!, eff_cld_frac_l
   !write(6,100) eff_cld_frac_l
   write(6,103) 
   write(6,100) refl(i,iLine), ai(i,iLine), eff_cld_frac(i,iLine)
   !write(6,100) rad_cld_frac(i,iLine), eff_cld_frac(i,iLine), refl(i,iLine)
  endif

102 format('get_f: rad_cld_frac, eff_cld_frac, refl')
103 format('get_f: refl, ai, eff_cld_frac')
104 format('get_f: refl_clr, refl_cld')
105 format('get_f: ratio_clr, ratio_obs')
100 format(4f12.3)
101 format(4e12.3)

end subroutine get_f

end module m_get_f
