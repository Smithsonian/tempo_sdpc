module m_cloud_pres_mod

INTEGER(KIND=4), EXTERNAL :: OMI_SMF_setmsg
integer :: status, NISE !, OMI_SMF_setmsg
real, parameter :: sz_max=86
real, parameter :: bad_thresh=0.02
real, parameter :: refl_cld_mask=0.4
real, parameter :: cld_mask_press_diff=-0.2
real, parameter :: sat_rad=6.5e13!6e13! this one worked! old value was 5.5e13
!7e13 too high
logical :: do_sat=.false. ! OPF8 fixed it!!!
integer :: nsat
integer :: nbad
integer, allocatable, dimension(:) :: bad_obs, good_obs
integer :: nterms=2 !2
integer :: nst !=nterms+3 ! + const, delta_lambda, cloud pres
integer :: ntm
real, dimension(:,:), allocatable :: x, x_fg
integer :: iter,nwave, nobs1 
integer :: nobs=0
integer :: ix1, ix2, ix1_old, ic1_old, ic1, ic2, ierr
integer :: i0x1, i0x2, i_np0
integer :: ip, i, ii, jj
integer :: i1_ler,i2_ler
real ::  int1_ler,  int2_ler
real :: ix1r
real :: delt_shift=-0.002
real :: noise=0.005
real :: noise2=0.10   
real :: delx=0.00001   
real :: delc=0.00001   
real :: diff_chi, shift_old
real :: chisq_old, chisq, bias, std, bias_old
real :: ixd, icd, np, np0, nc, ir, irc, irco, nt, j, l, nt_o, j_o, l_o
real :: refl_oc
integer :: ngood
integer, dimension(:), allocatable :: good
real, dimension(2) :: res
real, dimension(3) :: res3
real, dimension(:), allocatable :: chls
integer, dimension(:), allocatable :: ind
real, dimension(:), allocatable :: waves, y_obs, y_obs1, y_calc, y_frac, &
        y_calc_sh, y_calc_squeeze, sflx, sflx1, wave_diff
real, dimension(:), allocatable :: wavesd 
real, dimension(:,:), allocatable :: wavesp 
real, dimension(:,:), allocatable :: y_obs2, y_obs2b
real, dimension(:,:), allocatable :: y_resid, y_back
real, dimension(:), allocatable :: ycalc, rad_tot, rad_tot_noring, ref_rad, ref_ring
real, dimension(:), allocatable :: rad_clr, rad_cld, ring_cld2,ring_clr, ring_cld, jacob_dummy, jacob_rad, fit_rad
real, dimension(:), allocatable :: &
        r_i, r_i_sol, b_i, ring_oc, ring_sim, rad_sim, rad_clr_oc, radring_clr
real, dimension(:,:), allocatable :: h, htr, err_cov, corr
real, dimension(:,:), pointer :: ring_clds, rad_clds, ring_ocs, rad_clrs, ring_clrs, di_dr
real, dimension(:), allocatable :: w1p, f1p
real, dimension(:,:), allocatable :: i0a_clds, tra_clds, nra_clds, nia_clds, i01a_clds
real, dimension(:,:), allocatable :: z1_clds, z2_clds
real, dimension(:,:), allocatable :: ntot
real, dimension(:), allocatable :: o3_xsect
real, dimension(:,:,:,:), pointer :: table
real, dimension(:,:,:), pointer :: temp3D
real, dimension(:,:), pointer :: temp2D
real :: sz, satz, az

real :: wavetol=0.25 !0.22! nm tolerance for marking observations as bad due to bad qc flag
real :: var_inv_cp=1./0.5**2
real :: var_inv_chl=1./0.5**2 !1./5.0**2
real :: var_inv_big=1./1e8**2
real :: var_inv_small=1./1e-30**2
real :: var_inv_shift=1./0.05**2
real :: wave_min = 300.
real :: cld_frac, cld_frac_oc
real :: diff_chi_max=0.03
real :: diff_chi_iter_max=-0.05
real :: pres_sim=0.6
real :: i_obs_s, i_obs_l
real :: sflx_s, sflx_l
logical, dimension(:), allocatable :: computed, comp_clr, comp_all, comp_all_ring 
logical, dimension(:), allocatable :: computed_oc, comp_oc_clr
logical :: new_h, new_hcl, new_r
logical :: check_solar=.true.
logical :: check_rad=.true.
logical :: ret_chl
logical :: add_oc
logical :: comp_clear
logical :: bad_pix
real    :: chloro
real    :: adj
integer :: ino2
logical :: add_background_error = .true.
real    :: I0_l, I0_s, sb_l, sb_s, tr_l, tr_s
integer :: ind0, indt, indts, ind0s
integer :: nsh
real    :: reflec_cld
real    :: psurf, reflec
integer, parameter :: nbad_wavel=1
real, dimension(nbad_wavel) :: start_bad
real, dimension(nbad_wavel) :: end_bad
real, parameter :: small_frac=0.001
real, parameter :: ice_lat = 60.
real, parameter :: refl_ice = 0.7
integer :: niters
integer :: counts
integer, dimension(:), allocatable      :: indw
real :: i0_l1, i0_l2, sb_l1, sb_l2, tr_l1, tr_l2
real :: i0_s1, i0_s2, sb_s1, sb_s2, tr_s1, tr_s2
real :: i0_ls, sb_ls, tr_ls
real :: i0_ss, sb_ss, tr_ss
logical :: set_cld_frac
integer, parameter :: neta=1, nrefltm=0
integer, dimension(neta) :: order
real :: refl_land_fg=0.03, refl_oc_fg=0.08
real :: a
real :: contrast_thresh=0.002, contrast
real :: refl_clr_thresh=0.11
real :: refl_thresh=0.15
real :: a_thresh=5
real :: reflec1
integer :: eta_w
real :: cloud_term, dpdf, cld_pres2a, dring_dcld
real :: cldwgt, clrwgt, rad_tot_oc, ring_tot_oc
integer :: indf, i_r
logical :: do_o2_jacob=.false. !.true.
integer :: nchlr
integer :: ind1
real    :: elastic
end module m_cloud_pres_mod

