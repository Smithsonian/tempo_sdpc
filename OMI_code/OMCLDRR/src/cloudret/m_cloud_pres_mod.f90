module m_cloud_pres_mod

INTEGER(KIND=4), EXTERNAL :: OMI_SMF_setmsg
integer :: status, NISE !, OMI_SMF_setmsg
real (KIND=8), parameter :: sz_max=86d0
real (KIND=8), parameter :: bad_thresh=0.02d0
real (KIND=8), parameter :: refl_cld_mask=0.4d0
real (KIND=8), parameter :: cld_mask_press_diff=-0.2d0
real (KIND=8), parameter :: sat_rad=6.5e13!6e13! this one worked! old value was 5.5e13
!7e13 too high
logical :: do_sat=.false. ! OPF8 fixed it!!!
integer :: nsat
integer :: nbad
integer, allocatable, dimension(:) :: bad_obs, good_obs
integer :: nterms=2 !2
integer :: nst !=nterms+3 ! + const, delta_lambda, cloud pres
integer :: ntm
real (KIND=8), dimension(:,:), allocatable :: x, x_fg
integer :: iter,nwave, nobs1 
integer :: nobs=0
integer :: ix1, ix2, ix1_old, ic1_old, ic1, ic2, ierr
integer :: i0x1, i0x2, i_np0
integer :: ip, i, ii, jj
integer :: i1_ler,i2_ler
real (KIND=8) ::  int1_ler,  int2_ler
real (KIND=8) :: ix1r
real (KIND=8) :: delt_shift=-0.002d0
real (KIND=8) :: noise=0.005d0
real (KIND=8) :: noise2=0.10d0   
real (KIND=8) :: delx=0.00001d0   
real (KIND=8) :: delc=0.00001d0   
real (KIND=8) :: diff_chi, shift_old
real (KIND=8) :: chisq_old, chisq, bias, std, bias_old
real (KIND=8) :: ixd, icd, np, np0, nc, ir, irc, irco, nt, j, l, nt_o, j_o, l_o
real (KIND=8) :: refl_oc
integer :: ngood
integer, dimension(:), allocatable :: good
real (KIND=8), dimension(2) :: res
real (KIND=8), dimension(3) :: res3
real (KIND=8), dimension(:), allocatable :: chls
integer, dimension(:), allocatable :: ind
real (KIND=8), dimension(:), allocatable :: waves, y_obs, y_obs1, y_calc, y_frac, &
        y_calc_sh, y_calc_squeeze, sflx, sflx1, wave_diff
real (KIND=8), dimension(:), allocatable :: wavesd 
real (KIND=8), dimension(:,:), allocatable :: wavesp 
real (KIND=8), dimension(:,:), allocatable :: y_obs2, y_obs2b
real (KIND=8), dimension(:,:), allocatable :: y_resid, y_back
real (KIND=8), dimension(:), allocatable :: ycalc, rad_tot, rad_tot_noring, ref_rad, ref_ring
real (KIND=8), dimension(:), allocatable :: rad_clr, rad_cld, ring_cld2,ring_clr, ring_cld, jacob_dummy, jacob_rad, fit_rad
real (KIND=8), dimension(:), allocatable :: &
        r_i, r_i_sol, b_i, ring_oc, ring_sim, rad_sim, rad_clr_oc, radring_clr
real (KIND=8), dimension(:,:), allocatable :: h, htr, err_cov, corr
real (KIND=8), dimension(:,:), pointer :: ring_clds, rad_clds, ring_ocs, rad_clrs, ring_clrs, di_dr
real (KIND=8), dimension(:), allocatable :: w1p, f1p
real (KIND=8), dimension(:,:), allocatable :: i0a_clds, tra_clds, nra_clds, nia_clds, i01a_clds
real (KIND=8), dimension(:,:), allocatable :: z1_clds, z2_clds
real (KIND=8), dimension(:,:), allocatable :: ntot
real (KIND=8), dimension(:), allocatable :: o3_xsect
real (KIND=8), dimension(:,:,:,:), pointer :: table
!real (KIND=4), dimension(:,:,:,:), pointer :: table
real (KIND=8), dimension(:,:,:), pointer :: temp3D
!real (KIND=4), dimension(:,:,:), pointer :: temp3D
real (KIND=8), dimension(:,:), pointer :: temp2D
real (KIND=8) :: sz, satz, az

real (KIND=8) :: wavetol=0.25 !0.22! nm tolerance for marking observations as bad due to bad qc flag
real (KIND=8) :: var_inv_cp=1./0.5d0**2
real (KIND=8) :: var_inv_chl=1./0.5d0**2 !1./5.0**2
real (KIND=8) :: var_inv_big=1./1e8**2
real (KIND=8) :: var_inv_small=(1.d0/1e-30)**2 
real (KIND=8) :: var_inv_shift=1./0.05d0**2
real (KIND=8) :: wave_min = 300.d0
real (KIND=8) :: cld_frac, cld_frac_oc
real (KIND=8) :: diff_chi_max=0.03d0
real (KIND=8) :: diff_chi_iter_max=-0.05d0
real (KIND=8) :: pres_sim=0.6d0
real (KIND=8) :: i_obs_s, i_obs_l
real (KIND=8) :: sflx_s, sflx_l
logical, dimension(:), allocatable :: computed, comp_clr, comp_all, comp_all_ring 
logical, dimension(:), allocatable :: computed_oc, comp_oc_clr
logical :: new_h, new_hcl, new_r
logical :: check_solar=.true.
logical :: check_rad=.true.
logical :: ret_chl
logical :: add_oc
logical :: comp_clear
logical :: bad_pix
real (KIND=8)    :: chloro
real (KIND=8)   :: adj
integer :: ino2
logical :: add_background_error = .true.
real (KIND=8)   :: I0_l, I0_s, sb_l, sb_s, tr_l, tr_s
integer :: ind0, indt, indts, ind0s
integer :: nsh
real (KIND=8)   :: reflec_cld
real (KIND=8)   :: psurf, reflec
integer, parameter :: nbad_wavel=1
real (KIND=8), dimension(nbad_wavel) :: start_bad
real (KIND=8), dimension(nbad_wavel) :: end_bad
real (KIND=8), parameter :: small_frac=0.001d0
real (KIND=8), parameter :: ice_lat = 60.d0
real (KIND=8), parameter :: refl_ice = 0.7d0
integer :: niters
integer :: counts
integer, dimension(:), allocatable      :: indw
real (KIND=8) :: i0_l1, i0_l2, sb_l1, sb_l2, tr_l1, tr_l2
real (KIND=8) :: i0_s1, i0_s2, sb_s1, sb_s2, tr_s1, tr_s2
real (KIND=8) :: i0_ls, sb_ls, tr_ls
real (KIND=8) :: i0_ss, sb_ss, tr_ss
logical :: set_cld_frac
integer, parameter :: neta=1, nrefltm=0
integer, dimension(neta) :: order
real (KIND=8) :: refl_land_fg=0.03d0, refl_oc_fg=0.08d0
real (KIND=8) :: a
real (KIND=8) :: contrast_thresh=0.002d0, contrast
real (KIND=8) :: refl_clr_thresh=0.11d0
real (KIND=8) :: refl_thresh=0.15d0
real (KIND=8) :: a_thresh=5.d0
real (KIND=8) :: reflec1
integer :: eta_w
real (KIND=8) :: cloud_term, dpdf, cld_pres2a, dring_dcld
real (KIND=8) :: cldwgt, clrwgt, rad_tot_oc, ring_tot_oc
integer :: indf, i_r
logical :: do_o2_jacob=.false. !.true.
integer :: nchlr
integer :: ind1
real (KIND=8)    :: elastic
end module m_cloud_pres_mod

