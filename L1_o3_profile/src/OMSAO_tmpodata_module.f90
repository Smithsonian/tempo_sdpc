!
MODULE OMSAO_tmpodata_module
   USE OMSAO_precision_module
   USE OMSAO_parameters_module, ONLY: maxchlen, n_rad_winwav, maxwin, &
       max_fit_pts, max_ring_pts, mrefl
   USE OMSAO_indices_module,    ONLY: n_max_fitpars, max_rs_idx, &
       max_calfit_idx, o3_t1_idx, o3_t3_idx, spc_idx
   USE OMSAO_variables_module, ONLY:ring_group, refl_group, rad_group, & 
       irrad_group,cali_group, geo_group, o3p_group
   IMPLICIT NONE
   

   ! ---------------------------------
   ! Maximum data/swath dimensions
   ! ---------------------------------
   INTEGER, PARAMETER :: mswath_tmpo = 2
   INTEGER (KIND=i4), PARAMETER :: &
       ntimes_max     = 129, nxtrack_max  = 2048, nlines_max=10
   INTEGER, PARAMETER :: nwavel_ccd = 1026, nwavel_max=1026+1026


   ! ---------------------------------
   ! swath names (UV and VIS)
   ! --------------------------------
   CHARACTER (LEN=maxchlen), DIMENSION(mswath_tmpo), PARAMETER :: & 
            rad_swathname = (/ 'band_290_490_nm', 'band_540_740_nm'/)
   character (len=maxchlen), dimension(mswath_tmpo) :: &
           irrad_swathname = (/ 'band_290_490_nm', 'band_540_740_nm'/)

   ! ------------------------------------------------------------
   ! Boundary wavelengths (approximate) for UV and VIS channels
   ! ------------------------------------------------------------
   REAL (KIND=dp), DIMENSION(mswath_tmpo),  PARAMETER :: &
        lower_wvls_tmpo = (/290.0, 490.0/), upper_wvls_tmpo = (/540.0,740.0/) 
   REAL (KIND=dp), PARAMETER :: lower_spec_tmpo = 0.0, upper_spec_tmpo = 10.0E14 


   TYPE(ring_group) :: tmpo_ring
   TYPE(irrad_group) :: tmpo_irrad
   TYPE(rad_group) :: tmpo_rad
   TYPE(geo_group) :: tmpo_geo1, tmpo_geo2, tmpo_geo3, tmpo_geo4, tmpo_geo5
   TYPE(refl_group) :: tmpo_refl
   TYPE(cali_group):: tmpo_cali
   TYPE(o3p_group):: tmpo_o3p
END MODULE OMSAO_tmpodata_module
