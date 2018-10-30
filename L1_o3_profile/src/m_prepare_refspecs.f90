!
module m_prepare_refspecs

  public prepare_refspecs
  private

contains

  SUBROUTINE prepare_refspecs (n_radpts, curr_rad_wvl, pge_error_status )

    ! ***********************************************************
    !
    !   Calculate the splined fitting database.
    !   Note that the undersampled spectrum has just been done.
    !   Finish filling in database array.
    !
    ! ***********************************************************

    USE OMSAO_precision_module
    USE OMSAO_indices_module,    ONLY: solar_idx, wvl_idx, spc_idx!, &
    !us1_idx, us2_idx, ring_idx, ring1_idx, mxs_idx, max_rs_idx, max_calfit_idx
    !USE OMSAO_parameters_module, ONLY: maxchlen
    USE OMSAO_variables_module,  ONLY: curr_sol_spec, &
         database, n_irrad_wvl!, yn_doas, yn_smooth, up_radbnd, &
         !refsol_idx, rad_winwav_idx, n_refspec, n_refspec_pts, &
         !lo_radbnd, database_shiwf
    USE OMSAO_errstat_module
    USE m_ezspline_interpolation, only: interpolation

    IMPLICIT NONE

    ! *******************************************************************
    ! CAREFUL: Assumes that radiance and solar wavelength arrays have the
    ! same number of points. That must not be the case if we read in a
    ! general EL1 file. Examine and adjust! (tpk, note to himself)
    ! *******************************************************************

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                              INTENT (IN)    :: n_radpts
    REAL (KIND=dp), DIMENSION (n_radpts), INTENT (IN)    :: curr_rad_wvl
    INTEGER,                              INTENT (INOUT) :: pge_error_status

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: errstat!, i, j, npts, lu_rad, ll_rad
    REAL  (KIND=dp), DIMENSION (n_radpts) :: spline_sun
    real (kind=dp), dimension(n_irrad_wvl) :: tmp_x_in, tmp_y_in

    ! ------------------
    ! External functions
    ! ------------------
    !INTEGER :: OMI_SMF_setmsg

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    !CHARACTER (LEN=16), PARAMETER :: modulename = 'prepare_refspecs'


    ! Spline irradiance spectrum onto radiance grid

    ! PROBLEM: first wavlength in POS is smaller than first wavelength in
    !          curr_sol_spec(wvl_idx,*).
    ! SOLUTION: Don't include POS(1) in the interpolation, and assign 
    !           SPLINE_SUN(1) = SPEC_SUN(1)

    errstat = pge_errstat_ok

! FIXME - masking array temporaries
    tmp_x_in=curr_sol_spec(wvl_idx,1:n_irrad_wvl)
    tmp_y_in=curr_sol_spec(spc_idx,1:n_irrad_wvl)
    CALL interpolation (n_irrad_wvl, &
!         curr_sol_spec(wvl_idx,1:n_irrad_wvl), &
!         curr_sol_spec(spc_idx,1:n_irrad_wvl), &
         tmp_x_in, tmp_y_in, &
         n_radpts, curr_rad_wvl(1:n_radpts), &
         spline_sun(1:n_radpts), errstat )  
    database(solar_idx, 1:n_radpts) = spline_sun(1:n_radpts) 
    ! =================================================================
    ! Note that the UNDERSAMPLING spectrum has already been assigned to
    ! DATABASE(us1/2_idx,*) in the UNDERSPEC routine.
    ! =================================================================


    pge_error_status = MAX ( errstat, pge_error_status )

    RETURN
  END SUBROUTINE prepare_refspecs
end module m_prepare_refspecs
