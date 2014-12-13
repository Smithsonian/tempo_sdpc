module output_tools
  use terr_module
  use tio_module
  USE OMSAO_precision_module
  use ctrlvars, only: yn_diagnostic_run, yn_refseccor, yn_scat_weights

  implicit none
  private

  public create_template

contains

  subroutine create_template (filename, num_steps, num_xtrack, num_swlevels, &
                              errstat)
    USE OMSAO_parameters_module, ONLY: NWAVEL_MAX, nUTCdim
    USE OMSAO_indices_module,   ONLY: max_calfit_idx, max_rs_idx
    USE OMSAO_omidata_module,   ONLY: nclenfit, n_comm_wvl
    USE OMSAO_variables_module, ONLY: n_fitvar_rad
    use netcdf

    implicit none
    character (len=*), intent(in) :: filename
    integer (kind=i4), intent(in) :: num_steps, num_xtrack, num_swlevels
    integer, intent(inout) :: errstat

    type (tiof_object_type) :: obj
    type (tiof_dimlist_type) :: dimlist
    type (tiof_varlist_type) :: varlist
    integer, dimension(2) :: dimids_xtrack_step

    if (errstat < 0) return

    call tiof_dimlist_append (dimlist, tempo_dim_step, num_steps, errstat)
    call tiof_dimlist_append (dimlist, tempo_dim_xtrack, num_xtrack, errstat)
    call tiof_dimlist_append (dimlist, tempo_dim_swt_levels, num_swlevels, errstat)

    call tiof_open (filename, obj, nf90_write, errstat)
    call tiof_def_dims (obj, dimlist, errstat)

    call tiof_dimlist_lookup (dimlist, 2, &
                              [tempo_dim_xtrack, tempo_dim_step], &
                              dimids_xtrack_step, &
                              errstat)

    call tiof_varlist_append (varlist, &
                              tempo_var_column_amount, &  ! name
                              nf90_double,             &  ! xtype
                              dimids_xtrack_step,      &  ! dimids
                              errstat)

    call tiof_def_vars (obj, varlist, errstat)

    if (errstat < 0) then
      call terr_error (terr_io_error, &
                       "create_template: intializing file "//trim(filename), &
                       errstat)
      return
    endif

  end subroutine create_template

end module output_tools
