!> Interface module for ELSUNC optimizer
!! @file
module elsunc_interface_module
  use optimizer_interface_module
  use OMSAO_elsunc_fitting_module
  use tell_module
  implicit none

  public elsunc_optimizer
  private elsunc_objective
  private
  type(optimizer_type) :: this_optimizer

contains

  !> Objective function called by ELSUNC optimizer
  !! @details
  !! This function uses the parameters passed by the ELSUNC optimizer
  !! to call the objective function through the generic interface
  !! provided by \a optimizer_interface_module.
  !! @param[inout] params  Array of fit parameter values
  !! @param[in] num_params  Number of fit parameters
  !! @param[inout] residuals  Residuals from evaluating the objective function
  !! @param[in] num_residuals  Number of residuals
  !! @param[inout] elsunc_ctrl   ELSUNC integer control parameter
  !! @param[inout] cov_matrix   Covariance matrix
  !! @param[in]  dim1_cov_matrix  Leading dimension of covariance matrix array
  subroutine elsunc_objective (params, num_params, residuals, num_residuals, &
                               elsunc_ctrl, cov_matrix, dim1_cov_matrix)
    implicit none
    integer (kind=i4), intent(in) :: num_params, num_residuals
    real (kind=r8), dimension(num_params), intent(inout) :: params
    real (kind=r8), dimension(num_residuals), intent(inout) :: residuals
    integer (kind=i4), intent(inout) :: elsunc_ctrl
    integer (kind=i4), intent(in) :: dim1_cov_matrix
    real (kind=r8), dimension(dim1_cov_matrix,num_params), intent(inout) :: cov_matrix

    ! local variables
    integer (kind=i4) :: return_status, log_level
    character (len=1024) :: log_msg

    ! elsunc interprets the following return values of ctrl:
    integer (kind=i4), parameter :: UNCOMPUTABLE = -1
    integer (kind=i4), parameter :: JACOBIAN_NOT_AVAILABLE = 0

    ! 'this_optimizer' is a global

    if (elsunc_ctrl == 2) then
      this_optimizer%num_jac_calls = this_optimizer%num_jac_calls + 1
      if (this_optimizer%num_jac_calls < this_optimizer%max_num_fun_calls) then
        elsunc_ctrl = JACOBIAN_NOT_AVAILABLE
      else
        elsunc_ctrl = ELSUNC_INFLOOP_EVAL
      endif
      return
    endif

    this_optimizer%num_fun_calls = this_optimizer%num_fun_calls + 1
    if (this_optimizer%num_fun_calls > this_optimizer%max_num_fun_calls) then
      elsunc_ctrl = ELSUNC_INFLOOP_EVAL
      return
    endif

    if (.not.(allocated(this_optimizer%param_min) &
              .and.allocated(this_optimizer%param_max) &
              .and.allocated(this_optimizer%param_mask))) then
      call tell_error (tell_malloc_error, '*** elsunc_objective:  allocate failed', return_status)
      return
    endif

    if (any(params < this_optimizer%param_min) &
        .or. any(this_optimizer%param_max < params)) then
      elsunc_ctrl = UNCOMPUTABLE
      return
    endif

    call this_optimizer%objective (this_optimizer, params, num_params, &
                                   residuals, num_residuals, return_status)

    if (return_status < 0) then
      elsunc_ctrl = UNCOMPUTABLE
    endif

    log_level = tell_get_log_level()
    if (log_level > 4) then
      write(log_msg,'(1pe12.5,75(1x,1pe12.5))')sum(residuals(1:num_residuals)**2), &
        params(1:num_params)
      call tell_log (log_level, trim(log_msg))
    endif

  end subroutine elsunc_objective

  !> Generic Interface function to call the ELSUNC optimizer
  subroutine elsunc_optimizer (this, params, num_params, residuals, num_residuals, return_status, &
                               optional_cov_matrix)
   implicit none
   ! positional parameters
   type (optimizer_type) :: this
   real (kind=r8), dimension (:), intent(inout) :: params
   real (kind=r8), dimension (:),   intent(out) :: residuals
   integer (kind=i4),  intent(in) :: num_params, num_residuals
   integer (kind=i4), intent(out) :: return_status
   ! optional parameters
   real (kind=r8), dimension (:,:), intent(out), optional :: optional_cov_matrix

   !local

   ! Because of the way the elsunc code re-uses the covariance matrix array,
   ! the 2nd dimension of that array must have size >= 4. For details, look
   ! at the usage of the 'c' array in subroutine SECUC in the elsunc code.
   real (kind=r8), dimension(num_residuals,num_params + 4) :: cov_matrix

   ! elsunc specific objects
   integer (kind=i4), dimension (ELSUNC_NP) :: p
   real    (kind=r8), dimension (ELSUNC_NW) :: w
   integer (kind=i4) :: elbnd
   integer (kind=i4) :: elsunc_exval

   return_status = -1

   if ((0 < this%max_num_iterations) &
       .and. (this%max_num_iterations < (huge(1_i4)/num_params)/num_params)) then
     this%max_num_fun_calls = this%max_num_iterations * num_params * num_params
   else if ((0 < this%max_num_fun_calls) &
            .and. (this%max_num_fun_calls < huge(1_i4))) then
     this%max_num_iterations = (this%max_num_fun_calls/num_params)/num_params
   else
     this%max_num_fun_calls = huge(1_i4)
     this%max_num_iterations = huge(1_i4)
   endif

   ! Initialize elsunc control parameter arrays:
   p = -1;
   p(1) = 0
   p(3) = this%max_num_iterations
   w = -1.0
   w(1:4) = [real (kind=r8) :: this%tol, this%epsrel, this%epsabs, this%epsx]

   ! use a global to pass 'this' structure to elsunc_objective
   this_optimizer = this

   if (this%mode == opt_unbounded) then
     elbnd = 0
   else if (this%mode == opt_bounded) then
     elbnd = 2
   else
     call tell_error (tell_invalid_parm, "elsunc_optimizer: unsupported bounds type", return_status)
     return
   endif

   elsunc_exval = 0
   call elsunc (params, num_params, num_residuals, num_residuals, &
                elsunc_objective, elbnd, this%param_min, this%param_max, &
                p, w, elsunc_exval, residuals, cov_matrix)

   ! save the number of iterations
   this%num_iterations = p(6)

   ! map elsunc return code range onto the generic set
   if (elsunc_exval >= ELSUNC_LESS_IS_NOISE) then
     return_status = opt_convergence_good
   else if (0 <= elsunc_exval .and. elsunc_exval < ELSUNC_LESS_IS_NOISE) then
     return_status = opt_convergence_suspect
   else if (elsunc_exval == ELSUNC_MAXITER_EVAL) then
     return_status = opt_convergence_maxiter_exceeded
   else
     return_status = opt_convergence_failed
   endif

   if (present(optional_cov_matrix)) then
     optional_cov_matrix(1:num_residuals,1:num_params) = cov_matrix(1:num_residuals,1:num_params)
   endif

 end subroutine elsunc_optimizer
end module elsunc_interface_module
