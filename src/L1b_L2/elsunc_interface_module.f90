module elsunc_interface_module
  use optimizer_interface_module
  use OMSAO_elsunc_fitting_module
  implicit none

  public elsunc_optimizer
  private elsunc_objective
  private
  type(optimizer_type) :: this_optimizer

contains

  subroutine elsunc_objective (params, num_params, residuals, num_residuals, &
                               elsunc_ctrl, cov_matrix, dim1_cov_matrix)
    implicit none
    real (kind=r8), dimension(num_params), intent(inout) :: params
    real (kind=r8), dimension(num_residuals), intent(out) :: residuals
    integer (kind=i4), intent(in) :: num_params, num_residuals
    integer (kind=i4), intent(inout) :: elsunc_ctrl
    integer (kind=i4), intent(in) :: dim1_cov_matrix
    real (kind=r8), dimension(dim1_cov_matrix,num_params), intent(out) :: cov_matrix

    ! local variables
    integer (kind=i4) :: return_status

    ! elsunc interprets the following return values of ctrl:
    integer (kind=i4), parameter :: UNCOMPUTABLE = 1
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
      write(*,*)'*** elsunc_objective:  internal error'
      stop
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

  end subroutine elsunc_objective

 subroutine elsunc_optimizer (this, params, num_params, residuals, num_residuals, return_status, &
                              optional_cov_matrix)
   use OMSAO_indices_module, only: elsunc_userdef
   use OMSAO_variables_module, only: tol, epsrel, epsabs, epsx
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
   real (kind=r8), dimension(num_residuals,num_params) :: cov_matrix

   ! elsunc specific objects
   integer (kind=i4), dimension (ELSUNC_NP) :: p
   real    (kind=r8), dimension (ELSUNC_NW) :: w
   integer (kind=i4) :: elsunc_return_status

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
   w(1:4) = [real (kind=r8) :: tol, epsrel, epsabs, epsx]

   ! use a global to pass 'this' structure to elsunc_objective
   this_optimizer = this

   elsunc_return_status = 0
   call elsunc (params, num_params, num_residuals, num_residuals, &
                elsunc_objective, elsunc_userdef, &
                this%param_min, this%param_max, &
                p, w, elsunc_return_status, &
                residuals, cov_matrix)

   ! save the number of iterations
   this%num_iterations = p(6)

   return_status = elsunc_return_status

   if (present(optional_cov_matrix)) then
     optional_cov_matrix(1:num_residuals,1:num_params) = cov_matrix(1:num_residuals,1:num_params)
   endif

 end subroutine elsunc_optimizer

end module elsunc_interface_module
