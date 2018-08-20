!> Optimizer interface module
!! @file
!! @details
!! This module is intended to provide a generic optimizer interface
!!
!! Example:
!!
!! @code
!!  type (optimizer_type) :: opt
!!  call optimizer_open (opt, objective, num_params, errstat, &
!!                       optimizer_method=method, &
!!                       param_min=pmin, param_max=pmax, param_mask=mask, &
!!                       mode=mode, tol=tol, epsrel=epsrel, epsabs=epsabs, epsx=epsx, &
!!                       max_num_fun_calls=max_nfc,
!!                       max_num_iterations=max_itnum)
!!  call opt%optimize (opt, params, num_params, residuals, num_residuals, return_status, &
!!                     cov_matrix=c)
!!  call optimizer_close (opt, errstat)
!! @endcode
!
!!  For \a optimizer_open and \a optimizer_close:
!!  @verbatim
!!                errstat /= 0 means failure
!!                errstat = 0 means success
!!  @endverbatim
!!  The value of \a opt%optimize \a return_status depends on the particular
!!  optimizer's implementation, but should adhere to the same convention.
module optimizer_interface_module
  USE OMSAO_precision_module, only : i4, r8
  use tell_module
  implicit none

  public

  !> exit codes
  integer (kind=i4), parameter :: &
    opt_convergence_failed=-2, &
    opt_convergence_maxiter_exceeded=-1, &
    opt_convergence_suspect=0, &
    opt_convergence_good=1

  !> optimization modes
  integer (kind=i4), parameter :: opt_unbounded=0, opt_bounded=1

  !> optimizer object type definition
  type optimizer_type
    !> pointer to optimization procedure
    procedure(optimizer_interface), nopass, pointer :: optimize
    !> pointer to procedure to evaluate the objective function
    procedure(objective_interface), nopass, pointer :: objective
    !> Absolute and relative convergence tolerances
    real    (kind=r8) :: tol, epsrel, epsabs, epsx
    !> Allowed range for each fit parameter
    real    (kind=r8), dimension(:), allocatable :: param_min, param_max
    !> Variable fit parameter index array
    integer (kind=i4), dimension(:), allocatable :: param_mask
    integer (kind=i4) :: mode
    integer (kind=i4) :: num_params
    integer (kind=i4) :: num_iterations, max_num_iterations
    integer (kind=i4) :: num_fun_calls, max_num_fun_calls
    integer (kind=i4) :: num_jac_calls
  end type optimizer_type

  !> optimizer interface
  interface
    subroutine optimizer_interface (this, params, num_params, residuals, num_residuals, return_status, &
                                   cov_matrix)
     import i4, r8, optimizer_type
     ! 'implicit none' was commented out to work around a bug
     ! in ifort 12.0.2 20110112
     ! See https://software.intel.com/en-us/forums/topic/337047
     ! The bug was fixed in ifort-13.0.
     implicit none
     ! positional parameters
     type(optimizer_type) :: this
     real (kind=r8), dimension (:), intent(inout) :: params
     real (kind=r8), dimension (:),   intent(out) :: residuals
     integer (kind=i4),  intent(in) :: num_params, num_residuals
     integer (kind=i4), intent(out) :: return_status
     ! optional parameters
     real (kind=r8), dimension (:,:), intent(out), optional :: cov_matrix
   end subroutine optimizer_interface
 end interface

 !> objective function interface
 interface
   subroutine objective_interface (this, params, num_params, residuals, num_residuals, return_status)
     import i4, r8, optimizer_type
     implicit none
     type(optimizer_type) :: this
     real (kind=r8), dimension (:), intent(in) :: params
     real (kind=r8), dimension (:), intent(out) :: residuals
     integer (kind=i4), intent(in) :: num_params, num_residuals
     integer (kind=i4), intent(out) :: return_status
   end subroutine objective_interface
 end interface

 procedure(optimizer_interface), private, pointer :: default_optimizer_method => null()

contains

  !> Define the default optimization method
  !! @param[in] optimizer_method  Optimization procedure.
  subroutine optimizer_set_default_method (optimizer_method)
    implicit none
    procedure(optimizer_interface) :: optimizer_method
    default_optimizer_method => optimizer_method
  end subroutine optimizer_set_default_method

  !> Release any resources associated with this \a optimizer_type instance
  !! @param[inout] this  The current \a optimizer_type instance
  !! @param[inout] errstat Error status variable
  subroutine optimizer_close (this, errstat)
    implicit none
    type(optimizer_type) :: this
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    if (allocated(this%param_mask)) deallocate (this%param_mask, stat=errstat)
    if (allocated(this%param_min) .and. errstat /= 0) &
         deallocate (this%param_min, stat=errstat)
    if (allocated(this%param_max) .and. errstat /= 0) &
         deallocate (this%param_max, stat=errstat)

    if (errstat /= 0) then
      call tell_error(tell_malloc_error,"optimizer_close: deallocate failed",&
           errstat)
      return
    endif
  end subroutine optimizer_close

  !> Initialize an instance of \a optimizer_type
  !! @param[inout] this  The current \a optimizer_type instance
  !! @param[in] objective  The objective function to be optimized
  !! @param[in] num_params  Number of fit parameters
  !! @param[inout] errstat Error status variable
  !! @param[in]   optimizer_method  (Optional) optimization procedure.
  !!                 If not specified, the default optimization procedure
  !!                 will be used.
  !! @param[in]  mode  (Optional) Optimizer mode selection.
  !! @param[in]  tol   (Optional) Convergence tolerance.
  !!                   Precise usage depends on the optimizer selection.
  !! @param[in]  epsabs   (Optional) Absolute convergence tolerance.
  !!                   Precise usage depends on the optimizer selection.
  !! @param[in]  epsrel   (Optional) Relative convergence tolerance.
  !!                   Precise usage depends on the optimizer selection.
  !! @param[in]  epsx   (Optional) Parameter value relative convergence tolerance.
  !!                   Precise usage depends on the optimizer selection.
  !! @param[in]  param_min[] (Optional) Minimum allowed parameter value.
  !! @param[in]  param_max[] (Optional) Minimum allowed parameter value.
  !! @param[in]  param_mask[] (Optional) Array of non-frozen parameter indices.
  !! @param[in]  max_num_fun_calls  (Optional) Maximum allowed number of
  !!                        of times the objective function may be evaluated.
  !! @param[in]  max_num_iterations  (Optional) Maximum number of times the
  !!                       optimizer's inner loop may be executed.
  !!                   Precise usage depends on the optimizer selection.
  subroutine optimizer_open (this, objective, num_params, errstat, &
                             optimizer_method, &
                             mode, tol, epsabs, epsrel, epsx, &
                             param_min, param_max, param_mask, &
                             max_num_fun_calls, &
                             max_num_iterations)
    implicit none
    ! positional parameters
    type(optimizer_type) :: this
    procedure(objective_interface) :: objective
    integer (kind=i4), intent(in) :: num_params
    integer, intent(inout) :: errstat
    ! optional parameters
    procedure(optimizer_interface), optional :: optimizer_method
    integer (kind=i4), intent(in), optional :: mode
    real    (kind=r8), intent(in), optional :: tol, epsabs, epsrel, epsx
    real    (kind=r8), dimension(:), intent(inout), optional :: param_min, param_max
    integer (kind=i4), dimension(:), intent(in), optional :: param_mask
    integer (kind=i4), intent(in), optional :: max_num_fun_calls
    integer (kind=i4), intent(in), optional :: max_num_iterations

    ! local variables
    integer (kind=i4) :: i, status

    if (num_params < 1) then
      call tell_error (tell_invalid_parm, "optimizer_open:  invalid number of parameters", errstat)
      return
    endif

    if (present(optimizer_method)) then
      this%optimize => optimizer_method
    else
      this%optimize => default_optimizer_method
    endif

    this%objective => objective
    this%num_params = num_params

    if (present(mode)) then
      this%mode = mode
    else
      if (present(param_min).or.present(param_max)) then
        this%mode = opt_bounded
      else
        this%mode = opt_unbounded
      endif
    endif
    if (present(tol)) then
      this%tol = tol
    else
      this%tol = -1.0_r8
    endif
    if (present(epsrel)) then
      this%epsrel = epsrel
    else
      this%epsrel = -1.0_r8
    endif
    if (present(epsabs)) then
      this%epsabs = epsabs
    else
      this%epsabs = -1.0_r8
    endif
    if (present(epsx)) then
      this%epsx = epsx
    else
      this%epsx = -1.0_r8
    endif

    allocate (this%param_min(num_params), stat=status)
    if (status == 0) allocate (this%param_max(num_params), stat=status)
    if (status == 0) allocate (this%param_mask(num_params), stat=status)
    if (status /= 0) then
      call tell_error (tell_malloc_error, "optimizer_open:  allocate failed", errstat)
      return
    endif

    if (present(param_min)) then
      if (size(param_min) < num_params) then
        call tell_error (tell_invalid_parm, "optimizer_open:  invalid param_min array", errstat)
        return
      endif
      this%param_min = param_min
    else
      this%param_min = -huge(1.0_r8)
    endif

    if (present(param_max)) then
      if (size(param_max) < num_params) then
        call tell_error (tell_invalid_parm, "optimizer_open:  invalid param_max array", errstat)
        return
      endif
      this%param_max = param_max
    else
      this%param_max = huge(1.0_r8)
    endif

    if (present(param_mask)) then
      if (size(param_mask) < num_params) then
        call tell_error (tell_invalid_parm, "optimizer_open:  invalid param_mask array", errstat)
        return
      endif
      this%param_mask = param_mask
    else
      this%param_mask = [integer(kind=i4) :: (i, i=1_i4,num_params)]
    endif

    if (present(max_num_fun_calls)) then
      this%max_num_fun_calls = max_num_fun_calls
    else
      this%max_num_fun_calls = huge(1_i4)
    endif

    if (present(max_num_iterations)) then
      this%max_num_iterations = max_num_iterations
    else
      this%max_num_iterations = huge(1_i4)
    endif

    this%num_iterations = 0
    this%num_fun_calls = 0
    this%num_jac_calls = 0

  end subroutine optimizer_open

end module optimizer_interface_module
