module optimizer_interface_module
  USE OMSAO_precision_module, only : i4, r8
  implicit none

  ! This module is intended to provide a generic optimizer interface
  !
  ! Example:
  !
  !  type (optimizer_type) :: opt
  !  call optimizer_open (opt, optimizer_method, objective, num_params, return_status, &
  !                       param_min=pmin, param_max=pmax, param_mask=mask, &
  !                       mode=mode, tol=tol, epsrel=epsrel, epsabs=epsabs, epsx=epsx, &
  !                       max_num_fun_calls=max_nfc,
  !                       max_num_iterations=max_itnum)
  !  call opt%optimize (opt, params, num_params, residuals, num_residuals, return_status, &
  !                     cov_matrix=c)
  !  call optimizer_close (opt, return_status)
  !
  !  For optimizer_open and optimizer_close:
  !                return_status < 0 means failure
  !                return_status >= 0 means success
  !  opt%optimize depends on the particular implementation,
  !  but should adhere to the same convention.

  public

  ! exit codes
  integer (kind=i4), parameter :: &
    opt_convergence_failed=-2, &
    opt_convergence_maxiter_exceeded=-1, &
    opt_convergence_suspect=0, &
    opt_convergence_good=1

  ! modes
  integer (kind=i4), parameter :: opt_unbounded=0, opt_bounded=1

  type optimizer_type
    procedure(optimizer_interface), nopass, pointer :: optimize
    procedure(objective_interface), nopass, pointer :: objective
    real    (kind=r8) :: tol, epsrel, epsabs, epsx
    real    (kind=r8), dimension(:), allocatable :: param_min, param_max
    integer (kind=i4), dimension(:), allocatable :: param_mask
    integer (kind=i4) :: mode
    integer (kind=i4) :: num_params
    integer (kind=i4) :: num_iterations, max_num_iterations
    integer (kind=i4) :: num_fun_calls, max_num_fun_calls
    integer (kind=i4) :: num_jac_calls
  end type optimizer_type

  interface
    subroutine optimizer_interface (this, params, num_params, residuals, num_residuals, return_status, &
                                   cov_matrix)
     import i4, r8, optimizer_type
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

  subroutine optimizer_set_default_method (optimizer_method)
    implicit none
    procedure(optimizer_interface) :: optimizer_method
    default_optimizer_method => optimizer_method
  end subroutine optimizer_set_default_method

  subroutine optimizer_close (this, return_status)
    implicit none
    type(optimizer_type) :: this
    integer (kind=i4), intent(out) :: return_status
    if (allocated(this%param_mask)) deallocate (this%param_mask)
    if (allocated(this%param_min)) deallocate (this%param_min)
    if (allocated(this%param_max)) deallocate (this%param_max)
    return_status = 0
  end subroutine optimizer_close

  subroutine optimizer_open (this, objective, num_params, return_status, &
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
    integer (kind=i4), intent(out) :: return_status
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

    return_status = -1

    if (num_params < 1) then
      write(*,*)'optimizer_open:  invalid number of parameters'
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
      write(*,*)'optimizer_open:  allocate failed'
      return
    endif

    if (present(param_min)) then
      if (size(param_min) < num_params) then
        write(*,*)'optimizer_open:  invalid param_min array'
        return
      endif
      this%param_min = param_min
    else
      this%param_min = -huge(1.0_r8)
    endif

    if (present(param_max)) then
      if (size(param_max) < num_params) then
        write(*,*)'optimizer_open:  invalid param_max array'
        return
      endif
      this%param_max = param_max
    else
      this%param_max = huge(1.0_r8)
    endif

    if (present(param_mask)) then
      if (size(param_mask) < num_params) then
        write(*,*)'optimizer_open:  invalid param_mask array'
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

    return_status = 0
  end subroutine optimizer_open

end module optimizer_interface_module
