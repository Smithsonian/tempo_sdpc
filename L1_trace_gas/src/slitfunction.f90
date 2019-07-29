!> Generic slit-function object
!! @file
!! @details
!! This module defines a generic interface for convolving
!! a spectrum with a slit function.
module slitfunction
  use tell_module
  use slitfunction_asym_gaussian
  use slitfunction_super_gaussian
  use OMSAO_precision_module, only : i4, r8
  implicit none

  interface
    subroutine sf_init_interface (errstat)
      implicit none
      integer, intent(inout) :: errstat
    end subroutine sf_init_interface
  end interface

  interface
    subroutine sf_convolve_interface (pixel, num_wvl, wvl, spec, spec_conv, errstat)
      import i4, r8
      implicit none
      integer (kind=i4), intent(in) :: pixel, num_wvl
      real (kind=r8), dimension(:), intent(in) :: wvl, spec
      real (kind=r8), dimension(:), intent(out) :: spec_conv
      integer, intent(inout) :: errstat
    end subroutine sf_convolve_interface
  end interface

  procedure (sf_init_interface), private, pointer :: sf_init => null()
  procedure (sf_convolve_interface), private, pointer :: sf_convolve => null()

contains

  subroutine slitfunction_select (sf_init_method, sf_convolve_method)
    implicit none
    procedure (sf_init_interface) :: sf_init_method
    procedure (sf_convolve_interface) :: sf_convolve_method
    sf_init => sf_init_method
    sf_convolve => sf_convolve_method
  end subroutine slitfunction_select

  subroutine slitfunction_open (errstat, use_table)
    implicit none
    integer, intent(inout) :: errstat
    logical, optional, intent(in) :: use_table

    if (errstat /= 0) return

    if (present(use_table)) then
      if (.not.use_table) return
    endif

    call sf_init (errstat)
    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "slitfunction_open: failed", errstat)
      return
    endif

  end subroutine slitfunction_open

  subroutine slitfunction_convolve (num_wvl, wvl, spec, spec_conv, &
                                    pixel, params, num_params, errstat)
    use ctrlvars, only: yn_use_labslitfunc
    implicit none
    integer (kind=i4), intent(in) :: num_wvl
    real (kind=r8), dimension(:), intent(in) :: wvl, spec
    real (kind=r8), dimension(:), intent(out) :: spec_conv
    integer (kind=i4), intent(in) :: pixel, num_params
    real (kind=r8), dimension(:), intent(in) :: params
    integer, intent(inout) :: errstat

    if (errstat /= 0) return

    if (yn_use_labslitfunc) then
      call sf_convolve (pixel, num_wvl, wvl, spec, spec_conv, errstat)
    else
      if (num_params /= 2) then
        call tell_error (tell_invalid_parm, "slitfunction_convolve: analytic slit function requires num_params=2", errstat)
        return
      endif
      call asymmetric_gaussian_sf (num_wvl, params(1), params(2), wvl, spec, spec_conv)
! asymmetric_gaussian_sf ( npoints, hw1e, e_asym, wvlarr, specarr, specmod)
!      call super_gaussian_sf (num_wvl, params(1), 10.0_r8 , wvl, spec, spec_conv)
    endif

    if (errstat /= 0) then
      call tell_error (tell_runtime_error, "slitfunction_convolve: convolution failed", errstat)
      return
    endif

  end subroutine slitfunction_convolve

end module slitfunction
