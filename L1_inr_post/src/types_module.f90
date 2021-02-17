!> Definition of data types shared by polarization correction code
!! @file

module types_module
  use, intrinsic :: iso_c_binding
  use tio_module
  use ndinterp_module

  integer, public, parameter :: maxscene = 2

  integer, public, parameter :: &
    i1 = selected_int_kind (2**1), &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4), &
    r4 = kind(1.0), &
    r8 = selected_real_kind (2*precision(1.0_r4))

  !! Wavelength limits [nm] spanning the full TEMPO band,
  !! and a bit more to be on the safe side.
  real (kind=r8), parameter :: &
    lut_wav_min = 290.0, &
    lut_wav_max = 743.0

  !> Data type for storing diagnostic output
  type, public :: diag_output_type
    ! shape: (wave, xtrack) for specified mirror step
    real (kind=r8), allocatable, dimension(:,:) :: q, u
    ! shape: (xtrack) for specified mirror step
    real (kind=r8), allocatable, dimension(:) :: delta_irp
    logical :: active
  end type

  !> Radiance data structure
  type, public :: radiance_type
    type(tiof_file_type) :: obj
    integer :: this_step                        !< step index in the current granule
    integer :: num_step, num_xtrack, num_wave   !< granule dimensions
    ! shape: (wave, xtrack) for specified mirror step
    real (kind=r8), allocatable, dimension(:,:) :: radiance  !< radiance
    real (kind=r8), allocatable, dimension(:,:) :: wave      !< wavelength
    ! shape: (xtrack, step)
    real (kind=r8), allocatable, dimension(:,:) :: lon, lat  !< longitude, latitude
    real (kind=r8), allocatable, dimension(:,:) :: sza, saa  !< solar zenith, azimuth angles
    real (kind=r8), allocatable, dimension(:,:) :: vza, vaa  !< viewing zenith, azimuth angles
    real (kind=r8), allocatable, dimension(:,:) :: raa       !< relative azimuth angle
    real (kind=r8), allocatable, dimension(:,:) :: hgt       !< terrain height
    real (kind=r8), allocatable, dimension(:,:) :: pre       !< surface pressure
    integer, allocatable, dimension(:,:) :: inr_quality_flag !< INR quality flag
  end type

  !< data type to support subsetting of a single coordinate array
  type, public :: range_type
    real (kind=r8) :: min, max   !< floating point min/max coordinate values
    integer :: imin, imax        !< integer array indices corresponding to min/max
  end type

  !< data type to support subsetting the polarization lookup table on input
  type, public :: radiance_subset_type
    type(range_type) :: lon, lat   !< longitude, latitude subset
    type(range_type) :: sza        !< solar zenith angle subset
    type(range_type) :: vza        !< viewing zenith angle subset
    type(range_type) :: raa        !< relative azimuth angle subset
    type(range_type) :: wav        !< wavelength subset
  end type

  interface
    function spline (x,y,n, xs,ys,ns) bind (c, name='spline')
      use, intrinsic :: iso_c_binding, only : c_double, c_int
      implicit none
      integer (c_int), value :: n, ns
      real (kind=c_double), dimension(n), intent(in) :: x, y
      real (kind=c_double), dimension(ns), intent(in) :: xs
      real (kind=c_double), dimension(ns), intent(out) :: ys
      integer (c_int) :: spline
    end function spline
  end interface

end module types_module
