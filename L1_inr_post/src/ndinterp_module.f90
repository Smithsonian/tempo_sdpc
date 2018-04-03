module ndinterp_module
  use tell_module
  implicit none
  private

  ! dimension(1) varies fastest.
  ! Value array is indexed as val(d1,d2,...,dn)

  ! Recommend compiling this with -O3 to aggressively inline functions.
  ! With a gfortran test doing 10^6 trials interpolating on a large 6-D
  ! table, -O2 runs in 9 sec, -O3 runs in 5 sec.

  public ndi_find_index, ndi_find_indices
  public ndi_calc_weights, ndi_calc_term_weight
  public ndi_dims_alloc, ndi_dims_dealloc
  public ndi_table_interp

  include 'ndinterp_decl.inc'

  integer, parameter :: &
    i1 = selected_int_kind (2**1), &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4), &
    r4 = kind(1.0), &
    r8 = selected_real_kind (2*precision(1.0_r4))

  type, public :: ndi_dim
    real (kind=r8), allocatable, dimension(:) :: x
    integer :: dimlen
    character (len=32) :: name
  end type ndi_dim

contains

  integer function ndi_find_index (x0, x)
    implicit none
    real (kind=r8), intent(in) :: x0
    real (kind=r8), dimension(:), intent(in) :: x
    integer :: i, n

    ! Assume size(n) >= 2 (not tested)
    !
    ! On return, we want to ensure that 'i' satisfies
    !     x(i) <= x0 < x(i+1)
    ! such that both i and i+1 are valid array indices.
    !
    ! Extrapolation is not supported.
    n = size(x)

    if (x0 < x(1) .or. x(n) < x0) then
      ndi_find_index = -1
      return
    endif

    ! For long coordinate arrays (e.g. dimlen>64?), binary search
    ! will be faster. However, at the moment I'm most interested
    ! in support for interpolation in high-dimensional tables,
    ! where the coordinate arrays mostly have dimlen<64.

    i = minloc (abs(x0 - x), dim=1)

    if (x0 < x(i) .or. i == n) then
      i = i - 1
    endif

    ndi_find_index = i
  end function ndi_find_index

  subroutine ndi_find_indices (dims, x, indices, errstat)
    implicit none
    type (ndi_dim), dimension(:), intent(in) :: dims
    real (kind=r8), dimension(:), intent(in) :: x
    integer, dimension(:), intent(out) :: indices
    integer, intent(inout) :: errstat
    character (len=72) :: msg
    integer :: d

    if (errstat /= 0) return

    do d = 1, size(dims)
      indices(d) = ndi_find_index (x(d), dims(d) % x)
      if (indices(d) < 0) then
        write(msg, '(a,1pe12.4)')'x = ',x(d)
        call tell_error (tell_runtime_error, &
                         "ndi_find_indices: value out of range"//msg, &
                         errstat)
        return
      endif
    enddo
  end subroutine ndi_find_indices

  subroutine ndi_calc_weights (dims, x, indices, weights, errstat)
    implicit none
    type (ndi_dim), target, dimension(:), intent(in) :: dims
    real (kind=r8), dimension(:), intent(in) :: x
    integer, dimension(:), intent(in) :: indices
    real (kind=r8), dimension(:), intent(out) :: weights
    integer, intent(inout) :: errstat

    real (kind=r8), dimension(:), pointer :: tx
    integer :: d, m

    if (errstat /= 0) return

    ! Linear interpolation

    do d = 1, size(dims)
      tx => dims(d) % x
      m = indices(d)
      weights(d) = (tx(m+1) - x(d)) / (tx(m+1) - tx(m))

      if (weights(d) < 0 .or. weights(d) > 1) then
        call tell_error (tell_runtime_error, &
                         "ndi_calc_weights: extrapolation is not supported", &
                         errstat)
        return
      endif
    enddo
  end subroutine ndi_calc_weights

  subroutine ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
    implicit none
    integer, intent(in) :: k
    real (kind=r8), dimension(:), intent(in) :: weights
    integer, dimension(:), intent(in) :: indices
    integer, dimension(:), intent(out) :: id_k
    real (kind=r8), intent(out) :: wt_k
    integer :: d

    ! The interpolation point, x(:), lies inside an N-dimensional table
    ! cell with one corner at indices(d).
    ! The "diagonally opposite" corner of that cell is at indices(d)+1.
    ! In terms of the table coordinate arrays, X => tbl % dims(d) % x,
    ! we have
    !            X(index(d)) <= x(d) < X(index(d)+1).
    !
    ! The linearly interpolated value at x(:) may be viewed as a weighted
    ! sum of table values at that cell's 2^N corners.
    ! The coefficient of the kth term in that sum is given by the product
    ! of N weights, one weight factor for each dimension.
    ! The weights are defined so that the corner with index(d)
    ! carries weight, weights(d), (e.g. equality => weights(d)=1)
    ! and the corner with index(d)+1 carries weight, (1 - weights(d)).
    !
    ! For the kth term in the sum, this loop generates both
    ! the term's lookup table array indices, id_k(:), and
    ! the term's weight coefficient, 'wt_k'.

    wt_k = 1.0
    do d = 1, size(weights)
      if (iand (k-1, ishft(1, d-1)) == 0) then
        id_k(d) = indices(d) + 1
        wt_k = wt_k * (1.0 - weights(d))
      else
        id_k(d) = indices(d)
        wt_k = wt_k * weights(d)
      endif
    enddo
  end subroutine ndi_calc_term_weight

  subroutine ndi_dims_dealloc (dims)
    implicit none
    type (ndi_dim), dimension(:), intent(inout) :: dims
    integer :: i, err

    do i = 1,size(dims)
      if (allocated (dims(i) % x)) then
        deallocate (dims(i) % x, stat=err)
        if (err /= 0) return
      endif
    enddo
  end subroutine ndi_dims_dealloc

  subroutine ndi_dims_alloc (dims, dimlens, errstat)
    implicit none
    type (ndi_dim), dimension(:), intent(inout) :: dims
    integer, dimension(:), intent(in) :: dimlens
    integer, intent(inout) :: errstat
    integer :: i, err

    do i = 1, size(dims)
      if (dimlens(i) .lt. 2) then
        write(*,'(a,i0)')'*** Error: unsupported dimension size = ',dimlens(i)
        errstat = -1
        return
      endif

      dims(i) % dimlen = dimlens(i)

      allocate (dims(i) % x(dimlens(i)), stat=err)
      if (err /= 0) then
        errstat = -1
        return
      endif
    enddo
  end subroutine ndi_dims_alloc

  include 'ndinterp_code.inc'

end module ndinterp_module
