!> Multi-dimensional linear interpolation
!! @file
!! @sa ndinterp_code.inc  (Table interpolation functions for 1-D to 7-D tables)
!! @sa ndinterp_decl.inc  (Common interface definition for table interpolation functions)
!! @sa ndinterp_module_code.in  (Template for autogeneration of table interpolation functions)
!! @sa ndinterp_module_autogen.sh  (Shell script to autogenerate table interpolation functions)
!!
!! @details
!!
!! These functions support linear interpolation on lookup tables
!! implemented as a multidimensional array plus a 1-D array of
!! coordinates for each array dimension.  The coordinate arrays
!! are managed in an array of objects of type \a ndi_dim,
!! such that dimension 1 corresponds to the fastest varying
!! (leftmost) index of the multidimensional table value array,
!! and dimension N corresponds to the slowest varying (rightmost)
!! index of the multidimensional table value array.  Coordinate
!! arrays are assumed to be stored in ascending order.
!! Extrapolation is not supported.
!!
!! **Error handling:**
!!
!! Most member functions take an integer error code, \a errstat
!! When the input \a errstat is negative, the function returns
!! immediately, with \a errstat unchanged. On return, \a errstat<0
!! indicates that an error occured. Error messages and global error
!! status are handled using libtell. Additional information on the
!! type of error that occurred may be available from the global error
!! status, via \a tell_get_error and \a tell_copy_strerror.

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

  !> Linearly interpolate in a multidimensional array
  !! @fn subroutine ndi_table_interp (dims, tbl, x, val, errstat, input_indices, output_indices)
  !! @param[in]  dims  Array of N \a ndi_dim objects
  !! @param[in]  tbl   N-dimensional array of table values
  !! @param[in]  x     Coordinates of interpolation point
  !! @param[out] val   Interpolated table value at \a x.
  !! @param[inout] errstat   Error status code
  !! @param[optional, in]   input_indices   Indices of the table cell known to contain the interpolation point, \a x
  !! @param[optional, out]  output_indices  Indices of the table cell containing the interpolation point, \a x
  include 'ndinterp_decl.inc'

  integer, parameter :: &
    i1 = selected_int_kind (2**1), &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4), &
    r4 = kind(1.0), &
    r8 = selected_real_kind (2*precision(1.0_r4))

  !> Coordinate grid object
  type, public :: ndi_dim
    real (kind=r8), allocatable, dimension(:) :: x  !< coordinate array
    integer :: dimlen                               !< size of this array dimension
    character (len=32) :: name                      !< name of this coordinate
  end type ndi_dim

contains

  !> Find coordinate x0 in array x(:), assumed to be in ascending order
  !! @param[in]   x0   Coordinate value to search for
  !! @param[in]   x    Coordinate array to search. Must be in ascending order
  !! @return on error, returns -1.
  !!  on success, returns array index \a i, such that $x(i) \le x0 \lt x(i+1)$.
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

  !> Find the table cell containing coordinate vector x(:)
  !! @param[in]  dims     Array of N \a ndi_dim coordinate objects
  !! @param[in]  x        Vector of N coordinates of the interpolation point
  !! @param[out] indices  Vector of array indices of the table cell
  !!                      containing the interpolation point, \a x.
  !! @param[inout]  errstat  Error status
  !! @return on success, errstat=0, on error, errstat /=0
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

  !> Compute linear interpolation weight for each array dimension
  !! @param[in]  dims    Array of N \a ndi_dim coordinate objects
  !! @param[in]  x       Vector of N coordinates of the interpolation point
  !! @param[in] indices  Vector of array indices of the table cell
  !!                     containing the interpolation point, \a x.
  !! @param[out] weights Vector of N linear interpolation weights.
  !! @param[inout]  errstat  Error status
  !! @return on success, errstat=0, on error, errstat /=0
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
      if (m < 1 .or. (m+1) > dims(d) % dimlen) then
        call tell_error (tell_runtime_error, &
                         "ndi_calc_weights: invalid array index", errstat)
        return
      endif

      weights(d) = (tx(m+1) - x(d)) / (tx(m+1) - tx(m))

      if (weights(d) < 0 .or. weights(d) > 1) then
        call tell_error (tell_runtime_error, &
                         "ndi_calc_weights: extrapolation is not supported", &
                         errstat)
        return
      endif
    enddo
  end subroutine ndi_calc_weights

  !> Compute linear interpolation weight for each lookup table array element
  !! @param[in]  k      Index of a term in the sum of table entries,
  !!                    1 <= k <= 2**N, where N is the number of array dimensions
  !!                    in the lookup table
  !! @param[in] indices Vector of array indices of the table cell
  !!                    containing the interpolation point, \a x.
  !! @param[in] weights Vector of N linear interpolation weights for
  !!                    the interpolation point, \a x.
  !! @param[in] id_k    N multidimensional array indices of the table
  !!                    entry corresponding to term k
  !! @param[in] wt_k    Interpolation weight for term k
  !! @param[inout]  errstat  Error status
  !! @return on success, errstat=0, on error, errstat /=0
  !!
  !! Using this routine, N-dimensional interpolation requires accumulating
  !! a sum of 2^N weighted terms.  For example, a 3-D table interpolation
  !! would look like:
  !!@v+
  !!   num_terms = 2^3
  !!   s = 0.0
  !!   do k = 1, num_terms
  !!      call ndi_calc_term_weight (k, indices, weights, id_k, wt_k)
  !!      s = s + wt_k * value(id_k(1),id_k(2),id_k(3))
  !!   enddo
  !!@v-
  !!
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

  !> Deallocate array of \a ndi_dim objects
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

  !> Allocate an array of \a ndi_dim objects
  !! @param[in]  dims     Array of \a ndi_dim objects to be allocated
  !! @param[in]  dimlens  Array containing the size of each dimension
  !! @param[inout] errstat   Error status code
  !! @return on success, \a errstat=0, on error, \a errstat /= 0
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
