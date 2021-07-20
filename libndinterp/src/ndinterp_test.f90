program ndinterp_test
  use ndinterp_module
  implicit none
  integer, allocatable, dimension(:) :: seed
  integer :: errstat, n, i, pid

  call random_seed (size = n)
  allocate (seed(n))
  ! For a regression test, we want repeatability
  pid = 31415 ! getpid()
  do i = 1, n
    seed(i) = pid
  enddo
  call random_seed (put = seed)

  errstat = 0

  call test2d ((/5,7/), errstat)
  if (errstat /= 0) stop 1

  !call test6d ((/5,7,3,4,9, 2,10/), errstat)
  call test6d ((/601, 10, 10, 7, 9, 11/), errstat)
  if (errstat /= 0) stop 2

contains

  function ftest(x)
    implicit none
    real (kind=8), dimension(:), intent(in) :: x
    real (kind=8) :: ftest
    ftest = sum(x)
  end function ftest

  subroutine test2d (dimlens, errstat)
    implicit none
    integer, dimension(:), intent(in) :: dimlens
    integer, intent(inout) :: errstat

    type (ndi_dim), dimension(2) :: dims
    real, allocatable, dimension(:,:) :: tbl
    real (kind=8), dimension(size(dimlens)) :: x
    integer, dimension(size(dimlens)) :: indices
    integer :: i, d, i1, i2, nd
    real :: val

    if (errstat /= 0) return

    call ndi_dims_alloc (dims, dimlens, errstat)
    if (errstat /= 0) return

    allocate (tbl(dimlens(1),dimlens(2)))

    nd = size(dimlens)

    do d = 1, nd
      do i = 1, dimlens(d)
        dims(d) % x(i) = (i - 1) * 1.0 / (dimlens(d) - 1)
      enddo
    enddo

    do i2 = 1, dimlens(2)
      x(2) = dims(2) % x(i2)
      do i1 = 1, dimlens(1)
        x(1) = dims(1) % x(i1)
        tbl(i1,i2) = real(ftest (x), kind=4)
      enddo
    enddo

    call random_number (x)
    call ndi_table_interp (dims, tbl, x, val, errstat, output_indices=indices)
    if (errstat /= 0) stop 2

    write (*,'(a,6(e10.5,1x))')'x = ',x
    write (*,'(a,6(i3,1x))')' indices = ',indices
    write (*,'(a,e10.5)')'          val = ',val
    write (*,'(a,e10.5)')'expected f(x) = ',ftest(x)

  end subroutine test2d

  subroutine test6d (dimlens, errstat)
    implicit none
    integer, parameter :: nd=6
    integer, dimension(:), intent(in) :: dimlens
    integer, intent(inout) :: errstat

    type (ndi_dim), dimension(nd) :: dims
    real, allocatable, dimension(:,:,:,:,:,:) :: tbl
    real (kind=8), dimension(size(dimlens)) :: x
    integer :: i, d, i1, i2, i3, i4, i5, i6
    real :: val
    !real :: start, finish
    real (kind=8) :: sum_sq, fv, favg
    integer :: num_trials

    if (errstat /= 0) return

    call ndi_dims_alloc (dims, dimlens, errstat)
    if (errstat /= 0) return

    allocate (tbl(dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5),dimlens(6)))

    do d = 1, nd
      do i = 1, dimlens(d)
        dims(d) % x(i) = (i - 1) * 1.0 / (dimlens(d) - 1)
      enddo
    enddo

    do i6=1, dimlens(6)
      x(6) = dims(6) % x(i6)
      do i5=1, dimlens(5)
        x(5) = dims(5) % x(i5)
        do i4=1, dimlens(4)
          x(4) = dims(4) % x(i4)
          do i3=1, dimlens(3)
            x(3) = dims(3) % x(i3)
            do i2=1, dimlens(2)
              x(2) = dims(2) % x(i2)
              do i1=1, dimlens(1)
                x(1) = dims(1) % x(i1)
                tbl (i1,i2,i3,i4,i5,i6) = real (ftest (x), kind=4)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    call random_number (x)
    call ndi_table_interp (dims, tbl, x, val, errstat)
    if (errstat /= 0) stop 2

    write (*,'(a,6(e10.5,1x))')'x = ',x
    write (*,'(a,e10.5)')'          val = ',val
    write (*,'(a,e10.5)')'expected f(x) = ',ftest(x)

    num_trials = 1000
    write (*,*)'Starting ',num_trials,' trials'

    sum_sq = 0.0
    favg = 0.0

    !call cpu_time(start)
    do i = 1,num_trials
      call random_number (x)
      call ndi_table_interp (dims, tbl, x, val, errstat)
      if (errstat /= 0) stop 3
      fv = ftest(x)
      sum_sq = sum_sq + (val-fv)**2
      favg = favg + fv
    enddo
    !call cpu_time(finish)
    !write (*,*)'Elapsed time=',finish-start,' sec'
    write (*,'(a,e10.3)')'RMS error = ',sqrt(sum_sq/num_trials)
    write (*,'(a,e10.5)')'F mean = ',favg/num_trials

  end subroutine test6d

end program ndinterp_test
