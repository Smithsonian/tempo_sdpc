module tio_module
  use netcdf
  use tell_module
  implicit none

  include '_tempo_dims.inc'
  include '_tempo_grps.inc'
  include '_tempo_vars.inc'

  integer, private, parameter :: &
    i1 = selected_int_kind (2**1), &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4)

  integer :: tiof_get_var_section, tiof_put_att1
  external   tiof_get_var_section, tiof_put_att1

  integer, public, parameter :: &
    tiof_max_var_dims = 7, &
    tiof_max_name_len = 64

  type, public :: tiof_object_type
    integer :: fileid = -1
    integer :: groupid = -1
  end type tiof_object_type

  type, public :: tiof_dim_type
    character (len=tiof_max_name_len) :: name
    integer :: len = 0
    integer :: len_name = 0
    integer :: dimid = -1
    type (tiof_dim_type), pointer :: next => null()
  end type

  type, public :: tiof_dimlist_type
    integer :: num_items = 0
    type (tiof_dim_type), pointer :: head => null(), tail => null()
  end type

  type, public :: tiof_var_type
    character (len=tiof_max_name_len) :: name
    integer :: len_name = 0
    integer :: xtype = -1
    integer :: varid = -1
    character (len=tiof_max_name_len), pointer, dimension(:) :: dim_names => null()
    integer, dimension(tiof_max_var_dims) :: dimids
    integer :: rank = 0
    integer :: deflate_level=0
    logical :: contiguous = .true., shuffle=.false.
    integer, dimension(tiof_max_var_dims) :: chunksizes
    type (tiof_var_type), pointer :: next
  end type

  type, public :: tiof_varlist_type
    integer :: num_items = 0
    type (tiof_var_type), pointer :: head => null(), tail => null()
  end type

  private

  public tiof_open, tiof_close, tiof_inq_group, tiof_inq_dimlen, &
    tiof_get1d_r8, &
    tiof_get1d_r4, tiof_get2d_r4, tiof_get3d_r4, &
    tiof_get2d_i2, tiof_get3d_i2, &
    tiof_get1d_i1, tiof_get2d_i1, &
    tiof_dimlist_append, tiof_dimlist_lookup, tiof_def_dims, &
    tiof_varlist_append, tiof_varlist_lookup, tiof_def_vars

contains

  subroutine tiof_open (file, obj, open_mode, errstat)
    implicit none
    character (len=*), intent(in) :: file
    type (tiof_object_type), intent(out) :: obj
    integer, intent(in) :: open_mode
    integer, intent(inout) :: errstat

    integer :: fileid, status

    if (errstat < 0) return

    status = nf90_open (file, open_mode, fileid)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_open_error, "opening file "//file//" ("//trim(nf90_strerror(status))//")", errstat)
      obj % fileid = -1
      return
    endif

    obj % fileid = fileid
    obj % groupid = fileid

  end subroutine tiof_open

  subroutine tiof_close (obj, errstat)
    implicit none
    type (tiof_object_type), intent(inout) :: obj
    integer, intent(inout) :: errstat

    integer :: status

    if (errstat < 0) return

    if (obj%fileid >= 0) then
      status = nf90_close (obj % fileid)
      if (status /= nf90_noerr) then
        call tell_error (tell_io_error, "closing file ("//trim(nf90_strerror(status))//")", errstat)
      endif
      obj % fileid = -1
      obj % groupid = -1
    endif
  end subroutine tiof_close

  subroutine tiof_inq_group (obj, grpname, errstat)
    implicit none
    type (tiof_object_type), intent(inout) :: obj
    character (len=*), intent(in) :: grpname
    integer, intent(inout) :: errstat

    integer :: status, grp

    if (errstat < 0) return

    status = nf90_inq_ncid (obj % fileid, grpname, grp)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing group "//grpname//" ("//trim(nf90_strerror(status))//")", errstat)
    endif
    obj % groupid = grp
  end subroutine tiof_inq_group

  subroutine tiof_inq_dimlen (obj, name, dimlen, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(out) :: dimlen
    integer, intent(inout) :: errstat

    integer :: status, dimid

    if (errstat < 0) return

    status = nf90_inq_dimid (obj % groupid, name, dimid)
    if (status /= nf90_noerr) then
      call terr_error (terr_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    status = nf90_inquire_dimension (obj % groupid, dimid, len=dimlen)
    if (status /= nf90_noerr) then
      call terr_error (terr_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif
  end subroutine tiof_inq_dimlen

  subroutine tiof_get1d_r8 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=8), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_double, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get1d_r8

  subroutine tiof_get3d_r4 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get3d_r4

  subroutine tiof_get2d_r4 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get2d_r4

  subroutine tiof_get1d_r4 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    real (kind=4), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_float, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get1d_r4

  subroutine tiof_get3d_i2 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i2), dimension (:,:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get3d_i2

  subroutine tiof_get2d_i2 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i2), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_short, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get2d_i2

  subroutine tiof_get2d_i1 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i1), dimension (:,:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_byte, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get2d_i1

  subroutine tiof_get1d_i1 (obj, name, step0, numsteps, array, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: step0, numsteps
    integer (kind=i1), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, step0, numsteps, nf90_byte, array)

    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // name // " from file", errstat)
      return
    endif
  end subroutine tiof_get1d_i1

  subroutine tiof_dimlist_append (list, dim_name, dim_len, errstat)
    implicit none
    type (tiof_dimlist_type), intent(inout) :: list
    character (len=*), intent(in) :: dim_name
    integer, intent(in) :: dim_len
    integer, intent(inout) :: errstat

    ! local
    type (tiof_dim_type), pointer :: item
    integer :: status

    if (errstat < 0) return

    allocate (item, stat=status)
    if (status /= 0) then
      call terr_error (terr_malloc_error, &
                       "tiof_dimlist_append:  allocate failed", errstat)
      return
    endif
    item % next => null()
    item % len  = dim_len
    item % name = adjustl(dim_name)
    item % len_name = len_trim(item % name)

    if (associated(list%head)) then
      list % tail % next => item
    else
      list % head => item
    endif
    list % tail => item
    list % num_items = list % num_items + 1

  end subroutine tiof_dimlist_append

  subroutine tiof_dimlist_lookup (list, num, names, dimids, errstat)
    implicit none
    type (tiof_dimlist_type), intent(in) :: list
    integer, intent(in) :: num
    character (len=*), target, dimension(:), intent(in) :: names
    integer, dimension(:), intent(out) :: dimids
    integer, intent(inout) :: errstat

    type (tiof_dim_type), pointer :: item => null()
    character (len=tiof_max_name_len) :: name_i
    integer :: i, len_i

    if (errstat < 0) return

    if (.not.associated(list%head)) then
      call terr_error (terr_invalid_parm, &
                       "tiof_dimlist_lookup: null dimension list", &
                       errstat)
      return
    endif

    do i=1, num

      item => list % head
      name_i = adjustl(names(i))
      len_i = len_trim(name_i)

      dimids(i) = -1
      search: do while (associated(item))
        if (item%len_name /= len_i) cycle
        if (item % name(1:item%len_name) == name_i(1:len_i)) then
          dimids(i) = item % dimid
          exit search
        endif
        item => item % next
      enddo search

    enddo

  end subroutine tiof_dimlist_lookup

  subroutine tiof_def_dims (obj, list, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: list
    integer, intent(inout) :: errstat

    type (tiof_dim_type), pointer :: item => null()
    integer :: status

    if (errstat < 0) return

    if (.not.associated (list%head)) then
      call terr_error (terr_invalid_parm, &
                       "tiof_def_dims: null dimension list", &
                       errstat)
      return
    endif

    item => list % head
    do while (associated(item))
      status = nf90_def_dim (obj % groupid, item % name, item % len, &
                             item % dimid)
      if (status /= nf90_noerr) then
        call terr_error (terr_io_error, "defining dimension " &
                         //item % name(1:item%len_name)//" (" // &
                         trim(nf90_strerror(status))//") ", &
                         errstat)
        return
      endif
      item => item % next
    enddo

  end subroutine tiof_def_dims

  subroutine tiof_varlist_append (list, var_name, xtype, dimids, errstat, &
                                  shuffle, deflate_level, contiguous, chunksizes)
    implicit none
    type (tiof_varlist_type), intent(inout) :: list
    character (len=*), intent(in) :: var_name
    integer, intent(in) :: xtype
    integer, dimension(:), intent(in) :: dimids
    integer, intent(inout) :: errstat
    integer, optional, intent(in) :: deflate_level
    logical, optional, intent(in) :: contiguous, shuffle
    integer, optional, dimension(:), intent(in) :: chunksizes

    ! local
    type (tiof_var_type), pointer :: item
    integer :: status

    if (errstat < 0) return

    allocate (item, stat=status)
    if (status /= 0) then
      call terr_error (terr_malloc_error, &
                       "tiof_varlist_append:  allocate failed", errstat)
      return
    endif
    item % next => null()
    item % name = adjustl(var_name)
    item % len_name = len_trim(item % name)
    item % xtype = xtype
    item % rank = size(dimids)
    item%dimids(1:item%rank) = dimids(1:item%rank)

    if (present(shuffle)) then
      item % shuffle = shuffle
    endif

    if (present(deflate_level)) then
      item % deflate_level = deflate_level
    endif

    if (present(contiguous)) then
      item % contiguous = contiguous
      if (present(chunksizes)) then
        item % chunksizes(1:item%rank) = chunksizes(1:item%rank)
      endif
    endif

    if (associated(list%head)) then
      list % tail % next => item
    else
      list % head => item
    endif
    list % tail => item
    list % num_items = list % num_items + 1

  end subroutine tiof_varlist_append

  subroutine tiof_varlist_lookup (list, name, varid, errstat)
    implicit none
    type (tiof_varlist_type), intent(in) :: list
    character (len=*), intent(in) :: name
    integer, intent(out) :: varid
    integer, intent(inout) :: errstat

    type (tiof_var_type), pointer :: item => null()
    integer :: len_name

    varid = -1

    if (errstat < 0) return

    if (.not.associated(list%head)) then
      call terr_error (terr_invalid_parm, &
                       "tiof_dimlist_lookup: null dimension list", &
                       errstat)
      return
    endif

    len_name = len_trim(name)

    item => list % head

    do while (associated(item))
      if (item % len_name /= len_name) cycle
      if (item % name(1:item%len_name) == name (1:len_name)) then
        varid = item % varid
        return
      endif
      item => item % next
    enddo

  end subroutine tiof_varlist_lookup

  subroutine tiof_def_vars (obj, list, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    type (tiof_varlist_type), intent(in) :: list
    integer, intent(inout) :: errstat

    type (tiof_var_type), pointer :: item => null()
    integer :: status

    if (errstat < 0) return

    if (.not.associated (list%head)) then
      call terr_error (terr_invalid_parm, &
                       "tiof_def_vars: null variable list", &
                       errstat)
      return
    endif

    item => list % head
    do while (associated(item))
      status = nf90_def_var (obj % groupid, item % name, item % xtype, &
                             item % dimids(1:item%rank), &
                             item % varid, &
                             contiguous = item % contiguous, &
                             chunksizes = item % chunksizes(1:item%rank), &
                             deflate_level = item % deflate_level, &
                             shuffle = item % shuffle)
      if (status /= nf90_noerr) then
        call terr_error (terr_io_error, "defining variable " &
                         //item % name(1:item%len_name)//" (" // &
                         trim(nf90_strerror(status))//") ", &
                         errstat)
        return
      endif
      item => item % next
    enddo

  end subroutine tiof_def_vars

end module tio_module
