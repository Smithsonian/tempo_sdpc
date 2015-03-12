!> Fortran interface module
!! @sa get_put_code.inc
!!
!! @details
!! **Error handling:**
!!
!! All member functions take an integer error code, \a errstat.
!! When the input \a errstat is negative, the function returns
!! immediately, with \a errstat unchanged. On return, \a errstat<0
!! indicates that an error occured. Error messages and global error
!! status are handled using libtell. Additional information on the
!! type of error that occurred may be available from the global error
!! status, via \a tell_get_error and \a tell_copy_strerror.
!
module tio_module
  use netcdf
  use tell_module
  implicit none
  private

  integer, public, parameter :: &
    tiof_max_var_dims = 7, &
    tiof_max_name_len = 64, &
    tiof_max_att_len = 256

  include '_tempo_dims.inc'
  include '_tempo_grps.inc'
  include '_tempo_vars.inc'

  integer, private, parameter :: &
    i1 = selected_int_kind (2**1), &
    i2 = selected_int_kind (2**2), &
    i4 = selected_int_kind (2**3), &
    i8 = selected_int_kind (2**4)

  !> File object
  type, public :: tiof_file_type
    integer :: fileid
    integer :: groupid
    integer, private :: grp_stack(10), grp_stack_depth
  end type tiof_file_type

  !> Dimension object
  type, public :: tiof_dim_type
    character (len=tiof_max_name_len) :: name  !< dimension name
    integer :: len_name = 0
    integer :: len = 0
    integer :: dimid = -huge(1)   !< dimension id assigned by netCDF
    type (tiof_dim_type), private, pointer :: next => null()
  end type

  !> Dimension object list
  type, public :: tiof_dimlist_type
    integer :: num_items = 0
    type (tiof_dim_type), private, pointer :: head => null(), tail => null()
  end type

  !> Attribute object
  type, public :: tiof_att_type
    character (len=tiof_max_name_len) :: name  !< attribute name
    integer :: len_name = 0
    integer :: xtype = -1       !< attribute value data type
    character (len=tiof_max_att_len) :: att_text
    integer (kind=i4), allocatable, dimension(:) :: att_i4
    integer (kind=i8), allocatable, dimension(:) :: att_i8
    real (kind=4), allocatable, dimension(:) :: att_r4
    real (kind=8), allocatable, dimension(:) :: att_r8
    type (tiof_att_type), private, pointer :: next => null()
  end type

  !> Attribute object list
  type, public :: tiof_attlist_type
    integer :: num_items = 0
    type (tiof_att_type), private, pointer :: head => null(), tail => null()
  end type

  !> Variable object
  type, public :: tiof_var_type
    character (len=tiof_max_name_len) :: name  !< variable name
    integer :: len_name = 0
    integer :: xtype = -1          !< variable external data type
    integer :: varid = -huge(1)    !< variable id assigned by netCDF
    integer :: rank = 0            !< number of dimensions
    integer, dimension(tiof_max_var_dims) :: dimids  !< ordered list of dimension ids
    character (len=tiof_max_att_len) :: comment  !< attribute: comment
    character (len=tiof_max_att_len) :: units    !< attribute: units
    real (kind=8), dimension(2) :: valid_range = [0.0, 0.0]  !< attribute: valid_min, valid_max
    integer :: deflate_level = 0   !< attribute: compression deflate level
    logical :: shuffle = .false.   !< attribute: compress with shuffle?
    logical :: contiguous = .true. !< attribute: use contiguous storage?
    integer, dimension(tiof_max_var_dims) :: chunksizes = 0  !< attribute: chunk sizes
    integer :: no_fill = 0         !< attribute: if non-zero, turn off auto-fill
    real (kind=8) :: fillvalue     !< attribute: fill value
    logical, private :: have_comment=.false., have_units=.false., &
      have_valid_range=.false., have_fillvalue = .false.
    type (tiof_attlist_type), private, pointer :: attlist => null()
    type (tiof_var_type), private, pointer :: next => null()
  end type

  !> Variable object list
  type, public :: tiof_varlist_type
    integer :: num_items = 0
    type (tiof_var_type), private, pointer :: head => null(), tail => null()
  end type

  integer :: tiof_get_var_section, tiof_put_var_section, tio_f_put_git_hash, &
    tio_f_def_grp
  external   tiof_get_var_section, tiof_put_var_section, tio_f_put_git_hash, &
    tio_f_def_grp
  public     tiof_put_var_section, tiof_get_var_section, tio_f_put_git_hash, &
    tio_f_def_grp

  public tiof_create, tiof_open, tiof_close, &
    tiof_put_git_commit_hash, &
    tiof_def_group, tiof_push_group, tiof_pop_group, tiof_inq_group, &
    tiof_inq_dimlen, &
    tiof_dimlist_append, tiof_dimlist_free, tiof_def_dims, tiof_dimlist_lookup, &
    tiof_varlist_append, tiof_varlist_free, tiof_def_vars, tiof_varlist_lookup, &
    tiof_attlist_append, tiof_attlist_free, tiof_def_atts

  public tiof_put1d_text, tiof_get1d_text
  public tiof_put1d_string, tiof_get1d_string
  include 'get_put_decl.inc'

contains

  include 'get_put_code.inc'

  subroutine push_group (obj, groupid, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    integer, intent(in) :: groupid
    integer, intent(inout) :: errstat

    if (errstat < 0) return

    if (obj % grp_stack_depth == size (obj % grp_stack)) then
      call tell_error (tell_internal_error, "push_group:  stack depth exceeded", errstat)
      return
    endif

    obj % grp_stack_depth = obj % grp_stack_depth + 1
    obj % grp_stack (obj % grp_stack_depth) = obj % groupid
    obj % groupid = groupid

  end subroutine push_group

  subroutine pop_group (obj, groupid, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    integer, intent(out) :: groupid
    integer, intent(inout) :: errstat

    if (errstat < 0) return

    if (obj % grp_stack_depth == 0) then
      groupid = -1
      return
    endif

    groupid = obj % groupid
    obj % groupid = obj % grp_stack (obj % grp_stack_depth)
    obj % grp_stack_depth = obj % grp_stack_depth - 1

  end subroutine pop_group

  subroutine tiof_put_git_commit_hash (obj, errstat, name)
    use iso_c_binding, only : c_null_ptr, c_null_char
    implicit none
    type (tiof_file_type), intent(in) :: obj
    integer, intent(inout):: errstat
    character (len=*), optional, intent(in) :: name

    if (errstat < 0) return

    if (present(name)) then
      errstat = tio_f_put_git_hash (obj % groupid, trim(adjustl(name))//c_null_char)
    else
      errstat = tio_f_put_git_hash (obj % groupid, c_null_ptr)
    endif

  end subroutine

  !> write a 1d array of strings as a 2D array of characters
  subroutine tiof_put1d_text (obj, name, start, edge, array, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: start, edge
    character (len=*), dimension (:), intent(in) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_put_var_section (obj % groupid, name, [start,0], [edge, -1], &
                                nf90_char, array)
    if (err < 0) then
      call tell_error (tell_io_write_error, "Unable to write " // &
                       trim(name) // " to file", errstat)
      return
    endif
  end subroutine tiof_put1d_text

  !> read a 1d array of strings stored as a 2D array of characters
  subroutine tiof_get1d_text (obj, name, start, edge, array, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: start, edge
    character (len=*), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err

    if (errstat < 0) return

    err = tiof_get_var_section (obj % groupid, name, [start,0], [edge, -1], &
                                nf90_char, array)
    if (err < 0) then
      call tell_error (tell_io_read_error, "Unable to read " // &
                       trim(name) // " from file", errstat)
      return
    endif
  end subroutine tiof_get1d_text

  !> write a 1d array of strings as null-terminated strings
  subroutine tiof_put1d_string (obj, name, start, edge, array, errstat)
    use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
    implicit none
    type (tiof_file_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: start, edge
    character (len=*), dimension (:), intent(in) :: array
    integer, intent(inout) :: errstat

    integer :: err, i, num
    character(len=len(array)), dimension(:), allocatable, target :: copy
    type(c_ptr), dimension(:), allocatable :: ptrs

    if (errstat < 0) return

    num = size(array)
    allocate (ptrs(num), copy(num), stat=err)
    if (err /= 0) then
      errstat = tell_malloc_error
      call tell_set_error (tell_malloc_error)
      return
    endif

    do i=1,num
      copy(i) = trim(adjustl(array(i)))//c_null_char
      ptrs(i) = c_loc(copy(i))
    enddo

    err = tiof_put_var_section (obj % groupid, name, [start], [edge], &
                                nf90_string, ptrs)
    if (err < 0) then
      call tell_error (tell_io_write_error, "Unable to write " // &
                       trim(name) // " to file", errstat)
      return
    endif
  end subroutine tiof_put1d_string

  !> read a 1d array of strings stored as null-terminated strings
  subroutine tiof_get1d_string (obj, name, start, edge, array, errstat)
    use iso_c_binding, only: c_ptr, c_null_ptr, c_f_pointer, c_null_char
    implicit none
    type (tiof_file_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(in) :: start, edge
    character (len=*), dimension (:), intent(out) :: array
    integer, intent(inout) :: errstat

    integer :: err, num, i, k, alen
    type(c_ptr), dimension(size(array)) :: ptrs
    character (len=size(array)), pointer :: fptr

    if (errstat < 0) return

    fptr => null()
    ptrs(:) = c_null_ptr
    err = tiof_get_var_section (obj % groupid, name, [start], [edge], &
                                nf90_string, ptrs)
    if (err < 0) then
      call tell_error (tell_io_write_error, "Unable to read " // &
                       trim(name) // " from file", errstat)
      return
    endif

    num = size(array)
    alen = len(array)
    do i=1,num
      call c_f_pointer (ptrs(i), fptr)
      if (associated(fptr)) then
        array(i)(1:alen) = ' '
        k = min(index(fptr, c_null_char), alen)
        array(i)(1:k) = fptr(1:k)
        ! FIXME?: Netcdf allocated memory for these strings, and
        ! that memory must be freed here.  Unfortunately, netcdf
        ! provides no fortran interface for a function to free that
        ! memory. I could write a fortran wrapper for nc_free_string(),
        ! but deallocate seems to work just as well.
        deallocate(fptr)
      endif
    enddo

  end subroutine tiof_get1d_string

  !> Create a new netcdf4/HDF5 file
  subroutine tiof_create (obj, file, create_mode, errstat)
    implicit none
    type (tiof_file_type), intent(out) :: obj
    character (len=*), intent(in) :: file
    integer, intent(in) :: create_mode
    integer, intent(inout) :: errstat

    integer :: fileid, status

    if (errstat < 0) return

    status = nf90_create (file, ior (create_mode, nf90_netcdf4), fileid)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_open_error, "creating file "//file// &
                       " ("//trim(nf90_strerror(status))//")", errstat)
      obj % fileid = -1
      return
    endif

    obj % fileid = fileid
    obj % groupid = fileid
    obj % grp_stack_depth = 0

  end subroutine tiof_create

  !> Open an existing netcdf4/HDF5 file
  subroutine tiof_open (file, obj, open_mode, errstat)
    implicit none
    character (len=*), intent(in) :: file
    type (tiof_file_type), intent(out) :: obj
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
    obj % grp_stack_depth = 0

  end subroutine tiof_open

  !> Close a netcdf4/HDF5 file
  subroutine tiof_close (obj, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
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
      obj % grp_stack_depth = 0
    endif

  end subroutine tiof_close

  !> Define a new group
  subroutine tiof_def_group (obj, grpname, errstat, groupid)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    character (len=*), intent(in) :: grpname
    integer, intent(inout) :: errstat
    integer, optional, intent(out) :: groupid

    integer :: status, grp

    if (errstat < 0) return

    status = tio_f_def_grp (obj % groupid, grpname, grp)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_write_error, "creating group "//grpname//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    if (present(groupid)) groupid = grp

  end subroutine tiof_def_group

  !> Associate an open file object with a specific group
  subroutine tiof_inq_group (obj, grpname, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    character (len=*), intent(in) :: grpname
    integer, intent(inout) :: errstat

    integer :: status, grp

    if (errstat < 0) return

    if (len(grpname) == 0) then
      status = nf90_inq_ncid (obj % groupid, grpname, grp)
    else if (grpname(1:1) == '/') then
      status = nf90_inq_grp_full_ncid (obj % fileid, grpname, grp)
    else
      status = nf90_inq_ncid (obj % groupid, grpname, grp)
    endif

    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing group "//grpname//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    obj % groupid = grp
    obj % grp_stack_depth = 0  ! reset the stack

  end subroutine tiof_inq_group

  !> Associate an open file object with a specific group, saving the current group
  subroutine tiof_push_group (obj, grpname, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    character (len=*), intent(in) :: grpname
    integer, intent(inout) :: errstat

    integer :: status, grp

    if (errstat < 0) return

    if (len(grpname) == 0) then
      status = nf90_inq_ncid (obj % groupid, grpname, grp)
    else if (grpname(1:1) == '/') then
      status = nf90_inq_grp_full_ncid (obj % fileid, grpname, grp)
    else
      status = nf90_inq_ncid (obj % groupid, grpname, grp)
    endif

    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing group "//grpname//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    call push_group (obj, grp, errstat)

  end subroutine tiof_push_group

  !> Restore the groupid saved by the most recent call to tiof_push_group
  subroutine tiof_pop_group (obj, errstat)
    implicit none
    type (tiof_file_type), intent(inout) :: obj
    integer, intent(inout) :: errstat

    integer :: prev_grp

    !if (errstat < 0) return

    call pop_group (obj, prev_grp, errstat)

  end subroutine tiof_pop_group

  !> Inquire the size of a dimension
  subroutine tiof_inq_dimlen (obj, name, dimlen, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    character (len=*), intent(in) :: name
    integer, intent(out) :: dimlen
    integer, intent(inout) :: errstat

    integer :: status, dimid

    if (errstat < 0) return

    status = nf90_inq_dimid (obj % groupid, name, dimid)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    status = nf90_inquire_dimension (obj % groupid, dimid, len=dimlen)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif
  end subroutine tiof_inq_dimlen

  !> Append a new dimension object to a dimension list
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
      call tell_error (tell_malloc_error, &
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

  !> Free list of dimension objects (tiof_dimlist_type)
  subroutine tiof_dimlist_free (dimlist)
    implicit none
    type (tiof_dimlist_type), intent(inout) :: dimlist
    type (tiof_dim_type), pointer :: item

    do while (associated(dimlist % head))
      item => dimlist % head % next
      deallocate (dimlist % head)
      dimlist % head => item
    enddo
    dimlist % tail => null()
    dimlist % head => null()

  end subroutine tiof_dimlist_free

  !> Retrieve information about an array of dimension objects
  !! from a dimension list
  subroutine tiof_dimlist_lookup (list, names, dimids, errstat, &
                                  dimsizes)
    implicit none
    type (tiof_dimlist_type), intent(in) :: list
    character (len=*), dimension(:), intent(in) :: names
    integer, dimension(:), intent(out) :: dimids
    integer, intent(inout) :: errstat
    integer, optional, intent(inout), dimension(size(dimids)) :: dimsizes

    type (tiof_dim_type), pointer :: item
    character (len=tiof_max_name_len) :: name_i
    integer :: i, num, len_i

    if (errstat < 0) return

    if (.not.associated(list%head)) then
      call tell_error (tell_invalid_parm, &
                       "tiof_dimlist_lookup: null dimension list", &
                       errstat)
      return
    endif

    num = size(names)

    if (present (dimsizes)) then
      dimsizes(1:num) = -1
    endif

    do i=1, num

      item => list % head
      name_i = adjustl(names(i))
      len_i = len_trim(name_i)

      dimids(i) = -1
      search: do while (associated(item))
        if (item%len_name == len_i &
            .and. item % name(1:item%len_name) == name_i(1:len_i)) then
          dimids(i) = item % dimid
          if (present(dimsizes)) then
            dimsizes(i) = item % len
          endif
          exit search
        endif
        item => item % next
      enddo search

    enddo

  end subroutine tiof_dimlist_lookup

  !> Write a dimension list to a file
  !! @details
  !! As each dimension object (\a tiof_dim_type) is written
  !! to the file, the assigned id number is saved in the
  !! \a dimid field.
  !! @sa tiof_dimlist_lookup
  subroutine tiof_def_dims (obj, list, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: list
    integer, intent(inout) :: errstat

    type (tiof_dim_type), pointer :: item
    integer :: status

    if (errstat < 0) return

    if (.not.associated (list%head)) then
      call tell_error (tell_invalid_parm, &
                       "tiof_def_dims: null dimension list", &
                       errstat)
      return
    endif

    item => list % head
    do while (associated(item))
      status = nf90_def_dim (obj % groupid, item % name, item % len, &
                             item % dimid)
      if (status /= nf90_noerr) then
        call tell_error (tell_io_error, "defining dimension " &
                         //item % name(1:item%len_name)//" (" // &
                         trim(nf90_strerror(status))//") ", &
                         errstat)
        return
      endif
      item => item % next
    enddo

  end subroutine tiof_def_dims

  !> Append a new attribute object to an attribute list
  subroutine tiof_attlist_append (list, errstat, att_name, &
                                  att_i4, att_i8, att_r4, att_r8, att_text)
    implicit none
    type (tiof_attlist_type), intent(inout) :: list
    character (len=*), intent(in) :: att_name
    integer, intent(inout) :: errstat
    integer (kind=i4), optional, dimension(:) :: att_i4
    integer (kind=i8), optional, dimension(:) :: att_i8
    real (kind=4), optional, dimension(:) :: att_r4
    real (kind=8), optional, dimension(:) :: att_r8
    character (len=*), optional :: att_text

    ! local
    type (tiof_att_type), pointer :: item
    integer :: status

    if (errstat < 0) return

    allocate (item, stat=status)
    if (status /= 0) then
      call tell_error (tell_malloc_error, &
                       "tiof_attlist_append:  allocate failed", errstat)
      return
    endif
    item % next => null()
    item % name = adjustl(att_name)
    item % len_name = len_trim(item % name)

    ! each struct may contain only a single data type

    if (present(att_i4)) then
      allocate (item % att_i4(size(att_i4)), stat=status)
      if (status /= 0) then
        call tell_error (tell_malloc_error, &
                         "tiof_attlist_append:  allocate failed", errstat)
      endif
      item % att_i4(:) = att_i4(:)
      item % xtype = nf90_int
    else if (present(att_i8)) then
      allocate (item % att_i8(size(att_i8)), stat=status)
      if (status /= 0) then
        call tell_error (tell_malloc_error, &
                         "tiof_attlist_append:  allocate failed", errstat)
      endif
      item % att_i8(:) = att_i8(:)
      item % xtype = nf90_int64
    else if (present(att_r4)) then
      allocate (item % att_r4(size(att_r4)), stat=status)
      if (status /= 0) then
        call tell_error (tell_malloc_error, &
                         "tiof_attlist_append:  allocate failed", errstat)
      endif
      item % att_r4(:) = att_r4(:)
      item % xtype = nf90_float
    else if (present(att_r8)) then
      allocate (item % att_r8(size(att_r8)), stat=status)
      if (status /= 0) then
        call tell_error (tell_malloc_error, &
                         "tiof_attlist_append:  allocate failed", errstat)
      endif
      item % att_r8(:) = att_r8(:)
      item % xtype = nf90_double
    else if (present(att_text)) then
      item % att_text = adjustl(att_text)
      item % xtype = nf90_char
    endif

    if (associated(list%head)) then
      list % tail % next => item
    else
      list % head => item
    endif
    list % tail => item
    list % num_items = list % num_items + 1

  end subroutine tiof_attlist_append

  !> Free list of attribute objects (tiof_attlist_type)
  subroutine tiof_attlist_free (attlist)
    implicit none
    type (tiof_attlist_type), intent(inout) :: attlist
    type (tiof_att_type), pointer :: item

    do while (associated(attlist % head))
      item => attlist % head % next
      deallocate (attlist % head)
      attlist % head => item
    enddo
    attlist % tail => null()
    attlist % head => null()

  end subroutine tiof_attlist_free

  !> Write an attribute list to a file
  subroutine tiof_def_atts (obj, list, varid, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    type (tiof_attlist_type), intent(in) :: list
    integer, intent(in) :: varid
    integer, intent(inout) :: errstat

    type (tiof_att_type), pointer :: item
    integer :: status

    if (errstat < 0) return

    if (.not.associated (list%head)) then
      call tell_error (tell_invalid_parm, &
                       "tiof_def_atts: null attribute list", &
                       errstat)
      return
    endif

    item => list % head
    do while (associated(item))

      select case (item % xtype)
        case (nf90_char)
          status = nf90_put_att (obj % groupid, varid, item % name, item % att_text)
        case (nf90_double)
          status = nf90_put_att (obj % groupid, varid, item % name, item % att_r8)
        case (nf90_float)
          status = nf90_put_att (obj % groupid, varid, item % name, item % att_r4)
        case (nf90_int)
          status = nf90_put_att (obj % groupid, varid, item % name, item % att_i4)
        case (nf90_int64)
          status = nf90_put_att (obj % groupid, varid, item % name, item % att_i8)
        case default
          call tell_error (tell_invalid_parm, &
                           "tiof_def_atts: unsupported attribute type", &
                           errstat)
          return
      end select

      if (status /= nf90_noerr) then
        call tell_error (tell_io_error, "defining attribute " &
                         //item % name(1:item%len_name)//" (" // &
                         trim(nf90_strerror(status))//") ", &
                         errstat)
        return
      endif

      item => item % next
    enddo

  end subroutine tiof_def_atts

  !> Append a new variable object to a variable list
  subroutine tiof_varlist_append (list, errstat, var_name, xtype, dimids, &
                                  shuffle, deflate_level, contiguous, chunksizes, &
                                  comment, units, valid_range, &
                                  no_fill, fillvalue, &
                                  attlist)
    implicit none
    type (tiof_varlist_type), intent(inout) :: list
    integer, intent(inout) :: errstat
    character (len=*), intent(in) :: var_name
    integer, intent(in) :: xtype
    integer, optional, dimension(:), intent(in) :: dimids
    integer, optional, intent(in) :: deflate_level
    logical, optional, intent(in) :: contiguous, shuffle
    integer, optional, dimension(:), intent(in) :: chunksizes
    character (len=*), optional, intent(in) :: comment, units
    real (kind=8), optional, dimension(2), intent(in) :: valid_range
    integer, optional, intent(in) :: no_fill
    real (kind=8), optional, intent(in) :: fillvalue
    type (tiof_attlist_type), optional, target, intent(in) :: attlist

    ! local
    type (tiof_var_type), pointer :: item
    integer :: status

    if (errstat < 0) return

    allocate (item, stat=status)
    if (status /= 0) then
      call tell_error (tell_malloc_error, &
                       "tiof_varlist_append:  allocate failed", errstat)
      return
    endif
    item % next => null()
    item % name = adjustl(var_name)
    item % len_name = len_trim(item % name)
    item % xtype = xtype

    if (.not.present(dimids)) then
      item % rank = 0
    else
      item % rank = size(dimids)
      item % dimids(1:item%rank) = dimids(1:item%rank)

      if (present(contiguous)) then
        item % contiguous = contiguous
      endif

      if (present(chunksizes)) then
        item % chunksizes(1:item%rank) = chunksizes(1:item%rank)
        if (any (chunksizes(1:item%rank) > 1)) item % contiguous = .false.
      endif

      if (present(shuffle)) then
        item % shuffle = shuffle
        if (shuffle) item % contiguous = .false.
      endif

      if (present(deflate_level)) then
        item % deflate_level = deflate_level
        if (deflate_level > 0) item % contiguous = .false.
      endif
    endif

    if (present(comment)) then
      item % have_comment = .true.
      item % comment = adjustl(comment)
    endif

    if (present(units)) then
      item % have_units = .true.
      item % units = adjustl(units)
    endif

    if (present(valid_range)) then
      item % have_valid_range = .true.
      item % valid_range = valid_range(1:2)
    endif

    if (present(no_fill)) then
      item % no_fill = no_fill
      item % have_fillvalue = (no_fill /= 0)
    endif

    select case (xtype)
      case (nf90_double)
        item % fillvalue = nf90_fill_double
      case (nf90_float)
        item % fillvalue = nf90_fill_float
      case (nf90_int)
        item % fillvalue = nf90_fill_int
      case (nf90_short)
        item % fillvalue = nf90_fill_short
      case (nf90_byte)
        item % fillvalue = nf90_fill_byte
      case (nf90_char)
        item % fillvalue = nf90_fill_byte
      case default
        item % fillvalue = nf90_fill_double
    end select

    if (present(fillvalue)) then
      item % have_fillvalue = .true.
      item % fillvalue = fillvalue
    endif

    if (present(attlist)) then
      item % attlist => attlist
    endif

    if (associated(list%head)) then
      list % tail % next => item
    else
      list % head => item
    endif
    list % tail => item
    list % num_items = list % num_items + 1

  end subroutine tiof_varlist_append

  !> Free list of variable objects (tiof_varlist_type)
  subroutine tiof_varlist_free (varlist)
    implicit none
    type (tiof_varlist_type), intent(inout) :: varlist
    type (tiof_var_type), pointer :: item

    do while (associated(varlist % head))
      item => varlist % head % next
      deallocate (varlist % head)
      varlist % head => item
    enddo
    varlist % tail => null()
    varlist % head => null()

  end subroutine tiof_varlist_free

  !> Retrieve a variable object from a given variable list
  subroutine tiof_varlist_lookup (list, name, var_ptr, errstat)
    implicit none
    type (tiof_varlist_type), intent(in) :: list
    character (len=*), intent(in) :: name
    type (tiof_var_type), pointer, intent(out) :: var_ptr
    integer, intent(inout) :: errstat

    type (tiof_var_type), pointer :: item
    integer :: len_name

    var_ptr => null()

    if (.not.associated(list%head)) then
      call tell_error (tell_invalid_parm, &
                       "tiof_dimlist_lookup: null dimension list", &
                       errstat)
      return
    endif

    len_name = len_trim(name)

    item => list % head

    do while (associated(item))
      if (item % len_name == len_name &
          .and. item % name(1:item%len_name) == name (1:len_name)) then
        var_ptr => item
        return
      endif
      item => item % next
    enddo

  end subroutine tiof_varlist_lookup

  !> Write a variable list to a file
  !! @details
  !! As each variable object (\a tiof_var_type) is written
  !! to the file, the assigned id number is saved in the
  !! \a varid field.
  !! @sa tiof_varlist_lookup
  subroutine tiof_def_vars (obj, list, errstat)
    implicit none
    type (tiof_file_type), intent(in) :: obj
    type (tiof_varlist_type), intent(in) :: list
    integer, intent(inout) :: errstat

    type (tiof_var_type), pointer :: item
    integer :: status

    if (errstat < 0) return

    if (.not.associated (list%head)) then
      call tell_error (tell_invalid_parm, &
                       "tiof_def_vars: null variable list", &
                       errstat)
      return
    endif

    item => list % head
    do while (associated(item))
      if (item % rank == 0) then
        status = nf90_def_var (obj % groupid, item % name, &
                               item % xtype, item % varid)
      else if (item % contiguous) then
        status = nf90_def_var (obj % groupid, item % name, &
                               item % xtype, &
                               item % dimids(1:item%rank), &
                               item % varid)
      else if (any(item%chunksizes(1:item%rank) > 1)) then
        status = nf90_def_var (obj % groupid, item % name, &
                               item % xtype, &
                               item % dimids(1:item%rank), &
                               item % varid, &
                               chunksizes = item % chunksizes(1:item%rank), &
                               deflate_level = item % deflate_level, &
                               shuffle = item % shuffle)
      else
        status = nf90_def_var (obj % groupid, item % name, &
                               item % xtype, &
                               item % dimids(1:item%rank), &
                               item % varid, &
                               deflate_level = item % deflate_level, &
                               shuffle = item % shuffle)
      endif

      if (status /= nf90_noerr) then
        call tell_error (tell_io_error, "defining variable " &
                         //item % name(1:item%len_name)//" (" // &
                         trim(nf90_strerror(status))//") ", &
                         errstat)
        return
      endif

      if (item % have_fillvalue) then
        select case (item % xtype)
          case (nf90_double)
            status = nf90_def_var_fill (obj % groupid, item % varid, item % no_fill, item % fillvalue)
          case (nf90_float)
            status = nf90_def_var_fill (obj % groupid, item % varid, item % no_fill, real(item % fillvalue, kind=4))
          case (nf90_int64)
            status = nf90_def_var_fill (obj % groupid, item % varid, item % no_fill, int(item % fillvalue, kind=8))
          case (nf90_int)
            status = nf90_def_var_fill (obj % groupid, item % varid, item % no_fill, int(item % fillvalue, kind=4))
          case (nf90_short)
            status = nf90_def_var_fill (obj % groupid, item % varid, item % no_fill, int(item % fillvalue, kind=2))
          case (nf90_byte)
            status = nf90_def_var_fill (obj % groupid, item % varid, item % no_fill, int(item % fillvalue, kind=1))
          case (nf90_char)
            status = nf90_def_var_fill (obj % groupid, item % varid, item % no_fill, int(item % fillvalue, kind=1))
          case default
            call tell_error (tell_runtime_error, "tio_def_vars:  unsupported fill value type", errstat)
        end select
        if (status /= nf90_noerr) then
          errstat = tell_io_error
          call tell_set_error (errstat)
          return
        endif
      endif

      if (item % have_comment) then
        status = nf90_put_att (obj % groupid, item % varid, "comment", item % comment)
        if (status /= nf90_noerr) then
          errstat = tell_io_error
          call tell_set_error (errstat)
          return
        endif
      endif

      if (item % have_units) then
        status = nf90_put_att (obj % groupid, item % varid, "units", item % units)
        if (status /= nf90_noerr) then
          errstat = tell_io_error
          call tell_set_error (errstat)
          return
        endif
      endif

      if (item % have_valid_range) then
        status = nf90_put_att (obj % groupid, item % varid, "valid_min", item % valid_range(1))
        if (status /= nf90_noerr) then
          errstat = tell_io_error
          call tell_set_error (errstat)
          return
        endif
        status = nf90_put_att (obj % groupid, item % varid, "valid_max", item % valid_range(2))
        if (status /= nf90_noerr) then
          errstat = tell_io_error
          call tell_set_error (errstat)
          return
        endif
      endif

      if (associated(item % attlist)) then
        call tiof_def_atts (obj, item % attlist, item % varid, errstat)
        if (errstat < 0) return
      endif

      item => item % next
    enddo

  end subroutine tiof_def_vars

end module tio_module
