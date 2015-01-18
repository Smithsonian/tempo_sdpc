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

  type, public :: tiof_object_type
    integer :: fileid = -1
    integer :: groupid = -1
  end type tiof_object_type

  type, public :: tiof_dim_type
    character (len=tiof_max_name_len) :: name
    integer :: len_name = 0
    integer :: len = 0
    integer :: dimid = -huge(1)
    type (tiof_dim_type), pointer :: next => null()
  end type

  type, public :: tiof_dimlist_type
    integer :: num_items = 0
    type (tiof_dim_type), pointer :: head => null(), tail => null()
  end type

  type, public :: tiof_att_type
    character (len=tiof_max_name_len) :: name
    integer :: len_name = 0
    integer :: xtype = -1
    character (len=tiof_max_att_len) :: att_text
    integer (kind=i4), allocatable, dimension(:) :: att_i4
    real (kind=4), allocatable, dimension(:) :: att_r4
    real (kind=8), allocatable, dimension(:) :: att_r8
    type (tiof_att_type), pointer :: next => null()
  end type

  type, public :: tiof_attlist_type
    integer :: num_items = 0
    type (tiof_att_type), pointer :: head => null(), tail => null()
  end type

  type, public :: tiof_var_type
    character (len=tiof_max_name_len) :: name
    integer :: len_name = 0
    integer :: xtype = -1
    integer :: varid = -huge(1)
    integer :: rank = 0
    logical :: have_comment=.false., have_units=.false., have_valid_range=.false., &
      have_fillvalue = .false.
    character (len=tiof_max_att_len) :: comment, units
    real (kind=8), dimension(2) :: valid_range = [0.0, 0.0]
    integer, dimension(tiof_max_var_dims) :: dimids
    integer :: deflate_level = 0
    logical :: contiguous = .true., shuffle = .false.
    integer, dimension(tiof_max_var_dims) :: chunksizes = 0
    integer :: no_fill = 0
    real (kind=8) :: fillvalue
    type (tiof_attlist_type), pointer :: attlist => null()
    type (tiof_var_type), pointer :: next => null()
  end type

  type, public :: tiof_varlist_type
    integer :: num_items = 0
    type (tiof_var_type), pointer :: head => null(), tail => null()
  end type

  integer :: tiof_get_var_section, tiof_put_var_section
  external   tiof_get_var_section, tiof_put_var_section

  public tiof_create, tiof_open, tiof_close, &
    tiof_put_var_section, tiof_get_var_section, &
    tiof_inq_group, tiof_inq_dimlen, &
    tiof_dimlist_append, tiof_dimlist_lookup, tiof_def_dims, &
    tiof_varlist_append, tiof_varlist_lookup, tiof_def_vars, &
    tiof_attlist_append,                      tiof_def_atts

  include 'get_put_decl.inc'

contains

  include 'get_put_code.inc'

  subroutine tiof_create (obj, file, create_mode, errstat)
    implicit none
    type (tiof_object_type), intent(out) :: obj
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

  end subroutine tiof_create

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
      call tell_error (tell_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif

    status = nf90_inquire_dimension (obj % groupid, dimid, len=dimlen)
    if (status /= nf90_noerr) then
      call tell_error (tell_io_read_error, "accessing dimension "//trim(name)//" ("//trim(nf90_strerror(status))//")", errstat)
      return
    endif
  end subroutine tiof_inq_dimlen

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

  subroutine tiof_dimlist_lookup (list, names, dimids, errstat, &
                                  dimsizes)
    implicit none
    type (tiof_dimlist_type), intent(in) :: list
    character (len=*), dimension(:), intent(in) :: names
    integer, dimension(:), intent(out) :: dimids
    integer, intent(inout) :: errstat
    integer, optional, intent(inout), dimension(size(dimids)) :: dimsizes

    type (tiof_dim_type), pointer :: item => null()
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

  subroutine tiof_def_dims (obj, list, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    type (tiof_dimlist_type), intent(in) :: list
    integer, intent(inout) :: errstat

    type (tiof_dim_type), pointer :: item => null()
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

  subroutine tiof_attlist_append (list, errstat, att_name, &
                                  att_i4, att_r4, att_r8, att_text)
    implicit none
    type (tiof_attlist_type), intent(inout) :: list
    character (len=*), intent(in) :: att_name
    integer, intent(inout) :: errstat
    integer (kind=i4), optional, dimension(:) :: att_i4
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

  subroutine tiof_def_atts (obj, list, varid, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    type (tiof_attlist_type), intent(in) :: list
    integer, intent(in) :: varid
    integer, intent(inout) :: errstat

    type (tiof_att_type), pointer :: item => null()
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

  subroutine tiof_varlist_lookup (list, name, var_ptr, errstat)
    implicit none
    type (tiof_varlist_type), intent(in) :: list
    character (len=*), intent(in) :: name
    type (tiof_var_type), pointer, intent(out) :: var_ptr => null()
    integer, intent(inout) :: errstat

    type (tiof_var_type), pointer :: item => null()
    integer :: len_name

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

  subroutine tiof_def_vars (obj, list, errstat)
    implicit none
    type (tiof_object_type), intent(in) :: obj
    type (tiof_varlist_type), intent(in) :: list
    integer, intent(inout) :: errstat

    type (tiof_var_type), pointer :: item => null()
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
