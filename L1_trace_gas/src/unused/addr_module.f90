module addr_module
  use iso_c_binding, only: c_ptr, c_loc
  ! It would be nice to put these in an interface block to
  ! define a generic function.  However, this does not
  ! appear to be possible.
  public addr_i2, addr_i4, addr_r4, addr_r8
contains
  function addr_i2 (a) result (addr)
    use OMSAO_precision_module, only: i2
    integer(kind=i2), target :: a(*)
    type (c_ptr) :: addr
    addr = c_loc (a)
  end function addr_i2

  function addr_i4 (a) result (addr)
    use OMSAO_precision_module, only: i4
    integer(kind=i4), target :: a(*)
    type (c_ptr) :: addr
    addr = c_loc (a)
  end function addr_i4

  function addr_r4 (a) result (addr)
    use OMSAO_precision_module, only: r4
    real(kind=r4), target :: a(*)
    type (c_ptr) :: addr
    addr = c_loc (a)
  end function addr_r4

  function addr_r8 (a) result (addr)
    use OMSAO_precision_module, only: r8
    real(kind=r8), target :: a(*)
    type (c_ptr) :: addr
    addr = c_loc (a)
  end function addr_r8

end module
