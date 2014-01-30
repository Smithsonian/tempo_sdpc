module errormodule

  private
  public err_message_error, err_message_warn

contains
  subroutine err_message_error (msg, errcode)

    implicit none
    character (len=*), intent(in) :: msg
    integer, intent(inout) :: errcode

    write (*,*) "ERROR: ", trim(msg)
    if (errcode >= 0) errcode = -1
    return
  end subroutine

  subroutine err_message_warn (msg)

    implicit none
    character (len=*), intent(in) :: msg

    write (*,*) "WARNING: ", trim(msg)
    return
  end subroutine

end module
