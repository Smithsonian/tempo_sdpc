module errormodule

  private
  public err_message_error

contains
  subroutine err_message_error (msg, errcode)

    implicit none
    character (len=*), intent(in) :: msg
    integer, intent(inout) :: errcode

    write (*,*) "ERROR: ", trim(msg)
    if (errcode >= 0) errcode = -1
    return
  end subroutine

end module
