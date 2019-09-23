program test_met
  use, intrinsic :: iso_c_binding
  use met_module
  implicit none
  type (c_ptr) :: met
  integer :: flags
  real (kind=4) :: lon, lat, psurf, ptrop
  integer, parameter :: num_isobars = 3
  real (kind=4), dimension(num_isobars) :: isobars, temp_on_isobar
  integer :: i, status, errstat

  character (kind=c_char, len=1024) :: argbuf

  errstat = 0

  lon = -90.0
  lat = 36.0

  flags = 7

  met = met_list_new (flags)

  i = 1
  do
    call get_command_argument (i, argbuf)
    if (len_trim(argbuf) == 0) exit
    status = met_list_add_file (met, trim(argbuf)//c_null_char)
    if (status /= 0) stop 1
    i = i + 1
  enddo

  isobars = (/9.e4, 5.e4, 1.e4/)

  call met_list_interp_f (met, lon, lat, errstat, psurf)
  write(*,*)'psurf=',psurf

  call met_list_interp_f (met, lon, lat, errstat, ptrop=ptrop)
  write(*,*)'ptrop=',ptrop

  call met_list_interp_f (met, lon, lat, errstat, psurf, ptrop)
  write(*,*)'psurf=',psurf,' ptrop=',ptrop

  call met_list_interp_f (met, lon, lat, errstat, &
                          psurf, ptrop, isobars, temp_on_isobar)
  write(*,*)'psurf=',psurf,' ptrop=',ptrop
  write(*,*)'isobars=', isobars(:)
  write(*,*)'temp_on_isobar=', temp_on_isobar(:)

  call met_list_free (met)

end program
