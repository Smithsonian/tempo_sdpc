program test_met
  use, intrinsic :: iso_c_binding
  use met_module
  implicit none
  type (c_ptr) :: met
  integer :: flags
  real (kind=4) :: lon, lat, psurf, ptrop
  integer, parameter :: num_isobars = 3
  real (kind=4), dimension(num_isobars) :: isobars, temp_on_isobar
  ! netcdf4 test is hard-coded
  character (len=*), parameter :: ncfile = 'data/test_met.nc'
  logical :: ncfile_exists
  integer, parameter :: nlev = 72
  real (kind=8), dimension(nlev) :: tprof
  integer :: i, status, errstat
  type (synth_met_type) :: smt

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

  if (i == 1) then
    write(*,*)'Usage: ftest_met GRIB_FILE [GRIB_FILE ...]'
    call met_list_free (met)
    call exit(0)
  endif

  isobars = (/900.0, 500.0, 100.0/)

  call met_list_interp_f (met, lon, lat, errstat, psurf)
  write(*,'(f8.1)')psurf

  call met_list_interp_f (met, lon, lat, errstat, ptrop=ptrop)
  write(*,'(f8.1)')ptrop

  call met_list_interp_f (met, lon, lat, errstat, psurf, ptrop)
  write(*,'(f8.1)')psurf
  write(*,'(f8.1)')ptrop

  call met_list_interp_f (met, lon, lat, errstat, &
                          psurf, ptrop, isobars, temp_on_isobar)
  write(*,'(f8.1)')psurf
  write(*,'(f8.1)')ptrop
  write(*,'(f8.1)')temp_on_isobar(:)

  call met_list_free (met)

  inquire (file=ncfile, exist=ncfile_exists)
  if (ncfile_exists) then
    call open_synth_met_data (smt, ncfile, errstat)
    call read_synth_met_data (smt, lat, lon, ptrop, errstat, &
                              psurf, tprof)
    call close_synth_met_data (smt, errstat)
    write(*,'(f8.1)')psurf
    write(*,'(f8.1)')ptrop
    write(*,'(f8.1)')tprof
  endif

end program
