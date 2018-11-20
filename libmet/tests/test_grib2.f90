program test_grib2

  use met_module

  implicit none

  ! inputs
  real (kind=4), parameter :: lat = 42.0d0, lon = -98.0d0
  character(len=11), parameter :: metfile = 'small.grib2'
  real (kind=4), parameter :: tolerance = 1e-8

  ! correct outputs
  real (kind=4), parameter :: corr_troppres = 232.674576, &
       corr_surfpres = 950.577454
  real (kind=4), dimension(10), parameter :: corr_tprof = &
       (/282.000000, 275.600006, 277.825745, 276.299438, 267.428619, &
       257.895905, 244.879578, 227.800003, 219.300003, 214.706390 /)

  integer (kind=4) :: errstat
  real (kind=4) :: troppres, surfpres
  real (kind=4), dimension(10) :: tprof

  errstat = 0

  call read_met_data (metfile, lat, lon, troppres, surfpres, tprof, errstat)
  if (errstat /= 0) then
    print *, "*** test_grib2: failed to read GRIB2 file"
    stop 1
  endif

  if (abs(troppres-corr_troppres) .gt. tolerance) then
    print *, "*** test_grib2: failed: troposphere pressure incorrect"
    stop 1
  endif

  if (abs(surfpres-corr_surfpres) .gt. tolerance) then
    print *, "*** test_grib2: failed: surface pressure incorrect"
    stop 1
  endif

  if (any(abs(tprof-corr_tprof) .gt. tolerance)) then
    print *, "*** test_grib2: failed:  incorrect"
    stop 1
  endif


end program test_grib2
