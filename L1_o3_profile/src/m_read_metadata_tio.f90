!> Subroutines to read in metadata from L1 radiance netCDF file
module m_read_metadata_tio

  implicit none

  public read_date_tio
  private

contains

  !> Get year, month, day and julian day of observation start
  !---------------------------------------------------------------------
  !
  !> @param[in] l1file L1 netCDF radiance file name
  !> @param[out] year year of observation start
  !> @param[out] month month of observation start
  !> @param[out] day day of observation start
  !> @param[out] jday julian day of observation start
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan    July 2016
  !---------------------------------------------------------------------
  subroutine read_date_tio ( l1file, year, month, day, jday, errstat )

    use netcdf
    use tio_module
    use tell_module
    use OMSAO_parameters_module, only : maxchlen
    use m_utilities, only: day_of_year

    implicit none

    !input variables
    character (len=*), intent (in) :: l1file

    !ouput variables
    integer, intent (out) :: year, month, day, jday
    integer, intent (inout) :: errstat

    !local variables
    integer :: ncerr
    character (len=maxchlen) :: rbd_string
    type(tiof_file_type) :: tio_l1obj


    if (errstat /= 0) return

    call tiof_open (l1file, tio_l1obj, nf90_nowrite, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_date_tio: failed to open L1 radiance file", errstat)
      return
    endif

    ncerr = nf90_get_att (tio_l1obj % fileid, nf90_global, &
         "time_coverage_start", rbd_string)
    if (ncerr /= nf90_noerr) then
      call tell_error (tell_io_read_error, &
        "read_date_tio: failed to read global attribute time_coverage_start", &
           errstat)
      return
    endif

    call tiof_close (tio_l1obj, errstat)
    if (errstat /= 0) then
      call tell_error (tell_io_error, &
           "read_date_tio: failed to close L1 radiance file", errstat)
      return
    endif

    read (rbd_string, '(i4,1x,i2,1x,i2)') year, month, day
    jday = day_of_year ( year, month, day ) 

  end subroutine read_date_tio


end module m_read_metadata_tio
