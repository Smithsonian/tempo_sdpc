!>Memory allocation for matrices (covariance, jacobians, etc)
module m_alloc2

contains

  !-------------------------------------------------------------------------
  !
  ! !ROUTINE:  alloc2
  ! 
  ! !DESCRIPTION: 
  !> alloc2 allocates/deallocates memory for matrices (covariance, 
  !! jacobians, etc)
  !
  ! !CALLING SEQUENCE: 
  !
  !        call alloc2
  !     
  ! !INPUT PARAMETERS:   
  !> @param errstat error reporting integer, non-zero = failure
  !
  ! !REVISION HISTORY: 
  !
  !> @author   05Jan01   Joiner      original fortran 90
  !> @author   26Mar15   O'Sullivan  updated for TEMPO
  !
  !EOP
  !-------------------------------------------------------------------------
  subroutine alloc2(errstat)

    use m_cloud_pres_mod
    use tell_module
    implicit none

    integer, intent(inout) :: errstat

    if (errstat /= 0) return


    !deallocate memory
    !=================
    if (allocated(x)) then
      deallocate(x, x_fg, h, htr, err_cov, corr, b_i, stat=errstat)
    endif
    if(allocated(y_back)) deallocate(y_back, stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "alloc2: deallocation failure", &
           errstat)
      return
    endif


    !Allocate memory
    allocate(x(0:nst-1,1), &
         x_fg(0:nst-1,1), &
         h(0:nobs-1,0:nst-1), &
         htr(0:nst-1,0:nobs-1), &   
         err_cov(0:nst-1,0:nst-1), &   
         corr(0:nst-1,0:nst-1), &   
         b_i(0:nst-1), & 
         y_back(0:nst-1,1), &
         stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "alloc2: allocation failure", &
           errstat)
      return
    endif


  end subroutine alloc2

end module m_alloc2
