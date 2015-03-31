module m_alloc2

contains

  subroutine alloc2(errstat)

    use m_cloud_pres_mod
    use tell_module
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  alloc2
    ! 
    ! !DESCRIPTION: alloc2 allocates/deallocates memory for retrievals	    
    !
    ! !CALLING SEQUENCE: 
    !
    !        call alloc2
    !     
    ! !INPUT PARAMETERS:   
    !
    ! !OUTPUT PARAMETERS:  
    !
    ! !SEE ALSO:  
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner      original fortran 90
    !  26Mar15   O'Sullivan  update for TEMPO
    !
    !EOP
    !-------------------------------------------------------------------------

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
