!> Routines to allocate and deallocate arrays
module m_allocate

  use OMSAO_variables_module, only: yn_varyslit, slit_rad, scnwrt
  use OMSAO_omidata_module, only: omi_solslitfit, omi_radslitfit, &
       omi_radiance_spec, omi_radiance_prec, omi_radiance_wavl, &
       omi_radiance_qflg, omi_fitvar, radwind
  use tell_module

  implicit none

  public alloc, dealloc
  private

contains

  subroutine alloc(errstat)
    use OMSAO_parameters_module, only: max_fit_pts
    use OMSAO_omidata_module, only: nxtrack_max, nlines_max
    use OMSAO_indices_module, only: max_calfit_idx, n_max_fitpars

    implicit none

    !output variables
    integer (kind=4), intent(inout) :: errstat

    if (errstat /= 0) return

    ! Variables used in every run
    if (scnwrt) print *, 'Allocating variables'
    allocate (omi_radiance_spec(max_fit_pts,nxtrack_max,0:nlines_max-1), &
         omi_radiance_prec(max_fit_pts,nxtrack_max,0:nlines_max-1), &
         omi_radiance_wavl(max_fit_pts,nxtrack_max,0:nlines_max-1), &
         stat = errstat)
    allocate (radwind(max_fit_pts, nxtrack_max,0:nlines_max-1), &
         omi_radiance_qflg(max_fit_pts, nxtrack_max,0:nlines_max-1), &
         stat = errstat)
    allocate (omi_fitvar(nxtrack_max,0:nlines_max-1, n_max_fitpars), &
         stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "alloc: failed to allocate variables", errstat)
      return
    endif

    ! Variables used in variable slit fitting, see omi_cross_calibrate
    if (yn_varyslit) then
      if (scnwrt) print *, 'Allocating slit fitting variables'
      allocate (omi_solslitfit(max_fit_pts, max_calfit_idx, 2, nxtrack_max), &
           stat=errstat)
      omi_solslitfit = 0.0d0
      if (slit_rad) then
        allocate (omi_radslitfit(max_fit_pts, max_calfit_idx, 2, nxtrack_max),&
             stat=errstat)
        omi_radslitfit = 0.0d0
      endif
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "alloc: failed to allocate slit fitting variables", errstat)
        return
      endif
    endif

  end subroutine alloc


  subroutine dealloc (errstat)

    implicit none

    !output variables
    integer (kind=4), intent(inout) :: errstat

    ! used in every run
    if (errstat /= 0) return

    if (scnwrt) print *, "Deallocating variables"
    if (allocated(omi_radiance_spec)) deallocate (omi_radiance_spec, &
         stat=errstat)
    if (allocated(omi_radiance_prec)) deallocate (omi_radiance_prec, &
         stat=errstat)
    if (allocated(omi_radiance_wavl)) deallocate (omi_radiance_wavl, &
         stat=errstat)
    if (allocated(omi_radiance_qflg)) deallocate (omi_radiance_qflg, &
         stat=errstat)
    if (allocated(radwind)) deallocate (radwind, stat=errstat)
    if (allocated(omi_fitvar)) deallocate (omi_fitvar, stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "alloc: failed to deallocate variables", errstat)
      return
    endif


    ! Variables used in variable slit fitting, see omi_cross_calibrate
    if (yn_varyslit) then
      if (scnwrt) print *, 'Deallocating slit fitting variables'
      if (allocated(omi_solslitfit)) deallocate(omi_solslitfit, stat=errstat)
      if (allocated(omi_radslitfit)) deallocate(omi_radslitfit, stat=errstat)
      if (errstat /= 0) then
        call tell_error (tell_malloc_error, &
             "alloc: failed to deallocate slit fitting variables", errstat)
        return
      endif
    endif

  end subroutine dealloc


end module m_allocate
