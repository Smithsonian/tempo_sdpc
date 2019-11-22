module climatology
  use, intrinsic :: iso_c_binding
  use tell_module
  use tio_module
  use clim_module
  implicit none
  private

  public c_clim_species_vtrop

  ! This must be interoperable with Pixel_Grid_Param_Type
  type, bind(C) :: pixel_grid_param_type
    real (kind=c_double) :: xmin, xmax, ymin, ymax
    integer (c_int) :: nx, ny
    integer (c_int) :: num_extra_xpoints, num_extra_ypoints
  end type

contains

  function c_clim_species_vtrop (grid, taix_beg, taix_end, vtrop) &
      bind (C, name="c_clim_species_vtrop")
    implicit none
    type (pixel_grid_param_type) :: grid
    real (kind=c_double), value :: taix_beg, taix_end
    real (kind=c_double), dimension(*), intent(out) :: vtrop
    integer (c_int) :: c_clim_species_vtrop

    type (clim_pres_type) :: cpt
    type (clim_species_type) :: cst
    type (clim_pres_bounds_type) :: bounds
    real (kind=4), dimension(:), allocatable :: pres_z, vmr_z, partial_column_z
    real (kind=c_double) :: dx, dy
    real (kind=c_float) :: lon, lat
    real (kind=4) :: ptrop
    integer :: nx, ny, nz, nlayers, errstat, i, j, k, idx
    integer, dimension(2) :: year, month, day
    real (kind=8), dimension(2) :: hour

    errstat = 0

    call tiof_taix_time_to_utc_caldate (taix_beg, year(1), month(1), day(1), hour(1), errstat)
    call tiof_taix_time_to_utc_caldate (taix_end, year(2), month(2), day(2), hour(2), errstat)
    c_clim_species_vtrop = errstat
    if (errstat /= 0) return

    bounds % hour_beg = real (hour(1), kind=4)
    bounds % hour_end = real (hour(2), kind=4)
    bounds % lon_min  = real (grid % xmin, kind=4)
    bounds % lon_max  = real (grid % xmax, kind=4)
    bounds % lat_min  = real (grid % ymin, kind=4)
    bounds % lat_max  = real (grid % ymax, kind=4)

    call clim_pres_init (cpt, month(1), day(1), bounds, errstat)
    call clim_species_init (cst, cpt, 'NO2', errstat)
    c_clim_species_vtrop = errstat
    if (errstat /= 0) return

    dx = (grid % xmax - grid % xmin) / grid % nx
    dy = (grid % ymax - grid % ymin) / grid % ny

    nx = grid % nx
    ny = grid % ny
    nz = clim_pres_nz (cpt)
    nlayers = nz-1

    allocate (pres_z(0:nz-1), vmr_z(nlayers), partial_column_z(nlayers))

    do i = 1, nx
      lon = real (grid % xmin + (i-1) * dx, kind=4)
      do j = 1, ny
        lat = real (grid % ymin + (j-1) * dy, kind=4)

        call clim_pres (cpt, bounds % hour_beg, lon, lat, pres_z, errstat, p_trop=ptrop)
        call clim_species_vmr (cst, cpt, bounds % hour_beg, lon, lat, vmr_z, errstat)
        call clim_partial_column (pres_z, vmr_z, partial_column_z, errstat)
        c_clim_species_vtrop = errstat
        if (errstat /= 0) return

        do k = nlayers, 0, -1
          if (pres_z(k) > ptrop) exit
        enddo
        if (k == 0) then
          call tell_error (tell_runtime_error, &
                           "c_clim_species_vtrop: invalid troposphere pressure", errstat)
          c_clim_species_vtrop = errstat
          return
        endif

        idx = i + (j-1) * nx
        vtrop(idx) = sum(partial_column_z(1:k))
      enddo
    enddo

    c_clim_species_vtrop = 0

  end function

end module
