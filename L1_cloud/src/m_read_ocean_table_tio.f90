! Read in netCDF Ocean Ring effect reference data
module m_read_ocean_table_tio

  private
  public read_ocean_table_tio

contains

  subroutine read_ocean_table_tio (errstat)

    use m_vars, only: nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl, &
         oc_perms, nwave2, wgrid_out_oc, wgrid_oc, oc_table, w_grid, &
         theta_oc,scan_oc,phi_oc,ocrefl,chl,iprt
    use m_interpol
    use m_LUN_set
    use m_pgs_include
    use tio_module
    use tell_module
    use netcdf, only : nf90_nowrite

    implicit none

    !Input variables
    integer, intent(inout) :: errstat

    !Local variables
    real (KIND=8), allocatable, dimension(:) :: oc_perms2
    ! ocean Raman correction coefficient
    real (KIND=8), parameter :: coef = 1.5
    integer :: i, j, m, ext_index
    integer :: version, status, pgs_pc_getreference
    character (len=128) :: oc_fn_nc, oc_fn
    type (tiof_file_type) :: tio_oc_obj

    !temporary variables for checking data matches binary table
    real (kind=8), dimension(:), allocatable :: t_theta_oc, t_scan_oc
    real (kind=8), dimension(:), allocatable :: t_ocrefl, t_chl, t_phi_oc
    real (kind=4), dimension(:,:,:,:,:,:), allocatable :: t_oc_perms 
    real (kind=4), dimension(:), allocatable :: t_wgrid_out_oc 
    real (kind=4) :: tol
    tol=1e-8

    version = 1

    if (errstat /= 0) return

    !Get the name of the Ocean file
    status = pgs_pc_getreference ( oc_ram_id, version, oc_fn)
    if (status /= 0) then
      errstat = -1
      call tell_error (tell_io_read_error, &
           "read_ocean_table_tio: failed to get Ocean Raman filename", &
           errstat)
      return
    endif
    ext_index=index(oc_fn,'.dat')
    oc_fn_nc=oc_fn(1:ext_index-1)//'.nc'


    !Open the file and get dimension sizes
    call tiof_open (oc_fn_nc, tio_oc_obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (tio_oc_obj, "chl", nchl, errstat)
    call tiof_inq_dimlen (tio_oc_obj, "ocrefl", nocrefl, errstat)
    call tiof_inq_dimlen (tio_oc_obj, "phi_oc", nphi_oc, errstat)
    call tiof_inq_dimlen (tio_oc_obj, "scan_oc", nscan_oc, errstat)
    call tiof_inq_dimlen (tio_oc_obj, "theta_oc", nthet_oc, errstat)
    call tiof_inq_dimlen (tio_oc_obj, "wgrid_out_oc", nwave_oc, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_ocean_table_tio: failed to open Ocean Raman file", &
           errstat)
      return
    endif

    if (iprt >= 1) then
      print *, 'read_ocean_table_tio: nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl'
      print *, nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl
    endif

    ! Allocate memory - first arrays that should be unallocated
    allocate(oc_perms(nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl), &
         oc_perms2(nwave_oc), wgrid_out_oc(nwave_oc), stat=errstat)
    !Next arrays that may have been allocated if binary table read as well
    if (.not. allocated(chl)) then
      allocate (oc_table(nthet_oc,nscan_oc,nchl,nwave2), &
           theta_oc(nthet_oc), scan_oc(nscan_oc), phi_oc(nphi_oc), &
           wgrid_oc(nwave_oc), ocrefl(nocrefl), chl(nchl), stat=errstat)
    endif
    !Temporary arrays for value checking
    allocate (t_oc_perms(nwave_oc,nthet_oc,nscan_oc,nphi_oc,nocrefl,nchl), &
         t_wgrid_out_oc(nwave_oc), t_theta_oc(nthet_oc), t_scan_oc(nscan_oc), &
         t_phi_oc(nphi_oc), t_ocrefl(nocrefl), t_chl(nchl), &
         stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_ocean_table_tio: allocation failure", &
           errstat)
      return
    endif

    !Read in data
!    call tiof_get6d_r4 (tio_oc_obj, "oc_perms", [0,0,0,0,0,0], &
!       [nchl, nocrefl, nphi_oc, nscan_oc, nthet_oc, nwave_oc], &
!       t_oc_perms(1:nwave_oc,1:nthet_oc,1:nscan_oc,1:nphi_oc,1:nocrefl,1:nchl), &
!       errstat)
    call tiof_get1d_r8 (tio_oc_obj, "chl", [0], [nchl], &
         t_chl (1:nchl), errstat)
    call tiof_get1d_r8 (tio_oc_obj, "ocrefl", [0], [nocrefl], &
         t_ocrefl (1:nocrefl), errstat)
    call tiof_get1d_r8 (tio_oc_obj, "phi_oc", [0], [nphi_oc], &
         t_phi_oc (1:nphi_oc), errstat)
    call tiof_get1d_r8 (tio_oc_obj, "scan_oc", [0], [nscan_oc], &
         t_scan_oc (1:nscan_oc), errstat)
    call tiof_get1d_r8 (tio_oc_obj, "theta_oc", [0], [nthet_oc], &
         t_theta_oc (1:nthet_oc), errstat)
    call tiof_get1d_r4 (tio_oc_obj, "wgrid_out_oc", [0], [nwave_oc], &
         t_wgrid_out_oc (1:nwave_oc), errstat)

    !Close file
    call tiof_close (tio_oc_obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_ocean_table_tio: failed to read data", &
           errstat)
      return
    endif

    !Temporary: check values
    if (any(abs(t_chl - chl) > tol)) print *, "mismatch: chl"
    if (any(abs(t_ocrefl - ocrefl) > tol)) print *, "mismatch: ocrefl"
    if (any(abs(t_phi_oc - phi_oc) > tol)) print *, "mismatch: phi_oc"
    if (any(abs(t_scan_oc - scan_oc) > tol)) print *, "mismatch: scan_oc"
    if (any(abs(t_theta_oc - theta_oc) > tol)) print *, "mismatch: theta_oc"
!    if (any(abs(t_oc_perms - oc_perms) > tol)) print *, "mismatch: oc_perms"

    !Temporary: copy values
    chl = t_chl
    ocrefl = t_ocrefl
    phi_oc = t_phi_oc
    scan_oc = t_scan_oc
    theta_oc = t_theta_oc
    wgrid_out_oc = t_wgrid_out_oc
!    oc_perms = t_oc_perms

    ! Reformatting
    wgrid_oc=wgrid_out_oc
    do i=1,nthet_oc
      do j=1,nscan_oc
        do m=1,nchl
          ! use 10% reflectivity
          oc_perms2(:)=oc_perms(:,i,j,1,2,m)
          oc_table(i,j,m,:)= coef * interpol(oc_perms2(:),wgrid_oc,w_grid)
        enddo
      enddo
    enddo

    ! Debugging output
    if (iprt >= 6) then
      print *,'ocean_table wavelengths'
      write(6,'(6f12.3)') wgrid_out_oc
      print *,'ocean_table_new (%)'
      write(6,'(6f12.3)') oc_table(1,1,1,:)*100
    endif

    ! Deallocate unnecessary arrays
    deallocate(oc_perms2)
    deallocate(oc_perms)
    deallocate(wgrid_out_oc)

    end subroutine read_ocean_table_tio

end module m_read_ocean_table_tio
 
