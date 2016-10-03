!> Read in netCDF format Ring-effect reference data
!
!------------------------------------------------------------------------
!
!>   See Joiner et al. (1995), Appl. Optics, p.4513
!>   for more details of the approach used to create the table.
!
!   Variable descriptions, some guesswork involved:
!> @param errstat[inout] error reporting integer, non-zero indicates problem
!> @param   i0a: backscattered intensity at top of atmosphere
!> @param   i01a: ?
!> @param   z1: proportional to singly-scattered radiation?
!> @param   z2: proportional to twice-scattered radiation?
!> @param   k1bar: filling-in factor from singly-scattered photons
!> @param   nba, nia, nra: See Joiner et al eqns 26-30, avg number of 
!>      scatterings for reflected surface radiance scattered back
!>      to the groups (nba), radiation scattered back by the atmosphere 
!>      before it reaches the ground (nia) and radiation reflected by
!>      the ground and transmitted by the atmosphere (nra)?
!> @param   sba: component of reflected surface radiance that is reflected
!>        back to the ground by the atmosphere
!> @param   tra: surface radiance transmittance factor?
!
! !REVISION HISTORY:
!
!> @author  05Jan01   Joiner     original fortran 90
!> @author   12Aug02   Vasilkov   read filenames from PCF
!> @author  07Aug14   O'Sullivan added variable descriptions
!
!------------------------------------------------------------------------
module m_read_tables_tio

  private
  public read_tables_tio

contains

  subroutine read_tables_tio (errstat)

    use m_vars, ONLY: w_grid, nwave2, nwave, ntheta, nscan, nphi, &
         theta, scan, phi, sflx, wgrid_out, npres, pres,  k1bar, sba, nba, &
         i01a, i0a, tra, nia, nra, z1, z2
    use m_LUN_set
    use m_pgs_include
    use tio_module
    use tell_module
    use netcdf, only : nf90_nowrite


    implicit none

    ! input parameters
    integer, intent (inout) :: errstat

    !local variables
    integer :: ext_index
    integer :: status, version, pgs_pc_getreference
    character (len=128) :: ring_fn_nc, ring_fn, logmsg
    type (tiof_file_type) :: tio_ring_obj

    version = 1

    if (errstat /= 0) return

    !Get the name of the Ring effect table. If PCF file is out of date and
    ! still has .dat, change to .nc
    status = pgs_pc_getreference ( ring_id, version, ring_fn)
    if (status /= 0) then
      call tell_error (tell_io_read_error, &
           "read_tables_tio: failed to get Ring-effect table filename", &
           errstat)
      return
    endif
    ext_index=index(ring_fn,'.dat')
    if (ext_index > 0) then ! PCF still says .dat
      ring_fn_nc=ring_fn(1:ext_index-1)//'.nc'
    else ! PCF says .nc
      ring_fn_nc=ring_fn
    endif


    !Open the file and get dimension sizes
    call tiof_open (ring_fn_nc, tio_ring_obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (tio_ring_obj, "pres", npres, errstat)
    call tiof_inq_dimlen (tio_ring_obj, "phi", nphi, errstat)
    call tiof_inq_dimlen (tio_ring_obj, "scan", nscan, errstat)
    call tiof_inq_dimlen (tio_ring_obj, "theta", ntheta, errstat)
    call tiof_inq_dimlen (tio_ring_obj, "wave", nwave, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_tables_tio: failed to open Ring-effect table file", &
           errstat)
      return
    endif

    call tell_log(1,'read_tables_tio:')
    call tell_log(1,'nwave, ntheta, nscan, nphi, npres')
    write(logmsg,"(I5,I8,I7,I6,I7)") nwave, ntheta, nscan, nphi, npres
    call tell_log(1,logmsg)


    !Allocate arrays
    allocate(sflx(nwave), wgrid_out(nwave), theta(ntheta), scan(nscan), &
         phi(nphi), pres(npres), k1bar(nwave), sba(npres,nwave), &
         nba(npres,nwave), i01a(npres,nphi,nscan,ntheta,nwave), &
         i0a(npres,nscan,ntheta,nwave), z1(npres,nscan,ntheta,nwave), &
         z2(npres,nscan,ntheta,nwave), tra(npres,nscan,ntheta,nwave), &
         nia(npres,nphi,nscan,ntheta,nwave), nra(npres,nscan,ntheta,nwave), &
         w_grid(nwave), stat=errstat)

    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "read_tables_tio: allocation failure", &
           errstat)
      return
    endif


    !Read in data
    call tiof_get1d_r8 (tio_ring_obj, "theta", [0], [ntheta], &
         theta (1:ntheta), errstat)
    call tiof_get1d_r8 (tio_ring_obj, "scan", [0], [nscan], &
         scan (1:nscan), errstat)
    call tiof_get1d_r8 (tio_ring_obj, "phi", [0], [nphi], &
         phi (1:nphi), errstat)
    call tiof_get1d_r8 (tio_ring_obj, "pres", [0], [npres], &
         pres (1:npres), errstat)
    call tiof_get1d_r8 (tio_ring_obj, "wgrid_out", [0], [nwave], &
         wgrid_out (1:nwave), errstat)
    call tiof_get1d_r8 (tio_ring_obj, "sflx", [0], [nwave], &
         sflx (1:nwave), errstat)
    call tiof_get1d_r8 (tio_ring_obj, "k1bar", [0], [nwave], &
         k1bar (1:nwave), errstat)
    call tiof_get2d_r8 (tio_ring_obj, "nba", [0,0], [nwave, npres], &
         nba (1:npres, 1:nwave), errstat)
    call tiof_get2d_r8 (tio_ring_obj, "sba", [0,0], [nwave, npres], &
         sba (1:npres, 1:nwave), errstat)
    call tiof_get5d_r8 (tio_ring_obj, "i01a", [0,0,0,0,0], &
         [nwave, ntheta, nscan, nphi, npres], &
         i01a (1:npres, 1:nphi, 1:nscan, 1:ntheta, 1:nwave), errstat)
    call tiof_get5d_r8 (tio_ring_obj, "nia", [0,0,0,0,0], &
         [nwave, ntheta, nscan, nphi, npres], &
         nia (1:npres, 1:nphi, 1:nscan, 1:ntheta, 1:nwave), errstat)
    call tiof_get4d_r8 (tio_ring_obj, "nra", [0,0,0,0], &
         [nwave, ntheta, nscan, npres], &
         nra (1:npres, 1:nscan, 1:ntheta, 1:nwave), errstat)
    call tiof_get4d_r8 (tio_ring_obj, "tra", [0,0,0,0], &
         [nwave, ntheta, nscan, npres], &
         tra (1:npres, 1:nscan, 1:ntheta, 1:nwave), errstat)
    call tiof_get4d_r8 (tio_ring_obj, "i0a", [0,0,0,0], &
         [nwave, ntheta, nscan, npres], &
         i0a (1:npres, 1:nscan, 1:ntheta, 1:nwave), errstat)
    call tiof_get4d_r8 (tio_ring_obj, "z1", [0,0,0,0], &
         [nwave, ntheta, nscan, npres], &
         z1 (1:npres, 1:nscan, 1:ntheta, 1:nwave), errstat)
    call tiof_get4d_r8 (tio_ring_obj, "z2", [0,0,0,0], &
         [nwave, ntheta, nscan, npres], &
         z2 (1:npres, 1:nscan, 1:ntheta, 1:nwave), errstat)

    !Close file
    call tiof_close (tio_ring_obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, &
           "read_tables_tio: failed to read data", &
           errstat)
      return
    endif

    i0a=log(i0a)

    !Duplicate values used later
    nwave2=nwave
    w_grid=wgrid_out



  end subroutine read_tables_tio

end module m_read_tables_tio
