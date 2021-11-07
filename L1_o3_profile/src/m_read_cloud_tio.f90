!> Read cloud fraction, pressure, and quality flags from L2 cloud netCDF fil
module m_read_cloud_tio
  use o3p_names_module
  use tio_module
  use tell_module
  use netcdf, only: nf90_nowrite
  use OMSAO_tmpodata_module, ONLY: nxtrack_max, ntimes_max
  use OMSAO_variables_module, only:  scnwrt, nxbin, nybin, l2file=>l2_cld_filename

  implicit none
  ! Type declaration for cloud data
  type, public :: tempo_cloud_block
    real (kind=4) :: CFRmissing, CFRscale, CFRoffset, CTPmissing, &
         CTPscale, CTPoffset
    real (kind=4),   dimension (nxtrack_max,0:ntimes_max-1)  :: cfr
    real (kind=4),   dimension (nxtrack_max,0:ntimes_max-1)  :: ctp
    real (kind=4),   dimension (nxtrack_max,0:ntimes_max-1)  :: ai
    integer(kind=1), dimension (nxtrack_max,0:ntimes_max-1) :: qflags
  end type tempo_cloud_block
  type (tempo_cloud_block), public :: L2_cloud


  public read_cloud_tio
  private read_cloud_dims, read_cloud_data, fill_in_tempoctp


contains

  !> Top-level subroutine to read in an L2 netCDF cloud file
  !----------------------------------------------------------------------
  !
  !> @param[in] l2file filename for L2 netCDF cloud file
  !> @param[in] nstep number of mirror steps in use
  !> @param[in] nxtrack number of xtrack positions in use
  !> @param[in] nl number of binned lines
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan June 2016
  !----------------------------------------------------------------------

  subroutine read_cloud_tio (ntimes, nxtrack,sline,eline, errstat)

    use m_convert_coadd, only: convert_2bytes_to_16bits
    ! To allow comparison with values read in from he5

    implicit none

    !input variables
    integer (kind=4), intent (in) :: sline, eline, ntimes, nxtrack

    !output variables
    integer (kind=4), intent (out) :: errstat

    !local variables
    INTEGER :: ntimes_loc, nxtrack_loc
    integer (kind=4), parameter :: nbit = 16
    integer (kind=4) :: ix, iix, i, ii, j, k, &
                       nbx, nbin, nl
    real    (kind=4) :: scfr, scfr1, scfr0, tmpsum
    integer (kind=2), dimension(ntimes) :: tmp_byte_num
    integer (kind=2), dimension(ntimes, 0:nbit-1) :: tmp_bit_num
    type (tiof_file_type) :: tio_l2obj
    ! Save variables
    real    (kind=4), dimension (:,:),  POINTER :: cfr, ctp
    integer (kind=2), dimension (:,:),  POINTER :: qflag
    integer (kind=2), dimension (:,:,:),POINTER :: flgbits 
    !(nxtrack_max,0:ntimes_max-1, 0:nbit-1) :: flgbits

    allocate (flgbits(nxtrack_max,0:ntimes_max-1, 0:nbit-1))
    allocate ( cfr (nxtrack_max, 0:ntimes_max-1))
    allocate ( ctp (nxtrack_max, 0:ntimes_max-1))
    allocate ( qflag (nxtrack_max, 0:ntimes_max-1))

      !get dimensions of L2 cloud file
      call read_cloud_dims(l2file, tio_l2obj, ntimes_loc, nxtrack_loc,errstat)
      if (errstat /= 0) then 
      call tell_error (tell_io_read_error, &
           "read_cloud_dims: failed",  errstat)
      return
      endif

      !Check dimensions are consistent with input radiance data
      if (ntimes_loc /= ntimes .OR. nxtrack_loc /= nxtrack) then
        call tell_error (tell_io_error, &
           "inconsistent dimensions between radiance and cloud files", errstat)
        return
      endif
      
      !Read cloud fraction, pressure, quality flags for the swath
      call read_cloud_data (l2file, tio_l2obj, ntimes, nxtrack, &
           cfr, ctp, qflag, errstat)
      if (errstat /= 0) then 
        call tell_error (tell_io_read_error, "read_cloud_data: failed",  errstat)
      return
      endif

      flgbits = 0
      DO ix = 1, nxtrack
      tmp_byte_num=qflag(ix, 0:ntimes-1)
      CALL convert_2bytes_to_16bits (nbit,ntimes, tmp_byte_num, &
           tmp_bit_num)
      flgbits(ix, 0:ntimes-1, 0:nbit-1)=tmp_bit_num
      ENDDO

     ! Fill in cloud top pressure values for bad pixels (interpolation/extrapolation)
     ! 0   - failed convergence check
     ! 1   - solar zenith angle, lat., lon., out of range (SZA > 88 deg) 
     ! 2   - cloud pressure less than low range of table
     ! 3   - cloud pressure greater than surface pressure
     ! 4   - matrix inversion failed
     ! 5   - snow/ice (if second byte of GroundPixelQualityFlags is 50-130)
     ! 6   - reflectivity < 0 or > 1.0
     ! 7   - bad radiances detected
     ! 8   - aerosol index flag
     ! 9   - radiance PixelQuality error
     ! 10  - radiance PixelQuality warning
     ! 11  - irradiance PixelQuality error
     ! 12  - irradiance PixelQuality warning
     ! 13  - effective surface pressure retrieved because cloud fraction < 0.05
     ! 14  - missing data
     ! 15  - geolocation error
     qflag(1:nxtrack, :)  = &
         flgbits(1:nxtrack, :, 0) + flgbits(1:nxtrack, :, 1) + &
         flgbits(1:nxtrack, :, 2) + flgbits(1:nxtrack, :, 3) + &
         flgbits(1:nxtrack, :, 2) + flgbits(1:nxtrack, :, 3) + &
         flgbits(1:nxtrack, :, 4) + flgbits(1:nxtrack, :, 6) + &
         flgbits(1:nxtrack, :, 7) + flgbits(1:nxtrack, :, 9) + &
         flgbits(1:nxtrack, :, 11) + flgbits(1:nxtrack, :, 13) + &
         flgbits(1:nxtrack, :, 14) + flgbits(1:nxtrack, :, 15) 

    ! Fill in cloud top pressure values for bad pixels
     call fill_in_tempoctp(nxtrack, ntimes, ctp(1:nxtrack, 0:ntimes-1), &
         qflag(1:nxtrack, 0:ntimes-1))
     

    ! Move cloud arrays into cloud block and rebin
    nl = (eline-sline+1)/nybin
    nbin = nxbin
    nbx  = nxtrack / nbin
  
    L2_cloud%cfr   (1:nbx, 0:nl-1) = 0.0
    L2_cloud%ctp   (1:nbx, 0:nl-1) = 0.0
    L2_cloud%qflags(1:nbx, 0:nl-1) = 0

    do ix = 1, nbx
      do i = 0, nl-1
        iix = (ix - 1) * nbin + 1 
        ii  = i * nybin + sline-1

        scfr = 0.0
        scfr1 = 0.0
        scfr0 = 0.0
        tmpsum = 0.0
        do j = iix, iix + nbin - 1
          do k = ii, ii + nybin - 1
            if (cfr(j, k) >= 0.0) then
              L2_cloud%cfr(ix, i) = L2_cloud%cfr(ix, i) + cfr(j, k)
              scfr1 = scfr1 + 1.0
            endif

            if ( ctp(j, k) > 0.0 .and. cfr(j, k) >= 0.0 ) then                 
              L2_cloud%ctp(ix, i) = &
                   L2_cloud%ctp(ix, i) + log(ctp(j, k)) * cfr(j, k)
              tmpsum = tmpsum + log(ctp(j, k)) 
              scfr = scfr + cfr(j, k)
              scfr0 = scfr0 + 1.0
            endif
          enddo
        enddo

        if (scfr /= 0.0) then       ! Weighted by Cloud Fraction
          L2_cloud%ctp(ix, i) = exp(L2_cloud%ctp(ix, i) / scfr)
        else if (scfr0 > 0.0 ) then ! Simple average if cloud fraction is all zero
          L2_cloud%ctp(ix, i) = exp(tmpsum / scfr0)
        else
          L2_cloud%ctp(ix, i) = 0.0
        endif

        if (scfr1 /= 0.0) then
          L2_cloud%cfr(ix, i) = L2_cloud%cfr(ix, i) / scfr1
        else
          L2_cloud%cfr(ix, i) = 0.0
        endif
      enddo
    enddo

    do ix = 1, nbx
      do i = 0, nl-1
        if (L2_cloud%ctp(ix, i) == 0.0 ) then
          L2_cloud%qflags(ix, i) = 10  ! Bad results (should not be used)
        else
          L2_cloud%qflags(ix, i)  = 0  ! Good results
        endif
      enddo
    enddo

    deallocate(flgbits, cfr, ctp, qflag)
    RETURN
  end subroutine read_cloud_tio



  !>Open L2 netCDF cloud file and get dimensions
  !---------------------------------------------------------------------
  !
  !> @param[in] l2file filename for L2 netCDF cloud file
  !> @param tio_l2obj L2 cloud file object
  !> @param[out] nstep_loc mirror step dimension size in L2 cloud file
  !> @param[out] nxtrack_loc xtrack dimension size in L2 cloud file
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan June 2016
  !---------------------------------------------------------------------
  subroutine read_cloud_dims(l2file, tio_l2obj, nstep_loc, nxtrack_loc, errstat)

    implicit none

    !input variables
    character (len=*), intent (in) :: l2file 

    !output variables
    integer (kind=4), intent(out) :: nstep_loc, nxtrack_loc
    integer (kind=4), intent(inout) :: errstat

    type (tiof_file_type) :: tio_l2obj

    if (errstat /= 0) return

    call tiof_open (l2file, tio_l2obj, nf90_nowrite, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_xtrack, nxtrack_loc, errstat)
    call tiof_inq_dimlen (tio_l2obj, o3p_dim_step, nstep_loc, errstat)
    call tiof_close (tio_l2obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_open_error, &
           "read_cloud_dims: failed to open L2 cloud file", &
           errstat)
      return
    endif

  end subroutine read_cloud_dims



  !>Open L2 netCDF cloud file and get dimensions
  !---------------------------------------------------------------------
  !
  !> @param[in] l2file filename for L2 netCDF cloud file
  !> @param tio_l2obj L2 cloud file object
  !> @param[in] nstep number of mirror steps to read data from
  !> @param[in] nxtrack number of xtrack positions to read data from
  !> @param[out] cfr cloud fraction values read from L2 cloud file
  !> @param[out] ctp cloud pressure values read from L2 cloud file
  !> @param[out] qflag cloud processing quality values read from L2 cloud file
  !> @param errstat error handling integer, non-zero indicates failure
  !
  !> @author E. O'Sullivan June 2016
  !---------------------------------------------------------------------
  subroutine read_cloud_data (l2file, tio_l2obj, nstep, nxtrack, cfr, ctp, &
       qflag, errstat)

    implicit none

    !input variables
    character (len=*), intent (in) :: l2file 
    integer (kind=4), intent (in) :: nstep, nxtrack

    !output variables
    real    (kind=4), dimension (nxtrack, 0:nstep-1), intent(out) :: cfr, ctp
    integer (kind=2), dimension (nxtrack, 0:nstep-1), intent(out) :: qflag
    integer (kind=4), intent(inout) :: errstat

    type (tiof_file_type) :: tio_l2obj

    if (errstat /=0) return

    call tiof_open (l2file, tio_l2obj, nf90_nowrite, errstat)
!    call tiof_inq_group (tio_l2obj, "/product", errstat)
    call tiof_push_group (tio_l2obj, o3p_grp_product, errstat)
    call tiof_get2d_r4 (tio_l2obj, cld_var_cld_frac, [0,0], [nstep, nxtrack], &
         cfr(1:nxtrack, 0:nstep-1), errstat)
    call tiof_get2d_r4 (tio_l2obj, cld_var_cld_pres, [0,0], [nstep, nxtrack], &
         ctp(1:nxtrack, 0:nstep-1), errstat)
!rr    call tiof_pop_group (tio_l2obj, errstat)
!rr    call tiof_push_group (tio_l2obj, o3p_grp_qa_stats, errstat)
!    call tiof_inq_group (tio_l2obj, "/", errstat)
!    call tiof_inq_group (tio_l2obj, o3p_grp_qa_stats, errstat)
    call tiof_get2d_i2 (tio_l2obj, cld_var_qflag, [0,0], [nstep, nxtrack], &
         qflag(1:nxtrack, 0:nstep-1), errstat)
    call tiof_close (tio_l2obj, errstat)

    if (errstat /= 0) then
      call tell_error (tell_io_read_error, "read_coud_data: failed", errstat)
      return
    endif

  end subroutine read_cloud_data



  subroutine fill_in_tempoctp(nxtrack, nstep, ctp, qflag)

    implicit none

    ! input variables
    integer (kind=4), intent (in) :: nstep, nxtrack

    ! output variables
    integer (kind=2), dimension (nxtrack, 0:nstep-1), intent(inout) :: qflag
    real    (kind=4), dimension (nxtrack, 0:nstep-1), intent(inout) :: ctp

    ! local variables
    integer (kind=4) :: ix, i, j, fidx, lidx, sidx, eidx
    real    (kind=4) :: frac, dis

    ! If cloud top pressure is too small, then flag that retrieval
    where (ctp < 90.0 ) 
      qflag = int(qflag + 1 , kind=2)
    end where

    ! If cloud top pressure is too large, then reset it to 1013 mb 
    ! (i.e., surface pressure)
    where (ctp > 1013.0 ) 
      ctp = 1013.0
    end where

    ! Reset all flagged pixels to zero 
    where (qflag >= 1) 
      ctp = 0.0
    end where
    !PRINT *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nxtrack * nstep * 1.0)

    ! Fill in flagged pixels
    do i = 1, nstep - 2   ! Fill in along the track
      do ix = 1, nxtrack 
        if ( ctp(ix, i) == 0.0 .and. qflag(ix, i-1) == 0 .and. &
             qflag(ix, i+1) == 0 ) &
             ctp(ix, i) = ( ctp(ix, i-1) + ctp(ix, i+1) ) / 2.0
      enddo
    enddo

    do i = 0, nstep - 1   ! Fill in across the track 
      do ix = 2, nxtrack - 1 
        if ( ctp(ix, i) == 0.0 .and. qflag(ix - 1, i) == 0 .and. &
             qflag(ix + 1, i) == 0 ) &
             ctp(ix, i) = ( ctp(ix - 1, i) + ctp(ix + 1, i) ) / 2.0
      enddo
    enddo
    !PRINT *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nxtrack * nstep * 1.0)

    ! Linear interpolation along the track
    do ix = 1, nxtrack      
      do i = 0, nstep - 1
        if (ctp(ix, i) > 0.0) exit
      enddo
      fidx = i

      do i = nstep-1, 0, -1
        if (ctp(ix, i) > 0.0) exit
      enddo
      lidx = i

      if (fidx >= lidx ) cycle

      i = fidx + 1
      do while ( i <= lidx )
        if (ctp(ix, i) == 0.0 ) then
          sidx = i - 1
          i = i + 1

          eidx = sidx - 1
          do while (i <= lidx) 
            if ( ctp(ix, i) > 0.0 ) then
              eidx = i
              i = i + 1
              exit
            else
              i = i + 1
            endif
          enddo

          dis = real((eidx - sidx), kind=4)
          if (dis <= 12) then
            do j = sidx + 1, eidx - 1
              frac = 1.0 - real((j - sidx), kind=4) / dis
              ctp(ix, j) = frac * ctp(ix, sidx) + (1.0 - frac) * ctp(ix, eidx)
            enddo
          endif
        else
          i = i + 1
        endif
      enddo
    enddo
    !print *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nxtrack * nstep * 1.0)

    ! Linear interpolation across the track
    do i = 0, nstep - 1     
      do ix = 1, nxtrack
        if (ctp(ix, i) > 0.0) exit
      enddo
      fidx = ix

      do ix = nxtrack, 1, -1
        if (ctp(ix, i) > 0.0) exit
      enddo
      lidx = ix

      if (fidx >= lidx ) cycle

      ix = fidx + 1
      do while ( ix <= lidx )
        if (ctp(ix, i) == 0.0 ) then
          sidx = ix - 1
          ix = ix + 1

          eidx = sidx  - 1
          do while (ix <= lidx) 
            if ( ctp(ix, i) > 0.0 ) then
              eidx = ix
              ix = ix + 1
              exit
            else
              ix = ix + 1
            endif
          enddo

          dis = real((eidx - sidx), kind=4)
          if (dis <= 6) then
            do j = sidx + 1, eidx - 1
              frac = 1.0 - real((j - sidx), kind=4) / dis
              ctp(j, i) = frac * ctp(sidx, i) + (1.0 - frac) * ctp(eidx, i)
            enddo
          endif
        else
          ix = ix + 1
        endif
      enddo
    enddo
    !print *, 1.D0 * COUNT(ctp == 0.0) / (1.D0 * nxtrack * nstep)  

    return

  end subroutine fill_in_tempoctp

end module m_read_cloud_tio
