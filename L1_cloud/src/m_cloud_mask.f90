module m_cloud_mask

contains

  subroutine cld_mask (errstat)
    !!-------------------------------------------------------------------------
    !
    ! subtoutine cld_mask:
    !
    !  This subroutine reads in small pixel data
    !  SmPx_reader_class module to obtain the geolocation and data fields from 
    !  an OMI L1B data file in HE4 format and computes a cloud mask based on 
    !  spatial homogeneity
    !
    !  SmPx Reader Functions Used:
    !
    !   L1Br_open        - opens data block structure
    !   L1Br_getSWdims   - gets sizes of dimensions defined in swath
    !   L1Br_getGEOline  - gets geolocation information
    !   L1Br_getDATAline - gets data fields for one "scan" line
    !   L1Br_close       - closes data block structure
    !
    !  Inputs:
    !
    !   blk        - data block structure
    !   filename   - OMI L1B data file name (see PCF file)
    !   swathname  - name of swath ("Sun Daily VIS Swath" in this example)
    !   iLine      - "scan" line number
    !
    !  Outputs:
    !
    !    tim     - time
    !    nPix    - number of small pixel columns in this scan line
    !    valname - small pixel value name -- radiance, irradiance, or signal
    !
    !   Geolocation fields at center of each ground pixel for each "scan" line:
    !    latitude  - latitudes
    !    longitude - longitudes
    !    szenith   - solar zenith angles
    !    sazimuth  - solar azimuth angles
    !    vzenith   - view zenith angles
    !    vazimuth  - view azimuth angles
    !    height    - terrain heights
    !    geoflg    - geolocation flags
    !
    !    mflg - measurement quality flag
    !
    !   Data fields at center of each ground pixel for each "scan" line:
    !    smvaluesL     - small pixel values, either irradiance, radiance, or signal
    !    quality_flagL - pixel quality flags
    !    wavelengthL   - small pixel wavelength
    !
    ! Function Calls Using Keywords:
    !
    !  The use of Fortran 90 keywords allows subsets of the geolocation and 
    !  data fields to be obtained. For example, if just the small pixel values 
    !  and wavelengths for one "scan" line are needed, then the L1Br_getSIGline
    !  function can be called as follows to obtain these fields:
    !
    !   status = L1Br_getDATAline( blk, iLine, nPix, &
    !                             Data_k=smvaluesL, &
    !                             Wavelength_k=wavelengthL )
    !
    ! Author:
    !
    !  Jeremy Warner, SSAI, 25 April 2003, original driver
    !  Joanna Joiner, GSFC  08 July  2003, modifications for cloud mask
    !
    !!-------------------------------------------------------------------------

    !!-------------------------------------------------------------------
    !
    ! Ewan O'Sullivan 27Jul14
    !
    ! Commenting out all material relating to using_2_channels
    !   It appears the code was originally designed to have the option to
    !   use data from both the VIS and UV channels, but was then put into 
    !   OMCLDRR which only uses UV. Possible the same code exists in 
    !   OMCLDO2 and there uses only VIS? In any case, using_2_channels
    !   code uses Numerical Recipes routines, and appears to be never called
    !   since using_2_channels=.false.
    !
    ! 6Aug14
    !
    ! Separated processing step out into m_cloud_mask_proc. This module
    !   now performs only I/O functions, hopefully making it easier to 
    !   rewrite for TEMPO
    !
    !!--------------------------------------------------------------------

    use m_cloud_mask_proc
    use hdfeos4_parameters
    use L1B_Reader_class
    use m_strpos
    use m_vars, only: cloud_mask, smpx_mean, smpx_stddev, smpx_wavel, &
         filename, input_data_path, iprt, smpx_nPix, fill_value
    use m_swathnames
    use m_pgs_include
    use tell_module

    implicit none

    integer, intent (inout) :: errstat

    !Local variables
    integer (KIND = 4), parameter :: maxCoadd=5
    integer (KIND = 4) :: status, ierr, nTimes, nXtrack, iLine, nTimesSmPx
    character (LEN = 200) :: filenamen, swathname
    type (L1b_block_type) :: blk 
    integer (KIND = 2) :: nPix 
    real (KIND = 4), dimension(:,:), allocatable :: wavelengthL, wavelengthL2
    real (KIND = 4), dimension(:,:), allocatable :: smvaluesL, smvaluesL2

    if (errstat /= 0) return

    ! obtain name of swath
    filenamen=trim(input_data_path)//filename
    uvsz = strpos (filename, 'BRUZ') > 0

    if (uvsz) then
      swathname = uv2swathz
    else
      swathname = uv2swath
    endif
    if (iprt >= 2) print *,'cloud_mask: filename ',filenamen

    ! open data block structure with default size of 1 lines
    status = L1Br_open( blk, filenamen, swathname)!, valname )
    if( status .ne. OMI_S_SUCCESS ) then
      ierr = OMI_SMF_setmsg( status, &
           "L1Br_open failed,"//trim(filenamen), "cloud_mask", 0 )
      stop
    end if


    ! obtain sizes of dimensions defined in swath
    status = L1Br_getSWdims( blk,  NumTimes_k=nTimes, nXtrack_k=nXtrack, &
         NumTimesSmallPixel_k=nTimesSmPx )
    if( status .ne. OMI_S_SUCCESS ) then
      ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
           "L1Br_getSWdims failed.", "cloud_mask", 0 )
      errstat = -1
      return
    else if (iprt > 2) then
      print *,' cloud_mask: nTimes, nXtrack, nTimesSmPx', &
           nTimes, nXtrack, nTimesSmPx
    end if

    !Allocate arrays
    allocate(smvaluesL(nXtrack,maxCoadd), &
         smvaluesL2(nXtrack,maxCoadd), &
         wavelengthL(nXtrack,maxCoadd), &
         wavelengthL2(nXtrack,maxCoadd), &
         cloud_mask(nXtrack,nTimes), &
         smpx_mean(nXtrack,nTimes), &
         smpx_wavel(nXtrack,nTimes), &
         smpx_stddev(nXtrack,nTimes), &
         smpx_nPix(nXtrack,nTimes), stat=errstat )
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "cloud_mask: failed to allocate memory", &
           errstat)
      return
    endif
    cloud_mask=2
    smpx_mean=fill_value
    smpx_stddev=fill_value
    smpx_nPix=0

    !Loop over all scan lines
    do iLine = 0, nTimes-1

      if (iprt >= 7) print *,'cloud_mask: reading data from line ',iLine
      status = L1Br_getDATAline( blk, iLine, nPix, &
           Data_k=smvaluesL, &
           Wavelength_k=wavelengthL)!, quality_flagL )
      if( status .ne. OMI_S_SUCCESS ) then
        ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
             "L1Brd_getSIGline failed", "cloud_mask", 0 )
        errstat = -1
        return
      end if

      !Science data processing
      call cloud_mask_proc(nXtrack, nPix, iLine, maxCoadd, &
           wavelengthL, smvaluesL, errstat)


    end do ! loop over iLine

    ! deallocate memory
    deallocate (smvaluesL,wavelengthL, stat=errstat)
    if (errstat /= 0) then
      call tell_error (tell_malloc_error, &
           "cloud_mask: failed to deallocate memory", &
           errstat)
      return
    endif

    ! close data block structure
    status = L1Br_close( blk )
    if( status .ne. OMI_S_SUCCESS ) then
      ierr = OMI_SMF_setmsg( status, &
           "L1Br_close failed.", "cloud_mask", 0 )
      errstat = -1
      return
    end if
    ierr = OMI_SMF_setmsg(PGS_S_SUCCESS, "Test Done", "cld_mask, m_cloud_mask", 0 )

  end subroutine cld_mask

end module m_cloud_mask
