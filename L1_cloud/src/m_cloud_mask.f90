module m_cloud_mask

contains

  subroutine cld_mask
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

    USE m_cloud_mask_proc
    USE hdfeos4_parameters
    USE L1B_Reader_class
    USE m_strpos
    USE m_vars, ONLY: cloud_mask, smpx_mean, smpx_stddev, smpx_wavel, &
         stddev_thresh, filename, input_data_path, iprt, smpx_nPix, &
         npixels, thresholds, fill_value, noret, filename_cm
    USE m_swathnames
    USE m_pgs_include

    IMPLICIT NONE

    !Local variables
    INTEGER (KIND = 4), PARAMETER :: maxCoadd=5
    INTEGER (KIND = 4) :: status, ierr, nTimes, nXtrack, iLine, nTimesSmPx
    CHARACTER (LEN = 200) :: filenamen, swathname
    TYPE (L1b_block_type) :: blk 
    INTEGER (KIND = 4), PARAMETER :: zero = 0
    INTEGER (KIND = 2) :: nPix 
    REAL (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: wavelengthL, wavelengthL2
    REAL (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: smvaluesL, smvaluesL2

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
    IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( status, &
           "L1Br_open failed,"//trim(filenamen), "cloud_mask", 0 )
      STOP
    END IF


    ! obtain sizes of dimensions defined in swath
    status = L1Br_getSWdims( blk,  NumTimes_k=nTimes, nXtrack_k=nXtrack, &
         NumTimesSmallPixel_k=nTimesSmPx )
    IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
           "L1Br_getSWdims failed.", "cloud_mask", 0 )
      call exit(1)
    else if (iprt > 2) then
      print *,' cloud_mask: nTimes, nXtrack, nTimesSmPx', &
           nTimes, nXtrack, nTimesSmPx
    END IF

    ALLOCATE( smvaluesL(nXtrack,maxCoadd), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "smvaluesL allocation failure", "cloud_mask", 0 )
      call exit(1)
    END IF

    ALLOCATE( smvaluesL2(nXtrack,maxCoadd), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "smvaluesL2 allocation failure", "cloud_mask", 0 )
      call exit(1)
    END IF

    ALLOCATE( wavelengthL(nXtrack,maxCoadd), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "wavelengthL allocation failure", "cloud_mask", 0 )
      call exit(1)
    END IF

    ALLOCATE( wavelengthL2(nXtrack,maxCoadd), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "wavelengthL2 allocation failure", "cloud_mask", 0 )
      call exit(1)
    END IF

    ALLOCATE( cloud_mask(nXtrack,nTimes), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "cloud_mask allocation failure", "cloud_mask", 0 )
      call exit(1)
    ELSE
      cloud_mask=2
    END IF

    ALLOCATE( smpx_mean(nXtrack,nTimes), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "smpx_mean allocation failure", "cloud_mask", 0 )
      call exit(1)
    ELSE
      smpx_mean=fill_value
    END IF

    ALLOCATE( smpx_stddev(nXtrack,nTimes), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "smpx_stddev allocation failure", "cloud_mask", 0 )
      call exit(1)
    ELSE
      smpx_stddev=fill_value
    END IF

    ALLOCATE( smpx_nPix(nXtrack,nTimes), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "smpx_nPix allocation failure", "cloud_mask", 0 )
      call exit(1)
    ELSE
      smpx_nPix=0
    END IF

    ALLOCATE( smpx_wavel(nXtrack,nTimes), STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "smpx_wavel allocation failure", "cloud_mask", 0 )
      call exit(1)
    END IF

    !Loop over all scan lines
    DO iLine = 0, nTimes-1

      if (iprt >= 7) print *,'cloud_mask: reading data from line ',iLine
      status = L1Br_getDATAline( blk, iLine, nPix, &
           Data_k=smvaluesL, &
           Wavelength_k=wavelengthL)!, quality_flagL )
      IF( status .NE. OMI_S_SUCCESS ) THEN
        ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
             "L1Brd_getSIGline failed", "cloud_mask", 0 )
        call exit(1)
      END IF

      !Science data processing
      call cloud_mask_proc(nXtrack, nPix, iLine, maxCoadd, &
       wavelengthL, smvaluesL)


    END DO ! loop over iLine

    ! deallocate memory

    DEALLOCATE( smvaluesL, STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "smvaluesL deallocation failure", "cloud_mask", 0 )
      call exit(1)
    END IF

    DEALLOCATE( wavelengthL, STAT=ierr )
    IF( ierr .NE. zero ) THEN
      ierr = OMI_SMF_setmsg( OMI_E_MEM_ALLOC, &
           "wavelengthL deallocation failure", "cloud_mask", 0 )
      call exit(1)
    END IF

    ! close data block structure
    status = L1Br_close( blk )
    IF( status .NE. OMI_S_SUCCESS ) THEN
      ierr = OMI_SMF_setmsg( status, &
           "L1Br_close failed.", "cloud_mask", 0 )
      call exit(1)
    END IF
    ierr = OMI_SMF_setmsg(PGS_S_SUCCESS, "Test Done", "cld_mask, m_cloud_mask", 0 )

  end subroutine cld_mask

end module m_cloud_mask
