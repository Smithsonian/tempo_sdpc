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

    USE m_avg
    USE m_sortind
    USE m_eigen
    USE m_sigma
    USE m_matmul
    USE m_find
    USE hdfeos4_parameters
    USE L1B_Reader_class
    USE m_strpos
    USE m_vars, ONLY: cloud_mask, smpx_mean, smpx_stddev, smpx_wavel, &
         stddev_thresh, filename, input_data_path, iprt, smpx_nPix, &
         npixels, thresholds, fill_value, noret, filename_cm
    USE m_swathnames
    USE m_pgs_include

    IMPLICIT NONE

    !   INCLUDE 'PGS_PC.f'
    !   INCLUDE 'PGS_PC_9.f'
    !   INCLUDE 'PGS_SMF.f'
    !   INCLUDE 'PGS_OMI_1900.f'

    INTEGER :: nPixold=2
    INTEGER (KIND = 4), PARAMETER :: maxCoadd=5

    ! declaration of variables used in both examples
    INTEGER (KIND = 4) :: status, ierr, &
         nTimes, nXtrack, &
         iLine, nTimesSmPx
    CHARACTER (LEN = 200) :: filenamen, swathname
    TYPE (L1b_block_type) :: blk, blk2
    INTEGER (KIND = 4), PARAMETER :: zero = 0
    INTEGER (KIND = 2) :: nPix, nPix2
    CHARACTER (LEN = 8) :: fmt101="(5e12.3)",fmt102="(5f12.3)"


    ! declaration of variables specific to Example 1
    INTEGER (KIND = 4) :: itrack

    ! declaration of variables specific to Example 2
    REAL (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: wavelengthL, wavelengthL2
    REAL (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: smvaluesL, smvaluesL2
    integer :: indc, ind
    real (KIND=8), dimension(maxcoadd,2) :: r1_ri
    real (KIND=4), dimension(maxcoadd) :: D, E 
    integer, dimension(maxcoadd) :: order, sorted
    real (KIND=4), dimension(maxcoadd,maxcoadd) :: A
    integer :: ndim
    logical :: using_2_channels=.false.

    ! obtain name of swath
    filenamen=trim(input_data_path)//filename
    vis  = strpos (filename, 'BRVG') > 0
    visz = strpos (filename, 'BRVZ') > 0
    uvsz = strpos (filename, 'BRUZ') > 0
    if (visz) then
      swathname = visswathz
    else if (uvsz) then
      swathname = uv2swathz
    else if (vis) then
      swathname = visswath
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

    if (using_2_channels) then
      filenamen=trim(input_data_path)//filename_cm
      vis = strpos (filename_cm, 'BRUG') < 0 .and. strpos (filename_cm, 'BRUZ') < 0
      if (vis) then
        swathname = "Earth VIS Swath"
      else
        swathname = "Earth UV-2 Swath"
      endif
      if (iprt >= 2) print *,'cloud_mask: filename ',filenamen

      ! open data block structure with default size of 1 lines
      status = L1Br_open( blk2, filenamen, swathname)
      IF( status .NE. OMI_S_SUCCESS ) THEN
        ierr = OMI_SMF_setmsg( status, &
             "L1Br_open2 failed.", "cloud_mask", 0 )
        STOP
      END IF
    endif


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
      if (using_2_channels) then
        status = L1Br_getDATAline( blk2, iLine, nPix2, &
             Data_k=smvaluesL2, &
             Wavelength_k=wavelengthL2)!, quality_flagL )
        IF( status .NE. OMI_S_SUCCESS ) THEN
          ierr = OMI_SMF_setmsg( OMI_E_FAILURE, &
               "L1Brd_getSIGline2 failed", "cloud_mask", 0 )
          call exit(1)
        END IF
      endif

      ! science data processing

      do iTrack=1, nXtrack
        !if (iprt >= 7) print *,'cloud_mask: proc. data from track ',iTrack
        smpx_nPix(iTrack,iLine+1)=nPix
        smpx_wavel(iTrack,iLine+1)=avg(wavelengthL(iTrack,1:nPix))
        indc=count(npixels .eq. nPix)
        if (indc .eq. 1) then
          smpx_mean(iTrack,iLine+1)=sum(smvaluesL(iTrack,1:nPix))/nPix

          if (smpx_mean(iTrack,iLine+1) /= 0.) then

            if (using_2_channels) then

              if (nPix > 2) then
                sorted(1:nPix)=smvaluesL(iTrack,1:nPix)
                order(1:nPix)=i_sortind(sorted(1:nPix))
                sorted(1:nPix)=smvaluesL(iTrack,order(1:nPix))
                r1_ri(1:nPix-1,1)=(sorted(1)-sorted(2:nPix))/sorted(1)
                sorted(1:nPix)=smvaluesL2(iTrack,order(1:nPix))
                r1_ri(1:nPix-1,2)=(sorted(1)-sorted(2:nPix))/sorted(1)
                A=r1_ri .mm. transpose(r1_ri)
                ndim=nPix-1
                call TRED2(A,ndim,maxcoadd,D,E)
                call TQLI(D,E,ndim,maxcoadd,A)
                call EIGSRT(D,A,ndim,maxcoadd)
                smpx_stddev(iTrack,iLine+1) = d(1)/sum(d(1:npix))*100.
                smpx_mean(iTrack,iLine+1) = d(1)
                if (iprt >= 8) then
                  print *, iTrack, smvaluesL(iTrack,order(1))
                  write(6,fmt101) r1_ri(1:npix-1,1)/smvaluesL(iTrack,order(1))
                  write(6,fmt101) r1_ri(1:npix-1,2)/smvaluesL2(iTrack,order(1))
                  write(6,fmt102) d(1:npix-1)/sum(d(1:npix))*100. 
                  write(6,fmt101) A(1:npix-1,1)
                  write(6,fmt101) A(1:npix-1,2)
                endif ! iprt 

              endif ! npix > 2

            else ! use only vis
              smpx_stddev(iTrack,iLine+1)=abs(smvaluesL(iTrack,1)- &
                   smvaluesL(iTrack,nPix)) / &
                   smpx_mean(iTrack,iLine+1)

            endif ! use 2 channels
          else ! mean was zero
            smpx_stddev(iTrack,iLine+1)=0.
          endif
          ind=find1(npixels .eq. nPix)
          if (smpx_stddev(iTrack,iLine+1) >= thresholds(ind,iTrack)) then
            cloud_mask(iTrack,iLine+1)=1
          elseif (smpx_stddev(iTrack,iLine+1) > 0) then
            cloud_mask(iTrack,iLine+1)=0
          endif
        endif ! if npix found
      enddo ! loop over Xtrack

      if (iLine /= 0 .and. nPix /= nPixold) cloud_mask(:,iLine+1)=3
      nPixold=nPix

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
