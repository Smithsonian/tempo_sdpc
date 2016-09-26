!>Read in small pixel data and allocate memory preparatory to 
!>calculating cloud mask
!
!-------------------------------------------------------------------------
!
! subroutine cld_mask:
!
!>  This subroutine reads in small pixel data
!>  SmPx_reader_class module to obtain the geolocation and data fields from 
!>  an OMI L1B data file in HE4 format and computes a cloud mask based on 
!>  spatial homogeneity\n
!>
!>  SmPx Reader Functions Used:\n
!>
!>   L1Br_open        - opens data block structure\n
!>   L1Br_getSWdims   - gets sizes of dimensions defined in swath\n
!>   L1Br_getGEOline  - gets geolocation information\n
!>   L1Br_getDATAline - gets data fields for one "scan" line\n
!>   L1Br_close       - closes data block structure\n
!>
!>  Inputs:\n
!>
!>   blk        - data block structure\n
!>   filename   - OMI L1B data file name (see PCF file)\n
!>   swathname  - name of swath ("Sun Daily VIS Swath" in this example)\n
!>   iLine      - "scan" line number\n
!>
!>  Outputs:\n
!>
!>    tim     - time\n
!>    nPix    - number of small pixel columns in this scan line\n
!>    valname - small pixel value name -- radiance, irradiance, or signal\n
!>
!>   Geolocation fields at center of each ground pixel for each "scan" line:\n
!>    latitude  - latitudes\n
!>    longitude - longitudes\n
!>    szenith   - solar zenith angles\n
!>    sazimuth  - solar azimuth angles\n
!>    vzenith   - view zenith angles\n
!>    vazimuth  - view azimuth angles\n
!>    height    - terrain heights\n
!>    geoflg    - geolocation flags\n
!>
!>    mflg - measurement quality flag\n
!>
!>   Data fields at center of each ground pixel for each "scan" line:\n
!>    smvaluesL     - small pixel values, either irradiance, radiance, 
!>                     or signal\n
!>    quality_flagL - pixel quality flags\n
!>    wavelengthL   - small pixel wavelength\n
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
!> @author  Jeremy Warner, SSAI, 25 April 2003, original driver
!> @author Joanna Joiner, GSFC  08 July  2003, modifications for cloud mask
!
!!-------------------------------------------------------------------------
!!-------------------------------------------------------------------
!
!> @author Ewan O'Sullivan 27Jul14 Updated for TEMPO, separated processing
!! step into m_cloud_mask_proc
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
module m_cloud_mask

contains

  subroutine cld_mask (errstat)

    use m_cloud_mask_proc
    use hdfeos4_parameters
    use L1B_Reader_class
    use m_strpos
    use m_vars, only: cloud_mask, smpx_mean, smpx_stddev, smpx_wavel, &
         filename, input_data_path, smpx_nPix, fill_value
    use m_swathnames
    use m_pgs_include
    use tell_module

    implicit none

    integer, intent (inout) :: errstat

    !Local variables
    integer (KIND = 4), parameter :: maxCoadd=5
    integer (KIND = 4) :: status, nTimes, nXtrack, iLine, nTimesSmPx
    character (LEN = 200) :: filenamen, swathname, logmsg
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
    write(logmsg,"(A21,A)") 'cloud_mask: filename ',trim(filenamen)
    call tell_log(2,logmsg)

    ! open data block structure with default size of 1 lines
    status = L1Br_open( blk, filenamen, swathname)!, valname )
    if( status .ne. OMI_S_SUCCESS ) then
      errstat = -1
      call tell_error(tell_io_open_error,"cloud_mask: L1Br_open failed", &
           errstat)
      return
    end if


    ! obtain sizes of dimensions defined in swath
    status = L1Br_getSWdims( blk,  NumTimes_k=nTimes, nXtrack_k=nXtrack, &
         NumTimesSmallPixel_k=nTimesSmPx )
    if( status .ne. OMI_S_SUCCESS ) then
      errstat = -1
      call tell_error(tell_io_read_error,"cloud_mask: L1Br_getSWdims failed", &
           errstat)
      return
    else 
      call tell_log(3,'cloud_mask: nTimes, nXtrack, nTimesSmPx')
      write(logmsg,"(3I6)") nTimes, nXtrack, nTimesSmPx
      call tell_log(3,logmsg)
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

      write(logmsg,"(A35,I4)") 'cloud_mask: reading data from line ',iLine
      call tell_log(3,logmsg)
      status = L1Br_getDATAline( blk, iLine, nPix, &
           Data_k=smvaluesL, &
           Wavelength_k=wavelengthL)!, quality_flagL )
      if( status .ne. OMI_S_SUCCESS ) then
        errstat = -1
        call tell_error(tell_io_read_error, &
             "cloud_mask: L1Br_getDATAline failed", errstat)
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
      errstat = -1
      call tell_error(tell_io_error,"cloud_mask: L1Br_close failed", errstat)
      return
    end if
    call tell_log(1,'cloud_mask finished successfully')

  end subroutine cld_mask

end module m_cloud_mask
