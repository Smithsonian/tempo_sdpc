!Determine cloud status based on std.dev of small pixel values
module m_cloud_mask_proc

contains

  subroutine cloud_mask_proc(nXtrack, nPix, iLine, maxCoadd, &
       wavelengthL, smvaluesL)
    !!------------------------------------------------------------------
    !
    ! Carries out data processing stage of cloud mask generation,
    ! calculating std. dev. of small pixel values and comparing them
    ! to preset thresholds.
    !
    ! Taken out of m_cloud Mask to separate processing from I/O
    !
    ! Called by m_cloud_mask
    !
    ! INPUT VARIABLES
    ! maxCoadd: maximum number of small pixels to coadd, =5 in m_cloud_mask
    ! nXtrack: Number of cross-track pixels, read in from data file by 
    !           m_cloud_mask
    ! iLine: along-track "scan" line number currently under consideration, 
    !         m_cloud_mask_proc is called once for each line
    ! nPix: number of small pixel columns in this scan line
    ! wavelengthL: small pixel wavelength
    ! smvaluesL: small pixel values
    !
    ! INPUT/OUTPUT VARIABLES FROM M_VARS
    ! smpx_nPix: number of small pixels in each full-size pixel
    ! smpx_wavel: average of small pixel wavelengths in each full pixel
    ! smpx_mean: mean of smvaluesL in each full-size pixel
    ! smpx_stddev: standard deviataion of smvaluesL in each full pixel
    ! cloud_mask: mask value (0=clear, 1=cloud, 2=no data, 3=nPix changed)
    ! thresholds: Standard deviation thresholds to compare with smpx_stddev
    !             pixels above threshold contain cloud
    ! npixels: valid nPix values for which thresholds exist, read from
    !          header of threshold file by m_read_thresholds, =2,4,5
    !
    ! LOCAL VARIABLES
    ! iTrack: cross-track pixel index
    ! indc: is nPix value valid? 1=yes, continue with processing
    ! ind: index for accessing thresholds values.
    ! nPixold: used to flag pixels where nPix changes (it should be const)
    !
    ! Author: O'Sullivan, 05 August 2014
    !
    !!------------------------------------------------------------------

    use m_avg
    use m_find
    use m_vars, only: smpx_nPix, smpx_wavel, smpx_mean, smpx_stddev, &
         cloud_mask, thresholds, npixels

    implicit none

    !input variables
    integer (KIND = 4), intent(in) :: maxCoadd, nXtrack, iLine
    integer (KIND = 2), intent(in) :: nPix
    real (KIND = 4), dimension(nXtrack,maxCoadd), intent(in) :: wavelengthL, &
         smvaluesL
    !local variables
    integer (KIND = 4) :: ind, indc, iTrack
    integer (KIND = 2) :: nPixold=2

    !Data processing loop
    do iTrack=1, nXtrack
      smpx_nPix(iTrack,iLine+1)=nPix
      smpx_wavel(iTrack,iLine+1)=avg(wavelengthL(iTrack,1:nPix))
      indc=count(npixels .eq. nPix)
      if (indc .eq. 1) then
        smpx_mean(iTrack,iLine+1)=sum(smvaluesL(iTrack,1:nPix))/nPix
        if (smpx_mean(iTrack,iLine+1) /= 0.) then
          smpx_stddev(iTrack,iLine+1)=abs(smvaluesL(iTrack,1)- &
               smvaluesL(iTrack,nPix)) / &
               smpx_mean(iTrack,iLine+1)
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


  end subroutine cloud_mask_proc

end module m_cloud_mask_proc
