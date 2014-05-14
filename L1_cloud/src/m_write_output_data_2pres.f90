module m_write_output_data_2pres

contains

  subroutine write_output_data_2pres(outfile, swathname)
    use m_write_swath_field
    use m_vars
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  write_output_data
    ! 
    ! !DESCRIPTION: write_output_data writes OMI cloud level 2 data to HDF-EOS
    !		file
    !
    ! !CALLING SEQUENCE: 
    !
    !        call write_output_data
    !     
    ! !INPUT PARAMETERS:   
    !
    ! !OUTPUT PARAMETERS:  
    !
    ! !SEE ALSO:  
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner     original fortran 90
    !  15Jan03   Vasilkov   write GroundpixelQualityFlags
    !  30Aug03   Vasilkov   write MeasurementQualityFlags
    !
    !EOP
    !-------------------------------------------------------------------------
    !
    !inputs
    !------
    character(len=*), intent(in) :: outfile, swathname

    !ouputs
    !------
    !       integer, intent(out) :: ierr

    ! Declare the HDF-EOS file and swath identification numbers, and
    ! the status of the HDF-EOS functions calls.
    !-----------------------------------------------------------------
    integer (kind = 4) swfid, swid

    ! Declare the HDF-EOS functions.
    !-------------------------------
    !       integer (kind = 4) :: he5_swcreate, he5_swdefdim
    integer, parameter :: HE5_ACC_RDWR=100
    !       integer, parameter :: HE5_ACC_TRUNC=102
    integer (kind = 4) :: he5_swopen, he5_swattach, he5_swdetach, he5_swclose

    !Local variables
    !---------------
    character(len=255) :: nTimesstr="nTimes"
    character(len=255) :: nXtrackstr="nXtrack"
    character(len=255) :: dims2
    integer            :: nTime, iLine1
    !      integer            :: nWaveRes
    real (kind = 4) :: misval_r4 = fill_value
    integer (kind = 2) :: misval_i2 = fill_value_int 

    iLine1=iLine-start_line
    nTime=1
    dims2 = trim(nXtrackstr)//','//trim(nTimesstr)

    ! Open the OMI Level 2 HDF-EOS cloud output file.
    !-------------------------------------------------
    swfid = he5_swopen (outfile, HE5_ACC_RDWR)
    ! Attach to the swath
    !----------------------
    swid = he5_swattach(swfid, swathname)

    ! Write the fields
    !----------------------
    !      cld_pres2=cld_pres2*1013.25
    status = put_data (swid, "CloudPressureforO3", &
         dims2, cld_pres2, misval_r4, "Cloud Pressure for O3", "hPa", &
         offset=(/0,0/),iprt=iprt)
    status = put_data (swid, "CloudFractionforO3", &
         dims2, eff_cld_frac2, misval_r4, "Cloud Fraction for O3", "NoUnits", &
         offset=(/0,0/), iprt=iprt)
    status = put_data (swid, "ProcessingQualityFlagsforO3", &
         dims2, qc2, misval_i2, "Processing Quality Flags for O3", "NoUnits", &
         offset=(/0,0/), iprt=iprt)

    ! Detach from the swath interface.
    !-------------------------------------------------
    status = he5_swdetach (swid)
    if (iprt >= 2) &
         write (6, *) 'write_output_data: detached swath ',status

    ! Close the OMI Level2  HDF-EOS output file.
    !----------------------------------------------
    status = he5_swclose (swfid)
    if (iprt >= 2) &
         write (6, *) 'write_output_data: closed file ', status

  end subroutine write_output_data_2pres

end module m_write_output_data_2pres
