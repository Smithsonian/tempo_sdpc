!! F90
 ! 
 !!Description:
 ! MODULE L1B_smlpix_class contains the OMIL1B small pixel radiance or 
 ! small pixel irradiance data structure and the functions for the creation, 
 ! deletion and access of specific information stored in the data structure. 
 ! A user can use the following functions:  L1Bsmpx_open, L1Bsmpx_close, 
 ! L1Bsmpx_getSWdims, L1Bsmpx_getLine, L1Bsmpx_getLineWL.
 !
 !!Revision History:
 ! Revision 0.1  08/26/2003  Kai Yang/UMBC
 !!Team-unique Header:
 ! This software was developed by the OMI Science Team Support
 ! Group for the National Aeronautics and Space Administration, Goddard
 ! Space Flight Center, under NASA Task 916-003-1
 !
 !!References and Credits
 ! Written by
 ! Kai Yang
 ! University of Maryland Baltimore County 
 ! email: Kai.Yang-1@nasa.gov
 !
!!

MODULE L1B_smlpix_class
   USE OMI_SMF_class
   USE UTIL_tools_class
   IMPLICIT NONE
   PUBLIC  :: L1Bsmpx_open, L1Bsmpx_close
   PUBLIC  :: L1Bsmpx_getSWdims
   PUBLIC  :: L1Bsmpx_getLine
   PUBLIC  :: L1Bsmpx_getMeasurementQF
   PUBLIC  :: L1Bsmpx_getInstConfigId
   PUBLIC  :: L1Bsmpx_iLineSmPx2LgPx
   PUBLIC  :: L1Bsmpx_getAccuNcols

   PRIVATE :: fill_smlpix_blk
   INTEGER (KIND = 4), PARAMETER, PUBLIC :: MAX_NAME_LENGTH = 255
   INTEGER, PRIVATE :: ierr  !!error code returned from a function
   INTEGER, PARAMETER, PRIVATE :: zero = 0
   INTEGER, PARAMETER, PRIVATE :: one  = 1
   INTEGER, PARAMETER, PRIVATE :: two  = 2
!!
 ! L1B_smlpix_type is designed to store a block of L1B swath data.
 ! A block contains all the small pixel (ir)radiance and wavelength
 ! data fields for all the "scan" lines or "exposure times"
 ! in a L1B swath. 
 ! filename: the file name for the L1B file (including path)
 ! swathname: the swath name in the L1B file, either "Earth UV-1 Swath",
 !            "Earth UV-2 Swath", or "Earth VIS Swath"
 ! nTimes: the total number of lines contained in the L1B swath
 ! nTimesSmallPixel: the total number of lines contained in the small pixel data
 ! nXtrack: the number of cross-track pixels in the swath
 ! nWavel: the number of wavelengths in the swath
 ! initialized: the logical variable to indicate whether the
 !              block has been initialized.  The block data structure
 !              has to be created (initialized) before it can be used.
 !              
 ! Radiance and Irradiance fields:
 !   sxData: SmallPixelRadiance || SmallPixelIrradiance
 !   sxWl: SmallPixelWavelength
 !   PixelQualityFlags
 ! These are two-dimensional data arrays created for the number of
 ! lines and the number of pixels for the data block when the data structure 
 ! is created. 
 !
 ! NumberSmallPixelColumns
 ! SmallPixelColumn
 ! MeasurementQualityFlags
 ! InstrumentConfigurationId
 !
 ! These are one-dimensional data arrays created for the number of
 ! lines for the data block when the data structure is created. 
 !
 ! Again, a user cannot see these data fields directly but the 
 ! information can be retrieved only through functions L1Bsmpx_getLine.
 ! The "line" functions retrieve information for one 
 ! sepcific line defined by iLineSmPx.
!!

   TYPE, PUBLIC :: L1B_smlpix_type
      PRIVATE
      CHARACTER ( LEN = MAX_NAME_LENGTH ) :: filename, swathname
      INTEGER (KIND = 4) :: nTimes
      INTEGER (KIND = 4) :: nTimesSmallPixel
      INTEGER (KIND = 4) :: nXtrack
      INTEGER (KIND = 4) :: nWavel
      LOGICAL :: initialized

      ! Small Pixel Radiance, Irradiance, or Signal fields
      REAL (KIND = 4)   , DIMENSION(:,:), POINTER :: sxData, sxWl
      INTEGER (KIND = 2), DIMENSION(:,:), POINTER :: PixelQualityFlags
      INTEGER (KIND = 1), DIMENSION(:)  , POINTER :: NumberSmallPixelColumns,&
                                                     InstrumentConfigurationId
      INTEGER (KIND = 2), DIMENSION(:)  , POINTER :: SmallPixelColumn, &
                                                     MeasurementQualityFlags
      INTEGER (KIND = 4), DIMENSION(:)  , POINTER :: AccumNsmpxCols
   END TYPE L1B_smlpix_type

   CONTAINS

!! 1. L1Bsmpx_open
 !    This function should be called first to initiate the
 !    the interface with the L1B swath. 
 !    this: the block data structure
 !    fn  : the L1B file name
 !    swn : the swath name in the L1B file
 !    nL  : the number of lines the data block can store.  This is
 !          an optional input.  If it is not present, then a default
 !          value of 100 is used.
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bsmpx_open( this, fn, swn ) RESULT (status)
        USE HE4_class
        TYPE (L1B_smlpix_type), INTENT( INOUT ) :: this
        CHARACTER ( LEN = * ), INTENT( IN ) :: fn, swn
        INTEGER (KIND = 4) :: swfid, swid, status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg
        
        !! open the L1B swath file
        swfid = swopen( fn, DFACC_READ )
        IF( swfid < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_FILE_OPEN, fn, "L1Bsmpx_open", zero )
           RETURN
        ENDIF 

        !! attach to the swath
        swid = swattach( swfid, swn )
        IF( swid < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_SWATH_ATTACH, swn, &
                                 "L1Bsmpx_open", zero )
           RETURN
        ENDIF 
        
        !! retrieve the dimension info from the swath file.
        !! dimension names are obtained from file spec
        status = swrdattr( swid, "NumTimes", this%nTimes )
        IF( status < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get nTimes size failed", &
                                 "L1Bsmpx_open", zero )
           RETURN
        ENDIF 
        IF( this%nTimes <= zero ) THEN
           WRITE( msg,'(A,I3,A)' ) "nTimes =", this%nTimes, &
                 ", No scan lines in "//TRIM(swn)// ", " // TRIM(fn) 
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, "L1Bsmpx_open", zero )
           RETURN
        ENDIF 

        status = swrdattr( swid, "NumTimesSmallPixel", this%nTimesSmallPixel )
        IF( status < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS,  &
                                 "get NumTimesSmallPixel size failed", &
                                 "L1Bsmpx_open", zero )
           RETURN
        ENDIF 
        IF( this%nTimesSmallPixel <= zero ) THEN
           WRITE( msg,'(A,I3,A)' ) "nTimesSmallPixel =", this%nTimesSmallPixel,&
                 ", No scan lines in "//TRIM(swn)// ", " // TRIM(fn) 
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, "L1Bsmpx_open", zero )
           RETURN
        ENDIF 

        this%nXtrack    = swdiminfo( swid, "nXtrack" )
        IF( this%nXtrack < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get xTrack size failed", &
                                  "L1Bsmpx_open", zero )
           RETURN
        ENDIF 

        this%nWavel     = swdiminfo( swid, "nWavel" )
        IF( this%nWavel < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get nWavel size failed", &
                                  "L1Bsmpx_open", zero )
           RETURN
        ENDIF

        !! detach and close the L1B swath files.  No need for 
        !! error checking, for an error is unlikely to occur here.
        ierr = swdetach( swid )
        ierr = swclose( swfid )

        this%filename = fn   
        this%swathname = swn
        
        !! Allocate the memory for storage of radiance fields.  
        !! First make sure pointers are not associated with any 
        !! specific memory location.  If they are, deallocate the memory, 
        !! then allocate the the proper amount of memory for each data field,
        !! error checking for each memory allocation to make
        !! sure memory is allocated successfully. 
        IF( ASSOCIATED( this%NumberSmallPixelColumns   ) ) &
            DEALLOCATE( this%NumberSmallPixelColumns   )
        IF( ASSOCIATED( this%SmallPixelColumn          ) ) &
            DEALLOCATE( this%SmallPixelColumn          )
        IF( ASSOCIATED( this%MeasurementQualityFlags   ) ) &
            DEALLOCATE( this%MeasurementQualityFlags   )
        IF( ASSOCIATED( this%InstrumentConfigurationId ) ) &
            DEALLOCATE( this%InstrumentConfigurationId )
        IF( ASSOCIATED( this%PixelQualityFlags         ) ) &
            DEALLOCATE( this%PixelQualityFlags         )
        IF( ASSOCIATED( this%sxData                    ) ) &
            DEALLOCATE( this%sxData                    )
        IF( ASSOCIATED( this%sxWl                      ) ) &
            DEALLOCATE( this%sxWl                      )
        IF( ASSOCIATED( this%AccumNsmpxCols            ) ) &
            DEALLOCATE( this%AccumNsmpxCols            )

        ALLOCATE( this%NumberSmallPixelColumns(this%nTimes),        &
                  this%SmallPixelColumn(this%nTimes),               &      
                  this%InstrumentConfigurationId(this%nTimes),      &      
                  this%MeasurementQualityFlags(this%nTimes),        &
                  this%AccumNsmpxCols(this%nTimes+1),               &
                  this%PixelQualityFlags(this%nXtrack,this%nTimes), &
                  this%sxData(this%nXtrack,this%nTimesSmallPixel),  &
                  this%sxWl(this%nXtrack,this%nTimesSmallPixel),    & 
                  STAT = ierr )
        IF( ierr .NE. zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                 "this:small pixel fields data block", &
                                 "L1Bsmpx_open", two )
           RETURN
        ENDIF 

        this%initialized = .TRUE.
        status = fill_smlpix_blk( this ) 
        RETURN      
      END FUNCTION L1Bsmpx_open 

!! 2. L1Bsmpx_getSWdims
 !    This function retrieves the dimension sizes from the L1B swath
 !    this: the block data structure
 !    nTimes: total number of lines in the swath
 !    nXtrack: number of cross-track pixels in the swath
 !    nWavel: number of wavelengths for each line and pixel number
 !    nWavelCoef: number of wavelength coefficients
!!
      FUNCTION L1Bsmpx_getSWdims( this, nTimes_k, nXtrack_k, &
                                  nTimesSmallPixel_k ) RESULT( status )
        TYPE (L1B_smlpix_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nTimes_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nTimesSmallPixel_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nXtrack_k
        INTEGER (KIND = 4) :: status

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bsmpx_getSWdims", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( PRESENT( nTimes_k  ) ) nTimes_k  = this%nTimes
        IF( PRESENT( nXtrack_k ) ) nXtrack_k = this%nXtrack
        IF( PRESENT( nTimesSmallPixel_k ) ) &
            nTimesSmallPixel_k = this%nTimesSmallPixel
        status = OZT_S_SUCCESS
        RETURN
      END FUNCTION L1Bsmpx_getSWdims

!! Private function: fill_smlpix_blk
 !    This is a private function which is only used by other functions 
 !    in this MODULE to fill the block data structure. 
 !    this: the block data structure
 !    iLine: the starting line for the data block
 !    status: the return PGS_SMF status value
!!
      FUNCTION fill_smlpix_blk( this ) RESULT( status ) 
        USE HE4_class
        TYPE (L1B_smlpix_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 2), DIMENSION(this%nWavel,this%nXtrack,this%nTimes) :: &
                                                                           Lqau
        INTEGER (KIND = 4) :: status
        INTEGER (KIND = 4) :: getl1bdblk
        EXTERNAL              getl1bdblk
        INTEGER (KIND = 4) :: rank
        INTEGER (KIND = 4) :: iLine, nL, di
        INTEGER (KIND = 4), DIMENSION(1:3) :: dims
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "fill_smlpix_blk", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! read small pixel radiance

        iLine = 0
        nL    = this%nTimes
        status = getl1bdblk( this%filename, this%swathname, &
                             "PixelQualityFlags", DFNT_UINT16, &
                             iLine, nL, rank, dims, &
                             Lqau )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve PixelQualityFlags", &
                                 "fill_smlpix_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "MeasurementQualityFlags", DFNT_UINT16,& 
                             iLine, nL, rank, dims, &
                             this%MeasurementQualityFlags )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve MeasurementQualityFlags",&
                                 "fill_smlpix_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "InstrumentConfigurationId", DFNT_UINT8,& 
                             iLine, nL, rank, dims, &
                             this%InstrumentConfigurationId )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve InstrumentConfigurationId",&
                                 "fill_smlpix_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "NumberSmallPixelColumns", DFNT_INT8, & 
                             iLine, nL, rank, dims, &
                             this%NumberSmallPixelColumns )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve NumberSmallPixelColumns",&
                                 "fill_smlpix_blk", one )
           RETURN
        ENDIF

        this%AccumNsmpxCols(1) = 0
        DO di = 1, this%nTimes
           this%AccumNsmpxCols(di+1) = this%AccumNsmpxCols(di) &
                                     + this%NumberSmallPixelColumns(di) 
        ENDDO

        !! Check for consistency
        IF( this%AccumNsmpxCols(this%nTimes+1) /= this%nTimesSmallPixel ) THEN
           status = OZT_E_INPUT
           ierr = OMI_SMF_setmsg( status, "nTimesSmallPixel not consistent "//&
                                 "with the sum of NumberSmallPixelColumns",  &
                                 "fill_smlpix_blk", zero )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "SmallPixelColumn", DFNT_INT16,& 
                             iLine, nL, rank, dims, &
                             this%SmallPixelColumn )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve SmallPixelColumn",&
                                 "fill_smlpix_blk", one )
           RETURN
        ENDIF

        iLine = 0
        nL    = this%nTimesSmallPixel
        status = getl1bdblk( this%filename, this%swathname, &
                            "SmallPixelWavelength", DFNT_FLOAT32, & 
                             iLine, nL, rank, dims, &
                             this%sxWl )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve SmallPixelWavelength", &
                                 "fill_smlpix_blk", one )
           RETURN
        ENDIF
        IF( INDEX( this%swathname, "Earth" ) > 0 ) THEN
           status = getl1bdblk( this%filename, this%swathname, &
                               "SmallPixelRadiance", DFNT_FLOAT32, & 
                                iLine, nL, rank, dims, &
                                this%sxData )
           IF( status .NE. OZT_S_SUCCESS ) THEN
              ierr = OMI_SMF_setmsg( status, "retrieve SmallPixelRadiance", &
                                    "fill_smlpix_blk", one )
              RETURN
           ENDIF
        ELSE IF( INDEX( this%swathname, "Sun" ) > 0 ) THEN
           status = getl1bdblk( this%filename, this%swathname, &
                               "SmallPixelIrradiance", DFNT_FLOAT32, & 
                                iLine, nL, rank, dims, &
                                this%sxData )
           IF( status .NE. OZT_S_SUCCESS ) THEN
              ierr = OMI_SMF_setmsg( status, "retrieve SmallPixelIrradiance", &
                                    "fill_smlpix_blk", one )
              RETURN
           ENDIF
        ELSE
           WRITE( msg,'(A)' ) "Unknow swathname:" // TRIM( this%swathname )
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "fill_smlpix_blk", one )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! this%SmallPixelColumn(di) is in C convetions
        !! add 1 to become fortran index
        DO iLine = 1, this%nTimes
          di = this%SmallPixelColumn(iLine) + 1
            this%PixelQualityFlags(1:this%nXtrack,iLine) &
          = Lqau(di,1:this%nXtrack,iLine)
        ENDDO

        RETURN
      END FUNCTION fill_smlpix_blk

!! 3. L1Bsmpx_getLine
 !    This function gets one or more radiance data field values from
 !    the data block.
 !    this: the block data structure
 !    iLineSmPx: the line number in the L1B swath.  NOTE: this input is 0 based
 !           range from 0 to (nTimesSmallPixel-1) inclusive.
 !
 !    SmPxRadIrr_k, SmPxPixelQualityFlags_k, and SmPxWavelength_k
 !    are keyword arguments.  Only those present in the argument 
 !    list will be set by the the function.
 !
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bsmpx_getLine( this, iLineSmPx, SmPxRadIrr_k, &
                              SmPxWavelength_k, SmPxQualityFlags_k ) &
                              RESULT (status )
        TYPE (L1B_smlpix_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLineSmPx 
        REAL (KIND = 4), OPTIONAL, DIMENSION(:), INTENT( OUT ) :: & 
                                         SmPxRadIrr_k, SmPxWavelength_k
        INTEGER (KIND = 2), OPTIONAL, DIMENSION(:), INTENT( OUT ) :: &
                                              SmPxQualityFlags_k
        INTEGER :: i, j, jj
        INTEGER (KIND = 4) :: status

        status = OZT_S_SUCCESS
        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bsmpx_getLine", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLineSmPx < 0 .OR. iLineSmPx >= this%nTimesSmallPixel ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLineSmPx out of range", &
                                 "L1Bsmpx_getLine", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        j = iLineSmPx + 1

        IF( PRESENT( SmPxRadIrr_k ) ) THEN
           IF( SIZE( SmPxRadIrr_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input RadIrr_k array too small", &
                                    "L1Bsmpx_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           SmPxRadIrr_k(1:this%nXtrack) = this%sxData(1:this%nXtrack,j )
        ENDIF

        IF( PRESENT( SmPxQualityFlags_k ) ) THEN
           IF( SIZE( SmPxQualityFlags_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                "input SmPxQualityFlags_k array too small",&
                                "L1Bsmpx_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
            !! jj = iLineLgPx + 1 (add 1 for fortran scheme)
            jj = L1Bsmpx_iLineSmPx2LgPx( this, iLineSmPx ) + 1
            SmPxQualityFlags_k( 1:this%nXtrack )  &
          = this%PixelQualityFlags( 1:this%nXtrack,jj )
        ENDIF

        IF( PRESENT( SmPxWavelength_k ) ) THEN
           IF( SIZE( SmPxWavelength_k ) < this%nXtrack ) THEN 
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input SmPxWavelength_k array too small",&
                                    "L1Bsmpx_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           SmPxWavelength_k(1:this%nXtrack) = this%sxWl( 1:this%nXtrack,j )
        ENDIF

        RETURN
      END FUNCTION L1Bsmpx_getLine

!! 4.  L1Bsmpx_close
 !    This function should be called when the data block is no longer
 !    needed.  It deallocates all the allocated memory, and sets
 !    all the parameters to invalid values.
 !    this: the block data structure
 !    
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bsmpx_close( this ) RESULT( status )
        TYPE (L1B_smlpix_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4) :: status
        status = OZT_S_SUCCESS

        !! DEALLOCATE all the memory
        IF( ASSOCIATED( this%NumberSmallPixelColumns   ) ) &
            DEALLOCATE( this%NumberSmallPixelColumns   )
        IF( ASSOCIATED( this%SmallPixelColumn          ) ) &
            DEALLOCATE( this%SmallPixelColumn          )
        IF( ASSOCIATED( this%InstrumentConfigurationId ) ) &
            DEALLOCATE( this%InstrumentConfigurationId )
        IF( ASSOCIATED( this%MeasurementQualityFlags   ) ) &
            DEALLOCATE( this%MeasurementQualityFlags   )
        IF( ASSOCIATED( this%PixelQualityFlags         ) ) &
            DEALLOCATE( this%PixelQualityFlags         )
        IF( ASSOCIATED( this%sxData                    ) ) &
            DEALLOCATE( this%sxData                    )
        IF( ASSOCIATED( this%sxWl                      ) ) &
            DEALLOCATE( this%sxWl                      )
        IF( ASSOCIATED( this%AccumNsmpxCols            ) ) &
            DEALLOCATE( this%AccumNsmpxCols            )

        this%nTimes           = -1
        this%nTimesSmallPixel = -1
        this%nXtrack          = -1
        this%nWavel           = -1     
        this%filename         = ""
        this%swathname        = ""
        this%initialized      = .FALSE.
        RETURN
      END FUNCTION L1Bsmpx_close

!! 5. L1Bsmpx_getMeasurementQF
      FUNCTION L1Bsmpx_getMeasurementQF( this, iLineLgPx, mQF ) &
                                         RESULT (status)
        TYPE (L1B_smlpix_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLineLgPx
        INTEGER (KIND=2), INTENT( OUT ) :: mQF
        INTEGER :: jj
        INTEGER (KIND=4) :: status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bsmpx_getMeasurementQF", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLineLgPx < 0 .OR. iLineLgPx >= this%nTimes ) THEN
           WRITE( msg, '(A,I9,A)') "iLineLgPx = ", iLineLgPx, &
                 " out of range [0,", this%nTimes, ")"
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                 "L1Bsmpx_getMeasurementQF", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! add 1 for fortran array indexing scheme
        jj = iLineLgPx + 1 
        mQF = this%MeasurementQualityFlags(jj)
        status = OZT_S_SUCCESS

      END FUNCTION L1Bsmpx_getMeasurementQF

!! 6. L1Bsmpx_getInstConfigId
      FUNCTION L1Bsmpx_getInstConfigId( this, iLineLgPx, &
                                        instId, rebinningFlg ) RESULT (status)
        TYPE (L1B_smlpix_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLineLgPx
        INTEGER (KIND = 1), INTENT( OUT ) :: instId, rebinningFlg
        INTEGER (KIND=2) :: mQF
        INTEGER :: jj
        INTEGER (KIND=4) :: status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bsmpx_getInstConfigId", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLineLgPx < 0 .OR. iLineLgPx >= this%nTimes ) THEN
           WRITE( msg, '(A,I9,A)') "iLineLgPx = ", iLineLgPx, &
                 " out of range [0,", this%nTimes, ")"
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                  "L1Bsmpx_getInstConfigId", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF
        jj     = iLineLgPx  + 1
        instId = this%InstrumentConfigurationId(jj)
        mQF    = this%MeasurementQualityFlags(jj)
        rebinningFlg = IBITS( mQF, 7, 1 ) 
        status = OZT_S_SUCCESS

      END FUNCTION L1Bsmpx_getInstConfigId

!! 7. L1Bsmpx_iLineSmPx2LgPx
      FUNCTION L1Bsmpx_iLineSmPx2LgPx( this, iLineSmPx ) RESULT (iLineLgPx ) 
        TYPE (L1B_smlpix_type), INTENT( IN ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLineSmPx
        INTEGER (KIND = 4) :: iLineLgPx
        INTEGER :: j, jj

        !! both iLineSmPx, iLineLgPx starts with 0
        !! From small pixel array index j [1,nTimesSmallPixel] to
        !! global pixel jj [1, nTimes]
        j      = iLineSmPx + 1
        jj     = locate( this%AccumNsmpxCols, j )
        IF( j == this%AccumNsmpxCols(jj) ) jj = jj-1
        iLineLgPx = jj - 1
        RETURN
      END FUNCTION L1Bsmpx_iLineSmPx2LgPx

!! 8. L1Bsmpx_getAccuNcols
      FUNCTION L1Bsmpx_getAccuNcols( this, AccumNsmpxCols_k, & 
                                     NumberSmallPixelColumns_k, &
                                     SmallPixelColumn_k  ) &
                                     RESULT (status)
        TYPE (L1B_smlpix_type), INTENT( IN ) :: this
        INTEGER (KIND = 1), DIMENSION(:), INTENT(OUT), OPTIONAL ::  &
                                                  NumberSmallPixelColumns_k
        INTEGER (KIND = 2), DIMENSION(:), INTENT(OUT), OPTIONAL ::  &
                                                  SmallPixelColumn_k
        INTEGER (KIND = 4), DIMENSION(:), INTENT(OUT), OPTIONAL ::  &
                                                  AccumNsmpxCols_k
        INTEGER (KIND=4) :: status

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bsmpx_getAccuNcols", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( PRESENT( NumberSmallPixelColumns_k ) ) THEN
           IF( SIZE( NumberSmallPixelColumns_k ) < this%nTimes ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                            "input NumberSmallPixelColumns_k array too small", &
                            "L1Bsmpx_getAccuNcols", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
             NumberSmallPixelColumns_k(1:this%nTimes) &
           = this%NumberSmallPixelColumns(1:this%nTimes )
        ENDIF
           
        IF( PRESENT( SmallPixelColumn_k ) ) THEN
           IF( SIZE( SmallPixelColumn_k ) < this%nTimes ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                            "input SmallPixelColumn_k array too small", &
                            "L1Bsmpx_getAccuNcols", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
             SmallPixelColumn_k(1:this%nTimes) &
           = this%SmallPixelColumn(1:this%nTimes )
        ENDIF

        IF( PRESENT( AccumNsmpxCols_k ) ) THEN
           IF( SIZE( AccumNsmpxCols_k ) < this%nTimes+1 ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                   "input AccumNsmpxCols_k array too small", &
                                   "L1Bsmpx_getAccuNcols", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
             AccumNsmpxCols_k(1:this%nTimes+1) &
           = this%AccumNsmpxCols(1:this%nTimes+1 )
        ENDIF

        status = OZT_S_SUCCESS
        RETURN 
      END FUNCTION L1Bsmpx_getAccuNcols

END MODULE L1B_smlpix_class
