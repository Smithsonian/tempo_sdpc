!! F90
 ! 
 !!Description:
 ! MODULE L1B_radirr_class contains the OMIL1B radiance or irradiance data 
 ! structure and the functions for the creation, deletion and access of 
 ! specific information stored in the data structure. A user can use the 
 ! following functions:  L1Bri_open, L1Bri_close, 
 ! L1Bri_getSWdims, L1Bri_getLine, L1Bri_getLineWL.
 !
 !!Revision History:
 ! Revision 0.1  08/26/2003  Kai Yang (GEST/UMBC)
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

MODULE L1B_radirr_class
   USE PGS_PC_class
   USE OMI_SMF_class
   USE UTIL_tools_class
   USE HE4_class
   IMPLICIT NONE
   PUBLIC  :: L1Bri_open, L1Bri_close
   PUBLIC  :: L1Bri_getSWdims
   PUBLIC  :: L1Bri_getLine
   PUBLIC  :: L1Bri_getLineWL
   PUBLIC  :: L1Bri_getMeasurementQF
   PUBLIC  :: L1Bri_getInstConfigId
   PUBLIC  :: L1B_extractMqa
   PUBLIC  :: L1B_getIRRnmqa
   PUBLIC  :: L1B_selectIRR
!   PUBLIC  :: L1Bri_getLineWC

   PRIVATE :: fill_rad_blk
   PRIVATE :: fill_irr_blk
   INTEGER (KIND = 4), PARAMETER, PUBLIC :: MAX_NAME_LENGTH = 255
   INTEGER, PRIVATE :: ierr  !!error code returned from a function
   INTEGER, PARAMETER, PRIVATE :: zero  = 0
   INTEGER, PARAMETER, PRIVATE :: one   = 1
   INTEGER, PARAMETER, PRIVATE :: two   = 2
   INTEGER, PARAMETER, PRIVATE :: three = 3
   REAL(KIND=4), PRIVATE :: EPSILON10
!!
 ! L1B_radirr_type is designed to store a block of L1B swath data.
 ! A block contains all the (ir)radiance and wavelength
 ! data fields for a number of "scan" lines or "exposure times"
 ! in a L1B swath. 
 ! filename: the file name for the L1B file (including path)
 ! swathname: the swath name in the L1B file, either "Earth UV-1 Swath",
 !            "Earth UV-2 Swath", or "Earth VIS Swath"
 ! iLine: the starting line for the block.  It is a reference to the
 !        L1B file, where 0 is the first line, and (nTimes - 1)
 !        is the last line.
 ! eLine: the ending line for the block (inclusive).
 ! nLine: the number of lines containing in the block with iLine
 !        as start and eLine as the end (inclusive)
 ! nTimes: the total number of lines contained in the L1B swath
 ! nXtrack: the number of cross-track pixels in the swath
 ! nWavel: the number of wavelengths in the swath
 ! nWavelCoef: the number of wavelength coefficients 
 ! initialized: the logical variable to indicate whether the
 !              block has been initialized.  The block data structure
 !              has to be created (initialized) before it can be used.
 !              
 ! Radiance and Irradiance fields:
 !   Lmsa: RadianceMantissa          || IrradianceMantissa
 !   Lpsa: RadiancePrecisionMantissa || IrradiancePrecisionMantissa
 !   Lexp: RadianceExponent          || IrradianceExponent
 !   Lqau: PixelQualityFlags
 ! These are three-dimensional data arrays created for the number of
 ! lines, the number of pixels and the number of wavelengths for the
 ! data block when the data structure is created. 
 !
 ! Wavelength Coefficient fields:
 !   wCof: WavelengthCoefficient
 !   wPcf: WavelengthCoefficientPrecision
 ! These are threee-dimensional data arrays created for the number of
 ! lines, the number of pixels, and the number of wavelength coefficients 
 ! for the data block when the data structure is created. 
 !
 ! wRefCol : WavelengthReferenceColumn
 ! MeasurementQF: MeasurementQualityFlags
 ! instId: InstrumentConfigurationId
 !
 ! These are one-dimensional data arrays created for the number of
 ! lines for the data block when the data structure is created. 
 !
 ! Again, a user cannot see these data fields directly but the 
 ! information can be retrieved only through functions L1Bri_getLine, 
 ! L1Bri_getLineWL.  The "line" functions retrieve information for one 
 ! sepcific line defined by iLine for all the wavelengths, or a subset 
 ! of those wavelengths, or a set of wavelength values defined by the user.
!!

   TYPE, PUBLIC :: L1B_radirr_type
      PRIVATE
      CHARACTER ( LEN = MAX_NAME_LENGTH ) :: filename, swathname
      INTEGER (KIND = 4) :: iLine, eLine, nLine
      INTEGER (KIND = 4) :: nTimes
      INTEGER (KIND = 4) :: nTimesSmallPixel
      INTEGER (KIND = 4) :: nXtrack
      INTEGER (KIND = 4) :: nWavel
      INTEGER (KIND = 4) :: nWavelCoef
      LOGICAL :: initialized

      ! Radiance and Irradiance fields
      INTEGER (KIND = 2), DIMENSION(:,:,:), POINTER :: Lmsa, Lpsa, Lqau  
      INTEGER (KIND = 1), DIMENSION(:,:,:), POINTER :: Lexp
      REAL (KIND = 4), DIMENSION(:,:,:), POINTER :: wCof, wPcf
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: wRefCol
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: MeasurementQF
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: instId
   END TYPE L1B_radirr_type

   CONTAINS

!! 1. L1Bri_open
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
      FUNCTION L1Bri_open( this, fn, swn, nL ) RESULT (status)
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        CHARACTER ( LEN = * ), INTENT( IN ) :: fn, swn
        INTEGER (KIND = 4), OPTIONAL, INTENT( IN ) :: nL
        INTEGER (KIND = 4) :: swfid, swid, status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg
        
        !! nullify pointers so we don't run into allocation issues later
        nullify(this%Lmsa)
        nullify(this%Lpsa)
        nullify(this%Lqau)
        nullify(this%Lexp)
        nullify(this%wCof)
        nullify(this%wPcf)
        nullify(this%wRefCol)
        nullify(this%MeasurementQF)
        nullify(this%instId)
        

        !! open the L1B swath file
        swfid = swopen( fn, DFACC_READ )
        IF( swfid < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_FILE_OPEN, fn, "L1Bri_open", zero )
           RETURN
        ENDIF 

        !! attach to the swath
        swid = swattach( swfid, swn )
        IF( swid < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_SWATH_ATTACH, swn, &
                                 "L1Bri_open", zero )
           RETURN
        ENDIF 
        
        !! intitialize the start and end line to invalid values
        this%iLine      = -1
        this%eLine      = -1

        !! retrieve the dimension info from the swath file.
        !! dimension names are obtained from file spec
        status = swrdattr( swid, "NumTimes", this%nTimes )
        IF( status < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get nTimes size failed", &
                                 "L1Bri_open", zero )
           RETURN
        ENDIF 
        IF( this%nTimes <= zero ) THEN
           WRITE( msg,'(A,I3,A)' ) "nTimes =", this%nTimes, &
                 ", No scan lines in "//TRIM(swn)// ", " // TRIM(fn) 
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, "L1Bri_open", zero )
           RETURN
        ENDIF 

        status = swrdattr( swid, "NumTimesSmallPixel", this%nTimesSmallPixel )
        IF( status < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, &
                                 "get nTimesSmallPixel size failed", &
                                 "L1Bri_open", three )
           this%nTimesSmallPixel = -1
        ENDIF 

        this%nXtrack    = swdiminfo( swid, "nXtrack" )
        IF( this%nXtrack < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get xTrack size failed", &
                                  "L1Bri_open", zero )
           RETURN
        ENDIF 

        this%nWavel     = swdiminfo( swid, "nWavel" )
        IF( this%nWavel < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get nWavel size failed", &
                                  "L1Bri_open", zero )
           RETURN
        ENDIF 

        this%nWavelCoef = swdiminfo( swid, "nWavelCoef" )
        IF( this%nWavelCoef < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get nWavelCoef size failed", &
                                  "L1Bri_open", zero )
           RETURN
        ENDIF 

        !! detach and close the L1B swath files.  No need for 
        !! error checking, for an error is unlikely to occur here.
        ierr = swdetach( swid )
        ierr = swclose( swfid )
          
        !! set the block size according to input, if present
        !! if not set it to 100
        IF( .not. PRESENT( nL) ) THEN
           this%nLine = 100
        ELSE
           IF( nL <= 0 ) THEN
              this%nLine = 100
           ELSE 
              this%nLine = nL
           ENDIF
        ENDIF

        !! the maximum number of lines cannot exceed the total
        !! number of lines in the L1B swath.
        IF( this%nLine > this%nTimes ) this%nLine = this%nTimes
        
        this%filename = fn   
        this%swathname = swn
        
        !! Allocate the memory for storage of radiance fields.  
        !! First make sure pointers are not associated with any 
        !! specific memory location.  If they are, deallocate the memory, 
        !! then allocate the the proper amount of memory for each data field,
        !! error checking for each memory allocation to make
        !! sure memory is allocated successfully. 

        IF( ASSOCIATED( this%Lmsa ) ) deallocate( this%Lmsa )
        IF( ASSOCIATED( this%Lpsa ) ) deallocate( this%Lpsa )
        IF( ASSOCIATED( this%Lqau ) ) deallocate( this%Lqau )
        IF( ASSOCIATED( this%Lexp ) ) deallocate( this%Lexp )

        ALLOCATE( this%Lmsa(this%nWavel,this%nXtrack,this%nLine), &
                  this%Lpsa(this%nWavel,this%nXtrack,this%nLine), &
                  this%Lqau(this%nWavel,this%nXtrack,this%nLine), &
                  this%Lexp(this%nWavel,this%nXtrack,this%nLine), STAT = ierr )
        IF( ierr .NE. zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                 "this:radiance fields data block", &
                                 "L1Bri_open", two )
           RETURN
        ENDIF 

        IF( ASSOCIATED( this%wCof          ) ) deallocate( this%wCof          )
        IF( ASSOCIATED( this%wPcf          ) ) deallocate( this%wPcf          )
        IF( ASSOCIATED( this%wRefCol       ) ) deallocate( this%wRefCol       )
        IF( ASSOCIATED( this%MeasurementQF ) ) deallocate( this%MeasurementQF )
        IF( ASSOCIATED( this%instId        ) ) deallocate( this%instId        )

        ALLOCATE( this%wCof(this%nWavelCoef,this%nXtrack,this%nLine), &
                  this%wPcf(this%nWavelCoef,this%nXtrack,this%nLine), &
                  this%wRefCol(this%nLine), STAT = ierr )
        IF( ierr .NE. zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                 "this:wavelength coefficient data block", &
                                 "L1Bri_open", two )
           RETURN
        ENDIF 

        ALLOCATE( this%MeasurementQF(this%nLine), &
                  this%instId(this%nLine), STAT = ierr )
        IF( ierr .NE. zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                 "this:MeasurementQF allocation failure", &
                                 "L1Bri_open", two )
           RETURN
        ENDIF 

        this%initialized = .TRUE.
        status = OZT_S_SUCCESS
        RETURN      
      END FUNCTION L1Bri_open

      subroutine define_global_epsilon10
        EPSILON10 = 10*EPSILON(1.0)
      end subroutine

!! 2. L1Bri_getSWdims

 !    This function retrieves the dimension sizes from the L1B swath
 !    this: the block data structure
 !    nTimes: total number of lines in the swath
 !    nXtrack: number of cross-track pixels in the swath
 !    nWavel: number of wavelengths for each line and pixel number
 !    nWavelCoef: number of wavelength coefficients
!!
      FUNCTION L1Bri_getSWdims( this, nTimes_k, nXtrack_k, &
                               nWavel_k, nWavelCoef_k, &
                               nTimesSmallPixel_k ) RESULT( status )
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nTimes_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nXtrack_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nWavel_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nWavelCoef_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nTimesSmallPixel_k
        INTEGER (KIND = 4) :: status

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bri_getSWdims", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( PRESENT( nTimes_k     ) ) nTimes_k     = this%nTimes
        IF( PRESENT( nXtrack_k    ) ) nXtrack_k    = this%nXtrack
        IF( PRESENT( nWavel_k     ) ) nWavel_k     = this%nWavel
        IF( PRESENT( nWavelCoef_k ) ) nWavelCoef_k = this%nWavelCoef
        IF( PRESENT( nTimesSmallPixel_k ) ) &
                                      nTimesSmallPixel_k = this%nTimesSmallPixel
        status = OZT_S_SUCCESS
        RETURN
      END FUNCTION L1Bri_getSWdims

!! Private function: fill_rad_blk
 !    This is a private function which is only used by other functions 
 !    in this MODULE to fill the block data structure. 
 !    this: the block data structure
 !    iLine: the starting line for the data block
 !    status: the return PGS_SMF status value
!!
      FUNCTION fill_rad_blk( this, iLine ) RESULT( status ) 
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLine
        INTEGER (KIND = 4) :: status
        INTEGER (KIND = 4) :: getl1bdblk
        EXTERNAL              getl1bdblk
        INTEGER (KIND = 4) :: rank
        INTEGER (KIND = 4) :: nL
        INTEGER (KIND = 4), DIMENSION(1:3) :: dims

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "fill_rad_blk", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        this%iLine = iLine
        
        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "iLine out of range", "fill_rad_blk", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! If, from the starting line, the number of lines to be read into
        !! the block data structure goes past the final line in the L1B swath,
        !! then recalculate the number of lines to be read. 

        IF( (iLine + this%nLine) > this%nTimes ) THEN          
           nL = this%nTimes - iLine 
        ELSE
           nL = this%nLine
        ENDIF
        this%eLine = this%iLine + nL - 1

        !! read radiance
        status = getl1bdblk( this%filename, this%swathname, &
                            "RadianceMantissa", DFNT_INT16, & 
                             iLine, nL, rank, dims, &
                             this%Lmsa )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve RadianceMantissa", &
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "RadiancePrecisionMantissa", DFNT_INT16, & 
                             iLine, nL, rank, dims, &
                             this%Lpsa )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve RadiancePrecisionMantissa", &
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "RadianceExponent", DFNT_INT8, & 
                             iLine, nL, rank, dims, &
                             this%Lexp )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve RadianceExponent", &
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "PixelQualityFlags", DFNT_UINT16, & 
                             iLine, nL, rank, dims, &
                             this%Lqau )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve PixelQualityFlags", &
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        !! read wavelength coefficient
        status = getl1bdblk( this%filename, this%swathname, &
                            "WavelengthCoefficient", DFNT_FLOAT32, & 
                             iLine, nL, rank, dims, &
                             this%wCof )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve WavelengthCoefficient", &
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "WavelengthCoefficientPrecision", DFNT_FLOAT32,& 
                             iLine, nL, rank, dims, &
                             this%wPcf )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve WavelengthCoefficientPrecision",&
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "WavelengthReferenceColumn", DFNT_INT16,& 
                             iLine, nL, rank, dims, &
                             this%wRefCol )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve WavelengthReferenceColumn",&
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "MeasurementQualityFlags", DFNT_UINT16,& 
                             iLine, nL, rank, dims, &
                             this%MeasurementQF )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve MeasurementQualityFlags",&
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "InstrumentConfigurationId", DFNT_UINT8,& 
                             iLine, nL, rank, dims, &
                             this%instId )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve InstrumentConfigurationId",&
                                 "fill_rad_blk", one )
           RETURN
        ENDIF

        RETURN
      END FUNCTION fill_rad_blk

!! Private function: fill_irr_blk
 !    This is a private function which is only used by other functions 
 !    in this MODULE to fill the block data structure. 
 !    this: the block data structure
 !    status: the return PGS_SMF status value
!!
      FUNCTION fill_irr_blk( this, iLine ) RESULT( status ) 
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLine
        INTEGER (KIND = 4) :: status
        INTEGER (KIND = 4) :: getl1bdblk
        EXTERNAL              getl1bdblk
        INTEGER (KIND = 4) :: rank
        INTEGER (KIND = 4) :: nL
        INTEGER (KIND = 4), DIMENSION(1:2) :: dims

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "fill_irr_blk", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        this%iLine = iLine

        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "iLine out of range", "fill_irr_blk", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! If, from the starting line, the number of lines to be read into
        !! the block data structure goes past the final line in the L1B swath,
        !! then recalculate the number of lines to be read.

        IF( (iLine + this%nLine) > this%nTimes ) THEN
           nL = this%nTimes - iLine
        ELSE
           nL = this%nLine
        ENDIF
        this%eLine = this%iLine + nL - 1

        !! read irradiance
        status = getl1bdblk( this%filename, this%swathname, &
                            "IrradianceMantissa", DFNT_INT16, & 
                             iLine, nL, rank, dims, this%Lmsa )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve IrradianceMantissa", &
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "IrradiancePrecisionMantissa", DFNT_INT16, & 
                             iLine, nL, rank, dims, this%Lpsa )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve IrradiancePrecisionMantissa", &
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "IrradianceExponent", DFNT_INT8, & 
                             iLine, nL, rank, dims, this%Lexp )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve IrradianceExponent", &
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "PixelQualityFlags", DFNT_UINT16, & 
                             iLine, nL, rank, dims, this%Lqau )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve PixelQualityFlags", &
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        !! read wavelength coefficient
        status = getl1bdblk( this%filename, this%swathname, &
                            "WavelengthCoefficient", DFNT_FLOAT32, & 
                             iLine, nL, rank, dims, this%wCof )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve WavelengthCoefficient", &
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "WavelengthCoefficientPrecision", DFNT_FLOAT32,& 
                             iLine, nL, rank, dims, this%wPcf )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve WavelengthCoefficientPrecision",&
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "WavelengthReferenceColumn", DFNT_INT16,& 
                             iLine, nL, rank, dims, this%wRefCol )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve WavelengthReferenceColumn",&
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "MeasurementQualityFlags", DFNT_UINT16,& 
                             iLine, nL, rank, dims, this%MeasurementQF )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve MeasurementQualityFlags",&
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                            "InstrumentConfigurationId", DFNT_UINT8,& 
                             iLine, nL, rank, dims, &
                             this%instId )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve InstrumentConfigurationId",&
                                 "fill_irr_blk", one )
           RETURN
        ENDIF

        RETURN
      END FUNCTION fill_irr_blk

!! 3. L1Bri_getLine
 !    This function gets one or more radiance data field values from
 !    the data block.
 !    this: the block data structure
 !    iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !           range from 0 to (nTimes-1) inclusive.
 !    Wlmin_k: input wavelength value for the band beginning (in nm)   
 !    Wlmax_k: input wavelength value for the band end (in nm)
 !    (These two inputs are optional.  Note that Wlmin_k < Wlmax_k.)
 !
 !    RadIrr_k, RadIrrPrecision_k, PixelQualityFlags_k, Wavelength_k
 !    and Nwl_k are keyword arguments.  Only those present in the argument 
 !    list will be set by the the function.
 !
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bri_getLine( this, iLine, Wlmin_k, Wlmax_k, &
                              RadIrr_k, RadIrrPrecision_k, &
                              PixelQualityFlags_k, Wavelength_k, &
                              Nwl_k) RESULT (status)
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLine 
        INTEGER (KIND = 4), OPTIONAL, INTENT( OUT ) :: Nwl_k
        REAL (KIND = 4), OPTIONAL, INTENT( IN ) :: Wlmin_k, &
                                                   Wlmax_k    ! wavelength range
        REAL (KIND = 4), OPTIONAL, DIMENSION(:,:), INTENT( OUT ) :: & 
                                         RadIrr_k, RadIrrPrecision_k, &
                                         Wavelength_k
        INTEGER (KIND = 2), OPTIONAL, DIMENSION(:,:), INTENT( OUT ) :: &
                                              PixelQualityFlags_k
        REAL (KIND = 4), DIMENSION(1:this%nWavel,1:this%nXtrack) :: &
                                              wl_local
        INTEGER :: i, j, k, di, ic
        INTEGER (KIND = 4) :: Nwl_l
        INTEGER (KIND = 4) :: il, ih, i_foo(2)
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        status = OZT_S_SUCCESS
        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bri_getLine", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLine out of range", &
                                 "L1Bri_getLine", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < this%iLine .OR. iLine > this%eLine ) THEN
           IF( INDEX( this%swathname, "Earth" ) > 0 ) THEN
              status = fill_rad_blk( this, iLine )
           ELSE IF( INDEX( this%swathname, "Sun" ) > 0 ) THEN
              status = fill_irr_blk( this, iLine )
           ELSE
              WRITE( msg,'(A)' ) "Unknow swathname:" // TRIM( this%swathname )
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           IF( status .NE. OZT_S_SUCCESS ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, "retrieve data block", &
                                    "L1Bri_getLine", zero )
              RETURN
           ENDIF
        ENDIF

        !! j is the index into the block
        j = iLine - this%iLine + 1

        
        IF( PRESENT(Wlmin_k) .OR. PRESENT(Wlmax_k) .OR. PRESENT(Wavelength_k) ) THEN
           DO i = 1, this%nXtrack
           DO k = 1, this%nWavel
             di = k-1-this%wRefCol(j)   ! wavelength index is 0 based in HDF file.
                                        ! So subtract 1 from FORTRAN index k.
             wl_local(k,i) = di*this%wCof(this%nWavelCoef,i,j)
             DO ic = this%nWavelCoef - 1, 2, -1
               wl_local(k,i) = di*(wl_local(k,i) + this%wCof(ic,i,j)) 
             ENDDO
             wl_local(k,i) = wl_local(k,i) + this%wCof(1,i,j)
           ENDDO
           ENDDO
        ENDIF

        IF( PRESENT( Wlmin_k ) .AND. PRESENT(  Wlmax_k ) ) THEN
           IF( (Wlmin_k-Wlmax_k ) > 0.0 ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input Wlmin_k greater than Wlmax_k", &
                                    "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        IF( PRESENT( Wlmin_k ) ) THEN
           IF( Wlmin_k > MAXVAL( wl_local ) )THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input Wlmin_k out of bound", &
                                    "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
 
           IF( Wlmin_k < MINVAL(wl_local) ) THEN
              il = 1
           ELSE
              i_foo = MAXLOC( wl_local, MASK = wl_local-Wlmin_k <= 0.0 )
              il    = i_foo( 1 )
           ENDIF            
        ELSE
           il = 1
        ENDIF

        IF( PRESENT( Wlmax_k ) ) THEN
           IF( Wlmax_k < MINVAL( wl_local ) ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input Wlmax_k out of bound", &
                                    "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           IF( Wlmax_k > MAXVAL( wl_local ) ) THEN
              ih = this%nWavel
           ELSE
              i_foo = MINLOC( wl_local, MASK = wl_local-Wlmax_k >= 0.0 )
              ih    = i_foo( 1 )
           ENDIF       
        ELSE
           ih = this%nWavel
        ENDIF

        Nwl_l = ih-il+1

        IF( PRESENT( RadIrr_k ) ) THEN
           IF( SIZE( RadIrr_k, 1 ) < Nwl_l  .OR. &
               SIZE( RadIrr_k, 2 ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input RadIrr_k array too small", &
                                    "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           RadIrr_k(1:Nwl_l, 1:this%nXtrack) =  &
                             this%Lmsa(il:ih,1:this%nXtrack,j ) &
                           * 10.0**this%Lexp(il:ih,1:this%nXtrack,j)
        ENDIF

        IF( PRESENT( RadIrrPrecision_k ) ) THEN
           IF( SIZE( RadIrrPrecision_k, 1 ) < Nwl_l .OR. &
               SIZE( RadIrrPrecision_k, 2 ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                   "input RadIrrPrecision_k array too small",&
                                    "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           RadIrrPrecision_k(1:Nwl_l, 1:this%nXtrack) = &
                             this%Lpsa(il:ih,1:this%nXtrack,j) &
                           * 10.0**this%Lexp(il:ih,1:this%nXtrack,j)
        ENDIF

        IF( PRESENT( PixelQualityFlags_k ) ) THEN
           IF( SIZE( PixelQualityFlags_k, 1 ) < Nwl_l .OR. &
               SIZE( PixelQualityFlags_k, 2 ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                "input PixelQualityFlags_k array too small",&
                                "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           PixelQualityFlags_k(1:Nwl_l, 1:this%nXtrack ) = &
                                         this%Lqau(il:ih,1:this%nXtrack,j)
        ENDIF

        IF( PRESENT( Wavelength_k ) ) THEN
           IF( SIZE( Wavelength_k, 1 ) < Nwl_l .OR. &
               SIZE( Wavelength_k, 2 ) < this%nXtrack ) THEN 
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input Wavelength_k array too small",&
                                    "L1Bri_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           Wavelength_k(1:Nwl_l, 1:this%nXtrack) = &
                                         wl_local(il:ih,1:this%nXtrack)
        ENDIF

        IF( PRESENT( Nwl_k ) ) Nwl_k = Nwl_l

        RETURN
      END FUNCTION L1Bri_getLine

!! 4. L1Bri_getLineWL
 !    This function gets one or more radiance data field values from
 !    the data block.
 !    this: the block data structure
 !    iLine : the line number in the L1B swath.  NOTE: this input is 0 based
 !            range from 0 to (nTimes-1) inclusive.
 !    Wavelength_k is list of input wavelengths.  Only information for
 !    these specific wavelengths is returned.
 !    Nwl_k is optional.  If present, then it indicates the number of valid
 !    wavelengths in Wavelength_k.  If not, then all the values in the
 !    Wavelength_k are assumed to be valid.
 !
 !    RadIrr_k, RadIrrPrecision_k, and PixelQualityFlags_k
 !    are keyword arguments.  Only those present in the argument 
 !    list will be set by the function.
 !
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bri_getLineWL( this, iLine, &
                                Wavelength_k, Nwl_k, RadIrr_k, &
                                RadIrrPrecision_k, &
                                PixelQualityFlags_k ) &
                                RESULT (status)

        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLine
        REAL (KIND = 4), DIMENSION(:), INTENT( IN ) :: Wavelength_k
        INTEGER (KIND = 4), OPTIONAL, INTENT( IN ) :: Nwl_k
        REAL (KIND = 4), OPTIONAL, DIMENSION(:,:), INTENT( OUT ) :: &
                                         RadIrr_k, RadIrrPrecision_k
        INTEGER (KIND=2), OPTIONAL, DIMENSION(:,:), INTENT( OUT ) :: &
                                                PixelQualityFlags_k
        REAL (KIND = 4), DIMENSION(1:this%nWavel,1:this%nXtrack) :: &
                                                wl_local
        INTEGER (KIND = 4) :: Nwl_l
        INTEGER :: i, j, k, di, ic
        INTEGER :: iw
        INTEGER (KIND = 4) :: il, ih
        REAL (KIND = 4) :: ri_il, ri_ih, frac
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        call define_global_epsilon10

        status = OZT_S_SUCCESS
        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bri_getLineWL", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLine out of range", & 
                                 "L1Bri_getLineWL", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < this%iLine .OR. iLine > this%eLine ) THEN
           IF( INDEX( this%swathname, "Earth" ) > 0 ) THEN
              status = fill_rad_blk( this, iLine )
           ELSE IF( INDEX( this%swathname, "Sun" ) > 0 ) THEN
              status = fill_irr_blk( this, iLine )
           ELSE
              WRITE( msg,'(A)' ) "Unknow swathname:" // TRIM( this%swathname )
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L1Bri_getLineWL", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           IF( status .NE. OZT_S_SUCCESS ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, "retrieve data block", &
                                    "L1Bri_getLineWL", zero )
              RETURN
           ENDIF
        ENDIF

        j = iLine - this%iLine + 1
        DO i = 1, this%nXtrack
        DO k = 1, this%nWavel
          di = k-1-this%wRefCol(j)   ! wavelength index is 0 based in HDF file.
                                     ! So subtract 1 from FORTRAN index k.
          wl_local(k,i) = di*this%wCof(this%nWavelCoef,i,j)
          DO ic = this%nWavelCoef - 1, 2, -1
            wl_local(k,i) = di*(wl_local(k,i) + this%wCof(ic,i,j)) 
          ENDDO
          wl_local(k,i) = wl_local(k,i) + this%wCof(1,i,j)
        ENDDO
        ENDDO

        IF( PRESENT( Nwl_k ) ) THEN
           Nwl_l = Nwl_k
           IF( SIZE( Wavelength_k ) < Nwl_l ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, "input Nwl_l too large",&
                                    "L1Bri_getLineWL", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
        ELSE
           Nwl_l = SIZE( Wavelength_k )
        ENDIF

        IF( PRESENT( RadIrr_k) ) THEN
           IF( SIZE( RadIrr_k, 1 ) < Nwl_l .OR. &
               SIZE( RadIrr_k, 2 ) < this%nXtrack ) THEN
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                     "input RadIrr_k array too small", &
                                     "L1Bri_getLineWL", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        IF( PRESENT( RadIrrPrecision_k ) ) THEN
           IF( SIZE( RadIrrPrecision_k, 1 ) < Nwl_l .OR. &
               SIZE( RadIrrPrecision_k, 2 ) <  this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "input RadIrrPrecision_k array too small",&
                                    "L1Bri_getLineWL", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        IF( PRESENT( PixelQualityFlags_k ) ) THEN
           IF( SIZE( PixelQualityFlags_k, 1 ) < Nwl_l .OR. &
               SIZE( PixelQualityFlags_k, 2 ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input PixelQualityFlags_k array too small",&
                                 "L1Bri_getLineWL", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        DO iw = 1, Nwl_l
           il = iw
        DO i = 1, this%nXtrack
           il = hunt( wl_local(1:this%nWavel,i ), Wavelength_k(iw), il )
           IF( il == 0 .OR. il == this%nWavel ) THEN
              !! no measurement data for this wl, return 0 values
              !! and set PixelQualityFlags to be missing (bit 1 equal 1)
              IF( PRESENT(RadIrr_k) ) RadIrr_k(iw,i) =  0
              IF( PRESENT(RadIrrPrecision_k) ) &
                                       RadIrrPrecision_k(iw,i) = 0
              IF( PRESENT(PixelQualityFlags_k) ) PixelQualityFlags_k(iw,i) = 1
           ELSE
              ih = il+1
              frac = (Wavelength_k(iw) - wl_local(il,i))/ &
                     (wl_local(ih,i)   - wl_local(il,i))
              IF( ABS( frac ) < EPSILON10 ) THEN
                  ih = il
                  frac = 0.0
              ELSE IF( ABS( frac -1.0 ) < EPSILON10 ) THEN
                  il = ih
                  frac = 0.0
              ENDIF

              IF( ih == il ) THEN
                 IF( PRESENT(RadIrr_k) ) RadIrr_k(iw,i) &  
                    = this%Lmsa(il,i,j) * 10.0**this%Lexp(il,i,j)
                 IF( PRESENT(RadIrrPrecision_k) ) RadIrrPrecision_k(iw,i) &
                    = this%Lpsa(il,i,j) * 10.0**this%Lexp(il,i,j)
                 IF( PRESENT(PixelQualityFlags_k) ) PixelQualityFlags_k(iw,i) &
                    = this%Lqau(il,i,j)
              ELSE 
                 !! no need for fill value check on Lmsa and Lexp,  
                 !! ri_il and ri_ih will not go underflow or overflow
                 !! for all the possible values of Lmsa and Lexp.
                 !! PixelQualityFlags_k is used to indicate if there is
                 !! any eorror or warning associated with the values
                 !! (interpolated or not). 

                 IF( PRESENT(RadIrr_k) ) THEN
                    ri_il = this%Lmsa(il,i,j) * 10.0**this%Lexp(il,i,j)
                    ri_ih = this%Lmsa(ih,i,j) * 10.0**this%Lexp(ih,i,j)
                    RadIrr_k(iw,i) = ri_il + frac*( ri_ih - ri_il )
                 ENDIF
   
                 IF( PRESENT(RadIrrPrecision_k) ) THEN
                    ri_il = this%Lpsa(il,i,j) * 10.0**this%Lexp(il,i,j)
                    ri_ih = this%Lpsa(ih,i,j) * 10.0**this%Lexp(ih,i,j)
                    RadIrrPrecision_k(iw,i) = MAX( ri_il, ri_ih )
                 ENDIF
   
                 !! combine both the setted flags in il and ih, and return
                 !! it as the qa flag for the interpolated wl.
                 IF( PRESENT(PixelQualityFlags_k) ) THEN
                    PixelQualityFlags_k(iw,i) = IOR( this%Lqau(il,i,j), &
                                                     this%Lqau(ih,i,j) )
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        ENDDO

        RETURN
      END FUNCTION L1Bri_getLineWL

!! 5.  L1Bri_close
 !    This function should be called when the data block is no longer
 !    needed.  It deallocates all the allocated memory, and sets
 !    all the parameters to invalid values.
 !    this: the block data structure
 !    
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bri_close( this ) RESULT( status )
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4) :: status
        status = OZT_S_SUCCESS

        !! DEALLOCATE all the memory
        IF( ASSOCIATED( this%Lmsa          ) ) DEALLOCATE( this%Lmsa          )
        IF( ASSOCIATED( this%Lpsa          ) ) DEALLOCATE( this%Lpsa          )
        IF( ASSOCIATED( this%Lqau          ) ) DEALLOCATE( this%Lqau          )
        IF( ASSOCIATED( this%Lexp          ) ) DEALLOCATE( this%Lexp          )
        IF( ASSOCIATED( this%wCof          ) ) DEALLOCATE( this%wCof          )
        IF( ASSOCIATED( this%wPcf          ) ) DEALLOCATE( this%wPcf          )
        IF( ASSOCIATED( this%wRefCol       ) ) DEALLOCATE( this%wRefCol       )
        IF( ASSOCIATED( this%MeasurementQF ) ) DEALLOCATE( this%MeasurementQF )
        IF( ASSOCIATED( this%instId        ) ) DEALLOCATE( this%instId        )

        this%iLine       = -1
        this%nLine       = -1
        this%eLine       = -1
        this%nTimes      = -1
        this%nXtrack     = -1
        this%nWavel      = -1     
        this%nWavelCoef  = -1     
        this%filename    = ""
        this%swathname   = ""
        this%initialized = .FALSE.
        RETURN
      END FUNCTION L1Bri_close

!! 6. L1Bri_getMeasurementQF
      FUNCTION L1Bri_getMeasurementQF( this, iLine, mQF ) RESULT (status)
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLine
        INTEGER (KIND=2), INTENT( OUT ) :: mQF
        INTEGER :: j
        INTEGER (KIND=4) :: status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bri_getMeasurementQF", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLine out of range", & 
                                 "L1Bri_getMeasurementQF", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < this%iLine .OR. iLine > this%eLine ) THEN
           IF( INDEX( this%swathname, "Earth" ) > 0 ) THEN
              status = fill_rad_blk( this, iLine )
           ELSE IF( INDEX( this%swathname, "Sun" ) > 0 ) THEN
              status = fill_irr_blk( this, iLine )
           ELSE
              WRITE( msg,'(A)' ) "Unknow swathname:" // TRIM( this%swathname )
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                     "L1Bri_getMeasurementQF", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           IF( status .NE. OZT_S_SUCCESS ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, "retrieve data block", &
                                    "L1Bri_getMeasurementQF", zero )
              RETURN
           ENDIF
        ENDIF

        j = iLine - this%iLine + 1
        mQF = this%MeasurementQF(j)
        status = OZT_S_SUCCESS

      END FUNCTION L1Bri_getMeasurementQF

!! 7. L1Bri_getInstConfigId
      FUNCTION L1Bri_getInstConfigId( this, iLine, instId, &
               rebinningFlg ) RESULT (status)
        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLine
        INTEGER (KIND = 1), INTENT( OUT ) :: instId
        INTEGER (KIND = 1), INTENT( OUT ), OPTIONAL :: rebinningFlg
        INTEGER (KIND=2) :: mQF
        INTEGER :: j
        INTEGER (KIND=4) :: status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bri_getInstConfigId", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLine out of range", & 
                                 "L1Bri_getInstConfigId", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < this%iLine .OR. iLine > this%eLine ) THEN
           IF( INDEX( this%swathname, "Earth" ) > 0 ) THEN
              status = fill_rad_blk( this, iLine )
           ELSE IF( INDEX( this%swathname, "Sun" ) > 0 ) THEN
              status = fill_irr_blk( this, iLine )
           ELSE
              WRITE( msg,'(A)' ) "Unknow swathname:" // TRIM( this%swathname )
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                     "L1Bri_getInstConfigId", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           IF( status .NE. OZT_S_SUCCESS ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, "retrieve data block", &
                                    "L1Bri_getInstConfigId", zero )
              RETURN
           ENDIF
        ENDIF

        j = iLine - this%iLine + 1
        instId = this%instId(j)
        mQF    = this%MeasurementQF(j)
        IF( PRESENT( rebinningFlg ) ) rebinningFlg = int(IBITS( mQF, 7, 1 ),kind(rebinningFlg))
        status = OZT_S_SUCCESS

      END FUNCTION L1Bri_getInstConfigId

!! 8. L1B_extractMqa
      FUNCTION L1B_extractMqa( measurementQAflags, mqaL_warning, &
                               rebinF, saaF, maneuverF ) &
                               RESULT( mqaL_error )
        INTEGER (KIND = 2), INTENT(IN) :: measurementQAflags
        INTEGER (KIND = 1), INTENT(OUT) :: mqaL_warning
        INTEGER (KIND = 1), INTENT(OUT), OPTIONAL :: rebinF, saaF, maneuverF 
        INTEGER (KIND = 1) :: mqaL_error

        !!bit type is set according to the recommendation described in
        !!RP-OMIE-KNMI-396 Interpretation GDPS flags Issue 1, 22 November 2002.
        !!E for error, data should not be used
        !!W for warning, data can be used, but warning flags set
        !!R for rebinning
        !!S for SAA possibility
        !!M for spacecraft Maneuver
        !!U for unused flags

        CHARACTER( LEN = 1  ), DIMENSION(16) :: bittype = &
           (/'E','E','W','E','W','W','U','R','W','W','S','M','W','U','U','U'/)
        INTEGER (KIND=4) :: di

        mqaL_error   = 0
        mqaL_warning = 0
        rebinF       = 0
        saaF         = 0
        maneuverF    = 0
        DO di = 1, 16
          SELECT CASE ( bittype(di) )
            CASE( 'E' )
              mqaL_error   = mqaL_error + int(IBITS( measurementQAflags, di, 1),kind(mqaL_error))
            CASE( 'W' )
              mqaL_warning = mqaL_warning + int(IBITS( measurementQAflags, di, 1),kind(mqaL_warning))
            CASE( 'R' )
              IF( PRESENT( rebinF ) ) rebinF = int(IBITS( measurementQAflags,di,1),kind(rebinF))
            CASE( 'S' )
              IF( PRESENT( saaF   ) ) saaF   = int(IBITS( measurementQAflags,di,1),kind(saaF))
            CASE( 'M' )
              IF(PRESENT(maneuverF))maneuverF= int(IBITS( measurementQAflags,di,1 ),kind(maneuverF))
            CASE DEFAULT
              CYCLE
           END SELECT
        ENDDO
        RETURN
      END FUNCTION L1B_extractMqa

!! 9. L1B_getIRRnmqa
       FUNCTION L1B_getIRRnmqa( irrLUN, bandname, irrFN, irrSN, &
                                irrE, irrW, irrA ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN) :: irrLUN
         CHARACTER( LEN = * ), INTENT(IN) :: bandname
         CHARACTER( LEN = * ), INTENT(OUT) :: irrFN, irrSN
         INTEGER (KIND=1), INTENT(OUT) :: irrE, irrW, irrA
         INTEGER (KIND=4) :: numfiles, nswath, strbufsize
         INTEGER (KIND = 4) :: rank
         INTEGER (KIND = 4) :: iLine, nL
         INTEGER (KIND = 4), DIMENSION(1:3) :: dims
         INTEGER (KIND = 2), DIMENSION(1) :: MeasurementQF
         INTEGER (KIND = 1) :: rebinF, saaF, maneuverF
         INTEGER (KIND=4) :: version
         INTEGER (KIND=4) :: status, ierr
         INTEGER (KIND=4) :: di, ii
         CHARACTER( LEN = PGSd_PC_FILE_PATH_MAX ) :: swathlist
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND = 4) :: getl1bdblk
         EXTERNAL              getl1bdblk

         !! Get number of Irradiance files from PCF
         !! and make sure there is only one file
         status = pgs_pc_getnumberoffiles( irrLUN, numfiles )
         IF( status /= PGS_S_SUCCESS ) THEN
            WRITE( msg,'(A,I9)' ) "can't get numfiles from PCF file at LUN =", &
                                 irrLUN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L1B_getIRRnmqa", zero )
            status = OZT_E_FAILURE
            RETURN
         ELSE IF( numfiles /= one ) THEN
            WRITE( msg,'(A,I9,A,I9)' ) "numfiles =", numfiles, "not equal to 1"//&
                                     " at LUN =", irrLUN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L1B_getIRRnmqa", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         !! Get the L1B Irradiance file name from the PCF.
         version = 1
         status = pgs_pc_getreference( irrLUN, version, irrFN )
         IF( status /= PGS_S_SUCCESS ) THEN
            WRITE( msg,'(A,I9)' ) "get filename from PCF file at LUN =", &
                                  irrLUN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L1B_getIRRnmqa", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         nswath = SWinqswath( irrFN, swathlist, strbufsize )
         IF( INDEX( swathlist, "Sun" ) == 0  .OR. &
             INDEX( swathlist, TRIM(bandname) ) == 0 ) THEN 
            WRITE( msg,'(A)' ) "Can't find Sun "//TRIM(bandname)// &
                               " Swath in" // TRIM( swathlist) // "file:"// &
                                TRIM( irrFN )
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L1B_getIRRnmqa", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         IF( nswath == 1 ) THEN
            irrSN = TRIM(swathlist)
         ELSE
            DO ii = 1, nswath 
              di = INDEX( swathlist, ',' ) 
              IF( di == 0 ) di = strbufsize + 1
              irrSN = swathlist(1:di-1)
              IF( INDEX( irrSN, TRIM(bandname) ) > 0 ) EXIT 
              swathlist = swathlist(di+1:strbufsize)
              strbufsize = strbufsize - (di-1)
            ENDDO
         ENDIF

         iLine = 0
         nL    = 1
         status = getl1bdblk( irrFN, irrSN, &
                             "MeasurementQualityFlags", DFNT_UINT16,&
                             iLine, nL, rank, dims, MeasurementQF )
         IF( status .NE. OZT_S_SUCCESS ) THEN
            ierr = OMI_SMF_setmsg( status, &
                                  "retrieve MeasurementQualityFlags",&
                                  "L1B_getIRRnmqa", zero )
            RETURN
         ENDIF

         irrE = L1B_extractMqa( MeasurementQF(1), irrW, &
                                rebinF, saaF, maneuverF )
         irrA = rebinF + saaF + maneuverF 

         status = OZT_S_SUCCESS
         RETURN
       END FUNCTION L1B_getIRRnmqa

!!10. L1B_selectIRR
       FUNCTION L1B_selectIRR( bandname, usedLUN, IRR_filename, IRR_swathname, &
                               fileType, normalIRRmissing, irrE, irrW, irrA ) &
                               RESULT( status )
         USE OMI_LUN_set
         CHARACTER( LEN = * ), INTENT(IN) :: bandname
         CHARACTER( LEN = * ), INTENT(OUT) :: IRR_filename,  &
                                              IRR_swathname, &
                                              fileType
         CHARACTER( LEN = PGSd_PC_FILE_PATH_MAX ) :: IRR_filenameB, &
                                                     IRR_swathnameB
         INTEGER( KIND=4 ), INTENT(OUT) :: usedLUN
         LOGICAL, INTENT(OUT) ::  normalIRRmissing
         INTEGER (KIND=1), INTENT(OUT) :: irrE, irrW, irrA
         INTEGER (KIND=1) :: irrEB, irrWB, irrAB
         INTEGER (KIND=4) :: status, status1


         status = L1B_getIRRnmqa( L1B_IRR_FILE_LUN, bandname,  &
                                  IRR_filename, IRR_swathname, &
                                  irrE, irrW, irrA )

         IF( status /= OZT_S_SUCCESS ) THEN
            normalIRRmissing = .TRUE.
            ierr = OMI_SMF_setmsg( OZT_W_GENERAL, &
                                  "Normal Solar Product Missing",&
                                  "L1B_selectIRR", zero )
                              
            status1 = L1B_getIRRnmqa( BACKUP_L1BIRR_LUN, bandname, &
                                      IRR_filenameB, IRR_swathnameB, &
                                      irrEB, irrWB, irrAB )
            IF( status1 /= OZT_S_SUCCESS ) THEN
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                          "No valid (Noraml or Backup) solar IRR file found.",&
                          "L1B_selectIRR", zero )
               RETURN
            ELSE IF( irrEB > 0 ) THEN
               ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                    "Missing Noraml L1BIRR and Measurment error for Backup .",&
                    "L1B_selectIRR", zero )
               status = OZT_E_FAILURE
               RETURN
            ELSE
               fileType = "Backup"
               usedLUN  = BACKUP_L1BIRR_LUN
            ENDIF
         ELSE 
            normalIRRmissing = .FALSE.
            IF( irrE > 0 .OR. irrW > 0 ) THEN
               status1 = L1B_getIRRnmqa( BACKUP_L1BIRR_LUN, bandname, &
                                      IRR_filenameB, IRR_swathnameB,  &
                                      irrEB, irrWB, irrAB )
               IF( status1 /= OZT_S_SUCCESS ) THEN
                  IF( irrE == 0 ) THEN
                     fileType = "Normal"
                     usedLUN  = L1B_IRR_FILE_LUN
                  ELSE
                    ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                    "Missing Backup L1BIRR and Measurment error for Normal .",&
                    "L1B_selectIRR", zero )
                    status = OZT_E_FAILURE
                    RETURN
                  ENDIF
               ELSE !! either has error or warning in normal product
                  IF( irrE > 0 ) THEN  ! error normal
                     IF( irrEB > 0 ) THEN  ! no backup
                        ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                              "Measurement error of both "// &
                                              "Normal and Backup OML1BIRR", &
                                              "L1B_selectIRR", zero )
                        status = OZT_E_FAILURE
                        RETURN
                     ELSE                     ! error normal,backup error free
                        fileType = "Backup"
                        usedLUN  = BACKUP_L1BIRR_LUN
                     ENDIF
                  ELSE   ! error free normal
                     IF( irrW > 0 ) THEN   ! normal warning
                        IF( irrEB > 0 ) THEN ! backup error, use normal
                           fileType = "Normal"
                           usedLUN  = L1B_IRR_FILE_LUN
                        ELSE                      ! backup error free
                           IF( irrW == 0 ) THEN ! backup warning free
                              fileType = "Backup"
                              usedLUN  = BACKUP_L1BIRR_LUN
                           ELSE                     !backup warning, use noraml
                              fileType = "Normal"
                              usedLUN  = L1B_IRR_FILE_LUN
                           ENDIF
                        ENDIF
                     ELSE ! normal warning free
                        fileType = "Normal"
                        usedLUN  = L1B_IRR_FILE_LUN
                     ENDIF
                  ENDIF
               ENDIF
            ELSE   !! no error or warning, use normal
               fileType = "Normal"
               usedLUN  = L1B_IRR_FILE_LUN
            ENDIF
         ENDIF

         status = L1B_getIRRnmqa( usedLUN, bandname, &
                                  IRR_filename, IRR_swathname, &
                                  irrE, irrW, irrA )

         RETURN
       END FUNCTION L1B_selectIRR

!!11. L1Bri_interpWL
      !! This function is a simple scheme to set the Precision and
      !! pixel quality flag at the input wavelength (wl_com).
      FUNCTION L1Bri_interpWL( Wavelength, PixelQualityFlags, &
                               wl_com, PixelQualityFlagsWL,   &
                               RadIrrPrecision_k, PrecisionWL_k ) &
                               RESULT (status)
        REAL (KIND = 4), DIMENSION(:,:), INTENT( IN ) :: Wavelength
        REAL (KIND = 4), DIMENSION(:), INTENT( IN ) :: wl_com
        INTEGER (KIND = 2), DIMENSION(:,:), INTENT( IN ) :: &
                                              PixelQualityFlags
        INTEGER (KIND = 2), DIMENSION(:,:), INTENT( OUT ) :: &
                                              PixelQualityFlagsWL
        REAL (KIND = 4), OPTIONAL, DIMENSION(:,:), INTENT( IN ) :: & 
                                         RadIrrPrecision_k
        REAL (KIND = 4), OPTIONAL, DIMENSION(:,:), INTENT( OUT ) :: & 
                                         PrecisionWL_k
        INTEGER (KIND = 4 ) :: iw, Nwl, Nwl_com, nXtrack, i, il, ih
        INTEGER (KIND = 4) :: status
        REAL (KIND = 4) :: frac

        call define_global_epsilon10

        status = OZT_S_SUCCESS
        Nwl_com = SIZE( wl_com        )
        Nwl     = SIZE( Wavelength, 1 )
        nXtrack = SIZE( Wavelength, 2 )

        IF( SIZE( PixelQualityFlags, 1 ) /= Nwl ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Nwl size not equal to " // &
                                  " PixelQualityFlags(:,*)",&
                                  "L1Bri_interpWL", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( SIZE( PixelQualityFlagsWL, 1 ) /= Nwl_com ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "Nwl_com size not " // &
                                  " equal to PixelQualityFlagsWL(:,*)",&
                                  "L1Bri_interpWL", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF
                      
        IF( SIZE( PixelQualityFlagsWL, 2 ) /= nXtrack .OR. &
            SIZE( PixelQualityFlags,   2 ) /= nXtrack ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "nXtrack size not good:" // &
                                 " PixelQualityFlags or PixelQualityFlagsWL",&
                                 "L1Bri_interpWL", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        DO iw = 1, Nwl_com
           il = iw
        DO i = 1, nXtrack
           il = hunt( Wavelength(1:Nwl,i ), wl_com(iw), il )
           IF( il == 0 .OR. il == Nwl ) THEN
              !! no measurement data for this wl, return 0 values
              !! and set PixelQualityFlags to be missing (bit 1 equal 1)
              PixelQualityFlagsWL(iw,i) = 1
              IF( PRESENT(PrecisionWL_k) ) PrecisionWL_k(iw,i) = 0
           ELSE
              ih = il+1
              frac = (wl_com(iw)       - Wavelength(il,i))/ &
                     (Wavelength(ih,i) - Wavelength(il,i))
              IF( ABS( frac ) < EPSILON10 ) THEN
                  ih = il
                  frac = 0.0
              ELSE IF( ABS( frac -1.0 ) < EPSILON10 ) THEN
                  il = ih
                  frac = 0.0
              ENDIF
 
              IF( ih == il ) THEN
                 PixelQualityFlagsWL(iw,i) = PixelQualityFlags(il,i)
                 IF( PRESENT(PrecisionWL_k) .AND. &
                     PRESENT(RadIrrPrecision_k) ) PrecisionWL_k(iw,i)  &
                                                = RadIrrPrecision_k(il,i)
              ELSE
                 PixelQualityFlagsWL(iw,i) = IOR( PixelQualityFlags(il,i), &
                                                  PixelQualityFlags(ih,i) )
                 IF( PRESENT(PrecisionWL_k) .AND. & 
                     PRESENT(RadIrrPrecision_k) ) THEN
                    PrecisionWL_k(iw,i)  = MAX( RadIrrPrecision_k(il,i), &
                                                RadIrrPrecision_k(ih,i))
                 ENDIF
              ENDIF
           ENDIF
        ENDDO 
        ENDDO

        RETURN
      END FUNCTION L1Bri_interpWL

!! 12. L1Bri_getLineWC
 !    This function gets the wavelength coefficients the reference column from
 !    the data block.
 !    this: the block data structure
 !    iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !           range from 0 to (nTimes-1) inclusive.
 !    wlCoef, nWavelCoef, and RefCol will be set by the the function.
 !
 !    status: the return PGS_SMF status value
!!
!      FUNCTION L1Bri_getLineWC( this, iLine, WlCoef, nWavelCoef, RefCol ) &
!                                RESULT (status)
!        TYPE (L1B_radirr_type), INTENT( INOUT ) :: this
!        INTEGER (KIND = 4), INTENT( IN ) :: iLine 
!        INTEGER (KIND = 4), INTENT( OUT ) :: nWavelCoef
!        INTEGER (KIND = 2), INTENT( OUT ) :: RefCol
!        REAL (KIND = 4), DIMENSION(:,:), INTENT( OUT ) :: WlCoef 
!
!        INTEGER :: i, j, k, di, ic
!        INTEGER (KIND = 4) :: Nwl_l
!        INTEGER (KIND = 4) :: status
!        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg
!
!        status = OZT_S_SUCCESS
!        IF( .NOT. this%initialized ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
!                                 "input block not initialized", &
!                                 "L1Bri_getLineWC", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!
!        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLine out of range", &
!                                 "L1Bri_getLineWC", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!
!        IF( iLine < this%iLine .OR. iLine > this%eLine ) THEN
!           IF( INDEX( this%swathname, "Earth" ) > 0 ) THEN
!              status = fill_rad_blk( this, iLine )
!           ELSE IF( INDEX( this%swathname, "Sun" ) > 0 ) THEN
!              status = fill_irr_blk( this, iLine )
!           ELSE
!              WRITE( msg,'(A)' ) "Unknow swathname:" // TRIM( this%swathname )
!              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L1Bri_getLineWC", zero )
!              status = OZT_E_FAILURE
!              RETURN
!           ENDIF
!
!           IF( status .NE. OZT_S_SUCCESS ) THEN
!              ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, "retrieve data block", &
!                                    "L1Bri_getLineWC", zero )
!              RETURN
!           ENDIF
!        ENDIF
!
!        nWavelCoef = this%nWavelCoef
!        IF( SIZE( WlCoef, 2 ) /= this%nXtrack .OR. &
!            SIZE( WlCoef, 1 ) /= nWavelCoef  ) THEN
!           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "WlCoef size not good", &
!                                 "L1Bri_getLineWC", zero )
!           status = OZT_E_FAILURE
!           RETURN
!        ENDIF
!
!        j = iLine - this%iLine + 1
!        RefCol = this%wRefCol(j)   
!        WlCoef(:,:) = this%wCof(:,:,j)
!        RETURN
!      END FUNCTION L1Bri_getLineWC

END MODULE L1B_radirr_class
