!! F90
 ! 
 !!Description:
 ! MODULE L1B_geoang_class contains the OMIL1B geolcation
 ! data structure and the functions for the creation, deletion and access 
 ! of specific information stored in the data structure. A user can use the
 ! following functions:  L1Bga_open, L1Bga_close, L1Bga_getSWdims, 
 ! and L1Bga_getLine.
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

MODULE L1B_geoang_class
   USE OMI_SMF_class
   USE UTIL_tools_class
   IMPLICIT NONE
   PUBLIC  :: L1Bga_open, L1Bga_close
   PUBLIC  :: L1Bga_getSWdims
   PUBLIC  :: L1Bga_getLine

   PRIVATE :: fill_geoang_blk
   INTEGER (KIND = 4), PARAMETER, PUBLIC :: MAX_NAME_LENGTH = 255
   INTEGER, PRIVATE :: ierr  !!error code returned from a function
   INTEGER, PARAMETER, PRIVATE :: zero = 0
   INTEGER, PARAMETER, PRIVATE :: one  = 1
   INTEGER, PARAMETER, PRIVATE :: two  = 2
!!
 ! L1B_geoang_type is designed to store a block of L1B swath data.
 ! A block contains all the geolocation data fields for all the
 ! "scan" lines or "exposure times" in a L1B swath. 
 ! filename: the file name for the L1B file (including path)
 ! swathname: the swath name in the L1B file, either "Earth UV-1 Swath",
 !            "Earth UV-2 Swath", or "Earth VIS Swath"
 ! nTimes: the total number of lines contained in the L1B swath
 ! nXtrack: the number of cross-track pixels in the swath
 ! initialized: the logical variable to indicate whether the
 !              block has been initialized.  The block data structure
 !              has to be created (initialized) before it can be used.
 !              
 ! Geolocation fields:
 !    lat: Latitude
 !    lon: Longitude
 !    sza: SolarZenithAngle
 !    saz: SolarAzimuthAngle
 !    vza: ViewingZenithAngle
 !    vaz: ViewingAzimuthAngle
 !    hgt: TerrainHeight
 !    flg: GroundPixelQualityFlags
 ! These are two-dimensional data arrays created for all the
 ! lines and the number of pixels for the data block when the 
 ! data structure is created.  A user cannot see these data fields
 ! directly but the information can be retrieved only through
 ! functions L1Bga_getLine, which return one line of the data 
 ! fields.
 ! 
!!

   TYPE, PUBLIC :: L1B_geoang_type
      PRIVATE
      CHARACTER ( LEN = MAX_NAME_LENGTH ) :: filename, swathname
      INTEGER (KIND = 4) :: nTimes
      INTEGER (KIND = 4) :: nXtrack
      LOGICAL :: initialized

      ! Geolocation fields
      REAL (KIND = 8), DIMENSION(:), POINTER :: tim
      REAL (KIND = 4), DIMENSION(:), POINTER :: sid, sla, sln, sat
      REAL (KIND = 4), DIMENSION(:,:), POINTER :: lon, lat, sza, saz, vza, vaz
      INTEGER (KIND = 2), DIMENSION(:,:), POINTER :: hgt
      INTEGER (KIND = 2), DIMENSION(:,:), POINTER :: flg
      INTEGER (KIND = 1), DIMENSION(:,:), POINTER :: anomflg
   END TYPE L1B_geoang_type

   CONTAINS

!! 1. L1Bga_open
 !    This function should be called first to initiate the
 !    the interface with the L1B swath. 
 !    this: the block data structure
 !    fn  : the L1B file name
 !    swn : the swath name in the L1B file
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bga_open( this, fn, swn ) RESULT (status)
        USE HE4_class
        TYPE (L1B_geoang_type), INTENT( INOUT ) :: this
        CHARACTER ( LEN = * ), INTENT( IN ) :: fn, swn
        INTEGER (KIND = 4) :: swfid, swid, status
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg
        
        !! open the L1B swath file
        swfid = swopen( fn, DFACC_READ )
        IF( swfid < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_FILE_OPEN, fn, "L1Bga_open", zero )
           RETURN
        ENDIF 

        !! attach to the swath
        swid = swattach( swfid, swn )
        IF( swid < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_SWATH_ATTACH, swn, "L1Bga_open", zero )
           RETURN
        ENDIF 
        
        !! retrieve the dimension info from the swath file.
        !! dimension names are obtained from file spec
        status = swrdattr( swid, "NumTimes", this%nTimes )
        IF( status < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get nTimes size failed", &
                                 "L1Bga_open", zero )
           RETURN
        ENDIF 
        IF( this%nTimes <= zero ) THEN
           WRITE( msg,'(A,I3,A)' ) "nTimes =", this%nTimes, &
                 ", No scan lines in "//TRIM(swn)// ", " // TRIM(fn) 
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, "L1Bga_open", zero )
           RETURN
        ENDIF 

        this%nXtrack    = swdiminfo( swid, "nXtrack" )
        IF( this%nXtrack < zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, "get xTrack size failed", &
                                  "L1Bga_open", zero )
           RETURN
        ENDIF 

        !! detach and close the L1B swath files.  No need for 
        !! error checking, for an error is unlikely to occur here.
        ierr = swdetach( swid )
        ierr = swclose( swfid )
          
        this%filename = fn   
        this%swathname = swn
        
        !! Allocate the memory for storage of geolocation 
        !! data fields.  First make sure pointers are not
        !! associated with any specific memory location.  If
        !! they are, deallocate the memory, then allocate the
        !! the proper amount of memory for each data field,
        !! error checking for each memory allocation to make
        !! sure memory is allocated successfully. 

        IF( ASSOCIATED( this%tim ) ) DEALLOCATE( this%tim )
        IF( ASSOCIATED( this%sid ) ) DEALLOCATE( this%sid )
        IF( ASSOCIATED( this%sla ) ) DEALLOCATE( this%sla )
        IF( ASSOCIATED( this%sln ) ) DEALLOCATE( this%sln )
        IF( ASSOCIATED( this%sat ) ) DEALLOCATE( this%sat )
        IF( ASSOCIATED( this%lon ) ) DEALLOCATE( this%lon )
        IF( ASSOCIATED( this%lat ) ) DEALLOCATE( this%lat )
        IF( ASSOCIATED( this%sza ) ) DEALLOCATE( this%sza )
        IF( ASSOCIATED( this%saz ) ) DEALLOCATE( this%saz )
        IF( ASSOCIATED( this%vza ) ) DEALLOCATE( this%vza )
        IF( ASSOCIATED( this%vaz ) ) DEALLOCATE( this%vaz )
        IF( ASSOCIATED( this%hgt ) ) DEALLOCATE( this%hgt )
        IF( ASSOCIATED( this%flg ) ) DEALLOCATE( this%flg )
        IF( ASSOCIATED( this%anomflg ) ) DEALLOCATE( this%anomflg )

        ALLOCATE( this%tim(this%nTimes),              &
                  this%sid(this%nTimes),              &
                  this%sla(this%nTimes),              &
                  this%sln(this%nTimes),              &
                  this%sat(this%nTimes),              &
                  this%lon(this%nXtrack,this%nTimes), &
                  this%lat(this%nXtrack,this%nTimes), &
                  this%sza(this%nXtrack,this%nTimes), &
                  this%saz(this%nXtrack,this%nTimes), &
                  this%vza(this%nXtrack,this%nTimes), &
                  this%vaz(this%nXtrack,this%nTimes), &
                  this%hgt(this%nXtrack,this%nTimes), &
                  this%flg(this%nXtrack,this%nTimes), &
                  this%anomflg(this%nXtrack,this%nTimes), STAT = ierr )
        IF( ierr .NE. zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                 "this: geolocation data block", &
                                 "L1Bga_open", two )
           RETURN
        ENDIF 

        this%initialized = .TRUE.
        status = fill_geoang_blk( this )

        RETURN      
      END FUNCTION L1Bga_open 

!! 2. L1Bga_getSWdims
 !    This function retrieves the dimension sizes from the L1B swath
 !    this: the block data structure
 !    nTimes: total number of lines in the swath
 !    nXtrack: number of cross-track pixels in the swath
!!
      FUNCTION L1Bga_getSWdims( this, nTimes_k, nXtrack_k ) RESULT( status )
        TYPE (L1B_geoang_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nTimes_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nXtrack_k
        INTEGER (KIND = 4) :: status

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Bga_getSWdims", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( PRESENT( nTimes_k     ) ) nTimes_k     = this%nTimes
        IF( PRESENT( nXtrack_k    ) ) nXtrack_k    = this%nXtrack
        status = OZT_S_SUCCESS
        RETURN
      END FUNCTION L1Bga_getSWdims

!! Private function: fill_geoang_blk
 !    This is a private function which is only used by other functions 
 !    in this MODULE to fill the block data structure. 
 !    this: the block data structure
 !    status: the return PGS_SMF status value
!!
      FUNCTION fill_geoang_blk( this ) RESULT( status ) 
        USE HE4_class
        TYPE (L1B_geoang_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4) :: status
        INTEGER (KIND = 4) :: getl1bdblk
        EXTERNAL              getl1bdblk
        INTEGER (KIND = 4) :: rank
        INTEGER (KIND = 4), DIMENSION(1:2) :: dims

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "input block not initialized", &
                                 "fill_geoang_blk", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        !! read geolocation
        status = getl1bdblk( this%filename, this%swathname, &
                               "Time", DFNT_FLOAT64, & 
                                0, this%nTimes, rank, dims, &
                                this%tim )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve Time", "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "SecondsInDay", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%sid )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                               "retrieve SecondsInDay", "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "SpacecraftLatitude", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%sla )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                        "retrieve SpacecraftLatitude", "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "SpacecraftLongitude", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%sln )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                        "retrieve SpacecraftLongitude", "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "SpacecraftAltitude", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%sat )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                        "retrieve SpacecraftAltitude", "fill_geoang_blk", one )
           RETURN
        ENDIF
 
        status = getl1bdblk( this%filename, this%swathname, &
                               "Latitude", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%lat )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve Latitude", "fill_geoang_blk", one )
           RETURN
        ENDIF
 
        status = getl1bdblk( this%filename, this%swathname, &
                               "Longitude", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%lon )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, &
                                 "retrieve Longitude", "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "SolarZenithAngle", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%sza )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve SolarZenithAngle", &
                                 "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "SolarAzimuthAngle", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%saz )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve SolarAzimuthAngle", &
                                 "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "ViewingZenithAngle", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%vza )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve ViewingZenithAngle", & 
                                  "fill_geoang_blk", one )
           RETURN
        ENDIF
        
        status = getl1bdblk( this%filename, this%swathname, &
                               "ViewingAzimuthAngle", DFNT_FLOAT32, & 
                                0, this%nTimes, rank, dims, &
                                this%vaz )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve ViewingAzimuthAngle", &
                                 "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "TerrainHeight", DFNT_INT16, & 
                                0, this%nTimes, rank, dims, &
                                this%hgt )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve TerrainHeight", &
                                 "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "GroundPixelQualityFlags", DFNT_UINT16, & 
                                0, this%nTimes, rank, dims, &
                                this%flg )
        IF( status .NE. OZT_S_SUCCESS ) THEN
           ierr = OMI_SMF_setmsg( status, "retrieve GroundPixelQualityFlags", &
                                 "fill_geoang_blk", one )
           RETURN
        ENDIF

        status = getl1bdblk( this%filename, this%swathname, &
                               "XTrackQualityFlags", DFNT_UINT8, & 
                                0, this%nTimes, rank, dims, &
                                this%anomflg )

        ! status not equal to OZT_S_SUCCESS taken to indicate swath has no XTrackQualityFlags. don't abort.
        IF( status .NE. OZT_S_SUCCESS ) THEN

           ! print warning to log file
           ierr = OMI_SMF_setmsg( OZT_W_GENERAL, "Retrieval of XTrackQualityFlags from L1B failed.  All remain unset.", &
                                 "fill_geoang_blk", one )
           ! set XTrackAnomalyFlags all to zero and reset status to OZT_S_SUCCESS so code can continue
           this%anomflg(:,:) = 0
           status = OZT_S_SUCCESS
        ENDIF

        RETURN ! function return
      END FUNCTION fill_geoang_blk

!! 3. L1Bga_getLine
 !    This function gets one or more geolocation data field values from
 !    the data block.
 !    this: the block data structure
 !    iLine : the line number in the L1B swath.  NOTE: this input is 0 based
 !            range from 0 to (nTimes-1) inclusive.
 !    Latitude_k, Longitude_k, SolarZenith_k, SolarAzimuth_k, ViewZenith_k,
 !    ViewAzimuth_k, Height_k and GeoFlag_k are keyword arguments.
 !    Only those present in the argument list will be set by the
 !    the function.
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bga_getLine( this, iLine, Time_k, & 
                              Latitude_k, Longitude_k, &
                              solarZenith_k, SolarAzimuth_k, &
                              ViewZenith_k, ViewAzimuth_k, &
                              Height_k, Geoflag_k, XTrackQualityFlags_k, SecondsInDay_k, &
                              SpacecraftLat_k, SpacecraftLon_k, &
                              SpacecraftAltitude_k ) RESULT (status)

        TYPE (L1B_geoang_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4), INTENT( IN ) :: iLine
        REAL (KIND=8), OPTIONAL, INTENT( OUT ) :: Time_k
        REAL (KIND=4), OPTIONAL, INTENT( OUT ) :: SecondsInDay_k, &
               SpacecraftLat_k, SpacecraftLon_k, SpacecraftAltitude_k  
        REAL (KIND=4), OPTIONAL, DIMENSION(:), INTENT( OUT ) :: &
                                            Latitude_k, Longitude_k, &
                                            SolarZenith_k, SolarAzimuth_k, &
                                            ViewZenith_k, ViewAzimuth_k
        INTEGER (KIND=2), OPTIONAL, DIMENSION(:), INTENT( OUT ) :: Height_k
        INTEGER (KIND=2), OPTIONAL, DIMENSION(:), INTENT( OUT ) :: GeoFlag_k
        INTEGER (KIND=1), OPTIONAL, DIMENSION(:), INTENT( OUT ) :: XTrackQualityFlags_k
        INTEGER :: i
        INTEGER (KIND=4) :: status

        status = OZT_S_SUCCESS

        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "input block not initialized", &
                                 "L1Bga_getLine", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < 0 .OR. iLine >= this%nTimes ) THEN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLine out of range", &
                                 "L1Bga_getLine", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        i = iLine + 1

        IF( PRESENT( Time_k ) ) THEN
           Time_k = this%tim( i )
        ENDIF

        IF( PRESENT( SecondsInDay_k ) ) THEN
           SecondsInDay_k = this%sid( i )
        ENDIF

        IF( PRESENT( SpacecraftLat_k ) ) THEN
           SpacecraftLat_k = this%sla( i )
        ENDIF

        IF( PRESENT( SpacecraftLon_k ) ) THEN
           SpacecraftLon_k = this%sln( i )
        ENDIF

        IF( PRESENT( SpacecraftAltitude_k ) ) THEN
           SpacecraftAltitude_k = this%sat( i )
        ENDIF

        IF( PRESENT( Latitude_k ) ) THEN
           IF( SIZE( Latitude_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input latitude array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           Latitude_k = this%lat( 1:this%nXtrack, i )
        ENDIF
     
        IF( PRESENT( Longitude_k ) ) THEN
           IF( SIZE( Longitude_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input longitude array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           Longitude_k = this%lon( 1:this%nXtrack, i )
        ENDIF
            
        IF( PRESENT( SolarZenith_k ) ) THEN
           IF( SIZE( SolarZenith_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input solarZenith array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           SolarZenith_k = this%sza( 1:this%nXtrack, i )
        ENDIF

        IF( PRESENT( SolarAzimuth_k ) ) THEN
           IF( SIZE( SolarAzimuth_k ) <  this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input solarAzimuth array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           SolarAzimuth_k= this%saz( 1:this%nXtrack, i )
        ENDIF

        IF( PRESENT( ViewZenith_k ) ) THEN
           IF( SIZE( ViewZenith_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input viewZenith array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           ViewZenith_k  = this%vza( 1:this%nXtrack, i )
        ENDIF

        IF( PRESENT( ViewAzimuth_k ) ) THEN
           IF( SIZE( ViewAzimuth_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input viewAzimuth array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           ViewAzimuth_k = this%vaz( 1:this%nXtrack, i )
        ENDIF

        IF( PRESENT( Height_k ) ) THEN
           IF( SIZE( Height_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input height array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           Height_k  = this%hgt( 1:this%nXtrack, i )
        ENDIF

        IF( PRESENT( GeoFlag_k ) ) THEN
           IF( SIZE( GeoFlag_k ) < this%nXtrack  ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input geolocation flag array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           GeoFlag_k = this%flg( 1:this%nXtrack, i )
        ENDIF

        IF( PRESENT( XTrackQualityFlags_k )) THEN
           IF( SIZE( XTrackQualityFlags_k ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                    "input cross-track quality flags array too small", &
                                    "L1Bga_getLine", zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF
           XTrackQualityFlags_k = this%anomflg(1:this%nXtrack,i)
        ENDIF

        RETURN

      END FUNCTION L1Bga_getLine

!! 4. L1Bga_close
 !    This function should be called when the data block is no longer
 !    needed.  It deallocates all the allocated memory, and sets
 !    all the parameters to invalid values.
 !    this: the block data structure
 !    
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Bga_close( this ) RESULT( status )
        TYPE (L1B_geoang_type), INTENT( INOUT ) :: this
        INTEGER (KIND = 4) :: status
        status = OZT_S_SUCCESS

        !! DEALLOCATE all the memory
        IF( ASSOCIATED( this%tim ) ) DEALLOCATE( this%tim )
        IF( ASSOCIATED( this%sid ) ) DEALLOCATE( this%sid )
        IF( ASSOCIATED( this%sla ) ) DEALLOCATE( this%sla )
        IF( ASSOCIATED( this%sln ) ) DEALLOCATE( this%sln )
        IF( ASSOCIATED( this%sat ) ) DEALLOCATE( this%sat )
        IF( ASSOCIATED( this%lon ) ) DEALLOCATE( this%lon )
        IF( ASSOCIATED( this%lat ) ) DEALLOCATE( this%lat )
        IF( ASSOCIATED( this%sza ) ) DEALLOCATE( this%sza )
        IF( ASSOCIATED( this%saz ) ) DEALLOCATE( this%saz )
        IF( ASSOCIATED( this%vza ) ) DEALLOCATE( this%vza )
        IF( ASSOCIATED( this%vaz ) ) DEALLOCATE( this%vaz )
        IF( ASSOCIATED( this%hgt ) ) DEALLOCATE( this%hgt )
        IF( ASSOCIATED( this%flg ) ) DEALLOCATE( this%flg )
        IF( ASSOCIATED( this%anomflg ) ) DEALLOCATE( this%anomflg )

        this%nTimes      = -1
        this%nXtrack     = -1
        this%filename    = ""
        this%swathname   = ""
        this%initialized = .FALSE.
        RETURN
      END FUNCTION L1Bga_close

!! 5. L1Bga_EarthSunDistance
      FUNCTION L1Bga_EarthSunDistance( he4filename, swathname ) &
                                       RESULT( ESdistance )
        USE HE4_class
        CHARACTER (LEN = *), INTENT(IN) :: he4filename, swathname
        REAL (KIND = 4) :: ESdistance
        INTEGER (KIND = 4) :: swfid, swid, status, ierr 
        CHARACTER (LEN = MAX_NAME_LENGTH ) :: msg

        swfid = swopen( he4filename, DFACC_READ )
        IF( swfid /=  -1) THEN
           swid = swattach( swfid, swathname )
           IF( swid /= -1 ) THEN
              status = swrdattr( swid, "EarthSunDistance", ESdistance )
              IF( status == -1 ) THEN
                 WRITE(msg,'(A)') "Get swath attribute EarthSunDistance"// &
                                  "failed from "//TRIM(swathname)//","  // &
                                  TRIM( he4filename )  
                 ierr = OMI_SMF_setmsg( OZT_E_INPUT,  msg, &
                                       "L1Bga_EarthSunDistance", zero )
                 ESdistance = -1.0
              ENDIF
              status = swdetach(swid)
           ENDIF
           status = swclose(swfid)
        ENDIF

      END FUNCTION L1Bga_EarthSunDistance

END MODULE L1B_geoang_class
