!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_radgeo_class
! 
! This module contains the allocatable arrays used to store in GEOlocation
! inforamtion and Earth view RADiance information retrived from the input
! L1B file. 
!
!!Revision History:
! Initial version 03/26/2002  Kai Yang/UMBC
!
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
!!Design Notes
!
!!END
!!****************************************************************************
MODULE O3T_radgeo_class
    USE PGS_PC_class ! define PGSd_PC_FILE_PATH_MAX, and pgs_pc functions 
                     ! include PGS_SMF.f define PGS_SMF_MAX_MSG_SIZE
    USE L1B_class    ! L1B_class contains module to read OMI L1B
                     ! geolocation, radiance, and irradiance parameters
    IMPLICIT NONE
    REAL (KIND = 8) :: time
    REAL (KIND = 4) :: sInD, scLat, scLon, scHgt
    REAL (KIND = 4), DIMENSION(:), ALLOCATABLE :: latitude, longitude, &
                                                  szenith, sazimuth, &
                                                  vzenith, vazimuth, &
                                                  phiArray, ptArray, pcArray, &
                                                  snowIceArray
    REAL (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: lat_bounds, lon_bounds
    LOGICAL(KIND=4), DIMENSION(:), ALLOCATABLE :: PclimQ
    INTEGER (KIND = 2), DIMENSION(:), ALLOCATABLE :: height
    INTEGER (KIND = 2), DIMENSION(:), ALLOCATABLE :: geoflg
    INTEGER (KIND = 1), DIMENSION(:), ALLOCATABLE :: anomflg
    REAL    (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: radiance, &
                                                       radPrecision, &
                                                       radWavelength
    INTEGER (KIND = 2), DIMENSION(:,:), ALLOCATABLE :: radQAflags
    INTEGER (KIND = 2) :: mqa_rad
    INTEGER (KIND = 1) :: instID_rad
    INTEGER (KIND = 4) :: nTimes_rad, nXtrack_rad, nWavel_rad
    INTEGER (KIND = 4) :: nWavelCoef_rad, nTimesSmallPixel_rad
    REAL (KIND = 4) :: EarthSunDistance
    TYPE (L1B_geoang_type) :: geo_blk
    TYPE (L1B_radirr_type) :: rad_blk

    type (l1b_radgeo_type) :: radgeo_blk

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four =4, five = 5

    PUBLIC  :: O3T_initRAD
    PUBLIC  :: O3T_freeRAD

    CONTAINS

 subroutine o3t_tio_initrad (this, filename, groupname, errstat)
   use l1b_tio_class
   use tell_module
   implicit none
   type (l1b_tio_type), intent(inout) :: this
   character (len=*), intent(in) :: filename, groupname
   integer, intent(inout) :: errstat

   integer :: ext, ierr
   character (len=len(filename)) :: file_name

   if (errstat < 0) return

   ! FIXME hack filename
   ext = index (filename, ".he4", back=.true.)
   if (ext > 0) then
     file_name = filename(1:ext)//"nc"
   else
     file_name = filename
   endif
   call l1b_tio_open (this, trim(file_name), groupname, errstat)
   if (errstat < 0) then
     call tell_error (tell_io_open_error, &
                      "o3t_tio_initrad: opening "//trim(file_name), &
                      errstat)
     return
   endif

   call l1b_tio_getdims (this, nTimes_rad, nXtrack_rad, nWavel_rad, errstat)
   nTimesSmallPixel_rad = -1  ! OMI-legacy: code will fail if this is set to 0
   call l1b_tio_earthsun_distance (this, EarthSunDistance, errstat)
   call l1b_tio_init_rad (this, radgeo_blk, errstat)
   if (errstat /= 0) then
     call tell_error (tell_runtime_error, "o3t_tio_initrad: initializing radiances", errstat)
     return
   endif

   allocate (latitude (nXtrack_rad), &
             longitude (nXtrack_rad), &
             lon_bounds (4, nXtrack_rad), &
             lat_bounds (4, nXtrack_rad), &
             szenith (nXtrack_rad), &
             sazimuth (nXtrack_rad), &
             vzenith (nXtrack_rad), &
             vazimuth (nXtrack_rad), &
             phiArray( nXtrack_rad ), &
             ptArray( nXtrack_rad ), &
             pcArray( nXtrack_rad ), &
             PclimQ( nXtrack_rad ), &
             snowIceArray( nXtrack_rad ), &
             height( nXtrack_rad ), &
             geoflg( nXtrack_rad ), &
             anomflg( nXtrack_rad ), &
             radiance (nWavel_rad, nXtrack_rad), &
             radprecision (nWavel_rad, nXtrack_rad), &
             radWavelength (nWavel_rad, nXtrack_rad), &
             radQAflags (nWavel_rad, nXtrack_rad), stat=ierr)
   if (ierr /= 0) then
     call tell_error (tell_malloc_error, "o3t_tio_initrad: allocate failed", errstat)
     return
   endif
   radprecision (:,:) = 0.0

 end subroutine
      
       FUNCTION O3T_initRAD( L1B_filename, L1B_swathname, wl_com ) &
                RESULT( status )
         USE OMI_LUN_set      ! define Logical Unit for Input and Output
         USE OMI_SMF_class    ! include PGE specific messages and OMI_SMF_setmsg
         CHARACTER( LEN = * ), INTENT(IN) :: L1B_filename, L1B_swathname
         REAL (KIND = 4), DIMENSION(:), INTENT(IN), OPTIONAL :: wl_com
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         !INTEGER (KIND=4) :: iLine
         !INTEGER (KIND=4) :: version
         INTEGER (KIND=4) :: status, ierr

         !! clean up everything before proceed
         CALL O3T_freeRAD
         EarthSunDistance = L1Bga_EarthSunDistance( L1B_filename,L1B_swathname )
         status = L1Bga_open( geo_blk, L1B_filename, L1B_swathname )
         IF( status /= OZT_S_SUCCESS ) THEN
            WRITE( msg,'(A)' ) "L1Bga_open "// TRIM(L1B_swathname) //&
                               " in file " // TRIM(L1B_filename) // " failed."
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_initRAD", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         status = L1Bri_open( rad_blk, L1B_filename, L1B_swathname )
         IF( status /= OZT_S_SUCCESS ) THEN
            WRITE( msg,'(A)' ) "L1Bri_open "// TRIM(L1B_swathname) //&
                               " in file " // TRIM(L1B_filename) // " failed."
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "O3T_initRAD", zero )
            status = OZT_E_FAILURE
            RETURN 
         ELSE
            WRITE( msg,'(A)' ) "L1Bri_open "// TRIM(L1B_swathname) //&
                               " in file " //TRIM(L1B_filename)//" succesfully."
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, "O3T_initRAD", four )
         ENDIF

         ! obtain sizes of dimensions defined in swath
         status = L1Bri_getSWdims( rad_blk, nTimes_rad, nXtrack_rad, &
                                   nWavel_rad, nWavelCoef_rad, &
                                   nTimesSmallPixel_rad )
         IF( status /= OZT_S_SUCCESS ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_FAILURE, &
                                "L1Bri_getSWdims failed.", "O3T_initRAD", zero )
            status = OZT_E_FAILURE
            RETURN
         ELSE
            WRITE( msg,'(A,4I4)' ) "(nTimes, nXtract, nWavel, nWavelCoef) =", &
                         nTimes_rad, nXtrack_rad, nWavel_rad, nWavelCoef_rad
            ierr = OMI_SMF_setmsg( OZT_S_SUCCESS, msg, "O3T_initRAD", four )
         ENDIF

         IF( PRESENT( wl_com ) ) THEN
            nWavel_rad = SIZE( wl_com )
            WRITE( msg,'(A,I3)' ) "Common array of fixed wavelength used, " //&
                               "nWavel_rad = ", nWavel_rad
            ierr = OMI_SMF_setmsg( OZT_W_GENERAL, msg, "O3T_initRAD", four )
         ENDIF
 
         ! allocate memory for arrays
         ALLOCATE( latitude( nXtrack_rad ), &
                   longitude( nXtrack_rad ), &
                   szenith( nXtrack_rad ), &
                   sazimuth( nXtrack_rad ), &
                   vzenith( nXtrack_rad ), &
                   vazimuth( nXtrack_rad ), &
                   phiArray( nXtrack_rad ), &
                   ptArray( nXtrack_rad ), &
                   pcArray( nXtrack_rad ), &
                   PclimQ( nXtrack_rad ), &
                   snowIceArray( nXtrack_rad ), &
                   height( nXtrack_rad ), &
                   geoflg( nXtrack_rad ), &
                   anomflg( nXtrack_rad ), &
                   radiance( nWavel_rad, nXtrack_rad ), &
                   radPrecision( nWavel_rad, nXtrack_rad ), &
                   radQAflags( nWavel_rad, nXtrack_rad ), &
                   radWavelength( nWavel_rad, nXtrack_rad ), & 
                   STAT=ierr )
         IF( ierr /= zero ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                  "radiance and geo allocation failure", &
                                  "O3T_initRAD", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF
       END FUNCTION O3T_initRAD

       SUBROUTINE O3T_freeRAD
         INTEGER (KIND=4) :: status
         IF( ALLOCATED( latitude      ) ) DEALLOCATE( latitude      )
         IF( ALLOCATED( longitude     ) ) DEALLOCATE( longitude     )
         IF( ALLOCATED( lat_bounds    ) ) DEALLOCATE( lat_bounds    )
         IF( ALLOCATED( lon_bounds    ) ) DEALLOCATE( lon_bounds    )
         IF( ALLOCATED( szenith       ) ) DEALLOCATE( szenith       )
         IF( ALLOCATED( sazimuth      ) ) DEALLOCATE( sazimuth      )
         IF( ALLOCATED( vzenith       ) ) DEALLOCATE( vzenith       )
         IF( ALLOCATED( vazimuth      ) ) DEALLOCATE( vazimuth      )
         IF( ALLOCATED( phiArray      ) ) DEALLOCATE( phiArray      )
         IF( ALLOCATED( ptArray       ) ) DEALLOCATE( ptArray       )
         IF( ALLOCATED( pcArray       ) ) DEALLOCATE( pcArray       )
         IF( ALLOCATED( PclimQ        ) ) DEALLOCATE( PclimQ        )
         IF( ALLOCATED( snowIceArray  ) ) DEALLOCATE( snowIceArray  )
         IF( ALLOCATED( height        ) ) DEALLOCATE( height        )
         IF( ALLOCATED( geoflg        ) ) DEALLOCATE( geoflg        )
         IF( ALLOCATED( anomflg       ) ) DEALLOCATE( anomflg       )
         IF( ALLOCATED( radiance      ) ) DEALLOCATE( radiance      )
         IF( ALLOCATED( radPrecision  ) ) DEALLOCATE( radPrecision  )
         IF( ALLOCATED( radWavelength ) ) DEALLOCATE( radWavelength )
         IF( ALLOCATED( radQAflags    ) ) DEALLOCATE( radQAflags    )

         ! close data block structure
         status = L1Bga_close( geo_blk )
         status = L1Bri_close( rad_blk )
       END SUBROUTINE O3T_freeRAD

END MODULE O3T_radgeo_class
