!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE O3T_L2output_class
! 
!  contains functions to allocate memory space for storing L2 results and 
!  transfer these results the output block and set fill values in case there
!  are no L2 output results. Thses functions are specific to the OMTO3 and
!  when output parameter list are changed, these functions need to be updated
!  accordingly.
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
! 
!!Return
! None 
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
MODULE O3T_L2output_class
    USE OMI_SMF_class ! include PGE specific messages and OMI_SMF_setmsg
    USE O3T_radgeo_class
    USE OMI_L2writer_class
    USE O3T_omto3_fs
    USE O3T_pixel_class

    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four =4, five = 5

    REAL    (KIND = 4), DIMENSION(:), ALLOCATABLE :: xnvalm, res_stp1, &
                                                     res_stp2, res_stp3, &
                                                     dndomega_t, dNdT, &
                                                     dndr
    INTEGER (KIND = 2), DIMENSION(:,:), ALLOCATABLE :: radQAflags_com, &
                                                       irrQAflags_com
    REAL    (KIND = 4), DIMENSION(:,:), ALLOCATABLE :: radPrecision_com, &
                                                       irrPrecision_com

    PUBLIC :: O3T_initL2out
    PUBLIC :: O3T_L2setGeoLine
    PUBLIC :: O3T_L2setDataPix
    PUBLIC :: O3T_L2fillDataPix
    PUBLIC :: O3T_L2setMqaLine
    PUBLIC :: O3T_freeL2out

    CONTAINS
       FUNCTION O3T_initL2out( wl_com ) RESULT( status )
         REAL (KIND = 4), DIMENSION(:), INTENT(IN), OPTIONAL :: wl_com
         !CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND=4) :: status, ierr
         INTEGER (KIND=4) :: nwl_com

         status = OZT_S_SUCCESS
         IF( PRESENT( wl_com ) ) THEN
            nwl_com = SIZE( wl_com )
         ELSE
            nwl_com = nWavel_rad
         ENDIF

         CALL O3T_freeL2out
         ALLOCATE( radQAflags_com( nwl_com, nXtrack_rad ), &
                   irrQAflags_com( nwl_com, nXtrack_rad ), &
                   radPrecision_com( nwl_com, nXtrack_rad ), &
                   irrPrecision_com( nwl_com, nXtrack_rad ), &
                   xnvalm( nwl_com ), &
                   res_stp1( nwl_com ), &
                   res_stp2( nwl_com ), &
                   res_stp3( nwl_com ), &
                   dndomega_t( nwl_com ), &
                   dNdT( nwl_com ), &
                   dndr( nwl_com ), &
                   STAT=ierr )
        IF( ierr /= zero ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, &
                                  "L2out allocation failure", &
                                  "O3T_initL2out", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF
         RETURN
 
       END FUNCTION O3T_initL2out

       SUBROUTINE O3T_freeL2out
         IF( ALLOCATED( radQAflags_com   ) ) DEALLOCATE( radQAflags_com   )
         IF( ALLOCATED( irrQAflags_com   ) ) DEALLOCATE( irrQAflags_com   )
         IF( ALLOCATED( radPrecision_com ) ) DEALLOCATE( radPrecision_com )
         IF( ALLOCATED( irrPrecision_com ) ) DEALLOCATE( irrPrecision_com )
         IF( ALLOCATED( xnvalm           ) ) DEALLOCATE( xnvalm           )
         IF( ALLOCATED( res_stp1         ) ) DEALLOCATE( res_stp1         )
         IF( ALLOCATED( res_stp2         ) ) DEALLOCATE( res_stp2         )
         IF( ALLOCATED( res_stp3         ) ) DEALLOCATE( res_stp3         )
         IF( ALLOCATED( dndomega_t       ) ) DEALLOCATE( dndomega_t       )
         IF( ALLOCATED( dNdT             ) ) DEALLOCATE( dNdT             )
         IF( ALLOCATED( dndr             ) ) DEALLOCATE( dndr             )
       END SUBROUTINE O3T_freeL2out

       SUBROUTINE O3T_L2setGeoLine( iT, geoblk, datablk ) 
         TYPE (L2_generic_type), INTENT( INOUT ) :: geoblk, datablk
         INTEGER (KIND=4), INTENT(IN) :: iT
         INTEGER (KIND=4) :: ig, id, Ls, Le, bsize
         !REAL (KIND=4), DIMENSION(nXtrack_rad) :: R4Array

         DO ig = 1, geoblk%nFields
           bsize = geoblk%lineSize(ig)
           Ls    = geoblk%accuBlkSize(ig-1) + (iT-1)*bsize
           Le    = Ls + bsize
           Ls    = Ls + 1       !! fortran index scheme

           IF(      ig == 1 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( geoflg,    I1, bsize )
           ELSE IF( ig == 2 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( latitude,  I1, bsize )
           ELSE IF( ig == 3 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( longitude, I1, bsize )
           ELSE IF( ig == 4 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( szenith,   I1, bsize )
           ELSE IF( ig == 5 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( sazimuth,  I1, bsize )
           ELSE IF( ig == 6 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( vzenith,   I1, bsize )
           ELSE IF( ig == 7 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( vazimuth,  I1, bsize )
           ELSE IF( ig == 8 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( phiArray,  I1, bsize )
           ELSE IF( ig == 9 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( height,    I1, bsize )
           ELSE IF( ig == 10 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( time,      I1, bsize )
           ELSE IF( ig == 11 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( sInD,      I1, bsize )
           ELSE IF( ig == 12 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( scLat,     I1, bsize )
           ELSE IF( ig == 13 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( scLon,     I1, bsize )
           ELSE IF( ig == 14 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( scHgt,     I1, bsize )
           ELSE IF( ig == 15 ) THEN
              geoblk%data( Ls:Le ) = TRANSFER( anomflg,   I1, bsize)
           ENDIF
         ENDDO

         DO id = 1, 2
           bsize = datablk%lineSize(id)
           Ls    = datablk%accuBlkSize(id-1) + (iT-1)*bsize
           Le    = Ls + bsize
           Ls    = Ls + 1       !! fortran index scheme
           IF( id == 1 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( pcArray,  I1, bsize )
           ELSE IF( id == 2 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( ptArray,  I1, bsize )
           ENDIF
         ENDDO

       END SUBROUTINE O3T_L2setGeoLine
         
       SUBROUTINE O3T_L2setDataPix( iT, iX, nWavel, nLayers, &
                                    algflg, QAflags, radBadPixflgs, & 
                                    stp1oz, stp2oz, stp3oz, &
                                    oz_cld, aerind, so2ind, & 
                                    pixSURF, eff, aprfoz, datablk )
         TYPE (L2_generic_type), INTENT( INOUT ) :: datablk
         INTEGER (KIND=4), INTENT(IN) :: iT, iX, nWavel, nLayers
         TYPE (O3T_pixcover_type), INTENT(IN)  :: pixSURF
         REAL (KIND=4), DIMENSION(:), INTENT(IN) :: eff, aprfoz
         INTEGER (KIND=1), INTENT(IN) :: algflg
         INTEGER (KIND=2), INTENT(IN) :: QAflags, radBadPixflgs
         REAL (KIND=4), INTENT(IN) :: stp1oz, stp2oz, stp3oz, &
                                      oz_cld, aerind, so2ind 

         INTEGER (KIND=4) :: id, Ls, Le, bsize ! , is, ie
         INTEGER (KIND=4) :: LL
         !INTEGER (KIND=1) :: I1temp
         !INTEGER (KIND=2) :: I2temp
         REAL    (KIND=4) :: R4temp
         REAL    (KIND=4), DIMENSION( nWavel  ) :: R4wlArray

         DO id = 3, 24
           LL    = (iT-1)*nXtrack_rad + (iX-1)     !! 0 based index
           bsize = datablk%pixSize(id)
           Ls    = datablk%accuBlkSize(id-1) + LL*bsize
           Le    = Ls + bsize
           Ls    = Ls + 1       !! fortran index scheme
           IF(      id == 3 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( algflg,        I1, bsize )
           ELSE IF( id == 4 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( QAflags,       I1, bsize )
           ELSE IF( id == 5 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( radBadPixflgs, I1, bsize )
           ELSE IF( id == 6 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( pixSURF%clfrac,I1, bsize )
           ELSE IF( id == 7 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( pixSURF%rcf1  ,I1, bsize )
!!            datablk%data( Ls:Le ) = TRANSFER( pixSURF%rcf2  ,I1, bsize )
           ELSE IF( id == 8 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( stp3oz,        I1, bsize )
           ELSE IF( id == 9 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( oz_cld,        I1, bsize )
           ELSE IF( id == 10) THEN
              R4temp = 100.0*pixSURF%ref 
              datablk%data( Ls:Le ) = TRANSFER( R4temp,        I1, bsize )
           ELSE IF( id == 11 ) THEN
              R4temp = 100.0*pixSURF%ref360
              datablk%data( Ls:Le ) = TRANSFER( R4temp,        I1, bsize )
           ELSE IF( id == 12 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( so2ind,        I1, bsize )
           ELSE IF( id == 13 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( stp1oz,        I1, bsize )
           ELSE IF( id == 14 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( stp2oz,        I1, bsize )
           ELSE IF( id == 15 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( aerind,        I1, bsize )
           ELSE IF( id == 16 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( dndr(:),       I1, bsize )
           ELSE IF( id == 17 ) THEN
              R4wlArray(:) = 100.0*xnvalm(:)
              datablk%data( Ls:Le ) = TRANSFER( R4wlArray,     I1, bsize )
           ELSE IF( id == 18 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( res_stp3(:),   I1, bsize )
           ELSE IF( id == 19 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( res_stp1(:),   I1, bsize )
           ELSE IF( id == 20 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( res_stp2(:),   I1, bsize )
           ELSE IF( id == 21 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( dndomega_t(:), I1, bsize )
           ELSE IF( id == 22 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( dNdT(:),       I1, bsize )
           ELSE IF( id == 23 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( aprfoz(:),     I1, bsize )
           ELSE IF( id == 24 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( eff(:),        I1, bsize )
           ENDIF
         ENDDO 
         
       END SUBROUTINE O3T_L2setDataPix

       SUBROUTINE O3T_L2fillDataPix( iT, iX, nWavel, nLayers, &
                                     algflg, QAflags, radBadPixflgs, datablk )
         TYPE (L2_generic_type), INTENT( INOUT ) :: datablk
         INTEGER (KIND=4), INTENT(IN) :: iT, iX, nWavel, nLayers
         INTEGER (KIND=1), INTENT(IN) :: algflg
         INTEGER (KIND=2), INTENT(IN) :: QAflags, radBadPixflgs
         INTEGER (KIND=4) :: id, is, ie, Ls, Le, bsize
         INTEGER (KIND=4) :: LL
         !INTEGER (KIND=1) :: I1temp
         !INTEGER (KIND=2) :: I2temp
         REAL    (KIND=4) :: R4temp
         REAL    (KIND=4), DIMENSION( nWavel  ) :: R4wlArray
         REAL    (KIND=4), DIMENSION( nLayers ) :: R4lyrArray

         is = 3
         ie = 15
         DO id = is, ie
           LL    = (iT-1)*nXtrack_rad + (iX-1)     !! 0 based index
           bsize = datablk%pixSize(id)
           Ls    = datablk%accuBlkSize(id-1) + LL*bsize
           Le    = Ls + bsize
           Ls    = Ls + 1       !! fortran index scheme
           IF(      id == 3 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( algflg,  I1, bsize )
           ELSE IF( id == 4 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( QAflags, I1, bsize )
           ELSE IF( id == 5 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( radBadPixflgs, I1, bsize )
           ELSE 
              R4temp = fill_float32
              datablk%data( Ls:Le ) = TRANSFER( R4temp,  I1, bsize )
           ENDIF
         ENDDO

         is = 16
         ie = 22
         R4wlArray(:) = fill_float32
         DO id = is, ie   !! 3-D datafields: nWavel,nXtrack,nTimes
           LL    = (iT-1)*nXtrack_rad + (iX-1)     !! 0 based index
           bsize = datablk%pixSize(id)
           Ls    = datablk%accuBlkSize(id-1) + LL*bsize
           Le    = Ls + bsize
           Ls    = Ls + 1       !! fortran index scheme
           IF( id == is+1 ) THEN
              datablk%data( Ls:Le ) = TRANSFER( 100.0*xnvalm(:), I1, bsize )
           ELSE
              datablk%data( Ls:Le ) = TRANSFER( R4wlArray,  I1, bsize )
           ENDIF
         ENDDO 

         is = 23
         ie = 24
         R4lyrArray(:) = fill_float32
         DO id = is, ie   !! 3-D datafields: nLayers,nXtrack,nTimes
           LL    = (iT-1)*nXtrack_rad + (iX-1)     !! 0 based index
           bsize = datablk%pixSize(id)
           Ls    = datablk%accuBlkSize(id-1) + LL*bsize
           Le    = Ls + bsize
           Ls    = Ls + 1       !! fortran index scheme
           datablk%data( Ls:Le ) = TRANSFER( R4lyrArray,  I1, bsize )
         ENDDO 

       END SUBROUTINE O3T_L2fillDataPix

       SUBROUTINE O3T_L2setMqaLine( iT, mqaL2, datablk ) 
         TYPE (L2_generic_type), INTENT( INOUT ) :: datablk
         INTEGER (KIND=1), INTENT(IN) :: mqaL2
         INTEGER (KIND=4), INTENT(IN) :: iT
         INTEGER (KIND=4) :: id, Ls, Le, bsize
         DO id = 25, 25
           bsize = datablk%lineSize(id)
           Ls    = datablk%accuBlkSize(id-1) + (iT-1)*bsize
           Le    = Ls + bsize
           Ls    = Ls + 1       !! fortran index scheme
           datablk%data( Ls:Le ) = TRANSFER( mqaL2,  I1, bsize )
         ENDDO       
         
       END SUBROUTINE O3T_L2setMqaLine
END MODULE O3T_L2output_class
