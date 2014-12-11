MODULE UTIL_lh4_class
    USE HDF4_class
    USE OMI_SMF_class
    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: MAXRANK = 7, zero = 0
   
    TYPE, PUBLIC :: SDS_T
      INTEGER (KIND=4) :: id, rank
      INTEGER, DIMENSION(MAXRANK) :: dims
      CHARACTER (LEN = 80 ), DIMENSION(MAXRANK) :: dims_name
      CHARACTER (LEN = 80 ) :: name
      CHARACTER (LEN =160)  :: long_name 
      INTEGER (KIND=4) :: data_type
      CHARACTER (LEN = 80 ) :: label,unit,format,coordsys
      REAL (KIND = 8 ) :: cal, cal_err, offset, offset_err
      INTEGER (KIND=4) :: data_type_uncal
      REAL (KIND = 8 ) :: fill_value, max, min
    END TYPE SDS_T

    INTERFACE UTIL_lh4_get
      MODULE PROCEDURE UTIL_lh4_get1DIK2
      MODULE PROCEDURE UTIL_lh4_get2DIK2
      MODULE PROCEDURE UTIL_lh4_get3DIK2
      MODULE PROCEDURE UTIL_lh4_get4DIK2
      MODULE PROCEDURE UTIL_lh4_get5DIK2
      MODULE PROCEDURE UTIL_lh4_get6DIK2
    END INTERFACE

    INTERFACE UTIL_lh4_put
      MODULE PROCEDURE UTIL_lh4_put1DIK2
      MODULE PROCEDURE UTIL_lh4_put2DIK2
      MODULE PROCEDURE UTIL_lh4_put3DIK2
      MODULE PROCEDURE UTIL_lh4_put4DIK2
      MODULE PROCEDURE UTIL_lh4_put5DIK2
      MODULE PROCEDURE UTIL_lh4_put6DIK2
      MODULE PROCEDURE UTIL_lh4_put1DRK4
      MODULE PROCEDURE UTIL_lh4_put2DRK4
      MODULE PROCEDURE UTIL_lh4_put3DRK4
      MODULE PROCEDURE UTIL_lh4_put4DRK4
      MODULE PROCEDURE UTIL_lh4_put5DRK4
      MODULE PROCEDURE UTIL_lh4_put6DRK4
    END INTERFACE

    CONTAINS
      FUNCTION UTIL_lh4_selectSDS( sd_fid, sds ) RESULT (status)
        INTEGER , INTENT( IN )  :: sd_fid    !Science data file id
        TYPE (SDS_T), INTENT( INOUT ) :: sds !Data structure defining the SDS
        INTEGER (KIND = 4) :: status
        INTEGER          :: sds_index, rank, data_type, num_attrs 
        CHARACTER (LEN=80 ):: sds_name 
        INTEGER, DIMENSION(MAXRANK) :: dims
        INTEGER          :: ierr
        CHARACTER (LEN =255) :: msg

        status = OZT_S_SUCCESS

        sds_index = sfn2index( sd_fid, sds%name )
        IF( sds_index == -1 ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "sfn2index failed on dataset: ", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh4_selectSDS", zero )
           RETURN
        ENDIF

        sds%id = sfselect( sd_fid, sds_index ) 
        IF( sds%id == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "index-to-sds_id failed on dataset: ", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh4_selectSDS", zero )
           RETURN
        ENDIF

        !get the rank and dim_sizes
        ierr = sfginfo(sds%id, sds_name, rank, dims, data_type, num_attrs) 
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "get sds info failed on dataset: ", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh4_selectSDS", zero )
           RETURN
        ENDIF

        IF( sds%data_type .NE. -1 ) THEN
           IF( data_type .NE. sds%data_type ) THEN
              status = OZT_E_FAILURE
              WRITE( msg,* ) "data type not matched: ", sds%name
              ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                     "UTIL_lh4_selectSDS", zero )
              RETURN
           ENDIF
        ENDIF

        IF( sds%rank .NE. -1 ) THEN
           IF( rank .NE. sds%rank ) THEN
              status = OZT_E_FAILURE
              WRITE( msg,* ) "rank not matched: ", sds%name
              ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                     "UTIL_lh4_selectSDS", zero )
              RETURN
           ENDIF
        ENDIF

        sds%dims(1:rank) = dims(1:rank)

        RETURN
      END FUNCTION UTIL_lh4_selectSDS

      FUNCTION UTIL_lh4_createSDS( sd_fid, sds ) RESULT (status)
        INTEGER , INTENT( IN )  :: sd_fid    !Science data file id
        TYPE (SDS_T), INTENT( INOUT ) :: sds !Data structure defining the SDS
        INTEGER (KIND = 4) :: status
        INTEGER          :: rank, data_type, num_attrs
        INTEGER, DIMENSION(MAXRANK) :: dims 
        INTEGER          :: ierr, dim_id
        CHARACTER (LEN =255) :: msg
        INTEGER :: di, count 
        INTEGER(KIND=1) :: fill_int8, max_int8, min_int8
        INTEGER(KIND=2) :: fill_int16, max_int16, min_int16
        INTEGER(KIND=4) :: fill_int32, max_int32, min_int32
        REAL(KIND=4) :: fill_f32, max_f32, min_f32
        status = OZT_S_SUCCESS

        dims(1:sds%rank) = sds%dims(1:sds%rank)
        sds%id = sfcreate(sd_fid, sds%name, sds%data_type, sds%rank, dims ) 
        IF( sds%id == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "sfcreate failed on dataset: ", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_createSDS", zero )
           RETURN
        ENDIF

        DO di = 1, sds%rank 
          dim_id = sfdimid(sds%id, di-1)
          ierr   = sfsdmname(dim_id, sds%dims_name(di) )
        ENDDO

        ierr = sfscal( sds%id, sds%cal, sds%cal_err, sds%offset, &
                       sds%offset_err, sds%data_type_uncal ) 
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write calibration record failed on dataset:", &
                           sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_createSDS", zero )
           RETURN
        ENDIF

        count = LEN_TRIM( sds%long_name )
        IF( count > 0 ) THEN
           ierr = sfscatt( sds%id, "long_name", DFNT_CHAR8, &
                           count, TRIM(sds%long_name) ) 
           IF( ierr == FAIL ) THEN
              status = OZT_E_FAILURE
              WRITE( msg,* ) "write long name failed on dataset:", sds%name
              ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_createSDS", zero )
              RETURN
           ENDIF
        ENDIF

        ierr = sfsdtstr(sds%id, sds%label, sds%unit, sds%format, sds%coordsys) 
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write predefined attributes failed on dataset:", &
                          sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_createSDS", zero )
           RETURN
        ENDIF

        SELECT CASE( sds%data_type )
          CASE ( DFNT_CHAR8, DFNT_UCHAR8 )
            RETURN
          CASE ( DFNT_INT8, DFNT_UINT8 )
            fill_int8 = sds%fill_value
            max_int8  = sds%max  
            min_int8  = sds%min  
            ierr = sfsfill( sds%id, fill_int8 )
            ierr = sfsrange(sds%id, min_int8, max_int8)
            RETURN
          CASE ( DFNT_INT16, DFNT_UINT16 )
            fill_int16 = sds%fill_value
            max_int16  = sds%max  
            min_int16  = sds%min  
            ierr = sfsfill( sds%id, fill_int16 )
            ierr = sfsrange(sds%id, min_int16, max_int16)
            RETURN
          CASE ( DFNT_INT32, DFNT_UINT32 )
            fill_int32 = sds%fill_value
            max_int32  = sds%max  
            min_int32  = sds%min  
            ierr = sfsfill( sds%id, fill_int32 )
            ierr = sfsrange(sds%id, min_int32, max_int32)
            RETURN
          CASE ( DFNT_INT64, DFNT_UINT64 )
            RETURN
          CASE ( DFNT_FLOAT32 )
            fill_f32 = sds%fill_value
            max_f32  = sds%max  
            min_f32  = sds%min  
            ierr = sfsfill( sds%id, fill_f32 )
            ierr = sfsrange(sds%id, min_f32, max_f32)
            RETURN
          CASE DEFAULT
            ierr = sfsfill( sds%id, sds%fill_value )
            ierr = sfsrange(sds%id, sds%min, sds%max )
            RETURN
        END SELECT 
      END FUNCTION UTIL_lh4_createSDS

      FUNCTION UTIL_lh4_put1DIK2( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put1DIK2", zero )
           RETURN
        ENDIF

      END FUNCTION UTIL_lh4_put1DIK2

      FUNCTION UTIL_lh4_put2DIK2( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put2DIK2", zero )
           RETURN
        ENDIF

      END FUNCTION UTIL_lh4_put2DIK2

      FUNCTION UTIL_lh4_put3DIK2( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put3DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put3DIK2

      FUNCTION UTIL_lh4_put4DIK2( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put4DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put4DIK2

      FUNCTION UTIL_lh4_put5DIK2( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put5DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put5DIK2

      FUNCTION UTIL_lh4_put6DIK2( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put6DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put6DIK2

      FUNCTION UTIL_lh4_put1DRK4( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        REAL (KIND=4), DIMENSION(:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put1DRK4", zero )
           RETURN
        ENDIF

      END FUNCTION UTIL_lh4_put1DRK4

      FUNCTION UTIL_lh4_put2DRK4( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        REAL (KIND=4), DIMENSION(:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put2DRK4", zero )
           RETURN
        ENDIF

      END FUNCTION UTIL_lh4_put2DRK4

      FUNCTION UTIL_lh4_put3DRK4( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        REAL (KIND=4), DIMENSION(:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put3DRK4", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put3DRK4

      FUNCTION UTIL_lh4_put4DRK4( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        REAL (KIND=4), DIMENSION(:,:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put4DRK4", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put4DRK4

      FUNCTION UTIL_lh4_put5DRK4( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        REAL (KIND=4), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put5DRK4", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put5DRK4

      FUNCTION UTIL_lh4_put6DRK4( sds, data_in, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        REAL (KIND=4), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfwdata( sds%id, start_l, stride_l, count_l, data_in )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_put6DRK4", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_put6DRK4

      FUNCTION UTIL_lh4_get1DIK2( sds, data_out, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:), INTENT(OUT) :: data_out 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfrdata( sds%id, start_l, stride_l, count_l, data_out )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_get1DIK2", zero )
           RETURN
        ENDIF

      END FUNCTION UTIL_lh4_get1DIK2

      FUNCTION UTIL_lh4_get2DIK2( sds, data_out, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:), INTENT(OUT) :: data_out 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfrdata( sds%id, start_l, stride_l, count_l, data_out )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_get2DIK2", zero )
           RETURN
        ENDIF

      END FUNCTION UTIL_lh4_get2DIK2

      FUNCTION UTIL_lh4_get3DIK2( sds, data_out, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:), INTENT(OUT) :: data_out 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfrdata( sds%id, start_l, stride_l, count_l, data_out )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_get3DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_get3DIK2

      FUNCTION UTIL_lh4_get4DIK2( sds, data_out, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfrdata( sds%id, start_l, stride_l, count_l, data_out )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_get4DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_get4DIK2

      FUNCTION UTIL_lh4_get5DIK2( sds, data_out, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfrdata( sds%id, start_l, stride_l, count_l, data_out )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_get5DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_get5DIK2

      FUNCTION UTIL_lh4_get6DIK2( sds, data_out, start, count ) RESULT (status)
        TYPE (SDS_T), INTENT( IN ) :: sds
        INTEGER (KIND=2), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out 
        INTEGER (KIND=4), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
        INTEGER (KIND=4), DIMENSION(MAXRANK) :: start_l, count_l, stride_l 
        CHARACTER (LEN =255) :: msg
        INTEGER :: status 
        INTEGER :: ierr 
        stride_l = 1
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           start_l(1:sds%rank) = start(1:sds%rank) 
           count_l(1:sds%rank) = count(1:sds%rank) 
        ELSE
           start_l(1:sds%rank) = 0
           count_l(1:sds%rank) = sds%dims(1:sds%rank) 
        ENDIF

        ierr = sfrdata( sds%id, start_l, stride_l, count_l, data_out )
               
        IF( ierr == FAIL ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "write failed on dataset:", sds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh4_get6DIK2", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh4_get6DIK2
END MODULE UTIL_lh4_class

