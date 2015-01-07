MODULE UTIL_lh5_class
    USE HDF5
    USE OMI_SMF_class
    IMPLICIT NONE
    INTEGER (KIND=4), PARAMETER, PRIVATE :: MAXRANK = 7, zero = 0
   
    TYPE, PUBLIC :: DSh5_T
      CHARACTER (LEN = 80) :: name
      INTEGER (HID_T)  :: dataset_id 
      INTEGER (HID_T)  :: datatype  !datatype determines class, order, size, 
                                    !and vise versa 
      INTEGER          :: class
      INTEGER (SIZE_T) :: size
      INTEGER (HID_T)  :: dataspace !dataspace determines rank, dims 
      INTEGER          :: rank
      INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
    END TYPE DSh5_T

    INTERFACE UTIL_lh5_get
      MODULE PROCEDURE UTIL_lh5_get1DIK4
      MODULE PROCEDURE UTIL_lh5_get2DIK4
      MODULE PROCEDURE UTIL_lh5_get3DIK4
      MODULE PROCEDURE UTIL_lh5_get4DIK4
      MODULE PROCEDURE UTIL_lh5_get5DIK4
      MODULE PROCEDURE UTIL_lh5_get6DIK4
      MODULE PROCEDURE UTIL_lh5_get1DRK4
      MODULE PROCEDURE UTIL_lh5_get2DRK4
      MODULE PROCEDURE UTIL_lh5_get3DRK4
      MODULE PROCEDURE UTIL_lh5_get4DRK4
      MODULE PROCEDURE UTIL_lh5_get5DRK4
      MODULE PROCEDURE UTIL_lh5_get6DRK4
      MODULE PROCEDURE UTIL_lh5_get1DRK8
      MODULE PROCEDURE UTIL_lh5_get2DRK8
      MODULE PROCEDURE UTIL_lh5_get3DRK8
      MODULE PROCEDURE UTIL_lh5_get4DRK8
      MODULE PROCEDURE UTIL_lh5_get5DRK8
      MODULE PROCEDURE UTIL_lh5_get6DRK8
    END INTERFACE 
    
    INTERFACE UTIL_lh5_put
      MODULE PROCEDURE UTIL_lh5_put1DIK4
      MODULE PROCEDURE UTIL_lh5_put2DIK4
      MODULE PROCEDURE UTIL_lh5_put3DIK4
      MODULE PROCEDURE UTIL_lh5_put4DIK4
      MODULE PROCEDURE UTIL_lh5_put5DIK4
      MODULE PROCEDURE UTIL_lh5_put6DIK4
      MODULE PROCEDURE UTIL_lh5_put1DRK4
      MODULE PROCEDURE UTIL_lh5_put2DRK4
      MODULE PROCEDURE UTIL_lh5_put3DRK4
      MODULE PROCEDURE UTIL_lh5_put4DRK4
      MODULE PROCEDURE UTIL_lh5_put5DRK4
      MODULE PROCEDURE UTIL_lh5_put6DRK4
      MODULE PROCEDURE UTIL_lh5_put1DRK8
      MODULE PROCEDURE UTIL_lh5_put2DRK8
      MODULE PROCEDURE UTIL_lh5_put3DRK8
      MODULE PROCEDURE UTIL_lh5_put4DRK8
      MODULE PROCEDURE UTIL_lh5_put5DRK8
      MODULE PROCEDURE UTIL_lh5_put6DRK8
    END INTERFACE 

    PUBLIC :: UTIL_lh5_inquire
    PUBLIC :: UTIL_lh5_selectDS
!    PUBLIC :: UTIL_lh5_createDS
    PUBLIC :: UTIL_lh5_disposeDS

    CONTAINS
      FUNCTION UTIL_lh5_inquire( file_id, grpname, ds ) RESULT (status)
        INTEGER (HID_T), INTENT( IN )  :: file_id 
        CHARACTER(LEN=*), INTENT(IN) :: grpname        ! Name of the group 
        TYPE (DSh5_T), INTENT( IN ) :: ds 
        INTEGER :: status 
        INTEGER :: nmembers            ! Number of members in the group
        INTEGER :: hdferr              ! Error code 
        INTEGER :: idx                 ! Index of member object 
        INTEGER :: obj_type            ! Object type : 
        CHARACTER(LEN=256) :: obj_name ! Name of the object 
        status = OZT_E_FAILURE
        CALL h5gn_members_f(file_id, grpname, nmembers, hdferr)           
        IF( hdferr == -1 ) RETURN

        DO idx = 0, nmembers-1
          CALL h5gget_obj_info_idx_f( file_id, grpname, idx, & 
                                      obj_name, obj_type, hdferr)           
          IF( hdferr == -1 ) RETURN
        
          IF( INDEX(obj_name, ds%name) > 0 ) THEN
             status = OZT_S_SUCCESS
             RETURN   
          ENDIF
        ENDDO 
        RETURN   
      END FUNCTION UTIL_lh5_inquire


      FUNCTION UTIL_lh5_selectDS( file_id, ds ) RESULT (status)
        INTEGER (HID_T), INTENT( IN )  :: file_id 
        TYPE (DSh5_T), INTENT( INOUT ) :: ds 
        INTEGER (KIND = 4) :: status
        INTEGER          :: class
        INTEGER (SIZE_T) :: size
        INTEGER          :: irank, rank;
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims, maxdims
        INTEGER (HID_T)  :: datatype 
        INTEGER          :: error, ierr
        CHARACTER (LEN =255) :: msg

        status = UTIL_lh5_inquire( file_id, "/", ds ) 
        IF( status /= OZT_S_SUCCESS ) RETURN

        CALL h5dopen_f(file_id, ds%name, ds%dataset_id, error ) 
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dopen_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ENDIF

        CALL h5dget_type_f( ds%dataset_id, datatype, error )
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dget_type_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ENDIF

        !  make sure class and size match
        CALL h5tget_class_f(datatype, class, error )
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5tget_class_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ELSE IF( ds%class .NE. class ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "class_in=", ds%class, "not equal to", &
                           class, " in ", ds%name
           RETURN
        ENDIF

        CALL h5tget_size_f( datatype, size, error )
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5tget_size_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ELSE IF( ds%size .NE. size ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "size_in=", ds%size, "not equal to", &
                           size, " in ", ds%name
           RETURN
        ENDIF
        
        CALL h5dget_space_f(ds%dataset_id, ds%dataspace, error )
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5tget_space_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ENDIF
        
        ! get the rank and dim_sizes 
        CALL h5sget_simple_extent_ndims_f( ds%dataspace, rank, error )
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5sget_simple_extent_ndims_f failed on dataset: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ELSE IF( rank > MAXRANK ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "rank=",rank," exceeds maxrank=",MAXRANK, &
                          "for ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ELSE IF( rank .NE. ds%rank ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "rank_in=",ds%rank," not equal to rank=",rank, &
                          "for ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ENDIF

        CALL h5sget_simple_extent_dims_f(ds%dataspace, dims, maxdims, error ) 
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5sget_simple_extent_dims_f failed on dataset: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_selectDS", zero )
           RETURN
        ENDIF
        ds%dims(1:rank) = dims(1:rank)
      
        RETURN
      END FUNCTION UTIL_lh5_selectDS

!      FUNCTION UTIL_lh5_createDS( file_id, ds ) RESULT (status)
!        INTEGER (HID_T), INTENT( IN )  :: file_id 
!        TYPE (DSh5_T), INTENT( INOUT ) :: ds 
!        INTEGER(HID_T) :: dataspace     ! Dataspace identifier
!        INTEGER (KIND = 4) :: status
!        INTEGER          :: class
!        INTEGER (SIZE_T) :: size
!        INTEGER          :: irank, rank;
!        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims, maxdims
!        INTEGER (HID_T)  :: datatype 
!        INTEGER          :: error, ierr
!        CHARACTER (LEN =255) :: msg
!
!        status = OZT_S_SUCCESS
!
!        IF( ds%datatype < 0 ) THEN
!           status = OZT_E_FAILURE
!           WRITE( msg,* ) "valid datatype must be set before creation: ",&
!                           ds%name
!           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!                                  "UTIL_lh5_createDS", zero )
!           RETURN
!        ENDIF
!
!        !
!        ! Create the data space for the  dataset.
!        !
!        CALL h5screate_simple_f(ds%rank, ds%dims, ds%dataspace, error )
!        IF( error < zero ) THEN
!           status = OZT_E_FAILURE
!           WRITE( msg,* ) "h5screate_simple_f failed on dataset: ", ds%name
!           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!                                 "UTIL_lh5_createDS", zero )
!           RETURN
!        ENDIF
!
!        !
!        ! Create the dataset with default properties.
!        !
!        CALL h5dcreate_f(file_id, ds%name, ds%datatype, ds%dataspace, &
!                         ds%dataset_id, error)
!           
!        IF( error < zero ) THEN
!           status = OZT_E_FAILURE
!           WRITE( msg,* ) "h5dcreate_f dataspace failed on dataset: ", ds%name
!           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
!                                 "UTIL_lh5_createDS", zero )
!           RETURN
!        ENDIF
!
!      END FUNCTION UTIL_lh5_createDS

      FUNCTION UTIL_lh5_get1DIK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds 
        INTEGER (KIND = 4), DIMENSION(:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 1
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
        
        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get1DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank) 
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get1DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get1DIK4

      FUNCTION UTIL_lh5_get2DIK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 2
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get2DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get2DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get2DIK4

      FUNCTION UTIL_lh5_get3DIK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 3
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get3DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get3DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get3DIK4

      FUNCTION UTIL_lh5_get4DIK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 4
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get4DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get4DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get4DIK4

      FUNCTION UTIL_lh5_get5DIK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 5
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get5DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get5DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get5DIK4

      FUNCTION UTIL_lh5_get6DIK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 6
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get6DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get6DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get6DIK4

      FUNCTION UTIL_lh5_put1DIK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds 
        INTEGER (KIND = 4), DIMENSION(:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 1
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
        
        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put1DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank) 
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put1DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put1DIK4

      FUNCTION UTIL_lh5_put2DIK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 2
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put2DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put2DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put2DIK4

      FUNCTION UTIL_lh5_put3DIK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 3
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put3DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put3DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put3DIK4

      FUNCTION UTIL_lh5_put4DIK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 4
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put4DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put4DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put4DIK4

      FUNCTION UTIL_lh5_put5DIK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 5
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put5DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put5DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put5DIK4

      FUNCTION UTIL_lh5_put6DIK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        INTEGER (KIND = 4), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 6
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put6DIK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put6DIK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put6DIK4


      FUNCTION UTIL_lh5_get1DRK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds 
        REAL (KIND = 4), DIMENSION(:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 1
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
        
        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get1DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank) 
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get1DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get1DRK4

      FUNCTION UTIL_lh5_get2DRK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 2
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get2DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get2DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get2DRK4

      FUNCTION UTIL_lh5_get3DRK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 3
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get3DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get3DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get3DRK4

      FUNCTION UTIL_lh5_get4DRK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 4
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get4DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get4DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get4DRK4

      FUNCTION UTIL_lh5_get5DRK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 5
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get5DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get5DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get5DRK4

      FUNCTION UTIL_lh5_get6DRK4( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 6
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get6DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get6DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get6DRK4

      FUNCTION UTIL_lh5_put1DRK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds 
        REAL (KIND = 4), DIMENSION(:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 1
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
        
        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put1DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank) 
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put1DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put1DRK4

      FUNCTION UTIL_lh5_put2DRK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 2
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put2DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put2DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put2DRK4

      FUNCTION UTIL_lh5_put3DRK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 3
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put3DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put3DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put3DRK4

      FUNCTION UTIL_lh5_put4DRK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 4
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put4DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put4DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put4DRK4

      FUNCTION UTIL_lh5_put5DRK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 5
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put5DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put5DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put5DRK4

      FUNCTION UTIL_lh5_put6DRK4( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 4), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 6
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put6DRK4", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put6DRK4", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put6DRK4


      FUNCTION UTIL_lh5_get1DRK8( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds 
        REAL (KIND = 8), DIMENSION(:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 1
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
        
        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get1DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank) 
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get1DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get1DRK8

      FUNCTION UTIL_lh5_get2DRK8( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 2
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get2DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get2DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get2DRK8

      FUNCTION UTIL_lh5_get3DRK8( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 3
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get3DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get3DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get3DRK8

      FUNCTION UTIL_lh5_get4DRK8( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 4
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get4DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get4DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get4DRK8

      FUNCTION UTIL_lh5_get5DRK8( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 5
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get5DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get5DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get5DRK8

      FUNCTION UTIL_lh5_get6DRK8( ds, data_out, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 6
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before read: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_get6DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_out )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dread_f( ds%dataset_id, ds%datatype, data_out, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_get6DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_get6DRK8

      FUNCTION UTIL_lh5_put1DRK8( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds 
        REAL (KIND = 8), DIMENSION(:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 1
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
        
        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put1DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank) 
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF

        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put1DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put1DRK8

      FUNCTION UTIL_lh5_put2DRK8( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 2
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put2DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put2DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put2DRK8

      FUNCTION UTIL_lh5_put3DRK8( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 3
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put3DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put3DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put3DRK8

      FUNCTION UTIL_lh5_put4DRK8( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 4
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put4DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )
        
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put4DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put4DRK8

      FUNCTION UTIL_lh5_put5DRK8( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 5
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put5DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put5DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put5DRK8

      FUNCTION UTIL_lh5_put6DRK8( ds, data_in, start, count ) RESULT (status)
        TYPE (DSh5_T), INTENT( IN ) :: ds
        REAL (KIND = 8), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
        INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
        INTEGER          :: error, ierr
        INTEGER (HID_T)  :: dataspace 
        INTEGER (KIND = 4) :: status
        CHARACTER (LEN =255) :: msg
        INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
        INTEGER :: rankm = 6
        INTEGER(HID_T) :: memspace 
        INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

        IF( ds%dataset_id < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "dataset must be selected before write: ",&
                           ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                  "UTIL_lh5_put6DRK8", zero )
           RETURN
        ENDIF

        dims(1:ds%rank) = ds%dims(1:ds%rank)
        IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
           dataspace = ds%dataspace
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                      start, count, error)
           dimsm = SHAPE( data_in )
           CALL h5screate_simple_f( rankm, dimsm, memspace, error)
        ELSE
           dataspace = H5S_ALL_F
           memspace = H5S_ALL_F
        ENDIF
        
        CALL h5dwrite_f( ds%dataset_id, ds%datatype, data_in, dims, error, &
                        memspace, dataspace )

        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_put6DRK8", zero )
           RETURN
        ENDIF
        IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)
      END FUNCTION UTIL_lh5_put6DRK8

      FUNCTION UTIL_lh5_disposeDS( ds ) RESULT (status)
        TYPE (DSh5_T), INTENT( INOUT ) :: ds 
        INTEGER (KIND = 4) :: status
        INTEGER          :: error, ierr
        CHARACTER (LEN =255) :: msg

        status = OZT_S_SUCCESS
        !
        ! Close the dataspace for the dataset.
        !
        CALL h5sclose_f(ds%dataspace, error)
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5sclose_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_dispoaeDS", zero )
           RETURN
        ENDIF

        !
        ! Close the dataset.
        !
        CALL h5dclose_f(ds%dataset_id, error)
        IF( error < zero ) THEN
           status = OZT_E_FAILURE
           WRITE( msg,* ) "h5dclose_f failed on dataset: ", ds%name
           ierr = OMI_SMF_setmsg( OZT_E_HDF, TRIM(msg), &
                                 "UTIL_lh5_dispoaeDS", zero )
           RETURN
        ENDIF
      END FUNCTION UTIL_lh5_disposeDS

END MODULE UTIL_lh5_class

