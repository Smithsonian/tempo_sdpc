MODULE netcdf_module

  USE parameters_module
  USE error_module,     ONLY : ErrorType, RaiseFatalError, CheckError
  
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  TYPE NCDimType
      CHARACTER(LEN=maxChar)             :: NAME
      INTEGER                            :: DIMSIZE
      INTEGER                            :: NCID
      CHARACTER(LEN=maxChar)             :: DESCRIPTION
  ENDTYPE NCDimType
  
  ! For error checking
  CHARACTER(LEN=*), PARAMETER :: ModuleName = 'netcdf_module'
  PRIVATE :: ModuleName

  CONTAINS
  
  SUBROUTINE define_ncdf_dim( ncid, ncdim, Error )
    
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    INTEGER,           INTENT(IN   ) :: ncid
    TYPE( NCDimType ), INTENT(INOUT) :: ncdim
    TYPE(ErrorType),   INTENT(INOUT) :: Error ! Error variable
    
    ! ---------------
    ! Local variables
    ! ---------------
    
    INTEGER :: rcode
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'define_ncdf_dim'

    !==============================================================================
    ! define_ncdf_dim starts here
    !==============================================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
    
    ! Create the dimensions of the dataset:
    rcode = NF_DEF_DIM(ncid, &
                       TRIM(ADJUSTL(ncdim%NAME)), &
                       ncdim%DIMSIZE,             &
                       ncdim%NCID                 )
    
    ! Handle error
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    
  END SUBROUTINE define_ncdf_dim
  
  SUBROUTINE close_ncdf_file( ncid, error )
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    INTEGER,           INTENT(IN   ) :: ncid
    TYPE(ErrorType),   INTENT(INOUT) :: Error ! Error variable
    
    ! ---------------
    ! Local variables
    ! ---------------
    
    INTEGER :: rcode
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'close_ncdf_file'

    !==============================================================================
    ! close_ncdf_file starts here
    !==============================================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Create the dimensions of the dataset:
    rcode = NF_CLOSE(ncid)

    ! Handle error
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)

  END SUBROUTINE close_ncdf_file
  
  LOGICAL FUNCTION ncdf_var_exists(ncid,varname)

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    ! --------------------
    ! Subroutine arguments
    ! --------------------
    INTEGER,           INTENT(IN) :: ncid
    CHARACTER(LEN=*),  INTENT(IN) :: varname
    
    ! ---------------
    ! Local variables
    ! ---------------
    
    INTEGER :: rcode, vid
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'ncdf_var_exists'

    !==============================================================================
    ! ncdf_var_exists starts here
    !==============================================================================
    
    ! Initialize return
    ncdf_var_exists = .FALSE.

    ! Attempt to attach variable
    rcode = nf_inq_varid( ncid, TRIM(ADJUSTL(varname)), vid )

    ! Check for error
    IF(rcode .EQ. NF_NOERR) ncdf_var_exists = .TRUE.

  END FUNCTION ncdf_var_exists

  SUBROUTINE netcdf_handle_error(location,status)
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER, INTENT(IN) :: status
    CHARACTER(*)        :: location

    ! ---------------
    ! Local variables
    !----------------
    CHARACTER(NF_MAX_NAME) :: message
    
    ! Code starts here
    IF (status .NE. NF_NOERR) THEN
       message = 'Error in '//TRIM(location)//': '//NF_STRERROR(status)
       print*,message
       STOP 
    ENDIF

  END SUBROUTINE netcdf_handle_error
  
  SUBROUTINE CheckNetCDFErrorStatus(Error,status,CallingModuleName,CallingSubroutineName,OperationDescription)
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

    ! ---------------
    ! Input variables
    ! ---------------
    TYPE(ErrorType),          INTENT(INOUT) :: Error
    INTEGER,                  INTENT(IN)    :: status
    CHARACTER(LEN=*),         INTENT(IN)    :: CallingModuleName
    CHARACTER(LEN=*),         INTENT(IN)    :: CallingSubroutineName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: OperationDescription


    ! ---------------
    ! Local variables
    !----------------
    CHARACTER(Maxchar) :: Message

    ! =====================================================
    !  CheckNetCDFErrorStatus starts here
    ! =====================================================

    IF (status .NE. NF_NOERR) THEN
      Message = TRIM(ADJUSTL(NF_STRERROR(status)))
      IF(PRESENT(OperationDescription)) THEN
        Message = TRIM(ADJUSTL(Message)) // '(' // &
                  TRIM(ADJUSTL(OperationDescription)) // ')'
      ENDIF
      CALL RaiseFatalError(Error, ErrorCode_FileIO,                      &
                           CallingModuleName, CallingSubroutineName,     &
                           Message_in=Message)

    ENDIF
  END SUBROUTINE CheckNetCDFErrorStatus

  SUBROUTINE match_names_in_dimlist( nc_dimlist, nc_ndimlist, DIMNAMES, ND, DIM_ID, DIMS, Error, Varname )
    
    
    ! ---------------------------
    ! Subroutine arguments
    ! ---------------------------
    INTEGER,                                  INTENT(IN)    :: nc_ndimlist
    TYPE(NCDimType),  DIMENSION(nc_ndimlist), INTENT(IN)    :: nc_dimlist
    INTEGER,                                  INTENT(IN)    :: ND
    CHARACTER(LEN=*), DIMENSION(ND),          INTENT(IN)    :: DIMNAMES
    INTEGER,          DIMENSION(ND),          INTENT(OUT)   :: DIM_ID, DIMS
    TYPE(ErrorType),                          INTENT(INOUT) :: Error
    CHARACTER(LEN=*), OPTIONAL,               INTENT(IN)    :: Varname ! For debug logger

    ! ---------------------------
    ! Local variables
    ! ---------------------------
    INTEGER :: n,d
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'match_names_in_dimlist'
    CHARACTER(LEN=maxChar)      :: ErrorMessage

    ! =====================================================
    ! match_names_in_dimlist starts here
    ! =====================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
    
    ! Initialize dim_id
    DIM_ID(:) = -1
    
    ! Match dimension names
    DO n=1,nc_ndimlist
      DO d=1,ND
        IF( TRIM(ADJUSTL(DIMNAMES(d))) .EQ. TRIM(ADJUSTL(nc_dimlist(n)%NAME)) ) THEN
          DIM_ID(d) = nc_dimlist(n)%NCID
          DIMS(d)   = nc_dimlist(n)%DIMSIZE
        ENDIF
      ENDDO
    ENDDO
    
    ! Error check
    DO d=1,ND
      IF( DIM_ID(d) < 0 ) THEN
        ErrorMessage = 'Could match dimension for ' // TRIM(ADJUSTL(DIMNAMES(d))) 
        IF(PRESENT(Varname)) ErrorMessage = TRIM(ADJUSTL(ErrorMessage))  // &
                             ' for variable' // TRIM(ADJUSTL(Varname))
        CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,   &
                              Message_in= TRIM(ADJUSTL(ErrorMessage)),               &
                              Action_in='Ensure Dimension is defined in NetCDF file' )
      ENDIF
    ENDDO
    
  END SUBROUTINE match_names_in_dimlist
  
  SUBROUTINE create_ncvar( NCID, nc_dimlist, nc_ndimlist, VARNAME, DIMNAME, NDIM, NCTYPE, DO_XY, VAR_ID, ERROR  )
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
      
    ! ---------------------------
    ! Subroutine arguments
    ! ---------------------------
    INTEGER,                                    INTENT(IN)  :: nc_ndimlist
    TYPE(NCDimType),  DIMENSION(nc_ndimlist),   INTENT(IN)  :: nc_dimlist
    INTEGER,                                    INTENT(IN)  :: NCID, NDIM, NCTYPE
    CHARACTER(LEN=*),                           INTENT(IN)  :: VARNAME
    CHARACTER(LEN=*),          DIMENSION(NDIM), INTENT(IN)  :: DIMNAME
    LOGICAL,                                    INTENT(IN)  :: DO_XY
    INTEGER,                                    INTENT(OUT) :: VAR_ID
    TYPE(ErrorType),                          INTENT(INOUT) :: Error

    ! ---------------------------
    ! Local variables
    ! ---------------------------
    INTEGER,                                DIMENSION(NDIM) :: DIM_ID, DIMS
    INTEGER,                   ALLOCATABLE, DIMENSION(:)    :: DIM_ID_XY, DIMS_XY, CHUNK
    CHARACTER(LEN=maxChar), ALLOCATABLE, DIMENSION(:)       :: DIMNAME_XY
    INTEGER                                                 :: NDIM_XY
    INTEGER                                                 :: d, rcode
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'create_ncvar'

    ! =====================================================
    ! create_ncvar starts here
    ! =====================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
      
    IF( DO_XY ) THEN
      
      ! Append X-Y Dimensions
      NDIM_XY = NDIM+2
        
      ! Allocate arrays
      ALLOCATE(  DIM_ID_XY(NDIM_XY) )
      ALLOCATE(    DIMS_XY(NDIM_XY) )
      ALLOCATE(      CHUNK(NDIM_XY) )
      ALLOCATE( DIMNAME_XY(NDIM_XY) )
        
      ! Set dimension names
      DIMNAME_XY(1) = 'imx'; DIMNAME_XY(2) = 'jmx'
      DO d=1,NDIM
        DIMNAME_XY(d+2) = DIMNAME(d)
      ENDDO
      
      ! Match the dimension names
      CALL match_names_in_dimlist( nc_dimlist, nc_ndimlist, DIMNAME_XY, NDIM_XY, DIM_ID_XY, DIMS_XY, Error, VARNAME )
        
      ! Create variable
      rcode = NF_DEF_VAR( NCID, TRIM(ADJUSTL(varname)), NCTYPE, NDIM_XY, DIM_ID_XY, VAR_ID )
        
      ! Check error 
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
      
      ! Set up chunking array
      CHUNK(1:2)       = 1
      CHUNK(3:NDIM_XY) = DIMS_XY(3:NDIM_XY)
      
      ! Define chunking
      rcode = NF_DEF_VAR_CHUNKING(NCID, VAR_ID, NF_CHUNKED, CHUNK)
      
      ! Handle error
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
        
      ! Deallocate arrays
      DEALLOCATE(  DIM_ID_XY )
      DEALLOCATE(    DIMS_XY )
      DEALLOCATE(      CHUNK )
      DEALLOCATE( DIMNAME_XY )
        
    ELSE
        
      ! Match the dimension names
      CALL match_names_in_dimlist( nc_dimlist, nc_ndimlist, DIMNAME, NDIM, DIM_ID, DIMS, Error, VARNAME )
      
      ! Create variable
      rcode = NF_DEF_VAR( NCID, TRIM(ADJUSTL(varname)), NCTYPE, NDIM, DIM_ID, VAR_ID )
      
      ! Handle error
      CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
        
    ENDIF
    
  END SUBROUTINE create_ncvar
  
  SUBROUTINE write_ncvar( NCID, NC_XID, NC_YID,  VARNAME, START, SUBDIM, NDIM, VDIM, NCTYPE, DO_XY, ERROR, DATA_r8, DATA_i4 )
    
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
      
    ! ---------------------------
    ! Subroutine arguments
    ! ---------------------------
      
    INTEGER,                                    INTENT(IN)    :: NCID, NC_XID, NC_YID, NCTYPE, NDIM, VDIM
    CHARACTER(LEN=*),                           INTENT(IN)    :: VARNAME
    INTEGER,                   DIMENSION(NDIM), INTENT(IN)    :: START, SUBDIM
    LOGICAL,                                    INTENT(IN)    :: DO_XY
    TYPE(ErrorType),                            INTENT(INOUT) :: Error
    REAL(KIND=8),    OPTIONAL, DIMENSION(VDIM), INTENT(IN)    :: DATA_r8
    INTEGER(KIND=4), OPTIONAL, DIMENSION(VDIM), INTENT(IN)    :: DATA_i4
      
    ! ---------------------------
    ! Local variables
    ! ---------------------------
    INTEGER                            :: rcode, var_prec, NDIM_OUT, VAR_ID
    INTEGER, ALLOCATABLE, DIMENSION(:) :: START_OUT, SUBDIM_OUT
    INTEGER, DIMENSION(7)              :: ARDIM ! For unwrapping the vector
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'write_ncvar'
    CHARACTER(LEN=maxChar)      :: nctype_str

    ! =====================================================
    ! create_ncvar starts here
    ! =====================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
    
    ! Figure out variable precision
    IF(     NCTYPE .EQ. nf_byte  .OR. NCTYPE .EQ. nf_char      ) THEN
      var_prec = 1
    ELSEIF( NCTYPE .EQ. nf_short .OR. NCTYPE .EQ. nf_ushort    ) THEN
      var_prec = 2
    ELSEIF( NCTYPE .EQ. nf_int   .OR. NCTYPE .EQ. nf_uint   .OR. &
            NCTYPE .EQ. nf_float                               ) THEN
      var_prec = 4
    ELSEIF( NCTYPE .EQ. nf_int64 .OR. NCTYPE .EQ. nf_uint64 .OR. &
            NCTYPE .EQ. nf_double                              ) THEN
      var_prec = 8
    ELSE
      
      WRITE(nctype_str,'(I100)') nctype
      CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,   &
                            Message_in = 'Error writing field for variable ' // TRIM(ADJUSTL(VARNAME)) &
                                      // ' - Could not find precision for NC variable type:'           &
                                      // TRIM(ADJUSTL(nctype_str)),                                    &
                            Action_in='Ensure NCTYPE is one of nf_byte, nf_char, nf_short, nf_ushort,' &
                                      //'nf_int, nf_uint, nf_float, nf_int64, nf_uint64, or nf_double '&
                                      //'as defined in the netcdf.inc header file'                     )

    ENDIF
      
    ! Check if we are writing to XY
    IF( DO_XY ) THEN
      
      ! Dimensions of XY output
      NDIM_OUT = NDIM + 2
      
      ! Allocate start/count arrays
      ALLOCATE(  START_OUT(NDIM_OUT) )
      ALLOCATE( SUBDIM_OUT(NDIM_OUT) )
      
      ! Fill start/count arrays
      START_OUT(1)          = nc_xid; START_OUT(2) = nc_yid
      START_OUT(3:NDIM_OUT) = START(1:NDIM)
      SUBDIM_OUT(1:2) = 1
      SUBDIM_OUT(3:NDIM_OUT) = SUBDIM(1:NDIM)
      
    ELSE
      
      ! Kludge
      NDIM_OUT = NDIM
      
      ! Allocate start/count arrays
      ALLOCATE(  START_OUT(NDIM_OUT) )
      ALLOCATE( SUBDIM_OUT(NDIM_OUT) )
      
      ! Fill start/count arrays
      START_OUT = START ; SUBDIM_OUT = SUBDIM
      
    ENDIF
    
    ! Get variable index
    rcode = NF_INQ_VARID( NCID, TRIM(ADJUSTL(VARNAME)), VAR_ID)
      
    ! Handle error
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)
    
    ! Reshaping dimension
    ARDIM(:) = 1
    ARDIM( 1:NDIM ) = SUBDIM
    
    ! Write the variable
    IF( PRESENT(DATA_r8) ) THEN
      
      ! Write the variable
      SELECT CASE( NCTYPE )
        CASE( nf_float )
          rcode = NF_PUT_VARA_REAL(NCID, VAR_ID, START_OUT, SUBDIM_OUT,  &
                                   REAL( RESHAPE(DATA_r8,ARDIM), KIND=4) )
        CASE( nf_double )
          rcode = NF_PUT_VARA_DOUBLE(NCID, VAR_ID, START_OUT, SUBDIM_OUT, &
                                     RESHAPE(DATA_r8,ARDIM)               )
          
        CASE DEFAULT
          WRITE(nctype_str,'(I100)') nctype
          CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,   &
                                Message_in = 'Error writing field for variable ' // TRIM(ADJUSTL(VARNAME)) &
                                          // ' - No REAL output type for :' //  TRIM(ADJUSTL(nctype_str)), &
                                Action_in='Ensure NCTYPE is one of nf_float or nf_double '                 &
                                          //'as defined in the netcdf.inc header file'                     )
      END SELECT
        
      
     
      
    ELSEIF( PRESENT(DATA_i4) ) THEN
      
      ! Write the variable
      SELECT CASE( NCTYPE )
        CASE( nf_byte )
          rcode = NF_PUT_VARA_INT1( NCID, VAR_ID, START_OUT, SUBDIM_OUT,  &
                                    INT( RESHAPE(DATA_i4,ARDIM), KIND=1 ) )
        CASE( nf_int2 )
          rcode = NF_PUT_VARA_INT2( NCID, VAR_ID, START_OUT, SUBDIM_OUT,  &
                                    INT( RESHAPE(DATA_i4,ARDIM), KIND=2 ) )
        CASE( nf_int )
          rcode = NF_PUT_VARA_INT( NCID, VAR_ID, START_OUT, SUBDIM_OUT,  &
                                   RESHAPE(DATA_i4,ARDIM)                )
        CASE DEFAULT
          WRITE(nctype_str,'(I100)') nctype
          CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,   &
                                Message_in = 'Error writing field for variable ' // TRIM(ADJUSTL(VARNAME)) &
                                          // ' - No INT output type for :' //  TRIM(ADJUSTL(nctype_str)),  &
                                Action_in='Ensure NCTYPE is one of nf_byte, nf_int2, or nf_int '           &
                                          //'as defined in the netcdf.inc header file'                     )
      END SELECT

    ELSE
      
      CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,                         &
                            Message_in = 'Error: No data present for variable' // TRIM(ADJUSTL(VARNAME)),&
                            Action_in= 'One of Data_i4 or Data_r8 must be set!!!'                        )

    ENDIF
    
    ! Check for error writing field
    CALL CheckNetCDFErrorStatus(Error,rcode,ModuleName,SubroutineName)

    ! Deallocate arrays
    DEALLOCATE( START_OUT  )
    DEALLOCATE( SUBDIM_OUT )
    
  END SUBROUTINE write_ncvar
  
  SUBROUTINE nc_fld_1d( DIM_1D, LMX, PROFDATA, FLD_NAME, NCID, ACTION, DO_XY,             &
                        nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error, StartIn, LengthIn )
   
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    ! ---------------------------
    ! Subroutine arguments
    ! ---------------------------
    INTEGER,                INTENT(IN)     :: LMX
    CHARACTER(LEN=maxChar), INTENT(IN)    :: DIM_1D(1)
    CHARACTER(LEN=*),       INTENT(IN)    :: FLD_NAME
    REAL(KIND=8),           INTENT(IN)    :: PROFDATA(LMX)
    INTEGER,                INTENT(IN)    :: NCID
    INTEGER,                INTENT(IN)    :: ACTION
    LOGICAL,                INTENT(IN)    :: DO_XY
    INTEGER,                INTENT(IN)    :: nc_ndimlist
    TYPE(NCDimType),        INTENT(IN)    :: nc_dimlist(nc_ndimlist)
    INTEGER,                INTENT(IN)    :: NC_XID, NC_YID
    TYPE(ErrorType),        INTENT(INOUT) :: Error
    INTEGER, OPTIONAL,      INTENT(IN)    :: StartIn(1)
    INTEGER, OPTIONAL,      INTENT(IN)    :: LengthIn(1)

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: VAR_ID
    INTEGER :: Start(1)
    INTEGER :: Length(1)

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'nc_fld_1d'
    CHARACTER(LEN=maxChar)      :: ActStr

    ! =====================================================
    ! nc_fld_1d starts here
    ! =====================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN

    ! Optional Inputs - Write full dimension when not specified
    IF( PRESENT(StartIn) ) THEN
      Start(1) = StartIn(1)
    ELSE
      Start(1) = 1
    ENDIF
    IF( PRESENT(LengthIn) ) THEN
      Length(1) = LengthIn(1)
    ELSE
      Length(1) = lmx
    ENDIF
    
    ! Define field
    IF( ACTION == 1 ) THEN
      
      CALL create_ncvar( NCID, nc_dimlist, nc_ndimlist,                       &
                         TRIM(ADJUSTL(fld_name)),   DIM_1D,      1,           &
                         nf_float,                   DO_XY,   VAR_ID,  Error  )
      
      
    ! Write field
    ELSEIF( ACTION == 2 ) THEN
      
      ! Write 
      CALL write_ncvar( NCID, NC_XID, NC_YID,                                   &
                        TRIM(ADJUSTL(fld_name)),   Start, Length,     1,        &
                        lmx,        nf_float,  DO_XY,  Error,                   &
                        DATA_r8=profdata(1:Length(1))                           )
      
    ELSE
      WRITE(ActStr,'(I100)') ACTION
      CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,                 &
                            Message_in = ' No action defined for case ' // TRIM(ADJUSTL(ActStr)) &
                                         // ' and therefore could not update '                   &
                                         // TRIM(ADJUSTL(FLD_NAME)),                             &
                            Action_in= 'ACTION must be 1 (Create Field) or 2 (Write Field)'      )
      
    ENDIF
      
  END SUBROUTINE nc_fld_1d
  
  SUBROUTINE nc_fld_2d( DIM_2D, IMX, JMX, PROFDATA, FLD_NAME, NCID, ACTION, DO_XY,  &
                        nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error,             &
                        StartIn, LengthIn, NCType_in                                )
      
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    ! ---------------------------
    ! Subroutine arguments
    ! ---------------------------
    INTEGER,                                      INTENT(IN) :: IMX,JMX
    CHARACTER(LEN=maxChar), DIMENSION(2),         INTENT(IN) :: DIM_2D
    CHARACTER(LEN=*         ),                    INTENT(IN) :: FLD_NAME
    REAL(KIND=8),              DIMENSION(IMX,JMX),INTENT(IN) :: PROFDATA
    INTEGER, INTENT(IN) :: NCID
    INTEGER, INTENT(IN) :: ACTION
    LOGICAL, INTENT(IN) :: DO_XY
    TYPE(NCDimType),  DIMENSION(nc_ndimlist),   INTENT(IN)    :: nc_dimlist
    INTEGER,                                    INTENT(IN)    :: nc_ndimlist
    INTEGER,                                    INTENT(IN)    :: NC_XID, NC_YID
    TYPE(ErrorType),                            INTENT(INOUT) :: Error
    INTEGER, OPTIONAL,                          INTENT(IN)    :: StartIn(2)
    INTEGER, OPTIONAL,                          INTENT(IN)    :: LengthIn(2)
    INTEGER, OPTIONAL,                          INTENT(IN)    :: NCType_in
 
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                                 :: VAR_ID, npt
    INTEGER                                 :: Start(2)
    INTEGER                                 :: Length(2)
    INTEGER                                 :: NCType

    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'nc_fld_2d'
    CHARACTER(LEN=maxChar)      :: ActStr

    ! =====================================================
    ! nc_fld_2d starts here
    ! =====================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
    
    ! Default precision
    NCType = nf_float ; IF(PRESENT(NCType_in)) NCType = NCType_in

    ! Optional Inputs - Write full dimension when not specified
    IF( PRESENT(StartIn) ) THEN
      Start(:) = StartIn(:)
    ELSE
      Start(:) = 1
    ENDIF
    IF( PRESENT(LengthIn) ) THEN
      Length(:) = LengthIn(:)
    ELSE
      Length(:) = (/imx,jmx/)
    ENDIF


    ! Define field
    IF( ACTION == 1 ) THEN
      
      CALL create_ncvar( NCID, nc_dimlist, nc_ndimlist,       &
                         TRIM(ADJUSTL(fld_name)),  DIM_2D, 2, &
                         NCType ,   DO_XY,   VAR_ID,  ERROR   )
        
    ! Write field
    ELSEIF( ACTION == 2 ) THEN
      
      npt = imx*jmx
      
      ! Write 
      CALL write_ncvar( NCID, NC_XID, NC_YID,                      &
                        TRIM(ADJUSTL(fld_name)), Start, Length, 2, &
                        npt,    NCType ,  DO_XY,  ERROR,           &
                        DATA_r8= RESHAPE(PROFDATA,(/npt/))         )
      
    ELSE
      WRITE(ActStr,'(I100)') ACTION
      CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,                 &
                            Message_in = ' No action defined for case ' // TRIM(ADJUSTL(ActStr)) &
                                         // ' and therefore could not update '                   &
                                         // TRIM(ADJUSTL(FLD_NAME)),                             &
                            Action_in= 'ACTION must be 1 (Create Field) or 2 (Write Field)'      )
    ENDIF
    
    RETURN
      
  END SUBROUTINE nc_fld_2d
    
  SUBROUTINE nc_fld_3d( DIM_3D, IMX, JMX, LMX, PROFDATA, FLD_NAME, NCID, ACTION, DO_XY, &
                        nc_dimlist, nc_ndimlist, NC_XID, NC_YID, Error                  )
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    ! ---------------------------
    ! Subroutine arguments
    ! ---------------------------
    INTEGER,                                          INTENT(IN) :: IMX,JMX,LMX
    CHARACTER(LEN=maxChar), DIMENSION(3),             INTENT(IN) :: DIM_3D
    CHARACTER(LEN=*         ),                        INTENT(IN) :: FLD_NAME
    REAL(KIND=8),              DIMENSION(IMX,JMX,LMX),INTENT(IN) :: PROFDATA
    INTEGER, INTENT(IN) :: NCID
    INTEGER, INTENT(IN) :: ACTION
    LOGICAL, INTENT(IN) :: DO_XY
    TYPE(NCDimType),  DIMENSION(nc_ndimlist),   INTENT(IN)    :: nc_dimlist
    INTEGER,                                    INTENT(IN)    :: nc_ndimlist
    INTEGER,                                    INTENT(IN)    :: NC_XID, NC_YID
    TYPE(ErrorType),                            INTENT(INOUT) :: Error
    
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                                 :: VAR_ID, npt
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'nc_fld_3d'
    CHARACTER(LEN=maxChar)      :: ActStr

    ! =====================================================
    ! nc_fld_3d starts here
    ! =====================================================
      
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
      
    ! Define field
    IF( ACTION == 1 ) THEN
      
      CALL create_ncvar( NCID, nc_dimlist, nc_ndimlist,      &
                         TRIM(ADJUSTL(fld_name)), DIM_3D, 3, &
                         nf_float,  DO_XY,   VAR_ID,  ERROR  )
      
    ! Write field
    ELSEIF( ACTION == 2 ) THEN
        
      npt = imx*jmx*lmx
      
      ! Write field
      CALL write_ncvar( NCID, NC_XID, NC_YID,                                      &
                        TRIM(ADJUSTL(fld_name)), (/1,1,1/),    (/imx,jmx,lmx/), 3, &
                        npt,                nf_float,     DO_XY,  ERROR,           &
                        DATA_r8= RESHAPE(PROFDATA,(/npt/))                         )
        
    ELSE
      
      WRITE(ActStr,'(I100)') ACTION
      CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,                 &
                            Message_in = ' No action defined for case ' // TRIM(ADJUSTL(ActStr)) &
                                         // ' and therefore could not update '                   &
                                         // TRIM(ADJUSTL(FLD_NAME)),                             &
                            Action_in= 'ACTION must be 1 (Create Field) or 2 (Write Field)'      )
    ENDIF
      
    RETURN
    
  END SUBROUTINE nc_fld_3d
    
  SUBROUTINE nc_fld_4d( DIM_4D, IMX, JMX, LMX, NMX, PROFDATA, FLD_NAME, NCID, ACTION, DO_XY, &
                        nc_dimlist, nc_ndimlist, NC_XID, NC_YID,  Error                      )
    
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    
    ! ---------------------------
    ! Subroutine arguments
    ! ---------------------------
    INTEGER,                                              INTENT(IN) :: IMX,JMX,LMX,NMX
    CHARACTER(LEN=maxChar), DIMENSION(4),                 INTENT(IN) :: DIM_4D
    CHARACTER(LEN=*         ),                            INTENT(IN) :: FLD_NAME
    REAL(KIND=8),              DIMENSION(IMX,JMX,LMX,NMX),INTENT(IN) :: PROFDATA
    INTEGER, INTENT(IN) :: NCID
    INTEGER, INTENT(IN) :: ACTION
    LOGICAL, INTENT(IN) :: DO_XY
    TYPE(NCDimType),  DIMENSION(nc_ndimlist),   INTENT(IN)    :: nc_dimlist
    INTEGER,                                    INTENT(IN)    :: nc_ndimlist
    INTEGER,                                    INTENT(IN)    :: NC_XID, NC_YID
    TYPE(ErrorType),                            INTENT(INOUT) :: Error
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                                 :: VAR_ID, npt
    
    ! For error checking
    CHARACTER(LEN=*), PARAMETER :: SubroutineName = 'nc_fld_4d'
    CHARACTER(LEN=maxChar)      :: ActStr

    ! =====================================================
    ! nc_fld_4d starts here
    ! =====================================================
    
    ! Dont execute if error has been flagged
    IF(CheckError(Error)) RETURN
    
    ! Define field
    IF( ACTION == 1 ) THEN
      
      CALL create_ncvar( NCID, nc_dimlist, nc_ndimlist,      &
                         TRIM(ADJUSTL(fld_name)), DIM_4D, 4, &
                         nf_float, DO_XY,   VAR_ID,  Error   )
      
    ! Write field
    ELSEIF( ACTION == 2 ) THEN
      
      npt = imx*jmx*lmx*nmx
      
      ! Write field
      CALL write_ncvar( NCID, NC_XID, NC_YID,                                            &
                        TRIM(ADJUSTL(fld_name)), (/1,1,1,1/),    (/imx,jmx,lmx,nmx/), 4, &
                        npt,                nf_float,     DO_XY,  Error,                 &
                        DATA_r8= RESHAPE(PROFDATA,(/npt/))                               )
        
    ELSE
      
      WRITE(ActStr,'(I100)') ACTION
      CALL RaiseFatalError( Error, ErrorCode_FileIO, ModuleName, SubroutineName,                 &
                            Message_in = ' No action defined for case ' // TRIM(ADJUSTL(ActStr)) &
                                         // ' and therefore could not update '                   &
                                         // TRIM(ADJUSTL(FLD_NAME)),                             &
                            Action_in= 'ACTION must be 1 (Create Field) or 2 (Write Field)'      )

    ENDIF
    
    RETURN
      
  END SUBROUTINE nc_fld_4d
  


  ! ##################################################################################
  ! Some reading functions for convenience
  ! ##################################################################################
  
  ! Convenience routine for getting variable id and checking error in one step
  INTEGER  FUNCTION nc_inq_varid(ncid,varname,Error,            &
                                 ModuleName_in,SubroutineName_in)

    ! ------------------
    ! function arguments
    ! ------------------
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: varname
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: ModuleName_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: SubroutineName_in

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                :: rcode, vid
    CHARACTER(LEN=maxChar) :: SubroutineName,ModName

    ! =====================================================================
    ! nc_inq_varid starts here
    ! =====================================================================

    SubroutineName = 'nc_inq_varid'
    IF(PRESENT(SubroutineName_in)) SubroutineName = SubroutineName_in
    ModName = ModuleName
    IF(PRESENT(ModuleName_in)) ModName = ModuleName_in

    ! Attach
    rcode = nf_inq_varid(ncid,TRIM(ADJUSTL(varname)),nc_inq_varid)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)),   &
                                TRIM(ADJUSTL(SubroutineName)),        &
                               'nf_inq_varid:'//TRIM(ADJUSTL(varname)))

  END FUNCTION nc_inq_varid

  ! For reading the logic indices
  INTEGER(KIND=2) FUNCTION nc_read_i2_index(ncid,varname,Error,            &
                                            ModuleName_in,SubroutineName_in)

    ! ------------------
    ! function arguments
    ! ------------------
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: varname
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: ModuleName_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: SubroutineName_in

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                :: rcode, vid
    CHARACTER(LEN=maxChar) :: SubroutineName,ModName

    ! =====================================================================
    ! netcdf_read_i2_index starts here
    ! =====================================================================

    SubroutineName = 'nc_read_i2_index'
    IF(PRESENT(SubroutineName_in)) SubroutineName = SubroutineName_in
    ModName = ModuleName
    IF(PRESENT(ModuleName_in)) ModName = ModuleName_in

    ! 
    rcode = nf_inq_varid(ncid,TRIM(ADJUSTL(varname)),vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)),   &
                                TRIM(ADJUSTL(SubroutineName)),        &
                               'nf_inq_varid:'//TRIM(ADJUSTL(varname)))
    rcode = nf_get_vara_int2(ncid,vid,(/1/),(/1/),nc_read_i2_index)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)), &
                                TRIM(ADJUSTL(SubroutineName)),      &
                               'nf_get_vara:'//TRIM(ADJUSTL(varname)))

  END FUNCTION nc_read_i2_index

  ! For getting the i^th dimension of an array

  INTEGER FUNCTION nc_get_arr_dim(ncid,varname,dim_idx,Error,    &
                                  ModuleName_in,SubroutineName_in)

    ! ------------------
    ! function arguments
    ! ------------------
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: varname
    INTEGER,                    INTENT(IN)    :: dim_idx
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: ModuleName_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: SubroutineName_in

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                :: rcode, vid, ndim
    INTEGER, ALLOCATABLE   :: dimid(:)

    CHARACTER(LEN=maxChar) :: SubroutineName,ModName,tmpchar

    ! =====================================================================
    ! nc_get_arr_dim starts here
    ! =====================================================================

    SubroutineName = 'nc_get_arr_dim'
    IF(PRESENT(SubroutineName_in)) SubroutineName = SubroutineName_in
    ModName = ModuleName
    IF(PRESENT(ModuleName_in)) ModName = ModuleName_in
    
    ! Attach variable
    rcode = nf_inq_varid(ncid,TRIM(ADJUSTL(varname)), vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)),   &
                                TRIM(ADJUSTL(SubroutineName)),        &
                               'nf_inq_varid:'//TRIM(ADJUSTL(varname)))

    ! Get Number of dimensions
    rcode = nf_inq_varndims(ncid,vid,ndim)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)),   &
                                TRIM(ADJUSTL(SubroutineName)),        &
                               'nf_inq_varndims:'//TRIM(ADJUSTL(varname)))

    ! Get dimension IDs
    ALLOCATE(dimid(ndim))
    rcode = nf_inq_vardimid(ncid, vid, dimid)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)),   &
                                TRIM(ADJUSTL(SubroutineName)),        &
                               'nf_inq_vardimid:'//TRIM(ADJUSTL(varname)))

    ! Now get the index
    rcode = nf_inq_dim(ncid, dimid(dim_idx),tmpchar, nc_get_arr_dim)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)), &
                                TRIM(ADJUSTL(SubroutineName)),      &
                               'nf_inq_dim:'//TRIM(ADJUSTL(varname)))

  END FUNCTION nc_get_arr_dim

  INTEGER(KIND=2) FUNCTION nc_load_surf_i2(ncid,varname,i,j,Error,        &
                                           ModuleName_in,SubroutineName_in)

    ! ------------------
    ! function arguments
    ! ------------------
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: varname
    INTEGER,                    INTENT(IN)    :: i, j
    TYPE(ErrorType),            INTENT(INOUT) :: Error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: ModuleName_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: SubroutineName_in

    ! ---------------
    ! local variables
    ! ---------------
    INTEGER                :: rcode, vid, ndim

    CHARACTER(LEN=maxChar) :: SubroutineName,ModName,tmpchar

    ! =====================================================================
    ! nc_load_surf_i2 starts here
    ! =====================================================================

    ! Set names to pass to error module
    SubroutineName = 'nc_load_surf_i2'
    IF(PRESENT(SubroutineName_in)) SubroutineName = SubroutineName_in
    ModName = ModuleName
    IF(PRESENT(ModuleName_in)) ModName = ModuleName_in

    ! Get variable
    rcode = nf_inq_varid(ncid,TRIM(ADJUSTL(varname)), vid)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)),   &
                                TRIM(ADJUSTL(SubroutineName)),        &
                               'nf_inq_varid:'//TRIM(ADJUSTL(varname)))

    ! Get the variable
    rcode = nf_get_vara_int2(ncid,vid,(/i,j/),(/1,1/),nc_load_surf_i2)
    CALL CheckNetCDFErrorStatus(Error,rcode,TRIM(ADJUSTL(ModName)),  &
                                TRIM(ADJUSTL(SubroutineName)),       &
                               'nf_get_vara:'//TRIM(ADJUSTL(varname)))
  
  END FUNCTION nc_load_surf_i2

  
END MODULE netcdf_module
