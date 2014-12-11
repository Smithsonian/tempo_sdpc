MODULE OMI_copyHE4toHE5_class
   USE HE4_class
   USE HE5_class
   USE PGS_PC_class
   USE OMI_SMF_class
   USE ISO_C_BINDING, ONLY: C_LONG
   IMPLICIT NONE

   INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0
   INTEGER (KIND=4), PARAMETER, PRIVATE :: four = 4
   INTEGER (KIND=4), PARAMETER, PRIVATE :: NMAX = 100
   INTEGER (KIND=4 ), PRIVATE :: nattrs  = 0  
   INTEGER (KIND=4 ), PRIVATE :: ninflds = 0  
   CHARACTER ( LEN=PGSd_PC_VALUE_LENGTH_MAX ), DIMENSION(NMAX), PRIVATE :: &
                                               attrName, fieldName, dimName 
   LOGICAL (KIND=4 ), DIMENSION(NMAX), PRIVATE :: skipit 
   INTEGER (KIND=4 ), DIMENSION(NMAX), PRIVATE :: ntype, nbytes, count, rank
   INTEGER (KIND=4 ), DIMENSION(NMAX,7), PRIVATE :: dims 
   !INTEGER (KIND=4 ), DIMENSION(7), PRIVATE :: start, stride, edge
   INTEGER (KIND=C_LONG ), DIMENSION(7), PRIVATE :: start, stride, edge
   INTEGER (KIND=4 ), DIMENSION(0:NMAX), PRIVATE :: nBytesAccu
   INTEGER (KIND=1 ), DIMENSION(:), ALLOCATABLE, PRIVATE :: memBuf
   CHARACTER ( LEN=1 ), DIMENSION(NMAX), PRIVATE :: grp
   CHARACTER ( LEN=PGSd_PC_VALUE_LENGTH_MAX ), DIMENSION(NMAX), PRIVATE::dimlist

   PUBLIC :: OMI_readHE4swAttr
   PUBLIC :: OMI_cpwtHE5glAttr
   PUBLIC :: OMI_readHE4fields
   PUBLIC :: OMI_cpwtHE5fields
   PUBLIC :: OMI_4he2he5NumberType
   CONTAINS
      FUNCTION OMI_readHE4fields( he4swfLUN, he4swn, fieldList ) RESULT (status)
        INTEGER (KIND=4 ), INTENT( IN ) :: he4swfLUN
        CHARACTER ( LEN = * ), INTENT(IN) :: he4swn, fieldList
        INTEGER (KIND=4 ) :: status, ierr, version
        INTEGER (KIND=4 ) :: i, k, j, strbufsize
        INTEGER (KIND=4 ) :: swfid4, SWid4, nswath4, ngflds, ndflds
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: swathFN4
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: swathlist4, swathname4
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ), DIMENSION(NMAX) :: swN
        CHARACTER ( LEN = 1000 ) :: gfldlist, dfldlist
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

        version = 1
        status = PGS_PC_GetReference( he4swfLUN, version, swathFN4 )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I9)' ) "get filename from PCF file at LUN =", &
                                 he4swfLUN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "OMI_readHE4fields", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        ! open HDF-EOS swath file
        swfid4 = swopen( swathFN4, DFACC_READ )
        IF( swfid4 == -1 ) THEN
           WRITE( msg,'(A)' ) "open "// TRIM(swathFN4)//" for read failed"
           ierr = OMI_SMF_setmsg( OZT_E_FILE_OPEN, msg, &
                                 "OMI_readHE4fields", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        nswath4 = swinqswath( swathFN4, swathlist4, strbufsize )
        IF( nswath4 > 1 ) THEN
           IF( INDEX( swathlist4, TRIM(he4swn) ) == 0 ) THEN
              WRITE( msg,'(A)' ) "Can't find "// TRIM( he4swn ) // &
                                 "in " // TRIM( swathlist4)//","//TRIM(swathFN4)
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg,"OMI_readHE4fields",zero )
              status = OZT_E_FAILURE
              RETURN
           ELSE
              nswath4 = EH_parsestrF( swathlist4, ',', swN )
              DO j = 1, nswath4
                IF( INDEX( swN(j), TRIM(he4swn) ) > 0 ) THEN
                   swathname4 = swN(j)
                   EXIT
                ENDIF
              ENDDO
           ENDIF
        ELSE
           swathname4 = swathlist4(1:strbufsize)
        ENDIF

        ! Attach to the swath
        SWid4 = swattach( swfid4, swathname4 )
        IF( SWid4 == FAIL ) THEN
           WRITE( msg,'(A)' ) "Can't attach "// TRIM( swathname4 ) // &
                               "to file" //TRIM(swathFN4)
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "OMI_readHE4fields", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        ninflds = EH_parsestrF( fieldList, ',', fieldName )
        IF( ninflds <= 0 ) THEN
           status = swdetach( SWid4 )
           status = swclose( swfid4 )
           RETURN
        ENDIF
       
        ngflds = swinqgflds( SWid4, gfldlist, rank, ntype )
        ndflds = swinqdflds( SWid4, dfldlist, rank, ntype )

        nBytesAccu(0) = 0 
        DO i = 1, ninflds
          status = swfldinfo( SWid4, fieldName(i), rank(i), dims(i,:), &
                              ntype(i), dimlist(i) )
          IF( status == FAIL ) THEN
             WRITE( msg,'(A)' ) "Can't get info "// TRIM( fieldName(i) ) //&
                                " from "// TRIM(swathname4) // TRIM(swathFN4)//&
                                 ", skip this field."
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg,"OMI_readHE4fields",zero )
             skipit(i) = .TRUE.
             nbytes(i) = 0
          ELSE
             skipit(i) = .FALSE.
             nbytes(i) = DFKNTsize( ntype(i) )
             DO j = 1, rank(i)
               nbytes(i) = nbytes(i)*dims(i,j)
             ENDDO
          ENDIF
          nBytesAccu(i) = nBytesAccu(i-1) + nbytes(i)
          IF( INDEX( gfldlist, TRIM( fieldName(i) )) >  0 ) THEN
             grp(i) = 'g'
          ELSE IF( INDEX( dfldlist, TRIM( fieldName(i) )) >  0 ) THEN
             grp(i) = 'd'
          ENDIF
        ENDDO

        IF( ALLOCATED( memBuf ) ) DEALLOCATE( memBuf )
        ALLOCATE( memBuf( nBytesAccu(ninflds) ), STAT = ierr )
        IF( ierr /= zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, "memBuf", &
                                 "OMI_readHE4fields", zero )
           RETURN
        ENDIF
           
        start(:)  = 0
        stride(:) = 1
        DO i = 1, ninflds
          IF( skipit(i) ) CYCLE
          j = nBytesAccu(i-1)+1
          k = nBytesAccu(i)
          edge(1:rank(i)) = dims(i,1:rank(i))
          status = swrdfld( SWid4, fieldName(i), &
                            start, stride, edge, memBuf(j:k) ) 
          IF( status == FAIL ) THEN
             WRITE( msg,'(A)' ) "Read info "// TRIM( fieldName(i) ) //&
                                 " from "// TRIM(swathname4) // TRIM(swathFN4)
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg,"OMI_readHE4fields", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF
        ENDDO

        status = swdetach( SWid4 )
        status = swclose( swfid4 )

        status = OZT_S_SUCCESS
        RETURN

      END FUNCTION OMI_readHE4fields

      FUNCTION OMI_cpwtHE5Fields( he5swfLUN, he5swn ) RESULT (status)

        USE ISO_C_BINDING, ONLY: C_LONG

        INTEGER (KIND=4 ), INTENT( IN ) :: he5swfLUN 
        CHARACTER ( LEN = * ), INTENT(IN) :: he5swn
        INTEGER (KIND=4 ) :: status, ierr, version
        INTEGER (KIND=4 ) :: i, j, k, strbufsize, ndim
        INTEGER (KIND=4 ) :: swfid5, SWid5, nswath5
        INTEGER (KIND=4 ) :: numtype, nt, rankLocal
        INTEGER (KIND=4 ) :: size_numtype, size_nt
        INTEGER (KIND=4 ), DIMENSION(7) :: dimsLocal 
        !INTEGER (KIND=4 ), DIMENSION(NMAX) :: dimSize
        INTEGER (KIND=C_LONG ), DIMENSION(NMAX) :: dimSize
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: swathFN5 
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: swathlist5, swathname5
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ), DIMENSION(NMAX) :: swN
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: dimNamesList, &
                                                        dimens, maxdims 

        version = 1
        status = PGS_PC_GetReference( he5swfLUN, version, swathFN5 )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I9)' ) "get filename from PCF file at LUN =", &
                                 he5swfLUN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "OMI_cpwtHE5Fields", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        swfid5 = HE5_swopen( swathFN5, HE5F_ACC_RDWR )
        IF( swfid5 == FAIL ) THEN
           WRITE( msg,'(A)' ) "open "// TRIM(swathFN5)//" for read/write failed"
           ierr = OMI_SMF_setmsg( OZT_E_FILE_OPEN, msg, &
                                 "OMI_cpwtHE5fields", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        nswath5 = HE5_swinqswath( swathFN5, swathlist5, strbufsize )
        IF( nswath5 > 1 ) THEN
           IF( INDEX( swathlist5, TRIM(he5swn) ) == 0 ) THEN
              WRITE( msg,'(A)' ) "Can't find "// TRIM( he5swn ) // &
                                 "in " // TRIM( swathlist5)//","//TRIM(swathFN5)
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg,"OMI_cpwtHE5fields",zero )
              status = OZT_E_FAILURE
              RETURN
           ELSE
              nswath5 = EH_parsestrF( swathlist5, ',', swN )
              DO j = 1, nswath5
                IF( INDEX( swN(j), TRIM(he5swn) ) > 0 ) THEN
                   swathname5 = swN(j)
                   EXIT
                ENDIF
              ENDDO
           ENDIF
        ELSE
           swathname5 = swathlist5(1:strbufsize)
        ENDIF

        SWid5 = HE5_swattach( swfid5, swathname5 )
        IF( SWid5 == FAIL ) THEN
           WRITE( msg,'(A)' ) "Can't attach "// TRIM( swathname5 ) // &
                               "to file" //TRIM(swathFN5)
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg,"OMI_cpwtHE5fields", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        status = HE5_swinqdims( SWid5, dimNamesList, dimSize )
        IF( status == FAIL ) THEN
           WRITE( msg,'(A)' ) "HE5_swinqdims failed "// TRIM( swathname5 ) // &
                               "to file" //TRIM(swathFN5) 
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg,"OMI_cpwtHE5fields", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        DO i = 1, ninflds
          IF( skipit(i) ) CYCLE
 
          !! check to see if the field is already defined
          status = HE5_swfldinfo( SWid5,  fieldName(i), rankLocal, &
                                  dimsLocal, nt, dimens, maxdims )
          IF( status == 0 ) THEN
             !! if yes, now check data type, rank, and dim size to see 
             !! if they are the same, if yes, over-written
             !! whatever that is in it.
             numtype = OMI_4he2he5NumberType( ntype(i) )
             size_numtype = HE5Tget_size( numtype )
             size_nt      = HE5Tget_size( nt )
             IF( size_numtype /= size_nt .OR. rankLocal /= rank(i) ) THEN
                WRITE( msg,'(A,2I9,A,2I9,A)' ) "numberSize: (",size_numtype, &
                       size_nt, ") different or ranks: (", rankLocal, &
                       rank(i), ") mis-match between datafield defined "//&
                       "in HE4 and HE5 swaths."
                ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                       "OMI_cpwtHE5fields", zero )
                status = OZT_E_FAILURE
                RETURN
             ENDIF

             DO j = 1, rankLocal
               IF( dimsLocal(j) /= dims(i,j) ) THEN
                  WRITE( msg,'(A)' ) "dim size mis-mathch between "//&
                                     "datafield defined in HE4 and HE5 swaths."
                  ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, &
                                       "OMI_cpwtHE5fields", zero )
                  status = OZT_E_FAILURE
                  RETURN
               ENDIF
             ENDDO
          ELSE  
             !! if datafiled not already in the file, crate it.
             ndim = EH_parsestrF( dimlist(i), ',', dimName )
             DO j = 1, ndim
               IF( INDEX( dimNamesList, TRIM(dimName(j)) ) == 0 ) THEN
                  status = HE5_swdefdim( SWid5, dimName(j), dimSize(j) )
                  IF( status == FAIL ) THEN
                     WRITE( msg,'(A)' ) "HE5_swdefdim failed "//  &
                            TRIM(dimName(j))// "in swath:"//TRIM(swathname5)// &
                           "in file" //TRIM(swathFN5) 
                     ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                           "OMI_cpwtHE5fields", zero )
                     status = OZT_E_FAILURE
                     RETURN
                  ELSE
                     dimNamesList = TRIM(dimNamesList)//','//TRIM(dimName(j))
                  ENDIF
               ENDIF
             ENDDO

             numtype = OMI_4he2he5NumberType( ntype(i) )
             status = FAIL
             IF( grp(i) == 'g' ) THEN
                status = he5_swdefgfld( SWid5, fieldName(i), dimlist(i), " ", &
                                        numtype, HE5_HDFE_NOMERGE )
             ELSE IF( grp(i) == 'd' ) THEN
                status = he5_swdefdfld( SWid5, fieldName(i), dimlist(i), " ", &
                                        numtype, HE5_HDFE_NOMERGE )
             ENDIF
             IF( status == FAIL ) THEN
                WRITE( msg,'(A)') "Can't define field: "//TRIM(fieldName(i))//&
                                 "in "// TRIM(swathname5) // TRIM(swathFN5) 
                ierr = OMI_SMF_setmsg( OZT_E_HDFEOS,msg, "OMI_cpwtHE5fields", &
                                       zero )
                status = OZT_E_FAILURE
                RETURN
             ENDIF
          ENDIF
        ENDDO

        start(:)  = 0
        stride(:) = 1
        DO i = 1, ninflds
          IF( skipit(i) ) CYCLE
          j = nBytesAccu(i-1)+1
          k = nBytesAccu(i)
          edge(1:rank(i)) = dims(i,1:rank(i))
          status = HE5_swwrfld( SWid5, fieldName(i), &
                                start, stride, edge, memBuf(j:k) ) 
          IF( status == FAIL ) THEN
             WRITE( msg,'(A)' ) "Write filed "// TRIM( fieldName(i) ) //&
                                 "from "// TRIM(swathname5) // TRIM(swathFN5)
             ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg,"OMI_cpwtHE5fields",zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF
        ENDDO

        status = HE5_swdetach( SWid5 )
        status = HE5_swclose( swfid5 )
        IF( ALLOCATED( memBuf ) ) DEALLOCATE( memBuf )
        status = OZT_S_SUCCESS
      END FUNCTION OMI_cpwtHE5fields

      FUNCTION OMI_readHE4swAttr( he4swfLUN ) RESULT (status)
        INTEGER (KIND=4 ), INTENT( IN ) :: he4swfLUN 
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: swathFN 
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: swathlist, swathname
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: attrlist
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
        INTEGER (KIND=4 ) :: status, ierr, version
        INTEGER (KIND=4 ) :: i, k, j
        INTEGER (KIND=4 ) :: swfid, SWid, nswath, strbufsize

        version = 1
        status = PGS_PC_GetReference( he4swfLUN, version, swathFN )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I9)' ) "get filename from PCF file at LUN =", &
                                 he4swfLUN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "OMI_readHE4swAttr", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        ! open HDF-EOS swath file
        swfid = swopen( swathFN, DFACC_READ )
        IF( swfid == -1 ) THEN
           WRITE( msg,'(A)' ) "open "// TRIM(swathFN)//" for read failed"
           ierr = OMI_SMF_setmsg( OZT_E_FILE_OPEN, msg, &
                                 "OMI_readHE4swAttr", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF
 
        nswath = swinqswath( swathFN, swathlist, strbufsize )
        IF( nswath /= 1 ) THEN
           WRITE( msg,'(A,I9,A)' ) "only 1 swath is expected in "//TRIM(swathFN)&
                                  // " but ", nswath, " found."
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "OMI_readHE4swAttr", four )
           nattrs = 0
           status = OZT_E_FAILURE
           RETURN
        ELSE 
           swathname = swathlist(1:strbufsize)
        ENDIF
 
        ! Attach to the swath
        SWid = swattach( swfid, swathname )

        nattrs = swinqattrs( SWid, attrlist, strbufsize )
  
        IF( nattrs <= 0 ) RETURN
        IF( nattrs > NMAX ) THEN
           WRITE( msg,'(A,I9)' ) "too many attr in "//TRIM(swathname)//","// &
                                   TRIM(swathFN)//":", nattrs
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "OMI_readHE4swAttr", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        nattrs = EH_parsestrF( attrlist, ',', attrName )
        IF( nattrs <= 0 ) RETURN

        nBytesAccu(0) = 0 
        DO i = 1, nattrs
          status = swattrinfo( SWid, attrName(i), ntype(i), nbytes(i) )
          IF( status == FAIL ) THEN
             WRITE( msg,'(A)' ) " swattrinfo:"//TRIM(attrName(i))//","// &
                                  TRIM(swathFN) // ", skip this attribute."
             ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg,"OMI_readHE4swAttr",zero )
             skipit(i) = .TRUE.
             count(i) = 0
             nbytes(i) = 0 
             CYCLE
          ELSE
             skipit(i) = .FALSE.
             count(i) = nbytes(i)/DFKNTsize(ntype(i))
          ENDIF
          nBytesAccu(i) = nBytesAccu(i-1) + nbytes(i)
        ENDDO


        IF( ALLOCATED( memBuf ) ) DEALLOCATE( memBuf )
        ALLOCATE( memBuf( nBytesAccu(nattrs) ), STAT = ierr )
        IF( ierr /= zero ) THEN
           status = OZT_E_FAILURE
           ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, "memBuf", &
                                 "OMI_readHE4swAttr", zero )
           RETURN
        ENDIF

        DO i = 1, nattrs
          IF( skipit(i) ) CYCLE
          j = nBytesAccu(i-1)+1
          k = nBytesAccu(i)
          status = swrdattr( SWid, attrName(i), memBuf(j:k) )
        ENDDO

        status = swdetach( SWid )
        status = swclose( swfid )
        status = OZT_S_SUCCESS  
  
      END FUNCTION OMI_readHE4swAttr

      FUNCTION OMI_cpwtHE5glAttr( he5swfLUN ) RESULT (status)
        INTEGER (KIND=4 ), INTENT( IN ) :: he5swfLUN 
        CHARACTER ( LEN = PGSd_PC_VALUE_LENGTH_MAX ) :: swathFN 
        CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
        INTEGER (KIND=4 ) :: status, ierr, version, numtype
        INTEGER (KIND=4 ) :: swfid
        INTEGER (KIND=4 ) :: i, j, k

        !! Get the L2 file name from the PCF.
        version = 1
        status = PGS_PC_getreference( he5swfLUN, version, swathFN )
        IF( status /= PGS_S_SUCCESS ) THEN
           WRITE( msg,'(A,I9)' ) "get filename from PCF file at LUN =", &
                                 he5swfLUN
           ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "OMI_cpwtHE5glAttr", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF

        swfid = he5_swopen( swathFN, HE5F_ACC_RDWR )
        IF( swfid == -1 ) THEN
           WRITE( msg,'(A)' ) "he5_swopen:"// TRIM(swathFN) // " failed."
           ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, "OMI_cpwtHE5glAttr", zero )
           status = OZT_E_FAILURE
           RETURN
        ENDIF
                
        DO i = 1, nattrs
          IF( INDEX( attrName(i), "NumTimes" ) > 0 .OR. &
              INDEX( attrName(i), "EarthSunDistance" ) > 0 .OR. &
              skipit(i) ) CYCLE  !skip copy NumTimes
          j = nBytesAccu(i-1) + 1
          k = nBytesAccu(i)
          numtype = OMI_4he2he5NumberType( ntype(i) )
          IF( numtype == -1 ) THEN
             WRITE( msg,'(A, I9)' ) "Unknown HE4 number type:", ntype(i)
             ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                    "OMI_cpwtHE5glAttr", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF

          ierr = he5_ehwrglatt( swfid, attrName(i), &
                                numtype, count(i), memBuf(j:k) )
          IF( ierr == -1 ) THEN
             WRITE( msg,'(A)' ) "he5_ehwrglatt:"// TRIM(attrName(i)) //&
                                " failed in file " // TRIM(swathFN )
             ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
                                    "OMI_cpwtHE5glAttr", zero )
             status = OZT_E_FAILURE
             RETURN
          ENDIF
        ENDDO
        status = he5_swclose( swfid )
        IF( ALLOCATED( memBuf ) ) DEALLOCATE( memBuf )
        status = OZT_S_SUCCESS 
        RETURN
      END FUNCTION OMI_cpwtHE5glAttr

      FUNCTION OMI_4he2he5NumberType( nType_he4 ) RESULT (nType_he5)
        INTEGER (KIND=4 ), INTENT( IN ) :: nType_he4 
        INTEGER (KIND=4 ) :: nType_he5 
        
        SELECT CASE( nType_he4 )
          CASE ( DFNT_CHAR )
            nType_he5 = HE5T_NATIVE_CHAR
          CASE ( DFNT_UCHAR )
            nType_he5 = HE5T_NATIVE_UCHAR
          CASE ( DFNT_INT8 )
            nType_he5 = HE5T_NATIVE_INT8
          CASE ( DFNT_UINT8 )
            nType_he5 = HE5T_NATIVE_UINT8
          CASE ( DFNT_INT16 )
            nType_he5 = HE5T_NATIVE_INT16
          CASE ( DFNT_UINT16 )
            nType_he5 = HE5T_NATIVE_UINT16
          CASE ( DFNT_INT32 )
            nType_he5 = HE5T_NATIVE_INT32
          CASE ( DFNT_UINT32 )
            nType_he5 = HE5T_NATIVE_UINT32
          CASE ( DFNT_INT64 )
            nType_he5 = HE5T_NATIVE_INT64
          CASE ( DFNT_UINT64 )
            nType_he5 = HE5T_NATIVE_UINT64
          CASE ( DFNT_FLOAT32 )
            nType_he5 = HE5T_NATIVE_FLOAT
          CASE (DFNT_FLOAT64 )
            nType_he5 = HE5T_NATIVE_DOUBLE
          CASE ( DFNT_FLOAT128 )
            nType_he5 = HE5T_NATIVE_LDOUBLE
          CASE DEFAULT
            nType_he5 = -1
        END SELECT
        RETURN
      END FUNCTION OMI_4he2he5NumberType
      
END MODULE OMI_copyHE4toHE5_class
