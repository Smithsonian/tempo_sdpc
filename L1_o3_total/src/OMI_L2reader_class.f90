MODULE OMI_L2reader_class
    USE L2_data_structure
    USE HE5_class
    USE UTIL_tools_class
    USE PGS_PC_class  ! define PGSd_PC_FILE_PATH_MAX, and PGS_PC functions 
                      ! include PGS_SMF.f define PGS_SMF_MAX_MSG_SIZE
    USE OMI_SMF_class ! include PGE specific messages and OMI_SMF_setmsg

    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER, PRIVATE :: zero = 0, one = 1, two = 2
    INTEGER (KIND=4), PARAMETER, PRIVATE :: three = 3, four =4, five = 5

    PUBLIC  :: L2_getFileNames
    PUBLIC  :: L2_newBlock
    PUBLIC  :: L2_readBlock
    PUBLIC  :: L2_getSWdim
    PUBLIC  :: L2_getLine
!    PUBLIC  :: L2_getNvADJinfo
!    PUBLIC  :: L2_getNvADJ
!    PUBLIC  :: L2_getYMD
!    PUBLIC  :: L2_getEarthSunDistance
    PUBLIC  :: L2_disposeBlock

    CONTAINS
       FUNCTION L2_getFileNames( L2_file_LUN, numfiles, L2_filenames, &
                                 L2_swathlist ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN) :: L2_file_LUN
         CHARACTER( LEN = * ), DIMENSION(:), INTENT(OUT) :: L2_filenames
         CHARACTER( LEN = * ),               INTENT(OUT) :: L2_swathlist
         INTEGER (KIND=4), INTENT(OUT) :: numfiles
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         CHARACTER( LEN = LEN(L2_swathlist) ) :: swathlist
         INTEGER (KIND=4) :: di, version, nswath, strbufsize
         INTEGER (KIND=4) :: status, ierr
         INTEGER (KIND=4) :: SW_fileid 

         !! Get number of L2 output files from PCF
         !! and make sure there is only one 

         status = PGS_PC_getnumberoffiles( L2_file_LUN,  numfiles )
         IF( status /= PGS_S_SUCCESS ) THEN
            WRITE( msg,'(A,I9)' ) "can't get numfiles from PCF file at LUN =", &
                                  L2_file_LUN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT,msg,"L2_getFileNames",zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF
         !! Get the L2 file name from the PCF.
         DO di = 1, numfiles
           version = di
           status = PGS_PC_getreference( L2_file_LUN, version, &
                                         L2_filenames(di) )
           IF( status /= PGS_S_SUCCESS ) THEN
              WRITE( msg,'(2(A,I9))' ) "get filename from PCF file at LUN =", &
                                      L2_file_LUN, " Version =", version
              ierr = OMI_SMF_setmsg( OZT_E_INPUT,msg,"L2_getFileNames",zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

           SW_fileid = HE5_SWopen( TRIM(L2_filenames(di)), &
                                   HE5F_ACC_RDONLY )
           IF( SW_fileid == -1 ) THEN
              WRITE( msg,'(A)' ) "HE5_SWopen:"// TRIM(L2_filenames(di))//&
                                 " failed."
              ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getFileNames",zero)
              status = OZT_E_FAILURE
              RETURN 
           ENDIF

           nswath = HE5_SWinqswath( TRIM(L2_filenames(di)), swathlist, &
                                    strbufsize )
           IF( di == 1 ) THEN
              L2_swathlist = swathlist(1:strbufsize)
           ELSE IF( TRIM(L2_swathlist) /= swathlist(1:strbufsize) ) THEN
              WRITE( msg,'(A)' ) "file contains different swathlist:"// &
                                 TRIM( L2_swathlist ) // "," // &
                                 swathlist(1:strbufsize) 
              ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getFileNames",zero)
              status = OZT_E_FAILURE
              RETURN 
           ENDIF 

           status = HE5_SWclose( SW_fileid )
         ENDDO
         status = OZT_S_SUCCESS
         RETURN 
       END FUNCTION L2_getFileNames

!       FUNCTION L2_getYMD( L2_filename, Year, Month, Day ) RESULT( status )
!         CHARACTER( LEN = * ), INTENT(IN) :: L2_filename
!         INTEGER (KIND=4), INTENT(OUT) :: Year, Month, Day
!         INTEGER (KIND=4) :: status, ierr
!         INTEGER (KIND=4) :: SW_fileid 
!         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
!
!         SW_fileid = HE5_SWopen( L2_filename, HE5F_ACC_RDONLY )
!         IF( SW_fileid == -1 ) THEN
!            WRITE( msg,'(A)' ) "HE5_SWopen:"// TRIM(L2_filename) // " failed."
!            ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getYMD",zero)
!            status = OZT_E_FAILURE
!            RETURN 
!         ENDIF
!
!         ierr = he5_ehrdglatt( SW_fileid, "GranuleDay", Day )
!         IF( ierr == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleDay" //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
!                                  "L2_getYMD", zero )
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!
!         ierr = he5_ehrdglatt( SW_fileid, "GranuleMonth", Month )
!         IF( ierr == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleMonth" //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
!                                  "L2_getYMD", zero )
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!         
!         ierr = he5_ehrdglatt( SW_fileid, "GranuleYear", Year )
!         IF( ierr == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "GranuleYear" //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
!                                  "L2_getYMD", zero )
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!         ierr = HE5_SWclose( SW_fileid )
!         status = OZT_S_SUCCESS
!         RETURN
!       END FUNCTION L2_getYMD

!       FUNCTION L2_getNvADJinfo( L2_filename, nWvc, nXtc ) RESULT( status )
!         CHARACTER( LEN = * ), INTENT(IN) :: L2_filename
!         INTEGER (KIND=4), INTENT(OUT) :: nWvc, nXtc
!         INTEGER (KIND=4) :: status, ierr, numberType, count
!         INTEGER (KIND=4) :: SW_fileid 
!         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
!
!         SW_fileid = HE5_SWopen( L2_filename, HE5F_ACC_RDONLY )
!         IF( SW_fileid == -1 ) THEN
!            WRITE( msg,'(A)' ) "HE5_SWopen:"// TRIM(L2_filename) // " failed."
!            ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getNvADJinfo",zero)
!            status = OZT_E_FAILURE
!            RETURN 
!         ENDIF
!         ierr = he5_ehglattinf( SW_fileid, "WavelengthOfAdjustment", &
!                                numberType,  count )
!
!         IF( ierr == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_ehglattinf:"// "WavelengthOfAdjustment" //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getNvADJinfo",zero)
!            status = OZT_E_FAILURE
!            ierr = HE5_SWclose( SW_fileid )
!            RETURN
!         ELSE
!            nWvc = count
!         ENDIF
!
!         ierr = he5_ehglattinf( SW_fileid, "NVXAdjustment", numberType, count)
!         IF( ierr == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_ehglattinf:"// "NVXAdjustment" //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getNvADJinfo",zero)
!            status = OZT_E_FAILURE
!            ierr = HE5_SWclose( SW_fileid )
!            RETURN
!         ELSE
!            nXtc = count/nWvc
!         ENDIf
!         ierr = HE5_SWclose( SW_fileid )
!         status = OZT_S_SUCCESS
!         RETURN
!       END FUNCTION L2_getNvADJinfo

!       FUNCTION L2_getNvADJ( L2_filename, crwl, swpcr ) RESULT( status )
!         CHARACTER( LEN = * ), INTENT(IN) :: L2_filename
!         REAL(KIND = 4), DIMENSION(:,:), INTENT(OUT) :: swpcr
!         REAL(KIND = 4), DIMENSION(:), INTENT(OUT) :: crwl
!         INTEGER (KIND=4) :: status, ierr
!         INTEGER (KIND=4) :: SW_fileid 
!         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
!
!         SW_fileid = HE5_SWopen( L2_filename, HE5F_ACC_RDONLY )
!         IF( SW_fileid == -1 ) THEN
!            WRITE( msg,'(A)' ) "HE5_SWopen:"// TRIM(L2_filename) // " failed."
!            ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getNvADJ",zero)
!            status = OZT_E_FAILURE
!            RETURN 
!         ENDIF
!
!         ierr = he5_ehrdglatt( SW_fileid, "WavelengthOfAdjustment", crwl )
!         IF( ierr == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "WavelengthOfAdjustment" //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
!                                  "L2_getNvADJ", zero )
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!
!         ierr = he5_ehrdglatt( SW_fileid, "NVXAdjustment", swpcr )
!         IF( ierr == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_ehwrglatt:"// "NVXAdjustment" //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
!                                  "L2_getNvADJ", zero )
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!
!         ierr = HE5_SWclose( SW_fileid )
!         status = OZT_S_SUCCESS
!         RETURN
!       END FUNCTION L2_getNvADJ

       FUNCTION L2_newBlock( this, filename, swathname, fieldlist, nL, &
                             he5accessTag ) RESULT( status )
         INTEGER (KIND=4), INTENT(IN), OPTIONAL :: nL
         INTEGER (KIND=4), INTENT(IN), OPTIONAL :: he5accessTag
         CHARACTER( LEN = * ), INTENT(IN) :: filename, swathname, fieldlist
         TYPE (L2_generic_type), INTENT( OUT ) :: this
         INTEGER (KIND=4) :: id, k, rankID
         CHARACTER( LEN = MAX_STR_LEN ) :: dimlist, maxdimlist
         INTEGER (KIND=4) :: ntype
         INTEGER (KIND=4) :: he5accessTag_l
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND=4) :: status, ierr

         this%filename  = filename 
         this%swathname = swathname
         IF( PRESENT( he5accessTag ) ) THEN
            he5accessTag_l = he5accessTag
         ELSE
            he5accessTag_l = HE5F_ACC_RDONLY 
         ENDIF

         !! open the file
         this%sw_fid = HE5_SWopen( this%filename, he5accessTag_l )
         IF( this%sw_fid == -1 ) THEN
            WRITE( msg,'(A)' ) "HE5_SWopen:"// TRIM(this%filename)// " failed."
            ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_newBlock",zero)
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         !! attach the swath 
         this%swathID = HE5_SWattach( this%sw_fid, this%swathname )
         IF( this%swathID == -1 ) THEN
            WRITE( msg, '(A)' ) "HE5_SWattach:" // TRIM( this%swathname ) // &
                                ", ", TRIM(this%filename)// " failed."
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, "L2_newBlock", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF

         !! get the dimensions associated with this swath
         this%nDims = HE5_SWinqdims( this%swathID, this%dimnames, &
                                     this%dimSizes )
         IF( this%nDims == -1 ) THEN
            WRITE( msg, '(A)' ) "HE5_SWinqdims:" // TRIM( this%swathname ) // &
                                ", ", TRIM(this%filename)// " failed."
            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, "L2_newBlock", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF
       
         !! Extract the fieldnames from the input list
         this%nFields =  EH_parsestrF( fieldlist, ',', this%fieldname )
         IF( this%nFields == 0 .OR. this%nFields > NFLDS_MAX ) THEN
            WRITE( msg, '(A)' ) "fieldlist:" // TRIM( fieldlist ) // &
                                " contains no or too many names."
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_newBlock", zero )
            status = OZT_E_FAILURE
            RETURN 
         ENDIF
             
         !! Find out the dimension size, the rank, and the datatype and 
         !! and byte size for each of the field (geo or data).
         this%SumElmSize = 0
         DO id = 1, this%nFields
           ierr = HE5_SWfldinfo( this%swathID, this%fieldname(id), rankID, &
                                 this%dims(id,:), ntype, dimlist, maxdimlist )

           IF( ierr == - 1 ) THEN
              WRITE( msg, '(A)') "HE5_SWfldinfo failed for "// &
                                 TRIM( this%fieldname(id)) // ", "// &
                                 TRIM( this%swathname )
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_newBlock", zero )
              status = OZT_E_FAILURE
              RETURN 
           ENDIF

           IF( rankID > 3 ) THEN
              WRITE( msg,'(A,I9,A)' ) "rank= ", rankID, " exceed 3 for " // &
                                     TRIM( this%fieldname(id) )
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_newBlock", zero )
              status = OZT_E_FAILURE
              RETURN
           ELSE
              this%rank(id) = rankID
           ENDIF
           
           this%elmSize(id) = HE5Tget_size( ntype )
           IF( this%elmSize(id) <= 0 ) THEN
              WRITE( msg,'(A,I9,A)' ) "datatype: ", ntype, &
                                     " error for " // TRIM( this%fieldname(id) )
              ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_newBlock", zero )
              status = OZT_E_FAILURE
              RETURN
           ELSE
              this%SumElmSize = this%SumElmSize + this%elmSize(id) 
           ENDIF

           !! make sure every field has the same first (HDF5) dimension size 
           !! Note that first HDF5 or C dimension is the last Fortran 
           !! dimension.
           IF( id == 1 ) THEN
              this%nTotLine = this%dims(1,this%rank(1))
           ELSE 
              IF( this%nTotLine /= this%dims(id,this%rank(id)) ) THEN
                 WRITE( msg, '(A,2I5)') "First dimension size not equal:", &
                        this%dims(id,this%rank(id)), this%nTotLine
                 ierr = OMI_SMF_setmsg( OZT_E_INPUT, msg, "L2_newBlock", zero )
                 status = OZT_E_FAILURE
                 RETURN 
              ENDIF
           ENDIF
         ENDDO

         !! the first C (HDF5) dimesnion is usually referred to as the line
         !! The block size is determined by the number of Lines in the block.
         IF( PRESENT(nL) ) THEN
            IF( nL < this%nTotLine ) THEN
               this%nLine = nL      !! use input nL (if present) and if it is 
            ELSE                    !! smaller than the dim size.
               this%nLine = this%nTotLine
            ENDIF 
         ELSE
            IF( this%nTotLine < 100 ) THEN
               this%nLine = this%nTotLine
            ELSE
               this%nLine = 100   !! default nLine is 100
            ENDIF
         ENDIF 
           
         !!calculate the line size and block size and the accumulated block size
         this%SumLineSize = 0
         this%accuLineSize(0) = 0
         this%accuBlkSize(0) = 0
         DO id = 1, this%nFields
           rankID = this%rank(id)
           IF( rankID == 1 ) THEN
               this%pixSize(id) = this%elmSize(id)
              this%lineSize(id) = this%elmSize(id)
           ELSE IF( rankID == 2 ) THEN
               this%pixSize(id) = this%elmSize(id)
              this%lineSize(id) = this%elmSize(id) * this%dims(id,1)
           ELSE
               this%pixSize(id) = this%elmSize(id) * this%dims(id,1)
              this%lineSize(id) = this%pixSize(id) * this%dims(id,2)
           ENDIF
           this%SumLineSize = this%SumLineSize + this%lineSize(id)
           this%blkSize(id) = this%lineSize(id) * this%nLine
           this%accuLineSize(id) = this%accuLineSize(id-1) + this%lineSize(id)
            this%accuBlkSize(id) =  this%accuBlkSize(id-1) +  this%blkSize(id)
         ENDDO

         IF( this%accuBlkSize(this%nFields) > 0 ) THEN
            IF( ASSOCIATED( this%data ) ) DEALLOCATE( this%data )
            ALLOCATE( this%data( this%accuBlkSize(this%nFields) ), STAT = ierr )
            IF( ierr /= zero ) THEN
               status = OZT_E_FAILURE
               ierr = OMI_SMF_setmsg( OZT_E_MEM_ALLOC, "this%data", &
                                      "L2_newBlock", zero )
               RETURN
            ENDIF
         ENDIF
         this%iLine = -1
         this%eLine = -1
         this%initialized = .TRUE.
         status = OZT_S_SUCCESS

         RETURN 
       END FUNCTION L2_newBlock

       FUNCTION L2_readBlock( this, iLine ) RESULT( status )

         USE ISO_C_BINDING, ONLY: C_LONG

         TYPE (L2_generic_type), INTENT( INOUT ) :: this
         INTEGER (KIND = 4), INTENT( IN ) :: iLine
         INTEGER (KIND = 4) :: ierr, status
         INTEGER (KIND = 4) :: ii, jj, id, k, rank, nL
         !INTEGER (KIND=4 ), DIMENSION(3) :: start, stride, edge
         INTEGER (KIND= C_LONG), DIMENSION(3) :: start, stride, edge
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg

         IF( .NOT. this%initialized ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "input block not initialized", &
                                  "L2_readBlock", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         IF( iLine < 0 .OR. iLine >= this%nTotLine ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                 "iLine out of range", "L2_readBlock", zero )
            status = OZT_E_FAILURE
            RETURN
         ELSE
            this%iLine = iLine
         ENDIF

         IF( (iLine + this%nLine) > this%nTotLine ) THEN
            nL = this%nTotLine - iLine
         ELSE
            nL = this%nLine
         ENDIF
         this%eLine = this%iLine + nL - 1

         stride(:) = 1
         DO id = 1, this%nFields
           rank = this%rank(id) 
           DO k = 1, rank
             IF( k == rank ) THEN
                start(k) = iLine
                edge(k)  = nL
             ELSE
                start(k) = 0
                edge(k)  = this%dims(id,k)
             ENDIF
           ENDDO
                 
           ii = this%accuBlkSize(id-1)+1
           jj = this%accuBlkSize(id)
           status = HE5_swrdfld( this%swathID, this%fieldname(id), &
                                 start, stride, edge, this%data(ii:jj) )
           IF( status == -1 ) THEN
              WRITE( msg,'(A)' ) "Read "//TRIM(this%fieldname(id))//&
                                  " failed."
              ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg,"OMI_readBlock",zero )
              status = OZT_E_FAILURE
              RETURN
           ENDIF

         ENDDO
         status = OZT_S_SUCCESS
         RETURN
       END FUNCTION L2_readBlock

       FUNCTION L2_getLine( this, iLine, lineMem ) RESULT( status )
         TYPE (L2_generic_type), INTENT( INOUT ) :: this
         INTEGER (KIND = 4), INTENT( IN ) :: iLine
         INTEGER (KIND = 1), DIMENSION(:), INTENT( OUT ) :: lineMem
         INTEGER (KIND = 4) :: ierr, status
         INTEGER (KIND = 4) :: j, id, ii, jj, ll, nn

         status = OZT_S_SUCCESS
         IF( .NOT. this%initialized ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "input block not initialized", &
                                  "L2_getLine", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         IF( iLine < 0 .OR. iLine >= this%nTotLine ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, "iLine out of range", &
                                   "L2_getLine", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         IF( iLine < this%iLine .OR. iLine > this%eLine ) THEN
            status = L2_readBlock( this, iLine )
            IF( status /= OZT_S_SUCCESS ) THEN
               ierr = OMI_SMF_setmsg( OZT_E_DATA_BLOCK, "retrieve data block", &
                                     "L2_getLine", zero )
               RETURN
            ENDIF
         ENDIF

         IF( SIZE( lineMem ) < this%SumLineSize ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                   "input lineMem size too small", &
                                   "L2_getLine", zero )
            status = OZT_E_FAILURE
            RETURN
         ENDIF

         j = iLine - this%iLine 
         DO id = 1, this%nFields
           ii = this%accuBlkSize(id-1)+j*this%lineSize(id)
           jj = ii + this%lineSize(id)
           ii = ii + 1  !! Fortran index 
           ll = this%accuLineSize(id-1) + 1
           nn = this%accuLineSize(id)
           lineMem(ll:nn) = this%data(ii:jj)  
         ENDDO   
         RETURN
       END FUNCTION L2_getLine

       FUNCTION L2_getSWdim( this, dimname ) RESULT( dimsize )
         TYPE (L2_generic_type), INTENT( IN ) :: this
         CHARACTER( LEN = * ), INTENT(IN) :: dimname
         CHARACTER( LEN = 80 ), DIMENSION(NDIM_MAX) :: dimnamesArray
         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
         INTEGER (KIND=4) :: ierr, id, nd, dimsize 
 
         IF( .NOT. this%initialized ) THEN
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, &
                                  "input block not initialized", &
                                  "L2_getSWdim", zero )
            dimsize = -1
            RETURN
         ENDIF

         nd =  EH_parsestrF( this%dimnames, ',', dimnamesArray )
         IF( nd /= this%nDims ) THEN
            WRITE( msg, '(A,I3,A,I3)') "nd =", nd, "nDims=", this%nDims
            ierr = OMI_SMF_setmsg( OZT_E_INPUT, TRIM( msg ) // &
                                  ":number of dimensions not consistent", &
                                  "L2_getSWdim", zero )
            dimsize = -1
            RETURN
         ENDIF
      
         DO id = 1, nd
           IF( TRIM(dimname) == TRIM( dimnamesArray(id) ) ) THEN 
              dimsize = this%dimSizes(id) 
              RETURN
           ENDIF
         ENDDO

         ierr = OMI_SMF_setmsg( OZT_W_GENERAL, &
                               "Unknown dimname:" // TRIM( dimname ), &
                               "L2_getSWdim", two )
         dimsize = -1
         RETURN
       END FUNCTION L2_getSWdim

       SUBROUTINE L2_disposeBlock( this )
         TYPE (L2_generic_type), INTENT( INOUT ) :: this
         INTEGER (KIND=4) :: ierr
         this%nDims   = 0
         this%nFields = 0
         this%iLine   = -1
         this%eLine   = -1
         this%nLine   =  0
         IF( ASSOCIATED( this%data ) ) DEALLOCATE( this%data )
         IF( this%initialized ) THEN
            ierr = HE5_SWdetach( this%swathID )
            ierr = HE5_SWclose( this%sw_fid )
         ENDIF
         this%initialized = .FALSE.
       END SUBROUTINE L2_disposeBlock
       
!       FUNCTION L2_getEarthSunDistance( L2_filename, L2_swathname, &
!                                        EarthSunDistance) RESULT( status )
!         CHARACTER( LEN = * ), INTENT(IN) :: L2_filename, L2_swathname
!         REAL (KIND=4), INTENT(OUT)      :: EarthSunDistance
!         INTEGER (KIND=4) :: status, ierr
!         INTEGER (KIND=4) :: SW_fileid, SW_id
!         CHARACTER( LEN = PGS_SMF_MAX_MSG_SIZE  ) :: msg
!
!         SW_fileid = HE5_SWopen( L2_filename, HE5F_ACC_RDONLY )
!         IF( SW_fileid == -1 ) THEN
!            WRITE( msg,'(A)' ) "HE5_SWopen:"// TRIM(L2_filename) // " failed."
!            ierr = OMI_SMF_setmsg(OZT_E_HDFEOS,msg,"L2_getEarthSunDistance",zero)
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!
!         !! attach the swath
!         SW_id = he5_swattach( SW_fileid, L2_swathname )
!         IF( SW_id == -1 ) THEN
!            WRITE( msg,'(A)' ) "he5_swattach:"// TRIM(L2_swathname) //&
!                               " failed in file " // TRIM(L2_filename )
!            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, msg, &
!                                  "L2_getEarthSunDistance", zero )
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!
!         status = he5_swrdattr( SW_id, "EarthSunDistance", EarthSunDistance )
!         IF( status == -1 ) THEN
!            ierr = OMI_SMF_setmsg( OZT_E_HDFEOS, &
!                            "Read Swath Attribute EarthSunDistance failed.", &
!                            "L2_getEarthSunDistance", zero )
!            status = OZT_E_FAILURE
!            RETURN
!         ENDIF
!
!         ierr = HE5_SWdetach( SW_id )
!         ierr = HE5_SWclose( SW_fileid )
!         status = OZT_S_SUCCESS
!         RETURN
!       END FUNCTION L2_getEarthSunDistance

END MODULE OMI_L2reader_class
