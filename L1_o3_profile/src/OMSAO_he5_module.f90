MODULE OMSAO_he5_module
  ! ==============================
  ! Module for HDF-EOS5 parameters
  ! ==============================
  USE OMSAO_precision_module,  ONLY: i4, i8
  USE OMSAO_parameters_module, ONLY: &
       maxchlen, zero_r8, one_r8, &
       i2_missval, i4_missval, r4_missval, r8_missval, blank13, &
       valid_max_i2, valid_max_i4, valid_min_r8, valid_max_r8, &
       elsunc_usrstop_eval_r8, elsunc_highest_eval_r8, &
       main_qa_min_flag_r8, main_qa_max_flag_r8
  USE OMSAO_indices_module,    ONLY: &
       n_sao_pge, max_calfit_idx, sao_pge_min_idx, sao_pge_max_idx, &
       o3_t1_idx, o3_t3_idx
  USE ISO_C_BINDING

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  ! Integer variables for writing to swath data and geolocation fields
  ! ------------------------------------------------------------------
  INTEGER (KIND=C_LONG)  :: &
    he5_start_1d, he5_stride_1d, he5_edge_1d
  INTEGER (KIND=C_LONG), DIMENSION (2) :: &
    he5_start_2d, he5_stride_2d, he5_edge_2d
  INTEGER (KIND=C_LONG), DIMENSION (3) :: &
    he5_start_3d, he5_stride_3d, he5_edge_3d
  INTEGER (KIND=C_LONG), DIMENSION (6) :: &
    he5_start_6d, he5_stride_6d, he5_edge_6d
  INTEGER (KIND=C_LONG), DIMENSION (7) :: &
    he5_start_7d, he5_stride_7d, he5_edge_7d
  
  ! The following structures are based on Kai Yang OMTO3SHR routines
  ! HE5_class and L2_data_structure
  INTEGER (KIND=i4), PARAMETER :: HE5_DTSETRANKMAX = 8
  INTEGER (KIND=i4), PARAMETER :: HE5_FLDNUMBERMAX = 500
  INTEGER (KIND=i4), PARAMETER :: MAXRANK = 4
  TYPE, PUBLIC :: DFHE5_T
      REAL (KIND=i8)        :: ValidRange_l, ValidRange_h, ScaleFactor, Offset
      CHARACTER (LEN = 80)  :: name
      CHARACTER (LEN = 256) :: dimnames
      CHARACTER (LEN = 80)  :: Units
      CHARACTER (LEN = 256) :: LongName
      CHARACTER (LEN = 256) :: UniqueFieldDefinition
      INTEGER (KIND = i4)   :: swath_id
      INTEGER (KIND = i4)   :: datatype
      INTEGER (KIND = i4)   :: rank
      INTEGER (KIND = i4), DIMENSION(MAXRANK) :: dims
  END TYPE DFHE5_T

  INTEGER (KIND=4), PARAMETER :: NDIM_MAX = 12 ! 22 ?
  INTEGER (KIND=4), PARAMETER :: NFLDS_MAX = 75 ! 100  ?
  INTEGER (KIND=4), PARAMETER :: MAX_STR_LEN = 512

  TYPE, PUBLIC :: L2_generic_type
      CHARACTER ( LEN = MAX_STR_LEN )             :: filename, swathname
      CHARACTER ( LEN = MAX_STR_LEN )             :: dimnames
      CHARACTER ( LEN = 80 ), DIMENSION(NFLDS_MAX):: fieldname
      INTEGER (KIND = C_LONG)  :: sw_fid, swathID
      INTEGER (KIND = C_LONG)  :: nDims, nFields
      INTEGER (KIND = C_LONG)  :: iLine, eLine, nLine, nTotLine
      INTEGER (KIND = C_LONG), DIMENSION(NDIM_MAX)    :: dimSizes
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX)   :: rank
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX, MAXRANK) :: dims
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX)   :: elmSize,  pixSize
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX)   :: lineSize, blkSize
      INTEGER (KIND = C_LONG)                         :: SumElmSize, SumLineSize
      INTEGER (KIND = C_LONG), DIMENSION(0:NFLDS_MAX) :: accuLineSize, accuBlkSize
      LOGICAL                                         :: initialized

      ! geo and data fields
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: data
  END TYPE L2_generic_type

  INCLUDE 'hdfeos5.inc'
  ! ---------------------------------------------------------------
  ! Finally some system definitions that come with the HE5 Library.
  ! We include this here so that we don't have to worry about it
  ! inside the subroutines that use it.
  ! ---------------------------------------------------------------
  INTEGER (KIND = i4), EXTERNAL :: &
       he5_ehrdglatt, he5_ehwrglatt, he5_swattach,  he5_swclose,   he5_swcreate, &
       he5_swdefdfld, he5_swdefdim,  he5_swdefgfld, he5_swdetach,  he5_swopen,   &
       he5_swrdattr,  he5_swrdfld,   he5_swrdgattr, he5_swrdlattr, he5_swwrattr, &
       he5_swwrfld,   he5_swwrgattr, he5_swwrlattr, he5_swfldinfo, he5_ehglattinf, &
       he5_SWinqswath,he5_SWinqdims, he5_SWinqgflds,he5_SWinqdflds,he5_ehinqglatts, &
       he5_swdefchunk, he5_swdefcomp, &
       he5_gdopen, he5_gdattach, he5_gdrdfld, he5_gddetach, he5_gdclose

  PUBLIC :: EH_parsestrF
  CONTAINS
    FUNCTION EH_parsestrF( instring, delim, outstrs, strln ) RESULT( nstrout )
      CHARACTER ( LEN = * ), INTENT(IN) :: instring
      CHARACTER ( LEN = 1 ), INTENT(IN) :: delim
      CHARACTER ( LEN = * ), DIMENSION(:), INTENT(OUT) :: outstrs
      INTEGER (KIND=4 ), DIMENSION(:), INTENT(OUT), OPTIONAL :: strln
      CHARACTER ( LEN = LEN( instring) ):: localStr
      INTEGER (KIND=4 ) :: i, j, k, nstrout
      INTEGER (KIND=4 ) :: sOut, strlnS
      
      !! input string is empty
      IF( LEN_TRIM( instring ) == 0 ) THEN
         nstrout  = 0
         IF( PRESENT( strln ) ) strln(1) = 0
         RETURN
      ENDIF
      
      sOut = SIZE( outstrs )
      IF( PRESENT( strln ) ) strlnS = SIZE( strln )
      
      IF( LEN_TRIM( delim ) == 0 ) THEN   ! input string not empty
         outstrs(1) = instring            ! but delim is empty
         nstrout    = 1
         IF( PRESENT( strln ) ) strln(1) = LEN_TRIM( instring )
      ELSE                                ! delim is not empty
         localStr = instring
         i = 1
         DO WHILE( INDEX( localStr, delim ) > 0 )
            j = LEN_TRIM( localStr )
            k = INDEX( localStr, delim )
            outstrs(i) = localStr( 1:k-1 )
            localStr   = localStr( k+1:j )
            IF( PRESENT( strln ) ) strln(i) = k - 1
            i = i+1
            IF( i > sOut ) THEN
               nstrout = -1      !! error return when outstrs array is
               RETURN            !! not large enough to hold the results
            ENDIF
            
            IF( PRESENT( strln ) ) THEN
               strlnS = SIZE( strln )
               IF( i > strlnS ) THEN
                  nstrout = -1
                  RETURN
               ENDIF
            ENDIF
         ENDDO
         outstrs(i) = localStr
         nstrout    = i
         IF( PRESENT( strln ) ) strln(i) = LEN_TRIM( localStr )
      ENDIF
      
      RETURN
    END FUNCTION EH_parsestrF
    
END MODULE OMSAO_he5_module


