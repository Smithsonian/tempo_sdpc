MODULE HE5_class
  INCLUDE   'hdfeos5.inc'

  INTEGER (KIND=4), PARAMETER :: HE5_DTSETRANKMAX = 8
  INTEGER (KIND=4), PARAMETER :: HE5_FLDNUMBERMAX = 500
  INTEGER (KIND=4), PARAMETER, PRIVATE :: MAXRANK = 3
  TYPE, PUBLIC :: DFHE5_T
      REAL (KIND=8) :: ValidRange_l, ValidRange_h, ScaleFactor, Offset
      CHARACTER (LEN = 80)  :: name
      CHARACTER (LEN = 256) :: dimnames
      CHARACTER (LEN = 80)  :: Units
      CHARACTER (LEN = 256) :: LongName
      CHARACTER (LEN = 256) :: UniqueFieldDefinition
      INTEGER (KIND=4) :: swath_id
      INTEGER (KIND=4) :: datatype
      INTEGER (KIND=4) :: rank
      INTEGER (KIND=4), DIMENSION(MAXRANK) :: dims
  END TYPE DFHE5_T

  INTEGER (KIND = 4), EXTERNAL :: he5_swattach, &
                                  he5_swclose, &
                                  he5_swcreate, &
                                  he5_swinqswath, &
                                  he5_swinqdims, &
                                  he5_swinqgflds, &
                                  he5_swinqdflds, &
                                  he5_swdefdfld, &
                                  he5_swdefdim, &
                                  he5_swdefgfld, &
                                  he5_swdetach, &
                                  he5_swfldinfo, &
                                  he5_swopen, &
                                  he5_swrdattr, &
                                  he5_swrdfld, &
                                  he5_swrdgattr, &
                                  he5_swrdlattr, &
                                  he5_swwrattr, &
                                  he5_swwrfld, &
                                  he5_swwrgattr, &
                                  he5_swwrlattr, &
                                  he5_ehwrglatt, &
                                  he5_ehrdglatt, &
                                  he5_ehglattinf, &
                                  he5_ehinqglatts, &
                                  he5tget_size,    &
                                  he5_swdefchunk,  &
                                  he5_swdefcomp

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
        !IF( PRESENT( strln ) ) strlnS = SIZE( strln )

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

END MODULE HE5_class
