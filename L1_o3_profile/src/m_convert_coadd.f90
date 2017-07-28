!
module m_convert_coadd

  public convert_gpqualflag_info, convert_2bytes_to_16bits, &
       coadd_2bytes_qflgs, coadd_byte_qflgs, convert_xtrackqflag_info
  private convert_16bits_to_2bytes, convert_byte_to_8bits, &
       convert_8bits_to_byte

contains


  SUBROUTINE convert_gpqualflag_info ( &
       nxtrack, omi_geoflg, land_water_flg, glint_flg, snow_ice_flg )

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
!    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (IN) :: omi_geoflg
    INTEGER (KIND=i2), DIMENSION (:), INTENT (IN) :: omi_geoflg

    ! ----------------
    ! Output variables
    ! ----------------
!    INTEGER (KIND=i2), DIMENSION (nxtrack), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg
    INTEGER (KIND=i2), DIMENSION (:), INTENT (OUT) :: land_water_flg, glint_flg, snow_ice_flg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4),                PARAMETER      :: nbyte = 16
    INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = int((/ 1, 2, 4, 8, 16, 32, 64 /), kind=i2)
    INTEGER (KIND=i4)                                :: i
    INTEGER (KIND=i2), DIMENSION (nxtrack)           :: tmp_flg
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    land_water_flg = 0 
    glint_flg = 0 
    snow_ice_flg = 0

    ! -----------------------------------------------
    ! Save input variable in TMP_FLG for modification
    ! -----------------------------------------------
    tmp_flg(1:nxtrack) = int(omi_geoflg(1:nxtrack), kind=2)  ;  tmp_bytes = 0
    ! FIXME - TEMPO ground_pixel_flag is now int4, but all subroutines
    ! in this module assume int2

    CALL convert_2bytes_to_16bits ( &
         nbyte, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-1) )

    ! ------------------------------
    ! The Glint flag is easy: Byte 4
    ! ------------------------------
    glint_flg(1:nxtrack) = tmp_bytes(1:nxtrack,4)

    ! ------------------------------------------------------------------
    ! Land/Water and Ice require a bit more work. The BIT slices must be
    ! multiplied with the corresponding powers of 2. The sum over this
    ! product is the information we seek.
    ! ------------------------------------------------------------------
    DO i = 1, nxtrack
      land_water_flg(i) = SUM(tmp_bytes(i,0:3 )*seven_byte(1:4))
      snow_ice_flg  (i) = SUM(tmp_bytes(i,8:14)*seven_byte(1:7))
    END DO

    RETURN
  END SUBROUTINE convert_gpqualflag_info

  SUBROUTINE convert_2bytes_to_16bits ( nbits, ndim, byte_num, bit_num )

    ! ==========================================================
    ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
    ! converts it into an NDIM x 16 dimensional interger BIT_NUM
    ! ==========================================================

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                   INTENT (IN) :: ndim, nbits
    INTEGER (KIND=i2), DIMENSION (ndim), INTENT (IN) :: byte_num

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i2), DIMENSION (ndim,0:nbits-1), INTENT (OUT) :: bit_num

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                   :: i!, j
    REAL (KIND=dp)                      :: powval
    INTEGER (KIND=i2), DIMENSION (ndim) :: tmp_byte

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    bit_num = 0_i2

    ! ------------------------------------------------
    ! Save input variable in TMP_BYTE for modification
    ! ------------------------------------------------
    tmp_byte(1:ndim) = byte_num(1:ndim)

    WHERE ( tmp_byte(1:ndim) < 0)
      tmp_byte(1:ndim) = int(tmp_byte(1:ndim) + 65536 , kind=i2)
    ENDWHERE

    ! -------------------------------------------------------------------
    ! Starting with the highest power NBYTES-1, subtract powers of 2 and
    ! assign 1 whereever the power fits in the flag number. At the end we
    ! arrive at a 16 BIT binary representation from which we can extract
    ! the surface information.
    ! -------------------------------------------------------------------
    DO i = nbits-1, 0, -1
      powval = 2.0 ** i
      IF ( powval > 0 ) THEN
        WHERE ( tmp_byte(1:ndim) >= powval )
          bit_num(1:ndim,i) = 1_i2
          tmp_byte(1:ndim) = tmp_byte(1:ndim) - int(powval, kind=i2)
        ENDWHERE
      END IF
    END DO

    RETURN
  END SUBROUTINE convert_2bytes_to_16bits

  SUBROUTINE convert_16bits_to_2bytes ( nbits, ndim, bit_num, byte_num )

    ! ==========================================================
    ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
    ! converts it into an NDIM x 16 dimensional interger BIT_NUM
    ! ==========================================================

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ----------------------
    ! Input/Out variables
    ! ----------------------
    INTEGER (KIND=i4),                   INTENT (IN)           :: ndim, nbits
    INTEGER (KIND=i2), DIMENSION (ndim,0:nbits-1), INTENT (IN) :: bit_num
    INTEGER (KIND=i2), DIMENSION (ndim), INTENT (OUT)          :: byte_num

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: i!, j

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    byte_num = 0

    DO i = 0, nbits - 1 
      byte_num = int(byte_num + bit_num(:, i) * 2 ** i , kind=i2)
    ENDDO


    WHERE (byte_num > 32767)
      byte_num = int(byte_num - 65536 , kind=i2)
    ENDWHERE


    RETURN
  END SUBROUTINE convert_16bits_to_2bytes


  SUBROUTINE coadd_2bytes_qflgs(nbits, ndim, qflg1)!, qflg2)

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! --------------------------
    ! Input/Output variables
    ! --------------------------
    INTEGER (KIND=i4),                    INTENT (IN)  :: nbits, ndim
    INTEGER (KIND=i2), DIMENSION(ndim), INTENT (INOUT) :: qflg1
    !INTEGER (KIND=i2), DIMENSION(ndim), INTENT (IN)    :: qflg2

    INTEGER (KIND=i2), DIMENSION (ndim, 0:nbits-1) :: bit_num1, bit_num2

    CALL convert_2bytes_to_16bits ( nbits, ndim, qflg1, bit_num1 )
    CALL convert_2bytes_to_16bits ( nbits, ndim, qflg1, bit_num2 )
    bit_num1 = bit_num1 + bit_num2

    WHERE(bit_num1 > 1) 
      bit_num1 = 1
    ENDWHERE

    ! Convert 16 bits to 2bytes
    CALL convert_16bits_to_2bytes (nbits, ndim, bit_num1, qflg1)

    RETURN

  END SUBROUTINE coadd_2bytes_qflgs


  ! xliu, 03/25/2011, add several subroutines for bit-based operation of 8-bit unsigned integer
  SUBROUTINE convert_byte_to_8bits ( nbits, ndim, byte_num, bit_num )

    ! ==========================================================
    ! Takes an NDIM dimensional unsigned Byte integer BYTE_NUM and
    ! converts it into an NDIM x 8 dimensional interger BIT_NUM
    ! ==========================================================

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                   INTENT (IN) :: ndim, nbits
    INTEGER (KIND=i1), DIMENSION (ndim), INTENT (IN) :: byte_num

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i1), DIMENSION (ndim,0:nbits-1), INTENT (OUT) :: bit_num

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                   :: i!, j
    REAL (KIND=dp)                      :: powval
    INTEGER (KIND=i2), DIMENSION (ndim) :: tmp_byte

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    bit_num = 0_i1

    ! ------------------------------------------------
    ! Save input variable in TMP_BYTE for modification
    ! ------------------------------------------------
    ! Note: fortran does not have unsigned integer, so copy it to 2 bytes
    ! and 256 for negative nunmbers
    tmp_byte(1:ndim) = byte_num(1:ndim)
    WHERE ( tmp_byte(1:ndim) < 0 )
      tmp_byte(1:ndim) = int(tmp_byte(1:ndim) + 256 , kind=i2)
    ENDWHERE

    ! -------------------------------------------------------------------
    ! Starting with the highest power NBYTES-1, subtract powers of 2 and
    ! assign 1 whereever the power fits in the flag number. At the end we
    ! arrive at a 8 BIT binary representation from which we can extract
    ! the surface information.
    ! -------------------------------------------------------------------
    DO i = nbits-1, 0, -1
      powval = 2.0 ** i
      IF ( powval > 0 ) THEN
        WHERE ( tmp_byte(1:ndim) >= powval )
          bit_num(1:ndim,i) = 1_i1
          tmp_byte(1:ndim) = tmp_byte(1:ndim) - int(powval, kind=i2)
        ENDWHERE
      END IF
    END DO

    RETURN
  END SUBROUTINE convert_byte_to_8bits

  SUBROUTINE convert_8bits_to_byte ( nbits, ndim, bit_num, byte_num )

    ! ==========================================================
    ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
    ! converts it into an NDIM x 16 dimensional interger BIT_NUM
    ! ==========================================================

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ----------------------
    ! Input/Out variables
    ! ----------------------
    INTEGER (KIND=i4),                   INTENT (IN)           :: ndim, nbits
    INTEGER (KIND=i1), DIMENSION (ndim,0:nbits-1), INTENT (IN) :: bit_num
    INTEGER (KIND=i1), DIMENSION (ndim), INTENT (OUT)          :: byte_num

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: i

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    byte_num = 0

    DO i = 0, nbits - 1 
      byte_num = int(byte_num + bit_num(:, i) * 2 ** i , kind=i1)
    ENDDO

    WHERE (byte_num > 127)
      byte_num = int(byte_num - 256 , kind=i1)
    ENDWHERE

    RETURN
  END SUBROUTINE convert_8bits_to_byte



  SUBROUTINE coadd_byte_qflgs(nbits, ndim, qflg1)!, qflg2)

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! --------------------------
    ! Input/Output variables
    ! --------------------------
    INTEGER (KIND=i4),                    INTENT (IN)  :: nbits, ndim
    INTEGER (KIND=i1), DIMENSION(ndim), INTENT (INOUT) :: qflg1
    !INTEGER (KIND=i1), DIMENSION(ndim), INTENT (IN)    :: qflg2

    INTEGER (KIND=i1), DIMENSION (ndim, 0:nbits-1) :: bit_num1, bit_num2

    CALL convert_byte_to_8bits ( nbits, ndim, qflg1, bit_num1 )
    CALL convert_byte_to_8bits ( nbits, ndim, qflg1, bit_num2 )
    bit_num1 = bit_num1 + bit_num2

    WHERE(bit_num1 > 1) 
      bit_num1 = 1
    ENDWHERE

    ! Convert 8 bits to 1 byte
    CALL convert_8bits_to_byte (nbits, ndim, bit_num1, qflg1)

    RETURN

  END SUBROUTINE coadd_byte_qflgs


  SUBROUTINE convert_xtrackqflag_info ( nxtrack, omi_xtrackqflg, &
       rowanomaly_flg, waveshift_flg, blockage_flg, straysun_flg, strayearth_flg )

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN) :: nxtrack
    INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (IN) :: omi_xtrackqflg

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i1), DIMENSION (nxtrack), INTENT (OUT) :: rowanomaly_flg, &
         waveshift_flg, blockage_flg, straysun_flg, strayearth_flg

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4),                PARAMETER      :: nbyte = 8
    INTEGER (KIND=i2), DIMENSION (7), PARAMETER      :: seven_byte = int((/ 1, 2, 4, 8, 16, 32, 64 /), kind=i2)
    INTEGER (KIND=i4)                                :: i
    INTEGER (KIND=i1), DIMENSION (nxtrack)           :: tmp_flg
    INTEGER (KIND=i1), DIMENSION (nxtrack,0:nbyte-1) :: tmp_bytes

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    rowanomaly_flg = 0; waveshift_flg = 0; blockage_flg = 0
    straysun_flg = 0; strayearth_flg = 0

    ! -----------------------------------------------
    ! Save input variable in TMP_FLG for modification
    ! -----------------------------------------------
    tmp_flg(1:nxtrack) = omi_xtrackqflg(1:nxtrack)  ;  tmp_bytes = 0

    CALL convert_byte_to_8bits (nbyte, nxtrack, tmp_flg(1:nxtrack), tmp_bytes(1:nxtrack,0:nbyte-1))

    waveshift_flg(1:nxtrack)  = tmp_bytes(1:nxtrack,4)
    blockage_flg(1:nxtrack)   = tmp_bytes(1:nxtrack,5)
    straysun_flg(1:nxtrack)   = tmp_bytes(1:nxtrack,6)
    strayearth_flg(1:nxtrack) = tmp_bytes(1:nxtrack,7)

    ! ------------------------------------------------------------------
    ! Row anomaly require a bit more work. The BIT slices must be
    ! multiplied with the corresponding powers of 2. The sum over this
    ! product is the information we seek.
    ! ------------------------------------------------------------------
    DO i = 1, nxtrack
      rowanomaly_flg(i) = int(SUM(tmp_bytes(i,0:2)*seven_byte(1:3)), kind=i1)
    ENDDO

    RETURN
  END SUBROUTINE convert_xtrackqflag_info

end module m_convert_coadd
