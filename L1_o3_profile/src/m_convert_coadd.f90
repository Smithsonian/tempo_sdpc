!
module m_convert_coadd

  public convert_2bytes_to_16bits,convert_2byte_to_32bits,  &
         convert_byte_to_8bits,coadd_2bytes_qflgs, coadd_byte_qflgs,  &
         prespec_align, solwavcal_coadd, radwavcal_coadd, stray_coadd


  private convert_16bits_to_2bytes,  &
          convert_8bits_to_byte
           !prespec_align1 
contains



  SUBROUTINE convert_2bytes_to_32bits ( nbits, ndim, byte_num, bit_num)
    ! ==========================================================
    ! Takes an NDIM dimensional 2Byte integer BYTE_NUM and
    ! converts it into an NDIM x 32 dimensional interger BIT_NUM
    ! ==========================================================
    ! jbak please check this subroutine 
    USE OMSAO_precision_module
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                   INTENT (IN) :: ndim, nbits
    INTEGER (KIND=i4), DIMENSION (ndim), INTENT (IN) :: byte_num

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), DIMENSION (ndim,0:nbits-1), INTENT (OUT) :: bit_num

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                   :: i!, j
    REAL (KIND=dp)                      :: powval
    INTEGER (KIND=i4), DIMENSION (ndim) :: tmp_byte

    ! ----------------------------
    ! Initialize output quantities
    ! ----------------------------
    bit_num = 0_i4

    ! ------------------------------------------------
    ! Save input variable in TMP_BYTE for modification
    ! ------------------------------------------------
    tmp_byte(1:ndim) = byte_num(1:ndim)

    WHERE ( tmp_byte(1:ndim) < 0)
      tmp_byte(1:ndim) = int(tmp_byte(1:ndim) + 4294967296 , kind=i4)
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
          bit_num(1:ndim,i) = 1_i4
          tmp_byte(1:ndim) = tmp_byte(1:ndim) - int(powval, kind=i4)
        ENDWHERE
      END IF
    END DO

    RETURN
  END SUBROUTINE convert_2bytes_to_32bits


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


  SUBROUTINE coadd_2bytes_qflgs(nbits, ndim, qflg1, qflg2)

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! --------------------------
    ! Input/Output variables
    ! --------------------------
    INTEGER (KIND=i4),                    INTENT (IN)  :: nbits, ndim
    INTEGER (KIND=i2), DIMENSION(ndim), INTENT (INOUT) :: qflg1
    INTEGER (KIND=i2), DIMENSION(ndim), INTENT (IN)    :: qflg2

    INTEGER (KIND=i2), DIMENSION (ndim, 0:nbits-1) :: bit_num1, bit_num2

    CALL convert_2bytes_to_16bits ( nbits, ndim, qflg1, bit_num1 )
    CALL convert_2bytes_to_16bits ( nbits, ndim, qflg2, bit_num2 )
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



  SUBROUTINE coadd_byte_qflgs(nbits, ndim, qflg1, qflg2)

    USE OMSAO_precision_module
    IMPLICIT NONE

    ! --------------------------
    ! Input/Output variables
    ! --------------------------
    INTEGER (KIND=i4),                    INTENT (IN)  :: nbits, ndim
    INTEGER (KIND=i1), DIMENSION(ndim), INTENT (INOUT) :: qflg1
    INTEGER (KIND=i1), DIMENSION(ndim), INTENT (IN)    :: qflg2

    INTEGER (KIND=i1), DIMENSION (ndim, 0:nbits-1) :: bit_num1, bit_num2

    CALL convert_byte_to_8bits ( nbits, ndim, qflg1, bit_num1 )
    CALL convert_byte_to_8bits ( nbits, ndim, qflg2, bit_num2 )
    bit_num1 = bit_num1 + bit_num2

    WHERE(bit_num1 > 1) 
      bit_num1 = 1
    ENDWHERE

    ! Convert 8 bits to 1 byte
    CALL convert_8bits_to_byte (nbits, ndim, bit_num1, qflg1)

    RETURN

  END SUBROUTINE coadd_byte_qflgs

  

  ! Properly align several cross-track spectra (to be coadded) to within
  ! one pixel
  SUBROUTINE prespec_align(nw, nspec, wavl, spec, prec, qflg)

    USE OMSAO_precision_module
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ================================
    ! Input and Output variables
    ! =================================
    INTEGER, INTENT(IN)                                    :: nspec, nw
    REAL (KIND=r4),    DIMENSION(nw, nspec), INTENT(INOUT) :: wavl, spec, prec
    INTEGER (KIND=i2), DIMENSION(nw, nspec), INTENT(INOUT) :: qflg

    ! Local variables
    REAL (KIND=dp)                                         :: dwvl, mx_fwav
    REAL (KIND=dp), DIMENSION(nw)                          :: diff
    INTEGER                                                :: fidx, i, j, nsel
    mx_fwav = MAXVAL(wavl(1, 1:nspec))

    DO i = 1, nspec
      diff(1:nw) = ABS(wavl(1:nw, i) - mx_fwav)
      fidx       = MINVAL(MINLOC(diff(1:nw)))
      nsel = nw - fidx + 1

      ! Shift the spectrum
      wavl(1:nsel, i) = wavl(fidx:nw, i)
      spec(1:nsel, i) = spec(fidx:nw, i); spec(nsel+1:nw, i) = 0.0
      prec(1:nsel, i) = prec(fidx:nw, i)
      qflg(1:nsel, i) = qflg(fidx:nw, i)

      ! Make sure that wavelengths are increasing
      dwvl = wavl(2, i) - wavl(1, i)
      DO j = nsel + 1, nw
        wavl(j, i) = real(wavl(j-1, i) + dwvl , kind=r4)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE prespec_align
  
  SUBROUTINE stray_coadd(nspec, ncoadd, allspec)

    USE OMSAO_precision_module
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ================================
    ! Input and Output variables
    ! =================================
    INTEGER, INTENT(IN)                                        :: nspec, ncoadd
    REAL (KIND=dp), DIMENSION(ncoadd, 2, nspec), INTENT(INOUT) :: allspec

    ! ===============
    ! Local variables
    ! ===============
    INTEGER       :: i

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    !CHARACTER (LEN=11), PARAMETER :: modulename = 'stray_coadd'

    DO i = 2, ncoadd
      allspec(1, :, :)    = allspec(1, :, :) + allspec(i, :, :)
    ENDDO
    allspec(1, 1, :) = allspec(1, 1, :) / ncoadd
    allspec(1, 2, :) = allspec(1, 2, :) / ncoadd

    RETURN

  END SUBROUTINE stray_coadd

  SUBROUTINE solwavcal_coadd (wcal_bef_coadd, nspec, ncoadd, allspec, wshis, wsqus, error)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,     ONLY: max_calfit_idx, shi_idx, squ_idx, &
         wvl_idx, spc_idx, sig_idx, hwe_idx, asy_idx, vgl_idx, hwr_idx, spk_idx
    USE OMSAO_variables_module,   ONLY: n_fitvar_sol, fitvar_sol,   &
         mask_fitvar_sol, fitvar_sol_saved, lo_sunbnd, up_sunbnd, fixslitcal, fitwavs, &
         fitweights, currspec, which_slit, fitvar_sol_init, lo_sunbnd_init, &
         up_sunbnd_init, sol_wav_avg, xbin_decerr, correct_lamda
    USE OMSAO_errstat_module
    use m_cal_fit_one
    use m_ezspline_interpolation, only: interpolation

    IMPLICIT NONE

    ! ================================
    ! Input and Output variables
    ! =================================
    INTEGER, INTENT(IN)  :: nspec, ncoadd
    LOGICAL, INTENT(IN)  :: wcal_bef_coadd
    REAL (KIND=dp), DIMENSION(ncoadd, sig_idx, nspec), INTENT(INOUT) :: allspec
    REAL (KIND=dp), DIMENSION(ncoadd),  INTENT(OUT) :: wshis, wsqus
    LOGICAL,                            INTENT(OUT) :: error

    ! ===============
    ! Local variables
    ! ===============
    INTEGER, PARAMETER :: slit_unit =1000
    REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_varstd
    REAL (KIND=dp)                               :: tmpwave, solar_norm, dwvl
    INTEGER       :: i, solfit_exval, errstat!, j
    LOGICAL, SAVE :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=15), PARAMETER :: modulename = 'solwavcal_coadd'

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    error = .FALSE.
    wshis  = 0.0; wsqus = 0.0

    IF (wcal_bef_coadd) THEN

      ! find the locations of actually used fitting variables
      IF (first) THEN
        fixslitcal = .TRUE.    ; slitcal = .TRUE.
        wrt_to_screen = .FALSE.; wrt_to_file = .FALSE.

        IF (which_slit == 5) THEN
          fitvar_sol_init(hwe_idx:asy_idx) = 0_dp
          lo_sunbnd_init(hwe_idx:asy_idx)  = 0_dp
          up_sunbnd_init(hwe_idx:asy_idx)  = 0_dp

          fitvar_sol_init(vgl_idx:spk_idx) = 0_dp
          lo_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
          up_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
          slitcal = .FALSE.; fixslitcal = .FALSE.
        ENDIF
        lo_sunbnd=lo_sunbnd_init
        up_sunbnd=up_sunbnd_init    
        fitvar_sol_saved = fitvar_sol_init
        n_fitvar_sol = 0
        DO i = 1, max_calfit_idx
          IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
            n_fitvar_sol =  n_fitvar_sol + 1
            mask_fitvar_sol(n_fitvar_sol) = i
          END IF
        ENDDO
        first = .FALSE.
      ENDIF

      solar_norm = SUM(allspec(1, spc_idx, 1:nspec)) / nspec
      DO i = 1, ncoadd
        fitwavs   (1:nspec)  = allspec(i, wvl_idx, 1:nspec)
        currspec  (1:nspec)  = allspec(i, spc_idx, 1:nspec) / solar_norm
        fitweights(1:nspec)  = allspec(i, sig_idx, 1:nspec) / solar_norm
        fitvar_sol = fitvar_sol_saved

        CALL cal_fit_one (nspec, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
             slitcal, slit_unit, tmpwave, tmp_varstd, solfit_exval)

        IF (solfit_exval < 0) THEN
          WRITE(*, *) &
               'solwavcal_coadd: calibration does not converge for pixel: ', i
          error = .TRUE.; wshis(i) = 0.; wsqus(i) = 0.
        ELSE 

          ! Shift and squeeze earthshine spectrum        
          IF (correct_lamda == 1) THEN
            allspec(i, wvl_idx, 1:nspec) = (fitwavs(1:nspec) - fitvar_sol(shi_idx)) / & 
                                         (1.0 + fitvar_sol(squ_idx))
          ELSE
            allspec(i, wvl_idx, 1:nspec) = (fitwavs(1:nspec) - fitvar_sol(shi_idx) + & 
                       sol_wav_avg * fitvar_sol(squ_idx) ) / (1.0 + fitvar_sol(squ_idx))
          ENDIF

          ! Make sure that the correction is less than a pixel
          dwvl = fitwavs(2) - fitwavs(1)
          IF (ANY(ABS(allspec(i, wvl_idx, 1:nspec) - fitwavs(1:nspec)) > dwvl)) &
             THEN
            allspec(i, wvl_idx, 1:nspec) = fitwavs(1:nspec) ! Roll back
            wshis(i) = 0.0; wsqus(i) = 0.0
          ELSE
            wshis(i) = fitvar_sol(shi_idx); wsqus(i) = fitvar_sol(squ_idx)
          ENDIF
        ENDIF
      ENDDO ! end  coadd
    ENDIF

    IF (.NOT. wcal_bef_coadd) THEN
      ! simple coadding due to unsuccessful calibration
      DO i = 2, ncoadd
        allspec(1, :, :)    = allspec(1, :, :) + allspec(i, :, :)
      ENDDO
      allspec(1, wvl_idx, :) = allspec(1, wvl_idx, :) / ncoadd
      allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
      allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd  !/ SQRT(1.0 * ncoadd)  ! reduce S/N
    ELSE
      ! Interpolate all spectra to the wavelength grids of first spectrum
      DO i = 2, ncoadd
        CALL interpolation (nspec, allspec(i, wvl_idx, 1:nspec), &
             allspec(i, spc_idx, 1:nspec), nspec - 2,  &
             allspec(1, wvl_idx, 2:nspec-1), allspec(i, spc_idx, 2:nspec-1), &
             errstat )
        IF ( errstat > pge_errstat_warning ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
          STOP 1
        END IF

        allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) + &
             allspec(i, wvl_idx, 1)
        allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) + &
             allspec(i, wvl_idx, nspec)
        allspec(1, spc_idx, :) = allspec(1, spc_idx, :) + &
             allspec(i, spc_idx, :)
        allspec(1, sig_idx, :) = allspec(1, sig_idx, :) + &
             allspec(i, sig_idx, :)
      ENDDO
      allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) / ncoadd
      allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) / ncoadd
      allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
      allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd  !/SQRT(1.0 * ncoadd)  ! reduce S/N
    ENDIF
    ! Reduce measurement error after coadding
    IF (xbin_decerr) allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / SQRT(1.0 *ncoadd)
    RETURN
  END SUBROUTINE solwavcal_coadd

  SUBROUTINE radwavcal_coadd(wcal_bef_coadd, wavcal, & ! iwin, ix, &
       nspec, ncoadd, allspec, wshis, wsqus, error)

    USE OMSAO_precision_module
    USE OMSAO_indices_module,     ONLY: max_calfit_idx, shi_idx, squ_idx, &
         wvl_idx, spc_idx, sig_idx, hwe_idx, asy_idx, vgl_idx, hwr_idx, spk_idx
    USE OMSAO_variables_module,   ONLY: n_fitvar_sol, fitvar_sol,fitvar_sol_saved, &
         mask_fitvar_sol, lo_sunbnd, up_sunbnd, fixslitcal, fitwavs, &
         fitweights, currspec, which_slit, fitvar_sol_init, lo_sunbnd_init, &
         up_sunbnd_init, sol_wav_avg, xbin_decerr, correct_lamda
    USE OMSAO_errstat_module
    use m_cal_fit_one
    use m_ezspline_interpolation, only: interpolation

    IMPLICIT NONE

    ! ================================
    ! Input and Output variables
    ! =================================
    INTEGER, INTENT(IN)  :: nspec, ncoadd
    LOGICAL, INTENT(IN)  :: wcal_bef_coadd, wavcal
    REAL (KIND=dp), DIMENSION(ncoadd, sig_idx, nspec), INTENT(INOUT) :: allspec
    REAL (KIND=dp), DIMENSION(ncoadd),  INTENT(OUT) :: wshis, wsqus
    LOGICAL,                            INTENT(OUT) :: error

    ! ===============
    ! Local variables
    ! ===============
    INTEGER, PARAMETER :: slit_unit =1000
    REAL (KIND=dp), DIMENSION(max_calfit_idx, 2) :: tmp_varstd
    REAL (KIND=dp)                               :: tmpwave, solar_norm, dwvl
    INTEGER       :: i, solfit_exval, errstat!, j
    LOGICAL, SAVE :: wrt_to_screen, wrt_to_file, slitcal, first = .TRUE.

    ! ------------------------------
    ! Name of this subroutine/module
    ! ------------------------------
    CHARACTER (LEN=15), PARAMETER :: modulename = 'radwavcal_coadd'

    ! ------------------
    ! External functions
    ! ------------------
    INTEGER :: OMI_SMF_setmsg

    error = .FALSE.
    wshis  = 0.0; wsqus = 0.0

    IF (wcal_bef_coadd .AND. wavcal) THEN

      ! find the locations of actually used fitting variables
     ! IF (first) THEN
        fixslitcal = .TRUE.    ; slitcal = .TRUE.
        wrt_to_screen = .FALSE.; wrt_to_file = .FALSE.

        IF (which_slit == 5) THEN
          fitvar_sol_init(hwe_idx:asy_idx) = 0_dp
          lo_sunbnd_init(hwe_idx:asy_idx)  = 0_dp
          up_sunbnd_init(hwe_idx:asy_idx)  = 0_dp

          fitvar_sol_init(vgl_idx:spk_idx) = 0_dp
          lo_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
          up_sunbnd_init(vgl_idx:spk_idx)  = 0_dp
          slitcal = .FALSE.; fixslitcal = .FALSE.
        ENDIF
          fitvar_sol_saved=fitvar_sol_init
          lo_sunbnd=lo_sunbnd_init
          up_sunbnd=up_sunbnd_init
        n_fitvar_sol = 0
        DO i = 1, max_calfit_idx
          IF (lo_sunbnd(i) < up_sunbnd(i) ) THEN
            n_fitvar_sol =  n_fitvar_sol + 1
            mask_fitvar_sol(n_fitvar_sol) = i
          END IF
        ENDDO
        first = .FALSE.
     ! ENDIF

      solar_norm = SUM(allspec(1, spc_idx, 1:nspec)) / nspec
      DO i = 1, ncoadd
        fitwavs   (1:nspec)  = allspec(i, wvl_idx, 1:nspec)
        currspec  (1:nspec)  = allspec(i, spc_idx, 1:nspec) / solar_norm
        fitweights(1:nspec)  = allspec(i, sig_idx, 1:nspec) / solar_norm
        fitvar_sol = fitvar_sol_saved

        CALL cal_fit_one (nspec, n_fitvar_sol, wrt_to_screen, wrt_to_file,&
             slitcal, slit_unit, tmpwave, tmp_varstd, solfit_exval)

        IF (solfit_exval < 0) THEN
          WRITE(*, *) &
               'solwavcal_coadd: calibration does not converge for pixel: ', i
          error = .TRUE.; wshis(i) = 0.; wsqus(i) = 0.
        ELSE 

          ! Shift and squeeze earthshine spectrum        
          IF (correct_lamda == 1) THEN
            allspec(i, wvl_idx, 1:nspec) = (fitwavs(1:nspec) - fitvar_sol(shi_idx)) / & 
                                         (1.0 + fitvar_sol(squ_idx))
          ELSE
            allspec(i, wvl_idx, 1:nspec) = (fitwavs(1:nspec) - fitvar_sol(shi_idx) + & 
                       sol_wav_avg * fitvar_sol(squ_idx) ) / (1.0 + fitvar_sol(squ_idx))
          ENDIF

          ! Make sure that the correction is less than a pixel
          dwvl = fitwavs(2) - fitwavs(1)
          IF (ANY(ABS(allspec(i, wvl_idx, 1:nspec) - fitwavs(1:nspec)) > dwvl)) &
             THEN
            allspec(i, wvl_idx, 1:nspec) = fitwavs(1:nspec) ! Roll back
            wshis(i) = 0.0; wsqus(i) = 0.0
          ELSE
            wshis(i) = fitvar_sol(shi_idx); wsqus(i) = fitvar_sol(squ_idx)
          ENDIF
        ENDIF
      ENDDO ! end  coadd
   ELSE IF (wcal_bef_coadd .AND. .NOT. wavcal) THEN
      DO i = 1, ncoadd
        IF (correct_lamda == 1) THEN
        allspec(i, wvl_idx, 1:nspec) = (allspec(i, wvl_idx, 1:nspec) - wshis(i)) / (1.0 + wsqus(i))
        ELSE
        allspec(i, wvl_idx, 1:nspec) = (allspec(i, wvl_idx, 1:nspec) - wshis(i) + sol_wav_avg * wsqus(i) ) / (1.0 + wsqus(i))
        ENDIF
      ENDDO
    ENDIF

    IF (.NOT. wcal_bef_coadd) THEN
      ! simple coadding due to unsuccessful calibration
      DO i = 2, ncoadd
        allspec(1, :, :)    = allspec(1, :, :) + allspec(i, :, :)
      ENDDO
      allspec(1, wvl_idx, :) = allspec(1, wvl_idx, :) / ncoadd
      allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
      allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd  !/ SQRT(1.0 * ncoadd)  ! reduce S/N
    ELSE
      ! Interpolate all spectra to the wavelength grids of first spectrum
      DO i = 2, ncoadd
        CALL interpolation (nspec, allspec(i, wvl_idx, 1:nspec), &
             allspec(i, spc_idx, 1:nspec), nspec - 2,  &
             allspec(1, wvl_idx, 2:nspec-1), allspec(i, spc_idx, 2:nspec-1), &
             errstat )
        IF ( errstat > pge_errstat_warning ) THEN
          errstat = OMI_SMF_setmsg (omsao_e_interpol, modulename, '', 0)
          STOP 1
        END IF

        allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) + &
             allspec(i, wvl_idx, 1)
        allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) + &
             allspec(i, wvl_idx, nspec)
        allspec(1, spc_idx, :) = allspec(1, spc_idx, :) + &
             allspec(i, spc_idx, :)
        allspec(1, sig_idx, :) = allspec(1, sig_idx, :) + &
             allspec(i, sig_idx, :)
      ENDDO
      allspec(1, wvl_idx, 1) = allspec(1, wvl_idx, 1) / ncoadd
      allspec(1, wvl_idx, nspec) = allspec(1, wvl_idx, nspec) / ncoadd
      allspec(1, spc_idx, :) = allspec(1, spc_idx, :) / ncoadd  ! reduce S/N
      allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / ncoadd  !/SQRT(1.0 * ncoadd)  ! reduce S/N
    ENDIF
    ! Reduce measurement error after coadding
    IF (xbin_decerr) allspec(1, sig_idx, :) = allspec(1, sig_idx, :) / SQRT(1.0 *ncoadd)
    RETURN

  END SUBROUTINE radwavcal_coadd


  !  Unused?
  !
  ! Properly align several cross-track spectra (to be coadded) to within
  ! one pixel
!  SUBROUTINE prespec_align1(nw, nspec, wavl, spec, spec1, spec2, prec, qflg)
!
!    USE OMSAO_precision_module
!    USE OMSAO_errstat_module
!
!    IMPLICIT NONE
!
!    ! ================================
!    ! Input and Output variables
!    ! =================================
!    INTEGER, INTENT(IN)                                    :: nspec, nw
!    REAL (KIND=r4),    DIMENSION(nw, nspec), INTENT(INOUT) :: wavl, spec, &
!         prec, spec1, spec2
!    INTEGER (KIND=i2), DIMENSION(nw, nspec), INTENT(INOUT) :: qflg
!
!    ! Local variables
!    REAL (KIND=dp)                                         :: dwvl, mx_fwav
!    REAL (KIND=dp), DIMENSION(nw)                          :: diff
!    INTEGER                                               :: fidx, i, j, nsel
!
!    mx_fwav = MAXVAL(wavl(1, 1:nspec))
!
!    DO i = 1, nspec
!      diff(1:nw) = ABS(wavl(1:nw, i) - mx_fwav)
!      fidx       = MINVAL(MINLOC(diff(1:nw)))
!      nsel = nw - fidx + 1
!
!      ! Shift the spectrum
!      wavl(1:nsel, i)  = wavl(fidx:nw, i)
!      spec(1:nsel, i)  = spec(fidx:nw, i);  spec(nsel+1:nw, i)  = 0.0
!      spec1(1:nsel, i) = spec1(fidx:nw, i); spec1(nsel+1:nw, i) = 0.0
!      spec2(1:nsel, i) = spec2(fidx:nw, i); spec2(nsel+1:nw, i) = 0.0
!      prec(1:nsel, i)  = prec(fidx:nw, i)
!      qflg(1:nsel, i)  = qflg(fidx:nw, i)
!
!      ! Make sure that wavelengths are increasing
!      dwvl = wavl(2, i) - wavl(1, i)
!      DO j = nsel + 1, nw
!        wavl(j, i) = real(wavl(j-1, i) + dwvl , kind=r4)
!      ENDDO
!    ENDDO
!
!    RETURN
!  END SUBROUTINE prespec_align1
end module m_convert_coadd
