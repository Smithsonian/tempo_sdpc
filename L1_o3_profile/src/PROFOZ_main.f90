PROGRAM PROFOZ_main

  USE PGS_PC_class, only : pgs_pc_getconfigdata
  USE OMSAO_indices_module, ONLY: omi_idx, tempo_idx, instrument_idx, gome2_idx
  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  USE OMSAO_variables_module, ONLY:l2_filename
  USE OMI_LUN_set
  USE omi_pge_process
  USE tmpo_pge_process
  !USE gome2_pge_process
  USE m_read_fitting_controls, only: read_fitting_control_file
  USE m_read_reference_spectra

  IMPLICIT NONE

  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'PROFOZ_main'

  ! -------------------------
  ! OMI L1B related variables
  ! -------------------------
  INTEGER :: errstat, OMI_SMF_setmsg, status

  integer, parameter :: versionid_lun = 123456
  integer :: processing_version = 1 ! default value, over-ridden by VERSIONID from PCF file

  ! ---------------------------------------------------------------------------
  ! Some variables/parameters that are specific to the GOME way of doing things
  ! ---------------------------------------------------------------------------
  INTEGER, PARAMETER       :: fcunit   = FIT_CTRL_LUN, &
                              specunit = OZPROF_CTRL_LUN

  ! ------------------------------------------------------------------
  ! The general PGE error status variable. This is a relatively late
  ! addition and is not used consistently in all routines yet. However,
  ! it is envisaged that it will be used ubiquitously througout the
  ! PGE once the PGE developer gets around to implementing it as such.
  ! ------------------------------------------------------------------
  INTEGER :: pge_error_status
  REAL (kind=dp) :: e1, e2, runtime
  CHARACTER (100) :: message

  pge_error_status = pge_errstat_ok

  call cpu_time(e2)
  !----------------------------------------------------------------------------
  ! read fitting control / ozprof fitting control
  !----------------------------------------------------------------------------
  CALL read_fitting_control_file (fcunit, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) THEN 
       WRITE(*,*) www_message
       GOTO 666
  ENDIF
  !----------------------------------------------------------------------------
  ! Read reference spectra
  !----------------------------------------------------------------------------
  CALL read_reference_spectra ( specunit, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) THEN 
       WRITE(*,*) "Errors in read_reference_spectra"
       GOTO 666
  ENDIF

  status = pgs_pc_getconfigdata (versionid_lun, message)
  if (status == 0) then
    read (message, *)processing_version
  endif

  ! ---------------------------------------------------------------------------
  ! Algorithm running
  ! ---------------------------------------------------------------------------
  SELECT CASE ( instrument_idx )
  CASE ( omi_idx )
    CALL omi_fitting_process  (message, pge_error_status)
  CASE ( gome2_idx)
    !CALL gome2_fitting_process  (message, pge_error_status)
  CASE (tempo_idx)
    CALL tmpo_fitting_process (message, processing_version, pge_error_status)
  CASE DEFAULT
    pge_error_status = pge_errstat_error
  END SELECT
  IF (pge_error_status /= pge_errstat_ok) THEN 
    WRITE(*,'(A)') message
    STOP 1
  ENDIF

  WRITE(*,'(A)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  call cpu_time(e1) 
  runtime = e1-e2
  WRITE(*,'(A,f8.2)') '@ Finish:'//ADJUSTL(TRIM(l2_filename))//"time:",runtime

 666 continue
  write(*,'(a,i9)')' pge_error_status = ',pge_error_status
  ! ------------------------------------
  ! Write END_OF_RUN message to log file
  ! ------------------------------------
  SELECT CASE ( pge_error_status )
  CASE ( pge_errstat_ok )
    ! ----------------------------------------------------------------
    ! PGE execution completed successfully. All is well.
    ! ----------------------------------------------------------------
    errstat = OMI_SMF_setmsg ( OMSAO_S_ENDOFRUN, '', modulename, 0 )
    STOP 0
  CASE ( pge_errstat_warning )
    ! ----------------------------------------------------------------
    ! PGE execution raised non-terminal warnings. Nothing serious, we
    ! hope, so execution completed but with a non-zero exit status.
    ! ----------------------------------------------------------------
    errstat = OMI_SMF_setmsg ( OMSAO_W_ENDOFRUN, '', modulename, 0 )
    STOP 1
  CASE ( pge_errstat_error )

    ! ----------------------------------------------------------------
    ! PGE execution encountered an error that lead to termination.
    ! ----------------------------------------------------------------
    errstat = OMI_SMF_setmsg ( OMSAO_E_ENDOFRUN, '', modulename, 0 )
    STOP 1
  CASE DEFAULT
    ! ----------------------------------------------------------------
    ! If we ever reach here, then PGE_ERRSTAT has been set to a funny
    ! value. This should never happen, but we buffer this case anyway.
    ! ----------------------------------------------------------------
    errstat = OMI_SMF_setmsg ( OMSAO_U_ENDOFRUN, '', modulename, 0 )
    STOP 1
  END SELECT
END PROGRAM PROFOZ_main

