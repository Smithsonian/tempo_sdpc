PROGRAM PROFOZ_main

  USE OMSAO_precision_module
  USE OMSAO_indices_module,   ONLY: icf_idx
  USE OMSAO_variables_module, ONLY: static_input_fnames,    &
       l1b_rad_filename,                &
       l1b_irrad_filename, l2_filename, &
       l2_cld_filename, & !, pge_idx
       use_he5_in, use_he5_out, use_tio_in, use_tio_out
  USE OMSAO_gome_data_module, ONLY: omi_idx, instrument_idx!, &
  !scia_idx, gome_idx, gome2_idx
  USE OMSAO_errstat_module
  USE OMI_LUN_set
  use m_omi_fitting_process
  use m_read_fitting_controls, only: read_fitting_control_file
  use m_read_reference_spectra

  IMPLICIT NONE

  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'PROFOZ_main'

  ! -------------------------
  ! OMI L1B related variables
  ! -------------------------
  INTEGER :: errstat, OMI_SMF_setmsg

  ! ---------------------------------------------------------------------------
  ! Some variables/parameters that are specific to the GOME way of doing things
  ! ---------------------------------------------------------------------------
  INTEGER, PARAMETER       :: fcunit   = FIT_CTRL_LUN, &
       specunit = OZPROF_CTRL_LUN!, amfunit  = REFDB_DIR_LUN
  !INTEGER, PARAMETER       :: l1funit = 21, l2funit = 22
  CHARACTER (LEN=maxchlen) :: l1_inputs_fname_sol, l1_inputs_fname_rad, &
       l2_output_fname, arg
  INTEGER                  :: l2_hdf_flag, iarg

  ! ------------------------------------------------------------------
  ! The general PGE error status variable. This is a relatively late
  ! addition and is not used consistently in all routines yet. However,
  ! it is envisaged that it will be used ubiquitously througout the
  ! PGE once the PGE developer gets around to implementing it as such.
  ! ------------------------------------------------------------------
  INTEGER :: pge_error_status

  pge_error_status = pge_errstat_ok

  ! Assign fitting input control file
  !! static_input_fnames(icf_idx) = 'INP/SOMIPROF.inp'
  static_input_fnames(icf_idx) = ''

  !! Commented by Kai
  !!CALL unbufferSTDout()            ! Make PGE write STD/IO unbuffered
  !! End Commented by Kai

  !Check for command line args controlling he5/netCDF I/O
  ! default is netCDF input & output
  iarg = 0
  do
    call get_command_argument (iarg, arg)
    if (len_trim(arg) == 0) exit
    if (trim(arg) == "+he5_out") then
      use_he5_out = .true.
      ! currently not possible to produce he5 output
      ! without reading metadata from he5 input
      use_he5_in = .true.
    else if (trim(arg) == "-nc_in") then
      use_tio_in = .false.
      use_he5_in = .true.
    else if (trim(arg) == "-nc_out") then
      use_tio_out = .false.
    else if (trim(arg) == "+he5_in") then
      use_he5_in = .true.
    endif
    iarg = iarg + 1
  enddo


  ! Read fitting conrol parameters from input file
  ! CALL read_fitting_control_file ( pge_idx, pge_error_status )
  CALL read_fitting_control_file (fcunit, static_input_fnames(icf_idx), &
       instrument_idx, l1_inputs_fname_sol, l1_inputs_fname_rad,        &
       l2_output_fname, l2_cld_filename, l2_hdf_flag, pge_error_status )


  IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  ! Read reference spectra
  CALL read_reference_spectra ( specunit, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) GOTO 666


  ! -----------------------------------------------
  ! Different instruments require different actions
  ! -----------------------------------------------
  SELECT CASE ( instrument_idx )

  CASE ( omi_idx )

    l1b_irrad_filename = l1_inputs_fname_sol
    l1b_rad_filename   = l1_inputs_fname_rad
    l2_filename        = l2_output_fname
    CALL omi_fitting_process ( l2_hdf_flag, pge_error_status )
  CASE DEFAULT
    pge_error_status = pge_errstat_error
  END SELECT

  ! ------------------------------------
  ! Write END_OF_RUN message to log file
  ! ------------------------------------
666 SELECT CASE ( pge_error_status )
  CASE ( pge_errstat_ok )
    ! ----------------------------------------------------------------
    ! PGE execution completed successfully. All is well.
    ! ----------------------------------------------------------------
    errstat = OMI_SMF_setmsg ( OMSAO_S_ENDOFRUN, '', modulename, 0 )

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
