SUBROUTINE OMSAO_main ( exit_value )

  ! ************************************************************************
  !
  ! This is the main program of the SAO Product Generation Executives (PGEs)
  ! for  the Ozone Monitoring Instrument (OMI).
  !
  ! Authors: Thomas P. Kurosu, Kelly Chance
  !          Smithsonian Astrophysical Observatory
  !          60 Garden Street (MS 50)
  !          Cambridge, MA 02138 (USA)
  !
  !          EMail: tkurosu@cfa.harvard.edu
  !                 kchance@cfa.harvard.edu
  !
  ! ************************************************************************

  USE OMSAO_precision_module
  !USE OMSAO_indices_module,    ONLY: pge_bro_idx
  USE OMSAO_parameters_module
  !USE OMSAO_errstat_module
  use tell_module
  USE metadata_tools, ONLY: init_metadata
  USE omi_pge_fitting_aux, ONLY: set_input_pointer_and_versions
  USE pcf_file_module, ONLY: read_pcf_file
  USE read_reference_spectra, ONLY: read_ref_spectra
  USE omi_pge_fitting_process, ONLY: omi_pge_fitting
  use optimizer_interface_module, only : optimizer_set_default_method
  use elsunc_interface_module, only : elsunc_optimizer
  use slitfunction, only : slitfunction_select, slitfunction_open
  use slitfunction_omi, only : omi_slitfunc_read, omi_slitfunc_convolve
  use ctrlvars, only: yn_use_labslitfunc, yn_do_he5_output, yn_wrt_odl
  use OMSAO_omidata_module, only : initialize_omidata_structs, &
       deallocate_omidata_structs
  use OMSAO_variables_module, only : allocate_refspec_storage, &
    allocate_common_mode_storage, common_mode_spec, &
    deallocate_refspec_common_mode
  use irradiance_data, only: Irr_Data, deallocate_irr_data_type
  IMPLICIT NONE

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (inout) :: exit_value

  ! -------------------------
  ! Name of module/subroutine
  ! -------------------------
  !CHARACTER (LEN=10), PARAMETER :: modulename = 'OMSAO_main'

  ! ------------------------------------------------------------------
  ! The general PGE error status variable. This is a relatively late
  ! addition and is not used consistently in all routines yet. However,
  ! it is envisaged that it will be used ubiquitously througout the
  ! PGE once the PGE developer gets around to implementing it as such.
  ! ------------------------------------------------------------------
  INTEGER   (KIND=i4)      :: errstat, pge_idx !pge_error_status
  CHARACTER (LEN=12) :: pge_name
  !integer :: env_variable_status

  ! --------------------------------------------------------------------
  ! Maximum number of points in any reference spectrum. This is used for
  ! automatic (i.e., SUBROUTINE argument) memory allocation and hence is
  ! defined here instead of in a MODULE.
  ! --------------------------------------------------------------------
  INTEGER (KIND=i4) :: n_max_rspec

  !exit_value = -1   ! early return will indicate an error has occured
                     ! already set in main.f90 and passed to this subroutine

  ! command-line arguments
  integer(kind=4) :: iarg=0
  character (len=255) :: arg

  ! ----------------------------
  ! Set PGE_ERROR_STATUS to O.K.
  ! ----------------------------
  !pge_error_status = pge_errstat_ok
  errstat = 0

  call tell_set_log_level (1)

  ! Control he5 output
  ! Default for other codes is to not write he5 unless switch is set
  yn_do_he5_output = .false.
  ! Whether to write ODL metadata as text file alongside netCDF
  yn_wrt_odl = .false.
  do
    call get_command_argument (iarg, arg)
    if (len(trim(arg)) == 0) exit
    if (trim(arg) == "+he5_out") then
      yn_do_he5_output = .true.
    else if (trim(arg) == "-wrt_odl") then
      yn_wrt_odl = .true.
    else if (trim(arg) == "-h") then
      print *, "Usage: L1_trace_gas [options]"
      print *, ""
      print *, "   +he5_out    enable he5 output"
      print *, "   -wrt_odl    write ODL metadata with netCDF"
      return
    endif
    iarg = iarg + 1
  enddo

  ! Do he5 output unless this environment variable is set.
!  call get_environment_variable ("TG_NO_HE5_OUTPUT", status=env_variable_status)
!  yn_do_he5_output = (env_variable_status == 1 .or. env_variable_status == 2)

  ! ----------------------------------------------------------------------------
  CALL unbufferSTDout()                       ! Make PGE write STD/IO unbuffered
  ! ----------------------------------------------------------------------------
  call maybe_setenv_msgenv (errstat)
  if (errstat /= 0) return

  call allocate_refspec_storage (errstat)
  if (errstat /= 0) return

  ! allocate common mode here in case we read in common mode (see read_ref_spectra)
  call allocate_common_mode_storage (common_mode_spec, errstat)
  if (errstat /= 0) return

  call initialize_omidata_structs (errstat)
  if (errstat /= 0) return

  call optimizer_set_default_method (elsunc_optimizer)
  call slitfunction_select (omi_slitfunc_read, omi_slitfunc_convolve)

  !errstat = pge_errstat_ok
  CALL read_pcf_file (pge_idx, pge_name, errstat )
  if (errstat /= 0) return
  !CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
  !  modulename//f_sep//"READ_PCF_FILE.", vb_lev_default, pge_error_status )
  !IF ( pge_error_status >= pge_errstat_error ) GOTO 666

  CALL init_metadata (errstat)  ! Initialize MetaData
  if (errstat /= 0) return
  !CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
  !  modulename//f_sep//"INIT_METADATA.", vb_lev_default, pge_error_status )
  !IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  CALL read_ref_spectra ( pge_idx, n_max_rspec, errstat )     ! Read reference spectra
  if (errstat /= 0) return
  !CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
  !  modulename//f_sep//"READ_REFERENCE_SPECTRA.", vb_lev_default, pge_error_status )
  !IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  ! ------------------------------------------------------------------------------------
  !CALL omi_slitfunc_read ( errstat )                   ! Read OMI slit function
  ! ------------------------------------------------------------------------------------
  !CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
  !  modulename//f_sep//"OMI_SLITFUNC_READ.", vb_lev_default, pge_error_status )
  !IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666
  call slitfunction_open (errstat, use_table=yn_use_labslitfunc)
  if (errstat /= 0) return

  ! ---------------------------------------------
  ! Set number of InputPointers and InputVersions
  ! ---------------------------------------------
  CALL set_input_pointer_and_versions ( pge_idx )

  CALL omi_pge_fitting  ( pge_idx, n_max_rspec, errstat )   ! Where all the work is done
  if (errstat /= 0) return
  !CALL error_check ( errstat, pge_errstat_warning, errstat, OMSAO_A_SUBROUTINE, &
  !  modulename//f_sep//"OMI_PGE_FITTING_PROCESS.", vb_lev_default, pge_error_status )
  !IF ( pge_error_status >= pge_errstat_fatal ) GOTO 666

  !-------------------------------------
  call deallocate_omidata_structs (errstat)
  if (errstat /= 0) return

  call deallocate_refspec_common_mode (common_mode_spec, errstat)
  if (errstat /= 0) return

  call deallocate_irr_data_type (Irr_Data, errstat)
  if (errstat /= 0) return
  ! ------------------------------------
  ! Write END_OF_RUN message to log file
  ! ------------------------------------
 !666 CALL pge_error_status_exit ( pge_error_status, exit_value )
  exit_value = 0

  RETURN
END SUBROUTINE OMSAO_main
