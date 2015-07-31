PROGRAM main
  USE OMSAO_precision_module, ONLY: i4
  use tell_module
  IMPLICIT NONE
  INTEGER (KIND=i4) :: exit_value
  exit_value = -1
  call tell_open ("L1_trace_gas", 0)
  CALL OMSAO_main ( exit_value )
  call tell_close ()
  call c_exit (exit_value)
END PROGRAM main
