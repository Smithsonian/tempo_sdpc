PROGRAM main
  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE
  INTEGER (KIND=i4) :: exit_value
  exit_value = -1
  CALL OMSAO_main ( exit_value )
  call c_exit (exit_value)
END PROGRAM main
