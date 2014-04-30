PROGRAM main
  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  INTEGER (KIND=i4) :: exit_value

  CALL OMSAO_main ( exit_value )

  CALL EXIT (exit_value)
  STOP
END PROGRAM main
