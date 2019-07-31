module cache_module
  use OMSAO_precision_module, only :i4, r8
  implicit none

  ! ---------------------------------------------------------------------
  ! The following two quantities are used to determine whether we need to
  ! reconvolve the solar spectrum. Only if either SHIFT or SQUEEZE have
  ! changed from one iteration to the other is a reconvolution necessary.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8) :: saved_shift = -1.0E+30_r8, saved_squeeze = -1.0E+30_r8
  REAL (KIND=r8) :: saved_hwe = -1.0E+30_r8, saved_asy = -1.0E+30_r8, &
       saved_sgk = -1.0E+30_r8

end module cache_module
