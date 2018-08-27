MODULE omi_parameters_module

  USE OMSAO_precision_module, ONLY: i4
  ! ---------------------------------
  ! Maximum OMI data/swath dimensions
  ! ---------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
    omi_nxtrack_max    = 2048, & !60, &   ! JCH 2048 is required for TEMPO
    omi_nwavel_max     = 1026, &  ! JCH synthetic data uses 1026
    omi_nwavelcoef_max =   5

END MODULE
