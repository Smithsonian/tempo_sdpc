MODULE omi_parameters_module

  USE OMSAO_precision_module, ONLY: i4
  ! ---------------------------------
  ! Maximum OMI data/swath dimensions
  ! ---------------------------------
  INTEGER (KIND=i4), PARAMETER :: &
    omi_nxtrack_max    =  60, &
    omi_nwavel_max     = 1024, &
    omi_nwavelcoef_max =   5

END MODULE
