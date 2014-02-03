module ctrlvars

  ! The variables in this file are set via the control file and are not
  ! changed while the program is running.  Hence, they are effectively
  ! constants.

  logical, public :: yn_radiance_reference
  logical, public :: yn_common_iter
  logical, public :: yn_spectrum_norm

  ! Logical for newshift following Xiong comments -- gga
  logical, public :: yn_newshift

  logical, public :: yn_smooth
  logical, public :: yn_doas
  logical, public :: yn_use_labslitfunc
  logical, public :: yn_solar_i0
  logical, public :: yn_diagnostic_run
  logical, public :: yn_solar_comp  ! use solar composite spectrum
  logical, public :: yn_o3amf_cor
  logical, public :: yn_solmonthave
  logical, public :: yn_refseccor   !reference sector correction

  ! Scattering weights, gas profile, averaging kernels
  ! and albedo attributes. gga
  logical, public :: yn_scat_weights

  logical, public :: yn_remove_target

end module
