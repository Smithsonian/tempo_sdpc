module ctrlvars

  ! The variables in this file are set via the control file and are not
  ! changed while the program is running.  Hence, they are effectively
  ! constants.

  logical, public :: yn_radiance_reference

  logical, public :: yn_common_iter
  ! If TRUE, two passes will be made:  The first pass will perform fitting
  !    and produce a common-mode spectrum.  The second pass will perform the
  !    fitting using the derived common mode.
  ! If FALSE, A common mode spectrum will be read from a file and used.


  logical, public :: yn_spectrum_norm

  ! Logical for newshift following Xiong comments -- gga
  logical, public :: yn_newshift

  logical, public :: yn_smooth
  logical, public :: yn_doas
  logical, public :: yn_use_labslitfunc
  logical, public :: yn_solar_i0
  logical, public :: yn_diagnostic_run

  logical, public :: yn_solar_comp
  ! If TRUE, read a composite solar spectrum from a file.
  ! The solar_comp_typ determines the type to read.
  integer, public :: solar_comp_typ


  logical, public :: yn_o3amf_cor
  logical, public :: yn_solmonthave
  logical, public :: yn_refseccor   !reference sector correction

  ! Scattering weights, gas profile, averaging kernels
  ! and albedo attributes. gga
  logical, public :: yn_scat_weights

  ! Stratospheric and tropospheric AMF calculation
  logical, public :: yn_stratrop

  logical, public :: yn_remove_target

  ! FIXME JCH temporary switches for tempo development
  logical, public :: yn_disable_omi_features, yn_do_he5_output

end module
