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

  logical, public :: yn_I0


  logical, public :: yn_o3amf_cor
  logical, public :: yn_refseccor   !reference sector correction

  ! Scattering weights, gas profile, averaging kernels
  ! and albedo attributes. gga
  logical, public :: yn_scat_weights

  ! Stratospheric and tropospheric AMF calculation
  logical, public :: yn_stratrop

  logical, public :: yn_remove_target

  ! TEMPO/OMI data flag, set by command line switch, true=OMI, false=TEMPO
  logical, public :: yn_omi_data

  ! Write ODL-format metadata text file alongside netCDF
  logical, public :: yn_wrt_odl

  ! FIXME JCH temporary switches for tempo development
  logical, public :: yn_disable_omi_features, yn_do_he5_output

  !GEMS data flag, set by command line switch, true = GEMS
  logical, public :: yn_gems

  ! Logicals for saving and reading solar irradiance
  ! These are set from PCF inputs
  logical, public :: yn_do_solar_cal
  logical, public :: yn_exit_post_solar_cal
  logical, public :: yn_write_solar_cal
  logical, public :: yn_read_solar_cal

end module
