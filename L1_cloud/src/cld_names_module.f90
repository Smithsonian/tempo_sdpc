!>Inclusion of parameter name lists for netCDF files read/written by L1_cloud
!
!> @param cld_max_name_len maximum string length for names in include files
!
module cld_names_module
  implicit none
  integer, public, parameter :: cld_max_name_len = 64
  include 'cld_names_grp.inc'
  include 'cld_names_dim.inc'
  include 'cld_names_var.inc'
end module cld_names_module
