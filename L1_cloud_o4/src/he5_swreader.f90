module he5_swreader

integer(kind=4),external::he5_swopen,he5_swattach, &
                          he5_swrdfld,he5_swdetach, &
                          he5_swclose,he5_swcreate, &
                          he5_swdefdim,he5_swdefcomch, &
                          he5_swdefdfld,he5_swwrfld, &
                          he5_swrdattr,he5_swrdgattr, &
                          he5_swrdlattr,he5_swinqgflds, &
                          he5_swinqdflds,he5_swinqattrs,&
                          he5_swinqfldalias,he5_swinqgfldalias, &
                          he5_swinqgattrs,he5_inqlattrs, &
                          he5_swlattrinfo,he5_swattrinfo, &
                          he5_swfldinfo,he5_swinqswath, &
                          he5_swdefboxreg,he5_swextreg, &
                          he5_swdeftmeper,he5_swinqlattrs, &
                          he5_swdiminfo,he5_swidtype, &
                          he5_ehrdglatt,he5_ehglattinf, &
                          he5_ehinqglatts,he5_swwrgattr, &
                          he5_swdefgfld

! HDFEOS5 FILE ACCESS TAGS 
! ========================
  integer HE5F_ACC_RDWR
  integer HE5F_ACC_RDONLY
  integer HE5F_ACC_TRUNC
  parameter(HE5F_ACC_RDWR   = 100)
  parameter(HE5F_ACC_RDONLY = 101)
  parameter(HE5F_ACC_TRUNC  = 102)

end module he5_swreader
