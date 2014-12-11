module HE4_class
      USE HDF4_class
! he4 functions

      INTEGER( KIND = 4 ), EXTERNAL :: swcreate,  &
                                       swdefdim,  &
                                       swdetach,  &
                                       swclose,   &
                                       swopen,    &
                                       swattach,  &
                                       swfldinfo, &
                                       swdiminfo, &
                                       swrdfld,   &
                                       swwrfld,   &
                                       swdefgfld, &
                                       swdefdfld, &
                                       swrdattr,  &
                                       swinqgflds,&
                                       swinqdflds,&
                                       swinqswath,&
                                       swinqattrs,&
                                       swattrinfo

      INTEGER (KIND=4) :: DFKNTSize 
      EXTERNAL            DFKNTSize 

      integer HDFE_NOMERGE
      parameter (HDFE_NOMERGE=0)
      integer HDFE_AUTOMERGE
      parameter (HDFE_AUTOMERGE=1)
      integer HDFE_MIDPOINT
      parameter (HDFE_MIDPOINT=0)
      integer HDFE_ENDPOINT
      parameter (HDFE_ENDPOINT=1)
      integer HDFE_INTERNAL
      parameter (HDFE_INTERNAL=0)
      integer HDFE_COMP_SKPHUFF
      parameter (HDFE_COMP_SKPHUFF=3)
      integer HDFE_COMP_NONE
      parameter (HDFE_COMP_NONE=0)
      integer HDFE_NOPREVSUB
      parameter (HDFE_NOPREVSUB=-1)

end module HE4_class
