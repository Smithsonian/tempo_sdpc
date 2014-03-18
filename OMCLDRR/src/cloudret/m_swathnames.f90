MODULE m_swathnames

! visswath = swathname for VIS
! uv1swath = swathname for UV1
! uv2swath = swathname for UV2

   CHARACTER (LEN=200), PARAMETER ::  visswath = "Earth VIS Swath"
   CHARACTER (LEN=200), PARAMETER ::  uv1swath = "Earth UV-1 Swath"
   CHARACTER (LEN=200), PARAMETER ::  uv2swath = "Earth UV-2 Swath"
   CHARACTER (LEN=200), PARAMETER ::  sunvisswath = "Sun Volume VIS Swath" 
   CHARACTER (LEN=200), PARAMETER ::  sunuv2swath = "Sun Volume UV-2 Swath" 

   CHARACTER (LEN=200), PARAMETER ::  visswathz = "Earth VIS Swath (60x751x4)"
   CHARACTER (LEN=200), PARAMETER ::  uv1swathz = "Earth UV-1 Swath (60x557x4)" 
   CHARACTER (LEN=200), PARAMETER ::  uv2swathz = "Earth UV-2 Swath (60x557x4)"
   CHARACTER (LEN=200), PARAMETER ::  sunvisswathz = "Sun Volume VIS Swath (60x751x4)"
   CHARACTER (LEN=200), PARAMETER ::  sunuv2swathz = "Sun Volume UV-2 Swath (60x751x4)"

   logical :: visz=.false.
   logical :: uvsz=.false.
   logical :: vis=.false.

END MODULE m_swathnames
