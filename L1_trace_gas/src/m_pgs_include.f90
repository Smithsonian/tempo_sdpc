! include PGS code in one location to avoid warnings during build
!
!
module m_pgs_include
  include "PGS_MET.f"
  include "PGS_MET_13.f"
  include "PGS_PC.f"
  include 'PGS_PC_9.f'
  include 'PGS_SMF.f'
  include 'PGS_IO.f'
  include 'PGS_TD_3.f'
  include 'PGS_OMSAO_52500.f'
  include 'PGS_OMI_1900.f'
end module m_pgs_include
