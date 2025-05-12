!> include problematic PGS code via module to avoid warnings during build
!
module m_pgs_include
  include 'PGS_IO.f'
  include 'PGS_OMI_1900.f'
  include 'PGS_OMSAO_52510.f'
  include 'PGS_SMF.f'
  include 'PGS_MET_13.f'
end module m_pgs_include
