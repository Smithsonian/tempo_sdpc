module m_LUN_set
  integer, parameter :: L2_out = 520001
  integer, parameter :: chl_id = 510004
  integer, parameter :: refl_id = 510005
  integer, parameter :: thresh_id = 510006
  integer, parameter :: resid_id_early = 510007
  integer, parameter :: resid_id_late = 510008
  integer, parameter :: o3_id = 510009
  integer, parameter :: ref_id = 510007
  integer, parameter :: terr_prs_id = 510003
  integer, parameter :: oc_ram_id = 510011
  integer, parameter :: ler354_id = 510012
  integer, parameter :: IRR1B_file = 510001
  integer, parameter :: ring_id = 510010
  integer, parameter :: L1B_LUN=299001
  integer, parameter :: L1B_LUN_cm=299002
  integer, parameter :: MCF_LUN = 511001

  !run time parameter LUNs
  integer, parameter :: pgeversion_lun = 200105, processingcenter_lun = 200110, &
       instrumentname_lun = 200175, processlevel_lun = 200170, &
       processinghost_lun = 200115, authoraffiliation_lun = 200185, &
       authorname_lun = 200190, operationmode_lun = 200180, &
       versionid_lun = 200205, &
       OrbNum_LUN = 200200, ReprAct_LUN = 200120, using_resid_lun=200300, &
       write_resid_lun=200301, write_obs_lun=200302, &
       no_ret_ps_lun=200304, no_ret_lun=200305, transient_chk=200306, &
       wmin_LUN=200307, wmax_LUN=200308, do_o3_LUN=200309, ThreshOrbNum_LUN=200310
  integer, parameter :: test_solar_LUN = 200311, add_shift_LUN=200312
  integer, parameter :: using_spline_LUN = 200313
  integer, parameter :: mdlist_LUN=511002
end module m_LUN_set
