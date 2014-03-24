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
      integer, parameter :: cal_id = 510012
      integer, parameter :: L1B_LUN=299001
      integer, parameter :: L1B_LUN_cm=299002
      integer, parameter :: MCF_LUN = 511001

!run time parameter LUNs
      integer, parameter :: pgeversion_lun = 200105, processingcenter_lun = 200110, &
        instrumentname_lun = 200175, processlevel_lun = 200170, &
        processinghost_lun = 200115, authoraffiliation_lun = 200185, &
        authorname_lun = 200190, operationmode_lun = 200180, &
        OrbNum_LUN = 200200, ReprAct_LUN = 200120, using_resid_lun=200300, &
        write_resid_lun=200301, write_obs_lun=200302, using_cal_lun=200303, &
        no_ret_ps_lun=200304, no_ret_lun=200305, transient_chk=200306, &
        wmin_LUN=200307, wmax_LUN=200308, do_o3_LUN=200309, ThreshOrbNum_LUN=200310
 end module m_LUN_set
