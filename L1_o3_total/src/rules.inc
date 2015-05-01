$(ODIR)/HDF4_class.o: $(SDIR)/HDF4_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/HDF4_class.f90 -o $(ODIR)/HDF4_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/HE4_class.o: $(SDIR)/HE4_class.f90 \
  $(ODIR)/HDF4_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/HE4_class.f90 -o $(ODIR)/HE4_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/HE5_class.o: $(SDIR)/HE5_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/HE5_class.f90 -o $(ODIR)/HE5_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L1B_class.o: $(SDIR)/L1B_class.f90 \
  $(ODIR)/L1B_geoang_class.o \
  $(ODIR)/L1B_radirr_class.o \
  $(ODIR)/l1b_tio_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L1B_class.f90 -o $(ODIR)/L1B_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L1B_geoang_class.o: $(SDIR)/L1B_geoang_class.f90 \
  $(ODIR)/HE4_class.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L1B_geoang_class.f90 -o $(ODIR)/L1B_geoang_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L1B_getNames_m.o: $(SDIR)/L1B_getNames_m.f90 \
  $(ODIR)/HE4_class.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L1B_getNames_m.f90 -o $(ODIR)/L1B_getNames_m.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L1B_metaData_class.o: $(SDIR)/L1B_metaData_class.f90 \
  $(ODIR)/HDF4_class.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_MET_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L1B_metaData_class.f90 -o $(ODIR)/L1B_metaData_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L1B_radirr_class.o: $(SDIR)/L1B_radirr_class.f90 \
  $(ODIR)/HE4_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L1B_radirr_class.f90 -o $(ODIR)/L1B_radirr_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/l1b_tio_class.o: $(SDIR)/l1b_tio_class.f90 \
  $(ODIR)/o3t_names_module.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/l1b_tio_class.f90 -o $(ODIR)/l1b_tio_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L2_attr_class.o: $(SDIR)/L2_attr_class.f90 \
  $(ODIR)/HE5_class.o \
  $(ODIR)/L1B_metaData_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_MET_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/PGS_TD_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L2_attr_class.f90 -o $(ODIR)/L2_attr_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L2_data_structure.o: $(SDIR)/L2_data_structure.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L2_data_structure.f90 -o $(ODIR)/L2_data_structure.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L2_metaData_class.o: $(SDIR)/L2_metaData_class.f90 \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_MET_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/PGS_TD_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L2_metaData_class.f90 -o $(ODIR)/L2_metaData_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_anomflg.o: $(SDIR)/m_anomflg.f90 \
  $(ODIR)/HE4_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_anomflg.f90 -o $(ODIR)/m_anomflg.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_nvalc.o: $(SDIR)/m_nvalc.f90 \
  $(ODIR)/HE4_class.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_nvalc.f90 -o $(ODIR)/m_nvalc.o $(MODOUTFLG) $(ODIR)
$(ODIR)/o3t_names_module.o: $(SDIR)/o3t_names_module.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/o3t_names_module.f90 -o $(ODIR)/o3t_names_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_apriori_class.o: $(SDIR)/O3T_apriori_class.f90 \
  $(ODIR)/O3T_stnprof_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_apriori_class.f90 -o $(ODIR)/O3T_apriori_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_class.o: $(SDIR)/O3T_class.f90 \
  $(ODIR)/O3T_apriori_class.o \
  $(ODIR)/O3T_const.o \
  $(ODIR)/O3T_dndx_class.o \
  $(ODIR)/O3T_dndx_m.o \
  $(ODIR)/O3T_iztrsb_m.o \
  $(ODIR)/O3T_lpolycoef_class.o \
  $(ODIR)/O3T_nval_class.o \
  $(ODIR)/O3T_pixel_class.o \
  $(ODIR)/O3T_stnprof_class.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_class.f90 -o $(ODIR)/O3T_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_cloudPres_class.o: $(SDIR)/O3T_cloudPres_class.f90 \
  $(ODIR)/OMI_L2reader_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_MET_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_cloudPres_class.f90 -o $(ODIR)/O3T_cloudPres_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_const.o: $(SDIR)/O3T_const.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_const.f90 -o $(ODIR)/O3T_const.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_dndx_class.o: $(SDIR)/O3T_dndx_class.f90 \
  $(ODIR)/O3T_LUT_fs.o \
  $(ODIR)/O3T_nval_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_dndx_class.f90 -o $(ODIR)/O3T_dndx_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_dndx_m.o: $(SDIR)/O3T_dndx_m.f90 \
  $(ODIR)/O3T_dndx_class.o \
  $(ODIR)/O3T_lpolycoef_class.o \
  $(ODIR)/O3T_lpolyinterp_class.o \
  $(ODIR)/O3T_pixel_class.o \
  $(ODIR)/O3T_stnprof_class.o \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_dndx_m.f90 -o $(ODIR)/O3T_dndx_m.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_irrad_class.o: $(SDIR)/O3T_irrad_class.f90 \
  $(ODIR)/HE4_class.o \
  $(ODIR)/L1B_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_irrad_class.f90 -o $(ODIR)/O3T_irrad_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_iztrsb_m.o: $(SDIR)/O3T_iztrsb_m.f90 \
  $(ODIR)/O3T_lpolycoef_class.o \
  $(ODIR)/O3T_nval_class.o \
  $(ODIR)/O3T_pixel_class.o \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_iztrsb_m.f90 -o $(ODIR)/O3T_iztrsb_m.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_L2output_class.o: $(SDIR)/O3T_L2output_class.f90 \
  $(ODIR)/O3T_omto3_fs.o \
  $(ODIR)/O3T_pixel_class.o \
  $(ODIR)/O3T_radgeo_class.o \
  $(ODIR)/OMI_L2writer_class.o \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_L2output_class.f90 -o $(ODIR)/O3T_L2output_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_lpolycoef_class.o: $(SDIR)/O3T_lpolycoef_class.f90 \
  $(ODIR)/O3T_nval_class.o \
  $(ODIR)/O3T_pixel_class.o \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_lpolycoef_class.f90 -o $(ODIR)/O3T_lpolycoef_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_lpolyinterp_class.o: $(SDIR)/O3T_lpolyinterp_class.f90 \
  $(ODIR)/O3T_dndx_class.o \
  $(ODIR)/O3T_lpolycoef_class.o \
  $(ODIR)/O3T_nval_class.o \
  $(ODIR)/O3T_pixel_class.o \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_lpolyinterp_class.f90 -o $(ODIR)/O3T_lpolyinterp_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_LUT_fs.o: $(SDIR)/O3T_LUT_fs.f90 \
  $(ODIR)/UTIL_lh5_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_LUT_fs.f90 -o $(ODIR)/O3T_LUT_fs.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_nval_class.o: $(SDIR)/O3T_nval_class.f90 \
  $(ODIR)/O3T_LUT_fs.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_nval_class.f90 -o $(ODIR)/O3T_nval_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_omto3_fs.o: $(SDIR)/O3T_omto3_fs.f90 \
  $(ODIR)/HE5_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_omto3_fs.f90 -o $(ODIR)/O3T_omto3_fs.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_pixel_class.o: $(SDIR)/O3T_pixel_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_pixel_class.f90 -o $(ODIR)/O3T_pixel_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_QA_class.o: $(SDIR)/O3T_QA_class.f90 \
  $(ODIR)/L1B_class.o \
  $(ODIR)/PGS_MET_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_QA_class.f90 -o $(ODIR)/O3T_QA_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_radgeo_class.o: $(SDIR)/O3T_radgeo_class.f90 \
  $(ODIR)/L1B_class.o \
  $(ODIR)/l1b_tio_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_radgeo_class.f90 -o $(ODIR)/O3T_radgeo_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_so2_class.o: $(SDIR)/O3T_so2_class.f90 \
  $(ODIR)/O3T_const.o \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_so2_class.f90 -o $(ODIR)/O3T_so2_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3T_stnprof_class.o: $(SDIR)/O3T_stnprof_class.f90 \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/O3T_stnprof_class.f90 -o $(ODIR)/O3T_stnprof_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_copyHE4toHE5.o: $(SDIR)/OMI_copyHE4toHE5.f90 \
  $(ODIR)/HE4_class.o \
  $(ODIR)/HE5_class.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/OMI_copyHE4toHE5.f90 -o $(ODIR)/OMI_copyHE4toHE5.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_L2reader_class.o: $(SDIR)/OMI_L2reader_class.f90 \
  $(ODIR)/HE5_class.o \
  $(ODIR)/L2_data_structure.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/OMI_L2reader_class.f90 -o $(ODIR)/OMI_L2reader_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_L2writer_class.o: $(SDIR)/OMI_L2writer_class.f90 \
  $(ODIR)/HE5_class.o \
  $(ODIR)/L2_data_structure.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/OMI_L2writer_class.f90 -o $(ODIR)/OMI_L2writer_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_LUN_set.o: $(SDIR)/OMI_LUN_set.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/OMI_LUN_set.f90 -o $(ODIR)/OMI_LUN_set.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_SMF_class.o: $(SDIR)/OMI_SMF_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/OMI_SMF_class.f90 -o $(ODIR)/OMI_SMF_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_MET_class.o: $(SDIR)/PGS_MET_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/PGS_MET_class.f90 -o $(ODIR)/PGS_MET_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_PC_class.o: $(SDIR)/PGS_PC_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/PGS_PC_class.f90 -o $(ODIR)/PGS_PC_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_TD_class.o: $(SDIR)/PGS_TD_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/PGS_TD_class.f90 -o $(ODIR)/PGS_TD_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/UTIL_lh5_class.o: $(SDIR)/UTIL_lh5_class.f90 \
  $(ODIR)/OMI_SMF_class.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/UTIL_lh5_class.f90 -o $(ODIR)/UTIL_lh5_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/UTIL_tools_class.o: $(SDIR)/UTIL_tools_class.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/UTIL_tools_class.f90 -o $(ODIR)/UTIL_tools_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/DFKNTwrapper.o:	$(SDIR)/DFKNTwrapper.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/DFKNTwrapper.c -o $(ODIR)/DFKNTwrapper.o
$(ODIR)/HE5Tget_size.o:	$(SDIR)/HE5Tget_size.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/HE5Tget_size.c -o $(ODIR)/HE5Tget_size.o
$(ODIR)/L1B_BLK_reader.o:	$(SDIR)/L1B_BLK_reader.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/L1B_BLK_reader.c -o $(ODIR)/L1B_BLK_reader.o
$(ODIR)/lhdf.o:	$(SDIR)/lhdf.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/lhdf.c -o $(ODIR)/lhdf.o
$(ODIR)/OMI_localGranuleID.o:	$(SDIR)/OMI_localGranuleID.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_localGranuleID.c -o $(ODIR)/OMI_localGranuleID.o
$(ODIR)/OMI_mmddInterp.o:	$(SDIR)/OMI_mmddInterp.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_mmddInterp.c -o $(ODIR)/OMI_mmddInterp.o
$(ODIR)/OMI_pixGetCldPres.o:	$(SDIR)/OMI_pixGetCldPres.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_pixGetCldPres.c -o $(ODIR)/OMI_pixGetCldPres.o
$(ODIR)/OMI_pixGetSnowIce.o:	$(SDIR)/OMI_pixGetSnowIce.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_pixGetSnowIce.c -o $(ODIR)/OMI_pixGetSnowIce.o
$(ODIR)/OMI_pixGetTerPres.o:	$(SDIR)/OMI_pixGetTerPres.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_pixGetTerPres.c -o $(ODIR)/OMI_pixGetTerPres.o
$(ODIR)/OMI_SMF_setmsg.o:	$(SDIR)/OMI_SMF_setmsg.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_SMF_setmsg.c -o $(ODIR)/OMI_SMF_setmsg.o
$(ODIR)/rFill.o:	$(SDIR)/rFill.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/rFill.c -o $(ODIR)/rFill.o
