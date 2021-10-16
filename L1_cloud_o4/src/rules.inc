$(ODIR)/cal_ecf.o: $(SDIR)/cal_ecf.f90 \
  $(ODIR)/m_read_DEM.o \
  $(ODIR)/m_read_GMI.o \
  $(ODIR)/m_read_hdf5.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/cal_ecf.f90 -o $(ODIR)/cal_ecf.o $(MODOUTFLG) $(ODIR)
$(ODIR)/cal_ocp.o: $(SDIR)/cal_ocp.f90 \
  $(ODIR)/m_read_DEM.o \
  $(ODIR)/m_read_GMI.o \
  $(ODIR)/m_read_hdf5.o \
  $(ODIR)/m_scd_adjust.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/cal_ocp.f90 -o $(ODIR)/cal_ocp.o $(MODOUTFLG) $(ODIR)
$(ODIR)/cal_pscene.o: $(SDIR)/cal_pscene.f90 \
  $(ODIR)/m_read_DEM.o \
  $(ODIR)/m_read_GMI.o \
  $(ODIR)/m_read_hdf5.o \
  $(ODIR)/m_scd_adjust.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/cal_pscene.f90 -o $(ODIR)/cal_pscene.o $(MODOUTFLG) $(ODIR)
$(ODIR)/Constants_class.o: $(SDIR)/Constants_class.f90 \
  $(ODIR)/DataTypeDef.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/Constants_class.f90 -o $(ODIR)/Constants_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/DataTypeDef.o: $(SDIR)/DataTypeDef.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/DataTypeDef.f90 -o $(ODIR)/DataTypeDef.o $(MODOUTFLG) $(ODIR)
$(ODIR)/ErrorHandler_class.o: $(SDIR)/ErrorHandler_class.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/ErrorHandler_class.f90 -o $(ODIR)/ErrorHandler_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/H5Util_class.o: $(SDIR)/H5Util_class.f90 \
  $(ODIR)/Constants_class.o \
  $(ODIR)/DataTypeDef.o \
  $(ODIR)/ErrorHandler_class.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/H5Util_class.f90 -o $(ODIR)/H5Util_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/hdfeos4_parameters.o: $(SDIR)/hdfeos4_parameters.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/hdfeos4_parameters.f90 -o $(ODIR)/hdfeos4_parameters.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_swreader.o: $(SDIR)/he5_swreader.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_swreader.f90 -o $(ODIR)/he5_swreader.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_input_clim.o: $(SDIR)/m_read_input_clim.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_input_clim.f90 -o $(ODIR)/m_read_input_clim.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_input_gler.o: $(SDIR)/m_read_input_gler.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_input_gler.f90 -o $(ODIR)/m_read_input_gler.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_input_tio.o: $(SDIR)/m_read_input_tio.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_input_tio.f90 -o $(ODIR)/m_read_input_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_write_output_tio.o: $(SDIR)/m_write_output_tio.f90 \
  $(ODIR)/m_read_input_tio.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_write_output_tio.f90 -o $(ODIR)/m_write_output_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_DEM.o: $(SDIR)/m_read_DEM.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_DEM.f90 -o $(ODIR)/m_read_DEM.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_GMI.o: $(SDIR)/m_read_GMI.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_GMI.f90 -o $(ODIR)/m_read_GMI.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_hdf5.o: $(SDIR)/m_read_hdf5.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_hdf5.f90 -o $(ODIR)/m_read_hdf5.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_input_kleipool.o: $(SDIR)/m_read_input_kleipool.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/he5_swreader.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_input_kleipool.f90 -o $(ODIR)/m_read_input_kleipool.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_lut.o: $(SDIR)/m_read_lut.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/he5_swreader.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_lut.f90 -o $(ODIR)/m_read_lut.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_vars.o: $(SDIR)/m_vars.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_vars.f90 -o $(ODIR)/m_vars.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_scd_adjust.o: $(SDIR)/m_scd_adjust.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_scd_adjust.f90 -o $(ODIR)/m_scd_adjust.o $(MODOUTFLG) $(ODIR)
$(ODIR)/GetConfigDoubleArr.o:	$(SDIR)/GetConfigDoubleArr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/GetConfigDoubleArr.c -o $(ODIR)/GetConfigDoubleArr.o
$(ODIR)/GetConfigDouble.o:	$(SDIR)/GetConfigDouble.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/GetConfigDouble.c -o $(ODIR)/GetConfigDouble.o
$(ODIR)/GetConfigInteger.o:	$(SDIR)/GetConfigInteger.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/GetConfigInteger.c -o $(ODIR)/GetConfigInteger.o
$(ODIR)/GetConfigRealArr.o:	$(SDIR)/GetConfigRealArr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/GetConfigRealArr.c -o $(ODIR)/GetConfigRealArr.o
$(ODIR)/GetConfigReal.o:	$(SDIR)/GetConfigReal.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/GetConfigReal.c -o $(ODIR)/GetConfigReal.o
$(ODIR)/GetConfigString.o:	$(SDIR)/GetConfigString.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/GetConfigString.c -o $(ODIR)/GetConfigString.o
$(ODIR)/LoadControlFile.o:	$(SDIR)/LoadControlFile.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/LoadControlFile.c -o $(ODIR)/LoadControlFile.o
$(ODIR)/Messages.o:	$(SDIR)/Messages.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/Messages.c -o $(ODIR)/Messages.o
$(ODIR)/MetaConfig.o:	$(SDIR)/MetaConfig.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/MetaConfig.c -o $(ODIR)/MetaConfig.o
$(ODIR)/OMI_localGranuleID.o:	$(SDIR)/OMI_localGranuleID.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_localGranuleID.c -o $(ODIR)/OMI_localGranuleID.o
$(ODIR)/OMI_MET_bindFORTRAN.o:	$(SDIR)/OMI_MET_bindFORTRAN.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_MET_bindFORTRAN.c -o $(ODIR)/OMI_MET_bindFORTRAN.o
$(ODIR)/OMI_PC_bindFORTRAN.o:	$(SDIR)/OMI_PC_bindFORTRAN.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_PC_bindFORTRAN.c -o $(ODIR)/OMI_PC_bindFORTRAN.o
$(ODIR)/OMI_PC_GetReference.o:	$(SDIR)/OMI_PC_GetReference.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_PC_GetReference.c -o $(ODIR)/OMI_PC_GetReference.o
$(ODIR)/OMI_SMF_setmsg.o:	$(SDIR)/OMI_SMF_setmsg.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_SMF_setmsg.c -o $(ODIR)/OMI_SMF_setmsg.o
$(ODIR)/ResetConfig.o:	$(SDIR)/ResetConfig.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ResetConfig.c -o $(ODIR)/ResetConfig.o
$(ODIR)/SetTrueIfLittleEndian.o:	$(SDIR)/SetTrueIfLittleEndian.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/SetTrueIfLittleEndian.c -o $(ODIR)/SetTrueIfLittleEndian.o
$(ODIR)/ShareArgsWithFortran.o:	$(SDIR)/ShareArgsWithFortran.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ShareArgsWithFortran.c -o $(ODIR)/ShareArgsWithFortran.o
$(ODIR)/str2lwr.o:	$(SDIR)/str2lwr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/str2lwr.c -o $(ODIR)/str2lwr.o
$(ODIR)/str2upr.o:	$(SDIR)/str2upr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/str2upr.c -o $(ODIR)/str2upr.o
$(ODIR)/strechr.o:	$(SDIR)/strechr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/strechr.c -o $(ODIR)/strechr.o
$(ODIR)/strxchr.o:	$(SDIR)/strxchr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/strxchr.c -o $(ODIR)/strxchr.o
$(ODIR)/unbufferSTDOUT.o:	$(SDIR)/unbufferSTDOUT.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/unbufferSTDOUT.c -o $(ODIR)/unbufferSTDOUT.o
$(ODIR)/ymlAddKey.o:	$(SDIR)/ymlAddKey.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlAddKey.c -o $(ODIR)/ymlAddKey.o
$(ODIR)/ymlAddVal.o:	$(SDIR)/ymlAddVal.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlAddVal.c -o $(ODIR)/ymlAddVal.o
$(ODIR)/ymlCfort.o:	$(SDIR)/ymlCfort.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlCfort.c -o $(ODIR)/ymlCfort.o
$(ODIR)/ymlChkTxt.o:	$(SDIR)/ymlChkTxt.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlChkTxt.c -o $(ODIR)/ymlChkTxt.o
$(ODIR)/ymlClrUsdFlg.o:	$(SDIR)/ymlClrUsdFlg.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlClrUsdFlg.c -o $(ODIR)/ymlClrUsdFlg.o
$(ODIR)/ymlDmpStk.o:	$(SDIR)/ymlDmpStk.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlDmpStk.c -o $(ODIR)/ymlDmpStk.o
$(ODIR)/ymlDmpTbl.o:	$(SDIR)/ymlDmpTbl.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlDmpTbl.c -o $(ODIR)/ymlDmpTbl.o
$(ODIR)/ymlErrors.o:	$(SDIR)/ymlErrors.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlErrors.c -o $(ODIR)/ymlErrors.o
$(ODIR)/ymlEstNrE.o:	$(SDIR)/ymlEstNrE.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlEstNrE.c -o $(ODIR)/ymlEstNrE.o
$(ODIR)/ymlGetArrStr.o:	$(SDIR)/ymlGetArrStr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlGetArrStr.c -o $(ODIR)/ymlGetArrStr.o
$(ODIR)/ymlGetKeys.o:	$(SDIR)/ymlGetKeys.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlGetKeys.c -o $(ODIR)/ymlGetKeys.o
$(ODIR)/ymlGetOneStr.o:	$(SDIR)/ymlGetOneStr.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlGetOneStr.c -o $(ODIR)/ymlGetOneStr.o
$(ODIR)/ymlIniTbl.o:	$(SDIR)/ymlIniTbl.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlIniTbl.c -o $(ODIR)/ymlIniTbl.o
$(ODIR)/ymlMaxDtl.o:	$(SDIR)/ymlMaxDtl.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlMaxDtl.c -o $(ODIR)/ymlMaxDtl.o
$(ODIR)/ymlNrmKey.o:	$(SDIR)/ymlNrmKey.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlNrmKey.c -o $(ODIR)/ymlNrmKey.o
$(ODIR)/ymlOpen.o:	$(SDIR)/ymlOpen.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlOpen.c -o $(ODIR)/ymlOpen.o
$(ODIR)/ymlParTxt.o:	$(SDIR)/ymlParTxt.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlParTxt.c -o $(ODIR)/ymlParTxt.o
$(ODIR)/ymlRdFile.o:	$(SDIR)/ymlRdFile.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlRdFile.c -o $(ODIR)/ymlRdFile.o
$(ODIR)/ymlStkPsh.o:	$(SDIR)/ymlStkPsh.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/ymlStkPsh.c -o $(ODIR)/ymlStkPsh.o
