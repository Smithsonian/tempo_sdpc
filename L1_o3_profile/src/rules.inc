$(ODIR)/adj_measurement_data.o: $(SDIR)/adj_measurement_data.f90 \
  $(ODIR)/m_amf_conversion.o \
  $(ODIR)/m_angle_sat2toa.o \
  $(ODIR)/prepare_databases.o \
  $(ODIR)/m_rad_fit_vary.o \
  $(ODIR)/m_rad_shisqu.o \
  $(ODIR)/m_rad_wavcal_vary.o \
  $(ODIR)/m_radiance_wavcal.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_gome_data_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/adj_measurement_data.f90 -o $(ODIR)/adj_measurement_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/avg_band.o: $(SDIR)/avg_band.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/avg_band.f90 -o $(ODIR)/avg_band.o $(MODOUTFLG) $(ODIR)
$(ODIR)/bounded_nonlin_ls.o: $(SDIR)/bounded_nonlin_ls.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/bounded_nonlin_ls.f90 -o $(ODIR)/bounded_nonlin_ls.o $(MODOUTFLG) $(ODIR)
$(ODIR)/fitting_functions.o: $(SDIR)/fitting_functions.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/spectra.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/fitting_functions.f90 -o $(ODIR)/fitting_functions.o $(MODOUTFLG) $(ODIR)
$(ODIR)/get_cloud.o: $(SDIR)/get_cloud.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/prepare_atmosphere.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/get_cloud.f90 -o $(ODIR)/get_cloud.o $(MODOUTFLG) $(ODIR)
$(ODIR)/gsvd_data_module.o: $(SDIR)/gsvd_data_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/gsvd_data_module.f90 -o $(ODIR)/gsvd_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/gsvd_o3prof_utilities.o: $(SDIR)/gsvd_o3prof_utilities.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/gsvd_o3prof_utilities.f90 -o $(ODIR)/gsvd_o3prof_utilities.o $(MODOUTFLG) $(ODIR)
$(ODIR)/hdfeos4_parameters.o: $(SDIR)/hdfeos4_parameters.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/hdfeos4_parameters.f90 -o $(ODIR)/hdfeos4_parameters.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_input_tools.o: $(SDIR)/he5_input_tools.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_input_tools.f90 -o $(ODIR)/he5_input_tools.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_l2_fs.o: $(SDIR)/he5_l2_fs.f90 \
  $(ODIR)/OMSAO_he5_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_l2_fs.f90 -o $(ODIR)/he5_l2_fs.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_l2writer_class.o: $(SDIR)/he5_l2writer_class.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_omidata_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_l2writer_class.f90 -o $(ODIR)/he5_l2writer_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_output_module.o: $(SDIR)/he5_output_module.f90 \
  $(ODIR)/he5_l2_fs.o \
  $(ODIR)/he5_l2writer_class.o \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/L2_attr_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_output_module.f90 -o $(ODIR)/he5_output_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L1B_Reader_class.o: $(SDIR)/L1B_Reader_class.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/Szoom_Parameter_Module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/L1B_Reader_class.f90 -o $(ODIR)/L1B_Reader_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L2_attr_class.o: $(SDIR)/L2_attr_class.f90 \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/PGS_MET_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/PGS_TD_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/L2_attr_class.f90 -o $(ODIR)/L2_attr_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/lidort_prof_utilities.o: $(SDIR)/lidort_prof_utilities.f90 \
  $(ODIR)/avg_band.o \
  $(ODIR)/m_dataspline.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_gauss.o \
  $(ODIR)/m_lidort_prof_prep.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/lidort_prof_utilities.f90 -o $(ODIR)/lidort_prof_utilities.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_amf_conversion.o: $(SDIR)/m_amf_conversion.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_amf_conversion.f90 -o $(ODIR)/m_amf_conversion.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_angle_sat2toa.o: $(SDIR)/m_angle_sat2toa.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_angle_sat2toa.f90 -o $(ODIR)/m_angle_sat2toa.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_cal_fit_one.o: $(SDIR)/m_cal_fit_one.f90 \
  $(ODIR)/fitting_functions.o \
  $(ODIR)/m_specfit.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_cal_fit_one.f90 -o $(ODIR)/m_cal_fit_one.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_convert_coadd.o: $(SDIR)/m_convert_coadd.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_convert_coadd.f90 -o $(ODIR)/m_convert_coadd.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_dataspline.o: $(SDIR)/m_dataspline.f90 \
  $(ODIR)/avg_band.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_gauss.o \
  $(ODIR)/m_triangle.o \
  $(ODIR)/m_voigt.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_dataspline.f90 -o $(ODIR)/m_dataspline.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_ezspline_interpolation.o: $(SDIR)/m_ezspline_interpolation.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_ezspline_interpolation.f90 -o $(ODIR)/m_ezspline_interpolation.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_gauss.o: $(SDIR)/m_gauss.f90 \
  $(ODIR)/m_triangle.o \
  $(ODIR)/m_voigt.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_gauss.f90 -o $(ODIR)/m_gauss.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_apriori_covar.o: $(SDIR)/m_get_apriori_covar.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/prepare_atmosphere.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_apriori_covar.f90 -o $(ODIR)/m_get_apriori_covar.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_initial_albedo.o: $(SDIR)/m_get_initial_albedo.f90 \
  $(ODIR)/lidort_prof_utilities.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_initial_albedo.f90 -o $(ODIR)/m_get_initial_albedo.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_raman.o: $(SDIR)/m_get_raman.f90 \
  $(ODIR)/avg_band.o \
  $(ODIR)/lidort_prof_utilities.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_lidort_prof_prep.o \
  $(ODIR)/m_raman.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_raman.f90 -o $(ODIR)/m_get_raman.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_reg_matrix.o: $(SDIR)/m_get_reg_matrix.f90 \
  $(ODIR)/gsvd_o3prof_utilities.o \
  $(ODIR)/m_gsvd_lcurve_gcv.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_reg_matrix.f90 -o $(ODIR)/m_get_reg_matrix.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_gsvd_lcurve_gcv.o: $(SDIR)/m_gsvd_lcurve_gcv.f90 \
  $(ODIR)/gsvd_data_module.o \
  $(ODIR)/gsvd_o3prof_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_gsvd_lcurve_gcv.f90 -o $(ODIR)/m_gsvd_lcurve_gcv.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_int2string.o: $(SDIR)/m_int2string.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_int2string.f90 -o $(ODIR)/m_int2string.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_lidort_prof_prep.o: $(SDIR)/m_lidort_prof_prep.f90 \
  $(ODIR)/avg_band.o \
  $(ODIR)/m_dataspline.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_lidort_prof_prep.f90 -o $(ODIR)/m_lidort_prof_prep.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_make_atm.o: $(SDIR)/m_make_atm.f90 \
  $(ODIR)/get_cloud.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_read_aerosol_prof.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/prepare_atmosphere.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_make_atm.f90 -o $(ODIR)/m_make_atm.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_oe_inversion.o: $(SDIR)/m_oe_inversion.f90 \
  $(ODIR)/gsvd_o3prof_utilities.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_oe_inversion.f90 -o $(ODIR)/m_oe_inversion.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_ozprof_inverse.o: $(SDIR)/m_ozprof_inverse.f90 \
  $(ODIR)/bounded_nonlin_ls.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_raman.o \
  $(ODIR)/m_gsvd_lcurve_gcv.o \
  $(ODIR)/m_make_atm.o \
  $(ODIR)/m_oe_inversion.o \
  $(ODIR)/m_pseudo_model.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/prepare_atmosphere.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_ozprof_inverse.f90 -o $(ODIR)/m_ozprof_inverse.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_prepare_refspecs.o: $(SDIR)/m_prepare_refspecs.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_prepare_refspecs.f90 -o $(ODIR)/m_prepare_refspecs.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_rad_fit_vary.o: $(SDIR)/m_rad_fit_vary.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_subtract_poly.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_rad_fit_vary.f90 -o $(ODIR)/m_rad_fit_vary.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_radiance_wavcal.o: $(SDIR)/m_radiance_wavcal.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_radiance_wavcal.f90 -o $(ODIR)/m_radiance_wavcal.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_rad_wavcal_vary.o: $(SDIR)/m_rad_wavcal_vary.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_subtract_poly.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_rad_wavcal_vary.f90 -o $(ODIR)/m_rad_wavcal_vary.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_raman.o: $(SDIR)/m_raman.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_raman.f90 -o $(ODIR)/m_raman.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_solar_fit.o: $(SDIR)/m_solar_fit.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_solar_fit.f90 -o $(ODIR)/m_solar_fit.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_solar_wavcal_vary.o: $(SDIR)/m_solar_wavcal_vary.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_subtract_poly.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_solar_wavcal_vary.f90 -o $(ODIR)/m_solar_wavcal_vary.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_specfit.o: $(SDIR)/m_specfit.f90 \
  $(ODIR)/bounded_nonlin_ls.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_specfit.f90 -o $(ODIR)/m_specfit.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_subtract_poly.o: $(SDIR)/m_subtract_poly.f90 \
  $(ODIR)/bounded_nonlin_ls.o \
  $(ODIR)/fitting_functions.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_subtract_poly.f90 -o $(ODIR)/m_subtract_poly.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_triangle.o: $(SDIR)/m_triangle.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_triangle.f90 -o $(ODIR)/m_triangle.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_undersample.o: $(SDIR)/m_undersample.f90 \
  $(ODIR)/m_dataspline.o \
  $(ODIR)/m_gauss.o \
  $(ODIR)/m_triangle.o \
  $(ODIR)/m_voigt.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_undersample.f90 -o $(ODIR)/m_undersample.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_voigt.o: $(SDIR)/m_voigt.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_voigt.f90 -o $(ODIR)/m_voigt.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_cross_calibrate.o: $(SDIR)/omi_cross_calibrate.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_rad_fit_vary.o \
  $(ODIR)/m_rad_wavcal_vary.o \
  $(ODIR)/m_radiance_wavcal.o \
  $(ODIR)/m_solar_fit.o \
  $(ODIR)/solar_fit_vary.o \
  $(ODIR)/m_solar_wavcal_vary.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_cross_calibrate.f90 -o $(ODIR)/omi_cross_calibrate.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_fitting_aux.o: $(SDIR)/omi_fitting_aux.f90 \
  $(ODIR)/adj_measurement_data.o \
  $(ODIR)/avg_band.o \
  $(ODIR)/get_cloud.o \
  $(ODIR)/m_angle_sat2toa.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/prepare_databases.o \
  $(ODIR)/m_prepare_refspecs.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omicloud_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_fitting_aux.f90 -o $(ODIR)/omi_fitting_aux.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_omi_fitting_process.o: $(SDIR)/m_omi_fitting_process.f90 \
  $(ODIR)/he5_output_module.o \
  $(ODIR)/m_read_cloud_tio.o \
  $(ODIR)/m_read_geo_tio.o \
  $(ODIR)/m_read_l1_tio.o \
  $(ODIR)/m_read_metadata_tio.o \
  $(ODIR)/m_specfit_ozprof.o \
  $(ODIR)/m_write_final.o \
  $(ODIR)/omi_cross_calibrate.o \
  $(ODIR)/omi_fitting_aux.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMI_MET_class.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_omicloud_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/OZPROF_metaDef.o \
  $(ODIR)/tio_output_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_omi_fitting_process.f90 -o $(ODIR)/m_omi_fitting_process.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_LUN_set.o: $(SDIR)/OMI_LUN_set.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMI_LUN_set.f90 -o $(ODIR)/OMI_LUN_set.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_MET_class.o: $(SDIR)/OMI_MET_class.f90 \
  $(ODIR)/OMI_SMF_class.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/PGS_MET_class.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/PGS_TD_class.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMI_MET_class.f90 -o $(ODIR)/OMI_MET_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_read_l1b_data.o: $(SDIR)/omi_read_l1b_data.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/m_gauss.o \
  $(ODIR)/m_read_l1_tio.o \
  $(ODIR)/m_triangle.o \
  $(ODIR)/omi_cross_calibrate.o \
  $(ODIR)/omi_fitting_aux.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_read_l1b_data.f90 -o $(ODIR)/omi_read_l1b_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_SMF_class.o: $(SDIR)/OMI_SMF_class.f90 \
  $(ODIR)/OMSAO_errstat_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMI_SMF_class.f90 -o $(ODIR)/OMI_SMF_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_errstat_module.o: $(SDIR)/OMSAO_errstat_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_errstat_module.f90 -o $(ODIR)/OMSAO_errstat_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_gome_data_module.o: $(SDIR)/OMSAO_gome_data_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_gome_data_module.f90 -o $(ODIR)/OMSAO_gome_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_he5_module.o: $(SDIR)/OMSAO_he5_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_he5_module.f90 -o $(ODIR)/OMSAO_he5_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_indices_module.o: $(SDIR)/OMSAO_indices_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_indices_module.f90 -o $(ODIR)/OMSAO_indices_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_omicloud_module.o: $(SDIR)/OMSAO_omicloud_module.f90 \
  $(ODIR)/he5_input_tools.o \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_omicloud_module.f90 -o $(ODIR)/OMSAO_omicloud_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_omidata_module.o: $(SDIR)/OMSAO_omidata_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_omidata_module.f90 -o $(ODIR)/OMSAO_omidata_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_parameters_module.o: $(SDIR)/OMSAO_parameters_module.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_parameters_module.f90 -o $(ODIR)/OMSAO_parameters_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_pixelcorner_module.o: $(SDIR)/OMSAO_pixelcorner_module.f90 \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/m_angle_sat2toa.o \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_pixelcorner_module.f90 -o $(ODIR)/OMSAO_pixelcorner_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_precision_module.o: $(SDIR)/OMSAO_precision_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_precision_module.f90 -o $(ODIR)/OMSAO_precision_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_slitfunction_module.o: $(SDIR)/OMSAO_slitfunction_module.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_slitfunction_module.f90 -o $(ODIR)/OMSAO_slitfunction_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_variables_module.o: $(SDIR)/OMSAO_variables_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_variables_module.f90 -o $(ODIR)/OMSAO_variables_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/ozprof_data_module.o: $(SDIR)/ozprof_data_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/ozprof_data_module.f90 -o $(ODIR)/ozprof_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OZPROF_metaDef.o: $(SDIR)/OZPROF_metaDef.f90 \
  $(ODIR)/OMI_MET_class.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OZPROF_metaDef.f90 -o $(ODIR)/OZPROF_metaDef.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_MET_class.o: $(SDIR)/PGS_MET_class.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/PGS_MET_class.f90 -o $(ODIR)/PGS_MET_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_PC_class.o: $(SDIR)/PGS_PC_class.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/PGS_PC_class.f90 -o $(ODIR)/PGS_PC_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_TD_class.o: $(SDIR)/PGS_TD_class.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/PGS_TD_class.f90 -o $(ODIR)/PGS_TD_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/prepare_atmosphere.o: $(SDIR)/prepare_atmosphere.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/prepare_atmosphere.f90 -o $(ODIR)/prepare_atmosphere.o $(MODOUTFLG) $(ODIR)
$(ODIR)/prepare_databases.o: $(SDIR)/prepare_databases.f90 \
  $(ODIR)/avg_band.o \
  $(ODIR)/m_dataspline.o \
  $(ODIR)/m_prepare_refspecs.o \
  $(ODIR)/m_undersample.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/prepare_databases.f90 -o $(ODIR)/prepare_databases.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_pseudo_model.o: $(SDIR)/m_pseudo_model.f90 \
  $(ODIR)/adj_measurement_data.o \
  $(ODIR)/lidort_prof_utilities.o \
  $(ODIR)/m_spectra_reflectance.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_pseudo_model.f90 -o $(ODIR)/m_pseudo_model.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_radiance_fit.o: $(SDIR)/m_radiance_fit.f90 \
  $(ODIR)/fitting_functions.o \
  $(ODIR)/m_specfit.o \
  $(ODIR)/m_subtract_cubic.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_radiance_fit.f90 -o $(ODIR)/m_radiance_fit.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_rad_shisqu.o: $(SDIR)/m_rad_shisqu.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_rad_shisqu.f90 -o $(ODIR)/m_rad_shisqu.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_aerosol_prof.o: $(SDIR)/m_read_aerosol_prof.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/prepare_atmosphere.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_aerosol_prof.f90 -o $(ODIR)/m_read_aerosol_prof.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_fitting_controls.o: $(SDIR)/m_read_fitting_controls.f90 \
  $(ODIR)/m_read_ozprof_input.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_gome_data_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/UTIL_tools_class.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_fitting_controls.f90 -o $(ODIR)/m_read_fitting_controls.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_ozprof_input.o: $(SDIR)/m_read_ozprof_input.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_ozprof_input.f90 -o $(ODIR)/m_read_ozprof_input.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_reference_spectra.o: $(SDIR)/m_read_reference_spectra.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_reference_spectra.f90 -o $(ODIR)/m_read_reference_spectra.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_set_cldalb.o: $(SDIR)/m_set_cldalb.f90 \
  $(ODIR)/get_cloud.o \
  $(ODIR)/m_get_initial_albedo.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_gome_data_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_set_cldalb.f90 -o $(ODIR)/m_set_cldalb.o $(MODOUTFLG) $(ODIR)
$(ODIR)/solar_fit_vary.o: $(SDIR)/solar_fit_vary.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_subtract_poly.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/solar_fit_vary.f90 -o $(ODIR)/solar_fit_vary.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_specfit_ozprof.o: $(SDIR)/m_specfit_ozprof.f90 \
  $(ODIR)/m_get_apriori_covar.o \
  $(ODIR)/m_get_reg_matrix.o \
  $(ODIR)/m_make_atm.o \
  $(ODIR)/m_ozprof_inverse.o \
  $(ODIR)/m_set_cldalb.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/prepare_atmosphere.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_specfit_ozprof.f90 -o $(ODIR)/m_specfit_ozprof.o $(MODOUTFLG) $(ODIR)
$(ODIR)/spectra.o: $(SDIR)/spectra.f90 \
  $(ODIR)/m_gauss.o \
  $(ODIR)/m_triangle.o \
  $(ODIR)/m_voigt.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/spectra.f90 -o $(ODIR)/spectra.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_spectra_reflectance.o: $(SDIR)/m_spectra_reflectance.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_spectra_reflectance.f90 -o $(ODIR)/m_spectra_reflectance.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_subtract_cubic.o: $(SDIR)/m_subtract_cubic.f90 \
  $(ODIR)/bounded_nonlin_ls.o \
  $(ODIR)/fitting_functions.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_subtract_cubic.f90 -o $(ODIR)/m_subtract_cubic.o $(MODOUTFLG) $(ODIR)
$(ODIR)/Szoom_Parameter_Module.o: $(SDIR)/Szoom_Parameter_Module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/Szoom_Parameter_Module.f90 -o $(ODIR)/Szoom_Parameter_Module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/utilities.o: $(SDIR)/utilities.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/utilities.f90 -o $(ODIR)/utilities.o $(MODOUTFLG) $(ODIR)
$(ODIR)/UTIL_tools_class.o: $(SDIR)/UTIL_tools_class.f90 \
  $(ODIR)/OMSAO_errstat_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/UTIL_tools_class.f90 -o $(ODIR)/UTIL_tools_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_write_final.o: $(SDIR)/m_write_final.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_write_final.f90 -o $(ODIR)/m_write_final.o $(MODOUTFLG) $(ODIR)
$(ODIR)/tio_output_module.o: $(SDIR)/tio_output_module.f90 \
  $(ODIR)/m_o3p_params.o \
  $(ODIR)/m_read_geo_tio.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/tio_output_module.f90 -o $(ODIR)/tio_output_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/o3p_names_module.o: $(SDIR)/o3p_names_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/o3p_names_module.f90 -o $(ODIR)/o3p_names_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_cloud_tio.o: $(SDIR)/m_read_cloud_tio.f90 \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_omicloud_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_cloud_tio.f90 -o $(ODIR)/m_read_cloud_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_l1_tio.o: $(SDIR)/m_read_l1_tio.f90 \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_l1_tio.f90 -o $(ODIR)/m_read_l1_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_geo_tio.o: $(SDIR)/m_read_geo_tio.f90 \
  $(ODIR)/m_angle_sat2toa.o \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_geo_tio.f90 -o $(ODIR)/m_read_geo_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_metadata_tio.o: $(SDIR)/m_read_metadata_tio.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/utilities.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_metadata_tio.f90 -o $(ODIR)/m_read_metadata_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_o3p_params.o: $(SDIR)/m_o3p_params.f90 \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_o3p_params.f90 -o $(ODIR)/m_o3p_params.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_L2_o3p_tio.o: $(SDIR)/m_read_L2_o3p_tio.f90 \
  $(ODIR)/m_o3p_params.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_L2_o3p_tio.f90 -o $(ODIR)/m_read_L2_o3p_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/HE5Tget_size.o:	$(SDIR)/HE5Tget_size.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/HE5Tget_size.c -o $(ODIR)/HE5Tget_size.o
$(ODIR)/OMI_L1B_reader.o:	$(SDIR)/OMI_L1B_reader.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_L1B_reader.c -o $(ODIR)/OMI_L1B_reader.o
$(ODIR)/OMI_localGranuleID.o:	$(SDIR)/OMI_localGranuleID.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_localGranuleID.c -o $(ODIR)/OMI_localGranuleID.o
$(ODIR)/OMI_SMF_setmsg.o:	$(SDIR)/OMI_SMF_setmsg.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_SMF_setmsg.c -o $(ODIR)/OMI_SMF_setmsg.o
$(ODIR)/unpack.o:	$(SDIR)/unpack.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/unpack.c -o $(ODIR)/unpack.o
