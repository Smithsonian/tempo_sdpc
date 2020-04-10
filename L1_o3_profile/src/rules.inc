$(ODIR)/OMI_SMF_class.o: $(SDIR)/OMI_SMF_class.f90 \
  $(ODIR)/OMSAO_errstat_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMI_SMF_class.f90 -o $(ODIR)/OMI_SMF_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_errstat_module.o: $(SDIR)/OMSAO_errstat_module.f90 \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_errstat_module.f90 -o $(ODIR)/OMSAO_errstat_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_precision_module.o: $(SDIR)/OMSAO_precision_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_precision_module.f90 -o $(ODIR)/OMSAO_precision_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_indices_module.o: $(SDIR)/OMSAO_indices_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_indices_module.f90 -o $(ODIR)/OMSAO_indices_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_parameters_module.o: $(SDIR)/OMSAO_parameters_module.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_parameters_module.f90 -o $(ODIR)/OMSAO_parameters_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_variables_module.o: $(SDIR)/OMSAO_variables_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_variables_module.f90 -o $(ODIR)/OMSAO_variables_module.o $(MODOUTFLG) $(ODIR)
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
$(ODIR)/OMSAO_he5_module.o: $(SDIR)/OMSAO_he5_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_he5_module.f90 -o $(ODIR)/OMSAO_he5_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_pixelcorner_module.o: $(SDIR)/OMSAO_pixelcorner_module.f90 \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/m_angle_sat2toa.o \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_pixelcorner_module.f90 -o $(ODIR)/OMSAO_pixelcorner_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_slitfunction_module.o: $(SDIR)/OMSAO_slitfunction_module.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_slitfunction_module.f90 -o $(ODIR)/OMSAO_slitfunction_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_omidata_module.o: $(SDIR)/OMSAO_omidata_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_omidata_module.f90 -o $(ODIR)/OMSAO_omidata_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_gome_data_module.o: $(SDIR)/OMSAO_gome_data_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_gome_data_module.f90 -o $(ODIR)/OMSAO_gome_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_tmpodata_module.o: $(SDIR)/OMSAO_tmpodata_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_tmpodata_module.f90 -o $(ODIR)/OMSAO_tmpodata_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/ozprof_data_module.o: $(SDIR)/ozprof_data_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/ozprof_data_module.f90 -o $(ODIR)/ozprof_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/cloud_data_module.o: $(SDIR)/cloud_data_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/cloud_data_module.f90 -o $(ODIR)/cloud_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_fitting_controls.o: $(SDIR)/m_read_fitting_controls.f90 \
  $(ODIR)/m_read_ozprof_input.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_tmpodata_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/PGS_PC_class.o \
  $(ODIR)/UTIL_tools_class.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_fitting_controls.f90 -o $(ODIR)/m_read_fitting_controls.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_ozprof_input.o: $(SDIR)/m_read_ozprof_input.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_ozprof_input.f90 -o $(ODIR)/m_read_ozprof_input.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_reference_spectra.o: $(SDIR)/m_read_reference_spectra.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_reference_spectra.f90 -o $(ODIR)/m_read_reference_spectra.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_allocate.o: $(SDIR)/m_allocate.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_allocate.f90 -o $(ODIR)/m_allocate.o $(MODOUTFLG) $(ODIR)
$(ODIR)/tmpo_read_l1b_data.o: $(SDIR)/tmpo_read_l1b_data.f90 \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/m_fitting_util.o \
  $(ODIR)/m_gauss.o \
  $(ODIR)/m_read_l1_tio.o \
  $(ODIR)/m_read_metadata_tio.o \
  $(ODIR)/m_triangle.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_tmpodata_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/tmpo_read_l1b_data.f90 -o $(ODIR)/tmpo_read_l1b_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_bounded_nonlin_ls.o: $(SDIR)/m_bounded_nonlin_ls.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_bounded_nonlin_ls.f90 -o $(ODIR)/m_bounded_nonlin_ls.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_ezspline_interpolation.o: $(SDIR)/m_ezspline_interpolation.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_ezspline_interpolation.f90 -o $(ODIR)/m_ezspline_interpolation.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_utilities.o: $(SDIR)/m_utilities.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_utilities.f90 -o $(ODIR)/m_utilities.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_convert_coadd.o: $(SDIR)/m_convert_coadd.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_convert_coadd.f90 -o $(ODIR)/m_convert_coadd.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_angle_sat2toa.o: $(SDIR)/m_angle_sat2toa.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_angle_sat2toa.f90 -o $(ODIR)/m_angle_sat2toa.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_super_gauss.o: $(SDIR)/m_super_gauss.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_super_gauss.f90 -o $(ODIR)/m_super_gauss.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_gauss.o: $(SDIR)/m_gauss.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_gauss.f90 -o $(ODIR)/m_gauss.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_triangle.o: $(SDIR)/m_triangle.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_triangle.f90 -o $(ODIR)/m_triangle.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_voigt.o: $(SDIR)/m_voigt.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_voigt.f90 -o $(ODIR)/m_voigt.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_convol.o: $(SDIR)/m_convol.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_gauss.o \
  $(ODIR)/m_super_gauss.o \
  $(ODIR)/m_triangle.o \
  $(ODIR)/m_voigt.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_convol.f90 -o $(ODIR)/m_convol.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_avg_band.o: $(SDIR)/m_avg_band.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_avg_band.f90 -o $(ODIR)/m_avg_band.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_radiance_fit.o: $(SDIR)/m_radiance_fit.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_subtract_poly.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_radiance_fit.f90 -o $(ODIR)/m_radiance_fit.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_radiance_wavcal.o: $(SDIR)/m_radiance_wavcal.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_subtract_poly.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_radiance_wavcal.f90 -o $(ODIR)/m_radiance_wavcal.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_rad_shisqu.o: $(SDIR)/m_rad_shisqu.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_rad_shisqu.f90 -o $(ODIR)/m_rad_shisqu.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_solar_fit.o: $(SDIR)/m_solar_fit.f90 \
  $(ODIR)/m_cal_fit_one.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_fitting_util.o \
  $(ODIR)/m_subtract_poly.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
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
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_solar_wavcal_vary.f90 -o $(ODIR)/m_solar_wavcal_vary.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_subtract_poly.o: $(SDIR)/m_subtract_poly.f90 \
  $(ODIR)/m_bounded_nonlin_ls.o \
  $(ODIR)/m_fitting_util.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_subtract_poly.f90 -o $(ODIR)/m_subtract_poly.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_subtract_cubic.o: $(SDIR)/m_subtract_cubic.f90 \
  $(ODIR)/m_bounded_nonlin_ls.o \
  $(ODIR)/m_fitting_util.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_subtract_cubic.f90 -o $(ODIR)/m_subtract_cubic.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_spectra.o: $(SDIR)/m_spectra.f90 \
  $(ODIR)/m_convol.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_slitfunction_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_spectra.f90 -o $(ODIR)/m_spectra.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_fitting_util.o: $(SDIR)/m_fitting_util.f90 \
  $(ODIR)/m_bounded_nonlin_ls.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_fitting_util.f90 -o $(ODIR)/m_fitting_util.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_cal_fit_one.o: $(SDIR)/m_cal_fit_one.f90 \
  $(ODIR)/m_bounded_nonlin_ls.o \
  $(ODIR)/m_fitting_util.o \
  $(ODIR)/m_spectra.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_cal_fit_one.f90 -o $(ODIR)/m_cal_fit_one.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_prepare_databases.o: $(SDIR)/m_prepare_databases.f90 \
  $(ODIR)/m_avg_band.o \
  $(ODIR)/m_convol.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_prepare_databases.f90 -o $(ODIR)/m_prepare_databases.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_cross_calibrate.o: $(SDIR)/m_cross_calibrate.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_radiance_fit.o \
  $(ODIR)/m_radiance_wavcal.o \
  $(ODIR)/m_solar_fit.o \
  $(ODIR)/m_solar_wavcal_vary.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_cross_calibrate.f90 -o $(ODIR)/m_cross_calibrate.o $(MODOUTFLG) $(ODIR)
$(ODIR)/tmpo_adj_data.o: $(SDIR)/tmpo_adj_data.f90 \
  $(ODIR)/m_angle_sat2toa.o \
  $(ODIR)/m_avg_band.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_fitting_util.o \
  $(ODIR)/m_get_cloud.o \
  $(ODIR)/m_prepare_databases.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_tmpodata_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/tmpo_adj_data.f90 -o $(ODIR)/tmpo_adj_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/tmpo_pge_process.o: $(SDIR)/tmpo_pge_process.f90 \
  $(ODIR)/m_allocate.o \
  $(ODIR)/m_cross_calibrate.o \
  $(ODIR)/m_read_cloud_tio.o \
  $(ODIR)/m_read_geo_tio.o \
  $(ODIR)/m_specfit_ozprof.o \
  $(ODIR)/O3P_output_module.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_tmpodata_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/tio_output_module.o \
  $(ODIR)/tmpo_adj_data.o \
  $(ODIR)/tmpo_read_l1b_data.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/tmpo_pge_process.f90 -o $(ODIR)/tmpo_pge_process.o $(MODOUTFLG) $(ODIR)
$(ODIR)/gsvd_data_module.o: $(SDIR)/gsvd_data_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/gsvd_data_module.f90 -o $(ODIR)/gsvd_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_gsvd_o3prof_utilities.o: $(SDIR)/m_gsvd_o3prof_utilities.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_gsvd_o3prof_utilities.f90 -o $(ODIR)/m_gsvd_o3prof_utilities.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_gsvd_lcurve_gcv.o: $(SDIR)/m_gsvd_lcurve_gcv.f90 \
  $(ODIR)/gsvd_data_module.o \
  $(ODIR)/m_gsvd_o3prof_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_gsvd_lcurve_gcv.f90 -o $(ODIR)/m_gsvd_lcurve_gcv.o $(MODOUTFLG) $(ODIR)
$(ODIR)/hdfeos4_parameters.o: $(SDIR)/hdfeos4_parameters.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/hdfeos4_parameters.f90 -o $(ODIR)/hdfeos4_parameters.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_input_tools.o: $(SDIR)/he5_input_tools.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_input_tools.f90 -o $(ODIR)/he5_input_tools.o $(MODOUTFLG) $(ODIR)
$(ODIR)/L1B_Reader_class.o: $(SDIR)/L1B_Reader_class.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/m_pgs_include.o \
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
$(ODIR)/m_get_cloud.o: $(SDIR)/m_get_cloud.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_cloud.f90 -o $(ODIR)/m_get_cloud.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_ncep.o: $(SDIR)/m_get_ncep.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_ncep.f90 -o $(ODIR)/m_get_ncep.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_fnl.o: $(SDIR)/m_get_fnl.f90 \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_fnl.f90 -o $(ODIR)/m_get_fnl.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_met_tempo.o: $(SDIR)/m_get_met_tempo.f90 \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_met_tempo.f90 -o $(ODIR)/m_get_met_tempo.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_omigeos5.o: $(SDIR)/m_get_omigeos5.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_omigeos5.f90 -o $(ODIR)/m_get_omigeos5.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_bclayer.o: $(SDIR)/m_get_bclayer.f90 \
  $(ODIR)/m_get_fnl.o \
  $(ODIR)/m_get_met_tempo.o \
  $(ODIR)/m_get_ncep.o \
  $(ODIR)/m_get_omigeos5.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_bclayer.f90 -o $(ODIR)/m_get_bclayer.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_toz.o: $(SDIR)/m_get_toz.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_toz.f90 -o $(ODIR)/m_get_toz.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_o3prof.o: $(SDIR)/m_get_o3prof.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_o3prof.f90 -o $(ODIR)/m_get_o3prof.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_tracegas.o: $(SDIR)/m_get_tracegas.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_tracegas.f90 -o $(ODIR)/m_get_tracegas.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_xcrs.o: $(SDIR)/m_get_xcrs.f90 \
  $(ODIR)/m_avg_band.o \
  $(ODIR)/m_convol.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_xcrs.f90 -o $(ODIR)/m_get_xcrs.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_initial_albedo.o: $(SDIR)/m_get_initial_albedo.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_lidort_env_vv2p7.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_initial_albedo.f90 -o $(ODIR)/m_get_initial_albedo.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_reg_matrix.o: $(SDIR)/m_get_reg_matrix.f90 \
  $(ODIR)/m_gsvd_lcurve_gcv.o \
  $(ODIR)/m_gsvd_o3prof_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_reg_matrix.f90 -o $(ODIR)/m_get_reg_matrix.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_prepare_atmosphere.o: $(SDIR)/m_prepare_atmosphere.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_prepare_atmosphere.f90 -o $(ODIR)/m_prepare_atmosphere.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_aerosol.o: $(SDIR)/m_get_aerosol.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_aerosol.f90 -o $(ODIR)/m_get_aerosol.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_raman.o: $(SDIR)/m_raman.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_raman.f90 -o $(ODIR)/m_raman.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_raman.o: $(SDIR)/m_get_raman.f90 \
  $(ODIR)/m_avg_band.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_xcrs.o \
  $(ODIR)/m_raman.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_get_raman.f90 -o $(ODIR)/m_get_raman.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_make_atm.o: $(SDIR)/m_make_atm.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_aerosol.o \
  $(ODIR)/m_get_cloud.o \
  $(ODIR)/m_get_fnl.o \
  $(ODIR)/m_get_met_tempo.o \
  $(ODIR)/m_get_ncep.o \
  $(ODIR)/m_get_o3prof.o \
  $(ODIR)/m_get_omigeos5.o \
  $(ODIR)/m_get_tracegas.o \
  $(ODIR)/m_prepare_atmosphere.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_make_atm.f90 -o $(ODIR)/m_make_atm.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_set_tracegas.o: $(SDIR)/m_set_tracegas.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_set_tracegas.f90 -o $(ODIR)/m_set_tracegas.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_set_brdf.o: $(SDIR)/m_set_brdf.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_input_module.o \
  $(ODIR)/brdf_level1_def.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/brdf_surface_module.o \
  $(ODIR)/brdf_window_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_set_brdf.f90 -o $(ODIR)/m_set_brdf.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_parameters_module.o: $(SDIR)/brdf_parameters_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_parameters_module.f90 -o $(ODIR)/brdf_parameters_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_error_module.o: $(SDIR)/brdf_error_module.f90 \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_error_module.f90 -o $(ODIR)/brdf_error_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_netcdf_module.o: $(SDIR)/brdf_netcdf_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_netcdf_module.f90 -o $(ODIR)/brdf_netcdf_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_profilepar_module.o: $(SDIR)/brdf_profilepar_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_profilepar_module.f90 -o $(ODIR)/brdf_profilepar_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_covariance_module.o: $(SDIR)/brdf_covariance_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_covariance_module.f90 -o $(ODIR)/brdf_covariance_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_time_module.o: $(SDIR)/brdf_time_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_time_module.f90 -o $(ODIR)/brdf_time_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_level1_def.o: $(SDIR)/brdf_level1_def.f90 \
  $(ODIR)/brdf_covariance_module.o \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o \
  $(ODIR)/brdf_time_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_level1_def.f90 -o $(ODIR)/brdf_level1_def.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_interpolation_module.o: $(SDIR)/brdf_interpolation_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_level1_def.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_interpolation_module.f90 -o $(ODIR)/brdf_interpolation_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_string_utils_module.o: $(SDIR)/brdf_string_utils_module.f90 \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_string_utils_module.f90 -o $(ODIR)/brdf_string_utils_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_profile_module.o: $(SDIR)/brdf_profile_module.f90 \
  $(ODIR)/brdf_covariance_module.o \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_interpolation_module.o \
  $(ODIR)/brdf_level1_def.o \
  $(ODIR)/brdf_netcdf_module.o \
  $(ODIR)/brdf_parameters_module.o \
  $(ODIR)/brdf_profilepar_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_profile_module.f90 -o $(ODIR)/brdf_profile_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_rtm_def_module.o: $(SDIR)/brdf_rtm_def_module.f90 \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_rtm_def_module.f90 -o $(ODIR)/brdf_rtm_def_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_pca_bin_module.o: $(SDIR)/brdf_pca_bin_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o \
  $(ODIR)/brdf_string_utils_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_pca_bin_module.f90 -o $(ODIR)/brdf_pca_bin_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_isrf_module.o: $(SDIR)/brdf_isrf_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_interpolation_module.o \
  $(ODIR)/brdf_netcdf_module.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_isrf_module.f90 -o $(ODIR)/brdf_isrf_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_window_module.o: $(SDIR)/brdf_window_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_interpolation_module.o \
  $(ODIR)/brdf_isrf_module.o \
  $(ODIR)/brdf_level1_def.o \
  $(ODIR)/brdf_netcdf_module.o \
  $(ODIR)/brdf_parameters_module.o \
  $(ODIR)/brdf_pca_bin_module.o \
  $(ODIR)/brdf_rtm_def_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_window_module.f90 -o $(ODIR)/brdf_window_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_spatial_fa_module.o: $(SDIR)/brdf_spatial_fa_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_interpolation_module.o \
  $(ODIR)/brdf_level1_def.o \
  $(ODIR)/brdf_netcdf_module.o \
  $(ODIR)/brdf_parameters_module.o \
  $(ODIR)/brdf_profile_module.o \
  $(ODIR)/brdf_time_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_spatial_fa_module.f90 -o $(ODIR)/brdf_spatial_fa_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_chebyshev_module.o: $(SDIR)/brdf_chebyshev_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_chebyshev_module.f90 -o $(ODIR)/brdf_chebyshev_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_sif_module.o: $(SDIR)/brdf_sif_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_interpolation_module.o \
  $(ODIR)/brdf_netcdf_module.o \
  $(ODIR)/brdf_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_sif_module.f90 -o $(ODIR)/brdf_sif_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_surface_module.o: $(SDIR)/brdf_surface_module.f90 \
  $(ODIR)/brdf_chebyshev_module.o \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_interpolation_module.o \
  $(ODIR)/brdf_level1_def.o \
  $(ODIR)/brdf_netcdf_module.o \
  $(ODIR)/brdf_parameters_module.o \
  $(ODIR)/brdf_profile_module.o \
  $(ODIR)/brdf_sif_module.o \
  $(ODIR)/brdf_spatial_fa_module.o \
  $(ODIR)/brdf_window_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_surface_module.f90 -o $(ODIR)/brdf_surface_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/brdf_input_module.o: $(SDIR)/brdf_input_module.f90 \
  $(ODIR)/brdf_error_module.o \
  $(ODIR)/brdf_parameters_module.o \
  $(ODIR)/brdf_string_utils_module.o \
  $(ODIR)/brdf_surface_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/brdf_input_module.f90 -o $(ODIR)/brdf_input_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_set_cldalb.o: $(SDIR)/m_set_cldalb.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_initial_albedo.o \
  $(ODIR)/m_set_brdf.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_set_cldalb.f90 -o $(ODIR)/m_set_cldalb.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_spectra_reflectance.o: $(SDIR)/m_spectra_reflectance.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_spectra_reflectance.f90 -o $(ODIR)/m_spectra_reflectance.o $(MODOUTFLG) $(ODIR)
$(ODIR)/vlidort_data_module.o: $(SDIR)/vlidort_data_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/vlidort_data_module.f90 -o $(ODIR)/vlidort_data_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/GEMSTOOL_binset.o: $(SDIR)/GEMSTOOL_binset.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/GEMSTOOL_binset.f90 -o $(ODIR)/GEMSTOOL_binset.o $(MODOUTFLG) $(ODIR)
$(ODIR)/GEMSTOOL_CreateBins.o: $(SDIR)/GEMSTOOL_CreateBins.f90 \
  $(ODIR)/GEMSTOOL_binset.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/GEMSTOOL_CreateBins.f90 -o $(ODIR)/GEMSTOOL_CreateBins.o $(MODOUTFLG) $(ODIR)
$(ODIR)/vlidort_eofpc_module.o: $(SDIR)/vlidort_eofpc_module.f90 \
  $(ODIR)/GEMSTOOL_CreateBins.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/vlidort_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/vlidort_eofpc_module.f90 -o $(ODIR)/vlidort_eofpc_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/vlidort_exact_module.o: $(SDIR)/vlidort_exact_module.f90 \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/vlidort_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/vlidort_exact_module.f90 -o $(ODIR)/vlidort_exact_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_lidort_util.o: $(SDIR)/m_lidort_util.f90 \
  $(ODIR)/m_avg_band.o \
  $(ODIR)/m_convol.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_xcrs.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_lidort_util.f90 -o $(ODIR)/m_lidort_util.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_lidort_env_vv2p7_pca.o: $(SDIR)/m_lidort_env_vv2p7_pca.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_xcrs.o \
  $(ODIR)/m_lidort_util.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/vlidort_data_module.o \
  $(ODIR)/vlidort_eofpc_module.o \
  $(ODIR)/vlidort_exact_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_lidort_env_vv2p7_pca.f90 -o $(ODIR)/m_lidort_env_vv2p7_pca.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_lidort_env_vv2p7.o: $(SDIR)/m_lidort_env_vv2p7.f90 \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_xcrs.o \
  $(ODIR)/m_lidort_util.o \
  $(ODIR)/m_set_brdf.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/brdf_surface_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_lidort_env_vv2p7.f90 -o $(ODIR)/m_lidort_env_vv2p7.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_pseudo_model.o: $(SDIR)/m_pseudo_model.f90 \
  $(ODIR)/cloud_data_module.o \
  $(ODIR)/m_fitting_util.o \
  $(ODIR)/m_gsvd_o3prof_utilities.o \
  $(ODIR)/m_lidort_env_vv2p7.o \
  $(ODIR)/m_lidort_env_vv2p7_pca.o \
  $(ODIR)/m_lidort_util.o \
  $(ODIR)/m_set_brdf.o \
  $(ODIR)/m_spectra_reflectance.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_gome_data_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_pseudo_model.f90 -o $(ODIR)/m_pseudo_model.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_oe_inversion.o: $(SDIR)/m_oe_inversion.f90 \
  $(ODIR)/m_gsvd_o3prof_utilities.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_oe_inversion.f90 -o $(ODIR)/m_oe_inversion.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_ozprof_inverse.o: $(SDIR)/m_ozprof_inverse.f90 \
  $(ODIR)/m_bounded_nonlin_ls.o \
  $(ODIR)/m_ezspline_interpolation.o \
  $(ODIR)/m_get_o3prof.o \
  $(ODIR)/m_get_raman.o \
  $(ODIR)/m_gsvd_lcurve_gcv.o \
  $(ODIR)/m_make_atm.o \
  $(ODIR)/m_oe_inversion.o \
  $(ODIR)/m_pseudo_model.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_ozprof_inverse.f90 -o $(ODIR)/m_ozprof_inverse.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_specfit_ozprof.o: $(SDIR)/m_specfit_ozprof.f90 \
  $(ODIR)/m_get_bclayer.o \
  $(ODIR)/m_get_o3prof.o \
  $(ODIR)/m_get_reg_matrix.o \
  $(ODIR)/m_get_toz.o \
  $(ODIR)/m_lidort_util.o \
  $(ODIR)/m_make_atm.o \
  $(ODIR)/m_ozprof_inverse.o \
  $(ODIR)/m_set_brdf.o \
  $(ODIR)/m_set_cldalb.o \
  $(ODIR)/m_set_tracegas.o \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_specfit_ozprof.f90 -o $(ODIR)/m_specfit_ozprof.o $(MODOUTFLG) $(ODIR)
$(ODIR)/ascii_output_module.o: $(SDIR)/ascii_output_module.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/ascii_output_module.f90 -o $(ODIR)/ascii_output_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/tio_output_module.o: $(SDIR)/tio_output_module.f90 \
  $(ODIR)/m_o3p_params.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/tio_output_module.f90 -o $(ODIR)/tio_output_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_output_module.o: $(SDIR)/he5_output_module.f90 \
  $(ODIR)/he5_l2_fs.o \
  $(ODIR)/he5_l2writer_class.o \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/L2_attr_class.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_output_module.f90 -o $(ODIR)/he5_output_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/O3P_output_module.o: $(SDIR)/O3P_output_module.f90 \
  $(ODIR)/ascii_output_module.o \
  $(ODIR)/he5_output_module.o \
  $(ODIR)/OMI_LUN_set.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o \
  $(ODIR)/tio_output_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/O3P_output_module.f90 -o $(ODIR)/O3P_output_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_l2writer_class.o: $(SDIR)/he5_l2writer_class.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_omidata_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_l2writer_class.f90 -o $(ODIR)/he5_l2writer_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_l2_fs.o: $(SDIR)/he5_l2_fs.f90 \
  $(ODIR)/OMSAO_he5_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_l2_fs.f90 -o $(ODIR)/he5_l2_fs.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_LUN_set.o: $(SDIR)/OMI_LUN_set.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMI_LUN_set.f90 -o $(ODIR)/OMI_LUN_set.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_MET_class.o: $(SDIR)/PGS_MET_class.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/PGS_MET_class.f90 -o $(ODIR)/PGS_MET_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_PC_class.o: $(SDIR)/PGS_PC_class.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/PGS_PC_class.f90 -o $(ODIR)/PGS_PC_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/PGS_TD_class.o: $(SDIR)/PGS_TD_class.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/PGS_TD_class.f90 -o $(ODIR)/PGS_TD_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/Szoom_Parameter_Module.o: $(SDIR)/Szoom_Parameter_Module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/Szoom_Parameter_Module.f90 -o $(ODIR)/Szoom_Parameter_Module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/UTIL_tools_class.o: $(SDIR)/UTIL_tools_class.f90 \
  $(ODIR)/OMSAO_errstat_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/UTIL_tools_class.f90 -o $(ODIR)/UTIL_tools_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/o3p_names_module.o: $(SDIR)/o3p_names_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/o3p_names_module.f90 -o $(ODIR)/o3p_names_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_cloud_tio.o: $(SDIR)/m_read_cloud_tio.f90 \
  $(ODIR)/m_convert_coadd.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_tmpodata_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_cloud_tio.f90 -o $(ODIR)/m_read_cloud_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_l1_tio.o: $(SDIR)/m_read_l1_tio.f90 \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_tmpodata_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_l1_tio.f90 -o $(ODIR)/m_read_l1_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_geo_tio.o: $(SDIR)/m_read_geo_tio.f90 \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_tmpodata_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_geo_tio.f90 -o $(ODIR)/m_read_geo_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_metadata_tio.o: $(SDIR)/m_read_metadata_tio.f90 \
  $(ODIR)/m_utilities.o \
  $(ODIR)/OMSAO_parameters_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_metadata_tio.f90 -o $(ODIR)/m_read_metadata_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_L2_o3p_tio.o: $(SDIR)/m_read_L2_o3p_tio.f90 \
  $(ODIR)/m_o3p_params.o \
  $(ODIR)/o3p_names_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_read_L2_o3p_tio.f90 -o $(ODIR)/m_read_L2_o3p_tio.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_o3p_params.o: $(SDIR)/m_o3p_params.f90 \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/ozprof_data_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_o3p_params.f90 -o $(ODIR)/m_o3p_params.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_pgs_include.o: $(SDIR)/m_pgs_include.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/m_pgs_include.f90 -o $(ODIR)/m_pgs_include.o $(MODOUTFLG) $(ODIR)
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
