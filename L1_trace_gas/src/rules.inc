$(ODIR)/slitfunction.o: $(SDIR)/slitfunction.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/slitfunction_asym_gaussian.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/slitfunction.f90 -o $(ODIR)/slitfunction.o $(MODOUTFLG) $(ODIR)
$(ODIR)/slitfunction_omi.o: $(SDIR)/slitfunction_omi.f90 \
  $(ODIR)/arrayutils.o \
  $(ODIR)/ezspline_interpolation.o \
  $(ODIR)/integration_routines.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/sao_pge_utils.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/slitfunction_omi.f90 -o $(ODIR)/slitfunction_omi.o $(MODOUTFLG) $(ODIR)
$(ODIR)/slitfunction_asym_gaussian.o: $(SDIR)/slitfunction_asym_gaussian.f90 \
  $(ODIR)/integration_routines.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/SLATEC_davint.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/slitfunction_asym_gaussian.f90 -o $(ODIR)/slitfunction_asym_gaussian.o $(MODOUTFLG) $(ODIR)
$(ODIR)/cache_module.o: $(SDIR)/cache_module.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/cache_module.f90 -o $(ODIR)/cache_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/optimizer_interface_module.o: $(SDIR)/optimizer_interface_module.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/optimizer_interface_module.f90 -o $(ODIR)/optimizer_interface_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/elsunc_interface_module.o: $(SDIR)/elsunc_interface_module.f90 \
  $(ODIR)/OMSAO_elsunc_fitting_module.o \
  $(ODIR)/optimizer_interface_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/elsunc_interface_module.f90 -o $(ODIR)/elsunc_interface_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/angle_sat2toa.o: $(SDIR)/angle_sat2toa.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/angle_sat2toa.f90 -o $(ODIR)/angle_sat2toa.o $(MODOUTFLG) $(ODIR)
$(ODIR)/dataspline_module.o: $(SDIR)/dataspline_module.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/slitfunction.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/dataspline_module.f90 -o $(ODIR)/dataspline_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/ezspline_interpolation.o: $(SDIR)/ezspline_interpolation.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/ezspline_interpolation.f90 -o $(ODIR)/ezspline_interpolation.o $(MODOUTFLG) $(ODIR)
$(ODIR)/fitting_loops.o: $(SDIR)/fitting_loops.f90 \
  $(ODIR)/cache_module.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/he5_output_tools.o \
  $(ODIR)/irradiance_data.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_radiance_ref_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/output_tools.o \
  $(ODIR)/prepare_databases.o \
  $(ODIR)/radiance_fit.o \
  $(ODIR)/radiance_wavcal.o \
  $(ODIR)/sao_pge_utils.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/fitting_loops.f90 -o $(ODIR)/fitting_loops.o $(MODOUTFLG) $(ODIR)
$(ODIR)/he5_output_tools.o: $(SDIR)/he5_output_tools.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/metadata_tools.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/he5_output_tools.f90 -o $(ODIR)/he5_output_tools.o $(MODOUTFLG) $(ODIR)
$(ODIR)/output_tools.o: $(SDIR)/output_tools.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/tg_names_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/output_tools.f90 -o $(ODIR)/output_tools.o $(MODOUTFLG) $(ODIR)
$(ODIR)/integration_routines.o: $(SDIR)/integration_routines.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/integration_routines.f90 -o $(ODIR)/integration_routines.o $(MODOUTFLG) $(ODIR)
$(ODIR)/metadata_tools.o: $(SDIR)/metadata_tools.f90 \
  $(ODIR)/OMSAO_casestring_module.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/strutils.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/metadata_tools.f90 -o $(ODIR)/metadata_tools.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_pge_fitting_aux.o: $(SDIR)/omi_pge_fitting_aux.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/metadata_tools.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_pge_fitting_aux.f90 -o $(ODIR)/omi_pge_fitting_aux.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_pge_fitting_process.o: $(SDIR)/omi_pge_fitting_process.f90 \
  $(ODIR)/commonmode.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/fitting_loops.o \
  $(ODIR)/he5_output_tools.o \
  $(ODIR)/irradiance_data.o \
  $(ODIR)/l1bread_utils.o \
  $(ODIR)/metadata_tools.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/omi_pge_postprocessing.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_radiance_ref_module.o \
  $(ODIR)/OMSAO_solar_wavcal_module.o \
  $(ODIR)/OMSAO_solcomp_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/OMSAO_wfamf_module.o \
  $(ODIR)/output_tools.o \
  $(ODIR)/swathline_loop.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_pge_fitting_process.f90 -o $(ODIR)/omi_pge_fitting_process.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_pge_postprocessing.o: $(SDIR)/omi_pge_postprocessing.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/he5_output_tools.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_destriping_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_pixelcorner_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_Reference_sector_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/OMSAO_wfamf_module.o \
  $(ODIR)/output_tools.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_pge_postprocessing.f90 -o $(ODIR)/omi_pge_postprocessing.o $(MODOUTFLG) $(ODIR)
$(ODIR)/swathline_loop.o: $(SDIR)/swathline_loop.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/fitting_loops.o \
  $(ODIR)/he5_output_tools.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_radiance_ref_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/output_tools.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/swathline_loop.f90 -o $(ODIR)/swathline_loop.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_parameters_module.o: $(SDIR)/omi_parameters_module.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_parameters_module.f90 -o $(ODIR)/omi_parameters_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/omi_read_l1b_data.o: $(SDIR)/omi_read_l1b_data.f90 \
  $(ODIR)/angle_sat2toa.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/strutils.o \
  $(ODIR)/tg_names_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/omi_read_l1b_data.f90 -o $(ODIR)/omi_read_l1b_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_casestring_module.o: $(SDIR)/OMSAO_casestring_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_casestring_module.f90 -o $(ODIR)/OMSAO_casestring_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_destriping_module.o: $(SDIR)/OMSAO_destriping_module.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_median_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o \
  $(ODIR)/SLATEC_davint.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_destriping_module.f90 -o $(ODIR)/OMSAO_destriping_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_elsunc_fitting_module.o: $(SDIR)/OMSAO_elsunc_fitting_module.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_elsunc_fitting_module.f90 -o $(ODIR)/OMSAO_elsunc_fitting_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_errstat_module.o: $(SDIR)/OMSAO_errstat_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_errstat_module.f90 -o $(ODIR)/OMSAO_errstat_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_he5_module.o: $(SDIR)/OMSAO_he5_module.f90 \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_he5_module.f90 -o $(ODIR)/OMSAO_he5_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_indices_module.o: $(SDIR)/OMSAO_indices_module.f90 \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_indices_module.f90 -o $(ODIR)/OMSAO_indices_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_linterpolation_module.o: $(SDIR)/OMSAO_linterpolation_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_linterpolation_module.f90 -o $(ODIR)/OMSAO_linterpolation_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_main.o: $(SDIR)/OMSAO_main.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/elsunc_interface_module.o \
  $(ODIR)/irradiance_data.o \
  $(ODIR)/metadata_tools.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/omi_pge_fitting_process.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o \
  $(ODIR)/pcf_file_module.o \
  $(ODIR)/read_reference_spectra.o \
  $(ODIR)/slitfunction.o \
  $(ODIR)/slitfunction_omi.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_main.f90 -o $(ODIR)/OMSAO_main.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_median_module.o: $(SDIR)/OMSAO_median_module.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_median_module.f90 -o $(ODIR)/OMSAO_median_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_omidata_module.o: $(SDIR)/OMSAO_omidata_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_omidata_module.f90 -o $(ODIR)/OMSAO_omidata_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_parameters_module.o: $(SDIR)/OMSAO_parameters_module.f90 \
  $(ODIR)/omi_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_parameters_module.f90 -o $(ODIR)/OMSAO_parameters_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_pixelcorner_module.o: $(SDIR)/OMSAO_pixelcorner_module.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_pixelcorner_module.f90 -o $(ODIR)/OMSAO_pixelcorner_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_precision_module.o: $(SDIR)/OMSAO_precision_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_precision_module.f90 -o $(ODIR)/OMSAO_precision_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_prefitcol_module.o: $(SDIR)/OMSAO_prefitcol_module.f90 \
  $(ODIR)/datafields.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_prefitcol_module.f90 -o $(ODIR)/OMSAO_prefitcol_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_radiance_ref_module.o: $(SDIR)/OMSAO_radiance_ref_module.f90 \
  $(ODIR)/arrayutils.o \
  $(ODIR)/cache_module.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/irradiance_data.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_median_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_solcomp_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/radiance_fit.o \
  $(ODIR)/SLATEC_davint.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_radiance_ref_module.f90 -o $(ODIR)/OMSAO_radiance_ref_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_Reference_sector_module.o: $(SDIR)/OMSAO_Reference_sector_module.f90 \
  $(ODIR)/commonmode.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/ezspline_interpolation.o \
  $(ODIR)/l1bread_utils.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/omi_read_l1b_data.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/OMSAO_wfamf_module.o \
  $(ODIR)/output_tools.o \
  $(ODIR)/swathline_loop.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_Reference_sector_module.f90 -o $(ODIR)/OMSAO_Reference_sector_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_solar_wavcal_module.o: $(SDIR)/OMSAO_solar_wavcal_module.f90 \
  $(ODIR)/cache_module.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/ezspline_interpolation.o \
  $(ODIR)/irradiance_data.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o \
  $(ODIR)/wavecal.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_solar_wavcal_module.f90 -o $(ODIR)/OMSAO_solar_wavcal_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_solcomp_module.o: $(SDIR)/OMSAO_solcomp_module.f90 \
  $(ODIR)/arrayutils.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_solcomp_module.f90 -o $(ODIR)/OMSAO_solcomp_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_variables_module.o: $(SDIR)/OMSAO_variables_module.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_variables_module.f90 -o $(ODIR)/OMSAO_variables_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMSAO_wfamf_module.o: $(SDIR)/OMSAO_wfamf_module.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/ezspline_interpolation.o \
  $(ODIR)/omi_pge_fitting_aux.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_linterpolation_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/output_tools.o \
  $(ODIR)/sao_pge_utils.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/OMSAO_wfamf_module.f90 -o $(ODIR)/OMSAO_wfamf_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/prepare_databases.o: $(SDIR)/prepare_databases.f90 \
  $(ODIR)/arrayutils.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/dataspline_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/subtract_cubic.o \
  $(ODIR)/undersample.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/prepare_databases.f90 -o $(ODIR)/prepare_databases.o $(MODOUTFLG) $(ODIR)
$(ODIR)/radiance_fit.o: $(SDIR)/radiance_fit.f90 \
  $(ODIR)/arrayutils.o \
  $(ODIR)/cache_module.o \
  $(ODIR)/commonmode.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_solcomp_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/subtract_cubic.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/radiance_fit.f90 -o $(ODIR)/radiance_fit.o $(MODOUTFLG) $(ODIR)
$(ODIR)/radiance_wavcal.o: $(SDIR)/radiance_wavcal.f90 \
  $(ODIR)/commonmode.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o \
  $(ODIR)/wavecal.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/radiance_wavcal.f90 -o $(ODIR)/radiance_wavcal.o $(MODOUTFLG) $(ODIR)
$(ODIR)/control_module.o: $(SDIR)/control_module.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/OMSAO_casestring_module.o \
  $(ODIR)/OMSAO_destriping_module.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/OMSAO_wfamf_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/strutils.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/control_module.f90 -o $(ODIR)/control_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/pcf_file_module.o: $(SDIR)/pcf_file_module.f90 \
  $(ODIR)/control_module.o \
  $(ODIR)/metadata_tools.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_prefitcol_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/OMSAO_wfamf_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/strutils.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/pcf_file_module.f90 -o $(ODIR)/pcf_file_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/read_reference_spectra.o: $(SDIR)/read_reference_spectra.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/datafields.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_solcomp_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/read_reference_spectra.f90 -o $(ODIR)/read_reference_spectra.o $(MODOUTFLG) $(ODIR)
$(ODIR)/sao_pge_utils.o: $(SDIR)/sao_pge_utils.f90 \
  $(ODIR)/arrayutils.o \
  $(ODIR)/ezspline_interpolation.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/sao_pge_utils.f90 -o $(ODIR)/sao_pge_utils.o $(MODOUTFLG) $(ODIR)
$(ODIR)/strutils.o: $(SDIR)/strutils.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/strutils.f90 -o $(ODIR)/strutils.o $(MODOUTFLG) $(ODIR)
$(ODIR)/subtract_cubic.o: $(SDIR)/subtract_cubic.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/subtract_cubic.f90 -o $(ODIR)/subtract_cubic.o $(MODOUTFLG) $(ODIR)
$(ODIR)/SLATEC_davint.o: $(SDIR)/SLATEC_davint.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/SLATEC_davint.f90 -o $(ODIR)/SLATEC_davint.o $(MODOUTFLG) $(ODIR)
$(ODIR)/tg_names_module.o: $(SDIR)/tg_names_module.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/tg_names_module.f90 -o $(ODIR)/tg_names_module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/undersample.o: $(SDIR)/undersample.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/slitfunction.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/undersample.f90 -o $(ODIR)/undersample.o $(MODOUTFLG) $(ODIR)
$(ODIR)/datafields.o: $(SDIR)/datafields.f90 \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/OMSAO_he5_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/datafields.f90 -o $(ODIR)/datafields.o $(MODOUTFLG) $(ODIR)
$(ODIR)/irradiance_data.o: $(SDIR)/irradiance_data.f90 \
  $(ODIR)/arrayutils.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/l1bread_utils.o \
  $(ODIR)/OMSAO_errstat_module.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_solcomp_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/sao_pge_utils.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/irradiance_data.f90 -o $(ODIR)/irradiance_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/l1bread_utils.o: $(SDIR)/l1bread_utils.f90 \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/l1bread_utils.f90 -o $(ODIR)/l1bread_utils.o $(MODOUTFLG) $(ODIR)
$(ODIR)/wavecal.o: $(SDIR)/wavecal.f90 \
  $(ODIR)/cache_module.o \
  $(ODIR)/ctrlvars.o \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_parameters_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o \
  $(ODIR)/optimizer_interface_module.o \
  $(ODIR)/sao_pge_utils.o \
  $(ODIR)/slitfunction.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/wavecal.f90 -o $(ODIR)/wavecal.o $(MODOUTFLG) $(ODIR)
$(ODIR)/commonmode.o: $(SDIR)/commonmode.f90 \
  $(ODIR)/OMSAO_indices_module.o \
  $(ODIR)/OMSAO_omidata_module.o \
  $(ODIR)/OMSAO_precision_module.o \
  $(ODIR)/OMSAO_variables_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/commonmode.f90 -o $(ODIR)/commonmode.o $(MODOUTFLG) $(ODIR)
$(ODIR)/ctrlvars.o: $(SDIR)/ctrlvars.f90
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/ctrlvars.f90 -o $(ODIR)/ctrlvars.o $(MODOUTFLG) $(ODIR)
$(ODIR)/arrayutils.o: $(SDIR)/arrayutils.f90 \
  $(ODIR)/OMSAO_precision_module.o
	$(FC) -c $(FCFLAGS) $(IFLAGS) $(SDIR)/arrayutils.f90 -o $(ODIR)/arrayutils.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_SMF_setmsg.o:	$(SDIR)/OMI_SMF_setmsg.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(OMI_SMF_setmsg_INCS) $(OMI_SMF_setmsg_DEFINES) $(SDIR)/OMI_SMF_setmsg.c -o $(ODIR)/OMI_SMF_setmsg.o
$(ODIR)/unbufferSTDOUT.o:	$(SDIR)/unbufferSTDOUT.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(unbufferSTDOUT_INCS) $(unbufferSTDOUT_DEFINES) $(SDIR)/unbufferSTDOUT.c -o $(ODIR)/unbufferSTDOUT.o
