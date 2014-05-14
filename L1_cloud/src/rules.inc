$(ODIR)/L1B_Reader_class.o: $(SDIR)/L1B_Reader_class.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/Szoom_Parameter_Module.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/L1B_Reader_class.f90 -o $(ODIR)/L1B_Reader_class.o $(MODOUTFLG) $(ODIR)
$(ODIR)/MetadatOBJModule.o: $(SDIR)/MetadatOBJModule.f90 \
  $(ODIR)/m_pgs_include.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/MetadatOBJModule.f90 -o $(ODIR)/MetadatOBJModule.o $(MODOUTFLG) $(ODIR)
$(ODIR)/Szoom_Parameter_Module.o: $(SDIR)/Szoom_Parameter_Module.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/Szoom_Parameter_Module.f90 -o $(ODIR)/Szoom_Parameter_Module.o $(MODOUTFLG) $(ODIR)
$(ODIR)/hdfeos4_parameters.o: $(SDIR)/hdfeos4_parameters.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/hdfeos4_parameters.f90 -o $(ODIR)/hdfeos4_parameters.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_LUN_set.o: $(SDIR)/m_LUN_set.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_LUN_set.f90 -o $(ODIR)/m_LUN_set.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_pgs_include.o: $(SDIR)/m_pgs_include.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_pgs_include.f90 -o $(ODIR)/m_pgs_include.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_alloc1.o: $(SDIR)/m_alloc1.f90 \
  $(ODIR)/m_cloud_pres_mod.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_alloc1.f90 -o $(ODIR)/m_alloc1.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_alloc2.o: $(SDIR)/m_alloc2.f90 \
  $(ODIR)/m_cloud_pres_mod.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_alloc2.f90 -o $(ODIR)/m_alloc2.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_avg.o: $(SDIR)/m_avg.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_avg.f90 -o $(ODIR)/m_avg.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_bilinear.o: $(SDIR)/m_bilinear.f90 \
  $(ODIR)/m_findgen.o \
  $(ODIR)/m_interpol.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_bilinear.f90 -o $(ODIR)/m_bilinear.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_cloud_mask.o: $(SDIR)/m_cloud_mask.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/m_avg.o \
  $(ODIR)/m_eigen.o \
  $(ODIR)/m_find.o \
  $(ODIR)/m_matmul.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_sigma.o \
  $(ODIR)/m_sortind.o \
  $(ODIR)/m_strpos.o \
  $(ODIR)/m_swathnames.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_cloud_mask.f90 -o $(ODIR)/m_cloud_mask.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_cloud_pres_mod.o: $(SDIR)/m_cloud_pres_mod.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_cloud_pres_mod.f90 -o $(ODIR)/m_cloud_pres_mod.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_cloud_pres_ret.o: $(SDIR)/m_cloud_pres_ret.f90 \
  $(ODIR)/m_alloc1.o \
  $(ODIR)/m_alloc2.o \
  $(ODIR)/m_bilinear.o \
  $(ODIR)/m_cloud_pres_mod.o \
  $(ODIR)/m_find.o \
  $(ODIR)/m_findgen.o \
  $(ODIR)/m_get_ai_refl.o \
  $(ODIR)/m_get_f.o \
  $(ODIR)/m_indgen.o \
  $(ODIR)/m_interp_pres.o \
  $(ODIR)/m_interp_ring_rad.o \
  $(ODIR)/m_interpol.o \
  $(ODIR)/m_invert.o \
  $(ODIR)/m_lambda_qual.o \
  $(ODIR)/m_matmul.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_poly_fit.o \
  $(ODIR)/m_print_ret.o \
  $(ODIR)/m_sigma.o \
  $(ODIR)/m_spline.o \
  $(ODIR)/m_trilin.o \
  $(ODIR)/m_vars.o \
  $(ODIR)/mathcons.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_cloud_pres_ret.f90 -o $(ODIR)/m_cloud_pres_ret.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_eigen.o: $(SDIR)/m_eigen.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_eigen.f90 -o $(ODIR)/m_eigen.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_find.o: $(SDIR)/m_find.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_find.f90 -o $(ODIR)/m_find.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_findgen.o: $(SDIR)/m_findgen.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_findgen.f90 -o $(ODIR)/m_findgen.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_ai_refl.o: $(SDIR)/m_get_ai_refl.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_get_ai_refl.f90 -o $(ODIR)/m_get_ai_refl.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_get_f.o: $(SDIR)/m_get_f.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_get_f.f90 -o $(ODIR)/m_get_f.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_indgen.o: $(SDIR)/m_indgen.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_indgen.f90 -o $(ODIR)/m_indgen.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_initialize.o: $(SDIR)/m_initialize.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_initialize.f90 -o $(ODIR)/m_initialize.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_instr_config.o: $(SDIR)/m_instr_config.f90 \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_instr_config.f90 -o $(ODIR)/m_instr_config.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_interp_pres.o: $(SDIR)/m_interp_pres.f90 \
  $(ODIR)/m_cloud_pres_mod.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_interp_pres.f90 -o $(ODIR)/m_interp_pres.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_interp_ring_rad.o: $(SDIR)/m_interp_ring_rad.f90 \
  $(ODIR)/m_cloud_pres_mod.o \
  $(ODIR)/m_trilin.o \
  $(ODIR)/m_triquad.o \
  $(ODIR)/m_vars.o \
  $(ODIR)/mathcons.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_interp_ring_rad.f90 -o $(ODIR)/m_interp_ring_rad.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_interpol.o: $(SDIR)/m_interpol.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_interpol.f90 -o $(ODIR)/m_interpol.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_invert.o: $(SDIR)/m_invert.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_invert.f90 -o $(ODIR)/m_invert.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_lambda_qual.o: $(SDIR)/m_lambda_qual.f90 \
  $(ODIR)/m_cloud_pres_mod.o \
  $(ODIR)/m_find.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_lambda_qual.f90 -o $(ODIR)/m_lambda_qual.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_matmul.o: $(SDIR)/m_matmul.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_matmul.f90 -o $(ODIR)/m_matmul.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_num2string.o: $(SDIR)/m_num2string.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_num2string.f90 -o $(ODIR)/m_num2string.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_poly_fit.o: $(SDIR)/m_poly_fit.f90 \
  $(ODIR)/m_invert.o \
  $(ODIR)/m_matmul.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_poly_fit.f90 -o $(ODIR)/m_poly_fit.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_print_ret.o: $(SDIR)/m_print_ret.f90 \
  $(ODIR)/m_cloud_pres_mod.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_print_ret.f90 -o $(ODIR)/m_print_ret.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_rd_chl.o: $(SDIR)/m_rd_chl.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_rd_chl.f90 -o $(ODIR)/m_rd_chl.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_rd_terr.o: $(SDIR)/m_rd_terr.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_rd_terr.f90 -o $(ODIR)/m_rd_terr.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_rd_toms_refl.o: $(SDIR)/m_rd_toms_refl.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_rd_toms_refl.f90 -o $(ODIR)/m_rd_toms_refl.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_cal.o: $(SDIR)/m_read_cal.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_cal.f90 -o $(ODIR)/m_read_cal.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_ocean_table.o: $(SDIR)/m_read_ocean_table.f90 \
  $(ODIR)/m_interpol.o \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_ocean_table.f90 -o $(ODIR)/m_read_ocean_table.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_references.o: $(SDIR)/m_read_references.f90 \
  $(ODIR)/m_interp_pres.o \
  $(ODIR)/m_interp_ring_rad.o \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_references.f90 -o $(ODIR)/m_read_references.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_resid.o: $(SDIR)/m_read_resid.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_resid.f90 -o $(ODIR)/m_read_resid.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_tables.o: $(SDIR)/m_read_tables.f90 \
  $(ODIR)/m_find.o \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_tables.f90 -o $(ODIR)/m_read_tables.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_thresholds.o: $(SDIR)/m_read_thresholds.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_thresholds.f90 -o $(ODIR)/m_read_thresholds.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_sigma.o: $(SDIR)/m_sigma.f90 \
  $(ODIR)/m_avg.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_sigma.f90 -o $(ODIR)/m_sigma.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_sortind.o: $(SDIR)/m_sortind.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_sortind.f90 -o $(ODIR)/m_sortind.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_spline.o: $(SDIR)/m_spline.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_spline.f90 -o $(ODIR)/m_spline.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_str_replace.o: $(SDIR)/m_str_replace.f90 \
  $(ODIR)/m_strpos.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_str_replace.f90 -o $(ODIR)/m_str_replace.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_strpos.o: $(SDIR)/m_strpos.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_strpos.f90 -o $(ODIR)/m_strpos.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_swathnames.o: $(SDIR)/m_swathnames.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_swathnames.f90 -o $(ODIR)/m_swathnames.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_trilin.o: $(SDIR)/m_trilin.f90 \
  $(ODIR)/m_cloud_pres_mod.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_trilin.f90 -o $(ODIR)/m_trilin.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_triquad.o: $(SDIR)/m_triquad.f90 \
  $(ODIR)/m_cloud_pres_mod.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_triquad.f90 -o $(ODIR)/m_triquad.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_vars.o: $(SDIR)/m_vars.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_vars.f90 -o $(ODIR)/m_vars.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_write_output_data.o: $(SDIR)/m_write_output_data.f90 \
  $(ODIR)/m_vars.o \
  $(ODIR)/m_write_swath_field.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_write_output_data.f90 -o $(ODIR)/m_write_output_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_write_output_data_2pres.o: $(SDIR)/m_write_output_data_2pres.f90 \
  $(ODIR)/m_vars.o \
  $(ODIR)/m_write_swath_field.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_write_output_data_2pres.f90 -o $(ODIR)/m_write_output_data_2pres.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_write_swath_field.o: $(SDIR)/m_write_swath_field.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_write_swath_field.f90 -o $(ODIR)/m_write_swath_field.o $(MODOUTFLG) $(ODIR)
$(ODIR)/mathcons.o: $(SDIR)/mathcons.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/mathcons.f90 -o $(ODIR)/mathcons.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_HDFEOS_attr.o: $(SDIR)/m_HDFEOS_attr.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_HDFEOS_attr.f90 -o $(ODIR)/m_HDFEOS_attr.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_input_data.o: $(SDIR)/m_read_input_data.f90 \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/m_instr_config.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_read_solar_flux.o \
  $(ODIR)/m_strpos.o \
  $(ODIR)/m_swathnames.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_input_data.f90 -o $(ODIR)/m_read_input_data.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_read_solar_flux.o: $(SDIR)/m_read_solar_flux.f90 \
  $(ODIR)/hdfeos4_parameters.o \
  $(ODIR)/L1B_Reader_class.o \
  $(ODIR)/m_earth_sun_dist.o \
  $(ODIR)/m_lambda_qual.o \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_swathnames.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_read_solar_flux.f90 -o $(ODIR)/m_read_solar_flux.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_write_HDFEOS_attr.o: $(SDIR)/m_write_HDFEOS_attr.f90 \
  $(ODIR)/m_HDFEOS_attr.o \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_strpos.o \
  $(ODIR)/m_swathnames.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_write_HDFEOS_attr.f90 -o $(ODIR)/m_write_HDFEOS_attr.o $(MODOUTFLG) $(ODIR)
$(ODIR)/m_earth_sun_dist.o: $(SDIR)/m_earth_sun_dist.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_strpos.o \
  $(ODIR)/m_swathnames.o \
  $(ODIR)/m_vars.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/m_earth_sun_dist.f90 -o $(ODIR)/m_earth_sun_dist.o $(MODOUTFLG) $(ODIR)
$(ODIR)/MetadataModule.o: $(SDIR)/MetadataModule.f90 \
  $(ODIR)/m_LUN_set.o \
  $(ODIR)/m_pgs_include.o \
  $(ODIR)/m_swathnames.o \
  $(ODIR)/m_vars.o \
  $(ODIR)/MetadatOBJModule.o
	$(FC) -c $(FFLAGS) $(IFLAGS) $(SDIR)/MetadataModule.f90 -o $(ODIR)/MetadataModule.o $(MODOUTFLG) $(ODIR)
$(ODIR)/OMI_L1B_reader.o:	$(SDIR)/OMI_L1B_reader.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_L1B_reader.c -o $(ODIR)/OMI_L1B_reader.o
$(ODIR)/OMI_SMF_setmsg.o:	$(SDIR)/OMI_SMF_setmsg.c
	$(CC) -c $(CCFLAGS) $(CINCFLGS) $(SDIR)/OMI_SMF_setmsg.c -o $(ODIR)/OMI_SMF_setmsg.o
