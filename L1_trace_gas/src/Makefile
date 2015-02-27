# -*- sh -*-

EXECUTABLE = L1_trace_gas

include ../../Makefile.defines

srcdir     = $(config_dir)/L1_trace_gas/src
objdir     = $(srcdir)/$(archobjs)
msgdir     = $(srcdir)/msgs

bindir     = $(prefix)/bin
incdir     = $(prefix)/include
msgs_installdir = $(prefix)/msgs

LNFLAGS    = #
ADDLIBD    = #
ADDLIBS    = $(RPATH) -lsz -lm $(DL_LIB)

FCFLAGS  = $(MODFLG) ./ $(FFLAGS)
CCFLAGS = $(CFLAGS) -DLINUX -Df2cFortran

#---------------------------------------------------------------
# Ideally, you shouldn't need to change anything below this line
#---------------------------------------------------------------

PGSHOME     = $(OMIUTIL)/TOOLKIT

MV = /bin/mv
RM = /bin/rm
CP = /bin/cp

TIO_LIB  = $(prefix)/lib
TIO_SRC  = $(config_dir)/libtio/src
TIO_SRC_OBJDIR = $(TIO_SRC)/$(archobjs)

PSPLINELIBS  = -lpspline -lezcdf
TIOLIBS      = -Wl,-rpath=$(TIO_SRC_OBJDIR) -L$(TIO_SRC_OBJDIR) \
               -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB) -ltio -ltio_f \
               -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB) -ltell -ltell_f
NETCDFLIBS   = -lnetcdf -lnetcdff
HE4LIBS      = -lhdfeos -lmfhdf -lGctp -ldf -ljpeg -lz
HE5LIBS      = -lhe5_hdfeos -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lhdf5
PGSLIBS      = -lPGSTK
OMILIBS      =

OTS_INCDIR = $(OTS_ROOT)/include
OTS_MODDIR = $(OTS_ROOT)/mod
OTS_LIBDIR = $(OTS_ROOT)/lib

PSPLINEMOD = $(OTS_MODDIR)
PSPLINELIB = $(OTS_LIBDIR)

IFLAGS = $(MODFLG)$(HDFEOS_INC) \
         $(MODFLG)$(HDF5INC) \
         $(MODFLG)$(HDFINC) \
         $(MODFLG)$(PGSINC) \
         $(MODFLG)$(HDFEOS5_INC) \
         $(MODFLG)$(objdir) \
         $(MODFLG)$(TIO_SRC_OBJDIR) \
         $(MODFLG)$(TELL_INC) \
         $(MODFLG)$(NETCDF_INC) \
	 $(MODFLG)$(PSPLINEMOD)

LFLAGS = -L$(objdir) \
         -L$(HDFEOS_LIB) -L$(HDFLIB) -L$(PGSLIB) -L$(HDF5LIB) -L$(HDFEOS5_LIB) \
         -L$(NETCDF_LIB) -L$(PSPLINELIB) $(ADDLIBD)

LIBS   = $(PSPLINELIBS) $(TIOLIBS) $(NETCDFLIBS) \
         $(OMILIBS) $(PGSLIBS) $(HE5LIBS) $(HE4LIBS) $(ADDLIBS)

CINCFLGS = -I$(PGSINC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS5_INC) -I$(objdir)

unbufferSTDOUT_DEFINES = -DMSG_DIR='"'$(msgs_installdir)'"'
unbufferSTDOUT_INCS = -I$(TELL_INC)

# During development, make it easy to change the list of source files.
# For releases, use a static list of source files.
F90SOURCES = $(shell grep -v '\#' f90_files.lis)
CSOURCES   = $(shell grep -v '\#' c_files.lis)
LINKOBJS   = $(addprefix $(objdir)/,$(F90SOURCES:.f90=.o) $(CSOURCES:.c=.o))
# F90SOURCES = slitfunction.f90 slitfunction_omi.f90 slitfunction_asym_gaussian.f90 cache_module.f90 optimizer_interface_module.f90 elsunc_interface_module.f90 angle_sat2toa.f90 dataspline_module.f90 ezspline_interpolation.f90 fitting_loops.f90 hdfeos4_parameters.f90 he5_output_tools.f90 integration_routines.f90 metadata_tools.f90 omi_pge_fitting_aux.f90 omi_pge_fitting_process.f90 omi_pge_postprocessing.f90 swathline_loop.f90 omi_parameters_module.f90 omi_read_l1b_data.f90 OMSAO_casestring_module.f90 OMSAO_destriping_module.f90 OMSAO_elsunc_fitting_module.f90 OMSAO_errstat_module.f90 OMSAO_he5_module.f90 OMSAO_indices_module.f90 OMSAO_linterpolation_module.f90 OMSAO_main.f90 OMSAO_median_module.f90 OMSAO_omidata_module.f90 OMSAO_parameters_module.f90 OMSAO_pixelcorner_module.f90 OMSAO_precision_module.f90 OMSAO_prefitcol_module.f90 OMSAO_radiance_ref_module.f90 OMSAO_Reference_sector_module.f90 OMSAO_solar_wavcal_module.f90 OMSAO_solcomp_module.f90 OMSAO_variables_module.f90 OMSAO_wfamf_module.f90 prepare_databases.f90 radiance_fit.f90 radiance_wavcal.f90 control_module.f90 pcf_file_module.f90 read_reference_spectra.f90 sao_pge_utils.f90 strutils.f90 subtract_cubic.f90 Szoom_Parameter_Module.f90 SLATEC_davint.f90 undersample.f90 datafields.f90 irradiance_data.f90 errormodule.f90 l1bread.f90 wavecal.f90 commonmode.f90 ctrlvars.f90 arrayutils.f90
# CSOURCES   = OMI_SMF_setmsg.c unbufferSTDOUT.c l1breader.c
# LINKOBJS   = $(objdir)/arrayutils.o $(objdir)/ctrlvars.o $(objdir)/commonmode.o $(objdir)/wavecal.o $(objdir)/l1bread.o $(objdir)/errormodule.o $(objdir)/irradiance_data.o $(objdir)/datafields.o $(objdir)/undersample.o $(objdir)/SLATEC_davint.o $(objdir)/Szoom_Parameter_Module.o $(objdir)/subtract_cubic.o $(objdir)/strutils.o $(objdir)/sao_pge_utils.o $(objdir)/read_reference_spectra.o $(objdir)/pcf_file_module.o $(objdir)/control_module.o $(objdir)/radiance_wavcal.o $(objdir)/radiance_fit.o $(objdir)/prepare_databases.o $(objdir)/OMSAO_wfamf_module.o $(objdir)/OMSAO_variables_module.o $(objdir)/OMSAO_solcomp_module.o $(objdir)/OMSAO_solar_wavcal_module.o $(objdir)/OMSAO_Reference_sector_module.o $(objdir)/OMSAO_radiance_ref_module.o $(objdir)/OMSAO_prefitcol_module.o $(objdir)/OMSAO_precision_module.o $(objdir)/OMSAO_pixelcorner_module.o $(objdir)/OMSAO_parameters_module.o $(objdir)/OMSAO_omidata_module.o $(objdir)/OMSAO_median_module.o $(objdir)/OMSAO_main.o $(objdir)/OMSAO_linterpolation_module.o $(objdir)/OMSAO_indices_module.o $(objdir)/OMSAO_he5_module.o $(objdir)/OMSAO_errstat_module.o $(objdir)/OMSAO_elsunc_fitting_module.o $(objdir)/OMSAO_destriping_module.o $(objdir)/OMSAO_casestring_module.o $(objdir)/omi_read_l1b_data.o $(objdir)/omi_parameters_module.o $(objdir)/swathline_loop.o $(objdir)/omi_pge_postprocessing.o $(objdir)/omi_pge_fitting_process.o $(objdir)/omi_pge_fitting_aux.o $(objdir)/metadata_tools.o $(objdir)/integration_routines.o $(objdir)/he5_output_tools.o $(objdir)/hdfeos4_parameters.o $(objdir)/fitting_loops.o $(objdir)/ezspline_interpolation.o $(objdir)/dataspline_module.o $(objdir)/angle_sat2toa.o $(objdir)/elsunc_interface_module.o $(objdir)/optimizer_interface_module.o $(objdir)/cache_module.o $(objdir)/slitfunction_asym_gaussian.o $(objdir)/slitfunction_omi.o $(objdir)/slitfunction.o  $(objdir)/l1breader.o $(objdir)/unbufferSTDOUT.o $(objdir)/OMI_SMF_setmsg.o

#---------------------------------
# Primary 'make' targets for users
#   make [all|install|clean]
#---------------------------------

all:	dirs messages tg-names $(EXECUTABLE)

dirs:
	mkdir -p $(objdir) $(msgdir)

$(EXECUTABLE):	$(LINKOBJS)
	$(FC) main.f90 -o $(objdir)/$(EXECUTABLE) $(LINKOBJS) $(LNFLAGS) $(FCFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS)

install-executable:
	mkdir -p $(bindir)
	@echo installing $(EXECUTABLE) in $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/$(EXECUTABLE) $(bindir)/$(EXECUTABLE)

install-header:
	mkdir -p $(incdir)
	@echo installing $(NAME_INCLUDE) in $(incdir)/
	$(INSTALL_DATA) $(srcdir)/$(NAME_INCLUDE) $(incdir)/$(NAME_INCLUDE)

install-msgs:
	mkdir -p $(msgs_installdir)
	@echo installing message files in $(msgs_installdir)/
	for m in $(MESSAGE_FILES) ; do \
	   $(INSTALL_DATA) $(srcdir)/msgs/$$m $(msgs_installdir)/$$m ; \
	done

install:	install-executable install-header install-msgs

clean:	clean-tg-names
	-$(RM) $(objdir)/*.o $(objdir)/*.mod

distclean:	clean
	-$(RM) -rf $(objdir)

tags:
	ctags $(F90SOURCES) $(CSOURCES)
#-----------------------------------------------------------------
# The hdf.f90 file contains numerous non-conforming tabs.
# It's unlikely anyone will ever fix that in the HDF distribution,
# so we'll just fix it this way:
#-----------------------------------------------------------------
fix-hdf-tabs:
	sed --in-place=_ 's,\t, ,g' $(OTS_ROOT)/include/hdf.f90

# ---------------------------------------
# The PGE and TOOLKIT error message files
# ---------------------------------------

MESSAGE_FILES = PGS_9 PGS_1900 PGS_52500

PGEPGSSRC = OMSAO_52500.t OMI_1900.t
PGEPGSCMP = $(msgdir)/PGS_52500 $(msgdir)/PGS_1900
TKPGSCMP = PGS_9

messages: $(PGEPGSCMP) tkmessages

tkmessages:
	@for tkpgs in $(TKPGSCMP); do                                 \
	  if [ ! -e $(msgdir)/$$tkpgs ]; then                        \
	     echo 'Copying' $$tkpgs 'to' $(msgdir);                  \
	     touch $(msgdir)/$$tkpgs ; chmod u+w $(msgdir)/$$tkpgs;  \
	     $(CP) -p $(PGSHOME)/message/$$tkpgs $(msgdir)/;            \
	  fi;                                                        \
	done

# ----------------------------------------------------------------------
# Pattern rule for the creation of SMF messages and associated .f and .h
# ----------------------------------------------------------------------
$(msgdir)/PGS_1900:
	export PGSMSG="$(msgdir)" ; \
	$(PGSBIN)/smfcompile -r -f77 -f OMI_1900.t ; \
	$(PGSBIN)/smfcompile -r      -f OMI_1900.t
$(msgdir)/PGS_52500:
	export PGSMSG="$(msgdir)" ; \
	$(PGSBIN)/smfcompile -r -f77 -f OMSAO_52500.t ; \
	$(PGSBIN)/smfcompile -r      -f OMSAO_52500.t

#
# Name tables
#
NAME_FILE_BASENAMES = tg_names_dim tg_names_var
NAME_INCLUDE = tg_names.h

NAME_FILE_SOURCES = $(addsuffix .csv, $(NAME_FILE_BASENAMES))
NAME_FILE_TARGETS = $(addsuffix .inc, $(NAME_FILE_BASENAMES)) $(NAME_INCLUDE)
tg-names:	$(NAME_FILE_TARGETS)
$(NAME_FILE_TARGETS):	$(NAME_FILE_SOURCES)
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("#define %35s   %s\n",a[1],a[2])}' \
	        $${bn}.csv > .$${bn}_h.tmp ; \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("character (len=tg_max_name_len), public, parameter :: %35s = %s\n",tolower(a[1]),a[2])}' \
	        $${bn}.csv > $${bn}.inc ; \
	   sed -i -e "1i! Generated from $${bn}.csv by Makefile. Do not edit." \
	        $${bn}.inc ; \
	done
	sed -e '/@CSV_DIMS@/{r .tg_names_dim_h.tmp' -e 'd}' \
	    -e '/@CSV_VARS@/{r .tg_names_var_h.tmp' -e 'd}' \
	    $(srcdir)/$(NAME_INCLUDE)in > $(NAME_INCLUDE)
	sed -i -e '1i/* Generated by Makefile. Do not edit. */' \
	     $(NAME_INCLUDE)

clean-tg-names:
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f .$${bn}_h.tmp ; \
	done
distclean-tg-names:	clean-tg-names
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f $${bn}.inc $(NAME_INCLUDE) ; \
	done

# ---- build dependencies ----
# During development, the rules will be dynamically updated.
# For releases, a static set of rules can be appended.
Makefile:	rules.inc
rules.inc:	$(F90SOURCES) $(CSOURCES) f90_files.lis c_files.lis
	f90depmod.sl -r $(F90SOURCES) > rules.inc
	./mkrules.sh c $(CSOURCES) >> rules.inc

ODIR = $(objdir)
SDIR = $(srcdir)
-include rules.inc

.PHONY: tags
