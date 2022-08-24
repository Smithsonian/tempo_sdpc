# -*- sh -*-

include ../Makefile.defines

COMPILE_CMD = $(CC) $(CCFLAGS)

all:
	@echo "targets: install"

SCRIPT_FILES = plot_iru.py plot_smc.py plot_image.py plot_hk.py \
               plot_level2.py plot_level3.py sdpc_panoply.sh sdpc_logview.sh \
               check_time_series.py plot_trend.py asdc_upload_report.py \
               production_report.py update_leapsec_file.sh

SCRIPT_FILES_INSTALL_DIR = $(prefix)/bin

install_directories:
	$(MKINSDIR) $(SCRIPT_FILES_INSTALL_DIR)

install_scripts:
	@for X in $(SCRIPT_FILES); \
	do \
		echo $(INSTALL) $$X $(SCRIPT_FILES_INSTALL_DIR); \
		$(INSTALL) $$X $(SCRIPT_FILES_INSTALL_DIR); \
		if [ "$$?" != "0" ]; then \
		  exit 1; \
		fi; \
	done

install:	install_directories install_scripts

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

SLNOVAS_LIBPATHS = -Wl,-rpath=$(SLNOVAS_LIB) -L$(SLNOVAS_LIB)
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)

CONVERT_DEFINES = # -DDEBUG
CONVERT_INCS = -I$(SLNOVAS_INC) -I$(INR_INC) -I$(TIO_INC) -I$(NETCDF_INC)
CONVERT_LIBS = -Wl,-rpath=$(INR_LIB) -L$(INR_LIB) -lj2e -lextern \
               $(SLNOVAS_LIBPATHS) -llibnovas \
	       $(TIO_LIBPATHS) -ltio \
	       $(NETCDF_LIBPATHS) -lnetcdf

ecef_to_j2k:	convert_ecef_to_j2k.c
	$(COMPILE_CMD) -o $@ convert_ecef_to_j2k.c $(CONVERT_DEFINES) $(CONVERT_INCS) $(CONVERT_LIBS) -lm

check:

docs:

clean:
	$(RM) -f ecef_to_j2k

distclean:



