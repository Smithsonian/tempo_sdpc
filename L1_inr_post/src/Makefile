# -*- make -*-

include ../../Makefile.defines

THIS_MODULE = L1_inr_post

srcdir  = $(config_dir)/$(THIS_MODULE)/src
objdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

REGRID_INC = $(config_dir)/L2_regrid/src
REGRID_LIB = $(config_dir)/L2_regrid/src/$(archobjs)

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
REGRID_LIBPATHS  = -Wl,-rpath=$(REGRID_LIB) -L$(REGRID_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)
TIFF_LIBPATHS = -Wl,-rpath=$(TIFF_LIB) -L$(TIFF_LIB)
GEOTIFF_LIBPATHS = -Wl,-rpath=$(GEOTIFF_LIB) -L$(GEOTIFF_LIB)
PROJ_LIBPATHS = -Wl,-rpath=$(PROJ_LIB) -L$(PROJ_LIB)

INR_LIBPATHS = -Wl,-rpath=$(INR_LIB):$(INROTS_LIB) -L$(INR_LIB) -L$(INROTS_LIB) \
               -lTempoGeospatial -llapacke $(INROTS_EXTRA_LIBS)

EXTRA_INCS = -I$(PROJ_INC)
EXTRA_LIBS = $(PROJ_LIBPATHS) -lproj

INC_DIRS = $(REGRID_INC) $(TIO_INC) $(TELL_INC) $(LIBCONFIG_INC) $(NETCDF_INC) \
           $(GEOTIFF_INC) $(TIFF_INC) $(INR_INC) $(config_dir)/include
INCS = $(patsubst %,-I%, $(INC_DIRS))

LIBS = $(TIO_LIBPATHS) -ltio \
       $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell \
       $(REGRID_LIBPATHS) -lregrid \
       $(GEOTIFF_LIBPATHS) $(TIFF_LIBPATHS) -lgeotiff -ltiff \
       $(INR_LIBPATHS) \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l1_inr_post.o map.o elevation.o land_cover.o granule.o snow.o ezlhconv.o
OBJS     = $(patsubst %.o,$(objdir)/%.o, $(OFILES))

ELF_FCFLAGS = $(FCFLAGS) -fPIC
ELF_F_LINK_CMD = $(ELF_F_LINK)
ELF_F_COMPILE_CMD = $(ELF_FC) -c $(ELF_FCFLAGS)

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	 $(objdir) $(objdir)/$(THIS_MODULE) $(objdir)/L1_polcorr

polcorr:	$(objdir)/L1_polcorr

tests:	 $(objdir)/ndinterp_test

$(objdir):
	-mkdir -p $(objdir)

tags:
	ctags $(patsubst %.o,%.c,$(OFILES)) \
	      $(patsubst %.o,%.c,$(POLCORR_C_OFILES)) \
	      $(patsubst %.o,%.f90,$(POLCORR_F_OFILES))

.PHONY: clean distclean all install check tags polcorr

# --------------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------------

bindir = $(prefix)/bin
incdir = $(prefix)/include
libdir = $(prefix)/lib

install-executables:
	mkdir -p $(bindir)
	@echo installing executables in $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/$(THIS_MODULE) $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/L1_polcorr $(bindir)/

install:	install-executables

# --------------------------------------------------------------------------
# Housekeeping
# --------------------------------------------------------------------------

clean:
	-$(RM) -f .##* $(objdir)/*
distclean: clean
	-$(RM) -rf $(objdir)

# --------------------------------------------------------------------------
# Object Rules
# --------------------------------------------------------------------------

$(objdir)/$(THIS_MODULE):	$(OBJS)
	$(COMPILE_CMD) -o $(objdir)/$(THIS_MODULE) $(OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/l1_inr_post.o:	$(srcdir)/l1_inr_post.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l1_inr_post.c

$(objdir)/map.o:	$(srcdir)/map.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/map.c

$(objdir)/elevation.o:	$(srcdir)/elevation.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/elevation.c

$(objdir)/land_cover.o:	$(srcdir)/land_cover.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/land_cover.c

$(objdir)/granule.o:	$(srcdir)/granule.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/granule.c

$(objdir)/snow.o:	$(srcdir)/snow.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/snow.c

$(objdir)/ezlhconv.o:	$(srcdir)/ezlhconv.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/ezlhconv.c

# -----------------------------------------------------------------------

GSL_LIBPATHS = -lgsl -lgslcblas

POLCORR_F_OFILES = ndinterp_module.o types_module.o polpredict_module.o polcorrect_module.o
POLCORR_C_OFILES = l1_polcorr.o lps.o spline.o lps_apply.o
POLCORR_OFILES = $(POLCORR_C_OFILES) $(POLCORR_F_OFILES)

POLCORR_OBJS = $(patsubst %.o,$(objdir)/%.o, $(POLCORR_OFILES))

POLCORR_INC_DIRS = $(config_dir)/include $(PROJ_INC) \
                   $(TIO_INC) $(TELL_INC) $(LIBCONFIG_INC) $(NETCDF_INC)
POLCORR_INCS = $(patsubst %,-I%, $(POLCORR_INC_DIRS))

POLCORR_LIBS = $(LIBCONFIG_LIBPATHS) -lconfig \
               $(NETCDF_LIBPATHS) -lnetcdff -lnetcdf \
               $(TELL_LIBPATHS) -ltell_f -ltell \
               $(TIO_LIBPATHS) -ltio_f -ltio \
	       $(PROJ_LIBPATHS) -lproj \
	       $(GSL_LIBPATHS) \
               $(RPATH) $(INSTALL_LIBPATHS) -lgfortran

MODS_F = -I$(NETCDF_INC) -I$(TELL_INC) -I$(TIO_LIB) -I$(prefix)/include

$(objdir)/L1_polcorr:	$(POLCORR_OBJS)
	$(COMPILE_CMD) -o $(objdir)/L1_polcorr $(POLCORR_OBJS) $(POLCORR_LIBS) -lm

$(objdir)/l1_polcorr.o:	$(srcdir)/l1_polcorr.c
	cd $(objdir) && $(ELF_COMPILE_CMD) -o $(objdir)/l1_polcorr.o $(POLCORR_INCS) $(srcdir)/l1_polcorr.c

$(objdir)/lps.o:	$(srcdir)/lps.c
	cd $(objdir) && $(ELF_COMPILE_CMD) -o $(objdir)/lps.o $(POLCORR_INCS) $(srcdir)/lps.c

$(objdir)/lps_apply.o:	$(srcdir)/lps_apply.c
	cd $(objdir) && $(ELF_COMPILE_CMD) -o $(objdir)/lps_apply.o $(POLCORR_INCS) $(srcdir)/lps_apply.c

$(objdir)/spline.o:	$(srcdir)/spline.c
	cd $(objdir) && $(ELF_COMPILE_CMD) -o $(objdir)/spline.o $(POLCORR_INCS) $(srcdir)/spline.c

$(objdir)/types_module.o:	$(srcdir)/types_module.f90
	cd $(objdir) && $(ELF_F_COMPILE_CMD) $(MODS_F) $(srcdir)/types_module.f90

$(objdir)/polcorrect_module.o:	$(srcdir)/polcorrect_module.f90
	cd $(objdir) && $(ELF_F_COMPILE_CMD) $(MODS_F) $(srcdir)/polcorrect_module.f90

$(objdir)/polpredict_module.o:	$(srcdir)/polpredict_module.f90
	cd $(objdir) && $(ELF_F_COMPILE_CMD) $(MODS_F) $(srcdir)/polpredict_module.f90

NDINTERP_AUTOGEN_SRC = ndinterp_module_code.in ndinterp_module_autogen.sh

ndinterp-autogen-timestamp:	ndinterp_module.f90 $(NDINTERP_AUTOGEN_SRC)
	./ndinterp_module_autogen.sh
	touch $@

$(objdir)/ndinterp_module.o:	ndinterp_module.f90 ndinterp-autogen-timestamp
	$(ELF_F_COMPILE_CMD) -o $(objdir)/ndinterp_module.o -I$(TELL_INC) ndinterp_module.f90

NDI_LIBS = $(TELL_LIBPATHS) -ltell_f -ltell
$(objdir)/ndinterp_test:	$(objdir)/ndinterp_module.o ndinterp_test.f90
	$(ELF_FC) $(ELF_FCFLAGS) -o $(objdir)/ndinterp_test ndinterp_test.f90 $(objdir)/ndinterp_module.o $(NDI_LIBS)
