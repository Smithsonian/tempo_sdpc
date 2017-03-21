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

INC_DIRS = $(config_dir)/include $(REGRID_INC) \
           $(TIO_INC) $(TELL_INC) $(LIBCONFIG_INC) $(NETCDF_INC) \
           $(GEOTIFF_INC) $(TIFF_INC) $(INR_INC)
INCS = $(patsubst %,-I%, $(INC_DIRS))

LIBS = $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell \
       $(TIO_LIBPATHS) -ltio \
       $(REGRID_LIBPATHS) -lregrid \
       $(GEOTIFF_LIBPATHS) $(TIFF_LIBPATHS) -lgeotiff -ltiff \
       $(INR_LIBPATHS) \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l1_inr_post.o map.o elevation.o land_cover.o granule.o snow.o ezlhconv.o
OBJS     = $(patsubst %.o,$(objdir)/%.o, $(OFILES))

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	 $(objdir) $(objdir)/$(THIS_MODULE)

$(objdir):
	-mkdir -p $(objdir)

tags:
	ctags $(patsubst %.o,%.c,$(OFILES))

.PHONY: clean distclean all install check tags

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
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l1_inr_post.c $(LIBS)

$(objdir)/map.o:	$(srcdir)/map.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/map.c $(LIBS)

$(objdir)/elevation.o:	$(srcdir)/elevation.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/elevation.c $(LIBS)

$(objdir)/land_cover.o:	$(srcdir)/land_cover.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/land_cover.c $(LIBS)

$(objdir)/granule.o:	$(srcdir)/granule.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/granule.c $(LIBS)

$(objdir)/snow.o:	$(srcdir)/snow.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/snow.c $(LIBS)

$(objdir)/ezlhconv.o:	$(srcdir)/ezlhconv.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/ezlhconv.c $(LIBS)
