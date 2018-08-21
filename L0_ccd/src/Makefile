# -*- make -*-

include ../../Makefile.defines

THIS_MODULE = L0_ccd

srcdir  = $(config_dir)/$(THIS_MODULE)/src
objdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)
GSL_LIBPATHS     = -lgsl -lgslcblas

IOCSDPC_INC = $(prefix)/include
IOCSDPC_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -liocsdpc -lcsv

EXTRA_INCS =#
EXTRA_LIBS = $(GSL_LIBPATHS) #

INCS = -I$(config_dir)/include -I$(TIO_INC) -I$(IOCSDPC_INC) \
       -I$(TELL_INC) -I$(LIBCONFIG_INC) -I$(NETCDF_INC)
LIBS = $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf $(IOCSDPC_LIBPATHS) \
       $(TELL_LIBPATHS) -ltell \
       $(TIO_LIBPATHS) -ltio -lmpfit \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l0_ccd.o process.o image.o ccd.o granule.o instr.o pixelqf.o \
           bpix.o dark.o sensorcal.o output.o util.o \
	   wavecal.o interp.o lsq.o shapefun.o \
	   interp_cheb.o interp_poly.o interp_cspline.o
OBJS     = $(patsubst %.o,$(objdir)/%.o, $(OFILES))

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	 $(objdir) $(objdir)/$(THIS_MODULE) $(objdir)/wavecal_driver $(objdir)/wavecal_merge

tests:	$(objdir) $(objdir)/wavecal_driver $(objdir)/test_interp

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
	$(INSTALL_PROGRAM) $(objdir)/wavecal_driver $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/wavecal_merge $(bindir)/

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

TEST_INTERP_OFILES = test_interp.o lsq.o interp.o \
              interp_cheb.o interp_poly.o interp_cspline.o
TEST_INTERP_OBJS = $(patsubst %.o,$(objdir)/%.o, $(TEST_INTERP_OFILES))

$(objdir)/test_interp:	$(TEST_INTERP_OBJS)
	$(COMPILE_CMD) -o $(objdir)/test_interp $(TEST_INTERP_OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/test_interp.o:	$(srcdir)/test_interp.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/test_interp.c $(LIBS)

$(objdir)/wavecal_merge:	$(objdir)/wavecal_merge.o
	$(COMPILE_CMD) -o $(objdir)/wavecal_merge $(objdir)/wavecal_merge.o \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/wavecal_merge.o:	$(srcdir)/wavecal_merge.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal_merge.c $(LIBS)

WAVECAL_DRIVER_OFILES = wavecal_driver.o util.o wavecal.o interp.o shapefun.o \
     lsq.o interp_cheb.o interp_poly.o interp_cspline.o
WAVECAL_DRIVER_OBJS = $(patsubst %.o,$(objdir)/%.o, $(WAVECAL_DRIVER_OFILES))

$(objdir)/wavecal_driver:	$(WAVECAL_DRIVER_OBJS)
	$(COMPILE_CMD) -o $(objdir)/wavecal_driver $(WAVECAL_DRIVER_OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/wavecal_driver.o:	$(srcdir)/wavecal_driver.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal_driver.c $(LIBS)

$(objdir)/$(THIS_MODULE):	$(OBJS)
	$(COMPILE_CMD) -o $(objdir)/$(THIS_MODULE) $(OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/l0_ccd.o:	$(srcdir)/l0_ccd.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l0_ccd.c $(LIBS)

$(objdir)/process.o:	$(srcdir)/process.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/process.c $(LIBS)

$(objdir)/image.o:	$(srcdir)/image.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/image.c $(LIBS)

$(objdir)/ccd.o:	$(srcdir)/ccd.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/ccd.c $(LIBS)

$(objdir)/granule.o:	$(srcdir)/granule.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/granule.c $(LIBS)

$(objdir)/instr.o:	$(srcdir)/instr.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/instr.c $(LIBS)

$(objdir)/pixelqf.o:	$(srcdir)/pixelqf.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/pixelqf.c $(LIBS)

$(objdir)/bpix.o:	$(srcdir)/bpix.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/bpix.c $(LIBS)

$(objdir)/dark.o:	$(srcdir)/dark.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/dark.c $(LIBS)

$(objdir)/sensorcal.o:	$(srcdir)/sensorcal.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/sensorcal.c $(LIBS)

$(objdir)/output.o:	$(srcdir)/output.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/output.c $(LIBS)

$(objdir)/util.o:	$(srcdir)/util.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/util.c $(LIBS)

$(objdir)/wavecal.o:	$(srcdir)/wavecal.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal.c $(LIBS)

$(objdir)/interp.o:	$(srcdir)/interp.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp.c $(LIBS)

$(objdir)/shapefun.o:	$(srcdir)/shapefun.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/shapefun.c $(LIBS)

$(objdir)/lsq.o:	$(srcdir)/lsq.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/lsq.c $(LIBS)

$(objdir)/interp_cspline.o:	$(srcdir)/interp_cspline.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp_cspline.c $(LIBS)

$(objdir)/interp_poly.o:	$(srcdir)/interp_poly.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp_poly.c $(LIBS)

$(objdir)/interp_cheb.o:	$(srcdir)/interp_cheb.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp_cheb.c $(LIBS)

