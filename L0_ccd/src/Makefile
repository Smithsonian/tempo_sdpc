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

EXTRA_INCS =#
EXTRA_LIBS = $(GSL_LIBPATHS) #

INCS = -I$(config_dir)/include -I$(TIO_INC) \
       -I$(TELL_INC) -I$(LIBCONFIG_INC) -I$(NETCDF_INC)
LIBS = $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell \
       $(TIO_LIBPATHS) -ltio \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l0_ccd.o process.o image.o ccd.o granule.o instr.o pixelqf.o \
           bpix.o dark.o sensorcal.o output.o
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

