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

SLNOVAS_LIBPATHS = -Wl,-rpath=$(SLNOVAS_LIB) -L$(SLNOVAS_LIB) -llibnovas

IOCSDPC_INC = $(prefix)/include
IOCSDPC_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -liocsdpc -lcsv

EXTRA_INCS =#
EXTRA_LIBS = $(GSL_LIBPATHS) $(SLNOVAS_LIBPATHS) #

INCS = -I$(TIO_INC) -I$(IOCSDPC_INC) -I$(TELL_INC) \
        -I$(LIBCONFIG_INC) -I$(NETCDF_INC) -I$(SLNOVAS_INC) \
        -I$(config_dir)/include
LIBS = $(TIO_LIBPATHS) -ltio -lmpfit \
       $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf $(IOCSDPC_LIBPATHS) \
       $(TELL_LIBPATHS) -ltell \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

# Yes, this is an ugly hack.
BORROWED_OFILES = solar.o
OFILES   = l0_ccd.o process.o image.o ccd.o granule.o instr.o pixelqf.o \
           bpix.o dark.o sensorcal.o output.o util.o
OBJS     = $(patsubst %.o,$(objdir)/%.o, $(OFILES)) \
           $(patsubst %.o,$(objdir)/%.o, $(BORROWED_OFILES))

LINEARITY_OFILES = l0_linearity.o linearity.o image.o ccd.o granule.o util.o
LINEARITY_OBJS = $(patsubst %.o,$(objdir)/%.o, $(LINEARITY_OFILES))

EXECUTABLES = $(objdir)/$(THIS_MODULE) $(objdir)/L0_linearity

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	 $(objdir) $(EXECUTABLES)

$(objdir):
	-mkdir -p $(objdir)

tags:
	ctags $(patsubst %.o,%.c,$(OFILES) $(LINEARITY_OFILES))

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
	for pgm in $(EXECUTABLES) ; do \
	    $(INSTALL_PROGRAM) $$pgm $(bindir)/ ; \
	done

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

$(objdir)/L0_linearity:	$(LINEARITY_OBJS)
	$(COMPILE_CMD) -o $(objdir)/L0_linearity $(LINEARITY_OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

plan_srcdir = $(config_dir)/plan/src

$(objdir)/l0_ccd.o:	$(srcdir)/l0_ccd.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l0_ccd.c

$(objdir)/process.o:	$(srcdir)/process.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) -I$(plan_srcdir) $(srcdir)/process.c

$(objdir)/image.o:	$(srcdir)/image.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/image.c

$(objdir)/ccd.o:	$(srcdir)/ccd.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/ccd.c

$(objdir)/granule.o:	$(srcdir)/granule.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/granule.c

$(objdir)/instr.o:	$(srcdir)/instr.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/instr.c

$(objdir)/pixelqf.o:	$(srcdir)/pixelqf.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/pixelqf.c

$(objdir)/bpix.o:	$(srcdir)/bpix.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/bpix.c

$(objdir)/dark.o:	$(srcdir)/dark.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/dark.c

$(objdir)/sensorcal.o:	$(srcdir)/sensorcal.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/sensorcal.c

$(objdir)/output.o:	$(srcdir)/output.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/output.c

$(objdir)/util.o:	$(srcdir)/util.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/util.c

$(objdir)/solar.o:	$(plan_srcdir)/solar.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) -I$(plan_srcdir) $(plan_srcdir)/solar.c

$(objdir)/l0_linearity.o:	$(srcdir)/l0_linearity.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l0_linearity.c

$(objdir)/linearity.o:	$(srcdir)/linearity.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/linearity.c
