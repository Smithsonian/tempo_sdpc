# -*- make -*-

include ../../Makefile.defines

THIS_MODULE = L2_split

srcdir  = $(config_dir)/$(THIS_MODULE)/src
objdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

REGRID_INC = $(config_dir)/L2_regrid/src
REGRID_LIB = $(config_dir)/L2_regrid/src/$(archobjs)
REGRID_LIBPATHS  = -Wl,-rpath=$(REGRID_LIB) -L$(REGRID_LIB)

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
PROJ_LIBPATHS    = -Wl,-rpath=$(PROJ_LIB) -L$(PROJ_LIB)
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)
GSL_LIBPATHS     = -lgsl -lgslcblas

EXTRA_INCS = -I$(PROJ_INC)
EXTRA_LIBS = $(PROJ_LIBPATHS) -lproj $(GSL_LIBPATHS)

INCS = -I$(TIO_INC) -I$(REGRID_INC) -I$(TELL_INC) -I$(LIBCONFIG_INC) -I$(NETCDF_INC) \
       -I$(config_dir)/include
LIBS = $(TIO_LIBPATHS) -ltio \
       $(REGRID_LIBPATHS) -lregrid \
       $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l2_split.o proj.o scan.o process.o filter.o
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

$(objdir)/l2_split.o:	$(srcdir)/l2_split.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l2_split.c

$(objdir)/proj.o:	$(srcdir)/proj.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/proj.c

$(objdir)/scan.o:	$(srcdir)/scan.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/scan.c

$(objdir)/process.o:	$(srcdir)/process.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/process.c

$(objdir)/filter.o:	$(srcdir)/filter.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/filter.c
