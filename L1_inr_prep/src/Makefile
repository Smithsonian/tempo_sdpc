# -*- make -*-

include ../../Makefile.defines

THIS_MODULE = L1_inr_prep

srcdir  = $(config_dir)/$(THIS_MODULE)/src
objdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

IOCLIB_INC = $(prefix)/include
IOCLIB_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -lcsv

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)

EXTRA_INCS =
EXTRA_LIBS =

INC_DIRS = $(TIO_INC) $(IOCLIB_INC) $(TELL_INC) $(LIBCONFIG_INC) $(NETCDF_INC) \
           $(config_dir)/include
INCS = $(patsubst %,-I%, $(INC_DIRS))

LIBS = $(TIO_LIBPATHS) -ltio \
       $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell \
       $(IOCLIB_LIBPATHS) \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l1_inr_prep.o bsearch.o radiance.o row_select.o ephem.o
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

$(objdir)/l1_inr_prep.o:	$(srcdir)/l1_inr_prep.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l1_inr_prep.c $(LIBS)

$(objdir)/row_select.o:	$(srcdir)/row_select.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/row_select.c $(LIBS)

$(objdir)/radiance.o:	$(srcdir)/radiance.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/radiance.c $(LIBS)

$(objdir)/ephem.o:	$(srcdir)/ephem.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/ephem.c $(LIBS)

$(objdir)/bsearch.o:	$(srcdir)/bsearch.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/bsearch.c $(LIBS)

