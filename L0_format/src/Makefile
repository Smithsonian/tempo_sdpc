# -*- make -*-

include ../../Makefile.defines

THIS_MODULE = L0_format

srcdir  = $(config_dir)/$(THIS_MODULE)/src
objdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

LIBCSV_INC =#
LIBCSV_LIB =#

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -L$(TIO_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)
LIBCSV_LIBPATHS =# -Wl,-rpath=$(LIBCSV_LIB) -L$(LIBCSV_LIB)

IOCSDPC_INC = $(prefix)/include
IOCSDPC_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -liocsdpc

EXTRA_INCS =#
EXTRA_LIBS =#

INCS = -I$(TIO_INC) -I$(IOCSDPC_INC) $(LIBCSV_INC) \
       -I$(TELL_INC) -I$(LIBCONFIG_INC) -I$(NETCDF_INC) -I$(config_dir)/include
LIBS = $(TIO_LIBPATHS) -ltio \
       $(LIBCONFIG_LIBPATHS) -lconfig \
       $(LIBCSV_LIBPATHS) -lcsv \
       $(NETCDF_LIBPATHS) -lnetcdf $(IOCSDPC_LIBPATHS) \
       $(TELL_LIBPATHS) -ltell \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l0_format.o tpsec.o exprec_cache_disk.o exprec_cache_mem.o \
           exprec.o tpinfo.o smc.o iru.o enum.o maneuver.o ephemeris.o
OBJS     = $(patsubst %.o,$(objdir)/%.o, $(OFILES))

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	 $(objdir) $(objdir)/L0_format

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
	$(INSTALL_PROGRAM) $(objdir)/L0_format $(bindir)/

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

$(objdir)/L0_format:	$(OBJS)
	$(COMPILE_CMD) -o $(objdir)/L0_format $(OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/l0_format.o:	$(srcdir)/l0_format.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l0_format.c

$(objdir)/tpsec.o:	$(srcdir)/tpsec.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tpsec.c

$(objdir)/exprec.o:	$(srcdir)/exprec.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/exprec.c

$(objdir)/exprec_cache_disk.o:	$(srcdir)/exprec_cache_disk.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/exprec_cache_disk.c

$(objdir)/exprec_cache_mem.o:	$(srcdir)/exprec_cache_mem.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/exprec_cache_mem.c

$(objdir)/smc.o:	$(srcdir)/smc.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/smc.c

$(objdir)/iru.o:	$(srcdir)/iru.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/iru.c

$(objdir)/tpinfo.o:	$(srcdir)/tpinfo.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tpinfo.c

$(objdir)/enum.o:	$(srcdir)/enum.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/enum.c

$(objdir)/maneuver.o:	$(srcdir)/maneuver.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/maneuver.c

$(objdir)/ephemeris.o:	$(srcdir)/ephemeris.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/ephemeris.c
