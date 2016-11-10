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
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)
LIBCSV_LIBPATHS =# -Wl,-rpath=$(LIBCSV_LIB) -L$(LIBCSV_LIB)

IOCSDPC_INC = $(prefix)/include
IOCSDPC_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -liocsdpc

EXTRA_INCS =#
EXTRA_LIBS =#

INCS = -I$(config_dir)/include -I$(TIO_INC) -I$(IOCSDPC_INC) $(LIBCSV_INC) \
       -I$(TELL_INC) -I$(LIBCONFIG_INC) -I$(NETCDF_INC)
LIBS = $(LIBCONFIG_LIBPATHS) -lconfig \
       $(LIBCSV_LIBPATHS) -lcsv \
       $(NETCDF_LIBPATHS) -lnetcdf $(IOCSDPC_LIBPATHS) \
       $(TELL_LIBPATHS) -ltell \
       $(TIO_LIBPATHS) -ltio \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = l0_format.o tpsec.o exprec.o tpinfo.o enum.o daemon.o
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
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/l0_format.c $(LIBS)

$(objdir)/tpsec.o:	$(srcdir)/tpsec.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tpsec.c $(LIBS)

$(objdir)/exprec.o:	$(srcdir)/exprec.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/exprec.c $(LIBS)

$(objdir)/tpinfo.o:	$(srcdir)/tpinfo.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tpinfo.c $(LIBS)

$(objdir)/enum.o:	$(srcdir)/enum.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/enum.c $(LIBS)

$(objdir)/daemon.o:	$(srcdir)/daemon.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/daemon.c $(LIBS)
