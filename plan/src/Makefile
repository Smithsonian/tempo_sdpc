# -*- make -*-

include ../../Makefile.defines

THIS_MODULE = plan

srcdir  = $(config_dir)/$(THIS_MODULE)/src
objdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
TELL_LIBPATHS = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)
PROJ_LIBPATHS = -Wl,-rpath=$(PROJ_LIB) -L$(PROJ_LIB)
SLNOVAS_LIBPATHS = -Wl,-rpath=$(SLNOVAS_LIB) -L$(SLNOVAS_LIB)

INR_LIBPATHS = -Wl,-rpath=$(INR_LIB):$(INROTS_LIB) -L$(INR_LIB) -L$(INROTS_LIB) \
               -lTempoGeospatial -llapacke $(INROTS_EXTRA_LIBS)

EXTRA_INCS = -I$(PROJ_INC) -I$(SLNOVAS_INC)
EXTRA_LIBS = $(PROJ_LIBPATHS) -lproj $(SLNOVAS_LIBPATHS) -llibnovas

INC_DIRS = $(TELL_INC) $(TIO_INC) $(NETCDF_INC) $(LIBCONFIG_INC) $(INR_INC) $(SLNOVAS_INC) \
           $(config_dir)/include
INCS = $(patsubst %,-I%, $(INC_DIRS))

LIBS = $(TIO_LIBPATHS) -ltio \
       $(LIBCONFIG_LIBPATHS) -lconfig \
       $(TELL_LIBPATHS) -ltell \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(INR_LIBPATHS) \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = plan.o bisect.o solar.o scan.o plan_list.o scan_methods.o vis.o
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
	$(INSTALL_PROGRAM) $(srcdir)/szaplot.py $(bindir)/

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

$(objdir)/plan.o:	$(srcdir)/plan.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/plan.c $(LIBS)

$(objdir)/bisect.o:	$(srcdir)/bisect.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/bisect.c $(LIBS)

$(objdir)/solar.o:	$(srcdir)/solar.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/solar.c $(LIBS)

$(objdir)/scan.o:	$(srcdir)/scan.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/scan.c $(LIBS)

$(objdir)/plan_list.o:	$(srcdir)/plan_list.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/plan_list.c $(LIBS)

$(objdir)/scan_methods.o:	$(srcdir)/scan_methods.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/scan_methods.c $(LIBS)

$(objdir)/vis.o:	$(srcdir)/vis.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/vis.c $(LIBS)

