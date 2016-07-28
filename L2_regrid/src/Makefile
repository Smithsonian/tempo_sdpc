# -*- make -*-

EXECUTABLE = L2_regrid

include ../../Makefile.defines

THIS_LIB = regrid

srcdir  = $(config_dir)/L2_regrid/src
elfdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
REGRID_LIBPATHS  = -Wl,-rpath=$(elfdir) -L$(elfdir)
PROJ_LIBPATHS    = -Wl,-rpath=$(PROJ_LIB) -L$(PROJ_LIB)
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)

EXTRA_INCS = -I$(PROJ_INC)
EXTRA_LIBS = $(PROJ_LIBPATHS) -lproj

INCS = -I$(config_dir)/include -I$(TIO_INC) -I$(TELL_INC) -I$(LIBCONFIG_INC) -I$(NETCDF_INC)
LIBS = $(INSTALL_LIBPATHS) $(LIBCONFIG_LIBPATHS) -lconfig \
       $(TIO_LIBPATHS) -ltio \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell $(RPATH)

ELF_MAJOR_VERSION = 0#
ELF_MINOR_VERSION = 0#
ELF_MICRO_VERSION = 1#

# On most systems, these variables have values of the form:
#   libname.so, libname.so.X, libname.so.X.Y, libname.so.X.Y.Z
ELFLIB = lib$(THIS_LIB).so
ELFLIB_MAJOR = $(ELFLIB).$(ELF_MAJOR_VERSION)
ELFLIB_MAJOR_MINOR = $(ELFLIB_MAJOR).$(ELF_MINOR_VERSION)
ELFLIB_MAJOR_MINOR_MICRO = $(ELFLIB_MAJOR_MINOR).$(ELF_MICRO_VERSION)

# Set this to one of the above ELFLIB_MAJOR_* values
ELFLIB_BUILD_NAME = $(ELFLIB_MAJOR_MINOR_MICRO)
ELFDIR_ELF_LIB = $(elfdir)/$(ELFLIB_BUILD_NAME)#

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

OFILES   = poly.o pixel.o regrid.o var.o
ELFOBJS  = $(elfdir)/poly.o $(elfdir)/pixel.o $(elfdir)/regrid.o $(elfdir)/var.o

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	elf $(elfdir)/$(EXECUTABLE)
elf:	$(elfdir) $(ELFDIR_ELF_LIB)

tags:
	ctags $(patsubst %.o,%.c,$(OFILES) l2_regrid.o)

$(elfdir):
	-$(MKINSDIR) $(elfdir)

$(ELFDIR_ELF_LIB):	$(elfdir) $(ELFOBJS)
	-$(RM) -f $(ELFDIR_ELF_LIB) $(elfdir)/$(ELFLIB_MAJOR) $(elfdir)/$(ELFLIB)
	cd $(elfdir) && $(ELF_LINK_CMD) -o $(ELFLIB_BUILD_NAME) $(OFILES) $(LDFLAGS) $(LIBS) $(EXTRA_LIBS) $(ELF_DEP_LIBS)
	if [ $(ELFLIB_MAJOR) != $(ELFLIB_BUILD_NAME) ]; then \
	  cd $(elfdir) && $(LN) $(ELFLIB_BUILD_NAME) $(ELFLIB_MAJOR); \
	fi
	cd $(elfdir) && $(LN) $(ELFLIB_BUILD_NAME) $(ELFLIB)

.PHONY: clean distclean all install check tags

# --------------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------------

bindir = $(prefix)/bin
incdir = $(prefix)/include
libdir = $(prefix)/lib

install-executable:
	mkdir -p $(bindir)
	@echo installing $(EXECUTABLE) in $(bindir)/
	$(INSTALL_PROGRAM) $(elfdir)/$(EXECUTABLE) $(bindir)/$(EXECUTABLE)

install-elf:
	$(MKINSDIR) $(libdir)
	-$(RM) $(libdir)/$(ELFLIB)
	-$(RM) $(libdir)/$(ELFLIB_MAJOR)
	@echo installing $(ELFLIB_BUILD_NAME) in $(libdir)
	$(INSTALL_DATA) $(ELFDIR_ELF_LIB) $(libdir)
	-cd $(libdir) && $(LN) $(ELFLIB_BUILD_NAME) $(ELFLIB_MAJOR)
	-cd $(libdir) && $(LN) $(ELFLIB_MAJOR) $(ELFLIB)
	@echo ""
	@echo $(ELFLIB_BUILD_NAME) created in $(libdir)/.
	@echo The links $(ELFLIB) and $(ELFLIB_MAJOR) to $(ELFLIB_BUILD_NAME) were also created.
	@echo ""

install:	elf install-elf install-executable

# --------------------------------------------------------------------------
# Housekeeping
# --------------------------------------------------------------------------

clean:
	-$(RM) -f .##* $(elfdir)/*
distclean: clean
	-$(RM) -rf $(elfdir)

# --------------------------------------------------------------------------
# Object Rules
# --------------------------------------------------------------------------

$(elfdir)/$(EXECUTABLE):	$(ELFDIR_ELF_LIB) $(elfdir)/l2_regrid.o
	$(COMPILE_CMD) -o $(elfdir)/$(EXECUTABLE) $(elfdir)/l2_regrid.o \
	    $(INCS) $(EXTRA_INCS) $(REGRID_LIBPATHS) -lregrid $(LIBS) $(EXTRA_LIBS) -lm

$(elfdir)/poly.o:	$(srcdir)/poly.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/poly.c $(LIBS)

$(elfdir)/pixel.o:	$(srcdir)/pixel.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/pixel.c $(LIBS)

$(elfdir)/regrid.o:	$(srcdir)/regrid.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(EXTRA_INCS) $(srcdir)/regrid.c

$(elfdir)/var.o:	$(srcdir)/var.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(EXTRA_INCS) $(srcdir)/var.c

$(elfdir)/l2_regrid.o:	$(srcdir)/l2_regrid.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(EXTRA_INCS) $(srcdir)/l2_regrid.c
