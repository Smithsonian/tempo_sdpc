# -*- make -*-

include ../../Makefile.defines

THIS_LIB = gler

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/$(archobjs)

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
TIO_LIBPATHS     = -Wl,-rpath=$(config_dir)/libtio/src/$(archobjs) -L$(config_dir)/libtio/src/$(archobjs)
TIO_INC          = $(config_dir)/libtio/src
NDINTERP_LIBPATHS = -Wl,-rpath=$(config_dir)/libndinterp/src/$(archobjs) -L$(config_dir)/libndinterp/src/$(archobjs)
NDINTERP_INC     = $(config_dir)/libndinterp/src

IOCSDPC_INC = $(prefix)/include
IOCSDPC_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -liocsdpc -lcsv

GLER_LIBPATHS    = -Wl,-rpath=$(elfdir) -L$(elfdir)

INCS = -I$(TELL_INC) -I$(NETCDF_INC) -I$(TIO_INC)
LIBS = -Wl,-rpath=$(elfdir) -L$(elfdir) \
        $(NETCDF_LIBPATHS) -lnetcdf \
	$(TELL_LIBPATHS) -ltell \
	$(TIO_LIBPATHS) -ltio $(RPATH) \
	$(IOCSDPC_LIBPATHS) $(INSTALL_LIBPATHS)

MODS_F = -I$(NETCDF_INC) -I$(TELL_INC) -I$(TIO_INC)/$(archobjs) -I$(NDINTERP_INC)
LIBS_F = $(NETCDF_LIBPATHS) -lnetcdff \
         $(TELL_LIBPATHS) -ltell_f -ltell \
	 $(TIO_LIBPATHS) -ltio_f -ltio \
	 $(NDINTERP_LIBPATHS) -lndinterp_f

ELF_MAJOR_VERSION = 0#
ELF_MINOR_VERSION = 0#
ELF_MICRO_VERSION = 1#

# On most systems, these variables have values of the form:
#   libtio.so, libtio.so.X, libtio.so.X.Y, libtio.so.X.Y.Z
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

OFILES   = gler_util.o
ELFOBJS  = $(patsubst %.o,$(elfdir)/%.o, $(OFILES))

# On most systems, these variables have values of the form:
#   libtio.so, libtio.so.X, libtio.so.X.Y, libtio.so.X.Y.Z
# --- Fortran interface
ELFLIB_F = lib$(THIS_LIB)_f.so
ELFLIB_F_MAJOR = $(ELFLIB_F).$(ELF_MAJOR_VERSION)
ELFLIB_F_MAJOR_MINOR = $(ELFLIB_F_MAJOR).$(ELF_MINOR_VERSION)
ELFLIB_F_MAJOR_MINOR_MICRO = $(ELFLIB_F_MAJOR_MINOR).$(ELF_MICRO_VERSION)
ELFLIB_F_BUILD_NAME = $(ELFLIB_F_MAJOR_MINOR_MICRO)
ELFDIR_F_ELF_LIB = $(elfdir)/$(ELFLIB_F_BUILD_NAME)#

ELF_FCFLAGS = $(FCFLAGS) -fPIC
ELF_F_LINK_CMD = $(ELF_F_LINK)
ELF_F_COMPILE_CMD = $(ELF_FC) -c $(ELF_FCFLAGS)

OFILES_F = gler_module.o
ELFOBJS_F = $(patsubst %.o,$(elfdir)/%.o, $(OFILES_F))

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

EXECUTABLES = gler_test

all:	elf
elf:	$(elfdir) $(ELFDIR_F_ELF_LIB) # $(ELFDIR_ELF_LIB)

tags:
	ctags $(patsubst %.o,%.c,$(OFILES)) $(patsubst %.o,%.f90,$(OFILES_F))

$(elfdir):
	-$(MKINSDIR) $(elfdir)

# $(ELFDIR_ELF_LIB):	$(elfdir) $(ELFOBJS)
# 	-$(RM) -f $(ELFDIR_ELF_LIB) $(elfdir)/$(ELFLIB_MAJOR) $(elfdir)/$(ELFLIB)
# 	cd $(elfdir) && $(ELF_LINK_CMD) -o $(ELFLIB_BUILD_NAME) $(OFILES) $(LDFLAGS) $(LIBS) $(ELF_DEP_LIBS)
# 	if [ $(ELFLIB_MAJOR) != $(ELFLIB_BUILD_NAME) ]; then \
# 	  cd $(elfdir) && $(LN) $(ELFLIB_BUILD_NAME) $(ELFLIB_MAJOR); \
# 	fi
# 	cd $(elfdir) && $(LN) $(ELFLIB_BUILD_NAME) $(ELFLIB)

$(ELFDIR_F_ELF_LIB):	$(elfdir) $(ELFOBJS) $(ELFOBJS_F)
	-$(RM) -f $(ELFDIR_F_ELF_LIB) $(elfdir)/$(ELFLIB_F_MAJOR) $(elfdir)/$(ELFLIB_F)
	cd $(elfdir) && $(ELF_F_LINK_CMD) -o $(ELFLIB_F_BUILD_NAME) $(OFILES) $(OFILES_F) $(LDFLAGS) $(LIBS_F) $(ELF_DEP_LIBS)
	if [ $(ELFLIB_F_MAJOR) != $(ELFLIB_F_BUILD_NAME) ]; then \
	  cd $(elfdir) && $(LN) $(ELFLIB_F_BUILD_NAME) $(ELFLIB_F_MAJOR); \
	fi
	cd $(elfdir) && $(LN) $(ELFLIB_F_BUILD_NAME) $(ELFLIB_F)

gler_test:	$(ELFDIR_ELF_F_LIB) gler_test.f90
	$(FC) $(FCFLAGS) -o gler_test gler_test.f90 -I$(elfdir) -L$(elfdir) $(MODS_F) $(LIBS) $(GLER_LIBPATHS) -lgler_f $(LIBS_F)

.PHONY: clean distclean tags all install check \
        $(EXECUTABLES)

# --------------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------------

incdir = $(prefix)/include
libdir = $(prefix)/lib
bindir = $(prefix)/bin

install-header:
	$(MKINSDIR) $(incdir)
	@echo installing headers in $(incdir)/
	$(INSTALL_DATA) $(elfdir)/gler_module.mod $(incdir)

install-execs:
	$(MKINSDIR) $(bindir)
	@echo installing $(EXECUTABLES) in $(bindir)/
	for exec in $(EXECUTABLES) ; do \
	   $(INSTALL_PROGRAM) $${exec} $(bindir)/$${exec} ; \
	done

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

install-elf-f:
	$(MKINSDIR) $(libdir)
	-$(RM) $(libdir)/$(ELFLIB_F)
	-$(RM) $(libdir)/$(ELFLIB_F_MAJOR)
	@echo installing $(ELFLIB_F_BUILD_NAME) in $(libdir)
	$(INSTALL_DATA) $(ELFDIR_F_ELF_LIB) $(libdir)
	-cd $(libdir) && $(LN) $(ELFLIB_F_BUILD_NAME) $(ELFLIB_F_MAJOR)
	-cd $(libdir) && $(LN) $(ELFLIB_F_MAJOR) $(ELFLIB_F)
	@echo ""
	@echo $(ELFLIB_F_BUILD_NAME) created in $(libdir)/.
	@echo The links $(ELFLIB_F) and $(ELFLIB_F_MAJOR) to $(ELFLIB_F_BUILD_NAME) were also created.
	@echo ""

install:	elf install-elf-f install-header # install-elf install-execs

# --------------------------------------------------------------------------
# Housekeeping
# --------------------------------------------------------------------------

clean:
	-$(RM) -f .##* $(elfdir)/*
distclean: clean
	-$(RM) -f $(EXECUTABLES)
	-/bin/rmdir $(elfdir)

# --------------------------------------------------------------------------
# Object Rules
# --------------------------------------------------------------------------

$(elfdir)/gler_util.o:	$(srcdir)/gler_util.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/gler_util.c

$(elfdir)/gler_module.o:	$(srcdir)/gler_module.f90
	cd $(elfdir) && $(ELF_F_COMPILE_CMD) $(MODS_F) $(srcdir)/gler_module.f90

