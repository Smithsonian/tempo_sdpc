# -*- make -*-

include ../../Makefile.defines

THIS_LIB = met

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/$(archobjs)

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
TIO_LIBPATHS     = -Wl,-rpath=$(config_dir)/libtio/src/$(archobjs) -L$(config_dir)/libtio/src/$(archobjs)
TIO_INC          = $(config_dir)/libtio/src/$(archobjs)

ECCODES_LIBPATH  = -Wl,-rpath=/data/tempo_gs/sdpc/install/v1/ots/lib -L/data/tempo_gs/sdpc/install/v1/ots/lib

MET_LIBPATHS     = -Wl,-rpath=$(elfdir) -L$(elfdir)

INCS = -I$(TELL_INC) -I$(NETCDF_INC) -I$(TIO_INC)
LIBS = $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell \
       $(ECCODES_LIBPATH) -leccodes \
       $(TIO_LIBPATHS) -ltio $(RPATH) $(INSTALL_LIBPATHS)

MODS_F = -I$(NETCDF_INC) -I$(TELL_INC) -I($TIO_INC)
LIBS_F = $(NETCDF_LIBPATHS) -lnetcdff $(TELL_LIBPATHS) -ltell_f $(TIO_LIBPATHS) -ltio_f $(ECCODES_LIBPATH) -leccodes_f90

ELF_MAJOR_VERSION = 0#
ELF_MINOR_VERSION = 0#
ELF_MICRO_VERSION = 1#

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

OFILES_F = met_module.o
ELFOBJS_F = $(elfdir)/met_module.o

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	elf
elf:	$(elfdir) $(ELFDIR_F_ELF_LIB)

tags:
	ctags $(patsubst %.o,%.f90,$(OFILES_F))

$(elfdir):
	-$(MKINSDIR) $(elfdir)

$(ELFDIR_F_ELF_LIB):	$(elfdir) $(ELFOBJS_F)
	-$(RM) -f $(ELFDIR_F_ELF_LIB) $(elfdir)/$(ELFLIB_F_MAJOR) $(elfdir)/$(ELFLIB_F)
	cd $(elfdir) && $(ELF_F_LINK_CMD) -o $(ELFLIB_F_BUILD_NAME) $(OFILES_F) $(LDFLAGS) $(LIBS_F) $(ELF_DEP_LIBS)
	if [ $(ELFLIB_F_MAJOR) != $(ELFLIB_F_BUILD_NAME) ]; then \
	  cd $(elfdir) && $(LN) $(ELFLIB_F_BUILD_NAME) $(ELFLIB_F_MAJOR); \
	fi
	cd $(elfdir) && $(LN) $(ELFLIB_F_BUILD_NAME) $(ELFLIB_F)

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
	$(INSTALL_DATA) $(elfdir)/met_module.mod $(incdir)

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

install:	elf install-header install-elf-f

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


$(elfdir)/met_module.o:	$(srcdir)/met_module.f90
	git_commit="(unknown)" ; \
	if test -d $(config_dir)/.git ; then \
	   git_commit=`git rev-parse --verify HEAD` ; \
	fi ; \
	cd $(elfdir) && $(ELF_F_COMPILE_CMD) $(INCS) -DGIT_COMMIT_HASH='"'$${git_commit}'"' $(srcdir)/met_module.f90 $(LIBS)

#$(elfdir)/tio_ppc.o:	$(srcdir)/tio_ppc.c
#	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tio_ppc.c $(LIBS)
