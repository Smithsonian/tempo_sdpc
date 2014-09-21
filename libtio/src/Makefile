# -*- make -*-

include ../../Makefile.defines

THIS_LIB = tio

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/objs

INCS = -I$(config_dir)/include -I$(NETCDF_INC)
LIBS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB) -lnetcdf

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

OFILES   = tio.o tio_l1.o tio_err.o
ELFOBJS  = $(elfdir)/tio.o $(elfdir)/tio_l1.o $(elfdir)/tio_err.o

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

OFILES_F = tio_module.o
ELFOBJS_F = $(elfdir)/tio_module.o

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	elf
elf:	$(elfdir) $(ELFDIR_ELF_LIB) # $(ELFDIR_F_ELF_LIB)

$(elfdir):
	-$(MKINSDIR) $(elfdir)

$(ELFDIR_ELF_LIB):	$(elfdir) $(ELFOBJS)
	-$(RM) -f $(ELFDIR_ELF_LIB) $(elfdir)/$(ELFLIB_MAJOR) $(elfdir)/$(ELFLIB)
	cd $(elfdir) && $(ELF_LINK_CMD) -o $(ELFLIB_BUILD_NAME) $(OFILES) $(LDFLAGS) $(LIBS) $(ELF_DEP_LIBS)
	if [ $(ELFLIB_MAJOR) != $(ELFLIB_BUILD_NAME) ]; then \
	  cd $(elfdir) && $(LN) $(ELFLIB_BUILD_NAME) $(ELFLIB_MAJOR); \
	fi
	cd $(elfdir) && $(LN) $(ELFLIB_BUILD_NAME) $(ELFLIB)

$(ELFDIR_F_ELF_LIB):	$(elfdir) $(ELFOBJS_F)
	-$(RM) -f $(ELFDIR_F_ELF_LIB) $(elfdir)/$(ELFLIB_F_MAJOR) $(elfdir)/$(ELFLIB_F)
	cd $(elfdir) && $(ELF_F_LINK_CMD) -o $(ELFLIB_F_BUILD_NAME) $(OFILES_F) $(LDFLAGS) $(LIBS) $(ELF_DEP_LIBS)
	if [ $(ELFLIB_F_MAJOR) != $(ELFLIB_F_BUILD_NAME) ]; then \
	  cd $(elfdir) && $(LN) $(ELFLIB_F_BUILD_NAME) $(ELFLIB_F_MAJOR); \
	fi
	cd $(elfdir) && $(LN) $(ELFLIB_F_BUILD_NAME) $(ELFLIB_F)

.PHONY: clean distclean all install check

# --------------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------------

incdir = $(prefix)/include
libdir = $(prefix)/lib

install-header:
	$(MKINSDIR) $(incdir)
	@echo installing tio.h and tio_module.mod in $(incdir)/
	$(INSTALL_DATA) $(srcdir)/tio.h $(incdir)
#	$(INSTALL_DATA) $(elfdir)/tio_module.mod $(incdir)

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

install-elf:	elf install-header install-elf # install-elf-f
install:	install-elf

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

$(elfdir)/tio.o:	$(srcdir)/tio.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tio.c $(LIBS)

$(elfdir)/tio_l1.o:	$(srcdir)/tio_l1.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tio_l1.c $(LIBS)

$(elfdir)/tio_err.o:	$(srcdir)/tio_err.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tio_err.c $(LIBS)

# $(elfdir)/tio_module.o:	$(srcdir)/tio_module.f90
#	cd $(elfdir) && $(ELF_F_COMPILE_CMD) $(srcdir)/tio_module.f90
