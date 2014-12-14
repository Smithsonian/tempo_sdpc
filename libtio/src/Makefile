# -*- make -*-

include ../../Makefile.defines

THIS_LIB = tio

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/objs

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib

INCS = -I$(config_dir)/include -I$(TELL_INC) -I$(NETCDF_INC)
LIBS = $(INSTALL_LIBPATHS) \
       $(NETCDF_LIBPATHS) -lnetcdf \
       $(TELL_LIBPATHS) -ltell $(RPATH)

MODS_F = -I$(NETCDF_INC) -I$(TELL_INC)
LIBS_F = $(NETCDF_LIBPATHS) -lnetcdff $(TELL_LIBPATHS) -ltell_f

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

OFILES   = tio.o tio_template.o tio_l1_radiance.o tio_l1_irradiance.o
ELFOBJS  = $(elfdir)/tio.o $(elfdir)/tio_template.o $(elfdir)/tio_l1_radiance.o $(elfdir)/tio_l1_irradiance.o
TIO_TEMPLATE_INCLUDE = $(srcdir)/tio_template.h

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
elf:	tempo-names $(elfdir) $(ELFDIR_ELF_LIB) $(ELFDIR_F_ELF_LIB)

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
	cd $(elfdir) && $(ELF_F_LINK_CMD) -o $(ELFLIB_F_BUILD_NAME) $(OFILES_F) $(LDFLAGS) $(LIBS_F) $(ELF_DEP_LIBS)
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
	$(INSTALL_DATA) $(elfdir)/tio_module.mod $(incdir)

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

install:	elf install-header install-elf install-elf-f

# --------------------------------------------------------------------------
# Housekeeping
# --------------------------------------------------------------------------

clean:	clean-tempo-names
	-$(RM) -f .##* $(elfdir)/*
distclean: clean distclean-tempo-names
	-$(RM) -rf $(elfdir)

# --------------------------------------------------------------------------
# Object Rules
# --------------------------------------------------------------------------

NAME_FILE_BASENAMES = tempo_dims tempo_grps tempo_vars

tempo-names:
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("#define %25s   %s\n",a[1],a[2])}' \
	        $${bn}.csv > .$${bn}_h.tmp ; \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("character (len=tiof_max_name_len), public, parameter :: %25s = %s\n",tolower(a[1]),a[2])}' \
	        $${bn}.csv > _$${bn}.inc ; \
	   sed -i -e "1i! Generated from $${bn}.csv by Makefile. Do not edit." \
	        _$${bn}.inc ; \
	done
	sed -e '/@CSV_DIMS@/{r .tempo_dims_h.tmp' -e 'd}' \
	    -e '/@CSV_GRPS@/{r .tempo_grps_h.tmp' -e 'd}' \
	    -e '/@CSV_VARS@/{r .tempo_vars_h.tmp' -e 'd}' \
	    $(srcdir)/tio_template.hin > $(TIO_TEMPLATE_INCLUDE)
	sed -i -e '1i/* Generated from tio_template.hin by Makefile. Do not edit. */' \
	     $(TIO_TEMPLATE_INCLUDE)

clean-tempo-names:
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f .$${bn}_h.tmp ; \
	done
distclean-tempo-names:	clean-tempo-names
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f _$${bn}.inc $(TIO_TEMPLATE_INCLUDE) ; \
	done

$(elfdir)/tio.o:	$(srcdir)/tio.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tio.c $(LIBS)

$(elfdir)/tio_l1_radiance.o:	$(srcdir)/tio_l1_radiance.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tio_l1_radiance.c $(LIBS)

$(elfdir)/tio_l1_irradiance.o:	$(srcdir)/tio_l1_irradiance.c
	cd $(elfdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/tio_l1_irradiance.c $(LIBS)

$(elfdir)/tio_module.o:	$(srcdir)/tio_module.f90
	cd $(elfdir) && $(ELF_F_COMPILE_CMD) $(MODS_F) $(srcdir)/tio_module.f90
