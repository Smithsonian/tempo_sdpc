# -*- sh -*-

EXECUTABLE = L1_cloud_o4

include ../../Makefile.defines

srcdir     = $(config_dir)/L1_cloud_o4/src
objdir     = $(srcdir)/$(archobjs)
msgdir     = $(srcdir)/msgs

incdir     = $(prefix)/include
bindir     = $(prefix)/bin
libdir     = $(prefix)/lib

ADDLIBD    = #
ADDLIBS    = $(RPATH) -lsz -lm $(DL_LIB) # -llapack -lblas

#---------------------------------------------------------------
# Ideally, you shouldn't need to change anything below this line
#---------------------------------------------------------------

PGSHOME     = $(OMIUTIL)/TOOLKIT

MV = /bin/mv
RM = /bin/rm
CP = /bin/cp

TIO_LIB  = $(libdir)
TIO_SRC  = $(config_dir)/libtio/src
TIO_SRC_OBJDIR = $(TIO_SRC)/$(archobjs)

MD_LIB  = $(libdir)
MD_SRC  = $(config_dir)/libmd/src
MD_SRC_OBJDIR = $(MD_SRC)/$(archobjs)

CLIM_LIB = $(libdir)
CLIM_SRC = $(config_dir)/libclim/src
CLIM_SRC_OBJDIR = $(CLIM_SRC)/$(archobjs)

INTERP_LIB = $(libdir)
INTERP_SRC = $(config_dir)/libndinterp/src
INTERP_SRC_OBJDIR = $(INTERP_SRC)/$(archobjs)

GLER_LIB = $(libdir)
GLER_SRC = $(config_dir)/libgler/src
GLER_SRC_OBJDIR = $(GLER_SRC)/$(archobjs)

YAML_LIB = $(OTS_LIBDIR)
YAML_INC = $(OTS_INCDIR)

OTS_INCDIR = $(OTS_ROOT)/include
OTS_MODDIR = $(OTS_ROOT)/mod
OTS_LIBDIR = $(OTS_ROOT)/lib

LOCAL_OBJDIRS = -L$(TIO_SRC_OBJDIR) \
                -L$(MD_SRC_OBJDIR) \
                -L$(CLIM_SRC_OBJDIR) \
                -L$(INTERP_SRC_OBJDIR) \
                -L$(GLER_SRC_OBJDIR)

PSPLINELIBS  = -lpspline -lezcdf
TIOLIBS      = -L$(TIO_LIB) -ltio -ltio_f \
	       -L$(MD_LIB) -lmd_f \
               -L$(TELL_LIB) -ltell -ltell_f
NETCDFLIBS   = -lnetcdf -lnetcdff
HE4LIBS      = -lhdfeos -lmfhdf -lGctp -ldf -ljpeg -lz
HE5LIBS      = -lhe5_hdfeos -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lhdf5
PGSLIBS      = -lPGSTK
YAMLLIBS     = -lyaml
CLIM_LIBS    = -lcsv -lioclib -lclim_f
INTERPLIBS   = -lndinterp_f
GLERLIBS     = -lgler_f

PSPLINEMOD = $(OTS_MODDIR)
PSPLINELIB = $(OTS_LIBDIR)

IFLAGS = $(MODFLG)$(HDFEOS_INC) \
         $(MODFLG)$(HDF5INC) \
         $(MODFLG)$(HDFINC) \
         $(MODFLG)$(PGSINC) \
         $(MODFLG)$(HDFEOS5_INC) \
         $(MODFLG)$(objdir) \
         $(MODFLG)$(TELL_INC) \
         $(MODFLG)$(NETCDF_INC) \
	 $(MODFLG)$(PSPLINEMOD) \
         $(MODFLG)$(OTS_MODDIR) $(MODFLG)$(OTS_INCDIR) \
         $(MODFLG)$(TIO_SRC_OBJDIR) \
         $(MODFLG)$(MD_SRC_OBJDIR) \
         $(MODFLG)$(CLIM_SRC_OBJDIR) \
         $(MODFLG)$(GLER_SRC_OBJDIR) \
         $(MODFLG)$(INTERP_SRC_OBJDIR)

LFLAGS = $(LOCAL_OBJDIRS) -L$(objdir) \
         -Wl,-rpath=$(libdir) -L$(libdir) \
         -Wl,-rpath=$(OTS_LIBDIR) -L$(OTS_LIBDIR) \
         -L$(HDFEOS_LIB) -L$(HDFLIB) -L$(PGSLIB) -L$(HDF5LIB) -L$(HDFEOS5_LIB) \
         -L$(NETCDF_LIB) -L$(PSPLINELIB) $(ADDLIBD)

LIBS   = $(PSPLINELIBS) $(TIOLIBS) $(CLIM_LIBS) $(INTERPLIBS) $(GLERLIBS) $(NETCDFLIBS) \
         $(PGSLIBS) $(HE5LIBS) $(HE4LIBS) $(YAMLLIBS) $(ADDLIBS)

CINCFLGS = -I$(PGSINC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS5_INC) -I$(YAML_INC) \
           -I$(objdir) -Dgfortran

# During development, make it easy to change the list of source files.
# For releases, use a static list of source files.
F90SOURCES = $(shell grep -v '\#' f90_files.lis)
CSOURCES   = $(shell grep -v '\#' c_files.lis)
LINKOBJS   = $(addprefix $(objdir)/,$(F90SOURCES:.f90=.o) $(CSOURCES:.c=.o))

#---------------------------------
# Primary 'make' targets for users
#   make [all|install|clean]
#---------------------------------

all:	dirs messages $(EXECUTABLE)

dirs:
	mkdir -p $(objdir) $(msgdir)

$(EXECUTABLE):	$(LINKOBJS)
	$(FC) OMCDO2N.f90 -o $(objdir)/$(EXECUTABLE) $(LINKOBJS) $(LDFLAGS) $(FCFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS)

install-executable:
	mkdir -p $(bindir)
	@echo installing $(EXECUTABLE) in $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/$(EXECUTABLE) $(bindir)/$(EXECUTABLE)

install-header:
	mkdir -p $(incdir)
	@echo installing $(NAME_INCLUDE) in $(incdir)/
	$(INSTALL_DATA) $(srcdir)/$(NAME_INCLUDE) $(incdir)/$(NAME_INCLUDE)

install:	install-executable # install-header

clean:
	-$(RM) $(objdir)/*

distclean:	clean
	-$(RM) -rf $(objdir)

tags:
	ctags $(F90SOURCES) $(CSOURCES)
#-----------------------------------------------------------------
# The hdf.f90 file contains numerous non-conforming tabs.
# It's unlikely anyone will ever fix that in the HDF distribution,
# so we'll just fix it this way:
#-----------------------------------------------------------------
fix-hdf-tabs:
	sed --in-place=_ 's,\t, ,g' $(OTS_ROOT)/include/hdf.f90

# ---------------------------------------
# The PGE and TOOLKIT error message files
# ---------------------------------------

PGEPGSSRC = OMI_52090.t
PGEPGSCMP = $(msgdir)/OMI_52090
#TKPGSCMP = PGS_9 PGS_1

messages: $(PGEPGSCMP) #tkmessages

tkmessages:
	@for tkpgs in $(TKPGSCMP); do                                 \
	  if [ ! -e $(msgdir)/$$tkpgs ]; then                        \
	     echo 'Copying' $$tkpgs 'to' $(msgdir);                  \
	     touch $(msgdir)/$$tkpgs ; chmod u+w $(msgdir)/$$tkpgs;  \
	     $(CP) -p $(PGSHOME)/message/$$tkpgs $(msgdir)/;            \
	  fi;                                                        \
	done

# ----------------------------------------------------------------------
# Pattern rule for the creation of SMF messages and associated .f and .h
# ----------------------------------------------------------------------
$(msgdir)/OMI_52090:
	export PGSMSG="$(msgdir)" ; \
	$(PGSBIN)/smfcompile -r -f77 -f OMI_52090.t ; \
	$(PGSBIN)/smfcompile -r      -f OMI_52090.t

# ---- build dependencies ----
# During development, the rules will be dynamically updated.
# For releases, a static set of rules can be appended.
Makefile:	rules.inc
rules.inc:	$(F90SOURCES) $(CSOURCES) f90_files.lis c_files.lis
ifeq (, $(shell which slsh))
   $(error "slsh is not available on $$PATH, please install S-Lang")
endif
	$(DEPMOD) -r $(F90SOURCES) > rules.inc
	./mkrules.sh c $(CSOURCES) >> rules.inc

ODIR = $(objdir)
SDIR = $(srcdir)
-include rules.inc

.PHONY: tags $(NAME_FILE_TARGETS)
