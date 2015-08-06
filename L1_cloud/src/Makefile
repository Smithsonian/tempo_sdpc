# -*- sh -*-

EXECUTABLE = L1_cloud

include ../../Makefile.defines

srcdir     = $(config_dir)/L1_cloud/src
objdir     = $(srcdir)/$(archobjs)
msgdir     = $(srcdir)/msgs

incdir     = $(prefix)/include
bindir     = $(prefix)/bin

ADDLIBD    = #
ADDLIBS    = $(RPATH) -lsz -lm $(DL_LIB)

#---------------------------------------------------------------
# Ideally, you shouldn't need to change anything below this line
#---------------------------------------------------------------

PGSHOME     = $(OMIUTIL)/TOOLKIT

MV = /bin/mv
RM = /bin/rm
CP = /bin/cp

TIO_LIB  = $(prefix)/lib
TIO_SRC  = $(config_dir)/libtio/src
TIO_SRC_OBJDIR = $(TIO_SRC)/$(archobjs)

PSPLINELIBS  = -lpspline -lezcdf
TIOLIBS      = -Wl,-rpath=$(TIO_SRC_OBJDIR) -L$(TIO_SRC_OBJDIR) \
               -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB) -ltio -ltio_f \
               -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB) -ltell -ltell_f
NETCDFLIBS   = -lnetcdf -lnetcdff
HE4LIBS      = -lhdfeos -lmfhdf -lGctp -ldf -ljpeg -lz
HE5LIBS      = -lhe5_hdfeos -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lhdf5
PGSLIBS      = -lPGSTK
OMILIBS      =

OTS_INCDIR = $(OTS_ROOT)/include
OTS_MODDIR = $(OTS_ROOT)/mod
OTS_LIBDIR = $(OTS_ROOT)/lib

PSPLINEMOD = $(OTS_MODDIR)
PSPLINELIB = $(OTS_LIBDIR)

IFLAGS = $(MODFLG)$(HDFEOS_INC) \
         $(MODFLG)$(HDF5INC) \
         $(MODFLG)$(HDFINC) \
         $(MODFLG)$(PGSINC) \
         $(MODFLG)$(HDFEOS5_INC) \
         $(MODFLG)$(objdir) \
         $(MODFLG)$(TIO_SRC_OBJDIR) \
         $(MODFLG)$(TELL_INC) \
         $(MODFLG)$(NETCDF_INC) \
	 $(MODFLG)$(PSPLINEMOD)

LFLAGS = -L$(objdir) \
         -L$(HDFEOS_LIB) -L$(HDFLIB) -L$(PGSLIB) -L$(HDF5LIB) -L$(HDFEOS5_LIB) \
         -L$(NETCDF_LIB) -L$(PSPLINELIB) $(ADDLIBD)

LIBS   = $(PSPLINELIBS) $(TIOLIBS) $(NETCDFLIBS) \
         $(OMILIBS) $(PGSLIBS) $(HE5LIBS) $(HE4LIBS) $(ADDLIBS)

CINCFLGS = -I$(PGSINC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS5_INC) -I$(objdir)

# During development, make it easy to change the list of source files.
# For releases, use a static list of source files.
F90SOURCES = $(shell grep -v '\#' f90_files.lis)
CSOURCES   = $(shell grep -v '\#' c_files.lis)
LINKOBJS   = $(addprefix $(objdir)/,$(F90SOURCES:.f90=.o) $(CSOURCES:.c=.o))

#---------------------------------
# Primary 'make' targets for users
#   make [all|install|clean]
#---------------------------------

all:	dirs messages cld-names $(EXECUTABLE)

dirs:
	mkdir -p $(objdir) $(msgdir)

$(EXECUTABLE):	$(LINKOBJS)
	$(FC) OMCLDRR.f90 -o $(objdir)/$(EXECUTABLE) $(LINKOBJS) $(LDFLAGS) $(FCFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS)

install-executable:
	mkdir -p $(bindir)
	@echo installing $(EXECUTABLE) in $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/$(EXECUTABLE) $(bindir)/$(EXECUTABLE)

install-header:
	mkdir -p $(incdir)
	@echo installing $(NAME_INCLUDE) in $(incdir)/
	$(INSTALL_DATA) $(srcdir)/$(NAME_INCLUDE) $(incdir)/$(NAME_INCLUDE)

install:	install-executable install-header

clean:  clean-cld-names
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

PGEPGSSRC = OMCLDRR_52251.t OMI_1900.t
PGEPGSCMP = $(msgdir)/PGS_OMCLDRR_52251 $(msgdir)/PGS_OMI_1900
TKPGSCMP = PGS_9 PGS_1

messages: $(PGEPGSCMP) tkmessages

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
$(msgdir)/PGS_OMI_1900:
	export PGSMSG="$(msgdir)" ; \
	$(PGSBIN)/smfcompile -r -f77 -f OMI_1900.t ; \
	$(PGSBIN)/smfcompile -r      -f OMI_1900.t
$(msgdir)/PGS_OMCLDRR_52251:
	export PGSMSG="$(msgdir)" ; \
	$(PGSBIN)/smfcompile -r -f77 -f OMCLDRR_52251.t ; \
	$(PGSBIN)/smfcompile -r      -f OMCLDRR_52251.t

#
# Name tables
#
NAME_FILE_BASENAMES = cld_names_dim cld_names_var cld_names_grp
NAME_INCLUDE = cld_names.h

NAME_FILE_SOURCES = $(addsuffix .csv, $(NAME_FILE_BASENAMES))
NAME_FILE_TARGETS = $(addsuffix .inc, $(NAME_FILE_BASENAMES)) $(NAME_INCLUDE)
cld-names:	$(NAME_FILE_TARGETS)
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("#define %35s   %s\n",a[1],a[2])}' \
	        $${bn}.csv > .$${bn}_h.tmp ; \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("character (len=cld_max_name_len), public, parameter :: %35s = %s\n",tolower(a[1]),a[2])}' \
	        $${bn}.csv > $${bn}.inc ; \
	   sed -i -e "1i! Autogenerated from $${bn}.csv by Makefile. Do not edit." \
	        $${bn}.inc ; \
	done
	sed -e '/@CSV_DIMS@/{r .cld_names_dim_h.tmp' -e 'd}' \
	    -e '/@CSV_VARS@/{r .cld_names_var_h.tmp' -e 'd}' \
	    $(srcdir)/$(NAME_INCLUDE)in > $(NAME_INCLUDE)
	sed -i -e '1i/*! @brief Autogenerated by Makefile. Do not edit. \n *  @file */' \
	     $(NAME_INCLUDE)

clean-cld-names:
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f .$${bn}_h.tmp ; \
	done
distclean-cld-names:	clean-cld-names
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f $${bn}.inc $(NAME_INCLUDE) ; \
	done

# ---- build dependencies ----
# During development, the rules will be dynamically updated.
# For releases, a static set of rules can be appended.
Makefile:	rules.inc
rules.inc:	$(F90SOURCES) $(CSOURCES) f90_files.lis c_files.lis
	$(DEPMOD) -r $(F90SOURCES) > rules.inc
	./mkrules.sh c $(CSOURCES) >> rules.inc

ODIR = $(objdir)
SDIR = $(srcdir)
-include rules.inc

.PHONY: tags cld-names $(NAME_FILE_TARGETS)
