# -*- sh -*-

EXECUTABLE = L1_o3_profile

MERGE_EXEC = merge_o3p_files

include ../../Makefile.defines

srcdir     = $(config_dir)/L1_o3_profile/src
objdir     = $(srcdir)/$(archobjs)
msgdir     = $(srcdir)/msgs

incdir     = $(prefix)/include
bindir     = $(prefix)/bin
msgs_installdir = $(prefix)/msgs

ADDLIBD    = -Wl,-rpath=$(OTS_ROOT)/lib -L$(OTS_ROOT)/lib
ADDLIBS    = $(RPATH) -lsz -lm $(DL_LIB) -llidort_vv2p7_pcachrisbrdf -llapack -lblas

#---------------------------------------------------------------
# Ideally, you shouldn't need to change anything below this line
#---------------------------------------------------------------

MV = /bin/mv
RM = /bin/rm
CP = /bin/cp

TIO_LIB  = $(prefix)/lib
TIO_SRC  = $(config_dir)/libtio/src
TIO_SRC_OBJDIR = $(TIO_SRC)/$(archobjs)
TIOINCS     = -I$(TIOSRC_OBJS)
TIOLIBS     = -L$(TIO_LIB) -ltio_f -ltio

LUTSAO_DIR = $(OTS_ROOT)/lib/lutsao-1.0.0
LUTSAO_LIB = -L$(LUTSAO_DIR) -lLUT_SAO

GLER_LIB     = $(prefix)/lib
GLER_SRC     = $(config_dir)/libgler/src
GLERSRC_OBJS = $(GLER_SRC)/$(archobjs)
GLERINCS     = -I$(GLERSRC_OBJS)
GLERLIBS     = -L$(GLERSRC_OBJS) -L$(GLER_LIB) -lgler_f $(NDINTERP_LIB)

CLIM_LIB     = $(prefix)/lib
CLIM_SRC     = $(config_dir)/libclim/src
CLIMSRC_OBJS = $(CLIM_SRC)/$(archobjs)
CLIMINCS     = -I$(CLIMSRC_OBJS)
CLIMLIBS     = -L$(CLIMSRC_OBJS) -L$(CLIM_LIB) -lclim_f -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -lcsv

MET_LIB     = $(prefix)/lib
MET_SRC     = $(config_dir)/libmet/src
METSRC_OBJS = $(MET_SRC)/$(archobjs)
METINCS     = -I$(METSRC_OBJS)
METLIBS     = -L$(METSRC_OBJS) -L$(MET_LIB) -lmet_f -lmet

MD_LIB      = $(prefix)/lib
MD_SRC      = $(config_dir)/libmd/src
MDSRC_OBJS  = $(MD_SRC)/$(archobjs)
MDINCS      = -I$(MDSRC_OBJS)
MDLIBS      = -L$(MDSRC_OBJS) -L$(MD_LIB) -lmd_f

PSPLINELIBS  = -lpspline -lezcdf

TIOLIBS      = -L$(TIO_SRC_OBJDIR) \
               -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB) -ltio -ltio_f \
               -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB) -ltell -ltell_f \

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
	 $(MODFLG)$(PSPLINEMOD) \
         $(MODFLG)$(LUTSAO_DIR) \
	 $(CLIMINCS) \
	 $(GLERINCS) \
	 $(METINCS) \
         $(MDINCS)

LFLAGS = -L$(objdir) \
         -L$(HDFEOS_LIB) -L$(HDFLIB) -L$(PGSLIB) -L$(HDF5LIB) -L$(HDFEOS5_LIB) \
         -L$(NETCDF_LIB) -L$(PSPLINELIB) $(ADDLIBD)

LIBS   = $(PSPLINELIBS) $(TIOLIBS) $(NETCDFLIBS) \
         $(OMILIBS) $(PGSLIBS) $(HE5LIBS) $(HE4LIBS) $(ADDLIBS) \
	 $(CLIMLIBS) $(GLERLIBS) $(METLIBS) $(MDLIBS) $(LUTSAO_LIB)

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

all:	dirs messages o3p-names $(EXECUTABLE) $(MERGE_EXEC)
#all:	dirs messages $(EXECUTABLE)

dirs:
	mkdir -p $(objdir) $(msgdir)

$(EXECUTABLE):	$(LINKOBJS)
	$(FC) PROFOZ_main.f90 -o $(objdir)/$(EXECUTABLE) $(LINKOBJS) $(LDFLAGS) $(FCFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS)

$(MERGE_EXEC):	$(LINKOBJS)
	$(FC) merge_o3p_files.f90 -o $(objdir)/$(MERGE_EXEC) $(LINKOBJS) $(LDFLAGS) $(FCFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS)

install-executable:
	mkdir -p $(bindir)
	@echo installing $(EXECUTABLE) in $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/$(EXECUTABLE) $(bindir)/$(EXECUTABLE)

install-merge-exec:
	mkdir -p $(bindir)
	@echo installing $(MERGE_EXEC) in $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/$(MERGE_EXEC) $(bindir)/$(MERGE_EXEC)

install-header:
	mkdir -p $(incdir)
	@echo installing $(NAME_INCLUDE) in $(incdir)/
	$(INSTALL_DATA) $(srcdir)/$(NAME_INCLUDE) $(incdir)/$(NAME_INCLUDE)

install-msgs:
	mkdir -p $(msgs_installdir)
	@echo installing message files in $(msgs_installdir)/
	for m in $(MESSAGE_FILES) ; do \
	   $(INSTALL_DATA) $(srcdir)/msgs/$$m $(msgs_installdir)/$$m ; \
	done

install:	install-executable install-merge-exec install-header install-msgs

clean:  clean-o3p-names
	-$(RM) -f $(objdir)/*

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
PGSHOME = $(OMIUTIL)/TOOLKIT

MESSAGE_FILES = PGS_9 PGS_3 PGS_13 PGS_1900 PGS_52510

PGEPGSSRC = OMSAO_52510.t OMI_1900.t
PGEPGSCMP = $(msgdir)/PGS_OMSAO_52510 $(msgdir)/PGS_OMI_1900
TKPGSCMP = PGS_9 PGS_3 PGS_13

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
$(msgdir)/PGS_OMSAO_52510:
	export PGSMSG="$(msgdir)" ; \
	$(PGSBIN)/smfcompile -r -f77 -f OMSAO_52510.t ; \
	$(PGSBIN)/smfcompile -r      -f OMSAO_52510.t

#
# Name tables
#
NAME_FILE_BASENAMES = o3p_names_dim o3p_names_var o3p_names_grp
NAME_INCLUDE = o3p_names.h

NAME_FILE_SOURCES = $(addsuffix .csv, $(NAME_FILE_BASENAMES))
NAME_FILE_TARGETS = $(addsuffix .inc, $(NAME_FILE_BASENAMES)) $(NAME_INCLUDE)
o3p-names:	$(NAME_FILE_TARGETS)
$(NAME_FILE_TARGETS):	$(NAME_FILE_SOURCES)
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("#define %35s   %s\n",a[1],a[2])}' \
	        $${bn}.csv > .$${bn}_h.tmp ; \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("character (len=o3p_max_name_len), public, parameter :: %35s = %s\n",tolower(a[1]),a[2])}' \
	        $${bn}.csv > $${bn}.inc ; \
	   sed -i -e "1i! Autogenerated from $${bn}.csv by Makefile. Do not edit." \
	        $${bn}.inc ; \
	done
	sed -e '/@CSV_DIMS@/{r .o3p_names_dim_h.tmp' -e 'd}' \
	    -e '/@CSV_VARS@/{r .o3p_names_var_h.tmp' -e 'd}' \
	    $(srcdir)/$(NAME_INCLUDE)in > $(NAME_INCLUDE)
	sed -i -e '1i/*! @brief Autogenerated by Makefile. Do not edit. \n *  @file */' \
	     $(NAME_INCLUDE)

clean-o3p-names:
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f .$${bn}_h.tmp ; \
	done
distclean-o3p-names:	clean-o3p-names
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f $${bn}.inc $(NAME_INCLUDE) ; \
	done

# ---- build dependencies ----
# During development, the rules will be dynamically updated.
# For releases, a static set of rules can be appended.
Makefile:	rules.inc
rules.inc:	$(F90SOURCES) $(CSOURCES) f90_files.lis c_files.lis
#rules.inc:	f90_files.lis c_files.lis
ifeq (, $(shell which slsh))
   $(error "slsh is not available on $$PATH, please install S-Lang")
endif
	$(DEPMOD) -r $(F90SOURCES) > rules.inc
	./mkrules.sh c $(CSOURCES) >> rules.inc

ODIR = $(objdir)
SDIR = $(srcdir)
-include rules.inc

#.PHONY: tags
.PHONY: tags o3p-names $(NAME_FILE_TARGETS)
