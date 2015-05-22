# -*- sh -*-

EXECUTABLE = L1_o3_total

include ../../Makefile.defines

srcdir     = $(config_dir)/L1_o3_total/src
objdir     = $(srcdir)/$(archobjs)
msgdir     = $(srcdir)/msgs

bindir = $(prefix)/bin
incdir = $(prefix)/include
msgs_installdir = $(prefix)/msgs

ADDLIBD    = #
ADDLIBS    = $(RPATH) -lsz -lm $(DL_LIB)

#---------------------------------------------------------------
# Ideally, you shouldn't need to change anything below this line
#---------------------------------------------------------------

PGSHOME     = $(OMIUTIL)/TOOLKIT

MV = /bin/mv
RM = /bin/rm
CP = /bin/cp

LDPATH  = -Wl,-rpath=$(TIOSRC_OBJS):$(TIO_LIB):$(TELL_LIB) $(RPATH)

TIO_LIB     = $(prefix)/lib
TIO_SRC     = $(config_dir)/libtio/src
TIOSRC_OBJS = $(TIO_SRC)/$(archobjs)
TIOINCS     = -I$(TIOSRC_OBJS)
TIOLIBS     = -L$(TIO_LIB) -L$(TIOSRC_OBJS) -ltio_f -ltio

TELLINCS    = -I$(TELL_INC)
TELLLIBS    = -L$(TELL_LIB) -ltell_f -ltell

PSPLINEMODS = -I$(OTS_ROOT)/mod
PSPLINELIBS = -L$(OTS_ROOT)/lib -lpspline -lezcdf

NETCDFINCS  = -I$(NETCDF_INC)
NETCDFLIBS  = -L$(NETCDF_LIB) -lnetcdf -lnetcdff

HE5INCS     = -I$(HDF5INC) -I$(HDFEOS5_INC)
HE5LIBS     = -L$(HDF5LIB) -L$(HDFEOS5_LIB) -lhe5_hdfeos -lGctp -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lhdf5

HE4INCS     = -I$(HDFINC) -I$(HDFEOS_INC)
HE4LIBS     = -L$(HDFLIB) -L$(HDFEOS_LIB) -lhdfeos -lmfhdf -lGctp -ldf -ljpeg -lz

ZLIBS       = -L$(OTS_ROOT)/lib -lsz -lz

PGSINCS     = -I$(PGSINC)
PGSLIBS     = -L$(PGSLIB) -lPGSTK

IFLAGS = -I$(objdir) $(TIOINCS) $(TELLINCS) $(PSPLINEMODS) \
         $(NETCDFINCS) $(HE5INCS) $(PGSINCS) $(HE4INCS)

LIBS   = $(LDPATH) $(TIOLIBS) $(TELLLIBS) $(PSPLINELIBS) \
         $(NETCDFLIBS) $(HE5LIBS) $(HE4LIBS) $(ZLIBS) $(PGSLIBS) \
         -lm $(DL_LIB)

CINCFLGS = -I$(objdir) $(HE5INCS) $(PGSINCS) $(HE4INCS) $(HDFEOS_INC)

# During development, make it easy to change the list of source files.
# For releases, use a static list of source files.
F90SOURCES = $(shell grep -v '\#' f90_files.lis)
CSOURCES   = $(shell grep -v '\#' c_files.lis)
LINKOBJS   = $(addprefix $(objdir)/,$(F90SOURCES:.f90=.o) $(CSOURCES:.c=.o))

#---------------------------------
# Primary 'make' targets for users
#   make [all|install|clean]
#---------------------------------

all:	dirs messages o3t-names $(EXECUTABLE)

dirs:
	mkdir -p $(objdir) $(msgdir)

$(EXECUTABLE):	$(LINKOBJS)
	$(FC) O3T_main.f90 -o $(objdir)/$(EXECUTABLE) $(LINKOBJS) $(LDFLAGS) $(FCFLAGS) $(IFLAGS) $(LIBS)

install-executable:
	mkdir -p $(bindir)
	@echo installing $(EXECUTABLE) in $(bindir)/        
	$(INSTALL_PROGRAM) $(objdir)/$(EXECUTABLE) $(bindir)/$(EXECUTABLE)

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

install:	install-executable install-header install-msgs

clean:	clean-o3t-names
	-$(RM) $(objdir)/* PGS_OZT_52050.f PGS_OZT_52050.h

distclean:	clean distclean-o3t-names
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

MESSAGE_FILES = PGS_9 PGS_1 PGS_52050

PGEPGSSRC = OZT_52050.t
PGEPGSCMP = $(msgdir)/PGS_OZT_52050
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
$(msgdir)/PGS_OZT_52050:
	export PGSMSG="$(msgdir)" ; \
	$(PGSBIN)/smfcompile -r -f77 -f OZT_52050.t ; \
	$(PGSBIN)/smfcompile -r      -f OZT_52050.t

#
# Name tables
#
NAME_FILE_BASENAMES = o3t_names_dim o3t_names_var o3t_names_grp
NAME_INCLUDE = o3t_names.h

NAME_FILE_SOURCES = $(addsuffix .csv, $(NAME_FILE_BASENAMES))
NAME_FILE_TARGETS = $(addsuffix .inc, $(NAME_FILE_BASENAMES)) $(NAME_INCLUDE)
o3t-names:	$(NAME_FILE_TARGETS)
$(NAME_FILE_TARGETS):	$(NAME_FILE_SOURCES)
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("#define %35s   %s\n",a[1],a[2])}' \
	        $${bn}.csv > .$${bn}_h.tmp ; \
	   awk 'NF && !/^#/ {split($$0,a,","); printf("character (len=o3t_max_name_len), public, parameter :: %35s = %s\n",tolower(a[1]),a[2])}' \
	        $${bn}.csv > $${bn}.inc ; \
	   sed -i -e "1i! Generated from $${bn}.csv by Makefile. Do not edit." \
	        $${bn}.inc ; \
	done
	sed -e '/@CSV_DIMS@/{r .o3t_names_dim_h.tmp' -e 'd}' \
	    -e '/@CSV_VARS@/{r .o3t_names_var_h.tmp' -e 'd}' \
	    $(srcdir)/$(NAME_INCLUDE)in > $(NAME_INCLUDE)
	sed -i -e '1i/* Generated by Makefile. Do not edit. */' \
	     $(NAME_INCLUDE)

clean-o3t-names:
	for bn in $(NAME_FILE_BASENAMES) ; do \
	   $(RM) -f .$${bn}_h.tmp ; \
	done
distclean-o3t-names:	clean-o3t-names
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

.PHONY: tags
