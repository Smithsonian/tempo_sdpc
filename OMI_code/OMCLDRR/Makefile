#
SAOPGE = OMCLDRR

ifeq ($(HOST),vex)
   ROOT_DIR = /vex/d2/omi
else
   ROOT_DIR = /home/jedavis/sys/linux-x86_64
#   ROOT_DIR = /home/jhouck/omi/install
#    ROOT_DIR = /home/eosulliv/omi/hdfeos/install
endif

# COMPILER_FAMILY = intel | gfortran
#  OMIUTIL_SYSDIR = linux32 | linux64
#
# COMPILER_FAMILY = intel
COMPILER_FAMILY = gfortran
OMIUTIL_SYSDIR  = linux64
INTEL_ROOT = /opt/intel

OPTROOT = $(ROOT_DIR)/$(COMPILER_FAMILY)
OMIROOT = /home/eosulliv/omi/OMI_code

#----------------------------------------------------------------
# Ideally, you shouldn't need to change anything below this line,
# unless maybe you want to change the compiler options.
#----------------------------------------------------------------

#OMISAOSHARED = $(OMIROOT)/OMI_SAO_Shared_VOCs/src
OMCLDRR = $(OMIROOT)/$(SAOPGE)/src/cloudret
MKF_DIR = $(OMIROOT)/build_OMCLDRR

EXECUTABLE = $(SAOPGE).exe

PGEHOME    = $(OMIROOT)/$(SAOPGE)
PGEOUTDIR  = $(PGEHOME)/out
PGETMPDIR  = $(PGEHOME)/tmp
PGETBLDIR  = $(PGEHOME)/tbl
PGEBINDIR  = $(PGEHOME)/bin
PGESRCDIR  = $(PGEHOME)/src/cloudret/
PGSMSG     = $(PGEHOME)/bin

# Note that some things below depend on the
# environment configuration established here:
#
OMIUTIL = $(OPTROOT)/sdptk/SDPTK5.2.18v1.00
IGNORE := $(shell bash -c "mkdir -p $(PGETMPDIR); source $(OMIUTIL)/TOOLKIT/bin/$(OMIUTIL_SYSDIR)/pgs-dev-env.ksh ; env | sed 's/^/export /' > $(PGETMPDIR)/sdptk_env")
include $(PGETMPDIR)/sdptk_env

PGS_PC_INFO_FILE = $(PGEHOME)/src/$(SAOPGE).pcf
PGEVERSION = `grep PGEVERSION $(PGS_PC_INFO_FILE) | cut -d\| -f3 | tr -d \"`

VPATH = $(PGESRCDIR):$(PGETMPDIR)   #:$(OMISAOSHARED)

MV = /bin/mv
RM = /bin/rm
CP = /bin/cp

ifeq ($(COMPILER_FAMILY), gfortran)
   CC = gcc
   FC = gfortran

   # GFORT_FLAGS = -Wall -Wextra -W -pedantic -fimplicit-none -Warray-bounds -Wno-unused -fcheck-array-temporaries
   GFORT_FLAGS = -Wall -W -Wextra -Warray-bounds -Wno-unused-dummy-argument \
                 -fbounds-check -fmax-errors=0 -fbacktrace \
                 -fcheck-array-temporaries -fconvert=swap
   # GFORT_FLAGS = -Wall -W -Warray-bounds -Wno-unused -fprofile-arcs -ftest-coverage -pg
   # GFORT_OPT_FLAGS = -march=native -funroll-loops -O3
   # GFORT_OPT_FLAGS = -march=native -g -O2
   GFORT_OPT_FLAGS = -g
   # ALL_FFLAGS = $(GFORT_OPT_FLAGS) -fPIC -static $(GFORT_FLAGS)
   ALL_FFLAGS = $(GFORT_OPT_FLAGS) -fPIC $(GFORT_FLAGS)
   MODCMPFLG  =
   MODFLG     = -I#
   MODOUTFLG  = -J#
   ADDLIBD    =
   ADDLIBS    = -Wl,-rpath=$(OPTROOT)/lib -lsz -lm -ldl
   LNFLAGS    = # -lpg
   CFLAGS     += -DLINUX -Df2cFortran -fPIC
endif

ifeq ($(COMPILER_FAMILY), intel)
   FC        = $(INTEL_ROOT)/bin/ifort
   CC        = $(INTEL_ROOT)/bin/icc
   #FFLAGSopt = -xN -O3 -ip -no-prec-div -nowarn
   #FFLAGSopt  = -O3 -static  -warn interfaces
   FFLAGSopt  = -O3 -static
   #FFLAGSdbg  = -g -C -traceback -check bounds -check format -check output_conversion -check uninit \
   #             -warn -fp-stack-check -ftrapuv -debug -fpe0 -gen-interfaces -warn interfaces
   #FFLAGSdbg  = -g -traceback -fpe0 -C -check noarg_temp_created -static
   FFLAGSdbg  = -g -C -traceback -check bounds -check format -check output_conversion -check uninit \
                -warn -fp-stack-check -ftrapuv -debug -fpe0 -fPIC -static
   #FFLAGSprof = -O3 -p -prof-dir ../prof -opt-report-file ../prof/opti.log
   ALL_FFLAGS = $(FFLAGSopt)
   MODCMPFLG  =
   MODFLG     = -I#
   MODOUTFLG  = -module#
   ADDLIBD    =
   ADDLIBS    = -lsz -lm
   LNFLAGS    = -static
   CFLAGS     = -O3 -DLINUX -Df2cFortran
endif

PSPLINELIBS  = -lpspline -lezcdf
NETCDFLIBS   = # -lnetcdf
HE4LIBS      = -lhdfeos -lmfhdf -lGctp -ldf -ljpeg -lz
HE5LIBS      = -lhe5_hdfeos -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lhdf5
PGSLIBS      = -lPGSTK
OMILIBS      =

OPTINCDIR=$(OPTROOT)/include
OPTMODDIR=$(OPTROOT)/mod
OPTLIBDIR=$(OPTROOT)/lib

IFLAGS = -I$(OMIUTIL)/include -I$(OPTINCDIR) \
         $(MODFLG)$(OPTMODDIR) $(MODFLG)$(OMIUTIL)/mod/$(OMIUTIL_SYSDIR) \
         $(MODFLG)$(PGETMPDIR) -I$(PGETMPDIR) 
LFLAGS = -L$(OMIUTIL)/lib/$(OMIUTIL_SYSDIR) $(ADDLIBD) -L$(PGETMPDIR)
LIBS   = $(PSPLINELIBS) $(NETCDFLIBS)                \
         $(OMILIBS) $(PGSLIBS) $(HE5LIBS) $(HE4LIBS) $(ADDLIBS)

NETCDFINC  = $(OMIUTIL)/netcdf/include
NETCDFLIB  = $(OMIUTIL)/netcdf/lib
PSPLINEMOD = $(OPTMODDIR)
PSPLINELIB = $(OPTLIBDIR)

IFLAGS = -I$(OMIUTIL)/include                     \
         -I$(HDFEOS_INC)  $(MODFLG)$(HDFEOS_INC)  \
         -I$(HDF5INC)     $(MODFLG)$(HDF5INC)     \
         -I$(HDFINC)      $(MODFLG)$(HDFINC)      \
         -I$(PGSINC)      $(MODFLG)$(PGSINC)      \
         -I$(HDFEOS5_INC) $(MODFLG)$(HDFEOS5_INC) \
         -I$(PGETMPDIR)   $(MODFLG)$(PGETMPDIR)   \
         -I$(NETCDFINC)   $(MODFLG)$(NETCDFINC)   \
                          $(MODFLG)$(PSPLINEMOD)

IFLAGS = $(MODFLG)$(OMIUTIL)/include \
         $(MODFLG)$(HDFEOS_INC) \
         $(MODFLG)$(HDF5INC) \
         $(MODFLG)$(HDFINC) \
         $(MODFLG)$(PGSINC) \
         $(MODFLG)$(HDFEOS5_INC) \
         $(MODFLG)$(PGETMPDIR) \
         $(MODFLG)$(NETCDFINC) \
	 $(MODFLG)$(PSPLINEMOD) 

LFLAGS = -L$(PGETMPDIR) -L$(OMIUTIL)/lib                              \
         -L$(HDFEOS_LIB) -L$(HDFLIB) -L$(PGSLIB) -L$(HDF5LIB) -L$(HDFEOS5_LIB) \
         -L$(NETCDFLIB) -L$(PSPLINELIB) $(ADDLIBD)

LIBS   = $(PSPLINELIBS) $(NETCDFLIBS)                \
         $(OMILIBS) $(PGSLIBS) $(HE5LIBS) $(HE4LIBS) $(ADDLIBS)

FFLAGS   = $(MODFLG) ./ $(ALL_FFLAGS)
CINCFLGS = -I$(PGSINC) -I$(HDFINC) -I$(HDFEOS_INC) -I$(PGETMPDIR)

#---------------------------------------------------------------
#  Construct paths to most C and Fortran
#  source and object files, excluding the main program OMCLDRR.f90:
#
#  f90_files.lis:  one filename per line, including .f90 suffix
#    c_files.lis:  one filename per line, including .c suffix
#---------------------------------------------------------------
F90_FILE_LIST = $(MKF_DIR)/f90_files.lis
C_FILE_LIST   = $(MKF_DIR)/c_files.lis

F90SOURCES = $(shell cat $(F90_FILE_LIST))
CSOURCES = $(shell cat $(C_FILE_LIST))

F90SRCPATHS = $(addprefix $(OMCLDRR)/, $(F90SOURCES))
CSRCPATHS = $(addprefix $(OMCLDRR)/, $(CSOURCES))

F90OBJECTS = $(F90SOURCES:.f90=.o)
COBJECTS = $(CSOURCES:.c=.o)

LINKOBJS = $(addprefix $(PGETMPDIR)/,$(F90OBJECTS) $(COBJECTS))

#---------------------------------
# Primary 'make' targets for users
#   make [all|clean]
#---------------------------------

RULES_FILE = $(MKF_DIR)/rules.inc
#RULES_FILE = $(PGETMPDIR)/rules.inc

all:	dirs messages $(RULES_FILE) $(EXECUTABLE)

rules:
	$(MKF_DIR)/mkrules.sh f90 $(F90SRCPATHS) > $(RULES_FILE)
	$(MKF_DIR)/mkrules.sh c $(CSRCPATHS) >> $(RULES_FILE)

$(EXECUTABLE):	$(LINKOBJS)
	$(FC) $(PGEHOME)/src/cloudret/$(SAOPGE).f90 -o $(PGETMPDIR)/$(EXECUTABLE) $(LINKOBJS) $(LNFLAGS) $(FFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS)

install:
	$(CP) $(PGETMPDIR)/$(EXECUTABLE) $(PGEOUTDIR)

clean:
	-$(RM) $(PGETMPDIR)/*.o $(PGETMPDIR)/*.mod
	-$(RM) $(RULES_FILE)

# -----------------------------------------------------------------
# The hdf.f90 file contains numerous non-conforming tabs.
# It's unlikely anyone will ever fix that in the HDF distribution,
# so I'll just fix it myself:
# -----------------------------------------------------------------
fix-hdf-tabs:
	sed --in-place=_ 's,\t, ,g' $(OPTROOT)/include/hdf.f90

emptydirs = $(PGETMPDIR) $(PGSMSG) $(PGEBINDIR) $(PGEOUTDIR) $(PGETBLDIR)
SHELL=/bin/bash
export SHELL

dirs:
	@for d in $(emptydirs); do \
	  if [ ! -d $$d ]; then \
	      mkdir -p $$d ; \
	  fi; \
	done

.f.o:
	$(FC) -c $(FFLAGS) $(IFLAGS) $< -o $(PGETMPDIR)/$@
	#$(FC) -c $(FFLAGS) $(IFLAGS) $< -o $@
	#$(MV) -f $@ $(PGETMPDIR)/

.F.o:
	$(FC) -c $(FFLAGS) $(IFLAGS) $< -o $(PGETMPDIR)/$@
	#$(FC) -c $(FFLAGS) $(IFLAGS) $< -o $@
	#$(MV) -f $@ $(PGETMPDIR)/

# ---------------------------------------
# The PGE and TOOLKIT error message files
# ---------------------------------------

#PGEPGSSRC = OMSAO_52500.t OMI_1900.t
PGEPGSSRC = OMCLDRR_52251.t OMI_1900.t
#PGEPGSCMP = PGS_OMCLDRR_52251.f PGS_OMI_1900.f PGS_OMI_1900.h
PGEPGSCMP = PGS_OMCLDRR_52251.f PGS_OMI_1900.f PGS_OMI_1900.h
TKPGSCMP =  PGS_1 PGS_3 PGS_4 PGS_5 PGS_6 PGS_7 PGS_9 PGS_11 \
            PGS_12 PGS_13 PGS_14 PGS_15 PGS_10 PGS_101 

PGSHOME = $(OMIUTIL)/TOOLKIT/
PGSMSG     = $(PGEHOME)/bin

messages: pgemessages tkmessages

pgemessages: $(PGEPGSCMP)

tkmessages:
	 @for tkpgs in $(TKPGSCMP); do                                 \
	   if [ ! -e $(PGSMSG)/$$tkpgs ]; then                        \
	      echo 'Copying' $$tkpgs 'to' $(PGSMSG);                  \
	      touch $(PGSMSG)/$$tkpgs ; chmod u+w $(PGSMSG)/$$tkpgs;  \
	      cp -p $(PGSHOME)/message/$$tkpgs $(PGSMSG)/;            \
	   fi;                                                        \
	 done

# ---------------------------------------------
# Pattern rule for the creation of SMF messages
# ---------------------------------------------
PGS_%.f : %.t
	$(PGSBIN)/smfcompile -r -f77 -f $< -o $(PGETMPDIR)/$@
	#$(PGSBIN)/smfcompile -r -f77 -f $< -o $@
	#$(MV) $(OMCLDRR)/$@ $(PGETMPDIR)/
PGS_%.h : %.t
	$(PGSBIN)/smfcompile -r      -f $< -o $(PGETMPDIR)/$@
	#$(PGSBIN)/smfcompile -r      -f $< -o $@
	#$(MV) $(OMCLDRR)/$@ $(PGETMPDIR)/

#-------------------------------------------------------
# Auto-generate and include object file dependency rules
#-------------------------------------------------------

$(RULES_FILE):
	$(MKF_DIR)/mkrules.sh f90 $(F90SRCPATHS) > $(RULES_FILE)
#	$(MKF_DIR)/mkrules.sh c $(CSRCPATHS) >> $(RULES_FILE)

ODIR = $(PGETMPDIR)
SDIR = $(OMCLDRR)
-include $(RULES_FILE)

.PHONY:	rules.inc
