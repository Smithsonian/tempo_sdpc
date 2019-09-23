# -*- make -*-

include ../../Makefile.defines

THIS_LIB = met

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/$(archobjs)

MET_INC = $(srcdir)
TIO_INC = $(config_dir)/libtio/src/$(archobjs)
#TIO_INC = $(prefix)/include

MET_LIBPATHS    = -Wl,-rpath=$(elfdir) -L$(elfdir)
TELL_LIBPATHS   = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
TIO_LIBPATHS    = -Wl,-rpath=$(config_dir)/libtio/src/$(archobjs) -L$(config_dir)/libtio/src/$(archobjs)
#TIO_LIBPATHS   = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
NETCDF_LIBPATHS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)

INCS_F = -I$(elfdir) -I$(TELL_INC) -I$(TIO_INC)
LIBS_F = $(MET_LIBPATHS) -lmet_f -lmet \
         $(TELL_LIBPATHS) -ltell_f -ltell \
         $(TIO_LIBPATHS) -ltio_f -ltio \
         $(NETCDF_LIBPATHS) -lnetcdff -lnetcdf $(RPATH)

COMPILE_F_CMD = $(FC) $(FCFLAGS)

COMPILE_CMD = $(CC) $(CCFLAGS)

# ---------------------------------------------------------

EXECUTABLES = test_met test_met_list ftest_met

all:  $(EXECUTABLES)

# check:	$(EXECUTABLES)
# 	@./test_met && echo "running test_met...ok"

#GCDA_FILES_DIR = "../src/gnuobjs"

# coverage:	$(EXECUTABLES)
#	# zero all execution counters
#	lcov -z -d $(GCDA_FILES_DIR) -d .
#	# create baseline coverage data file
#	lcov -c -i -d $(GCDA_FILES_DIR) -d . -o base.info
#	# perform test
#	$(MAKE) check
#	# create test coverage data file
#	lcov -c -d $(GCDA_FILES_DIR) -d . -o test_cov.info
#	# combine baseline and test coverage data
#	lcov -a base.info -a test_cov.info -o total.info
#	# remove selected paths from coverage statistics
#	lcov -r total.info $(OTS_ROOT)/include/cfortran.h -o filtered.info
#	# generate summary web pages
#	genhtml filtered.info -o lcov.dir

test_met:	$(ELFDIR_ELF_LIB) test_met.c
	$(COMPILE_CMD) -o test_met test_met.c -I$(MET_INC) $(MET_LIBPATHS) -lmet $(LIBS)

test_met_list:	$(ELFDIR_ELF_LIB) test_met_list.c
	$(COMPILE_CMD) -o test_met_list test_met_list.c -I$(MET_INC) $(MET_LIBPATHS) -lmet $(LIBS)

ftest_met:	$(ELFDIR_ELF_F_LIB) ftest_met.f90
	$(FC) $(FCFLAGS) -o ftest_met ftest_met.f90 -I$(elfdir) -L$(elfdir) $(MET_LIBPATHS) -lmet_f -lmet $(LIBS_F)

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o *.mod
covclean:
	-$(RM) -f *.gcda *.gcno *.info
	-$(RM) -rf lcov.dir
distclean: clean covclean
	-$(RM) -f $(EXECUTABLES)

.PHONY: clean covclean distclean check coverage $(EXECUTABLES)
