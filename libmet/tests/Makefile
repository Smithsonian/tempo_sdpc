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
LIBS_F = $(MET_LIBPATHS) -lmet_f \
         $(TELL_LIBPATHS) -ltell_f -ltell \
         $(TIO_LIBPATHS) -ltio_f -ltio \
         $(NETCDF_LIBPATHS) -lnetcdff -lnetcdf $(RPATH)

COMPILE_F_CMD = $(FC) $(FCFLAGS)

# ---------------------------------------------------------

EXECUTABLES = test_nc

all:  $(EXECUTABLES)

check:	$(EXECUTABLES)
	@./test_nc && echo "running test_nc...ok"

GCDA_FILES_DIR = "../src/gnuobjs"

coverage:	$(EXECUTABLES)
	# zero all execution counters
	lcov -z -d $(GCDA_FILES_DIR) -d .
	# create baseline coverage data file
	lcov -c -i -d $(GCDA_FILES_DIR) -d . -o base.info
	# perform test
	$(MAKE) check
	# create test coverage data file
	lcov -c -d $(GCDA_FILES_DIR) -d . -o test_cov.info
	# combine baseline and test coverage data
	lcov -a base.info -a test_cov.info -o total.info
	# remove selected paths from coverage statistics
	lcov -r total.info $(OTS_ROOT)/include/cfortran.h -o filtered.info
	# generate summary web pages
	genhtml filtered.info -o lcov.dir

test_nc:
	$(COMPILE_F_CMD) -o test_nc $(INCS_F) test_nc.f90 $(LIBS_F)

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o *.mod
covclean:
	-$(RM) -f *.gcda *.gcno *.info
	-$(RM) -rf lcov.dir
distclean: clean covclean
	-$(RM) -f $(EXECUTABLES) *.nc

.PHONY: clean covclean distclean check coverage $(EXECUTABLES)
