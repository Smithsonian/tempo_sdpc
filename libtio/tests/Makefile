# -*- make -*-

include ../../Makefile.defines

THIS_LIB = tio

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/$(archobjs)

TIO_INC = -I$(srcdir)

TIO_LIBPATHS    = -Wl,-rpath=$(elfdir) -L$(elfdir)
TELL_LIBPATHS   = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)

INCS   = $(TIO_INC) -I$(NETCDF_INC)
LIBS   = $(TIO_LIBPATHS) -ltio \
         $(TELL_LIBPATHS) -ltell \
         $(NETCDF_LIBPATHS) -lnetcdf $(RPATH) -lm

INCS_F = -I$(elfdir) -I$(NETCDF_INC) -I$(TELL_INC)
LIBS_F = $(TIO_LIBPATHS) -ltio -ltio_f \
         $(TELL_LIBPATHS) -ltell_f \
	 $(NETCDF_LIBPATHS) -lnetcdff $(RPATH)

COMPILE_CMD = $(CC) $(CCFLAGS)
COMPILE_F_CMD = $(FC) $(FCFLAGS)

# ---------------------------------------------------------

EXECUTABLES = tio_test tio_test_f ftest1

all:  $(EXECUTABLES)

check:	$(EXECUTABLES)
	@./tio_test 2> c.test.log && echo "testing C...ok"
#	@./tio_test_f 2> f90.test.log && echo "testing Fortran...ok"
	@./ftest1 2> ftest1.test.log && echo "running ftest1...ok"
	-/bin/rm delete_*.nc ftest1.nc *.test.log

GCDA_FILES_DIR = "../src/gnuobjs"

coverage:	$(EXECUTABLES)
	# zero all execution counters
	lcov -z -d $(GCDA_FILES_DIR)
	# create baseline coverage data file
	lcov -c -i -d $(GCDA_FILES_DIR) -o base.info
	# perform test
	$(MAKE) check
	# create test coverage data file
	lcov -c -d $(GCDA_FILES_DIR) -o test_cov.info
	# combine baseline and test coverage data
	lcov -a base.info -a test_cov.info -o total.info
	# remove selected paths from coverage statistics
	lcov -r total.info $(OTS_ROOT)/include/cfortran.h -o filtered.info
	# generate summary web pages
	genhtml filtered.info -o lcov.dir

tio_test:
	$(COMPILE_CMD) -o tio_test $(INCS) tio_test.c $(LIBS)

tio_test_f:
	$(COMPILE_F_CMD) -o tio_test_f $(INCS_F) tio_test_f.f90 $(LIBS_F)

ftest1:
	$(COMPILE_F_CMD) -o ftest1 $(INCS_F) ftest1.f90 $(LIBS_F)

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o
covclean:
	-$(RM) -f *.gcda *.gcno *.info
	-$(RM) -rf lcov.dir
distclean: clean covclean
	-$(RM) -f $(EXECUTABLES) *.nc

.PHONY: clean covclean distclean check coverage $(EXECUTABLES)
