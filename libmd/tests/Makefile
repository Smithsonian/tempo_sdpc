# -*- make -*-

include ../../Makefile.defines

THIS_LIB = md

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/$(archobjs)

libtio_srcdir = $(config_dir)/libtio/src
libtio_objdir = $(libtio_srcdir)/$(archobjs)

MD_INC = $(srcdir)
TIO_INC = $(prefix)/include

MD_LIBPATHS    = -Wl,-rpath=$(elfdir) -L$(elfdir)
TELL_LIBPATHS   = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
TIO_LIBPATHS    = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
NETCDF_LIBPATHS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)

INCS_F = -I$(libtio_srcdir) -I$(libtio_objdir) -I$(elfdir) -I$(TELL_INC) -I$(TIO_INC)
LIBS_F = $(MD_LIBPATHS) -lmd_f \
         -L$(libtio_objdir) -Wl,-rpath=$(libtio_objdir) $(TIO_LIBPATHS) -ltio_f -ltio \
         $(TELL_LIBPATHS) -ltell_f -ltell \
         $(NETCDF_LIBPATHS) -lnetcdff -lnetcdf $(RPATH)

COMPILE_F_CMD = $(FC) $(FCFLAGS)

# ---------------------------------------------------------

EXECUTABLES = test_attr

all:  $(EXECUTABLES)

check:	$(EXECUTABLES)
	@./test_attr && ./compare_attr.sh && echo "running test_nc...ok"

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

test_attr:
	$(COMPILE_F_CMD) -o test_attr $(INCS_F) test_attr.f90 $(LIBS_F)

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o *.mod
covclean:
	-$(RM) -f *.gcda *.gcno *.info
	-$(RM) -rf lcov.dir
distclean: clean covclean
	-$(RM) -f $(EXECUTABLES) *.nc

.PHONY: clean covclean distclean check coverage $(EXECUTABLES)
