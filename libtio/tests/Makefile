# -*- make -*-

include ../../Makefile.defines

THIS_LIB = tio

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/$(archobjs)

TIO_INC = $(srcdir)

TIO_LIBPATHS    = -Wl,-rpath=$(elfdir) -L$(elfdir)
TELL_LIBPATHS   = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)

INCS   = -I$(TIO_INC) -I$(TELL_INC) -I$(NETCDF_INC)
LIBS   = $(TIO_LIBPATHS) -ltio \
         $(TELL_LIBPATHS) -ltell \
         $(NETCDF_LIBPATHS) -lnetcdf $(RPATH) -lm

INCS_F = -I$(elfdir) -I$(NETCDF_INC) -I$(TELL_INC)
LIBS_F = $(TIO_LIBPATHS) -ltio -ltio_f \
         $(TELL_LIBPATHS) -ltell_f -ltell \
	 $(NETCDF_LIBPATHS) -lnetcdff $(RPATH)

COMPILE_CMD = $(CC) $(CCFLAGS)
COMPILE_F_CMD = $(FC) $(FCFLAGS)

# ---------------------------------------------------------

EXECUTABLES = tio_test time_test meta_test epoch_test tio_test_f ftest1 ppc_test history_test

all:  $(EXECUTABLES)

check:	$(EXECUTABLES)
	@./time_test < times.in > times.out && diff times.out times.chk && echo "running time_test...ok"
	@./epoch_test && echo "running epoch_test...ok"
	@./tio_test 2> c.test.log && echo "testing C...ok"
	@./tio_test_f 2> f90.test.log && echo "testing Fortran...ok"
	@./ftest1 2> ftest1.test.log && echo "running ftest1...ok"
	@./ppc_test > ppc_test.out && diff ppc_test.out ppc_test.chk && echo "running ppc_test...ok"
	@./history_test && echo "history update [0 args]...ok"
	@./history_test arg1 arg2 arg3 && echo "history update [3 args]...ok"
	@./meta_test < meta.in > meta.out && diff meta.out meta.chk && echo "running meta_test...ok"
	-/bin/rm delete_*.nc ftest1.nc *.test.log ppc_test.out times.out meta.out history_test.nc

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

tio_test:
	$(COMPILE_CMD) -o tio_test $(INCS) tio_test.c $(LIBS)

meta_test:
	$(COMPILE_CMD) -o meta_test $(INCS) meta_test.c $(LIBS)

time_test:
	$(COMPILE_CMD) -o time_test $(INCS) time_test.c $(LIBS)

epoch_test:
	$(COMPILE_CMD) -o epoch_test $(INCS) epoch_test.c $(LIBS)

ppc_test:
	$(COMPILE_CMD) -o ppc_test $(INCS) -DSELF_TEST ../src/tio_ppc.c $(LIBS)

history_test:
	$(COMPILE_CMD) -o history_test $(INCS) history_test.c $(LIBS)

tio_test_f:
	$(COMPILE_F_CMD) -o tio_test_f $(INCS_F) tio_test_f.f90 $(LIBS_F)

ftest1:
	$(COMPILE_F_CMD) -o ftest1 $(INCS_F) ftest1_module.f90 ftest1.f90 $(LIBS_F)

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o *.mod
covclean:
	-$(RM) -f *.gcda *.gcno *.info
	-$(RM) -rf lcov.dir
distclean: clean covclean
	-$(RM) -f $(EXECUTABLES) *.nc

.PHONY: clean covclean distclean check coverage $(EXECUTABLES)
