# -*- make -*-

include ../../Makefile.defines

srcdir  = $(config_dir)/L2_regrid/src
elfdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src

TIO_LIBPATHS    = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
TELL_LIBPATHS   = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
REGRID_LIBPATHS = -Wl,-rpath=$(elfdir) -L$(elfdir)

INCS   = -I$(TIO_INC) -I$(TELL_INC) -I$(NETCDF_INC) -I$(srcdir)
LIBS   = $(REGRID_LIBPATHS) -lregrid \
         $(TIO_LIBPATHS) -ltio \
         $(TELL_LIBPATHS) -ltell \
         $(NETCDF_LIBPATHS) -lnetcdf $(RPATH) -lm

PROJ_DIR = $(HOME)/sys/$(SYS_CARCH)
PROJ_INC = $(PROJ_DIR)/include
PROJ_LIB = $(PROJ_DIR)/lib
PROJ_LIBPATHS = -Wl,-rpath=$(PROJ_LIB) -L$(PROJ_LIB)

EXTRA_INCS = -I$(PROJ_INC)
EXTRA_LIBS = $(PROJ_LIBPATHS) -lproj

COMPILE_CMD = $(CC) $(CCFLAGS)

# ---------------------------------------------------------

EXECUTABLES = test_poly test_pixel mktestdata_l2 test_area

all:  $(EXECUTABLES)

check:	$(EXECUTABLES)
	@./test_poly && echo "tested poly...ok"
	@./test_pixel && echo "tested pixel...ok"
	@./test_run.sh && echo "tested L2_regrid...ok"

check-area:   test_area
	./test_area

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

test_poly:
	$(COMPILE_CMD) -o test_poly $(INCS) test_poly.c $(LIBS) -lm

test_pixel:
	$(COMPILE_CMD) -o test_pixel $(INCS) test_pixel.c $(LIBS) -lm

test_area:
	$(COMPILE_CMD) -o test_area $(INCS) $(EXTRA_INCS) test_area.c $(LIBS) $(EXTRA_LIBS) -lm

mktestdata_l2:
	$(COMPILE_CMD) -o mktestdata_l2 mktestdata_l2.c $(EXTRA_INCS) $(INCS) $(EXTRA_LIBS) $(LIBS) -lm

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o *.mod
	-$(RM) -f tmp/*
	-rmdir tmp
covclean:
	-$(RM) -f *.gcda *.gcno *.info
	-$(RM) -rf lcov.dir
distclean: clean covclean
	-$(RM) -f $(EXECUTABLES) *.nc

.PHONY: clean covclean distclean check coverage $(EXECUTABLES)
