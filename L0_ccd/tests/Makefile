# -*- make -*-

include ../../Makefile.defines

srcdir  = $(config_dir)/L0_ccd/src
elfdir  = $(srcdir)/$(archobjs)

TIO_LIB = $(config_dir)/libtio/src/$(archobjs)
TIO_INC = $(config_dir)/libtio/src

TIO_LIBPATHS    = -Wl,-rpath=$(TIO_LIB) -L$(TIO_LIB)
TELL_LIBPATHS   = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)
L0CCD_LIBPATHS = -Wl,-rpath=$(elfdir) -L$(elfdir)

INCS   = -I$(TIO_INC) -I$(TELL_INC) -I$(NETCDF_INC) -I$(srcdir)
LIBS   = $(TIO_LIBPATHS) -ltio \
         $(TELL_LIBPATHS) -ltell \
	 $(LIBCONFIG_LIBPATHS) -lconfig \
         $(NETCDF_LIBPATHS) -lnetcdf $(RPATH) -lm

EXTRA_INCS =
EXTRA_LIBS =

COMPILE_CMD = $(CC) $(CCFLAGS)
ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)

# ---------------------------------------------------------

EXECUTABLES = test_image test_pixelqf # test_bpix test_ccd test_dark

all:  $(EXECUTABLES)

TEST_OUTPUT_FILES = img_dup_pixels.arr img_dup_pqf.arr dark_table_test.nc \
                    bpix_test.nc

check:	all
	@./test_image && echo "test_image ok"
#	@./test_ccd && echo "test_ccd ok"
#	@./test_dark sdc && echo "test_dark sdc ok"
#	@./test_dark fptemp && echo "test_dark fptemp ok"
#	@./test_dark exptime && echo "test_dark exptime ok"
	@./test_pixelqf && echo "test_pixelqf ok"
#	@./test_bpix && echo "test_bpix ok"
	-$(RM) -f $(TEST_OUTPUT_FILES)

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

LIBOBJS = util.o $(elfdir)/ccd.o $(elfdir)/image.o $(elfdir)/dark.o \
          $(elfdir)/pixelqf.o $(elfdir)/bpix.o $(elfdir)/util.o

test_image:	util.o
	$(COMPILE_CMD) -o test_image $(INCS) test_image.c $(LIBOBJS) $(LIBS) -lm

test_ccd:	util.o
	$(COMPILE_CMD) -o test_ccd $(INCS) test_ccd.c $(LIBOBJS) $(LIBS) -lm

test_dark:	util.o
	$(COMPILE_CMD) -o test_dark $(INCS) test_dark.c $(LIBOBJS) $(LIBS) -lm

test_pixelqf:	util.o
	$(COMPILE_CMD) -o test_pixelqf $(INCS) test_pixelqf.c $(LIBOBJS) $(LIBS) -lm

test_bpix:
	$(COMPILE_CMD) -o test_bpix $(INCS) test_bpix.c $(LIBOBJS) $(LIBS) -lm

util.o:       util.c
	$(ELF_COMPILE_CMD) $(INCS) util.c

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o $(TEST_OUTPUT_FILES)
covclean:
	-$(RM) -f *.gcda *.gcno *.info
	-$(RM) -rf lcov.dir
distclean: clean covclean
	-$(RM) -f $(EXECUTABLES)

.PHONY: clean covclean distclean check coverage $(EXECUTABLES)
