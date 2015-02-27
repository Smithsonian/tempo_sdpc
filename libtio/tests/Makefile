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

INCS_F = -I$(elfdir) -I$(NETCDF_INC)
LIBS_F = $(TIO_LIBPATHS) -ltio -ltio_f \
         $(TELL_LIBPATHS) -ltell_f \
	 $(NETCDF_LIBPATHS) -lnetcdff $(RPATH)

COMPILE_CMD = $(CC) $(CFLAGS)
COMPILE_F_CMD = $(FC) $(FCFLAGS)

# ---------------------------------------------------------

all:  tio_test tio_test_f
check:	tio_test tio_test_f
	@./tio_test 2> tio_test.log && echo "testing C...ok"
	@./tio_test_f 2>> tio_test.log && echo "testing Fortran...ok"
	/bin/rm delete_*.nc tio_test.log

tio_test:
	$(COMPILE_CMD) -o tio_test $(INCS) tio_test.c $(LIBS)

tio_test_f:
	$(COMPILE_F_CMD) -o tio_test_f $(INCS_F) tio_test_f.f90 $(LIBS_F)

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o
distclean: clean
	-$(RM) -f tio_test tio_test_f

.PHONY: clean distclean check tio_test tio_test_f
