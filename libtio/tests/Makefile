# -*- make -*-

include ../../Makefile.defines

THIS_LIB = tio

srcdir  = $(config_dir)/lib$(THIS_LIB)/src
elfdir  = $(srcdir)/objs

TIO_LIB = -Wl,-rpath=$(elfdir) -L$(elfdir) -ltio
TIO_INC = -I$(srcdir)

INCS   = $(TIO_INC) -I$(NETCDF_INC)
LIBS   = $(TIO_LIB) -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB) -lnetcdf -lm
# INCS_F = -I$(elfdir)
# LIBS_F = $(TIO_LIB) -lrefdb_f

COMPILE_CMD = $(CC) $(CFLAGS)
# COMPILE_F_CMD = $(FC) $(FCFLAGS)

# ---------------------------------------------------------

all:  tio_test
check:	tio_test
	@./tio_test 2> tio_test.log && echo ok && /bin/rm delete_me.nc tio_test.log

tio_test:	tio_test.c
	$(COMPILE_CMD) -o tio_test $(INCS) tio_test.c $(LIBS)

# tio_test_f:	tio_test_f.f90
# 	$(COMPILE_F_CMD) -o tio_test_f $(INCS_F) tio_test_f.f90 $(LIBS_F)

# ---------------------------------------------------------

clean:
	-$(RM) -f .##* *.o
distclean: clean
	-$(RM) -f tio_test tio_test_f

.PHONY: clean distclean check
