include ../../../Makefile.defines

srcdir  = $(config_dir)/libtio/src
elfdir  = $(srcdir)/$(archobjs)

TIO_LIBPATHS    = -Wl,-rpath=$(elfdir) -L$(elfdir)
NETCDF_LIBPATHS = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)

INCS_F = -I$(elfdir) -I$(NETCDF_INC)
LIBS_F = $(TIO_LIBPATHS) -ltio -ltio_f \
	 $(NETCDF_LIBPATHS) -lnetcdff

COMPILE_F_CMD = $(FC) $(FCFLAGS)

example1:	example1.f90
	$(COMPILE_F_CMD) -o example1 example1.f90 $(INCS_F) $(LIBS_F)

clean:
	$(RM) ./example1
