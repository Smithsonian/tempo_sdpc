# -*- make -*-

include ../../Makefile.defines

srcdir  = $(config_dir)/wavecal/src
objdir  = $(srcdir)/$(archobjs)

TIO_INC = $(config_dir)/libtio/src
TIO_LIB = $(config_dir)/libtio/src/$(archobjs)

TELL_LIBPATHS    = -Wl,-rpath=$(TELL_LIB) -L$(TELL_LIB)
NETCDF_LIBPATHS  = -Wl,-rpath=$(NETCDF_LIB) -L$(NETCDF_LIB)
TIO_LIBPATHS     = -L$(TIO_LIB)
INSTALL_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib
LIBCONFIG_LIBPATHS = -Wl,-rpath=$(LIBCONFIG_LIB) -L$(LIBCONFIG_LIB)

IOCSDPC_INC = $(prefix)/include
IOCSDPC_LIBPATHS = -Wl,-rpath=$(prefix)/lib -L$(prefix)/lib -lioclib -lcsv

EXTRA_INCS =#
EXTRA_LIBS = $(GSL_LIBPATHS) #

INCS = -I$(TIO_INC) -I$(IOCSDPC_INC) -I$(TELL_INC) -I$(LIBCONFIG_INC) -I$(NETCDF_INC) \
        -I$(config_dir)/include
LIBS = $(TIO_LIBPATHS) -ltio -lmpfit \
       $(LIBCONFIG_LIBPATHS) -lconfig \
       $(NETCDF_LIBPATHS) -lnetcdf $(IOCSDPC_LIBPATHS) \
       $(TELL_LIBPATHS) -ltell \
       $(RPATH) $(INSTALL_LIBPATHS)

ELF_COMPILE_CMD = $(ELF_CC) -c $(ELF_CFLAGS)
ELF_LINK_CMD = $(ELF_LINK)

COMPILE_CMD = $(CC) $(CCFLAGS)

# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

all:	 $(objdir) $(objdir)/wavecal_driver $(objdir)/wavecal_merge $(objdir)/wavecal_dump

TEST_EXECS = test_asg test_convolve test_table # test_interp

tests:	$(objdir) $(TEST_EXECS)

$(objdir):
	-mkdir -p $(objdir)

tags:
	ctags $(WAVECAL_DRIVER_SRCFILES)

.PHONY: clean distclean all install check tags

# --------------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------------

bindir = $(prefix)/bin
incdir = $(prefix)/include
libdir = $(prefix)/lib

install-executables:
	mkdir -p $(bindir)
	@echo installing executables in $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/wavecal_driver $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/wavecal_merge $(bindir)/
	$(INSTALL_PROGRAM) $(objdir)/wavecal_dump $(bindir)/

install:	install-executables

# --------------------------------------------------------------------------
# Housekeeping
# --------------------------------------------------------------------------

clean:
	-$(RM) -f .##* $(objdir)/* $(TEST_EXECS)
distclean: clean
	-$(RM) -rf $(objdir)

# --------------------------------------------------------------------------
# Object Rules
# --------------------------------------------------------------------------

TEST_INTERP_OFILES = test_interp.o lsq.o interp.o \
              interp_cheb.o interp_poly.o interp_cspline.o
TEST_INTERP_OBJS = $(patsubst %.o,$(objdir)/%.o, $(TEST_INTERP_OFILES))

test_interp:	$(TEST_INTERP_OBJS)
	$(COMPILE_CMD) -o test_interp $(TEST_INTERP_OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/test_interp.o:	$(srcdir)/test_interp.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/test_interp.c

$(objdir)/wavecal_merge:	$(objdir)/wavecal_merge.o
	$(COMPILE_CMD) -o $(objdir)/wavecal_merge $(objdir)/wavecal_merge.o \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/wavecal_merge.o:	$(srcdir)/wavecal_merge.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal_merge.c

$(objdir)/wavecal_dump:	$(objdir)/wavecal_dump.o
	$(COMPILE_CMD) -o $(objdir)/wavecal_dump $(objdir)/wavecal_dump.o \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/wavecal_dump.o:	$(srcdir)/wavecal_dump.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal_dump.c

WAVECAL_DRIVER_SRCFILES = wavecal_driver.c wavecal_adj.c util.c wavecal.c interp.c shapefun.c \
     lsq.c interp_cheb.c interp_poly.c interp_cspline.c \
     slit_function_asg.c slit_function_apply.c slit_function_table.c
WAVECAL_DRIVER_OFILES = $(patsubst %.c,%.o, $(WAVECAL_DRIVER_SRCFILES))
WAVECAL_DRIVER_OBJS = $(patsubst %.o,$(objdir)/%.o, $(WAVECAL_DRIVER_OFILES))

$(objdir)/wavecal_driver:	$(WAVECAL_DRIVER_OBJS)
	$(COMPILE_CMD) -o $(objdir)/wavecal_driver $(WAVECAL_DRIVER_OBJS) \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

$(objdir)/wavecal_driver.o:	$(srcdir)/wavecal_driver.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal_driver.c

$(objdir)/util.o:	$(srcdir)/util.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/util.c

$(objdir)/wavecal.o:	$(srcdir)/wavecal.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal.c

$(objdir)/wavecal_adj.o:	$(srcdir)/wavecal_adj.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/wavecal_adj.c

$(objdir)/interp.o:	$(srcdir)/interp.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp.c

$(objdir)/shapefun.o:	$(srcdir)/shapefun.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/shapefun.c

$(objdir)/lsq.o:	$(srcdir)/lsq.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/lsq.c

$(objdir)/interp_cspline.o:	$(srcdir)/interp_cspline.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp_cspline.c

$(objdir)/interp_poly.o:	$(srcdir)/interp_poly.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp_poly.c

$(objdir)/interp_cheb.o:	$(srcdir)/interp_cheb.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/interp_cheb.c

$(objdir)/slit_function_asg.o:	$(srcdir)/slit_function_asg.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/slit_function_asg.c

$(objdir)/slit_function_apply.o:	$(srcdir)/slit_function_apply.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/slit_function_apply.c

$(objdir)/slit_function_table.o:	$(srcdir)/slit_function_table.c
	cd $(objdir) && $(ELF_COMPILE_CMD) $(INCS) $(srcdir)/slit_function_table.c

test_asg:	$(srcdir)/slit_function_asg.c
	$(COMPILE_CMD) -o test_asg -DUNIT_TEST $(srcdir)/slit_function_asg.c $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

test_convolve:	$(srcdir)/slit_function_apply.c $(objdir)/slit_function_asg.o
	$(COMPILE_CMD) -o test_convolve -DUNIT_TEST $(srcdir)/slit_function_apply.c \
	    $(objdir)/slit_function_asg.o \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm

test_table:	$(srcdir)/slit_function_table.c $(objdir)/util.o
	$(COMPILE_CMD) -o test_table -DUNIT_TEST $(srcdir)/slit_function_table.c \
	    $(objdir)/util.o \
	    $(INCS) $(EXTRA_INCS) $(LIBS) $(EXTRA_LIBS) -lm
