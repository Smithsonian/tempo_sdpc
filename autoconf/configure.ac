AC_INIT(L1_trace_gas/src/OMSAO_main.f90)

CONFIG_DIR=`pwd`
AC_SUBST(CONFIG_DIR)

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_CANONICAL_HOST

dnl JD_ANSI_CC
dnl JD_ELF_COMPILER

dnl AC_SUBST(CC)
dnl AC_SUBST(CFLAGS)

dnl AC_PROG_FC
dnl AC_SUBST(FC)
dnl AC_SUBST(FCFLAGS)

dnl JD_WITH_LIBRARY(sqlite3)

CFLAGS_COV=""
FCFLAGS_COV=""
LDFLAGS_COV=""
AC_ARG_ENABLE([gcov],
  AS_HELP_STRING([--enable-gcov], [Do not use for production code!]))
AS_IF([test "x$enable_gcov" = "xyes"], [
  CFLAGS_COV="--coverage"
  FCFLAGS_COV="--coverage"
  LDFLAGS_COV="--coverage"
])
AC_SUBST(CFLAGS_COV)
AC_SUBST(FCFLAGS_COV)
AC_SUBST(LDFLAGS_COV)

JH_SDPTK_SETUP

AC_ARG_WITH(hdfeos,
  [  --with-hdfeos[=DIR]       Use HDFEOS libraries installed in DIR],
  [jh_with_hdfeos=$withval],
  [jh_with_hdfeos=no])
if ! test "x$jh_with_hdfeos" = "xno"; then
   HDFEOS_LIB="$jh_with_hdfeos/lib/$OMIUTIL_SYSDIR"
   HDFEOS_INC="$jh_with_hdfeos/include"
fi
AC_SUBST(HDFEOS_LIB)
AC_SUBST(HDFEOS_INC)

AC_ARG_WITH(hdfeos5,
  [  --with-hdfeos5[=DIR]      Use HDFEOS5 libraries installed in DIR],
  [jh_with_hdfeos5=$withval],
  [jh_with_hdfeos5=no])
if ! test "x$jh_with_hdfeos5" = "xno"; then
   HDFEOS5_LIB="$jh_with_hdfeos5/lib/$OMIUTIL_SYSDIR"
   HDFEOS5_INC="$jh_with_hdfeos5/include"
fi
AC_SUBST(HDFEOS5_LIB)
AC_SUBST(HDFEOS5_INC)

AC_ARG_WITH(netcdf,
  [  --with-netcdf[=DIR]       Use netCDF library installed in DIR],
  [jh_with_netcdf=$withval],
  [jh_with_netcdf=no])
if test "x$jh_with_netcdf" = "xno"; then
   NETCDF_LIB="\$(OTS_ROOT)/lib"
   NETCDF_INC="\$(OTS_ROOT)/include"
else
   NETCDF_LIB="$jh_with_netcdf/lib"
   NETCDF_INC="$jh_with_netcdf/include"
fi
AC_SUBST(NETCDF_LIB)
AC_SUBST(NETCDF_INC)

AC_ARG_WITH(proj,
  [  --with-proj[=DIR]         Use proj4 library installed in DIR],
  [jh_with_proj=$withval],
  [jh_with_proj=no])
if test "x$jh_with_proj" = "xno"; then
   PROJ_LIB="\$(OTS_ROOT)/lib"
   PROJ_INC="\$(OTS_ROOT)/include"
else
   PROJ_LIB="$jh_with_proj/lib"
   PROJ_INC="$jh_with_proj/include"
fi
AC_SUBST(PROJ_LIB)
AC_SUBST(PROJ_INC)

AC_ARG_WITH(slnovas,
  [  --with-slnovas[=DIR]      Use slnovas library installed in DIR],
  [jh_with_slnovas=$withval],
  [jh_with_slnovas=no])
if test "x$jh_with_slnovas" = "xno"; then
   SLNOVAS_LIB="\$(OTS_ROOT)/lib"
   SLNOVAS_INC="\$(OTS_ROOT)/include"
else
   SLNOVAS_LIB="$jh_with_proj/lib"
   SLNOVAS_INC="$jh_with_proj/include"
fi
AC_SUBST(SLNOVAS_LIB)
AC_SUBST(SLNOVAS_INC)

AC_ARG_WITH(libconfig,
  [  --with-libconfig[=DIR]    Use libconfig library installed in DIR],
  [jh_with_libconfig=$withval],
  [jh_with_libconfig=no])
if test "x$jh_with_libconfig" = "xno"; then
   LIBCONFIG_LIB="\$(OTS_ROOT)/lib"
   LIBCONFIG_INC="\$(OTS_ROOT)/include"
else
   LIBCONFIG_LIB="$jh_with_libconfig/lib"
   LIBCONFIG_INC="$jh_with_libconfig/include"
fi
AC_SUBST(LIBCONFIG_LIB)
AC_SUBST(LIBCONFIG_INC)

AC_ARG_WITH(eccodes,
  [  --with-eccodes[=DIR]      Use eccodes library installed in DIR],
  [jh_with_eccodes=$withval],
  [jh_with_eccodes=no])
if test "x$jh_with_eccodes" = "xno"; then
   ECCODES_LIB="\$(OTS_ROOT)/lib"
   ECCODES_INC="\$(OTS_ROOT)/include"
else
   ECCODES_LIB="$jh_with_eccodes/lib"
   ECCODES_INC="$jh_with_eccodes/include"
fi
AC_SUBST(ECCODES_LIB)
AC_SUBST(ECCODES_INC)

AC_ARG_WITH(inrlib,
  [  --with-inrlib[=DIR]       Use INR library installed in DIR],
  [jh_with_inrlib=$withval],
  [jh_with_inrlib=no])
if test "x$jh_with_inrlib" = "xno"; then
   INR_LIB="\$(OTS_ROOT)/lib"
   INR_INC="\$(OTS_ROOT)/include"
else
   INR_LIB="$jh_with_inrlib/lib"
   INR_INC="$jh_with_inrlib/include"
fi
AC_SUBST(INR_LIB)
AC_SUBST(INR_INC)

AC_ARG_WITH(inrots,
  [  --with-inrots[=DIR]       Use INR OTS software installed in DIR],
  [jh_with_inrots=$withval],
  [jh_with_inrots=no])
if test "x$jh_with_inrots" = "xno"; then
   INROTS_LIB="\$(OTS_ROOT)/lib"
   INROTS_INC="\$(OTS_ROOT)/include"
else
   INROTS_LIB="$jh_with_inrots/lib"
   INROTS_INC="$jh_with_inrots/include"
   if test -f "${INROTS_LIB}/libtmglib.a" ; then
     INROTS_EXTRA_LIBS="-ltmglib"
   else
     INROTS_EXTRA_LIBS=""
   fi
fi
AC_SUBST(INROTS_LIB)
AC_SUBST(INROTS_INC)
AC_SUBST(INROTS_EXTRA_LIBS)

AC_ARG_WITH(tell,
  [  --with-tell[=DIR]         Use tell library installed in DIR (default:$prefix)],
  [jh_with_tell=$withval],
  [jh_with_tell=no])
if test "x$jh_with_tell" = "xno"; then
   TELL_LIB="\$(OTS_ROOT)/lib"
   TELL_INC="\$(OTS_ROOT)/include"
else
   TELL_LIB="$jh_with_tell/lib"
   TELL_INC="$jh_with_tell/include"
fi
AC_SUBST(TELL_LIB)
AC_SUBST(TELL_INC)

AC_ARG_WITH(tiff,
  [  --with-tiff[=DIR]         Use tiff library installed in DIR (default:$prefix)],
  [jh_with_tiff=$withval],
  [jh_with_tiff=no])
if test "x$jh_with_tiff" = "xno"; then
   TIFF_LIB="\$(OTS_ROOT)/lib"
   TIFF_INC="\$(OTS_ROOT)/include"
else
   TIFF_LIB="$jh_with_tiff/lib"
   TIFF_INC="$jh_with_tiff/include"
fi
AC_SUBST(TIFF_LIB)
AC_SUBST(TIFF_INC)

AC_ARG_WITH(geotiff,
  [  --with-geotiff[=DIR]      Use geotiff library installed in DIR (default:$prefix)],
  [jh_with_geotiff=$withval],
  [jh_with_geotiff=no])
if test "x$jh_with_geotiff" = "xno"; then
   GEOTIFF_LIB="\$(OTS_ROOT)/lib"
   GEOTIFF_INC="\$(OTS_ROOT)/include"
else
   GEOTIFF_LIB="$jh_with_geotiff/lib"
   GEOTIFF_INC="$jh_with_geotiff/include"
fi
AC_SUBST(GEOTIFF_LIB)
AC_SUBST(GEOTIFF_INC)

AC_ARG_WITH(leapsec,
  [  --with-leapsec[=FILE]     Use this leap second file],
  [jh_with_leapsec=$withval],
  [jh_with_leapsec=no])
if test "x$jh_with_leapsec" = "xno"; then
   LEAP_SECONDS_LIST="/usr/share/zoneinfo/leap-seconds.list"
else
   LEAP_SECONDS_LIST="$jh_with_leapsec"
fi
if ! test -f "$LEAP_SECONDS_LIST" ; then
   AC_MSG_ERROR(Leap seconds file not found: $LEAP_SECONDS_LIST)
elif ! test -r "$LEAP_SECONDS_LIST" ; then
   AC_MSG_ERROR(Cannot read leap seconds file: $LEAP_SECONDS_LIST)
fi
AC_SUBST(LEAP_SECONDS_LIST)

AC_ARG_ENABLE(imgcomp,
  [  --disable-imgcomp       Don't compress image data on output],
  [ jh_imgcomp_enable=$enableval ], [jh_imgcomp_enable=yes])
if test X"$jh_imgcomp_enable" = X"yes" ; then
  AC_DEFINE(ENABLE_IMAGE_COMPRESSION)
  echo 'Enable image compression on output'
fi

JD_WITH_LIBRARY(slang)
JD_SLANG_MODULE_INSTALL_DIR
INSTALL_MODULE="\$(INSTALL_DATA)"
AC_SUBST(INSTALL_MODULE)

dnl Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

AC_CHECK_SIZEOF(short, 2)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(double, 8)

JD_SET_RPATH($libdir)

JH_LIST_SOURCES(L1_trace_gas/src)

AC_CONFIG_HEADER(\
 libtio/src/config.h:libtio/src/config.hin\
 L0_format/src/config.h:L0_format/src/config.hin\
 pipe/src/config.h:pipe/src/config.hin\
)

AC_CONFIG_FILES([\
 Makefile:autoconf/Makefile.in \
 Makefile.defines:autoconf/Makefile.defines.in \
 pipe/src/Makefile:pipe/src/Makefile.in \
 pipe/slurm/Makefile:pipe/slurm/Makefile.in \
])
AC_OUTPUT


