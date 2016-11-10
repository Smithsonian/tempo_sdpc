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

JH_LIST_SOURCES(L1_trace_gas/src)

dnl AC_CONFIG_HEADER(include/config.h:include/config.hin)

AC_CONFIG_FILES([\
 Makefile:autoconf/Makefile.in \
 Makefile.defines:autoconf/Makefile.defines.in \
])
AC_OUTPUT


