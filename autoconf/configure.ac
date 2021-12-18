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

AC_ARG_WITH(atlasblas,
  [  --with-atlasblas[=DIR]     Use atlas library installed in DIR],
  [jh_with_atlasblas=$withval],
  [jh_with_atlasblas=no])
if test "x$jh_with_atlasblas" = "xno"; then
   GSL_BLAS="-lgslcblas"
else
   if test x"$jh_with_atlasblas" = xyes ; then
      GSL_BLAS="-L/usr/lib64/atlas -lsatlas"
   else
      GSL_BLAS="-L$jh_with_atlasblas -lsatlas"
   fi
fi
AC_SUBST(GSL_BLAS)

AC_ARG_WITH(archive-home,
  [  --with-archive-home[=DIR]  Put archive directory trees in DIR],
  [jh_with_archive_home=$withval],
  [jh_with_archive_home=no])
if test "x$jh_with_archive_home" = "xno"; then
   ARCHIVE_HOME="/tempo/nas0/sdpc_archive"
else
   ARCHIVE_HOME="$jh_with_archive_home"
fi
AC_SUBST(ARCHIVE_HOME)

AC_ARG_WITH(pipe-home,
  [  --with-pipe-home[=DIR]     Put pipeline directories in DIR],
  [jh_with_pipe_home=$withval],
  [jh_with_pipe_home=no])
if test "x$jh_with_pipe_home" = "xno"; then
   PIPE_HOME="/tempo/nas0/sdpc/proc"
else
   PIPE_HOME="$jh_with_pipe_home"
fi
AC_SUBST(PIPE_HOME)

AC_ARG_WITH(refdata,
  [  --with-refdata[=DIR]       Use refdata directory=DIR],
  [jh_with_refdata=$withval],
  [jh_with_refdata=no])
if test "x$jh_with_refdata" = "xno"; then
   REFDATA_DIR="/tempo/nas0/sdpc_soft/refdata"
else
   REFDATA_DIR="$jh_with_refdata"
fi
AC_SUBST(REFDATA_DIR)

AC_ARG_WITH(ancdata,
  [  --with-ancdata[=DIR]       Use ancillary data directory=DIR],
  [jh_with_ancdata=$withval],
  [jh_with_ancdata=no])
if test "x$jh_with_ancdata" = "xno"; then
   ANCILLARY_DATA_DIR="/tempo/nas0/sdpc/ancillary"
else
   ANCILLARY_DATA_DIR="$jh_with_ancdata"
fi
AC_SUBST(ANCILLARY_DATA_DIR)

AC_ARG_WITH(s6root,
  [  --with-s6root[=DIR]        Use s6 supervisor installed in DIR],
  [jh_with_s6root=$withval],
  [jh_with_s6root=no])
if test "x$jh_with_s6root" = "xno"; then
   S6_ROOT_DIR="/dev/null"
else
   S6_ROOT_DIR="$jh_with_s6root"
fi
AC_SUBST(S6_ROOT_DIR)

AC_ARG_WITH(inrroot,
  [  --with-inrroot[=DIR]       Use INR software installed in DIR],
  [jh_with_inrroot=$withval],
  [jh_with_inrroot=no])
if test "x$jh_with_inrroot" = "xno"; then
   INRSW_ROOT="$OTS_ROOT"
   INR_LIB="\$(OTS_ROOT)/lib"
   INR_INC="\$(OTS_ROOT)/include"
else
   INRSW_ROOT="$jh_with_inrroot"
   INR_LIB="$jh_with_inrroot/lib"
   INR_INC="$jh_with_inrroot/include"
fi
AC_SUBST(INR_LIB)
AC_SUBST(INR_INC)
AC_SUBST(INRSW_ROOT)

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
fi
if test -f "${INROTS_LIB}/libtmglib.a" ; then
  INROTS_EXTRA_LIBS="-ltmglib"
else
  INROTS_EXTRA_LIBS=""
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
   LEAP_SECONDS_LIST="$REFDATA_DIR/zoneinfo/leap-seconds.list"
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

AC_ARG_WITH(sc-timezone,
  [  --with-sc-timezone[=HH]   Set the host spacecraft time zone (integer hours relative to UTC)],
  [jh_with_timezone=$withval],
  [jh_with_timezone=no])
if test "x$jh_with_timezone" != "xno"; then
   AC_DEFINE_UNQUOTED(SPACECRAFT_TIMEZONE,[($jh_with_timezone)])
   echo "Setting spacecraft time zone to $jh_with_timezone"
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
)

AC_CONFIG_FILES([\
 Makefile:autoconf/Makefile.in \
 Makefile.defines:autoconf/Makefile.defines.in \
 pipe/src/Makefile:pipe/src/Makefile.in \
 pipe/slurm/Makefile:pipe/slurm/Makefile.in \
 pipe/etc/sdpc_env.sh:pipe/etc/sdpc_env.sh.in \
 cron/Makefile:cron/Makefile.in \
 cron/crontab:cron/crontab.in \
 cron/logrotate.conf:cron/logrotate.conf.in \
 utils/Makefile:utils/Makefile.in
])
AC_OUTPUT

echo ""
echo "Configuration summary:"
echo "     Install dir: $prefix"
echo "    INR software: $INRSW_ROOT"
echo "  Reference data: $REFDATA_DIR"
echo "  Ancillary data: $ANCILLARY_DATA_DIR"
echo "   Pipeline home: $PIPE_HOME"
echo "    Archive home: $ARCHIVE_HOME"
echo ""
if test "x$jh_with_timezone" != "xno"; then
   echo "   SPACECRAFT_TIMEZONE: $jh_with_timezone"
fi
