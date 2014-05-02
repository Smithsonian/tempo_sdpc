AC_INIT(L1_trace_gas/src/OMSAO_main.f90)

CONFIG_DIR=`pwd`
AC_SUBST(CONFIG_DIR)

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_CANONICAL_HOST

JH_ANSI_CC
AC_SUBST(CC)
AC_SUBST(CFLAGS)

AC_PROG_FC
AC_SUBST(FC)
AC_SUBST(FCFLAGS)

JH_SDPTK_SETUP

AC_ARG_WITH(hdfeos,
  [  --with-hdfeos[=DIR]      Use HDFEOS libraries installed in DIR],
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
if ! test "x$jh_with_hdfeos" = "xno"; then
   HDFEOS5_LIB="$jh_with_hdfeos5/lib/$OMIUTIL_SYSDIR"
   HDFEOS5_INC="$jh_with_hdfeos5/include"
fi
AC_SUBST(HDFEOS5_LIB)
AC_SUBST(HDFEOS5_INC)

JH_LIST_SOURCES(L1_trace_gas/src)

dnl L1_trace_gas/src/Makefile:L1_trace_gas/src/Makefile.in \

AC_CONFIG_FILES([\
 Makefile:autoconf/Makefile.in \
 Makefile.defines:autoconf/Makefile.defines.in \
])
AC_OUTPUT


