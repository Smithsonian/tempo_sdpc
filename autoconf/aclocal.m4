dnl# -*- mode: sh; mode: fold -*-

AC_DEFUN([JH_ANSI_CC],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_CPP])
AC_ISC_POSIX
])

AC_DEFUN([JH_SDPTK_SETUP],
[
AC_ARG_WITH(otsroot,
  [  --with-otsroot[=DIR]      Use OTS libraries installed in DIR],
  [jh_with_otsroot=$withval],
  [jh_with_otsroot=no])
if test "x$jh_with_otsroot" = "xno"; then
   if test X"$GCC" = X"yes" ; then
      COMPILER_FAMILY=gfortran
   else
      COMPILER_FAMILY=intel
   fi
   OTS_ROOT="/nfs/$HOST/d2/omi/install/$COMPILER_FAMILY"
else
   OTS_ROOT="$jh_with_otsroot"
fi
AC_SUBST(OTS_ROOT)

AC_CHECK_SIZEOF(long)
if test "$ac_cv_sizeof_long" -eq 4 ; then
   OMIUTIL_SYSDIR=linux32
elif test "$ac_cv_sizeof_long" -eq 8 ; then
   OMIUTIL_SYSDIR=linux64
fi

dnl # Let's try to avoid running this stupid ksh setup script:
dnl OMIUTIL="$OTS_ROOT/sdptk/SDPTK5.2.18v1.00"
dnl ksh -c "source $OMIUTIL/TOOLKIT/bin/$OMIUTIL_SYSDIR/pgs-dev-env.ksh; env > .sdptk_env"
dnl hdf_envs=`grep HDF .sdptk_env`
dnl eval `echo $hdf_envs`
dnl pgs_envs=`grep PGS .sdptk_env`
dnl eval `echo $pgs_envs`
dnl # Instead of running the ksh setup script,
dnl # we'll set by hand the few symbols we actually need:
HDFLIB="\$(OTS_ROOT)/lib"
HDFINC="\$(OTS_ROOT)/include"
HDF5LIB="\$(OTS_ROOT)/lib"
HDF5INC="\$(OTS_ROOT)/include"
PGSLIB="\$(OMIUTIL)/TOOLKIT/lib/$OMIUTIL_SYSDIR"
PGSINC="\$(OMIUTIL)/TOOLKIT/include"
HDFEOS_LIB="\$(OMIUTIL)/TOOLKIT/hdfeos/lib/$OMIUTIL_SYSDIR"
HDFEOS_INC="\$(OMIUTIL)/TOOLKIT/hdfeos/include"
HDFEOS5_LIB="\$(OMIUTIL)/TOOLKIT/hdfeos5/lib/$OMIUTIL_SYSDIR"
HDFEOS5_INC="\$(OMIUTIL)/TOOLKIT/hdfeos5/include"
dnl #
AC_SUBST(HDFLIB)
AC_SUBST(HDFINC)
AC_SUBST(HDFEOS_LIB)
AC_SUBST(HDFEOS_INC)
AC_SUBST(HDF5LIB)
AC_SUBST(HDF5INC)
AC_SUBST(HDFEOS5_LIB)
AC_SUBST(HDFEOS5_INC)
AC_SUBST(PGSLIB)
AC_SUBST(PGSINC)
])

AC_DEFUN([JH_LIST_SOURCES],
[
dir="$1"
F90SOURCES=`grep -v '\#' $dir/f90_files.lis|tr '\n' ' '`
CSOURCES=`cat $dir/c_files.lis|tr '\n' ' '`
AC_SUBST(F90SOURCES)
AC_SUBST(CSOURCES)

F90OBJS=""
for file in $F90SOURCES ; do
   bn=`basename $file .f90`
   F90OBJS="\$(objdir)/$bn.o $F90OBJS"
done
COBJS=""
for file in $CSOURCES ; do
   bn=`basename $file .c`
   COBJS="\$(objdir)/$bn.o $COBJS"
done
LINKOBJS="$F90OBJS $COBJS"
AC_SUBST(LINKOBJS)
])
