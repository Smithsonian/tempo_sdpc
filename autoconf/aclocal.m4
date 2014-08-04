dnl# -*- mode: sh; mode: fold -*-

AC_DEFUN([JD_ANSI_CC], dnl#{{{
[
AC_AIX
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_CPP])
AC_REQUIRE([AC_PROG_GCC_TRADITIONAL])
AC_ISC_POSIX

dnl #This stuff came from Yorick config script
dnl
dnl # HPUX needs special stuff
dnl
AC_EGREP_CPP(yes,
[#ifdef hpux
  yes
#endif
], [
AC_DEFINE(_HPUX_SOURCE,1,[Special define needed for HPUX])
if test "$CC" = cc; then CC="cc -Ae"; fi
])dnl
dnl
dnl #Be sure we've found compiler that understands prototypes
dnl
AC_MSG_CHECKING(C compiler that understands ANSI prototypes)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[ ]], [[
 extern int silly (int);]])],[
 AC_MSG_RESULT($CC looks ok.  Good.)],[
 AC_MSG_RESULT($CC is not a good enough compiler)
 AC_MSG_ERROR(Set env variable CC to your ANSI compiler and rerun configure.)
 ])dnl
])dnl

dnl#}}}

AC_DEFUN([JD_ELF_COMPILER], dnl#{{{
[
dnl #-------------------------------------------------------------------------
dnl # Check for dynamic linker
dnl #-------------------------------------------------------------------------
DYNAMIC_LINK_LIB=""

dnl# AH_TEMPLATE([HAVE_DLOPEN],1,[Define if you have dlopen])

AC_CHECK_HEADER(dlfcn.h,[
  AC_DEFINE(HAVE_DLFCN_H,1,[Define if you have the dlfcn.h header])
  AC_CHECK_LIB(dl,dlopen,[
    DYNAMIC_LINK_LIB="-ldl"
    AC_DEFINE(HAVE_DLOPEN,1,[Define if you have dlopen])
   ],[
    AC_CHECK_FUNC(dlopen,AC_DEFINE(HAVE_DLOPEN,[Define if you have dlopen]))
    if test "$ac_cv_func_dlopen" != yes
    then
      AC_MSG_WARN(cannot perform dynamic linking)
    fi
   ])])
AC_SUBST(DYNAMIC_LINK_LIB)

if test "$GCC" = yes
then
  if test X"$CFLAGS" = X
  then
     CFLAGS="-O2"
  fi
fi

dnl #Some defaults
ELFLIB="lib\$(THIS_LIB).so"
ELFLIB_MAJOR="\$(ELFLIB).\$(ELF_MAJOR_VERSION)"
ELFLIB_MAJOR_MINOR="\$(ELFLIB_MAJOR).\$(ELF_MINOR_VERSION)"
ELFLIB_MAJOR_MINOR_MICRO="\$(ELFLIB_MAJOR_MINOR).\$(ELF_MICRO_VERSION)"
ELFLIB_F="lib\$(THIS_LIB)_f.so"
ELFLIB_F_MAJOR="\$(ELFLIB_F).\$(ELF_MAJOR_VERSION)"

dnl# This specifies the target to use in the makefile to install the shared library
INSTALL_ELFLIB_TARGET="install-elf-and-links"
ELFLIB_BUILD_NAME="\$(ELFLIB_MAJOR_MINOR_MICRO)"
INSTALL_MODULE="\$(INSTALL_DATA)"
M_LIB="-lm"

case "$host_os" in
  *linux*|*gnu*|k*bsd*-gnu )
    DYNAMIC_LINK_FLAGS="-Wl,-export-dynamic"
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -fPIC"
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-O1 -Wl,-soname,\$(ELFLIB_MAJOR)"
    ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-O1 -Wl,-soname,\$(ELFLIB_F_MAJOR)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
    CC_SHARED_FLAGS="-shared -fPIC"
    CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
    ;;
  *solaris* )
    if test "$GCC" = yes
    then
      DYNAMIC_LINK_FLAGS=""
      ELF_CC="\$(CC)"
      ELF_CFLAGS="\$(CFLAGS) -fPIC"
      ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-ztext -Wl,-h,\$(ELFLIB_MAJOR)"
      ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-ztext -Wl,-h,\$(ELFLIB_F_MAJOR)"
      ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
      CC_SHARED_FLAGS="-G -fPIC"
      CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
    else
      DYNAMIC_LINK_FLAGS=""
      ELF_CC="\$(CC)"
      ELF_CFLAGS="\$(CFLAGS) -K PIC"
      ELF_LINK="\$(CC) \$(LDFLAGS) -G -h\$(ELFLIB_MAJOR)"
      ELF_F_LINK="\$(CC) \$(LDFLAGS) -G -h\$(ELFLIB_F_MAJOR)"
      ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
      CC_SHARED_FLAGS="-G -K PIC"
      CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
    fi
    ;;
   # osr5 or unixware7 with current or late autoconf
  *sco3.2v5* | *unixware-5* | *sco-sysv5uw7*)
     if test "$GCC" = yes
     then
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS) -fPIC"
       ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-h,\$(ELFLIB_MAJOR)"
       ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-h,\$(ELFLIB_F_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED_FLAGS="-G -fPIC"
       CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
     else
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS) -K pic"
       # ELF_LINK="ld -G -z text -h#"
       ELF_LINK="\$(CC) \$(LDFLAGS) -G -z text -h\$(ELFLIB_MAJOR)"
       ELF_F_LINK="\$(CC) \$(LDFLAGS) -G -z text -h\$(ELFLIB_F_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED_FLAGS="-G -K pic"
       CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
     fi
     ;;
  *irix6.5* )
     echo "Note: ELF compiler for host_os=$host_os may not be correct"
     echo "double-check: 'mode_t', 'pid_t' may be wrong!"
     if test "$GCC" = yes
     then
       # not tested
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS) -fPIC"
       ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-h,\$(ELFLIB_MAJOR)"
       ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-h,\$(ELFLIB_F_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED_FLAGS="-shared -fPIC"
       CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
     else
       DYNAMIC_LINK_FLAGS=""
       ELF_CC="\$(CC)"
       ELF_CFLAGS="\$(CFLAGS)"     # default anyhow
       ELF_LINK="\$(CC) \$(LDFLAGS) -shared -o \$(ELFLIB_MAJOR)"
       ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -o \$(ELFLIB_F_MAJOR)"
       ELF_DEP_LIBS=
       CC_SHARED_FLAGS="-shared"
       CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
     fi
     ;;
  *darwin* )
     DYNAMIC_LINK_FLAGS=""
     ELF_CC="\$(CC)"
     ELF_CFLAGS="\$(CFLAGS) -fno-common"
     ELF_LINK="\$(CC) \$(LDFLAGS) -dynamiclib -install_name \$(install_lib_dir)/\$(ELFLIB_MAJOR) -compatibility_version \$(ELF_MAJOR_VERSION) -current_version \$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION)"
     ELF_F_LINK="\$(CC) \$(LDFLAGS) -dynamiclib -install_name \$(install_lib_dir)/\$(ELFLIB_F_MAJOR) -compatibility_version \$(ELF_MAJOR_VERSION) -current_version \$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION)"
     ELF_DEP_LIBS="\$(LDFLAGS) \$(DL_LIB)"
     CC_SHARED_FLAGS="-bundle -flat_namespace -undefined suppress -fno-common"
     CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
     ELFLIB="lib\$(THIS_LIB).dylib"
     ELFLIB_MAJOR="lib\$(THIS_LIB).\$(ELF_MAJOR_VERSION).dylib"
     ELFLIB_MAJOR_MINOR="lib\$(THIS_LIB).\$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION).dylib"
     ELFLIB_MAJOR_MINOR_MICRO="lib\$(THIS_LIB).\$(ELF_MAJOR_VERSION).\$(ELF_MINOR_VERSION).\$(ELF_MICRO_VERSION).dylib"
     ;;
  *freebsd* )
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -fPIC"
    #if test "X$PORTOBJFORMAT" = "Xelf" ; then
    #  ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-soname,\$(ELFLIB_MAJOR)"
    #else
    #  ELF_LINK="ld -Bshareable -x"
    #fi
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-soname,\$(ELFLIB_MAJOR)"
    ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-soname,\$(ELFLIB_F_MAJOR)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm"
    CC_SHARED_FLAGS="-shared -fPIC"
    CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
    ;;
  *cygwin* )
    DYNAMIC_LINK_FLAGS=""
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -DBUILD_DLL=1"
    #ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-O1 -Wl,-soname,\$(ELFLIB_MAJOR) -Wl,-export-all-symbols -Wl,-enable-auto-import"
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-O1 -Wl,-soname,\$(ELFLIB_MAJOR)"
    ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-O1 -Wl,-soname,\$(ELFLIB_F_MAJOR)"
    ELF_DEP_LIBS="\$(DL_LIB) -lm"
    CC_SHARED_FLAGS="-shared"
    CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
    dnl# CYGWIN prohibits undefined symbols when linking shared libs
    INSTALL_MODULE="\$(INSTALL)"
    INSTALL_ELFLIB_TARGET="install-elf-cygwin"
    ELFLIB="lib\$(THIS_LIB).dll"
    ELFLIB_MAJOR="cyg\$(THIS_LIB)-\$(ELF_MAJOR_VERSION).dll"
    ELFLIB_MAJOR_MINOR="cyg\$(THIS_LIB)-\$(ELF_MAJOR_VERSION)_\$(ELF_MINOR_VERSION).dll"
    ELFLIB_MAJOR_MINOR_MICRO="cyg\$(THIS_LIB)-\$(ELF_MAJOR_VERSION)_\$(ELF_MINOR_VERSION)_\$(ELF_MICRO_VERSION).dll"
    ELFLIB_BUILD_NAME="\$(ELFLIB_MAJOR)"
    ;;
  *haiku* )
    M_LIB=""
    DYNAMIC_LINK_FLAGS="-Wl,-export-dynamic"
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -fPIC"
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-O1 -Wl,-soname,\$(ELFLIB_MAJOR)"
    ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared -Wl,-O1 -Wl,-soname,\$(ELFLIB_F_MAJOR)"
    ELF_DEP_LIBS="\$(DL_LIB)"
    CC_SHARED_FLAGS="-shared -fPIC"
    CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
    ;;
  * )
    echo "Note: ELF compiler for host_os=$host_os may be wrong"
    ELF_CC="\$(CC)"
    ELF_CFLAGS="\$(CFLAGS) -fPIC"
    ELF_LINK="\$(CC) \$(LDFLAGS) -shared"
    ELF_F_LINK="\$(CC) \$(LDFLAGS) -shared"
    ELF_DEP_LIBS="\$(DL_LIB) -lm -lc"
    CC_SHARED_FLAGS="-shared -fPIC"
    CC_SHARED="\$(CC) $CC_SHARED_FLAGS \$(CFLAGS)"
esac

AC_SUBST(ELF_CC)
AC_SUBST(ELF_CFLAGS)
AC_SUBST(ELF_LINK)
AC_SUBST(ELF_F_LINK)
AC_SUBST(ELF_DEP_LIBS)
AC_SUBST(DYNAMIC_LINK_FLAGS)
AC_SUBST(CC_SHARED_FLAGS)
AC_SUBST(CC_SHARED)
AC_SUBST(ELFLIB)
AC_SUBST(ELFLIB_MAJOR)
AC_SUBST(ELFLIB_MAJOR_MINOR)
AC_SUBST(ELFLIB_MAJOR_MINOR_MICRO)
AC_SUBST(INSTALL_MODULE)
AC_SUBST(INSTALL_ELFLIB_TARGET)
AC_SUBST(ELFLIB_BUILD_NAME)
AC_SUBST(M_LIB)
])

dnl#}}}

RPATH=""
AC_DEFUN([JD_INIT_RPATH], dnl#{{{
[
dnl# determine whether or not -R or -rpath can be used
case "$host_os" in
  *linux*|*solaris* )
    if test "X$GCC" = Xyes
    then
      if test "X$ac_R_nospace" = "Xno"
      then
        RPATH="-Wl,-R,"
      else
        RPATH="-Wl,-R"
      fi
    else
      if test "X$ac_R_nospace" = "Xno"
      then
        RPATH="-R "
      else
	RPATH="-R"
      fi
    fi
  ;;
  *osf*|*openbsd*)
    if test "X$GCC" = Xyes
    then
      RPATH="-Wl,-rpath,"
    else
      RPATH="-rpath "
    fi
  ;;
  *netbsd*)
    if test "X$GCC" = Xyes
    then
      RPATH="-Wl,-R"
    fi
  ;;
esac
])

dnl#}}}

AC_DEFUN([JD_SET_RPATH], dnl#{{{
[
if test "X$1" != "X"
then
  if test "X$RPATH" = "X"
  then
    JD_INIT_RPATH
    if test "X$RPATH" != "X"
    then
      RPATH="$RPATH$1"
    fi
  else
    _already_there=0
    for X in `echo $RPATH | sed 's/:/ /g'`
    do
      if test "$X" = "$1"
      then
        _already_there=1
	break
      fi
    done
    if test $_already_there = 0
    then
      RPATH="$RPATH:$1"
    fi
  fi
fi
])
AC_SUBST(RPATH)dnl

dnl#}}}

AC_DEFUN([JD_UPPERCASE], dnl#{{{
[
changequote(<<, >>)dnl
define(<<$2>>, translit($1, [a-z], [A-Z]))dnl
changequote([, ])dnl
])
#}}}

dnl# This function expand the "prefix variables.  For example, it will expand
dnl# values such as ${exec_prefix}/foo when ${exec_prefix} itself has a
dnl# of ${prefix}.  This function produces the shell variables:
dnl# jd_prefix_libdir, jd_prefix_incdir
AC_DEFUN([JD_EXPAND_PREFIX], dnl#{{{
[
  if test "X$jd_prefix" = "X"
  then
    jd_prefix=$ac_default_prefix
    if test "X$prefix" != "XNONE"
    then
      jd_prefix="$prefix"
    fi
    jd_exec_prefix="$jd_prefix"
    if test "X$exec_prefix" != "XNONE"
    then
      jd_exec_prefix="$exec_prefix"
    fi

    dnl#Unfortunately, exec_prefix may have a value like ${prefix}, etc.
    dnl#Let the shell expand those.  Yuk.
    eval `sh <<EOF
      prefix=$jd_prefix
      exec_prefix=$jd_exec_prefix
      libdir=$libdir
      includedir=$includedir
      echo jd_prefix_libdir="\$libdir" jd_prefix_incdir="\$includedir"
EOF
`
  fi
])
#}}}

AC_DEFUN([JD_GET_SYS_INCLIBS], dnl#{{{
[
  if test -x $ac_aux_dir/scripts/getsyslibs.sh
  then
    JD_SYS_INCLIBS=`$ac_aux_dir/scripts/getsyslibs.sh`
  else
    JD_SYS_INCLIBS=""
  fi
])
dnl#}}}

dnl# This macro process the --with-xxx, --with-xxxinc, and --with-xxxlib
dnl# command line arguments and returns the values as shell variables
dnl# jd_xxx_include_dir and jd_xxx_library_dir.  It does not perform any
dnl# substitutions, nor check for the existence of the supplied values.
AC_DEFUN([JD_WITH_LIBRARY_PATHS], dnl#{{{
[
 JD_UPPERCASE($1,JD_ARG1)
 jd_$1_include_dir=""
 jd_$1_library_dir=""
 if test X"$jd_with_$1_library" = X
 then
   jd_with_$1_library=""
 fi

 AC_ARG_WITH($1,
  [  --with-$1=DIR      Use DIR/lib and DIR/include for $1],
  [jd_with_$1_arg=$withval], [jd_with_$1_arg=unspecified])

 case "x$jd_with_$1_arg" in
   xno)
     jd_with_$1_library="no"
    ;;
   x)
    dnl# AC_MSG_ERROR(--with-$1 requires a value-- try yes or no)
    jd_with_$1_library="yes"
    ;;
   xunspecified)
    ;;
   xyes)
    jd_with_$1_library="yes"
    ;;
   *)
    jd_with_$1_library="yes"
    jd_$1_include_dir="$jd_with_$1_arg"/include
    jd_$1_library_dir="$jd_with_$1_arg"/lib
    ;;
 esac

 AC_ARG_WITH($1lib,
  [  --with-$1lib=DIR   $1 library in DIR],
  [jd_with_$1lib_arg=$withval], [jd_with_$1lib_arg=unspecified])
 case "x$jd_with_$1lib_arg" in
   xunspecified)
    ;;
   xno)
    ;;
   x)
    AC_MSG_ERROR(--with-$1lib requres a value)
    ;;
   *)
    jd_with_$1_library="yes"
    jd_$1_library_dir="$jd_with_$1lib_arg"
    ;;
 esac

 AC_ARG_WITH($1inc,
  [  --with-$1inc=DIR   $1 include files in DIR],
  [jd_with_$1inc_arg=$withval], [jd_with_$1inc_arg=unspecified])
 case "x$jd_with_$1inc_arg" in
   x)
     AC_MSG_ERROR(--with-$1inc requres a value)
     ;;
   xunspecified)
     ;;
   xno)
     ;;
   *)
    jd_with_$1_library="yes"
    jd_$1_include_dir="$jd_with_$1inc_arg"
   ;;
 esac
])
dnl#}}}

dnl# This function checks for the existence of the specified library $1 with
dnl# header file $2.  If the library exists, then the shell variables will
dnl# be created:
dnl#  jd_with_$1_library=yes/no,
dnl#  jd_$1_inc_file
dnl#  jd_$1_include_dir
dnl#  jd_$1_library_dir
dnl# If $3 is present, then also look in $3/include+$3/lib
AC_DEFUN([JD_CHECK_FOR_LIBRARY], dnl#{{{
[
  AC_REQUIRE([JD_EXPAND_PREFIX])dnl
  AC_REQUIRE([JD_GET_SYS_INCLIBS])dnl
  dnl JD_UPPERCASE($1,JD_ARG1)
  JD_WITH_LIBRARY_PATHS($1)
  AC_MSG_CHECKING(for the $1 library and header files $2)
  if test X"$jd_with_$1_library" != Xno
  then
    jd_$1_inc_file=$2
    dnl# jd_with_$1_library="yes"

    if test "X$jd_$1_inc_file" = "X"
    then
       jd_$1_inc_file=$1.h
    fi

    if test X"$jd_$1_include_dir" = X
    then
      inc_and_lib_dirs="\
         $jd_prefix_incdir,$jd_prefix_libdir \
	 /usr/local/$1/include,/usr/local/$1/lib \
	 /usr/local/include/$1,/usr/local/lib \
	 /usr/local/include,/usr/local/lib \
	 $JD_SYS_INCLIBS \
	 /usr/include/$1,/usr/lib \
	 /usr/$1/include,/usr/$1/lib \
	 /usr/include,/usr/lib \
	 /opt/include/$1,/opt/lib \
	 /opt/$1/include,/opt/$1/lib \
	 /opt/include,/opt/lib"

      if test X$3 != X
      then
        inc_and_lib_dirs="$3/include,$3/lib $inc_and_lib_dirs"
      fi

      case "$host_os" in
         *darwin* )
	   exts="dylib so a"
	   ;;
	 *cygwin* )
	   exts="dll.a so a"
	   ;;
	 * )
	   exts="so a"
      esac

      xincfile="$jd_$1_inc_file"
      xlibfile="lib$1"
      jd_with_$1_library="no"

      for include_and_lib in $inc_and_lib_dirs
      do
        # Yuk.  Is there a better way to set these variables??
        xincdir=`echo $include_and_lib | tr ',' ' ' | awk '{print [$]1}'`
	xlibdir=`echo $include_and_lib | tr ',' ' ' | awk '{print [$]2}'`
	found=0
	if test -r $xincdir/$xincfile
	then
	  for E in $exts
	  do
	    if test -r "$xlibdir/$xlibfile.$E"
	    then
	      jd_$1_include_dir="$xincdir"
	      jd_$1_library_dir="$xlibdir"
	      jd_with_$1_library="yes"
	      found=1
	      break
	    fi
	  done
	fi
	if test $found -eq 1
	then
	  break
	fi
      done
    fi
  fi

  if test X"$jd_$1_include_dir" != X -a X"$jd_$1_library_dir" != X
  then
    AC_MSG_RESULT(yes: $jd_$1_library_dir and $jd_$1_include_dir)
    jd_with_$1_library="yes"
    dnl#  Avoid using /usr/lib and /usr/include because of problems with
    dnl#  gcc on some solaris systems.
    JD_ARG1[]_LIB=-L$jd_$1_library_dir
    JD_ARG1[]_LIB_DIR=$jd_$1_library_dir
    if test "X$jd_$1_library_dir" = "X/usr/lib" -o "X$jd_$1_include_dir" = "X/usr/include"
    then
      JD_ARG1[]_LIB=""
    else
      JD_SET_RPATH($jd_$1_library_dir)
    fi

    JD_ARG1[]_INC=-I$jd_$1_include_dir
    JD_ARG1[]_INC_DIR=$jd_$1_include_dir
    if test "X$jd_$1_include_dir" = "X/usr/include"
    then
      JD_ARG1[]_INC=""
    fi
  else
    AC_MSG_RESULT(no)
    jd_with_$1_library="no"
    JD_ARG1[]_INC=""
    JD_ARG1[]_LIB=""
    JD_ARG1[]_INC_DIR=""
    JD_ARG1[]_LIB_DIR=""
  fi
  AC_SUBST(JD_ARG1[]_LIB)
  AC_SUBST(JD_ARG1[]_INC)
  AC_SUBST(JD_ARG1[]_LIB_DIR)
  AC_SUBST(JD_ARG1[]_INC_DIR)
])
dnl#}}}

AC_DEFUN([JD_WITH_LIBRARY], dnl#{{{
[
  JD_CHECK_FOR_LIBRARY($1, $2, $3)
  if test "$jd_with_$1_library" = "no"
  then
    AC_MSG_ERROR(unable to find the $1 library and header file $jd_$1_inc_file)
  fi
])
dnl#}}}

AC_DEFUN([JH_SDPTK_SETUP], #{{{
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
PGSBIN="\$(OMIUTIL)/TOOLKIT/bin/$OMIUTIL_SYSDIR"
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
AC_SUBST(PGSBIN)
AC_SUBST(PGSLIB)
AC_SUBST(PGSINC)
])

#}}}

AC_DEFUN([JH_LIST_SOURCES], #{{{
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

#}}}
