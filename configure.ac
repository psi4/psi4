# Process this file with autoconf to produce a configure script.
#
#
#
#

define([psi_version],[4.0.0-beta1])
define([psi_buildid],[alpha])
define([psi_bugreport],[psicode@users.sourceforge.net])
AC_PREREQ(2.57)
AC_INIT([psi],[psi_version],[psi_bugreport])

AC_CONFIG_HEADER(include/psiconfig.h)

PSI_VERSION=psi_version
PSI_BUILDID=psi_buildid
PSI_BUGREPORT=psi_bugreport
AC_SUBST(PSI_VERSION)
AC_SUBST(PSI_BUILDID)
AC_SUBST(PSI_BUGREPORT)

AC_DEFINE(PSI_VERSION,"psi_version")
AC_DEFINE(PSI_BUILDID,"psi_buildid")
AC_DEFINE(PSI_BUGREPORT,"psi_bugreport")

AC_SUBST(abs_top_srcdir)
AC_SUBST(abs_top_builddir)

if test "$srcdir" = "."; then
  echo "You are not allowed to compile in the source directory."
  echo "Create an objdir or something similar as stated in the installation manual."
  exit 1
fi

PSI_TOP_SRCDIR=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED([PSI_TOP_SRCDIR],["$PSI_TOP_SRCDIR"])
PSI_TOP_OBJDIR=`pwd`
AC_DEFINE_UNQUOTED([PSI_TOP_OBJDIR],["$PSI_TOP_OBJDIR"])

AC_CONFIG_AUX_DIR(bin)
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_TARGET
AC_PREFIX_DEFAULT(/usr/local/psi4)

AC_PROG_MAKE_SET([])

srcdir=`(cd $srcdir; pwd)`

if test -z $docdir; then
       docdir='$(datadir)/doc/$(PACKAGE)'
       AC_SUBST([docdir])
fi
if test -z $txtdir; then
       txtdir='$(docdir)/txt'
       AC_SUBST([txtdir])
fi
if test -z $htmldir; then
       htmldir='$(docdir)/html'
       AC_SUBST([htmldir])
fi
if test -z $dvidir; then
       dvidir='$(docdir)/dvi'
       AC_SUBST([dvidir])
fi
pkgdatadir='$(datadir)/psi'
AC_SUBST([pkgdatadir])

INSTALL_INCLUDE='$(top_srcdir)/bin/installinc.sh' AC_SUBST(INSTALL_INCLUDE)
MKDIRS='$(top_srcdir)/bin/mkdirs.sh' AC_SUBST(MKDIRS)

#default compiler flags
FDBG='-g' AC_SUBST(FDBG)
FOPT='-O' AC_SUBST(FOPT)
FDEF='' AC_SUBST(FDEF)
FOTH='' AC_SUBST(FOTH)
FSTATIC='' AC_SUBST(FSTATIC)
FINC='$(CPPFLAGS) $(PSICPPFLAGS)' AC_SUBST(FINC)
CDBG='-g' AC_SUBST(CDBG)
COPT='-O2' AC_SUBST(COPT)
CDEF='' AC_SUBST(CDEF)
COTH='' AC_SUBST(COTH)
CINC='$(CPPFLAGS) $(PSICPPFLAGS)' AC_SUBST(CINC)
CXXDBG='-g' AC_SUBST(CXXDBG)
CXXOPT='-O2' AC_SUBST(CXXOPT)
CXXDEF='' AC_SUBST(CXXDEF)
CXXOTH='' AC_SUBST(CXXOTH)
CXXINC='$(CPPFLAGS) $(PSICPPFLAGS) $(PYTHON_INCLUDE)' AC_SUBST(CXXINC)

########
## PSI-specific options
########

# Options for libint
LIBINT_NEW_AM=8
LIBINT_OPT_AM=4
AC_ARG_WITH(max-am-eri,[  --with-max-am-eri       maximum AM in ERIs],[
  LIBINT_NEW_AM=$((2*($withval)))
  LIBINT_OPT_AM=$(($withval))
])
AC_SUBST(LIBINT_NEW_AM)
AC_SUBST(LIBINT_OPT_AM)

# Options for libderiv
LIBDERIV_NEW_AM1=6
AC_ARG_WITH(max-am-deriv1,[  --with-max-am-deriv1    maximum AM in ERI first derivatives],[
  LIBDERIV_NEW_AM1=$((2*($withval)))
])
AC_SUBST(LIBDERIV_NEW_AM1)
LIBDERIV_NEW_AM2=6
AC_ARG_WITH(max-am-deriv2,[  --with-max-am-deriv2    maximum AM in ERI second derivatives],[
  LIBDERIV_NEW_AM2=$((2*($withval)))
])
AC_SUBST(LIBDERIV_NEW_AM2)

PSI_INTEGRALS_STANDARD=1
AC_ARG_WITH(integrals-standard,[  --with-integrals-standard     specifies the standard for the basis normalization, etc.
                                The allowed values are cca and psi (default).],[
case $withval in
  psi)
    PSI_INTEGRALS_STANDARD=0
    ;;
  cca)
    PSI_INTEGRALS_STANDARD=1
    ;;
  *)
    AC_MSG_WARN([did not recognize value $withval. will assume the default])
    ;;
esac
])
AC_DEFINE_UNQUOTED(PSI_INTEGRALS_STANDARD,$PSI_INTEGRALS_STANDARD)

########
## Check for programs
########

case $target_os in
  aix*)
    CC_LIST="cc_r gcc cc"
    CXX_LIST="xlC_r g++ c++"
    F77_LIST="xlf_r gfortran g77 ifort f77"
    ;;
  irix*)
    CC_LIST="cc gcc"
    CXX_LIST="CC g++ c++"
    F77_LIST="f77 gfortran g77 ifort"
    ;;
  *)
    CC_LIST="gcc icc cc"
    CXX_LIST="g++ c++ icpc cxx"
    F77_LIST="gfortran g77 ifort f77"
    ;;
esac


AC_ARG_WITH(cc,[  --with-cc               which C compiler to use],[
CC=$withval
])
AC_PROG_CC([$CC_LIST])

AC_ARG_WITH(cxx,[  --with-cxx              which C++ compiler to use],[
CXX=$withval
])
AC_PROG_CXX([$CXX_LIST])

AC_ARG_WITH(fc,[  --with-fc               which Fortran compiler to use],[
case $withval in
  no)
    F77=="no"
    ;;
  *)
    F77=$withval
    ;;
esac
])

if test -z "$F77"; then
  AC_PROG_F77
elif test "$F77" = "no"; then
  F77=""
fi
AC_SUBST(FC,[$F77])

## Fortran-C linking convention tests

AC_ARG_WITH(f77symbol,[  --with-f77symbol        FORTRAN77 symbol convention
                          Allowed values are:
                            lc  : lower-case
                            lcu : lower-case with underscore (default)
                            uc  : upper-case
                            ucu : upper-case with underscore],[
case $withval in
  lc)
    FC_SYMBOL=1
    ;;
  lcu)
    FC_SYMBOL=2
    ;;
  uc)
    FC_SYMBOL=3
    ;;
  ucu)
    FC_SYMBOL=4
    ;;
  *)
    AC_MSG_WARN([did not recognize value $withval. will assume lower-case with underscore])
    FC_SYMBOL=2
    ;;
esac
])

if test -z "$FC_SYMBOL"; then
  AC_MSG_CHECKING(Fortran symbols)
if test -n "$F77" ; then
  /bin/rm -f ffunc.f flink.c
  echo "      subroutine ffunc()" > ffunc.f
  echo "      return" >> ffunc.f
  echo "      end" >> ffunc.f
  $F77 -c ffunc.f 1>/dev/null 2>/dev/null
  echo "main(){ FF(); return 0; }" > flink.c
  if $CC -o flink -DFF=ffunc flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(same as C)
    FC_SYMBOL=1
  elif $CC -o flink -DFF=ffunc_ flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(lowercase with underscore)
    FC_SYMBOL=2
  elif $CC -o flink -DFF=FFUNC flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(uppercase)
    FC_SYMBOL=3
  elif $CC -o flink -DFF=FFUNC_ flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(uppercase with underscore)
    FC_SYMBOL=4
  else
    AC_MSG_RESULT(giving up)
    AC_MSG_ERROR(could not determine fortran symbol names)
  fi
  /bin/rm -f ffunc.f ffunc.o flink flink.c flink.o ffunc
else
  FC_SYMBOL=2
  AC_MSG_RESULT(No Fortran compiler provided; guessing lowercase with underscore)
fi
fi
AC_SUBST(FC_SYMBOL)

#ACX_MPI()

AC_PROG_RANLIB
AC_PROG_CPP
AC_PROG_CXXCPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_CHECK_PROGS([QHELPGENERATOR],[qhelpgenerator])
AC_CHECK_PROGS([DSYMUTIL],[dsymutil])
AC_CHECK_PROGS(AR,ar,ar)
AC_ARG_WITH(ar,
[  --with-ar               Names the archive creator.],
AR=$withval
)
AC_SUBST(AR)

AC_LANG([C++])

# I need the restrict keyword, this test handles the language set by AC_LANG
AC_C_RESTRICT

ARFLAGS=r
AC_ARG_WITH(ar-flags,
[  --with-ar-flags         Flags for the the archive creator.],
ARFLAGS=$withval
)
case $os in
  aix*)
    ARFLAGS="-X32_64 $ARFLAGS"
    ;;
esac
AC_SUBST(ARFLAGS)

AC_ARG_WITH(ld,
[  --with-ld               Names the object linker.],
LD=$withval
)
AC_SUBST(LD)

EXTRAINCDIRS=
CPPFLAGS=
AC_ARG_WITH(incdirs,
[  --with-incdirs          Specifies include directories (-Idir1 -Idir2).],
EXTRAINCDIRS=$withval
CPPFLAGS=$withval
AC_MSG_RESULT([Using extra include directories: $withval])
)
AC_SUBST(EXTRAINCDIRS)
AC_SUBST(CPPFLAGS)

AC_ARG_WITH(libs,
[  --with-libs             Specifies libraries (-llib1 -llib2).],
LDLIBS=$withval
AC_MSG_RESULT([Using extra libraries: $withval])
)
AC_SUBST(LDLIBS)



LDFLAGS=
LIBDIRS=
AC_ARG_WITH(libdirs,
[  --with-libdirs          Specifies library directories (-Ldir1 -Ldir2).],
LIBDIRS=$withval
LDFLAGS=$withval
AC_MSG_RESULT([Using extra library directories: $withval])
)
AC_SUBST(LDFLAGS)



AC_CHECK_PROGS(PERL,perl,perl) AC_SUBST(PERL)

# Check for python interpreter and libraries
#  Version 2.4 is required because we call Py_InitializeEx(0)
AC_PYTHON([2.4])



AC_PROG_LATEX
AC_PROG_LATEX2HTML
AC_PROG_DVIPS
AC_PROG_PS2PDF
AC_PROG_BIBTEX



####
# Check for doxygen and it buddies
####
AC_CHECK_PROGS(DOXYGEN,doxygen,NO)
DOT=yes
AC_ARG_WITH(dot,
[  --with-dot              Gives the path to the dot graph generator.],
DOT=$withval
)
DOT_PATH=""
if test X$DOT = Xno; then
  HAVE_DOT=NO
elif test X$DOT = Xyes; then
  AC_CHECK_PROG(HAVE_DOT,dot,YES,NO)
else
  HAVE_DOT=YES
  DOT_PATH=$DOT
fi
AC_SUBST(DOT_PATH)
AC_SUBST(HAVE_DOT)



########
## Check for libraries
########

if test ! -z "$F77"; then
  AC_F77_LIBRARY_LDFLAGS
fi
# This just implicitly checks that FLIBS can co-exist with the C linker
AC_CHECK_LIB(c,exit,[],[
FLIBS=
AC_MSG_WARN([F77 libraries don't seem to link with C programs])
],$FLIBS)

AC_ARG_WITH(blas,[  --with-blas             which blas library to use],[
case $withval in
  yes)
    BLAS='-lblas'
    ;;
  no)
    BLAS=''
    ;;
  *)
    BLAS=$withval
    ;;
esac
])

AC_ARG_WITH(lapack,[  --with-lapack           which lapack library to use],[
case $withval in
  yes)
    LAPACK='-llapack'
    ;;
  no)
    LAPACK=''
    ;;
  *)
    LAPACK=" $withval"
    ;;
esac
])

# Set strings for particular BLAS and LAPACK functions for later testing
case $FC_SYMBOL in
  1)
    DGEMM_SYMBOL=dgemm
    DGEEV_SYMBOL=dgeev
    ;;
  2)
    DGEMM_SYMBOL=dgemm_
    DGEEV_SYMBOL=dgeev_
    ;;
  3)
    DGEMM_SYMBOL=DGEMM
    DGEEV_SYMBOL=DGEEV
    ;;
  4)
    DGEMM_SYMBOL=DGEMM_
    DGEEV_SYMBOL=DGEEV_
    ;;
esac



SAVE_LIBS=$LIBS
LIBS="$LIBS $CLIBS $FLIBS"

### If the user didn't give us a BLAS lib, see if we can find one

# If were using the Intel compilers, we prefer MKL as well
if test -z "$BLAS"; then
  if test $CC = icc; then
    AC_CHECK_LIB(mkl, $DGEMM_SYMBOL, BLAS="-lmkl")
  fi
fi
# Otherwise, run through a list of possibilities
if test -z "$BLAS"; then
  AC_CHECK_LIB(essl, $DGEMM_SYMBOL, BLAS="-lessl")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(scs, $DGEMM_SYMBOL, BLAS="-lscs")
fi

if test -z "$BLAS"; then
  SAVE_LIBS2=$LIBS
  LIBS="-lpthread $LIBS"
  AC_CHECK_LIB(goto, $DGEMM_SYMBOL, BLAS="-lgoto -lpthread")
  LIBS=$SAVE_LIBS2
fi
if test -z "$BLAS"; then
  SAVE_LIBS2=$LIBS
  LIBS="-latlas $LIBS"
  AC_CHECK_LIB(f77blas, $DGEMM_SYMBOL, BLAS="-lf77blas -latlas")
  LIBS=$SAVE_LIBS2
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(mkl, $DGEMM_SYMBOL, BLAS="-lmkl")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(blas, $DGEMM_SYMBOL, BLAS="-lblas")
fi
if test -z "$BLAS"; then
case $target_vendor in
  apple)
    BLAS="-altivec -framework vecLib"
    AUTOAPPLE="yes"
    AC_MSG_RESULT([No BLAS specified and vendor is Apple: Assuming vecLib is available!])
    ;;
esac
fi
if test -z "$BLAS"; then
  AC_MSG_ERROR("Did not find a BLAS library.")
fi
AC_SUBST(BLAS)
AC_MSG_RESULT([BLAS library is $BLAS])
LIBS=$SAVE_LIBS



# Similar checks for LAPACK, except that sometimes BLAS includes LAPACK
# Also note that LAPACK often (always?) depends on BLAS
SAVE_LIBS=$LIBS
LIBS="$BLAS $LIBS $CLIBS $FLIBS"
if test -z "$LAPACK"; then
  if test ! -z "$AUTOAPPLE"; then
    AC_MSG_RESULT([Assuming vecLib for LAPACK])
    LAPACK="$BLAS"
    LAPACK_IN_BLAS=yes
  else
    AC_CHECK_FUNC($DGEEV_SYMBOL,LAPACK_IN_BLAS=yes)
  fi
fi
if test -z "$LAPACK_IN_BLAS"; then
  if test -z "$LAPACK"; then
    if test $CC = icc; then
      # Intel's libmkl_lapack64 contains the double-precision routines
      AC_CHECK_LIB(mkl_lapack64,$DGEEV_SYMBOL,LAPACK="-lmkl_lapack64")
    fi
  fi
  if test -z "$LAPACK"; then
    AC_CHECK_LIB(lapack,$DGEEV_SYMBOL,LAPACK="-llapack")
  fi
  if test -z "$LAPACK"; then
    AC_MSG_ERROR("Did not find a LAPACK library.");
  fi
  AC_MSG_RESULT([LAPACK library is $LAPACK])
else
  AC_MSG_RESULT([LAPACK library is included in BLAS+LIBS])
fi
AC_SUBST(LAPACK)
LIBS=$SAVE_LIBS


AC_CHECK_LIB([m], MKL_Free_Buffers, AC_DEFINE([HAVE_MKL]),[], $BLAS)
# The following check is disable until issues arising because of it are resolved.
#AC_CHECK_LIB([m], mkl_malloc, AC_DEFINE([HAVE_MKL_MALLOC]), [], $BLAS)


########
## Check for header files
########

# check for header files
AC_CHECK_HEADERS([errno.h fcntl.h unistd.h stdint.h sys/types.h sys/param.h sys/times.h sys/time.h utime.h time.h])
# check for C++ headers
AC_CXX_PROCESS_CHECK([cmath],[#include <cmath>],,[HAVE_CMATH=1],[HAVE_CMATH=0])
#CXXDEF="$CXXDEF -DHAVE_CMATH=$HAVE_CMATH"
AC_DEFINE_UNQUOTED(HAVE_CMATH,$HAVE_CMATH)
AC_CXX_PROCESS_CHECK([cstdio],[#include <cstdio>],,[HAVE_CSTDIO=1],[HAVE_CSTDIO=0])
#CXXDEF="$CXXDEF -DHAVE_CSTDIO=$HAVE_CSTDIO"
AC_DEFINE_UNQUOTED(HAVE_CSTDIO,$HAVE_CSTDIO)
AC_CXX_PROCESS_CHECK([cstdlib],[#include <cstdlib>],,[HAVE_CSTDLIB=1],[HAVE_CSTDLIB=0])
#CXXDEF="$CXXDEF -DHAVE_CSTDLIB=$HAVE_CSTDLIB"
AC_DEFINE_UNQUOTED(HAVE_CSTDLIB,$HAVE_CSTDLIB)
AC_CXX_PROCESS_CHECK([cstring],[#include <cstring>],,[HAVE_CSTRING=1],[HAVE_CSTRING=0])
#CXXDEF="$CXXDEF -DHAVE_CSTRING=$HAVE_CSTRING"
AC_DEFINE_UNQUOTED(HAVE_CSTRING,$HAVE_CSTRING)
AC_CXX_PROCESS_CHECK([cstddef],[#include <cstddef>],,[HAVE_CSTDDEF=1],[HAVE_CSTDDEF=0])
#CXXDEF="$CXXDEF -DHAVE_CSTDDEF=$HAVE_CSTDDEF"
AC_DEFINE_UNQUOTED(HAVE_CSTDDEF,$HAVE_CSTDDEF)

# Check for dlfcn.h
AC_CXX_PROCESS_CHECK([dlfcn.h],[#include <dlfcn.h>],,[HAVE_DLFCN_H=1],[HAVE_DLFCN_H=0])
AC_DEFINE_UNQUOTED(HAVE_DLFCN_H,$HAVE_DLFCN_H)

# Check for mpi.h
AC_CXX_PROCESS_CHECK([mpi.h],[#include <mpi.h>],,[
  # Check if we can link against mpi
  AC_MSG_CHECKING([if we can link to MPI])
  AC_TRY_LINK([#include <mpi.h>], [ MPI_Finalize(); ],
    [
      HAVE_MPI=1
      AC_MSG_RESULT([yes])
      AC_DEFINE_UNQUOTED(HAVE_MPI,$HAVE_MPI)
      AC_SUBST(HAVE_MPI)
    ],[
    AC_MSG_RESULT([no])]
  )
],[AC_MSG_RESULT([no])])


# Set up the location of the madness libraries.  The default
# is to use the psi4 internal madness.  The user can also specify
# not to use madness, or the user can specify a different madness
# library to link with.
HAVE_MADNESS=0
BUILD_MADNESS=0
if test "$HAVE_MPI" = "1"; then
  AC_ARG_WITH([madness],
  [  --with-madness          Specifies the location of the MADNESS libraries. You must be able to link to an MPI library. (default is internal).], [
  case $withval in
    yes)
      BUILD_MADNESS=1
      HAVE_MADNESS=1
      AC_DEFINE_UNQUOTED(HAVE_MADNESS,1)
      madlibs='$(top_objdir)/madness/src/lib'
      madinc='-I$(top_srcdir)/madness/src/lib -I$(top_srcdir)/madness/include -I$(top_objdir)/madness/include'
      ;;
    no)
      AC_DEFINE_UNQUOTED(HAVE_MADNESS,0)
      ;;
    *)
      HAVE_MADNESS=1
      AC_DEFINE_UNQUOTED(HAVE_MADNESS,1)
      madlibs="$withval/lib"
      madinc="-I$withval/include"
      ;;
  esac
  ],[
  BUILD_MADNESS=1
  HAVE_MADNESS=1
  AC_DEFINE_UNQUOTED(HAVE_MADNESS,1)
  madlibs='$(top_objdir)/madness/src/lib'
  madinc='-I$(top_srcdir)/madness/src/lib -I$(top_srcdir)/madness/include -I$(top_objdir)/madness/include'
  ])
fi
AC_SUBST(BUILD_MADNESS)

/bin/rm -f a.out


########
## Check for types
########



#check to see if void* is allowed
AC_MSG_CHECKING(for void pointers)
psi_ac_voidptr=yes
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
;}
void *func(void*ptr) { return ptr;
]])],[],[psi_ac_voidptr=no])
AC_MSG_RESULT($psi_ac_voidptr)
if test $psi_ac_voidptr = "no"; then
  CDEF="$CDEF -DNO_VOID_PTR"
fi

########
## Check for compiler characteristics
########

test -z "$F77SUF" && F77SUF=f AC_SUBST(F77SUF)
test -z "$LIBSUF" && LIBSUF=a AC_SUBST(LIBSUF)
test -z "$OBJSUF" && OBJSUF=o AC_SUBST(OBJSUF)

dnl -- check how dependency information is built --
# The GNU compilers work with:
CCDEPENDSUF=none
CXXDEPENDSUF=none
CCDEPENDFLAGS=-M
CXXDEPENDFLAGS=-M
/bin/rm -f depcheck.u depcheck.c depcheck.cc depcheck.o
# Check for an IBM visual age C compiler
echo "#include <stdio.h>" > depcheck.c
$CC $CPPFLAGS $CFLAGS -M depcheck.c > /dev/null 2>&1
if test -f depcheck.u; then
  CCDEPENDSUF=u
fi
/bin/rm -f depcheck.u depcheck.c depcheck.cc depcheck.o
# Check for an IBM visual age C++ compiler
echo "#include <iostream>" > depcheck.cc
$CXX $CPPFLAGS $CXXFLAGS -M -E depcheck.cc > /dev/null 2>&1
if test -f depcheck.u; then
  CXXDEPENDSUF=u
  CXXDEPENDFLAGS="-M -E"
fi
/bin/rm -f depcheck.u depcheck.c depcheck.cc depcheck.o

AC_SUBST(CCDEPENDSUF)
AC_SUBST(CXXDEPENDSUF)
AC_SUBST(CCDEPENDFLAGS)
AC_SUBST(CXXDEPENDFLAGS)

AC_CHECK_LIB(m,exit,CLIBS="$CLIBS -lm")
AC_CHECK_LIB(bsd,exit,CLIBS="$CLIBS -lbsd")
AC_MSG_RESULT(c libs: $CLIBS)
AC_SUBST(CLIBS)

AC_CHECK_LIB(m,exit,CXXLIBS="$CXXLIBS -lm")
AC_MSG_RESULT(c++ libs: $CXXLIBS)
AC_SUBST(CXXLIBS)

BOOST_ROOT='$(top_objdir)/boost/boost_1_48_0'
BOOST_INC='$(BOOST_ROOT)'
BOOST_LIB='$(BOOST_ROOT)/bin.v2/libs'
AC_SUBST(BOOST_ROOT)
AC_SUBST(BOOST_INC)
AC_SUBST(BOOST_LIB)

if test "$HAVE_MADNESS" = "1"; then
  PSICPPFLAGS='-I. -I$(BOOST_INC) -I$(top_objdir)/include -I$(top_srcdir)/include -I$(top_srcdir)/src/lib -I$(top_objdir)/src/lib $(madinc)'
  AC_SUBST(PSICPPFLAGS)
  AC_SUBST(madlibs)
  AC_SUBST(madinc)
else
  PSICPPFLAGS='-I. -I$(BOOST_INC) -I$(top_objdir)/include -I$(top_srcdir)/include -I$(top_srcdir)/src/lib -I$(top_objdir)/src/lib'
  AC_SUBST(PSICPPFLAGS)
fi
AC_SUBST(CPPFLAGS)
AC_SUBST(top_objdir)

AC_C_CONST



AC_CHECK_DECLS([putenv],[AC_DEFINE_UNQUOTED(HAVE_DECL_PUTENV,1)])
AC_CHECK_FUNC(putenv,[AC_DEFINE_UNQUOTED(HAVE_PUTENV,1)])

AC_CHECK_DECLS([setenv],[AC_DEFINE_UNQUOTED(HAVE_DECL_SETENV,1)])
AC_CHECK_FUNC(setenv,[AC_DEFINE_UNQUOTED(HAVE_SETENV,1)])

AC_CHECK_LIB(m,isinf,[AC_DEFINE_UNQUOTED(HAVE_FUNC_ISINF,1)],[
  AC_CHECK_LIB(m,finite,[AC_DEFINE_UNQUOTED(HAVE_FUNC_FINITE,1)],[
    AC_MSG_ERROR([Did not find isinf() or finite() functions.])
  ])
])

AC_CHECK_LIB(m,erf,[AC_DEFINE_UNQUOTED(HAVE_FUNC_ERF,1)])

# set up cpu specific flags
AC_MSG_CHECKING([for architecture/compiler specific optimization flags])
case $target_cpu in
  rs6000)
    if test "$FC" = xlf_r ; then
      FOPT="-qarch=com -qxflag=xalias -qdpc -qfloat=fold -qfloat=maf -qmaxmem=-1"
    fi
    if test "$CC" = xlc_r -o "$CC" = cc_r ; then
      COPT="-qarch=com -qfloat=maf -qfloat=fold -qmaxmem=-1"
    fi
    if test "$CXX" = xlC_r ; then
      CXXOPT="-qarch=com -qfloat=maf -qfloat=fold -qmaxmem=-1"
    fi
    ;;
  powerpc)
    if test "$FC" = xlf_r ; then
      FOPT="-qarch=ppc -qxflag=xalias -qdpc -qfloat=fold -qmaxmem=-1 -qfloat=maf"
    fi
    if test "$CC" = xlc_r -o "$CC" = cc_r ; then
      COPT="-qarch=ppc -qfloat=maf -qfloat=fold -qmaxmem=-1"
    fi
    if test "$CXX" = xlC_r ; then
      CXXOPT="-qarch=ppc -qfloat=maf -qfloat=fold -qmaxmem=-1"
    fi
    ;;
  i686)
    if test "$FC" = g77 ; then
      FOPT=""
    fi
    if test "$CC" = gcc -o "$CC" = cc ; then
      COPT=""
    fi
    if test "$CXX" = g++ -o "$CXX" = c++ -o "$CXX" = gcc ; then
      CXXOPT=""
    fi
    ;;
  x86_64)
    if test "$CC" = gcc -o "$CC" = cc ; then
      COPT="-m64"
    fi
    if test "$CXX" = g++ -o "$CXX" = c++ -o "$CXX" = gcc ; then
      CXXOPT="-m64"
    fi
    ;;
esac
AC_MSG_RESULT([C optimization flags: $COPT])
AC_MSG_RESULT([C++ optimization flags: $CXXOPT])

AC_MSG_CHECKING([for OS specific flags])
case $target_os in
  aix*)
# _ALL_SOURCE is required because some include files on AIX don't seem to work with pure POSIX
    EXTRADEFS="-DAIX -D_LARGE_FILES -D_ALL_SOURCE"
    FDEF="$FDEF $EXTRADEFS"
    CDEF="$CDEF $EXTRADEFS"
    CXXDEF="$CXXDEF $EXTRADEFS"
    ;;
  linux*)
    EXTRADEFS="-DLinux -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE"
    FDEF="$FDEF $EXTRADEFS"
    CDEF="$CDEF $EXTRADEFS"
    CXXDEF="$CXXDEF $EXTRADEFS"
    ;;
esac
AC_MSG_RESULT([C preprocessor directives: $CDEF])
AC_MSG_RESULT([C++ preprocessor directives: $CXXDEF])

# Adjust the debug flags.
# The default is to not debug.  --with-debug will turn on debugging.
AC_ARG_WITH(debug,[  --with-debug            will turn on debugging],[
case $withval in
  yes)
        CDBG='-g -DPSIDEBUG'
        CXXDBG='-g -DPSIDEBUG'
    ;;
  no)
    CDBG=''
    CXXDBG=''
    ;;
  *)
    CDBG="$withval"
    CXXDBG="$withval"
    ;;
esac
],CDBG='' CXXDBG='')

# Adjust the optimize flags.
# The default is to optimize.  --with-opt=no or --without-opt will turn
# off optimization.
AC_ARG_WITH(opt,[  --with-opt              --with-opt=no to turn off optimization],[
case $withval in
  yes)
    ;;
  no)
    COPT='-O0'
    CXXOPT='-O0'
    ;;
  *)
    COPT="$withval"
    CXXOPT="$withval"
    ;;
esac
])

# Check for OpenMP, the test reports for C++ only, because we called AC_LANG([C++]) near the top
AX_OPENMP([
    CXXOTH="$CXXOTH $OPENMP_CXXFLAGS"
    LDFLAGS="$LDFLAGS $OPENMP_CXXFLAGS"
])

# Check to see if the compiler supports __builtin_expect and __builtin_prefetch
AX_BUILTINS()

# Adjust the flags dealing with strictness of the compilers.
# The default is to disable strict language checks.
AC_ARG_WITH(strict,[  --with-strict           --with-strict=yes to turn on strict compiler warnings],[
case $withval in
  yes)
    # "Strict" flags only known for GNU compilers
    if test "$G77" = yes ; then
      FOTH="$FOTH -Wall"
    fi
    if test "$GCC" = yes ; then
      COTH="$COTH -Wall -ansi -pedantic -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-long-long"
    fi
    if test "$GXX" = yes ; then
      CXXOTH="$CXXOTH -Wall -ansi -pedantic -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-long-long"
    fi
    ;;
  no)
    ;;
esac
])

AC_SUBST(TARGET_OS, $target_os)
# Handle plugins
AC_ARG_WITH(plugins,[  --with-plugins          Compiles PSI4 to be compatible with plugins.],[
case $withval in
  yes)
     AC_SUBST(HAVE_PLUGINS, yes)
     case $target_os in
       linux*)
         CXXOTH="$CXXOTH -fpic"
         LDFLAGS="$LDFLAGS -Wl,-export-dynamic"
         ;;
       darwin*)
         CXXOTH="$CXXOTH -fpic"
         LDFLAGS="$LDFLAGS -all_load"
         ;;
     esac
     ;;
  no)
     AC_SUBST(HAVE_PLUGINS, no)
     ;;
esac
])

# Adjust the flags dealing with strictness of the compilers.
# The default is to disable strict language checks.
AC_ARG_WITH(insanely-strict,[  --with-insanely-strict  --with-insanely-strict=yes to turn on strict compiler warnings (cover your eyes, it ain't pretty, for now)],[
case $withval in
  yes)
    # "Strict" flags only known for GNU compilers
    if test "$GXX" = yes ; then
      CXXOTH="$CXXOTH -pipe -Wall -ansi -pedantic -W -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wno-long-long -Wctor-dtor-privacy -Woverloaded-virtual -Wextra -Wwrite-strings -Winit-self -Wold-style-cast -Wstrict-aliasing -Wformat=2 -Wuninitialized -Wnon-virtual-dtor"
    fi
    ;;
  no)
    ;;
esac
])
AC_MSG_RESULT([C extra flags: $COTH])
AC_MSG_RESULT([C++ extra flags: $CXXOTH])

# Set up CFLAGS and FFLAGS.
CFLAGS='$(CDBG) $(COPT) $(CDEF) $(COTH) $(CINC)' AC_SUBST(CFLAGS)
CXXFLAGS='$(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(CXXINC)' AC_SUBST(CXXFLAGS)
FFLAGS='$(FDBG) $(FOPT) $(FOTH) $(FSTATIC)' AC_SUBST(FFLAGS)

########
## Check for library functions
########


########
## Output
######## 

AC_CONFIG_FILES([
          boost/Makefile
          include/Makefile
          Makefile lib/Makefile src/Makefile
          src/lib/Makefile src/lib/MakeRules src/lib/MakeVars
          src/util/Makefile src/util/MakeRules src/util/MakeVars
          src/bin/Makefile src/bin/MakeRules src/bin/MakeVars

          doc/MakeVars doc/MakeRules doc/Makefile
          doc/macros.tex doc/progman/svn/svn.eps
          doc/installation/Makefile
          doc/progman/Makefile doc/userman/Makefile
          doc/userman/prsty.bst
          doc/doxygen/Makefile
          doc/doxygen/doxygen.cfg

          src/lib/libdpd/Makefile
          src/lib/libutil/Makefile
          src/lib/libiwl/Makefile
          src/lib/libpsio/Makefile
          src/lib/libqt/Makefile
          src/lib/libciomr/Makefile
          src/lib/libchkpt/Makefile
          src/lib/libderiv/Makefile
          src/lib/libderiv/MakeVars
          src/lib/libderiv/libderiv_config.h
          src/lib/libdiis/Makefile
          src/lib/libint/Makefile
          src/lib/libint/MakeVars
          src/lib/libint/libint_config.h
          src/lib/liboptions/Makefile
          src/lib/libmoinfo/Makefile
          src/lib/libmints/Makefile
          src/lib/libdist_matrix/Makefile
          src/lib/libscf_solver/Makefile
          src/lib/liblmp2_solver/Makefile
          src/lib/libtrans/Makefile
          src/lib/libparallel/Makefile
          src/lib/libsapt_solver/Makefile
          src/lib/libfunctional/Makefile
          src/lib/libplugin/Makefile
          src/lib/lib3index/Makefile
          src/lib/libsmartptr/Makefile
          src/lib/libyeti/Makefile
          src/lib/libdfcc_solver/Makefile
          src/util/tocprint/Makefile

          src/bin/psi4/Makefile
          src/bin/response/Makefile
          src/bin/dcft/Makefile
          src/bin/mcscf/Makefile
          src/bin/transqt2/Makefile
          src/bin/transqt/Makefile
          src/bin/detci/Makefile
          src/bin/mp2/Makefile
          src/bin/clag/Makefile
          src/bin/stable/Makefile
          src/bin/ccsort/Makefile
          src/bin/ccenergy/Makefile
          src/bin/cctriples/Makefile
          src/bin/ccresponse/Makefile
          src/bin/cchbar/Makefile
          src/bin/cclambda/Makefile
          src/bin/ccdensity/Makefile
          src/bin/cceom/Makefile
          src/bin/cis/Makefile
          src/bin/scf/Makefile
          src/bin/lmp2/Makefile
          src/bin/psimrcc/Makefile
          src/bin/mrcc/Makefile
          src/bin/df-mp2/Makefile
          src/bin/sapt/Makefile
          src/bin/mints/Makefile
          src/bin/deriv/Makefile
          src/bin/optking/Makefile
          src/bin/df-cc/Makefile
          src/bin/findif/Makefile
          src/bin/omp2/Makefile
          src/bin/adc/Makefile

          tests/Makefile
          tests/MakeVars
          tests/MakeRules
          tests/runtest.pl
          tests/tu1-h2o-energy/Makefile
          tests/tu2-ch2-energy/Makefile
          tests/tu3-h2o-opt/Makefile
          tests/tu4-h2o-freq/Makefile
          tests/tu5-sapt/Makefile
          tests/tu6-cp-ne2/Makefile
          tests/dcft1/Makefile
          tests/dcft2/Makefile
          tests/dcft3/Makefile
          tests/dcft4/Makefile
          tests/dcft5/Makefile
          tests/scf1/Makefile
          tests/scf2/Makefile
          tests/scf3/Makefile
          tests/scf4/Makefile
          tests/scf5/Makefile
          tests/dfscf-bz2/Makefile
          tests/scf-bz2/Makefile
          tests/opt1/Makefile
          tests/opt1-fd/Makefile
          tests/opt2/Makefile
          tests/opt2-fd/Makefile
          tests/opt3/Makefile
          tests/opt4/Makefile
          tests/opt5/Makefile
          tests/sad1/Makefile
          tests/fd-gradient/Makefile
          tests/fd-freq-energy/Makefile
          tests/castup1/Makefile
          tests/sapt1/Makefile
          tests/dfmp2_1/Makefile
          tests/dfmp2_2/Makefile
          tests/pywrap_db1/Makefile
          tests/pywrap_db2/Makefile
          tests/pywrap_cbs1/Makefile
          tests/pywrap_all/Makefile
          tests/pywrap_alias/Makefile
          tests/cisd-h2o+-0/Makefile
          tests/cisd-h2o+-1/Makefile
          tests/cisd-h2o+-2/Makefile
          tests/cisd-h2o-clpse/Makefile
          tests/cisd-opt-fd/Makefile
          tests/cisd-sp/Makefile
          tests/cisd-sp-2/Makefile
          tests/fci-h2o/Makefile
          tests/fci-h2o-2/Makefile
          tests/fci-h2o-fzcv/Makefile
          tests/fci-dipole/Makefile
          tests/fci-tdm/Makefile
          tests/fci-tdm-2/Makefile
          tests/rasci-c2-active/Makefile
          tests/rasci-h2o/Makefile
          tests/rasci-ne/Makefile
          tests/zaptn-nh2/Makefile
          tests/mpn-bh/Makefile
          tests/cc1/Makefile
          tests/cc2/Makefile
          tests/cc4/Makefile
          tests/cc4a/Makefile
          tests/cc5a/Makefile
          tests/cc8/Makefile
          tests/cc8a/Makefile
          tests/cc8b/Makefile
          tests/cc8c/Makefile
          tests/cc9/Makefile
          tests/cc9a/Makefile
          tests/cc10/Makefile
          tests/cc11/Makefile
          tests/cc12/Makefile
          tests/cc13/Makefile
          tests/cc13a/Makefile
          tests/cc15/Makefile
          tests/cc16/Makefile
          tests/cc17/Makefile
          tests/cc18/Makefile
          tests/cc19/Makefile
          tests/cc22/Makefile
          tests/cc24/Makefile
          tests/cc25/Makefile
          tests/cc26/Makefile
          tests/cc27/Makefile
          tests/cc28/Makefile
          tests/cc29/Makefile
          tests/cc30/Makefile
          tests/cc31/Makefile
          tests/cc32/Makefile
          tests/cc33/Makefile
          tests/cc34/Makefile
          tests/cc35/Makefile
          tests/cc36/Makefile
          tests/cc37/Makefile
          tests/cc38/Makefile
          tests/cc39/Makefile
          tests/cc40/Makefile
          tests/cc41/Makefile
          tests/cc42/Makefile
          tests/cc43/Makefile
          tests/cc44/Makefile
          tests/cc45/Makefile
          tests/cc49/Makefile
          tests/cc50/Makefile
          tests/cc51/Makefile
          tests/cc52/Makefile
          tests/mcscf1/Makefile
          tests/mcscf2/Makefile
          tests/mcscf3/Makefile
          tests/psimrcc-sp1/Makefile
          tests/mints1/Makefile
          tests/mints2/Makefile
          tests/mints3/Makefile
          tests/mints4/Makefile
          tests/mp2_1/Makefile
          tests/matrix1/Makefile
          tests/props1/Makefile
          tests/props2/Makefile
          tests/dft1/Makefile
          tests/mom/Makefile
          tests/plugin_aointegrals/Makefile
          tests/plugin_aointegrals/input.dat
          tests/plugin_aointegrals/aointegrals.cc
          tests/plugin_sointegrals/Makefile
          tests/plugin_sointegrals/input.dat
          tests/plugin_sointegrals/sointegrals.cc
          tests/plugin_backtrans/Makefile
          tests/plugin_backtrans/input.dat
          tests/plugin_backtrans/backtrans.cc
          tests/plugin_mp2/Makefile
          tests/plugin_mp2/input.dat tests/plugin_mp2/globals.h
          tests/plugin_mp2/plugin_mp2.cc tests/plugin_mp2/plugin_mp2_restricted.cc
          tests/plugin_mp2/plugin_mp2_unrestricted.cc
          tests/plugin_mointegrals/Makefile
          tests/plugin_mointegrals/input.dat
          tests/plugin_mointegrals/mointegrals.cc

          tests/plugin_testparse/Makefile
          tests/plugin_testparse/input.dat
          tests/plugin_testparse/plugin_testparse.cc

          bin/psi4-config
])

AC_CONFIG_COMMANDS([psi4-config], [chmod +x bin/psi4-config])

# If we are using madness, and we are linking the internal madness
# library then set up the internal madness to be built
if test "$HAVE_MADNESS" = "1"; then
  if test "$BUILD_MADNESS" = "1"; then
    AC_CONFIG_SUBDIRS([madness])
  fi
fi

AC_OUTPUT
