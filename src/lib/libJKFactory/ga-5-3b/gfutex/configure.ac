# Process this file with autoconf to produce a configure script.

###############################################################################
## Init autoconf
################################################################################
AC_PREREQ([2.68])
AC_INIT([Global Futures], [0.1], [hpctools@pnnl.gov], [gfutex])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/gfutex.C])

###############################################################################
# Must occur before automake init
###############################################################################
GA_TARGET
GA_CROSS_COMPILING

###############################################################################
## Init automake
################################################################################
AM_INIT_AUTOMAKE([color-tests foreign parallel-tests silent-rules subdir-objects])

## Don't emit "rebuild rules" for configure, Makefile.ins, etc.
AM_MAINTAINER_MODE

###############################################################################
# Misc. information and package setup.
###############################################################################
GA_WITH_HELP

GFUTEX_TOP_BUILDDIR="`pwd`"
cd "$srcdir"
GFUTEX_TOP_SRCDIR="`pwd`"
cd "$GFUTEX_TOP_BUILDDIR"
AS_IF([test "$GFUTEX_TOP_BUILDDIR" != "GFUTEX_TOP_SRCDIR"],
      [AC_MSG_NOTICE([Detected VPATH build])])

AC_CACHE_CHECK([whether gfutex is built outside of GA], [gfutex_cv_standalone],
    [AS_IF([test -f "$GFUTEX_TOP_BUILDDIR/../armci/config.status"],
        [gfutex_cv_standalone=no],
        [gfutex_cv_standalone=yes])])
AM_CONDITIONAL([GFUTEX_STANDALONE], [test x$gfutex_cv_standalone = xyes])

# If MPI is desired we use the MPI compiler wrappers instead of the standard
# compilers.
GA_MSG_COMMS([no-tcgmsg])

# F77 and C++ are always enabled for gfutex
enable_f77=yes
enable_cxx=yes

GA_DISABLE_MPI_TESTS

###############################################################################
# C compiler
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C compiler])
AC_MSG_NOTICE

tascel_save_CFLAGS="$CFLAGS"
GA_PROG_MPICC
CFLAGS="$tascel_save_CFLAGS"
AS_IF([test x$with_mpi_wrappers = xyes],
    [GA_MPI_UNWRAP],
    [GA_ARG_PARSE([with_mpi], [GA_MP_LIBS], [GA_MP_LDFLAGS], [GA_MP_CPPFLAGS])])

GA_MPICC_TEST_COMPILE

# Hack to append .x to executables.
AC_SUBST([EXEEXT], [.x$EXEEXT])

GA_TARGET64

# Checks for C libraries.

# Checks for C header files.
AC_HEADER_ASSERT
GA_CHECK_HEADERS([limits.h math.h omp.h pthread.h sched.h stdint.h stdio.h stdlib.h string.h sys/time.h unistd.h])

# Checks for C typedefs, structures, and compiler characteristics.

# Checks for C library functions.

###############################################################################
# C++ compiler
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C++ compiler])
AC_MSG_NOTICE
AC_LANG_PUSH([C++])

gfutex_save_CXXFLAGS="$CXXFLAGS"
GA_PROG_MPICXX
CXXFLAGS="$gfutex_save_CXXFLAGS"
AS_IF([test x$with_mpi_wrappers = xyes], [GA_MPI_UNWRAP])
GA_MPICXX_TEST

# Checks for libraries.
GA_TBB
GA_BOOST

# Checks for C++ header files.

# Checks for C++ typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for C++ type sizes.
AC_CHECK_SIZEOF([void*])
AC_CHECK_SIZEOF([char])
AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AS_IF([test x$ac_cv_type_long_long_int = xyes], [AC_CHECK_SIZEOF([long long])])
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])
AS_IF([test x$ac_cv_type_long_double = xyes], [AC_CHECK_SIZEOF([long double])])

# Checks for C++ library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday memset pow sqrt strerror])

AC_LANG_POP([C++])

###############################################################################
# Look for Global Arrays
###############################################################################
GA_REQUIRE

###############################################################################
# Libtool setup -- no compiler/linker tests after this
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Libtool setup])
AC_MSG_NOTICE
# temporarily restore unwrapped compilers
# this works around a bug where libtool sadly relies on matching compiler
# names in order to determine features (Fortran only, I think)
# libtool doesn't recognize MPI compiler names, nor should it
AS_IF([test x$with_mpi_wrappers = xyes], [GA_MPI_UNWRAP_PUSH])
GA_AR
LT_INIT([disable-shared])
# and now that that's over, put the MPI compilers back
# also, the above hack incorrectly sets the base compiler as the linker
AS_IF([test x$with_mpi_wrappers = xyes],
    [GA_MPI_UNWRAP_POP
     compiler="$CC"
     LTCC="$CC"
     lt_save_CC="$CC"
     compiler_DEFAULT="$CC"
     compiler_CXX="$CXX"
     compiler_F77="$F77"])

###############################################################################
# The End
###############################################################################
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# Report on what we found.
AC_MSG_NOTICE([])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([ $PACKAGE_NAME configured as follows:])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([    GA_MP_LIBS=$GA_MP_LIBS])
AC_MSG_NOTICE([ GA_MP_LDFLAGS=$GA_MP_LDFLAGS])
AC_MSG_NOTICE([GA_MP_CPPFLAGS=$GA_MP_CPPFLAGS])
AC_MSG_NOTICE([            CC=$CC])
AC_MSG_NOTICE([        CFLAGS=$CFLAGS])
AC_MSG_NOTICE([           CXX=$CXX])
AC_MSG_NOTICE([      CXXFLAGS=$CXXFLAGS])
AC_MSG_NOTICE([           CPP=$CPP])
AC_MSG_NOTICE([      CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([       LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([          LIBS=$LIBS])
AC_MSG_NOTICE([            AR=$AR])
AC_MSG_NOTICE([      AR_FLAGS=$AR_FLAGS])
AC_MSG_NOTICE([   GA_CPPFLAGS=$GA_CPPFLAGS])
AC_MSG_NOTICE([    GA_LDFLAGS=$GA_LDFLAGS])
AC_MSG_NOTICE([       GA_LIBS=$GA_LIBS $GA_FLIBS])
AC_MSG_NOTICE([          DEFS=$DEFS])
AC_MSG_NOTICE([         SHELL=$SHELL])
AC_MSG_NOTICE([])
