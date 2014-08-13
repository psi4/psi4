# Process this file with autoconf to produce a configure script.

###############################################################################
# Init autoconf
###############################################################################
AC_PREREQ([2.67])
AC_INIT([Communication Runtime for Exascale (comex)],
        [1.0],
        [hpctools@pnl.gov],
        [comex],
        [http://www.emsl.pnl.gov/docs/global/])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_LIBOBJ_DIR([compat])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src-common/comex.h])

###############################################################################
# Init automake
###############################################################################
AM_INIT_AUTOMAKE([color-tests foreign parallel-tests silent-rules subdir-objects])

# Don't emit "rebuild rules" for configure, Makefile.ins, etc.
AM_MAINTAINER_MODE

###############################################################################
# Misc. information and package setup.
###############################################################################
COMEX_WITH_HELP

COMEX_TOP_BUILDDIR="`pwd`"
cd "$srcdir"
COMEX_TOP_SRCDIR="`pwd`"
cd "$COMEX_TOP_BUILDDIR"
AS_IF([test "$COMEX_TOP_BUILDDIR" != "$COMEX_TOP_SRCDIR"],
      [AC_MSG_NOTICE([Detected VPATH build])])

# We use the MPI compiler wrappers instead of the standard compilers.
COMEX_WITH_MPI

# In case users don't care whether MPI actually works during configure.
COMEX_DISABLE_MPI_TESTS

# ARMCI profiling layer
COMEX_ENABLE_PROFILING_ARMCI

#########################################
# C compiler
########################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C compiler])
AC_MSG_NOTICE

COMEX_PROG_MPICC
AC_USE_SYSTEM_EXTENSIONS

AS_IF([test x$with_mpi_wrappers = xyes],
    [COMEX_MPI_UNWRAP],
    [COMEX_ARG_PARSE([with_mpi], [MPI_LIBS], [MPI_LDFLAGS], [MPI_CPPFLAGS])])
AS_CASE([$enable_mpi_tests],
    [yes],[COMEX_MPICC_TEST_LINK],
    [no], [COMEX_MPICC_TEST_COMPILE])

AM_CONDITIONAL([CROSS_COMPILING], [test "x$cross_compiling" = xyes])

# Establish the underlying network infrastructure (MPI, OFA, DMAPP, etc)
COMEX_NETWORK_SETUP

# Checks for C header files.
AC_HEADER_ASSERT
AC_HEADER_DIRENT
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
COMEX_CHECK_HEADERS([assert.h])
COMEX_CHECK_HEADERS([errno.h])
COMEX_CHECK_HEADERS([getopt.h])
COMEX_CHECK_HEADERS([math.h])
COMEX_CHECK_HEADERS([pthread.h])
COMEX_CHECK_HEADERS([semaphore.h])
COMEX_CHECK_HEADERS([stdint.h])
COMEX_CHECK_HEADERS([stdio.h])
COMEX_CHECK_HEADERS([stdlib.h])
COMEX_CHECK_HEADERS([string.h])
COMEX_CHECK_HEADERS([strings.h])
COMEX_CHECK_HEADERS([sys/time.h])
COMEX_CHECK_HEADERS([sys/types.h])
COMEX_CHECK_HEADERS([unistd.h])

# Checks for C typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_C_VOLATILE
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_CHECK_TYPES([ptrdiff_t])
COMEX_FUNCTION
COMEX_DISABLE_SYS_WEAK_ALIAS
COMEX_SYS_WEAK_ALIAS

# Checks for C type sizes.
comex_save_LIBS="$LIBS"; LIBS="$LIBS $MPI_LIBS"
comex_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$LDFLAGS $MPI_LDFLAGS"
comex_save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $MPI_CPPFLAGS"
AC_CHECK_SIZEOF([MPI_Aint], [], [[#include <mpi.h>]])
LIBS="$comex_save_LIBS"
LDFLAGS="$comex_save_LDFLAGS"
CPPFLAGS="$comex_save_CPPFLAGS"
AC_CHECK_SIZEOF([void*])
AC_CHECK_SIZEOF([char])
AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])

# Checks for C library functions.
AC_FUNC_MMAP
COMEX_SEARCH_LIBS([sqrt], [m])
COMEX_SEARCH_LIBS([sem_open], [rt])
COMEX_CHECK_FUNCS([bzero])
COMEX_CHECK_FUNCS([gettimeofday])
COMEX_CHECK_FUNCS([memset])
COMEX_CHECK_FUNCS([munmap])
COMEX_CHECK_FUNCS([strchr])
COMEX_CHECK_FUNCS([strdup])
COMEX_CHECK_FUNCS([strncasecmp])
COMEX_CHECK_FUNCS([strstr])

# Checks for C libraries.
COMEX_BLAS

###############################################################################
# Checks for programs.
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Checks for additional programs])
AC_MSG_NOTICE

AC_PROG_GREP
AC_PROG_SED

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
AS_IF([test x$with_mpi_wrappers = xyes], [COMEX_MPI_UNWRAP_PUSH])
COMEX_AR
LT_INIT([disable-shared])
# and now that that's over, put the MPI compilers back
# also, the above hack incorrectly sets the base compiler as the linker
AS_IF([test x$with_mpi_wrappers = xyes],
    [COMEX_MPI_UNWRAP_POP
     compiler="$CC"
     LTCC="$CC"
     lt_save_CC="$CC"
     compiler_DEFAULT="$CC"])

###############################################################################
# Test suite setup
###############################################################################
AC_ARG_VAR([NPROCS], [number of procs to use for parallel tests (default 4)])
AS_IF([test "x$NPROCS" = x], [NPROCS=4])
AC_SUBST([NPROCS])
AC_ARG_VAR([MPIEXEC],
           [how to run parallel tests if built with MPI e.g. "mpiexec -np %NP%"])
AS_IF([test "x$MPIEXEC" = x],
    [AC_PATH_PROGS([MPIEXEC], [mpirun mpiexec])
     MPIEXEC="$MPIEXEC -n %NP%"])
AC_SUBST([MPIEXEC])

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
AC_MSG_NOTICE([              MPI_LIBS=$MPI_LIBS])
AC_MSG_NOTICE([           MPI_LDFLAGS=$MPI_LDFLAGS])
AC_MSG_NOTICE([          MPI_CPPFLAGS=$MPI_CPPFLAGS])
AC_MSG_NOTICE([         COMEX_NETWORK=$comex_network])
AC_MSG_NOTICE([ COMEX_NETWORK_LDFLAGS=$COMEX_NETWORK_LDFLAGS])
AC_MSG_NOTICE([    COMEX_NETWORK_LIBS=$COMEX_NETWORK_LIBS])
AC_MSG_NOTICE([COMEX_NETWORK_CPPFLAGS=$COMEX_NETWORK_CPPFLAGS])
AC_MSG_NOTICE([                    CC=$CC])
AS_IF([test "x$with_mpi_wrappers" = xyes], [
AC_MSG_NOTICE([          unwrapped CC=$comex_cv_mpic_naked])
])
AC_MSG_NOTICE([                CFLAGS=$CFLAGS])
AC_MSG_NOTICE([                   CPP=$CPP])
AC_MSG_NOTICE([              CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([               LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([                  LIBS=$LIBS])
AC_MSG_NOTICE([                 FLIBS=$FLIBS])
AC_MSG_NOTICE([                    AR=$AR])
AC_MSG_NOTICE([              AR_FLAGS=$AR_FLAGS])
AC_MSG_NOTICE([                  DEFS=$DEFS])
AC_MSG_NOTICE([                 SHELL=$SHELL])
AC_MSG_NOTICE([               MPIEXEC=$MPIEXEC])
AC_MSG_NOTICE([                NPROCS=$NPROCS])
AC_MSG_NOTICE([])
