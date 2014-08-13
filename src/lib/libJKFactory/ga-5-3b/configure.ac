# Process this file with autoconf to produce a configure script.

###############################################################################
# Init autoconf
###############################################################################
AC_PREREQ([2.67])
AC_INIT([Global Arrays (GA)], [5.1],
        [hpctools@pnl.gov], [ga],
        [http://www.emsl.pnl.gov/docs/global/])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_LIBOBJ_DIR([compat])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([global/src/ga.h])

###############################################################################
# Must occur before automake init
###############################################################################
GA_TARGET
GA_CROSS_COMPILING

###############################################################################
# Init automake
###############################################################################
AM_INIT_AUTOMAKE([color-tests foreign parallel-tests silent-rules subdir-objects])

# Don't emit "rebuild rules" for configure, Makefile.ins, etc.
AM_MAINTAINER_MODE

###############################################################################
# Misc. information and package setup.
###############################################################################
GA_WITH_HELP

GA_TOP_BUILDDIR="`pwd`"
cd "$srcdir"
GA_TOP_SRCDIR="`pwd`"
cd "$GA_TOP_BUILDDIR"
AS_IF([test "$GA_TOP_BUILDDIR" != "$GA_TOP_SRCDIR"],
      [AC_MSG_NOTICE([Detected VPATH build])])

# Determine messaging library up front because if MPI is desired we use the
# MPI compiler wrappers instead of the standard compilers.
GA_MSG_COMMS

# Whether to disable/enable Fortran 77 or C++ bindings.
# Note that compiling the Fortran 77 source is enabled by default.
GA_F77_DISABLE
GA_CXX_ENABLE

GA_DISABLE_MPI_TESTS
GA_ENABLE_WARNINGS
GA_ENABLE_CHECKPOINT
GA_ENABLE_OPT
GA_ENABLE_PEIGS
GA_ENABLE_PROFILING
GA_ENABLE_TRACE
GA_THREAD_SAFE

# Unless we rename all preprocessor symbols of MPI and MPIPP
# in the source code, this hack is necessary.
AH_BOTTOM([/* DO NOT MODIFY THE FOLLOWING PREPROCESSOR DIRECTIVES */
#ifdef MSG_COMMS_MPI
#   ifdef __cplusplus
#       define MPIPP
#   else
#       define MPI
#   endif /* __cplusplus */
#endif /* MSG_COMMS_MPI */])

#########################################
# C compiler
#########################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C compiler])
AC_MSG_NOTICE

ga_save_CFLAGS="$CFLAGS"
GA_PROG_MPICC
CFLAGS="$ga_save_CFLAGS"
AC_USE_SYSTEM_EXTENSIONS
GA_COMPILER_VENDOR
GA_WARN_FLAGS
ga_save_CFLAGS="$CFLAGS"
CFLAGS="$ga_cv_c_warning_flags $CFLAGS"

AS_IF([test x$with_mpi_wrappers = xyes],
    [GA_MPI_UNWRAP],
    [GA_ARG_PARSE([with_mpi], [GA_MP_LIBS], [GA_MP_LDFLAGS], [GA_MP_CPPFLAGS])])
# If F77 is enabled, only perform a simple compilation test.
AS_IF([test "x$with_mpi" != xno], [
    AS_CASE([$enable_f77:$enable_mpi_tests],
        [yes:*], [GA_MPICC_TEST_COMPILE],
        [no:yes],[GA_MPICC_TEST_LINK],
        [no:no], [GA_MPICC_TEST_COMPILE])])
# Hack to append .x to executables.
AC_SUBST([EXEEXT], [.x$EXEEXT])

GA_TARGET64

# Establish the underlying network infrastructure (SOCKETS, OPENIB, etc)
GA_ARMCI_NETWORK

# Checks for C header files.
AC_HEADER_ASSERT
AC_HEADER_DIRENT
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
GA_CHECK_HEADERS([assert.h])
GA_CHECK_HEADERS([errno.h])
GA_CHECK_HEADERS([fcntl.h])
GA_CHECK_HEADERS([ffio.h])
GA_CHECK_HEADERS([limits.h])
GA_CHECK_HEADERS([linux/limits.h])
GA_CHECK_HEADERS([malloc.h])
GA_CHECK_HEADERS([math.h])
GA_CHECK_HEADERS([stddef.h])
GA_CHECK_HEADERS([stdint.h])
GA_CHECK_HEADERS([stdio.h])
GA_CHECK_HEADERS([stdlib.h])
GA_CHECK_HEADERS([string.h])
GA_CHECK_HEADERS([strings.h])
GA_CHECK_HEADERS([sys/file.h])
GA_CHECK_HEADERS([sys/fsid.h])
GA_CHECK_HEADERS([sys/fstyp.h])
GA_CHECK_HEADERS([sys/ioctl.h])
GA_CHECK_HEADERS([sys/iosw.h])
GA_CHECK_HEADERS([sys/ipc.h])
GA_CHECK_HEADERS([sys/mman.h])
GA_CHECK_HEADERS([sys/mount.h])
GA_CHECK_HEADERS([sys/param.h])
GA_CHECK_HEADERS([sys/sem.h])
GA_CHECK_HEADERS([sys/shm.h])
GA_CHECK_HEADERS([sys/statfs.h])
GA_CHECK_HEADERS([sys/stat.h])
GA_CHECK_HEADERS([sys/statvfs.h])
GA_CHECK_HEADERS([sys/time.h])
GA_CHECK_HEADERS([sys/types.h])
GA_CHECK_HEADERS([sys/vfs.h])
GA_CHECK_HEADERS([time.h])
GA_CHECK_HEADERS([unistd.h])
GA_CHECK_HEADERS([windows.h])

# Checks for C typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_C_VOLATILE
AC_TYPE_INT64_T
AC_TYPE_LONG_DOUBLE
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIGNAL
AC_TYPE_SIZE_T
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
GA_FUNCTION
GA_PROGNAME
GA_DISABLE_SYS_WEAK_ALIAS
GA_SYS_WEAK_ALIAS
GA_UNION_SEMUN

# Checks for C type sizes.
AC_CHECK_SIZEOF([void*])
AC_CHECK_SIZEOF([char])
AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])
AS_IF([test x$ac_cv_type_long_double = xyes], [AC_CHECK_SIZEOF([long double])])
GA_C_POINTER_AS_INTEGER

# Checks for C library functions.
AC_FUNC_FORK
AC_FUNC_MMAP
AC_FUNC_SELECT_ARGTYPES
AC_SEARCH_LIBS([floor], [m], [have_floor=1], [have_floor=0])
AC_DEFINE_UNQUOTED([HAVE_FLOOR], [$have_floor],
    [Define to 1 if you have the 'floor' function.])
AC_SEARCH_LIBS([modf], [m], [have_modf=1], [have_modf=0])
AC_DEFINE_UNQUOTED([HAVE_MODF], [$have_modf],
    [Define to 1 if you have the 'modf' function.])
AC_SEARCH_LIBS([pow], [m], [have_pow=1], [have_pow=0])
AC_DEFINE_UNQUOTED([HAVE_POW], [$have_pow],
    [Define to 1 if you have the 'pow' function.])
AC_SEARCH_LIBS([sqrt], [m], [have_sqrt=1], [have_sqrt=0])
AC_DEFINE_UNQUOTED([HAVE_SQRT], [$have_sqrt],
    [Define to 1 if you have the 'sqrt' function.])
AC_SEARCH_LIBS([hypot], [m], [have_hypot=1], [have_hypot=0])
AC_DEFINE_UNQUOTED([HAVE_HYPOT], [$have_hypot],
    [Define to 1 if you have the 'hypot' function.])
GA_CHECK_FUNCS([bcopy])
GA_CHECK_FUNCS([bzero])
GA_CHECK_FUNCS([create_shared_region])
GA_CHECK_FUNCS([ftruncate])
GA_CHECK_FUNCS([gethostbyname])
GA_CHECK_FUNCS([gethostname])
GA_CHECK_FUNCS([getpagesize])
GA_CHECK_FUNCS([gettimeofday])
GA_CHECK_FUNCS([memalign])
GA_CHECK_FUNCS([memcpy])
GA_CHECK_FUNCS([memset])
GA_CHECK_FUNCS([munmap])
GA_CHECK_FUNCS([pause])
GA_CHECK_FUNCS([posix_memalign])
GA_CHECK_FUNCS([putenv])
GA_CHECK_FUNCS([select])
GA_CHECK_FUNCS([setenv])
GA_CHECK_FUNCS([shared_malloc])
GA_CHECK_FUNCS([shmalloc])
GA_CHECK_FUNCS([socket])
GA_CHECK_FUNCS([srand])
GA_CHECK_FUNCS([srandom])
GA_CHECK_FUNCS([strchr])
GA_CHECK_FUNCS([strdup])
GA_CHECK_FUNCS([strerror])
GA_CHECK_FUNCS([strstr])
GA_CHECK_FUNCS([strtol])

# Checks for C libraries.
GA_CHECK_LIBPLOT

GA_C_OPT

CFLAGS="$ga_save_CFLAGS"

#########################################
# C++ compiler
#########################################
AC_MSG_NOTICE
AS_IF([test "x$enable_cxx" = xyes],
    [AC_MSG_NOTICE([C++ compiler])],
    [AC_MSG_NOTICE([C++ compiler (disabled, but some tests still required)])])
AC_MSG_NOTICE
AC_LANG_PUSH([C++])

# GA_PROG_MPICXX is required to silence complaints that C++ source exists
# even if C++ support has been disabled.
ga_save_CXXFLAGS="$CXXFLAGS"
GA_PROG_MPICXX
CXXFLAGS="$ga_save_CXXFLAGS"
AS_IF([test x$with_mpi_wrappers = xyes], [GA_MPI_UNWRAP])

AS_IF([test "x$enable_cxx" = xyes], [

GA_COMPILER_VENDOR
GA_WARN_FLAGS
ga_save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$ga_cv_cxx_warning_flags $CXXFLAGS"

AS_IF([test "x$with_mpi" != xno], [GA_MPICXX_TEST])

# Checks for C++ header files.

# Checks for C++ typedefs, structures, and compiler characteristics.
GA_CXX_NAMESPACES

# Checks for C++ type sizes.

# Checks for C++ library functions.

GA_CXX_OPT

CXXFLAGS="$ga_save_CXXFLAGS"

])

# conditionals must be defined unconditionally
AM_CONDITIONAL([CXX_NAMESPACES], [test x$ga_cv_cxx_namespaces = xyes])

AC_LANG_POP([C++])

#########################################
# Assembler
#########################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Assembler])
AC_MSG_NOTICE

GA_AS
AM_PROG_AS

#########################################
# Fortran compiler
#########################################
AC_MSG_NOTICE
AS_IF([test "x$enable_f77" = xyes],
    [AC_MSG_NOTICE([Fortran compiler])],
    [AC_MSG_NOTICE([Fortran compiler (disabled, but some tests still required)])])
AC_MSG_NOTICE
AC_LANG_PUSH([Fortran 77])

ga_save_FFLAGS="$FFLAGS"
GA_PROG_MPIF77
FFLAGS="$ga_save_FFLAGS"

# If Fortran 77 compiler was not found, we disable support for it.
AS_IF([test "x$F77" = x], [enable_f77=no])

# If Fortran 77 compiler does not work, we disable support for it.
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([],[])], [],
    [AC_MSG_WARN([Fortran compiler does not work, disabling])
     enable_f77=no
     F77=])

GA_F77_DISABLE_RESULTS

# Perform Fortran 77 tests only if enabled.
if test "x$enable_f77" = xyes; then :

AS_IF([test x$with_mpi_wrappers = xyes], [GA_MPI_UNWRAP])

GA_COMPILER_VENDOR
GA_WARN_FLAGS
ga_save_FFLAGS="$FFLAGS"
FFLAGS="$ga_cv_f77_warning_flags $FFLAGS"

# If F77 is enabled, perform a simple link test for MPI.
AS_IF([test "x$with_mpi" != xno],
    [AS_IF([test "x$enable_mpi_tests" = xyes],
        [GA_MPIF77_TEST_LINK],
        [GA_MPIF77_TEST_COMPILE])])

GA_SUPPRESS_MESSAGE
GA_F77_LIBRARY_LDFLAGS
GA_F77_MAYBE_UNDERSCORING
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS
GA_F77_CPP_SYMBOL
GA_F77_FIXED
GA_F77_MISMATCH_TYPE
GA_F77_INTEGER_SIZE

# Checks for Fortran typedefs, structures, and compiler characteristics.

# Checks for Fortran type sizes.
ga_save_FFLAGS="$FFLAGS"; FFLAGS="$FFLAGS $FFLAG_INT"
GA_F77_CHECK_SIZEOF([INTEGER], [$ga_cv_f77_integer_size])
FFLAGS="$ga_save_FFLAGS"
GA_F77_CHECK_SIZEOF([REAL], [4])
GA_F77_CHECK_SIZEOF([DOUBLE PRECISION], [8])

# Checks for Fortran library functions.
GA_CHECK_FUNCS([dtime etime])
GA_F77_FLUSH

# The rest of the Fortran 77 to C compatibility checks.
GA_F2C_NOMAIN
GA_F2C_CMDARGS
GA_F2C_HIDDEN_STRING_LENGTH_CONVENTION
GA_F2C_MATCH_TYPES([INTEGER], [short, int, long, long long])
GA_F2C_MATCH_TYPES([REAL], [float, double, long double])
GA_F2C_MATCH_TYPES([DOUBLE PRECISION], [float, double, long double])

GA_F77_OPT

FFLAGS="$ga_save_FFLAGS"

else

# enable_f77=no version of AC_F77_WRAPPERS
AC_DEFINE([F77_FUNC(name,NAME)], [name @%:@@%:@ _])
AC_DEFINE([F77_FUNC_(name,NAME)],[name @%:@@%:@ _])
# enable_f77=no version of GA_F77_CHECK_SIZEOF
AC_DEFINE_UNQUOTED(AS_TR_CPP(sizeof_f77_INTEGER), [$ac_cv_sizeof_voidp],
    [The size of 'INTEGER' as computed by C's sizeof.])
AC_DEFINE_UNQUOTED(AS_TR_CPP(sizeof_f77_REAL), [4],
    [The size of 'REAL' as computed by C's sizeof.])
AC_DEFINE_UNQUOTED(AS_TR_CPP(sizeof_f77_DOUBLE PRECISION), [8],
    [The size of 'DOUBLE PRECISION' as computed by C's sizeof.])
# enable_f77=no version of GA_F2C_NOMAIN
AC_SUBST([FLD_NOMAIN], [])
# enable_f77=no version of GA_F2C_CMDARGS
F2C_GETARG="F2C_GETARG"
F2C_IARGC="F2C_IARGC"
AC_SUBST([F2C_GETARG])
AC_SUBST([F2C_IARGC])
# enable_f77=no version of GA_F2C_HIDDEN_STRING_LENGTH_CONVENTION
AC_DEFINE([F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS], [1],
    [whether the hidden string length comes after all other args])
# enable_f77=no version of GA_F2C_MATCH_TYPES
AC_SUBST([F2C_INTEGER_C_TYPE], [long])
AC_SUBST([F2C_REAL_C_TYPE], [float])
AC_SUBST([F2C_DOUBLE_PRECISION_C_TYPE], [double])

fi

# This is required to override what automake does, if needed.
GA_F77_LD_OVERRIDE

AC_LANG_POP([Fortran 77])

###############################################################################
# Checks for BLAS/LAPACK/SCALAPACK require both C and Fortran features to
# have been detected. Further, they assume the current AC_LANG is C.
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Checks for BLAS,LAPACK,ScaLAPACK])
AC_MSG_NOTICE

GA_BLAS
GA_LAPACK
GA_SCALAPACK

###############################################################################
# Checks for programs.
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Checks for additional programs])
AC_MSG_NOTICE

AC_PROG_GREP
AC_PROG_SED

AM_PATH_PYTHON([],[],[:])
AM_CONDITIONAL([HAVE_PYTHON], [test "x$PYTHON" != x])

AC_PATH_PROG([M4], [m4], [notfound])
AM_CONDITIONAL([HAVE_M4], [test "x$M4" != xnotfound])

AC_PATH_PROG([PERL], [perl], [notfound])
AM_CONDITIONAL([HAVE_PERL], [test "x$PERL" != xnotfound])

# X needed for xregion program, but only if desired.
AC_ARG_ENABLE([xregion],
    [AS_HELP_STRING([--enable-xregion], [build global/X/xregion])],
    [enable_xregion=yes],
    [enable_xregion=no])
AS_IF([test "x$enable_xregion" = xyes],
    [AC_PATH_X
     AS_IF([test "x$x_includes" != x],
        [X_CPPFLAGS="-I$x_includes"],
        [X_CPPFLAGS=])
     AC_SUBST([X_CPPFLAGS])
     AS_IF([test "x$x_libraries" != x],
        [X_LDFLAGS="-L$x_libraries"],
        [X_LDFLAGS=])
     AC_SUBST([X_LDFLAGS])])
AM_CONDITIONAL([ENABLE_XREGION], [test "x$enable_xregion" = xyes])
AM_CONDITIONAL([NO_X], [test "$no_x" = xyes])

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
# Remaining setup -- some tests, some individual components
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Miscellaneous remaining setup])
AC_MSG_NOTICE

# Establish some preprocessor symbols.
ARMCI_SETUP

# Set up MA.
MA_ENABLE_ARMCI_MEM_OPTION
MA_LONG_DOUBLE_TYPEDEF
MA_STATS
MA_VERIFY
AS_IF([test "x$ga_cv_target" = xLINUX64],
    [AS_CASE([$host_cpu],
        [x86_64|ppc64],
            [AC_DEFINE([NOUSE_MMAP], [1], [Set some mallopt options])])])

# Create proper types for our access functions.
AS_IF([test "x$ac_cv_sizeof_voidp" = "x8"],
    [ga_access_index_type="integer*8"
     ga_access_index_type_c="long"],
    [ga_access_index_type="integer"
     ga_access_index_type_c="int"])
AC_SUBST([GA_ACCESS_INDEX_TYPE], [$ga_access_index_type])
AC_SUBST([MA_ACCESS_INDEX_TYPE], [$ga_access_index_type])
AC_SUBST([MA_ACCESS_INDEX_TYPE_C], [$ga_access_index_type_c])

# Make GP arrays optional.
AC_ARG_ENABLE([gparrays],
    [AS_HELP_STRING([--enable-gparrays], [enable Global Pointer Arrays])],
    [],
    [enable_gparrays=no])
AM_CONDITIONAL([ENABLE_GPARRAYS], [test "x$enable_gparrays" = xyes])

###############################################################################
# Test suite setup
###############################################################################
AC_ARG_VAR([NPROCS], [number of procs to use for parallel tests (default 4)])
AS_IF([test "x$NPROCS" = x], [NPROCS=4])
AC_SUBST([NPROCS])
AC_ARG_VAR([MPIEXEC],
           [how to run parallel tests if built with MPI e.g. "mpiexec -np %NP%"])
AS_CASE([$ga_msg_comms],
    [TCGMSG],  [TCGEXEC="`pwd`/armci/parallel.x"],
    [TCGMSG5], [TCGEXEC=],
    [MPI|TCGMSGMPI],
        [AS_IF([test "x$MPIEXEC" = x],
            [AC_PATH_PROGS([MPIEXEC], [mpirun mpiexec])
             MPIEXEC="$MPIEXEC -n %NP%"])])
AC_SUBST([MPIEXEC])
AC_SUBST([TCGEXEC])

###############################################################################
# The End
###############################################################################
AC_CONFIG_FILES([Makefile
                 gaf2c/farg.h
                 gaf2c/typesf2c.h
                 global/src/global.fh
                 ma/mafdecls.fh
                 ma/matypes.h])

AC_CONFIG_FILES([tools/ga-config], [chmod +x tools/ga-config])
AC_CONFIG_SUBDIRS([armci])
AC_CONFIG_SUBDIRS([comex])
AC_OUTPUT

# Report on what we found.
AC_MSG_NOTICE([])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([ $PACKAGE_NAME configured as follows:])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([                TARGET=$ga_cv_target])
AC_MSG_NOTICE([             MSG_COMMS=$ga_msg_comms])
AC_MSG_NOTICE([            GA_MP_LIBS=$GA_MP_LIBS])
AC_MSG_NOTICE([         GA_MP_LDFLAGS=$GA_MP_LDFLAGS])
AC_MSG_NOTICE([        GA_MP_CPPFLAGS=$GA_MP_CPPFLAGS])
AC_MSG_NOTICE([         ARMCI_NETWORK=$ga_armci_network])
AC_MSG_NOTICE([ ARMCI_NETWORK_LDFLAGS=$ARMCI_NETWORK_LDFLAGS])
AC_MSG_NOTICE([    ARMCI_NETWORK_LIBS=$ARMCI_NETWORK_LIBS])
AC_MSG_NOTICE([ARMCI_NETWORK_CPPFLAGS=$ARMCI_NETWORK_CPPFLAGS])
AS_IF([test "x$enable_f77" = xyes], [
AC_MSG_NOTICE([                   F77=$F77])
AS_IF([test "x$with_mpi_wrappers" = xyes], [
AC_MSG_NOTICE([         unwrapped F77=$ga_cv_mpif77_naked])
])
AC_MSG_NOTICE([                FFLAGS=$FFLAGS])
AC_MSG_NOTICE([             FFLAG_INT=$FFLAG_INT])
AC_MSG_NOTICE([               GA_FOPT=$GA_FOPT])
])
AC_MSG_NOTICE([                    CC=$CC])
AS_IF([test "x$with_mpi_wrappers" = xyes], [
AC_MSG_NOTICE([          unwrapped CC=$ga_cv_mpic_naked])
])
AC_MSG_NOTICE([                CFLAGS=$CFLAGS])
AC_MSG_NOTICE([               GA_COPT=$GA_COPT])
AS_IF([test "x$enable_cxx" = xyes], [
AC_MSG_NOTICE([                   CXX=$CXX])
AS_IF([test "x$with_mpi_wrappers" = xyes], [
AC_MSG_NOTICE([         unwrapped CXX=$ga_cv_mpicxx_naked])
])
AC_MSG_NOTICE([              CXXFLAGS=$CXXFLAGS])
AC_MSG_NOTICE([             GA_CXXOPT=$GA_CXXOPT])
])
AC_MSG_NOTICE([                   CPP=$CPP])
AC_MSG_NOTICE([              CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([               LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([                  LIBS=$LIBS])
AC_MSG_NOTICE([                 FLIBS=$FLIBS])
AC_MSG_NOTICE([          BLAS_LDFLAGS=$BLAS_LDFLAGS])
AC_MSG_NOTICE([             BLAS_LIBS=$BLAS_LIBS])
AC_MSG_NOTICE([         BLAS_CPPFLAGS=$BLAS_CPPFLAGS])
AC_MSG_NOTICE([                    AR=$AR])
AC_MSG_NOTICE([              AR_FLAGS=$AR_FLAGS])
AC_MSG_NOTICE([                  CCAS=$CCAS])
AC_MSG_NOTICE([            CCAS_FLAGS=$CCAS_FLAGS])
AC_MSG_NOTICE([                  DEFS=$DEFS])
AC_MSG_NOTICE([                 SHELL=$SHELL])
AS_CASE([$ga_msg_comms], [MPI|TCGMSGMPI], [
AC_MSG_NOTICE([               MPIEXEC=$MPIEXEC])
])
AC_MSG_NOTICE([                NPROCS=$NPROCS])
AC_MSG_NOTICE([])
