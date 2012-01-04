AC_PREREQ(2.59)
# was 2.61 but jaguar is not up to date so backported which
# involved replacing AC_PROG_SED and UINT*_T type checks.
AC_INIT([madness], [0.9], [http://code.google.com/p/m-a-d-n-e-s-s])
AC_CONFIG_SRCDIR([configure.ac])

# This is where autoconf automatically generated files go
AC_CONFIG_AUX_DIR([config])

# This is where local macros will be stored (by us)
AC_CONFIG_MACRO_DIR([config])

# Default options for automake ??
AM_INIT_AUTOMAKE([-Wall -Werror nostdinc foreign])

# We want #defines output into a header file
AC_CONFIG_HEADERS([include/madness_config.h])

# Detect specific systems so can provide their whacky defaults
ACX_CRAYXT
ACX_MAC
ACX_IBMBGP
#AC_AIX

# Set default thread processor-affinity
BIND=${BIND-"-1 -1 -1"}
AC_DEFINE_UNQUOTED([MAD_BIND_DEFAULT], ["$BIND"], [The default binding for threads])


# This sets host, host_cpu, host_vendor, host_os
AC_CANONICAL_HOST
echo "HOST INFORMATION:   $host    $host_cpu    $host_vendor    $host_os"
AC_DEFINE_UNQUOTED([HOST_CPU], ["$host_cpu"], [Defines the host cpu (x86, x86_64, ...).])
AC_DEFINE_UNQUOTED([HOST_SYSTEM], ["$host_os"], [Defines the host os (linux-gnu, ...).])

# Pick out CPUs for which we have special source files (assembler)
use_x86_64_asm=no
use_x86_32_asm=no
case $host_cpu in
     x86_64)   
         AC_DEFINE([X86_64], [1], [If set indicates we are using a 64-bit Intel/AMD cpu.]) 
         use_x86_64_asm=yes
         ;;
     i?86*) 
         AC_DEFINE([X86_32], [1], [If set indicates we are using a 32-bit Intel/AMD cpu.]) 
         use_x86_32_asm=yes
         ;;
esac

# Define AM conditional for use of aforementioned assembler
AM_CONDITIONAL([USE_X86_64_ASM], [test $use_x86_64_asm = yes])
AM_CONDITIONAL([USE_X86_32_ASM], [test $use_x86_32_asm = yes])

# All source will be built with the MPI C++ and C compilers
ACX_MPI

initial_CFLAGS="$CFLAGS"
initial_CXXFLAGS="$CXXFLAGS"

AC_PROG_CC([$MPICC])
AC_PROG_CXX([$MPICXX])
AC_PROG_F77([$MPIF77]) # f90 is a superset of f77
AC_LANG_CPLUSPLUS
echo "Compiler choices after 'AC LANG' ... $CXX and $CC"

CFLAGS="$initial_CFLAGS"
CXXFLAGS="$initial_CXXFLAGS"

ACX_DETECT_CXX

# Set debugging
ACX_ENABLE_DEBUGGING

# Set optimization
ACX_ENABLE_OPTIMIZATION

# Set warning flags by compiler vendor
ACX_ENABLE_WARN
ACX_ENABLE_OPTIMAL

# Check for boost
ACX_WITH_BOOST
ACX_CHECK_TR1

# Get optional external libraries inplace so that building will partially check them
ACX_WITH_LIBUNWIND
ACX_WITH_GOOGLE_PERF
ACX_WITH_GOOGLE_TEST

# Check for basic build functionality
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AM_PROG_AS
#AC_PROG_SED
AC_CHECK_PROG([HAVE_SED],[sed],[YES],[NO])
if test $HAVE_SED = YES; then
   AC_PATH_PROG(SED,sed)
fi

# We need these to build the documentation
AC_CHECK_PROG([HAVE_DOXYGEN],[doxygen],[YES],[NO])
if test $HAVE_DOXYGEN = YES; then
   AC_PATH_PROG(DOXYGEN,doxygen)
fi
AC_CHECK_PROG([HAVE_DOT],[dot],[YES],[NO])
if test $HAVE_DOT = YES; then
   AC_PATH_PROG(DOT,dot)
fi
AC_CHECK_PROG([HAVE_PDFLATEX],[pdflatex],[YES],[NO])
if test $HAVE_PDFLATEX = YES; then
   AC_PATH_PROG(PDFLATEX,pdflatex)
fi
AC_CHECK_PROG([HAVE_LATEX],[latex],[YES],[NO])
if test $HAVE_LATEX = YES; then
   AC_PATH_PROG(LATEX,latex)
fi
AM_CONDITIONAL([CAN_BUILD_DOC], [test $HAVE_DOXYGEN = YES])

AC_HEADER_TIME
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stddef.h stdint.h stdlib.h string.h sys/param.h sys/time.h unistd.h])
AC_CHECK_HEADERS([bits/atomicity.h ext/atomicity.h])
AC_HEADER_STDBOOL
# No longer used
#ACX_GNU_HASHMAP

# Types
AC_TYPE_PID_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_HEADER_TIME
#AC_TYPE_UINT64_T
#AC_TYPE_INT64_T
AC_CHECK_TYPES([uint64_t, int64_t, uint32_t, int32_t, uint16_t, uint8_t])
AC_C_VOLATILE
AC_C_CONST
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])

# Library functions
AC_CHECK_FUNCS([sleep random execv perror gettimeofday memmove memset pow sqrt strchr strdup getenv])
AC_FUNC_FORK
# EFV: when using Intel compiler on OS X malloc is found to be non-GNU-compatible, which breaks gcc's cstdlib
#AC_FUNC_MALLOC
AC_FUNC_ERROR_AT_LINE
ACX_POSIX_MEMALIGN

# Check for EIGEN3
ACX_WITH_EIGEN3

# Checks for Fortran linking conventions and BLAS+LAPACK at the same time
#ACX_FORTRAN_SYMBOLS

# Compiler quirks
ACX_UNQUALIFIED_STATIC_DECL
ACX_STD_ABS

# Check for xterm
AC_CHECK_PROG(XTERM,xterm,xterm,no)
if test $XTERM = 'xterm'; then
  AC_DEFINE(HAVE_XTERM,[1],[If set indicates xterm is in path])
fi

# Check for gdb
AC_CHECK_PROG(GDB,gdb,gdb,no)
if test $GDB = 'gdb'; then
  AC_DEFINE(HAVE_GDB,[1],[If set indicates gdb is in path])
fi

# How to handle MADNESS exceptions.
AC_ARG_ENABLE(madex,
              [AC_HELP_STRING([--enable-madex=arg],
                [Specifies behavior of MADNESS assertions: throw (default), abort, assert, disable])],
              [madex=$enableval],
              [madex=throw])
case $madex in
  abort)
    echo MADNESS assertions will abort
    AC_DEFINE([MADNESS_ASSERTIONS_ABORT], [1], [Set if MADNESS assertions abort])
  ;;
  assert)
    echo MADNESS assertions will assert
    AC_DEFINE([MADNESS_ASSERTIONS_ASSERT], [1], [Set if MADNESS assertions assert])
  ;;
  disable)
    echo MADNESS assertions are disabled
    AC_DEFINE([MADNESS_ASSERTIONS_DISABLE], [1], [Set if MADNESS assertions disabled])
  ;;
  *)
    echo MADNESS assertions will throw
    AC_DEFINE([MADNESS_ASSERTIONS_THROW], [1], [Set if MADNESS assertions throw])
  ;;
esac

# Determine if we should use the instrumented new/delete
AC_ARG_ENABLE([memstats],
              [AC_HELP_STRING([--enable-memstats], [If yes, gather memory statistics (expensive)])],
              [memstats=$enableval],
              [memstats=no])
if test "$memstats" = "yes"; then
  echo 'MADNESS will gather memory statistics (expensive)'
  AC_DEFINE([WORLD_GATHER_MEM_STATS], [1], [Set if MADNESS gathers memory statistics])
fi

AC_ARG_ENABLE([profile], 
              [AC_HELP_STRING([--enable-profile], [Enables profiling])], 
              [echo "Enabling profiling"; AC_DEFINE(WORLD_PROFILE_ENABLE, [1], [Define if should enable profiling])], 
              [])

AC_ARG_ENABLE([tensor-bounds-checking], 
              [AC_HELP_STRING([--enable-tensor-bounds-checking],
                [Enable checking of bounds in tensors ... slow but useful for debugging])], 
              [echo "Enabling tensor bounds checking"; AC_DEFINE(TENSOR_BOUNDS_CHECKING, [1], [Define if should enable bounds checking in tensors])], 
              [])

AC_ARG_ENABLE([tensor-instance-count], 
              [AC_HELP_STRING([--enable-tensor-instance-count],
                [Enable counting of allocated tensors for memory leak detection])], 
              [echo "Enabling tensor instance counting"; AC_DEFINE(TENSOR_INSTANCE_COUNT, [1], [Define if should enable instance counting in tensors])], 
              [])

AC_ARG_ENABLE([spinlocks], 
              [AC_HELP_STRING([--enable-spinlocks],
                [Enables use of spinlocks instead of mutexs (faster unless over subscribing processors)])], 
              [echo "Enabling use of spinlocks"; AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])], 
              [])

AC_ARG_ENABLE([never-spin], 
              [AC_HELP_STRING([--enable-never-spin],
                [Disables use of spinlocks (notably for use inside virtual machines)])], 
              [echo "Disabling use of spinlocks"; AC_DEFINE(NEVER_SPIN, [1], [Define if should use never use spinlocks])], 
              [])

AC_ARG_WITH([papi], 
            [AC_HELP_STRING([--with-papi], [Enables use of PAPI])], 
            [echo "Enabling use of PAPI"; AC_DEFINE(HAVE_PAPI,[1], [Define if have PAPI])], 
            [])

# Capture configuration info for user by compiled code
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_USER],["$USER"],[User that configured the code])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_DATE],["`date`"],[Date of configuration])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_CXX], ["$CXX"],[Configured C++ compiler])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_CXXFLAGS], ["$CXXFLAGS"],[Configured C++ compiler flags])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_HOST], ["$HOSTNAME"],[Configured on this machine])

# Not yet enabled the user to configure this stuff but probably should
AC_DEFINE([WORLD_WATCHDOG],[1],[Enables timeout for dead lock detection])
AC_DEFINE([WATCHDOG_BARK_INTERVAL],[300],[Interval for warning about possible deadlock])
AC_DEFINE([WATCHDOG_TIMEOUT],[1200],[Abort if idle for this long])
AC_DEFINE([MADNESS_FORINT],[long],[Equivalent C type for Fortran integers ... should always be long?])


# Stuff that either does not work or is ill-advised
######AC_DEFINE([WORLD_TAU_TRACE],[1],[Enable tracing using tau])
#####AC_DEFINE([WORLD_SCHED_BACKOFF],[1],[Yield during busy spin])

# Display options
AC_MSG_NOTICE([CPPFLAGS = $CPPFLAGS])
AC_MSG_NOTICE([CFLAGS   = $CFLAGS])
AC_MSG_NOTICE([CXXFLAGS = $CXXFLAGS])
AC_MSG_NOTICE([LDFLAGS  = $LDFLAGS])
AC_MSG_NOTICE([LIBS     = $LIBS])

# These are the files that must be generated/preprocessed ... as you add files
# please maintain the directory tree file order and indenting.
AC_CONFIG_FILES([
 Makefile
    config/Makefile.sample
    include/Makefile
    src/Makefile
    src/lib/Makefile
       src/lib/world/Makefile
       src/lib/muParser/Makefile
       src/lib/misc/Makefile
       src/lib/tinyxml/Makefile
       src/lib/tensor/Makefile
       src/lib/linalg/Makefile
       src/lib/mra/Makefile
    src/apps/Makefile
       src/apps/examples/Makefile
       src/apps/tdse/Makefile
       src/apps/hf/Makefile
       src/apps/exciting/Makefile
       src/apps/moldft/Makefile
       src/apps/nick/Makefile
       src/apps/jacob/Makefile
    doc/Makefile
    doc/doxygen.cfg
])

AC_OUTPUT
