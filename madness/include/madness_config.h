/* include/madness_config.h.  Generated from madness_config.h.in by configure.  */
/* include/madness_config.h.in.  Generated from configure.ac by autoheader.  */

/* "Target for tuning mtxmq kernels" */
/* #undef AMD_QUADCORE_TUNE */

/* Fortran-C linking convention lower case (no underscore) */
/* #undef FORTRAN_LINKAGE_LC */

/* Fortran-C linking convention lower case with single underscore */
#define FORTRAN_LINKAGE_LCU 1

/* Fortran-C linking convention lower case with double underscore */
/* #undef FORTRAN_LINKAGE_LCUU */

/* Fortran-C linking convention upper case */
/* #undef FORTRAN_LINKAGE_UC */

/* Fortran-C linking convention upper case with single underscore */
/* #undef FORTRAN_LINKAGE_UCU */

/* Define if AMD math library available - ACML */
/* #undef HAVE_ACML */

/* Define to 1 if you have the <bits/atomicity.h> header file. */
/* #undef HAVE_BITS_ATOMICITY_H */

/* Defined if we are running on a Cray-XT */
/* #undef HAVE_CRAYXT */

/* Define to 1 if you have the `execv' function. */
#define HAVE_EXECV 1

/* Define to 1 if you have the <ext/atomicity.h> header file. */
/* #undef HAVE_EXT_ATOMICITY_H */

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* If set indicates gdb is in path */
#define HAVE_GDB 1

/* Define to 1 if you have the `getenv' function. */
#define HAVE_GETENV 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Defined if we are running on an IBMBGP */
/* #undef HAVE_IBMBGP */

/* Define to 1 if the system has the type `int32_t'. */
#define HAVE_INT32_T 1

/* Define to 1 if the system has the type `int64_t'. */
#define HAVE_INT64_T 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if have std::labs(long) */
/* #undef HAVE_LABS */

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the `memmove' function. */
#define HAVE_MEMMOVE 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define if have PAPI */
/* #undef HAVE_PAPI */

/* Define to 1 if you have the `perror' function. */
#define HAVE_PERROR 1

/* Set if have posix_memalign */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the `random' function. */
#define HAVE_RANDOM 1

/* Define to 1 if you have the `sleep' function. */
#define HAVE_SLEEP 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if stdbool.h conforms to C99. */
/* #undef HAVE_STDBOOL_H */

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define if have std::abs(long) */
#define HAVE_STD_ABS_LONG 1

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if the system has the type `uint16_t'. */
#define HAVE_UINT16_T 1

/* Define to 1 if the system has the type `uint32_t'. */
#define HAVE_UINT32_T 1

/* Define to 1 if the system has the type `uint64_t'. */
#define HAVE_UINT64_T 1

/* Define to 1 if the system has the type `uint8_t'. */
#define HAVE_UINT8_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Set if compiler will instantiate static templates */
#define HAVE_UNQUALIFIED_STATIC_DECL 1

/* Define to 1 if you have the `vfork' function. */
#define HAVE_VFORK 1

/* Define to 1 if you have the <vfork.h> header file. */
/* #undef HAVE_VFORK_H */

/* Define to 1 if `fork' works. */
#define HAVE_WORKING_FORK 1

/* Define to 1 if `vfork' works. */
#define HAVE_WORKING_VFORK 1

/* If set indicates xterm is in path */
#define HAVE_XTERM 1

/* Define to 1 if the system has the type `_Bool'. */
#define HAVE__BOOL 1

/* Defines the host cpu (x86, x86_64, ...). */
#define HOST_CPU "x86_64"

/* Defines the host os (linux-gnu, ...). */
#define HOST_SYSTEM "darwin11.0.0"

/* define if array has fill member function. */
/* #undef MADNESS_ARRAY_HAS_FILL */

/* Set if MADNESS assertions abort */
/* #undef MADNESS_ASSERTIONS_ABORT */

/* Set if MADNESS assertions assert */
/* #undef MADNESS_ASSERTIONS_ASSERT */

/* Set if MADNESS assertions disabled */
/* #undef MADNESS_ASSERTIONS_DISABLE */

/* Set if MADNESS assertions throw */
#define MADNESS_ASSERTIONS_THROW 1

/* Configured C++ compiler */
#define MADNESS_CONFIGURATION_CXX "mpicxx"

/* Configured C++ compiler flags */
#define MADNESS_CONFIGURATION_CXXFLAGS " -O3 -Wall -Wno-strict-aliasing -Wno-deprecated"

/* Date of configuration */
#define MADNESS_CONFIGURATION_DATE "Mon Aug  8 02:39:30 JST 2011"

/* Configured on this machine */
#define MADNESS_CONFIGURATION_HOST "Osiris.local"

/* User that configured the code */
#define MADNESS_CONFIGURATION_USER "masaaki"

/* Equivalent C type for Fortran integers ... should always be long? */
#define MADNESS_FORINT long

/* Madness has Boost make_shared and allocate_shared available. */
/* #undef MADNESS_HAS_BOOST_MAKE_SHARED */

/* Madness will use Boost.TR1 where the compiler lacks support for TR1. */
/* #undef MADNESS_HAS_BOOST_TR1 */

/* Madness will use Eigen3 for linear algebra operations */
/* #undef MADNESS_HAS_EIGEN3 */

/* Define if using Google PerformanceTools */
/* #undef MADNESS_HAS_GOOGLE_PERF */

/* Define if using Google PerformanceTools without libunwind */
/* #undef MADNESS_HAS_GOOGLE_PERF_MINIMAL */

/* Define if should use Google unit testing */
/* #undef MADNESS_HAS_GOOGLE_TEST */

/* Define if should use libunwind for Google performance tools */
/* #undef MADNESS_HAS_LIBUNWIND */

/* define if std::tr1::array is available. */
/* #undef MADNESS_HAS_STD_ARRAY */

/* define if std::tr1::hash is available. */
/* #undef MADNESS_HAS_STD_HASH */

/* define if std::make_shared and std::allocate_shared are available. */
/* #undef MADNESS_HAS_STD_MAKE_SHARED */

/* define if std::tr1::shared_ptr is available. */
/* #undef MADNESS_HAS_STD_SHARED_PTR */

/* define if std::tr1::array is available. */
#define MADNESS_HAS_STD_TR1_ARRAY 1

/* define if std::tr1::hash is available. */
#define MADNESS_HAS_STD_TR1_HASH 1

/* define if std::tr1::shared_ptr is available. */
#define MADNESS_HAS_STD_TR1_SHARED_PTR 1

/* define if std::tr1 type traits are available. */
#define MADNESS_HAS_STD_TR1_TYPE_TRAITS 1

/* define if std::tr1 type traits are available. */
/* #undef MADNESS_HAS_STD_TYPE_TRAITS */

/* define if MADNESS is using <array>. */
/* #undef MADNESS_USE_ARRAY */

/* define if MADNESS is using <boost/tr1/array.hpp>. */
/* #undef MADNESS_USE_BOOST_TR1_ARRAY_HPP */

/* define if MADNESS is using <boost/tr1/functional.hpp>. */
/* #undef MADNESS_USE_BOOST_TR1_FUNCTIONAL_HPP */

/* define if MADNESS is using <boost/tr1/memory.hpp>. */
/* #undef MADNESS_USE_BOOST_TR1_MEMORY_HPP */

/* define if MADNESS is using <boost/tr1/type_traits.hpp>. */
/* #undef MADNESS_USE_BOOST_TR1_TYPE_TRAITS_HPP */

/* define if MADNESS is using <functional>. */
/* #undef MADNESS_USE_FUNCTIONAL */

/* define if MADNESS is using <memory>. */
/* #undef MADNESS_USE_MEMORY */

/* define if MADNESS is using <tr1/array>. */
#define MADNESS_USE_TR1_ARRAY 1

/* define if MADNESS is using <tr1/functional>. */
#define MADNESS_USE_TR1_FUNCTIONAL 1

/* define if MADNESS is using <tr1/memory>. */
#define MADNESS_USE_TR1_MEMORY 1

/* define if MADNESS is using <tr1/type_traits>. */
#define MADNESS_USE_TR1_TYPE_TRAITS 1

/* define if MADNESS is using <type_traits>. */
/* #undef MADNESS_USE_TYPE_TRAITS */

/* The default binding for threads */
#define MAD_BIND_DEFAULT "-1 -1 -1"

/* Set if the posix_memalign prototype is missing */
/* #undef MISSING_POSIX_MEMALIGN_PROTO */

/* Define if should use never use spinlocks */
/* #undef NEVER_SPIN */

/* Set if building on a mac */
#define ON_A_MAC 1

/* Name of package */
#define PACKAGE "madness"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "http://code.google.com/p/m-a-d-n-e-s-s"

/* Define to the full name of this package. */
#define PACKAGE_NAME "madness"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "madness 0.9"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "madness"

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.9"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define if should enable bounds checking in tensors */
/* #undef TENSOR_BOUNDS_CHECKING */

/* Define if should enable instance counting in tensors */
/* #undef TENSOR_INSTANCE_COUNT */

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Define if should use spinlocks */
/* #undef USE_SPINLOCKS */

/* Version number of package */
#define VERSION "0.9"

/* Interval for warning about possible deadlock */
#define WATCHDOG_BARK_INTERVAL 300

/* Abort if idle for this long */
#define WATCHDOG_TIMEOUT 1200

/* Set if MADNESS gathers memory statistics */
/* #undef WORLD_GATHER_MEM_STATS */

/* Define if should enable profiling */
/* #undef WORLD_PROFILE_ENABLE */

/* Enables timeout for dead lock detection */
#define WORLD_WATCHDOG 1

/* If set indicates we are using a 32-bit Intel/AMD cpu. */
/* #undef X86_32 */

/* If set indicates we are using a 64-bit Intel/AMD cpu. */
#define X86_64 1

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to equivalent of C99 restrict keyword, or to nothing if this is not
   supported. Do not define if restrict is supported directly. */
#define restrict __restrict

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define as `fork' if `vfork' does not work. */
/* #undef vfork */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */
