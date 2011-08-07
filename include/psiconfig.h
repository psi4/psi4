/* include/psiconfig.h.  Generated from psiconfig.h.in by configure.  */

#ifndef _psi_include_psiconfig_h
#define _psi_include_psiconfig_h

/* Define if you have <stdint.h>.  */
#define HAVE_STDINT_H 1

/* PSI version information */
#define PSI_VERSION "4.0.0-alpha"
#define PSI_BUILDID "alpha"
#define PSI_BUGREPORT "psicode@users.sourceforge.net"

/* Top source code directory */
#define PSI_TOP_SRCDIR "/Users/masaaki/psi4"

/* MPI? */
#define HAVE_MPI 1

/* MADNESS? */
#define HAVE_MADNESS 1

/* Defined if we are using MKL */
/* #undef HAVE_MKL */

/* Defined if we are using an MKL with mkl_malloc (MKL version 10+ supposedly) */
/* #undef HAVE_MKL_MALLOC */

#define HAVE_CMATH 1
#define HAVE_CSTDIO 1
#define HAVE_CSTDLIB 1
#define HAVE_CSTRING 1
#define HAVE_CSTDDEF 1
#define HAVE_DECL_PUTENV 1
#define HAVE_PUTENV 1
#define HAVE_DECL_SETENV 1
#define HAVE_SETENV 1
#define HAVE_FUNC_ISINF 1
#define HAVE_FUNC_ERF 1

/* Have dlfcn.h for dlopen and friends */
#define HAVE_DLFCN_H 1

/* Compiler supports __builtin_expected */
#define HAVE_BUILTIN_EXPECT 1

/* Compiler supports __builtin_prefetch */
#define HAVE_BUILTIN_PREFETCH 1

/* Compiler supports __builtin_constant_p */
#define HAVE_BUILTIN_CONSTANT_P 1

/* Restrict keyword definition */
#define restrict __restrict

#endif /* _psi_src_psiconfig_h */

