#ifndef _psi_include_psiconfig_h
#define _psi_include_psiconfig_h

/* Define if you have <stdint.h>.  */
#undef HAVE_STDINT_H

/* PSI version information */
#undef PSI_VERSION
#undef PSI_BUILDID
#undef PSI_BUGREPORT

/* Top source code directory */
#cmakedefine PSI_TOP_SRCDIR "@PSI_TOP_SRCDIR@"

/* Which integrals standard is used? */
#undef PSI_INTEGRALS_STANDARD

/* MPI? */
#cmakedefine HAVE_MPI @MPI_FOUND@

/* MADNESS? */
#undef HAVE_MADNESS

/* Defined if we are using MKL */
#undef HAVE_MKL

/* Defined if we are using an MKL with mkl_malloc (MKL version 10+ supposedly) */
#undef HAVE_MKL_MALLOC

#cmakedefine HAVE_CMATH
#cmakedefine HAVE_CSTDIO
#cmakedefine HAVE_CSTDLIB
#cmakedefine HAVE_CSTRING
#cmakedefine HAVE_CSTDDEF
#cmakedefine HAVE_FUNC_ISINF

/* Have dlfcn.h for dlopen and friends */
#cmakedefine HAVE_DLFCN_H

/* Compiler supports __builtin_expected */
#undef HAVE_BUILTIN_EXPECT

/* Compiler supports __builtin_prefetch */
#undef HAVE_BUILTIN_PREFETCH

/* Compiler supports __builtin_constant_p */
#undef HAVE_BUILTIN_CONSTANT_P

#endif /* _psi_src_psiconfig_h */

