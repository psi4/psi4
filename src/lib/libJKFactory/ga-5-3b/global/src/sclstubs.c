#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "globalp.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#define _MAX_PROB_SIZE_ 10000 /* 100x100 */

#if ENABLE_F77
#   define gai_lu_solve_alt_ F77_FUNC_(gai_lu_solve_alt,GAI_LU_SOLVE_ALT)
#   define gai_llt_solve_ F77_FUNC_(gai_llt_solve,GAI_LLT_SOLVE)
#   define gai_solve_ F77_FUNC_(gai_solve,GAI_SOLVE)
#   define gai_spd_invert_ F77_FUNC_(gai_spd_invert,GAI_SPD_INVERT)
#endif

#ifdef SCALAPACK_I8
#   if   SIZEOF_SHORT == 8
#       define SL_INT short
#   elif SIZEOF_INT == 8
#       define SL_INT int
#   elif SIZEOF_LONG == 8
#       define SL_INT long
#   elif SIZEOF_LONG_LONG == 8
#       define SL_INT long long
#   else
#       error
#   endif
#else
#   if   SIZEOF_SHORT == 4
#       define SL_INT short
#   elif SIZEOF_INT == 4
#       define SL_INT int
#   elif SIZEOF_LONG == 4
#       define SL_INT long
#   elif SIZEOF_LONG_LONG == 4
#       define SL_INT long long
#   else
#       error
#   endif
#endif

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_lu_solve_alt = pnga_lu_solve_alt
#endif
void pnga_lu_solve_alt(Integer tran, Integer g_a, Integer g_b) {
#if HAVE_SCALAPACK
#   if ENABLE_F77
    gai_lu_solve_alt_(&tran, &g_a, &g_b);
#   else
    pnga_error("ga_lu_solve:scalapack interfaced, need configure --enable-f77",0L);
#   endif
#else
    pnga_error("ga_lu_solve:scalapack not interfaced",0L);
#endif
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_lu_solve = pnga_lu_solve
#endif
void pnga_lu_solve(char *tran, Integer g_a, Integer g_b) {

  Integer dimA1, dimA2, typeA;
  Integer dimB1, dimB2, typeB;
  Integer dimsA[2], dimsB[2], ndim;

  /** check GA info for input arrays */
  pnga_check_handle(g_a, "ga_lu_solve: a");
  pnga_check_handle(g_b, "ga_lu_solve: b");
  pnga_inquire (g_a, &typeA, &ndim, dimsA);
  pnga_inquire (g_b, &typeB, &ndim, dimsB);
  dimA1 = dimsA[0];
  dimA2 = dimsA[1];
  dimB1 = dimsB[0];
  dimB2 = dimsB[1];
  
  if( (dimA1*dimA2 > _MAX_PROB_SIZE_) || (dimB1*dimB2 > _MAX_PROB_SIZE_) )
    pnga_error("ga_lu_solve:Array size too large. Use scalapack for optimum performance. configure --with-scalapack or --with-scalapack-i8 for ga_lu_solve to use Scalapack interface",0L);

  pnga_lu_solve_seq(tran, g_a, g_b);
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_llt_solve = pnga_llt_solve
#endif
Integer pnga_llt_solve(Integer g_a, Integer g_b) {
#if HAVE_SCALAPACK
#   if ENABLE_F77
    return gai_llt_solve_(&g_a, &g_b);
#   else
    pnga_error("ga_lu_solve:scalapack interfaced, need configure --enable-f77",0L);
    return FALSE;
#   endif
#else
    pnga_error("ga_lu_solve:scalapack not interfaced",0L);
    return FALSE;
#endif
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_solve = pnga_solve
#endif
Integer pnga_solve(Integer g_a, Integer g_b) {
#if HAVE_SCALAPACK
#   if ENABLE_F77
    return gai_solve_(&g_a, &g_b);
#   else
    pnga_error("ga_lu_solve:scalapack interfaced, need configure --enable-f77",0L);
    return FALSE;
#   endif
#else
    pnga_error("ga_lu_solve:scalapack not interfaced",0L);
    return FALSE;
#endif
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_spd_invert = pnga_spd_invert
#endif
Integer pnga_spd_invert(Integer g_a) {
#if HAVE_SCALAPACK
#   if ENABLE_F77
    return gai_spd_invert_(&g_a);
#   else
    pnga_error("ga_lu_solve:scalapack interfaced, need configure --enable-f77",0L);
    return FALSE;
#   endif
#else
    pnga_error("ga_lu_solve:scalapack not interfaced",0L);
    return FALSE;
#endif
}

