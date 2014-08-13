#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "globalp.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#if ENABLE_PEIGS
#   if ENABLE_F77
#       define gai_diag_       F77_FUNC_(gai_diag,GAI_DIAG)
#       define gai_diag_std_   F77_FUNC_(gai_diag_std,GAI_DIAG_STD)
#       define gai_diag_reuse_ F77_FUNC_(gai_diag_reuse,GAI_DIAG_REUSE)
extern gai_diag_(Integer*,Integer*,Integer*,DoublePrecision*);
extern gai_diag_std_(Integer*,Integer*,DoublePrecision*);
extern gai_diag_reuse_(Integer*,Integer*,Integer*,Integer*,DoublePrecision*);
#   else
#   endif
#else
#endif

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_diag = pnga_diag
#endif
void pnga_diag(Integer g_a, Integer g_s, Integer g_v, DoublePrecision *eval) {
#if ENABLE_PEIGS
#   if ENABLE_F77
    gai_diag_(&g_a, &g_s, &g_v, eval);
#   else
    pnga_error("ga_diag:peigs interfaced, need to configure --enable-f77",0L);
#   endif
#else
    pnga_error("ga_diag:peigs not interfaced",0L);
#endif
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_diag_std = pnga_diag_std
#endif
void pnga_diag_std(Integer g_a, Integer g_v, DoublePrecision *eval) {
#if ENABLE_PEIGS
#   if ENABLE_F77
    gai_diag_std_(&g_a, &g_v, eval);
#   else
    pnga_error("ga_diag:peigs interfaced, need to configure --enable-f77",0L);
#   endif
#else
    pnga_error("ga_diag:peigs not interfaced",0L);
#endif
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_diag_reuse = pnga_diag_reuse
#endif
void pnga_diag_reuse(Integer reuse, Integer g_a, Integer g_s, 
		   Integer g_v, DoublePrecision *eval) {
#if ENABLE_PEIGS
#   if ENABLE_F77
    gai_diag_reuse_(&reuse, &g_a, &g_s, &g_v, eval);
#   else
    pnga_error("ga_diag:peigs interfaced, need to configure --enable-f77",0L);
#   endif
#else
    pnga_error("ga_diag:peigs not interfaced",0L);
#endif
}
