#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "globalp.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#if ENABLE_F77
extern void gai_diag_seq_(Integer*, Integer*, Integer*, DoublePrecision*);
extern void gai_diag_std_seq_(Integer*, Integer*, DoublePrecision*);
#   define gai_diag_seq_     F77_FUNC_(gai_diag_seq,GAI_DIAG_SEQ)
#   define gai_diag_std_seq_ F77_FUNC_(gai_diag_std_seq,GAI_DIAG_STD_SEQ)
#endif

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_diag_seq = pnga_diag_seq
#endif
void pnga_diag_seq(Integer g_a, Integer g_s, Integer g_v, 
		       DoublePrecision *eval) {
#if ENABLE_F77
    gai_diag_seq_(&g_a, &g_s, &g_v, eval);
#else
    pnga_error("ga_diag_seq: you must configure --enable-f77", 0L);
#endif
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_diag_std_seq = pnga_diag_std_seq
#endif
void pnga_diag_std_seq(Integer g_a, Integer g_v, 
			   DoublePrecision *eval) {
#if ENABLE_F77
    gai_diag_std_seq_(&g_a, &g_v, eval);
#else
    pnga_error("ga_diag_std_seq: you must configure --enable-f77", 0L);
#endif
}
