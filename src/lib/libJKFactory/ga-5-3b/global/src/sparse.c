#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif

#include "abstract_ops.h"
#include "globalp.h"
#include "macdecls.h"
#include "message.h"
#include "ga-papi.h"
#include "ga-wapi.h"

/*\ sets values for specified array elements by enumerating with stride
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_patch_enum = pnga_patch_enum
#endif
void pnga_patch_enum(Integer g_a, Integer lo, Integer hi, void* start, void* stride)
{
Integer dims[1],lop,hip;
Integer ndim, type, me, off;

   pnga_sync();
   me = pnga_nodeid();

   pnga_check_handle(g_a, "ga_patch_enum");

   ndim = pnga_ndim(g_a);
   if (ndim > 1) pnga_error("ga_patch_enum:applicable to 1-dim arrays",ndim);

   pnga_inquire(g_a, &type, &ndim, dims);
   pnga_distribution(g_a, me, &lop, &hip);

   if ( lop > 0 ){ /* we get 0 if no elements stored on this process */

      /* take product of patch owned and specified by the user */ 
      if(hi <lop || hip <lo); /* we got no elements to update */
      else{
        void *ptr;
        Integer ld;
        register Integer i;
        register Integer nelem;

        if(lop < lo)lop = lo;
        if(hip > hi)hip = hi;
        nelem = hip-lop+1;
        off = lop - lo;
        pnga_access_ptr(g_a, &lop, &hip, &ptr, &ld);
        
        switch (type) {
#define TYPE_CASE(MT,T,AT)                                                  \
            case MT:                                                        \
                {                                                           \
                    T *aptr = (T*)ptr;                                      \
                    T astart = *((T*)start);                                \
                    T astride = *((T*)stride);                              \
                    for (i=0; i<nelem; i++) {                               \
                        T offset;                                           \
                        assign_mul_constant_##AT(offset,off+i,astride);     \
                        assign_add_##AT(aptr[i],astart,offset);             \
                    }                                                       \
                    break;                                                  \
                }
#include "types.xh"
#undef TYPE_CASE
            default: pnga_error("ga_patch_enum:wrong data type ",type);
        }

        pnga_release_update(g_a, &lop, &hip);
      }
   }
   
   pnga_sync();
}



#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_scan_copy = pnga_scan_copy
#endif
void pnga_scan_copy(Integer g_src, Integer g_dst, Integer g_msk,
                           Integer lo, Integer hi)
{       
    long *lim=NULL;
    Integer i, nproc, me, elems, ioff, lop, hip, ndim, dims, ld;
    Integer type_src, type_dst, type_msk, combined_type;
    void *ptr_src=NULL;
    void *ptr_dst=NULL;
    void *ptr_msk=NULL;

    nproc = pnga_nnodes();
    me = pnga_nodeid();

    pnga_check_handle(g_src, "ga_scan_copy 1");
    pnga_check_handle(g_dst, "ga_scan_copy 2");
    pnga_check_handle(g_msk, "ga_scan_copy 3");

    if(!pnga_compare_distr(g_src, g_msk))
        pnga_error("ga_scan_copy: different distribution src",0);
    if(!pnga_compare_distr(g_dst, g_msk))
        pnga_error("ga_scan_copy: different distribution dst",0);

    pnga_inquire(g_src, &type_src, &ndim, &dims);
    pnga_inquire(g_dst, &type_dst, &ndim, &dims);
    pnga_inquire(g_msk, &type_msk, &ndim, &dims);
    if(ndim>1)pnga_error("ga_scan_copy: applicable to 1-dim arrays",ndim);
    if(g_src == g_dst) {
        pnga_error("ga_scan_copy: src and dst must be different arrays", 0);
    }
    if(type_src != type_dst) {
        pnga_error("ga_scan_copy: src and dst arrays must be same type", 0);
    }

    pnga_sync();

    pnga_distribution(g_msk, me, &lop, &hip);

    /* create arrays to hold last bit set on a given process */
    lim = (long *) ga_malloc(nproc, MT_C_LONGINT, "ga scan buf");
    bzero(lim,sizeof(long)*nproc);
    lim[me] = -1;

    /* find last bit set on given process (store as global index) */
    if ( lop > 0 ){ /* we get 0 if no elements stored on this process */ 
        elems = hip - lop + 1;
        pnga_access_ptr(g_msk, &lop, &hip, &ptr_msk, &ld);
        switch (type_msk) {
#define TYPE_CASE(MT,T,AT)                                      \
            case MT:                                            \
                {                                               \
                    T * restrict buf = (T*)ptr_msk;             \
                    for(i=0; i<elems; i++) {                    \
                        if (neq_zero_##AT(buf[i])) {            \
                            ioff = i + lop;                     \
                            if (ioff >= lo && ioff <= hi) {     \
                                lim[me]= ioff;                  \
                            }                                   \
                        }                                       \
                    }                                           \
                    break;                                      \
                }
#include "types.xh"
#undef TYPE_CASE
        }
        pnga_release(g_msk, &lop, &hip);
    }
    pnga_gop(C_LONG,lim,nproc,"+");

    if(hi <lop || hip <lo) {
        /* we have no elements to update */
    }
    else {
        Integer rmt_idx, start, stop;

        if (lop < lo) {
            start = lo - lop;
        } else {
            start = 0;
        }
        if (hip > hi) {
            stop = hi-lop+1;
        } else {
            stop = hip-lop+1;
        }

        /* If start bit is not first local bit, find last remote bit.
         * We must scan the entire lim to find the next-highest index.
         * Otherwise, this algorithm won't work with restricted arrays? */
        rmt_idx = -1;
        for (i=0; i<nproc; i++) {
            if (-1 != lim[i] && lim[i] > rmt_idx && lim[i] < lop) {
                rmt_idx = lim[i];
            }
        }

        pnga_access_ptr(g_src, &lop, &hip, &ptr_src, &ld);
        pnga_access_ptr(g_dst, &lop, &hip, &ptr_dst, &ld);
        pnga_access_ptr(g_msk, &lop, &hip, &ptr_msk, &ld);
        combined_type = MT_NUMTYPES*type_src + type_msk;
        switch (combined_type) {
#define TYPE_CASE(MT,T,AT,MT_MSK,T_MSK,AT_MSK)                              \
            case (MT_NUMTYPES*MT) + MT_MSK:                                 \
                {                                                           \
                    T * restrict src = (T*)ptr_src;                         \
                    T * restrict dst = (T*)ptr_dst;                         \
                    T_MSK * restrict msk = (T_MSK*)ptr_msk;                 \
                    T last_val;                                             \
                    assign_zero_##AT(last_val);                             \
                    if (-1 != rmt_idx) {                                    \
                        pnga_get(g_src, &rmt_idx, &rmt_idx, &last_val, &ld);\
                    }                                                       \
                    for (i=start; i<stop; i++) {                            \
                        if (neq_zero_##AT_MSK(msk[i])) {                    \
                            assign_##AT(last_val, src[i]);                  \
                        }                                                   \
                        assign_##AT(dst[i], last_val);                      \
                    }                                                       \
                    break;                                                  \
                }
#include "types2.xh"
#undef TYPE_CASE
            default: pnga_error("ga_scan_copy:wrong data type",combined_type);
        }
        /* release local access to arrays */
        pnga_release(g_src, &lop, &hip);
        pnga_release(g_msk, &lop, &hip);
        pnga_release_update(g_dst, &lop, &hip);
    }

    pnga_sync();
    ga_free(lim);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_scan_add = pnga_scan_add
#endif
void pnga_scan_add(Integer g_src, Integer g_dst, Integer g_msk,
                           Integer lo, Integer hi, Integer excl)
{
    long *lim=NULL;
    Integer i, nproc, me, elems, ioff, lop, hip, ndim, dims, ld;
    Integer type_src, type_dst, type_msk, combined_type;
    void *ptr_src=NULL;
    void *ptr_dst=NULL;
    void *ptr_msk=NULL;

    nproc = pnga_nnodes();
    me = pnga_nodeid();

    pnga_check_handle(g_src, "ga_scan_add 1");
    pnga_check_handle(g_dst, "ga_scan_add 2");
    pnga_check_handle(g_msk, "ga_scan_add 3");

    if(!pnga_compare_distr(g_src, g_msk))
        pnga_error("ga_scan_add: different distribution src",0);
    if(!pnga_compare_distr(g_dst, g_msk))
        pnga_error("ga_scan_add: different distribution dst",0);

    pnga_inquire(g_src, &type_src, &ndim, &dims);
    pnga_inquire(g_dst, &type_dst, &ndim, &dims);
    pnga_inquire(g_msk, &type_msk, &ndim, &dims);
    if(ndim>1)pnga_error("ga_scan_add: applicable to 1-dim arrays",ndim);
    if(g_src == g_dst) {
        pnga_error("ga_scan_add: src and dst must be different arrays", 0);
    }
    if(type_src != type_dst) {
        pnga_error("ga_scan_add: src and dst arrays must be same type", 0);
    }

    pnga_sync();

    pnga_distribution(g_msk, me, &lop, &hip);

    /* create arrays to hold last bit set on a given process */
    lim = (long *) ga_malloc(nproc, MT_C_LONGINT, "ga scan buf");
    bzero(lim,sizeof(long)*nproc);
    lim[me] = -1;

    /* find last bit set on given process (store as global index) */
    if ( lop > 0 ){ /* we get 0 if no elements stored on this process */ 
        elems = hip - lop + 1;
        pnga_access_ptr(g_msk, &lop, &hip, &ptr_msk, &ld);
        switch (type_msk) {
#define TYPE_CASE(MT,T,AT)                                      \
            case MT:                                            \
                {                                               \
                    T * restrict buf = (T*)ptr_msk;             \
                    for(i=0; i<elems; i++) {                    \
                        if (neq_zero_##AT(buf[i])) {            \
                            ioff = i + lop;                     \
                            if (ioff >= lo && ioff <= hi) {     \
                                lim[me]= ioff;                  \
                            }                                   \
                        }                                       \
                    }                                           \
                    break;                                      \
                }
#include "types.xh"
#undef TYPE_CASE
        }
        pnga_release(g_msk, &lop, &hip);
    }
    pnga_gop(C_LONG,lim,nproc,"+");
#if VERBOSE_DEBUG
    if (0==me) printf("(%d) lo=%d hi=%d\n", me, lo, hi);
    for (i=0; i<nproc; i++) {
        if (i==me) printf("(%d) lop=%d hip=%d\n", me, lop, hip);
        GA_Sync();
    }
#endif

    if(hi <lop || hip <lo) {
        /* we have no elements to update, but there are sync's in the else */
#if VERBOSE_DEBUG
        for (i=0; i<nproc; i++) {
            if (i==me) printf("(%d) rmt_idx=N/A\n", me, lop, hip);
            GA_Sync();
        }
#endif
        pnga_sync();
        pnga_sync();
    }
    else {
        Integer rmt_idx, start, stop;

        /* find the nearest set bit on another process */
        rmt_idx = -1;
        for (i=0; i<nproc; i++) {
            if (-1 != lim[i] && lim[i] > rmt_idx && lim[i] < lop) {
                rmt_idx = lim[i];
            }
        }
#if VERBOSE_DEBUG
        for (i=0; i<nproc; i++) {
            if (i==me) printf("(%d) rmt_idx=%d\n", me, rmt_idx);
            GA_Sync();
        }
#endif

        if (lop < lo) {
            start = lo - lop;
        } else {
            start = 0;
        }
        if (hip > hi) {
            stop = hi-lop+1;
        } else {
            stop = hip-lop+1;
        }

        /* first, perform local scan add */
        pnga_access_ptr(g_src, &lop, &hip, &ptr_src, &ld);
        pnga_access_ptr(g_dst, &lop, &hip, &ptr_dst, &ld);
        pnga_access_ptr(g_msk, &lop, &hip, &ptr_msk, &ld);
        combined_type = MT_NUMTYPES*type_src + type_msk;
        switch (combined_type) {
#define TYPE_CASE(MT,T,AT,MT_MSK,T_MSK,AT_MSK)                              \
            case (MT_NUMTYPES*MT) + MT_MSK:                                 \
                {                                                           \
                    T * restrict src = (T*)ptr_src;                         \
                    T * restrict dst = (T*)ptr_dst;                         \
                    T_MSK * restrict msk = (T_MSK*)ptr_msk;                 \
                    int found_bit = 0;                                      \
                    /* find first set bit on first segment */               \
                    if (lop <= lo) {                                        \
                        while (eq_zero_##AT_MSK(msk[start]) && start<stop) {\
                            start++;                                        \
                        }                                                   \
                    }                                                       \
                    /* set first index then use for subsequent indices */   \
                    if (start < stop) {                                     \
                        i = start;                                          \
                        if (neq_zero_##AT_MSK(msk[i])) {                    \
                            found_bit = 1;                                  \
                            if (excl != 0) {                                \
                                assign_zero_##AT(dst[i]);                   \
                            } else {                                        \
                                assign_##AT(dst[i], src[i]);                \
                            }                                               \
                        } else if (-1 != rmt_idx) {                         \
                            found_bit = 1;                                  \
                            if (excl != 0) {                                \
                                Integer loc = i+lop - 1;                    \
                                if (loc > 0) {                              \
                                    pnga_get(g_src, &loc, &loc, &dst[i], NULL);\
                                }                                           \
                            } else {                                        \
                                assign_##AT(dst[i], src[i]);                \
                            }                                               \
                        }                                                   \
                    }                                                       \
                    if (excl != 0) {                                        \
                        for (i=start+1; i<stop; i++) {                      \
                            if (neq_zero_##AT_MSK(msk[i])) {                \
                                assign_zero_##AT(dst[i]);                   \
                                found_bit = 1;                              \
                            } else if (1 == found_bit) {                    \
                                assign_add_##AT(dst[i], dst[i-1], src[i-1]);\
                            }                                               \
                        }                                                   \
                    } else {                                                \
                        for (i=start+1; i<stop; i++) {                      \
                            if (neq_zero_##AT_MSK(msk[i])) {                \
                                assign_##AT(dst[i], src[i]);                \
                                found_bit = 1;                              \
                            } else if (1 == found_bit) {                    \
                                assign_add_##AT(dst[i], dst[i-1], src[i]);  \
                            }                                               \
                        }                                                   \
                    }                                                       \
                    pnga_sync();                                            \
                    /* lastly, reconcile segment boundaries on other procs */\
                    if (eq_zero_##AT_MSK(msk[start]) && -1 != rmt_idx) {    \
                        Integer np, *map, *proclist, *subs, rmt_hi;         \
                        T *v, sum;                                          \
                        rmt_hi = lop-1;                                     \
                        pnga_locate_nnodes(g_dst, &rmt_idx, &rmt_hi, &np);  \
                        map = ga_malloc(4*np, MT_F_INT, "ga scan add locate");\
                        v = ga_malloc(np, MT, "ga scan add gather values"); \
                        proclist = map+(2*np);                              \
                        subs = map+(3*np);                                  \
                        pnga_locate_region(g_dst, &rmt_idx, &rmt_hi, map, proclist, &np);\
                        for (i=0; i<np; i++) {                              \
                            subs[i] = map[i*2+1];                           \
                        }                                                   \
                        pnga_gather(g_dst, v, subs, 0, np);                 \
                        pnga_sync();                                        \
                        assign_zero_##AT(sum);                              \
                        for (i=0; i<np; i++) {                              \
                            add_assign_##AT(sum, v[i]);                     \
                        }                                                   \
                        for (i=start; i<stop; i++) {                        \
                            if (eq_zero_##AT_MSK(msk[i])) {                 \
                                add_assign_##AT(dst[i], sum);               \
                            } else {                                        \
                                break;                                      \
                            }                                               \
                        }                                                   \
                        ga_free(v);                                         \
                        ga_free(map);                                       \
                    } else {                                                \
                        pnga_sync();                                        \
                    }                                                       \
                    break;                                                  \
                }
#include "types2.xh"
#undef TYPE_CASE
            default: pnga_error("ga_scan_add:wrong data type",combined_type);
        }
        /* release local access to arrays */
        pnga_release(g_src, &lop, &hip);
        pnga_release(g_msk, &lop, &hip);
        pnga_release_update(g_dst, &lop, &hip);
    }

    pnga_sync();
    ga_free(lim);
}


static void sga_pack(Integer first, long lim, Integer elems,
                     Integer type_src, Integer type_msk,
                     void *ptr_src, void *ptr_dst, void *ptr_msk)
{
    Integer combined_type, i, pck_idx=0;
    combined_type = MT_NUMTYPES*type_src + type_msk;
    switch (combined_type) {
#define TYPE_CASE(MT,T,AT,MT_MSK,T_MSK,AT_MSK)                          \
        case (MT_NUMTYPES*MT) + MT_MSK:                                 \
            {                                                           \
                T * restrict pck = (T*)ptr_dst;                         \
                T * restrict src = (T*)ptr_src;                         \
                T_MSK * restrict msk = (T_MSK*)ptr_msk;                 \
                for (i=first; i<elems&&pck_idx<lim; i++) {              \
                    if (neq_zero_##AT_MSK(msk[i])) {                    \
                        assign_##AT(pck[pck_idx], src[i]);              \
                        ++pck_idx;                                      \
                    }                                                   \
                }                                                       \
                break;                                                  \
            }
#include "types2.xh"
#undef TYPE_CASE
    }
}


static void sga_unpack(Integer first, long lim, Integer elems,
                       Integer type_src, Integer type_msk,
                       void *ptr_src, void *ptr_dst, void *ptr_msk)
{
    Integer combined_type, i, pck_idx=0;
    combined_type = MT_NUMTYPES*type_src + type_msk;
    switch (combined_type) {
#define ga_unpack_case(MT,T,AT,MT_MSK,T_MSK,AT_MSK)                     \
        case (MT_NUMTYPES*MT) + MT_MSK:                                 \
            {                                                           \
                T * restrict pck = (T*)ptr_src;                         \
                T * restrict dst = (T*)ptr_dst;                         \
                T_MSK * restrict msk = (T_MSK*)ptr_msk;                 \
                for (i=first; i<elems&&pck_idx<lim; i++) {              \
                    if (neq_zero_##AT_MSK(msk[i])) {                    \
                        assign_##AT(dst[i], pck[pck_idx]);              \
                        ++pck_idx;                                      \
                    }                                                   \
                }                                                       \
                break;                                                  \
            }
        ga_unpack_case(C_INT,int,reg,           C_INT,int,reg)
        ga_unpack_case(C_LONG,long,reg,         C_INT,int,reg)
        ga_unpack_case(C_LONGLONG,long long,reg,C_INT,int,reg)
        ga_unpack_case(C_FLOAT,float,reg,       C_INT,int,reg)
        ga_unpack_case(C_DBL,double,reg,        C_INT,int,reg)
        ga_unpack_case(C_SCPL,SingleComplex,cpl,C_INT,int,reg)
        ga_unpack_case(C_DCPL,DoubleComplex,cpl,C_INT,int,reg)
        ga_unpack_case(C_INT,int,reg,           C_LONG,long,reg)
        ga_unpack_case(C_LONG,long,reg,         C_LONG,long,reg)
        ga_unpack_case(C_LONGLONG,long long,reg,C_LONG,long,reg)
        ga_unpack_case(C_FLOAT,float,reg,       C_LONG,long,reg)
        ga_unpack_case(C_DBL,double,reg,        C_LONG,long,reg)
        ga_unpack_case(C_SCPL,SingleComplex,cpl,C_LONG,long,reg)
        ga_unpack_case(C_DCPL,DoubleComplex,cpl,C_LONG,long,reg)
        ga_unpack_case(C_INT,int,reg,           C_LONGLONG,long long,reg)
        ga_unpack_case(C_LONG,long,reg,         C_LONGLONG,long long,reg)
        ga_unpack_case(C_LONGLONG,long long,reg,C_LONGLONG,long long,reg)
        ga_unpack_case(C_FLOAT,float,reg,       C_LONGLONG,long long,reg)
        ga_unpack_case(C_DBL,double,reg,        C_LONGLONG,long long,reg)
        ga_unpack_case(C_SCPL,SingleComplex,cpl,C_LONGLONG,long long,reg)
        ga_unpack_case(C_DCPL,DoubleComplex,cpl,C_LONGLONG,long long,reg)
        ga_unpack_case(C_INT,int,reg,           C_FLOAT,float,reg)
        ga_unpack_case(C_LONG,long,reg,         C_FLOAT,float,reg)
        ga_unpack_case(C_LONGLONG,long long,reg,C_FLOAT,float,reg)
        ga_unpack_case(C_FLOAT,float,reg,       C_FLOAT,float,reg)
        ga_unpack_case(C_DBL,double,reg,        C_FLOAT,float,reg)
        ga_unpack_case(C_SCPL,SingleComplex,cpl,C_FLOAT,float,reg)
        ga_unpack_case(C_DCPL,DoubleComplex,cpl,C_FLOAT,float,reg)
        ga_unpack_case(C_INT,int,reg,           C_DBL,double,reg)
        ga_unpack_case(C_LONG,long,reg,         C_DBL,double,reg)
        ga_unpack_case(C_LONGLONG,long long,reg,C_DBL,double,reg)
        ga_unpack_case(C_FLOAT,float,reg,       C_DBL,double,reg)
        ga_unpack_case(C_DBL,double,reg,        C_DBL,double,reg)
        ga_unpack_case(C_SCPL,SingleComplex,cpl,C_DBL,double,reg)
        ga_unpack_case(C_DCPL,DoubleComplex,cpl,C_DBL,double,reg)
        ga_unpack_case(C_INT,int,reg,           C_SCPL,SingleComplex,cpl)
        ga_unpack_case(C_LONG,long,reg,         C_SCPL,SingleComplex,cpl)
        ga_unpack_case(C_LONGLONG,long long,reg,C_SCPL,SingleComplex,cpl)
        ga_unpack_case(C_FLOAT,float,reg,       C_SCPL,SingleComplex,cpl)
        ga_unpack_case(C_DBL,double,reg,        C_SCPL,SingleComplex,cpl)
        ga_unpack_case(C_SCPL,SingleComplex,cpl,C_SCPL,SingleComplex,cpl)
        ga_unpack_case(C_DCPL,DoubleComplex,cpl,C_SCPL,SingleComplex,cpl)
        ga_unpack_case(C_INT,int,reg,           C_DCPL,DoubleComplex,cpl)
        ga_unpack_case(C_LONG,long,reg,         C_DCPL,DoubleComplex,cpl)
        ga_unpack_case(C_LONGLONG,long long,reg,C_DCPL,DoubleComplex,cpl)
        ga_unpack_case(C_FLOAT,float,reg,       C_DCPL,DoubleComplex,cpl)
        ga_unpack_case(C_DBL,double,reg,        C_DCPL,DoubleComplex,cpl)
        ga_unpack_case(C_SCPL,SingleComplex,cpl,C_DCPL,DoubleComplex,cpl)
        ga_unpack_case(C_DCPL,DoubleComplex,cpl,C_DCPL,DoubleComplex,cpl)
    }
#undef ga_unpack_case
}


static void sga_pack_unpack(Integer g_src, Integer g_dst, Integer g_msk,
                            Integer lo, Integer hi, Integer* icount,
                            Integer pack)
{
    void *ptr;
    long *lim=NULL;
    Integer nproc, me, i, myplace, np=0, first=-1;
    Integer lop, hip, dims;
    Integer ndim_src, ndim_dst, ndim_msk;
    Integer type_src, type_dst, type_msk;

    nproc = pnga_nnodes();
    me = pnga_nodeid();

    pnga_check_handle(g_src, "ga_pack src");
    pnga_check_handle(g_dst, "ga_pack dst");
    pnga_check_handle(g_msk, "ga_pack msk");
    pnga_inquire(g_src, &type_src, &ndim_src, &dims);
    pnga_inquire(g_dst, &type_dst, &ndim_dst, &dims);
    pnga_inquire(g_msk, &type_msk, &ndim_msk, &dims);
    if (1 != ndim_src) {
        pnga_error("ga_pack: supports 1-dim arrays only: src", ndim_src);
    }
    if (1 != ndim_dst) {
        pnga_error("ga_pack: supports 1-dim arrays only: dst", ndim_dst);
    }
    if (1 != ndim_msk) {
        pnga_error("ga_pack: supports 1-dim arrays only: msk", ndim_msk);
    }
    if (type_src != type_dst) {
        pnga_error("ga_pack: src and dst must be same type", 0);
    }
    if (1 == pack) {
        if(!pnga_compare_distr(g_src, g_msk)) {
            pnga_error("ga_pack: src and msk distributions differ",0);
        }
    } else if (0 == pack) {
        if(!pnga_compare_distr(g_dst, g_msk)) {
            pnga_error("ga_unpack: dst and msk distributions differ",0);
        }
    } else {
        pnga_error("ga_pack/unpack bad pack flag",0);
    }

    pnga_sync();

    lim = (long*) ga_malloc(nproc, MT_C_LONGINT, "ga_pack lim buf");
    bzero(lim,sizeof(long)*nproc);
    pnga_distribution(g_msk, me, &lop, &hip);

    /* how many elements do we have to copy? */
    if ( lop > 0 ){ /* we get 0 if no elements stored on this process */
        Integer lop_1 = lop-1;
        /* adjust the range of elements to be within <lo,hi> */
        if(lop < lo) lop = lo;
        if(hip > hi) hip = hi;
        if(hi <lop || hip <lo); /* we have no elements to update */
        else{
            Integer elems, ONE=1;
            pnga_locate_nnodes(g_msk, &ONE, &lop_1, &np);
            pnga_access_ptr(g_msk, &lop, &hip, &ptr, &elems);
            elems = hip-lop+1;
            switch (type_msk) {
#define TYPE_CASE(MT,T,AT) case MT:                                         \
                {                                                           \
                    T * restrict aptr = (T*)ptr;                            \
                    for (i=0; i<elems; i++) {                               \
                        if (neq_zero_##AT(aptr[i])) {                       \
                            if (first<0) first=i;                           \
                            lim[np]++;                                      \
                        }                                                   \
                    }                                                       \
                    break;                                                  \
                }
#include "types.xh"
#undef TYPE_CASE
            }
        }
    }

    /* find number of elements everybody else is contributing */
    pnga_gop(C_LONG,lim,nproc,"+");

    for(i= myplace= *icount= 0; i<nproc; i++){
        if(i<np) myplace += lim[i];
        *icount += lim[i];
    }

    if(hi<lop || hip<lo || lim[np]==0 ); /* we have no elements to update */
    else{
        Integer ignore;
        void *buf=NULL, *msk=NULL;
        buf = ga_malloc(lim[np], type_dst, "ga pack buf");
        pnga_access_ptr(g_msk, &lop, &hip, &msk, &ignore);
        if (1 == pack) {
            void *src=NULL;
            Integer dst_lo=myplace+1, dst_hi=myplace+lim[np];
            pnga_access_ptr(g_src, &lop, &hip, &src, &ignore);
            sga_pack(first, lim[np], hip-lop+1,
                     type_src, type_msk, src, buf, msk);
            pnga_put(g_dst, &dst_lo, &dst_hi, buf, &ignore);
            pnga_release(g_src, &lop, &hip);
        } else if (0 == pack) {
            void *dst=NULL;
            Integer src_lo=myplace+1, src_hi=myplace+lim[np];
            pnga_access_ptr(g_dst, &lop, &hip, &dst, &ignore);
            pnga_get(g_src, &src_lo, &src_hi, buf, &ignore);
            sga_unpack(first, lim[np], hip-lop+1,
                       type_src, type_msk, buf, dst, msk);
            pnga_release_update(g_dst, &lop, &hip);
        }
        ga_free(buf);
    }
    ga_free(lim);
    pnga_sync();
}



#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_pack = pnga_pack
#endif
void pnga_pack(Integer g_src, Integer g_dst, Integer g_msk,
              Integer lo, Integer hi, Integer* icount)
{
     sga_pack_unpack(g_src, g_dst, g_msk, lo, hi, icount, 1);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_unpack = pnga_unpack
#endif
void pnga_unpack(Integer g_src, Integer g_dst, Integer g_msk,
              Integer lo, Integer hi, Integer* icount)
{
     sga_pack_unpack(g_src, g_dst, g_msk, lo, hi, icount, 0);
}



#define NWORK 2000
int workR[NWORK], workL[NWORK];

/*\ compute offset for each of n bins for the given processor to contribute its
 *  elements, number of which for each bin is specified in x
\*/ 
static void sgai_bin_offset(int scope, int *x, int n, int *offset)
{
int root, up, left, right;
int len, lenmes, tag=32100, i, me=armci_msg_me();

    if(!x)pnga_error("sgai_bin_offset: NULL pointer", n);
    if(n>NWORK)pnga_error("sgai_bin_offset: >NWORK", n);
    len = sizeof(int)*n;

    armci_msg_bintree(scope, &root, &up, &left, &right);

    /* up-tree phase: collect number of elements */
    if (left > -1) armci_msg_rcv(tag, workL, len, &lenmes, left);
    if (right > -1) armci_msg_rcv(tag, workR, len, &lenmes, right);

    /* add number of elements in each bin */
    if((right > -1) && left>-1) for(i=0;i<n;i++)workL[i] += workR[i] +x[i];
    else if(left > -1) for(i=0;i<n;i++)workL[i] += x[i];
    else for(i=0;i<n;i++)workL[i] = x[i]; 

    /* now, workL on root contains the number of elements in each bin*/
         
    if (me != root && up!=-1) armci_msg_snd(tag, workL, len, up);

    /* down-tree: compute offset subtracting elements for self and right leaf*/
    if (me != root && up!=-1){
             armci_msg_rcv(tag, workL, len, &lenmes, up);
    }
    for(i=0; i<n; i++) offset[i] = workL[i]-x[i];

    if (right > -1) armci_msg_snd(tag, offset, len, right);
    if (left > -1) {
            /* we saved num elems for right subtree to adjust offset for left*/
            for(i=0; i<n; i++) workR[i] = offset[i] -workR[i]; 
            armci_msg_snd(tag, workR, len, left);
    }
/*    printf("%d:left=%d right=%d up=%d root=%d off=%d\n",me,left, right,up,root,offset[0]);
    fflush(stdout);
*/
}

static 
Integer sgai_match_bin2proc(Integer blo, Integer bhi, Integer plo, Integer phi)
{
int rc=0;
       if(blo == plo) rc=1;
       if(bhi == phi) rc+=2; 
       return rc; /* 1 - first 2-last 3-last+first */
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_create_bin_range = pnga_create_bin_range
#endif
logical pnga_create_bin_range(Integer g_bin, Integer g_cnt, Integer g_off, Integer *g_range)
{
Integer type, ndim, nbin, lobin, hibin, me=pnga_nodeid(),crap;
Integer dims[2], nproc=pnga_nnodes(),chunk[2];
Integer tlo[2], thi[2];

    pnga_inquire(g_bin, &type, &ndim, &nbin);
    if(ndim !=1) pnga_error("ga_bin_index: 1-dim array required",ndim);
    if(type!= C_INT && type!=C_LONG && type!=C_LONGLONG)
       pnga_error("ga_bin_index: not integer type",type);

    chunk[0]=dims[0]=2; dims[1]=nproc; chunk[1]=1;
    if(!pnga_create(MT_F_INT, 2, dims, "bin_proc",chunk,g_range)) return FALSE;

    pnga_distribution(g_off,me, &lobin,&hibin);

    if(lobin>0){ /* enter this block when we have data */
      Integer first_proc, last_proc, p;
      Integer first_off, last_off;
      Integer *myoff, bin;

      /* get offset values stored on my processor to first and last bin */
      pnga_access_ptr(g_off, &lobin, &hibin, &myoff, &crap);
      first_off = myoff[0]; last_off = myoff[hibin-lobin];
/*
      pnga_get(g_off,&lobin,&lobin,&first_off,&lo);
      pnga_get(g_off,&hibin,&hibin,&last_off,&hi);
*/

      /* since offset starts at 0, add 1 to get index to g_bin */
      first_off++; last_off++;

      /* find processors on which these bins are located */
      if(!pnga_locate(g_bin, &first_off, &first_proc))
          pnga_error("ga_bin_sorter: failed to locate region f",first_off);
      if(!pnga_locate(g_bin, &last_off, &last_proc))
          pnga_error("ga_bin_sorter: failed to locate region l",last_off);

      /* inspect range of indices to bin elements stored on these processors */
      for(p=first_proc, bin=lobin; p<= last_proc; p++){
          Integer lo, hi, buf[2], off, cnt; 
          buf[0] =-1; buf[1]=-1;

          pnga_distribution(g_bin,p,&lo,&hi);

          for(/* start from current bin */; bin<= hibin; bin++, myoff++){ 
              Integer blo,bhi,stat;

              blo = *myoff +1;
              if(bin == hibin){
                 pnga_get(g_cnt, &hibin, &hibin, &cnt, &hibin); /* local */
                 bhi = blo + cnt-1; 
              }else
                 bhi = myoff[1]; 

              stat= sgai_match_bin2proc(blo, bhi, lo, hi);

              switch (stat) {
              case 0:  /* bin in a middle */ break;
              case 1:  /* first bin on that processor */
                       buf[0] =bin; break;
              case 2:  /* last bin on that processor */
                       buf[1] =bin; break;
              case 3:  /* first and last bin on that processor */
                       buf[0] =bin; buf[1] =bin; break;
              }

              if(stat>1)break; /* found last bin on that processor */
          }
          
          /* set range of bins on processor p */
          cnt =0; off=1;
          if(buf[0]!=-1){cnt=1; off=0;} 
          if(buf[1]!=-1)cnt++; 
          if(cnt){
                 Integer p1 = p+1;
                 lo = 1+off; hi = lo+cnt-1;
                 tlo[0] = lo;
                 tlo[1] = p1;
                 thi[0] = hi;
                 thi[1] = p1;
                 pnga_put(*g_range, tlo, thi, buf+off, &cnt);
          }
      }
   }
   return TRUE;
}


extern void gai_hsort(Integer *list, int n);


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_bin_sorter = pnga_bin_sorter
#endif
void pnga_bin_sorter(Integer g_bin, Integer g_cnt, Integer g_off)
{
Integer nbin,totbin,type,ndim,lo,hi,me=pnga_nodeid(),crap;
Integer g_range;

    if(FALSE==pnga_create_bin_range(g_bin, g_cnt, g_off, &g_range))
        pnga_error("ga_bin_sorter: failed to create temp bin range array",0); 

    pnga_inquire(g_bin, &type, &ndim, &totbin);
    if(ndim !=1) pnga_error("ga_bin_sorter: 1-dim array required",ndim);
     
    pnga_distribution(g_bin, me, &lo, &hi);
    if (lo > 0 ){ /* we get 0 if no elements stored on this process */
        Integer bin_range[2], rlo[2],rhi[2];
        Integer *bin_cnt, *ptr, i;

        /* get and inspect range of bins stored on current processor */
        rlo[0] = 1; rlo[1]= me+1; rhi[0]=2; rhi[1]=rlo[1];
        pnga_get(g_range, rlo, rhi, bin_range, rhi); /* local */
        nbin = bin_range[1]-bin_range[0]+1;
        if(nbin<1 || nbin> totbin || nbin>(hi-lo+1))
           pnga_error("ga_bin_sorter:bad nbin",nbin);

        /* get count of elements in each bin stored on this task */
        if(!(bin_cnt = (Integer*)malloc(nbin*sizeof(Integer))))
           pnga_error("ga_bin_sorter:memory allocation failed",nbin);
        pnga_get(g_cnt,bin_range,bin_range+1,bin_cnt,&nbin);

        /* get access to local bin elements */
        pnga_access_ptr(g_bin, &lo, &hi, &ptr, &crap);
        
        for(i=0;i<nbin;i++){ 
            int elems =(int) bin_cnt[i];
            gai_hsort(ptr, elems);
            ptr+=elems;
        }
        pnga_release_update(g_bin, &lo, &hi);             
    }

    pnga_sync();
}


/*\ note that subs values must be sorted; bins numbered from 1
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_bin_index = pnga_bin_index
#endif
void pnga_bin_index(Integer g_bin, Integer g_cnt, Integer g_off, 
                   Integer *values, Integer *subs, Integer n, Integer sortit)
{
int i, my_nbin=0;
int *all_bin_contrib, *offset;
Integer type, ndim, nbin;

    pnga_inquire(g_bin, &type, &ndim, &nbin);
    if(ndim !=1) pnga_error("ga_bin_index: 1-dim array required",ndim);
    if(type!= C_INT && type!=C_LONG && type!=C_LONGLONG)
       pnga_error("ga_bin_index: not integer type",type);

    all_bin_contrib = (int*)calloc(nbin,sizeof(int));
    if(!all_bin_contrib)pnga_error("ga_binning:calloc failed",nbin);
    offset = (int*)malloc(nbin*sizeof(int));
    if(!offset)pnga_error("ga_binning:malloc failed",nbin);

    /* count how many elements go to each bin */
    for(i=0; i< n; i++){
       int selected = subs[i];
       if(selected <1 || selected> nbin) pnga_error("wrong bin",selected);

       if(all_bin_contrib[selected-1] ==0) my_nbin++; /* new bin found */
       all_bin_contrib[selected-1]++;
    }

    /* process bins in chunks to match available buffer space */
    for(i=0; i<nbin; i+=NWORK){
        int cnbin = ((i+NWORK)<nbin) ? NWORK: nbin -i;
        sgai_bin_offset(SCOPE_ALL, all_bin_contrib+i, cnbin, offset+i);
    }

    for(i=0; i< n; ){
       Integer lo, hi;
       Integer selected = subs[i];
       int elems = all_bin_contrib[selected-1];

       pnga_get(g_off,&selected,&selected, &lo, &selected);
       lo += offset[selected-1]+1;
       hi = lo + elems -1;
/*
       printf("%d: elems=%d lo=%d sel=%d off=%d contrib=%d nbin=%d\n",pnga_nodeid(), elems, lo, selected,offset[selected-1],all_bin_contrib[0],nbin);
*/
       if(lo > nbin) {
	      printf("Writing off end of bins array: index=%d elems=%d lo=%ld hi=%ld values=%ld nbin=%ld\n",
                i,elems,(long)lo,(long)hi,(long)values+i,(long)nbin);
         break;   
       }else{
          pnga_put(g_bin, &lo, &hi, values+i, &selected); 
       }
       i+=elems;
    }
    
    free(offset);
    free(all_bin_contrib);

    if(sortit)pnga_bin_sorter(g_bin, g_cnt, g_off);
    else pnga_sync();
}

