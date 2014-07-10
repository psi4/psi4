#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdlib.h>

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

#define VERBOSE 0
#define NELEM 200000

#define test(MT,T,F) \
void test_##MT() \
{ \
    int g_a, lo[1], hi[1], ld[1], dims[1], chunk[1], i; \
    T *buf, *get, start, inc; \
    start = 1; \
    inc = 1; \
    lo[0] = 0; \
    hi[0] = NELEM-1; \
    dims[0] = NELEM; \
    chunk[0] = 0; \
    buf = (T*)malloc(sizeof(T) * NELEM); \
    get = (T*)malloc(sizeof(T) * NELEM); \
    g_a = NGA_Create(MT, 1, dims, "g_a", NULL); \
    for (i=0; i<100; i++ ){ \
    GA_Patch_enum(g_a, lo[0], hi[0], &start, &inc); \
    } \
    NGA_Get(g_a, lo, hi, get, ld); \
    for (i=0; i<NELEM; i++) { \
        buf[i] = start + i*inc; \
    } \
    for (i=0; i<NELEM; i++) { \
        if (buf[i] != get[i]) { \
            fprintf(stderr, "at %d value mismatch " F " != " F "\n", \
                    i, buf[i], get[i]); \
            GA_Error("error", i); \
        } \
    } \
    GA_Destroy(g_a); \
    free(buf); \
    free(get); \
}
test(C_INT,int,"%d")
test(C_LONG,long,"%ld")
test(C_LONGLONG,long long,"%lld")
test(C_FLOAT,float,"%f")
test(C_DBL,double,"%f")
#undef test

#define test(MT,T) \
void test_##MT() \
{ \
    int g_a, lo[1], hi[1], ld[1], dims[1], chunk[1], i; \
    T *buf, *get, start, inc; \
    start.real = 1; \
    start.imag = 1; \
    inc.real = 1; \
    inc.imag = 1; \
    lo[0] = 0; \
    hi[0] = NELEM-1; \
    dims[0] = NELEM; \
    chunk[0] = 0; \
    buf = (T*)malloc(sizeof(T) * NELEM); \
    get = (T*)malloc(sizeof(T) * NELEM); \
    g_a = NGA_Create(MT, 1, dims, "g_a", NULL); \
    GA_Patch_enum(g_a, lo[0], hi[0], &start, &inc); \
    NGA_Get(g_a, lo, hi, get, ld); \
    for (i=0; i<NELEM; i++) { \
        buf[i].real = start.real + i*inc.real; \
        buf[i].imag = start.imag + i*inc.imag; \
    } \
    for (i=0; i<NELEM; i++) { \
        if (buf[i].real != get[i].real || buf[i].imag != get[i].imag) { \
            fprintf(stderr, "at %d value mismatch {%f,%f} != {%f,%f}\n", \
                    i, buf[i].real, buf[i].imag, get[i].real, get[i].imag); \
            GA_Error("error", i); \
        } \
    } \
    GA_Destroy(g_a); \
    free(buf); \
    free(get); \
}
test(C_SCPL,SingleComplex)
test(C_DCPL,DoubleComplex)
#undef test


int main(int argc, char **argv)
{
    int me, nproc;
    MP_INIT(argc,argv);
    GA_INIT(argc,argv);

    me = GA_Nodeid();
    nproc = GA_Nnodes();

#if VERBOSE
    if (0 == me) printf("TESTING INT\n");
#endif
    test_C_INT();
#if VERBOSE
    if (0 == me) printf("TESTING LONG\n");
#endif
    test_C_LONG();
#if VERBOSE
    if (0 == me) printf("TESTING LONG LONG\n");
#endif
    test_C_LONGLONG();
#if VERBOSE
    if (0 == me) printf("TESTING FLOAT\n");
#endif
    test_C_FLOAT();
#if VERBOSE
    if (0 == me) printf("TESTING DOUBLE\n");
#endif
    test_C_DBL();
#if VERBOSE
    if (0 == me) printf("TESTING SINGLECOMPLEX\n");
#endif
    test_C_SCPL();
#if VERBOSE
    if (0 == me) printf("TESTING DOUBLECOMPLEX\n");
#endif
    test_C_DCPL();

    GA_Terminate();
    MP_FINALIZE();

    return 0;
}
