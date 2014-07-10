/**
 * Tests the scan_add function in GA.
 *
 * Each test will locally perform the same functionality, then compare
 * local buffers against global buffers.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#define NELEM 200000
#define HEAP 200*200*4
#define FUDGE 100
#define STACK 200*200

#include "ga.h"
#include "macdecls.h"
#include "mp3.h" 
#include <stdlib.h>
#include <string.h>

static int me;
static int nproc;

#define test_scan_copy_reg1 local_src[i] = i+1
#define test_scan_copy_cpl1 local_src[i].real = i+1; local_src[i].imag = i+2
#define test_scan_copy_reg2 local_msk[i] = rand()%q
#define test_scan_copy_cpl2 local_msk[i].real = rand()%q; \
                            local_msk[i].imag = 0
#define test_scan_copy_reg3 local_msk[i] != 0
#define test_scan_copy_cpl3 local_msk[i].real != 0
#define test_scan_copy_reg4 local_dst[i] != buf_dst[i]
#define test_scan_copy_cpl4 local_dst[i].real != buf_dst[i].real || \
                            local_dst[i].imag != buf_dst[i].imag
#define test_scan_copy(MT,T,INNER,MT_MSK,T_MSK,INNER_MSK) \
static void test_scan_copy_##MT##_##MT_MSK(int llo, int lhi, int q) \
{ \
    int g_src, g_dst, g_msk; \
    int ndim = 1; \
    int dims[] = {NELEM}; \
    int lo[1]; \
    int hi[1]; \
    int alo[] = {0}; \
    int ahi[] = {NELEM-1}; \
    int i; \
    T *local_src, *local_dst, *buf_dst, last_val; \
    T_MSK *local_msk; \
 \
    lo[0] = llo; \
    hi[0] = lhi; \
 \
    g_src = NGA_Create(MT,    ndim, dims, "g_src", NULL); \
    g_dst = NGA_Create(MT,    ndim, dims, "g_dst", NULL); \
    g_msk = NGA_Create(MT_MSK,ndim, dims, "g_msk", NULL); \
    local_src = malloc(sizeof(T)*NELEM); \
    local_dst = malloc(sizeof(T)*NELEM); \
    local_msk = malloc(sizeof(T_MSK)*NELEM); \
    buf_dst   = malloc(sizeof(T)*NELEM); \
 \
    (void)memset(&last_val, 0, sizeof(T)); \
    (void)memset(local_src, 0, sizeof(T)*NELEM); \
    (void)memset(local_dst, 0, sizeof(T)*NELEM); \
    (void)memset(local_msk, 0, sizeof(T_MSK)*NELEM); \
    (void)memset(buf_dst,   0, sizeof(T)*NELEM); \
 \
    /* process 0 initializes all buffers */ \
    if (0 == me) { \
        for (i=0; i<NELEM; i++) { \
            test_scan_copy_##INNER##1; \
            test_scan_copy_##INNER_MSK##2; \
        } \
        NGA_Put(g_src, alo, ahi, local_src, NULL); \
        NGA_Put(g_msk, alo, ahi, local_msk, NULL); \
    } \
    GA_Zero(g_dst); \
    GA_Sync(); \
    if (0 != me) { \
        NGA_Get(g_src, alo, ahi, local_src, NULL); \
        NGA_Get(g_msk, alo, ahi, local_msk, NULL); \
    } \
    GA_Sync(); \
 \
    /* perform local scan copy on all procs */ \
    for (i=llo; i<=lhi; i++) { \
        if (test_scan_copy_##INNER_MSK##3) { \
            last_val = local_src[i]; \
        } \
        local_dst[i] = last_val; \
    } \
 \
    /* perform global scan copy, get result, compare result */ \
    GA_Scan_copy(g_src, g_dst, g_msk, *lo, *hi); \
    GA_Sync(); \
    NGA_Get(g_dst, alo, ahi, buf_dst, NULL); \
    GA_Sync(); \
    for (i=0; i<NELEM; i++) { \
        if (test_scan_copy_##INNER##4) { \
            GA_Error("scan_copy mismatch in " #MT, i); \
        } \
    } \
 \
    free(local_src); \
    free(local_dst); \
    free(local_msk); \
    free(buf_dst); \
    GA_Destroy(g_src); \
    GA_Destroy(g_dst); \
    GA_Destroy(g_msk); \
}
test_scan_copy(C_INT,int,reg,           C_INT,int,reg)
test_scan_copy(C_LONG,long,reg,         C_INT,int,reg)
test_scan_copy(C_LONGLONG,long long,reg,C_INT,int,reg)
test_scan_copy(C_FLOAT,float,reg,       C_INT,int,reg)
test_scan_copy(C_DBL,double,reg,        C_INT,int,reg)
test_scan_copy(C_SCPL,SingleComplex,cpl,C_INT,int,reg)
test_scan_copy(C_DCPL,DoubleComplex,cpl,C_INT,int,reg)
test_scan_copy(C_INT,int,reg,           C_LONG,long,reg)
test_scan_copy(C_LONG,long,reg,         C_LONG,long,reg)
test_scan_copy(C_LONGLONG,long long,reg,C_LONG,long,reg)
test_scan_copy(C_FLOAT,float,reg,       C_LONG,long,reg)
test_scan_copy(C_DBL,double,reg,        C_LONG,long,reg)
test_scan_copy(C_SCPL,SingleComplex,cpl,C_LONG,long,reg)
test_scan_copy(C_DCPL,DoubleComplex,cpl,C_LONG,long,reg)
test_scan_copy(C_INT,int,reg,           C_LONGLONG,long long,reg)
test_scan_copy(C_LONG,long,reg,         C_LONGLONG,long long,reg)
test_scan_copy(C_LONGLONG,long long,reg,C_LONGLONG,long long,reg)
test_scan_copy(C_FLOAT,float,reg,       C_LONGLONG,long long,reg)
test_scan_copy(C_DBL,double,reg,        C_LONGLONG,long long,reg)
test_scan_copy(C_SCPL,SingleComplex,cpl,C_LONGLONG,long long,reg)
test_scan_copy(C_DCPL,DoubleComplex,cpl,C_LONGLONG,long long,reg)
test_scan_copy(C_INT,int,reg,           C_FLOAT,float,reg)
test_scan_copy(C_LONG,long,reg,         C_FLOAT,float,reg)
test_scan_copy(C_LONGLONG,long long,reg,C_FLOAT,float,reg)
test_scan_copy(C_FLOAT,float,reg,       C_FLOAT,float,reg)
test_scan_copy(C_DBL,double,reg,        C_FLOAT,float,reg)
test_scan_copy(C_SCPL,SingleComplex,cpl,C_FLOAT,float,reg)
test_scan_copy(C_DCPL,DoubleComplex,cpl,C_FLOAT,float,reg)
test_scan_copy(C_INT,int,reg,           C_DBL,double,reg)
test_scan_copy(C_LONG,long,reg,         C_DBL,double,reg)
test_scan_copy(C_LONGLONG,long long,reg,C_DBL,double,reg)
test_scan_copy(C_FLOAT,float,reg,       C_DBL,double,reg)
test_scan_copy(C_DBL,double,reg,        C_DBL,double,reg)
test_scan_copy(C_SCPL,SingleComplex,cpl,C_DBL,double,reg)
test_scan_copy(C_DCPL,DoubleComplex,cpl,C_DBL,double,reg)
test_scan_copy(C_INT,int,reg,           C_SCPL,SingleComplex,cpl)
test_scan_copy(C_LONG,long,reg,         C_SCPL,SingleComplex,cpl)
test_scan_copy(C_LONGLONG,long long,reg,C_SCPL,SingleComplex,cpl)
test_scan_copy(C_FLOAT,float,reg,       C_SCPL,SingleComplex,cpl)
test_scan_copy(C_DBL,double,reg,        C_SCPL,SingleComplex,cpl)
test_scan_copy(C_SCPL,SingleComplex,cpl,C_SCPL,SingleComplex,cpl)
test_scan_copy(C_DCPL,DoubleComplex,cpl,C_SCPL,SingleComplex,cpl)
test_scan_copy(C_INT,int,reg,           C_DCPL,DoubleComplex,cpl)
test_scan_copy(C_LONG,long,reg,         C_DCPL,DoubleComplex,cpl)
test_scan_copy(C_LONGLONG,long long,reg,C_DCPL,DoubleComplex,cpl)
test_scan_copy(C_FLOAT,float,reg,       C_DCPL,DoubleComplex,cpl)
test_scan_copy(C_DBL,double,reg,        C_DCPL,DoubleComplex,cpl)
test_scan_copy(C_SCPL,SingleComplex,cpl,C_DCPL,DoubleComplex,cpl)
test_scan_copy(C_DCPL,DoubleComplex,cpl,C_DCPL,DoubleComplex,cpl)


int main(int argc, char **argv)
{
    int i=0,lo=0,hi=0,q=0;

    MP_INIT(argc,argv);
    GA_INIT(argc,argv);

    me = GA_Nodeid();
    nproc = GA_Nnodes();
    MA_init(MT_DCPL, STACK, HEAP/nproc + FUDGE);

    if (0 == me) {
        printf("NELEM=%d\n", NELEM);
        fflush(stdout);
    }

    for (q=1; q<=6; q++) {
    for (i=0; i<5; i++) {
        switch (i) {
            case 0: lo=0; hi=NELEM-1; break;
            case 1: lo=1; hi=NELEM-1; break;
            case 2: lo=0; hi=NELEM-2; break;
            case 3: lo=NELEM/4; hi=NELEM/2; break;
            case 4: lo=NELEM/3; hi=NELEM/3+10; break;
            default: GA_Error("oops",1); break;
        }
#define tests(MT,MT_MSK) \
        if (0 == me) { \
            printf("testing lo=%d hi=%d q=%d " #MT "\t" #MT_MSK "\n", \
                    lo, hi, q); \
            fflush(stdout); \
        } \
        test_scan_copy_##MT##_##MT_MSK(lo,hi,q)

        tests(C_INT,     C_INT);
        tests(C_LONG,    C_INT);
        tests(C_LONGLONG,C_INT);
        tests(C_FLOAT,   C_INT);
        tests(C_DBL,     C_INT);
        tests(C_SCPL,    C_INT);
        tests(C_DCPL,    C_INT);

        tests(C_INT,     C_LONG);
        tests(C_LONG,    C_LONG);
        tests(C_LONGLONG,C_LONG);
        tests(C_FLOAT,   C_LONG);
        tests(C_DBL,     C_LONG);
        tests(C_SCPL,    C_LONG);
        tests(C_DCPL,    C_LONG);

        tests(C_INT,     C_LONGLONG);
        tests(C_LONG,    C_LONGLONG);
        tests(C_LONGLONG,C_LONGLONG);
        tests(C_FLOAT,   C_LONGLONG);
        tests(C_DBL,     C_LONGLONG);
        tests(C_SCPL,    C_LONGLONG);
        tests(C_DCPL,    C_LONGLONG);

        tests(C_INT,     C_FLOAT);
        tests(C_LONG,    C_FLOAT);
        tests(C_LONGLONG,C_FLOAT);
        tests(C_FLOAT,   C_FLOAT);
        tests(C_DBL,     C_FLOAT);
        tests(C_SCPL,    C_FLOAT);
        tests(C_DCPL,    C_FLOAT);

        tests(C_INT,     C_DBL);
        tests(C_LONG,    C_DBL);
        tests(C_LONGLONG,C_DBL);
        tests(C_FLOAT,   C_DBL);
        tests(C_DBL,     C_DBL);
        tests(C_SCPL,    C_DBL);
        tests(C_DCPL,    C_DBL);

        tests(C_INT,     C_SCPL);
        tests(C_LONG,    C_SCPL);
        tests(C_LONGLONG,C_SCPL);
        tests(C_FLOAT,   C_SCPL);
        tests(C_DBL,     C_SCPL);
        tests(C_SCPL,    C_SCPL);
        tests(C_DCPL,    C_SCPL);

        tests(C_INT,     C_DCPL);
        tests(C_LONG,    C_DCPL);
        tests(C_LONGLONG,C_DCPL);
        tests(C_FLOAT,   C_DCPL);
        tests(C_DBL,     C_DCPL);
        tests(C_SCPL,    C_DCPL);
        tests(C_DCPL,    C_DCPL);
    }
    }

    GA_Terminate();
    MP_FINALIZE();

    return 0;
}
