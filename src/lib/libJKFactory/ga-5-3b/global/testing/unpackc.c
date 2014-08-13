/**
 * Tests the unpack function in GA.
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

#include <stdlib.h>
#include <string.h>

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

static int me;
static int nproc;

#define assign_reg(a,b)       (a) = (b)
#define assign_cpl(a,b)       (a).real = (b).real; (a).imag = (b).imag
#define assign_val_reg(a,b,c) (a) = (b)
#define assign_val_cpl(a,b,c) (a).real = (b); (a).imag = (c)
#define eq_zero_reg(a)        (0 == (a))
#define eq_zero_cpl(a)        (0 == (a).real && 0 == (a).imag)
#define neq_zero_reg(a)       (0 != (a))
#define neq_zero_cpl(a)       (0 != (a).real || 0 != (a).imag)
#define neq_reg(a,b)          ((a) != (b))
#define neq_cpl(a,b)          ((a).real != (b).real || (a).imag != (b).imag)
#define assign_add_reg(a,b,c) (a) = (b) + (c)
#define assign_add_cpl(a,b,c) (a).real = (b).real + (c).real; \
                              (a).imag = (b).imag + (c).imag
#define cast_reg(a)           ((float)(a))
#define cast_cpl(a)           ((float)(a).real)

#if 0
#   define PRINT(AT,AT_MSK) do {                                              \
    int p;                                                                    \
    for (p=0; p<nproc; p++) {                                                 \
        if (p == me) {                                                        \
            printf("================= %3d ==========================\n", me); \
            for (i=0; i<NELEM; i++) {                                         \
                printf("[%3d] %6.0f %6.0f %6.0f %6.0f", i,                    \
                        cast_##AT(local_src[i]), cast_##AT_MSK(local_msk[i]), \
                        cast_##AT(local_dst[i]), cast_##AT(buf_dst[i]));      \
                if (neq_##AT(local_dst[i], buf_dst[i])) {                     \
                    printf("<--------------------------\n");                  \
                } else {                                                      \
                    printf("\n");                                             \
                }                                                             \
            }                                                                 \
        }                                                                     \
        GA_Sync();                                                            \
    }                                                                         \
} while (0)
#else
#   define PRINT(AT,AT_MSK) do { } while(0)
#endif

#define test_unpack(MT,T,AT,MT_MSK,T_MSK,AT_MSK)                              \
static void test_unpack_##MT##_##MT_MSK(int llo, int lhi, int q)              \
{                                                                           \
    int g_src, g_dst, g_msk;                                                \
    int ndim = 1;                                                           \
    int dims[] = {NELEM};                                                   \
    int alo[] = {0};                                                        \
    int ahi[] = {NELEM-1};                                                  \
    int i, local_count=0, remote_count=0;                                   \
    T *local_src=NULL, *local_dst=NULL, *buf_dst=NULL;                      \
    T_MSK *local_msk=NULL;                                                  \
                                                                            \
    g_src = NGA_Create(MT,     ndim, dims, "g_src", NULL);                  \
    g_dst = NGA_Create(MT,     ndim, dims, "g_dst", NULL);                  \
    g_msk = NGA_Create(MT_MSK, ndim, dims, "g_msk", NULL);                  \
    buf_dst   = malloc(sizeof(T)*NELEM);                                    \
    local_src = malloc(sizeof(T)*NELEM);                                    \
    local_dst = malloc(sizeof(T)*NELEM);                                    \
    local_msk = malloc(sizeof(T_MSK)*NELEM);                                \
                                                                            \
    (void)memset(buf_dst,   0, sizeof(T)*NELEM);                            \
    (void)memset(local_src, 0, sizeof(T)*NELEM);                            \
    (void)memset(local_dst, 0, sizeof(T)*NELEM);                            \
    (void)memset(local_msk, 0, sizeof(T_MSK)*NELEM);                        \
                                                                            \
    /* process 0 initializes all buffers */                                 \
                                                                            \
    if (0 == me) {                                                          \
        for (i=0; i<NELEM; i++) {                                           \
            assign_val_##AT(local_src[i], i+1, 0);                          \
            assign_val_##AT_MSK(local_msk[i], rand()%q==0?1:0, 0);          \
        }                                                                   \
        NGA_Put(g_src, alo, ahi, local_src, NULL);                          \
        NGA_Put(g_msk, alo, ahi, local_msk, NULL);                          \
    }                                                                       \
    GA_Zero(g_dst);                                                         \
    GA_Sync();                                                              \
    if (0 != me) {                                                          \
        NGA_Get(g_src, alo, ahi, local_src, NULL);                          \
        NGA_Get(g_msk, alo, ahi, local_msk, NULL);                          \
    }                                                                       \
    GA_Sync();                                                              \
                                                                            \
    /* preform local unpack on all procs */                                 \
                                                                            \
    for (i=llo; i<=lhi; i++) {                                              \
        if (neq_zero_##AT_MSK(local_msk[i])) {                              \
            assign_##AT(local_dst[i], local_src[local_count]);              \
            ++local_count;                                                  \
        }                                                                   \
    }                                                                       \
                                                                            \
    /* perform global unpack, get result, compare result */                 \
                                                                            \
    GA_Unpack(g_src, g_dst, g_msk, llo, lhi, &remote_count);                \
    GA_Sync();                                                              \
    NGA_Get(g_dst, alo, ahi, buf_dst, NULL);                                \
    GA_Sync();                                                              \
    PRINT(AT,AT_MSK);                                                       \
    for (i=0; i<NELEM; i++) {                                               \
        if (neq_##AT(local_dst[i], buf_dst[i])) {                           \
            GA_Error("unpack mismatch in " #MT "/" #MT_MSK, i);             \
        }                                                                   \
    }                                                                       \
    if (local_count != remote_count) {                                      \
        GA_Error("unpack count mismatch in " #MT "/" #MT_MSK, i);           \
    }                                                                       \
                                                                            \
    free(local_src);                                                        \
    free(local_dst);                                                        \
    free(local_msk);                                                        \
    free(buf_dst);                                                          \
    GA_Destroy(g_src);                                                      \
    GA_Destroy(g_dst);                                                      \
    GA_Destroy(g_msk);                                                      \
}
test_unpack(C_INT,int,reg,           C_INT,int,reg)
test_unpack(C_LONG,long,reg,         C_INT,int,reg)
test_unpack(C_LONGLONG,long long,reg,C_INT,int,reg)
test_unpack(C_FLOAT,float,reg,       C_INT,int,reg)
test_unpack(C_DBL,double,reg,        C_INT,int,reg)
test_unpack(C_SCPL,SingleComplex,cpl,C_INT,int,reg)
test_unpack(C_DCPL,DoubleComplex,cpl,C_INT,int,reg)
test_unpack(C_INT,int,reg,           C_LONG,long,reg)
test_unpack(C_LONG,long,reg,         C_LONG,long,reg)
test_unpack(C_LONGLONG,long long,reg,C_LONG,long,reg)
test_unpack(C_FLOAT,float,reg,       C_LONG,long,reg)
test_unpack(C_DBL,double,reg,        C_LONG,long,reg)
test_unpack(C_SCPL,SingleComplex,cpl,C_LONG,long,reg)
test_unpack(C_DCPL,DoubleComplex,cpl,C_LONG,long,reg)
test_unpack(C_INT,int,reg,           C_LONGLONG,long long,reg)
test_unpack(C_LONG,long,reg,         C_LONGLONG,long long,reg)
test_unpack(C_LONGLONG,long long,reg,C_LONGLONG,long long,reg)
test_unpack(C_FLOAT,float,reg,       C_LONGLONG,long long,reg)
test_unpack(C_DBL,double,reg,        C_LONGLONG,long long,reg)
test_unpack(C_SCPL,SingleComplex,cpl,C_LONGLONG,long long,reg)
test_unpack(C_DCPL,DoubleComplex,cpl,C_LONGLONG,long long,reg)
test_unpack(C_INT,int,reg,           C_FLOAT,float,reg)
test_unpack(C_LONG,long,reg,         C_FLOAT,float,reg)
test_unpack(C_LONGLONG,long long,reg,C_FLOAT,float,reg)
test_unpack(C_FLOAT,float,reg,       C_FLOAT,float,reg)
test_unpack(C_DBL,double,reg,        C_FLOAT,float,reg)
test_unpack(C_SCPL,SingleComplex,cpl,C_FLOAT,float,reg)
test_unpack(C_DCPL,DoubleComplex,cpl,C_FLOAT,float,reg)
test_unpack(C_INT,int,reg,           C_DBL,double,reg)
test_unpack(C_LONG,long,reg,         C_DBL,double,reg)
test_unpack(C_LONGLONG,long long,reg,C_DBL,double,reg)
test_unpack(C_FLOAT,float,reg,       C_DBL,double,reg)
test_unpack(C_DBL,double,reg,        C_DBL,double,reg)
test_unpack(C_SCPL,SingleComplex,cpl,C_DBL,double,reg)
test_unpack(C_DCPL,DoubleComplex,cpl,C_DBL,double,reg)
test_unpack(C_INT,int,reg,           C_SCPL,SingleComplex,cpl)
test_unpack(C_LONG,long,reg,         C_SCPL,SingleComplex,cpl)
test_unpack(C_LONGLONG,long long,reg,C_SCPL,SingleComplex,cpl)
test_unpack(C_FLOAT,float,reg,       C_SCPL,SingleComplex,cpl)
test_unpack(C_DBL,double,reg,        C_SCPL,SingleComplex,cpl)
test_unpack(C_SCPL,SingleComplex,cpl,C_SCPL,SingleComplex,cpl)
test_unpack(C_DCPL,DoubleComplex,cpl,C_SCPL,SingleComplex,cpl)
test_unpack(C_INT,int,reg,           C_DCPL,DoubleComplex,cpl)
test_unpack(C_LONG,long,reg,         C_DCPL,DoubleComplex,cpl)
test_unpack(C_LONGLONG,long long,reg,C_DCPL,DoubleComplex,cpl)
test_unpack(C_FLOAT,float,reg,       C_DCPL,DoubleComplex,cpl)
test_unpack(C_DBL,double,reg,        C_DCPL,DoubleComplex,cpl)
test_unpack(C_SCPL,SingleComplex,cpl,C_DCPL,DoubleComplex,cpl)
test_unpack(C_DCPL,DoubleComplex,cpl,C_DCPL,DoubleComplex,cpl)



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
#define tests(MT,MT_MSK)                                                 \
            if (0 == me) {                                               \
                printf("testing lo=%d hi=%d q=%d " #MT " " #MT_MSK "\n", \
                        lo, hi, q);                                      \
                fflush(stdout);                                          \
            }                                                            \
            test_unpack_##MT##_##MT_MSK(lo,hi,q)

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
