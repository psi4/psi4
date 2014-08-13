#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* C and/or system headers */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 3rd party headers */
#include <mpi.h>

/* our headers */
#include "comex.h"
#include "comex_impl.h"


/* needed for complex accumulate */
typedef struct {
    double real;
    double imag;
} DoubleComplex;

typedef struct {
    float real;
    float imag;
} SingleComplex;


int comex_nbputs(void *src_ptr, int src_stride_ar[/*stride_levels*/],
                 void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                 int count[/*stride_levels+1*/], int stride_levels,
                 int proc, comex_group_t group,
                 comex_request_t *hdl)
{
    int i, j;
    long src_idx, dst_idx;    /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++) {
        n1dim *= count[i];
    }

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++)
    {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    /* index mangling */
    for(i=0; i<n1dim; i++)
    {
        src_idx = 0;
        dst_idx = 0;
        for(j=1; j<=stride_levels; j++)
        {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0)
                src_bvalue[j]++;
            if(src_bvalue[j] > (count[j]-1))
                src_bvalue[j] = 0;
        }

        for(j=1; j<=stride_levels; j++)
        {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0)
                dst_bvalue[j]++;
            if(dst_bvalue[j] > (count[j]-1))
                dst_bvalue[j] = 0;
        }
        COMEXD_put_nbi((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc);

    }

    comex_wait_proc(proc, group);

    return COMEX_SUCCESS;
}


int comex_nbgets(void *src_ptr, int src_stride_ar[/*stride_levels*/],
                 void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                 int count[/*stride_levels+1*/], int stride_levels,
                 int proc, comex_group_t group,
                 comex_request_t *hdl)
{
    int i, j;
    long src_idx, dst_idx;    /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    
    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;
    
    for(i=2; i<=stride_levels; i++)
    {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    for(i=0; i<n1dim; i++)
    {
        src_idx = 0;
        for(j=1; j<=stride_levels; j++)
        {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0)
                src_bvalue[j]++;
            if(src_bvalue[j] > (count[j]-1))
                src_bvalue[j] = 0;
        }

        dst_idx = 0;

        for(j=1; j<=stride_levels; j++)
        {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0)
                dst_bvalue[j]++;
            if(dst_bvalue[j] > (count[j]-1))
                dst_bvalue[j] = 0;
        }

        COMEXD_get_nbi((char *)src_ptr + src_idx,
                (char *)dst_ptr + dst_idx, count[0], proc);
    }

    COMEXD_waitproc(proc);
    return COMEX_SUCCESS;
}

int comex_nbaccs(int datatype, void *scale,
                 void *src_ptr, int src_stride_ar[/*stride_levels*/],
                 void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                 int count[/*stride_levels+1*/], int stride_levels,
                 int proc, comex_group_t group,
                 comex_request_t *hdl)
{
    int i, j;
    long src_idx, dst_idx;    /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int sizetogetput;
    void *get_buf;

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++)
    {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    sizetogetput = count[0];

    /* TODO: Can we allocate a temporary buffer like we did for the
     * gemini port? */
    if (sizetogetput <= l_state.acc_buf_len) {
        get_buf = l_state.acc_buf;
    }
    else {
        get_buf = comex_malloc_local(sizetogetput);
    }

    assert(get_buf);

    /* grab the atomics lock */
    COMEXD_network_lock(proc);

    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        dst_idx = 0;

        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }

        comex_get((char *)dst_ptr + dst_idx, l_state.acc_buf, 
                count[0], proc, group);

#define EQ_ONE_REG(A) ((A) == 1.0)
#define EQ_ONE_CPL(A) ((A).real == 1.0 && (A).imag == 0.0)
#define IADD_REG(A,B) (A) += (B)
#define IADD_CPL(A,B) (A).real += (B).real; (A).imag += (B).imag
#define IADD_SCALE_REG(A,B,C) (A) += (B) * (C)
#define IADD_SCALE_CPL(A,B,C) (A).real += ((B).real*(C).real) - ((B).imag*(C).imag);\
                              (A).imag += ((B).real*(C).imag) + ((B).imag*(C).real);
#define ACC(WHICH, COMEX_TYPE, C_TYPE)                                      \
        if (datatype == COMEX_TYPE) {                                       \
            int m;                                                          \
            int m_lim = count[0]/sizeof(C_TYPE);                            \
            C_TYPE *iterator = (C_TYPE *)get_buf;                           \
            C_TYPE *value = (C_TYPE *)((char *)src_ptr + src_idx);          \
            C_TYPE calc_scale = *(C_TYPE *)scale;                           \
            if (EQ_ONE_##WHICH(calc_scale)) {                               \
                for (m = 0 ; m < m_lim; ++m) {                              \
                    IADD_##WHICH(iterator[m], value[m]);                    \
                }                                                           \
            }                                                               \
            else {                                                          \
                for (m = 0 ; m < m_lim; ++m) {                              \
                    IADD_SCALE_##WHICH(iterator[m], value[m], calc_scale);  \
                }                                                           \
            }                                                               \
        } else
        ACC(REG, COMEX_ACC_DBL, double)
        ACC(REG, COMEX_ACC_FLT, float)
        ACC(REG, COMEX_ACC_INT, int)
        ACC(REG, COMEX_ACC_LNG, long)
        ACC(CPL, COMEX_ACC_DCP, DoubleComplex)
        ACC(CPL, COMEX_ACC_CPL, SingleComplex)
        {
            assert(0);
        }
#undef ACC
#undef EQ_ONE_REG
#undef EQ_ONE_CPL
#undef IADD_REG
#undef IADD_CPL
#undef IADD_SCALE_REG
#undef IADD_SCALE_CPL

        comex_put((char *)l_state.acc_buf, 
                (char *)dst_ptr + dst_idx, count[0], proc, group);
    }

    COMEXD_waitproc(proc);
    COMEXD_network_unlock(proc);

    if (sizetogetput > l_state.acc_buf_len) {
        comex_free_local(get_buf);
    }

    return COMEX_SUCCESS;
}


int comex_puts(void *src_ptr, int src_stride_ar[/*stride_levels*/],
        void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
        int count[/*stride_levels+1*/], int stride_levels,
        int proc, comex_group_t group)
{
    comex_nbputs(src_ptr, src_stride_ar, dst_ptr,
            dst_stride_ar, count,stride_levels,proc, group, NULL);
    COMEXD_waitproc(proc);

    return COMEX_SUCCESS;
}


int comex_gets(void *src_ptr, int src_stride_ar[/*stride_levels*/],
                   void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                   int count[/*stride_levels+1*/], int stride_levels,
                   int proc, comex_group_t group)
{   
    comex_nbgets(src_ptr, src_stride_ar, dst_ptr,
            dst_stride_ar, count,stride_levels, proc, group, NULL);
    COMEXD_waitproc(proc);

    return COMEX_SUCCESS;
}
    
int comex_accs(int datatype, void *scale,
                   void *src_ptr, int src_stride_ar[/*stride_levels*/],
                   void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                   int count[/*stride_levels+1*/],
                   int stride_levels, int proc, comex_group_t group)
{
    comex_nbaccs(datatype, scale, src_ptr, src_stride_ar, 
            dst_ptr, dst_stride_ar, count, stride_levels, proc, group, NULL);
    COMEXD_waitproc(proc);

    return COMEX_SUCCESS;
}   

