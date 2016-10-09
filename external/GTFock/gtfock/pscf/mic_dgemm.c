#include "mic_dgemm.h"


#ifdef __INTEL_OFFLOAD
static int max_pcl_matrix_size = -1;
static int usemic = 0;
#pragma offload_attribute(push, target(mic))
double *pcl_a_mic;
double *pcl_b_mic;
double *pcl_c_mic;
void cblas_dgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                  CBLAS_TRANSPOSE TransB, MKL_INT M,
                  MKL_INT N, MKL_INT K,
                  double alpha, double *A, MKL_INT lda,
                  double *B, MKL_INT ldb,
                  double beta, double *C, MKL_INT ldc);
#pragma offload_attribute(pop)
#endif /* __INTEL_OFFLOAD */


void init_pcl_dgemm (int matrix_size, int usemic_)
{
#ifdef __INTEL_OFFLOAD
    usemic = usemic_;
    if (!usemic)
        return;
    max_pcl_matrix_size = matrix_size * 2;
    pcl_a_mic =
        (double *) _mm_malloc (max_pcl_matrix_size * max_pcl_matrix_size *
                               sizeof (double), 64);
    pcl_b_mic =
        (double *) _mm_malloc (max_pcl_matrix_size * max_pcl_matrix_size *
                               sizeof (double), 64);
    pcl_c_mic =
        (double *) _mm_malloc (max_pcl_matrix_size * max_pcl_matrix_size *
                               sizeof (double), 64);
    assert (pcl_a_mic != NULL &&
            pcl_b_mic != NULL &&
            pcl_c_mic != NULL);
    #pragma offload target(mic : 0) \
        in(pcl_a_mic: length(max_pcl_matrix_size*max_pcl_matrix_size) ALLOC align(64)) \
        in(pcl_b_mic: length(max_pcl_matrix_size*max_pcl_matrix_size) ALLOC align(64)) \
        in(pcl_c_mic: length(max_pcl_matrix_size*max_pcl_matrix_size) ALLOC align(64))
    {
    }
#else
    matrix_size = usemic_ = 0;
#endif
}


void deinit_pcl_dgemm (void)
{
#ifdef __INTEL_OFFLOAD
    if (!usemic)
        return;
    
    #pragma offload target(mic : 0) \
        in(pcl_a_mic: length(max_pcl_matrix_size*max_pcl_matrix_size) FREE align(64)) \
        in(pcl_b_mic: length(max_pcl_matrix_size*max_pcl_matrix_size) FREE align(64)) \
        in(pcl_c_mic: length(max_pcl_matrix_size*max_pcl_matrix_size) FREE align(64))

    _mm_free (pcl_a_mic);
    _mm_free (pcl_b_mic);
    _mm_free (pcl_c_mic);
#endif
}


__attribute__ ((noinline))
void pcl_dgemm (CBLAS_ORDER order, CBLAS_TRANSPOSE transA,
                CBLAS_TRANSPOSE transB, int M, int N, int K,
                double alpha, double *A, int lda, double *B, int ldb,
                double beta, double *C, int ldc)
{
    assert (order == CblasRowMajor);

#ifdef __INTEL_OFFLOAD
    if (usemic && M >= 800)
    {
        assert ((M * lda) < max_pcl_matrix_size * max_pcl_matrix_size);
        assert ((K * ldb) < max_pcl_matrix_size * max_pcl_matrix_size);
        assert ((M * ldc) < max_pcl_matrix_size * max_pcl_matrix_size);
        assert (beta == 0.0);
        #pragma offload target(mic : 0) \
            in(order, M,N,K,lda,ldb,ldc,alpha,beta) \
            in(A[0:(M*lda)]: into(pcl_a_mic[0:(M*lda)]) REUSE align(64))\
            in(B[0:(K*ldb)]: into(pcl_b_mic[0:(K*ldb)]) REUSE align(64))\
            out(pcl_c_mic[0:(M*ldc)]: into(C[0:(M*ldc)]) REUSE align(64))
        {
            cblas_dgemm (order, transA, transB, M, N, K, alpha, pcl_a_mic,
                         lda, pcl_b_mic, ldb, beta, pcl_c_mic, ldc);
        }        
    }
    else
#endif        
    {
        cblas_dgemm (order, transA, transB, M, N, K, alpha, A, lda,
                     B, ldb, beta, C, ldc);
    }    
}
