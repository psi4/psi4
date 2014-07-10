#ifndef XGEMM_H_
#define XGEMM_H_

void xb_sgemm(char *transa, char *transb, int *M, int *N, int *K,
              float *alpha, const float *a, int *p_lda, const float *b,
              int *p_ldb, float *beta, float *c, int *p_ldc);
void xb_dgemm(char *transa, char *transb, int *M, int *N, int *K,
              double *alpha, const double *a, int *p_lda, const double *b,
              int *p_ldb, double *beta, double *c, int *p_ldc);
void xb_zgemm(char * transa, char *transb, int *M, int *N, int *K,
              const void *alpha, const void *a, int *p_lda, const void *b,
              int *p_ldb, const void *beta, void *c, int *p_ldc);
void xb_cgemm(char * transa, char *transb, int *M, int *N, int *K,
              const void *alpha, const void *a, int *p_lda, const void *b,
              int *p_ldb, const void *beta, void *c, int *p_ldc);

#endif /* XGEMM_H_ */
