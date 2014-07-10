#ifndef _ACC_H_
#define _ACC_H_

typedef struct {
    float real;
    float imag;
} complex_t;

typedef struct {
    double real;
    double imag;
} dcomplex_t;

void c_d_accumulate_1d_(const double* const restrict alpha,
                              double* const restrict A,
                        const double* const restrict B,
                        const int*    const restrict rows);
void c_f_accumulate_1d_(const float* const restrict alpha,
                              float* const restrict A,
                        const float* const restrict B,
                        const int*   const restrict rows);
void c_c_accumulate_1d_(const complex_t* const restrict alpha,
                              complex_t* const restrict A,
                        const complex_t* const restrict B,
                        const int*       const restrict rows);
void c_z_accumulate_1d_(const dcomplex_t* const restrict alpha,
                              dcomplex_t* const restrict A,
                        const dcomplex_t* const restrict B,
                        const int*        const restrict rows);
void c_i_accumulate_1d_(const int* const restrict alpha,
                              int* const restrict A,
                        const int* const restrict B,
                        const int* const restrict rows);
void c_l_accumulate_1d_(const long* const restrict alpha,
                              long* const restrict A,
                        const long* const restrict B,
                        const int*  const restrict rows);
void c_ll_accumulate_1d_(const long long* const restrict alpha,
                               long long* const restrict A,
                         const long long* const restrict B,
                         const int*       const restrict rows);

void c_d_accumulate_2d_(const double* const restrict alpha,
                        const int*    const restrict rows,
                        const int*    const restrict cols,
                              double* const restrict A,
                        const int*    const restrict ald,
                        const double* const restrict B,
                        const int*    const restrict bld);
void c_f_accumulate_2d_(const float*  const restrict alpha,
                        const int*    const restrict rows,
                        const int*    const restrict cols,
                              float*  const restrict A,
                        const int*    const restrict ald,
                        const float*  const restrict B,
                        const int*    const restrict bld);
void c_c_accumulate_2d_(const complex_t* const restrict alpha,
                        const int*       const restrict rows,
                        const int*       const restrict cols,
                              complex_t* const restrict A,
                        const int*       const restrict ald,
                        const complex_t* const restrict B,
                        const int*       const restrict bld);
void c_z_accumulate_2d_(const dcomplex_t* const restrict alpha,
                        const int*        const restrict rows,
                        const int*        const restrict cols,
                              dcomplex_t* const restrict A,
                        const int*        const restrict ald,
                        const dcomplex_t* const restrict B,
                        const int*        const restrict bld);
void c_i_accumulate_2d_(const int* const restrict alpha,
                        const int* const restrict rows,
                        const int* const restrict cols,
                              int* const restrict A,
                        const int* const restrict ald,
                        const int* const restrict B,
                        const int* const restrict bld);
void c_l_accumulate_2d_(const long* const restrict alpha,
                        const int*  const restrict rows,
                        const int*  const restrict cols,
                              long* const restrict A,
                        const int*  const restrict ald,
                        const long* const restrict B,
                        const int*  const restrict bld);
void c_ll_accumulate_2d_(const long long* const restrict alpha,
                         const int*       const restrict rows,
                         const int*       const restrict cols,
                               long long* const restrict A,
                         const int*       const restrict ald,
                         const long long* const restrict B,
                         const int*       const restrict bld);

void c_d_accumulate_2d_u_(const double* const restrict alpha,
                          const int*    const restrict rows,
                          const int*    const restrict cols,
                                double* const restrict A,
                          const int*    const restrict ald,
                          const double* const restrict B,
                          const int*    const restrict bld);
void c_f_accumulate_2d_u_(const float*  const restrict alpha,
                          const int*    const restrict rows,
                          const int*    const restrict cols,
                                float*  const restrict A,
                          const int*    const restrict ald,
                          const float*  const restrict B,
                          const int*    const restrict bld);
void c_c_accumulate_2d_u_(const complex_t* const restrict alpha,
                          const int*       const restrict rows,
                          const int*       const restrict cols,
                                complex_t* const restrict A,
                          const int*       const restrict ald,
                          const complex_t* const restrict B,
                          const int*       const restrict bld);
void c_z_accumulate_2d_u_(const dcomplex_t* const restrict alpha,
                          const int*        const restrict rows,
                          const int*        const restrict cols,
                                dcomplex_t* const restrict A,
                          const int*        const restrict ald,
                          const dcomplex_t* const restrict B,
                          const int*        const restrict bld);
void c_i_accumulate_2d_u_(const int* const restrict alpha,
                          const int* const restrict rows,
                          const int* const restrict cols,
                                int* const restrict A,
                          const int* const restrict ald,
                          const int* const restrict B,
                          const int* const restrict bld);
void c_l_accumulate_2d_u_(const long* const restrict alpha,
                          const int*  const restrict rows,
                          const int*  const restrict cols,
                                long* const restrict A,
                          const int*  const restrict ald,
                          const long* const restrict B,
                          const int*  const restrict bld);
void c_ll_accumulate_2d_u_(const long long* const restrict alpha,
                           const int*       const restrict rows,
                           const int*       const restrict cols,
                                 long long* const restrict A,
                           const int*       const restrict ald,
                           const long long* const restrict B,
                           const int*       const restrict bld);

void RA_ACCUMULATE_2D_(long*, int*, int*, long*, int*, long*, int*);

void c_dadd_(const int*    const restrict n,
                   double* const restrict x,
             const double* const restrict work);
void c_dadd2_(const int*    const restrict n,
                    double* const restrict x,
              const double* const restrict work,
              const double* const restrict work2);
void c_dmult_(const int*    const restrict n,
                    double* const restrict x,
              const double* const restrict work);
void c_dmult2_(const int*    const restrict n,
                     double* const restrict x,
               const double* const restrict work,
               const double* const restrict work2);

#define I_ACCUMULATE_1D   c_i_accumulate_1d_
#define L_ACCUMULATE_1D   c_l_accumulate_1d_
#define LL_ACCUMULATE_1D c_ll_accumulate_1d_
#define D_ACCUMULATE_1D   c_d_accumulate_1d_
#define C_ACCUMULATE_1D   c_c_accumulate_1d_
#define Z_ACCUMULATE_1D   c_z_accumulate_1d_
#define F_ACCUMULATE_1D   c_f_accumulate_1d_
#define I_ACCUMULATE_2D   c_i_accumulate_2d_
#define L_ACCUMULATE_2D   c_l_accumulate_2d_
#define LL_ACCUMULATE_2D c_ll_accumulate_2d_
#define D_ACCUMULATE_2D   c_d_accumulate_2d_
#define C_ACCUMULATE_2D   c_c_accumulate_2d_
#define Z_ACCUMULATE_2D   c_z_accumulate_2d_
#define F_ACCUMULATE_2D   c_f_accumulate_2d_
#define FORT_DADD   c_dadd_
#define FORT_DADD2  c_dadd2_
#define FORT_DMULT  c_dmult_
#define FORT_DMULT2 c_dmult2_

#endif /* _ACC_H_ */
