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

#if ENABLE_F77
#   ifdef WIN32
#       define ATR __stdcall
#   else
#       define ATR
#   endif
#   define i_accumulate_1d_    F77_FUNC_(i_accumulate_1d,I_ACCUMULATE_2D)
#   define l_accumulate_1d_              c_l_accumulate_1d_
#   define ll_accumulate_1d_             c_ll_accumulate_1d_
#   define f_accumulate_1d_    F77_FUNC_(f_accumulate_1d,F_ACCUMULATE_2D)
#   define d_accumulate_1d_    F77_FUNC_(d_accumulate_1d,D_ACCUMULATE_2D)
#   define c_accumulate_1d_    F77_FUNC_(c_accumulate_1d,C_ACCUMULATE_2D)
#   define z_accumulate_1d_    F77_FUNC_(z_accumulate_1d,Z_ACCUMULATE_2D)
#   define i_accumulate_2d_    F77_FUNC_(i_accumulate_2d,I_ACCUMULATE_2D)
#   define l_accumulate_2d_              c_l_accumulate_2d_
#   define ll_accumulate_2d_             c_ll_accumulate_2d_
#   define f_accumulate_2d_    F77_FUNC_(f_accumulate_2d,F_ACCUMULATE_2D)
#   define d_accumulate_2d_    F77_FUNC_(d_accumulate_2d,D_ACCUMULATE_2D)
#   define c_accumulate_2d_    F77_FUNC_(c_accumulate_2d,C_ACCUMULATE_2D)
#   define z_accumulate_2d_    F77_FUNC_(z_accumulate_2d,Z_ACCUMULATE_2D)
#   define i_accumulate_2d_u_  F77_FUNC_(i_accumulate_2d_u,I_ACCUMULATE_2D_U)
#   define l_accumulate_2d_u_            c_l_accumulate_2d_u_
#   define ll_accumulate_2d_u_           c_ll_accumulate_2d_u_
#   define f_accumulate_2d_u_  F77_FUNC_(f_accumulate_2d_u,F_ACCUMULATE_2D_U)
#   define d_accumulate_2d_u_  F77_FUNC_(d_accumulate_2d_u,D_ACCUMULATE_2D_U)
#   define c_accumulate_2d_u_  F77_FUNC_(c_accumulate_2d_u,C_ACCUMULATE_2D_U)
#   define z_accumulate_2d_u_  F77_FUNC_(z_accumulate_2d_u,Z_ACCUMULATE_2D_U)
#   define fort_dadd_          F77_FUNC_(fort_dadd,FORT_DADD)
#   define fort_dadd2_         F77_FUNC_(fort_dadd2,FORT_DADD2)
#   define fort_dmult_         F77_FUNC_(fort_dmult,FORT_DMULT)
#   define fort_dmult2_        F77_FUNC_(fort_dmult2,FORT_DMULT2)
void ATR d_accumulate_1d_(const double* const restrict alpha,
                                double* const restrict A,
                          const double* const restrict B,
                          const int*    const restrict rows);
void ATR f_accumulate_1d_(const float* const restrict alpha,
                                float* const restrict A,
                          const float* const restrict B,
                          const int*   const restrict rows);
void ATR c_accumulate_1d_(const complex_t* const restrict alpha,
                                complex_t* const restrict A,
                          const complex_t* const restrict B,
                          const int*       const restrict rows);
void ATR z_accumulate_1d_(const dcomplex_t* const restrict alpha,
                                dcomplex_t* const restrict A,
                          const dcomplex_t* const restrict B,
                          const int*        const restrict rows);
void ATR i_accumulate_1d_(const int* const restrict alpha,
                                int* const restrict A,
                          const int* const restrict B,
                          const int* const restrict rows);
void ATR l_accumulate_1d_(const long* const restrict alpha,
                                long* const restrict A,
                          const long* const restrict B,
                          const int*  const restrict rows);
void ATR ll_accumulate_1d_(const long long* const restrict alpha,
                                 long long* const restrict A,
                           const long long* const restrict B,
                           const int*       const restrict rows);

void ATR d_accumulate_2d_(const double* const restrict alpha,
                          const int*    const restrict rows,
                          const int*    const restrict cols,
                                double* const restrict A,
                          const int*    const restrict ald,
                          const double* const restrict B,
                          const int*    const restrict bld);
void ATR f_accumulate_2d_(const float*  const restrict alpha,
                          const int*    const restrict rows,
                          const int*    const restrict cols,
                                float*  const restrict A,
                          const int*    const restrict ald,
                          const float*  const restrict B,
                          const int*    const restrict bld);
void ATR c_accumulate_2d_(const complex_t* const restrict alpha,
                          const int*       const restrict rows,
                          const int*       const restrict cols,
                                complex_t* const restrict A,
                          const int*       const restrict ald,
                          const complex_t* const restrict B,
                          const int*       const restrict bld);
void ATR z_accumulate_2d_(const dcomplex_t* const restrict alpha,
                          const int*        const restrict rows,
                          const int*        const restrict cols,
                                dcomplex_t* const restrict A,
                          const int*        const restrict ald,
                          const dcomplex_t* const restrict B,
                          const int*        const restrict bld);
void ATR i_accumulate_2d_(const int* const restrict alpha,
                          const int* const restrict rows,
                          const int* const restrict cols,
                                int* const restrict A,
                          const int* const restrict ald,
                          const int* const restrict B,
                          const int* const restrict bld);
void ATR l_accumulate_2d_(const long* const restrict alpha,
                          const int*  const restrict rows,
                          const int*  const restrict cols,
                                long* const restrict A,
                          const int*  const restrict ald,
                          const long* const restrict B,
                          const int*  const restrict bld);
void ATR ll_accumulate_2d_(const long long* const restrict alpha,
                           const int*       const restrict rows,
                           const int*       const restrict cols,
                                 long long* const restrict A,
                           const int*       const restrict ald,
                           const long long* const restrict B,
                           const int*       const restrict bld);

void ATR d_accumulate_2d_u_(const double* const restrict alpha,
                            const int*    const restrict rows,
                            const int*    const restrict cols,
                                  double* const restrict A,
                            const int*    const restrict ald,
                            const double* const restrict B,
                            const int*    const restrict bld);
void ATR f_accumulate_2d_u_(const float*  const restrict alpha,
                            const int*    const restrict rows,
                            const int*    const restrict cols,
                                  float*  const restrict A,
                            const int*    const restrict ald,
                            const float*  const restrict B,
                            const int*    const restrict bld);
void ATR c_accumulate_2d_u_(const complex_t* const restrict alpha,
                            const int*       const restrict rows,
                            const int*       const restrict cols,
                                  complex_t* const restrict A,
                            const int*       const restrict ald,
                            const complex_t* const restrict B,
                            const int*       const restrict bld);
void ATR z_accumulate_2d_u_(const dcomplex_t* const restrict alpha,
                            const int*        const restrict rows,
                            const int*        const restrict cols,
                                  dcomplex_t* const restrict A,
                            const int*        const restrict ald,
                            const dcomplex_t* const restrict B,
                            const int*        const restrict bld);
void ATR i_accumulate_2d_u_(const int* const restrict alpha,
                            const int* const restrict rows,
                            const int* const restrict cols,
                                  int* const restrict A,
                            const int* const restrict ald,
                            const int* const restrict B,
                            const int* const restrict bld);
void ATR l_accumulate_2d_u_(const long* const restrict alpha,
                            const int*  const restrict rows,
                            const int*  const restrict cols,
                                  long* const restrict A,
                            const int*  const restrict ald,
                            const long* const restrict B,
                            const int*  const restrict bld);
void ATR ll_accumulate_2d_u_(const long long* const restrict alpha,
                             const int*       const restrict rows,
                             const int*       const restrict cols,
                                   long long* const restrict A,
                             const int*       const restrict ald,
                             const long long* const restrict B,
                             const int*       const restrict bld);

void ATR fort_dadd_(const int*    const restrict n,
                          double* const restrict x,
                    const double* const restrict work);
void ATR fort_dadd2_(const int*    const restrict n,
                           double* const restrict x,
                     const double* const restrict work,
                     const double* const restrict work2);
void ATR fort_dmult_(const int*    const restrict n,
                           double* const restrict x,
                     const double* const restrict work);
void ATR fort_dmult2_(const int*    const restrict n,
                            double* const restrict x,
                      const double* const restrict work,
                      const double* const restrict work2);
#endif

#if NOFORT
#   define I_ACCUMULATE_1D   c_i_accumulate_1d_
#   define L_ACCUMULATE_1D   c_l_accumulate_1d_
#   define LL_ACCUMULATE_1D c_ll_accumulate_1d_
#   define D_ACCUMULATE_1D   c_d_accumulate_1d_
#   define C_ACCUMULATE_1D   c_c_accumulate_1d_
#   define Z_ACCUMULATE_1D   c_z_accumulate_1d_
#   define F_ACCUMULATE_1D   c_f_accumulate_1d_
#   define I_ACCUMULATE_2D   c_i_accumulate_2d_
#   define L_ACCUMULATE_2D   c_l_accumulate_2d_
#   define LL_ACCUMULATE_2D c_ll_accumulate_2d_
#   define D_ACCUMULATE_2D   c_d_accumulate_2d_
#   define C_ACCUMULATE_2D   c_c_accumulate_2d_
#   define Z_ACCUMULATE_2D   c_z_accumulate_2d_
#   define F_ACCUMULATE_2D   c_f_accumulate_2d_
#   define FORT_DADD   c_dadd_
#   define FORT_DADD2  c_dadd2_
#   define FORT_DMULT  c_dmult_
#   define FORT_DMULT2 c_dmult2_
#else
#   if defined(AIX) || defined(BGML) || defined(SGI_)
#       define I_ACCUMULATE_2D     i_accumulate_2d_u_
#       define L_ACCUMULATE_2D   c_l_accumulate_2d_u_
#       define LL_ACCUMULATE_2D c_ll_accumulate_2d_u_
#       define D_ACCUMULATE_2D     d_accumulate_2d_u_
#       define C_ACCUMULATE_2D     c_accumulate_2d_u_
#       define Z_ACCUMULATE_2D     z_accumulate_2d_u_
#       define F_ACCUMULATE_2D     f_accumulate_2d_u_
#   else
#       define I_ACCUMULATE_2D     i_accumulate_2d_
#       define L_ACCUMULATE_2D   c_l_accumulate_2d_
#       define LL_ACCUMULATE_2D c_ll_accumulate_2d_
#       define D_ACCUMULATE_2D     d_accumulate_2d_
#       define C_ACCUMULATE_2D     c_accumulate_2d_
#       define Z_ACCUMULATE_2D     z_accumulate_2d_
#       define F_ACCUMULATE_2D     f_accumulate_2d_
#   endif
#   if defined(CRAY) && !defined(__crayx1)
#       undef  D_ACCUMULATE_2D 
#       define D_ACCUMULATE_2D F77_FUNC_(daxpy_2d,DAXPY_2D)
#   endif
#   define I_ACCUMULATE_1D     i_accumulate_1d_
#   define L_ACCUMULATE_1D   c_l_accumulate_1d_
#   define LL_ACCUMULATE_1D c_ll_accumulate_1d_
#   define D_ACCUMULATE_1D     d_accumulate_1d_
#   define C_ACCUMULATE_1D     c_accumulate_1d_
#   define Z_ACCUMULATE_1D     z_accumulate_1d_
#   define F_ACCUMULATE_1D     f_accumulate_1d_
#   define FORT_DADD   fort_dadd_
#   define FORT_DADD2  fort_dadd2_
#   define FORT_DMULT  fort_dmult_
#   define FORT_DMULT2 fort_dmult2_
#endif /* !NOFORT */

// specific to src-gemini
#if defined(AIX) || defined(NOUNDERSCORE)
#   define RA_ACCUMULATE_2D ra_accumulate_2d_u
#elif defined(BGML)
#   define RA_ACCUMULATE_2D ra_accumulate_2d_u__
#elif defined(SGI_)
#   define RA_ACCUMULATE_2D RA_ACCUMULATE_2D_
#elif !defined(CRAY) && !defined(WIN32) && !defined(HITACHI) ||defined(__crayx1)
#   define RA_ACCUMULATE_2D RA_ACCUMULATE_2D_
#endif

#ifndef CRAY_T3E
void ATR RA_ACCUMULATE_2D(long*, int*, int*, long*, int*, long*, int*);
#else
#define RA_ACCUMULATE_2D RA_ACCUMULATE_2D_
void RA_ACCUMULATE_2D_(long*, int*, int*, long*, int*, long*, int*);
#endif

#endif /* _ACC_H_ */
