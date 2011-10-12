/*!
  \file
  \brief Matrix multiplication routine (deprecated)
  \ingroup QT
*/

namespace psi {

/*
** NEWMM():  A general matrix multiplication routine with nearly ideal
** memory access patterns and low overhead.  In the near future, loops will
** need to be unrolled in each of the four algorithms below before this
** code can be used for production.
**
** Perform C = alpha A * B + beta C
**
** T. Daniel Crawford
** November, 1996
**
** As fast as I can make it using negligible additional memory.
** For A*B and A*Bt it beats mmult and is very nice.
** For At*B newmm is somewhat slower than mmult.
** For At*Bt the code is very slow.. (See timings.)
**
** Rollin A. King
** July, 1997
**
** Note: Deprecated, C_DGEMM is preferred
**
** \param A         = matrix A
** \param transa    = 1 if A should be transposed
** \param B         = matrix B
** \param transb    = 1 if B should be transposed
** \param C         = matrix for result of A*B
** \param num_rows  = number of rows of A
** \param num_links = columns of A or rows of B
** \param num_cols  = number of columns of A
** \param alpha     = coefficient multiplying A
** \param beta      = coefficient multiplying C
**
** \ingroup QT
*/

void newmm_rking(double **A, int transa, double **B, int transb, double **C,
           int num_rows, int num_links, int num_cols,
           double alpha, double beta)
{
  int i,j,k;
  double sum11,sum12,sum13,sum21,sum22,sum23,sum31,sum32,sum33;
  double tmp11,tmp12,tmp13,tmp21,tmp22,tmp23,tmp31,tmp32,tmp33;
  double *tmpA1,*tmpA2,*tmpA3;
  double *tmpB1,*tmpB2,*tmpB3,*tmpB4,*tmpB5;
  double *tmpC1,*tmpC2,*tmpC3;
  double tmp1,tmp2,tmp3,tmp4;

  if(!num_rows || !num_links || !num_cols) return;

  if(beta == 0.0) {
      for(i=0;i<num_rows;++i)
          for(j=0;j<num_cols;++j)
              C[i][j] = 0.0;
  }
  else if(beta != 1.0)
      for(i=0; i < num_rows; i++)
          for(j=0; j < num_cols; j++)
              C[i][j] *= beta;

  if(!transa) {
      if (alpha != 1.0) {
        for (i=0;i<num_rows;++i)
          for(k=0;k<num_links;++k)
            A[i][k] *= alpha;
      }
      if(!transb) {
          /* C = alpha A * B + beta C */
          for(i=0; i < (num_rows-2); i += 3) {
              tmpC1 = C[i];
              tmpC2 = C[i+1];
              tmpC3 = C[i+2];
              tmpA1 = A[i];
              tmpA2 = A[i+1];
              tmpA3 = A[i+2];
              for(k=0; k < (num_links-2); k += 3) {
                  tmp11 = *(tmpA1+k);
                  tmp12 = *(tmpA1+k+1);
                  tmp13 = *(tmpA1+k+2);
                  tmp21 = *(tmpA2+k);
                  tmp22 = *(tmpA2+k+1);
                  tmp23 = *(tmpA2+k+2);
                  tmp31 = *(tmpA3+k);
                  tmp32 = *(tmpA3+k+1);
                  tmp33 = *(tmpA3+k+2);
                  tmpB1 = B[k];
                  tmpB2 = B[k+1];
                  tmpB3 = B[k+2];
                  for(j=0; j < (num_cols-2); j += 3) {
                    sum11 = tmp11 * *(tmpB1+j);
                    sum12 = tmp11 * *(tmpB1+j+1);
                    sum13 = tmp11 * *(tmpB1+j+2);
                    sum21 = tmp21 * *(tmpB1+j);
                    sum22 = tmp21 * *(tmpB1+j+1);
                    sum23 = tmp21 * *(tmpB1+j+2);
                    sum31 = tmp31 * *(tmpB1+j);
                    sum32 = tmp31 * *(tmpB1+j+1);
                    sum33 = tmp31 * *(tmpB1+j+2);

                    sum11 += tmp12 * *(tmpB2+j);
                    sum12 += tmp12 * *(tmpB2+j+1);
                    sum13 += tmp12 * *(tmpB2+j+2);
                    sum21 += tmp22 * *(tmpB2+j);
                    sum22 += tmp22 * *(tmpB2+j+1);
                    sum23 += tmp22 * *(tmpB2+j+2);
                    sum31 += tmp32 * *(tmpB2+j);
                    sum32 += tmp32 * *(tmpB2+j+1);
                    sum33 += tmp32 * *(tmpB2+j+2);

                    sum11 += tmp13 * *(tmpB3+j);
                    sum12 += tmp13 * *(tmpB3+j+1);
                    sum13 += tmp13 * *(tmpB3+j+2);
                    sum21 += tmp23 * *(tmpB3+j);
                    sum22 += tmp23 * *(tmpB3+j+1);
                    sum23 += tmp23 * *(tmpB3+j+2);
                    sum31 += tmp33 * *(tmpB3+j);
                    sum32 += tmp33 * *(tmpB3+j+1);
                    sum33 += tmp33 * *(tmpB3+j+2);

                    *(tmpC1+j)   += sum11;
                    *(tmpC1+j+1) += sum12;
                    *(tmpC1+j+2) += sum13;
                    *(tmpC2+j)   += sum21;
                    *(tmpC2+j+1) += sum22;
                    *(tmpC2+j+2) += sum23;
                    *(tmpC3+j)   += sum31;
                    *(tmpC3+j+1) += sum32;
                    *(tmpC3+j+2) += sum33;
                  }
                  /* picking up residual j cols */
                  for( ; j < num_cols; ++j) {
                    sum11 = tmp11 * *(tmpB1+j);
                    sum21 = tmp21 * *(tmpB1+j);
                    sum31 = tmp31 * *(tmpB1+j);
                    sum11 += tmp12 * *(tmpB2+j);
                    sum21 += tmp22 * *(tmpB2+j);
                    sum31 += tmp32 * *(tmpB2+j);
                    sum11 += tmp13 * *(tmpB3+j);
                    sum21 += tmp23 * *(tmpB3+j);
                    sum31 += tmp33 * *(tmpB3+j);

                    *(tmpC1+j)   += sum11;
                    *(tmpC2+j)   += sum21;
                    *(tmpC3+j)   += sum31;
                  }
              }
              /* getting rest of k links */
              for ( ; k<num_links;++k) {
                  tmp11 = *(tmpA1+k);
                  tmp21 = *(tmpA2+k);
                  tmp31 = *(tmpA3+k);
                  tmpB1 = B[k];
                  for(j=0; j < (num_cols-2); j += 3) {
                    sum11 = tmp11 * *(tmpB1+j);
                    sum12 = tmp11 * *(tmpB1+j+1);
                    sum13 = tmp11 * *(tmpB1+j+2);
                    sum21 = tmp21 * *(tmpB1+j);
                    sum22 = tmp21 * *(tmpB1+j+1);
                    sum23 = tmp21 * *(tmpB1+j+2);
                    sum31 = tmp31 * *(tmpB1+j);
                    sum32 = tmp31 * *(tmpB1+j+1);
                    sum33 = tmp31 * *(tmpB1+j+2);

                    *(tmpC1+j)   += sum11;
                    *(tmpC1+j+1) += sum12;
                    *(tmpC1+j+2) += sum13;
                    *(tmpC2+j)   += sum21;
                    *(tmpC2+j+1) += sum22;
                    *(tmpC2+j+2) += sum23;
                    *(tmpC3+j)   += sum31;
                    *(tmpC3+j+1) += sum32;
                    *(tmpC3+j+2) += sum33;
                  }
                  /* picking up residual j cols */
                  for( ; j < num_cols; ++j) {
                    sum11 = tmp11 * *(tmpB1+j);
                    sum21 = tmp21 * *(tmpB1+j);
                    sum31 = tmp31 * *(tmpB1+j);

                    *(tmpC1+j)   += sum11;
                    *(tmpC2+j)   += sum21;
                    *(tmpC3+j)   += sum31;
                  }
              }
          }
          /* getting rest of i rows */
          for ( ; i<num_rows;++i) {
              tmpC1 = C[i];
              tmpA1 = A[i];
              for(k=0; k < (num_links-2); k += 3) {
                  tmp11 = *(tmpA1+k);
                  tmp12 = *(tmpA1+k+1);
                  tmp13 = *(tmpA1+k+2);
                  tmpB1 = B[k];
                  tmpB2 = B[k+1];
                  tmpB3 = B[k+2];
                  for(j=0; j < (num_cols-2); j += 3) {
                    sum11 = tmp11 * *(tmpB1+j);
                    sum12 = tmp11 * *(tmpB1+j+1);
                    sum13 = tmp11 * *(tmpB1+j+2);
                    sum11 += tmp12 * *(tmpB2+j);
                    sum12 += tmp12 * *(tmpB2+j+1);
                    sum13 += tmp12 * *(tmpB2+j+2);
                    sum11 += tmp13 * *(tmpB3+j);
                    sum12 += tmp13 * *(tmpB3+j+1);
                    sum13 += tmp13 * *(tmpB3+j+2);
                    *(tmpC1+j)   += sum11;
                    *(tmpC1+j+1) += sum12;
                    *(tmpC1+j+2) += sum13;
                  }
                  /* picking up residual j cols */
                  for( ; j < num_cols; ++j) {
                    sum11 = tmp11 * *(tmpB1+j);
                    sum11 += tmp12 * *(tmpB2+j);
                    sum11 += tmp13 * *(tmpB3+j);
                    *(tmpC1+j)   += sum11;
                  }
              }
              /* getting rest of k links */
              for ( ; k<num_links;++k) {
                  tmp11 = *(tmpA1+k);
                  tmpB1 = B[k];
                  for(j=0; j < (num_cols-2); j += 3) {
                    sum11 = tmp11 * *(tmpB1+j);
                    sum12 = tmp11 * *(tmpB1+j+1);
                    sum13 = tmp11 * *(tmpB1+j+2);
                    *(tmpC1+j)   += sum11;
                    *(tmpC1+j+1) += sum12;
                    *(tmpC1+j+2) += sum13;
                  }
                  /* picking up residual j cols */
                  for( ; j < num_cols; ++j) {
                    sum11 = tmp11 * *(tmpB1+j);
                    *(tmpC1+j)   += sum11;
                  }
              }
          }
      }
      else { /* begin A*Bt */
          for(i=0;i<(num_rows-1);i+=2) {
              tmpA1 = A[i];
              tmpA2 = A[i+1];
              for(j=0;j<(num_cols-1);j+=2) {
                  tmpB1 = B[j];
                  tmpB2 = B[j+1];
                  tmp1  = 0.0;
                  tmp2  = 0.0;
                  tmp3  = 0.0;
                  tmp4  = 0.0;
                  for(k=0;k<(num_links-1);k+=2) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                    tmp1  += *(tmpA1+k+1) * *(tmpB1+k+1);
                    tmp2  += *(tmpA2+k)   * *(tmpB2+k);
                    tmp2  += *(tmpA2+k+1) * *(tmpB2+k+1);
                    tmp3  += *(tmpA1+k)   * *(tmpB2+k);
                    tmp3  += *(tmpA1+k+1) * *(tmpB2+k+1);
                    tmp4  += *(tmpA2+k)   * *(tmpB1+k);
                    tmp4  += *(tmpA2+k+1) * *(tmpB1+k+1);
                  }
                  /* cleaning up rest of k links */
                  for( ; k<num_links;++k) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                    tmp2  += *(tmpA2+k)   * *(tmpB2+k);
                    tmp3  += *(tmpA1+k)   * *(tmpB2+k);
                    tmp4  += *(tmpA2+k)   * *(tmpB1+k);
                  }
                  C[i][j]     += tmp1;
                  C[i][j+1]   += tmp3;
                  C[i+1][j]   += tmp4;
                  C[i+1][j+1] += tmp2;
              }
              /* cleaning up left of j columns */
              for ( ; j<num_cols; ++j) {
                  tmpB1 = B[j];
                  tmp1  = 0.0;
                  tmp4  = 0.0;
                  for(k=0;k<(num_links-1);k+=2) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                    tmp1  += *(tmpA1+k+1) * *(tmpB1+k+1);
                    tmp4  += *(tmpA2+k)   * *(tmpB1+k);
                    tmp4  += *(tmpA2+k+1) * *(tmpB1+k+1);
                  }
                  /* cleaning up rest of k links */
                  for( ; k<num_links;++k) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                    tmp4  += *(tmpA2+k)   * *(tmpB1+k);
                  }
                  C[i][j]     += tmp1;
                  C[i+1][j]   += tmp4;
              }
          }
          /* cleaning up rest of i rows */
          for( ;i<num_rows;i++) {
              tmpA1 = A[i];
              for(j=0;j<(num_cols-1);j+=2) {
                  tmpB1 = B[j];
                  tmpB2 = B[j+1];
                  tmp1  = 0.0;
                  tmp3  = 0.0;
                  for(k=0;k<(num_links-1);k+=2) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                    tmp1  += *(tmpA1+k+1) * *(tmpB1+k+1);
                    tmp3  += *(tmpA1+k)   * *(tmpB2+k);
                    tmp3  += *(tmpA1+k+1) * *(tmpB2+k+1);
                  }
                  /* cleaning up rest of k links */
                  for( ; k<num_links;++k) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                    tmp3  += *(tmpA1+k)   * *(tmpB2+k);
                  }
                  C[i][j]     += tmp1;
                  C[i][j+1]   += tmp3;
              }
              /* cleaning up left of j columns */
              for ( ; j<num_cols; ++j) {
                  tmpB1 = B[j];
                  tmp1  = 0.0;
                  for(k=0;k<(num_links-1);k+=2) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                    tmp1  += *(tmpA1+k+1) * *(tmpB1+k+1);
                  }
                  /* cleaning up rest of k links */
                  for( ; k<num_links;++k) {
                    tmp1  += *(tmpA1+k)   * *(tmpB1+k);
                  }
                  C[i][j] += tmp1;
              }
          }
      }
      if (alpha != 1.0) {
        for (i=0;i<num_rows;++i)
          for(k=0;k<num_links;++k)
            A[i][k] /= alpha;
        }
  }
      else { /* A is transposed */
        if (alpha != 1.0) {
          for (k=0;k<num_links;++k)
            for (i=0;i<num_rows;++i)
              A[k][i] *= alpha;
        }
        if(!transb) { /* At * B begins */
            for(i=0; i<(num_rows-1);i+=2) {
               tmpC1 = C[i];
               tmpC2 = C[i+1];
               for (k=0;k<(num_links-1);k+=2) {
                  tmpB1 = B[k];
                  tmpB2 = B[k+1];
                  tmp11 = A[k][i];
                  tmp12 = A[k+1][i];
                  tmp21 = A[k][i+1];
                  tmp22 = A[k+1][i+1];
                  for (j=0;j<(num_cols-3);j+=4) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                     *(tmpC1+j+1) += tmp11 * *(tmpB1+j+1);
                     *(tmpC1+j+2) += tmp11 * *(tmpB1+j+2);
                     *(tmpC1+j+3) += tmp11 * *(tmpB1+j+3);
                     *(tmpC1+j)   += tmp12 * *(tmpB2+j);
                     *(tmpC1+j+1) += tmp12 * *(tmpB2+j+1);
                     *(tmpC1+j+2) += tmp12 * *(tmpB2+j+2);
                     *(tmpC1+j+3) += tmp12 * *(tmpB2+j+3);
                     *(tmpC2+j)   += tmp21 * *(tmpB1+j);
                     *(tmpC2+j+1) += tmp21 * *(tmpB1+j+1);
                     *(tmpC2+j+2) += tmp21 * *(tmpB1+j+2);
                     *(tmpC2+j+3) += tmp21 * *(tmpB1+j+3);
                     *(tmpC2+j)   += tmp22 * *(tmpB2+j);
                     *(tmpC2+j+1) += tmp22 * *(tmpB2+j+1);
                     *(tmpC2+j+2) += tmp22 * *(tmpB2+j+2);
                     *(tmpC2+j+3) += tmp22 * *(tmpB2+j+3);
                  }
                  /* get rest of j cols */
                  for ( ; j<num_cols;++j) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                     *(tmpC1+j)   += tmp12 * *(tmpB2+j);
                     *(tmpC2+j)   += tmp21 * *(tmpB1+j);
                     *(tmpC2+j)   += tmp22 * *(tmpB2+j);
                  }
               }
               /* get rest of k links */
               for ( ; k<num_links;++k) {
                  tmpB1 = B[k];
                  tmp11 = A[k][i];
                  tmp21 = A[k][i+1];
                  for (j=0;j<(num_cols-3);j+=4) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                     *(tmpC1+j+1) += tmp11 * *(tmpB1+j+1);
                     *(tmpC1+j+2) += tmp11 * *(tmpB1+j+2);
                     *(tmpC1+j+3) += tmp11 * *(tmpB1+j+3);
                     *(tmpC2+j)   += tmp21 * *(tmpB1+j);
                     *(tmpC2+j+1) += tmp21 * *(tmpB1+j+1);
                     *(tmpC2+j+2) += tmp21 * *(tmpB1+j+2);
                     *(tmpC2+j+3) += tmp21 * *(tmpB1+j+3);
                  }
                  /* get rest of j cols */
                  for ( ; j<num_cols;++j) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                     *(tmpC2+j)   += tmp21 * *(tmpB1+j);
                  }
               }
            }
            /* get rest of i rows */
            for ( ; i<num_rows; ++i) {
               tmpC1 = C[i];
               for (k=0;k<(num_links-1);k+=2) {
                  tmpB1 = B[k];
                  tmpB2 = B[k+1];
                  tmp11 = A[k][i];
                  tmp12 = A[k+1][i];
                  for (j=0;j<(num_cols-3);j+=4) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                     *(tmpC1+j+1) += tmp11 * *(tmpB1+j+1);
                     *(tmpC1+j+2) += tmp11 * *(tmpB1+j+2);
                     *(tmpC1+j+3) += tmp11 * *(tmpB1+j+3);
                     *(tmpC1+j)   += tmp12 * *(tmpB2+j);
                     *(tmpC1+j+1) += tmp12 * *(tmpB2+j+1);
                     *(tmpC1+j+2) += tmp12 * *(tmpB2+j+2);
                     *(tmpC1+j+3) += tmp12 * *(tmpB2+j+3);
                  }
                  /* get rest of j cols */
                  for ( ; j<num_cols;++j) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                     *(tmpC1+j)   += tmp12 * *(tmpB2+j);
                  }
               }
               /* get rest of k links */
               for ( ; k<num_links;++k) {
                  tmpB1 = B[k];
                  tmp11 = A[k][i];
                  for (j=0;j<(num_cols-3);j+=4) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                     *(tmpC1+j+1) += tmp11 * *(tmpB1+j+1);
                     *(tmpC1+j+2) += tmp11 * *(tmpB1+j+2);
                     *(tmpC1+j+3) += tmp11 * *(tmpB1+j+3);
                  }
                  /* get rest of j cols */
                  for ( ; j<num_cols;++j) {
                     *(tmpC1+j)   += tmp11 * *(tmpB1+j);
                  }
               }
            } /* end rest of i rows */
        }

        else {   /* if At and Bt */
            for(i=0;i<(num_rows-2);i+=3) {
               tmpC1 = C[i];
               tmpC2 = C[i+1];
               tmpC3 = C[i+2];
               for(k=0;k<(num_links-2);k+=3) {
                 tmp11 = A[k][i];
                 tmp12 = A[k][i+1];
                 tmp13 = A[k][i+2];
                 tmp21 = A[k+1][i];
                 tmp22 = A[k+1][i+1];
                 tmp23 = A[k+1][i+2];
                 tmp31 = A[k+2][i];
                 tmp32 = A[k+2][i+1];
                 tmp33 = A[k+2][i+2];
                 for(j=0;j<(num_cols-4);j+=5) {
                    tmpB1 = B[j];
                    tmpB2 = B[j+1];
                    tmpB3 = B[j+2];
                    tmpB4 = B[j+3];
                    tmpB5 = B[j+4];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                    *(tmpC1+j+1) += tmp11 * *(tmpB2+k);
                    *(tmpC1+j+2) += tmp11 * *(tmpB3+k);
                    *(tmpC1+j+3) += tmp11 * *(tmpB4+k);
                    *(tmpC1+j+4) += tmp11 * *(tmpB5+k);
                    *(tmpC1+j)   += tmp21 * *(tmpB1+k+1);
                    *(tmpC1+j+1) += tmp21 * *(tmpB2+k+1);
                    *(tmpC1+j+2) += tmp21 * *(tmpB3+k+1);
                    *(tmpC1+j+3) += tmp21 * *(tmpB4+k+1);
                    *(tmpC1+j+4) += tmp21 * *(tmpB5+k+1);
                    *(tmpC1+j)   += tmp31 * *(tmpB1+k+2);
                    *(tmpC1+j+1) += tmp31 * *(tmpB2+k+2);
                    *(tmpC1+j+2) += tmp31 * *(tmpB3+k+2);
                    *(tmpC1+j+3) += tmp31 * *(tmpB4+k+2);
                    *(tmpC1+j+4) += tmp31 * *(tmpB5+k+2);
                    *(tmpC2+j)   += tmp12 * *(tmpB1+k);
                    *(tmpC2+j+1) += tmp12 * *(tmpB2+k);
                    *(tmpC2+j+2) += tmp12 * *(tmpB3+k);
                    *(tmpC2+j+3) += tmp12 * *(tmpB4+k);
                    *(tmpC2+j+4) += tmp12 * *(tmpB5+k);
                    *(tmpC2+j)   += tmp22 * *(tmpB1+k+1);
                    *(tmpC2+j+1) += tmp22 * *(tmpB2+k+1);
                    *(tmpC2+j+2) += tmp22 * *(tmpB3+k+1);
                    *(tmpC2+j+3) += tmp22 * *(tmpB4+k+1);
                    *(tmpC2+j+4) += tmp22 * *(tmpB5+k+1);
                    *(tmpC2+j)   += tmp32 * *(tmpB1+k+2);
                    *(tmpC2+j+1) += tmp32 * *(tmpB2+k+2);
                    *(tmpC2+j+2) += tmp32 * *(tmpB3+k+2);
                    *(tmpC2+j+3) += tmp32 * *(tmpB4+k+2);
                    *(tmpC2+j+4) += tmp32 * *(tmpB5+k+2);
                    *(tmpC3+j)   += tmp13 * *(tmpB1+k);
                    *(tmpC3+j+1) += tmp13 * *(tmpB2+k);
                    *(tmpC3+j+2) += tmp13 * *(tmpB3+k);
                    *(tmpC3+j+3) += tmp13 * *(tmpB4+k);
                    *(tmpC3+j+4) += tmp13 * *(tmpB5+k);
                    *(tmpC3+j)   += tmp23 * *(tmpB1+k+1);
                    *(tmpC3+j+1) += tmp23 * *(tmpB2+k+1);
                    *(tmpC3+j+2) += tmp23 * *(tmpB3+k+1);
                    *(tmpC3+j+3) += tmp23 * *(tmpB4+k+1);
                    *(tmpC3+j+4) += tmp23 * *(tmpB5+k+1);
                    *(tmpC3+j)   += tmp33 * *(tmpB1+k+2);
                    *(tmpC3+j+1) += tmp33 * *(tmpB2+k+2);
                    *(tmpC3+j+2) += tmp33 * *(tmpB3+k+2);
                    *(tmpC3+j+3) += tmp33 * *(tmpB4+k+2);
                    *(tmpC3+j+4) += tmp33 * *(tmpB5+k+2);
                 }
                 /* get rest of j cols */
                 for ( ; j<num_cols;++j) {
                    tmpB1 = B[j];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                    *(tmpC1+j)   += tmp21 * *(tmpB1+k+1);
                    *(tmpC1+j)   += tmp31 * *(tmpB1+k+2);
                    *(tmpC2+j)   += tmp12 * *(tmpB1+k);
                    *(tmpC2+j)   += tmp22 * *(tmpB1+k+1);
                    *(tmpC2+j)   += tmp32 * *(tmpB1+k+2);
                    *(tmpC3+j)   += tmp13 * *(tmpB1+k);
                    *(tmpC3+j)   += tmp23 * *(tmpB1+k+1);
                    *(tmpC3+j)   += tmp33 * *(tmpB1+k+2);
                 }
               }
               /* get rest of k links */
               for ( ; k<num_links; ++k) {
                 tmp11 = A[k][i];
                 tmp12 = A[k][i+1];
                 tmp13 = A[k][i+2];
                 for(j=0;j<(num_cols-4);j+=5) {
                    tmpB1 = B[j];
                    tmpB2 = B[j+1];
                    tmpB3 = B[j+2];
                    tmpB4 = B[j+3];
                    tmpB5 = B[j+4];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                    *(tmpC1+j+1) += tmp11 * *(tmpB2+k);
                    *(tmpC1+j+2) += tmp11 * *(tmpB3+k);
                    *(tmpC1+j+3) += tmp11 * *(tmpB4+k);
                    *(tmpC1+j+4) += tmp11 * *(tmpB5+k);
                    *(tmpC2+j)   += tmp12 * *(tmpB1+k);
                    *(tmpC2+j+1) += tmp12 * *(tmpB2+k);
                    *(tmpC2+j+2) += tmp12 * *(tmpB3+k);
                    *(tmpC2+j+3) += tmp12 * *(tmpB4+k);
                    *(tmpC2+j+4) += tmp12 * *(tmpB5+k);
                    *(tmpC3+j)   += tmp13 * *(tmpB1+k);
                    *(tmpC3+j+1) += tmp13 * *(tmpB2+k);
                    *(tmpC3+j+2) += tmp13 * *(tmpB3+k);
                    *(tmpC3+j+3) += tmp13 * *(tmpB4+k);
                    *(tmpC3+j+4) += tmp13 * *(tmpB5+k);
                 }
                 /* get rest of j cols */
                 for ( ; j<num_cols;++j) {
                    tmpB1 = B[j];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                    *(tmpC2+j)   += tmp12 * *(tmpB1+k);
                    *(tmpC3+j)   += tmp13 * *(tmpB1+k);
                 }
               }
            }
            /* get rest of i rows */
            for ( ; i<num_rows;++i) {
               tmpC1 = C[i];
               for(k=0;k<(num_links-2);k+=3) {
                 tmp11 = A[k][i];
                 tmp21 = A[k+1][i];
                 tmp31 = A[k+2][i];
                 for(j=0;j<(num_cols-4);j+=5) {
                    tmpB1 = B[j];
                    tmpB2 = B[j+1];
                    tmpB3 = B[j+2];
                    tmpB4 = B[j+3];
                    tmpB5 = B[j+4];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                    *(tmpC1+j+1) += tmp11 * *(tmpB2+k);
                    *(tmpC1+j+2) += tmp11 * *(tmpB3+k);
                    *(tmpC1+j+3) += tmp11 * *(tmpB4+k);
                    *(tmpC1+j+4) += tmp11 * *(tmpB5+k);
                    *(tmpC1+j)   += tmp21 * *(tmpB1+k+1);
                    *(tmpC1+j+1) += tmp21 * *(tmpB2+k+1);
                    *(tmpC1+j+2) += tmp21 * *(tmpB3+k+1);
                    *(tmpC1+j+3) += tmp21 * *(tmpB4+k+1);
                    *(tmpC1+j+4) += tmp21 * *(tmpB5+k+1);
                    *(tmpC1+j)   += tmp31 * *(tmpB1+k+2);
                    *(tmpC1+j+1) += tmp31 * *(tmpB2+k+2);
                    *(tmpC1+j+2) += tmp31 * *(tmpB3+k+2);
                    *(tmpC1+j+3) += tmp31 * *(tmpB4+k+2);
                    *(tmpC1+j+4) += tmp31 * *(tmpB5+k+2);
                 }
                 /* get rest of j cols */
                 for ( ; j<num_cols;++j) {
                    tmpB1 = B[j];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                    *(tmpC1+j)   += tmp21 * *(tmpB1+k+1);
                    *(tmpC1+j)   += tmp31 * *(tmpB1+k+2);
                 }
               }
               /* get rest of k links */
               for ( ; k<num_links; ++k) {
                 tmp11 = A[k][i];
                 for(j=0;j<(num_cols-4);j+=5) {
                    tmpB1 = B[j];
                    tmpB2 = B[j+1];
                    tmpB3 = B[j+2];
                    tmpB4 = B[j+3];
                    tmpB5 = B[j+4];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                    *(tmpC1+j+1) += tmp11 * *(tmpB2+k);
                    *(tmpC1+j+2) += tmp11 * *(tmpB3+k);
                    *(tmpC1+j+3) += tmp11 * *(tmpB4+k);
                    *(tmpC1+j+4) += tmp11 * *(tmpB5+k);
                 }
                 /* get rest of j cols */
                 for ( ; j<num_cols;++j) {
                    tmpB1 = B[j];
                    *(tmpC1+j)   += tmp11 * *(tmpB1+k);
                 }
               }
            }
        }

        if (alpha != 1.0) {
          for (k=0;k<num_links;++k)
            for (i=0;i<num_rows;++i)
              A[k][i] /= alpha;
        }
      }
}

}

