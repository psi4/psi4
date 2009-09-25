#include "algebra_interface.h"
#include <libmoinfo/libmoinfo.h>

extern FILE *infile, *outfile;

namespace psi{ namespace mcscf{

/*
** C_DGEMM_12()
**
** This function calculates C(m,n)=alpha*A(k,m)*B(n,k)+ beta*C(m,n)
**
** nra = number of rows in A
** ncb = number of columns in B
** ncc = number of columns in C
*/
void C_DGEMM_12(int m, int n, int k, double alpha,
           double *A, int nra, double *B, int ncb, double beta, double *C,
           int ncc)
{
  //  the only strange thing we need to do is reverse everything
  //  since the stride runs differently in C vs. Fortran

  /* also, do nothing if a dimension is 0 */
  if (m == 0 || n == 0 || k == 0) return;

  F_DGEMM("t","t",&n,&m,&k,&alpha,B,&ncb,A,&nra,&beta,C,&ncc);
}

/*
** C_DGEMM_22()
**
** This function calculates C(m,n)=alpha*A(m,k)*B(n,k)+ beta*C(m,n)
**
** nra = number of columns in A
** ncb = number of columns in B
** ncc = number of columns in C
*/
void C_DGEMM_22(int m, int n, int k, double alpha,
           double *A, int nca, double *B, int ncb, double beta, double *C,
           int ncc)
{
  //  the only strange thing we need to do is reverse everything
  //  since the stride runs differently in C vs. Fortran

  /* also, do nothing if a dimension is 0 */
  if (m == 0 || n == 0 || k == 0) return;

  F_DGEMM("t","n",&n,&m,&k,&alpha,B,&ncb,A,&nca,&beta,C,&ncc);
}


}} /* End Namespaces */
