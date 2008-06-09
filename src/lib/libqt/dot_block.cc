/*!
  \file
  \brief Take dot product of two block matrices
  \ingroup QT
*/

namespace psi {
	
/*!
** dot_block(): Find dot product of two block matrices
**
** \param A     = block matrix A
** \param B     = block matrix B 
** \param nrows = number of rows of A and B
** \param ncols = number of columns of A and B
** \param alpha = scale factor by which the dot product is multiplied
**
** Returns: dot product
** \ingroup QT
*/
double dot_block(double **A, double **B, int rows, int cols, double alpha)
{
  register long int i;
  double *a, *b;
  double value;
  long int size;

  size = ((long) rows) * ((long) cols);

  if(!size) return 0.0;

  a = A[0]; b = B[0];

  value = 0.0;
  for(i=0; i < size; i++,a++,b++) value += (*a) * (*b);

  return alpha*value;
}

}

