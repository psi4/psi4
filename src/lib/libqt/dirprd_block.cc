/*!
  \file dirprd_block.c
  \ingroup (QT)
*/

extern "C" {
	
/*!

  dirprd_block()

  This function takes two block matrices A and B and multiplies
  each element of B by the corresponding element of A

  \param **A: block matrix A
  \param **B: block matrix B 
  \param nrows: number of rows of A and B
  \param ncols: number of columns of A and B

  \ingroup (QT)
*/
void dirprd_block(double **A, double **B, int rows, int cols)
{
  register long int i;
  double *a, *b;
  long size;

  size = ((long) rows) * ((long) cols);

  if(!size) return;

  a = A[0]; b= B[0];

  for(i=0; i < size; i++, a++, b++) (*b) = (*a) * (*b);
}

} /* extern "C" */