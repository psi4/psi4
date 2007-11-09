/*!
  \file fill_sym_matrix.c
  \ingroup (QT)
*/

extern "C" {
	
/*!
** fill_sym_matrix():
** This function fills a symmetric matrix by placing the elements of 
** the lower triangle into the upper triangle.
**
** \param double A: matrix to symmetrize
** \param int size: number of rows or columns (assume square)
**
** Returns: 
**    none
** \ingroup (QT)
*/
void fill_sym_matrix(double **A, int size)
{
   double **row, *col; 
   int rc, cc;

   row = A ;
   for (rc = 0; rc < (size-1); rc++) {
      col = *row;
      for (cc = 0; cc < size; cc++) {
         if (cc > rc) {
            *col = A[cc][rc] ;
            }
         col++ ;
         }
      row++ ;
      }
}

} /* extern "C" */