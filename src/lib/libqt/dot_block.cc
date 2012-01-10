/*!
  \file
  \brief Take dot product of two block matrices
  \ingroup QT
*/

#include <libqt/qt.h>

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
    double value;
    long int size;
        
    size = ((long) rows) * ((long) cols);
        
    if(!size) return 0.0;
        
    C_DGEMM('T', 'N', 1, 1, size, alpha, A[0], 1, B[0], 1, 0.0, &value, 1);
        
    return value;
}

}

