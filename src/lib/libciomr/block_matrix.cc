/*!
\file
\brief Allocate a blocked (memory-contiguous) 2D matrix of doubles
\ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <strings.h>
#include <psifiles.h>

namespace psi {

/*!
** block_matrix(): Allocate a 2D array of doubles using contiguous memory
**
** Allocates a contiguous block of memory for an array of
** doubles, allocates an array of pointers to the beginning of each row and
** returns the pointer to the first row pointer.  This allows transparent
** 2d-array style access, but keeps memory together such that the matrix 
** could be used in conjunction with FORTRAN matrix routines.
**
** Allocates memory for an n x m matrix and returns a pointer to the
** first row. 
**
** \param n = number of rows (unsigned long to allow large matrices)
** \param m = number of columns (unsigned long to allow large matrices)
**
** Returns: double star pointer to newly allocated matrix
**
** T. Daniel Crawford
** Sometime in 1994
**
** Based on init_matrix() from libciomr
** \ingroup CIOMR
*/

double ** block_matrix(unsigned long int n, unsigned long int m)
{
    double **A=NULL;
    double *B=NULL;
    unsigned long int i;

    if(!m || !n) return(static_cast<double **>(0));

//  if ((A = (double **) malloc(n * (unsigned long int)sizeof(double *)))==NULL) {
    if ((A = new double*[n])==NULL) {
        fprintf(stderr,"block_matrix: trouble allocating memory \n");
        fprintf(stderr,"n = %ld\n",n);
        exit(PSI_RETURN_FAILURE);
    }

//  if ((B = (double *) malloc(m*n * (unsigned long int)sizeof(double)))==NULL) {
    if ((B = new double[n*m])==NULL) {
        fprintf(stderr,"block_matrix: trouble allocating memory \n");
        fprintf(stderr,"m = %ld\n",m);
        exit(PSI_RETURN_FAILURE);
    }

    // bzero is not in the C standard, use memset instead.
    //bzero(B, m*n*(unsigned long int)sizeof(double));
    memset(static_cast<void*>(B), 0, m*n*sizeof(double));
    
    for (i = 0; i < n; i++) {
        A[i] = &(B[i*m]);
    }

    return(A);
}


/*!
** free_block(): Free a block matrix
**
** \param array = pointer to matrix to be freed
**
** Returns: none
**
** \ingroup CIOMR
*/
void free_block(double **array)
{
    if(array == NULL) return;
    delete [] array[0];
    delete [] array;
}

}
