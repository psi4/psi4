/*!
** \file
** \brief Initialize a matrix of doubles
** \ingroup CIOMR
*/

#include <psifiles.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <strings.h>

namespace psi {

  /**
  *  WARNING: Psi 3 init/free_matrix routines deprecated
  *  by Robert Parrish, robparrish@gmail.com
  *
  *  the signature of this method remains the same
  *
  *  June 22, 2010
  **/
/*!
** init_matrix(): Initialize an nxm matrix of doubles and return a pointer to
** the first row.  Note that this does not form a matrix which is
** necessarily contiguous in memory.  Use block_matrix() for that.
**
** \param n = number of rows (unsigned long to allow large matrices)
** \param m = number of columns (unsigned long to allow large matrices)
**
** Returns: pointer to first row
**
** \ingroup CIOMR
*/
double ** init_matrix(unsigned long int n, unsigned long int m)
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
//<<<<<<<<<<<<<<<<<<<<<
// BEGIN DEPRECATED CODE
//<<<<<<<<<<<<<<<<<<<<<<

  /**
  double **array=NULL;
  unsigned long int i;

  if ((array = (double **) malloc(n*(unsigned long int)sizeof(double *)))
    ==NULL) {
    fprintf(stderr,"init_matrix: trouble allocating memory \n");
    fprintf(stderr,"n = %ld\n",n);
    exit(PSI_RETURN_FAILURE);
  }

  for (i = 0; i < n; i++) {
    if ((array[i] = (double *) malloc(m*(unsigned long int)sizeof(double)))
      ==NULL) {
      fprintf(stderr,"init_matrix: trouble allocating memory \n");
      fprintf(stderr,"i = %ld m = %ld\n",i,m);
      exit(PSI_RETURN_FAILURE);
    }
    bzero(array[i],m*(unsigned long int)sizeof(double));
  }
  return(array);
  **/
}


  /**
  *  WARNING: Psi 3 init/free_matrix routines deprecated
  *  by Robert Parrish, robparrish@gmail.com
  *
  *  the signature of this method remains the same
  *
  *  June 22, 2010
  **/
/*!
** free_matrix(): Free a 2D matrix allocated with init_matrix().
**
** \param array = matrix to free
** \param size = number of rows (unsigned long to allow large matrices)
**
** Returns: none
**
** \ingroup CIOMR
*/
void free_matrix(double **array, unsigned long int /*size*/)
{
    if(array == NULL) return;
    delete [] array[0];
    delete [] array;
//<<<<<<<<<<<<<<<<<<<<<
// BEGIN DEPRECATED CODE
//<<<<<<<<<<<<<<<<<<<<<<

  /**
  unsigned long int i;

  for (i=0; i < size ; i++) {
    free(array[i]);
  }

  free(array);
**/
}
}

