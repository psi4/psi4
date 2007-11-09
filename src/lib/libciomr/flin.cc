#include "includes.h"

extern "C" {

extern void ludcmp(double **, int, int *, double *);
extern void lubksb(double **, int, int *, double *);

/*!
** \file flin.cc
** \ingroup (CIOMR)
*/ 

/*
** flin(): solves linear equations A * x = b.
**
** \param a   = coefficient matrix
** \param b   = known vectors
** \param in  = dimension of a(in*in)
** \param im  = number of b vectors
** \param det = pointer to hold determinant of matrix a
**
** \ingroup (CIOMR)
*/

void flin(double **a, double *b, int in, int im, double *det)
{
    int i,j,k,*indx;

    indx = (int *) init_array(in);

   ludcmp(a,in,indx,det);

   for (i=0; i < in ; i++) *det *= a[i][i];

   for (j=0; j<im; j++)
      lubksb(a,in,indx,b+j*in);

   free(indx);
}

} /* extern "C" */
