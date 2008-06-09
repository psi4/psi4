/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.1  2000/02/04 22:52:32  evaleev
 * Initial revision
 *
/* Revision 1.2  1999/08/17 19:04:18  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 */

static char *rcsid = "$Id: sdot.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void sdot(double** a, double** b, int n, double* value)
   {
      register int i,j;
      double *ta, *tb, tval;

      tval = 0.0;
      for (i=0; i < n; i++) {
         ta = a[i];
         tb = b[i];
         for (j=0; j <= i; j++,ta++,tb++) {
            tval += (*ta) * (*tb);
            }
         }
      *value = tval;
      }

void vdot(double* a, double* b, int n, double* value)
{
  int i;
  double tval=0.0;

  for(i=0; i < n ; i++) tval += a[i]*b[i];

  *value = tval;
  }

}} // namespace psi::cscf
