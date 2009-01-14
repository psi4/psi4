/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.3  2002/12/06 15:50:32  crawdad
 * Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
 * necessary.  This is new for the PSI3 execution driver.
 * -TDC
 *
 * Revision 1.2  2000/10/13 19:51:22  evaleev
 * Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
 *
 * Revision 1.1.1.1  2000/02/04 22:52:32  evaleev
 * Started PSI 3 repository
 *
 * Revision 1.2  1999/08/17 19:04:18  evaleev
 * Changed the default symmetric orthogonalization to the canonical
 * orthogonalization. Now, if near-linear dependencies in the basis are found,
 * eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
 * left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
 * longer a square matrix. Had to rework some routines in libfile30, and add some.
 * The progrem prints out a warning if near-linear dependencies are found. TRANSQT
 * and a whole bunch of other codes has to be fixed to work with such basis sets.
 *
 * Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 *
 * Revision 1.1  1991/06/15  20:22:42  seidl
 * Initial revision
 * */

static char *rcsid = "$Id: schmit.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void schmit(int all)
{
   int i,j,ij,nn,num_mo;
   int n,m,ncol;
   double *v,**ctmp,vtmp;
   struct symm *s;

   v = (double *) init_array(nsfmax);
   ctmp = (double **) init_matrix(nsfmax,nsfmax);

   for(n=0; n < num_ir ; n++) {
      s= &scf_info[n];
      if(nn=s->num_so) {
	 num_mo = s->num_mo;
         for (i=0; i < nn ; i++)
            for (j=0; j < num_mo ; j++)
               ctmp[j][i] = s->cmat[i][j];

         ncol = s->nclosed + s->nopen;
         if(s->nhalf) ncol++;
	 if(all) ncol = num_mo;
         if(!ncol) continue;
         for(m=0; m < ncol ; m++) {
            v[0]=ctmp[m][0]*s->smat[0];
            for(i=1; i < nn ; i++) {
               for(j=0,vtmp=0.0; j < i ; j++) {
                  ij=ioff[i]+j;
                  vtmp += ctmp[m][j]*s->smat[ij];
                  v[j] += ctmp[m][i]*s->smat[ij];
                  }
               v[i] = vtmp+ctmp[m][i]*s->smat[ioff[i]+j];
               }
            for(i=0,vtmp=0.0; i < nn ; i++) vtmp += v[i]*ctmp[m][i];
            if(!vtmp) {
		exit(PSI_RETURN_FAILURE);
	    }
            if(vtmp < 10.0e-20) vtmp = 10.0e-20;
            vtmp = 1.0/sqrt(vtmp);

            for(i=0; i < nn ; i++) {
               v[i] *= vtmp;
               ctmp[m][i] *= vtmp;
               }

            if(m < ncol-1) {
               for(i=m+1,vtmp=0.0; i < ncol ; i++) {
                  for(j=0,vtmp=0.0; j<nn ;j++) vtmp += v[j]*ctmp[i][j];
                  for(j=0; j < nn ; j++) ctmp[i][j] -=
                                                      vtmp*ctmp[m][j];
                  }
               }
            }

         for (i=0; i < nn ; i++)
            for (j=0; j < num_mo ; j++)
               s->cmat[i][j] = ctmp[j][i];

         }
      }
   free(v); v = NULL;
   free_matrix(ctmp,nsfmax); ctmp = NULL;
   }

}} // namespace psi::cscf
