/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.1  2000/02/04 22:52:32  evaleev
 * Initial revision
 *
 * Revision 1.2  1999/08/17 19:04:17  evaleev
 * Changed the default symmetric orthogonalization to the canonical
 * orthogonalization. Now, if near-linear dependencies in the basis are found,
 * eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
 * left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
 * longer a square matrix. Had to rework some routines in libfile30, and add some.
 * The progrem prints out a warning if near-linear dependencies are found. TRANSQT
 * and a whole bunch of other codes has to be fixed to work with such basis sets.
 *
 * Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 * */

static char *rcsid = "$Id: rotate_vector.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void check_rot(int nn, int num_mo, double **cold, double **cnew, 
  double *smat_pac, double *fock_evals, int irrep);
extern void formg_two(int iju, int* optest);
extern void formg_open();

void rotate_vector()
{
   int i,j,ij,k,nn,num_mo;
   int nc,no;
   int joff,oj;
   int *optest,iju=0;
   double **scr1,**scr2,**scr3;
   double tol=1.0e-15;
   struct symm *s;

   scr1 = (double **) init_matrix(nsfmax,nsfmax);
   scr2 = (double **) init_matrix(nsfmax,nsfmax);
   scr3 = (double **) init_matrix(nsfmax,nsfmax);

   dmat();

   if(twocon) {
      /* optest = (int *) init_array(ioff[nbasis]/2); */
      optest = (int *) init_int_array(ioff[nbasis]);

/* find open shells */

      for (i=0; i < num_ir ; i++) {
         iju += scf_info[i].num_so;
         if(scf_info[i].nopen) {
            iju = ioff[iju];
            break;
            }
         }

/* set up array of flags indicating open shells */

      for (k=0,joff=0; k < num_ir ; k++) {
         s = &scf_info[k];
         if (nn=s->num_so) {
            for (i=0; i < nn ; i++)
               for (j=0; j <= i ; j++)
                  if(s->nopen) optest[ioff[i+joff]+j+joff] = 1;
            joff += nn;
            }
         }

      formg_two(iju,optest);
      }
   else formg_open();

   for(i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if(nn=s->num_so) {
	 num_mo = s->num_mo;
         nc=s->nclosed;
         no=s->nopen;
         add_arr(s->hmat,s->gmat,s->fock_pac,ioff[nn]);

         tri_to_sq(s->fock_pac,scr1,nn);
         mmult(s->cmat,1,scr1,0,scr2,0,num_mo,nn,nn,0);
         mmult(scr2,0,s->cmat,0,scr1,0,num_mo,nn,num_mo,0);

         zero_mat(scr2,nn,nn);
         for(j=0; j < nc ; j++)
            for(k=0; k < nc ; k++)
               scr2[j][k]=scr1[j][k];
         for(j=nc; j < nc+no ; j++)
            for(k=nc; k < nc+no ; k++)
               scr2[j][k]=scr1[j][k];
         for(j=nc+no; j < num_mo ; j++)
            for(k=nc+no; k < num_mo ; k++)
               scr2[j][k]=scr1[j][k];

         sq_rsp(num_mo,num_mo,scr2,s->fock_evals,1,scr3,tol);

         mmult(s->cmat,0,scr3,0,scr2,0,nn,num_mo,num_mo,0);

         if (icheck_rot)
	   check_rot(nn, num_mo, s->cmat, scr2, s->smat, s->fock_evals, i);

         for(j=0; j < nn ; j++)
            for(k=0; k < num_mo ; k++)
               s->cmat[j][k]=scr2[j][k];

         }
      }
   free_matrix(scr1,nsfmax);
   free_matrix(scr2,nsfmax);
   free_matrix(scr3,nsfmax);
   scr1 = scr2 = scr3 = NULL;
   }

}} // namespace psi::cscf
