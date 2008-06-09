/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.1  2000/02/04 22:52:29  evaleev
 * Initial revision
 *
/* Revision 1.2  1999/08/17 19:04:14  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: dmat_2.cc 3815 2008-02-13 21:50:07Z sherrill $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void dmat_2(int opblk)
{
   int i,j,k,l,ij,n;
   int ndocc;
   double ptempc,ctmp;
   struct symm *s;

   for (l=0; l < num_ir ; l++) {
      s = &scf_info[l];
      if (n=s->num_so) {
         ndocc = s->nclosed;
         if (l == opblk) ndocc++;
         
         ij=0;
         for (i=0; i < n ; i++) {
            for (j=0; j < i; j++,ij++) {
               ptempc=0.0;
               for (k=0; k < ndocc ; k++) {
                  ptempc += 4.0*s->cmat[i][k]*s->cmat[j][k];
                  }
               if(opblk==opblk1) s->pmato2[ij] = ptempc;
               else s->pmat2[ij] = ptempc;
               }
            ptempc=0.0;
            for (k=0; k < ndocc ; k++) {
               ctmp=s->cmat[i][k];
               ptempc += 2.0*ctmp*ctmp;
               }
            if(opblk==opblk1) s->pmato2[ij] = ptempc;
            else s->pmat2[ij] = ptempc;
            ij++;
            }
         if(print & 4) {
            if(opblk==opblk1) {
               fprintf(outfile,
                       "\ndensity matrix 1 for irrep %s",s->irrep_label);
               print_array(s->pmato2,n,outfile);
               }
            else {
               fprintf(outfile,
                       "\ndensity matrix 1 for irrep %s",s->irrep_label);
               print_array(s->pmat2,n,outfile);
               }
            }
         }
      }
   }

}} // namespace psi::cscf
