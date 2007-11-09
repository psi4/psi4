/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.2  2001/06/29 20:39:29  evaleev
 * Modified cscf to use libpsio to store supermatrix files.
 *
/* Revision 1.1.1.1  2000/02/04 22:52:30  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.3  1999/10/22 19:47:18  evaleev
/* A direct SCF-enabled version (set DIRECT_SCF=TRUE in input.dat).
/*
/* Revision 1.2  1999/08/17 19:04:15  evaleev
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

static char *rcsid = "$Id: formgc.cc 3592 2007-09-28 13:01:33Z evaleev $";

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

static double *gtmp,*ptmp;
extern struct c_pkints {
         int ij;
         int kl;
         double pval;
         } *c_outbuf;
extern int last;
static int where=0;
static int *int_nums;

void formg_closed()

{
   register int i,j,k,joff,nn;
   register int ij,kl;
   int ilast,num;
   int tmpsiz;
   double tmpval;
   struct c_pkints *ctmp;

   tmpsiz = ioff[nbasis];

   if(gtmp == NULL) {
      gtmp = (double *) init_array(tmpsiz);
      ptmp = (double *) init_array(tmpsiz);
      }
   else bzero(gtmp,sizeof(double)*tmpsiz);

   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++)
            for(j=0; j <= i ; j++)
               ptmp[ioff[i+joff]+j+joff] = scf_info[k].dpmat[ioff[i]+j];
         joff += nn;
         }
      }

   if(!where) {
      /* int_nums = (int *) init_array(num_bufs+1); */
      int_nums = (int *) init_int_array(num_bufs+1);
      for(i=1; i < num_bufs ; i++) int_nums[i]=maxbuf;
      int_nums[num_bufs]=last;
      where=num_bufs;
      }

   num=int_nums[where];
   for (j=0; j < num_bufs; j++) {
      ctmp = c_outbuf;

      for (i=num; i ; i--,ctmp++) {
         ij = (*ctmp).ij;
         kl = (*ctmp).kl;
         tmpval = (*ctmp).pval;

         gtmp[ij] += ptmp[kl]*tmpval;
         gtmp[kl] += ptmp[ij]*tmpval;
         }

      if (readflg && j < num_bufs-1) {
         if(where==num_bufs) {
            where=0;
	    Pmat.bufpos = PSIO_ZERO;
            }
         where++;
         num=int_nums[where];
	 psio_read(Pmat.unit, Pmat.key, (char *) c_outbuf, sizeof(struct c_pkints)*num,
		   Pmat.bufpos, &(Pmat.bufpos));
         }
      }
   for(k=joff=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
         for(i=0; i < nn ; i++)
            for(j=0; j <= i ; j++)
               scf_info[k].gmat[ioff[i]+j] += gtmp[ioff[i+joff]+j+joff];
         joff += nn;
         }
      }
}

}} // namespace psi::cscf
