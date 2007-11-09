/*! \file 
    \ingroup (CSCF)
    \brief Enter brief description of file here 
*/
/* $Log$
 * Revision 1.5  2007/04/05 15:45:25  crawdad
 * Fixed a few memory leaks identified by valgrind. -TDC
 *
/* Revision 1.4  2004/05/03 04:32:40  crawdad
/* Major mods based on merge with stable psi-3-2-1 release.  Note that this
/* version has not been fully tested and some scf-optn test cases do not run
/* correctly beccause of changes in mid-March 2004 to optking.
/* -TDC
/*
/* Revision 1.3.6.4  2004/04/21 15:45:07  evaleev
/* Modified DIIS algorithm for RHF and ROHF to work in OSO basis rather than in
/* AO basis, to avoid difficulties of transforming between MO and AO bases
/* when linear dependencies are present.
/*
/* Revision 1.3.6.3  2004/04/10 19:41:32  crawdad
/* Fixed the DIIS code for UHF cases.  The new version uses the Pulay scheme of
/* building the error vector in the AO basis as FDS-SDF, followed by xformation
/* into the orthogonal AO basis.   This code converges faster for test cases
/* like cc8, but fails for linearly dependent basis sets for unknown reasons.
/* -TDC
/*
/* Revision 1.3.6.2  2004/04/09 00:16:29  evaleev
/* Added messages explaining why DGETRF and DGETRI most likely fail.
/*
/* Revision 1.3.6.1  2004/04/06 21:29:05  crawdad
/* Corrections to the RHF/ROHF DIIS algorithm, which was simply incorrect.
/* The backtransformation of the DIIS error vectors to the AO basis was not
/* mathematically right.
/* -TDC and EFV
/*
/* Revision 1.3  2002/11/24 22:52:17  crawdad
/* Merging the gbye-file30 branch into the main trunk.
/* -TDC
/*
/* Revision 1.2.6.1  2002/11/23 21:54:45  crawdad
/* Removal of mxcoef stuff for chkpt runs.
/* -TDC
/*
/* Revision 1.2  2000/10/13 19:51:22  evaleev
/* Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
/*
/* Revision 1.1.1.1  2000/02/04 22:52:33  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.3  1999/08/17 19:04:18  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.2  1999/08/11 18:39:03  evaleev
/* Added some checks on the lowest eigenvalue of the overlap matrix.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id: shalf.cc 3592 2007-09-28 13:01:33Z evaleev $";

/* construct S-1/2 matrix 'sahalf' using CANONICAL orthogonalization  */

#define EXTERN
#include "includes.h"
#include "common.h"

namespace psi { namespace cscf {

void shalf(void)
{
  int i,nn,num_mo;
  int ii,jj,kk;
  double *eig_vals, **eig_vecs;
  double **shalf;
  double tol = 1.0e-20;
  double min_eval = 100000.0;
  struct symm *s;
  int info, lwork, *ipiv;
  double *work, **P, **T, **X;
     
  eig_vals = (double *) init_array(nsfmax);
  eig_vecs = (double **) init_matrix(nsfmax,nsfmax);
  shalf = block_matrix(nsfmax,nsfmax);

#if !USE_LIBCHKPT
  mxcoef = 0;
#endif
  mxcoef2 = 0;
  nmo = 0;
     
  /*  diagonalize smat to get eigenvalues and eigenvectors  */

  for (i=0; i < num_ir ; i++) {
    s = &scf_info[i];
    if (nn=s->num_so) {

      rsp(nn,nn,ioff[nn],s->smat,eig_vals,1,eig_vecs,tol);

      /*--- Find the lowest eigenvalue ---*/
      for(ii=0; ii < nn ; ii++)
	if (eig_vals[ii] < min_eval)
	  min_eval = eig_vals[ii];

      if(print & 64) {
	fprintf(outfile,"\noverlap eigenstuff\n");
	eivout(eig_vecs,eig_vals,nn,nn,outfile);
      }

      /*--- Go through the eigenvalues and "throw away"
	the dangerously small ones ---*/
      num_mo = 0;
      for(ii=0; ii < nn ; ii++)
	if (eig_vals[ii] >= LINDEP_CUTOFF) {
	  eig_vals[ii] = 1.0/sqrt(eig_vals[ii]);
	  num_mo++;
	}
      s->num_mo = num_mo;
#if !USE_LIBCHKPT
      mxcoef += num_mo * nn;
#endif
      mxcoef2 += ioff[nn];
      nmo += num_mo;
      if (num_mo < nn)
	fprintf(outfile,"\n  In symblk %d %d eigenvectors of S with eigenvalues < %lf are thrown away",
		i,nn-num_mo,LINDEP_CUTOFF);

      /* form 'sahalf' matrix sahalf = U*(s-1/2)  */

      for (ii=0; ii < nn ; ii++)
	for (jj=0; jj < num_mo ; jj++)
	  s->sahalf[ii][jj] += eig_vecs[ii][jj+nn-num_mo]*eig_vals[jj+nn-num_mo];

      /* form 'shalf' matrix shalf = U*(s^1/2)  */

      for (ii=0; ii < nn ; ii++)
	for (jj=0; jj < num_mo ; jj++)
	  shalf[ii][jj] = eig_vecs[ii][jj+nn-num_mo]*(1.0/eig_vals[jj+nn-num_mo]);

      /*
      fprintf(outfile, "S^-1/2 Original:\n");
      print_mat(s->sahalf, nn, num_mo, outfile);
      */

      /* new transformation matrix for the DIIS error vectors: P^-1 = S^{1/2} * (S^{1/2})^+ */
      C_DGEMM('n', 't', nn, nn, num_mo, 1.0, shalf[0], nsfmax, shalf[0], nsfmax,
	      0.0, s->pinv[0], nn);

    }
  }

  free(eig_vals);
  free_matrix(eig_vecs,nsfmax);

  fprintf(outfile,"\n  The lowest eigenvalue of the overlap matrix was %e\n\n",
	  min_eval);
   
}

}} // namespace psi::cscf
