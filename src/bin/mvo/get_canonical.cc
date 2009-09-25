/*! \file
    \ingroup MVO
    \brief Enter brief description of file here 
*/
/*
**  get_canonical.c
**  
**  Code to rediagonalize the Fock matrix in each RAS subspace to
**  obtain (semi) canonical orbitals again
**
**  C. David Sherrill
**  Georgia Institute of Technology
**  January 2003
**
** Let me try to explain what this is for.  I forsee two uses:
**  (1) Improve convergence of CI based on NO's
**  (2) Allow proper computation of (T) correction in NO basis
**
** The code will attempt to re-diagonalize the Fock matrix in each
** "RAS" space defined by moinfo.ras_opi.  This makes it nearly arbitrary.
** However, one must be careful how one uses this routine.
**
** For the purpose of improving CI convergence, this probably works 
** ok, or not, but either way nothing is truly wrong.
**
** For the purpose of computing the (T) correction, one must be very
** careful.  The (T) correction assumes the Fock matrix is diagonal.
** It will *not* be diagonal if we blindly diagonalize it in a series
** of MO subspaces!  This routine simply ignores the coupling of the
** Fock matrix blocks among each other, which is surely wrong and does
** not truly make the Fock matrix diagonal!  Thus, we must make sure
** that the *full* Fock matrix of all active orbitals is indeed diagonal.
** This can be accomplished by calling this routine with only two "RAS"
** spaces: (a) All active orbitals, (b) the frozen virtual orbitals. 
** This ignores coupling between active and frozen virtuals, but since
** the virtuals are deleted, any coupling is ignored, too!  (Should be...)
**
** Revisions:
**  15 May 2006 by C. David Sherrill 
**  - Cleaned up this routine and added ability to handle restricted orbitals
**  - Tested on canonical SCF orbitals and appear to just get them back,
**    which is what should happen
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "params.h"
#include "globals.h"


/* First definitions of globals */
extern "C" {
  extern FILE *infile, *outfile;
}

namespace psi { namespace mvo {

extern int *ioff;
extern struct MOInfo moinfo;
extern struct Params params;

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))


void get_canonical(void)
{
  int ras_space, irrep, size, ntri;
  int i,j,ij,i2,j2,offset;
  double *evals, *new_evals;
  double **evecs;
  double **C, **C_new, **C_sub, **C_sub_new;
  double **F, *F_sub;
  double **tras;
  void form_fock_full(double **);

  /* Form the full Fock matrix */
  F = block_matrix(moinfo.nmo,moinfo.nmo);
  form_fock_full(F);
  if (params.print_lvl) {
    fprintf(outfile, "Fock Matrix (MO):\n");
    print_mat(F, moinfo.nmo, moinfo.nmo, outfile);
  }

  /* Read current MO's from checkpoint file into C matrix */
  chkpt_init(PSIO_OPEN_OLD);
  C = chkpt_rd_scf();
  chkpt_close();

  if (params.print_lvl) {
    fprintf(outfile, "\nFull C matrix\n");
    print_mat(C, moinfo.nso, moinfo.nmo, outfile);
  } 

  /* allocate space */
  C_new = block_matrix(moinfo.nso,moinfo.nmo);

  /* go ahead and allocate subblocks big enough for any case */
  C_sub = block_matrix(moinfo.nso,moinfo.nmo);
  C_sub_new = block_matrix(moinfo.nso,moinfo.nmo);
  ntri = (moinfo.nmo*(moinfo.nmo+1))/2;
  F_sub = init_array(ntri);
  evals = init_array(moinfo.nmo);
  evecs = block_matrix(moinfo.nmo,moinfo.nmo);
  new_evals = init_array(moinfo.nmo);

  /* put the orbs per subspace into a big array "tras".
     frozen docc, then restricted docc, then the RAS spaces,
     then restricted uocc, then frozen uocc */
  tras = block_matrix(MAX_RAS_SPACES+3,moinfo.nirreps);
  for (i=0; i<moinfo.nirreps; i++) {
    tras[0][i] = params.fzc ? 0 : moinfo.frdocc[i];
    tras[1][i] = moinfo.rstrdocc[i];
    tras[2+MAX_RAS_SPACES][i] = moinfo.rstruocc[i];
    /* don't think we need frozen uocc actually */
  }
  for (ras_space=0; ras_space<MAX_RAS_SPACES; ras_space++) {
    for (i=0; i<moinfo.nirreps; i++) {
      tras[2+ras_space][i] = moinfo.ras_opi[ras_space][i];
    }
  }

  for (ras_space=0,offset=0; ras_space<MAX_RAS_SPACES+3; ras_space++) {
    for (i=0,j=0; i<moinfo.nirreps; i++) {
      j += tras[ras_space][i];
    }
    if (j<1) continue; /* no orbitals in this ras space */
    for (irrep=0; irrep<moinfo.nirreps; 
		  offset+=tras[ras_space][irrep],irrep++) {
      size = tras[ras_space][irrep];
      if (size < 1) continue;

      /* copy the sub-block of the full Fock matrix into F_sub */
      ntri = (size * (size+1))/2;
      for (i=0,ij=0; i<size; i++) {
        for (j=0; j<=i; j++,ij++) {
          F_sub[ij] = F[i+offset][j+offset];
	}
      } 
      if (params.print_lvl) {
        fprintf(outfile, "\nOld F_sub matrix for subspace %d, irrep %d\n",
                ras_space, irrep);
	print_array(F_sub, size, outfile);
      }
      rsp(size, size, ntri, F_sub, evals, 1, evecs, 1E-12);
      if (params.print_lvl) {
        fprintf(outfile, "\nDiagonalization of matrix F_sub\n");
	eivout(evecs, evals, size, size, outfile);
      }
      /* place the eigenvalues where they go */
      for (j=0; j<size; j++) {
        j2 = moinfo.corr2pitz[j+offset];
        new_evals[j2] = evals[j];
      }

      /* load up the corresponding block of the SCF matrix C */
      for (i=0; i<moinfo.sopi[irrep]; i++) {
        i2 = moinfo.first_so[irrep] + i;
        for (j=0; j<size; j++) {
          j2 = moinfo.corr2pitz[j+offset]; 
          C_sub[i][j] = C[i2][j2];
	}
      }
      if (params.print_lvl) {
        fprintf(outfile, "Old C_sub for subspace %d, irrep %d\n", ras_space,
                irrep);
	print_mat(C_sub, moinfo.sopi[irrep], size, outfile);
      }
      /* multiply by the Fock eigenvectors to get new C_sub */
      mmult(C_sub, 0, evecs, 0, C_sub_new, 0, moinfo.sopi[irrep], size,
            size, 0);

      if (params.print_lvl) {
        fprintf(outfile, "C_sub_new for subspace %d, irrep %d\n", ras_space,
                irrep);
	print_mat(C_sub_new, moinfo.sopi[irrep], size, outfile);
      }

      /* load up this contribution to the new C matrix */
      for (i=0; i<moinfo.sopi[irrep]; i++) {
        i2 = moinfo.first_so[irrep] + i;
        for (j=0; j<size; j++) {
          j2 = moinfo.corr2pitz[j+offset]; 
	  C_new[i2][j2] = C_sub_new[i][j];
	}
      }

    } /* end loop over irreps */
  } /* end loop over RAS spaces */

  /* write out the new C matrix */
  if (params.print_lvl) {
    fprintf(outfile, "\nFull new C matrix\n");
    print_mat(C_new, moinfo.nso, moinfo.nmo, outfile);
  } 

  if (params.print_lvl) {
    fprintf(outfile, "\nNew orbital eigenvalues\n");
    for (i=0; i<moinfo.nmo; i++) {
      fprintf(outfile, "%12.6lf  ", new_evals[i]);
      if ((i+1)%5 == 0) fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_scf(C_new);
  chkpt_wt_evals(new_evals);
  chkpt_close();

  free_block(F);
  free(F_sub);
  free(evals);
  free(new_evals);
  free_block(evecs);
  free_block(C_sub);
  free_block(C_sub_new);
  free_block(C_new);
  free_block(C);
  free_block(tras);
}

}} // end namespace psi::mvo

