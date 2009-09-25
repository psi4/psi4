/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi { namespace cphf {

/* cphf_X(): Solve the first-order CPHF equations for a nuclear
** perturbation.
**
** The CPHF equations are:
**
** A U^a = B0^a
**
** where A is the MO hessian [computed in mohess()], U^a is the
** orbital response (CPHF) coefficient, and B0^a is the
** perturbation-dependent inhomogenous factor. B0^a is given by (MO
** basis):
**
** (B0^a)_ai = (F^a)_ai - (S^a)_ai eps_i - (S^a)_jk [ 2 (ai|jk) -
** (aj|ik) ]
**
** where (F^a)_pq is a Fock derivative and (S^a)_pq is an overlap
** derivative.
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.10,
** esp. pp. 128-132.
**
** Also see N.C. Handy, R.D. Amos, J.F. Gaw, J.E. Rice, and
** E.D. Simandiras, "The elimination of singularities in derivative
** calculations," Chem. Phys. Lett. 120, 151 (1985) for an explanation
** of the treatment of dependent-pairs in the U coefficients.
** 
** TDC, December 2001 (revised October 2002)
*/

void cphf_X(double ***S, double **Baijk, double **Aaibj, double ***U)
{
  int coord, coord_a, coord_b, AI, *ipiv, error;
  double **Acopy;
  double *inbuf, value;
  char *label;
  int p, q, k, l, a, i, j, m, pq, kl, pk, qk, ql, pl, pj, jq, jl, ij;
  int pqkl, pkql, plqk, pjql, pqjl, pljq;
  int pfirst, plast, qfirst, qlast, kfirst, klast, lfirst, llast;
  int afirst, alast, ifirst, ilast, jfirst, jlast;
  int psym, qsym, ksym, lsym, asym, isym, jsym;

  ipiv = init_int_array(num_ai);

  /* Solve the CPHF equations */
  Acopy = block_matrix(num_ai, num_ai); /* keep a copy of A */
  memcpy(Acopy[0], Aaibj[0], num_ai*num_ai*sizeof(double));

  for(coord=0; coord < natom*3; coord++) {
    error = C_DGESV(num_ai, 1, &(Aaibj[0][0]), num_ai, &(ipiv[0]), 
            &(Baijk[coord][0]), num_ai);

    /* Recopy A because DGESV corrupts it */
    memcpy(Aaibj[0], Acopy[0], num_ai*num_ai*sizeof(double));
  }

  /* Sort the U matrices to matrix form */
  for(coord=0; coord < natom*3; coord++) {
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

        for(isym=0; isym < nirreps; isym++) {

	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    U[coord][a][i] = Baijk[coord][AI];
          }
        }
      }
    }
  } 

  if (print_lvl > 5) {
    for(coord=0; coord < natom*3; coord++) {
      fprintf(outfile, "\nU[%d] Matrix (Baikj):\n", coord);
      print_mat(U[coord], nmo, nmo, outfile);
    }
  }
  
  /* Add the dependent pairs on the lower triangle */
  for(isym=0; isym < nirreps; isym++) {
    ifirst = ofirst[isym];
    ilast = olast[isym];

    for(i=ifirst; i <= ilast; i++) {

      for(jsym=0; jsym < nirreps; jsym++) {

	jfirst = ofirst[jsym];
	jlast = olast[jsym];

	for(j=jfirst; (j <= jlast) && (j <= i); j++) {

	  if(i==j) continue; /* apply only to non-diagonal terms */

	  /* this is the Handy -1/2 S trick for dependent pairs */
	  for(coord=0; coord < natom*3; coord++) {
            U[coord][i][j] = -0.5 * S[coord][i][j];
	  }
        }
      }
    }
  }

  /* Add the upper triangle */
  for(coord=0; coord < natom*3; coord++) {

    /* dependent pairs */
    for(isym=0; isym < nirreps; isym++) {
      ifirst = ofirst[isym];
      ilast = olast[isym];
      for(i=ifirst; i <= ilast; i++) {

	U[coord][i][i] -= 0.5 * S[coord][i][i];

	for(jsym=0; jsym < nirreps; jsym++) {

	  jfirst = ofirst[jsym];
	  jlast = olast[jsym];

	  for(j=jfirst; (j <= jlast) && (j < i); j++)
	    U[coord][j][i] = -U[coord][i][j] - S[coord][i][j];
	}

      }
    }

    /* independent pairs */
    for(asym=0; asym < nirreps; asym++) {
      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym=0; isym < nirreps; isym++) {
	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++) {
	    U[coord][i][a] = -U[coord][a][i] - S[coord][a][i];
	  }
	}
      }
    }
  }


  /* Write the U matrices out to disk */
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  psio_open(PSIF_CPHF, 1);
  for(coord=0; coord < natom*3; coord++) {
    sprintf(label, "UX(%d)", coord);
    psio_write_entry(PSIF_CPHF, label, (char *) &(U[coord][0][0]), nmo*nmo*sizeof(double));
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
    if(print_lvl > 5) {
      fprintf(outfile, "\nU[%d] Matrix (MO):\n", coord);
      print_mat(U[coord], nmo, nmo, outfile);
    }
  }
  psio_close(PSIF_CPHF, 1);

  free_block(Acopy);
}

}} // namespace psi::cphf
