/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <physconst.h>
#define EXTERN
#include "globals.h"

namespace psi { namespace cphf {

/* polarize(): Compute the total polarizability from the dipole moment
** integrals computed by "cints --oeprop" and the field-perturbation
** CPHF coefficients computed in cphf_F().  The equation for the
** polarizability is (MO basis):
**
** alpha_fg = - 4 (U^g)_ai (mu^f)_ai
**
** where (U^g)_pq is a CPHF coefficient and (mu^f)_pq is a dipole
** moment integral.  Summation over all indices is implied.
**
** In addition, the dipole moment is computed, just because we have
** the integrals and it's easy:
**
** mu^f = (mu^f)_elec + (mu^f)_nuc = 2 (mu^f)_ii + R_Af Z_A
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.17,
** esp. pp. 314-317.
**
** TDC, October 2002
*/

void polarize(double ***UF)
{
  double **alpha;
  double ***MU, *mu;
  double *scratch;
  int stat, f, g, j, ij;
  int a, asym, afirst, alast;
  int i, isym, ifirst, ilast;

  MU = (double ***) malloc(3 * sizeof(double **));
  for(i=0; i < 3; i++) MU[i] = block_matrix(nmo,nmo);

  /* Grab the MO-basis dipole integrals from disk */
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_MO_MX, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      MU[0][i][j] = MU[0][j][i] = scratch[ij];
  zero_arr(scratch,ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_MO_MY, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      MU[1][i][j] = MU[1][j][i] = scratch[ij];
  zero_arr(scratch,ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_MO_MZ, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      MU[2][i][j] = MU[2][j][i] = scratch[ij];
  if(print_lvl & 2) {
    fprintf(outfile, "MO Mu-X Ints:\n");
    print_mat(MU[0], nmo, nmo, outfile);
    fprintf(outfile, "MO Mu-Y Ints:\n");
    print_mat(MU[1], nmo, nmo, outfile);
    fprintf(outfile, "MO Mu-Z Ints:\n");
    print_mat(MU[2], nmo, nmo, outfile);
  }
  free(scratch);

  /* Compute the dipole moment */
  mu = init_array(3);
  for(f=0; f < 3; f++) {
    /* electronic contribution */
    for(isym=0; isym < nirreps; isym++) {
      ifirst = ofirst[isym];
      ilast = olast[isym];
      for(i=ifirst; i <= ilast; i++) 
	mu[f] += 2.0 * MU[f][i][i];
    }
    /* nuclear contribution */
    for(a=0; a < natom; a++) {
      mu[f] += geom[a][f] * zvals[a];
    }
  }
 
  fprintf(outfile, "\n\tHartree-Fock Electric Dipole Moments:\n");
  fprintf(outfile, "\t------------------------------------\n");
  for(f=0; f < 3; f++) {
    fprintf(outfile, "\tmu[%d] = %10.6f (e a0) = %10.6f (debye)\n", 
            f, mu[f], mu[f]*_dipmom_au2debye);
  }
  fprintf(outfile, "\n\tTotal dipole moment (debye) = %20.10f\n", 
          sqrt(mu[0]*mu[0] + mu[1]*mu[1] + mu[2]*mu[2])*_dipmom_au2debye);
  
  free(mu);

  /* Compute the polarizability tensor */
  alpha = block_matrix(3,3);
  for(f=0; f < 3; f++) {
    for(g=0; g < 3; g++) {
      for(asym=0; asym < nirreps; asym++) {
	afirst = vfirst[asym];
	alast = vlast[asym];
	for(a=afirst; a <= alast; a++) {
	  for(isym=0; isym < nirreps; isym++) {
	    ifirst = ofirst[isym];
	    ilast = olast[isym];
	    for(i=ifirst; i <= ilast; i++) {
	      alpha[f][g] -= 4.0 * UF[g][a][i] * MU[f][a][i];
	    }
	  }
	}
      }
    }
  }
  
  fprintf(outfile, "\n\tHartree-Fock Electric Polarizability Tensor [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "\t---------------------------------------------------------------\n");
  print_mat(alpha,3,3,outfile);

  for(i=0; i < 3; i++) free_block(MU[i]);
  free(MU);
  free_block(alpha);
}

}} // namespace psi::cphf
