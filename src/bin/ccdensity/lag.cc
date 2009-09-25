/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* lag(): Build the orbital Lagrangian, I'pq, defined in spin-orbitals
** as:
**
** I'pq = sum_r fpr (Dqr + Drq) + sum_rs <pr||qs> (Drs + Dsr) (q,occ)
**        + sum_rst <pr||st> Gqrst + 2 fpq (q,occ)
**
** The orbital-response component of the gradient (for non-correlated
** orbitals) is defined as
**
** dE/dx <--- sum_pq I'pq U(x)pq
**
** where U(x)pq is the usual CPHF coefficient.  Note, however, that
** the final expression we want involves not CPHF coefficients, but
** overlap derivatives.  For example, in the occupied-occupied and
** virtual-virtual blocks, a choice of non-canonical perturbed
** orbitals allows the assignments
**
** U(x)ij = -1/2 S(x)ij      and      U(x)ab = -1/2 S(x)ab
**
** to be made.  We also choose to incorporate the -1/2 prefactor
** into the Largrangian itself so the final orbital response
** expression will appear as
**
** dE/dx <--- sum_pq Ipq S(x)pq
**
** where Ipq is the "relaxed" Lagrangian (see relax_I.c).
**
** The final set of loops force the appropriate open-shell terms to
** zero for ROHF refs. (See the description of the treatment of
** open-shells in the ROHF-CCSD code as discussed in CCSORT for an
** explanation of why this is necessary.) */

void Iij(struct RHO_Params rho_params);
void Iab(struct RHO_Params rho_params);
void Iai(struct RHO_Params rho_params);
void Iia(struct RHO_Params rho_params);

void lag(struct RHO_Params rho_params)
{
  int h, nirreps, i, j, a, b;
  int *occpi, *virtpi, *openpi;
  dpdfile2 I;
  
  Iij(rho_params);
  Iab(rho_params);
  Iai(rho_params);
  Iia(rho_params);

  /* Multiply all I'pq components by -1/2 for compatibility with the
     final gradient expression */

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 2, 2, "I'ij");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 3, 3, "I'ab");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);
    dpd_file2_init(&I, CC_OEI, 0, 3, 2, "I'ai");
    dpd_file2_scm(&I, -0.5);
    dpd_file2_close(&I);

  }

  /* Now go through all terms involving open-shell orbitals and force
     the appropriate spin cases to zero. */

  if(params.ref == 1) { /** ROHF **/
    nirreps = moinfo.nirreps;
    occpi = moinfo.occpi; virtpi = moinfo.virtpi; openpi = moinfo.openpi;


    dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");
    dpd_file2_mat_init(&I);
    dpd_file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(i=(occpi[h]-openpi[h]); i < occpi[h]; i++) {
	for(j=(occpi[h]-openpi[h]); j < occpi[h]; j++) {
	  I.matrix[h][i][j] = 0.0;
	}
      }
      for(i=(occpi[h]-openpi[h]); i < occpi[h]; i++) {
	for(j=0; j < occpi[h]; j++) {
	  I.matrix[h][i][j] = 0.0;
	}
      }
      for(i=0; i < occpi[h]; i++) {
	for(j=(occpi[h]-openpi[h]); j < occpi[h]; j++) {
	  I.matrix[h][i][j] = 0.0;
	}
      }
    }
    dpd_file2_mat_wrt(&I);
    dpd_file2_mat_close(&I);
    dpd_file2_close(&I);

    dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'AB");
    dpd_file2_mat_init(&I);
    dpd_file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	for(b=(virtpi[h]-openpi[h]); b < virtpi[h]; b++) {
	  I.matrix[h][a][b] = 0.0;
	}
      }
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	for(b=0; b < virtpi[h]; b++) {
	  I.matrix[h][a][b] = 0.0;
	}
      }
      for(a=0; a < virtpi[h]; a++) {
	for(b=(virtpi[h]-openpi[h]); b < virtpi[h]; b++) {
	  I.matrix[h][a][b] = 0.0;
	}
      }
    }
    dpd_file2_mat_wrt(&I);
    dpd_file2_mat_close(&I);
    dpd_file2_close(&I);

    dpd_file2_init(&I, CC_OEI, 0, 1, 1, "I'ab");
    dpd_file2_close(&I);

    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");
    dpd_file2_mat_init(&I);
    dpd_file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	for(i=0; i < occpi[h]; i++) {
	  I.matrix[h][a][i] = 0.0;
	}
      }
    }
    dpd_file2_mat_wrt(&I);
    dpd_file2_mat_close(&I);
    dpd_file2_close(&I);

    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");
    dpd_file2_mat_init(&I);
    dpd_file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(a=0; a < virtpi[h]; a++) {
	for(i=(occpi[h] - openpi[h]); i < occpi[h]; i++) {
	  I.matrix[h][a][i] = 0.0;
	}
      }
    }
    dpd_file2_mat_wrt(&I);
    dpd_file2_mat_close(&I);
    dpd_file2_close(&I);

    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
    dpd_file2_mat_init(&I);
    dpd_file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
	for(a=(virtpi[h] - openpi[h]); a < virtpi[h]; a++) {
	  I.matrix[h][i][a] = 0.0;
	}
      }
    }
    dpd_file2_mat_wrt(&I);
    dpd_file2_mat_close(&I);
    dpd_file2_close(&I);

    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");
    dpd_file2_mat_init(&I);
    dpd_file2_mat_rd(&I);
    for(h=0; h < nirreps; h++) {
      for(i=(occpi[h] - openpi[h]); i < occpi[h]; i++) {
	for(a=0; a < virtpi[h]; a++) {
	  I.matrix[h][i][a] = 0.0;
	}
      }
    }
    dpd_file2_mat_wrt(&I);
    dpd_file2_mat_close(&I);
    dpd_file2_close(&I);

  }
}

}} // namespace psi::ccdensity
