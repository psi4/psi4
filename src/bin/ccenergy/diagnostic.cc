/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

double diagnostic(void)
{
  int h, nirreps, Gi, Ga;
  int i, a, I, A, row, col;
  int num_elec, num_elec_a, num_elec_b;
  int *occpi, *virtpi;
  int *occ_sym, *vir_sym;
  int *clsdpi, *uoccpi;
  int *openpi;
  double t1diag, t1diag_a, t1diag_b;
  dpdfile2 T1A, T1B;

  nirreps = moinfo.nirreps;
  clsdpi = moinfo.clsdpi; 
  uoccpi = moinfo.uoccpi;
  openpi = moinfo.openpi;

  if(params.ref != 2) {
    occpi = moinfo.occpi; virtpi = moinfo.virtpi;
    occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  }

  /* Compute the number of electrons */
  for(h=0,num_elec_a=0,num_elec_b=0; h < nirreps; h++) {
    num_elec_a += clsdpi[h] + openpi[h];
    num_elec_b += clsdpi[h];
  }
  num_elec = num_elec_a + num_elec_b;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
    t1diag = dpd_file2_dot_self(&T1A);
    dpd_file2_close(&T1A);

    t1diag /= num_elec;
    t1diag = sqrt(t1diag);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1A);
    dpd_file2_mat_rd(&T1A);
    dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_mat_init(&T1B);
    dpd_file2_mat_rd(&T1B);

    t1diag = 0.0;
    for(h=0; h < nirreps; h++) {

      for(i=0; i < (occpi[h] - openpi[h]); i++) {
	for(a=0; a < (virtpi[h] - openpi[h]); a++) {

	  t1diag += (T1A.matrix[h][i][a] + T1B.matrix[h][i][a]) *
	    (T1A.matrix[h][i][a] + T1B.matrix[h][i][a]);
	}
      }

      for(i=0; i < (occpi[h] - openpi[h]); i++) {
	for(a=0; a < openpi[h]; a++) {
	  A = a + uoccpi[h];

	  t1diag += 2 * T1B.matrix[h][i][A] * T1B.matrix[h][i][A];
	}
      }

      for(i=0; i < openpi[h]; i++) {
	I = i + clsdpi[h];
	for(a=0; a < (virtpi[h] - openpi[h]); a++) {
		  
	  t1diag += 2 * T1A.matrix[h][I][a] * T1A.matrix[h][I][a];
	}
      }
    }

    t1diag /= num_elec;
    t1diag = sqrt(t1diag);
    t1diag *= 0.5;

    dpd_file2_mat_close(&T1A);
    dpd_file2_close(&T1A);
    dpd_file2_mat_close(&T1B);
    dpd_file2_close(&T1B);

  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1A);
    dpd_file2_mat_rd(&T1A);
    dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&T1B);
    dpd_file2_mat_rd(&T1B);

    t1diag_a = 0.0;
    t1diag_b = 0.0;
    for(h=0; h < nirreps; h++) {

      for(row=0; row < T1A.params->rowtot[h]; row++) 
	for(col=0; col < T1A.params->coltot[h]; col++) 
	  t1diag_a += (T1A.matrix[h][row][col] * T1A.matrix[h][row][col]);

      for(row=0; row < T1B.params->rowtot[h]; row++) 
	for(col=0; col < T1B.params->coltot[h]; col++) 
	  t1diag_b += (T1B.matrix[h][row][col] * T1B.matrix[h][row][col]);
    }


    t1diag = sqrt((t1diag_a + t1diag_b)/(num_elec_a + num_elec_b));

    dpd_file2_mat_close(&T1A);
    dpd_file2_mat_close(&T1B);
    dpd_file2_close(&T1A);
    dpd_file2_close(&T1B);
	      
  }

  return t1diag; 
}
}} // namespace psi::ccenergy
