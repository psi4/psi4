/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/* Computes the D1 diagnostic as defined in:
 *
 * C.L. Janssen and I.M.B. Nielsen, Chem. Phys. Lett. 290, 423 (1998). [RHF
 * version]
 *
 * I.M.B. Nielsen and C.L. Janssen, Chem. Phys. Lett. 310, 568 (1999).
 *
 * M.L. Leininger, I.M.B. Nielsen, T.D. Crawford, and C.L. Janssen, 
 *   Chem. Phys. Lett. 328, 431-436 (2000).  [ROHF version]
 *
 * */

double d1diag_t1_rhf(void)
{
  int h, nirreps, i;
  double **T, **C, *E, max;
  dpdfile2 T1;

  nirreps = moinfo.nirreps;
  max = 0.0;

  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);

  for(h=0; h < nirreps; h++) {
    if(T1.params->rowtot[h]) {
      T = block_matrix(T1.params->rowtot[h], T1.params->rowtot[h]);

      if(T1.params->rowtot[h] && T1.params->coltot[h]) {
	C_DGEMM('n','t',T1.params->rowtot[h],T1.params->rowtot[h],
		T1.params->coltot[h],1.0,T1.matrix[h][0],T1.params->coltot[h],
		T1.matrix[h][0],T1.params->coltot[h], 0.0,
		T[0], T1.params->rowtot[h]);
      }
      /*
	newmm(T1.matrix[h], 0, T1.matrix[h], 1, T, T1.params->rowtot[h],
	T1.params->coltot[h], T1.params->rowtot[h], 1.0, 0.0);
      */

      E = init_array(T1.params->rowtot[h]);
      C = block_matrix(T1.params->rowtot[h], T1.params->rowtot[h]);
      sq_rsp(T1.params->rowtot[h], T1.params->rowtot[h], T, E, 0, C, 1e-12);

      /* Find maximum eigenvalue of T */
      for(i=0; i < T1.params->rowtot[h]; i++) if(E[i] > max) max = E[i];
	     
      free_block(T);
      free_block(C);
      free(E);
    }
  }

  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);

  max = sqrt(max);

  return max;
}

static double
d1diag_subblock(double **Tave, int row0, int rown, int col0, int coln)
{
  int i,j;
  int nrow = rown - row0;
  int ncol = coln - col0;
  double max = 0.;
  double **Tsub;
  double **Tsq;
  double *E, **C;

  if (nrow && ncol) {
      Tsub  = block_matrix(nrow, ncol);
      Tsq  = block_matrix(nrow, nrow);

      for (i=row0; i<rown; i++) {
          for (j=col0; j<coln; j++) {
              Tsub[i-row0][j-col0] = Tave[i][j];
            }
        }

      C_DGEMM('n','t',nrow,nrow,ncol,1.0,Tsub[0],ncol,Tsub[0],ncol,
	  0.0,Tsq[0],nrow);
      /* newmm(Tsub, 0, Tsub, 1, Tsq, nrow, ncol, nrow, 1.0, 0.0); */

      E = init_array(nrow);
      C = block_matrix(nrow, nrow);
      sq_rsp(nrow, nrow, Tsq, E, 0, C, 1e-12);

      /* Find maximum eigenvalue of T */
      for(i=0; i < nrow; i++) if(E[i] > max) max = E[i];
	     
      free_block(C);
      free(E);
      free_block(Tsq);
      free_block(Tsub);
    }

  return max;
}

static double
d1diag_t1_rohf()
{
  int h, nirreps, i, j;
  double **Tave, tmp, max;
  double max_ph=0.0, max_xp=0.0, max_hx=0.0;
  dpdfile2 T1_a, T1_b;

  nirreps = moinfo.nirreps;

  dpd_file2_init(&T1_a, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&T1_a);
  dpd_file2_mat_rd(&T1_a);
      
  dpd_file2_init(&T1_b, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1_b);
  dpd_file2_mat_rd(&T1_b);

  for(h=0; h < nirreps; h++) {
      int nrow = T1_a.params->rowtot[h];
      int ncol = T1_a.params->coltot[h];
      int nopen = moinfo.openpi[h];
      if(nrow && ncol) {
         Tave = block_matrix(nrow, ncol);

         for (i=0; i<nrow; i++) {
             for (j=0; j<ncol; j++) {
                 Tave[i][j] = (T1_a.matrix[h][i][j] + T1_b.matrix[h][i][j])/2.;
               }
           }

         tmp = d1diag_subblock(Tave, 0, nrow-nopen, 0, ncol-nopen);
         if (tmp > max_ph) max_ph = tmp;

         tmp = d1diag_subblock(Tave, 0, nrow-nopen, ncol-nopen, ncol);
         if (tmp > max_hx) max_hx = tmp;

         tmp = d1diag_subblock(Tave, nrow-nopen, nrow, 0, ncol-nopen);
         if (tmp > max_xp) max_xp = tmp;

         free_block(Tave);
       }
    }

  dpd_file2_mat_close(&T1_a);
  dpd_file2_close(&T1_a);

  dpd_file2_mat_close(&T1_b);
  dpd_file2_close(&T1_b);

  max_ph = sqrt(max_ph);
  max_hx = sqrt(max_hx);
  max_xp = sqrt(max_xp);

  /*
  fprintf(outfile, "D1:  hp=%8.6f hx=%8.6f xp=%8.6f\n", max_ph, max_hx, max_xp);
  */

  max = max_ph;
  if (max_hx > max) max = max_hx;
  if (max_xp > max) max = max_xp;

  return max;
}

double d1diag(void)
{
  double norm = 0.0;

  if(params.ref == 0) { /** RHF **/
    norm = d1diag_t1_rhf();
  }
  else if (params.ref == 1) { /** ROHF **/
    norm = d1diag_t1_rohf();
  }
  return norm;
}
}} // namespace psi::ccenergy
