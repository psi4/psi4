/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* BUILD_X(): Construct the orbital rotation gradient, XAI, for CC
** gradient calculations:
**
**  Xai = I'ia - I'ai
** */

void build_X(void)
{
  dpdfile2 X, I, I1, I2, I3, F, F1, X2;
  dpdbuf4 E, Fi;
  

  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    global_dpd_->file2_copy(&I, PSIF_CC_OEI, "XAI");
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 1, 0, "XAI");
    global_dpd_->file2_scm(&X, -1.0);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    global_dpd_->file2_axpy(&I, &X, 1.0, 1);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&X);

    /* Add the orbital response contributions of dependent pairs 
     (i,j) and (a,b) to the X_AI due to the use of canonical 
     perturbed orbitals */
    if(params.wfn == "CCSD_T" && params.dertype ==1){
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);

    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&F1);
    global_dpd_->file2_mat_rd(&F1);

    /* delta_I/delta_f_IJ = (I'ij - I'ji)/(fii - fjj)  if (I'ij - I'ji) > 1e-8 */
    /* delta_I/delta_f_AB = (I'ab - I'ba)/(faa - fbb)  if (I'ab - I'ba) > 1e-8 */
    global_dpd_->file2_init(&I, PSIF_CC_TMP, 0, 0, 0, "delta_I/delta_f_IJ");
    global_dpd_->file2_init(&I1, PSIF_CC_TMP, 0, 1, 1, "delta_I/delta_f_AB");
    global_dpd_->file2_init(&I2, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
    global_dpd_->file2_init(&I3, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    global_dpd_->file2_mat_init(&I);
    global_dpd_->file2_mat_init(&I1);
    global_dpd_->file2_mat_init(&I2);
    global_dpd_->file2_mat_init(&I3);
    global_dpd_->file2_mat_rd(&I2);
    global_dpd_->file2_mat_rd(&I3);

    for(int h=0; h < moinfo.nirreps; h++){
      for(int i=0; i < moinfo.occpi[h]; i++)
          for(int j=0; j < moinfo.occpi[h]; j++){
             double diff = fabs(I2.matrix[h][i][j] - I2.matrix[h][j][i]); 
             if (diff > 1e-8) 
		  I.matrix[h][i][j] = (I2.matrix[h][i][j] - I2.matrix[h][j][i])/(F.matrix[h][i][i] - F.matrix[h][j][j]) ;
             else 
		  I.matrix[h][i][j] = 0.0 ;
      }
      for(int a=0; a < moinfo.virtpi[h]; a++)
          for(int b=0; b < moinfo.virtpi[h]; b++) {
             double diff = fabs(I3.matrix[h][a][b] - I3.matrix[h][b][a]); 
             if (diff > 1e-8) 
		  I1.matrix[h][a][b] = (I3.matrix[h][a][b] - I3.matrix[h][b][a]) /(F1.matrix[h][a][a] - F1.matrix[h][b][b]) ;
             else 
		  I1.matrix[h][a][b] = 0.0 ;
    }
  }
    global_dpd_->file2_mat_wrt(&I);
    global_dpd_->file2_mat_wrt(&I1);
    global_dpd_->file2_mat_close(&I);
    global_dpd_->file2_mat_close(&I1);
    global_dpd_->file2_mat_close(&I2);
    global_dpd_->file2_mat_close(&I3);
    global_dpd_->file2_mat_close(&F);
    global_dpd_->file2_mat_close(&F1);

    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&I1);
    global_dpd_->file2_close(&I2);
    global_dpd_->file2_close(&I3);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&F1);

    global_dpd_->file2_init(&I, PSIF_CC_TMP, 0, 0, 0, "delta_I/delta_f_IJ");
    global_dpd_->file2_init(&I1, PSIF_CC_TMP, 0, 1, 1, "delta_I/delta_f_AB");
    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 1, 0, "XAI");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
    global_dpd_->buf4_scmcopy(&E, PSIF_CC_EINTS, "4 <ka|ji> - <ka|ij> - <ki|aj>", 4);
    global_dpd_->buf4_sort_axpy(&E, PSIF_CC_EINTS, pqsr, 10, 0, "4 <ka|ji> - <ka|ij> - <ki|aj>", -1);
    global_dpd_->buf4_sort_axpy(&E, PSIF_CC_EINTS, rqsp, 10, 0, "4 <ka|ji> - <ka|ij> - <ki|aj>", -1);
    global_dpd_->buf4_close(&E);

    /* Xai += 1/2 sum_kj (4<ka|ji> - <ka|ij> - <ki|aj>) * (I'kj - I'jk)/(fkk - fjj) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "4 <ka|ji> - <ka|ij> - <ki|aj>");
    global_dpd_->dot13(&I, &E, &X, 0, 0, 0.5, 1.0);
    global_dpd_->buf4_close(&E);
    global_dpd_->file2_close(&X);
    global_dpd_->file2_close(&I);


    global_dpd_->buf4_init(&Fi, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_scmcopy(&Fi, PSIF_CC_FINTS, "4 <bi|ca> - <bi|ac> - <ci|ab>", 4);
    global_dpd_->buf4_sort_axpy(&Fi, PSIF_CC_FINTS, pqsr, 11, 5, "4 <bi|ca> - <bi|ac> - <ci|ab>", -1);
    global_dpd_->buf4_sort_axpy(&Fi, PSIF_CC_FINTS, rqsp, 11, 5, "4 <bi|ca> - <bi|ac> - <ci|ab>", -1);
    global_dpd_->buf4_close(&Fi);

    /* Xia(tmp) = 1/2 sum_bc (4<bi|ca> - <bi|ac> - <ci|ab>) * (I'bc - I'cb)/(fbb - fcc) */
    global_dpd_->buf4_init(&Fi, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "4 <bi|ca> - <bi|ac> - <ci|ab>");
    global_dpd_->file2_init(&X2, PSIF_CC_OEI, 0, 0, 1, "XIA_tmp");
    global_dpd_->dot13(&I1, &Fi, &X2, 0, 0, 0.5, 1.0);
    global_dpd_->buf4_close(&Fi);
    global_dpd_->file2_close(&I1);
    global_dpd_->file2_close(&X2);

    /* Xai +=  transpose(XIA_tmp) */
    global_dpd_->file2_init(&X2, PSIF_CC_OEI, 0, 0, 1, "XIA_tmp");
    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 1, 0, "XAI");
    global_dpd_->file2_axpy(&X2, &X, 1, 1);
    global_dpd_->file2_close(&X2);
    global_dpd_->file2_close(&X);

    }
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    global_dpd_->file2_copy(&I, PSIF_CC_OEI, "XAI");
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 1, 0, "XAI");
    global_dpd_->file2_scm(&X, -1.0);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    global_dpd_->file2_axpy(&I, &X, 1.0, 1);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&X);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");
    global_dpd_->file2_copy(&I, PSIF_CC_OEI, "Xai");
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 1, 0, "Xai");
    global_dpd_->file2_scm(&X, -1.0);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'ia");
    global_dpd_->file2_axpy(&I, &X, 1.0, 1);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&X);

    /* Build spatial orbital version of X for Zvector:
       X(A,I) = 1/2[X(A,I)+X(a,i)]  */
    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 1, 0, "XAI");
    global_dpd_->file2_copy(&X, PSIF_CC_MISC, "X(A,I)");
    global_dpd_->file2_close(&X);
    global_dpd_->file2_init(&X, PSIF_CC_MISC, 0, 1, 0, "X(A,I)");
    global_dpd_->file2_init(&X2, PSIF_CC_OEI, 0, 1, 0, "Xai");
    global_dpd_->file2_axpy(&X2, &X, 1.0, 0);
    global_dpd_->file2_close(&X2);
    global_dpd_->file2_scm(&X, 0.5);
    global_dpd_->file2_close(&X);

  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    global_dpd_->file2_copy(&I, PSIF_CC_OEI, "XAI");
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 1, 0, "XAI");
    global_dpd_->file2_scm(&X, -1.0);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    global_dpd_->file2_axpy(&I, &X, 1.0, 1);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&X);

    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 3, 2, "I'ai");
    global_dpd_->file2_copy(&I, PSIF_CC_OEI, "Xai");
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&X, PSIF_CC_OEI, 0, 3, 2, "Xai");
    global_dpd_->file2_scm(&X, -1.0);
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 3, "I'ia");
    global_dpd_->file2_axpy(&I, &X, 1.0, 1);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&X);

  }

}

}} // namespace psi::ccdensity
