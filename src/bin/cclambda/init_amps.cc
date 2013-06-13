/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void init_amps(struct L_Params L_params)
{
  double norm;
  dpdfile2 T1, R1, LIA, Lia, dIA, dia, XIA, Xia;
  dpdbuf4 T2, R2, LIJAB, Lijab, LIjAb, dIJAB, dijab, dIjAb, XIJAB, Xijab, XIjAb;
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];
  int L_irr;
  L_irr = L_params.irrep;

  /* if solving zeta equations, initial guess is Xi * denom */
  if (params.zeta) {
    if (params.ref == 0) { /* RHF */
      dpd_->file2_init(&XIA, PSIF_EOM_XI, L_irr, 0, 1, "XIA");
      dpd_->file2_copy(&XIA, PSIF_CC_LAMBDA, "LIA");
      dpd_->file2_close(&XIA);
      dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
      dpd_->file2_dirprd(&dIA, &LIA);
      dpd_->file2_close(&dIA);
      dpd_->file2_close(&LIA);

      dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
      dpd_->buf4_copy(&XIjAb, PSIF_CC_LAMBDA, "LIjAb");
      dpd_->buf4_close(&XIjAb);
      dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
      dpd_->buf4_dirprd(&dIjAb, &LIjAb);
      dpd_->buf4_close(&dIjAb);
      dpd_->buf4_close(&LIjAb);
    }
    else if (params.ref == 1) { /* ROHF */
      dpd_->file2_init(&XIA, PSIF_EOM_XI, L_irr, 0, 1, "XIA");
      dpd_->file2_copy(&XIA, PSIF_CC_LAMBDA, "LIA");
      dpd_->file2_close(&XIA);
      dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
      dpd_->file2_dirprd(&dIA, &LIA);
      dpd_->file2_close(&dIA);
      dpd_->file2_close(&LIA);

      dpd_->file2_init(&Xia, PSIF_EOM_XI, L_irr, 0, 1, "Xia");
      dpd_->file2_copy(&Xia, PSIF_CC_LAMBDA, "Lia");
      dpd_->file2_close(&Xia);
      dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
      dpd_->file2_init(&dia, PSIF_CC_DENOM, L_irr, 0, 1, "dia");
      dpd_->file2_dirprd(&dia, &Lia);
      dpd_->file2_close(&dia);
      dpd_->file2_close(&Lia);

      dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
      dpd_->buf4_copy(&XIJAB, PSIF_CC_LAMBDA, "LIJAB");
      dpd_->buf4_close(&XIJAB);
      dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, L_irr, 2, 7, 2, 7, 0, "dIJAB");
      dpd_->buf4_dirprd(&dIJAB, &LIJAB);
      dpd_->buf4_close(&dIJAB);
      dpd_->buf4_close(&LIJAB);

      dpd_->buf4_init(&Xijab, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "Xijab");
      dpd_->buf4_copy(&Xijab, PSIF_CC_LAMBDA, "Lijab");
      dpd_->buf4_close(&Xijab);
      dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_->buf4_init(&dijab, PSIF_CC_DENOM, L_irr, 2, 7, 2, 7, 0, "dijab");
      dpd_->buf4_dirprd(&dijab, &Lijab);
      dpd_->buf4_close(&dijab);
      dpd_->buf4_close(&Lijab);

      dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
      dpd_->buf4_copy(&XIjAb, PSIF_CC_LAMBDA, "LIjAb");
      dpd_->buf4_close(&XIjAb);
      dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
      dpd_->buf4_dirprd(&dIjAb, &LIjAb);
      dpd_->buf4_close(&dIjAb);
      dpd_->buf4_close(&LIjAb);
    }
    else if(params.ref == 2) { /** UHF **/
      dpd_->file2_init(&XIA, PSIF_EOM_XI, L_irr, 0, 1, "XIA");
      dpd_->file2_copy(&XIA, PSIF_CC_LAMBDA, "LIA");
      dpd_->file2_close(&XIA);
      dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_->file2_init(&dIA, PSIF_CC_DENOM, L_irr, 0, 1, "dIA");
      dpd_->file2_dirprd(&dIA, &LIA);
      dpd_->file2_close(&dIA);
      dpd_->file2_close(&LIA);

      dpd_->file2_init(&Xia, PSIF_EOM_XI, L_irr, 2, 3, "Xia");
      dpd_->file2_copy(&Xia, PSIF_CC_LAMBDA, "Lia");
      dpd_->file2_close(&Xia);
      dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
      dpd_->file2_init(&dia, PSIF_CC_DENOM, L_irr, 2, 3, "dia");
      dpd_->file2_dirprd(&dia, &Lia);
      dpd_->file2_close(&dia);
      dpd_->file2_close(&Lia);

      dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
      dpd_->buf4_copy(&XIJAB, PSIF_CC_LAMBDA, "LIJAB");
      dpd_->buf4_close(&XIJAB);
      dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, L_irr, 2, 7, 2, 7, 0, "dIJAB");
      dpd_->buf4_dirprd(&dIJAB, &LIJAB);
      dpd_->buf4_close(&dIJAB);
      dpd_->buf4_close(&LIJAB);

      dpd_->buf4_init(&Xijab, PSIF_EOM_XI, L_irr, 12, 17, 12, 17, 0, "Xijab");
      dpd_->buf4_copy(&Xijab, PSIF_CC_LAMBDA, "Lijab");
      dpd_->buf4_close(&Xijab);
      dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      dpd_->buf4_init(&dijab, PSIF_CC_DENOM, L_irr, 12, 17, 12, 17, 0, "dijab");
      dpd_->buf4_dirprd(&dijab, &Lijab);
      dpd_->buf4_close(&dijab);
      dpd_->buf4_close(&Lijab);

      dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, L_irr, 22, 28, 22, 28, 0, "XIjAb");
      dpd_->buf4_copy(&XIjAb, PSIF_CC_LAMBDA, "LIjAb");
      dpd_->buf4_close(&XIjAb);
      dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, L_irr, 22, 28, 22, 28, 0, "dIjAb");
      dpd_->buf4_dirprd(&dIjAb, &LIjAb);
      dpd_->buf4_close(&dIjAb);
      dpd_->buf4_close(&LIjAb);
    }
    return;
  }

  /* ground state guess L <= T */
  /* excited state guess L <= R0 * T + R */
  if (L_params.ground || L_params.irrep == 0) {
    if(params.ref == 0) { /** RHF **/
      if(!params.restart || !psio_tocscan(PSIF_CC_LAMBDA, "LIA")) {
        dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        dpd_->file2_copy(&T1, PSIF_CC_LAMBDA, "LIA");
        dpd_->file2_copy(&T1, PSIF_CC_LAMBDA, "Lia");
        dpd_->file2_close(&T1);
      }
      else fprintf(outfile, "\tUsing old L1 amplitudes.\n");

      if(!params.restart || !psio_tocscan(PSIF_CC_LAMBDA, "LIjAb")) {
        dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "LIjAb");
        dpd_->buf4_close(&T2);
      }
      else fprintf(outfile, "\tUsing old L2 amplitudes.\n");

      dpd_->buf4_init(&T2, PSIF_CC_LAMBDA, 0, 2, 7, 0, 5, 1, "LIjAb");
      dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "LIJAB");
      dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "Lijab");
      dpd_->buf4_close(&T2);
    }
    else if(params.ref == 1) { /** ROHF **/
      if(!params.restart || !psio_tocscan(PSIF_CC_LAMBDA, "LIA") ||
         !psio_tocscan(PSIF_CC_LAMBDA, "Lia")) {
        dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        dpd_->file2_copy(&T1, PSIF_CC_LAMBDA, "LIA");
        dpd_->file2_close(&T1);
  
        dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
        dpd_->file2_copy(&T1, PSIF_CC_LAMBDA, "Lia");
        dpd_->file2_close(&T1);
      }
      else fprintf(outfile, "\tUsing old L1 amplitudes.\n");
  
      if(!params.restart || !psio_tocscan(PSIF_CC_LAMBDA, "LIjAb") ||
         !psio_tocscan(PSIF_CC_LAMBDA, "LIJAB") || 
         !psio_tocscan(PSIF_CC_LAMBDA, "Lijab")) {
        dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
        dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "LIJAB");
        dpd_->buf4_close(&T2);
  
        dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
        dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "Lijab");
        dpd_->buf4_close(&T2);

        dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "LIjAb");
        dpd_->buf4_close(&T2);
      }
      else fprintf(outfile, "\tUsing old L2 amplitudes.\n");
    }
    else if(params.ref == 2) { /** UHF **/
      if(!params.restart || !psio_tocscan(PSIF_CC_LAMBDA, "LIA") ||
         !psio_tocscan(PSIF_CC_LAMBDA, "Lia")) {
        dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        dpd_->file2_copy(&T1, PSIF_CC_LAMBDA, "LIA");
        dpd_->file2_close(&T1);
  
        dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
        dpd_->file2_copy(&T1, PSIF_CC_LAMBDA, "Lia");
        dpd_->file2_close(&T1);
      }
      else fprintf(outfile, "\tUsing old L1 amplitudes.\n");
  
      if(!params.restart || !psio_tocscan(PSIF_CC_LAMBDA, "LIjAb") ||
         !psio_tocscan(PSIF_CC_LAMBDA, "LIJAB") || 
         !psio_tocscan(PSIF_CC_LAMBDA, "Lijab")) {
        dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
        dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "LIJAB");
        dpd_->buf4_close(&T2);
  
        dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
        dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "Lijab");
        dpd_->buf4_close(&T2);
  
        dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        dpd_->buf4_copy(&T2, PSIF_CC_LAMBDA, "LIjAb");
        dpd_->buf4_close(&T2);
      }
      else fprintf(outfile, "\tUsing old L2 amplitudes.\n");
    }
  }

  if (!L_params.ground) {
    sprintf(R1A_lbl, "RIA %d %d", L_params.irrep, L_params.root);
    sprintf(R1B_lbl, "Ria %d %d", L_params.irrep, L_params.root);
    sprintf(R2AA_lbl, "RIJAB %d %d", L_params.irrep, L_params.root);
    sprintf(R2BB_lbl, "Rijab %d %d", L_params.irrep, L_params.root);
    sprintf(R2AB_lbl, "RIjAb %d %d", L_params.irrep, L_params.root);

    /* multiply by R0 and create nonsymmetric L files */
    dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    if (params.ref <= 1) {
      dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
      dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    }
    else {
      dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
      dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    }

    dpd_->file2_scm(&LIA, L_params.R0);
    dpd_->file2_scm(&Lia, L_params.R0);
    dpd_->buf4_scm(&LIJAB, L_params.R0);
    dpd_->buf4_scm(&Lijab, L_params.R0);
    dpd_->buf4_scm(&LIjAb, L_params.R0);
  
      /* add R1 and R2 */
    dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    dpd_->file2_axpy(&R1, &LIA, 1.0, 0);
    dpd_->file2_close(&R1);
    dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    dpd_->buf4_axpy(&R2, &LIJAB, 1.0);
    dpd_->buf4_close(&R2);

    if (params.ref <= 1) {
      dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      dpd_->file2_axpy(&R1, &Lia, 1.0, 0);
      dpd_->file2_close(&R1);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      dpd_->buf4_axpy(&R2, &Lijab, 1.0);
      dpd_->buf4_close(&R2);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      dpd_->buf4_axpy(&R2, &LIjAb, 1.0);
      dpd_->buf4_close(&R2);
    }
    else {
      dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      dpd_->file2_axpy(&R1, &Lia, 1.0, 0);
      dpd_->file2_close(&R1);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      dpd_->buf4_axpy(&R2, &Lijab, 1.0);
      dpd_->buf4_close(&R2);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      dpd_->buf4_axpy(&R2, &LIjAb, 1.0);
      dpd_->buf4_close(&R2);
    }
  
    /* dot L and R together */
    dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    norm = dpd_->file2_dot(&LIA, &R1);
    dpd_->file2_close(&R1);
    dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += dpd_->buf4_dot(&LIJAB, &R2);
    dpd_->buf4_close(&R2);
    if (params.ref <= 1) {
      dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      norm += dpd_->file2_dot(&Lia, &R1);
      dpd_->file2_close(&R1);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      norm += dpd_->buf4_dot(&Lijab, &R2);
      dpd_->buf4_close(&R2);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      norm += dpd_->buf4_dot(&LIjAb, &R2);
      dpd_->buf4_close(&R2);
    }
    else {
      dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      norm += dpd_->file2_dot(&Lia, &R1);
      dpd_->file2_close(&R1);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      norm += dpd_->buf4_dot(&Lijab, &R2);
      dpd_->buf4_close(&R2);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      norm += dpd_->buf4_dot(&LIjAb, &R2);
      dpd_->buf4_close(&R2);
    }
  
    fprintf(outfile,"\tInitial overlap of initial guess <L|R> = %15.10lf\n", norm);
  
    dpd_->file2_scm(&LIA, 1.0/norm);
    dpd_->file2_scm(&Lia, 1.0/norm);
    dpd_->buf4_scm(&LIJAB, 1.0/norm);
    dpd_->buf4_scm(&Lijab, 1.0/norm);
    dpd_->buf4_scm(&LIjAb, 1.0/norm);
  
    dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    norm = dpd_->file2_dot(&LIA, &R1);
    dpd_->file2_close(&R1);
    dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += dpd_->buf4_dot(&LIJAB, &R2);
    dpd_->buf4_close(&R2);

    if (params.ref <= 1) {
      dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      norm += dpd_->file2_dot(&Lia, &R1);
      dpd_->file2_close(&R1);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      norm += dpd_->buf4_dot(&Lijab, &R2);
      dpd_->buf4_close(&R2);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      norm += dpd_->buf4_dot(&LIjAb, &R2);
      dpd_->buf4_close(&R2);
    }
    else {
      dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      norm += dpd_->file2_dot(&Lia, &R1);
      dpd_->file2_close(&R1);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      norm += dpd_->buf4_dot(&Lijab, &R2);
      dpd_->buf4_close(&R2);
      dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      norm += dpd_->buf4_dot(&LIjAb, &R2);
      dpd_->buf4_close(&R2);
    }
    fprintf(outfile,"\tChecking overlap of initial guess <L|R> = %15.10lf\n", norm);
  
    dpd_->file2_close(&LIA);
    dpd_->file2_close(&Lia);
    dpd_->buf4_close(&LIJAB);
    dpd_->buf4_close(&Lijab);
    dpd_->buf4_close(&LIjAb);
  }

#ifdef EOM_DEBUG
  fprintf(outfile,"initial guess\n");
  dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
  dpd_file2_print(&LIA,outfile);
  dpd_file2_close(&LIA);
#endif
}

}} // namespace psi::cclambda
