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
      dpd_file2_init(&XIA, EOM_XI, L_irr, 0, 1, "XIA");
      dpd_file2_copy(&XIA, CC_LAMBDA, "LIA");
      dpd_file2_close(&XIA);
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_file2_init(&dIA, CC_DENOM, L_irr, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &LIA);
      dpd_file2_close(&dIA);
      dpd_file2_close(&LIA);

      dpd_buf4_init(&XIjAb, EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
      dpd_buf4_copy(&XIjAb, CC_LAMBDA, "LIjAb");
      dpd_buf4_close(&XIjAb);
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_init(&dIjAb, CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
      dpd_buf4_dirprd(&dIjAb, &LIjAb);
      dpd_buf4_close(&dIjAb);
      dpd_buf4_close(&LIjAb);
    }
    else if (params.ref == 1) { /* ROHF */
      dpd_file2_init(&XIA, EOM_XI, L_irr, 0, 1, "XIA");
      dpd_file2_copy(&XIA, CC_LAMBDA, "LIA");
      dpd_file2_close(&XIA);
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_file2_init(&dIA, CC_DENOM, L_irr, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &LIA);
      dpd_file2_close(&dIA);
      dpd_file2_close(&LIA);

      dpd_file2_init(&Xia, EOM_XI, L_irr, 0, 1, "Xia");
      dpd_file2_copy(&Xia, CC_LAMBDA, "Lia");
      dpd_file2_close(&Xia);
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
      dpd_file2_init(&dia, CC_DENOM, L_irr, 0, 1, "dia");
      dpd_file2_dirprd(&dia, &Lia);
      dpd_file2_close(&dia);
      dpd_file2_close(&Lia);

      dpd_buf4_init(&XIJAB, EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
      dpd_buf4_copy(&XIJAB, CC_LAMBDA, "LIJAB");
      dpd_buf4_close(&XIJAB);
      dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_init(&dIJAB, CC_DENOM, L_irr, 2, 7, 2, 7, 0, "dIJAB");
      dpd_buf4_dirprd(&dIJAB, &LIJAB);
      dpd_buf4_close(&dIJAB);
      dpd_buf4_close(&LIJAB);

      dpd_buf4_init(&Xijab, EOM_XI, L_irr, 2, 7, 2, 7, 0, "Xijab");
      dpd_buf4_copy(&Xijab, CC_LAMBDA, "Lijab");
      dpd_buf4_close(&Xijab);
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_init(&dijab, CC_DENOM, L_irr, 2, 7, 2, 7, 0, "dijab");
      dpd_buf4_dirprd(&dijab, &Lijab);
      dpd_buf4_close(&dijab);
      dpd_buf4_close(&Lijab);

      dpd_buf4_init(&XIjAb, EOM_XI, L_irr, 0, 5, 0, 5, 0, "XIjAb");
      dpd_buf4_copy(&XIjAb, CC_LAMBDA, "LIjAb");
      dpd_buf4_close(&XIjAb);
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_init(&dIjAb, CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
      dpd_buf4_dirprd(&dIjAb, &LIjAb);
      dpd_buf4_close(&dIjAb);
      dpd_buf4_close(&LIjAb);
    }
    else if(params.ref == 2) { /** UHF **/
      dpd_file2_init(&XIA, EOM_XI, L_irr, 0, 1, "XIA");
      dpd_file2_copy(&XIA, CC_LAMBDA, "LIA");
      dpd_file2_close(&XIA);
      dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
      dpd_file2_init(&dIA, CC_DENOM, L_irr, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &LIA);
      dpd_file2_close(&dIA);
      dpd_file2_close(&LIA);

      dpd_file2_init(&Xia, EOM_XI, L_irr, 2, 3, "Xia");
      dpd_file2_copy(&Xia, CC_LAMBDA, "Lia");
      dpd_file2_close(&Xia);
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");
      dpd_file2_init(&dia, CC_DENOM, L_irr, 2, 3, "dia");
      dpd_file2_dirprd(&dia, &Lia);
      dpd_file2_close(&dia);
      dpd_file2_close(&Lia);

      dpd_buf4_init(&XIJAB, EOM_XI, L_irr, 2, 7, 2, 7, 0, "XIJAB");
      dpd_buf4_copy(&XIJAB, CC_LAMBDA, "LIJAB");
      dpd_buf4_close(&XIJAB);
      dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_init(&dIJAB, CC_DENOM, L_irr, 2, 7, 2, 7, 0, "dIJAB");
      dpd_buf4_dirprd(&dIJAB, &LIJAB);
      dpd_buf4_close(&dIJAB);
      dpd_buf4_close(&LIJAB);

      dpd_buf4_init(&Xijab, EOM_XI, L_irr, 12, 17, 12, 17, 0, "Xijab");
      dpd_buf4_copy(&Xijab, CC_LAMBDA, "Lijab");
      dpd_buf4_close(&Xijab);
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      dpd_buf4_init(&dijab, CC_DENOM, L_irr, 12, 17, 12, 17, 0, "dijab");
      dpd_buf4_dirprd(&dijab, &Lijab);
      dpd_buf4_close(&dijab);
      dpd_buf4_close(&Lijab);

      dpd_buf4_init(&XIjAb, EOM_XI, L_irr, 22, 28, 22, 28, 0, "XIjAb");
      dpd_buf4_copy(&XIjAb, CC_LAMBDA, "LIjAb");
      dpd_buf4_close(&XIjAb);
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      dpd_buf4_init(&dIjAb, CC_DENOM, L_irr, 22, 28, 22, 28, 0, "dIjAb");
      dpd_buf4_dirprd(&dIjAb, &LIjAb);
      dpd_buf4_close(&dIjAb);
      dpd_buf4_close(&LIjAb);
    }
    return;
  }

  /* ground state guess L <= T */
  /* excited state guess L <= R0 * T + R */
  if (L_params.ground || L_params.irrep == 0) {
    if(params.ref == 0) { /** RHF **/
      if(!params.restart || !psio_tocscan(CC_LAMBDA, "LIA")) {
        dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
        dpd_file2_copy(&T1, CC_LAMBDA, "LIA");
        dpd_file2_copy(&T1, CC_LAMBDA, "Lia");
        dpd_file2_close(&T1);
      }
      else fprintf(outfile, "\tUsing old L1 amplitudes.\n");

      if(!params.restart || !psio_tocscan(CC_LAMBDA, "LIjAb")) {
        dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        dpd_buf4_copy(&T2, CC_LAMBDA, "LIjAb");
        dpd_buf4_close(&T2);
      }
      else fprintf(outfile, "\tUsing old L2 amplitudes.\n");

      dpd_buf4_init(&T2, CC_LAMBDA, 0, 2, 7, 0, 5, 1, "LIjAb");
      dpd_buf4_copy(&T2, CC_LAMBDA, "LIJAB");
      dpd_buf4_copy(&T2, CC_LAMBDA, "Lijab");
      dpd_buf4_close(&T2);
    }
    else if(params.ref == 1) { /** ROHF **/
      if(!params.restart || !psio_tocscan(CC_LAMBDA, "LIA") ||
         !psio_tocscan(CC_LAMBDA, "Lia")) {
        dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
        dpd_file2_copy(&T1, CC_LAMBDA, "LIA");
        dpd_file2_close(&T1);
  
        dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
        dpd_file2_copy(&T1, CC_LAMBDA, "Lia");
        dpd_file2_close(&T1);
      }
      else fprintf(outfile, "\tUsing old L1 amplitudes.\n");
  
      if(!params.restart || !psio_tocscan(CC_LAMBDA, "LIjAb") ||
         !psio_tocscan(CC_LAMBDA, "LIJAB") || 
         !psio_tocscan(CC_LAMBDA, "Lijab")) {
        dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
        dpd_buf4_copy(&T2, CC_LAMBDA, "LIJAB");
        dpd_buf4_close(&T2);
  
        dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
        dpd_buf4_copy(&T2, CC_LAMBDA, "Lijab");
        dpd_buf4_close(&T2);

        dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        dpd_buf4_copy(&T2, CC_LAMBDA, "LIjAb");
        dpd_buf4_close(&T2);
      }
      else fprintf(outfile, "\tUsing old L2 amplitudes.\n");
    }
    else if(params.ref == 2) { /** UHF **/
      if(!params.restart || !psio_tocscan(CC_LAMBDA, "LIA") ||
         !psio_tocscan(CC_LAMBDA, "Lia")) {
        dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
        dpd_file2_copy(&T1, CC_LAMBDA, "LIA");
        dpd_file2_close(&T1);
  
        dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
        dpd_file2_copy(&T1, CC_LAMBDA, "Lia");
        dpd_file2_close(&T1);
      }
      else fprintf(outfile, "\tUsing old L1 amplitudes.\n");
  
      if(!params.restart || !psio_tocscan(CC_LAMBDA, "LIjAb") ||
         !psio_tocscan(CC_LAMBDA, "LIJAB") || 
         !psio_tocscan(CC_LAMBDA, "Lijab")) {
        dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
        dpd_buf4_copy(&T2, CC_LAMBDA, "LIJAB");
        dpd_buf4_close(&T2);
  
        dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
        dpd_buf4_copy(&T2, CC_LAMBDA, "Lijab");
        dpd_buf4_close(&T2);
  
        dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        dpd_buf4_copy(&T2, CC_LAMBDA, "LIjAb");
        dpd_buf4_close(&T2);
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
    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    if (params.ref <= 1) {
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    }
    else {
      dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");
      dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    }

    dpd_file2_scm(&LIA, L_params.R0);
    dpd_file2_scm(&Lia, L_params.R0);
    dpd_buf4_scm(&LIJAB, L_params.R0);
    dpd_buf4_scm(&Lijab, L_params.R0);
    dpd_buf4_scm(&LIjAb, L_params.R0);
  
      /* add R1 and R2 */
    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    dpd_file2_axpy(&R1, &LIA, 1.0, 0);
    dpd_file2_close(&R1);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    dpd_buf4_axpy(&R2, &LIJAB, 1.0);
    dpd_buf4_close(&R2);

    if (params.ref <= 1) {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      dpd_file2_axpy(&R1, &Lia, 1.0, 0);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      dpd_buf4_axpy(&R2, &Lijab, 1.0);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      dpd_buf4_axpy(&R2, &LIjAb, 1.0);
      dpd_buf4_close(&R2);
    }
    else {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      dpd_file2_axpy(&R1, &Lia, 1.0, 0);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      dpd_buf4_axpy(&R2, &Lijab, 1.0);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      dpd_buf4_axpy(&R2, &LIjAb, 1.0);
      dpd_buf4_close(&R2);
    }
  
    /* dot L and R together */
    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    norm = dpd_file2_dot(&LIA, &R1);
    dpd_file2_close(&R1);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += dpd_buf4_dot(&LIJAB, &R2);
    dpd_buf4_close(&R2);
    if (params.ref <= 1) {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
    else {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
  
    fprintf(outfile,"\tInitial overlap of initial guess <L|R> = %15.10lf\n", norm);
  
    dpd_file2_scm(&LIA, 1.0/norm);
    dpd_file2_scm(&Lia, 1.0/norm);
    dpd_buf4_scm(&LIJAB, 1.0/norm);
    dpd_buf4_scm(&Lijab, 1.0/norm);
    dpd_buf4_scm(&LIjAb, 1.0/norm);
  
    dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    norm = dpd_file2_dot(&LIA, &R1);
    dpd_file2_close(&R1);
    dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    norm += dpd_buf4_dot(&LIJAB, &R2);
    dpd_buf4_close(&R2);

    if (params.ref <= 1) {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 0, 1, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
    else {
      dpd_file2_init(&R1, CC_RAMPS, L_irr, 2, 3, R1B_lbl);
      norm += dpd_file2_dot(&Lia, &R1);
      dpd_file2_close(&R1);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
      norm += dpd_buf4_dot(&Lijab, &R2);
      dpd_buf4_close(&R2);
      dpd_buf4_init(&R2, CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
      norm += dpd_buf4_dot(&LIjAb, &R2);
      dpd_buf4_close(&R2);
    }
    fprintf(outfile,"\tChecking overlap of initial guess <L|R> = %15.10lf\n", norm);
  
    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&LIjAb);
  }

#ifdef EOM_DEBUG
  fprintf(outfile,"initial guess\n");
  dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
  dpd_file2_print(&LIA,outfile);
  dpd_file2_close(&LIA);
#endif
}

}} // namespace psi::cclambda
