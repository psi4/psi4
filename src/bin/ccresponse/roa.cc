/*! \file
    \ingroup ccresponse
    \brief Compute the three tensors needed for Raman Optical Activity.

    ROA requires the following polarizability tensors:
      (1) electric-dipole/electric-dipole; 
      (2) electric-dipole/electric-quadrupole; and 
      (3) electric-dipole/magnetic-dipole.

  -TDC, August 2009
*/
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void pertbar(const char *pert, int irrep, int anti);
void compute_X(const char *pert, int irrep, double omega);
void linresp(double *tensor, double A, double B,
	     const char *pert_x, int x_irrep, double omega_x,
	     const char *pert_y, int y_irrep, double omega_y);

void roa(void)
{
  double ***tensor_rl, ***tensor_pl, **tensor0, ***tensor_rr; 
  double ****tensor_rQ, ***tensor_rQ0, ***tensor_rQ1;
  double **tensor_rl0, **tensor_rl1, **tensor_pl0, **tensor_pl1;
  char **cartcomp, pert[32], pert_x[32], pert_y[32];
  int alpha, beta, gamma, i, j, k, l, irrep;
  double omega_nm, omega_ev, omega_cm;
  char lbl1[32], lbl2[32], lbl3[32], lbl4[32];
  int compute_rl=0, compute_pl=0;
  psio_address next;
  double value;

  /* Booleans for convenience */
  if(params.gauge == "LENGTH" || params.gauge == "BOTH") 
    compute_rl=1;
  if(params.gauge == "VELOCITY" || params.gauge == "BOTH") 
    compute_pl=1;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor_rQ = (double ****) malloc(params.nomega * sizeof(double ***));
  for(i=0; i < params.nomega; i++) {
    tensor_rQ[i] = (double ***) malloc(3 * sizeof(double **));
    for(j=0; j < 3; j++) tensor_rQ[i][j] = block_matrix(3,3);
  }
  tensor_rl = (double ***) malloc(params.nomega * sizeof(double **));
  tensor_pl = (double ***) malloc(params.nomega * sizeof(double **));
  tensor_rr = (double ***) malloc(params.nomega * sizeof(double **));
  for(i=0; i < params.nomega; i++) {
    tensor_rl[i] = block_matrix(3,3);
    tensor_pl[i] = block_matrix(3,3);
    tensor_rr[i] = block_matrix(3,3);
  }
  tensor0 = block_matrix(3,3);
  tensor_rl0 = block_matrix(3,3);
  tensor_rl1 = block_matrix(3,3);
  tensor_pl0 = block_matrix(3,3);
  tensor_pl1 = block_matrix(3,3);

  tensor_rQ0 = (double ***) malloc(3 * sizeof(double **));
  tensor_rQ1 = (double ***) malloc(3 * sizeof(double **));
  for(i=0; i < 3; i++) {
    tensor_rQ0[i] = block_matrix(3,3);
    tensor_rQ1[i] = block_matrix(3,3);
  }

  if(compute_pl) {
    sprintf(lbl1, "<<P;L>>_(%5.3f)", 0.0);
    if(!params.restart || !psio_tocscan(CC_INFO, lbl1)) {
      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "P_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.mu_irreps[alpha], 1);
        compute_X(pert, moinfo.mu_irreps[alpha], 0);

        sprintf(pert, "L_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.l_irreps[alpha], 1);
        compute_X(pert, moinfo.l_irreps[alpha], 0);
      }

      fprintf(outfile, "\n\tComputing %s tensor.\n", lbl1); fflush(outfile);
      for(alpha=0; alpha < 3; alpha ++) {
        for(beta=0; beta < 3; beta++) {
          sprintf(pert_x, "P_%1s", cartcomp[alpha]);
          sprintf(pert_y, "L_%1s", cartcomp[beta]);
          linresp(&tensor0[alpha][beta], +1.0, 0.0,
                  pert_x, moinfo.mu_irreps[alpha], 0.0,
                  pert_y, moinfo.l_irreps[beta], 0.0);
        }
      }
      psio_write_entry(CC_INFO, lbl1, (char *) tensor0[0], 9*sizeof(double));
      psio_close(CC_LR, 0);
      psio_open(CC_LR, 0);
      for(j=CC_TMP; j <= CC_TMP11; j++) {
        psio_close(j,0);
        psio_open(j,0);
      }
    }    
    else {
      fprintf(outfile, "Using %s tensor found on disk.\n", lbl1);
      psio_read_entry(CC_INFO, lbl1, (char *) tensor0[0], 9*sizeof(double));
    }

    if (params.wfn == "CC2")
      fprintf(outfile, "\n     CC2 Optical Rotation Tensor (Velocity Gauge): %s\n", lbl1);
    else if(params.wfn == "CCSD")
      fprintf(outfile, "\n    CCSD Optical Rotation Tensor (Velocity Gauge): %s\n", lbl1);

    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    fprintf(outfile,   "   Evaluated at omega = 0.00 E_h (Inf nm, 0.0 eV, 0.0 cm-1)\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor0, 3, 3, outfile);
  }

  for(i=0; i < params.nomega; i++) {
    zero_mat(tensor_rl0, 3, 3);
    zero_mat(tensor_rl1, 3, 3);
    zero_mat(tensor_pl0, 3, 3);
    zero_mat(tensor_pl1, 3, 3);
    for(alpha=0; alpha < 3; alpha++) {
      zero_mat(tensor_rQ0[alpha], 3, 3);
      zero_mat(tensor_rQ1[alpha], 3, 3);
    }

    sprintf(lbl1, "1/2 <<Mu;L>>_(%5.3f)", params.omega[i]);
    sprintf(lbl2, "1/2 <<P;L>>_(%5.3f)", params.omega[i]);
    sprintf(lbl3, "1/2 <<Mu;Mu>>_(%5.3f)", params.omega[i]);
    sprintf(lbl4, "1/2 <<Mu;Q>>_(%5.3f)", params.omega[i]);

    if(!params.restart || 
       ( (compute_rl && !psio_tocscan(CC_INFO,lbl1)) ||
         (compute_pl && !psio_tocscan(CC_INFO,lbl2)) || 
         !psio_tocscan(CC_INFO,lbl3) || 
         !psio_tocscan(CC_INFO,lbl4) 
       )
      ) {

      /* prepare the dipole-length and/or dipole-velocity integrals */
      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "Mu_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.mu_irreps[alpha], 0);
      }

      if(compute_pl) {
        for(alpha=0; alpha < 3; alpha++) {
          sprintf(pert, "P_%1s", cartcomp[alpha]);
	  pertbar(pert, moinfo.mu_irreps[alpha], 1);
        }
      }

      /* prepare the magnetic-dipole integrals */
      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "L_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.l_irreps[alpha], 1);
      }

      /* electric quadrupole integrals */
      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          sprintf(pert, "Q_%1s%1s", cartcomp[alpha], cartcomp[beta]);
          irrep = moinfo.mu_irreps[alpha]^moinfo.mu_irreps[beta];
          pertbar(pert, irrep, 0);
        }
      }

      for(alpha=0; alpha < 3; alpha++) {
        /* -omega electric-dipole CC wave functions */
        sprintf(pert, "Mu_%1s", cartcomp[alpha]);
        compute_X(pert, moinfo.mu_irreps[alpha], -params.omega[i]);

        /* +omega electric-dipole CC wave functions */
        sprintf(pert, "Mu_%1s", cartcomp[alpha]);
        compute_X(pert, moinfo.mu_irreps[alpha], +params.omega[i]);

	if(compute_pl) {
          /* -omega velocity electric-dipole CC wave functions */
          sprintf(pert, "P_%1s", cartcomp[alpha]);
	  compute_X(pert, moinfo.mu_irreps[alpha], -params.omega[i]);
        }

        /* +omega magnetic-dipole CC wave functions */
        sprintf(pert, "L_%1s", cartcomp[alpha]);
	compute_X(pert, moinfo.l_irreps[alpha], +params.omega[i]);
      }

      /* +omega electric-quadrupole CC wave functions */
      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          sprintf(pert, "Q_%1s%1s", cartcomp[alpha], cartcomp[beta]);
          irrep = moinfo.mu_irreps[alpha]^moinfo.mu_irreps[beta];
	  compute_X(pert, irrep, params.omega[i]);
        }
      }

      fprintf(outfile, "\n");
      fprintf(outfile, "\tComputing %s tensor.\n", lbl3); fflush(outfile);
      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
          sprintf(pert_y, "Mu_%1s", cartcomp[beta]);
          linresp(&tensor_rr[i][alpha][beta], -1.0, 0.0, 
                  pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
                  pert_y, moinfo.mu_irreps[beta], +params.omega[i]);
        }
      }
      psio_write_entry(CC_INFO, lbl3, (char *) tensor_rr[0], 9*sizeof(double));

      if(compute_rl) {
        fprintf(outfile, "\tComputing %s tensor.\n", lbl1); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L_%1s", cartcomp[beta]);
            linresp(&tensor_rl0[alpha][beta], -0.5, 0.0, 
                    pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
		    pert_y, moinfo.l_irreps[beta], params.omega[i]);
          }
        }
	psio_write_entry(CC_INFO, lbl1, (char *) tensor_rl0[0], 9*sizeof(double));
      }
      if(compute_pl) {
        fprintf(outfile, "\tComputing %s tensor.\n", lbl2); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "P_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L_%1s", cartcomp[beta]);
	    linresp(&tensor_pl0[alpha][beta], +0.5, 0.0, 
                    pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
		    pert_y, moinfo.l_irreps[beta], params.omega[i]);
          }
        }
	psio_write_entry(CC_INFO, lbl2, (char *) tensor_pl0[0], 9*sizeof(double));
      }

      fprintf(outfile, "\tComputing %s tensor.\n", lbl4); fflush(outfile);
      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          for(gamma=0; gamma < 3; gamma++) {
            sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
            sprintf(pert_y, "Q_%1s%1s", cartcomp[beta], cartcomp[gamma]);
            linresp(&tensor_rQ0[alpha][beta][gamma], -0.5, 0.0, 
                  pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
                  pert_y, moinfo.mu_irreps[beta]^moinfo.mu_irreps[gamma], 
                  params.omega[i]);
          }
        }
      }
      next = PSIO_ZERO;
      for(alpha=0; alpha < 3; alpha++) 
        psio_write(CC_INFO, lbl4, (char *) tensor_rQ0[alpha][0], 
                   9*sizeof(double), next, &next);

      /* Clean up disk space */
      for(j=CC_TMP; j <= CC_TMP11; j++) {
	psio_close(j,0);
	psio_open(j,0);
      }
    }
    else {
      fprintf(outfile, "\n");
      fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl3); fflush(outfile);
      psio_read_entry(CC_INFO, lbl1, (char *) tensor_rr[i][0], 9*sizeof(double));

      if(compute_rl) {
	fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl1); fflush(outfile);
	psio_read_entry(CC_INFO, lbl1, (char *) tensor_rl0[0], 9*sizeof(double));
      }
      if(compute_pl) {
	fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl2); fflush(outfile);
	psio_read_entry(CC_INFO, lbl2, (char *) tensor_pl0[0], 9*sizeof(double));
      }

      fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl4); fflush(outfile);
      next = PSIO_ZERO;
      for(alpha=0; alpha < 3; alpha++)
        psio_read(CC_INFO, lbl4, (char *) tensor_rQ0[alpha][0], 
                  9*sizeof(double), next, & next);

    }

    sprintf(lbl1, "1/2 <<Mu;L*>>_(%5.3f)", params.omega[i]);
    sprintf(lbl2, "1/2 <<P*;L*>>_(%5.3f)", params.omega[i]);
    sprintf(lbl3, "<<Mu;Q>>_(%5.3f)", -params.omega[i]);
    if(!params.restart || 
       ( (compute_rl && !psio_tocscan(CC_INFO,lbl1)) || 
         (compute_pl && !psio_tocscan(CC_INFO,lbl2)) ||
         !psio_tocscan(CC_INFO,lbl3) )
      ) {

      if(compute_pl) {
        for(alpha=0; alpha < 3; alpha++) {
          sprintf(pert, "P*_%1s", cartcomp[alpha]);
	  pertbar(pert, moinfo.mu_irreps[alpha], 1);
        }
      }

      /* prepare the complex-conjugate of the magnetic-dipole integrals */
      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "L*_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.l_irreps[alpha], 1);
      }

      /* +omega velocity electric-dipole CC wave functions */
      for(alpha=0; alpha < 3; alpha++) {
	if(compute_pl) {
          sprintf(pert, "P*_%1s", cartcomp[alpha]);
	  compute_X(pert, moinfo.mu_irreps[alpha], params.omega[i]);
        }

        /* -omega magnetic-dipole CC wave functions */
        sprintf(pert, "L*_%1s", cartcomp[alpha]);
	compute_X(pert, moinfo.l_irreps[alpha], -params.omega[i]);
      }

      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          sprintf(pert, "Q_%1s%1s", cartcomp[alpha], cartcomp[beta]);
          compute_X(pert, moinfo.mu_irreps[alpha]^moinfo.mu_irreps[beta], -params.omega[i]);
        }
      }

      fprintf(outfile, "\n");
      if(compute_rl) {
	fprintf(outfile, "\tComputing %s tensor.\n", lbl1); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L*_%1s", cartcomp[beta]);
	    linresp(&tensor_rl1[alpha][beta], -0.5, 0.0, 
                    pert_x, moinfo.mu_irreps[alpha], params.omega[i],
		    pert_y, moinfo.l_irreps[beta], -params.omega[i]);
          }
        }
	psio_write_entry(CC_INFO, lbl1, (char *) tensor_rl1[0], 9*sizeof(double));
      }
      if(compute_pl) {
	fprintf(outfile, "\tComputing %s tensor.\n", lbl2); fflush(outfile);
        for(alpha=0; alpha < 3; alpha++) {
          for(beta=0; beta < 3; beta++) {
            sprintf(pert_x, "P*_%1s", cartcomp[alpha]);
            sprintf(pert_y, "L*_%1s", cartcomp[beta]);
	    linresp(&tensor_pl1[alpha][beta], +0.5, 0.0, 
                    pert_x, moinfo.mu_irreps[alpha], params.omega[i],
		    pert_y, moinfo.l_irreps[beta], -params.omega[i]);
          }
        }
	psio_write_entry(CC_INFO, lbl2, (char *) tensor_pl1[0], 9*sizeof(double));
      }

      fprintf(outfile, "\tComputing %s tensor.\n", lbl3); fflush(outfile);
      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          for(gamma=0; gamma < 3; gamma++) {
            sprintf(pert_x, "Mu_%1s", cartcomp[alpha]);
            sprintf(pert_y, "Q_%1s%1s", cartcomp[beta], cartcomp[gamma]);
            linresp(&tensor_rQ1[alpha][beta][gamma], -0.5, 0.0, 
                  pert_x, moinfo.mu_irreps[alpha], +params.omega[i],
                  pert_y, moinfo.mu_irreps[beta]^moinfo.mu_irreps[gamma], 
                  -params.omega[i]);
          }
        }
      }

      /* Clean up disk space */
      psio_close(CC_LR, 0);
      psio_open(CC_LR, 0);

      for(j=CC_TMP; j <= CC_TMP11; j++) {
	psio_close(j,0);
	psio_open(j,0);
      }

    }
    else {
      fprintf(outfile, "\n");
      if(compute_rl) {
	fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl1); fflush(outfile);
	psio_read_entry(CC_INFO, lbl1, (char *) tensor_rl1[0], 9*sizeof(double));
      }
      if(compute_pl) {
	fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl2); fflush(outfile);
	psio_read_entry(CC_INFO, lbl2, (char *) tensor_pl1[0], 9*sizeof(double));
      }

      fprintf(outfile, "\tUsing %s tensor found on disk.\n", lbl3); fflush(outfile);
      next = PSIO_ZERO;
      for(alpha=0; alpha < 3; alpha++)
        psio_read(CC_INFO, lbl3, (char *) tensor_rQ1[alpha][0], 
                  9*sizeof(double), next, &next);
    }

    /* sum the two 1/2 tensors for the mixed perturbations */
    for(j=0; j < 3; j++)
      for(k=0; k < 3; k++) {
	if(compute_rl) tensor_rl[i][j][k] = tensor_rl0[j][k] + tensor_rl1[j][k];
	if(compute_pl) tensor_pl[i][j][k] = tensor_pl0[j][k] + tensor_pl1[j][k];
      }

    for(j=0; j < 3; j++)
      for(k=0; k < 3; k++)
        for(l=0; l < 3; l++)
          tensor_rQ[i][j][k][l] = tensor_rQ0[j][k][l] + tensor_rQ1[j][k][l];

    /* Also symmetrize the rr tensor */
    for(j=0; j < 3; j++)
      for(k=0; k < j; k++) {
        if(k!=j) {
          value = 0.5 * (tensor_rr[i][j][k] + tensor_rr[i][k][j]);
          tensor_rr[i][j][k] = tensor_rr[i][k][j] = value;
        }
      }

    if(params.wfn == "CC2")
      fprintf(outfile, "\n                 CC2 Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
    else
      fprintf(outfile, "\n                 CCSD Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");

      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i], _hartree2wavenumbers*params.omega[i]);
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor_rr[i], 3, 3, outfile);

    if(compute_rl) {
      if (params.wfn == "CC2") 
	fprintf(outfile, "\n            CC2 Optical Rotation Tensor (Length Gauge):\n");
      else if(params.wfn == "CCSD")
	fprintf(outfile, "\n           CCSD Optical Rotation Tensor (Length Gauge):\n");

      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i], _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      mat_print(tensor_rl[i], 3, 3, outfile);
    }

    if(compute_pl) {

      if (params.wfn == "CC2") 
	fprintf(outfile, "\n          CC2 Optical Rotation Tensor (Velocity Gauge):\n");
      else if(params.wfn == "CCSD")
	fprintf(outfile, "\n         CCSD Optical Rotation Tensor (Velocity Gauge):\n");

      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i], _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      mat_print(tensor_pl[i], 3, 3, outfile);

      /* subtract the zero-frequency beta tensor */
      for(j=0; j < 3; j++)
	for(k=0; k < 3; k++)
	  tensor_pl[i][j][k] -= tensor0[j][k];

      if (params.wfn == "CC2")
	fprintf(outfile, "\n        CC2 Optical Rotation Tensor (Modified Velocity Gauge):\n");
      else if(params.wfn == "CCSD")
	fprintf(outfile, "\n        CCSD Optical Rotation Tensor (Modified Velocity Gauge):\n");

      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i], _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      mat_print(tensor_pl[i], 3, 3, outfile);

    }

    if(params.wfn == "CC2")
      fprintf(outfile, "\n    CC2 Electric-Dipole/Quadrupole Polarizability [(e^2 a0^2)/E_h]:\n");
    else
      fprintf(outfile, "\n    CCSD Electric-Dipole/Quadrupole Polarizability [(e^2 a0^2)/E_h]:\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i], _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
    for(alpha=0; alpha < 3; alpha++)
      mat_print(tensor_rQ[i][alpha], 3, 3, outfile);

  } /* loop i over nomega */

  for(i=0; i < params.nomega; i++) {
    for(j=0; j < 3; j++) free_block(tensor_rQ[i][j]);
    free(tensor_rQ[i]);
  }
  free(tensor_rQ);
  for(i=0; i < 3; i++) {
    free_block(tensor_rQ0[i]);
    free_block(tensor_rQ1[i]);
  }
  free(tensor_rQ0);
  free(tensor_rQ1);

  for(i=0; i < params.nomega; i++) {
    free_block(tensor_rl[i]);
    free_block(tensor_pl[i]);
    free_block(tensor_rr[i]);
  }
  free(tensor_rl);
  free(tensor_pl);
  free(tensor_rr);
  free_block(tensor0);
  free_block(tensor_rl0);
  free_block(tensor_rl1);
  free_block(tensor_pl0);
  free_block(tensor_pl1);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}

}} // namespace psi::ccresponse
