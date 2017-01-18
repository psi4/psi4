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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
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

void polar(void)
{
  double ***tensor;
  char **cartcomp, pert[32], pert_x[32], pert_y[32];
  int alpha, beta, i;
  double omega_nm, omega_ev, omega_cm, *trace;
  char lbl[32];
  double value;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor = (double ***) malloc(params.nomega * sizeof(double **));
  for(i=0; i < params.nomega; i++)
    tensor[i] = block_matrix(3,3);

  trace = init_array(params.nomega);

  for(i=0; i < params.nomega; i++) {

    sprintf(lbl, "<<Mu;Mu>_(%5.3f)", params.omega[i]);
    if(!params.restart || !psio_tocscan(PSIF_CC_INFO, lbl)) {

      for(alpha=0; alpha < 3; alpha++) {
        sprintf(pert, "Mu_%1s", cartcomp[alpha]);
        pertbar(pert, moinfo.mu_irreps[alpha], 0);
	compute_X(pert, moinfo.mu_irreps[alpha], params.omega[i]);
	if(params.omega[i] != 0.0) compute_X(pert, moinfo.mu_irreps[alpha], -params.omega[i]);
      }

      outfile->Printf( "\n\tComputing %s tensor.\n", lbl);
      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          sprintf(pert_x,"Mu_%1s", cartcomp[alpha]);
          sprintf(pert_y,"Mu_%1s", cartcomp[beta]);
          linresp(&tensor[i][alpha][beta], -1.0, 0.0,
                  pert_x, moinfo.mu_irreps[alpha], -params.omega[i],
                  pert_y, moinfo.mu_irreps[beta], params.omega[i]);
        }
      }

      psio_write_entry(PSIF_CC_INFO, lbl, (char *) tensor[i][0], 9*sizeof(double));

      psio_close(PSIF_CC_LR, 0);
      psio_open(PSIF_CC_LR, 0);
    }
    else {
      outfile->Printf( "Using %s tensor found on disk.\n", lbl);
      psio_read_entry(PSIF_CC_INFO, lbl, (char *) tensor[i], 9*sizeof(double));
    }

    /* symmetrize the polarizability */
    for(alpha=0; alpha < 3; alpha++)
      for(beta=0; beta < alpha; beta++) {
        if(alpha!=beta) {
          value = 0.5 * (tensor[i][alpha][beta] + tensor[i][beta][alpha]);
          tensor[i][alpha][beta] = value;
        }
      }

    if (params.wfn == "CC2")
      outfile->Printf( "\n                 CC2 Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
    else
      outfile->Printf( "\n                 CCSD Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
    outfile->Printf( "  -------------------------------------------------------------------------\n");
    if(params.omega[i] != 0.0)
      omega_nm = (pc_c*pc_h*1e9)/(pc_hartree2J*params.omega[i]);
    omega_ev = pc_hartree2ev*params.omega[i];
    omega_cm = pc_hartree2wavenumbers*params.omega[i];
    if(params.omega[i] != 0.0)
      outfile->Printf(   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n",
	      params.omega[i], omega_nm, omega_ev, omega_cm);
    else
      outfile->Printf(   "   Evaluated at omega = %8.6f E_h (Inf nm, %5.3f eV, %8.2f cm-1)\n",
	      params.omega[i], omega_ev, omega_cm);
    outfile->Printf( "  -------------------------------------------------------------------------\n");
    mat_print(tensor[i], 3, 3, "outfile");
    trace[i] = tensor[i][0][0] + tensor[i][1][1] + tensor[i][2][2];
    outfile->Printf( "\n\talpha_(%5.3f) = %20.12f a.u.\n", params.omega[i], trace[i]/3.0);

    // Need to generalize this for multi-wavelength calcs
    if (params.wfn == "CC2")
      Process::environment.globals["CC2 DIPOLE POLARIZABILITY"] = trace[i]/3.0;
    else
      Process::environment.globals["CCSD DIPOLE POLARIZABILITY"] = trace[i]/3.0;
  }

  if(params.nomega > 1) {  /* print a summary table for multi-wavelength calcs */

    outfile->Printf( "\n\t-------------------------------\n");
    if (params.wfn == "CC2")
      outfile->Printf(   "\t      CC2 Polarizability\n");
    else
      outfile->Printf(   "\t      CCSD Polarizability\n");
    outfile->Printf(   "\t-------------------------------\n");
    outfile->Printf(   "\t    Omega          alpha\n");
    outfile->Printf(   "\t E_h      nm        a.u.        \n");
    outfile->Printf(   "\t-----   ------ ----------------\n");
    for(i=0; i < params.nomega; i++)
      outfile->Printf( "\t%5.3f   %6.2f      %10.5f\n", params.omega[i], (pc_c*pc_h*1e9)/(pc_hartree2J*params.omega[i]), trace[i]/3.0);
  }

  for(i=0; i < params.nomega; i++)
    free_block(tensor[i]);
  free(tensor);

  free(trace);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}

}} // namespace psi::ccresponse
