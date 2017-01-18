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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libipv1/ip_lib.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* counts the number of states and Cs for each irrep based on the orbital numbers given
 * in input -- overrides the "roots_per_irrep" keyword */

void read_guess_init(void)
{
  int i, a, k, l, spin, errcod, C_irrep;
  int num_vectors, vector_len, this_irrep, useit;
  char lbl[32];
  double norm, value;
  dpdfile2 CME;

  for (i=0;i<moinfo.nirreps;++i) {
    eom_params.cs_per_irrep[i] = 0;
    eom_params.states_per_irrep[i] = 0;
  }

  /* Read number of guess = number of final states to solve for */
  errcod = ip_count("EOM_GUESS_VECTORS",&num_vectors,0);
  if(errcod != IPE_OK) {
    outfile->Printf( "\nread_guess(): Unable to read number of guesses from input.\n");
    exit(2);
  }

  for(k=0; k < num_vectors; k++) {
    ip_count("EOM_GUESS_VECTORS", &vector_len, 1, k);

    for(l=0; l < vector_len; l++) {
      errcod = ip_data("EOM_GUESS_VECTORS", "%d", &i, 3, k, l, 0);
      errcod = ip_data("EOM_GUESS_VECTORS", "%d", &a, 3, k, l, 1);
      errcod = ip_data("EOM_GUESS_VECTORS", "%lf", &value, 3, k, l, 2);
      errcod = ip_data("EOM_GUESS_VECTORS", "%d", &spin, 4, k, l, 3);

      if ((spin != 0) && (params.eom_ref == 0)) {
        outfile->Printf("only alpha guesses allowed for EOM_REF = RHF\n");
        exit(1);
      }

      if (params.eom_ref == 0) {
        if(l==0) { /* check symmetry of first excitation */
          this_irrep = moinfo.occ_sym[i]^moinfo.vir_sym[a];
          eom_params.cs_per_irrep[this_irrep] += 1;
          eom_params.states_per_irrep[this_irrep^moinfo.sym] += 1;
        }
        else { /* check consistency of other excitations */
          if (moinfo.occ_sym[i]^moinfo.vir_sym[a] != this_irrep) {
            outfile->Printf( "\nInconsisent symmetries in components of guess %d.\n", k);
            exit(2);
          }
        }
      }
      else {
        if(l==0) { /* check symmetry of first excitation */
          if (spin == 0)
            this_irrep = moinfo.aocc_sym[i]^moinfo.avir_sym[a];
          else
            this_irrep = moinfo.bocc_sym[i]^moinfo.bvir_sym[a];

	      eom_params.cs_per_irrep[this_irrep] += 1;
	      eom_params.states_per_irrep[this_irrep^moinfo.sym] += 1;
        }
        else { /* check consistency of other excitations */
          if (spin == 0) {
            if (moinfo.aocc_sym[i]^moinfo.avir_sym[a] != this_irrep) {
              outfile->Printf( "\nInconsisent symmetries in components of guess %d.\n", k);
              exit(2);
            }
          }
          else {
            if (moinfo.bocc_sym[i]^moinfo.bvir_sym[a] != this_irrep) {
              outfile->Printf( "\nInconsisent symmetries in components of guess %d.\n", k);
              exit(2);
            }
          }
        }
      }
    }
  }

  outfile->Printf("EOM_GUESS_VECTORS implies roots_per_irrep: \n\t");
  for (i=0;i<moinfo.nirreps;++i)
    outfile->Printf("%s %d, ",moinfo.irr_labs[i], eom_params.states_per_irrep[i]);
  outfile->Printf("\n");
  outfile->Printf("and Rs_per_irrep: \n\t");
  for (i=0;i<moinfo.nirreps;++i)
    outfile->Printf("%s %d, ",moinfo.irr_labs[i], eom_params.cs_per_irrep[i]);
  outfile->Printf("These numbers should match those given by the roots_per_irrep keyword\n");
  outfile->Printf("\n\n");

  return;
}

void read_guess(int C_irr)
{
  int i, a, k, l, spin, errcod;
  int num_vectors, vector_len, this_irrep;
  char lbl[32];
  double norm, value;
  dpdfile2 CME, Cme;

  errcod = ip_count("EOM_GUESS_VECTORS",&num_vectors,0);
  if(errcod != IPE_OK) {
    outfile->Printf( "\nread_guess(): Unable to read number of guesses from input.\n");
    exit(2);
  }

  /* loop over number of initial guess of this symmetry */
  for(k=0; k < eom_params.cs_per_irrep[C_irr]; k++) {
    sprintf(lbl, "%s %d", "CME", k);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
    global_dpd_->file2_scm(&CME, 0);
    global_dpd_->file2_mat_init(&CME);

    if(params.eom_ref <= 1) {
      sprintf(lbl, "%s %d", "Cme", k);
      global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
      global_dpd_->file2_scm(&Cme, 0);
      global_dpd_->file2_mat_init(&Cme);
    }
    else if (params.eom_ref == 2) {
      sprintf(lbl, "%s %d", "Cme", k);
      global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
      global_dpd_->file2_scm(&Cme, 0);
      global_dpd_->file2_mat_init(&Cme);
    }

    norm = 0.0;
    ip_count("EOM_GUESS_VECTORS", &vector_len, 1, k);
    for(l=0; l < vector_len; l++) {
      errcod = ip_data("EOM_GUESS_VECTORS", "%d", &i, 3, k, l, 0);
      errcod = ip_data("EOM_GUESS_VECTORS", "%d", &a, 3, k, l, 1);
      errcod = ip_data("EOM_GUESS_VECTORS", "%lf", &value, 3, k, l, 2);
      errcod = ip_data("EOM_GUESS_VECTORS", "%d", &spin, 3, k, l, 3);

      if (params.ref == 0) {
        if(l==0) { /* check symmetry of this state */
          this_irrep = moinfo.occ_sym[i]^moinfo.vir_sym[a];
        }
        else{ /* check other excitations for consistency */
          if (moinfo.occ_sym[i]^moinfo.vir_sym[a] != this_irrep) {
            outfile->Printf( "\nInconsisent symmetries in components of guess %d.\n", k);
            exit(2);
          }
        }
      }
      else {
        if(l==0) { /* check symmetry of this state */
          if (spin == 0)
            this_irrep = moinfo.aocc_sym[i]^moinfo.avir_sym[a];
          else
            this_irrep = moinfo.bocc_sym[i]^moinfo.bvir_sym[a];
        }
        else{ /* check other excitations for consistency */
          if (spin == 0) {
            if (moinfo.aocc_sym[i]^moinfo.avir_sym[a] != this_irrep) {
              outfile->Printf( "\nInconsisent symmetries in components of guess %d.\n", k);
              exit(2);
            }
          }
          else {
            if (moinfo.bocc_sym[i]^moinfo.bvir_sym[a] != this_irrep) {
              outfile->Printf( "\nInconsisent symmetries in components of guess %d.\n", k);
              exit(2);
            }
          }
        }
      }

      if (spin == 0)
        CME.matrix[C_irr][i][a] = value;
      else
        Cme.matrix[C_irr][i][a] = value;

      norm += value * value;
    }
    global_dpd_->file2_mat_wrt(&CME);
    global_dpd_->file2_mat_wrt(&Cme);
    global_dpd_->file2_mat_close(&CME);
    global_dpd_->file2_mat_close(&Cme);

    global_dpd_->file2_scm(&CME,1.0/sqrt(norm));
    global_dpd_->file2_scm(&Cme,1.0/sqrt(norm));

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }

  return;
}

}} // namespace psi::cceom
