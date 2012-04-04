/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libdpd/dpd.h>
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
    fprintf(outfile, "\nread_guess(): Unable to read number of guesses from input.\n");
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
        fprintf(outfile,"only alpha guesses allowed for EOM_REF = RHF\n");
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
            fprintf(outfile, "\nInconsisent symmetries in components of guess %d.\n", k);
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
              fprintf(outfile, "\nInconsisent symmetries in components of guess %d.\n", k);
              exit(2);
            }
          }
          else {
            if (moinfo.bocc_sym[i]^moinfo.bvir_sym[a] != this_irrep) {
              fprintf(outfile, "\nInconsisent symmetries in components of guess %d.\n", k);
              exit(2);
            }
          }
        }
      }
    }
  }

  fprintf(outfile,"EOM_GUESS_VECTORS implies roots_per_irrep: \n\t");
  for (i=0;i<moinfo.nirreps;++i)
    fprintf(outfile,"%s %d, ",moinfo.irr_labs[i], eom_params.states_per_irrep[i]);
  fprintf(outfile,"\n");
  fprintf(outfile,"and Rs_per_irrep: \n\t");
  for (i=0;i<moinfo.nirreps;++i)
    fprintf(outfile,"%s %d, ",moinfo.irr_labs[i], eom_params.cs_per_irrep[i]);
  fprintf(outfile,"These numbers should match those given by the roots_per_irrep keyword\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);
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
    fprintf(outfile, "\nread_guess(): Unable to read number of guesses from input.\n");
    exit(2);
  }

  /* loop over number of initial guess of this symmetry */
  for(k=0; k < eom_params.cs_per_irrep[C_irr]; k++) {
    sprintf(lbl, "%s %d", "CME", k);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_scm(&CME, 0);
    dpd_file2_mat_init(&CME);

    if(params.eom_ref <= 1) {
      sprintf(lbl, "%s %d", "Cme", k);
      dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
      dpd_file2_scm(&Cme, 0);
      dpd_file2_mat_init(&Cme);
    }
    else if (params.eom_ref == 2) {
      sprintf(lbl, "%s %d", "Cme", k);
      dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
      dpd_file2_scm(&Cme, 0);
      dpd_file2_mat_init(&Cme);
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
            fprintf(outfile, "\nInconsisent symmetries in components of guess %d.\n", k);
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
              fprintf(outfile, "\nInconsisent symmetries in components of guess %d.\n", k);
              exit(2);
            }
          }
          else {
            if (moinfo.bocc_sym[i]^moinfo.bvir_sym[a] != this_irrep) {
              fprintf(outfile, "\nInconsisent symmetries in components of guess %d.\n", k);
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
    dpd_file2_mat_wrt(&CME);
    dpd_file2_mat_wrt(&Cme);
    dpd_file2_mat_close(&CME);
    dpd_file2_mat_close(&Cme);

    dpd_file2_scm(&CME,1.0/sqrt(norm));
    dpd_file2_scm(&Cme,1.0/sqrt(norm));

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }

  return;
}

}} // namespace psi::cceom
