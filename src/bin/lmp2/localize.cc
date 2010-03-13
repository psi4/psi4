/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
*/

#include <iostream>
#include <fstream>              // file I/O support
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
#include <cmath>
#include <libiwl/iwl.hpp>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
//#include <libpsio/psio.hpp>
//#include <libqt/qt.h>
#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>
#include <libmints/wavefunction.h>
#include <libparallel/parallel.h>
#include <libint/libint.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{

namespace lmp2{

void LMP2::localize() {

  using namespace psi;

  if(Communicator::world->me() == 0) {
    fprintf(outfile, "\n********************* Entering Localization Scope ********************************\n");
  }

  int iter, s, t, A, k, l, iold, max;
  int i, j, ij, am, atom, shell_length, offset;
  int ntri, puream, *soccpi, *stype, *snuc, nfzc;;
  double  **LCtmp, **F_occ;
  double *scratch, *evals;
  int *orb_order, *orb_boolean;
  double P, PiiA, Pst, Pss, Ptt, Ast, Bst, AB;
  double Uss, Utt, Ust, Uts, LCks, LCkt, **U, **V, **VV;
  double cos4a, alpha, alphamax, alphalast, conv;

  if(print > 2 && Communicator::world->me() == 0) {
    fprintf(outfile, "\nC Matrix in the AO basis:\n");
    print_mat(C, nso, nso, outfile);
  }

  evals = get_evals();
  puream = get_puream();
  snuc = get_snuc();

  ntri = nso*(nso+1)/2;
  scratch = init_array(ntri);

  // **** Read in the overlap matrix ****
  aoovlp = block_matrix(nso, nso);

//  shared_ptr<BasisSet> basis(new BasisSet(chkpt));
//  shared_ptr<IntegralFactory> integral(new IntegralFactory(basis, basis, basis, basis));
//  shared_ptr<OneBodyInt> S(integral->overlap());
//  shared_ptr<Matrix> overlap(factory->create_matrix("Overlap"));

//  S->compute(overlap);

    if(Communicator::world->me() == 0) {
        IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_S, scratch, ntri, 0, 0, outfile);
    }
    Communicator::world->bcast(scratch, ntri, 0);
    for(i=0, ij=0; i < nso; i++)
        for(j=0; j <= i; j++, ij++)
            aoovlp[i][j] = aoovlp[j][i] = scratch[ij];
  
  //send_overlap(aoovlp);

  free(scratch);

  if(print > 2 && Communicator::world->me() == 0) {
    fprintf(outfile, "Overlap Matrix");
    print_mat(aoovlp, nso, nso, outfile);
  }

  soccpi = get_soccpi();
  if(soccpi[0]) {
    throw PsiException("Local MP2 is only valid for closed-shell molecules.", __FILE__, __LINE__);
  }
  free(soccpi);

  // Frozen orbital info
  nfzc = get_frdocc();

  // Compute the length of each AM block
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=1; l < (LIBINT_MAX_AM); l++) {
    if(puream) l_length[l] = 2 * l + 1;
    else l_length[l] = l_length[l-1] + l + 1;
  }

  // Set up the atom->AO and AO->atom lookup arrays
  aostart = init_int_array(natom);
  aostop = init_int_array(natom);
  stype = get_stype();
  for(i=0,atom=-1,offset=0; i < nshell; i++) {
    am = stype[i] - 1;                  // am is the angular momentum of the orbital
    shell_length = l_length[am];        // shell_length is the number of obritals in each shell

    if(atom != snuc[i]-1) {             // snuc is the nucleus that the shell belongs to
      if(atom != -1) aostop[atom] = offset-1;
      atom = snuc[i]-1;
      aostart[atom] = offset;
    }

    offset += shell_length;
  }
  aostop[atom] = offset-1;

  ao2atom = init_int_array(nso);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) {
      ao2atom[j] = i;                   // ao2atom is the atom number that the AO is located on
    }

  if(Communicator::world->me() == 0) {
    fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n\n", nocc);

    fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv\n");
    fprintf(outfile, "\t------------------------------------------------------------\n");
  }

  U = block_matrix(nocc, nocc);
  V = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++) V[i][i] = 1.0;
  VV = block_matrix(nocc, nocc);

  for(iter=0; iter < 100; iter++) {

  P = 0.0;
  for(i=nfzc; i < nocc; i++) {
    for(A=0; A < natom; A++) {
      PiiA = 0.0;

      for(l=aostart[A]; l <= aostop[A]; l++)
        for(k=0; k < nso; k++)
          PiiA += C[k][i] * C[l][i] * aoovlp[k][l];

      P += PiiA * PiiA;
    }
  }

  // Compute 2x2 rotations for Pipek-Mezey lo.zation
  alphamax = 0.0;

  for(s=nfzc; s < nocc; s++) {
    for(t=nfzc; t < s; t++) {

      Ast = Bst = 0.0;

      for(A=0; A < natom; A++) {

        Pst = Pss = Ptt = 0.0;

        for(l=aostart[A]; l <= aostop[A]; l++) {
          for(k=0; k < nso; k++) {
            Pst += 0.5 * (C[k][s] * C[l][t] +
                          C[l][s] * C[k][t]) * aoovlp[k][l];            // Eqn 31 (JCP 90, 4916)

            Pss += C[k][s] * C[l][s] * aoovlp[k][l];                    // Eqn 31 (JCP 90, 4916)

            Ptt += C[k][t] * C[l][t] * aoovlp[k][l];                    // Eqn 31 (JCP 90, 4916)
          }
        }

        Ast += Pst * Pst - 0.25 * (Pss - Ptt) * (Pss - Ptt);               // Eqn 29A (JCP 90, 4916)
        Bst += Pst * (Pss - Ptt);                                  // Eqn 29B (JCP 90, 4916)

      } // A-loop

      // Compute the rotation angle
      AB = Ast * Ast + Bst * Bst;

        if(fabs(AB) > 0.0) {
          cos4a = -Ast/sqrt(AB);                                     // Eqn 13b (JCP 90, 4916)
          alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
        }
        else alpha = 0.0;

        // Keep up with the maximum 2x2 rotation angle
        alphamax = (fabs(alpha) > alphamax ? alpha : alphamax);


        Uss = cos(alpha);                                            // Eqn 10a/b (JCP 90, 4916)
        Utt = cos(alpha);                                               // Eqn 10a/b (JCP 90, 4916)
        Ust = sin(alpha);                                               // Eqn 10a/b (JCP 90, 4916)
        Uts = -Ust;                                                  // Eqn 10a/b (JCP 90, 4916)

        // Now do the rotation
        for(k=0; k < nso; k++) {
          LCks = C[k][s];
          LCkt = C[k][t];
          C[k][s] = Uss * LCks + Ust * LCkt;
          C[k][t] = Uts * LCks + Utt * LCkt;
        }

        zero_mat(U, nocc, nocc);
        for(i=0; i < nocc; i++) U[i][i] = 1.0;

        U[s][s] = Uss;
        U[t][t] = Utt;
        U[s][t] = Ust;
        U[t][s] = Uts;

        zero_mat(VV, nocc, nocc);
        for(i=0; i < nocc; i++) {
          for(j=0; j < nocc; j++) {
            for(k=0; k < nocc; k++) {
              VV[i][j] += V[i][k] * U[j][k];
            }
          }
        }

        for(i=0; i < nocc; i++)
          for(j=0; j < nocc; j++)
            V[i][j] = VV[i][j];

      } // t-loop
    } // s-loop

    conv = fabs(alphamax) - fabs(alphalast);
    if(Communicator::world->me() == 0) {
      fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, P, alphamax, conv);
    }
    if((iter > 2) && ((fabs(conv) < 1e-12) || alphamax == 0.0)) break;
    alphalast = alphamax;

    fflush(outfile);

  } // iter-loop

  // Transform occupied orbital eigenvalues
  F_occ = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nocc; j++)
      for(k=0; k < nocc; k++)
        F_occ[i][j] += V[k][i] * evals[k] * V[k][j];

  // Compute a reordering array based on the diagonal elements of F
  orb_order = init_int_array(nocc);
  orb_boolean = init_int_array(nocc);
  for(i=0; i < nocc; i++) { orb_order[i] = 0;  orb_boolean[i] = 0; }

  for(i=0,max=0; i < nocc; i++) // First, find the overall maximum
    if(fabs(F_occ[i][i]) > fabs(F_occ[max][max])) max = i;

  orb_order[0] = max;  orb_boolean[max] = 1;

  for(i=1; i < nocc; i++) {
    max = 0;
    while(orb_boolean[max]) max++; // Find an unused max
    for(j=0; j < nocc; j++)
      if((fabs(F_occ[j][j]) >= fabs(F_occ[max][max])) && !orb_boolean[j]) max = j;
    orb_order[i] = max; orb_boolean[max] = 1;
  }

  // Now reorder the localized MO's according to F
  LCtmp = block_matrix(nso,nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nso; j++) LCtmp[j][i] = C[j][i];

  for(i=0; i < nocc; i++) {
    iold = orb_order[i];
    for(j=0; j < nso; j++) C[j][i] = LCtmp[j][iold];
    evals[i] = F_occ[iold][iold];
  }
  free_block(LCtmp);

  if(print > 3 && Communicator::world->me() == 0) {
    fprintf(outfile, "\nC Matrix in the LO basis:\n");
    print_mat(C, nso, nso, outfile);
  }

  free(evals);

  if(Communicator::world->me() == 0) {
    fprintf(outfile, "\n********************* Exiting Localization Scope ********************************\n");
  }
}

}} // namespace psi::lmp2

