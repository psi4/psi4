/*! \file
    \ingroup LMP2
    \brief Construct MO integrals from density-fitted 3-center quantities
 */

#include <psi4-dec.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <iostream>
#include <fstream>              // file I/O support
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libparallel/parallel.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

#define TIME_DF_LMP2 1

namespace psi {

namespace lmp2 {

void LMP2::direct_df_transformation() {
#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0){
  timer_init();
  timer_on("Compute DF-LMP2");
}
#endif

  using namespace std;

  // Create a basis set object and initialize it using the checkpoint file.
    // Create a basis set object and initialize it.
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, Process::environment.molecule(), orbital_basis);
    shared_ptr<BasisSet> ribasis = BasisSet::construct(parser, Process::environment.molecule(), ri_basis);
  //ribasis->print();

  shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

  // Create integral factory
  IntegralFactory rifactory(ribasis, zero, basis, basis);
  IntegralFactory rifactory_J(ribasis, zero, ribasis, zero);

  // Create an integral object for ERIs

  TwoBodyInt* eri = rifactory.eri();
  TwoBodyInt* Jint = rifactory_J.eri();
  double **J = block_matrix(ribasis->nbf(), ribasis->nbf());
  double **J_mhalf = block_matrix(ribasis->nbf(), ribasis->nbf());
  const double *Jbuffer = Jint->buffer();

#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0)  timer_on("Form J");
#endif

  int index = 0;
  for(int MU = 0; MU < ribasis->nshell(); ++MU) {
    int nummu = ribasis->shell(MU)->nfunction();

    for(int NU = 0; NU < ribasis->nshell(); ++NU) {
      int numnu = ribasis->shell(NU)->nfunction();

      Jint->compute_shell(MU, 0, NU, 0);

      for(int mu = 0, index=0; mu < nummu; ++mu) {
        int omu = ribasis->shell(MU)->function_index() + mu;

        for(int nu = 0; nu < numnu; ++nu, ++index) {
          int onu = ribasis->shell(NU)->function_index() + nu;

          J[omu][onu] = Jbuffer[index];
        }
      }
    }
  }

  // Form J^-1/2
  // First, diagonalize J
  // the C_DSYEV call replaces the original matrix J with its eigenvectors
  double* eigval = init_array(ribasis->nbf());
  int lwork = ribasis->nbf() * 3;
  double* work = init_array(lwork);
  int stat = C_DSYEV('v', 'u', ribasis->nbf(), J[0], ribasis->nbf(), eigval,
          work, lwork);
  if(stat != 0) {
    fprintf(outfile, "C_DSYEV failed\n");
    exit(PSI_RETURN_FAILURE);
  }
  free(work);

  // Now J contains the eigenvectors of the original J
  // Copy J to J_copy
  double **J_copy = block_matrix(ribasis->nbf(), ribasis->nbf());
  C_DCOPY(ribasis->nbf() * ribasis->nbf(), J[0], 1, J_copy[0], 1);

  // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
  // where j^{-1/2} is the diagonal matrix of the inverse square roots
  // of the eigenvalues, and U is the matrix of eigenvectors of J
  for(int i = 0; i < ribasis->nbf(); i++) {
    if(eigval[i] < 1.0E-10)
      eigval[i] = 0.0;
    else {
      eigval[i] = 1.0 / sqrt(eigval[i]);
      }
    // scale one set of eigenvectors by the diagonal elements j^{-1/2}
    C_DSCAL(ribasis->nbf(), eigval[i], J[i], 1);
  }
  free(eigval);

  // J_mhalf = J_copy(T) * J
  C_DGEMM('t', 'n', ribasis->nbf(), ribasis->nbf(), ribasis->nbf(), 1.0,
          J_copy[0], ribasis->nbf(), J[0], ribasis->nbf(), 0.0, J_mhalf[0], ribasis->nbf());

  free_block(J);
  free_block(J_copy);

#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0)  timer_off("Form J");
#endif

  int i, j, ij, l, k, a, b, t, u, v;

  int *ij_owner, *ij_local;
  int **ij_map, *abs_ij_map;
  int **pairdomain, *pairdom_len;

  ij_owner = get_ij_owner();
  ij_local = get_ij_local();

  ij_map = get_ij_map();
  abs_ij_map = original_ij_map();

  pairdomain = compute_pairdomain(ij_map);
  pairdom_len = compute_pairdomlen(ij_map);

  /* Build the "united pair domains" (Werner JCP 118, 8149 (2003) */
  /* This will only work for now, where we have assumed that all
     pairs are strong pairs.  This assumption will change, then
     we may need to change this code --CDS 12/09
   */
  int **uniteddomain = init_int_matrix(nocc, natom);
  int *uniteddomain_len = init_int_array(nocc);
  for(i = 0; i < nocc; i++) {
    for(k = 0; k < natom; k++) {
      for(j = 0; j < nocc; j++) {
        ij = INDEX(i, j);
        ij = abs_ij_map[ij]; // map org ij to new ij
        if(pairdomain[ij][k] && uniteddomain[i][k] == 0) {
          uniteddomain[i][k] = 1;
          uniteddomain_len[i] += aostop[k] - aostart[k] + 1;
        }
      }
    }
  }

  int **uniteddomain_abs2rel = init_int_matrix(nocc, nso);
  for(i = 0; i < nocc; i++) {
    for(k = 0, a = 0; k < natom; k++) {
      if(uniteddomain[i][k]) {
        for(t = aostart[k]; t <= aostop[k]; t++, a++) {
          uniteddomain_abs2rel[i][t] = a;
        }
      }
    }
  }

  /* This suppose to map local i and j for given local ij
   * Here *local* mean local to the given node
   */

  int **proc_has_i = init_int_matrix(Communicator::world->nproc(),nocc);
  int **i_local = init_int_matrix(Communicator::world->nproc(),nocc);
  int *i_local_count = init_int_array(Communicator::world->nproc());
  int proc, count;

  for(ij = 0; ij < ij_pairs; ij++) {
    proc = ij_owner[ij];
    i = ij_map[ij][0];
    j = ij_map[ij][1];
    proc_has_i[proc][i] = 1;
    proc_has_i[proc][j] = 1;
  }

  for(proc = 0; proc < Communicator::world->nproc(); proc++) {
    for(i = 0, count = 0; i < nocc; i++, count++) {
      if(proc_has_i[proc][i]) {
        //i_list[proc][count] = i;
        i_local[proc][i] = count;
      }
    }
    i_local_count[proc] = count;
  }



  double ***mo_p_ir = (double ***) malloc(sizeof (double **) * nocc);
  for(i = 0; i < nocc; i++)
    mo_p_ir[i] = (double **) malloc(sizeof (double *) * uniteddomain_len[i]);
  for(i = 0; i < nocc; i++) {
    for(j = 0; j < uniteddomain_len[i]; j++) {
      mo_p_ir[i][j] = (double *) malloc(sizeof (double) * ribasis->nbf());
      memset(mo_p_ir[i][j], '\0', sizeof (double) * ribasis->nbf());
    }
  }

  double ***B_ir_p = (double ***) malloc(sizeof (double **) * nocc);
  for(i = 0; i < nocc; i++)
    B_ir_p[i] = (double **) malloc(sizeof (double *) * uniteddomain_len[i]);
  for(i = 0; i < nocc; i++){
    for(j = 0; j < uniteddomain_len[i]; j++) {
      B_ir_p[i][j] = (double *) malloc(sizeof (double) * ribasis->nbf());
      memset(B_ir_p[i][j], '\0', sizeof (double) * ribasis->nbf());
    }
  }

  IntegralFactory ao_eri_factory(basis, basis, basis, basis);
  TwoBodyInt* ao_eri = ao_eri_factory.eri();
  const double *ao_buffer = ao_eri->buffer();

  double *Schwartz = init_array(basis->nshell() * (basis->nshell()+1) / 2);
  double *DFSchwartz = init_array(ribasis->nshell());

  for(int P=0,PQ=0;P<basis->nshell();P++) {
    int numw = basis->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basis->shell(Q)->nfunction();
      double tei, max=0.0;

      ao_eri->compute_shell(P, Q, P, Q);

      for(int w=0;w<numw;w++) {
        for(int x=0;x<numx;x++) {
          index = ( ( (w*numx + x) * numw + w) * numx + x);
          tei = ao_buffer[index];
          if(fabs(tei) > max) max = fabs(tei);
        }
      }
      Schwartz[PQ] = max;
    }
  }

  for(int P=0;P<ribasis->nshell();P++) {
    int numw = ribasis->shell(P)->nfunction();
    double tei, max=0.0;

    Jint->compute_shell(P, 0, P, 0);

    for(int w=0;w<numw;w++) {
      tei = Jbuffer[w];
      if(fabs(tei) > max) max = fabs(tei);
    }
    DFSchwartz[P] = max;
  }



  double **half = block_matrix(nocc, nso);

  int numPshell, Pshell, MU, NU, P, oP, Q, oQ, mu, nu, nummu, numnu, omu, onu, mn;

  // find out the max number of P's in a P shell
  int screened=0;
  int maxPshell = 0;
  for(Pshell = 0; Pshell < ribasis->nshell(); ++Pshell) {
    numPshell = ribasis->shell(Pshell)->nfunction();
    if (numPshell > maxPshell) maxPshell = numPshell;
  }

  double*** temp = new double**[maxPshell];
  for (P = 0; P < maxPshell; P++) temp[P] = block_matrix(nso, nso);

#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0) timer_on("Form mo_p_ir");
#endif

  const double *buffer = eri->buffer();
  for(Pshell = 0; Pshell < ribasis->nshell(); ++Pshell) {
    numPshell = ribasis->shell(Pshell)->nfunction();

    for(P = 0; P < numPshell; ++P) {
      zero_mat(temp[P], nso, nso);
    }

    for(MU = 0, mn=0; MU < basis->nshell(); ++MU) {
      nummu = basis->shell(MU)->nfunction();
      for(NU = 0; NU <= MU; ++NU,++mn) {
        numnu = basis->shell(NU)->nfunction();

        if (sqrt(Schwartz[mn]*DFSchwartz[Pshell])>tol) {
        eri->compute_shell(Pshell, 0, MU, NU);

        for(P = 0, index = 0; P < numPshell; ++P) {
          for(mu = 0; mu < nummu; ++mu) {
            omu = basis->shell(MU)->function_index() + mu;
            for(nu = 0; nu < numnu; ++nu, ++index) {
              onu = basis->shell(NU)->function_index() + nu;
              temp[P][omu][onu] = buffer[index]; // (oP | omu onu) integral
              temp[P][onu][omu] = buffer[index]; // (oP | onu omu) integral
            }
          }
        } // end loop over P in Pshell

        } // end Schwartz inequality
        else screened++;
      } // end loop over NU shell
    } // end loop over MU shell
    // now we've gone through all P, mu, nu for a given Pshell
    // transform the integrals for all P in the given P shell

    for(P = 0; P < numPshell; ++P) {
      oP = ribasis->shell(Pshell)->function_index() + P;

      // Do transform
      // (inu|A) = sum_mu (munu|A)C_mui
      C_DGEMM('T', 'N', nocc, nso, nso, 1.0, &(C[0][0]),
              nso, &(temp[P][0][0]), nso, 0.0, &(half[0][0]), nso);

      for(i = 0; i < nocc; i++) {
        for(k = 0, a = 0; k < natom; k++) {
          if(uniteddomain[i][k]) {
            for(t = aostart[k]; t <= aostop[k]; t++, a++) {
              for(mu = 0; mu < nso; ++mu) {
                 mo_p_ir[i][a][oP] += Rt_full[mu][t] * half[i][mu];
              }
            }
          }
        }
      }
    } // loop over the functions within the P shell
  } // end loop over P shells; done with forming MO basis (P|ia)'s

  fprintf(outfile,"  %d shell triplets screened via Schwartz inequality\n\n",
    screened);

  // should free temp here
  for(P = 0; P < maxPshell; P++) free_block(temp[P]);
  // destruct temp[] itself?

#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0) timer_off("Form mo_p_ir");
#endif

// fprintf(outfile, "mo_p_ia:\n");
// print_mat(mo_p_ia, ribasis->nbf(), nact_docc*nact_virt, outfile);

#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0) timer_on("Form B_ir^P");
#endif

// mo_p_ir has integrals
// B_ir^P = Sum_Q (i r | Q) (J^-1/2)_QP

  for (i = 0; i < nocc; i++) {
      for (a = 0; a < uniteddomain_len[i]; a++) {
          for (oP = 0; oP < ribasis->nbf(); ++oP) {
              B_ir_p[i][a][oP] = C_DDOT(ribasis->nbf(), &(mo_p_ir[i][a][0]), 1, &(J_mhalf[oP][0]), 1);
          }
      }
  }

  for(i = 0; i < nocc; i++)
    for(j = 0; j < uniteddomain_len[i]; j++)
      free(mo_p_ir[i][j]);
  for(i = 0; i < nocc; i++)
    free(mo_p_ir[i]);
  free(mo_p_ir);

#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0) timer_off("Form B_ir^P");
#endif

  //Allocating the memory needed by Ktilde
  if(ij_pairs % Communicator::world->nproc() == 0) {
    Ktilde = (double ***) malloc((ij_pairs / Communicator::world->nproc()) * sizeof (double **));
    for(ij = 0; ij < ij_pairs; ij++) {
      if(Communicator::world->me() == ij_owner[ij])
        Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
    }
  }
  else {
    if(Communicator::world->me() < ij_pairs % Communicator::world->nproc()) {
      Ktilde = (double ***) malloc(((ij_pairs / Communicator::world->nproc()) + 1) * sizeof (double **));
      for(ij = 0; ij < ij_pairs; ij++) {
        if(Communicator::world->me() == ij_owner[ij])
          Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
      }
    }
    else {
      Ktilde = (double ***) malloc(ij_pairs / Communicator::world->nproc() * sizeof (double **));
      for(ij = 0; ij < ij_pairs; ij++) {
        if(Communicator::world->me() == ij_owner[ij])
          Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
      }
    }
  }

  v = 0;
  int a2, b2;
  /* This now loops over all ij_pairs even if we are neglecting
   * distant pairs */
  for(ij = 0; ij < ij_pairs; ij++, v++) {
    i = ij_map[ij][0];
    j = ij_map[ij][1];
    if(v % Communicator::world->nproc() == Communicator::world->me()) {
      for(k = 0, a = 0; k < natom; k++) {
        if(pairdomain[ij][k]) {
          for(t = aostart[k]; t <= aostop[k]; t++, a++) {
            for(l = 0, b = 0; l < natom; l++) {
              if(pairdomain[ij][l]) {
                for(u = aostart[l]; u <= aostop[l]; u++, b++) {
                  a2 = uniteddomain_abs2rel[i][t];
                  b2 = uniteddomain_abs2rel[j][u];
                  Ktilde[ij_local[ij]][a][b] = C_DDOT(ribasis->nbf(), &(B_ir_p[i][a2][0]), 1, &(B_ir_p[j][b2][0]), 1);
                }
              }
            }
          }
        }
      }
    }
  } // End of ij_pair loop

  free_block(J_mhalf);

  for (i = 0; i < nocc; i++)
      for (j = 0; j < uniteddomain_len[i]; j++)
          free(B_ir_p[i][j]);
  for (i = 0; i < nocc; i++)
      free(B_ir_p[i]);
  free(B_ir_p);

#ifdef TIME_DF_LMP2
if(Communicator::world->me() == 0) timer_off("Compute DF-LMP2");
#endif

}

}} // end namespaces

