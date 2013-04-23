/*! \file
    \ingroup CLAG
    \brief Compute orbital relaxation contributions to gradient
*/


#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include "indpairs.h"

extern "C" {
  extern FILE *infile;
  extern FILE *outfile;
}

namespace psi { namespace clag {

#define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)
#define INDEX2(i,j,n) ((i)*(n) + (j))



extern int *ioff;
extern int print_lvl;

/*
** relax_pdms()
**
** Add in the orbital relaxation from the orbital z vector into the
** onepdm, twopdm, and effective Lagrangian
**
** Note: using equations from "Analytic Gradients of Configuration
** Interaction Energies", C. D. Sherrill, 2008
**
** C. David Sherrill
** Georgia Institute of Technology
** July 2009
**
** NOTE: I have assumed that all occ/vir rotations are listed among the
** CI independent pair list.  This is not guaranteed if we run a 
** calculation with both occ and vir inside RAS II.  Need to go back
** and put in a flag for IndepPairs to force all occ/vir to be included
** (it's not problematic to do so because the terms involving those at
** the CI level will cancel out anyway).  See CI gradient notes p. 10.
*/
void relax_pdms(double **onepdm, double *tpdm, double *tei,
  double **X_tilde, double *epsilon, IndepPairs &IndPairs, 
  int nmo, int nocc, int npop, double *Zvec, int opdm_file, 
  int lag_file, int tpdm_file)
{
  int *p_arr = IndPairs.get_p_ptr();
  int *q_arr = IndPairs.get_q_ptr();

  for (int pair=0; pair<IndPairs.get_num_pairs(); pair++) {
    int p = p_arr[pair];  int q = q_arr[pair];
    if (p<q) printf("relax_pdms: found p<q!\n"); // debug test
    int pq = INDEX(p,q);

    // oops, I think the PDM's maybe need to be bigger for gradients,
    // not just size of populated orbs...
    if (p >= npop || q >= npop) {
      fprintf(outfile, "relax_pdms: ERROR: oops, need bigger pdms!\n");
      continue;
    }
    onepdm[p][q] += 2.0*Zvec[pq];

    // halve the 2PDM contributions because the backtransform assumes
    // bra-ket symmetry, which effectively doubles whatever we add/subtract
    // here
    for (int k=0; k<nocc; k++) {
      int kk = k*npop + k;
      int pqkk = INDEX(INDEX2(p,q,npop),kk);
      tpdm[pqkk] += 2.0 * Zvec[pq];  // would have been 4.0

      int pk = INDEX2(p,k,npop);
      int qk = INDEX2(q,k,npop);
      int pkqk = INDEX(pk,qk);
      tpdm[pkqk] -= 1.0 * Zvec[pq];  // would have been 2.0
    }

    X_tilde[p][q] += 2.0 * Zvec[pq] * epsilon[q];
    for (int i=0; i<nocc; i++) {
      for (int j=0; j<nocc; j++) { 
        int pqij = INDEX(pq,INDEX(i,j));
        int piqj = INDEX(INDEX(p,i),INDEX(q,j));
        X_tilde[i][j] += 2.0 * Zvec[pq] * (2.0 * tei[pqij] - tei[piqj]);
      }
    }
        
  }

  // write out modified matrices
  psio_open(opdm_file, PSIO_OPEN_OLD);
  psio_write_entry(opdm_file, "MO-basis OPDM", (char *) onepdm[0],
    npop * npop * sizeof(double));
  psio_close(opdm_file,1); 

  if (print_lvl > 1) {
    fprintf(outfile,"\nRelaxed one-particle density matrix\n\n");
    print_mat(onepdm, nmo, nmo, outfile);
  }

  psio_open(lag_file, PSIO_OPEN_OLD);
  psio_write_entry(lag_file, "MO-basis Lagrangian", (char *) X_tilde[0],
    nmo*nmo*sizeof(double));
  psio_close(lag_file,1);

  if (print_lvl > 1) {
    fprintf(outfile,"\nEffective Lagrangian Matrix\n\n");
    print_mat(X_tilde, nmo, nmo, outfile);
  }

  if (print_lvl > 1) {
    fprintf(outfile,"\nRelaxed two-particle density matrix\n\n");
    print_array(tpdm, nmo*nmo, outfile);
  }

  struct iwlbuf TBuff;
  iwl_buf_init(&TBuff, tpdm_file, 0.0, 0, 0);

  for (int i=0; i<npop; i++) {
    for (int j=0; j<npop; j++) {
      for (int k=0; k<=i; k++) {
        int lmax = (k==i) ? j+1 : npop;
        for (int l=0; l<lmax; l++) {
          int ij = i * npop + j;
          int kl = k * npop + l;
          int ijkl = ioff[ij] + kl;
          iwl_buf_wrt_val(&TBuff,i,j,k,l,tpdm[ijkl],0,outfile,0);
        }
      }
    }
  }

  iwl_buf_flush(&TBuff, 1);
  iwl_buf_close(&TBuff, 1);
 
}

}} // end namespace psi::clag

