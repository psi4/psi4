/*! \file
    \ingroup CLAG
    \brief Compute orbital Z vector
*/


#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "indpairs.h"

extern "C" {
  extern FILE *infile;
  extern FILE *outfile;
}

namespace psi { namespace clag {

#define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)

extern int *ioff;
extern int print_lvl;

/*
** compute_zorb()
**
** Compute the orbital response z vector and the X_tilde matrix
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
void compute_zorb (double *tei,
  double **lag, double *epsilon, IndepPairs &IndPairs, 
  int nmo, int nocc, int nvir, double *Zvec, double **X_tilde)
{

  // form X_tilde_pq = X_{pq} p \leq q, or X_{qp} p>q (eq 41)
  for (int p=0; p<nmo; p++) {
    for (int q=0; q<nmo; q++) {
      if (p <= q) X_tilde[p][q] = lag[p][q];
      else X_tilde[p][q] = lag[q][p];
    }
  }

  // form Delta_X_{pq} = X_{pq} - X_{qp} where X is the Lagrangian (eq 42)
  double **Delta_X = block_matrix(nmo, nmo);
  for (int p=0; p<nmo; p++) {
    for (int q=0; q<nmo; q++) {
      Delta_X[p][q] = lag[p][q] - lag[q][p];
    }
  }

  fprintf(outfile, "Delta_X matrix:\n");
  print_mat(Delta_X,nmo,nmo,outfile);

  // form Delta_X_tilde[ai] = Delta_X[a][i] + Delta_X_prime[ai] +
  //   Delta_X_prime_prime[ai] (eq 60), where
  // Delta_X_prime[ai] =  \sum_{j>k,CI-IP} \Delta_X[j][k] * A_{jk,ai} /
  //   (\epsilon_k - \epsilon_j) (eq 55), and
  // Delta_X_prime_prime[ai] = \sum_{b>c,CI-IP} \Delta X[b][c] * A_{bc,ai} / 
  //   (\epsilon_c - \epsilon_b) (eq 56)

  double *Delta_X_tilde = init_array(nocc*nvir);
  for (int a=nocc,idx=0; a<nmo; a++) {
    for (int i=0; i<nocc; i++,idx++) {
      int pq, pqai, pa, qi, paqi, pi, qa, piqa;
      int ai = ioff[a] + i;
      Delta_X_tilde[idx] = Delta_X[a][i];
      // loop over CI independent pairs
      int *p_arr = IndPairs.get_p_ptr();
      int *q_arr = IndPairs.get_q_ptr();
      for (int pair=0; pair<IndPairs.get_num_pairs(); pair++) {
        int p = p_arr[pair];  int q = q_arr[pair];
        if ((p < nocc && q < nocc) ||
            (p >= nocc && q >= nocc)) {
          pq = INDEX(p,q);
          pqai = INDEX(pq,ai);
          pa = INDEX(p,a);
          qi = INDEX(q,i);
          paqi = INDEX(pa,qi);
          pi = INDEX(p,i);
          qa = INDEX(q,a);
          piqa = INDEX(pi,qa);
          // diagonal term never enters 
          // CDS: CHANGED THE SIGN BELOW FOR FUN...
          Delta_X_tilde[idx] += Delta_X[p][q] * (4.0 * tei[pqai] -
            tei[paqi] - tei[piqa]) / (epsilon[q] - epsilon[p]);
        }
      }
      fprintf(outfile, "Delta_X_tilde[%d][%d] = %10.6lf\n", Delta_X_tilde[idx],
        a, i);
    }
  }
 
  double **A = block_matrix(nocc*nvir,nocc*nvir);
  // form A^T_{bj,ai} = A_{bj,ai} = A_{ai,bj} I believe (should be symmetric)
  for (int a=nocc,row=0; a<nmo; a++) {
    for (int i=0; i<nocc; i++, row++) {
      for (int b=nocc,col=0; b<nmo; b++) {
        for (int j=0; j<nocc; j++,col++) {
          int bj = ioff[b] + j;
          int ai = ioff[a] + i;
          int aibj = INDEX(ai,bj);
          int ab = INDEX(a,b);
          int ij = INDEX(i,j);
          int abij = ioff[ab]+ij;
          int aj = ioff[a]+j;
          int ib = ioff[b]+i;
          int ajib = INDEX(aj,ib); 
          if (b==a && j==i) A[row][col] = epsilon[i] - epsilon[a];
          A[row][col] -= (4.0 * tei[aibj] - tei[abij] - tei[ajib]);
        }
      }
    }
  } 

  fprintf(outfile, "A matrix:\n");
  print_mat(A,nocc*nvir,nocc*nvir,outfile);

  // Solve \sum_{ai} A^T_{bj,ai} Z_{ai} = \Delta_X_tilde_{bj}
  double det;
  double *Zai = init_array(nocc*nvir);
  int *tmpi = init_int_array(nocc*nvir);
  C_DCOPY(nocc*nvir,Delta_X_tilde,1,Zai,1);
  // flin(A,Zai,nocc*nvir,1, &det);
  if (C_DGESV(nocc*nvir,1,A[0],nocc*nvir,tmpi,Zai,nocc*nvir) != 0) {
    fprintf(outfile, "compute_zorb: C_DGESV returned an error\n");
    exit(1);
  }
  free_block(A);
  free(tmpi);

  // now construct Z vector.  First the SCF-IP parts
  for (int a=nocc,idx=0; a<nmo; a++) {
    for (int i=0; i<nocc; i++,idx++) {
      int ai = ioff[a] + i;
      Zvec[ai] = Zai[idx];
      // fprintf(outfile, "Zvec[%d][%d] = %10.6lf\n", a, i, Zai[idx]);
    }
  }   

  // get CI IP's which aren't SCF IP's
  int *p_arr = IndPairs.get_p_ptr();
  int *q_arr = IndPairs.get_q_ptr();
  for (int pair=0; pair<IndPairs.get_num_pairs(); pair++) {
    int p = p_arr[pair];  int q = q_arr[pair];
    if ((p < nocc && q < nocc) ||
        (p >= nocc && q >= nocc)) {
        if (p<q) printf("Error: p<q in compute_zorb!\n");
        int pq = ioff[p]+q;
        Zvec[pq] = Delta_X[p][q] / (epsilon[q] - epsilon[p]);
     }
   }
        
  // we got Z!
  free(Zai);
  free(Delta_X_tilde);
  free_block(Delta_X);

}

}} // end namespace psi::clag

