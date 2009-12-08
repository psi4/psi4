/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/
#ifndef _psi_src_bin_lmp2_local_h_
#define _psi_src_bin_lmp2_local_h_
struct local {
  int **domain;
  int *domain_len;
  int **pairdomain;
  int *pairdom_len;
  int *aostart;
  int *aostop;
  double ***W;
  int *pairdom_nrlen;
  double **Rt_full;
  double ****T;
  double **evals;
  int **MR_shell;
  int num_unique_shells;
  double ***eri_1;
  double *B1;

  double **C;            /* SCF MO coefficients in the LO basis */
  double **F;            /* Fock matrix in the LO basis */
  double ***Ktilde;      /* ERI transformed to (ai|bj) where ij are occupied localized orbitals and ab are projected ao's */
  double Emp2;           /* LMP2 energy */
  double Emp2_old;       /* LMP2 energy from previous iteration */
  int iter;              /* iteration number */
  int conv;
  double DEmp2;
  double Drms;
  int ij_pairs;
};

struct trans_thread_data {
  int P, Q;
  int nump, numq;
  int my_t;
};

struct comm_thread_data {
  int t;
  int root;
  int div;
};

#endif /* Header guard */

