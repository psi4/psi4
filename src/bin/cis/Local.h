/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/

namespace psi { namespace cis {

struct Local {
  int natom;
  int nso;
  int nocc;
  int nvir;
  int *aostart;
  int *aostop;
  int **domain;
  int **pairdomain;
  int *pairdom_len;
  int *pairdom_nrlen;
  int *weak_pairs;
  int ghost;
  double ***V;
  double ***W;
  double *eps_occ;
  double **eps_vir;
  double cutoff;
  char *method;
  char *weakp;
  char *precon;
  double weak_pair_energy;
  double **U;
  double **WW;
  double amp_print_cutoff;
};

}} // namespace psi::cis
