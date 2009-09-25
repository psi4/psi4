/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/

namespace psi { namespace cceom {

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
  int do_singles;
  double ***V;
  double ***W;
  double *eps_occ;
  double **eps_vir;
  double cutoff;
  char *method;
  char *weakp;
  char *precon;
  int filter_singles;
  double weak_pair_energy;
};

}} // namespace psi::cceom
