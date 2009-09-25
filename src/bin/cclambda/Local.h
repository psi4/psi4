/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/

namespace psi { namespace cclambda {

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
  double ***V;
  double ***W;
  double *eps_occ;
  double **eps_vir;
  double cutoff;
  char *method;
  char *weakp;
  int filter_singles;
  double cphf_cutoff;
  char *freeze_core;
  char *pairdef;
};

}} // namespace psi::cclambda
