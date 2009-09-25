/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/

namespace psi { namespace ccresponse {

struct Local {
  int natom;
  int nso;
  int nocc;
  int nvir;
  int cart;
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
  int filter_singles;
  char *method;
  char *weakp;
  double cphf_cutoff;
  char *freeze_core;
  char *pairdef;
};

}} // namespace psi::ccresponse
