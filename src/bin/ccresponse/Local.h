/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <string>

namespace psi { namespace ccresponse {

struct Local {
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
  std::string method;
  std::string weakp;
  double cphf_cutoff;
  std::string freeze_core;
  std::string pairdef;
};

}} // namespace psi::ccresponse
