/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_ccenergy_local_h
#define _psi_src_bin_ccenergy_local_h

namespace psi { namespace ccenergy {

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
  std::string method;
  std::string weakp;
  int filter_singles;
  double weak_pair_energy;
  double cphf_cutoff;
  int freeze_core;
  std::string pairdef;
};

}} // namespace psi::ccenergy

#endif // _psi_src_bin_ccenergy_local_h
