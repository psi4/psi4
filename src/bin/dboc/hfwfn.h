/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_bin_dboc_hfwfn_h_
#define _psi3_bin_dboc_hfwfn_h_

#include <libchkpt/chkpt.h>

namespace psi { namespace dboc {

class HFWavefunction {

  int num_mo_;
  int num_so_;
  int num_ao_;
  reftype refnum_;

  int ndocc_;
  int nsocc_;

  int nirreps_;
  int *clsdpi_;
  int *openpi_;
  int *orbspi_;

  double **alpha_evec_;
  double **beta_evec_;
  double **aotoso_;
  double **rref_;

  public:
  HFWavefunction();
  ~HFWavefunction();

  int num_mo();
  int num_so();
  int ndocc();
  int nsocc();
  int nalpha();
  int nbeta();
  int num_ao();
  int nirreps();
  int* clsdpi();
  int* openpi();
  int* orbspi();
  double** alpha_evec();
  double** beta_evec();
  double** aotoso();
  double** rref();
  void set_rref(double**);
};

}} /* namespace psi::dboc */

#endif

