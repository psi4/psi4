#ifndef LIBPE_CPPE_PE_CALC_HANDLER_H
#define LIBPE_CPPE_PE_CALC_HANDLER_H

#include "pe_calc_handler_i.hh"
#include <cppe/core/cppe_state.hh>

namespace psi {

class CppePeCalcHandler : public PeCalcHandler {
private:
  // TODO: this will go into PeStateHandler/PeState
  int iteration;
  std::vector<libcppe::Potential> potentials;
  libcppe::CppeState cppe_state;
  int dm_count;

public:
  ~CppePeCalcHandler() {};
  CppePeCalcHandler(libcppe::PeOptions options);
  bool initialize();
  void fock_contribution(arma::mat& Ptot,
    arma::mat& out, double* energy);
  void print_energy();
  // void perturbative_exc_energy_correction(double* p_exc, double* energy);
};

}
#endif // LIBPE_CPPE_PE_CALC_HANDLER_H
