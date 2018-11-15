#ifndef LIBPSIPE_PE_CALC_HANDLER_H
#define LIBPSIPE_PE_CALC_HANDLER_H

#include <string>

#include <armadillo>

#include <cppe/libcppe.hh>

using libcppe::CPPE;

namespace psi {

class PeCalcHandler {
protected:
  CPPE cppe_interface; //!< CPPE instance for calls to the CPPE library
  arma::mat P_gs; //!< copy of the most recent density matrix
  int m_nbas = 0; //!< number of basis functions
  libcppe::PeOptions m_options;

public:
  PeCalcHandler(libcppe::PeOptions options) : m_options(options) {};
  virtual ~PeCalcHandler () = default;
  virtual bool initialize() = 0;
  virtual void fock_contribution(arma::mat& Ptot, arma::mat& out, double* energy) = 0;
  virtual void print_energy() = 0;
  // virtual void perturbative_exc_energy_correction(double* p_exc, double* energy) = 0;
};

}
#endif // LIBPE_PE_CALC_HANDLER_H
