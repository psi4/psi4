#include <liboptions/liboptions.hpp>

#include "mrccsd_t.h"

namespace psi{ namespace psimrcc{

MRCCSD_T::MRCCSD_T(Hamiltonian* h_eff_) : h_eff(h_eff_)
{
  startup();
  check_intruders();
  if(triples_algorithm == SpinAdaptedTriples)
    compute_spin_adapted();
  else if(triples_algorithm == RestrictedTriples)
    compute_restricted();
  else
    compute();
}

MRCCSD_T::~MRCCSD_T()
{
  cleanup();
}

}} /* End Namespaces */
