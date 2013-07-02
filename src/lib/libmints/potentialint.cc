#include "mints.h"
#include "potentialint.h"

namespace psi{

PCMPotentialInt::PCMPotentialInt(std::vector<SphericalTransform>& trans,
boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv):
    PotentialInt(trans, bs1, bs1)
{
    // We don't want to transform the integrals from Cartesian (6d, 10f, ...) to Pure (5d, 7f, ...)
    // for each external charge.  It'll be better to backtransform the density / Fock matrices to 
    // the Cartesian basis, if necessary.
    force_cartesian_ = true;
}

} //Namespace
