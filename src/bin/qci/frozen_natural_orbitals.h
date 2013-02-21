#ifndef FROZENNO_H
#define FROZENNO_H
#include"psi4-dec.h"
#include<libmints/wavefunction.h>

namespace boost {
template<class T> class shared_ptr;
}
namespace psi{
  class Wavefunction;
}

namespace psi{namespace qci{
class FrozenNO : public Wavefunction{
  public:
    FrozenNO(Options&options);
    ~FrozenNO();

    double compute_energy();
    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }

  protected:

    Options & options_;

    // mp2 energy in full basis
    double emp2;
    long int nirreps,nso,nmo,ndocc,nvirt,nfzc,nfzv,ndoccact,nvirt_no;

    void NaturalOrbitals();
    void TransformIntegrals(double*Dab);
    void TransformOVOV();
};
}}

#endif
