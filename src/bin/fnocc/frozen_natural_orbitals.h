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

namespace psi{namespace fnocc{

// base class
class FrozenNO : public Wavefunction {
  public:
    FrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options);
    ~FrozenNO();

    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
    void ComputeNaturalOrbitals();

  protected:

    // mp2 energy in full basis
    double emp2;
    long int nso,nmo,ndocc,nvirt,nfzc,nfzv,ndoccact,nvirt_no;

    void common_init();
};

class DFFrozenNO : public FrozenNO {
  public:
    DFFrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options);
    ~DFFrozenNO();

    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
    void ComputeNaturalOrbitals();
    void ThreeIndexIntegrals();

  protected:

    void ModifyCa(double*Dab);
    void ModifyCa_occ(double*Dij);
    void BuildFock(long int nQ,double*Qso,double*F);
    void TransformQ(long int nQ,double*Qso);

};

}}

#endif
