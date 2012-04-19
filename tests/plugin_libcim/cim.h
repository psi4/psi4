#ifndef CIM_H
#define CIM_H
#include"psi4-dec.h"
#include"boys.h"
#include<libmints/matrix.h>
#include<libmints/wavefunction.h>

#define isempty -999

namespace boost {
template<class T> class shared_ptr;
}
namespace psi{
  class Wavefunction;
}

namespace psi{
class CIM : public Wavefunction{
  public:
    CIM();
    ~CIM();
    void BuildClusters();
    double *epsSave;
    double **Fock;
    int nirreps,nso,nmo,ndocc,nvirt,nfzc,nfzv,ndoccact,maxndomains;
    // local variables:
    int nvirt_,nfzc_,ndoccact_,nfzv_;
    double thresh1,thresh2,thresh3;
    Options options_;

    // boys localization
    boost::shared_ptr<Boys> boys;

    /*
     * apparently if i want to inherit wavefunction, i need the following:
     */
    double compute_energy();
    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }

    /*
     * build occupied domains:
     * central orbitals | central mo domain | environmental domain
     */
    int*skip;
    int*domainsize,**domain,ndomains;
    int**central,*ncentral;
    int**env,*nenv;
    int**modomain,*nmodomain;
    void OccupiedDomains();
    void SECIM();

    /*
     * cluster lmo->mo transformation matrix
     */
    SharedMatrix localClmo;
    /*
     * cluster Fock matrix
     */
    SharedMatrix localFock;

    /*
     * build virtual spaces
     */
    void VirtualSpaces(int cluster,int clusternum);
    void MP2Density(double*Dab,int cluster);
    void OppositeSpinMP2Density(double*Dab,int cluster);

    /*
     * transform integrals
     */
    void TransformIntegrals(int cluster,SharedMatrix Clmo,int clusternum);

    /*
     * customize wavefunction for cluster, {i}
     */
    void ClusterWavefunction(int cluster);

    /*
     * build quasicanonical orbitals for cluster, {i}
     */
    void QuasiCanonicalOrbitals(int cluster);

    /*
     * df integrals to approximate mp2 opdm
     */
    double**Qov;
    int nQ;
};
};
#endif
