#ifndef CIM_H
#define CIM_H
namespace boost {
template<class T> class shared_ptr;
}
namespace psi{
class CIM{
  public:
    CIM(Options&options);
    ~CIM();
    void BuildClusters();
    double **Fock;
    int nirreps,nso,nmo,ndocc,nvirt,nfzc,nfzv,ndoccact;
    double thresh1,thresh2,thresh3;
    Options options_;

    /*
     * build occupied domains:
     * central orbitals | central mo domain | environmental domain
     */
    int*cdomainsizes,**cdomains,ndomains;
    int**central,*ncentral;
    int**env,*nenv;
    int**modomain,*nmodomain;
    void OccupiedDomains();
};
};
#endif
