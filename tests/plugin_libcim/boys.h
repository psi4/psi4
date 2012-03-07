#ifndef BOYS_H
#define BOYS_H
namespace boost {
template<class T> class shared_ptr;
}
namespace psi{
class Boys{
  public:
    Boys(Options&options);
    ~Boys();
    void Localize();
    double LocalizationSum();
    void ShuffleOrbitals();
    void GetAngle(int i,int j,double&angle);
    void Rotate(int i,int j,double angle);
    void SoToMo(int nso,int nmo,double**mat,double**trans,double**temp);
    int nirreps,nso,nmo,ndocc,nvirt,nfzc,nfzv,ndoccact,*reorder;
    double ***mu,**Clmo_pointer,conv;
    SharedMatrix Clmo;
    SharedMatrix Fock;
    Options options_;
};
};
#endif
