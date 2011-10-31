#ifndef APPS_H 
#define APPS_H

#include <libmints/wavefunction.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class BasisSet;
class Matrix;
class TwoBodyAOInt;

// => BASE CLASSES <= //

class RBase : public Wavefunction {

protected:

    int print_;

    SharedMatrix C_;

    SharedMatrix Cfocc_;
    SharedMatrix Cfvir_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;

    boost::shared_ptr<Vector> eps_focc_;
    boost::shared_ptr<Vector> eps_fvir_;
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;

    SharedMatrix AO2USO_;

    double Eref_;

    void common_init();
    
public:

    RBase();
    virtual ~RBase();

    virtual bool restricted() const { return true; }

    void set_print(int print) { print_ = print; }

}; 

// => APPLIED CLASSES <= //

class RCIS : public RBase {

protected:

    std::vector<SharedMatrix > singlets_;  
    std::vector<SharedMatrix > triplets_;  
    std::vector<double> E_singlets_;  
    std::vector<double> E_triplets_;  

    virtual void print_header();
    virtual void print_wavefunctions();
    virtual void print_amplitudes();
    virtual void print_transitions();

    virtual SharedMatrix TDmo(SharedMatrix T1, bool singlet = true);
    virtual SharedMatrix TDso(SharedMatrix T1, bool singlet = true);
    virtual SharedMatrix TDao(SharedMatrix T1, bool singlet = true);

    virtual SharedMatrix Dmo(SharedMatrix T1, bool diff = false);
    virtual SharedMatrix Dso(SharedMatrix T1, bool diff = false);
    virtual SharedMatrix Dao(SharedMatrix T1, bool diff = false);

    virtual std::pair<SharedMatrix, boost::shared_ptr<Vector> > Nmo(SharedMatrix T1, bool diff = false);
    virtual std::pair<SharedMatrix, boost::shared_ptr<Vector> > Nso(SharedMatrix T1, bool diff = false);
    virtual std::pair<SharedMatrix, boost::shared_ptr<Vector> > Nao(SharedMatrix T1, bool diff = false);

public:
    RCIS();
    virtual ~RCIS();

    virtual double compute_energy();

};


}
#endif
