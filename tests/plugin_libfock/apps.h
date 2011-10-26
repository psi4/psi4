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

    boost::shared_ptr<Matrix> C_;

    boost::shared_ptr<Matrix> Cfocc_;
    boost::shared_ptr<Matrix> Cfvir_;
    boost::shared_ptr<Matrix> Caocc_;
    boost::shared_ptr<Matrix> Cavir_;

    boost::shared_ptr<Vector> eps_focc_;
    boost::shared_ptr<Vector> eps_fvir_;
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;

    boost::shared_ptr<Matrix> AO2USO_;

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

    std::vector<boost::shared_ptr<Matrix> > singlets_;  
    std::vector<boost::shared_ptr<Matrix> > triplets_;  
    std::vector<double> E_singlets_;  
    std::vector<double> E_triplets_;  

    virtual void print_header();
    virtual void print_wavefunctions();
    virtual void print_amplitudes();
    virtual void print_transitions();

    virtual boost::shared_ptr<Matrix> TDmo(boost::shared_ptr<Matrix> T1, bool singlet = true);
    virtual boost::shared_ptr<Matrix> TDso(boost::shared_ptr<Matrix> T1, bool singlet = true);
    virtual boost::shared_ptr<Matrix> TDao(boost::shared_ptr<Matrix> T1, bool singlet = true);

    virtual boost::shared_ptr<Matrix> Dmo(boost::shared_ptr<Matrix> T1, bool diff = false);
    virtual boost::shared_ptr<Matrix> Dso(boost::shared_ptr<Matrix> T1, bool diff = false);
    virtual boost::shared_ptr<Matrix> Dao(boost::shared_ptr<Matrix> T1, bool diff = false);

    virtual std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > Nmo(boost::shared_ptr<Matrix> T1, bool diff = false);
    virtual std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > Nso(boost::shared_ptr<Matrix> T1, bool diff = false);
    virtual std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > Nao(boost::shared_ptr<Matrix> T1, bool diff = false);

public:
    RCIS();
    virtual ~RCIS();

    virtual double compute_energy();

};


}
#endif
