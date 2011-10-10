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

    boost::shared_ptr<Matrix> Cfocc_;
    boost::shared_ptr<Matrix> Cfvir_;
    boost::shared_ptr<Matrix> Caocc_;
    boost::shared_ptr<Matrix> Cavir_;

    boost::shared_ptr<Vector> eps_focc_;
    boost::shared_ptr<Vector> eps_fvir_;
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;

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

    std::vector<boost::shared_ptr<Vector> > singlets_;  
    std::vector<boost::shared_ptr<Vector> > triplets_;  
    std::vector<double> E_singlets_;  
    std::vector<double> E_triplets_;  

    virtual void print_header();

public:
    RCIS();
    virtual ~RCIS();

    virtual double compute_energy();

};


}
#endif
