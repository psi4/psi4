#ifndef LIBSCF_SAD_H
#define LIBSCF_SAD_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class BasisSet;
class Molecule;
class Matrix;

namespace scf {

class SADGuess {

protected:

    int print_;
    int debug_;

    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
    SharedMatrix AO2SO_;

    int nalpha_;
    int nbeta_;

    Options& options_;

    SharedMatrix Da_;
    SharedMatrix Db_;
    SharedMatrix Ca_;
    SharedMatrix Cb_;

    void common_init();

    SharedMatrix form_D_AO();
    void getUHFAtomicDensity(boost::shared_ptr<BasisSet> atomic_basis, int n_electrons, int multiplicity, double** D);
    void atomicUHFHelperFormCandD(int nelec, int norbs,double** Shalf, double**F, double** C, double** D);

    void form_D();
    void form_C();

public:

    SADGuess(boost::shared_ptr<BasisSet> basis, int nalpha, int nbeta, Options& options);
    virtual ~SADGuess();

    void compute_guess();
    
    SharedMatrix Da() const { return Da_; } 
    SharedMatrix Db() const { return Db_; } 
    SharedMatrix Ca() const { return Ca_; } 
    SharedMatrix Cb() const { return Cb_; } 

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

}; 

}}

#endif
