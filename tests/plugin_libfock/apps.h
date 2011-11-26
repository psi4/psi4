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

    std::vector<boost::tuple<double, int, int, int> > states_;
    std::vector<SharedMatrix > singlets_;  
    std::vector<SharedMatrix > triplets_;  
    std::vector<double> E_singlets_;  
    std::vector<double> E_triplets_;  

    void sort_states();    

    virtual void print_header();
    virtual void print_wavefunctions();
    virtual void print_amplitudes();
    virtual void print_transitions();
    virtual void print_densities();

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

class RTDHF : public RBase {

protected:

    std::vector<SharedMatrix > singlets_X_;  
    std::vector<SharedMatrix > triplets_X_;  
    std::vector<SharedMatrix > singlets_Y_;  
    std::vector<SharedMatrix > triplets_Y_;  
    std::vector<double> E_singlets_;  
    std::vector<double> E_triplets_;  

    virtual void print_header();

public:
    RTDHF();
    virtual ~RTDHF();

    virtual double compute_energy();

};


class RCPHF : public RBase {

protected:

    // OV-Rotations
    std::vector<SharedMatrix> x_;  
    // OV-Perturbations
    std::vector<SharedMatrix> b_;  

    virtual void print_header();

    void add_named_tasks();
    void analyze_named_tasks();

    void add_polarizability();
    void analyze_polarizability();

    std::set<std::string> tasks_;

public:
    RCPHF();
    virtual ~RCPHF();

    /// Solve for all perturbations currently in b 
    virtual double compute_energy();

    /// Perturbation vector queue, shove tasks onto this guy before compute_energy
    std::vector<SharedMatrix>& b() { return b_; }
    /// Resultant solution vectors, available after compute_energy is called
    std::vector<SharedMatrix>& x() { return x_; }

    /// Add a named task
    void add_task(const std::string& task);
        
};

}
#endif
