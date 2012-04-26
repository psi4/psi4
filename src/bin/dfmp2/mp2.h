#ifndef TOTAL_DFMP2_H
#define TOTAL_DFMP2_H

#include <libmints/wavefunction.h>
#include <map>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class PSIO;
class Chkpt;

namespace dfmp2 {

class DFMP2 : public Wavefunction {

protected:

    // Auxiliary basis
    boost::shared_ptr<BasisSet> ribasis_;
    // Energy map
    std::map<std::string, double> energies_;
    // Gradients map
    std::map<std::string, SharedMatrix> gradients_;

    // Same-spin scale
    double sss_;
    // Opposite-spin scale
    double oss_;

    void common_init();
    // Common printing of energies/SCS
    virtual void print_energies();
    virtual void print_gradients();

    // Print header/reference information
    virtual void print_header() = 0;
    // Form the (A|ia) = (A|mn) C_mi C_na tensor(s)
    virtual void form_Aia() = 0;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia)
    virtual void form_Qia() = 0;
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia) and J_QA^-1 (A|ia)
    virtual void form_Qia_grad() = 0;
    // Form the energy contributions
    virtual void form_energy() = 0;
    // Form the energy contributions and gradients
    virtual void form_gradient() = 0;
    // Form the small gamma
    virtual void form_gamma() = 0;
    // Form the (A|B)^x contribution to the gradient
    virtual void form_AB_x_terms() = 0;
    // Form the (A|mn)^x contribution to the gradient
    virtual void form_Amn_x_terms() = 0;

    // Compute singles correction [for ROHF-MBPT(2) or dual-basis]
    virtual void form_singles();
    // Apply the fitting and transposition to a given disk entry Aia tensor
    virtual void apply_fitting(SharedMatrix Jm12, unsigned int file, unsigned long int naux, unsigned long int nia);
    // Apply the fitting again to a given disk entry Qia tensor
    virtual void apply_fitting_grad(SharedMatrix Jm12, unsigned int file, unsigned long int naux, unsigned long int nia);
    // Form the inverse square root of the fitting metric, or read it off disk
    virtual SharedMatrix form_inverse_metric();
    // Form an abstract gamma
    virtual void apply_gamma(unsigned int file, unsigned long int naux, unsigned long int nia);

public:
    DFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    virtual ~DFMP2();

    double compute_energy();
    SharedMatrix compute_gradient();

}; 

class RDFMP2 : public DFMP2 {

protected:

    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    
    SharedVector eps_aocc_;
    SharedVector eps_avir_;

    void common_init();

    // Print additional header
    virtual void print_header();
    // Form the (A|ia) = (A|mn) C_mi C_na tensor(s)
    virtual void form_Aia();
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia)
    virtual void form_Qia();
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia) and J_QA^-1 (A|ia)
    virtual void form_Qia_grad();
    // Form the energy contributions
    virtual void form_energy();
    // Form the energy contributions and gradients
    virtual void form_gradient();
    // Form the small gamma
    virtual void form_gamma();
    // Form the (A|B)^x contribution to the gradient
    virtual void form_AB_x_terms();
    // Form the (A|mn)^x contribution to the gradient
    virtual void form_Amn_x_terms();

public:
    RDFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    virtual ~RDFMP2();

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
};

class UDFMP2 : public DFMP2 {

protected:

    SharedMatrix Caocc_a_;
    SharedMatrix Cavir_a_;
    SharedMatrix Caocc_b_;
    SharedMatrix Cavir_b_;
    
    SharedVector eps_aocc_a_;
    SharedVector eps_avir_a_;
    SharedVector eps_aocc_b_;
    SharedVector eps_avir_b_;

    void common_init();

    // Print additional header
    virtual void print_header();
    // Form the (A|ia) = (A|mn) C_mi C_na tensor(s)
    virtual void form_Aia();
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia)
    virtual void form_Qia();
    // Apply the fitting (Q|ia) = J_QA^-1/2 (A|ia) and J_QA^-1 (A|ia)
    virtual void form_Qia_grad();
    // Form the energy contributions
    virtual void form_energy();
    // Form the energy contributions and gradients
    virtual void form_gradient();
    // Form the small gamma
    virtual void form_gamma();
    // Form the (A|B)^x contribution to the gradient
    virtual void form_AB_x_terms();
    // Form the (A|mn)^x contribution to the gradient
    virtual void form_Amn_x_terms();

public:
    UDFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    virtual ~UDFMP2();

    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }

};

class RODFMP2 : public UDFMP2 {

protected:

    void common_init();

    // Print additional header
    virtual void print_header();

public:
    RODFMP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    virtual ~RODFMP2();

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return false; }
};

}}

#endif
