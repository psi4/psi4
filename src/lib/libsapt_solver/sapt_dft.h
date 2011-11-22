#ifndef SAPTDFT_H
#define SAPTDFT_H

#include "sapt0.h"

namespace psi { 

class Quadrature;

namespace sapt {

class MP2C : public SAPT0 {

private:

protected:

    /// The quadrature over i*omega
    boost::shared_ptr<Quadrature> quad_;
  
    /// The metric
    SharedMatrix J_; 
    /// The inverse of the metric
    SharedMatrix Jinv_; 
    /// The coupling matrix for monomer A
    SharedMatrix W_A_;
    /// The coupling matrix for monomer B
    SharedMatrix W_B_;
    /// The vector of X matrices for monomer A
    std::vector<SharedMatrix > X_A_;
    /// The vector of X matrices for monomer B
    std::vector<SharedMatrix > X_B_;
    
    /// Print the header
    virtual void print_header(); 
    /// Print the results
    virtual void print_results();
    
    /// Form the omega quadrature
    void form_quadrature();
    /// Form the J matrix inverse (Jinv_)
    void form_J();
    /// Form the W coupling matrices
    void form_W();
    /// Form the X0 response tensors
    void form_X0();
    /// Form the XC response tensors
    void form_XC();

    /// Do a Casimir-Polder Quadrature with whatever is in X 
    std::vector<double> casimirPolder();


public:
    MP2C(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    virtual ~MP2C();

    virtual double compute_energy();
};

}}

#endif
