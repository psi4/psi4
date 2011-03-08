/*
 *  ks.h
 *  matrix
 *
 *  Created by Rob Parrish on 3/7/2011  
 *
 */

#ifndef KS_H
#define KS_H

#include <psi4-dec.h>
#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include <libfunctional/superfunctional.h>

using namespace psi;
using namespace psi::functional;

namespace psi {

class Matrix;
class Integrator;
class Properties;

namespace scf{

class KS {

protected:
    // Superfunctional object
    shared_ptr<SuperFunctional> functional_;
    // Properties object
    shared_ptr<Properties> properties_;
    // Integrator object
    shared_ptr<Integrator> integrator_;
    // Values (E_xc, <\rho>, <\rho x>, etc)
    std::map<std::string, double> quad_values_;
    // primary basis set (might get fancy later)
    shared_ptr<BasisSet> basisset_;
    // Options object
    Options& options_; 
    // Molecule object  
    shared_ptr<Molecule> molecule_;
    // PSIO object
    shared_ptr<PSIO> psio_;

    // Compute E_xc and the V matrix
    virtual void form_V() = 0;
    // Build functional, grid, etc
    void common_init();
public:
    KS(Options & options, shared_ptr<PSIO> psio); 
    ~KS();
};

class RKS : public RHF, public KS {

protected:
    // Alpha/Beta spin Kohn-Sham Potential (identical)
    shared_ptr<Matrix> V_;
    // Compute E_xc and the V matrix
    virtual void form_V();
    virtual void form_G();
    virtual double compute_E();
    
    void common_init();
public:
    RKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    RKS(Options& options, shared_ptr<PSIO> psio);
    ~RKS();
};

class UKS : public UHF, public KS {

protected:
    // Alpha spin Kohn-Sham Potential
    shared_ptr<Matrix> Va_;
    // Beta spin Kohn-Sham Potential
    shared_ptr<Matrix> Vb_;
    // Compute E_xc and the V matrices
    virtual void form_V();
    virtual void form_G();
    virtual double compute_E();

    void common_init();
public:
    UKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    UKS(Options& options, shared_ptr<PSIO> psio);
    ~UKS();
};

}} // Namespaces

#endif

