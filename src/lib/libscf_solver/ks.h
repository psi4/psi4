/*
 *  ks.h
 *  matrix
 *
 *  Created by Rob Parrish on 3/7/2011
 *
 */

#ifndef KS_H
#define KS_H

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include <libfunctional/superfunctional.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Matrix;
class Integrator;
class Properties;

namespace scf{

class KS {

protected:
    /// Superfunctional object
    boost::shared_ptr<functional::SuperFunctional> functional_;
    /// Properties object
    boost::shared_ptr<Properties> properties_;
    /// Integrator object
    boost::shared_ptr<Integrator> integrator_;
    /// Values (E_xc, <\rho>, <\rho x>, etc)
    std::map<std::string, double> quad_values_;
    /// primary basis set (might get fancy later)
    boost::shared_ptr<BasisSet> basisset_;
    /// Options object
    Options& options_;
    /// Molecule object
    boost::shared_ptr<Molecule> molecule_;
    /// PSIO object
    boost::shared_ptr<PSIO> psio_;

    /// Compute E_xc and the V matrix
    virtual void form_V() = 0;
    /// Build functional, grid, etc
    void common_init();
public:
    KS(Options & options, boost::shared_ptr<PSIO> psio);
    ~KS();
};

class RKS : public RHF, public KS {

protected:
    /// Alpha/Beta spin Kohn-Sham Potential (identical)
    boost::shared_ptr<Matrix> V_;
    /// Compute E_xc and the V matrix
    virtual void form_V();
    virtual void form_G();
    virtual double compute_E();

    void common_init();
public:
    RKS(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    RKS(Options& options, boost::shared_ptr<PSIO> psio);
    ~RKS();
};

class UKS : public UHF, public KS {

protected:
    /// Alpha spin Kohn-Sham Potential
    boost::shared_ptr<Matrix> Va_;
    /// Beta spin Kohn-Sham Potential
    boost::shared_ptr<Matrix> Vb_;
    /// Compute E_xc and the V matrices
    virtual void form_V();
    virtual void form_G();
    virtual double compute_E();

    void common_init();
public:
    UKS(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    UKS(Options& options, boost::shared_ptr<PSIO> psio);
    ~UKS();
};

}} // Namespaces

#endif

