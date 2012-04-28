#ifndef KS_H
#define KS_H
/*
 *  ks.h
 *  matrix
 *
 *  Created by Rob Parrish on 3/7/2011
 *
 */


#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include <libfunctional/superfunctional.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Matrix;
class Properties;
class TwoBodySOInt;
class ErfERI;
class VBase;

namespace scf{

class KS {

protected:
    /// KS Potential (the heart of the algorithm)
    boost::shared_ptr<VBase> potential_;
    /// Pointer to potential's functional
    boost::shared_ptr<SuperFunctional> functional_;
    /// primary basis set (might get fancy later)
    boost::shared_ptr<BasisSet> basisset_;
    /// primary so basis set
    boost::shared_ptr<SOBasisSet> sobasisset_;
    /// Options object
    Options& options_;
    /// Molecule object
    boost::shared_ptr<Molecule> molecule_;
    /// PSIO object
    boost::shared_ptr<PSIO> psio_;
    /// ERI object for omega integrals
    boost::shared_ptr<TwoBodySOInt> omega_eri_;
    /// Factory (for Spherical Harmonics)
    boost::shared_ptr<IntegralFactory> omega_factory_;

    /// Compute E_xc and the V matrix
    virtual void form_V() = 0;
    /// Build functional, grid, etc
    void common_init();
    /// Analog to HF::process_tei for omega K
    template <class OmegaKFunctor> void process_omega_tei(OmegaKFunctor & functor);

public:
    KS(Options & options, boost::shared_ptr<PSIO> psio);
    virtual ~KS();
};

class RKS : public RHF, public KS {

protected:
    /// Alpha/Beta spin Kohn-Sham Potential (identical)
    SharedMatrix V_;
    /// Omega K
    SharedMatrix wK_;
    /// Compute E_xc and the V matrix
    virtual void form_V();
    virtual void form_G();
    virtual double compute_E();
    virtual void stability_analysis();
    virtual void integrals();
    virtual void finalize();

    void common_init();
public:
    RKS(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    RKS(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~RKS();
};

class UKS : public UHF, public KS {

protected:
    /// Alpha spin Kohn-Sham Potential
    SharedMatrix Va_;
    /// Beta spin Kohn-Sham Potential
    SharedMatrix Vb_;
    /// Omega Ka
    SharedMatrix wKa_;
    /// Omega Kb
    SharedMatrix wKb_;
    /// Compute E_xc and the V matrices
    virtual void form_V();
    virtual void form_G();
    virtual double compute_E();
    virtual void stability_analysis();
    virtual void integrals();
    virtual void finalize();

    void common_init();
public:
    UKS(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    UKS(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~UKS();
};

}} // Namespaces

#endif

