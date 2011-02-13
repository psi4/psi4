#ifndef SCF_PSEUDO_H
#define SCF_PSEUDO_H
#include <psi4-dec.h>

using namespace psi;

namespace psi {

class BasisSet;
class Matrix;
class Options;
class FittingMetric;
class SchwarzSieve;
class TwoBodyAOInt;
class PseudopotentialInt;
class PSIO;

namespace scf {

class PseudospectralHF {
    
    protected:

        shared_ptr<PSIO> psio_;
        // Threadsafe (A|mn) integral objects
        std::vector<shared_ptr<TwoBodyAOInt> > eri_;
        // Threadsafe (m|V|n) integral objects
        std::vector<shared_ptr<PseudopotentialInt> > pot_;
        // Inverse of fitting metric (Coulomb side)
        shared_ptr<FittingMetric> Jinv_;
        // X pseudospectral matrix
        shared_ptr<Matrix> X_;
        // Shared pointer to density matrix
        shared_ptr<Matrix> D_;
        // Shared pointer to Coulomb matrix
        shared_ptr<Matrix> J_;
        // Shared pointer to Exchange matrix
        shared_ptr<Matrix> K_;
        // Constant reference to the options object
        Options& options_;
        // Primary basis set
        shared_ptr<BasisSet> primary_;
        // Auxiliary basis set (coulomb)
        shared_ptr<BasisSet> auxiliary_;
        // Quadrature Points (P x 3)
        shared_ptr<Matrix> points_;
        // Schwarz Sieve object 
        shared_ptr<SchwarzSieve> schwarz_;

        // memory in doubles
        unsigned long int memory_;
        // number of pseudospectral points
        unsigned long int P_;
    
        // Helper methods
        void form_J_DF_RHF();
        void form_K_PS_RHF();

    public:
        // Constructor for RHF/RKS
        PseudospectralHF(shared_ptr<BasisSet> basis, shared_ptr<Matrix> D,
        shared_ptr<Matrix> J, shared_ptr<Matrix> K, shared_ptr<PSIO> psio, Options& opt);
        // Destructor
        ~PseudospectralHF();
        // RAII, builds grids, sieves, fitting metric
        void common_init();
        
        // form J and K
        void form_G_RHF(); 

};

}}

#endif
