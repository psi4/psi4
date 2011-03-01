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
class PseudospectralInt;
class PSIO;

namespace scf {

class PseudospectralHF {
    
    protected:
        // Whether the alpha and beta orbitals are equal or not
        bool restricted_;
        // The psio object
        shared_ptr<PSIO> psio_;
        // Threadsafe (A|mn) integral objects
        std::vector<shared_ptr<TwoBodyAOInt> > eri_;
        // Threadsafe (m|V|n) integral objects
        std::vector<shared_ptr<PseudospectralInt> > pot_;
        // Inverse of fitting metric (Coulomb side)
        shared_ptr<FittingMetric> Jinv_;
        // X pseudospectral matrix
        shared_ptr<Matrix> X_;
        // Shared pointer to alpha density matrix
        shared_ptr<Matrix> Da_;
        // Shared pointer to beta density matrix
        shared_ptr<Matrix> Db_;
        // Shared pointer to alpha Coulomb matrix
        shared_ptr<Matrix> Ja_;
        // Shared pointer to beta Coulomb matrix
        shared_ptr<Matrix> Jb_;
        // Shared pointer to alpha Exchange matrix
        shared_ptr<Matrix> Ka_;
        // Shared pointer to beta Exchange matrix
        shared_ptr<Matrix> Kb_;
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
        void form_J_DF_UHF();
        void form_K_PS_UHF();

    public:
        // Constructor for RHF/RKS
        PseudospectralHF(shared_ptr<BasisSet> basis, shared_ptr<Matrix> Da,
        shared_ptr<Matrix> Ja, shared_ptr<Matrix> Ka, shared_ptr<PSIO> psio, Options& opt);
        // Constructor for generic HF, with J, K, D to be set later
        PseudospectralHF(shared_ptr<BasisSet> basis, shared_ptr<PSIO> psio, Options& opt);
        // Destructor
        ~PseudospectralHF();
        // RAII, builds grids, sieves, fitting metric
        void common_init();
        
        // Setter methods, to be called from the JK functors
        void set_restricted(bool y_n) { restricted_ = y_n; }
        void set_Ja(shared_ptr<Matrix> Ja) {Ja_ = Ja;}
        void set_Jb(shared_ptr<Matrix> Jb) {Jb_ = Jb;}
        void set_Ka(shared_ptr<Matrix> Ka) {Ka_ = Ka;}
        void set_Kb(shared_ptr<Matrix> Kb) {Kb_ = Kb;}
        void set_Da(shared_ptr<Matrix> Da) {Da_ = Da;}
        void set_Db(shared_ptr<Matrix> Db) {Db_ = Db;}
        // form J and K
        void form_J_DF() { restricted_ ? form_J_DF_RHF() : form_J_DF_UHF(); }
        void form_K_PS() { restricted_ ? form_K_PS_RHF() : form_K_PS_UHF(); }
        void form_G_RHF();

};

}}

#endif
