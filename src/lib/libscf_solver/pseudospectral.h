#ifndef SCF_PSEUDO_H
#define SCF_PSEUDO_H

namespace boost {
template<class T> class shared_ptr;
}

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
        // Print or not?
        int print_;
        // Whether the alpha and beta orbitals are equal or not
        bool restricted_;
        // The psio object
        boost::shared_ptr<PSIO> psio_;
        // Threadsafe (A|mn) integral objects
        std::vector<boost::shared_ptr<TwoBodyAOInt> > eri_;
        // Threadsafe (m|V|n) integral objects
        std::vector<boost::shared_ptr<PseudospectralInt> > pot_;
        // Inverse of fitting metric (Coulomb side)
        boost::shared_ptr<FittingMetric> Jinv_;
        // X pseudospectral matrix
        SharedMatrix X_;
        // Shared pointer to alpha density matrix
        SharedMatrix Da_;
        // Shared pointer to beta density matrix
        SharedMatrix Db_;
        // Shared pointer to total Coulomb matrix
        SharedMatrix Ja_;
        // Shared pointer to alpha Exchange matrix
        SharedMatrix Ka_;
        // Shared pointer to beta Exchange matrix
        SharedMatrix Kb_;
        // Constant reference to the options object
        Options& options_;
        // Primary basis set
        boost::shared_ptr<BasisSet> primary_;
        // Auxiliary basis set (coulomb)
        boost::shared_ptr<BasisSet> auxiliary_;
        // Quadrature Points (P x 3)
        SharedMatrix points_;
        // Schwarz Sieve object
        boost::shared_ptr<SchwarzSieve> schwarz_;

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
        PseudospectralHF(boost::shared_ptr<BasisSet> basis, SharedMatrix Da,
        SharedMatrix Ja, SharedMatrix Ka, boost::shared_ptr<PSIO> psio, Options& opt);
        // Constructor for generic HF, with J, K, D to be set later
        PseudospectralHF(boost::shared_ptr<BasisSet> basis, boost::shared_ptr<PSIO> psio, Options& opt);
        // Destructor
        ~PseudospectralHF();
        // RAII, builds grids, sieves, fitting metric
        void common_init();

        // Setter methods, to be called from the JK functors
        void set_restricted(bool y_n) { restricted_ = y_n; }
        void set_J(SharedMatrix Ja) {Ja_ = Ja;}
        void set_Ka(SharedMatrix Ka) {Ka_ = Ka;}
        void set_Kb(SharedMatrix Kb) {Kb_ = Kb;}
        void set_Da(SharedMatrix Da) {Da_ = Da;}
        void set_Db(SharedMatrix Db) {Db_ = Db;}
        // form J and K
        void form_J_DF() { restricted_ ? form_J_DF_RHF() : form_J_DF_UHF(); }
        void form_K_PS() { restricted_ ? form_K_PS_RHF() : form_K_PS_UHF(); }
        void form_G_RHF();

};

}}

#endif
