#ifndef SCF_DF_H
#define SCF_DF_H

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
class IntegralFactory;
class PSIO;
class AIOHandler;

namespace scf {

class DFHFDiskIterator;
class DFHFLocalizer; 

class DFHF {

    private:

        // JK J buffers
        // Triangular D matrix
        double* Dtri_;
        // Triangular J matrix
        double* Jtri_;
        // d vector
        double* dQ_;

        // JK K buffers
        // Exchange tensor
        double** Ep_;
        // QS tensor chunk
        double*** QS_;
        // Ctemp matrix
        double*** Ctemp_;
        // Indexing
        int** m_ij_indices_;
        // More indexing
        int** n_indices_;
        // Even more indexing
        int* index_sizes_;

    protected:
        // Constant reference to the options object
        Options& options_;
       
        // Print or not?
        int print_;
        // Is this disk or core?
        bool is_disk_;
        // Is this object RI-J or RI-JK?
        bool is_jk_;
        // Is this object initialized?
        bool is_initialized_;
        // Is this object's disk iterator initialized?
        bool is_disk_initialized_;
        // Whether the alpha and beta orbitals are equal or not
        bool restricted_;
        
        // memory in doubles
        unsigned long int memory_;
        
        // Omega (0.0 if not used)
        double omega_;

        // Inverse of fitting metric
        boost::shared_ptr<FittingMetric> Jinv_;
        // The three index tensor (or a chunk if on disk)
        SharedMatrix Qmn_;

        // AO to USO transform matrix
        SharedMatrix AO2USO_;

        // SO Basis quantities
        // Shared pointer to alpha density matrix
        SharedMatrix Da_;
        // Shared pointer to beta density matrix
        SharedMatrix Db_;
        // Shared pointer to alpha occupation matrix
        SharedMatrix Ca_;
        // Shared pointer to beta occupation matrix
        SharedMatrix Cb_;
        // Shared pointer to total Coulomb matrix
        SharedMatrix Ja_;
        // Shared pointer to alpha Exchange matrix
        SharedMatrix Ka_;
        // Shared pointer to beta Exchange matrix
        SharedMatrix Kb_;
        // Number of alpha electrons
        const int* nalphapi_;
        // Number of beta electrons
        const int* nbetapi_;
       
        // AO Basis quantities
        // Shared pointer to alpha density matrix
        SharedMatrix Da_ao_;
        // Shared pointer to beta density matrix
        SharedMatrix Db_ao_;
        // Shared pointer to alpha occupation matrix
        SharedMatrix Ca_ao_;
        // Shared pointer to beta occupation matrix
        SharedMatrix Cb_ao_;
        // Shared pointer to total Coulomb matrix
        SharedMatrix Ja_ao_;
        // Shared pointer to alpha Exchange matrix
        SharedMatrix Ka_ao_;
        // Shared pointer to beta Exchange matrix
        SharedMatrix Kb_ao_;
        // Number of alpha electrons
        int nalpha_;
        // Number of beta electrons
        int nbeta_;
 
        // Zero basis set
        boost::shared_ptr<BasisSet> zero_;
        // Primary basis set
        boost::shared_ptr<BasisSet> primary_;
        // Auxiliary basis set
        boost::shared_ptr<BasisSet> auxiliary_;
        // Schwarz Sieve object
        boost::shared_ptr<SchwarzSieve> schwarz_;

        // The psio object
        boost::shared_ptr<PSIO> psio_;
        // The asynchronous, cyclic disk iterator
        boost::shared_ptr<DFHFDiskIterator> disk_iter_;
        // The unit number
        unsigned int unit_;

        // Setup
        void common_init();
        // Initialize methods
        void initialize();
        void initialize_JK_disk();
        void initialize_JK_core();
        void initialize_J_core();

        // Fast core J
        void compute_J_core();
        // Block methods
        void compute_JK_block_J(double** Qmnp, int nrows, int max_rows);
        void compute_JK_block_K(double** Qmnp, int nrows, int max_rows, bool is_alpha);

        // Transform C/D back to AO
        void USO2AO();
        // Transform J/K forward to USO
        void AO2USO();

    public:
        // Constructor for generic HF, with J, K, D to be set later
        DFHF(boost::shared_ptr<BasisSet> basis, boost::shared_ptr<PSIO> psio, Options& opt);
        // Erf Constructor for generic HF, with J, K, D to be set later
        DFHF(boost::shared_ptr<BasisSet> basis, boost::shared_ptr<PSIO> psio, Options& opt, double omega);
        // Destructor
        ~DFHF();

        // Set unit number
        void set_unit(unsigned int unit); 
        // Set memory
        void set_memory(unsigned long int memory);

        // Setter methods, to be called from the JK functors
        void set_jk(bool jk) { is_jk_ = jk; }
        void set_restricted(bool y_n) { restricted_ = y_n; }
        void set_J(SharedMatrix Ja) {Ja_ = Ja;}
        void set_Ka(SharedMatrix Ka) {Ka_ = Ka;}
        void set_Kb(SharedMatrix Kb) {Kb_ = Kb;}
        void set_Da(SharedMatrix Da) {Da_ = Da;}
        void set_Db(SharedMatrix Db) {Db_ = Db;}
        void set_Ca(SharedMatrix Ca) {Ca_ = Ca;}
        void set_Cb(SharedMatrix Cb) {Cb_ = Cb;}
        void set_Na(const int* Na) {nalphapi_ = Na;}
        void set_Nb(const int* Nb) {nbetapi_ = Nb;}

        // form J only
        void form_J_DF();
        // form J and K
        void form_JK_DF();

};

/*!
 * Disk iterator object for use with DFHF
 * Uses a cyclic algorithm to provide blocks
 *
 */
class DFHFDiskIterator {

    // The PSIO object
    boost::shared_ptr<PSIO> psio_;
    // The AIOHandler object
    boost::shared_ptr<AIOHandler> aio_;

    // The File
    unsigned int unit_;

    // Buffer A
    SharedMatrix A_;
    // Buffer B
    SharedMatrix B_;    

    // Fast index size
    int ntri_;
    // Slow index total size
    int naux_;
    // Max rows in each block (memory purposes)
    int max_rows_;
    // Currnet rows
    int current_rows_;

    // How many times has next_block been called?
    int iteration_;
  
    // number of blocks
    int nblocks_; 
    // Start row of each block
    std::vector<int> block_starts_;
    // Size in rows of each block
    std::vector<int> block_sizes_;
    // Block address order for this cycle
    std::vector<int> blocks_;
 
    // Initialize the block
    void common_init();
    // reset (not user called typically)
    void reset();
    // Post the read for a block
    void read(SharedMatrix A, int start, int rows);

public: 
    // Opens disk
    DFHFDiskIterator(boost::shared_ptr<PSIO>, int ntri, int naux, int max_rows);
    // Closes disk
    ~DFHFDiskIterator();

    // Set unit number
    void set_unit(unsigned int unit) { unit_ = unit; }

    // Call to get a block/Post theread for the next, if needed 
    double** next_block();
    // How many rows are there in the current block 
    int current_rows() const { return current_rows_; }
    // are we finished (not idempotent, resets)
    bool finished();
    // How many blocks are there?
    int nblock() const { return nblocks_; }
    // What is max_rows?
    int max_rows() const { return max_rows_; }
    // Sychronize the AIO
    void synchronize();

};


class DFHFLocalizer {

protected:
    /// The number of orbitals to work with 
    int nocc_;
    /// Debug? Defaults false
    bool debug_;

    /// The alpha parameter
    double alpha_;

    /// Speed tricks and cutoffs
    /// Variable threshold schwarz sieve
    boost::shared_ptr<SchwarzSieve> schwarz_;
    /// Charge cutoff (a.u.) for primary domain selection
    double charge_cutoff_;
    /// Charge leak allowed a.u.
    double charge_leak_cutoff_;
    /// Radial cutoff (bohr) for extended domain selection
    double R_ext_;
    /// number of am cardinal numbers to drop in the fitting bases
    int l_drop_;

    /// The primary basis 
    boost::shared_ptr<BasisSet> primary_;
    /// The auxiliary basis 
    boost::shared_ptr<BasisSet> auxiliary_;
    /// The zero basis
    boost::shared_ptr<BasisSet> zero_;

    /// The (mm|00) one-body integral constructor (note: not for conventional TEIs)
    boost::shared_ptr<IntegralFactory> factory_; 
    /// The (A0|mm) three-index integral constructor
    boost::shared_ptr<IntegralFactory> auxfactory_; 
    /// The (A0|A0) J matrix constructor
    boost::shared_ptr<IntegralFactory> Jfactory_; 

    /// The J matrix ERI
    boost::shared_ptr<TwoBodyAOInt> Jint_;
    /// The Amn matrix ERI
    boost::shared_ptr<TwoBodyAOInt> Amnint_;

    /// The localized orbitals (C1 obviously)
    SharedMatrix C_;
    /// The C1 S^1/2 matrix (for Lowdin charges)
    SharedMatrix Sp12_;
    /// The Lowdin charges (basis functions in rows, orbitals in cols)
    SharedMatrix Itemp_;
    /// The Lowdin charges (Atoms in rows, orbitals in cols)
    SharedMatrix I_;
    /// The custom int array for selection of domains (0-no, 1-charge, 2-distance, 3-unification, 4-alpha delocalized, 5-leak delocalized) (atoms x orbitals)
    int** domains_;
    /// The custom old int array for selection of domains (0-no, 1-charge, 2-distance, 3-unification, 4-alpha delocalized, 5-leak delocalized) (atoms x orbitals)
    int** old_domains_;
    /// Have domains changed/require recomputation?
    std::vector<bool> domains_changed_;
    /// Have superdomains changed/require recomputation? 
    std::vector<bool> superdomains_changed_; 

    /// Vector of Cholesky decompositions of J (or pointers where superdomains exist)
    std::vector<SharedMatrix > JCholesky_;

    /// Vector of significant shell indices (fast) for each orbital (slow) 
    std::vector<std::vector<int> > primary_shells_;
    /// Vector of significant function indices (fast) for each orbital (slow) 
    std::vector<std::vector<int> > primary_funs_;
    /// Vecotr of significant function starts per significant shell (fast) for each orbital (slow)
    std::vector<std::vector<int> > primary_starts_;
    /// Vector of significant shell indices (fast) for each orbital (slow) 
    std::vector<std::vector<int> > auxiliary_shells_;
    /// Vector of significant function indices (fast) for each orbital (slow) 
    std::vector<std::vector<int> > auxiliary_funs_;
    /// Vecotr of significant function starts per significant shell (fast) for each orbital (slow)
    std::vector<std::vector<int> > auxiliary_starts_;

    /// The vector of superdomains, which is indexes as superdomains[superdomain][member] = orbital member
    /// Superdomains are for localizable orbitals only 
    std::vector<std::vector<int> > superdomains_;
    /// The vector of orbitals in the completely delocalized set 
    std::set<int> delocal_set_;
    /// The vector of orbitals in the localizeable set
    std::set<int> local_set_;

    /// The common init routine
    void common_init();
    /// Computes the J matrix and associated Cholesky decomposition for orbital i
    SharedMatrix computeJCholesky(int i);
    /// Computes Amn for orbital i
    SharedMatrix computeAmn(int i);
    /// Computes the local K matrix for orbital i
    SharedMatrix computeK(int i, SharedMatrix J, SharedMatrix Amn);
    /// Find the primary/extended domains due to Lowdin charges
    void lowdinDomains();
    /// Identify completely delocalized domains 
    void delocalDomains();
    /// Populate the superdomain maps due to synergetic interdomain overlap
    void superDomains();
    /// Find the basis functions for each domain based on atoms and l_drop_ 
    void basisDomains();
    /// Check for changes in the atomic domains
    void checkForChanges();
    /// A silent global chnage has occurred, flag all the domains as stale 
    void globalChange();
    /// We've just computed everything, no need to repeat unless stuff changes 
    void globalReset();
    

public:
    DFHFLocalizer(int nocc, boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary);
    ~DFHFLocalizer();

    /// The local C matrix
    SharedMatrix C() const { return C_; }    
    /// The delocal orbital set
    std::set<int> local_set() const { return local_set_; }
    /// The local orbital set
    std::set<int> delocal_set() const { return delocal_set_; }

    /// Set the debug flag 
    void set_debug(bool d) { debug_ = d; }
    /// Set the Schwarz cutoff 
    void set_schwarz_cutoff(double schwarz); 
    /// Set the Lowdin charge cutoff (a.u.)
    void set_charge_cutoff(double cut) { charge_cutoff_ = cut; } 
    /// Set the Lowdin charge leak cutoff (a.u.)
    void set_charge_leak_cutoff(double cut) { charge_leak_cutoff_ = cut; } 
    /// Set the extended domains cutoff radius (bohr)
    void set_R_ext(double R) { R_ext_ = R; } 
    /// The fraction of atoms in a domain to declare the orbital delocalized
    void set_alpha(double a) { alpha_ = a; }
    /// Set the number of angular momenta to drop in the fitting basis (1 or 0 usually)
    void set_l_drop(int l) { l_drop_ = l; globalChange(); }

    /*!
    * Localize based on C1 square SPSD matrix D 
    * @param D density matrix in C1 AO basis
    */
    void choleskyLocalize(SharedMatrix D); 
    /// Computes the local part of matrix K, as determind by indices
    void computeKLocal(SharedMatrix Kglobal);

};

}}

#endif
