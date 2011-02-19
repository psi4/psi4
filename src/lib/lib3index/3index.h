#ifndef three_index_H
#define three_index_H

#include <psi4-dec.h>
#include <psiconfig.h>

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class GridBlock;


class FittingMetric {

protected:
    /// Pointer to the auxiliary basis set
    shared_ptr<BasisSet> aux_;
    /// Pointer to the poisson basis set
    shared_ptr<BasisSet> pois_;
    /// Is the metric poisson?
    bool is_poisson_;

    /// The fitting metric or symmetric inverse
    shared_ptr<Matrix> metric_;
    /// The indices (per irrep) of pivots
    shared_ptr<IntVector> pivots_;
    /// The indices (per irrep) of reverse pivots 
    shared_ptr<IntVector> rev_pivots_;

    /// The fitting algorithm selected 
    std::string algorithm_;
    /// Is the metric inverted or just a J matrix?
    bool is_inverted_;

    /// Fully pivot the fitting metric
    void pivot();

public:

    /// Default constructor, for python
    FittingMetric();
    /// DF Fitting Metric
    FittingMetric(shared_ptr<BasisSet> aux);
    /// Poisson Fitting Metric
    FittingMetric(shared_ptr<BasisSet> aux, shared_ptr<BasisSet> pois);

    /// Destructor
    ~FittingMetric();
   
    /// What algorithm to use for symmetric inverse? 
    std::string get_algorithm() const {return algorithm_; }
    /// Are poisson functions used?
    bool is_poisson() const {return is_poisson_; }
    /// Is the metric inverted? 
    bool is_inverted() const {return is_inverted_; }

    /// The fitting metric or symmetric inverse
    shared_ptr<Matrix> get_metric() const {return metric_; }
    /// The vector of pivots (for stability) (pivoted->global)
    shared_ptr<IntVector> get_pivots() const {return pivots_; }
    /// The vector of back pivots (for stability) (global->pivoted)
    shared_ptr<IntVector> get_reverse_pivots() const {return rev_pivots_; }

    /// The gaussian fitting basis
    shared_ptr<BasisSet> get_auxiliary_basis() const {return aux_; }
    /// The poisson fitting basis
    shared_ptr<BasisSet> get_poisson_basis() const {return pois_; }

    /// Build the raw fitting metric (sets up indices to canonical)
    void form_fitting_metric();
    /// Build the Cholesky half inverse metric (calls form_fitting_metric)
    void form_cholesky_inverse();
    /// Build the QR half inverse metric (calls form_fitting_metric)
    void form_QR_inverse(double tol = 1.0E-10);
    /// Build the eigendecomposed half inverse metric (calls form_fitting_metric)
    void form_eig_inverse(double tol = 1.0E-10);
    /// Build the full inverse metric. NOT RECOMMENDED: Numerical stability (calls form_fitting_metric)
    void form_full_inverse();
    /// Build the full metric's Cholesky factor. RECOMMENDED: Numerical stability
    void form_cholesky_factor();
};

class PseudoGrid {

protected:
    /// The grid (in a Jeff Bridges voice)
    shared_ptr<GridBlock> grid_;
    /// Name (ie: cc-pVQZ-ultrafine)
    std::string name_;
    /// The molecule the grid is build on
    shared_ptr<Molecule> molecule_;
    
public:
    /// Construct a grid object from a molecule and name, then parse (or maybe autogen)
    PseudoGrid(shared_ptr<Molecule> mol, const std::string& name);
    /// Destructor, frees grid memory if created
    ~PseudoGrid();

    /// The GridBlock associated with this grid
    shared_ptr<GridBlock> getBlock() const { return grid_; } 
    /// Parse grid given full file path to G94-style grid file
    void parse(const std::string& file);
    
};

class SchwarzSieve {

protected:

    // The schwarz cutoff 
    double schwarz_;

    // Basis set for this schwarz
    shared_ptr<BasisSet> basis_;  
 
    // number of significant shell pairs 
    unsigned long int nshell_pairs_;
    // number of significant function pairs
    unsigned long int nfun_pairs_;

    double max_global_val_;
    int* schwarz_shells_;
    int* schwarz_funs_;

    long int* schwarz_shells_reverse_;
    long int* schwarz_funs_reverse_;

    double* schwarz_shell_vals_;
    double* schwarz_fun_vals_;

    void form_schwarz_ints();

public:
    SchwarzSieve(shared_ptr<BasisSet>, double cutoff);
    virtual ~SchwarzSieve();

    void form_schwarz_sieve(double cutoff);
    // Sizes of the significant bra/ket pairs
    unsigned long int get_nshell_pairs() const { return nshell_pairs_; }
    unsigned long int get_nfun_pairs() const { return nfun_pairs_; }
    // I_global = arr[2*I_local], J_global = arr[2*I_local + 1]
    // These are only defined up to nshell_pairs_ and nfun_pairs_, respectively
    int* get_schwarz_shells() const { return schwarz_shells_; }
    int* get_schwarz_funs() const { return schwarz_funs_; }
    // Canonical compound indexing, -1 if not present
    long int* get_schwarz_shells_reverse() const { return schwarz_shells_reverse_; }
    long int* get_schwarz_funs_reverse() const { return schwarz_funs_reverse_; }

};

class ThreeIndexChunk {
    
    protected:
        // Name
        std::string name_;
        // PSIO (for disk-based) 
        shared_ptr<PSIO> psio_;
        // PSIO address (for disk-based)
        psio_address address_;
        // Core tensor [0][0][0] pointer
        double* core_tensor_;

        // memory (doubles)
        unsigned long int memory_;
        // tensor size (doubles)
        unsigned long int tensor_size_;
        // max_rows permitted by memory
        int max_rows_;  
        // slow index total size
        int slow_size_; 
        // middle index total size
        int middle_size_; 
        // fast index total size
        int fast_size_; 
 
        // is this tensor core or disk?
        bool is_core_;
        // if disk, can we fully cache this tensor?
        bool is_cached_;

        // is this tensor finished?
        bool is_done_;
        // Starting slow index
        int current_row_;
        // Number of slow indices
        int current_rows_;
        // Tensor chunk
        double*chunk_;


    public:
        // Disk algorithm constructor
        ThreeIndexChunk(shared_ptr<PSIO> psio, 
                        const std::string& name,
                        int slow_size,
                        int middle_size,
                        int fast_size,
                        unsigned long int memory 
                        );
        // Core algorithm constructor
        ThreeIndexChunk(double* core_tensor, 
                        const std::string& name,
                        int slow_size,
                        int middle_size,
                        int fast_size,
                        unsigned long int memory 
                        );
        ~ThreeIndexChunk();

        void reset();
        bool isDone() const { return is_done_; }       
        void next();        
        double* pointer() const { return chunk_; } 
        int current_index() const { return current_row_; }
        int current_rows() const { return current_rows_; }

        bool is_core() const { return is_core_; }
        std::string name() const { return name_; }
        int slow_size() const { return slow_size_; }        
        int middle_size() const { return middle_size_; }        
        int fast_size() const { return fast_size_; }        
        unsigned long int memory() const { return memory_; }
        int max_rows() const { return max_rows_; } 
};


class DFTensor {

protected:

    // The fitting metric (if still needed)
    shared_ptr<FittingMetric> metric_;
    // The Schwarz sieve 
    shared_ptr<SchwarzSieve> schwarz_;
    
    // The primary basis
    shared_ptr<BasisSet> primary_;
    // The auxiliary basis
    shared_ptr<BasisSet> auxiliary_;

    // The PSIO object
    shared_ptr<PSIO> psio_;
    // is this tensor core or disk
    bool is_core_;

    // Number of occupied orbitals
    int nocc_;
    // Number of virtual orbitals
    int nvir_;

public:
    DFTensor(shared_ptr<PSIO>, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> aux, double schwarz = 0.0);
    virtual ~DFTensor();

    void common_init();     
    
    shared_ptr<SchwarzSieve> get_schwarz() const { return schwarz_; }
    shared_ptr<FittingMetric> get_metric() const { return metric_; }

    // Form all MO integrals, disk algorithm, default striping Qov
    void form_MO_disk(shared_ptr<Matrix> Cocc, shared_ptr<Matrix> Cvir, const std::string& algorithm, double cond);
    
    // Iterators to the various blocks of the MO DF integrals, 
    // memory in doubles  
    shared_ptr<ThreeIndexChunk> get_oo_iterator(unsigned long int memory);
    shared_ptr<ThreeIndexChunk> get_vv_iterator(unsigned long int memory);
    shared_ptr<ThreeIndexChunk> get_ov_iterator(unsigned long int memory);

};

class Pseudospectral {

protected:
    shared_ptr<BasisSet> primary_;
    shared_ptr<BasisSet> dealias_;
    shared_ptr<PseudoGrid> grid_;
    int npoints_;
    shared_ptr<PSIO> psio_;
public:
    Pseudospectral(shared_ptr<PSIO>, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> dealias, shared_ptr<PseudoGrid> grid);
    ~Pseudospectral();
    
    shared_ptr<Matrix> form_X();
    shared_ptr<Matrix> form_X_dealias();
    shared_ptr<Matrix> form_S();
    shared_ptr<Matrix> form_S_dealias();
    shared_ptr<Matrix> form_Q();
    shared_ptr<Matrix> form_A();
    void form_A_disk();
    int npoints() const { return npoints_; }

    shared_ptr<Matrix> form_I(); 

};
class Denominator {

protected:
    // Denominator (ia in columns, w in rows)
    shared_ptr<Matrix> denominator_;

    // Pointer to active occupied orbitals
    shared_ptr<Vector> eps_occ_;
    // Pointer to active virtual orbitals
    shared_ptr<Vector> eps_vir_;

    // Gauge reference (HOMO-LUMO splitting)
    double gauge_;

    virtual void decompose() = 0; 
public:
    Denominator(shared_ptr<Vector> eps_occ_, shared_ptr<Vector> eps_vir);
    virtual ~Denominator();
};

class LaplaceDenominator : public Denominator {

protected:
    int nvector_;    
    void decompose();
public:
    LaplaceDenominator(shared_ptr<Vector> eps_occ_, shared_ptr<Vector> eps_vir, int nvector);
    ~LaplaceDenominator();

};

class CholeskyDenominator : public Denominator {

protected:
    double delta_;    
    double degeneracy_multiplier_;    
    void decompose();
    static bool criteria(std::pair<int, double> a, std::pair<int, double> b);
public:
    CholeskyDenominator(shared_ptr<Vector> eps_occ_, shared_ptr<Vector> eps_vir, double delta, double degeneracy_multiplier = 100.0);
    ~CholeskyDenominator();

};

}
#endif
