#ifndef three_index_H
#define three_index_H

#include <psi4-dec.h>
#include <psiconfig.h>

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class IntVector;


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
    unsigned long int get_nshell_pairs() const { return nshell_pairs_; }
    unsigned long int get_nfun_pairs() const { return nfun_pairs_; }
    // I_global = arr[2*I_local], J_global = arr[2*I_local + 1]
    // These are only defined up to nshell_pairs_ and nfun_pairs_, respectively
    int* get_schwarz_shells() const { return schwarz_shells_; }
    int* get_schwarz_funs() const { return schwarz_funs_; }
    // Canonical compound indexingi, -1 if not present
    long int* get_schwarz_shells_reverse() const { return schwarz_shells_reverse_; }
    long int* get_schwarz_funs_reverse() const { return schwarz_funs_reverse_; }

};



}
#endif
