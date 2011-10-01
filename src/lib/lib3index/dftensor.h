#ifndef three_index_df_H
#define three_index_df_H

#include <libmints/mints.h>

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class Vector3;

class FittingMetric {

protected:
    /// Pointer to the auxiliary basis set
    boost::shared_ptr<BasisSet> aux_;
    /// Pointer to the poisson basis set
    boost::shared_ptr<BasisSet> pois_;
    /// Is the metric poisson?
    bool is_poisson_;
    /// Should we force C1?
    bool force_C1_;
    /// Range separation omega (0.0 if not used)
    double omega_;

    /// The fitting metric or symmetric inverse
    boost::shared_ptr<Matrix> metric_;
    /// The indices (per irrep) of pivots
    boost::shared_ptr<IntVector> pivots_;
    /// The indices (per irrep) of reverse pivots
    boost::shared_ptr<IntVector> rev_pivots_;

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
    FittingMetric(boost::shared_ptr<BasisSet> aux, bool force_C1 = false);
    /// DF Fitting Metric
    FittingMetric(boost::shared_ptr<BasisSet> aux, double omega, bool force_C1 = false);
    /// Poisson Fitting Metric
    FittingMetric(boost::shared_ptr<BasisSet> aux, boost::shared_ptr<BasisSet> pois, bool force_C1 = false);

    /// Destructor
    ~FittingMetric();

    /// What algorithm to use for symmetric inverse?
    std::string get_algorithm() const {return algorithm_; }
    /// Are poisson functions used?
    bool is_poisson() const {return is_poisson_; }
    /// Is the metric inverted?
    bool is_inverted() const {return is_inverted_; }

    /// The fitting metric or symmetric inverse
    boost::shared_ptr<Matrix> get_metric() const {return metric_; }
    /// The vector of pivots (for stability) (pivoted->global)
    boost::shared_ptr<IntVector> get_pivots() const {return pivots_; }
    /// The vector of back pivots (for stability) (global->pivoted)
    boost::shared_ptr<IntVector> get_reverse_pivots() const {return rev_pivots_; }

    /// The gaussian fitting basis
    boost::shared_ptr<BasisSet> get_auxiliary_basis() const {return aux_; }
    /// The poisson fitting basis
    boost::shared_ptr<BasisSet> get_poisson_basis() const {return pois_; }

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
    /// Build the full inverse metric.
    void form_full_eig_inverse(double tol = 1.0E-10);
    /// Build the full metric's Cholesky factor. RECOMMENDED: Numerical stability
    void form_cholesky_factor();
};

class DFTensor {

protected:

    /// Debug level
    int debug_;
    /// Print level
    int print_;

    /// Molecule (fo convenience)
    boost::shared_ptr<Molecule> molecule_;
    /// Primary basis set
    boost::shared_ptr<BasisSet> primary_;
    /// Dealias basis set
    boost::shared_ptr<BasisSet> auxiliary_;
    /// options reference
    Options& options_;

    /// Symmetric inverse fitting metric
    boost::shared_ptr<Matrix> metric_;

    /// Full C matrix (must provide orthonormal MO basis)
    boost::shared_ptr<Matrix> C_;
    /// Active occupied C Matrix (for convenience)
    boost::shared_ptr<Matrix> Caocc_;
    /// Active virtual C Matrix (for convenience)
    boost::shared_ptr<Matrix> Cavir_;

    /// Number of AO primary functions
    int nso_;
    /// Number of MO primary functions
    int nmo_;

    /// Number of grid points
    int naux_;

    /// Number of frozen occupieds
    int nfocc_;
    /// Total number of occupieds
    int nocc_;
    /// Number of active occupieds
    int naocc_;
    /// Number of frozen virtuals
    int nfvir_;
    /// Total number of virtuals
    int nvir_;
    /// Number of active virtuals
    int navir_;

    void common_init();
    void build_metric();
    void print_header();

public:

    DFTensor(boost::shared_ptr<BasisSet> primary,
             boost::shared_ptr<BasisSet> auxiliary,
             boost::shared_ptr<Matrix> C,
             int nocc,
             int nvir,
             int naocc,
             int navir,
             Options& options);
    ~DFTensor();

    boost::shared_ptr<Matrix> Qso();
    boost::shared_ptr<Matrix> Qmo();
    boost::shared_ptr<Matrix> Qoo();
    boost::shared_ptr<Matrix> Qov();
    boost::shared_ptr<Matrix> Qvv();

    boost::shared_ptr<Matrix> Imo();
    boost::shared_ptr<Matrix> Idfmo();
};

}
#endif
