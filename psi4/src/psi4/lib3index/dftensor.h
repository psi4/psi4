/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef three_index_df_H
#define three_index_df_H

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class Molecule;
class IntVector;
class Vector3;

class FittingMetric {

protected:
    /// Pointer to the auxiliary basis set
    std::shared_ptr<BasisSet> aux_;
    /// Pointer to the poisson basis set
    std::shared_ptr<BasisSet> pois_;
    /// Is the metric poisson?
    bool is_poisson_;
    /// Should we force C1?
    bool force_C1_;
    /// Range separation omega (0.0 if not used)
    double omega_;

    /// The fitting metric or symmetric inverse
    SharedMatrix metric_;
    /// The indices (per irrep) of pivots
    std::shared_ptr<IntVector> pivots_;
    /// The indices (per irrep) of reverse pivots
    std::shared_ptr<IntVector> rev_pivots_;

    /// The fitting algorithm selected
    std::string algorithm_;
    /// Is the metric inverted or just a J matrix?
    bool is_inverted_;

    /// Fully pivot the fitting metric
    void pivot();

public:

    /// DF Fitting Metric
    FittingMetric(std::shared_ptr<BasisSet> aux, bool force_C1 = false);
    /// DF Fitting Metric
    FittingMetric(std::shared_ptr<BasisSet> aux, double omega, bool force_C1 = false);
    /// Poisson Fitting Metric
    FittingMetric(std::shared_ptr<BasisSet> aux, std::shared_ptr<BasisSet> pois, bool force_C1 = false);

    /// Destructor
    ~FittingMetric();

    /// What algorithm to use for symmetric inverse?
    std::string get_algorithm() const {return algorithm_; }
    /// Are poisson functions used?
    bool is_poisson() const {return is_poisson_; }
    /// Is the metric inverted?
    bool is_inverted() const {return is_inverted_; }

    /// The fitting metric or symmetric inverse
    SharedMatrix get_metric() const {return metric_; }
    /// The vector of pivots (for stability) (pivoted->global)
    std::shared_ptr<IntVector> get_pivots() const {return pivots_; }
    /// The vector of back pivots (for stability) (global->pivoted)
    std::shared_ptr<IntVector> get_reverse_pivots() const {return rev_pivots_; }

    /// The gaussian fitting basis
    std::shared_ptr<BasisSet> get_auxiliary_basis() const {return aux_; }
    /// The poisson fitting basis
    std::shared_ptr<BasisSet> get_poisson_basis() const {return pois_; }

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
    std::shared_ptr<Molecule> molecule_;
    /// Primary basis set
    std::shared_ptr<BasisSet> primary_;
    /// Dealias basis set
    std::shared_ptr<BasisSet> auxiliary_;
    /// options reference
    Options& options_;

    /// Symmetric inverse fitting metric
    SharedMatrix metric_;

    /// Full C matrix (must provide orthonormal MO basis)
    SharedMatrix C_;
    /// Active occupied C Matrix (for convenience)
    SharedMatrix Caocc_;
    /// Active virtual C Matrix (for convenience)
    SharedMatrix Cavir_;

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

    DFTensor(std::shared_ptr<BasisSet> primary,
             std::shared_ptr<BasisSet> auxiliary,
             SharedMatrix C,
             int nocc,
             int nvir,
             int naocc,
             int navir,
             Options& options);

    /**
    * Assumes all orbitals are active and pull options from enviroment
    **/
    DFTensor(std::shared_ptr<BasisSet> primary,
             std::shared_ptr<BasisSet> auxiliary,
             SharedMatrix C,
             int nocc,
             int nvir);
    ~DFTensor();

    SharedMatrix Qso();
    SharedMatrix Qmo();
    SharedMatrix Qoo();
    SharedMatrix Qov();
    SharedMatrix Qvv();

    SharedMatrix Imo();
    SharedMatrix Idfmo();
};

}
#endif
