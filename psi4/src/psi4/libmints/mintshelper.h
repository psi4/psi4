/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_helper_h
#define _psi_src_lib_libmints_helper_h

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/multipolesymmetry.h"
#include "psi4/libpsi4util/process.h"

#include <vector>

namespace psi {

class CdSalcList;
class CorrelationFactor;
class TwoBodyAOInt;
class PetiteList;
class ThreeCenterOverlapInt;
class OneBodyAOInt;

/**
 * The MintsHelper object, places molecular integrals
 * (and later derivative integrals) on disk
 **/
class PSI_API MintsHelper {
   private:
    /// The Options reference for basis sets and things
    Options& options_;
    std::shared_ptr<PSIO> psio_;
    std::shared_ptr<MatrixFactory> factory_;
    std::shared_ptr<Molecule> molecule_;
    std::shared_ptr<IntegralFactory> integral_;

    /// The ORBITAL basis
    std::shared_ptr<BasisSet> basisset_;

    /// DF/RI/F12/etc basis sets
    std::map<std::string, std::shared_ptr<BasisSet>> basissets_;

    std::shared_ptr<SOBasisSet> sobasis_;
    std::shared_ptr<TwoBodyAOInt> eriInts_;

    int print_;
    int nthread_;

    std::map<std::pair<std::string, bool>, SharedMatrix> cached_oe_ints_;

    /// Value which any two-electron integral is below is discarded
    double cutoff_;

    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix Co, SharedMatrix Cv);
    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// In-core builds spin eri's
    SharedMatrix mo_spin_eri_helper(SharedMatrix Iso, int n1, int n2);

    SharedMatrix ao_helper(const std::string& label, std::shared_ptr<TwoBodyAOInt> ints);
    SharedMatrix ao_shell_getter(const std::string& label, std::shared_ptr<TwoBodyAOInt> ints, int M, int N, int P,
                                 int Q);

    SharedMatrix ao_3coverlap_helper(const std::string& label, std::shared_ptr<ThreeCenterOverlapInt> ints);

    void common_init();

    /**
     * @brief Compute a one-body operator matrix.
     *
     * This function takes in a vector of OneBodyAOInt objects and outputs a matrix
     * representation of the one-body operator.
     *
     * @param[in] ints Vector of OneBodyAOInt integrals
     * @param[out] out Matrix containing the one-body operator
     * @param[in] symm Use symmetry flag
     */
    void one_body_ao_computer(std::vector<std::shared_ptr<OneBodyAOInt>> ints, SharedMatrix out, bool symm);
    void grad_two_center_computer(std::vector<std::shared_ptr<OneBodyAOInt>> ints, SharedMatrix D, SharedMatrix out);
    /// Helper function to convert ao integrals to so and cache them
    void cache_ao_to_so_ints(SharedMatrix ao_ints, const std::string& label, bool include_perturbation);
    /// Returns true if an integral type is already computed and cached
    bool are_ints_cached(const std::string& label, bool include_perturbation);

    /// Computes X2C overlap, kinetic, and potential integrals
    void compute_so_x2c_ints(bool include_perturbations = true);
    /// Add dipole perturbation to the potential integrals
    void add_dipole_perturbation(SharedMatrix potential_mat);
    /// Returns the non-relativistic overlap integrals in the so basis
    SharedMatrix so_overlap_nr();
    /// Returns the non-relativistic kinetic integrals in the so basis
    SharedMatrix so_kinetic_nr();
    /// Returns the non-relativistic potential integrals in the so basis
    SharedMatrix so_potential_nr(bool include_perturbations = true);

   public:
    /// Class initialization from Wavefunction
    void init_helper(std::shared_ptr<Wavefunction> wavefunction = std::shared_ptr<Wavefunction>());
    /// Class initialization from orbital BasisSet and a map of BasisSet objects
    void init_helper(std::shared_ptr<BasisSet> basis, std::map<std::string, std::shared_ptr<psi::BasisSet>> basissets =
                                                          std::map<std::string, std::shared_ptr<psi::BasisSet>>());

    /// Constructor, using basisset
    MintsHelper(std::shared_ptr<BasisSet> basis, Options& options = Process::environment.options, int print = 0);
    MintsHelper(std::shared_ptr<BasisSet> basis, std::map<std::string, std::shared_ptr<psi::BasisSet>> basissets,
                Options& options = Process::environment.options, int print = 0);

    /// Constructor, using wavefunction
    MintsHelper(std::shared_ptr<Wavefunction> wavefunction);
    /// Destructor, does nothing
    ~MintsHelper();

    OperatorSymmetry operator_symmetry(int order) { return OperatorSymmetry(order, molecule_, integral_, factory_); }

    /// Returns the number of basis functions
    int nbf() const;

    /// Sets the print level
    void set_print(int print) { print_ = print; }
    void set_nthread(int nthread) { nthread_ = nthread; }

    /// Returns petite list that is capable of transforming basis functions (nbf) to SO's.
    std::shared_ptr<PetiteList> petite_list() const;

    enum { kFromCartesianAO = true, kFromBF = false };
    /** Returns petite list that is capable of transforming AO basis functions (nbf) to SO's.
     *  \param include_pure_transform Is either kFromCartesianAO or kFromBF.
     */
    std::shared_ptr<PetiteList> petite_list(bool include_pure_transform) const;
    /// Basis set being used.
    std::shared_ptr<BasisSet> basisset() const;
    /// SO basis set being used.
    std::shared_ptr<SOBasisSet> sobasisset() const;
    /// Matrix factory being used
    std::shared_ptr<MatrixFactory> factory() const;
    /// Integral factory being used
    std::shared_ptr<IntegralFactory> integral() const;

    /// Getters and setters for other basis sets
    std::map<std::string, std::shared_ptr<BasisSet>> basissets() const { return basissets_; };
    std::shared_ptr<BasisSet> get_basisset(std::string label);
    void set_basisset(std::string label, std::shared_ptr<BasisSet> basis);
    bool basisset_exists(std::string label);

    /// Molecular integrals (just like cints used to do)
    void integrals();
    void integrals_erf(double w = -1.0);
    void integrals_erfc(double w = -1.0);

    /// Standard one electron integrals (just like oeints used to do)
    void one_electron_integrals();
    /// Derivative integrals (not implemented)
    void integral_gradients();
    /// Hessian integrals (not implemented)
    void integral_hessians();

    /// AO ERI Integrals (Full matrix, not recommended for large systems)
    SharedMatrix ao_eri(std::shared_ptr<IntegralFactory> = nullptr);
    SharedMatrix ao_eri(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3,
                        std::shared_ptr<BasisSet> bs4);
    /// AO ERI Shell
    SharedMatrix ao_eri_shell(int M, int N, int P, int Q);

    // Derivatives of OEI in AO and MO basis
    std::vector<SharedMatrix> ao_oei_deriv1(const std::string& oei_type, int atom);
    std::vector<SharedMatrix> ao_oei_deriv2(const std::string& oei_type, int atom1, int atom2);
    std::vector<SharedMatrix> ao_overlap_half_deriv1(const std::string& half_der_side, int atom);
    std::vector<SharedMatrix> ao_overlap_kinetic_deriv1_helper(const std::string& type, int atom);
    std::vector<SharedMatrix> ao_overlap_half_deriv1_helper(const std::string& half_der_side, int atom);
    std::vector<SharedMatrix> ao_potential_deriv1_helper(int atom);
    std::vector<SharedMatrix> ao_overlap_kinetic_deriv2_helper(const std::string& type, int atom1, int atom2);
    std::vector<SharedMatrix> ao_potential_deriv2_helper(int atom1, int atom2);
    std::vector<SharedMatrix> mo_oei_deriv1(const std::string& oei_type, int atom, SharedMatrix C1, SharedMatrix C2);
    std::vector<SharedMatrix> mo_oei_deriv2(const std::string& oei_type, int atom1, int atom2, SharedMatrix C1,
                                            SharedMatrix C2);
    std::vector<SharedMatrix> mo_overlap_half_deriv1(const std::string& half_der_side, int atom, SharedMatrix C1,
                                                     SharedMatrix C2);

    // Derivatives of electric dipole moment integrals in AO and MO basis
    std::vector<SharedMatrix> ao_elec_dip_deriv1(int atom);
    std::vector<SharedMatrix> ao_elec_dip_deriv1_helper(int atom);
    std::vector<SharedMatrix> mo_elec_dip_deriv1(int atom, SharedMatrix C1, SharedMatrix C2);

    // Derivatives of TEI in AO and MO basis
    std::vector<SharedMatrix> ao_tei_deriv1(int atom, double omega = 0.0, std::shared_ptr<IntegralFactory> = nullptr);
    std::vector<SharedMatrix> ao_tei_deriv2(int atom1, int atom2);
    std::vector<SharedMatrix> mo_tei_deriv1(int atom, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3,
                                            SharedMatrix C4);
    std::vector<SharedMatrix> mo_tei_deriv2(int atom1, int atom2, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3,
                                            SharedMatrix C4);
    std::vector<SharedMatrix> ao_metric_deriv1(int atom, const std::string& aux_name);
    std::vector<SharedMatrix> ao_3center_deriv1(int atom, const std::string& aux_name);

    /// AO ERF Integrals
    SharedMatrix ao_erf_eri(double omega, std::shared_ptr<IntegralFactory> = nullptr);
    /// AO ERFC Omega Integrals
    SharedMatrix ao_erfc_eri(double omega);
    /// AO F12 Integrals
    SharedMatrix ao_f12(std::vector<std::pair<double, double>> exp_coeff);
    SharedMatrix ao_f12(std::vector<std::pair<double, double>> exp_coeff, std::shared_ptr<BasisSet> bs1,
                        std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);
    /// AO F12 Integrals
    SharedMatrix ao_f12_scaled(std::vector<std::pair<double, double>> exp_coeff);
    SharedMatrix ao_f12_scaled(std::vector<std::pair<double, double>> exp_coeff, std::shared_ptr<BasisSet> bs1,
                               std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3,
                               std::shared_ptr<BasisSet> bs4);
    /// AO F12 squared Integrals
    SharedMatrix ao_f12_squared(std::vector<std::pair<double, double>> exp_coeff);
    SharedMatrix ao_f12_squared(std::vector<std::pair<double, double>> exp_coeff, std::shared_ptr<BasisSet> bs1,
                                std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3,
                                std::shared_ptr<BasisSet> bs4);
    /// AO F12G12 Integrals
    SharedMatrix ao_f12g12(std::vector<std::pair<double, double>> exp_coeff);
    SharedMatrix ao_f12g12(std::vector<std::pair<double, double>> exp_coeff, std::shared_ptr<BasisSet> bs1,
                        std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);
    /// AO F12 double commutator Integrals
    SharedMatrix ao_f12_double_commutator(std::vector<std::pair<double, double>> exp_coeff);
    SharedMatrix ao_f12_double_commutator(std::vector<std::pair<double, double>> exp_coeff, std::shared_ptr<BasisSet> bs1,
                        std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);
    /// F12 Fitted Slater Correlation Factor
    std::vector<std::pair<double, double>> f12_cgtg(double exponent = 1.0);

    /// 3Center overlap integrals
    SharedMatrix ao_3coverlap();
    SharedMatrix ao_3coverlap(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                              std::shared_ptr<BasisSet> bs3);

    /// Erf-attenuated Coulomb potential on origin
    SharedMatrix ao_potential_erf(const std::vector<double> &origin, double omega = 0.0, int deriv = 0);
    /// Erfc-attenuated Coulomb potential on origin
    SharedMatrix ao_potential_erf_complement(const std::vector<double> &origin, double omega = 0.0, int deriv = 0);

    /// Symmetric MO ERI Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    SharedMatrix mo_eri(SharedMatrix Cocc, SharedMatrix Cvir);
    /// Non Symmetric MO ERI Omega Integrals, (12|34) type  (Full matrix, N^5, not recommended for large systems)
    SharedMatrix mo_eri(SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO ERI Omega Integrals (Full matrix, not recommended for large systems)
    SharedMatrix mo_erf_eri(double omega, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO ERFC Omega Integrals
    SharedMatrix mo_erfc_eri(double omega, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO F12 Integrals
    SharedMatrix mo_f12(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1, SharedMatrix C2,
                        SharedMatrix C3, SharedMatrix C4);
    /// MO F12 squared Integrals
    SharedMatrix mo_f12_squared(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1, SharedMatrix C2,
                                SharedMatrix C3, SharedMatrix C4);
    /// MO F12G12 Integrals
    SharedMatrix mo_f12g12(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1, SharedMatrix C2,
                           SharedMatrix C3, SharedMatrix C4);
    /// MO F12 double commutator Integrals
    SharedMatrix mo_f12_double_commutator(std::vector<std::pair<double, double>> exp_coeff, SharedMatrix C1,
                                          SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);

    /// Symmetric MO ERI Omega Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    SharedMatrix mo_erf_eri(double omega, SharedMatrix Cocc, SharedMatrix Cvir);

    /// Symmetric MO Spin ERI Integrals, <oo|vv> type (Full matrix (16x larger than MO ERI), N^5,
    /// most definitely not recommended for large systems)
    /// Pass C_ C_ for <aa|aa> type, Cocc_, Cocc_ for <oo|oo> type, or Cvir_, Cvir_ for <vv|vv> type
    SharedMatrix mo_spin_eri(SharedMatrix Co, SharedMatrix Cv);

    /// AO Overlap Integrals
    SharedMatrix ao_overlap();
    SharedMatrix ao_overlap(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    /// AO Kinetic Integrals
    SharedMatrix ao_kinetic();
    SharedMatrix ao_kinetic(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    /// AO Potential Integrals
    SharedMatrix ao_potential();
    SharedMatrix ao_potential(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    /// AO ECP Integrals
#ifdef USING_ecpint
    SharedMatrix ao_ecp();
    SharedMatrix ao_ecp(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
#endif
    /// AO pVp Integrals
    SharedMatrix ao_pvp();
    /// AO DKH Integrals
    SharedMatrix ao_dkh(int dkh_order = -1);
    /// SO DKH Integrals
    SharedMatrix so_dkh(int dkh_order = -1);
    /// Vector AO Dipole Integrals
    std::vector<SharedMatrix> ao_dipole();
    /// Vector AO Quadrupole Integrals
    std::vector<SharedMatrix> ao_quadrupole();
    /// Vector AO Traceless Quadrupole Integrals
    std::vector<SharedMatrix> ao_traceless_quadrupole();
    /// Vector AO Multipole Integrals up to given order (in CCA lexicographic order)
    std::vector<SharedMatrix> ao_multipoles(int order, const std::vector<double>& origin);
    /// AO EFP Multipole Potential Integrals
    std::vector<SharedMatrix> ao_efp_multipole_potential(const std::vector<double>& origin,
                                                         int deriv = 0);
    // AO Multipole Potential Integrals up to given order (in CCA lexicographic order)
    std::vector<SharedMatrix> ao_multipole_potential(int order, const std::vector<double>& origin, int deriv = 0);
    /// Electric Field Integrals
    std::vector<SharedMatrix> electric_field(const std::vector<double>& origin, int deriv = 0);
    /// Induction Operator for dipole moments at given sites
    SharedMatrix induction_operator(SharedMatrix coords, SharedMatrix moment);
    /// Electric Field Value at given sites
    SharedMatrix electric_field_value(SharedMatrix coords, SharedMatrix D);
    /// Vector AO Angular Momentum Integrals
    std::vector<SharedMatrix> ao_angular_momentum();
    /// Vector AO Nabla Integrals
    std::vector<SharedMatrix> ao_nabla();
    /// SO Overlap Integrals
    SharedMatrix so_overlap(bool include_perturbations = true);
    /// SO Kinetic Integrals
    SharedMatrix so_kinetic(bool include_perturbations = true);
    /// SO ECP Integrals
#ifdef USING_ecpint
    SharedMatrix so_ecp();
#endif
    /// SO Potential Integrals
    SharedMatrix so_potential(bool include_perturbations = true);
    /// Vector SO Dipole Integrals
    std::vector<SharedMatrix> so_dipole() const;
    /// Vector SO Nabla Integrals
    std::vector<SharedMatrix> so_nabla() const;
    /// Vector SO Angular Momentum Integrals
    std::vector<SharedMatrix> so_angular_momentum() const;
    /// Vector SO Quadrupole Integrals
    std::vector<SharedMatrix> so_quadrupole();
    /// Vector SO Traceless Quadrupole Integrals
    std::vector<SharedMatrix> so_traceless_quadrupole();

    /// Electrostatic potential values at given sites with associated charge, specified as an (n_sites, 4) matrix
    SharedVector electrostatic_potential_value(SharedVector charges, SharedMatrix coords, SharedMatrix D);

    /// Returns a CdSalcList object
    std::shared_ptr<CdSalcList> cdsalcs(int needed_irreps = 0xF, bool project_out_translations = true,
                                        bool project_out_rotations = true);

    /// N^5 ao->mo transform, in memory, smart indexing
    SharedMatrix mo_transform(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);

    /// Gradient Integrals
    SharedMatrix overlap_grad(SharedMatrix D);

    // Computes all "core" gradient terms T + V + perturb
    SharedMatrix core_hamiltonian_grad(SharedMatrix D);
    // Computes the metric derivative gradient terms for DF methods
    // Uses a vector of "densities" for methods that decompose the "densities"
    std::map<std::string, SharedMatrix> metric_grad(std::map<std::string, SharedMatrix>& D,
                                                    const std::string& aux_name);
    SharedMatrix three_idx_grad(const std::string& aux_name, const std::string& intermed_name,
                                const std::string& gradient_name);

    // Computes the gradient due to the kinetic energy term of the core Hamiltonian
    SharedMatrix kinetic_grad(SharedMatrix D);
    // Computes the gradient due to the potential energy term of the core Hamiltonian
    SharedMatrix potential_grad(SharedMatrix D);
    // Computes the electric dipole derivatives
    SharedMatrix dipole_grad(SharedMatrix D);
    SharedMatrix multipole_grad(SharedMatrix D, int order, const std::vector<double>& origin = {0.0, 0.0, 0.0});
    // Computes the gradient due to external dipole perturbation
    SharedMatrix perturb_grad(SharedMatrix D);
    // Computes the gradient due to ECPs in the basis set
#ifdef USING_ecpint
    SharedMatrix effective_core_potential_grad(SharedMatrix D);
#endif

    /// Play function
    void play();
};
}  // namespace psi

#endif
