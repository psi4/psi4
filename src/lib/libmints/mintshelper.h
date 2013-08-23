/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmints_helper_h
#define _psi_src_lib_libmints_helper_h

#include <vector>
#include "wavefunction.h"
#include "multipolesymmetry.h"

namespace psi {

class Options;
class CdSalcList;

/**
* The MintsHelper object, places molecular integrals
* (and later derivative integrals) on disk
**/
class MintsHelper {

private:
    /// The Options reference for basis sets and things
    Options& options_;
    boost::shared_ptr<PSIO> psio_;
    boost::shared_ptr<MatrixFactory> factory_;
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<IntegralFactory> integral_;
    boost::shared_ptr<BasisSet> basisset_;
    boost::shared_ptr<SOBasisSet> sobasis_;
    boost::shared_ptr<TwoBodyAOInt> eriInts_;
    int print_;

    /// Value which any two-electron integral is below is discarded
    double cutoff_;

    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix Co, SharedMatrix Cv);
    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2,
                                                 SharedMatrix C3, SharedMatrix C4);

    SharedMatrix ao_helper(const std::string& label, boost::shared_ptr<TwoBodyAOInt> ints);
    SharedMatrix ao_shell_getter(const std::string& label, boost::shared_ptr<TwoBodyAOInt> ints, int M, int N, int P, int Q);


    void common_init();

public:

    void init_helper(boost::shared_ptr<Wavefunction> wavefunction = boost::shared_ptr<Wavefunction>());
    void init_helper(boost::shared_ptr<BasisSet> basis);

    /// Constructor, just lines references up
    MintsHelper(Options&, int print = 1);
    /// Constructor, uses a basis set
    MintsHelper(boost::shared_ptr<BasisSet> basis);
    /// Constructor, uses globals
    MintsHelper();
    /// Constructor, using wavefunction
    MintsHelper(boost::shared_ptr<Wavefunction> wavefunction);
    /// Destructor, does nothing
    ~MintsHelper();

    OperatorSymmetry operator_symmetry(int order) {
        return OperatorSymmetry(order, molecule_, integral_, factory_);
    }

    /// Returns the number of basis functions
    int nbf() const;

    /// Returns petite list that is capable of transforming basis functions (nbf) to SO's.
    boost::shared_ptr<PetiteList> petite_list() const;

    enum {
        kFromCartesianAO = true,
        kFromBF = false
    };
    /** Returns petite list that is capable of transforming AO basis functions (nbf) to SO's.
     *  \param include_pure_transform Is either kFromCartesianAO or kFromBF.
     */
    boost::shared_ptr<PetiteList> petite_list(bool include_pure_transform) const;
    /// Basis set being used.
    boost::shared_ptr<BasisSet> basisset() const;
    /// SO basis set being used.
    boost::shared_ptr<SOBasisSet> sobasisset() const;
    /// Matrix factory being used
    boost::shared_ptr<MatrixFactory> factory() const;
    /// Integral factory being used
    boost::shared_ptr<IntegralFactory> integral() const;

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
    SharedMatrix ao_eri();
    /// AO ERI Shell
    SharedMatrix ao_eri_shell(int M, int N, int P, int Q);
    /// AO ERF Integrals
    SharedMatrix ao_erf_eri(double omega);
    /// MO ERFC Omega Integrals
    SharedMatrix ao_erfc_eri(double omega);
    /// MO F12 Integrals
    SharedMatrix ao_f12(boost::shared_ptr<CorrelationFactor> corr);
    /// MO F12 squared Integrals
    SharedMatrix ao_f12_squared(boost::shared_ptr<CorrelationFactor> corr);
    /// MO F12G12 Integrals
    SharedMatrix ao_f12g12(boost::shared_ptr<CorrelationFactor> corr);
    /// MO F12 double commutator Integrals
    SharedMatrix ao_f12_double_commutator(boost::shared_ptr<CorrelationFactor> corr);
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
    SharedMatrix mo_f12(boost::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO F12 squared Integrals
    SharedMatrix mo_f12_squared(boost::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO F12G12 Integrals
    SharedMatrix mo_f12g12(boost::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO F12 double commutator Integrals
    SharedMatrix mo_f12_double_commutator(boost::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);

    /// Symmetric MO ERI Omega Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    SharedMatrix mo_erf_eri(double omega, SharedMatrix Cocc, SharedMatrix Cvir);

    /// AO Overlap Integrals
    SharedMatrix ao_overlap();
    /// AO Kinetic Integrals
    SharedMatrix ao_kinetic();
    /// AO Potential Integrals
    SharedMatrix ao_potential();
    /// Vector AO Dipole Integrals
    std::vector<SharedMatrix> ao_dipole();
    /// Vector AO Angular Momentum Integrals
    std::vector<SharedMatrix > ao_angular_momentum();
    /// Vector AO Nabla Integrals
    std::vector<SharedMatrix > ao_nabla();
    /// SO Overlap Integrals
    SharedMatrix so_overlap();
    /// SO Kinetic Integrals
    SharedMatrix so_kinetic();
    /// SO Potential Integrals
    SharedMatrix so_potential(bool include_perturbations = true);
    /// Vector SO Dipole Integrals
    std::vector<SharedMatrix > so_dipole();
    /// Vector SO Nabla Integrals
    std::vector<SharedMatrix > so_nabla();
    /// Vector SO Angular Momentum Integrals
    std::vector<SharedMatrix > so_angular_momentum();
    /// Vector SO Quadrupole Integrals
    std::vector<SharedMatrix > so_quadrupole();
    /// Vector SO Traceless Quadrupole Integrals
    std::vector<SharedMatrix > so_traceless_quadrupole();

    /// Returns a CdSalcList object
    boost::shared_ptr<CdSalcList> cdsalcs(int needed_irreps=0xF,
                                          bool project_out_translations=true,
                                          bool project_out_rotations=true);

    /// Play function
    void play();
};

}

#endif
