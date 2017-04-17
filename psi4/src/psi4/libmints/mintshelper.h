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

#ifndef _psi_src_lib_libmints_helper_h
#define _psi_src_lib_libmints_helper_h

#include <vector>
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/multipolesymmetry.h"
#include "psi4/libmints/ecp.h"

namespace psi {

class Options;
class CdSalcList;
class CorrelationFactor;
class TwoBodyAOInt;
class PetiteList;
class ThreeCenterOverlapInt;
class OneBodyAOInt;
class OneBodySOInt;

/**
* The MintsHelper object, places molecular integrals
* (and later derivative integrals) on disk
**/
class MintsHelper {

private:
    /// The Options reference for basis sets and things
    Options& options_;
    std::shared_ptr<PSIO> psio_;
    std::shared_ptr<MatrixFactory> factory_;
    std::shared_ptr<Molecule> molecule_;
    std::shared_ptr<IntegralFactory> integral_;
    std::shared_ptr<BasisSet> basisset_;
    std::shared_ptr<SOBasisSet> sobasis_;
    std::shared_ptr<BasisSet> ecpbasis_;
    std::shared_ptr<TwoBodyAOInt> eriInts_;
    std::shared_ptr<BasisSet> rel_basisset_;
    int print_;
    int nthread_;

    /// Value which any two-electron integral is below is discarded
    double cutoff_;

    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix Co, SharedMatrix Cv);
    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2,
                                                 SharedMatrix C3, SharedMatrix C4);
    /// In-core builds spin eri's
    SharedMatrix mo_spin_eri_helper(SharedMatrix Iso, int n1, int n2);


    SharedMatrix ao_helper(const std::string& label, std::shared_ptr<TwoBodyAOInt> ints);
    SharedMatrix ao_shell_getter(const std::string& label, std::shared_ptr<TwoBodyAOInt> ints, int M, int N, int P, int Q);

    SharedMatrix ao_3coverlap_helper(const std::string &label, std::shared_ptr<ThreeCenterOverlapInt> ints);

    void common_init();

    void one_body_ao_computer(std::vector<std::shared_ptr<OneBodyAOInt>> obv, SharedMatrix out, bool symm);

public:

    void init_helper(std::shared_ptr<Wavefunction> wavefunction = std::shared_ptr<Wavefunction>());
    void init_helper(std::shared_ptr<BasisSet> basis, std::shared_ptr<BasisSet> ecpbasis = nullptr);

    /// Constructor, using basisset
    MintsHelper(std::shared_ptr<BasisSet> basis,
                Options& options = Process::environment.options,
                int print = 0,  std::shared_ptr<BasisSet> ecpbasis = nullptr);

    /// Constructor, using wavefunction
    MintsHelper(std::shared_ptr<Wavefunction> wavefunction);
    /// Destructor, does nothing
    ~MintsHelper();

    OperatorSymmetry operator_symmetry(int order) {
        return OperatorSymmetry(order, molecule_, integral_, factory_);
    }

    /// Returns the number of basis functions
    int nbf() const;

    /// Sets the print level
    void set_print(int print) {print_ = print; }
    void set_nthread(int nthread) {nthread_ = nthread; }

    /// Returns petite list that is capable of transforming basis functions (nbf) to SO's.
    std::shared_ptr<PetiteList> petite_list() const;

    enum {
        kFromCartesianAO = true,
        kFromBF = false
    };
    /** Returns petite list that is capable of transforming AO basis functions (nbf) to SO's.
     *  \param include_pure_transform Is either kFromCartesianAO or kFromBF.
     */
    std::shared_ptr<PetiteList> petite_list(bool include_pure_transform) const;
    /// Basis set being used.
    std::shared_ptr<BasisSet> basisset() const;
    /// SO basis set being used.
    std::shared_ptr<SOBasisSet> sobasisset() const;
    /// ECP basis set being used.
    std::shared_ptr<BasisSet> ecpbasisset() const;
    /// Matrix factory being used
    std::shared_ptr<MatrixFactory> factory() const;
    /// Integral factory being used
    std::shared_ptr<IntegralFactory> integral() const;

    void set_rel_basisset(std::shared_ptr<BasisSet> rel_basis) { rel_basisset_ = rel_basis; }

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
    SharedMatrix ao_eri(std::shared_ptr<BasisSet> bs1,
                        std::shared_ptr<BasisSet> bs2,
                        std::shared_ptr<BasisSet> bs3,
                        std::shared_ptr<BasisSet> bs4);
    /// AO ERI Shell
    SharedMatrix ao_eri_shell(int M, int N, int P, int Q);
    /// AO ERF Integrals
    SharedMatrix ao_erf_eri(double omega);
    /// MO ERFC Omega Integrals
    SharedMatrix ao_erfc_eri(double omega);
    /// MO F12 Integrals
    SharedMatrix ao_f12(std::shared_ptr<CorrelationFactor> corr);
    SharedMatrix ao_f12(std::shared_ptr<CorrelationFactor> corr,
                        std::shared_ptr<BasisSet> bs1,
                        std::shared_ptr<BasisSet> bs2,
                        std::shared_ptr<BasisSet> bs3,
                        std::shared_ptr<BasisSet> bs4);
    /// MO F12 Integrals
    SharedMatrix ao_f12_scaled(std::shared_ptr<CorrelationFactor> corr);
    SharedMatrix ao_f12_scaled(std::shared_ptr<CorrelationFactor> corr,
                               std::shared_ptr<BasisSet> bs1,
                               std::shared_ptr<BasisSet> bs2,
                               std::shared_ptr<BasisSet> bs3,
                               std::shared_ptr<BasisSet> bs4);
    /// MO F12 squared Integrals
    SharedMatrix ao_f12_squared(std::shared_ptr<CorrelationFactor> corr);
    SharedMatrix ao_f12_squared(std::shared_ptr<CorrelationFactor> corr,
                                std::shared_ptr<BasisSet> bs1,
                                std::shared_ptr<BasisSet> bs2,
                                std::shared_ptr<BasisSet> bs3,
                                std::shared_ptr<BasisSet> bs4);
    /// MO F12G12 Integrals
    SharedMatrix ao_f12g12(std::shared_ptr<CorrelationFactor> corr);
    /// MO F12 double commutator Integrals
    SharedMatrix ao_f12_double_commutator(std::shared_ptr<CorrelationFactor> corr);
    /// 3Center overlap integrals
    SharedMatrix ao_3coverlap();
    SharedMatrix ao_3coverlap(std::shared_ptr<BasisSet> bs1,
                              std::shared_ptr<BasisSet> bs2,
                              std::shared_ptr<BasisSet> bs3);

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
    SharedMatrix mo_f12(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO F12 squared Integrals
    SharedMatrix mo_f12_squared(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO F12G12 Integrals
    SharedMatrix mo_f12g12(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);
    /// MO F12 double commutator Integrals
    SharedMatrix mo_f12_double_commutator(std::shared_ptr<CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4);

    /// Symmetric MO ERI Omega Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    SharedMatrix mo_erf_eri(double omega, SharedMatrix Cocc, SharedMatrix Cvir);

    /// Symmetric MO Spin ERI Integrals, <oo|vv> type (Full matrix (16x larger than MO ERI), N^5,
    /// most definitely not recommended for large systems)
    /// Pass C_ C_ for <aa|aa> type, Cocc_, Cocc_ for <oo|oo> type, or Cvir_, Cvir_ for <vv|vv> type
    SharedMatrix mo_spin_eri(SharedMatrix Co, SharedMatrix Cv);

    /// AO Overlap Integrals
    SharedMatrix ao_overlap();
    // JWM 4/3/2015
    SharedMatrix ao_overlap(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    /// AO Kinetic Integrals
    SharedMatrix ao_kinetic();
    SharedMatrix ao_kinetic(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    /// AO Potential Integrals
    SharedMatrix ao_potential();
    SharedMatrix ao_potential(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    /// AO ECP Integrals
    SharedMatrix ao_ecp();
    SharedMatrix ao_ecp(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    /// AO pVp Integrals
    SharedMatrix ao_pvp();
    /// AO DKH Integrals
    SharedMatrix ao_dkh(int dkh_order = -1);
    /// SO DKH Integrals
    SharedMatrix so_dkh(int dkh_order = -1);
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
    std::shared_ptr<CdSalcList> cdsalcs(int needed_irreps=0xF,
                                          bool project_out_translations=true,
                                          bool project_out_rotations=true);

    /// N^5 ao->mo transform, in memory, smart indexing
    SharedMatrix mo_transform(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2,
                                                SharedMatrix C3, SharedMatrix C4);
    /// Play function
    void play();
};

}

#endif
