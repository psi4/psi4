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

#include "mp2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"

#include <Einsums/TensorAlgebra.hpp>
#include <Einsums/Tensor/DiskTensor.hpp>

namespace psi {
namespace f12 {

MP2F12::MP2F12(SharedWavefunction ref_wfn, Options& options) : Wavefunction(options) {
    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;

    common_init();
}

void MP2F12::common_init() {
    if (options_.get_str("REFERENCE") != "RHF") {
        throw PsiException("Only a restricted reference may be used", __FILE__, __LINE__);
    }

    options_.set_global_str("SCREENING", "NONE");

    print_ = options_.get_int("PRINT");
    singles_ = options_.get_bool("CABS_SINGLES");

    f12_type_ = options_.get_str("MP2_TYPE");
    f12_subtype_ = options_.get_str("F12_SUBTYPE");
    f12_read_ints_ = options_.get_bool("F12_READ_INTS");

    nobs_ = reference_wavefunction_->basisset()->nbf();
    nocc_ = reference_wavefunction_->doccpi()[0];
    nvir_ = nobs_ - nocc_;
    nfrzn_ = reference_wavefunction_->frzcpi()[0];
    nact_ = nocc_ - nfrzn_;
    naux_ = 0;

    if (f12_type_.find("DF") != std::string::npos) {
        use_df_ = true;
        DFBS_ = reference_wavefunction_->get_basisset("DF_BASIS_MP2");
        naux_ = DFBS_->nbf();
    }

    beta_ = options_.get_double("F12_BETA");
    cgtg_ = reference_wavefunction_->mintshelper()->f12_cgtg(beta_);

    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = Process::environment.get_n_threads();
#endif
}

void MP2F12::print_header() {
    outfile->Printf("\n -----------------------------------------------------------\n");
    if (use_df_) {
        outfile->Printf("                      DF-MP2-F12/3C(FIX)                    \n");
        outfile->Printf("             Density-Fitted Explicitly Correlated           \n");
        outfile->Printf("               2nd Order Moeller-Plesset Theory             \n");
        outfile->Printf("                RMP2 Wavefunction, %2d Threads              \n\n", nthreads_);
        outfile->Printf("                        Erica Mitchell                      \n");
    } else {
        outfile->Printf("                        MP2-F12/3C(FIX)                     \n");
        outfile->Printf("                     Explicitly Correlated                  \n");
        outfile->Printf("               2nd Order Moeller-Plesset Theory             \n");
        outfile->Printf("                RMP2 Wavefunction, %2d Threads              \n\n", nthreads_);
        outfile->Printf("                        Erica Mitchell                      \n");
    }
    outfile->Printf(" -----------------------------------------------------------\n\n");
    outfile->Printf(" Using %s %s algorithm \n\n", f12_type_.c_str(), f12_subtype_.c_str());
}

void MP2F12::form_basissets() {
    outfile->Printf(" ===> Forming the OBS and CABS <===\n\n");

    outfile->Printf("  Orbital Basis Set (OBS)\n");
    OrbitalSpace OBS = reference_wavefunction_->alpha_orbital_space("p", "SO", "ALL");
    OBS.basisset()->print();

    outfile->Printf("  Complementary Auxiliary Basis Set (CABS)\n");
    OrbitalSpace RI = OrbitalSpace::build_ri_space(reference_wavefunction_->get_basisset("CABS"), 1.0e-8);
    OrbitalSpace CABS = OrbitalSpace::build_cabs_space(OBS, RI, 1.0e-6);
    CABS.basisset()->print();

    if (use_df_) {
        outfile->Printf("  Auxiliary Basis Set\n");
        DFBS_->print();
    }

    nri_ = CABS.dim().max() + nobs_;
    ncabs_ = nri_ - nobs_;

    if (nfrzn_ != 0) {
        outfile->Printf("  Frozen Core Orbitals: %3d \n\n", nfrzn_);
    }

    outfile->Printf("  ----------------------------------------\n");
    outfile->Printf("     %5s  %5s   %5s  %5s  %5s   \n", "NOCC", "NOBS", "NCABS", "NRI", "NAUX");
    outfile->Printf("  ----------------------------------------\n");
    outfile->Printf("     %5d  %5d   %5d  %5d  %5d   \n", nocc_, nobs_, ncabs_, nri_, naux_);
    outfile->Printf("  ----------------------------------------\n");

    bs_ = {OBS, CABS};
}

void MP2F12::form_D(einsums::Tensor<double, 4>* D, einsums::Tensor<double, 2>* f) {
    using namespace einsums;

    std::vector<double> e_vir(nvir_);
    for (size_t a = 0; a < nvir_; a++) e_vir[a] = (*f)(nocc_ + a, nocc_ + a);

#pragma omp parallel for collapse(2) num_threads(nthreads_)
    for (size_t i = nfrzn_; i < nocc_; i++) {
        for (size_t j = nfrzn_; j < nocc_; j++) {
            const double e_ij = (*f)(i, i) + (*f)(j, j);
            for (size_t a = 0; a < nvir_; a++) {
                const double e_a = e_vir[a];
                for (size_t b = 0; b < nvir_; b++) {
                    (*D)(i - nfrzn_, j - nfrzn_, a, b) = 1.0 / (e_vir[b] + e_a - e_ij);
                }
            }
        }
    }
}

void MP2F12::form_f12_energy(einsums::Tensor<double, 4>* V, einsums::Tensor<double, 4>* X,
                             einsums::Tensor<double, 4>* C, einsums::Tensor<double, 4>* B,
                             einsums::Tensor<double, 2>* f, einsums::Tensor<double, 4>* G,
                             einsums::Tensor<double, 4>* D) {
    using namespace einsums;
    using namespace einsums::tensor_algebra;
    using namespace einsums::index;

    // Pre-build triangular pair list for load-balanced parallel dispatch.
    struct Pair { size_t i, j; };
    std::vector<Pair> pairs;
    pairs.reserve(nact_ * (nact_ + 1) / 2);
    for (size_t i = 0; i < nact_; i++)
        for (size_t j = i; j < nact_; j++)
            pairs.push_back({i, j});

    struct PairResult { size_t i, j; double E_s, E_t; };
    std::vector<PairResult> results(pairs.size());

    double E_f12_s = 0.0, E_f12_t = 0.0;

#pragma omp parallel for schedule(dynamic) reduction(+:E_f12_s,E_f12_t) num_threads(nthreads_)
    for (size_t p = 0; p < pairs.size(); p++) {
        const size_t i = pairs[p].i;
        const size_t j = pairs[p].j;

        // B_tilde scalars: only the two elements consumed by the energy expression.
        const double f_scale = (*f)(i + nfrzn_, i + nfrzn_) + (*f)(j + nfrzn_, j + nfrzn_);
        double B_val_ijij = (*B)(i, j, i, j) - f_scale * (*X)(i, j, i, j);
        double B_val_ijji = (*B)(i, j, j, i) - f_scale * (*X)(i, j, j, i);

        // V_tilde scalars: only the two elements consumed by the energy expression.
        double V_val_ij = (*V)(i, j, i, j);
        double V_val_ji = (*V)(i, j, j, i);

        // Apply C/D/G corrections in O(nvir^2) per pair, no rank-4 temporaries.
        {
            const auto G_ij = TensorView<double, 2>{*G, Dim<2>{nvir_, nvir_}, Offset<4>{i, j, 0, 0},
                                                     Stride<2>{(*G).stride(2), (*G).stride(3)}};
            const auto D_ij = TensorView<double, 2>{*D, Dim<2>{nvir_, nvir_}, Offset<4>{i, j, 0, 0},
                                                     Stride<2>{(*D).stride(2), (*D).stride(3)}};
            const auto C_ij = TensorView<double, 2>{*C, Dim<2>{nvir_, nvir_}, Offset<4>{i, j, 0, 0},
                                                     Stride<2>{(*C).stride(2), (*C).stride(3)}};

            // GD = G_ij .* D_ij  (for V correction)
            Tensor<double, 2> GD{"GD", nvir_, nvir_};
            einsum(Indices{a, b}, &GD, Indices{a, b}, G_ij, Indices{a, b}, D_ij);

            // CD_ij = C_ij .* D_ij  (for B correction)
            Tensor<double, 2> CD_ij{"CD_ij", nvir_, nvir_};
            einsum(Indices{a, b}, &CD_ij, Indices{a, b}, C_ij, Indices{a, b}, D_ij);

            Tensor<double, 0> tmp;

            // V(i,j,i,j) correction: sum_{a,b} C_ij(a,b) * G_ij(a,b) * D_ij(a,b)
            einsum(Indices{}, &tmp, Indices{a, b}, C_ij, Indices{a, b}, GD);
            V_val_ij -= tmp;

            // B(i,j,i,j) correction: sum_{a,b} C_ij(a,b)^2 * D_ij(a,b)
            einsum(Indices{}, &tmp, Indices{a, b}, C_ij, Indices{a, b}, CD_ij);
            B_val_ijij -= tmp;

            const auto C_ji = TensorView<double, 2>{*C, Dim<2>{nvir_, nvir_}, Offset<4>{j, i, 0, 0},
                                                     Stride<2>{(*C).stride(2), (*C).stride(3)}};

            // V(i,j,j,i) correction: sum_{a,b} C_ji(a,b) * G_ij(a,b) * D_ij(a,b)
            einsum(Indices{}, &tmp, Indices{a, b}, C_ji, Indices{a, b}, GD);
            V_val_ji -= tmp;

            // B(i,j,j,i) correction: sum_{a,b} C_ji(a,b) * C_ij(a,b) * D_ij(a,b)
            einsum(Indices{}, &tmp, Indices{a, b}, C_ji, Indices{a, b}, CD_ij);
            B_val_ijji -= tmp;
        }

        // Fixed-amplitude energy expression (matches original formula exactly).
        const int kd = (i == j) ? 1 : 2;
        const double ts = t_(i, j, i, j) + t_(i, j, j, i);
        const double V_s = 0.25 * ts * kd * (V_val_ij + V_val_ji);
        const double B_s = 0.125 * ts * kd * (B_val_ijij + B_val_ijji) * ts * kd;
        const double E_s = kd * (2.0 * V_s + B_s);
        E_f12_s += E_s;

        double E_t = 0.0;
        if (i != j) {
            const double tt = t_(i, j, i, j) - t_(i, j, j, i);
            const double V_t = 0.25 * tt * kd * (V_val_ij - V_val_ji);
            const double B_t = 0.125 * tt * kd * (B_val_ijij - B_val_ijji) * tt * kd;
            E_t = 3.0 * kd * (2.0 * V_t + B_t);
            E_f12_t += E_t;
        }

        results[p] = {i, j, E_s, E_t};
    }

    // Print pair energies in canonical (i,j) order.
    outfile->Printf("  \n");
    outfile->Printf("  %1s   %1s  |     %14s     %14s     %12s \n", "i", "j", "E_F12(Singlet)", "E_F12(Triplet)",
                    "E_F12");
    outfile->Printf(" ----------------------------------------------------------------------\n");
    for (const auto& r : results) {
        outfile->Printf("%3d %3d  |   %16.12f   %16.12f     %16.12f \n",
                        r.i + nfrzn_ + 1, r.j + nfrzn_ + 1, r.E_s, r.E_t, r.E_s + r.E_t);
    }

    set_scalar_variable("MP2-F12 OPPOSITE-SPIN CORRELATION ENERGY", E_f12_s + scalar_variable("MP2 OPPOSITE-SPIN CORRELATION ENERGY"));
    set_scalar_variable("MP2-F12 SAME-SPIN CORRELATION ENERGY", E_f12_t + scalar_variable("MP2 SAME-SPIN CORRELATION ENERGY"));

    E_f12_ = E_f12_s + E_f12_t;
}

void MP2F12::form_cabs_singles(einsums::Tensor<double, 2>* f) {
    using namespace einsums;
    using namespace einsums::linear_algebra;

    int all_vir = nvir_ + ncabs_;

    // Diagonalize f_ij and f_AB
    Tensor<double, 2> C_ij{"occupied e-vecs", nocc_, nocc_};
    Tensor<double, 2> C_AB{"vir and CABS e-vecs", all_vir, all_vir};

    Tensor<double, 1> e_ij{"occupied e-vals", nocc_};
    Tensor<double, 1> e_AB{"vir and CABS e-vals", all_vir};
    {
        C_ij = (*f)(Range{0, nocc_}, Range{0, nocc_});
        C_AB = (*f)(Range{nocc_, nri_}, Range{nocc_, nri_});

        syev(&C_ij, &e_ij);
        syev(&C_AB, &e_AB);
    }

    // Form f_iA
    Tensor<double, 2> f_iA{"Fock Occ-All_vir", nocc_, all_vir};
    {
        Tensor f_view = (*f)(Range{0, nocc_}, Range{nocc_, nri_});

        gemm<false, false>(1.0, C_ij, gemm<false, true>(1.0, f_view, C_AB), 0.0, &f_iA);
    }

    double E_s = 0.0;
#pragma omp parallel for collapse(2) num_threads(nthreads_) reduction(+ : E_s)
    for (size_t A = 0; A < all_vir; A++) {
        for (size_t i = 0; i < nocc_; i++) {
            E_s += 2 * f_iA(i, A) * f_iA(i, A) / (e_ij(i) - e_AB(A));
        }
    }

    set_scalar_variable("F12 CABS CORRECTION ENERGY", E_s);
    E_singles_ = E_s;
}

double MP2F12::compute_energy() {
    timer_on("MP2-F12 Compute Energy");
    using namespace einsums;
    einsums::profile::initialize();

    print_header();

    /* Form the orbital spaces */
    timer_on("OBS and CABS");
    form_basissets();
    timer_off("OBS and CABS");

    outfile->Printf("\n ===> Forming the Integrals <===");
    outfile->Printf("\n No screening will be used to compute integrals\n");

    /* Form the Fock Matrix */
    outfile->Printf("   Fock Matrix\n");
    auto f = std::make_unique<Tensor<double, 2>>("Fock Matrix", nri_, nri_);
    auto k = std::make_unique<Tensor<double, 2>>("Exchange MO Integral", nri_, nri_);
    timer_on("Fock Matrix");
    if (use_df_) {
        form_df_fock(f.get(), k.get());
    } else {
        form_fock(f.get(), k.get());
    }
    timer_off("Fock Matrix");

    /* Form the F12 Intermediates */
    outfile->Printf("\n ===> Forming the F12 Intermediate Tensors <===\n");
    auto V = std::make_unique<Tensor<double, 4>>("V Intermediate Tensor", nact_, nact_, nact_, nact_);
    auto X = std::make_unique<Tensor<double, 4>>("X Intermediate Tensor", nact_, nact_, nact_, nact_);
    auto C = std::make_unique<Tensor<double, 4>>("C Intermediate Tensor", nact_, nact_, nvir_, nvir_);
    auto B = std::make_unique<Tensor<double, 4>>("B Intermediate Tensor", nact_, nact_, nact_, nact_);
    auto D = std::make_unique<Tensor<double, 4>>("D Tensor", nact_, nact_, nvir_, nvir_);
    auto G_ijab = std::make_unique<Tensor<double, 4>>("ERI <ij|ab>", nact_, nact_, nvir_, nvir_);

    if (use_df_) {
        outfile->Printf("   [J_AB]^(-1)\n");
        auto J_inv_AB = std::make_unique<Tensor<double, 3>>("Metric MO ([J_AB]^{-1})", naux_, nact_, nri_);
        timer_on("Metric Integrals");
        form_metric_ints(J_inv_AB.get(), false);
        timer_off("Metric Integrals");

        outfile->Printf("   V Intermediate\n");
        outfile->Printf("   X Intermediate\n");
        timer_on("V and X Intermediate");
        form_df_V_X(V.get(), X.get(), J_inv_AB.get());
        timer_off("V and X Intermediate");

        outfile->Printf("   C Intermediate\n");
        timer_on("C Intermediate");
        form_df_C(C.get(), f.get(), J_inv_AB.get());
        timer_off("C Intermediate");

        outfile->Printf("   B Intermediate\n");
        timer_on("B Intermediate");
        form_df_B(B.get(), f.get(), k.get(), J_inv_AB.get());
        k.reset();
        timer_off("B Intermediate");

        timer_on("Energy Denom");
        form_D(D.get(), f.get());
        timer_off("Energy Denom");

        auto G = Tensor<double, 4>{"MO G Tensor", nact_, nact_, nobs_, nobs_};
        timer_on("ERI <ij|ab>");
        form_df_teints("G", &G, J_inv_AB.get(), {'o', 'O', 'o', 'O'});
        (*G_ijab) = G(Range{0, nact_}, Range{0, nact_}, Range{nocc_, nobs_}, Range{nocc_, nobs_});
        timer_off("ERI <ij|ab>");
    } else {
        outfile->Printf("   V Intermediate\n");
        outfile->Printf("   X Intermediate\n");
        timer_on("V and X Intermediate");
        form_V_X(V.get(), X.get());
        timer_off("V and X Intermediate");

        outfile->Printf("   C Intermediate\n");
        timer_on("C Intermediate");
        form_C(C.get(), f.get());
        timer_off("C Intermediate");

        outfile->Printf("   B Intermediate\n");
        timer_on("B Intermediate");
        form_B(B.get(), f.get(), k.get());
        k.reset();
        timer_off("B Intermediate");

        timer_on("Energy Denom");
        form_D(D.get(), f.get());
        timer_off("Energy Denom");

        auto G = Tensor<double, 4>{"MO G Tensor", nact_, nact_, nobs_, nobs_};
        timer_on("ERI <ij|ab>");
        form_teints("G", &G, {'o', 'O', 'o', 'O'});
        (*G_ijab) = G(All, All, Range{nocc_, nobs_}, Range{nocc_, nobs_});
        timer_off("ERI <ij|ab>");
    }

    /* Compute the MP2F12/3C Energy */
    outfile->Printf("\n ===> Computing F12/3C(FIX) Energy Correction <===\n");
    timer_on("F12 Energy Correction");
    form_f12_energy(V.get(), X.get(), C.get(), B.get(), f.get(), G_ijab.get(), D.get());
    V.reset();
    X.reset();
    C.reset();
    B.reset();
    G_ijab.reset();
    D.reset();
    timer_off("F12 Energy Correction");

    if (singles_ == true) {
        timer_on("CABS Singles Correction");
        form_cabs_singles(f.get());
        timer_off("CABS Singles Correction");
    }

    print_results();

    if (print_ > 1) {
        einsums::profile::report("timer_mp2f12.dat", false);
    }
    einsums::profile::finalize();
    timer_off("MP2-F12 Compute Energy");

    // Typically you would build a new wavefunction and populate it with data
    return E_mp2f12_;
}

void MP2F12::print_results() {
    if (use_df_) {
        outfile->Printf("\n ===> DF-MP2-F12/3C(FIX) Energies <===\n\n");
    } else {
        outfile->Printf("\n ===> MP2-F12/3C(FIX) Energies <===\n\n");
    }

    auto E_rhf = scalar_variable("CURRENT REFERENCE ENERGY");
    auto E_mp2 = scalar_variable("MP2 CORRELATION ENERGY");

    E_mp2f12_ = E_rhf + E_mp2 + E_f12_ + E_singles_;

    outfile->Printf("     RHF Reference Energy:              %16.12f \n", E_rhf);
    if (singles_ == true) {
        outfile->Printf("     CABS Singles Correction:           %16.12f \n", E_singles_);
    }
    outfile->Printf("     CABS-corrected Reference Energy:   %16.12f \n", E_rhf + E_singles_);

    outfile->Printf("     MP2 Correlation Energy:            %16.12f \n", E_mp2);
    outfile->Printf("     F12/3C(FIX) Correction Energy:     %16.12f \n", E_f12_);
    outfile->Printf("     MP2-F12/3C(FIX) Correlation Energy:%16.12f \n", E_f12_ + E_mp2);
    if (use_df_) {
        outfile->Printf("  Total DF-MP2-F12/3C(FIX) Energy:      %16.12f \n", E_mp2f12_);
    } else {
        outfile->Printf("  Total MP2-F12/3C(FIX) Energy:         %16.12f \n", E_mp2f12_);
    }

    set_scalar_variable("HF-CABS TOTAL ENERGY", E_rhf + E_singles_);
    set_scalar_variable("MP2-F12 CORRECTION ENERGY", E_f12_);
    set_scalar_variable("MP2-F12 CORRELATION ENERGY", E_mp2 + E_f12_);
    set_scalar_variable("MP2-F12 TOTAL ENERGY", E_mp2f12_);

    set_scalar_variable("MP2-F12 SINGLES ENERGY", 0.0);  // RHF
    set_scalar_variable("MP2-F12 DOUBLES ENERGY", E_mp2 + E_f12_);
}

double MP2F12::t_(int p, int q, int r, int s) {
    if (p == q && p == r && p == s) {
        return 0.5;
    } else if (p == r && q == s) {
        return 0.375;
    } else if (q == r && p == s) {
        return 0.125;
    }
    return 0.0;
}

////////////////////////////////
//* Disk Algorithm (CONV/DF) *//
////////////////////////////////

DiskMP2F12::DiskMP2F12(SharedWavefunction reference_wavefunction, Options& options)
    : MP2F12(reference_wavefunction, options) {
}

void DiskMP2F12::form_D(einsums::DiskTensor<double, 4>* D, einsums::DiskTensor<double, 2>* f) {
    using namespace einsums;

    auto f_view = (*f)(Range{0, nobs_}, Range{0, nobs_});
    f_view.set_read_only(true);

    // HDF5 concurrent writes to the same dataset are not thread-safe; the
    // O(nact^2 * nvir^2) arithmetic cost is negligible, so compute serially.
    for (size_t i = nfrzn_; i < nocc_; i++) {
        for (size_t j = nfrzn_; j < nocc_; j++) {
            auto D_view = (*D)(i - nfrzn_, j - nfrzn_, All, All);
            double e_ij = -1.0 * (f_view(i, i) + f_view(j, j));

            for (size_t a = nocc_; a < nobs_; a++) {
                for (size_t b = nocc_; b < nobs_; b++) {
                    D_view(a - nocc_, b - nocc_) = 1.0 / (e_ij + f_view(a, a) + f_view(b, b));
                }
            }
        }
    }
}

void DiskMP2F12::form_f12_energy(einsums::DiskTensor<double, 4>* V, einsums::DiskTensor<double, 4>* X,
                                 einsums::DiskTensor<double, 4>* C, einsums::DiskTensor<double, 4>* B,
                                 einsums::DiskTensor<double, 2>* f, einsums::DiskTensor<double, 4>* G,
                                 einsums::DiskTensor<double, 4>* D) {
    using namespace einsums;
    using namespace einsums::tensor_algebra;
    using namespace einsums::index;

    auto E_f12_s = 0.0;
    auto E_f12_t = 0.0;

    auto f_act = (*f)(Range{nfrzn_, nocc_}, Range{nfrzn_, nocc_});
    f_act.set_read_only(true);

    // Load B, X, V once from HDF5 — they are nact^4 and fit in memory by
    // assumption (the disk algorithm is chosen only when MO integrals fit).
    Tensor<double, 4> B_mem{"B_mem", nact_, nact_, nact_, nact_};
    Tensor<double, 4> X_mem{"X_mem", nact_, nact_, nact_, nact_};
    Tensor<double, 4> V_mem{"V_mem", nact_, nact_, nact_, nact_};
    {
        auto Bv = (*B)(All, All, All, All);
        Bv.set_read_only(true);
        B_mem = Bv.get();

        auto Xv = (*X)(All, All, All, All);
        Xv.set_read_only(true);
        X_mem = Xv.get();

        auto Vv = (*V)(All, All, All, All);
        Vv.set_read_only(true);
        V_mem = Vv.get();
    }

    outfile->Printf("  \n");
    outfile->Printf("  %1s   %1s  |     %14s     %14s     %12s \n", "i", "j", "E_F12(Singlet)", "E_F12(Triplet)",
                    "E_F12");
    outfile->Printf(" ----------------------------------------------------------------------\n");
    for (size_t i = 0; i < nact_; i++) {
        for (size_t j = i; j < nact_; j++) {
            const double f_scale = f_act(i, i) + f_act(j, j);

            // B_tilde scalars: only the two elements consumed by the energy expression.
            double B_val_ijij = B_mem(i, j, i, j) - f_scale * X_mem(i, j, i, j);
            double B_val_ijji = B_mem(i, j, j, i) - f_scale * X_mem(i, j, j, i);

            // V_tilde scalars.
            double V_val_ij = V_mem(i, j, i, j);
            double V_val_ji = V_mem(i, j, j, i);

            // Read G and D for this pair from disk (nvir^2 each).
            Tensor<double, 2> G_ij{"G_ij", nvir_, nvir_};
            {
                auto Gv = (*G)(i, j, Range{nocc_, nobs_}, Range{nocc_, nobs_});
                Gv.set_read_only(true);
                G_ij = Gv.get();
            }
            Tensor<double, 2> D_ij{"D_ij", nvir_, nvir_};
            {
                auto Dv = (*D)(i, j, All, All);
                Dv.set_read_only(true);
                D_ij = Dv.get();
            }

            // Apply C corrections reading only C(i,j) and C(j,i) from disk.
            {
                // GD = G_ij .* D_ij  (for V correction)
                Tensor<double, 2> GD{"GD", nvir_, nvir_};
                einsum(Indices{a, b}, &GD, Indices{a, b}, G_ij, Indices{a, b}, D_ij);

                auto C_ij_view = (*C)(i, j, All, All);
                C_ij_view.set_read_only(true);
                Tensor<double, 2> C_ij = C_ij_view.get();

                // CD_ij = C_ij .* D_ij  (for B correction)
                Tensor<double, 2> CD_ij{"CD_ij", nvir_, nvir_};
                einsum(Indices{a, b}, &CD_ij, Indices{a, b}, C_ij, Indices{a, b}, D_ij);

                Tensor<double, 0> tmp;

                einsum(Indices{}, &tmp, Indices{a, b}, C_ij, Indices{a, b}, GD);
                V_val_ij -= tmp;

                einsum(Indices{}, &tmp, Indices{a, b}, C_ij, Indices{a, b}, CD_ij);
                B_val_ijij -= tmp;

                auto C_ji_view = (*C)(j, i, All, All);
                C_ji_view.set_read_only(true);
                Tensor<double, 2> C_ji = C_ji_view.get();

                einsum(Indices{}, &tmp, Indices{a, b}, C_ji, Indices{a, b}, GD);
                V_val_ji -= tmp;

                einsum(Indices{}, &tmp, Indices{a, b}, C_ji, Indices{a, b}, CD_ij);
                B_val_ijji -= tmp;
            }

            // Fixed-amplitude energy expression (matches original formula exactly).
            const int kd = (i == j) ? 1 : 2;
            const double ts = t_(i, j, i, j) + t_(i, j, j, i);
            const double V_s = 0.25 * ts * kd * (V_val_ij + V_val_ji);
            const double B_s = 0.125 * ts * kd * (B_val_ijij + B_val_ijji) * ts * kd;
            const double E_s = kd * (2.0 * V_s + B_s);
            E_f12_s += E_s;

            double E_t = 0.0;
            if (i != j) {
                const double tt = t_(i, j, i, j) - t_(i, j, j, i);
                const double V_t = 0.25 * tt * kd * (V_val_ij - V_val_ji);
                const double B_t = 0.125 * tt * kd * (B_val_ijij - B_val_ijji) * tt * kd;
                E_t = 3.0 * kd * (2.0 * V_t + B_t);
                E_f12_t += E_t;
            }
            outfile->Printf("%3d %3d  |   %16.12f   %16.12f     %16.12f \n", i + 1, j + 1, E_s, E_t, E_s + E_t);
        }
    }

    set_scalar_variable("MP2-F12 OPPOSITE-SPIN CORRELATION ENERGY", E_f12_s + scalar_variable("MP2 OPPOSITE-SPIN CORRELATION ENERGY"));
    set_scalar_variable("MP2-F12 SAME-SPIN CORRELATION ENERGY", E_f12_t + scalar_variable("MP2 SAME-SPIN CORRELATION ENERGY"));

    E_f12_ = E_f12_s + E_f12_t;
}

void DiskMP2F12::form_cabs_singles(einsums::DiskTensor<double, 2>* f) {
    using namespace einsums;
    using namespace einsums::linear_algebra;

    int all_vir = nvir_ + ncabs_;

    // Diagonalize f_ij and f_AB
    Tensor<double, 1> e_ij{"occupied e-vals", nocc_};
    Tensor<double, 1> e_AB{"vir and CABS e-vals", all_vir};

    auto C_ij = (*f)(Range{0, nocc_}, Range{0, nocc_});
    auto C_AB = (*f)(Range{nocc_, nri_}, Range{nocc_, nri_});

    syev(&(C_ij.get()), &e_ij);
    syev(&(C_AB.get()), &e_AB);

    // Form f_iA
    Tensor<double, 2> f_iA{"Fock Occ-All_vir", nocc_, all_vir};
    {
        auto f_view = (*f)(Range{0, nocc_}, Range{nocc_, nri_});
        f_view.set_read_only(true);

        gemm<false, false>(1.0, C_ij.get(), gemm<false, true>(1.0, f_view.get(), C_AB.get()), 0.0, &f_iA);
    }

    double E_s = 0.0;
#pragma omp parallel for collapse(2) num_threads(nthreads_) reduction(+ : E_s)
    for (size_t A = 0; A < all_vir; A++) {
        for (size_t i = 0; i < nocc_; i++) {
            E_s += 2 * f_iA(i, A) * f_iA(i, A) / (e_ij(i) - e_AB(A));
        }
    }

    set_scalar_variable("F12 CABS CORRECTION ENERGY", E_s);
    E_singles_ = E_s;
}

double DiskMP2F12::compute_energy() {
    timer_on("MP2-F12 Compute Energy");
    using namespace einsums;
    einsums::profile::initialize();

    print_header();

    /* Form the orbital spaces */
    timer_on("OBS and CABS");
    form_basissets();
    timer_off("OBS and CABS");

    // Disable HDF5 diagnostic reporting.
    H5Eset_auto(0, nullptr, nullptr);
    std::string file_name = "Data_" + std::to_string(nocc_) + "_" + std::to_string(ncabs_);
    if (use_df_) {
        file_name += "_" + std::to_string(naux_);
    }
    file_name += ".h5";

    if (f12_read_ints_) {
        // Reads existing file
        // formerly lhs einsums::state::data() . this is likely to change again in Einsums in the near future
        ein_state_data_ = h5::open(file_name, H5F_ACC_RDWR);
    } else {
        // Creates new file
        ein_state_data_ = h5::create(file_name, H5F_ACC_TRUNC);
    }

    outfile->Printf("\n ===> Forming the Integrals <===");
    outfile->Printf("\n No screening will be used to compute integrals\n\n");

    /* Form the two-electron integrals */
    // Two-Electron Integrals
    auto G = std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "MO G Tensor", nact_, nact_, nobs_, nri_);
    auto F = std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "MO F12 Tensor", nact_, nact_, nri_, nri_);
    auto F2 =
        std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "MO F12_Squared Tensor", nact_, nact_, nact_, nri_);
    auto FG = std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "MO F12G12 Tensor", nact_, nact_, nact_, nact_);
    auto Uf = std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "MO F12_DoubleCommutator Tensor", nact_, nact_,
                                                      nact_, nact_);

    std::vector<std::string> teint = {};
    if (!(*FG).existed()) teint.push_back("FG");
    if (!(*Uf).existed()) teint.push_back("Uf");
    if (!(*G).existed()) teint.push_back("G");
    if (!(*F).existed()) teint.push_back("F");
    if (!(*F2).existed()) teint.push_back("F2");
    if (teint.size() == 0) outfile->Printf("   Two-Electron Integrals\n");

    // Fock Matrices
    auto f = std::make_unique<DiskTensor<double, 2>>(ein_state_data_, "Fock Matrix", nri_, nri_);
    auto k = std::make_unique<DiskTensor<double, 2>>(ein_state_data_, "Exchange Matrix", nri_, nri_);
    auto fk = std::make_unique<DiskTensor<double, 2>>(ein_state_data_, "Fock-Exchange Matrix", nri_, nri_);

    if (use_df_) {
        outfile->Printf("   Fock Matrix\n");
        if (!(*f).existed() && !(*k).existed() && !(*fk).existed()) {
            timer_on("Fock Matrix");
            form_df_fock(f.get(), k.get(), fk.get());
            timer_off("Fock Matrix");
        }

        timer_on("Metric Integrals");
        auto J_inv_AB = std::make_unique<Tensor<double, 3>>("Metric MO ([J_AB]^{-1})", naux_, nact_, nri_);
        form_metric_ints(J_inv_AB.get(), false);
        timer_off("Metric Integrals");

        for (int i = 0; i < teint.size(); i++) {
            if (teint[i] == "F") {
                outfile->Printf("   F Integral\n");
                timer_on("F_12 Integral");
                form_df_teints(teint[i], F.get(), J_inv_AB.get());
                timer_off("F_12 Integral");
            } else if (teint[i] == "FG") {
                outfile->Printf("   FG Integral\n");
                timer_on("FG_12 Integral");
                form_df_teints(teint[i], FG.get(), J_inv_AB.get());
                timer_off("FG_12 Integral");
            } else if (teint[i] == "F2") {
                outfile->Printf("   F Squared Integral\n");
                timer_on("F^2_12 Integral");
                form_df_teints(teint[i], F2.get(), J_inv_AB.get());
                timer_off("F^2_12 Integral");
            } else if (teint[i] == "Uf") {
                outfile->Printf("   F Double Commutator Integral\n");
                timer_on("U^F_12 Integral");
                form_df_teints(teint[i], Uf.get(), J_inv_AB.get());
                timer_off("U^F_12 Integral");
            } else {
                outfile->Printf("   G Integral\n");
                timer_on("G Integral");
                form_df_teints(teint[i], G.get(), J_inv_AB.get());
                timer_off("G Integral");
            }
        }
    } else {
        outfile->Printf("   Fock Matrix\n");
        if (!(*f).existed() && !(*k).existed() && !(*fk).existed()) {
            timer_on("Fock Matrix");
            form_fock(f.get(), k.get(), fk.get());
            timer_off("Fock Matrix");
        }

        for (int i = 0; i < teint.size(); i++) {
            if (teint[i] == "F") {
                outfile->Printf("   F Integral\n");
                timer_on("F_12 Integral");
                form_teints(teint[i], F.get());
                timer_off("F_12 Integral");
            } else if (teint[i] == "FG") {
                outfile->Printf("   FG Integral\n");
                timer_on("FG_12 Integral");
                form_teints(teint[i], FG.get());
                timer_off("FG_12 Integral");
            } else if (teint[i] == "F2") {
                outfile->Printf("   F Squared Integral\n");
                timer_on("F^2_12 Integral");
                form_teints(teint[i], F2.get());
                timer_off("F^2_12 Integral");
            } else if (teint[i] == "Uf") {
                outfile->Printf("   F Double Commutator Integral\n");
                timer_on("U^F_12 Integral");
                form_teints(teint[i], Uf.get());
                timer_off("U^F_12 Integral");
            } else {
                outfile->Printf("   G Integral\n");
                timer_on("G Integral");
                form_teints(teint[i], G.get());
                timer_off("G Integral");
            }
        }
    }

    /* Form the F12 Matrices */
    outfile->Printf("\n ===> Forming the F12 Intermediate Tensors <===\n");
    auto V =
        std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "V Intermediate Tensor", nact_, nact_, nact_, nact_);
    auto X =
        std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "X Intermediate Tensor", nact_, nact_, nact_, nact_);
    auto C =
        std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "C Intermediate Tensor", nact_, nact_, nvir_, nvir_);
    auto B =
        std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "B Intermediate Tensor", nact_, nact_, nact_, nact_);
    auto D = std::make_unique<DiskTensor<double, 4>>(ein_state_data_, "D Tensor", nact_, nact_, nvir_, nvir_);

    outfile->Printf("   V Intermediate\n");
    if (!(*V).existed()) {
        timer_on("V Intermediate");
        form_V_X(V.get(), F.get(), G.get(), FG.get());
        timer_off("V Intermediate");
    }

    outfile->Printf("   X Intermediate\n");
    if (!(*X).existed()) {
        timer_on("X Intermediate");
        form_V_X(X.get(), F.get(), F.get(), F2.get());
        timer_off("X Intermediate");
    }

    outfile->Printf("   C Intermediate\n");
    if (!(*C).existed()) {
        timer_on("C Intermediate");
        form_C(C.get(), F.get(), f.get());
        timer_off("C Intermediate");
    }

    outfile->Printf("   B Intermediate\n");
    if (!(*B).existed()) {
        timer_on("B Intermediate");
        form_B(B.get(), Uf.get(), F2.get(), F.get(), f.get(), fk.get(), k.get());
        timer_off("B Intermediate");
    }

    if (!(*D).existed()) {
        timer_on("Energy Denom");
        form_D(D.get(), f.get());
        timer_off("Energy Denom");
    }

    /* Compute the MP2F12/3C Energy */
    outfile->Printf("\n ===> Computing F12/3C(FIX) Energy Correction <===\n");
    timer_on("F12 Energy Correction");
    form_f12_energy(V.get(), X.get(), C.get(), B.get(), f.get(), G.get(), D.get());
    timer_off("F12 Energy Correction");

    if (singles_ == true) {
        timer_on("CABS Singles Correction");
        form_cabs_singles(f.get());
        timer_off("CABS Singles Correction");
    }

    print_results();

    if (print_ > 1) {
        einsums::profile::report("timer_mp2f12.dat", false);
    }
    einsums::profile::finalize();
    timer_off("MP2-F12 Compute Energy");

    // Typically you would build a new wavefunction and populate it with data
    return E_mp2f12_;
}

}  // namespace f12
}  // namespace psi
