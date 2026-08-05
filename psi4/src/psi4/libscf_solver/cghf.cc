/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "cghf.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/complexmatrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libfock/ComplexJK.h"
#include "psi4/libfock/v.h"

#ifdef USING_OpenOrbitalOptimizer
#include <openorbitaloptimizer/scfsolver.hpp>
#include <any>
#endif

// TODO: REMOVE THESE
#include <Einsums/LinearAlgebra.hpp>
#include <Einsums/TensorAlgebra.hpp>

namespace {

// Takes an nsopi_-shaped square SharedMatrix and copies to 2 (two) diagonal
// blocks **per irrep** into each tile of the provided ComplexMatrix.
void copy_matrix_to_complex(const psi::Matrix& A, psi::ComplexMatrix& B) {
    const int nirrep = A.nirrep();
    psi::Dimension row_dim(nirrep);
    psi::Dimension col_dim(nirrep);

    for (int h = 0; h < nirrep; h++) {
        row_dim[h] = A.rowspi(h) * 2;
        col_dim[h] = A.colspi(h) * 2;
    }

    B = psi::ComplexMatrix{B.name(), row_dim, col_dim};
    B.zero();

    for (int h = 0; h < nirrep; h++) {
        const int r_ = A.rowspi(h);
        const int c_ = A.colspi(h);
        for (int i = 0; i < r_; i++) {
            for (int j = 0; j < c_; j++) {
				B.set(h, i, j, A(h, i, j));
				B.set(h, i + r_, j + c_, A(h, i, j));
            }
        }
    }
}

}

namespace psi {
namespace scf {

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func)
    : CGHF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {}

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : ComplexWavefunction(options), BaseHF(func) {
    // Analogous to HF::HF's shallow_copy(ref_wfn); full ComplexWavefunction::shallow_copy comes later.
    molecule_ = ref_wfn->molecule();
    basisset_ = ref_wfn->basisset();
    psio_ = psio;
    // Options-only base ctor skipped init; run it now that molecule/basis are set.
    ComplexWavefunction::common_init();
    common_init();
}

CGHF::~CGHF() {}

void CGHF::common_init() {
    BaseHF::common_init(options_, module_, molecule_, dipole_field_strength_);
    name_ = "CGHF";
    // Prefer BaseHF::nelectron_ over BaseWavefunction's; copy from ComplexWavefunction.
    nelectron_ = nelec_;

    if (molecule_->schoenflies_symbol() != "c1" || nirrep_ != 1) {
        throw PSIEXCEPTION("CGHF currently supports only C1 symmetry. Set symmetry c1 in the molecule block.");
    }

    // DFT stuff (would typically go in subclass_init)
    setup_potential();

    Dimension irrep_dim(nirrep_);
    for (int h = 0; h < nirrep_; h++) {
        irrep_dim[h] = nsopi_[h] * 2;
    }

    T_ = std::make_shared<ComplexMatrix>("T", irrep_dim);
    V_ = std::make_shared<ComplexMatrix>("V", irrep_dim);
    H_ = std::make_shared<ComplexMatrix>("H", irrep_dim);
    D_ = std::make_shared<ComplexMatrix>("SCF density", irrep_dim);
    F_ = std::make_shared<ComplexMatrix>("F", irrep_dim);
    G_ = std::make_shared<ComplexMatrix>("G", irrep_dim);
    J_ = std::make_shared<ComplexMatrix>("J", irrep_dim);
    K_ = std::make_shared<ComplexMatrix>("K", irrep_dim);

    T_->zero();
    V_->zero();
    H_->zero();
    D_->zero();
    F_->zero();
    G_->zero();
    J_->zero();
    K_->zero();

    // We don't know the sizes of these until nmopi_ fills in form_Shalf();
    S_ = std::make_shared<ComplexMatrix>(); S_->set_name("Overlap");
    X_ = std::make_shared<ComplexMatrix>(); X_->set_name("Orthogonalization");
    // C_ is resized in form_C once X_ is known
    C_ = std::make_shared<ComplexMatrix>(); C_->set_name("MO coefficients");

    // How much stuff shall we echo to the user?
    if (options_["PRINT"].has_changed()) print_ = options_.get_int("PRINT");

    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("         ---------------------------------------------------------\n");
        outfile->Printf("                                COMPLEX SCF\n");
        outfile->Printf("                            by Nathan Gillispie\n");
        outfile->Printf("                              %4s Reference\n", options_.get_str("REFERENCE").c_str());
        outfile->Printf("                           %6ld MiB Core\n", memory_ / 1048576L);
        outfile->Printf("         ---------------------------------------------------------\n");
        outfile->Printf("\n");
        outfile->Printf("  ==> Geometry <==\n\n");

        molecule_->print();

        outfile->Printf("  Running in %s symmetry.\n\n", molecule_->point_group()->symbol().c_str());

        molecule_->print_rotational_constants();

        outfile->Printf("  Nuclear repulsion = %20.15f\n\n", nuclearrep_);
        outfile->Printf("  Charge       = %d\n", charge_);
        outfile->Printf("  Multiplicity = %d\n", multiplicity_);
        outfile->Printf("  Electrons    = %d\n", nelectron_);

        outfile->Printf("  ==> Algorithm <==\n\n");
        outfile->Printf("  SCF Algorithm Type is %s.\n", options_.get_str("SCF_TYPE").c_str());
        outfile->Printf("  DIIS %s.\n", options_.get_bool("DIIS") ? "enabled" : "disabled");
        if ((options_.get_int("MOM_START") != 0) && (options_["MOM_OCC"].size() != 0))  // TROUBLE, NOT SET YET?
            outfile->Printf("  Excited-state MOM enabled.\n");
        else
            outfile->Printf("  MOM %s.\n", (options_.get_int("MOM_START") == 0) ? "disabled" : "enabled");
        outfile->Printf("  Fractional occupation %s.\n", (options_.get_int("FRAC_START") == 0) ? "disabled" : "enabled");
        outfile->Printf("  Guess Type is %s.\n", options_.get_str("GUESS").c_str());
        outfile->Printf("  Energy threshold   = %3.2e\n", options_.get_double("E_CONVERGENCE"));
        outfile->Printf("  Density threshold  = %3.2e\n", options_.get_double("D_CONVERGENCE"));
        outfile->Printf("  Integral threshold = %3.2e\n\n", integral_threshold_);

        outfile->Printf("  ==> Primary Basis <==\n\n");

        basisset_->print_by_level("outfile", print_);
    }
}

void CGHF::print_orbitals() {
    std::vector<std::string> labels = molecule_->irrep_labels();

    outfile->Printf("    Orbital Energies [Eh]\n    ---------------------\n\n");

    std::vector<std::pair<double, std::pair<std::string, int>>> occ;
    std::vector<std::pair<double, std::pair<std::string, int>>> vir;

    for (int h = 0; h < nirrep_; h++) {
        std::vector<std::pair<double, int> > orb_e;
        for (int a = 0; a < nmopi_[h]; a++) orb_e.push_back(std::make_pair(epsilon_->get(h, a), a));
        std::sort(orb_e.begin(), orb_e.end());

        std::vector<int> orb_order(nmopi_[h]);
        for (int a = 0; a < nmopi_[h]; a++) orb_order[orb_e[a].second] = a;

        for (int a = 0; a < nelecpi_[h]; a++)
            occ.push_back(std::make_pair(epsilon_->get(h, a), std::make_pair(labels[h], orb_order[a] + 1)));
        for (int a = nelecpi_[h]; a < nmopi_[h]; a++)
            vir.push_back(std::make_pair(epsilon_->get(h, a), std::make_pair(labels[h], orb_order[a] + 1)));
    }
    std::sort(occ.begin(), occ.end());
    std::sort(vir.begin(), vir.end());

    outfile->Printf("    Occupied:\n\n    ");
    int count = 0;
    for (int i = 0; i < occ.size(); i++) {
        outfile->Printf("%4d%-4s%11.6f  ", occ[i].second.second, occ[i].second.first.c_str(), occ[i].first);
        if (count++ % 3 == 2 && count != occ.size()) outfile->Printf("\n    ");
    }
    outfile->Printf("\n\n    Virtual:\n\n    ");
    count = 0;
    for (int i = 0; i < vir.size(); i++) {
        outfile->Printf("%4d%-4s%11.6f  ", vir[i].second.second, vir[i].second.first.c_str(), vir[i].first);
        if (count++ % 3 == 2 && count != vir.size()) outfile->Printf("\n    ");
    }

    outfile->Printf("\n\n    Final Occupation by Irrep:\n          ");
    for (int h = 0; h < nirrep_; ++h) outfile->Printf(" %4s ", labels[h].c_str());
    outfile->Printf("\n");

    // TODO: print something like amount of alpha beta per electron? Or whatever
    // is most physically meaningful.

    outfile->Printf("\n");
}

void CGHF::setup_potential() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("DFT functionals are not supported with CGHF.");
    } else {
        potential_ = nullptr;
    }
}

void CGHF::finalize() {
    // Clean memory off, handle diis closeout, etc

    // This will be the only one
    if (!options_.get_bool("SAVE_JK")) {
        jk_.reset();
    }

    // Clean up after DIIS
    if (initialized_diis_manager_) diis_manager_.attr("delete_diis_file")();
    diis_manager_ = py::none();
    initialized_diis_manager_ = false;

    // No frozen virtual and frozen core irrep vars to set

    energy_ = energies_["Total Energy"];

    X_.reset();
    T_.reset();
}

void CGHF::set_jk(std::shared_ptr<ComplexJK> jk) {
    // Cheap basis check
    int jk_nbf = jk->basisset()->nbf();
    int hf_nbf = basisset_->nbf();
    if (hf_nbf != jk_nbf) {
        throw PSIEXCEPTION("Tried setting a JK object whos number of basis functions does not match HF's!");
    }

    jk_ = jk;
}

void CGHF::form_H() {
    SharedMatrix T_real = mintshelper()->so_kinetic();
    SharedMatrix V_real = mintshelper()->so_potential();

    copy_matrix_to_complex(*T_real, *T_);
    copy_matrix_to_complex(*V_real, *V_);

    if (debug_ > 2) T_->print("outfile");
    if (debug_ > 2) V_->print("outfile");

    if (perturb_h_)
        throw PSIEXCEPTION("Perturbed Hamiltonian not supported with CGHF.");

    if (external_pot_)
        throw PSIEXCEPTION("External potential not supported with CGHF.");

    H_->zero();
    H_->add(*T_);
    H_->add(*V_);

    if (print_ > 3) H_->print("outfile");
}

void CGHF::form_Shalf() {
    BasisSetOrthogonalization::OrthogonalizationMethod method;
    if (options_.get_str("S_ORTHOGONALIZATION") == "SYMMETRIC")
        method = BasisSetOrthogonalization::Symmetric;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "CANONICAL")
        method = BasisSetOrthogonalization::Canonical;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "PARTIALCHOLESKY")
        method = BasisSetOrthogonalization::PartialCholesky;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "AUTO")
        method = BasisSetOrthogonalization::Automatic;
    else
        throw PSIEXCEPTION("Unrecognized S_ORTHOGONALIZATION method\n");


    double lindep_tolerance = options_.get_double("S_TOLERANCE");
    double cholesky_tolerance = options_.get_double("S_CHOLESKY_TOLERANCE");

    SharedMatrix S_temp = mintshelper()->so_overlap();
    BasisSetOrthogonalization orthog(method, S_temp, lindep_tolerance, cholesky_tolerance, print_);

    // Transform
    SharedMatrix X_temp = orthog.basis_to_orthog_basis();

    // Update nmo_
    nmopi_ = X_temp->colspi();

    // Double for spin blocks
    for (auto& h : nmopi_) { h *= 2; }
    nmo_ = nmopi_.sum();

    // This check may not apply in most cases.
    for (int h = 0; h < X_temp->nirrep(); ++h) {
        if (nelecpi_[h] > nmopi_[h]) {
            throw PSIEXCEPTION("Not enough molecular orbitals to satisfy requested occupancies");
        }
    }

    // Create the correct sized quantities now.
    copy_matrix_to_complex(*S_temp, *S_);
    copy_matrix_to_complex(*X_temp, *X_);
}

void CGHF::guess() {
    double guess_E;
    std::string guess_type = options_.get_str("GUESS");

    if (guess_type == "CORE") {
        if (print_) outfile->Printf("  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");

        (*F_) = (*H_); // Try the core Hamiltonian as the Fock Matrix
        form_initial_C(); // calls only CGHF::form_C()
        form_D();
        guess_E = compute_initial_E();
    } else {
        throw PSIEXCEPTION("CGHF '" + guess_type + "' GUESS not implemented. Use 'CORE'.");
    }

    energies_["Total Energy"] = 0.0;  // don't use this guess in our convergence checks
}

double CGHF::compute_initial_E() {
    /*  \sum_{ij}h_{ij} \gamma_{ji} */
    std::complex<double> E = H_->product_trace(*D_);

    if (E.imag() > 1e-12) {
        outfile->Printf("WARNING: CGHF::compute_initial_E found large imaginary %+fi\n"
                        "  Is D Hermitian?\n", E.imag());
    }

    return nuclearrep_ + E.real();
}

double CGHF::compute_E() {
    std::complex<double> one_electron_E = H_->product_trace(*D_);
    std::complex<double> kinetic_E = T_->product_trace(*D_);
    std::complex<double> coulomb_E = J_->product_trace(*D_);
    std::complex<double> exchange_E = K_->product_trace(*D_);

    double two_electron_E = 0.5 * (coulomb_E.real() - exchange_E.real());

    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = kinetic_E.real();
    energies_["One-Electron"] = one_electron_E.real();
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"] = 0.0;
    energies_["VV10"] = 0.0;
    energies_["-D"] = 0.0;

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E.real();
    Etotal += two_electron_E;

    return Etotal;
}


void CGHF::form_C(double shift) {
    if (shift != 0.0) throw PSIEXCEPTION("Level shifting not available for CGHF.");

    // Form F' = X'FX for canonical orthogonalization
    auto Forth = linalg::triplet<true, false, false>(X_, F_, X_);
    Forth->set_name("Orthogonalized Fock");

    // Form C' = eig(F')
    epsilon_ = std::make_shared<Vector>("Orbital energies", nmopi_);

    for (int h = 0; h < nirrep_; h++) {
        // Do not diagonalize 0x0 matrix
        if (nmopi_[h] == 0) continue;

        auto evals = einsums::Tensor<double, 1>("Fock evals", nmopi_[h]);
        evals.zero();

        // Hermitian eigensolver one block at a time
        einsums::linear_algebra::heev<true>(&Forth->get(h), &evals);

        double last_value = - std::numeric_limits<double>::infinity();
        for (int m = 0; m < nmopi_[h]; m++) {
            const double& current_value = evals(m);
            if (last_value > current_value + 1e-16) throw PSIEXCEPTION("CGHF Orbitals are not ordered!");
            epsilon_->set(h, m, current_value);
            last_value = current_value;
        }
    }

    // heev retuns the wrong side, so we need to take the conjugate transpose for the proper eigenvectors
    auto temp = Forth->conjT();

    // Form C_ := X_ @ C' (temp)
    C_ = linalg::doublet<false, false>(X_, temp);
    C_->set_name("MO coefficients");

    find_occupation();
}

void CGHF::find_occupation() {
    if (debug_) epsilon_->print("outfile");

	if (nirrep_ > 1) throw PSIEXCEPTION("CGHF::find_occupation() C1 symmetry only!");

    // No changing C_ because we assume the orbitals are sorted within their irreps

    std::vector<std::pair<double, int>> pvec;
    pvec.reserve(nmo_);
    for (int h = 0; h < nirrep_; h++) {
        for (int m = 0; m < nmopi_[h]; m++) {
            pvec.emplace_back(epsilon_->get(h, m), h);
        }
    }

    std::sort(pvec.begin(), pvec.end());

    for (auto& e : nelecpi_) { e = 0; }

    // Aufbau principle
    for (int i = 0; i < nelectron_; i++) {
        const auto& [energy, irrep] = pvec[i];
        nelecpi_[irrep] += 1;
    }
}

void CGHF::form_D() {
    D_->zero();
    for (int h = 0; h < nirrep_; ++h) {
        int nso = C_->rowdim(h);
        int nocc = static_cast<int>(nelecpi_[h]);
        if (!nso || !nocc) continue;

        auto const& C_h = C_->get(h);
        auto& D_h = D_->get(h);

        // D_h = C_occ * C_occ^H using the leading occupied columns (lda = full MO stride)
        einsums::blas::gemm('n', 'c', nso, nso, nocc, std::complex<double>{1.0}, C_h.data(),
                            static_cast<int>(C_h.stride(0)), C_h.data(), static_cast<int>(C_h.stride(0)),
                            std::complex<double>{0.0}, D_h.data(), static_cast<int>(D_h.stride(0)));
    }

    if (debug_ > 0) {
        outfile->Printf("in CGHF::form_D:\n");
        D_->print("outfile");
    }
}

void CGHF::form_F() {
    (*F_) = (*H_);
    (*F_) += (*G_);

    F_->set_name("Fock");

    if (debug_) {
        F_->print("outfile");
        J_->print("outfile");
        K_->print("outfile");
        G_->print("outfile");
    }
}

void CGHF::form_G() {
    G_->zero();

    // Push the C matrix
    std::vector<SharedComplexMatrix>& C = jk_->C_left();
    C.clear();

    // Grab the occupied columns
    auto C2 = std::make_shared<ComplexMatrix>("C SO OCC", C_->rowspi(), nelecpi_);
    for (int h = 0; h < nirrep_; h++) {
        if (nelecpi_[h] == 0) continue;
        // Einsums analog of Matrix's per-column C_DCOPY
        C2->get(h) = C_->get(h)(einsums::All, einsums::Range{0, nelecpi_[h]});
    }
    C.push_back(C2);

    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedComplexMatrix>& J = jk_->J();
    const std::vector<SharedComplexMatrix>& K = jk_->K();
    J_ = J[0];
    K_ = K[0];

    double alpha = functional_->x_alpha();
    if (alpha != 1) throw PSIEXCEPTION("Who let the DFT in?");
    if (!functional_->is_x_hybrid()) throw PSIEXCEPTION("Who let the DFT in?");

    G_->add(*J_);
    G_->subtract(*K_);
}

std::shared_ptr<ComplexJK> CGHF::build_jk(size_t memory) const {
    auto jk = ComplexJK::build_JK(get_basisset("ORBITAL"), get_basisset("DF_BASIS_SCF"),
                                  Process::environment.options, false, memory);
    return jk;
}

void CGHF::openorbital_scf() {
#ifndef USING_OpenOrbitalOptimizer
    throw PSIEXCEPTION(
        "OpenOrbitalOptimizer support has not been enabled in this Psi4 build! Reconfigure with `-D "
        "ENABLE_OpenOrbitalOptimizer=ON`.\n");
#else
    // Store AO-basis DIIS error to use for convergence instead of orthogonal-basis error
    double ao_basis_diis_error = 1.0;

    std::function<OpenOrbitalOptimizer::FockBuilderReturn<std::complex<double>, double>(
        const OpenOrbitalOptimizer::DensityMatrix<std::complex<double>, double>&)>
        fock_builder =
            [&](const OpenOrbitalOptimizer::DensityMatrix<std::complex<double>, double>& dm) {
                // Grab the orbitals and occupations
                std::vector<arma::cx_mat> orbitals = dm.first;
                std::vector<arma::vec> occupations = dm.second;
                assert(orbitals.size() == static_cast<size_t>(nirrep_));
                assert(occupations.size() == static_cast<size_t>(nirrep_));

                // Throw away zero occupations and calculate size of the matrix
                Dimension nmopi(nirrep_);
                Dimension nsopi_spin(nirrep_);
                for (int h = 0; h < nirrep_; h++) {
                    nsopi_spin[h] = X_->rowdim(h);
                    if (nsopi_spin[h] == 0) {
                        nmopi[h] = 0;
                        continue;
                    }
                    arma::uvec idx(arma::find(occupations[h] != 0.0));
                    orbitals[h] = orbitals[h].cols(idx);
                    occupations[h] = occupations[h](idx);
                    nmopi[h] = idx.n_elem;
                    // This interface can't handle negative occupations
                    if (idx.n_elem > 0 and arma::min(occupations[h]) < 0.0) {
                        throw PSIEXCEPTION("Negative orbital occupations not supported in Psi4 interface!\n");
                    }
                }

                // Form the dummy orbitals for the jk object (weighted by sqrt(occupation))
                auto Cdummy = std::make_shared<ComplexMatrix>("Dummy orbitals", nsopi_spin, nmopi);
                for (int h = 0; h < nirrep_; h++) {
                    if (nmopi[h] == 0) continue;
                    const arma::cx_mat Xblock(X_->to_armadillo_matrix(h));
                    arma::cx_mat Cblock = Xblock * orbitals[h] * arma::diagmat(arma::sqrt(occupations[h]));
                    Cdummy->from_armadillo_matrix(Cblock, h);
                }

                std::vector<SharedComplexMatrix>& jkC = jk_->C_left();
                jkC.clear();
                jkC.push_back(Cdummy);
                jk_->compute();
                const std::vector<SharedComplexMatrix>& Jvec = jk_->J();
                const std::vector<SharedComplexMatrix>& Kvec = jk_->K();

                // CGHF is Hartree–Fock only; no DFT exchange-correlation.
                std::vector<arma::cx_mat> fock(nirrep_);
                double Ecore = 0.0, Ecoul = 0.0, Eexch = 0.0;
                double diis_sum = 0.0;
                long diis_terms = 0;
                for (int h = 0; h < nirrep_; h++) {
                    if (nsopi_spin[h] == 0) continue;
                    const arma::cx_mat Xblock(X_->to_armadillo_matrix(h));
                    const arma::cx_mat Sblock(S_->to_armadillo_matrix(h));
                    const arma::cx_mat J_AO(Jvec[0]->to_armadillo_matrix(h));
                    const arma::cx_mat K_AO(Kvec[0]->to_armadillo_matrix(h));
                    const arma::cx_mat coreH(H_->to_armadillo_matrix(h));
                    const arma::cx_mat Cblock(Cdummy->to_armadillo_matrix(h));

                    // Orthogonal-basis Fock: X^H (H + J - K) X
                    arma::cx_mat F_AO = coreH + J_AO - K_AO;
                    fock[h] = Xblock.t() * F_AO * Xblock;

                    arma::cx_mat P_AO(Cblock * Cblock.t());
                    Ecore += std::real(arma::trace(P_AO * coreH));
                    Ecoul += 0.5 * std::real(arma::trace(J_AO * P_AO));
                    Eexch += -0.5 * std::real(arma::trace(K_AO * P_AO));

                    // AO-basis DIIS error: FDS - SDF, projected like HF::form_FDSmSDF
                    arma::cx_mat FDSmSDF = F_AO * P_AO * Sblock - Sblock * P_AO * F_AO;
                    FDSmSDF = Xblock.t() * FDSmSDF * Xblock;
                    diis_sum += std::real(arma::accu(FDSmSDF % arma::conj(FDSmSDF)));
                    diis_terms += FDSmSDF.n_elem;
                }
                double Etot = Ecore + Ecoul + Eexch + nuclearrep_;
                ao_basis_diis_error = (diis_terms > 0) ? std::sqrt(diis_sum / static_cast<double>(diis_terms)) : 0.0;

                return std::make_pair(Etot, fock);
            };

    arma::uword nirrep(nirrep_);
    // One particle type (generalized spinors), with nirrep blocks
    arma::uvec number_of_blocks_per_particle_type({nirrep});
    // Each spinor has maximal occupation 1
    arma::vec maximum_occupation(nirrep, arma::fill::value(1.0));
    arma::vec number_of_particles({(double)nelectron_});

    std::vector<std::string> block_descriptions(nirrep);
    CharacterTable ct = molecule_->point_group()->char_table();
    for (int h = 0; h < nirrep_; h++) block_descriptions[h] = ct.gamma(h).symbol();

    double start_diis = options_.get_double("SCF_INITIAL_START_DIIS_TRANSITION");
    double finish_diis = options_.get_double("SCF_INITIAL_FINISH_DIIS_TRANSITION");
    int maxvecs = options_.get_double("DIIS_MAX_VECS");
    int maxiter = options_.get_int("MAXITER");
    bool fail_on_maxiter = options_.get_bool("FAIL_ON_MAXITER");

    // Get the orbital guess
    std::vector<arma::cx_mat> orbitals(nirrep);
    std::vector<arma::vec> occupations(nirrep);
    for (int h = 0; h < nirrep_; h++) {
        if (X_->rowdim(h) == 0) continue;

        auto Xblock(X_->to_armadillo_matrix(h));
        auto Sblock(S_->to_armadillo_matrix(h));
        auto Cblock(C_->to_armadillo_matrix(h));

        if (Cblock.n_cols) {
            orbitals[h] = Xblock.t() * Sblock * Cblock;
            occupations[h].zeros(Cblock.n_cols);
            if (nelecpi_[h] > 0) occupations[h].subvec(0, nelecpi_[h] - 1).ones();
        }
    }

    std::function<bool(const std::map<std::string, std::any>&)> callback_convergence_function =
        [&](const std::map<std::string, std::any>& data) -> bool {
        double e_conv = options_.get_double("E_CONVERGENCE");
        double d_conv = options_.get_double("D_CONVERGENCE");
        double e_delta = std::any_cast<double>(data.at("dE"));
        // Like UHF: no closed-shell 0.5 factor on the AO DIIS error
        double d_rms = ao_basis_diis_error;

        bool converged = (fabs(e_delta) < e_conv) && (d_rms < d_conv);
        if (iteration_ == options_.get_int("MAXITER"))
            printf("%d = |%e| < %e && %e < %e\n", converged, e_delta, e_conv, d_rms, d_conv);
        return converged;
    };

    std::function<void(const std::map<std::string, std::any>&)> callback_function =
        [&](const std::map<std::string, std::any>& data) {
            std::string reference = options_.get_str("REFERENCE");
            if (options_.get_str("SCF_TYPE").ends_with("DF")) reference = "DF-" + reference;

            iteration_++;
            double E = std::any_cast<double>(data.at("E"));
            double dE = std::any_cast<double>(data.at("dE"));
            double Dnorm = ao_basis_diis_error;
            std::string step = std::any_cast<std::string>(data.at("step"));

            outfile->Printf("   @%s iter %3i: %20.14f   %12.5e   %-11.5e %s\n", reference.c_str(), iteration_, E, dE,
                            Dnorm, step.c_str());
        };

    OpenOrbitalOptimizer::SCFSolver<std::complex<double>, double> scfsolver(
        number_of_blocks_per_particle_type, maximum_occupation, number_of_particles, fock_builder, block_descriptions);
    scfsolver.maximum_iterations(maxiter);
    scfsolver.verbosity(options_.get_int("OOO_PRINT"));
    scfsolver.maximum_history_length(maxvecs);
    scfsolver.oda_restart_steps(options_.get_int("OOO_ODA_RESTART_STEPS"));
    scfsolver.callback_function(callback_function);
    scfsolver.callback_convergence_function(callback_convergence_function);
    scfsolver.diis_epsilon(start_diis);
    scfsolver.diis_threshold(finish_diis);
    scfsolver.diis_restart_factor(options_.get_double("OOO_DIIS_RESTART_FACTOR"));
    scfsolver.optimal_damping_threshold(options_.get_double("OOO_OPTIMAL_DAMPING_THRESHOLD"));
    scfsolver.initialize_with_orbitals(orbitals, occupations);
    scfsolver.run();
    if (fail_on_maxiter and not scfsolver.converged())
        throw PSIEXCEPTION("SCF did not converge and FAIL_ON_MAXITER is set to true.\n");

    // Update the orbitals with OOO's solution
    auto solution = scfsolver.get_solution();
    orbitals = solution.first;
    occupations = solution.second;
    for (int h = 0; h < nirrep_; h++) {
        if (X_->rowdim(h) == 0) continue;
        const arma::cx_mat Xblock(X_->to_armadillo_matrix(h));
        arma::cx_mat Cblock(Xblock * orbitals[h]);
        C_->from_armadillo_matrix(Cblock, h);

        // Occupations should be integer 0/1 for CGHF spinors
        arma::uvec nonzero(arma::find(occupations[h] > 0.0));
        double diff(arma::norm(occupations[h](nonzero) - arma::ones<arma::vec>(nonzero.n_elem), 2));
        if (diff != 0.0) {
            fprintf(stderr, "Non-integer occupations in symmetry block %i\n", h);
            occupations[h].print("occs");
        }
        nelecpi_[h] = (int)nonzero.n_elem;
    }

    // Form the density matrix
    form_D();
    // Form the two-electron part
    form_G();
    // Compute the energy
    compute_E();
#endif
}

}  // namespace scf
}  // namespace psi
