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
#include "sad.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/complexmatrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libfock/ComplexJK.h"
#include "psi4/libfock/v.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#ifdef USING_OpenOrbitalOptimizer
#include <openorbitaloptimizer/scfsolver.hpp>
#include <any>
#endif

#include <cmath>
#include <complex>

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

    if (options_.get_int("MOM_START") != 0) throw PSIEXCEPTION("MOM not available for CGHF.");

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
    S_ = std::make_shared<ComplexMatrix>("Overlap");
    X_ = std::make_shared<ComplexMatrix>("Orthogonalization");
    // C_ is resized in form_C once X_ is known
    C_ = std::make_shared<ComplexMatrix>("MO coefficients");

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
        outfile->Printf("  MOM not available for Complex SCF.\n");
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

    T_ = ComplexMatrix::spin_blocked_from(T_real);
    V_ = ComplexMatrix::spin_blocked_from(V_real);

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

    // Create the correctly sized complex quantities now.
    S_ = ComplexMatrix::spin_blocked_from(S_temp);
    X_ = ComplexMatrix::spin_blocked_from(X_temp);
}

void CGHF::guess() {
    double guess_E;
    std::string guess_type = options_.get_str("GUESS");

    if ((guess_type == "READ") && !guess_C_) {
        outfile->Printf("\nWarning! Guess was READ without C set, switching to CORE!\n");
        outfile->Printf("           This option should have been configured at the driver level.\n\n");
        guess_type = "CORE";
    }

    if (guess_C_) {
        if (print_) outfile->Printf("  SCF Guess: Orbitals guess was supplied from a previous computation.\n\n");

        if (guess_C_->nirrep() != nirrep_) {
            throw PSIEXCEPTION(
                "Number of irreps of the input orbitals do not match number of irreps of the wavefunction.");
        }

        // Full C_ sized to current orthog basis; copy leading occupied columns from guess
        C_ = std::make_shared<ComplexMatrix>("MO coefficients", X_->rowspi(), nmopi_);
        C_->zero();

        for (int h = 0; h < nirrep_; h++) {
            const int nso_spin = X_->rowdim(h);
            const int nocc = guess_C_->coldim(h);
            if (guess_C_->rowdim(h) != nso_spin) {
                throw PSIEXCEPTION("Nso of the guess orbitals do not match Nso of the wavefunction.");
            }
            if (nocc > nmopi_[h]) {
                throw PSIEXCEPTION("Guess has more occupied orbitals than available MOs.");
            }
            nelecpi_[h] = nocc;
            for (int i = 0; i < nso_spin; i++) {
                for (int j = 0; j < nocc; j++) {
                    C_->set(h, i, j, guess_C_->get(h, i, j));
                }
            }
        }

        form_D();
        iteration_ = -1;
        guess_E = compute_initial_E();

    } else if (guess_type == "CORE") {
        if (print_) outfile->Printf("  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");

        F_->add(*H_); // Try the core Hamiltonian as the Fock Matrix
        form_initial_C(); // calls only CGHF::form_C()
        form_D();
        guess_E = compute_initial_E();
    } else if (guess_type == "SAD") {
        if (print_)
            outfile->Printf(
                "  SCF Guess: Superposition of Atomic Densities via on-the-fly atomic UHF (no occupation "
                "information).\n\n");

        compute_SAD_guess();
        // Guess iteration: occupations must be reset in SCF; SAD has no usable MOs yet.
        iteration_ = -1;
        sad_ = true;
        guess_E = compute_initial_E();
    } else {
        throw PSIEXCEPTION("CGHF '" + guess_type + "' GUESS not implemented. Use 'CORE', 'SAD', or 'READ'.");
    }

    energies_["Total Energy"] = 0.0;  // don't use this guess in our convergence checks
}

void CGHF::compute_SAD_guess() {
    if (sad_basissets_.empty()) {
        throw PSIEXCEPTION("  SCF guess was set to SAD, but sad_basissets_ was empty!\n\n");
    }

    auto guess = std::make_shared<SADGuess>(basisset_, sad_basissets_, options_);
    // Driver only fills fitting bases when SAD_SCF_TYPE is DF-like.
    if (!sad_fitting_basissets_.empty()) {
        guess->set_atomic_fit_bases(sad_fitting_basissets_);
    }

    guess->compute_guess();

    // Spin-restricted SAD: Da == Db → block-diagonal spinor density [[Da,0],[0,Da]]
    D_ = ComplexMatrix::spin_blocked_from(guess->Da());

    // Embed Cholesky factors as temporary occupied spinors for form_G on SAD iter 0:
    // Ca in alpha spatial block (cols [0, nchol)), Cb in beta block (cols [nchol, 2*nchol)).
    SharedMatrix Ca_sad = guess->Ca();
    SharedMatrix Cb_sad = guess->Cb();
    C_ = std::make_shared<ComplexMatrix>("MO coefficients", X_->rowspi(), nmopi_);
    C_->zero();

    for (int h = 0; h < nirrep_; h++) {
        int nso = Ca_sad->rowspi()[h];
        int nchol = Ca_sad->colspi()[h];
        // nmopi_ is already spin-doubled; clamp Cholesky columns to spatial MO count
        int nmo_spatial = nmopi_[h] / 2;
        if (nchol > nmo_spatial) nchol = nmo_spatial;

        nelecpi_[h] = 2 * nchol;

        if (!nso || !nchol) continue;

        for (int i = 0; i < nso; i++) {
            for (int j = 0; j < nchol; j++) {
                C_->set(h, i, j, Ca_sad->get(h, i, j));
                C_->set(h, i + nso, j + nchol, Cb_sad->get(h, i, j));
            }
        }
    }

    energies_["Total Energy"] = 0.0;  // This is the -1th iteration
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

    auto [evals, evecs] = linalg::diagonalize(*F_, *X_);
    epsilon_ = evals;

    // Form C_ := X_ @ C' (evecs)
    C_ = linalg::doublet(X_, evecs);
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
    // Build C_occ: the occupied-column submatrix of the spinor coefficient matrix.
    // C_ is (2*nsopi_) × (2*nsopi_), nelecpi_ occupied spinor columns per irrep.
    auto C_occ = std::make_shared<ComplexMatrix>("C_occ", C_->rowspi(), nelecpi_);

    for (int h = 0; h < nirrep_; ++h) {
        int nso = C_->rowdim(h);
        int nocc = nelecpi_[h];
        if (!nso || !nocc) continue;

        // Copy occupied columns: C_occ(h) = C(h)[:, :nocc]
        C_occ->get(h) = C_->get(h)(einsums::All, einsums::Range{0, nocc});
    }

    // D = C_occ * C_occ^H   (conjugate-transpose via doublet(..., false, true))
    D_ = linalg::doublet(C_occ, C_occ, false, true);
    D_->set_name("Density");

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

std::tuple<double, double> CGHF::spin_square() const {
    // Spatial SO overlap (not spin-blocked). Alpha/beta blocks of C_ share this metric.
    ComplexMatrix S{mintshelper()->so_overlap()};
    S.set_name("S (spatial)");

    // Occupied alpha / beta MO coefficients from the spinor C_ (rows: [α; β]).
    ComplexMatrix Ca("Ca occ", nsopi_, nelecpi_);
    ComplexMatrix Cb("Cb occ", nsopi_, nelecpi_);
    for (int h = 0; h < nirrep_; h++) {
        const int nso = nsopi_[h];
        const int nocc = nelecpi_[h];
        if (!nso || !nocc) continue;
        for (int i = 0; i < nso; i++) {
            for (int j = 0; j < nocc; j++) {
                Ca.set(h, i, j, C_->get(h, i, j));
                Cb.set(h, i, j, C_->get(h, i + nso, j));
            }
        }
    }

    // Sσστ = Cσ^H S Cτ
    auto Saa = linalg::triplet(Ca, S, Ca, true, false, false);
    auto Sbb = linalg::triplet(Cb, S, Cb, true, false, false);
    auto Sab = linalg::triplet(Ca, S, Cb, true, false, false);
    auto Sba = linalg::triplet(Cb, S, Ca, true, false, false);

    const auto nocc_a = Saa.trace();
    const auto nocc_b = Sbb.trace();

    // ⟨S+S- + S-S+⟩/2 contribution
    std::complex<double> ssxy = (nocc_a + nocc_b) * 0.5;
    ssxy += Sba.trace() * Sab.trace() - Sba.product_trace(Sab);

    // ⟨Sz²⟩ contribution
    std::complex<double> ssz = (nocc_a + nocc_b) * 0.25;
    ssz += (nocc_a - nocc_b) * (nocc_a - nocc_b) * 0.25;
    ComplexMatrix tmp = Saa;
    tmp.subtract(Sbb);
    ssz -= tmp.product_trace(tmp) * 0.25;

    const double ss = (ssxy + ssz).real();
    const double s = std::sqrt(ss + 0.25) - 0.5;
    const double multiplicity = 2.0 * s + 1.0;

    if (ss < 5e-10) {
        outfile->Printf("   @S^2:              0\n");
        outfile->Printf("   @2S+1:             1\n\n");
    } else {
        outfile->Printf("   @S^2:              %17.9f\n", ss);
        outfile->Printf("   @2S+1:             %17.9f\n\n", multiplicity);
    }

    return std::make_tuple(ss, multiplicity);
}

void CGHF::check_phases() {
    // Complex counterpart of HF::check_phases: for each MO column, find the first
    // significant AO coefficient and multiply the column by conj(c)/|c| so that
    // element is real and positive. For real orbitals this reduces to a ±1 flip.
    for (int h = 0; h < nirrep_; ++h) {
        const int nso = C_->rowdim(h);
        const int nmo = C_->coldim(h);
        for (int p = 0; p < nmo; ++p) {
            for (int mu = 0; mu < nso; ++mu) {
                const auto c = C_->get(h, mu, p);
                if (std::abs(c) > 1.0E-3) {
                    const auto phase = std::conj(c) / std::abs(c);
                    for (int nu = 0; nu < nso; ++nu) {
                        C_->set(h, nu, p, C_->get(h, nu, p) * phase);
                    }
                    break;
                }
            }
        }
    }
}

SharedComplexMatrix CGHF::form_FDSmSDF(SharedComplexMatrix Fso, SharedComplexMatrix Dso) {
    // FDS - SDF in the (spin-blocked) SO basis, then project with X like HF::form_FDSmSDF /
    // the OOO DIIS error (X^H e X). For Hermitian F, D, S this matches (FDS) - (FDS)^H.
    auto FDSmSDF = linalg::triplet(Fso, Dso, S_);
    auto SDF = linalg::triplet(S_, Dso, Fso);
    FDSmSDF->subtract(*SDF);
    FDSmSDF->set_name("FDS-SDF");

    auto orthog = linalg::triplet(X_, FDSmSDF, X_, true, false, false);
    orthog->set_name("Orthogonal FDS-SDF");
    return orthog;
}

}  // namespace scf
}  // namespace psi
