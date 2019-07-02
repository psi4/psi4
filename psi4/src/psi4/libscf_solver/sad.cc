/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

/*
 *  sad.cc
 *
 * Routines for the high-maintenance SAD guess
 * and dual-basis projections
 *
 */
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/factory.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libdiis/diisentry.h"
#include "psi4/libfock/jk.h"
#include "psi4/lib3index/dfhelper.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "hf.h"
#include "sad.h"

using namespace psi;

namespace psi {
namespace scf {
// Parse options: either use DF or exact integrals.
static bool SAD_use_fitting(const Options& opt) {
    std::string jk_type(opt.get_str("SAD_SCF_TYPE"));
    if ((jk_type == "DIRECT") || (jk_type == "PK") || (jk_type == "OUT_OF_CORE") || (jk_type == "CD") ||
        (jk_type == "GTFOCK")) {
        return false;
    }
    if ((jk_type == "DF") || (jk_type == "MEM_DF") || (jk_type == "DISK_DF")) {
        return true;
    }
    throw PSIEXCEPTION("SAD_SCF_TYPE " + opt.get_str("SAD_SCF_TYPE") + " not implemented.\n");
}

SADGuess::SADGuess(std::shared_ptr<BasisSet> basis, std::vector<std::shared_ptr<BasisSet>> atomic_bases,
                   Options& options)
    : basis_(basis), atomic_bases_(atomic_bases), options_(options) {
    common_init();
}
SADGuess::~SADGuess() {}
void SADGuess::common_init() {
    molecule_ = basis_->molecule();

    auto ints = std::make_shared<IntegralFactory>(basis_);
    auto petite = std::make_shared<PetiteList>(basis_, ints);
    AO2SO_ = petite->aotoso();

    print_ = options_.get_int("SAD_PRINT");
    debug_ = options_.get_int("DEBUG");
    if (options_["SOCC"].size() > 0 || options_["DOCC"].size() > 0)
        PSIEXCEPTION("SAD guess not implemented for user-specified SOCCs and/or DOCCs yet");
}
void SADGuess::compute_guess() {
    timer_on("SAD Guess");
    start_skip_timers();
    form_D();
    form_C();
    stop_skip_timers();
    timer_off("SAD Guess");
}
void SADGuess::form_D() {
    // Build Neutral D in AO basis (block diagonal)
    SharedMatrix DAO;
    // Huckel matrices
    SharedMatrix HuckelC;
    SharedVector HuckelE;
    run_atomic_calculations(DAO, HuckelC, HuckelE);

    // Transform Neutral D from AO to SO basis
    Da_ = std::make_shared<Matrix>("Da SAD", AO2SO_->colspi(), AO2SO_->colspi());
    Da_->apply_symmetry(DAO, AO2SO_);

    // Set Db to Da
    Db_ = Da_;

    if (debug_) {
        Da_->print();
        Db_->print();
    }
}
void SADGuess::form_C() {
    Ca_ = Da_->partial_cholesky_factorize(options_.get_double("SAD_CHOL_TOLERANCE"));
    Ca_->set_name("Ca SAD");
    Cb_ = Ca_;

    if (debug_) {
        Ca_->print();
        Cb_->print();
    }
}
void SADGuess::run_atomic_calculations(SharedMatrix& DAO, SharedMatrix& HuckelC, SharedVector& HuckelE) {
    if (print_ > 6) {
        for (int A = 0; A < molecule_->natom(); A++) {
            outfile->Printf("  SAD: Atomic Basis Set %d\n", A);
            atomic_bases_[A]->molecule()->print();
            outfile->Printf("\n");
            atomic_bases_[A]->print("outfile");
            outfile->Printf("\n");
        }
    }

    // Spin occupations per atom, to be determined by Hund's Rules
    // or user input
    std::vector<double> nalpha(molecule_->natom(), 0);
    std::vector<double> nbeta(molecule_->natom(), 0);
    std::vector<int> nelec(molecule_->natom(), 0);

    // Ground state high spin occupency array, atoms 0 to 36 (see Giffith's Quantum Mechanics, pp. 217)
    // For 37 to 86, save for f-block: Atomic, Molecular, & Optical Physics Handbook, Ed. Gordon W. F. Drake, American
    // Institute of Physics, Woodbury, New York, USA, 1996.

    // clang-format off
    const std::vector<int> reference_S = { 0,
                                           1,                                                                                           0,
                                           1, 0,                                                                         1, 2, 3, 2, 1, 0,
                                           1, 0,                                                                         1, 2, 3, 2, 1, 0,
                                           1, 0,                                           1, 2, 3, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0,
                                           1, 0,                                           1, 2, 5, 6, 5, 4, 3, 0, 1, 0, 1, 2, 3, 2, 1, 0,
                                           1, 0, 1, 0, 3, 4, 5, 6, 7, 8, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0 };
    // clang-format on

    if (print_ > 1) outfile->Printf("  Determining Atomic Occupations\n");

    for (int A = 0; A < molecule_->natom(); A++) {
        int Z = std::round(molecule_->Z(A));
        // Number of ECP electrons on center
        int ECP = basis_->n_ecp_core(molecule_->label(A));
        // Assuming neutral atoms for now
        nelec[A] = Z;

        if (options_.get_bool("SAD_SPIN_AVERAGE")) {
            // Spin-averaged occupations
            nalpha[A] = nbeta[A] = 0.5 * nelec[A];
        } else {
            // Target ground spin state
            if (Z + ECP >= reference_S.size()) {
                std::ostringstream err;
                err << " Only atoms up to Z = " << reference_S.size() - 1
                    << " are currently supported with SAD_SPIN_AVERAGE = false\n";
                throw std::domain_error(err.str());
            }
            int nhigh = reference_S[Z + ECP];
            nalpha[A] = 0.5 * (nelec[A] + nhigh);
            nbeta[A] = 0.5 * (nelec[A] - nhigh);
        }

        if (print_ > 1) outfile->Printf("  Atom %d, Z = %d, nalpha = %.1f, nbeta = %.1f\n", A, Z, nalpha[A], nbeta[A]);
    }

    // Determine redundant atoms
    std::vector<int> unique_indices(molecule_->natom(), 0);  // All atoms to representative unique atom
    std::vector<int> atomic_indices(molecule_->natom(), 0);  // unique atom to first representative atom
    std::vector<int> offset_indices(molecule_->natom(), 0);  // unique atom index to rank
    int nunique = 0;
    for (int l = 0; l < molecule_->natom(); l++) {
        unique_indices[l] = l;
        atomic_indices[l] = l;
    }

    // This is all overkill for now, but lets leave it incase we do something oddball in the future
    for (int l = 0; l < molecule_->natom() - 1; l++) {
        for (int m = l + 1; m < molecule_->natom(); m++) {
            if (unique_indices[m] != m) continue;  // Already assigned
            if (molecule_->Z(l) != molecule_->Z(m)) continue;
            if (nalpha[l] != nalpha[m]) continue;
            if (nbeta[l] != nbeta[m]) continue;
            if (atomic_bases_[l]->nbf() != atomic_bases_[m]->nbf()) continue;
            if (atomic_bases_[l]->nshell() != atomic_bases_[m]->nshell()) continue;
            if (atomic_bases_[l]->nprimitive() != atomic_bases_[m]->nprimitive()) continue;
            if (atomic_bases_[l]->max_am() != atomic_bases_[m]->max_am()) continue;
            if (atomic_bases_[l]->max_nprimitive() != atomic_bases_[m]->max_nprimitive()) continue;
            if (atomic_bases_[l]->has_puream() != atomic_bases_[m]->has_puream()) continue;
            if (atomic_bases_[l]->n_ecp_core() != atomic_bases_[m]->n_ecp_core()) continue;

            // Semi-Rigorous match obtained
            unique_indices[m] = l;
        }
    }
    for (int l = 0; l < molecule_->natom(); l++) {
        if (unique_indices[l] == l) {
            atomic_indices[nunique] = l;
            offset_indices[l] = nunique;
            nunique++;
        }
    }

    // Atomic density matrices
    std::vector<SharedMatrix> atomic_D(nunique);
    // Atomic orbitals for Huckel
    std::vector<SharedMatrix> atomic_Chu(nunique);
    // Atomic orbital energies for Huckel
    std::vector<SharedVector> atomic_Ehu(nunique);

    if (print_ > 1) outfile->Printf("\n  Performing Atomic UHF Computations:\n");
    for (int uniA = 0; uniA < nunique; uniA++) {
        int index = atomic_indices[uniA];
        int nbf = atomic_bases_[index]->nbf();
        int Z = molecule_->Z(index);
        if (nelec[index] == 0) {
            // No electrons on atom!
            continue;
        }

        if (print_ > 1) {
            outfile->Printf("\n  UHF Computation for Unique Atom %d which is Atom %d:\n", uniA, index);
            outfile->Printf("  Occupation: nalpha = %.1f, nbeta = %.1f, nbf = %d\n", nalpha[uniA], nbeta[uniA], nbf);
        }

        // Occupation numbers
        SharedVector occ_a, occ_b;
        // Number of orbitals occupied, partially or fully
        int nocc_a, nocc_b;
        if (options_.get_bool("SAD_FRAC_OCC")) {
            // The density is spread over the whole of the possible valence shell.
            // Only the noble gas core is doubly occupied
            static const std::vector<int> magic_values = {0, 2, 10, 18, 36, 54, 86, 118};
            // Find the noble gas core.
            auto imagic = std::lower_bound(magic_values.begin(), magic_values.end(), Z);
            if (imagic == magic_values.end()) {
                throw PSIEXCEPTION("SAD: Fractional occupations are not supported beyond Oganesson");
            }
            // lower_bound gives a value that is equal or greater than the
            // wanted value, which is handled next.
            if (*imagic > Z) imagic--;

            // Number of frozen and active orbitals
            int nfzc, nact;
            if ((*imagic) == Z) {
                // Special case: we can hit the boundary at the end of the array
                nfzc = 0;
                nact = (*imagic) / 2;
            } else {
                nfzc = (*imagic) / 2;
                nact = (*(++imagic)) / 2 - nfzc;
            }

            // Sanity check: can't have more active orbitals than basis functions
            if (nact > nbf - nfzc) {
                nact = nbf - nfzc;
            }

            // Number of occupied orbitals is
            nocc_a = nocc_b = nfzc + nact;

            // Fractional alpha and beta occupation. Occupations are
            // squared in the density calculation, so take the root
            double frac_a = std::sqrt((nalpha[index] - nfzc) / nact);
            double frac_b = std::sqrt((nbeta[index] - nfzc) / nact);

            occ_a = std::make_shared<Vector>("Alpha fractional occupation", nocc_a);
            for (size_t x = 0; x < nfzc; x++) occ_a->set(x, 1.0);
            for (size_t x = nfzc; x < nocc_a; x++) occ_a->set(x, frac_a);

            occ_b = std::make_shared<Vector>("Beta fractional occupation", nocc_b);
            for (size_t x = 0; x < nfzc; x++) occ_b->set(x, 1.0);
            for (size_t x = nfzc; x < nocc_b; x++) occ_b->set(x, frac_b);

            if (print_ > 1) {
                outfile->Printf(
                    "  %d fully and %d partially occupied orbitals with % .3f alpha and % .3f beta electrons "
                    "each\n",
                    nfzc, nact, (nalpha[index] - nfzc) / nact, (nbeta[index] - nfzc) / nact);
            }
        } else {
            // Conventional occupations
            nocc_a = std::round(nalpha[index]);
            occ_a = std::make_shared<Vector>("Alpha occupation", nocc_a);
            for (size_t x = 0; x < nocc_a; x++) occ_a->set(x, 1.0);

            nocc_b = std::round(nbeta[index]);
            occ_b = std::make_shared<Vector>("Beta occupation", nocc_b);
            for (size_t x = 0; x < nocc_b; x++) occ_b->set(x, 1.0);
        }

        int nhu = occ_a->dim();
        atomic_D[uniA] = std::make_shared<Matrix>("Atomic D_AO", nbf, nbf);
        atomic_Chu[uniA] = std::make_shared<Matrix>("Atomic Huckel C", nbf, nhu);
        atomic_Ehu[uniA] = std::make_shared<Vector>("Atomic Huckel E", nhu);

        if (SAD_use_fitting(options_)) {
            get_uhf_atomic_density(atomic_bases_[index], atomic_fit_bases_[index], occ_a, occ_b, atomic_D[uniA],
                                   atomic_Chu[uniA], atomic_Ehu[uniA]);
        } else {
            std::shared_ptr<BasisSet> zbas = BasisSet::zero_ao_basis_set();
            get_uhf_atomic_density(atomic_bases_[index], zbas, occ_a, occ_b, atomic_D[uniA], atomic_Chu[uniA],
                                   atomic_Ehu[uniA]);
        }
        if (print_ > 1) outfile->Printf("Finished UHF Computation!\n");
    }
    if (print_) outfile->Printf("\n");

    // Add atomic_D into D (scale by 1/2, we like effective pairs)
    DAO = std::make_shared<Matrix>("D_SAD (AO)", basis_->nbf(), basis_->nbf());
    DAO->zero();
    for (int A = 0, offset = 0; A < molecule_->natom(); A++) {
        int nbf = atomic_bases_[A]->nbf();
        // Handle ghost atoms
        if (nelec[A] > 0) {
            int back_index = unique_indices[A];
            for (int m = 0; m < nbf; m++)
                for (int n = 0; n < nbf; n++)
                    DAO->set(0, m + offset, n + offset, 0.5 * atomic_D[offset_indices[back_index]]->get(m, n));
        }
        offset += nbf;
    }

    // Total number of Huckel orbitals
    int nhuckel = 0;
    for (int A = 0, offset = 0; A < molecule_->natom(); A++) {
        // Handle ghost atoms
        if (nelec[A] > 0) {
            int back_index = unique_indices[A];
            int uniA = offset_indices[back_index];
            nhuckel += atomic_Chu[uniA]->coldim();
        }
    }

    // Collect Huckel orbital coefficients
    HuckelC = std::make_shared<Matrix>("C_Huckel (MINAO)", basis_->nbf(), nhuckel);
    HuckelE = std::make_shared<Vector>("E_Huckel (MINAO)", nhuckel);
    HuckelC->zero();
    HuckelE->zero();
    for (int A = 0, offset = 0, ioffset = 0; A < molecule_->natom(); A++) {
        int nbf = atomic_bases_[A]->nbf();
        // Handle ghost atoms
        if (nelec[A] > 0) {
            int back_index = unique_indices[A];
            int uniA = offset_indices[back_index];
            int nhu = atomic_Chu[uniA]->coldim();
            assert(atomic_Chu[uniA]->rowdim() == nbf);
            for (int ibf = 0; ibf < nbf; ibf++)
                for (int io = 0; io < nhu; io++)
                    HuckelC->set(0, ibf + offset, io + ioffset, atomic_Chu[uniA]->get(ibf, io));
            for (int io = 0; io < nhu; io++) HuckelE->set(io + ioffset, atomic_Ehu[uniA]->get(io));
            ioffset += nhu;
        }
        offset += nbf;
    }

    if (debug_) {
        DAO->print();
        HuckelC->print();
        HuckelE->print();
    }
}
void SADGuess::get_uhf_atomic_density(std::shared_ptr<BasisSet> bas, std::shared_ptr<BasisSet> fit, SharedVector occ_a,
                                      SharedVector occ_b, SharedMatrix D, SharedMatrix Chuckel, SharedVector Ehuckel) {
    std::shared_ptr<Molecule> mol = bas->molecule();
    mol->update_geometry();
    if (print_ > 1) {
        mol->print();
    }

    int natom = mol->natom();
    int nbf = bas->nbf();
    int Z = bas->molecule()->Z(0);

    if (occ_a->dim() > nbf || occ_b->dim() > nbf) throw PSIEXCEPTION("Atom has more electrons than basis functions.");

    if (print_ > 1) {
        outfile->Printf("\n");
        bas->print("outfile");
        outfile->Printf("\n  Atom:\n");
        mol->print();
    }

    if (natom != 1) {
        throw std::domain_error("SAD Atomic UHF has been given a molecule, not an atom");
    }

    IntegralFactory integral(bas, bas, bas, bas);
    MatrixFactory mat;
    mat.init_with(1, &nbf, &nbf);
    std::unique_ptr<OneBodyAOInt> S_ints = std::unique_ptr<OneBodyAOInt>(integral.ao_overlap());
    std::unique_ptr<OneBodyAOInt> T_ints = std::unique_ptr<OneBodyAOInt>(integral.ao_kinetic());
    std::unique_ptr<OneBodyAOInt> V_ints = std::unique_ptr<OneBodyAOInt>(integral.ao_potential());
    std::unique_ptr<OneBodyAOInt> ECP_ints = std::unique_ptr<OneBodyAOInt>(integral.ao_ecp());

    // Compute overlap S and orthogonalizer X;
    SharedMatrix S(mat.create_matrix("Overlap Matrix"));
    S_ints->compute(S);

    SharedMatrix X = S->clone();
    X->power(-0.5, 1.e-10);
    X->set_name("Orthogonalizer X^-1/2 Matrix");

    if (print_ > 6) {
        S->print();
        X->print();
    }

    // Compute H
    SharedMatrix T(mat.create_matrix("T"));
    T_ints->compute(T);
    SharedMatrix V(mat.create_matrix("V"));
    V_ints->compute(V);
    SharedMatrix ECP(mat.create_matrix("ECP"));
    ECP_ints->compute(ECP);
    SharedMatrix H(mat.create_matrix("Core Hamiltonian Matrix H"));
    H->zero();
    H->add(T);
    H->add(V);
    H->add(ECP);

    T.reset();
    V.reset();
    ECP.reset();

    if (print_ > 6) {
        H->print();
    }

    // Init temps
    SharedMatrix Ca(mat.create_matrix("Ca"));
    SharedMatrix Cb(mat.create_matrix("Cb"));

    SharedMatrix Da(mat.create_matrix("Da"));
    SharedMatrix Db(mat.create_matrix("Db"));

    SharedVector Ea = std::make_shared<Vector>("Ea", nbf);
    SharedVector Eb = std::make_shared<Vector>("Eb", nbf);

    SharedMatrix gradient_a(mat.create_matrix("gradient_a"));
    SharedMatrix gradient_b(mat.create_matrix("gradient_b"));

    SharedMatrix Fa(mat.create_matrix("Fa"));
    SharedMatrix Fb(mat.create_matrix("Fb"));

    auto Ca_occ = std::make_shared<Matrix>("Ca occupied", nbf, occ_a->dim());
    auto Cb_occ = std::make_shared<Matrix>("Cb occupied", nbf, occ_b->dim());

    // Compute initial Cx, Dx, and D from core guess
    form_C_and_D(X, H, Ca, Ea, Ca_occ, occ_a, Da);
    form_C_and_D(X, H, Cb, Eb, Cb_occ, occ_b, Db);

    D->zero();
    D->add(Da);
    D->add(Db);

    if (print_ > 6) {
        Ca->print();
        Cb->print();
        Ca_occ->print();
        Cb_occ->print();
        Da->print();
        Db->print();
        D->print();
    }

    // Compute initial E for reference
    double E = D->vector_dot(H);
    E *= 0.5;

    double E_tol = options_.get_double("SAD_E_CONVERGENCE");
    double D_tol = options_.get_double("SAD_D_CONVERGENCE");
    int sad_maxiter = options_.get_int("SAD_MAXITER");
    bool diis_rms = options_.get_bool("DIIS_RMS_ERROR");

    double E_old = E;
    int iteration = 0;

    // Setup DIIS
    DIISManager diis_manager(6, "SAD DIIS", DIISManager::LargestError, DIISManager::InCore);
    diis_manager.set_error_vector_size(2, DIISEntry::Matrix, gradient_a.get(), DIISEntry::Matrix, gradient_b.get());
    diis_manager.set_vector_size(2, DIISEntry::Matrix, Fa.get(), DIISEntry::Matrix, Fb.get());

    // Setup JK
    std::unique_ptr<JK> jk;

    // Need a very special auxiliary basis here
    if (SAD_use_fitting(options_)) {
        MemDFJK* dfjk = new MemDFJK(bas, fit);
        if (options_["DF_INTS_NUM_THREADS"].has_changed())
            dfjk->set_df_ints_num_threads(options_.get_int("DF_INTS_NUM_THREADS"));
        dfjk->dfh()->set_print_lvl(0);
        jk = std::unique_ptr<JK>(dfjk);
    } else {
        DirectJK* directjk(new DirectJK(bas));
        if (options_["DF_INTS_NUM_THREADS"].has_changed())
            directjk->set_df_ints_num_threads(options_.get_int("DF_INTS_NUM_THREADS"));
        jk = std::unique_ptr<JK>(directjk);
    }

    jk->set_memory((size_t)(0.5 * (Process::environment.get_memory() / 8L)));
    jk->initialize();
    if (print_ > 1) jk->print_header();

    // These are static so lets just grab them now
    std::vector<SharedMatrix>& jkC = jk->C_left();
    jkC.push_back(Ca_occ);
    jkC.push_back(Cb_occ);
    const std::vector<SharedMatrix>& Jvec = jk->J();
    const std::vector<SharedMatrix>& Kvec = jk->K();

    // Print a header
    bool converged = false;
    if (print_ > 1) {
        std::string measure = diis_rms ? "RMS |[F,P]|  " : "MAX |[F,P]|  ";
        outfile->Printf("\n  Initial Atomic UHF Energy:    %14.10f\n\n", E);
        outfile->Printf("  %33s %20s    %20s %20s\n", "", "Total Energy   ", "Delta E   ", measure.c_str());
    }

    // Run the iterations
    do {
        iteration++;

        // Copy the old values over for error analysis
        E_old = E;

        // Compute JK matrices
        jk->compute();

        // Form Fa and Fb
        Fa->copy(H);
        Fa->add(Jvec[0]);
        Fa->add(Jvec[1]);

        Fb->copy(Fa);

        Fa->subtract(Kvec[0]);
        Fb->subtract(Kvec[1]);

        // Compute E
        E = H->vector_dot(D);
        E += Da->vector_dot(Fa);
        E += Db->vector_dot(Fb);
        E *= 0.5;

        double deltaE = std::fabs(E - E_old);

        // Build Gradient
        form_gradient(gradient_a, Fa, Da, S, X);
        form_gradient(gradient_b, Fb, Db, S, X);
        double Dnorm = diis_rms ? std::sqrt(0.5 * (std::pow(gradient_a->rms(), 2) + std::pow(gradient_b->rms(), 2)))
                                : std::max(gradient_a->absmax(), gradient_b->absmax());

        // Add and extrapolate DIIS
        diis_manager.add_entry(4, gradient_a.get(), gradient_b.get(), Fa.get(), Fb.get());
        diis_manager.extrapolate(2, Fa.get(), Fb.get());

        // Diagonalize Fa and Fb to from Ca and Cb and Da and Db
        form_C_and_D(X, Fa, Ca, Ea, Ca_occ, occ_a, Da);
        form_C_and_D(X, Fb, Cb, Eb, Cb_occ, occ_b, Db);

        // Form D
        D->copy(Da);
        D->add(Db);

        if (print_ > 6) {
            H->print();
            Fa->print();
            Fb->print();
            Ca->print();
            Cb->print();
            Da->print();
            Db->print();
            D->print();
        }
        if (print_ > 1)
            outfile->Printf("  @Atomic UHF iteration %3d energy: %20.14f    %20.14f %20.14f\n", iteration, E, E - E_old,
                            Dnorm);

        // Check convergence
        if (iteration > 1) {
            converged = (deltaE < E_tol && Dnorm < D_tol);
        }

        if (iteration > sad_maxiter) {
            outfile->Printf(
                "\n WARNING: Atomic UHF is not converging! Try casting from a smaller basis or call Rob at CCMST.\n");
            break;
        }

    } while (!converged);

    if (converged && print_ > 1)
        outfile->Printf("  @Atomic UHF Final Energy for atom %s: %20.14f\n", mol->symbol(0).c_str(), E);

    // Copy Huckel coefficients and energies
    double** Coccp = Chuckel->pointer();
    double** Cp = Ca->pointer();
    for (int i = 0; i < nbf; i++) {
        C_DCOPY(occ_a->dim(), Cp[i], 1, Coccp[i], 1);
    }
    double* Eoccp = Ehuckel->pointer();
    double* Ep = Ea->pointer();
    for (int i = 0; i < occ_a->dim(); i++) {
        Eoccp[i] = Ep[i];
    }
}
void SADGuess::form_gradient(SharedMatrix grad, SharedMatrix F, SharedMatrix D, SharedMatrix S, SharedMatrix X) {
    int nbf = X->rowdim();
    auto Scratch1 = std::make_shared<Matrix>("Scratch1", nbf, nbf);
    auto Scratch2 = std::make_shared<Matrix>("Scratch2", nbf, nbf);

    // FDS
    Scratch1->gemm(false, false, 1.0, F, D, 0.0);
    Scratch2->gemm(false, false, 1.0, Scratch1, S, 0.0);

    // SDF
    Scratch1->copy(Scratch2);
    Scratch1->transpose_this();

    // FDS - SDF
    grad->copy(Scratch2);
    grad->subtract(Scratch1);

    // Level it out X(FDS - SDF)X
    Scratch1->gemm(false, false, 1.0, X, grad, 0.0);
    grad->gemm(false, false, 1.0, Scratch1, X, 0.0);

    Scratch1.reset();
    Scratch2.reset();
}

void SADGuess::form_C_and_D(SharedMatrix X, SharedMatrix F, SharedMatrix C, SharedVector E, SharedMatrix Cocc,
                            SharedVector occ, SharedMatrix D) {
    int nbf = X->rowdim();
    int nocc = occ->dim();
    if (nocc == 0) return;

    // Forms C in the AO basis for SAD Guesses
    auto Scratch1 = std::make_shared<Matrix>("Scratch1", nbf, nbf);
    auto Scratch2 = std::make_shared<Matrix>("Scratch2", nbf, nbf);

    // Form Fp = XFX
    Scratch1->gemm(true, false, 1.0, X, F, 0.0);
    Scratch2->gemm(false, false, 1.0, Scratch1, X, 0.0);

    // Diagonalize
    Scratch2->diagonalize(Scratch1, E);

    // Form C = XC'
    C->gemm(false, false, 1.0, X, Scratch1, 0.0);

    // Copy over Cocc
    double** Coccp = Cocc->pointer();
    double** Cp = C->pointer();
    for (int i = 0; i < nbf; i++) {
        C_DCOPY(nocc, Cp[i], 1, Coccp[i], 1);
    }
    // Scale by occ
    for (int i = 0; i < nocc; i++) {
        C_DSCAL(nbf, occ->get(i), &Coccp[0][i], nocc);
    }
    // Form D = Cocc*Cocc'
    D->gemm(false, true, 1.0, Cocc, Cocc, 0.0);

    Scratch1.reset();
    Scratch2.reset();
}
SharedMatrix SADGuess::huckel_guess() {
    // Build Neutral D in AO basis (block diagonal)
    SharedMatrix DAO;
    // Huckel matrices
    SharedMatrix Chu;
    SharedVector Ehu;
    run_atomic_calculations(DAO, Chu, Ehu);

    IntegralFactory integral(basis_, basis_, basis_, basis_);
    MatrixFactory mat;

    int nbf = basis_->nbf();
    mat.init_with(1, &nbf, &nbf);
    std::unique_ptr<OneBodyAOInt> S_ints = std::unique_ptr<OneBodyAOInt>(integral.ao_overlap());

    // Compute overlap S
    SharedMatrix S(mat.create_matrix("Overlap Matrix"));
    S_ints->compute(S);

    // Compute Huckel basis overlap S*Chu
    int nhu = Chu->coldim();

    // Compute S*Chu
    auto SChu = std::make_shared<Matrix>("SChu", nbf, nhu);
    SChu->gemm(false, false, 1.0, S, Chu, 0.0);

    // Compute Chu^T*S*Chu
    auto ChuSChu = std::make_shared<Matrix>("ChuSChu", nhu, nhu);
    ChuSChu->gemm(true, false, 1.0, Chu, SChu, 0.0);

    // Huckel matrix in Huckel basis
    auto huckelmo = std::make_shared<Matrix>("Huckel MO matrix", nhu, nhu);
    double** huckelmop = huckelmo->pointer();
    for (int i = 0; i < nhu; i++) {
        huckelmo->set(i, i, Ehu->get(i));
        for (int j = 0; j < nhu; j++) {
            huckelmo->set(i, j, 0.875 * ChuSChu->get(i, j) * (Ehu->get(i) + Ehu->get(j)));
        }
    }

    // Half-transform
    auto scratch = std::make_shared<Matrix>("Scratch memory", nbf, nhu);
    scratch->gemm(false, false, 1.0, SChu, huckelmo, 0.0);

    // Full transform to AO basis
    auto huckelao = std::make_shared<Matrix>("Huckel AO matrix", nbf, nbf);
    huckelao->gemm(false, true, 1.0, SChu, scratch, 0.0);

    // Now, transform from AO to SO basis
    auto huckel = std::make_shared<Matrix>("Huckel SO matrix", AO2SO_->colspi(), AO2SO_->colspi());
    huckel->apply_symmetry(huckelao, AO2SO_);

    return huckel;
}
void HF::compute_SAD_guess(bool natorb) {
    if (sad_basissets_.empty()) {
        throw PSIEXCEPTION("  SCF guess was set to SAD, but sad_basissets_ was empty!\n\n");
    }

    auto guess = std::make_shared<SADGuess>(basisset_, sad_basissets_, options_);
    if (SAD_use_fitting(options_)) {
        if (sad_fitting_basissets_.empty()) {
            throw PSIEXCEPTION("  SCF guess was set to SAD with DiskDFJK, but sad_fitting_basissets_ was empty!\n\n");
        }
        guess->set_atomic_fit_bases(sad_fitting_basissets_);
    }

    guess->compute_guess();

    if (natorb) {
        // SAD natural orbitals (doi:10.1021/acs.jctc.8b01089)

        // Number of basis functions
        auto nbf(basisset_->nbf());
        // Grab the density matrix
        auto Dhelp = std::make_shared<Matrix>("Helper density matrix", AO2SO_->colspi(), AO2SO_->colspi());
        Dhelp->copy(guess->Da());
        // Take its negative so that most strongly occupied orbitals come first
        Dhelp->scale(-1.0);

        // Transfrom to an orthonormal basis
        Dhelp->transform(S_);
        Dhelp->transform(X_);
        if (debug_) {
            outfile->Printf("SAD Density Matrix (orthonormal basis):\n");
            Dhelp->print();
        }

        // Diagonalize the SAD density and form the natural orbitals
        auto Cno_temp_ = SharedMatrix(factory_->create_matrix("SAD NO temp"));
        Dhelp->diagonalize(Cno_temp_, epsilon_a_);
        if (debug_) {
            outfile->Printf("SAD Natural Orbital Occupations:\n");
            epsilon_a_->print();
        }
        Ca_->gemm(false, false, 1.0, X_, Cno_temp_, 0.0);

        // Same orbitals for beta
        Cb_->copy(Ca_);
        epsilon_b_->copy(*epsilon_a_);

    } else {
        SharedMatrix Ca_sad = guess->Ca();
        SharedMatrix Cb_sad = guess->Cb();
        Da_->copy(guess->Da());
        Db_->copy(guess->Db());
        Dimension sad_dim(Da_->nirrep(), "SAD Dimensions");

        for (int h = 0; h < Da_->nirrep(); h++) {
            int nso = Ca_sad->rowspi()[h];
            int nmo = Ca_sad->colspi()[h];
            if (nmo > X_->colspi()[h]) nmo = X_->colspi()[h];

            sad_dim[h] = nmo;

            if (!nso || !nmo) continue;

            double** Cap = Ca_->pointer(h);
            double** Cbp = Cb_->pointer(h);
            double** Ca2p = Ca_sad->pointer(h);
            double** Cb2p = Cb_sad->pointer(h);
            for (int i = 0; i < nso; i++) {
                C_DCOPY(nmo, Ca2p[i], 1, Cap[i], 1);
                C_DCOPY(nmo, Cb2p[i], 1, Cbp[i], 1);
            }
        }

        nalphapi_ = sad_dim;
        nbetapi_ = sad_dim;
        nalpha_ = sad_dim.sum();
        nbeta_ = sad_dim.sum();
        doccpi_ = sad_dim;
        soccpi_ = Dimension(Da_->nirrep(), "SAD SOCC dim (0's)");
        energies_["Total Energy"] = 0.0;  // This is the -1th iteration
    }
}
void HF::compute_huckel_guess() {
    if (sad_basissets_.empty()) {
        throw PSIEXCEPTION("  SCF guess was set to SAD, but sad_basissets_ was empty!\n\n");
    }

    auto guess = std::make_shared<SADGuess>(basisset_, sad_basissets_, options_);
    if (SAD_use_fitting(options_)) {
        if (sad_fitting_basissets_.empty()) {
            throw PSIEXCEPTION("  SCF guess was set to SAD with DiskDFJK, but sad_fitting_basissets_ was empty!\n\n");
        }
        guess->set_atomic_fit_bases(sad_fitting_basissets_);
    }

    SharedMatrix Fhuckel = guess->huckel_guess();
    Fa_->copy(Fhuckel);
    Fb_->copy(Fhuckel);

    energies_["Total Energy"] = 0.0;  // This is the -1th iteration
}
}  // namespace scf
}  // namespace psi
