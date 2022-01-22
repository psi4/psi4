/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "local.h"

#include <algorithm>

#include "psi4/libqt/qt.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

using namespace psi;

namespace psi {

Localizer::Localizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C) : primary_(primary), C_(C) {
    if (C->nirrep() != 1) {
        throw PSIEXCEPTION("Localizer: C matrix is not C1");
    }
    if (C->rowspi()[0] != primary->nbf()) {
        throw PSIEXCEPTION("Localizer: C matrix does not match basis");
    }
    common_init();
}
Localizer::~Localizer() {}
void Localizer::common_init() {
    print_ = 0;
    debug_ = 0;
    bench_ = 0;
    convergence_ = 1.0E-8;
    maxiter_ = 50;
    converged_ = false;
}
std::shared_ptr<Localizer> Localizer::build(const std::string& type, std::shared_ptr<BasisSet> primary,
                                            std::shared_ptr<Matrix> C, Options& options) {
    std::shared_ptr<Localizer> local;

    if (type == "BOYS") {
        local = std::make_shared<BoysLocalizer>(primary, C);
    } else if (type == "PIPEK_MEZEY") {
        local = std::make_shared<PMLocalizer>(primary, C);
    } else {
        throw PSIEXCEPTION("Localizer: Unrecognized localization algorithm");
    }

    local->set_print(options.get_int("PRINT"));
    local->set_debug(options.get_int("DEBUG"));
    local->set_bench(options.get_int("BENCH"));
    local->set_convergence(options.get_double("LOCAL_CONVERGENCE"));
    local->set_maxiter(options.get_int("LOCAL_MAXITER"));

    return local;
}
std::shared_ptr<Localizer> Localizer::build(const std::string& type, std::shared_ptr<BasisSet> primary,
                                            std::shared_ptr<Matrix> C) {
    return Localizer::build(type, primary, C, Process::environment.options);
}
std::shared_ptr<Localizer> Localizer::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C,
                                            Options& options) {
    return Localizer::build(options.get_str("LOCAL_TYPE"), primary, C, options);
}
std::shared_ptr<Matrix> Localizer::fock_update(std::shared_ptr<Matrix> Fc) {
    if (!L_ || !U_) {
        throw PSIEXCEPTION("Localizer: run compute() first");
    }

    int nso = L_->rowspi()[0];
    int nmo = L_->colspi()[0];

    if (nmo < 1) return Fc;

    std::shared_ptr<Matrix> Fl = linalg::triplet(U_, Fc, U_, true, false, false);
    double** Fp = Fl->pointer();
    double** Lp = L_->pointer();
    double** Up = U_->pointer();

    std::vector<std::pair<double, int> > order;
    for (int i = 0; i < nmo; i++) {
        order.push_back(std::pair<double, int>(Fp[i][i], i));
    }
    std::sort(order.begin(), order.end());

    std::shared_ptr<Matrix> Fl2(Fl->clone());
    Fl2->copy(Fl);
    double** F2p = Fl2->pointer();
    for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < nmo; j++) {
            Fp[i][j] = F2p[order[i].second][order[j].second];
        }
    }

    std::shared_ptr<Matrix> L2(L_->clone());
    L2->copy(L_);
    double** L2p = L2->pointer();
    std::shared_ptr<Matrix> U2(U_->clone());
    U2->copy(U_);
    double** U2p = U2->pointer();
    for (int i = 0; i < nmo; i++) {
        C_DCOPY(nso, &L2p[0][order[i].second], nmo, &Lp[0][i], nmo);
        C_DCOPY(nmo, &U2p[0][order[i].second], nmo, &Up[0][i], nmo);
    }

    return Fl;
}

BoysLocalizer::BoysLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C) : Localizer(primary, C) {
    common_init();
}
BoysLocalizer::~BoysLocalizer() {}
void BoysLocalizer::common_init() {}
void BoysLocalizer::print_header() const {
    outfile->Printf("  ==> Boys Localizer <==\n\n");
    outfile->Printf("    Convergence = %11.3E\n", convergence_);
    outfile->Printf("    Maxiter     = %11d\n", maxiter_);
    outfile->Printf("\n");
}
void BoysLocalizer::localize() {
    print_header();

    // => Sizing <= //

    int nso = C_->rowspi()[0];
    int nmo = C_->colspi()[0];

    // => Dipole Integrals <= //

    auto fact = std::make_shared<IntegralFactory>(primary_);
    std::shared_ptr<OneBodyAOInt> Dint(fact->ao_dipole());
    std::vector<std::shared_ptr<Matrix> > D;
    for (int xyz = 0; xyz < 3; xyz++) {
        D.push_back(std::make_shared<Matrix>("D", nso, nso));
    }
    Dint->compute(D);
    Dint.reset();
    fact.reset();

    std::vector<std::shared_ptr<Matrix> > Dmo;
    auto T = std::make_shared<Matrix>("T", nso, nmo);
    for (int xyz = 0; xyz < 3; xyz++) {
        Dmo.push_back(std::make_shared<Matrix>("D", nmo, nmo));
        C_DGEMM('N', 'N', nso, nmo, nso, 1.0, D[xyz]->pointer()[0], nso, C_->pointer()[0], nmo, 0.0, T->pointer()[0],
                nmo);
        C_DGEMM('T', 'N', nmo, nmo, nso, 1.0, C_->pointer()[0], nmo, T->pointer()[0], nmo, 0.0, Dmo[xyz]->pointer()[0],
                nmo);
    }
    D.clear();
    T.reset();

    // => Targets <= //

    L_ = std::make_shared<Matrix>("L", nso, nmo);
    U_ = std::make_shared<Matrix>("U", nmo, nmo);
    L_->copy(C_);
    U_->identity();
    converged_ = false;

    if (nmo < 1) return;

    // => Pointers <= //

    std::vector<double**> Dp;
    for (int xyz = 0; xyz < 3; xyz++) {
        Dp.push_back(Dmo[xyz]->pointer());
    }
    double** Cp = C_->pointer();
    double** Lp = L_->pointer();
    double** Up = U_->pointer();

    // => Seed the random idempotently <= //

    srand(0L);

    // => Metric <= //

    double metric = 0.0;
    for (int xyz = 0; xyz < 3; xyz++) {
        metric += C_DDOT(nmo, Dp[xyz][0], nmo + 1, Dp[xyz][0], nmo + 1);
    }
    double old_metric = metric;

    // => Iteration Print <= //
    outfile->Printf("    Iteration %24s %14s\n", "Metric", "Residual");
    outfile->Printf("    @Boys %4d %24.16E %14s\n", 0, metric, "-");

    // ==> Master Loop <== //

    double Ad, Ao, a, b, c, Hd, Ho, theta, cc, ss;

    for (int iter = 1; iter <= maxiter_; iter++) {
        // => Random Permutation <= //

        std::vector<int> order;
        for (int i = 0; i < nmo; i++) {
            order.push_back(i);
        }
        std::vector<int> order2;
        for (int i = 0; i < nmo; i++) {
            int pivot = (1L * (nmo - i) * rand()) / RAND_MAX;
            int i2 = order[pivot];
            order[pivot] = order[nmo - i - 1];
            order2.push_back(i2);
        }

        // => Jacobi sweep <= //

        for (int i2 = 0; i2 < nmo - 1; i2++) {
            for (int j2 = i2 + 1; j2 < nmo; j2++) {
                int i = order2[i2];
                int j = order2[j2];

                // > Compute the rotation < //

                // H elements
                a = 0.0;
                b = 0.0;
                c = 0.0;
                for (int xyz = 0; xyz < 3; xyz++) {
                    double** Ak = Dp[xyz];
                    Ad = (Ak[i][i] - Ak[j][j]);
                    Ao = 2.0 * Ak[i][j];
                    a += Ad * Ad;
                    b += Ao * Ao;
                    c += Ad * Ao;
                }

                // Theta
                Hd = a - b;
                Ho = 2.0 * c;
                theta = 0.5 * atan2(Ho, Hd + sqrt(Hd * Hd + Ho * Ho));

                // Check for trivial (maximal) rotation, which might be better with theta = pi/4
                if (std::fabs(theta) < 1.0E-8) {
                    double O0 = 0.0;
                    double O1 = 0.0;
                    ;
                    for (int xyz = 0; xyz < 3; xyz++) {
                        double** Ak = Dp[xyz];
                        O0 += Ak[i][j] * Ak[i][j];
                        O1 += 0.25 * (Ak[j][j] - Ak[i][i]) * (Ak[j][j] - Ak[i][i]);
                    }
                    if (O1 < O0) {
                        theta = M_PI / 4.0;
                        if (debug_ > 3) {
                            outfile->Printf("@Break\n");
                        }
                    }
                }

                // Givens rotation
                cc = cos(theta);
                ss = sin(theta);

                if (debug_ > 3) {
                    outfile->Printf("@Rotation, i = %4d, j = %4d, Theta = %24.16E\n", i, j, theta);
                    outfile->Printf("@Info, a = %24.16E, b = %24.16E, c = %24.16E\n", a, b, c);
                }

                // > Apply the rotation < //

                // rows and columns of A^k
                for (int xyz = 0; xyz < 3; xyz++) {
                    double** Ak = Dp[xyz];
                    C_DROT(nmo, &Ak[i][0], 1, &Ak[j][0], 1, cc, ss);
                    C_DROT(nmo, &Ak[0][i], nmo, &Ak[0][j], nmo, cc, ss);
                }

                // Q
                C_DROT(nmo, Up[i], 1, Up[j], 1, cc, ss);
            }
        }

        // => Metric <= //

        metric = 0.0;
        for (int xyz = 0; xyz < 3; xyz++) {
            metric += C_DDOT(nmo, Dp[xyz][0], nmo + 1, Dp[xyz][0], nmo + 1);
        }

        double conv = std::fabs(metric - old_metric) / std::fabs(old_metric);
        old_metric = metric;

        // => Iteration Print <= //

        outfile->Printf("    @Boys %4d %24.16E %14.6E\n", iter, metric, conv);

        // => Convergence Check <= //

        if (conv < convergence_) {
            converged_ = true;
            break;
        }
    }

    outfile->Printf("\n");
    if (converged_) {
        outfile->Printf("    Boys Localizer converged.\n\n");
    } else {
        outfile->Printf("    Boys Localizer failed.\n\n");
    }

    U_->transpose_this();
    C_DGEMM('N', 'N', nso, nmo, nmo, 1.0, Cp[0], nmo, Up[0], nmo, 0.0, Lp[0], nmo);
}

PMLocalizer::PMLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C) : Localizer(primary, C) {
    common_init();
}
PMLocalizer::~PMLocalizer() {}
void PMLocalizer::common_init() {}
void PMLocalizer::print_header() const {
    outfile->Printf("  ==> Pipek-Mezey Localizer <==\n\n");
    outfile->Printf("    Convergence = %11.3E\n", convergence_);
    outfile->Printf("    Maxiter     = %11d\n", maxiter_);
    outfile->Printf("\n");
}
void PMLocalizer::localize() {
    print_header();

    // => Sizing <= //

    int nso = C_->rowspi()[0];
    int nmo = C_->colspi()[0];
    int nA = primary_->molecule()->natom();

    // => Overlap Integrals <= //

    auto fact = std::make_shared<IntegralFactory>(primary_);
    std::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());
    auto S = std::make_shared<Matrix>("S", nso, nso);
    Sint->compute(S);
    Sint.reset();
    fact.reset();
    double** Sp = S->pointer();

    // => Targets <= //

    L_ = std::make_shared<Matrix>("L", nso, nmo);
    U_ = std::make_shared<Matrix>("U", nmo, nmo);
    L_->copy(C_);
    U_->identity();
    converged_ = false;

    if (nmo < 1) return;

    // => Pointers <= //

    double** Lp = L_->pointer();
    double** Up = U_->pointer();

    // => LS product (avoids GEMV) <= //

    auto LS = std::make_shared<Matrix>("LS", nso, nmo);
    double** LSp = LS->pointer();
    C_DGEMM('N', 'N', nso, nmo, nso, 1.0, Sp[0], nso, Lp[0], nmo, 0.0, LSp[0], nmo);

    // => Starting functions on each atomic center <= //

    std::vector<int> Astarts;
    int Aoff = 0;
    for (int m = 0; m < primary_->nbf(); m++) {
        if (primary_->function_to_center(m) == Aoff) {
            Astarts.push_back(m);
            Aoff++;
        }
    }
    Astarts.push_back(primary_->nbf());

    // => Seed the random idempotently <= //

    srand(0L);

    // => Metric <= //

    double metric = 0.0;
    for (int i = 0; i < nmo; i++) {
        for (int A = 0; A < nA; A++) {
            int nm = Astarts[A + 1] - Astarts[A];
            int off = Astarts[A];
            double PA = C_DDOT(nm, &LSp[off][i], nmo, &Lp[off][i], nmo);
            metric += PA * PA;
        }
    }
    double old_metric = metric;

    // => Iteration Print <= //
    outfile->Printf("    Iteration %24s %14s\n", "Metric", "Residual");
    outfile->Printf("    @PM %4d %24.16E %14s\n", 0, metric, "-");

    // ==> Master Loop <== //

    double Aii, Ajj, Aij, Ad, Ao, a, b, c, Hd, Ho, theta, cc, ss;

    for (int iter = 1; iter <= maxiter_; iter++) {
        // => Random Permutation <= //

        std::vector<int> order;
        for (int i = 0; i < nmo; i++) {
            order.push_back(i);
        }
        std::vector<int> order2;
        for (int i = 0; i < nmo; i++) {
            int pivot = (1L * (nmo - i) * rand()) / RAND_MAX;
            int i2 = order[pivot];
            order[pivot] = order[nmo - i - 1];
            order2.push_back(i2);
        }

        // => Jacobi sweep <= //

        for (int i2 = 0; i2 < nmo - 1; i2++) {
            for (int j2 = i2 + 1; j2 < nmo; j2++) {
                int i = order2[i2];
                int j = order2[j2];

                // > Compute the rotation < //

                // H elements
                a = 0.0;
                b = 0.0;
                c = 0.0;
                for (int A = 0; A < nA; A++) {
                    int nm = Astarts[A + 1] - Astarts[A];
                    int off = Astarts[A];
                    Aii = C_DDOT(nm, &LSp[off][i], nmo, &Lp[off][i], nmo);
                    Ajj = C_DDOT(nm, &LSp[off][j], nmo, &Lp[off][j], nmo);
                    Aij = 0.5 * C_DDOT(nm, &LSp[off][i], nmo, &Lp[off][j], nmo) +
                          0.5 * C_DDOT(nm, &LSp[off][j], nmo, &Lp[off][i], nmo);

                    Ad = (Aii - Ajj);
                    Ao = 2.0 * Aij;
                    a += Ad * Ad;
                    b += Ao * Ao;
                    c += Ad * Ao;
                }

                // Theta
                Hd = a - b;
                Ho = 2.0 * c;
                theta = 0.5 * atan2(Ho, Hd + sqrt(Hd * Hd + Ho * Ho));

                // Check for trivial (maximal) rotation, which might be better with theta = pi/4
                if (std::fabs(theta) < 1.0E-8) {
                    double O0 = 0.0;
                    double O1 = 0.0;
                    ;
                    for (int A = 0; A < nA; A++) {
                        int nm = Astarts[A + 1] - Astarts[A];
                        int off = Astarts[A];
                        Aii = C_DDOT(nm, &LSp[off][i], nmo, &Lp[off][i], nmo);
                        Ajj = C_DDOT(nm, &LSp[off][j], nmo, &Lp[off][j], nmo);
                        Aij = 0.5 * C_DDOT(nm, &LSp[off][i], nmo, &Lp[off][j], nmo) +
                              0.5 * C_DDOT(nm, &LSp[off][j], nmo, &Lp[off][i], nmo);
                        O0 += Aij * Aij;
                        O1 += 0.25 * (Ajj - Aii) * (Ajj - Aii);
                    }
                    if (O1 < O0) {
                        theta = M_PI / 4.0;
                        if (debug_ > 3) {
                            outfile->Printf("@Break\n");
                        }
                    }
                }

                // Givens rotation
                cc = cos(theta);
                ss = sin(theta);

                if (debug_ > 3) {
                    outfile->Printf("@Rotation, i = %4d, j = %4d, Theta = %24.16E\n", i, j, theta);
                    outfile->Printf("@Info, a = %24.16E, b = %24.16E, c = %24.16E\n", a, b, c);
                }

                // > Apply the rotation < //

                // columns of LS and L
                C_DROT(nso, &LSp[0][i], nmo, &LSp[0][j], nmo, cc, ss);
                C_DROT(nso, &Lp[0][i], nmo, &Lp[0][j], nmo, cc, ss);

                // Q
                C_DROT(nmo, Up[i], 1, Up[j], 1, cc, ss);
            }
        }

        // => Metric <= //

        metric = 0.0;
        for (int i = 0; i < nmo; i++) {
            for (int A = 0; A < nA; A++) {
                int nm = Astarts[A + 1] - Astarts[A];
                int off = Astarts[A];
                double PA = C_DDOT(nm, &LSp[off][i], nmo, &Lp[off][i], nmo);
                metric += PA * PA;
            }
        }

        double conv = std::fabs(metric - old_metric) / std::fabs(old_metric);
        old_metric = metric;

        // => Iteration Print <= //

        outfile->Printf("    @PM %4d %24.16E %14.6E\n", iter, metric, conv);

        // => Convergence Check <= //

        if (conv < convergence_) {
            converged_ = true;
            break;
        }
    }

    outfile->Printf("\n");
    if (converged_) {
        outfile->Printf("    PM Localizer converged.\n\n");
    } else {
        outfile->Printf("    PM Localizer failed.\n\n");
    }

    U_->transpose_this();
}

IBOLocalizer::IBOLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> minao, std::shared_ptr<Matrix> C,
                            std::shared_ptr<Matrix> F, const std::vector<int>& ranges)
    : Localizer(primary, C), minao_(minao), Focc_(F), ranges_(ranges) {
    common_init();
}
IBOLocalizer::~IBOLocalizer() {}
void IBOLocalizer::common_init() {
    print_ = 0;
    debug_ = 0;
    bench_ = 0;
    convergence_ = 1.0E-12;
    maxiter_ = 50;
    use_ghosts_ = false;
    power_ = 4;
    condition_ = 1.0E-7;
    use_stars_ = false;
    stars_completeness_ = 0.9;
    stars_.clear();
}
std::shared_ptr<IBOLocalizer> IBOLocalizer::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> minao,
                                                std::shared_ptr<Matrix> C, std::shared_ptr<Matrix> F, 
                                                const std::vector<int>& ranges, Options& options) {

    auto local = std::make_shared<IBOLocalizer>(primary, minao, C, F, ranges);

    local->set_print(options.get_int("PRINT"));
    local->set_debug(options.get_int("DEBUG"));
    local->set_bench(options.get_int("BENCH"));
    local->set_convergence(options.get_double("LOCAL_CONVERGENCE"));
    local->set_maxiter(options.get_int("LOCAL_MAXITER"));
    local->set_use_ghosts(options.get_bool("LOCAL_USE_GHOSTS"));
    local->set_condition(options.get_double("LOCAL_IBO_CONDITION"));
    local->set_power(options.get_double("LOCAL_IBO_POWER"));
    local->set_use_stars(options.get_bool("LOCAL_IBO_USE_STARS"));
    local->set_stars_completeness(options.get_double("LOCAL_IBO_STARS_COMPLETENESS"));

    std::vector<int> stars;
    for (int ind = 0; ind < options["LOCAL_IBO_STARS"].size(); ind++) {
        stars.push_back(options["LOCAL_IBO_STARS"][ind].to_integer() - 1);
    }
    local->set_stars(stars);

    return local;
}

std::shared_ptr<IBOLocalizer> IBOLocalizer::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> minao, std::shared_ptr<Matrix> C, 
                                                std::shared_ptr<Matrix> F, const std::vector<int>& ranges) {
    return IBOLocalizer::build(primary, minao, C, F, ranges, Process::environment.options);
}

void IBOLocalizer::print_header() const {
    outfile->Printf("  ==> IBO Localizer <==\n\n");
    outfile->Printf("    MinAO Basis = %11s\n", minao_->name().c_str());
    outfile->Printf("    Use Ghosts  = %11s\n", (use_ghosts_ ? "TRUE" : "FALSE"));
    outfile->Printf("    Use Stars   = %11s\n", (use_stars_ ? "TRUE" : "FALSE"));
    outfile->Printf("    Condition   = %11.3E\n", condition_);
    outfile->Printf("    Power       = %11d\n", power_);
    outfile->Printf("    Convergence = %11.3E\n", convergence_);
    outfile->Printf("    Maxiter     = %11d\n", maxiter_);
    outfile->Printf("\n");
}

void IBOLocalizer::build_iaos() {

    // => Ghosting <= //
    std::shared_ptr<Molecule> mol = minao_->molecule();
    true_atoms_.clear();
    true_iaos_.clear();
    iaos_to_atoms_.clear();
    for (int A = 0; A < mol->natom(); A++) {
        if (!use_ghosts_ && mol->Z(A) == 0.0) continue;
        int Atrue = true_atoms_.size();
        int nPshells = minao_->nshell_on_center(A);
        int sPshells = minao_->shell_on_center(A, 0);
        for (int P = sPshells; P < sPshells + nPshells; P++) {
            int nP = minao_->shell(P).nfunction();
            int oP = minao_->shell(P).function_index();
            for (int p = 0; p < nP; p++) {
                true_iaos_.push_back(p + oP);
                iaos_to_atoms_.push_back(Atrue);
            }
        }
        true_atoms_.push_back(A);
    }

    // => Overlap Integrals <= //

    auto fact11 = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    auto fact12 = std::make_shared<IntegralFactory>(primary_, minao_, primary_, minao_);
    auto fact22 = std::make_shared<IntegralFactory>(minao_, minao_, minao_, minao_);

    std::unique_ptr<OneBodyAOInt> ints11(fact11->ao_overlap());
    std::unique_ptr<OneBodyAOInt> ints12(fact12->ao_overlap());
    std::unique_ptr<OneBodyAOInt> ints22(fact22->ao_overlap());

    auto S11 = std::make_shared<Matrix>("S11", primary_->nbf(), primary_->nbf());
    auto S12f = std::make_shared<Matrix>("S12f", primary_->nbf(), minao_->nbf());
    auto S22f = std::make_shared<Matrix>("S22f", minao_->nbf(), minao_->nbf());

    ints11->compute(S11);
    ints12->compute(S12f);
    ints22->compute(S22f);

    ints11.reset();
    ints12.reset();
    ints22.reset();

    fact11.reset();
    fact12.reset();
    fact22.reset();

    // => Ghosted Overlap Integrals <= //

    auto S12 = std::make_shared<Matrix>("S12", primary_->nbf(), true_iaos_.size());
    auto S22 = std::make_shared<Matrix>("S22", true_iaos_.size(), true_iaos_.size());

    double** S12p = S12->pointer();
    double** S12fp = S12f->pointer();
    for (int m = 0; m < primary_->nbf(); m++) {
        for (int p = 0; p < true_iaos_.size(); p++) {
            S12p[m][p] = S12fp[m][true_iaos_[p]];
        }
    }

    double** S22p = S22->pointer();
    double** S22fp = S22f->pointer();
    for (int p = 0; p < true_iaos_.size(); p++) {
        for (int q = 0; q < true_iaos_.size(); q++) {
            S22p[p][q] = S22fp[true_iaos_[p]][true_iaos_[q]];
        }
    }

    // => Metric Inverses <= //

    std::shared_ptr<Matrix> S11_m12(S11->clone());
    std::shared_ptr<Matrix> S22_m12(S22->clone());
    S11_m12->copy(S11);
    S22_m12->copy(S22);
    S11_m12->power(-1.0 / 2.0, condition_);
    S22_m12->power(-1.0 / 2.0, condition_);

    // => Tilde C <= //

    auto C = C_;
    auto T1 = linalg::doublet(S22_m12, S12, false, true);
    auto T2 = linalg::doublet(S11_m12, linalg::triplet(T1, T1, C, true, false, false), false, false);
    auto T3 = linalg::doublet(T2, T2, true, false);
    T3->power(-1.0 / 2.0, condition_);
    auto Ctilde = linalg::triplet(S11_m12, T2, T3, false, false, false);

    // => D and Tilde D <= //

    auto D = linalg::doublet(C, C, false, true);
    auto Dtilde = linalg::doublet(Ctilde, Ctilde, false, true);

    // => A (Before Orthogonalization) <= //

    auto DSDtilde = linalg::triplet(D, S11, Dtilde, false, false, false);
    DSDtilde->scale(2.0);

    auto L = linalg::doublet(S11_m12, S11_m12, false, false);  // TODO: Possibly Unstable
    L->add(DSDtilde);
    L->subtract(D);
    L->subtract(Dtilde);

    auto AN = linalg::doublet(L, S12, false, false);

    // => A (After Orthogonalization) <= //

    auto V = linalg::triplet(AN, S11, AN, true, false, false);
    V->power(-1.0 / 2.0, condition_);

    auto A = linalg::doublet(AN, V, false, false);

    // => Assignment <= //

    S_ = S11;
    A_ = A;
}

std::map<std::string, std::shared_ptr<Matrix> > IBOLocalizer::localize_task(
    std::shared_ptr<Matrix> L, const std::vector<std::vector<int> >& minao_inds,
    const std::vector<std::pair<int, int> >& rot_inds, double convergence, int maxiter, int power) {
    int nmin = L->colspi()[0];
    int nocc = L->rowspi()[0];

    std::shared_ptr<Matrix> L2(L->clone());
    L2->copy(L);
    double** Lp = L2->pointer();

    auto U = std::make_shared<Matrix>("U", nocc, nocc);
    U->identity();
    double** Up = U->pointer();

    bool converged = false;

    if (power != 2 && power != 4) throw PSIEXCEPTION("IAO: Invalid metric power.");

    outfile->Printf("    @IBO %4s: %24s %14s\n", "Iter", "Metric", "Gradient");

    for (int iter = 1; iter <= maxiter; iter++) {
        double metric = 0.0;
        for (int i = 0; i < nocc; i++) {
            for (const auto& minao_ind_list : minao_inds) {
                double Lval = 0.0;
                for (const int mind : minao_ind_list) {
                    Lval += Lp[i][mind] * Lp[i][mind];
                }
                metric += pow(Lval, power);
            }
        }
        metric = pow(metric, 1.0 / power);

        double gradient = 0.0;
        for (int ind = 0; ind < rot_inds.size(); ind++) {
            int i = rot_inds[ind].first;
            int j = rot_inds[ind].second;

            double Aij = 0.0;
            double Bij = 0.0;
            for (const auto& minao_ind_list : minao_inds) {
                double Qii = 0.0;
                double Qij = 0.0;
                double Qjj = 0.0;
                for (const int mind : minao_ind_list) {
                    Qii += Lp[i][mind] * Lp[i][mind];
                    Qij += Lp[i][mind] * Lp[j][mind];
                    Qjj += Lp[j][mind] * Lp[j][mind];
                }
                if (power == 2) {
                    Aij += 4.0 * Qij * Qij - (Qii - Qjj) * (Qii - Qjj);
                    Bij += 4.0 * Qij * (Qii - Qjj);
                } else {
                    Aij += (-1.0) * Qii * Qii * Qii * Qii - Qjj * Qjj * Qjj * Qjj +
                           6.0 * (Qii * Qii + Qjj * Qjj) * Qij * Qij + Qii * Qii * Qii * Qjj + Qii * Qjj * Qjj * Qjj;
                    Bij += 4.0 * Qij * (Qii * Qii * Qii - Qjj * Qjj * Qjj);
                }
            }

            double phi = 0.25 * atan2(Bij, -Aij);
            double c = cos(phi);
            double s = sin(phi);

            C_DROT(nmin, Lp[i], 1, Lp[j], 1, c, s);
            C_DROT(nocc, Up[i], 1, Up[j], 1, c, s);

            gradient += Bij * Bij;
        }
        gradient = sqrt(gradient);

        outfile->Printf("    @IBO %4d: %24.16E %14.6E\n", iter, metric, gradient);

        if (gradient < convergence) {
            converged = true;
            break;
        }
    }

    outfile->Printf("\n");
    if (converged) {
        outfile->Printf("    IBO Localizer converged.\n\n");
    } else {
        outfile->Printf("    IBO Localizer failed.\n\n");
    }

    U->transpose_this();

    std::map<std::string, std::shared_ptr<Matrix> > ret;
    ret["U"] = U;
    ret["L"] = L2;

    ret["U"]->set_name("U");
    ret["L"]->set_name("L");

    return ret;
}

std::shared_ptr<Matrix> IBOLocalizer::reorder_orbitals(std::shared_ptr<Matrix> F, const std::vector<int>& ranges) {
    int nmo = F->rowspi()[0];
    double** Fp = F->pointer();

    auto U = std::make_shared<Matrix>("U", nmo, nmo);
    double** Up = U->pointer();

    for (int ind = 0; ind < ranges.size() - 1; ind++) {
        int start = ranges[ind];
        int stop = ranges[ind + 1];
        std::vector<std::pair<double, int> > fvals;
        for (int i = start; i < stop; i++) {
            fvals.push_back(std::pair<double, int>(Fp[i][i], i));
        }
        std::sort(fvals.begin(), fvals.end());
        for (int i = start; i < stop; i++) {
            Up[i][fvals[i - start].second] = 1.0;
        }
    }

    return U;
}

void IBOLocalizer::localize() {
    print_header();

    if (!A_) build_iaos();

    if (!ranges_.size()) {
        ranges_.push_back(0);
        ranges_.push_back(C_->colspi()[0]);
    }

    std::vector<std::vector<int> > minao_inds;
    for (int A = 0; A < true_atoms_.size(); A++) {
        std::vector<int> vec;
        for (int m = 0; m < iaos_to_atoms_.size(); m++) {
            if (iaos_to_atoms_[m] == A) {
                vec.push_back(m);
            }
        }
        minao_inds.push_back(vec);
    }

    std::vector<std::pair<int, int> > rot_inds;
    for (int ind = 0; ind < ranges_.size() - 1; ind++) {
        int start = ranges_[ind];
        int stop = ranges_[ind + 1];
        for (int i = start; i < stop; i++) {
            for (int j = start; j < i; j++) {
                rot_inds.push_back(std::pair<int, int>(i, j));
            }
        }
    }

    auto L = linalg::triplet(C_, S_, A_, true, false, false);
    L->set_name("L");

    auto ret1 = IBOLocalizer::localize_task(L, minao_inds, rot_inds, convergence_, maxiter_, power_);
    L = ret1["L"];
    auto U = ret1["U"];

    if (use_stars_) {
        auto Q = orbital_charges(L);
        double** Qp = Q->pointer();
        int nocc = Q->colspi()[0];
        int natom = Q->rowspi()[0];

        std::vector<int> pi_orbs;
        for (int i = 0; i < nocc; i++) {
            std::vector<double> Qs;
            for (int A = 0; A < natom; A++) {
                Qs.push_back(std::fabs(Qp[A][i]));
            }
            std::sort(Qs.begin(), Qs.end(), std::greater<double>());
            double Qtot = 0.0;
            for (int A = 0; A < natom && A < 2; A++) {
                Qtot += Qs[A];
            }
            if (Qtot < stars_completeness_) {
                pi_orbs.push_back(i);
            }
        }

        std::vector<std::pair<int, int> > rot_inds2;
        for (int iind = 0; iind < pi_orbs.size(); iind++) {
            for (int jind = 0; jind < iind; jind++) {
                rot_inds2.push_back(std::pair<int, int>(pi_orbs[iind], pi_orbs[jind]));
            }
        }

        std::vector<std::vector<int> > minao_inds2;
        for (int Aind = 0; Aind < stars_.size(); Aind++) {
            int A = -1;
            for (int A2 = 0; A2 < true_atoms_.size(); A2++) {
                if (stars_[Aind] == true_atoms_[A2]) {
                    A = A2;
                    break;
                }
            }
            if (A == -1) continue;
            std::vector<int> vec;
            for (int m = 0; m < iaos_to_atoms_.size(); m++) {
                if (iaos_to_atoms_[m] == A) {
                    vec.push_back(m);
                }
            }
            minao_inds2.push_back(vec);
        }

        outfile->Printf("    *** Stars Procedure ***\n\n");
        outfile->Printf("    Pi Completeness = %11.3f\n", stars_completeness_);
        outfile->Printf("    Number of Pis   = %11zu\n", pi_orbs.size());
        outfile->Printf("    Number of Stars = %11zu\n", stars_.size());
        outfile->Printf("    Star Centers: ");
        for (int ind = 0; ind < stars_.size(); ind++) {
            outfile->Printf("%3d ", stars_[ind] + 1);
        }
        outfile->Printf("\n\n");

        auto ret2 = IBOLocalizer::localize_task(L, minao_inds2, rot_inds2, convergence_, maxiter_, power_);
        L = ret2["L"];
        auto U3 = ret2["U"];
        U = linalg::doublet(U, U3, false, false);

        auto ret3 = IBOLocalizer::localize_task(L, minao_inds, rot_inds, convergence_, maxiter_, power_);
        L = ret3["L"];
        auto U4 = ret3["U"];
        U = linalg::doublet(U, U4, false, false);

        // => Analysis <= //

        Q = orbital_charges(L);
        Qp = Q->pointer();

        pi_orbs.clear();
        for (int i = 0; i < nocc; i++) {
            std::vector<double> Qs;
            for (int A = 0; A < natom; A++) {
                Qs.push_back(std::fabs(Qp[A][i]));
            }
            std::sort(Qs.begin(), Qs.end(), std::greater<double>());
            double Qtot = 0.0;
            for (int A = 0; A < natom && A < 2; A++) {
                Qtot += Qs[A];
            }
            if (Qtot < stars_completeness_) {
                pi_orbs.push_back(i);
            }
        }

        std::vector<int> centers;
        for (int i2 = 0; i2 < pi_orbs.size(); i2++) {
            int i = pi_orbs[i2];
            int ind = 0;
            for (int A = 0; A < natom; A++) {
                if (std::fabs(Qp[A][i]) >= std::fabs(Qp[ind][i])) {
                    ind = A;
                }
            }
            centers.push_back(ind);
        }
        std::sort(centers.begin(), centers.end());

        outfile->Printf("    *** Stars Analysis ***\n\n");
        outfile->Printf("    Pi Centers: ");
        for (int ind = 0; ind < centers.size(); ind++) {
            outfile->Printf("%3d ", centers[ind] + 1);
        }
        outfile->Printf("\n\n");
        Q_ = Q;
    }

    auto Focc2 = linalg::triplet(U, Focc_, U, true, false, false);
    auto U2 = IBOLocalizer::reorder_orbitals(Focc2, ranges_);

    // Compute Orbital Charges
    L = linalg::doublet(U2, L, true, false);
    Q_ = orbital_charges(L);

    // Form new U and L
    U_ = linalg::doublet(U, U2, false, false);
    L_ = linalg::doublet(C_, U_, false, false);
}

std::shared_ptr<Matrix> IBOLocalizer::orbital_charges(std::shared_ptr<Matrix> L) {
    double** Lp = L->pointer();
    int nocc = L->rowspi()[0];
    int nmin = L->colspi()[0];
    int natom = true_atoms_.size();

    auto Q = std::make_shared<Matrix>("Q", natom, nocc);
    double** Qp = Q->pointer();

    for (int i = 0; i < nocc; i++) {
        for (int m = 0; m < nmin; m++) {
            Qp[iaos_to_atoms_[m]][i] += Lp[i][m] * Lp[i][m];
        }
    }

    return Q;
}

void IBOLocalizer::print_charges(double scale) {
    if (!A_) build_iaos();

    auto L = linalg::triplet(C_, S_, A_, true, false, false);

    int nocc = L->rowspi()[0];
    int natom = true_atoms_.size();

    auto Q = orbital_charges(L);
    double** Qp = Q->pointer();

    auto N = Vector("N", natom);
    double* Np = N.pointer();

    for (int A = 0; A < natom; A++) {
        for (int i = 0; i < nocc; i++) {
            Np[A] += Qp[A][i];
        }
    }

    std::shared_ptr<Molecule> mol = minao_->molecule();

    outfile->Printf("   > Atomic Charges <\n\n");
    outfile->Printf("    %4s %3s %11s %11s %11s\n", "N", "Z", "Nuclear", "Electronic", "Atomic");
    double Ztot = 0.0;
    double Qtot = 0.0;
    for (int A = 0; A < natom; A++) {
        int Afull = true_atoms_[A];
        double Z = mol->Z(Afull);
        double Q = -scale * Np[A];
        outfile->Printf("    %4d %3s %11.3E %11.3E %11.3E\n", Afull + 1, mol->symbol(Afull).c_str(), Z, Q, Z + Q);
        Ztot += Z;
        Qtot += Q;
    }
    outfile->Printf("    %8s %11.3E %11.3E %11.3E\n", "Total", Ztot, Qtot, Ztot + Qtot);
    outfile->Printf("\n");

    outfile->Printf("    True Molecular Charge: %11.3E\n", (double)mol->molecular_charge());
    outfile->Printf("    IBO  Molecular Charge: %11.3E\n", Ztot + Qtot);
    outfile->Printf("    IBO  Error:            %11.3E\n", Ztot + Qtot - (double)mol->molecular_charge());
    outfile->Printf("\n");
}

}  // Namespace psi
