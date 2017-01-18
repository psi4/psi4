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
#include "psi4/libqt/qt.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"

;
using namespace psi;

namespace psi {

Localizer::Localizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C) :
    primary_(primary), C_(C)
{
    if (C->nirrep() != 1) {
        throw PSIEXCEPTION("Localizer: C matrix is not C1");
    }
    if (C->rowspi()[0] != primary->nbf()) {
        throw PSIEXCEPTION("Localizer: C matrix does not match basis");
    }
    common_init();
}
Localizer::~Localizer()
{
}
void Localizer::common_init()
{
    print_ = 0;
    debug_ = 0;
    bench_ = 0;
    convergence_ = 1.0E-8;
    maxiter_ = 50;
    converged_ = false;
}
std::shared_ptr<Localizer> Localizer::build(const std::string& type, std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C, Options& options)
{
    std::shared_ptr<Localizer> local;

    if (type == "BOYS") {
        BoysLocalizer* l = new BoysLocalizer(primary,C);
        local = std::shared_ptr<Localizer>(l);
    } else if (type == "PIPEK_MEZEY") {
        PMLocalizer* l = new PMLocalizer(primary,C);
        local = std::shared_ptr<Localizer>(l);
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
std::shared_ptr<Localizer> Localizer::build(const std::string& type, std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C)
{
    return Localizer::build(type, primary, C, Process::environment.options);
}
std::shared_ptr<Localizer> Localizer::build(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C, Options& options)
{
    return Localizer::build(options.get_str("LOCAL_TYPE"), primary, C, options);
}
std::shared_ptr<Matrix> Localizer::fock_update(std::shared_ptr<Matrix> Fc)
{
    if (!L_ || !U_) {
        throw PSIEXCEPTION("Localizer: run compute() first");
    }

    int nso = L_->rowspi()[0];
    int nmo = L_->colspi()[0];

    if (nmo < 1) return Fc;

    std::shared_ptr<Matrix> Fl = Matrix::triplet(U_, Fc, U_, true, false, false);
    double** Fp = Fl->pointer();
    double** Lp = L_->pointer();
    double** Up = U_->pointer();

    std::vector<std::pair<double, int> > order;
    for (int i = 0; i < nmo; i++) {
        order.push_back(std::pair<double, int>(Fp[i][i],i));
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
        C_DCOPY(nso,&L2p[0][order[i].second],nmo,&Lp[0][i],nmo);
        C_DCOPY(nmo,&U2p[0][order[i].second],nmo,&Up[0][i],nmo);
    }

    return Fl;
}

BoysLocalizer::BoysLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C) :
    Localizer(primary,C)
{
    common_init();
}
BoysLocalizer::~BoysLocalizer()
{
}
void BoysLocalizer::common_init()
{
}
void BoysLocalizer::print_header() const
{
    outfile->Printf( "  ==> Boys Localizer <==\n\n");
    outfile->Printf( "    Convergence = %11.3E\n", convergence_);
    outfile->Printf( "    Maxiter     = %11d\n", maxiter_);
    outfile->Printf( "\n");

}
void BoysLocalizer::localize()
{
    print_header();

    // => Sizing <= //

    int nso = C_->rowspi()[0];
    int nmo = C_->colspi()[0];

    // => Dipole Integrals <= //

    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_));
    std::shared_ptr<OneBodyAOInt> Dint(fact->ao_dipole());
    std::vector<std::shared_ptr<Matrix> > D;
    for (int xyz = 0; xyz < 3; xyz++) {
        D.push_back(std::shared_ptr<Matrix>(new Matrix("D", nso, nso)));
    }
    Dint->compute(D);
    Dint.reset();
    fact.reset();

    std::vector<std::shared_ptr<Matrix> > Dmo;
    std::shared_ptr<Matrix> T(new Matrix("T",nso,nmo));
    for (int xyz = 0; xyz < 3; xyz++) {
        Dmo.push_back(std::shared_ptr<Matrix>(new Matrix("D", nmo, nmo)));
        C_DGEMM('N','N',nso,nmo,nso,1.0,D[xyz]->pointer()[0],nso,C_->pointer()[0],nmo,0.0,T->pointer()[0],nmo);
        C_DGEMM('T','N',nmo,nmo,nso,1.0,C_->pointer()[0],nmo,T->pointer()[0],nmo,0.0,Dmo[xyz]->pointer()[0],nmo);
    }
    D.clear();
    T.reset();

    // => Targets <= //

    L_ = std::shared_ptr<Matrix>(new Matrix("L",nso,nmo));
    U_ = std::shared_ptr<Matrix>(new Matrix("U",nmo,nmo));
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
        metric += C_DDOT(nmo,Dp[xyz][0],nmo+1,Dp[xyz][0],nmo+1);
    }
    double old_metric = metric;

    // => Iteration Print <= //
    outfile->Printf( "    Iteration %24s %14s\n", "Metric", "Residual");
    outfile->Printf( "    @Boys %4d %24.16E %14s\n", 0, metric, "-");

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

        for (int i2 = 0; i2 < nmo-1; i2++) {
            for (int j2 = i2+1; j2 < nmo; j2++) {

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
                if (fabs(theta) < 1.0E-8) {
                    double O0 = 0.0;
                    double O1 = 0.0;;
                    for (int xyz = 0; xyz < 3; xyz++) {
                        double** Ak = Dp[xyz];
                        O0 += Ak[i][j] * Ak[i][j];
                        O1 += 0.25 * (Ak[j][j] - Ak[i][i]) * (Ak[j][j] - Ak[i][i]);
                    }
                    if (O1 < O0) {
                        theta = M_PI / 4.0;
                        if (debug_ > 3) {
                            outfile->Printf( "@Break\n");
                        }
                    }
                }

                // Givens rotation
                cc = cos(theta);
                ss = sin(theta);

                if (debug_ > 3) {
                    outfile->Printf( "@Rotation, i = %4d, j = %4d, Theta = %24.16E\n", i, j,theta);
                    outfile->Printf( "@Info, a = %24.16E, b = %24.16E, c = %24.16E\n", a, b, c);
                }

                // > Apply the rotation < //

                // rows and columns of A^k
                for (int xyz = 0; xyz < 3; xyz++) {
                    double** Ak = Dp[xyz];
                    C_DROT(nmo,&Ak[i][0],1,&Ak[j][0],1,cc,ss);
                    C_DROT(nmo,&Ak[0][i],nmo,&Ak[0][j],nmo,cc,ss);
                }

                // Q
                C_DROT(nmo,Up[i],1,Up[j],1,cc,ss);
            }
        }

        // => Metric <= //

        metric = 0.0;
        for (int xyz = 0; xyz < 3; xyz++) {
            metric += C_DDOT(nmo,Dp[xyz][0],nmo+1,Dp[xyz][0],nmo+1);
        }

        double conv = fabs(metric - old_metric) / fabs(old_metric);
        old_metric = metric;

        // => Iteration Print <= //

        outfile->Printf( "    @Boys %4d %24.16E %14.6E\n", iter, metric, conv);

        // => Convergence Check <= //

        if (conv < convergence_) {
            converged_ = true;
            break;
        }

    }

    outfile->Printf( "\n");
    if (converged_) {
        outfile->Printf( "    Boys Localizer converged.\n\n");
    } else {
        outfile->Printf( "    Boys Localizer failed.\n\n");
    }

    U_->transpose_this();
    C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Up[0],nmo,0.0,Lp[0],nmo);
}

PMLocalizer::PMLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C) :
    Localizer(primary,C)
{
    common_init();
}
PMLocalizer::~PMLocalizer()
{
}
void PMLocalizer::common_init()
{
}
void PMLocalizer::print_header() const
{
    outfile->Printf( "  ==> Pipek-Mezey Localizer <==\n\n");
    outfile->Printf( "    Convergence = %11.3E\n", convergence_);
    outfile->Printf( "    Maxiter     = %11d\n", maxiter_);
    outfile->Printf( "\n");

}
void PMLocalizer::localize()
{
    print_header();

    // => Sizing <= //

    int nso = C_->rowspi()[0];
    int nmo = C_->colspi()[0];
    int nA  = primary_->molecule()->natom();

    // => Overlap Integrals <= //

    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_));
    std::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());
    std::shared_ptr<Matrix> S(new Matrix("S",nso,nso));
    Sint->compute(S);
    Sint.reset();
    fact.reset();
    double** Sp = S->pointer();

    // => Targets <= //

    L_ = std::shared_ptr<Matrix>(new Matrix("L",nso,nmo));
    U_ = std::shared_ptr<Matrix>(new Matrix("U",nmo,nmo));
    L_->copy(C_);
    U_->identity();
    converged_ = false;

    if (nmo < 1) return;

    // => Pointers <= //

    double** Lp = L_->pointer();
    double** Up = U_->pointer();

    // => LS product (avoids GEMV) <= //

    std::shared_ptr<Matrix> LS(new Matrix("LS",nso,nmo));
    double** LSp = LS->pointer();
    C_DGEMM('N','N',nso,nmo,nso,1.0,Sp[0],nso,Lp[0],nmo,0.0,LSp[0],nmo);

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
            int nm = Astarts[A+1] - Astarts[A];
            int off = Astarts[A];
            double PA = C_DDOT(nm,&LSp[off][i],nmo,&Lp[off][i],nmo);
            metric += PA * PA;
        }
    }
    double old_metric = metric;

    // => Iteration Print <= //
    outfile->Printf( "    Iteration %24s %14s\n", "Metric", "Residual");
    outfile->Printf( "    @PM %4d %24.16E %14s\n", 0, metric, "-");

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

        for (int i2 = 0; i2 < nmo-1; i2++) {
            for (int j2 = i2+1; j2 < nmo; j2++) {

                int i = order2[i2];
                int j = order2[j2];

                // > Compute the rotation < //

                // H elements
                a = 0.0;
                b = 0.0;
                c = 0.0;
                for (int A = 0; A < nA; A++) {

                    int nm = Astarts[A+1] - Astarts[A];
                    int off = Astarts[A];
                    Aii = C_DDOT(nm,&LSp[off][i],nmo,&Lp[off][i],nmo);
                    Ajj = C_DDOT(nm,&LSp[off][j],nmo,&Lp[off][j],nmo);
                    Aij = 0.5 * C_DDOT(nm,&LSp[off][i],nmo,&Lp[off][j],nmo) +
                          0.5 * C_DDOT(nm,&LSp[off][j],nmo,&Lp[off][i],nmo);


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
                if (fabs(theta) < 1.0E-8) {
                    double O0 = 0.0;
                    double O1 = 0.0;;
                    for (int A = 0; A < nA; A++) {
                        int nm = Astarts[A+1] - Astarts[A];
                        int off = Astarts[A];
                        Aii = C_DDOT(nm,&LSp[off][i],nmo,&Lp[off][i],nmo);
                        Ajj = C_DDOT(nm,&LSp[off][j],nmo,&Lp[off][j],nmo);
                        Aij = 0.5 * C_DDOT(nm,&LSp[off][i],nmo,&Lp[off][j],nmo) +
                              0.5 * C_DDOT(nm,&LSp[off][j],nmo,&Lp[off][i],nmo);
                        O0 += Aij * Aij;
                        O1 += 0.25 * (Ajj - Aii) * (Ajj - Aii);
                    }
                    if (O1 < O0) {
                        theta = M_PI / 4.0;
                        if (debug_ > 3) {
                            outfile->Printf( "@Break\n");
                        }
                    }
                }

                // Givens rotation
                cc = cos(theta);
                ss = sin(theta);

                if (debug_ > 3) {
                    outfile->Printf( "@Rotation, i = %4d, j = %4d, Theta = %24.16E\n", i, j,theta);
                    outfile->Printf( "@Info, a = %24.16E, b = %24.16E, c = %24.16E\n", a, b, c);
                }

                // > Apply the rotation < //

                // columns of LS and L
                C_DROT(nso,&LSp[0][i],nmo,&LSp[0][j],nmo,cc,ss);
                C_DROT(nso,&Lp[0][i],nmo,&Lp[0][j],nmo,cc,ss);

                // Q
                C_DROT(nmo,Up[i],1,Up[j],1,cc,ss);
            }
        }

        // => Metric <= //

        metric = 0.0;
        for (int i = 0; i < nmo; i++) {
            for (int A = 0; A < nA; A++) {
                int nm = Astarts[A+1] - Astarts[A];
                int off = Astarts[A];
                double PA = C_DDOT(nm,&LSp[off][i],nmo,&Lp[off][i],nmo);
                metric += PA * PA;
            }
        }

        double conv = fabs(metric - old_metric) / fabs(old_metric);
        old_metric = metric;

        // => Iteration Print <= //

        outfile->Printf( "    @PM %4d %24.16E %14.6E\n", iter, metric, conv);

        // => Convergence Check <= //

        if (conv < convergence_) {
            converged_ = true;
            break;
        }

    }

    outfile->Printf( "\n");
    if (converged_) {
        outfile->Printf( "    PM Localizer converged.\n\n");
    } else {
        outfile->Printf( "    PM Localizer failed.\n\n");
    }

    U_->transpose_this();
}

} // Namespace psi
