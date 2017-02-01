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
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/fisapt/local2.h"



namespace psi {

namespace fisapt {

IBOLocalizer2::IBOLocalizer2(
    std::shared_ptr<BasisSet> primary,
    std::shared_ptr<BasisSet> minao,
    std::shared_ptr<Matrix> C) :
    primary_(primary),
    minao_(minao),
    C_(C)
{
    if (C->nirrep() != 1) {
        throw PSIEXCEPTION("Localizer: C matrix is not C1");
    }
    if (C->rowspi()[0] != primary->nbf()) {
        throw PSIEXCEPTION("Localizer: C matrix does not match basis");
    }
    common_init();
}
IBOLocalizer2::~IBOLocalizer2()
{
}
void IBOLocalizer2::common_init()
{
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
std::shared_ptr<IBOLocalizer2> IBOLocalizer2::build(
    std::shared_ptr<BasisSet> primary,
    std::shared_ptr<BasisSet> minao,
    std::shared_ptr<Matrix> C,
    Options& options)
{
//    Options& options = Process::environment.options;

    std::shared_ptr<IBOLocalizer2> local(new IBOLocalizer2(primary, minao, C));

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
void IBOLocalizer2::print_header() const
{
    outfile->Printf( "  ==> IBO Localizer 2 <==\n\n");
    outfile->Printf( "    MinAO Basis = %11s\n", minao_->name().c_str());
    outfile->Printf( "    Use Ghosts  = %11s\n", (use_ghosts_ ? "TRUE" : "FALSE"));
    outfile->Printf( "    Use Stars   = %11s\n", (use_stars_  ? "TRUE" : "FALSE"));
    outfile->Printf( "    Condition   = %11.3E\n", condition_);
    outfile->Printf( "    Power       = %11d\n", power_);
    outfile->Printf( "    Convergence = %11.3E\n", convergence_);
    outfile->Printf( "    Maxiter     = %11d\n", maxiter_);
    outfile->Printf( "\n");
    //fflush(outfile);
}
void IBOLocalizer2::build_iaos()
{
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

    std::shared_ptr<IntegralFactory> fact11(new IntegralFactory(primary_,primary_));
    std::shared_ptr<IntegralFactory> fact12(new IntegralFactory(primary_,minao_));
    std::shared_ptr<IntegralFactory> fact22(new IntegralFactory(minao_,minao_));

    std::shared_ptr<OneBodyAOInt> ints11(fact11->ao_overlap());
    std::shared_ptr<OneBodyAOInt> ints12(fact12->ao_overlap());
    std::shared_ptr<OneBodyAOInt> ints22(fact22->ao_overlap());

    std::shared_ptr<Matrix> S11(new Matrix("S11", primary_->nbf(), primary_->nbf()));
    std::shared_ptr<Matrix> S12f(new Matrix("S12f", primary_->nbf(), minao_->nbf()));
    std::shared_ptr<Matrix> S22f(new Matrix("S22f", minao_->nbf(), minao_->nbf()));

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

    std::shared_ptr<Matrix> S12(new Matrix("S12", primary_->nbf(), true_iaos_.size()));
    std::shared_ptr<Matrix> S22(new Matrix("S22", true_iaos_.size(), true_iaos_.size()));

    double** S12p  = S12->pointer();
    double** S12fp = S12f->pointer();
    for (int m = 0; m < primary_->nbf(); m++) {
        for (int p = 0; p < true_iaos_.size(); p++) {
            S12p[m][p] = S12fp[m][true_iaos_[p]];
        }
    }

    double** S22p  = S22->pointer();
    double** S22fp = S22f->pointer();
    for (int p = 0; p < true_iaos_.size(); p++) {
        for (int q = 0; q < true_iaos_.size(); q++) {
            S22p[p][q] = S22fp[true_iaos_[p]][true_iaos_[q]];

        }
    }

    // => Metric Inverses <= //

    std::shared_ptr<Matrix>S11_m12(S11->clone());
    std::shared_ptr<Matrix>S22_m12(S22->clone());
    S11_m12->copy(S11);
    S22_m12->copy(S22);
    S11_m12->power(-1.0/2.0, condition_);
    S22_m12->power(-1.0/2.0, condition_);

    // => Tilde C <= //

    std::shared_ptr<Matrix> C = C_;
    std::shared_ptr<Matrix> T1 = Matrix::doublet(S22_m12, S12, false, true);
    std::shared_ptr<Matrix> T2 = Matrix::doublet(S11_m12, Matrix::triplet(T1, T1, C, true, false, false), false, false);
    std::shared_ptr<Matrix> T3 = Matrix::doublet(T2, T2, true, false);
    T3->power(-1.0/2.0, condition_);
    std::shared_ptr<Matrix> Ctilde = Matrix::triplet(S11_m12, T2, T3, false, false, false);

    // => D and Tilde D <= //

    std::shared_ptr<Matrix> D = Matrix::doublet(C, C, false, true);
    std::shared_ptr<Matrix> Dtilde = Matrix::doublet(Ctilde, Ctilde, false, true);

    // => A (Before Orthogonalization) <= //

    std::shared_ptr<Matrix> DSDtilde = Matrix::triplet(D, S11, Dtilde,false, false, false);
    DSDtilde->scale(2.0);

    std::shared_ptr<Matrix> L = Matrix::doublet(S11_m12, S11_m12, false, false); // TODO: Possibly Unstable
    L->add(DSDtilde);
    L->subtract(D);
    L->subtract(Dtilde);

    std::shared_ptr<Matrix> AN = Matrix::doublet(L, S12, false, false);

    // => A (After Orthogonalization) <= //

    std::shared_ptr<Matrix> V = Matrix::triplet(AN, S11, AN, true, false, false);
    V->power(-1.0/2.0, condition_);

    std::shared_ptr<Matrix> A = Matrix::doublet(AN, V, false, false);

    // => Assignment <= //

    S_ = S11;
    A_ = A;
}
std::map<std::string, std::shared_ptr<Matrix> > IBOLocalizer2::localize_task(
    std::shared_ptr<Matrix> L,
    const std::vector<std::vector<int> >& minao_inds,
    const std::vector<std::pair<int, int> >& rot_inds,
    double convergence,
    int maxiter,
    int power
    )
{
    int nmin = L->colspi()[0];
    int nocc = L->rowspi()[0];

    std::shared_ptr<Matrix> L2(L->clone());
    L2->copy(L);
    double** Lp = L2->pointer();

    std::shared_ptr<Matrix> U(new Matrix("U", nocc, nocc));
    U->identity();
    double** Up = U->pointer();

    bool converged = false;

    if (power != 2 && power != 4) throw PSIEXCEPTION("IAO: Invalid metric power.");

    outfile->Printf( "    @IBO %4s: %24s %14s\n", "Iter", "Metric", "Gradient");

    for (int iter = 1; iter <= maxiter; iter++) {

        double metric = 0.0;
        for (int i = 0; i < nocc; i++) {
            for (int A = 0; A < minao_inds.size(); A++) {
                double Lval = 0.0;
                for (int m = 0; m < minao_inds[A].size(); m++) {
                    int mind = minao_inds[A][m];
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
            for (int A = 0; A < minao_inds.size(); A++) {
                double Qii = 0.0;
                double Qij = 0.0;
                double Qjj = 0.0;
                for (int m = 0; m < minao_inds[A].size(); m++) {
                    int mind = minao_inds[A][m];
                    Qii += Lp[i][mind] * Lp[i][mind];
                    Qij += Lp[i][mind] * Lp[j][mind];
                    Qjj += Lp[j][mind] * Lp[j][mind];
                }
                if (power == 2) {
                    Aij += 4.0 * Qij * Qij - (Qii - Qjj) * (Qii - Qjj);
                    Bij += 4.0 * Qij * (Qii - Qjj);
                } else {
                    Aij += (-1.0) * Qii * Qii * Qii * Qii - Qjj * Qjj * Qjj * Qjj + 6.0 * (Qii * Qii + Qjj * Qjj) * Qij * Qij + Qii * Qii * Qii * Qjj + Qii * Qjj * Qjj * Qjj;
                    Bij += 4.0 * Qij * (Qii * Qii * Qii - Qjj * Qjj * Qjj);
                }
            }

            double phi = 0.25 * atan2(Bij, -Aij);
            double c = cos(phi);
            double s = sin(phi);

            C_DROT(nmin,Lp[i],1,Lp[j],1,c,s);
            C_DROT(nocc,Up[i],1,Up[j],1,c,s);

            gradient += Bij * Bij;

        }
        gradient = sqrt(gradient);

        outfile->Printf( "    @IBO %4d: %24.16E %14.6E\n", iter, metric, gradient);

        if (gradient < convergence) {
            converged = true;
            break;
        }

    }

    outfile->Printf( "\n");
    if (converged) {
        outfile->Printf( "    IBO Localizer 2 converged.\n\n");
    } else {
        outfile->Printf( "    IBO Localizer 2 failed.\n\n");
    }

    U->transpose_this();

    std::map<std::string, std::shared_ptr<Matrix> > ret;
    ret["U"] = U;
    ret["L"] = L2;

    ret["U"]->set_name("U");
    ret["L"]->set_name("L");

    return ret;
}
std::shared_ptr<Matrix> IBOLocalizer2::reorder_orbitals(
    std::shared_ptr<Matrix> F,
    const std::vector<int>& ranges)
{
    int nmo = F->rowspi()[0];
    double** Fp = F->pointer();

    std::shared_ptr<Matrix> U(new Matrix("U", nmo, nmo));
    double** Up = U->pointer();

    for (int ind = 0; ind < ranges.size() - 1; ind++) {
        int start = ranges[ind];
        int stop = ranges[ind+1];
        std::vector<std::pair<double,int> > fvals;
        for (int i = start; i < stop; i++) {
            fvals.push_back(std::pair<double,int>(Fp[i][i],i));
        }
        std::sort(fvals.begin(),fvals.end());
        for (int i = start; i < stop; i++) {
            Up[i][fvals[i-start].second] = 1.0;
        }
    }

    return U;
}
std::map<std::string, std::shared_ptr<Matrix> > IBOLocalizer2::localize(
        std::shared_ptr<Matrix> Cocc,
        std::shared_ptr<Matrix> Focc,
        const std::vector<int>& ranges2
        )
{
    if (!A_) build_iaos();

    std::vector<int> ranges = ranges2;
    if (!ranges.size()) {
        ranges.push_back(0);
        ranges.push_back(Cocc->colspi()[0]);
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

    std::vector<std::pair<int,int> > rot_inds;
    for (int ind = 0; ind < ranges.size() - 1; ind++) {
        int start = ranges[ind];
        int stop  = ranges[ind+1];
        for (int i = start; i < stop; i++) {
            for (int j = start; j < i; j++) {
                rot_inds.push_back(std::pair<int,int>(i,j));
            }
        }
    }

    std::shared_ptr<Matrix> L = Matrix::triplet(Cocc,S_,A_,true,false,false);
    L->set_name("L");

    std::map<std::string, std::shared_ptr<Matrix> > ret1 = IBOLocalizer2::localize_task(L,minao_inds,rot_inds,convergence_,maxiter_,power_);
    L = ret1["L"];
    std::shared_ptr<Matrix> U = ret1["U"];

    if (use_stars_) {
        std::shared_ptr<Matrix> Q = orbital_charges(L);
        double** Qp = Q->pointer();
        int nocc  = Q->colspi()[0];
        int natom = Q->rowspi()[0];

        std::vector<int> pi_orbs;
        for (int i = 0; i < nocc; i++) {
            std::vector<double> Qs;
            for (int A = 0; A < natom; A++) {
                Qs.push_back(fabs(Qp[A][i]));
            }
            std::sort(Qs.begin(),Qs.end(),std::greater<double>());
            double Qtot = 0.0;
            for (int A = 0; A < natom && A < 2; A++) {
                Qtot += Qs[A];
            }
            if (Qtot < stars_completeness_) {
                pi_orbs.push_back(i);
            }
        }

        std::vector<std::pair<int,int> > rot_inds2;
        for (int iind = 0; iind < pi_orbs.size(); iind++) {
            for (int jind = 0; jind < iind; jind++) {
                rot_inds2.push_back(std::pair<int,int>(pi_orbs[iind],pi_orbs[jind]));
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

        outfile->Printf( "    *** Stars Procedure ***\n\n");
        outfile->Printf( "    Pi Completeness = %11.3f\n", stars_completeness_);
        outfile->Printf( "    Number of Pis   = %11zu\n", pi_orbs.size());
        outfile->Printf( "    Number of Stars = %11zu\n", stars_.size());
        outfile->Printf( "    Star Centers: ");
        for (int ind = 0; ind < stars_.size(); ind++) {
            outfile->Printf( "%3d ", stars_[ind]+1);
        }
        outfile->Printf( "\n\n");

        std::map<std::string, std::shared_ptr<Matrix> > ret2 = IBOLocalizer2::localize_task(L,minao_inds2,rot_inds2,convergence_,maxiter_,power_);
        L = ret2["L"];
        std::shared_ptr<Matrix> U3 = ret2["U"];
        U = Matrix::doublet(U,U3,false,false);

        std::map<std::string, std::shared_ptr<Matrix> > ret3 = IBOLocalizer2::localize_task(L,minao_inds,rot_inds,convergence_,maxiter_,power_);
        L = ret3["L"];
        std::shared_ptr<Matrix> U4 = ret3["U"];
        U = Matrix::doublet(U,U4,false,false);

        // => Analysis <= //

        Q = orbital_charges(L);
        Qp = Q->pointer();

        pi_orbs.clear();
        for (int i = 0; i < nocc; i++) {
            std::vector<double> Qs;
            for (int A = 0; A < natom; A++) {
                Qs.push_back(fabs(Qp[A][i]));
            }
            std::sort(Qs.begin(),Qs.end(),std::greater<double>());
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
                if (fabs(Qp[A][i]) >= fabs(Qp[ind][i])) {
                    ind = A;
                }
            }
            centers.push_back(ind);
        }
        std::sort(centers.begin(), centers.end());

        outfile->Printf( "    *** Stars Analysis ***\n\n");
        outfile->Printf( "    Pi Centers: ");
        for (int ind = 0; ind < centers.size(); ind++) {
            outfile->Printf( "%3d ", centers[ind]+1);
        }
        outfile->Printf( "\n\n");
    }

    std::shared_ptr<Matrix> Focc2 = Matrix::triplet(U,Focc,U,true,false,false);
    std::shared_ptr<Matrix> U2 = IBOLocalizer2::reorder_orbitals(Focc2, ranges);

    std::shared_ptr<Matrix> Uocc3 = Matrix::doublet(U,U2,false,false);
    std::shared_ptr<Matrix> Focc3 = Matrix::triplet(Uocc3,Focc,Uocc3,true,false,false);
    std::shared_ptr<Matrix> Locc3 = Matrix::doublet(Cocc,Uocc3,false,false);

    L = Matrix::doublet(U2,L,true,false);
    std::shared_ptr<Matrix> Q = orbital_charges(L);

    std::map<std::string, std::shared_ptr<Matrix> > ret;
    ret["L"] = Locc3;
    ret["U"] = Uocc3;
    ret["F"] = Focc3;
    ret["Q"] = Q;

    ret["L"]->set_name("L");
    ret["U"]->set_name("U");
    ret["F"]->set_name("F");
    ret["Q"]->set_name("Q");

    return ret;
}
std::shared_ptr<Matrix> IBOLocalizer2::orbital_charges(
    std::shared_ptr<Matrix> L)
{
    double** Lp = L->pointer();
    int nocc = L->rowspi()[0];
    int nmin = L->colspi()[0];
    int natom = true_atoms_.size();

    std::shared_ptr<Matrix> Q(new Matrix("Q", natom, nocc));
    double** Qp = Q->pointer();

    for (int i = 0; i < nocc; i++) {
        for (int m = 0; m < nmin; m++) {
            Qp[iaos_to_atoms_[m]][i] += Lp[i][m] * Lp[i][m];
        }
    }

    return Q;
}
void IBOLocalizer2::print_charges(double scale)
{
    if (!A_) build_iaos();

    std::shared_ptr<Matrix> L = Matrix::triplet(C_, S_, A_, true, false, false);

    int nocc = L->rowspi()[0];
    int natom = true_atoms_.size();

    std::shared_ptr<Matrix> Q = orbital_charges(L);
    double** Qp = Q->pointer();

    std::shared_ptr<Vector> N(new Vector("N", natom));
    double* Np = N->pointer();

    for (int A = 0; A < natom; A++) {
        for (int i = 0; i < nocc; i++) {
            Np[A] += Qp[A][i];
        }
    }

    std::shared_ptr<Molecule> mol = minao_->molecule();

    outfile->Printf("   > Atomic Charges <\n\n");
    outfile->Printf("    %4s %3s %11s %11s %11s\n",
        "N", "Z", "Nuclear", "Electronic", "Atomic");
    double Ztot = 0.0;
    double Qtot = 0.0;
    for (int A = 0; A < natom; A++) {
        int Afull = true_atoms_[A];
        double Z = mol->Z(Afull);
        double Q = -scale * Np[A];
        outfile->Printf("    %4d %3s %11.3E %11.3E %11.3E\n",
            Afull+1, mol->symbol(Afull).c_str(), Z, Q, Z + Q);
        Ztot += Z;
        Qtot += Q;
    }
    outfile->Printf("    %8s %11.3E %11.3E %11.3E\n",
            "Total", Ztot, Qtot, Ztot + Qtot);
    outfile->Printf("\n");

    outfile->Printf("    True Molecular Charge: %11.3E\n", (double) mol->molecular_charge());
    outfile->Printf("    IBO  Molecular Charge: %11.3E\n", Ztot + Qtot);
    outfile->Printf("    IBO  Error:            %11.3E\n", Ztot + Qtot - (double) mol->molecular_charge());
    outfile->Printf("\n");
}

} // Namespace fisapt

} // Namespace psi
