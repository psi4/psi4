/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>

namespace psi {

ExternalPotential::ExternalPotential() : debug_(0), print_(1) {}

ExternalPotential::~ExternalPotential() {}

void ExternalPotential::clear() {
    charges_.clear();
    bases_.clear();
    exchange_bases_.clear();
}

void ExternalPotential::addCharge(double Z, double x, double y, double z) {
    charges_.push_back(std::make_tuple(Z, x, y, z));
}

void ExternalPotential::addBasis(std::shared_ptr<BasisSet> basis, SharedVector coefs) {
    bases_.push_back(std::make_pair(basis, coefs));
}

void ExternalPotential::addExchangeBasis(std::shared_ptr<BasisSet> basis, SharedVector coefs) {
    exchange_bases_.push_back(std::make_pair(basis, coefs));
}

void ExternalPotential::addMultipoles(std::shared_ptr<BasisSet> basis, SharedVector coefs) {
    // This is where the addBasis and addMultipoles routines differ. Add basis is for a general basis set that has been
    // fit to a given density, so the normalization is included in the fit.  When specifying some arbitrary multipoles
    // from a classical force field, we want each function |A) defining the multipole to be normalized to 1.  However,
    // standard basis sets are normalized such that (A|A) = 1 so the user would not get the expected pertubation.
    // Instead we work out the renormalization factor needed to convert the condition (A|A) = 1 into ||A|) = 1.
    std::vector<double> scale_facs;
    for (int Q = 0; Q < basis->nshell(); Q++) {
        const auto &shell = basis->shell(Q);
        const auto &l = shell.am();
        const auto *exps = shell.exps();
        const double PI_3_2 = M_PI * sqrt(M_PI);
        double N2 = sqrt((pow(2.0, l) * pow(2.0 * exps[0], l + 1.5)) / (PI_3_2 * df[2 * l]));
        double N1;
        if (shell.nprimitive() > 1) throw PSIEXCEPTION("Multipoles must be specified as uncontracted basis functions");
        if (l % 2 == 1) {
            // clang-format off
            //                  2
            //    /inf  2n+1 -ax           n!
            //    |    x    e    dx  =  -------
            //    /0                       n+1
            //                           2a
            // clang-format on
            N1 = pow(2.0, 0.5 * l + 0.5) * pow(exps[0], 0.5 * l + 2) / (PI_3_2 * df[l]);
        } else {
            // clang-format off
            //                2
            //    /inf  2n -ax        (2n-1)!!     ---
            //    |    x  e    dx  =  --------    / pi
            //    /0                    n  n+1   / ---
            //                        a   2    \/   a
            // clang-format on
            N1 = pow(2.0, 0.5 * l) * pow(exps[0], 0.5 * l + 1.5) / (PI_3_2 * df[l]);
        }
        // When a user specifies a positive coefficient, we assume that function represents a positive charge.  Because
        // the diffuse bases are modeling negative charges, we flip the sign of the coefficient to account for this.
        double normfac = -N1 / N2;
        if (shell.is_pure()) {
            for (int m = 0; m < 2 * l + 1; ++m) {
                scale_facs.push_back(normfac);
            }
        } else {
            // Cartesian functions in Psi4 are normalized according to the CCA standard, so only axis aligned functions
            // (e.g. xxx, yyy, zzz) are normalized to unity.  Here we provide the extra angular momentum normalization
            // to ensure that all functions are unit normalized.
            for (int lx = l; lx >= 0; lx--) {
                for (int lz = 0; lz <= l - lx; lz++) {
                    int ly = l - lx - lz;
                    double angmomfac = sqrt(df[lx] * df[ly] * df[lz] / df[l]);
                    scale_facs.push_back(normfac * angmomfac);
                }
            }
        }
    }
    if (coefs->dim() != scale_facs.size()) {
        throw PSIEXCEPTION(
            "The number of coefficients specified is inconsistent with the number of multipoles provided");
    }
    auto scaled = std::shared_ptr<Vector>(new Vector(coefs->dim()));
    std::transform(scale_facs.cbegin(), scale_facs.cend(), coefs->pointer(), scaled->pointer(), std::multiplies<>{});
    bases_.push_back(std::make_pair(basis, scaled));
}

void ExternalPotential::print(std::string out) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => External Potential Field: %s <= \n\n", name_.c_str());

    // Charges
    if (charges_.size()) {
        printer->Printf("    > Charges [a.u.] < \n\n");
        printer->Printf("     %10s %10s %10s %10s\n", "Z", "x", "y", "z");
        for (size_t i = 0; i < charges_.size(); i++) {
            printer->Printf("     %10.5f %10.5f %10.5f %10.5f\n", std::get<0>(charges_[i]), std::get<1>(charges_[i]),
                            std::get<2>(charges_[i]), std::get<3>(charges_[i]));
        }
        printer->Printf("\n");
    }

    // Bases
    if (bases_.size()) {
        printer->Printf("    > Diffuse Bases < \n\n");
        for (size_t i = 0; i < bases_.size(); i++) {
            printer->Printf("    Molecule %d\n\n", i + 1);
            bases_[i].first->molecule()->print();
            printer->Printf("    Basis %d\n\n", i + 1);
            bases_[i].first->print_by_level(out, print_);
            if (print_ > 2) {
                printer->Printf("    Density Coefficients %d\n\n", i + 1);
                bases_[i].second->print();
            }
        }
    }

    // Exchange Bases
    if (exchange_bases_.size()) {
        printer->Printf("    > Diffuse Exchange Bases < \n\n");
        for (size_t i = 0; i < exchange_bases_.size(); i++) {
            printer->Printf("    Molecule %d\n\n", i + 1);
            exchange_bases_[i].first->molecule()->print();
            printer->Printf("    Basis %d\n\n", i + 1);
            exchange_bases_[i].first->print_by_level(out, print_);
            if (print_ > 2) {
                printer->Printf("    Exchange Overlap Coefficients %d\n\n", i + 1);
                exchange_bases_[i].second->print();
            }
        }
    }
}

SharedMatrix ExternalPotential::computePotentialMatrix(std::shared_ptr<BasisSet> basis) {
    int n = basis->nbf();
    auto V = std::make_shared<Matrix>("External Potential", n, n);
    auto fact = std::make_shared<IntegralFactory>(basis, basis, basis, basis);

    // Thread count
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    double convfac = 1.0;
    if (basis->molecule()->units() == Molecule::Angstrom) convfac /= pc_bohr2angstroms;

    // Monopoles
    auto Zxyz = std::make_shared<Matrix>("Charges (Z,x,y,z)", charges_.size(), 4);
    double **Zxyzp = Zxyz->pointer();
    for (size_t i = 0; i < charges_.size(); ++i) {
        Zxyzp[i][0] = std::get<0>(charges_[i]);
        Zxyzp[i][1] = convfac * std::get<1>(charges_[i]);
        Zxyzp[i][2] = convfac * std::get<2>(charges_[i]);
        Zxyzp[i][3] = convfac * std::get<3>(charges_[i]);
    }

    std::vector<SharedMatrix> V_charge;
    std::vector<std::shared_ptr<PotentialInt> > pot;
    for (size_t t = 0; t < nthreads; ++t) {
        V_charge.push_back(std::make_shared<Matrix>("External Potential (Charges)", n, n));
        V_charge[t]->zero();
        pot.push_back(std::shared_ptr<PotentialInt>(static_cast<PotentialInt *>(fact->ao_potential())));
        pot[t]->set_charge_field(Zxyz);
    }

    // Monopole potential is symmetric, so generate unique pairs of shells
    std::vector<std::pair<size_t, size_t> > ij_pairs;
    for (size_t i = 0; i < basis->nshell(); ++i) {
        for (size_t j = 0; j <= i; ++j) {
            ij_pairs.push_back(std::pair<size_t, size_t>(i, j));
        }
    }

    // Calculate monopole potential
#pragma omp parallel for schedule(guided) num_threads(nthreads)
    for (size_t p = 0; p < ij_pairs.size(); ++p) {
        size_t i = ij_pairs[p].first;
        size_t j = ij_pairs[p].second;
        size_t ni = basis->shell(i).nfunction();
        size_t nj = basis->shell(j).nfunction();
        size_t index_i = basis->shell(i).function_index();
        size_t index_j = basis->shell(j).function_index();

        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        const double *buffer = pot[rank]->buffer();
        double **Vp = V_charge[rank]->pointer();
        pot[rank]->compute_shell(i, j);

        size_t index = 0;
        for (size_t ii = index_i; ii < (index_i + ni); ++ii) {
            for (size_t jj = index_j; jj < (index_j + nj); ++jj) {
                Vp[ii][jj] = Vp[jj][ii] = buffer[index++];
            }
        }
    } // p

    for (size_t t = 0; t < nthreads; ++t) {
        V->add(V_charge[t]);
        V_charge[t].reset();
        pot[t].reset();
    }

    // Diffuse Bases
    for (size_t ind = 0; ind < bases_.size(); ind++) {
        std::shared_ptr<BasisSet> aux = bases_[ind].first;
        SharedVector d = bases_[ind].second;

        // TODO thread this
        auto fact2 = std::make_shared<IntegralFactory>(aux, BasisSet::zero_ao_basis_set(), basis, basis);
        std::shared_ptr<TwoBodyAOInt> eri(fact2->eri());

        double **Vp = V->pointer();
        double *dp = d->pointer();

        for (int Q = 0; Q < aux->nshell(); Q++) {
            int numQ = aux->shell(Q).nfunction();
            int Qstart = aux->shell(Q).function_index();
            for (int M = 0; M < basis->nshell(); M++) {
                int numM = basis->shell(M).nfunction();
                int Mstart = basis->shell(M).function_index();
                for (int N = 0; N < basis->nshell(); N++) {
                    int numN = basis->shell(N).nfunction();
                    int Nstart = basis->shell(N).function_index();

                    eri->compute_shell(Q, 0, M, N);
                    const double *buffer = eri->buffers()[0];

                    for (int oq = 0, index = 0; oq < numQ; oq++) {
                        for (int om = 0; om < numM; om++) {
                            for (int on = 0; on < numN; on++, index++) {
                                Vp[om + Mstart][on + Nstart] += dp[oq + Qstart] * buffer[index];
                            }
                        }
                    }
                }
            }
        }
    }
    // Exchange overlap Bases
    for (size_t ind = 0; ind < exchange_bases_.size(); ind++) {
        std::shared_ptr<BasisSet> aux = exchange_bases_[ind].first;
        SharedVector d = exchange_bases_[ind].second;

        // TODO thread this
        auto fact3 = std::make_shared<IntegralFactory>(aux, basis, basis, BasisSet::zero_ao_basis_set());
        std::shared_ptr<ThreeCenterOverlapInt> ov3c(fact3->overlap_3c());

        double **Vp = V->pointer();
        double *dp = d->pointer();

        for (int Q = 0; Q < aux->nshell(); Q++) {
            int numQ = aux->shell(Q).nfunction();
            int Qstart = aux->shell(Q).function_index();
            for (int M = 0; M < basis->nshell(); M++) {
                int numM = basis->shell(M).nfunction();
                int Mstart = basis->shell(M).function_index();
                for (int N = 0; N < basis->nshell(); N++) {
                    int numN = basis->shell(N).nfunction();
                    int Nstart = basis->shell(N).function_index();

                    ov3c->compute_shell(Q, M, N);
                    const double *buffer = ov3c->buffer();

                    for (int oq = 0, index = 0; oq < numQ; oq++) {
                        for (int om = 0; om < numM; om++) {
                            for (int on = 0; on < numN; on++, index++) {
                                Vp[om + Mstart][on + Nstart] += dp[oq + Qstart] * buffer[index];
                            }
                        }
                    }
                }
            }
        }
    }
    return V;
}

SharedMatrix ExternalPotential::computePotentialGradients(std::shared_ptr<BasisSet> basis, std::shared_ptr<Matrix> Dt) {
    SharedMolecule mol = basis->molecule();
    int natom = mol->natom();
    int nextc = charges_.size();
    auto grad = std::make_shared<Matrix>("External Potential Gradient", natom, 3);
    double **Gp = grad->pointer();

    auto Zxyz = std::make_shared<Matrix>("Charges (Z,x,y,z)", charges_.size(), 4);
    double **Zxyzp = Zxyz->pointer();

    double convfac = 1.0;
    if (mol->units() == Molecule::Angstrom) convfac /= pc_bohr2angstroms;

    for (size_t i = 0; i < charges_.size(); i++) {
        Zxyzp[i][0] = std::get<0>(charges_[i]);
        Zxyzp[i][1] = convfac * std::get<1>(charges_[i]);
        Zxyzp[i][2] = convfac * std::get<2>(charges_[i]);
        Zxyzp[i][3] = convfac * std::get<3>(charges_[i]);
    }

    // Start with the nuclear contribution
    for (int cen = 0; cen < natom; ++cen) {
        double xc = mol->x(cen);
        double yc = mol->y(cen);
        double zc = mol->z(cen);
        double cencharge = mol->Z(cen);
        for (int ext = 0; ext < nextc; ++ext) {
            double charge = cencharge * Zxyzp[ext][0];
            double x = Zxyzp[ext][1] - xc;
            double y = Zxyzp[ext][2] - yc;
            double z = Zxyzp[ext][3] - zc;
            double r2 = x * x + y * y + z * z;
            double r = sqrt(r2);
            Gp[cen][0] += charge * x / (r * r2);
            Gp[cen][1] += charge * y / (r * r2);
            Gp[cen][2] += charge * z / (r * r2);
        }
    }

    // Thread count
    int threads = 1;
#ifdef _OPENMP
    threads = Process::environment.get_n_threads();
#endif
    // Potential derivative matrices for each thread
    std::vector<SharedMatrix> Vtemps;
    for (int t = 0; t < threads; t++) {
        Vtemps.push_back(std::make_shared<Matrix>("Gradient thread temp", natom, 3));
    }

    // Lower Triangular shell pairs
    std::vector<std::pair<int, int> > PQ_pairs;
    for (int P = 0; P < basis->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }
    long int nPQ = PQ_pairs.size();

    if (charges_.size()) {
        // Now the electronic contribution.
        auto fact = std::make_shared<IntegralFactory>(basis, basis, basis, basis);
        std::vector<std::shared_ptr<PotentialInt> > Vint;
        for (int t = 0; t < threads; t++) {
            Vint.push_back(std::shared_ptr<PotentialInt>(dynamic_cast<PotentialInt *>(fact->ao_potential(1))));
            Vint[t]->set_charge_field(Zxyz);
        }
#pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (long int PQ = 0L; PQ < nPQ; PQ++) {
            int P = PQ_pairs[PQ].first;
            int Q = PQ_pairs[PQ].second;

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            Vint[thread]->compute_shell_deriv1_no_charge_term(P, Q);
            const double *buffer = Vint[thread]->buffer();
            const auto &shellP = basis->shell(P);
            const auto &shellQ = basis->shell(Q);

            int aP = shellP.ncenter();
            int nP = shellP.nfunction();
            int oP = shellP.function_index();

            int aQ = shellQ.ncenter();
            int nQ = shellQ.nfunction();
            int oQ = shellQ.function_index();

            double perm = (P == Q ? 1.0 : 2.0);

            double **Vp = Vtemps[thread]->pointer();
            double **Dp = Dt->pointer();

            int npq = nP * nQ;
            const double *intPx = buffer + 3 * aP * npq + 0 * npq;
            const double *intPy = buffer + 3 * aP * npq + 1 * npq;
            const double *intPz = buffer + 3 * aP * npq + 2 * npq;
            const double *intQx = buffer + 3 * aQ * npq + 0 * npq;
            const double *intQy = buffer + 3 * aQ * npq + 1 * npq;
            const double *intQz = buffer + 3 * aQ * npq + 2 * npq;

            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    double Vval = perm * Dp[p + oP][q + oQ];
                    Vp[aP][0] += Vval * (*intPx);
                    Vp[aP][1] += Vval * (*intPy);
                    Vp[aP][2] += Vval * (*intPz);
                    if (aP != aQ) {
                        Vp[aQ][0] += Vval * (*intQx);
                        Vp[aQ][1] += Vval * (*intQy);
                        Vp[aQ][2] += Vval * (*intQz);
                    }
                    ++intPx;
                    ++intPy;
                    ++intPz;
                    ++intQx;
                    ++intQy;
                    ++intQz;
                }
            }
        }
    }

    if (bases_.size()) {
        //                       x              x
        // Add the contribution V  <- d_A (A|PQ)
        //

        // Potential derivatives
        for (size_t ind = 0; ind < bases_.size(); ind++) {
            std::shared_ptr<BasisSet> aux = bases_[ind].first;
            SharedVector d = bases_[ind].second;
            const double *pD = d->pointer();

            const auto &zero = BasisSet::zero_ao_basis_set();
            auto Afact = std::make_shared<IntegralFactory>(aux, zero, zero, zero);
            std::shared_ptr<PotentialInt> pot(static_cast<PotentialInt *>(Afact->ao_potential(1)));
            auto Zxyz = std::make_shared<Matrix>("Charges (Z,x,y,z)", 1, 4);
            double **Zxyzp = Zxyz->pointer();
            for (int atom = 0; atom < mol->natom(); atom++) {
                Zxyzp[0][0] = mol->Z(atom);
                Zxyzp[0][1] = mol->x(atom);
                Zxyzp[0][2] = mol->y(atom);
                Zxyzp[0][3] = mol->z(atom);
                pot->set_charge_field(Zxyz);

                for (int A = 0; A < aux->nshell(); A++) {
                    pot->compute_shell_deriv1_no_charge_term(A, 0);
                    const auto *buffer = pot->buffer();

                    const auto &shellA = aux->shell(A);
                    int aA = shellA.ncenter();
                    int nA = shellA.nfunction();
                    int oA = shellA.function_index();

                    const double *intAx = buffer + 3 * aA * nA + 0 * nA;
                    const double *intAy = buffer + 3 * aA * nA + 1 * nA;
                    const double *intAz = buffer + 3 * aA * nA + 2 * nA;

                    for (int a = 0; a < nA; a++) {
                        double prefac = pD[a + oA];
                        Gp[atom][0] -= prefac * (*intAx);
                        Gp[atom][1] -= prefac * (*intAy);
                        Gp[atom][2] -= prefac * (*intAz);
                        ++intAx;
                        ++intAy;
                        ++intAz;
                    }
                }
            }

            auto APQfact = std::make_shared<IntegralFactory>(aux, zero, basis, basis);
            std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
            for (int t = 0; t < threads; t++) {
                eri.push_back(std::shared_ptr<TwoBodyAOInt>(APQfact->eri(1)));
            }
            long int nAPQ = aux->nshell() * nPQ;
#pragma omp parallel for schedule(dynamic) num_threads(threads)
            for (long int APQ = 0L; APQ < nAPQ; APQ++) {
                int thread = 0;
#ifdef _OPENMP
                thread = omp_get_thread_num();
#endif
                size_t A = APQ / nPQ;
                size_t PQ = APQ % nPQ;
                int P = PQ_pairs[PQ].first;
                int Q = PQ_pairs[PQ].second;

                eri[thread]->compute_shell_deriv1(A, 0, P, Q);

                double **Vp = Vtemps[thread]->pointer();
                double **Dp = Dt->pointer();

                // We don't need any derivatives w.r.t. the external potential center (A| here
                const auto &shellA = aux->shell(A);
                int nA = shellA.nfunction();
                int oA = shellA.function_index();

                const auto &shellP = basis->shell(P);
                int nP = shellP.nfunction();
                int aP = shellP.ncenter();
                int oP = shellP.function_index();

                const auto &shellQ = basis->shell(Q);
                int nQ = shellQ.nfunction();
                int aQ = shellQ.ncenter();
                int oQ = shellQ.function_index();

                const auto &buffers = eri[thread]->buffers();
                const double *Px = buffers[3];
                const double *Py = buffers[4];
                const double *Pz = buffers[5];
                const double *Qx = buffers[6];
                const double *Qy = buffers[7];
                const double *Qz = buffers[8];

                double perm = (P == Q ? 1.0 : 2.0);

                for (int a = 0; a < nA; a++) {
                    double scale_fac = perm * pD[a + oA];
                    for (int p = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++) {
                            double val = scale_fac * Dp[p + oP][q + oQ];
                            Vp[aP][0] += val * (*Px);
                            Vp[aP][1] += val * (*Py);
                            Vp[aP][2] += val * (*Pz);
                            Vp[aQ][0] += val * (*Qx);
                            Vp[aQ][1] += val * (*Qy);
                            Vp[aQ][2] += val * (*Qz);

                            Px++;
                            Py++;
                            Pz++;
                            Qx++;
                            Qy++;
                            Qz++;
                        }
                    }
                }
            }
        }
    }
    for (int t = 0; t < threads; t++) {
        grad->add(Vtemps[t]);
    }
    return grad;
}

double ExternalPotential::computeNuclearEnergy(std::shared_ptr<Molecule> mol) {
    double E = 0.0;
    double convfac = 1.0;
    if (mol->units() == Molecule::Angstrom) convfac /= pc_bohr2angstroms;

    // Nucleus-charge interaction
    for (int A = 0; A < mol->natom(); A++) {
        double xA = mol->x(A);
        double yA = mol->y(A);
        double zA = mol->z(A);
        double ZA = mol->Z(A);

        if (ZA > 0) {  // skip Ghost interaction
            for (size_t B = 0; B < charges_.size(); B++) {
                double ZB = std::get<0>(charges_[B]);
                double xB = convfac * std::get<1>(charges_[B]);
                double yB = convfac * std::get<2>(charges_[B]);
                double zB = convfac * std::get<3>(charges_[B]);
                double dx = xA - xB;
                double dy = yA - yB;
                double dz = zA - zB;
                double R = sqrt(dx * dx + dy * dy + dz * dz);
                E += ZA * ZB / R;
            }
        }
    }

    if (bases_.size()) {
        auto Zxyz = std::make_shared<Matrix>("Nuclear Charges (Z,x,y,z)", mol->natom(), 4);
        double **Zxyzp = Zxyz->pointer();
        for (int A = 0; A < mol->natom(); A++) {
            Zxyzp[A][0] = mol->Z(A);
            Zxyzp[A][1] = mol->x(A);
            Zxyzp[A][2] = mol->y(A);
            Zxyzp[A][3] = mol->z(A);
        }
        for (size_t ind = 0; ind < bases_.size(); ind++) {
            // QM nucleus-diffuse interaction
            std::shared_ptr<BasisSet> aux = bases_[ind].first;
            SharedVector d = bases_[ind].second;
            auto V = std::make_shared<Matrix>("(Q|Z|0) Integrals", aux->nbf(), 1);
            std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
            auto fact = std::make_shared<IntegralFactory>(aux, zero, zero, zero);
            std::shared_ptr<PotentialInt> pot(static_cast<PotentialInt *>(fact->ao_potential()));
            pot->set_charge_field(Zxyz);
            pot->compute(V);
            E += C_DDOT(aux->nbf(), d->pointer(), 1, V->pointer()[0], 1);
        }
    }

    return E;
}

double ExternalPotential::computeExternExternInteraction(std::shared_ptr<ExternalPotential> other_extern, bool in_angstrom) {
    double E = 0.0;
    double convfac = 1.0; // assume the geometry of the point charges is in Bohr.
    if (in_angstrom) {
        convfac /= pc_bohr2angstroms; // Here we convert to Angstrom
    }

    // charge-charge interaction
    for (auto self_charge: charges_) {
        double ZA = std::get<0>(self_charge);
        double xA = convfac * std::get<1>(self_charge);
        double yA = convfac * std::get<2>(self_charge);
        double zA = convfac * std::get<3>(self_charge);

        for (auto other_charge: other_extern->charges_) {
            double ZB = std::get<0>(other_charge);
            double xB = convfac * std::get<1>(other_charge);
            double yB = convfac * std::get<2>(other_charge);
            double zB = convfac * std::get<3>(other_charge);

            double dx = xA - xB;
            double dy = yA - yB;
            double dz = zA - zB;
            double R = sqrt(dx * dx + dy * dy + dz * dz);

            E += ZA * ZB / R;
        }
    }

    return E;
}


}  // namespace psi
