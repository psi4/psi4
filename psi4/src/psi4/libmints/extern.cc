/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include "psi4/libmints/extern.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/potential.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

ExternalPotential::ExternalPotential() : debug_(0), print_(1) {}

ExternalPotential::~ExternalPotential() {}

void ExternalPotential::clear() {
    charges_.clear();
    bases_.clear();
}

void ExternalPotential::addCharge(double Z, double x, double y, double z) {
    charges_.push_back(std::make_tuple(Z, x, y, z));
}

void ExternalPotential::addBasis(std::shared_ptr<BasisSet> basis, SharedVector coefs) {
    bases_.push_back(std::make_pair(basis, coefs));
}

SharedMatrix ExternalPotential::gradient_on_charges() { return gradient_on_charges_; }

void ExternalPotential::print(const std::string& out) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => External Potential Field: %s <= \n\n", name_.c_str());

    // Charges
    if (charges_.size()) {
        printer->Printf("    > Charges [e] [a0] < \n\n");
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
            printer->Printf("    Molecule %zu\n\n", i + 1);
            bases_[i].first->molecule()->print();
            printer->Printf("    Basis %zu\n\n", i + 1);
            bases_[i].first->print_by_level(out, print_);
            if (print_ > 2) {
                printer->Printf("    Density Coefficients %zu\n\n", i + 1);
                bases_[i].second->print();
            }
        }
    }
    
    // Matrix
    if (matrix_) {
        printer->Printf("    > One-Electron Potential Matrix < \n\n");
        matrix_->print();
        printer->Printf("\n");
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

    // Monopoles
    if (charges_.size()) {
        std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
        for (size_t i=0; i< charges_.size(); ++i) {
            Zxyz.push_back({std::get<0>(charges_[i]),{{std::get<1>(charges_[i]),
                                                       std::get<2>(charges_[i]),
                                                       std::get<3>(charges_[i])}}});
        }

        std::vector<SharedMatrix> V_charge;
        std::vector<std::shared_ptr<PotentialInt> > pot;
        for (size_t t = 0; t < nthreads; ++t) {
            V_charge.push_back(std::make_shared<Matrix>("External Potential (Charges)", n, n));
            V_charge[t]->zero();
            pot.push_back(std::shared_ptr<PotentialInt>(static_cast<PotentialInt *>(fact->ao_potential().release())));
            pot[t]->set_charge_field(Zxyz);
        }

        // Monopole potential is symmetric, so generate unique pairs of shells
        const auto& ij_pairs = pot[0]->shellpairs();

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

            double **Vp = V_charge[rank]->pointer();
            pot[rank]->compute_shell(i, j);
            const auto* buffer = pot[rank]->buffers()[0];

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
            for (int M = 0; M < basis->nshell(); M++) {
                for (int N = 0; N < basis->nshell(); N++) {
                    int numQ = aux->shell(Q).nfunction();
                    int numM = basis->shell(M).nfunction();
                    int numN = basis->shell(N).nfunction();
                    int Qstart = aux->shell(Q).function_index();
                    int Mstart = basis->shell(M).function_index();
                    int Nstart = basis->shell(N).function_index();

                    eri->compute_shell(Q, 0, M, N);
                    const double *buffer = eri->buffer();

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

    // User-provided one-electron potential matrix
    if (matrix_) {
        if (matrix_->rowdim() != matrix_->coldim())
            throw PSIEXCEPTION("ExternalPotential: user-provided 1e potential matrix must be a square matrix.");
        if (V->coldim() != matrix_->coldim())
            throw PSIEXCEPTION("ExternalPotential: user-provided 1e potential matrix must be nbf x nbf.");
        V->add(matrix_);
    }

    return V;
}

SharedMatrix ExternalPotential::computePotentialGradients(std::shared_ptr<BasisSet> basis, std::shared_ptr<Matrix> Dt) {
    // This will be easy to implement, I think, but just throw for now.
    if (bases_.size()) throw PSIEXCEPTION("Gradients with blurred external charges are not implemented yet.");

    SharedMolecule mol = basis->molecule();
    int natom = mol->natom();
    int nextc = charges_.size();

    auto grad_on_atoms = std::make_shared<Matrix>("External Potential Gradient", natom, 3);
    auto grad_on_charges = std::make_shared<Matrix>("Gradient on External Charges", nextc, 3);
    grad_on_atoms->zero();
    grad_on_charges->zero();

    double **Gp = grad_on_atoms->pointer();
    double **EGp = grad_on_charges->pointer();

    if (charges_.size()) {
        std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
        for (size_t i=0; i< charges_.size(); ++i) {
            Zxyz.push_back({std::get<0>(charges_[i]),{{std::get<1>(charges_[i]),
                                                       std::get<2>(charges_[i]),
                                                       std::get<3>(charges_[i])}}});
        }

        // Start with the nuclear contribution
        for (int cen = 0; cen < natom; ++cen) {
            double xc = mol->x(cen);
            double yc = mol->y(cen);
            double zc = mol->z(cen);
            double cencharge = mol->Z(cen);
            for (int ext = 0; ext < nextc; ++ext) {
                double charge = cencharge * Zxyz[ext].first;
                double x = Zxyz[ext].second[0] - xc;
                double y = Zxyz[ext].second[1] - yc;
                double z = Zxyz[ext].second[2] - zc;
                double r2 = x * x + y * y + z * z;
                double r = sqrt(r2);
                Gp[cen][0] += charge * x / (r * r2);
                Gp[cen][1] += charge * y / (r * r2);
                Gp[cen][2] += charge * z / (r * r2);
                EGp[ext][0] -= charge * x / (r * r2);
                EGp[ext][1] -= charge * y / (r * r2);
                EGp[ext][2] -= charge * z / (r * r2);
            }
        }

        // Now the electronic contribution.
        auto fact = std::make_shared<IntegralFactory>(basis, basis, basis, basis);

        // Thread count
        int threads = 1;
#ifdef _OPENMP
        threads = Process::environment.get_n_threads();
#endif

        // Potential derivatives
        std::vector<std::shared_ptr<PotentialInt> > Vint;
        std::vector<SharedMatrix> Vtemps;
        std::vector<SharedMatrix> EVtemps;
        for (int t = 0; t < threads; t++) {
            Vint.push_back(std::shared_ptr<PotentialInt>(dynamic_cast<PotentialInt *>(fact->ao_potential(1).release())));
            Vint[t]->set_charge_field(Zxyz);
            Vtemps.push_back(SharedMatrix(grad_on_atoms->clone()));
            Vtemps[t]->zero();
            EVtemps.push_back(SharedMatrix(grad_on_charges->clone()));
            EVtemps[t]->zero();
        }

        // Lower Triangle
        const auto &PQ_pairs = Vint[0]->shellpairs();

#pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
            int P = PQ_pairs[PQ].first;
            int Q = PQ_pairs[PQ].second;

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif
            Vint[thread]->compute_shell_deriv1(P, Q);
            const auto& buffers = Vint[thread]->buffers();

            int cP = basis->shell(P).ncenter();
            int nP = basis->shell(P).nfunction();
            int oP = basis->shell(P).function_index();

            int cQ = basis->shell(Q).ncenter();
            int nQ = basis->shell(Q).nfunction();
            int oQ = basis->shell(Q).function_index();

            double perm = (P == Q ? 1.0 : 2.0);

            double **Vp = Vtemps[thread]->pointer();
            double **EVp = EVtemps[thread]->pointer();
            double **Dp = Dt->pointer();
            const double *ref0 = buffers[0];
            const double *ref1 = buffers[1];
            const double *ref2 = buffers[2];
            const double *ref3 = buffers[3];
            const double *ref4 = buffers[4];
            const double *ref5 = buffers[5];
            std::vector<const double*> refx(3*nextc);
            for (int ext = 0; ext < nextc; ext++) {
                refx[ext*3] = buffers[ext*3 + 6];
                refx[ext*3 + 1] = buffers[ext*3 + 7];
                refx[ext*3 + 2] = buffers[ext*3 + 8];
            }
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    double Vval = perm * Dp[p + oP][q + oQ];
                    Vp[cP][0] += Vval * (*ref0++);
                    Vp[cP][1] += Vval * (*ref1++);
                    Vp[cP][2] += Vval * (*ref2++);
                    Vp[cQ][0] += Vval * (*ref3++);
                    Vp[cQ][1] += Vval * (*ref4++);
                    Vp[cQ][2] += Vval * (*ref5++);
                    const double** refxp = &refx[0];
                    for (int ext = 0; ext < nextc; ext++) {
                        EVp[ext][0] += Vval * (*refxp[ext*3]++);
                        EVp[ext][1] += Vval * (*refxp[ext*3 + 1]++);
                        EVp[ext][2] += Vval * (*refxp[ext*3 + 2]++);
                    }
                }
            }
        }

        for (int t = 0; t < threads; t++) {
            grad_on_atoms->add(Vtemps[t]);
            grad_on_charges->add(EVtemps[t]);
        }
    }
    gradient_on_charges_ = grad_on_charges;
    return grad_on_atoms;
}

double ExternalPotential::computeNuclearEnergy(std::shared_ptr<Molecule> mol) {
    double E = 0.0;

    // Nucleus-charge interaction
    for (int A = 0; A < mol->natom(); A++) {
        double xA = mol->x(A);
        double yA = mol->y(A);
        double zA = mol->z(A);
        double ZA = mol->Z(A);

        if (ZA > 0) { // skip Ghost interaction
            for (size_t B = 0; B < charges_.size(); B++) {
                double ZB = std::get<0>(charges_[B]);
                double xB = std::get<1>(charges_[B]);
                double yB = std::get<2>(charges_[B]);
                double zB = std::get<3>(charges_[B]);

                double dx = xA - xB;
                double dy = yA - yB;
                double dz = zA - zB;
                double R = sqrt(dx * dx + dy * dy + dz * dz);

                E += ZA * ZB / R;
            }
        }
    }

    if (bases_.size()) {
        // Nucleus-diffuse interaction
        std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
        for (int A = 0; A < mol->natom(); A++) {
            Zxyz.push_back({mol->Z(A),{{mol->x(A),
                                        mol->y(A),
                                        mol->z(A)}}});
        }

        for (size_t ind = 0; ind < bases_.size(); ind++) {
            std::shared_ptr<BasisSet> aux = bases_[ind].first;
            SharedVector d = bases_[ind].second;

            auto V = std::make_shared<Matrix>("(Q|Z|0) Integrals", aux->nbf(), 1);

            std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
            auto fact = std::make_shared<IntegralFactory>(aux, zero, zero, zero);
            std::shared_ptr<PotentialInt> pot(static_cast<PotentialInt *>(fact->ao_potential().release()));
            pot->set_charge_field(Zxyz);
            pot->compute(V);

            E += C_DDOT(aux->nbf(), d->pointer(), 1, V->pointer()[0], 1);
        }
    }

    return E;
}

double ExternalPotential::computeExternExternInteraction(std::shared_ptr<ExternalPotential> other_extern) {
    double E = 0.0;

    // charge-charge interaction
    for (auto self_charge: charges_) {
        double ZA = std::get<0>(self_charge);
        double xA = std::get<1>(self_charge);
        double yA = std::get<2>(self_charge);
        double zA = std::get<3>(self_charge);

        for (auto other_charge: other_extern->charges_) {
            double ZB = std::get<0>(other_charge);
            double xB = std::get<1>(other_charge);
            double yB = std::get<2>(other_charge);
            double zB = std::get<3>(other_charge);

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
