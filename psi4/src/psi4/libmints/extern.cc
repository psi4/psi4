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

const std::vector<SharedMatrix> ExternalPotential::gradients_on_bases() { return gradients_on_bases_; }

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
        printer->Printf("  Hellmann-Feynman contributions to the gradient from user-provided 1e potential"
                        "\n  matrices cannot be calculated by Psi4 and must be calculated by the user.\n\n");
        matrix_->print();
    }
}

SharedMatrix ExternalPotential::computePotentialMatrix(std::shared_ptr<BasisSet> basis) {
    auto nbf = basis->nbf();
    auto V = std::make_shared<Matrix>("External Potential", nbf, nbf);
    auto fact = std::make_shared<IntegralFactory>(basis, basis, basis, basis);

    // Thread count
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // Monopoles
    if (charges_.size()) {
        std::vector<std::pair<double, std::array<double, 3>>> Qxyz;
        for (const auto& [Q, x, y, z] : charges_) {Qxyz.emplace_back(Q, std::array{x, y, z});}

        std::vector<SharedMatrix> V_charge;
        std::vector<std::shared_ptr<PotentialInt> > pot;
        for (size_t t = 0; t < nthreads; ++t) {
            V_charge.push_back(std::make_shared<Matrix>("External Potential (Charges)", nbf, nbf));
            V_charge[t]->zero();
            pot.push_back(std::shared_ptr<PotentialInt>(static_cast<PotentialInt *>(fact->ao_potential().release())));
            pot[t]->set_charge_field(Qxyz);
        }

        // Monopole potential is symmetric, so generate unique pairs of shells
        const auto& ij_pairs = pot[0]->shellpairs();

        // Calculate monopole potential
#pragma omp parallel for schedule(guided) num_threads(nthreads)
        for (size_t p = 0; p < ij_pairs.size(); ++p) {
            auto [i, j] = ij_pairs[p];
            auto ni = basis->shell(i).nfunction();
            auto nj = basis->shell(j).nfunction();
            auto index_i = basis->shell(i).function_index();
            auto index_j = basis->shell(j).function_index();

            size_t rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif

            auto **Vp = V_charge[rank]->pointer();
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
    auto zero = BasisSet::zero_ao_basis_set();
    for (const auto& [aux, d] : bases_) {

        auto fact2 = std::make_shared<IntegralFactory>(aux, zero, basis, basis);
        std::vector<SharedMatrix> V_diffuse;
        std::vector<std::shared_ptr<TwoBodyAOInt> > eri(nthreads);
        eri[0] = std::shared_ptr<TwoBodyAOInt>(fact2->eri());
        for (int t = 1; t < nthreads; t++) eri[t] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
        for (size_t t = 0; t < nthreads; t++) {
            V_diffuse.push_back(std::make_shared<Matrix>("External Potential (Diffuses)", nbf, nbf));
            V_diffuse[t]->zero();
        }
        const auto& shell_pairs = eri[0]->shell_pairs();
        auto npairs = shell_pairs.size();
        auto *dp = d->pointer();

        // Lower Triangle
        for (int P = 0; P < aux->nshell(); P++) {

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (long int MN = 0L; MN < static_cast<long int>(npairs); MN++) {
                int rank = 0;
#ifdef _OPENMP
                rank = omp_get_thread_num();
#endif
                auto pair = shell_pairs[MN];
                auto [M, N] = pair;

                auto nP = aux->shell(P).nfunction();
                auto oP = aux->shell(P).function_index();

                auto nM = basis->shell(M).nfunction();
                auto oM = basis->shell(M).function_index();

                auto nN = basis->shell(N).nfunction();
                auto oN = basis->shell(N).function_index();

                eri[rank]->compute_shell(P, 0, M, N);
                const auto* buffer = eri[rank]->buffer();
                double **Vp = V_diffuse[rank]->pointer();

                for (int p = 0, index = 0; p < nP; p++) {
                    for (int m = 0; m < nM; m++) {
                        for (int n = 0; n < nN; n++, ++index) {
                            Vp[m + oM][n + oN] += dp[p + oP] * buffer[index];
                        }
                    }
                }
            }
        }
        for (size_t t = 0; t < nthreads; t++) {
            for (size_t i = 0; i < nbf; i++) {
                for (size_t j = i+1; j < nbf; j++) {
                    auto **Vp = V_diffuse[t]->pointer();
                    Vp[i][j] = Vp[j][i];
                }
            }
            V->add(V_diffuse[t]);
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

    SharedMolecule mol = basis->molecule();
    auto natom = mol->natom();
    auto nextc = charges_.size();

    auto grad_on_atoms = std::make_shared<Matrix>("External Potential Gradient", natom, 3);
    auto grad_on_charges = std::make_shared<Matrix>("Gradient on External Charges", nextc, 3);
    grad_on_atoms->zero();
    grad_on_charges->zero();

    auto **Gp = grad_on_atoms->pointer();
    auto **EGp = grad_on_charges->pointer();
    auto **Dp = Dt->pointer();

    std::vector<std::pair<double, std::array<double, 3>>> Qxyz;
    for (const auto& [Q, x, y, z] : charges_) {Qxyz.emplace_back(Q, std::array{x, y, z});}
    std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
    for (int A = 0; A < natom; A++) {
        Zxyz.emplace_back(mol->Z(A), std::array{mol->x(A), mol->y(A), mol->z(A)});
    }

    // Thread count
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // Evaluate contributions from point charges, if they exist.
    if (nextc > 0) {
        // Start with the nuclear contribution from point charges.
        grad_on_atoms->zero();
        grad_on_charges->zero();
        for (size_t cen = 0; cen < natom; ++cen) {
            auto xc = mol->x(cen);
            auto yc = mol->y(cen);
            auto zc = mol->z(cen);
            auto cencharge = mol->Z(cen);
            for (size_t ext = 0; const auto& [Q, Qx, Qy, Qz] : charges_) {
                auto charge = cencharge * Q;
                auto x = Qx - xc;
                auto y = Qy - yc;
                auto z = Qz - zc;
                auto r2 = x * x + y * y + z * z;
                auto r = sqrt(r2);
                Gp[cen][0] += charge * x / (r * r2);
                Gp[cen][1] += charge * y / (r * r2);
                Gp[cen][2] += charge * z / (r * r2);
                EGp[ext][0] -= charge * x / (r * r2);
                EGp[ext][1] -= charge * y / (r * r2);
                EGp[ext][2] -= charge * z / (r * r2);
                ext++;
            }
        }

        // Now the electronic contribution from point charges.
        auto fact = std::make_shared<IntegralFactory>(basis, basis, basis, basis);

        // Potential derivatives
        std::vector<std::shared_ptr<PotentialInt> > Vint;
        std::vector<SharedMatrix> Vtemps;
        std::vector<SharedMatrix> EVtemps;
        for (int t = 0; t < nthreads; t++) {
            Vint.push_back(std::shared_ptr<PotentialInt>(dynamic_cast<PotentialInt *>(fact->ao_potential(1).release())));
            Vint[t]->set_charge_field(Qxyz);
            Vtemps.push_back(SharedMatrix(grad_on_atoms->clone()));
            Vtemps[t]->zero();
            EVtemps.push_back(SharedMatrix(grad_on_charges->clone()));
            EVtemps[t]->zero();
        }

        // Lower Triangle
        const auto &PQ_pairs = Vint[0]->shellpairs();

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
            auto [P, Q] = PQ_pairs[PQ];

            int rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            Vint[rank]->compute_shell_deriv1(P, Q);
            const auto& buffers = Vint[rank]->buffers();

            auto cP = basis->shell(P).ncenter();
            auto nP = basis->shell(P).nfunction();
            auto oP = basis->shell(P).function_index();

            auto cQ = basis->shell(Q).ncenter();
            auto nQ = basis->shell(Q).nfunction();
            auto oQ = basis->shell(Q).function_index();

            double perm = (P == Q ? 1.0 : 2.0);

            auto **Vp = Vtemps[rank]->pointer();
            auto **EVp = EVtemps[rank]->pointer();
            const auto *ref0 = buffers[0];
            const auto *ref1 = buffers[1];
            const auto *ref2 = buffers[2];
            const auto *ref3 = buffers[3];
            const auto *ref4 = buffers[4];
            const auto *ref5 = buffers[5];
            std::vector<const double*> refx(3*nextc);
            for (int ext = 0; ext < nextc; ext++) {
                refx[ext*3] = buffers[ext*3 + 6];
                refx[ext*3 + 1] = buffers[ext*3 + 7];
                refx[ext*3 + 2] = buffers[ext*3 + 8];
            }
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    auto Vval = perm * Dp[p + oP][q + oQ];
                    Vp[cP][0] += Vval * (*ref0++);
                    Vp[cP][1] += Vval * (*ref1++);
                    Vp[cP][2] += Vval * (*ref2++);
                    Vp[cQ][0] += Vval * (*ref3++);
                    Vp[cQ][1] += Vval * (*ref4++);
                    Vp[cQ][2] += Vval * (*ref5++);
                    const auto** refxp = &refx[0];
                    for (int ext = 0; ext < nextc; ext++) {
                        EVp[ext][0] += Vval * (*refxp[ext*3]++);
                        EVp[ext][1] += Vval * (*refxp[ext*3 + 1]++);
                        EVp[ext][2] += Vval * (*refxp[ext*3 + 2]++);
                    }
                }
            }
        }
        for (int t = 0; t < nthreads; t++) {
            grad_on_atoms->add(Vtemps[t]);
            grad_on_charges->add(EVtemps[t]);
        }
    }

    auto zero = BasisSet::zero_ao_basis_set();
    // Evaluate contributions from diffuse charges, if they exist.
    for (const auto& [aux, d] : bases_) {

        auto grad_on_diffuses = std::make_shared<Matrix>("Gradient on External Diffuse Charges", aux->nshell(), 3);
        grad_on_diffuses->zero();
        auto **EdGp = grad_on_diffuses->pointer();

        auto fact = std::make_shared<IntegralFactory>(aux, zero, zero, zero);
        auto fact2 = std::make_shared<IntegralFactory>(aux, zero, basis, basis);
        std::shared_ptr<PotentialInt> pot(dynamic_cast<PotentialInt *>(fact->ao_potential(1).release()));
        pot->set_charge_field(Zxyz);
        auto *dp = d->pointer();

        // Start with the interaction between diffuse charges and nuclei.
        for (int P = 0; P < aux->nshell(); P++) {
            pot->compute_shell_deriv1(P, 0);
            const auto& buffers = pot->buffers();

            auto cP = aux->shell(P).ncenter();
            auto nP = aux->shell(P).nfunction();
            auto oP = aux->shell(P).function_index();

            const auto *ref0 = buffers[0];
            const auto *ref1 = buffers[1];
            const auto *ref2 = buffers[2];
            std::vector<const double*> refa(3*natom);
            for (int A = 0; A < natom; A++) {
                refa[A*3] = buffers[A*3 + 6];
                refa[A*3 + 1] = buffers[A*3 + 7];
                refa[A*3 + 2] = buffers[A*3 + 8];
            }
            for (int p = 0; p < nP; p++) {
                auto coef = dp[p + oP];
                EdGp[cP][0] += coef * (*ref0++);
                EdGp[cP][1] += coef * (*ref1++);
                EdGp[cP][2] += coef * (*ref2++);
                const auto** refap = &refa[0];
                for (int A = 0; A < natom; A++) {
                    Gp[A][0] += coef * (*refap[A*3]++);
                    Gp[A][1] += coef * (*refap[A*3 + 1]++);
                    Gp[A][2] += coef * (*refap[A*3 + 2]++);
                }
            }
        }

        std::vector<std::shared_ptr<TwoBodyAOInt> > eri(nthreads);
        eri[0] = std::shared_ptr<TwoBodyAOInt>(fact2->eri(1));
        for (int t = 1; t < nthreads; t++) eri[t] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
        std::vector<SharedMatrix> Vtemps;
        std::vector<SharedMatrix> EVtemps;
        for (int t = 0; t < nthreads; t++) {
            Vtemps.push_back(SharedMatrix(grad_on_atoms->clone()));
            Vtemps[t]->zero();
            EVtemps.push_back(SharedMatrix(grad_on_diffuses->clone()));
            EVtemps[t]->zero();
        }
        const auto& shell_pairs = eri[0]->shell_pairs();
        auto npairs = shell_pairs.size();

        // Now the interaction between diffuse charges and electrons.
        for (int P = 0; P < aux->nshell(); P++) {

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (long int MN = 0L; MN < static_cast<long int>(npairs); MN++) {
                int rank = 0;
#ifdef _OPENMP
                rank = omp_get_thread_num();
#endif
                auto pair = shell_pairs[MN];
                auto [M, N] = pair;

                eri[rank]->compute_shell_deriv1(P, 0, M, N);
                const auto & buffers = eri[rank]->buffers();
                auto **Vp = Vtemps[rank]->pointer();
                auto **EdVp = EVtemps[rank]->pointer();

                auto cP = aux->shell(P).ncenter();
                auto nP = aux->shell(P).nfunction();
                auto oP = aux->shell(P).function_index();

                auto cM = basis->shell(M).ncenter();
                auto nM = basis->shell(M).nfunction();
                auto oM = basis->shell(M).function_index();

                auto cN = basis->shell(N).ncenter();
                auto nN = basis->shell(N).nfunction();
                auto oN = basis->shell(N).function_index();

                double perm = (M == N ? 1.0 : 2.0);

                const auto *ref0 = buffers[0];
                const auto *ref1 = buffers[1];
                const auto *ref2 = buffers[2];
                const auto *ref3 = buffers[3];
                const auto *ref4 = buffers[4];
                const auto *ref5 = buffers[5];
                const auto *ref6 = buffers[6];
                const auto *ref7 = buffers[7];
                const auto *ref8 = buffers[8];

                for (int p = 0; p < nP; p++) {
                    for (int m = 0; m < nM; m++) {
                        for (int n = 0; n < nN; n++) {
                            auto coef = perm * dp[p + oP] * Dp[m + oM][n + oN];
                            EdVp[cP][0] += coef * (*ref0++);
                            EdVp[cP][1] += coef * (*ref1++);
                            EdVp[cP][2] += coef * (*ref2++);
                            Vp[cM][0] += coef * (*ref3++);
                            Vp[cM][1] += coef * (*ref4++);
                            Vp[cM][2] += coef * (*ref5++);
                            Vp[cN][0] += coef * (*ref6++);
                            Vp[cN][1] += coef * (*ref7++);
                            Vp[cN][2] += coef * (*ref8++);
                        }
                    }
                }
            }
        }
        for (int t = 0; t < nthreads; t++) {
            grad_on_atoms->add(Vtemps[t]);
            grad_on_diffuses->add(EVtemps[t]);
        }
        gradients_on_bases_.push_back(grad_on_diffuses);
    }

    gradient_on_charges_ = grad_on_charges;
    return grad_on_atoms;
}

double ExternalPotential::computeNuclearEnergy(std::shared_ptr<Molecule> mol) {
    double E = 0.0;
    auto natom = mol->natom();

    // Nucleus-charge interaction
    for (int A = 0; A < natom; A++) {
        auto xA = mol->x(A);
        auto yA = mol->y(A);
        auto zA = mol->z(A);
        auto ZA = mol->Z(A);

        if (ZA > 0) { // skip Ghost interaction
            for (size_t B = 0; B < charges_.size(); B++) {
                auto ZB = std::get<0>(charges_[B]);
                auto xB = std::get<1>(charges_[B]);
                auto yB = std::get<2>(charges_[B]);
                auto zB = std::get<3>(charges_[B]);

                auto dx = xA - xB;
                auto dy = yA - yB;
                auto dz = zA - zB;
                auto R = sqrt(dx * dx + dy * dy + dz * dz);

                E += ZA * ZB / R;
            }
        }
    }

    // Nucleus-diffuse interaction
    std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
    for (int A = 0; A < natom; A++) {
        Zxyz.emplace_back(mol->Z(A), std::array{mol->x(A), mol->y(A), mol->z(A)});
    }

    auto zero = BasisSet::zero_ao_basis_set();
    for (const auto& [aux, d] : bases_) {

        auto V = std::make_shared<Matrix>("(Q|Z|0) Integrals", aux->nbf(), 1);

        auto fact = std::make_shared<IntegralFactory>(aux, zero, zero, zero);
        std::shared_ptr<PotentialInt> pot(static_cast<PotentialInt *>(fact->ao_potential().release()));
        pot->set_charge_field(Zxyz);
        pot->compute(V);

        E += C_DDOT(aux->nbf(), d->pointer(), 1, V->pointer()[0], 1);
    }

    return E;
}

double ExternalPotential::computeExternExternInteraction(std::shared_ptr<ExternalPotential> other_extern) {
    double E = 0.0;

    // charge-charge interaction
    for (auto self_charge: charges_) {
        auto ZA = std::get<0>(self_charge);
        auto xA = std::get<1>(self_charge);
        auto yA = std::get<2>(self_charge);
        auto zA = std::get<3>(self_charge);

        for (auto other_charge: other_extern->charges_) {
            auto ZB = std::get<0>(other_charge);
            auto xB = std::get<1>(other_charge);
            auto yB = std::get<2>(other_charge);
            auto zB = std::get<3>(other_charge);

            auto dx = xA - xB;
            auto dy = yA - yB;
            auto dz = zA - zB;
            auto R = sqrt(dx * dx + dy * dy + dz * dz);

            E += ZA * ZB / R;
        }
    }

    auto zero = BasisSet::zero_ao_basis_set();
    // diffuse-charge interaction
    if (other_extern->charges_.size()) {
        std::vector<std::pair<double, std::array<double, 3>>> Qxyz;
        for (const auto& [Q, x, y, z] : other_extern->charges_) {Qxyz.emplace_back(Q, std::array{x, y, z});}
        for (const auto& [aux, d] : bases_) {
       
            auto V = std::make_shared<Matrix>("(Q|Z|0) Integrals", aux->nbf(), 1);
       
            auto fact = std::make_shared<IntegralFactory>(aux, zero, zero, zero);
            std::shared_ptr<PotentialInt> pot(static_cast<PotentialInt *>(fact->ao_potential().release()));
            pot->set_charge_field(Qxyz);
            pot->compute(V);
       
            E += C_DDOT(aux->nbf(), d->pointer(), 1, V->pointer()[0], 1);
        }
    }

    // charge-diffuse interaction
    if (charges_.size()) {
        std::vector<std::pair<double, std::array<double, 3>>> Qxyz;
        for (const auto& [Q, x, y, z] : charges_) {Qxyz.emplace_back(Q, std::array{x, y, z});}
        for (const auto& [aux, d] : other_extern->bases_) {
       
            auto V = std::make_shared<Matrix>("(Q|Z|0) Integrals", aux->nbf(), 1);
       
            auto fact = std::make_shared<IntegralFactory>(aux, zero, zero, zero);
            std::shared_ptr<PotentialInt> pot(static_cast<PotentialInt *>(fact->ao_potential().release()));
            pot->set_charge_field(Qxyz);
            pot->compute(V);
       
            E += C_DDOT(aux->nbf(), d->pointer(), 1, V->pointer()[0], 1);
        }
    }

    // diffuse-diffuse interaction
    for (const auto& [aux0, d0] : bases_) {
        auto *dp0 = d0->pointer();

        for (const auto& [aux1, d1] : other_extern->bases_) {
            auto *dp1 = d1->pointer();
        
            auto fact2 = std::make_shared<IntegralFactory>(aux0, zero, aux1, zero);
            std::shared_ptr<TwoBodyAOInt> eri(fact2->eri());

            for (int P = 0; P < aux0->nshell(); P++) {
                for (int Q = 0; Q < aux1->nshell(); Q++) {
                    eri->compute_shell(P, 0, Q, 0);
                    const auto *buffer = eri->buffer();

                    auto nP = aux0->shell(P).nfunction();
                    auto oP = aux0->shell(P).function_index();

                    auto nQ = aux1->shell(Q).nfunction();
                    auto oQ = aux1->shell(Q).function_index();

                    for (int p = 0, index = 0; p < nP; p++) {
                        for (int q = 0; q < nQ; q++, ++index) {
                                E += dp0[p + oP] * dp1[q + oQ] * buffer[index];
                        }
                    }
                }
            }
        }
    }

    return E;
}

}  // namespace psi
