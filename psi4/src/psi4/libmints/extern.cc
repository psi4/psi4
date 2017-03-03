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
#include "psi4/libparallel/ParallelPrinter.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

ExternalPotential::ExternalPotential() :
        debug_(0), print_(1)
{
}

ExternalPotential::~ExternalPotential()
{
}

void ExternalPotential::clear()
{
    charges_.clear();
    bases_.clear();
}

void ExternalPotential::addCharge(double Z, double x, double y, double z)
{
    charges_.push_back(std::make_tuple(Z, x, y, z));
}

void ExternalPotential::addBasis(std::shared_ptr <BasisSet> basis, SharedVector coefs)
{
    bases_.push_back(std::make_pair(basis, coefs));
}

void ExternalPotential::print(std::string out) const
{
    std::shared_ptr <psi::PsiOutStream> printer = (out == "outfile" ? outfile :
                                                   std::shared_ptr<OutFile>(new OutFile(out)));
    printer->Printf("   => External Potential Field: %s <= \n\n", name_.c_str());

    // Charges
    if (charges_.size()) {
        printer->Printf("    > Charges [a.u.] < \n\n");
        printer->Printf("     %10s %10s %10s %10s\n", "Z", "x", "y", "z");
        for (size_t i = 0; i < charges_.size(); i++) {
            printer->Printf("     %10.5f %10.5f %10.5f %10.5f\n",
                            std::get<0>(charges_[i]),
                            std::get<1>(charges_[i]),
                            std::get<2>(charges_[i]),
                            std::get<3>(charges_[i]));
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
}

SharedMatrix ExternalPotential::computePotentialMatrix(std::shared_ptr <BasisSet> basis)
{
    int n = basis->nbf();
    SharedMatrix V(new Matrix("External Potential", n, n));
    std::shared_ptr <IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));

    double convfac = 1.0;
    if (basis->molecule()->units() == Molecule::Angstrom)
        convfac /= pc_bohr2angstroms;

    // Monopoles
    SharedMatrix V_charge(new Matrix("External Potential (Charges)", n, n));

    SharedMatrix Zxyz(new Matrix("Charges (Z,x,y,z)", charges_.size(), 4));
    double **Zxyzp = Zxyz->pointer();
    for (size_t i = 0; i < charges_.size(); i++) {
        Zxyzp[i][0] = std::get<0>(charges_[i]);
        Zxyzp[i][1] = convfac * std::get<1>(charges_[i]);
        Zxyzp[i][2] = convfac * std::get<2>(charges_[i]);
        Zxyzp[i][3] = convfac * std::get<3>(charges_[i]);
    }

    std::shared_ptr <PotentialInt> pot(static_cast<PotentialInt *>(fact->ao_potential()));
    pot->set_charge_field(Zxyz);
    pot->compute(V_charge);

    V->add(V_charge);
    V_charge.reset();
    pot.reset();

    // Diffuse Bases
    for (size_t ind = 0; ind < bases_.size(); ind++) {

        std::shared_ptr <BasisSet> aux = bases_[ind].first;
        SharedVector d = bases_[ind].second;

        // TODO thread this
        std::shared_ptr <IntegralFactory> fact2(new IntegralFactory(aux, BasisSet::zero_ao_basis_set(), basis, basis));
        std::shared_ptr <TwoBodyAOInt> eri(fact2->eri());

        const double *buffer = eri->buffer();

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


SharedMatrix ExternalPotential::computePotentialGradients(std::shared_ptr <BasisSet> basis, std::shared_ptr <Matrix> Dt)
{
    // This will be easy to implement, I think, but just throw for now.
    if (bases_.size())
        throw PSIEXCEPTION("Gradients with blurred external charges are not implemented yet.");

    SharedMolecule mol = basis->molecule();
    int natom = mol->natom();
    int nextc = charges_.size();
    SharedMatrix grad(new Matrix("External Potential Gradient", natom, 3));
    double **Gp = grad->pointer();

    SharedMatrix Zxyz(new Matrix("Charges (Z,x,y,z)", charges_.size(), 4));
    double **Zxyzp = Zxyz->pointer();

    double convfac = 1.0;
    if (mol->units() == Molecule::Angstrom)
        convfac /= pc_bohr2angstroms;

    for (size_t i = 0; i < charges_.size(); i++) {
        Zxyzp[i][0] = std::get<0>(charges_[i]);
        Zxyzp[i][1] = convfac * std::get<1>(charges_[i]);
        Zxyzp[i][2] = convfac * std::get<2>(charges_[i]);
        Zxyzp[i][3] = convfac * std::get<3>(charges_[i]);
    }

    // Start with the nuclear contribution
    grad->zero();
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

    // Now the electronic contribution.
    std::shared_ptr <IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
#if 0
    // Slow, but correct, memory hog version
    std::shared_ptr<PotentialInt> potential_deriv_ints(dynamic_cast<PotentialInt*>(fact->ao_potential(1)));

    int nbf = basis->nbf();
    std::vector<SharedMatrix> intmats;
    for(int i = 0; i < 3*basis->molecule()->natom(); ++i)
        intmats.push_back(SharedMatrix(new Matrix("V derivative integrals", nbf, nbf)));


    potential_deriv_ints->set_charge_field(Zxyz);
    potential_deriv_ints->compute_deriv1_no_charge_term(intmats);
    SharedMatrix nucgrad = grad->clone();
    nucgrad->set_name("Nuclear grad");
    grad->zero();
    for (int i = 0; i < basis->molecule()->natom(); ++i){
        Gp[i][0] += Dt->vector_dot(intmats[3*i+0]);
        Gp[i][1] += Dt->vector_dot(intmats[3*i+1]);
        Gp[i][2] += Dt->vector_dot(intmats[3*i+2]);
    }
    nucgrad->print();
    grad->print();
    grad->add(nucgrad);
    grad->print();
    return grad;
#else

    // Thread count
    int threads = 1;
#ifdef _OPENMP
    threads = Process::environment.get_n_threads();
#endif

    // Potential derivatives
    std::vector <std::shared_ptr<PotentialInt>> Vint;
    std::vector <SharedMatrix> Vtemps;
    for (int t = 0; t < threads; t++) {
        Vint.push_back(std::shared_ptr<PotentialInt>(dynamic_cast<PotentialInt *>(fact->ao_potential(1))));
        Vint[t]->set_charge_field(Zxyz);
        Vtemps.push_back(SharedMatrix(grad->clone()));
        Vtemps[t]->zero();
    }

    // Lower Triangle
    std::vector <std::pair<int, int>> PQ_pairs;
    for (int P = 0; P < basis->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }

#pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        Vint[thread]->compute_shell_deriv1_no_charge_term(P, Q);
        const double *buffer = Vint[thread]->buffer();

        int nP = basis->shell(P).nfunction();
        int oP = basis->shell(P).function_index();

        int nQ = basis->shell(Q).nfunction();
        int oQ = basis->shell(Q).function_index();

        double perm = (P == Q ? 1.0 : 2.0);

        double **Vp = Vtemps[thread]->pointer();
        double **Dp = Dt->pointer();

        for (int A = 0; A < basis->molecule()->natom(); A++) {
            const double *ref0 = &buffer[3 * A * nP * nQ + 0 * nP * nQ];
            const double *ref1 = &buffer[3 * A * nP * nQ + 1 * nP * nQ];
            const double *ref2 = &buffer[3 * A * nP * nQ + 2 * nP * nQ];
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    double Vval = perm * Dp[p + oP][q + oQ];
                    Vp[A][0] += Vval * (*ref0++);
                    Vp[A][1] += Vval * (*ref1++);
                    Vp[A][2] += Vval * (*ref2++);
                }
            }
        }
    }

    for (int t = 0; t < threads; t++) {
        grad->add(Vtemps[t]);
    }
    return grad;
#endif
}

double ExternalPotential::computeNuclearEnergy(std::shared_ptr <Molecule> mol)
{
    double E = 0.0;
    double convfac = 1.0;

    if (mol->units() == Molecule::Angstrom)
        convfac /= pc_bohr2angstroms;


    // Nucleus-charge interaction
    for (int A = 0; A < mol->natom(); A++) {

        double xA = mol->x(A);
        double yA = mol->y(A);
        double zA = mol->z(A);
        double ZA = mol->Z(A);

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

    if (bases_.size()) {
        // Nucleus-diffuse interaction
        SharedMatrix Zxyz(new Matrix("Charges (Z,x,y,z)", mol->natom(), 4));
        double **Zxyzp = Zxyz->pointer();
        for (int A = 0; A < mol->natom(); A++) {
            Zxyzp[A][0] = mol->Z(A);
            Zxyzp[A][1] = mol->x(A);
            Zxyzp[A][2] = mol->y(A);
            Zxyzp[A][3] = mol->z(A);
        }

        for (size_t ind = 0; ind < bases_.size(); ind++) {

            std::shared_ptr <BasisSet> aux = bases_[ind].first;
            SharedVector d = bases_[ind].second;

            std::shared_ptr <Matrix> V(new Matrix("(Q|Z|0) Integrals", aux->nbf(), 1));

            std::shared_ptr <BasisSet> zero = BasisSet::zero_ao_basis_set();
            std::shared_ptr <IntegralFactory> fact(new IntegralFactory(aux, zero, zero, zero));
            std::shared_ptr <PotentialInt> pot(static_cast<PotentialInt *>(fact->ao_potential()));
            pot->set_charge_field(Zxyz);
            pot->compute(V);

            E += C_DDOT(aux->nbf(), d->pointer(), 1, V->pointer()[0], 1);
        }
    }

    return E;
}

}
