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

/**********************************************************
* dispersion.cc: definitions for -D for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/wavefunction.h"
#include "dispersion.h"
#include "dispersion_defines.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

namespace psi {

Dispersion::Dispersion() {}

Dispersion::~Dispersion() {}

std::shared_ptr<Dispersion> Dispersion::build(const std::string &name, double s6, double alpha6, double sr6) {
    // Options &options = Process::environment.options;
    // if (options["DFT_DISPERSION_PARAMETERS"].has_changed()) {
    //     int temp = options["DFT_DISPERSION_PARAMETERS"].size();
    //     if (temp > 0) { s6 = options["DFT_DISPERSION_PARAMETERS"][0].to_double(); }
    //     if (temp > 1) { p1 = options["DFT_DISPERSION_PARAMETERS"][1].to_double(); }
    //     if (temp > 2) { p2 = options["DFT_DISPERSION_PARAMETERS"][2].to_double(); }
    //     if (temp > 3) { p3 = options["DFT_DISPERSION_PARAMETERS"][3].to_double(); }
    //     if (temp > 4) {
    //         throw PSIEXCEPTION("DFT_DISPERSION_PARAMETERS takes no more than four elements.");
    //     }
    // }

    if (to_upper_copy(name) == "D1") {
        auto disp = std::make_shared<Dispersion>();
        disp->name_ = "-D1";
        disp->description_ = "    Grimme's -D1 Dispersion Correction\n";
        disp->citation_ = "    Grimme, S. (2004), J. Comp. Chem., 25: 1463-1473\n";
        disp->bibtex_ = "Grimme:2004:1463";
        disp->s6_ = s6;
        disp->d_ = 23.0;
        disp->C6_ = C6_D1_;
        disp->RvdW_ = RvdW_D1_;
        disp->C6_type_ = C6_arit;
        disp->Damping_type_ = Damping_D1;
        return disp;
    } else if (to_upper_copy(name) == "D2") {
        auto disp = std::make_shared<Dispersion>();
        disp->name_ = "-D2";
        disp->description_ = "    Grimme's -D2 Dispersion Correction\n";
        disp->citation_ = "    Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799\n";
        disp->bibtex_ = "Grimme:2006:1787";
        disp->s6_ = s6;
        disp->d_ = alpha6;  // 20.0
        disp->sr6_ = sr6;  // 1.1
        disp->C6_ = C6_D2_;
        disp->RvdW_ = RvdW_D2_;
        disp->C6_type_ = C6_geom;
        disp->Damping_type_ = Damping_D1;
        return disp;
    } else if (to_upper_copy(name) == "CHG") {
        auto disp = std::make_shared<Dispersion>();
        disp->name_ = "-CHG";
        disp->description_ = "    Chai and Head-Gordon Dispersion Correction\n";
        disp->citation_ = "    Chai, J.-D.; Head-Gordon, M. (2010), J. Chem. Phys., 132: 6615-6620\n";
        disp->bibtex_ = "Chai:2010:6615";
        disp->s6_ = s6;
        disp->d_ = 6.0;
        disp->C6_ = C6_D2_;
        disp->RvdW_ = RvdW_D2_;
        disp->C6_type_ = C6_geom;
        disp->Damping_type_ = Damping_CHG;
        return disp;
    } else if (to_upper_copy(name) == "DAS2009") {
        auto disp = std::make_shared<Dispersion>();
        disp->name_ = "-DAS2009";
        disp->description_ = "    Podeszwa and Szalewicz Dispersion Correction\n";
        disp->citation_ =
            "    Pernal, K.; Podeszwa, R.; Patkowski, K.; Szalewicz, K. (2009), Phys. Rev. Lett., 103: 263201\n";
        disp->bibtex_ = "Pernal:2009:263201";
        disp->s6_ = s6;
        disp->C6_ = C6_Das2009_;
        disp->C8_ = C8_Das2009_;
        disp->A_ = A_Das2009_;
        disp->Beta_ = Beta_Das2009_;
        disp->C6_type_ = C6_geom;
        disp->C8_type_ = C8_geom;
        disp->Damping_type_ = Damping_TT;
        disp->Spherical_type_ = Spherical_Das;
        return disp;
    } else if (to_upper_copy(name) == "DAS2010") {
        auto disp = std::make_shared<Dispersion>();
        disp->name_ = "-DAS2010";
        disp->description_ = "    Podeszwa and Szalewicz Dispersion Correction\n";
        disp->citation_ =
            "    Podeszwa, R.; Pernal, K.; Patkowski, K.; Szalewicz, K. (2010), J. Phys. Chem. Lett., 1: 550\n";
        disp->bibtex_ = "Podeszwa:2010:550";
        disp->s6_ = s6;
        disp->C6_ = C6_Das2010_;
        disp->C8_ = C8_Das2010_;
        disp->Beta_ = Beta_Das2010_;
        disp->C6_type_ = C6_geom;
        disp->C8_type_ = C8_geom;
        disp->Damping_type_ = Damping_TT;
        disp->Spherical_type_ = Spherical_zero;
        return disp;
    } else {
        printf("can't find %s", to_upper_copy(name).c_str());
        throw PSIEXCEPTION("Dispersion: Unknown -D type specified");
    }
}

void Dispersion::print(std::string out, int level) const {
    if (level < 1) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => %s: Empirical Dispersion <=\n\n", name_.c_str());

    printer->Printf("%s", description_.c_str());
    printer->Printf("\n");

    printer->Printf("%s", citation_.c_str());
    printer->Printf("\n");

    printer->Printf("    S6  = %14.6E\n", s6_);
    if ((name_ == "-D1") || (name_ == "-D2") || (name_ == "-CHG") || (name_ == "-D2GR"))
        printer->Printf("    A6  = %14.6E\n", d_);
    printer->Printf("\n");
}

std::string Dispersion::print_energy(std::shared_ptr<Molecule> m) {
    double e = compute_energy(m);
    std::stringstream s;
    s.setf(std::ios::scientific);
    s.precision(11);

    s << "   " << name_ << " Dispersion Energy: " << e << " [Eh]" << std::endl;

    return s.str();
}

std::string Dispersion::print_gradient(std::shared_ptr<Molecule> m) {
    SharedMatrix G = compute_gradient(m);
    double *g = G->pointer()[0];
    std::stringstream s;
    s.setf(std::ios::scientific);
    s.precision(11);

    s << "   " << name_ << " Dispersion Gradient ([a.u.]): " << std::endl << std::endl;
    s << "    Atom #:       E_x                E_y                 E_z" << std::endl;
    s << "   -----------------------------------------------------------------" << std::endl;

    for (int k = 1; k <= m->natom(); k++) {
        // clang-format off
        s << "  " << std::setw(5) << k <<
          std::setw(20) << g[(k - 1) * 3 + 0] <<
          std::setw(20) << g[(k - 1) * 3 + 1] <<
          std::setw(20) << g[(k - 1) * 3 + 2] << std::endl;
        // clang-format on
    }
    return s.str();
}

std::string Dispersion::print_hessian(std::shared_ptr<Molecule> m) {
    SharedMatrix H = compute_hessian(m);
    double **h = H->pointer();

    std::stringstream s;
    s.setf(std::ios::scientific);
    s.precision(11);

    // print_mat(h, 3*m->natom(), 3*m->natom(), outfile);

    s << "   " << name_ << " Dispersion Hessian ([a.u.]): " << std::endl << std::endl;
    for (int k = 1; k <= m->natom(); k++) {
        for (int j = 1; j <= m->natom(); j++) {
            // clang-format off
            s << "    Atom Pair A = " << k << " B = " << j << ":" << std::endl << std::endl;
            s << "                   xB                 yB                  zB" << std::endl;
            s << "   -----------------------------------------------------------------" << std::endl;

            s << "  " << std::setw(5) << "xA" <<
              std::setw(20) << h[(k - 1) * 3 + 0][(j - 1) * 3 + 0] <<
              std::setw(20) << h[(k - 1) * 3 + 0][(j - 1) * 3 + 1] <<
              std::setw(20) << h[(k - 1) * 3 + 0][(j - 1) * 3 + 2] << std::endl;
            s << "  " << std::setw(5) << "yA" <<
              std::setw(20) << h[(k - 1) * 3 + 1][(j - 1) * 3 + 0] <<
              std::setw(20) << h[(k - 1) * 3 + 1][(j - 1) * 3 + 1] <<
              std::setw(20) << h[(k - 1) * 3 + 1][(j - 1) * 3 + 2] << std::endl;
            s << "  " << std::setw(5) << "zA" <<
              std::setw(20) << h[(k - 1) * 3 + 2][(j - 1) * 3 + 0] <<
              std::setw(20) << h[(k - 1) * 3 + 2][(j - 1) * 3 + 1] <<
              std::setw(20) << h[(k - 1) * 3 + 2][(j - 1) * 3 + 2] << std::endl;
            s << std::endl;
            // clang-format on
        }
    }
    return s.str();
}

double Dispersion::compute_energy(std::shared_ptr<Molecule> m) {
    double E = 0.0;

    if (Damping_type_ == Damping_TT) {
        // -DAS dispersion only involves inter-fragment terms
        if (m->nactive_fragments() == 1) {
            // Just in case, since auto fragment is not called
            outfile->Printf("\n    Only one fragment provided, no empirical dispersion will be added.\n\n");
            return 0.0;
        }

        // need a check if there is only one active fragment ...

        // list of atoms in monomer A
        std::vector<int> realsA;
        realsA.push_back(0);
        std::vector<int> ghostsA;
        ghostsA.push_back(1);
        std::shared_ptr<Molecule> monoA = m->extract_subsets(realsA, ghostsA);
        std::shared_ptr<Vector> alist = set_atom_list(monoA);
        double *alist_p = alist->pointer();

        // list of atoms in monomer B
        std::vector<int> realsB;
        realsB.push_back(1);
        std::vector<int> ghostsB;
        ghostsB.push_back(0);
        std::shared_ptr<Molecule> monoB = m->extract_subsets(realsB, ghostsB);
        std::shared_ptr<Vector> blist = set_atom_list(monoB);
        double *blist_p = blist->pointer();

        for (int i = 0; i < monoA->natom(); i++) {
            if ((int)monoA->Z(i) == 0) continue;
            for (int j = 0; j < monoB->natom(); j++) {
                if ((int)monoB->Z(j) == 0) continue;

                double C6, C8, Rm6, Rm8, f_6, f_8, g, beta;

                double dx = monoB->x(j) - monoA->x(i);
                double dy = monoB->y(j) - monoA->y(i);
                double dz = monoB->z(j) - monoA->z(i);

                double R2 = dx * dx + dy * dy + dz * dz;

                double R = sqrt(R2);
                double R6 = R2 * R2 * R2;
                double R8 = R2 * R2 * R2 * R2;
                Rm6 = 1.0 / R6;
                Rm8 = 1.0 / R8;

                // Compute geometric mean of atomic C6 coefficients
                C6 = sqrt(C6_[(int)alist_p[i]] * C6_[(int)blist_p[j]]);

                // Compute geometric mean of atomic C8 coefficients
                C8 = sqrt(C8_[(int)alist_p[i]] * C8_[(int)blist_p[j]]);

                // Tang-Toennies Damping function
                double f_6_sum = 1.0;
                double f_8_sum = 1.0;
                beta = sqrt(Beta_[(int)alist_p[i]] * Beta_[(int)blist_p[j]]);
                for (int n = 1; n <= 6; n++) {
                    f_6_sum += pow(R * beta, n) / fac[n];
                }
                for (int n = 1; n <= 8; n++) {
                    f_8_sum += pow(R * beta, n) / fac[n];
                }

                f_6 = 1.0 - exp(-R * beta) * f_6_sum;
                f_8 = 1.0 - exp(-R * beta) * f_8_sum;

                if (Spherical_type_ == Spherical_Das) {
                    g = sqrt(A_[(int)alist_p[i]] * A_[(int)blist_p[j]]) * exp(-R * beta);
                } else if (Spherical_type_ == Spherical_zero) {
                    g = 0.0;
                } else {
                    throw PSIEXCEPTION("Unrecognized Spherical Type");
                }

                E += C6 * Rm6 * f_6;
                E += C8 * Rm8 * f_8;
                E += g;
            }
        }
    } else {
        std::shared_ptr<Vector> atom_list = set_atom_list(m);
        double *atom_list_p = atom_list->pointer();
        for (int i = 0; i < m->natom(); i++) {
            for (int j = 0; j < i; j++) {
                double C6, Rm6, f;

                double dx = m->x(j) - m->x(i);
                double dy = m->y(j) - m->y(i);
                double dz = m->z(j) - m->z(i);

                double R2 = dx * dx + dy * dy + dz * dz;
                double R = sqrt(R2);
                double R6 = R2 * R2 * R2;
                Rm6 = 1.0 / R6;

                if (C6_type_ == C6_arit) {
                    C6 = 2.0 * C6_[(int)atom_list_p[i]] * C6_[(int)atom_list_p[j]] /
                         (C6_[(int)atom_list_p[i]] + C6_[(int)atom_list_p[j]]);
                } else if (C6_type_ == C6_geom) {
                    C6 = sqrt(C6_[(int)atom_list_p[i]] * C6_[(int)atom_list_p[j]]);
                } else {
                    throw PSIEXCEPTION("Unrecognized C6 Type");
                }

                if (Damping_type_ == Damping_D1) {
                    double RvdW = (RvdW_[(int) atom_list_p[i]] + RvdW_[(int) atom_list_p[j]]) / 1.1;
                    f = 1.0 / (1.0 + exp(-d_ * (R / (sr6_ * RvdW) - 1)));
                } else if (Damping_type_ == Damping_CHG) {
                    double RvdW = RvdW_[(int)atom_list_p[i]] + RvdW_[(int)atom_list_p[j]];
                    f = 1.0 / (1.0 + d_ * pow((R / RvdW), -12.0));
                } else {
                    throw PSIEXCEPTION("Unrecognized Damping Function");
                }

                E += C6 * Rm6 * f;
            }
        }
    }
    E *= -s6_;

    return E;
}

SharedMatrix Dispersion::compute_gradient(std::shared_ptr<Molecule> m) {
    auto G = std::make_shared<Matrix>("Dispersion Gradient", m->natom(), 3);
    double **Gp = G->pointer();

    if (Damping_type_ == Damping_TT) {
        throw PSIEXCEPTION("+Das Gradients not yet implemented");
    }

    for (int i = 0; i < m->natom(); i++) {
        for (int j = 0; j < i; j++) {
            double C6, Rm6, f;
            double C6_R, Rm6_R, f_R;

            double R;
            double R_xi, R_yi, R_zi;
            double R_xj, R_yj, R_zj;

            double dx = m->x(j) - m->x(i);
            double dy = m->y(j) - m->y(i);
            double dz = m->z(j) - m->z(i);

            double R2 = dx * dx + dy * dy + dz * dz;
            R = sqrt(R2);

            R_xi = -dx / R;
            R_xj = dx / R;
            R_yi = -dy / R;
            R_yj = dy / R;
            R_zi = -dz / R;
            R_zj = dz / R;

            double R6 = R2 * R2 * R2;
            Rm6 = 1.0 / R6;
            Rm6_R = -6.0 * Rm6 / R;

            double RvdW = RvdW_[(int)m->Z(i)] + RvdW_[(int)m->Z(j)];

            if (C6_type_ == C6_arit) {
                C6 = 2.0 * C6_[(int)m->Z(i)] * C6_[(int)m->Z(j)] / (C6_[(int)m->Z(i)] + C6_[(int)m->Z(j)]);
                C6_R = 0.0;
            } else if (C6_type_ == C6_geom) {
                C6 = sqrt(C6_[(int)m->Z(i)] * C6_[(int)m->Z(j)]);
                C6_R = 0.0;
            } else {
                throw PSIEXCEPTION("Unrecognized C6 Type");
            }
            if (Damping_type_ == Damping_D1) {
                f = 1.0 / (1.0 + exp(-d_ * (R / RvdW - 1.0)));
                f_R = -f * f * exp(-d_ * (R / RvdW - 1.0)) * (-d_ / RvdW);
            } else if (Damping_type_ == Damping_CHG) {
                f = 1.0 / (1.0 + d_ * pow((R / RvdW), -12.0));
                f_R = -f * f * d_ * (-12.0) * pow((R / RvdW), -13.0) * (1.0 / RvdW);
            } else if (Damping_type_ == Damping_TT) {
                throw PSIEXCEPTION("+Das Gradients not yet implemented");
            } else {
                throw PSIEXCEPTION("Unrecognized Damping Function");
            }

            double E_R = C6_R * Rm6 * f + C6 * Rm6_R * f + C6 * Rm6 * f_R;

            Gp[i][0] += E_R * R_xi;
            Gp[i][1] += E_R * R_yi;
            Gp[i][2] += E_R * R_zi;
            Gp[j][0] += E_R * R_xj;
            Gp[j][1] += E_R * R_yj;
            Gp[j][2] += E_R * R_zj;
        }
    }

    G->scale(-s6_);
    return G;
}

SharedMatrix Dispersion::compute_hessian(std::shared_ptr<Molecule> m) {
    throw PSIEXCEPTION("Dispersion: Hessians not implemented");
}

std::shared_ptr<Vector> Dispersion::set_atom_list(std::shared_ptr<Molecule> mol) {
    auto atom_list = std::make_shared<Vector>(mol->natom());
    double *atom_list_p = atom_list->pointer();

    // look for hydrogens:
    for (int a = 0; a < mol->natom(); a++) {
        atom_list_p[a] = mol->Z(a);
        if (name_ != "-DAS2010") continue;

        if ((int)atom_list_p[a] > 54) {
            throw PsiException("libdisp does not currently support atoms with Z > 54", __FILE__, __LINE__);
        }

        if ((int)atom_list_p[a] != 1) continue;

        double Ax = mol->x(a);
        double Ay = mol->y(a);
        double Az = mol->z(a);
        double minr = 9.0e99;
        ;
        int minb = a;
        for (int b = 0; b < mol->natom(); b++) {
            if (a == b) continue;
            double r = 0.0;
            double Bx = mol->x(b);
            double By = mol->y(b);
            double Bz = mol->z(b);
            r = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
            r = sqrt(r);
            if (r < minr) {
                minr = r;
                minb = b;
            }
        }
        // what is the h bonded to?
        int atom = (int)mol->Z(minb);
        if (atom == 6)
            atom_list_p[a] = 55.0;
        else if (atom == 7)
            atom_list_p[a] = 56.0;
        else if (atom == 8)
            atom_list_p[a] = 57.0;
        else if (atom == 9)
            atom_list_p[a] = 58.0;
        else if (atom == 16)
            atom_list_p[a] = 59.0;
        else if (atom == 17)
            atom_list_p[a] = 60.0;
        else
            throw PsiException("libdisp did not find an appropriate neighbor for h", __FILE__, __LINE__);
    }
    return atom_list;
}

}  // end namespace
