/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2016-2017 Robert A. Shaw.
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

/* Implements ecpint.hpp */

#include "psi4/libmints/ecpint.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/molecule.h"

#include "psi4/libciomr/libciomr.h"

#include <libint2/shell.h>

#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <map>
#include <vector>

namespace psi {

ECPInt::ECPInt(std::vector<SphericalTransform> &st, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
               int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv), engine_(bs1->max_am(), bs1->max_ecp_am(), deriv) {
    int maxam1 = bs1->max_am();
    int maxam2 = bs2->max_am();

    // It will be easy to lift this restriction by duplicating the code below for bs1 and bs2
    if (bs1 != bs2)
        throw PSIEXCEPTION("Mixed basis sets are not supported for ECP integrals yet.");

    // Make LibECP Gaussian basis set objects
    for (int shell = 0; shell < bs1->nshell(); ++shell){
        const GaussianShell &psi_shell = bs1->shell(shell);
        const double *center = psi_shell.center();
        std::array<double,3> C{center[0], center[1], center[2]};
        libecpint::GaussianShell newshell(C, psi_shell.am());
        for (int prim = 0; prim < psi_shell.nprimitive(); ++prim){
            newshell.addPrim(psi_shell.exp(prim), psi_shell.coef(prim));
        }
        // Use the location of the shell's first bf as a unique ID
        libecp_shell_lookup_[psi_shell.start()] = libecp_shells_.size();
        libecp_shells_.push_back(newshell);
    }

    double oldCx, oldCy, oldCz;
    libecpint::ECP ecp;
    // Make LibECP ECP objects, grouping all functions on a given center into a single ECP object
    for (int ecp_shell = 0; ecp_shell < bs1->n_ecp_shell(); ++ecp_shell){
        const GaussianShell &psi_ecp_shell = bs1->ecp_shell(ecp_shell);
        const double *center = psi_ecp_shell.center();
        double Cx = center[0];
        double Cy = center[1];
        double Cz = center[2];
        if (ecp_shell == 0) {
            // Initialize the first ECP object
            oldCx = Cx;
            oldCy = Cy;
            oldCz = Cz;
            ecp = libecpint::ECP(psi_ecp_shell.center());
        }
        if (oldCx != Cx && oldCy != Cy && oldCz != Cz) {
            // We're on a different center now; add the current ECP and make a new one
            oldCx = Cx;
            oldCy = Cy;
            oldCz = Cz;
            centers_and_libecp_ecps_.push_back(std::make_pair(psi_ecp_shell.ncenter(),ecp));
            ecp = libecpint::ECP(center);
        }
        int nprim = psi_ecp_shell.nprimitive();
        for (int prim = 0; prim < nprim; ++prim) {
            ecp.addPrimitive(psi_ecp_shell.nval(prim), psi_ecp_shell.am(), psi_ecp_shell.exp(prim), psi_ecp_shell.coef(prim), prim==nprim-1);
        }
        if (ecp_shell == bs1->n_ecp_shell()-1){
            // Make sure the last one gets pushed back in!
            centers_and_libecp_ecps_.push_back(std::make_pair(psi_ecp_shell.ncenter(),ecp));
        }
    }

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (deriv == 1) {
        // We set chunk count for normalize_am and pure_transform
        // We can't use the trick of using less memory that I implemented in overlap & kinetic
        // since potential integral derivatives also have a contribution to center c...which is
        // over all atoms.
        set_chunks(3 * natom_);

        maxnao1 *= 3 * natom_;
    } else if (deriv == 2) {
        set_chunks(27 * natom_);
        maxnao1 *= 27 * natom_;
    }

    buffer_ = new double[maxnao1 * maxnao2];
    buffers_.resize(1);
    buffers_[0] = buffer_;
}

ECPInt::~ECPInt() { delete[] buffer_; }

void ECPInt::compute_pair(const libint2::Shell &shellA, const libint2::Shell &shellB) {
    // Start by finding the LibECP shells, using the lookup table
    const size_t size = s1.ncartesian() * s2.ncartesian();
    memset(buffer_, 0, size * sizeof(double));
    //int idx1 = libecp_shell_lookup_[s1.start()];
    //int idx2 = libecp_shell_lookup_[s2.start()];
    int idx1 = 0;
    int idx2 = 0;
    const libecpint::GaussianShell &LibECPShell1 = libecp_shells_[idx1];
    const libecpint::GaussianShell &LibECPShell2 = libecp_shells_[idx2];
    for (const auto &center_and_ecp : centers_and_libecp_ecps_){
        libecpint::TwoIndex<double> results;
        engine_.compute_shell_pair(center_and_ecp.second, LibECPShell1, LibECPShell2, results);
        // Accumulate the results into buffer_
        std::transform (results.data.begin(), results.data.end(), buffer_, buffer_, std::plus<double>());
    }
}

void ECPInt::compute_pair_deriv1(const GaussianShell &s1, const GaussianShell &s2) {
    const size_t size = s1.ncartesian() * s2.ncartesian();
    memset(buffer_, 0, 3 * natom_ * size * sizeof(double));
    int idx1 = libecp_shell_lookup_[s1.start()];
    int idx2 = libecp_shell_lookup_[s2.start()];
    const libecpint::GaussianShell &LibECPShell1 = libecp_shells_[idx1];
    const libecpint::GaussianShell &LibECPShell2 = libecp_shells_[idx2];
    int center1 = s1.ncenter();
    int center2 = s2.ncenter();
    for (const auto &center_and_ecp : centers_and_libecp_ecps_){
        int center3 = center_and_ecp.first;
        std::array<libecpint::TwoIndex<double>, 9> results;
        engine_.compute_shell_pair_derivative(center_and_ecp.second, LibECPShell1, LibECPShell2, results);
        // Accumulate the results into buffer_
        const size_t offsets[9] = { center1*3*size + 0*size, center1 * 3*size + 1*size, center1 * 3*size + 2*size,
                                    center2*3*size + 0*size, center2 * 3*size + 1*size, center2 * 3*size + 2*size,
                                    center3*3*size + 0*size, center3 * 3*size + 1*size, center3 * 3*size + 2*size };
        for (int i = 0; i < 9; ++i){
            const size_t offset = offsets[i];
            std::transform (results[i].data.begin(), results[i].data.end(), buffer_ + offset, buffer_ + offset, std::plus<double>());
        }
    }
}


ECPSOInt::ECPSOInt(const std::shared_ptr<OneBodyAOInt> &aoint, const std::shared_ptr<IntegralFactory> &fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

ECPSOInt::ECPSOInt(const std::shared_ptr<OneBodyAOInt> &aoint, const IntegralFactory *fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

}  // namespace psi
