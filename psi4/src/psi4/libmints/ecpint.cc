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
            libecp_ecps_.push_back(ecp);
            ecp = libecpint::ECP(center);
        }
        int nprim = psi_ecp_shell.nprimitive();
        for (int prim = 0; prim < nprim; ++prim) {
            ecp.addPrimitive(psi_ecp_shell.nval(prim), psi_ecp_shell.am(), psi_ecp_shell.exp(prim), psi_ecp_shell.coef(prim), prim==nprim-1);
        }
        if (ecp_shell == bs1->n_ecp_shell()-1){
            // Make sure the last one gets pushed back in!
            libecp_ecps_.push_back(ecp);
        }
    }

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    buffer_ = new double[maxnao1 * maxnao2];
    buffers_.resize(1);
    buffers_[0] = buffer_;
}

ECPInt::~ECPInt() { delete[] buffer_; }


void ECPInt::compute_pair(const libint2::Shell &shellA, const libint2::Shell &shellB) {
    // Start by finding the LibECP shells, using the lookup table
    memset(buffer_, 0, shellA.cartesian_size() * shellB.cartesian_size() * sizeof(double));
    //int idxA = libecp_shell_lookup_[shellA.start()];
    //int idxB = libecp_shell_lookup_[shellB.start()];
    int idxA = 0;
    int idxB = 0;
    const libecpint::GaussianShell &LibECPShellA = libecp_shells_[idxA];
    const libecpint::GaussianShell &LibECPShellB = libecp_shells_[idxB];
    for (const auto &ecp : libecp_ecps_){
        libecpint::TwoIndex<double> results;
        engine_.compute_shell_pair(ecp, LibECPShellA, LibECPShellB, results);
        // Accumulate the results into buffer_
        std::transform (results.data.begin(), results.data.end(), buffer_, buffer_, std::plus<double>());
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
