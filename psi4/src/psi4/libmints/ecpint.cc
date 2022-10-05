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

    // Make LibECP Gaussian basis set objects for the bra...
    for (int shell = 0; shell < bs1->nshell(); ++shell){
        const GaussianShell &psi_shell = bs1->shell(shell);
        const double *center = psi_shell.center();
        std::array<double,3> C{center[0], center[1], center[2]};
        libecpint::GaussianShell newshell(C, psi_shell.am());
        for (int prim = 0; prim < psi_shell.nprimitive(); ++prim){
            newshell.addPrim(psi_shell.exp(prim), psi_shell.coef(prim));
        }
        libecp_shells1_.push_back(newshell);
    }
    // ... and the ket
    for (int shell = 0; shell < bs2->nshell(); ++shell){
        const GaussianShell &psi_shell = bs2->shell(shell);
        const double *center = psi_shell.center();
        std::array<double,3> C{center[0], center[1], center[2]};
        libecpint::GaussianShell newshell(C, psi_shell.am());
        for (int prim = 0; prim < psi_shell.nprimitive(); ++prim){
            newshell.addPrim(psi_shell.exp(prim), psi_shell.coef(prim));
        }
        libecp_shells2_.push_back(newshell);
    }

    double oldCx, oldCy, oldCz;
    std::pair<int,libecpint::ECP> ecp;
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
            ecp = std::make_pair(psi_ecp_shell.ncenter(), libecpint::ECP(psi_ecp_shell.center()));
        }
        if (oldCx != Cx || oldCy != Cy || oldCz != Cz) {
            // We're on a different center now; add the current ECP and make a new one
            oldCx = Cx;
            oldCy = Cy;
            oldCz = Cz;
            centers_and_libecp_ecps_.push_back(ecp);
            ecp = std::make_pair(psi_ecp_shell.ncenter(), libecpint::ECP(psi_ecp_shell.center()));
        }
        int nprim = psi_ecp_shell.nprimitive();
        for (int prim = 0; prim < nprim; ++prim) {
            ecp.second.addPrimitive(psi_ecp_shell.nval(prim), psi_ecp_shell.am(), psi_ecp_shell.exp(prim), psi_ecp_shell.coef(prim), prim==nprim-1);
        }
        if (ecp_shell == bs1->n_ecp_shell()-1){
            // Make sure the last one gets pushed back in!
            centers_and_libecp_ecps_.push_back(ecp);
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
        set_chunks(45);
        maxnao1 *= 45;
    }

    buffer_ = new double[nchunk_ * maxnao1 * maxnao2];
    buffers_.resize(nchunk_);
    buffers_[0] = buffer_;
}

ECPInt::~ECPInt() { delete[] buffer_; }

void ECPInt::compute_shell(int s1, int s2) {
    const libecpint::GaussianShell &LibECPShell1 = libecp_shells1_[s1];
    const libecpint::GaussianShell &LibECPShell2 = libecp_shells2_[s2];
    const size_t size = LibECPShell1.ncartesian() * LibECPShell2.ncartesian();
    memset(buffer_, 0, size * sizeof(double));
    for (const auto &center_and_ecp : centers_and_libecp_ecps_){
        libecpint::TwoIndex<double> results;
        engine_.compute_shell_pair(center_and_ecp.second, LibECPShell1, LibECPShell2, results);
        // Accumulate the results into buffer_
        std::transform (results.data.begin(), results.data.end(), buffer_, buffer_, std::plus<double>());
    }
    pure_transform(bs1_->l2_shell(s1), bs2_->l2_shell(s2));
    buffers_[0] = buffer_;
}

void ECPInt::compute_shell_deriv1(int s1, int s2) {
    const libecpint::GaussianShell &LibECPShell1 = libecp_shells1_[s1];
    const libecpint::GaussianShell &LibECPShell2 = libecp_shells2_[s2];
    const size_t size = LibECPShell1.ncartesian() * LibECPShell2.ncartesian();
    memset(buffer_, 0, 3 * natom_ * size * sizeof(double));
    int center1 = bs1_->shell(s1).ncenter();
    int center2 = bs2_->shell(s2).ncenter();
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
            std::transform(results[i].data.begin(), results[i].data.end(), buffer_ + offset, buffer_ + offset, std::plus<double>());
        }
    }
    pure_transform(bs1_->l2_shell(s1), bs2_->l2_shell(s2), nchunk_);
    for (int chunk = 0; chunk < nchunk_; ++chunk) {
        buffers_[chunk] = buffer_ + chunk * bs1_->shell(s1).nfunction() * bs2_->shell(s2).nfunction();
    }
}

void ECPInt::compute_shell_deriv2(int s1, int s2) {
    // The derivative ordering in the LibECPInt buffers
    //    0     1     2     3     4     5
    //  AxAx, AxAy, AxAz, AyAy, AyAz, AzAz,
    //    6     7     8     9    10    11    12    13    14
    //  AxBx, AxBy, AxBz, AyBx, AyBy, AyBz, AzBx, AzBy, AzBz,
    //   15    16    17    18    19    20    21    22    23
    //  AxCx, AxCy, AxCz, AyCx, AyCy, AyCz, AzCx, AzCy, AzCz,
    //   24    25    26    27    28    29
    //  BxBx, BxBy, BxBz, ByBy, ByBz, BzBz,
    //   30    31    32    33    34    35    36    37    38
    //  BxCx, BxCy, BxCz, ByCx, ByCy, ByCz, BzCx, BzCy, BzCz,
    //   39    40    41    42    43    44
    //  CxCx, CxCy, CxCz, CyCy, CyCz, CzCz
    const libecpint::GaussianShell &LibECPShell1 = libecp_shells1_[s1];
    const libecpint::GaussianShell &LibECPShell2 = libecp_shells2_[s2];
    const size_t size = LibECPShell1.ncartesian() * LibECPShell2.ncartesian();
    memset(buffer_, 0, 45 * size * sizeof(double));
    for (const auto &center_and_ecp : centers_and_libecp_ecps_){
        std::array<libecpint::TwoIndex<double>, 45> results;
        engine_.compute_shell_pair_second_derivative(center_and_ecp.second, LibECPShell1, LibECPShell2, results);
        // Accumulate the results into buffer_
        for (int i = 0; i < 45; ++i){
            const size_t offset = i * size;
            std::transform (results[i].data.begin(), results[i].data.end(), buffer_ + offset, buffer_ + offset, std::plus<double>());
        }
    }
    pure_transform(bs1_->l2_shell(s1), bs2_->l2_shell(s2), nchunk_);
    for (int chunk = 0; chunk < nchunk_; ++chunk) {
        buffers_[chunk] = buffer_ + chunk * bs1_->shell(s1).nfunction() * bs2_->shell(s2).nfunction();
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
