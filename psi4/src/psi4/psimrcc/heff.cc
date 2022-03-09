/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include <cstdio>
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "heff.h"

#include <algorithm>
#include <functional>
#include <utility>

namespace psi {

namespace psimrcc {

Hamiltonian::Hamiltonian(std::shared_ptr<PSIMRCCWfn> wfn) : wfn_(wfn) { startup(); }

Hamiltonian::~Hamiltonian() { cleanup(); }

void Hamiltonian::startup() {}

void Hamiltonian::cleanup() {}

void Hamiltonian::print_matrix() {
    if (ndets < 8) {
        outfile->Printf("\n\n  Hamiltonian Matrix\n");
        for (int mu = 0; mu < ndets; ++mu) {
            outfile->Printf("\n  ");
            for (int nu = 0; nu < ndets; ++nu) outfile->Printf(" %22.15f", matrix[mu][nu]);
        }
    }
}

void Hamiltonian::print() {
    print_matrix();

    std::vector<std::pair<double, int> > eigenvector_index_pair;
    for (int mu = 0; mu < ndets; ++mu) {
        eigenvector_index_pair.push_back(std::make_pair(right_eigenvector[mu] * right_eigenvector[mu], mu));
    }
    std::sort(eigenvector_index_pair.begin(), eigenvector_index_pair.end(), std::greater<std::pair<double, int> >());
    int max_size_list = std::min(10, static_cast<int>(eigenvector_index_pair.size()));
    outfile->Printf("\n\n  Most important determinants in the wave function");
    outfile->Printf("\n\n  determinant  eigenvector   eigenvector^2\n");
    for (int i = 0; i < max_size_list; ++i) {
        outfile->Printf("\n  %11d   %9.6f    %9.6f  %s", eigenvector_index_pair[i].second,
                        right_eigenvector[eigenvector_index_pair[i].second], eigenvector_index_pair[i].first,
                        wfn_->moinfo()->get_determinant_label(eigenvector_index_pair[i].second).c_str());
    }
}

void Hamiltonian::set_matrix(double** M, int ndets_) {
    ndets = ndets_;

    matrix.clear();
    for (int mu = 0; mu < ndets; ++mu) {
        std::vector<double> row(ndets, 0.0);
        matrix.push_back(row);
    }

    for (int mu = 0; mu < ndets; ++mu) {
        for (int nu = 0; nu < ndets; ++nu) {
            matrix[mu][nu] = M[mu][nu];
        }
    }
}

void Hamiltonian::set_left_eigenvector(const std::vector<double>& v) {
    ndets = v.size();
    left_eigenvector = std::vector<double>(v);
}

void Hamiltonian::set_right_eigenvector(const std::vector<double>& v) {
    ndets = v.size();
    right_eigenvector = std::vector<double>(v);
}

void Hamiltonian::set_zeroth_order_eigenvector(const std::vector<double>& v) {
    ndets = v.size();
    zeroth_order_eigenvector = std::vector<double>(v);
}

double Hamiltonian::expectation_value() {
    double value = 0.0;
    for (int mu = 0; mu < ndets; ++mu) {
        for (int nu = 0; nu < ndets; ++nu) {
            value += left_eigenvector[mu] * matrix[mu][nu] * right_eigenvector[nu];
        }
    }
    return (value);
}

double Hamiltonian::trace() {
    double value = 0.0;
    for (int mu = 0; mu < ndets; ++mu) {
        value += left_eigenvector[mu] * matrix[mu][mu] * right_eigenvector[mu];
    }
    return (value);
}

}  // namespace psimrcc
}  // namespace psi
