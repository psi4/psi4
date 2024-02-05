/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

//
// Created by Justin Turney on 12/17/15.
//

#include <stdexcept>

#include "integrals.h"
#include <ambit/tensor.h>
//#include <tensor/core/core.h>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"

namespace ambit {

namespace helpers {

namespace psi4 {

void PSI_API integrals(psi::OneBodyAOInt &integral, ambit::Tensor *target) {
    // One-electron integrals are generally small enough to compute
    // on a single core and broadcast them out.

    if (settings::rank == 0) {
        size_t row = static_cast<size_t>(integral.basis1()->nbf());
        size_t col = static_cast<size_t>(integral.basis2()->nbf());

        auto tmp = std::make_shared<psi::Matrix>(row, col);
        integral.compute(tmp);

        Tensor local_tensor = Tensor::build(CoreTensor, "Local Data", {row, col});

        // copy data from SharedMatrix to local_tensor
        std::copy(tmp->pointer()[0], tmp->pointer()[0] + (row * col), local_tensor.data().begin());

        // Splice data into the target tensor
        (*target)() = local_tensor();
    } else {
        Dimension zero;
        IndexRange zero_range;

        for (size_t i = 0; i < target->rank(); ++i) {
            zero.push_back(0);
            zero_range.push_back({0, 0});
        }
        Tensor local_data = Tensor::build(CoreTensor, "Local Data", zero);

        (*target)(zero_range) = local_data(zero_range);
    }
}

void PSI_API integrals(psi::TwoBodyAOInt &integral, Tensor *target) {
    if (target->type() != CoreTensor)
        throw std::runtime_error(
            "integrals(TwoBodyAOInt, Tensor) is only "
            "compatible with CoreTensor type");

    Dimension max_quartet;

    const psi::BasisSet &basis1 = *integral.basis1().get();
    const psi::BasisSet &basis2 = *integral.basis2().get();
    const psi::BasisSet &basis3 = *integral.basis3().get();
    const psi::BasisSet &basis4 = *integral.basis4().get();

    int centers_dim[4] = {-1, -1, -1, -1};

    {
        int count = 0;
        if (basis1.nbf() > 1) {
            max_quartet.push_back(static_cast<size_t>(basis1.max_function_per_shell()));
            centers_dim[0] = count++;
        }
        if (basis2.nbf() > 1) {
            max_quartet.push_back(static_cast<size_t>(basis2.max_function_per_shell()));
            centers_dim[1] = count++;
        }
        if (basis3.nbf() > 1) {
            max_quartet.push_back(static_cast<size_t>(basis3.max_function_per_shell()));
            centers_dim[2] = count++;
        }
        if (basis4.nbf() > 1) {
            max_quartet.push_back(static_cast<size_t>(basis4.max_function_per_shell()));
            centers_dim[3] = count++;
        }
    }

    if (max_quartet.size() != target->rank())
        throw std::runtime_error("TwoBodyAOInt and Tensor do not have same rank.");

    // Allocate local CoreTensor that can hold a quartets worth
    // of integrals.
    Tensor local_tensor = Tensor::build(CoreTensor, "Local Data", max_quartet);

    IndexRange target_range(target->rank());
    IndexRange local_range(target->rank());

    for (int i = 0; i < target->rank(); i++) {
        target_range[i] = {0L, 0L};
        local_range[i] = {0L, 0L};
    }

    const double *buffer = integral.buffer();
    for (int P = 0; P < basis1.nshell(); P++) {
        int nP = basis1.shell(P).nfunction();
        int startP = basis1.shell(P).function_index();
        if (centers_dim[0] != -1) {
            max_quartet[centers_dim[0]] = static_cast<size_t>(nP);
            target_range[centers_dim[0]][0] = static_cast<size_t>(startP);
            target_range[centers_dim[0]][1] = static_cast<size_t>(startP + nP);
            local_range[centers_dim[0]][1] = static_cast<size_t>(nP);
        }

        for (int Q = 0; Q < basis2.nshell(); Q++) {
            int nQ = basis2.shell(Q).nfunction();
            int startQ = basis2.shell(Q).function_index();
            if (centers_dim[1] != -1) {
                max_quartet[centers_dim[1]] = static_cast<size_t>(nQ);
                target_range[centers_dim[1]][0] = static_cast<size_t>(startQ);
                target_range[centers_dim[1]][1] = static_cast<size_t>(startQ + nQ);
                local_range[centers_dim[1]][1] = static_cast<size_t>(nQ);
            }

            for (int R = 0; R < basis3.nshell(); R++) {
                int nR = basis3.shell(R).nfunction();
                int startR = basis3.shell(R).function_index();
                if (centers_dim[2] != -1) {
                    max_quartet[centers_dim[2]] = static_cast<size_t>(nR);
                    target_range[centers_dim[2]][0] = static_cast<size_t>(startR);
                    target_range[centers_dim[2]][1] = static_cast<size_t>(startR + nR);
                    local_range[centers_dim[2]][1] = static_cast<size_t>(nR);
                }

                for (int S = 0; S < basis4.nshell(); S++) {
                    int nS = basis4.shell(S).nfunction();
                    int startS = basis4.shell(S).function_index();
                    if (centers_dim[3] != -1) {
                        max_quartet[centers_dim[3]] = static_cast<size_t>(nS);
                        target_range[centers_dim[3]][0] = static_cast<size_t>(startS);
                        target_range[centers_dim[3]][1] = static_cast<size_t>(startS + nS);
                        local_range[centers_dim[3]][1] = static_cast<size_t>(nS);
                    }

                    // Have Psi4 compute the integral
                    integral.compute_shell(P, Q, R, S);

                    // Unfortunately we have to perform a memcpy :(
                    // from Psi4 integral buffer to local_tensor
                    std::copy(buffer, buffer + (nP * nQ * nR * nS), local_tensor.data().begin());

                    // "reshape" the local tensor
                    local_tensor.reshape(max_quartet);

                    // Slice the data from local_tensor into the target.
                    // For CoreTensor this isn't the optimal procedure to
                    // follow.
                    // However, this is similar to the parallel scheme that will
                    // follow.
                    (*target)(target_range) = local_tensor(local_range);
                }
            }
        }
    }
}

}  // namespace psi4

}  // namespace helpers

}  // namespace ambit
