/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "occwave.h"

#include <cmath>

namespace psi {
namespace occwave {

void OCCWave::update_mo() {
    if (do_diis_) {
        // Convert the Array1d of the orbital amplitude update to a Vector, then add it to the total orbital amplitude
        // Vector, and DIIS extrapolate.
        const auto wogA_vec = std::make_shared<Vector>(idp_dimensions_[SpinType::Alpha], *wogA);
        if (reference_ == "RESTRICTED") {
            orbitalDiis->add_entry(2, wogA_vec.get(), kappa_bar_[SpinType::Alpha].get());
            // Choosing to DIIS extrapolate ONLY when you have all error vectors seems un-wise, but this will be changing in
            // the next commit.
            if (orbitalDiis->subspace_size() >= num_vecs) {
                orbitalDiis->extrapolate(1, kappa_bar_[SpinType::Alpha].get());
            }
        } else if (reference_ == "UNRESTRICTED") {
            const auto wogB_vec = std::make_shared<Vector>(idp_dimensions_[SpinType::Beta], *wogB);
            orbitalDiis->add_entry(4, wogA_vec.get(), wogB_vec.get(), kappa_bar_[SpinType::Alpha].get(),
                                   kappa_bar_[SpinType::Beta].get());
            if (orbitalDiis->subspace_size() >= num_vecs) {
                orbitalDiis->extrapolate(2, kappa_bar_[SpinType::Alpha].get(), kappa_bar_[SpinType::Beta].get());
            }
        }
    }

    // Now that we have DIIS'd amplitudes (if needed), perform the actual updates.
    update_mo_spincase(SpinType::Alpha);
    if (reference_ == "UNRESTRICTED") update_mo_spincase(SpinType::Beta);
}

void OCCWave::update_mo_spincase(const SpinType spin) {
    const auto& kappa_bar = kappa_bar_[spin];

    // Form the linearized orbital rotation matrix, K, from the amplitudes, kappa_bar
    auto Korb = std::make_shared<Matrix>("K MO rotation", nirrep_, nmopi_, nmopi_);
    const auto idpS = idp_dimensions_[spin].sum();
    const auto& idprowS = idprow_[spin];
    const auto& idpcolS = idpcol_[spin];
    const auto& idpirrS = idpirr_[spin];
    const auto& occpiS = occpi_[spin];
    for (int x = 0; x < idpS; x++) {
        int a = idprowS[x];
        int i = idpcolS[x];
        int h = idpirrS[x];
        Korb->set(h, a + occpiS[h], i, kappa_bar->get(x));
        Korb->set(h, i, a + occpiS[h], -kappa_bar->get(x));
    }

    // Approximate orbital rotation matrix U = exp(K) by I + K + K^2/2.
    auto Uorb = std::make_shared<Matrix>("MO rotation matrix", nirrep_, nmopi_, nmopi_);
    Uorb->identity();
    Uorb->add(Korb);
    auto Ksqr = linalg::doublet(Korb, Korb);
    Ksqr->scale(0.5);
    Uorb->add(Ksqr);

    // An MO rotation had better be orthogonal, so orthogonalize U.
    if (orth_type == "MGS") {
        double rmgs1, rmgs2;

        // loop-over nirrep_
        for (int h = 0; h < nirrep_; h++) {
            // loop-1
            for (int k = 0; k < nmopi_[h]; k++) {
                rmgs1 = 0.0;

                // loop-1a
                for (int i = 0; i < nmopi_[h]; i++) {
                    rmgs1 += Uorb->get(h, i, k) * Uorb->get(h, i, k);
                }  // end 1a

                rmgs1 = std::sqrt(rmgs1);

                // loop-1b
                for (int i = 0; i < nmopi_[h]; i++) {
                    Uorb->set(h, i, k, Uorb->get(h, i, k) / rmgs1);
                }  // end 1b

                // loop-2
                for (int j = (k + 1); j < nmopi_[h]; j++) {
                    rmgs2 = 0;

                    // loop-2a
                    for (int i = 0; i < nmopi_[h]; i++) {
                        rmgs2 += Uorb->get(h, i, k) * Uorb->get(h, i, j);
                    }  // end 2a

                    // loop-2b
                    for (int i = 0; i < nmopi_[h]; i++) {
                        Uorb->set(h, i, j, Uorb->get(h, i, j) - (rmgs2 * Uorb->get(h, i, k)));
                    }  // end 2b

                }  // end 2
            }      // end 1
        }          // end loop-over nirrep_
    }              // end main if

    else if (orth_type == "GS") {
        Uorb->schmidt();
    }

    // Now, finally update the orbitals!
    auto& CS = C_[spin];
    const auto& C_refS = C_ref_[spin];
    CS->gemm(false, false, 1.0, C_refS, Uorb, 0.0);

    if (print_ > 1) {
        Uorb->print();
        CS->print();
    }
}  // end main
}
}  // End Namespaces
