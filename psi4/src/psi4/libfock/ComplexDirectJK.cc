/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "psi4/libfock/ComplexJK.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"

namespace psi {

ComplexDirectJK::ComplexDirectJK(std::shared_ptr<BasisSet> primary, Options& options) :
    ComplexJK(primary), options_(options) {

}

ComplexDirectJK::~ComplexDirectJK() {}

size_t ComplexDirectJK::memory_estimate() {
    return 0; // Effectively?
}

void ComplexDirectJK::preiterations() { /*no-op*/ }

void ComplexDirectJK::compute_JK() {
    if (!do_J_ && !do_K_) {
        outfile->Printf("\n  WARNING: ComplexDirectJK tried to compute JK, but "
                        "do_J_ and do_K_ were both set to false.\n");
        return;
    }

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    for (int N = 0; N < D_.size(); N++) {
        if (D_[N]->grid_size(0) != 1)
            throw PSIEXCEPTION("Non-C1 symmetries not allowed with ComplexJK and SCF_TYPE DIRECT");

        const auto& D_ref = D_[N]->tile(0, 0);

        auto ints = std::shared_ptr<TwoBodyAOInt>(factory->eri());
        if (do_J_ && do_K_) {
            build_JK_naive(ints, D_ref, J_[N]->tile(0, 0), K_[N]->tile(0, 0));
        } else {
            // TODO: figure out later
            throw PSIEXCEPTION("Both J and K must be computed with ComplexJK and SCF_TYPE DIRECT");
        }
    }
}

void ComplexDirectJK::build_JK_naive(std::shared_ptr<TwoBodyAOInt> ints, const ComplexT& D, ComplexT& J, ComplexT& K) {
    // TODO: Something like below
    // bool build_J = (!J.empty());
    // bool build_K = (!K.empty());
    timer_on("build_JK_naive");

    J.zero();
    K.zero();
    const int nshell = primary_->nshell();
    for (int P = 0; P < nshell; P++) {
        for (int Q = 0; Q < nshell; Q++) {
            for (int R = 0; R < nshell; R++) {
                for (int S = 0; S < nshell; S++) {
                    if (ints->compute_shell(P, Q, R, S) == 0) continue;
                    const double* buf = ints->buffer();
                    const int Psize = primary_->shell(P).nfunction();
                    const int Qsize = primary_->shell(Q).nfunction();
                    const int Rsize = primary_->shell(R).nfunction();
                    const int Ssize = primary_->shell(S).nfunction();
                    const int Poff  = primary_->shell(P).function_index();
                    const int Qoff  = primary_->shell(Q).function_index();
                    const int Roff  = primary_->shell(R).function_index();
                    const int Soff  = primary_->shell(S).function_index();
                    for (int p = 0; p < Psize; p++) {
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
                                for (int s = 0; s < Ssize; s++) {
                                    const double I = *buf++;
                                    // J_pq += D_rs * (pq|rs)
                                    J(p + Poff, q + Qoff) += D(r + Roff, s + Soff) * I;
                                    // K_pr += D_qs * (pq|rs)
                                    K(p + Poff, r + Roff) += D(q + Qoff, s + Soff) * I;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("build_JK_naive");
}


void ComplexDirectJK::postiterations() { /*no-op*/ }

void ComplexDirectJK::common_init() {
    // Whatever is needed here. Below is just an example.

    if (options_.get_int("INCFOCK_FULL_FOCK_EVERY") <= 0) {
        throw PSIEXCEPTION("Invalid input for option INCFOCK_FULL_FOCK_EVERY (<= 0)");
    }
    // other options
    auto screening_type = options_.get_str("SCREENING");
}

void ComplexDirectJK::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> DirectJK: Integral-Direct J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    Screening Type:    %11s\n", screen_type.c_str());
        outfile->Printf("    Screening Cutoff:  %11.0E\n", cutoff_);
        outfile->Printf("\n");
    }
}

} // namespace psi
