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
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <libint2/engine.h>
#include <libint2/shell.h>

#include <stdexcept>

namespace psi {

namespace {
static void transform1e_1(int am, SphericalTransformIter &sti, double *s, double *t, int nl) {
    memset(t, 0, INT_NPURE(am) * nl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex() * nl;
        double *tptr = t + sti.pureindex() * nl;
        double coef = sti.coef();

        //        outfile->Printf( "1e_1: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(),
        //        sti.coef());

        for (int l = 0; l < nl; l++) {
            //            outfile->Printf( "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
            *(tptr++) += coef * *(sptr++);
        }
    }
}

static void transform1e_2(int am, SphericalTransformIter &sti, double *s, double *t, int nk, int nl) {
    const int sl = nl;
    const int tl = INT_NPURE(am);

    memset(t, 0, nk * tl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        // In the standard spherical transform, cartindex() inverses the order
        // per shell (s + (L*(L-1)/2), s + (L*(L-1)/2) -1, ...), and
        // pureindex() is just a sequence 0, 1, ..., 2*L + 1
        double *sptr = s + sti.cartindex();
        double *tptr = t + sti.pureindex();
        double coef = sti.coef();

        //        outfile->Printf( "1e_2: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(),
        //        sti.coef());

        for (int k = 0; k < nk; k++, sptr += sl, tptr += tl) {
            //            outfile->Printf( "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
            *(tptr) += coef * *(sptr);
        }
    }
}
}  // namespace

std::vector<std::pair<int, int>> build_shell_pair_list_no_spdata(std::shared_ptr<BasisSet> bs1,
                                                                 std::shared_ptr<BasisSet> bs2) {
    const int nsh1 = bs1->nshell();
    const int nsh2 = bs2->nshell();
    const bool bs1_equiv_bs2 = (bs1 == bs2);

    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(bs1_equiv_bs2 ? static_cast<size_t>(nsh1) * (nsh1 + 1) / 2 : static_cast<size_t>(nsh1) * nsh2);
    for (int s1 = 0; s1 < nsh1; ++s1) {
        const int s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
        for (int s2 = 0; s2 <= s2_max; ++s2) pairs.emplace_back(s1, s2);
    }
    return pairs;
}

OneBodyAOInt::OneBodyAOInt(std::vector<SphericalTransform> &spherical_transforms, std::shared_ptr<BasisSet> bs1,
                           std::shared_ptr<BasisSet> bs2, int deriv)
    : bs1_(bs1), bs2_(bs2), spherical_transforms_(spherical_transforms), deriv_(deriv), nchunk_(1) {
    buffer_ = nullptr;
    natom_ = bs1_->molecule()->natom();

    size_t buffsize = INT_NCART(bs1->max_am()) * INT_NCART(bs2->max_am());

    tformbuf_ = new double[buffsize];
    target_ = new double[buffsize];

    shellpairs_ = build_shell_pair_list_no_spdata(bs1, bs2);
}

void OneBodyAOInt::require_deriv_range(const char *name, int max_deriv) const {
    if (deriv_ < 0 || deriv_ > max_deriv)
        throw PSIEXCEPTION(std::string(name) + ": derivative order " + std::to_string(deriv_) +
                           " is not in [0, " + std::to_string(max_deriv) + "].");
}

void OneBodyAOInt::require_deriv(int order) const {
    if (deriv_ < order)
        throw PSIEXCEPTION("OneBodyAOInt: the integral object was created for derivatives of order " +
                           std::to_string(deriv_) + ", not " + std::to_string(order) + ".");
}

void OneBodyAOInt::set_two_center_chunks() {
    set_chunks(deriv_ == 2 ? 21 : (deriv_ == 1 ? 6 : 1));
    buffers_.assign(nchunk_, nullptr);
}

OneBodyAOInt::~OneBodyAOInt() {
    delete[] tformbuf_;
    delete[] target_;
}

std::shared_ptr<BasisSet> OneBodyAOInt::basis() { return bs1_; }

std::shared_ptr<BasisSet> OneBodyAOInt::basis1() { return bs1_; }

std::shared_ptr<BasisSet> OneBodyAOInt::basis2() { return bs2_; }

bool OneBodyAOInt::cloneable() const { return false; }

OneBodyAOInt *OneBodyAOInt::clone() const {
    throw FeatureNotImplemented("libmints", "OneBodyInt::clone()", __FILE__, __LINE__);
}

void OneBodyAOInt::pure_transform(const libint2::Shell &s1, const libint2::Shell &s2, int chunks) {
    for (int chunk = 0; chunk < chunks; ++chunk) {
        const int am1 = s1.contr[0].l;
        const int is_pure1 = s1.contr[0].pure && am1 > 0;
        const int ncart1 = s1.cartesian_size();
        const int nbf1 = s1.size();

        const int am2 = s2.contr[0].l;
        const int is_pure2 = s2.contr[0].pure && am2 > 0;
        const int ncart2 = s2.cartesian_size();
        const int nbf2 = s2.size();

        // spherical_transforms_ only reaches IntegralFactory::init_spherical_harmonics's max_am
        const int max_pure_am = static_cast<int>(spherical_transforms_.size()) - 1;
        if ((is_pure1 && am1 > max_pure_am) || (is_pure2 && am2 > max_pure_am))
            throw PSIEXCEPTION("OneBodyAOInt: pure (spherical) angular momentum " + std::to_string(std::max(am1, am2)) +
                               " exceeds the supported maximum of " + std::to_string(max_pure_am));

        int ncart12 = ncart1 * ncart2;
        int nbf12 = nbf1 * nbf2;

        // Memory pointers that aid in transform
        double *source1, *target1;
        double *source2, *target2;
        double *source = buffer_ + (chunk * ncart12);
        double *target = target_;
        double *tmpbuf = tformbuf_;

        int transform_index = 2 * is_pure1 + is_pure2;
        switch (transform_index) {
            case 0:
                break;
            case 1:
                source2 = source;
                target2 = target;
                break;
            case 2:
                source1 = source;
                target1 = target;
                break;
            case 3:
                source2 = source;
                target2 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;
        }

        if (is_pure2) {
            SphericalTransformIter stiter(spherical_transforms_[am2]);
            transform1e_2(am2, stiter, source2, target2, ncart1, ncart2);
        }
        if (is_pure1) {
            SphericalTransformIter stiter(spherical_transforms_[am1]);
            transform1e_1(am1, stiter, source1, target1, nbf2);
        }

        if (transform_index) {
            memcpy(buffer_ + (chunk * nbf12), target_, sizeof(double) * nbf12);
        }
    }
}

void OneBodyAOInt::compute_pair(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine0_->compute(s1, s2);
    // An object created for derivatives counts their chunks (6 or 21 for overlap); computing the integrals
    // themselves gives only the engine's results. Operators without derivatives keep the count their
    // constructor chose, which can be fewer than the engine computes (NablaInt keeps 3 of 6).
    if (deriv_ > 0) set_chunks(static_cast<int>(engine0_->results().size()));
    buffers_.resize(nchunk_);
    for (int chunk = 0; chunk < nchunk_; chunk++) {
        buffers_[chunk] = engine0_->results()[chunk];
    }
    buffer_ = const_cast<double *>(buffers_[0]);
}

void OneBodyAOInt::compute_pair_deriv1(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine1_->compute(s1, s2);
    set_chunks(engine1_->results().size());
    buffers_.resize(nchunk_);
    for (int chunk = 0; chunk < nchunk_; chunk++) {
        buffers_[chunk] = engine1_->results()[chunk];
    }
}

void OneBodyAOInt::compute_pair_deriv2(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine2_->compute(s1, s2);
    set_chunks(engine2_->results().size());
    buffers_.resize(nchunk_);
    for (int chunk = 0; chunk < nchunk_; chunk++) {
        buffers_[chunk] = engine2_->results()[chunk];
    }
}

void OneBodyAOInt::compute(SharedMatrix &result) {
    const auto bs1_equiv_bs2 = (bs1_ == bs2_);

    double sign = is_antisymmetric() ? -1 : 1;
    for (auto pair : shellpairs_) {
        int p1 = pair.first;
        int p2 = pair.second;

        int ni = bs1_->shell(p1).nfunction();
        int nj = bs2_->shell(p2).nfunction();
        int i_offset = bs1_->shell_to_basis_function(p1);
        int j_offset = bs2_->shell_to_basis_function(p2);

        compute_shell(p1, p2);

        const double *location = buffers_[0];
        for (int p = 0; p < ni; ++p) {
            for (int q = 0; q < nj; ++q) {
                result->add(0, i_offset + p, j_offset + q, *location);
                if (bs1_equiv_bs2 && p1 != p2) {
                    result->add(0, j_offset + q, i_offset + p, (*location) * sign);
                }
                location++;
            }
        }
    }
}

void OneBodyAOInt::compute(std::vector<SharedMatrix> &result) {
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    const auto bs1_equiv_bs2 = (bs1_ == bs2_);

    // Check the length of result, must be chunk
    // There not an easy way of checking the size now.
    if (result.size() != (size_t)nchunk_) {
        outfile->Printf("result length = %ld, nchunk = %d\n", result.size(), nchunk_);
        throw SanityCheckError("OneBodyInt::compute(result): result incorrect length.", __FILE__, __LINE__);
    }

    // Check the individual matrices, we can only handle nirrep() == 1
    for (int chunk = 0; chunk < result.size(); ++chunk) {
        const auto a = result[chunk];
        if (a->nirrep() != 1) {
            throw SanityCheckError("OneBodyInt::compute(result): one or more of the matrices given has symmetry.",
                                   __FILE__, __LINE__);
        }
    }

    double sign = is_antisymmetric() ? -1 : 1;
    for (const auto &pair : shellpairs_) {
        int p1 = pair.first;
        int p2 = pair.second;

        int ni = bs1_->shell(p1).nfunction();
        int nj = bs2_->shell(p2).nfunction();
        int i_offset = bs1_->shell_to_basis_function(p1);
        int j_offset = bs2_->shell_to_basis_function(p2);

        // Compute the shell
        compute_shell(p1, p2);

        // For each integral that we got put in its contribution
        for (int r = 0; r < nchunk_; ++r) {
            const double *location = buffers_[r];
            for (int p = 0; p < ni; ++p) {
                for (int q = 0; q < nj; ++q) {
                    result[r]->add(0, i_offset + p, j_offset + q, *location);
                    if (bs1_equiv_bs2 && p1 != p2) {
                        result[r]->add(0, j_offset + q, i_offset + p, *location * sign);
                    }
                    location++;
                }
            }
        }
    }
}

void OneBodyAOInt::compute_shell(int sh1, int sh2) {
    const auto &s1 = bs1_->l2_shell(sh1);
    const auto &s2 = bs2_->l2_shell(sh2);
    compute_pair(s1, s2);
}

void OneBodyAOInt::compute_shell_deriv1(int sh1, int sh2) {
    require_deriv(1);
    const auto &s1 = bs1_->l2_shell(sh1);
    const auto &s2 = bs2_->l2_shell(sh2);
    compute_pair_deriv1(s1, s2);
}

void OneBodyAOInt::compute_shell_deriv2(int sh1, int sh2) {
    require_deriv(2);
    const auto &s1 = bs1_->l2_shell(sh1);
    const auto &s2 = bs2_->l2_shell(sh2);
    compute_pair_deriv2(s1, s2);
}

}  // namespace psi
