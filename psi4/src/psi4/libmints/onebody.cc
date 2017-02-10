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
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"

#include <stdexcept>

namespace psi {

namespace {
static void transform1e_1(int am, SphericalTransformIter &sti, double *s, double *t, int nl)
{
    memset(t, 0, INT_NPURE(am) * nl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex() * nl;
        double *tptr = t + sti.pureindex() * nl;
        double coef = sti.coef();

//        outfile->Printf( "1e_1: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for (int l = 0; l < nl; l++) {
//            outfile->Printf( "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
            *(tptr++) += coef * *(sptr++);
        }
    }
}

static void transform1e_2(int am, SphericalTransformIter &sti, double *s, double *t, int nk, int nl)
{
    const int sl = nl;
    const int tl = INT_NPURE(am);

    memset(t, 0, nk * tl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex();
        double *tptr = t + sti.pureindex();
        double coef = sti.coef();

//        outfile->Printf( "1e_2: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for (int k = 0; k < nk; k++, sptr += sl, tptr += tl) {
//            outfile->Printf( "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
            *(tptr) += coef * *(sptr);
        }
    }
}
} // namespace anonmyous

OneBodyAOInt::OneBodyAOInt(std::vector <SphericalTransform> &spherical_transforms, std::shared_ptr <BasisSet> bs1, std::shared_ptr <BasisSet> bs2, int deriv)
        : bs1_(bs1), bs2_(bs2), spherical_transforms_(spherical_transforms), deriv_(deriv), nchunk_(1)
{
    force_cartesian_ = false;
    buffer_ = 0;
    natom_ = bs1_->molecule()->natom();

    size_t buffsize = INT_NCART(bs1->max_am()) * INT_NCART(bs2->max_am());

    tformbuf_ = new double[buffsize];
    target_ = new double[buffsize];
}

OneBodyAOInt::~OneBodyAOInt()
{
    delete[] tformbuf_;
    delete[] target_;
}

std::shared_ptr <BasisSet> OneBodyAOInt::basis()
{
    return bs1_;
}

std::shared_ptr <BasisSet> OneBodyAOInt::basis1()
{
    return bs1_;
}

std::shared_ptr <BasisSet> OneBodyAOInt::basis2()
{
    return bs2_;
}

const double *OneBodyAOInt::buffer() const
{
    return buffer_;
}

bool OneBodyAOInt::cloneable() const
{
    return false;
}

OneBodyAOInt *OneBodyAOInt::clone() const
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::clone()", __FILE__, __LINE__);
}

void OneBodyAOInt::normalize_am(const GaussianShell & /*s1*/, const GaussianShell & /*s2*/, int /*nchunk*/)
{
    // ACS removed this; the normalize function just returns 1.0
//    // Integrals are done. Normalize for angular momentum
//    int am1 = s1.am();
//    int am2 = s2.am();
//    int length = INT_NCART(am1) * INT_NCART(am2);

//    int ao12 = 0;
//    for(int ii = 0; ii <= am1; ii++) {
//        int l1 = am1 - ii;
//        for(int jj = 0; jj <= ii; jj++) {
//            int m1 = ii - jj;
//            int n1 = jj;
//            /*--- create all am components of sj ---*/
//            for(int kk = 0; kk <= am2; kk++) {
//                int l2 = am2 - kk;
//                for(int ll = 0; ll <= kk; ll++) {
//                    int m2 = kk - ll;
//                    int n2 = ll;

//                    for (int chunk=0; chunk<nchunk; ++chunk) {
//                        buffer_[ao12+(chunk*length)] *= GaussianShell::normalize(l1, m1, n1) * GaussianShell::normalize(l2, m2, n2);
//                    }
//                    ao12++;
//                }
//            }
//        }
//    }
}

void OneBodyAOInt::pure_transform(const GaussianShell &s1,
                                  const GaussianShell &s2, int chunks)
{
    for (int chunk = 0; chunk < chunks; ++chunk) {
        const int am1 = s1.am();
        const int is_pure1 = s1.is_pure() && am1 > 0;
        const int ncart1 = s1.ncartesian();
        const int nbf1 = s1.nfunction();

        const int am2 = s2.am();
        const int is_pure2 = s2.is_pure() && am2 > 0;
        const int ncart2 = s2.ncartesian();
        const int nbf2 = s2.nfunction();

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

void OneBodyAOInt::compute_shell(int sh1, int sh2)
{
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);

    // Call the child's compute_pair method, results better be in buffer_.
    compute_pair(s1, s2);

    // Normalize for angular momentum
    normalize_am(s1, s2, nchunk_);
    if (!force_cartesian_) {
        // Pure angular momentum (6d->5d, ...) transformation
        pure_transform(s1, s2, nchunk_);
        buffer_size_ = nchunk_ * s1.nfunction() * s2.nfunction();
    } else {
        buffer_size_ = nchunk_ * s1.ncartesian() * s2.ncartesian();
    }
}

void OneBodyAOInt::compute_pair_deriv1(const GaussianShell &, const GaussianShell &)
{
    throw PSIEXCEPTION("OneBodyAOInt::compute_pair_deriv1: Not implemented.");
}

void OneBodyAOInt::compute_pair_deriv2(const GaussianShell &, const GaussianShell &)
{
    throw PSIEXCEPTION("OneBodyAOInt::compute_pair_deriv1: Not implemented.");
}


void OneBodyAOInt::compute_shell_deriv1(int sh1, int sh2)
{
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);

    // Call the child's compute_pair method, results better be in buffer_.
    compute_pair_deriv1(s1, s2);

    // Normalize for angular momentum
    normalize_am(s1, s2, nchunk_);

    // Pure angular momentum (6d->5d, ...) transformation
    if (!force_cartesian_)
        pure_transform(s1, s2, nchunk_);
}

void OneBodyAOInt::compute(SharedMatrix &result)
{
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();

    int i_offset = 0;
    double *location;

    // Leave as this full double for loop. We could be computing nonsymmetric integrals
    for (int i = 0; i < ns1; ++i) {
        int ni = force_cartesian_ ? bs1_->shell(i).ncartesian() : bs1_->shell(i).nfunction();
        int j_offset = 0;
        for (int j = 0; j < ns2; ++j) {
            int nj = force_cartesian_ ? bs2_->shell(j).ncartesian() : bs2_->shell(j).nfunction();

            // Compute the shell (automatically transforms to pure am in needed)
            compute_shell(i, j);

            // For each integral that we got put in its contribution
            location = buffer_;
            for (int p = 0; p < ni; ++p) {
                for (int q = 0; q < nj; ++q) {
                    result->add(0, i_offset + p, j_offset + q, *location);
                    location++;
                }
            }

            j_offset += nj;
        }
        i_offset += ni;
    }
}

void OneBodyAOInt::compute(std::vector <SharedMatrix> &result)
{
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    int i_offset = 0;
    double *location = 0;

    // Check the length of result, must be chunk
    // There not an easy way of checking the size now.
    if (result.size() != (size_t) nchunk_) {
        outfile->Printf("result length = %ld, nchunk = %d\n", result.size(), nchunk_);
        throw SanityCheckError("OneBodyInt::compute(result): result incorrect length.", __FILE__, __LINE__);
    }

    // Check the individual matrices, we can only handle nirrep() == 1
    for (SharedMatrix a : result) {
        if (a->nirrep() != 1) {
            throw SanityCheckError("OneBodyInt::compute(result): one or more of the matrices given has symmetry.", __FILE__, __LINE__);
        }
    }

    for (int i = 0; i < ns1; ++i) {
        int ni = force_cartesian_ ? bs1_->shell(i).ncartesian() : bs1_->shell(i).nfunction();
        int j_offset = 0;
        for (int j = 0; j < ns2; ++j) {
            int nj = force_cartesian_ ? bs2_->shell(j).ncartesian() : bs2_->shell(j).nfunction();

            // Compute the shell
            compute_shell(i, j);

            // For each integral that we got put in its contribution
            location = buffer_;
            for (int r = 0; r < nchunk_; ++r) {
                for (int p = 0; p < ni; ++p) {
                    for (int q = 0; q < nj; ++q) {
                        result[r]->add(0, i_offset + p, j_offset + q, *location);
                        location++;
                    }
                }
            }
            j_offset += nj;
        }
        i_offset += ni;
    }
}

void OneBodyAOInt::compute_deriv1(std::vector <SharedMatrix> &result)
{
    if (deriv_ < 1)
        throw SanityCheckError("OneBodyInt::compute_deriv1(result): integral object not created to handle derivatives.", __FILE__, __LINE__);

    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    int i_offset = 0;
    double *location = 0;

    // Check the length of result, must be 3*natom_
    if (result.size() != (size_t) 3 * natom_)
        throw SanityCheckError("OneBodyInt::compute_deriv1(result): result must be 3 * natom in length.", __FILE__, __LINE__);

    if (result[0]->nirrep() != 1)
        throw SanityCheckError("OneBodyInt::compute_deriv1(result): results must be C1 symmetry.", __FILE__, __LINE__);

    for (int i = 0; i < ns1; ++i) {
        int ni = force_cartesian_ ? bs1_->shell(i).ncartesian() : bs1_->shell(i).nfunction();
        int center_i3 = 3 * bs1_->shell(i).ncenter();
        int j_offset = 0;
        for (int j = 0; j < ns2; ++j) {
            int nj = force_cartesian_ ? bs2_->shell(i).ncartesian() : bs2_->shell(j).nfunction();
            int center_j3 = 3 * bs2_->shell(j).ncenter();

            if (center_i3 != center_j3) {

                // Compute the shell
                compute_shell_deriv1(i, j);

                // Center i
                location = buffer_;
                for (int r = 0; r < 3; ++r) {
                    for (int p = 0; p < ni; ++p) {
                        for (int q = 0; q < nj; ++q) {
                            result[center_i3 + r]->add(0, i_offset + p, j_offset + q, *location);
                            location++;
                        }
                    }
                }

                // Center j -- only if center i != center j
                for (int r = 0; r < 3; ++r) {
                    for (int p = 0; p < ni; ++p) {
                        for (int q = 0; q < nj; ++q) {
                            result[center_j3 + r]->add(0, i_offset + p, j_offset + q, *location);
                            location++;
                        }
                    }
                }

            }

            j_offset += nj;
        }
        i_offset += ni;
    }
}

void OneBodyAOInt::compute_deriv2(std::vector <SharedMatrix> &result)
{
    if (deriv_ < 2)
        throw SanityCheckError("OneBodyInt::compute_deriv2(result): integral object not created to handle derivatives.", __FILE__, __LINE__);

    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    int i_offset = 0;

    // Check the length of result, must be 3*natom_
    if (result.size() != (size_t) 9 * natom_ * natom_)
        throw SanityCheckError("OneBodyInt::compute_deriv2(result): result must be 9 * natom^2 in length.", __FILE__, __LINE__);

    if (result[0]->nirrep() != 1)
        throw SanityCheckError("OneBodyInt::compute_deriv2(result): results must be C1 symmetry.", __FILE__, __LINE__);

    for (int i = 0; i < ns1; ++i) {
        int ni = force_cartesian_ ? bs1_->shell(i).ncartesian() : bs1_->shell(i).nfunction();
        int center_i3 = 3 * bs1_->shell(i).ncenter();
        int j_offset = 0;
        for (int j = 0; j < ns2; ++j) {
            int nj = force_cartesian_ ? bs2_->shell(i).ncartesian() : bs2_->shell(j).nfunction();
            int center_j3 = 3 * bs2_->shell(j).ncenter();

            if (center_i3 != center_j3) {

                // Compute the shell
                compute_shell_deriv2(i, j);

//                // Center i
//                location = buffer_;
//                for (int r=0; r<3; ++r) {
//                    for (int p=0; p<ni; ++p) {
//                        for (int q=0; q<nj; ++q) {
//                            result[center_i3+r]->add(0, i_offset+p, j_offset+q, *location);
//                            location++;
//                        }
//                    }
//                }

//                // Center j -- only if center i != center j
//                for (int r=0; r<3; ++r) {
//                    for (int p=0; p<ni; ++p) {
//                        for (int q=0; q<nj; ++q) {
//                            result[center_j3+r]->add(0, i_offset+p, j_offset+q, *location);
//                            location++;
//                        }
//                    }
//                }

            }

            j_offset += nj;
        }
        i_offset += ni;
    }
}

void OneBodyAOInt::compute_shell_deriv2(int i, int j)
{
    compute_pair_deriv2(bs1_->shell(i), bs2_->shell(j));
    // Normalize for angular momentum
    normalize_am(bs1_->shell(i), bs2_->shell(j), nchunk_);

    // Pure angular momentum (6d->5d, ...) transformation
    if (!force_cartesian_)
        pure_transform(bs1_->shell(i), bs2_->shell(j), nchunk_);
}

}
