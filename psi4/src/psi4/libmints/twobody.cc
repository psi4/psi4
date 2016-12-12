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

#include <stdexcept>
#include "psi4/libqt/qt.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"

;
using namespace psi;

static void transform2e_1(int, SphericalTransformIter&, double*, double*, int);
static void transform2e_2(int, SphericalTransformIter&, double*, double*, int, int, int);
static void transform2e_3(int, SphericalTransformIter&, double*, double*, int, int, int);
static void transform2e_4(int, SphericalTransformIter&, double*, double*, int, int);

TwoBodyAOInt::TwoBodyAOInt(const IntegralFactory* intsfactory, int deriv) :
    integral_(intsfactory),
    original_bs1_(integral_->basis1()),
    original_bs2_(integral_->basis2()),
    original_bs3_(integral_->basis3()),
    original_bs4_(integral_->basis4()),
    bs1_(original_bs1_),
    bs2_(original_bs2_),
    bs3_(original_bs3_),
    bs4_(original_bs4_),
    target_full_(nullptr), target_(nullptr),
    source_full_(nullptr), source_(nullptr),
    deriv_(deriv)
{
    //outfile->Printf( "TwoBodyAOInt object created with: %s, %s, %s, %s\n",
    //        original_bs1_->name().c_str(),
    //        original_bs2_->name().c_str(),
    //        original_bs3_->name().c_str(),
    //        original_bs4_->name().c_str());
    // The derived classes allocate this memory.
    force_cartesian_ = false;
    tformbuf_ = nullptr;
    natom_ = original_bs1_->molecule()->natom();  // This assumes the 4 bases come from the same molecule.
}

TwoBodyAOInt::TwoBodyAOInt(const TwoBodyAOInt & rhs)
    : TwoBodyAOInt(rhs.integral_, rhs.deriv_)
{
    blocks12_ = rhs.blocks12_;
    blocks34_ = rhs.blocks34_;
}

TwoBodyAOInt::~TwoBodyAOInt()
{
}

std::shared_ptr<BasisSet> TwoBodyAOInt::basis()
{
    return original_bs1_;
}

std::shared_ptr<BasisSet> TwoBodyAOInt::basis1()
{
    return original_bs1_;
}

std::shared_ptr<BasisSet> TwoBodyAOInt::basis2()
{
    return original_bs2_;
}

std::shared_ptr<BasisSet> TwoBodyAOInt::basis3()
{
    return original_bs3_;
}

std::shared_ptr<BasisSet> TwoBodyAOInt::basis4()
{
    return original_bs4_;
}

bool TwoBodyAOInt::cloneable() const
{
    return false;
}

TwoBodyAOInt* TwoBodyAOInt::clone() const
{
    throw FeatureNotImplemented("libmints", "TwoBodyInt::clone()", __FILE__, __LINE__);
}

std::vector<ShellPairBlock> TwoBodyAOInt::get_blocks12(void) const
{
    return blocks12_;
}

std::vector<ShellPairBlock> TwoBodyAOInt::get_blocks34(void) const
{
    return blocks34_;
}


// Expected to be overridden by derived classes
void TwoBodyAOInt::create_blocks(void)
{
    // Default implementation : do no blocking.
    // Each ShellPairBlock will only contain one shell pair
    blocks12_.clear();
    blocks34_.clear();

    const auto nshell1 = basis1()->nshell();
    const auto nshell2 = basis2()->nshell();
    const auto nshell3 = basis3()->nshell();
    const auto nshell4 = basis4()->nshell();

    blocks12_.reserve(nshell1*nshell2);
    blocks34_.reserve(nshell3*nshell4);

    for(int ishell = 0; ishell < basis1()->nshell(); ishell++)
    for(int jshell = 0; jshell < basis2()->nshell(); jshell++)
        blocks12_.push_back({{ishell, jshell}});

    for(int kshell = 0; kshell < basis3()->nshell(); kshell++)
    for(int lshell = 0; lshell < basis4()->nshell(); lshell++)
        blocks34_.push_back({{kshell, lshell}});
}


void TwoBodyAOInt::compute_shell_blocks(int shellpair12, int shellpair34,
                                        int npair12, int npair34)
{
    // Default implementation - go through the blocks and do each quartet
    // one at a time

    // reset the target & source pointers
    target_ = target_full_;
    source_ = source_full_;

    auto vsh12 = blocks12_[shellpair12];
    auto vsh34 = blocks34_[shellpair34];

    for(const auto sh12 : vsh12)
    {
        const auto & shell1 = original_bs1_->shell(sh12.first);
        const auto & shell2 = original_bs2_->shell(sh12.second);
   
        int n1, n2;
        if (force_cartesian_)
        {
            n1 = shell1.ncartesian();
            n2 = shell2.ncartesian();
        }
        else
        {
            n1 = shell1.nfunction();
            n2 = shell2.nfunction();
        }


        for(const auto sh34 : vsh34)
        {
            const auto & shell3 = original_bs3_->shell(sh34.first);
            const auto & shell4 = original_bs4_->shell(sh34.second);

            int n3, n4;
            if (force_cartesian_)
            {
                n3 = shell3.ncartesian();
                n4 = shell4.ncartesian();
            }
            else
            {
                n3 = shell3.nfunction();
                n4 = shell4.nfunction();
            }

            const int n1234 = n1 * n2 * n3 * n4;

            // actually compute the eri
            // this will put the results in target_
            auto ret = compute_shell(sh12.first, sh12.second,
                                     sh34.first, sh34.second);

            // advance the target pointer
            target_ += n1234; 
            
            // Since we are only doing one at a time we don't need to
            // move the source_ pointer
        }
    }
}

void TwoBodyAOInt::normalize_am(std::shared_ptr<GaussianShell> s1, std::shared_ptr<GaussianShell> s2, std::shared_ptr<GaussianShell> s3, std::shared_ptr<GaussianShell> s4, int nchunk)
{
    // Integrals assume this normalization is 1.0.
    return;
#if 0
#ifdef MINTS_TIMER
    timer_on("Angular momentum normalization");
#endif
    int am1 = s1->am();
    int am2 = s2->am();
    int am3 = s3->am();
    int am4 = s4->am();
    int length = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

    // Need to go through and grab all the integrals for this given shell and add them
    // to the running totals.
    int nprim = 0;
    for (int i = 0; i <= am1; i++) {
        int l1 = am1 - i;
        for (int j = 0; j <= i; j++) {
            int m1 = i - j;
            int n1 = j;
            double norm_a = s1->normalize(l1, m1, n1);

            for (int k = 0; k <= am2; k++) {
                int l2 = am2 - k;
                for (int l = 0; l <= k; l++) {
                    int m2 = k - l;
                    int n2 = l;
                    double norm_b = s2->normalize(l2, m2, n2);

                    for (int m = 0; m <= am3; m++) {
                        int l3 = am3 - m;
                        for (int n = 0; n <= m; n++) {
                            int m3 = m - n;
                            int n3 = n;
                            double norm_c = s3->normalize(l3, m3, n3);

                            for (int o = 0; o <= am4; o++) {
                                int l4 = am4 - o;
                                for (int p = 0; p <= o; p++) {
                                    int m4 = o - p;
                                    int n4 = p;
                                    double norm_d = s4->normalize(l4, m4, n4);

                                    // printf("normalization %f %f %f %f\n", norm_a, norm_b, norm_c, norm_d);
                                    for (int chunk=0; chunk < nchunk; ++chunk) {
                                        source_[nprim+(chunk*length)] *= norm_a * norm_b * norm_c * norm_d;
                                    }
                                    nprim++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("Angular momentum normalization");
#endif
#endif
}

void TwoBodyAOInt::permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24)
{
#ifdef MINTS_TIMER
    timer_on("Permute target");
#endif
    const GaussianShell& s1 = bs1_->shell(sh1);
    const GaussianShell& s2 = bs2_->shell(sh2);
    const GaussianShell& s3 = bs3_->shell(sh3);
    const GaussianShell& s4 = bs4_->shell(sh4);

    int nbf1, nbf2, nbf3, nbf4;
    if(force_cartesian_){
        nbf1 = s1.ncartesian();
        nbf2 = s2.ncartesian();
        nbf3 = s3.ncartesian();
        nbf4 = s4.ncartesian();
    }
    else{
        nbf1 = s1.nfunction();
        nbf2 = s2.nfunction();
        nbf3 = s3.nfunction();
        nbf4 = s4.nfunction();
    }

    if (!p13p24) {
        if (p12) {
            if (p34) {
                permute_1234_to_2143(s, t, nbf1, nbf2, nbf3, nbf4);
            } else {
                permute_1234_to_2134(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        } else {
            permute_1234_to_1243(s, t, nbf1, nbf2, nbf3, nbf4);
        }
    } else {
        if (p12) {
            if (p34) {
                permute_1234_to_4321(s, t, nbf1, nbf2, nbf3, nbf4);
            } else {
                permute_1234_to_4312(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        } else {
            if (p34) {
                permute_1234_to_3421(s, t, nbf1, nbf2, nbf3, nbf4);
            } else {
                permute_1234_to_3412(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("Permute target");
#endif
}

void TwoBodyAOInt::permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf1;
    int f2=nbf2;
    int f3=nbf4;
    int f4=nbf3;
    for (int bf1=0; bf1<f1; bf1++) {
        for (int bf2=0; bf2<f2; bf2++) {
            for (int bf4=0; bf4<f4; bf4++) {
                for (int bf3=0; bf3<f3; bf3++) {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf2;
    int f2=nbf1;
    int f3=nbf3;
    int f4=nbf4;
    for (int bf2=0; bf2<f2; bf2++) {
        for (int bf1=0; bf1<f1; bf1++) {
            for (int bf3=0; bf3<f3; bf3++) {
                for (int bf4=0; bf4<f4; bf4++) {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf2;
    int f2=nbf1;
    int f3=nbf4;
    int f4=nbf3;
    for (int bf2=0; bf2<f2; bf2++) {
        for (int bf1=0; bf1<f1; bf1++) {
            for (int bf4=0; bf4<f4; bf4++) {
                for (int bf3=0; bf3<f3; bf3++) {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf3;
    int f2=nbf4;
    int f3=nbf1;
    int f4=nbf2;
    for (int bf3=0; bf3<f3; bf3++) {
        for (int bf4=0; bf4<f4; bf4++) {
            for (int bf1=0; bf1<f1; bf1++) {
                for (int bf2=0; bf2<f2; bf2++) {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf4;
    int f2=nbf3;
    int f3=nbf1;
    int f4=nbf2;
    for (int bf3=0; bf3<f3; bf3++) {
        for (int bf4=0; bf4<f4; bf4++) {
            for (int bf2=0; bf2<f2; bf2++) {
                for (int bf1=0; bf1<f1; bf1++) {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf3;
    int f2=nbf4;
    int f3=nbf2;
    int f4=nbf1;
    for (int bf4=0; bf4<f4; bf4++) {
        for (int bf3=0; bf3<f3; bf3++) {
            for (int bf1=0; bf1<f1; bf1++) {
                for (int bf2=0; bf2<f2; bf2++) {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf4;
    int f2=nbf3;
    int f3=nbf2;
    int f4=nbf1;
    for (int bf4=0; bf4<f4; bf4++) {
        for (int bf3=0; bf3<f3; bf3++) {
            for (int bf2=0; bf2<f2; bf2++) {
                for (int bf1=0; bf1<f1; bf1++) {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::pure_transform(int sh1, int sh2, int sh3, int sh4, int nchunk,
                                  bool copy_to_source)
{
#ifdef MINTS_TIMER
    timer_on("Pure transformation");
#endif
    const GaussianShell& s1 = bs1_->shell(sh1);
    const GaussianShell& s2 = bs2_->shell(sh2);
    const GaussianShell& s3 = bs3_->shell(sh3);
    const GaussianShell& s4 = bs4_->shell(sh4);

    // Get the transforms from the basis set
    SphericalTransformIter trans1(*integral()->spherical_transform(s1.am()));
    SphericalTransformIter trans2(*integral()->spherical_transform(s2.am()));
    SphericalTransformIter trans3(*integral()->spherical_transform(s3.am()));
    SphericalTransformIter trans4(*integral()->spherical_transform(s4.am()));

    // Get the angular momentum for each shell
    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();

    // Get number of Cartesian functions for each shell
    int nao1 = s1.ncartesian();
    int nao2 = s2.ncartesian();
    int nao3 = s3.ncartesian();
    int nao4 = s4.ncartesian();

    int nbf1 = s1.nfunction();
    int nbf2 = s2.nfunction();
    int nbf3 = s3.nfunction();
    int nbf4 = s4.nfunction();

    // Get if each shell has pure functions
    bool is_pure1 = s1.is_pure();
    bool is_pure2 = s2.is_pure();
    bool is_pure3 = s3.is_pure();
    bool is_pure4 = s4.is_pure();

    for (int ichunk=0; ichunk < nchunk; ++ichunk) {
        // Compute the offset in source_, and target
        size_t sourcechunkoffset = ichunk * (nao1 * nao2 * nao3 * nao4);
        double *source1, *target1;
        double *source2, *target2;
        double *source3, *target3;
        double *source4, *target4;
        double *source = source_+sourcechunkoffset;
        double *target = target_+sourcechunkoffset;
        double *tmpbuf = tformbuf_;

        int transform_index = 8*is_pure1 + 4*is_pure2 + 2*is_pure3 + is_pure4;
        switch (transform_index) {
            case 0:
                break;

            case 1:
                source4 = source;
                target4 = target;
                break;

            case 2:
                source3 = source;
                target3 = target;
                break;

            case 3:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = target;
                break;

            case 4:
                source2 = source;
                target2 = target;
                break;

            case 5:
                source4 = source;
                target4 = tmpbuf;
                source2 = tmpbuf;
                target2 = target;
                break;

            case 6:
                source3 = source;
                target3 = tmpbuf;
                source2 = tmpbuf;
                target2 = target;
                break;

            case 7:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = source;
                source2 = source;
                target2 = target;
                break;

            case 8:
                source1 = source;
                target1 = target;
                break;

            case 9:
                source4 = source;
                target4 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;

            case 10:
                source3 = source;
                target3 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;

            case 11:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = source;
                source1 = source;
                target1 = target;
                break;

            case 12:
                source2 = source;
                target2 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;

            case 13:
                source4 = source;
                target4 = tmpbuf;
                source2 = tmpbuf;
                target2 = source;
                source1 = source;
                target1 = target;
                break;

            case 14:
                source3 = source;
                target3 = tmpbuf;
                source2 = tmpbuf;
                target2 = source;
                source1 = source;
                target1 = target;
                break;

            case 15:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = source;
                source2 = source;
                target2 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;
        }

        size_t size = nbf1 * nbf2 * nbf3 * nbf4;
        if (is_pure4) {
            transform2e_4(am4, trans4, source4, target4, nao1*nao2*nao3,nao4);
//            size *= nbf4;
        }
        if (is_pure3) {
            transform2e_3(am3, trans3, source3, target3, nao1*nao2,nao3,nbf4);
//            size *= nbf3;
        }
        if (is_pure2) {
            transform2e_2(am2, trans2, source2, target2, nao1,nao2,nbf3*nbf4);
//            size *= nbf2;
        }
        if (is_pure1) {
            transform2e_1(am1, trans1, source1, target1, nbf2*nbf3*nbf4);
//            size *= nbf1;
        }

        // The permute indices routines depend on the integrals being in source_
        if (copy_to_source && (is_pure1 || is_pure2 || is_pure3 || is_pure4))
            memcpy(source, target, size * sizeof(double));
    }
#ifdef MINTS_TIMER
    timer_off("Pure transformation");
#endif
}

static void transform2e_1(int am, SphericalTransformIter& sti, double *s, double *t, int njkl)
{
    memset(t,0,INT_NPURE(am)*njkl*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex()*njkl;
        double *tptr = t + sti.pureindex()*njkl;
        double coef = sti.coef();

//        outfile->Printf( "2e_1: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int jkl=0; jkl<njkl; jkl++)
            *(tptr++) += coef * *(sptr++);
    }
}

static void transform2e_2(int am, SphericalTransformIter& sti, double *s, double *t, int ni, int nj, int nkl)
{
    int sj = INT_NPURE(am);
    const int sjkl = nj*nkl;
    const int tjkl = sj*nkl;

    memset(t,0,ni*tjkl*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex()*nkl;
        double *tptr = t + sti.pureindex()*nkl;
        double coef = sti.coef();

//        outfile->Printf( "2e_2: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int i=0; i<ni; i++,sptr+=sjkl,tptr+=tjkl) {
            for(int kl=0; kl<nkl; kl++)
                tptr[kl] += coef * sptr[kl];
        }
    }
}

static void transform2e_3(int am, SphericalTransformIter& sti, double *s, double *t, int nij, int nk, int nl)
{
    int sk = INT_NPURE(am);
    const int skl = nk*nl;
    const int tkl = sk*nl;

    memset(t,0,nij*tkl*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex()*nl;
        double *tptr = t + sti.pureindex()*nl;
        // printf("cartindex = %d, pureindex = %d\n", sti.cartindex(), sti.pureindex());

//        outfile->Printf( "2e_3: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        double coef = sti.coef();
        for(int ij=0; ij<nij; ij++,sptr+=skl,tptr+=tkl) {
            for(int l=0; l<nl; l++)
                tptr[l] += coef * sptr[l];
        }
    }
}

// am => angular momentum of l
// sti => spherical tranformation iterator
// s => source integrals buffer
// t => target buffer
// nijk => how many i, j, k combinations are there?
// nl => how man l's are there?
static void transform2e_4(int am, SphericalTransformIter& sti, double *s, double *t, int nijk, int nl)
{
    // Protect ourselves
    const int sl = nl;
    const int tl = INT_NPURE(am);

    // Clear out target memory
    memset(t, 0, nijk*tl*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        // Starting point in source and target buffers
        double *sptr = s + sti.cartindex();
        double *tptr = t + sti.pureindex();

//        outfile->Printf( "2e_4: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        // What's the coefficient we're using
        double coef = sti.coef();
        for (int ijk=0; ijk<nijk; ++ijk) {
            // Add contribution of the source to the target
            *(tptr) += coef * *(sptr);

            // skip ahead to the next ijk
            sptr += sl;
            tptr += tl;
        }
    }
}

