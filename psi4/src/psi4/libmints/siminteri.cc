/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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


#include <algorithm>

#include "psi4/libmints/siminteri.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/cartesianiter.h"
#include "psi4/libmints/integral.h"

#define SIMINT_SCREEN_TOL 0.0
#define SIMINT_SCREEN 0

using ShellVec = psi::SimintTwoElectronInt::ShellVec;
using ShellPairVec = psi::SimintTwoElectronInt::ShellPairVec;

namespace psi {


// Some helpers
static void shell_vector_deleter_(ShellVec * mpv)
{
    for(auto & mp : *mpv)
       simint_free_shell(&mp);
    mpv->clear();
}


static void multishell_vector_deleter_(ShellPairVec * mpv)
{
    for(auto & mp : *mpv)
       simint_free_multi_shellpair(&mp);
    mpv->clear();
}


static simint_shell psi_shell_to_simint_(const GaussianShell & s)
{
    size_t nprim = s.nprimitive();
    simint_shell newgs;

    simint_initialize_shell(&newgs);
    simint_allocate_shell(nprim, &newgs);
    newgs.am = s.am();
    std::copy(s.exps(), s.exps()+nprim, newgs.alpha);
    std::copy(s.coefs(), s.coefs()+nprim, newgs.coef);
    newgs.nprim = static_cast<int>(nprim);
    newgs.x = s.center()[0];
    newgs.y = s.center()[1];
    newgs.z = s.center()[2];

    return newgs;
}


static std::shared_ptr<ShellVec>
create_shell_vec_(const BasisSet & bs)
{
    std::shared_ptr<ShellVec> shells = std::shared_ptr<ShellVec>(new ShellVec, &shell_vector_deleter_);
    for(size_t i = 0; i < bs.nshell(); i++)
        shells->push_back(psi_shell_to_simint_(bs.shell(i)));

    return shells;
} 

static std::shared_ptr<ShellPairVec>
create_shell_pair_(const ShellVec & bs1, const ShellVec & bs2)
{
    auto spairs = std::shared_ptr<ShellPairVec>(new ShellPairVec, &multishell_vector_deleter_);

    size_t nshell1 = bs1.size();
    size_t nshell2 = bs2.size();
    spairs->resize(nshell1 * nshell2);

    for(size_t i = 0; i < nshell1; i++)
    {
        const simint_shell & sh_i = bs1[i];
        for(size_t j = 0; j < nshell2; j++)
        {
            const simint_shell & sh_j = bs2[j];
            simint_multi_shellpair P;
            simint_initialize_multi_shellpair(&P);
            simint_create_multi_shellpair(1, &sh_i, 1, &sh_j, &P, SIMINT_SCREEN);
            (*spairs)[i*nshell2+j] = P;
        }
    }

    return spairs;
}

static ShellPairVec
create_multi_shellpair_(const std::vector<ShellPairBlock> & vsh,
                        const ShellVec & shell1,
                        const ShellVec & shell2)
{
    ShellPairVec ret;

    // copying simint shells will copy the pointers. That is ok - we won't
    // do anything special to delete them
    for(const auto & spairs : vsh)
    {
        std::vector<simint_shell> simint_shells;
        
        for(const auto & s : spairs)
        {
            simint_shells.push_back(shell1[s.first]);
            simint_shells.push_back(shell2[s.second]);
        }

        // spairs.size() should be simint_shells.size()/2
        simint_multi_shellpair P;
        simint_initialize_multi_shellpair(&P);
        simint_create_multi_shellpair2(spairs.size(), simint_shells.data(), &P, SIMINT_SCREEN);
        ret.push_back(P);
    }

    return ret;
}

static std::shared_ptr<ShellPairVec>
create_shared_multi_shellpair_(const std::vector<ShellPairBlock> & vsh,
                               const ShellVec & shell1,
                               const ShellVec & shell2)
{
    auto vec = create_multi_shellpair_(vsh, shell1, shell2);
    return std::shared_ptr<ShellPairVec>(new ShellPairVec(std::move(vec)));
}




////////////////////////////////////////////////
// Actual class code starts here
////////////////////////////////////////////////

SimintTwoElectronInt::SimintTwoElectronInt(const IntegralFactory * integral, int deriv, bool use_shell_pairs)
    : TwoBodyAOInt(integral, deriv)
{
    maxam_ = std::max( std::max(basis1()->max_am(), basis2()->max_am()),
                       std::max(basis3()->max_am(), basis4()->max_am()) );

    if(maxam_ > SIMINT_OSTEI_MAXAM)
        throw PSIEXCEPTION("ERROR - SIMINT CANNOT HANDLE AM THIS HIGH\n");
    if(deriv_ > SIMINT_OSTEI_MAXDER)
        throw PSIEXCEPTION("ERROR - SIMINT CANNOT HANDLE DERIVATIVES THIS HIGH\n");

    // initialize the library
    // It is safe to call this multiple times (it is only thread unsafe)
    simint_init();

    batchsize_ = 32;

    size_t size = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                  INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    size_t fullsize = size * batchsize_;
    allwork_size_ = sizeof(double) * (fullsize * 2 + size);

    const size_t simint_workmem = simint_ostei_workmem(deriv_, maxam_);
    sharedwork_ = (double *)SIMINT_ALLOC(simint_workmem);
    allwork_ = (double *)SIMINT_ALLOC(allwork_size_);

    target_full_ = allwork_;
    source_full_ = allwork_ + fullsize;
    tformbuf_ = allwork_ + 2*fullsize;

    target_ = target_full_;
    source_ = source_full_;

    // build plain shells
    shells1_ = create_shell_vec_(*original_bs1_);

    if(original_bs2_ == original_bs1_)
        shells2_ = shells1_;
    else
        shells2_ = create_shell_vec_(*original_bs2_);

    if(original_bs3_ == original_bs1_)
        shells3_ = shells1_;
    else if(original_bs3_ == original_bs2_)
        shells3_ = shells2_;
    else
        shells3_ = create_shell_vec_(*original_bs3_);

    if(original_bs4_ == original_bs1_)
        shells4_ = shells1_;
    else if(original_bs4_ == original_bs2_)
        shells4_ = shells2_;
    else if(original_bs4_ == original_bs3_)
        shells4_ = shells3_;
    else
        shells4_ = create_shell_vec_(*original_bs4_);
    

    bra_same_ = (original_bs1_ == original_bs2_);
    ket_same_ = (original_bs3_ == original_bs4_);
    braket_same_ = (original_bs1_ == original_bs3_ &&
                    original_bs2_ == original_bs4_);

    single_spairs_bra_ = create_shell_pair_(*shells1_, *shells2_);

    if(braket_same_)
        single_spairs_ket_ = single_spairs_bra_;
    else
        single_spairs_ket_ = create_shell_pair_(*shells3_, *shells4_);

    create_blocks();
}

SimintTwoElectronInt::SimintTwoElectronInt(const SimintTwoElectronInt & rhs)
    : TwoBodyAOInt(rhs),
      maxam_(rhs.maxam_),
      batchsize_(rhs.batchsize_),
      allwork_size_(rhs.allwork_size_),
      bra_same_(rhs.bra_same_),
      ket_same_(rhs.ket_same_),
      braket_same_(rhs.braket_same_),
      shells1_(rhs.shells1_),
      shells2_(rhs.shells2_),
      shells3_(rhs.shells3_),
      shells4_(rhs.shells4_),
      single_spairs_bra_(rhs.single_spairs_bra_),
      single_spairs_ket_(rhs.single_spairs_ket_),
      multi_spairs_bra_(rhs.multi_spairs_bra_),
      multi_spairs_ket_(rhs.multi_spairs_ket_)
{
    // allocate and reconfigure workspace
    const size_t simint_workmem = simint_ostei_workmem(deriv_, maxam_);
    sharedwork_ = (double *)SIMINT_ALLOC(simint_workmem);
    allwork_ = (double *)SIMINT_ALLOC(allwork_size_);

    // we can do this without storing all the various
    // sizes for the components by simply using the offsets
    // from the object we are copying
    target_full_ = allwork_;
    source_full_ = allwork_ + (rhs.source_full_ - rhs.allwork_);
    tformbuf_ = allwork_ + (rhs.tformbuf_ - rhs.allwork_);
    target_ = target_full_;
    source_ = source_full_;
}

SimintTwoElectronInt::~SimintTwoElectronInt()
{
    simint_finalize();
    SIMINT_FREE(allwork_);
    SIMINT_FREE(sharedwork_);
}


size_t SimintTwoElectronInt::compute_shell(const AOShellCombinationsIterator & shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t SimintTwoElectronInt::compute_shell_deriv1(int, int, int, int)
{
    throw PSIEXCEPTION("ERROR - SIMINT CANNOT HANDLE DERIVATIVES\n");
}

size_t SimintTwoElectronInt::compute_shell_deriv2(int, int, int, int)
{
    throw PSIEXCEPTION("ERROR - SIMINT CANNOT HANDLE DERIVATIVES\n");
}


size_t SimintTwoElectronInt::compute_shell(int sh1, int sh2, int sh3, int sh4)
{
    target_ = target_full_;
    source_ = source_full_;

    const auto nsh2 = original_bs2_->nshell();
    const auto nsh4 = original_bs4_->nshell();

    const auto & shell1 = original_bs1_->shell(sh1);
    const auto & shell2 = original_bs2_->shell(sh2);
    const auto & shell3 = original_bs3_->shell(sh3);
    const auto & shell4 = original_bs4_->shell(sh4);

    bool do_cart = force_cartesian_ || (shell1.is_cartesian() &&
                                        shell2.is_cartesian() &&
                                        shell3.is_cartesian() &&
                                        shell4.is_cartesian());

    int n1, n2, n3, n4;

    if (force_cartesian_) {
        n1 = shell1.ncartesian();
        n2 = shell2.ncartesian();
        n3 = shell3.ncartesian();
        n4 = shell4.ncartesian();
    }
    else {
        n1 = shell1.nfunction();
        n2 = shell2.nfunction();
        n3 = shell3.nfunction();
        n4 = shell4.nfunction();
    }

    curr_buff_size_ = n1 * n2 * n3 * n4;

    // get the precomputed simint_multi_shellpair
    const simint_multi_shellpair * P = &(*single_spairs_bra_)[sh1*nsh2 + sh2];
    const simint_multi_shellpair * Q = &(*single_spairs_ket_)[sh3*nsh4 + sh4];

    // actually compute
    // if we are doing cartesian, put directly in target. Otherwise, put in source
    // and let pure_transform put it in target
    size_t ncomputed = 0;

    if(do_cart)
        ncomputed = simint_compute_eri(P, Q, SIMINT_SCREEN_TOL, sharedwork_, target_);
    else
    {
        ncomputed = simint_compute_eri(P, Q, SIMINT_SCREEN_TOL, sharedwork_, source_);
        pure_transform(sh1, sh2, sh3, sh4, 1, false);
    }

    return ncomputed;
}




void
SimintTwoElectronInt::compute_shell_blocks(int shellpair1, int shellpair2,
                                           int npair1, int npair2)
{
    // set the pointers to the beginning of the work spaces
    target_ = target_full_;
    source_ = source_full_;

    auto vsh12 = blocks12_[shellpair1];
    auto vsh34 = blocks34_[shellpair2];

    // get the basic info
    std::pair<int, int> sh12 = *vsh12.begin();
    std::pair<int, int> sh34 = *vsh34.begin();

    const auto & shell1 = original_bs1_->shell(sh12.first);
    const auto & shell2 = original_bs2_->shell(sh12.second);
    const auto & shell3 = original_bs3_->shell(sh34.first);
    const auto & shell4 = original_bs4_->shell(sh34.second);

    bool do_cart = force_cartesian_ || (shell1.is_cartesian() &&
                                        shell2.is_cartesian() &&
                                        shell3.is_cartesian() &&
                                        shell4.is_cartesian());


    const int ncart1 = shell1.ncartesian();
    const int ncart2 = shell2.ncartesian();
    const int ncart3 = shell3.ncartesian();
    const int ncart4 = shell4.ncartesian();
    const int n1 = shell1.nfunction();
    const int n2 = shell2.nfunction();
    const int n3 = shell3.nfunction();
    const int n4 = shell4.nfunction();

    const int ncart1234 = ncart1 * ncart2 * ncart3 * ncart4;
    const int n1234 = n1 * n2 * n3 * n4;
    curr_buff_size_ = n1234;

    // get P and Q from the precomputed values
    simint_multi_shellpair P = (*multi_spairs_bra_)[shellpair1];
    simint_multi_shellpair Q = (*multi_spairs_ket_)[shellpair2];

    // notice we made a copy. This involves copies
    // of pointers.
    // We can change the number of shell pair without causing
    // problems.
    if(npair1 < 0) npair1 = P.nshell12;
    if(npair2 < 0) npair2 = Q.nshell12;

    // Check that we have enough space now
    size_t nshell1234 = npair1 * npair2;

    if(nshell1234 == 0)
        throw PSIEXCEPTION("No shells passed to calculate");
    if(nshell1234 > batchsize_)
        throw PSIEXCEPTION("Not enough space allocated for that many shells\n");

    P.nshell12_clip = npair1;
    Q.nshell12_clip = npair2;

    // actually compute
    // if we are doing cartesian, put directly in target. Otherwise, put in source
    // and let pure_transform put it in target
    size_t ncomputed = 0;

    if(do_cart)
        ncomputed = simint_compute_eri(&P, &Q, SIMINT_SCREEN_TOL, sharedwork_, target_);
    else
    {
        ncomputed = simint_compute_eri(&P, &Q, SIMINT_SCREEN_TOL, sharedwork_, source_);
        if(!do_cart)
        {
            for(int i = 0; i < npair1; i++)
            {
                auto sh12 = vsh12[i];

                for(int j = 0; j < npair2; j++)
                {
                    auto sh34 = vsh34[j];
                    pure_transform(sh12.first, sh12.second,
                                   sh34.first, sh34.second, 1, false);
    
                    source_ += ncart1234;
                    target_ += n1234;
                }
            }
        }
    }
}

void SimintTwoElectronInt::create_blocks(void)
{
    blocks12_.clear();
    blocks34_.clear();

    // make the significant shell pairs 
    // (form blocks of MU,NU)
    const auto am1 = basis1()->max_am();
    const auto am2 = basis2()->max_am();
    const auto am3 = basis3()->max_am();
    const auto am4 = basis4()->max_am();

    // sort the basis set AM
    std::vector<std::vector<int>> sorted_shells1(am1+1),
                                  sorted_shells2(am2+1),
                                  sorted_shells3(am3+1),
                                  sorted_shells4(am4+1);

    for(int ishell = 0; ishell < basis1()->nshell(); ishell++)
    {
        int am = basis1()->shell(ishell).am();
        sorted_shells1[am].push_back(ishell);
    }

    for(int ishell = 0; ishell < basis2()->nshell(); ishell++)
    {
        int am = basis2()->shell(ishell).am();
        sorted_shells2[am].push_back(ishell);
    }

    for(int ishell = 0; ishell < basis3()->nshell(); ishell++)
    {
        int am = basis3()->shell(ishell).am();
        sorted_shells3[am].push_back(ishell);
    }

    for(int ishell = 0; ishell < basis4()->nshell(); ishell++)
    {
        int am = basis4()->shell(ishell).am();
        sorted_shells4[am].push_back(ishell);
    }


    // form pairs for the bra
    // these aren't batched

    for(int iam = 0; iam <= am1; iam++)
    for(int jam = 0; jam <= am2; jam++)
    {
        for(int ishell : sorted_shells1[iam])
        for(int jshell : sorted_shells2[jam])
        {
            if(!bra_same_ || (bra_same_ && jshell <= ishell) )
                blocks12_.push_back( {{ishell, jshell}} );
        }
    }

    // form pairs for the ket
    for(int iam = 0; iam <= am3; iam++)
    for(int jam = 0; jam <= am4; jam++)
    {
        ShellPairBlock curblock;

        for(int ishell : sorted_shells3[iam])
        for(int jshell : sorted_shells4[jam])
        {
            if(!ket_same_ || (ket_same_ && jshell <= ishell) )
            {
                curblock.push_back({ishell, jshell});
                if(curblock.size() == batchsize_)
                {
                    blocks34_.push_back(curblock);
                    curblock.clear();
                }
            }
        }

        if(curblock.size())
            blocks34_.push_back(std::move(curblock));
    }

    multi_spairs_bra_ = create_shared_multi_shellpair_(blocks12_, *shells1_, *shells2_); 
    multi_spairs_ket_ = create_shared_multi_shellpair_(blocks34_, *shells3_, *shells4_); 

}

SimintERI::SimintERI(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : SimintTwoElectronInt(integral, deriv, use_shell_pairs)
{ }

SimintERI::SimintERI(const SimintERI & rhs)
    : SimintTwoElectronInt(rhs)
{ }





} // close namespace psi
