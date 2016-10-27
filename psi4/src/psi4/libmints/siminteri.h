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

#ifndef _psi4_libmints_siminteri_h
#define _psi4_libmints_siminteri_h

#include "psi4/libmints/twobody.h"

#include <simint/simint.h>

namespace psi {


class SimintTwoElectronInt : public TwoBodyAOInt
{
    private:
        typedef std::vector<simint_multi_shellpair> ShellPairVec;
        typedef std::vector<simint_shell> ShellVec;

        static simint_shell PsiShellToSimint(const GaussianShell & s);

        static void MultishellVectorDeleter(ShellPairVec * mpv);
        static void ShellVectorDeleter(ShellVec * mpv);

        static std::shared_ptr<ShellPairVec> 
               CreateShellPairs(const ShellVec & bs1, const ShellVec & bs2);

        std::shared_ptr<ShellVec> CreateShellVec(const BasisSet & bs);

        std::shared_ptr<ShellPairVec>
        initialize_shell_pair_single_(const std::vector<ShellPairBlock> & vsh,
                                      const ShellVec & shell1, const ShellVec & shell2);


        size_t multi_batchsize_;
        double * allwork_;
        double * sharedwork_;

        // plain simint shells
        std::shared_ptr<ShellVec> shells1_;
        std::shared_ptr<ShellVec> shells2_;
        std::shared_ptr<ShellVec> shells3_;
        std::shared_ptr<ShellVec> shells4_;

        // my shell pairs
        std::shared_ptr<ShellPairVec> single_spairs_bra_;
        std::shared_ptr<ShellPairVec> single_spairs_ket_;
        std::shared_ptr<ShellPairVec> multi_spairs_bra_;
        std::shared_ptr<ShellPairVec> multi_spairs_ket_;


        // a buffer for shell pairs
        simint_multi_shellpair P_, Q_;

        // number of shells in each basis set
        int nsh1_, nsh2_, nsh3_, nsh4_;

        void create_blocks(void);

    public:
        SimintTwoElectronInt(const IntegralFactory * integral, int deriv=0, bool use_shell_pairs=false);
        virtual ~SimintTwoElectronInt();

        virtual size_t compute_shell(const AOShellCombinationsIterator&) override;

        virtual size_t compute_shell(int, int, int, int) override;

        virtual void
        compute_shell_blocks(int shellpair1, int shellpair2,
                             int npair1 = -1, int npair2 = -1) override;

        virtual size_t compute_shell_deriv1(int, int, int, int);

        virtual size_t compute_shell_deriv2(int, int, int, int);
};


class SimintERI : public SimintTwoElectronInt
{
    public:
        SimintERI(const IntegralFactory * integral, int deriv=0, bool use_shell_pairs=false);
};



}

#endif
