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

#include "pointgrp.h"
#include "petitelist.h"
#include "sobasis.h"
#include "molecule.h"
#include "basisset.h"
#include "gshell.h"
#include "dimension.h"

#include "matrix.h"
#include "psi4/psi4-dec.h"
#include <cstdio>
#include "psi4/libparallel/ParallelPrinter.h"
;
using namespace psi;

///////////////////////////////////////////////////////////////////////////////

SOTransform::SOTransform()
{
    naoshell_allocated = 0;
    naoshell = 0;
    aoshell = 0;
}

SOTransform::~SOTransform()
{
    delete[] aoshell;
}

void SOTransform::set_naoshell(int n)
{
    naoshell = 0;
    delete[] aoshell;
    naoshell_allocated = n;
    aoshell = new SOTransformShell[n];
}

void SOTransform::add_transform(int aoshellnum, int irrep,
                                double coef, int aofunc, int sofunc)
{
//    outfile->Printf( "SOTransform::add_transform(aoshellnum = %d, irrep = %d, coef = %lf, aofunc = %d, sofunc = %d)\n", aoshellnum, irrep, coef, aofunc, sofunc);

    int i;
    for (i=0; i<naoshell; i++) {
        if (aoshell[i].aoshell == aoshellnum) break;
    }
    if (i>=naoshell_allocated) {
        throw PSIEXCEPTION("SOTransform::add_transform: allocation too small");
    }
    aoshell[i].add_func(irrep,coef,aofunc,sofunc);
    aoshell[i].aoshell = aoshellnum;
    if (i==naoshell) naoshell++;
}

///////////////////////////////////////////////////////////////////////////////

AOTransform::AOTransform()
{
    for (int h=0; h<8; h++)
        nfuncpi[h] = 0;
}

AOTransform::~AOTransform()
{
}

void AOTransform::add_transform(int irrep,
                                double coef, int aofunc, int sofunc)
{
    soshell.push_back(AOTransformFunction(coef, aofunc, sofunc, irrep));
    soshellpi[irrep].push_back(AOTransformFunction(coef, aofunc, sofunc, irrep));
    nfuncpi[irrep]++;
}

///////////////////////////////////////////////////////////////////////////////

SOTransformShell::SOTransformShell()
{
    nfunc = 0;
    func = 0;
}

SOTransformShell::~SOTransformShell()
{
    if (func)
        delete[] func;
}

void SOTransformShell::add_func(int irrep, double coef, int aofunc, int sofunc)
{
    SOTransformFunction *newfunc = new SOTransformFunction[nfunc+1];
    for (int i=0; i<nfunc; i++) newfunc[i] = func[i];
    delete[] func;
    func = newfunc;
    func[nfunc].irrep = irrep;
    func[nfunc].coef = coef;
    func[nfunc].aofunc = aofunc;
    func[nfunc].sofunc = sofunc;
    nfunc++;
}

///////////////////////////////////////////////////////////////////////////////

SOBasisSet::SOBasisSet(const std::shared_ptr<BasisSet> &basis, const IntegralFactory *integral)
    : basis_(basis), integral_(integral)
{
    init();
}

SOBasisSet::SOBasisSet(const std::shared_ptr<BasisSet> &basis, const std::shared_ptr<IntegralFactory> &integral)
    : basis_(basis), integral_(integral.get())
{
    init();
}

void SOBasisSet::init()
{
    int i,j,k;

    std::shared_ptr<Molecule> mol = basis_->molecule();

    CharacterTable ct = mol->point_group()->char_table();
    nirrep_ = ct.nirrep();

    // count the number of so shells
    nshell_ = 0;
    for (i=0; i<mol->nunique(); i++) {
        nshell_ += basis_->nshell_on_center(mol->unique(i));
    }

    //=----- Begin debug printing -----=
//    outfile->Printf( "SOBasis:\n");
//    outfile->Printf( "nshell_ = %d\n", nshell_);
    //=-----  End debug printing  -----=

    // Allocate memory for unique shell to am
    ushell_am_ = new int[nshell_];

    // map each ao shell to an so shell
    int *aoshell_to_soshell = new int[basis_->nshell()];
    int soshell = 0;
    for (i=0; i<mol->nunique(); i++) {
        for (j=0; j<basis_->nshell_on_center(mol->unique(i)); j++) {
            for (k=0; k<mol->nequivalent(i); k++) {
                int aoshell = basis_->shell_on_center(mol->equivalent(i,k),j);
                aoshell_to_soshell[aoshell] = soshell;
//                outfile->Printf( "i = %d j = %d k = %d aoshell = %d soshell = %d, mol->equivalent = %d\n",
//                        i, j, k, aoshell, soshell, mol->equivalent(i,k));
            }

            // For each so shell obtain its angular momentum
            ushell_am_[soshell] = basis_->shell(mol->unique(i), j).am();

            soshell++;
        }
    }

    //=----- Begin debug printing -----=
//    outfile->Printf( "Final aoshell_to_soshell:\n");
//    for (i = 0; i < basis_->nshell(); ++i) {
//        outfile->Printf( "aoshell_to_soshell[%d] = %d\n", i, aoshell_to_soshell[i]);
//    }
    //=-----  End debug printing  -----=

    ncomp_ = new int[nirrep_];
    for (i=0; i<nirrep_; i++) {
        ncomp_[i] = ct.gamma(i).degeneracy();
        if (ncomp_[i] != 1) {
            throw PSIEXCEPTION("SOBasis::SOBasis: not tested for degenerate point groups");
        }
    }

    naofunc_ = new int[nshell_];
    memset(naofunc_, 0, sizeof(int)*nshell_);

    nfunc_ = new int*[nshell_];
    funcoff_ = new int*[nshell_];
    for (i=0; i<nshell_; i++) {
        nfunc_[i] = new int[nirrep_];
        funcoff_[i] = new int[nirrep_];
        for (j=0; j<nirrep_; j++) {
            nfunc_[i][j] = 0;
        }
    }

    bool include_pure_transform = true;

    petite_ = std::shared_ptr<PetiteList>(new PetiteList(basis_, integral_, include_pure_transform));

//    petite_->print();


    int nblocks = petite_->nblocks();
    SO_block *soblocks(petite_->compute_aotoso_info());

//    for (i=0; i<nblocks; ++i) {
//        outfile->Printf( "soblock[%d]\n", i);
//        soblocks[i].print("");
//    }

    // == Begin forming (A|S)OTransform array ==
    sotrans_ = new SOTransform[nshell_];   // nshell_ is symmetry unique shells
    aotrans_ = new AOTransform[basis_->nshell()]; // we need the ao shell number here

    for (i=0; i<nblocks; i++) {
        for (j=0; j<soblocks[i].len; j++) {
            if (soblocks[i].so[j].length == 0) continue;
            int bfn0 = soblocks[i].so[j].cont[0].bfn;
            int aoshell0 = include_pure_transform ?
                        basis_->ao_to_shell(bfn0) : basis_->function_to_shell(bfn0);
            int soshell0 = aoshell_to_soshell[aoshell0];
            int atom0 = basis_->shell_to_center(aoshell0);
            int nequiv0 = mol->nequivalent(mol->atom_to_unique(atom0));
            sotrans_[soshell0].set_naoshell(nequiv0);
//            outfile->Printf( "i = %d j = %d bfn0 = %d aoshell0 = %d soshell0 = %d atom0 = %d nequiv0 = %d\n", i, j, bfn0, aoshell0, soshell0, atom0, nequiv0);
        }
    }

    int nfuncall = 0;
    for (i=0; i<nblocks; i++) {
        int irrep = ct.which_irrep(i);
        for (j=0; j<soblocks[i].len; j++) {
            if (soblocks[i].so[j].length == 0) continue;
            int bfn0 = soblocks[i].so[j].cont[0].bfn;
            int aoshell0 = include_pure_transform ?
                        basis_->ao_to_shell(bfn0) : basis_->function_to_shell(bfn0);
            int soshell0 = aoshell_to_soshell[aoshell0];
            int sofunc = nfunc_[soshell0][irrep];

            int naofunc = include_pure_transform ? basis_->shell(aoshell0).ncartesian() : basis_->shell(aoshell0).nfunction();
            if (naofunc_[soshell0] && (naofunc_[soshell0] != naofunc)) {
                throw PSIEXCEPTION("SOBasis::SOBasis: mismatch in naofunc");
            }
            naofunc_[soshell0] = naofunc;

            nfunc_[soshell0][irrep]++;
            nfuncall++;

            for (k=0; k<soblocks[i].so[j].length; k++) {
                int bfn = soblocks[i].so[j].cont[k].bfn;
                double coef = soblocks[i].so[j].cont[k].coef;
                int aoshell = include_pure_transform ? basis_->ao_to_shell(bfn) : basis_->function_to_shell(bfn);
                int aoshellfunc = bfn - (include_pure_transform ?
                            basis_->shell_to_ao_function(aoshell) : basis_->shell_to_basis_function(aoshell));
                int soshell = aoshell_to_soshell[aoshell];

                if (soshell != soshell0) {
                    throw PSIEXCEPTION("SOBasis::SOBasis: shell changed");
                }

                sotrans_[soshell].add_transform(aoshell, irrep, coef, aoshellfunc, sofunc);
                aotrans_[aoshell].add_transform(irrep, coef, aoshellfunc, sofunc);
            }
        }
    }

    // == End forming (A|S)OTransform array ==

    if (nfuncall != basis_->nbf()) {
        throw PSIEXCEPTION("SOBasis::SOBasis: miscounted number of functions");
    }

    for (i=0; i<nshell_; i++) {
        funcoff_[i][0] = 0;
        for (j=1; j<nirrep_; j++) {
            funcoff_[i][j] = funcoff_[i][j-1] + nfunc_[i][j-1];
//            outfile->Printf( "funcoff_[%d][%d] = %d\n", i, j, funcoff_[i][j]);
        }
    }

    for(int i=0; i < basis_->nshell(); ++i){
        int usoshell = aoshell_to_soshell[i];
        aotrans_[i].add_offsets(nirrep_, funcoff_[usoshell]);
    }

    delete[] aoshell_to_soshell;
    delete[] soblocks;

    func_ = new int[nshell_];
    irrep_ = new int[basis_->nbf()];
    func_within_irrep_ = new int[basis_->nbf()];
    nfunc_in_irrep_ = new int[nirrep_];

    for (i=0; i<nirrep_; i++) nfunc_in_irrep_[i] = 0;

    if (nshell_) {
        func_[0] = 0;
        for (i=1; i<nshell_; i++) {
            func_[i] = func_[i-1] + nfunction(i-1);
//            outfile->Printf( "func_[%d] = %d\n", i, func_[i]);
        }
        int ibasis_ = 0;
        for (i=0; i<nshell_; i++) {
            for (j=0; j<nirrep_; j++) {
                for (k=0; k<nfunc_[i][j]; k++,ibasis_++) {
                    irrep_[ibasis_] = j;
                    func_within_irrep_[ibasis_] = nfunc_in_irrep_[j]++;
//                    outfile->Printf( "irrep_[%d] = %d func_within_irrep_[%d] = %d\n", ibasis_, j, ibasis_, func_within_irrep_[ibasis_]);
                }
            }
        }
    }

    // Create a map that has a key/value pair
    // The key is the angular momentum function of the SO shell arranged in decending order
    // The value is the actual shell number
    typedef std::pair<int, int> am_to_so_shell_pair;
    std::multimap< int, int, std::less<int> > am_to_so_shell_list;
    for(int i=0; i < nshell_; i++) {
        am_to_so_shell_list.insert(am_to_so_shell_pair(naofunction(i), i));
        //std::cout << "naofunctions(" << i << ") = " << naofunction(i) << std::endl;
    }
    // This puts the sorted SO shell values into the sorted_so_shell_list_ vector,
    // which can be used by the integral iterator to look up the value of the sorted shells
    std::multimap< int, int, std::less<int> >::iterator it;
    for (it=am_to_so_shell_list.begin(); it != am_to_so_shell_list.end(); it++) {
        //std::cout << "sorted shell size = " << it->first <<
        //        "\t, which belongs to shell number " << it->second << std::endl;
        sorted_so_shell_list_.push_back(it->second);
    }
//    print();
}

SOBasisSet::~SOBasisSet()
{
    for (int i=0; i<nshell_; i++) {
        delete[] nfunc_[i];
        delete[] funcoff_[i];
    }
    delete[] nfunc_;
    delete[] funcoff_;
    delete[] naofunc_;
    delete[] ncomp_;
    delete[] sotrans_;
    delete[] aotrans_;
    delete[] func_;
    delete[] irrep_;
    delete[] func_within_irrep_;
    delete[] nfunc_in_irrep_;
    delete[] ushell_am_;
}

int SOBasisSet::max_nfunction_in_shell() const
{
    int maxn = 0;
    for (int i=0; i<nshell_; i++) {
        int n = nfunction(i);
        if (n > maxn) maxn = n;
    }
    return maxn;
}

int SOBasisSet::nfunction(int ishell) const
{
    int n=0;
    for (int i=0; i<nirrep_; i++) {
        n += nfunc_[ishell][i];
    }
    return n;
}

void SOBasisSet::print(std::string out) const
{
    int i,j,k;
    std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
          std::shared_ptr<OutFile>(new OutFile(out)));
    printer->Printf( "  SOBasis:\n");
    printer->Printf( "    nshell(SO) = %d\n", nshell_);
    printer->Printf( "    nirrep = %d\n", nirrep_);

    printer->Printf( "    ncomp = [");
    for (i=0; i<nirrep_; i++)
        printer->Printf( " %3d", ncomp_[i]);
    printer->Printf( " ]\n");

    printer->Printf( "    nfunc:\n");
    for (i=0; i<nshell_; i++) {
        printer->Printf( "      %3d:", i);
        for (j=0; j<nirrep_; j++)
            printer->Printf( "  %3d", nfunc_[i][j]);
        printer->Printf( "\n");
    }

    printer->Printf( "    irrep             = [");
    for (i=0; i<basis_->nbf(); ++i) {
        printer->Printf( " %4d", irrep_[i]);
    }
    printer->Printf( "]\n");

    printer->Printf( "    func              = [");
    for (i=0; i<nshell_; ++i) {
        printer->Printf( " %4d", func_[i]);
    }
    printer->Printf( "]\n");

    printer->Printf( "    func_within_irrep = [");
    for (i=0; i<basis_->nbf(); ++i) {
        printer->Printf( " %4d", func_within_irrep_[i]);
    }
    printer->Printf( "]\n");

    printer->Printf( "    nfunc_in_irrep    = [");
    for (i=0; i<nirrep_; ++i) {
        printer->Printf( " %4d", nfunc_in_irrep_[i]);
    }
    printer->Printf( "]\n");

    printer->Printf( "    funcoff           = [\n");
    for (i=0; i<nshell_; i++) {
        printer->Printf( "      %3d:", i);
        for (j=0; j<nirrep_; j++)
            printer->Printf( "  %3d", funcoff_[i][j]);
        printer->Printf( "\n");
    }

    printer->Printf( "    sotransform:\n");
    for (i=0; i<nshell_; i++) {
        if (i>0) printer->Printf( "\n");
        for (j=0; j<sotrans_[i].naoshell; j++) {
            for (k=0; k<sotrans_[i].aoshell[j].nfunc; k++) {
                printer->Printf( "      SO(%3d %2d %d [%2d]) += %12.8f * AO(%3d %2d)\n",
                        i,
                        sotrans_[i].aoshell[j].func[k].sofunc,
                        sotrans_[i].aoshell[j].func[k].irrep,
                        function_offset_within_shell(
                            i, sotrans_[i].aoshell[j].func[k].irrep)
                        + sotrans_[i].aoshell[j].func[k].sofunc,
                        sotrans_[i].aoshell[j].func[k].coef,
                        sotrans_[i].aoshell[j].aoshell,
                        sotrans_[i].aoshell[j].func[k].aofunc);
            }
        }
    }

    printer->Printf( "    aotransform:\n");
    for (i=0; i<basis_->nshell(); ++i) {
        if (i>0) printer->Printf( "\n");
        for (j=0; j<(int)aotrans_[i].soshell.size(); ++j) {
            printer->Printf( "      AO(%3d) sofunc %d aofunc %d irrep %d coef %12.8f\n",
                    i,
                    aotrans_[i].soshell[j].sofunc,
                    aotrans_[i].soshell[j].aofunc,
                    aotrans_[i].soshell[j].irrep,
                    aotrans_[i].soshell[j].coef);
        }
    }
}

std::shared_ptr<BasisSet> SOBasisSet::basis() const
{
    return basis_;
}

Dimension SOBasisSet::dimension() const
{
    std::shared_ptr<PetiteList> petite = std::shared_ptr<PetiteList>(new PetiteList(basis_, integral_));
    return petite->SO_basisdim();
}

const std::shared_ptr<PetiteList> SOBasisSet::petite_list() const
{
    return petite_;
}
