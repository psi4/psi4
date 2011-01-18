#include "pointgrp.h"
#include "petitelist.h"
#include "sobasis.h"
#include "molecule.h"
#include "basisset.h"
#include "gshell.h"

#include <psi4-dec.h>
#include <cstdio>

using namespace psi;

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

SOTransformShell::SOTransformShell()
{
    nfunc = 0;
    func = 0;
}

SOTransformShell::~SOTransformShell()
{
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

SOBasis::SOBasis(const boost::shared_ptr<BasisSet> &basis, const boost::shared_ptr<IntegralFactory> &integral)
    : basis_(basis)
{
    int i,j,k;

    basis_ = basis;

    boost::shared_ptr<Molecule> mol = basis_->molecule();

    CharacterTable ct = mol->point_group()->char_table();
    nirrep_ = ct.nirrep();

    // count the number of so shells
    nshell_ = 0;
    for (i=0; i<mol->nunique(); i++) {
        nshell_ += basis_->nshell_on_center(mol->unique(i));
    }

    // map each ao shell to an so shell
    int *aoshell_to_soshell = new int[basis_->nshell()];
    int soshell = 0;
    for (i=0; i<mol->nunique(); i++) {
        for (j=0; j<basis_->nshell_on_center(mol->unique(i)); j++) {
            for (k=0; k<mol->nequivalent(i); k++) {
                int aoshell = basis_->shell_on_center(mol->equivalent(i,k),j);
                aoshell_to_soshell[aoshell] = soshell;
            }
            soshell++;
        }
    }

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

    shared_ptr<PetiteList> petite = shared_ptr<PetiteList>(new PetiteList(basis_, integral));

    int nblocks = petite->nblocks();
    SO_block *soblocks = petite->aotoso_info();

    trans_ = new SOTransform[nshell_];
    for (i=0; i<nblocks; i++) {
        for (j=0; j<soblocks[i].len; j++) {
            if (soblocks[i].so[j].length == 0) continue;
            int bfn0 = soblocks[i].so[j].cont[0].bfn;
            int aoshell0 = basis_->function_to_shell(bfn0);
            int soshell0 = aoshell_to_soshell[aoshell0];
            int atom0 = basis_->shell_to_center(aoshell0);
            int nequiv0 = mol->nequivalent(mol->atom_to_unique(atom0));
            trans_[soshell0].set_naoshell(nequiv0);
        }
    }

    int nfuncall = 0;
    for (i=0; i<nblocks; i++) {
        int irrep = ct.which_irrep(i);
        for (j=0; j<soblocks[i].len; j++) {
            if (soblocks[i].so[j].length == 0) continue;
            int bfn0 = soblocks[i].so[j].cont[0].bfn;
            int aoshell0 = basis_->function_to_shell(bfn0);
            int soshell0 = aoshell_to_soshell[aoshell0];
            int sofunc = nfunc_[soshell0][irrep];

            int naofunc = basis_->shell(aoshell0)->nfunction();
            if (naofunc_[soshell0] && (naofunc_[soshell0] != naofunc)) {
                throw PSIEXCEPTION("SOBasis::SOBasis: mismatch in naofunc");
            }
            naofunc_[soshell0] = naofunc;

            nfunc_[soshell0][irrep]++;
            nfuncall++;

            for (k=0; k<soblocks[i].so[j].length; k++) {
                int bfn = soblocks[i].so[j].cont[k].bfn;
                double coef = soblocks[i].so[j].cont[k].coef;
                int aoshell = basis_->function_to_shell(bfn);
                int aoshellfunc = bfn - basis_->shell_to_function(aoshell);
                int soshell = aoshell_to_soshell[aoshell];

                if (soshell != soshell0) {
                    throw PSIEXCEPTION("SOBasis::SOBasis: shell changed");
                }

                trans_[soshell].add_transform(aoshell,irrep, coef,aoshellfunc,sofunc);
            }
        }
    }

    if (nfuncall != basis_->nbf()) {
        throw PSIEXCEPTION("SOBasis::SOBasis: miscounted number of functions");
    }

    delete[] soblocks;
    delete[] aoshell_to_soshell;

    for (i=0; i<nshell_; i++) {
        funcoff_[i][0] = 0;
        for (j=1; j<nirrep_; j++) {
            funcoff_[i][j] = funcoff_[i][j-1] + nfunc_[i][j-1];
        }
    }

    func_ = new int[nshell_];
    irrep_ = new int[basis_->nbf()];
    func_within_irrep_ = new int[basis_->nbf()];
    nfunc_in_irrep_ = new int[nirrep_];

    for (i=0; i<nirrep_; i++) nfunc_in_irrep_[i] = 0;

    if (nshell_) {
        func_[0] = 0;
        for (i=1; i<nshell_; i++) {
            func_[i] = func_[i-1] + nfunction(i-1);
        }
        int ibasis_ = 0;
        for (i=0; i<nshell_; i++) {
            for (j=0; j<nirrep_; j++) {
                for (k=0; k<nfunc_[i][j]; k++,ibasis_++) {
                    irrep_[ibasis_] = j;
                    func_within_irrep_[ibasis_] = nfunc_in_irrep_[j]++;
                }
            }
        }
    }
}

SOBasis::~SOBasis()
{
    for (int i=0; i<nshell_; i++) {
        delete[] nfunc_[i];
        delete[] funcoff_[i];
    }
    delete[] nfunc_;
    delete[] funcoff_;
    delete[] naofunc_;
    delete[] ncomp_;
    delete[] trans_;
    delete[] func_;
    delete[] irrep_;
    delete[] func_within_irrep_;
    delete[] nfunc_in_irrep_;
}

int SOBasis::max_nfunction_in_shell() const
{
    int maxn = 0;
    for (int i=0; i<nshell_; i++) {
        int n = nfunction(i);
        if (n > maxn) maxn = n;
    }
    return maxn;
}

int SOBasis::nfunction(int ishell) const
{
    int n=0;
    for (int i=0; i<nirrep_; i++) {
        n += nfunc_[ishell][i];
    }
    return n;
}


