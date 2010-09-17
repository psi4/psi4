#include <cstdio>

#include "mints.h"

#include <psi4-dec.h>

using namespace psi;

SOTransformComponent::SOTransformComponent()
{
    aofunc_ = 0;
    sofunc_ = 0;
    irrep_ = 0;
    coef_ = 0.0;
}

void SOTransformComponent::init(int aofunc, int sofunc, int sofuncirrep,
                                int irrep, double coef)
{
    aofunc_      = aofunc;
    sofunc_      = sofunc;
    sofuncirrep_ = sofuncirrep;
    irrep_       = irrep;
    coef_        = coef;
}

SOTransformShell::SOTransformShell()
{
    aoshell_ = 0;
}

void SOTransformShell::add_function(int irrep, double coef, int aofunc,
                                    int sofunc, int sofuncirrep)
{
    SOTransformComponent newfunc;
    newfunc.init(aofunc, sofunc, sofuncirrep, irrep, coef);
    funcs_.push_back(newfunc);
}

void SOTransform::init(int nshells)
{
    for (int i=0; i<nshells; ++i) {
        SOTransformShell shell;
        shell.aoshell(i);
        aoshell_.push_back(shell);
    }
}

void SOTransform::add_transform(int aoshell, int irrep, int sofuncirrep,
                                double coef, int aofunc, int sofunc)
{
//    printf("aoshell %d, irrep %d, sofuncirrep %d, coef %14.10f, aofunc %d, sofunc %d\n",aoshell,irrep,sofuncirrep,coef,aofunc,sofunc);

    unsigned int i;

    // Find the aoshell
    for (i=0; i<aoshell_.size(); ++i) {
        if (aoshell_[i].aoshell() == aoshell)
            break;
    }

    aoshell_[i].add_function(irrep, coef, aofunc, sofunc, sofuncirrep);
    aoshell_[i].aoshell(aoshell);
}

///////////////////////////////////////////////////////////////////////////////

SOBasis::SOBasis(const boost::shared_ptr<BasisSet> &basis, const boost::shared_ptr<IntegralFactory> &integral)
    : basis_(basis)
{
    int i, j, k;

    Molecule& mol = *basis_->molecule().get();

    CharacterTable ct = mol.point_group()->char_table();
    nirrep_ = ct.nirrep();

    // Count the number of SO shells.
    nshell_ = 0;
    for (i=0; i<mol.nunique(); ++i) {
        nshell_ += basis_->nshell_on_center(mol.unique(i));
    }

    // map each ao shell to an so shell
    int *aoshell_to_soshell = new int[basis_->nshell()];
    int soshell = 0;
    for (i=0; i<mol.nunique(); ++i) {
        for (j=0; j<basis_->nshell_on_center(mol.unique(i)); ++j) {
            for (k=0; k<mol.nequivalent(i); ++k) {
                int aoshell = basis_->shell_on_center(mol.equivalent(i, j), j);
                aoshell_to_soshell[aoshell] = soshell;
            }
            soshell++;
        }
    }

    ncomp_ = new int[nirrep_];
    for (i=0; i<nirrep_; ++i) {
        ncomp_[i] = ct.gamma(i).degeneracy();
        if (ncomp_[i] != 1) {
            throw PSIEXCEPTION("WARNING: SOBasis not tested for degenerate point groups");
        }
    }

    naofunc_ = new int[nshell_];
    memset(naofunc_, 0, sizeof(int) *nshell_);

    nfunc_ = new int*[nshell_];
    funcoff_ = new int*[nshell_];
    for (i=0; i<nshell_; ++i) {
        nfunc_[i] = new int[nirrep_];
        funcoff_[i] = new int[nirrep_];
        for (j=0; j<nirrep_; ++j) {
            nfunc_[i][j] = 0;
        }
    }

    PetiteList petite(basis_, integral);

    int nblocks = petite.nblocks();
    SO_block *soblocks = petite.aotoso_info();

    trans_ = new SOTransform[nshell_];

    int nfuncall = 0;
    for (i=0; i<nblocks; ++i) {
        int irrep = ct.which_irrep(i);
        for (j=0; j<soblocks[i].len; ++j) {
            if (soblocks[i].so[j].length == 0) continue;

            int bfn0 = soblocks[i].so[j].cont[0].bfn;
            int aoshell0 = basis_->function_to_shell(bfn0);
            int soshell0 = aoshell_to_soshell[aoshell0];
            int sofunc = nfunc_[aoshell0][irrep];

            int naofunc = basis_->shell(aoshell0)->nfunction();
            if (naofunc_[soshell0] && (naofunc_[soshell0] != naofunc)) {
                throw PSIEXCEPTION("Error: SOBasis: mismatch in naofunc");
            }
            naofunc_[soshell0] = naofunc;

            nfunc_[soshell0][irrep]++;
            nfuncall++;

            for (k=0; k<soblocks[i].so[j].length; ++k) {
                int bfn = soblocks[i].so[j].cont[k].bfn;
                double coef = soblocks[i].so[j].cont[k].coef;
                int aoshell = basis_->function_to_shell(bfn);
                int aoshellfunc = bfn - basis_->shell_to_function(aoshell);
                int soshell = aoshell_to_soshell[aoshell];

                if (soshell != soshell0) {
                    throw PSIEXCEPTION("Error: SOBasis: shell changed");
                }

                // completely different function call. Function signature:
//                void SOTransform::add_transform(int aoshell, int irrep, int sofuncirrep,
//                                                double coef, int aofunc, int sofunc)
                //                                            should not be zero
                trans_[soshell].add_transform(aoshell, irrep, 0, coef, aoshellfunc, sofunc);
            }
        }
    }

    if (nfuncall != basis_->nbf()) {
        throw PSIEXCEPTION("Error: SOBasis: miscounted number of functions");
    }

    delete[] soblocks;
    delete[] aoshell_to_soshell;

    for (i=0; i<nshell_; ++i) {
        funcoff_[i][0] = 0;
        for (j=1; j<nirrep_; ++j) {
            funcoff_[i][j] = funcoff_[i][j-1] + nfunc_[i][j-1];
        }
    }

    func_ = new int[nshell_];
    irrep_ = new int[basis_->nbf()];
    func_within_irrep_ = new int[basis_->nbf()];
    nfunc_in_irrep_ = new int[nirrep_];

    for (i=0; i<nirrep_; ++i) nfunc_in_irrep_[i] = 0;

    if (nshell_) {
        func_[0] = 0;
        for (i=1; i<nshell_; ++i) {
            func_[i] = func_[i-1] + nfunction(i-1);
        }
        int ibasis_ = 0;
        for (i=0; i<nshell_; ++i) {
            for (j=0; j<nirrep_; ++j) {
                for (k=0; k<nfunc_[i][j]; ++k,++ibasis_) {
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


