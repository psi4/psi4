#include <cstdio>

#include <libmints/sobasis.h>

using namespace psi;

extern FILE *outfile;

SOTransformComponent::SOTransformComponent()
{
    aofunc_ = 0;
    sofunc_ = 0;
    irrep_ = 0;
    coef_ = 0.0;
}

void SOTransformComponent::init(int aofunc, int sofunc, int sofuncirrep, int irrep, double coef)
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

void SOTransformShell::add_function(int irrep, double coef, int aofunc, int sofunc, int sofuncirrep)
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

void SOTransform::add_transform(int aoshell, int irrep, int sofuncirrep, double coef, int aofunc, int sofunc)
{
    unsigned int i;
    
    // Find the aoshell
    for (i=0; i<aoshell_.size(); ++i) {
        if (aoshell_[i].aoshell() == aoshell)
            break;
    }

    aoshell_[i].add_function(irrep, coef, aofunc, sofunc, sofuncirrep);
    aoshell_[i].aoshell(aoshell);
}
