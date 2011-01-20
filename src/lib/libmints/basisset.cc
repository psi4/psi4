/*!
    \defgroup MINTS libmints: Integral library
    \ingroup MINTS
*/

#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libparallel/parallel.h>
#include <psifiles.h>

#include "vector3.h"
#include "molecule.h"
#include "basisset.h"
#include "dimension.h"
#include "sobasis.h"
#include "integral.h"
#include "symmetry.h"
#include "gshell.h"
#include "factory.h"
#include "basisset_parser.h"
#include "pointgrp.h"
#include "wavefunction.h"

using namespace std;
using namespace psi;
using namespace boost;

once_flag BasisSet::initialized_shared_ = BOOST_ONCE_INIT;

std::vector<Vector3> BasisSet::exp_ao[LIBINT_MAX_AM];

BasisSet::BasisSet()
{
    call_once(initialize_singletons, initialized_shared_);
}

BasisSet::~BasisSet()
{

}

void BasisSet::initialize_singletons()
{
    // Populate the exp_ao arrays
    for (int l=0; l<LIBINT_MAX_AM; ++l) {
        for (int i=0; i<=l; ++i) {
            int x = l-i;
            for (int j=0; j<=i; ++j) {
                int y = i-j;
                int z = j;

                Vector3 xyz_ao(x, y, z);
                BasisSet::exp_ao[l].push_back(xyz_ao);
            }
        }
    }
}

void BasisSet::print(FILE *out) const
{
    fprintf(out, "  Basis Set\n");
    fprintf(out, "    Number of shells: %d\n", nshell());
    fprintf(out, "    Number of basis function: %d\n", nbf());
    fprintf(out, "    Number of Cartesian functions: %d\n", nao());
    fprintf(out, "    Spherical Harmonics?: %s\n", has_puream() ? "true" : "false");
    fprintf(out, "    Max angular momentum: %d\n\n", max_am());

//    fprintf(out, "    Shells:\n\n");
//    for (int s=0; s<nshell(); ++s)
//        shells_[s]->print(out);
}

shared_ptr<GaussianShell> BasisSet::shell(int si) const
{
    if (si < 0 || si > nshell())
        throw PSIEXCEPTION("BasisSet::shell: requested shell is out-of-bounds.");
    return shells_[si];
}

shared_ptr<GaussianShell> BasisSet::shell(int center, int si) const
{
    return shell(center_to_shell_[center] + si);
}

shared_ptr<BasisSet> BasisSet::zero_ao_basis_set()
{
    shared_ptr<BasisSet> new_basis(new BasisSet());

    // Setup all the parameters needed for a zero basis set
//    new_basis->shell_first_basis_function_ = NULL;
//    new_basis->shell_first_ao_ = NULL;
//    new_basis->shell_center_ = new int[1];
//    new_basis->shell_center_[0] = 0;
    new_basis->shell_center_.push_back(0);
    new_basis->max_nprimitive_ = 1;
    new_basis->max_am_ = 0;

    new_basis->puream_ = false;

    // Add "ghost" atom to the molecule for this basis
    new_basis->molecule_ = shared_ptr<Molecule>(new Molecule);
    // Ghost atoms are now handled differently, they are not added to the normal xyz information array,
    // but to the fxyz array.
    new_basis->molecule_->add_atom(0, 0.0, 0.0, 0.0);
    Vector3 center = new_basis->molecule_->fxyz(0);

    new_basis->nshell_ = 1;
    new_basis->nprimitive_ = 1;
    new_basis->nao_ = 1;
    new_basis->nbf_ = 1;

    // Create shell array
    new_basis->shells_.push_back(shared_ptr<GaussianShell>(new GaussianShell));

    // Spherical and SO-transforms are expected even if not used.
//    new_basis->sphericaltransforms_.push_back(SphericalTransform(0));

    // Create out basis set arrays
    // null basis set
    int am   = 0;
    double e = 0.0;
    double c = 1.0;

    // Add the null-s-function
    new_basis->shells_[0]->init(1, &e, am, Cartesian, &c, 0, center, 0, Normalized);

    return new_basis;
}

shared_ptr<SOBasisSet> BasisSet::zero_so_basis_set(const shared_ptr<IntegralFactory>& factory)
{
    shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    shared_ptr<SOBasisSet> sozero(new SOBasisSet(zero, factory));
    return sozero;
}

shared_ptr<BasisSet> BasisSet::construct(const shared_ptr<BasisSetParser>& parser,
        const shared_ptr<Molecule>& mol,
        const std::string& basisname)
{
    // Construct vector with the same basis name for each element and pass to
    // the other version of construct
    vector<string> basisnames;

    // Update geometry in molecule, if there is a problem an exception is thrown
    mol->update_geometry();

    for (int i=0; i<mol->natom(); ++i) {
        basisnames.push_back(basisname);
    }

    return construct(parser, mol, basisnames);
}

shared_ptr<BasisSet> BasisSet::construct(const shared_ptr<BasisSetParser>& parser,
        const shared_ptr<Molecule>& mol,
        const std::vector<std::string>& basisnames)
{
    shared_ptr<BasisSet> basisset(new BasisSet);

    // Update geometry in molecule, if there is a problem an exception is thrown.
    mol->update_geometry();

    basisset->molecule_ = mol;
    parser->parse(basisset, basisnames);

    return basisset;
}

void BasisSet::refresh()
{
    // Reset data to initial values
    nshell_ = shells_.size();
    nprimitive_ = 0;
    nao_ = 0;
    nbf_ = 0;
    max_am_ = 0;
    max_nprimitive_ = 0;
    puream_ = false;

    shell_first_basis_function_.clear(); shell_first_basis_function_.resize(nshell_, 0);
    shell_first_ao_.clear();             shell_first_ao_.resize(nshell_, 0);
    shell_center_.clear();               shell_center_.resize(nshell_, 0);
    center_to_nshell_.clear();           center_to_nshell_.resize(molecule_->natom(), 0);
    center_to_shell_.clear();            center_to_shell_.resize(molecule_->natom(), 0);
    center_to_shell_[0] = 0;

    int current_center = 0;

    for (int i=0; i<nshell(); ++i) {
        shell_center_[i]   = shells_[i]->ncenter();
        shell_first_ao_[i] = nao_;
        shell_first_basis_function_[i] = nbf_;
        shells_[i]->set_function_index(nbf_);

        center_to_nshell_[shell_center_[i]]++;
        if (current_center != shell_center_[i]) {
            center_to_shell_[shell_center_[i]] = i;
            current_center = shell_center_[i];
        }

        nprimitive_ += shells_[i]->nprimitive();
        nao_        += shells_[i]->ncartesian();
        nbf_        += shells_[i]->nfunction();

        if (max_am_ < shells_[i]->am())
            max_am_ = shells_[i]->am();

        if (max_nprimitive_ < shells_[i]->nprimitive())
            max_nprimitive_ = shells_[i]->nprimitive();

        if (puream_ == false && shells_[i]->is_pure())
            puream_ = true;
    }

    function_to_shell_.resize(nbf());
    int ifunc = 0;
    for (int i=0; i<nshell_; ++i) {
        int nfun = shell(i)->nfunction();
        for (int j=0; j<nfun; ++j) {
            function_to_shell_[ifunc] = i;
            ifunc++;
        }
    }
}

shared_ptr<BasisSet> BasisSet::atomic_basis_set(int center)
{
    //May only work in C1!!!!
    //Construct a blank BasisSet on the heap
    shared_ptr<BasisSet> bas(new BasisSet);
    //Construct a blank Molecule on the heap
    shared_ptr<Molecule> mol(new Molecule);

    int Z = molecule_->Z(center);
    double x = molecule_->x(center);
    double y = molecule_->y(center);
    double z = molecule_->z(center);
    double mass = molecule_->fmass(center);
    double charge = molecule_->fcharge(center);
    std::string lab = molecule_->flabel(center);
    char* label = new char[lab.length() + 1];
    strcpy(label,lab.c_str());

    //Put the atomic info into mol
    mol->add_atom(Z, 0.0, 0.0, 0.0, label, mass, charge);
    Vector3 v(0.0,0.0,0.0);

    //Assign the atomic molecule to bas
    bas->molecule_ = mol;

    //Go through shells in current basis set
    //Push shells that belong to center
    //to bas
    int current_shells = 0;
    int shell_start = -1;
    int ao_start = 0;
    int so_start = 0;
    for (int i = 0; i<nshell(); i++) {
        if (shell_center_[i] < center) {
            ao_start += shells_[i]->ncartesian();
            so_start += shells_[i]->nfunction();
        }
        if (shell_center_[i] == center) {
            if (shell_start == -1)
                shell_start = i;
            shared_ptr<GaussianShell> shell(new GaussianShell);
            int nprm = shells_[i]->nprimitive();
            int am = shells_[i]->am();
            GaussianType harmonics = (shells_[i]->is_pure() ? Pure : Cartesian);
            int nc = 0; // In the atomic basis, always on the 0th atom
            int start = 0; //Will be reset later
            double* e = shells_[i]->exps();
            double* c = shells_[i]->coefs();

            shell->init(nprm,e,am, harmonics, c,nc,
              v,start);

            bas->shells_.push_back(shell);
            current_shells++;
        }
        if (shell_center_[i] > center)
            break;
    }

    //Populate SOTransform
#warning Sorry Rob, I probably broke SAD again.
#if 0
    for (int i = 0; i< current_shells; i++) {
        //OK, this is the SOTransformShell in the full basis
        SOTransformShell* full_so = sotransform_->aoshell(shell_start+i);
        //and it must go into the SOTransform in the atomic basis, and the offset in shell, ao, and so must be respected
        for (int ao = 0; ao < full_so->nfunc; ao++) {
            //printf("AO %d, AOStart %d, SOStart%d, i %d, shell_start %d.\n",ao,ao_start,so_start,i,shell_start);
            bas->sotransform_->add_transform(i,full_so->func[ao].irrep(),full_so->func[ao].sofuncirrep()-so_start,full_so->func[ao].coef(), \
                full_so->func[ao].aofunc(),full_so->func[ao].sofunc()-so_start);
        }
        //ao_start += shells_[i+shell_start]->ncartesian();
        //so_start += shells_[i+shell_start]->nfunction();
    }
#endif

    //Setup the indexing in the atomic basis
    bas->refresh();

    //And ... return
    return bas;
}

