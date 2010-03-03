/*!
    \defgroup MINTS libmints: Integral library
    \ingroup MINTS
*/

#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <psifiles.h>

#include "basisset.h"
#include "integral.h"
#include "symmetry.h"
#include "factory.h"
#include "vector3.h"
#include "basisset_parser.h"

using namespace psi;
using namespace boost;

once_flag BasisSet::initialized_shared_ = BOOST_ONCE_INIT;

std::vector<Vector3> BasisSet::exp_ao[LIBINT_MAX_AM];

BasisSet::BasisSet() :
    max_nprimitives_(0), shell_first_basis_function_(NULL), shell_first_ao_(NULL), shell_center_(NULL),
    max_stability_index_(0), uso2ao_(NULL), uso2bf_(NULL)
{
    call_once(initialize_singletons, initialized_shared_);
}

BasisSet::BasisSet(shared_ptr<Chkpt> chkpt, std::string basiskey) :
    max_nprimitives_(0), shell_first_basis_function_(NULL), shell_first_ao_(NULL), shell_center_(NULL),
    max_stability_index_(0), uso2ao_(NULL), uso2bf_(NULL),
    molecule_(new Molecule)
{
    call_once(initialize_singletons, initialized_shared_);

    // This requirement holds no matter what.
    puream_ = chkpt->rd_puream(basiskey.c_str()) ? true : false;

    // Initialize molecule, retrieves number of center and geometry
    molecule_->init_with_chkpt(chkpt);

    // Obtain symmetry information from the molecule.
    shared_ptr<PointGroup> pg = molecule_->find_point_group();
    molecule_->set_point_group(pg);
    molecule_->form_symmetry_information();

    // Initialize the shells
    initialize_shells(chkpt, basiskey);
}

BasisSet::~BasisSet()
{
    if (shell_first_basis_function_)
        Chkpt::free(shell_first_basis_function_);
    if (shell_first_ao_)
        Chkpt::free(shell_first_ao_);
    if (shell_center_)
        Chkpt::free(shell_center_);
    if (uso2ao_)
        Chkpt::free(uso2ao_);
    if (uso2bf_)
        Chkpt::free(uso2bf_);
}

void BasisSet::initialize_singletons()
{
//    fprintf(outfile, "In BasisSet::initialize_singletons:\n");
    // Populate the exp_ao arrays
    int ao;
    for (int l=0; l<LIBINT_MAX_AM; ++l) {
//    fprintf(outfile, "  l = %d\n", l);
    ao = 0;
    for (int i=0; i<=l; ++i) {
        int x = l-i;
        for (int j=0; j<=i; ++j) {
        int y = i-j;
        int z = j;

        Vector3 xyz_ao(x, y, z);
        BasisSet::exp_ao[l].push_back(xyz_ao);

//        fprintf(outfile, "    ao = %d %s\n", ao, xyz_ao.to_string().c_str());
        ao++;
        }
    }
    }
}

void BasisSet::initialize_shells(shared_ptr<Chkpt> chkpt, std::string& basiskey)
{
    // Initialize some data from checkpoint.
    nshells_      = chkpt->rd_nshell(basiskey.c_str());
    nprimitives_  = chkpt->rd_nprim(basiskey.c_str());
    nao_          = chkpt->rd_nao(basiskey.c_str());

    // Psi3 only allows either all Cartesian or all Spherical harmonic
    nbf_          = chkpt->rd_nso(basiskey.c_str());
    max_am_       = chkpt->rd_max_am(basiskey.c_str());
    uso2ao_       = chkpt->rd_usotao(basiskey.c_str());
    uso2bf_       = chkpt->rd_usotbf(basiskey.c_str());

    simple_mat_uso2ao_ = shared_ptr<SimpleMatrix>(new SimpleMatrix("Unique SO to AO transformation matrix", nbf_, nao_));
    simple_mat_uso2ao_->set(uso2ao_);
    // simple_mat_uso2ao_.print();

    simple_mat_uso2bf_ = shared_ptr<SimpleMatrix>(new SimpleMatrix("Unique SO to BF transformation matrix", nbf_, nbf_));
    simple_mat_uso2bf_->set(uso2bf_);
    // simple_mat_uso2bf_.print();

    // Retrieve angular momentum of each shell (1=s, 2=p, ...)
    int *shell_am = chkpt->rd_stype(basiskey.c_str());

    // Retrieve number of primitives per shell
    int *shell_num_prims = chkpt->rd_snumg(basiskey.c_str());

    // Retrieve exponents of primitive Gaussians
    double *exponents = chkpt->rd_exps(basiskey.c_str());

    // Retrieve coefficients of primitive Gaussian
    double **ccoeffs = chkpt->rd_contr_full(basiskey.c_str());

    // Retrieve pointer to first primitive in shell
    int *shell_fprim = chkpt->rd_sprim(basiskey.c_str());

    // Retrieve pointer to first basis function in shell
    shell_first_basis_function_ = chkpt->rd_sloc_new(basiskey.c_str());

    // Retrieve pointer to first AO in shell
    shell_first_ao_ = chkpt->rd_sloc(basiskey.c_str());

    // Retrieve location of shells (which atom it's centered on)
    shell_center_ = chkpt->rd_snuc(basiskey.c_str());

    // Initialize SphericalTransform
    for (int i=0; i<=max_am_; ++i) {
        sphericaltransforms_.push_back(SphericalTransform(i));
    }

    // Initialize SOTransform
    sotransform_ = shared_ptr<SOTransform>(new SOTransform);
    sotransform_->init(nshells_);

    int *so2symblk = new int[nbf_];
    int *so2index  = new int[nbf_];
    int *sopi = chkpt->rd_sopi(basiskey.c_str());
    int nirreps = chkpt->rd_nirreps();

    // Create so2symblk and so2index
    int ij = 0; int offset = 0;
    for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<sopi[h]; ++i) {
            so2symblk[ij] = h;
            so2index[ij] = ij-offset;

            ij++;
        }
        offset += sopi[h];
    }

    // Currently all basis sets are treated as segmented contractions
    // even though GaussianShell is generalized (well not really).
    int ao_start = 0;
    int puream_start = 0;

    for (int i=0; i<nshells_; ++i) {
        int am;
        am = shell_am[i] - 1;
        int fprim = shell_fprim[i] - 1;
        int nprims = shell_num_prims[i];
        Vector3 center = molecule_->xyz(shell_center_[i] - 1);
        double *cc = new double[nprims];
        for (int p=0; p<nprims; ++p) {
            cc[p] = ccoeffs[fprim+p][am];
        }

        // Construct a new shell. GaussianShell copies the data to new memory
        shells_.push_back(shared_ptr<GaussianShell>(new GaussianShell));
        shells_[i]->init(nprims, &(exponents[fprim]), am,
            puream_ ? GaussianShell::Pure : GaussianShell::Cartesian, cc, shell_center_[i]-1, center,
            puream_start);

        if (nprims > max_nprimitives_)
            max_nprimitives_ = nprims;

        delete[] cc;

        // OK, for a given number of AO functions in a shell INT_NCART(am)
        // beginning at column ao_start go through all rows finding where this
        // AO function contributes to an SO.
        for (int ao = 0; ao < INT_NCART(am); ++ao) {
            int aooffset = ao_start + ao;
            for (int so = 0; so < nbf_; ++so) {
                if (fabs(uso2ao_[so][aooffset]) >= 1.0e-14)
                    sotransform_->add_transform(i, so2symblk[so], so2index[so], uso2ao_[so][aooffset], ao, so);
            }
        }

        // Shift the ao starting index over to the next shell
        ao_start += INT_NCART(am);
        puream_start += INT_NFUNC(puream_, am);
    }

    delete[] so2symblk;
    delete[] so2index;
    Chkpt::free(sopi);
    Chkpt::free(ccoeffs);
    Chkpt::free(exponents);
    Chkpt::free(shell_am);
    Chkpt::free(shell_num_prims);
    Chkpt::free(shell_fprim);
}

void BasisSet::print(FILE *out) const
{
    fprintf(out, "  Basis Set\n");
    fprintf(out, "    Number of shells: %d\n", nshell());
    fprintf(out, "    Number of basis function: %d\n", nbf());
    fprintf(out, "    Number of Cartesian functions: %d\n", nao());
    fprintf(out, "    Spherical Harmonics?: %s\n", has_puream() ? "true" : "false");
    fprintf(out, "    Max angular momentum: %d\n\n", max_am());

    fprintf(out, "    Shells:\n\n");
    for (int s=0; s<nshell(); ++s)
        shells_[s]->print(out);
}

shared_ptr<GaussianShell> BasisSet::shell(int si) const
{
    #ifdef DEBUG
    assert(si < nshell());
    #endif
    return shells_[si];
}

shared_ptr<BasisSet> BasisSet::zero_basis_set()
{
    shared_ptr<BasisSet> new_basis(new BasisSet());

    // Setup all the parameters needed for a zero basis set
    new_basis->shell_first_basis_function_ = NULL;
    new_basis->shell_first_ao_ = NULL;
    new_basis->shell_center_ = new int[1];
    new_basis->shell_center_[0] = 0;
    new_basis->uso2ao_ = NULL;
    new_basis->max_nprimitives_ = 1;
    new_basis->max_stability_index_ = 0;
    new_basis->max_am_ = 0;

    new_basis->puream_ = false;

    // Add "ghost" atom to the molecule for this basis
    new_basis->molecule_ = shared_ptr<Molecule>(new Molecule);
    // Ghost atoms are now handled differently, they are not added to the normal xyz information array,
    // but to the fxyz array.
    new_basis->molecule_->add_atom(0, 0.0, 0.0, 0.0);
    Vector3 center = new_basis->molecule_->fxyz(0);

    new_basis->nshells_ = 1;
    new_basis->nprimitives_ = 1;
    new_basis->nao_ = 1;
    new_basis->nbf_ = 1;
    new_basis->uso2ao_ = Chkpt::matrix<double>(1, 1);
    new_basis->uso2ao_[0][0] = 1.0;
    new_basis->uso2bf_ = Chkpt::matrix<double>(1, 1);
    new_basis->uso2bf_[0][0] = 1.0;

    // Create shell array
    new_basis->shells_.push_back(shared_ptr<GaussianShell>(new GaussianShell));

    // Spherical and SO-transforms are expected even if not used.
    new_basis->sphericaltransforms_.push_back(SphericalTransform(0));
    new_basis->sotransform_ = shared_ptr<SOTransform>(new SOTransform);
    new_basis->sotransform_->init(1);

    // Create out basis set arrays
    // null basis set
    int am   = 0;
    double e = 0.0;
    double c = 1.0;

    // Add the null-s-function
    new_basis->shells_[0]->init(1, &e, am, GaussianShell::Cartesian, &c, 0, center, 0, GaussianShell::Normalized);

    // Add s-function SO transform.
    new_basis->sotransform_->add_transform(0, 0, 0, 1.0, 0, 0);

    return new_basis;
}

shared_ptr<BasisSet> BasisSet::construct(const shared_ptr<BasisSetParser>& parser,
        const shared_ptr<Molecule>& mol,
        const std::string& basisname)
{
    // Construct vector with the same basis name for each element and pass to
    // the other version of construct
    vector<string> basisnames;

    for (int i=0; i<mol->natom(); ++i)
    basisnames.push_back(basisname);

    return construct(parser, mol, basisnames);
}

shared_ptr<BasisSet> BasisSet::construct(const shared_ptr<BasisSetParser>& parser,
        const shared_ptr<Molecule>& mol,
        const std::vector<std::string>& basisnames)
{
    shared_ptr<BasisSet> basisset(new BasisSet);
    basisset->molecule_ = mol;
    parser->parse(basisset, basisnames);

    return basisset;
}
