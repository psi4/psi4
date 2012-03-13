/*!
    \defgroup MINTS libmints: Integral library
    \ingroup MINTS
*/
#include <boost/shared_ptr.hpp>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>

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
#include "coordentry.h"

using namespace std;
using namespace psi;
using namespace boost;

once_flag BasisSet::initialized_shared_ = BOOST_ONCE_INIT;

std::vector<Vector3> BasisSet::exp_ao[LIBINT_MAX_AM];

namespace {
bool has_ending (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}
}

BasisSet::BasisSet()
{
    call_once(initialize_singletons, initialized_shared_);
}

boost::shared_ptr<BasisSet> BasisSet::build(boost::shared_ptr<Molecule> molecule,
                                            const std::vector<GaussianShell>& shells)
{
    boost::shared_ptr<BasisSet> basis(new BasisSet());
    basis->molecule_ = molecule;
    basis->shells_ = shells;
    basis->refresh();

    return basis;
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

boost::shared_ptr<Molecule> BasisSet::molecule() const
{
    return molecule_;
}

void BasisSet::print(FILE *out) const
{
    if (Communicator::world->me() == 0) {
        fprintf(out, "  Basis Set: %s\n", name_.c_str());
        fprintf(out, "    Number of shells: %d\n", nshell());
        fprintf(out, "    Number of basis function: %d\n", nbf());
        fprintf(out, "    Number of Cartesian functions: %d\n", nao());
        fprintf(out, "    Spherical Harmonics?: %s\n", has_puream() ? "true" : "false");
        fprintf(out, "    Max angular momentum: %d\n\n", max_am());
    }
}

void BasisSet::print_by_level(FILE* out, int level) const
{
    if (level < 1)
        return;
    else if (level == 1)
        print(out);
    else if (level == 2)
        print_summary(out);
    else if (level > 2)
        print_detail(out);
}

void BasisSet::print_summary(FILE* out) const
{
    if (Communicator::world->me() == 0) {
        fprintf(out, "  -AO BASIS SET INFORMATION:\n");
        fprintf(out, "    Name                   = %s\n", name_.c_str());
        fprintf(out, "    Total number of shells = %d\n", nshell());
        fprintf(out, "    Number of primitives   = %d\n", nprimitive_);
        fprintf(out, "    Number of AO           = %d\n", nao_);
        fprintf(out, "    Number of SO           = %d\n", nbf_);
        fprintf(out, "    Maximum AM             = %d\n", max_am_);
        fprintf(out, "    Spherical Harmonics    = %s\n", (puream_ ? "TRUE" : "FALSE"));
        fprintf(out, "\n");

        fprintf(out,"  -Contraction Scheme:\n");
        fprintf(out, "    Atom   Type   All Primitives // Shells:\n");
        fprintf(out, "   ------ ------ --------------------------\n");
    }

    int *nprims = new int[max_am_ + 1];
    int *nunique = new int[max_am_ + 1];
    int *nshells = new int[max_am_ + 1];
    char *amtypes = new char[max_am_ + 1];

    for (int A = 0; A < molecule_->natom(); A++) {

        memset((void*) nprims , '\0', (max_am_ + 1) * sizeof(int));
        memset((void*) nunique, '\0', (max_am_ + 1) * sizeof(int));
        memset((void*) nshells, '\0', (max_am_ + 1) * sizeof(int));

        if (Communicator::world->me() == 0) {
            fprintf(out, "    %4d    ", A+1);
            fprintf(out, "%2s     ", molecule_->symbol(A).c_str());
        }

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++) {
            const GaussianShell& shell = shells_[Q + first_shell];
            nshells[shell.am()]++;
            nunique[shell.am()]+= shell.nprimitive();
            nprims [shell.am()]+= shell.nprimitive();
            amtypes[shell.am()] = shell.amchar();
        }

        // All Primitives
        for (int l = 0; l < max_am_ + 1; l++) {
            if (nprims[l] == 0)
                continue;
            if (Communicator::world->me() == 0)
                fprintf(out, "%d%c ", nprims[l], amtypes[l]);
        }
        // Shells
        if (Communicator::world->me() == 0)
            fprintf(out, "// ");
        for (int l = 0; l < max_am_ + 1; l++) {
            if (nshells[l] == 0)
                continue;
            if (Communicator::world->me() == 0)
                fprintf(out, "%d%c ", nshells[l], amtypes[l]);
        }
        if (Communicator::world->me() == 0)
            fprintf(out, "\n");
    }
    if (Communicator::world->me() == 0)
        fprintf(out, "\n");

    delete[] nprims;
    delete[] nunique;
    delete[] nshells;
    delete[] amtypes;
}

void BasisSet::print_detail(FILE* out) const
{
    print_summary(out);

    //TODO: Use unique atoms (C1 for now)
    for (int A = 0; A < molecule_->natom(); A++) {
        if (Communicator::world->me() == 0)
            fprintf(out, "  -Basis set on unique center %d: %s\n", A+1,molecule_->symbol(A).c_str());

        int first_shell = center_to_shell_[A];
        int n_shell = center_to_nshell_[A];

        for (int Q = 0; Q < n_shell; Q++) {
            const GaussianShell& shell = shells_[Q + first_shell];

            for (int K = 0; K < shell.nprimitive(); K++) {
                if (Communicator::world->me() == 0) {
                    if (K == 0)
                        fprintf(outfile, "     %c ", shell.AMCHAR());
                    else
                        fprintf(outfile, "       ");
                    fprintf(outfile, "(%20.8f %20.8f)\n",shell.exp(K), shell.coef(K));
                }

            }
        }
        if (Communicator::world->me() == 0)
            fprintf(out, "\n");
    }
}

const GaussianShell& BasisSet::shell(int si) const
{
    if (si < 0 || si > nshell())
        throw PSIEXCEPTION("BasisSet::shell: requested shell is out-of-bounds.");
    return shells_[si];
}

const GaussianShell& BasisSet::shell(int center, int si) const
{
    return shell(center_to_shell_[center] + si);
}

boost::shared_ptr<BasisSet> BasisSet::zero_ao_basis_set()
{
    boost::shared_ptr<BasisSet> new_basis(new BasisSet());

    // Setup all the parameters needed for a zero basis set
    new_basis->shell_center_.push_back(0);
    new_basis->max_nprimitive_ = 1;
    new_basis->max_am_ = 0;

    new_basis->puream_ = false;

    // Add "ghost" atom to the molecule for this basis
    new_basis->molecule_ = boost::shared_ptr<Molecule>(new Molecule);
    // Ghost atoms are now handled differently, they are not added to the normal xyz information array,
    // but to the fxyz array.
    new_basis->molecule_->add_atom(0, 0.0, 0.0, 0.0);
    Vector3 center = new_basis->molecule_->fxyz(0);

    new_basis->nprimitive_ = 1;
    new_basis->nao_ = 1;
    new_basis->nbf_ = 1;

    // Create shell array
    std::vector<double> e(1);
    std::vector<double> c(1);
    e[0] = 0.0;
    c[0] = 1.0;

    new_basis->shells_.push_back(GaussianShell(0, c, e, Cartesian, 0, center, 0));

    new_basis->refresh();

    return new_basis;
}

//boost::shared_ptr<SOBasisSet> BasisSet::zero_so_basis_set(const boost::shared_ptr<IntegralFactory>& factory)
//{
//    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
//    boost::shared_ptr<SOBasisSet> sozero(new SOBasisSet(zero, factory));
//    return sozero;
//}

boost::shared_ptr<BasisSet> BasisSet::construct(const boost::shared_ptr<BasisSetParser>& parser,
        const boost::shared_ptr<Molecule>& mol,
        const std::string& type)
{
    boost::shared_ptr<BasisSet> basisset(new BasisSet);

    // Update geometry in molecule, if there is a problem an exception is thrown.
    mol->update_geometry();

    // Set the default name to just the type
    basisset->set_name(type);

    // Assign the molecule to the basis set
    basisset->molecule_ = mol;

    // For each one try to load the basis set
    const list<string>& user_list = Process::environment.user_basis_files;

    // Map of GaussianShells
    //  basis           atom        gaussian shells
    typedef map<string, map<string, vector<GaussianShell> > > map_ssv;
    typedef map<string, vector<GaussianShell> > map_sv;
    map_ssv basis_atom_shell;

    map<string, int> names;

    for (int atom=0; atom<mol->natom(); ++atom) {

        const string& symbol = mol->atom_entry(atom)->symbol();
        const string& basisname = mol->atom_entry(atom)->basisset(type);

        if (basisname.empty())
            throw PSIEXCEPTION("BasisSet::construct: No basis set specified for "
                               +symbol+ " and " +type+" type.");

        names[basisname] = 1;

        // Add basisname, symbol to the list by clearing the vector.
        basis_atom_shell[basisname][symbol].clear();
    }

    BOOST_FOREACH(map_ssv::value_type& basis, basis_atom_shell)
    {
//fprintf(outfile, "Working on basis %s\n", basis.first.c_str());
         bool not_found = true;

        BOOST_FOREACH(string user_file, user_list)
        {
            boost::filesystem::path bf_path;
            bf_path = boost::filesystem::system_complete(user_file);
//fprintf(outfile, "Working on file %s\n", user_file.c_str());
            // Load in the basis set and remove it from atomsymbol_to_basisname
            vector<string> file = parser->load_file(bf_path.string(), basis.first);

            BOOST_FOREACH(map_sv::value_type& atom, basis.second) {
//fprintf(outfile, "Working on atom %s\n", atom.first.c_str());
                string symbol = atom.first;

                // Don't even look, if this has already been found
                if(!basis_atom_shell[basis.first][symbol].empty()) continue;

                try {
                    // Need to wrap this is a try catch block
                    basis_atom_shell[basis.first][symbol] = parser->parse(symbol, file);

                    if (Communicator::world->me() == 0)
                        fprintf(outfile, "  Basis set %s for %s read from %s\n",
                                basis.first.c_str(), symbol.c_str(), user_file.c_str());
                    not_found = false;
                }
                catch (BasisSetNotFound& e) {
                    // This is thrown when load_file fails
                    if (Communicator::world->me() == 0)
                        fprintf(outfile, "  Unable to find %s for %s in %s.\n",
                                basis.first.c_str(), symbol.c_str(), user_file.c_str());
                    not_found = true;
                }
            }
        }

        string filename = make_filename(basis.first);
        string path = Process::environment("PSIDATADIR");
        vector<string> file;

        try {
            // Don't even look, if this has already been found
            BOOST_FOREACH(map_sv::value_type& atom, basis.second) {
                string symbol = atom.first;
                // Don't bother looking if we've already found this
                if (atom.second.empty()){
                    if(file.empty()) file = parser->load_file(path + "/basis/" + filename);
                    // If not found this will throw...let it.
                    basis_atom_shell[basis.first][symbol] = parser->parse(symbol, file);
                }
            }
        }
        catch (BasisSetFileNotFound& e) {
            fprintf(outfile, " Unable to load %s from the default Psi4 basis set library.\n", filename.c_str());
            throw PSIEXCEPTION("  Unable to load "+ filename + " from the default Psi4 basis set library.");
        }
    }

    // Go through the atoms and copy the shells to the basis set.
    for (int atom=0; atom<mol->natom(); ++atom) {
        string basis = mol->atom_entry(atom)->basisset(type);
        string symbol = mol->atom_entry(atom)->symbol();
        Vector3 center = mol->xyz(atom);

        vector<GaussianShell>& shells = basis_atom_shell[basis][symbol];

        for (int i=0; i<shells.size(); ++i) {
            basisset->shells_.push_back(shells[i].copy(atom, center));
        }
    }

    // This step is very important. Without it the basis set is useless.
    basisset->refresh();

    basisset->name_.clear();
    for (map<string, int>::iterator iter = names.begin(), end = names.end();
         iter != end;
         ++iter) {
        basisset->name_ += (*iter).first + " + ";
    }

    if (has_ending(basisset->name_, " + "))
        basisset->name_ = basisset->name_.substr(0, basisset->name_.length()-3);

    return basisset;
}

std::string BasisSet::make_filename(const std::string& name)
{
    // Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
    string basisname = name;

    // First make it lower case
    transform(basisname.begin(), basisname.end(), basisname.begin(), ::tolower);

    string format_underscore("_"); // empty string
    // Replace all '(' with '_'
    xpressive::sregex match_format = xpressive::as_xpr("(");
    basisname = regex_replace(basisname, match_format, format_underscore);

    // Replace all ')' with '_'
    match_format = xpressive::as_xpr(")");
    basisname = regex_replace(basisname, match_format, format_underscore);

    // Replace all ',' with '_'
    match_format = xpressive::as_xpr(",");
    basisname = regex_replace(basisname, match_format, format_underscore);

    // Replace all '*' with 's'
    match_format = xpressive::as_xpr("*");
    string format_star("s");
    basisname = regex_replace(basisname, match_format, format_star);

    // Replace all '+' with 'p'
    match_format = xpressive::as_xpr("+");
    string format_plus("p");
    basisname = regex_replace(basisname, match_format, format_plus);

    // Add file extension
    basisname += ".gbs";

    return basisname;
}

void BasisSet::refresh()
{
    // Reset data to initial values
    nprimitive_ = 0;
    nao_ = 0;
    nbf_ = 0;
    max_am_ = 0;
    max_nprimitive_ = 0;
    puream_ = false;

    shell_first_basis_function_.clear(); shell_first_basis_function_.resize(nshell(), 0);
    shell_first_ao_.clear();             shell_first_ao_.resize(nshell(), 0);
    shell_center_.clear();               shell_center_.resize(nshell(), 0);
    center_to_nshell_.clear();           center_to_nshell_.resize(molecule_->natom(), 0);
    center_to_shell_.clear();            center_to_shell_.resize(molecule_->natom(), 0);
    center_to_shell_[0] = 0;

    int current_center = 0;

    for (int i=0; i<nshell(); ++i) {
        shell_center_[i]   = shells_[i].ncenter();
        shell_first_ao_[i] = nao_;
        shell_first_basis_function_[i] = nbf_;
        shells_[i].set_function_index(nbf_);

        center_to_nshell_[shell_center_[i]]++;
        if (current_center != shell_center_[i]) {
            center_to_shell_[shell_center_[i]] = i;
            current_center = shell_center_[i];
        }

        nprimitive_ += shells_[i].nprimitive();
        nao_        += shells_[i].ncartesian();
        nbf_        += shells_[i].nfunction();

        if (max_am_ < shells_[i].am())
            max_am_ = shells_[i].am();

        if (max_nprimitive_ < shells_[i].nprimitive())
            max_nprimitive_ = shells_[i].nprimitive();

        if (puream_ == false && shells_[i].is_pure())
            puream_ = true;
    }

    function_to_shell_.resize(nbf());
    int ifunc = 0;
    for (int i=0; i<nshell(); ++i) {
        int nfun = shells_[i].nfunction();
        for (int j=0; j<nfun; ++j) {
            function_to_shell_[ifunc] = i;
            ifunc++;
        }
    }
    ao_to_shell_.resize(nao());
    ifunc = 0;
    for (int i=0; i<nshell(); ++i) {
        int nfun = shells_[i].ncartesian();
        for (int j=0; j<nfun; ++j) {
            ao_to_shell_[ifunc] = i;
            ifunc++;
        }
    }

    // Create a map that has a key/value pair
    // The key is the angular momentum function of the shell arranged in decending order
    // The value is the actual shell number
    typedef std::pair<int, int> am_to_shell_pair;
    std::multimap< int, int, std::less<int> > am_to_shell_list;
    for (int i=0; i < shells_.size(); i++) {
        am_to_shell_list.insert(am_to_shell_pair(shells_[i].nfunction(), i));
    }
    // This puts the sorted shell values into the sorted_shell_list_ vector
    // This can be used by the integral iterator to look up the value of the sorted shells
    std::multimap< int, int, std::less<int> >::iterator it;
    sorted_ao_shell_list_.clear();
    for (it=am_to_shell_list.begin(); it != am_to_shell_list.end(); it++) {
        //std::cout << "sorted shell size = " << it->first <<
        //        "\t, which belongs to shell number " << it->second << std::endl;
        sorted_ao_shell_list_.push_back(it->second);
    }
}

std::pair<std::vector<std::string>, boost::shared_ptr<BasisSet> > BasisSet::test_basis_set(int max_am)
{
    int max_centers = 4;
    int max_primitives = 10;
    int max_shells;

    std::vector<int> nprim;
    nprim.push_back(10);
    nprim.push_back(1);
    nprim.push_back(6);
    nprim.push_back(1);
    nprim.push_back(2);
    nprim.push_back(1);
    nprim.push_back(1);
    nprim.push_back(1);
    nprim.push_back(1);
    nprim.push_back(1);

    std::vector<int> am;
    am.push_back(0);
    am.push_back(0);
    am.push_back(1);
    am.push_back(1);
    am.push_back(2);
    am.push_back(2);
    am.push_back(3);
    am.push_back(4);
    am.push_back(5);
    am.push_back(2);

    std::vector<std::vector<double> > c;
    c.push_back(std::vector<double>(10));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(6));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(2));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));
    c.push_back(std::vector<double>(1));

    std::vector<std::vector<double> > e;
    e.push_back(std::vector<double>(10));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(6));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(2));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));
    e.push_back(std::vector<double>(1));

    c[0][0] = 0.458878E-03;
    c[0][1] = 0.355070E-02;
    c[0][2] = 0.182618E-01;
    c[0][3] = 0.716650E-01;
    c[0][4] = 0.212346E+00;
    c[0][5] = 0.416203E+00;
    c[0][6] = 0.373020E+00;
    c[0][7] = 0.625054E-01;
    c[0][8] = 0.624532E-02;
    c[0][9] = 0.243374E-02;
    c[1][0] =     1.0;
    c[2][0] = 0.458878E-03;
    c[2][1] = 0.355070E-02;
    c[2][2] = 0.182618E-01;
    c[2][3] = 0.716650E-01;
    c[2][4] = 0.212346E+00;
    c[2][5] = 0.416203E+00;
    c[3][0] =     1.0;
    c[4][0] = 0.458878E-03;
    c[4][1] = 0.355070E-02;
    c[5][0] =     1.0;
    c[6][0] =     1.0;
    c[7][0] =     1.0;
    c[8][0] =     1.0;
    c[9][0] =     1.0;

    e[0][0] = 31700.0;
    e[0][1] =  4755.0;
    e[0][2] =  1082.0;
    e[0][3] =   306.0;
    e[0][4] =    99.0;
    e[0][5] =    33.0;
    e[0][6] =    13.0;
    e[0][7] =     4.0;
    e[0][8] =     2.0;
    e[0][9] =     0.5;
    e[1][0] =     1.0;
    e[2][0] = 31700.0;
    e[2][1] =  4755.0;
    e[2][2] =  1082.0;
    e[2][3] =   306.0;
    e[2][4] =    99.0;
    e[2][5] =    33.0;
    e[3][0] =     1.0;
    e[4][0] = 31700.0;
    e[4][1] =  4755.0;
    e[5][0] =     1.0;
    e[6][0] =     1.0;
    e[7][0] =     1.0;
    e[8][0] =     1.0;
    e[9][0] =     1.0;

    std::vector<std::string> labels;
    if (max_am > -1) {
        labels.push_back("S");
        labels.push_back("s");
        max_shells = 2;
    }
        labels.push_back("P");
    if (max_am > 0) {
        labels.push_back("p");
        max_shells = 4;
    }
    if (max_am > 1) {
        labels.push_back("D");
        labels.push_back("d");
        max_shells = 6;
    }
    if (max_am > 2) {
        labels.push_back("f");
        max_shells = 7;
    }
    if (max_am > 3) {
        labels.push_back("g");
        max_shells = 8;
    }
    if (max_am > 4) {
        labels.push_back("h");
        max_shells = 9;
    }
    if (max_am > 5) {
        labels.push_back("i");
        max_shells = 10;
    }

    boost::shared_ptr<BasisSet> new_basis(new BasisSet());

    // Add 4 atoms to the molecule for this basis (max integal centers is 4 at the moment)
    new_basis->molecule_ = boost::shared_ptr<Molecule>(new Molecule);
    // Ghost atoms are now handled differently, they are not added to the normal xyz information array,
    // but to the fxyz array.
    double x = 0.0;
    for (int A = 0; A < max_centers; A++) {
        new_basis->molecule_->add_atom(0, x, x, x);
        x += 1.0;
    }

    // Setup all the parameters needed for a zero basis set
    new_basis->shell_center_.resize(max_shells * max_centers);
    for (int A = 0; A < max_centers; A++)
        for (int Q = 0; Q < max_shells; Q++)
            new_basis->shell_center_[A*max_shells + Q] = A;
    new_basis->max_nprimitive_ = max_primitives;
    new_basis->max_am_ = max_am;

    // We'll time puream for now
    new_basis->puream_ = true;

    // Add shells
    for (int A = 0; A < max_centers; A++) {
        Vector3 center = new_basis->molecule_->fxyz(A);
        for (int Q = 0; Q < max_shells; Q++) {
            new_basis->shells_.push_back(GaussianShell(am[Q], c[Q], e[Q], Pure, A, center, 0));
        }
    }

    new_basis->refresh();

    return make_pair(labels, new_basis);
}

boost::shared_ptr<BasisSet> BasisSet::atomic_basis_set(int center)
{
    //May only work in C1!!!!
    //Construct a blank BasisSet on the heap
    boost::shared_ptr<BasisSet> bas(new BasisSet);
    //Construct a blank Molecule on the heap
    boost::shared_ptr<Molecule> mol(new Molecule);

    int Z = molecule_->Z(center);
    double x = molecule_->x(center);
    double y = molecule_->y(center);
    double z = molecule_->z(center);
    double mass = molecule_->mass(center);
    double charge = molecule_->charge(center);
    std::string lab = molecule_->label(center);
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
            ao_start += shells_[i].ncartesian();
            so_start += shells_[i].nfunction();
        }
        if (shell_center_[i] == center) {
            if (shell_start == -1)
                shell_start = i;
            int nprm = shells_[i].nprimitive();
            int am = shells_[i].am();
            GaussianType harmonics = (shells_[i].is_pure() ? Pure : Cartesian);
            int nc = 0; // In the atomic basis, always on the 0th atom
            int start = 0; //Will be reset later
            const std::vector<double>& e = shells_[i].exps();
            const std::vector<double>& c = shells_[i].coefs();

            bas->shells_.push_back(GaussianShell(am, c, e, harmonics, nc, v, start));
            current_shells++;
        }
        if (shell_center_[i] > center)
            break;
    }

    //Setup the indexing in the atomic basis
    bas->refresh();

    //And ... return
    return bas;
}
