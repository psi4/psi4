#include "aobasis.h"
#include "index.h"

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define MAX_AM 10

static ShellBasis* null_shell_basis = 0;

ShellBasis::ShellBasis(usi ang_mom)
    : ang_mom_(ang_mom)
{

}

void
ShellBasis::add_shell(uli size)
{
    shell_sizes_.push_back(size);
}

uli
ShellBasis::nshell() const
{
    return shell_sizes_.size();
}

ShellBasis::iterator
ShellBasis::begin() const
{
    return shell_sizes_.begin();
}

ShellBasis::iterator
ShellBasis::end() const
{
    return shell_sizes_.end();
}

AtomBasis::AtomBasis(
    const std::string &symbol,
    const std::string &name,
    double x, double y, double z
) :
    symbol_(symbol),
    name_(name),
    nbasis_(0),
    nshell_(0),
    x_(x), y_(y), z_(z),
    shells_(MAX_AM, null_shell_basis)
{
}

uli
AtomBasis::nbasis() const
{
    return nbasis_;
}

uli
AtomBasis::nshell() const
{
    return nshell_;
}

AtomBasis::iterator
AtomBasis::begin() const
{
    return shells_.begin();
}

const std::string&
AtomBasis::name() const
{
    return name_;
}

const std::string&
AtomBasis::symbol() const
{
    return symbol_;
}

AtomBasis::iterator
AtomBasis::end() const
{
    return shells_.end();
}

void
AtomBasis::add_shell(uli size, usi ang_mom)
{
    //ShellBasisPtr& shell = shells_[ang_mom];
    //if (!shell)
    ShellBasisPtr shell = new ShellBasis(ang_mom);
    shell->add_shell(size);
    shells_.push_back(shell);
    ++nshell_;
    nbasis_ += size;
}


AOBasis::AOBasis(
    const std::string& id,
    const std::string& name
) : id_(id), name_(name), nbasis_(0)
{
}

uli
AOBasis::nbasis() const
{
    return nbasis_;
}

void
AOBasis::add_atom(const AtomBasisPtr &atom)
{
    nbasis_ += atom->nbasis();
    atoms_.push_back(atom);
}

IndexRange*
AOBasis::get_index_range(uli nranges) const
{
    //create a vector of atom shells
    vector<AtomShellPtr> atom_shells;




    //loop atoms
    vector<AtomBasisPtr>::const_iterator it_atom(atoms_.begin());
    vector<AtomBasisPtr>::const_iterator stop_atom(atoms_.end());
    for ( ; it_atom != stop_atom; ++it_atom)
    {
        AtomBasisPtr atom(*it_atom);
        AtomBasis::iterator it_shell(atom->begin());
        AtomBasis::iterator stop_shell(atom->end());
        for ( ; it_shell != stop_shell; ++it_shell)
        {
            ShellBasisPtr shell(*it_shell);
            if (!shell)
                continue;

            AtomShellPtr atom_shell = new AtomShell(atom, shell);
            atom_shells.push_back(atom_shell);
        }
    }


    uli atom_shell_start = 0;
    uli shell_start = 0;
    uli atom_shell_idx = 0;
    SubindexTuplePtr atom_shell_tuple(new SubindexTuple(atom_shells.size()));
    vector<AtomShellPtr>::const_iterator it(atom_shells.begin());
    vector<AtomShellPtr>::const_iterator stop(atom_shells.end());
    for ( ; it != stop; ++it)
    {
        AtomShellPtr atom_shell = *it;

        SubindexTuplePtr shell_tuple(new SubindexTuple(atom_shell->nshell()));

        //loop shells
        AtomShell::iterator it_shell(atom_shell->begin());
        AtomShell::iterator stop_shell(atom_shell->end());
        uli shell_idx = 0;
        for (uli nfxn_shell = 0; it_shell != stop_shell; ++it_shell, ++shell_idx)
        {
            nfxn_shell = *it_shell;
            IndexRange* shell_range(
                new IndexRange(
                    shell_start,
                    nfxn_shell
                )
            );
            shell_tuple->set(shell_idx, shell_range);
            shell_start += nfxn_shell;
        }

        IndexRange* atom_shell_range(
            new IndexRange(
                atom_shell_start,
                shell_tuple
           )
        );
        atom_shell_tuple->set(atom_shell_idx, atom_shell_range);
        atom_shell_start += atom_shell->nshell();
        ++atom_shell_idx;
    }

    uli ndesired = nbasis_ / nranges;

    uli ntuple = atom_shell_tuple->size();
    uli nfinal_tuple = 0;
    uli sizes[1000];
    uli nranges_on_tuple = 0;
    uli nfxns = 0;
    for (uli i=0; i < ntuple; ++i)
    {
        uli ntot = atom_shell_tuple->get(i)->ntot();
        ++nranges_on_tuple;
        nfxns += ntot;
        if (nfxns >= ndesired)
        {
            sizes[nfinal_tuple] = nranges_on_tuple;
            ++nfinal_tuple;
            nranges_on_tuple = 0;
            nfxns = 0;
        }
    }

    if (nranges_on_tuple)
    {
        sizes[nfinal_tuple] = nranges_on_tuple;
        ++nfinal_tuple;
    }

    SubindexTuplePtr final_tuple = new SubindexTuple(nfinal_tuple);

    uli itot = 0;
    uli start = 0;
    for (uli i=0; i < nfinal_tuple; ++i)
    {
        uli tuple_size = sizes[i];
        SubindexTuplePtr subtuple = new SubindexTuple(tuple_size);
        for (uli isub = 0; isub < tuple_size; ++isub, ++itot)
        {
            subtuple->set(isub, atom_shell_tuple->get(itot)->get_subranges()->get(0));
        }
        IndexRange* subrange = new IndexRange(start, subtuple);
        start += tuple_size;
        final_tuple->set(i, subrange);
    }

    IndexRange* range(new IndexRange(0, final_tuple));

    return range;
}

AtomShell::AtomShell(
    const AtomBasisPtr &atom_basis,
    const ShellBasisPtr &shell_basis
)
    : shell_basis_(shell_basis),
    atom_basis_(atom_basis)
{


}

AtomShell::iterator
AtomShell::begin() const
{
    return shell_basis_->begin();
}

AtomShell::iterator
AtomShell::end() const
{
    return shell_basis_->end();
}

uli
AtomShell::nshell() const
{
    return shell_basis_->nshell();
}


