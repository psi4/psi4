#include "aobasis.h"
#include "index.h"

using namespace yeti;
using namespace std;

AtomBasis::AtomBasis(
    const std::string &symbol,
    const std::string &name
) : symbol_(symbol), name_(name), nbasis_(0), nshell_(0)
{
}

size_t
AtomBasis::nbasis() const
{
    return nbasis_;
}

size_t
AtomBasis::nshell() const
{
    return nshell_;
}

AtomBasis::iterator
AtomBasis::begin() const
{
    return shell_sizes_.begin();
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
    return shell_sizes_.end();
}

void
AtomBasis::add_shell(size_t size)
{
    shell_sizes_.push_back(size);
    ++nshell_;
    nbasis_ += size;
}


AOBasis::AOBasis(
    const std::string& id,
    const std::string& name
) : id_(id), name_(name), nbasis_(0)
{
}

size_t
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
AOBasis::get_index_range() const
{


    SubindexTuplePtr atom_tuple(new SubindexTuple(atoms_.size()));

    //loop atoms
    vector<AtomBasisPtr>::const_iterator it_atom(atoms_.begin());
    vector<AtomBasisPtr>::const_iterator stop_atom(atoms_.end());
    uli atom_start = 0;
    uli shell_start = 0;
    uli atom_idx = 0;
    for ( ; it_atom != stop_atom; ++it_atom, ++atom_idx)
    {
        AtomBasisPtr atom(*it_atom);
        SubindexTuplePtr shell_tuple(new SubindexTuple(atom->nshell()));

        //loop shells
        AtomBasis::iterator it_shell(atom->begin());
        AtomBasis::iterator stop_shell(atom->end());
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

        IndexRange* atom_range(
            new IndexRange(
                atom_start,
                shell_tuple
           )
        );
        atom_tuple->set(atom_idx, atom_range);
        atom_start += atom->nshell();
    }
    IndexRange* range(new IndexRange(0, atom_tuple));
    return range;
}


