#include "aobasis.h"
#include "index.h"
#include "runtime.h"
#include "gigmatrix.h"

#include <algorithm>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define MAX_AM 10

static Shell* null_shell = 0;

bool Shell::statics_done_ = false;
std::map<std::string, std::map<usi, double> > Shell::diffuse_cutoffs_;

int myupper(int c)
{
    return std::toupper((unsigned char)c);
}

bool
_lt_shells_am_diffusetype_atomgroup_atom(
    const ShellPtr& sh1,
    const ShellPtr& sh2
)
{
    usi am1 = sh1->angular_momentum();
    if (am1 == 1)
        am1 = 0;
    usi am2 = sh2->angular_momentum();
    if (am2 == 1)
        am2 = 0;
    if (am1 != am2)
        return am1 < am2;

    Shell::diffuse_type_t diff1 = sh1->diffuse_type();
    Shell::diffuse_type_t diff2 = sh2->diffuse_type();
    if (diff1 != diff2)
        return diff1 < diff2;

    uli group1 = sh1->get_atom()->group_number();
    uli group2 = sh2->get_atom()->group_number();
    if (group1 != group2)
        return group1 < group2;

    uli at1 = sh1->get_atom()->number();
    uli at2 = sh2->get_atom()->number();
    if (at1 != at2)
        return at1 < at2;

    //exactly the same... so... false...
    return false;
}

Shell::Shell(
    const AtomPtr& atom,
    usi ang_mom,
    uli nsets,
    uli nfxn,
    double min_exponent,
    uli program_shell_number,
    uli program_fxn_start
) :  atom_(atom),
    ang_mom_(ang_mom),
    nsets_(nsets),
    nfxn_(nfxn),
    min_exponent_(min_exponent),
    program_shell_number_(program_shell_number),
    program_fxn_start_(program_fxn_start),
    diffuse_type_(Shell::valence)
{
    if (!statics_done_)
        init_statics();

    double diff_cutoff = diffuse_cutoffs_[atom->name()][ang_mom];
    if (min_exponent < diff_cutoff)
        diffuse_type_ = Shell::diffuse;
}

void
Shell::init_statics()
{
    diffuse_cutoffs_["H"][0] = 0.05;
    diffuse_cutoffs_["H"][1] = 0.15;

    diffuse_cutoffs_["He"][0] = 0.08;
    diffuse_cutoffs_["He"][1] = 0.28;

    diffuse_cutoffs_["C"][0] = 0.05;
    diffuse_cutoffs_["C"][1] = 0.05;
    diffuse_cutoffs_["C"][2] = 0.17;

    diffuse_cutoffs_["N"][0] = 0.08;
    diffuse_cutoffs_["N"][1] = 0.07;
    diffuse_cutoffs_["N"][2] = 0.3;

    diffuse_cutoffs_["O"][0] = 0.1;
    diffuse_cutoffs_["O"][1] = 0.08;
    diffuse_cutoffs_["O"][2] = 0.4;

    diffuse_cutoffs_["F"][0] = 0.12;
    diffuse_cutoffs_["F"][1] = 0.10;
    diffuse_cutoffs_["F"][2] = 0.5;

    diffuse_cutoffs_["NE"][0] = 0.15;
    diffuse_cutoffs_["NE"][1] = 0.13;
    diffuse_cutoffs_["NE"][2] = 0.7;

    statics_done_ = true;
}

Shell::diffuse_type_t
Shell::diffuse_type() const
{
    return diffuse_type_;
}

uli
Shell::nfxn() const
{
    return nfxn_;
}

void
Shell::print(std::ostream& os) const
{
    os << stream_printf("Shell %d: atom=%d am=%d ncxn=%d",
                    this->program_shell_number_,
                    atom_->number(),
                    ang_mom_,
                    nsets_);
    if (diffuse_type_ == valence)
        os << " valence";
    else
        os << " diffuse";
}

const AtomPtr&
Shell::get_atom() const
{
    return atom_;
}

usi
Shell::angular_momentum() const
{
    return ang_mom_;
}

uli
Shell::program_shell_number() const
{
    return program_shell_number_;
}

uli
Shell::program_fxn_start() const
{
    return program_fxn_start_;
}

Atom::Atom(
    const std::string& symbol,
    const std::string& name,
    double x,
    double y,
    double z,
    uli atom_number
) : symbol_(symbol),
    name_(name),
    x_(x),
    y_(y),
    z_(z),
    number_(atom_number),
    group_number_(0)
{
    std::transform(symbol_.begin(), symbol_.end(), symbol_.begin(), myupper);
    std::transform(name_.begin(), name_.end(), name_.begin(), myupper);
}

uli
Atom::group_number() const
{
    return group_number_;
}

uli
Atom::number() const
{
    return number_;
}

const std::string&
Atom::name() const
{
    return name_;
}

const std::string&
Atom::symbol() const
{
    return symbol_;
}

void
Atom::set_group_number(uli num)
{
    group_number_ = num;
}

AOBasis::AOBasis(const std::string& name, const std::string& descr)
: name_(name),
    descr_(descr),
    nbasis_(0),
    range_(0),
    index_descr_(0),
    nfxn_min_per_tile_(0),
    multi_shell_map_(0),
    yeti_to_program_shell_map_(0),
    yeti_to_program_fxn_map_(0),
    program_to_yeti_shell_map_(0),
    program_to_yeti_fxn_map_(0)
{
}

void
AOBasis::add_atom(const AtomPtr& atom)
{
    atoms_.push_back(atom);
}

void
AOBasis::add_shell(const ShellPtr& shell)
{
    nbasis_ += shell->nfxn();
    shells_.push_back(shell);
}

void
AOBasis::configure_atom_groups_as_atoms()
{
    uli natoms = atoms_.size();
    for (uli i=0; i < natoms; ++i)
        atoms_[i]->set_group_number(i);
}

void
AOBasis::configure_atom_groups_by_number(uli natoms_per_group)
{
    uli natoms = atoms_.size();

    uli ngroups = natoms / natoms_per_group;
    if (ngroups == 0)
        ngroups = 1;

    uli nper = natoms / ngroups;
    uli nextra = natoms % nper;
    uli groupnum = 0;
    uli atomnum = 0;
    for (uli i=0; i < nextra; ++i, ++groupnum)
    {
        uli natoms_in_group = nper + 1;
        for (uli j=0; j < natoms_in_group; ++j, ++atomnum)
        {
            atoms_[atomnum]->set_group_number(groupnum);
        }
    }
    for (uli i=nextra; i < ngroups; ++i, ++groupnum)
    {
        uli natoms_in_group = nper;
        for (uli j=0; j < natoms_in_group; ++j, ++atomnum)
        {
            atoms_[atomnum]->set_group_number(groupnum);
        }
    }
}

void
AOBasis::configure_many_atoms_large_basis(uli natoms_per_group)
{
    configure_atom_groups_by_number(natoms_per_group);

    //create the bottom level group - determine the number groups
    std::vector<ShellPtr>::const_iterator it = shells_.begin() + 1;
    std::vector<ShellPtr>::const_iterator stop = shells_.end();

    std::vector<uli> data_partitions;
    std::vector<uli> node_partitions;
    std::vector<uli> thread_partitions;

    uli thread_start = 0;
    uli node_start = 0;
    uli nelements_data = shells_[0]->nfxn();
    usi am = shells_[0]->angular_momentum();
    if (am == 1) //treat s and p equivalently
        am = 0;
    Shell::diffuse_type_t difftype = shells_[0]->diffuse_type();
    uli groupnum = shells_[0]->get_atom()->group_number();
    uli atomnum = shells_[0]->get_atom()->number();
    for ( ; it != stop; ++it)
    {
        const ShellPtr& sh = *it;
        usi am_shell = sh->angular_momentum();
        if (am_shell == 1)
            am_shell = 0;

        Shell::diffuse_type_t difftype_shell = sh->diffuse_type();

        uli groupnum_shell = sh->get_atom()->group_number();

        uli atomnum_shell = sh->get_atom()->number();

        if (am_shell != am || difftype_shell != difftype || groupnum_shell != groupnum)
        {
            data_partitions.push_back(nelements_data);

            uli nranges_thread = data_partitions.size()  - thread_start;
            thread_partitions.push_back(nranges_thread);
            thread_start += nranges_thread;

            uli nranges_node = thread_partitions.size() - node_start;
            node_partitions.push_back(nranges_node);
            node_start += nranges_node;

            nelements_data = sh->nfxn();
            am = am_shell;
            difftype = difftype_shell;
            groupnum = groupnum_shell;
            atomnum = atomnum_shell;
        }
        else if (atomnum != atomnum_shell)
        {
            data_partitions.push_back(nelements_data);

            uli nranges_thread = data_partitions.size()  - thread_start;
            thread_partitions.push_back(nranges_thread);
            thread_start += nranges_thread;

            nelements_data = sh->nfxn();
            atomnum = atomnum_shell;
        }
        else
        {
            nelements_data += sh->nfxn();
        }
    }
    //add the last one
    data_partitions.push_back(nelements_data);

    uli nranges_thread = data_partitions.size()  - thread_start;
    thread_partitions.push_back(nranges_thread);

    uli nranges_node = thread_partitions.size() - node_start;
    node_partitions.push_back(nranges_node);

    uli start = 0;
    IndexRangePtr data_range = new IndexRange(start, data_partitions);

    IndexRangePtr thread_range = new IndexRange(start, thread_partitions);
    thread_range->acquire_subranges(data_range->get_subranges());

    range_ = new IndexRange(start, node_partitions);
    range_->acquire_subranges(thread_range->get_subranges());

    if (nfxn_min_per_tile_)
    {
        range_ = range_->squeeze_together_bottom_ranges(nfxn_min_per_tile_);
    }


}

void
AOBasis::configure_many_atoms_large_basis()
{
    uli nshells = shells_.size();
    for (uli i=0; i < nshells; ++i)
    {
        //shells_[i]->print(cout); cout << endl;
    }

    std::sort(shells_.begin(), shells_.end(), _lt_shells_am_diffusetype_atomgroup_atom);

    //cout << endl;
    for (uli i=0; i < nshells; ++i)
    {
        //shells_[i]->print(cout); cout << endl;
    }

    if (yeti_to_program_shell_map_ == 0)
        yeti_to_program_shell_map_ = new uli[nshells];
    if (yeti_to_program_fxn_map_ == 0)
        yeti_to_program_fxn_map_ = new uli[nbasis_];
    if (program_to_yeti_shell_map_ == 0)
        program_to_yeti_shell_map_ = new uli[nshells];
    if (program_to_yeti_fxn_map_ == 0)
        program_to_yeti_fxn_map_ = new uli[nbasis_];


    uli fxn = 0;
    for (uli i=0; i < nshells; ++i)
    {
        uli program_number = shells_[i]->program_shell_number();
        yeti_to_program_shell_map_[i] = program_number;
        program_to_yeti_shell_map_[program_number] = i;

        uli nfxn_in_shell = shells_[i]->nfxn();
        uli fxn_start = shells_[i]->program_fxn_start();
        for (uli j=0; j < nfxn_in_shell; ++j, ++fxn)
        {
            yeti_to_program_fxn_map_[fxn] = fxn_start + j;
            program_to_yeti_fxn_map_[fxn_start + j] = fxn;
        }
    }


    for (uli i=0; i < nbasis_; ++i)
    {
        //cout << stream_printf("yeti %d -> mpqc %d", i, yeti_to_program_fxn_map_[i]) << endl;
    }

    uli min_nelements = 5;
    nshells = 1;
    uli atomnum = shells_[0]->get_atom()->number();
    for (uli i=1; i < shells_.size(); ++i)
    {
        uli atomnum_shell = shells_[i]->get_atom()->number();
        if (atomnum_shell != atomnum)
        {
            atomnum = atomnum_shell;
            ++nshells;
        }
    }
    if (min_nelements > nshells)
        min_nelements = nshells;

    uli natoms_per_group = 5;
    configure_many_atoms_large_basis(natoms_per_group);


    while (range_->nelements() < min_nelements)
    {
        range_ = 0;
        --natoms_per_group;
        configure_many_atoms_large_basis(natoms_per_group);
    }




    YetiRuntime::register_index_range(range_, descr_, name_);
    index_descr_ = YetiRuntime::get_descr(name_);
    multi_shell_map_ = new MultiShellMap(shells_, index_descr_);
}

void
AOBasis::configure_few_atoms_large_basis()
{
    configure_many_atoms_large_basis();
    std::sort(shells_.begin(), shells_.end(), _lt_shells_am_diffusetype_atomgroup_atom);
}

void
AOBasis::configure_many_atoms_small_basis()
{
    configure_many_atoms_large_basis();
    std::sort(shells_.begin(), shells_.end(), _lt_shells_am_diffusetype_atomgroup_atom);
}

void
AOBasis::configure_few_atoms_small_basis()
{
    configure_many_atoms_large_basis();
    std::sort(shells_.begin(), shells_.end(), _lt_shells_am_diffusetype_atomgroup_atom);
}

void
AOBasis::configure_min_nfxn_per_tile(uli nmin)
{
    nfxn_min_per_tile_ = nmin;

    if (range_)
        range_->squeeze_together_bottom_ranges(nfxn_min_per_tile_);
}

uli
AOBasis::program_shell_number(uli yeti_shell_number)
{
    return yeti_to_program_shell_map_[yeti_shell_number];
}

const std::string&
AOBasis::name() const
{
    return name_;
}

const std::string&
AOBasis::descr() const
{
    return descr_;
}

const MultiShellMapPtr&
AOBasis::get_multishell_map() const
{
    return multi_shell_map_;
}

IndexDescr*
AOBasis::get_index_descr() const
{
    return index_descr_.get();
}

void
AOBasis::convert_to_yeti_numbering(RectMatrixPtr& coefs)
{
    double* Ctmp = new double[nbasis_];
    for (uli mo=0; mo < coefs.ncol(); ++mo)
    {
        for (uli yetifxn=0; yetifxn < nbasis_; ++yetifxn)
        {
            uli progfxn = yeti_to_program_fxn_map_[yetifxn];
            double Cxi = coefs.get_element(progfxn, mo);
            Ctmp[yetifxn] = Cxi;
        }

        for (uli yetifxn=0; yetifxn < nbasis_; ++yetifxn)
            coefs.set_element(yetifxn, mo, Ctmp[yetifxn]);
    }
    delete[] Ctmp;
}

MultiShellMap::MultiShellMap(
    const std::vector<ShellPtr>& shells,
    const IndexDescrPtr& ao_descr
) :
    shell_sizes_(0),
    shell_start_on_index_(0),
    nshells_per_index_(0),
    nfxn_per_index_(0)
{
    uli nshell = shells.size();
    usi metadata_depth = 1;
    usi nranges = ao_descr->nelements(metadata_depth);
    shell_start_on_index_ = new uli[nranges];
    nshells_per_index_ = new uli[nranges];
    shell_sizes_ = new uli[nshell];
    nfxn_per_index_ = new uli[nranges];

    uli idx = 0;
    uli bsfxn = 0;
    uli nsh_on_idx = 0;
    uli nfxn_on_idx = ao_descr->nelements(idx);
    uli start = 0;
    for (uli sh=0; sh < nshell; ++sh)
    {
        uli nfxn_on_shell = shells[sh]->nfxn();
        bsfxn += nfxn_on_shell;
        ++nsh_on_idx;
        shell_sizes_[sh] = nfxn_on_shell;
        if (bsfxn == nfxn_on_idx)
        {
            nfxn_per_index_[idx] = nfxn_on_idx;
            shell_start_on_index_[idx] = start;
            start = sh + 1;
            nshells_per_index_[idx] = nsh_on_idx;
            bsfxn = 0;
            nsh_on_idx = 0;
            ++idx;
            if (idx < nranges)
                nfxn_on_idx = ao_descr->nelements(idx);
        }
    }
}

uli
MultiShellMap::nshell(uli idx)
{
    return nshells_per_index_[idx];
}

uli
MultiShellMap::nfxn(uli idx)
{
    return nfxn_per_index_[idx];
}

uli
MultiShellMap::shell_size(uli sh)
{
    return shell_sizes_[sh];
}

uli
MultiShellMap::shell_start(uli idx)
{
    return shell_start_on_index_[idx];
}
