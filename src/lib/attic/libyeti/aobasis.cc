#include "aobasis.h"
#include "index.h"
#include "runtime.h"
#include "gigmatrix.h"
#include "exception.h"
#include "env.h"

#include <algorithm>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define MAX_AM 10

#define PRINT_SHELLS 0

static Shell* null_shell = 0;

bool Shell::statics_done_ = false;
std::map<std::string, std::map<usi, double> > Shell::diffuse_cutoffs_;

int myupper(int c)
{
    return std::toupper((unsigned char)c);
}


#define fix_intra_idx(am,x) \
if (am > 1) x += (am - 1);
//else if (am == 2 && x== 0) x = 1;

#define OLD_SORT 0

#if OLD_SORT
#define fix_am(x) if (x <= 1) x = 0
DEPRECATED bool
_lt_shells_amgroups_diffusetype_atomgroup_atom(
    const ShellPtr& sh1,
    const ShellPtr& sh2
)
{
    usi am1 = sh1->angular_momentum();
    usi am2 = sh2->angular_momentum();

    // Group s and p together...
    if (am1 == 1) am1 = 0;
    if (am2 == 1) am2 = 0;

    // Group everything else together
    if (am1 > 2) am1 = 2;
    if (am2 > 2) am2 = 2;

    if (am1 != am2)
        return am1 < am2;

    Shell::diffuse_type_t diff1 = sh1->diffuse_type();
    Shell::diffuse_type_t diff2 = sh2->diffuse_type();
    if (diff1 != diff2)
        return diff1 < diff2;

    uli mod1 = sh1->get_atom()->modulo_number();
    uli mod2 = sh2->get_atom()->modulo_number();
    if (mod1 != mod2)
        return mod1 < mod2;

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
#else
#define fix_am(x) if (x <= 2) x = 0
// This sorts stuff
DEPRECATED bool
_lt_shells_am_diffusetype_atomgroup_atom(
    const ShellPtr& sh1,
    const ShellPtr& sh2
)
{
    usi am1 = sh1->angular_momentum();
    usi am2 = sh2->angular_momentum();
    uli idx1 = sh1->intra_am_index();
    uli idx2 = sh2->intra_am_index();

    fix_intra_idx(am1,idx1);
    fix_intra_idx(am2,idx2);
    //fix_am(am1);
    //fix_am(am2);

    if (idx1 != idx2)
        return idx1 < idx2;

    Shell::am_type_t amtype1 = sh1->am_type();
    Shell::am_type_t amtype2 = sh2->am_type();
    if (amtype1 != amtype2)
        return amtype1 < amtype2;


    Shell::diffuse_type_t diff1 = sh1->diffuse_type();
    Shell::diffuse_type_t diff2 = sh2->diffuse_type();
    if (diff1 != diff2)
        return diff1 < diff2;

    uli mod1 = sh1->get_atom()->modulo_number();
    uli mod2 = sh2->get_atom()->modulo_number();
    if (mod1 != mod2)
        return mod1 < mod2;

    uli group1 = sh1->get_atom()->group_number();
    uli group2 = sh2->get_atom()->group_number();
    if (group1 != group2)
        return group1 < group2;

    uli at1 = sh1->get_atom()->number();
    uli at2 = sh2->get_atom()->number();
    if (at1 != at2)
        return at1 < at2;

    if (am1 != am2)
        return am1 < am2;

    //exactly the same... so... false...
    return false;
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////
//  Shell class
////////////////////////////////////////////////////////////////////////////////////////////////////


Shell::Shell(
    const AtomPtr& atom,
    usi ang_mom,
    uli nsets,
    uli nfxn,
    double min_exponent,
    uli program_shell_number,
    uli program_fxn_start,
    uli intra_am_idx
) :  atom_(atom),
    ang_mom_(ang_mom),
    nsets_(nsets),
    nfxn_(nfxn),
    min_exponent_(min_exponent),
    program_shell_number_(program_shell_number),
    program_fxn_start_(program_fxn_start),
    diffuse_type_(Shell::valence),
    intra_am_index_(intra_am_idx),
    am_type_(Shell::low_am)
{
    if (!statics_done_)
        init_statics();

    double diff_cutoff = diffuse_cutoffs_[atom->symbol()][ang_mom];
    if (min_exponent < diff_cutoff)
        diffuse_type_ = Shell::diffuse;

    if (ang_mom >= 3)
        am_type_ = Shell::high_am;
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

    //diffuse_cutoffs_["NE"][0] = 0.15;
    //diffuse_cutoffs_["NE"][1] = 0.13;
    //diffuse_cutoffs_["NE"][2] = 0.7;
    //diffuse_cutoffs_["NE"][3] = 0.7;
    //diffuse_cutoffs_["NE"][4] = 0.7;

    diffuse_cutoffs_["NE"][0] = 0.45;
    diffuse_cutoffs_["NE"][1] = 0.36;
    diffuse_cutoffs_["NE"][2] = 0.7;
    diffuse_cutoffs_["NE"][3] = 0.7;


    statics_done_ = true;
}

Shell::am_type_t
Shell::am_type() const
{
    return am_type_;
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
    os << stream_printf("Shell %2d: atom=%d am=%d ncxn=%d mod=%d grp=%d, intra=%d, nfxn=%2d",
                    this->program_shell_number_,
                    atom_->number(),
                    ang_mom_,
                    nsets_,
                    atom_->modulo_number(),
                    atom_->group_number(),
                    intra_am_index_,
                    nfxn_
                );
    if (am_type_ == high_am)
        os << " high";
    else
        os << " low";

    if (diffuse_type_ == valence)
        os << " valence";
    else
        os << " diffuse";
    os << endl;
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

uli
Shell::intra_am_index() const
{
    return intra_am_index_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//  Atom class
////////////////////////////////////////////////////////////////////////////////////////////////////

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
    group_number_(0),
    modulo_number_(0)
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
Atom::modulo_number() const
{
    return modulo_number_;
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

void
Atom::set_modulo_number(uli num)
{
    modulo_number_ = num;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//  AOBasis class
////////////////////////////////////////////////////////////////////////////////////////////////////

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


std::vector<uli>
_repartition_smaller(
    std::vector<uli> shell_sizes,
    uli max_elements
)
{
    // TODO: This could actually be done much better by rearranging shells to create more even partitions, but the general solution to that problem is NP-complete, so it would take exponential time.  Still, could be fun...
    std::vector<uli> newparts(0);
    uli nelements_current = 0;
    for(int i = 0; i < shell_sizes.size(); ++i) {
        if(nelements_current + shell_sizes[i] > max_elements) {
            newparts.push_back(nelements_current);
            nelements_current = 0;
        }

        nelements_current += shell_sizes[i];

    }
    // Now push on the last one...
    if(nelements_current >= max_elements) {
        newparts.push_back(nelements_current);
    }
    else {
        // Oh well, we tried...
        newparts[newparts.size() - 1] += nelements_current;
    }

    return newparts;

}

void
AOBasis::configure(
    PartitioningPolicyPtr policy
)
{
    uli min_nelements = 2;
    uli nshells = 1;
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

    uli natoms_per_group = 2;

    // This will need to be replaced with a more reasonble grouping
    // scheme at some point.
    configure_atom_groups_by_number(natoms_per_group);

    // Sort the shells so that like groups are together
    policy->sort_shells(shells_);

    //create the bottom level group - determine the number groups
    std::vector<ShellPtr>::const_iterator it = shells_.begin() + 1;
    std::vector<ShellPtr>::const_iterator stop = shells_.end();

    std::vector< std::vector<uli> > partitions(0);
    std::vector< uli > nfxns_so_far(0);
    std::vector< std::vector<IndexRangePtr> > ranges(0);
    std::vector< uli > current_shell_sizes(0);
    int nlevels = policy->nlevels();
    uli starts[nlevels], subrngstarts[nlevels];
    for(int i = 0; i < nlevels; ++i) {
        starts[i] = 0;
        subrngstarts[i] = 0;
        std::vector<uli> ipt(0);
        nfxns_so_far.push_back(0);
        std::vector<IndexRangePtr> rng(0);
        partitions.push_back(ipt);
        ranges.push_back(rng);
        Env::out0() << Env::indent << "{" << endl;
        ++Env::indent;
    }
    uli nelements_data = shells_[0]->nfxn();
    current_shell_sizes.push_back(shells_[0]->nfxn());

    ShellPtr prev = shells_[0];
    for ( ; it != stop; ++it)
    {
        const ShellPtr& sh = *it;

        if(policy->should_partition_at_level(prev, sh, 0, nelements_data, nfxns_so_far)) {

            // Check if the partition is too big.  If so, split it.
            if(nelements_data > policy->max_fxn_per_tile()) {
                // Repartition smaller
                std::vector<uli> newparts = _repartition_smaller(current_shell_sizes, policy->max_fxn_per_tile());

                Env::out0() << --Env::indent << "} (" << nelements_data << " fxns split into " << newparts.size() << " partitions: ";

                // And add the partitions
                for(int i = 0; i < newparts.size(); ++i) {
                    partitions[0].push_back(newparts[i]);
                    ranges[0].push_back(new IndexRange(starts[0], newparts[i]));
                    starts[0] += newparts[i];
                    Env::out0() << newparts[i] << (i == newparts.size() - 1 ? " fxns, " : " fxns");
                }

                Env::out0() << ")" << endl;

            }
            else {
                partitions[0].push_back(nelements_data);
                ranges[0].push_back(new IndexRange(starts[0], nelements_data));
                starts[0] += nelements_data;
                --Env::indent;
                Env::out0() << Env::indent << "}" << endl;

            }

            nfxns_so_far[0] = 0;
            current_shell_sizes.clear();
            nelements_data = 0;


            int nunindents = 0;
            for(int i = 1; i < nlevels; ++i) {
                uli nrangesi = partitions[i-1].size() - starts[i];

                if(policy->should_partition_at_level(prev, sh, i, nrangesi, nfxns_so_far)) {

                    SubindexTuplePtr subranges = new SubindexTuple(nrangesi);
                    for(int j = starts[i]; j < starts[i] + nrangesi; ++j) {
                        subranges->insert(j - starts[i], ranges[i-1][j].get());
                    }
                    ranges[i].push_back(new IndexRange(starts[i], subranges));

                    nfxns_so_far[i] = 0;
                    partitions[i].push_back(nrangesi);
                    starts[i] += nrangesi;
                    nunindents++;
                    --Env::indent;
                    Env::out0() << Env::indent << "}" << endl;
                }
            }
            for(int a = 0; a < nunindents; ++a) {
                Env::out0() << Env::indent << "{" << endl;
                ++Env::indent;
            }

            Env::out0() << Env::indent << "{" << endl;
            ++Env::indent;

        }



        nelements_data += sh->nfxn();
        current_shell_sizes.push_back(sh->nfxn());
        for(int i = 0; i < nlevels; ++i)
            nfxns_so_far[i] += sh->nfxn();


#if PRINT_SHELLS
        sh->print(Env::out0() << Env::indent);
#endif
        prev = sh;

    }

    //Now add the last ones
    if(nelements_data > policy->max_fxn_per_tile()) {
        // Repartition smaller
        std::vector<uli> newparts = _repartition_smaller(current_shell_sizes, policy->max_fxn_per_tile());

        Env::out0() << --Env::indent << "} (" << nelements_data << " fxns split into " << newparts.size() << " partitions: ";

        // And add the partitions
        for(int i = 0; i < newparts.size(); ++i) {
            partitions[0].push_back(newparts[i]);
            ranges[0].push_back(new IndexRange(starts[0], newparts[i]));
            starts[0] += newparts[i];
            Env::out0() << newparts[i] << (i == newparts.size() - 1 ? " fxns" : " fxns, ");
        }
        Env::out0() << ")" << endl;
    }
    else {
        partitions[0].push_back(nelements_data);
        ranges[0].push_back(new IndexRange(starts[0], nelements_data));
        starts[0] += nelements_data;
        Env::out0() << --Env::indent << "}" << endl;
    }

    for(int i = 1; i < nlevels; ++i) {
        uli nrangesi = partitions[i-1].size() - starts[i];
        partitions[i].push_back(nrangesi);

        SubindexTuplePtr subranges = new SubindexTuple(nrangesi);
        for(int j = starts[i]; j < starts[i] + nrangesi; ++j) {
            subranges->insert(j - starts[i], ranges[i-1][j].get());
        }
        ranges[i].push_back(new IndexRange(starts[i], subranges));

        Env::out0() << --Env::indent << "}" << endl;
    }


    // Now split the tiles with more than the maximum number of functions
    /*for(int parent = 0; parent < partitions[1].size(); ++parent) {
        bool repartition = false;
        for(int i = 0; i < partitions[1][parent]; ++i) {

        }
    }*/
    /*
    for(int i = 0; i < partitions[0].size(); ++i) {
        while(ranges[i]->n)
    }
    */


    /*
      Old code that might be useful again soon.


    // Do the top-level combining first
    int lvl = nlevels - 1;
    int minfxn = policy->min_fxns_at_level(lvl);
    for(int i = 0; i < partitions[lvl].size() - 1; ++i) {
        if(nfxns_in_ranges[lvl][i] < minfxn || nfxns_in_ranges[lvl][i + 1] < minfxn) {
            uli nidx_next = partitions[lvl][i+1];
            partitions[lvl].erase(i+1);
            partitions[lvl][i] += nidx_next;

            uli nfxn_next = nfxns_in_ranges[lvl][i+1];
            nfxns_in_ranges[lvl].erase(i+1);
            nfxns_in_ranges[lvl][i] += nfxn_next;

            i--;  // Decrement to account for the increment causing the size to be checked again.
        }
    }

    // Now check the subranges of each range for the need to combine
    for(int lvl = nlevels-1; lvl > 0; --lvl) {
        int start_pos = 0;
        int minfxn = policy->min_fxns_at_level(lvl-1);
        for(int i = 0; i < partitions[lvl].size(); ++i) {
            for(int j = start_pos; j < start_pos + partitions[lvl][i] - 1; ++j) {
                if(nfxns_in_ranges[lvl-1][j] < minfxn || nfxns_in_ranges[lvl-1][j + 1] < minfxn) {
                    uli nidx_next = partitions[lvl-1][j+1];
                    partitions[lvl-1].erase(j+1);
                    partitions[lvl-1][j] += nidx_next;

                    uli nfxn_next = nfxns_in_ranges[lvl][j+1];
                    nfxns_in_ranges[lvl-1].erase(j+1);
                    nfxns_in_ranges[lvl-1][j] += nfxn_next;

                    j--;  // Decrement to account for the increment causing the size to be checked again.
                    partitions[lvl][i]--;
                }
            }
            start_pos += partitions[lvl][i];
        }
    }


    // Now check all the way back up until we don't need to do any combining anymore
    for(int lvl = 1; lvl < nlevels; ++lvl) {
        int start_pos = 0;
        int minfxn = policy->min_fxns_at_level(lvl-1);
        for(int ipart = 0; ipart < partitions[lvl].size() - 1; ++i) {
            // If one of the first two is not true, then all of the subranges are already happy
            if(partitions[lvl][ipart] == 1) {
                if(partitions[lvl-1][start_pos] < min)
            }

            start_pos += partitions[lvl][ipart];
        }
    }

    IndexRangePtr old_range = new IndexRange(0, partitions[0]);
    for(int i = 1; i < nlevels - 1; ++i) {
        IndexRangePtr i_range = new IndexRange(0, partitions[i]);
        i_range->acquire_subranges(old_range->get_subranges());
        old_range = i_range;
    }
    range_ = new IndexRange(0, partitions[nlevels - 1]);
    range_->acquire_subranges(old_range->get_subranges());
    */

    SubindexTuplePtr subranges = new SubindexTuple(partitions[nlevels - 1].size());
    for(int i = 0; i < ranges[nlevels - 1].size(); ++i) {
        subranges->insert(i, ranges[nlevels-1][i].get());
    }

    range_ = new IndexRange(0, subranges);

    ////////////////////////////////////////////////////////////////////////////////////////////////////


    nshells = shells_.size();

    //for (uli i=0; i < nshells; ++i)
    //{
        //shells_[i]->print(cout); cout << endl;
    //}
    //cout << endl;

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

    YetiRuntime::register_index_range(range_, descr_, name_);
    index_descr_ = YetiRuntime::get_descr(name_);
    multi_shell_map_ = new MultiShellMap(shells_, index_descr_);

}


void
AOBasis::configure_atom_groups_as_atoms()
{
    uli natoms = atoms_.size();
    for (uli i=0; i < natoms; ++i)
    {
        atoms_[i]->set_group_number(i);
        atoms_[i]->set_modulo_number(0);
    }
}

void
AOBasis::configure_atom_groups_by_number(uli natoms_per_group)
{
    uli natoms = atoms_.size();

    uli ngroups = natoms / natoms_per_group;
    if (ngroups == 0)
        ngroups = 1;

    if (natoms % ngroups)
        ++ngroups;


    uli nper = natoms / ngroups;
    uli nextra = natoms % ngroups;
    uli groupnum = 0;
    uli atomnum = 0;
    for (uli i=0; i < nextra; ++i, ++groupnum)
    {
        uli natoms_in_group = nper + 1;
        for (uli j=0; j < natoms_in_group; ++j, ++atomnum)
        {
            atoms_[atomnum]->set_group_number(groupnum);
            atoms_[atomnum]->set_modulo_number(j);
        }
    }
    for (uli i=nextra; i < ngroups; ++i, ++groupnum)
    {
        uli natoms_in_group = nper;
        for (uli j=0; j < natoms_in_group; ++j, ++atomnum)
        {
            atoms_[atomnum]->set_group_number(groupnum);
            atoms_[atomnum]->set_modulo_number(j);
        }
    }
}

void
AOBasis::configure_many_atoms_large_basis(uli natoms_per_group)
{
    configure_atom_groups_by_number(natoms_per_group);

    // Sort the shells so that like groups are together
    std::sort(shells_.begin(), shells_.end(), _lt_shells_am_diffusetype_atomgroup_atom);

    //create the bottom level group - determine the number groups
    std::vector<ShellPtr>::const_iterator it = shells_.begin() + 1;
    std::vector<ShellPtr>::const_iterator stop = shells_.end();


    std::vector<uli> data_partitions;
    std::vector<uli> thread_partitions;
    std::vector<uli> node_partitions;

    uli node_start = 0;
    uli thread_start = 0;
    uli nelements_data = shells_[0]->nfxn();
    uli intra_idx = shells_[0]->intra_am_index();
    usi am = shells_[0]->angular_momentum();
    Shell::diffuse_type_t difftype = shells_[0]->diffuse_type();
    Shell::am_type_t amtype = shells_[0]->am_type();
    uli groupnum = shells_[0]->get_atom()->group_number();
    uli atomnum = shells_[0]->get_atom()->number();
    uli modulonum = shells_[0]->get_atom()->modulo_number();
    fix_intra_idx(am,intra_idx);
    fix_am(am);
    //shells_[0]->print(cout); cout << endl;
    for ( ; it != stop; ++it)
    {
        const ShellPtr& sh = *it;
        //sh->print(cout); cout << endl;

        usi am_shell = sh->angular_momentum();

        Shell::diffuse_type_t difftype_shell = sh->diffuse_type();

        Shell::am_type_t amtype_shell = sh->am_type();

        uli groupnum_shell = sh->get_atom()->group_number();

        uli modulonum_shell = sh->get_atom()->modulo_number();

        uli atomnum_shell = sh->get_atom()->number();

        uli intra_idx_shell = sh->intra_am_index();

        fix_intra_idx(am_shell, intra_idx_shell);
        fix_am(am_shell);

#if OLD_SORT
        if (am_shell != am || difftype_shell != difftype || modulonum_shell != modulonum)
        {

            // We need to check if any of the other levels will be split.  This should be done with a little more foresight in the future...
            if (nelements_data >= nfxn_min_per_tile_ || difftype_shell != difftype || am_shell != am)
            {
                data_partitions.push_back(nelements_data);

                nelements_data = sh->nfxn();

                if(difftype_shell != difftype || am_shell != am) {
                    uli nranges_thread = data_partitions.size()  - diffuse_start;
                    diffuse_partitions.push_back(nranges_thread);
                    diffuse_start += nranges_thread;

                    if(am_shell != am)
                    {
                        uli nranges_node = diffuse_partitions.size() - angmom_start;
                        angmom_partitions.push_back(nranges_node);
                        angmom_start += nranges_node;
                        //cout << "-\n--\n---\n--" << endl;
                    }
                    else {
                        //cout << "-\n--" << endl;
                    }


                }

                //cout << "-" << endl;
            }
            else
            {
                nelements_data += sh->nfxn();
            }

            // Now set all of the values for the next check
            am = am_shell;
            difftype = difftype_shell;
            groupnum = groupnum_shell;
            modulonum = modulonum_shell;
            atomnum = atomnum_shell;
        }
#else
        if (difftype_shell != difftype 
            //|| am_shell != am
            || intra_idx_shell != intra_idx
            || modulonum_shell != modulonum
            || amtype_shell != amtype
        )
        {
            if (nelements_data >= nfxn_min_per_tile_)
            {
                data_partitions.push_back(nelements_data);

                uli nranges_thread = data_partitions.size()  - thread_start;
                thread_partitions.push_back(nranges_thread);
                thread_start += nranges_thread;

                uli nranges_node = thread_partitions.size() - node_start;
                node_partitions.push_back(nranges_node);
                node_start += nranges_node;

                nelements_data = sh->nfxn();
            }
            else
            {
                //cout << stream_printf("n=%d nmin=%d", nelements_data, nfxn_min_per_tile_) << endl;
                nelements_data += sh->nfxn();
            }



            am = am_shell;
            difftype = difftype_shell;
            amtype = amtype_shell;
            intra_idx = intra_idx_shell;
            groupnum = groupnum_shell;
            modulonum = modulonum_shell;
            atomnum = atomnum_shell;

        }
        else if (groupnum != groupnum_shell)
        {
            if (nelements_data >= nfxn_min_per_tile_)
            {
                data_partitions.push_back(nelements_data);
                nelements_data = sh->nfxn();
            }
            else
            {
                nelements_data += sh->nfxn();
            }

            groupnum = groupnum_shell;
            atomnum = atomnum_shell;

        }
#endif
        else
        {
#if PRINT_SHELLS
            sh->print(Env::out0());
#endif
            nelements_data += sh->nfxn();
        }
    }
    //cout << "-\n--\n---" << endl;
    //add the last one
    data_partitions.push_back(nelements_data);

    uli nranges_thread = data_partitions.size()  - thread_start;
    thread_partitions.push_back(nranges_thread);

    uli nranges_node = thread_partitions.size() - node_start;
    node_partitions.push_back(nranges_node);

#if 0
    for(int i = 0; i < data_partitions.size(); ++i) {
        cout << "data_partitions[" << i << "] = " << data_partitions[i] << endl;
    }

    for(int i = 0; i < thread_partitions.size(); ++i) {
        cout << "thread_partitions[" << i << "] = " << thread_partitions[i] << endl;
    }
    for(int i = 0; i < node_partitions.size(); ++i) {
        cout << "node_partitions[" << i << "] = " << node_partitions[i] << endl;
    }
#endif

    uli start = 0;
    IndexRangePtr data_range = new IndexRange(start, data_partitions);

    IndexRangePtr thread_range = new IndexRange(start, thread_partitions);
    thread_range->acquire_subranges(data_range->get_subranges());

    range_ = new IndexRange(start, node_partitions);
    range_->acquire_subranges(thread_range->get_subranges());
}

void
AOBasis::configure_many_atoms_large_basis()
{

    uli min_nelements = 2;
    uli nshells = 1;
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
    //if (min_nelements < nshells)
        min_nelements = 0; //just forget it

    uli natoms_per_group = 3;
    configure_many_atoms_large_basis(natoms_per_group);

#if 0
    while (range_->nelements() < min_nelements)
    {
        range_ = 0;
        ++natoms_per_group; //more atoms per group generates more blocks
        if (natoms_per_group > atoms_.size())
            yeti_throw(SanityCheckError, "Yeti AO Basis configuration error");
        configure_many_atoms_large_basis(natoms_per_group);
    }
#endif


    nshells = shells_.size();

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

IndexRange*
AOBasis::get_index_range() const
{
    return range_.get();
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


////////////////////////////////////////////////////////////////////////////////////////////////////
//  MultiShellMap class
////////////////////////////////////////////////////////////////////////////////////////////////////

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
        else if(bsfxn > nfxn_on_idx) {
            yeti_throw(SanityCheckError, "Shell misalignment or something else weird.");
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



////////////////////////////////////////////////////////////////////////////////////////////////////
//  PartitioningPolicy class
////////////////////////////////////////////////////////////////////////////////////////////////////


PartitioningPolicy::PartitioningPolicy()
    : levels_(0),
      sort_importance_(0),
      sizecheck_traits_(0),
      nfxn_min_per_tile_(0),
      nidx_min_per_tile_(0)
{
    // This is the original policy implemented by Jeremy.
    levels_.push_back(NoPartition);
    levels_.push_back(NoPartition);
    levels_.push_back(SmallOrLargeAngularMomentum|DiffuseType|ModuloNumber|GroupNumber);

    sort_importance_.push_back(SmallOrLargeAngularMomentum);
    sort_importance_.push_back(DiffuseType);
    sort_importance_.push_back(ModuloNumber);
    sort_importance_.push_back(GroupNumber);
    sort_importance_.push_back(AtomNumber);


    sizecheck_traits_.push_back(NoPartition);
    sizecheck_traits_.push_back(NoPartition);
    sizecheck_traits_.push_back(GroupNumber);

#if YETI_SANITY_CHECK || 1
    check_sort_importance();
#endif

    init();

}


PartitioningPolicy::PartitioningPolicy(
    std::vector<int> levels
) : levels_(levels),
    sort_importance_(0),
    sizecheck_traits_(0),
    nfxn_min_per_tile_(0),
    nidx_min_per_tile_(0)
{
    build_sort_importance();

#if YETI_SANITY_CHECK || 1
    check_sort_importance();
#endif

    init();
}


PartitioningPolicy::PartitioningPolicy(
    int d1,
    int d2,
    int d3,
    int d4,
    int d5,
    int d6,
    int d7,
    int d8,
    int d9
) : levels_(0),
    sort_importance_(0),
    sizecheck_traits_(0),
    nfxn_min_per_tile_(0),
    nidx_min_per_tile_(0)
{
    levels_.push_back(d1);
    if(d2 != -1) levels_.push_back(d2);
    if(d3 != -1) levels_.push_back(d3);
    if(d4 != -1) levels_.push_back(d4);
    if(d5 != -1) levels_.push_back(d5);
    if(d6 != -1) levels_.push_back(d6);
    if(d7 != -1) levels_.push_back(d7);
    if(d8 != -1) levels_.push_back(d8);
    if(d9 != -1) levels_.push_back(d9);

    build_sort_importance();

#if YETI_SANITY_CHECK || 1
    check_sort_importance();
#endif

    init();
}


PartitioningPolicy::~PartitioningPolicy()
{

}


void
PartitioningPolicy::init()
{
    min_fxn_at_level_ = new uli[levels_.size()];
    for(uli i = 0; i < levels_.size(); ++i) {
        min_fxn_at_level_[i] = 0;
    }
    max_fxn_per_tile_ = REALLY_BIG_INT;
}


void
PartitioningPolicy::build_sort_importance()
{
    for(int i = levels_.size() - 1; i >= 0; --i) {

        if(int(levels_[i]) / int(AngularMomentum) % 2 == 1)
            sort_importance_.push_back(AngularMomentum);

        if(int(levels_[i]) / int(DiffuseType) % 2 == 1)
            sort_importance_.push_back(DiffuseType);

        if(int(levels_[i]) / int(ModuloNumber) % 2 == 1)
            sort_importance_.push_back(ModuloNumber);

        if(int(levels_[i]) / int(GroupNumber) % 2 == 1)
            sort_importance_.push_back(GroupNumber);

        if(int(levels_[i]) / int(AtomNumber) % 2 == 1)
            sort_importance_.push_back(AtomNumber);

        if(int(levels_[i]) / int(SmallOrLargeAngularMomentum) % 2 == 1)
            sort_importance_.push_back(SmallOrLargeAngularMomentum);

    }
}


std::vector<PartitionableTrait>
PartitioningPolicy::split_traits(int item) {

    std::vector<PartitionableTrait> ret_val(0);

    if(item / int(AngularMomentum) % 2 == 1)
        ret_val.push_back(AngularMomentum);

    if(item / int(DiffuseType) % 2 == 1)
        ret_val.push_back(DiffuseType);

    if(item / int(ModuloNumber) % 2 == 1)
        ret_val.push_back(ModuloNumber);

    if(item / int(GroupNumber) % 2 == 1)
        ret_val.push_back(GroupNumber);

    if(item / int(AtomNumber) % 2 == 1)
        ret_val.push_back(AtomNumber);

    if(item / int(SmallOrLargeAngularMomentum) % 2 == 1)
        ret_val.push_back(SmallOrLargeAngularMomentum);

    return ret_val;
}


void
PartitioningPolicy::check_sort_importance()
{
    return;
    int idx = 0;
    for(int i = levels_.size() - 1; i >= 0; --i) {
        int traits_so_far = 0;
        // Bitwise AND checks to see if all of the traits in the
        // current level have been accounted for.
        while((traits_so_far & levels_[i]) != levels_[i] && idx < sort_importance_.size()) {
            traits_so_far |= sort_importance_[idx];
            ++idx;
        }

        if(idx >= sort_importance_.size())
            yeti_throw(SanityCheckError, "Invalid PartitioningPolicy sort importance.  Check your code and think about it...");

        int extra_sort = traits_so_far ^ levels_[i];  // Exclusive or to get the extra sort traits at this level
        for(int j = i-1; j >= 0; --j) {
            if((extra_sort & levels_[j]) != 0) // See if these traits exist at a lower level, which would mess things up
                yeti_throw(SanityCheckError, "Invalid PartitioningPolicy sort importance.  Check your code and think about it...");
        }
    }
}


bool
PartitioningPolicy::LessThanComparator::operator ()(
    const ShellPtr &sh1,
    const ShellPtr &sh2
)
{
    for(int i = 0; i < sort_importance_.size(); ++i) {
        switch(sort_importance_[i]) {
        case AngularMomentum:
        {
            int val1 = sh1->angular_momentum();
            int val2 = sh2->angular_momentum();
            if(val1 != val2) return val1 < val2;
            break;
        }
        case DiffuseType:
        {
            Shell::diffuse_type_t val1 = sh1->diffuse_type();
            Shell::diffuse_type_t val2 = sh2->diffuse_type();
            if(val1 != val2) return val1 < val2;
            break;
        }
        case ModuloNumber:
        {
            int val1 = sh1->get_atom()->modulo_number();
            int val2 = sh2->get_atom()->modulo_number();
            if(val1 != val2) return val1 < val2;
            break;
        }
        case GroupNumber:
        {
            int val1 = sh1->get_atom()->group_number();
            int val2 = sh2->get_atom()->group_number();
            if(val1 != val2) return val1 < val2;
            break;
        }
        case AtomNumber:
        {
            int val1 = sh1->get_atom()->number();
            int val2 = sh2->get_atom()->number();
            if(val1 != val2) return val1 < val2;
            break;
        }
        case SmallOrLargeAngularMomentum:
        {
            int val1 = sh1->angular_momentum() >= 2 ? 1 : 0;
            int val2 = sh2->angular_momentum() >= 2 ? 1 : 0;
            if(val1 != val2) return val1 < val2;
            break;
        }
        }
    }

    // If everything is the same, return false
    return false;
}


bool
PartitioningPolicy::should_partition_at_level(
    const ShellPtr &sh1,
    const ShellPtr &sh2,
    int level,
    uli numSoFar,
    std::vector<uli> fxnSoFar
)
{
    while(sizecheck_traits_.size() < levels_.size()) sizecheck_traits_.push_back(NoPartition);

    for(int i = level; i < levels_.size(); ++i) {

        if(int(sizecheck_traits_[i]) / int(AngularMomentum) % 2 == 1) {
            int val1 = sh1->angular_momentum();
            int val2 = sh2->angular_momentum();
            if(val1 != val2) {
                if(i == 0 && numSoFar >= nfxn_min_per_tile_) return true;
                if(i > 0 && (numSoFar >= nidx_min_per_tile_ && fxnSoFar[i] >= min_fxn_at_level_[i])) return true;
            }

        }
        else if(int(levels_[i]) / int(AngularMomentum) % 2 == 1) {
            int val1 = sh1->angular_momentum();
            int val2 = sh2->angular_momentum();
            if(val1 != val2) return true;
        }


        if(int(sizecheck_traits_[i]) / int(DiffuseType) % 2 == 1) {
            Shell::diffuse_type_t val1 = sh1->diffuse_type();
            Shell::diffuse_type_t val2 = sh2->diffuse_type();
            if(val1 != val2) {
                if(i == 0 && numSoFar >= nfxn_min_per_tile_) return true;
                if(i > 0 && (numSoFar >= nidx_min_per_tile_ && fxnSoFar[i] >= min_fxn_at_level_[i])) return true;
            }
        }
        else if(int(levels_[i]) / int(DiffuseType) % 2 == 1) {
            Shell::diffuse_type_t val1 = sh1->diffuse_type();
            Shell::diffuse_type_t val2 = sh2->diffuse_type();
            if(val1 != val2) return true;
        }



        if(int(sizecheck_traits_[i]) / int(ModuloNumber) % 2 == 1) {
            int val1 = sh1->get_atom()->modulo_number();
            int val2 = sh2->get_atom()->modulo_number();
            if(val1 != val2) {
                if(i == 0 && numSoFar >= nfxn_min_per_tile_) return true;
                if(i > 0 && (numSoFar >= nidx_min_per_tile_ && fxnSoFar[i] >= min_fxn_at_level_[i])) return true;
            }
        }
        else if(int(levels_[i]) / int(ModuloNumber) % 2 == 1) {
            int val1 = sh1->get_atom()->modulo_number();
            int val2 = sh2->get_atom()->modulo_number();
            if(val1 != val2) return true;
        }


        if(int(sizecheck_traits_[i]) / int(GroupNumber) % 2 == 1) {
            int val1 = sh1->get_atom()->group_number();
            int val2 = sh2->get_atom()->group_number();
            if(val1 != val2) {
                if(i == 0 && numSoFar >= nfxn_min_per_tile_) return true;
                if(i > 0 && (numSoFar >= nidx_min_per_tile_ && fxnSoFar[i] >= min_fxn_at_level_[i])) return true;
            }
        }
        else if(int(levels_[i]) / int(GroupNumber) % 2 == 1) {
            int val1 = sh1->get_atom()->group_number();
            int val2 = sh2->get_atom()->group_number();
            if(val1 != val2) return true;
        }


        if(int(sizecheck_traits_[i]) / int(AtomNumber) % 2 == 1) {
            int val1 = sh1->get_atom()->number();
            int val2 = sh2->get_atom()->number();
            if(val1 != val2) {
                if(i == 0 && numSoFar >= nfxn_min_per_tile_) return true;
                if(i > 0 && (numSoFar >= nidx_min_per_tile_ && fxnSoFar[i] >= min_fxn_at_level_[i])) return true;
            }
        }
        else if(int(levels_[i]) / int(AtomNumber) % 2 == 1) {
            int val1 = sh1->get_atom()->number();
            int val2 = sh2->get_atom()->number();
            if(val1 != val2) return true;
        }


        if(int(sizecheck_traits_[i]) / int(SmallOrLargeAngularMomentum) % 2 == 1) {
            int val1 = sh1->angular_momentum() >= 2 ? 1 : 0;
            int val2 = sh2->angular_momentum() >= 2 ? 1 : 0;
            if(val1 != val2) {
                if(i == 0 && numSoFar >= nfxn_min_per_tile_) return true;
                if(i > 0 && (numSoFar >= nidx_min_per_tile_ && fxnSoFar[i] >= min_fxn_at_level_[i])) return true;
            }
        }
        else if(int(levels_[i]) / int(SmallOrLargeAngularMomentum) % 2 == 1) {
            int val1 = sh1->angular_momentum() >= 2 ? 1 : 0;
            int val2 = sh2->angular_momentum() >= 2 ? 1 : 0;
            if(val1 != val2) return true;
        }


    }

    // No differences or the size isn't big enough, so don't insert a
    // partition between these to shells at this level
    return false;
}

void
PartitioningPolicy::sort_shells(std::vector<ShellPtr> &shells) const
{
    std::sort(shells.begin(), shells.end(), PartitioningPolicy::LessThanComparator(sort_importance_));
}


void
PartitioningPolicy::set_sizecheck_all(bool szchk)
{
    if(szchk) {
        sizecheck_traits_ = levels_;
    }
    else {
        sizecheck_traits_.clear();
        while(sizecheck_traits_.size() < levels_.size()) sizecheck_traits_.push_back(NoPartition);
    }
}

void
PartitioningPolicy::set_sizecheck(int trait)
{
    while(sizecheck_traits_.size() < levels_.size()) sizecheck_traits_.push_back(NoPartition);
    int i;
    for(i = 0; i < levels_.size(); ++i) {
        if((levels_[i] & trait) == trait) {
            sizecheck_traits_[i] |= trait;
            break;
        }
    }
    if(i == levels_.size()) {
        yeti_throw(SanityCheckError, "Can't set a policy to size-check a trait that is not part of its level hierarchy.");
    }
}


void
PartitioningPolicy::unset_sizecheck(int trait)
{
    while(sizecheck_traits_.size() < levels_.size()) sizecheck_traits_.push_back(NoPartition);
    int i;
    for(i = 0; i < levels_.size(); ++i) {
        if(levels_[i] & int(trait) == trait) {
            sizecheck_traits_[i] &= ~trait;
            break;
        }
    }
    if(i == levels_.size()) {
        yeti_throw(SanityCheckError, "Can't unset a policy to size-check a trait that is not part of its level hierarchy.");
    }
}


void
PartitioningPolicy::set_sort_priority(
    int item0,
    int item1,
    int item2,
    int item3,
    int item4,
    int item5,
    int item6,
    int item7,
    int item8,
    int item9
)
{
    sort_importance_.clear();

    int item = item0;
    std::vector<PartitionableTrait> items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    item = item1;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    item = item2;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    item = item3;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    item = item4;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    item = item5;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }


    item = item6;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }


    item = item7;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    item = item8;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    item = item9;
    if(item == -1) return;
    items.clear();
    items = split_traits(item);
    for(int i = 0; i < items.size(); ++i) {
        sort_importance_.push_back(items[i]);
    }

    check_sort_importance();
}


void
PartitioningPolicy::set_min_fxn_per_tile(int n)
{
    nfxn_min_per_tile_ = n;
    min_fxn_at_level_[0] = n;
    for(int i = 0; i < levels_.size(); ++i) {
        if(min_fxn_at_level_[i] < n) {
            min_fxn_at_level_[i] = n;
        }
    }
}


void
PartitioningPolicy::set_min_fxn_at_level(int level, uli minfxn)
{
    // TODO remove redundancy here and in the should_partition function
    min_fxn_at_level_[level] = minfxn;
    if(level == 0) {
        nfxn_min_per_tile_ = minfxn;
    }
}
