#ifndef yeti_aobasis_h
#define yeti_aobasis_h

#include "class.h"

#include "index.hpp"
#include "aobasis.hpp"

#include "gigmatrix.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class Shell :
    public smartptr::Countable
{

    public:
        typedef enum { valence, diffuse } diffuse_type_t;

        typedef enum { high_am, low_am} am_type_t;

    private:
        AtomPtr atom_;

        usi ang_mom_;

        uli intra_am_index_;

        uli nfxn_;

        uli nsets_;

        uli program_shell_number_;

        uli program_fxn_start_;

        double min_exponent_;

        diffuse_type_t diffuse_type_;

        am_type_t am_type_;

        static bool statics_done_;

        void init_statics();

        static std::map<std::string, std::map<usi, double> > diffuse_cutoffs_;

    public:
        /**
          @param atom
          @param ang_mom
          @param nsets  The number of sets of functions.  For example, you could have a DZ basis with
                        with 2 sets of p-functions combined in a single shell
          @param nfxn   The total number of basis functions.  For example, you could have a DZ basis
                        with 2 sets of p-functions or 6 total functions
          @param min_exponent   The minimum Gaussian exponent (or minimum ``average'' exponent if
                                a contracted function) for all sets of functions in the shell
          @param program_shell_number   The shell number as stored by the interfacing integral program
        */
        Shell(
            const AtomPtr& atom,
            usi ang_mom,
            uli nsets,
            uli nfxn,
            double min_exponent,
            uli program_shell_number,
            uli program_fxn_start,
            uli intra_am_index
        );


        diffuse_type_t diffuse_type() const;

        double min_exponent() const;

        uli intra_am_index() const;

        usi angular_momentum() const;

        am_type_t am_type() const;

        uli nfxn() const;

        const AtomPtr& get_atom() const;

        uli program_shell_number() const;

        uli program_fxn_start() const;

        void print(std::ostream& os = std::cout) const;

};

class Atom :
    public smartptr::Countable
{

    private:
        std::string symbol_;

        std::string name_;

        double x_;

        double y_;

        double z_;

        uli number_;

        uli group_number_;

        uli modulo_number_;

    public:

        Atom(
            const std::string& symbol,
            const std::string& name,
            double x,
            double y,
            double z,
            uli atom_number
        );

        double x() const;

        double y() const;

        double z() const;

        const std::string& symbol() const;

        const std::string& name() const;

        uli number() const;

        uli group_number() const;

        uli modulo_number() const;

        void set_group_number(uli num);
        
        void set_modulo_number(uli num);

};

class AOBasis :
    public smartptr::Countable
{

    private:
        std::vector<ShellPtr> shells_;

        std::vector<AtomPtr> atoms_;

        uli nbasis_;

        uli nfxn_min_per_tile_;

        uli* yeti_to_program_shell_map_;

        uli* yeti_to_program_fxn_map_;

        uli* program_to_yeti_shell_map_;

        uli* program_to_yeti_fxn_map_;

        std::string name_;

        std::string descr_;

        IndexRangePtr range_;

        IndexDescrPtr index_descr_;

        MultiShellMapPtr multi_shell_map_;

        void configure_atom_groups_as_atoms();

        void configure_atom_groups_by_number(uli natoms_per_group);

        void configure_many_atoms_large_basis(uli natoms_per_group);

    public:
        AOBasis(const std::string& name, const std::string& descr);

        void add_shell(const ShellPtr& shell);

        void add_atom(const AtomPtr& atom);

        void configure(PartitioningPolicyPtr policy);

        DEPRECATED void configure_many_atoms_small_basis();

        DEPRECATED void configure_few_atoms_large_basis();

        DEPRECATED void configure_many_atoms_large_basis();

        DEPRECATED void configure_few_atoms_small_basis();

        void configure_min_nfxn_per_tile(uli nmin);

        void convert_to_yeti_numbering(RectMatrixPtr& coefs);

        IndexRange* get_index_range() const;

        IndexDescr* get_index_descr() const;

        const MultiShellMapPtr& get_multishell_map() const;

        const std::string& name() const;

        const std::string& descr() const;

        uli program_shell_number(uli yeti_shell_number);

        uli program_fxn_number(uli yeti_fxn_number);

        uli yeti_shell_number(uli program_shell_number);

        uli yeti_fxn_number(uli program_fxn_number);

        uli nbasis() const;
};

class MultiShellMap :
    public smartptr::Countable
{

    private:
        uli* shell_sizes_;

        uli* shell_start_on_index_;

        uli* nshells_per_index_;

        uli* nfxn_per_index_;

    public:
        MultiShellMap(
            const std::vector<ShellPtr>& shells,
            const IndexDescrPtr& ao_descr
        );

        uli shell_size(uli shell);

        uli shell_start(uli idx);

        uli nshell(uli idx);

        uli nfxn(uli idx);

};

enum PartitionableTrait {
    NoPartition = 0, // Not necessarily no partitioning.  It still adopts all of the higher levels' partitioning traits
    AngularMomentum = 1,
    DiffuseType = 2,
    ModuloNumber = 4,
    GroupNumber = 8,
    AtomNumber = 16,
    SmallOrLargeAngularMomentum = 32  // Two possible values:  0 for s and p, 1 for d, f, g, etc.
};

class PartitioningPolicy :
    public smartptr::Countable
{
    protected:

        /**
          Describes the various levels of partitioning.
          levels_[0] corresponds to depth = 1 (so the tree
          structure is inverted, as in the rest of the Yeti suite)
          */
        std::vector<int> levels_;

        /**
          This defaults to the same as levels_, but offers some finer-grained
          control in cases where the user wishes to combine traits in one level
          but still specify an order of importance in sorting that level.  Note that
          these values cannot be bitwise OR'd
          */
        std::vector<PartitionableTrait> sort_importance_;


        /**
          This is a list of traits (which can be bitwise OR'd) at each level that
          only cause partitioning to occur if there
          are at least nfxn_min_per_tile_ functions in the tile (at depth = 1) or
          nidx_min_per_tile indices in a metadata tile (at depth > 1)
          */
        std::vector<int> sizecheck_traits_;

        /**
          After the splitting [see AOBasis::configure(PartitioningPolicyPtr)] is completed,
          ranges at each level are checked to see if they contain the minimum number of
          functions specified for that level.  If not, yeti will \em try to combine ranges at that level
          to get above the minimum.  However, it will not combine ranges with different parents
          unless the parents' level's minimum function value is less that the functions (not indices or shells,
          basis functions) in that level, in which case the neighboring parents will be merged and
          the merge will be propagated down.  The merge will also only take place if the merged
          data tile's size is less than max_fxn_per_tile_ (or if max_fxn_per_tile_ is < 1, which
          turns off that particular feature)
          */
        uli* min_fxn_at_level_;

        /**
          The maximum number of functions per data tile.  After the splitting [see AOBasis::configure(PartitioningPolicyPtr)] is completed,
          ranges at the data level (level = 0, a.k.a. depth = 1) are checked to see if they contain more than the maximum number
          of basis functions specified.  If so, and if splitting would not cause the number to fall below min_fxn_at_level_[0],
          the data block is split into blocks as evenly sized as possible (with as few splits as possible to get the number
          below the max_fxn_per_tile_ value.
          */
        uli max_fxn_per_tile_;

        /**
          This value affects the partitioning or lack there of with respect to the sizecheck_traits_ vector
          at the data level (level = 0, a.k.a. depth = 1).  At this point it is exactly equivalent to
          min_fxn_at_level_[0].
          */
        int nfxn_min_per_tile_;

        /**
          The same concept as PartitioningPolicy::nfxn_min_per_tile_, except this counts indices at levels
          higher than the data level.
          */
        int nidx_min_per_tile_;

        void build_sort_importance();

        void check_sort_importance();

        class LessThanComparator {
            private:
                std::vector<PartitionableTrait> sort_importance_;

            public:

                LessThanComparator(std::vector<PartitionableTrait> si) : sort_importance_(si) { }

                bool operator() (const ShellPtr& sh1, const ShellPtr& sh2);
        };

        static std::vector<PartitionableTrait> split_traits(int item);

        void init();

    public:

        /// A default constructor that uses some reasonable defaults
        PartitioningPolicy();

        PartitioningPolicy(
            std::vector<int> levels
        );

        PartitioningPolicy(
            int d1,
            int d2 = -1,
            int d3 = -1,
            int d4 = -1,
            int d5 = -1,
            int d6 = -1,
            int d7 = -1,
            int d8 = -1,
            int d9 = -1
        );

        ~PartitioningPolicy();

        int nlevels() const { return levels_.size(); }

        /// This function takes into account sizecheck_traits_
        bool should_partition_at_level(const ShellPtr &sh1, const ShellPtr &sh2, int level, uli numSoFar, std::vector<uli> fxnSoFar);

        /**
          Sorts the shells according to the specified sort priority associated with the PartitioningPolicy instance
          for which the method is called.  See also PartitioningPolicy::set_sort_priority and PartitioningPolicy::sort_importance_
          */
        void sort_shells(std::vector<ShellPtr> &shells) const;

        /** Setter function for PartitioningPolicy::nfxn_min_per_tile_.
          This function also sets min_fxn_at_level_[0] and min_fxn_at_level[i] if the previously specified value
          (or if unspecified) is less than the parameter given.
          */
        void set_min_fxn_per_tile(int n);

        /// Getter function for PartitioningPolicy::nfxn_min_per_tile_.
        int min_fxn_per_tile() const { return nfxn_min_per_tile_; }

        /// Setter function for PartitioningPolicy::nidx_min_per_tile_
        void set_min_idx_per_tile(int n) { nidx_min_per_tile_ = n; }

        /// Getter function for PartitioningPolicy::nidx_min_per_tile_
        int min_idx_per_tile() const { return nidx_min_per_tile_; }

        uli min_fxn_at_level(int level) const { return min_fxn_at_level_[level]; }

        void set_min_fxn_at_level(int level, uli minfxn);

        uli max_fxn_per_tile() const { return max_fxn_per_tile_; }

        void set_max_fxn_per_tile(uli n) { max_fxn_per_tile_ = n; }

        /**
          If the parameter given is true or no parameter is given, all traits are set to be size-checked.
          If the parameter given is false, no parameters are size-checked.
          Either way, any previous size-checked settings are obliterated.
          */
        void set_sizecheck_all(bool szchk = true);

        /**
          Sets the trait to be size-checked.  Note that you can use bitwise OR'd values here, but only if
          the traits are on the same level.  Make multiple to set multiple levels.
          */
        void set_sizecheck(int trait);

        /**
          Similar to set_sizecheck
          */
        void unset_sizecheck(int trait);

        /**
          See set_sort_priority(int, int, int, int, int, int, int, int, int, int)
          */
        void set_sort_priority(std::vector<PartitionableTrait> sort_importance) {
            sort_importance_ = sort_importance;
            check_sort_importance();
        }

        /**
          Set the contents of thesort priority of the traits.  Traits not in the levels_ vector may be
          included and bitwise OR'd values are allowed, however the sort priority must at least include all
          traits from the levels_ vector, with the traits from higher levels coming <i>before</i> the traits
          from lower levels (failure to do so will result in a runtime error).  If a bitwise OR'd value is given,
          the order of the OR'd values is not specified (though it is deterministic; see source
          code) and should not be relied upon (i. e. this ordering may change in the future).
          Only give bitwise OR'd values if it really does not matter what order the sorting is
          carried out in relative to these traits.
          The sort priority will cause the PartitioningPolicy::LessThanComparator functor to first compare
          based on item0, then item1, then item2, and so on (up to 9 traits may be specified here, which is
          more than the number of traits currently available in the enumerated type PartitionableTrait, so
          you should be okay; if you need more for whatever reason, there is also a version that takes a
          std::vector of PartitionableTrait values).
          */
        void set_sort_priority(
            int item0,
            int item1 = -1,
            int item2 = -1,
            int item3 = -1,
            int item4 = -1,
            int item5 = -1,
            int item6 = -1,
            int item7 = -1,
            int item8 = -1,
            int item9 = -1
        );
};

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif


