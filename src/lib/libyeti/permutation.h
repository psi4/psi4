#ifndef yeti_permutation_h
#define yeti_permutation_h

#include <vector>
#include <functional>

#include "class.h"
#include "yetiobject.h"

#include "tuple.h"
#include "permutation.hpp"
#include "index.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

#define NPERMUTATIONS 1500
#define NINDEX_PLUS_ONE 7

/**
  @class Permutation
  Encapsulates a permutation of indices
*/
class Permutation :
    public smartptr::Countable
{
    
    private:
        friend class PermutationTest;
        friend struct permutation_less;

        /** 
            The permutation mapping. pmap_[i] = j means the jth element gets mapped to index i in the permutaiton.
            e.g. The permutation {1 2 0} permutes the elements {1 2 3} to {1 3 2}
        */
        usi pmap_[NINDEX];

        /** The number of indices in the permutation */
        usi nindex_;

        /** The number of indices actually permuted */
        usi rank_;

        /** The order of the element, i.e. at what power does p^n = 1 */
        usi order_;

        /**
          The sign of the permutation, either plus or minus
        */
        short sign_;

        uli id_index_;

        void init();

        usi compute_rank(const usi* arr);

        usi compute_order();

        static void build_permutation(usi nindex, const usi* index_map);

        static void next_permutation_index(
            usi nindex,
            usi index_number,
            const usi* available_indices,
            usi* index_map
        );

        static void build_permutations(usi nindex);

        static void add_permutation(Permutation* p);

        static uli factorial_[NINDEX_PLUS_ONE];

        static Permutation* minus_perms_[NINDEX_PLUS_ONE][NPERMUTATIONS];

        static Permutation* plus_perms_[NINDEX_PLUS_ONE][NPERMUTATIONS];

        static uli inverses_[NINDEX_PLUS_ONE][NPERMUTATIONS];

        static uli products_[NINDEX_PLUS_ONE][NPERMUTATIONS][NPERMUTATIONS];

        static Permutation* identities_[NINDEX_PLUS_ONE];

        static uli id_index(usi nindex, const usi* indexmap);

    public:
        /** How to expand a permutation. See #expand */
        typedef enum { ExpandForward, ExpandBackward } expand_t;
        
        typedef enum { ValidPermutationMalloc } permutation_malloc_flag_t;


        /**
            @param n      The number of items in the complete permutation group
            @param sign   The sign of the permutation for antisymmetric or symmetric quantities
            @param cycle An integer sequence defining the permutation.  See #pmap_;
        */
        Permutation(
            usi n,
            short sign,
            const usi* cycle
        );

        virtual ~Permutation();

        static void init_permutations();

        static Permutation* get_permutation(short sign, usi nindex, uli id_index);

        static Permutation* get_permutation(short sign, usi nindex, const usi* indexmap);

        static Permutation* get_permutation(
            short sign,
            usi nindex,
            usi nsubindex,
            const usi* indexmap,
            const usi* subset
       );

        static Permutation* get_permutation(
                short sign,
                usi nsubindex,
                usi nindex,
                const usi* indexmap,
                expand_t type
        );


        static Permutation* get_identity(usi nindex);

        /**
            Template method for permuting indices
            @param indices 
        */
        template <class T>
        T* permute(const T* indices) const;

        /**
            Template method for permuting indices
            @param indices
        */
        template <class T>
        void permute(const T* src, T* dest) const;

        /**
            @param indices map the given set of indices based on the permutation
        */
        uli* permute(const uli* indices) const;

        void permute(const uli* src, uli* dest) const;

        void image(const usi* src, usi* img, usi nindex) const;

        bool image_equals(const uli* src, const uli* img) const;

        void validate();

        bool is_inverse(Permutation* p) const;

        /**
            @param indices map the given set of indices based on the permutation
        */
        int* permute(const int* indices) const;

        uli id_index() const;

        /**
            @param indices map the given set of indices based on the permutation
        */
        usi* permute(const usi* indices) const;

        /**
            @return The order of the permutation group this belongs to
        */
        usi nindex() const;

        /**
            @return The number of indices that are actually permuted
        */
        usi rank() const;

        /**
            @return The power of the element that brings it to unity
        */
        usi order() const;

        /**
         * @param index
         * @return Whether the given index is shifted by the permutation
         */
        bool permutes(usi index) const;

        /**
         * @param index
         * @return Whether the given index is unaffected by the permutation
         */
        bool fixes(usi index) const;

        usi image(usi index) const;

        /**
            Determines whether the permutation leaves the values unchanged. For example,
            given the values (0 0 1 2) the permutation (1 0 2 3) leaves the
            values unchanged.  The permutation (0 1 3 2) would have the image
            (0 0 2 1), making the values different.
            @return
        */
        bool fixes_arr(const uli* vals) const;

        template <typename data_t>
        bool fixes_arr(const data_t* vals) const;

        bool lt(Permutation* p, const uli* vals) const;

        /**
            @return The sign of the permutation, either symmetric or anti-symmetric
        */
        short sign() const;

        const usi* indexmap() const;

        /**
            Returns the product permutation correponding to this * p
            If p is { 1 0 2 } and this is { 2 0 1 } the product is
            { 2 1 0 }.  The order here is very important.  p * this
            would be { 0 2 1 }.
            @param p The right permutation in the product
            @return The product permutation
        */
        Permutation* product(Permutation* p) const;

        Permutation* inverse() const;

        /**
            No error check is done on subset or n.
            Returns whether this permutation is valid for a 
            permutational subgroup which only permutes the subset. For example,
            consder the permutations {0 1 3 2}.  If subset is {1 2 3}, then the
            permutation is valid because it "acts" on indices 2, 3 and they are included
            in the subset.  In constrast, if the subset is {1 2} then the permutation is
            not valid because the permutation acts on index 3, but it is not included in the subset.
            @param subset
            @param n  The number of elements in the subset
            @return Validity of subset
        */
        bool valid_subgroup(const usi* subset, usi n) const;

        /**
            @return Whether a permutation of subindices is equivalent if expanded
        */
        bool is_equivalent(Permutation* p, expand_t type) const;

        /**
            Determines whether the indices are permutationally improvable.
            For example, given the transposition (1 2), the index set {0 3 2}
            would be improved because it maps to {0 2 3}, which is "less than."
            Given the transposition (0 1), the index set {0 3 2} would not be
            improved because it maps to {3 0 2} which is "greater than"
            @param indexset
            @return If the permutation would improve the sort of the indices
        */
        bool improves_sort(const uli* indices) const;

        void lowest_image(const uli* src, uli* dst) const;

        /**
          Compress the permutation to give a permutation for
          the n indices defined by subset
          @param subset The indices to pick out for the new permutation
          @param n The number of subset indices
          */
        Permutation* compress(const usi* subset, usi n) const;

        /**
          Whether the two permutations are equal. This requires the
          permutations to be of both the same length.
          @param p
          @return Whether the permutations are equal
          */
        bool equals(Permutation* p) const;

        /**
          Whether the permutation is a transposition of i,j
          and only a transposition. No other indices move.
          @param i The first index in the transpostion
          @param j The second index in the transposition
          @return Whether the permutation is a transposition of indices i,j
          */
        bool is_transposition(usi i, usi j) const;

        void print(std::ostream& os = std::cout) const;

        bool is_identity() const;

        /**
          Expands the permutation to a larger number of elements. The
          expanded indices are fixed in the new permutation.  For example,
          the permutation (1 0) could be expanded to 4 elements, giving the
          permutation (1 0 2 3).
          @param n The number of indices in the new permutation
          @param type Whether the padded indices are appended to the
                        end or beginning of the permutation
          @return The expanded permutation
         */
        Permutation* expand(usi n, expand_t type) const;

        bool maps_to(const uli* src, const uli* dst) const;

        void* operator new(size_t size, permutation_malloc_flag_t flag);

        void operator delete(void* ptr, permutation_malloc_flag_t flag);

        void* operator new(size_t size);

        void operator delete(void* ptr);
};

/** @class PermutationGroup
    Encapsulates a group of permutations
*/
class PermutationSet :
        public smartptr::Countable
{

    public:
        typedef std::vector<Permutation*>::const_iterator iterator;

    protected:
        void sanity_check(const PermutationSetPtr& grp) const;

        bool closed_;

    protected:
        /** The set of permutations in the group */
        std::vector<Permutation*> perms_;

        /** The number of indices each permutation in the group permutes */
        usi nindex_;

        /**
            Order the permutations such that the lower rank permutations
            come first.
        */
        void sort();

        void print_elements(
            std::ostream& os,
            const std::string& title
        );

    public:
        /**
            @param n The order of the permutation group
        */
        PermutationSet(usi n);

        PermutationSet(const PermutationSetPtr& set);

        virtual ~PermutationSet();

        /**
            @return The number of elements the permutation group is meant to permute
        */
        usi nindex() const;

        /**
            @return The number of permutations in the group
        */
        uli order() const;

        /**
            Adds a permutation to the group. Duplicates are skipped.
            @param p
            @throw ValuesNotEqual if the permutation corresponds to a different number of elements
        */
        void add(Permutation* p);

        void add(const PermutationSetPtr& pgrp);

        /**
            Add a permutation to the group which permutes less than #n_ indices
            and expand it.  See #_Permutation::expand. Duplicates are not added.
            @param pgrp The set of permutations to add
            @param type Whether to append indices to the front or back
        */
        void add_and_expand(const PermutationSetPtr& pgrp,
                             Permutation::expand_t type);

        /**
            Check if a permutation if already contained in the group. This does not throw if p
            correponds to a different number of n elements.
            @param p
        */
        bool contains(Permutation* p) const;

        /**
            Check if a permutation if already contained in the group. This does not throw if p
            correponds to a different number of n elements.
            @param p
        */
        bool contains(const PermutationSetPtr& p) const;

        bool equals(const PermutationSetPtr& set) const;

        Permutation* get(uli index) const;

        /**
            @param pgrp
            @return A permutation set containing the permutations in both this and pgrp
            @throw ValuesNotEqual if the groups correspons to a different number of n elements
        */
        PermutationSetPtr intersection_set(const PermutationSetPtr& pgrp) const;

        Permutation* lowest_image(const uli* src, uli* dst) const;

        static Permutation* lowest_image(Permutation** perm_list, usi nperms, const uli* src, uli* dst, usi nindex);

        /**
            @param pset
            @return A permutation set containing the permutations in either this or pset
            @throw ValuesNotEqual if the groups correspons to a different number of n elements
        */
        PermutationSetPtr union_set(const PermutationSetPtr& pgrp) const;

        PermutationSetPtr orbit(Permutation* p) const;

        bool is_product(Permutation* p) const;

        /**
            Builds a permutation group where the permutations involve only the subset
            @param subset
            @return The permutations involving the subset
        */
        PermutationSetPtr compressed_subset(const usi* subset, usi n) const;

        /**
            Builds a permutation group where the permutations involve only the subset
            @param subset
            @return The permutations involving the subset
        */
        PermutationSetPtr subset(const usi* subset, usi n) const;

        PermutationSetPtr complement_set(const PermutationSetPtr& set) const;

        /**
         * The set of permutations that fixes the vals
         * @param vals the values
         * @return the group
         */
        PermutationSetPtr isotropy_set(const uli* vals) const;
        
        bool nontrivial_isotropy(const uli* indices) const;

        /**
            Whether the group has the transposition of i and j
            @param i
            @param j
            @return If group contains (i j)
        */
        bool has_transposition(usi i, usi j) const;

        iterator begin() const;

        iterator end() const;

        void print(std::ostream& os = std::cout) const;

        /**
            Whether the set contains any permutations the improve the sort.
            See Permutation::improves_sort.
            @param indices
            @return If the indices can be sorted by the permutations set
        */
        virtual bool improves_sort(const uli* indices) const;

        /**
            Removes all elements form the group
        */
        void clear();

        Permutation* get_lowest_permutation(const uli* vals) const;
        
        PermutationSetPtr get_generator_set() const;

        /**
            Given a set of generating permutations, create all possible
            products, making a "closed" group
        */
        void close();

        bool closed() const;

        void set_closed();


};


class PermutationGroup :
        public PermutationSet
{

    private:
        void sanity_check() const;

        void sanity_check(const PermutationGroupPtr& grp) const;

    public:
        PermutationGroup(usi n);

        PermutationGroup(const PermutationSetPtr& p);

        PermutationGroup(
            usi n,
            const PermutationGroupPtr& grp,
            const usi* subset
        );

        virtual ~PermutationGroup();


        /**
            @param pgrp
            @return A permutation set containing the permutations in both this and pgrp
            @throw ValuesNotEqual if the groups correspons to a different number of n elements
        */
        PermutationGroupPtr intersection_grp(const PermutationGroupPtr& pgrp) const;

        /**
            @param pset
            @return A permutation set containing the permutations in either this or pset
            @throw ValuesNotEqual if the groups correspons to a different number of n elements
        */
        PermutationGroupPtr union_grp(const PermutationSetPtr& pset) const;

        /**
            Builds a permutation group where the permutations involve only the subset
            @param subset
            @return The permutations involving the subset
        */
        PermutationGroupPtr compressed_subgrp(const usi* subset, usi n) const;

        /**
            Builds a permutation group where the permutations involve only the subset
            @param subset
            @return The permutations involving the subset
        */
        PermutationGroupPtr subgrp(const usi* subset, usi n) const;

        PermutationGroupPtr conjugate(Permutation* p) const;

        /**
         * The set of permutations that fixes the vals
         * @param vals the values
         * @return the group
         */
        PermutationGroupPtr isotropy_grp(const uli* vals) const;

        template <typename data_t>
        PermutationGroupPtr isotropy_grp(const data_t* indices);

        bool improves_sort(const uli* indices) const;

        Permutation* get_identity() const;

        /**
            Builds a group of "unrelated" permutations where the relation is based
            on the denominator group pgrp. If we have two elements e,f in the current
            permutation group and there exists p in pgrp such that e*p = f,
            then f is a non-unique element. Only e will be placed in the new permutation group.
            @param pgrp
            @return The group of unique elements based on quotient relation
        */
        PermutationSetPtr quotient_set(const PermutationGroupPtr& pgrp) const;

};

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
