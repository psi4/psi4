#ifndef yeti_sort_h
#define yeti_sort_h

#include "class.h"

#include "permutation.hpp"
#include "sort.hpp"
#include "index.hpp"

#include "mallocimpl.h"
#include "yetiobject.h"

namespace yeti {

class Sort :
    public YetiRuntimeSerializable,
    public Malloc<Sort>
{

    private:
        friend class SortTest;
        friend class ThreadedSort;

        /**
         * The permutation generating the sort
         */
        PermutationPtr p_;

        /** Array of length #nindex_ giving the length of each stride */
        uli* lengths_;

        /** Array of length #nindex_ giving the number of strides on each index */
        uli* nstrides_;

        /** The total number of values in target/src arrays*/
        uli ntot_;

        /** The number of indices in the sort */
        usi nindex_;

        Sort(const PermutationPtr& p);

    public:
        /**
         * Constructor
         * @param p The permutation defining the sort
         * @param sizes Array of size p->n() giving the sizes
         * @param stepsize The number of entries to skip on each step.  This is mostly
                            intended for sorting, e.g. an array of doubles treated as
                            char* in which case the stepsize is sizeof(double), usually 8.
                            This should always be 1 except in rare cases.
         */
        Sort(
            const PermutationPtr& p,
            const size_t* sizes
        );

        Sort(
            const PermutationPtr& p,
            const IndexRangeTuplePtr& tuple
        );

        virtual ~Sort();

        /*
         * @return The number of indices in the sort AFTER compression
        */
        usi n() const;

        uli ntot() const;

        /*
         * @return The permutation object defining the sort
         */
        Permutation* get_permutation() const;

        void configure(const uli* sizes);

        void configure(const IndexRangeTuplePtr& tuple);

        void configure(const PermutationPtr& p, const IndexRangeTuplePtr& tuple);

        template <class T>
        T* sort(const T* src) const;

         /**
         * Allocate new array and sort values
         * @param src The source array
         * @return Sorted values
         */
        template <class T, class U>
        void sort(const T* src, U* dst) const;

        template <class T, class U>
        void sort_noscale(const T* src, U* dst) const;

        void print(std::ostream &os = std::cout) const;

         /**
         * Allocate new array and sort values
         * @param src The source array
         * @return Sorted values
         */
        template <class T, class U>
        void accumulate(const T* src, U* dst, U scale = 0) const;

    private:
        template <class T, class U>
        void stride_sort(usi sortlevel, const T* src, U*& target) const;

        template <class T, class U>
        void stride_accumulate(usi sortlevel, const T* src, U*& target, U scale) const;


};


class ThreadedSort :
    public YetiRuntimeCountable
{

    private:
        std::vector<SortPtr> sorters_;

        PermutationPtr perm_;

    public:
        ThreadedSort(const PermutationPtr& p);

        Sort* get_sorter(uli threadnum) const;

        Permutation* get_permutation() const;
};


}

#endif

