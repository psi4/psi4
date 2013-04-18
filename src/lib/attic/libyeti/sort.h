#ifndef yeti_sort_h
#define yeti_sort_h

#include "class.h"

#include "permutation.hpp"
#include "sort.hpp"
#include "index.hpp"

#include "mallocimpl.h"
#include "yetiobject.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

#define SORT_BUFFER_SIZE 5000000 //5 MB
class Sort :
    public smartptr::Countable,
    public Malloc<Sort>
{

    private:
        friend class SortTest;
        friend class ThreadedSort;

        /**
         * The permutation generating the sort
         */
        Permutation* p_;

        Permutation* inv_p_;

        /** Array of length #nindex_ giving the length of each stride */
        uli lengths_[NINDEX];

        /** Array of length #nindex_ giving the number of strides on each index */
        uli nstrides_[NINDEX];

        char metadata_buffer_[SORT_BUFFER_SIZE];

        char data_buffer_[SORT_BUFFER_SIZE];

        /** The total number of values in target/src arrays*/
        uli ntot_;

        /** The number of indices in the sort */
        usi nindex_;

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
        Sort(Permutation* p);

        virtual ~Sort();

        /*
         * @return The number of indices in the sort AFTER compression
        */
        usi nindex() const;

        const uli* lengths() const;

        const uli* nstrides() const;

        uli ntot() const;

        uli nelements() const;

        uli total_stride() const;

        char* data_buffer() const;

        char* metadata_buffer(usi depth) const;

        /*
         * @return The permutation object defining the sort
         */
        Permutation* get_permutation() const;

        Permutation* get_inverse_permutation() const;

        void configure(const uli* sizes);
        
        void configure(Permutation* p);

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

        template <class T, class U, class Functor>
        void sort_functor(
            T* src,
            U* dst,
            Functor& functor
        ) const;

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

        template <class T, class U, class Functor>
        void stride_functor(
            usi sortlevel,
            T* src,
            U*& target,
            Functor& functor
        ) const;

        template <class T, class U>
        void stride_accumulate(usi sortlevel, const T* src, U*& target, U scale) const;


};



}

#ifdef redefine_size_t
#undef size_t
#endif

#endif

