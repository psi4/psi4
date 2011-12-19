#ifndef yeti_MAPIMPL_H
#define yeti_MAPIMPL_H

#include "class.h"
#include "tuple.h"
#include "index.h"
#include "malloc.h"
#include "sortimpl.h"
#include "runtime.h"

namespace yeti {

class Indexer {

    protected:
        uli sizes_[NINDEX];

        uli offsets_[NINDEX];

        uli cumulative_sizes_[NINDEX];

        usi nindex_;

        usi depth_;

        uli size_;

    public:
        Indexer(
            const uli* sizes,
            const uli* offsets,
            usi nindex
        ) : nindex_(nindex),
            depth_(0),
            size_(1)
        {
            for (usi i=0; i < nindex_; ++i)
            {
                sizes_[i] = sizes[i];
                offsets_[i] = offsets[i];
            }

            init_sizes();
        }

        Indexer(
            TensorIndexDescr* tuple,
            const uli* indices = 0,
            usi depth = 0
        ) :
            nindex_(tuple->nindex()),
            size_(1),
            depth_(depth)
        {
            usi mapdepth = indices ? depth + 1 : tuple->depth();
            for (usi i=0; i < nindex_; ++i)
            {
                IndexDescr* descr = tuple->get(i);
                IndexRange* range = descr->get_range(mapdepth);
                if (indices)
                    range = range->get_subindex(indices[i]);
                sizes_[i] = range->nelements();
                offsets_[i] = range->start();
            }

            init_sizes();

        }

        void init_sizes()
        {
            size_ = 1;
            for (short i=nindex_ - 1; i >= 0; --i)
            {
                cumulative_sizes_[i] = size_;
                size_ *= sizes_[i];
            }

        }

        const uli* sizes() const {return sizes_;}

        const uli* cumulative_sizes() const {return cumulative_sizes_;}

        const uli* offsets() const {return offsets_;}

        uli size() const { return size_; }

        uli index(const uli* indices) const
        {
            uli idx = 0;
            for (usi i=0; i < nindex_; ++i)
            {
                idx += (indices[i] - offsets_[i]) * cumulative_sizes_[i];
            }
            return idx;
        }

        void indices(uli index, uli* indexset) const
        {
            uli remainder = index;
            for (usi i=0; i < nindex_; ++i)
            {
                uli ni = remainder / cumulative_sizes_[i];
                indexset[i] = ni + offsets_[i];
                remainder -= ni * cumulative_sizes_[i];
            }
        }

};

template <class T>
class IndexableMap :
    public Indexer
{

    private:
        bool mallocd_;

    protected:
        T** items_;

    public:
        typedef T** iterator;

        IndexableMap(
            TensorIndexDescr* tuple,
            const uli* indices = 0,
            usi depth = 0,
            MemoryPool* mempool = 0
        ) :
            Indexer(tuple, indices, depth),
            mallocd_(!mempool),
            items_(0)
        {

            if (mempool)
                items_ = reinterpret_cast<T**>(mempool->get(size_ * sizeof(T*)));
            else
                items_ = reinterpret_cast<T**>(YetiRuntime::malloc(size_ * sizeof(T*)));
            ::memset(items_, 0, size_ * sizeof(T*));
        }

        virtual ~IndexableMap()
        {
            if (mallocd_) YetiRuntime::free(items_, size_ * sizeof(T*));
        }

        /**
          @param p P is a sorting permutation that has been applied to the index.
                   this indexing operation "undoes" the effect of the permutation.
                   Thus indexset[pmap[i]] = ni + offsets[i].
        */
        static void permuted_indices(
            uli index,
            usi nindex,
            uli* indexset,
            const uli* sizes,
            const uli* offsets,
            Permutation* p
        )
        {
            uli remainder = index;
            const usi* indexmap = p->indexmap();
            for (usi i=0; i < nindex; ++i)
            {
                uli ni = remainder / sizes[i];
                usi imap = indexmap[i];
                indexset[imap] = ni + offsets[imap];
                remainder -= ni * sizes[i];
            }
        }

        T* get(const uli* indices) const
        {
            uli idx = index(indices);
#if YETI_SANITY_CHECK
            if (idx >= size_)
            {
                std::cerr << "Invalid index " << idx << " passed to map of size " << size_ << std::endl;
            }
#endif
            return items_[idx];
        }

        T* get(uli idx) const
        {
#if YETI_SANITY_CHECK
            if (idx >= size_)
            {
                std::cerr << "Invalid index " << idx << " passed to map of size " << size_ << std::endl;
            }
#endif
            return items_[idx];
        }

        iterator begin() const {return items_;}

        iterator end() const {return items_ + size_;}

        void set(const uli* indices, T* node)
        {
            uli idx = index(indices);
            items_[idx] = node;
        }

        void set(uli index, T* node)
        {
#if YETI_SANITY_CHECK
            if (items_[index])
            {
                std::cerr << "cannot set on already existing node location" << std::endl;
                abort();
            }
#endif
            items_[index] = node;
        }


        void sort(const SortPtr& sort)
        {
            sort->configure(sizes_);
            T** buffer = reinterpret_cast<T**>(sort->metadata_buffer(0));
            sort->sort_noscale<T*>(items_, buffer);
            ::memcpy(items_, buffer, size_ * sizeof(T*));
            //sort the nodes

            metadata_sort(sort);
        }

        void metadata_sort(const SortPtr& sort)
        {
            uli* tmp = reinterpret_cast<uli*>(sort->metadata_buffer(0));
            sort->get_permutation()->template permute<uli>(sizes_, tmp);
            ::memcpy(sizes_, tmp, nindex_ * sizeof(uli));

            sort->get_permutation()->template permute<uli>(offsets_, tmp);
            ::memcpy(offsets_, tmp, nindex_ * sizeof(uli));

            init_sizes();
        }

};

template <class T>
class BlockMap :
    public IndexableMap<T>
{

    protected:
        using IndexableMap<T>::size_;
        using IndexableMap<T>::items_;

    public:
        using IndexableMap<T>::begin;
        using IndexableMap<T>::end;

    public:
         BlockMap(
            TensorIndexDescr* tuple,
            const uli* indices = 0,
            usi depth = 0
        ) :
            IndexableMap<T>(
                tuple,
                indices,
                depth
            )
        {
        }

         void set(const uli* indices, T* node)
         {
             node->incref();
             IndexableMap<T>::set(indices, node);
         }

         void set(uli index, T* node)
         {
            node->incref();
            IndexableMap<T>::set(index, node);
         }

         ~BlockMap()
         {
             T** it = begin();
             T** stop = end();
             for ( ; it != stop; ++it)
             {
                 T* t = *it;
                 if (t)
                 {
                     boost::intrusive_ptr_release(t);
                 }
             }
         }

         void erase(uli index)
         {
             T* t = items_[index];
             if (t)
                 boost::intrusive_ptr_release(t);
             items_[index] = 0;
         }
};

template <class T>
class NodeMap :
    public IndexableMap<T>,
    public MempoolVirtualAddressMalloc
{

    protected:
        using IndexableMap<T>::size_;
        using IndexableMap<T>::items_;

    public:
        using IndexableMap<T>::begin;
        using IndexableMap<T>::end;

    public:
         NodeMap(
            MemoryPool* mempool,
            TensorIndexDescr* tuple,
            const uli* indices,
            usi depth
        ) :
            IndexableMap<T>(
                tuple,
                indices,
                depth,
                mempool
            )
        {
        }

        ~NodeMap()
        {
             //no need to delete anything
        }

        void realign_memory_pool(void* oldptr, void* newptr)
        {
            realign_pointer(items_, oldptr, newptr);
            T** it = begin();
            T** stop = end();
            for ( ; it != stop; ++it)
            {
                T* next = *it;
                if (next)
                {
                    realign_pointer(next, oldptr, newptr);
                    *it = next;
                }
            }
        }


};

}

#endif // MAPIMPL_H
