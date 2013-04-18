#ifndef yeti_tuple_h
#define yeti_tuple_h

#include <vector>
#include "class.h"
#include "mallocimpl.h"
#include "yetiobject.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define foreach_nonnull(item, container, cls, ...) \
{ \
cls** item##_it(container->begin()); \
cls** item##_stop(container->end()); \
cls* item; \
for ( ; item##_it != item##_stop; ++item##_it) \
{ \
    item = *item##_it; \
    if (item) \
    { \
        __VA_ARGS__ \
    } \
 } \
}

#define foreach(item, container, cls, ...) \
{ \
cls** item##_it(container->begin()); \
cls** item##_stop(container->end()); \
cls* item; \
for ( ; item##_it != item##_stop; ++item##_it) \
{ \
    item = *item##_it; \
    __VA_ARGS__ \
 } \
}

namespace yeti {


template <class T, uli size_>
class StorageArray
{

    private:
        uli n_;

        T* items_[size_];

    public:
        StorageArray()
            :
            n_(0)
        {
            ::memset(items_, 0, size_ * sizeof(T*));
        }

        ~StorageArray()
        {
            clear();
        }

        void
        append(T* item)
        {
            if (n_ >= size_)
            {
                std::cerr << "Insufficient space to append in storage array" << std::endl;
                abort();
            }
            items_[n_] = item;
            ++n_;
        }

        T*
        get(uli idx)
        {
            return items_[idx];
        }

        T* last()
        {
            return n_ == 0 ? 0 : items_[n_ - 1];
        }

        void clear()
        {
        }

        bool full() const
        {
            return n_ == size_;
        }

        uli nelements() const
        {
            return n_;
        }

        T** begin() const
        {
            return const_cast<T**>(items_);
        }

        T** end() const
        {
            T** stop = const_cast<T**>(items_);
            return stop + n_;
        }

        void realign_memory_pool(void* oldptr, void* newptr)
        {
            T** it(items_);
            T** end(items_ + size_);
            for ( ; it != end; ++it)
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


template <
    class T
>
class CountableArray :
    public smartptr::Countable
{

    protected:
        T** items_;

        uli size_;

        bool mallocd_;

        CountableArray(
            uli size,
            T** items
        ) : size_(size),
            items_(items),
            mallocd_(false)
        {
            init();
        }

    public:
        typedef T** iterator;

        CountableArray(
            uli size
        )
            :
            size_(size),
            items_(new T*[size]),
            mallocd_(true)
        {
            init();
        }

        void
        init()
        {
            for (unsigned long i=0; i < size_; ++i)
                items_[i] = 0;
        }

        ~CountableArray()
        {
            clear();
            if (mallocd_)
                delete[] items_;
        }

        T* get(uli n) const
        {
            return items_[n];
        }

        void erase(uli index)
        {
            T* ptr = items_[index];
            items_[index] = 0;
            if (ptr)
                boost::intrusive_ptr_release(ptr);
        }

        void insert(uli n, T* item)
        {
            items_[n] = item;
            item->incref();
        }

        void set(uli n, T* item)
        {
            insert(n, item);
        }

        T** begin() const
        {
            return items_;
        }

        T** end() const
        {
            return items_ + size_;
        }

        void clear()
        {
            T** it(items_);
            T** end(items_ + size_);
            for ( ; it != end; ++it)
            {
                T* next(*it);
                if (next)
                    boost::intrusive_ptr_release(next);
                *it = 0;
            }
        }

        uli n_nonnull() const
        {
            uli ntot = 0;
            T** it(items_);
            T** end(items_ + size_);
            for ( ; it != end; ++it)
            {
                T* next(*it);
                if (next)
                    ++ntot;
            }
            return ntot;
        }

        uli size() const
        {
            return size_;
        }


};

template <class T>
class MetaDataArray :
    public CountableArray<T>,
    public MempoolVirtualAddressMalloc
{
    protected:
        using CountableArray<T>::items_;
        using CountableArray<T>::size_;

    private:
        uli n_;

    public:
        MetaDataArray(
            MemoryPool* mem,
            uli size
        )
            :
            CountableArray<T>(
                size,
                reinterpret_cast<T**>(
                    mem->get(size * sizeof(T*))
                )
            ),
            n_(0)
        {
        }

        void realign_memory_pool(void* oldptr, void* newptr)
        {
            realign_pointer(items_, oldptr, newptr);
            T** it(items_);
            T** end(items_ + n_);
            for ( ; it != end; ++it)
            {
                T* next = *it;
                realign_pointer(next, oldptr, newptr);
                *it = next;
            }
        }

        void
        append(T* item)
        {
            if (n_ >= size_)
            {
                std::cerr << "Insufficient space to append in metadata array" << std::endl;
                abort();
            }
            CountableArray<T>::set(n_, item);
            ++n_;
        }

        T* pop_back()
        {
            return items_[n_-1];
        }

        uli nelements() const
        {
            return n_;
        }

        T** begin() const
        {
            return items_;
        }

        T** end() const
        {
            return items_ + n_;
        }

        /**
            Clear all data, but do not attempt to inc/dec ref count.
        */
        void memset()
        {
            ::memset(items_, 0, n_ * sizeof (T*));
            n_ = 0;
        }

};

template <class T>
class SparseMetadataArray :
    public CountableArray<T>,
    public MempoolVirtualAddressMalloc
{
    protected:
        using CountableArray<T>::items_;
        using CountableArray<T>::size_;


    public:
        SparseMetadataArray(
            MemoryPool* mem,
            uli size
        )
            :
            CountableArray<T>(
                size,
                reinterpret_cast<T**>(
                    mem->get(size * sizeof(T*))
                )
            )
        {
        }

        void realign_memory_pool(void* oldptr, void* newptr)
        {
            realign_pointer(items_, oldptr, newptr);
            T** it(items_);
            T** end(items_ + size_);
            for ( ; it != end; ++it)
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

template <class T>
class Tuple :
    public smartptr::Countable
{
    
    public:
        typedef typename std::vector<T>::const_iterator iterator;

    private:
        std::vector<T> items_;

    public:
        Tuple(uli n, T tmpl)
            : items_(n, tmpl)
        {
        }

        Tuple()
        {
        }

        ~Tuple()
        {
        }

        iterator begin() const {return items_.begin();}
        iterator end() const {return items_.end();}

        int size() const {return items_.size();}

        void append(T val)
        {
            items_.push_back(val);
        }

        const T& get(uli n) const {return items_[n];}

        void set(uli n, const T& item) { items_[n] = item; }
};

template <class T>
std::ostream&
operator<<(std::ostream& os, const boost::intrusive_ptr<Tuple<T> >& tuple)
{
    if (tuple->size() == 0)
        return os;

    typename Tuple<T>::iterator it(tuple->begin());

    os << *it;
    ++it;

    for ( ; it != tuple->end(); ++it)
    {
        os << std::endl << *it;
    }

    return os;
}

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
