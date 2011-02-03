#ifndef yeti_tuple_h
#define yeti_tuple_h

#include <vector>
#include "class.h"

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


template <class T>
class CountableArray : public smartptr::Serializable {

    private:
        T** items_;

        uli size_;

        uli n_;

    public:
        typedef T** iterator;

        CountableArray(uli size)
            : items_(new T*[size]), size_(size), n_(0)
        {
            for (unsigned long i=0; i < size_; ++i)
                items_[i] = 0;
        }

        CountableArray(uli size, T* tmpl)
            : items_(new T*[size]), size_(size), n_(0)
        {
            for (unsigned long i=0; i < size_; ++i)
            {
                items_[i] = tmpl;
                if (tmpl) tmpl->incref();
            }
        }

        ~CountableArray()
        {
            T** it(items_);
            T** end(items_ + size_);
            for ( ; it != end; ++it)
            {
                T* next(*it);
                if (next)
                {
                    //std::cout << "releasing item=" << next << std::endl;
                    boost::intrusive_ptr_release(next);
                }
            }
            delete[] items_;
        }

        T* get(size_t n) const 
        {
            return items_[n];
        }

        void insert(uli n, T* item)
        {
            if (n >= size_)
            {
                std::cerr << "index too large" << std::endl;
                abort();
            }

            if (items_[n] == 0)
                ++n_;

            //std::cout << "inserted item=" << item << " at n=" << n << std::endl;

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
            T** end(items_ + n_);
            for ( ; it != end; ++it)
            {
                T* next(*it);
                if (next)
                    boost::intrusive_ptr_release(next);
                *it = 0;
            }
            n_ = 0;
        }

        uli n() const
        {
            return n_;
        }

        uli size() const
        {
            return size_;
        }


};

template <class T>
class Tuple : public smartptr::Serializable {
    
    public:
        typedef typename std::vector<T>::const_iterator iterator;

    private:
        std::vector<T> items_;

    public:
        Tuple(size_t n, T tmpl)
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

#endif
