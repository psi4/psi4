#ifndef smartptr_set_h
#define smartptr_set_h

#include <vector>
#include "ref.h"

namespace smartptr {

template <class T>
class SetType {

    public:
        typedef T type;
};

template <class T>
class SetType<const T>
{
    public:
        typedef T type;
};

/**
    @class Set
    Container that automatically converts types. The vector constructor does not
    allow automatic conversion of dynamic types, e.g.  A vector of SimpleInternalCoordinates
    cannot be directly converted to a vector of InternalCoordinates.  The implicit
    template constructor of Set automatically performs type conversions.
*/
template <class T>
class Set {
    
    public:
        typedef typename SetType<T>::type val_t;
        typedef typename std::vector<val_t>::const_iterator const_iterator;
        typedef typename std::vector<val_t>::iterator iterator;

    private:
        /**
            The underlying vector containing the items
        */
        std::vector<val_t> items_;

    public:
        /**
            Constructor for initializing empty set
        */
        Set()
        {
        }

        /**
            Template constructor for initializing Set from vector.
            Template Y must be a derived type of T or compile-time
            error will be thrown.
            @param v
        */
        template <class Y>
        Set(const std::vector<Y>& v) 
        {
            //we need to recreate the vector as the new elements
            //inefficient, but otherwise no type safety
            typename std::vector<Y>::const_iterator it;
            for (it = v.begin(); it != v.end(); ++it)
                items_.push_back(*it); //do the type conversion
        }

        /**
            Template constructor for initializing Set.
            Template Y must be a derived type of T or compile-time
            error will be thrown.
             @param v
        */
        template <class Y>
        Set(const Set<Y>& set) 
        {
            //we need to recreate the vector as the new elements
            //inefficient, but otherwise no type safety
            typename Set<Y>::const_iterator it;
            for (it = set.begin(); it != set.end(); ++it)
                items_.push_back(*it);
        }

        /**
            Template method for appending methods to set.
            Template Y must be a derived type of T or compile-time
            error will be thrown.
            @param element
        */
        template <class Y>
        void
        append(const Y& element)
        {
            items_.push_back(element);
        }

        iterator begin() {return items_.begin();}

        iterator end() {return items_.end();}

        const_iterator begin() const {return items_.begin();}

        const_iterator end() const {return items_.end();}

        int size() const {return items_.size();}

        const T& operator[](int n) const {return items_[n];}

        
};

}

#endif
