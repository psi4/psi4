#ifndef smartptr_ref_h
#define smartptr_ref_h

//do not use default Boost IO stream behavior on smart ptrs
#define BOOST_NO_IOSTREAM

namespace smartptr {

/**
    @class Countable The base class from which all reference counted
                     classes must derive
*/
class Countable {
    
    private:
        unsigned int refcount_;

    public:
        Countable();

        /**
            Increments the reference count by one.
            @return The new reference count
        */
        unsigned int incref();

        /**
            Decrements the reference count by one.  Aborts
            if reference count is already zero.
            @return The new reference count
        */
        unsigned int decref();

        /**
            @return The current reference count
        */
        unsigned int nref() const;

        virtual ~Countable(){}

};

}

namespace boost {

/**
    Boost function to increase the reference count of const objects
    @param c Reference counted object
*/
void
intrusive_ptr_add_ref(const smartptr::Countable* c);

/**
    Boost function to increase the reference count of objects
    @param c Reference counted object
*/
void
intrusive_ptr_add_ref(smartptr::Countable* c);

/**
    Boost function to decrease the reference count of const objects
    @param c Reference counted object
*/
void
intrusive_ptr_release(const smartptr::Countable* c);

/**
    Boost function to decrease the reference count of objects
    @param c Reference counted object
*/
void
intrusive_ptr_release(smartptr::Countable* c);

}

#include "boost/intrusive_ptr.hpp"

#endif
