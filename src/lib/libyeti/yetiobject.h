#ifndef YETIOBJECT_H
#define YETIOBJECT_H

#include "class.h"
#include "thread.h"

namespace yeti {

/**
    @class YetiRuntimeObject
    Class which governs thread safety and retrieval.  A thread lock and reference count
    is used for objects that have to be retrieved before us.  The reference count ensures
    that multiple threads can use the same object and do not release objeccts being used
    by another thread.
*/
class YetiRuntimeObject {

    private:
        uli runtime_count_;

    protected:
        virtual void _release(uli threadnum);

        virtual void _retrieve(uli threadnum);

    public:
        YetiRuntimeObject();

        virtual ~YetiRuntimeObject();

        /**
            @return Whether the class is currently retrieved by a thread
        */
        bool is_retrieved() const;

        /**
            Perform the retrieval operation if necessary and increment the reference count
        */
        virtual void retrieve(uli threadnum = 0);

        /**
            Perform the release operation if necessary and increment the reference count
        */
        virtual void release(uli threadnum = 0);

};

class YetiThreadedRuntimeObject : public YetiRuntimeObject {

    protected:
        ThreadLock* lock_;

    public:
        YetiThreadedRuntimeObject();

        virtual ~YetiThreadedRuntimeObject();

        void retrieve(uli threadnum);

        void release(uli threadnum);

        /**
            Lock the object and wait until released by all other threads
        */
        void lock();

        /**
            Try to lock the object, but do not block.
            @return True if lock was successful and current thread now holds lock
        */
        bool trylock();

        /**
            Unlock the object.  This should only be called after a successful lock
        */
        void unlock();

};

/**
    @class YetiRuntimeSerializable
    Class for combining yeti runtime operations with serializable attributes
*/
class YetiRuntimeSerializable :
    public smartptr::Serializable,
    public YetiThreadedRuntimeObject
{
};

/**
    @class YetiRuntimeSerializable
    Class for combining yeti runtime operations with countable attributes
*/
class YetiRuntimeCountable :
    public smartptr::Countable,
    public YetiThreadedRuntimeObject
{
};

}

#define USE_YETI_REFCOUNT 1
#if USE_YETI_REFCOUNT

namespace boost {
    void intrusive_ptr_add_ref(yeti::YetiRuntimeSerializable* obj);
    void intrusive_ptr_release(yeti::YetiRuntimeSerializable* obj);
    void intrusive_ptr_add_ref(yeti::YetiRuntimeCountable* obj);
    void intrusive_ptr_release(yeti::YetiRuntimeCountable* obj);
}

#endif


#endif // YETIOBJECT_H
