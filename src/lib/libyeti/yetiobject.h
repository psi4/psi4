#ifndef YETIOBJECT_H
#define YETIOBJECT_H

#include "class.h"
#include "thread.h"
#include "tensor.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

/**
    @class YetiRuntimeObject
    Class which governs thread safety and retrieval.  A thread lock and reference count
    is used for objects that have to be retrieved before us.  The reference count ensures
    that multiple threads can use the same object and do not release objeccts being used
    by another thread.
*/
class YetiRuntimeObject {

    protected:
        uli runtime_count_;

        bool initialized_;

        bool finalized_;

        virtual void _initialize();

        virtual void _release();

        virtual void _retrieve();

        virtual void _obsolete();

        virtual void _renew();

        void set_initialized(bool flag);

        void set_finalized(bool flag);

    public:
        typedef enum {thread_safe, not_thread_safe} thread_safety_flag_t;

        YetiRuntimeObject();

        virtual ~YetiRuntimeObject();

        /**
            @return Whether the class is currently retrieved by a thread
        */
        bool is_retrieved() const;

        bool is_initialized() const;

        bool is_finalized() const;

        /**
            Perform the retrieval operation if necessary
            and increment the reference count. The reference
            count BEFORE retrieval is returned.  Usually, this is
            0 if not retrieved providing essentially a boolean
            flag
        */
        virtual uli retrieve();

        /**
            Perform the release operation if necessary and increment the reference count
        */
        virtual void release();

        virtual void obsolete();

        virtual void initialize();

};

class YetiThreadedRuntimeObject :
    public YetiRuntimeObject
{

    protected:
        ThreadLock* lock_;

        bool locked_;

    public:
        YetiThreadedRuntimeObject();

        YetiThreadedRuntimeObject(YetiRuntimeObject::thread_safety_flag_t flag);

        virtual ~YetiThreadedRuntimeObject();

        void initialize();

        uli retrieve();

        uli retrieve_lock();

        uli retrieve_nolock();

        bool is_locked() const;

        void release();

        void release_lock();

        void release_nolock();

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

        void obsolete();

        void set_lock(ThreadLock* lock);

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
    protected:
        YetiRuntimeCountable();

        YetiRuntimeCountable(YetiRuntimeObject::thread_safety_flag_t flag);
};

}

namespace boost {
    void intrusive_ptr_add_ref(yeti::YetiRuntimeSerializable* obj);
    void intrusive_ptr_release(yeti::YetiRuntimeSerializable* obj);
    void intrusive_ptr_add_ref(yeti::YetiRuntimeCountable* obj);
    void intrusive_ptr_release(yeti::YetiRuntimeCountable* obj);
}

#ifdef redefine_size_t
#undef size_t
#endif

#endif // YETIOBJECT_H
