#include "yetiobject.h"
#include "runtime.h"
#include "thread.h"
#include "tensor.h"
#include "tensorblock.h"
#include "contraction.h"
#include "exception.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

using namespace yeti;
using namespace std;

YetiRuntimeObject::YetiRuntimeObject()
    :
   runtime_count_(0),
   initialized_(false),
   finalized_(true)
{
}

YetiRuntimeObject::~YetiRuntimeObject()
{
}

void
YetiRuntimeObject::_retrieve()
{
}

void
YetiRuntimeObject::obsolete()
{
   if (!initialized_ || finalized_) //already obsolete
        return;
   _obsolete();
}

void
YetiRuntimeObject::_obsolete()
{
}

uli
YetiRuntimeObject::retrieve()
{
    if (!initialized_)
    {
        if (runtime_count_ != 0)
            yeti_throw(SanityCheckError, "retrieved object was not initialized");
        YetiRuntimeObject::initialize();
    }

    if (runtime_count_ > 0)
    {
        ++runtime_count_;
        return runtime_count_;
    }


    if (finalized_)
    {
        _retrieve();
        finalized_ = false;
    }
    else
    {
        _renew();
    }

    ++runtime_count_;
    return runtime_count_ - 1;
}

void
YetiRuntimeObject::_renew()
{
}

void
YetiRuntimeObject::_release()
{
}

void
YetiRuntimeObject::release()
{
    if (runtime_count_ == 0)
    {
        cerr << "invalid negative runtime count" << endl;
        abort();
    }

    if (runtime_count_ == 1)
    {
        _release();
    }
    --runtime_count_;
}

void
YetiRuntimeObject::_initialize()
{
}

void
YetiRuntimeObject::initialize()
{
    if (initialized_)
        return;

    _initialize();
    initialized_ = true;
}


bool
YetiRuntimeObject::is_initialized() const
{
    return initialized_;
}

bool
YetiRuntimeObject::is_finalized() const
{
    return finalized_;
}

void
YetiRuntimeObject::set_initialized(bool flag)
{
    initialized_ = flag;
}

void
YetiRuntimeObject::set_finalized(bool flag)
{
    finalized_ = flag;
}

bool
YetiRuntimeObject::is_retrieved() const
{
    return runtime_count_;
}

YetiThreadedRuntimeObject::YetiThreadedRuntimeObject()
    :
   lock_(0),
   locked_(false)
{
    lock_ = YetiRuntime::get_thread_grp()->get_lock();
}

YetiThreadedRuntimeObject::YetiThreadedRuntimeObject(
    YetiRuntimeObject::thread_safety_flag_t flag
)
    :
   lock_(0),
   locked_(false)
{
    if (flag == YetiRuntimeObject::thread_safe)
    {
        lock_ = new NullThreadLock;
    }
    else
    {
        lock_ = YetiRuntime::get_thread_grp()->get_lock();
    }
}


YetiThreadedRuntimeObject::~YetiThreadedRuntimeObject()
{
#if YETI_SANITY_CHECK
    if (is_retrieved())
    {
        yeti_throw(SanityCheckError, "cannot delete retrieved object");
    }
#endif
    
    //lock_->unlock();
    delete lock_;
}

void
YetiThreadedRuntimeObject::set_lock(ThreadLock* lock)
{
    delete lock_;
    lock_ = lock;
}

void
YetiThreadedRuntimeObject::lock()
{
    lock_->lock();
    locked_ = true;
}

void
YetiThreadedRuntimeObject::unlock()
{
#if YETI_SANITY_CHECK
    if (!locked_)
        yeti_throw(SanityCheckError, "You unlocked a mutex that isn't locked");
#endif
    locked_ = false;
    lock_->unlock();
}

bool
YetiThreadedRuntimeObject::is_locked() const
{
    return locked_;
}

bool
YetiThreadedRuntimeObject::trylock()
{
    //if (locked_)
    //{
    //    cout << "failed lock on " << this << endl;
    //    return false;
    //}

    bool trylock = lock_->trylock();
    if (trylock) 
    {
        locked_ = true;
    }

    return trylock;
}

void
YetiThreadedRuntimeObject::initialize()
{
    lock();
    YetiRuntimeObject::initialize();
    unlock();
}


void
YetiThreadedRuntimeObject::obsolete()
{
    lock();
    if (runtime_count_ != 0)
    {
        cerr << "cannot make retrieved object obsolete" << endl;
        abort();
    }
    YetiRuntimeObject::obsolete();
    unlock();
}

uli
YetiThreadedRuntimeObject::retrieve()
{
    lock();
    uli count = YetiRuntimeObject::retrieve();
    unlock();
    return count;
}

uli
YetiThreadedRuntimeObject::retrieve_lock()
{
    /** multiple compute threads OR multiple processors */
   if (YetiRuntime::nthread() > 1 || YetiRuntime::nproc() > 1)
   {
       pThreadLock* test = dynamic_cast<pThreadLock*>(lock_);
       if (!test)
       {
           cerr << "no thread lock" << endl;
           abort();
       }
   }


    lock();
    uli count = YetiRuntimeObject::retrieve();
    return count;
}

uli
YetiThreadedRuntimeObject::retrieve_nolock()
{
    return YetiRuntimeObject::retrieve();
}

void
YetiThreadedRuntimeObject::release()
{
    lock();
    YetiRuntimeObject::release();
    unlock();
}

void
YetiThreadedRuntimeObject::release_lock()
{
    YetiRuntimeObject::release();
    unlock();
}

void
YetiThreadedRuntimeObject::release_nolock()
{
    return YetiRuntimeObject::release();
}

template <class T>
void
intrusive_ptr_incref(T* obj)
{
    obj->lock();
    obj->incref();
    obj->unlock();
}

template <class T>
void
intrusive_ptr_decref(T* obj)
{
    obj->lock();
    //there is no need to "unlock" as the obj is getting deleted 
    //it's destructor will handle the unlock
    if (obj->decref() == 0)
        delete obj;
    else
        obj->unlock();
}

void
boost::
intrusive_ptr_add_ref(yeti::YetiRuntimeCountable *obj)
{
    intrusive_ptr_incref<YetiRuntimeCountable>(obj);
}

void
boost::
intrusive_ptr_release(yeti::YetiRuntimeCountable *obj)
{
    intrusive_ptr_decref<YetiRuntimeCountable>(obj);
}

void
boost::
intrusive_ptr_add_ref(yeti::YetiRuntimeSerializable *obj)
{
    intrusive_ptr_incref<YetiRuntimeSerializable>(obj);
}

void
boost::
intrusive_ptr_release(yeti::YetiRuntimeSerializable *obj)
{
    intrusive_ptr_decref<YetiRuntimeSerializable>(obj);
}


YetiRuntimeCountable::YetiRuntimeCountable()
{
}

YetiRuntimeCountable::YetiRuntimeCountable(YetiRuntimeObject::thread_safety_flag_t flag)
 : YetiThreadedRuntimeObject(flag)
{
}
