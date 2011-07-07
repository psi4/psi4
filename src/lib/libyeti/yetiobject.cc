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
   runtime_count_(0)
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
YetiRuntimeObject::_obsolete()
{
}

uli
YetiRuntimeObject::retrieve()
{
    ++runtime_count_;
    if (runtime_count_ == 1)
    {
        _retrieve();
    }
    return runtime_count_ - 1;
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

bool
YetiRuntimeObject::is_retrieved() const
{
    return runtime_count_;
}

YetiThreadedRuntimeObject::YetiThreadedRuntimeObject()
    :
   lock_(0)
{
    lock_ = YetiRuntime::get_thread_grp()->get_lock();
}

YetiThreadedRuntimeObject::YetiThreadedRuntimeObject(
    YetiRuntimeObject::thread_safety_flag_t flag
)
    :
   lock_(0)
{
    if (flag == YetiRuntimeObject::thread_safe)
        lock_ = new NullThreadLock;
    else
        lock_ = YetiRuntime::get_thread_grp()->get_lock();
}


YetiThreadedRuntimeObject::~YetiThreadedRuntimeObject()
{
    if (YetiRuntime::is_threaded_runtime())
        lock_->unlock();


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
    if (YetiRuntime::is_threaded_runtime())
    {
        lock_->lock();
    }
}

void
YetiThreadedRuntimeObject::unlock()
{
    if (YetiRuntime::is_threaded_runtime())
    {
        lock_->unlock();
    }
}

bool
YetiThreadedRuntimeObject::trylock()
{
    if (YetiRuntime::is_threaded_runtime())
    {
        return lock_->trylock();
    }
    else
    {
        return true;
    }
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
    _obsolete();
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
   if (YetiRuntime::nthread() > 1)
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
    if (YetiRuntime::is_threaded_runtime())
    {
        obj->lock();
        obj->incref();
        obj->unlock();
    }
    else
    {
        obj->incref();
    }
}

template <class T>
void
intrusive_ptr_decref(T* obj)
{
    if (YetiRuntime::is_threaded_runtime())
    {
        obj->lock();
        //there is no need to "unlock" as the obj is getting deleted anyway
        if (obj->decref() == 0)
            delete obj;
        else
            obj->unlock();
    }
    else
    {
        //there is no need to "unlock" as the obj is getting deleted anyway
        if (obj->decref() == 0)
            delete obj;
    }
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
