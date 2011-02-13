#include "yetiobject.h"
#include "runtime.h"
#include "thread.h"

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
YetiRuntimeObject::_retrieve(uli threadnum)
{
}

void
YetiRuntimeObject::retrieve(uli threadnum)
{
    if (runtime_count_ == 0)
        _retrieve(threadnum);
    ++runtime_count_;
}

void
YetiRuntimeObject::_release(uli threadnum)
{
}

void
YetiRuntimeObject::release(uli threadnum)
{
    if (runtime_count_ == 0)
    {
        cerr << "invalid negative runtime count" << endl;
        abort();
    }

    if (runtime_count_ == 1)
    {
        _release(threadnum);
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

YetiThreadedRuntimeObject::~YetiThreadedRuntimeObject()
{
    if (YetiRuntime::is_threaded_runtime())
        lock_->unlock();
    delete lock_;
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
YetiThreadedRuntimeObject::retrieve(uli threadnum)
{
    lock();
    YetiRuntimeObject::retrieve(threadnum);
    unlock();
}

void
YetiThreadedRuntimeObject::release(uli threadnum)
{
    lock();
    YetiRuntimeObject::release(threadnum);
    unlock();
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

