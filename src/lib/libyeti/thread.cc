#include "threadimpl.h"
#include "exception.h"

#include <errno.h>

using namespace yeti;
using namespace std;


DECLARE_SUBMALLOC(ThreadLock,pThreadLock);

list<ThreadWorkspaceAllocator*>* ThreadEnvironment::allocators_ = 0;


void*
pthread_run(void* thrptr)
{
    Thread* thr = static_cast<Thread*>(thrptr);
    thr->run();
    return 0;
}

void
ThreadEnvironment::allocate()
{
    list<ThreadWorkspaceAllocator*>::const_iterator it(allocators_->begin());
    list<ThreadWorkspaceAllocator*>::const_iterator stop(allocators_->end());
    for ( ; it != stop; ++it)
    {
        (*it)->allocate();
    }
}

void
ThreadEnvironment::deallocate()
{
    list<ThreadWorkspaceAllocator*>::const_iterator it(allocators_->begin());
    list<ThreadWorkspaceAllocator*>::const_iterator stop(allocators_->end());
    for ( ; it != stop; ++it)
    {
        (*it)->deallocate();
    }
}

void
ThreadEnvironment::add_allocator(ThreadWorkspaceAllocator *allocator)
{
    if (allocators_ == 0) //allocate
        allocators_ = new list<ThreadWorkspaceAllocator*>;
    allocators_->push_back(allocator);
}

ThreadLock::~ThreadLock()
{
}

pThreadLock::pThreadLock()
{
    int signal = pthread_mutexattr_init(&attr_);
    if (signal != 0)
    {
        cerr << "mutex attr init error " << signal << endl;
        abort();
    }

    signal = pthread_mutex_init(&mutex_, &attr_);
    if (signal != 0)
    {
        cerr << "mutex init error " << signal << endl;
        abort();
    }
}

pThreadLock::~pThreadLock()
{
    int signal = pthread_mutexattr_destroy(&attr_);
    if (signal != 0)
    {
        cerr << "mutex attr destroy error " << signal << endl;
        abort();
    }

    signal = pthread_mutex_destroy(&mutex_);
    if (signal != 0)
    {
        cerr << "mutex destroy error " << signal << endl;
        abort();
    }
}

void
pThreadLock::lock()
{
    int signal = pthread_mutex_lock(&mutex_);
    if (signal != 0)
    {
        cerr << "mutex lock error " << signal << " on " << &mutex_ << endl;
        abort();
    }
}

bool
pThreadLock::trylock()
{
    int signal = pthread_mutex_trylock(&mutex_);
    return signal == 0; //signal of 0 means lock successful
}

void
pThreadLock::unlock()
{
    pthread_mutex_unlock(&mutex_);
}

pthread_mutex_t*
pThreadLock::mutex()
{
    return &mutex_;
}

void
pThreadLock::print(std::ostream &os) const
{
    os << "pThread Lock" << endl;
}

ThreadGroup::ThreadGroup(uli nthread)
  :
    nthread_max_(nthread),
    threads_(new Thread*[nthread]),
    nthread_(0)
{
}

void
ThreadGroup::clear()
{
    for (uli i=0; i < nthread_; ++i)
    {
        delete threads_[i];
        threads_[i] = 0;
    }

    nthread_ = 0;
}

pThreadGroup::pThreadGroup(uli nthread)
    : ThreadGroup(nthread),
    pthreads_(new pthread_t[nthread]),
    attrs_(new pthread_attr_t[nthread])
{
    ::memset(threads_, 0, nthread * sizeof(void*));
    ::memset(pthreads_, 0, nthread * sizeof(void*));
    for (uli i=0; i < nthread; ++i)
    {
        pthread_attr_init(&attrs_[i]);

        size_t stacksize = 1e8;
        pthread_attr_setstacksize (&attrs_[i], stacksize);
    }
}

void
pThreadGroup::clear()
{
    for (uli i=0; i < nthread_; ++i)
    {
    }
    ThreadGroup::clear();
}

void
pThreadGroup::add(Thread* thr)
{
    if (nthread_ >= nthread_max_)
    {
        raise(SanityCheckError, "too many threads added to thread group");
    }

    threads_[nthread_] = thr;
    ++nthread_;
}

void
pThreadGroup::run()
{
    if (nthread_ > 1)
        YetiRuntime::set_threaded_runtime(true);

    for (uli t=1; t < nthread_; ++t)
    {
        int status = pthread_create(&pthreads_[t], &attrs_[t], pthread_run, reinterpret_cast<void*>(threads_[t]));
        if (status != 0)
        {
            cerr << "thread " << t << " not created properly" << endl;
            abort();
        }
    }

    threads_[0]->run();
}

void
pThreadGroup::wait()
{
    void* retval;
    for (uli t=1; t < nthread_; ++t)
    {
        int status = pthread_join(pthreads_[t], &retval);
        if (status != 0)
        {
            cerr << "thread " << t << " returned error status" << endl;
            abort();
        }

        if (retval != 0)
        {
            cerr << "thread " << t << " returned error status" << endl;
            abort();
        }
    }


    YetiRuntime::set_threaded_runtime(false);
}


ThreadLock*
pThreadGroup::get_lock()
{
    ThreadLock* lock;
    if (nthread_max_ > 1)
        lock = new pThreadLock;
    else
        lock = new NullThreadLock;

    return lock;
}

void
NullThreadLock::lock()
{
    //do nothing
}

void
NullThreadLock::unlock()
{
}

bool
NullThreadLock::trylock()
{
    return true; //always free
}

void
NullThreadLock::print(std::ostream &os) const
{
    os << "NullThreadLock" << endl;
}

Thread::Thread(uli threadnum)
    : threadnum_(threadnum)
{
}

Thread::~Thread()
{
}


