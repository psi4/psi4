#include "malloc.h"
#include "thread.h"
#include "runtime.h"

#include <iostream>

using namespace yeti;
using namespace std;

std::map<std::string, int>* FastMalloc::registration_queue_ = 0;

#define MALLOC_FAST_MALLOC 1

FastMalloc::FastMalloc(
    size_t size,
    size_t n,
    const std::string& name
)
    :
    data_(0),
    size_(0),
    n_(0),
    name_(name)
{

    init(size, n);
}

FastMalloc::FastMalloc(
    const std::string& name
)
    :
    data_(0),
    size_(0),
    n_(0),
    name_(name),
    refcount_(0)
{
}


void
FastMalloc::init(uli size, uli n)
{
    if (data_)
    {
        cerr << "malloc for " << name_ << " has already been initialized" << endl;
        abort();
    }

    n_ = n;
    size_ = size;

    MallocOverride o;
    lock_ = new (o) pThreadLock;

    data_ = new char[size * n];
    mallocd_ = new char[n + 1];
    memset(mallocd_, 0, n + 1);
    dataptr_ = data_;
    mallocptr_ = mallocd_;
    index_ = 0;
}

FastMalloc::~FastMalloc()
{
    char* mallocptr = mallocd_;
    uli nmalloc = 0;
    for (uli i=0; i < n_; ++i, ++mallocptr)
    {
        if (*mallocptr)
        {
            ++nmalloc;
        }
    }

    ::free(lock_);

    if (nmalloc)
    {
        cout << stream_printf("There are still %ld unfreed entries on malloc for %s",
                              nmalloc, this->name_.c_str()) << endl;
    }
    else //everything is clear for deletion
    {
        if (data_)
            delete[] data_;
        if (mallocd_)
            delete[] mallocd_;
    }


}

void
FastMalloc::search()
{
    for ( ; *mallocptr_ && index_ < n_; ++mallocptr_, dataptr_ += size_, ++index_);
}

void*
FastMalloc::malloc()
{
#if MALLOC_FAST_MALLOC
    YetiRuntime::lock_malloc();
    void* ptr = ::malloc(size_);
    YetiRuntime::unlock_malloc();
    return ptr;
#else
    if (YetiRuntime::is_threaded_runtime())
    {
        lock_->lock();
    }
    ++refcount_;

    //maybe an empty spot
    search();
    if (index_ == n_) //try again
    {
        index_ = 0;
        dataptr_ = data_;
        mallocptr_ = mallocd_;
        search();
        if (index_ == n_)
        {
            //everything must be allocated...
            cerr << "Exceeded the " << n_ << " entries in malloc for class " << name_ << endl;
            abort();
        }

        if (refcount_ == 2)
            abort();
    }

    *mallocptr_ = 1;

    void* ptr = dataptr_;

    --refcount_;
    if (YetiRuntime::is_threaded_runtime())
        lock_->unlock();


    return ptr;
#endif
}

void
FastMalloc::free(const void* entry)
{
#if MALLOC_FAST_MALLOC
    YetiRuntime::lock_malloc();
    ::free(const_cast<void*>(entry));
    YetiRuntime::unlock_malloc();
#else
    if (YetiRuntime::is_threaded_runtime())
        lock_->lock();

    size_t offset = (reinterpret_cast<const char*>(entry) - data_) / size_;
    *(mallocd_ + offset) = 0;

    if (YetiRuntime::is_threaded_runtime())
        lock_->unlock();

#if YETI_DEBUG_MALLOC
    std::map<const void*,std::string>::iterator it = tracker_.find(entry);
    if (it != tracker_.end())
    {
        tracker_.erase(it);
    }
#endif

#endif
}

void
FastMalloc::queue(const std::string &name)
{
    if (!registration_queue_)
        registration_queue_ = new std::map<std::string, int>;

    std::map<std::string, int>::iterator it(registration_queue_->find(name));
    if (it != registration_queue_->end())
    {
        cerr << name << " has already been registered for malloc" << endl;
        abort();
    }

    (*registration_queue_)[name] = 1;
}

void
FastMalloc::unqueue(const std::string &name)
{
    std::map<std::string, int>::iterator it(registration_queue_->find(name));
    if (it == registration_queue_->end())
    {
        cerr << name << " has not been registered for malloc" << endl;
        abort();
    }

    registration_queue_->erase(it);
}

void
FastMalloc::check_queue()
{
    if (registration_queue_->size() != 0)
    {
        cerr << "Not all malloc objects registered. Failure on" << endl;
        std::map<std::string, int>::iterator it(registration_queue_->begin());
        std::map<std::string, int>::iterator stop(registration_queue_->end());
        for ( ; it != stop; ++it)
            cerr << it->first << endl;
        abort();
    }
}

MemoryAllocation::MemoryAllocation(size_t size)
    : data_(::malloc(size)),
        size_(size),
        remaining_(size),
        ptr_(reinterpret_cast<char*>(data_))
{
}

MemoryAllocation::~MemoryAllocation()
{
    ::free(data_);
}

void*
MemoryAllocation::get(uli size)
{
    if (size > remaining_)
    {
        cerr << "insufficient memory for allocation" << endl;
        abort();
    }

    void* tmp = ptr_;

    remaining_ -= size;
    ptr_ += size;

    return tmp;
}

uli
MemoryAllocation::size() const
{
    return size_;
}



