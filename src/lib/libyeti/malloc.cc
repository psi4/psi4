#include "malloc.h"
#include "thread.h"
#include "runtime.h"

#include <iostream>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

std::map<std::string, int>* FastMalloc::registration_queue_ = 0;

DECLARE_MALLOC(MemoryPool);

#define MALLOC_FAST_MALLOC 0

FastMalloc::FastMalloc(
    uli size,
    uli n,
    const std::string& name
)
    :
    data_(0),
    size_(0),
    n_(0),
    name_(name),
    thread_offsets_(0),
    thread_starts_(0),
    thread_stops_(0)
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
    refcount_(0),
    thread_offsets_(0),
    thread_starts_(0),
    thread_stops_(0)
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

    data_ = YetiRuntime::malloc(size_ * n_);
    mallocd_ = YetiRuntime::malloc(n_ + 1);
    memset(mallocd_, 0, n + 1);
    offset_ = 0;

    uli nthread = YetiRuntime::nthread();
    thread_offsets_ = new uli[nthread];
    thread_starts_ = new uli[nthread];
    thread_stops_ = new uli[nthread];

    uli increment = n / nthread;
    uli offset = 0;
    for (uli i=0; i < nthread; ++i)
    {
        thread_offsets_[i] = offset;
        thread_starts_[i] = offset;
        offset += increment;
        thread_stops_[i] = offset;
    }
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
        cerr << stream_printf("There are still %ld unfreed entries on malloc for %s",
                              nmalloc, this->name_.c_str()) << endl;
    }
    else //everything is clear for deletion
    {
        YetiRuntime::free(data_, size_ * n_);
        YetiRuntime::free(mallocd_, n_ + 1);

        if (thread_offsets_)
            delete[] thread_offsets_;
        if (thread_starts_)
            delete[] thread_starts_;
        if (thread_stops_)
            delete[] thread_stops_;
    }


}

uli
FastMalloc::search(
    uli start,
    uli stop
)
{
    char* ptr = mallocd_ + start;
    uli index = start;
    for ( ; *ptr && index != stop; ++index, ++ptr);
    return index;
}

void*
FastMalloc::malloc(uli threadnum)
{
#if MALLOC_FAST_MALLOC
    YetiRuntime::lock_malloc();
    void* ptr = ::malloc(size_);
    YetiRuntime::unlock_malloc();
    return ptr;
#else
    uli offset = thread_offsets_[threadnum];
    uli stop = thread_stops_[threadnum];
    //maybe an empty spot
    uli index = search(offset, stop);
    if (index == stop) //try again
    {
        offset = thread_starts_[threadnum];
        index = search(offset, stop);
        if (index == stop)
        {
            //everything must be allocated...
            cerr << "Exceeded the " << n_ << " entries in malloc for class " << name_ << endl;
            abort();
        }
    }

    thread_offsets_[threadnum] = index + 1;
    mallocd_[index] = 1;
    char* ptr = data_ + index * size_;
    return ptr;
#endif
}

void*
FastMalloc::malloc()
{
    if (YetiRuntime::is_threaded_runtime())
        return malloc(YetiRuntime::get_thread_number());

#if MALLOC_FAST_MALLOC
    YetiRuntime::lock_malloc();
    void* ptr = ::malloc(size_);
    if (ptr == 0)
    {
        cerr << "malloc failed to allocate block of size " << size_ << endl;
        abort();
    }
    YetiRuntime::unlock_malloc();
    return ptr;
#else
    if (YetiRuntime::is_threaded_runtime())
    {
        lock_->lock();
    }

    //maybe an empty spot
    uli index = search(offset_, n_);
    if (index == n_) //try again
    {
        offset_ = 0;
        index = search(offset_, n_);
        if (index == n_)
        {
            //everything must be allocated...
            cerr << "Exceeded the " << n_ << " entries in malloc for class " << name_ << endl;
            abort();
        }
    }

    offset_ = index + 1;
    mallocd_[index] = 1;
    char* ptr = data_ + index * size_;

    if (YetiRuntime::is_threaded_runtime())
        lock_->unlock();


    return ptr;
#endif
}

void
FastMalloc::free(const void* entry)
{
#if MALLOC_FAST_MALLOC
    ::free(const_cast<void*>(entry));
#else
    uli offset = (reinterpret_cast<const char*>(entry) - data_) / size_;
    *(mallocd_ + offset) = 0;
#endif
}

void
FastMalloc::queue(const std::string &name)
{
    if (!registration_queue_)
        registration_queue_ = new std::map<std::string, int>;

    std::map<std::string, int>::iterator it(registration_queue_->find(name));
    if (it == registration_queue_->end())
    {
        (*registration_queue_)[name] = 1;
    }


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

MemoryPool::MemoryPool(
    size_t size,
    char* data
)
    : data_(data),
        size_(size),
        remaining_(size),
        ptr_(data),
        mallocd_(false)
{
}


MemoryPool::MemoryPool(size_t size)
    : data_(0),
        size_(size),
        remaining_(size),
        ptr_(data_),
        mallocd_(true)
{
    data_ = YetiRuntime::malloc(size);
}

MemoryPool::~MemoryPool()
{
    if (mallocd_) YetiRuntime::free(data_, size_);
}

char*
MemoryPool::get(uli size)
{
    if (size > remaining_)
    {
        cerr << "Insufficient memory in pool. Need " << size << " bytes"
                " but only have " << remaining_ << " bytes out of total "
            << size_ << endl;
        abort();
    }

#if 0
    dout << stream_printf("allocated %d bytes on pool of size %d with %d remaining",
                          size, size_, remaining_
                          ) << endl;
#endif

    char* tmp = ptr_;

    remaining_ -= size;
    ptr_ += size;

    return tmp;
}

char*
MemoryPool::data() const
{
    return data_;
}

size_t
MemoryPool::remaining() const
{
    return remaining_;
}

size_t
MemoryPool::size() const
{
    return size_;
}

void
MemoryPool::memset()
{
    ::memset(data_, 0, size_);
}

void
MemoryPool::memcpy(MemoryPool* pool)
{
    size_t cpysize = pool->size_ - pool->remaining_;
    if (cpysize > size_)
    {
        cerr << "Memory pool is too large for memcpy!" << endl;
        cerr << "Parent pool is copying " << cpysize << " of " << pool->size_
             << " bytes, but destination only has "
             << size_ << " bytes"
             << endl;
        abort();
    }

    ::memcpy(data_, pool->data_, cpysize);

    remaining_ = size_ - cpysize;
}

void
MemoryPool::set(char* data)
{
    data_ = data;
    ptr_ = data_ + size_ - remaining_;
}

void
MemoryPool::reset()
{
    ptr_ = data_;
    remaining_ = size_;
}

void
MempoolVirtualAddressMalloc::operator delete(void* ptr)
{
    //do nothing
}

void*
MempoolVirtualAddressMalloc::operator new(size_t size, MemoryPool* mem)
{
    //return ::malloc(size);
    return mem->get(size);
}

void
MempoolVirtualAddressMalloc::operator delete(void* ptr, MemoryPool* mem)
{
    cerr << "Something failed in memory pool allocation" << endl;
    abort();
}
