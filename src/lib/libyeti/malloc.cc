#include "malloc.h"
#include "thread.h"
#include "runtime.h"
#include "exception.h"

#include <iostream>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

std::map<std::string, int>* FastMalloc::registration_queue_ = 0;


DECLARE_MALLOC(MemoryPool);

#define MALLOC_FAST_MALLOC 0

#if MALLOC_FAST_MALLOC
static std::map<void*, uli> malloc_numbers_;
static std::map<size_t, std::map<uli, void*> > malloc_objects_;
static std::map<size_t, uli> malloc_counts_;
#endif

FastMalloc::FastMalloc(
    uli size,
    uli n,
    const std::string& name
)
    :
    data_(0),
    size_(0),
    n_(0),
    offset_(0),
    name_(name),
    static_mem_(false)
{
    init(size, n);
}

FastMalloc::FastMalloc(
    char* data,
    char* mallocd,
    uli size,
    uli n,
    const std::string& name
)
    :
    data_(data),
    size_(size),
    n_(n),
    offset_(0),
    name_(name),
    mallocd_(mallocd),
    lock_(0),
    static_mem_(true)
{
    memset(mallocd_, 0, n + 1);
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
    offset_(0),
    lock_(0),
    static_mem_(false)
{
}


void
FastMalloc::init(uli size, uli n)
{
    if (data_)
    {
        cerr << "Malloc for " << name_ << " has already been initialized" << endl;
        abort();
    }

    n_ = n;
    size_ = size;

    MallocOverride o;
    lock_ = new (o) pThreadLock;

#if MALLOC_FAST_MALLOC
    return;
#endif

    data_ = YetiRuntime::malloc(size_ * n_);
    mallocd_ = YetiRuntime::malloc(n_ + 1);
    memset(mallocd_, 0, n + 1);
    offset_ = 0;
}

FastMalloc::~FastMalloc()
{
#if MALLOC_FAST_MALLOC
    return;
#endif
    char* mallocptr = mallocd_;
    uli nmalloc = 0;
    for (uli i=0; i < n_; ++i, ++mallocptr)
    {
        if (*mallocptr)
        {
            ++nmalloc;
        }
    }

    if (lock_) ::free(lock_);

    if (nmalloc)
    {
        cerr << stream_printf("There are still %ld unfreed entries on malloc for %s on node %ld",
                              nmalloc, this->name_.c_str(), YetiRuntime::me()) << endl;
    }
    else if (!static_mem_) //everything is clear for deletion
    {
        YetiRuntime::free(data_, size_ * n_);
        YetiRuntime::free(mallocd_, n_ + 1);
    }


}

uli
FastMalloc::nfree() const
{
    char* ptr = mallocd_;
    uli nblocks_free = 0;
    for (uli i=0; i < n_; ++i, ++ptr)
    {
        if (*ptr)
            ++nblocks_free;
    }
    return nblocks_free;
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
FastMalloc::malloc()
{
#if MALLOC_FAST_MALLOC
    YetiRuntime::lock_malloc();
    void* ptr = ::malloc(size_);
    if (ptr == 0)
    {
        cerr << "malloc failed to allocate block of size " << size_ << endl;
        abort();
    }

    uli malloc_number = malloc_counts_[size_] + 1;
    malloc_counts_[size_] = malloc_number;
    malloc_numbers_[ptr] = malloc_number;
    malloc_objects_[size_][malloc_number] = ptr;


    YetiRuntime::unlock_malloc();
    return ptr;
#else
    if (lock_) lock_->lock();


    void* ptr = malloc_no_lock();

    if (lock_) lock_->unlock();


    return ptr;
#endif
}

void*
FastMalloc::malloc_no_lock()
{
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
    return ptr;
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

uli
FastMalloc::get_malloc_number(void* obj)
{
#if MALLOC_FAST_MALLOC
    return malloc_numbers_[obj];
#else
    size_t diff = (size_t) obj - (size_t) data_;
    uli n = diff / size_;
    if (n >= n_)
    {
	cerr << "get malloc number " << n << " is not valid for "
             << name_ << " since there are only " << n_ << " entries" << endl;
        yeti_throw(SanityCheckError, "malloc number exceeds size");
    }
    return n;
#endif
}

void*
FastMalloc::get_object(uli malloc_number)
{
#if MALLOC_FAST_MALLOC
    std::map<uli,void*>::const_iterator it = malloc_objects_[size_].find(malloc_number);
    if (it == malloc_objects_[size_].end())
    {
        cerr << "Object number " << malloc_number
            << " is not valid for " << name_ 
            << endl;
        abort();
    }
    return malloc_objects_[size_][malloc_number];
#else
#if YETI_SANITY_CHECK
    if (malloc_number >= n_)
    {
	cerr << "get object " << malloc_number << " is not valid for "
             << name_ << " since there are only " << n_ << " entries" << endl;
        yeti_throw(SanityCheckError, "malloc number exceeds size");
    }
#endif
    return data_ + malloc_number * size_;
#endif
}

MemoryPool::MemoryPool(
    size_t size,
    char* data
)
    : data_(data),
        total_size_(size),
        remaining_(size),
        ptr_(data),
        mallocd_(false)
{
}


MemoryPool::MemoryPool(size_t size)
    : data_(0),
        total_size_(size),
        remaining_(size),
        ptr_(data_),
        mallocd_(true)
{
    data_ = YetiRuntime::malloc(size);
}

MemoryPool::~MemoryPool()
{
    if (mallocd_) YetiRuntime::free(data_, total_size_);
}

char*
MemoryPool::get(uli size)
{
    if (size > remaining_)
    {
        cerr << "Insufficient memory in pool. Need " << size << " bytes"
                " but only have " << remaining_ << " bytes out of total "
            << total_size_ << endl;
        abort();
    }

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
MemoryPool::total_size() const
{
    return total_size_;
}

size_t
MemoryPool::data_size() const
{
    return total_size_ - remaining_;
}

void
MemoryPool::memset()
{
    ::memset(data_, 0, total_size_);
}

void
MemoryPool::memcpy(MemoryPool* pool)
{
    size_t cpysize = pool->total_size_ - pool->remaining_;
    if (cpysize > total_size_)
    {
        cerr << "Memory pool is too large for memcpy!" << endl;
        cerr << "Parent pool is copying " << cpysize << " of " << pool->total_size_
             << " bytes, but destination only has "
             << total_size_ << " bytes"
             << endl;
        abort();
    }

    ::memcpy(data_, pool->data_, cpysize);

    remaining_ = total_size_ - cpysize;
    ptr_ = data_ + cpysize;
}

void
MemoryPool::set_data_size(size_t size)
{
    ptr_ = data_ + size;
    remaining_ = total_size_ - size;
}

void
MemoryPool::set(char* data)
{
    data_ = data;
    ptr_ = data_ + total_size_ - remaining_;
}

void
MemoryPool::reset()
{
    ptr_ = data_;
    remaining_ = total_size_;
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
