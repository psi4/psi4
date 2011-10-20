#ifndef yeti_MALLOCIMPL_H
#define yeti_MALLOCIMPL_H

#include "malloc.h"


namespace yeti {

/**
    @class FastMallocTemplate
    Class used for providing a memory pool for a specific class. This allocates
    blocks of the appropriate size for the given template parameter.
*/
template <class T>
class FastMallocTemplate {

    private:
        static uli nblocks_;

        static uli blocksize_;

        static std::string classname_;

    public:
        static FastMalloc* malloc;

        /**
            Constructor for registering a particular memory pool.
            Template parameter T correponds to the parent class.  All
            derived classes will have a special memory allocation.  The block
            sizes, though, will be determined by the constructor template
            parameter U.  For example, DataBlock is T and RecomputedBlock is U.
            This ensures that allocated block sizes for all derived classes of
            DataBlock will be large enough since the largest derived class is
            RecomputedBlock.
            @param u Pointer used for implicit template instantiaion
            @param classname
        */
        template <class U>
        FastMallocTemplate(U* u, const std::string& classname)
        {
            if (sizeof(U) > blocksize_)
                blocksize_ = sizeof(U);
            classname_ = classname;
            FastMalloc::queue(classname);
        }

        ~FastMallocTemplate();

        /**
            Allocate the memory pool

        */
        static void init(uli nblocks)
        {
            FastMalloc::unqueue(classname_);
            nblocks_ = nblocks;
            malloc = new FastMalloc(blocksize_, nblocks_, classname_);
        }

        static size_t blocksize()
        {
            return blocksize_;
        }


};

/**
    @class MallocOverride
    Data type used for placement new. This actually calls malloc in lieu of memory pool
    for particular objects.  This is needed in rare exceptions.
*/
class MallocOverride {
};


template <class T>
class Malloc {

    private:
        static bool initialized_;

        size_t malloc_number_;

    public:

        virtual ~Malloc()
        {
        }

        void operator delete(void* ptr)
        {
            if (initialized_)
                FastMallocTemplate<T>::malloc->free(ptr);
            else; //leak it
        }

        void* operator new(size_t size)
        {
#if YETI_SANITY_CHECK
            if (size > FastMallocTemplate<T>::blocksize())
            {
                std::cerr << "Blocksize requested is too large" << std::endl;
                abort();
            }
#endif
            if (initialized_)
            {
                return FastMallocTemplate<T>::malloc->malloc();
            }
            else
            {
                std::cerr << "Malloc class not initialized" << std::endl;
                abort();
            }
        }

        void* operator new(size_t size, MallocOverride& o)
        {
            return ::malloc(size);
        }

        void operator delete(void* ptr, MallocOverride& o)
        {
        }

        static void init_malloc(uli nblocks)
        {
            FastMallocTemplate<T>::init(nblocks);
            initialized_ = true;
        }

        static void uninit_malloc()
        {
            initialized_ = false;
        }

        size_t get_malloc_number() const
        {
            return malloc_number_;
        }

        static T* get_object(uli malloc_number)
        {
            return static_cast<T*>(FastMallocTemplate<T>::malloc->get_object(malloc_number));
        }

        static uli get_malloc_number(T* object)
        {
            return FastMallocTemplate<T>::malloc->get_malloc_number(object);
        }

};


/**
    @class MemoryAllocation
    Class encapsulating a memory pool which hands out blocks of data
*/
class MemoryPool :
    public smartptr::Countable,
    public Malloc<MemoryPool>
{

    private:
        size_t total_size_;

        size_t remaining_;

        char* data_;

        char* ptr_;

        /**
            Whether this memory allocation
            mallocd the data itself
        */
        bool mallocd_;

    public:
        MemoryPool(
            size_t size,
            char* data
        );

        MemoryPool(size_t size);

        ~MemoryPool();

        char* get(size_t size);

        void memset();

        void set(char* data);

        void set_data_size(size_t size);

        void reset();

        void memcpy(MemoryPool* pool);

        char* data() const;

        size_t total_size() const;

        size_t data_size() const;

        size_t remaining() const;

};

class MempoolVirtualAddressMalloc
{

    public:
        void operator delete(void* ptr);

        void* operator new(size_t size, MemoryPool* mem);

        void operator delete(void* ptr, MemoryPool* mem);
};

template <class T>
void
realign_pointer(
    T*& ptr,
    void* oldptr,
    void* newptr
)
{
    char* final_ptr = newptr > oldptr ?
        reinterpret_cast<char*>(const_cast<T*>(ptr))
            + ((char*) newptr - (char*) oldptr) :
        reinterpret_cast<char*>(const_cast<T*>(ptr))
            - ((char*) oldptr - (char*) newptr) ;
    ptr = reinterpret_cast<T*>(final_ptr);
}

template <class T>
T*
get_realigned_pointer(
    const T* ptr,
    void* oldptr,
    void* newptr
)
{
    char* final_ptr = newptr > oldptr ?
        reinterpret_cast<char*>(const_cast<T*>(ptr))
            + ((char*) newptr - (char*) oldptr) :
        reinterpret_cast<char*>(const_cast<T*>(ptr))
            - ((char*) oldptr - (char*) newptr) ;
    return reinterpret_cast<T*>(final_ptr);
}

template <class T>
FastMallocTemplate<T>::~FastMallocTemplate()
{
    Malloc<T>::uninit_malloc();
    if (malloc)
        delete malloc;
    malloc = 0;
}

#define DECLARE_PARENT_MALLOC(cls) \
typedef FastMallocTemplate<cls> Malloc##cls; \
template<> bool Malloc<cls>::initialized_ = false; \
template<> uli FastMallocTemplate<cls>::blocksize_ = 0; \
template<> uli FastMallocTemplate<cls>::nblocks_ = 0; \
template<> std::string FastMallocTemplate<cls>::classname_ = ""; \
template<> FastMalloc* FastMallocTemplate<cls>::malloc = 0;

#define DECLARE_SUB_MALLOC(cls, mcls) \
mcls* malloc_init_##mcls = 0; \
static Malloc##cls malloc_decl_##mcls (malloc_init_##mcls, #cls)

#define DECLARE_MALLOC(cls) \
    DECLARE_PARENT_MALLOC(cls); \
    DECLARE_SUB_MALLOC(cls,cls)

}

#endif // MALLOCIMPL_H
